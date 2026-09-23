"""ScrubField — the precise numeric knob (docs/ui_redesign/decisions.md D3).

One row per knob: the label carries its help on hover, the value is a
spinbutton you can drag, type into or step with the keyboard, and a 3 px track
underneath gives the coarse sweep a slider used to give.

Division of labour: the browser owns the gesture and the Python side owns the
number. A typed entry travels as raw text and `gabes_ui.units` parses it, so
unit handling has one implementation and one test suite. Anything the runtime
cannot mount (an older Streamlit, or `AppTest`, which stubs the component
manager) falls back to `st.slider` with the same session-state key.
"""
import contextlib
import json
import logging
from html import escape

import streamlit as st

from gabes_ui.units import clamp, decimals_for, display_text, parse_entry

NOTE_SUFFIX = "__scrub_note"

_LOG = logging.getLogger(__name__)


def state_key(key):
    """Widget key for the component instance behind the knob ``key``.

    Streamlit reserves ``__`` inside a bidirectional component's id, and every
    knob key already contains one (``scheme__param``), so it is folded out.
    """
    return "scrub-" + key.replace("__", "-")


_HTML = """<div class="sf">
  <div class="sf-head">
    <span class="sf-label"></span>
    <button class="sf-reset" type="button" tabindex="-1" hidden></button>
    <span class="sf-value" role="spinbutton" tabindex="0"><span
      class="sf-num"></span><span class="sf-unit"></span></span>
    <input class="sf-input" type="text" inputmode="decimal" autocomplete="off"
           spellcheck="false" hidden>
  </div>
  <div class="sf-rail">
    <span class="sf-end sf-end-l"></span>
    <span class="sf-track"><span class="sf-fill"></span><span class="sf-knob"></span></span>
    <span class="sf-end sf-end-r"></span>
  </div>
  <div class="sf-note" role="status" hidden></div>
  <span class="sf-help" id="sf-help" hidden></span>
</div>"""

_CSS = """
.sf{font:400 .8rem/1.3 var(--g-font-sans,sans-serif);color:var(--g-ink-2,#1F2D3D);
  padding:.18rem 0 .1rem;}
.sf-head{display:flex;align-items:baseline;gap:.35rem;}
.sf-label{flex:1 1 auto;min-width:0;color:var(--g-ink-3,#3B4A5E);
  overflow:hidden;text-overflow:ellipsis;white-space:nowrap;}
.sf-label[title]{cursor:help;text-decoration:underline dotted transparent;
  text-underline-offset:3px;}
.sf-label[title]:hover{text-decoration-color:var(--g-line-strong,#CBD5E1);}
.sf-value{flex:0 0 auto;display:inline-flex;align-items:baseline;gap:.15rem;
  font:500 .82rem/1.3 var(--g-font-mono,monospace);color:var(--g-ink,#0E1726);
  padding:0 .2rem;border-radius:4px;cursor:ew-resize;touch-action:none;
  user-select:none;-webkit-user-select:none;white-space:nowrap;}
.sf-value:hover{background:var(--g-fill-soft,#EDF1F5);}
.sf-value:focus-visible{outline:2px solid var(--g-accent,#0284C7);outline-offset:1px;}
.sf-unit{font-size:.72rem;color:var(--g-muted,#5B6B80);}
.sf-input{flex:0 0 7.5rem;width:7.5rem;box-sizing:border-box;
  font:500 .82rem/1.3 var(--g-font-mono,monospace);color:var(--g-ink,#0E1726);
  text-align:right;padding:.1rem .3rem;border:1px solid var(--g-accent,#0284C7);
  border-radius:4px;background:var(--g-surface,#fff);}
.sf-input:focus{outline:none;}
.sf-reset{flex:0 0 auto;width:12px;height:12px;padding:0;margin:0;border:0;
  border-radius:50%;background:transparent;cursor:pointer;position:relative;
  align-self:center;}
.sf-reset::before{content:"";position:absolute;inset:3px;border-radius:50%;
  background:var(--g-accent,#0284C7);}
.sf-reset:hover::before{inset:1px;}
.sf-reset:focus-visible{outline:2px solid var(--g-accent,#0284C7);outline-offset:1px;}
.sf-rail{display:flex;align-items:center;gap:.35rem;height:12px;}
.sf-end{flex:0 0 auto;font-size:.65rem;line-height:1;color:var(--g-muted,#5B6B80);
  white-space:nowrap;}
.sf-end:empty{display:none;}
.sf-track{position:relative;flex:1 1 auto;height:3px;border-radius:2px;
  background:var(--g-line,#E3E8EF);cursor:pointer;touch-action:none;}
.sf-fill{position:absolute;left:0;top:0;height:100%;border-radius:2px;
  background:var(--g-line-strong,#CBD5E1);}
.sf-knob{position:absolute;top:50%;width:9px;height:9px;border-radius:50%;
  background:var(--g-accent,#0284C7);transform:translate(-50%,-50%);}
.sf.is-dragging .sf-knob{box-shadow:0 0 0 4px var(--g-accent-soft,#E6F2FA);}
.sf.no-range .sf-rail{display:none;}
.sf-note{font-size:.7rem;line-height:1.3;color:var(--g-warn-ink,#92400E);
  padding-top:.1rem;}
.sf-help{position:absolute;width:1px;height:1px;overflow:hidden;clip:rect(0 0 0 0);
  white-space:nowrap;}
@media (pointer:coarse){
  .sf-value{padding:.25rem .35rem;}
  .sf-rail{height:22px;}
  .sf-knob{width:14px;height:14px;}
}
"""

_JS = r"""
export default function (component) {
  const { parentElement, data, setStateValue } = component;
  const root = parentElement.querySelector(".sf");
  if (!root) return;
  const q = (c) => root.querySelector(c);
  const el = {
    label: q(".sf-label"), value: q(".sf-value"), num: q(".sf-num"),
    unit: q(".sf-unit"), input: q(".sf-input"), reset: q(".sf-reset"),
    track: q(".sf-track"), fill: q(".sf-fill"), knob: q(".sf-knob"),
    endL: q(".sf-end-l"), endR: q(".sf-end-r"), note: q(".sf-note"),
    help: q(".sf-help"),
  };
  const S = root._sf || (root._sf = { seq: 0, editing: false });
  const min = data.min, max = data.max, step = data.step || 0;
  const span = (min != null && max != null) ? max - min : 0;
  const hasRange = span > 0;
  root.classList.toggle("no-range", !hasRange);

  const tidy = (v) => Number(v.toFixed(12));
  const snap = (v) => {
    if (!step) return tidy(v);
    const base = min != null ? min : 0;
    return tidy(base + Math.round((v - base) / step) * step);
  };
  const bound = (v) => {
    if (min != null && v < min) v = min;
    if (max != null && v > max) v = max;
    return tidy(v);
  };
  const show = (v) => data.decimals > 0 ? v.toFixed(data.decimals) : String(tidy(v));

  // `text` comes from Python, which spells out a typed value finer than the
  // format (1.234567 under "%.2f"); a value the browser just moved is on the
  // step grid, so the format is exact for it.
  function paint(v, text) {
    S.shown = v;
    if (!S.editing) el.num.textContent = text != null ? text : show(v);
    el.value.setAttribute("aria-valuenow", String(v));
    el.value.setAttribute("aria-valuetext", show(v) + (data.unit ? " " + data.unit : ""));
    if (hasRange) {
      const pct = Math.max(0, Math.min(1, (v - min) / span)) * 100;
      el.fill.style.width = pct + "%";
      el.knob.style.left = pct + "%";
    }
    const off = data.default != null && Math.abs(v - data.default) > (step ? step / 2 : 1e-12);
    el.reset.hidden = !off;
  }

  function send(payload) {
    // A JSON string, not an object: the sequence number makes two identical
    // gestures two distinct values, so neither is swallowed as "no change".
    S.seq += 1;
    setStateValue("value", JSON.stringify(Object.assign({ n: S.seq }, payload)));
  }
  function stream(v) {           // navigate-only knobs follow the drag
    if (!data.live) return;
    const now = Date.now();
    if (S.sentAt && now - S.sentAt < 125) return;
    S.sentAt = now;
    send({ v: v });
  }

  el.label.textContent = data.label || "";
  if (data.help) {
    el.label.title = data.help;
    el.help.textContent = data.help;
    el.help.hidden = false;
    el.value.setAttribute("aria-describedby", "sf-help");
  } else {
    // A label too long for the rail is clipped; hovering must still reveal it.
    const clipped = el.label.scrollWidth > el.label.clientWidth + 1;
    if (clipped) el.label.title = data.label || "";
    else el.label.removeAttribute("title");
    el.help.hidden = true;
  }
  el.unit.textContent = data.unit || "";
  el.endL.textContent = (data.ends && data.ends[0]) || "";
  el.endR.textContent = (data.ends && data.ends[1]) || "";
  el.value.setAttribute("aria-label", data.label + (data.unit ? ", " + data.unit : ""));
  if (min != null) el.value.setAttribute("aria-valuemin", String(min));
  if (max != null) el.value.setAttribute("aria-valuemax", String(max));
  el.reset.title = "Reset to " + (data.defaultText || "default");
  el.reset.setAttribute("aria-label", el.reset.title);
  el.note.textContent = data.note || "";
  el.note.hidden = !data.note;
  paint(data.value, data.text);

  // ---- typing ---------------------------------------------------------
  function beginEdit() {
    S.editing = true;
    el.input.value = el.num.textContent;
    el.value.hidden = true;
    el.input.hidden = false;
    el.input.focus();
    el.input.select();
  }
  function endEdit(commit) {
    if (!S.editing) return;
    S.editing = false;
    el.input.hidden = true;
    el.value.hidden = false;
    if (commit) send({ t: el.input.value });
    else el.value.focus();
  }
  el.input.onkeydown = (e) => {
    if (e.key === "Enter") { e.preventDefault(); endEdit(true); }
    else if (e.key === "Escape") { e.preventDefault(); endEdit(false); }
  };
  el.input.onblur = () => endEdit(true);

  // ---- dragging -------------------------------------------------------
  function fromTrack(clientX) {
    const r = el.track.getBoundingClientRect();
    if (!r.width || !hasRange) return S.shown;
    return bound(snap(min + ((clientX - r.left) / r.width) * span));
  }
  el.value.onpointerdown = (e) => {
    if (S.editing || e.button !== 0) return;
    e.preventDefault();
    el.value.setPointerCapture(e.pointerId);
    S.drag = { x: e.clientX, v: S.shown, moved: false };
    root.classList.add("is-dragging");
  };
  el.value.onpointermove = (e) => {
    const d = S.drag;
    if (!d) return;
    const dx = e.clientX - d.x;
    if (!d.moved && Math.abs(dx) < 3) return;
    d.moved = true;
    let per = hasRange ? Math.max(step || 0, span / 400) : (step || 1);
    if (e.shiftKey) per = Math.max(per * 0.1, step ? step : per * 0.1);
    const v = bound(snap(d.v + dx * per));
    paint(v);
    stream(v);
  };
  el.value.onpointerup = (e) => {
    const d = S.drag;
    S.drag = null;
    root.classList.remove("is-dragging");
    if (el.value.hasPointerCapture(e.pointerId)) el.value.releasePointerCapture(e.pointerId);
    if (!d) return;
    if (!d.moved) beginEdit();
    else if (S.shown !== d.v) send({ v: S.shown });
  };
  el.value.onpointercancel = () => { S.drag = null; root.classList.remove("is-dragging"); };

  el.track.onpointerdown = (e) => {
    if (!hasRange || e.button !== 0) return;
    e.preventDefault();
    el.track.setPointerCapture(e.pointerId);
    S.track = { v: S.shown };
    root.classList.add("is-dragging");
    const v = fromTrack(e.clientX);
    paint(v);
    stream(v);
  };
  el.track.onpointermove = (e) => {
    if (!S.track) return;
    const v = fromTrack(e.clientX);
    paint(v);
    stream(v);
  };
  el.track.onpointerup = (e) => {
    const t = S.track;
    S.track = null;
    root.classList.remove("is-dragging");
    if (el.track.hasPointerCapture(e.pointerId)) el.track.releasePointerCapture(e.pointerId);
    if (t && S.shown !== t.v) send({ v: S.shown });
  };
  el.track.onpointercancel = () => { S.track = null; root.classList.remove("is-dragging"); };

  // ---- keyboard -------------------------------------------------------
  el.value.onkeydown = (e) => {
    if (S.editing) return;
    if (e.key === "Enter" || e.key === " ") { e.preventDefault(); beginEdit(); return; }
    const unit = (step || (hasRange ? span / 100 : 1)) * (e.shiftKey ? 0.1 : 1);
    let v = null;
    if (e.key === "ArrowUp" || e.key === "ArrowRight") v = S.shown + unit;
    else if (e.key === "ArrowDown" || e.key === "ArrowLeft") v = S.shown - unit;
    else if (e.key === "PageUp") v = S.shown + unit * 10;
    else if (e.key === "PageDown") v = S.shown - unit * 10;
    else if (e.key === "Home" && min != null) v = min;
    else if (e.key === "End" && max != null) v = max;
    else return;
    e.preventDefault();
    paint(bound(tidy(v)));
    clearTimeout(S.keyTimer);            // a key repeat is one gesture
    S.keyTimer = setTimeout(() => send({ v: S.shown }), 160);
  };

  el.reset.onclick = (e) => { e.preventDefault(); send({ r: true }); };
}
"""

try:
    _COMPONENT = st.components.v2.component(
        "gabes_scrub", html=_HTML, css=_CSS, js=_JS)
except Exception:                                       # pragma: no cover
    _COMPONENT = None

_MOUNTABLE = _COMPONENT is not None


def available():
    """False once mounting has failed in this process (AppTest, old runtime)."""
    return _MOUNTABLE


def _within(container):
    """`with` context for a container that may be the `st` module itself."""
    if container is None or container is st:
        return contextlib.nullcontext()
    return container


def decode_gesture(raw):
    """One gesture as the component sends it, or None if there is nothing to do.

    ``{"n": seq}`` plus either ``v`` (a number the browser already snapped),
    ``t`` (text the person typed) or ``r`` (reset to the default).
    """
    if isinstance(raw, str):
        try:
            raw = json.loads(raw)
        except ValueError:                              # pragma: no cover
            return None
    return raw if isinstance(raw, dict) else None


def resolve_gesture(payload, *, unit, vmin, vmax, integer, default):
    """What one gesture means for the knob: ``(value, note)``.

    ``value`` is None when the entry could not be trusted — the knob keeps
    what it had and ``note`` says why. A note without a refusal (a clamp)
    still moves the knob.
    """
    if not isinstance(payload, dict):
        return None, None
    if payload.get("r"):
        value = default
    elif "t" in payload:
        value, error = parse_entry(payload.get("t"), unit)
        if error:
            return None, error
    else:
        try:
            value = float(payload.get("v"))
        except (TypeError, ValueError):                 # pragma: no cover
            return None, None
    if value is None:                                   # pragma: no cover
        return None, None
    value, note = clamp(value, vmin, vmax)
    if integer:
        value = int(round(value))
    return value, note


def _commit(key, comp_key, unit, vmin, vmax, integer, default, on_change):
    """Turn one gesture into the knob's value. Runs before the script rerun."""
    state = st.session_state.get(comp_key)
    raw = None
    if state is not None:
        try:
            raw = state["value"]
        except (KeyError, TypeError):                   # pragma: no cover
            raw = None
    payload = decode_gesture(raw)
    if payload is None:
        return
    value, note = resolve_gesture(payload, unit=unit, vmin=vmin, vmax=vmax,
                                  integer=integer, default=default)
    if value is None:
        if note:
            st.session_state[key + NOTE_SUFFIX] = note
        return
    st.session_state[key] = value
    st.session_state[key + NOTE_SUFFIX] = note or ""
    if on_change is not None:
        on_change()


def _fallback(container, key, label, unit, vmin, vmax, step, fmt, help_, endpoints,
              on_change):
    """The plain slider, for runtimes that cannot mount the component."""
    shown = label + (f"  [{unit}]" if unit else "")
    kwargs = {"format": fmt} if fmt else {}
    value = container.slider(shown, vmin, vmax, step=step, key=key, help=help_,
                             on_change=on_change, **kwargs)
    if endpoints:
        left, right = endpoints
        container.markdown(
            "<div class='gabes-endpoints'>"
            f"<span>{escape(str(left))}</span><span>{escape(str(right))}</span></div>",
            unsafe_allow_html=True)
    return value


def scrub_field(container, *, key, label, value, vmin, vmax, step, unit="",
                default=None, fmt=None, help=None, endpoints=None, live=False,
                integer=False, on_change=None):
    """Render one precise numeric knob and return its value.

    `st.session_state[key]` is the value — a plain entry, not widget state, so
    presets and default sets write to it and a knob that stops being rendered
    keeps what it had.
    """
    global _MOUNTABLE
    if key not in st.session_state:
        st.session_state[key] = value
    current = st.session_state[key]
    if not _MOUNTABLE:
        return _fallback(container, key, label, unit, vmin, vmax, step, fmt, help,
                         endpoints, on_change)
    # Read and clear: a notice explains the gesture that just happened, and
    # would be stale by the next one.
    note = st.session_state.pop(key + NOTE_SUFFIX, "") or ""
    data = {
        "label": label,
        "unit": unit or "",
        "value": current,
        "text": display_text(current, fmt),
        "min": vmin,
        "max": vmax,
        "step": step,
        "decimals": decimals_for(fmt),
        "default": default,
        "defaultText": display_text(default, fmt),
        "help": help or "",
        "ends": list(endpoints) if endpoints else None,
        "live": bool(live),
        "note": note,
    }
    comp_key = state_key(key)
    try:
        with _within(container):
            _COMPONENT(
                key=comp_key,
                data=data,
                default={"value": None},
                width="stretch",
                on_value_change=lambda: _commit(
                    key, comp_key, unit, vmin, vmax, integer, default,
                    on_change),
            )
    except Exception:
        # One failure settles it for the process: AppTest stubs the component
        # manager, so every later mount would fail the same way.
        _MOUNTABLE = False
        _LOG.warning("ScrubField unavailable, falling back to sliders",
                     exc_info=True)
        return _fallback(container, key, label, unit, vmin, vmax, step, fmt, help,
                         endpoints, on_change)
    return st.session_state[key]
