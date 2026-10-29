"""Static Matplotlib background plus browser-local SVG curves/cursors.

Only gesture completion is sent to Streamlit. The owned figure is restored
before returning, so exported arrays and figures retain the committed state.
"""
from __future__ import annotations

import base64
import hashlib
from io import BytesIO
import json
import logging
from pathlib import Path

import numpy as np
import streamlit as st
import matplotlib as mpl
from matplotlib.colors import to_hex

from gabes.plot_style import apply_gabes_plot_style
from .component_runtime import missing_apptest_component

_ASSETS = Path(__file__).resolve().parents[1] / "assets" / "ui"
_COMPONENT = st.components.v2.component(
    "gabes_live_plot",
    html='<div class="lp"><output class="lp-status" aria-live="polite"></output>'
         '<div class="lp-stage"><img alt="Spectrum" />'
         '<svg aria-label="Interactive spectrum"></svg></div>'
         '<div class="lp-controls"></div></div>',
    css="""
    .lp{width:100%;font:13px var(--g-font-sans,sans-serif);color:var(--g-ink-2,#334155)}
    .lp-stage{position:relative;line-height:0;margin-inline:auto}.lp-stage img{width:100%;display:block}
    .lp-stage svg{position:absolute;inset:0;width:100%;height:100%;touch-action:none;cursor:ew-resize}
    .lp-controls{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:8px 20px;margin-top:8px}
    .lp-controls label{display:grid;grid-template-columns:1fr 90px;gap:4px 8px;align-items:center}
    .lp-controls input[type=range]{grid-column:1/-1;width:100%;accent-color:#0284c7;margin:0}
    .lp-controls input[type=number]{min-width:0;width:100%;box-sizing:border-box;color:inherit;
      border:1px solid #cbd5e1;border-radius:4px;background:transparent;padding:3px;font:inherit}
    .lp-status{display:block;min-height:2.6em;margin:0 0 6px;line-height:1.3;
      font-variant-numeric:tabular-nums;color:var(--g-ink-2,#334155)}
    .lp[aria-busy=true] .lp-status{color:var(--g-accent,#0284c7)}
    @media(max-width:600px){.lp-controls{grid-template-columns:1fr}}
    """,
    js=(_ASSETS / "liveplot.js").read_text(encoding="utf-8"),
)
_MOUNTABLE = True


def _numbers(values):
    """JSON-safe arrays keep unavailable intervals as gaps, never as zero."""
    return [float(v) if np.isfinite(v) else None for v in values]


def _axis_data(axis, width, height):
    box = axis.get_position()
    return dict(left=box.x0 * width, top=(1 - box.y1) * height,
                width=box.width * width, height=box.height * height,
                xlim=list(map(float, axis.get_xlim())),
                ylim=list(map(float, axis.get_ylim())),
                log=axis.get_yscale() == "log")


def preview_data(fig, scheme_name, comparison=None, *, canonical_noise=False):
    """Extract only the mutable artists; no physics is solved here."""
    moving = []
    if scheme_name == "sas" and comparison:
        alignment = comparison["alignment"]
        axis = fig.axes[comparison["descriptor"]["axis_index"]]
        traces = []
        for line in axis.lines:
            if line.get_label().startswith(("Experimental CSV", "CSV ·")):
                moving.append(line)
                traces.append(dict(x=_numbers(line.get_xdata()), y=_numbers(line.get_ydata()),
                                   color=to_hex(line.get_color()), alpha=line.get_alpha() or 1,
                                   dashed=line.get_linestyle() == ":"))
        if not traces:
            return None, []
        xlim = axis.get_xlim()
        span = abs(xlim[1] - xlim[0])
        shift = alignment["shift"]
        controls = [
            dict(name="shift", label=f"X shift [{comparison['descriptor']['axis_x_unit']}]",
                 value=shift, min=min(xlim[0] - 2 * span, shift),
                 max=max(xlim[1] + 2 * span, shift), step=span / 1000),
            dict(name="multiplier", label="X fine scale ×", value=alignment["manual_multiplier"],
                 min=0.1, max=3.0, step=0.005),
        ]
        if comparison["calibration_mode"] == "relative_extrema":
            controls += [
                dict(name="signal_scale", label="Y scale ×", value=alignment["signal_scale"],
                     min=min(0.01, alignment["signal_scale"] / 10),
                     max=max(3., alignment["signal_scale"] * 2),
                     step=min(0.01, alignment["signal_scale"] / 100)),
                dict(name="signal_offset", label="Y offset", value=alignment["signal_offset"],
                     min=min(-1., alignment["signal_offset"] - 1),
                     max=max(1., alignment["signal_offset"] + 1), step=0.005),
            ]
        data = dict(kind="alignment", traces=traces, controls=controls,
                    original=dict(shift=shift, multiplier=alignment["manual_multiplier"],
                                  signal_scale=alignment["signal_scale"],
                                  signal_offset=alignment["signal_offset"]),
                    pivot=(min(traces[-1]["x"]) + max(traces[-1]["x"])) / 2 - shift,
                    axis_indices=[fig.axes.index(axis)])
    elif scheme_name == "fwm" and len(fig.axes) == 2:
        # Seeded FWM's two fixed curves each carry one selected-point marker.
        if not all(len(ax.lines) >= 3 and ax.collections for ax in fig.axes):
            return None, []
        curves = []
        for axis in fig.axes:
            line = axis.lines[0]
            curves.append(dict(x=_numbers(line.get_xdata()), y=_numbers(line.get_ydata())))
            moving.extend([axis.lines[1], axis.collections[0]])
        value = float(fig.axes[0].lines[1].get_xdata()[0])
        data = dict(kind="cursor", curves=curves, linear_noise=canonical_noise,
                    original=dict(tpd=value), controls=[], axis_indices=[0, 1],
                    event_key="fwm__tpd")
    else:
        return None, []
    return data, moving


def _commit(component_key, scheme_name, context):
    state = st.session_state.get(component_key, {})
    payload = state.get("gesture") if hasattr(state, "get") else getattr(state, "gesture", None)
    try:
        payload = json.loads(payload)
        if payload["context"] != context:
            return
        sequence = payload.get("n", 0)
        if isinstance(sequence, bool) or not isinstance(sequence, int) or sequence < 0:
            return
        previous = st.session_state.get(component_key + "-ack", {})
        if previous.get("context") == context and sequence <= previous.get("n", -1):
            return
        values = {k: float(v) for k, v in payload["values"].items()}
        if not all(np.isfinite(v) for v in values.values()):
            return
    except (ValueError, TypeError, KeyError, AttributeError):
        return
    if scheme_name == "sas" and any(st.session_state.get(f"sas__{name}", False)
                                    for name in ("_csv_apply", "_csv_frame", "_csv_reset", "_csv_fit")):
        # Native explicit alignment actions win if Streamlit receives their
        # widget states in the same rerun as an older live gesture.
        st.session_state[component_key + "-ack"] = {"context": context, "n": sequence}
        return
    if scheme_name == "fwm":
        if "tpd" not in values:
            return
        intent = payload.get("intent", 0)
        if isinstance(intent, bool) or not isinstance(intent, int) or intent < 0:
            return
        intent_key = "fwm__tpd-input-intent"
        if intent < st.session_state.get(intent_key, 0):
            # A newer explicit knob edit wins, regardless of callback order.
            # Acknowledge the discarded request so the browser can drain it.
            st.session_state[component_key + "-ack"] = {"context": context, "n": sequence}
            return
        st.session_state[intent_key] = intent
    if scheme_name == "sas":
        names = dict(shift="_csv_x_shift", multiplier="_csv_x_multiplier",
                     signal_scale="_csv_y_scale", signal_offset="_csv_y_offset")
        if not 0.1 <= values.get("multiplier", 1) <= 3 or values.get("signal_scale", 1) <= 0:
            return
        for name, value in values.items():
            if name in names:
                st.session_state[f"sas__{names[name]}"] = value
        st.session_state.pop("sas___csv_alignment_fit", None)
    elif scheme_name == "fwm" and "tpd" in values:
        from gabes.schemes.fwm import TPD_LIMIT_MHZ
        st.session_state["fwm__tpd"] = float(np.clip(values["tpd"], -TPD_LIMIT_MHZ, TPD_LIMIT_MHZ))
    else:
        return
    st.session_state[component_key + "-ack"] = {"context": context, "n": sequence}


def _preview_context(data, comparison, interaction_context=None):
    """Identity of the data being edited, independent of cursor and viewport.

    An acknowledgement for one spectrum/file must not replay queued gestures
    onto a different spectrum or a differently calibrated import.
    """
    if comparison:
        alignment = comparison.get("alignment", {})
        identity = dict(source=comparison["source_fingerprint_sha256"],
                        revision=comparison.get("alignment_revision", 0),
                        calibration=comparison.get("calibration_mode"),
                        inputs=comparison.get("calibration_inputs"),
                        descriptor=comparison.get("descriptor"),
                        fixed_alignment={k: alignment.get(k) for k in
                                         ("base_scale", "reverse", "invert")})
    elif interaction_context is not None:
        # FWM figures insert the exact selected point into their sampled
        # curves. Identify the physical settings, not that display-only knot.
        identity = dict(settings=interaction_context, linear_noise=data.get("linear_noise"))
    else:
        identity = dict(curves=data.get("curves"), linear_noise=data.get("linear_noise"))
    return hashlib.sha256(json.dumps(identity, sort_keys=True, allow_nan=False).encode()).hexdigest()


def render_live_plot(fig, scheme_name, comparison=None, *, canonical_noise=False, gain_digits=6,
                     interaction_context=None):
    global _MOUNTABLE
    if not _MOUNTABLE:
        return False
    data, moving = preview_data(fig, scheme_name, comparison, canonical_noise=canonical_noise)
    if data is None:
        return False
    apply_gabes_plot_style(fig)
    visibility = [artist.get_visible() for artist in moving]
    try:
        for artist in moving:
            artist.set_visible(False)
        stream = BytesIO()
        # Fixed canvas bounds are essential: bbox_inches='tight' changes the
        # origin and breaks alignment between the raster and SVG coordinates.
        with mpl.rc_context({"savefig.bbox": None}):
            fig.savefig(stream, format="png", dpi=120)
    finally:
        for artist, visible in zip(moving, visibility):
            artist.set_visible(visible)
    width, height = fig.get_size_inches() * 120
    data.update(image="data:image/png;base64," + base64.b64encode(stream.getvalue()).decode(),
                width=float(width), height=float(height),
                axes=[_axis_data(fig.axes[i], width, height) for i in data.pop("axis_indices")])
    context = _preview_context(data, comparison, interaction_context)
    data["context"] = context
    data["gain_digits"] = gain_digits
    component_key = f"live-plot-{scheme_name}"
    acknowledgement = st.session_state.get(component_key + "-ack", {})
    data["ack"] = (acknowledgement.get("n", 0)
                   if acknowledgement.get("context") == context else 0)
    data["scheme"] = scheme_name
    try:
        _COMPONENT(key=component_key, data=data, default={"gesture": None}, width="stretch",
                   on_gesture_change=lambda: _commit(component_key, scheme_name, context))
    except TypeError:
        # AppTest stubs the v2 component registry; keep its native-widget path.
        _MOUNTABLE = False
        logging.getLogger(__name__).warning("Live plot unavailable; using static plot", exc_info=True)
        return False
    except ValueError as error:
        if not missing_apptest_component(error, "gabes_live_plot"):
            raise
        # AppTest owns a fresh registry while this module's callable survives.
        return False
    return True
