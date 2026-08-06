"""
Canvas spike: proves the click round-trip before anything is built on it.

Reached at `?app=sabes&dev=canvas`. It stays out of the normal SABES page on
purpose but rides the same deployment, because the thing that actually needs
proving is that a hand-rolled component works on Streamlit Community Cloud, not
just on a laptop.

What it must demonstrate, in order:
  1. the component loads at all (no npm build, served as package data);
  2. clicking an optic returns its id;
  3. clicking empty table returns a coordinate in spec units, not screen pixels;
  4. clicking the SAME optic twice still triggers a rerun (the `seq` guard);
  5. the SVG scales with the browser window without breaking (3).

Delete this module once the real table view lands.
"""
import streamlit as st

from . import canvas as canvas_module   # the module, not the function

_SPEC_WIDTH = 900
_SPEC_HEIGHT = 300


def _demo_spec(probe=None):
    """A miniature beamline: source, two optics, a detector, one beam."""
    shapes = [
        canvas_module.beam(
            [(60, 150), (300, 150), (300, 90), (620, 90)],
            stroke="#E8873A", width=7),
        canvas_module.beam(
            [(300, 150), (620, 150)],
            stroke="#D9534F", width=3, dash="7 5"),
        canvas_module.rect("source", 60, 150, 46, 58, label="ECDL",
                           fill="#111827", stroke="#111827"),
        canvas_module.rect("pbs", 300, 150, 40, 40, angle=0, label="PBS",
                           fill="#F3D9F5"),
        canvas_module.circle("lens", 470, 90, 9, label="f=150",
                             fill="#E9D5FF"),
        canvas_module.rect("cell", 470, 150, 76, 34, label="Rb cell",
                           fill="#FDE9C9"),
        canvas_module.rect("pd", 640, 90, 26, 44, label="PD",
                           fill="#111827", stroke="#111827"),
        # A deliberately non-clickable element: steering mirrors are geometry,
        # not parameters, and the table view must be able to say so.
        canvas_module.rect("mirror_decor", 300, 90, 30, 8, angle=45,
                           label="M", fill="#BFDBFE", clickable=False),
    ]
    return canvas_module.make_spec(shapes, width=_SPEC_WIDTH,
                                   height=_SPEC_HEIGHT, probe=probe)


def render():
    st.markdown("### Canvas spike")
    st.caption(
        "Stage B of the optical-table plan. Click an optic, then click empty "
        "table. Both must come back to Python, and clicking the same optic "
        "twice must still register."
    )

    history = st.session_state.setdefault("_sabes_spike_log", [])
    probe = st.session_state.get("_sabes_spike_probe")

    event = canvas_module.svg_canvas(_demo_spec(probe=probe), key="sabes_spike")

    if event and (not history or event.get("seq") != history[-1].get("seq")):
        history.append(event)
        if event.get("kind") == canvas_module.KIND_POINT:
            st.session_state["_sabes_spike_probe"] = [event["x"], event["y"]]
            st.rerun()

    left, right = st.columns([2, 3])
    with left:
        st.markdown("**Latest event**")
        st.json(event if event else {"note": "nothing clicked yet"})
    with right:
        st.markdown("**Checks**")
        ids = {e.get("id") for e in history if e.get("kind") == canvas_module.KIND_OPTIC}
        # Scan the whole log, not just its tail: the interesting event is any
        # adjacent pair naming the same optic, and a later click on empty table
        # must not un-prove it.
        repeated = any(
            a.get("kind") == canvas_module.KIND_OPTIC
            and a.get("id") == b.get("id")
            for a, b in zip(history, history[1:]))
        rows = [
            ("component loaded", event is not None or bool(history)),
            ("optic click returns an id", bool(ids)),
            ("empty click returns a coordinate",
             any(e.get("kind") == "point" for e in history)),
            ("repeat click on one optic re-fires", repeated),
        ]
        st.markdown(
            "| Check | Result |\n|---|---|\n"
            + "".join(f"| {name} | {'PASS' if ok else 'not yet'} |\n"
                      for name, ok in rows))

    if history:
        st.markdown("**Event log** (newest first)")
        st.markdown(
            "| seq | kind | payload |\n|---|---|---|\n"
            + "".join(f"| {e.get('seq')} | {e.get('kind')} | {_payload(e)} |\n"
                      for e in reversed(history[-12:])))
        st.button("Clear log", on_click=_clear)


def _payload(event):
    if event.get("kind") == canvas_module.KIND_OPTIC:
        return event.get("id", "?")
    return f"({event.get('x')}, {event.get('y')})"


def _clear():
    st.session_state["_sabes_spike_log"] = []
    st.session_state.pop("_sabes_spike_probe", None)
