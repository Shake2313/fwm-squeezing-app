"""
Clickable SVG canvas as a Streamlit component.

Streamlit cannot return a click coordinate on its own, which blocks every part of
the optical-table plan: selecting an optic, and dropping an instrument head on an
arbitrary point. This is the substrate that unblocks it.

Deliberately built without npm. `declare_component(path=...)` will serve any
static folder, and the component protocol underneath the official helper library
is three postMessage shapes, so `frontend/index.html` implements them directly.
The cost is writing ~200 lines of vanilla JS once; the benefit is that the
component is package data with no toolchain, no lockfile and no build step in the
deployment.

Shapes are plain dicts so the caller owns the drawing vocabulary:

    {"id": "pbs_1", "kind": "rect", "x": 120, "y": 140, "w": 44, "h": 44,
     "angle": 45, "label": "PBS", "clickable": True}
    {"kind": "path", "points": [[0, 162], [120, 162]], "stroke": "#E8873A",
     "width": 6, "dash": "6 5"}

The return value is one of

    {"kind": "optic", "id": "pbs_1", "seq": 3}
    {"kind": "point", "x": 412.0, "y": 88.5, "seq": 4}
    None                                        (nothing clicked yet)

`seq` exists because Streamlit only reruns when a component's value *changes*;
without it, clicking the same optic twice would be silently swallowed.
"""
from pathlib import Path

import streamlit.components.v1 as components

_FRONTEND_DIR = Path(__file__).resolve().parent / "frontend"

#: Returned `kind` values, so callers can branch without magic strings.
KIND_OPTIC = "optic"
KIND_POINT = "point"

_component = None


def _instance():
    global _component
    if _component is None:
        if not (_FRONTEND_DIR / "index.html").exists():
            raise FileNotFoundError(
                f"canvas frontend missing at {_FRONTEND_DIR}; it ships as "
                f"package data and must be present next to this module")
        _component = components.declare_component(
            "sabes_canvas", path=str(_FRONTEND_DIR))
    return _component


def svg_canvas(spec, *, selected=None, key=None, default=None):
    """Render `spec` and return the latest click event, or `default`.

    `selected` highlights one shape id without waiting for a round trip, so the
    caller can keep the canvas in step with state changed elsewhere (a dialog
    closing, a preset being applied).
    """
    return _instance()(spec=spec, selected=selected, key=key, default=default)


def make_spec(shapes, *, width=900, height=360, background="#F6F8FB",
              probe=None):
    """Assemble a canvas spec. Thin, but it keeps key names in one place."""
    return {
        "width": float(width),
        "height": float(height),
        "background": background,
        "shapes": list(shapes),
        "probe": list(probe) if probe is not None else None,
    }


def rect(shape_id, x, y, w, h, *, angle=0.0, label=None, title=None,
         fill="#E2E8F0", stroke="#1F4E79", stroke_width=1.6, clickable=True,
         label_dy=None):
    """A boxy optic: cube, plate, mount, cell."""
    return {
        "id": shape_id, "kind": "rect", "x": float(x), "y": float(y),
        "w": float(w), "h": float(h), "angle": float(angle),
        "label": label, "title": title, "fill": fill, "stroke": stroke,
        "stroke_width": stroke_width, "clickable": bool(clickable),
        "label_dy": label_dy,
    }


def circle(shape_id, x, y, r, *, label=None, title=None, fill="#E2E8F0",
           stroke="#1F4E79", stroke_width=1.6, clickable=True, label_dy=None):
    """A round optic: lens, iris, fiber tip."""
    return {
        "id": shape_id, "kind": "circle", "x": float(x), "y": float(y),
        "r": float(r), "label": label, "title": title, "fill": fill,
        "stroke": stroke, "stroke_width": stroke_width,
        "clickable": bool(clickable), "label_dy": label_dy,
    }


def beam(points, *, stroke="#E8873A", width=3.0, dash=None, opacity=1.0):
    """A beam segment. Never clickable: clicking a beam means clicking the table
    underneath it, which is how an instrument head gets placed on it."""
    return {
        "kind": "path", "points": [[float(x), float(y)] for x, y in points],
        "stroke": stroke, "width": float(width), "dash": dash,
        "opacity": float(opacity), "clickable": False,
    }
