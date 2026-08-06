"""
Draw a layout: optics as glyphs, beams as strokes carrying live power.

Presentation lives here and nowhere else. `Node` knows what a thing *is* and
which parameters it owns; how it looks is this module's business, so the physics
model never accumulates styling.

Two conventions worth stating out loud, because both are places a drawing can
quietly lie:

  Beam thickness is NOT to scale. A 530 um beam on a metre of table is 0.05 % of
  the width -- drawn honestly it would be invisible. Thickness encodes *power*,
  logarithmically, and the legend says so.

  Colour follows the setup drawings themselves: orange for the high-power pump,
  red for the seed and the probe it becomes, blue for the conjugate. Reusing the
  lab's own key means nobody has to translate.

Modelled optics are drawn saturated and are clickable; hardware the physics does
not know about is drawn muted and inert. `Node.clickable` decides both, so the
picture cannot claim more than the model delivers.
"""
import math

from ..components import canvas as canvas_module
from .table import (BEAM_COMMON, BEAM_CONJUGATE, BEAM_PROBE, BEAM_PUMP,
                    BEAM_REFERENCE, BEAM_SEED)

# --- palette -------------------------------------------------------------
TABLE_BG = "#F1F5F9"
GRID_COLOR = "#CBD5E1"
INK = "#1E293B"
MUTED_INK = "#94A3B8"

BEAM_STYLE = {
    BEAM_COMMON: {"stroke": "#E8873A", "dash": None},
    BEAM_PUMP: {"stroke": "#E8873A", "dash": None},
    BEAM_SEED: {"stroke": "#D9534F", "dash": "7 5"},
    BEAM_PROBE: {"stroke": "#D9534F", "dash": None},
    BEAM_CONJUGATE: {"stroke": "#3B82F6", "dash": None},
    BEAM_REFERENCE: {"stroke": "#F0A868", "dash": "3 4"},
}

GLASS = "#E9D5FF"
CUBE = "#F3D9F5"
METAL = "#BFDBFE"
VAPOUR = "#FDE9C9"
DARK = "#1F2937"

# --- beam thickness ------------------------------------------------------
#: Power that maps to `WIDTH_AT_REF`, and pixels gained per decade above it.
#: Chosen so a microwatt seed stays visible and a watt of pump reads as heavy
#: without swamping the optics it passes through.
WIDTH_REF_W = 1e-6
WIDTH_AT_REF = 2.6
WIDTH_PER_DECADE = 1.25
WIDTH_MIN = 1.4
WIDTH_MAX = 10.0


def beam_width(power_w):
    """Stroke width for an optical power. Logarithmic, and not a beam size."""
    power = max(float(power_w), 1e-18)
    width = WIDTH_AT_REF + WIDTH_PER_DECADE * math.log10(power / WIDTH_REF_W)
    return min(max(width, WIDTH_MIN), WIDTH_MAX)


# --- glyphs --------------------------------------------------------------
def glyph_for(node):
    """Which glyph draws a node.

    Matched on the id's own vocabulary rather than a per-node styling field, so
    adding an optic to the layout does not require touching the renderer, and a
    node that matches nothing still draws as a plain box.
    """
    node_id, kind = node.id, node.kind
    if kind == "mirror" or node_id.startswith("m_"):
        return "mirror"
    if node_id.startswith("dm_"):
        return "dmirror"
    if node_id.startswith("pbs"):
        return "cube"
    if node_id.startswith(("hwp", "qwp")):
        return "waveplate"
    if node_id.startswith("etalon"):
        return "etalon"
    if node_id.startswith("iris"):
        return "iris"
    if node_id.startswith(("gt", "pol")):
        return "polarizer"
    if node_id.startswith("oi_"):
        return "isolator"
    if kind == "cell":
        return "cell"
    if kind == "fiber":
        return "collimator"
    if kind == "detector":
        return "detector"
    if kind == "source":
        return "box"
    if node.shape == "circle":
        return "lens"
    return "box"


def _parts_for(glyph, node, active):
    """Sub-primitives inside a node's clickable group."""
    ink = INK if active else MUTED_INK
    half_w, half_h = node.w / 2.0, node.h / 2.0
    if glyph == "cube":
        return [{"kind": "line", "x1": -half_w, "y1": -half_h,
                 "x2": half_w, "y2": half_h, "stroke": ink,
                 "stroke_width": 1.3}]
    if glyph == "waveplate":
        return [{"kind": "line", "x1": 0, "y1": -half_h + 2,
                 "x2": 0, "y2": half_h - 2, "stroke": ink, "stroke_width": 1.0}]
    if glyph == "etalon":
        return [{"kind": "line", "x1": dx, "y1": -half_h, "x2": dx,
                 "y2": half_h, "stroke": ink, "stroke_width": 1.1}
                for dx in (-half_w * 0.45, half_w * 0.45)]
    if glyph == "iris":
        return [{"kind": "circle", "r": max(node.w * 0.22, 2.0),
                 "fill": TABLE_BG, "stroke": ink, "stroke_width": 1.1}]
    if glyph == "polarizer":
        return [{"kind": "line", "x1": -half_w * 0.5, "y1": dy,
                 "x2": half_w * 0.5, "y2": dy, "stroke": ink,
                 "stroke_width": 0.8}
                for dy in (-half_h * 0.4, 0.0, half_h * 0.4)]
    if glyph == "isolator":
        return [{"kind": "line", "x1": -4, "y1": 0, "x2": 4, "y2": 0,
                 "stroke": "#F8FAFC", "stroke_width": 1.6}]
    if glyph == "detector":
        return [{"kind": "rect", "dx": -node.w * 0.28, "w": node.w * 0.3,
                 "h": node.h * 0.55, "fill": "#64748B", "rx": 1}]
    if glyph == "collimator":
        return [{"kind": "rect", "dx": node.w * 0.34, "w": node.w * 0.28,
                 "h": node.h * 0.34, "fill": "#4ADE80", "rx": 1}]
    if glyph == "dmirror":
        return [{"kind": "line", "x1": -node.w / 2, "y1": 0,
                 "x2": 0, "y2": 0, "stroke": "#F8FAFC", "stroke_width": 2.2}]
    return []


_FILL = {
    "cube": CUBE, "waveplate": GLASS, "etalon": "#E2E8F0", "lens": GLASS,
    "iris": "#CBD5E1", "polarizer": "#DDD6FE", "mirror": METAL,
    "dmirror": METAL, "cell": VAPOUR, "collimator": DARK, "detector": DARK,
    "box": DARK, "isolator": "#334155",
}


def node_shape(node, *, selected=False):
    """One optic as a canvas shape, glyph and styling included."""
    glyph = glyph_for(node)
    active = node.clickable
    fill = _FILL.get(glyph, "#E2E8F0")
    stroke = INK if active else MUTED_INK
    common = dict(
        label=node.label or None,
        title=(node.note or node.label or node.id) + (
            "" if active else "  (drawn, not modelled)"),
        fill=fill,
        stroke="#0F766E" if selected else stroke,
        stroke_width=2.6 if selected else (1.5 if active else 1.0),
        clickable=active,
        label_dy=-(node.h / 2.0) - 5.0,
    )
    if node.shape == "circle":
        shape = canvas_module.circle(node.id, node.x, node.y,
                                     max(node.w, node.h) / 2.0, **common)
    else:
        shape = canvas_module.rect(node.id, node.x, node.y, node.w, node.h,
                                   angle=node.angle, **common)
    shape["parts"] = _parts_for(glyph, node, active)
    if not active:
        shape["opacity"] = 0.55
    return shape


# --- assembly ------------------------------------------------------------
def _link_power(chain, link, twin_states, resolve):
    state = resolve(chain, link, twin_states)
    return None if state is None else state.total_power_w


def build_spec(layout, chain, *, selected=None, probe=None, twin_states=None):
    """A canvas spec for one layout, with beam widths from the live chain.

    Beams are drawn first so optics sit on top of them, and each beam gets a
    translucent underlay a couple of pixels wider -- the cheapest way to read as
    light rather than as wire, with no SVG filters involved.
    """
    from .probe import resolve_link_beam

    shapes = []
    for link in layout.links:
        style = BEAM_STYLE.get(link.beam, BEAM_STYLE[BEAM_COMMON])
        power = _link_power(chain, link, twin_states, resolve_link_beam)
        width = beam_width(power) if power is not None else WIDTH_MIN
        points = layout.points(link)
        shapes.append(canvas_module.beam(
            points, stroke=style["stroke"], width=width + 2.5,
            dash=style["dash"], opacity=0.18))
        shapes.append(canvas_module.beam(
            points, stroke=style["stroke"], width=width, dash=style["dash"],
            opacity=1.0 if power is not None else 0.45))

    for node in layout.nodes:
        shapes.append(node_shape(node, selected=(node.id == selected)))

    spec = canvas_module.make_spec(
        shapes, width=layout.width, height=layout.height,
        background=TABLE_BG, probe=probe)
    spec["grid"] = 24.0
    spec["grid_color"] = GRID_COLOR
    return spec


def legend_rows():
    """(swatch colour, dashed, label) for the caption under the drawing."""
    return (
        (BEAM_STYLE[BEAM_PUMP]["stroke"], False, "pump"),
        (BEAM_STYLE[BEAM_SEED]["stroke"], True, "seed"),
        (BEAM_STYLE[BEAM_PROBE]["stroke"], False, "probe"),
        (BEAM_STYLE[BEAM_CONJUGATE]["stroke"], False, "conjugate"),
        (BEAM_STYLE[BEAM_REFERENCE]["stroke"], True, "locking reference"),
    )
