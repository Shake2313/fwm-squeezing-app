"""
The optical table as data: where things sit, and how far the light travels.

Two coordinate systems, deliberately independent.

  display   `Node.x/y/angle` and `Link.via`, in the arbitrary units of the setup
            drawings. Used to render the table. NEVER enters the physics.
  path      `Link.path_m`, the optical path length of a segment in metres, with a
            provenance tag. Used by the physics. NEVER derived from display.

They are allowed to disagree, because the setup drawings are schematics: in
part 2 the pump telescope lenses must sit `f1 + f2 = 450 mm` apart, and nothing
makes the pixel distance agree. Deriving path lengths from drawing coordinates
would quietly corrupt every beam size downstream. Keeping the two apart means a
schematic stays a schematic and the numbers stay honest.

Only segments where the beam meaningfully diverges need a real `path_m`; the
rest carry `None` and are geometry for the eye alone. That is why this does not
require surveying the table with a tape measure.
"""
from dataclasses import dataclass
from typing import Optional, Tuple

import math

#: Which beam a segment carries. Doubles as the drawing's colour key.
BEAM_PUMP = "pump"
BEAM_SEED = "seed"
BEAM_PROBE = "probe"
BEAM_CONJUGATE = "conjugate"
BEAM_COMMON = "common"
BEAM_REFERENCE = "reference"

#: Where a `path_m` came from. Same vocabulary as `sabes.calibration`, so the UI
#: can badge a distance the same way it badges a coefficient.
PROVENANCE_DESIGN = "design"
PROVENANCE_MEASURED = "measured"
PROVENANCE_NOMINAL = "nominal"


@dataclass(frozen=True)
class Node:
    """One thing on the table.

    Two independent honesty flags, and they answer different questions.

    `modelled`  does the physics know what this element *does*? A half-wave
                plate rotates polarization; an optical isolator, as far as this
                model is concerned, does nothing but attenuate.
    `lumped`    is a single transmission the whole story? Anything on the table
                that is not detail-modelled is drawn with "(lumped)" after its
                label, so the picture never implies more than the model has.

    Every element downstream of the amplifier carries `transmission_key`, the
    calibration coefficient holding its per-pass transmission. That is what makes
    alignment loss a variable you can set per optic rather than one lumped number
    per arm, and it is why a lumped optic is still worth clicking.
    """
    id: str
    label: str
    x: float
    y: float
    kind: str = "optic"              # optic | source | detector | cell | fiber
    shape: str = "rect"              # rect | circle
    w: float = 26.0
    h: float = 26.0
    angle: float = 0.0
    modelled: bool = False
    lumped: bool = False
    params: Tuple[str, ...] = ()     # settings fields this optic owns
    transmission_key: str = ""       # calibration coefficient, "" = not tracked
    note: str = ""

    @property
    def clickable(self):
        """Anything with something to edit: a physics knob or a transmission."""
        return bool(self.params) or bool(self.transmission_key)

    @property
    def display_label(self):
        return f"{self.label} (lumped)" if self.lumped and self.label else self.label


@dataclass(frozen=True)
class Link:
    """A beam segment between two nodes.

    `carries` names the source-chain stage whose beam state travels this
    segment, which is how a point on the drawing is resolved to real physics.
    Segments downstream of the cell name a twin-beam arm instead.
    """
    a: str
    b: str
    beam: str = BEAM_COMMON
    path_m: Optional[float] = None
    provenance: str = PROVENANCE_NOMINAL
    via: Tuple[Tuple[float, float], ...] = ()
    carries: str = ""
    note: str = ""

    def __post_init__(self):
        if self.path_m is not None and self.path_m < 0.0:
            raise ValueError(f"{self.a}->{self.b}: negative path length")
        if self.provenance not in (PROVENANCE_DESIGN, PROVENANCE_MEASURED,
                                   PROVENANCE_NOMINAL):
            raise ValueError(f"{self.a}->{self.b}: bad provenance "
                             f"{self.provenance!r}")


@dataclass(frozen=True)
class Layout:
    """One drawing's worth of table: part 1 or part 2."""
    key: str
    title: str
    width: float
    height: float
    nodes: Tuple[Node, ...] = ()
    links: Tuple[Link, ...] = ()

    def __post_init__(self):
        ids = [n.id for n in self.nodes]
        duplicates = {i for i in ids if ids.count(i) > 1}
        if duplicates:
            raise ValueError(f"{self.key}: duplicate node ids {sorted(duplicates)}")
        known = set(ids)
        for link in self.links:
            missing = {link.a, link.b} - known
            if missing:
                raise ValueError(
                    f"{self.key}: link {link.a}->{link.b} references unknown "
                    f"node(s) {sorted(missing)}")

    # ---- lookup ----
    def node(self, node_id) -> Node:
        for candidate in self.nodes:
            if candidate.id == node_id:
                return candidate
        raise KeyError(f"{self.key}: no node {node_id!r}")

    def clickable_nodes(self):
        return tuple(n for n in self.nodes if n.clickable)

    def transmission_nodes(self):
        """Everything downstream of the amplifier, in declaration order."""
        return tuple(n for n in self.nodes if n.transmission_key)

    def transmission_along(self, *node_ids, calibration=None, default=1.0):
        """Product of the per-optic transmissions over a route.

        Replaces the one lumped loss coefficient an arm used to carry, so that
        "where did the power go" has an answer per optic instead of per arm.
        """
        total = 1.0
        for node_id in node_ids:
            key = self.node(node_id).transmission_key
            if not key:
                continue
            total *= (default if calibration is None
                      else calibration.value(key, default))
        return total

    def links_for(self, beam):
        return tuple(l for l in self.links if l.beam == beam)

    # ---- geometry ----
    def points(self, link) -> Tuple[Tuple[float, float], ...]:
        """Display polyline for a link: start, waypoints, end."""
        start = self.node(link.a)
        end = self.node(link.b)
        return ((start.x, start.y),) + tuple(link.via) + ((end.x, end.y),)

    def path_length_m(self, *node_ids, default=None):
        """Summed optical path along a chain of nodes.

        Raises rather than guessing when a segment on the requested route has no
        measured or designed length, because silently substituting zero is how a
        beam size goes wrong without anyone noticing.
        """
        total = 0.0
        for a, b in zip(node_ids, node_ids[1:]):
            link = self.link_between(a, b)
            if link.path_m is None:
                if default is None:
                    raise ValueError(
                        f"{self.key}: segment {a}->{b} has no path length; give "
                        f"it one or pass an explicit default")
                total += float(default)
            else:
                total += link.path_m
        return total

    def link_between(self, a, b) -> Link:
        for link in self.links:
            if (link.a, link.b) == (a, b) or (link.a, link.b) == (b, a):
                return link
        raise KeyError(f"{self.key}: no link between {a!r} and {b!r}")

    def unmeasured_links(self):
        """Physics-carrying segments still resting on a nominal number."""
        return tuple(l for l in self.links
                     if l.path_m is not None
                     and l.provenance == PROVENANCE_NOMINAL)


# ---------------------------------------------------------------------------
# Polyline geometry, used to resolve a click into a place on a beam
# ---------------------------------------------------------------------------


def segment_projection(point, start, end):
    """(closest point, distance, fraction along) for a point onto a segment."""
    px, py = point
    ax, ay = start
    bx, by = end
    dx, dy = bx - ax, by - ay
    span = dx * dx + dy * dy
    if span <= 0.0:
        return (ax, ay), math.hypot(px - ax, py - ay), 0.0
    t = ((px - ax) * dx + (py - ay) * dy) / span
    t = min(max(t, 0.0), 1.0)
    cx, cy = ax + t * dx, ay + t * dy
    return (cx, cy), math.hypot(px - cx, py - cy), t


def polyline_projection(point, points):
    """Closest approach of `point` to a polyline.

    Returns (closest point, distance, fraction of total display length). The
    fraction is what maps a click onto a physical distance along the segment.
    """
    if len(points) < 2:
        raise ValueError("a polyline needs at least two points")
    spans = [math.dist(points[i], points[i + 1]) for i in range(len(points) - 1)]
    total = sum(spans)
    best = None
    travelled = 0.0
    for index in range(len(points) - 1):
        closest, distance, t = segment_projection(point, points[index],
                                                  points[index + 1])
        along = travelled + t * spans[index]
        if best is None or distance < best[1]:
            best = (closest, distance, along)
        travelled += spans[index]
    closest, distance, along = best
    return closest, distance, (along / total if total > 0 else 0.0)


def polyline_length(points):
    return sum(math.dist(points[i], points[i + 1])
               for i in range(len(points) - 1))
