"""
Resolve a point on the drawing into a real beam.

`beam_at(layout, x, y)` is the whole reason the layout exists: it is what turns a
click into physics, and it is the single input every virtual instrument takes. A
power meter, a beam profiler and a wavemeter differ only in what they do with its
result; an oscilloscope and a spectrum analyser take a photodiode's output, and
the photodiode takes this.

The mode is propagated to the probe point rather than reported at the segment's
start, so a beam profiler dropped halfway down a long arm reads the size the beam
actually has there. That is the difference between a diagram and a simulator.
"""
from dataclasses import dataclass
from typing import Optional

from ..beamstate import BeamState
from .table import Layout, Link, polyline_length, polyline_projection

#: How close a click must land to a beam, in display units, to count as being on
#: it. Generous on purpose: the drawn stroke is a visual encoding rather than the
#: beam's true width, so demanding pixel accuracy would just be annoying.
PICK_TOLERANCE = 14.0


@dataclass(frozen=True)
class ProbeReading:
    """What sits at a probed point on the table."""
    link: Link
    beam: BeamState
    distance_along_m: float
    fraction_along: float
    display_point: tuple
    pick_distance: float

    @property
    def beam_name(self):
        return self.link.beam

    @property
    def radius_m(self):
        return self.beam.mode.radius_m

    @property
    def power_w(self):
        return self.beam.total_power_w


#: A link's beam maps onto the source chain's path. Both arms have a stage named
#: "at cell", so keying stages by name alone silently hands back the wrong beam.
_BEAM_TO_CHAIN_PATH = {
    "pump": "pump",
    "seed": "seed",
    "common": "common",
}


def beam_states_by_stage(chain):
    """{(path, stage name): BeamState} for the source chain."""
    return {(stage.path, stage.name): stage.beam for stage in chain.stages}


def resolve_link_beam(chain, link, twin_states=None):
    """The BeamState travelling a link, or None if nothing is known about it.

    `twin_states` lets the caller supply probe/conjugate states once the solve
    has produced gains; without it the twins fall back to the seed as it leaves
    the cell, which is the right *mode* but not the right power.
    """
    if link.beam in ("probe", "conjugate"):
        if twin_states and link.beam in twin_states:
            return twin_states[link.beam]
        return chain.seed
    if not link.carries:
        return None
    path = _BEAM_TO_CHAIN_PATH.get(link.beam)
    if path is None:
        return None
    states = beam_states_by_stage(chain)
    state = states.get((path, link.carries))
    if state is None:
        # Common-path stages are shared, so an arm may legitimately point at one.
        state = states.get(("common", link.carries))
    return state


def beam_at(layout: Layout, x, y, chain, *, twin_states=None,
            tolerance=PICK_TOLERANCE) -> Optional[ProbeReading]:
    """Nearest beam to a display point, with its mode propagated to that point.

    Returns None when the click is not on any beam -- which is a legitimate
    answer, not an error: most of an optical table is empty.
    """
    best = None
    for link in layout.links:
        points = layout.points(link)
        closest, distance, fraction = polyline_projection((x, y), points)
        if distance > tolerance:
            continue
        if best is None or distance < best[1]:
            best = (link, distance, closest, fraction)
    if best is None:
        return None

    link, distance, closest, fraction = best
    state = resolve_link_beam(chain, link, twin_states)
    if state is None:
        return None

    # Display fraction maps onto physical distance only when the segment has a
    # real length. Without one there is nothing to propagate along, and the
    # honest answer is the beam as it entered the segment.
    along_m = 0.0 if link.path_m is None else link.path_m * fraction
    if along_m:
        state = state.with_mode(state.mode.propagated(along_m))

    return ProbeReading(link=link, beam=state, distance_along_m=along_m,
                        fraction_along=fraction, display_point=closest,
                        pick_distance=distance)


def optic_at(layout: Layout, x, y):
    """Node whose drawn body contains a point, preferring clickable ones."""
    hits = []
    for node in layout.nodes:
        half_w = max(node.w, 10.0) / 2.0
        half_h = max(node.h, 10.0) / 2.0
        if abs(x - node.x) <= half_w and abs(y - node.y) <= half_h:
            hits.append(node)
    if not hits:
        return None
    clickable = [n for n in hits if n.clickable]
    return (clickable or hits)[0]


def display_scale_m_per_unit(layout: Layout):
    """Metres of optical path per display unit, per link that has both.

    Diagnostic, not a conversion: a schematic has no single scale, and the spread
    this reports is exactly why display coordinates are never used as distances.
    """
    ratios = {}
    for link in layout.links:
        if link.path_m is None:
            continue
        span = polyline_length(layout.points(link))
        if span > 0:
            ratios[(link.a, link.b)] = link.path_m / span
    return ratios
