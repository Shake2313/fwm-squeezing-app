"""
Where things sit on the table, and how far the light travels between them.

Two coordinate systems that never mix: `display` for drawing, `path` for physics.
See `table.py` for why, and `parts.py` for the Sim et al. layouts.
"""
from .parts import LAYOUTS, PART1, PART2, cell_to_dmirror_m, cell_to_lens_m, get
from .probe import ProbeReading, beam_at, optic_at
from .render import beam_width, build_spec, legend_rows
from .table import (BEAM_COMMON, BEAM_CONJUGATE, BEAM_PROBE, BEAM_PUMP,
                    BEAM_REFERENCE, BEAM_SEED, Layout, Link, Node)

__all__ = [
    "Layout", "Link", "Node",
    "LAYOUTS", "PART1", "PART2", "get",
    "cell_to_dmirror_m", "cell_to_lens_m",
    "beam_at", "optic_at", "ProbeReading",
    "build_spec", "beam_width", "legend_rows",
    "BEAM_COMMON", "BEAM_PUMP", "BEAM_SEED", "BEAM_PROBE", "BEAM_CONJUGATE",
    "BEAM_REFERENCE",
]
