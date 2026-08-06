"""
Table-drawing checks (Stage C of the optical-table plan).

    python tests/test_sabes_render.py   # or: pytest tests/test_sabes_render.py

The things worth pinning here are the ones where a drawing can quietly lie:
that stroke width encodes power rather than beam size, that "clickable" and
"modelled" cannot drift apart, and that a spec still crosses the iframe boundary
as plain JSON.
"""
import json
import sys
from dataclasses import replace
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from sabes import layout  # noqa: E402
from sabes.beamline import SetupSettings, build_source_chain  # noqa: E402
from sabes.layout import parts, render, table  # noqa: E402


@pytest.fixture(scope="module")
def chain():
    return build_source_chain()


def _shapes(spec, kind=None):
    return [s for s in spec["shapes"]
            if kind is None or s.get("kind") == kind]


# ------------------------------------------------------------ beam width

def test_beam_width_encodes_power_logarithmically():
    """Not beam size: a 530 um beam on a metre of table would be invisible."""
    decade = render.beam_width(1e-3) - render.beam_width(1e-4)
    assert decade == pytest.approx(render.WIDTH_PER_DECADE, abs=1e-9)
    assert render.beam_width(0.6) > render.beam_width(8e-6)


def test_beam_width_is_clamped_at_both_ends():
    assert render.beam_width(1e-15) == render.WIDTH_MIN
    assert render.beam_width(1e6) == render.WIDTH_MAX
    assert render.beam_width(0.0) == render.WIDTH_MIN     # no log domain error


def test_the_pump_is_drawn_heavier_than_the_seed(chain):
    spec = layout.build_spec(parts.PART2, chain)
    def width(a, b):
        link = parts.PART2.link_between(a, b)
        points = [[float(x), float(y)] for x, y in parts.PART2.points(link)]
        strokes = [s["width"] for s in _shapes(spec, "path")
                   if s["points"] == points]
        return max(strokes)
    assert width("gtp_pump", "pbs_in") > width("gtp_seed", "pbs_in")


# ------------------------------------------------------------ honesty

def test_clickable_and_modelled_cannot_drift_apart(chain):
    for part in layout.LAYOUTS:
        spec = layout.build_spec(part, chain)
        for shape in _shapes(spec):
            if shape.get("kind") == "path":
                continue
            assert shape["clickable"] == part.node(shape["id"]).clickable


def test_decorative_hardware_recedes(chain):
    spec = layout.build_spec(parts.PART1, chain)
    by_id = {s["id"]: s for s in _shapes(spec) if s.get("id")}
    assert by_id["noise_eater"]["opacity"] < 1.0
    assert "not modelled" in by_id["noise_eater"]["title"]
    assert "opacity" not in by_id["eom"] or by_id["eom"]["opacity"] == 1.0


def test_beams_are_never_clickable(chain):
    """Clicking a beam means clicking the table under it -- that is how an
    instrument head gets placed on the beam."""
    for part in layout.LAYOUTS:
        for shape in _shapes(layout.build_spec(part, chain), "path"):
            assert shape["clickable"] is False


def test_segments_with_no_known_beam_are_drawn_faint(chain):
    """The locking arm carries light the source chain does not model."""
    spec = layout.build_spec(parts.PART1, chain)
    reference = [s for s in _shapes(spec, "path")
                 if s["stroke"] == render.BEAM_STYLE[table.BEAM_REFERENCE]["stroke"]]
    assert reference
    assert all(s["opacity"] < 1.0 for s in reference)
    assert all(s["width"] in (render.WIDTH_MIN, render.WIDTH_MIN + 2.5)
               for s in reference)


# ------------------------------------------------------------ glyphs

def test_every_node_gets_a_glyph_and_unknown_ones_fall_back():
    seen = set()
    for part in layout.LAYOUTS:
        for node in part.nodes:
            glyph = render.glyph_for(node)
            assert glyph
            seen.add(glyph)
    for expected in ("cube", "waveplate", "mirror", "etalon", "lens", "cell",
                     "collimator", "detector", "iris", "polarizer", "dmirror"):
        assert expected in seen, expected
    stranger = table.Node("nothing_matches", "?", 0, 0)
    assert render.glyph_for(stranger) == "box"


def test_glyph_parts_stay_inside_the_clickable_group(chain):
    """Sub-primitives make a polarizing cube look like one without splitting it
    into separate click targets."""
    spec = layout.build_spec(parts.PART2, chain)
    by_id = {s["id"]: s for s in _shapes(spec) if s.get("id")}
    assert by_id["pbs_in"]["parts"], "a PBS needs its diagonal"
    assert by_id["cell"]["parts"] == []
    for shape in _shapes(spec):
        for part in shape.get("parts", []):
            assert part["kind"] in ("rect", "circle", "line")


def test_selection_is_visible_in_the_spec(chain):
    plain = layout.build_spec(parts.PART2, chain)
    picked = layout.build_spec(parts.PART2, chain, selected="cell")
    def cell(spec):
        return next(s for s in _shapes(spec) if s.get("id") == "cell")
    assert cell(picked)["stroke"] != cell(plain)["stroke"]
    assert cell(picked)["stroke_width"] > cell(plain)["stroke_width"]


# ------------------------------------------------------------ transport

def test_spec_survives_the_json_round_trip(chain):
    """It crosses an iframe boundary, so a numpy scalar would break rendering."""
    for part in layout.LAYOUTS:
        spec = layout.build_spec(part, chain, probe=(10, 20))
        restored = json.loads(json.dumps(spec))
        assert restored["width"] == part.width
        assert restored["probe"] == [10.0, 20.0]
        assert len(restored["shapes"]) == len(spec["shapes"])


def test_spec_stays_small_enough_to_ship_every_rerun(chain):
    for part in layout.LAYOUTS:
        payload = json.dumps(layout.build_spec(part, chain))
        assert len(payload) < 60_000, f"{part.key}: {len(payload)} bytes"


def test_grid_is_declared_but_is_not_a_scale(chain):
    """The breadboard dots are atmosphere: display units are not millimetres."""
    spec = layout.build_spec(parts.PART2, chain)
    assert spec["grid"] > 0
    ratios = list(layout.probe.display_scale_m_per_unit(parts.PART2).values())
    assert max(ratios) / min(ratios) > 10.0


# ------------------------------------------------------------ live values

def test_widths_follow_the_chain_rather_than_being_baked_in(chain):
    """Turn the pump down and the drawing has to thin out."""
    quiet = build_source_chain(replace(SetupSettings(), ta_current_a=0.9))
    loud = layout.build_spec(parts.PART2, chain)
    soft = layout.build_spec(parts.PART2, quiet)
    def pump_width(spec):
        link = parts.PART2.link_between("cell", "pbs_out")
        points = [[float(x), float(y)] for x, y in parts.PART2.points(link)]
        return max(s["width"] for s in _shapes(spec, "path")
                   if s["points"] == points)
    assert pump_width(soft) < pump_width(loud)


def test_twin_states_light_up_the_post_cell_arms(chain):
    seed = chain.seed
    twin = {"probe": seed.scaled(14.0), "conjugate": seed.scaled(13.0)}
    plain = layout.build_spec(parts.PART2, chain)
    gained = layout.build_spec(parts.PART2, chain, twin_states=twin)
    def probe_width(spec):
        link = parts.PART2.link_between("dm_probe", "iris_probe")
        points = [[float(x), float(y)] for x, y in parts.PART2.points(link)]
        return max(s["width"] for s in _shapes(spec, "path")
                   if s["points"] == points)
    assert probe_width(gained) > probe_width(plain)


def test_legend_covers_every_beam_colour_actually_drawn(chain):
    drawn = set()
    for part in layout.LAYOUTS:
        drawn.update(s["stroke"] for s in _shapes(layout.build_spec(part, chain),
                                                  "path"))
    assert drawn <= {colour for colour, _dash, _name in layout.legend_rows()}


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-q"]))
