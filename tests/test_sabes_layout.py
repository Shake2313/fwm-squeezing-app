"""
Layout-model checks (Stage A of the optical-table plan).

    python tests/test_sabes_layout.py   # or: pytest tests/test_sabes_layout.py

The load-bearing property is the separation of the two coordinate systems:
`display` coordinates draw the table and `path_m` lengths do the physics, and
nothing is allowed to derive one from the other. The setup drawings are
schematics, so a pixel distance is not a metre distance, and letting them mix
would corrupt every beam size downstream without anything looking wrong.
"""
import sys
from dataclasses import replace
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from sabes import detection, layout  # noqa: E402
from sabes.beamline import SetupSettings, build_source_chain  # noqa: E402
from sabes.calibration import default_calibration  # noqa: E402
from sabes.detection import DetectionSettings  # noqa: E402
from sabes.layout import parts, table  # noqa: E402


@pytest.fixture(scope="module")
def chain():
    return build_source_chain()


# ------------------------------------------------------- structural integrity

def test_both_layouts_are_well_formed():
    assert {l.key for l in layout.LAYOUTS} == {"part1", "part2"}
    for part in layout.LAYOUTS:
        assert part.nodes and part.links
        for link in part.links:
            part.node(link.a)      # raises if a link dangles
            part.node(link.b)


def test_duplicate_node_ids_are_rejected():
    node = table.Node("dup", "x", 0, 0)
    with pytest.raises(ValueError, match="duplicate node ids"):
        table.Layout("bad", "t", 10, 10, nodes=(node, node))


def test_links_to_unknown_nodes_are_rejected():
    with pytest.raises(ValueError, match="unknown node"):
        table.Layout("bad", "t", 10, 10,
                     nodes=(table.Node("a", "", 0, 0),),
                     links=(table.Link("a", "ghost"),))


def test_negative_and_mistyped_path_lengths_are_rejected():
    with pytest.raises(ValueError, match="negative path length"):
        table.Link("a", "b", path_m=-1.0)
    with pytest.raises(ValueError, match="bad provenance"):
        table.Link("a", "b", path_m=1.0, provenance="vibes")


# --------------------------------------------- the two coordinate systems

def test_display_coordinates_do_not_determine_path_lengths():
    """A schematic has no single scale, and this proves it does not need one.

    If display distance and optical path were the same thing, every link would
    share one metres-per-unit ratio. They do not, by more than an order of
    magnitude, which is exactly why `path_m` is stored rather than measured off
    the drawing.
    """
    from sabes.layout.probe import display_scale_m_per_unit

    ratios = list(display_scale_m_per_unit(parts.PART2).values())
    assert len(ratios) >= 6
    assert max(ratios) / min(ratios) > 10.0


def test_path_lengths_carry_provenance_and_the_nominal_ones_are_visible():
    physics_links = [l for l in parts.PART2.links if l.path_m is not None]
    assert physics_links
    assert {l.provenance for l in physics_links} <= {
        table.PROVENANCE_DESIGN, table.PROVENANCE_MEASURED,
        table.PROVENANCE_NOMINAL}
    # The post-cell path is still partly guessed, and the model must say so.
    unmeasured = {(l.a, l.b) for l in parts.PART2.unmeasured_links()}
    assert ("dm_probe", "iris_probe") in unmeasured
    assert ("iris_probe", "lens_probe") in unmeasured


def test_telescope_separations_are_design_values_not_guesses():
    """An afocal pair has its lenses f1 + f2 apart by definition."""
    for a, b, expected in (("seed_l1", "seed_l2", 0.400 + 0.150),
                           ("pump_l1", "pump_l2", 0.300 + 0.150)):
        link = parts.PART2.link_between(a, b)
        assert link.path_m == pytest.approx(expected)
        assert link.provenance == table.PROVENANCE_DESIGN


def test_missing_path_length_raises_instead_of_being_treated_as_zero():
    with pytest.raises(ValueError, match="no path length"):
        parts.PART2.path_length_m("m_seed_a", "gtp_seed")
    assert parts.PART2.path_length_m("m_seed_a", "gtp_seed", default=0.0) == 0.0


# ------------------------------------------------- distances drive the physics

def test_layout_reproduces_the_distances_the_detection_model_used():
    assert parts.cell_to_dmirror_m() == pytest.approx(0.5842, abs=1e-9)
    assert parts.cell_to_lens_m() == pytest.approx(0.9000, abs=1e-9)


def test_detection_geometry_reads_distances_from_the_layout(chain):
    settings = SetupSettings()
    base = detection.geometry(chain, settings)
    moved = detection.geometry(
        chain, settings,
        layout=parts.with_detection_distances(cell_to_dmirror=1.4))
    assert moved.twin_radius_at_dmirror_m > base.twin_radius_at_dmirror_m
    assert moved.pump_separation_m > base.pump_separation_m


def test_moving_the_focusing_lens_changes_the_photodiode_spot(chain):
    settings = SetupSettings()
    base = detection.geometry(chain, settings)
    further = detection.geometry(
        chain, settings, layout=parts.with_detection_distances(iris_to_lens=0.9))
    assert further.arms[0].spot_radius_m < base.arms[0].spot_radius_m


def test_crossing_angle_uses_the_measured_table_distance():
    assert parts.CELL_TO_DMIRROR_M == pytest.approx(0.5842)
    assert SetupSettings().crossing_angle_deg == pytest.approx(0.32, abs=1e-4)


# ------------------------------------------------------------- beam_at

def test_beam_at_returns_the_right_arm_not_just_the_right_stage_name(chain):
    """Both arms have a stage called 'at cell'; keying on the name alone would
    hand the pump's segment the seed's beam."""
    on_pump = layout.beam_at(parts.PART2, 300, 256, chain)
    assert on_pump is not None and on_pump.beam_name == layout.BEAM_PUMP
    assert on_pump.power_w > 0.1                      # hundreds of mW
    on_seed = layout.beam_at(parts.PART2, 500, 60, chain)
    assert on_seed is not None and on_seed.beam_name == layout.BEAM_SEED
    assert on_seed.power_w < 0.01


def test_beam_at_propagates_the_mode_to_the_probed_point(chain):
    """A profiler halfway along an arm must read the size the beam has there."""
    link = parts.PART2.link_between("pump_l1", "pump_l2")
    points = parts.PART2.points(link)
    near = layout.beam_at(parts.PART2, *_lerp(points[0], points[-1], 0.05), chain)
    far = layout.beam_at(parts.PART2, *_lerp(points[0], points[-1], 0.95), chain)
    assert near is not None and far is not None
    assert far.distance_along_m > near.distance_along_m
    assert near.radius_m != far.radius_m


def test_beam_at_returns_nothing_on_empty_table(chain):
    """Most of a table has no beam on it, and that is an answer, not a failure."""
    assert layout.beam_at(parts.PART2, 620, 380, chain) is None


def test_beam_at_respects_the_pick_tolerance(chain):
    # Take a point that is genuinely on the probe arm, so the test keeps meaning
    # something if the drawing is nudged in Stage C.
    points = parts.PART2.points(parts.PART2.link_between("pbs_out", "dm_probe"))
    on_beam_x, on_beam_y = _lerp(points[0], points[-1], 0.5)
    hit = layout.beam_at(parts.PART2, on_beam_x, on_beam_y, chain, tolerance=6.0)
    assert hit is not None and hit.link.b == "dm_probe"
    assert hit.pick_distance <= 6.0

    # Step off by more than the tolerance and this link is no longer the answer.
    # It may not be *nothing*: the conjugate arm runs parallel below, and picking
    # up the neighbouring beam is correct behaviour on a table this dense.
    stepped = layout.beam_at(parts.PART2, on_beam_x, on_beam_y + 40, chain,
                             tolerance=6.0)
    assert stepped is None or stepped.link.b != "dm_probe"


def test_segments_without_a_length_do_not_fake_propagation(chain):
    """No length means nothing to propagate along; report the entry state."""
    reading = layout.beam_at(parts.PART2, 250, 200, chain)
    if reading is not None and reading.link.path_m is None:
        assert reading.distance_along_m == 0.0


# ---------------------------------------------------------- clickability

def test_clickable_means_there_is_something_to_edit():
    """Anything with a physics knob OR a transmission opens; nothing else does."""
    for part in layout.LAYOUTS:
        for node in part.nodes:
            assert node.clickable == bool(node.params or node.transmission_key)
    # A steering mirror after the amplifier owns its transmission and opens...
    assert parts.PART2.node("m_pump_a").clickable
    # ...while one before it has nothing to edit.
    assert not parts.PART1.node("m_ecdl_1").clickable


def test_clickable_optics_name_real_settings_fields():
    from dataclasses import fields
    known = {f.name for f in fields(SetupSettings)}
    known |= {f.name for f in fields(DetectionSettings)}
    known |= {"eom_offset_mhz"}
    known |= {f"etalon_detune_ghz_{n}" for n in (1, 2, 3)}
    for part in layout.LAYOUTS:
        for node in part.clickable_nodes():
            for param in node.params:
                assert param in known, f"{part.key}/{node.id}: {param}"


def test_every_settings_knob_has_an_optic_that_owns_it():
    """Stage D moves the sidebar into popups, so nothing may be orphaned."""
    owned = set()
    for part in layout.LAYOUTS:
        for node in part.nodes:
            owned.update(node.params)
    from dataclasses import fields
    knobs = {f.name for f in fields(SetupSettings)}
    knobs |= {f.name for f in fields(DetectionSettings)}
    # Recorded-only annotations and internal numerics are deliberately unowned.
    knobs -= {"ecdl_temp_c", "ecdl_current_ma", "ecdl_pzt_v", "ta_temp_c",
              "eom_frequency_hz", "bpd_gain_v_per_a", "analysis_bandwidth_hz"}
    # Exposed per stage rather than as one tuple, so the owner check is on the
    # per-stage names instead.
    knobs.discard("etalon_detune_ghz")
    assert {f"etalon_detune_ghz_{n}" for n in (1, 2, 3)} <= owned
    missing = knobs - owned
    assert not missing, f"no optic owns: {sorted(missing)}"


# ------------------------------------------- part 1 corrections from the drawing

def test_part1_has_the_glan_taylor_after_the_etalon_chain():
    """The drawing has several Glan-Taylors; the linear model had one."""
    parts.PART1.node("gt_seed")
    assert parts.PART1.link_between("etalon_3", "gt_seed")
    assert parts.PART1.link_between("gt_seed", "m_seed_out")
    # And part 2 still has one on each arm ahead of the cell PBS.
    parts.PART2.node("gtp_seed")
    parts.PART2.node("gtp_pump")


def test_part1_has_the_half_wave_plate_between_etalons_one_and_two():
    parts.PART1.node("hwp_etalon")
    assert parts.PART1.link_between("m_et_1", "hwp_etalon")
    assert parts.PART1.link_between("hwp_etalon", "m_et_2")


def test_the_two_collimators_are_distinct_parts():
    assert parts.PART1.node("f230_pump").label == "F230"
    assert parts.PART1.node("f220_seed").label == "F220"


def test_unmodelled_hardware_is_drawn_and_labelled_lumped():
    """Drawn so the table looks like the table, and labelled so nothing is
    implied about how much of it the physics actually knows."""
    for node_id in ("oi_ecdl", "oi_mopa", "cl", "noise_eater", "rb_ref"):
        node = parts.PART1.node(node_id)
        assert node.lumped
        if node.label:
            assert node.display_label.endswith("(lumped)")
    # A detail-modelled optic must not be labelled that way.
    assert not parts.PART2.node("cell").display_label.endswith("(lumped)")
    assert "locking" in parts.PART1.node("ecdl").note.lower() or \
        "reference arm" in parts.PART1.node("ecdl").note.lower()


def _lerp(a, b, t):
    return (a[0] + (b[0] - a[0]) * t, a[1] + (b[1] - a[1]) * t)


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-q"]))
