"""Dated apparatus inputs must affect the correct optical path and plane."""
from dataclasses import replace

import pytest

from gabes.schemes import fwm
from sabes import bridge
from sabes.beamline import SetupSettings, build_source_chain
from sabes.calibration import default_calibration
from sabes.detection import DetectionSettings, geometry, readout
from sabes.layout import beam_at, parts
from sabes.layout.probe import resolve_link_beam
from sabes.layout.table import BEAM_REFERENCE


def test_nel_is_a_connected_post_filter_transmission_only_element():
    assert "noise_eater" not in parts.ROUTE_SEED_TO_EOM
    assert "noise_eater" in parts.ROUTE_SEED_FILTERS
    parts.PART1.link_between("m_seed_out", "noise_eater")
    parts.PART1.link_between("noise_eater", "m_seed_down")

    calibration = default_calibration()
    baseline = build_source_chain(calibration=calibration)
    attenuated = build_source_chain(calibration=calibration.with_value(
        "t_noise_eater", calibration.value("t_noise_eater") * 0.5))
    assert attenuated.stage_power("seed", "into fiber EOM") == pytest.approx(
        baseline.stage_power("seed", "into fiber EOM"))
    assert attenuated.eom_output == baseline.eom_output
    assert attenuated.seed_power_w == pytest.approx(0.5 * baseline.seed_power_w)
    assert attenuated.seed_purity == pytest.approx(baseline.seed_purity)
    assert baseline.stage_power("seed", "noise eater") == pytest.approx(
        baseline.stage_power("seed", "before noise eater")
        * calibration.value("t_noise_eater"))


def test_sql_bypass_is_carrier_reference_and_not_generated_twin_light():
    chain = build_source_chain()
    links = [link for link in parts.PART2.links
             if link.carries.startswith("SQL")]
    assert links
    for link in links:
        assert link.beam == BEAM_REFERENCE
        assert link.path_m is None
        beam = resolve_link_beam(chain, link)
        assert beam is not None
        assert beam.total_power_w > 0.0
        assert all(line.offset_hz == 0.0 for line in beam.lines)
    probe, conjugate = chain.sql_reference_arms
    bypass_w = chain.stage_power("reference", "SQL bypass")
    assert probe.total_power_w + conjugate.total_power_w == pytest.approx(
        bypass_w)
    assert probe.total_power_w == pytest.approx(conjugate.total_power_w)
    assert sum(left.conjugate() * right for left, right in
               zip(probe.jones, conjugate.jones)) == pytest.approx(0.0)
    on_bypass = beam_at(parts.PART2, 400, 420, chain)
    assert on_bypass is not None
    assert on_bypass.beam_name == BEAM_REFERENCE
    assert on_bypass.power_w == pytest.approx(bypass_w)
    # Drawing the alternative route must not select SQL at a shared normal
    # collection stroke when no detector switching state has been supplied.
    on_collection = beam_at(parts.PART2, 850, 214, chain)
    assert on_collection is not None
    assert on_collection.beam_name == "probe"


def test_exit_alignment_triangle_is_separate_from_nominal_mode_propagation():
    settings = SetupSettings()
    chain = build_source_chain(settings)
    calibration = default_calibration()
    geom = geometry(chain, settings, calibration=calibration)
    assert geom.exit_to_dmirror_m == pytest.approx(0.5842)
    assert geom.mode_to_dmirror_m == pytest.approx(
        0.5842 + 0.5 * calibration.value("cell_length_mm") * 1e-3)
    assert geom.pump_separation_m == pytest.approx(
        settings.dmirror_separation_mm * 1e-3)
    assert geom.twin_radius_at_dmirror_m == pytest.approx(
        chain.seed.mode.propagated(geom.mode_to_dmirror_m).radius_m)
    assert "unmeasured" in geom.mode_origin


def test_recorded_cell_radii_override_telescope_estimates_only_when_requested():
    baseline = build_source_chain()
    settings = replace(SetupSettings(), pump_waist_override_um=530.0,
                       seed_waist_override_um=330.0)
    recorded = build_source_chain(settings)
    assert recorded.pump_waist_m == pytest.approx(530e-6)
    assert recorded.seed_waist_m == pytest.approx(330e-6)
    assert recorded.pump_power_w == pytest.approx(baseline.pump_power_w)
    assert recorded.seed_power_w == pytest.approx(baseline.seed_power_w)
    assert recorded.seed.mode.z_m == pytest.approx(baseline.seed.mode.z_m)


@pytest.mark.parametrize("value", (0.0, -1.0, float("nan"), float("inf")))
def test_recorded_cell_radii_reject_invalid_values(value):
    with pytest.raises(ValueError, match="finite and positive"):
        replace(SetupSettings(), seed_waist_override_um=value)


def test_aggregate_efficiency_inputs_apply_consistently_to_both_detector_arms():
    settings = SetupSettings()
    chain = build_source_chain(settings)
    detector = DetectionSettings(post_cell_transmission_override=0.9709,
                                 qe_pct_override=90.5)
    geom = geometry(chain, settings, detector)
    assert geom.optics_transmission == pytest.approx(0.9709)
    assert geom.qe_pct == pytest.approx(90.5)
    assert geom.eta == pytest.approx(0.9709 * 0.905)
    assert all(arm.optics_transmission == pytest.approx(0.9709)
               for arm in geom.arms)
    result = readout(geom, chain, (15.0, 14.0), settings, detector)
    expected_responsivity = 0.905 * 795.0 / 1239.841984
    assert result.shot_noise_a_per_rthz ** 2 == pytest.approx(
        2.0 * 1.602176634e-19 * expected_responsivity * result.total_power_w)


@pytest.mark.parametrize("kwargs", (
    {"post_cell_transmission_override": -0.1},
    {"post_cell_transmission_override": 1.1},
    {"qe_pct_override": 101.0},
    {"qe_pct_override": 0.0},
    {"qe_pct_override": float("nan")},
))
def test_aggregate_efficiency_inputs_reject_unphysical_values(kwargs):
    with pytest.raises(ValueError, match="must be"):
        DetectionSettings(**kwargs)


@pytest.mark.parametrize("leakage_dbm", (-60.0, -float("inf")))
def test_zero_transmission_has_no_shot_noise_or_detector_clearance(leakage_dbm):
    settings = SetupSettings()
    chain = build_source_chain(settings)
    detector = DetectionSettings(post_cell_transmission_override=0.0,
                                 qe_pct_override=90.5,
                                 pump_leakage_dbm=leakage_dbm)
    geom = geometry(chain, settings, detector)
    result = readout(geom, chain, (15.0, 14.0), settings, detector)
    assert result.total_power_w == 0.0
    assert result.shot_noise_a_per_rthz == 0.0
    assert result.electronic_noise_a_per_rthz > 0.0
    assert result.clearance_db == -float("inf")
    assert not result.measurable
    assert result.eom_noise.sql_psd_a2_per_hz == 0.0
    assert result.eom_noise.modes == ()
    assert all(arm.power_w == 0.0 for arm in result.arms)
    assert any("normalization is unavailable" in message
               for message in result.warnings)


def test_dark_detector_never_reports_normalized_noise_from_legacy_baseline():
    settings = SetupSettings()
    chain = build_source_chain(settings)
    detector = DetectionSettings(post_cell_transmission_override=0.0)
    geom = geometry(chain, settings, detector)
    actual_readout = readout(geom, chain, (15.0, 14.0), settings, detector)
    params = bridge.to_gabes_params(chain, settings, geom, quantize=False)
    center = fwm.branch_center_GHz(settings.opd_ghz, -1)
    raw = {
        "probe_axis_GHz": [center + settings.tpd_mhz * 1e-3],
        "gain_referred_noise_dB": [-8.0],
        "G_s": [15.0], "G_c": [14.0],
    }
    result = bridge.SabesResult(settings, chain, geom, params, raw=raw,
                               readout=actual_readout)
    assert result.gain_referred_noise_db == -8.0
    assert result.readout.eom_noise.sql_psd_a2_per_hz == 0.0
    assert result.eom_shot_noise_only_db is None
    assert result.eom_rin_loaded_noise_db is None
    assert result.eom_rin_penalty_db is None
    assert result.eom_coherent_mixture_noise_db is None
    assert result.eom_loaded_noise_db is None
    assert result.eom_noise_penalty_db is None


@pytest.mark.parametrize("arm,suffix", (("probe", "probe"),
                                       ("conjugate", "conj")))
def test_profiler_at_dmirror_matches_detector_geometry(arm, suffix):
    settings = replace(SetupSettings(), seed_waist_override_um=330.0)
    chain = build_source_chain(settings)
    calibration = default_calibration().with_value("cell_length_mm", 20.0)
    geom = geometry(chain, settings, calibration=calibration)
    mirror = parts.PART2.node(f"dm_{suffix}")
    reading = beam_at(parts.PART2, mirror.x, mirror.y, chain,
                      calibration=calibration)
    assert reading.beam_name == arm
    assert reading.radius_m == pytest.approx(geom.twin_radius_at_dmirror_m)
    assert "cell-centre" in reading.mode_scope


@pytest.mark.parametrize("arm", ("probe", "conjugate"))
def test_profiler_at_detector_matches_focused_spot_and_active_defocus(arm):
    settings = SetupSettings()
    chain = build_source_chain(settings)
    detector = DetectionSettings(probe_lens_focal_mm=140.0,
                                 conjugate_lens_focal_mm=120.0,
                                 pd_defocus_mm=7.0)
    calibration = default_calibration().with_value("cell_length_mm", 20.0)
    geom = geometry(chain, settings, detector, calibration)
    # Both eyes share one detector glyph. Suppress the other twin stroke to
    # inspect each arm at that exact endpoint without a geometric tie.
    other_arm = "conjugate" if arm == "probe" else "probe"
    arm_layout = replace(parts.PART2, links=tuple(
        link for link in parts.PART2.links if link.beam != other_arm))
    bpd = arm_layout.node("bpd")
    reading = beam_at(arm_layout, bpd.x, bpd.y, chain,
                      calibration=calibration, detection=detector)
    assert reading.beam_name == arm
    expected = next(value for value in geom.arms if value.name == arm)
    assert reading.radius_m == pytest.approx(expected.spot_radius_m)
    assert "post-lens waist" in reading.mode_scope
