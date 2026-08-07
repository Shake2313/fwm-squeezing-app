"""
SABES source-chain checks.

    python tests/test_sabes_optics.py   # or: pytest tests/test_sabes_optics.py

Three things are pinned here:
  1. device laws against closed-form values that do not depend on calibration;
  2. the chain reproducing the Sim et al. operating point (600 mW / 8 uW /
     530 / 330 um / -8 MHz TPD / << 0.5 % carrier);
  3. the honesty machinery -- provenance tags and the etalon leakage floor,
     because a silently optimistic filter model would inflate every squeezing
     number downstream.
"""
import sys
from dataclasses import replace
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from gabes import constants  # noqa: E402
from sabes import beamstate, calibration as calib, devices  # noqa: E402
from sabes.beamline import (SetupSettings, build_source_chain,  # noqa: E402
                            nominal_settings, solve_seed_polarizer_deg,
                            solve_seed_trim_deg, solve_split_angle_deg)

WAVELENGTH = constants.WAVELENGTH_D1_85RB


# ---------------------------------------------------------------- primitives

def test_bessel_matches_known_values():
    assert devices.bessel_j(0, 1.0) == pytest.approx(0.7651976866, abs=1e-9)
    assert devices.bessel_j(1, 1.0) == pytest.approx(0.4400505857, abs=1e-9)
    assert devices.bessel_j(2, 2.0) == pytest.approx(0.3528340286, abs=1e-9)
    assert devices.bessel_j(-1, 2.0) == pytest.approx(-devices.bessel_j(1, 2.0))
    assert devices.bessel_j(0, 0.0) == 1.0
    assert devices.bessel_j(3, 0.0) == 0.0


def test_phase_modulation_conserves_power_across_the_ladder():
    """Sum of J_n^2 is 1, so a lossless modulator moves power but never adds it."""
    eom = devices.PhaseModulator(frequency_hz=3.036e9, v_pi_v=4.4,
                                 rf_power_dbm=18.0, max_order=12)
    assert sum(eom.sideband_fractions().values()) == pytest.approx(1.0, abs=1e-9)

    beam = beamstate.monochromatic(1.0, WAVELENGTH, 500e-6)
    assert eom.apply(beam).total_power_w == pytest.approx(1.0, rel=1e-9)


def test_rf_drive_reaches_the_j1_maximum_at_18_dbm():
    """18 dBm is not arbitrary: it puts beta near the J1 peak at 1.841."""
    eom = devices.PhaseModulator(frequency_hz=3.036e9, v_pi_v=4.4,
                                 rf_power_dbm=18.0)
    assert devices.rf_peak_volts(18.0) == pytest.approx(2.5119, abs=1e-3)
    assert eom.modulation_index == pytest.approx(1.7935, abs=1e-3)
    assert eom.sideband_fractions()[-1] == pytest.approx(0.338, abs=2e-3)


def test_gaussian_free_space_and_lens_are_self_consistent():
    mode = beamstate.GaussianMode(waist_m=330e-6, wavelength_m=WAVELENGTH)
    assert mode.rayleigh_m == pytest.approx(0.4300, abs=2e-3)
    assert mode.propagated(0.5842).radius_m == pytest.approx(556e-6, abs=2e-6)
    # A lens one focal length after a waist re-collimates to a new waist at f.
    relayed = mode.propagated(0.1).through_thin_lens(0.1).propagated(0.1)
    assert abs(relayed.z_m) < 1e-9


def test_m2_raises_divergence_without_moving_the_current_radius():
    mode = beamstate.GaussianMode(waist_m=500e-6, wavelength_m=WAVELENGTH,
                                  z_m=0.2)
    degraded = mode.with_m2(2.0)
    assert degraded.radius_m == pytest.approx(mode.radius_m, rel=1e-12)
    assert degraded.divergence_rad == pytest.approx(2.0 * mode.divergence_rad,
                                                    rel=1e-9)


def test_telescope_magnification_is_the_focal_ratio():
    """This is what reproduces the paper waists from the setup drawing."""
    for diameter_mm, f1, f2, expected_um in ((2.0, 300.0, 150.0, 500.0),
                                             (1.6, 400.0, 150.0, 300.0)):
        beam = beamstate.monochromatic(1.0, WAVELENGTH, 1e-3)
        beam = devices.FiberCollimator(diameter_mm * 1e-3).apply(beam)
        beam = devices.Telescope(f1 * 1e-3, f2 * 1e-3).apply(beam)
        assert beam.mode.waist_m * 1e6 == pytest.approx(expected_um, rel=1e-6)
        assert abs(beam.mode.z_m) < 1e-9


def test_half_wave_plate_and_polarizer_follow_malus():
    beam = beamstate.monochromatic(1.0, WAVELENGTH, 500e-6)
    for angle in (0.0, 10.0, 22.5, 45.0):
        rotated = devices.HalfWavePlate(angle).apply(beam)
        assert rotated.total_power_w == pytest.approx(1.0, rel=1e-12)
        passed = devices.Polarizer(0.0, extinction_ratio=1e12).apply(rotated)
        expected = pytest.approx(
            devices.math.cos(devices.math.radians(2 * angle)) ** 2, abs=1e-9)
        assert passed.total_power_w == expected


def test_pbs_ports_are_complementary_and_leak_by_the_extinction():
    beam = beamstate.monochromatic(1.0, WAVELENGTH, 500e-6)
    beam = devices.HalfWavePlate(45.0).apply(beam)      # rotate H -> V
    transmitted, reflected = devices.PolarizingBeamSplitter(3000.0).split(beam)
    assert reflected.total_power_w == pytest.approx(1.0, abs=1e-3)
    assert transmitted.total_power_w == pytest.approx(1.0 / 3000.0, rel=1e-6)


def test_amplifier_saturates_rather_than_scaling_linearly():
    amp = devices.TaperedAmplifier(current_a=2.5, threshold_a=0.5,
                                   slope_w_per_a=1.05, saturation_w=0.020)
    assert amp.max_power_w == pytest.approx(2.1)
    doubled_in = amp.gain_output_w(0.080) / amp.gain_output_w(0.040)
    assert 1.0 < doubled_in < 2.0        # saturated: less than proportional
    assert amp.gain_output_w(10.0) == pytest.approx(2.1, rel=1e-2)


def test_fiber_cleans_beam_quality_at_a_power_cost():
    beam = beamstate.monochromatic(1.0, WAVELENGTH, 500e-6)
    beam = beam.with_mode(beam.mode.with_m2(2.0))
    out = devices.FiberCoupling(1.0).apply(beam)
    assert out.mode.m2 == 1.0
    assert out.total_power_w == pytest.approx(devices.FiberCoupling.mode_penalty(2.0))
    assert devices.FiberCoupling.mode_penalty(1.0) == pytest.approx(1.0)


def test_aperture_clipping_is_monotone_and_bounded():
    mode = beamstate.GaussianMode(waist_m=500e-6, wavelength_m=WAVELENGTH)
    assert mode.clipping_transmission(1e9) == pytest.approx(1.0, abs=1e-12)
    assert mode.clipping_transmission(0.0) == pytest.approx(0.0, abs=1e-12)
    # A beam pushed well off the aperture is rejected.
    assert mode.clipping_transmission(1e-3, offset_m=1e-2) < 1e-12


# ---------------------------------------------------------------- etalon model

def test_etalon_peak_is_the_stated_transmission_and_stopband_is_the_floor():
    et = devices.EtalonFilter(fsr_hz=17e9, fwhm_hz=950e6,
                              peak_transmission=0.80, leak=0.0765)
    assert et.finesse == pytest.approx(17.895, abs=1e-3)
    assert et.transmission(0.0) == pytest.approx(0.80, rel=1e-12)
    assert et.transmission(17e9) == pytest.approx(0.80, rel=1e-12)   # next order
    # Mid-stopband sits just above the floor: the Airy tail adds ~1/F of the peak.
    assert et.transmission(8.5e9) == pytest.approx(0.0820, abs=5e-4)
    assert et.transmission(8.5e9) > et.leak


def test_etalon_leak_floor_reproduces_both_published_suppressions():
    """An ideal Airy cascade is far too optimistic; the floor is why.

    Sim et al. report ~0.5 % carrier:seed with two etalons and much less with
    three. Fitting the floor to the two-stage number must then predict the
    three-stage number, or the model is not describing the same filter.
    """
    eom = devices.PhaseModulator(frequency_hz=constants.NU_GROUND_HF,
                                 v_pi_v=4.4, rf_power_dbm=18.0)
    fractions = eom.sideband_fractions()
    ratio_in = fractions[0] / fractions[-1]

    bare = devices.EtalonFilter(fsr_hz=17e9, fwhm_hz=950e6,
                                peak_transmission=0.80)
    leak = devices.leak_from_carrier_ratio(
        bare, 2, constants.NU_GROUND_HF, ratio_in, 0.005)
    assert leak == pytest.approx(
        calib.default_calibration().value("etalon_leak"), rel=1e-3)

    tuned = replace(bare, leak=leak)
    def ratio(stages):
        return ratio_in * (tuned.transmission(constants.NU_GROUND_HF)
                           / tuned.transmission(0.0)) ** stages

    assert ratio(2) == pytest.approx(0.005, rel=1e-6)     # the fitted point
    assert ratio(3) < 0.005 / 5.0                         # "much less", predicted
    # The ideal cascade would claim an order of magnitude better than observed.
    assert ratio_in * (bare.transmission(constants.NU_GROUND_HF)
                       / bare.transmission(0.0)) ** 2 < 0.005 / 10.0


def test_etalon_passband_offset_cannot_be_used_to_solve_a_floor():
    bare = devices.EtalonFilter(fsr_hz=17e9, fwhm_hz=950e6)
    with pytest.raises(ValueError):
        devices.leak_from_carrier_ratio(bare, 2, 0.0, 0.35, 0.005)


# ---------------------------------------------------------------- provenance

def test_every_shipped_coefficient_declares_a_known_provenance():
    cal = calib.default_calibration()
    assert len(cal) > 20
    for coefficient in cal:
        assert coefficient.provenance in calib.PROVENANCE_ORDER


def test_unverified_coefficients_are_reported_for_the_ui_badge():
    cal = calib.default_calibration()
    unverified = {c.name for c in cal.unverified()}
    # These are placeholders and must not masquerade as measurements.
    assert "eom_v_pi_v" in unverified
    assert "ta_slope_w_per_a" in unverified
    # These are read off a datasheet or the paper and must not be flagged.
    assert "pd_responsivity_a_per_w" not in unverified
    assert "etalon_fsr_hz" not in unverified


def test_anchoring_promotes_a_coefficient_to_fitted():
    cal = calib.default_calibration()
    assert cal.get("eom_v_pi_v").provenance == "nominal"
    anchored = cal.anchored_to("eom_v_pi_v", 4.1, note="measured 2026-08-04")
    assert anchored.get("eom_v_pi_v").value == 4.1
    assert anchored.get("eom_v_pi_v").provenance == "fitted"
    assert anchored.get("eom_v_pi_v").verified
    # The original is untouched, so a session can branch safely.
    assert cal.get("eom_v_pi_v").value == 4.4


def test_calibration_round_trips_through_json(tmp_path):
    cal = calib.default_calibration()
    path = tmp_path / "roundtrip.json"
    cal.save(path)
    again = calib.Calibration.load(path)
    assert again.names() == cal.names()
    for name in cal.names():
        assert again.get(name) == cal.get(name)


def test_unknown_provenance_is_rejected():
    with pytest.raises(ValueError):
        calib.Coefficient("x", 1.0, provenance="probably-fine")


def test_missing_coefficient_error_points_at_the_calibration_file():
    with pytest.raises(KeyError, match="calibration file"):
        calib.default_calibration().get("no_such_coefficient")


# ---------------------------------------------------------------- full chain

def test_default_chain_reproduces_the_sim_operating_point():
    chain = build_source_chain()
    assert chain.pump_power_w == pytest.approx(0.600, rel=2e-3)
    assert chain.seed_power_w == pytest.approx(8.0e-6, rel=2e-3)
    # Waists come from collimator diameter x telescope magnification x one
    # shared scale factor, so they land near but not exactly on 530 / 330 um.
    assert chain.pump_waist_m * 1e6 == pytest.approx(530.0, rel=0.03)
    assert chain.seed_waist_m * 1e6 == pytest.approx(330.0, rel=0.03)
    assert chain.carrier_ratio < 0.005          # paper: << 0.5 % with 3 etalons
    assert chain.seed_purity > 0.99
    assert chain.warnings == ()


def test_modulator_frequency_sets_the_two_photon_detuning_exactly():
    """No calibration coefficient stands between the dial and delta."""
    assert SetupSettings().tpd_mhz == pytest.approx(-8.0, abs=1e-6)
    resonant = replace(SetupSettings(), eom_frequency_hz=constants.NU_GROUND_HF)
    assert resonant.tpd_mhz == pytest.approx(0.0, abs=1e-9)
    # The "3.036 GHz" shorthand is NOT the -8 MHz operating point.
    shorthand = replace(SetupSettings(), eom_frequency_hz=3.036e9)
    assert shorthand.tpd_mhz == pytest.approx(-0.2676, abs=1e-3)


def test_crossing_angle_follows_the_dmirror_separation():
    assert SetupSettings().crossing_angle_deg == pytest.approx(0.32, abs=2e-3)
    wider = replace(SetupSettings(), dmirror_separation_mm=6.52)
    assert wider.crossing_angle_deg == pytest.approx(0.64, abs=5e-3)


def test_seed_arm_carries_far_more_power_than_the_fwm_wants():
    """The trim throws most of it away -- so the seed budget is not a constraint."""
    chain = build_source_chain()
    before_trim = chain.stage_power("seed", "delivery fiber + collimator")
    assert before_trim > 10.0 * chain.seed_power_w


def test_split_and_polarizer_knobs_solve_independently():
    settings = nominal_settings()
    chain = build_source_chain(settings)
    assert chain.pump_power_w == pytest.approx(0.600, rel=1e-4)
    assert chain.stage_power("seed", "into fiber EOM") == pytest.approx(7.0e-3,
                                                                       rel=1e-4)
    assert chain.seed_power_w == pytest.approx(8.0e-6, rel=1e-4)


def test_pushing_the_pump_up_overdrives_the_fiber_eom():
    """The 10 mW modulator rating is a real constraint, not a footnote.

    Raising the pump by turning the split half-wave plate dumps *more* light
    into the seed arm, which is why the setup needs a separate polarizer in
    front of the EOM rather than one shared knob.
    """
    at_nominal = build_source_chain()
    assert at_nominal.warnings == ()

    # Same split, but with the seed polarizer wide open: the arm alone carries
    # far more than the modulator tolerates.
    opened = build_source_chain(replace(SetupSettings(), seed_polarizer_deg=90.0))
    assert opened.stage_power("seed", "into fiber EOM") > 0.010
    assert any("exceeds" in w for w in opened.warnings)

    # And turning the split towards more pump makes it worse, not better.
    more_pump = build_source_chain(replace(SetupSettings(), hwp_split_deg=20.0,
                                           seed_polarizer_deg=90.0))
    assert (more_pump.stage_power("seed", "into fiber EOM")
            > opened.stage_power("seed", "into fiber EOM"))


def test_out_of_reach_target_is_reported_rather_than_silently_clipped():
    with pytest.raises(ValueError, match="out of reach"):
        solve_split_angle_deg(50.0)


def test_more_etalon_stages_buy_purity_almost_for_free():
    """The trade-off the seed path does NOT have.

    Each stage costs 20 % of a power budget with a 25x surplus, and buys about
    an order of magnitude of carrier suppression -- so the stage count should be
    chosen on purity alone.
    """
    cal = calib.default_calibration()
    ratios, powers = [], []
    for stages in (2, 3, 4):
        chain = build_source_chain(
            SetupSettings(), cal.with_value("etalon_stages", stages))
        ratios.append(chain.carrier_ratio)
        powers.append(chain.stage_power("seed", f"etalon {stages}"))
    assert ratios[1] < ratios[0] / 5.0
    assert ratios[2] < ratios[1] / 5.0
    assert powers[2] > 0.5 * powers[0]


def test_waist_follows_the_telescope_lenses():
    """Lenses are swappable, so the cell waist is a design variable."""
    base = build_source_chain()
    longer = build_source_chain(replace(SetupSettings(),
                                        seed_telescope_f1_mm=250.0))
    assert longer.seed_waist_m > base.seed_waist_m
    assert longer.seed_waist_m / base.seed_waist_m == pytest.approx(400.0 / 250.0,
                                                                   rel=1e-6)


if __name__ == "__main__":
    import pytest as _pytest
    raise SystemExit(_pytest.main([__file__, "-q"]))
