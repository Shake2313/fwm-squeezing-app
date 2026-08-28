"""
SABES bridge and detection checks.

    python tests/test_sabes_bridge.py   # or: pytest tests/test_sabes_bridge.py

The load-bearing test here is
`test_reference_geometry_reproduces_the_gabes_default_diagnostic`: fed the paper's
own geometry, the whole lab-settings-to-diagnostic pipeline must land on the number
GABES produces when driven by hand. Everything else in SABES is only meaningful if
that holds.
"""
import math
import sys
from dataclasses import replace
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from gabes.schemes import fwm  # noqa: E402
from sabes import bridge, detection  # noqa: E402
from sabes.beamline import SetupSettings, build_source_chain  # noqa: E402
from sabes.calibration import default_calibration  # noqa: E402
from sabes.detection import DetectionSettings  # noqa: E402
from sabes.layout import parts as layout_parts  # noqa: E402


def _reference_overrides():
    """Swap SABES's own geometry for the paper values the GABES defaults use."""
    defaults = fwm.FWMScheme().defaults()
    return dict(pump_waist_um=defaults["pump_waist_um"],
                probe_waist_um=defaults["probe_waist_um"],
                qe_pct=defaults["qe_pct"],
                loss_pct=defaults["loss_pct"],
                seeded_angle_deg=defaults["seeded_angle_deg"])


# ------------------------------------------------------------------ the bridge

def test_bridge_derives_every_knob_it_claims_to_own():
    result = bridge.run(solve=False)
    for key in bridge.DERIVED_KEYS:
        assert key in result.params
        assert np.isfinite(float(result.params[key])), key
    # And it leaves the rest of the scheme alone, so runs stay comparable.
    defaults = fwm.FWMScheme().defaults()
    untouched = set(defaults) - set(bridge.DERIVED_KEYS) - {"mode"}
    for key in untouched:
        assert result.params[key] == defaults[key], key


def test_reference_geometry_reproduces_the_gabes_default_diagnostic():
    """The completion criterion for the bridge: same inputs, same answer."""
    scheme = fwm.FWMScheme()
    defaults = scheme.defaults()
    reference = scheme.compute(defaults)
    reference_db = float(np.asarray(reference["gain_referred_noise_dB"])[
        bridge._operating_index(reference, defaults)])

    result = bridge.run(**_reference_overrides())
    for key in ("opd", "tpd", "temp_c", "cell_mm", "pump_mw", "probe_uw"):
        assert result.params[key] == pytest.approx(defaults[key], rel=1e-5), key
    assert result.gain_referred_noise_db == pytest.approx(reference_db, abs=1e-6)
    # Not bit-identical: the split half-wave plate is solved numerically, so the
    # pump lands on 600.0008 mW rather than exactly 600. That residual is the
    # only thing between the two runs.
    assert result.gains[0] == pytest.approx(
        float(np.asarray(reference["G_s"])[
            bridge._operating_index(reference, defaults)]), rel=1e-5)


def test_native_geometry_differs_only_through_qe_and_waists():
    """The SABES-native answer must be explainable, not merely close."""
    native = bridge.run()
    assert native.params["qe_pct"] == pytest.approx(90.45, abs=0.01)
    assert native.params["pump_waist_um"] == pytest.approx(540.0, abs=1.0)
    assert native.params["probe_waist_um"] == pytest.approx(324.0, abs=1.0)

    # Putting only the QE back recovers most of the gap, because eta sets the
    # floor 10 log10(1 - eta) that every result is measured against.
    qe_only = bridge.run(qe_pct=fwm.FWMScheme().defaults()["qe_pct"])
    reference = bridge.run(**_reference_overrides())
    assert abs(qe_only.gain_referred_noise_db - reference.gain_referred_noise_db) < \
        abs(native.gain_referred_noise_db - reference.gain_referred_noise_db)


def test_two_photon_detuning_needs_no_calibration_coefficient():
    from gabes import constants
    for tpd_mhz in (-8.0, 0.0, +12.5):
        settings = replace(
            SetupSettings(),
            eom_frequency_hz=constants.NU_GROUND_HF - tpd_mhz * 1e6)
        result = bridge.run(settings, solve=False)
        assert result.params["tpd"] == pytest.approx(tpd_mhz, abs=1e-6)


def test_crossing_angle_comes_from_millimetres_at_the_mirrors():
    result = bridge.run(solve=False)
    assert result.params["seeded_angle_deg"] == pytest.approx(0.32, abs=1e-4)


def test_geometry_is_computable_without_the_solve():
    """loss_pct is a solve *input*, so it cannot depend on the solve."""
    result = bridge.run(solve=False)
    assert result.raw is None and result.readout is None
    assert result.params["loss_pct"] > 0.0
    assert result.geometry.eta > 0.0


# --------------------------------------------------------------- detection

def test_loss_pct_is_the_per_optic_product_times_the_measured_residual():
    """The per-optic numbers are guesses; 5.5 % is the one measured total."""
    result = bridge.run(solve=False)
    assert result.params["loss_pct"] == pytest.approx(5.5, abs=0.01)

    calibration = default_calibration()
    product = layout_parts.transmission_of(layout_parts.ROUTE_POST_CELL,
                                           calibration)
    residual = calibration.value("loss_post_cell_residual")
    assert result.geometry.optics_transmission == pytest.approx(
        product * residual, rel=1e-9)


def test_dirtying_one_post_cell_optic_shows_up_in_the_loss():
    """Which is the point of per-optic transmissions over one lumped number."""
    calibration = default_calibration()
    dirty = calibration.with_value("t_lens_probe", 0.80)
    clean = bridge.run(solve=False)
    smudged = bridge.run(calibration=dirty, solve=False)
    assert smudged.params["loss_pct"] > clean.params["loss_pct"] + 15.0
    assert dirty.get("t_lens_probe").provenance == "hand"


def test_separation_margin_saturates_with_distance():
    """The margin rises with distance but is capped at angle / divergence.

    Near the cell the beam barely grows, so separation wins and the margin climbs
    linearly. Far away both grow together and the ratio locks onto
    `theta * pi * w0 / lambda`. Moving the D-shaped mirrors back therefore buys a
    bounded amount, and only a larger waist or a larger angle raises the ceiling.
    """
    settings = SetupSettings()
    chain = build_source_chain(settings)
    margins = []
    for distance in (0.2, 0.5842, 2.0, 8.0):
        # Swept through the layout, not a coefficient: moving the mirrors is a
        # change to the table, and that is now where distances live.
        geom = detection.geometry(
            chain, settings, DetectionSettings(),
            layout=layout_parts.with_detection_distances(
                cell_to_dmirror=distance))
        margins.append(geom.separation_margin)
    assert margins == sorted(margins)                 # monotone increasing
    limit = (math.radians(settings.crossing_angle_deg)
             / chain.seed.mode.divergence_rad)
    assert margins[-1] == pytest.approx(limit, rel=0.02)
    # At the actual 23 inch the margin is already most of the way to the ceiling,
    # so there is little left to win by lengthening the path.
    assert 0.75 * limit < margins[1] < 0.90 * limit


def test_a_much_larger_waist_breaks_the_pump_separation():
    """Waist is a design variable, so the overlap limit has to be checkable."""
    settings = SetupSettings()
    base = detection.geometry(build_source_chain(settings), settings)
    assert base.separation_margin > detection.SEPARATION_MARGIN_WARN
    assert not any("twin radii" in w for w in base.warnings)

    fat = replace(settings, seed_telescope_f1_mm=100.0)   # 4x the seed waist
    wide = detection.geometry(build_source_chain(fat), fat)
    assert wide.separation_margin < detection.SEPARATION_MARGIN_WARN
    assert any("twin radii" in w for w in wide.warnings)


def test_pump_leakage_is_observed_rather_than_derived():
    """A Gaussian tail claimed 1e-9 rejection: useless, and contradicted.

    The paper identifies pump scatter as a real experimental limit, so leakage is
    a measured detector input rather than a geometry-derived noise coefficient.
    """
    settings = SetupSettings()
    chain = build_source_chain(settings)
    geom = detection.geometry(chain, settings)

    quiet = detection.readout(geom, chain, (14.0, 13.0), settings,
                              detection=DetectionSettings(pump_leakage_dbm=-80.0))
    loud = detection.readout(geom, chain, (14.0, 13.0), settings,
                             detection=DetectionSettings(pump_leakage_dbm=-50.0))
    assert loud.arms[0].residual_pump_w == pytest.approx(
        1000.0 * quiet.arms[0].residual_pump_w, rel=1e-6)
    # And it no longer depends on beam geometry at all.
    assert not hasattr(geom.arms[0], "pump_leak_fraction")


def test_clearance_matches_the_closed_form_and_scales_with_power():
    cal = default_calibration()
    responsivity = cal.value("pd_responsivity_a_per_w")
    nep = cal.value("bpd_nep_w_per_rthz")
    settings = SetupSettings()
    chain = build_source_chain(settings)
    geom = detection.geometry(chain, settings)

    read = detection.readout(geom, chain, (14.0, 13.0), settings)
    expected = 20.0 * math.log10(
        math.sqrt(2.0 * detection.ELEMENTARY_CHARGE * responsivity
                  * read.total_power_w) / (responsivity * nep))
    assert read.clearance_db == pytest.approx(expected, abs=1e-9)

    louder = detection.readout(geom, chain, (140.0, 130.0), settings)
    assert louder.clearance_db == pytest.approx(read.clearance_db + 10.0, abs=0.1)


def test_power_for_clearance_inverts_the_clearance_formula():
    target = 20.0
    power = detection.power_for_clearance(target)
    assert power == pytest.approx(2.0e-3, rel=0.1)     # ~2 mW, well under 9 mW
    cal = default_calibration()
    responsivity = cal.value("pd_responsivity_a_per_w")
    nep = cal.value("bpd_nep_w_per_rthz")
    got = 20.0 * math.log10(
        math.sqrt(2.0 * detection.ELEMENTARY_CHARGE * responsivity * power)
        / (responsivity * nep))
    assert got == pytest.approx(target, abs=1e-9)


def test_corrected_absolute_normalization_exposes_paper_power_mismatch():
    """Do not silently retune the inherited 0.74 scale to the paper power.

    Removing the duplicate ground-manifold population factor raises the current
    approximate model's absolute gain far above the measured operating point.
    Until the remaining microscopic response is rebuilt, this disagreement is a
    model diagnostic rather than detector-clearance evidence for the experiment.
    """
    result = bridge.run(resolution=fwm.FIDELITY_ULTRA)
    assert result.readout.total_power_w > 1.0e-3
    assert result.readout.clearance_db > 15.0
    assert result.readout.margin_above_electronic_db(
        result.gain_referred_noise_db) > 8.0
    assert result.readout.measurable
    assert any("linearity limit" in w for w in result.warnings)


def test_more_seed_power_buys_clearance_the_atoms_do_not_notice():
    """The free lunch: seed power is inert for the physics but sets clearance."""
    settings = SetupSettings()
    chain = build_source_chain(settings)
    geom = detection.geometry(chain, settings)
    base = detection.readout(geom, chain, (14.0, 13.0), settings)

    # Trim wide open: no rotation, so the Glan-Taylor passes everything.
    louder = replace(settings, seed_trim_hwp_deg=0.0)
    chain_louder = build_source_chain(louder)
    assert chain_louder.seed_power_w > 10.0 * chain.seed_power_w
    read = detection.readout(detection.geometry(chain_louder, louder),
                             chain_louder, (14.0, 13.0), louder)
    assert read.clearance_db > base.clearance_db + 9.0
    assert read.measurable


def test_hard_focus_onto_the_diode_breaks_the_linearity_limit():
    """f=75/100 mm puts the twins far inside the Thorlabs 4 W/cm^2 limit."""
    result = bridge.run(resolution=fwm.FIDELITY_ULTRA)
    probe = result.readout.arms[0]
    assert probe.spot_radius_m < 50e-6
    assert probe.power_density_w_per_cm2 > 4.0
    assert any("linearity limit" in w for w in result.warnings)


def test_a_longer_focal_length_relieves_the_density_limit():
    # The corrected, still-uncalibrated absolute model predicts much more power
    # than the paper; 600 mm is the first standard-scale choice below the limit.
    gentle = DetectionSettings(probe_lens_focal_mm=600.0,
                               conjugate_lens_focal_mm=600.0)
    result = bridge.run(detection=gentle, resolution=fwm.FIDELITY_ULTRA)
    assert result.readout.arms[0].power_density_w_per_cm2 < 4.0
    assert not any("linearity limit" in w for w in result.warnings)
    # And the spot still fits the S3883 active area comfortably.
    active = 0.5 * default_calibration().value("pd_active_diameter_mm") * 1e-3
    assert result.readout.arms[0].spot_radius_m < active


def test_defocus_also_relieves_the_density_limit():
    result = bridge.run(detection=DetectionSettings(pd_defocus_mm=40.0),
                        resolution=fwm.FIDELITY_ULTRA)
    assert result.readout.arms[0].power_density_w_per_cm2 < 4.0


def test_focal_length_helper_hits_the_requested_density():
    wavelength = 795e-9
    radius_at_lens = 700e-6
    focal = detection.focal_length_for_density(4.0, 100e-6, radius_at_lens,
                                              wavelength)
    spot = wavelength * focal / (math.pi * radius_at_lens)
    assert 2.0 * 100e-6 / (math.pi * spot ** 2) / 1e4 == pytest.approx(4.0,
                                                                      rel=1e-9)


# --------------------------------------------------------------- reporting

def test_derived_table_pairs_every_setting_with_its_gabes_quantity():
    result = bridge.run(solve=False)
    rows = bridge.derived_table(result)
    assert len(rows) >= 8
    for setting_label, setting_value, derived_label, derived_value in rows:
        assert setting_label and setting_value and derived_label and derived_value
    joined = " ".join(r[2] for r in rows)
    for expected in ("TPD", "Pump power", "Seed power", "Crossing angle",
                     "Detector QE"):
        assert expected in joined


def test_warnings_aggregate_from_every_stage():
    """A fat seed beam makes the pump-separation warning fire in geometry, and
    the readout adds its own on top."""
    fat = replace(SetupSettings(), seed_telescope_f1_mm=100.0)
    result = bridge.run(fat, resolution=fwm.FIDELITY_ULTRA)
    assert any("twin radii" in w for w in result.geometry.warnings)
    assert len(result.warnings) >= len(result.geometry.warnings)


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-q"]))
