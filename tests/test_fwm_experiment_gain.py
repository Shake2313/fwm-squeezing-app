"""Experimental model calibration continues outside support and tracks controls."""
from copy import deepcopy
import json
from pathlib import Path

import numpy as np
import pytest

from gabes import fwm_experiment_gain as practical
from gabes import fwm_measured_gain as measured


RECORD = Path(__file__).resolve().parents[1]/"analysis/fwm_experiment_20261022/gain_profile.json"


@pytest.fixture(scope="module")
def nominal_raw():
    record = json.loads(RECORD.read_text(encoding="utf-8"))
    context = measured.nominal_context()
    return {
        "probe_axis_GHz": np.array(record["reference_model_internal_delta_mhz"])*1e-3,
        "raman_center_minus_GHz": 0.0,
        "G_s": np.array(record["reference_model_G_s"]),
        "G_c": np.array(record["reference_model_G_c"]),
        "D_GHz": context["D_GHz"], "temperature_K": context["T"],
        "input_pump_power_W": context["P_pump"], "input_seed_power_W": context["P_probe"],
        "w_pump_m": context["w_pump"], "w_probe_m": context["w_probe"],
        "cell_length_m": context["L"], "pump_probe_angle_deg": context["pump_probe_angle_deg"],
        "branch": -1,
    }


def test_nominal_model_calibration_recovers_raw_measurements_and_moves_peaks(nominal_raw):
    observed = measured.observations()
    nodes = measured.internal_delta_mhz_from_eom(observed["eom_MHz"])
    result = practical.evaluate(nominal_raw, nodes)
    for arm, key in (("probe", "G_s"), ("conjugate", "G_c")):
        np.testing.assert_allclose(result[key], observed[arm+"_ratio"], rtol=2e-13, atol=1e-13)
    dense = np.linspace(-12., 18., 3001)
    corrected = practical.evaluate(nominal_raw, dense)
    assert np.argmax(nominal_raw["G_s"]) == np.argmax(nominal_raw["G_c"])
    for key in ("G_s", "G_c"):
        peak = dense[np.argmax(corrected[key])]
        assert 1.7 <= peak <= 3.1
    assert not result["context_changed"]


def test_affine_fit_and_frozen_coefficients_have_reproducible_source_metadata():
    record = json.loads(RECORD.read_text(encoding="utf-8"))
    profile = practical.profile()
    assert profile["measurement_date"] == "2026-10-13"
    assert profile["source_workbook_sha256"] == measured.SOURCE_WORKBOOK_SHA256
    assert not profile["independent_validation"]
    assert not profile["noise_covariance_calibrated"]
    for arm in ("probe", "conjugate"):
        fitted = profile["arms"][arm]
        assert fitted["gain_scale"] > 0
        assert fitted["frequency_scale"] > 0
        for key, value in fitted.items():
            np.testing.assert_allclose(value, record["arms"][arm][key], rtol=1e-13, atol=1e-13)


def test_wide_frequency_continuation_uses_current_model_not_a_clamped_gain(nominal_raw):
    delta = np.array([-1e6, -550., -100., -25., 25., 100., 550., 1e6])
    first = practical.evaluate(nominal_raw, delta)
    assert first["applied"]
    assert not np.any(first["within_measured_support"])
    assert np.all(first["point_status"] == "EXTRAPOLATED_FREQUENCY")
    for key in ("G_s", "G_c"):
        assert np.all(np.isfinite(first[key]))
        assert np.all(first[key] >= 0)
        assert np.ptp(first[key]) > .01
    changed = deepcopy(nominal_raw)
    changed["G_s"] *= 1.4
    changed["G_c"] *= .6
    changed["input_pump_power_W"] *= 1.1
    changed["temperature_K"] += 5
    second = practical.evaluate(changed, [-100., -25., 25., 100.])
    original = practical.evaluate(nominal_raw, [-100., -25., 25., 100.])
    np.testing.assert_allclose(second["G_s"], 1.4*original["G_s"], rtol=1e-13)
    np.testing.assert_allclose(second["G_c"], .6*original["G_c"], rtol=1e-13)
    assert second["applied"] and second["context_changed"]
    assert "input_pump_power_W" in second["changed_context_inputs"]
    assert "temperature_K" in second["changed_context_inputs"]


@pytest.mark.parametrize("context_key,value", [
    ("D_GHz", -1.4), ("temperature_K", 395.), ("input_pump_power_W", .6),
    ("input_seed_power_W", 8e-6), ("w_pump_m", 700e-6), ("w_probe_m", 450e-6),
    ("cell_length_m", .02), ("pump_probe_angle_deg", .5), ("branch", 1),
])
def test_operating_conditions_never_disable_requested_correction(nominal_raw, context_key, value):
    changed = dict(nominal_raw, **{context_key: value})
    result = practical.calibrate_spectrum(changed, enabled=True, selected_delta_mhz=150.)
    assert result["applied"] and result["context_changed"]
    point = result["selected_point"]
    assert point["available"] and point["applied"]
    assert point["G_s"] >= 0 and point["G_c"] >= 0
    assert "EXTRAPOLATED_FREQUENCY" in point["status"]
    assert "CHANGED_OPERATING_CONDITIONS" in point["status"]


def test_disabled_identity_nonmutation_and_profile_copy(nominal_raw):
    before = deepcopy(nominal_raw)
    result = practical.calibrate_spectrum(nominal_raw, enabled=False)
    assert not result["applied"]
    for key in ("G_s", "G_c"):
        np.testing.assert_array_equal(result[key], nominal_raw[key])
        result[key][:] = -1
        np.testing.assert_array_equal(nominal_raw[key], before[key])
    practical.calibrate_spectrum(nominal_raw, enabled=True)
    for key in nominal_raw:
        np.testing.assert_array_equal(nominal_raw[key], before[key])
    profile = practical.profile()
    profile["arms"]["probe"]["gain_scale"] = -1
    assert practical.profile()["arms"]["probe"]["gain_scale"] > 0


def test_zero_conjugate_remains_zero_instead_of_affine_background(nominal_raw):
    raw = dict(nominal_raw, G_c=np.zeros_like(nominal_raw["G_c"]))
    result = practical.evaluate(raw, [-1e6, -10., 2., 15., 1e6])
    np.testing.assert_array_equal(result["G_c"], np.zeros(5))


def test_spectrum_order_does_not_change_calibrated_evaluation(nominal_raw):
    reversed_raw = dict(nominal_raw)
    for key in ("probe_axis_GHz", "G_s", "G_c"):
        reversed_raw[key] = nominal_raw[key][::-1]
    first = practical.evaluate(nominal_raw, [-15., 2.2, 20.])
    second = practical.evaluate(reversed_raw, [-15., 2.2, 20.])
    np.testing.assert_array_equal(first["G_s"], second["G_s"])
    np.testing.assert_array_equal(first["G_c"], second["G_c"])


def test_invalid_inputs_are_reported_instead_of_generating_finite_fake_gains(nominal_raw):
    with pytest.raises(ValueError, match="finite"):
        practical.evaluate(nominal_raw, [np.nan])
    with pytest.raises(ValueError, match="nonnegative"):
        practical.calibrate_spectrum(dict(nominal_raw, G_s=-nominal_raw["G_s"]), enabled=True)


def test_live_atomic_solver_matches_calibration_and_changed_physics_stays_active():
    from gabes.schemes.fwm import FIDELITY_BALANCED, FWMScheme

    scheme = FWMScheme()
    context = measured.nominal_context()
    params = dict(scheme.defaults(), opd=context["D_GHz"], temp_c=context["T"]-273.15,
                  pump_mw=context["P_pump"]*1e3, probe_uw=context["P_probe"]*1e6,
                  cell_mm=context["L"]*1e3, pump_waist_um=context["w_pump"]*1e6,
                  probe_waist_um=context["w_probe"]*1e6,
                  seeded_angle_deg=context["pump_probe_angle_deg"], resolution=FIDELITY_BALANCED)
    observed = measured.observations()
    nodes = measured.internal_delta_mhz_from_eom(observed["eom_MHz"])
    raw = scheme.compute(params)
    nominal = practical.evaluate(raw, nodes)
    # This integration check deliberately uses a fresh current atomic solve,
    # rather than only the archived nominal arrays used by the unit tests.
    np.testing.assert_allclose(nominal["G_s"], observed["probe_ratio"], rtol=.01)
    np.testing.assert_allclose(nominal["G_c"], observed["conjugate_ratio"], rtol=.01)
    changed_raw = scheme.compute(dict(params, pump_mw=350., temp_c=118.))
    changed = practical.calibrate_spectrum(changed_raw, enabled=True, selected_delta_mhz=150.)
    assert changed["applied"] and changed["context_changed"]
    assert changed["selected_point"]["available"]
    assert not np.allclose(changed["G_s"], practical.calibrate_spectrum(raw, enabled=True)["G_s"])
    assert np.all(np.isfinite(changed["G_s"]))
    assert np.all(np.isfinite(changed["G_c"]))
