"""Measured-reference units, support and explicit apparatus-context controls."""

import hashlib
import json
from pathlib import Path

import numpy as np
import pytest
from scipy.interpolate import PchipInterpolator

from gabes import fwm_measured_gain as measured
from gabes_ui.export import to_json_primitive


ROOT = Path(__file__).resolve().parents[1]
EVIDENCE = ROOT/"docs/grand_challenge/gain_calibration_2026_10_16"


def test_observations_are_raw_workbook_voltages_with_common_pump_off_denominator():
    document = json.loads((EVIDENCE/"workbook_ooxml.json").read_text(encoding="utf-8"))
    cells = {cell["cell"]: cell["value"] for row in document["sheets"][0]["rows"] for cell in row}
    actual = measured.observations()
    np.testing.assert_array_equal(actual["eom_MHz"], [cells[f"A{row}"] for row in range(5, 28)])
    for arm, column in (("probe", "C"), ("conjugate", "D")):
        voltages = np.array([cells[f"{column}{row}"] for row in range(5, 28)])
        np.testing.assert_array_equal(actual[f"{arm}_voltage_mV"], voltages)
        np.testing.assert_array_equal(actual[f"{arm}_ratio"], voltages/cells["B1"])
    assert actual["probe_ratio"][-1] == .86/35.
    assert actual["probe_ratio"][-1] != .8046*.86/35.
    actual["probe_ratio"][:] = -1
    assert np.all(measured.observations()["probe_ratio"] > 0)


def test_frequency_sign_and_unrounded_hyperfine_origin():
    eom = np.array([3022., 3033., 3034., 3044.])
    delta = measured.internal_delta_mhz_from_eom(eom)
    np.testing.assert_allclose(delta, [13.732439, 2.732439, 1.732439, -8.267561], atol=1e-9)
    np.testing.assert_allclose(measured.eom_mhz_from_internal_delta(delta), eom, atol=1e-12)
    # Workbook's spacing coordinate is EOM-3036, not -delta_internal.
    np.testing.assert_allclose(eom-3036., -delta-.267561, atol=1e-9)


def test_strict_support_no_extrapolation_or_passive_gain_floor():
    frequency = np.array([3021.999999, 3022., 3033., 3044., 3044.000001])
    actual = measured.reference_curve(frequency)
    np.testing.assert_array_equal(actual["within_support"], [False, True, True, True, False])
    for arm in ("probe", "conjugate"):
        values = actual[f"measured_gain_{arm}"]
        assert np.all(np.isnan(values[[0, 4]]))
        assert np.all(np.isfinite(values[1:4]))
        assert 0 < values[1] < 1 < values[2]
    assert actual["status"][0] == "OUTSIDE_MEASURED_RF_SUPPORT"


def test_both_arm_arrays_keep_scalar_and_multidimensional_shapes():
    scalar = measured.reference_curve(3033.)
    assert np.shape(scalar["measured_gain_probe"]) == ()
    assert np.shape(scalar["measured_gain_conjugate"]) == ()
    frequency = np.array([[3022., 3023.], [3043., 3044.]])
    actual = measured.reference_curve(frequency)
    assert actual["measured_gain_probe"].shape == frequency.shape
    assert actual["measured_gain_conjugate"].shape == frequency.shape
    assert np.all(actual["measured_gain_probe"] > 0)
    assert np.all(actual["measured_gain_conjugate"] > 0)


def test_interpolator_reproduces_raw_observations_without_full_cell_scaling():
    data = measured.observations()
    actual = measured.reference_curve(data["eom_MHz"])
    for arm in ("probe", "conjugate"):
        np.testing.assert_allclose(actual[f"measured_gain_{arm}"], data[f"{arm}_ratio"],
                                   rtol=5e-15, atol=0.)
    assert np.all(actual["status"] == "SAMPLED_EMPIRICAL_REFERENCE_INTERPOLATION")
    assert not actual["profile"]["predictive_calibration_accepted"]
    assert not actual["profile"]["independent_validation"]


def test_shape_preserving_intervals_do_not_overshoot_either_measured_arm():
    data = measured.observations()
    frequencies = data["eom_MHz"][::-1]
    for index in range(len(frequencies)-1):
        sample = np.linspace(frequencies[index], frequencies[index+1], 33)
        actual = measured.reference_curve(sample)
        for arm in ("probe", "conjugate"):
            left, right = data[f"{arm}_ratio"][::-1][index:index+2]
            values = actual[f"measured_gain_{arm}"]
            tolerance = 16*np.finfo(float).eps*max(left, right)
            assert np.min(values) >= min(left, right)-tolerance
            assert np.max(values) <= max(left, right)+tolerance


def test_probe_measured_peak_plateau_and_conjugate_peak_are_preserved():
    actual = measured.reference_curve(np.linspace(3033., 3034., 101))
    np.testing.assert_allclose(actual["measured_gain_probe"], 870./35., rtol=5e-15)
    assert np.max(actual["measured_gain_conjugate"]) == pytest.approx(890./35., rel=5e-15)
    assert measured.reference_curve(3033.)["measured_gain_conjugate"] > (
        measured.reference_curve(3034.)["measured_gain_conjugate"])


def test_reference_record_identity_and_observation_contract_are_preserved():
    profile = measured.profile()
    body = (ROOT/profile["reference_record"]).read_bytes()
    assert hashlib.sha256(body).hexdigest() == measured.REFERENCE_RECORD_SHA256
    record = json.loads(body)
    assert record["source_workbook_sha256"] == measured.SOURCE_WORKBOOK_SHA256
    assert record["model_id"] == measured.MODEL_ID
    data = measured.observations()
    for key in ("eom_MHz", "probe_ratio", "conjugate_ratio"):
        np.testing.assert_array_equal(record["source_points"][key], data[key][::-1])
    assert record["quality"]["sampled_interpolation_validation_passed"]
    assert not record["quality"]["development_acceptance_passed"]
    assert profile["sampled_interpolation_validation_passed"]
    assert not profile["predictive_calibration_accepted"]


def test_fixed_alternating_diagnostic_is_distinct_from_exact_node_agreement():
    # Fit the interpolator only to the original even rows, then evaluate the
    # excluded odd rows. This checks the reported diagnostic independently of
    # the production reference, which contains all 23 observations.
    data = measured.observations()
    train, validation = np.arange(0, 23, 2), np.arange(1, 23, 2)
    profile = measured.profile()
    for arm in ("probe", "conjugate"):
        frequency, ratios = data["eom_MHz"], data[f"{arm}_ratio"]
        interpolation = PchipInterpolator(frequency[train][::-1], np.log(ratios[train][::-1]),
                                           extrapolate=False)
        prediction = np.exp(interpolation(frequency[validation]))
        maximum = float(np.max(abs(prediction/ratios[validation]-1)))
        assert maximum == pytest.approx(
            profile["alternating_interpolation_max_relative_error"][arm], rel=1e-12)
        assert maximum > .1
    assert profile["alternating_interpolation_diagnostic_passed"]
    assert not profile["development_acceptance_passed"]
    assert profile["failed_predictive_trials"]["chebyshev_best_joint_max_relative_error"] > .2


def test_invalid_frequencies_remain_unavailable_in_strict_json():
    result = measured.reference_curve([np.nan, np.inf, -np.inf, 3033.])
    payload = json.loads(json.dumps(to_json_primitive(result), allow_nan=False))
    assert payload["measured_gain_probe"][:3] == [None, None, None]
    assert payload["measured_gain_conjugate"][:3] == [None, None, None]
    assert payload["status"][:3] == ["INVALID_RF_FREQUENCY"]*3
    assert payload["profile"]["pump_off_atomic_transmission"] is None


def test_explicit_conditional_profile_is_disabled_by_default():
    context = measured.nominal_context()
    assert not measured.context_eligibility(context)["applied"]
    enabled = measured.context_eligibility(context, enabled=True)
    assert enabled["applied"] and enabled["eligible"]
    assert enabled["status"] == "CONDITIONAL_MEASURED_REFERENCE_PROFILE"
    profile = enabled["profile"]
    assert not profile["independently_calibrated_vapor_input"]
    assert not profile["full_cell_transmission_is_atomic_A0"]
    assert not profile["atomic_gain_modified"]
    assert not profile["canonical_transfer_modified"]
    assert not profile["NF_modified"]
    assert not profile["SABES_detector_power_gain_modified"]


def test_measured_seed_power_is_not_accepted_as_the_assumed_atomic_input_proxy():
    context = measured.nominal_context()
    expected = measured.REFERENCE_SEED_POWER_W/(.8046*.9709)
    assert context["P_probe"] == pytest.approx(expected)
    context["P_probe"] = measured.REFERENCE_SEED_POWER_W
    assert not measured.context_eligibility(context, enabled=True)["applied"]


def test_seed_proxy_assumption_must_be_explicit_even_when_numbers_match():
    context = measured.nominal_context()
    context.pop("seed_input_assumption")
    ledger = measured.context_eligibility(context, enabled=True)
    assert ledger["context_matches"]
    assert not ledger["seed_proxy_assumption_adopted"]
    assert not ledger["eligible"] and not ledger["applied"]


@pytest.mark.parametrize("key", ("D_GHz", "T", "P_pump", "P_probe", "w_pump", "w_probe",
                                  "L", "pump_probe_angle_deg", "branch"))
def test_context_gate_rejects_changed_physical_inputs(key):
    context = measured.nominal_context()
    context[key] *= 1.001
    ledger = measured.context_eligibility(context, enabled=True)
    assert not ledger["eligible"] and not ledger["applied"]
    assert not ledger["numeric_context_matches"][key]


def test_roundoff_allowance_is_not_a_physical_opd_acceptance_interval():
    context = measured.nominal_context()
    context["P_pump"] = .2948040000000001
    assert measured.context_eligibility(context, enabled=True)["eligible"]
    context["D_GHz"] += 1e-6
    assert not measured.context_eligibility(context, enabled=True)["eligible"]


@pytest.mark.parametrize("flag", ("legacy_mean_field_fixture", "legacy_gain_calibration_fixture"))
def test_current_profile_cannot_compose_with_historical_fixtures(flag):
    context = dict(measured.nominal_context(), **{flag: True})
    assert not measured.context_eligibility(context, enabled=True)["applied"]


def test_detection_and_retired_switches_do_not_define_empirical_reference_context():
    context = dict(measured.nominal_context(), after_cell_eff_pct=0., qe_pct=0.,
                   seed_classical_rin_scaled=100., gain_closure_enabled=True)
    ledger = measured.context_eligibility(context, enabled=True)
    assert ledger["applied"]
    assert not measured.context_eligibility(context)["applied"]


def test_missing_and_nonfinite_context_fail_closed_without_metadata_repair():
    missing = measured.nominal_context()
    missing.pop("T")
    assert not measured.context_eligibility(missing, enabled=True)["eligible"]
    invalid = dict(measured.nominal_context(), P_probe=np.nan)
    assert not measured.context_eligibility(invalid, enabled=True)["eligible"]
