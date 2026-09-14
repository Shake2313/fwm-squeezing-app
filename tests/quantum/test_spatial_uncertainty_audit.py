"""Independent geometry and evidence checks for spatial input uncertainty."""

from dataclasses import asdict, replace
import hashlib
import json
from types import SimpleNamespace

import numpy as np
import pytest

from analysis.grand_challenge import spatial_uncertainty_audit as audit
from analysis.grand_challenge.normalization_audit import conditional_inputs
from gabes import constants as c, core
from gabes.fwm_quantum.normalization import optical_carriers
from gabes.fwm_quantum.kinetic import CarrierGeometry
from gabes.quantum.contracts import AnalysisFrequencyAxis, ParameterEvidence
from gabes.quantum.readout import DetectorResponse, IntensityDifferenceSpectrum
from gabes.quantum.spectrum_analysis import SpectralTrace
from gabes.quantum.uncertainty import JointInputUncertainty, propagate_first_order


@pytest.fixture(autouse=True)
def single_thread():
    with core.blas_single_thread():
        yield


def detector():
    return DetectorResponse(AnalysisFrequencyAxis.from_hz([.1e6, 1e6, 4e6]),
        [.85, .85], np.ones((3, 2)), 1., np.zeros(3),
        "conditional detector fixture; not a measured calibration")


def evaluator(*, inputs=None, solver=None):
    return audit.SpatialUncertaintyEvaluator(
        conditional_inputs() if inputs is None else inputs, detector(), solver=solver)


def test_angle_samples_rebuild_physical_carriers_with_the_analytic_vector_derivative():
    original = conditional_inputs()
    provider = evaluator(inputs=original)
    theta, step = -.004, 1e-5
    center = provider.medium([original.pump_power_W, theta], order_x=4)
    plus = provider.medium([original.pump_power_W, theta+step], order_x=4)
    minus = provider.medium([original.pump_power_W, theta-step], order_x=4)
    k_conjugate = optical_carriers(original.detunings)[2]/c.C_LIGHT
    expected = k_conjugate*np.array([np.cos(theta), 0., -np.sin(theta)])
    actual = (plus.geometry.wavevectors_rad_m[2]-minus.geometry.wavevectors_rad_m[2])/(2*step)
    np.testing.assert_allclose(actual, expected, rtol=2e-8, atol=2e-4)
    for changed in (plus, minus):
        np.testing.assert_array_equal(changed.geometry.wavevectors_rad_m[:2], center.geometry.wavevectors_rad_m[:2])
        assert changed.inputs == original
        np.testing.assert_allclose(np.linalg.norm(changed.geometry.wavevectors_rad_m, axis=1),
                                   optical_carriers(original.detunings)/c.C_LIGHT, rtol=2e-15)
    # The spatial loop is kc+kp-2*k0. Its angle derivative must survive the
    # adapter; a scalar mismatch replacement or a frozen geometry would fail.
    actual_loop = (-plus.geometry.mismatch_rad_m+minus.geometry.mismatch_rad_m)/(2*step)
    np.testing.assert_allclose(actual_loop, expected, rtol=2e-8, atol=2e-4)
    assert not center.geometry.closure_audit()["passed"]
    assert original.phase_mismatch_rad_m == 0.


def test_power_samples_rebuild_atomic_drive_and_leave_fixed_inputs_and_modes_bound():
    original = conditional_inputs()
    provider = evaluator(inputs=original)
    center = provider.medium([.6, -.004], order_x=4)
    changed = provider.medium([.606, -.004], order_x=8)
    assert changed.inputs == replace(original, pump_power_W=.606)
    assert original.pump_power_W == .6
    assert not np.array_equal(center.h0, changed.h0)
    np.testing.assert_array_equal(center.geometry.wavevectors_rad_m, changed.geometry.wavevectors_rad_m)
    assert len(center.modes.positions_m) == 4
    assert len(changed.modes.positions_m) == 8
    for medium in (center, changed):
        np.testing.assert_allclose(medium.modes.area_weights_m2@abs(medium.modes.mode_values_m_inverse)**2,
                                   [1., 1.], atol=1e-15)
        assert medium.mean_order == 4 and medium.response_order == 3


def test_scalar_mismatch_cannot_be_added_to_rebuilt_vector_geometry():
    with pytest.raises(ValueError, match="scalar mismatch"):
        evaluator(inputs=replace(conditional_inputs(), phase_mismatch_rad_m=1.)).medium([.6, -.004], order_x=2)


def test_real_no_atom_spatial_solver_preserves_shot_noise_for_both_uncertain_inputs():
    provider = evaluator(inputs=replace(conditional_inputs(), number_density_m3=0.))
    first = provider.evaluate([.6, -.004], order_x=2)
    second = provider.evaluate([.61, -.0039], order_x=3)
    for row in (first, second):
        np.testing.assert_allclose(row["outputs"], [1., 0., 1., 1., 1.], atol=1e-13, rtol=0)
        vectors = np.asarray(row["wavevectors_rad_m"])
        np.testing.assert_allclose(row["loop_wavevector_rad_m"], vectors[1]+vectors[2]-2*vectors[0],
                                   atol=4e-9, rtol=0)
    assert not np.array_equal(first["wavevectors_rad_m"], second["wavevectors_rad_m"])
    np.testing.assert_array_equal(provider.callback(2)([.6, -.004]), first["outputs"])


@pytest.mark.parametrize("correlation", [-.6, .6])
def test_shared_gain_and_ratio_callback_propagates_the_full_correlated_input_covariance(monkeypatch, correlation):
    provider = evaluator()
    # The fixture reads the angle from the actual rebuilt optical vector. Its
    # response is analytic; no atomic or experimental noise is inferred here.
    derivative = np.array([[2., 1000.], [.5, -2000.], [.3, 20.], [.4, 30.], [.2, -10.]])
    nominal = np.array([2., 1., .8, .7, .9])
    calls = []
    def analytic_evaluate(values, *, order_x, propagation_rtol=2e-9):
        medium = provider.medium(values, order_x=order_x)
        kc = medium.geometry.wavevectors_rad_m[2]
        angle = np.arctan2(kc[0], kc[2])
        offsets = np.array([medium.inputs.pump_power_W-.6, angle+.004])
        calls.append((order_x, propagation_rtol))
        return {"outputs": nominal + derivative@offsets}
    monkeypatch.setattr(provider, "evaluate", analytic_evaluate)
    parameters = (
        ParameterEvidence("pump_power_W", .6, "W", "assumed", uncertainty=.006),
        ParameterEvidence("conjugate_angle_rad", -.004, "rad", "assumed", uncertainty=1e-5),
    )
    model = JointInputUncertainty(parameters, [[1., correlation], [correlation, 1.]])
    result = propagate_first_order(model, provider.callback(8, propagation_rtol=2e-11),
        steps=[.0006, 1e-6], output_ids=("Gp", "Gc", "R0", "R1", "R2"), output_units=("1",)*5)
    expected_covariance = derivative@model.covariance@derivative.T
    np.testing.assert_allclose(result.nominal, nominal, atol=1e-13, rtol=0)
    np.testing.assert_allclose(result.jacobian, derivative, rtol=2e-10, atol=2e-10)
    np.testing.assert_allclose(result.covariance, expected_covariance, rtol=5e-10, atol=1e-18)
    assert calls == [(8, 2e-11)]*5
    assert np.any(result.covariance[:2, 2:] != 0)
    assert not result.model.audit_evidence().passed
    assert not result.experimental_validation


def toy_solver(*, failure=None):
    """Synthetic adapter records, explicitly not independent atomic validation."""
    calls = []
    def solve(medium, detector, *, propagation_rtol):
        calls.append((medium, propagation_rtol))
        ratio = np.array([.8, .7, .9])
        sql = np.full(3, 2e-24)
        spectrum = IntensityDifferenceSpectrum(
            detector.analysis_axis, ratio*sql, sql, np.zeros(3), ratio, ratio,
            [1e12, 1e12], 0., 0., "synthetic adapter fixture")
        result = {"transfer": SimpleNamespace(audit=lambda: {"passed": True}),
            "spectrum": spectrum, "probe_power_gain": 2., "conjugate_power_gain": 1.,
            "minimum_channel_cp_eigenvalue": 0., "quantum_propagation": {"fixture": True},
            "mean_ode_evaluations": 1}
        if failure == "transfer":
            result["transfer"] = SimpleNamespace(audit=lambda: {"passed": False})
        elif failure == "channel":
            result["minimum_channel_cp_eigenvalue"] = -1e-3
        elif failure == "nan_channel":
            result["minimum_channel_cp_eigenvalue"] = np.nan
        elif failure == "covariance":
            result["spectrum"] = replace(spectrum, minimum_covariance_uncertainty_eigenvalue=-1e-3)
        elif failure == "gain":
            result["probe_power_gain"] = -1.
        elif failure == "sql":
            result["spectrum"] = replace(spectrum, sql_psd_A2_Hz=[0., 2e-24, 2e-24])
        elif failure == "negative_psd":
            result["spectrum"] = replace(spectrum, quantum_psd_A2_Hz=-ratio*sql)
        elif failure == "axis":
            result["spectrum"] = replace(spectrum, analysis_axis=AnalysisFrequencyAxis.from_hz([2e5, 1e6, 4e6]))
        return result
    return solve, calls


def test_cached_records_are_owned_and_separate_grids_and_tolerances_recompute():
    solver, calls = toy_solver()
    provider = evaluator(solver=solver)
    first = provider.evaluate([.6, -.004], order_x=4)
    first["outputs"][0] = 999.
    first["band_analysis"]["bands"][0]["width_hz"] = 123.
    replay = provider.evaluate([.6, -.004], order_x=4)
    assert len(calls) == 1 and replay["outputs"][0] == 2.
    assert replay["band_analysis"]["bands"][0]["width_hz"] is None
    provider.evaluate([.6, -.004], order_x=8)
    provider.evaluate([.6, -.004], order_x=8, propagation_rtol=2e-11)
    provider.evaluate([.60000000001, -.004], order_x=8, propagation_rtol=2e-11)
    assert len(calls) == 4
    assert len(provider.evaluations) == 4
    assert provider.evaluations[0]["outputs"][0] == 2.
    band = replay["band_analysis"]
    assert band["complete_domain"] and len(band["bands"]) == 1
    assert band["bands"][0]["lower"]["kind"] == "domain"
    assert band["bands"][0]["upper"]["kind"] == "domain"
    assert band["bands"][0]["observed_span_hz"] == pytest.approx(3.9e6)


@pytest.mark.parametrize("failure", ["transfer", "channel", "nan_channel", "covariance", "gain", "sql", "negative_psd", "axis"])
def test_failed_solver_or_readout_is_not_cached_as_a_valid_uncertainty_sample(failure):
    solver, calls = toy_solver(failure=failure)
    provider = evaluator(solver=solver)
    for _ in range(2):
        with pytest.raises(ValueError):
            provider.evaluate([.6, -.004], order_x=4)
    assert len(calls) == 2 and provider.evaluations == []


@pytest.mark.parametrize("mask_name", ["valid", "converged"])
def test_unavailable_or_unconverged_rf_masks_cannot_be_promoted_to_a_complete_callback(monkeypatch, mask_name):
    solver, _ = toy_solver()
    original_adapter = SpectralTrace.from_intensity_difference
    def masked(spectrum, **kwargs):
        return replace(original_adapter(spectrum, **kwargs), **{mask_name: [True, False, True]})
    monkeypatch.setattr(SpectralTrace, "from_intensity_difference", masked)
    provider = evaluator(solver=solver)
    with pytest.raises(ValueError, match="physical/current-spectrum"):
        provider.evaluate([.6, -.004], order_x=4)
    assert provider.evaluations == []


def test_declared_uncertainty_model_retains_units_correlations_and_assumed_provenance():
    model = audit.uncertainty_model(correlation=.6)
    assert model.parameter_ids == ("pump_power_W", "conjugate_angle_rad")
    assert model.units == ("W", "rad")
    np.testing.assert_array_equal(model.values, [.6, -.004])
    np.testing.assert_array_equal(model.standard_uncertainties, [.006, 1e-5])
    np.testing.assert_array_equal(model.correlation, [[1., .6], [.6, 1.]])
    assert all(item.status == "assumed" for item in model.parameters)
    assert model.correlation_evidence.status == "assumed"
    assert not model.audit_evidence().passed


def parent_fixture():
    inputs = conditional_inputs()
    vectors = CarrierGeometry.vacuum_beams(inputs, probe_angle_rad=.005, conjugate_angle_rad=-.004).wavevectors_rad_m
    source = "gabes/fwm_quantum/spatial_cell.py"
    return {"schema": "gabes.stationary_spatial_field_audit.v2",
        "expected_controls_passed": True, "source_stable_during_run": True,
        "velocity_m_s": [0., 0., 0.], "detector_frequencies_Hz": [1e5, 1e6, 4e6],
        "dipole_convention": "uniform-zeeman-rms", "inputs": asdict(inputs),
        "cells": [{"name": f"unequal_Nx{order}", "passed": True,
                   "wavevectors_rad_m": vectors.tolist()} for order in (8, 12)],
        "source_sha256": {source: hashlib.sha256((audit.ROOT/source).read_bytes()).hexdigest()},
        "consumed_scalar_ledger": inputs.consumed_inputs(detector(), convention="uniform-zeeman-rms")["scalars"],
        "detector_transmissions": [.85, .85], "detector_balance": 1.}


def test_historical_parent_hash_mismatches_remain_explicit_without_rewriting_evidence(tmp_path):
    parent = parent_fixture()
    key = next(iter(parent["source_sha256"]))
    parent["source_sha256"][key] = "0"*64
    path = tmp_path/"parent.json"
    path.write_text(json.dumps(parent), encoding="utf-8")
    before = path.read_bytes()
    loaded, provenance = audit.load_parent(path)
    assert loaded == json.loads(before) and path.read_bytes() == before
    assert provenance["sha256"] == hashlib.sha256(before).hexdigest()
    assert not provenance["parent_code_matches_current_sources"]
    assert len(provenance["changed_parent_sources"]) == 1
    mismatch = provenance["changed_parent_sources"][0]
    assert mismatch["path"] == key and mismatch["parent_sha256"] == "0"*64
    assert mismatch["current_sha256"] == hashlib.sha256((audit.ROOT/key).read_bytes()).hexdigest()
    assert "current nominal recomputed" in provenance["scope"]


@pytest.mark.parametrize("change", ["schema", "failed", "unstable", "moving", "frequency", "convention", "scalar", "geometry", "source_escape"])
def test_parent_must_match_stationary_geometry_and_provenance_scope(tmp_path, change):
    parent = parent_fixture()
    if change == "schema":
        parent["schema"] = "gabes.stationary_spatial_field_audit.v1"
    elif change == "failed":
        parent["expected_controls_passed"] = False
    elif change == "unstable":
        parent["source_stable_during_run"] = False
    elif change == "moving":
        parent["velocity_m_s"] = [0., 0., 1.]
    elif change == "frequency":
        parent["detector_frequencies_Hz"][0] = 2e5
    elif change == "convention":
        parent["dipole_convention"] = "legacy-reciprocal"
    elif change == "scalar":
        parent["inputs"]["phase_mismatch_rad_m"] = 1.
    elif change == "geometry":
        parent["cells"][0]["wavevectors_rad_m"][2][0] += 10.
    elif change == "source_escape":
        parent["source_sha256"] = {"../outside.py": "0"*64}
    path = tmp_path/"invalid.json"
    path.write_text(json.dumps(parent), encoding="utf-8")
    with pytest.raises(ValueError):
        audit.load_parent(path)


def test_parent_detector_keeps_complex_response_electronics_and_balance_from_ledger():
    parent = parent_fixture()
    parent["detector_balance"] = 1.2
    parent["consumed_scalar_ledger"]["probe_response_imag_1"] = [.25, "A/A"]
    parent["consumed_scalar_ledger"]["conjugate_response_real_2"] = [.9, "A/A"]
    parent["consumed_scalar_ledger"]["electronics_psd_0"] = [3e-26, "A^2/Hz"]
    actual = audit.default_detector(parent)
    assert actual.balance == 1.2
    assert actual.current_response[1, 0] == 1+.25j
    assert actual.current_response[2, 1] == .9
    np.testing.assert_array_equal(actual.electronics_difference_psd_A2_Hz, [3e-26, 0., 0.])
    np.testing.assert_allclose(actual.analysis_axis.frequency_hz, [1e5, 1e6, 4e6], rtol=1e-15)


def test_real_historical_v2_parent_loads_and_preserves_its_declared_detector():
    original = audit.PARENT.read_bytes()
    parent, provenance = audit.load_parent()
    actual = audit.default_detector(parent)
    assert provenance["sha256"] == hashlib.sha256(original).hexdigest()
    assert parent["schema"] == "gabes.stationary_spatial_field_audit.v2"
    assert parent["expected_controls_passed"]
    np.testing.assert_array_equal(actual.transmissions, parent["detector_transmissions"])
    np.testing.assert_allclose(actual.analysis_axis.frequency_hz, parent["detector_frequencies_Hz"], rtol=1e-15)
    assert audit.PARENT.read_bytes() == original


@pytest.mark.parametrize("conflict", ["report", "plot", "same_path"])
def test_cli_checks_immutable_output_conflicts_before_reading_or_computing_parent(tmp_path, monkeypatch, conflict):
    target, plot = tmp_path/"report.json", tmp_path/"plot.png"
    def forbidden(*args, **kwargs):
        raise AssertionError("artifact conflicts must be checked before any expensive build")
    monkeypatch.setattr(audit, "build_report", forbidden)
    if conflict == "same_path":
        plot = target
    else:
        (target if conflict == "report" else plot).write_bytes(b"prior immutable evidence")
    args = ["--parent", str(tmp_path/"absent_parent.json"), "--output", str(target), "--plot", str(plot)]
    with pytest.raises(ValueError if conflict == "same_path" else FileExistsError):
        audit.main(args)
    if conflict != "same_path":
        assert (target if conflict == "report" else plot).read_bytes() == b"prior immutable evidence"
