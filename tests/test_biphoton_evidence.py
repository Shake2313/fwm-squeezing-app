"""Source anchors, numerical scope and quantum claims stay separate in SFWM."""
import json
from types import SimpleNamespace
from unittest.mock import patch

import numpy as np
import pytest

from gabes.biphoton_evidence import biphoton_output_evidence
from gabes import observables
from gabes.schemes import fwm


def _evidence(*, reference=True, predictive=True, converged=True, filter_mhz=300.0):
    topology = SimpleNamespace(
        name="reference" if reference else "generic diamond",
        reference_od=112.0 if reference else None,
        reference_fwhm_ns=0.56 if reference else None,
        reference_bandwidth_mhz=300.0 if reference else None,
        target_g2_peak=44.0 if reference else None)
    return biphoton_output_evidence(
        {"topology": topology, "predictive": predictive,
         "velocity_converged": converged}, {}, {"filter_bandwidth_mhz": filter_mhz})


def test_template_coefficient_cannot_inherit_a_literature_anchor():
    reference, template = _evidence(), _evidence(reference=False)
    assert reference["observables"]["source_pair_rate_cps"]["origin"] == "reference_anchored"
    assert template["source_model"] == "assumed_template_model"
    assert template["observables"]["source_pair_rate_cps"]["origin"] == "assumed_template"
    assert template["observables"]["pair_rate_cps"]["status"] == "assumed_template"
    for evidence in (reference, template):
        assert not evidence["measured_data_available"]
        assert not evidence["independent_validation"]
        assert not evidence["seeded_fwm_validation_applicable"]
        assert not any(evidence["claims"][key] for key in (
            "physical_nonclassicality_available", "physical_squeezing_available",
            "absolute_rate_prediction_available"))


@pytest.mark.parametrize("reference", [False, True])
def test_fwhm_refinement_never_certifies_whole_waveform_spectrum_or_rate(reference):
    evidence = _evidence(reference=reference)
    assert evidence["numerical_scope"]["velocity_fwhm_status"] == "passed"
    assert evidence["numerical_scope"]["velocity_fwhm_relative_tolerance"] == 0.03
    for key in ("whole_waveform_certified", "time_grid_and_support_certified",
                "spectral_bandwidth_certified", "absolute_rate_certified"):
        assert evidence["numerical_scope"][key] is False
    assert all(record["independent_validation"] is False
               for record in evidence["observables"].values())


def test_velocity_failure_propagates_by_dependency_not_to_unfiltered_anchor():
    outputs = _evidence(converged=False)["observables"]
    for key in ("source_fwhm_ns", "source_bandwidth_mhz", "detected_fwhm_ns",
                "gate_capture_fraction", "g2_peak", "CAR", "coincidence_cps",
                "heralding_signal", "pair_rate_cps", "filter_transmission"):
        assert outputs[key]["status"] == "qualitative_unconverged"
        assert not outputs[key]["display_numeric"]
    assert outputs["source_pair_rate_cps"]["display_numeric"]
    assert outputs["phase_match_weight"]["display_numeric"]
    no_filter = _evidence(converged=False, filter_mhz=0.0)["observables"]
    assert no_filter["pair_rate_cps"]["display_numeric"]
    assert no_filter["raw_accidental_cps"]["display_numeric"]
    assert not no_filter["coincidence_cps"]["display_numeric"]


def test_reference_timing_and_added_background_are_not_validation():
    evidence = _evidence(predictive=False)
    assert evidence["numerical_scope"]["velocity_fwhm_status"] == "not_assessed"
    outputs = evidence["observables"]
    assert outputs["detected_fwhm_ns"]["origin"] == "detector_convolution"
    assert outputs["added_accidental_cps"]["origin"] == "reference_added_background"
    assert outputs["g2_peak"]["origin"] == "reference_anchored"
    assert "target_g2_peak" in outputs["added_accidental_cps"]["dependencies"]
    assert outputs["cauchy_schwarz_R"]["origin"] == "assumed_thermal_arms"
    assert not outputs["physical_cauchy_schwarz_R"]["display_numeric"]


@pytest.mark.parametrize("channel", [fwm.CS_CHANNEL_917, fwm.CS_CHANNEL_795])
def test_cs_reference_model_does_not_turn_fallback_bandwidth_into_a_reference(channel):
    topology = fwm._cs_btw_spec(channel)
    assert topology.reference_fwhm_ns is not None
    assert topology.reference_bandwidth_mhz is None
    outputs = biphoton_output_evidence(
        {"topology": topology, "predictive": False, "velocity_converged": True},
        {}, {})["observables"]
    assert outputs["source_v"]["origin"] == "reduced_model"
    assert outputs["psi_tau"]["origin"] == "reference_anchored"
    assert "emission_decay_ns" in outputs["psi_tau"]["dependencies"]
    assert outputs["source_bandwidth_mhz"]["origin"] == "assumed_input"
    assert outputs["source_bandwidth_mhz"]["status"] == "assumed_input"
    assert "psi_tau" not in outputs["source_bandwidth_mhz"]["dependencies"]
    assert "300 MHz fallback" in " ".join(outputs["source_bandwidth_mhz"]["assumptions"])


def test_reference_bandwidth_and_lifetime_are_separate_input_provenance():
    outputs = _evidence(predictive=False)["observables"]
    assert outputs["source_v"]["origin"] == "reduced_model"
    assert outputs["source_bandwidth_mhz"]["origin"] == "reference_anchored"
    assert outputs["source_bandwidth_mhz"]["dependencies"] == ["reference_bandwidth_mhz"]


def test_raw_correlations_precede_target_background_and_histogram_is_gate_independent():
    topology = SimpleNamespace(name="reference", reference_od=112.0,
                               reference_fwhm_ns=.56, reference_bandwidth_mhz=300.,
                               target_g2_peak=44.0)
    tau = np.linspace(0.0, 2.0, 21)
    wave = np.exp(-tau)
    kwargs = dict(signal_eff=.5, idler_eff=.5, coincidence_window_ns=.5)
    raw_stats = observables.biphoton_stats(tau, wave, 1000., **kwargs)
    calibrated = observables.biphoton_stats(tau, wave, 1000., target_g2_peak=44., **kwargs)
    assert calibrated["g2_peak"] < raw_stats["g2_peak"]
    assert calibrated["raw_g2_peak"] == pytest.approx(raw_stats["g2_peak"], rel=2e-15)
    assert calibrated["raw_CAR"] == pytest.approx(raw_stats["CAR"], rel=2e-15)
    outputs = biphoton_output_evidence(
        {"topology": topology, "predictive": False}, calibrated, kwargs)["observables"]
    for key in ("raw_g2_peak", "raw_CAR"):
        assert outputs[key]["origin"] == "detector_derived"
        assert "target_g2_peak" not in outputs[key]["dependencies"]
        assert "accidental_cps" not in outputs[key]["dependencies"]
    assert "raw_accidental_cps" in outputs["raw_CAR"]["dependencies"]
    assert "target_g2_peak" in outputs["g2_peak"]["dependencies"]
    assert "coincidence_window_ns" in outputs["raw_accidental_cps"]["dependencies"]
    assert "gate_captured_coincidence_cps" in outputs["CAR"]["dependencies"]
    assert "gate_captured_coincidence_cps" not in outputs["g2_peak"]["dependencies"]


@pytest.fixture(scope="module")
def rb_source():
    scheme = fwm.FWMScheme()
    params = dict(scheme.defaults(), mode=fwm.MODE_BIPHOTON,
                  topology=fwm.TOPOLOGY_RB87_TELECOM,
                  biphoton_model=fwm.BIPHOTON_PREDICTIVE)
    params.update(scheme.recommended_defaults(params)[fwm.MODE_BIPHOTON])
    return scheme, params, scheme.compute(params)


def test_readout_has_current_evidence_without_source_resolve_or_array_mutation(rb_source):
    scheme, params, raw = rb_source
    saved = {key: value.copy() for key, value in raw.items() if isinstance(value, np.ndarray)}
    keys = set(raw)
    with patch.object(fwm.GenericFWMSolver, "compute_biphoton", side_effect=AssertionError("source re-solved")):
        base = scheme.headless_observables(raw, params)
        changed = scheme.headless_observables(
            raw, dict(params, signal_eff_pct=40.0, coincidence_window_ns=0.05))
    assert set(raw) == keys
    for key, original in saved.items():
        np.testing.assert_array_equal(raw[key], original)
    assert (base["biphoton_evidence"]["observables"]["source_pair_rate_cps"]
            == changed["biphoton_evidence"]["observables"]["source_pair_rate_cps"])
    assert changed["biphoton_readout"]["coincidence_cps"] != base["biphoton_readout"]["coincidence_cps"]
    stats = base["biphoton_readout"]
    assert stats["thermal_arm_comparison_R"] == stats["cauchy_schwarz_R"]
    assert stats["physical_cauchy_schwarz_R"] is None
    detection = next(t["markdown"] for t in base["tables"] if t["title"] == "Detection estimates")
    assert "| Cauchy-Schwarz R |" not in detection
    assert "assumed g²_SS = g²_II = 2" in detection
    assert "Physical nonclassicality | unavailable" in detection


def test_failed_waveform_check_withholds_affected_heroes_and_preserves_raw_diagnostics(rb_source):
    scheme, params, original = rb_source
    raw = dict(original, velocity_converged=False)
    view = scheme.headless_observables(raw, params)
    metrics = {m["label"]: m for m in view["metrics"]}
    for label in ("Peak signal–idler correlation", "CAR", "Pair rate", "Source bandwidth"):
        assert metrics[label]["value"] == "unconverged"
        assert metrics[label]["kind"] == "status"
    assert np.isfinite(view["biphoton_readout"]["g2_peak"])
    assert np.isfinite(view["biphoton_readout"]["coincidence_cps"])
    literature = next(t["markdown"] for t in view["tables"] if t["title"] == "Literature comparison")
    assert "| g2 peak | unconverged |" in literature
    assert "| Source pair rate / pump |" in literature
    assert "| Filter-passed pair rate / pump | unconverged |" in literature


def test_strict_export_preserves_evidence_and_readout_status(rb_source):
    from gabes_ui.export import build_export_payload
    scheme, params, raw = rb_source
    view = scheme.headless_observables(raw, params)
    payload = build_export_payload(
        scheme=scheme, raw=raw, params=params, view=view,
        plot_payloads=[], comparison_payload=None, cache_version=scheme.cache_version,
        readout_cache_version="evidence-test")
    decoded = json.loads(json.dumps(payload, allow_nan=False))
    assert decoded["view"]["biphoton_evidence"] == view["biphoton_evidence"]
    assert decoded["view"]["biphoton_readout"]["physical_cauchy_schwarz_R"] is None
    assert decoded["view"]["biphoton_readout"]["coincidence_cps"] == view["biphoton_readout"]["coincidence_cps"]
