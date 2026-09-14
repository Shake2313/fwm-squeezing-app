"""Independent angular sums, dimensional reciprocity and consumed-input evidence."""

from dataclasses import replace
import json

import numpy as np
import pytest

from analysis.grand_challenge.normalization_audit import conditional_inputs, main
from analysis.grand_challenge.reference.uncoupled_d1 import d1_dipole_reference
from gabes import constants as c, observables, species
from gabes.fwm_quantum.field import reduced_readout_operators
from gabes.fwm_quantum.inputs import audit_consumed_inputs, power_normalized_readout
from gabes.fwm_quantum.normalization import (
    d1_dipole_from_decay, gaussian_peak_field, normalization_audit,
    optical_carriers, reduced_dipoles, validated_transition_scales,
)
from gabes.quantum.contracts import AnalysisFrequencyAxis, ParameterEvidence
from gabes.quantum.readout import DetectorResponse


def test_dipole_decay_convention_and_gaussian_power_integral():
    d = d1_dipole_from_decay(c.GAMMA, c.NU_D1_85RB)
    gamma = c.OMEGA_D1**3*d*d/(3*np.pi*c.EPS_0*c.HBAR*c.C_LIGHT**3)
    assert gamma == pytest.approx(c.GAMMA, rel=3e-15)
    E = gaussian_peak_field(.6, 530e-6)
    assert c.EPS_0*c.C_LIGHT*E**2/2*np.pi*(530e-6)**2/2 == pytest.approx(.6)
    assert gaussian_peak_field(0, 1e-3) == 0
    assert gaussian_peak_field(4*.6, 2*530e-6) == pytest.approx(E)


@pytest.mark.parametrize("q", [0, 1, 2])
def test_uncoupled_electronic_reference_matches_manifold_average_for_each_polarization(q):
    ref = d1_dipole_reference()
    assert ref["basis_unitarity_error"] < 1e-14
    assert ref["emission_closure_error"] < 1e-14
    coupling = reduced_dipoles("uniform-zeeman-rms")
    for row in ref["transition_rows"]:
        g, e = row["Fg"]-2, row["Fe"]
        assert coupling.transition_scales[g, e]**2 == pytest.approx(row["mean_strength_by_q"][q], rel=2e-14)
    # Each normalized, unpolarized manifold couples with total dJ^2/3.
    np.testing.assert_allclose(np.sum(coupling.transition_scales[:2, 2:]**2, axis=1), 1/3, atol=1e-15)


def test_no_single_scalar_after_sublevel_sum_can_match_both_manifolds():
    audit = normalization_audit()
    assert audit["pump_to_weak_dipole_ratio"] == pytest.approx(1.99937341410737)
    for row in audit["transition_rows"]:
        assert row["sum_rule_error"] < 1e-14
        assert row["historical_to_uniform_strength_ratio"] == pytest.approx((2*row["Fg"]+1)/12)
    assert not audit["single_scalar_matches_both_manifolds"]
    assert 1/5 != 1/7


@pytest.mark.parametrize("convention", ["legacy-reciprocal", "uniform-zeeman-rms"])
def test_same_dipole_drives_pump_and_photon_flux_field(convention):
    coupling = reduced_dipoles(convention)
    area, power = 1.2e-7, .6
    # A=pi*w^2/2 gives the same local E at a Gaussian pump peak.
    waist = np.sqrt(2*area/np.pi)
    Q = observables.photon_flux_mode_matrix(c.OMEGA_D1, c.OMEGA_D1, area, area)[0, 0]
    photon_flux = power/(c.HBAR*c.OMEGA_D1)
    g = coupling.base_dipole_C_m*Q/(2*c.HBAR)
    assert 2*g*np.sqrt(photon_flux) == pytest.approx(coupling.pump_rabi_rad_s(power, waist))
    ops = reduced_readout_operators(coupling.transition_scales)
    np.testing.assert_array_equal(ops[0, 1, 2:], coupling.transition_scales[1, 2:])
    np.testing.assert_array_equal(ops[1, 2:, 0], coupling.transition_scales[0, 2:])


def test_carrier_anchor_uses_same_hyperfine_energies_as_hamiltonian():
    inputs = conditional_inputs()
    pump, probe, conjugate = optical_carriers(inputs.detunings)
    # Independent Casimir energies using the exact stored splittings.
    eg = species.hf_energy_mhz(c.NU_GROUND_HF/3e6, 0, 2.5, .5, 2)*1e6
    ee = species.hf_energy_mhz(c.NU_EXCITED_HF_D1/3e6, 0, 2.5, .5, 3)*1e6
    assert (pump-c.OMEGA_D1-inputs.one_photon_rad_s)/(2*np.pi) == pytest.approx(ee-eg, abs=.1)
    assert (probe+conjugate-2*pump) == pytest.approx(0, abs=1.)
    assert (probe-pump)/(2*np.pi) == pytest.approx(-c.NU_GROUND_HF-8e6, abs=.1)


def _detector():
    return DetectorResponse(AnalysisFrequencyAxis.from_hz([1e6]), [.85, .85], [[1., 1.]], 1., [0.], "fixture")


def test_evidence_is_bound_to_actual_small_values_units_and_heldout_ids():
    ledger = conditional_inputs().consumed_inputs(_detector(), convention="legacy-reciprocal")
    evidence = [ParameterEvidence(k, v, unit, "independent", "synthetic test record", 0.,
                    "test only", ("calibration",), ("test observable",), "unit test")
                for k, (v, unit) in ledger["scalars"].items()]
    assert audit_consumed_inputs(ledger, evidence).passed
    k = next(i for i, e in enumerate(evidence) if e.parameter_id == "stored_D1_dipole")
    bad = list(evidence)
    bad[k] = replace(bad[k], value=bad[k].value*2)
    assert not audit_consumed_inputs(ledger, bad).passed
    bad[k] = replace(evidence[k], unit="debye")
    assert not audit_consumed_inputs(ledger, bad).passed
    assert not audit_consumed_inputs(ledger, evidence, target_dataset_ids=["calibration"]).passed
    assert not audit_consumed_inputs(ledger, evidence[:-1]).passed
    assert not audit_consumed_inputs(ledger, [replace(e, status="assumed") for e in evidence]).passed


def test_power_adapter_vacuum_density_limit_and_input_validation():
    inputs = replace(conditional_inputs(), number_density_m3=0.)
    rf = AnalysisFrequencyAxis.from_hz([-1e6, 0., 1e6])
    result = power_normalized_readout(inputs, rf, _detector())
    assert result.probe_power_gain == pytest.approx(1.)
    np.testing.assert_allclose(result.spectrum.quantum_ratio, 1., atol=1e-12)
    with pytest.raises(ValueError):
        replace(inputs, pump_waist_m=0)
    with pytest.raises(ValueError):
        validated_transition_scales(np.ones((4, 4)))
    with pytest.raises(ValueError):
        validated_transition_scales(np.zeros((2, 2)))
    coupling = reduced_dipoles("uniform-zeeman-rms")
    with pytest.raises(ValueError):
        coupling.transition_scales.setflags(write=True)


def test_normalization_artifact_checks_joint_readout_and_refuses_overwrite(tmp_path):
    path = tmp_path/"normalization.json"
    assert main(["--output", str(path)]) == 0
    original = path.read_bytes()
    report = json.loads(original)
    assert report["expected_controls_passed"]
    assert not report["independent_input_audit"]["passed"]
    assert not report["absolute_hot_vapor_prediction"]
    assert all(case["passed"] for case in report["configurations"])
    # Pin the already-audited historical fixture while new paths stay explicit.
    old = report["configurations"][0]
    assert old["probe_power_gain"] == pytest.approx(1.0517029257110075, rel=1e-10)
    assert old["detected_S_minus_db"][0] == pytest.approx(-.3658522796429018, abs=1e-9)
    with pytest.raises(FileExistsError):
        main(["--output", str(path)])
    assert path.read_bytes() == original
