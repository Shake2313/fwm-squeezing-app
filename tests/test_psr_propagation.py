"""Longitudinal pump/noise propagation against independent physical limits."""
from dataclasses import replace

import numpy as np
import pytest

from gabes import beam, constants, psr_propagation as propagation
from gabes.psr_conventions import jones_from_angles
from gabes.psr_full_obe import solve_psr_pump
from gabes.psr_homodyne import psr_output_homodyne_readout
from gabes.psr_spectrum import psr_local_generator, psr_uniform_transfer

RF = np.array([-4e5, 4e5])
COMMON = dict(intensity_mw_cm2=12.5, density_m3=1e17, ground_relaxation_hz=1e5)


def test_frozen_uniform_pump_reproduces_s2_without_extra_loss_noise():
    pump = solve_psr_pump(400., **COMMON)
    local = psr_local_generator(pump, RF, method="optimized")
    reference = psr_uniform_transfer(local, pump, .025)
    result = propagation.solve_psr_propagation(400., RF, .025, freeze_pump=True, **COMMON)
    assert result["convergence"]["passed"]
    for key in ("transfer", "noise_greater", "noise_lesser"):
        np.testing.assert_allclose(getattr(result["transfer"], key), getattr(reference, key), rtol=2e-13, atol=2e-13)
    np.testing.assert_allclose(result["profile"]["intensity_mw_cm2"], [12.5, 12.5], rtol=0, atol=3e-15)
    assert result["segments"] == 1
    assert result["transfer"].audit()["passed"]


@pytest.mark.parametrize("length,density", [(0., 1e17), (.025, 0.)])
def test_zero_length_or_density_preserves_pump_and_vacuum(length, density):
    jones = jones_from_angles(.31, .17)*np.exp(.47j)
    result = propagation.solve_psr_propagation(400., RF, length, jones_xy=jones,
                                              **dict(COMMON, density_m3=density))
    assert result["convergence"]["passed"]
    np.testing.assert_allclose(result["output_jones_xy"], jones, atol=1e-14)
    assert result["output_intensity_mw_cm2"] == pytest.approx(12.5, abs=1e-13)
    noise = psr_output_homodyne_readout(result["transfer"], result["output_jones_xy"])
    np.testing.assert_allclose(noise["minimum_db"], 0., atol=1e-12)
    np.testing.assert_allclose(noise["maximum_db"], 0., atol=1e-12)


@pytest.fixture(scope="module")
def elliptic_case():
    return propagation.solve_psr_propagation(400., RF, .025, jones_xy=jones_from_angles(.31, .17), **COMMON)


def test_actual_depletion_jones_change_and_longitudinal_convergence(elliptic_case):
    result = elliptic_case
    assert result["convergence"]["passed"]
    assert result["segments"] >= 8
    assert result["convergence"]["history"][-1]["passed"]
    assert all(np.diff(result["profile"]["intensity_mw_cm2"]) < 0)
    assert 7.69 < result["output_intensity_mw_cm2"] < 7.70
    incident = jones_from_angles(.31, .17)
    assert abs(np.vdot(incident, result["output_jones_xy"]))**2 < .99
    assert .55 < result["carrier_phase_advance_rad"] < .57
    # Independent DOP853 physical-laboratory field ODE, rtol=2e-9.
    assert result["output_intensity_mw_cm2"] == pytest.approx(7.69358247976, rel=7e-6)
    lab_field = (beam.field_amplitude(10*result["output_intensity_mw_cm2"])*result["output_jones_xy"]
                 *np.exp(1j*result["carrier_phase_advance_rad"]))
    reference_field = np.array([207.56757748949113+109.59162846165461j,
                                14.391970168506031+51.63682559357033j])
    np.testing.assert_allclose(lab_field, reference_field, rtol=5e-5)
    noise = psr_output_homodyne_readout(result["transfer"],result["output_jones_xy"],detection_efficiency=.9025)
    np.testing.assert_allclose([noise["minimum_db"][0],noise["maximum_db"][0]],
                               [.6121504807454617,1.124004845801523],rtol=0,atol=1e-3)
    assert len(result["segment_diagnostics"]) == result["segments"]
    assert result["segment_diagnostics"][0]["intensity_mw_cm2"] > result["segment_diagnostics"][-1]["intensity_mw_cm2"]
    assert result["transfer"].audit()["passed"]


def test_equal_column_density_is_equal_in_this_independent_atom_model(elliptic_case):
    scaled = propagation.solve_psr_propagation(400., RF, .0125, jones_xy=jones_from_angles(.31, .17),
                                              **dict(COMMON, density_m3=2e17))
    assert scaled["convergence"]["passed"]
    assert scaled["segments"] == elliptic_case["segments"]
    for key in ("transfer", "noise_greater", "noise_lesser"):
        np.testing.assert_allclose(getattr(scaled["transfer"], key), getattr(elliptic_case["transfer"], key),
                                   rtol=2e-11, atol=2e-11)
    np.testing.assert_allclose(scaled["output_jones_xy"], elliptic_case["output_jones_xy"], atol=1e-13)
    assert scaled["output_intensity_mw_cm2"] == pytest.approx(elliptic_case["output_intensity_mw_cm2"], rel=1e-13)


def test_thin_cell_restores_uniform_limit():
    true = propagation.solve_psr_propagation(400., RF, 1e-6, **COMMON)
    fixed = propagation.solve_psr_propagation(400., RF, 1e-6, freeze_pump=True, **COMMON)
    assert true["convergence"]["passed"]
    assert true["output_intensity_mw_cm2"] < 12.5
    for key in ("transfer", "noise_greater", "noise_lesser"):
        np.testing.assert_allclose(getattr(true["transfer"], key), getattr(fixed["transfer"], key), atol=3e-8, rtol=1e-4)


def test_analytic_passive_absorber_preserves_vacuum_and_updates_pump(monkeypatch):
    base = solve_psr_pump(400., **COMMON)
    template = psr_local_generator(base, RF, method="optimized")
    alpha, beta = 40., 17.
    k = 2*np.pi*base["conditions"]["optical_frequency_hz"]/constants.C_LIGHT
    def pump(detuning, *, intensity_mw_cm2, jones_xy, **kwargs):
        field = beam.field_amplitude(10*intensity_mw_cm2)
        return dict(base, conditions=dict(base["conditions"], intensity_mw_cm2=intensity_mw_cm2,
                                          peak_field_v_m=field, jones_xy=np.asarray(jones_xy)),
                    chi_probe=(2*beta+1j*alpha)/k, alpha_per_m=alpha,
                    polarization_xy_c_m2=constants.EPS_0*field*(2*beta+1j*alpha)/k*np.asarray(jones_xy))
    def local(*args, **kwargs):
        drift=np.broadcast_to(-alpha/2*np.eye(4)+1j*beta*np.diag(template.signs),(2,4,4))
        greater=np.broadcast_to(alpha*np.diag((1+template.signs)/2),(1,2,4,4))
        lesser=np.broadcast_to(alpha*np.diag((1-template.signs)/2),(1,2,4,4))
        return replace(template,drift=drift,noise_greater_by_reservoir=greater,
                       noise_lesser_by_reservoir=lesser,reservoir_names=("passive",))
    monkeypatch.setattr(propagation,"solve_psr_pump",pump)
    monkeypatch.setattr(propagation,"psr_local_generator",local)
    result=propagation.solve_psr_propagation(400.,RF,.03,max_segments=128,rtol=1e-4,**COMMON)
    assert result["convergence"]["passed"]
    assert result["output_intensity_mw_cm2"]==pytest.approx(12.5*np.exp(-alpha*.03),rel=3e-5)
    assert result["carrier_phase_advance_rad"]==pytest.approx(beta*.03,abs=1e-13)
    noise=psr_output_homodyne_readout(result["transfer"],result["output_jones_xy"])
    np.testing.assert_allclose(noise["minimum_db"],0.,atol=1e-11)
    np.testing.assert_allclose(noise["maximum_db"],0.,atol=1e-11)
    np.testing.assert_allclose(result["transfer"].transfer,np.broadcast_to(np.eye(4)*np.exp(-alpha*.03/2),(2,4,4)),atol=1e-12)


def test_budget_and_unconverged_grid_are_explicit():
    one = propagation.solve_psr_propagation(400., RF, .025, segments=1, max_segments=1, **COMMON)
    assert not one["convergence"]["passed"]
    budget = propagation.solve_psr_propagation(400., RF, .025, max_atomic_evaluations=2, **COMMON)
    assert not budget["convergence"]["passed"]
    assert "budget" in budget["convergence"]["status"]
    assert budget["performance"]["atomic_evaluations"] <= 2
    frozen = propagation.solve_psr_propagation(400., RF, .025, max_atomic_evaluations=1, freeze_pump=True, **COMMON)
    assert frozen["convergence"]["passed"]


@pytest.mark.parametrize("changes", [dict(intensity_mw_cm2=0.),dict(segments=0),dict(max_segments=1),
                                      dict(velocity_model="unknown"),dict(jones_xy=(2.,0.)),
                                      dict(thermal_options={"reference":True})])
def test_invalid_propagation_inputs(changes):
    with pytest.raises(ValueError):
        propagation.solve_psr_propagation(400., RF, .025, **dict(COMMON, **changes))


def test_thermal_local_recomputation_and_accumulated_velocity_checks(monkeypatch):
    from gabes import magneto_multilevel
    atom = magneto_multilevel.build_multilevel_atom("Rb87", "D1")
    shifts = np.array([-30.,80.])*atom.frequency_hz/constants.C_LIGHT/1e6
    monkeypatch.setattr(magneto_multilevel,"velocity_quadrature",
                        lambda *args,**kwargs:(shifts,np.array([.25,.75]),200.))
    result = propagation.solve_psr_propagation(400.,RF,.01,velocity_model="thermal",segments=1,max_segments=8,
        jones_xy=jones_from_angles(.31,.17),**dict(COMMON,density_m3=1e16))
    assert result["convergence"]["passed"]
    assert result["velocity_convergence"]["passed"]
    assert len(result["velocity_convergence"]["checks"])==2
    assert result["performance"]["atomic_evaluations"]==1+7*result["performance"]["local_evaluations"]
    assert "velocity_m_s" not in result["conditions"]
    assert all(d["velocity_convergence"]["converged"] for d in result["segment_diagnostics"])
    assert result["output_intensity_mw_cm2"]<12.5
    exhausted = propagation.solve_psr_propagation(400.,RF,.01,velocity_model="thermal",max_atomic_evaluations=8,**COMMON)
    assert not exhausted["convergence"]["passed"]
    assert "budget" in exhausted["convergence"]["status"]
    assert exhausted["performance"]["atomic_evaluations"]==8
    assert not exhausted["velocity_convergence"]["passed"]


def test_accumulated_velocity_error_can_reject_locally_accepted_kernels(monkeypatch):
    from gabes import magneto_multilevel
    atom = magneto_multilevel.build_multilevel_atom("Rb87", "D1")
    shifts = np.array([-30.,80.])*atom.frequency_hz/constants.C_LIGHT/1e6
    monkeypatch.setattr(magneto_multilevel,"velocity_quadrature",
                        lambda *args,**kwargs:(shifts,np.array([.25,.75]),200.))
    original = propagation.solve_psr_thermal
    def inaccurate_history(*args, **kwargs):
        result = original(*args, **kwargs)
        assert result["convergence"]["converged"]
        # Model a local error estimator which accepted a materially different
        # previous-order generator. Every individual generator stays physical.
        old = result["history"][-3]["local"]
        result["history"][-3]["local"] = replace(old,drift=1.5*old.drift,
            noise_greater_by_reservoir=1.5*old.noise_greater_by_reservoir,
            noise_lesser_by_reservoir=1.5*old.noise_lesser_by_reservoir)
        return result
    monkeypatch.setattr(propagation,"solve_psr_thermal",inaccurate_history)
    result = propagation.solve_psr_propagation(400.,RF,.01,velocity_model="thermal",segments=1,max_segments=2,**COMMON)
    assert not result["convergence"]["passed"]
    assert not result["velocity_convergence"]["passed"]
    assert "accumulated velocity" in result["convergence"]["status"]
    assert result["transfer"] is not None  # Diagnostic only; acceptance is false.
