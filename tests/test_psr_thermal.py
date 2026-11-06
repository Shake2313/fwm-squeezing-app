"""Thermal PSR uses one common drift/noise measure before propagation."""
from copy import deepcopy

import numpy as np
import pytest

from gabes import constants, magneto_multilevel
from gabes.psr_full_obe import solve_psr_pump
from gabes.psr_spectrum import psr_local_generator, psr_uniform_transfer
from gabes.psr_homodyne import psr_homodyne_readout
from gabes.psr_thermal import solve_psr_thermal, thermal_convergence, _kernel_difference


RF = np.array([-4e5, 4e5])
INPUTS = dict(intensity_mw_cm2=12.5, temp_c=70., ground_relaxation_hz=1e5, density_m3=1e16)


def test_shared_grid_default_and_explicit_options():
    atom = magneto_multilevel.build_multilevel_atom("Rb87", "D1")
    laser = magneto_multilevel.component_centroid_mhz(atom)+400.
    original = magneto_multilevel.velocity_quadrature(atom, laser, 343.15, order=4)
    explicit = magneto_multilevel.velocity_quadrature(
        atom, laser, 343.15, order=4, cutoff_sigma=6., extra_shifts_mhz=())
    for first, second in zip(original, explicit):
        np.testing.assert_array_equal(first, second)
    assert len(original[0]) == 120
    assert np.sum(original[1]) == pytest.approx(1., abs=1e-15)
    # Wider cutoff changes the domain without changing physical temperature.
    widened = magneto_multilevel.velocity_quadrature(
        atom, laser, 343.15, order=4, cutoff_sigma=7., extra_shifts_mhz=(-20., 20.))
    assert widened[0][0] < original[0][0] and widened[0][-1] > original[0][-1]
    assert widened[2] == original[2]


def test_two_class_linear_drift_noise_and_mean_field(monkeypatch):
    atom = magneto_multilevel.build_multilevel_atom("Rb87", "D1")
    velocities = np.array([-30., 80.])
    shifts = velocities*atom.frequency_hz/constants.C_LIGHT/1e6
    weights = np.array([.25, .75])
    monkeypatch.setattr(magneto_multilevel, "velocity_quadrature",
                        lambda *a, **kw: (shifts, weights, 200.))
    result = solve_psr_thermal(400., RF, **INPUTS, velocity_order=2, max_order=2, reference=True)
    states = [solve_psr_pump(400., velocity_m_s=v, **INPUTS) for v in velocities]
    local = [psr_local_generator(state, RF, method="direct") for state in states]
    for name in ("drift", "noise_greater_by_reservoir", "noise_lesser_by_reservoir"):
        expected = sum(w*getattr(item, name) for w, item in zip(weights, local))
        np.testing.assert_allclose(getattr(result["local"], name), expected, rtol=2e-11, atol=1e-12)
    assert result["pump_state"]["chi_probe"] == pytest.approx(sum(w*s["chi_probe"] for w, s in zip(weights, states)))
    assert result["pump_state"]["alpha_per_m"] == pytest.approx(sum(w*s["alpha_per_m"] for w, s in zip(weights, states)))
    assert result["local"].audit()["passed"]
    assert "rho" not in result["pump_state"] and "generator_s_inverse" not in result["pump_state"]
    assert "velocity_m_s" not in result["pump_state"]["conditions"]
    assert "effective_detuning_mhz" not in result["pump_state"]["conditions"]
    assert result["pump_state"]["model_id"] != states[0]["model_id"]
    assert not thermal_convergence(result, .025)["passed"]


def test_zero_thermal_width_recovers_single_velocity_and_passive_vacuum():
    # Positive T approaching zero, with density independently specified.
    inputs = dict(INPUTS, temp_c=-273.15+1e-9, intensity_mw_cm2=0.)
    result = solve_psr_thermal(400., RF, **inputs, velocity_order=2, max_order=4, rtol=.001)
    state = solve_psr_pump(400., **inputs)
    direct = psr_local_generator(state, RF, method="direct")
    np.testing.assert_allclose(result["local"].drift, direct.drift, rtol=2e-8, atol=1e-10)
    np.testing.assert_allclose(result["local"].noise_greater_by_reservoir,
                               direct.noise_greater_by_reservoir, rtol=2e-8, atol=1e-10)
    assert result["convergence"]["converged"]
    assert thermal_convergence(result, .025)["passed"]
    noise = psr_homodyne_readout(psr_uniform_transfer(result["local"], result["pump_state"], .025))
    np.testing.assert_allclose(noise["minimum_db"], 0., atol=2e-11)
    np.testing.assert_allclose(noise["maximum_db"], 0., atol=2e-11)


def test_explicit_budget_and_noise_sensitive_refinement():
    result = solve_psr_thermal(400., RF, **INPUTS, velocity_order=2, max_order=10, max_classes=150)
    assert result["performance"]["class_evaluations"] <= 150
    assert not result["convergence"]["converged"]
    assert "budget" in result["convergence"]["status"]
    assert not thermal_convergence(result, .025)["passed"]
    # A noise-only change cannot hide behind unchanged mean susceptibility.
    before = result["history"][0]
    after = deepcopy(before)
    after["local"].noise_greater_by_reservoir[0, 0, 1, 1] += .1
    check = _kernel_difference(before, after, rtol=.005, atol=1e-6)
    assert not check["passed"]
    assert check["checks"]["mean_polarization"]["max_absolute"] == 0.


def test_propagated_noise_convergence_cannot_be_overruled_by_kernel_flag(monkeypatch):
    atom = magneto_multilevel.build_multilevel_atom("Rb87", "D1")
    monkeypatch.setattr(magneto_multilevel, "velocity_quadrature",
                        lambda *a, **kw: (np.array([0.]), np.array([1.]), 1.))
    result = solve_psr_thermal(400., RF, **INPUTS, velocity_order=2, max_order=4)
    assert thermal_convergence(result, .025)["passed"]
    original = psr_homodyne_readout
    count = 0
    def mismatched(*args, **kwargs):
        nonlocal count
        out = original(*args, **kwargs)
        count += 1
        if count > 1:
            out["maximum_db"] = out["maximum_db"]+1.
        return out
    monkeypatch.setattr("gabes.psr_thermal.psr_homodyne_readout", mismatched)
    checked = thermal_convergence(result, .025)
    assert not checked["passed"]
    assert "extrema" in checked["status"]


@pytest.mark.parametrize("kwargs", [dict(rtol=0), dict(atol_per_m=0), dict(cutoff_sigma=7),
                                      dict(max_classes=0), dict(velocity_order=True),
                                      dict(velocity_order=6, max_order=4), dict(method="other")])
def test_invalid_quadrature_controls_fail_early(kwargs):
    with pytest.raises(ValueError):
        solve_psr_thermal(400., RF, **INPUTS, **kwargs)


def test_unsymmetric_rf_and_too_small_budget_rejected():
    with pytest.raises(ValueError, match="symmetric"):
        solve_psr_thermal(400., [4e5], **INPUTS)
    with pytest.raises(ValueError, match="budget"):
        solve_psr_thermal(400., RF, **INPUTS, max_classes=1)


def test_room_temperature_refinement_tail_and_independent_reference():
    result = solve_psr_thermal(400., RF, **INPUTS)
    audit = thermal_convergence(result, .025, detection_efficiency=.9025)
    assert audit["passed"]
    assert audit["order"]["passed"] and audit["cutoff"]["passed"]
    assert audit["omitted_maxwell_probability"] < 3e-12
    assert {item["detection_efficiency"] for item in audit["output_checks"]} == {1., .9025}
    noise = psr_homodyne_readout(
        psr_uniform_transfer(result["local"], result["pump_state"], .025), detection_efficiency=.9025)
    # Independently evaluated original optical-only grid, interval order 20,
    # direct jump sums, at the same pump and detector conditions.
    np.testing.assert_allclose(noise["minimum_db"], .02499930, atol=1e-4)
    np.testing.assert_allclose(noise["maximum_db"], .14212376, atol=1e-4)
    assert result["performance"]["retained_integral_bytes"] < 100_000


def test_elliptic_two_class_mean_polarization_and_kernels(monkeypatch):
    from gabes.psr_conventions import jones_from_angles

    atom = magneto_multilevel.build_multilevel_atom("Rb87", "D1")
    velocities = np.array([-35., 65.])
    weights = np.array([.3, .7])
    shifts = velocities*atom.frequency_hz/constants.C_LIGHT/1e6
    monkeypatch.setattr(magneto_multilevel, "velocity_quadrature",
                        lambda *a, **kw: (shifts, weights, 200.))
    jones = np.exp(.39j)*jones_from_angles(.23, .14)
    result = solve_psr_thermal(400., RF, **INPUTS, jones_xy=jones,
                               velocity_order=2, max_order=4)
    states = [solve_psr_pump(400., velocity_m_s=v, **INPUTS, jones_xy=jones) for v in velocities]
    local = [psr_local_generator(state, RF, method="direct") for state in states]
    for key in ("polarization_xy_c_m2", "chi_probe", "alpha_per_m"):
        expected = sum(w*state[key] for w, state in zip(weights, states))
        np.testing.assert_allclose(result["pump_state"][key], expected, rtol=2e-12, atol=0.)
    for key in ("drift", "noise_greater_by_reservoir", "noise_lesser_by_reservoir"):
        np.testing.assert_allclose(getattr(result["local"], key),
            sum(w*getattr(part, key) for w, part in zip(weights, local)), rtol=2e-11, atol=1e-12)
    np.testing.assert_array_equal(result["pump_state"]["conditions"]["jones_xy"], jones)
    assert result["convergence"]["converged"]
    assert "polarization_vector" in result["convergence"]["order"]["checks"]

    # A transverse mean-polarization error must not be hidden by an unchanged
    # longitudinal susceptibility or field noise kernel.
    before = result["history"][0]
    after = deepcopy(before)
    field = after["pump_state"]["conditions"]["peak_field_v_m"]
    after["pump_state"]["polarization_xy_c_m2"][1] += constants.EPS_0*field*1e-6
    audit = _kernel_difference(before, after, rtol=.005, atol=1e-6)
    assert not audit["passed"]
    assert audit["checks"]["mean_polarization"]["max_absolute"] == 0.
    assert audit["checks"]["polarization_vector"]["max_scaled"] > 1.


def test_thermal_invalid_jones_fails_before_class_integration(monkeypatch):
    def forbidden(*args, **kwargs):
        raise AssertionError("invalid Jones input must fail before velocity integration")
    monkeypatch.setattr(magneto_multilevel, "velocity_quadrature", forbidden)
    with pytest.raises(ValueError, match="Jones"):
        solve_psr_thermal(400., RF, **INPUTS, jones_xy=(1., 1.))
