"""Finite-seed mean-state diagnostics, with independent dense Floquet checks."""

import numpy as np
import pytest

from analysis.grand_challenge.normalization_audit import conditional_inputs
from gabes import constants as c, core, observables
from gabes.fwm_quantum.field import reduced_field_pair, reduced_readout_operators
from gabes.fwm_quantum.model import reduced_pump_noise
from gabes.fwm_quantum.normalization import reduced_dipoles, optical_carriers
from gabes.fwm_quantum.seed_validity import seed_harmonic_hamiltonian, local_seed_backaction
from gabes.quantum.contracts import AnalysisFrequencyAxis


def fixture():
    inputs = conditional_inputs()
    dipoles = reduced_dipoles('uniform-zeeman-rms')
    atom = reduced_pump_noise(dipoles.pump_rabi_rad_s(inputs.pump_power_W, inputs.pump_waist_m), inputs.one_photon_rad_s,
        transit_rate_s_inverse=inputs.transit_rate_s_inverse, transition_scales=dipoles.transition_scales)
    omega = optical_carriers(inputs.detunings)[1:]
    beta = np.array([np.sqrt(inputs.seed_power_W/(c.HBAR*omega[0])), .2j*np.sqrt(inputs.seed_power_W/(c.HBAR*omega[1]))])
    v = seed_harmonic_hamiltonian(dipoles, omega, inputs.uniform_area_m2, beta)
    return atom, v, -c.OMEGA_HF+inputs.two_photon_rad_s


def test_local_floquet_matches_dense_and_converges_with_order():
    atom, v, beat = fixture()
    last = None
    for nf in (2, 3, 4):
        result = local_seed_backaction(atom, v, beat, n_f=nf)
        dense = core.floquet_solve_direct(atom.generator, core.comm_super(v), core.comm_super(v.conj().T), beat,
            [0.], np.zeros((16, 16)), 4, n_f=nf, return_harmonics=True)[0]
        np.testing.assert_allclose(result['harmonics'], dense, rtol=1e-10, atol=1e-13)
        assert result['diagnostics']['phase_sampled_minimum_state_eigenvalue'] > 0
        assert result['diagnostics']['cancellation_aware_residual_passed']
        if last is not None:
            assert np.linalg.norm(result['harmonics'][1:-1]-last) < 1e-10
        last = result['harmonics']


def test_seed_zero_limit_and_perturbative_power_scaling():
    atom, v, beat = fixture()
    zero = local_seed_backaction(atom, v*0, beat, n_f=2)
    np.testing.assert_allclose(zero['harmonics'][2], atom.stationary_state, atol=1e-13)
    np.testing.assert_allclose(zero['harmonics'][[0, 1, 3, 4]], 0, atol=1e-14)
    small = local_seed_backaction(atom, v*.1, beat, n_f=2)
    large = local_seed_backaction(atom, v*.2, beat, n_f=2)
    assert large['mean_state_trace_distance']/small['mean_state_trace_distance'] == pytest.approx(4, rel=1e-4)
    assert large['first_harmonic_relative_error']/small['first_harmonic_relative_error'] == pytest.approx(4, rel=1e-4)
    assert local_seed_backaction(atom, v*1e-5, beat, n_f=2)['diagnostics']['cancellation_aware_residual_passed']
    with pytest.raises(ValueError):
        local_seed_backaction(atom, v, 0, n_f=2)


def test_first_order_floquet_drive_matches_microscopic_field_drift_and_phase_rotation():
    inputs = conditional_inputs()
    dipoles = reduced_dipoles('uniform-zeeman-rms')
    atom, _, beat = fixture()
    omega = optical_carriers(inputs.detunings)[1:]
    beta = np.array([2e4+1e4j, -1e4j])
    v = seed_harmonic_hamiltonian(dipoles, omega, inputs.uniform_area_m2, beta)
    response = local_seed_backaction(atom, v, beat, n_f=3)
    field = reduced_field_pair(atom, inputs.detunings, AnalysisFrequencyAxis.from_hz([0.]),
        number_density_m3=inputs.number_density_m3, uniform_area_m2=inputs.uniform_area_m2,
        optical_omega_rad_s=omega, effective_dipole_C_m=dipoles.base_dipole_C_m,
        transition_scales=dipoles.transition_scales)
    g = dipoles.base_dipole_C_m*np.diag(observables.photon_flux_mode_matrix(*omega, inputs.uniform_area_m2, inputs.uniform_area_m2))/(2*c.HBAR)
    pol = np.einsum('iab,ba->i', reduced_readout_operators(dipoles.transition_scales), response['weak_harmonic'])
    expected = -1j*np.array([1, -1])*g*inputs.number_density_m3*inputs.uniform_area_m2*pol
    actual = field.main.drift[0]@np.array([beta[0], beta[1].conjugate()])
    np.testing.assert_allclose(actual, expected, rtol=1e-10, atol=1e-7)
    phase = .7
    rotated = local_seed_backaction(atom, v*np.exp(1j*phase), beat, n_f=3)
    np.testing.assert_allclose(rotated['harmonics'], response['harmonics']*np.exp(1j*np.arange(-3, 4)*phase)[:, None, None], atol=1e-13)
