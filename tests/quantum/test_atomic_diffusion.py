"""Atomic Einstein diffusion, independent QRT and analytic spectral limits."""

from dataclasses import replace

import numpy as np
import pytest

from analysis.grand_challenge.reference.atomic_qrt import (
    adjoint_einstein_diffusion, qrt_ordered_spectrum,
)
from gabes import atoms, constants
from gabes.fwm_quantum.model import reduced_pump_noise, minus_branch_atomic_frequencies
from gabes.quantum.contracts import AnalysisFrequencyAxis, GeneratorFrequencyAxis, OpticalDetunings
from gabes.quantum.diffusion import (
    ordered_jump_diffusion, stationary_atomic_noise, traceless_hermitian_basis,
)
from gabes.quantum.reservoirs import CollapseChannel, ExplicitReservoirs, thermal_reset_channels
from gabes.schemes import fwm


def _two_level(gamma=2.4, resonance=3.7, drive=0.0):
    h = np.array([[0., drive/2], [drive/2, resonance]])
    reservoirs = ExplicitReservoirs(2, (CollapseChannel(
        "emission", [[0., np.sqrt(gamma)], [0., 0.]], "radiative", "two-level fixture"),))
    return h, reservoirs


@pytest.mark.parametrize("n", [2, 3, 4])
def test_complete_basis_closes_all_trace_zero_atomic_operators(n):
    f = traceless_hermitian_basis(n)
    assert f.shape == (n*n-1, n, n)
    np.testing.assert_allclose(np.trace(f, axis1=1, axis2=2), 0, atol=2e-16)
    np.testing.assert_allclose(f, f.conj().swapaxes(1, 2), atol=0)
    b = f.reshape(n*n-1, n*n).T
    np.testing.assert_allclose(b.conj().T@b, np.eye(n*n-1), atol=4e-16)
    identity = np.eye(n).reshape(-1)/np.sqrt(n)
    np.testing.assert_allclose(b@b.conj().T+np.outer(identity, identity),
                               np.eye(n*n), atol=4e-16)


def test_jump_diffusion_matches_independent_adjoint_product_rule_in_a_nonstationary_state():
    rng = np.random.default_rng(408)
    n = 3
    raw = rng.normal(size=(n, n))+1j*rng.normal(size=(n, n))
    rho = raw@raw.conj().T; rho /= np.trace(rho)
    h = (raw+raw.conj().T)/2
    jumps = tuple(CollapseChannel(
        f"j{k}", rng.normal(size=(n, n))+1j*rng.normal(size=(n, n)),
        "test", "independent random Markov fixture") for k in range(3))
    reservoirs = ExplicitReservoirs(n, jumps)
    f = traceless_hermitian_basis(n)
    by_channel = ordered_jump_diffusion(reservoirs, rho, f)
    reference = adjoint_einstein_diffusion(reservoirs.generator(h), rho, f)
    np.testing.assert_allclose(by_channel.sum(axis=0), reference, rtol=2e-13, atol=2e-14)
    assert np.linalg.eigvalsh(by_channel).min() > -2e-13
    # Hamiltonian contributions cancel in the Einstein product rule.
    np.testing.assert_allclose(reference, adjoint_einstein_diffusion(
        reservoirs.dissipator(), rho, f), rtol=2e-13, atol=2e-14)
    shifted = f+np.arange(len(f))[:, None, None]*np.eye(n)
    np.testing.assert_allclose(ordered_jump_diffusion(reservoirs, rho, shifted),
                               by_channel, rtol=2e-13, atol=2e-14)


def test_ground_state_ordered_spectrum_has_correct_sign_width_and_vacuum_asymmetry():
    gamma, resonance = 2.4, 3.7
    h, reservoirs = _two_level(gamma, resonance)
    model = stationary_atomic_noise(h, reservoirs)
    frequencies = GeneratorFrequencyAxis(np.linspace(-10, 10, 401), "two-level-lab-frame")
    spectra = model.spectrum(frequencies)
    w, a = frequencies.omega_rad_s, gamma/2
    expected = a/(a*a+(w-resonance)**2)
    expected_sym = (expected+a/(a*a+(w+resonance)**2))/2
    # Basis operator 1 is (|g><e|+|e><g|)/sqrt(2).
    np.testing.assert_allclose(spectra.ordered[:, 1, 1], expected, rtol=2e-14, atol=2e-15)
    np.testing.assert_allclose(spectra.symmetrized[:, 1, 1], expected_sym, rtol=2e-14, atol=2e-15)
    assert w[np.argmax(spectra.ordered[:, 1, 1].real)] == pytest.approx(resonance, abs=0.025)
    # Integral over Omega/(2*pi) recovers C_xx=1/2, including the analytic tail.
    bound = 500*gamma
    wide = GeneratorFrequencyAxis(np.linspace(-bound, bound, 30001), "two-level-lab-frame")
    xx = model.spectrum(wide).ordered[:, 1, 1].real
    area = np.trapezoid(xx, wide.omega_rad_s)/(2*np.pi)
    tail = 0.5-(np.arctan((bound-resonance)/a)-np.arctan((-bound-resonance)/a))/(2*np.pi)
    assert area+tail == pytest.approx(model.ordered_covariance[1, 1].real, abs=1e-9)
    assert model.diagnostics.lyapunov_relative_residual < 1e-13


@pytest.mark.parametrize("drive", [0.0, 2.1, 8.0])
def test_qrt_agrees_with_diffusion_including_dc_and_complex_cross_spectra(drive):
    h, reservoirs = _two_level(drive=drive)
    model = stationary_atomic_noise(h, reservoirs)
    axis = GeneratorFrequencyAxis(np.array([-9., -3., 0., 3., 9.]), "two-level-lab-frame")
    result = model.spectrum(axis)
    reference = qrt_ordered_spectrum(model.generator, model.stationary_state,
                                     model.operators, axis.omega_rad_s)
    np.testing.assert_allclose(result.ordered, reference, rtol=1e-12, atol=2e-15)
    np.testing.assert_allclose(result.symmetrized,
                               (reference+reference[::-1].transpose(0, 2, 1))/2,
                               rtol=1e-12, atol=2e-15)
    assert np.linalg.eigvalsh(result.ordered).min() > -2e-14
    assert np.linalg.eigvalsh(result.symmetrized).min() > -2e-14
    assert np.max(np.abs(result.ordered.imag)) > 1e-3
    assert model.diagnostics.commutator_relative_residual < 1e-12


@pytest.mark.parametrize("pump_power,reset_rate", [
    (0., constants.GAMMA_GG), (0.6, 0.), (0.6, constants.GAMMA_GG),
])
def test_rb85_pump_state_parity_and_full_atomic_noise_at_rf_and_beat_frequencies(pump_power, reset_rate):
    pump = fwm.rabi_freq(pump_power, fwm.W_PUMP)
    optical = OpticalDetunings(2*np.pi*0.9e9, -2*np.pi*8e6)
    model = reduced_pump_noise(pump, optical.one_photon_rad_s,
                               transit_rate_s_inverse=reset_rate)
    # The existing pump reference receives the SAME declared physics, not its
    # default collisional model. The independently assembled factory keeps the
    # original radiative decay list and adds only reset collapse operators.
    reference_atom = replace(atoms.double_lambda_rb85(gamma_gg=0), collapse_ops=tuple(
        c.operator for c in thermal_reset_channels(
            reset_rate, [5/12, 7/12, 0, 0], source="parity fixture")))
    reference = fwm.pump_only_weak_response_reference(
        pump, pump, [optical.two_photon_rad_s], [optical.one_photon_rad_s], atom=reference_atom)
    np.testing.assert_allclose(model.stationary_state, reference.pump_state[0],
                               rtol=1e-10, atol=3e-13)
    rf = AnalysisFrequencyAxis.from_hz([-4e6, 0, 4e6])
    shifted = minus_branch_atomic_frequencies(optical, rf)
    for axis in (shifted, GeneratorFrequencyAxis(rf.omega_rad_s, "static-pump-near-zero")):
        spectrum = model.spectrum(axis)
        qrt = qrt_ordered_spectrum(model.generator, model.stationary_state,
                                   model.operators, axis.omega_rad_s)
        for calculated, expected in zip(spectrum.ordered, qrt):
            relative = np.linalg.norm(calculated-expected)/np.linalg.norm(expected)
            assert relative < 2e-8
    assert model.diagnostics.lyapunov_relative_residual < 1e-9
    assert model.diagnostics.decay_gap_s_inverse > model.diagnostics.stability_resolution_s_inverse


def test_nondecaying_modes_and_nonstationary_state_fail_instead_of_losing_elastic_noise():
    h, reservoirs = _two_level()
    with pytest.raises(ValueError, match="not stationary"):
        stationary_atomic_noise(h, reservoirs, stationary_state=np.eye(2)/2)
    with pytest.raises(ValueError, match="nondecaying"):
        stationary_atomic_noise(h, ExplicitReservoirs(2, ()))
    with pytest.raises(ValueError, match="nondecaying"):
        reduced_pump_noise(0, 2*np.pi*0.9e9)  # unconnected ground populations
    with pytest.raises(ValueError, match="nondecaying"):
        qrt_ordered_spectrum(np.zeros((4, 4)), np.eye(2)/2,
                            traceless_hermitian_basis(2), [0.])


def test_atomic_covariance_is_not_assigned_a_canonical_bosonic_commutator():
    h, reservoirs = _two_level()
    model = stationary_atomic_noise(h, reservoirs)
    c = model.ordered_covariance
    # Three SU(2) generators, with state-dependent expected commutators.
    assert c.shape == (3, 3)
    for i, fi in enumerate(model.operators):
        for j, fj in enumerate(model.operators):
            expected = np.trace(model.stationary_state@(fi@fj-fj@fi))
            assert c[i, j]-c[j, i] == pytest.approx(expected, abs=1e-14)
    with pytest.raises(TypeError, match="generator-frame"):
        model.spectrum(AnalysisFrequencyAxis.from_hz([0, 1]))
    with pytest.raises(TypeError):
        minus_branch_atomic_frequencies(AnalysisFrequencyAxis.from_hz([0]),
                                       OpticalDetunings(0, 0))
    with pytest.raises(ValueError):
        model.ordered_covariance.setflags(write=True)
