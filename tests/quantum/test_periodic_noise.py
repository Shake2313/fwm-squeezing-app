"""Dynamic Einstein identity, cyclic orderings and independent periodic QRT."""

import json

import numpy as np
import pytest

from analysis.grand_challenge.reference.atomic_qrt import adjoint_einstein_diffusion
from analysis.grand_challenge.reference.periodic_qrt import time_domain_qrt
from gabes import core
from gabes.quantum.contracts import GeneratorFrequencyAxis, AnalysisFrequencyAxis
from gabes.quantum.diffusion import ordered_jump_diffusion, stationary_atomic_noise
from gabes.quantum.periodic import periodic_atomic_noise
from gabes.quantum.reservoirs import ExplicitReservoirs, CollapseChannel


def fixture():
    h = np.array([[0., .3], [.3, 1.]])
    v = np.array([[.2, .25j], [.12, 0.]])
    reservoirs = ExplicitReservoirs(2, (CollapseChannel('decay', [[0., 1.2], [0., 0.]], 'radiative', 'fixture'),))
    return h, v, reservoirs, 4.


def test_time_dependent_jump_noise_matches_product_rule_and_dynamic_balance():
    h, v, reservoirs, beat = fixture()
    model = periodic_atomic_noise(h, v, reservoirs, beat, mean_order=7)
    phases = np.array([.17, .83, 2.13, 4.21])
    sampled = model.at_phase(phases)
    for k, phi in enumerate(phases):
        jump = ordered_jump_diffusion(reservoirs, sampled['state'][k], model.operators)
        np.testing.assert_allclose(jump, sampled['diffusion_by_reservoir'][:, k], atol=2e-15)
        H = h+v*np.exp(-1j*phi)+v.conj().T*np.exp(1j*phi)
        einstein = adjoint_einstein_diffusion(reservoirs.generator(H), sampled['state'][k], model.operators)
        np.testing.assert_allclose(einstein, jump.sum(axis=0), atol=3e-15)
    assert model.diagnostics['maximum_dynamic_einstein_relative_residual'] < 1e-10
    assert model.diagnostics['maximum_dynamic_commutator_relative_residual'] < 1e-10
    assert np.linalg.norm(sampled['diffusion_by_reservoir'][:, 0]-sampled['diffusion_by_reservoir'][:, 1]) > .01


@pytest.mark.parametrize('beat', [4., -4.])
def test_zero_modulation_reduces_to_independent_stationary_frequency_blocks(beat):
    h, v, reservoirs, _ = fixture()
    model = periodic_atomic_noise(h, v*0, reservoirs, beat, mean_order=3)
    stationary = stationary_atomic_noise(h, reservoirs)
    axis = GeneratorFrequencyAxis([-.3, .2], 'Floquet test')
    spectrum = model.lift(2).spectrum(axis)
    for i in range(-2, 3):
        ref = stationary.spectrum(GeneratorFrequencyAxis(axis.omega_rad_s+i*beat, 'static test'))
        np.testing.assert_allclose(spectrum.block(i, i), ref.ordered, atol=1e-14)
        for j in range(-2, 3):
            if i != j:
                np.testing.assert_allclose(spectrum.block(i, j), 0, atol=1e-14)


def test_lesser_ordering_requires_atomic_transpose_and_frequency_harmonic_reflection():
    h, v, reservoirs, beat = fixture()
    model = periodic_atomic_noise(h, v, reservoirs, beat, mean_order=7)
    spectrum = model.lift(3).spectrum(GeneratorFrequencyAxis([-.4, .4], 'Floquet'))
    reflect = np.concatenate([np.arange(i*3, (i+1)*3) for i in range(6, -1, -1)])
    expected = spectrum.greater[0].T[reflect][:, reflect]
    np.testing.assert_allclose(spectrum.lesser[1], expected, rtol=2e-12, atol=1e-14)
    assert np.linalg.norm(spectrum.greater-spectrum.lesser) > .1


def test_time_origin_rotates_cross_harmonics_without_changing_diagonal_spectra():
    h, v, reservoirs, beat = fixture()
    phase = .73
    a = periodic_atomic_noise(h, v, reservoirs, beat, mean_order=7)
    b = periodic_atomic_noise(h, v*np.exp(1j*phase), reservoirs, beat, mean_order=7)
    axis = GeneratorFrequencyAxis([.2], 'Floquet')
    left, right = a.lift(3).spectrum(axis), b.lift(3).spectrum(axis)
    u = np.diag(np.repeat(np.exp(1j*np.arange(-3, 4)*phase), 3))
    np.testing.assert_allclose(right.greater, u@left.greater@u.conj().T, atol=1e-13)


def test_negative_beat_with_conjugate_drive_describes_the_same_physical_model():
    h, v, reservoirs, beat = fixture()
    a = periodic_atomic_noise(h, v, reservoirs, beat, mean_order=7)
    b = periodic_atomic_noise(h, v.conj().T, reservoirs, -beat, mean_order=7)
    axis = GeneratorFrequencyAxis([.2], 'Floquet')
    left, right = a.lift(2).spectrum(axis), b.lift(2).spectrum(axis)
    for i in (-1, 0, 1):
        for j in (-1, 0, 1):
            np.testing.assert_allclose(left.block(i, j), right.block(-i, -j), atol=1e-13)


def test_cyclic_spectrum_agrees_with_time_domain_full_liouville_qrt_and_converges():
    h, v, reservoirs, beat = fixture()
    model = periodic_atomic_noise(h, v, reservoirs, beat, mean_order=7)
    axis = GeneratorFrequencyAxis([-.3, 0., .6], 'Floquet')
    with core.blas_single_thread():
        reference = time_domain_qrt(reservoirs.generator(h), core.comm_super(v), core.comm_super(v.conj().T),
            beat, model.operators, axis.omega_rad_s, [-1, 0, 1], phase_samples=24)
    errors = []
    for order in (1, 2, 3, 4):
        spectrum = model.lift(order).spectrum(axis)
        indices = np.arange((order-1)*3, (order+2)*3)
        actual = spectrum.greater[:, indices][:, :, indices]
        errors.append(np.linalg.norm(actual-reference['ordered_spectrum'])/np.linalg.norm(reference['ordered_spectrum']))
    assert errors[-1] < 2e-9
    assert errors[0] > 1e-5 and errors[-1] < errors[0]/1e4
    np.testing.assert_allclose(model.at_phase(reference['physical_phases_rad'])['state'], reference['phase_states'], atol=5e-10)


def test_averaged_noise_is_not_the_microscopic_periodic_noise():
    h, v, reservoirs, beat = fixture()
    model = periodic_atomic_noise(h, v, reservoirs, beat, mean_order=7)
    lift = model.lift(3)
    total = lift.greater_by_reservoir.sum(axis=0)
    average = np.kron(np.eye(7), model.diffusion_harmonics_by_reservoir[:, 7].sum(axis=0))
    assert np.linalg.norm(total-average)/np.linalg.norm(total) > .01
    assert lift.audit()['passed']
    # The nonzero Fourier coefficient is not itself Hermitian/PSD; only the
    # full phase diffusion and harmonic Toeplitz matrix have that requirement.
    d1 = model.diffusion_harmonics_by_reservoir[:, 8].sum(axis=0)
    assert np.linalg.norm(d1-d1.conj().T) > .01


def test_invalid_mean_truncation_unresolved_modes_and_frequency_coordinates_fail():
    h, v, reservoirs, beat = fixture()
    with pytest.raises(ValueError, match='consistency'):
        periodic_atomic_noise(h, v, reservoirs, beat, mean_order=1)
    unitary = ExplicitReservoirs(2, (CollapseChannel('zero', np.zeros((2, 2)), 'test', 'test'),))
    with pytest.raises(ValueError, match='nondecaying'):
        periodic_atomic_noise(h, v, unitary, beat, mean_order=7)
    model = periodic_atomic_noise(h, v, reservoirs, beat, mean_order=7)
    with pytest.raises(TypeError):
        model.lift(2).spectrum(AnalysisFrequencyAxis([0.]))
    with pytest.raises(ValueError, match='base frequencies'):
        model.lift(2).spectrum(GeneratorFrequencyAxis([2.], 'edge'))
    with pytest.raises(ValueError):
        model.lift(-1)
    with pytest.raises(ValueError):
        model.state_harmonics.setflags(write=True)


def test_reduced_rb_periodic_report_independent_controls_and_artifact_preservation(tmp_path):
    from analysis.grand_challenge.periodic_noise_audit import main

    path = tmp_path/'periodic.json'
    assert main(['--output', str(path)]) == 0
    original = path.read_bytes()
    report = json.loads(original)
    assert report['expected_controls_passed']
    assert report['periodic_atomic_quantum_noise_implemented']
    assert not report['finite_seed_field_quantum_noise_implemented']
    assert not report['absolute_hot_vapor_prediction']
    assert not report['experimental_validation']
    assert report['zero_seed_maximum_relative_spectral_error'] < 1e-10
    zero, entrance, exit_case, strong = report['configurations']
    for case in report['configurations']:
        assert case['passed']
        assert case['independent_time_qrt_relative_error'] < 1e-8
        assert case['mean_diagnostics']['maximum_dynamic_einstein_relative_residual'] < 1e-10
        assert case['mean_diagnostics']['maximum_dynamic_commutator_relative_residual'] < 1e-10
    assert zero['cross_harmonic_spectrum_relative_norm'] < 1e-12
    assert exit_case['averaged_noise_spectrum_relative_change'] > 1e-4
    assert strong['averaged_noise_spectrum_relative_change'] > .01
    assert strong['adjacent_response_interior_relative_changes'][0] > 1e-6
    assert strong['adjacent_response_interior_relative_changes'][-1] < 1e-10
    # A small mean back-action does not bound a selected atomic noise component.
    assert exit_case['probe_dipole_lesser_spectrum_s'][5]/zero['probe_dipole_lesser_spectrum_s'][5] > 1.2
    assert entrance['time_dependent_diffusion_harmonics_relative_norm'] > 1e-3
    with pytest.raises(FileExistsError):
        main(['--output', str(path)])
    assert path.read_bytes() == original
