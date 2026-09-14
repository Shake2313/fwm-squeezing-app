"""Two-phase convection, microscopic diffusion and independent trajectory QRT."""

from dataclasses import replace

import numpy as np
import pytest

from analysis.grand_challenge.normalization_audit import conditional_inputs
from analysis.grand_challenge.reference.atomic_qrt import adjoint_einstein_diffusion
from analysis.grand_challenge.reference.torus_qrt import trajectory_qrt, static_grating_qrt
from gabes import core
from gabes.fwm_quantum.kinetic import CarrierGeometry, periodic_velocity_atom, _pump_system
from gabes.fwm_quantum.spatial import spatial_velocity_atom
from gabes.quantum.contracts import AnalysisFrequencyAxis, GeneratorFrequencyAxis
from gabes.quantum.diffusion import stationary_atomic_noise
from gabes.quantum.periodic import periodic_atomic_noise
from gabes.quantum.reservoirs import ExplicitReservoirs, thermal_reset_channels
from gabes.quantum.spatial import torus_atomic_noise


@pytest.fixture(autouse=True)
def single_thread():
    with core.blas_single_thread():
        yield


def fixture():
    h = np.array([[0., .2], [.2, .7]])
    drives = np.array([[[.08, .09j], [.02, 0]], [[.03, .06], [.04j, -.02]]])
    reservoirs = ExplicitReservoirs(2, thermal_reset_channels(1., [.75, .25], source='synthetic contraction fixture'))
    return h, [(1, 0), (-1, 1)], drives, reservoirs, np.array([2.3, np.sqrt(.5)])


def make_atom(**kwargs):
    return torus_atomic_noise(*fixture(), mean_orders=(5, 5), **kwargs)


def test_jump_diffusion_independently_matches_full_adjoint_product_rule():
    atom = make_atom()
    h, labels, drives, reservoirs, _ = fixture()
    phases = np.array([[.17, 1.2], [.81, -.5], [2.2, 4.1]])
    samples = atom.at_phase(phases)
    for i, theta in enumerate(phases):
        v = np.einsum('q,qab->ab', np.exp(-1j*np.asarray(labels)@theta), drives)
        reference = adjoint_einstein_diffusion(reservoirs.generator(h+v+v.conj().T), samples['state'][i], atom.operators)
        np.testing.assert_allclose(samples['diffusion_by_reservoir'][:, i].sum(axis=0), reference, atol=2e-15)
    assert atom.diagnostics['passed']


def test_two_irrational_frequencies_match_independent_history_and_qrt():
    atom = make_atom()
    h, labels, drives, reservoirs, frequencies = fixture()
    # Assemble full density generators independently, never use atomic A or D.
    reference_labels = [(0, 0), *labels, *[tuple(-np.array(q)) for q in labels]]
    generators = [reservoirs.generator(h), *[core.comm_super(v) for v in drives],
                  *[core.comm_super(v.conj().T) for v in drives]]
    outputs = [(0, 0), (1, 0), (-1, 1)]
    reference = trajectory_qrt(generators, reference_labels, frequencies, atom.operators, [.23], outputs,
        phase_samples=(8, 8), history=28., delay=28.)
    actual = atom.spectrum(GeneratorFrequencyAxis([.23], 'two-phase reference'), response_orders=(4, 4), output_labels=outputs)
    np.testing.assert_allclose(atom.at_phase(reference['phases_rad'])['state'], reference['phase_states'], atol=5e-10)
    np.testing.assert_allclose(actual['greater'], reference['greater'], rtol=2e-8, atol=1e-10)
    assert reference['maximum_regression_endpoint_norm'] < 1e-10


def test_zero_modulation_matches_shifted_stationary_spectra_without_folding():
    h, labels, drives, reservoirs, frequencies = fixture()
    atom = torus_atomic_noise(h, labels, drives*0, reservoirs, frequencies, mean_orders=(2, 2))
    outputs = [(0, 0), (1, 0), (-1, 1)]
    axis = GeneratorFrequencyAxis([.13, 7.2], 'unfolded offsets')
    actual = atom.spectrum(axis, response_orders=(2, 2), output_labels=outputs)
    static = stationary_atomic_noise(h, reservoirs)
    for i, q in enumerate(outputs):
        ref = static.spectrum(GeneratorFrequencyAxis(axis.omega_rad_s+np.array(q)@frequencies, 'shifted atomic frequencies'))
        np.testing.assert_allclose(actual['greater'][:, i*3:(i+1)*3, i*3:(i+1)*3], ref.ordered, atol=2e-15)
        for j in range(len(outputs)):
            if i != j:
                np.testing.assert_allclose(actual['greater'][:, i*3:(i+1)*3, j*3:(j+1)*3], 0, atol=1e-15)


def test_lesser_reflects_both_lattice_indices_and_frequency():
    atom = make_atom()
    labels = [(0, 0), (1, 0), (-1, 1), (-1, 0), (1, -1)]
    actual = atom.spectrum(GeneratorFrequencyAxis([-.23, .23], 'reflection'), response_orders=(3, 3), output_labels=labels)
    perm = np.concatenate([np.arange(labels.index(tuple(-np.array(q)))*3, labels.index(tuple(-np.array(q)))*3+3) for q in labels])
    np.testing.assert_allclose(actual['lesser'][1], actual['greater'][0].T[perm][:, perm], atol=2e-14)


def test_independent_phase_origins_rotate_mean_and_cross_spectra():
    h, labels, drives, reservoirs, frequencies = fixture()
    shift = np.array([.38, -.72])
    a = make_atom()
    b = torus_atomic_noise(h, labels, drives*np.exp(1j*np.array(labels)@shift)[:, None, None], reservoirs,
        frequencies, mean_orders=(5, 5))
    np.testing.assert_allclose(b.state_coefficients, a.state_coefficients*np.exp(1j*a.state_labels@shift)[:, None, None], atol=1e-15)
    outputs = [(0, 0), (1, 0), (-1, 1)]
    axis = GeneratorFrequencyAxis([.23], 'phase covariance')
    left, right = [x.spectrum(axis, response_orders=(3, 3), output_labels=outputs) for x in (a, b)]
    u = np.diag(np.repeat(np.exp(1j*np.array(outputs)@shift), 3))
    np.testing.assert_allclose(right['greater'], u@left['greater']@u.conj().T, atol=2e-14)


def test_static_grating_requires_phase_selection_not_torus_average():
    h, labels, drives, reservoirs, _ = fixture()
    atom = torus_atomic_noise(h, labels, drives, reservoirs, [2.3, 0.], mean_orders=(6, 6))
    periodic = periodic_atomic_noise(h, drives[0]+drives[1].conj().T, reservoirs, 2.3, mean_order=9)
    theta = np.linspace(0, 2*np.pi, 13)
    actual = atom.at_phase(np.column_stack([theta, np.zeros(len(theta))]))
    expected = periodic.at_phase(theta)
    np.testing.assert_allclose(actual['state'], expected['state'], atol=2e-11)
    np.testing.assert_allclose(actual['diffusion_by_reservoir'], expected['diffusion_by_reservoir'], atol=2e-11)
    zero_grating = atom.state_coefficients.copy()
    zero_grating[atom.state_labels[:, 1] != 0] = 0
    averaged = replace(atom, state_coefficients=zero_grating).at_phase(np.column_stack([theta, theta*0]))['state']
    assert np.max(abs(averaged-expected['state'])) > 1e-3


def test_insufficient_mean_cutoff_is_rejected_instead_of_psd_repaired():
    with pytest.raises(ValueError, match='consistency failed'):
        torus_atomic_noise(*fixture(), mean_orders=(1, 1))


def test_bad_coordinates_and_noninteger_cutoffs_are_rejected():
    atom = make_atom()
    with pytest.raises(ValueError, match='integer'):
        torus_atomic_noise(*fixture(), mean_orders=(2.5, 2))
    with pytest.raises(TypeError, match='generator-frequency'):
        atom.spectrum(AnalysisFrequencyAxis([.1]), response_orders=(2, 2), output_labels=[(0, 0)])
    for labels in ([(0, 0), (0, 0)], [(3, 0)], [(0.5, 0)]):
        with pytest.raises(ValueError):
            atom.spectrum(GeneratorFrequencyAxis([.1], 'bad labels'), response_orders=(2, 2), output_labels=labels)


def test_nonclosed_rb_geometry_keeps_both_moving_frequencies_and_lab_coordinates():
    inputs = conditional_inputs()
    geometry = CarrierGeometry.vacuum_beams(inputs, probe_angle_rad=.005, conjugate_angle_rad=-.005)
    velocity = np.array([150., 0., 380.])
    atom = spatial_velocity_atom(inputs, geometry, velocity, [5e6+1e6j, .5e6j])
    _, probe, conjugate = geometry.atomic_frequencies(inputs, velocity[None])
    np.testing.assert_allclose(atom.convective_frequencies_rad_s[0], probe[0], rtol=0, atol=1e-6)
    np.testing.assert_allclose(-atom.convective_frequencies_rad_s[0]+atom.convective_frequencies_rad_s[1], conjugate[0], rtol=0, atol=4e-6)
    np.testing.assert_allclose(atom.convective_frequencies_rad_s[1], geometry.mismatch_rad_m@velocity, atol=1e-6)
    assert atom.coordinate_kind == 'two-phase-convective-lattice'
    t = np.array([0., 1e-9]); r0 = np.array([.001, 0., .003])
    expected = atom.atom.at_phase(t[:, None]*atom.convective_frequencies_rad_s-r0@atom.wavevectors_rad_m.T)
    np.testing.assert_allclose(atom.at_events(t, r0+t[:, None]*velocity)['state'], expected['state'], atol=1e-14)
    # v.Q=0 alone does not make the nonzero spatial wavevector disappear.
    static = spatial_velocity_atom(inputs, geometry, [0, 0, 0], [5e6, .5e6j])
    assert static.coordinate_kind == 'two-phase-convective-lattice'
    assert static.convective_frequencies_rad_s[1] == 0


def test_closed_rb_quotient_reuses_the_existing_solver_without_duplicate_modes():
    inputs = conditional_inputs(); geometry = CarrierGeometry.vacuum_beams(inputs)
    velocity = [100., 0., 200.]; beta = [5e6, .5e6j]
    a = spatial_velocity_atom(inputs, geometry, velocity, beta)
    b = periodic_velocity_atom(inputs, geometry, velocity, beta, mean_order=3)
    assert a.coordinate_kind == 'closed-single-phase-quotient'
    np.testing.assert_array_equal(a.atom.state_harmonics, b.state_harmonics)
    np.testing.assert_array_equal(a.atom.diffusion_harmonics_by_reservoir, b.diffusion_harmonics_by_reservoir)
    with pytest.raises(ValueError, match='scalar mismatch'):
        spatial_velocity_atom(replace(inputs, phase_mismatch_rad_m=1.), geometry, velocity, beta)


def test_rb_spatial_grating_matches_independent_period_qrt_projection():
    inputs = conditional_inputs()
    geometry = CarrierGeometry.vacuum_beams(inputs, probe_angle_rad=.005, conjugate_angle_rad=-.005)
    model = spatial_velocity_atom(inputs, geometry, [0, 0, 0], [5e6, .5e6j])
    h0, reservoirs = _pump_system(inputs, inputs.one_photon_rad_s, 'uniform-zeeman-rms')
    pairs = [(core.comm_super(v), core.comm_super(v.conj().T)) for v in model.drive_operators_rad_s]
    labels = [(0, 0), (1, 0), (-1, 1)]
    axis = GeneratorFrequencyAxis([2*np.pi*1e6], 'Rb independent grating reference')
    ref = static_grating_qrt(reservoirs.generator(h0), *pairs, model.convective_frequencies_rad_s[0],
        model.atom.operators, axis.omega_rad_s, labels, grating_samples=4, phase_samples=16)
    actual = model.atom.spectrum(axis, response_orders=(2, 2), output_labels=labels)
    assert np.linalg.norm(actual['greater']-ref['greater'])/np.linalg.norm(ref['greater']) < 1e-9
    np.testing.assert_allclose(model.atom.at_phase(ref['phases_rad'])['state'], ref['phase_states'], atol=1e-10)


def test_nondecaying_response_lattice_is_rejected():
    atom = make_atom()
    with pytest.raises(ValueError, match='nondecaying'):
        replace(atom, drift_coefficients=atom.drift_coefficients*0).spectrum(
            GeneratorFrequencyAxis([.23], 'unstable control'), response_orders=(2, 2), output_labels=[(0, 0)])


def test_audit_will_not_overwrite_an_existing_scientific_artifact(tmp_path, monkeypatch):
    from analysis.grand_challenge import spatial_audit
    output = tmp_path/'existing.json'
    output.write_text('preserved evidence', encoding='utf-8')
    def forbidden():
        raise AssertionError('expensive audit must not start before artifact checks')
    monkeypatch.setattr(spatial_audit, 'build_report', forbidden)
    with pytest.raises(FileExistsError):
        spatial_audit.main(['--output', str(output)])
    assert output.read_text(encoding='utf-8') == 'preserved evidence'
