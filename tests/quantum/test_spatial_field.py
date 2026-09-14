"""Normalized physical mode projection, microscopic noise and nonlinear tangent."""

from dataclasses import replace

import numpy as np
import pytest

from analysis.grand_challenge.normalization_audit import conditional_inputs
from analysis.grand_challenge.reference.periodic_field import ordered_current_psd
from analysis.grand_challenge.reference.spatial_field import spatial_field_reference
from gabes import core, constants as c
from gabes.fwm_quantum.kinetic import CarrierGeometry
from gabes.fwm_quantum.normalization import optical_carriers
from gabes.fwm_quantum.periodic import reduced_periodic_noise
from gabes.fwm_quantum.periodic_cell import reduced_periodic_ports, integrate_reduced_carriers
from gabes.fwm_quantum.spatial_cell import StationarySpatialMedium, stationary_spatial_cell
from gabes.fwm_quantum.spatial import spatial_velocity_atom
from gabes.quantum.contracts import AnalysisFrequencyAxis, GeneratorFrequencyAxis
from gabes.quantum.periodic_field import eliminate_periodic_atom
from gabes.quantum.periodic import periodic_atomic_noise
from gabes.quantum.readout import DetectorResponse
from gabes.quantum.spatial_modes import TransverseModeGrid


@pytest.fixture(autouse=True)
def single_thread():
    with core.blas_single_thread():
        yield


def setup(order=2, closed=False):
    inputs = conditional_inputs()
    geometry = (CarrierGeometry.vacuum_beams(inputs) if closed else
        CarrierGeometry.vacuum_beams(inputs, probe_angle_rad=.005, conjugate_angle_rad=-.004))
    modes = TransverseModeGrid.rectangle(.0004, .0003, order_x=order)
    return StationarySpatialMedium(inputs, geometry, modes)


def detector():
    return DetectorResponse(AnalysisFrequencyAxis.from_hz([.1e6, 1e6]), [.85, .72],
        [[1., .9], [1+.2j, .9-.1j]], 1.1, [0., 0.], 'declared conditional detector')


def test_profiles_are_normalized_without_hidden_rescaling():
    modes = TransverseModeGrid.rectangle(.0004, .0003, order_x=3, order_y=2)
    np.testing.assert_allclose(modes.area_weights_m2@abs(modes.mode_values_m_inverse)**2, 1, atol=1e-15)
    with pytest.raises(ValueError, match='no silent normalization'):
        replace(modes, mode_values_m_inverse=modes.mode_values_m_inverse*1.01)
    with pytest.raises(ValueError, match='positive area'):
        replace(modes, area_weights_m2=-modes.area_weights_m2)
    with pytest.raises(ValueError, match='integer'):
        TransverseModeGrid.rectangle(.001, .001, order_x=1.5)


def test_original_optical_phases_and_independent_time_qrt_match_projected_M_D():
    medium = setup()
    u = medium.modes.mode_values_m_inverse.copy()
    u[:, 0] *= np.sqrt([1.3, .7])
    u[:, 1] *= np.sqrt([.8, 1.2])*np.exp(1j*np.array([.3, -.4]))
    medium = StationarySpatialMedium(medium.inputs, medium.geometry, replace(medium.modes, mode_values_m_inverse=u))
    z, beta = .006, np.array([5e6+1e6j, .5e6j])
    axis = GeneratorFrequencyAxis([0., 2*np.pi*1e6], 'independent optical phase reference')
    actual = medium.local_field(z, beta, axis)
    k0, kp, kc = medium.geometry.wavevectors_rad_m
    ref = spatial_field_reference(medium.h0, medium.reservoirs, medium.beat_rad_s,
        medium.ports.lowering_operators, medium.ports.coupling_s_inverse_sqrt_flux, beta,
        medium.modes.positions_m, medium.modes.area_weights_m2, u, medium.inputs.uniform_area_m2,
        np.array([kp-k0, kc-k0]), z, medium.inputs.number_density_m3, axis.omega_rad_s)
    for name, expected in (('drift', ref['drift']), ('noise_greater', ref['greater']), ('noise_lesser', ref['lesser'])):
        assert np.linalg.norm(getattr(actual, name)-expected)/np.linalg.norm(expected) < 1e-8
    np.testing.assert_allclose(medium.mean_rate(z, beta), ref['mean_rate'], rtol=2e-11, atol=1e-4)


def test_reference_area_is_a_normalization_coordinate_not_a_fit_parameter():
    a = setup()
    b = StationarySpatialMedium(replace(a.inputs, uniform_area_m2=a.inputs.uniform_area_m2*7), a.geometry, a.modes)
    beta = [5e6, .5e6j]; axis = GeneratorFrequencyAxis([0., 2*np.pi*1e6], 'area invariance')
    left, right = [m.local_field(.004, beta, axis) for m in (a, b)]
    for name in ('drift', 'noise_greater', 'noise_lesser'):
        np.testing.assert_allclose(getattr(left, name), getattr(right, name), rtol=2e-10, atol=1e-10)
    np.testing.assert_allclose(a.mean_rate(.004, beta), b.mean_rate(.004, beta), rtol=2e-10)


def test_splitting_the_same_area_element_does_not_add_correlated_noise():
    a = setup()
    grid = a.modes
    b = StationarySpatialMedium(a.inputs, a.geometry, TransverseModeGrid(
        np.repeat(grid.positions_m, 2, axis=0), np.repeat(grid.area_weights_m2/2, 2),
        np.repeat(grid.mode_values_m_inverse, 2, axis=0), 'same quadrature cells split in half'))
    assert len(a.sections(.004)) == len(b.sections(.004))
    axis = GeneratorFrequencyAxis([2*np.pi*1e6], 'area splitting')
    left, right = [m.local_field(.004, [5e6, .5e6j], axis) for m in (a, b)]
    for name in ('drift', 'noise_greater', 'noise_lesser'):
        np.testing.assert_array_equal(getattr(left, name), getattr(right, name))


def test_collinear_quotient_matches_existing_mean_and_field():
    medium = setup(order=5, closed=True)
    assert len(medium.sections(.004)) == 1
    beta = [5e6, .5e6j]; axis = GeneratorFrequencyAxis([0., 2*np.pi*1e6], 'closed limit')
    a = medium.local_field(.004, beta, axis)
    atom = reduced_periodic_noise(medium.inputs, beta, mean_order=4)
    b = eliminate_periodic_atom(atom, reduced_periodic_ports(medium.inputs), axis, response_order=3,
        linear_density_m_inverse=medium.inputs.number_density_m3*medium.inputs.uniform_area_m2)
    for name in ('drift', 'noise_greater', 'noise_lesser'):
        np.testing.assert_allclose(getattr(a, name), getattr(b, name), rtol=2e-10, atol=2e-10)
    lhs, scale = medium.integrate_carriers()
    rhs, old_scale = integrate_reduced_carriers(medium.inputs)
    np.testing.assert_allclose(lhs.y[:, -1]*scale, rhs.y[:, -1]*old_scale, rtol=2e-10)


def test_dc_local_generator_is_the_derivative_of_the_projected_nonlinear_mean():
    medium = setup()
    beta = np.array([5e6+1e6j, .5e6j])
    local = medium.local_field(.004, beta, GeneratorFrequencyAxis([0.], 'DC tangent'))
    for direction in (np.array([1., .2j]), np.array([1j, -.3])):
        epsilon = 100.
        actual = (medium.mean_rate(.004, beta+epsilon*direction)-medium.mean_rate(.004, beta-epsilon*direction))/(2*epsilon)
        expected = (local.drift[0]@np.r_[direction, direction.conj()])[:2]
        assert np.linalg.norm(actual-expected)/np.linalg.norm(expected) < 1e-6


@pytest.fixture(scope='module')
def cell():
    with core.blas_single_thread():
        medium = setup()
        result = stationary_spatial_cell(medium, detector())
    return medium, result


def test_nonclosed_cell_preserves_field_commutator_and_physical_readout(cell):
    medium, result = cell
    assert not medium.geometry.closure_audit()['passed']
    assert result['transfer'].audit()['passed']
    assert result['minimum_channel_cp_eigenvalue'] >= -1e-10
    assert result['spectrum'].minimum_covariance_uncertainty_eigenvalue >= -1e-10
    transfer, beta = result['transfer'], result['output_carriers_sqrt_flux']
    vac = np.diag([1., 1., 0., 0.])
    for j, w in enumerate(detector().analysis_axis.omega_rad_s):
        i = np.flatnonzero(transfer.frequency_axis.omega_rad_s == w)[0]
        greater = transfer.transfer[i]@vac@transfer.transfer[i].conj().T+transfer.noise_greater[i]
        ref = ordered_current_psd(greater, beta, detector().transmissions, detector().current_response[j], detector().balance, c.ELEMENTARY_CHARGE)
        assert abs(ref/result['spectrum'].quantum_psd_A2_Hz[j]-1) < 1e-11


def test_dc_cell_transfer_is_nonlinear_mean_map_tangent(cell):
    medium, result = cell
    scale = np.sqrt(medium.inputs.seed_power_W/(c.HBAR*optical_carriers(medium.inputs.detunings)[1]))
    initial = np.array([scale, 0.]); direction = np.array([.4+.2j, -.3+.7j])
    outputs = []
    epsilon = 100.
    for sign in (-1, 1):
        trajectory, norm = medium.integrate_carriers(initial_amplitudes=initial+sign*epsilon*direction)
        outputs.append(trajectory.y[:, -1]*norm)
    transfer = result['transfer']; index = np.flatnonzero(transfer.frequency_axis.omega_rad_s == 0)[0]
    expected = (transfer.transfer[index]@np.r_[direction, direction.conj()])[:2]
    assert np.linalg.norm((outputs[1]-outputs[0])/(2*epsilon)-expected)/np.linalg.norm(expected) < 1e-6


@pytest.mark.parametrize('no_atoms', [True, False])
def test_no_atoms_and_zero_length_return_vacuum_shot_noise(no_atoms):
    medium = setup()
    inputs = replace(medium.inputs, **({'number_density_m3': 0.} if no_atoms else {'length_m': 0.}))
    medium = StationarySpatialMedium(inputs, medium.geometry, medium.modes)
    result = stationary_spatial_cell(medium, detector())
    np.testing.assert_allclose(result['transfer'].transfer, np.broadcast_to(np.eye(4), result['transfer'].transfer.shape), atol=1e-13)
    np.testing.assert_allclose(result['spectrum'].quantum_ratio, 1., atol=1e-13)


def test_transport_and_double_counted_mismatch_are_explicitly_rejected():
    medium = setup()
    for velocity in ([1., 0., 0.], [0., 0., 1.]):
        with pytest.raises(ValueError, match='transport noise'):
            StationarySpatialMedium(medium.inputs, medium.geometry, medium.modes, velocity_m_s=velocity)
    with pytest.raises(ValueError, match='scalar mismatch'):
        StationarySpatialMedium(replace(medium.inputs, phase_mismatch_rad_m=1.), medium.geometry, medium.modes)


def test_Q_to_zero_limit_is_continuous_at_fixed_physical_aperture():
    closed = setup(closed=True)
    axis = GeneratorFrequencyAxis([2*np.pi*1e6], 'fixed-aperture geometry limit')
    beta = [5e6, .5e6j]
    reference = closed.local_field(.006, beta, axis)
    errors = []
    for theta in (.001, .0003, .0001):
        geometry = CarrierGeometry.vacuum_beams(closed.inputs, probe_angle_rad=theta, conjugate_angle_rad=-theta)
        local = StationarySpatialMedium(closed.inputs, geometry, closed.modes).local_field(.006, beta, axis)
        errors.append(np.linalg.norm(local.drift-reference.drift)/np.linalg.norm(reference.drift))
    assert errors[1] < .15*errors[0]
    assert errors[2] < .15*errors[1]


def test_physical_mode_section_is_the_stationary_section_of_the_two_phase_atom():
    medium = setup()
    beta = np.array([5e6+1e6j, .5e6j])
    torus = spatial_velocity_atom(medium.inputs, medium.geometry, [0, 0, 0], beta).atom
    z, theta = .006, .17
    for xy, _area, values in zip(medium.modes.positions_m, medium.modes.area_weights_m2, medium.modes.mode_values_m_inverse):
        loop = -medium.geometry.mismatch_rad_m@np.r_[xy, z]
        eta = np.sqrt(medium.inputs.uniform_area_m2)*values.copy()
        eta[1] *= np.exp(1j*loop)
        periodic = periodic_atomic_noise(medium.h0, medium.harmonic_drive(beta, eta), medium.reservoirs,
            medium.beat_rad_s, mean_order=4)
        a, b = torus.at_phase([[theta, -loop]]), periodic.at_phase([theta])
        np.testing.assert_allclose(a['state'], b['state'], atol=2e-12)
        np.testing.assert_allclose(a['diffusion_by_reservoir'], b['diffusion_by_reservoir'], rtol=1e-8, atol=2e-7)


def test_spatial_field_audit_preserves_existing_artifacts(tmp_path, monkeypatch):
    from analysis.grand_challenge import spatial_field_audit
    target = tmp_path/'kept.json'
    target.write_text('unchanged scientific evidence', encoding='utf-8')
    def forbidden():
        raise AssertionError('must check artifact paths before computation')
    monkeypatch.setattr(spatial_field_audit, 'build_report', forbidden)
    with pytest.raises(FileExistsError):
        spatial_field_audit.main(['--output', str(target)])
    assert target.read_text(encoding='utf-8') == 'unchanged scientific evidence'


def test_positive_noise_does_not_validate_a_frozen_moving_grating():
    from analysis.grand_challenge.spatial_field_audit import frozen_moving_grating_control
    result = frozen_moving_grating_control()
    assert result['passed']
    assert result['candidate_audit']['maximum_commutator_relative_residual'] > 1e-6
