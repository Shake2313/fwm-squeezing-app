"""Carrier momentum, laboratory/atomic frequencies and independent velocity noise."""

from dataclasses import replace
from math import erf, exp, pi, sqrt

import numpy as np
import pytest

from analysis.grand_challenge.normalization_audit import conditional_inputs
from analysis.grand_challenge.reference.kinetic import pump_state_field_reference
from gabes import core, constants as c, doppler
from gabes.fwm_quantum.inputs import power_normalized_readout
from gabes.fwm_quantum.kinetic import CarrierGeometry, VelocityQuadrature, kinetic_local_field, kinetic_ports, periodic_velocity_atom, thermal_pump_cell
from gabes.fwm_quantum.model import reduced_pump_system
from gabes.fwm_quantum.normalization import optical_carriers, reduced_dipoles
from gabes.fwm_quantum.periodic import reduced_periodic_noise
from gabes.quantum.contracts import AnalysisFrequencyAxis, GeneratorFrequencyAxis
from gabes.quantum.periodic_field import eliminate_periodic_atom
from gabes.quantum.readout import DetectorResponse


@pytest.fixture(autouse=True)
def single_thread():
    with core.blas_single_thread():
        yield


def closed_geometry(inputs, theta=.005):
    vacuum = CarrierGeometry.vacuum_beams(inputs, probe_angle_rad=theta, conjugate_angle_rad=-theta)
    k0, kp, _ = vacuum.wavevectors_rad_m
    return CarrierGeometry([k0, kp, 2*k0-kp], 'synthetic explicitly non-vacuum conjugate; test only')


def det():
    return DetectorResponse(AnalysisFrequencyAxis.from_hz([.1e6, 1e6, 4e6]), [.85, .85], np.ones((3, 2)), 1., np.zeros(3), 'assumed test')


def test_lab_carriers_stay_fixed_and_doppler_follows_direct_moving_plane_wave_phases():
    inputs = conditional_inputs()
    geometry = closed_geometry(inputs)
    velocities = np.array([[50., 0., 130.], [-130., 0., -220.]])
    delta, probe, conjugate = geometry.atomic_frequencies(inputs, velocities)
    omega = optical_carriers(inputs.detunings)
    moving_frequencies = omega[None, :]-velocities@geometry.wavevectors_rad_m.T
    np.testing.assert_allclose(probe, moving_frequencies[:, 1]-moving_frequencies[:, 0], atol=.5, rtol=0)
    np.testing.assert_allclose(conjugate, moving_frequencies[:, 2]-moving_frequencies[:, 0], atol=.5, rtol=0)
    np.testing.assert_allclose(probe+conjugate, 0., atol=1e-5)
    expected_delta, expected_two = doppler.noncollinear_atomic_detunings_rad_s(
        inputs.one_photon_rad_s, inputs.two_photon_rad_s, velocities[:, 0], velocities[:, 2],
        np.linalg.norm(geometry.wavevectors_rad_m[0]), np.linalg.norm(geometry.wavevectors_rad_m[1]), .005)
    np.testing.assert_allclose(delta, expected_delta, atol=1e-6)
    np.testing.assert_allclose(probe, -c.OMEGA_HF+expected_two, atol=1e-5)
    q = VelocityQuadrature(velocities, [.3, .7], 'two velocities')
    r = kinetic_local_field(inputs, geometry, q, AnalysisFrequencyAxis.from_hz([-1e6, 0., 1e6]))
    np.testing.assert_array_equal(r['lab_optical_omega_rad_s'], omega)
    np.testing.assert_array_equal(r['generator'].frequency_axis.omega_rad_s, 2*np.pi*np.array([-1e6, 0., 1e6]))


def test_vacuum_equal_opposite_angles_fail_single_phase_closure_with_measured_defect():
    inputs = conditional_inputs()
    geometry = CarrierGeometry.vacuum_beams(inputs, probe_angle_rad=.005, conjugate_angle_rad=-.005)
    audit = geometry.closure_audit(temperature_K=394.15, length_m=inputs.length_m)
    assert not audit['passed'] and abs(audit['axial_mismatch_phase_rad']) > 2.
    velocities = np.array([[100., 0., 150.]])
    _, p, cc = geometry.atomic_frequencies(inputs, velocities)
    np.testing.assert_allclose(p+cc, velocities@geometry.mismatch_rad_m, atol=1e-5)
    with pytest.raises(ValueError, match='momentum closure'):
        kinetic_ports(inputs, geometry)
    assert CarrierGeometry.vacuum_beams(inputs).closure_audit()['passed']


def test_maxwell_grid_moments_match_conditional_truncated_gaussian_and_tail():
    q = VelocityQuadrature.maxwell_xz(394.15, longitudinal_order=64, transverse_order=48, cutoff_sigma=5.)
    variance = c.KB*394.15/c.MASS_85RB
    retained = erf(5/sqrt(2))
    expected_variance = variance*(1-2*5*exp(-25/2)/sqrt(2*pi)/retained)
    np.testing.assert_allclose(q.probabilities@q.velocities_m_s, 0, atol=1e-12)
    actual = q.velocities_m_s.T@(q.probabilities[:, None]*q.velocities_m_s)
    np.testing.assert_allclose(actual.diagonal()[[0, 2]], expected_variance, rtol=1e-12)
    assert q.omitted_probability == pytest.approx(1-retained**2)
    assert q.integrated_velocity_axes == (1,)


def test_zero_velocity_full_periodic_and_stationary_thermal_limits_match_previous_core():
    inputs = conditional_inputs()
    geometry = CarrierGeometry.vacuum_beams(inputs)
    q = VelocityQuadrature([[0., 0., 0.]], [1.], 'cold limit')
    axis = AnalysisFrequencyAxis.from_hz([-1e6, 0., 1e6])
    beta = [np.sqrt(inputs.seed_power_W/(c.HBAR*optical_carriers(inputs.detunings)[1])), .2j*1e6]
    atom = reduced_periodic_noise(inputs, beta, mean_order=4)
    expected = eliminate_periodic_atom(atom, kinetic_ports(inputs, geometry),
        GeneratorFrequencyAxis(axis.omega_rad_s, 'reference'), response_order=3,
        linear_density_m_inverse=inputs.number_density_m3*inputs.uniform_area_m2)
    actual = kinetic_local_field(inputs, geometry, q, axis, carrier_amplitudes=beta)['generator']
    np.testing.assert_allclose(actual.drift, expected.drift, rtol=1e-10, atol=1e-12)
    np.testing.assert_allclose(actual.noise_greater, expected.noise_greater, rtol=1e-10, atol=1e-12)
    rf = AnalysisFrequencyAxis.from_hz([-4e6, -1e6, -.1e6, 0., .1e6, 1e6, 4e6])
    old = power_normalized_readout(inputs, rf, det())
    new = thermal_pump_cell(inputs, geometry, q, det())
    assert new['probe_power_gain'] == pytest.approx(old.probe_power_gain, rel=1e-11)
    np.testing.assert_allclose(new['spectrum'].quantum_db, old.spectrum.quantum_db, atol=1e-10)


def test_fast_stationary_path_is_exact_zero_seed_periodic_limit_at_nonzero_velocity():
    inputs = conditional_inputs(); geometry = closed_geometry(inputs)
    q = VelocityQuadrature([[90., 0., 230.], [-160., 0., -120.]], [.4, .6], 'two-class test')
    axis = AnalysisFrequencyAxis.from_hz([-1e6, 0., 1e6])
    a = kinetic_local_field(inputs, geometry, q, axis)['generator']
    b = kinetic_local_field(inputs, geometry, q, axis, carrier_amplitudes=[0., 0.])['generator']
    for name in ('drift', 'noise_greater', 'noise_lesser'):
        np.testing.assert_allclose(getattr(a, name), getattr(b, name), rtol=1e-10, atol=2e-11)


def test_velocity_split_and_streaming_sum_preserve_covariances_not_noise_amplitudes():
    inputs = conditional_inputs(); geometry = closed_geometry(inputs)
    axis = AnalysisFrequencyAxis.from_hz([-1e6, 0., 1e6])
    one = VelocityQuadrature([[90., 0., 230.]], [1.], 'one class')
    split = VelocityQuadrature([[90., 0., 230.], [90., 0., 230.]], [.25, .75], 'same class split')
    a = kinetic_local_field(inputs, geometry, one, axis)['generator']
    b = kinetic_local_field(inputs, geometry, split, axis)['generator']
    streaming = kinetic_local_field(inputs, geometry, split, axis, retain_classes=False)['generator']
    for name in ('drift', 'noise_greater', 'noise_lesser'):
        np.testing.assert_allclose(getattr(a, name), getattr(b, name), rtol=1e-11, atol=1e-12)
        np.testing.assert_allclose(getattr(b, name), getattr(streaming, name), rtol=1e-11, atol=1e-12)
    wrong_factor = (sqrt(.25)+sqrt(.75))**2
    wrong = replace(a, noise_greater_by_reservoir=a.noise_greater_by_reservoir*wrong_factor,
                    noise_lesser_by_reservoir=a.noise_lesser_by_reservoir*wrong_factor)
    assert not wrong.audit()['passed']


def test_finite_seed_weighted_mean_derivative_matches_summed_dc_drift():
    inputs = conditional_inputs(); geometry = closed_geometry(inputs)
    q = VelocityQuadrature([[90., 0., 230.], [-160., 0., -120.]], [.4, .6], 'two classes')
    beta = np.array([3e6+1e6j, .2e6j])
    direction = np.array([1+.2j, -.3+.7j])
    axis = AnalysisFrequencyAxis.from_hz([0.])
    reference = kinetic_local_field(inputs, geometry, q, axis, carrier_amplitudes=beta)
    derivative = []
    for sign in (-1, 1):
        result = kinetic_local_field(inputs, geometry, q, axis, carrier_amplitudes=beta+sign*100*direction)
        derivative.append(-1j*inputs.number_density_m3*inputs.uniform_area_m2*
            result['ports'].coupling_s_inverse_sqrt_flux*result['mean_polarization'])
    finite_difference = (derivative[1]-derivative[0])/200
    actual = (reference['generator'].drift[0]@np.r_[direction, direction.conj()])[:2]
    np.testing.assert_allclose(actual, finite_difference, rtol=1e-6, atol=1e-6)


def test_transverse_reflection_symmetry_and_zero_density_limit():
    inputs = conditional_inputs()
    q = VelocityQuadrature.maxwell_xz(394.15, longitudinal_order=8, transverse_order=8)
    a = thermal_pump_cell(inputs, closed_geometry(inputs), q, det())
    b = thermal_pump_cell(inputs, closed_geometry(inputs, -.005), q, det())
    np.testing.assert_allclose(a['spectrum'].quantum_db, b['spectrum'].quantum_db, atol=1e-11)
    empty = thermal_pump_cell(replace(inputs, number_density_m3=0), closed_geometry(inputs), q, det())
    assert empty['probe_power_gain'] == pytest.approx(1., abs=1e-12)
    np.testing.assert_allclose(empty['spectrum'].quantum_ratio, 1., atol=1e-12)


def test_z_plane_flux_and_invalid_quadrature_axes_are_explicit():
    inputs = conditional_inputs(); geometry = closed_geometry(inputs)
    collinear = kinetic_ports(inputs, CarrierGeometry.vacuum_beams(inputs))
    ports = kinetic_ports(inputs, geometry)
    cosine = geometry.wavevectors_rad_m[1:, 2]/np.linalg.norm(geometry.wavevectors_rad_m[1:], axis=1)
    np.testing.assert_allclose(ports.coupling_s_inverse_sqrt_flux, collinear.coupling_s_inverse_sqrt_flux/np.sqrt(cosine), rtol=1e-14)
    bad = geometry.wavevectors_rad_m.copy();bad[1,1]=1.;bad[2,1]=-1.
    with pytest.raises(ValueError, match='integrated out'):
        kinetic_local_field(inputs, CarrierGeometry(bad, 'noncoplanar'),
            VelocityQuadrature.maxwell_xz(394.15,longitudinal_order=8,transverse_order=8), AnalysisFrequencyAxis.from_hz([0.]))
    with pytest.raises(ValueError):
        VelocityQuadrature([[0.,0.,0.]], [.5], 'invalid normalization')
    with pytest.raises(ValueError):
        VelocityQuadrature.maxwell_xz(394.15,longitudinal_order=2.5,transverse_order=8)
    with pytest.raises(ValueError):
        geometry.wavevectors_rad_m.setflags(write=True)


def test_velocity_averaged_M_and_both_D_match_independent_full_density_response_and_qrt():
    inputs = conditional_inputs(); geometry = closed_geometry(inputs)
    q = VelocityQuadrature([[150., 0., 380.], [-110., 0., -220.], [0., 0., 700.]], [.4, .5, .1], 'off-resonant and resonant velocity fixture')
    axis = AnalysisFrequencyAxis.from_hz([-1e6, 0., 1e6])
    actual = kinetic_local_field(inputs, geometry, q, axis, retain_classes=False)
    hs, ops, gs, signs, _ = actual['ports'].nambu_coordinates()
    dipoles = reduced_dipoles('uniform-zeeman-rms')
    totals = {name: np.zeros_like(actual['generator'].drift) for name in ('drift', 'greater', 'lesser')}
    for weight, delta, beats in zip(q.probabilities, actual['atomic_one_photon_rad_s'], actual['atomic_signed_beats_rad_s']):
        h, reservoirs = reduced_pump_system(dipoles.pump_rabi_rad_s(inputs.pump_power_W, inputs.pump_waist_m), delta,
            transit_rate_s_inverse=inputs.transit_rate_s_inverse, transition_scales=dipoles.transition_scales)
        reference = pump_state_field_reference(reservoirs.generator(h), ops, gs, signs, hs, beats[0],
            axis.omega_rad_s, inputs.number_density_m3*inputs.uniform_area_m2*weight)
        for name in totals:
            totals[name] += reference[name]
    for expected, name in ((totals['drift'], 'drift'), (totals['greater'], 'noise_greater'), (totals['lesser'], 'noise_lesser')):
        assert np.linalg.norm(getattr(actual['generator'], name)-expected)/np.linalg.norm(expected) < 1e-9
