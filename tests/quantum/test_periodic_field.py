"""Periodic atom -> fields: independent references, physical modes and mean tangent."""

from dataclasses import replace
import json

import numpy as np
import pytest

from analysis.grand_challenge.normalization_audit import conditional_inputs
from analysis.grand_challenge.reference.periodic_field import forced_liouville_response, ordered_current_psd
from analysis.grand_challenge.reference.periodic_qrt import time_domain_qrt
from gabes import core, constants as c
from gabes.fwm_quantum.field import reduced_field_pair
from gabes.fwm_quantum.model import reduced_pump_noise
from gabes.fwm_quantum.normalization import optical_carriers, reduced_dipoles
from gabes.fwm_quantum.periodic import reduced_periodic_noise
from gabes.fwm_quantum.periodic_cell import reduced_periodic_ports, reduced_periodic_cell, integrate_reduced_carriers
from gabes.quantum.contracts import AnalysisFrequencyAxis, GeneratorFrequencyAxis
from gabes.quantum.periodic_field import PeriodicFieldPorts, eliminate_periodic_atom, paired_frequency_channel, adaptive_field_propagation
from gabes.quantum.readout import DetectorResponse
from gabes.quantum.sidebands import sideband_channels
from gabes.quantum.traveling import constant_segment


@pytest.fixture(autouse=True)
def single_thread():
    with core.blas_single_thread():
        yield


def detector():
    return DetectorResponse(AnalysisFrequencyAxis.from_hz([.1e6, 1e6]), [.85, .72],
        [[1., .9], [1+.2j, .9-.1j]], 1.1, [0., 0.], 'declared test detector')


def setup_atom(power=8e-6, zero=False):
    inputs = replace(conditional_inputs(), seed_power_W=power)
    omega = optical_carriers(inputs.detunings)[1:]
    beta = np.array([np.sqrt(power/(c.HBAR*omega[0])), .2j*np.sqrt(power/(c.HBAR*omega[1]))])
    atom = reduced_periodic_noise(inputs, beta*0 if zero else beta, mean_order=4)
    return inputs, atom, reduced_periodic_ports(inputs), beta


def test_zero_seed_full_nambu_reproduces_both_sectors_and_four_mode_channel():
    inputs, atom, ports, _ = setup_atom(zero=True)
    rf = AnalysisFrequencyAxis.from_hz([-1e6, 0., 1e6])
    full = eliminate_periodic_atom(atom, ports, GeneratorFrequencyAxis(rf.omega_rad_s, 'Floquet'),
        response_order=3, linear_density_m_inverse=inputs.number_density_m3*inputs.uniform_area_m2)
    dip = reduced_dipoles('uniform-zeeman-rms')
    static = reduced_pump_noise(dip.pump_rabi_rad_s(inputs.pump_power_W, inputs.pump_waist_m), inputs.one_photon_rad_s,
        transit_rate_s_inverse=inputs.transit_rate_s_inverse, transition_scales=dip.transition_scales)
    pair = reduced_field_pair(static, inputs.detunings, rf, number_density_m3=inputs.number_density_m3,
        uniform_area_m2=inputs.uniform_area_m2, optical_omega_rad_s=optical_carriers(inputs.detunings)[1:],
        effective_dipole_C_m=dip.base_dipole_C_m, transition_scales=dip.transition_scales)
    for indices, old in (([0, 3], pair.main), ([2, 1], pair.companion)):
        for attr in ('drift', 'noise_greater', 'noise_lesser'):
            np.testing.assert_allclose(getattr(full, attr)[:, indices][:, :, indices], getattr(old, attr), rtol=2e-10, atol=2e-11)
    np.testing.assert_allclose(full.drift[:, [0, 3]][:, :, [2, 1]], 0, atol=1e-12)
    t = constant_segment(full, inputs.length_m)
    new = paired_frequency_channel(t, ports.mode_labels, 2, 0)
    old = sideband_channels(constant_segment(pair.main, inputs.length_m), constant_segment(pair.companion, inputs.length_m), rf).channels[0]
    perm = [0, 1, 6, 7, 4, 5, 2, 3]
    np.testing.assert_allclose(new.transfer[perm][:, perm], old.transfer, atol=2e-12)
    np.testing.assert_allclose(new.added_covariance[perm][:, perm], old.added_covariance, atol=2e-12)


def test_microscopic_field_drift_matches_full_liouville_and_noise_matches_time_qrt():
    inputs, atom, ports, _ = setup_atom()
    axis = GeneratorFrequencyAxis([0., 2*np.pi*1e6], 'Floquet')
    density = inputs.number_density_m3*inputs.uniform_area_m2
    local = eliminate_periodic_atom(atom, ports, axis, response_order=3, linear_density_m_inverse=density)
    hs, ops, g, signs, _ = ports.nambu_coordinates()
    args = (atom.reservoirs.generator(atom.hamiltonian_zero_rad_s), core.comm_super(atom.hamiltonian_plus_rad_s),
            core.comm_super(atom.hamiltonian_plus_rad_s.conj().T), atom.beat_rad_s)
    reference = forced_liouville_response(*args, atom.state_harmonics, ops, hs, g, signs, axis.omega_rad_s,
        response_order=3, linear_density=density)
    assert np.linalg.norm(reference-local.drift)/np.linalg.norm(reference) < 1e-10
    qrt = time_domain_qrt(*args, atom.operators, axis.omega_rad_s, [-1, 1], phase_samples=16, rtol=2e-12, atol=2e-14)
    readout = np.zeros((4, 30), complex)
    for j, (h, op, coupling, sign) in enumerate(zip(hs, ops, g, signs)):
        row = 0 if h == -1 else 15
        readout[j, row:row+15] = -1j*sign*coupling*np.einsum('ab,kba->k', op, atom.operators)
    d = density*(readout@qrt['ordered_spectrum']@readout.conj().T)
    assert np.linalg.norm(d-local.noise_greater)/np.linalg.norm(d) < 1e-8


def test_spatial_mean_tangent_equals_dc_quantum_transfer_not_transfer_times_mean():
    inputs = conditional_inputs()
    result = reduced_periodic_cell(inputs, detector(), segments=8)
    full = result['transfer']
    t0 = full.transfer[np.flatnonzero(full.frequency_axis.omega_rad_s == 0)[0]]
    beta = np.array([np.sqrt(inputs.seed_power_W/(c.HBAR*optical_carriers(inputs.detunings)[1])), 0.])
    for j in range(2):
        for phase in (1., 1j):
            direction = np.zeros(2, complex); direction[j] = phase
            eps = np.linalg.norm(beta)*1e-4
            outputs = []
            for sign in (-1, 1):
                solution, scale = integrate_reduced_carriers(inputs, initial_amplitudes=beta+sign*eps*direction)
                outputs.append(solution.y[:, -1]*scale)
            derivative = (outputs[1]-outputs[0])/(2*eps)
            predicted = t0[:2]@np.r_[direction, direction.conj()]
            np.testing.assert_allclose(predicted, derivative, rtol=2e-6, atol=1e-7)
    wrong_mean = t0[:2]@np.r_[beta, beta.conj()]
    assert np.linalg.norm(wrong_mean-result['output_carrier_amplitudes_sqrt_flux'])/np.linalg.norm(beta) > 1e-6


def test_full_quadrature_current_matches_ordered_nambu_with_complex_detector_response():
    det = detector()
    result = reduced_periodic_cell(conditional_inputs(), det, segments=4)
    transfer = result['transfer']
    greater, _ = transfer.vacuum_output()
    for j, w in enumerate(det.analysis_axis.omega_rad_s):
        index = np.flatnonzero(transfer.frequency_axis.omega_rad_s == w)[0]
        reference = ordered_current_psd(greater[index], result['output_carrier_amplitudes_sqrt_flux'],
            det.transmissions, det.current_response[j], det.balance, c.ELEMENTARY_CHARGE)
        assert reference == pytest.approx(result['spectrum'].quantum_psd_A2_Hz[j], rel=1e-12)
    assert min(ch.audit().minimum_cp_eigenvalue for ch in result['channels']) > 0


@pytest.mark.parametrize('zero', ['density', 'length'])
def test_empty_cell_gives_vacuum_sql_and_unchanged_mean(zero):
    inputs = replace(conditional_inputs(), seed_phase_rad=.37,
                     **({'number_density_m3': 0.} if zero == 'density' else {'length_m': 0.}))
    result = reduced_periodic_cell(inputs, detector(), segments=2)
    assert result['probe_power_gain'] == pytest.approx(1., abs=1e-12)
    assert result['conjugate_power_gain'] == pytest.approx(0., abs=1e-20)
    np.testing.assert_allclose(result['spectrum'].quantum_ratio, 1., atol=1e-12)
    np.testing.assert_allclose(np.angle(result['midpoint_carriers_sqrt_flux'][:, 0]), .37, atol=1e-12)


def test_periodic_time_origin_rotates_full_fields_and_additional_ports_stay_explicit():
    inputs, atom, ports, beta = setup_atom()
    phase = .73
    changed = reduced_periodic_noise(inputs, beta*np.exp(1j*np.array([1, -1])*phase), mean_order=4)
    axis = GeneratorFrequencyAxis([2*np.pi*1e6], 'Floquet')
    density = inputs.number_density_m3*inputs.uniform_area_m2
    a = eliminate_periodic_atom(atom, ports, axis, response_order=3, linear_density_m_inverse=density)
    b = eliminate_periodic_atom(changed, ports, axis, response_order=3, linear_density_m_inverse=density)
    u = np.diag(np.exp(1j*ports.nambu_coordinates()[0]*phase))
    np.testing.assert_allclose(b.drift, u@a.drift@u.conj().T, rtol=1e-9, atol=1e-11)
    np.testing.assert_allclose(b.noise_greater, u@a.noise_greater@u.conj().T, rtol=1e-9, atol=1e-11)
    # A declared extra port is a separate model, not a larger atomic cutoff.
    extra = PeriodicFieldPorts(('probe', 'conjugate', 'test-extra-band'), [1, -1, 2],
        np.concatenate([ports.lowering_operators, ports.lowering_operators[:1]]),
        np.r_[ports.coupling_s_inverse_sqrt_flux, ports.coupling_s_inverse_sqrt_flux[0]])
    expanded = eliminate_periodic_atom(atom, extra, axis, response_order=4, linear_density_m_inverse=density)
    assert len(expanded.signs) == 6 and expanded.audit()['passed']


def test_insufficient_atomic_cutoff_and_invalid_physical_coordinates_fail():
    inputs, atom, ports, _ = setup_atom()
    axis = GeneratorFrequencyAxis([2*np.pi*1e6], 'Floquet')
    with pytest.raises(ValueError, match='commutator/PSD failed'):
        eliminate_periodic_atom(atom, ports, axis, response_order=1, linear_density_m_inverse=1e11)
    with pytest.raises(ValueError, match='carrier harmonic'):
        eliminate_periodic_atom(atom, ports, axis, response_order=0, linear_density_m_inverse=1e11)
    with pytest.raises(ValueError):
        replace(ports, carrier_harmonics=[1.5, -1])
    with pytest.raises(ValueError):
        replace(ports, mode_labels=('probe', 'probe'))
    with pytest.raises(ValueError, match='phase mismatch'):
        reduced_periodic_ports(replace(inputs, phase_mismatch_rad_m=1.))
    with pytest.raises(ValueError):
        ports.lowering_operators.setflags(write=True)


def test_adaptive_covariance_propagation_matches_exact_constant_generator_integral():
    inputs, atom, ports, _ = setup_atom()
    local = eliminate_periodic_atom(atom, ports, GeneratorFrequencyAxis([-2*np.pi*1e6, 2*np.pi*1e6], 'Floquet'),
        response_order=3, linear_density_m_inverse=inputs.number_density_m3*inputs.uniform_area_m2)
    reference = constant_segment(local, inputs.length_m)
    actual, diagnostics = adaptive_field_propagation(lambda z: local, inputs.length_m, rtol=2e-11, atol=2e-13)
    for attribute in ('transfer', 'noise_greater', 'noise_lesser'):
        np.testing.assert_allclose(getattr(actual, attribute), getattr(reference, attribute), rtol=2e-10, atol=2e-12)
    assert diagnostics['maximum_local_commutator_relative_residual'] < 1e-10


def test_periodic_cell_report_checks_convergence_refs_and_preserves_existing_output(tmp_path):
    from analysis.grand_challenge.periodic_field_audit import main
    path = tmp_path/'field.json'
    assert main(['--output', str(path)]) == 0
    original = path.read_bytes()
    report = json.loads(original)
    assert report['expected_controls_passed']
    assert report['finite_seed_two_band_field_implemented']
    assert report['fixed_classical_pump']
    assert not report['absolute_hot_vapor_prediction']
    assert not report['experimental_validation']
    assert not report['additional_physical_optical_ports_propagated']
    assert not report['initial_midpoint_convergence_passed']
    assert report['adaptive_tolerance_change_db_8uW'] < 1e-7
    assert report['adaptive_tolerance_change_db_1mW'] < 1e-7
    assert report['independent_full_liouville_drift_relative_error'] < 1e-9
    assert report['independent_time_qrt_field_noise_relative_error'] < 1e-8
    assert report['nonlinear_mean_tangent_maximum_relative_error'] < 1e-6
    assert report['incorrect_transfer_times_mean_relative_error'] > 1e-6
    assert report['configurations'][1]['maximum_cross_sector_pinching_change_db'] > .01
    with pytest.raises(FileExistsError):
        main(['--output', str(path)])
    assert path.read_bytes() == original
