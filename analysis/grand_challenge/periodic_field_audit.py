"""Conditional finite-seed two-band microscopic field propagation and readout.

python -m analysis.grand_challenge.periodic_field_audit --output NEW.json --plot NEW.png
"""

import argparse
from dataclasses import replace
import hashlib
import json
from pathlib import Path
import platform

import numpy as np
import scipy

from gabes import core, constants as c
from gabes.fwm_quantum.inputs import INPUT_UNITS, power_normalized_readout
from gabes.fwm_quantum.normalization import optical_carriers
from gabes.fwm_quantum.periodic import reduced_periodic_noise
from gabes.fwm_quantum.periodic_cell import integrate_reduced_carriers, reduced_periodic_cell, reduced_periodic_ports
from gabes.quantum.contracts import AnalysisFrequencyAxis, GeneratorFrequencyAxis
from gabes.quantum.periodic_field import eliminate_periodic_atom, paired_frequency_channel, adaptive_field_propagation
from gabes.quantum.readout import DetectorResponse, intensity_difference_spectrum
from gabes.quantum.sidebands import SIDEBAND_MODES
from gabes.quantum.traveling import constant_segment, compose_segments
from .normalization_audit import conditional_inputs
from .reference.periodic_field import forced_liouville_response, ordered_current_psd
from .reference.periodic_qrt import time_domain_qrt


ROOT = Path(__file__).resolve().parents[2]


def _complex(a):
    a = np.asarray(a)
    return {'real': a.real.tolist(), 'imag': a.imag.tolist()}


def _relative(a, b):
    return float(np.linalg.norm(a-b)/max(np.linalg.norm(b), np.finfo(float).tiny))


def _pinched_readout(result, inputs, detector):
    """Diagnostic: remove cross-sector couplings and noise, preserving each block."""
    mask = np.equal.outer([0, 1, 1, 0], [0, 1, 1, 0])
    total = None
    for local in result['local_generators']:
        altered = replace(local, drift=local.drift*mask,
            noise_greater_by_reservoir=local.noise_greater_by_reservoir*mask,
            noise_lesser_by_reservoir=local.noise_lesser_by_reservoir*mask)
        segment = constant_segment(altered, inputs.length_m/len(result['local_generators']))
        total = segment if total is None else compose_segments(total, segment)
    if result['quantum_propagation']['method'] == 'DOP853':
        trajectory, scale = integrate_reduced_carriers(inputs)
        ports = reduced_periodic_ports(inputs)
        axis = result['transfer'].frequency_axis
        def local_at_z(z):
            atom = reduced_periodic_noise(inputs, trajectory.sol(z)*scale, mean_order=4)
            local = eliminate_periodic_atom(atom, ports, axis, response_order=3,
                linear_density_m_inverse=inputs.number_density_m3*inputs.uniform_area_m2)
            return replace(local, drift=local.drift*mask,
                noise_greater_by_reservoir=local.noise_greater_by_reservoir*mask,
                noise_lesser_by_reservoir=local.noise_lesser_by_reservoir*mask)
        total, _ = adaptive_field_propagation(local_at_z, inputs.length_m, rtol=2e-11, atol=2e-13)
    covariances = []
    perm = [0, 1, 6, 7, 4, 5, 2, 3]
    for w in detector.analysis_axis.omega_rad_s:
        index, mirror = [int(np.flatnonzero(total.frequency_axis.omega_rad_s == sign*w)[0]) for sign in (1, -1)]
        channel = paired_frequency_channel(total, ('probe', 'conjugate'), index, mirror)
        covariance = channel.apply_covariance(np.eye(8)/2)
        covariances.append(covariance[perm][:, perm])
    return intensity_difference_spectrum(covariances, result['output_carrier_amplitudes_sqrt_flux'], detector,
        mode_labels=SIDEBAND_MODES, analysis_axis=detector.analysis_axis)


def build_report():
    with core.blas_single_thread():
        return _build_report()


def _build_report():
    base = conditional_inputs()
    det = DetectorResponse(AnalysisFrequencyAxis.from_hz([.1e6, .5e6, 1e6, 2e6, 3e6, 4e6]),
        [.85, .85], np.ones((6, 2)), 1., np.zeros(6), 'assumed flat balanced detector, no fitted loss')
    signed_rf = AnalysisFrequencyAxis(np.unique(np.r_[-det.analysis_axis.omega_rad_s, 0., det.analysis_axis.omega_rad_s]))
    grids = [reduced_periodic_cell(base, det, segments=n) for n in (4, 8, 16)]
    spatial_changes = [float(np.max(abs(b['spectrum'].quantum_db-a['spectrum'].quantum_db))) for a, b in zip(grids[:-1], grids[1:])]
    adaptive = [reduced_periodic_cell(base, det, segments=4, propagation_rtol=tol) for tol in (2e-9, 2e-11)]
    adaptive_change = float(np.max(abs(adaptive[-1]['spectrum'].quantum_db-adaptive[0]['spectrum'].quantum_db)))
    results, cases = [], []
    for power in (1e-8, 8e-6, 1e-3):
        inputs = replace(base, seed_power_W=power)
        result = adaptive[-1] if power == 8e-6 else reduced_periodic_cell(inputs, det, segments=4, propagation_rtol=2e-11)
        results.append(result)
        old = power_normalized_readout(inputs, signed_rf, det)
        pinched = _pinched_readout(result, inputs, det)
        greater, _ = result['transfer'].vacuum_output()
        reference_current = []
        for j, w in enumerate(det.analysis_axis.omega_rad_s):
            idx = np.flatnonzero(result['transfer'].frequency_axis.omega_rad_s == w)[0]
            reference_current.append(ordered_current_psd(greater[idx], result['output_carrier_amplitudes_sqrt_flux'],
                det.transmissions, det.current_response[j], det.balance, c.ELEMENTARY_CHARGE))
        current_error = _relative(reference_current, result['spectrum'].quantum_psd_A2_Hz)
        cases.append({'seed_power_W': power, 'probe_power_gain': result['probe_power_gain'],
            'conjugate_power_gain': result['conjugate_power_gain'],
            'pump_state_probe_power_gain': old.probe_power_gain,
            'detected_S_minus_db': result['spectrum'].quantum_db.tolist(),
            'pump_state_S_minus_db': old.spectrum.quantum_db.tolist(),
            'cross_sector_pinched_S_minus_db': pinched.quantum_db.tolist(),
            'maximum_change_from_pump_state_db': float(np.max(abs(result['spectrum'].quantum_db-old.spectrum.quantum_db))),
            'maximum_cross_sector_pinching_change_db': float(np.max(abs(pinched.quantum_db-result['spectrum'].quantum_db))),
            'quantum_psd_A2_Hz': result['spectrum'].quantum_psd_A2_Hz.tolist(),
            'sql_psd_A2_Hz': result['spectrum'].sql_psd_A2_Hz.tolist(),
            'output_coherent_powers_W': result['output_coherent_powers_W'].tolist(),
            'output_carriers_sqrt_flux': _complex(result['output_carrier_amplitudes_sqrt_flux']),
            'midpoint_z_m': result['midpoint_z_m'].tolist(),
            'midpoint_carriers_sqrt_flux': _complex(result['midpoint_carriers_sqrt_flux']),
            'quantum_propagation': result['quantum_propagation'],
            'maximum_local_commutator_relative_residual': max(g.audit()['maximum_commutator_relative_residual'] for g in result['local_generators']),
            'global_audit': result['transfer'].audit(),
            'minimum_channel_cp_eigenvalue': min(ch.audit().minimum_cp_eigenvalue for ch in result['channels']),
            'minimum_detected_covariance_uncertainty_eigenvalue': result['spectrum'].minimum_covariance_uncertainty_eigenvalue,
            'independent_ordered_photocurrent_relative_error': current_error,
            'local_drift_m_inverse': _complex([g.drift for g in result['local_generators']]),
            'local_noise_greater_by_reservoir_m_inverse': _complex([g.noise_greater_by_reservoir for g in result['local_generators']]),
            'local_noise_lesser_by_reservoir_m_inverse': _complex([g.noise_lesser_by_reservoir for g in result['local_generators']]),
            'transfer': _complex(result['transfer'].transfer),
            'added_noise_greater': _complex(result['transfer'].noise_greater),
            'added_noise_lesser': _complex(result['transfer'].noise_lesser),
            'passed': bool(current_error < 1e-10 and result['transfer'].audit()['passed'])})
    # Independently verify M and D at a finite-seed midpoint, including RF=0.
    selected = results[1]
    atom = reduced_periodic_noise(base, selected['midpoint_carriers_sqrt_flux'][2], mean_order=4)
    ports = reduced_periodic_ports(base)
    axis = GeneratorFrequencyAxis([0., 2*np.pi*1e6], 'reference-Floquet-quasifrequency')
    density = base.number_density_m3*base.uniform_area_m2
    hs, ops, g, signs, _ = ports.nambu_coordinates()
    args = (atom.reservoirs.generator(atom.hamiltonian_zero_rad_s), core.comm_super(atom.hamiltonian_plus_rad_s),
            core.comm_super(atom.hamiltonian_plus_rad_s.conj().T), atom.beat_rad_s)
    candidates = [eliminate_periodic_atom(atom, ports, axis, response_order=n, linear_density_m_inverse=density) for n in (2, 3, 4)]
    response_changes = [_relative(b.drift, a.drift) for a, b in zip(candidates[:-1], candidates[1:])]
    noise_changes = [_relative(b.noise_greater, a.noise_greater) for a, b in zip(candidates[:-1], candidates[1:])]
    try:
        eliminate_periodic_atom(atom, ports, axis, response_order=1, linear_density_m_inverse=density)
        insufficient_order = {'rejected': False}
    except ValueError as exc:
        insufficient_order = {'rejected': True, 'message': str(exc)}
    ref_m = forced_liouville_response(*args, atom.state_harmonics, ops, hs, g, signs, axis.omega_rad_s,
        response_order=4, linear_density=density)
    qrt = time_domain_qrt(*args, atom.operators, axis.omega_rad_s, [-1, 1], phase_samples=16, rtol=2e-12, atol=2e-14)
    readout = np.zeros((4, 30), complex)
    for j, (h, op, coupling, sign) in enumerate(zip(hs, ops, g, signs)):
        start = 0 if h == -1 else 15
        readout[j, start:start+15] = -1j*sign*coupling*np.einsum('ab,kba->k', op, atom.operators)
    ref_d = density*(readout@qrt['ordered_spectrum']@readout.conj().T)
    m_error, d_error = _relative(candidates[-1].drift, ref_m), _relative(candidates[-1].noise_greater, ref_d)
    # DC quantum transfer must be the derivative of the nonlinear mean map.
    beta0 = np.array([np.sqrt(base.seed_power_W/(c.HBAR*optical_carriers(base.detunings)[1])), 0.])
    t0 = selected['transfer'].transfer[np.flatnonzero(selected['transfer'].frequency_axis.omega_rad_s == 0)[0]]
    tangent_errors = []
    for j in range(2):
        for phase in (1., 1j):
            u = np.zeros(2, complex); u[j] = phase
            eps = np.linalg.norm(beta0)*1e-4
            outputs = []
            for sign in (-1, 1):
                ode, scale = integrate_reduced_carriers(base, initial_amplitudes=beta0+sign*eps*u)
                outputs.append(ode.y[:, -1]*scale)
            fd = (outputs[1]-outputs[0])/(2*eps)
            tangent_errors.append(_relative(t0[:2]@np.r_[u, u.conj()], fd))
    wrong_mean_error = _relative(t0[:2]@np.r_[beta0, beta0.conj()], selected['output_carrier_amplitudes_sqrt_flux'])
    sources = [Path(__file__), ROOT/'analysis/grand_challenge/normalization_audit.py',
        ROOT/'analysis/grand_challenge/reference/periodic_field.py', ROOT/'analysis/grand_challenge/reference/periodic_qrt.py']
    sources += list((ROOT/'gabes/quantum').glob('*.py'))+list((ROOT/'gabes/fwm_quantum').glob('*.py'))
    sources += [ROOT/'gabes'/name for name in ('core.py', 'atoms.py', 'constants.py', 'observables.py', 'species.py', 'hyperfine.py', 'zeeman.py', 'plot_style.py', 'schemes/fwm.py')]
    strong_loose = reduced_periodic_cell(replace(base, seed_power_W=1e-3), det, segments=4, propagation_rtol=2e-9)
    strong_change = float(np.max(abs(strong_loose['spectrum'].quantum_db-results[-1]['spectrum'].quantum_db)))
    passed = all(case['passed'] for case in cases) and max(adaptive_change, strong_change) < 1e-7 and m_error < 1e-9 and d_error < 1e-8 and max(tangent_errors) < 1e-6 and insufficient_order['rejected']
    return {'stage': 'S1 finite-seed microscopic field, nonlinear two-carrier mean and bright readout',
        'expected_controls_passed': bool(passed), 'finite_seed_two_band_field_implemented': True,
        'absolute_hot_vapor_prediction': False, 'experimental_validation': False,
        'fixed_classical_pump': True, 'additional_physical_optical_ports_propagated': False,
        'rf_hz': det.analysis_axis.frequency_hz.tolist(),
        'signed_quasifrequencies_rad_s': selected['transfer'].frequency_axis.omega_rad_s.tolist(),
        'field_nambu_order': list(selected['transfer'].mode_labels),
        'field_carrier_harmonics': hs.tolist(), 'field_nambu_signs': signs.tolist(),
        'mean_order': 4, 'response_order': 3, 'quantum_propagation_method': 'adaptive-DOP853', 'local_diagnostic_midpoint_count': 4,
        'base_input_ledger': {k: {'value': getattr(base, k), 'unit': unit, 'status': 'assumed'} for k, unit in INPUT_UNITS.items()},
        'detector': {'transmissions': det.transmissions.tolist(), 'balance': det.balance,
                     'current_response': _complex(det.current_response), 'source': det.source},
        'configurations': cases,
        'spatial_orders': [4, 8, 16], 'adjacent_spatial_maximum_S_minus_change_db': spatial_changes,
        'initial_midpoint_convergence_passed': bool(spatial_changes[-1] < 1e-7),
        'adaptive_tolerances': [2e-9, 2e-11],
        'adaptive_tolerance_change_db_8uW': adaptive_change,
        'adaptive_tolerance_change_db_1mW': strong_change,
        'midpoint_to_adaptive_maximum_S_minus_errors_db': [float(np.max(abs(g['spectrum'].quantum_db-adaptive[-1]['spectrum'].quantum_db))) for g in grids],
        'atomic_response_orders': [2, 3, 4], 'adjacent_response_relative_drift_changes': response_changes,
        'adjacent_response_relative_noise_changes': noise_changes,
        'insufficient_atomic_order_control': insufficient_order,
        'independent_full_liouville_drift_relative_error': m_error,
        'independent_time_qrt_field_noise_relative_error': d_error,
        'nonlinear_mean_tangent_maximum_relative_error': max(tangent_errors),
        'incorrect_transfer_times_mean_relative_error': wrong_mean_error,
        'source_sha256': {str(p.relative_to(ROOT)).replace('\\', '/'): hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(set(sources))},
        'environment': {'python': platform.python_version(), 'numpy': np.__version__, 'scipy': scipy.__version__},
        'sources': ['https://arxiv.org/abs/2005.08249', 'https://arxiv.org/abs/2301.11993'],
        'limits': ['Same assumed reduced RMS dipoles and reservoirs, single velocity, zero phase mismatch, uniform area.',
            'Only the declared probe/conjugate optical bands propagate. Atomic harmonic convergence does not establish optical-mode closure.',
            'Probe/conjugate means evolve nonlinearly while the classical pump remains externally prescribed; pump depletion and quantum pump ports remain open.',
            'Bright linear photocurrent: previous pump-state Gaussian quadratic correction is not silently transplanted to this periodic field.',
            'Cross-sector pinching is a diagnostic alteration, not an accepted microscopic prediction or fitted noise model.',
            'No full Zeeman, angular-Doppler/transverse modes, independent measured input uncertainties or held-out experimental validation.']}


def save_plot(report, path):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from gabes.plot_style import PALETTE, apply_gabes_plot_style
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.6), layout='constrained')
    f = np.asarray(report['rf_hz'])/1e6
    baseline = report['configurations'][1]
    axes[0].plot(f, baseline['pump_state_S_minus_db'], '--', color=PALETTE['muted'], label='Pump-state baseline')
    for case, label, color in zip(report['configurations'], ('10 nW', '8 uW', '1 mW'), (PALETTE['cyan'], PALETTE['rose'], PALETTE['warm'])):
        axes[0].plot(f, case['detected_S_minus_db'], 'o-', color=color, label=label)
    axes[0].set(xlabel='Analysis frequency [MHz]', ylabel='Detected S_minus [dB]', title='Full four-coordinate Nambu propagation')
    axes[0].legend(fontsize=8)
    for name, label, color in (('detected_S_minus_db', 'Full periodic field', PALETTE['rose']),
                               ('cross_sector_pinched_S_minus_db', 'Cross sectors removed', PALETTE['cyan'])):
        axes[1].plot(f, np.array(baseline[name])-baseline['pump_state_S_minus_db'], 'o-', color=color, label=label)
    axes[1].set(xlabel='Analysis frequency [MHz]', ylabel='Change from pump-state result [dB]', title='8 uW: effect on intensity difference')
    axes[1].legend(fontsize=8)
    fig.suptitle('Conditional two-band finite-seed field\nFixed pump; bright readout; no hot-vapor validation')
    apply_gabes_plot_style(fig)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('xb') as stream:
        fig.savefig(stream, format='png', dpi=170)
    plt.close(fig)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    parser.add_argument('--plot', type=Path)
    args = parser.parse_args(argv)
    for path in (args.output, args.plot):
        if path is not None and path.exists():
            raise FileExistsError(path)
    report = build_report()
    payload = json.dumps(report, ensure_ascii=False, indent=2, allow_nan=False)+'\n'
    if args.output is None:
        print(payload, end='')
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        with args.output.open('x', encoding='utf-8') as stream:
            stream.write(payload)
        print(f'Wrote {args.output}')
    if args.plot is not None:
        save_plot(report, args.plot)
        print(f'Wrote {args.plot}')
    return 0 if report['expected_controls_passed'] else 1


if __name__ == '__main__':
    raise SystemExit(main())
