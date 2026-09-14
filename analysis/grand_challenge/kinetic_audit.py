"""Conditional carrier-closure, velocity-noise and thermal pump-state audit.

python -m analysis.grand_challenge.kinetic_audit --output NEW.json --plot NEW.png
"""

import argparse
from dataclasses import replace
import hashlib
import json
from pathlib import Path
import platform
import time

import numpy as np
import scipy

from gabes import core, constants as c
from gabes.fwm_quantum.inputs import INPUT_UNITS
from gabes.fwm_quantum.kinetic import CarrierGeometry, VelocityQuadrature, kinetic_local_field, thermal_pump_cell
from gabes.fwm_quantum.model import reduced_pump_system
from gabes.fwm_quantum.normalization import optical_carriers, reduced_dipoles
from gabes.quantum.contracts import AnalysisFrequencyAxis
from gabes.quantum.readout import DetectorResponse
from .normalization_audit import conditional_inputs
from .reference.kinetic import pump_state_field_reference


ROOT = Path(__file__).resolve().parents[2]
TEMPERATURE_K = 394.15
# Tolerances apply to the declared reduced model, not an experimental error bar.
TOLERANCES = {'relative_M': 1e-5, 'relative_D_greater': 1e-5,
    'relative_D_lesser': 1e-5, 'relative_gain': 1e-7, 'maximum_S_minus_change_db': 1e-6}


def _complex(a):
    a = np.asarray(a)
    return {'real': a.real.tolist(), 'imag': a.imag.tolist()}


def _relative(a, b):
    return float(np.linalg.norm(a-b)/max(np.linalg.norm(b), np.finfo(float).tiny))


def _worst_frequency(a, b):
    return float(np.max(np.linalg.norm(a-b, axis=(-2, -1))/
        np.maximum(np.linalg.norm(b, axis=(-2, -1)), np.finfo(float).tiny)))


def _comparison(a, b):
    ga, gb = a['local']['generator'], b['local']['generator']
    values = dict(zip(('relative_M', 'relative_D_greater', 'relative_D_lesser'),
        (_worst_frequency(getattr(ga, key), getattr(gb, key)) for key in ('drift', 'noise_greater', 'noise_lesser'))))
    values['relative_gain'] = abs(a['probe_power_gain']/b['probe_power_gain']-1)
    values['maximum_S_minus_change_db'] = float(np.max(abs(a['spectrum'].quantum_db-b['spectrum'].quantum_db)))
    values['passed'] = bool(all(values[key] < limit for key, limit in TOLERANCES.items()))
    return values


def _case(result):
    local = result['local']
    return {'probe_power_gain': result['probe_power_gain'], 'conjugate_power_gain': result['conjugate_power_gain'],
        'detected_S_minus_db': result['spectrum'].quantum_db.tolist(),
        'quantum_psd_A2_Hz': result['spectrum'].quantum_psd_A2_Hz.tolist(),
        'sql_psd_A2_Hz': result['spectrum'].sql_psd_A2_Hz.tolist(),
        'local_audit': local['generator'].audit(), 'global_audit': result['transfer'].audit(),
        'maximum_class_commutator_relative_residual': local['maximum_class_commutator_relative_residual'],
        'minimum_channel_cp_eigenvalue': result['minimum_channel_cp_eigenvalue'],
        'minimum_detected_covariance_uncertainty_eigenvalue': result['spectrum'].minimum_covariance_uncertainty_eigenvalue}


def _local_controls(inputs, geometry):
    axis = AnalysisFrequencyAxis.from_hz([-1e6, 0., 1e6])
    q = VelocityQuadrature([[150., 0., 380.], [-110., 0., -220.], [0., 0., 700.]], [.4, .5, .1],
        'three discrete reference velocities; includes a near-resonant atom, not a Maxwell approximation')
    actual = kinetic_local_field(inputs, geometry, q, axis, retain_classes=False)
    hs, ops, gs, signs, _ = actual['ports'].nambu_coordinates()
    dipoles = reduced_dipoles('uniform-zeeman-rms')
    totals = {key: np.zeros_like(actual['generator'].drift) for key in ('drift', 'greater', 'lesser')}
    for weight, delta, beats in zip(q.probabilities, actual['atomic_one_photon_rad_s'], actual['atomic_signed_beats_rad_s']):
        h, reservoirs = reduced_pump_system(dipoles.pump_rabi_rad_s(inputs.pump_power_W, inputs.pump_waist_m), delta,
            transit_rate_s_inverse=inputs.transit_rate_s_inverse, transition_scales=dipoles.transition_scales)
        ref = pump_state_field_reference(reservoirs.generator(h), ops, gs, signs, hs, beats[0],
            axis.omega_rad_s, inputs.number_density_m3*inputs.uniform_area_m2*weight)
        for key in totals:
            totals[key] += ref[key]
    errors = {name: _relative(getattr(actual['generator'], key), totals[name])
        for name, key in (('drift', 'drift'), ('greater', 'noise_greater'), ('lesser', 'noise_lesser'))}
    one = VelocityQuadrature([[90., 0., 230.]], [1.], 'one velocity class')
    split = VelocityQuadrature([[90., 0., 230.], [90., 0., 230.]], [.25, .75], 'identical class split')
    a = kinetic_local_field(inputs, geometry, one, axis)['generator']
    b = kinetic_local_field(inputs, geometry, split, axis, retain_classes=False)['generator']
    split_errors = {key: _relative(getattr(a, key), getattr(b, key)) for key in ('drift', 'noise_greater', 'noise_lesser')}
    wrong_factor = (np.sqrt(.25)+np.sqrt(.75))**2
    wrong = replace(a, noise_greater_by_reservoir=a.noise_greater_by_reservoir*wrong_factor,
        noise_lesser_by_reservoir=a.noise_lesser_by_reservoir*wrong_factor).audit()
    finite_q = VelocityQuadrature([[90., 0., 230.], [-160., 0., -120.]], [.4, .6], 'two discrete finite-seed reference velocities')
    beta = np.array([3e6+1e6j, .2e6j])
    direction = np.array([1+.2j, -.3+.7j])
    finite = kinetic_local_field(inputs, geometry, finite_q, axis, carrier_amplitudes=beta, retain_classes=False)
    derivatives = []
    for sign in (-1, 1):
        shifted = kinetic_local_field(inputs, geometry, finite_q, axis, carrier_amplitudes=beta+sign*100*direction,
            retain_classes=False)
        derivatives.append(-1j*inputs.number_density_m3*inputs.uniform_area_m2*
            shifted['ports'].coupling_s_inverse_sqrt_flux*shifted['mean_polarization'])
    dc = np.flatnonzero(axis.omega_rad_s == 0)[0]
    tangent = (finite['generator'].drift[dc]@np.r_[direction, direction.conj()])[:2]
    mean_error = _relative(tangent, (derivatives[1]-derivatives[0])/200)
    return {'independent_full_liouville_qrt_relative_errors': errors,
        'reference_velocities_m_s': q.velocities_m_s.tolist(), 'reference_probabilities': q.probabilities.tolist(),
        'split_class_relative_errors': split_errors, 'incorrect_noise_amplitude_sum_audit': wrong,
        'finite_seed_local_mean_tangent_relative_error': mean_error,
        'finite_seed_local_velocities_m_s': finite_q.velocities_m_s.tolist(),
        'finite_seed_local_probabilities': finite_q.probabilities.tolist(),
        'finite_seed_local_carriers_sqrt_flux': _complex(beta),
        'finite_seed_local_mean_polarization': _complex(finite['mean_polarization']),
        'finite_seed_local_audit': finite['generator'].audit(),
        'passed': bool(max(errors.values()) < 1e-9 and max(split_errors.values()) < 1e-11
            and not wrong['passed'] and mean_error < 1e-6 and finite['generator'].audit()['passed'])}


def build_report():
    with core.blas_single_thread():
        return _build_report()


def _build_report():
    inputs = conditional_inputs()
    detector = DetectorResponse(AnalysisFrequencyAxis.from_hz([.1e6, 1e6, 4e6]), [.85, .85], np.ones((3, 2)),
        1., np.zeros(3), 'assumed flat balanced detector; no fitted loss')
    vacuum = CarrierGeometry.vacuum_beams(inputs, probe_angle_rad=.005, conjugate_angle_rad=-.005)
    k0, kp, _ = vacuum.wavevectors_rad_m
    geometry = CarrierGeometry([k0, kp, 2*k0-kp],
        'synthetic non-vacuum conjugate defined explicitly as 2*k0-kp; diagnostic fixture, no measured dispersion')
    grids = [(64, 4, 5), (128, 16, 5), (256, 32, 5), (512, 32, 7), (768, 32, 7),
             (1024, 32, 7), (1024, 48, 7), (1024, 48, 6), (1024, 48, 8)]
    results, cases = [], []
    for nz, nx, cut in grids:
        started = time.monotonic()
        quadrature = VelocityQuadrature.maxwell_xz(TEMPERATURE_K, longitudinal_order=nz, transverse_order=nx, cutoff_sigma=cut)
        result = thermal_pump_cell(inputs, geometry, quadrature, detector)
        results.append(result)
        cases.append({'grid': {'longitudinal_order': nz, 'transverse_order': nx, 'cutoff_sigma': cut},
            'omitted_Maxwell_probability': quadrature.omitted_probability,
            'seconds': time.monotonic()-started, **_case(result)})
        print(f'Completed velocity grid {nz} x {nx}, cutoff {cut} sigma', flush=True)
    chosen_index = 6
    selected = results[chosen_index]
    refinements = [{'name': name, 'from_case': a, 'to_case': b, **_comparison(results[a], results[b])}
        for name, a, b in (('longitudinal', 4, 5), ('transverse', 5, 6), ('tail_6_to_7_sigma', 7, 6), ('tail_7_to_8_sigma', 6, 8))]
    collinear = CarrierGeometry.vacuum_beams(inputs)
    cold = thermal_pump_cell(inputs, collinear, VelocityQuadrature([[0., 0., 0.]], [1.], 'single atom at rest'), detector)
    hot = thermal_pump_cell(inputs, collinear,
        VelocityQuadrature.maxwell_xz(TEMPERATURE_K, longitudinal_order=1024, transverse_order=2, cutoff_sigma=7), detector)
    controls = _local_controls(inputs, geometry)
    closure = vacuum.closure_audit(temperature_K=TEMPERATURE_K, length_m=inputs.length_m)
    rejected = False
    try:
        vacuum.require_closed()
    except ValueError:
        rejected = True
    omega = optical_carriers(inputs.detunings)
    gen = selected['local']['generator']
    sources = [Path(__file__), ROOT/'analysis/grand_challenge/normalization_audit.py',
        ROOT/'analysis/grand_challenge/reference/kinetic.py', ROOT/'analysis/grand_challenge/reference/atomic_qrt.py']
    sources += list((ROOT/'gabes/quantum').glob('*.py'))+list((ROOT/'gabes/fwm_quantum').glob('*.py'))
    sources += [ROOT/'gabes'/name for name in ('core.py', 'atoms.py', 'constants.py', 'observables.py', 'species.py',
        'hyperfine.py', 'zeeman.py', 'doppler.py', 'plot_style.py', 'schemes/fwm.py')]
    return {'stage': 'S1 kinetic local mean/noise and conditional thermal pump-state cell',
        'expected_controls_passed': bool(controls['passed'] and rejected and all(row['passed'] for row in refinements)
            and all(row['local_audit']['passed'] and row['global_audit']['passed']
                and row['minimum_channel_cp_eigenvalue'] >= -1e-10
                and row['minimum_detected_covariance_uncertainty_eigenvalue'] >= -1e-10 for row in cases)),
        'finite_seed_kinetic_local_model_implemented': True, 'thermal_pump_state_cell_implemented': True,
        'finite_seed_thermal_propagation_implemented': False, 'general_noncollinear_geometry_implemented': False,
        'absolute_hot_vapor_prediction': False, 'experimental_validation': False,
        'rf_hz': detector.analysis_axis.frequency_hz.tolist(), 'rf_scope': 'three samples, not a bandwidth or minimum search',
        'base_input_ledger': {key: {'value': getattr(inputs, key), 'unit': unit, 'status': 'assumed'} for key, unit in INPUT_UNITS.items()},
        'kinetic_input_ledger': {'temperature_K': TEMPERATURE_K, 'probe_angle_rad': .005, 'status': 'assumed',
            'density_rule': 'independent supplied number density; kinetic T does not infer density or collision rate'},
        'detector': {'transmissions': detector.transmissions.tolist(), 'balance': detector.balance,
            'current_response': _complex(detector.current_response), 'source': detector.source},
        'vacuum_opposite_angle_geometry': {'wavevectors_rad_m': vacuum.wavevectors_rad_m.tolist(),
            'closure_audit': closure, 'rejected_by_single_phase_solver': rejected},
        'synthetic_closed_geometry': {'wavevectors_rad_m': geometry.wavevectors_rad_m.tolist(), 'source': geometry.source,
            'closure_audit': geometry.closure_audit(temperature_K=TEMPERATURE_K, length_m=inputs.length_m),
            'wavevector_magnitude_over_vacuum': (np.linalg.norm(geometry.wavevectors_rad_m, axis=1)*c.C_LIGHT/omega).tolist(),
            'normalization': 'vacuum photon energy/flux convention through common z-plane area; g_j divided by sqrt(cos(theta_j)); refractive energy normalization not certified'},
        'lab_optical_omega_rad_s': omega.tolist(), 'lab_signed_beat_rad_s': selected['local']['lab_signed_beat_rad_s'],
        'mean_order_for_finite_seed_local_fixture': 4, 'response_order_for_finite_seed_local_fixture': 3,
        'velocity_grid_cases': cases, 'selected_case': chosen_index,
        'coarse_to_selected_comparisons': [_comparison(result, selected) for result in results[:4]],
        'numerical_tolerances': TOLERANCES, 'refinements': refinements, 'local_controls': controls,
        'single_velocity_collinear_pump_state': _case(cold), 'thermal_collinear_pump_state': _case(hot),
        'selected_generator': {'signed_lab_rf_rad_s': gen.frequency_axis.omega_rad_s.tolist(),
            'mode_labels': list(gen.mode_labels), 'signs': gen.signs.tolist(), 'drift_m_inverse': _complex(gen.drift),
            'reservoir_names': list(gen.reservoir_names),
            'noise_greater_by_reservoir_m_inverse': _complex(gen.noise_greater_by_reservoir),
            'noise_lesser_by_reservoir_m_inverse': _complex(gen.noise_lesser_by_reservoir)},
        'source_sha256': {str(p.relative_to(ROOT)).replace('\\', '/'): hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(set(sources))},
        'environment': {'python': platform.python_version(), 'numpy': np.__version__, 'scipy': scipy.__version__},
        'sources': ['https://arxiv.org/abs/1303.7187', 'https://arxiv.org/abs/1007.1610', 'https://arxiv.org/abs/2005.08249'],
        'limits': ['Nonzero carrier-loop mismatch requires additional spatial/convective phases; scalar mismatch or shifted lab frequencies do not close this solver.',
            'Synthetic closed wavevectors are a declared numerical fixture, not a measured refractive-index or Maxwell solution.',
            'Finite seed is validated locally on discrete velocity fixtures; thermal cell uses the exact pump-state weak-field limit and fixed pump.',
            'Local ballistic response neglects convection of slow spatial envelopes and inter-slice atom transport correlations; no velocity-changing collisions.',
            'Conditional Maxwell truncation probability is not a rigorous bound on a resonantly weighted spectrum. Refinement is a numerical comparison, not an experimental error bar.',
            'Two optical bands, RMS reduced dipoles and explicit radiative/internal-reset reservoirs; no full Zeeman, transverse modes, pump depletion or periodic higher-order photocurrent.',
            'No independently measured geometry/input uncertainty or out-of-sample gain/squeezing validation.']}


def save_plot(report, path):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from gabes.plot_style import PALETTE, apply_gabes_plot_style
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.6), layout='constrained')
    f = np.asarray(report['rf_hz'])/1e6
    chosen = report['velocity_grid_cases'][report['selected_case']]
    for case, label, color in ((report['single_velocity_collinear_pump_state'], 'Atom at rest', PALETTE['muted']),
        (report['thermal_collinear_pump_state'], '394 K, collinear', PALETTE['cyan']),
        (chosen, '394 K, synthetic closed k', PALETTE['rose'])):
        axes[0].plot(f, case['detected_S_minus_db'], 'o-', color=color, label=label)
    axes[0].set(xlabel='Analysis frequency [MHz]', ylabel='Detected S_minus [dB]', title='Pump-state weak-field samples')
    axes[0].legend(fontsize=8)
    rows = report['refinements']
    x = np.arange(len(rows))
    for key, label, color in (('relative_M', 'M / tolerance', PALETTE['cyan']),
        ('relative_D_greater', 'D greater / tolerance', PALETTE['rose']),
        ('maximum_S_minus_change_db', 'S_minus / tolerance', PALETTE['warm'])):
        axes[1].semilogy(x, [max(row[key]/report['numerical_tolerances'][key], 1e-12) for row in rows], 'o-', color=color, label=label)
    axes[1].axhline(1., ls='--', color=PALETTE['muted'])
    axes[1].set(xticks=x, xticklabels=['N_z', 'N_x', '6 to 7 sigma', '7 to 8 sigma'],
        ylabel='Numerical change / declared tolerance', title='Independent quadrature refinements')
    axes[1].legend(fontsize=8)
    fig.suptitle('Conditional kinetic audit\nFixed pump; synthetic closed wavevectors; no experimental validation')
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
