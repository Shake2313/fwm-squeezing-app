"""Immutable thermal boundary-flux and common-phase Poisson audit.

python -m analysis.grand_challenge.inflow_audit --output NEW.json --plot NEW.png
"""

import argparse
import hashlib
import json
from pathlib import Path
import platform
import time

import numpy as np
import scipy
from scipy.special import ndtri

from gabes import core
from gabes.constants import KB, MASS_85RB
from gabes.quantum.contracts import GeneratorFrequencyAxis
from gabes.quantum.inflow import (maxwell_box_inflow, uniform_entry_phases,
    entry_phase_factors, average_marked_poisson)
from gabes.quantum.reservoirs import ExplicitReservoirs
from gabes.quantum.transport import integrate_characteristic


ROOT = Path(__file__).resolve().parents[2]
POWERS = (10, 12, 14, 16, 18)
SEEDS = (11, 211, 811)
TEMPERATURE_K = 373.
DENSITY_M3 = 1e17
SHAPES_M = {'cube': np.full(3, np.cbrt(2)*1e-3), 'thin_box': np.array([4., 2., .25])*1e-3}
OMEGA = 2*np.pi*np.array([0., 5e4, 1.5e5, 5e5])
TOLERANCES = {'final_phase_space_moment_absolute': .002,
    'identity_spectrum_last_refinement_relative': .015,
    'phase_raw_identity_relative': 2e-12, 'actual_transport_identity_relative': 2e-9}


def hashes():
    files = ['gabes/core.py', 'gabes/constants.py', 'gabes/quantum/__init__.py']
    files += ['gabes/quantum/'+name+'.py' for name in
              ('inflow', 'transport', 'contracts', 'diffusion', 'reservoirs', 'periodic')]
    files += ['analysis/grand_challenge/inflow_audit.py', 'tests/quantum/test_inflow.py']
    return {name: hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in files}


def relative(actual, expected):
    return float(np.linalg.norm(actual-expected)/max(np.linalg.norm(expected), np.finfo(float).tiny))


def make_inflow(lengths, power, seed):
    return maxwell_box_inflow(-lengths/2, lengths/2, temperature_K=TEMPERATURE_K,
        mass_kg=MASS_85RB, density_m3=DENSITY_M3, points_per_face_power=power, seed=seed,
        source='Synthetic box dimensions, 373 K and n=1e17 m^-3; explicit existing MASS_85RB; no measured apparatus')


def pulse(tau, omega=OMEGA):
    return tau[..., None]*np.sinc(tau[..., None]*omega/(2*np.pi))*np.exp(.5j*tau[..., None]*omega)


def moments_row(inflow):
    m = inflow.occupation_moments()
    sigma = np.sqrt(KB*inflow.temperature_K/inflow.mass_kg)
    lengths = inflow.upper_corner_m-inflow.lower_corner_m
    position_scale = lengths/np.sqrt(12)
    errors = {'occupancy': abs(m['occupancy_over_nV']-1),
        'velocity_mean': float(np.max(np.abs(m['velocity_mean_m_s']))/sigma),
        'velocity_second': float(np.max(np.abs(m['velocity_second_m2_s2']/sigma**2-np.eye(3)))),
        'mean_speed': abs(m['speed_mean_m_s']/(np.sqrt(8/np.pi)*sigma)-1),
        'speed_fourth': abs(m['speed_fourth_m4_s4']/(15*sigma**4)-1),
        'position_mean': float(np.max(np.abs(m['position_centered_mean_m']/position_scale))),
        'position_second': float(np.max(np.abs(m['position_centered_second_m2']
            /np.outer(position_scale, position_scale)-np.eye(3))))}
    face_rates = np.bincount(inflow.face_index, weights=inflow.rate_s_inverse, minlength=6)
    analytic_face = np.repeat(DENSITY_M3*np.prod(lengths)/lengths*sigma/np.sqrt(2*np.pi), 2)
    spectra = np.sum(inflow.rate_s_inverse[:, None]*np.abs(pulse(inflow.residence_time_s))**2, axis=0)
    return {'power': inflow.points_per_face_power, 'seed': inflow.seed,
        'path_count': len(inflow.rate_s_inverse), 'mean_occupancy': inflow.mean_occupancy,
        'equilibrium_nV': inflow.equilibrium_atom_number, 'moment_errors': errors,
        'max_moment_error': max(errors.values()),
        'face_rates_s_inverse': face_rates.tolist(), 'analytic_face_rates_s_inverse': analytic_face.tolist(),
        'face_rate_max_relative_error': float(np.max(np.abs(face_rates/analytic_face-1))),
        'total_rate_s_inverse': inflow.total_arrival_rate_s_inverse,
        'mean_residence_s': inflow.mean_occupancy/inflow.total_arrival_rate_s_inverse,
        'analytic_mean_residence_s': inflow.equilibrium_atom_number/analytic_face.sum(),
        'identity_number_spectrum_s': spectra.tolist()}


def volume_velocity_negative_control(inflow):
    """Replace only Rayleigh normal speeds by half-normal volume speeds.

    Retain exact face fluxes, tangential velocities and entry points. Recover
    the underlying uniform coordinate by Rayleigh CDF, then inverse half-normal
    CDF. This intentionally wrong measure cannot hide behind a fitted rate.
    """
    velocity = np.array(inflow.velocity_m_s)
    rows = np.arange(len(velocity))
    normal = inflow.face_index//2
    sigma = np.sqrt(KB*inflow.temperature_K/inflow.mass_kg)
    speed = np.abs(velocity[rows, normal])/sigma
    u = -np.expm1(-speed**2/2)
    velocity[rows, normal] = np.sign(velocity[rows, normal])*sigma*ndtri((1+u)/2)
    # Independent slab intersection expressed as entry-to-exit coordinates.
    target = np.where(velocity > 0, inflow.upper_corner_m, inflow.lower_corner_m)
    tau = np.min((target-inflow.entry_position_m)/velocity, axis=1)
    weight = inflow.rate_s_inverse*tau/inflow.equilibrium_atom_number
    second = np.einsum('n,ni->i', weight, velocity**2)/sigma**2
    return {'wrong_occupancy_over_nV': float(weight.sum()),
        'wrong_velocity_second_over_sigma2': second.tolist(),
        'same_exact_face_rates': True, 'rejected': bool(abs(weight.sum()-1) > .05)}


def phase_poisson_control():
    inflow = make_inflow(SHAPES_M['thin_box'], 4, 11)
    axis = GeneratorFrequencyAxis(OMEGA, 'laboratory age Fourier frequency for synthetic identity pulses')
    wavevectors = np.array([[240., 100., 0.], [400., -50., 60.]])
    offsets = np.array([0., .35])
    def packet(path, phase):
        f = pulse(np.asarray(path.residence_time_s))
        amplitudes = entry_phase_factors(path, phase, wavevectors, [1, 1], offsets_rad=offsets).real
        diagonal = np.abs(f[:, None, None])**2*np.eye(2)
        return {'frequency_axis': axis, 'residence_time_s': path.residence_time_s,
            'mean_pulse': f[:, None]*amplitudes, 'greater': .2*diagonal, 'lesser': .07*diagonal}
    delta = inflow.entry_position_m@(wavevectors[0]-wavevectors[1])+offsets[0]-offsets[1]
    phase_outer = np.empty((len(delta), 2, 2))
    phase_outer[:, 0, 0] = phase_outer[:, 1, 1] = .5
    phase_outer[:, 0, 1] = phase_outer[:, 1, 0] = .5*np.cos(delta)
    squared = np.abs(pulse(inflow.residence_time_s))**2
    reference_number = np.einsum('n,nf,nij->fij', inflow.rate_s_inverse, squared, phase_outer)
    reference_internal = np.einsum('n,nf,ij->fij', inflow.rate_s_inverse, squared, np.eye(2))
    refinements = []
    for count in (1, 2, 4, 8, 16):
        result = average_marked_poisson(inflow, uniform_entry_phases(count, offset_rad=.173), packet,
            source='Synthetic cosine-modulated identity means, two dimensionless readouts, fixed wavevectors; declared connected matrices')
        refinements.append({'phase_count': count,
            'number_relative_error': relative(result['poisson_number'], reference_number),
            'greater_relative_error': relative(result['greater'], .2*reference_internal+reference_number),
            'lesser_relative_error': relative(result['lesser'], .07*reference_internal+reference_number)})
    mean = result['arrival_weighted_conditional_mean_pulse']
    wrong_mean_outer = inflow.total_arrival_rate_s_inverse*mean[:, :, None]*mean[:, None, :].conj()
    wrong_independent = np.array(reference_number)
    wrong_independent[:, 0, 1] = wrong_independent[:, 1, 0] = 0
    return {'phase_refinement': refinements,
        'ordering_difference_relative_error': relative(result['greater']-result['lesser'], .13*reference_internal),
        'wrong_average_mean_before_outer_relative_error': relative(wrong_mean_outer, reference_number),
        'wrong_independent_two_phase_relative_error': relative(wrong_independent, reference_number),
        'wrong_drop_poisson_mean_relative_error': relative(.2*reference_internal, result['greater']),
        'minimum_greater_eigenvalue_s': float(np.linalg.eigvalsh(result['greater']).min()),
        'dc_number_matrix_s': reference_number[0].tolist(),
        'passed': bool(max(max(row[key] for key in ('number_relative_error', 'greater_relative_error',
            'lesser_relative_error')) for row in refinements[2:]) < TOLERANCES['phase_raw_identity_relative']
            and relative(result['greater']-result['lesser'], .13*reference_internal) < TOLERANCES['phase_raw_identity_relative']
            and relative(wrong_mean_outer, reference_number) > .9
            and relative(wrong_independent, reference_number) > .5
            and relative(.2*reference_internal, result['greater']) > .5)}


def actual_transport_identity_control():
    """Use actual finite-atom transport.wavepacket on one chord per face."""
    inflow = make_inflow(SHAPES_M['thin_box'], 0, 11)
    axis = GeneratorFrequencyAxis(OMEGA[:2], 'laboratory age Fourier frequency for identity readout')
    cache = {}
    def packet(path, phase):
        if path.source not in cache:
            characteristic = integrate_characteristic(path, lambda a, r: np.zeros((2, 2)),
                ExplicitReservoirs(2, ()), np.diag([1., 0.]), times_s=[0, path.residence_time_s])
            cache[path.source] = characteristic.wavepacket(axis, lambda a, r: np.eye(2)[None])
        return cache[path.source]
    result = average_marked_poisson(inflow, uniform_entry_phases(4), packet,
        source='Actual transport.wavepacket with O=I, two-level constant state, H=0, no jumps; phase independent exact control')
    expected = np.sum(inflow.rate_s_inverse[:, None]*np.abs(pulse(inflow.residence_time_s, axis.omega_rad_s))**2, axis=0)
    error = relative(result['greater'][:, 0, 0], expected)
    return {'characteristics_solved': len(cache), 'frequency_hz': (axis.omega_rad_s/(2*np.pi)).tolist(),
        'expected_spectrum_s': expected.tolist(), 'actual_greater_s': result['greater'][:, 0, 0].real.tolist(),
        'relative_error': error, 'internal_greater_norm': float(np.linalg.norm(result['internal_greater'])),
        'passed': bool(error < TOLERANCES['actual_transport_identity_relative']
            and np.linalg.norm(result['internal_greater']) == 0)}


def build_report():
    started = time.monotonic()
    before = hashes()
    geometry = {}
    negatives = {}
    with core.blas_single_thread():
        for name, lengths in SHAPES_M.items():
            rows = []
            for power in POWERS:
                for seed in SEEDS:
                    inflow = make_inflow(lengths, power, seed)
                    rows.append(moments_row(inflow))
                print(f'{name}: six faces x 2**{power}, three independent scrambles complete', flush=True)
            negatives[name] = volume_velocity_negative_control(inflow)
            last = [row for row in rows if row['power'] == POWERS[-1]]
            previous = [row for row in rows if row['power'] == POWERS[-2]]
            spectral_changes = [float(np.max(np.abs(np.array(a['identity_number_spectrum_s'])
                /np.array(b['identity_number_spectrum_s'])-1))) for a, b in zip(previous, last)]
            independent_spectra = np.array([row['identity_number_spectrum_s'] for row in last])
            seed_spread = np.ptp(independent_spectra, axis=0)/independent_spectra.mean(axis=0)
            envelopes = {str(p): max(row['max_moment_error'] for row in rows if row['power'] == p) for p in POWERS}
            geometry[name] = {'lengths_m': lengths.tolist(), 'volume_m3': float(np.prod(lengths)),
                'surface_area_m2': float(2*np.sum(np.prod(lengths)/lengths)), 'rows': rows,
                'max_moment_error_by_power': envelopes,
                'identity_spectrum_last_refinement_relative': spectral_changes,
                'identity_spectrum_independent_seed_relative_spread_by_frequency': seed_spread.tolist(),
                'passed': bool(envelopes[str(POWERS[-1])] < TOLERANCES['final_phase_space_moment_absolute']
                    and envelopes[str(POWERS[-1])] < envelopes[str(POWERS[0])]/4
                    and max(spectral_changes) < TOLERANCES['identity_spectrum_last_refinement_relative']
                    and seed_spread.max() < TOLERANCES['identity_spectrum_last_refinement_relative'])}
        phase = phase_poisson_control()
        transport = actual_transport_identity_control()
    shape_flux_ratio = geometry['thin_box']['rows'][-1]['total_rate_s_inverse']/geometry['cube']['rows'][-1]['total_rate_s_inverse']
    area_ratio = geometry['thin_box']['surface_area_m2']/geometry['cube']['surface_area_m2']
    area = {'equal_volume_m3': float(np.prod(SHAPES_M['cube'])), 'arrival_rate_ratio': shape_flux_ratio,
        'surface_area_ratio': area_ratio, 'occupancy_ratio': geometry['thin_box']['rows'][-1]['mean_occupancy']
            /geometry['cube']['rows'][-1]['mean_occupancy'],
        'passed': bool(abs(shape_flux_ratio/area_ratio-1) < 1e-12
            and abs(geometry['thin_box']['rows'][-1]['mean_occupancy']
                /geometry['cube']['rows'][-1]['mean_occupancy']-1) < .002)}
    after = hashes()
    return {'schema': 'gabes-thermal-inflow-audit-v1',
        'scope': 'Synthetic collisionless Maxwell reservoir boundaries and finite atom pulses; no Rb optical gain or squeezing validation',
        'inputs': {'temperature_K': TEMPERATURE_K, 'mass_kg': MASS_85RB, 'density_m3': DENSITY_M3,
            'input_status': 'explicit assumed thermal gas and geometries, not independently measured data'},
        'quadrature': {'powers_per_face': list(POWERS), 'independent_scramble_seeds': list(SEEDS),
            'sobol_dimensions': 5, 'sobol_bits': 30, 'open_cube_midpoint_shift': 2.**-31,
            'face_order': ['x-', 'x+', 'y-', 'y+', 'z-', 'z+'],
            'occupancy_renormalized': False, 'fitted_arrival_rate': False},
        'frequency_hz': (OMEGA/(2*np.pi)).tolist(), 'tolerances': TOLERANCES,
        'geometry_convergence': geometry, 'same_volume_shape_control': area,
        'volume_velocity_negative_controls': negatives, 'common_phase_poisson_control': phase,
        'actual_transport_identity_control': transport,
        'source_hashes_sha256': before, 'source_hashes_stable': before == after,
        'expected_controls_passed': bool(before == after and all(row['passed'] for row in geometry.values())
            and all(row['rejected'] for row in negatives.values()) and area['passed'] and phase['passed'] and transport['passed']),
        'environment': {'python': platform.python_version(), 'numpy': np.__version__, 'scipy': scipy.__version__,
            'elapsed_s': time.monotonic()-started, 'blas_threads': 1},
        'limitations': ['Moment and identity-pulse convergence do not certify a nonlinear Rb wavepacket integral.',
            'Independent Sobol scrambles are numerical controls, not experimental uncertainty or a certified error bound.',
            'The RF phase is one common lab-entry phase; spatial offsets are fixed by each entry position.',
            'For lab-periodic drive the result is the zero-cyclic time-averaged connected stream PSD; coherent periodic mean lines and other cyclic sectors are excluded.',
            'Box inflow assumes an independent equilibrium reservoir, no wall return memory or atom-atom correlation.',
            'Collision/depletion, beam profiles, self-consistent nonlocal Maxwell coupling, optical SQL and held-out data remain downstream.']}


def save_plot(report, path):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 3, figsize=(13, 4.2), layout='constrained')
    for name, row in report['geometry_convergence'].items():
        x = np.array(POWERS)
        axes[0].loglog(6*2**x, [row['max_moment_error_by_power'][str(p)] for p in POWERS], 'o-', label=name)
        values = row['rows'][-1]
        axes[1].plot(report['frequency_hz'], values['identity_number_spectrum_s'], 'o-', label=name)
    axes[0].axhline(TOLERANCES['final_phase_space_moment_absolute'], color='grey', linestyle='--', linewidth=1)
    axes[0].set(xlabel='Total boundary paths', ylabel='Largest dimensionless moment error', title='Independent phase-space controls')
    axes[0].legend()
    axes[1].set(xlabel='Fourier frequency (Hz)', ylabel='Atomic identity pulse spectrum (s)', title='Poisson number term at equal nV', yscale='log')
    axes[1].legend()
    phase = report['common_phase_poisson_control']
    axes[2].bar(['Average mean\nbefore outer', 'Independent\ntwo phases', 'Drop Poisson\nmean term'],
        [phase['wrong_average_mean_before_outer_relative_error'], phase['wrong_independent_two_phase_relative_error'],
            phase['wrong_drop_poisson_mean_relative_error']], color=['#cd5c5c', '#e79a41', '#bd699c'])
    axes[2].set(ylabel='Relative matrix error', title='Deliberately wrong phase/noise controls')
    for ax in axes:
        ax.grid(alpha=.2, axis='y')
    fig.suptitle('Thermal boundary flux and common entry phase | Atomic foundation; no optical squeezing claim', fontsize=12)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('xb') as stream:
        fig.savefig(stream, format='png', dpi=160)
    plt.close(fig)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--plot', type=Path)
    args = parser.parse_args(argv)
    if args.plot is not None and args.output.resolve() == args.plot.resolve():
        raise ValueError('report and plot require different paths')
    for path in (args.output, args.plot):
        if path is not None and path.exists():
            raise FileExistsError(path)
    report = build_report()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open('x', encoding='utf-8') as stream:
        stream.write(json.dumps(report, ensure_ascii=False, indent=2, allow_nan=False)+'\n')
    if args.plot is not None:
        save_plot(report, args.plot)
    print(f'Wrote {args.output}; controls={report["expected_controls_passed"]}', flush=True)
    return 0 if report['expected_controls_passed'] else 1


if __name__ == '__main__':
    raise SystemExit(main())
