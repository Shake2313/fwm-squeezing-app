"""Audit nonclosed carrier phases without promoting them to a field spectrum.

python -m analysis.grand_challenge.spatial_audit --output NEW.json --plot NEW.png
"""

import argparse
from dataclasses import asdict
import hashlib
import json
from pathlib import Path
import platform
import time

import numpy as np
import scipy

from gabes import core, constants as c
from gabes.fwm_quantum.kinetic import CarrierGeometry, _pump_system
from gabes.fwm_quantum.normalization import optical_carriers
from gabes.fwm_quantum.spatial import spatial_velocity_atom
from gabes.quantum.contracts import GeneratorFrequencyAxis
from gabes.quantum.reservoirs import ExplicitReservoirs, thermal_reset_channels
from gabes.quantum.spatial import torus_atomic_noise
from .normalization_audit import conditional_inputs
from .reference.torus_qrt import trajectory_qrt, static_grating_qrt


ROOT = Path(__file__).resolve().parents[2]
OUTPUT_LABELS = [(0, 0), (1, 0), (-1, 1)]
TOLERANCES = {'toy_state_absolute': 1e-8, 'reference_spectrum_relative': 1e-8,
    'mean_refinement_absolute': 1e-8, 'response_refinement_relative': 1e-8,
    'rb_state_reference_absolute': 1e-9}


def complex_array(a):
    a = np.asarray(a)
    return {'real': a.real.tolist(), 'imag': a.imag.tolist()}


def relative(a, b):
    return float(np.linalg.norm(a-b)/max(np.linalg.norm(b), np.finfo(float).tiny))


def noise_comparison(a, b):
    return {key: relative(a[key], b[key]) for key in ('greater', 'lesser')}


def synthetic_controls():
    h = np.array([[0., .2], [.2, .7]])
    labels = [(1, 0), (-1, 1)]
    drives = np.array([[[.08, .09j], [.02, 0]], [[.03, .06], [.04j, -.02]]])
    reservoirs = ExplicitReservoirs(2, thermal_reset_channels(1., [.75, .25], source='assumed contracting toy fixture'))
    frequencies = np.array([2.3, np.sqrt(.5)])
    atom = torus_atomic_noise(h, labels, drives, reservoirs, frequencies, mean_orders=(5, 5))
    refined = torus_atomic_noise(h, labels, drives, reservoirs, frequencies, mean_orders=(6, 6))
    # Independent full-density assembly; reference consumes neither atomic A/D
    # nor the mean-state coefficients of the implementation under test.
    ref_labels = [(0, 0), *labels, *[tuple(-np.array(q)) for q in labels]]
    generators = [reservoirs.generator(h), *[core.comm_super(v) for v in drives],
        *[core.comm_super(v.conj().T) for v in drives]]
    references = []
    for phases, tail, rtol in (((8, 8), 28., 2e-10), ((12, 12), 28., 2e-10),
                              ((12, 12), 36., 2e-10), ((12, 12), 36., 2e-12)):
        started = time.monotonic()
        ref = trajectory_qrt(generators, ref_labels, frequencies, atom.operators, [.23], OUTPUT_LABELS,
            phase_samples=phases, history=tail, delay=tail, rtol=rtol, atol=rtol/100)
        references.append(ref)
        print(f'Toy direct QRT phases={phases}, history/delay={tail}: {time.monotonic()-started:.2f}s', flush=True)
    axis = GeneratorFrequencyAxis([.23], 'unfolded auxiliary phase-kernel base frequency')
    spectra = [atom.spectrum(axis, response_orders=o, output_labels=OUTPUT_LABELS) for o in ((3, 3), (4, 3), (4, 4))]
    phase = references[-1]['phases_rad']
    state_error = float(np.max(abs(atom.at_phase(phase)['state']-references[-1]['phase_states'])))
    mean_change = float(np.max(abs(atom.at_phase(phase)['state']-refined.at_phase(phase)['state'])))
    spectrum_error = relative(spectra[-1]['greater'], references[-1]['greater'])
    refinement = [noise_comparison(a, b) for a, b in zip(spectra[:-1], spectra[1:])]
    reference_changes = [relative(a['greater'], b['greater']) for a, b in zip(references[:-1], references[1:])]
    rejected = False
    try:
        torus_atomic_noise(h, labels, drives, reservoirs, frequencies, mean_orders=(1, 1))
    except ValueError:
        rejected = True
    return {'hamiltonian_rad_s': complex_array(h), 'drive_labels': labels, 'drives_rad_s': complex_array(drives),
        'frequencies_rad_s': frequencies.tolist(), 'reset_rate_s_inverse': 1., 'reset_populations': [.75, .25],
        'mean_orders': [5, 5], 'refined_mean_orders': [6, 6],
        'mean_state_reference_maximum_absolute_error': state_error, 'mean_refinement_maximum_absolute_change': mean_change,
        'spectrum_reference_relative_error': spectrum_error,
        'response_orders': [[3, 3], [4, 3], [4, 4]], 'adjacent_response_changes': refinement,
        'reference_controls': [{key: ref[key] for key in ('phase_samples', 'history', 'delay', 'rtol', 'atol',
            'maximum_regression_endpoint_norm')} for ref in references],
        'reference_phase_tail_tolerance_relative_changes': reference_changes,
        'insufficient_mean_order_rejected': rejected, 'diagnostics': atom.diagnostics,
        'passed': bool(state_error < TOLERANCES['toy_state_absolute'] and mean_change < TOLERANCES['mean_refinement_absolute']
            and spectrum_error < TOLERANCES['reference_spectrum_relative']
            and max(reference_changes) < TOLERANCES['reference_spectrum_relative']
            and max(v for row in refinement for v in row.values()) < TOLERANCES['response_refinement_relative'] and rejected)}


def build_report():
    with core.blas_single_thread():
        return _build_report()


def _build_report():
    toy = synthetic_controls()
    inputs = conditional_inputs()
    geometry = CarrierGeometry.vacuum_beams(inputs, probe_angle_rad=.005, conjugate_angle_rad=-.005)
    amplitude = np.sqrt(inputs.seed_power_W/(c.HBAR*optical_carriers(inputs.detunings)[1]))
    beta = np.array([amplitude, .25j*amplitude])
    axis = GeneratorFrequencyAxis([0., 2*np.pi*1e6], 'atomic phase-kernel offsets; no optical-mode/SQL projection')
    cases = []
    static_selected = None
    psi = 2*np.pi*np.arange(65)/64
    for velocity in ([0., 0., 0.], [150., 0., 380.], [0., 0., 700.]):
        started = time.monotonic()
        model = spatial_velocity_atom(inputs, geometry, velocity, beta, mean_orders=(3, 3))
        refined = spatial_velocity_atom(inputs, geometry, velocity, beta, mean_orders=(4, 4))
        phases = np.column_stack([np.full(len(psi), .17), psi])
        mean_change = float(np.max(abs(model.atom.at_phase(phases)['state']-refined.atom.at_phase(phases)['state'])))
        spectra = [model.atom.spectrum(axis, response_orders=o, output_labels=OUTPUT_LABELS)
                   for o in ((1, 1), (2, 2), (3, 2), (3, 3))]
        comparisons = [noise_comparison(a, b) for a, b in zip(spectra[:-1], spectra[1:])]
        # Negative control changes the convective equation and recomputes mean,
        # A and D together. It is never used to repair the physical geometry.
        h0, reservoirs = _pump_system(inputs, inputs.one_photon_rad_s-geometry.wavevectors_rad_m[0]@velocity,
                                      'uniform-zeeman-rms')
        erased = torus_atomic_noise(h0, [(1, 0), (-1, 1)], model.drive_operators_rad_s, reservoirs,
            [model.convective_frequencies_rad_s[0], 0.], mean_orders=(3, 3))
        erased_spectrum = erased.spectrum(axis, response_orders=(3, 3), output_labels=OUTPUT_LABELS)
        row = {'velocity_m_s': velocity, 'coordinate_kind': model.coordinate_kind,
            'convective_frequencies_rad_s': model.convective_frequencies_rad_s.tolist(),
            'optical_drive_frequencies_rad_s': [model.convective_frequencies_rad_s[0],
                -model.convective_frequencies_rad_s[0]+model.convective_frequencies_rad_s[1]],
            'mean_orders': [3, 3], 'refined_mean_orders': [4, 4], 'mean_refinement_maximum_absolute_change': mean_change,
            'response_orders': [[1, 1], [2, 2], [3, 2], [3, 3]], 'adjacent_response_changes': comparisons,
            'mean_diagnostics': model.atom.diagnostics, 'spectrum_audit': spectra[-1]['audit'],
            'zero_loop_convection_negative_control': {'spectrum_relative_change': noise_comparison(erased_spectrum, spectra[-1]),
                'state_maximum_absolute_change': float(np.max(abs(erased.at_phase(phases)['state']-model.atom.at_phase(phases)['state']))),
                'scope': 'different equation, retained for sensitivity; not an approximation certified for readout'},
            'state_labels': model.atom.state_labels.tolist(), 'state_coefficients': complex_array(model.atom.state_coefficients),
            'selected_greater_s': complex_array(spectra[-1]['greater']), 'selected_lesser_s': complex_array(spectra[-1]['lesser']),
            'grating_phase_rad': psi.tolist(), 'g2_population_at_theta1_0p17': model.atom.at_phase(phases)['state'][:, 1, 1].real.tolist(),
            'seconds': time.monotonic()-started,
            'passed': bool(mean_change < TOLERANCES['mean_refinement_absolute']
                and max(v for r in comparisons[1:] for v in r.values()) < TOLERANCES['response_refinement_relative'])}
        cases.append(row)
        if static_selected is None:
            static_selected = (model, h0, reservoirs, spectra[-1])
        print(f'Rb spatial velocity={velocity}: {row["seconds"]:.2f}s', flush=True)
    model, h0, reservoirs, static_spectrum = static_selected
    pairs = [(core.comm_super(v), core.comm_super(v.conj().T)) for v in model.drive_operators_rad_s]
    references = [static_grating_qrt(reservoirs.generator(h0), *pairs, model.convective_frequencies_rad_s[0],
        model.atom.operators, axis.omega_rad_s, OUTPUT_LABELS, grating_samples=ng, phase_samples=nt)
        for ng, nt in ((4, 16), (8, 16), (8, 24))]
    reference = references[-1]
    rb_error = relative(static_spectrum['greater'], reference['greater'])
    rb_state_error = float(np.max(abs(model.atom.at_phase(reference['phases_rad'])['state']-reference['phase_states'])))
    ref_changes = [relative(a['greater'], b['greater']) for a, b in zip(references[:-1], references[1:])]
    rb_ref = {'spectrum_relative_error': rb_error, 'state_maximum_absolute_error': rb_state_error,
        'phase_grid_refinements': [[4, 16], [8, 16], [8, 24]], 'adjacent_reference_relative_changes': ref_changes,
        'scope': reference['scope'], 'passed': bool(rb_error < TOLERANCES['reference_spectrum_relative']
            and rb_state_error < TOLERANCES['rb_state_reference_absolute']
            and max(ref_changes) < TOLERANCES['reference_spectrum_relative'])}
    closed = spatial_velocity_atom(inputs, CarrierGeometry.vacuum_beams(inputs), [150, 0, 380], beta)
    previous = ROOT/'docs/grand_challenge/s1_kinetic_report_v2.json'
    sources = [ROOT/'gabes/core.py', ROOT/'gabes/doppler.py', ROOT/'gabes/observables.py', ROOT/'gabes/constants.py',
        ROOT/'gabes/atoms.py', ROOT/'gabes/hyperfine.py', ROOT/'gabes/schemes/fwm.py', ROOT/'gabes/species.py',
        *sorted((ROOT/'gabes/quantum').glob('*.py')), *sorted((ROOT/'gabes/fwm_quantum').glob('*.py')),
        Path(__file__), ROOT/'analysis/grand_challenge/normalization_audit.py',
        ROOT/'analysis/grand_challenge/reference/torus_qrt.py', ROOT/'analysis/grand_challenge/reference/periodic_qrt.py',
        ROOT/'tests/quantum/test_spatial.py']
    return {'stage': 'S1 general two-phase local convective atomic mean and ordered noise',
        'expected_controls_passed': bool(toy['passed'] and all(r['passed'] for r in cases) and rb_ref['passed']
            and not geometry.closure_audit()['passed'] and closed.coordinate_kind == 'closed-single-phase-quotient'),
        'numerical_tolerances': TOLERANCES, 'synthetic_independent_controls': toy,
        'input_ledger': asdict(inputs), 'input_provenance': 'assumed reduced RMS-dipole fixture, no gain/noise fit or independent-input promotion',
        'geometry': {'wavevectors_rad_m': geometry.wavevectors_rad_m.tolist(), 'source': geometry.source,
            'single_phase_audit': geometry.closure_audit(temperature_K=394.15, length_m=inputs.length_m),
            'phase_wavevectors_rad_m': model.wavevectors_rad_m.tolist(), 'phase_labels': ['probe beat', 'loop grating']},
        'prescribed_carriers_sqrt_photons_per_second': complex_array(beta),
        'carrier_scope': 'local prescribed probe and conjugate; conjugate=.25i*probe, not a self-consistent exit solution',
        'signed_base_angular_frequencies_rad_s': axis.omega_rad_s.tolist(), 'output_atomic_labels': OUTPUT_LABELS,
        'rb_cases': cases, 'rb_independent_static_grating_qrt': rb_ref,
        'closed_geometry_coordinate_kind': closed.coordinate_kind,
        'historical_predecessor': {'path': str(previous.relative_to(ROOT)).replace('\\', '/'),
            'sha256': hashlib.sha256(previous.read_bytes()).hexdigest(),
            'scope': 'historical context only; no earlier numerical controls reused or current-source parity asserted'},
        'source_sha256': {str(p.resolve().relative_to(ROOT)).replace('\\', '/'): hashlib.sha256(p.read_bytes()).hexdigest() for p in sources},
        'environment': {'python': platform.python_version(), 'numpy': np.__version__, 'scipy': scipy.__version__, 'blas_threads': 1},
        'limits': ['Local independent ballistic atom with fixed plane-wave envelopes and pump; no spatial-envelope transport or class-changing collisions.',
            'Two-phase atomic coordinates are not independent optical modes. No general-geometry canonical field M/D, gain or S_minus is reported.',
            'The one-velocity direct moving-trajectory reference uses a synthetic two-level atom; the four-level Rb time-QRT comparison fixes v.Q=0 and retains spatial phase.',
            'Finite rectangular harmonic cutoffs and sampled state positivity are numerical controls, not a rigorous uniform error bound.',
            'Full Zeeman/polarization, measured collection, pump depletion, higher-order currents and held-out no-fit experimental validation remain open.']}


def save_plot(report, path):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2), constrained_layout=True)
    colors = ['#355c9a', '#c55a31', '#438565']
    for row, color in zip(report['rb_cases'], colors):
        population = np.asarray(row['g2_population_at_theta1_0p17'])
        label = str(tuple(int(v) for v in row['velocity_m_s']))
        axes[0].plot(np.asarray(row['grating_phase_rad'])/np.pi, 1e6*(population-population.mean()), color=color, label=label)
    axes[0].set(xlabel='Loop grating phase / pi', ylabel='g2 population minus phase mean (ppm)', title='Local spatial modulation, fixed probe phase')
    axes[0].legend(title='Velocity (vx, vy, vz), m/s', fontsize=8)
    changes = [100*r['zero_loop_convection_negative_control']['spectrum_relative_change']['greater'] for r in report['rb_cases']]
    axes[1].bar(range(3), changes, color=colors)
    axes[1].set_xticks(range(3), ['(0, 0, 0)', '(150, 0, 380)', '(0, 0, 700)'])
    axes[1].set(xlabel='Velocity (vx, vy, vz), m/s', ylabel='Atomic spectrum relative change (%)', title='Removing loop convection changes the equation')
    for i, value in enumerate(changes):
        axes[1].annotate(f'{value:.3g}%', (i, value), xytext=(0, 5), textcoords='offset points', ha='center')
    axes[1].margins(y=.18)
    for ax in axes:
        ax.grid(alpha=.2, axis='y'); ax.set_axisbelow(True)
    fig.suptitle('Nonclosed vacuum geometry: local reduced atom, no field squeezing prediction', fontsize=12)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('xb') as stream:
        fig.savefig(stream, format='png', dpi=160)
    plt.close(fig)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--plot', type=Path)
    args = parser.parse_args(argv)
    for path in (args.output, args.plot):
        if path is not None and path.exists():
            raise FileExistsError(path)
    report = build_report()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open('x', encoding='utf-8') as stream:
        stream.write(json.dumps(report, ensure_ascii=False, indent=2, allow_nan=False)+'\n')
    print(f'Wrote {args.output}', flush=True)
    if args.plot is not None:
        save_plot(report, args.plot)
    return 0 if report['expected_controls_passed'] else 1


if __name__ == '__main__':
    raise SystemExit(main())
