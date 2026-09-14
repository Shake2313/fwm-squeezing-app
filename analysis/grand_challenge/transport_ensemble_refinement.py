"""Refine the existing constant-atom boundary ensemble without changing physics.

python -m analysis.grand_challenge.transport_ensemble_refinement --output NEW.json --plot NEW.png

Exact fixture factorization preserves every physical path and its arrival rate.
It is not a new physical-Rb solver or a production callback certification.
"""

import argparse
import hashlib
import itertools
import json
from pathlib import Path
import platform
import sys
import time

import numpy as np
import scipy

from gabes import core
from . import transport_ensemble_audit as a


POWERS = (10, 12, 14, 16)
SEEDS = (11, 211, 811)
CHUNK_SIZE = 4096
PATH_BUDGET = a.PATH_BUDGET
STREAM_BUDGET = a.STREAM_BUDGET
PARENT_PATH = a.ROOT/'docs/grand_challenge/transport_ensemble_report_v1.json'
PARENT_SHA256 = 'b2571c4777456b0a53ce528bfefecbc99c9dc013f475e79719f29948b9d9141c'


def source_hashes():
    """Hash imported project code and the tests of the consumed calculation.

    Concurrent, unrelated new reference files are intentionally not globbed.
    """
    paths = {Path(__file__).resolve(), Path(a.__file__).resolve()}
    for name, module in tuple(sys.modules.items()):
        if name == 'gabes' or name.startswith('gabes.'):
            filename = getattr(module, '__file__', None)
            if filename and filename.endswith('.py'):
                paths.add(Path(filename).resolve())
    paths.update(a.ROOT/name for name in (
        'tests/quantum/test_transport_ensemble_refinement.py',
        'tests/quantum/test_transport_ensemble.py', 'tests/quantum/test_inflow.py'))
    return {p.relative_to(a.ROOT).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in sorted(paths)}


def _digest(value):
    return hashlib.sha256(json.dumps(a.json_value(value), sort_keys=True,
        allow_nan=False, separators=(',', ':')).encode()).hexdigest()


def row_digest(row):
    return _digest({key: value for key, value in row.items() if key != 'content_digest'})


def batch_pulse_integrals(taus, *, order=None):
    """Evaluate the unchanged finite-pulse integrals for all supplied ages."""
    taus = np.asarray(taus, dtype=float)
    if taus.ndim != 1 or not len(taus) or not np.all(np.isfinite(taus)) or np.any(taus <= 0):
        raise ValueError('Residence times must be a nonempty positive finite vector')
    tau = taus[:, None]
    z = 1j*tau*a.AXIS.omega_rad_s[None, :]
    if order is not None:
        if not isinstance(order, (int, np.integer)) or isinstance(order, bool) or order < 1:
            raise ValueError('Gauss order must be a positive integer')
        nodes, weights = a.gauss_rule(int(order))
        oscillation = np.exp(z[:, :, None]*nodes)
        return tau*(oscillation@weights), tau**2*(oscillation@(weights*(1-nodes)))
    f = tau*np.sinc(tau*a.AXIS.omega_rad_s/(2*np.pi))*np.exp(z/2)
    response = np.empty_like(z)
    small = np.abs(z) < .1
    term = np.full(np.count_nonzero(small), .5, dtype=complex)
    series = term.copy()
    for power in range(1, 16):
        term = term*z[small]/(power+2)
        series += term
    response[small] = series
    response[~small] = (np.expm1(z[~small])-z[~small])/z[~small]**2
    return f, tau**2*response


def batch_path_errors(taus, candidate, reference, atom):
    """The original RF Frobenius errors, using ||scalar*A||=|scalar| ||A||."""
    taus = np.asarray(taus, dtype=float)
    expected = (len(taus), len(a.AXIS.omega_rad_s))
    fc, tc = (np.asarray(item) for item in candidate)
    fr, tr = (np.asarray(item) for item in reference)
    if any(item.shape != expected for item in (fc, tc, fr, tr)):
        raise ValueError('Pulse integrals must have shape (paths, laboratory RF)')
    outer = np.outer(atom['mean'], atom['mean'].conj())
    base = {'greater': atom['greater'], 'lesser': atom['lesser'],
        'greater_by_source': atom['greater'], 'lesser_by_source': atom['lesser'],
        'mean_pulse': atom['mean'], 'mean_outer': outer,
        'retarded_response': atom['response']}
    scales = a.metric_scales(taus[:, None], atom)
    result = {}
    for key in a.PATH_KEYS:
        c, r = ((fc, fr) if key == 'mean_pulse' else (tc, tr)
                if key == 'retarded_response' else (np.abs(fc)**2, np.abs(fr)**2))
        norm = np.linalg.norm(base[key])
        denominator = np.maximum(np.abs(r)*norm, 1e-12*scales[key])
        result[key] = np.max(np.abs(c-r)*norm/denominator, axis=1)
    return result


def _spectra(pulse_weight, response_weight, atom):
    # Only this fixed atom, zero offsets, common RF and uniform common mark.
    charge = np.array([1, -1, -1, 1])
    mask = charge[:, None] == charge[None, :]
    outer = mask*np.outer(atom['mean'], atom['mean'].conj())
    number = pulse_weight[:, None, None]*outer
    result = {'poisson_number': number,
        'retarded_response': response_weight[:, None, None]*(mask*atom['response'])}
    for name in ('greater', 'lesser'):
        connected = pulse_weight[:, None, None]*(mask*atom[name])
        result['internal_'+name] = connected
        result[name+'_by_source'] = connected[None]
        result[name] = connected+number
    return result


def factorized_grid(power, seed, atom=None, *, density=a.DENSITY_M3, chunk_size=CHUNK_SIZE):
    """Return actual grid spectra and measured evidence, without packet callbacks."""
    started = time.monotonic()
    if not isinstance(chunk_size, int) or isinstance(chunk_size, bool) or chunk_size < 1:
        raise ValueError('chunk_size must be a positive integer')
    atom = a.toy_atom() if atom is None else atom
    # This routine is deliberately limited to the exact original atom.
    if _digest(atom) != _digest(a.toy_atom()):
        raise ValueError('Factorization is restricted to the unchanged constant-atom fixture')
    inflow = a.make_inflow(power, seed, density=density)
    count = len(inflow.residence_time_s)
    maximum = {kind: {key: 0. for key in a.PATH_KEYS}
               for kind in ('path_refinement', 'independent_reference')}
    worst = None
    failure_count = 0
    pulse_weight = np.zeros(len(a.AXIS.omega_rad_s))
    response_weight = np.zeros(len(pulse_weight), dtype=complex)
    exact_pulse_weight = np.zeros_like(pulse_weight)
    exact_response_weight = np.zeros_like(response_weight)
    chain = hashlib.sha256()
    # Bind the actual complete boundary realization, not only its residence times.
    for array in (inflow.entry_position_m, inflow.velocity_m_s,
                  inflow.residence_time_s, inflow.rate_s_inverse, inflow.face_index):
        chain.update(np.ascontiguousarray(array).tobytes())
    for begin in range(0, count, chunk_size):
        end = min(begin+chunk_size, count)
        tau = inflow.residence_time_s[begin:end]
        rate = inflow.rate_s_inverse[begin:end]
        coarse = batch_pulse_integrals(tau, order=a.GAUSS_ORDERS[0])
        fine = batch_pulse_integrals(tau, order=a.GAUSS_ORDERS[1])
        exact = batch_pulse_integrals(tau)
        failed = np.zeros(len(tau), dtype=bool)
        for kind, reference in (('path_refinement', coarse), ('independent_reference', exact)):
            errors = batch_path_errors(tau, fine, reference, atom)
            for key, values in errors.items():
                failed |= ~np.isfinite(values) | (values > PATH_BUDGET)
                index = int(np.argmax(values))
                value = float(values[index])
                maximum[kind][key] = max(maximum[kind][key], value)
                if worst is None or value > worst['error']:
                    worst = {'kind': kind, 'metric': key, 'error': value,
                        'path_index': begin+index, 'residence_time_s': float(tau[index])}
        failure_count += int(np.count_nonzero(failed))
        for pair in (coarse, fine, exact):
            for array in pair:
                chain.update(np.ascontiguousarray(array).tobytes())
        pulse_weight += np.einsum('p,pf->f', rate, np.abs(fine[0])**2)
        response_weight += np.einsum('p,pf->f', rate, fine[1])
        exact_pulse_weight += np.einsum('p,pf->f', rate, np.abs(exact[0])**2)
        exact_response_weight += np.einsum('p,pf->f', rate, exact[1])
    spectra = _spectra(pulse_weight, response_weight, atom)
    scales = a.metric_scales(float(np.max(a.LENGTHS_M)/np.sqrt(
        a.c.KB*a.TEMPERATURE_K/a.c.MASS_85RB)), atom,
        rate=inflow.total_arrival_rate_s_inverse)
    exact_errors = a.observed_errors(spectra,
        _spectra(exact_pulse_weight, exact_response_weight, atom), a.STREAM_KEYS, scales)
    row = {'power': power, 'seed': seed, 'path_count': count,
        'evaluated_path_count': count, 'packet_callback_count': 0,
        'evaluation_method': 'Exact fixed-atom factorization; Gauss24/Gauss48/closed form evaluated on every actual path',
        'density_m3': density, 'mean_occupancy': inflow.mean_occupancy,
        'equilibrium_nV': inflow.equilibrium_atom_number,
        'occupancy_over_nV': inflow.mean_occupancy/inflow.equilibrium_atom_number,
        'total_arrival_rate_s_inverse': inflow.total_arrival_rate_s_inverse,
        'residence_time_range_s': [float(inflow.residence_time_s.min()), float(inflow.residence_time_s.max())],
        'path_evidence_passed': failure_count == 0, 'failed_path_count': failure_count,
        'max_path_errors': maximum, 'worst_path_error': worst,
        'all_boundary_and_integral_arrays_sha256': chain.hexdigest(),
        'closed_form_stream_errors': exact_errors, 'comparison_scales': scales,
        'pulse_weight_s': pulse_weight, 'response_weight_s': response_weight,
        'spectra': spectra, 'elapsed_s': time.monotonic()-started}
    row['content_digest'] = row_digest(row)
    return spectra, row


def convergence_gate(rows, *, powers=POWERS, seeds=SEEDS):
    """Require all seeds, two adjacent refinements, and all directed seed pairs.

    Values are remeasured from the actual matrices. Occupancy is never a gate.
    """
    reasons, refinements, scrambles = [], [], []
    powers, seeds = tuple(powers), tuple(seeds)
    if len(powers) < 3 or len(set(powers)) != len(powers) or tuple(sorted(powers)) != powers:
        return {'passed': False, 'reasons': ['Plan requires at least three increasing unique powers'],
                'refinements': [], 'scrambles': []}
    if len(seeds) < 3 or len(set(seeds)) != len(seeds):
        return {'passed': False, 'reasons': ['Plan requires at least three distinct seeds'],
                'refinements': [], 'scrambles': []}
    indexed = {}
    densities = set()
    try:
        for row in rows:
            key = (row['power'], row['seed'])
            if key in indexed or key[0] not in powers or key[1] not in seeds:
                reasons.append('Duplicate or unplanned power/seed')
            indexed[key] = row
            density = float(row['density_m3'])
            if not np.isfinite(density) or density <= 0:
                reasons.append('A positive finite physical density is required')
            densities.add(density)
            if row['content_digest'] != row_digest(row):
                reasons.append('Grid content digest mismatch')
            if not row['path_evidence_passed'] or row['failed_path_count'] != 0:
                reasons.append('A path failed the unchanged numerical evidence budget')
            if row['path_count'] != 6*2**row['power'] or row['evaluated_path_count'] != row['path_count']:
                reasons.append('Incomplete physical boundary path count')
            for kind in ('path_refinement', 'independent_reference'):
                values = [row['max_path_errors'][kind][key] for key in a.PATH_KEYS]
                if not np.all(np.isfinite(values)) or min(values) < 0 or max(values) > PATH_BUDGET:
                    reasons.append('Incomplete or failed seven-metric path evidence')
            stream_errors = [row['closed_form_stream_errors'][key] for key in a.STREAM_KEYS]
            if (not np.all(np.isfinite(stream_errors)) or min(stream_errors) < 0
                    or max(stream_errors) > a.PHASE_BUDGET):
                reasons.append('Closed-form stream comparison exceeds original 2e-12 budget')
        present = sorted({p for p, seed in indexed})
        if len(present) < 3 or tuple(present) != powers[:len(present)]:
            reasons.append('At least three consecutive planned grids required')
        if set(indexed) != set(itertools.product(present, seeds)):
            reasons.append('Every planned seed required at every evaluated grid')
        if len(densities) != 1:
            reasons.append('All refinements and scrambles must use the same physical density')
        if reasons:
            return {'passed': False, 'reasons': sorted(set(reasons)), 'refinements': [], 'scrambles': []}

        def compare(kind, candidate, reference):
            errors = a.observed_errors(candidate['spectra'], reference['spectra'],
                a.STREAM_KEYS, candidate['comparison_scales'])
            passed = bool(np.all(np.isfinite(list(errors.values()))) and max(errors.values()) <= STREAM_BUDGET)
            return {'kind': kind, 'candidate_power': candidate['power'], 'candidate_seed': candidate['seed'],
                'reference_power': reference['power'], 'reference_seed': reference['seed'],
                'candidate_digest': candidate['content_digest'], 'reference_digest': reference['content_digest'],
                'errors': errors, 'budget': STREAM_BUDGET, 'passed': passed}

        for previous, power in zip(present[-3:-1], present[-2:]):
            for seed in seeds:
                refinements.append(compare('ensemble_refinement', indexed[power, seed], indexed[previous, seed]))
        for power in present[-2:]:
            for candidate_seed, reference_seed in itertools.permutations(seeds, 2):
                scrambles.append(compare('independent_scramble', indexed[power, candidate_seed], indexed[power, reference_seed]))
    except (KeyError, ValueError, TypeError, IndexError) as error:
        reasons.append('Incomplete or invalid measured grid data: '+str(error))
    if not all(row['passed'] for row in refinements+scrambles):
        reasons.append('Actual spectrum/source/number/response convergence exceeds 5%')
    return {'passed': not reasons and len(refinements) == 2*len(seeds)
            and len(scrambles) == 2*len(seeds)*(len(seeds)-1),
        'reasons': sorted(set(reasons)), 'refinements': refinements, 'scrambles': scrambles}


def _from_json_arrays(value):
    if isinstance(value, dict) and set(value) == {'real', 'imag'}:
        return np.asarray(value['real'])+1j*np.asarray(value['imag'])
    if isinstance(value, dict):
        return {key: _from_json_arrays(item) for key, item in value.items()}
    return value


def adapter_parity(atom):
    rows = []
    for power, seed in ((0, 11), (2, 211)):
        candidate, original, independent_phase = a.run_grid(power, seed, atom)
        spectra, factorized = factorized_grid(power, seed, atom)
        errors = a.observed_errors(spectra, candidate['spectra'], a.STREAM_KEYS, original['comparison_scales'])
        phase_errors = a.observed_errors(spectra, independent_phase, a.STREAM_KEYS, original['comparison_scales'])
        rows.append({'power': power, 'seed': seed, 'path_count': original['path_count'],
            'actual_callback_count': original['callback_count'], 'original_candidate_digest': original['candidate_digest'],
            'factorized_digest': factorized['content_digest'], 'stream_errors': errors,
            'independent_literal_phase_errors': phase_errors,
            'passed': bool(candidate['path_evidence_passed'] and factorized['path_evidence_passed']
                and max(errors.values()) <= a.PHASE_BUDGET and max(phase_errors.values()) <= a.PHASE_BUDGET)})
    return {'grids': rows, 'budget': a.PHASE_BUDGET, 'passed': all(row['passed'] for row in rows)}


def build_report():
    started = time.monotonic()
    raw_parent = PARENT_PATH.read_bytes()
    if hashlib.sha256(raw_parent).hexdigest() != PARENT_SHA256:
        raise ValueError('Historical parent report hash changed')
    parent = json.loads(raw_parent)
    atom = a.toy_atom()
    rows, history = [], []
    with core.blas_single_thread():
        parity = adapter_parity(atom)
        constant = a.actual_constant_transport_control(atom)
        for power in POWERS:
            for seed in SEEDS:
                _, row = factorized_grid(power, seed, atom)
                rows.append(row)
                print(f'p{power} seed{seed}: {row["path_count"]} paths; path evidence={row["path_evidence_passed"]}; '
                      f'nV ratio={row["occupancy_over_nV"]:.9f}', flush=True)
            gate = convergence_gate(rows)
            comparisons = gate['refinements']+gate['scrambles']
            maximum = max((max(row['errors'].values()) for row in comparisons), default=None)
            history.append({'power': power, 'gate': gate, 'maximum_required_error': maximum})
            print(f'p{power}: two-refinement / all-seed gate={gate["passed"]}; maximum={maximum}', flush=True)
            if gate['passed'] or not all(row['path_evidence_passed'] for row in rows):
                break
    nominal = next(row for row in rows if (row['power'], row['seed']) == (10, 11))
    parent_errors = a.observed_errors(nominal['spectra'],
        _from_json_arrays(parent['final_candidate']['spectra']), a.STREAM_KEYS, nominal['comparison_scales'])
    parent_passed = max(parent_errors.values()) <= a.PHASE_BUDGET
    final = next(row for row in rows if row['power'] == history[-1]['power'] and row['seed'] == SEEDS[0])
    closed_form_passed = all(np.all(np.isfinite(list(row['closed_form_stream_errors'].values())))
        and max(row['closed_form_stream_errors'].values()) <= a.PHASE_BUDGET for row in rows)
    implementation_passed = bool(parity['passed'] and constant['passed'] and parent_passed
        and closed_form_passed and all(row['path_evidence_passed'] for row in rows))
    return {'schema': 'gabes-constant-atom-ensemble-refinement-v1',
        'scope': 'Same constant two-level phase-mark fixture as parent, evaluated by exact factorization; no physical-Rb or optical certification',
        'historical_parent': {'path': PARENT_PATH.relative_to(a.ROOT).as_posix(), 'sha256': PARENT_SHA256,
            'certified': parent['final_candidate']['certified'], 'p10_reproduction_errors': parent_errors,
            'reproduction_budget': a.PHASE_BUDGET, 'reproduction_passed': parent_passed},
        'toy_model': parent['toy_model'], 'assumed_boundary_inputs': parent['assumed_boundary_inputs'],
        'frequency_hz': a.AXIS.frequency_hz, 'declared_plan': {
            'powers_per_face': POWERS, 'seeds': SEEDS, 'gauss_orders': a.GAUSS_ORDERS,
            'chunk_size': CHUNK_SIZE, 'path_budget': PATH_BUDGET, 'stream_budget': STREAM_BUDGET,
            'closed_form_stream_budget': a.PHASE_BUDGET,
            'required_refinement_edges_per_seed': 2, 'required_scramble_grids': 2,
            'scramble_comparisons': 'All directed pairs of three predeclared independent seeds',
            'stopping_rule': 'First passing complete gate from p14; otherwise evaluate through p16; any path failure prevents certification',
            'error_definition': a.ERROR_DEFINITION,
            'density_count': 'Once in physical boundary arrival rates; every path retained; no occupancy renormalization',
            'factorization': 'A=sum(rate*abs(F48)^2), B=sum(rate*T48); full raw/source/mean-outer/complex response matrices reconstructed exactly for this fixed atom'},
        'grids': rows, 'convergence_history': history, 'final_convergence': history[-1]['gate'],
        'final_result': {'power': final['power'], 'nominal_seed': SEEDS[0], 'path_count_per_seed': final['path_count'],
            'spectra': final['spectra'], 'numerical_units': parent['final_candidate']['numerical_units'],
            'constant_atom_ensemble_converged': bool(implementation_passed and history[-1]['gate']['passed']),
            'production_adapter_certified': False,
            'production_adapter_reason': 'High-resolution calculation uses exact fixture factorization; no per-path production packet-digest chain was created',
            'physical_Rb_ensemble_certified': False},
        'same_grid_production_adapter_parity': parity, 'actual_constant_transport_control': constant,
        'all_grid_closed_form_stream_control_passed': bool(closed_form_passed),
        'implementation_controls_passed': implementation_passed,
        'expected_controls_passed': bool(implementation_passed and history[-1]['gate']['passed']),
        'environment': {'python': platform.python_version(), 'numpy': np.__version__, 'scipy': scipy.__version__,
            'blas_threads': 1, 'elapsed_s': time.monotonic()-started},
        'limitations': [
            'Exact factorization requires the fixed state/readouts/drives, H=0, no jumps and zero port offsets; it cannot replace physical Rb path solves.',
            'Only atomic_inflow exists here. Nonzero microscopic jump-source thermal convergence remains outside this fixture.',
            'Two refinements and independent scrambles are observed numerical controls, not rigorous statistical confidence intervals.',
            'Retarded response remains complex and age-integrated; it is not a spatial Maxwell matrix or a field squeezing prediction.',
            'The historical failed report remains unchanged; its production candidate has not been retroactively certified.']}


def save_plot(report, path):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.2), layout='constrained')
    indexed = {(row['power'], row['seed']): row for row in report['grids']}
    powers = sorted({power for power, seed in indexed})
    for seed in SEEDS:
        errors = [max(a.observed_errors(indexed[p, seed]['spectra'], indexed[previous, seed]['spectra'],
            a.STREAM_KEYS, indexed[p, seed]['comparison_scales']).values())
            for previous, p in zip(powers[:-1], powers[1:])]
        axes[0].plot(powers[1:], 100*np.asarray(errors), 'o-', label=f'Seed {seed}')
    axes[0].axhline(100*STREAM_BUDGET, color='#a5443b', ls='--', label='Unchanged 5% budget')
    axes[0].set(xlabel='Sobol power per face', ylabel='Maximum matrix error (%)', title='Every seed: adjacent grid refinement')
    axes[0].legend(fontsize=8)
    spread = []
    for power in powers:
        spread.append(max(max(a.observed_errors(indexed[power, i]['spectra'], indexed[power, j]['spectra'],
            a.STREAM_KEYS, indexed[power, i]['comparison_scales']).values()) for i, j in itertools.permutations(SEEDS, 2)))
    axes[1].plot(powers, 100*np.asarray(spread), 'o-', color='#247c77')
    axes[1].axhline(100*STREAM_BUDGET, color='#a5443b', ls='--')
    axes[1].set(xlabel='Sobol power per face', ylabel='Maximum matrix error (%)', title='All directed independent seed pairs')
    final = report['final_result']['spectra']
    for key, label in (('greater', 'Raw greater'), ('internal_greater', 'Connected atom'), ('poisson_number', 'Poisson mean outer')):
        axes[2].plot(a.AXIS.frequency_hz/1e3, np.asarray(final[key]).real[:, 0, 0], 'o-', label=label)
    axes[2].set(xlabel='Laboratory RF (kHz)', ylabel='Toy atomic spectrum (s)', title='Mean outer retained; density once')
    axes[2].legend(fontsize=8)
    for ax in axes:
        ax.grid(alpha=.2)
    fig.suptitle('Constant-atom boundary refinement | Same physics and 5% budget | No physical-Rb certification', fontsize=11)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('xb') as output:
        fig.savefig(output, format='png', dpi=160)
    plt.close(fig)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--plot', type=Path)
    args = parser.parse_args(argv)
    if args.plot is not None and args.output.resolve() == args.plot.resolve():
        raise ValueError('Report and plot require different paths')
    for path in (args.output, args.plot):
        if path is not None and path.exists():
            raise FileExistsError(path)
    before = source_hashes()
    report = build_report()
    after = source_hashes()
    if before != after:
        raise RuntimeError('Source files changed during the audit; no immutable artifact written')
    report['source_hashes_sha256'] = before
    report['source_hashes_stable'] = True
    payload = json.dumps(a.json_value(report), ensure_ascii=False, indent=2, allow_nan=False)+'\n'
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open('x', encoding='utf-8') as output:
        output.write(payload)
    if args.plot is not None:
        save_plot(report, args.plot)
    print(f'Wrote {args.output}; controls={report["expected_controls_passed"]}', flush=True)
    return 0 if report['expected_controls_passed'] else 1


if __name__ == '__main__':
    raise SystemExit(main())
