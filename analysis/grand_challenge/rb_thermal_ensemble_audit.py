"""Conditional thermal-Rb workload and six-boundary-path evidence pilot.

Uncomputed grid/seed combinations remain explicit. A complete pilot grid is
not an ensemble convergence certificate. Immutable references are verified
before any reuse; the prescribed selected 2 us path is never substituted.
"""

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
import hashlib
import json
import math
import time

import numpy as np
import threadpoolctl

from gabes import core
from . import rb_thermal_ensemble as r
from . import rb_thermal_reference_jobs as refs


WORKLOAD_POWERS = (0, 1, 2, 3)
SEEDS = (11, 211, 811)
PRIMARY_METHOD = 'CF4-exact-spectral-microscopic-source-propagation-v1'
REFERENCE_METHOD = 'independent-backward-observable-Lindblad-products-v1'
DEFAULT_SEGMENTS = (256, 512, 1024)
DEFAULT_MAX_STEPS_S = (2e-6/4096, 2e-6/8192, 2e-6/16384)


def plan_for(segments=DEFAULT_SEGMENTS):
    segments = tuple(segments)
    if (len(segments) != 3 or any(isinstance(n, (bool, np.bool_))
            or not isinstance(n, (int, np.integer)) or n <= 0 for n in segments)
            or any(a >= b for a, b in zip(segments[:-1], segments[1:]))):
        raise ValueError('Three strictly increasing positive integer segment counts required')
    return r.PathPlan(
        tuple(r.SolverSpec(PRIMARY_METHOD, int(n), {'segments': int(n), 'order': 4}) for n in segments),
        tuple(r.SolverSpec(REFERENCE_METHOD, 1/rtol, {'rtol': rtol, 'atol': atol, 'max_step_fraction': 1/64})
              for rtol, atol in refs.REFERENCE_TOLERANCES))


def step_plan(max_steps_s=DEFAULT_MAX_STEPS_S):
    steps = tuple(max_steps_s)
    if (len(steps) != 3 or any(isinstance(value, (bool, np.bool_))
            or not isinstance(value, (float, int, np.floating, np.integer))
            or not np.isfinite(value) or value <= 0 for value in steps)
            or any(a <= b for a, b in zip(steps[:-1], steps[1:]))):
        raise ValueError('Three strictly decreasing positive finite maximum time steps required')
    return r.PathPlan(
        tuple(r.SolverSpec(PRIMARY_METHOD, 1/step, {'max_step_s': step, 'order': 4}) for step in steps),
        plan_for().reference)


def solver_parameters(path, spec):
    parameters = dict(spec.parameters)
    if 'max_step_s' in parameters:
        maximum = parameters.pop('max_step_s')
        parameters['segments'] = int(math.ceil(path.residence_time_s/maximum))
    parameters['sample_count'] = 17
    return parameters


def workload(model):
    rows = []
    seen = set()
    for power in WORKLOAD_POWERS:
        for seed in SEEDS:
            inflow = model.inflow(power, seed)
            keys = {r.digest([inflow.entry_position_m[i], inflow.velocity_m_s[i], inflow.residence_time_s[i]])
                    for i in range(len(inflow.residence_time_s))}
            long = inflow.residence_time_s > 5e-6
            rows.append({'power': power, 'seed': seed, 'path_count': len(keys),
                'new_unique_physical_paths': len(keys-seen),
                'reused_physical_paths': len(keys & seen),
                'residence_range_s': [float(inflow.residence_time_s.min()), float(inflow.residence_time_s.max())],
                'speed_range_m_s': [float(np.linalg.norm(inflow.velocity_m_s, axis=1).min()),
                                    float(np.linalg.norm(inflow.velocity_m_s, axis=1).max())],
                'mean_occupancy': inflow.mean_occupancy, 'equilibrium_nV': inflow.equilibrium_atom_number,
                'occupancy_over_nV': inflow.mean_occupancy/inflow.equilibrium_atom_number,
                'total_arrival_rate_s_inverse': inflow.total_arrival_rate_s_inverse,
                'tau_over_5us_arrival_fraction': float(inflow.rate_s_inverse[long].sum()/inflow.total_arrival_rate_s_inverse),
                'tau_over_5us_occupancy_fraction': float(np.dot(inflow.rate_s_inverse[long], inflow.residence_time_s[long])/inflow.mean_occupancy),
                'status': 'NOT_EVALUATED'})
            seen.update(keys)
    return {'powers': WORKLOAD_POWERS, 'seeds': SEEDS, 'grids': rows,
        'grid_path_entries': sum(row['path_count'] for row in rows),
        'unique_physical_paths': len(seen), 'calculations_per_path': 5,
        'scope': 'Actual physical boundary nodes, no ODE or spectrum convergence inferred from these moments'}


def verified_references(directory, model):
    """Resolve the twelve independently solved jobs by exact physical inputs."""
    snapshot, jobs = refs.prepare_plan()
    output = []
    inflow = model.inflow(0, 11)
    for identity in jobs:
        path_index, level = identity['path_index'], identity['reference_level']
        path = inflow.path(path_index)
        p = model.problem(path)
        data = r.decode(identity)
        if r.encode(model.identity()['inputs_SI']) != identity['inputs_SI']:
            raise ValueError('Reference inputs differ from the declared thermal model')
        for key in ('lower_corner_m', 'upper_corner_m', 'temperature_K', 'mass_kg'):
            if not np.array_equal(data['box'][key], getattr(model, key)):
                raise ValueError('Reference boundary model differs')
        for key in ('entry_position_m', 'velocity_m_s', 'residence_time_s'):
            if not np.array_equal(data['path'][key], getattr(path, key)):
                raise ValueError('Reference is not this exact physical boundary path')
        for key in ('h0', 'h1', 'boundary_state', 'readouts', 'drives'):
            if not np.array_equal(data['atomic_problem'][key], p[key]):
                raise ValueError('Reference atomic problem differs: '+key)
        if (not np.array_equal(data['port_frequencies_rad_s'], p['frequencies_rad_s'])
                or r.encode(data['physical_metadata']) != r.encode(p['metadata'])):
            raise ValueError('Reference phases, RF or optical metadata differ')
        filename = Path(directory)/refs.job_filename(path_index, level)
        record = refs.validate_job(filename, identity, snapshot)
        packet = {**r.decode(record['values']), 'metadata': p['metadata'],
            'source_names': record['source_names'], 'frequencies_rad_s': np.asarray(record['frequencies_rad_s']),
            'residence_time_s': record['residence_time_s'], 'numerics': record['numerics']}
        output.append({'path_index': path_index, 'level': level, 'file': str(filename.resolve()),
            'file_sha256': hashlib.sha256(filename.read_bytes()).hexdigest(),
            'record_sha256': record['record_sha256'], 'identity_sha256': record['full_identity_sha256'],
            'runtime_blas_threads': record['runtime_blas_threads'], 'packet': packet})
    return output


def _primary_worker(job):
    from .reference.exponential_transport import exponential_wavepacket
    path_index, spec_identity, cache_directory, sources = job
    model = r.default_model()
    path = model.inflow(0, 11).path(path_index)
    spec = r.SolverSpec(**spec_identity)
    cache = r.PathCache(cache_directory, sources)

    def provider(model, path, spec):
        p = model.problem(path)
        with core.blas_single_thread():
            pools = [row for row in threadpoolctl.threadpool_info() if row['user_api'] == 'blas']
            if not pools or any(row['num_threads'] != 1 for row in pools):
                raise RuntimeError('Actual one-thread BLAS execution required')
            parameters = solver_parameters(path, spec)
            packet = exponential_wavepacket(**{k: v for k, v in p.items() if k != 'metadata'}, **parameters)
        packet['numerics']['runtime_blas_threads'] = pools
        packet['numerics']['cache_requested_numerical_parameters'] = r.encode(spec.parameters)
        return packet

    _, record = cache.get(model, path, spec, provider)
    return {'path_index': path_index, 'segments': solver_parameters(path, spec)['segments'], **record}


def build_report(reference_directory, cache_directory, *, segments=DEFAULT_SEGMENTS,
                 max_steps_s=None, compute_primary=False, workers=2):
    started = time.monotonic()
    if type(workers) is not int or not 1 <= workers <= 4:
        raise ValueError('workers must be an integer from one to four')
    plan = plan_for(segments) if max_steps_s is None else step_plan(max_steps_s)
    model = r.default_model()
    # Import before the source snapshot even for cache-only runs, so the same
    # mathematical request has one dependency identity in either mode.
    from .reference import exponential_transport
    extra = [Path(__file__), r.ROOT/'tests/quantum/test_rb_thermal_ensemble.py']
    extra += [p for p in (r.ROOT/'tests/quantum/test_rb_thermal_ensemble_audit.py',) if p.exists()]
    sources = r.consumed_hashes(extra)
    references = verified_references(reference_directory, model)
    cache = r.PathCache(cache_directory, sources)
    inflow = model.inflow(0, 11)
    reference_cache = []
    for item in references:
        path = inflow.path(item['path_index'])
        cached, record = cache.get(model, path, plan.reference[item['level']],
            lambda model, path, spec, packet=item['packet']: packet)
        binding_keys = r.METRICS+('metadata', 'numerics', 'frequencies_rad_s', 'residence_time_s', 'source_names')
        if r.digest({key: cached[key] for key in binding_keys}) != r.digest({key: item['packet'][key] for key in binding_keys}):
            raise ValueError('Cached reference differs from the verified immutable independent job')
        reference_cache.append(record)
    primary_jobs = []
    if compute_primary:
        jobs = [(i, spec.identity(), str(Path(cache_directory).resolve()), dict(sources))
                for i in range(6) for spec in plan.primary]
        with ProcessPoolExecutor(max_workers=workers) as pool:
            futures = [pool.submit(_primary_worker, job) for job in jobs]
            for future in as_completed(futures):
                record = future.result()
                primary_jobs.append(record)
                print(f'primary path{record["path_index"]}, segments={record["segments"]}, cache_hit={record["hit"]}', flush=True)
    per_path = []
    for i in range(6):
        path = inflow.path(i)
        p = model.problem(path)
        ref = [references[2*i+level]['packet'] for level in range(2)]
        refinement = r.path_errors(ref[-1], ref[0], p)
        try:
            _, ledger = r.build_path_evidence(model, path, plan, cache)
            status = 'PATH_VERIFIED' if ledger['passed'] else 'PATH_FAILED'
        except (ValueError, KeyError, OSError, RuntimeError) as error:
            ledger = {'passed': False, 'reasons': [type(error).__name__+': '+str(error)]}
            status = 'INCOMPLETE'
        per_path.append({'path_index': i, 'entry_position_m': path.entry_position_m,
            'velocity_m_s': path.velocity_m_s, 'residence_time_s': path.residence_time_s,
            'source_names': model.convention().source_names,
            'independent_reference_refinement': refinement,
            'independent_reference_refinement_passed': all(v <= r.BUDGETS['independent_refinement'] for v in refinement.values()),
            'status': status, 'ledger': ledger})
    candidate, grid = r.run_grid(model, 0, 11, plan, cache)
    gate = r.convergence_gate([grid], powers=WORKLOAD_POWERS, seeds=SEEDS)
    table = workload(model)
    table['grids'][0]['status'] = ('UNCONVERGED' if candidate['path_evidence_passed'] else
        'PATH_FAILED' if any(row['status'] == 'PATH_FAILED' for row in per_path) else 'INCOMPLETE')
    stable = sources == r.consumed_hashes(extra)
    cache._check_sources()
    if not stable:
        raise RuntimeError('Sources changed during pilot; no final artifact may be written')
    return {'schema': 'gabes-thermal-Rb-evidence-pilot-v1', 'physical_model': model.identity(),
        'source_hashes_sha256': sources, 'source_stable_during_run': stable,
        'reference_job_directory': str(Path(reference_directory).resolve()),
        'references': [{k: v for k, v in row.items() if k != 'packet'} for row in references],
        'reference_cache': reference_cache, 'primary_cache_jobs': primary_jobs,
        'declared_plan': {'workload_powers': WORKLOAD_POWERS, 'seeds': SEEDS,
            'executed_pilot': {'power': 0, 'seed': 11, 'all_boundary_paths': 6},
            'path_plan': plan.identity(), 'path_budgets': r.BUDGETS, 'ensemble_budget': r.STREAM_BUDGET,
            'error_definition': r.ERROR_DEFINITION,
            'restriction': 'No other grid/seed is evaluated. The six pilot paths are fixed before solving; no clipping or selected-path substitution.'},
        'workload': table, 'pilot_paths': per_path, 'pilot_grid': grid,
        'ensemble_gate': gate,
        'pilot_all_path_evidence_passed': all(row['ledger']['passed'] for row in per_path),
        'reference_refinements_passed': all(row['independent_reference_refinement_passed'] for row in per_path),
        'thermal_ensemble_converged': False,
        'production_adapter_certified': False, 'physical_optical_prediction': False,
        'environment': {**refs.environment(), 'workers': workers, 'elapsed_s': time.monotonic()-started},
        'limitations': [
            'Square area equals the optical normalization area by assumption; it is not a measured vapor-cell cross-section.',
            'Fresh unpolarized side inflow omits pump history outside this open column; the Gaussian field is already nonzero there.',
            'The workload table is actual boundary geometry, not evaluated atomic spectrum or tail convergence.',
            'A verified six-path grid would still lack the two grid refinements and all independent seed comparisons.',
            'The retarded response is age-integrated atomic response, not nonlocal Maxwell propagation or optical SQL/squeezing.']}


def save_plot(report, filename):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.3), layout='constrained')
    for seed in SEEDS:
        rows = [row for row in report['workload']['grids'] if row['seed'] == seed]
        axes[0].plot([row['power'] for row in rows], [row['occupancy_over_nV'] for row in rows], 'o-', label=f'Seed {seed}')
    axes[0].axhline(1., color='grey', ls='--')
    axes[0].set(xlabel='Sobol power per face', ylabel='Computed occupancy / nV', title='Workload geometry only; no spectrum claim')
    axes[0].legend(fontsize=8)
    axes[1].bar(range(6), [row['residence_time_s']*1e6 for row in report['pilot_paths']], color='#257e81')
    axes[1].set(xlabel='Fixed pilot boundary path', ylabel='Full residence (microseconds)', title='All six p0 / seed11 paths retained')
    for key in r.METRICS:
        axes[2].semilogy(range(6), [max(row['independent_reference_refinement'][key], 1e-17) for row in report['pilot_paths']], 'o-', label=key)
    axes[2].axhline(r.BUDGETS['independent_refinement'], color='black', ls='--')
    axes[2].set(xlabel='Fixed pilot boundary path', ylabel='Max relative error per RF/source', title='Independent adjoint refinement')
    axes[2].legend(fontsize=6, ncols=2)
    for ax in axes:
        ax.grid(alpha=.2)
    fig.suptitle('Conditional open-column Rb pilot | Thermal ensemble NOT CERTIFIED', fontsize=11)
    with Path(filename).open('xb') as handle:
        fig.savefig(handle, format='png', dpi=160)
    plt.close(fig)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True, type=Path)
    parser.add_argument('--plot', required=True, type=Path)
    parser.add_argument('--reference-dir', type=Path, default=refs.DEFAULT_OUTPUT)
    parser.add_argument('--cache-dir', type=Path, required=True)
    numerical = parser.add_mutually_exclusive_group()
    numerical.add_argument('--segments', type=int, nargs=3, default=DEFAULT_SEGMENTS)
    numerical.add_argument('--max-steps-s', type=float, nargs=3)
    parser.add_argument('--compute-primary', action='store_true')
    parser.add_argument('--workers', type=int, default=2)
    args = parser.parse_args(argv)
    if args.output.resolve() == args.plot.resolve():
        raise ValueError('Distinct immutable report and plot paths required')
    for path in (args.output, args.plot):
        if path.exists():
            raise FileExistsError(path)
    report = build_report(args.reference_dir, args.cache_dir, segments=args.segments,
                          max_steps_s=args.max_steps_s, compute_primary=args.compute_primary, workers=args.workers)
    payload = json.dumps(r.encode(report), ensure_ascii=False, indent=2, allow_nan=False)+'\n'
    for path in (args.output, args.plot):
        path.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open('x', encoding='utf-8') as handle:
        handle.write(payload)
    save_plot(report, args.plot)
    print(json.dumps({key: report[key] for key in ('pilot_all_path_evidence_passed',
        'reference_refinements_passed', 'thermal_ensemble_converged')}), flush=True)
    return 0 if report['pilot_all_path_evidence_passed'] and report['reference_refinements_passed'] else 1


if __name__ == '__main__':
    raise SystemExit(main())
