"""Resume all native paths of a declared grid in one frozen-source worker pool.

This controller is outside the numerical source inventory. Its captured bytes
are recorded separately; it never changes the campaign's equations or budgets.
"""

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import time

# Direct script launch needs the checkout (or captured capsule) on sys.path.
# Windows spawn reimports this captured file and the same captured modules.
ROOT = Path(__file__).resolve().parent
if ROOT.name == 'tools':
    ROOT = ROOT.parent
sys.path.insert(0, str(ROOT))
from analysis.grand_challenge import thermal_campaign as c


def declared_selections(plan, power, seed):
    if type(power) is not int or type(seed) is not int:
        raise ValueError('Integer declared grid required')
    grids = [g for g in plan['grids'] if (g['power'], g['seed']) == (power, seed)]
    if len(grids) != 1:
        raise ValueError('Select one declared grid')
    return [(power, seed, i) for i in range(len(grids[0]['paths']))]


def controller_identity():
    raw = Path(__file__).read_bytes()
    return {'name': 'thermal_campaign_batch.py',
            'raw_sha256': hashlib.sha256(raw).hexdigest(), 'raw_bytes': len(raw)}


def prepare_jobs(directory, plan, selections):
    """Validate all available records before scheduling any missing calculation."""
    cache = c.CampaignCache(Path(directory)/'cache', plan)
    numerical = c.pilot.step_plan(plan['max_steps_s'])
    specs = numerical.primary+numerical.reference
    existing, missing, keys = {}, [], set()
    for selection in selections:
        model, path = c.selected_path(plan, *selection)
        for level, spec in enumerate(specs):
            key = cache.key(model, path, spec)
            if key in keys:
                raise ValueError('Duplicate physical jobs in selected grid')
            keys.add(key)
            # Exact mathematical reuse: these records solve the same unweighted
            # path equations at identical tolerances. Grid labels/arrival rates
            # enter only the subsequent stream sum. Re-solving is unnecessary;
            # all five numerical comparisons are recomputed before publication.
            if (cache.directory/(key+'.json')).exists():
                packet, record = cache.get(model, path, spec)
                existing[(selection, level)] = {'selection': selection, 'level': level,
                    'elapsed_s': packet['numerics']['elapsed_seconds'], **record}
            else:
                # Longest estimated work first reduces the pool's idle tail.
                cost = path.residence_time_s * (1., 2., 4., 1.5, 1.6)[level]
                missing.append((cost, selection, level))
    missing.sort(key=lambda item: (-item[0], item[1], item[2]))
    return cache, existing, missing


def publish_path(directory, plan, selection, results, cache):
    model, path = c.selected_path(plan, *selection)
    _, ledger = c.r.build_path_evidence(model, path,
                                       c.pilot.step_plan(plan['max_steps_s']), cache)
    destination = Path(directory)/('path-p%d-s%d-i%d.json' % selection)
    if destination.exists():
        old = c.records.read_sealed(destination)
        if (old.get('schema') != 'gabes-frozen-thermal-path-audit-v1'
                or old.get('campaign_sha256') != plan['record_sha256']
                or old.get('selection') != list(selection)
                or old.get('physical_path') != c.path_identity(path)
                or old.get('path_passed') is not ledger['passed']
                or old.get('thermal_ensemble_converged') is not False
                or old.get('physical_optical_prediction') is not False
                or c.provenance.source_identity(old['source_manifest_after']) != plan['source_identity']):
            raise ValueError('Existing path report disagrees with current evidence')
        for name in ('target_digest', 'primary_refinement', 'independent_reference',
                     'independent_refinement', 'passed', 'reasons'):
            if old['ledger'][name] != ledger[name]:
                raise ValueError('Existing path comparison differs: '+name)
        report = old
        reused = True
    else:
        report = c.write_record(destination, {
            'schema': 'gabes-frozen-thermal-path-audit-v1', 'selection': selection,
            'campaign_sha256': plan['record_sha256'], 'physical_path': c.path_identity(path),
            'jobs': [results[(selection, level)] for level in range(5)], 'ledger': ledger,
            'path_passed': ledger['passed'], 'thermal_ensemble_converged': False,
            'physical_optical_prediction': False, 'source_manifest_after': cache._check_sources(),
            'controller': controller_identity()})
        reused = False
    return {'selection': selection, 'path_passed': ledger['passed'],
            'report': destination.name, 'record_sha256': report['record_sha256'],
            'existing_report_preserved': reused}


def execute(directory, power, seed, output, workers=4):
    if type(workers) is not int or not 1 <= workers <= 4:
        raise ValueError('One to four worker processes required')
    directory, output = Path(directory).resolve(), Path(output).resolve()
    if output.exists():
        raise FileExistsError(output)
    plan = c.load_plan(directory)
    c.require_capsule(plan)
    control = controller_identity()
    started = time.monotonic()
    selections = declared_selections(plan, power, seed)
    cache, results, missing = prepare_jobs(directory, plan, selections)
    initial_hits = len(results)
    print(json.dumps({'event': 'prepared', 'grid': [power, seed],
                      'cached_jobs': initial_hits, 'new_jobs': len(missing),
                      'workers': workers}), flush=True)
    published = {}

    def completed_paths():
        for selection in selections:
            if selection not in published and all((selection, level) in results for level in range(5)):
                published[selection] = publish_path(directory, plan, selection, results, cache)
                print(json.dumps({'event': 'path', **published[selection]}), flush=True)

    completed_paths()
    with ProcessPoolExecutor(max_workers=workers) as pool:
        pending = {pool.submit(c.execute_job, str(directory), selection, level): (selection, level)
                   for _, selection, level in missing}
        for future in as_completed(pending):
            key = pending[future]
            result = future.result()
            if (tuple(result['selection']), result['level']) != key:
                raise ValueError('Worker returned a different job')
            results[key] = result
            completed_paths()
    if len(published) != len(selections) or controller_identity() != control:
        raise ValueError('Incomplete grid or changed execution controller')
    payload = {'schema': 'gabes-frozen-thermal-batch-v1', 'grid': [power, seed],
        'campaign_sha256': plan['record_sha256'], 'source_identity': plan['source_identity'],
        'source_manifest_after': cache._check_sources(), 'controller': control,
        'workers': workers, 'initial_cache_hits': initial_hits,
        'new_solve_count': sum(not row['hit'] for row in results.values()),
        'job_count': len(results), 'wall_seconds': time.monotonic()-started,
        'paths': [published[s] for s in selections],
        'all_paths_passed': all(row['path_passed'] for row in published.values()),
        'thermal_ensemble_converged': False, 'physical_optical_prediction': False}
    return c.write_record(output, payload)


def launch(directory, power, seed, output, workers=4):
    directory, output = Path(directory).resolve(), Path(output).resolve()
    if output.exists():
        raise FileExistsError(output)
    scratch = ROOT/'.git'/'grand-challenge-runs'
    scratch.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='thermal-batch-', dir=scratch) as folder:
        c.extract_bundle(directory, folder)
        controller = Path(folder)/'thermal_campaign_batch.py'
        shutil.copyfile(__file__, controller)
        env = os.environ.copy()
        env.pop('PYTHONPATH', None)
        env.update(OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1', OMP_NUM_THREADS='1')
        subprocess.run([sys.executable, str(controller), '--execute', str(directory),
                        '--grid', str(power), str(seed), '--output', str(output),
                        '--workers', str(workers)], cwd=folder, env=env, check=True)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    action = parser.add_mutually_exclusive_group(required=True)
    action.add_argument('--run', type=Path)
    action.add_argument('--execute', type=Path, help=argparse.SUPPRESS)
    parser.add_argument('--grid', nargs=2, type=int, required=True, metavar=('POWER', 'SEED'))
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--workers', type=int, default=4)
    args = parser.parse_args(argv)
    if args.run:
        launch(args.run, *args.grid, args.output, args.workers)
        return 0 if c.records.read_sealed(args.output)['all_paths_passed'] else 1
    report = execute(args.execute, *args.grid, args.output, args.workers)
    return 0 if report['all_paths_passed'] else 1


if __name__ == '__main__':
    raise SystemExit(main())
