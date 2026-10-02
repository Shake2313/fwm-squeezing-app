"""Frozen-source thermal Rb jobs, independent of historical selected-path runs.

Create a campaign once; execute any declared boundary path in an isolated source
bundle. Missing paths stay missing. Source equivalence never certifies physics.
"""

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import zipfile

import numpy as np
import threadpoolctl

from gabes import core
from . import rb_thermal_ensemble as r
from . import rb_thermal_ensemble_audit as pilot
from . import rb_thermal_reference_jobs as records
from . import source_provenance as provenance
from .reference.adjoint_transport import adjoint_wavepacket
from .reference.exponential_transport import exponential_wavepacket


ROOT = Path(__file__).resolve().parents[2]
MAX_STEPS = (2e-6/16384, 2e-6/32768, 2e-6/65536)
SCHEMA = 'gabes-frozen-thermal-campaign-v1'


def source_paths(root):
    # Deliberately broad, explicit inventory. A frozen bundle avoids concurrent
    # app edits without pretending those edits share historical source hashes.
    return sorted(p.relative_to(root).as_posix()
                  for directory in ('gabes', 'analysis/grand_challenge')
                  for p in (Path(root)/directory).rglob('*.py')
                  if '__pycache__' not in p.parts)


def write_record(path, payload):
    record = records.sealed(r.encode(payload))
    text = json.dumps(record, ensure_ascii=False, allow_nan=False, indent=2)+'\n'
    destination = Path(path)
    if destination.exists():
        raise FileExistsError(destination)
    # A killed/interrupted writer must never leave a partial immutable record
    # at its final pathname. The same-directory temporary also keeps the final
    # publication on one filesystem. Existing evidence must not be replaced,
    # including a record published by another writer after the preflight check.
    descriptor, name = tempfile.mkstemp(prefix='.'+destination.name+'.',
                                        suffix='.tmp', dir=destination.parent)
    temporary = Path(name)
    try:
        with os.fdopen(descriptor, 'w', encoding='utf-8', newline='\n') as handle:
            handle.write(text)
            handle.flush()
            os.fsync(handle.fileno())
        if os.name == 'nt':
            # Windows rename refuses an existing destination. Drive-backed
            # Windows folders need not support hardlinks.
            os.rename(temporary, destination)
        else:
            # POSIX rename would overwrite; linking atomically reserves the
            # final name without replacing another writer's complete record.
            os.link(temporary, destination)
    finally:
        temporary.unlink(missing_ok=True)
    return record


def path_identity(path):
    return r.encode({'entry_position_m': path.entry_position_m,
                     'velocity_m_s': path.velocity_m_s,
                     'residence_time_s': path.residence_time_s})


def make_plan(root=ROOT, *, powers=(0, 1, 2), seeds=(11, 211, 811), max_steps=MAX_STEPS):
    powers, seeds = tuple(powers), tuple(seeds)
    if (len(powers) < 3 or any(type(p) is not int or not 0 <= p <= 10 for p in powers)
            or tuple(sorted(set(powers))) != powers
            or len(seeds) < 3 or len(set(seeds)) != len(seeds)
            or any(type(s) is not int or s < 0 for s in seeds)):
        raise ValueError('At least three increasing bounded powers and three distinct nonnegative seeds required')
    numerical = pilot.step_plan(max_steps)
    model = r.default_model()
    paths = source_paths(root)
    snapshot = provenance.capture_sources(root, paths)
    grids = []
    for power in powers:
        for seed in seeds:
            inflow = model.inflow(power, seed)
            grids.append({'power': power, 'seed': seed,
                'paths': [{'index': i, 'physical_path': path_identity(inflow.path(i)),
                           'rate_s_inverse': float(inflow.rate_s_inverse[i])}
                          for i in range(len(inflow.rate_s_inverse))],
                'mean_occupancy': inflow.mean_occupancy,
                'equilibrium_nV': inflow.equilibrium_atom_number})
    provenance.require_sources(root, snapshot)
    return {'schema': SCHEMA, 'source_manifest': snapshot, 'source_paths': paths,
        'source_identity': provenance.source_identity(snapshot),
        'model': r.encode(model.identity()), 'environment': records.environment(),
        'powers': powers, 'seeds': seeds, 'max_steps_s': tuple(max_steps),
        'path_plan': numerical.identity(), 'grids': grids,
        'path_budgets': r.BUDGETS, 'ensemble_budget': r.STREAM_BUDGET,
        'historical_numerical_results_reused': False,
        'scope': 'Conditional open-column pump-only atomic stream; no optical squeezing certification'}


def create_campaign(directory, *, root=ROOT, **kwargs):
    directory = Path(directory)
    if directory.exists():
        raise FileExistsError(directory)
    plan = make_plan(root, **kwargs)
    directory.mkdir(parents=True)
    archive = directory/'sources.zip'
    with zipfile.ZipFile(archive, 'x', compression=zipfile.ZIP_DEFLATED) as bundle:
        for name in plan['source_paths']:
            bundle.writestr(name, (Path(root)/name).read_bytes())
    provenance.require_sources(root, plan['source_manifest'])
    if source_paths(root) != plan['source_paths']:
        raise ValueError('Source inventory changed while freezing campaign')
    plan['source_bundle_sha256'] = hashlib.sha256(archive.read_bytes()).hexdigest()
    write_record(directory/'plan.json', plan)
    (directory/'cache').mkdir()
    return directory/'plan.json'


def load_plan(directory, *, root=ROOT):
    plan = records.read_sealed(Path(directory)/'plan.json')
    if (plan.get('schema') != SCHEMA
            or plan.get('source_identity') != provenance.source_identity(plan['source_manifest'])
            or plan.get('environment') != records.environment()
            or plan.get('model') != r.encode(r.default_model().identity())
            or plan.get('path_plan') != r.encode(pilot.step_plan(plan['max_steps_s']).identity())
            or plan.get('path_budgets') != r.BUDGETS
            or plan.get('ensemble_budget') != r.STREAM_BUDGET):
        raise ValueError('Campaign model, numerical plan, environment or source identity differs')
    provenance.require_sources(root, plan['source_manifest'])
    if source_paths(root) != plan['source_paths']:
        raise ValueError('Frozen source inventory differs')
    return plan


def selected_path(plan, power, seed, index):
    rows = [g for g in plan['grids'] if (g['power'], g['seed']) == (power, seed)]
    if len(rows) != 1 or type(index) is not int or not 0 <= index < len(rows[0]['paths']):
        raise ValueError('Select an exact declared grid and native path index')
    model = r.default_model()
    inflow = model.inflow(power, seed)
    path = inflow.path(index)
    expected = rows[0]['paths'][index]
    if (expected['index'] != index or expected['physical_path'] != path_identity(path)
            or expected['rate_s_inverse'] != float(inflow.rate_s_inverse[index])):
        raise ValueError('Declared boundary path or arrival rate changed')
    return model, path


def numerical_parameters(path, spec):
    if spec.method == pilot.PRIMARY_METHOD:
        return pilot.solver_parameters(path, spec)
    if spec.method == pilot.REFERENCE_METHOD:
        return {'rtol': spec.parameters['rtol'], 'atol': spec.parameters['atol'],
                'max_step_s': path.residence_time_s*spec.parameters['max_step_fraction']}
    raise ValueError('Unknown native mathematical solver')


def native_provider(model, path, spec):
    p = model.problem(path)
    solver = exponential_wavepacket if spec.method == pilot.PRIMARY_METHOD else adjoint_wavepacket
    parameters = numerical_parameters(path, spec)
    with core.blas_single_thread():
        pools = [row for row in threadpoolctl.threadpool_info() if row['user_api'] == 'blas']
        if not pools or any(row['num_threads'] != 1 for row in pools):
            raise RuntimeError('Measured single-thread BLAS required')
        raw = solver(**{k: v for k, v in p.items() if k != 'metadata'}, **parameters)
    # Store measured runtime evidence, never an assumed thread-count label.
    raw['numerics']['runtime_blas_threads'] = pools
    raw['numerics']['requested_parameters'] = parameters
    return raw


def validate_native(packet, path, spec):
    numbers = packet['numerics']
    expected = numerical_parameters(path, spec)
    if numbers.get('requested_parameters') != expected or numbers.get('synthetic_test_fixture'):
        raise ValueError('Native numerical parameters differ or synthetic record supplied')
    pools = numbers.get('runtime_blas_threads', [])
    if not pools or any(p.get('user_api') != 'blas' or p.get('num_threads') != 1 for p in pools):
        raise ValueError('Measured single-thread BLAS evidence required')
    if spec.method == pilot.PRIMARY_METHOD:
        if any(numbers.get(k) != expected[k] for k in ('segments', 'order')):
            raise ValueError('Primary solver work differs from requested resolution')
        if numbers.get('substeps') != 2*expected['segments']:
            raise ValueError('Incomplete CF4 substep work')
    else:
        if any(numbers.get(k) != v for k, v in expected.items()):
            raise ValueError('Independent solver tolerance or native residence differs')
        for k in ('forward_density_evaluations', 'backward_evaluations', 'density_mesh_points'):
            if type(numbers.get(k)) is not int or numbers[k] <= 0:
                raise ValueError('Positive independent solver work counters required')
    elapsed = numbers.get('elapsed_seconds')
    if isinstance(elapsed, bool) or not isinstance(elapsed, (float, int)) or not np.isfinite(elapsed) or elapsed < 0:
        raise ValueError('Finite nonnegative measured solve time required')


class CampaignCache:
    def __init__(self, directory, plan, *, root=ROOT):
        self.directory, self.root = Path(directory), Path(root)
        # JSON round trip prevents caller mutation of a frozen request.
        self.plan = json.loads(json.dumps(r.encode(plan), allow_nan=False))
        if (self.plan['source_identity'] != provenance.source_identity(self.plan['source_manifest'])
                or self.plan['environment'] != records.environment()):
            raise ValueError('Cache source identity or runtime environment differs')
        self._check_sources()

    def _check_sources(self):
        if source_paths(self.root) != self.plan['source_paths']:
            raise ValueError('Frozen dependency inventory differs')
        return provenance.require_sources(self.root, self.plan['source_manifest'])

    def identity(self, model, path, spec):
        # Identical physical paths on nested grids solve identical equations.
        # Grid labels/rates are omitted only here; stream rates and packet
        # digests are rebound and verified by build_path_evidence/run_grid.
        return {'schema': 'gabes-portable-native-thermal-path-v1',
            'source_identity': self.plan['source_identity'],
            'model': model.identity(), 'path': path_identity(path),
            'solver': spec.identity(), 'environment': self.plan['environment']}

    def key(self, model, path, spec):
        return r.digest(self.identity(model, path, spec))

    def get(self, model, path, spec, provider=None):
        before = self._check_sources()
        identity = r.encode(self.identity(model, path, spec))
        key = r.digest(identity)
        filename = self.directory/(key+'.json')
        hit = filename.exists()
        if hit:
            record = records.read_sealed(filename)
            if (record.get('schema') != 'gabes-native-thermal-packet-v1'
                    or record.get('identity') != identity
                    or record.get('key') != key
                    or provenance.source_identity(record['source_before']) != self.plan['source_identity']
                    or provenance.source_identity(record['source_after']) != self.plan['source_identity']
                    or record.get('source_stable_during_run') is not True):
                raise ValueError('Native cache identity or execution provenance differs')
            packet = r._validated_packet(model, path, r.decode(record['packet']))
            validate_native(packet, path, spec)
        else:
            if provider is None:
                raise FileNotFoundError('Native campaign path not computed: '+key)
            if provider is not native_provider:
                raise ValueError('Campaign publication requires the native primary or independent solver')
            require_capsule(self.plan)
            packet = r._validated_packet(model, path, provider(model, path, spec))
            validate_native(packet, path, spec)
            after = self._check_sources()
            self.directory.mkdir(parents=True, exist_ok=True)
            record = write_record(filename, {'schema': 'gabes-native-thermal-packet-v1',
                'identity': identity, 'key': key, 'original_path_source': path.source,
                'packet': packet, 'source_before': before, 'source_after': after,
                'source_stable_during_run': provenance.source_identity(before) == provenance.source_identity(after)})
        self._check_sources()
        return packet, {'key': key, 'hit': hit, 'record_sha256': record['record_sha256'],
                        'path': str(filename), 'payload_digest': r.digest(packet)}


def execute_job(directory, selection, level):
    plan = load_plan(directory)
    require_capsule(plan)
    if type(level) is not int or level not in range(5):
        raise ValueError('Native job level must be 0 through 4')
    model, path = selected_path(plan, *selection)
    numerical = pilot.step_plan(plan['max_steps_s'])
    spec = (numerical.primary+numerical.reference)[level]
    cache = CampaignCache(Path(directory)/'cache', plan)
    packet, record = cache.get(model, path, spec, native_provider)
    result = {'selection': selection, 'level': level, 'elapsed_s': packet['numerics']['elapsed_seconds'], **record}
    print(json.dumps(result), flush=True)
    return result


def execute(directory, selection, workers=4):
    if type(workers) is not int or not 1 <= workers <= 4:
        raise ValueError('One to four worker processes required')
    directory = Path(directory).resolve()
    plan = load_plan(directory)
    require_capsule(plan)
    model, path = selected_path(plan, *selection)
    name = 'path-p%d-s%d-i%d.json' % tuple(selection)
    destination = directory/name
    if destination.exists():
        raise FileExistsError(destination)
    with ProcessPoolExecutor(max_workers=workers) as pool:
        jobs = [pool.submit(execute_job, str(directory), selection, level) for level in range(5)]
        results = [job.result() for job in as_completed(jobs)]
    cache = CampaignCache(directory/'cache', plan)
    _, ledger = r.build_path_evidence(model, path, pilot.step_plan(plan['max_steps_s']), cache)
    # Finishing a selected path cannot imply that its full grid was integrated.
    report = {'schema': 'gabes-frozen-thermal-path-audit-v1', 'selection': selection,
        'campaign_sha256': plan['record_sha256'], 'physical_path': path_identity(path),
        'jobs': sorted(results, key=lambda row: row['level']), 'ledger': ledger,
        'path_passed': ledger['passed'], 'thermal_ensemble_converged': False,
        'physical_optical_prediction': False, 'source_manifest_after': cache._check_sources()}
    write_record(destination, report)
    return report


def extract_bundle(directory, target):
    directory, target = Path(directory), Path(target)
    plan = records.read_sealed(directory/'plan.json')
    archive = directory/'sources.zip'
    if hashlib.sha256(archive.read_bytes()).hexdigest() != plan['source_bundle_sha256']:
        raise ValueError('Frozen source bundle changed')
    with zipfile.ZipFile(archive) as bundle:
        if sorted(bundle.namelist()) != plan['source_paths'] or len(bundle.namelist()) != len(set(bundle.namelist())):
            raise ValueError('Frozen source bundle inventory differs')
        for name in bundle.namelist():
            destination = (target/name).resolve()
            if not destination.is_relative_to(target.resolve()):
                raise ValueError('Frozen source escapes execution directory')
            destination.parent.mkdir(parents=True, exist_ok=True)
            with destination.open('xb') as handle:
                handle.write(bundle.read(name))
    provenance.require_sources(target, plan['source_manifest'])
    with (target/'.campaign-source-identity').open('x', encoding='utf-8') as handle:
        handle.write(plan['record_sha256'])


def require_capsule(plan):
    marker = ROOT/'.campaign-source-identity'
    if not marker.exists() or marker.read_text(encoding='utf-8') != plan['record_sha256']:
        raise ValueError('Execute through --run in a fresh captured-source interpreter')


def launch(directory, selection, workers=4):
    directory = Path(directory).resolve()
    # New interpreter/process pool imports the captured sources, never the live
    # application's already-imported modules. Each job also verifies its files.
    scratch = ROOT/'.git'/'grand-challenge-runs'
    scratch.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='thermal-', dir=scratch) as folder:
        extract_bundle(directory, folder)
        env = os.environ.copy()
        env.pop('PYTHONPATH', None)
        command = [sys.executable, '-m', 'analysis.grand_challenge.thermal_campaign',
                   '--execute', str(directory), '--path', *map(str, selection), '--workers', str(workers)]
        subprocess.run(command, cwd=folder, env=env, check=True)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    action = parser.add_mutually_exclusive_group(required=True)
    action.add_argument('--create', type=Path)
    action.add_argument('--run', type=Path)
    action.add_argument('--execute', type=Path, help=argparse.SUPPRESS)
    parser.add_argument('--path', nargs=3, type=int, metavar=('POWER', 'SEED', 'INDEX'))
    parser.add_argument('--workers', type=int, default=4)
    args = parser.parse_args(argv)
    if args.create:
        print(create_campaign(args.create), flush=True)
    else:
        if args.path is None:
            parser.error('--path POWER SEED INDEX required')
        if args.run:
            launch(args.run, args.path, args.workers)
        else:
            result = execute(args.execute, args.path, args.workers)
            print(json.dumps({'path_passed': result['path_passed'], 'thermal_ensemble_converged': False}), flush=True)
            return 0 if result['path_passed'] else 1
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
