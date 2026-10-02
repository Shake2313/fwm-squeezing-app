"""Declare a finer frozen thermal campaign and transfer exact raw cache records.

Parent files are read-only. Publication requires a second fresh interpreter
using the staged campaign's own capsule. No atomic provider is ever supplied.
"""

import argparse
import ctypes
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parent
if ROOT.name == 'tools':
    ROOT = ROOT.parent
sys.path.insert(0, str(ROOT))
from tools import thermal_campaign_grid as g

c = g.c

SCHEMA = 'gabes-frozen-thermal-transfer-v1'


def sha(raw):
    return hashlib.sha256(raw).hexdigest()


def controller():
    if Path(g.__file__).resolve().parent != Path(__file__).resolve().parent:
        raise ValueError('Grid helper must originate beside the extension controller')
    raw = Path(__file__).read_bytes()
    token = sha(raw)
    return {'name': 'thermal_campaign_extend.py', 'raw_sha256': token,
        'raw_bytes': len(raw), 'archive': 'controllers/thermal_campaign_extend-'+token+'.py',
        'dependencies': [g.controller_identity()]}, raw


def controller_files(control, raw, *, stage=None):
    dependency = control['dependencies'][0]
    helper = (Path(stage)/dependency['archive']).read_bytes() if stage else Path(g.__file__).read_bytes()
    if sha(helper) != dependency['raw_sha256'] or len(helper) != dependency['raw_bytes']:
        raise ValueError('Captured grid helper changed')
    return [(control, raw), (dependency, helper)]


def copy_controllers(folder, captured):
    directory = Path(folder)/'tools'
    directory.mkdir()
    for identity, raw in captured:
        with (directory/identity['name']).open('xb') as handle:
            handle.write(raw)
    return directory/'thermal_campaign_extend.py'


def destinations(parent, target):
    parent, requested = Path(parent).resolve(), Path(target)
    if os.path.lexists(requested):
        raise FileExistsError(requested)
    target = requested.resolve()
    if target.is_relative_to(parent) or parent.is_relative_to(target):
        raise ValueError('Extension target must be separate from the immutable parent')
    if not target.parent.is_dir():
        raise FileNotFoundError('Target parent directory must already exist: '+str(target.parent))
    return parent, target


def source_origins(plan):
    expected = {row['path']: row for row in plan['source_manifest']['files']}
    rows = []
    for name, module in sorted(tuple(sys.modules.items())):
        if name != 'gabes' and not name.startswith(('gabes.', 'analysis.grand_challenge.')):
            continue
        filename = getattr(module, '__file__', None)
        if filename is None:
            continue
        path = Path(filename).resolve()
        if not path.is_relative_to(c.ROOT.resolve()):
            raise ValueError('Loaded numerical source originates outside the captured capsule')
        relative = path.relative_to(c.ROOT.resolve()).as_posix()
        if relative not in expected or sha(path.read_bytes()) != expected[relative]['raw_sha256']:
            raise ValueError('Loaded numerical source differs from the frozen raw inventory')
        rows.append({'module': name, 'source_path': relative, 'raw_sha256': expected[relative]['raw_sha256']})
    for filename in (Path(__file__), Path(g.__file__)):
        if not filename.resolve().is_relative_to(c.ROOT.resolve()):
            raise ValueError('Extension controller loaded outside the captured capsule')
    return rows


def snapshot(parent, plan):
    plan_raw, bundle = (parent/'plan.json').read_bytes(), (parent/'sources.zip').read_bytes()
    if c.records.read_sealed(parent/'plan.json') != plan or sha(bundle) != plan['source_bundle_sha256']:
        raise ValueError('Parent declaration or source archive changed')
    return {'campaign_sha256': plan['record_sha256'], 'plan_file_sha256': sha(plan_raw),
            'source_bundle_sha256': sha(bundle)}


def compatible(parent, child):
    for name in ('source_manifest', 'source_paths', 'source_identity', 'model', 'environment',
                 'max_steps_s', 'path_plan', 'path_budgets', 'ensemble_budget'):
        if c.r.encode(parent[name]) != c.r.encode(child[name]):
            raise ValueError('Extension changes immutable equations or execution contract: '+name)
    if tuple(child['seeds']) != tuple(parent['seeds']):
        raise ValueError('Extension must preserve every declared independent seed')


def job_index(plan, cache):
    """Enumerate and validate native declarations before trusting filenames."""
    powers, seeds = tuple(plan['powers']), tuple(plan['seeds'])
    if (len(powers) < 3 or any(type(p) is not int or not 0 <= p <= 10 for p in powers)
            or tuple(sorted(set(powers))) != powers
            or len(seeds) < 3 or len(set(seeds)) != len(seeds)
            or any(type(seed) is not int or seed < 0 for seed in seeds)):
        raise ValueError('Invalid complete campaign powers or seeds')
    declared = {(row['power'], row['seed']): row for row in plan['grids']}
    if len(declared) != len(plan['grids']) or set(declared) != {(p, s) for p in powers for s in seeds}:
        raise ValueError('Campaign must declare every unique grid and seed')
    numerical = c.pilot.step_plan(plan['max_steps_s'])
    indexed = {}
    for power in powers:
        for seed in seeds:
            # Validate ALL parent and target counts/rates/coordinates/occupancy,
            # including grids with no cached jobs. load_plan alone does not.
            model, inflow = g.selected_grid(plan, power, seed)
            for index in range(len(inflow.rate_s_inverse)):
                path = inflow.path(index)
                for level, spec in enumerate(numerical.primary+numerical.reference):
                    key = cache.key(model, path, spec)
                    indexed.setdefault(key, (model, path, spec, [power, seed, index, level]))
    return indexed


def regular_cache_file(path):
    # Virtual filesystems may report zero when link-count metadata is unavailable;
    # zero is not evidence of a hardlink. Reject observed multiple links. Imports
    # still create independent files exclusively, never links, and verify bytes.
    return not path.is_symlink() and path.is_file() and path.stat().st_nlink in (0, 1)


def cache_files(directory):
    files = {}
    for path in sorted(Path(directory).iterdir()):
        key = path.stem
        if (not regular_cache_file(path) or path.suffix != '.json'
                or len(key) != 64 or any(char not in '0123456789abcdef' for char in key)):
            raise ValueError('Unexpected or nonregular file in immutable packet cache: '+path.name)
        files[key] = path
    return files


def present_snapshot(parent, origin, control):
    rows = []
    for key, filename in cache_files(parent/'cache').items():
        raw = filename.read_bytes()
        rows.append({'key': key, 'file_sha256': sha(raw), 'raw_bytes': len(raw)})
    return c.records.sealed({'schema': 'gabes-thermal-cache-snapshot-v1',
                            'parent': origin, 'controller': control, 'records': rows})


def read_parent_records(parent, parent_index, child_index, cache, captured):
    imported, excluded = [], []
    rows = captured['records']
    if len({row['key'] for row in rows}) != len(rows):
        raise ValueError('Duplicate captured cache key')
    for captured_row in rows:
        key = captured_row['key']
        if key not in parent_index:
            raise ValueError('Parent cache key is outside its declared native jobs: '+key)
        filename = parent/'cache'/(key+'.json')
        raw = filename.read_bytes()
        if sha(raw) != captured_row['file_sha256'] or len(raw) != captured_row['raw_bytes']:
            raise ValueError('Captured parent packet changed after the launch snapshot')
        model, path, spec, selection = parent_index[key]
        _, record = cache.get(model, path, spec)  # provider=None: immutable read only.
        if not record['hit'] or filename.read_bytes() != raw:
            raise ValueError('Parent packet changed during validation')
        row = {'key': key, 'record_sha256': record['record_sha256'],
            'payload_digest': record['payload_digest'], 'file_sha256': sha(raw), 'raw_bytes': len(raw),
            'parent_selection': selection, 'cache_file': 'cache/'+key+'.json'}
        if key in child_index:
            row['target_selection'] = child_index[key][3]
            imported.append(row)
        else:
            excluded.append(row)
    return imported, excluded


def verify_stage(stage, *, expected_controller_sha256, expected_grid_sha256=None):
    """Run only inside the staged declaration's newly extracted interpreter."""
    stage = Path(stage).resolve()
    plan = c.load_plan(stage)
    c.require_capsule(plan)
    control, raw = controller()
    if control['raw_sha256'] != expected_controller_sha256:
        raise ValueError('Verification controller differs from captured bytes')
    if expected_grid_sha256 is not None and control['dependencies'][0]['raw_sha256'] != expected_grid_sha256:
        raise ValueError('Verification grid helper differs from captured bytes')
    for identity, contents in controller_files(control, raw):
        if (stage/identity['archive']).read_bytes() != contents:
            raise ValueError('Staged controller archive changed')
    request = c.records.read_sealed(stage/'transfer-request.json')
    parent = Path(request['parent_location']).resolve()
    old = c.load_plan(parent)
    compatible(old, plan)
    if (request['schema'] != SCHEMA or request['controller'] != control
            or request['target_campaign_sha256'] != plan['record_sha256']
            or request['parent'] != snapshot(parent, old)
            or plan['extension']['parent'] != request['parent']
            or plan['extension']['controller'] != control
            or plan['extension']['schema'] != 'gabes-frozen-thermal-extension-v1'
            or plan['extension']['path_and_grid_evidence_imported'] is not False
            or request['path_and_grid_evidence_imported'] is not False
            or type(request['new_solve_count']) is not int or request['new_solve_count'] != 0):
        raise ValueError('Transfer request disagrees with parent or staged declaration')
    cache = c.CampaignCache(stage/'cache', plan)
    indexed = job_index(plan, cache)
    old_cache = c.CampaignCache(parent/'cache', old)
    old_index = job_index(old, old_cache)
    imports = request['imports']
    keys = [row['key'] for row in imports]
    captured = request['parent_cache_snapshot']
    captured_payload = {key: value for key, value in captured.items() if key != 'record_sha256'}
    if (c.records.sealed(captured_payload) != captured or captured['parent'] != request['parent']
            or captured['controller'] != control or captured['schema'] != 'gabes-thermal-cache-snapshot-v1'):
        raise ValueError('Parent cache snapshot seal or provenance differs')
    captured_keys = [row['key'] for row in captured['records']]
    considered = imports+request['excluded_parent_records']
    if (len(set(captured_keys)) != len(captured_keys) or len(considered) != len(captured_keys)
            or {row['key'] for row in considered} != set(captured_keys)):
        raise ValueError('Transfer omits or duplicates a captured parent record')
    captured_by_key = {row['key']: row for row in captured['records']}
    for row in considered:
        observed = captured_by_key[row['key']]
        if row['file_sha256'] != observed['file_sha256'] or row['raw_bytes'] != observed['raw_bytes']:
            raise ValueError('Transfer file metadata differs from the launch snapshot')
    for row in request['excluded_parent_records']:
        key = row['key']
        if (key in indexed or key not in old_index or row['cache_file'] != 'cache/'+key+'.json'
                or row['parent_selection'] != old_index[key][3]):
            raise ValueError('Excluded record belongs to the target or is misbound')
        model, path, spec, _ = old_index[key]
        _, record = old_cache.get(model, path, spec)
        if record['record_sha256'] != row['record_sha256'] or record['payload_digest'] != row['payload_digest']:
            raise ValueError('Excluded record metadata differs from its validated parent record')
    for row in captured['records']:
        filename = parent/'cache'/(row['key']+'.json')
        original = filename.read_bytes()
        if (row['key'] not in old_index or not regular_cache_file(filename)
                or sha(original) != row['file_sha256'] or len(original) != row['raw_bytes']):
            raise ValueError('Captured parent record changed during fresh verification')
    if (len(set(keys)) != len(keys) or set(cache_files(stage/'cache')) != set(keys)
            or request['expected_native_keys'] != sorted(indexed)
            or request['requested_powers'] != plan['powers'] or request['required_seeds'] != plan['seeds']
            or request['missing_native_keys'] != sorted(set(indexed)-set(keys))
            or plan['extension']['imported_raw_record_count'] != len(keys)
            or plan['historical_numerical_results_reused'] is not bool(keys)):
        raise ValueError('Staged cache inventory differs from the requested transfer')
    validated = []
    for row in imports:
        key = row['key']
        if (key not in indexed or key not in old_index
                or row['cache_file'] != 'cache/'+key+'.json'
                or row['parent_selection'] != old_index[key][3]
                or row['target_selection'] != indexed[key][3]):
            raise ValueError('Transfer cache reference is not bound to declared native jobs')
        raw_source = (parent/'cache'/(key+'.json')).read_bytes()
        copied = (stage/'cache'/(key+'.json')).read_bytes()
        if raw_source != copied or sha(copied) != row['file_sha256'] or len(copied) != row['raw_bytes']:
            raise ValueError('Transferred record bytes differ from the immutable parent')
        model, path, spec, _ = indexed[key]
        _, record = cache.get(model, path, spec)  # no provider, repair, or new solve.
        if (not record['hit'] or record['record_sha256'] != row['record_sha256']
                or record['payload_digest'] != row['payload_digest']):
            raise ValueError('Transferred record fails its new-campaign native identity')
        validated.append({'key': key, 'record_sha256': record['record_sha256'], 'file_sha256': sha(copied)})
    if (stage/'sources.zip').read_bytes() != (parent/'sources.zip').read_bytes():
        raise ValueError('New campaign source ZIP is not byte-identical to its parent')
    cache._check_sources()
    origins = source_origins(plan)
    if controller()[0] != control or snapshot(parent, old) != request['parent']:
        raise ValueError('Controller or parent changed during staged verification')
    return {'schema': 'gabes-frozen-thermal-transfer-verification-v1', 'request_sha256': request['record_sha256'],
        'target_campaign_sha256': plan['record_sha256'], 'parent': request['parent'],
        'controller': control, 'source_identity': plan['source_identity'],
        'loaded_numerical_sources': origins, 'validated_records': validated,
        'imports_passed': True, 'new_solve_count': 0,
        'path_and_grid_evidence_imported': False, 'thermal_ensemble_converged': False,
        'physical_optical_prediction': False}


def verify_in_fresh_capsule(stage, control, raw):
    scratch = ROOT/'.git'/'grand-challenge-runs'
    scratch.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='thermal-import-check-', dir=scratch) as folder:
        c.extract_bundle(stage, folder)
        script = copy_controllers(folder, controller_files(control, raw, stage=stage))
        output = Path(folder)/'verification.json'
        env = os.environ.copy()
        env.pop('PYTHONPATH', None)
        subprocess.run([sys.executable, str(script), '--verify', str(stage),
            '--verification-output', str(output), '--controller-sha256', control['raw_sha256'],
            '--grid-sha256', control['dependencies'][0]['raw_sha256']],
            cwd=folder, env=env, check=True)
        return c.records.read_sealed(output)


def publish(stage, target):
    """Atomic no-replace directory rename; never overwrite a raced-in target."""
    if os.path.lexists(target):
        raise FileExistsError(target)
    if os.name == 'nt':
        os.rename(stage, target)  # Windows rename refuses any existing target.
    elif sys.platform.startswith('linux'):
        libc = ctypes.CDLL(None, use_errno=True)
        rename = getattr(libc, 'renameat2', None)
        if rename is None:
            raise RuntimeError('Atomic no-replace publication is unavailable on this platform')
        rename.argtypes = [ctypes.c_int, ctypes.c_char_p, ctypes.c_int, ctypes.c_char_p, ctypes.c_uint]
        rename.restype = ctypes.c_int
        if rename(-100, os.fsencode(stage), -100, os.fsencode(target), 1):
            number = ctypes.get_errno()
            raise OSError(number, os.strerror(number), str(target))
    else:
        raise RuntimeError('Atomic no-replace publication is supported on Windows and Linux only')


def execute(parent, target, powers, *, expected_controller_sha256=None, expected_grid_sha256=None, cache_snapshot=None):
    parent, target = destinations(parent, target)
    old = c.load_plan(parent)
    c.require_capsule(old)
    control, raw = controller()
    if expected_controller_sha256 is not None and control['raw_sha256'] != expected_controller_sha256:
        raise ValueError('Extension controller differs from captured launch bytes')
    if expected_grid_sha256 is not None and control['dependencies'][0]['raw_sha256'] != expected_grid_sha256:
        raise ValueError('Grid helper differs from captured launch bytes')
    captured_controllers = controller_files(control, raw)
    source_origins(old)
    origin = snapshot(parent, old)
    captured = present_snapshot(parent, origin, control) if cache_snapshot is None else cache_snapshot
    payload = {key: value for key, value in captured.items() if key != 'record_sha256'}
    if (c.records.sealed(payload) != captured or captured['parent'] != origin
            or captured['controller'] != control or captured['schema'] != 'gabes-thermal-cache-snapshot-v1'):
        raise ValueError('Captured launch snapshot differs from the parent or controller')
    plan = c.make_plan(c.ROOT, powers=tuple(powers), seeds=tuple(old['seeds']), max_steps=tuple(old['max_steps_s']))
    if max(plan['powers']) <= max(old['powers']):
        raise ValueError('Extension must declare a finer maximum power than its parent')
    compatible(old, plan)
    parent_cache = c.CampaignCache(parent/'cache', old)
    parent_index = job_index(old, parent_cache)
    new_index = job_index(plan, c.CampaignCache(parent/'cache', plan))
    imports, excluded = read_parent_records(parent, parent_index, new_index, parent_cache, captured)
    # Exact mathematical reuse: the ODE key binds model/source/path/spec/env,
    # not a campaign label or boundary-grid weight. Matching keys solve exactly
    # the same unweighted equations. Copy bytes, never refit or re-solve them.
    # Current-path gates, rates and grid evidence must be recomputed afterwards.
    plan.update(source_bundle_sha256=origin['source_bundle_sha256'],
        historical_numerical_results_reused=bool(imports),
        extension={'schema': 'gabes-frozen-thermal-extension-v1', 'parent': origin,
                   'controller': control, 'imported_raw_record_count': len(imports),
                   'path_and_grid_evidence_imported': False})
    with tempfile.TemporaryDirectory(prefix='.'+target.name+'-staging-', dir=target.parent) as folder:
        stage = Path(folder)/'campaign'
        stage.mkdir()
        (stage/'cache').mkdir()
        (stage/'controllers').mkdir()
        for identity, contents in captured_controllers:
            with (stage/identity['archive']).open('xb') as handle:
                handle.write(contents)
        (stage/'sources.zip').write_bytes((parent/'sources.zip').read_bytes())
        sealed_plan = c.write_record(stage/'plan.json', plan)
        for row in imports:
            copied = (parent/row['cache_file']).read_bytes()
            if sha(copied) != row['file_sha256']:
                raise ValueError('Parent cache record changed before copying')
            with (stage/row['cache_file']).open('xb') as handle:
                handle.write(copied)
        request = c.write_record(stage/'transfer-request.json', {'schema': SCHEMA,
            'parent_location': str(parent), 'parent': origin, 'controller': control,
            'target_campaign_sha256': sealed_plan['record_sha256'], 'requested_powers': plan['powers'],
            'parent_cache_snapshot': captured,
            'required_seeds': plan['seeds'], 'expected_native_keys': sorted(new_index),
            'imports': imports, 'excluded_parent_records': excluded,
            'missing_native_keys': sorted(set(new_index)-{row['key'] for row in imports}),
            'new_solve_count': 0, 'path_and_grid_evidence_imported': False})
        verified = verify_in_fresh_capsule(stage, control, raw)
        if (verified['request_sha256'] != request['record_sha256']
                or verified['target_campaign_sha256'] != sealed_plan['record_sha256']
                or verified['controller'] != control or verified['parent'] != origin
                or verified['imports_passed'] is not True or verified['new_solve_count'] != 0):
            raise ValueError('Fresh staged verification does not bind this exact transfer')
        for row in imports+excluded:
            if sha((parent/row['cache_file']).read_bytes()) != row['file_sha256']:
                raise ValueError('Parent cache changed before publication')
        for row in imports:
            if sha((stage/row['cache_file']).read_bytes()) != row['file_sha256']:
                raise ValueError('Copied cache changed after verification')
        if set(cache_files(stage/'cache')) != {row['key'] for row in imports}:
            raise ValueError('Staged cache inventory changed after verification')
        if (snapshot(parent, old) != origin or controller()[0] != control
                or c.records.read_sealed(stage/'plan.json') != sealed_plan
                or sha((stage/'sources.zip').read_bytes()) != origin['source_bundle_sha256']
                or any((stage/identity['archive']).read_bytes() != contents for identity, contents in captured_controllers)):
            raise ValueError('Declaration, archive or controller changed before publication')
        payload = {key: value for key, value in request.items() if key != 'record_sha256'}
        result = c.write_record(stage/'transfer.json', {**payload,
            'request_sha256': request['record_sha256'], 'verification': verified,
            'imported_raw_record_count': len(imports), 'transfer_passed': True,
            'thermal_ensemble_converged': False, 'physical_optical_prediction': False})
        (stage/'transfer-request.json').unlink()
        publish(stage, target)
    return result


def launch(parent, target, powers):
    parent, target = destinations(parent, target)
    control, raw = controller()
    captured_controllers = controller_files(control, raw)
    declaration = c.records.read_sealed(parent/'plan.json')
    captured = present_snapshot(parent, snapshot(parent, declaration), control)
    scratch = ROOT/'.git'/'grand-challenge-runs'
    scratch.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='thermal-extend-', dir=scratch) as folder:
        c.extract_bundle(parent, folder)
        script = copy_controllers(folder, captured_controllers)
        request = Path(folder)/'launch-snapshot.json'
        with request.open('x', encoding='utf-8') as handle:
            json.dump(captured, handle, ensure_ascii=False, allow_nan=False)
        env = os.environ.copy()
        env.pop('PYTHONPATH', None)
        subprocess.run([sys.executable, str(script), '--execute', str(parent),
            '--target', str(target), '--powers', *map(str, powers), '--controller-sha256', control['raw_sha256'],
            '--grid-sha256', control['dependencies'][0]['raw_sha256'], '--cache-snapshot', str(request)],
            cwd=folder, env=env, check=True)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    action = parser.add_mutually_exclusive_group(required=True)
    action.add_argument('--run', type=Path)
    action.add_argument('--execute', type=Path, help=argparse.SUPPRESS)
    action.add_argument('--verify', type=Path, help=argparse.SUPPRESS)
    parser.add_argument('--target', type=Path)
    parser.add_argument('--powers', nargs='+', type=int)
    parser.add_argument('--controller-sha256', help=argparse.SUPPRESS)
    parser.add_argument('--grid-sha256', help=argparse.SUPPRESS)
    parser.add_argument('--cache-snapshot', type=Path, help=argparse.SUPPRESS)
    parser.add_argument('--verification-output', type=Path, help=argparse.SUPPRESS)
    args = parser.parse_args(argv)
    if args.verify:
        if args.controller_sha256 is None or args.grid_sha256 is None or args.verification_output is None:
            parser.error('Captured verification requires controller hash and output')
        c.write_record(args.verification_output, verify_stage(args.verify,
                       expected_controller_sha256=args.controller_sha256, expected_grid_sha256=args.grid_sha256))
        return 0
    if args.target is None or args.powers is None:
        parser.error('--target NEW_DIRECTORY and --powers P0 P1 P2 [...] required')
    if args.run:
        launch(args.run, args.target, args.powers)
        return 0
    if args.controller_sha256 is None or args.grid_sha256 is None or args.cache_snapshot is None:
        parser.error('Captured --execute requires both controller hashes and the cache snapshot')
    result = execute(args.execute, args.target, args.powers,
                     expected_controller_sha256=args.controller_sha256, expected_grid_sha256=args.grid_sha256,
                     cache_snapshot=c.records.read_sealed(args.cache_snapshot))
    print(json.dumps({'target': str(args.target), 'imported_raw_record_count': result['imported_raw_record_count'],
                      'new_solve_count': 0, 'thermal_ensemble_converged': False}), flush=True)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
