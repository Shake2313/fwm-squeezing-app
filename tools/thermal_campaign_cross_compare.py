"""Read-only nested comparison across a direct frozen campaign extension.

Each input is first audited in its own marker-bound source capsule. Only the
new output is published; captured external controllers are embedded in it.
"""

import argparse
import base64
from contextlib import ExitStack
import hashlib
import os
from pathlib import Path
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from tools import thermal_campaign_compare as comparison
from tools import thermal_campaign_extend as extension

c, g = comparison.c, comparison.g
SCHEMA = 'gabes-frozen-thermal-cross-grid-comparison-v1'
SIDE_SCHEMA = 'gabes-frozen-thermal-cross-grid-side-audit-v1'


def sha(raw):
    return hashlib.sha256(raw).hexdigest()


def captured_controllers():
    captured = []
    for filename in (Path(__file__), Path(comparison.__file__), Path(g.__file__), Path(extension.__file__)):
        if filename.resolve().parent != Path(__file__).resolve().parent:
            raise ValueError('Cross comparison helpers must originate beside the controller')
        raw = filename.read_bytes()
        captured.append(({'name': filename.name, 'raw_sha256': sha(raw), 'raw_bytes': len(raw)}, raw))
    return captured


class ByteWatch:
    """Freeze only accessed files; unrelated append-only cache arrivals are fine."""

    def __init__(self):
        self.files = {}

    def read(self, filename):
        filename = Path(filename)
        if not extension.regular_cache_file(filename):
            raise ValueError('Audit input must be a regular independent file: '+str(filename))
        filename = filename.resolve()
        raw = filename.read_bytes()
        reference = {'path': str(filename), 'file_sha256': sha(raw), 'raw_bytes': len(raw)}
        previous = self.files.setdefault(str(filename), reference)
        if reference != previous:
            raise ValueError('Audited file bytes changed: '+str(filename))
        return raw

    def record(self, filename):
        # Parse the trusted byte snapshot itself. A separate native file read
        # could consume B between matching A reads and silently change a bound
        # launch request, parent proof or transfer record.
        result = g.sealed_snapshot(self.read(filename))
        self.read(filename)
        return result

    def include(self, references):
        for row in references:
            if set(row) != {'path', 'file_sha256', 'raw_bytes'}:
                raise ValueError('Invalid byte-watch reference')
            self.read(row['path'])
            if self.files[str(Path(row['path']).resolve())] != row:
                raise ValueError('Captured input byte reference differs')

    def finish(self):
        for filename in tuple(self.files):
            self.read(filename)
        return [self.files[name] for name in sorted(self.files)]


class ReadOnlyCache:
    def __init__(self, directory, plan, watch):
        self.directory, self.watch = Path(directory)/'cache', watch
        # The watch carries raw references to final publication; the shared
        # native guard also binds each returned seal to the bytes validated.
        self.cache = g.CacheByteGuard(c.CampaignCache(self.directory, plan))

    def key(self, model, path, spec):
        return self.cache.key(model, path, spec)

    def identity(self, model, path, spec):
        return self.cache.identity(model, path, spec)

    def get(self, model, path, spec, provider=None):
        if provider is not None:
            raise ValueError('Cross comparison forbids native providers')
        filename = self.directory/(self.key(model, path, spec)+'.json')
        self.watch.read(filename)
        packet, record = self.cache.get(model, path, spec, provider=None)
        self.watch.read(filename)
        if record['hit'] is not True:
            raise ValueError('Cross comparison requires existing immutable cache records')
        return packet, record

    def _check_sources(self):
        return self.cache._check_sources()


def capsule_files(root, plan, controllers, watch):
    root = Path(root).resolve()
    if watch.read(root/'.campaign-source-identity').decode('utf-8') != plan['record_sha256']:
        raise ValueError('Source capsule marker is not bound to this campaign')
    if c.source_paths(root) != plan['source_paths']:
        raise ValueError('Source capsule inventory differs')
    for row in plan['source_manifest']['files']:
        if sha(watch.read(root/row['path'])) != row['raw_sha256']:
            raise ValueError('Source capsule raw bytes differ')
    for identity in controllers:
        raw = watch.read(root/'tools'/identity['name'])
        if sha(raw) != identity['raw_sha256'] or len(raw) != identity['raw_bytes']:
            raise ValueError('Source capsule controller bytes differ')


def guard(plan, controllers, watch):
    c.require_capsule(plan)
    if [row[0] for row in captured_controllers()] != controllers:
        raise ValueError('Captured controller or dependency identity differs')
    capsule_files(c.ROOT, plan, controllers, watch)
    origins = extension.source_origins(plan)
    for filename in (Path(__file__), Path(comparison.__file__)):
        if not filename.resolve().is_relative_to(c.ROOT.resolve()):
            raise ValueError('Cross comparison controller loaded outside the capsule')
    return origins


def load_plan(directory, watch):
    declaration = watch.record(Path(directory)/'plan.json')
    if c.load_plan(directory) != declaration:
        raise ValueError('Campaign declaration changed while loading')
    if sha(watch.read(Path(directory)/'sources.zip')) != declaration['source_bundle_sha256']:
        raise ValueError('Frozen source ZIP differs from its declaration')
    return declaration


def read_grid(filename, directory, plan, watch):
    watch.read(filename)
    saved, reference, model, inflow = comparison.read_grid_report(filename, directory, plan)
    watch.read(filename)
    archived_controller(directory, saved['controller'], 'thermal_campaign_grid.py', watch)
    reference['campaign_sha256'] = plan['record_sha256']
    return saved, reference, model, inflow


def sealed(value):
    return c.records.sealed({key: item for key, item in value.items() if key != 'record_sha256'}) == value


def archived_controller(directory, identity, name, watch):
    token = identity['raw_sha256']
    if (identity['name'] != name or not isinstance(token, str) or len(token) != 64
            or any(char not in '0123456789abcdef' for char in token)
            or type(identity['raw_bytes']) is not int or identity['raw_bytes'] <= 0
            or identity['archive'] != 'controllers/'+Path(name).stem+'-'+token+'.py'):
        raise ValueError('Extension controller archive reference differs')
    raw = watch.read(Path(directory)/identity['archive'])
    if sha(raw) != token or len(raw) != identity['raw_bytes']:
        raise ValueError('Extension controller archive bytes differ')


def extension_contract(parent, directory, old, plan, old_cache, cache, watch):
    """Revalidate the published transfer, allowing later appended native jobs."""
    extension.compatible(old, plan)
    if (old['record_sha256'] == plan['record_sha256'] or max(plan['powers']) <= max(old['powers'])
            or watch.read(parent/'sources.zip') != watch.read(directory/'sources.zip')):
        raise ValueError('Comparison requires a distinct direct extension with the identical source ZIP')
    transfer = watch.record(directory/'transfer.json')
    origin = {'campaign_sha256': old['record_sha256'],
              'plan_file_sha256': sha(watch.read(parent/'plan.json')),
              'source_bundle_sha256': old['source_bundle_sha256']}
    provenance = plan['extension']
    control = transfer['controller']
    if (provenance['schema'] != 'gabes-frozen-thermal-extension-v1'
            or provenance['parent'] != origin or provenance['controller'] != control
            or provenance['path_and_grid_evidence_imported'] is not False
            or transfer['schema'] != extension.SCHEMA or transfer['parent'] != origin
            or transfer['target_campaign_sha256'] != plan['record_sha256']
            or transfer['requested_powers'] != plan['powers'] or transfer['required_seeds'] != plan['seeds']
            or transfer['transfer_passed'] is not True or transfer['new_solve_count'] != 0
            or type(transfer['new_solve_count']) is not int
            or any(transfer[name] is not False for name in
                   ('path_and_grid_evidence_imported', 'thermal_ensemble_converged', 'physical_optical_prediction'))):
        raise ValueError('Direct parent-to-child extension provenance differs')
    archived_controller(directory, control, 'thermal_campaign_extend.py', watch)
    if len(control['dependencies']) != 1:
        raise ValueError('Extension controller dependencies differ')
    archived_controller(directory, control['dependencies'][0], 'thermal_campaign_grid.py', watch)
    extra = {'record_sha256', 'request_sha256', 'verification', 'imported_raw_record_count',
             'transfer_passed', 'thermal_ensemble_converged', 'physical_optical_prediction'}
    request = {key: value for key, value in transfer.items() if key not in extra}
    if c.records.sealed(request)['record_sha256'] != transfer['request_sha256']:
        raise ValueError('Extension transfer request seal differs')
    verified = transfer['verification']
    if (not sealed(verified) or verified['schema'] != 'gabes-frozen-thermal-transfer-verification-v1'
            or verified['request_sha256'] != transfer['request_sha256'] or verified['parent'] != origin
            or verified['target_campaign_sha256'] != plan['record_sha256']
            or verified['controller'] != control or verified['source_identity'] != plan['source_identity']
            or verified['imports_passed'] is not True or type(verified['new_solve_count']) is not int
            or verified['new_solve_count'] != 0
            or any(verified[name] is not False for name in
                   ('path_and_grid_evidence_imported', 'thermal_ensemble_converged', 'physical_optical_prediction'))):
        raise ValueError('Extension transfer verification provenance differs')
    # The original verifier and this interpreter import the same numerical
    # module inventory from byte-identical source bundles, not the live tree.
    if verified['loaded_numerical_sources'] != extension.source_origins(plan):
        raise ValueError('Extension loaded numerical source provenance differs')
    captured = transfer['parent_cache_snapshot']
    if (not sealed(captured) or captured['schema'] != 'gabes-thermal-cache-snapshot-v1'
            or captured['parent'] != origin or captured['controller'] != control):
        raise ValueError('Extension parent snapshot provenance differs')
    old_index, index = extension.job_index(old, old_cache), extension.job_index(plan, cache)
    snapshots = {row['key']: row for row in captured['records']}
    imports = {row['key']: row for row in transfer['imports']}
    excluded = {row['key']: row for row in transfer['excluded_parent_records']}
    if (len(snapshots) != len(captured['records']) or len(imports) != len(transfer['imports'])
            or len(excluded) != len(transfer['excluded_parent_records']) or imports.keys() & excluded.keys()
            or imports.keys() | excluded.keys() != snapshots.keys()
            or transfer['expected_native_keys'] != sorted(index)
            or transfer['missing_native_keys'] != sorted(set(index)-imports.keys())
            or transfer['imported_raw_record_count'] != len(imports)
            or provenance['imported_raw_record_count'] != len(imports)
            or plan['historical_numerical_results_reused'] is not bool(imports)):
        raise ValueError('Extension transfer coverage differs from its declared native jobs')
    validated = []
    for key, row in {**imports, **excluded}.items():
        if (key not in old_index or (key in imports) != (key in index)
                or row['parent_selection'] != old_index[key][3] or row['cache_file'] != 'cache/'+key+'.json'):
            raise ValueError('Extension import/exclusion is not bound to its declared path')
        raw = watch.read(parent/row['cache_file'])
        if (sha(raw) != row['file_sha256'] or len(raw) != row['raw_bytes']
                or snapshots[key] != {'key': key, 'file_sha256': sha(raw), 'raw_bytes': len(raw)}):
            raise ValueError('Extension original record bytes differ from the transfer snapshot')
        model, path, spec, _ = old_index[key]
        _, reference = old_cache.get(model, path, spec)
        if comparison._reference(reference) != comparison._reference(row):
            raise ValueError('Extension parent native record reference differs')
        if key in imports:
            if row['target_selection'] != index[key][3] or watch.read(directory/row['cache_file']) != raw:
                raise ValueError('Extension imported record bytes or target binding differ')
            model, path, spec, _ = index[key]
            _, reference = cache.get(model, path, spec)
            if comparison._reference(reference) != comparison._reference(row):
                raise ValueError('Extension child native record reference differs')
            validated.append({'key': key, 'record_sha256': row['record_sha256'], 'file_sha256': sha(raw)})
    if verified['validated_records'] != validated:
        raise ValueError('Extension validated-record inventory differs')
    # parent_location records where the transfer happened. Portable lineage is
    # the exact sealed/raw parent plan and source ZIP, never its historical path.
    return {'record_sha256': transfer['record_sha256'], 'file_sha256': sha(watch.read(directory/'transfer.json')),
            'parent': origin, 'target_campaign_sha256': plan['record_sha256'],
            'historical_parent_location': transfer['parent_location'],
            'imported_record_keys': sorted(imports), 'validated': True}


def nested_reuse_audit(model, old_inflow, inflow, coarse, fine, numerical, old_cache, cache):
    """Read each side explicitly: some nested paths have identical source labels."""
    pairs = comparison.nested_indices(coarse['power'], fine['power'])
    if (coarse['seed'] != fine['seed'] or len(old_inflow.rate_s_inverse) != len(pairs)
            or len(inflow.rate_s_inverse) != 2*len(pairs)):
        raise ValueError('Cross reuse requires the same seed and doubled complete grids')
    rows, mapped = [], set()
    for old, new in pairs:
        before, after = old_inflow.path(old), inflow.path(new)
        face, offset = divmod(old, 2**coarse['power'])
        if (new in mapped or int(old_inflow.face_index[old]) != face or int(inflow.face_index[new]) != face
                or any(not comparison._bitwise_equal(getattr(before, name), getattr(after, name))
                       for name in ('entry_position_m', 'velocity_m_s', 'residence_time_s'))
                or inflow.rate_s_inverse[new] != .5*old_inflow.rate_s_inverse[old]):
            raise ValueError('Nested map does not preserve exact physical paths and half rates')
        mapped.add(new)
        lo, hi = coarse['path_ledgers'][old], fine['path_ledgers'][new]
        records = []
        for level, spec in enumerate(numerical.primary+numerical.reference):
            one, ref_one = comparison._checked_read(model, before, spec, old_cache, lo, level)
            two, ref_two = comparison._checked_read(model, after, spec, cache, hi, level)
            if comparison._reference(ref_one) != comparison._reference(ref_two):
                raise ValueError('Cross-campaign shared record identities differ')
            if any(not comparison._bitwise_equal(one[name], two[name]) for name in c.r.METRICS):
                raise ValueError('Cross-campaign raw eight-metric packets are not bitwise identical')
            records.append({'level': level, **comparison._reference(ref_one), 'all_eight_metrics_bitwise_equal': True})
        if (before.source == after.source) != (lo['target_digest'] == hi['target_digest']):
            raise ValueError('Nested target digests do not reflect current source rebinding')
        rows.append({'coarse_index': old, 'fine_index': new, 'face': face, 'offset': offset,
            'physical_path': c.path_identity(before), 'coarse_rate_s_inverse': float(old_inflow.rate_s_inverse[old]),
            'fine_rate_s_inverse': float(inflow.rate_s_inverse[new]), 'coarse_source': before.source,
            'fine_source': after.source, 'coarse_target_digest': lo['target_digest'], 'fine_target_digest': hi['target_digest'],
            'cache_records': records})
    return {'passed': True, 'paired_path_count': len(rows), 'shared_unique_jobs': 5*len(rows),
        'metric_names': c.r.METRICS, 'pairs': rows,
        'new_fine_indices': sorted(set(range(len(inflow.rate_s_inverse)))-mapped), 'new_solve_count': 0,
        'historical_hit_flags_used_as_execution_evidence': False,
        'scope': 'Both distinct stores reread; identical unweighted equations and five raw calculations per pair'}


def read_request(filename, expected_hash, watch):
    if sha(watch.read(filename)) != expected_hash:
        raise ValueError('Cross comparison launch request bytes differ')
    request = watch.record(filename)
    if request['schema'] != SCHEMA+'-request':
        raise ValueError('Wrong cross comparison launch schema')
    watch.include(request['input_bytes'])
    return request


def audit_parent(request_file, expected_hash, output):
    watch = ByteWatch()
    request = read_request(request_file, expected_hash, watch)
    directory = Path(request['parent_run'])
    plan = load_plan(directory, watch)
    origins = guard(plan, request['controllers'], watch)
    saved, reference, model, inflow = read_grid(request['coarse'], directory, plan, watch)
    cache = ReadOnlyCache(directory, plan, watch)
    row, direct = comparison.reread_grid(saved, model, inflow, c.pilot.step_plan(plan['max_steps_s']), cache, plan)
    cache._check_sources()
    guard(plan, request['controllers'], watch)
    return c.write_record(output, {'schema': SIDE_SCHEMA, 'request_sha256': request['record_sha256'],
        'campaign_sha256': plan['record_sha256'], 'capsule_root': str(c.ROOT.resolve()),
        'input_report': reference, 'row': row, 'aggregation_audit': direct,
        'loaded_numerical_sources': origins, 'controllers': request['controllers'],
        'input_bytes': watch.finish(), 'audit_passed': True, 'new_solve_count': 0})


def execute(request_file, expected_hash, parent_audit_file, parent_audit_hash):
    watch = ByteWatch()
    request = read_request(request_file, expected_hash, watch)
    parent, directory = Path(request['parent_run']), Path(request['run'])
    output = new_output(request['output'], parent, directory)
    plan, old = load_plan(directory, watch), load_plan(parent, watch)
    origins = guard(plan, request['controllers'], watch)
    if sha(watch.read(parent_audit_file)) != parent_audit_hash:
        raise ValueError('Parent capsule audit bytes differ')
    proof = c.r.decode(watch.record(parent_audit_file))
    if (proof['schema'] != SIDE_SCHEMA or proof['request_sha256'] != request['record_sha256']
            or proof['campaign_sha256'] != old['record_sha256'] or proof['controllers'] != request['controllers']
            or proof['audit_passed'] is not True or type(proof['new_solve_count']) is not int
            or proof['new_solve_count'] != 0 or Path(proof['capsule_root']).resolve() == c.ROOT.resolve()):
        raise ValueError('Parent report lacks a separate campaign-bound capsule audit')
    watch.include(proof['input_bytes'])
    capsule_files(proof['capsule_root'], old, request['controllers'], watch)
    if proof['loaded_numerical_sources'] != origins:
        raise ValueError('Parent capsule loaded source provenance differs')
    saved_old, ref_old, model, old_inflow = read_grid(request['coarse'], parent, old, watch)
    saved, reference, _, inflow = read_grid(request['fine'], directory, plan, watch)
    comparison.nested_indices(saved_old['selection'][0], saved['selection'][0])
    if saved_old['selection'][1] != saved['selection'][1] or proof['input_report'] != ref_old:
        raise ValueError('Cross comparison requires bound reports with the same seed')
    old_cache, cache = ReadOnlyCache(parent, old, watch), ReadOnlyCache(directory, plan, watch)
    transfer = extension_contract(parent, directory, old, plan, old_cache, cache, watch)
    numerical = c.pilot.step_plan(plan['max_steps_s'])
    coarse, coarse_audit = comparison.reread_grid(saved_old, model, old_inflow, numerical, old_cache, old)
    fine, fine_audit = comparison.reread_grid(saved, model, inflow, numerical, cache, plan)
    if (c.r.digest(comparison._row_semantics(proof['row'])) != c.r.digest(comparison._row_semantics(coarse))
            or c.r.digest(comparison._audit_semantics(proof['aggregation_audit'])) !=
               c.r.digest(comparison._audit_semantics(coarse_audit))):
        raise ValueError('Parent capsule result differs from the actual current evidence')
    # Exact mathematical reuse: identical model/source/path/spec/environment
    # defines the same unweighted ODE problem across campaign labels. The five
    # original records must be byte-identical; no repeat solve is justified.
    # Current-path source/digests, half rates, and all path gates are rebuilt.
    reuse = nested_reuse_audit(model, old_inflow, inflow, coarse, fine, numerical, old_cache, cache)
    imported = set(transfer['imported_record_keys'])
    for pair in reuse['pairs']:
        for row in pair['cache_records']:
            key = row['key']
            if key not in imported or watch.read(parent/'cache'/(key+'.json')) != watch.read(directory/'cache'/(key+'.json')):
                raise ValueError('Shared raw computation is not an exact declared imported record')
    split = comparison.nested_sum_audit(model, inflow, coarse, fine, numerical, cache, reuse)
    errors = c.r.stream_errors(fine['candidate']['spectra'], coarse['candidate']['spectra'], coarse['comparison_scales'])
    diagnostic = {'errors': errors, 'budget': plan['ensemble_budget'],
        'passed': all(value <= plan['ensemble_budget'] for value in errors.values()),
        'error_definition': c.r.ERROR_DEFINITION, 'reference_scales': 'coarse grid comparison_scales',
        'scope': 'One refinement edge; not a campaign convergence certificate'}
    declarations = []
    for current, row in ((old, coarse), (plan, fine)):
        declarations.append({'campaign_sha256': current['record_sha256'], 'required_powers': current['powers'],
            'required_seeds': current['seeds'], 'missing_grids': [[p, s] for p in current['powers'] for s in current['seeds']
                if [p, s] != [row['power'], row['seed']]],
            'convergence_gate': c.r.convergence_gate([row], powers=current['powers'], seeds=current['seeds'])})
    old_cache._check_sources()
    cache._check_sources()
    guard(plan, request['controllers'], watch)
    capsule_files(proof['capsule_root'], old, request['controllers'], watch)
    captured = captured_controllers()
    return c.write_record(output, {'schema': SCHEMA, 'input_reports': [ref_old, reference],
        'source_identity': plan['source_identity'], 'source_bundle_sha256': plan['source_bundle_sha256'],
        'extension_transfer': transfer, 'controllers': request['controllers'],
        'controller_archives': [dict(identity, encoding='base64', raw_base64=base64.b64encode(raw).decode('ascii'))
                                for identity, raw in captured],
        'capsule_audits': [{'campaign_sha256': old['record_sha256'], 'record_sha256': proof['record_sha256'],
                            'capsule_root': proof['capsule_root'], 'loaded_numerical_sources': proof['loaded_numerical_sources']},
                           {'campaign_sha256': plan['record_sha256'], 'capsule_root': str(c.ROOT.resolve()),
                            'loaded_numerical_sources': origins}],
        'rows': [coarse, fine], 'fresh_grid_aggregation_audits': [coarse_audit, fine_audit],
        'nested_reuse': reuse, 'nested_sum': split, 'refinement_diagnostic': diagnostic,
        'campaign_declarations': declarations, 'input_byte_checks': watch.finish(),
        'audit_passed': bool(reuse['passed'] and split['passed']), 'new_solve_count': 0,
        'thermal_ensemble_converged': False, 'physical_optical_prediction': False,
        'scope': 'Direct extension and exact nested atomic stream reuse; no whole-campaign or optical certification'})


def new_output(output, parent, directory):
    output = g._new_output(output, directory)
    g._new_output(output, parent)
    for root in (parent, directory):
        if output.is_relative_to(Path(root)/'controllers'):
            raise ValueError('Comparison output must not alter campaign controller archives')
    return output


def launch(parent, directory, coarse, fine, output):
    parent, directory = Path(parent).resolve(), Path(directory).resolve()
    output = new_output(output, parent, directory)
    if parent == directory:
        raise ValueError('Cross comparison requires distinct parent and child campaigns')
    coarse, fine = Path(coarse).resolve(), Path(fine).resolve()
    captured = captured_controllers()
    watch = ByteWatch()
    for filename in (parent/'plan.json', parent/'sources.zip', directory/'plan.json', directory/'sources.zip',
                     directory/'transfer.json', coarse, fine):
        watch.read(filename)
    for identity, raw in captured:
        if watch.read(Path(__file__).parent/identity['name']) != raw:
            raise ValueError('Controller changed during launch capture')
    scratch = ROOT/'.git/grand-challenge-runs'
    scratch.mkdir(parents=True, exist_ok=True)
    with ExitStack() as stack:
        roots = [Path(stack.enter_context(tempfile.TemporaryDirectory(prefix='thermal-cross-', dir=scratch)))
                 for _ in range(2)]
        for campaign, root in zip((parent, directory), roots):
            c.extract_bundle(campaign, root)
            (root/'tools').mkdir()
            for identity, raw in captured:
                with (root/'tools'/identity['name']).open('xb') as handle:
                    handle.write(raw)
        request = roots[0]/'request.json'
        c.write_record(request, {'schema': SCHEMA+'-request', 'parent_run': str(parent), 'run': str(directory),
            'coarse': str(coarse), 'fine': str(fine), 'output': str(output),
            'controllers': [item[0] for item in captured], 'input_bytes': watch.finish()})
        request_hash = sha(request.read_bytes())
        proof = roots[0]/'parent-audit.json'
        env = os.environ.copy()
        env.pop('PYTHONPATH', None)
        subprocess.run([sys.executable, str(roots[0]/'tools'/Path(__file__).name), '--audit-parent',
            '--request', str(request), '--request-sha256', request_hash, '--output', str(proof)],
            cwd=roots[0], env=env, check=True)
        completed = subprocess.run([sys.executable, str(roots[1]/'tools'/Path(__file__).name), '--execute',
            '--request', str(request), '--request-sha256', request_hash,
            '--parent-audit', str(proof), '--parent-audit-sha256', sha(proof.read_bytes())],
            cwd=roots[1], env=env, check=False)
        return completed.returncode


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--parent-run', type=Path)
    parser.add_argument('--run', type=Path)
    parser.add_argument('--coarse', type=Path)
    parser.add_argument('--fine', type=Path)
    parser.add_argument('--output', type=Path)
    internal = parser.add_mutually_exclusive_group()
    internal.add_argument('--audit-parent', action='store_true', help=argparse.SUPPRESS)
    internal.add_argument('--execute', action='store_true', help=argparse.SUPPRESS)
    parser.add_argument('--request', type=Path, help=argparse.SUPPRESS)
    parser.add_argument('--request-sha256', help=argparse.SUPPRESS)
    parser.add_argument('--parent-audit', type=Path, help=argparse.SUPPRESS)
    parser.add_argument('--parent-audit-sha256', help=argparse.SUPPRESS)
    args = parser.parse_args(argv)
    if args.audit_parent or args.execute:
        if not args.request or not args.request_sha256:
            parser.error('Captured execution requires the launch request and its raw hash')
        if args.audit_parent:
            if not args.output:
                parser.error('Parent capsule audit requires an exclusive output')
            audit_parent(args.request, args.request_sha256, args.output)
        else:
            if not args.parent_audit or not args.parent_audit_sha256:
                parser.error('Cross comparison requires the captured parent capsule audit')
            result = execute(args.request, args.request_sha256, args.parent_audit, args.parent_audit_sha256)
            return 0 if result['audit_passed'] else 1
    else:
        if any(value is None for value in (args.parent_run, args.run, args.coarse, args.fine, args.output)):
            parser.error('--parent-run PARENT --run CHILD --coarse GRID --fine GRID --output NEW_JSON required')
        return launch(args.parent_run, args.run, args.coarse, args.fine, args.output)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
