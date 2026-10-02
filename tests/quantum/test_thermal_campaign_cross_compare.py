"""Cross-campaign evidence contracts; synthetic packets, never native solves."""

from copy import deepcopy
import json
from pathlib import Path
import shutil
import subprocess

import numpy as np
import pytest

from tools import thermal_campaign_cross_compare as x
from test_thermal_campaign_grid import synthetic_packet
from test_thermal_campaign import native_numbers

c, g, e, t = x.c, x.g, x.extension, x.comparison
REAL_ROOT, REAL_CACHE, REAL_LOAD = c.ROOT, c.CampaignCache, c.load_plan
REAL_VALIDATE, REAL_INDEX = c.validate_native, e.job_index
ORIGINS = [{'module': 'gabes.fixture', 'source_path': 'gabes/fixture.py',
            'raw_sha256': x.sha(b'fixture = 1\n')}]


def test_request_uses_trusted_snapshot_when_a_separate_read_would_consume_another_record(
        tmp_path, monkeypatch):
    filename = tmp_path/'request.json'
    saved = c.write_record(filename, {'schema': x.SCHEMA+'-request',
        'input_bytes': [], 'output': 'declared-output.json'})
    raw = filename.read_bytes()
    alternate = c.records.sealed({'schema': x.SCHEMA+'-request',
        'input_bytes': [], 'output': 'alternate-output.json'})
    original_read = c.records.read_sealed
    alternate_reads = []

    def transient_record(path):
        alternate_reads.append(path)
        filename.write_text(json.dumps(alternate), encoding='utf-8')
        try:
            return original_read(path)
        finally:
            filename.write_bytes(raw)

    monkeypatch.setattr(c.records, 'read_sealed', transient_record)
    request = x.read_request(filename, x.sha(raw), x.ByteWatch())
    assert request['output'] == 'declared-output.json'
    assert request['record_sha256'] == saved['record_sha256']
    assert request == saved and filename.read_bytes() == raw
    assert not alternate_reads


@pytest.mark.parametrize('damage', ['seal', 'duplicate'])
def test_trusted_request_hash_does_not_accept_an_invalid_record(tmp_path, damage):
    filename = tmp_path/'request.json'
    record = c.records.sealed({'schema': x.SCHEMA+'-request',
        'input_bytes': [], 'output': 'declared-output.json'})
    if damage == 'seal':
        record['record_sha256'] = '0'*64
        raw = json.dumps(record).encode('utf-8')
    else:
        raw = (json.dumps(record)[:-1]+', "output": "alternate-output.json"}').encode('utf-8')
    filename.write_bytes(raw)
    with pytest.raises(ValueError, match='content hash mismatch|Duplicate audit input key'):
        x.read_request(filename, x.sha(raw), x.ByteWatch())


def test_cross_cache_rejects_a_transient_native_record_even_when_bytes_revert(
        tmp_path, monkeypatch):
    directory = tmp_path/'cache'
    directory.mkdir()
    key = 'a'*64
    filename = directory/(key+'.json')
    saved = c.write_record(filename, {'value': 1})
    raw = filename.read_bytes()
    alternate = c.records.sealed({'value': 2})

    class NativeCache:
        def __init__(self, directory, plan):
            self.directory = directory

        def key(self, *args):
            return key

        def get(self, *args, provider=None):
            assert provider is None
            filename.write_text(json.dumps(alternate), encoding='utf-8')
            try:
                native = c.records.read_sealed(filename)
            finally:
                filename.write_bytes(raw)
            return {}, {'key': key, 'hit': True, 'record_sha256': native['record_sha256']}

    monkeypatch.setattr(c, 'CampaignCache', NativeCache)
    cache = x.ReadOnlyCache(tmp_path, {}, x.ByteWatch())
    with pytest.raises(ValueError, match='differs from consumed cache byte snapshot'):
        cache.get(None, None, None)
    assert filename.read_bytes() == raw
    assert c.records.read_sealed(filename) == saved


def reseal(path, change):
    record = c.records.read_sealed(path)
    del record['record_sha256']
    change(record)
    path.write_text(json.dumps(c.records.sealed(record), ensure_ascii=False), encoding='utf-8')


def hashes(root):
    return {p.relative_to(root).as_posix(): x.sha(p.read_bytes()) for p in root.rglob('*') if p.is_file()}


def validate_synthetic_native(packet, path, spec):
    # Only the in-process fixture removes its explicit synthetic marker. Real
    # subprocesses retain the unmodified native validator and cannot accept it.
    numbers = dict(packet['numerics'])
    assert numbers.pop('synthetic_test_fixture') is True
    REAL_VALIDATE({'numerics': numbers}, path, spec)


def install_fixture_cache(patch, source):
    patch.setattr(c, 'CampaignCache', lambda directory, plan: REAL_CACHE(directory, plan, root=source))
    patch.setattr(c, 'load_plan', lambda directory: REAL_LOAD(directory, root=source))
    patch.setattr(c, 'validate_native', validate_synthetic_native)
    patch.setattr(c, 'native_provider', lambda *args: pytest.fail('cross comparison must never solve'))


def fixture_records(directory, plan, power):
    cache = c.CampaignCache(directory/'cache', plan)
    model = c.r.default_model()
    inflow = model.inflow(power, 211)
    numerical = c.pilot.step_plan(plan['max_steps_s'])
    for index in range(len(inflow.rate_s_inverse)):
        path = inflow.path(index)
        for spec in numerical.primary+numerical.reference:
            key = cache.key(model, path, spec)
            file = directory/'cache'/(key+'.json')
            if file.exists():
                continue
            packet = synthetic_packet(model, path, spec)
            packet['numerics'] = dict(native_numbers(path, spec), synthetic_test_fixture=True)
            packet = c.r._validated_packet(model, path, packet)
            c.write_record(file, {'schema': 'gabes-native-thermal-packet-v1',
                'identity': c.r.encode(cache.identity(model, path, spec)), 'key': key,
                'original_path_source': path.source, 'packet': packet,
                'source_before': plan['source_manifest'], 'source_after': plan['source_manifest'],
                'source_stable_during_run': True, 'synthetic_contract_fixture': True})


def transfer_record(parent, child, old, plan, indexes):
    origin = e.snapshot(parent, old)
    control = plan['extension']['controller']
    captured = e.present_snapshot(parent, origin, control)
    imported, excluded = e.read_parent_records(parent, indexes[old['record_sha256']],
        indexes[plan['record_sha256']], c.CampaignCache(parent/'cache', old), captured)
    request = c.records.sealed({'schema': e.SCHEMA, 'parent_location': str(parent), 'parent': origin,
        'controller': control, 'target_campaign_sha256': plan['record_sha256'],
        'requested_powers': plan['powers'], 'required_seeds': plan['seeds'], 'parent_cache_snapshot': captured,
        'expected_native_keys': sorted(indexes[plan['record_sha256']]), 'imports': imported,
        'excluded_parent_records': excluded, 'missing_native_keys': sorted(set(indexes[plan['record_sha256']])-
            {row['key'] for row in imported}), 'new_solve_count': 0, 'path_and_grid_evidence_imported': False})
    verified = c.records.sealed({'schema': 'gabes-frozen-thermal-transfer-verification-v1',
        'request_sha256': request['record_sha256'], 'target_campaign_sha256': plan['record_sha256'],
        'parent': origin, 'controller': control, 'source_identity': plan['source_identity'],
        'loaded_numerical_sources': ORIGINS,
        'validated_records': [{key: row[key] for key in ('key', 'record_sha256', 'file_sha256')} for row in imported],
        'imports_passed': True, 'new_solve_count': 0, 'path_and_grid_evidence_imported': False,
        'thermal_ensemble_converged': False, 'physical_optical_prediction': False})
    return {**{key: value for key, value in request.items() if key != 'record_sha256'},
        'request_sha256': request['record_sha256'], 'verification': verified,
        'imported_raw_record_count': len(imported), 'transfer_passed': True,
        'thermal_ensemble_converged': False, 'physical_optical_prediction': False}


@pytest.fixture(scope='module')
def base(tmp_path_factory):
    root = tmp_path_factory.mktemp('cross-base')
    source, parent, child = root/'source', root/'parent', root/'child'
    (source/'gabes').mkdir(parents=True)
    (source/'gabes/fixture.py').write_bytes(b'fixture = 1\n')
    c.create_campaign(parent, root=source)
    old = c.records.read_sealed(parent/'plan.json')
    child.mkdir()
    (child/'cache').mkdir()
    shutil.copyfile(parent/'sources.zip', child/'sources.zip')
    control, raw = e.controller()
    t.archive_controllers(child, e.controller_files(control, raw))
    plan = c.make_plan(source, powers=(1, 2, 3))
    plan.update(source_bundle_sha256=old['source_bundle_sha256'], historical_numerical_results_reused=True,
        extension={'schema': 'gabes-frozen-thermal-extension-v1', 'parent': e.snapshot(parent, old),
                   'controller': control, 'imported_raw_record_count': 30, 'path_and_grid_evidence_imported': False})
    plan = c.write_record(child/'plan.json', plan)
    with pytest.MonkeyPatch.context() as patch:
        install_fixture_cache(patch, source)
        patch.setattr(c, 'require_capsule', lambda *args: None)
        patch.setattr(g, '_loaded_from_capsule', lambda: None)
        fixture_records(parent, old, 0)
        for file in (parent/'cache').iterdir():
            shutil.copyfile(file, child/'cache'/file.name)
        fixture_records(child, plan, 1)
        g.execute(parent, 0, 211, parent/'grid.json')
        g.execute(child, 1, 211, child/'grid.json')
        indexes = {value['record_sha256']: REAL_INDEX(value, c.CampaignCache(directory/'cache', value))
                   for directory, value in ((parent, old), (child, plan))}
        c.write_record(child/'transfer.json', transfer_record(parent, child, old, plan, indexes))
    return root, indexes


@pytest.fixture
def campaign(base, tmp_path, monkeypatch):
    base_root, indexes = base
    parent, child = tmp_path/'parent', tmp_path/'child'
    shutil.copytree(base_root/'parent', parent)
    shutil.copytree(base_root/'child', child)
    source = base_root/'source'
    install_fixture_cache(monkeypatch, source)
    monkeypatch.setattr(e, 'source_origins', lambda plan: deepcopy(ORIGINS))
    # The declarations were completely enumerated by the real helper once.
    # Cache that pure index to keep repeated negative tests bounded.
    monkeypatch.setattr(e, 'job_index', lambda plan, cache: indexes[plan['record_sha256']])
    old, plan = c.load_plan(parent), c.load_plan(child)
    (child/'transfer.json').unlink()
    c.write_record(child/'transfer.json', transfer_record(parent, child, old, plan, indexes))
    captured = x.captured_controllers()
    roots = [tmp_path/'parent-capsule', tmp_path/'child-capsule']
    for directory, root in zip((parent, child), roots):
        root.mkdir()
        c.extract_bundle(directory, root)
        (root/'tools').mkdir()
        for identity, raw in captured:
            (root/'tools'/identity['name']).write_bytes(raw)
    def in_process_guard(value, identities, watch):
        # Numeric imports remain live only in these synthetic unit tests. Keep
        # real marker, inventory, raw source and external-controller checks.
        c.require_capsule(value)
        assert [row[0] for row in x.captured_controllers()] == identities
        x.capsule_files(c.ROOT, value, identities, watch)
        return deepcopy(ORIGINS)
    monkeypatch.setattr(x, 'guard', in_process_guard)
    watch = x.ByteWatch()
    for file in (parent/'plan.json', child/'plan.json', parent/'sources.zip', child/'sources.zip',
                 parent/'grid.json', child/'grid.json', child/'transfer.json'):
        watch.read(file)
    request_file, output, proof_file = tmp_path/'request.json', tmp_path/'comparison.json', tmp_path/'proof.json'
    request = {'schema': x.SCHEMA+'-request', 'parent_run': str(parent), 'run': str(child),
        'coarse': str(parent/'grid.json'), 'fine': str(child/'grid.json'), 'output': str(output),
        'controllers': [row[0] for row in captured], 'input_bytes': watch.finish()}
    c.write_record(request_file, request)
    def audit():
        monkeypatch.setattr(c, 'ROOT', roots[0])
        result = x.audit_parent(request_file, x.sha(request_file.read_bytes()), proof_file)
        monkeypatch.setattr(c, 'ROOT', roots[1])
        return result
    def run():
        audit()
        return x.execute(request_file, x.sha(request_file.read_bytes()), proof_file, x.sha(proof_file.read_bytes()))
    return {'parent': parent, 'child': child, 'old': old, 'plan': plan, 'roots': roots, 'output': output,
            'request': request_file, 'proof': proof_file, 'audit': audit, 'run': run, 'indexes': indexes}


def refresh_request(fixture):
    def update(request):
        for row in request['input_bytes']:
            raw = Path(row['path']).read_bytes()
            row.update(file_sha256=x.sha(raw), raw_bytes=len(raw))
    reseal(fixture['request'], update)


def execute_after_proof(fixture):
    return x.execute(fixture['request'], x.sha(fixture['request'].read_bytes()),
                     fixture['proof'], x.sha(fixture['proof'].read_bytes()))


def test_cross_comparison_revalidates_separate_stores_and_preserves_campaigns(campaign):
    before = [hashes(campaign[name]) for name in ('parent', 'child')]
    report = c.r.decode(campaign['run']())
    assert report['schema'] == 'gabes-frozen-thermal-cross-grid-comparison-v1'
    assert report['audit_passed'] and report['nested_sum']['passed'] and report['new_solve_count'] == 0
    assert not report['thermal_ensemble_converged'] and not report['physical_optical_prediction']
    assert [row['campaign_sha256'] for row in report['input_reports']] == [campaign['old']['record_sha256'], campaign['plan']['record_sha256']]
    assert [row['selection'] for row in report['input_reports']] == [[0, 211], [1, 211]]
    for reference in report['input_reports']:
        raw = Path(reference['path']).read_bytes()
        assert reference['file_sha256'] == x.sha(raw)
        assert reference['record_sha256'] == c.records.read_sealed(reference['path'])['record_sha256']
    assert report['nested_reuse']['shared_unique_jobs'] == 30
    assert len(report['nested_reuse']['new_fine_indices']) == 6
    assert max(report['nested_sum']['errors'].values()) < 1e-13
    assert report['refinement_diagnostic']['errors'] == c.r.stream_errors(
        report['rows'][1]['candidate']['spectra'], report['rows'][0]['candidate']['spectra'], report['rows'][0]['comparison_scales'])
    assert all(not row['convergence_gate']['passed'] and row['missing_grids'] for row in report['campaign_declarations'])
    assert report['capsule_audits'][0]['capsule_root'] != report['capsule_audits'][1]['capsule_root']
    for archive in report['controller_archives']:
        raw = x.base64.b64decode(archive['raw_base64'])
        assert x.sha(raw) == archive['raw_sha256'] and len(raw) == archive['raw_bytes']
    assert [hashes(campaign[name]) for name in ('parent', 'child')] == before


def test_unchanged_campaigns_can_move_without_rewriting_historical_lineage(campaign):
    parent, child = campaign['parent'], campaign['child']
    before = [hashes(parent), hashes(child)]
    relocated = campaign['request'].parent/'relocated'
    relocated.mkdir()
    new_parent, new_child = relocated/'parent', relocated/'child'
    # Owned temporary fixtures only; both old historical locations disappear.
    parent.rename(new_parent)
    child.rename(new_child)
    def move_request(request):
        request.update(parent_run=str(new_parent), run=str(new_child),
                       coarse=str(new_parent/'grid.json'), fine=str(new_child/'grid.json'))
        for row in request['input_bytes']:
            file = Path(row['path'])
            for old, new in ((parent, new_parent), (child, new_child)):
                if file.is_relative_to(old):
                    row['path'] = str(new/file.relative_to(old))
                    break
    reseal(campaign['request'], move_request)
    report = campaign['run']()
    assert report['audit_passed'] and not parent.exists() and not child.exists()
    assert report['extension_transfer']['historical_parent_location'] == str(parent)
    assert [hashes(new_parent), hashes(new_child)] == before


def test_cli_preserves_failed_split_evidence_and_returns_nonzero(campaign, monkeypatch):
    campaign['audit']()
    original = t._phase_terms
    def missing_number(packet):
        terms = original(packet)
        terms['poisson_number'] = np.zeros_like(terms['poisson_number'])
        return terms
    monkeypatch.setattr(t, '_phase_terms', missing_number)
    code = x.main(['--execute', '--request', str(campaign['request']),
        '--request-sha256', x.sha(campaign['request'].read_bytes()),
        '--parent-audit', str(campaign['proof']), '--parent-audit-sha256', x.sha(campaign['proof'].read_bytes())])
    report = c.records.read_sealed(campaign['output'])
    assert code == 1 and report['audit_passed'] is False
    assert report['nested_sum']['errors']['poisson_number'] > report['nested_sum']['budget']
    assert report['new_solve_count'] == 0 and report['thermal_ensemble_converged'] is False


@pytest.mark.parametrize('side', ['parent', 'child'])
@pytest.mark.parametrize('damage', ['missing', 'changed_bytes'])
def test_each_store_is_required_despite_identical_other_copy(campaign, side, damage):
    campaign['audit']()
    filename = next((campaign['parent']/'cache').glob('*.json')).name
    file = campaign[side]/'cache'/filename
    other = campaign['child' if side == 'parent' else 'parent']/'cache'/filename
    original = other.read_bytes()
    if damage == 'missing':
        file.unlink()
    else:
        file.write_bytes(file.read_bytes()+b'\n')
    with pytest.raises((ValueError, FileNotFoundError)):
        execute_after_proof(campaign)
    assert other.read_bytes() == original and not campaign['output'].exists()


@pytest.mark.parametrize('field', ['parent', 'controller', 'verification', 'coverage', 'sources'])
def test_resealed_transfer_provenance_cannot_substitute_for_contract(campaign, field):
    def change(record):
        if field == 'parent':
            record['parent']['campaign_sha256'] = '0'*64
        elif field == 'controller':
            record['controller']['raw_sha256'] = '0'*64
        elif field == 'verification':
            record['verification']['target_campaign_sha256'] = '0'*64
        elif field == 'coverage':
            record['missing_native_keys'].pop()
        else:
            record['verification']['loaded_numerical_sources'][0]['raw_sha256'] = '0'*64
    reseal(campaign['child']/'transfer.json', change)
    refresh_request(campaign)
    with pytest.raises(ValueError, match='provenance|archive|seal'):
        campaign['run']()
    assert not campaign['output'].exists()


@pytest.mark.parametrize('field', ['coverage', 'source_origins', 'validated_records'])
def test_transfer_semantics_are_checked_even_after_every_nested_seal_is_recomputed(campaign, field):
    def change(record):
        if field == 'coverage':
            record['missing_native_keys'].pop()
        elif field == 'source_origins':
            record['verification']['loaded_numerical_sources'][0]['raw_sha256'] = '0'*64
        else:
            record['verification']['validated_records'].pop()
        extra = {'request_sha256', 'verification', 'imported_raw_record_count', 'transfer_passed',
                 'thermal_ensemble_converged', 'physical_optical_prediction'}
        request = {key: value for key, value in record.items() if key not in extra}
        record['request_sha256'] = c.records.sealed(request)['record_sha256']
        verified = record['verification']
        del verified['record_sha256']
        verified['request_sha256'] = record['request_sha256']
        record['verification'] = c.records.sealed(verified)
    reseal(campaign['child']/'transfer.json', change)
    refresh_request(campaign)
    with pytest.raises(ValueError, match='coverage differs|source provenance differs|validated-record inventory differs'):
        campaign['run']()
    assert not campaign['output'].exists()


@pytest.mark.parametrize('field', ['model', 'environment', 'source_identity', 'source_manifest', 'path_budgets', 'max_steps_s'])
def test_mixed_equations_or_environment_fail_closed(campaign, field):
    def change(plan):
        if field == 'max_steps_s':
            plan[field][0] *= 2
        else:
            plan[field] = 'foreign'
    reseal(campaign['child']/'plan.json', change)
    refresh_request(campaign)
    with pytest.raises((ValueError, TypeError)):
        campaign['run']()
    assert not campaign['output'].exists()


@pytest.mark.parametrize('field', ['campaign_sha256', 'source_bundle_sha256', 'model_digest', 'current_path_source', 'spectrum'])
def test_resealed_unbound_grid_or_fabricated_evidence_is_rejected(campaign, field):
    def change(report):
        if field == 'model_digest':
            report['row'][field] = '0'*64
        elif field == 'current_path_source':
            report['row']['path_ledgers'][0][field] = 'wrong-grid-source'
        elif field == 'spectrum':
            # Recompute inner digests too: only reread evidence can reject this.
            row = c.r.decode(report['row'])
            row['candidate']['spectra']['retarded_response'] *= 1.01
            row['candidate']['candidate_digest'] = c.r.stream._spectrum_digest(row['candidate'])
            report['row'] = c.r.encode(row)
        else:
            report[field] = '0'*64
        report['row']['content_digest'] = c.r.row_digest(c.r.decode(report['row']))
    reseal(campaign['child']/'grid.json', change)
    refresh_request(campaign)
    with pytest.raises(ValueError, match='reference differs|identity|freshly validated'):
        campaign['run']()
    assert not campaign['output'].exists()


@pytest.mark.parametrize('side', [0, 1])
@pytest.mark.parametrize('damage', ['marker', 'source', 'controller'])
def test_capsule_marker_raw_sources_and_controller_bytes_are_bound(campaign, side, damage):
    root = campaign['roots'][side]
    file = {'marker': root/'.campaign-source-identity', 'source': root/'gabes/fixture.py',
            'controller': root/'tools/thermal_campaign_compare.py'}[damage]
    file.write_bytes(file.read_bytes()+b'\n')
    with pytest.raises(ValueError, match='fresh captured-source|capsule'):
        campaign['run']()
    assert not campaign['output'].exists()


def test_parent_capsule_audit_cannot_be_relabelled_as_child(campaign):
    campaign['audit']()
    reseal(campaign['proof'], lambda proof: proof.update(campaign_sha256=campaign['plan']['record_sha256']))
    with pytest.raises(ValueError, match='separate campaign-bound'):
        execute_after_proof(campaign)


@pytest.mark.parametrize('what', ['packet', 'plan', 'report', 'zip', 'lineage', 'controller'])
def test_final_byte_checks_reject_changes_after_all_numerical_reads(campaign, monkeypatch, what):
    original = t.nested_sum_audit
    def changed(*args):
        result = original(*args)
        child = campaign['child']
        filename = {'packet': next((child/'cache').glob('*.json')), 'plan': child/'plan.json',
            'report': child/'grid.json', 'zip': child/'sources.zip', 'lineage': child/'transfer.json',
            'controller': child/campaign['plan']['extension']['controller']['archive']}[what]
        filename.write_bytes(filename.read_bytes()+b'\n')
        return result
    monkeypatch.setattr(t, 'nested_sum_audit', changed)
    with pytest.raises(ValueError, match='bytes changed'):
        campaign['run']()
    assert not campaign['output'].exists()


def test_wrong_nested_mapping_is_rejected_across_stores(campaign, monkeypatch):
    monkeypatch.setattr(t, 'nested_indices', lambda *args: [(i, 2*((i+1) % 6)) for i in range(6)])
    with pytest.raises(ValueError, match='Nested map'):
        campaign['run']()


def pairing_inputs(campaign):
    coarse = c.r.decode(c.records.read_sealed(campaign['parent']/'grid.json'))['row']
    fine = c.r.decode(c.records.read_sealed(campaign['child']/'grid.json'))['row']
    model = c.r.default_model()
    stores = [x.ReadOnlyCache(campaign[name], campaign[plan], x.ByteWatch())
              for name, plan in (('parent', 'old'), ('child', 'plan'))]
    return (model, model.inflow(coarse['power'], coarse['seed']), model.inflow(fine['power'], fine['seed']),
            coarse, fine, c.pilot.step_plan(campaign['plan']['max_steps_s']), *stores)


@pytest.mark.parametrize('metric', c.r.METRICS)
def test_cross_pair_loop_independently_compares_all_eight_raw_metrics(campaign, monkeypatch, metric):
    args = pairing_inputs(campaign)
    child_cache, numerical = args[-1], args[5]
    original = child_cache.get
    def changed(model, path, spec, provider=None):
        packet, reference = original(model, path, spec, provider)
        if spec == numerical.reference[0]:
            packet = dict(packet)
            packet[metric] = packet[metric].copy()
            packet[metric].flat[0] += max(abs(packet[metric].flat[0]), 1e-30)*1e-6
        return packet, reference
    monkeypatch.setattr(child_cache, 'get', changed)
    with pytest.raises(ValueError, match='eight-metric packets'):
        x.nested_reuse_audit(*args)


def test_pair_loop_reads_both_stores_when_current_path_sources_are_identical(campaign, monkeypatch):
    args = pairing_inputs(campaign)
    assert args[1].path(0).source == args[2].path(0).source
    calls = [[], []]
    for side, cache in enumerate(args[-2:]):
        original = cache.get
        def tracked(model, path, spec, provider=None, original=original, side=side):
            calls[side].append((path.source, spec.method))
            return original(model, path, spec, provider)
        monkeypatch.setattr(cache, 'get', tracked)
    result = x.nested_reuse_audit(*args)
    assert result['passed'] and len(calls[0]) == len(calls[1]) == 30
    key = result['pairs'][0]['cache_records'][0]['key']
    (campaign['parent']/'cache'/(key+'.json')).unlink()
    assert (campaign['child']/'cache'/(key+'.json')).exists()
    with pytest.raises(ValueError, match='regular independent file'):
        x.nested_reuse_audit(*args)


def test_later_nested_map_preserves_within_face_offset_and_rejects_nonconsecutive():
    assert t.nested_indices(1, 2)[:4] == [(0, 0), (1, 1), (2, 4), (3, 5)]
    with pytest.raises(ValueError, match='consecutive'):
        t.nested_indices(0, 2)


def test_no_overwrite_and_no_campaign_archive_mutations(campaign):
    campaign['output'].write_bytes(b'owner data')
    before = [hashes(campaign[name]) for name in ('parent', 'child')]
    with pytest.raises(FileExistsError):
        x.launch(campaign['parent'], campaign['child'], campaign['parent']/'grid.json',
                 campaign['child']/'grid.json', campaign['output'])
    assert campaign['output'].read_bytes() == b'owner data'
    assert [hashes(campaign[name]) for name in ('parent', 'child')] == before


def test_real_fresh_subprocess_reaches_native_cache_reread_before_rejecting_fixture(base, tmp_path, monkeypatch, capfd):
    # Full captured source bundle, deliberately absent native cache. No monkey-
    # patched numeric implementation reaches the fresh child interpreter.
    parent, child = tmp_path/'parent', tmp_path/'child'
    c.create_campaign(parent, root=REAL_ROOT)
    shutil.copytree(parent, child)
    old = c.records.read_sealed(parent/'plan.json')
    coarse = c.records.read_sealed(base[0]/'parent/grid.json')
    del coarse['record_sha256']
    coarse.update(campaign_sha256=old['record_sha256'], source_identity=old['source_identity'],
        source_bundle_sha256=old['source_bundle_sha256'], source_manifest_after=old['source_manifest'])
    t.archive_controllers(parent, [(g.controller_identity(), Path(g.__file__).read_bytes())])
    c.write_record(parent/'grid.json', coarse)
    c.write_record(child/'grid.json', coarse)
    c.write_record(child/'transfer.json', {'synthetic_unreached': True})
    monkeypatch.setattr(x, 'ROOT', tmp_path)
    with pytest.raises(subprocess.CalledProcessError):
        x.launch(parent, child, parent/'grid.json', child/'grid.json', tmp_path/'output.json')
    stderr = capfd.readouterr().err
    assert 'Fresh grid requires complete passing path evidence' in stderr
    assert 'Audit input must be a regular independent file' in stderr
    assert not (tmp_path/'output.json').exists()
    assert not list((tmp_path/'.git/grand-challenge-runs').iterdir())
