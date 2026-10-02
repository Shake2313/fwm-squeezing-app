"""Nested sidecar checks from temporary synthetic grids, never native solves."""

import hashlib
import json
from pathlib import Path
import subprocess

import numpy as np
import pytest

from tools import thermal_campaign_compare as t
from test_thermal_campaign_grid import FixtureCache, DiskFixtureCache, install_fixture_guard, synthetic_packet

c, g = t.c, t.g


def prepare_campaign(tmp_path, monkeypatch, *, full_sources=False):
    root = c.ROOT if full_sources else tmp_path/'source'
    if not full_sources:
        (root/'gabes').mkdir(parents=True)
        (root/'gabes/fixture.py').write_text('x = 1\n', encoding='utf-8')
    directory = tmp_path/'campaign'
    c.create_campaign(directory, root=root)
    plan = c.records.read_sealed(directory/'plan.json')
    cache = FixtureCache(plan)
    monkeypatch.setattr(c, 'load_plan', lambda *args, **kwargs: plan)
    monkeypatch.setattr(c, 'require_capsule', lambda *args: None)
    monkeypatch.setattr(c, 'CampaignCache', lambda *args, **kwargs: cache)
    install_fixture_guard(monkeypatch)
    monkeypatch.setattr(c, 'native_provider', lambda *args: pytest.fail('comparison must not solve'))
    coarse_file, fine_file = tmp_path/'coarse.json', tmp_path/'fine.json'
    g.execute(directory, 0, 11, coarse_file)
    g.execute(directory, 1, 11, fine_file)
    cache.calls.clear()
    return directory, plan, cache, coarse_file, fine_file, tmp_path/'comparison.json'


@pytest.fixture
def campaign(tmp_path, monkeypatch):
    return prepare_campaign(tmp_path, monkeypatch)


def rewrite_report(filename, change):
    report = c.r.decode(c.records.read_sealed(filename))
    del report['record_sha256']
    change(report)
    report['row']['content_digest'] = c.r.row_digest(report['row'])
    sealed = c.records.sealed(c.r.encode(report))
    filename.write_text(json.dumps(sealed, ensure_ascii=False, allow_nan=False), encoding='utf-8')


def decoded_rows(campaign):
    return [c.r.decode(c.records.read_sealed(path))['row'] for path in campaign[3:5]]


def test_complete_comparison_rereads_evidence_and_preserves_all_input_bytes(campaign):
    directory, plan, cache, coarse_file, fine_file, output = campaign
    originals = {p: p.read_bytes() for p in (coarse_file, fine_file, directory/'plan.json', directory/'sources.zip')}
    keys = set(cache.packets)
    report = c.r.decode(t.execute(directory, coarse_file, fine_file, output))
    assert report['audit_passed'] and report['nested_sum']['passed']
    assert report['new_solve_count'] == report['nested_reuse']['new_solve_count'] == 0
    assert not report['thermal_ensemble_converged'] and not report['physical_optical_prediction']
    assert report['convergence_gate']['passed'] is False and report['convergence_gate']['reasons']
    assert report['nested_reuse']['paired_path_count'] == 6
    assert report['nested_reuse']['shared_unique_jobs'] == 30
    assert report['nested_reuse']['new_fine_indices'] == [1, 3, 5, 7, 9, 11]
    assert len(cache.calls) == 6*(6+12)+10*6+6  # fresh grids, direct audits, pairs, new nodes
    assert set(cache.packets) == keys
    assert max(report['nested_sum']['errors'].values()) < 1e-13
    expected = c.r.stream_errors(report['rows'][1]['candidate']['spectra'],
        report['rows'][0]['candidate']['spectra'], report['rows'][0]['comparison_scales'])
    assert report['refinement_diagnostic']['errors'] == expected
    assert report['refinement_diagnostic']['budget'] == .05
    assert report['refinement_diagnostic']['passed'] == all(v <= .05 for v in expected.values())
    for identity in report['controllers']:
        raw = (directory/identity['archive']).read_bytes()
        assert hashlib.sha256(raw).hexdigest() == identity['raw_sha256']
        assert len(raw) == identity['raw_bytes']
    for path, raw in originals.items():
        assert path.read_bytes() == raw


@pytest.mark.parametrize('artifact', ['packet', 'plan', 'zip', 'controller', 'helper'])
def test_comparison_refuses_inputs_changed_after_last_cache_read(campaign, monkeypatch, artifact):
    directory, plan, original_cache, _, _, output = campaign
    cache = DiskFixtureCache(plan, directory/'cache', original_cache.packets)
    monkeypatch.setattr(c, 'CampaignCache', lambda *args, **kwargs: cache)
    coarse_file, fine_file = output.with_name('disk-coarse.json'), output.with_name('disk-fine.json')
    g.execute(directory, 0, 11, coarse_file)
    g.execute(directory, 1, 11, fine_file)
    original = t.nested_sum_audit

    def mutate(*args, **kwargs):
        result = original(*args, **kwargs)
        key = cache.key(cache.model, cache.model.inflow(0, 11).path(0), cache.numerical.reference[0])
        captured = t.captured_controllers()
        filename = {'packet': cache.directory/(key+'.json'), 'plan': directory/'plan.json',
                    'zip': directory/'sources.zip', 'controller': directory/captured[0][0]['archive'],
                    'helper': directory/captured[1][0]['archive']}[artifact]
        filename.write_bytes(filename.read_bytes()+b' ')
        return result

    monkeypatch.setattr(t, 'nested_sum_audit', mutate)
    with pytest.raises(ValueError, match='bytes changed|controller archive differs'):
        t.execute(directory, coarse_file, fine_file, output)
    assert not output.exists()


@pytest.mark.parametrize('controller_index', [0, 1])
def test_comparison_archive_capture_binds_original_controller_bytes(campaign, monkeypatch, controller_index):
    directory, _, cache, coarse_file, fine_file, output = campaign
    original = t.archive_controllers

    def mutate(folder, captured):
        original(folder, captured)
        identity = captured[controller_index][0]
        (Path(folder)/identity['archive']).write_bytes(b'changed before initial archive snapshot')

    monkeypatch.setattr(t, 'archive_controllers', mutate)
    with pytest.raises(ValueError, match='Archived controller bytes differ'):
        t.execute(directory, coarse_file, fine_file, output)
    assert not output.exists() and not cache.calls


def test_grid_report_seal_is_bound_to_initial_raw_snapshot_during_a_b_a_race(campaign, monkeypatch):
    directory, plan, _, coarse_file, _, _ = campaign
    raw = coarse_file.read_bytes()
    saved = c.records.read_sealed(coarse_file)
    substituted = c.records.sealed({name: value for name, value in saved.items() if name != 'record_sha256'}
                                  | {'scope': 'Different sealed record observed only by a second native read'})
    native_reads = []

    def raced_read(filename):
        assert Path(filename) == coarse_file
        native_reads.append(filename)
        # Native read observes B while both raw reads observe A again.
        return substituted

    monkeypatch.setattr(c.records, 'read_sealed', raced_read)
    report, reference, _, _ = t.read_grid_report(coarse_file, directory, plan)
    assert not native_reads
    assert report['record_sha256'] == reference['record_sha256'] == saved['record_sha256']
    assert reference['file_sha256'] == hashlib.sha256(raw).hexdigest()
    assert coarse_file.read_bytes() == raw


def test_later_refinement_keeps_within_face_offsets_instead_of_doubling_global_indices(campaign):
    directory, _, cache, _, coarse_file, output = campaign
    inflow = cache.model.inflow(2, 11)
    for index in range(len(inflow.rate_s_inverse)):
        path = inflow.path(index)
        for spec in cache.numerical.primary+cache.numerical.reference:
            cache.packets.setdefault(cache.key(cache.model, path, spec), synthetic_packet(cache.model, path, spec))
    fine_file = output.with_name('p2.json')
    g.execute(directory, 2, 11, fine_file)
    report = t.execute(directory, coarse_file, fine_file, output)
    pairs = report['nested_reuse']['pairs']
    assert [(row['coarse_index'], row['fine_index']) for row in pairs[:4]] == [(0, 0), (1, 1), (2, 4), (3, 5)]
    assert report['nested_sum']['passed'] and not report['thermal_ensemble_converged']


@pytest.mark.parametrize('damage', ['missing', 'changed_packet'])
def test_saved_success_does_not_replace_actual_cache_evidence(campaign, damage):
    directory, _, cache, coarse_file, fine_file, output = campaign
    path = cache.model.inflow(0, 11).path(2)
    if damage == 'missing':
        del cache.packets[cache.key(cache.model, path, cache.numerical.reference[-1])]
    else:
        # Change all five equally: path convergence still passes. The saved row
        # must nevertheless be rejected because the actual physical data changed.
        for spec in cache.numerical.primary+cache.numerical.reference:
            cache.packets[cache.key(cache.model, path, spec)]['retarded_response'] *= 1.01
    with pytest.raises(ValueError, match='Fresh grid|freshly validated'):
        t.execute(directory, coarse_file, fine_file, output)
    assert not output.exists()


@pytest.mark.parametrize('location', ['ledger', 'aggregation'])
@pytest.mark.parametrize('field', t.REFERENCE_KEYS)
def test_resealed_fabricated_cache_references_cannot_pass(campaign, location, field):
    directory, _, _, coarse_file, fine_file, output = campaign
    def change(report):
        records = (report['row']['path_ledgers'][0]['cache_records'] if location == 'ledger'
                   else report['aggregation_audit']['cache_records'])
        records[0][field] = 'a'*64
    rewrite_report(fine_file, change)
    with pytest.raises(ValueError, match='freshly validated|references differ'):
        t.execute(directory, coarse_file, fine_file, output)
    assert not output.exists()


@pytest.mark.parametrize('metric', c.r.stream.STREAM_METRICS)
def test_resealed_rehashed_spectrum_is_checked_against_actual_packets(campaign, metric):
    directory, _, _, coarse_file, fine_file, output = campaign
    def change(report):
        candidate = report['row']['candidate']
        candidate['spectra'][metric] *= 1.001
        candidate['candidate_digest'] = c.r.stream._spectrum_digest(candidate)
    rewrite_report(fine_file, change)
    with pytest.raises(ValueError, match='freshly validated'):
        t.execute(directory, coarse_file, fine_file, output)
    assert not output.exists()


@pytest.mark.parametrize('field', ['current_path_source', 'target_digest', 'physical_path_key', 'independent_reference'])
def test_resealed_wrong_source_digest_or_gate_cannot_be_reused(campaign, field):
    directory, _, _, coarse_file, fine_file, output = campaign
    def change(report):
        ledger = report['row']['path_ledgers'][2]
        if field == 'independent_reference':
            ledger[field]['mean_outer'] = 1e-12
        else:
            ledger[field] = 'fabricated-current-path'
    rewrite_report(fine_file, change)
    with pytest.raises(ValueError, match='freshly validated'):
        t.execute(directory, coarse_file, fine_file, output)
    assert not output.exists()


def test_historical_hit_flags_and_old_storage_paths_do_not_count_as_new_solves(campaign):
    directory, _, _, coarse_file, fine_file, output = campaign
    def change(report):
        for ledger in report['row']['path_ledgers']:
            for record in ledger['cache_records']:
                record.update(path='old-location/'+record['key']+'.json', hit=False)
        for record in report['aggregation_audit']['cache_records']:
            record.update(path='old-location/'+record['key']+'.json', hit=False)
        report['row']['cache_misses'] = report['row']['cache_hits']
        report['row']['cache_hits'] = 0
    rewrite_report(coarse_file, change)
    rewrite_report(fine_file, change)
    report = t.execute(directory, coarse_file, fine_file, output)
    assert report['audit_passed'] and report['new_solve_count'] == 0
    assert [row['cache_hits'] for row in report['rows']] == [30, 60]
    assert not report['nested_reuse']['historical_hit_flags_used_as_execution_evidence']


@pytest.mark.parametrize('field', ['campaign_sha256', 'source_identity', 'source_bundle_sha256', 'model_digest', 'plan_digest'])
def test_wrong_campaign_source_or_model_is_rejected_before_cache_reads(campaign, field):
    directory, _, cache, coarse_file, fine_file, output = campaign
    def change(report):
        target = report['row'] if field.endswith('_digest') else report
        target[field] = '0'*64
    rewrite_report(fine_file, change)
    with pytest.raises(ValueError, match='reference differs|identity'):
        t.execute(directory, coarse_file, fine_file, output)
    assert not output.exists() and not cache.calls


def test_wrong_nested_map_is_detected_from_physical_paths(campaign, monkeypatch):
    _, _, cache, _, _, _ = campaign
    coarse, fine = decoded_rows(campaign)
    monkeypatch.setattr(t, 'nested_indices', lambda *args: [(i, 2*((i+1) % 6)) for i in range(6)])
    with pytest.raises(ValueError, match='Nested map'):
        t.nested_reuse_audit(cache.model, cache.model.inflow(0, 11), cache.model.inflow(1, 11),
                            coarse, fine, cache.numerical, cache)


@pytest.mark.parametrize('metric', c.r.METRICS)
def test_all_eight_raw_metrics_are_compared_bitwise(campaign, monkeypatch, metric):
    _, _, cache, _, _, _ = campaign
    coarse, fine = decoded_rows(campaign)
    old, new = cache.model.inflow(0, 11), cache.model.inflow(1, 11)
    original = cache.get
    def changed(model, path, spec, provider=None):
        packet, record = original(model, path, spec, provider)
        if path.source == new.path(2).source and spec == cache.numerical.reference[0]:
            packet = dict(packet)
            packet[metric] = packet[metric].copy()
            packet[metric].flat[0] += max(abs(packet[metric].flat[0]), 1e-30)*1e-6
        return packet, record
    monkeypatch.setattr(cache, 'get', changed)
    with pytest.raises(ValueError, match='eight-metric'):
        t.nested_reuse_audit(cache.model, old, new, coarse, fine, cache.numerical, cache)


def test_signed_zero_and_dtype_changes_are_not_bitwise_identity():
    assert not t._bitwise_equal(np.array([0.]), np.array([-0.]))
    assert not t._bitwise_equal(np.array([1.], dtype='float32'), np.array([1.], dtype='float64'))


def test_split_identity_detects_missing_poisson_number_contribution(campaign):
    _, _, cache, _, _, _ = campaign
    coarse, fine = decoded_rows(campaign)
    inflow = cache.model.inflow(1, 11)
    reuse = {'new_fine_indices': [1, 3, 5, 7, 9, 11]}
    fine['candidate']['spectra']['poisson_number'] *= 1.01
    result = t.nested_sum_audit(cache.model, inflow, coarse, fine, cache.numerical, cache, reuse)
    assert not result['passed'] and result['errors']['poisson_number'] > .001


def test_output_and_controller_archives_are_exclusive(campaign):
    directory, _, cache, coarse_file, fine_file, output = campaign
    output.write_bytes(b'previous comparison')
    for action in (t.execute, t.launch):
        with pytest.raises(FileExistsError):
            action(directory, coarse_file, fine_file, output)
    assert output.read_bytes() == b'previous comparison' and not cache.calls
    captured = t.captured_controllers()
    t.archive_controllers(directory, captured)
    before = [(directory/item[0]['archive']).stat().st_mtime_ns for item in captured]
    t.archive_controllers(directory, captured)
    assert before == [(directory/item[0]['archive']).stat().st_mtime_ns for item in captured]
    archive = directory/captured[0][0]['archive']
    archive.write_bytes(b'corrupt captured controller')
    with pytest.raises(ValueError, match='archive changed'):
        t.archive_controllers(directory, captured)
    assert archive.read_bytes() == b'corrupt captured controller'


def test_fresh_subprocess_cannot_accept_synthetic_reports_without_actual_native_cache(tmp_path, monkeypatch):
    directory, plan, _, coarse_file, fine_file, output = prepare_campaign(tmp_path, monkeypatch, full_sources=True)
    source_raw = (directory/'sources.zip').read_bytes()
    monkeypatch.setattr(t, 'ROOT', tmp_path)
    monkeypatch.setenv('PYTHONPATH', str(tmp_path/'untrusted-modules'))
    with pytest.raises(subprocess.CalledProcessError):
        t.launch(directory, coarse_file, fine_file, output)
    assert not output.exists()
    assert list((directory/'cache').iterdir()) == []
    assert (directory/'sources.zip').read_bytes() == source_raw
    assert c.records.read_sealed(directory/'plan.json') == plan
    assert list((tmp_path/'.git/grand-challenge-runs').iterdir()) == []
    archived = list((directory/'controllers').glob('*.py'))
    assert {path.stem.split('-')[0] for path in archived} == {'thermal_campaign_compare', 'thermal_campaign_grid'}
