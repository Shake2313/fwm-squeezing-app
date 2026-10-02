"""Read-only strata diagnostics: independent arithmetic and rejected evidence."""

import hashlib
import json
from pathlib import Path
import subprocess
from types import SimpleNamespace

import numpy as np
import pytest

from tools import thermal_campaign_diagnose as d
import test_thermal_campaign_compare as comparison_fixtures
from test_thermal_campaign_compare import rewrite_report
from test_thermal_campaign_grid import FixtureCache, synthetic_packet

c, g = d.c, d.g


class DiskFixtureCache(FixtureCache):
    """Honest sealed test files; the real native cache rejects their schema."""

    def __init__(self, plan, directory):
        super().__init__(plan)
        self.directory = Path(directory)
        self.directory.mkdir(exist_ok=True)
        for key in self.packets:
            self.persist(key)

    def persist(self, key):
        payload = {'synthetic_test_fixture': True, 'packet': c.r.encode(self.packets[key])}
        (self.directory/(key+'.json')).write_text(json.dumps(c.records.sealed(payload), allow_nan=False), encoding='utf-8')

    def get(self, model, path, spec, provider=None):
        assert provider is None, 'diagnostics must never solve'
        key = self.key(model, path, spec)
        self.calls.append((path.source, key))
        filename = self.directory/(key+'.json')
        record = c.records.read_sealed(filename)
        packet = c.r._validated_packet(model, path, c.r.decode(record['packet']))
        return packet, {'key': key, 'hit': True, 'record_sha256': record['record_sha256'],
                        'payload_digest': c.r.digest(packet), 'path': str(filename)}


def prepare_campaign(tmp_path, monkeypatch, *, full_sources=False):
    monkeypatch.setattr(comparison_fixtures, 'FixtureCache',
                        lambda plan: DiskFixtureCache(plan, tmp_path/'campaign'/'cache'))
    return comparison_fixtures.prepare_campaign(tmp_path, monkeypatch, full_sources=full_sources)


@pytest.fixture
def campaign(tmp_path, monkeypatch):
    return prepare_campaign(tmp_path, monkeypatch)


def _add_grid(directory, cache, power, seed, filename):
    inflow = cache.model.inflow(power, seed)
    for index in range(len(inflow.rate_s_inverse)):
        path = inflow.path(index)
        for spec in cache.numerical.primary+cache.numerical.reference:
            key = cache.key(cache.model, path, spec)
            cache.packets.setdefault(key, synthetic_packet(cache.model, path, spec))
            cache.persist(key)
    g.execute(directory, power, seed, filename)


def test_complete_diagnostic_reconstructs_native_sums_and_preserves_evidence(campaign):
    directory, plan, cache, reference, candidate, output = campaign
    originals = {path: path.read_bytes() for path in (reference, candidate, directory/'plan.json', directory/'sources.zip')}
    old_keys = set(cache.packets)
    report = c.r.decode(d.execute(directory, reference, candidate, output))
    assert report['audit_passed'] and report['new_solve_count'] == 0
    assert not report['thermal_ensemble_converged'] and not report['physical_optical_prediction']
    assert not report['historical_hit_flags_used_as_execution_evidence']
    diagnostic = report['directed_diagnostic']
    assert diagnostic['kind'] == 'ensemble_refinement' and diagnostic['diagnostic_only']
    assert diagnostic['candidate_grid'] == [1, 11] and diagnostic['reference_grid'] == [0, 11]
    left, right = [c.r.decode(c.records.read_sealed(path))['row'] for path in (reference, candidate)]
    expected = c.r.stream_errors(right['candidate']['spectra'], left['candidate']['spectra'], right['comparison_scales'])
    assert diagnostic['errors'] == expected
    assert diagnostic['passed'] == all(value <= .05 for value in expected.values())
    for partition in diagnostic['partitions'].values():
        groups = partition['groups']
        assert sorted(index for group in groups for index in group['reference_indices']) == list(range(6))
        assert sorted(index for group in groups for index in group['candidate_indices']) == list(range(12))
        assert sum(group['candidate_rate_s_inverse'] for group in groups) == pytest.approx(cache.model.inflow(1, 11).total_arrival_rate_s_inverse)
        for key in ('reference_reconstruction_errors', 'candidate_reconstruction_errors', 'difference_reconstruction_errors'):
            assert max(partition[key].values()) < 1e-12
        for metric in c.r.stream.STREAM_METRICS:
            expected_delta = right['candidate']['spectra'][metric]-left['candidate']['spectra'][metric]
            actual = sum(group['difference_spectra'][metric] for group in groups)
            np.testing.assert_allclose(actual, expected_delta, rtol=1e-12, atol=1e-28)
            information = partition['metrics'][metric]
            assert information['worst_relative_error'] == expected[metric]
            location = tuple(information['worst_location']['array_index'])
            projections = [np.asarray(group['metrics'][metric]['signed_projection']) for group in groups]
            error_array = np.asarray(information['relative_error'])
            np.testing.assert_allclose(sum(projections), error_array, rtol=1e-11, atol=1e-12)
            assert error_array[location] == expected[metric]
            assert len(information['groups_at_worst']) == len(groups)
            if metric.endswith('_by_source'):
                assert error_array.shape == (len(cache.model.convention().source_names), 3)
                assert information['worst_location']['source'] in cache.model.convention().source_names
    assert len(cache.calls) == 7*(6+12)  # 5 native gate + 1 direct + 1 strata; final byte hashes only
    assert set(cache.packets) == old_keys
    assert report['consumed_cache_byte_audit']['passed']
    assert len(report['consumed_cache_byte_audit']['records']) == 60
    for item in report['consumed_cache_byte_audit']['records']:
        raw = Path(item['path']).read_bytes()
        assert hashlib.sha256(raw).hexdigest() == item['file_sha256']
        assert len(raw) == item['raw_bytes']
    for path, raw in originals.items():
        assert path.read_bytes() == raw
    for identity in report['controllers']:
        raw = (directory/identity['archive']).read_bytes()
        assert hashlib.sha256(raw).hexdigest() == identity['raw_sha256']
    assert report['campaign_sha256'] == plan['record_sha256']


def test_seed_direction_uses_reference_norm_and_candidate_floor_without_averaging(campaign):
    directory, _, cache, reference, _, output = campaign
    candidate = output.with_name('independent.json')
    _add_grid(directory, cache, 0, 211, candidate)
    forward = c.r.decode(d.execute(directory, reference, candidate, output))['directed_diagnostic']
    reverse = c.r.decode(d.execute(directory, candidate, reference, output.with_name('reverse.json')))['directed_diagnostic']
    assert forward['kind'] == reverse['kind'] == 'independent_scramble'
    assert forward['candidate_grid'] == [0, 211] and reverse['candidate_grid'] == [0, 11]
    assert forward['errors']['greater'] != reverse['errors']['greater']
    for metric in c.r.stream.STREAM_METRICS:
        np.testing.assert_array_equal(forward['total_difference_spectra'][metric], -reverse['total_difference_spectra'][metric])
    # A tailored dark matrix independently pins which grid provides the floor.
    delta = np.array([np.eye(2)*3e-14], dtype=complex)
    dark = np.zeros_like(delta)
    attributed = d.metric_attribution([delta], delta, dark, 2.)
    assert attributed['denominator'][0] == 256*np.finfo(float).eps
    assert attributed['relative_error'][0] == pytest.approx(np.sqrt(2)*3e-14/(256*np.finfo(float).eps))


def test_signed_projection_tracks_cancellation_and_zero_difference():
    first = np.array([np.eye(2)*2], dtype=complex)
    second = np.array([-np.eye(2)], dtype=complex)
    result = d.metric_attribution([first, second], first+second, np.array([np.eye(2)]), 1.)
    assert result['relative_error'][0] == pytest.approx(1.)
    assert result['sum_relative_magnitudes'][0] == pytest.approx(3.)
    assert result['groups'][0]['signed_projection'][0] == pytest.approx(2.)
    assert result['groups'][1]['signed_projection'][0] == pytest.approx(-1.)
    zero = d.metric_attribution([first, -first], np.zeros_like(first), np.array([np.eye(2)]), 1.)
    assert zero['relative_error'][0] == 0 and zero['sum_relative_magnitudes'][0] == pytest.approx(4.)
    assert all(group['signed_projection'][0] == 0 for group in zero['groups'])


def test_complex_nonhermitian_attribution_keeps_real_inner_product():
    first = np.array([[[1j, 2+3j], [4-2j, -1j]]])
    second = np.array([[[2., -1j], [3j, .5]]])
    delta = first+second
    reference = np.array([np.eye(2)])
    result = d.metric_attribution([first, second], delta, reference, 1.)
    expected = np.vdot(delta[0], first[0]).real/np.linalg.norm(delta[0])/np.sqrt(2)
    assert result['groups'][0]['signed_projection'][0] == pytest.approx(expected)
    assert sum(item['signed_projection'][0] for item in result['groups']) == pytest.approx(np.linalg.norm(delta)/np.sqrt(2))


def test_fixed_partition_boundaries_keep_tails_and_put_equal_threshold_in_upper_bin():
    model = c.r.default_model()
    sigma = np.sqrt(c.r.c.KB*model.temperature_K/model.mass_kg)
    inflow = SimpleNamespace(rate_s_inverse=np.ones(6), face_index=np.arange(6),
        residence_time_s=np.array([.1, .25, .5, 1., 2., 40.])*1e-6,
        velocity_m_s=sigma*np.column_stack(([0., .5, 1., 2., 3., 100.], np.zeros(6), [-10., -2., -1., 0., 1., 2.])))
    result = d.partitions(model, inflow)
    assert result['inlet_face']['assignment'] == list(range(6))
    assert result['residence_time_us']['assignment'] == list(range(6))
    assert result['longitudinal_velocity_sigma']['assignment'] == list(range(6))
    assert result['transverse_speed_sigma']['assignment'] == [0, 1, 2, 3, 4, 4]
    assert result['transverse_speed_sigma']['groups'][-1]['indices'] == [4, 5]
    assert result['residence_time_us']['groups'][-1]['upper_exclusive'] is None


@pytest.mark.parametrize('metric', c.r.stream.STREAM_METRICS)
def test_changed_spectrum_resealed_to_match_its_digest_is_rejected(campaign, metric):
    directory, _, _, reference, candidate, output = campaign
    def change(report):
        report['row']['candidate']['spectra'][metric] *= 1.001
        report['row']['candidate']['candidate_digest'] = c.r.stream._spectrum_digest(report['row']['candidate'])
    rewrite_report(candidate, change)
    with pytest.raises(ValueError, match='freshly validated'):
        d.execute(directory, reference, candidate, output)
    assert not output.exists()


@pytest.mark.parametrize('damage', ['missing', 'changed', 'bad_source'])
def test_actual_packet_evidence_is_required(campaign, damage):
    directory, _, cache, reference, candidate, output = campaign
    path = cache.model.inflow(0, 11).path(2)
    if damage == 'missing':
        key = cache.key(cache.model, path, cache.numerical.reference[-1])
        del cache.packets[key]
        (cache.directory/(key+'.json')).unlink()
    elif damage == 'changed':
        for spec in cache.numerical.primary+cache.numerical.reference:
            key = cache.key(cache.model, path, spec)
            cache.packets[key]['greater'] *= 1.1
            cache.persist(key)
    else:
        rewrite_report(candidate, lambda report: report['row']['path_ledgers'][0].update(current_path_source='wrong'))
    with pytest.raises(ValueError, match='Fresh grid|freshly validated'):
        d.execute(directory, reference, candidate, output)
    assert not output.exists()


@pytest.mark.parametrize('field', ['source_identity', 'source_bundle_sha256', 'campaign_sha256'])
def test_mismatched_frozen_campaign_is_rejected_before_packets(campaign, field):
    directory, _, cache, reference, candidate, output = campaign
    rewrite_report(candidate, lambda report: report.update({field: '0'*64}))
    with pytest.raises(ValueError, match='reference differs'):
        d.execute(directory, reference, candidate, output)
    assert not cache.calls and not output.exists()


@pytest.mark.parametrize('before,after', [([0, 11], [0, 11]), ([1, 11], [0, 11]), ([0, 11], [2, 11]), ([0, 11], [1, 211])])
def test_ambiguous_or_reversed_refinement_comparisons_are_rejected(before, after):
    with pytest.raises(ValueError, match='distinct|adjacent'):
        d.comparison_kind({'selection': before}, {'selection': after})


def test_omitted_number_noise_cannot_reconstruct_groups(campaign, monkeypatch):
    directory, _, _, reference, candidate, output = campaign
    original = d.comparison._phase_terms
    def broken(packet):
        terms = original(packet)
        terms['poisson_number'] *= 0
        return terms
    monkeypatch.setattr(d.comparison, '_phase_terms', broken)
    with pytest.raises(ValueError, match='reconstruct the validated grid'):
        d.execute(directory, reference, candidate, output)
    assert not output.exists()


@pytest.mark.parametrize('damage', ['edit', 'remove', 'rename', 'replace_sealed'])
def test_cache_change_after_attribution_is_detected_by_final_hash(campaign, monkeypatch, damage):
    directory, _, cache, reference, candidate, output = campaign
    original = d.difference_diagnostic
    def mutate(*args, **kwargs):
        result = original(*args, **kwargs)
        path = cache.model.inflow(0, 11).path(0)
        key = cache.key(cache.model, path, cache.numerical.reference[0])
        filename = cache.directory/(key+'.json')
        if damage == 'edit':
            filename.write_bytes(filename.read_bytes()+b' ')
        elif damage == 'remove':
            filename.unlink()
        elif damage == 'rename':
            filename.rename(filename.with_suffix('.moved'))
        else:
            cache.packets[key]['retarded_response'] *= 1.001
            cache.persist(key)
        return result
    monkeypatch.setattr(d, 'difference_diagnostic', mutate)
    with pytest.raises((ValueError, FileNotFoundError)):
        d.execute(directory, reference, candidate, output)
    assert not output.exists()


@pytest.mark.parametrize('artifact', ['plan.json', 'sources.zip', 'candidate'])
def test_control_input_changes_during_diagnostic_are_rejected(campaign, monkeypatch, artifact):
    directory, _, _, reference, candidate, output = campaign
    original = d.difference_diagnostic
    def mutate(*args, **kwargs):
        result = original(*args, **kwargs)
        path = candidate if artifact == 'candidate' else directory/artifact
        path.write_bytes(path.read_bytes()+b' ')
        return result
    monkeypatch.setattr(d, 'difference_diagnostic', mutate)
    with pytest.raises(ValueError, match='changed during diagnostic'):
        d.execute(directory, reference, candidate, output)
    assert not output.exists()


def test_existing_output_and_corrupt_controller_archive_are_never_overwritten(campaign):
    directory, _, cache, reference, candidate, output = campaign
    output.write_bytes(b'existing evidence')
    for action in (d.execute, d.launch):
        with pytest.raises(FileExistsError):
            action(directory, reference, candidate, output)
    assert output.read_bytes() == b'existing evidence' and not cache.calls
    captured = d.captured_controllers()
    d.comparison.archive_controllers(directory, captured)
    archive = directory/captured[0][0]['archive']
    archive.write_bytes(b'changed')
    with pytest.raises(ValueError, match='archive changed'):
        d.execute(directory, reference, candidate, output.with_name('other.json'))
    assert archive.read_bytes() == b'changed'


@pytest.mark.parametrize('controller_index', [0, 1, 2])
def test_all_archived_controllers_are_rechecked_after_attribution(campaign, monkeypatch, controller_index):
    directory, _, _, reference, candidate, output = campaign
    original = d.difference_diagnostic
    archive = directory/d.captured_controllers()[controller_index][0]['archive']
    def mutate(*args, **kwargs):
        result = original(*args, **kwargs)
        archive.write_bytes(b'archive changed after attribution')
        return result
    monkeypatch.setattr(d, 'difference_diagnostic', mutate)
    with pytest.raises(ValueError, match='archive changed'):
        d.execute(directory, reference, candidate, output)
    assert not output.exists()
    assert archive.read_bytes() == b'archive changed after attribution'


def test_fresh_capsule_rejects_synthetic_reports_without_native_packets(tmp_path, monkeypatch, capfd):
    directory, _, _, reference, candidate, output = prepare_campaign(tmp_path, monkeypatch, full_sources=True)
    monkeypatch.setattr(d, 'ROOT', tmp_path)
    monkeypatch.setenv('PYTHONPATH', str(tmp_path/'wrong-source'))
    original = {path: path.read_bytes() for path in (directory/'cache').iterdir()}
    with pytest.raises(subprocess.CalledProcessError):
        d.launch(directory, reference, candidate, output)
    captured = capfd.readouterr()
    assert 'Fresh grid requires complete passing path evidence' in captured.err
    assert not output.exists()
    assert {path: path.read_bytes() for path in (directory/'cache').iterdir()} == original


def test_execute_cli_requires_all_controller_hashes(tmp_path):
    with pytest.raises(SystemExit) as error:
        d.main(['--execute', str(tmp_path), '--reference', 'a.json', '--candidate', 'b.json', '--output', 'out.json'])
    assert error.value.code == 2


@pytest.fixture
def byte_cache(tmp_path):
    key = 'a'*64
    filename = tmp_path/(key+'.json')
    filename.write_text(json.dumps(c.records.sealed({'value': 1})), encoding='utf-8')
    class Cache:
        directory = tmp_path
        calls = 0

        def key(self, *args):
            return key

        def get(self, *args, provider=None):
            assert provider is None
            self.calls += 1
            record = c.records.read_sealed(filename)
            return {}, {'key': key, 'hit': True, 'record_sha256': record['record_sha256']}

    cache = Cache()
    return cache, d.CacheByteGuard(cache), filename


def test_byte_guard_preserves_each_initial_get_and_final_check_is_only_bytes(byte_cache):
    cache, guard, filename = byte_cache
    guard.get(None, None, None)
    guard.get(None, None, None)
    assert cache.calls == 2
    (filename.parent/'unrelated-new-batch-record.json').write_bytes(b'new work outside consumed set')
    records = guard.verify()
    assert cache.calls == 2 and len(records) == 1
    assert records[0]['file_sha256'] == hashlib.sha256(filename.read_bytes()).hexdigest()


@pytest.mark.parametrize('damage', ['edit', 'remove', 'rename'])
def test_byte_guard_rejects_changed_bytes_before_repeated_native_get(byte_cache, damage):
    cache, guard, filename = byte_cache
    guard.get(None, None, None)
    if damage == 'edit':
        filename.write_bytes(filename.read_bytes()+b' ')
    elif damage == 'remove':
        filename.unlink()
    else:
        filename.rename(filename.with_suffix('.moved'))
    with pytest.raises((ValueError, FileNotFoundError)):
        guard.get(None, None, None)
    assert cache.calls == 1


def test_byte_guard_rejects_mutation_during_native_get(byte_cache, monkeypatch):
    cache, guard, filename = byte_cache
    original = cache.get
    def mutate(*args, **kwargs):
        result = original(*args, **kwargs)
        filename.write_bytes(filename.read_bytes()+b' ')
        return result
    monkeypatch.setattr(cache, 'get', mutate)
    with pytest.raises(ValueError, match='during native validation'):
        guard.get(None, None, None)
    assert not guard.snapshots


def test_byte_guard_binds_initial_snapshot_to_native_record_even_for_a_b_a_race(byte_cache, monkeypatch):
    cache, guard, filename = byte_cache
    original = cache.get
    before = filename.read_bytes()
    def replace_and_restore(*args, **kwargs):
        filename.write_text(json.dumps(c.records.sealed({'value': 2})), encoding='utf-8')
        result = original(*args, **kwargs)
        filename.write_bytes(before)
        return result
    monkeypatch.setattr(cache, 'get', replace_and_restore)
    with pytest.raises(ValueError, match='differs from exact cache byte snapshot'):
        guard.get(None, None, None)
    assert filename.read_bytes() == before and not guard.snapshots


@pytest.mark.parametrize('damage', ['bad_seal', 'duplicate_key', 'no_seal'])
def test_invalid_initial_sealed_snapshot_is_rejected_before_native_get(byte_cache, damage):
    cache, guard, filename = byte_cache
    if damage == 'bad_seal':
        record = json.loads(filename.read_text(encoding='utf-8'))
        record['value'] = 7
        filename.write_text(json.dumps(record), encoding='utf-8')
    elif damage == 'duplicate_key':
        filename.write_text('{"value": 1, "value": 2, "record_sha256": "bad"}', encoding='utf-8')
    else:
        filename.write_text('{"value": 1}', encoding='utf-8')
    with pytest.raises(ValueError, match='hash|duplicate|seal'):
        guard.get(None, None, None)
    assert cache.calls == 0


def test_byte_guard_never_accepts_a_provider(byte_cache):
    cache, guard, _ = byte_cache
    with pytest.raises(ValueError, match='forbids an atomic provider'):
        guard.get(None, None, None, provider=lambda *args: None)
    assert cache.calls == 0
