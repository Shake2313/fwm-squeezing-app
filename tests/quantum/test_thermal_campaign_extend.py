"""Transfer/publication contracts; temporary synthetic packets, no native solve."""

from copy import deepcopy
import json
import os
from pathlib import Path
import shutil
import subprocess

import pytest

from tools import thermal_campaign_extend as e
from test_thermal_campaign_grid import synthetic_packet
from test_thermal_campaign import native_numbers

c = e.c


def fixture_record(parent, plan, path_index=1, level=0):
    model = c.r.default_model()
    path = model.inflow(0, 11).path(path_index)
    numerical = c.pilot.step_plan(plan['max_steps_s'])
    spec = (numerical.primary+numerical.reference)[level]
    cache = c.CampaignCache(parent/'cache', plan)
    packet = synthetic_packet(model, path, spec)
    packet['numerics'] = dict(native_numbers(path, spec), synthetic_test_fixture=True)
    packet = c.r._validated_packet(model, path, packet)
    key = cache.key(model, path, spec)
    filename = parent/'cache'/(key+'.json')
    c.write_record(filename, {'schema': 'gabes-native-thermal-packet-v1',
        'identity': c.r.encode(cache.identity(model, path, spec)), 'key': key,
        'original_path_source': path.source, 'packet': packet,
        'source_before': plan['source_manifest'], 'source_after': plan['source_manifest'],
        'source_stable_during_run': True, 'synthetic_contract_fixture': True})
    return filename


def hashes(root):
    return {p.relative_to(root).as_posix(): e.sha(p.read_bytes()) for p in root.rglob('*') if p.is_file()}


def reseal(path, change):
    record = c.records.read_sealed(path)
    del record['record_sha256']
    change(record)
    path.write_text(json.dumps(c.records.sealed(record), ensure_ascii=False, allow_nan=False), encoding='utf-8')


@pytest.fixture(scope='module')
def base_parent(tmp_path_factory):
    parent = tmp_path_factory.mktemp('extension-base')/'parent'
    c.create_campaign(parent, root=c.ROOT)
    plan = c.records.read_sealed(parent/'plan.json')
    for level in range(5):
        filename = fixture_record(parent, plan, level=level)
        filename.write_bytes(filename.read_bytes().replace(b'\n', b'\r\n'))
    (parent/'path-old.json').write_bytes(b'historical report must not be imported')
    (parent/'controllers').mkdir()
    (parent/'controllers/original.py').write_bytes(b'# parent archive stays untouched\r\n')
    return parent


@pytest.fixture
def campaign(base_parent, tmp_path, monkeypatch):
    parent, target = tmp_path/'parent', tmp_path/'extension'
    shutil.copytree(base_parent, parent)
    plan = c.records.read_sealed(parent/'plan.json')
    required = []
    monkeypatch.setattr(c, 'require_capsule', lambda record: required.append(record['record_sha256']))
    native_validation = c.validate_native
    def test_validation(packet, path, spec):
        # Production rejects these explicitly synthetic records. Only this test
        # process removes the fixture marker, then exercises every native-number
        # contract against the real validator; subprocess tests keep it intact.
        numbers = dict(packet['numerics'])
        assert numbers.pop('synthetic_test_fixture') is True
        native_validation({'numerics': numbers}, path, spec)
    monkeypatch.setattr(c, 'validate_native', test_validation)
    monkeypatch.setattr(c, 'native_provider', lambda *args: pytest.fail('transfer must not solve'))
    def local_verifier(stage, control, raw):
        result = e.verify_stage(stage, expected_controller_sha256=control['raw_sha256'],
                                expected_grid_sha256=control['dependencies'][0]['raw_sha256'])
        return c.records.sealed(c.r.encode(result))
    monkeypatch.setattr(e, 'verify_in_fresh_capsule', local_verifier)
    return parent, target, plan, required


def patch_link_count(monkeypatch, root, count):
    original = Path.stat
    def changed(path, *args, **kwargs):
        value = original(path, *args, **kwargs)
        if path.parent.name == 'cache' and path.is_relative_to(root):
            fields = list(value)
            fields[3] = count
            return os.stat_result(fields)
        return value
    monkeypatch.setattr(Path, 'stat', changed)
    return original


@pytest.mark.parametrize('count', [0, 1, 2])
def test_link_count_policy_preserves_files_and_exclusive_destination(tmp_path, monkeypatch, count):
    parent, target = tmp_path/'parent', tmp_path/'target'
    cache = parent/'cache'
    cache.mkdir(parents=True)
    key = 'a'*64
    file = cache/(key+'.json')
    file.write_bytes(b'original packet bytes')
    target.mkdir()
    (target/'owner.txt').write_bytes(b'existing output')
    before = hashes(tmp_path)
    patch_link_count(monkeypatch, tmp_path, count)
    assert e.regular_cache_file(file) is (count in (0, 1))
    if count > 1:
        with pytest.raises(ValueError, match='nonregular file'):
            e.cache_files(cache)
    else:
        assert e.cache_files(cache) == {key: file}
    # Destination protection precedes reading even a malformed parent plan.
    with pytest.raises(FileExistsError):
        e.launch(parent, target, (2, 3, 4))
    assert hashes(tmp_path) == before


@pytest.mark.parametrize('kind', ['symlink', 'directory'])
def test_link_count_zero_does_not_allow_symlinks_or_nonregular_files(tmp_path, monkeypatch, kind):
    cache = tmp_path/'cache'
    cache.mkdir()
    file = cache/('a'*64+'.json')
    if kind == 'directory':
        file.mkdir()
    else:
        file.write_bytes(b'unchanged')
        original = Path.is_symlink
        monkeypatch.setattr(Path, 'is_symlink', lambda path: path == file or original(path))
    patch_link_count(monkeypatch, tmp_path, 0)
    assert not e.regular_cache_file(file)
    with pytest.raises(ValueError, match='nonregular file'):
        e.cache_files(cache)
    assert file.exists()


@pytest.mark.parametrize('virtual_count', [None, 0], ids=['native', 'link_count_zero'])
def test_exact_transfer_preserves_parent_bytes_and_rebinds_only_new_declaration(campaign, monkeypatch, virtual_count):
    parent, target, old, required = campaign
    actual_stat = Path.stat
    if virtual_count is not None:
        patch_link_count(monkeypatch, target.parent, virtual_count)
    before = hashes(parent)
    report = e.execute(parent, target, (2, 3, 4))
    new = c.load_plan(target)
    assert hashes(parent) == before
    assert new['record_sha256'] != old['record_sha256']
    assert required == [old['record_sha256'], new['record_sha256']]
    assert new['powers'] == [2, 3, 4] and new['seeds'] == old['seeds']
    assert new['historical_numerical_results_reused'] is True
    assert new['extension']['imported_raw_record_count'] == 5
    assert (parent/'sources.zip').read_bytes() == (target/'sources.zip').read_bytes()
    assert report == c.records.read_sealed(target/'transfer.json')
    assert report['transfer_passed'] and report['imported_raw_record_count'] == 5
    assert len(report['expected_native_keys']) == 1440 and len(report['missing_native_keys']) == 1435
    assert report['new_solve_count'] == 0 and not report['path_and_grid_evidence_imported']
    assert not report['thermal_ensemble_converged'] and not report['physical_optical_prediction']
    assert set(p.name for p in target.iterdir()) == {'plan.json', 'sources.zip', 'cache', 'controllers', 'transfer.json'}
    for row in report['imports']:
        one, two = parent/row['cache_file'], target/row['cache_file']
        assert one.read_bytes() == two.read_bytes()
        assert not os.path.samefile(one, two) and actual_stat(two).st_nlink == 1
        assert row['parent_selection'][:3] == [0, 11, 1]
        assert row['target_selection'][:3] == [2, 11, 4]
        assert c.records.read_sealed(two)['record_sha256'] == row['record_sha256']
    assert not list(target.parent.glob('.extension-staging-*'))


def test_append_only_parent_arrivals_after_snapshot_are_not_mutation_or_extra_imports(campaign):
    parent, target, plan, _ = campaign
    captured = e.present_snapshot(parent, e.snapshot(parent, plan), e.controller()[0])
    extra = fixture_record(parent, plan, path_index=2, level=0)
    before = hashes(parent)
    report = e.execute(parent, target, (2, 3, 4), cache_snapshot=captured)
    assert report['imported_raw_record_count'] == 5
    assert not (target/'cache'/extra.name).exists()
    assert hashes(parent) == before


@pytest.mark.parametrize('field', ['requested_powers', 'required_seeds', 'missing_native_keys'])
def test_resealed_request_coverage_must_match_staged_declaration(campaign, monkeypatch, field):
    parent, target, _, _ = campaign
    original = e.verify_in_fresh_capsule
    def changed(stage, control, raw):
        reseal(stage/'transfer-request.json', lambda request: request[field].pop())
        return original(stage, control, raw)
    monkeypatch.setattr(e, 'verify_in_fresh_capsule', changed)
    with pytest.raises(ValueError, match='Staged cache inventory differs'):
        e.execute(parent, target, (2, 3, 4))
    assert not target.exists()


def test_resealed_exclusion_cannot_hide_a_compatible_captured_record(campaign, monkeypatch):
    parent, target, _, _ = campaign
    original = e.verify_in_fresh_capsule
    def changed(stage, control, raw):
        def exclude(request):
            row = request['imports'].pop()
            request['excluded_parent_records'].append(row)
            (stage/row['cache_file']).unlink()
            request['missing_native_keys'] = sorted(request['missing_native_keys']+[row['key']])
        reseal(stage/'transfer-request.json', exclude)
        return original(stage, control, raw)
    monkeypatch.setattr(e, 'verify_in_fresh_capsule', changed)
    with pytest.raises(ValueError, match='Excluded record belongs to the target'):
        e.execute(parent, target, (2, 3, 4))
    assert not target.exists()


@pytest.mark.parametrize('kind', ['corrupt', 'missing_after_snapshot', 'changed_after_snapshot'])
def test_captured_present_records_cannot_be_repaired_or_silently_skipped(campaign, kind):
    parent, target, plan, _ = campaign
    file = next((parent/'cache').glob('*.json'))
    if kind == 'corrupt':
        file.write_bytes(b'{invalid')
    captured = e.present_snapshot(parent, e.snapshot(parent, plan), e.controller()[0])
    if kind == 'missing_after_snapshot':
        file.unlink()
    elif kind == 'changed_after_snapshot':
        file.write_bytes(file.read_bytes()+b'\n')
    before = hashes(parent)
    with pytest.raises((ValueError, FileNotFoundError)):
        e.execute(parent, target, (2, 3, 4), cache_snapshot=captured)
    assert not target.exists() and hashes(parent) == before


@pytest.mark.parametrize('kind', ['environment', 'model', 'source_identity', 'grid_count', 'rate', 'coordinates', 'occupancy'])
def test_resealed_parent_inconsistency_is_rejected_including_uncached_grids(campaign, kind):
    parent, target, _, _ = campaign
    def change(plan):
        grid = plan['grids'][-1]  # No cached records on this seed; still validate.
        if kind == 'environment':
            plan['environment']['python'] = 'incompatible'
        elif kind in ('model', 'source_identity'):
            plan[kind] = 'forged'
        elif kind == 'grid_count':
            grid['paths'].pop()
        elif kind == 'rate':
            grid['paths'][-1]['rate_s_inverse'] *= 2
        elif kind == 'coordinates':
            grid['paths'][-1]['physical_path']['velocity_m_s'][0] += 1
        else:
            grid['mean_occupancy'] *= 2
    reseal(parent/'plan.json', change)
    with pytest.raises((ValueError, TypeError)):
        e.execute(parent, target, (2, 3, 4))
    assert not target.exists()


def test_resealed_target_grid_damage_is_detected_before_publication(campaign, monkeypatch):
    parent, target, _, _ = campaign
    original = c.make_plan
    def damaged(*args, **kwargs):
        plan = original(*args, **kwargs)
        plan['grids'][-1]['paths'][-1]['rate_s_inverse'] *= 2
        return plan
    monkeypatch.setattr(c, 'make_plan', damaged)
    with pytest.raises(ValueError, match='arrival rate'):
        e.execute(parent, target, (2, 3, 4))
    assert not target.exists()


@pytest.mark.parametrize('kind', ['native_parameters', 'model_identity', 'source_manifest'])
def test_resealed_cache_inconsistency_cannot_be_imported(campaign, kind):
    parent, target, _, _ = campaign
    file = next((parent/'cache').glob('*.json'))
    def change(record):
        if kind == 'native_parameters':
            record['packet']['numerics']['requested_parameters'] = {}
        elif kind == 'model_identity':
            record['identity']['model'] = 'different-model'
        else:
            record['source_after']['source_sha256'] = '0'*64
    reseal(file, change)
    with pytest.raises(ValueError):
        e.execute(parent, target, (2, 3, 4))
    assert not target.exists()


def test_source_zip_and_launch_controller_hashes_are_bound(campaign):
    parent, target, _, _ = campaign
    with pytest.raises(ValueError, match='controller'):
        e.execute(parent, target, (2, 3, 4), expected_controller_sha256='0'*64)
    with pytest.raises(ValueError, match='Grid helper'):
        e.execute(parent, target, (2, 3, 4), expected_grid_sha256='0'*64)
    with (parent/'sources.zip').open('ab') as handle:
        handle.write(b'changed')
    with pytest.raises(ValueError, match='source archive'):
        e.execute(parent, target, (2, 3, 4))
    assert not target.exists()


@pytest.mark.parametrize('powers', [(0, 1, 2), (2, 3), (2, 3, 3), (True, 3, 4)])
def test_declaration_must_be_valid_and_finer(campaign, powers):
    parent, target, _, _ = campaign
    with pytest.raises(ValueError):
        e.execute(parent, target, powers)
    assert not target.exists()


def test_existing_and_parent_nested_targets_are_rejected_without_changes(campaign):
    parent, target, _, _ = campaign
    with pytest.raises(ValueError, match='separate'):
        e.execute(parent, parent/'nested', (2, 3, 4))
    target.mkdir()
    (target/'owner.txt').write_bytes(b'user data')
    before = hashes(parent)
    for action in (e.execute, e.launch):
        with pytest.raises(FileExistsError):
            action(parent, target, (2, 3, 4))
    assert hashes(parent) == before and (target/'owner.txt').read_bytes() == b'user data'


def test_failure_during_second_packet_copy_never_publishes_partial_campaign(campaign, monkeypatch):
    parent, target, _, _ = campaign
    before = hashes(parent)
    original, copies = Path.open, []
    def failing(path, mode='r', *args, **kwargs):
        if mode == 'xb' and path.parent.name == 'cache' and '.extension-staging-' in str(path):
            copies.append(path)
            if len(copies) == 2:
                raise OSError('injected second-copy failure')
        return original(path, mode, *args, **kwargs)
    monkeypatch.setattr(Path, 'open', failing)
    with pytest.raises(OSError, match='second-copy'):
        e.execute(parent, target, (2, 3, 4))
    assert len(copies) == 2 and not target.exists() and hashes(parent) == before
    assert not list(target.parent.glob('.extension-staging-*'))


def test_destination_race_is_exclusive_and_preserves_external_target(campaign, monkeypatch):
    parent, target, _, _ = campaign
    original = e.publish
    def race(stage, destination):
        Path(destination).mkdir()
        (Path(destination)/'owner.txt').write_bytes(b'concurrent owner')
        return original(stage, destination)
    monkeypatch.setattr(e, 'publish', race)
    with pytest.raises(FileExistsError):
        e.execute(parent, target, (2, 3, 4))
    assert {p.name for p in target.iterdir()} == {'owner.txt'}
    assert (target/'owner.txt').read_bytes() == b'concurrent owner'
    assert not list(target.parent.glob('.extension-staging-*'))


def test_real_fresh_parent_subprocess_rejects_synthetic_native_evidence(campaign, monkeypatch, capfd):
    parent, target, _, _ = campaign
    monkeypatch.setattr(e, 'ROOT', target.parent)
    with pytest.raises(subprocess.CalledProcessError):
        e.launch(parent, target, (2, 3, 4))
    assert 'Native numerical parameters differ or synthetic record supplied' in capfd.readouterr().err
    assert not target.exists()


def test_empty_campaign_extension_passes_both_real_capsules_without_any_native_solve(campaign, monkeypatch):
    parent, target, old, _ = campaign
    for path in (parent/'cache').glob('*.json'):
        path.unlink()
    before = hashes(parent)
    monkeypatch.setattr(e, 'ROOT', target.parent)
    e.launch(parent, target, (2, 3, 4))
    new = c.load_plan(target)
    report = c.records.read_sealed(target/'transfer.json')
    assert report['transfer_passed'] and report['imported_raw_record_count'] == 0
    assert report['verification']['target_campaign_sha256'] == new['record_sha256'] != old['record_sha256']
    assert new['historical_numerical_results_reused'] is False
    assert report['new_solve_count'] == 0 and not report['thermal_ensemble_converged']
    assert hashes(parent) == before
    assert not list((target.parent/'.git/grand-challenge-runs').iterdir())


def test_verification_must_use_target_marker_not_parent_marker(campaign, monkeypatch, capfd):
    parent, target, old, _ = campaign
    for path in (parent/'cache').glob('*.json'):
        path.unlink()
    original_extract = c.extract_bundle
    def wrong_marker(directory, folder):
        original_extract(directory, folder)
        (Path(folder)/'.campaign-source-identity').write_text(old['record_sha256'], encoding='utf-8')
    monkeypatch.setattr(c, 'extract_bundle', wrong_marker)
    # Restore the real second-process verifier; the first-process capsule guard
    # is explicitly stubbed by the fixture, but the fresh verifier is untouched.
    verifier = e.verify_stage.__globals__['subprocess']
    assert verifier is subprocess
    def fresh(stage, control, raw):
        with pytest.MonkeyPatch.context() as patch:
            patch.setattr(e, 'ROOT', target.parent)
            # Reuse the original function captured before fixture patching.
            return REAL_VERIFY(stage, control, raw)
    monkeypatch.setattr(e, 'verify_in_fresh_capsule', fresh)
    with pytest.raises(subprocess.CalledProcessError):
        e.execute(parent, target, (2, 3, 4))
    assert 'Execute through --run in a fresh captured-source interpreter' in capfd.readouterr().err
    assert not target.exists()


REAL_VERIFY = e.verify_in_fresh_capsule
