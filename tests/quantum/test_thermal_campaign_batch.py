"""Controller tests use synthetic packets; they are never atomic evidence."""

from concurrent.futures import ThreadPoolExecutor
import hashlib
from pathlib import Path

import pytest

from analysis.grand_challenge import thermal_campaign as c
from tools import thermal_campaign_batch as batch
from tests.quantum.test_rb_thermal_ensemble import synthetic_packet
from tests.quantum.test_thermal_campaign import native_numbers


@pytest.fixture
def campaign(tmp_path, monkeypatch):
    root = tmp_path/'source'
    (root/'gabes').mkdir(parents=True)
    (root/'gabes/test.py').write_text('# synthetic test dependency\n')
    directory = tmp_path/'campaign'
    directory.mkdir()
    plan = c.write_record(directory/'plan.json', c.make_plan(root))
    original_cache = c.CampaignCache
    monkeypatch.setattr(c, 'CampaignCache',
                        lambda path, plan: original_cache(path, plan, root=root))
    monkeypatch.setattr(c, 'load_plan', lambda path: plan)
    monkeypatch.setattr(c, 'require_capsule', lambda plan: None)
    monkeypatch.setattr(batch, 'ProcessPoolExecutor', ThreadPoolExecutor)
    calls = []

    def provider(model, path, spec):
        calls.append((c.path_identity(path), spec.identity()))
        raw = synthetic_packet(model, path, spec)
        raw['numerics'] = native_numbers(path, spec)
        return raw

    monkeypatch.setattr(c, 'native_provider', provider)
    return directory, plan, calls


def test_resume_completes_all_paths_and_preserves_old_records(campaign):
    directory, plan, calls = campaign
    # Simulate interruption after one successful native calculation.
    first = c.execute_job(str(directory), (0, 11, 0), 0)
    before = Path(first['path']).read_bytes()
    result = batch.execute(directory, 0, 11, directory/'batch-first.json', workers=2)
    assert len(calls) == 30
    assert result['initial_cache_hits'] == 1
    assert result['new_solve_count'] == 29
    assert result['job_count'] == 30
    assert len(result['paths']) == 6 and result['all_paths_passed']
    assert not result['thermal_ensemble_converged']
    assert not result['physical_optical_prediction']
    assert Path(first['path']).read_bytes() == before
    preserved = {p: p.read_bytes() for p in directory.glob('path-*.json')}
    resumed = batch.execute(directory, 0, 11, directory/'batch-second.json', workers=2)
    assert len(calls) == 30
    assert resumed['initial_cache_hits'] == 30 and resumed['new_solve_count'] == 0
    assert all(row['existing_report_preserved'] for row in resumed['paths'])
    assert all(p.read_bytes() == raw for p, raw in preserved.items())
    with pytest.raises(FileExistsError):
        batch.execute(directory, 0, 11, directory/'batch-second.json')
    assert len(calls) == 30


def test_corrupt_cache_refused_before_missing_jobs_run(campaign):
    directory, plan, calls = campaign
    first = c.execute_job(str(directory), (0, 11, 4), 0)
    path = Path(first['path'])
    raw = path.read_text(encoding='utf-8')
    path.write_text(raw.replace('gabes-native-thermal-packet-v1', 'changed'), encoding='utf-8')
    with pytest.raises(ValueError):
        batch.execute(directory, 0, 11, directory/'failed.json')
    assert len(calls) == 1
    assert not (directory/'failed.json').exists()
    assert not list(directory.glob('path-*.json'))


def test_existing_path_report_mismatch_is_never_overwritten(campaign):
    directory, plan, calls = campaign
    for level in range(5):
        c.execute_job(str(directory), (0, 11, 0), level)
    old = directory/'path-p0-s11-i0.json'
    c.write_record(old, {'schema': 'unrelated-sealed-report'})
    before = old.read_bytes()
    with pytest.raises(ValueError, match='Existing path report'):
        batch.execute(directory, 0, 11, directory/'failed.json')
    assert len(calls) == 5 and old.read_bytes() == before


@pytest.mark.parametrize('damage', ['empty_jobs', 'cache_hash', 'job_key', 'physical_key', 'source'])
def test_resealed_bad_report_references_fail_resume(campaign, damage):
    directory, plan, calls = campaign
    selection = (0, 11, 0)
    jobs = {(selection, level): c.execute_job(str(directory), selection, level) for level in range(5)}
    cache = c.CampaignCache(directory/'cache', plan)
    batch.publish_path(directory, plan, selection, jobs, cache)
    path = directory/'path-p0-s11-i0.json'
    record = c.records.read_sealed(path)
    if damage == 'empty_jobs':
        record['jobs'] = []
    elif damage == 'cache_hash':
        record['ledger']['cache_records'][0]['record_sha256'] = '0'*64
    elif damage == 'job_key':
        record['jobs'][1]['key'] = '0'*64
    elif damage == 'physical_key':
        record['ledger']['physical_path_key'] = '0'*64
    else:
        record['ledger']['current_path_source'] = 'different-node'
    # Re-sealing makes this a semantic corruption test, not just hash damage.
    record.pop('record_sha256')
    path.unlink()
    c.write_record(path, record)
    before = path.read_bytes()
    with pytest.raises(ValueError, match='Existing'):
        batch.execute(directory, 0, 11, directory/'bad.json')
    assert len(calls) == 5 and path.read_bytes() == before


def test_worker_failure_does_not_dispatch_the_rest_of_the_grid(campaign, monkeypatch):
    directory, plan, calls = campaign
    attempted = []
    def fail(directory, selection, level):
        attempted.append((selection, level))
        raise RuntimeError('deliberate worker failure')
    monkeypatch.setattr(c, 'execute_job', fail)
    with pytest.raises(RuntimeError, match='worker failure'):
        batch.execute(directory, 0, 11, directory/'failed.json', workers=2)
    assert 1 <= len(attempted) <= 2
    assert not (directory/'failed.json').exists() and not calls


def test_failed_path_retains_failed_evidence_not_success(campaign, monkeypatch):
    directory, plan, calls = campaign
    provider = c.native_provider

    def failed_provider(model, path, spec):
        packet = provider(model, path, spec)
        packet['audit']['passed'] = False
        return packet

    monkeypatch.setattr(c, 'native_provider', failed_provider)
    result = batch.execute(directory, 0, 11, directory/'failed-paths.json', workers=2)
    assert len(calls) == 30 and result['job_count'] == 30
    assert not result['all_paths_passed']
    assert all(not row['path_passed'] for row in result['paths'])
    assert not result['thermal_ensemble_converged']


@pytest.mark.parametrize('grid', [(False, 11), (1, True), (9, 11), (1, 12)])
def test_undeclared_grid_rejected(campaign, grid):
    directory, plan, calls = campaign
    with pytest.raises(ValueError, match='declared grid'):
        batch.execute(directory, *grid, directory/'bad.json')
    assert not calls


@pytest.mark.parametrize('damage', ['truncated', 'empty', 'occupancy'])
def test_declared_grid_cannot_drop_uncomputed_paths(campaign, damage):
    directory, plan, calls = campaign
    grid = next(row for row in plan['grids'] if (row['power'], row['seed']) == (0, 11))
    if damage == 'occupancy':
        grid['mean_occupancy'] *= 2
    else:
        grid['paths'] = grid['paths'][:-1] if damage == 'truncated' else []
    with pytest.raises(ValueError, match='grid count or occupancy'):
        batch.execute(directory, 0, 11, directory/'bad.json')
    assert not calls


def test_report_cannot_be_written_inside_native_cache(campaign):
    directory, plan, calls = campaign
    with pytest.raises(ValueError, match='outside.*cache'):
        batch.execute(directory, 0, 11, directory/'cache'/'bad.json')
    assert not calls


@pytest.mark.parametrize('workers', [0, 5, True, 1.5])
def test_invalid_process_count_rejected_before_loading_plan(tmp_path, workers):
    with pytest.raises(ValueError, match='worker processes'):
        batch.execute(tmp_path, 0, 11, tmp_path/'out.json', workers)


def test_live_interpreter_cannot_dispatch_batch(campaign, monkeypatch):
    directory, plan, calls = campaign
    def forbidden(plan):
        raise ValueError('fresh captured-source interpreter required')
    monkeypatch.setattr(c, 'require_capsule', forbidden)
    with pytest.raises(ValueError, match='captured-source'):
        batch.execute(directory, 0, 11, directory/'bad.json')
    assert not calls


def test_controller_archive_hash_uses_the_same_single_read_as_contents(tmp_path, monkeypatch):
    original = Path.read_bytes
    captured = []
    def changing_read(path):
        if path == Path(batch.__file__):
            captured.append(b'first version' if not captured else b'concurrent edit')
            return captured[-1]
        return original(path)
    monkeypatch.setattr(Path, 'read_bytes', changing_read)
    control = batch.archive_controller(tmp_path)
    assert captured == [b'first version']
    archive = tmp_path/'controllers'/('thermal_campaign_batch-'+control['raw_sha256']+'.py')
    assert archive.read_bytes() == captured[0]
    assert hashlib.sha256(archive.read_bytes()).hexdigest() == control['raw_sha256']


def test_launcher_copies_controller_outside_frozen_inventory(tmp_path, monkeypatch):
    directory = tmp_path/'campaign'
    root = tmp_path/'repo'
    root.mkdir()
    monkeypatch.setattr(batch, 'ROOT', root)
    seen = []
    monkeypatch.setenv('PYTHONPATH', 'unrelated-source-tree')
    monkeypatch.setattr(c, 'extract_bundle', lambda directory, folder: seen.append(Path(folder)))

    def run(command, *, cwd, env, check):
        assert Path(command[1]).read_bytes() == Path(batch.__file__).read_bytes()
        assert Path(command[1]).parent == Path(cwd) == seen[0]
        assert 'PYTHONPATH' not in env and env['OPENBLAS_NUM_THREADS'] == '1'
        assert check and '--execute' in command

    monkeypatch.setattr(batch.subprocess, 'run', run)
    batch.launch(directory, 1, 11, tmp_path/'out.json')
    assert not seen[0].exists()
