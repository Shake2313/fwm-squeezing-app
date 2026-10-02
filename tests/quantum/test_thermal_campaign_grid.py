"""Sidecar contracts with synthetic packets; no expensive native solves."""

import hashlib
import json
from pathlib import Path
import subprocess

import numpy as np
import pytest

from tools import thermal_campaign_grid as g


c = g.c
REAL_BYTE_GUARD = g.CacheByteGuard


def install_fixture_guard(patch):
    """Only in-memory test doubles lack persisted native records to bind."""
    class MemoryGuard:
        def __init__(self, cache):
            self.cache = cache

        def __getattr__(self, name):
            return getattr(self.cache, name)

        def verify(self):
            return []

    patch.setattr(g, 'CacheByteGuard', lambda cache:
                  REAL_BYTE_GUARD(cache) if hasattr(cache, 'directory') else MemoryGuard(cache))


def synthetic_packet(model, path, spec):
    """Finite Gram matrices, nonzero complex mean/response, and native metadata."""
    problem = model.problem(path)
    nf = len(model.analysis_axis.omega_rad_s)
    ns = len(model.convention().source_names)
    rng = np.random.default_rng(734)
    amplitude = rng.normal(size=(ns, nf, 4, 4))+1j*rng.normal(size=(ns, nf, 4, 4))
    factor = 1+.08*np.sin(path.entry_position_m@[700., 900., 300.])
    greater = path.residence_time_s**2*factor*(amplitude@amplitude.conj().swapaxes(-2, -1))/100
    mean = path.residence_time_s*factor*np.tile([.1+.2j, .3-.1j, .1-.2j, .3+.1j], (nf, 1))
    return {'greater': greater.sum(axis=0), 'lesser': .6*greater.sum(axis=0),
        'greater_by_source': greater, 'lesser_by_source': .6*greater,
        'mean_pulse': mean, 'retarded_response': path.residence_time_s**2*amplitude[0],
        'exit_state': model.boundary_state, 'metadata': problem['metadata'],
        'numerics': {'synthetic_test_fixture': True}, 'audit': {'passed': True}}


class FixtureCache:
    """Read-only test double, never accepted by the real CampaignCache writer."""

    def __init__(self, plan):
        self.plan, self.calls, self.packets = plan, [], {}
        self.model = c.r.default_model()
        self.numerical = c.pilot.step_plan(plan['max_steps_s'])
        inflow = self.model.inflow(1, 11)
        for index in range(len(inflow.rate_s_inverse)):
            path = inflow.path(index)
            for spec in self.numerical.primary+self.numerical.reference:
                self.packets[self.key(self.model, path, spec)] = synthetic_packet(self.model, path, spec)

    def identity(self, model, path, spec):
        return {'model': model.identity(), 'path': c.path_identity(path), 'spec': spec.identity()}

    def key(self, model, path, spec):
        return c.r.digest(self.identity(model, path, spec))

    def get(self, model, path, spec, provider=None):
        assert provider is None, 'aggregation may never invoke a solve'
        key = self.key(model, path, spec)
        self.calls.append((path.source, key))
        if key not in self.packets:
            raise FileNotFoundError('synthetic fixture packet missing')
        packet = c.r._validated_packet(model, path, self.packets[key])
        token = c.r.digest(packet)
        return packet, {'key': key, 'hit': True, 'record_sha256': token,
                        'payload_digest': token, 'path': 'synthetic-fixture:'+key}

    def _check_sources(self):
        return self.plan['source_manifest']


class DiskFixtureCache(FixtureCache):
    """Sealed disk test packets; real native validators reject their schema."""

    def __init__(self, plan, directory, packets=None):
        super().__init__(plan)
        self.directory = Path(directory)
        self.directory.mkdir(exist_ok=True)
        if packets is not None:
            self.packets = packets
        for key in self.packets:
            self.persist(key)

    def persist(self, key):
        value = c.records.sealed({'synthetic_test_fixture': True, 'packet': c.r.encode(self.packets[key])})
        (self.directory/(key+'.json')).write_text(json.dumps(value, allow_nan=False), encoding='utf-8')

    def get(self, model, path, spec, provider=None):
        assert provider is None
        key = self.key(model, path, spec)
        self.calls.append((path.source, key))
        filename = self.directory/(key+'.json')
        record = c.records.read_sealed(filename)
        packet = c.r._validated_packet(model, path, c.r.decode(record['packet']))
        return packet, {'key': key, 'hit': True, 'record_sha256': record['record_sha256'],
                        'payload_digest': c.r.digest(packet), 'path': str(filename)}


@pytest.fixture
def campaign(tmp_path, monkeypatch):
    root = tmp_path/'minimal-source'
    (root/'gabes').mkdir(parents=True)
    (root/'gabes/fixture.py').write_text('value = 1\n', encoding='utf-8')
    directory = tmp_path/'campaign'
    c.create_campaign(directory, root=root)
    plan = c.records.read_sealed(directory/'plan.json')
    cache = FixtureCache(plan)
    monkeypatch.setattr(c, 'load_plan', lambda *args, **kwargs: plan)
    monkeypatch.setattr(c, 'require_capsule', lambda plan: None)
    monkeypatch.setattr(c, 'CampaignCache', lambda *args, **kwargs: cache)
    install_fixture_guard(monkeypatch)
    return directory, plan, cache, tmp_path/'grid.json'


def test_complete_grid_is_sealed_read_only_and_cannot_claim_ensemble_convergence(campaign):
    directory, plan, cache, output = campaign
    original_keys = set(cache.packets)
    report = g.execute(directory, 1, 11, output)
    assert report == c.records.read_sealed(output)
    assert report['grid_passed'] and report['path_evidence_passed']
    assert not report['thermal_ensemble_converged'] and not report['physical_optical_prediction']
    assert report['row']['candidate']['certified'] is False
    assert report['row']['cache_hits'] == 60 and report['row']['cache_misses'] == 0
    assert len(cache.calls) == 72  # 5 comparisons + one independent read per path
    assert set(cache.packets) == original_keys
    assert max(report['aggregation_audit']['errors'].values()) < 1e-13
    assert report['convergence_gate']['reasons'] and not report['convergence_gate']['comparisons']
    assert report['campaign_sha256'] == plan['record_sha256']
    assert report['source_identity'] == plan['source_identity']
    assert report['source_bundle_sha256'] == plan['source_bundle_sha256']
    assert report['controller']['raw_sha256'] == hashlib.sha256(Path(g.__file__).read_bytes()).hexdigest()
    assert (directory/report['controller']['archive']).read_bytes() == Path(g.__file__).read_bytes()
    coverage = report['coverage_inventory']
    assert (coverage['declared_grid_count'], coverage['declared_path_entries'],
            coverage['unique_physical_paths'], coverage['declared_unique_jobs']) == (9, 126, 72, 360)
    assert coverage['revalidated_current_grid_jobs'] == 60
    # Synthetic in-memory packets are deliberately not masqueraded as persisted
    # native files. Inventory presence and validated reads are separate facts.
    assert coverage['present_job_files'] == 0 and coverage['missing_job_files'] == 360
    assert coverage['file_presence_is_validated_evidence'] is False


def test_nested_cache_reuse_rebinds_path_sources_and_keeps_current_rates(campaign):
    directory, plan, cache, output = campaign
    fine = c.r.decode(g.execute(directory, 1, 11, output))
    coarse = c.r.decode(g.execute(directory, 0, 11, output.with_name('coarse.json')))
    assert coarse['grid_passed'] and coarse['row']['cache_hits'] == 30
    for face in range(6):
        lo = coarse['row']['path_ledgers'][face]
        hi = fine['row']['path_ledgers'][2*face]
        assert lo['cache_records'][2]['key'] == hi['cache_records'][2]['key']
        if face:
            assert lo['current_path_source'] != hi['current_path_source']
            assert lo['target_digest'] != hi['target_digest']
    assert not np.allclose(coarse['row']['candidate']['spectra']['greater'],
                           fine['row']['candidate']['spectra']['greater'], atol=0)
    inflow = cache.model.inflow(0, 11)
    literal_number = np.zeros_like(coarse['row']['candidate']['spectra']['poisson_number'])
    mask = np.equal.outer([1, -1, -1, 1], [1, -1, -1, 1])
    for index, rate in enumerate(inflow.rate_s_inverse):
        packet = synthetic_packet(cache.model, inflow.path(index), cache.numerical.primary[-1])
        mean = packet['mean_pulse']
        literal_number += rate*np.einsum('fi,fj->fij', mean, mean.conj())*mask
    np.testing.assert_allclose(coarse['row']['candidate']['spectra']['poisson_number'], literal_number, rtol=1e-14)
    assert np.linalg.norm(literal_number) > 0
    gate = c.r.convergence_gate([coarse['row'], fine['row']], powers=plan['powers'], seeds=plan['seeds'])
    assert not gate['passed']  # Same scramble, only two grids, no fabricated rows.


@pytest.mark.parametrize('metric', c.r.stream.STREAM_METRICS)
def test_every_rf_source_number_and_response_output_controls_direct_check(campaign, metric):
    _, _, cache, _ = campaign
    _, row = c.r.run_grid(cache.model, 0, 11, cache.numerical, cache)
    row = c.r.decode(c.r.encode(row))
    changed = row['candidate']['spectra'][metric]
    changed[-1, -1] *= 1.01 if metric.endswith('_by_source') else 1.2
    audit = g.direct_weighted_audit(cache.model, cache.model.inflow(0, 11), cache.numerical, cache, row)
    assert audit['performed'] and not audit['passed']
    assert audit['errors'][metric] > g.AGGREGATION_BUDGET


@pytest.mark.parametrize('damage', ['missing', 'primary', 'independent', 'local_audit'])
def test_incomplete_or_failed_five_calculation_evidence_withholds_entire_grid(campaign, damage):
    directory, _, cache, output = campaign
    path = cache.model.inflow(0, 11).path(2)
    spec = cache.numerical.reference[-1] if damage == 'independent' else cache.numerical.primary[0]
    key = cache.key(cache.model, path, spec)
    if damage == 'missing':
        del cache.packets[key]
    elif damage == 'local_audit':
        cache.packets[key]['audit']['passed'] = False
    else:
        cache.packets[key]['retarded_response'] *= 1.1
    report = g.execute(directory, 0, 11, output)
    assert not report['grid_passed']
    assert report['row']['candidate']['spectra'] is None
    assert report['row']['candidate']['failed_path_index'] == 2
    assert not report['aggregation_audit']['performed']
    assert not report['thermal_ensemble_converged']


def test_changed_packet_between_path_gate_and_independent_read_cannot_pass(campaign):
    _, _, cache, _ = campaign
    _, row = c.r.run_grid(cache.model, 0, 11, cache.numerical, cache)
    key = cache.key(cache.model, cache.model.inflow(0, 11).path(0), cache.numerical.primary[-1])
    cache.packets[key]['retarded_response'] *= 1.0001
    with pytest.raises(ValueError, match='changed between'):
        g.direct_weighted_audit(cache.model, cache.model.inflow(0, 11), cache.numerical, cache, row)


def test_direct_audit_failure_withholds_spectrum_even_after_passing_path_gates(campaign, monkeypatch):
    directory, _, _, output = campaign
    run = c.r.run_grid
    def damaged(*args, **kwargs):
        candidate, row = run(*args, **kwargs)
        row = c.r.decode(c.r.encode(row))
        row['candidate']['spectra']['poisson_number'] *= 2
        return row['candidate'], row
    monkeypatch.setattr(c.r, 'run_grid', damaged)
    report = g.execute(directory, 0, 11, output)
    assert report['path_evidence_passed']
    assert not report['grid_passed'] and report['row']['candidate']['spectra'] is None
    assert report['withheld_candidate_digest'] is not None
    assert not report['convergence_gate']['passed']


@pytest.mark.parametrize('damage', ['count', 'rate', 'path', 'index', 'occupancy', 'duplicate'])
def test_selected_grid_declaration_is_checked_before_any_cache_read(campaign, damage):
    directory, plan, cache, output = campaign
    row = plan['grids'][0]
    if damage == 'count':
        row['paths'].pop()
    elif damage == 'rate':
        row['paths'][-1]['rate_s_inverse'] *= .5
    elif damage == 'path':
        row['paths'][-1]['physical_path']['velocity_m_s'][0] += 1
    elif damage == 'index':
        row['paths'][-1]['index'] = True
    elif damage == 'occupancy':
        row['mean_occupancy'] *= 2
    else:
        plan['grids'].append(row)
    # Keep the stored seal consistent so this still exercises the semantic
    # grid-declaration checks beyond the new immutable-input byte binding.
    resealed = c.records.sealed({name: value for name, value in plan.items() if name != 'record_sha256'})
    plan.clear()
    plan.update(resealed)
    (directory/'plan.json').write_text(json.dumps(plan), encoding='utf-8')
    with pytest.raises(ValueError, match='[Dd]eclared'):
        g.execute(directory, 0, 11, output)
    assert cache.calls == [] and not output.exists()


@pytest.mark.parametrize('power,seed', [(False, 11), (1., 11), (3, 11), (0, 99)])
def test_undeclared_or_noninteger_selection_rejected(campaign, power, seed):
    directory, _, cache, output = campaign
    with pytest.raises(ValueError, match='declared'):
        g.execute(directory, power, seed, output)
    assert not cache.calls and not output.exists()


def test_existing_outputs_are_never_replaced_and_reports_do_not_pollute_cache(campaign, monkeypatch):
    directory, _, cache, output = campaign
    output.write_bytes(b'prior user report')
    monkeypatch.setattr(c, 'load_plan', lambda *args: pytest.fail('preflight must reject first'))
    for action in (g.execute, g.launch):
        with pytest.raises(FileExistsError):
            action(directory, 0, 11, output)
    assert output.read_bytes() == b'prior user report'
    with pytest.raises(ValueError, match='outside'):
        g.execute(directory, 0, 11, directory/'cache/new.json')
    assert not cache.calls


def test_wrong_controller_hash_cannot_publish(campaign):
    directory, _, cache, output = campaign
    with pytest.raises(ValueError, match='controller bytes'):
        g.execute(directory, 0, 11, output, expected_controller_sha256='0'*64)
    assert not output.exists() and not cache.calls


def test_matching_controller_archive_preserved_and_corrupt_archive_never_replaced(tmp_path):
    control = g.archive_controller(tmp_path)
    archive = tmp_path/control['archive']
    before = archive.stat().st_mtime_ns
    assert g.archive_controller(tmp_path) == control
    assert archive.stat().st_mtime_ns == before
    archive.write_bytes(b'corrupt captured controller')
    with pytest.raises(ValueError, match='archive changed'):
        g.archive_controller(tmp_path)
    assert archive.read_bytes() == b'corrupt captured controller'


def test_fresh_interpreter_launch_keeps_frozen_inventory_and_missing_cache_read_only(tmp_path, monkeypatch):
    # A new temporary campaign contains current source files, but no numerical
    # result. The real captured subprocess must report the missing first path;
    # no solver is invoked, and no historical evidence fixture is required.
    directory = tmp_path/'campaign'
    c.create_campaign(directory, root=c.ROOT)
    plan = c.records.read_sealed(directory/'plan.json')
    before = {name: (directory/name).read_bytes() for name in ('plan.json', 'sources.zip')}
    monkeypatch.setattr(g, 'ROOT', tmp_path)
    monkeypatch.setenv('PYTHONPATH', str(tmp_path/'unrelated-code'))
    output = tmp_path/'missing.json'
    with pytest.raises(subprocess.CalledProcessError) as exc:
        g.launch(directory, 0, 11, output)
    assert exc.value.returncode == 1
    report = c.records.read_sealed(output)
    assert report['source_identity'] == plan['source_identity']
    assert report['row']['candidate']['spectra'] is None and not report['grid_passed']
    assert report['row']['candidate']['failed_path_index'] == 0
    assert report['row']['cache_misses'] == 0
    assert list((directory/'cache').iterdir()) == []
    assert list((tmp_path/'.git/grand-challenge-runs').iterdir()) == []
    for name, raw in before.items():
        assert (directory/name).read_bytes() == raw


def test_corrupt_real_cache_record_is_reported_without_repair(tmp_path, monkeypatch):
    root = tmp_path/'source'
    (root/'gabes').mkdir(parents=True)
    (root/'gabes/fixture.py').write_text('x = 1\n', encoding='utf-8')
    directory = tmp_path/'campaign'
    c.create_campaign(directory, root=root)
    plan = c.records.read_sealed(directory/'plan.json')
    real_cache = c.CampaignCache(directory/'cache', plan, root=root)
    model = c.r.default_model()
    numerical = c.pilot.step_plan(plan['max_steps_s'])
    key = real_cache.key(model, model.inflow(0, 11).path(0), numerical.primary[0])
    filename = directory/'cache'/(key+'.json')
    filename.write_text('{"record_sha256":"invalid"}', encoding='utf-8')
    raw = filename.read_bytes()
    monkeypatch.setattr(c, 'load_plan', lambda *args: plan)
    monkeypatch.setattr(c, 'require_capsule', lambda plan: None)
    monkeypatch.setattr(c, 'CampaignCache', lambda *args, **kwargs: real_cache)
    report = g.execute(directory, 0, 11, tmp_path/'corrupt.json')
    assert not report['grid_passed'] and report['row']['candidate']['spectra'] is None
    assert 'hash mismatch' in report['row']['candidate']['reasons'][0]
    assert filename.read_bytes() == raw


@pytest.mark.parametrize('artifact', ['packet', 'plan', 'zip', 'controller'])
def test_final_grid_byte_checks_refuse_inputs_changed_after_numerical_reads(campaign, monkeypatch, artifact):
    directory, plan, original_cache, output = campaign
    cache = DiskFixtureCache(plan, directory/'cache', original_cache.packets)
    monkeypatch.setattr(c, 'CampaignCache', lambda *args, **kwargs: cache)
    original = c.r.convergence_gate

    def mutate(*args, **kwargs):
        result = original(*args, **kwargs)
        key = cache.key(cache.model, cache.model.inflow(0, 11).path(0), cache.numerical.reference[0])
        filename = {'packet': cache.directory/(key+'.json'), 'plan': directory/'plan.json',
                    'zip': directory/'sources.zip', 'controller': directory/g.controller_identity()['archive']}[artifact]
        filename.write_bytes(filename.read_bytes()+b' ')
        return result

    monkeypatch.setattr(c.r, 'convergence_gate', mutate)
    with pytest.raises(ValueError, match='bytes changed'):
        g.execute(directory, 0, 11, output)
    assert not output.exists()


def test_grid_archive_capture_binds_original_controller_bytes(campaign, monkeypatch):
    directory, _, cache, output = campaign
    original = g.archive_controller

    def mutate(folder):
        identity = original(folder)
        (Path(folder)/identity['archive']).write_bytes(b'changed before initial archive snapshot')
        return identity

    monkeypatch.setattr(g, 'archive_controller', mutate)
    with pytest.raises(ValueError, match='Archived controller bytes differ'):
        g.execute(directory, 0, 11, output)
    assert not output.exists() and not cache.calls


def test_cache_snapshot_is_bound_to_native_seal_even_when_bytes_revert(tmp_path):
    key = 'a'*64
    filename = tmp_path/(key+'.json')
    saved = c.records.sealed({'value': 1})
    other = c.records.sealed({'value': 2})
    filename.write_text(json.dumps(saved), encoding='utf-8')

    class Cache:
        directory = tmp_path

        def key(self, *args):
            return key

        def get(self, *args, provider=None):
            assert provider is None
            # A different sealed record was validated between the guard's
            # two reads, even though the bytes then reverted to the first one.
            return {}, {'key': key, 'hit': True, 'record_sha256': other['record_sha256']}

    with pytest.raises(ValueError, match='differs from consumed cache byte snapshot'):
        REAL_BYTE_GUARD(Cache()).get(None, None, None)
