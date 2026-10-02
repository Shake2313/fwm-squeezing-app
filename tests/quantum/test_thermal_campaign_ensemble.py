"""Hermetic ensemble audit contracts; synthetic caches never native evidence."""

from copy import deepcopy
import hashlib
import itertools
import json
from pathlib import Path
import shutil
import subprocess
import sys
from types import SimpleNamespace

import numpy as np
import pytest

from tools import thermal_campaign_ensemble as e
from test_thermal_campaign_grid import FixtureCache, DiskFixtureCache, install_fixture_guard, synthetic_packet
from test_thermal_campaign_compare import rewrite_report

c, g = e.c, e.g


def constant_packet(model, path, spec):
    """Constant integrands allow a real 3x3 gate success without atomic solves."""
    packet = synthetic_packet(model, path, spec)
    factor = 1+.08*np.sin(path.entry_position_m@[700., 900., 300.])
    for name in ('greater', 'lesser', 'greater_by_source', 'lesser_by_source'):
        packet[name] *= 1e-18/(path.residence_time_s**2*factor)
    packet['mean_pulse'] *= 1e-9/(path.residence_time_s*factor)
    packet['retarded_response'] *= 1e-18/path.residence_time_s**2
    return packet


@pytest.fixture(scope='module')
def base_campaign(tmp_path_factory):
    folder = tmp_path_factory.mktemp('ensemble-fixture')
    directory = folder/'campaign'
    c.create_campaign(directory, root=c.ROOT)
    plan = c.records.read_sealed(directory/'plan.json')
    cache = FixtureCache(plan)
    cache.packets.clear()
    for power, seed in itertools.product(plan['powers'], plan['seeds']):
        inflow = cache.model.inflow(power, seed)
        for index in range(len(inflow.rate_s_inverse)):
            path = inflow.path(index)
            for spec in cache.numerical.primary+cache.numerical.reference:
                key = cache.key(cache.model, path, spec)
                if key not in cache.packets:
                    cache.packets[key] = constant_packet(cache.model, path, spec)
    with pytest.MonkeyPatch.context() as patch:
        patch.setattr(c, 'load_plan', lambda *args, **kwargs: plan)
        patch.setattr(c, 'require_capsule', lambda *args: None)
        patch.setattr(c, 'CampaignCache', lambda *args, **kwargs: cache)
        install_fixture_guard(patch)
        patch.setattr(c, 'native_provider', lambda *args: pytest.fail('fixture must never solve'))
        for power, seed in itertools.product(plan['powers'], plan['seeds']):
            g.execute(directory, power, seed, folder/f'grid-p{power}-s{seed}.json')
    return folder, cache


@pytest.fixture
def campaign(base_campaign, tmp_path, monkeypatch):
    folder, original = base_campaign
    directory = tmp_path/'campaign'
    shutil.copytree(folder/'campaign', directory)
    for source in folder.glob('grid-*.json'):
        shutil.copyfile(source, tmp_path/source.name)
    plan = c.records.read_sealed(directory/'plan.json')
    cache = object.__new__(FixtureCache)
    cache.plan, cache.calls = plan, []
    cache.model, cache.numerical = original.model, original.numerical
    cache.packets = deepcopy(original.packets)
    monkeypatch.setattr(c, 'load_plan', lambda *args, **kwargs: plan)
    monkeypatch.setattr(c, 'require_capsule', lambda *args: None)
    monkeypatch.setattr(c, 'CampaignCache', lambda *args, **kwargs: cache)
    install_fixture_guard(monkeypatch)
    monkeypatch.setattr(c, 'native_provider', lambda *args: pytest.fail('ensemble must never solve'))
    paths = [tmp_path/f'grid-p{power}-s11.json' for power in plan['powers']]
    return directory, plan, cache, paths, tmp_path/'ensemble.json'


def test_three_power_one_seed_subset_is_valid_but_cannot_certify_thermal_convergence(campaign):
    directory, plan, cache, paths, output = campaign
    before = {path: path.read_bytes() for path in [*paths, directory/'plan.json', directory/'sources.zip']}
    report = c.r.decode(e.execute(directory, list(reversed(paths)), output))
    assert report['audit_passed'] and not report['thermal_ensemble_converged']
    assert not report['physical_optical_prediction'] and report['new_solve_count'] == 0
    assert report['required_powers'] == plan['powers'] and report['required_seeds'] == plan['seeds']
    assert report['coverage']['provided_grids'] == [[0, 11], [1, 11], [2, 11]]
    assert report['coverage']['missing_grids'] == [[power, seed] for power in (0, 1, 2) for seed in (211, 811)]
    assert report['coverage']['required_grid_count'] == 9 and not report['coverage']['complete']
    assert not report['convergence_gate']['passed'] and report['convergence_gate']['comparisons'] == []
    diagnostics = report['pairwise_diagnostics']
    assert diagnostics['diagnostic_only'] and not diagnostics['independent_scrambles']
    assert [(r['reference_grid'], r['candidate_grid']) for r in diagnostics['refinements']] == [
        ([0, 11], [1, 11]), ([1, 11], [2, 11])]
    for item in diagnostics['refinements']:
        current, previous = report['rows'][item['candidate_grid'][0]], report['rows'][item['reference_grid'][0]]
        assert item['errors'] == c.r.stream_errors(current['candidate']['spectra'],
            previous['candidate']['spectra'], current['comparison_scales'])
    assert len(cache.calls) == 6*(6+12+24)
    assert report['loaded_source_provenance']['numerical_modules']
    for identity in report['controllers']:
        raw = (directory/identity['archive']).read_bytes()
        assert hashlib.sha256(raw).hexdigest() == identity['raw_sha256']
    assert {item['name'] for item in report['controllers']} == {
        'thermal_campaign_ensemble.py', 'thermal_campaign_compare.py', 'thermal_campaign_grid.py'}
    for path, raw in before.items():
        assert path.read_bytes() == raw


def test_complete_same_plan_three_by_three_can_pass_the_real_gate(campaign):
    directory, plan, _, _, output = campaign
    paths = sorted(output.parent.glob('grid-*.json'), reverse=True)
    report = c.r.decode(e.execute(directory, paths, output))
    assert report['audit_passed'] and report['coverage']['complete']
    assert report['coverage']['missing_grids'] == []
    assert report['thermal_ensemble_converged'] and not report['physical_optical_prediction']
    expected = c.r.convergence_gate(report['rows'], powers=plan['powers'], seeds=plan['seeds'])
    assert report['convergence_gate'] == expected and len(expected['comparisons']) == 18
    assert len(report['pairwise_diagnostics']['refinements']) == 6
    assert len(report['pairwise_diagnostics']['independent_scrambles']) == 18


def test_failed_available_edge_is_reported_while_missing_seeds_withhold_formal_gate(campaign):
    directory, _, cache, paths, output = campaign
    inflow = cache.model.inflow(1, 11)
    # Change the actual unweighted packets of p1's six added physical paths,
    # equally at all five resolutions. Path gates remain valid; the integrated
    # response changes from R0 to 1.5*R0 at p1 and 1.25*R0 at p2. No saved
    # spectrum, gate flag, or comparison error is fabricated for this test.
    for index in (1, 3, 5, 7, 9, 11):
        path = inflow.path(index)
        for spec in cache.numerical.primary+cache.numerical.reference:
            cache.packets[cache.key(cache.model, path, spec)]['retarded_response'] *= 2
    fresh = [paths[0]]
    for power in (1, 2):
        filename = output.parent/f'changed-p{power}.json'
        assert g.execute(directory, power, 11, filename)['grid_passed']
        fresh.append(filename)
    report = e.execute(directory, fresh, output)
    edges = report['pairwise_diagnostics']['refinements']
    assert len(edges) == 2 and all(not edge['passed'] for edge in edges)
    assert edges[0]['errors']['retarded_response'] == pytest.approx(.5, rel=1e-12)
    assert edges[1]['errors']['retarded_response'] == pytest.approx(1/6, rel=1e-12)
    assert report['audit_passed'] and not report['thermal_ensemble_converged']
    assert not report['convergence_gate']['passed'] and report['convergence_gate']['comparisons'] == []
    assert any('seed' in reason for reason in report['convergence_gate']['reasons'])
    assert report['coverage']['missing_grids'] == [[p, s] for p in (0, 1, 2) for s in (211, 811)]


def test_missing_middle_power_does_not_invent_a_consecutive_refinement(campaign):
    directory, _, _, paths, output = campaign
    report = e.execute(directory, [paths[0], paths[2]], output)
    assert report['audit_passed'] and not report['thermal_ensemble_converged']
    assert report['pairwise_diagnostics']['refinements'] == []
    assert [1, 11] in report['coverage']['missing_grids']


def test_subset_with_two_seeds_measures_both_directed_scramble_comparisons(campaign):
    directory, _, _, paths, output = campaign
    second = output.parent/'grid-p1-s211.json'
    report = e.execute(directory, [paths[1], second], output)
    comparisons = report['pairwise_diagnostics']['independent_scrambles']
    assert len(comparisons) == 2
    assert {(tuple(r['candidate_grid']), tuple(r['reference_grid'])) for r in comparisons} == {
        ((1, 11), (1, 211)), ((1, 211), (1, 11))}
    assert not report['thermal_ensemble_converged']


@pytest.mark.parametrize('kind', ['empty', 'missing_file', 'duplicate_file', 'duplicate_grid'])
def test_missing_or_duplicate_inputs_cannot_create_an_audit(campaign, kind):
    directory, _, cache, paths, output = campaign
    if kind == 'empty':
        inputs, error = [], ValueError
    elif kind == 'missing_file':
        inputs, error = [*paths, output.parent/'missing.json'], FileNotFoundError
    elif kind == 'duplicate_file':
        inputs, error = [*paths, paths[0]], ValueError
    else:
        alias = output.parent/'same-grid.json'
        shutil.copyfile(paths[0], alias)
        inputs, error = [*paths, alias], ValueError
    with pytest.raises(error):
        e.execute(directory, inputs, output)
    assert not output.exists() and not cache.calls


@pytest.mark.parametrize('field', ['campaign_sha256', 'source_identity', 'source_bundle_sha256',
    'model_digest', 'plan_digest', 'target_digest', 'current_path_source'])
def test_resealed_misbound_report_is_rejected(campaign, field):
    directory, _, _, paths, output = campaign
    def change(report):
        target = (report['row'] if field in ('model_digest', 'plan_digest')
                  else report['row']['path_ledgers'][0] if field in ('target_digest', 'current_path_source') else report)
        target[field] = '0'*64
    rewrite_report(paths[0], change)
    with pytest.raises(ValueError):
        e.execute(directory, paths, output)
    assert not output.exists()


@pytest.mark.parametrize('field', ['key', 'record_sha256', 'payload_digest'])
@pytest.mark.parametrize('location', ['ledger', 'aggregation'])
def test_forged_cache_metadata_is_not_trusted_even_with_valid_report_seals(campaign, field, location):
    directory, _, _, paths, output = campaign
    def change(report):
        records = (report['row']['path_ledgers'][0]['cache_records'] if location == 'ledger'
                   else report['aggregation_audit']['cache_records'])
        records[0][field] = 'f'*64
    rewrite_report(paths[0], change)
    with pytest.raises(ValueError, match='freshly validated|references differ'):
        e.execute(directory, paths, output)
    assert not output.exists()


@pytest.mark.parametrize('field', ['archive', 'raw_bytes', 'raw_sha256'])
def test_forged_grid_controller_metadata_is_rejected(campaign, field):
    directory, _, cache, paths, output = campaign
    def change(report):
        report['controller'][field] = 1 if field == 'raw_bytes' else 'not-the-archived-controller'
    rewrite_report(paths[0], change)
    with pytest.raises(ValueError, match='controller'):
        e.execute(directory, paths, output)
    assert not output.exists() and not cache.calls


@pytest.mark.parametrize('kind', ['claimed_convergence', 'forged_gate', 'incomplete_ledger', 'altered_spectrum'])
def test_incomplete_or_fabricated_evidence_cannot_claim_convergence(campaign, kind):
    directory, _, _, paths, output = campaign
    def change(report):
        if kind == 'claimed_convergence':
            report['thermal_ensemble_converged'] = True
        elif kind == 'forged_gate':
            report['convergence_gate'] = {'passed': True, 'comparisons': [], 'reasons': []}
        elif kind == 'incomplete_ledger':
            report['row']['path_ledgers'].pop()
        else:
            candidate = report['row']['candidate']
            candidate['spectra']['poisson_number'] *= 1.01
            candidate['candidate_digest'] = c.r.stream._spectrum_digest(candidate)
    rewrite_report(paths[0], change)
    with pytest.raises(ValueError):
        e.execute(directory, paths, output)
    assert not output.exists()


@pytest.mark.parametrize('damage', ['missing', 'changed', 'failed_reference'])
def test_actual_cache_is_revalidated_without_repair_or_native_solve(campaign, damage):
    directory, _, cache, paths, output = campaign
    path = cache.model.inflow(0, 11).path(0)
    spec = cache.numerical.reference[-1]
    key = cache.key(cache.model, path, spec)
    if damage == 'missing':
        del cache.packets[key]
    elif damage == 'failed_reference':
        cache.packets[key]['retarded_response'] *= 1.1
    else:
        for spec in cache.numerical.primary+cache.numerical.reference:
            cache.packets[cache.key(cache.model, path, spec)]['retarded_response'] *= 1.001
    before = set(cache.packets)
    with pytest.raises(ValueError, match='Fresh grid|freshly validated'):
        e.execute(directory, paths, output)
    assert not output.exists() and set(cache.packets) == before


def test_diagnostics_use_current_scales_without_rewriting_the_formal_gate(campaign):
    _, plan, _, paths, _ = campaign
    rows = [c.r.decode(c.records.read_sealed(path))['row'] for path in paths[:2]]
    # Dark test response isolates the SI floor; different comparison scales
    # make candidate-vs-coarse normalization observable without a native solve.
    rows[0]['candidate']['spectra']['retarded_response'][:] = 0
    rows[1]['candidate']['spectra']['retarded_response'][:] = 1e-30
    rows[0]['comparison_scales']['retarded_response'] = 1e-12
    rows[1]['comparison_scales']['retarded_response'] = 1e-9
    diagnostics = e.pairwise_diagnostics(rows, plan)['refinements'][0]
    expected = c.r.stream_errors(rows[1]['candidate']['spectra'], rows[0]['candidate']['spectra'], rows[1]['comparison_scales'])
    other = c.r.stream_errors(rows[1]['candidate']['spectra'], rows[0]['candidate']['spectra'], rows[0]['comparison_scales'])
    assert diagnostics['errors'] == expected
    assert expected['retarded_response'] != other['retarded_response']


def test_loaded_source_outside_capsule_or_manifest_is_rejected(campaign, tmp_path, monkeypatch):
    directory, _, _, paths, output = campaign
    foreign = tmp_path/'foreign.py'
    foreign.write_text('value = 1\n', encoding='utf-8')
    monkeypatch.setitem(sys.modules, 'gabes.injected_test_module', SimpleNamespace(__file__=str(foreign)))
    with pytest.raises(ValueError, match='outside the captured source capsule'):
        e.execute(directory, paths, output)
    assert not output.exists()


def test_existing_output_and_changed_controller_archives_are_never_overwritten(campaign):
    directory, _, cache, paths, output = campaign
    output.write_bytes(b'previous audit')
    for action in (e.execute, e.launch):
        with pytest.raises(FileExistsError):
            action(directory, paths, output)
    assert output.read_bytes() == b'previous audit' and not cache.calls
    captured = e.captured_controllers()
    e.comparison.archive_controllers(directory, captured)
    archive = directory/captured[0][0]['archive']
    before = archive.stat().st_mtime_ns
    e.comparison.archive_controllers(directory, captured)
    assert archive.stat().st_mtime_ns == before
    archive.write_bytes(b'changed controller')
    with pytest.raises(ValueError, match='archive changed'):
        e.execute(directory, paths, output.with_name('new.json'))
    assert archive.read_bytes() == b'changed controller'


def test_controller_hash_request_must_bind_all_three_dependencies(campaign):
    directory, _, cache, paths, output = campaign
    with pytest.raises(ValueError, match='launch request'):
        e.execute(directory, paths, output, expected_controller_hashes=['0'*64]*3)
    assert not output.exists() and not cache.calls


@pytest.mark.parametrize('artifact', ['packet', 'plan', 'zip', 'controller', 'comparison_helper', 'grid_helper'])
def test_ensemble_refuses_inputs_changed_after_numerical_gate(campaign, monkeypatch, artifact):
    directory, plan, original_cache, _, output = campaign
    cache = DiskFixtureCache(plan, directory/'cache', original_cache.packets)
    monkeypatch.setattr(c, 'CampaignCache', lambda *args, **kwargs: cache)
    paths = []
    for power in plan['powers']:
        filename = output.with_name(f'disk-grid-p{power}.json')
        g.execute(directory, power, 11, filename)
        paths.append(filename)
    original = c.r.convergence_gate

    def mutate(*args, **kwargs):
        result = original(*args, **kwargs)
        key = cache.key(cache.model, cache.model.inflow(0, 11).path(0), cache.numerical.reference[0])
        captured = e.captured_controllers()
        filename = {'packet': cache.directory/(key+'.json'), 'plan': directory/'plan.json',
                    'zip': directory/'sources.zip', 'controller': directory/captured[0][0]['archive'],
                    'comparison_helper': directory/captured[1][0]['archive'],
                    'grid_helper': directory/captured[2][0]['archive']}[artifact]
        filename.write_bytes(filename.read_bytes()+b' ')
        return result

    # The helper also calls the gate while rereading individual grids. Inject
    # mutation only after all requested grids reached the ensemble-level gate.
    def final_gate(rows, **kwargs):
        return mutate(rows, **kwargs) if len(rows) > 1 else original(rows, **kwargs)

    monkeypatch.setattr(c.r, 'convergence_gate', final_gate)
    with pytest.raises(ValueError, match='bytes changed|controller archive differs'):
        e.execute(directory, paths, output)
    assert not output.exists()


@pytest.mark.parametrize('controller_index', [0, 1, 2])
def test_ensemble_archive_capture_binds_original_controller_bytes(campaign, monkeypatch, controller_index):
    directory, _, cache, paths, output = campaign
    original = e.comparison.archive_controllers

    def mutate(folder, captured):
        original(folder, captured)
        identity = captured[controller_index][0]
        (Path(folder)/identity['archive']).write_bytes(b'changed before initial archive snapshot')

    monkeypatch.setattr(e.comparison, 'archive_controllers', mutate)
    with pytest.raises(ValueError, match='Archived controller bytes differ'):
        e.execute(directory, paths, output)
    assert not output.exists() and not cache.calls


def test_fresh_frozen_subprocess_rejects_synthetic_reports_without_native_cache(campaign, tmp_path, monkeypatch, capfd):
    directory, _, _, paths, output = campaign
    original = (directory/'sources.zip').read_bytes()
    monkeypatch.setattr(e, 'ROOT', tmp_path)
    monkeypatch.setenv('PYTHONPATH', str(tmp_path/'unrelated-code'))
    with pytest.raises(subprocess.CalledProcessError):
        e.launch(directory, paths, output)
    stderr = capfd.readouterr().err
    assert 'Fresh grid requires complete passing path evidence' in stderr
    assert 'Native campaign path not computed' in stderr
    assert not output.exists() and list((directory/'cache').iterdir()) == []
    assert (directory/'sources.zip').read_bytes() == original
    assert list((tmp_path/'.git/grand-challenge-runs').iterdir()) == []
    assert len(list((directory/'controllers').glob('thermal_campaign_*.py'))) == 3
