"""Fresh campaign contracts; synthetic storage fixtures never atomic evidence."""

from contextlib import contextmanager

import numpy as np
import pytest

from analysis.grand_challenge import thermal_campaign as c


@pytest.fixture
def source_root(tmp_path):
    root = tmp_path/'source'
    (root/'gabes').mkdir(parents=True)
    (root/'gabes/example.py').write_text('coefficient = 1\n', encoding='utf-8')
    return root


@pytest.fixture
def minimal_plan(source_root):
    paths = c.source_paths(source_root)
    manifest = c.provenance.capture_sources(source_root, paths)
    return {'source_manifest': manifest, 'source_paths': paths,
            'source_identity': c.provenance.source_identity(manifest),
            'environment': c.records.environment()}


def native_numbers(path, spec):
    parameters = c.numerical_parameters(path, spec)
    nums = {'requested_parameters': parameters, 'elapsed_seconds': 0.1,
            'runtime_blas_threads': [{'user_api': 'blas', 'num_threads': 1}],
            **{k: v for k, v in parameters.items() if k != 'sample_count'}}
    if spec.method == c.pilot.PRIMARY_METHOD:
        nums['substeps'] = 2*parameters['segments']
    else:
        nums.update(forward_density_evaluations=1, backward_evaluations=1, density_mesh_points=2)
    return nums


def test_native_paths_are_declared_before_any_solve(source_root, monkeypatch):
    monkeypatch.setattr(c, 'native_provider', lambda *args: pytest.fail('planning must not solve'))
    plan = c.make_plan(source_root)
    assert len(plan['grids']) == 9
    assert sum(len(g['paths']) for g in plan['grids']) == 126
    grid0 = plan['grids'][0]
    grid1 = next(g for g in plan['grids'] if (g['power'], g['seed']) == (1, 11))
    for i in range(6):
        assert grid0['paths'][i]['physical_path'] == grid1['paths'][2*i]['physical_path']
        assert grid0['paths'][i]['rate_s_inverse'] == 2*grid1['paths'][2*i]['rate_s_inverse']
    assert plan['historical_numerical_results_reused'] is False


@pytest.mark.parametrize('power', (0, 1, 2))
def test_nested_mapping_preserves_the_within_face_sobol_index(power):
    model = c.r.default_model()
    a, b = model.inflow(power, 11), model.inflow(power+1, 11)
    face_size = 2**power
    for i in range(len(a.rate_s_inverse)):
        face, offset = divmod(i, face_size)
        j = face*(2*face_size)+offset
        assert c.path_identity(a.path(i)) == c.path_identity(b.path(j))
        assert a.rate_s_inverse[i] == 2*b.rate_s_inverse[j]


@pytest.mark.parametrize('powers,seeds', [((0, 1), (11, 211, 811)),
    ((0, 1, 1), (11, 211, 811)), ((0, 1, 2), (11, 11, 811)),
    ((False, 1, 2), (11, 211, 811)), ((0, 1, 2), (-1, 211, 811))])
def test_invalid_campaign_grid_declarations(source_root, powers, seeds):
    with pytest.raises(ValueError):
        c.make_plan(source_root, powers=powers, seeds=seeds)


@pytest.mark.parametrize('damage', ['rate', 'velocity', 'index'])
def test_native_selection_cannot_reinterpret_a_declared_path(source_root, damage):
    plan = c.make_plan(source_root)
    row = plan['grids'][3]['paths'][1]
    if damage == 'rate':
        row['rate_s_inverse'] *= .5
    elif damage == 'velocity':
        row['physical_path']['velocity_m_s'][0] += 1
    else:
        row['index'] = 0
    with pytest.raises(ValueError, match='boundary path'):
        c.selected_path(plan, 1, 11, 1)


def test_cache_identity_reuses_native_path_but_not_different_equations(source_root, minimal_plan, tmp_path):
    model = c.r.default_model()
    coarse, fine = model.inflow(0, 11).path(2), model.inflow(1, 11).path(4)
    spec = c.pilot.step_plan(c.MAX_STEPS).primary[0]
    cache = c.CampaignCache(tmp_path/'cache', minimal_plan, root=source_root)
    assert cache.key(model, coarse, spec) == cache.key(model, fine, spec)
    assert cache.key(model, coarse, spec) != cache.key(model, model.inflow(1, 11).path(5), spec)
    assert cache.key(model, coarse, spec) != cache.key(model, coarse, c.pilot.step_plan(c.MAX_STEPS).primary[1])
    with pytest.raises(FileNotFoundError):
        cache.get(model, coarse, spec)
    with pytest.raises(ValueError, match='native'):
        cache.get(model, coarse, spec, lambda *args: {})


def test_cache_cannot_dispatch_from_live_interpreter(source_root, minimal_plan, tmp_path, monkeypatch):
    model = c.r.default_model()
    path = model.inflow(1, 11).path(1)
    spec = c.pilot.step_plan(c.MAX_STEPS).primary[0]
    cache = c.CampaignCache(tmp_path/'cache', minimal_plan, root=source_root)
    def forbidden(*args):
        pytest.fail('live provider must be rejected before any solve')
    monkeypatch.setattr(c, 'native_provider', forbidden)
    with pytest.raises(ValueError, match='fresh captured-source interpreter'):
        cache.get(model, path, spec, c.native_provider)


def test_packet_axis_and_complex_arrays_survive_sealed_round_trip(tmp_path):
    axis = c.r.default_model().analysis_axis
    payload = {'packet': {'analysis_axis': axis,
                         'greater': np.eye(4, dtype=complex)*(1+0.0j),
                         'retarded_response': np.array([1+2j])}}
    path = tmp_path/'packet.json'
    record = c.write_record(path, payload)
    restored = c.r.decode(c.records.read_sealed(path))
    assert restored['record_sha256'] == record['record_sha256']
    np.testing.assert_array_equal(restored['packet']['analysis_axis'].omega_rad_s, axis.omega_rad_s)
    np.testing.assert_array_equal(restored['packet']['retarded_response'], np.array([1+2j]))


@pytest.mark.parametrize('failure', [KeyboardInterrupt, OSError])
def test_interrupted_record_write_leaves_no_partial_final_and_can_retry(
        tmp_path, monkeypatch, failure):
    destination = tmp_path/'packet.json'
    payload = {'schema': 'storage-fixture', 'values': list(range(30))}
    original_fdopen = c.os.fdopen
    partial_writes = []

    @contextmanager
    def interrupted_stream(*args, **kwargs):
        with original_fdopen(*args, **kwargs) as handle:
            class PartialWriter:
                def write(self, text):
                    partial_writes.append(handle.write(text[:12]))
                    handle.flush()
                    raise failure('simulated interruption after a partial write')
            yield PartialWriter()

    with monkeypatch.context() as patch:
        patch.setattr(c.os, 'fdopen', interrupted_stream)
        with pytest.raises(failure, match='partial write'):
            c.write_record(destination, payload)
    assert partial_writes == [12]
    assert not destination.exists()
    assert not list(tmp_path.iterdir())
    written = c.write_record(destination, payload)
    assert c.records.read_sealed(destination) == written
    assert list(tmp_path.iterdir()) == [destination]


def test_record_publication_preserves_destination_created_during_write(
        tmp_path, monkeypatch):
    destination = tmp_path/'packet.json'
    existing = b'previous immutable evidence\n'
    original_fsync = c.os.fsync

    def concurrent_publication(descriptor):
        original_fsync(descriptor)
        destination.write_bytes(existing)

    monkeypatch.setattr(c.os, 'fsync', concurrent_publication)
    with pytest.raises(FileExistsError):
        c.write_record(destination, {'schema': 'storage-fixture'})
    assert destination.read_bytes() == existing
    assert list(tmp_path.iterdir()) == [destination]


def test_record_publication_never_replaces_existing_evidence(tmp_path):
    destination = tmp_path/'packet.json'
    original = c.write_record(destination, {'schema': 'original-fixture'})
    original_bytes = destination.read_bytes()
    with pytest.raises(FileExistsError):
        c.write_record(destination, {'schema': 'replacement-fixture'})
    assert destination.read_bytes() == original_bytes
    assert c.records.read_sealed(destination) == original
    assert list(tmp_path.iterdir()) == [destination]


def test_cache_rejects_source_edits_and_new_dependency_inventory(source_root, minimal_plan, tmp_path):
    cache = c.CampaignCache(tmp_path/'cache', minimal_plan, root=source_root)
    p = source_root/'gabes/example.py'
    p.write_text('coefficient = 2\n')
    with pytest.raises(ValueError):
        cache._check_sources()
    p.write_text('coefficient = 1\n')
    (source_root/'gabes/extra.py').write_text('other = 3\n')
    with pytest.raises(ValueError, match='inventory'):
        cache._check_sources()


@pytest.mark.parametrize('level', range(5))
def test_native_parameters_match_actual_solver_work(level):
    path = c.r.default_model().inflow(1, 11).path(1)
    plan = c.pilot.step_plan(c.MAX_STEPS)
    spec = (plan.primary+plan.reference)[level]
    c.validate_native({'numerics': native_numbers(path, spec)}, path, spec)


@pytest.mark.parametrize('damage', ['segments', 'substeps', 'requested', 'threads', 'synthetic', 'elapsed'])
def test_primary_cannot_be_labeled_with_a_different_solve(damage):
    path = c.r.default_model().inflow(1, 11).path(1)
    spec = c.pilot.step_plan(c.MAX_STEPS).primary[0]
    n = native_numbers(path, spec)
    if damage in ('segments', 'substeps'):
        n[damage] -= 1
    elif damage == 'requested':
        n['requested_parameters'] = {'segments': 1}
    elif damage == 'threads':
        n['runtime_blas_threads'][0]['num_threads'] = 2
    elif damage == 'synthetic':
        n['synthetic_test_fixture'] = True
    else:
        n['elapsed_seconds'] = float('nan')
    with pytest.raises(ValueError):
        c.validate_native({'numerics': n}, path, spec)


def test_reference_result_cannot_fill_a_primary_key():
    path = c.r.default_model().inflow(1, 11).path(1)
    plan = c.pilot.step_plan(c.MAX_STEPS)
    with pytest.raises(ValueError):
        c.validate_native({'numerics': native_numbers(path, plan.reference[0])}, path, plan.primary[0])


def test_source_bundle_survives_relocation_and_live_source_changes(source_root, tmp_path):
    directory = tmp_path/'campaign'
    c.create_campaign(directory, root=source_root)
    (source_root/'gabes/example.py').write_text('coefficient = 999\n')
    relocated = tmp_path/'relocated'
    c.extract_bundle(directory, relocated)
    assert (relocated/'gabes/example.py').read_text() == 'coefficient = 1\n'
    assert c.load_plan(directory, root=relocated)['schema'] == c.SCHEMA
    with pytest.raises(FileExistsError):
        c.create_campaign(directory, root=source_root)
    with pytest.raises(FileExistsError):
        c.extract_bundle(directory, relocated)


def test_corrupt_archive_refused_before_extract(source_root, tmp_path):
    directory = tmp_path/'campaign'
    c.create_campaign(directory, root=source_root)
    with (directory/'sources.zip').open('ab') as handle:
        handle.write(b'changed')
    with pytest.raises(ValueError, match='bundle changed'):
        c.extract_bundle(directory, tmp_path/'relocated')


def test_duplicate_json_keys_not_accepted(tmp_path):
    file = tmp_path/'bad.json'
    file.write_text('{"record_sha256":"a","record_sha256":"b"}')
    with pytest.raises(ValueError, match='duplicate'):
        c.records.read_sealed(file)


def test_live_interpreter_cannot_publish_as_captured_execution():
    with pytest.raises(ValueError, match='fresh captured-source interpreter'):
        c.require_capsule({'record_sha256': '0'*64})
