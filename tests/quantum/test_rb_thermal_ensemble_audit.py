"""Pilot-workload, native reference identity and immutable output contracts.

All stored numerical values below are explicit synthetic test data in temporary
directories. Native physical metadata and production record validation are used;
no independent or primary atomic ODE is executed by these tests.
"""

from copy import deepcopy
from dataclasses import replace
import hashlib
import json
from pathlib import Path

import numpy as np
import pytest

from analysis.grand_challenge import rb_thermal_ensemble_audit as audit
from analysis.grand_challenge import rb_thermal_reference_jobs as refs
from gabes.fwm_quantum.kinetic import CarrierGeometry
from gabes.quantum.contracts import AnalysisFrequencyAxis


@pytest.fixture(autouse=True)
def forbid_atomic_solves(monkeypatch):
    def forbidden(*args, **kwargs):
        pytest.fail('audit contract tests must not start atomic ODEs or primary workers')
    monkeypatch.setattr(refs, 'adjoint_wavepacket', forbidden)
    monkeypatch.setattr(audit, '_primary_worker', forbidden)
    monkeypatch.setattr(audit, 'ProcessPoolExecutor', forbidden)


@pytest.fixture(scope='module')
def model():
    return audit.r.default_model()


@pytest.fixture(scope='module')
def reference_plan():
    return refs.prepare_plan()


def synthetic_record(identity, snapshot):
    """Valid storage fixture; synthetic work counters are never real ODE evidence."""
    nf, ports = np.asarray(identity['port_frequencies_rad_s']).shape
    ns = len(identity['source_names'])
    tau = identity['path']['residence_time_s']
    scale = tau**2*(1+1e-8*identity['reference_level'])
    by_source = np.broadcast_to(np.eye(ports, dtype=complex), (ns, nf, ports, ports)).copy()
    by_source *= scale*np.arange(1, ns+1)[:, None, None, None]/(ns*100)
    mean = tau*np.broadcast_to(np.array([.04+.01j, .02-.03j, .04-.01j, .02+.03j]), (nf, ports)).copy()
    values = {'greater_by_source': by_source, 'lesser_by_source': .7*by_source,
        'greater': by_source.sum(axis=0), 'lesser': (.7*by_source).sum(axis=0),
        'mean_pulse': mean, 'mean_outer': mean[:, :, None]*mean[:, None, :].conj(),
        'retarded_response': np.broadcast_to((.2+.3j)*scale*np.eye(ports), (nf, ports, ports)).copy(),
        'exit_state': audit.r.decode(identity['atomic_problem']['boundary_state'])}
    numerics = {**identity['numerical_settings'], 'forward_density_evaluations': 1,
        'backward_evaluations': 1, 'density_mesh_points': 2,
        'backward_complex_variables': 1, 'elapsed_seconds': 0.,
        'synthetic_test_fixture': True}
    return refs.sealed({'schema': 'gabes-rb-thermal-reference-job-v1',
        'identity': identity, 'full_identity_sha256': refs.digest(identity),
        'values': refs.encode(values), 'source_names': identity['source_names'],
        'numerics': numerics, 'frequencies_rad_s': identity['port_frequencies_rad_s'],
        'residence_time_s': tau, 'environment': identity['environment'],
        'source_before': snapshot, 'source_after': snapshot, 'source_stable_during_run': True,
        'scope': refs.SCOPE, 'runtime_blas_threads': [
            {'user_api': 'blas', 'num_threads': 1, 'internal_api': 'synthetic storage test fixture'}]})


@pytest.fixture
def reference_directory(tmp_path, reference_plan):
    snapshot, identities = reference_plan
    directory = tmp_path/'references'
    directory.mkdir()
    for identity in identities:
        record = synthetic_record(identity, snapshot)
        # Exercise the real complete native identity and storage validator.
        refs.validate_record(record, identity, snapshot)
        filename = directory/refs.job_filename(identity['path_index'], identity['reference_level'])
        with filename.open('x', encoding='utf-8') as output:
            output.write(json.dumps(record, ensure_ascii=False, allow_nan=False))
    return directory


def reseal(record):
    return refs.sealed({key: value for key, value in record.items() if key != 'record_sha256'})


def test_workload_counts_actual_nested_boundary_paths_without_solving(model):
    table = audit.workload(model)
    assert table['grid_path_entries'] == 270
    assert table['unique_physical_paths'] == 144
    assert table['calculations_per_path'] == 5
    assert len(table['grids']) == 12
    assert all(row['status'] == 'NOT_EVALUATED' for row in table['grids'])
    by_power = {power: [row for row in table['grids'] if row['power'] == power] for power in range(4)}
    for power, rows in by_power.items():
        assert {row['seed'] for row in rows} == {11, 211, 811}
        assert all(row['path_count'] == 6*2**power for row in rows)
        expected_new = 6 if power == 0 else 6*2**(power-1)
        assert all(row['new_unique_physical_paths'] == expected_new for row in rows)
        assert all(row['reused_physical_paths']+row['new_unique_physical_paths'] == row['path_count'] for row in rows)
        assert all(0 <= row['tau_over_5us_arrival_fraction'] <= 1 for row in rows)
        assert all(0 <= row['tau_over_5us_occupancy_fraction'] <= 1 for row in rows)


def test_reference_jobs_and_pilot_share_exact_physical_paths_and_canonical_phase_metadata(model, reference_plan):
    _, jobs = reference_plan
    inflow = model.inflow(0, 11)
    assert len(jobs) == 12
    for identity in jobs:
        path = inflow.path(identity['path_index'])
        problem = model.problem(path)
        decoded = audit.r.decode(identity)
        for key in ('entry_position_m', 'velocity_m_s'):
            assert np.asarray(decoded['path'][key]).tobytes() == getattr(path, key).tobytes()
        assert decoded['path']['residence_time_s'] == path.residence_time_s
        for key in ('h0', 'h1', 'boundary_state', 'readouts', 'drives'):
            np.testing.assert_array_equal(decoded['atomic_problem'][key], problem[key])
        assert audit.r.encode(decoded['physical_metadata']) == audit.r.encode(problem['metadata'])
        assert isinstance(identity['physical_metadata']['mode_labels'], list)
        assert isinstance(problem['metadata']['mode_labels'], tuple)
        assert identity['path']['residence_time_s'] != 2e-6


def test_verified_references_use_actual_sealed_record_validation_and_file_hashes(model, reference_directory):
    values = audit.verified_references(reference_directory, model)
    assert [(row['path_index'], row['level']) for row in values] == [(path, level) for path in range(6) for level in range(2)]
    for row in values:
        filename = reference_directory/refs.job_filename(row['path_index'], row['level'])
        assert row['file_sha256'] == hashlib.sha256(filename.read_bytes()).hexdigest()
        assert len(row['record_sha256']) == len(row['identity_sha256']) == 64
        packet = row['packet']
        assert set(refs.METRICS) <= set(packet)
        mean = packet['mean_pulse']
        np.testing.assert_allclose(packet['mean_outer'], mean[:, :, None]*mean[:, None, :].conj(), rtol=1e-15)
        assert np.linalg.norm(packet['retarded_response'].imag) > 0.


@pytest.mark.parametrize('damage', ['missing', 'unsealed_payload', 'identity', 'source_hash', 'frequency',
    'missing_metric', 'mean_outer', 'source_closure', 'blas_threads', 'work_counter', 'duplicate_json_key'])
def test_verified_references_reject_missing_tampered_or_incomplete_jobs(model, reference_directory, damage):
    filename = reference_directory/refs.job_filename(2, 1)
    if damage == 'missing':
        filename.unlink()
    elif damage == 'duplicate_json_key':
        filename.write_text('{"record_sha256":"first","record_sha256":"second"}', encoding='utf-8')
    else:
        record = json.loads(filename.read_text(encoding='utf-8'))
        if damage == 'unsealed_payload':
            record['values']['mean_pulse']['real'][0][0] += 1e-9
        elif damage == 'identity':
            record['identity']['path']['residence_time_s'] *= 1.01
            record['full_identity_sha256'] = refs.digest(record['identity'])
        elif damage == 'source_hash':
            name = next(iter(record['source_before']['source_sha256']))
            record['source_before']['source_sha256'][name] = '0'*64
        elif damage == 'frequency':
            record['frequencies_rad_s'][0][0] += 1e6
        elif damage == 'missing_metric':
            del record['values']['retarded_response']
        elif damage == 'mean_outer':
            record['values']['mean_outer']['real'][0][0][0] *= 2
        elif damage == 'source_closure':
            record['values']['greater']['real'][0][0][0] *= 2
        elif damage == 'blas_threads':
            record['runtime_blas_threads'][0]['num_threads'] = 2
        elif damage == 'work_counter':
            record['numerics']['backward_evaluations'] = 0
        if damage != 'unsealed_payload':
            record = reseal(record)
        filename.write_text(json.dumps(record), encoding='utf-8')
    with pytest.raises((ValueError, FileNotFoundError)):
        audit.verified_references(reference_directory, model)


@pytest.mark.parametrize('change', ['inputs', 'box', 'temperature', 'geometry', 'rf', 'state', 'phase'])
def test_references_cannot_be_attached_to_a_different_physical_model(model, reference_directory, change):
    if change == 'inputs':
        changed = replace(model, inputs=replace(model.inputs, pump_power_W=model.inputs.pump_power_W*1.001))
    elif change == 'box':
        changed = replace(model, upper_corner_m=model.upper_corner_m+[0., 0., 1e-5])
    elif change == 'temperature':
        changed = replace(model, temperature_K=model.temperature_K+1)
    elif change == 'geometry':
        changed = replace(model, geometry=CarrierGeometry.vacuum_beams(model.inputs, probe_angle_rad=.0061, conjugate_angle_rad=-.005))
    elif change == 'rf':
        changed = replace(model, analysis_axis=AnalysisFrequencyAxis(model.analysis_axis.omega_rad_s*1.001))
    elif change == 'state':
        changed = replace(model, boundary_state=np.diag([.4, .6, 0., 0.]))
    else:
        changed = replace(model, entry_phase_rad=.01)
    with pytest.raises(ValueError):
        audit.verified_references(reference_directory, changed)


@pytest.mark.parametrize('segments', [(256.5, 512, 1024), (True, 512, 1024), ('256', 512, 1024),
    (0, 512, 1024), (-1, 512, 1024), (512, 256, 1024), (256, 256, 1024),
    (256, 512), (128, 256, 512, 1024)])
def test_primary_plan_requires_three_positive_increasing_integers(segments):
    with pytest.raises((ValueError, TypeError)):
        audit.plan_for(segments)


def test_primary_plan_preserves_exact_counts_and_independent_reference_tolerances():
    plan = audit.plan_for((128, 256, 512))
    assert [item.parameters['segments'] for item in plan.primary] == [128, 256, 512]
    assert [item.resolution for item in plan.primary] == [128, 256, 512]
    assert [(item.parameters['rtol'], item.parameters['atol']) for item in plan.reference] == list(refs.REFERENCE_TOLERANCES)
    assert all(item.parameters['max_step_fraction'] == 1/64 for item in plan.reference)


def test_maximum_time_step_plan_uses_each_full_physical_residence_time(model):
    maximum = (1e-6, 4e-7, 1e-7)
    plan = audit.step_plan(maximum)
    inflow = model.inflow(0, 11)
    counts = []
    for index in range(6):
        path = inflow.path(index)
        current = []
        for spec, limit in zip(plan.primary, maximum):
            parameters = audit.solver_parameters(path, spec)
            count = parameters['segments']
            assert count == int(np.ceil(path.residence_time_s/limit))
            assert path.residence_time_s/count <= limit*(1+1e-15)
            assert (count-1)*limit < path.residence_time_s
            assert parameters['order'] == 4 and parameters['sample_count'] == 17
            assert 'max_step_s' not in parameters
            assert spec.parameters['max_step_s'] == limit
            current.append(count)
        counts.append(tuple(current))
    assert len(set(counts)) > 1


@pytest.mark.parametrize('steps', [(1e-9, 2e-9, 3e-9), (3e-9, 3e-9, 1e-9),
    (3e-9, 2e-9, 0.), (3e-9, 2e-9, np.nan), (np.inf, 2e-9, 1e-9),
    (True, 2e-9, 1e-9), ('3e-9', 2e-9, 1e-9), (2e-9, 1e-9)])
def test_time_step_plan_requires_three_finite_positive_decreasing_limits(steps):
    with pytest.raises((ValueError, TypeError)):
        audit.step_plan(steps)


def test_reference_only_pilot_cannot_certify_paths_grid_or_thermal_ensemble(reference_directory, tmp_path):
    report = audit.build_report(reference_directory, tmp_path/'cache', compute_primary=False)
    assert len(report['references']) == len(report['reference_cache']) == 12
    assert not report['primary_cache_jobs']
    assert not report['pilot_all_path_evidence_passed']
    assert report['reference_refinements_passed']
    assert all(row['status'] == 'INCOMPLETE' for row in report['pilot_paths'])
    assert report['pilot_grid']['candidate']['spectra'] is None
    assert not report['ensemble_gate']['passed']
    assert not report['thermal_ensemble_converged']
    assert not report['production_adapter_certified']
    assert not report['physical_optical_prediction']
    assert report['workload']['grids'][0]['status'] == 'INCOMPLETE'
    assert all(row['status'] == 'NOT_EVALUATED' for row in report['workload']['grids'][1:])


def test_even_a_complete_synthetic_pilot_cache_cannot_certify_the_thermal_ensemble(
        model, reference_directory, tmp_path):
    # Match the runner's source inventory before seeding test-only cached data.
    from analysis.grand_challenge.reference import exponential_transport
    assert exponential_transport is not None
    extra = [Path(audit.__file__), audit.r.ROOT/'tests/quantum/test_rb_thermal_ensemble.py', Path(__file__)]
    cache_directory = tmp_path/'complete-synthetic-pilot'
    cache = audit.r.PathCache(cache_directory, audit.r.consumed_hashes(extra))
    references = audit.verified_references(reference_directory, model)
    plan = audit.plan_for()
    inflow = model.inflow(0, 11)
    for index in range(6):
        packet = references[2*index+1]['packet']
        for spec in plan.primary:
            cache.get(model, inflow.path(index), spec,
                provider=lambda model, path, spec, packet=packet: deepcopy(packet))
    report = audit.build_report(reference_directory, cache_directory, compute_primary=False)
    assert report['pilot_all_path_evidence_passed']
    assert report['reference_refinements_passed']
    assert report['pilot_grid']['candidate']['spectra'] is not None
    assert report['workload']['grids'][0]['status'] == 'UNCONVERGED'
    assert not report['ensemble_gate']['passed']
    assert not report['thermal_ensemble_converged']
    assert not report['production_adapter_certified']
    assert not report['physical_optical_prediction']
    assert not report['primary_cache_jobs']


@pytest.mark.parametrize('changed', ['matrices', 'numerics'])
def test_other_valid_cached_reference_cannot_replace_the_verified_immutable_job(
        model, reference_directory, tmp_path, changed):
    from analysis.grand_challenge.reference import exponential_transport
    assert exponential_transport is not None
    extra = [Path(audit.__file__), audit.r.ROOT/'tests/quantum/test_rb_thermal_ensemble.py', Path(__file__)]
    directory = tmp_path/'mismatched-reference-cache'
    cache = audit.r.PathCache(directory, audit.r.consumed_hashes(extra))
    reference = audit.verified_references(reference_directory, model)[0]
    altered = deepcopy(reference['packet'])
    if changed == 'matrices':
        altered['greater'] *= 1.01
        altered['greater_by_source'] *= 1.01
    else:
        altered['numerics']['backward_evaluations'] += 1
    # This alternate packet passes the actual cache's PSD/identity checks.
    _, stored = cache.get(model, model.inflow(0, 11).path(0), audit.plan_for().reference[0],
        provider=lambda model, path, spec: altered)
    before = Path(stored['path']).read_bytes()
    with pytest.raises(ValueError, match='Cached reference differs'):
        audit.build_report(reference_directory, directory, compute_primary=False)
    assert Path(stored['path']).read_bytes() == before


@pytest.mark.parametrize('existing', ['report', 'plot'])
def test_cli_rejects_existing_artifacts_before_reference_reads_or_build(tmp_path, monkeypatch, existing):
    output, plot = tmp_path/'report.json', tmp_path/'plot.png'
    protected = output if existing == 'report' else plot
    protected.write_bytes(b'original bytes')
    monkeypatch.setattr(audit, 'build_report', lambda *args, **kwargs: pytest.fail('preflight must reject before build'))
    with pytest.raises(FileExistsError):
        audit.main(['--output', str(output), '--plot', str(plot), '--cache-dir', str(tmp_path/'cache')])
    assert protected.read_bytes() == b'original bytes'
    assert not (plot if existing == 'report' else output).exists()
    assert not (tmp_path/'cache').exists()


def test_cli_rejects_same_new_report_and_plot_path_before_build(tmp_path, monkeypatch):
    path = tmp_path/'unwritten'
    monkeypatch.setattr(audit, 'build_report', lambda *args, **kwargs: pytest.fail('preflight must reject before build'))
    with pytest.raises(ValueError):
        audit.main(['--output', str(path), '--plot', str(path), '--cache-dir', str(tmp_path/'cache')])
    assert not path.exists()


def test_cli_cannot_mix_fixed_segment_counts_and_maximum_time_steps(tmp_path, monkeypatch):
    monkeypatch.setattr(audit, 'build_report', lambda *args, **kwargs: pytest.fail('argument conflict must reject before build'))
    with pytest.raises(SystemExit) as error:
        audit.main(['--output', str(tmp_path/'unused.json'), '--plot', str(tmp_path/'unused.png'),
            '--cache-dir', str(tmp_path/'unused-cache'), '--segments', '128', '256', '512',
            '--max-steps-s', '3e-9', '2e-9', '1e-9'])
    assert error.value.code == 2
    assert not list(tmp_path.iterdir())


def test_source_change_during_build_prevents_report_and_plot_creation(reference_directory, tmp_path, monkeypatch):
    hashes = audit.r.consumed_hashes
    calls = 0

    def changed_at_final_snapshot(extra=()):
        nonlocal calls
        calls += 1
        result = hashes(extra)
        if calls >= 2:
            result = {**result, 'simulated_changed_dependency.py': '0'*64}
        return result

    monkeypatch.setattr(audit.r, 'consumed_hashes', changed_at_final_snapshot)
    output, plot = tmp_path/'never-written.json', tmp_path/'never-written.png'
    with pytest.raises(RuntimeError, match='[Ss]ource|changed'):
        audit.main(['--output', str(output), '--plot', str(plot), '--reference-dir', str(reference_directory),
                    '--cache-dir', str(tmp_path/'cache')])
    assert calls >= 2
    assert not output.exists() and not plot.exists()
