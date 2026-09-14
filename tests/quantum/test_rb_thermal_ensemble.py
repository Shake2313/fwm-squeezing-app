"""Cache and ensemble contracts using native Rb geometry and synthetic packets.

The injected provider exercises the executor without large optical-frequency
ODEs. Its path-dependent Gram matrices are contract fixtures, not Rb predictions.
"""

from copy import deepcopy
from dataclasses import replace
import hashlib
import itertools

import numpy as np
import pytest

from analysis.grand_challenge import rb_thermal_ensemble as thermal
from gabes.fwm_quantum.kinetic import CarrierGeometry
from gabes.fwm_quantum.transport_ensemble import packet_digest
from gabes.quantum.contracts import AnalysisFrequencyAxis


PATH_KEYS = ('greater', 'lesser', 'greater_by_source', 'lesser_by_source',
             'mean_pulse', 'mean_outer', 'retarded_response')
STREAM_KEYS = ('greater', 'lesser', 'greater_by_source', 'lesser_by_source',
               'poisson_number', 'retarded_response')
CHARGES = np.array([1, -1, -1, 1])


@pytest.fixture
def source_manifest(tmp_path):
    source = tmp_path/'provider_source.py'
    source.write_text('# Explicit synthetic packet provider source fixture\n', encoding='utf-8')
    return source, {str(source): hashlib.sha256(source.read_bytes()).hexdigest()}


@pytest.fixture
def cache(tmp_path, source_manifest):
    return thermal.PathCache(tmp_path/'cache', source_manifest[1])


@pytest.fixture(scope='module')
def model():
    return thermal.default_model()


@pytest.fixture(scope='module')
def plan():
    return thermal.PathPlan(
        primary=tuple(thermal.SolverSpec('synthetic-primary', float(level), {'algorithm': 'test Gram packets'})
                      for level in (8, 16, 32)),
        reference=tuple(thermal.SolverSpec('synthetic-independent', float(level), {'algorithm': 'independent test family'})
                        for level in (16, 32)))


def synthetic_packet(model, path, spec):
    """Each source, RF and physical path has a distinct finite full matrix."""
    problem = model.problem(path)
    names = ('atomic_inflow',)+tuple('jump:'+item.name for item in problem['reservoirs'].channels)
    nf, ns = len(model.analysis_axis.omega_rad_s), len(names)
    rng = np.random.default_rng(723)
    a = rng.normal(size=(ns, nf, 4, 4))+1j*rng.normal(size=(ns, nf, 4, 4))
    b = rng.normal(size=(ns, nf, 4, 4))+1j*rng.normal(size=(ns, nf, 4, 4))
    path_factor = 1+.08*np.sin(path.entry_position_m@[700., 900., 300.])
    path_factor += .06*np.cos(path.velocity_m_s@[.004, -.006, .002])
    signed_coefficient = 2e-4 if spec.method == 'synthetic-primary' else -1e-4
    correction = signed_coefficient/spec.resolution**2
    scale = path.residence_time_s**2*path_factor*(1+correction)/ns/100
    greater = scale*(a@a.conj().swapaxes(-1, -2))
    lesser = .6*scale*(b@b.conj().swapaxes(-1, -2))
    amplitudes = np.array([.13+.04j, .09-.07j, .13-.04j, .09+.07j])
    mean = path.residence_time_s*(1+correction)*path_factor*np.arange(1, nf+1)[:, None]*amplitudes/2
    response = path.residence_time_s**2*(1+correction)*path_factor*(
        rng.normal(size=(nf, 4, 4))+1j*rng.normal(size=(nf, 4, 4)))
    return {'analysis_axis': model.analysis_axis, 'frequencies_rad_s': problem['frequencies_rad_s'],
        'source_names': names, 'greater': greater.sum(axis=0), 'lesser': lesser.sum(axis=0),
        'greater_by_source': greater, 'lesser_by_source': lesser,
        'mean_pulse': mean, 'retarded_response': response,
        'residence_time_s': path.residence_time_s, 'exit_state': model.boundary_state,
        'audit': {'passed': True}, 'metadata': problem['metadata']}


class CountingProvider:
    def __init__(self):
        self.calls = []

    def __call__(self, model, path, spec):
        self.calls.append((path.entry_position_m.copy(), path.velocity_m_s.copy(),
                           path.residence_time_s, spec.method, spec.resolution))
        return synthetic_packet(model, path, spec)


def literal_phase_average(packet):
    names = ('greater', 'lesser', 'greater_by_source', 'lesser_by_source', 'retarded_response')
    result = {key: np.zeros_like(packet[key]) for key in names}
    result['poisson_number'] = np.zeros_like(packet['greater'])
    result['mean_pulse'] = np.zeros_like(packet['mean_pulse'])
    for phase in .217+np.arange(4)*np.pi/2:
        d = np.diag(np.exp(1j*CHARGES*phase))
        for key in names:
            result[key] += d@packet[key]@d.conj().T/4
        mean = packet['mean_pulse']@d.T
        result['poisson_number'] += np.einsum('fi,fj->fij', mean, mean.conj())/4
        result['mean_pulse'] += mean/4
    for key in ('greater', 'lesser'):
        result[key] += result['poisson_number']
    return result


def test_default_model_uses_declared_area_length_and_explicit_boundary_transport(model):
    lengths = model.upper_corner_m-model.lower_corner_m
    np.testing.assert_allclose(lengths[:2], np.sqrt(model.inputs.uniform_area_m2), rtol=1e-15)
    assert lengths[2] == pytest.approx(model.inputs.length_m, rel=1e-15)
    assert model.inputs.transit_rate_s_inverse == 0.
    assert model.inputs.phase_mismatch_rad_m == 0.
    contract = model.convention()
    assert contract.number_density_m3 == model.inputs.number_density_m3
    assert len(contract.source_names) > 1
    assert np.linalg.norm(model.geometry.mismatch_rad_m) > 0


def test_cache_reuses_one_computation_and_returns_complete_complex_packet(model, plan, cache):
    path = model.inflow(0, 11).path(0)
    provider = CountingProvider()
    first, _ = cache.get(model, path, plan.primary[-1], provider=provider)
    loaded, _ = cache.get(model, path, plan.primary[-1])
    assert len(provider.calls) == 1
    for key in ('greater', 'lesser', 'greater_by_source', 'lesser_by_source',
                'mean_pulse', 'retarded_response', 'exit_state'):
        np.testing.assert_array_equal(loaded[key], first[key])
    assert tuple(loaded['source_names']) == tuple(model.convention().source_names)
    assert loaded['audit']['passed']
    assert packet_digest(path, loaded, model.convention()) == packet_digest(path, first, model.convention())


def test_existing_cache_never_invokes_a_replacement_provider_or_overwrites_values(model, plan, cache):
    path = model.inflow(0, 11).path(0)
    first, _ = cache.get(model, path, plan.primary[-1], provider=synthetic_packet)

    def forbidden(*args):
        pytest.fail('an existing immutable cache entry must be reused')

    loaded, _ = cache.get(model, path, plan.primary[-1], provider=forbidden)
    np.testing.assert_array_equal(loaded['mean_pulse'], first['mean_pulse'])
    original = loaded['mean_pulse'].copy()
    try:
        loaded['mean_pulse'][0, 0] *= 2
    except ValueError:
        pass
    again, _ = cache.get(model, path, plan.primary[-1])
    np.testing.assert_array_equal(again['mean_pulse'], original)


def test_cache_identity_ignores_grid_node_label_but_evidence_rebinds_current_path(model, plan, cache):
    path = model.inflow(0, 11).path(3)
    changed_source = replace(path, source='Same physical chord with a different grid/node provenance label')
    provider = CountingProvider()
    first, _ = thermal.build_path_evidence(model, path, plan, cache, provider=provider)
    second, _ = thermal.build_path_evidence(model, changed_source, plan, cache)
    assert len(provider.calls) == 5
    assert first.evidence.target_digest == packet_digest(path, first.packet, first.convention)
    assert second.evidence.target_digest == packet_digest(changed_source, second.packet, second.convention)
    assert second.evidence.target_digest != first.evidence.target_digest
    for key in STREAM_KEYS:
        if key != 'poisson_number':
            np.testing.assert_array_equal(second.packet[key], first.packet[key])


def test_missing_cache_requires_explicit_provider_and_never_becomes_successful_evidence(model, plan, cache):
    path = model.inflow(0, 11).path(0)
    with pytest.raises(FileNotFoundError):
        cache.get(model, path, plan.primary[-1])
    candidate, row = thermal.run_grid(model, 0, 11, plan, cache)
    assert candidate['spectra'] is None
    assert not candidate['certified']
    assert candidate['reasons']


@pytest.mark.parametrize('field', ['entry', 'velocity', 'duration', 'power', 'waist', 'density',
    'geometry', 'rf', 'boundary_state', 'pump_center', 'entry_phase', 'solver_resolution', 'solver_parameters'])
def test_each_consumed_physical_or_solver_change_causes_cache_miss(model, plan, cache, field):
    path = model.inflow(0, 11).path(0)
    spec = plan.primary[-1]
    cache.get(model, path, spec, provider=synthetic_packet)
    changed_model, changed_path, changed_spec = model, path, spec
    if field == 'entry':
        changed_path = replace(path, entry_position_m=path.entry_position_m+[0., 1e-8, 0.])
    elif field == 'velocity':
        changed_path = replace(path, velocity_m_s=path.velocity_m_s+[.01, 0., 0.])
    elif field == 'duration':
        changed_path = replace(path, residence_time_s=path.residence_time_s*1.00001)
    elif field in ('power', 'waist', 'density'):
        name = {'power': 'pump_power_W', 'waist': 'pump_waist_m', 'density': 'number_density_m3'}[field]
        inputs = replace(model.inputs, **{name: getattr(model.inputs, name)*1.001})
        changed_model = replace(model, inputs=inputs)
    elif field == 'geometry':
        changed_model = replace(model, geometry=CarrierGeometry.vacuum_beams(
            model.inputs, probe_angle_rad=.0061, conjugate_angle_rad=-.005))
    elif field == 'rf':
        changed_model = replace(model, analysis_axis=AnalysisFrequencyAxis(model.analysis_axis.omega_rad_s*1.001))
    elif field == 'boundary_state':
        changed_model = replace(model, boundary_state=np.diag([.4, .6, 0., 0.]))
    elif field == 'pump_center':
        changed_model = replace(model, pump_center_xy_m=[1e-8, 0.])
    elif field == 'entry_phase':
        changed_model = replace(model, entry_phase_rad=.001)
    elif field == 'solver_resolution':
        changed_spec = replace(spec, resolution=spec.resolution*2)
    elif field == 'solver_parameters':
        changed_spec = replace(spec, parameters={**spec.parameters, 'new_control': 1.})
    with pytest.raises(FileNotFoundError):
        cache.get(changed_model, changed_path, changed_spec)


def test_source_change_invalidates_existing_cache_instance_and_new_manifest_has_new_identity(
        model, plan, cache, source_manifest, tmp_path):
    path = model.inflow(0, 11).path(0)
    cache.get(model, path, plan.primary[-1], provider=synthetic_packet)
    source = source_manifest[0]
    source.write_text('# Changed numerical implementation\n', encoding='utf-8')
    with pytest.raises((ValueError, RuntimeError)):
        cache.get(model, path, plan.primary[-1])
    updated = {str(source): hashlib.sha256(source.read_bytes()).hexdigest()}
    other = thermal.PathCache(tmp_path/'cache', updated)
    with pytest.raises(FileNotFoundError):
        other.get(model, path, plan.primary[-1])


def test_wrong_source_digest_is_rejected_before_any_cache_use(source_manifest, tmp_path):
    source, _ = source_manifest
    with pytest.raises((ValueError, RuntimeError)):
        thermal.PathCache(tmp_path/'wrong-cache', {str(source): '0'*64})


def test_missing_source_file_is_not_accepted_as_provenance(tmp_path):
    with pytest.raises((ValueError, RuntimeError, FileNotFoundError)):
        thermal.PathCache(tmp_path/'cache', {str(tmp_path/'missing.py'): 'a'*64})


def test_source_modified_during_provider_call_cannot_create_a_stale_cache(
        model, plan, cache, source_manifest, tmp_path):
    path = model.inflow(0, 11).path(0)

    def changes_its_source(model, path, spec):
        packet = synthetic_packet(model, path, spec)
        source_manifest[0].write_text('# Changed while the calculation was running\n', encoding='utf-8')
        return packet

    with pytest.raises((ValueError, RuntimeError)):
        cache.get(model, path, plan.primary[-1], provider=changes_its_source)
    assert not [path for path in (tmp_path/'cache').rglob('*') if path.is_file()]


def test_cache_hit_rechecks_sources_after_packet_validation(model, plan, cache, source_manifest, monkeypatch):
    path = model.inflow(0, 11).path(0)
    cache.get(model, path, plan.primary[-1], provider=synthetic_packet)
    validate = thermal._validated_packet

    def source_changes_during_read(model, path, raw):
        packet = validate(model, path, raw)
        source_manifest[0].write_text('# Changed during cache validation\n', encoding='utf-8')
        return packet

    monkeypatch.setattr(thermal, '_validated_packet', source_changes_during_read)
    with pytest.raises((ValueError, RuntimeError)):
        cache.get(model, path, plan.primary[-1])


@pytest.mark.parametrize('failure', ['wrong_rf', 'missing_source', 'source_total_mismatch',
    'negative_covariance', 'nonhermitian_covariance', 'bad_exit_trace', 'bad_exit_psd', 'nan_response'])
def test_audit_true_cannot_override_invalid_packet_data(model, plan, cache, failure):
    path = model.inflow(0, 11).path(0)

    def invalid(model, path, spec):
        packet = synthetic_packet(model, path, spec)
        if failure == 'wrong_rf':
            packet['frequencies_rad_s'] = packet['frequencies_rad_s']+1e5
        elif failure == 'missing_source':
            packet['source_names'] = packet['source_names'][:-1]
        elif failure == 'source_total_mismatch':
            packet['greater_by_source'][0] *= 2
        elif failure == 'negative_covariance':
            packet['greater_by_source'][0] *= -100
            packet['greater'] = packet['greater_by_source'].sum(axis=0)
        elif failure == 'nonhermitian_covariance':
            packet['lesser'][0, 0, 1] += path.residence_time_s**2
        elif failure == 'bad_exit_trace':
            packet['exit_state'] = 2*packet['exit_state']
        elif failure == 'bad_exit_psd':
            packet['exit_state'] = np.diag([1.1, -.1, 0., 0.])
        elif failure == 'nan_response':
            packet['retarded_response'][0, 0, 0] = np.nan
        return packet

    with pytest.raises((ValueError, TypeError)):
        cache.get(model, path, plan.primary[-1], provider=invalid)
    with pytest.raises(FileNotFoundError):
        cache.get(model, path, plan.primary[-1])


def test_absent_audit_is_checked_but_a_lone_cached_packet_is_not_path_evidence(model, plan, cache):
    path = model.inflow(0, 11).path(0)

    def without_audit(model, path, spec):
        packet = synthetic_packet(model, path, spec)
        del packet['audit']
        return packet

    packet, _ = cache.get(model, path, plan.primary[-1], provider=without_audit)
    assert packet['audit']['passed']
    candidate, _ = thermal.run_grid(model, 0, 11, plan, cache)
    assert not candidate['certified'] and candidate['spectra'] is None


def test_explicit_failed_atomic_audit_cannot_be_promoted_by_cache_reload(model, plan, cache):
    path = model.inflow(0, 11).path(0)

    def failed(model, path, spec):
        packet = synthetic_packet(model, path, spec)
        packet['audit']['passed'] = False
        return packet

    packet, _ = cache.get(model, path, plan.primary[-1], provider=failed)
    loaded, _ = cache.get(model, path, plan.primary[-1])
    assert not packet['audit']['passed'] and not loaded['audit']['passed']


@pytest.mark.parametrize('damage', ['partial_json', 'changed_payload', 'changed_identity'])
def test_corrupt_cache_records_are_rejected_and_never_recomputed_or_overwritten(model, plan, cache, damage):
    import json
    from pathlib import Path
    path = model.inflow(0, 11).path(0)
    _, record = cache.get(model, path, plan.primary[-1], provider=synthetic_packet)
    filename = Path(record['path'])
    if damage == 'partial_json':
        filename.write_text('{"key":', encoding='utf-8')
    else:
        saved = json.loads(filename.read_text(encoding='utf-8'))
        if damage == 'changed_payload':
            saved['packet']['mean_pulse']['real'][0][0] += 1e-8
        else:
            saved['identity']['physical_path']['residence_time_s'] *= 2
        filename.write_text(json.dumps(saved), encoding='utf-8')
    corrupt = filename.read_bytes()

    def forbidden(*args):
        pytest.fail('an invalid immutable record must not be replaced by a new solve')

    with pytest.raises((ValueError, TypeError)):
        cache.get(model, path, plan.primary[-1], provider=forbidden)
    assert filename.read_bytes() == corrupt


def test_path_evidence_contains_actual_nonzero_comparisons_for_all_ordered_outputs(model, plan, cache):
    path = model.inflow(0, 11).path(0)
    supplied, ledger = thermal.build_path_evidence(model, path, plan, cache, provider=synthetic_packet)
    evidence = supplied.evidence
    assert evidence.target_digest == packet_digest(path, supplied.packet, supplied.convention)
    primary_rows = [row for row in evidence.comparisons if row.kind == 'path_refinement']
    independent_rows = [row for row in evidence.comparisons if row.kind == 'independent_reference']
    assert len(primary_rows) == 2 and independent_rows
    assert primary_rows[0].candidate_id == primary_rows[1].reference_id
    assert primary_rows[-1].candidate_id == evidence.target_digest
    assert independent_rows[-1].candidate_id == evidence.target_digest
    for comparison in primary_rows+independent_rows:
        assert set(PATH_KEYS) <= set(comparison.errors)
        assert comparison.passed
        assert max(comparison.errors.values()) > 0.
        assert min(comparison.errors.values()) >= 0.


@pytest.fixture
def completed_grid(model, plan, cache):
    provider = CountingProvider()
    candidate, row = thermal.run_grid(model, 0, 11, plan, cache, provider=provider)
    return candidate, row, provider


def test_native_geometry_grid_matches_literal_phase_sum_of_every_source_rf_and_mean_outer(
        model, plan, completed_grid):
    candidate, _, provider = completed_grid
    assert candidate['spectra'] is not None and candidate['path_evidence_passed']
    assert not candidate['certified']
    inflow = model.inflow(0, 11)
    assert len(provider.calls) == 5*len(inflow.rate_s_inverse)
    expected = None
    for index, rate in enumerate(inflow.rate_s_inverse):
        packet = synthetic_packet(model, inflow.path(index), plan.primary[-1])
        phase = literal_phase_average(packet)
        if expected is None:
            expected = {key: np.zeros_like(value) for key, value in phase.items()}
        for key, value in phase.items():
            expected[key] += rate*value
    for key in STREAM_KEYS:
        np.testing.assert_allclose(candidate['spectra'][key], expected[key],
            rtol=4e-13, atol=4e-15*np.linalg.norm(expected[key]))
    for key in ('greater', 'lesser'):
        np.testing.assert_allclose(candidate['spectra'][key],
            candidate['spectra'][key+'_by_source'].sum(axis=0)+candidate['spectra']['poisson_number'],
            rtol=1e-14)
    assert np.linalg.norm(candidate['spectra']['poisson_number']) > 0.
    assert np.linalg.norm(candidate['spectra']['retarded_response'].imag) > 0.


def test_nested_sobol_reuses_physical_paths_but_applies_each_new_grid_rate(model, plan, cache):
    provider = CountingProvider()
    previous_paths = 0
    for power in (0, 1, 2):
        candidate, _ = thermal.run_grid(model, power, 11, plan, cache, provider=provider)
        inflow = model.inflow(power, 11)
        count = len(inflow.rate_s_inverse)
        assert candidate['spectra'] is not None
        assert len(provider.calls) == 5*count
        if power:
            coarse = model.inflow(power-1, 11)
            for face in range(6):
                coarse_indices = np.flatnonzero(coarse.face_index == face)
                fine_indices = np.flatnonzero(inflow.face_index == face)[:len(coarse_indices)]
                np.testing.assert_array_equal(coarse.entry_position_m[coarse_indices], inflow.entry_position_m[fine_indices])
                np.testing.assert_array_equal(coarse.velocity_m_s[coarse_indices], inflow.velocity_m_s[fine_indices])
                np.testing.assert_allclose(inflow.rate_s_inverse[fine_indices],
                                           coarse.rate_s_inverse[coarse_indices]/2, rtol=1e-15)
        direct = None
        for index, rate in enumerate(inflow.rate_s_inverse):
            packet = synthetic_packet(model, inflow.path(index), plan.primary[-1])
            phase = literal_phase_average(packet)
            if direct is None:
                direct = {key: np.zeros_like(phase[key]) for key in STREAM_KEYS}
            for key in STREAM_KEYS:
                direct[key] += rate*phase[key]
        for key in STREAM_KEYS:
            np.testing.assert_allclose(candidate['spectra'][key], direct[key],
                rtol=5e-13, atol=5e-15*np.linalg.norm(direct[key]))
        assert count > previous_paths
        previous_paths = count


def test_changing_density_scales_the_complete_stream_once(model, plan, cache):
    first, _ = thermal.run_grid(model, 0, 11, plan, cache, provider=synthetic_packet)
    doubled = replace(model, inputs=replace(model.inputs, number_density_m3=2*model.inputs.number_density_m3))
    second, _ = thermal.run_grid(doubled, 0, 11, plan, cache, provider=synthetic_packet)
    assert first['spectra'] is not None and second['spectra'] is not None
    for key in first['spectra']:
        np.testing.assert_allclose(second['spectra'][key], 2*first['spectra'][key], rtol=3e-14)
    assert second['total_arrival_rate_s_inverse'] == pytest.approx(2*first['total_arrival_rate_s_inverse'], rel=1e-15)


def test_one_failed_path_hides_all_partial_spectra_and_never_renormalizes(model, plan, cache):
    inflow = model.inflow(0, 11)
    failed_position = inflow.path(2).entry_position_m

    def provider(model, path, spec):
        packet = synthetic_packet(model, path, spec)
        if np.array_equal(path.entry_position_m, failed_position):
            packet['audit']['passed'] = False
        return packet

    candidate, row = thermal.run_grid(model, 0, 11, plan, cache, provider=provider)
    assert not candidate['certified']
    assert candidate['spectra'] is None
    assert candidate['failed_path_index'] == 2
    assert candidate['total_arrival_rate_s_inverse'] == inflow.total_arrival_rate_s_inverse
    assert candidate['failed_path_rate_s_inverse'] == inflow.rate_s_inverse[2]
    assert candidate['reasons']


def test_final_aggregation_exception_reports_no_fictitious_out_of_range_failed_path(model, plan, cache, monkeypatch):
    aggregate = thermal.stream.aggregate_pump_only_stream

    def fail_after_last_packet(*args, **kwargs):
        aggregate(*args, **kwargs)
        raise ValueError('Deliberate final accumulation failure')

    monkeypatch.setattr(thermal.stream, 'aggregate_pump_only_stream', fail_after_last_packet)
    candidate, row = thermal.run_grid(model, 0, 11, plan, cache, provider=synthetic_packet)
    assert candidate['spectra'] is None and not candidate['certified']
    assert candidate['consumed_path_count'] == candidate['path_count'] == 6
    assert candidate['failed_path_index'] is None
    assert candidate['failed_path_rate_s_inverse'] is None
    assert candidate['last_attempted_path_index'] == 5
    assert 'final' in candidate['failure_stage']


def test_solver_parameters_are_defensively_copied_and_deeply_immutable():
    parameters = {'controls': {'rtol': 1e-9, 'steps': [8, 16]}}
    spec = thermal.SolverSpec('explicit solver', 32., parameters)
    identity = thermal.digest(spec.identity())
    parameters['controls']['rtol'] = 1e-3
    parameters['controls']['steps'][0] = 1
    assert thermal.digest(spec.identity()) == identity
    with pytest.raises(TypeError):
        spec.parameters['controls']['rtol'] = 1e-2
    with pytest.raises(TypeError):
        spec.parameters['controls']['steps'][0] = 2
    assert thermal.digest(spec.identity()) == identity


@pytest.mark.parametrize('method,resolution', [('', 1.), ('method', 0.), ('method', -1.),
    ('method', np.nan), ('method', np.inf)])
def test_solver_spec_requires_named_finite_positive_resolution(method, resolution):
    with pytest.raises((ValueError, TypeError)):
        thermal.SolverSpec(method, resolution, {})


@pytest.mark.parametrize('failure', ['short_primary', 'short_reference', 'same_method', 'mixed_primary',
    'mixed_reference', 'unordered_primary', 'unordered_reference'])
def test_path_plan_requires_two_successive_primary_controls_and_independent_refinement(plan, failure):
    primary, reference = plan.primary, plan.reference
    if failure == 'short_primary':
        primary = primary[:2]
    elif failure == 'short_reference':
        reference = reference[:1]
    elif failure == 'same_method':
        reference = tuple(replace(row, method=primary[0].method) for row in reference)
    elif failure == 'mixed_primary':
        primary = (*primary[:2], replace(primary[-1], method='different-primary'))
    elif failure == 'mixed_reference':
        reference = (reference[0], replace(reference[-1], method='different-reference'))
    elif failure == 'unordered_primary':
        primary = primary[::-1]
    elif failure == 'unordered_reference':
        reference = reference[::-1]
    with pytest.raises((ValueError, TypeError)):
        thermal.PathPlan(primary=primary, reference=reference)


def constant_packet_for_gate(model, path, spec):
    """Explicit constant packet fixture makes grid gates testable at tiny N.

    Each actual physical path still receives five separate provider calculations
    and the native Rb metadata. No atomic solver or factorization is reused.
    """
    packet = synthetic_packet(model, path, spec)
    path_factor = 1+.08*np.sin(path.entry_position_m@[700., 900., 300.])
    path_factor += .06*np.cos(path.velocity_m_s@[.004, -.006, .002])
    for key in ('greater', 'lesser', 'greater_by_source', 'lesser_by_source', 'retarded_response'):
        packet[key] = packet[key]*1e-18/(path.residence_time_s**2*path_factor)
    packet['mean_pulse'] = packet['mean_pulse']*1e-9/(path.residence_time_s*path_factor)
    return packet


@pytest.fixture(scope='module')
def converged_rows(model, plan, tmp_path_factory):
    directory = tmp_path_factory.mktemp('rb-gate-cache')
    source = directory/'constant_packet_provider.py'
    source.write_text('# Explicit constant packet contract fixture, not a physical Rb prediction\n', encoding='utf-8')
    manifest = {str(source): hashlib.sha256(source.read_bytes()).hexdigest()}
    local_cache = thermal.PathCache(directory/'cache', manifest)
    rows = []
    for power, seed in itertools.product((0, 1, 2), (11, 211, 811)):
        candidate, row = thermal.run_grid(model, power, seed, plan, local_cache, provider=constant_packet_for_gate)
        assert candidate['spectra'] is not None and candidate['path_evidence_passed']
        rows.append(row)
    return rows


def clone_rows(rows):
    # Recreate arrays and the typed analysis axis, including mapping-proxy data.
    return thermal.decode(thermal.encode(rows))


def rebind_synthetic_row(row):
    """Bind an explicitly modified unit-test table, never an audit artifact."""
    row['candidate']['candidate_digest'] = thermal.stream._spectrum_digest(row['candidate'])
    row['content_digest'] = thermal.row_digest(row)


def scale_synthetic_row(row, scale):
    row['candidate']['spectra'] = {key: scale*value for key, value in row['candidate']['spectra'].items()}
    rebind_synthetic_row(row)


def test_grid_gate_remeasures_two_edges_and_all_directed_pairs_at_both_recent_grids(converged_rows):
    result = thermal.convergence_gate(list(reversed(converged_rows)))
    assert result['passed'] and not result['reasons']
    assert len(result['comparisons']) == 18
    primary = [row for row in result['comparisons'] if row['kind'] == 'ensemble_refinement']
    independent = [row for row in result['comparisons'] if row['kind'] == 'independent_scramble']
    assert len(primary) == 6 and len(independent) == 12
    assert {(row['reference_grid'][0], row['candidate_grid'][0], row['candidate_grid'][1]) for row in primary} == {
        (lo, hi, seed) for lo, hi in ((0, 1), (1, 2)) for seed in (11, 211, 811)}
    assert {(row['candidate_grid'][0], row['candidate_grid'][1], row['reference_grid'][1]) for row in independent} == {
        (power, one, two) for power in (1, 2) for one, two in itertools.permutations((11, 211, 811), 2)}
    assert all(set(row['errors']) == set(STREAM_KEYS) for row in result['comparisons'])


def test_passing_last_refinement_does_not_erase_an_earlier_failed_required_edge(converged_rows):
    rows = clone_rows(converged_rows)
    for row in rows:
        if row['power'] == 0:
            scale_synthetic_row(row, .8)
    result = thermal.convergence_gate(rows)
    assert not result['passed']
    assert all(row['passed'] for row in result['comparisons']
               if row['kind'] == 'ensemble_refinement' and row['candidate_grid'][0] == 2)


def test_independent_non_nominal_seed_pair_is_checked_even_when_each_matches_nominal(converged_rows):
    rows = clone_rows(converged_rows)
    for row in rows:
        scale_synthetic_row(row, {11: 1., 211: .97, 811: 1.03}[row['seed']])
    result = thermal.convergence_gate(rows)
    assert not result['passed']
    assert all(row['passed'] for row in result['comparisons'] if row['kind'] == 'ensemble_refinement')
    scrambles = [row for row in result['comparisons'] if row['kind'] == 'independent_scramble']
    assert all(row['passed'] for row in scrambles if 11 in (row['candidate_grid'][1], row['reference_grid'][1]))
    assert any(not row['passed'] for row in scrambles if {row['candidate_grid'][1], row['reference_grid'][1]} == {211, 811})


@pytest.mark.parametrize('key', STREAM_KEYS)
def test_each_source_rf_number_and_response_output_controls_the_gate(converged_rows, key):
    rows = clone_rows(converged_rows)
    changed = rows[-1]
    changed['candidate']['spectra'][key] *= 1.2
    rebind_synthetic_row(changed)
    result = thermal.convergence_gate(rows)
    assert not result['passed']
    assert any(row['errors'][key] > .05 for row in result['comparisons'])


def test_dim_source_error_cannot_hide_under_a_bright_source_or_other_rf(converged_rows):
    rows = clone_rows(converged_rows)
    for row in rows:
        for key in ('greater_by_source', 'lesser_by_source'):
            row['candidate']['spectra'][key][0] *= 1e8
        for key in ('greater', 'lesser'):
            row['candidate']['spectra']['internal_'+key] = row['candidate']['spectra'][key+'_by_source'].sum(axis=0)
            row['candidate']['spectra'][key] = row['candidate']['spectra']['internal_'+key]+row['candidate']['spectra']['poisson_number']
        rebind_synthetic_row(row)
    changed = rows[-1]
    changed['candidate']['spectra']['greater_by_source'][-1, -1] *= 1.2
    changed['candidate']['spectra']['internal_greater'] = changed['candidate']['spectra']['greater_by_source'].sum(axis=0)
    changed['candidate']['spectra']['greater'] = changed['candidate']['spectra']['internal_greater']+changed['candidate']['spectra']['poisson_number']
    rebind_synthetic_row(changed)
    result = thermal.convergence_gate(rows)
    assert not result['passed']
    assert max(row['errors']['greater_by_source'] for row in result['comparisons']) > .05


@pytest.mark.parametrize('failure', ['missing_seed', 'duplicate', 'missing_middle_grid', 'only_two_grids',
    'wrong_plan', 'wrong_model', 'incomplete_path_count', 'failed_ledger', 'missing_primary_metric',
    'failed_independent_metric', 'negative_reference_metric'])
def test_missing_or_failed_grid_or_path_evidence_is_rejected(converged_rows, failure):
    rows = clone_rows(converged_rows)
    if failure == 'missing_seed':
        rows.pop()
    elif failure == 'duplicate':
        rows.append(rows[-1])
    elif failure == 'missing_middle_grid':
        rows = [row for row in rows if row['power'] != 1]
    elif failure == 'only_two_grids':
        rows = [row for row in rows if row['power'] != 2]
    else:
        changed = rows[-1]
        if failure == 'wrong_plan':
            changed['plan_digest'] = 'e'*64
        elif failure == 'wrong_model':
            changed['model_digest'] = 'f'*64
        elif failure == 'incomplete_path_count':
            changed['path_ledgers'].pop()
        elif failure == 'failed_ledger':
            changed['path_ledgers'][0]['passed'] = False
        elif failure == 'missing_primary_metric':
            del changed['path_ledgers'][0]['primary_refinement'][0]['mean_outer']
        elif failure == 'failed_independent_metric':
            changed['path_ledgers'][0]['independent_reference']['mean_pulse'] = 1e-3
        elif failure == 'negative_reference_metric':
            changed['path_ledgers'][0]['independent_refinement']['retarded_response'] = -1.
        changed['content_digest'] = thermal.row_digest(changed)
    result = thermal.convergence_gate(rows)
    assert not result['passed'] and result['reasons']


def test_raw_spectrum_tampering_is_rejected_even_if_the_outer_row_digest_is_rebound(converged_rows):
    rows = clone_rows(converged_rows)
    rows[-1]['candidate']['spectra']['retarded_response'][0, 0, 0] *= 1.01
    rows[-1]['content_digest'] = thermal.row_digest(rows[-1])
    result = thermal.convergence_gate(rows)
    assert not result['passed'] and result['reasons']


def test_same_candidate_cannot_be_relabelled_as_an_independent_scramble(converged_rows):
    rows = clone_rows(converged_rows)
    nominal = next(row for row in rows if (row['power'], row['seed']) == (2, 11))
    duplicate = next(row for row in rows if (row['power'], row['seed']) == (2, 211))
    duplicate['candidate'] = nominal['candidate']
    duplicate['content_digest'] = thermal.row_digest(duplicate)
    result = thermal.convergence_gate(rows)
    assert not result['passed']
    assert any('same candidate' in reason or 'Independent' in reason for reason in result['reasons'])


@pytest.mark.parametrize('value', [0., -1., np.nan, np.inf])
def test_stream_comparison_requires_positive_finite_dimensional_scales(converged_rows, value):
    row = converged_rows[-1]
    scales = {**row['comparison_scales'], 'retarded_response': value}
    with pytest.raises(ValueError):
        thermal.stream_errors(row['candidate']['spectra'], row['candidate']['spectra'], scales)


def test_overflowing_matrix_error_is_rejected_instead_of_being_accepted_as_finite(converged_rows):
    row = converged_rows[-1]
    spectra = dict(row['candidate']['spectra'])
    spectra['retarded_response'] = np.full_like(spectra['retarded_response'], 1e308)
    with pytest.raises(ValueError):
        thermal.stream_errors(spectra, row['candidate']['spectra'], row['comparison_scales'])
