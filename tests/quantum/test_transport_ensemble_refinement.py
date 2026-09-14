"""Independent numerical controls for the constant-atom refinement audit."""

from copy import deepcopy
import itertools

import numpy as np
import pytest
from scipy.integrate import quad

from analysis.grand_challenge import transport_ensemble_audit as original
from analysis.grand_challenge import transport_ensemble_refinement as refinement


PATH_KEYS = ('greater', 'lesser', 'greater_by_source', 'lesser_by_source',
             'mean_pulse', 'mean_outer', 'retarded_response')
STREAM_KEYS = ('greater', 'lesser', 'greater_by_source', 'lesser_by_source',
               'poisson_number', 'retarded_response')


def independent_pulse_integrals(tau, omega):
    """Adaptive real integrals on unit age, independent of batched GL/series."""
    pulse, response = [], []
    for frequency in omega:
        phase = frequency*tau
        f = quad(lambda x: np.cos(phase*x), 0., 1., epsabs=2e-13, epsrel=2e-13)[0]
        f += 1j*quad(lambda x: np.sin(phase*x), 0., 1., epsabs=2e-13, epsrel=2e-13)[0]
        r = quad(lambda x: (1-x)*np.cos(phase*x), 0., 1., epsabs=2e-13, epsrel=2e-13)[0]
        r += 1j*quad(lambda x: (1-x)*np.sin(phase*x), 0., 1., epsabs=2e-13, epsrel=2e-13)[0]
        pulse.append(tau*f)
        response.append(tau**2*r)
    return np.array(pulse), np.array(response)


def literal_phase_sum(packet):
    """Average raw ordered moments directly; never import the production mask."""
    keys = ('greater', 'lesser', 'greater_by_source', 'lesser_by_source', 'retarded_response')
    result = {key: np.zeros_like(packet[key]) for key in keys}
    result['poisson_number'] = np.zeros_like(packet['greater'])
    result['mean_pulse'] = np.zeros_like(packet['mean_pulse'])
    for phase in .319+np.arange(4)*np.pi/2:
        d = np.diag(np.exp(1j*np.array([1, -1, -1, 1])*phase))
        for key in keys:
            result[key] += d@packet[key]@d.conj().T/4
        mean = packet['mean_pulse']@d.T
        result['mean_pulse'] += mean/4
        result['poisson_number'] += np.einsum('fi,fj->fij', mean, mean.conj())/4
    for key in ('greater', 'lesser'):
        result[key] += result['poisson_number']
    return result


@pytest.mark.parametrize('order', [None, 48])
def test_batched_age_and_causal_integrals_match_independent_adaptive_reference(order):
    taus = np.array([1e-12, 1e-8, 1e-6, 3e-5])
    pulse, response = refinement.batch_pulse_integrals(taus, order=order)
    assert pulse.shape == response.shape == (len(taus), 3)
    for index, tau in enumerate(taus):
        expected_f, expected_t = independent_pulse_integrals(tau, original.AXIS.omega_rad_s)
        np.testing.assert_allclose(pulse[index]/tau, expected_f/tau, rtol=3e-12, atol=3e-13)
        np.testing.assert_allclose(response[index]/tau**2, expected_t/tau**2, rtol=3e-12, atol=3e-13)
    np.testing.assert_allclose(pulse[:, 0], taus, rtol=3e-15)
    np.testing.assert_allclose(response[:, 0], taus**2/2, rtol=3e-15)


def test_exact_batch_integrals_retain_long_life_oscillations_and_near_dark_pulses():
    taus = np.array([1e-4, 1e-3])
    pulse, response = refinement.batch_pulse_integrals(taus)
    for index, tau in enumerate(taus):
        f, t = independent_pulse_integrals(tau, original.AXIS.omega_rad_s)
        np.testing.assert_allclose(pulse[index]/tau, f/tau, rtol=2e-12, atol=2e-13)
        np.testing.assert_allclose(response[index]/tau**2, t/tau**2, rtol=2e-12, atol=2e-13)
    assert abs(pulse[0, 1]/taus[0]) < 1e-13
    assert abs(response[0, 1]/taus[0]**2) > .1


def test_highly_oscillatory_path_fails_the_unchanged_path_budget_instead_of_being_hidden():
    taus = np.array([1e-3])
    coarse = refinement.batch_pulse_integrals(taus, order=24)
    fine = refinement.batch_pulse_integrals(taus, order=48)
    exact = refinement.batch_pulse_integrals(taus)
    errors = refinement.batch_path_errors(taus, fine, coarse, original.toy_atom())
    analytic = refinement.batch_path_errors(taus, fine, exact, original.toy_atom())
    assert set(errors) == set(analytic) == set(PATH_KEYS)
    assert max(float(value[0]) for value in errors.values()) > original.PATH_BUDGET
    assert max(float(value[0]) for value in analytic.values()) > original.PATH_BUDGET


@pytest.mark.parametrize('taus', [[0.], [-1.], [np.nan], [np.inf], [[1e-6]], []])
def test_invalid_or_empty_residence_arrays_are_rejected(taus):
    with pytest.raises((ValueError, TypeError)):
        refinement.batch_pulse_integrals(taus)


def test_batch_path_errors_match_each_of_seven_full_matrix_metrics():
    atom = original.toy_atom()
    inflow = original.make_inflow(1, 11)
    contract = original.toy_convention(inflow, atom)
    taus = inflow.residence_time_s
    reference = refinement.batch_pulse_integrals(taus)
    # Deliberate finite discrepancy avoids a comparison of roundoff with itself.
    change = np.linspace(2e-5, 6e-5, len(taus))[:, None]
    candidate = (reference[0]*(1+change*(1+.3j)), reference[1]*(1-change*(.2+.7j)))
    actual = refinement.batch_path_errors(taus, candidate, reference, atom)
    assert set(actual) == set(PATH_KEYS)
    for index, tau in enumerate(taus):
        exact = original.toy_packet(inflow.path(index), contract, atom)
        f, t = candidate[0][index], candidate[1][index]
        packet = dict(exact)
        for key in ('greater', 'lesser'):
            packet[key] = np.abs(f[:, None, None])**2*atom[key]
            packet[key+'_by_source'] = packet[key][None]
        packet['mean_pulse'] = f[:, None]*atom['mean']
        packet['retarded_response'] = t[:, None, None]*atom['response']
        expected = original.observed_errors(original.metric_arrays(packet), original.metric_arrays(exact),
            PATH_KEYS, original.metric_scales(tau, atom))
        for key in PATH_KEYS:
            assert actual[key].shape == (len(taus),)
            assert actual[key][index] == pytest.approx(expected[key], rel=3e-10, abs=3e-15)


@pytest.fixture(scope='module')
def same_grid():
    atom = original.toy_atom()
    spectra, row = refinement.factorized_grid(2, 11, atom, chunk_size=5)
    candidate, original_row, _ = original.run_grid(2, 11, atom)
    return atom, spectra, row, candidate, original_row


def test_scalar_factorization_matches_original_complete_complex_stream_matrices(same_grid):
    _, spectra, row, candidate, previous_row = same_grid
    assert row['path_count'] == candidate['path_count'] == 24
    assert row['path_evidence_passed'] and previous_row['path_evidence_passed']
    assert set(spectra) == set(candidate['spectra'])
    for key, expected in candidate['spectra'].items():
        np.testing.assert_allclose(spectra[key], expected, rtol=5e-13, atol=5e-15*np.linalg.norm(expected))
    for kind in ('path_refinement', 'independent_reference'):
        assert set(row['max_path_errors'][kind]) == set(PATH_KEYS)
        for key in PATH_KEYS:
            assert row['max_path_errors'][kind][key] == pytest.approx(
                previous_row['max_path_errors'][kind][key], rel=2e-3, abs=8e-15)


def test_factorized_stream_matches_independent_literal_raw_common_phase_average(same_grid):
    atom, spectra, _, _, _ = same_grid
    inflow = original.make_inflow(2, 11)
    contract = original.toy_convention(inflow, atom)
    expected = None
    for index, rate in enumerate(inflow.rate_s_inverse):
        packet = original.toy_packet(inflow.path(index), contract, atom)
        phase = literal_phase_sum(packet)
        if expected is None:
            expected = {key: np.zeros_like(value) for key, value in phase.items()}
        for key, value in phase.items():
            expected[key] += rate*value
    for key in STREAM_KEYS:
        np.testing.assert_allclose(spectra[key], expected[key], rtol=8e-13, atol=8e-15*np.linalg.norm(expected[key]))
    mean_scale = inflow.mean_occupancy*np.linalg.norm(atom['mean'])*np.sqrt(3)
    assert np.linalg.norm(expected['mean_pulse']) < 1e-14*mean_scale
    assert np.linalg.norm(spectra['poisson_number']) > .1*np.linalg.norm(spectra['greater'])
    np.testing.assert_allclose(spectra['greater']-spectra['internal_greater'], spectra['poisson_number'], rtol=1e-14)


def test_density_enters_exactly_once_and_chunk_boundaries_preserve_the_integral(same_grid):
    atom, base, row, _, _ = same_grid
    doubled, other = refinement.factorized_grid(2, 11, atom, density=2*original.DENSITY_M3, chunk_size=7)
    for key in base:
        np.testing.assert_allclose(doubled[key], 2*base[key], rtol=4e-14, atol=4e-15*np.linalg.norm(base[key]))
    assert other['path_count'] == row['path_count']
    for kind in ('path_refinement', 'independent_reference'):
        for key in PATH_KEYS:
            assert other['max_path_errors'][kind][key] == pytest.approx(row['max_path_errors'][kind][key], abs=8e-15)


def synthetic_gate_rows(same_grid, *, powers=(0, 1, 2), seeds=(11, 211, 811)):
    """Explicitly synthetic equal-integral rows isolate certification guards.

    These rows do not claim that low-order physical Sobol grids converged.
    Actual numerical parity is established by the separate integration tests.
    """
    template = same_grid[2]
    rows = []
    for power, seed in itertools.product(powers, seeds):
        row = deepcopy(template)
        row.update(power=power, seed=seed, path_count=6*2**power,
            evaluated_path_count=6*2**power,
            evaluation_method='Explicit synthetic equal-integral table for unit-testing evidence gates')
        row['content_digest'] = refinement.row_digest(row)
        rows.append(row)
    return rows


def rescale_synthetic_row(row, scale):
    row['spectra'] = {key: scale*value for key, value in row['spectra'].items()}
    row['content_digest'] = refinement.row_digest(row)


def test_gate_requires_two_refinements_per_seed_and_every_directed_pair_on_both_recent_grids(same_grid):
    rows = synthetic_gate_rows(same_grid)
    result = refinement.convergence_gate(rows[::-1], powers=(0, 1, 2, 3))
    assert result['passed'] and not result['reasons']
    assert len(result['refinements']) == 6
    assert len(result['scrambles']) == 12
    edges = {(row['reference_power'], row['candidate_power'], row['candidate_seed'])
             for row in result['refinements']}
    assert edges == {(lo, hi, seed) for lo, hi in ((0, 1), (1, 2)) for seed in (11, 211, 811)}
    pairs = {(row['candidate_power'], row['candidate_seed'], row['reference_seed'])
             for row in result['scrambles']}
    assert pairs == {(power, first, second) for power in (1, 2)
                     for first, second in itertools.permutations((11, 211, 811), 2)}
    assert all(set(row['errors']) == set(STREAM_KEYS) for row in result['refinements']+result['scrambles'])
    assert all(row['budget'] == .05 for row in result['refinements']+result['scrambles'])


def test_recent_two_refinements_cannot_be_replaced_by_only_a_passing_last_edge(same_grid):
    rows = synthetic_gate_rows(same_grid)
    for row in rows:
        if row['power'] == 0:
            rescale_synthetic_row(row, .8)
    result = refinement.convergence_gate(rows, powers=(0, 1, 2))
    assert not result['passed']
    assert all(row['passed'] for row in result['refinements'] if row['candidate_power'] == 2)
    assert any(not row['passed'] for row in result['refinements'] if row['candidate_power'] == 1)


def test_gate_uses_last_three_grids_of_a_complete_plan_prefix(same_grid):
    rows = synthetic_gate_rows(same_grid, powers=(0, 1, 2, 3))
    for row in rows:
        if row['power'] == 0:
            rescale_synthetic_row(row, .5)
    result = refinement.convergence_gate(rows, powers=(0, 1, 2, 3))
    assert result['passed']
    assert {row['reference_power'] for row in result['refinements']} == {1, 2}


def test_independent_comparisons_are_required_at_the_penultimate_grid_too(same_grid):
    rows = synthetic_gate_rows(same_grid)
    for row in rows:
        if row['power'] < 2:
            rescale_synthetic_row(row, {11: 1., 211: .97, 811: 1.03}[row['seed']])
    result = refinement.convergence_gate(rows, powers=(0, 1, 2))
    assert not result['passed']
    assert all(row['passed'] for row in result['refinements'])
    assert all(row['passed'] for row in result['scrambles'] if row['candidate_power'] == 2)
    assert any(not row['passed'] for row in result['scrambles'] if row['candidate_power'] == 1)


def test_non_nominal_seed_pair_cannot_hide_behind_individual_agreement_with_nominal(same_grid):
    rows = synthetic_gate_rows(same_grid)
    for row in rows:
        rescale_synthetic_row(row, {11: 1., 211: .97, 811: 1.03}[row['seed']])
    result = refinement.convergence_gate(rows, powers=(0, 1, 2))
    assert not result['passed']
    assert all(row['passed'] for row in result['refinements'])
    nominal_pairs = [row for row in result['scrambles'] if 11 in (row['candidate_seed'], row['reference_seed'])]
    assert nominal_pairs and all(row['passed'] for row in nominal_pairs)
    assert any(not row['passed'] for row in result['scrambles']
               if {row['candidate_seed'], row['reference_seed']} == {211, 811})


@pytest.mark.parametrize('key', STREAM_KEYS)
def test_each_actual_stream_metric_is_required_even_when_number_moments_are_exact(same_grid, key):
    rows = synthetic_gate_rows(same_grid)
    changed = rows[-1]
    changed['spectra'][key] = 1.2*changed['spectra'][key]
    changed['occupancy_over_nV'] = 1.
    changed['mean_occupancy'] = changed['equilibrium_nV']
    changed['content_digest'] = refinement.row_digest(changed)
    result = refinement.convergence_gate(rows, powers=(0, 1, 2))
    assert not result['passed'] and result['reasons']
    assert any(row['errors'][key] > .05 for row in result['refinements']+result['scrambles'])


@pytest.mark.parametrize('failure', ['missing_seed', 'duplicate', 'unplanned_seed', 'unplanned_power',
    'missing_middle_grid', 'two_grids_only', 'missing_path_metric', 'failed_path', 'incomplete_count',
    'negative_path_error', 'failed_path_error', 'missing_closed_form_metric', 'failed_closed_form_error'])
def test_incomplete_or_failed_numerical_evidence_cannot_certify(same_grid, failure):
    rows = synthetic_gate_rows(same_grid)
    if failure == 'missing_seed':
        rows.pop()
    elif failure == 'duplicate':
        rows.append(deepcopy(rows[-1]))
    elif failure == 'missing_middle_grid':
        rows = [row for row in rows if row['power'] != 1]
    elif failure == 'two_grids_only':
        rows = [row for row in rows if row['power'] < 2]
    else:
        row = rows[-1]
        if failure == 'unplanned_seed':
            row['seed'] = 812
        elif failure == 'unplanned_power':
            row['power'] = 3
        elif failure == 'missing_path_metric':
            del row['max_path_errors']['independent_reference']['mean_outer']
        elif failure == 'failed_path':
            row['failed_path_count'] = 1
        elif failure == 'incomplete_count':
            row['evaluated_path_count'] -= 1
        elif failure == 'negative_path_error':
            row['max_path_errors']['path_refinement']['mean_pulse'] = -1.
        elif failure == 'failed_path_error':
            row['max_path_errors']['path_refinement']['mean_pulse'] = 2e-9
        elif failure == 'missing_closed_form_metric':
            del row['closed_form_stream_errors']['retarded_response']
        elif failure == 'failed_closed_form_error':
            row['closed_form_stream_errors']['retarded_response'] = 4e-12
        row['content_digest'] = refinement.row_digest(row)
    result = refinement.convergence_gate(rows, powers=(0, 1, 2))
    assert not result['passed'] and result['reasons']


def test_content_digest_binds_spectra_evidence_and_path_metadata(same_grid):
    for kind in ('spectra', 'evidence', 'path_count'):
        rows = synthetic_gate_rows(same_grid)
        changed = rows[-1]
        if kind == 'spectra':
            changed['spectra']['retarded_response'] = 1.01*changed['spectra']['retarded_response']
        elif kind == 'evidence':
            changed['max_path_errors']['independent_reference']['mean_outer'] = 1e-13
        else:
            changed['path_count'] += 1
        result = refinement.convergence_gate(rows, powers=(0, 1, 2))
        assert not result['passed']
        assert any('digest' in reason.lower() for reason in result['reasons'])


def test_mixed_physical_densities_cannot_pass_even_when_all_matrix_differences_are_small(same_grid):
    rows = synthetic_gate_rows(same_grid)
    rows[-1]['density_m3'] *= 1.001
    rows[-1]['content_digest'] = refinement.row_digest(rows[-1])
    result = refinement.convergence_gate(rows, powers=(0, 1, 2))
    assert not result['passed']
    assert any('density' in reason.lower() or 'densities' in reason.lower() for reason in result['reasons'])


@pytest.mark.parametrize('powers,seeds', [((0, 1), (11, 211, 811)), ((0, 1, 1), (11, 211, 811)),
    ((0, 2, 1), (11, 211, 811)), ((0, 1, 2), (11, 211)), ((0, 1, 2), (11, 211, 211))])
def test_plan_itself_requires_three_distinct_seeds_and_three_increasing_grids(same_grid, powers, seeds):
    result = refinement.convergence_gate(synthetic_gate_rows(same_grid), powers=powers, seeds=seeds)
    assert not result['passed'] and result['reasons']


def test_one_failed_path_is_counted_and_retained_in_unrenormalized_diagnostic(same_grid, monkeypatch):
    atom = same_grid[0]
    inflow = original.make_inflow(0, 11)
    unmodified = refinement.batch_pulse_integrals

    def inaccurate(taus, *, order=None):
        f, t = unmodified(taus, order=order)
        if order == 48:
            f = f.copy()
            f[0] *= 1.3
        return f, t

    monkeypatch.setattr(refinement, 'batch_pulse_integrals', inaccurate)
    spectra, row = refinement.factorized_grid(0, 11, atom, chunk_size=10)
    f, _ = inaccurate(inflow.residence_time_s, order=48)
    expected_weight = np.sum(inflow.rate_s_inverse[:, None]*abs(f)**2, axis=0)
    assert not row['path_evidence_passed']
    assert row['failed_path_count'] == 1
    assert row['path_count'] == row['evaluated_path_count'] == 6
    assert row['total_arrival_rate_s_inverse'] == inflow.total_arrival_rate_s_inverse
    np.testing.assert_allclose(row['pulse_weight_s'], expected_weight, rtol=3e-15)
    # The diagnostic includes the failed path with its unchanged physical rate.
    np.testing.assert_allclose(spectra['poisson_number'][:, 0, 0],
                              expected_weight*abs(atom['mean'][0])**2, rtol=3e-15)


def test_factorization_rejects_a_changed_atomic_model(same_grid):
    atom = deepcopy(same_grid[0])
    atom['response'] *= 1.01
    with pytest.raises(ValueError, match='fixture|atom|Factorization'):
        refinement.factorized_grid(0, 11, atom)


@pytest.mark.parametrize('existing', ['report', 'plot'])
def test_cli_rejects_existing_artifacts_before_expensive_calculation(tmp_path, monkeypatch, existing):
    report, plot = tmp_path/'report.json', tmp_path/'plot.png'
    artifact = report if existing == 'report' else plot
    artifact.write_bytes(b'user-owned original')

    def forbidden():
        pytest.fail('preflight must reject before building a report')

    monkeypatch.setattr(refinement, 'build_report', forbidden)
    with pytest.raises(FileExistsError):
        refinement.main(['--output', str(report), '--plot', str(plot)])
    assert artifact.read_bytes() == b'user-owned original'
    assert not (plot if existing == 'report' else report).exists()


def test_cli_rejects_identical_new_report_and_plot_targets(tmp_path, monkeypatch):
    artifact = tmp_path/'unwritten'

    def forbidden():
        pytest.fail('matching output paths must reject before building a report')

    monkeypatch.setattr(refinement, 'build_report', forbidden)
    with pytest.raises(ValueError):
        refinement.main(['--output', str(artifact), '--plot', str(artifact)])
    assert not artifact.exists()


def test_unstable_source_hashes_prevent_any_immutable_output(tmp_path, monkeypatch):
    artifact = tmp_path/'unwritten.json'
    snapshots = iter([{'source.py': 'before'}, {'source.py': 'after'}])
    monkeypatch.setattr(refinement, 'source_hashes', lambda: next(snapshots))
    monkeypatch.setattr(refinement, 'build_report', lambda: {
        'expected_controls_passed': True, 'implementation_controls_passed': True})
    with pytest.raises(RuntimeError, match='[Ss]ource|[Hh]ash|changed|stable'):
        refinement.main(['--output', str(artifact)])
    assert not artifact.exists()
