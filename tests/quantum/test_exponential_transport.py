"""Opt-in exponential transport: exact maps and finite-atom convergence.

These controls test prescribed atomic paths, not physical Rb/Maxwell convergence.
The matrix ODE and raw-density QRT references do not use the spectral source map.
"""

from dataclasses import replace

import numpy as np
import pytest
from scipy.integrate import quad, solve_ivp
from scipy.linalg import expm

from gabes import core
from gabes.quantum.diffusion import traceless_hermitian_basis
from gabes.quantum.reservoirs import CollapseChannel, ExplicitReservoirs
from gabes.quantum.segmented_transport import segmented_wavepacket
from gabes.quantum.smooth_transport import smooth_wavepacket
from analysis.grand_challenge.reference.exponential_transport import (
    _apply_source_map,
    _density_modes,
    _divided_exponential,
    _drift_modes,
    _source_map,
    exponential_wavepacket,
)
from analysis.grand_challenge.reference.smooth_qrt import smooth_qrt


@pytest.fixture(autouse=True)
def single_blas_thread():
    with core.blas_single_thread():
        yield


def _physical_problem(n=2):
    rng = np.random.default_rng(831+n)
    h = rng.normal(size=(2, n, n))+1j*rng.normal(size=(2, n, n))
    h = (h+h.conj().swapaxes(-1, -2))/3
    x = rng.normal(size=(n, n))+1j*rng.normal(size=(n, n))
    rho = x@x.conj().T
    rho /= rho.trace()
    lower = np.diag(np.ones(n-1), 1)
    reservoirs = ExplicitReservoirs(n, (
        CollapseChannel('decay', lower, 'decay', 'declared test atom'),
        CollapseChannel('pump', .3*lower.T, 'excitation', 'declared test atom'),
        CollapseChannel('dephasing', .2*np.diag(np.arange(n)),
                        'dephasing', 'declared test atom'),
    ))
    ops = np.array([lower+.17j*np.eye(n), x/3])
    return dict(h0=h[0], h1=h[1],
                envelope=lambda t: np.exp(-((t-.3)/.7)**2),
                reservoirs=reservoirs, boundary_state=rho, duration_s=1.1,
                readouts=ops, frequencies_rad_s=np.array([[.2, -.7], [.5, .3]]),
                drives=ops.conj().swapaxes(-1, -2),
                drive_frequencies_rad_s=np.array([[.6, -.4], [.1, .9]]))


def _matrix_problem():
    # Distinct stable eigenvalues, complex nonnormal eigenvectors, two sources
    # and deliberately unequal orderings detect transpose/conjugation mistakes.
    rng = np.random.default_rng(1957)
    aug = np.diag([-.4+.2j, -.9-.3j, -1.7+.6j])
    aug[0, 1], aug[1, 2] = .4+.3j, -.2+.1j
    density = np.diag([-.2+.1j, -.6-.4j, -1.2+.5j, -1.9-.2j])
    density[0, 1], density[1, 3] = .2-.3j, .1+.2j
    dmap = rng.normal(size=(4, 9, 4))+1j*rng.normal(size=(4, 9, 4))
    cov = rng.normal(size=(2, 3, 3, 3))+1j*rng.normal(size=(2, 3, 3, 3))
    rho = np.array([[.7, .1+.08j], [.1-.08j, .3]])
    return aug, density, dmap, cov, rho


def _matrix_ode(aug, density, dmap, cov, rho, dt):
    """Direct dC/dt=A C+C A†+D(rho), with no Kronecker/eigenbasis lift."""
    nr, dim, ns = rho.size, len(aug), cov.shape[1]

    def rhs(t, y):
        state = y[:nr]
        current = y[nr:].reshape(cov.shape)
        rate = aug@current+current@aug.conj().T
        for ordering in range(2):
            for source in range(1, ns):
                index = ordering*(ns-1)+source-1
                rate[ordering, source] += (dmap[index]@state).reshape(dim, dim)
        return np.r_[density@state, rate.ravel()]

    result = solve_ivp(rhs, (0., dt), np.r_[rho.ravel(), cov.ravel()],
                       method='DOP853', rtol=2e-13, atol=2e-15)
    assert result.success, result.message
    return result.y[nr:, -1].reshape(cov.shape)


@pytest.mark.parametrize('cached_density', [False, True])
@pytest.mark.parametrize('cached_drift', [False, True])
def test_spectral_and_forced_block_preserve_every_source_and_order(cached_density, cached_drift):
    aug, density, dmap, cov, rho = _matrix_problem()
    original_cov, original_rho = cov.copy(), rho.copy()
    modes = _density_modes(density, dmap) if cached_density else None
    drift = _drift_modes(aug) if cached_drift else None
    spectral = _source_map(aug, density, dmap, .73, 1e6, modes, drift)
    block = _source_map(aug, density, dmap, .73, 1, modes, drift)
    assert len(spectral) == len(block) == 6
    assert spectral[0] == 'spectral'
    assert block[0] == 'block'
    assert 1 < spectral[-1] < 1e6
    # Reuse one map for different entering densities, as the public solver does.
    for state in (rho, np.array([[.4, -.13j], [.13j, .6]])):
        expected = _matrix_ode(aug, density, dmap, cov, state, .73)
        for mapping in (spectral, block):
            actual = _apply_source_map(mapping, cov, state)
            np.testing.assert_allclose(actual, expected, rtol=3e-12, atol=3e-13)
            homogeneous = expm(aug*.73)
            np.testing.assert_allclose(actual[:, 0],
                                       homogeneous@cov[:, 0]@homogeneous.conj().T,
                                       rtol=3e-12, atol=3e-13)
    np.testing.assert_array_equal(cov, original_cov)
    np.testing.assert_array_equal(rho, original_rho)


@pytest.mark.parametrize('defective', ['augmented', 'density', 'both'])
def test_defective_jordan_generators_use_exact_block_fallback(defective):
    aug, density, dmap, cov, rho = _matrix_problem()
    if defective in ('augmented', 'both'):
        aug = (-.6+.2j)*np.eye(3)+np.diag([1.+.3j, .7-.2j], 1)
    if defective in ('density', 'both'):
        density = (-.4-.1j)*np.eye(4)+np.diag([1., .7j, -.3], 1)
    mapping = _source_map(aug, density, dmap, .9, 1e6)
    assert mapping[0] == 'block'
    assert mapping[-1] > 1e6
    actual = _apply_source_map(mapping, cov, rho)
    assert np.isfinite(actual).all()
    np.testing.assert_allclose(actual, _matrix_ode(aug, density, dmap, cov, rho, .9),
                               rtol=3e-12, atol=3e-13)


@pytest.mark.parametrize('max_condition', [1, 1e6])
def test_source_map_retains_signed_noise_without_clipping(max_condition):
    aug = np.array([[-.4+.2j, .3j], [0., -.8-.1j]])
    density = np.diag([0., -.2, -.5, -.9]).astype(complex)
    # A tiny negative diagnostic must stay signed, even with no jump sources.
    cov = np.broadcast_to(-1e-14*np.eye(2), (2, 1, 2, 2)).astype(complex).copy()
    mapping = _source_map(aug, density, np.empty((0, 4, 4)), .7, max_condition)
    assert mapping[0] == ('block' if max_condition == 1 else 'spectral')
    actual = _apply_source_map(mapping, cov, np.eye(2)/2)
    propagator = expm(aug*.7)
    np.testing.assert_allclose(actual, propagator@cov@propagator.conj().T,
                               rtol=3e-14, atol=0.)
    assert np.linalg.eigvalsh(actual).max() < 0.


@pytest.mark.parametrize('dt', [0., .73, 2.])
def test_divided_exponential_broadcast_exact_and_near_resonance(dt):
    left = np.array([0., -.4+.7j, -2.-.3j])[:, None]
    offsets = np.array([0., 1e-13+2e-13j, -1e-8j, .000999, .001001, .5+.3j])
    right = left+offsets
    actual = _divided_exponential(left, right, dt)
    expected = np.empty(right.shape, complex)
    for i, j in np.ndindex(expected.shape):
        # Even scipy's triangular expm can lose digits for almost repeated
        # complex diagonals. Integrate the defining smooth scalar convolution.
        integrand = lambda s: np.exp(left[i, 0]*(dt-s)+right[i, j]*s)
        expected[i, j] = (
            quad(lambda s: integrand(s).real, 0., dt, epsabs=2e-14, epsrel=2e-14)[0]
            +1j*quad(lambda s: integrand(s).imag, 0., dt, epsabs=2e-14, epsrel=2e-14)[0])
    np.testing.assert_allclose(actual, expected, rtol=5e-13, atol=2e-15)
    np.testing.assert_allclose(actual[:, 0], dt*np.exp(left[:, 0]*dt),
                               rtol=2e-15, atol=0.)
    scalar = _divided_exponential(-.4+.7j, -.4+.7j, dt)
    np.testing.assert_allclose(scalar, dt*np.exp((-.4+.7j)*dt), rtol=2e-15)


def test_divided_exponential_strong_decay_never_forms_zero_times_infinity():
    left = np.array([-1200., -.25, -1400., -1200.+.3j, -1600.])
    right = np.array([-.25, -1200., -1300., -1200.+.3j+1e-10, 0.])
    with np.errstate(over='raise', invalid='raise', divide='raise'):
        actual = _divided_exponential(left, right, 2.)
    expected = np.array([np.exp(-.5)/1199.75, np.exp(-.5)/1199.75,
                         0., 0., 1/1600])
    assert np.isfinite(actual).all()
    np.testing.assert_allclose(actual, expected, rtol=2e-15, atol=0.)


def _assert_quantum_consistency(result):
    states = result['states']
    np.testing.assert_allclose(states, states.conj().swapaxes(-1, -2), atol=3e-12)
    np.testing.assert_allclose(np.trace(states, axis1=-2, axis2=-1), 1., atol=3e-12)
    assert np.linalg.eigvalsh(states).min() > -3e-12
    atomic = result['atomic_covariance_by_source']
    np.testing.assert_allclose(atomic[:, 1], atomic[:, 0].swapaxes(-1, -2),
                               rtol=2e-10, atol=3e-12)
    basis = traceless_hermitian_basis(states.shape[-1])
    for state, sources in zip(states, atomic):
        means = np.array([np.trace(state@f) for f in basis])
        raw = np.array([[np.trace(state@f@g) for g in basis] for f in basis])
        expected = raw-means[:, None]*means[None, :]
        np.testing.assert_allclose(sources[0].sum(axis=0), expected, atol=4e-12)
        # Direct <[Fi,Fj]> fixes the unsymmetrized ordering and its sign.
        np.testing.assert_allclose((sources[0]-sources[1]).sum(axis=0),
                                   raw-raw.T, atol=4e-12)
    for order in ('greater', 'lesser'):
        sources = result[order+'_by_source']
        np.testing.assert_allclose(sources, sources.conj().swapaxes(-1, -2), atol=4e-12)
        assert np.linalg.eigvalsh(sources).min() > -4e-12
        np.testing.assert_allclose(sources.sum(axis=0), result[order], atol=3e-13)
    assert result['audit']['passed']


@pytest.mark.parametrize('order', [2, 4])
@pytest.mark.parametrize('max_condition', [1, 1e6])
def test_constant_physical_wavepacket_segmented_parity_drives_rf_and_seconds(order, max_condition):
    kw = _physical_problem()
    kw['envelope'] = lambda t: .4
    options = dict(segments=4, order=order, max_condition=max_condition)
    result = exponential_wavepacket(**kw, **options)
    reference = segmented_wavepacket(
        [kw['h0']+.4*kw['h1']], kw['reservoirs'], kw['boundary_state'],
        [kw['duration_s']], [kw['readouts']], kw['frequencies_rad_s'],
        drives=[kw['drives']], drive_frequencies_rad_s=kw['drive_frequencies_rad_s'])
    factor = 1e-6
    scaled = {**kw, 'h0': kw['h0']/factor, 'h1': kw['h1']/factor,
              'duration_s': kw['duration_s']*factor,
              'frequencies_rad_s': kw['frequencies_rad_s']/factor,
              'drive_frequencies_rad_s': kw['drive_frequencies_rad_s']/factor,
              'drives': kw['drives']/np.sqrt(factor),
              'reservoirs': ExplicitReservoirs(2, tuple(
                  replace(c, operator=c.operator/np.sqrt(factor))
                  for c in kw['reservoirs'].channels))}
    seconds = exponential_wavepacket(**scaled, **options)
    for key, power in [('greater', 2), ('lesser', 2), ('greater_by_source', 2),
                       ('lesser_by_source', 2), ('mean_pulse', 1),
                       ('retarded_response', 1.5), ('exit_state', 0)]:
        np.testing.assert_allclose(result[key], reference[key], rtol=3e-10, atol=3e-12)
        # Compare after undoing SI scaling: a loose s^2 absolute floor hides bugs.
        np.testing.assert_allclose(seconds[key]/factor**power, result[key],
                                   rtol=3e-10, atol=3e-12)
    assert result['source_names'] == reference['source_names'] == (
        'atomic_inflow', 'jump:decay', 'jump:pump', 'jump:dephasing')
    np.testing.assert_allclose(result['mean_outer'],
                               result['mean_pulse'][:, :, None]*result['mean_pulse'][:, None, :].conj())
    np.testing.assert_allclose(result['sample_ages_s'], np.linspace(0., kw['duration_s'], 5))
    assert result['states'].shape == (5, 2, 2)
    _assert_quantum_consistency(result)
    numerics = result['numerics']
    assert numerics['order'] == order and numerics['segments'] == 4
    assert numerics['substeps'] == (8 if order == 4 else 4)
    maps = numerics['spectral_maps']+numerics['block_fallback_maps']
    assert maps > 0
    assert maps == numerics['spectral_response_maps']+numerics['small_response_exponentials']
    assert numerics['cached_exponential_reuses'] > 0
    if max_condition == 1:
        assert numerics['spectral_maps'] == 0
        assert numerics['block_fallback_maps'] > 0
        assert numerics['small_response_exponentials'] > 0
        assert numerics['spectral_response_maps'] == 0
    else:
        assert numerics['spectral_maps'] > 0
        assert numerics['block_fallback_maps'] == 0
        assert numerics['small_response_exponentials'] == 0
        assert numerics['spectral_response_maps'] > 0


def test_rf_batch_shares_density_and_preserves_row_order():
    kw = _physical_problem()
    batch = exponential_wavepacket(**kw, segments=4)
    for row in (1, 0):
        single = exponential_wavepacket(**{
            **kw, 'frequencies_rad_s': kw['frequencies_rad_s'][row:row+1],
            'drive_frequencies_rad_s': kw['drive_frequencies_rad_s'][row:row+1]}, segments=4)
        for key in ('greater', 'lesser', 'mean_pulse', 'retarded_response'):
            np.testing.assert_allclose(single[key][0], batch[key][row], rtol=2e-11, atol=2e-12)
        for key in ('greater_by_source', 'lesser_by_source'):
            np.testing.assert_allclose(single[key][:, 0], batch[key][:, row], rtol=2e-11, atol=2e-12)
        np.testing.assert_allclose(single['states'], batch['states'], atol=2e-12)


@pytest.mark.parametrize('max_condition', [1, 1e6])
def test_named_jump_split_phase_and_permutation_keep_each_order(max_condition):
    kw = _physical_problem()
    channel = kw['reservoirs'].channels[0]
    kw['reservoirs'] = ExplicitReservoirs(2, (channel,))
    options = dict(segments=2, max_condition=max_condition)
    base = exponential_wavepacket(**kw, **options)
    channels = (
        replace(channel, name='b', operator=.8*np.exp(-.29j)*channel.operator),
        replace(channel, name='a', operator=.6*np.exp(.37j)*channel.operator))
    split = exponential_wavepacket(**{**kw, 'reservoirs': ExplicitReservoirs(2, channels)}, **options)
    assert split['source_names'] == ('atomic_inflow', 'jump:b', 'jump:a')
    for key in ('greater', 'lesser', 'mean_pulse', 'retarded_response', 'exit_state'):
        np.testing.assert_allclose(split[key], base[key], rtol=3e-11, atol=3e-12)
    for key in ('greater_by_source', 'lesser_by_source'):
        for index, weight in ((0, 1.), (1, .64), (2, .36)):
            np.testing.assert_allclose(split[key][index], weight*base[key][min(index, 1)],
                                       rtol=3e-11, atol=3e-12)
    _assert_quantum_consistency(split)


def test_ground_state_vacuum_ordering_matches_analytic_source_integrals():
    gamma, resonance, frequency, duration = .8, 1.3, .6, 1.7
    lower = np.array([[0., 1.], [0., 0.]])
    reservoirs = ExplicitReservoirs(2, (
        CollapseChannel('radiation', np.sqrt(gamma)*lower, 'decay', 'analytic control'),))
    result = exponential_wavepacket(np.diag([0., resonance]), np.zeros((2, 2)),
        lambda t: 0., reservoirs, np.diag([1., 0.]), duration, lower[None],
        [[frequency]], segments=4)
    z = -gamma/2+1j*(frequency-resonance)
    pulse = lambda t: np.expm1(z*t)/z
    boundary = abs(pulse(duration))**2
    bath = gamma*quad(lambda t: abs(pulse(t))**2, 0., duration, epsabs=1e-13)[0]
    np.testing.assert_allclose(result['greater_by_source'][:, 0, 0, 0],
                               [boundary, bath], rtol=3e-12, atol=2e-14)
    np.testing.assert_allclose(result['greater'],
                               2*np.real((np.expm1(z*duration)-z*duration)/z**2), rtol=3e-12)
    np.testing.assert_allclose(result['lesser_by_source'], 0., atol=3e-14)
    np.testing.assert_allclose(result['mean_pulse'], 0., atol=3e-14)
    assert result['retarded_response'].shape == (1, 1, 0)
    _assert_quantum_consistency(result)


@pytest.mark.parametrize('max_condition', [1, 1e6])
def test_constant_complex_response_with_unequal_rf_and_resonant_drive(max_condition):
    lower = np.array([[0., 1.], [0., 0.]])
    resonance, duration = 1.3, 1.1
    readout_coefficients = np.array([1., .4+.3j])
    drive_coefficients = np.array([.6-.2j, .3+.8j])
    frequencies = np.array([[.4, -.3], [0., 1.3], [1.3+1e-11, -.7]])
    drive_frequencies = np.array([[-.7, .9], [1.3, 0.], [1.3+1e-12, -.2]])
    kw = dict(h0=np.diag([0., resonance]), h1=np.zeros((2, 2)),
              envelope=lambda t: 0., reservoirs=ExplicitReservoirs(2, ()),
              boundary_state=np.diag([1., 0.]), duration_s=duration,
              readouts=readout_coefficients[:, None, None]*lower,
              frequencies_rad_s=frequencies,
              drives=drive_coefficients[:, None, None]*lower.T,
              drive_frequencies_rad_s=drive_frequencies)
    result = exponential_wavepacket(**kw, segments=4, max_condition=max_condition)
    expected = np.empty((3, 2, 2), complex)
    for row, port, drive in np.ndindex(expected.shape):
        w, nu = frequencies[row, port], drive_frequencies[row, drive]

        def integrand(t):
            # delta rho_eg(t)=-i*d*exp(-i*resonance*t)*integral exp(i*(resonance-nu)*s) ds.
            difference = resonance-nu
            inner = t if difference == 0 else np.expm1(1j*difference*t)/(1j*difference)
            return np.exp(1j*(w-resonance)*t)*inner

        integral = (
            quad(lambda t: integrand(t).real, 0., duration, epsabs=1e-13)[0]
            +1j*quad(lambda t: integrand(t).imag, 0., duration, epsabs=1e-13)[0])
        expected[row, port, drive] = -1j*readout_coefficients[port]*drive_coefficients[drive]*integral
    np.testing.assert_allclose(result['retarded_response'], expected, rtol=3e-11, atol=3e-13)
    assert np.linalg.norm(expected.imag) > .1
    assert result['source_names'] == ('atomic_inflow',)
    np.testing.assert_allclose(result['lesser'], 0., atol=3e-13)
    _assert_quantum_consistency(result)


@pytest.mark.parametrize('order', [2, 4])
def test_sampled_history_retains_macro_endpoints_and_full_streaming_audit(order):
    kw = {**_physical_problem(), 'segments': 32, 'order': order}
    full = exponential_wavepacket(**kw, sample_count=None)
    assert full['states'].shape == (33, 2, 2)
    for options, count in (({}, 17), ({'sample_count': 2}, 2),
                           ({'sample_count': 5}, 5), ({'sample_count': 100}, 33)):
        sampled = exponential_wavepacket(**kw, **options)
        ages = sampled['sample_ages_s']
        assert ages.shape == (count,)
        assert ages[0] == 0. and ages[-1] == kw['duration_s']
        assert np.all(np.diff(ages) > 0.)
        indices = np.searchsorted(full['sample_ages_s'], ages)
        np.testing.assert_array_equal(ages, full['sample_ages_s'][indices])
        for key in ('states', 'atomic_covariance_by_source'):
            np.testing.assert_array_equal(sampled[key], full[key][indices])
        for key in ('greater_by_source', 'lesser_by_source', 'mean_pulse',
                    'retarded_response', 'exit_state'):
            np.testing.assert_array_equal(sampled[key], full[key])
        # Audit statistics cover all substeps, regardless of retained history.
        assert sampled['audit'] == full['audit']
        assert sampled['numerics']['stored_time_points'] == count
        _assert_quantum_consistency(sampled)


def test_three_level_smooth_path_has_fourth_order_convergence_to_reference():
    kw = _physical_problem(3)
    reference = smooth_wavepacket(**kw, rtol=2e-12, atol=2e-14, max_step_s=kw['duration_s']/96)
    qrt_kw = {key: value for key, value in kw.items()
              if key not in ('drives', 'drive_frequencies_rad_s')}
    qrt = smooth_qrt(**qrt_kw, rtol=2e-12, atol=2e-14, max_step_s=kw['duration_s']/96)
    for key in ('greater', 'lesser', 'mean_pulse', 'exit_state'):
        np.testing.assert_allclose(reference[key], qrt[key], rtol=3e-10, atol=3e-12)
    keys = ('greater_by_source', 'lesser_by_source', 'mean_pulse', 'retarded_response', 'exit_state')
    errors = {key: [] for key in keys}
    for segments in (8, 16, 32):
        result = exponential_wavepacket(**kw, segments=segments, order=4, sample_count=None)
        _assert_quantum_consistency(result)
        assert result['states'].shape == (segments+1, 3, 3)
        assert result['atomic_covariance_by_source'].shape[0] == segments+1
        assert result['source_names'] == reference['source_names']
        np.testing.assert_allclose(result['sample_ages_s'], np.linspace(0., kw['duration_s'], segments+1))
        for key in keys:
            errors[key].append(np.linalg.norm(result[key]-reference[key])/np.linalg.norm(reference[key]))
    for key, values in errors.items():
        ratios = np.asarray(values[:-1])/values[1:]
        assert np.all((ratios > 12.) & (ratios < 20.)), (key, values, ratios)
        assert values[-1] < 2e-7, (key, values)
    midpoint = exponential_wavepacket(**kw, segments=32, order=2)
    for key in keys:
        midpoint_error = np.linalg.norm(midpoint[key]-reference[key])/np.linalg.norm(reference[key])
        assert midpoint_error > 30*errors[key][-1], (key, midpoint_error, errors[key])


@pytest.mark.parametrize('field,value', [
    ('segments', 0), ('segments', -2), ('segments', True), ('segments', 2.5),
    ('order', 0), ('order', 3), ('order', True),
    ('sample_count', 0), ('sample_count', 1), ('sample_count', True),
    ('sample_count', 2.5), ('sample_count', np.nan), ('sample_count', np.inf),
    ('max_condition', .9), ('max_condition', np.inf), ('max_condition', np.nan),
    ('duration_s', 0.), ('duration_s', -1.), ('duration_s', np.inf), ('duration_s', np.nan),
    ('envelope', None), ('envelope', lambda t: [1.]),
    ('envelope', lambda t: 1j), ('envelope', lambda t: np.nan), ('envelope', lambda t: np.inf),
    ('h0', np.array([[0., 1.], [0., 0.]])), ('h1', np.full((2, 2), np.nan)),
    ('h0', np.eye(3)),
    ('boundary_state', np.diag([1.2, -.2])), ('boundary_state', np.eye(2)),
    ('boundary_state', np.array([[.5, .2j], [0., .5]])),
    ('boundary_state', np.eye(3)/3), ('boundary_state', np.full((2, 2), np.nan)),
    ('readouts', np.eye(2)), ('readouts', np.empty((0, 2, 2))),
    ('readouts', np.full((2, 2, 2), np.nan)),
    ('frequencies_rad_s', [[.2]]), ('frequencies_rad_s', [[1j, 0.]]),
    ('frequencies_rad_s', np.empty((0, 2))), ('frequencies_rad_s', [[np.inf, 0.]]),
    ('drives', np.eye(2)), ('drives', np.empty((0, 2, 2))),
    ('drives', None), ('drive_frequencies_rad_s', None),
    ('drive_frequencies_rad_s', [[0.]]), ('drive_frequencies_rad_s', [[1j, 0.], [0., 0.]]),
])
def test_invalid_wavepacket_input_is_rejected(field, value):
    kw = {**_physical_problem(), 'segments': 2, field: value}
    with pytest.raises(ValueError):
        exponential_wavepacket(**kw)


def test_reservoirs_must_be_explicit():
    with pytest.raises(TypeError, match='explicit reservoirs'):
        exponential_wavepacket(**{**_physical_problem(), 'reservoirs': None}, segments=2)
