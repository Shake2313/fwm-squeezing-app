"""Backward source attribution versus forward moments and physical drives."""

from dataclasses import replace

import numpy as np
import pytest
from scipy.integrate import quad, solve_ivp

from gabes import core
from gabes.quantum.reservoirs import ExplicitReservoirs, CollapseChannel
from gabes.quantum.segmented_transport import segmented_wavepacket
from gabes.quantum.smooth_transport import smooth_wavepacket
from analysis.grand_challenge.reference.adjoint_transport import adjoint_wavepacket
from analysis.grand_challenge.reference.smooth_qrt import smooth_qrt


@pytest.fixture(autouse=True)
def single_thread():
    with core.blas_single_thread():
        yield


def fixture(n=2):
    rng = np.random.default_rng(831+n)
    h = rng.normal(size=(2, n, n))+1j*rng.normal(size=(2, n, n))
    h = (h+h.conj().swapaxes(-1, -2))/3
    x = rng.normal(size=(n, n))+1j*rng.normal(size=(n, n))
    rho = x@x.conj().T
    rho /= rho.trace()
    lower = np.diag(np.ones(n-1), 1)
    reservoirs = ExplicitReservoirs(n, (
        CollapseChannel('decay', lower, 'decay', 'test fixture'),
        CollapseChannel('pump', .3*lower.T, 'excitation', 'test fixture'),
        CollapseChannel('dephasing', .2*np.diag(np.arange(n)), 'dephasing', 'test fixture')))
    ops = np.array([lower+.17j*np.eye(n), x/3])
    return dict(h0=h[0], h1=h[1], envelope=lambda a: np.exp(-((a-.3)/.7)**2),
                reservoirs=reservoirs, boundary_state=rho, duration_s=1.1,
                readouts=ops, frequencies_rad_s=np.array([[.2, -.7], [.5, .3]]),
                drives=ops.conj().swapaxes(-1, -2),
                drive_frequencies_rad_s=np.array([[.6, -.4], [.1, .9]]), rtol=1e-10, atol=1e-12)


@pytest.mark.parametrize('n', [2, 3])
def test_all_named_sources_and_unequal_complex_response(n):
    kw = fixture(n)
    adj = adjoint_wavepacket(**kw)
    forward = smooth_wavepacket(**kw)
    qrt = smooth_qrt(**{k: v for k, v in kw.items() if k not in ('drives', 'drive_frequencies_rad_s')})
    for key in ('greater', 'lesser', 'greater_by_source', 'lesser_by_source',
                'mean_pulse', 'retarded_response', 'exit_state'):
        np.testing.assert_allclose(adj[key], forward[key], rtol=3e-9, atol=2e-11)
    for key in ('greater', 'lesser', 'mean_pulse', 'exit_state'):
        np.testing.assert_allclose(adj[key], qrt[key], rtol=3e-9, atol=2e-11)
    for key in ('greater_by_source', 'lesser_by_source'):
        assert np.linalg.eigvalsh(adj[key]).min() > -1e-11
    assert adj['source_names'] == forward['source_names']


def test_named_jump_splitting_phase_and_permutation():
    kw = fixture()
    channel = kw['reservoirs'].channels[0]
    kw['reservoirs'] = ExplicitReservoirs(2, (channel,))
    base = adjoint_wavepacket(**kw)
    split = (replace(channel, name='a', operator=.6*np.exp(.37j)*channel.operator),
             replace(channel, name='b', operator=.8*np.exp(-.29j)*channel.operator))
    divided = adjoint_wavepacket(**{**kw, 'reservoirs': ExplicitReservoirs(2, split[::-1])})
    for key in ('greater', 'lesser', 'mean_pulse', 'retarded_response', 'exit_state'):
        np.testing.assert_allclose(divided[key], base[key], rtol=2e-10, atol=2e-12)
    for key in ('greater_by_source', 'lesser_by_source'):
        np.testing.assert_allclose(divided[key][0], base[key][0], rtol=2e-10, atol=2e-12)
        np.testing.assert_allclose(divided[key][1], .64*base[key][1], rtol=2e-10, atol=2e-12)
        np.testing.assert_allclose(divided[key][2], .36*base[key][1], rtol=2e-10, atol=2e-12)
    assert divided['source_names'] == ('atomic_inflow', 'jump:b', 'jump:a')


def test_exact_constant_segments_and_seconds_units():
    kw = fixture()
    kw['envelope'] = lambda a: .4
    p = adjoint_wavepacket(**kw)
    ref = segmented_wavepacket([kw['h0']+.4*kw['h1']], kw['reservoirs'], kw['boundary_state'],
                               [kw['duration_s']], [kw['readouts']], kw['frequencies_rad_s'],
                               drives=[kw['drives']], drive_frequencies_rad_s=kw['drive_frequencies_rad_s'])
    factor = 1e-6
    scaled = {**kw, 'h0': kw['h0']/factor, 'h1': kw['h1']/factor,
              'duration_s': kw['duration_s']*factor,
              'frequencies_rad_s': kw['frequencies_rad_s']/factor,
              'drive_frequencies_rad_s': kw['drive_frequencies_rad_s']/factor,
              'drives': kw['drives']/np.sqrt(factor),
              'reservoirs': ExplicitReservoirs(2, tuple(replace(c, operator=c.operator/np.sqrt(factor))
                                                      for c in kw['reservoirs'].channels))}
    si = adjoint_wavepacket(**scaled)
    for key in ('greater', 'lesser', 'greater_by_source', 'lesser_by_source',
                'mean_pulse', 'retarded_response', 'exit_state'):
        power = 0 if key == 'exit_state' else 1 if key == 'mean_pulse' else 1.5 if key == 'retarded_response' else 2
        np.testing.assert_allclose(p[key], ref[key], rtol=3e-10, atol=3e-12)
        np.testing.assert_allclose(si[key]/factor**power, p[key], rtol=3e-10, atol=3e-12)


def test_rf_batching_and_identity_invariance():
    kw = fixture()
    batch = adjoint_wavepacket(**kw)
    for j in range(2):
        single = adjoint_wavepacket(**{**kw, 'frequencies_rad_s': kw['frequencies_rad_s'][j:j+1],
                                       'drive_frequencies_rad_s': kw['drive_frequencies_rad_s'][j:j+1]})
        for key in ('greater', 'lesser', 'mean_pulse', 'retarded_response'):
            np.testing.assert_allclose(single[key][0], batch[key][j], rtol=3e-9, atol=2e-11)
    identity = adjoint_wavepacket(**{**kw, 'readouts': np.array([np.eye(2), np.eye(2)]),
                                     'frequencies_rad_s': np.zeros((2, 2)), 'drives': None,
                                     'drive_frequencies_rad_s': None})
    np.testing.assert_allclose(identity['mean_pulse'], kw['duration_s'], atol=2e-12)
    np.testing.assert_allclose(identity['greater_by_source'], 0., atol=2e-12)
    assert identity['retarded_response'].shape == (2, 2, 0)


def test_no_jump_closed_atom_has_only_boundary_noise():
    kw = fixture()
    kw['reservoirs'] = ExplicitReservoirs(2, ())
    adj = adjoint_wavepacket(**kw)
    ref = smooth_wavepacket(**kw)
    assert adj['greater_by_source'].shape[0] == 1
    for key in ('greater', 'lesser', 'mean_pulse', 'retarded_response'):
        np.testing.assert_allclose(adj[key], ref[key], rtol=2e-9, atol=2e-11)


def test_analytic_ground_state_boundary_and_vacuum_jump_noise():
    gamma, resonance, frequency, duration = .8, 1.3, .6, 1.7
    lower = np.array([[0., 1.], [0., 0.]])
    reservoirs = ExplicitReservoirs(2, (CollapseChannel('radiation', np.sqrt(gamma)*lower, 'decay', 'analytic control'),))
    a = adjoint_wavepacket(np.diag([0., resonance]), np.zeros((2, 2)), lambda t: 0., reservoirs,
                          np.diag([1., 0.]), duration, lower[None], [[frequency]], rtol=1e-11, atol=1e-13)
    z = -gamma/2+1j*(frequency-resonance)
    pulse = lambda x: np.expm1(z*x)/z
    boundary = abs(pulse(duration))**2
    bath = gamma*quad(lambda x: abs(pulse(x))**2, 0., duration, epsabs=1e-13)[0]
    total = 2*np.real((np.expm1(z*duration)-z*duration)/z**2)
    np.testing.assert_allclose(a['greater_by_source'][:, 0, 0, 0], [boundary, bath], rtol=2e-10)
    np.testing.assert_allclose(a['greater'], total, rtol=2e-10)
    np.testing.assert_allclose(a['lesser_by_source'], 0., atol=1e-13)
    np.testing.assert_allclose(a['mean_pulse'], 0., atol=1e-13)


def test_analytic_complex_response_sign_and_mixed_frequencies():
    lower = np.array([[0., 1.], [0., 0.]])
    duration = 1.3
    w, v = np.array([.4, -.3, 0.]), np.array([-.7, .9, 0.])
    a = adjoint_wavepacket(np.zeros((2, 2)), np.zeros((2, 2)), lambda t: 0., ExplicitReservoirs(2, ()),
                          np.diag([1., 0.]), duration, lower[None], w[:, None],
                          drives=lower.T[None], drive_frequencies_rad_s=v[:, None])
    phi = lambda x: duration if x == 0 else np.expm1(1j*x*duration)/(1j*x)
    exact = [-1j*(phi(x)-phi(x-y))/(1j*y) if y else -1j*duration**2/2 for x, y in zip(w, v)]
    np.testing.assert_allclose(a['retarded_response'][:, 0, 0], exact, rtol=2e-11, atol=2e-13)


def test_identity_readout_shift_changes_only_mean():
    kw = fixture()
    c = np.array([.7+.3j, -.4j])
    base = adjoint_wavepacket(**kw)
    shifted = adjoint_wavepacket(**{**kw, 'readouts': kw['readouts']+c[:, None, None]*np.eye(2)})
    for key in ('greater_by_source', 'lesser_by_source', 'retarded_response'):
        np.testing.assert_allclose(shifted[key], base[key], rtol=2e-9, atol=2e-11)
    w = kw['frequencies_rad_s']
    np.testing.assert_allclose(shifted['mean_pulse']-base['mean_pulse'],
                               c*np.expm1(1j*w*kw['duration_s'])/(1j*w), rtol=2e-11, atol=2e-13)


def test_retarded_response_matches_physical_cosine_drive_difference():
    kw = fixture()
    v = np.array([[.2, .4j], [-.4j, -.1]])
    frequency = .8
    kw.update(drives=np.array([v, v]), drive_frequencies_rad_s=np.tile([frequency, -frequency], (2, 1)))
    reference = adjoint_wavepacket(**kw)
    n = 2
    def pulse(epsilon):
        def rhs(t, y):
            rho = y[:4].reshape(n, n)
            h = kw['h0']+kw['envelope'](t)*kw['h1']+epsilon*np.cos(frequency*t)*v
            rate = kw['reservoirs'].generator(h)@y[:4]
            observed = np.einsum('jab,ba->j', kw['readouts'], rho)*np.exp(1j*kw['frequencies_rad_s']*t)
            return np.r_[rate, observed.ravel()]
        ode = solve_ivp(rhs, (0., kw['duration_s']), np.r_[kw['boundary_state'].ravel(), np.zeros(4, complex)],
                        method='DOP853', rtol=1e-12, atol=1e-14)
        return ode.y[4:, -1].reshape(2, 2)
    epsilon = 1e-4
    difference = (pulse(epsilon)-pulse(-epsilon))/(2*epsilon)
    np.testing.assert_allclose(difference, reference['retarded_response'].mean(axis=-1), rtol=3e-7, atol=2e-9)


@pytest.mark.parametrize('change', ['hamiltonian', 'state', 'duration', 'readouts', 'frequency', 'empty_frequency',
                                  'drives', 'drive_frequency', 'orphan_frequency', 'rtol', 'atol', 'step', 'envelope', 'nan_envelope'])
def test_invalid_contracts(change):
    kw = fixture()
    if change == 'hamiltonian': kw['h1'] = np.array([[0, 1], [0, 0]])
    if change == 'state': kw['boundary_state'] = np.diag([1.2, -.2])
    if change == 'duration': kw['duration_s'] = 0.
    if change == 'readouts': kw['readouts'] = np.eye(2)
    if change == 'frequency': kw['frequencies_rad_s'] = [[1j, 0]]
    if change == 'empty_frequency': kw['frequencies_rad_s'] = np.empty((0, 2))
    if change == 'drives': kw['drives'] = np.eye(2)
    if change == 'drive_frequency': kw['drive_frequencies_rad_s'] = [[0]]
    if change == 'orphan_frequency': kw['drives'] = None
    if change == 'rtol': kw['rtol'] = 1.
    if change == 'atol': kw['atol'] = 0.
    if change == 'step': kw['max_step_s'] = -1.
    if change == 'envelope': kw['envelope'] = lambda a: [1.]
    if change == 'nan_envelope': kw['envelope'] = lambda a: np.nan
    with pytest.raises(ValueError):
        adjoint_wavepacket(**kw)
