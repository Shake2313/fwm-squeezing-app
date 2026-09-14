"""Boundary, quantum ordering, causal memory and finite residence controls."""

from dataclasses import replace

import numpy as np
import pytest

from gabes import core
from gabes.quantum.contracts import GeneratorFrequencyAxis
from gabes.quantum.diffusion import ordered_jump_diffusion, stationary_atomic_noise
from gabes.quantum.reservoirs import CollapseChannel, ExplicitReservoirs
from gabes.quantum.transport import BallisticPath, integrate_characteristic, poisson_beam_spectrum
from analysis.grand_challenge.transport_fixtures import moving_fixture
from analysis.grand_challenge.reference.transport_qrt import characteristic_qrt


@pytest.fixture(autouse=True)
def single_thread():
    with core.blas_single_thread():
        yield


@pytest.fixture(scope='module')
def moving():
    with core.blas_single_thread():
        f = moving_fixture()
        a = integrate_characteristic(f['path'], f['hamiltonian'], f['reservoirs'], f['boundary_state'])
        axis = GeneratorFrequencyAxis(2*np.pi*np.array([-.3e6, 0., .3e6]), 'physical Fourier frequency')
        packet = a.wavepacket(axis, f['readout'], drive_operators=f['drives'])
        reference = characteristic_qrt(f['path'], f['hamiltonian'], f['reservoirs'], f['boundary_state'],
            f['readout'], axis.omega_rad_s, drives=f['drives'])
    return f, a, packet, reference


def test_path_declares_the_actual_entry_velocity_and_exit():
    p = BallisticPath([0, 1, 2], [3, 4, 5], 2., 'declared')
    np.testing.assert_array_equal(p.position([0, 1, 2]), [[0,1,2], [3,5,7], [6,9,12]])
    with pytest.raises(ValueError, match='outside'):
        p.position(2.01)
    for kwargs in ({'residence_time_s': 0}, {'velocity_m_s': [1, 2]}, {'source': ''}):
        with pytest.raises(ValueError):
            replace(p, **kwargs)


def test_boundary_plus_jumps_reproduce_the_actual_nonstationary_atomic_moments(moving):
    f, a, _, _ = moving
    assert a.audit['passed']
    assert a.audit['covariance_relative_residual'] < 1e-10
    for index in (0, 8, 16):
        _, _, cached = a.coefficients(a.times_s[index], a.states[index])
        direct = ordered_jump_diffusion(f['reservoirs'], a.states[index], a.operators)
        np.testing.assert_allclose(cached, direct, rtol=1e-11, atol=1e-8)
    assert np.linalg.norm(a.covariance_by_source[0, 4]) > .1


@pytest.mark.parametrize('key', ['greater', 'lesser', 'mean_pulse', 'retarded_response'])
def test_full_density_time_qrt_and_drive_reference(moving, key):
    _, _, actual, ref = moving
    assert np.linalg.norm(actual[key]-ref[key])/np.linalg.norm(ref[key]) < 1e-9


def test_disjoint_positions_retain_cross_covariance_and_global_positivity(moving):
    f, a, _, _ = moving
    k = a.two_time_covariance(f['readout'])
    full = k.transpose(0,2,1,3).reshape(34,34)
    assert np.linalg.eigvalsh(full)[0] > -1e-10
    assert np.linalg.norm(k[9, 6]) > .05
    pieces = sum(a.two_time_covariance(f['readout'], source=i) for i in range(len(a.source_names)))
    np.testing.assert_allclose(k, pieces, rtol=1e-11, atol=1e-13)
    np.testing.assert_allclose(a.transition(16, 0), a.transition(16, 8)@a.transition(8, 0), rtol=1e-11, atol=1e-14)
    with pytest.raises(ValueError):
        a.transition(0, 2)


def test_retarded_kernel_never_responds_before_the_source(moving):
    f, a, _, _ = moving
    k = a.retarded_kernel(f['readout'], f['drives'])
    assert np.all(k[np.triu_indices(len(a.times_s), 1)] == 0)
    assert np.linalg.norm(k[9, 6]) > 1e3


def test_discarding_inflow_can_leave_psd_but_break_the_atomic_commutator(moving):
    _, a, packet, _ = moving
    correct = a.covariance_by_source.sum(axis=0)
    erased = a.covariance_by_source[1:].sum(axis=0)
    assert np.linalg.eigvalsh(erased).min() >= -1e-10
    defect = erased-erased.swapaxes(-1,-2)-correct+correct.swapaxes(-1,-2)
    assert np.linalg.norm(defect) > .1
    erased_spectrum = packet['greater_by_source'][1:].sum(axis=0)
    assert np.linalg.eigvalsh(erased_spectrum).min() > 0
    assert np.linalg.norm(packet['greater']-erased_spectrum)/np.linalg.norm(packet['greater']) > .05


def test_poisson_number_noise_is_fixed_by_the_mean_pulse(moving):
    _, _, p, _ = moving
    a = poisson_beam_spectrum(p, 2e6, source='declared independent Poisson entries')
    expected = 2e6*p['mean_pulse'][:,:,None]*p['mean_pulse'][:,None,:].conj()
    np.testing.assert_array_equal(a['poisson_number'], expected)
    np.testing.assert_allclose(a['greater']-a['lesser'], 2e6*(p['greater']-p['lesser']), rtol=1e-12, atol=1e-20)
    assert a['mean_occupancy'] == 8
    with pytest.raises(ValueError):
        poisson_beam_spectrum(p, -1., source='invalid')


def test_identity_readout_has_no_internal_noise_but_has_poisson_count_noise():
    path = BallisticPath([0,0,0],[1,0,0],2.,'identity-count test')
    a = integrate_characteristic(path, lambda t,r:np.zeros((2,2)), ExplicitReservoirs(2,()), np.diag([.6,.4]))
    axis = GeneratorFrequencyAxis([0., .7], 'finite residence')
    p = a.wavepacket(axis, lambda t,r:np.eye(2)[None])
    np.testing.assert_allclose(p['greater'], 0., atol=1e-25)
    expected = np.array([2., np.expm1(1.4j)/(.7j)])
    np.testing.assert_allclose(p['mean_pulse'][:,0], expected, rtol=1e-11)
    psd = poisson_beam_spectrum(p, 3., source='Poisson counting control')
    np.testing.assert_allclose(psd['greater'][:,0,0],3*abs(expected)**2,rtol=1e-11)


def test_finite_residence_matches_analytic_lorentzian_window_and_stationary_limit():
    lower = np.array([[0.,1.],[0.,0.]])
    res = ExplicitReservoirs(2,[CollapseChannel('emission',lower,'decay','unit toy decay')])
    h = np.diag([0., .7]); rho = np.diag([1., 0.])
    axis = GeneratorFrequencyAxis([.2,1.1], 'finite window')
    a = integrate_characteristic(BallisticPath([0,0,0],[0,0,1],3.,'toy'),lambda t,r:h,res,rho)
    p = a.wavepacket(axis,lambda t,r:lower[None])
    s = -.5+1j*(axis.omega_rad_s-.7)
    exact = 2*((np.expm1(3*s)-3*s)/s**2).real
    np.testing.assert_allclose(p['greater'][:,0,0],exact,rtol=1e-10)
    np.testing.assert_allclose(p['lesser'],0.,atol=1e-13)
    atom = stationary_atomic_noise(h,res)
    read = np.einsum('ab,iba->i',lower,atom.operators)
    stationary = np.real(read@atom.spectrum(axis).ordered@read.conj())
    errors=[]
    for duration in (20.,100.,500.):
        finite = 2*((np.expm1(duration*s)-duration*s)/s**2).real/duration
        errors.append(np.linalg.norm(finite-stationary))
    assert errors[1] < .21*errors[0]
    assert errors[2] < .21*errors[1]


def test_unphysical_boundary_and_missing_endpoints_are_rejected():
    f=moving_fixture()
    for state,times in ((np.diag([1.1,-.1]),None),(f['boundary_state'],[0.,1e-6])):
        with pytest.raises(ValueError):
            integrate_characteristic(f['path'],f['hamiltonian'],f['reservoirs'],state,times_s=times)


def test_seconds_and_rescaled_time_are_exactly_the_same_wavepacket(moving):
    f,_,packet,_=moving
    duration=f['path'].residence_time_s
    path=replace(f['path'],velocity_m_s=f['path'].velocity_m_s*duration,residence_time_s=1.)
    res=ExplicitReservoirs(2,tuple(replace(c,operator=np.sqrt(duration)*c.operator) for c in f['reservoirs'].channels))
    a=integrate_characteristic(path,lambda u,r:duration*f['hamiltonian'](u*duration,r),res,f['boundary_state'])
    p=a.wavepacket(GeneratorFrequencyAxis(packet['frequency_axis'].omega_rad_s*duration,'dimensionless time'),
        lambda u,r:f['readout'](u*duration,r),drive_operators=lambda u,r:duration*f['drives'](u*duration,r))
    for key,power in (('greater',2),('lesser',2),('mean_pulse',1),('retarded_response',1)):
        assert np.linalg.norm(p[key]*duration**power-packet[key])/np.linalg.norm(packet[key])<1e-9


def test_three_level_transient_does_not_require_a_stationary_mixing_state():
    path=BallisticPath([-.4,0,0],[.5,0,.2],1.6,'synthetic finite Lambda atom')
    l0=np.zeros((3,3),complex);l0[0,2]=np.sqrt(.7)
    l1=np.zeros((3,3),complex);l1[1,2]=np.sqrt(.3)
    res=ExplicitReservoirs(3,[CollapseChannel('to0',l0,'decay','toy'),CollapseChannel('to1',l1,'decay','toy')])
    def h(t,r):
        out=np.diag([0.,.2,.6]).astype(complex)
        out[0,2]=.5*np.exp(-r[0]**2)*np.exp(.7j*t)
        out[1,2]=.3*np.exp(-.5j*t)
        out[2,0],out[2,1]=out[0,2].conjugate(),out[1,2].conjugate()
        return out
    read=lambda t,r:np.array([l0+.2*l1.conj().T,l1+.1*np.eye(3)])
    rho=np.diag([.6,.3,.1]);axis=GeneratorFrequencyAxis([-.3,.4],'toy finite-time Fourier')
    a=integrate_characteristic(path,h,res,rho)
    p=a.wavepacket(axis,read)
    ref=characteristic_qrt(path,h,res,rho,read,axis.omega_rad_s)
    assert a.audit['passed']
    for key in ('greater','lesser','mean_pulse'):
        assert np.linalg.norm(p[key]-ref[key])/np.linalg.norm(ref[key])<1e-9


def test_transport_audit_will_not_overwrite_evidence(tmp_path,monkeypatch):
    from analysis.grand_challenge import transport_audit
    p=tmp_path/'kept.json';p.write_text('preserve',encoding='utf-8')
    def forbidden():
        raise AssertionError('must reject before calculating')
    monkeypatch.setattr(transport_audit,'build_report',forbidden)
    with pytest.raises(FileExistsError):
        transport_audit.main(['--output',str(p)])
    assert p.read_text(encoding='utf-8')=='preserve'
