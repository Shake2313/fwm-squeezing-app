"""Continuous microscopic moments versus independent finite-density QRT."""

from dataclasses import replace
import numpy as np
import pytest

from gabes import core
from gabes.quantum.contracts import GeneratorFrequencyAxis,AnalysisFrequencyAxis
from gabes.quantum.reservoirs import ExplicitReservoirs,CollapseChannel
from gabes.quantum.transport import BallisticPath,integrate_characteristic
from gabes.quantum.smooth_transport import smooth_wavepacket
from gabes.quantum.segmented_transport import segmented_wavepacket
from gabes.fwm_quantum.smooth_transport import smooth_rb_problem,smooth_rb_wavepacket
from gabes.fwm_quantum.kinetic import CarrierGeometry
from analysis.grand_challenge.normalization_audit import conditional_inputs
from analysis.grand_challenge.reference.transport_qrt import characteristic_qrt


@pytest.fixture(autouse=True)
def single_thread():
    with core.blas_single_thread(): yield


def toy(n=2):
    rng=np.random.default_rng(371+n)
    hs=rng.normal(size=(2,n,n))+1j*rng.normal(size=(2,n,n))
    hs=(hs+hs.conj().swapaxes(-1,-2))/3
    lower=np.diag(np.ones(n-1),1)
    r=ExplicitReservoirs(n,(CollapseChannel('cascade',lower,'decay','synthetic'),))
    ops=np.array([lower+.15*np.eye(n),lower.T+.2j*np.diag(np.arange(n))])
    return hs,r,np.eye(n)/n,ops


@pytest.mark.parametrize('n',[2,3])
def test_smooth_qrts_source_covariance_and_complex_response(n):
    hs,r,rho,ops=toy(n)
    t=1.2
    envelope=lambda age:np.exp(-((age-.4)/.7)**2)
    offsets=np.array([.9,-.4])
    axis=GeneratorFrequencyAxis([-.3,.5],'declared toy physical frequency')
    freq=axis.omega_rad_s[:,None]+offsets
    drives=ops.conj().swapaxes(-1,-2)
    p=smooth_wavepacket(*hs,envelope,r,rho,t,ops,freq,drives=drives,
        drive_frequencies_rad_s=freq,rtol=1e-11,atol=1e-13)
    path=BallisticPath([0,0,0],[0,0,0],t,'finite fixed-center control')
    h=lambda a,pos:hs[0]+envelope(a)*hs[1]
    read=lambda a,pos:ops*np.exp(1j*offsets*a)[:,None,None]
    drive=lambda a,pos:drives*np.exp(-1j*offsets*a)[:,None,None]
    ref=characteristic_qrt(path,h,r,rho,read,axis.omega_rad_s,drives=drive)
    for key in ('greater','lesser','mean_pulse','retarded_response','exit_state'):
        np.testing.assert_allclose(p[key],ref[key],rtol=2e-10,atol=2e-12)
    old=integrate_characteristic(path,h,r,rho,rtol=1e-11,atol=1e-13)
    q=old.wavepacket(axis,read)
    for key in ('greater_by_source','lesser_by_source'):
        np.testing.assert_allclose(p[key],q[key],rtol=2e-10,atol=2e-12)
    assert p['audit']['passed']


def test_constant_map_seconds_units_and_shared_rf_density():
    hs,r,rho,ops=toy()
    w=np.array([[.2,-.4],[.7,.3]])
    h=hs[0]+.4*hs[1]
    p=smooth_wavepacket(*hs,lambda t:.4,r,rho,1.,ops,w,rtol=1e-11,atol=1e-13)
    q=segmented_wavepacket([h],r,rho,[1.],[ops],w)
    factor=1e-6
    rr=ExplicitReservoirs(2,tuple(replace(ch,operator=ch.operator/np.sqrt(factor)) for ch in r.channels))
    si=smooth_wavepacket(hs[0]/factor,hs[1]/factor,lambda t:.4,rr,rho,factor,ops,w/factor,rtol=1e-11,atol=1e-13)
    for key in ('greater','lesser','greater_by_source','lesser_by_source','mean_pulse','exit_state'):
        np.testing.assert_allclose(p[key],q[key],rtol=2e-11,atol=2e-13)
        power=0 if key=='exit_state' else 1 if key=='mean_pulse' else 2
        np.testing.assert_allclose(si[key]/factor**power,p[key],rtol=2e-11,atol=2e-13)


def test_source_atomic_transpose_identity_and_frequency_batch_invariance():
    hs,r,rho,ops=toy()
    kwargs=dict(h0=hs[0],h1=hs[1],envelope=lambda t:1-t/2,reservoirs=r,
        boundary_state=rho,duration_s=1.,readouts=ops,rtol=1e-11,atol=1e-13)
    p=smooth_wavepacket(**kwargs,frequencies_rad_s=[[.2,-.3],[.7,-.6]])
    c=p['atomic_covariance_by_source']
    np.testing.assert_array_equal(c[:,1],c[:,0].swapaxes(-1,-2))
    for index,w in enumerate(([.2,-.3],[.7,-.6])):
        q=smooth_wavepacket(**kwargs,frequencies_rad_s=[w])
        for key in ('greater','lesser','mean_pulse'):
            np.testing.assert_allclose(p[key][index],q[key][0],rtol=2e-10,atol=2e-12)


def test_identity_and_closed_system_do_not_need_stationary_state():
    hs,_,rho,_=toy()
    p=smooth_wavepacket(*hs,lambda t:np.sin(t),ExplicitReservoirs(2,()),rho,1.,
        np.array([np.eye(2)]),[[.7]],rtol=1e-11,atol=1e-13)
    np.testing.assert_allclose(p['greater'],0.,atol=1e-14)
    np.testing.assert_allclose(p['lesser'],0.,atol=1e-14)
    np.testing.assert_allclose(p['mean_pulse'][0,0],np.expm1(.7j)/(.7j),rtol=1e-12)


@pytest.mark.parametrize('change',['nonhermitian','duration','frequency','envelope_complex','envelope_nan','drive','rtol','atol','samples','step'])
def test_contracts(change):
    hs,r,rho,ops=toy()
    kw=dict(h0=hs[0],h1=hs[1],envelope=lambda t:.2,reservoirs=r,
        boundary_state=rho,duration_s=1.,readouts=ops,frequencies_rad_s=[[.2,-.3]])
    if change=='nonhermitian':kw['h1']=np.array([[0,1],[0,0]])
    if change=='duration':kw['duration_s']=0
    if change=='frequency':kw['frequencies_rad_s']=[[.2]]
    if change=='envelope_complex':kw['envelope']=lambda t:1j
    if change=='envelope_nan':kw['envelope']=lambda t:np.nan
    if change=='drive':kw['drive_frequencies_rad_s']=[[0.]]
    if change=='rtol':kw['rtol']=0.
    if change=='atol':kw['atol']=-1.
    if change=='samples':kw['sample_count']=True
    if change=='step':kw['max_step_s']=0.
    with pytest.raises(ValueError): smooth_wavepacket(**kw)


def test_rb_continuous_profile_matches_frozen_midpoint_only_at_midpoint():
    i=replace(conditional_inputs(),transit_rate_s_inverse=0.)
    g=CarrierGeometry.vacuum_beams(i,probe_angle_rad=.006,conjugate_angle_rad=-.005)
    path=BallisticPath([-1e-4,2e-5,0],[150,10,100],2e-6,'declared moving control')
    p=smooth_rb_problem(i,g,path,AnalysisFrequencyAxis.from_hz([1e6]))
    for t in (0.,.25e-6,1e-6,2e-6):
        r=path.position(t)
        np.testing.assert_allclose(p['envelope'](t),np.exp(-sum(r[:2]**2)/i.pump_waist_m**2),rtol=2e-14)
    assert abs(p['envelope'](0)-p['envelope'](2e-6))>.05
    assert len(p['reservoirs'].channels)==4
    assert abs(p['metadata']['loop_convective_frequency_rad_s'])>1e5


def test_short_physical_rb_phase_theorem_survives_smooth_envelope():
    i=replace(conditional_inputs(),transit_rate_s_inverse=0.)
    g=CarrierGeometry.vacuum_beams(i,probe_angle_rad=.006,conjugate_angle_rad=-.005)
    path=BallisticPath([-1e-4,2e-5,0],[150,10,100],2e-9,'short validation path')
    axis=AnalysisFrequencyAxis.from_hz([1e6])
    p=smooth_rb_wavepacket(i,g,path,axis,rtol=1e-10,atol=1e-13)
    phi=.73
    q=smooth_rb_wavepacket(i,g,path,axis,entry_phase_rad=phi,rtol=1e-10,atol=1e-13)
    u=np.exp(1j*phi*np.array([1,-1,-1,1]))
    for key in ('greater','lesser','greater_by_source','lesser_by_source','retarded_response'):
        np.testing.assert_allclose(q[key],p[key]*u[:,None]*u.conj()[None,:],rtol=2e-8,atol=1e-28)


def test_smooth_audit_refuses_overwrite_before_worker_creation(tmp_path):
    from analysis.grand_challenge.smooth_transport_audit import write_audit
    p=tmp_path/'keep.json'; p.write_text('preserve',encoding='utf-8')
    with pytest.raises(FileExistsError):write_audit(p,tmp_path/'figure.png')
    assert p.read_text(encoding='utf-8')=='preserve'
