"""Fast carrier integration must preserve inflow, source noise and response."""

from dataclasses import replace
import numpy as np
import pytest

from gabes import core
from gabes.quantum.contracts import GeneratorFrequencyAxis, AnalysisFrequencyAxis
from gabes.quantum.reservoirs import CollapseChannel, ExplicitReservoirs
from gabes.quantum.segmented_transport import segmented_wavepacket
from gabes.quantum.transport import BallisticPath, integrate_characteristic
from gabes.fwm_quantum.kinetic import CarrierGeometry
from gabes.fwm_quantum.transport import reduced_ballistic_problem, reduced_ballistic_wavepacket, average_pump_only_entry_phase
from analysis.grand_challenge.normalization_audit import conditional_inputs
from analysis.grand_challenge.reference.transport_qrt import characteristic_qrt
from analysis.grand_challenge.reference.segmented_qrt import segmented_qrt


@pytest.fixture(autouse=True)
def one_thread():
    with core.blas_single_thread():
        yield


def fixture():
    h = np.array([[0., .31+.1j], [.31-.1j, .7]])
    lower = np.array([[0, 1], [0, 0]])
    r = ExplicitReservoirs(2, (CollapseChannel('decay', lower, 'radiative', 'synthetic'),))
    rho = np.array([[.6, .1j], [-.1j, .4]])
    ops = np.array([lower+.1*np.eye(2), lower.conj().T+.2*np.diag([1,-1])])
    return h, r, rho, ops


@pytest.mark.parametrize('offsets', [[0., 0.], [.9, -1.3]])
def test_independent_density_qrt_and_driven_response(offsets):
    h, r, rho, ops = fixture()
    t = 1.1
    axis = GeneratorFrequencyAxis([-.2, .4], 'synthetic physical frequency')
    w = axis.omega_rad_s[:, None]+offsets
    v = ops.conj().swapaxes(-1,-2)
    result = segmented_wavepacket([h], r, rho, [t], [ops], w, drives=[v], drive_frequencies_rad_s=w)
    path = BallisticPath([0,0,0], [0,0,0], t, 'fixed center finite observation')
    read = lambda a, pos: ops*np.exp(1j*np.array(offsets)*a)[:,None,None]
    drive = lambda a, pos: v*np.exp(-1j*np.array(offsets)*a)[:,None,None]
    ref = characteristic_qrt(path, lambda a,pos:h, r, rho, read, axis.omega_rad_s, drives=drive)
    for key in ('greater', 'lesser', 'mean_pulse', 'retarded_response'):
        np.testing.assert_allclose(result[key], ref[key], rtol=2e-10, atol=2e-12)
    old = integrate_characteristic(path, lambda a,pos:h, r, rho)
    packet = old.wavepacket(axis, read, drive_operators=drive)
    for key in ('greater_by_source', 'lesser_by_source'):
        np.testing.assert_allclose(result[key], packet[key], rtol=2e-10, atol=2e-12)


def test_dark_ordering_keeps_signed_roundoff_and_analytic_finite_residence():
    _, r, _, ops = fixture()
    t, w = 3., .2
    p = segmented_wavepacket([np.diag([0,.7])], r, np.diag([1.,0.]), [t], [[ops[0]-.1*np.eye(2)]], [[w]])
    s = -.5+1j*(w-.7)
    expected = 2*((np.exp(s*t)-1-s*t)/s**2).real
    np.testing.assert_allclose(p['greater'][0,0,0], expected, rtol=2e-13)
    assert abs(p['lesser'][0,0,0]) < 2e-14
    assert p['audit']['passed']


def test_exact_segment_composition_and_seconds_rescaling():
    h,r,rho,ops = fixture()
    one = segmented_wavepacket([h], r, rho, [1.2], [ops], [[.3,-.5]])
    many = segmented_wavepacket([h]*7, r, rho, [1.2/7]*7, [ops]*7, [[.3,-.5]])
    assert many['matrix_exponentials'] == 2
    assert many['identical_exponentials_reused'] == 12
    factor = 1e-6
    rr = ExplicitReservoirs(2, tuple(replace(c, operator=c.operator/np.sqrt(factor)) for c in r.channels))
    scaled = segmented_wavepacket([h/factor], rr, rho, [1.2*factor], [ops], np.array([[.3,-.5]])/factor)
    for key in ('greater', 'lesser', 'greater_by_source', 'lesser_by_source', 'mean_pulse', 'exit_state'):
        np.testing.assert_allclose(many[key], one[key], rtol=2e-12, atol=2e-14)
        power = 0 if key=='exit_state' else 1 if key=='mean_pulse' else 2
        np.testing.assert_allclose(scaled[key]/factor**power, one[key], rtol=2e-12, atol=2e-14)


def test_identity_has_no_internal_noise_and_no_jumps_is_valid():
    h,_,rho,_ = fixture()
    p = segmented_wavepacket([h], ExplicitReservoirs(2, ()), rho, [2.], [[np.eye(2)]], [[.7]])
    np.testing.assert_allclose(p['greater'], 0., atol=1e-14)
    np.testing.assert_allclose(p['lesser'], 0., atol=1e-14)
    np.testing.assert_allclose(p['mean_pulse'][0,0], np.expm1(1.4j)/(.7j), rtol=1e-13)


@pytest.mark.parametrize('n', [2, 3])
def test_independent_raw_qrt_after_hamiltonian_and_readout_changes(n):
    rng = np.random.default_rng(483+n)
    hs = rng.normal(size=(3,n,n))+1j*rng.normal(size=(3,n,n))
    hs = (hs+hs.conj().swapaxes(-1,-2))/2
    ops = rng.normal(size=(3,2,n,n))+1j*rng.normal(size=(3,2,n,n))
    jump = np.diag(np.ones(n-1), 1)
    reservoirs = ExplicitReservoirs(n,(CollapseChannel('cascade', jump, 'decay', 'synthetic'),))
    rho = np.eye(n)/n
    args = dict(hamiltonians=hs, reservoirs=reservoirs, boundary_state=rho,
                durations_s=[.31,.29,.4], readouts=ops, frequencies_rad_s=[[.8,-.3],[.1,.2]])
    p, q = segmented_wavepacket(**args), segmented_qrt(**args)
    for key in ('greater','lesser','mean_pulse','exit_state'):
        np.testing.assert_allclose(p[key], q[key], rtol=1e-11, atol=1e-13)


@pytest.mark.parametrize('change', ['duration', 'nonhermitian', 'readout', 'frequency', 'state', 'drive'])
def test_bad_contracts_rejected(change):
    h,r,rho,ops = fixture()
    args = dict(hamiltonians=[h], reservoirs=r, boundary_state=rho, durations_s=[1.], readouts=[ops], frequencies_rad_s=[[.2,.3]])
    if change=='duration': args['durations_s']=[0.]
    if change=='nonhermitian': args['hamiltonians']=[h+np.array([[0,1],[0,0]])]
    if change=='readout': args['readouts']=[ops[0]]
    if change=='frequency': args['frequencies_rad_s']=[[.2]]
    if change=='state': args['boundary_state']=np.diag([2.,-1.])
    if change=='drive': args['drive_frequencies_rad_s']=[[.2]]
    with pytest.raises(ValueError): segmented_wavepacket(**args)


def rb_fixture():
    inputs = replace(conditional_inputs(), transit_rate_s_inverse=0.)
    geometry = CarrierGeometry.vacuum_beams(inputs, probe_angle_rad=.006, conjugate_angle_rad=-.005)
    path = BallisticPath([-1e-4, 2e-5, 0.], [150., 10., 100.], 2e-6, 'declared Rb control path')
    axis = AnalysisFrequencyAxis.from_hz([1e6])
    return inputs, geometry, path, axis


def test_rb_retains_loop_doppler_and_does_not_double_count_transit():
    i,g,p,axis = rb_fixture()
    f = reduced_ballistic_problem(i,g,p,axis,segments=2)
    q = g.wavevectors_rad_m[1:]-g.wavevectors_rad_m[0]
    np.testing.assert_allclose(f['metadata']['loop_convective_frequency_rad_s'], -q.sum(axis=0)@p.velocity_m_s, atol=2e-6)
    assert abs(f['metadata']['loop_convective_frequency_rad_s']) > 1e5
    assert len(f['reservoirs'].channels)==4
    assert not any('reset' in x.name for x in f['reservoirs'].channels)
    assert np.ptp(f['metadata']['pump_amplitudes_rad_s']) > 1e6
    with pytest.raises(ValueError,match='transit'): reduced_ballistic_problem(replace(i,transit_rate_s_inverse=1.),g,p,axis,segments=2)
    with pytest.raises(ValueError,match='scalar'): reduced_ballistic_problem(replace(i,phase_mismatch_rad_m=1.),g,p,axis,segments=2)


def test_rb_phase_theorem_against_actual_resolve_and_wrong_independent_phase():
    i,g,path,axis = rb_fixture()
    p = reduced_ballistic_wavepacket(i,g,path,axis,segments=2)
    phi = .73
    shifted = reduced_ballistic_wavepacket(i,g,path,axis,segments=2,entry_phase_rad=phi)
    s = np.array([1,-1,-1,1])
    u = np.exp(1j*s*phi)
    for key in ('greater', 'lesser', 'greater_by_source', 'lesser_by_source', 'retarded_response'):
        np.testing.assert_allclose(shifted[key], p[key]*u[:,None]*u.conj()[None,:], rtol=2e-8, atol=2e-26)
    np.testing.assert_allclose(shifted['mean_pulse'], p['mean_pulse']*u, rtol=2e-8, atol=2e-19)
    avg = average_pump_only_entry_phase(p)
    discrete = sum(p['greater']*np.exp(1j*s*phi)[:,None]*np.exp(-1j*s*phi)[None,:]
                   for phi in 2*np.pi*np.arange(8)/8)/8
    np.testing.assert_allclose(avg['conditional_greater'], discrete, atol=2e-34)
    assert np.linalg.norm(avg['mean_outer_phase_average']) > 1e-26
    assert avg['conditional_greater'][0,0,3] != 0  # common beat preserves anomalous pair correlation
    assert avg['conditional_greater'][0,0,1] == 0


def test_boundary_omission_cannot_be_certified_by_noise_positivity():
    i,g,path,axis = rb_fixture()
    p = reduced_ballistic_wavepacket(i,g,path,axis,segments=2)
    assert p['audit']['passed']
    rest = p['greater_by_source'][1:].sum(axis=0)
    assert np.linalg.eigvalsh((rest+rest.conj().swapaxes(-1,-2))/2).min() >= -1e-29
    assert np.linalg.norm(p['greater']-rest)/np.linalg.norm(p['greater']) > .01
    atomic = p['atomic_covariance_by_source'][-1, 0]
    k = atomic.sum(axis=0)-atomic.sum(axis=0).T
    missing = atomic[1:].sum(axis=0)-atomic[1:].sum(axis=0).T
    defect = np.linalg.norm(k-missing)/np.linalg.norm(k)
    assert defect > 1e-6
    assert defect > 100*p['audit']['atomic_commutator_relative_residual']


def test_rb_audit_refuses_overwrite_before_computing(tmp_path):
    from analysis.grand_challenge.rb_transport_audit import write_audit
    output = tmp_path/'existing.json'
    output.write_text('preserve',encoding='utf-8')
    with pytest.raises(FileExistsError): write_audit(output,tmp_path/'new.png')
    assert output.read_text(encoding='utf-8')=='preserve'
