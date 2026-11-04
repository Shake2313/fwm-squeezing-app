"""Physical regression gates for CP dephasing and independent optical scans."""
import sys
from pathlib import Path

import numpy as np
import pytest
from scipy.linalg import expm

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from gabes import atoms, core
from gabes import rydberg_electrometry as electrometry
from gabes.schemes import rydberg
from gabes.schemes.rydberg import MHZ, RydbergEITScheme


def parameters(view="AT electrometry", **updates):
    scheme = RydbergEITScheme()
    return dict(scheme.recommended_defaults(scheme.defaults())[view], **updates)


def choi(matrix, n=4):
    return matrix.reshape(n, n, n, n).transpose(2, 0, 3, 1).reshape(n*n, n*n)


@pytest.mark.parametrize("gd,grf", [(0., 0.), (.25, 0.), (0., 1.), (.25, 1.), (2., .01)])
def test_correlated_dephasing_generator_is_conditionally_completely_positive(gd, grf):
    atom = RydbergEITScheme()._atom(gd*MHZ, grf*MHZ)
    identity = np.eye(4).reshape(16)/2
    orthogonal = np.eye(16)-np.outer(identity, identity)
    conditional = orthogonal @ choi(atom.lindblad/MHZ) @ orthogonal
    assert np.linalg.eigvalsh((conditional+conditional.conj().T)/2).min() >= -2e-12
    # A finite channel check exercises the full decay+phase generator too.
    channel = expm(atom.lindblad * (1./MHZ))
    assert np.linalg.eigvalsh((choi(channel)+choi(channel).conj().T)/2).min() >= -2e-12
    trace = np.eye(4).reshape(16)
    np.testing.assert_allclose(trace @ channel, trace, atol=2e-12)


def test_nested_projectors_impose_consistent_optical_and_rf_coherence_rates():
    atom = RydbergEITScheme()._atom(.3*MHZ, .8*MHZ)
    phase = sum(core.dissipator(operator) for operator in atom.collapse_ops)
    for i, j, expected in [(0, 2, .3), (2, 3, .8), (0, 3, 1.1), (1, 2, .3), (1, 3, 1.1)]:
        np.testing.assert_allclose(-phase[i*4+j, i*4+j]/MHZ, expected, atol=1e-13)
    assert atom.dephasing == ()


@pytest.mark.parametrize("rate", [-1., float("nan"), float("inf")])
def test_invalid_dephasing_rates_rejected(rate):
    scheme = RydbergEITScheme()
    with pytest.raises(ValueError, match="finite and nonnegative"):
        scheme._atom(rate, MHZ)
    with pytest.raises(ValueError, match="finite and nonnegative"):
        scheme._atom(MHZ, rate)


def test_legacy_negative_stationary_state_is_physical_after_cp_correction():
    scheme = RydbergEITScheme()
    params = parameters("EIT", probe_power_uw=.5)
    raw = scheme.compute(params)
    gd = (raw["rydberg_dephasing_mhz"]+raw["transit_mhz"])*MHZ
    grf = raw["rf_dephasing_mhz"]*MHZ
    sk = rydberg._cascade_skeleton()
    legacy = atoms.AtomModel(name="test_legacy_non_CP", n_levels=4, labels=sk["labels"],
        ground=(0,), excited=(1, 2, 3), decay=sk["decay"],
        dephasing=((0,2,gd),(2,0,gd),(0,3,gd),(3,0,gd),(2,3,grf),(3,2,grf)),
        doppler_levels=sk["doppler_levels"], doppler_ratios=sk["doppler_ratios"])
    hamiltonian = scheme._hamiltonian_at_scan_detuning(raw, params, -3*MHZ)
    bad = core.steady_state_from_liouvillian(core.build_liouvillian(hamiltonian, legacy), 4)
    good = core.steady_state_from_liouvillian(core.build_liouvillian(hamiltonian, scheme._atom(gd, grf)), 4)
    assert np.linalg.eigvalsh((bad+bad.conj().T)/2).min() < -1e-6
    assert np.linalg.eigvalsh((good+good.conj().T)/2).min() >= -1e-12
    np.testing.assert_allclose(np.trace(good), 1., atol=1e-12)


@pytest.mark.parametrize("axis", ["Probe detuning", "Coupling detuning"])
def test_stationary_states_remain_positive_across_optical_axis_and_velocity(axis):
    scheme = RydbergEITScheme()
    params = parameters(scan_axis=axis, probe_detuning_mhz=4., coupling_detuning_mhz=-4., probe_power_uw=.5)
    raw = scheme.compute(params)
    atom = scheme._atom((raw["rydberg_dephasing_mhz"]+raw["transit_mhz"])*MHZ,
                        raw["rf_dephasing_mhz"]*MHZ)
    center = scheme._resonance_scan_mhz(params)
    matrices = [core.build_liouvillian(scheme._hamiltonian_at_scan_detuning(raw, params, (center+delta)*MHZ), atom)
                + raw["k_vec"]*velocity*atom.S_v
                for delta in (-6., -3., 0., 3., 6.) for velocity in (-300., 0., 300.)]
    states = core.steady_state_from_liouvillian(np.array(matrices), 4)
    assert np.linalg.eigvalsh((states+states.conj().transpose(0,2,1))/2).min() >= -1e-11


@pytest.mark.parametrize("axis,expected", [
    ("Probe detuning", [2., -3., -3.7]),
    ("Coupling detuning", [-5., -3., -3.7]),
])
def test_optical_detunings_enter_one_and_two_photon_diagonals(axis, expected):
    fields = dict(probe_rabi_mhz=2., coupling_rabi_mhz=3., lo_rabi_mhz=3.7)
    params = parameters(scan_axis=axis, probe_detuning_mhz=5., coupling_detuning_mhz=5., mw_detuning_mhz=.7)
    h = RydbergEITScheme._hamiltonian_at_scan_detuning(fields, params, -2.*MHZ)
    np.testing.assert_allclose(np.diag(h)[1:]/MHZ, expected, atol=1e-12)
    np.testing.assert_allclose(h, h.conj().T)


def test_default_scan_grid_retained_and_nonzero_fixed_counterpart_recenters_it():
    scheme = RydbergEITScheme()
    params = parameters()
    half = max(10., 3*max(params["lo_rabi_mhz"], scheme._coupling_rabi(params)))
    np.testing.assert_array_equal(scheme._scan(params), np.linspace(-half, half, 801)*MHZ)
    shifted = scheme._scan(dict(params, coupling_detuning_mhz=12.))/MHZ
    np.testing.assert_allclose((shifted[0]+shifted[-1])/2, -12., atol=1e-12)
    coupling = scheme._scan(dict(params, scan_axis="Coupling detuning", probe_detuning_mhz=7.))/MHZ
    np.testing.assert_allclose((coupling[0]+coupling[-1])/2, -7., atol=1e-12)
    with pytest.raises(ValueError, match="scan_axis"):
        scheme._scan(dict(params, scan_axis="ambiguous"))


class PointScan(RydbergEITScheme):
    def __init__(self, center):
        self.point = center

    def _scan(self, params):
        return (self.point+np.array([-.01,0.,.01]))*MHZ


@pytest.mark.parametrize("doppler", ["off", "on"])
def test_same_physical_point_has_same_stationary_and_finite_if_response_on_either_axis(monkeypatch, doppler):
    probe = PointScan(2.)
    coupling = PointScan(-1.)
    pp = parameters(scan_axis="Probe detuning", coupling_detuning_mhz=-1., residual_zeeman_mhz=0.)
    pc = parameters(scan_axis="Coupling detuning", probe_detuning_mhz=2., residual_zeeman_mhz=0.)
    # Stationary comparison uses v=0; the finite-IF comparison below separately
    # exercises coherent velocity classes through the selected optical axis.
    rawp, rawc = probe.compute(pp), coupling.compute(pc)
    np.testing.assert_allclose(rawp["chi_bar"][1], rawc["chi_bar"][1], rtol=2e-10, atol=1e-20)
    velocity = np.array([-40.,0.,40.])
    weights = np.array([.2,.6,.2])
    monkeypatch.setattr(rydberg.affine_doppler.AffineScan, "quadrature_classes", lambda *a: (velocity,weights))
    pp, pc = dict(pp,doppler=doppler),dict(pc,doppler=doppler)
    rawp, rawc = dict(rawp,scan=np.array([2.*MHZ])),dict(rawc,scan=np.array([-1.*MHZ]))
    ap = probe._superheterodyne_atomic_phasors(rawp,pp,np.array([2.]),np.ones(1),static_indices=np.array([0]))
    ac = coupling._superheterodyne_atomic_phasors(rawc,pc,np.array([-1.]),np.ones(1),static_indices=np.array([0]))
    np.testing.assert_allclose(ap[1],ac[1],rtol=1e-11,atol=1e-21)
    atom = probe._atom((rawp["rydberg_dephasing_mhz"]+rawp["transit_mhz"])*MHZ, rawp["rf_dephasing_mhz"]*MHZ)
    h = probe._hamiltonian_at_scan_detuning(rawp,pp,2.*MHZ)
    velocities = velocity if doppler=="on" else np.array([0.])
    w = weights if doppler=="on" else np.array([1.])
    direct = []
    for v in velocities:
        L = core.build_liouvillian(h,atom)+rawp["k_vec"]*v*atom.S_v
        direct.append(electrometry.weak_signal_response_from_liouvillian(L,4,40e3*2*np.pi,
            signal_transition=(2,3)).real_observable_phasor_per_angular_rabi(probe._absorption_quadrature_operator()))
    np.testing.assert_allclose(ap[1][0],np.dot(w,direct),rtol=1e-10,atol=1e-21)


def test_atomic_cache_key_includes_axis_and_both_fixed_detunings(monkeypatch):
    scheme = RydbergEITScheme()
    params = parameters()
    raw = scheme.compute(params)
    x,transmission,_=scheme._transmission(raw['chi_bar'],raw,params)
    calls=[]
    def sentinel(*args,**kwargs):
        calls.append(1)
        return np.array([0]),np.array([1j]),"test"
    monkeypatch.setattr(scheme,"_superheterodyne_atomic_phasors",sentinel)
    scheme._cached_superheterodyne_atomic_phasors(raw,params,x,transmission)
    scheme._cached_superheterodyne_atomic_phasors(raw,params,x,transmission)
    assert len(calls)==1
    for changed in (dict(params,probe_detuning_mhz=1.),dict(params,coupling_detuning_mhz=1.),
                    dict(params,scan_axis="Coupling detuning")):
        scheme._cached_superheterodyne_atomic_phasors(raw,changed,x,transmission)
    assert len(calls)==4


def test_features_and_resonance_readout_use_nonzero_center():
    x=np.linspace(-8,8,1601)
    transmission=.5+.4/(1+(x/.7)**2)
    width,contrast=RydbergEITScheme._eit_features(x,transmission)
    shifted=RydbergEITScheme._eit_features(x+12.,transmission,center_mhz=12.)
    np.testing.assert_allclose(shifted,[width,contrast],atol=1e-12)
    scheme=RydbergEITScheme()
    params=parameters("EIT",scan_axis="Coupling detuning",probe_detuning_mhz=12.)
    raw=scheme.compute(params)
    ro=scheme._readout(raw,params)
    index=int(np.argmin(abs(ro['x']+12.)))
    assert ro['resonance_scan_mhz']==-12.
    np.testing.assert_allclose(ro['rydberg_readout']['resonance_transmission'],ro['T_trans'][index])
    assert raw['dephasing_model']=='correlated_phase_lindblad'
    assert ro['rydberg_readout']['scan_axis']=='Coupling detuning'
    doublet=.5+.2/(1+((x-2)/.2)**2)+.2/(1+((x+2)/.2)**2)
    shifted_peaks=scheme._transparency_maxima(x+12.,doublet,3.,center_mhz=12.)
    np.testing.assert_allclose([point[0] for point in shifted_peaks],[10.,14.])


def test_finite_if_export_separates_scan_coordinate_from_physical_detunings():
    scheme=RydbergEITScheme()
    params=parameters(scan_axis="Coupling detuning",probe_detuning_mhz=12.)
    ro=scheme._readout(scheme.compute(params),params)
    result=ro['superhet']
    assert result is not None
    assert result['scan_axis']=='Coupling detuning'
    assert result['optimum_probe_detuning_mhz']==12.
    assert result['optimum_scan_detuning_mhz']==result['optimum_detuning_mhz']
    assert result['optimum_coupling_detuning_mhz']==result['optimum_scan_detuning_mhz']
    np.testing.assert_allclose(result['optimum_two_photon_detuning_mhz'],12.+result['optimum_coupling_detuning_mhz'])
    exported=ro['rydberg_readout']['finite_if']
    assert exported['optimum_probe_detuning_mhz']==12.
    assert exported['optimum_coupling_detuning_mhz']==result['optimum_coupling_detuning_mhz']


def test_live_probe_anchor_scales_power_and_diameter_and_is_exported():
    scheme=RydbergEITScheme()
    params=parameters(probe_rabi_ref_mhz=3.)
    assert scheme._probe_rabi(params)==3.
    assert scheme._probe_rabi(dict(params,probe_power_uw=24.))==6.
    assert scheme._probe_rabi(dict(params,beam_diameter_mm=.3))==1.5
    raw=scheme.compute(params)
    assert raw['probe_rabi_ref_mhz']==3.
    assert raw['probe_rabi_mhz']==3.
    declared={p.name:p for p in scheme.param_schema()}
    assert declared['probe_rabi_ref_mhz'].recompute
    assert declared['scan_axis'].recompute
    assert declared['coupling_detuning_mhz'].recompute
    assert declared['probe_detuning_mhz'].recompute
