"""Exact mini re-implementation of the FWM pipeline (reference for validating the
closed-form theory). Conventions ported 1:1 from the model, seeded-beat corrected
(Omega_beat = nu_HF + branch*delta).
Target: beat-corrected optimum  Delta=-1.50 GHz, T=110 C  ->
        xi_finite=-8.10 dB, xi_ideal=-15.62 dB, G_s=22.6, gap=0.623, delta=-280 MHz
Reproduced here: G_s=22.59, gap=0.6226, xi=-8.101/-15.619 dB, delta*=-280.0 MHz.

Companion of docs: squeezing_analytic_reconstruction_v1.tex
"""
import numpy as np
from scipy.special import wofz

# ---------- constants (gabes/constants.py, hyperfine.py) ----------
HBAR = 1.054571817e-34; KB = 1.380649e-23; C = 299792458.0
EPS0 = 8.8541878128e-12
NU_D1 = 377.107385690e12; LAM = C / NU_D1
K_VEC = 2*np.pi/LAM; OMEGA_D1 = 2*np.pi*NU_D1
GAMMA = 2*np.pi*5.746e6
NU_HF = 3.035732439e9; OMEGA_HF = 2*np.pi*NU_HF
NU_EHF = 361.58e6; OMEGA_EHF = 2*np.pi*NU_EHF
MASS = 1.4100e-25; I_SAT = 44.84  # W/m^2
DIP = 2.5377e-29
GAMMA_GG0 = 2*np.pi*100e3; XSEC = 1.9e-18
BETA_SELF = 2*np.pi*0.69e-13     # rad/s * m^3
CF2 = {(2,2):10/81, (2,3):35/81, (3,2):35/81, (3,3):28/81}
GPOP = {2:5/12, 3:7/12}; NSUB = 12
DIPSQ = 3.0*DIP**2
DE = 1065.646e6
SHIFT = {(3,3):DE, (3,2):DE-NU_EHF, (2,3):DE+NU_HF, (2,2):DE+NU_HF-NU_EHF}

def rabi(P, w):
    I = 2*P/(np.pi*w*w)
    return GAMMA*np.sqrt(I/(2*I_SAT))

def density(T):
    return 10**(9.318-4040.0/T)/(KB*T)

# ---------- 4-level atom ----------
G1,G2,E2,E3 = 0,1,2,3; NL = 4
A_SC = np.zeros((4,4))
GF = {G1:2, G2:3}; EF = {E2:2, E3:3}
for g in (G1,G2):
    for e in (E2,E3):
        A_SC[g,e] = np.sqrt(3.0*CF2[(GF[g],EF[e])])

def comm_super(H):
    I4 = np.eye(NL)
    return -1j*(np.kron(H, I4) - np.kron(I4, H.T))

def lindblad(gamma_gg, gamma_opt):
    D = np.zeros((16,16), complex); I4 = np.eye(NL)
    kets = np.eye(NL)
    for e in (E2,E3):
        wts = {g: CF2[(GF[g],EF[e])] for g in (G1,G2)}
        tot = sum(wts.values())
        for g in (G1,G2):
            rate = GAMMA*wts[g]/tot
            L = np.sqrt(rate)*np.outer(kets[g], kets[e])
            LdL = L.conj().T@L
            D += np.kron(L, L.conj()) - 0.5*np.kron(LdL,I4) - 0.5*np.kron(I4,LdL.T)
    deph = [(G1,G2,gamma_gg),(G2,G1,gamma_gg)]
    for g in (G1,G2):
        for e in (E2,E3):
            deph += [(g,e,gamma_opt),(e,g,gamma_opt)]
    for (i,j,r) in deph:
        D[i*NL+j, i*NL+j] -= r
    return D

S_V = comm_super(np.diag([0,0,1.0,1.0]).astype(complex))

def H0_build(OpA, Os, delta):
    H = np.zeros((4,4), complex)
    H[G2,G2] = delta; H[E2,E2] = -OMEGA_EHF
    for e in (E2,E3):
        H[G1,e] += OpA*A_SC[G1,e]/2; H[e,G1] += OpA*A_SC[G1,e]/2
        H[G2,e] += Os*A_SC[G2,e]/2;  H[e,G2] += Os*A_SC[G2,e]/2
    return H

def Hp_build(OpB, Oc):
    H = np.zeros((4,4), complex)
    for e in (E2,E3):
        H[e,G1] += Oc*A_SC[G1,e]/2
        H[e,G2] += OpB*A_SC[G2,e]/2
    return H

def floquet(L0, Cp, Cm, Ob, deff_axis):
    n = deff_axis.size
    L0b = L0[None,:,:] - deff_axis[:,None,None]*S_V[None,:,:]
    iO = 1j*Ob*np.eye(16)
    Am = L0b - iO[None]; Ap = L0b + iO[None]
    Cmb = np.broadcast_to(Cm,(n,16,16)); Cpb = np.broadcast_to(Cp,(n,16,16))
    AmiCm = np.linalg.solve(Am, Cmb); ApiCp = np.linalg.solve(Ap, Cpb)
    Aeff = L0b - Cpb@AmiCm - Cmb@ApiCp
    Aeff[:,0,:] = 0
    for s in range(NL): Aeff[:,0,s*NL+s] = 1
    rhs = np.zeros((n,16,1), complex); rhs[:,0,0] = 1
    r0 = np.linalg.solve(Aeff, rhs)
    rp = -(ApiCp@r0)
    return r0[:,:,0].reshape(n,4,4), rp[:,:,0].reshape(n,4,4)

def pol(rho, g):
    return sum(A_SC[g,e]*rho[:,e,g] for e in (E2,E3))

def expm2(M, L):
    s = 0.5*(M[...,0,0]+M[...,1,1])
    q00 = M[...,0,0]-s; q01 = M[...,0,1]; q10 = M[...,1,0]; q11 = M[...,1,1]-s
    c = np.sqrt(-(q00*q11-q01*q10)+0j)
    big = np.abs(c) > 1e-30
    sc = np.where(big, c, 1.0)
    cL = c*L; sL = s*L
    soc = np.where(big, np.sinh(cL)/sc, L*np.ones_like(c))
    out = np.empty_like(M)
    e = np.exp(sL); ch = np.cosh(cL)
    out[...,0,0] = e*(ch+soc*q00); out[...,0,1] = e*soc*q01
    out[...,1,0] = e*soc*q10;      out[...,1,1] = e*(ch+soc*q11)
    return out

def voigt_alpha(beam_GHz, T, only_Fg=None):
    """alpha(beam) [1/m], analytic Voigt == code's unit-area OBE Voigt."""
    N = density(T); geff = GAMMA + BETA_SELF*N; gL = geff/2
    sig = K_VEC*np.sqrt(KB*T/MASS)
    ref = SHIFT[(2,3)]
    x = 2*np.pi*(np.asarray(beam_GHz,float)*1e9 + ref)
    K = np.pi*K_VEC*DIPSQ*N/(HBAR*EPS0)/NSUB
    a = np.zeros_like(x)
    for (Fg,Fe),sh in SHIFT.items():
        if only_Fg is not None and Fg != only_Fg: continue
        z = ((x-2*np.pi*sh)+1j*gL)/(sig*np.sqrt(2))
        V = np.real(wofz(z))/(sig*np.sqrt(2*np.pi))
        a = a + K*GPOP[Fg]*CF2[(Fg,Fe)]*V
    return a

# ---------- pipeline ----------
def run(D_GHz=-1.5, T=383.15, P_pump=0.6, P_probe=8e-6, ls=0.74,
        loss=0.055, qe=0.92, L=12.5e-3, wp=530e-6, ws=330e-6,
        coarse=81, window=0.7, vstep=5.0, vcut=3.0, theta_deg=0.32,
        kappa=0.1, save=None):
    branch = -1
    OpA = rabi(P_pump, wp); Os = rabi(P_probe, ws)
    N = density(T)
    vbar = np.sqrt(8*KB*T/(np.pi*MASS/2))
    ggg = GAMMA_GG0 + N*XSEC*vbar
    gopt = 0.5*BETA_SELF*N
    Ldiss = lindblad(ggg, gopt)
    cnorm = GPOP[3]/NSUB; cls = ls*cnorm
    eta = qe*(1-loss)
    Delta = 2*np.pi*D_GHz*1e9

    center = D_GHz - NU_HF/1e9
    probe = np.linspace(center-window, center+window, coarse)
    delta_ax = 2*np.pi*(probe-center)*1e9

    sig_v = np.sqrt(KB*T/MASS)
    vlim = np.ceil(vcut*sig_v/vstep)*vstep
    v = np.arange(-vlim, vlim+0.5*vstep, vstep)
    wts = np.exp(-v**2/(2*sig_v**2)); wts /= wts.sum()
    lo = Delta - K_VEC*v.max(); hi = Delta - K_VEC*v.min()
    step = K_VEC*vstep
    nde = int(np.ceil((hi-lo)/step))+1
    deff = np.linspace(lo, hi, nde)

    stp = (deff[-1]-deff[0])/(nde-1)
    dv = Delta - K_VEC*v
    idxf = np.clip((dv-deff[0])/stp, 0, nde-1)
    ilo = np.clip(np.floor(idxf).astype(int), 0, nde-2)
    fr = idxf-ilo

    chi = {k: np.zeros((coarse,nde), complex) for k in ("ss","cs","sc","cc")}
    for i,d in enumerate(delta_ax):
        Ob = OMEGA_HF + branch*d   # seeded_sideband_beat: nu_HF + branch*delta
        L0 = comm_super(H0_build(OpA, Os, d)) + Ldiss
        Hp = Hp_build(OpA, 0.0)
        r0, rp = floquet(L0, comm_super(Hp), comm_super(Hp.conj().T), Ob, deff)
        chi["ss"][i] = pol(r0, G2)/Os; chi["cs"][i] = pol(rp, G1)/Os
        L0 = comm_super(H0_build(OpA, 0.0, d)) + Ldiss
        Hp = Hp_build(OpA, Os)
        r0, rp = floquet(L0, comm_super(Hp), comm_super(Hp.conj().T), Ob, deff)
        chi["sc"][i] = pol(r0, G2)/Os; chi["cc"][i] = pol(rp, G1)/Os

    avg = {}
    for k in chi:
        tab = chi[k]
        interp = tab[:,ilo]*(1-fr)[None,:] + tab[:,ilo+1]*fr[None,:]
        avg[k] = interp@wts

    coupling = -2.0*N*cls*DIP**2/(EPS0*HBAR)
    chiph_ss = coupling*avg["ss"]; chiph_cc = coupling*avg["cc"]

    th = np.radians(theta_deg)
    pump_off = 2*np.pi*D_GHz*1e9
    seed_off = 2*np.pi*probe*1e9
    conj_off = 2*pump_off - seed_off
    kp = (OMEGA_D1+pump_off)/C
    ks = (OMEGA_D1+seed_off)/C
    kc = (OMEGA_D1+conj_off)/C
    ns = 1+np.clip(0.5*np.real(chiph_ss), -1e-5, 1e-5)
    nc = 1+np.clip(0.5*np.real(chiph_cc), -1e-5, 1e-5)
    dk = 2*kp - (ks*ns + kc*nc)*np.cos(th)
    kps = ks*ns; kcs = kc*nc

    alpha_ss = np.maximum(K_VEC*np.imag(chiph_ss), 0.0)
    seg_od = float(np.clip(np.nanmedian(alpha_ss)*L, 0, 2))
    nseg = 64
    zf = (np.arange(nseg)+0.5)/nseg
    segp = np.exp(-0.5*seg_od*zf)
    z = (zf-0.5)*L
    wsq = wp**2+ws**2
    spat = np.exp(-(np.abs(z*np.tan(th)))**2/wsq); spat /= spat.max()

    M = np.zeros((coarse,2,2), complex)
    M[:,0,0] = 0.5j*kps*coupling*avg["ss"] + 0.5j*dk
    M[:,0,1] = 0.5j*kps*coupling*avg["sc"]
    M[:,1,0] = -0.5j*kcs*coupling*np.conj(avg["cs"])
    M[:,1,1] = -0.5j*kcs*coupling*np.conj(avg["cc"]) - 0.5j*dk

    dz = L/nseg
    amp = np.zeros((coarse,2), complex); amp[:,0] = np.sqrt(P_probe)
    prem = np.full(coarse, P_pump)
    for cs_, sp_ in zip(segp, spat):
        psc = np.sqrt(np.clip(prem/P_pump, 0, 1))
        Mz = M.copy(); sc = cs_*sp_*psc
        Mz[:,0,1] *= sc; Mz[:,1,0] *= sc
        Tm = expm2(Mz, dz)
        amp = np.einsum("nij,nj->ni", Tm, amp)
        prem = np.maximum(P_pump - np.maximum(np.abs(amp[:,0])**2-P_probe,0)
                          - np.abs(amp[:,1])**2, 0)
    Gs = np.abs(amp[:,0])**2/P_probe; Gc = np.abs(amp[:,1])**2/P_probe

    gp = np.maximum(Gs-1,0); Pc = 0.5*P_pump
    Gs = 1+gp/(1+gp*P_probe/Pc); Gc = np.maximum(Gc,0)/(1+np.maximum(Gc,0)*P_probe/Pc)

    conj_GHz = 2*D_GHz - probe
    od_c = np.clip(voigt_alpha(conj_GHz, T, only_Fg=2)*L, 0, 5)
    od_p = np.clip(voigt_alpha(probe, T, only_Fg=2)*L, 0, 5)
    ps_alpha = voigt_alpha(np.array([D_GHz]), T)[0]
    od_pump = float(np.clip(ps_alpha*L, 0, 50))
    pump_sc = kappa*(1-np.exp(-od_pump))

    def xi(eta_):
        es = eta_*np.exp(-(seg_od+od_p)); ec = eta_*np.exp(-od_c)
        ms = es*Gs; mc = ec*Gc
        cov0 = np.clip(0.5*((Gs+Gc)-(Gs-Gc)**2), 0, np.sqrt(np.maximum(Gs*Gc,0)))
        cov = es*ec*cov0
        w = ms/np.maximum(mc,1e-30)
        var = ms + w*w*mc - 2*w*cov
        S = np.maximum(var/np.maximum(ms+w*w*mc,1e-30) + pump_sc, 1e-30)
        return 10*np.log10(S)

    Sf = xi(eta); Si = xi(1.0)
    gap = Gs-Gc
    valid = np.isfinite(Sf) & (gap>=0.5) & (gap<=1.5)
    Ss = np.where(valid, Sf, np.inf)
    i = int(np.nanargmin(Ss))
    out = dict(i=i, delta_MHz=(probe[i]-center)*1e3, Gs=Gs[i], Gc=Gc[i],
               gap=gap[i], xi_finite=Sf[i], xi_ideal=Si[i], seg_od=seg_od,
               od_conj=od_c[i], od_probe=od_p[i], od_pump=od_pump,
               pump_scatter=pump_sc, N=N, OpA_2pi_MHz=OpA/2/np.pi/1e6,
               ggg_2pi_kHz=ggg/2/np.pi/1e3, gopt_2pi_MHz=gopt/2/np.pi/1e6)
    if save:
        np.savez(save, probe=probe, delta_ax=delta_ax, deff=deff, v=v, wts=wts,
                 **{f"chi_{k}": chi[k] for k in chi},
                 **{f"avg_{k}": avg[k] for k in avg},
                 Gs=Gs, Gc=Gc, Sf=Sf, Si=Si, dk=dk, seg_od=seg_od,
                 od_c=od_c, od_p=od_p, params=np.array([D_GHz,T,OpA,Os,N,ggg,gopt,cls]))
    return out

import tempfile, os
NPZ_PATH = os.path.join(tempfile.gettempdir(), "ref_optimum.npz")

if __name__ == "__main__":
    import json, time
    t0 = time.time()
    r = run(save=NPZ_PATH)
    r = {k: (float(v) if np.isscalar(v) or isinstance(v,(np.floating,np.integer)) else v)
         for k,v in r.items()}
    print(json.dumps(r, indent=1))
    print("elapsed", round(time.time()-t0,1), "s")
