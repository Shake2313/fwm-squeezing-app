"""Final numbers for squeezing_analytic_reconstruction_v1.tex:
 (1) full analytic pipeline -> G_s, G_c, gap, xi_finite, xi_ideal,
 (2) reference (exact chi) numbers,
 (3) Faddeeva closed-form demo for the linear chi_cc direct term,
 (4) gap-evolution identity decomposition,
 (5) noise-budget decomposition at the reference gains.

The main routine refreshes the reference NPZ before reading it.
Beat-corrected optimum (Delta=-1.5 GHz, T=110 C), results (2026-07-06):
  exact (full Floquet) @ delta*=-280 MHz : Gs=22.59, gap=0.623, xi=-8.10 / -15.62 dB
  internal simulation                    : Gs=22.6,  gap=0.623, xi=-8.10 / -15.62 dB
  leading pencil truncation @ -280 MHz   : Gs=3.5,  gap=1.25, xi=-4.29 / -5.57 dB
    (right optimum + phase, but 0.5-0.8x magnitude -> gain undershoots)
"""
import json
import numpy as np
from scipy.special import wofz
if __package__:
    from . import theory_final as TH
    from .ref_solver import (
        rabi, density, voigt_alpha, expm2, run as ref_run,
        KB, MASS, K_VEC, OMEGA_D1, C, EPS0, HBAR, DIP,
        GAMMA_GG0, XSEC, BETA_SELF, NU_HF, GPOP, NSUB, GAMMA, NPZ_PATH,
    )
else:  # Direct script execution.
    import theory_final as TH
    from ref_solver import (
        rabi, density, voigt_alpha, expm2, run as ref_run,
        KB, MASS, K_VEC, OMEGA_D1, C, EPS0, HBAR, DIP,
        GAMMA_GG0, XSEC, BETA_SELF, NU_HF, GPOP, NSUB, GAMMA, NPZ_PATH,
    )

def run_pipeline(chi_fn, D_GHz=-1.5, T=383.15, P_pump=0.6, P_probe=8e-6,
                 ls=0.74, loss=0.055, qe=0.92, L=12.5e-3, wp=530e-6, ws=330e-6,
                 coarse=81, window=0.7, vstep=5.0, vcut=3.0, theta_deg=0.32,
                 kappa=0.1):
    OpA = rabi(P_pump, wp); Os = rabi(P_probe, ws)
    N = density(T)
    vbar = np.sqrt(8*KB*T/(np.pi*MASS/2))
    ggg = GAMMA_GG0 + N*XSEC*vbar; gopt = 0.5*BETA_SELF*N
    Go = GAMMA/2 + gopt
    cls = ls*GPOP[3]/NSUB
    eta = qe*(1-loss)
    Delta = 2*np.pi*D_GHz*1e9
    center = D_GHz - NU_HF/1e9
    probe = np.linspace(center-window, center+window, coarse)
    delta_ax = 2*np.pi*(probe-center)*1e9
    sig_v = np.sqrt(KB*T/MASS)
    vlim = np.ceil(vcut*sig_v/vstep)*vstep
    v = np.arange(-vlim, vlim+0.5*vstep, vstep)
    wts = np.exp(-v**2/(2*sig_v**2)); wts /= wts.sum()
    deffs = Delta - K_VEC*v
    avg = {k: np.zeros(coarse, complex) for k in ("ss","cs","sc","cc")}
    for i,dlt in enumerate(delta_ax):
        th = np.zeros((4, v.size), complex)
        for j,De in enumerate(deffs):
            th[:,j] = chi_fn(OpA, Os, dlt, De, Go, ggg)
        for k,name in enumerate(("ss","cs","sc","cc")):
            avg[name][i] = th[k]@wts
    coupling = -2.0*N*cls*DIP**2/(EPS0*HBAR)
    chiph_ss = coupling*avg["ss"]; chiph_cc = coupling*avg["cc"]
    th_r = np.radians(theta_deg)
    pump_off = 2*np.pi*D_GHz*1e9; seed_off = 2*np.pi*probe*1e9
    conj_off = 2*pump_off-seed_off
    kp=(OMEGA_D1+pump_off)/C; ks=(OMEGA_D1+seed_off)/C; kc=(OMEGA_D1+conj_off)/C
    ns = 1+np.clip(0.5*np.real(chiph_ss),-1e-5,1e-5)
    nc = 1+np.clip(0.5*np.real(chiph_cc),-1e-5,1e-5)
    dk = 2*kp-(ks*ns+kc*nc)*np.cos(th_r)
    alpha_ss = np.maximum(K_VEC*np.imag(chiph_ss),0.0)
    seg_od = float(np.clip(np.nanmedian(alpha_ss)*L,0,2))
    nseg=64; zf=(np.arange(nseg)+0.5)/nseg
    segp=np.exp(-0.5*seg_od*zf)
    z=(zf-0.5)*L; wsq=wp**2+ws**2
    spat=np.exp(-(np.abs(z*np.tan(th_r)))**2/wsq); spat/=spat.max()
    M = np.zeros((coarse,2,2), complex)
    M[:,0,0]=0.5j*ks*ns*coupling*avg["ss"]+0.5j*dk
    M[:,0,1]=0.5j*ks*ns*coupling*avg["sc"]
    M[:,1,0]=-0.5j*kc*nc*coupling*np.conj(avg["cs"])
    M[:,1,1]=-0.5j*kc*nc*coupling*np.conj(avg["cc"])-0.5j*dk
    dz=L/nseg
    amp=np.zeros((coarse,2),complex); amp[:,0]=np.sqrt(P_probe)
    prem=np.full(coarse,P_pump)
    for cs_,sp_ in zip(segp,spat):
        psc=np.sqrt(np.clip(prem/P_pump,0,1))
        Mz=M.copy(); sc=cs_*sp_*psc
        Mz[:,0,1]*=sc; Mz[:,1,0]*=sc
        Tm=expm2(Mz,dz)
        amp=np.einsum("nij,nj->ni",Tm,amp)
        prem=np.maximum(P_pump-np.maximum(np.abs(amp[:,0])**2-P_probe,0)
                        -np.abs(amp[:,1])**2,0)
    Gs=np.abs(amp[:,0])**2/P_probe; Gc=np.abs(amp[:,1])**2/P_probe
    gp=np.maximum(Gs-1,0); Pc=0.5*P_pump
    Gs=1+gp/(1+gp*P_probe/Pc); Gc=np.maximum(Gc,0)/(1+np.maximum(Gc,0)*P_probe/Pc)
    conj_GHz=2*D_GHz-probe
    od_c=np.clip(voigt_alpha(conj_GHz,T,only_Fg=2)*L,0,5)
    od_p=np.clip(voigt_alpha(probe,T,only_Fg=2)*L,0,5)
    od_pump=float(np.clip(voigt_alpha(np.array([D_GHz]),T)[0]*L,0,50))
    pump_sc=kappa*(1-np.exp(-od_pump))
    def xi(eta_):
        es=eta_*np.exp(-(seg_od+od_p)); ec=eta_*np.exp(-od_c)
        ms=es*Gs; mc=ec*Gc
        cov0=np.clip(0.5*((Gs+Gc)-(Gs-Gc)**2),0,np.sqrt(np.maximum(Gs*Gc,0)))
        cov=es*ec*cov0
        w=ms/np.maximum(mc,1e-30)
        var=ms+w*w*mc-2*w*cov
        return 10*np.log10(np.maximum(var/np.maximum(ms+w*w*mc,1e-30)+pump_sc,1e-30))
    Sf=xi(eta); Si=xi(1.0)
    return dict(probe=probe, center=center, Gs=Gs, Gc=Gc, Sf=Sf, Si=Si,
                seg_od=seg_od, od_c=od_c, od_p=od_p, od_pump=od_pump,
                pump_sc=pump_sc, M=M, dk=dk, segp=segp, spat=spat, L=L,
                avg=avg, coupling=coupling, ks=ks, kc=kc, ns=ns, nc=nc)

def main():
    out = {}
    r = run_pipeline(TH.response)
    i0 = 24  # delta = -280 MHz on the +-0.7 GHz / 81-point grid
    gap = r["Gs"] - r["Gc"]
    valid = np.isfinite(r["Sf"]) & (gap >= 0.5) & (gap <= 1.5)
    scores = np.where(valid, r["Sf"], np.inf)
    iam = int(np.nanargmin(scores)) if valid.any() else int(np.nanargmin(r["Sf"]))
    out["theory_at_297"] = dict(
        Gs=float(r["Gs"][i0]), Gc=float(r["Gc"][i0]), gap=float(gap[i0]),
        xi_f=float(r["Sf"][i0]), xi_i=float(r["Si"][i0]))
    out["theory_argmin"] = dict(
        delta=float((r["probe"][iam] - r["center"]) * 1e3),
        Gs=float(r["Gs"][iam]), gap=float(gap[iam]),
        xi_f=float(r["Sf"][iam]), xi_i=float(r["Si"][iam]))
    out["noise_pieces"] = dict(
        seg_od=r["seg_od"], od_c=float(r["od_c"][i0]),
        od_p=float(r["od_p"][i0]), od_pump=r["od_pump"], pump_sc=r["pump_sc"])

    ref = ref_run(save=NPZ_PATH)
    out["ref"] = {k: float(ref[k]) for k in
                  ("Gs", "Gc", "gap", "xi_finite", "xi_ideal", "delta_MHz",
                   "seg_od", "od_conj", "od_pump", "pump_scatter")}

    Gs0, Gc0 = ref["Gs"], ref["Gc"]
    odc, odp_lin, sod = ref["od_conj"], ref["od_probe"], ref["seg_od"]
    pump_scatter = ref["pump_scatter"]

    def budget(eta):
        eta_s = eta * np.exp(-(sod + odp_lin))
        eta_c = eta * np.exp(-odc)
        mean_s, mean_c = eta_s * Gs0, eta_c * Gc0
        cov0 = min(0.5 * ((Gs0 + Gc0) - (Gs0 - Gc0) ** 2), np.sqrt(Gs0 * Gc0))
        cov = eta_s * eta_c * cov0
        weight = mean_s / mean_c
        variance = mean_s + weight * weight * mean_c - 2 * weight * cov
        balanced = variance / (mean_s + weight * weight * mean_c)
        return dict(
            S_bal=float(balanced), S_tot=float(balanced + pump_scatter),
            xi=float(10 * np.log10(balanced + pump_scatter)))

    out["budget_finite"] = budget(0.92 * (1 - 0.055))
    out["budget_ideal"] = budget(1.0)
    out["budget_ideal_noscatter_xi"] = float(
        10 * np.log10(out["budget_ideal"]["S_bal"]))

    T = 383.15
    N = density(T)
    vbar = np.sqrt(8 * KB * T / (np.pi * MASS / 2))
    ggg = GAMMA_GG0 + N * XSEC * vbar
    gopt = 0.5 * BETA_SELF * N
    Go = GAMMA / 2 + gopt
    OpA = rabi(0.6, 530e-6)
    sig_v = np.sqrt(KB * T / MASS)
    sig_w = K_VEC * sig_v
    branch = -1
    Delta = 2 * np.pi * (-1.5e9)
    delta = 2 * np.pi * (-280e6)
    beat = 2 * np.pi * NU_HF + branch * delta
    vlim = np.ceil(3 * sig_v / 5.0) * 5.0
    velocity = np.arange(-vlim, vlim + 2.5, 5.0)
    weights = np.exp(-velocity ** 2 / (2 * sig_v ** 2))
    weights /= weights.sum()
    quad = 0.0
    fadd = 0.0
    carrier = TH.carrier(OpA, delta, Delta, Go, ggg)
    for excited in (2, 3):
        amplitude = np.sqrt(3 * TH.CF2[(2, excited)])
        x0 = Delta + TH.EHF[excited] + beat
        population = carrier["w1"][excited]
        quad += (amplitude ** 2 * population
                 / (2 * ((x0 - K_VEC * velocity) + 1j * Go))) @ weights
        z = (x0 + 1j * Go) / (sig_w * np.sqrt(2))
        fadd -= (amplitude ** 2 * population / 2
                 * (1j * np.sqrt(np.pi / 2) / sig_w) * wofz(z))
    out["faddeeva_demo"] = dict(
        quad=[float(quad.real), float(quad.imag)],
        faddeeva=[float(fadd.real), float(fadd.imag)])

    with np.load(NPZ_PATH) as data:
        avg_ref = {k: data[f"avg_{k}"][i0] for k in ("ss", "cs", "sc", "cc")}
        probe0 = float(data["probe"][i0])
        dk_v = float(data["dk"][i0])
    line_strength = 0.74 * GPOP[3] / NSUB
    coupling = -2.0 * N * line_strength * DIP ** 2 / (EPS0 * HBAR)
    pump_off = 2 * np.pi * (-1.5e9)
    seed_off = 2 * np.pi * probe0 * 1e9
    conj_off = 2 * pump_off - seed_off
    ks = (OMEGA_D1 + seed_off) / C
    kc = (OMEGA_D1 + conj_off) / C
    a = 0.5j * ks * coupling * avg_ref["ss"] + 0.5j * dk_v
    d = -0.5j * kc * coupling * np.conj(avg_ref["cc"]) - 0.5j * dk_v
    b = 0.5j * ks * coupling * avg_ref["sc"]
    c = -0.5j * kc * coupling * np.conj(avg_ref["cs"])
    out["gap_identity"] = dict(
        Re_a=float(np.real(a)), Re_d=float(np.real(d)), abs_b=float(abs(b)),
        abs_c=float(abs(c)),
        b_minus_cstar=[float(np.real(b - np.conj(c))),
                       float(np.imag(b - np.conj(c)))],
        dk=dk_v, L=12.5e-3)

    print(json.dumps(out, indent=1))
    return out


if __name__ == "__main__":
    main()
