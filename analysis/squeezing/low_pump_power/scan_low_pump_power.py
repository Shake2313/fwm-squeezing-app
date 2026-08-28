"""
Why 100-200 mW of pump is enough for 85Rb D1 double-Lambda FWM squeezing.

Motivated by U. Jain, J. Choi, C. Hull, A. M. Marino, Opt. Lett. 50, 5165 (2025)
("Compact fiber-coupled narrowband two-mode squeezed light source"), which
reports -7.2 dB of intensity-difference squeezing before the output fibers
(-4.4 dB after) from a 135 mW pump -- roughly a quarter of the 550-600 mW used
in the earlier 85Rb FWM squeezing experiments -- by shrinking the pump/seed
1/e^2 waist DIAMETER to ~0.6 mm (w = 300 um radius) in a 12 mm cell at ~99 C.

THE ARGUMENT HAS THREE INDEPENDENT LEGS, all computed here.

(1) PUMP POWER IS NOT A PARAMETER OF THE THEORY.
    In the closed-form chain of `docs/FWM physics and analytic reconstruction/
    squeezing_analytic_reconstruction_v2.tex`, the pump enters the Doppler-
    averaged Floquet response only as the Rabi frequency
        Omega_pump = Gamma*sqrt(I / 2 I_sat),   I = 2P/(pi w^2)
    (`gabes.constants.rabi_freq`), so it enters ONLY as P/w^2. The classical
    coupled-mode matrix M_cl, its trace/discriminant pair (s, q), the transfer
    T_cl = exp(M_cl L) and hence G_s, G_c inherit that dependence.
    => ISO-SQUEEZING CURVES IN THE (P, w) PLANE ARE P proportional to w^2.
    Section `invariance` verifies this BITWISE inside GABES: the small-signal
    gain array is float-identical for (600 mW, 530 um) and (192.2 mW, 300 um).
    The single exception in the engine is the Manley-Rowe pump-photon budget
    (`observables.pump_depletion_saturation`), which depends on the absolute
    pump flux and only ever LIMITS gain -- quantified in `depletion`.

(2) THE PUBLISHED RECORD ALREADY SITS AT ONE INTENSITY.
    Section `published`: Dowran 2018 (550 mW, w=500 um) = 140 W/cm^2,
    Sim 2025 (600 mW, w=530 um) = 136 W/cm^2, Jain 2025 (135-250 mW, w=300 um)
    = 95-177 W/cm^2. The Jain power window BRACKETS the high-power experiments
    in pump Rabi frequency. "Low power" is a smaller beam, not a new regime.

(3) THE GAIN THAT SQUEEZING NEEDS IS SMALL, AND GAIN IS DEEPLY SATURATED IN P.
    Section `noise_law`: the twin-beam law xi = 10log10[eta/(2G-1) + (1-eta)]
    (`observables.intensity_difference_squeezing_dB`) reproduces BOTH published
    squeezing numbers from their measured gains, and says -7 dB needs only
    G ~ 8. Section `power_scan`: at Omega/Gamma ~ 100 the engine's gain obeys
    dln(G-1)/dlnP ~ 0.1, so a 4.4x power cut costs a few percent of gain, i.e.
    < 0.1 dB of squeezing.

WHAT THIS SCRIPT DOES *NOT* CLAIM. The engine's ABSOLUTE gain is not
predictive: v2 Sec. "Model point and experimental comparison" records +7500 %
on G_p at the Sim operating point, and G = cosh^2(qL) turns a ~2x error in the
gain-length product into a ~50x error in G. Every headline number below is
therefore either (a) an exact structural invariance, (b) a logarithmic
sensitivity dln/dln, or (c) the measured-gain-in / squeezing-out noise law --
none of which need the absolute gain scale. Absolute engine gains are printed
as context, flagged as such, and never used as a prediction.

MEAN-FIELD LAYER = GABES legacy-fidelity `G_s_smallsignal`, peak over the
two-photon detuning delta (located coarsely, then refined). Legacy is used on
purpose: it is the layer in which the P/w^2 invariance is exact, so the
sensitivities are not contaminated by Ultra's segmented-overlap and dynamic
depletion corrections -- which are measured separately in `invariance`.
NOISE LAYER = the closed-form twin-beam law, per v2's layer table (the
microscopic quantum-Langevin layer is still pending there).

Usage:
    python analysis/squeezing/low_pump_power/scan_low_pump_power.py
    python analysis/squeezing/low_pump_power/scan_low_pump_power.py --quick
    python analysis/squeezing/low_pump_power/scan_low_pump_power.py --out DIR
"""
import argparse
import json
import math
import os
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", str(os.cpu_count() or 1))

# Korean Windows consoles are cp949; keep stdout ASCII and force UTF-8 on files.
for _stream in (sys.stdout, sys.stderr):
    try:
        _stream.reconfigure(encoding="utf-8", errors="replace")
    except Exception:
        pass

import numpy as np

ROOT = Path(__file__).resolve()
for cand in [ROOT, *ROOT.parents]:
    if (cand / "gabes").is_dir() and (cand / "tests").is_dir():
        REPO = cand
        break
else:  # pragma: no cover
    raise SystemExit("Could not locate the GABES repo root (needs gabes/ + tests/).")
sys.path.insert(0, str(REPO))

from gabes import atoms, constants, hyperfine, observables      # noqa: E402
from gabes.core import blas_single_thread                       # noqa: E402
from gabes.schemes import fwm                                   # noqa: E402

# --------------------------------------------------------------------------
# Operating points.  Beam sizes are 1/e^2 intensity RADII (v2 Sec. "Field,
# beam, and susceptibility normalization" flags the radius/diameter trap).
# --------------------------------------------------------------------------
SIM = dict(name="Sim 2025 (Sci. Rep. 15, 7727)", D_GHz=0.9, T=394.15,
           P_pump=0.600, P_probe=8e-6, w_pump=530e-6, w_probe=330e-6,
           L=12.5e-3, angle=0.32, G_measured=111.0 / 8.0, xi_measured=-7.8)
JAIN = dict(name="Jain 2025 (Opt. Lett. 50, 5165)", D_GHz=0.8, T=372.15,
            P_pump=0.135, P_probe=10e-6, w_pump=300e-6, w_probe=300e-6,
            L=12.0e-3, angle=0.40, G_measured=None, xi_measured=-7.2)
DOWRAN = dict(name="Dowran 2018 (Optica 5, 628)", P_pump=0.550, w_pump=500e-6,
              xi_measured=-9.0)

ETA_EXT = 0.8694          # v2 Table "model layer comparison": QE 0.92 x (1-5.5%)
LS_RESIDUAL = 0.74        # inherited reference residual (tests/test_regression)
LAMBDA_D1 = 794.98e-9     # 85Rb D1

# transit-time floor anchor: GAMMA_GG is coded as the residual ground-Raman
# rate "~ v_mp / w for the FWM beams" at the Sim geometry.
W_TRANSIT_REF = math.sqrt(530e-6 * 330e-6)
T_TRANSIT_REF = 394.15


def v_mp(T):
    """Most probable speed [m/s]."""
    return math.sqrt(2.0 * constants.KB * float(T) / constants.MASS_85RB)


def intensity_w_cm2(P, w):
    return 2.0 * float(P) / (math.pi * float(w) ** 2) * 1e-4


def rayleigh_mm(w, wavelength=LAMBDA_D1):
    return math.pi * float(w) ** 2 / float(wavelength) * 1e3


def xi_from_gain(G, eta=ETA_EXT):
    """Twin-beam intensity-difference squeezing [dB] for a parametric gain G.

    Delegated to `observables.intensity_difference_squeezing_dB` with the
    lossless twin-beam relation G_c = G_s - 1, i.e. S_ideal = 1/(2G-1)."""
    G = np.asarray(G, dtype=float)
    return float(observables.intensity_difference_squeezing_dB(
        G, np.maximum(G - 1.0, 0.0), eta)) if G.ndim == 0 else \
        observables.intensity_difference_squeezing_dB(
            G, np.maximum(G - 1.0, 0.0), eta)


def gain_for_xi(xi_db, eta=ETA_EXT):
    """Invert the twin-beam law: parametric gain needed for a target xi [dB]."""
    S = 10.0 ** (float(xi_db) / 10.0)
    S_ideal = (S - (1.0 - eta)) / eta
    if S_ideal <= 0.0:
        return float("inf")           # below the detection floor 10log10(1-eta)
    return 0.5 * (1.0 / S_ideal + 1.0)


# --------------------------------------------------------------------------
# Mean-field layer
# --------------------------------------------------------------------------
FID = {"legacy": dict(phase_detail=fwm.PHASE_LEGACY),
       "ultra": dict(phase_detail=fwm.PHASE_ULTRA,
                     model_fidelity=fwm.FIDELITY_ULTRA)}


def spectrum(D_GHz, *, T, P_pump, P_probe, w_pump, w_probe, L, angle,
             fid="legacy", ls=LS_RESIDUAL, dmin=-100.0, dmax=450.0, npts=751,
             vstep=5.0, vcut=3.0):
    """One probe scan; `delta_mhz` is the two-photon detuning off branch center."""
    center = fwm.branch_center_GHz(D_GHz, -1)
    with blas_single_thread():
        s = fwm.compute_spectrum(
            D_GHz, T=T, P_pump=P_pump, P_probe=P_probe,
            w_pump=w_pump, w_probe=w_probe, L=L, pump_probe_angle_deg=angle,
            line_strength=ls, loss_frac=0.055, qe=fwm.QE_DETECTOR,
            coarse_points=npts, fine_points=0,
            scan_min=center + dmin * 1e-3, scan_max=center + dmax * 1e-3,
            velocity_step=vstep, velocity_cutoff=vcut, branch=-1, **FID[fid])
    s["delta_mhz"] = (s["probe_axis_GHz"] - center) * 1e3
    return s


def gain_peak(cfg, *, fid="legacy", ls=LS_RESIDUAL, refine_half=6.0,
              refine_n=241, **over):
    """Peak small-signal seed gain over delta, coarse-located then refined.

    `G_s_smallsignal` is used rather than `G_s` because it is the layer in
    which the P/w^2 invariance is exact -- the Manley-Rowe cap that separates
    the two is reported on its own in `depletion`.
    """
    kw = {k: cfg[k] for k in
          ("T", "P_pump", "P_probe", "w_pump", "w_probe", "L", "angle")}
    kw.update(over)
    D = kw.pop("D_GHz", cfg["D_GHz"])
    s = spectrum(D, fid=fid, ls=ls, **kw)
    j = int(np.nanargmax(s["G_s_smallsignal"]))
    d0 = float(s["delta_mhz"][j])
    edge = bool(j <= 1 or j >= s["G_s_smallsignal"].size - 2)
    r = spectrum(D, fid=fid, ls=ls, dmin=d0 - refine_half, dmax=d0 + refine_half,
                 npts=refine_n, **kw)
    k = int(np.nanargmax(r["G_s_smallsignal"]))
    G_ss = float(r["G_s_smallsignal"][k])
    return dict(G_ss=G_ss, G_sat=float(r["G_s"][k]), G_c=float(r["G_c"][k]),
                delta_mhz=float(r["delta_mhz"][k]),
                xi_law_db=xi_from_gain(G_ss),
                coarse_edge=edge, no_resonance=bool(G_ss <= 1.0))


def transit_scaled_atom(w_eff):
    """`fwm.collisional_atom` with the ground-Raman floor rescaled as v_mp/w.

    GABES codes GAMMA_GG as a fixed 2pi x 100 kHz residual anchored at the Sim
    beam geometry, so shrinking the beam in the engine buys the intensity
    without paying the transit-time cost.  This monkeypatch restores the cost
    (analysis only; the engine is not modified)."""
    def collisional_atom(T, density=None):
        if density is None:
            density = hyperfine.number_density(T)
        floor = (constants.GAMMA_GG * (v_mp(T) / v_mp(T_TRANSIT_REF))
                 * (W_TRANSIT_REF / float(w_eff)))
        gamma_gg = constants.ground_coherence_dephasing(T, density, floor=floor)
        gamma_opt = 0.5 * (hyperfine.self_broadened_gamma(density) - constants.GAMMA)
        return atoms.double_lambda_rb85(gamma_gg=gamma_gg,
                                        gamma_opt=max(gamma_opt, 0.0))
    return collisional_atom


def transit_floor_khz(T, w_eff):
    return (constants.GAMMA_GG * (v_mp(T) / v_mp(T_TRANSIT_REF))
            * (W_TRANSIT_REF / float(w_eff)) / (2 * math.pi) / 1e3)


# ==========================================================================
# Section 1 -- the published operating points, in intensity
# ==========================================================================
def sec_published():
    w_ref = JAIN["w_pump"]
    rows = []
    for label, P, w, T, L, xi in (
            ("Dowran 2018, 550 mW", DOWRAN["P_pump"], DOWRAN["w_pump"],
             None, None, DOWRAN["xi_measured"]),
            ("Sim 2025, 600 mW", SIM["P_pump"], SIM["w_pump"], SIM["T"],
             SIM["L"], SIM["xi_measured"]),
            ("Jain 2025, 135 mW (used)", JAIN["P_pump"], JAIN["w_pump"],
             JAIN["T"], JAIN["L"], JAIN["xi_measured"]),
            ("Jain 2025, 250 mW (window top)", 0.250, JAIN["w_pump"],
             JAIN["T"], JAIN["L"], None)):
        Om = constants.rabi_freq(P, w)
        rows.append(dict(
            label=label, P_mW=P * 1e3, w_um=w * 1e6,
            I_W_cm2=intensity_w_cm2(P, w),
            Omega_2pi_MHz=Om / (2 * math.pi) / 1e6,
            Omega_over_Gamma=Om / constants.GAMMA,
            T_C=None if T is None else T - 273.15,
            L_mm=None if L is None else L * 1e3,
            N_m3=None if T is None else float(hyperfine.number_density(T)),
            xi_measured_dB=xi,
            # what the same pump INTENSITY costs at the Jain waist
            P_equiv_at_300um_mW=P * (w_ref / w) ** 2 * 1e3,
            rayleigh_mm=rayleigh_mm(w)))
    I = [r["I_W_cm2"] for r in rows]
    Om = [r["Omega_2pi_MHz"] for r in rows]
    return dict(rows=rows, I_span_W_cm2=[min(I), max(I)],
                Omega_span_2pi_MHz=[min(Om), max(Om)],
                Omega_ratio_jain_over_sim=Om[2] / Om[1],
                I_ratio_jain_over_sim=I[2] / I[1])


# ==========================================================================
# Section 2 -- the P/w^2 invariance, verified inside the engine
# ==========================================================================
def sec_invariance():
    base = dict(SIM)
    ref = dict(T=base["T"], P_probe=base["P_probe"], w_probe=base["w_probe"],
               L=base["L"], angle=base["angle"])
    out = {"note": "pump waist scaled at fixed seed waist/power, so any residual "
                   "difference is a pump-side effect only"}
    for fid in ("legacy", "ultra"):
        s0 = spectrum(base["D_GHz"], P_pump=base["P_pump"],
                      w_pump=base["w_pump"], fid=fid, **ref)
        rows = []
        for w in (300e-6, 400e-6, 700e-6):
            P = base["P_pump"] * (w / base["w_pump"]) ** 2
            s1 = spectrum(base["D_GHz"], P_pump=P, w_pump=w, fid=fid, **ref)
            g0, g1 = s0["G_s_smallsignal"], s1["G_s_smallsignal"]
            rows.append(dict(
                w_um=w * 1e6, P_mW=P * 1e3,
                Omega_2pi_MHz=constants.rabi_freq(P, w) / (2 * math.pi) / 1e6,
                smallsignal_bitwise_identical=bool(np.array_equal(g0, g1)),
                smallsignal_peak_ref=float(np.nanmax(g0)),
                smallsignal_peak=float(np.nanmax(g1)),
                smallsignal_peak_rel_dev=float(np.nanmax(g1) / np.nanmax(g0) - 1.0),
                saturated_peak_ref=float(np.nanmax(s0["G_s"])),
                saturated_peak=float(np.nanmax(s1["G_s"])),
                saturated_peak_rel_dev=float(
                    np.nanmax(s1["G_s"]) / np.nanmax(s0["G_s"]) - 1.0)))
        out[fid] = rows
    # Ultra's residual w-dependence is the crossed-Gaussian overlap.  Scaling
    # the crossing angle with the waist would keep the walk-off/waist ratio
    # fixed, but theta ALSO sets Delta_k_z at Ultra, so this is not a clean
    # overlap compensator -- recorded here as a caveat, not as a correction.
    w, lam = 300e-6, 300e-6 / SIM["w_pump"]
    s0 = spectrum(SIM["D_GHz"], P_pump=SIM["P_pump"], w_pump=SIM["w_pump"],
                  fid="ultra", **ref)
    ref_a = dict(ref)
    ref_a["angle"] = SIM["angle"] * lam
    s1 = spectrum(SIM["D_GHz"], P_pump=SIM["P_pump"] * lam ** 2, w_pump=w,
                  fid="ultra", **ref_a)
    out["angle_compensated"] = dict(
        w_um=w * 1e6, P_mW=SIM["P_pump"] * lam ** 2 * 1e3,
        angle_deg=SIM["angle"] * lam,
        smallsignal_peak_ref=float(np.nanmax(s0["G_s_smallsignal"])),
        smallsignal_peak=float(np.nanmax(s1["G_s_smallsignal"])),
        smallsignal_peak_rel_dev=float(
            np.nanmax(s1["G_s_smallsignal"]) / np.nanmax(s0["G_s_smallsignal"]) - 1.0))
    # geometric walk-off, closed form
    out["walk_off"] = [
        dict(label=lab, w_um=w_ * 1e6, angle_deg=a, L_mm=L_ * 1e3,
             end_separation_um=0.5 * L_ * math.tan(math.radians(a)) * 1e6,
             end_overlap=math.exp(-(0.5 * L_ * math.tan(math.radians(a))) ** 2
                                  / (w_ ** 2 + wp ** 2)))
        for lab, w_, wp, a, L_ in (
            ("Sim 2025", SIM["w_pump"], SIM["w_probe"], SIM["angle"], SIM["L"]),
            ("Jain 2025", JAIN["w_pump"], JAIN["w_probe"], JAIN["angle"], JAIN["L"]))]
    return out


# ==========================================================================
# Section 3 -- the twin-beam noise law, checked against both experiments
# ==========================================================================
def sec_noise_law():
    knee = [dict(xi_target_dB=x,
                 G_required_eta0p8694=gain_for_xi(x, 0.8694),
                 G_required_eta0p90=gain_for_xi(x, 0.90),
                 G_required_eta0p95=gain_for_xi(x, 0.95))
            for x in (-3.0, -4.0, -4.4, -5.0, -6.0, -7.0, -7.2, -7.8, -8.0, -8.5)]
    sim_pred = xi_from_gain(SIM["G_measured"], ETA_EXT)
    jain_G = gain_for_xi(JAIN["xi_measured"], ETA_EXT)
    # sensitivity of xi to the gain, around the operating gain
    sens = []
    for G in (4.0, 6.0, 8.0, 10.0, 13.875, 20.0, 40.0):
        sens.append(dict(G=G, xi_dB=xi_from_gain(G, ETA_EXT),
                         d_xi_per_octave_dB=xi_from_gain(2 * G, ETA_EXT)
                         - xi_from_gain(G, ETA_EXT)))
    return dict(
        eta_ext=ETA_EXT,
        detection_floor_dB=10.0 * math.log10(1.0 - ETA_EXT),
        gain_knee=knee,
        sim_check=dict(G_measured=SIM["G_measured"], xi_law_dB=sim_pred,
                       xi_measured_dB=SIM["xi_measured"],
                       residual_dB=sim_pred - SIM["xi_measured"]),
        jain_implied=dict(xi_measured_dB=JAIN["xi_measured"],
                          G_implied=jain_G,
                          seed_uW_for_200uW_output=200.0 / jain_G),
        sensitivity=sens,
        gain_ratio_jain_over_sim=jain_G / SIM["G_measured"])


# ==========================================================================
# Section 4 -- pump-power scans (saturation of the gain)
# ==========================================================================
def _pool(n):
    return ThreadPoolExecutor(max_workers=min(n, (os.cpu_count() or 4)))


def sec_power_scan(quick=False):
    Ps = np.array([1, 3, 10, 20, 35, 50, 80, 110, 135, 160, 200, 250, 300,
                   400, 500, 600, 800] if not quick else
                  [3, 20, 50, 135, 300, 600]) * 1e-3
    out = {}
    for label, cfg in (("jain_geometry", JAIN), ("sim_geometry", SIM)):
        with _pool(len(Ps)) as ex:
            res = list(ex.map(lambda P: gain_peak(cfg, P_pump=P), Ps))
        rows = []
        for P, r in zip(Ps, res):
            Om = constants.rabi_freq(P, cfg["w_pump"])
            rows.append(dict(P_mW=P * 1e3, I_W_cm2=intensity_w_cm2(P, cfg["w_pump"]),
                             Omega_2pi_MHz=Om / (2 * math.pi) / 1e6,
                             Omega_over_Gamma=Om / constants.GAMMA,
                             G_ss=r["G_ss"], delta_mhz=r["delta_mhz"],
                             xi_law_dB=r["xi_law_db"],
                             no_resonance=r["no_resonance"]))
        # local log-slope dln(G-1)/dlnP
        good = [r for r in rows if r["G_ss"] > 1.0005]
        for i, r in enumerate(good):
            if i == 0:
                r["dlnG1_dlnP"] = None
                continue
            p, q = good[i - 1], r
            r["dlnG1_dlnP"] = float(
                (math.log(q["G_ss"] - 1) - math.log(p["G_ss"] - 1))
                / (math.log(q["P_mW"]) - math.log(p["P_mW"])))
        out[label] = rows
    # the headline chain: 600 mW @ 530 um  ->  135 mW @ 300 um
    j = [r for r in out["jain_geometry"] if 80.0 <= r["P_mW"] <= 250.0]
    slopes = [r["dlnG1_dlnP"] for r in j if r.get("dlnG1_dlnP") is not None]
    slope = float(np.mean(slopes)) if slopes else float("nan")
    I_ratio = (intensity_w_cm2(JAIN["P_pump"], JAIN["w_pump"])
               / intensity_w_cm2(SIM["P_pump"], SIM["w_pump"]))
    G0 = SIM["G_measured"]
    G1 = 1.0 + (G0 - 1.0) * I_ratio ** slope
    out["transfer"] = dict(
        note="measured Sim gain propagated to the Jain intensity using only the "
             "engine's local logarithmic sensitivity (no absolute gain used)",
        slope_dlnG1_dlnP_80_250mW=slope,
        intensity_ratio=I_ratio,
        G_sim_measured=G0, G_jain_predicted=G1,
        xi_sim_dB=xi_from_gain(G0, ETA_EXT),
        xi_jain_predicted_dB=xi_from_gain(G1, ETA_EXT),
        penalty_dB=xi_from_gain(G1, ETA_EXT) - xi_from_gain(G0, ETA_EXT))
    return out


# ==========================================================================
# Section 5 -- (P, T) map at the Jain waist
# ==========================================================================
def sec_pt_map(quick=False):
    Ps = np.array([20, 35, 60, 100, 135, 185, 250, 400, 600] if not quick
                  else [35, 135, 600]) * 1e-3
    Ts = (np.arange(88.0, 128.0, 3.0) if not quick
          else np.array([95.0, 105.0, 115.0])) + 273.15
    jobs = [(iP, iT, P, T) for iT, T in enumerate(Ts) for iP, P in enumerate(Ps)]
    G = np.zeros((Ps.size, Ts.size))
    with _pool(len(jobs)) as ex:
        res = list(ex.map(lambda j: gain_peak(JAIN, P_pump=j[2], T=j[3]), jobs))
    for (iP, iT, _, _), r in zip(jobs, res):
        G[iP, iT] = r["G_ss"]
    XI = np.array([[xi_from_gain(g, ETA_EXT) for g in row] for row in G])
    # minimum temperature reaching a target squeezing, per pump power
    targets = {}
    for xi_t in (-6.0, -7.0, -7.2, -8.0):
        need = gain_for_xi(xi_t, ETA_EXT)
        tmin = []
        for iP in range(Ps.size):
            idx = np.where(G[iP, :] >= need)[0]
            tmin.append(float(Ts[idx[0]] - 273.15) if idx.size else None)
        targets[f"{xi_t:.1f}"] = dict(G_required=need,
                                      T_min_C=tmin, P_mW=list(Ps * 1e3))
    # Where does more pump stop helping?  At high density the AC-Stark/power-
    # broadening penalty overtakes the coupling gain and dG/dP turns negative --
    # the mean-field half of Jain's "powers outside 135-250 mW reduce gain or
    # increase noise".  (The noise half is outside this layer.)
    rolloff = []
    for iT, T in enumerate(Ts):
        iP = int(np.argmax(G[:, iT]))
        rolloff.append(dict(T_C=float(T - 273.15), P_best_mW=float(Ps[iP] * 1e3),
                            G_best=float(G[iP, iT]),
                            G_at_max_P=float(G[-1, iT]),
                            interior=bool(0 < iP < Ps.size - 1)))
    return dict(P_mW=list(Ps * 1e3), T_C=list(Ts - 273.15),
                G_ss=G.tolist(), xi_law_dB=XI.tolist(), targets=targets,
                rolloff=rolloff)


# ==========================================================================
# Section 6 -- waist scan, transit-time cost, diffraction limit, depletion
# ==========================================================================
def sec_waist(quick=False):
    ws = np.array([100, 130, 160, 200, 250, 300, 380, 450, 530, 650, 800]
                  if not quick else [150, 300, 530]) * 1e-6
    with _pool(len(ws)) as ex:
        res = list(ex.map(lambda w: gain_peak(JAIN, w_pump=w, w_probe=w), ws))
    I_ref = intensity_w_cm2(SIM["P_pump"], SIM["w_pump"])   # 136 W/cm^2
    rows = []
    for w, r in zip(ws, res):
        rows.append(dict(w_um=w * 1e6,
                         I_W_cm2=intensity_w_cm2(JAIN["P_pump"], w),
                         Omega_2pi_MHz=constants.rabi_freq(JAIN["P_pump"], w)
                         / (2 * math.pi) / 1e6,
                         G_ss=r["G_ss"], xi_law_dB=r["xi_law_db"],
                         rayleigh_mm=rayleigh_mm(w),
                         rayleigh_over_L=rayleigh_mm(w) / (JAIN["L"] * 1e3),
                         transit_floor_kHz=transit_floor_khz(JAIN["T"], w),
                         # pump power that reproduces the Sim reference intensity
                         P_for_ref_intensity_mW=0.5 * math.pi * w ** 2
                         * I_ref * 1e4 * 1e3))
    w_min = math.sqrt(LAMBDA_D1 * JAIN["L"] / math.pi)
    return dict(rows=rows, I_reference_W_cm2=I_ref,
                P_pump_mW=JAIN["P_pump"] * 1e3, T_C=JAIN["T"] - 273.15,
                w_confocal_um=w_min * 1e6,
                note="w for which the Rayleigh range equals the cell length is "
                     "sqrt(lambda L / pi); keeping z_R >= 4L needs w >= 2x that.")


def sec_transit():
    """Cost of the smaller beam through the ground-Raman transit floor."""
    original = fwm.collisional_atom
    rows = []
    try:
        for label, w in (("engine default (fixed 100 kHz floor)", None),
                         ("floor scaled to w = 530 um", 530e-6),
                         ("floor scaled to w = 300 um", 300e-6),
                         ("floor scaled to w = 200 um", 200e-6),
                         ("floor scaled to w = 130 um", 130e-6)):
            fwm.collisional_atom = (original if w is None
                                    else transit_scaled_atom(w))
            r = gain_peak(JAIN)
            rows.append(dict(label=label, w_um=None if w is None else w * 1e6,
                             floor_kHz=(constants.GAMMA_GG / (2 * math.pi) / 1e3
                                        if w is None
                                        else transit_floor_khz(JAIN["T"], w)),
                             G_ss=r["G_ss"], xi_law_dB=r["xi_law_db"]))
    finally:
        fwm.collisional_atom = original
    base = next(r for r in rows if r["w_um"] == 530.0)
    for r in rows:
        r["G_rel_to_530um"] = r["G_ss"] / base["G_ss"]
        r["xi_penalty_dB"] = r["xi_law_dB"] - base["xi_law_dB"]
    return dict(rows=rows,
                note="the engine's GAMMA_GG floor is anchored at the Sim beam "
                     "geometry and does not track the waist, so this patch adds "
                     "the transit cost the engine would otherwise not charge")


def sec_depletion():
    """Manley-Rowe photon budget: the engine's only absolute-power dependence."""
    rows = []
    for label, cfg, G in (("Sim 2025 (measured G)", SIM, SIM["G_measured"]),
                          ("Jain 2025 (implied G)", JAIN,
                           gain_for_xi(JAIN["xi_measured"], ETA_EXT))):
        P_cap = 0.5 * cfg["P_pump"]
        drawn = (G - 1.0) * cfg["P_probe"]
        rows.append(dict(label=label, P_pump_mW=cfg["P_pump"] * 1e3,
                         P_seed_uW=cfg["P_probe"] * 1e6, G=G,
                         power_drawn_uW=drawn * 1e6,
                         conversion_fraction=drawn / P_cap,
                         gain_reduction_factor=1.0 / (1.0 + drawn / P_cap)))
    # Jain's own numbers: ~200 uW per bright output beam
    G_j = gain_for_xi(JAIN["xi_measured"], ETA_EXT)
    rows.append(dict(label="Jain 2025 from reported ~200 uW/beam output",
                     P_pump_mW=JAIN["P_pump"] * 1e3,
                     P_seed_uW=200.0 / G_j, G=G_j, power_drawn_uW=200.0,
                     conversion_fraction=200e-6 / (0.5 * JAIN["P_pump"]),
                     gain_reduction_factor=1.0 / (1.0 + 200e-6 / (0.5 * JAIN["P_pump"]))))
    return dict(rows=rows,
                note="G_s_sat = 1 + (G_s-1)/(1 + (G_s-1) P_seed / (P_pump/2)) "
                     "(observables.pump_depletion_saturation)")


# ==========================================================================
# Reporting
# ==========================================================================
def make_figure(res, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from gabes.plot_style import apply_gabes_plot_style, PALETTE

    fig, axes = plt.subplots(2, 2, figsize=(12.4, 9.2))
    (axA, axB), (axC, axD) = axes

    # (a) iso-intensity lines in the (w, P) plane
    w = np.linspace(80.0, 800.0, 300)
    for I0, col, lab in ((95.5, PALETTE["cyan"], "95 W/cm2  (Jain, 135 mW)"),
                         (136.0, PALETTE["teal"], "136 W/cm2  (Sim, 600 mW)"),
                         (177.0, PALETTE["violet"], "177 W/cm2  (Jain window top)")):
        axA.plot(w, 0.5 * math.pi * (w * 1e-6) ** 2 * I0 * 1e4 * 1e3,
                 color=col, label=lab)
    axA.axhspan(100, 200, color=PALETTE["warm"], alpha=0.12, zorder=0)
    for lab, P, ww, m, off in (("Dowran'18", 550, 500, "s", (-58, -16)),
                               ("Sim'25", 600, 530, "o", (10, 4)),
                               ("Jain'25", 135, 300, "D", (10, -4))):
        axA.plot(ww, P, m, color=PALETTE["ink"], ms=7, zorder=5)
        axA.annotate(lab, (ww, P), textcoords="offset points", xytext=off,
                     fontsize=9, color=PALETTE["ink"])
    axA.set_xlabel("pump 1/e^2 waist radius w  [um]")
    axA.set_ylabel("pump power P  [mW]")
    axA.set_yscale("log")
    axA.set_title("(a) iso-intensity: P ~ w^2  (shaded = 100-200 mW)")
    axA.legend(fontsize=8, loc="upper left")

    # (b) gain vs pump power at w = 300 um, several temperatures
    pt = res["pt_map"]
    Ps = np.array(pt["P_mW"]); Ts = np.array(pt["T_C"]); G = np.array(pt["G_ss"])
    for T_want, col in ((95.0, PALETTE["cyan"]), (100.0, PALETTE["teal"]),
                        (106.0, PALETTE["warm"]), (112.0, PALETTE["rose"])):
        iT = int(np.argmin(np.abs(Ts - T_want)))
        axB.plot(Ps, G[:, iT], "o-", ms=3.5, color=col, label=f"T = {Ts[iT]:.0f} C")
    for xi_t, style in ((-6.0, ":"), (-7.2, "--")):
        axB.axhline(gain_for_xi(xi_t, ETA_EXT), ls=style, lw=1.1,
                    color=PALETTE["muted"])
        axB.annotate(f"G for {xi_t:.1f} dB", (Ps[0], gain_for_xi(xi_t, ETA_EXT)),
                     fontsize=8, color=PALETTE["muted"],
                     textcoords="offset points", xytext=(2, 3))
    axB.axvspan(100, 200, color=PALETTE["warm"], alpha=0.12, zorder=0)
    axB.set_xscale("log"); axB.set_yscale("log")
    axB.set_xticks([20, 50, 100, 200, 400, 600])
    axB.get_xaxis().set_major_formatter(
        matplotlib.ticker.FuncFormatter(lambda v, _: f"{v:g}"))
    axB.get_xaxis().set_minor_formatter(matplotlib.ticker.NullFormatter())
    axB.set_xlabel("pump power P  [mW]   (w = 300 um)")
    axB.set_ylabel("small-signal FWM gain  G")
    axB.set_title("(b) gain is saturated in P, steep in T")
    axB.legend(fontsize=8)

    # (c) the twin-beam noise law
    Gg = np.logspace(0.05, 2.2, 300)
    axC.plot(Gg, [xi_from_gain(g, ETA_EXT) for g in Gg], color=PALETTE["cyan"],
             label=f"eta = {ETA_EXT:.4f}")
    axC.plot(Gg, [xi_from_gain(g, 1.0) for g in Gg], color=PALETTE["muted"],
             ls=":", label="eta = 1 (lossless)")
    axC.axhline(10 * math.log10(1 - ETA_EXT), color=PALETTE["border"], lw=1)
    axC.plot(SIM["G_measured"], SIM["xi_measured"], "o", color=PALETTE["ink"], ms=7)
    axC.annotate("Sim'25 measured\n(G=13.9, -7.8 dB)", (SIM["G_measured"],
                 SIM["xi_measured"]), fontsize=8, textcoords="offset points",
                 xytext=(8, 6), color=PALETTE["ink"])
    Gj = res["noise_law"]["jain_implied"]["G_implied"]
    axC.plot(Gj, JAIN["xi_measured"], "D", color=PALETTE["rose"], ms=7)
    axC.annotate(f"Jain'25 measured\n(-7.2 dB => G={Gj:.1f})", (Gj,
                 JAIN["xi_measured"]), fontsize=8, textcoords="offset points",
                 xytext=(-6, -28), color=PALETTE["rose"])
    axC.set_xscale("log")
    axC.set_xlabel("parametric gain  G")
    axC.set_ylabel("intensity-difference squeezing  xi  [dB]")
    axC.set_title("(c) xi = 10log10[eta/(2G-1) + 1-eta]")
    axC.legend(fontsize=8, loc="lower left")

    # (d) waist scan at fixed 135 mW
    wr = res["waist"]["rows"]
    ww = np.array([r["w_um"] for r in wr])
    axD.plot(ww, [r["G_ss"] for r in wr], "o-", ms=4, color=PALETTE["cyan"])
    axD.axvline(300, color=PALETTE["ink"], lw=1, ls="--")
    axD.annotate("Jain w = 300 um", (300, max(r["G_ss"] for r in wr)),
                 fontsize=8, rotation=90, va="top", textcoords="offset points",
                 xytext=(-12, -4), color=PALETTE["ink"])
    axD.axvline(2 * res["waist"]["w_confocal_um"], color=PALETTE["rose"], lw=1, ls=":")
    axD.annotate("z_R = 4L", (2 * res["waist"]["w_confocal_um"],
                 max(r["G_ss"] for r in wr)), fontsize=8, rotation=90, va="top",
                 textcoords="offset points", xytext=(-12, -4), color=PALETTE["rose"])
    axD.set_yscale("log")
    axD.set_xlabel("beam 1/e^2 waist radius w  [um]   (P = 135 mW, T = 99 C)")
    axD.set_ylabel("small-signal FWM gain  G")
    axD.set_title("(d) focusing buys gain at constant power")

    apply_gabes_plot_style(fig)
    fig.suptitle("85Rb D1 double-Lambda FWM: low-pump gain/ideal-law diagnostic",
                 fontsize=13)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(path, dpi=150)
    print(f"figure -> {path}")


def write_markdown(res, path):
    L = []
    a = L.append
    pub = res["published"]
    a("# 저출력 펌프(100–200 mW) FWM gain/ideal-law 진단\n")
    a("자동 생성 파일 — `analysis/squeezing/low_pump_power/scan_low_pump_power.py`\n")
    a("**Claim gate — `MEAN_FIELD_IDEAL_LAW_DIAGNOSTIC`.** GABES의 절대 이득은 "
      "예측값이 아니며, 아래 ξ 값은 측정 이득 또는 mean-field gain에 조건부 "
      "ideal twin-beam 법칙을 적용한 진단이다. microscopic Langevin "
      "diffusion/covariance와 동일 조건의 측정 SQL이 없으므로 물리적 squeezing, "
      "대역폭 또는 SQL 상하를 독립적으로 예측하지 않는다.\n")
    a(f"검출 효율 η = {ETA_EXT:.4f} (v2 문서 η_ext), 검출 바닥 "
      f"10·log₁₀(1−η) = {10*math.log10(1-ETA_EXT):.3f} dB. "
      f"평균장 이득은 GABES legacy fidelity의 `G_s_smallsignal` (δ 최대점), "
      f"잡음은 폐형식 twin-beam 법칙.\n")

    a("\n## 1. 발표된 동작점은 모두 같은 펌프 세기에 있다\n")
    a("| 동작점 | P [mW] | w [µm] | I [W/cm²] | Ω/2π [MHz] | Ω/Γ | T [°C] | "
      "측정 ξ [dB] | w=300 µm 환산 P [mW] |")
    a("|---|---:|---:|---:|---:|---:|---:|---:|---:|")
    def fmt(x, spec_="{:.1f}"):
        return "—" if x is None else spec_.format(x)
    for r in pub["rows"]:
        a(f"| {r['label']} | {r['P_mW']:.0f} | {r['w_um']:.0f} | "
          f"{r['I_W_cm2']:.1f} | {r['Omega_2pi_MHz']:.0f} | "
          f"{r['Omega_over_Gamma']:.0f} | {fmt(r['T_C'], '{:.0f}')} | "
          f"{fmt(r['xi_measured_dB'])} | {r['P_equiv_at_300um_mW']:.0f} |")
    a(f"\n세 실험의 펌프 세기 범위: **{pub['I_span_W_cm2'][0]:.0f}–"
      f"{pub['I_span_W_cm2'][1]:.0f} W/cm²**, 펌프 라비 "
      f"**{pub['Omega_span_2pi_MHz'][0]:.0f}–{pub['Omega_span_2pi_MHz'][1]:.0f} MHz**. "
      f"Jain 135 mW의 Ω는 Sim 600 mW의 "
      f"{pub['Omega_ratio_jain_over_sim']*100:.0f}%에 불과하다.\n")

    inv = res["invariance"]
    a("\n## 2. 펌프 파워는 이론의 변수가 아니다 — P/w² 불변성\n")
    a("| fidelity | w [µm] | P [mW] | Ω/2π [MHz] | small-signal 비트동일 | "
      "small-signal 편차 | 포화이득 편차 |")
    a("|---|---:|---:|---:|:--:|---:|---:|")
    for fid in ("legacy", "ultra"):
        for r in inv[fid]:
            a(f"| {fid} | {r['w_um']:.0f} | {r['P_mW']:.1f} | "
              f"{r['Omega_2pi_MHz']:.1f} | "
              f"{'예' if r['smallsignal_bitwise_identical'] else '아니오'} | "
              f"{r['smallsignal_peak_rel_dev']*100:+.2f}% | "
              f"{r['saturated_peak_rel_dev']*100:+.1f}% |")
    ac = inv["angle_compensated"]
    a("\nlegacy에서 small-signal 이득 배열은 **비트 단위로 동일**하다 — 즉 "
      "엔진에서 펌프 파워는 독립 변수가 아니다. legacy의 포화이득 편차는 전부 "
      "Manley–Rowe 광자 예산(§6)이고, ultra의 small-signal 편차(≈1%)는 교차 "
      "가우시안 겹침이다. ultra 포화이득 편차가 큰 것은 엔진의 (비예측적) 최대 "
      "이득점 ~5×10³에서 읽었기 때문이며, 실제 동작 이득에서의 고갈 보정은 "
      "§6 표의 0.1% 수준이다.")
    a(f"\n교차각을 허리와 같은 비율로 줄이면 ({SIM['angle']:.2f}° → "
      f"{ac['angle_deg']:.3f}°) 편차는 {ac['smallsignal_peak_rel_dev']*100:+.1f}%로 "
      "오히려 커진다 — θ는 Ultra에서 겹침뿐 아니라 Δk_z 위상정합도 함께 "
      "정하므로 순수한 겹침 보상자가 아니다. 겹침 자체의 크기는 아래 폐형식 "
      "walk-off 표로 읽는 것이 맞다.")
    a("\n| 기하 | 셀 끝 walk-off [µm] | 셀 끝 겹침 |")
    a("|---|---:|---:|")
    for r in inv["walk_off"]:
        a(f"| {r['label']} (w={r['w_um']:.0f} µm, θ={r['angle_deg']:.2f}°, "
          f"L={r['L_mm']:.1f} mm) | {r['end_separation_um']:.1f} | "
          f"{r['end_overlap']:.4f} |")

    nl = res["noise_law"]
    a("\n## 3. 스퀴징이 요구하는 이득은 작다\n")
    a("| 목표 ξ [dB] | 필요 G (η=0.8694) | η=0.90 | η=0.95 |")
    a("|---:|---:|---:|---:|")
    for r in nl["gain_knee"]:
        f = lambda x: "∞" if not math.isfinite(x) else f"{x:.1f}"
        a(f"| {r['xi_target_dB']:.1f} | {f(r['G_required_eta0p8694'])} | "
          f"{f(r['G_required_eta0p90'])} | {f(r['G_required_eta0p95'])} |")
    sc = nl["sim_check"]
    ji = nl["jain_implied"]
    a(f"\n- Sim 2025: 측정 이득 G = {sc['G_measured']:.3f} → 법칙이 주는 ξ = "
      f"{sc['xi_law_dB']:.2f} dB, 측정 {sc['xi_measured_dB']:.1f} dB "
      f"(차이 {sc['residual_dB']:+.2f} dB).")
    a(f"- Jain 2025: 측정 ξ = {ji['xi_measured_dB']:.1f} dB ⇒ 필요 이득 "
      f"G = {ji['G_implied']:.2f}. 보고된 출력 ~200 µW/beam이면 시드 "
      f"≈ {ji['seed_uW_for_200uW_output']:.0f} µW.")
    a(f"- 두 실험의 이득 비 = {nl['gain_ratio_jain_over_sim']:.3f}. "
      f"즉 펌프 파워 4.4배 차이가 이득 1.8배, 스퀴징 0.6 dB 차이로 끝난다.")
    a("\n| G | ξ [dB] | G를 2배 했을 때 얻는 dB |")
    a("|---:|---:|---:|")
    for r in nl["sensitivity"]:
        a(f"| {r['G']:.1f} | {r['xi_dB']:.2f} | {r['d_xi_per_octave_dB']:+.2f} |")

    ps = res["power_scan"]
    a("\n## 4. 이득은 펌프 파워에 대해 깊게 포화되어 있다\n")
    a(f"Jain 기하 (w=300 µm, T={JAIN['T']-273.15:.0f} °C, L={JAIN['L']*1e3:.0f} mm, "
      f"Δ=+{JAIN['D_GHz']:.1f} GHz, θ={JAIN['angle']:.2f}°):\n")
    a("| P [mW] | I [W/cm²] | Ω/Γ | G | δ* [MHz] | ξ(법칙) [dB] | dln(G−1)/dlnP |")
    a("|---:|---:|---:|---:|---:|---:|---:|")
    for r in ps["jain_geometry"]:
        s = r.get("dlnG1_dlnP")
        flag = " *" if r.get("no_resonance") else ""
        a(f"| {r['P_mW']:.0f}{flag} | {r['I_W_cm2']:.1f} | "
          f"{r['Omega_over_Gamma']:.0f} | {r['G_ss']:.3f} | "
          f"{r['delta_mhz']:+.1f} | {r['xi_law_dB']:.2f} | "
          + ("—" if s is None else f"{s:.3f}") + " |")
    a("\n`*` = 스캔 창 안에 FWM 공명이 서지 않은 점(이득 없음). 그 아래로는 "
      "이득이 끊긴다.")
    sim_rows = [r for r in ps["sim_geometry"]
                if r.get("dlnG1_dlnP") is not None and 300.0 <= r["P_mW"] <= 800.0]
    if sim_rows:
        sl = float(np.mean([r["dlnG1_dlnP"] for r in sim_rows]))
        a(f"\n같은 기울기를 Sim 기하(w=530 µm, T=121 °C, L=12.5 mm)에서 재면 "
          f"300–800 mW 구간에서 dln(G−1)/dlnP = {sl:.3f} "
          + ("**로 부호가 음이다** — 즉 그 동작점에서는 펌프를 더 넣을수록 "
             "이득이 오히려 떨어진다(§5의 고밀도 roll-off와 같은 현상). "
             "600 mW라는 값 자체에는 이득상의 이유가 없다."
             if sl < 0 else
             "로, 포화는 기하가 아니라 세기가 만드는 효과다."))
    tr = ps["transfer"]
    a(f"\n**전이 계산** (모델의 절대 이득을 쓰지 않는다): 80–250 mW 구간의 국소 "
      f"기울기 dln(G−1)/dlnP = {tr['slope_dlnG1_dlnP_80_250mW']:.3f}, "
      f"세기 비 {tr['intensity_ratio']:.3f} ⇒ Sim의 측정 이득 "
      f"{tr['G_sim_measured']:.2f} → Jain 조건에서 {tr['G_jain_predicted']:.2f}, "
      f"즉 ξ {tr['xi_sim_dB']:.2f} dB → {tr['xi_jain_predicted_dB']:.2f} dB "
      f"(**대가 {tr['penalty_dB']:+.2f} dB**).")

    pt = res["pt_map"]
    a("\n## 5. (P, T) 지도 — 온도가 지렛대, 파워는 아니다\n")
    a("w = 300 µm에서의 small-signal 이득 G:\n")
    a("| T [°C] | " + " | ".join(f"{p:.0f} mW" for p in pt["P_mW"]) + " |")
    a("|---:|" + "---:|" * len(pt["P_mW"]))
    Gm = np.array(pt["G_ss"])
    for iT, T in enumerate(pt["T_C"]):
        a(f"| {T:.0f} | " + " | ".join(f"{Gm[iP, iT]:.2f}" for iP in
                                       range(len(pt["P_mW"]))) + " |")
    for key, tgt in pt["targets"].items():
        pairs = [f"{p:.0f} mW→" + ("—" if t is None else f"{t:.0f} °C")
                 for p, t in zip(tgt["P_mW"], tgt["T_min_C"])]
        a(f"\n- ξ = {key} dB (G ≥ {tgt['G_required']:.1f}) 최소 온도: "
          + ", ".join(pairs))
    # lever comparison, read straight off the grid
    Pv, Tv = np.array(pt["P_mW"]), np.array(pt["T_C"])
    iP0 = int(np.argmin(np.abs(Pv - 135.0)))
    iT0 = int(np.argmin(np.abs(Tv - 100.0)))
    dT = 15.0
    iT1 = int(np.argmin(np.abs(Tv - (Tv[iT0] + dT))))
    a(f"\n**지렛대 비교** — {Tv[iT0]:.0f} °C, {Pv[iP0]:.0f} mW를 기준으로: "
      f"파워를 {Pv[0]:.0f} → {Pv[-1]:.0f} mW ({Pv[-1]/Pv[0]:.0f}배) 올리면 "
      f"이득은 {Gm[0, iT0]:.2f} → {Gm[-1, iT0]:.2f} (×{Gm[-1, iT0]/Gm[0, iT0]:.1f}), "
      f"온도를 {Tv[iT0]:.0f} → {Tv[iT1]:.0f} °C ({Tv[iT1]-Tv[iT0]:.0f} °C) 올리면 "
      f"{Gm[iP0, iT0]:.2f} → {Gm[iP0, iT1]:.2f} (×{Gm[iP0, iT1]/Gm[iP0, iT0]:.0f}). "
      f"목표 스퀴징에 필요한 최소 온도는 20–400 mW 전 구간에서 사실상 같은 값으로 "
      f"나오며, 이것이 '펌프 파워는 자유 변수가 아니다'의 실무적 형태다.")
    ro = [r for r in pt["rolloff"] if r["interior"]]
    a("\n게다가 밀도가 높아지면 **파워를 더 넣는 것이 이득을 깎는다**: ")
    if ro:
        a("| T [°C] | 최적 P [mW] | 그때 G | 최대 P에서의 G |")
        a("|---:|---:|---:|---:|")
        for r in ro:
            a(f"| {r['T_C']:.0f} | {r['P_best_mW']:.0f} | {r['G_best']:.2f} | "
              f"{r['G_at_max_P']:.2f} |")
        a("\n이것이 Jain이 관측한 \"135–250 mW 창 바깥에서는 이득이 줄거나 "
          "잡음이 는다\"의 평균장 쪽 절반이다 (잡음 쪽 절반은 이 층 바깥이다).")
    else:
        a("(이 격자 범위에서는 내부 최적 P가 잡히지 않았다.)")
    a("\n> 주의: 엔진의 절대 온도 눈금은 보정되어 있지 않다. 여기서 읽을 것은 "
      "**∂lnG/∂T ≫ ∂lnG/∂lnP** 라는 지렛대 비교이지, 특정 목표를 위한 절대 "
      "온도값이 아니다.")

    wa = res["waist"]
    a("\n## 6. 허리를 줄이는 값 — 회절·통과시간·펌프 고갈\n")
    a(f"P = {wa['P_pump_mW']:.0f} mW 고정, T = {wa['T_C']:.0f} °C:\n")
    a("| w [µm] | I [W/cm²] | Ω/2π [MHz] | G | ξ(법칙) [dB] | z_R [mm] | "
      "z_R/L | 통과시간 바닥 [kHz] | Sim 세기 재현에 필요한 P [mW] |")
    a("|---:|---:|---:|---:|---:|---:|---:|---:|---:|")
    for r in wa["rows"]:
        a(f"| {r['w_um']:.0f} | {r['I_W_cm2']:.0f} | {r['Omega_2pi_MHz']:.0f} | "
          f"{r['G_ss']:.2f} | {r['xi_law_dB']:.2f} | {r['rayleigh_mm']:.0f} | "
          f"{r['rayleigh_over_L']:.0f} | {r['transit_floor_kHz']:.0f} | "
          f"{r['P_for_ref_intensity_mW']:.0f} |")
    a(f"\n마지막 열은 Sim 기준 세기 {wa['I_reference_W_cm2']:.0f} W/cm²를 그 "
      f"허리에서 재현하는 데 드는 펌프 파워다. z_R = L이 되는 허리는 "
      f"{wa['w_confocal_um']:.0f} µm이므로 z_R ≥ 4L을 지키려면 "
      f"w ≳ {2*wa['w_confocal_um']:.0f} µm이고, 그 허리에서도 "
      f"{0.5*math.pi*(2*wa['w_confocal_um']*1e-6)**2*wa['I_reference_W_cm2']*1e4*1e3:.0f} mW면 "
      f"Sim과 같은 펌프 세기에 도달한다.")

    tt = res["transit"]
    a("\n통과시간 바닥을 v̄/w로 실제 스케일했을 때의 대가:\n")
    a("| 설정 | 바닥 [kHz] | G | 530 µm 대비 | ξ 대가 [dB] |")
    a("|---|---:|---:|---:|---:|")
    for r in tt["rows"]:
        a(f"| {r['label']} | {r['floor_kHz']:.0f} | {r['G_ss']:.4f} | "
          f"{r['G_rel_to_530um']:.4f} | {r['xi_penalty_dB']:+.3f} |")

    dep = res["depletion"]
    a("\nManley–Rowe 펌프 광자 예산 (엔진의 유일한 절대 파워 의존성):\n")
    a("| 동작점 | P_pump [mW] | 시드 [µW] | G | 뽑아간 파워 [µW] | 변환율 | "
      "이득 감소 |")
    a("|---|---:|---:|---:|---:|---:|---:|")
    for r in dep["rows"]:
        a(f"| {r['label']} | {r['P_pump_mW']:.0f} | {r['P_seed_uW']:.1f} | "
          f"{r['G']:.2f} | {r['power_drawn_uW']:.1f} | "
          f"{r['conversion_fraction']*100:.3f}% | "
          f"{(1-r['gain_reduction_factor'])*100:.3f}% |")

    a("\n## 7. 결론\n")
    a("1. 펌프는 Ω² ∝ P/w²로만 이론에 들어간다. GABES legacy에서 "
      "(600 mW, 530 µm)와 (192 mW, 300 µm)의 small-signal 이득 배열은 "
      "**비트 단위로 동일**하다. 유일한 절대 파워 의존성은 Manley–Rowe "
      "광자 예산이고, 실제 변환율은 0.1% 수준이라 구속되지 않는다.")
    a("2. 이미 발표된 세 실험은 모두 같은 펌프 세기 대역에 있다. "
      "Jain의 135–250 mW 창은 Sim/Dowran의 550–600 mW 점을 라비 주파수에서 "
      "감싼다.")
    a("3. −7 dB에 필요한 이득은 G ≈ 8뿐이고, 이 영역에서 이득은 파워에 대해 "
      "dln(G−1)/dlnP ≈ 0.1로 깊게 포화되어 있다. 두 효과를 합치면 600 mW → "
      "135 mW의 대가는 0.1 dB 미만이다.")
    a("4. 허리를 줄이는 비용(통과시간 확장, 회절, 교차 겹침)은 300 µm에서 "
      "모두 1% 수준이다. 100–200 mW는 물리적 문턱이 아니라 단순한 재초점 "
      "문제다.")
    wmin4 = 2 * res["waist"]["w_confocal_um"]
    a(f"5. 여유도 남아 있다. 12 mm 셀에서 z_R ≥ 4L을 지키는 최소 허리는 "
      f"w ≈ {wmin4:.0f} µm이고, 그 허리에서 Sim과 같은 펌프 세기를 만드는 데 "
      f"필요한 파워는 "
      f"{0.5*math.pi*(wmin4*1e-6)**2*res['waist']['I_reference_W_cm2']*1e4*1e3:.0f} mW다. "
      f"즉 100–200 mW는 하한이 아니라 편안한 동작점이다. 대신 줄어드는 것은 "
      f"밝은 twin-beam의 절대 출력(모드 면적에 비례)과 공간 모드 수이지, "
      f"스퀴징 수준이 아니다.")
    a("\n### 이 문서가 주장하지 않는 것\n")
    a("엔진의 **절대 이득은 예측값이 아니다** (v2 「Model point and experimental "
      "comparison」: Sim 동작점에서 G_p가 +7500%). G = cosh²(qL)이므로 qL의 "
      "2배 오차가 G의 50배 오차가 된다. 위의 모든 헤드라인 수치는 (a) 구조적 "
      "불변성, (b) 로그 민감도, (c) 측정 이득을 입력으로 받는 잡음 법칙 중 "
      "하나이며, 절대 이득 눈금을 필요로 하지 않는다.")

    text = "\n".join(L) + "\n"
    Path(path).write_text(text, encoding="utf-8")
    print(f"report -> {path}")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default=str(Path(__file__).resolve().parent / "generated"))
    ap.add_argument("--quick", action="store_true", help="coarse grids, for smoke tests")
    ap.add_argument("--no-figure", action="store_true")
    args = ap.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    t0 = time.time()
    res = {}
    for name, fn in (("published", lambda: sec_published()),
                     ("noise_law", lambda: sec_noise_law()),
                     ("invariance", lambda: sec_invariance()),
                     ("power_scan", lambda: sec_power_scan(args.quick)),
                     ("pt_map", lambda: sec_pt_map(args.quick)),
                     ("waist", lambda: sec_waist(args.quick)),
                     ("transit", lambda: sec_transit()),
                     ("depletion", lambda: sec_depletion())):
        t = time.time()
        res[name] = fn()
        print(f"  {name:12s} {time.time()-t:7.1f}s")

    res["meta"] = dict(
        eta_ext=ETA_EXT, line_strength_residual=LS_RESIDUAL,
        mean_field_layer="GABES legacy fidelity G_s_smallsignal, peak over delta",
        noise_layer="observables.intensity_difference_squeezing_dB with G_c=G_s-1",
        diagnostic_claim_gate=dict(
            status="MEAN_FIELD_IDEAL_LAW_DIAGNOSTIC",
            physical_squeezing_prediction=False,
            bandwidth_prediction=False,
            sql_crossing_prediction=False,
            reason=("No microscopic Langevin diffusion/covariance or "
                    "same-condition measured SQL is included."),
        ),
        sim=dict(SIM), jain=dict(JAIN), dowran=dict(DOWRAN),
        runtime_s=time.time() - t0)

    jpath = out / "low_pump_power.json"
    jpath.write_text(json.dumps(res, indent=1, ensure_ascii=False,
                                default=float), encoding="utf-8")
    print(f"json   -> {jpath}")
    write_markdown(res, out / "low_pump_power.md")
    if not args.no_figure:
        make_figure(res, out / "low_pump_power.png")

    # console summary (ASCII only)
    nl, ps = res["noise_law"], res["power_scan"]["transfer"]
    print("\n=== HEADLINE ===")
    print(f"published pump intensities  : {res['published']['I_span_W_cm2'][0]:.0f}"
          f"-{res['published']['I_span_W_cm2'][1]:.0f} W/cm2 across 135-600 mW")
    print(f"Jain/Sim pump Rabi ratio    : "
          f"{res['published']['Omega_ratio_jain_over_sim']:.3f}")
    print(f"G needed for -7.2 dB        : {nl['jain_implied']['G_implied']:.2f}"
          f"   (Sim measured G = {nl['sim_check']['G_measured']:.2f})")
    print(f"twin-beam law on Sim        : {nl['sim_check']['xi_law_dB']:.2f} dB "
          f"vs measured {nl['sim_check']['xi_measured_dB']:.1f} dB")
    print(f"gain saturation dlnG/dlnP   : {ps['slope_dlnG1_dlnP_80_250mW']:.3f}")
    print(f"600 mW/530um -> 135 mW/300um: {ps['penalty_dB']:+.2f} dB")
    print(f"total {time.time()-t0:.1f}s")


if __name__ == "__main__":
    main()
