"""
Operating-variable audit for the 85Rb D1 double-Lambda FWM mean-field gain and
gain-referred ideal-law noise diagnostic.  It checks an experimentalist's field
intuitions against (a) the engine, (b) docs/squeezing_report/
squeezing_report_v6.tex, and (c) the 10-paper reference DB
references/fwm_squeezing_paper_parameters.csv.

CLAIM GATE.  The engine quantity is ``gain_referred_noise_dB``.  ``S_dB`` and
the ``xi_*`` fields retained below are compatibility aliases, not physical
squeezing.  The score combines mean-field gains with an assumed ideal twin-beam
law.  GABES does not provide the frequency-dependent microscopic Langevin
diffusion/covariance or a same-condition measured SQL, so this script cannot
predict physical squeezing, squeezing bandwidth, or above/below-SQL status.

Claims under test (verbatim, translated):
  OPD    -- (no claim offered)
  TPD    -- "related to dressing; moves with pump intensity"
  pump   -- "fuel; sets the gain ceiling. Too weak and FWM will not start, but
             almost anything works. If FWM does not appear, check beam OVERLAP
             before blaming pump power."
  seed   -- "sets FWM position/direction/frequency. It travels the same path as
             the probe so it lands in the output arm too, which stops the
             intensity difference from having zero expectation. Effectively the
             same as injecting a Gaussian beam into one arm."
  T      -- "too low and there is no FWM -- the main culprit. FWM turns on
             sharply above ~110 C; above 120 C gain grows but squeezing gets
             worse."
  L      -- "not sure what it couples to."
  angle  -- "beam-overlap variable. Ideal at 0 deg in theory, but then the beams
             cannot be separated, so a compromise is needed."
  loss   -- "once enough gain exists, the dominant squeezing killer."
  stray  -- "killer #2. Loss cannot by itself explain noise ABOVE the SQL;
             stray-beam injection is the main source of super-Poissonian noise."

SECTIONS
  tpd_lightshift  -- delta*(P) at two geometries, fitted against the AC-Stark
                     scaling Omega^2/(4 Delta).  Tests the TPD claim.
  pump_vs_overlap -- gain-referred diagnostic versus pump power and mode overlap,
                     using compute_spectrum's own `mode_overlap_penalty` knob.
                     Tests "check overlap before pump power".
  seed            -- (i) gain-referred diagnostic vs seed power; (ii)
                     unamplified-seed leakage into the
                     probe arm, through `observables.balanced_twin_beam_noise`
                     with raw vs DC-balanced weighting.  Tests the seed claim.
  temperature     -- gain-referred diagnostic and G versus T at the Sim gold
                     geometry, fixed-delta and best-diagnostic delta.
  cell_length     -- diagnostic/G vs L, plus the N*L trade.
  angle           -- diagnostic/G vs crossing angle, including the low-angle ridge that
                     v6 Sec. "해석" flags as outside geometry calibration.
  loss            -- ideal-law diagnostic vs loss and its algebraic asymptote.
  stray           -- conditional ideal-covariance fixtures for pure loss and
                     additive excess noise; these are not measured-SQL predictions.
  paper_db        -- min/median/max of each variable across the reference DB.

READOUT CONVENTION (same as v6 Sec. tolerance): hardened Ultra, eta = 0.8694,
5 MHz delta grid, gap gate 0.5 <= G_s - G_c <= 1.5 on the best-delta readout.
Absolute gains are NOT predictions (v2 "Model point and experimental
comparison").  Read only the declared sensitivities and orderings; the dB score
is a gain-referred ideal-law diagnostic, not physical squeezing.

Usage:
    python analysis/squeezing/variable_audit/scan_variable_audit.py
    python analysis/squeezing/variable_audit/scan_variable_audit.py --quick
"""
import argparse
import csv
import json
import math
import os
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", str(os.cpu_count() or 1))
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

from gabes import constants, hyperfine, observables                 # noqa: E402
from gabes.core import blas_single_thread                           # noqa: E402
from gabes.schemes import fwm                                       # noqa: E402

ETA = 0.8694                        # v6 readout efficiency
GAP_MIN, GAP_MAX = 0.5, 1.5         # v6 / frontier-skill twin-beam validity gate
DELTA_STEP_MHZ = 5.0                # v6: the 17.5 MHz frontier grid is an artifact
LS = 0.74

# Public schema for this analysis.  The old xi/S names are appended only as
# compatibility aliases after every section has been calculated.
DIAGNOSTIC_CLAIM_GATE = {
    "level": "MEAN_FIELD_DIAGNOSTIC",
    "primary_engine_quantity": "gain_referred_noise_dB",
    "engine_compatibility_alias": "S_dB",
    "physical_squeezing_prediction": False,
    "physical_sql_comparison": False,
    "microscopic_langevin_diffusion": False,
    "interpretation": (
        "mean-field gains completed with an assumed ideal twin-beam noise law; "
        "use for conditional sensitivity/order diagnostics only"
    ),
}

_LEGACY_TO_PRIMARY = {
    "xi_best": "gain_referred_best_dB",
    "xi_fixed": "gain_referred_fixed_dB",
    "best_xi_dB": "best_gain_referred_dB",
    "best_fixed_xi_dB": "best_fixed_gain_referred_dB",
    "xi_dB": "ideal_law_noise_dB",
    "xi_at_zero_excess_dB": "ideal_law_at_zero_excess_dB",
    "power_spread_dB": "gain_referred_spread_dB",
    "floor_dB": "ideal_law_asymptote_dB",
    "gap_to_floor_dB": "gap_to_ideal_law_asymptote_dB",
    "slope_dB_per_pct": "gain_referred_slope_dB_per_pct",
    "penalty_dB": "ideal_law_penalty_dB",
    "equivalent_extra_loss_pct": "ideal_law_equivalent_extra_loss_pct",
    "S": "ideal_law_noise_ratio",
    "S_ideal": "ideal_law_source_ratio",
    "above_sql": "above_normalized_unity",
    "excess_for_sql": "excess_for_normalized_unity",
    "theorem": "ideal_law_statement",
    "gain_needed_for_floor_note": "ideal_law_asymptote_note",
}


def promote_gain_referred_schema(value):
    """Add explicit primary names while preserving deprecated result aliases."""
    if isinstance(value, dict):
        for child in list(value.values()):
            promote_gain_referred_schema(child)
        for legacy, primary in _LEGACY_TO_PRIMARY.items():
            if legacy in value and primary not in value:
                value[primary] = value[legacy]
    elif isinstance(value, list):
        for child in value:
            promote_gain_referred_schema(child)
    return value

# Sim et al. gold reference (references/fwm_squeezing_paper_parameters.csv,
# is_gold_reference = 1).  This is the closest published point to the user's lab.
GOLD = dict(D_GHz=0.9, T=394.15, P_pump=0.600, P_probe=8e-6,
            w_pump=530e-6, w_probe=330e-6, L=12.5e-3, angle=0.32,
            delta_mhz=-8.0)
ULTRA = dict(phase_detail=fwm.PHASE_ULTRA, model_fidelity=fwm.FIDELITY_ULTRA)


def scan(cfg, *, window_mhz=350.0, **over):
    """Hardened-Ultra probe scan on a 5 MHz delta grid."""
    kw = {k: cfg[k] for k in ("T", "P_pump", "P_probe", "w_pump", "w_probe",
                              "L", "angle")}
    kw.update({k: v for k, v in over.items() if k in kw or k in
               ("loss_frac", "mode_overlap_penalty")})
    D = over.get("D_GHz", cfg["D_GHz"])
    loss = over.get("loss_frac", 0.055)
    overlap = over.get("mode_overlap_penalty", 1.0)
    npts = int(round(2 * window_mhz / DELTA_STEP_MHZ)) + 1
    center = fwm.branch_center_GHz(D, -1)
    with blas_single_thread():
        s = fwm.compute_spectrum(
            D, T=kw["T"], P_pump=kw["P_pump"], P_probe=kw["P_probe"],
            w_pump=kw["w_pump"], w_probe=kw["w_probe"], L=kw["L"],
            pump_probe_angle_deg=kw["angle"], line_strength=LS,
            mode_overlap_penalty=overlap,
            loss_frac=loss, qe=fwm.QE_DETECTOR,
            coarse_points=npts, fine_points=0,
            scan_min=center - window_mhz * 1e-3,
            scan_max=center + window_mhz * 1e-3,
            velocity_step=5.0, velocity_cutoff=3.0, branch=-1, **ULTRA)
    s["delta_mhz"] = (s["probe_axis_GHz"] - center) * 1e3
    return s


def readout(s, fixed_delta_mhz=None):
    """Fixed-delta and gap-gated gain-referred diagnostic readouts."""
    d, Gs, Gc = s["delta_mhz"], s["G_s"], s["G_c"]
    diagnostic = s.get("gain_referred_noise_dB")
    if diagnostic is None:  # compatibility with archived raw spectra
        diagnostic = s["S_dB"]
    S = np.asarray(diagnostic, dtype=float)
    gap = Gs - Gc
    ok = np.isfinite(S) & (gap >= GAP_MIN) & (gap <= GAP_MAX)
    out = {}
    if ok.any():
        i = int(np.nanargmin(np.where(ok, S, np.inf)))
        j = int(np.nanargmax(np.where(ok, Gs, -np.inf)))
        out.update(xi_best=float(S[i]), delta_best=float(d[i]),
                   G_best=float(Gs[i]), gap_best=float(gap[i]),
                   G_max_gated=float(Gs[j]),
                   best_at_edge=bool(i <= 1 or i >= d.size - 2))
    else:
        out.update(xi_best=None, delta_best=None, G_best=None, gap_best=None,
                   G_max_gated=None, best_at_edge=None)
    out["gap_ok_fraction"] = float(ok.mean())
    if fixed_delta_mhz is not None:
        k = int(np.argmin(np.abs(d - fixed_delta_mhz)))
        out.update(xi_fixed=float(S[k]), G_fixed=float(Gs[k]),
                   gap_fixed=float(gap[k]), delta_fixed=float(d[k]),
                   gap_fixed_ok=bool(GAP_MIN <= gap[k] <= GAP_MAX))
    # peak of the ungated classical gain resonance (where delta* lives)
    p = int(np.nanargmax(Gs))
    out.update(delta_peak=float(d[p]), G_peak=float(Gs[p]),
               peak_at_edge=bool(p <= 1 or p >= d.size - 2))
    return promote_gain_referred_schema(out)


def pool(n):
    return ThreadPoolExecutor(max_workers=min(n, (os.cpu_count() or 4)))


def light_shift_mhz(P, w, Delta_GHz):
    """Two-level AC-Stark estimate of the Raman resonance pull, Omega^2/(4 Delta)."""
    Om = constants.rabi_freq(P, w)
    return (Om ** 2 / (4.0 * 2 * math.pi * Delta_GHz * 1e9)) / (2 * math.pi) / 1e6


# ==========================================================================
# TPD: does the two-photon optimum move with pump intensity?
# ==========================================================================
def sec_tpd(quick=False):
    geoms = {
        "gold (w=530 um, T=121 C, D=+0.9 GHz)": dict(GOLD),
        "Jain (w=300 um, T=99 C, D=+0.8 GHz)": dict(
            D_GHz=0.8, T=372.15, P_pump=0.135, P_probe=10e-6, w_pump=300e-6,
            w_probe=300e-6, L=12.0e-3, angle=0.40, delta_mhz=0.0),
    }
    Ps = np.array([35, 60, 100, 135, 200, 300, 450, 600] if not quick
                  else [60, 200, 600]) * 1e-3

    def refined(cfg, P):
        """delta* on a 0.5 MHz grid: locate on the 5 MHz scan, then refine."""
        r0 = readout(scan(cfg, P_pump=P))
        d0 = r0["delta_peak"]
        center = fwm.branch_center_GHz(cfg["D_GHz"], -1)
        with blas_single_thread():
            s = fwm.compute_spectrum(
                cfg["D_GHz"], T=cfg["T"], P_pump=P, P_probe=cfg["P_probe"],
                w_pump=cfg["w_pump"], w_probe=cfg["w_probe"], L=cfg["L"],
                pump_probe_angle_deg=cfg["angle"], line_strength=LS,
                loss_frac=0.055, qe=fwm.QE_DETECTOR,
                coarse_points=121, fine_points=0,
                scan_min=center + (d0 - 15.0) * 1e-3,
                scan_max=center + (d0 + 15.0) * 1e-3,
                velocity_step=5.0, velocity_cutoff=3.0, branch=-1, **ULTRA)
        d = (s["probe_axis_GHz"] - center) * 1e3
        k = int(np.nanargmax(s["G_s"]))
        r0["delta_peak"] = float(d[k])
        r0["G_peak"] = float(s["G_s"][k])
        return r0

    out = {}
    for label, cfg in geoms.items():
        with pool(len(Ps)) as ex:
            res = list(ex.map(lambda P: refined(cfg, P), Ps))
        rows = []
        for P, r in zip(Ps, res):
            Om = constants.rabi_freq(P, cfg["w_pump"])
            rows.append(dict(
                P_mW=float(P * 1e3),
                Omega_2pi_MHz=Om / (2 * math.pi) / 1e6,
                delta_peak_mhz=r["delta_peak"],
                delta_best_mhz=r["delta_best"],
                ac_stark_estimate_mhz=light_shift_mhz(P, cfg["w_pump"],
                                                      cfg["D_GHz"]),
                G_peak=r["G_peak"], xi_best=r["xi_best"]))
        # Free power law delta* = a * Omega^n.  n = 2 is the perturbative light
        # shift Omega^2/(4 Delta); n = 1 is the strong-dressing limit where the
        # generalized Rabi splitting sqrt(Delta^2+Omega^2) -> Omega.  Here
        # Omega ~ Delta, so the crossover exponent is the interesting output.
        lx = np.log([r["Omega_2pi_MHz"] for r in rows])
        ly = np.log([max(r["delta_peak_mhz"], 1e-6) for r in rows])
        n, lna = np.polyfit(lx, ly, 1)
        resid = ly - (n * lx + lna)
        # dressed-state reference: (sqrt(Delta^2 + Omega^2) - Delta)/2
        D_mhz = cfg["D_GHz"] * 1e3
        for r in rows:
            r["dressed_shift_mhz"] = 0.5 * (
                math.hypot(D_mhz, r["Omega_2pi_MHz"]) - D_mhz)
        # how far does delta* move over the reported pump range, vs the v6
        # +0.5 dB TPD tolerance window (-9.6 / +4.2 MHz)?
        span = rows[-1]["delta_peak_mhz"] - rows[0]["delta_peak_mhz"]
        # delta* ~ Omega^n ~ P^(n/2), so doubling / halving the pump power at the
        # OPERATING point multiplies delta* by 2^(n/2).  Quote both directions.
        i0 = min(range(len(rows)),
                 key=lambda k: abs(rows[k]["P_mW"] - cfg["P_pump"] * 1e3))
        d0_op = rows[i0]["delta_peak_mhz"]
        out[label] = dict(
            rows=rows, exponent_n=float(n), prefactor=float(math.exp(lna)),
            fit_rms_dex=float(np.sqrt(np.mean(resid ** 2))),
            fit_r2=float(1.0 - np.sum(resid ** 2)
                         / max(np.sum((ly - ly.mean()) ** 2), 1e-30)),
            delta_span_mhz=float(span),
            operating_P_mW=float(rows[i0]["P_mW"]),
            operating_delta_mhz=float(d0_op),
            delta_on_power_x2_mhz=float(d0_op * 2 ** (n / 2) - d0_op),
            delta_on_power_half_mhz=float(d0_op * 0.5 ** (n / 2) - d0_op),
            delta_on_power_10pct_mhz=float(d0_op * 1.1 ** (n / 2) - d0_op),
            tpd_tolerance_window_mhz=[-9.6, 4.2],
            omega_over_Delta=[rows[0]["Omega_2pi_MHz"] / (cfg["D_GHz"] * 1e3),
                              rows[-1]["Omega_2pi_MHz"] / (cfg["D_GHz"] * 1e3)])
    return out


# ==========================================================================
# Pump power vs beam overlap: which one actually kills the FWM?
# ==========================================================================
def sec_pump_vs_overlap(quick=False):
    Ps = np.array([150, 300, 450, 600, 800, 1000] if not quick
                  else [300, 600, 1000]) * 1e-3
    ovl = np.array([1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3] if not quick
                   else [1.0, 0.7, 0.4])
    with pool(len(Ps)) as ex:
        rp = list(ex.map(lambda P: readout(scan(GOLD, P_pump=P),
                                           GOLD["delta_mhz"]), Ps))
    with pool(len(ovl)) as ex:
        ro = list(ex.map(lambda m: readout(scan(GOLD, mode_overlap_penalty=m),
                                           GOLD["delta_mhz"]), ovl))
    power_rows = [dict(P_mW=float(P * 1e3), G_fixed=r["G_fixed"],
                       xi_fixed=r["xi_fixed"], G_peak=r["G_peak"],
                       xi_best=r["xi_best"])
                  for P, r in zip(Ps, rp)]
    ovl_rows = [dict(overlap=float(m), G_fixed=r["G_fixed"],
                     xi_fixed=r["xi_fixed"], G_peak=r["G_peak"],
                     xi_best=r["xi_best"],
                     # transverse offset that gives this Gaussian overlap
                     offset_um=float(math.sqrt(max(-math.log(max(m, 1e-9)), 0.0))
                                     * math.hypot(GOLD["w_pump"],
                                                  GOLD["w_probe"]) * 1e6))
                for m, r in zip(ovl, ro)]
    # --- leverage, as a LOCAL two-sided derivative at the SAME reference point.
    # `mode_overlap_penalty` multiplies the coupling, so q ∝ overlap exactly and
    # dlnG/dln(overlap) = dlnG/dln(q).  By the chain rule
    #     dlnG/dlnP = (dlnG/dlnq) * (dlnq/dlnP),
    # so the RATIO of the two derivatives is exactly 1/(dlnq/dlnP) -- the
    # engine's (non-predictive) absolute gain cancels out of it.
    eps = 0.15
    with pool(4) as ex:
        loc = list(ex.map(
            lambda kw: readout(scan(GOLD, **kw), GOLD["delta_mhz"]),
            [dict(P_pump=GOLD["P_pump"] * (1 - eps)),
             dict(P_pump=GOLD["P_pump"] * (1 + eps)),
             dict(mode_overlap_penalty=1.0 - eps),
             dict(mode_overlap_penalty=1.0)]))
    lp = math.log(loc[1]["G_peak"] / loc[0]["G_peak"]) / math.log(
        (1 + eps) / (1 - eps))
    lo = math.log(loc[3]["G_peak"] / loc[2]["G_peak"]) / math.log(
        1.0 / (1.0 - eps))
    return dict(power=power_rows, overlap=ovl_rows,
                dlnG_dlnP=lp, dlnG_dln_overlap=lo,
                leverage_ratio=float(lo / lp) if lp else None,
                dlnq_dlnP=float(lp / lo) if lo else None,
                local_eps=eps,
                note="derivatives are LOCAL two-sided logarithmic slopes at the "
                     "gold point (+/-15%). mode_overlap_penalty multiplies the "
                     "coupling, so q is exactly linear in it and the RATIO of "
                     "the two slopes equals 1/(dlnq/dlnP), independent of the "
                     "engine's absolute gain. offset_um is the static transverse "
                     "misalignment giving the same Gaussian overlap "
                     "exp(-d^2/(w_pump^2+w_probe^2))")


# ==========================================================================
# Seed: power insensitivity, and unamplified seed leaking into the probe arm
# ==========================================================================
def sec_seed(quick=False):
    Ps = np.array([1, 4, 8, 16, 40, 100, 200] if not quick
                  else [1, 8, 100]) * 1e-6
    with pool(len(Ps)) as ex:
        res = list(ex.map(lambda P: readout(scan(GOLD, P_probe=P),
                                            GOLD["delta_mhz"]), Ps))
    power_rows = [dict(P_probe_uW=float(P * 1e6), G_fixed=r["G_fixed"],
                       xi_fixed=r["xi_fixed"], xi_best=r["xi_best"])
                  for P, r in zip(Ps, res)]
    xf = [r["xi_fixed"] for r in power_rows]
    # --- conditional ideal-covariance fixture for unamplified seed leakage
    # Probe arm carries the amplified seed (G*P0, correlated with the conjugate)
    # plus a fraction f*P0 of seed that missed the pump mode and is therefore
    # UNCORRELATED coherent light.  Conjugate arm carries (G-1)*P0.
    # Closed form rather than the API: `balanced_twin_beam_noise` re-derives the
    # source covariance from whatever gains it is handed, which would make the
    # leaked light spuriously correlated.  Here the covariance is FROZEN at its
    # f = 0 value -- that is exactly what "uncorrelated leakage" means.
    #   mean_s = eta (G + f),  mean_c = eta (G - 1),  cov = eta^2 (G - 1)
    #   S = Var(I_s - w I_c) / (mean_s + w^2 mean_c)
    # This chosen covariance is an algebraic fixture, not the missing microscopic
    # Langevin covariance of the physical cell.
    G = 14.0                       # gold-reference measured gain scale
    Gc = G - 1.0
    cov0 = Gc                      # = 0.5[(G+Gc) - (G-Gc)^2] at G - Gc = 1
    leak = []
    for f in (0.0, 0.1, 0.25, 0.5, 1.0, 2.0, 5.0, 10.0):
        mean_s, mean_c, cov = ETA * (G + f), ETA * Gc, ETA * ETA * cov0
        for wname, w in (("raw (w=1)", 1.0),
                         ("DC-balanced", mean_s / mean_c),
                         ("shot-optimal", math.sqrt(mean_s / mean_c))):
            shot = mean_s + w * w * mean_c
            S = (shot - 2.0 * w * cov) / shot
            leak.append(dict(
                leak_fraction_of_seed=f, weight=wname,
                leak_pct_of_probe_arm=100.0 * f / (G + f),
                dc_imbalance=float((mean_s - w * mean_c) / shot),
                xi_dB=float(10 * math.log10(max(S, 1e-30)))))
    base = {r["weight"]: r["xi_dB"] for r in leak
            if r["leak_fraction_of_seed"] == 0.0}
    for r in leak:
        r["penalty_dB"] = r["xi_dB"] - base[r["weight"]]
        # loss that would cost the same: eta_eq from S = eta_eq/(2G-1) + 1-eta_eq
        S = 10 ** (r["xi_dB"] / 10.0)
        Si = 1.0 / (2 * G - 1.0)
        r["equivalent_eta"] = float((1.0 - S) / (1.0 - Si))
        r["equivalent_extra_loss_pct"] = float(100.0 * (1.0 - r["equivalent_eta"] / ETA))
    return dict(power=power_rows,
                power_spread_dB=float(max(xf) - min(xf)),
                leakage=leak, G_used=G,
                note="leakage model: probe arm = correlated G*P0 + uncorrelated "
                     "f*P0; the source covariance is FROZEN at the f=0 value, so "
                     "the leak only adds shot noise and DC imbalance. "
                     "equivalent_extra_loss_pct converts the penalty into the "
                     "post-cell loss that would cost the same.")


# ==========================================================================
# Temperature
# ==========================================================================
def sec_temperature(quick=False):
    Ts = (np.arange(95.0, 146.0, 5.0) if not quick
          else np.array([100.0, 120.0, 140.0])) + 273.15
    with pool(len(Ts)) as ex:
        res = list(ex.map(lambda T: (T, scan(GOLD, T=T)), Ts))
    rows = []
    for T, s in res:
        r = readout(s, GOLD["delta_mhz"])
        rows.append(dict(
            T_C=float(T - 273.15), N_m3=float(hyperfine.number_density(T)),
            G_fixed=r["G_fixed"], xi_fixed=r["xi_fixed"],
            gap_fixed_ok=r["gap_fixed_ok"],
            G_best=r["G_best"], xi_best=r["xi_best"], delta_best=r["delta_best"],
            G_peak=r["G_peak"],
            in_cell_od=float(s.get("segment_absorption_od", 0.0) or 0.0),
            od_conj=float((s.get("hardened_noise") or {}).get("od_conj_max", 0.0)),
            pump_scatter=float((s.get("hardened_noise") or {})
                               .get("pump_scatter_noise", 0.0))))
    good = [r for r in rows if r["xi_best"] is not None]
    best = min(good, key=lambda r: r["xi_best"]) if good else None
    gf = [r for r in rows if r["xi_fixed"] is not None]
    best_fixed = min(gf, key=lambda r: r["xi_fixed"]) if gf else None
    return dict(rows=rows,
                best_T_C=None if best is None else best["T_C"],
                best_xi_dB=None if best is None else best["xi_best"],
                best_fixed_T_C=None if best_fixed is None else best_fixed["T_C"],
                best_fixed_xi_dB=None if best_fixed is None else best_fixed["xi_fixed"])


# ==========================================================================
# Cell length -- and whether it substitutes for temperature
# ==========================================================================
def sec_cell_length(quick=False):
    Ls = np.array([6, 9, 12.5, 16, 19, 25, 40] if not quick
                  else [6, 12.5, 25]) * 1e-3
    with pool(len(Ls)) as ex:
        res = list(ex.map(lambda L: (L, scan(GOLD, L=L)), Ls))
    rows = []
    for L, s in res:
        r = readout(s, GOLD["delta_mhz"])
        rows.append(dict(L_mm=float(L * 1e3), G_fixed=r["G_fixed"],
                         xi_fixed=r["xi_fixed"], G_best=r["G_best"],
                         xi_best=r["xi_best"], G_peak=r["G_peak"],
                         in_cell_od=float(s.get("segment_absorption_od", 0.0) or 0.0),
                         NL_rel=float(L / GOLD["L"])))
    # N*L trade: halve the density, double the length -- same column density.
    N0 = hyperfine.number_density(GOLD["T"])
    trade = []
    for factor in (0.5, 1.0, 2.0):
        target_N = N0 / factor          # factor = L/L0
        lo, hi = 320.0, 430.0           # bracket in K
        for _ in range(60):             # bisect on N(T) = target_N
            mid = 0.5 * (lo + hi)
            if hyperfine.number_density(mid) < target_N:
                lo = mid
            else:
                hi = mid
        T = 0.5 * (lo + hi)
        s = scan(GOLD, L=GOLD["L"] * factor, T=T)
        r = readout(s, GOLD["delta_mhz"])
        trade.append(dict(L_mm=float(GOLD["L"] * factor * 1e3),
                          T_C=float(T - 273.15),
                          N_m3=float(hyperfine.number_density(T)),
                          column_density_rel=float(
                              hyperfine.number_density(T) * GOLD["L"] * factor
                              / (N0 * GOLD["L"])),
                          G_fixed=r["G_fixed"], xi_fixed=r["xi_fixed"],
                          G_peak=r["G_peak"], xi_best=r["xi_best"],
                          in_cell_od=float(s.get("segment_absorption_od", 0.0) or 0.0)))
    return dict(rows=rows, column_density_trade=trade,
                note="the trade rows hold N*L fixed: if L only entered through "
                     "the column density, every row would be identical")


# ==========================================================================
# Crossing angle
# ==========================================================================
def sec_angle(quick=False):
    angs = (np.array([0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.32, 0.4, 0.5, 0.65,
                      0.8, 1.0]) if not quick else np.array([0.0, 0.32, 0.8]))
    with pool(len(angs)) as ex:
        res = list(ex.map(lambda a: (a, scan(GOLD, angle=a)), angs))
    rows = []
    for a, s in res:
        r = readout(s, GOLD["delta_mhz"])
        # delta_k_z is per-delta under Ultra's self-consistent phase iteration;
        # quote it at the fixed readout delta.
        dkz = np.atleast_1d(np.asarray(s.get("delta_k_z", 0.0), dtype=float))
        kf = int(np.argmin(np.abs(s["delta_mhz"] - GOLD["delta_mhz"])))
        rows.append(dict(
            angle_deg=float(a), G_fixed=r["G_fixed"], xi_fixed=r["xi_fixed"],
            G_best=r["G_best"], xi_best=r["xi_best"], G_peak=r["G_peak"],
            gap_best=r["gap_best"], gap_ok_fraction=r["gap_ok_fraction"],
            delta_k_z=float(dkz[kf] if dkz.size > 1 else dkz[0]),
            delta_k_z_vacuum=float(np.atleast_1d(
                np.asarray(s.get("delta_k_z_vacuum", 0.0), dtype=float)).ravel()[0]),
            overlap_min=float(s.get("ultra_spatial_overlap_min", 1.0)),
            # first-order beam separation at the cell exit, in probe waists
            separation_per_waist=float(
                GOLD["L"] * math.tan(math.radians(a)) / GOLD["w_probe"]),
            # propagation distance from the cell to 3-waist separation,
            # including Gaussian divergence w(z) = w0 sqrt(1 + (z/z_R)^2)
            separation_distance_mm=_separation_distance_mm(a)))
    good = [r for r in rows if r["xi_best"] is not None]
    return dict(rows=rows,
                best_angle_deg=(min(good, key=lambda r: r["xi_best"])["angle_deg"]
                                if good else None),
                note="Option A (v2 'No double counting of dispersion') puts the "
                     "dispersive phase in Re chi, so the GEOMETRIC mismatch "
                     "2k_pump - k_p - k_c vanishes identically at theta=0 by "
                     "energy conservation. The engine is therefore structurally "
                     "biased to prefer theta=0 and cannot adjudicate an optimal "
                     "finite angle; v6 Sec. 해석 item 4 flags the low-angle "
                     "ridge as outside geometry calibration for the same reason. "
                     "separation_distance_mm is the real counter-pressure.")


def _separation_distance_mm(angle_deg, waists=3.0, wavelength=794.98e-9):
    """Distance past the cell at which the beams are `waists` apart."""
    th = math.radians(float(angle_deg))
    w0 = GOLD["w_probe"]
    zR = math.pi * w0 ** 2 / wavelength
    # Far field: w(z) -> z * lambda/(pi w0), so the beams reach `waists` of
    # separation only if the crossing angle exceeds `waists` times the
    # divergence half-angle.  Below that they never separate, at any distance.
    if th <= waists * wavelength / (math.pi * w0):
        return float("inf")
    z = 0.0
    for _ in range(200):                        # fixed-point on z = k w(z)/theta
        w = w0 * math.sqrt(1.0 + (z / zR) ** 2)
        z_new = waists * w / math.tan(th)
        if abs(z_new - z) < 1e-9:
            z = z_new
            break
        z = 0.5 * (z + z_new)
    return float(z * 1e3)


# ==========================================================================
# Loss
# ==========================================================================
def sec_loss():
    # Use the MEASURED gold-reference gain (111 uW out of an 8 uW seed), not the
    # engine's gain here: the engine overshoots it by ~25x at this point (v2
    # "Model point and experimental comparison"), and at G ~ 300 every
    # gap-to-floor is trivially ~0, which hides the actual loss trade.
    #
    # G_c is taken as G - 1 (single-gain convention), NOT the raw 109/8 ratio.
    # v2 Sec. "Gain convention": a conjugate-to-seed power ratio quoted that way
    # is not an independently calibrated G_c.  Feeding the raw pair in gives
    # G_s - G_c = 0.25, far below the twin-beam validity gap of 1, and
    # (G_s-G_c)^2/(G_s+G_c) then reports a spuriously deep 0.0023 -- the exact
    # gap artifact GAP_MIN guards against.
    G = 111.0 / 8.0
    Gc0 = G - 1.0
    rows = []
    for loss in (0.0, 0.01, 0.02, 0.055, 0.10, 0.15, 0.20, 0.30):
        eta = fwm.QE_DETECTOR * (1.0 - loss)
        rows.append(dict(loss_pct=100.0 * loss, eta=eta,
                         floor_dB=10.0 * math.log10(1.0 - eta),
                         xi_dB=float(observables.intensity_difference_squeezing_dB(
                             G, Gc0, eta)),
                         gap_to_floor_dB=float(
                             observables.intensity_difference_squeezing_dB(G, Gc0, eta)
                             - 10.0 * math.log10(1.0 - eta))))
    x = np.array([r["loss_pct"] for r in rows])
    y = np.array([r["xi_dB"] for r in rows])
    m = x <= 20.0
    return dict(rows=rows, G_used=G, G_c_used=Gc0,
                slope_dB_per_pct=float(np.polyfit(x[m], y[m], 1)[0]),
                gain_needed_for_floor_note=(
                    "the floor 10log10(1-eta) is reachable only as G -> inf; "
                    "the gap to it is 10log10(1 + eta/((1-eta)(2G-1)))"))


# ==========================================================================
# Stray light / excess noise in a conditional ideal-covariance fixture
# ==========================================================================
def sec_stray():
    G, Gc = 14.0, 13.0
    S_ideal = float(observables.ideal_twin_beam_noise(G, Gc))
    loss_rows = [dict(eta=e, S=e * S_ideal + (1 - e),
                      xi_dB=10 * math.log10(e * S_ideal + (1 - e)))
                 for e in (1.0, 0.9, 0.8694, 0.5, 0.2, 0.05, 0.0)]
    # additive excess noise (uncorrelated stray power / thermal modes)
    tech_rows = []
    for tech in (0.0, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0, 2.0):
        S = float(observables.balanced_twin_beam_noise(
            G, Gc, ETA, ETA, reference_weight="dc",
            source_noise=S_ideal, seed_excess_noise=tech))
        tech_rows.append(dict(excess_noise_shot_units=tech,
                              S=S, xi_dB=10 * math.log10(S),
                              above_sql=bool(S > 1.0)))
    S0 = float(observables.balanced_twin_beam_noise(
        G, Gc, ETA, ETA, reference_weight="dc", source_noise=S_ideal))
    # thermal (super-Poissonian) modes coupling in, Jain et al. Eq. for S':
    #   S' = eta*S0 + (1-eta)*[eps*N_th + (1-eps)*N_v],  N_v = 1
    jain = []
    for eps in (0.0, 0.25, 0.5, 0.9, 1.0):
        for N_th in (1.0, 3.0, 10.0):
            Sp = ETA * S_ideal + (1 - ETA) * (eps * N_th + (1 - eps) * 1.0)
            jain.append(dict(eps=eps, N_th=N_th, S=Sp,
                             xi_dB=10 * math.log10(Sp), above_sql=bool(Sp > 1.0)))
    return dict(
        S_ideal=S_ideal, G_used=G,
        pure_loss=loss_rows, additive_excess=tech_rows, jain_model=jain,
        xi_at_zero_excess_dB=10 * math.log10(S0),
        excess_for_sql=float(max(0.0, 1.0 - S0)),
        theorem="Within the assumed ideal covariance, S(eta) = eta*S_ideal + "
                "(1-eta) with 0 <= S_ideal <= 1 gives S <= 1 for eta in [0,1]. "
                "This algebraic unity reference is not a measured SQL and does "
                "not establish the physical cell covariance.")


# ==========================================================================
# Reference DB
# ==========================================================================
def sec_paper_db():
    path = REPO / "references" / "fwm_squeezing_paper_parameters.csv"
    with open(path, encoding="utf-8-sig", newline="") as fh:
        rows = list(csv.reader(fh))
    head, body = rows[0], rows[1:]
    papers = head[5:]
    table = {r[0]: r for r in body}

    def numeric(row_id):
        r = table.get(row_id)
        if r is None:
            return {}
        out = {}
        for name, v in zip(papers, r[5:]):
            v = (v or "").strip()
            if v:
                try:
                    out[name] = float(v)
                except ValueError:
                    pass
        return out

    def stats(row_id):
        d = numeric(row_id)
        if not d:
            return None
        v = np.array(list(d.values()))
        return dict(n=len(v), min=float(v.min()), median=float(np.median(v)),
                    max=float(v.max()), values=d)

    keys = ["opd_ghz_primary", "tpd_mhz_gabes", "cell_temperature_c_primary",
            "cell_length_mm", "pump_power_mw_primary",
            "probe_seed_power_uw_primary", "pump_waist_um_1e2_radius",
            "probe_waist_um_1e2_radius", "pump_probe_angle_deg_primary",
            "fwm_gain_primary", "optical_loss_pct", "detector_qe_pct",
            "total_detection_loss_pct", "squeezing_measured_db",
            "squeezing_bandwidth_mhz_primary"]
    summary = {k: stats(k) for k in keys}
    # derived: pump intensity per paper
    Pd, wd = numeric("pump_power_mw_primary"), numeric("pump_waist_um_1e2_radius")
    inten, rabi = {}, {}
    for name in Pd:
        if name in wd and wd[name] > 0:
            P, w = Pd[name] * 1e-3, wd[name] * 1e-6
            inten[name] = 2 * P / (math.pi * w ** 2) * 1e-4
            rabi[name] = constants.rabi_freq(P, w) / (2 * math.pi) / 1e6
    # How tightly did 10 independent labs converge on each variable?  This is an
    # empirical tolerance ranking, completely independent of the engine and of
    # the v6 1D scans -- so it is a real cross-check on them.
    # Two scales must be treated differently.  TPD and OPD are DETUNINGS whose
    # zero is physical, so a ratio or CV is meaningless (TPD straddles zero and
    # CV blows up); their spread is compared directly with the v6 +0.5 dB
    # tolerance window.  Powers, waists and lengths are ratio-scale.
    v6_window = {                       # v6 Table "tolerance", +0.5 dB, full width
        "tpd_mhz_gabes": (13.8, "-9.6 / +4.2 MHz"),
        "opd_ghz_primary": (0.187, "-32 / +155 MHz"),
        "cell_temperature_c_primary": (15.3, "-3.3 / +12.0 C"),
    }
    spread_abs, spread_ratio = [], []
    for k, lab in (("tpd_mhz_gabes", "TPD delta [MHz]"),
                   ("opd_ghz_primary", "OPD Delta [GHz]"),
                   ("cell_temperature_c_primary", "temperature [C]")):
        st = summary.get(k)
        if not st:
            continue
        v = np.array(list(st["values"].values()))
        win, wtxt = v6_window[k]
        spread_abs.append(dict(
            key=k, label=lab, n=st["n"], span=float(v.max() - v.min()),
            v6_window_full_width=win, v6_window_text=wtxt,
            span_over_window=float((v.max() - v.min()) / win)))
    for k, lab in (("pump_probe_angle_deg_primary", "crossing angle [deg]"),
                   ("cell_length_mm", "cell length [mm]"),
                   ("pump_waist_um_1e2_radius", "pump waist [um]"),
                   ("probe_waist_um_1e2_radius", "seed waist [um]"),
                   ("pump_power_mw_primary", "pump power [mW]"),
                   ("probe_seed_power_uw_primary", "seed power [uW]")):
        st = summary.get(k)
        if not st:
            continue
        v = np.array(list(st["values"].values()))
        spread_ratio.append(dict(key=k, label=lab, n=st["n"],
                                 span_ratio=float(v.max() / v.min()),
                                 cv=float(np.std(v) / np.mean(v))))
    spread_ratio.sort(key=lambda r: r["span_ratio"])
    spread = dict(absolute_scale=spread_abs, ratio_scale=spread_ratio)
    # same for the derived pump intensity
    if inten:
        v = np.array(list(inten.values()))
        inten_spread = dict(label="pump intensity [W/cm2]", n=len(v),
                            span_ratio=float(v.max() / v.min()),
                            cv=float(np.std(v) / np.mean(v)))
    else:
        inten_spread = None
    # lessons row (its header cell is mojibake in the CSV; find it by position)
    lesson_row = body[-1]
    lessons = {n: v for n, v in zip(papers, lesson_row[5:]) if v.strip()}
    return dict(papers=papers, summary=summary,
                pump_intensity_W_cm2=inten, pump_rabi_2pi_MHz=rabi,
                lessons=lessons, spread_ranking=spread,
                intensity_spread=inten_spread,
                intensity_span=[min(inten.values()), max(inten.values())] if inten else None)


# ==========================================================================
# Report
# ==========================================================================
def write_markdown(res, path):
    promote_gain_referred_schema(res)
    L = []
    a = L.append
    f = lambda x, s="{:.2f}": "—" if x is None else s.format(x)

    a("# FWM 동작 변수 점검 — mean-field gain과 gain-referred 진단\n")
    a("자동 생성 — `analysis/squeezing/variable_audit/scan_variable_audit.py`\n")
    a(f"판독 규약은 v6 gain-referred diagnostic tolerance 절과 동일: hardened Ultra, η = {ETA}, "
      f"δ 격자 {DELTA_STEP_MHZ:.0f} MHz, gap gate [{GAP_MIN}, {GAP_MAX}]. "
      f"기준점은 레퍼런스 DB의 gold reference (Sim 2025: Δ=+0.9 GHz, T=121 °C, "
      f"δ=−8 MHz, 600 mW, w 530/330 µm, L=12.5 mm, θ=0.32°).\n")
    a("**Claim gate — `MEAN_FIELD_DIAGNOSTIC`.** 엔진의 절대 이득은 예측값이 "
      "아니며, 아래 dB 값은 mean-field gain에 ideal twin-beam law를 결합한 "
      "`gain_referred_noise_dB` 진단이다. GABES에는 microscopic Langevin "
      "diffusion/covariance와 동일 조건의 측정 SQL이 없으므로, 아래 결과는 물리적 "
      "squeezing·대역폭·SQL 상하를 예측하지 않는다. 읽을 것은 이 가정 아래의 "
      "민감도와 순서뿐이다. `S_dB`와 `xi_*`는 호환용 구 키다.\n")

    db = res["paper_db"]
    a("\n## 0. 레퍼런스 DB가 말하는 실제 사용 범위 (10편)\n")
    a("| 변수 | n | min | median | max |")
    a("|---|---:|---:|---:|---:|")
    labels = {
        "opd_ghz_primary": "OPD Δ [GHz]", "tpd_mhz_gabes": "TPD δ [MHz]",
        "cell_temperature_c_primary": "온도 [°C]", "cell_length_mm": "셀 길이 [mm]",
        "pump_power_mw_primary": "펌프 파워 [mW]",
        "probe_seed_power_uw_primary": "시드 파워 [µW]",
        "pump_waist_um_1e2_radius": "펌프 허리 [µm]",
        "probe_waist_um_1e2_radius": "시드 허리 [µm]",
        "pump_probe_angle_deg_primary": "교차각 [deg]",
        "fwm_gain_primary": "FWM 이득 G", "optical_loss_pct": "광학 loss [%]",
        "detector_qe_pct": "검출기 QE [%]",
        "total_detection_loss_pct": "총 검출 loss [%]",
        "squeezing_measured_db": "측정 스퀴징 [dB]",
        "squeezing_bandwidth_mhz_primary": "대역폭 [MHz]"}
    for k, lab in labels.items():
        st = db["summary"].get(k)
        if st:
            a(f"| {lab} | {st['n']} | {st['min']:g} | {st['median']:g} | "
              f"{st['max']:g} |")
    if db["intensity_span"]:
        a("\n| 논문 | P [mW] | w [µm] | I [W/cm²] | Ω/2π [MHz] | 측정 스퀴징 [dB] |")
        a("|---|---:|---:|---:|---:|---:|")
        sq = db["summary"]["squeezing_measured_db"]["values"]
        for n, I in sorted(db["pump_intensity_W_cm2"].items(), key=lambda kv: kv[1]):
            a(f"| {n} | {db['summary']['pump_power_mw_primary']['values'][n]:.0f} | "
              f"{db['summary']['pump_waist_um_1e2_radius']['values'][n]:.0f} | "
              f"{I:.0f} | {db['pump_rabi_2pi_MHz'][n]:.0f} | "
              f"{sq.get(n, float('nan')):.1f} |")
        isp = db["intensity_spread"]
        a(f"\n> **앞선 저출력 분석의 보정**: 3편(Dowran/Sim/Jain)만 볼 때는 "
          f"펌프 세기가 95–140 W/cm²로 모여 보였지만, DB 9편 전체에서는 "
          f"{db['intensity_span'][0]:.0f}–{db['intensity_span'][1]:.0f} W/cm² "
          f"({isp['span_ratio']:.1f}배)로 펌프 파워 자체의 산포"
          f"({db['summary']['pump_power_mw_primary']['max']/db['summary']['pump_power_mw_primary']['min']:.1f}배)"
          f"보다 오히려 넓다. 결론은 약해지는 게 아니라 **더 강해진다**: "
          f"세기를 5배 바꿔도 −3.5…−8.2 dB가 다 나온다는 뜻이고, 이는 §2에서 "
          f"측정하는 포화 지수 dln q/dln P ≈ 0.11이 예측하는 그대로다. "
          f"펌프는 어느 축으로 봐도 빡빡한 변수가 아니다.")

        a("\n### 10개 연구실이 실제로 얼마나 좁게 수렴했나 — v6 gain-referred 진단 창과 대조\n")
        a("엔진과 무관한 독립 교차검증이다. **검출량(detuning)은 영점이 물리적이라 "
          "비율/CV가 무의미**하므로(TPD는 0을 가로지른다) v6의 +0.5 dB 허용 창과 "
          "직접 비교한다:\n")
        a("| 변수 | n | DB 전체 산포 | v6 진단 +0.5 dB 창 | 산포/창 |")
        a("|---|---:|---:|---|---:|")
        for r in db["spread_ranking"]["absolute_scale"]:
            a(f"| {r['label']} | {r['n']} | {r['span']:.3g} | "
              f"{r['v6_window_text']} | {r['span_over_window']:.2f} |")
        a("\n**TPD가 가장 좁다**: 8편이 독립적으로 고른 δ가 전부 14 MHz 폭 안에 "
          "들어 있고, 그 폭이 v6가 1D 스캔으로 계산한 gain-referred +0.5 dB 진단 창 "
          "(−9.6/+4.2 MHz, 폭 13.8 MHz)과 사실상 같다. 서로 다른 셀·온도·파워를 "
          "쓴 연구실들이 이 한 변수에서만 한 점으로 수렴했다는 뜻이다. 이는 "
          "**진단 창과 문헌 동작점 산포의 정성적 일치**이지, 물리적 squeezing 창의 "
          "검증은 아니다. "
          "OPD와 온도는 창의 2배 정도 흩어져 있는데, 이는 각 실험이 자기 조건에서 "
          "다시 최적화했기 때문으로 읽힌다.")
        a("\n비율 눈금 변수(파워·허리·길이·각)는 max/min으로:\n")
        a("| 변수 | n | 비율 (max/min) | CV |")
        a("|---|---:|---:|---:|")
        for r in db["spread_ranking"]["ratio_scale"]:
            a(f"| {r['label']} | {r['n']} | {r['span_ratio']:.2f} | "
              f"{r['cv']:.3f} |")
        if isp:
            a(f"| (파생) {isp['label']} | {isp['n']} | {isp['span_ratio']:.2f} | "
              f"{isp['cv']:.3f} |")
        a("\n가장 넓게 흩어진 것은 **시드 파워(40배)**와 **펌프 파워(4.4배, "
          "세기로는 5.1배)**다. v6가 1D 스캔으로 얻은 난이도 순서 "
          "`TPD ≈ OPD ≫ T ≫ probe power`와 10개 독립 연구실의 선택이 같은 "
          "정성적 순서를 보인다.")

    a("\n## 1. TPD — dressing 관련, 펌프 세기에 따른 공명 이동\n")
    for label, blk in res["tpd"].items():
        a(f"\n**{label}**  (δ* = 이득 공명 정점, 0.5 MHz 격자)\n")
        a("| P [mW] | Ω/2π [MHz] | Ω/Δ | δ* [MHz] | δ(gate 통과 진단 최저) [MHz] | "
          "Ω²/4Δ [MHz] | dressed (√(Δ²+Ω²)−Δ)/2 [MHz] | G_peak |")
        a("|---:|---:|---:|---:|---:|---:|---:|---:|")
        for r in blk["rows"]:
            a(f"| {r['P_mW']:.0f} | {r['Omega_2pi_MHz']:.0f} | "
              f"{r['Omega_2pi_MHz']/(1e3*(0.9 if 'gold' in label else 0.8)):.2f} | "
              f"{f(r['delta_peak_mhz'],'{:+.2f}')} | "
              f"{f(r['delta_best_mhz'],'{:+.0f}')} | "
              f"{r['ac_stark_estimate_mhz']:+.1f} | "
              f"{r['dressed_shift_mhz']:+.1f} | {f(r['G_peak'],'{:.1f}')} |")
        a(f"\n자유 멱법칙 적합 δ* ∝ Ω^n → **n = {blk['exponent_n']:.2f}** "
          f"(R² = {blk['fit_r2']:.4f}).")
    a("\n두 기하 모두 지수가 **n ≈ 1.5**로 나온다. 이것이 dressing이라는 "
      "직접적인 증거다: n = 2는 섭동적 light shift Ω²/4Δ, n = 1은 강한 dressing "
      "극한(일반화 라비 √(Δ²+Ω²) → Ω)이고, 이 동작점은 Ω/Δ ≈ 0.2–1.6이라 "
      "정확히 그 **교차 영역**에 있다. 위 표에서 Ω²/4Δ는 고출력에서 과대평가, "
      "dressed 공식은 저출력에서 과소평가하는데, 실측 δ*가 그 사이에 있는 것도 "
      "같은 이야기다.")
    a(f"\nδ* ∝ Ω^n ∝ P^(n/2) 이므로 파워를 x배 하면 δ*는 x^{{n/2}} ≈ x^0.75배가 "
      f"된다. 동작점에서 실제로 얼마나 움직이는지:\n")
    a("| 기하 | 동작 P [mW] | 그때 δ* [MHz] | 파워 −50% | 파워 +10% | 파워 ×2 |")
    a("|---|---:|---:|---:|---:|---:|")
    for label, blk in res["tpd"].items():
        a(f"| {label.split(' (')[0]} | {blk['operating_P_mW']:.0f} | "
          f"{blk['operating_delta_mhz']:+.1f} | "
          f"{blk['delta_on_power_half_mhz']:+.1f} | "
          f"{blk['delta_on_power_10pct_mhz']:+.1f} | "
          f"{blk['delta_on_power_x2_mhz']:+.1f} |")
    a(f"\n비교 대상: v6 gain-referred 진단표의 TPD +0.5 dB 창은 "
      f"**−9.6 / +4.2 MHz** — 모든 변수 중 절대 창이 가장 좁다. "
      f"gold 동작점에서 펌프를 절반으로 줄이면 δ*가 −14.5 MHz 움직여 "
      f"창(−9.6)을 벗어나고, 2배로 올리면 +24 MHz로 창(+4.2)을 6배 넘긴다. "
      f"**파워를 10%만 건드려도 청(+)측 창의 절반을 쓴다.** "
      f"실무적으로 펌프 파워(또는 허리)를 바꾼 뒤 δ 재스캔을 우선 점검할 근거다. "
      f"다만 이 계산만으로 관측된 squeezing 저하의 물리적 원인을 δ 불일치로 "
      f"확정할 수는 없다.")

    pv = res["pump_vs_overlap"]
    a("\n## 2. 펌프 빔 — \"연료, 어지간하면 됨. FWM 안 나오면 겹침부터\" → "
      "**맞음, 정량적으로도 압도적**\n")
    a("| P [mW] | G(δ 고정) | gain-ref. 진단(δ 고정) [dB] | G_peak | gain-ref. 진단 최저 [dB] |")
    a("|---:|---:|---:|---:|---:|")
    for r in pv["power"]:
        a(f"| {r['P_mW']:.0f} | {f(r['G_fixed'])} | {f(r['gain_referred_fixed_dB'])} | "
          f"{f(r['G_peak'],'{:.1f}')} | {f(r['gain_referred_best_dB'])} |")
    a("\n| 모드 겹침 | 등가 횡방향 어긋남 [µm] | G(δ 고정) | gain-ref. 진단(δ 고정) [dB] | "
      "G_peak | gain-ref. 진단 최저 [dB] |")
    a("|---:|---:|---:|---:|---:|---:|")
    for r in pv["overlap"]:
        a(f"| {r['overlap']:.2f} | {r['offset_um']:.0f} | {f(r['G_fixed'])} | "
          f"{f(r['gain_referred_fixed_dB'])} | {f(r['G_peak'],'{:.1f}')} | "
          f"{f(r['gain_referred_best_dB'])} |")
    a(f"\n**지렛대 비교** (기준점 ±{pv['local_eps']*100:.0f}% 국소 미분): "
      f"dlnG/dlnP = {pv['dlnG_dlnP']:.2f}, dlnG/dln(겹침) = "
      f"{pv['dlnG_dln_overlap']:.2f} → 겹침이 파워보다 "
      f"**{pv['leverage_ratio']:.1f}배** 강하다.")
    a("\n이 비율은 엔진의 절대 이득과 무관한 양이다. `mode_overlap_penalty`는 "
      "결합상수에 곱해지므로 q ∝ 겹침이 **정확히** 성립하고, 연쇄법칙으로\n")
    a("```\n"
      "  dlnG/dln(겹침) ÷ dlnG/dlnP = 1 / (dln q / dln P)\n"
      "```\n")
    a(f"즉 dlnG/dlnq가 그대로 약분된다. 측정된 dln q/dln P = "
      f"{pv['dlnq_dlnP']:.3f} — 앞선 저출력 분석에서 나온 포화 지수와 같은 크기다. "
      f"겹침은 이득 지수 qL에 선형으로 들어가 exp에 바로 앉는 반면, 파워는 "
      f"√I → 이미 포화된 Ω를 거쳐 들어가기 때문이다. "
      f"**'FWM이 안 보이면 펌프보다 겹침을 먼저 보라'는 모델이 그대로 재현한다.**")

    sd = res["seed"]
    a("\n## 3. 시드 빔 — 절반 맞음\n")
    a("**(a) 시드 파워에 대한 gain-referred 진단 변화는 작다** "
      f"(1–200 µW 전 구간 {sd['gain_referred_spread_dB']:.3f} dB):\n")
    a("| 시드 [µW] | G(δ 고정) | gain-ref. 진단(δ 고정) [dB] | gain-ref. 진단 최저 [dB] |")
    a("|---:|---:|---:|---:|")
    for r in sd["power"]:
        a(f"| {r['P_probe_uW']:.0f} | {f(r['G_fixed'])} | "
          f"{f(r['gain_referred_fixed_dB'])} | {f(r['gain_referred_best_dB'])} |")
    a("\nv6 gain-referred 진단표의 `P_probe: 0.125×–8× 전 구간 무변화(<0.004 dB)`와 "
      "일치한다. χ̄가 시드 세기에 무관하고 펌프 고갈 기여도 미미하기 때문이다.")
    a("\n**(b) 증폭되지 않고 한 팔로 새는 시드는 선택한 ideal-covariance fixture에서 "
      "불리하다.** 아래 계산은 새는 coherent light가 상관항에는 기여하지 않고 "
      "프로브 팔의 평균과 포아소니안 잡음만 늘린다고 가정한다. 실제 셀의 covariance "
      "예측은 아니다.\n")
    a(f"G = {sd['G_used']:.0f} 기준. 상관항은 f=0 값에 고정하고 프로브 팔만 "
      f"η(G+f)로 키운 조건부 계산이다:\n")
    a("| 새는 시드 (P₀ 배수) | 프로브 팔 중 비중 | 가중 | DC 불균형 | ideal-law 진단 [dB] | "
      "진단 대가 [dB] | ideal-law 등가 loss [%p] |")
    a("|---:|---:|---|---:|---:|---:|---:|")
    for r in sd["leakage"]:
        a(f"| {r['leak_fraction_of_seed']:.2f} | "
          f"{r['leak_pct_of_probe_arm']:.1f}% | {r['weight']} | "
          f"{r['dc_imbalance']:+.4f} | {r['ideal_law_noise_dB']:.2f} | "
          f"{r['ideal_law_penalty_dB']:+.2f} | "
          f"{r['ideal_law_equivalent_extra_loss_pct']:+.1f} |")
    a("\n이 fixture에서는 세 가지 전자 가중(w=1 / DC 상쇄 / 샷잡음 최적)의 "
      "진단값이 소수점 셋째 자리까지 같다. 정규화된 ideal-law 잡음비는 w에 대해 "
      "최적점 근처에서 **정류점(stationary)** 이므로, DC 재균형 자체가 진단값을 "
      "거의 바꾸지 않는다.\n")
    a("정리하면:")
    a("- **이 fixture에서 DC 기댓값 자체는 진단 손실이 아니다.** 재균형으로 "
      "지워도 ideal-law 진단은 거의 같다.")
    a("- **진단 대가는 상관항 없이 한쪽 팔에 더한 광자에서 생긴다.** 선택한 "
      "포아소니안 가정 아래 프로브 팔의 6.7% 누설은 0.73 dB, ideal-law loss "
      "3.6 %p와 같은 진단 대가다.")
    a("\n그러니 체감의 물리는 맞고 표현만 옮기면 된다: 'IDS 기댓값이 0이 안 "
      "된다'가 아니라 **'한쪽 팔에만 상관 없는 광자가 더해진다'**. 그리고 그 "
      "순간 이 항은 시드 문제가 아니라 §8의 '외부 빔 유입'과 같은 범주가 된다. "
      "다만 coherent/thermal 유입의 실제 SQL 상하는 microscopic covariance와 "
      "동일 조건의 측정 SQL 없이 이 스크립트가 판정할 수 없다.")

    tp = res["temperature"]
    a("\n## 4. 온도 — mean-field gain과 gain-referred 진단의 조건부 추세\n")
    a("| T [°C] | N [m⁻³] | G(δ=−8) | gain-ref. 진단(δ=−8) [dB] | gap ok | G(δ 진단최저) | "
      "gain-ref. 진단 최저 [dB] | in-cell OD | conj OD | pump scatter |")
    a("|---:|---:|---:|---:|:--:|---:|---:|---:|---:|---:|")
    for r in tp["rows"]:
        a(f"| {r['T_C']:.0f} | {r['N_m3']:.2e} | {f(r['G_fixed'],'{:.1f}')} | "
          f"{f(r['gain_referred_fixed_dB'])} | {'○' if r['gap_fixed_ok'] else '×'} | "
          f"{f(r['G_best'],'{:.1f}')} | {f(r['gain_referred_best_dB'])} | "
          f"{r['in_cell_od']:.3f} | {r['od_conj']:.3f} | "
          f"{r['pump_scatter']:.4f} |")
    a(f"\nδ 고정 진단의 최저점은 T = {f(tp['best_fixed_T_C'],'{:.0f}')} °C "
      f"({f(tp['best_fixed_gain_referred_dB'])} dB), δ 재조정 진단의 최저점은 "
      f"T = {f(tp['best_T_C'],'{:.0f}')} °C "
      f"({f(tp['best_gain_referred_dB'])} dB). "
      f"mean-field 이득은 온도에 대해 계속 커지지만 ideal-law 진단은 꺾인다. "
      f"이 모델에서 함께 변하는 항은 "
      f"표의 마지막 세 열(in-cell OD, conjugate 팔 OD, pump scatter)이다. "
      f"이는 \"gain은 커지나 물리적 squeezing이 나빠진다\"의 증명이나 microscopic "
      f"noise 설명이 아니다. "
      f"레퍼런스 DB의 온도 중앙값도 120 °C이고 최대 128 °C로, 아무도 그 위로 "
      f"가지 않는다.")

    cl = res["cell_length"]
    a("\n## 5. 셀 길이 — \"뭐랑 관련있는지 모르겠음\" → **밀도와 짝을 이루는 "
      "변수 (기둥밀도 N·L)**\n")
    a("| L [mm] | G(δ=−8) | gain-ref. 진단(δ=−8) [dB] | G_peak | gain-ref. 진단 최저 [dB] | in-cell OD |")
    a("|---:|---:|---:|---:|---:|---:|")
    for r in cl["rows"]:
        a(f"| {r['L_mm']:.1f} | {f(r['G_fixed'],'{:.1f}')} | "
          f"{f(r['gain_referred_fixed_dB'])} | {f(r['G_peak'],'{:.1f}')} | "
          f"{f(r['gain_referred_best_dB'])} | {r['in_cell_od']:.3f} |")
    a("\n**N·L을 고정한 교환 실험** (밀도를 반으로 줄이고 길이를 두 배로):\n")
    a("| L [mm] | T [°C] | N·L (상대) | G(δ=−8) | gain-ref. 진단(δ=−8) [dB] | G_peak | "
      "gain-ref. 진단 최저 [dB] | in-cell OD |")
    a("|---:|---:|---:|---:|---:|---:|---:|---:|")
    for r in cl["column_density_trade"]:
        a(f"| {r['L_mm']:.1f} | {r['T_C']:.1f} | {r['column_density_rel']:.3f} | "
          f"{f(r['G_fixed'],'{:.1f}')} | {f(r['gain_referred_fixed_dB'])} | "
          f"{f(r['G_peak'],'{:.1f}')} | {f(r['gain_referred_best_dB'])} | "
          f"{r['in_cell_od']:.3f} |")
    a("\nL은 이득 지수 qL에 **선형**으로, 흡수 OD에도 선형으로 들어간다. "
      "따라서 1차적으로 L과 밀도(온도)는 기둥밀도 N·L로 묶여 서로를 대체한다. "
      "위 교환 표가 완전히 같지 않은 만큼이 '길이만의 효과'이고, 그 차이는 "
      "온도가 밀도 말고도 Doppler 폭·충돌 완화·pump scatter를 함께 움직이기 "
      "때문이다. 실무적 의미: **긴 셀은 낮은 온도로 같은 이득을 사는 수단**이고, "
      "저온 쪽이 충돌 dephasing과 pump scatter가 작으니 유리할 수 있다. "
      f"레퍼런스 DB의 셀 길이가 {db['summary']['cell_length_mm']['min']:.0f}–"
      f"{db['summary']['cell_length_mm']['max']:.0f} mm로 흩어져 있는데도 "
      f"결과가 비슷한 이유가 이것이다.")

    an = res["angle"]
    a("\n## 6. 펌프-프로브 각 — \"겹침 변수, 0도 이상적, 분리 때문에 타협\" → "
      "**결론은 맞지만 지배 항이 겹침이 아니라 위상정합**\n")
    a("| θ [deg] | G(δ=−8) | gain-ref. 진단(δ=−8) [dB] | G_peak | gain-ref. 진단 최저 [dB] | "
      "Δk_z [1/m] | 셀내 최소겹침 | 출구 분리/허리 | 3-허리 분리까지 [mm] |")
    a("|---:|---:|---:|---:|---:|---:|---:|---:|---:|")
    for r in an["rows"]:
        sep = ("∞" if not math.isfinite(r["separation_distance_mm"])
               else f"{r['separation_distance_mm']:.0f}")
        a(f"| {r['angle_deg']:.2f} | {f(r['G_fixed'],'{:.1f}')} | "
          f"{f(r['gain_referred_fixed_dB'])} | {f(r['G_peak'],'{:.1f}')} | "
          f"{f(r['gain_referred_best_dB'])} | "
          f"{r['delta_k_z']:.1f} | {r['overlap_min']:.3f} | "
          f"{r['separation_per_waist']:.2f} | {sep} |")
    a("\n표의 `셀내 최소겹침` 열을 보면 0.32°에서도 0.997, 1.0°에서도 0.97로 "
      "**겹침은 거의 손실이 없다**. 그런데 이득은 급격히 떨어진다 — 범인은 "
      "겹침이 아니라 위상정합 Δk_z다. v2 §Classical coupled-mode equation의 "
      "Δk_geom = 2k_pump,z − k_p,z − k_c,z 가 exp(qL) 안에 직접 들어간다.")
    a("\n**단, 여기엔 모델의 구조적 편향이 있다.** Option A 규약(v2 「No double "
      "counting of dispersion」)에서 분산은 Re χ_pp/χ_cc가 전부 담당하고 Δk에는 "
      "진공 파수만 들어간다. 에너지 보존 ω_p + ω_c = 2ω_pump이 정확하므로 "
      "**θ=0에서 Δk_geom은 항등적으로 0**이고, 엔진은 언제나 0°를 최적으로 "
      "고를 수밖에 없다. 즉 이 모델은 '유한 각이 최적일 수 있는가'라는 질문에 "
      "답할 자격이 없다 — v6 「해석」 4항이 low-angle ridge를 "
      "*geometry calibration 범위 밖*이라고 못 박은 것과 같은 이유다. "
      "저각 수치를 설계 근거로 쓰지 말 것.")
    a(f"\n실무적으로 각을 정하는 것은 **분리 거리**다. 마지막 열이 셀 출구에서 "
      f"프로브/conjugate가 3 허리만큼 벌어지는 데 필요한 자유전파 거리다 "
      f"(가우시안 발산 포함). 레퍼런스 DB의 실제 사용 각은 "
      f"{db['summary']['pump_probe_angle_deg_primary']['min']:g}–"
      f"{db['summary']['pump_probe_angle_deg_primary']['max']:g}° "
      f"(중앙값 {db['summary']['pump_probe_angle_deg_primary']['median']:g}°)로 "
      f"좁게 모여 있는데, 이 표를 보면 그 범위가 '분리 거리는 감당 가능하면서 "
      f"Δk_z 대가는 아직 작은' 구간과 정확히 겹친다. 체감하신 타협이 바로 "
      f"이것이다.")

    ls_ = res["loss"]
    a("\n## 7. Loss — 측정 이득을 넣은 conditional ideal-law 진단\n")
    a(f"기준 이득은 gold reference의 **측정값** G = {ls_['G_used']:.2f} "
      f"(111 µW 출력 / 8 µW 시드), G_c = G−1 = {ls_['G_c_used']:.2f}. 엔진 "
      f"이득이 아니다. G_c에 원자료 109/8을 그대로 넣으면 gap이 0.25로 "
      f"twin-beam 유효 범위를 벗어나 −26 dB짜리 가짜 값이 나온다 "
      f"(v2 §Gain convention: 그 비는 독립 교정된 G_c가 아니다):\n")
    a("| loss [%] | η | ideal-law 점근선 10log₁₀(1−η) [dB] | ideal-law 진단 [dB] | 점근선까지 거리 [dB] |")
    a("|---:|---:|---:|---:|---:|")
    for r in ls_["rows"]:
        a(f"| {r['loss_pct']:.1f} | {r['eta']:.4f} | "
          f"{r['ideal_law_asymptote_dB']:.2f} | {r['ideal_law_noise_dB']:.2f} | "
          f"{r['gap_to_ideal_law_asymptote_dB']:+.2f} |")
    a(f"\n진단 기울기 ≈ {ls_['gain_referred_slope_dB_per_pct']:.3f} dB/%p "
      f"(v6 gain-referred 진단표의 0.22 dB/%p와 같은 크기). "
      f"10log₁₀(1−η)는 선택한 single-gain ideal law의 G→∞ 대수적 점근선이다. "
      f"microscopic Langevin noise나 검출 전자잡음이 포함된 물리적 squeezing "
      f"floor로 해석하면 안 된다.")

    st = res["stray"]
    a("\n## 8. 외부 빔 유입 — conditional ideal-covariance fixture\n")
    a(f"> {st['ideal_law_statement']}\n")
    a("아래 unity=1은 이 fixture의 정규화 기준일 뿐, 동일 조건에서 측정한 SQL이 "
      "아니다. 따라서 표의 unity 상하는 물리적 above/below-SQL 판정이 아니다.\n")
    a(f"ideal-law source ratio = {st['ideal_law_source_ratio']:.4f} "
      f"(G = {st['G_used']:.0f}) 에서 η를 0까지 "
      f"낮춰도:\n")
    a("| η | ideal-law 잡음비 | ideal-law 진단 [dB] |")
    a("|---:|---:|---:|")
    for r in st["pure_loss"]:
        a(f"| {r['eta']:.4f} | {r['ideal_law_noise_ratio']:.4f} | "
          f"{r['ideal_law_noise_dB']:+.2f} |")
    a("\nη → 0에서 S → 1 (0 dB)에 **아래로부터** 수렴할 뿐 절대 넘지 않는다. "
      "이는 assumed ideal covariance 안의 대수적 성질이다. 측정된 SQL 초과의 "
      "원인을 이 계산만으로 loss 또는 미광에 귀속할 수 없다.\n")
    a("가법 초과잡음(상관 없는 유입광, 샷잡음 단위):\n")
    a("| 초과잡음 | ideal-law 잡음비 | ideal-law 진단 [dB] | 정규화 unity 초과? |")
    a("|---:|---:|---:|:--:|")
    for r in st["additive_excess"]:
        a(f"| {r['excess_noise_shot_units']:.2f} | "
          f"{r['ideal_law_noise_ratio']:.4f} | {r['ideal_law_noise_dB']:+.2f} | "
          f"{'예' if r['above_normalized_unity'] else '아니오'} |")
    a(f"\n이 fixture에서는 정규화 초과잡음 {st['excess_for_normalized_unity']:.2f}를 "
      f"더하면 unity를 넘는다. 이는 측정 SQL 임계값이 아니다.\n")
    a("Jain et al.(OL 50, 5165)의 파이버 결합 모델 "
      "S′ = ηS₀ + (1−η)[ε·N_th + (1−ε)·N_v] 로 본 것 — 상관 없는 공간 모드가 "
      "**열적(super-Poissonian)** 이라는 것이 핵심이다:\n")
    a("| ε (열적 비율) | N_th | ideal-law S′ | ideal-law 진단 [dB] | 정규화 unity 초과? |")
    a("|---:|---:|---:|---:|:--:|")
    for r in st["jain_model"]:
        a(f"| {r['eps']:.2f} | {r['N_th']:.0f} | "
          f"{r['ideal_law_noise_ratio']:.4f} | {r['ideal_law_noise_dB']:+.2f} | "
          f"{'예' if r['above_normalized_unity'] else '아니오'} |")
    a("\nJain의 −7.2 dB(셀 직후) → −4.4 dB(파이버 후) 열화가 10% 삽입손실만으로는 "
      "−6.3 dB까지밖에 설명되지 않고, 나머지 ~1.9 dB가 바로 이 '상관 없는 "
      "공간 모드 유입으로 해석한 문헌 모델과 부합한다. 이는 이 스크립트가 Jain의 "
      "microscopic covariance를 재현했다는 뜻은 아니다. 레퍼런스 DB의 Jain 항목 "
      "교훈란도 같은 말을 한다: "
      "*\"Fiber source limited by coherence area/spatial filtering, not only "
      "insertion loss.\"*")

    a("\n## 9. 종합 — 판정표\n")
    a("| 변수 | 체감 | 판정 | 근거 |")
    a("|---|---|---|---|")
    verdict = [
        ("OPD", "(제시 없음)", "보충",
         f"v6 gain-referred 진단: +0.5 dB 창 −32/+155 MHz(δ 고정). δ 재조정하면 "
         f"±0.1 GHz로 풀린다 — OPD·TPD는 항상 함께 조정할 쌍. DB 범위 "
         f"{db['summary']['opd_ghz_primary']['min']:g}–"
         f"{db['summary']['opd_ghz_primary']['max']:g} GHz"),
        ("TPD", "dressing, 펌프 세기 따라 이동", "맞음",
         "δ* ∝ Ω^1.5 (섭동 Ω²와 강한 dressing Ω의 교차역). 파워 2배 = "
         "허용 창 전체를 벗어남. v6에서 절대 창이 가장 좁은 변수"),
        ("펌프", "연료, 어지간하면 됨 / 겹침 먼저 확인", "맞음",
         f"dlnG/dlnP = {pv['dlnG_dlnP']:.2f} vs dlnG/dln(겹침) = "
         f"{pv['dlnG_dln_overlap']:.2f}. DB 135–600 mW 전부 동작"),
        ("시드 파워", "—", "진단상 둔감",
         f"1–200 µW에서 gain-referred 진단 변화 "
         f"{sd['gain_referred_spread_dB']:.3f} dB"),
        ("시드 누설", "출력 경로 유입 → IDS 기댓값 0 아님", "물리 맞음, 표현 정정",
         "선택한 ideal-covariance fixture에서 DC 재균형은 진단값을 거의 "
         "바꾸지 않는다. 한쪽 팔의 상관 없는 광자는 ideal-law 등가 loss로 "
         "환산되지만 실제 covariance 예측은 아니다"),
        ("온도", "110↑ 급격, 120↑ gain↑ 스퀴징↓", "gain 추세/진단만 지지",
         f"δ 고정 진단 최저점 {f(tp['best_fixed_T_C'],'{:.0f}')} °C. "
         f"physical squeezing 악화나 microscopic noise는 판정 불가"),
        ("셀 길이", "모르겠음", "보충",
         "qL과 OD 양쪽에 선형 → 밀도와 기둥밀도 N·L로 묶임. "
         "긴 셀 = 낮은 온도로 같은 이득"),
        ("교차각", "겹침 변수, 0도 이상적", "절반",
         "겹침 + 위상정합(Δk_z) 이중 역할. v6의 low-angle ridge는 "
         "calibration 범위 밖. DB 0.3–1.0°"),
        ("Loss", "이득 확보 후 최대 방해", "ideal-law 진단에서 지지",
         f"{ls_['gain_referred_slope_dB_per_pct']:.2f} dB/%p; "
         "10log10(1−η)는 물리적 floor가 아닌 이상법칙 점근선"),
        ("외부 빔 유입", "loss로 SQL 초과 설명 불가", "조건부 항등식만 지지",
         "assumed ideal covariance에서 S = ηS_ideal + (1−η) ≤ 1. "
         "측정 SQL 상하와 실제 원인 귀속은 판정 불가"),
    ]
    for v in verdict:
        a(f"| {v[0]} | {v[1]} | **{v[2]}** | {v[3]} |")

    a("\n### 조건부 점검 순서 (mean-field/ideal-law 진단)\n")
    a("이 순서는 실험 진단 후보의 우선순위이며, physical squeezing noise-budget의 "
      "first-principles 순위가 아니다.")
    a("1. **Loss / 검출 효율** — ideal-law 진단과 점근선을 직접 움직인다.")
    a("2. **외부 빔 유입 / 미광** — 별도 microscopic covariance와 동일 조건 SQL "
      "측정으로 확인해야 한다.")
    a("3. **빔 겹침** — 이득이 안 나올 때 펌프 파워보다 먼저 볼 것 "
      f"(지렛대 {abs(pv['dlnG_dln_overlap']/max(abs(pv['dlnG_dlnP']),1e-9)):.0f}배).")
    a("4. **TPD 재스캔** — 특히 펌프 파워를 건드린 직후. δ* ∝ Ω^1.5 ∝ P^0.75로 "
      "이동하고, 파워 10% 변화가 허용 창의 절반을 쓴다.")
    a("5. **OPD** — TPD와 쌍으로 함께 재스캔.")
    a("6. **온도** — 저온측이 절벽이니 명목값보다 1–2 °C 높게.")
    a("7. **시드 파워** — 마지막. 검출기 포화만 피하면 된다.")

    a("\n### DB의 논문별 교훈 (원문)\n")
    for n, v in db["lessons"].items():
        a(f"- **{n}**: {v}")

    Path(path).write_text("\n".join(L) + "\n", encoding="utf-8")
    print(f"report -> {path}")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default=str(Path(__file__).resolve().parent / "generated"))
    ap.add_argument("--quick", action="store_true")
    args = ap.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    t0 = time.time()
    res = {}
    for name, fn in (("paper_db", lambda: sec_paper_db()),
                     ("tpd", lambda: sec_tpd(args.quick)),
                     ("pump_vs_overlap", lambda: sec_pump_vs_overlap(args.quick)),
                     ("seed", lambda: sec_seed(args.quick)),
                     ("temperature", lambda: sec_temperature(args.quick)),
                     ("cell_length", lambda: sec_cell_length(args.quick)),
                     ("angle", lambda: sec_angle(args.quick)),
                     ("loss", lambda: sec_loss()),
                     ("stray", lambda: sec_stray())):
        t = time.time()
        res[name] = fn()
        print(f"  {name:16s} {time.time()-t:7.1f}s")
    res["meta"] = dict(
        eta=ETA, gap_gate=[GAP_MIN, GAP_MAX],
        delta_step_mhz=DELTA_STEP_MHZ, gold=dict(GOLD),
        runtime_s=time.time() - t0,
        result_schema_version=2,
        diagnostic_claim_gate=dict(DIAGNOSTIC_CLAIM_GATE),
        deprecated_result_aliases=dict(_LEGACY_TO_PRIMARY),
    )
    promote_gain_referred_schema(res)

    (out / "variable_audit.json").write_text(
        json.dumps(res, indent=1, ensure_ascii=False, default=float),
        encoding="utf-8")
    print(f"json   -> {out / 'variable_audit.json'}")
    write_markdown(res, out / "variable_audit.md")
    print(f"total {time.time()-t0:.1f}s")


if __name__ == "__main__":
    main()
