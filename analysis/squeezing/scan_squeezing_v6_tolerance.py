"""
v6 tolerance scan: how far each experimentally tunable knob can drift from the
v6 finite-loss hardened-Ultra optimum before the squeezing degrades past a
threshold.

Center (v6 global optimum, squeezing_report_v6.tex):
  Delta = -1.50 GHz, T = 110 C, delta(TPD) = -280 MHz,
  P_pump = 0.6 W, P_probe = 8 uW, line_strength = 0.74,
  QE = 0.92, loss = 0.055 (eta = 0.8694)  ->  xi = -8.102 dB, G_s = 22.6.

Scanned knobs (one at a time, everything else pinned at the center):
  1. OPD   (one-photon detuning Delta)      - re-solve per point
  2. TPD   (two-photon detuning delta)      - one fine-delta solve at the center
  3. T     (cell temperature)               - re-solve per point
  4. P_probe (seed power)                   - re-solve per point
  5. loss  (detection loss fraction)        - analytic eta remap, no re-solve

Two readouts for the re-solve scans:
  * fixed-delta: xi at delta = -280 MHz exactly (pure drift: nothing re-tuned).
    This is the physical drift scenario when probe is AOM/EOM-derived from the
    pump (TPD stays locked while OPD/T/power drift).
  * best-delta:  xi at the gap-gated deepest delta (experimenter re-tunes TPD).

Trust gate: same twin-beam gap gate as the frontier scan
  (GAP_MIN <= G_s - G_c <= GAP_MAX), delta-window edges excluded for best-delta.

Outputs:
  analysis/squeezing/squeezing_v6_tolerance/squeezing_v6_tolerance.npz
  docs/squeezing_report/squeezing_v6_tolerance.png

Run from the repo root:
  python analysis/squeezing/scan_squeezing_v6_tolerance.py
"""
from __future__ import annotations

import os
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

for _stream in (sys.stdout, sys.stderr):
    try:
        _stream.reconfigure(encoding="utf-8", errors="replace")
    except Exception:
        pass

from gabes.core import blas_single_thread  # noqa: E402
from gabes.schemes import fwm  # noqa: E402
from gabes import observables  # noqa: E402

OUT = ROOT / "analysis" / "squeezing" / "squeezing_v6_tolerance"
DOC_OUT = ROOT / "docs" / "squeezing_report" / "squeezing_v6_tolerance.png"

# ---- v6 finite-loss optimum (center point) ---------------------------------
D0_GHZ = -1.50
T0_C = 110.0
TPD0_MHZ = -280.0
P_PUMP = 0.6
P_PROBE0 = 8e-6
LINE_STRENGTH = 0.74
LOSS0 = 0.055
QE = fwm.QE_DETECTOR          # 0.92
BRANCH = -1

# Same delta window as the v6 frontier scan, but a 5 MHz delta grid throughout:
# the frontier's 17.5 MHz spacing misses valid-gap delta points next to the dip
# (observed at Delta=-1.45 / T=114-116), which fakes a collapse of the
# best-delta readout near the gap-gate boundary.
WINDOW_GHZ = 0.7
COARSE = 281                   # 5 MHz spacing; TPD0 = -280 sits on-grid
COARSE_FINE = 281              # same grid for the TPD tolerance curve
VELOCITY_STEP = 5.0
VELOCITY_CUTOFF = 3.0

GAP_MIN = 0.5
GAP_MAX = 1.5

# ---- scan axes --------------------------------------------------------------
OPD_AXIS_GHZ = np.round(np.arange(-2.30, -0.70 + 1e-9, 0.05), 4)       # 33 pts
TEMP_AXIS_C = np.round(np.arange(88.0, 132.0 + 1e-9, 2.0), 4)          # 23 pts
PPROBE_FACTORS = np.array([2.0 ** (k / 2.0) for k in range(-6, 7)])    # 13 pts
LOSS_AXIS = np.round(np.arange(0.0, 0.351 + 1e-9, 0.005), 4)           # 71 pts

THRESHOLDS_DB = (0.25, 0.5, 1.0)   # degradation above the center xi


def _solve(D_GHz, T_C, P_probe, coarse=COARSE):
    """One hardened-Ultra solve; returns the full delta spectrum + hardened arrays."""
    center = fwm.branch_center_GHz(D_GHz, BRANCH)
    with blas_single_thread():
        spec = fwm.compute_spectrum(
            D_GHz, T=T_C + 273.15, P_pump=P_PUMP, P_probe=P_probe,
            line_strength=LINE_STRENGTH, loss_frac=LOSS0, qe=QE,
            coarse_points=coarse, fine_points=0,
            scan_min=center - WINDOW_GHZ, scan_max=center + WINDOW_GHZ,
            velocity_step=VELOCITY_STEP, velocity_cutoff=VELOCITY_CUTOFF,
            branch=BRANCH,
            phase_detail=fwm.PHASE_ULTRA, model_fidelity=fwm.FIDELITY_ULTRA)
    tpd_mhz = (spec["probe_axis_GHz"] - center) * 1e3
    hn = spec.get("hardened_noise") or {}
    n = tpd_mhz.size
    return dict(
        tpd_mhz=np.asarray(tpd_mhz, float),
        xi=np.asarray(spec["S_dB"], float),
        Gs=np.asarray(spec["G_s"], float),
        Gc=np.asarray(spec["G_c"], float),
        seg_od=float(spec.get("segment_absorption_od", 0.0) or 0.0),
        od_conj=np.asarray(hn.get("od_conj_arr", np.zeros(n)), float),
        od_probe=np.asarray(hn.get("od_probe_lin_arr", np.zeros(n)), float),
        pump_scatter=float(hn.get("pump_scatter_noise", 0.0)),
    )


def _readouts(sol):
    """(fixed-delta xi/gap, best-delta xi/gap/tpd) from one solved spectrum."""
    tpd, xi, gap = sol["tpd_mhz"], sol["xi"], sol["Gs"] - sol["Gc"]
    i0 = int(np.argmin(np.abs(tpd - TPD0_MHZ)))
    fixed = dict(xi=float(xi[i0]), gap=float(gap[i0]),
                 trusted=bool(GAP_MIN <= gap[i0] <= GAP_MAX),
                 on_grid=bool(abs(tpd[i0] - TPD0_MHZ) < 1e-6), idx=i0)
    n = xi.size
    interior = np.ones(n, bool)
    interior[:2] = interior[-2:] = False
    valid = np.isfinite(xi) & interior & (gap >= GAP_MIN) & (gap <= GAP_MAX)
    if valid.any():
        j = int(np.nanargmin(np.where(valid, xi, np.inf)))
        best = dict(xi=float(xi[j]), gap=float(gap[j]), tpd=float(tpd[j]),
                    trusted=True, idx=j)
    else:
        j = int(np.nanargmin(xi))
        best = dict(xi=float(xi[j]), gap=float(gap[j]), tpd=float(tpd[j]),
                    trusted=False, idx=j)
    return fixed, best


def _xi_of_loss(sol, idx, loss_axis):
    """Analytic eta remap at one delta index (hardened balanced formula).

    Reproduces compute_spectrum's hardened S_dB: probe arm eta_s =
    eta*exp(-(seg_od+od_probe)), conjugate arm eta_c = eta*exp(-od_conj),
    pump-scatter excess added to the normalized noise ratio.
    """
    Gs, Gc = sol["Gs"][idx], sol["Gc"][idx]
    out = np.empty(loss_axis.size)
    for k, ell in enumerate(loss_axis):
        eta = QE * (1.0 - ell)
        eta_s = eta * np.exp(-(max(sol["seg_od"], 0.0)
                               + max(sol["od_probe"][idx], 0.0)))
        eta_c = eta * np.exp(-max(sol["od_conj"][idx], 0.0))
        S = observables.balanced_twin_beam_noise(
            Gs, Gc, eta_s, eta_c, reference_weight="dc")
        S = S + max(sol["pump_scatter"], 0.0)
        out[k] = 10.0 * np.log10(max(S, 1e-30))
    return out


def _crossings(x, xi, x_center, xi_center, thr_db):
    """(x_lo, x_hi) where xi first exceeds xi_center + thr on each side of the
    center (linear interp). np.nan if never crossed within the scan range."""
    target = xi_center + thr_db
    ic = int(np.argmin(np.abs(x - x_center)))
    lo = hi = np.nan
    for i in range(ic, x.size - 1):                      # upward side
        a, b = xi[i], xi[i + 1]
        if np.isfinite(a) and np.isfinite(b) and a <= target < b:
            hi = x[i] + (x[i + 1] - x[i]) * (target - a) / (b - a)
            break
    for i in range(ic, 0, -1):                           # downward side
        a, b = xi[i], xi[i - 1]
        if np.isfinite(a) and np.isfinite(b) and a <= target < b:
            lo = x[i] + (x[i - 1] - x[i]) * (target - a) / (b - a)
            break
    return lo, hi


def _fmt(v, unit=""):
    return "  beyond scan" if not np.isfinite(v) else f"{v:+.4g}{unit}"


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    t0 = time.time()

    # ---- center point (fine delta axis: TPD curve + loss remap + validation)
    print("center solve (fine delta axis)...", flush=True)
    sol_c = _solve(D0_GHZ, T0_C, P_PROBE0, coarse=COARSE_FINE)
    fixed_c, best_c = _readouts(sol_c)
    xi_c = fixed_c["xi"]
    print(f"  v6 center: xi(fixed delta -280 MHz) = {xi_c:.4f} dB, "
          f"gap = {fixed_c['gap']:.3f}")
    print(f"  best-delta readout: xi = {best_c['xi']:.4f} dB @ "
          f"TPD = {best_c['tpd']:.1f} MHz, gap = {best_c['gap']:.3f} "
          f"(v6 report: -8.102 dB @ -280.0 MHz, gap 0.623)")
    # loss remap must reproduce the solved S_dB at the reference loss
    xi_remap_ref = _xi_of_loss(sol_c, fixed_c["idx"], np.array([LOSS0]))[0]
    print(f"  eta-remap check @ loss={LOSS0}: remap {xi_remap_ref:.4f} dB "
          f"vs solve {xi_c:.4f} dB (diff {abs(xi_remap_ref - xi_c):.2e})")

    # ---- re-solve scans (OPD, T, P_probe) in one thread pool ---------------
    jobs = ([("opd", i, float(d), T0_C, P_PROBE0)
             for i, d in enumerate(OPD_AXIS_GHZ)]
            + [("temp", i, D0_GHZ, float(t), P_PROBE0)
               for i, t in enumerate(TEMP_AXIS_C)]
            + [("ppr", i, D0_GHZ, T0_C, float(P_PROBE0 * f))
               for i, f in enumerate(PPROBE_FACTORS)])
    res = {"opd": [None] * OPD_AXIS_GHZ.size,
           "temp": [None] * TEMP_AXIS_C.size,
           "ppr": [None] * PPROBE_FACTORS.size}
    workers = min(os.cpu_count() or 1, 8)
    print(f"\n1D scans: {len(jobs)} hardened-Ultra solves on {workers} workers",
          flush=True)
    with ThreadPoolExecutor(max_workers=workers) as ex:
        futs = {ex.submit(_solve, d, t, p): (kind, i)
                for (kind, i, d, t, p) in jobs}
        for n, fut in enumerate(as_completed(futs), start=1):
            kind, i = futs[fut]
            res[kind][i] = _readouts(fut.result())
            if n % 10 == 0 or n == len(jobs):
                print(f"  {n}/{len(jobs)} ({time.time() - t0:.0f}s)", flush=True)

    def _unpack(kind):
        fixed = np.array([r[0]["xi"] for r in res[kind]])
        fgap = np.array([r[0]["gap"] for r in res[kind]])
        ftrust = np.array([r[0]["trusted"] for r in res[kind]])
        best = np.array([r[1]["xi"] for r in res[kind]])
        btrust = np.array([r[1]["trusted"] for r in res[kind]])
        btpd = np.array([r[1]["tpd"] for r in res[kind]])
        return fixed, fgap, ftrust, best, btrust, btpd

    xi_opd_f, gap_opd, ftr_opd, xi_opd_b, tr_opd, tpd_opd = _unpack("opd")
    xi_t_f, gap_t, ftr_t, xi_t_b, tr_t, tpd_t = _unpack("temp")
    xi_p_f, gap_p, ftr_p, xi_p_b, tr_p, tpd_p = _unpack("ppr")

    # ---- TPD curve (from the fine center solve) + loss remap ---------------
    tpd_axis = sol_c["tpd_mhz"]
    xi_tpd = sol_c["xi"]
    gap_tpd = sol_c["Gs"] - sol_c["Gc"]
    xi_loss = _xi_of_loss(sol_c, fixed_c["idx"], LOSS_AXIS)

    # ---- tolerance table ----------------------------------------------------
    scans = [
        ("OPD [GHz]", OPD_AXIS_GHZ, xi_opd_f, D0_GHZ, " GHz"),
        ("TPD [MHz]", tpd_axis, xi_tpd, TPD0_MHZ, " MHz"),
        ("T [C]", TEMP_AXIS_C, xi_t_f, T0_C, " C"),
        ("P_probe [x ref]", PPROBE_FACTORS, xi_p_f, 1.0, " x"),
        ("loss [frac]", LOSS_AXIS, xi_loss, LOSS0, ""),
    ]
    gate_trust = {"OPD [GHz]": ftr_opd, "TPD [MHz]": (gap_tpd >= GAP_MIN) & (gap_tpd <= GAP_MAX),
                  "T [C]": ftr_t, "P_probe [x ref]": ftr_p, "loss [frac]": None}
    table = {}
    print("\n================  TOLERANCE (fixed-delta readout)  ================")
    print(f"center xi = {xi_c:.3f} dB; thresholds = center + "
          + " / ".join(f"{t}" for t in THRESHOLDS_DB) + " dB")
    for name, x, xi, xc, unit in scans:
        row = {}
        print(f"\n  {name}  (center {xc:+g})")
        for thr in THRESHOLDS_DB:
            lo, hi = _crossings(np.asarray(x, float), np.asarray(xi, float),
                                xc, xi_c, thr)
            row[thr] = (lo, hi)
            dlo = lo - xc if np.isfinite(lo) else np.nan
            dhi = hi - xc if np.isfinite(hi) else np.nan
            print(f"    +{thr:4.2f} dB: [{_fmt(lo, unit)}, {_fmt(hi, unit)}]"
                  f"   (drift {_fmt(dlo, unit)} / {_fmt(dhi, unit)})")
        tr = gate_trust[name]
        if tr is not None:
            x_arr = np.asarray(x, float)
            ic = int(np.argmin(np.abs(x_arr - xc)))
            i_lo = i_hi = ic
            while i_lo > 0 and tr[i_lo - 1]:
                i_lo -= 1
            while i_hi < x_arr.size - 1 and tr[i_hi + 1]:
                i_hi += 1
            print(f"    gap-gate trusted (fixed delta): "
                  f"[{x_arr[i_lo]:+g}, {x_arr[i_hi]:+g}]{unit}")
        table[name] = row

    print("\n================  BEST-DELTA (TPD re-tuned) readout  ============")
    print(f"center best-delta xi = {best_c['xi']:.3f} dB")
    for name, x, xc, xi_b, tr, unit in [
            ("OPD [GHz]", OPD_AXIS_GHZ, D0_GHZ, xi_opd_b, tr_opd, " GHz"),
            ("T [C]", TEMP_AXIS_C, T0_C, xi_t_b, tr_t, " C"),
            ("P_probe [x ref]", PPROBE_FACTORS, 1.0, xi_p_b, tr_p, " x")]:
        print(f"  {name}:")
        for thr in THRESHOLDS_DB:
            lo, hi = _crossings(np.asarray(x, float),
                                np.where(tr, xi_b, np.nan), xc,
                                best_c["xi"], thr)
            print(f"    +{thr:4.2f} dB: [{_fmt(lo, unit)}, {_fmt(hi, unit)}]")

    np.savez(
        OUT / "squeezing_v6_tolerance.npz",
        center=np.array([D0_GHZ, T0_C, TPD0_MHZ, P_PROBE0, LOSS0, xi_c]),
        thresholds_db=np.array(THRESHOLDS_DB),
        opd_ghz=OPD_AXIS_GHZ, xi_opd_fixed=xi_opd_f, xi_opd_best=xi_opd_b,
        gap_opd=gap_opd, trust_opd_fixed=ftr_opd, trust_opd_best=tr_opd,
        tpd_opd_best=tpd_opd,
        temp_c=TEMP_AXIS_C, xi_temp_fixed=xi_t_f, xi_temp_best=xi_t_b,
        gap_temp=gap_t, trust_temp_fixed=ftr_t, trust_temp_best=tr_t,
        tpd_temp_best=tpd_t,
        pprobe_factor=PPROBE_FACTORS, xi_pprobe_fixed=xi_p_f,
        xi_pprobe_best=xi_p_b, gap_pprobe=gap_p, trust_pprobe_fixed=ftr_p,
        trust_pprobe_best=tr_p,
        tpd_mhz=tpd_axis, xi_tpd=xi_tpd, gap_tpd=gap_tpd,
        loss_frac=LOSS_AXIS, xi_loss=xi_loss,
        qe=QE, loss0=LOSS0, p_pump=P_PUMP, p_probe0=P_PROBE0,
        line_strength=LINE_STRENGTH, coarse=COARSE, coarse_fine=COARSE_FINE,
        velocity_step=VELOCITY_STEP, velocity_cutoff=VELOCITY_CUTOFF,
        gap_min=GAP_MIN, gap_max=GAP_MAX,
    )
    print(f"\nsaved {OUT / 'squeezing_v6_tolerance.npz'}")

    _plot(xi_c, scans, table, gap_tpd, gate_trust,
          (xi_opd_b, tr_opd), (xi_t_b, tr_t), (xi_p_b, tr_p))
    print(f"total {time.time() - t0:.0f}s")


def _plot(xi_c, scans, table, gap_tpd, gate_trust, opd_b, t_b, p_b):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    best_curves = {"OPD [GHz]": opd_b, "T [C]": t_b, "P_probe [x ref]": p_b}
    fig, axes = plt.subplots(2, 3, figsize=(16.5, 8.6))
    axs = axes.ravel()
    for ax, (name, x, xi, xc, unit) in zip(axs, scans):
        x = np.asarray(x, float)
        xi = np.asarray(xi, float)
        ax.plot(x, xi, "o-", color="#1f77b4", ms=3, lw=1.5,
                label="fixed delta (pure drift)")
        tr = gate_trust.get(name)
        if tr is not None and not np.all(tr):
            ax.plot(x[~tr], xi[~tr], "x", color="#d62728", ms=5, mew=1.4,
                    label="gap gate FAILED (untrusted)")
        if name in best_curves:
            xb, trb = best_curves[name]
            ax.plot(x, np.where(trb, xb, np.nan), "s--", color="#2ca02c",
                    ms=3, lw=1.1, label="best delta (TPD re-tuned)")
        for thr, color in zip(THRESHOLDS_DB, ("#bbbbbb", "#ff7f0e", "#d62728")):
            ax.axhline(xi_c + thr, color=color, ls="--", lw=0.9)
            lo, hi = table[name][thr]
            if np.isfinite(lo) and np.isfinite(hi) and thr == 0.5:
                ax.axvspan(lo, hi, color="#ff7f0e", alpha=0.10)
        ax.axvline(xc, color="k", ls=":", lw=0.9)
        ax.scatter([xc], [xi_c], color="gold", edgecolor="k", s=70, zorder=5)
        if name.startswith("P_probe"):
            ax.set_xscale("log", base=2)
        if name.startswith("TPD"):
            ax2 = ax.twinx()
            ax2.plot(x, gap_tpd, "-", color="0.6", lw=0.8, alpha=0.7)
            ax2.axhspan(GAP_MIN, GAP_MAX, color="0.85", alpha=0.3)
            ax2.set_ylabel("G_s - G_c (gray)", color="0.45", fontsize=8)
            ax2.tick_params(axis="y", labelcolor="0.45", labelsize=7)
            ax.set_xlim(-450, -100)
        ax.set_xlabel(name)
        ax.set_ylabel("xi [dB]")
        ax.set_title(f"xi vs {name.split(' [')[0]}")
        ax.grid(alpha=0.3)
        ax.legend(fontsize=7, loc="upper right")

    # summary panel: +-0.5 dB drift half-widths as text
    ax = axs[5]
    ax.axis("off")
    lines = [f"v6 optimum: xi = {xi_c:.3f} dB",
             "tolerance windows (fixed-delta drift)", ""]
    for name, x, xi, xc, unit in scans:
        for thr in (0.5,):
            lo, hi = table[name][thr]
            lo_s = f"{lo - xc:+.3g}" if np.isfinite(lo) else "<scan"
            hi_s = f"{hi - xc:+.3g}" if np.isfinite(hi) else ">scan"
            lines.append(f"+{thr} dB  {name:<16s} {lo_s} .. {hi_s}{unit}")
    ax.text(0.02, 0.95, "\n".join(lines), va="top", ha="left", fontsize=9,
            family="monospace", transform=ax.transAxes)

    fig.suptitle("v6 optimum tolerance scan (hardened Ultra, eta = 0.8694)",
                 fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    DOC_OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(DOC_OUT, dpi=150)
    print(f"wrote {DOC_OUT}")


if __name__ == "__main__":
    main()
