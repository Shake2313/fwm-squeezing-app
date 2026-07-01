"""
Extended, high-resolution search for the *theoretical* squeezing maximum.

Differences from scan_squeezing_delta_T.py (the validated/loss=5.5% map):
  * loss = 0           -> detection floor moves to 10log10(1-QE) = -10.2 dB
  * also reports eta=1  -> intrinsic source squeezing S_ideal (no floor)
  * much wider OPD      -> Δ in [-4, +8] GHz
  * higher resolution   -> finer δ scan, finer velocity grid, more points
  * cap sweep           -> how the intrinsic maximum scales with the
                           Manley-Rowe pump-depletion ceiling 1 + P_pump/(2 P_seed)

Key efficiency: G_s, G_c do not depend on the detection efficiency η, and the
δ that minimises S is η-independent (S(η)=η·S_ideal+(1-η) is monotone in
S_ideal). So we run the (Δ,T) grid ONCE, store the deepest-squeezing operating
point's gains, and read off squeezing for any η afterwards.

    python analysis/scan_squeezing_theory_max.py
"""
import os
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", str(os.cpu_count() or 1))

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from gabes.schemes import fwm                       # noqa: E402
from gabes.core import blas_single_thread           # noqa: E402
from gabes import observables as ob                 # noqa: E402

# ---- Fixed reference values (Sim et al.), EXCEPT loss=0 ---------------------
P_PUMP = 0.6
P_SEED = 8e-6
LINE_STRENGTH = 0.74  # residual (Sim et al. anchored); physical p_F/[2(2I+1)] in-engine
QE = fwm.QE_DETECTOR          # 0.9047
BRANCH = -1
WINDOW_GHZ = 0.60             # half-width of the focused δ (probe) scan

# Detection efficiencies we read off the same gains:
ETA_IDEAL = 1.0                       # intrinsic source squeezing (no floor)
ETA_QE = QE                           # lossless detection, QE-limited
ETA_REF = QE * (1.0 - 0.055)          # original 5.5%-loss reference

def floor_dB(eta):
    return 10.0 * np.log10(max(1.0 - eta, 1e-300))

# ---- High-resolution scan grid ---------------------------------------------
DELTA_GHZ = np.round(np.arange(-4.0, 8.0 + 1e-9, 0.2), 3)     # wide OPD
TEMP_C = np.round(np.arange(70.0, 185.0 + 1e-9, 5.0), 3)      # up to 185 C
COARSE_POINTS = 201          # δ resolution ~6 MHz over the window
VELOCITY_STEP = 3.0          # finer Doppler grid


def deepest_point(D_GHz, T_K, coarse=COARSE_POINTS, vstep=VELOCITY_STEP,
                  P_pump=P_PUMP, P_seed=P_SEED, ls=LINE_STRENGTH):
    """At the δ that gives the *most* squeezing (min S_ideal), return the gains.

    Returns the operating point that minimises S_ideal (= maximises balanced
    gain) on the (-) branch; squeezing for any η follows from (G_s, G_c).
    """
    center = fwm.branch_center_GHz(D_GHz, BRANCH)
    with blas_single_thread():
        spec = fwm.compute_spectrum(
            D_GHz, T=T_K, P_pump=P_pump, P_probe=P_seed,
            line_strength=ls, loss_frac=0.0, qe=QE,
            coarse_points=coarse, fine_points=0,
            scan_min=center - WINDOW_GHZ, scan_max=center + WINDOW_GHZ,
            velocity_step=vstep, velocity_cutoff=3.0, branch=BRANCH)
    Gs, Gc = spec["G_s"], spec["G_c"]
    # S_ideal per δ (eta=1) -> argmin is the deepest squeezing point.
    s_ideal_dB = ob.intensity_difference_squeezing_dB(Gs, Gc, 1.0)
    i = int(np.nanargmin(s_ideal_dB))
    n = Gs.size
    return dict(
        S_ideal_dB=float(s_ideal_dB[i]),
        G_s=float(Gs[i]), G_c=float(Gc[i]),
        delta_mhz=float((spec["probe_axis_GHz"][i] - center) * 1e3),
        at_edge=bool(i <= 1 or i >= n - 2),
        smallsignal_peak=float(spec["G_s_smallsignal_peak"]),
        cap=float(spec["pump_depletion_cap"]),
    )


def run_grid():
    jobs = [(iD, iT, float(D), float(T) + 273.15)
            for iD, D in enumerate(DELTA_GHZ)
            for iT, T in enumerate(TEMP_C)]
    nD, nT = DELTA_GHZ.size, TEMP_C.size
    Sid = np.full((nD, nT), np.nan)
    Gs = np.full((nD, nT), np.nan)
    Gc = np.full((nD, nT), np.nan)
    Dlt = np.full((nD, nT), np.nan)
    Edge = np.zeros((nD, nT), dtype=bool)

    t0 = time.time()
    workers = min(os.cpu_count() or 1, 8)
    print(f"  {len(jobs)} points on {workers} workers ...", flush=True)
    with ThreadPoolExecutor(max_workers=workers) as ex:
        futs = {ex.submit(deepest_point, D, T): (iD, iT)
                for (iD, iT, D, T) in jobs}
        done = 0
        for fut in as_completed(futs):
            iD, iT = futs[fut]
            r = fut.result()
            Sid[iD, iT] = r["S_ideal_dB"]
            Gs[iD, iT] = r["G_s"]
            Gc[iD, iT] = r["G_c"]
            Dlt[iD, iT] = r["delta_mhz"]
            Edge[iD, iT] = r["at_edge"]
            done += 1
            if done % 100 == 0 or done == len(jobs):
                print(f"    {done}/{len(jobs)}  ({time.time()-t0:.0f}s)", flush=True)
    return Sid, Gs, Gc, Dlt, Edge


def squeezing_dB(Gs, Gc, eta):
    """Vectorised S_dB over the whole map for one η."""
    flat = ob.intensity_difference_squeezing_dB(
        np.asarray(Gs).ravel(), np.asarray(Gc).ravel(), eta)
    return flat.reshape(np.asarray(Gs).shape)


def report_eta(name, eta, Gs, Gc, Dlt):
    S = squeezing_dB(Gs, Gc, eta)
    flat = int(np.nanargmin(S))
    iD, iT = np.unravel_index(flat, S.shape)
    fl = floor_dB(eta)
    print(f"\n--- η = {eta:.4f}  ({name}) ---")
    print(f"  floor 10log10(1-η)      = {fl:8.3f} dB"
          if eta < 1 else "  floor                   = -inf (no detection floor)")
    print(f"  best squeezing ξ_max    = {S[iD,iT]:8.3f} dB")
    print(f"    at Δ = {DELTA_GHZ[iD]:+.2f} GHz, T = {TEMP_C[iT]:.0f} °C, "
          f"δ = {Dlt[iD,iT]:.1f} MHz")
    print(f"    G_s = {Gs[iD,iT]:.3e}, G_c = {Gc[iD,iT]:.3e}")
    if eta < 1:
        n_floor = int(np.sum(S < fl + 0.01))
        print(f"  points within 0.01 dB of floor: {n_floor}/{S.size}")
        # efficient frontier: lowest T that reaches the floor, and its Δ
        reached = S < fl + 0.05
        if reached.any():
            iTs = np.where(reached.any(axis=0))[0]
            jT = iTs.min()
            iD_best = int(np.nanargmin(S[:, jT]))
            print(f"  efficient frontier: floor first reached at "
                  f"T = {TEMP_C[jT]:.0f} °C, Δ = {DELTA_GHZ[iD_best]:+.2f} GHz "
                  f"(G_s = {Gs[iD_best,jT]:.1f})")
    return S


def cap_sweep():
    """Intrinsic (η=1) deepest squeezing vs the pump-depletion ceiling.

    cap = 1 + P_pump/(2 P_seed). Pushed to a high T so the small-signal gain
    saturates against the cap. The conserved-(N_s−N_c) ideal squeezing is
    S_ideal = (G_s−G_c)²/(G_s+G_c); in the LOSSLESS limit the twin-beam relation
    G_c = G_s−1 holds and this reduces to the reference law S_ideal ≈ 1/(2G−1).

    Caveat (collisional decoherence, `fwm.collisional_atom`): at the high T needed
    to saturate the cap, Rb self-broadening adds real in-cell absorption, so the
    propagation no longer preserves G_s−G_c = 1 (the twin gains are compressed
    toward each other, gap < 1). The lossless formula then reports S_ideal *below*
    the 1/(2G−1) reference — this is the ideal expression applied outside its
    lossless domain, not extra physical squeezing. A faithful lossy source needs
    the quantum-Langevin noise channel (docs/checklist.json
    `fwm-quantum-langevin-noise`); until then treat the η=1 column as a lossless
    idealization and the `gap` as the diagnostic of how far it is trusted.
    """
    print("\n================  CAP SWEEP (η=1, intrinsic)  ================")
    print("  using Δ = +2.0 GHz, T = 180 °C (high gain), loss = 0")
    print(f"  {'P_pump[W]':>9} {'P_seed[uW]':>10} {'cap':>12} {'G_s':>12} "
          f"{'gap Gs-Gc':>10} {'S_ideal[dB]':>12} {'1/(2G-1)[dB]':>13}")
    combos = [
        (0.6, 8e-6), (0.6, 1e-6), (1.2, 1e-6),
        (1.2, 1e-7), (2.0, 1e-8),
    ]
    rows = []
    for P_pump, P_seed in combos:
        r = deepest_point(2.0, 180.0 + 273.15, coarse=161, vstep=3.0,
                          P_pump=P_pump, P_seed=P_seed)
        cap = 1.0 + P_pump / (2.0 * P_seed)
        gap = r["G_s"] - r["G_c"]
        pred = 10.0 * np.log10(1.0 / (2.0 * r["G_s"] - 1.0))
        print(f"  {P_pump:9.2f} {P_seed*1e6:10.4g} {cap:12.3e} "
              f"{r['G_s']:12.3e} {gap:10.3f} {r['S_ideal_dB']:12.3f} {pred:13.3f}")
        rows.append((cap, r["G_s"], r["S_ideal_dB"]))
    print("  (gap Gs-Gc → 1 in the lossless limit; gap < 1 flags collisional "
          "absorption\n   breaking the twin-beam relation, so S_ideal runs below "
          "the 1/(2G-1) reference.)")
    return rows


def main():
    print("Extended theoretical-maximum search (loss = 0)")
    print(f"  QE = {QE},  floor(η=QE) = {floor_dB(ETA_QE):.3f} dB")
    print(f"  grid: {DELTA_GHZ.size} OPD ({DELTA_GHZ[0]}..{DELTA_GHZ[-1]} GHz) "
          f"x {TEMP_C.size} T ({TEMP_C[0]}..{TEMP_C[-1]} °C) "
          f"= {DELTA_GHZ.size*TEMP_C.size} points")
    print(f"  δ window ±{WINDOW_GHZ} GHz, coarse={COARSE_POINTS}, "
          f"velocity_step={VELOCITY_STEP} m/s\n")

    Sid, Gs, Gc, Dlt, Edge = run_grid()

    print(f"\n[edge check] deepest-δ hit scan-window edge at "
          f"{int(Edge.sum())}/{Edge.size} points "
          f"(max |δ| = {np.nanmax(np.abs(Dlt)):.0f} MHz)")

    S_qe = report_eta("lossless, QE-limited", ETA_QE, Gs, Gc, Dlt)
    S_ref = report_eta("original 5.5% loss", ETA_REF, Gs, Gc, Dlt)
    # intrinsic (η=1)
    flat = int(np.nanargmin(Sid))
    iD, iT = np.unravel_index(flat, Sid.shape)
    print(f"\n--- η = 1 (intrinsic source S_ideal, no floor) ---")
    print(f"  deepest S_ideal = {Sid[iD,iT]:.3f} dB at "
          f"Δ = {DELTA_GHZ[iD]:+.2f} GHz, T = {TEMP_C[iT]:.0f} °C")
    print(f"    G_s = {Gs[iD,iT]:.3e} (cap-limited); "
          f"S_ideal scales as 1/(2G-1) -> no interior maximum, only the gain "
          f"ceiling bounds it.")

    rows = cap_sweep()

    out = ROOT / "analysis"
    np.savez(out / "squeezing_theory_max.npz",
             delta_ghz=DELTA_GHZ, temp_c=TEMP_C, S_ideal_dB=Sid,
             G_s=Gs, G_c=Gc, delta_mhz=Dlt,
             S_qe_dB=S_qe, S_ref_dB=S_ref,
             eta_qe=ETA_QE, eta_ref=ETA_REF,
             cap_sweep=np.array(rows))
    _plot(Sid, S_qe, Gs, out)
    print(f"\nsaved: {out/'squeezing_theory_max.npz'} and squeezing_theory_max.png")


def _plot(Sid, S_qe, Gs, out):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 3, figsize=(17, 4.8))
    extent = [TEMP_C[0], TEMP_C[-1], DELTA_GHZ[0], DELTA_GHZ[-1]]

    # (1) lossless QE-limited squeezing map
    im0 = axes[0].imshow(S_qe, origin="lower", aspect="auto", extent=extent,
                         cmap="viridis_r")
    fig.colorbar(im0, ax=axes[0], label="ξ [dB]")
    axes[0].axhline(0.9, color="white", ls=":", lw=0.8)
    axes[0].set_title(f"Lossless (η=QE): floor {floor_dB(ETA_QE):.1f} dB")
    axes[0].set_xlabel("T [°C]"); axes[0].set_ylabel("Δ [GHz]")

    # (2) intrinsic S_ideal map (clipped for visibility)
    Sid_clip = np.clip(Sid, -40, 0)
    im1 = axes[1].imshow(Sid_clip, origin="lower", aspect="auto", extent=extent,
                         cmap="magma_r")
    fig.colorbar(im1, ax=axes[1], label="S_ideal [dB] (clip −40)")
    axes[1].set_title("Intrinsic source (η=1): no floor → gain-limited")
    axes[1].set_xlabel("T [°C]"); axes[1].set_ylabel("Δ [GHz]")

    # (3) log10 gain map
    im2 = axes[2].imshow(np.log10(np.clip(Gs, 1e-3, None)), origin="lower",
                         aspect="auto", extent=extent, cmap="cividis")
    fig.colorbar(im2, ax=axes[2], label="log10 G_s")
    axes[2].set_title("Balanced gain log10 G_s (cap-clamped at high T)")
    axes[2].set_xlabel("T [°C]"); axes[2].set_ylabel("Δ [GHz]")

    fig.tight_layout()
    fig.savefig(out / "squeezing_theory_max.png", dpi=130)


if __name__ == "__main__":
    main()
