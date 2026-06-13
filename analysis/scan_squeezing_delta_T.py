"""
Maximise intensity-difference squeezing ξ over the (Δ, T) plane.

Δ  = one-photon detuning (OPD, pump position relative to F=2 -> F'=3) [GHz]
T  = cell temperature [degC]
ξ  = best (most negative) twin-beam squeezing [dB] on the (-) Raman branch,
     obtained by scanning the two-photon detuning δ (probe-frequency scan)
     for each (Δ, T) and taking the minimum S_dB.

Everything except Δ and T is fixed at the Sim et al. reference operating point
(G. Sim, H. Kim, H. S. Moon, Sci. Rep. 15, 7727 (2025)) -- the same numbers the
regression test calls "sim_optimum":

    P_pump = 600 mW,  P_seed = 8 uW,  line_strength = 0.05,  loss = 5.5 %,
    QE = 0.9047,  L_cell = 12.5 mm,  w_pump = 530 um,  w_seed = 330 um.

Pump / seed power and the line-strength calibration do not change the qualitative
squeezing result (gain saturates against the pump-depletion bound); the detection
efficiency η = QE·(1-loss) does, because it sets the hard squeezing floor

    S_floor(dB) = 10·log10(1 - η).

    python analysis/scan_squeezing_delta_T.py
"""
import os
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", str(os.cpu_count() or 1))

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from gabes import kernels  # noqa: E402
from gabes.schemes import fwm  # noqa: E402
from gabes.core import blas_single_thread  # noqa: E402

# ---- Fixed Sim et al. reference values (everything except Δ, T) -------------
P_PUMP = 0.6        # W   (600 mW)  -- adjustable, gain saturates anyway
P_SEED = 8e-6       # W   (8 uW)
LINE_STRENGTH = 0.05
LOSS_FRAC = 0.055   # 5.5 % loss after the cell
QE = fwm.QE_DETECTOR          # 0.9047
ETA = QE * (1.0 - LOSS_FRAC)
BRANCH = -1                   # single Raman branch (the (-) line)
WINDOW_GHZ = 0.55             # half-width of the focused δ (probe) scan window

# η-limited theoretical squeezing floor: S(η)=η·S_ideal+(1-η) >= 1-η.
S_FLOOR_DB = 10.0 * np.log10(1.0 - ETA)

# ---- Scan grid --------------------------------------------------------------
DELTA_GHZ = np.round(np.arange(-0.5, 3.0 + 1e-9, 0.1), 3)      # OPD
TEMP_C = np.round(np.arange(60.0, 150.0 + 1e-9, 5.0), 3)       # temperature

# Map resolution (fast but converged for the *minimum* of a smooth curve).
COARSE_POINTS = 81
VELOCITY_STEP = 5.0


def best_squeezing(D_GHz, T_K, coarse=COARSE_POINTS, vstep=VELOCITY_STEP):
    """Min S_dB over the δ scan on the (-) branch, plus the gains there."""
    center = fwm.branch_center_GHz(D_GHz, BRANCH)
    with blas_single_thread():
        spec = fwm.compute_spectrum(
            D_GHz, T=T_K, P_pump=P_PUMP, P_probe=P_SEED,
            line_strength=LINE_STRENGTH, loss_frac=LOSS_FRAC, qe=QE,
            coarse_points=coarse, fine_points=0,
            scan_min=center - WINDOW_GHZ, scan_max=center + WINDOW_GHZ,
            velocity_step=vstep, velocity_cutoff=3.0, branch=BRANCH)
    S = spec["S_dB"]
    i = int(np.nanargmin(S))
    delta_mhz = (spec["probe_axis_GHz"][i] - center) * 1e3
    return dict(
        S_dB=float(S[i]),
        delta_mhz=float(delta_mhz),
        G_s=float(spec["G_s"][i]),
        G_c=float(spec["G_c"][i]),
        G_s_smallsignal_peak=float(spec["G_s_smallsignal_peak"]),
        pump_depletion_cap=float(spec["pump_depletion_cap"]),
    )


def run_grid():
    jobs = [(iD, iT, float(D), float(T) + 273.15)
            for iD, D in enumerate(DELTA_GHZ)
            for iT, T in enumerate(TEMP_C)]
    nD, nT = DELTA_GHZ.size, TEMP_C.size
    Xi = np.full((nD, nT), np.nan)
    Gs = np.full((nD, nT), np.nan)
    Gc = np.full((nD, nT), np.nan)
    Dlt = np.full((nD, nT), np.nan)

    t0 = time.time()
    # The numba kernel is internally parallel (saturates all cores per point),
    # so thread the grid only on the NumPy fallback path.
    workers = 1 if kernels.available() else min(os.cpu_count() or 1, 8)
    with ThreadPoolExecutor(max_workers=workers) as ex:
        futs = {ex.submit(best_squeezing, D, T): (iD, iT)
                for (iD, iT, D, T) in jobs}
        done = 0
        from concurrent.futures import as_completed
        for fut in as_completed(futs):
            iD, iT = futs[fut]
            r = fut.result()
            Xi[iD, iT] = r["S_dB"]
            Gs[iD, iT] = r["G_s"]
            Gc[iD, iT] = r["G_c"]
            Dlt[iD, iT] = r["delta_mhz"]
            done += 1
            if done % 50 == 0 or done == len(jobs):
                print(f"  {done}/{len(jobs)} points "
                      f"({time.time()-t0:.0f}s)", flush=True)
    return Xi, Gs, Gc, Dlt


def main():
    print(f"eta = QE*(1-loss) = {QE:.4f}*(1-{LOSS_FRAC}) = {ETA:.5f}")
    print(f"theoretical squeezing floor 10log10(1-eta) = {S_FLOOR_DB:.4f} dB")
    print(f"grid: {DELTA_GHZ.size} OPD x {TEMP_C.size} T = "
          f"{DELTA_GHZ.size*TEMP_C.size} points\n")

    Xi, Gs, Gc, Dlt = run_grid()

    # Global maximum squeezing (most negative S_dB).
    flat = int(np.nanargmin(Xi))
    iD, iT = np.unravel_index(flat, Xi.shape)
    Dbest, Tbest = DELTA_GHZ[iD], TEMP_C[iT]
    Sbest = Xi[iD, iT]

    print("\n================  RESULT  ================")
    print(f"Best squeezing  ξ_max = {Sbest:.4f} dB")
    print(f"  at  Δ = {Dbest:.2f} GHz,  T = {Tbest:.0f} °C")
    print(f"  G_s = {Gs[iD,iT]:.2f},  G_c = {Gc[iD,iT]:.2f},  "
          f"δ = {Dlt[iD,iT]:.1f} MHz")
    print(f"Theoretical floor      = {S_FLOOR_DB:.4f} dB")
    print(f"Gap to floor           = {Sbest - S_FLOOR_DB:+.4f} dB")
    below_floor = Xi[np.isfinite(Xi)] < S_FLOOR_DB - 1e-6
    print(f"points below the η-floor (should be 0): {int(below_floor.sum())}")

    # Sim reference point for context.
    sim = best_squeezing(0.9, 394.15, coarse=121, vstep=4.0)
    print(f"\nSim optimum (Δ=0.9 GHz, T=121 °C): "
          f"ξ = {sim['S_dB']:.4f} dB,  G_s = {sim['G_s']:.2f}")

    out = ROOT / "analysis"
    np.savez(out / "squeezing_map.npz",
             delta_ghz=DELTA_GHZ, temp_c=TEMP_C, xi_dB=Xi,
             G_s=Gs, G_c=Gc, delta_mhz=Dlt, eta=ETA, floor_dB=S_FLOOR_DB)

    _plot(Xi, Gs, Dbest, Tbest, Sbest, out)
    print(f"\nsaved: {out/'squeezing_map.npz'}  and  squeezing_map.png")


def _plot(Xi, Gs, Dbest, Tbest, Sbest, out):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, (axM, axT, axD) = plt.subplots(1, 3, figsize=(16, 4.6))

    extent = [TEMP_C[0], TEMP_C[-1], DELTA_GHZ[0], DELTA_GHZ[-1]]
    im = axM.imshow(Xi, origin="lower", aspect="auto", extent=extent,
                    cmap="viridis_r")
    fig.colorbar(im, ax=axM, label="best squeezing ξ  [dB]")
    axM.contour(TEMP_C, DELTA_GHZ, Xi,
                levels=[S_FLOOR_DB + 0.2, S_FLOOR_DB + 0.5, S_FLOOR_DB + 1.0],
                colors="white", linewidths=0.7, alpha=0.6)
    axM.scatter([Tbest], [Dbest], color="crimson", marker="*", s=160,
                edgecolor="k", zorder=5, label=f"ξ_max = {Sbest:.2f} dB")
    axM.scatter([121], [0.9], color="white", marker="o", s=40,
                edgecolor="k", zorder=5, label="Sim optimum")
    axM.set_xlabel("Temperature T  [°C]")
    axM.set_ylabel("One-photon detuning Δ  [GHz]")
    axM.set_title("ξ(Δ, T)  (min over δ, − branch)")
    axM.legend(loc="upper left", fontsize=8)

    # Slice vs T at the best Δ -> approach to the η floor.
    iD = int(np.argmin(np.abs(DELTA_GHZ - Dbest)))
    axT.plot(TEMP_C, Xi[iD, :], "o-", color="#2ca02c", lw=1.6, ms=3)
    axT.axhline(S_FLOOR_DB, color="crimson", ls="--",
                label=f"floor {S_FLOOR_DB:.2f} dB")
    axT.set_xlabel("Temperature T  [°C]")
    axT.set_ylabel("ξ  [dB]")
    axT.set_title(f"ξ vs T  (Δ = {Dbest:.2f} GHz)")
    axT.grid(alpha=0.3)
    axT.legend(fontsize=8)

    # Slice vs Δ at the best T.
    iT = int(np.argmin(np.abs(TEMP_C - Tbest)))
    axD.plot(DELTA_GHZ, Xi[:, iT], "o-", color="#1f77b4", lw=1.6, ms=3)
    axD.axhline(S_FLOOR_DB, color="crimson", ls="--",
                label=f"floor {S_FLOOR_DB:.2f} dB")
    axD.set_xlabel("One-photon detuning Δ  [GHz]")
    axD.set_ylabel("ξ  [dB]")
    axD.set_title(f"ξ vs Δ  (T = {Tbest:.0f} °C)")
    axD.grid(alpha=0.3)
    axD.legend(fontsize=8)

    fig.tight_layout()
    fig.savefig(out / "squeezing_map.png", dpi=130)


if __name__ == "__main__":
    main()
