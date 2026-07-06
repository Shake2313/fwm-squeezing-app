"""
Detailed v6 ideal-detection scan around the beat-sign-corrected loss=0 optimum.

This extends the existing (Delta, T) ideal scan by treating two-photon detuning
(TPD, the probe-axis offset from the branch center) and pump-probe beam angle as
explicit scan variables. The engine is unchanged: each (Delta, T, angle) point
uses the hardened Ultra path, and the returned TPD axis is sampled directly.

Outputs:
  analysis/squeezing_frontier_ideal_v6_tpd_angle_detail/squeezing_tpd_angle_detail.npz
  docs/squeezing_report/squeezing_frontier_ideal_v6_tpd_angle_detail.png

Run from the repo root:
  python analysis/scan_squeezing_v6_loss0_tpd_angle_detail.py
"""
from __future__ import annotations

import os
import sys
import time
import warnings
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from gabes.core import blas_single_thread  # noqa: E402
from gabes.schemes import fwm  # noqa: E402


OUT = ROOT / "analysis" / "squeezing_frontier_ideal_v6_tpd_angle_detail"
DOC_OUT = ROOT / "docs" / "squeezing_report" / "squeezing_frontier_ideal_v6_tpd_angle_detail.png"

P_PUMP = 0.6
P_SEED = 8e-6
LINE_STRENGTH = 0.74
BRANCH = -1
QE = 1.0
LOSS = 0.0

# Local ranges around the v6 loss=0 optimum after fixing the seeded sideband
# beat sign: Delta=-1.50 GHz, T=110 C, TPD about -280 MHz at the default
# 0.32 deg geometry. Keep the lower-T/low-angle side open enough to diagnose
# whether the geometry ridge persists.
D_AXIS_GHZ = np.round(np.arange(-1.56, -1.44 + 1e-9, 0.02), 3)
TEMP_AXIS_C = np.round(np.arange(96.0, 114.0 + 1e-9, 1.0), 3)
ANGLE_AXIS_DEG = np.round(np.arange(0.00, 0.36 + 1e-9, 0.02), 3)
TPD_AXIS_MHZ = np.round(np.linspace(-420.0, -160.0, 81), 3)

VELOCITY_STEP = 4.0
VELOCITY_CUTOFF = 3.0

GAP_MIN = 0.5
GAP_MAX = 1.5


def _compute_one(iD: int, iT: int, iA: int) -> dict:
    D = float(D_AXIS_GHZ[iD])
    T_C = float(TEMP_AXIS_C[iT])
    angle = float(ANGLE_AXIS_DEG[iA])
    center = fwm.branch_center_GHz(D, BRANCH)
    scan_min = center + float(TPD_AXIS_MHZ[0]) * 1e-3
    scan_max = center + float(TPD_AXIS_MHZ[-1]) * 1e-3

    with blas_single_thread():
        spec = fwm.compute_spectrum(
            D,
            T=T_C + 273.15,
            P_pump=P_PUMP,
            P_probe=P_SEED,
            line_strength=LINE_STRENGTH,
            loss_frac=LOSS,
            qe=QE,
            coarse_points=TPD_AXIS_MHZ.size,
            fine_points=0,
            scan_min=scan_min,
            scan_max=scan_max,
            velocity_step=VELOCITY_STEP,
            velocity_cutoff=VELOCITY_CUTOFF,
            branch=BRANCH,
            phase_detail=fwm.PHASE_ULTRA,
            pump_probe_angle_deg=angle,
        )

    tpd = (spec["probe_axis_GHz"] - center) * 1e3
    if not np.allclose(tpd, TPD_AXIS_MHZ, atol=1e-9, rtol=0.0):
        raise RuntimeError("TPD axis mismatch")

    hardened = spec["hardened_noise"] or {}
    return {
        "iD": iD,
        "iT": iT,
        "iA": iA,
        "xi": np.asarray(spec["S_dB"], dtype=float),
        "Gs": np.asarray(spec["G_s"], dtype=float),
        "Gc": np.asarray(spec["G_c"], dtype=float),
        "delta_k_z": np.asarray(spec["delta_k_z"], dtype=float),
        "od_conj": np.asarray(hardened.get("od_conj_arr", np.full(TPD_AXIS_MHZ.size, np.nan)), dtype=float),
        "od_probe": np.asarray(hardened.get("od_probe_lin_arr", np.full(TPD_AXIS_MHZ.size, np.nan)), dtype=float),
        "segment_od": float(spec.get("segment_absorption_od", np.nan)),
        "pump_scatter": float(hardened.get("pump_scatter_noise", np.nan)),
        "od_pump": float(hardened.get("od_pump", np.nan)),
        "spatial_overlap_min": float(spec.get("ultra_spatial_overlap_min", np.nan)),
        "phase_max_change": float(spec.get("ultra_phase_max_change", np.nan)),
    }


def _best(masked_xi: np.ndarray) -> tuple[int, ...]:
    flat = int(np.nanargmin(masked_xi))
    return np.unravel_index(flat, masked_xi.shape)


def _nanmin(a: np.ndarray, axis=None) -> np.ndarray:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        return np.nanmin(a, axis=axis)


def _row(label: str, idx: tuple[int, int, int, int], xi, Gs, Gc, od_conj, od_probe, delta_k_z) -> dict:
    iD, iT, iA, iP = idx
    return {
        "label": label,
        "Delta_GHz": float(D_AXIS_GHZ[iD]),
        "T_C": float(TEMP_AXIS_C[iT]),
        "angle_deg": float(ANGLE_AXIS_DEG[iA]),
        "TPD_MHz": float(TPD_AXIS_MHZ[iP]),
        "xi_dB": float(xi[idx]),
        "G_s": float(Gs[idx]),
        "G_c": float(Gc[idx]),
        "gap": float(Gs[idx] - Gc[idx]),
        "od_conj": float(od_conj[idx]),
        "od_probe": float(od_probe[idx]),
        "delta_k_z": float(delta_k_z[idx]),
    }


def _print_row(row: dict) -> None:
    print(
        f"{row['label']}: xi={row['xi_dB']:.4f} dB, "
        f"Delta={row['Delta_GHz']:+.3f} GHz, T={row['T_C']:.1f} C, "
        f"TPD={row['TPD_MHz']:+.1f} MHz, angle={row['angle_deg']:.3f} deg, "
        f"Gs={row['G_s']:.2f}, gap={row['gap']:.3f}, "
        f"od_c={row['od_conj']:.4f}, od_p={row['od_probe']:.4f}, "
        f"dkz={row['delta_k_z']:.2e} 1/m"
    )


def _plot(xi, trusted, summary_rows) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    txi = np.where(trusted, xi, np.nan)
    angle_tpd = _nanmin(txi, axis=(0, 1))  # angle x TPD
    dt_map = _nanmin(txi, axis=(2, 3))  # Delta x T
    best_vs_angle = _nanmin(angle_tpd, axis=1)
    best_tpd_idx = np.nanargmin(angle_tpd, axis=1)

    best = summary_rows[0]
    default = summary_rows[1]

    fig, axes = plt.subplots(1, 3, figsize=(17.5, 4.8))
    cmap_angle = plt.get_cmap("magma_r").copy()
    cmap_angle.set_bad("#eeeeee")
    cmap_dt = plt.get_cmap("viridis_r").copy()
    cmap_dt.set_bad("#eeeeee")

    im0 = axes[0].imshow(
        angle_tpd.T,
        origin="lower",
        aspect="auto",
        extent=[ANGLE_AXIS_DEG[0], ANGLE_AXIS_DEG[-1], TPD_AXIS_MHZ[0], TPD_AXIS_MHZ[-1]],
        cmap=cmap_angle,
        vmin=np.nanmin(angle_tpd),
        vmax=min(-13.0, np.nanmax(angle_tpd)),
    )
    fig.colorbar(im0, ax=axes[0], label="best trusted xi [dB]")
    axes[0].scatter([best["angle_deg"]], [best["TPD_MHz"]], s=70, c="cyan", edgecolor="k", zorder=5)
    axes[0].scatter([default["angle_deg"]], [default["TPD_MHz"]], s=55, c="white", edgecolor="k", zorder=5)
    axes[0].axvline(fwm.SEEDED_PHASE_ANGLE_DEG, color="white", ls=":", lw=1.0)
    axes[0].set_xlabel("pump-probe angle [deg]")
    axes[0].set_ylabel("TPD [MHz]")
    axes[0].set_title("min over Delta,T: new knobs")
    axes[0].text(0.01, 0.02, "gray: gap/edge rejected", transform=axes[0].transAxes,
                 fontsize=7, color="0.25")

    axes[1].plot(ANGLE_AXIS_DEG, best_vs_angle, "o-", color="#1f77b4", ms=3, lw=1.6)
    axes[1].axvline(fwm.SEEDED_PHASE_ANGLE_DEG, color="0.25", ls=":", lw=1.0, label="default 0.32 deg")
    axes[1].scatter([best["angle_deg"]], [best["xi_dB"]], c="cyan", edgecolor="k", zorder=5, label="scan best")
    axes[1].scatter([default["angle_deg"]], [default["xi_dB"]], c="white", edgecolor="k", zorder=5, label="default-angle best")
    for a, tpd in zip(ANGLE_AXIS_DEG, TPD_AXIS_MHZ[best_tpd_idx]):
        if int(round(a * 100)) % 4 == 0:
            axes[1].annotate(f"{tpd:.0f}", (a, best_vs_angle[int(np.argmin(np.abs(ANGLE_AXIS_DEG - a)))]),
                             textcoords="offset points", xytext=(0, 6), ha="center", fontsize=6)
    axes[1].set_xlabel("pump-probe angle [deg]")
    axes[1].set_ylabel("best trusted xi [dB]")
    axes[1].set_title("angle sensitivity (labels: best TPD MHz)")
    axes[1].grid(alpha=0.3)
    axes[1].legend(fontsize=7)

    im2 = axes[2].imshow(
        dt_map,
        origin="lower",
        aspect="auto",
        extent=[TEMP_AXIS_C[0], TEMP_AXIS_C[-1], D_AXIS_GHZ[0], D_AXIS_GHZ[-1]],
        cmap=cmap_dt,
    )
    fig.colorbar(im2, ax=axes[2], label="best trusted xi [dB]")
    axes[2].scatter([best["T_C"]], [best["Delta_GHz"]], c="cyan", edgecolor="k", s=70, zorder=5)
    axes[2].scatter([110.0], [-1.5], c="white", edgecolor="k", s=45, zorder=5, label="coarse v6 ideal")
    axes[2].set_xlabel("T [C]")
    axes[2].set_ylabel("Delta [GHz]")
    axes[2].set_title("min over TPD,angle")
    axes[2].legend(fontsize=7)

    fig.tight_layout()
    DOC_OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(DOC_OUT, dpi=150)
    print(f"wrote {DOC_OUT}")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    shape4 = (D_AXIS_GHZ.size, TEMP_AXIS_C.size, ANGLE_AXIS_DEG.size, TPD_AXIS_MHZ.size)
    shape3 = shape4[:3]
    xi = np.full(shape4, np.nan)
    Gs = np.full(shape4, np.nan)
    Gc = np.full(shape4, np.nan)
    delta_k_z = np.full(shape4, np.nan)
    od_conj = np.full(shape4, np.nan)
    od_probe = np.full(shape4, np.nan)
    segment_od = np.full(shape3, np.nan)
    pump_scatter = np.full(shape3, np.nan)
    od_pump = np.full(shape3, np.nan)
    spatial_overlap_min = np.full(shape3, np.nan)
    phase_max_change = np.full(shape3, np.nan)

    jobs = [(iD, iT, iA)
            for iD in range(D_AXIS_GHZ.size)
            for iT in range(TEMP_AXIS_C.size)
            for iA in range(ANGLE_AXIS_DEG.size)]
    workers = min(os.cpu_count() or 1, 8)
    print("Detailed v6 loss=0 TPD/angle scan")
    print(f"  grid: Delta {D_AXIS_GHZ.size} x T {TEMP_AXIS_C.size} x angle {ANGLE_AXIS_DEG.size} "
          f"x TPD {TPD_AXIS_MHZ.size} = {np.prod(shape4)} sampled points")
    print(f"  compute calls: {len(jobs)} on {workers} workers")
    print(f"  gap trust gate: {GAP_MIN} <= Gs-Gc <= {GAP_MAX}, TPD-edge excluded")

    t0 = time.time()
    with ThreadPoolExecutor(max_workers=workers) as ex:
        futs = {ex.submit(_compute_one, *job): job for job in jobs}
        for n, fut in enumerate(as_completed(futs), start=1):
            r = fut.result()
            iD, iT, iA = r["iD"], r["iT"], r["iA"]
            xi[iD, iT, iA, :] = r["xi"]
            Gs[iD, iT, iA, :] = r["Gs"]
            Gc[iD, iT, iA, :] = r["Gc"]
            delta_k_z[iD, iT, iA, :] = r["delta_k_z"]
            od_conj[iD, iT, iA, :] = r["od_conj"]
            od_probe[iD, iT, iA, :] = r["od_probe"]
            segment_od[iD, iT, iA] = r["segment_od"]
            pump_scatter[iD, iT, iA] = r["pump_scatter"]
            od_pump[iD, iT, iA] = r["od_pump"]
            spatial_overlap_min[iD, iT, iA] = r["spatial_overlap_min"]
            phase_max_change[iD, iT, iA] = r["phase_max_change"]
            if n % 50 == 0 or n == len(jobs):
                print(f"  {n}/{len(jobs)} calls done ({time.time() - t0:.1f}s)", flush=True)

    gap = Gs - Gc
    tpd_edge = np.zeros(shape4, dtype=bool)
    tpd_edge[..., 0] = True
    tpd_edge[..., -1] = True
    gap_bad = (gap < GAP_MIN) | (gap > GAP_MAX)
    trusted = np.isfinite(xi) & ~gap_bad & ~tpd_edge
    trusted_xi = np.where(trusted, xi, np.nan)

    best_idx = _best(trusted_xi)
    default_iA = int(np.argmin(np.abs(ANGLE_AXIS_DEG - fwm.SEEDED_PHASE_ANGLE_DEG)))
    default_xi = np.full_like(trusted_xi, np.nan)
    default_xi[:, :, default_iA, :] = trusted_xi[:, :, default_iA, :]
    default_idx = _best(default_xi)

    raw_idx = _best(np.where(np.isfinite(xi) & ~tpd_edge, xi, np.nan))
    rows = [
        _row("trusted 4D best", best_idx, xi, Gs, Gc, od_conj, od_probe, delta_k_z),
        _row("trusted default-angle best", default_idx, xi, Gs, Gc, od_conj, od_probe, delta_k_z),
        _row("raw no-gap-gate best", raw_idx, xi, Gs, Gc, od_conj, od_probe, delta_k_z),
    ]
    print("\n================ RESULT ================")
    for row in rows:
        _print_row(row)

    out_npz = OUT / "squeezing_tpd_angle_detail.npz"
    np.savez(
        out_npz,
        delta_ghz=D_AXIS_GHZ,
        temp_c=TEMP_AXIS_C,
        angle_deg=ANGLE_AXIS_DEG,
        tpd_mhz=TPD_AXIS_MHZ,
        xi_dB=xi,
        G_s=Gs,
        G_c=Gc,
        gap=gap,
        gap_bad=gap_bad,
        tpd_edge=tpd_edge,
        trusted=trusted,
        delta_k_z=delta_k_z,
        od_conj=od_conj,
        od_probe=od_probe,
        segment_od=segment_od,
        pump_scatter=pump_scatter,
        od_pump=od_pump,
        spatial_overlap_min=spatial_overlap_min,
        phase_max_change=phase_max_change,
        best_idx=np.array(best_idx, dtype=int),
        default_angle_idx=np.array(default_idx, dtype=int),
        raw_idx=np.array(raw_idx, dtype=int),
        summary_labels=np.array([r["label"] for r in rows]),
        summary_values=np.array([
            [r["Delta_GHz"], r["T_C"], r["angle_deg"], r["TPD_MHz"], r["xi_dB"],
             r["G_s"], r["G_c"], r["gap"], r["od_conj"], r["od_probe"], r["delta_k_z"]]
            for r in rows
        ]),
        qe=QE,
        loss=LOSS,
        line_strength=LINE_STRENGTH,
        velocity_step=VELOCITY_STEP,
        velocity_cutoff=VELOCITY_CUTOFF,
        gap_min=GAP_MIN,
        gap_max=GAP_MAX,
    )
    print(f"saved {out_npz}")
    _plot(xi, trusted, rows[:2])


if __name__ == "__main__":
    main()
