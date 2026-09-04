"""Shared runner for detailed two-photon-detuning and crossing-angle scans."""
from __future__ import annotations

import os
import sys
import threading
import time
import warnings
from concurrent.futures import ThreadPoolExecutor, as_completed
from contextlib import contextmanager
from dataclasses import dataclass
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from gabes.core import blas_single_thread  # noqa: E402
from gabes.schemes import fwm  # noqa: E402


P_PUMP = 0.6
P_SEED = 8e-6
LINE_STRENGTH = 0.74
BRANCH = -1
QE = 1.0
LOSS = 0.0
VELOCITY_STEP = 4.0
VELOCITY_CUTOFF = 3.0
GAP_MIN = 0.5
GAP_MAX = 1.5


@dataclass(frozen=True)
class ScanConfig:
    version: str
    delta_axis_ghz: np.ndarray
    temperature_axis_c: np.ndarray
    angle_axis_deg: np.ndarray
    two_photon_axis_mhz: np.ndarray
    reference_delta_ghz: float
    reference_temperature_c: float

    @property
    def output_dir(self) -> Path:
        return (ROOT / "analysis" / "squeezing"
                / f"squeezing_frontier_ideal_{self.version}_tpd_angle_detail")

    @property
    def document_figure(self) -> Path:
        return (ROOT / "docs" / "squeezing_report"
                / f"squeezing_frontier_ideal_{self.version}_tpd_angle_detail.png")


_TABLE_STATE = threading.local()


@contextmanager
def _reuse_atomic_tables():
    """Reuse angle-independent susceptibility tables within each worker."""
    original = fwm.chi_matrix_table

    def cached(*args, **kwargs):
        cache = getattr(_TABLE_STATE, "cache", None)
        if cache is None:
            return original(*args, **kwargs)
        order = int(kwargs.get("n_f", 1))
        if order not in cache:
            cache[order] = original(*args, **kwargs)
        return cache[order]

    fwm.chi_matrix_table = cached
    try:
        yield
    finally:
        fwm.chi_matrix_table = original


def _compute_one(config: ScanConfig, i_delta: int, i_temp: int, i_angle: int) -> dict:
    delta_ghz = float(config.delta_axis_ghz[i_delta])
    temperature_c = float(config.temperature_axis_c[i_temp])
    angle = float(config.angle_axis_deg[i_angle])
    center = fwm.branch_center_GHz(delta_ghz, BRANCH)
    scan_min = center + float(config.two_photon_axis_mhz[0]) * 1e-3
    scan_max = center + float(config.two_photon_axis_mhz[-1]) * 1e-3

    spec = fwm.compute_spectrum(
        delta_ghz,
        T=temperature_c + 273.15,
        P_pump=P_PUMP,
        P_probe=P_SEED,
        line_strength=LINE_STRENGTH,
        loss_frac=LOSS,
        qe=QE,
        coarse_points=config.two_photon_axis_mhz.size,
        fine_points=0,
        scan_min=scan_min,
        scan_max=scan_max,
        velocity_step=VELOCITY_STEP,
        velocity_cutoff=VELOCITY_CUTOFF,
        branch=BRANCH,
        phase_detail=fwm.PHASE_ULTRA,
        pump_probe_angle_deg=angle,
    )

    two_photon = (spec["probe_axis_GHz"] - center) * 1e3
    if not np.allclose(two_photon, config.two_photon_axis_mhz,
                       atol=1e-9, rtol=0.0):
        raise RuntimeError("two-photon detuning axis mismatch")

    hardened = spec["hardened_noise"] or {}
    n_points = config.two_photon_axis_mhz.size
    return {
        "i_delta": i_delta,
        "i_temp": i_temp,
        "i_angle": i_angle,
        "xi": np.asarray(spec["S_dB"], dtype=float),
        "Gs": np.asarray(spec["G_s"], dtype=float),
        "Gc": np.asarray(spec["G_c"], dtype=float),
        "delta_k_z": np.asarray(spec["delta_k_z"], dtype=float),
        "od_conj": np.asarray(
            hardened.get("od_conj_arr", np.full(n_points, np.nan)), dtype=float),
        "od_probe": np.asarray(
            hardened.get("od_probe_lin_arr", np.full(n_points, np.nan)), dtype=float),
        "segment_od": float(spec.get("segment_absorption_od", np.nan)),
        "pump_scatter": float(hardened.get("pump_scatter_noise", np.nan)),
        "od_pump": float(hardened.get("od_pump", np.nan)),
        "spatial_overlap_min": float(spec.get("ultra_spatial_overlap_min", np.nan)),
        "phase_max_change": float(spec.get("ultra_phase_max_change", np.nan)),
    }


def _compute_delta_temperature(config: ScanConfig, i_delta: int, i_temp: int):
    _TABLE_STATE.cache = {}
    try:
        return [
            _compute_one(config, i_delta, i_temp, i_angle)
            for i_angle in range(config.angle_axis_deg.size)
        ]
    finally:
        del _TABLE_STATE.cache


def _best(masked_xi: np.ndarray) -> tuple[int, ...]:
    flat = int(np.nanargmin(masked_xi))
    return np.unravel_index(flat, masked_xi.shape)


def _nanmin(array: np.ndarray, axis=None) -> np.ndarray:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        return np.nanmin(array, axis=axis)


def _row(config, label, idx, xi, Gs, Gc, od_conj, od_probe, delta_k_z):
    i_delta, i_temp, i_angle, i_two_photon = idx
    return {
        "label": label,
        "Delta_GHz": float(config.delta_axis_ghz[i_delta]),
        "T_C": float(config.temperature_axis_c[i_temp]),
        "angle_deg": float(config.angle_axis_deg[i_angle]),
        "TPD_MHz": float(config.two_photon_axis_mhz[i_two_photon]),
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
        f"{row['label']}: diagnostic={row['xi_dB']:.4f} dB, "
        f"Δ={row['Delta_GHz']:+.3f} GHz, T={row['T_C']:.1f} C, "
        f"δ={row['TPD_MHz']:+.1f} MHz, angle={row['angle_deg']:.3f} deg, "
        f"Gs={row['G_s']:.2f}, gap={row['gap']:.3f}, "
        f"od_c={row['od_conj']:.4f}, od_p={row['od_probe']:.4f}, "
        f"dkz={row['delta_k_z']:.2e} 1/m")


def _plot(config: ScanConfig, xi, trusted, summary_rows) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    trusted_xi = np.where(trusted, xi, np.nan)
    angle_detuning = _nanmin(trusted_xi, axis=(0, 1))
    delta_temperature = _nanmin(trusted_xi, axis=(2, 3))
    best_vs_angle = _nanmin(angle_detuning, axis=1)
    best_detuning_idx = np.nanargmin(angle_detuning, axis=1)
    best, default = summary_rows

    fig, axes = plt.subplots(1, 3, figsize=(17.5, 4.8))
    cmap_angle = plt.get_cmap("magma_r").copy()
    cmap_angle.set_bad("#eeeeee")
    cmap_dt = plt.get_cmap("viridis_r").copy()
    cmap_dt.set_bad("#eeeeee")

    im0 = axes[0].imshow(
        angle_detuning.T, origin="lower", aspect="auto",
        extent=[config.angle_axis_deg[0], config.angle_axis_deg[-1],
                config.two_photon_axis_mhz[0], config.two_photon_axis_mhz[-1]],
        cmap=cmap_angle, vmin=np.nanmin(angle_detuning),
        vmax=min(-13.0, np.nanmax(angle_detuning)))
    fig.colorbar(im0, ax=axes[0], label="best gain-derived diagnostic [dB]")
    axes[0].scatter([best["angle_deg"]], [best["TPD_MHz"]],
                    s=70, c="cyan", edgecolor="k", zorder=5)
    axes[0].scatter([default["angle_deg"]], [default["TPD_MHz"]],
                    s=55, c="white", edgecolor="k", zorder=5)
    axes[0].axvline(fwm.SEEDED_PHASE_ANGLE_DEG, color="white", ls=":", lw=1.0)
    axes[0].set_xlabel("pump-probe angle [deg]")
    axes[0].set_ylabel("two-photon detuning δ [MHz]")
    axes[0].set_title("minimum over one-photon detuning and T")

    axes[1].plot(config.angle_axis_deg, best_vs_angle, "o-",
                 color="#1f77b4", ms=3, lw=1.6)
    axes[1].axvline(fwm.SEEDED_PHASE_ANGLE_DEG, color="0.25", ls=":", lw=1.0,
                    label="default 0.32 deg")
    axes[1].scatter([best["angle_deg"]], [best["xi_dB"]], c="cyan",
                    edgecolor="k", zorder=5, label="scan best")
    axes[1].scatter([default["angle_deg"]], [default["xi_dB"]], c="white",
                    edgecolor="k", zorder=5, label="default-angle best")
    for angle, detuning in zip(
            config.angle_axis_deg, config.two_photon_axis_mhz[best_detuning_idx]):
        if int(round(angle * 100)) % 4 == 0:
            index = int(np.argmin(np.abs(config.angle_axis_deg - angle)))
            axes[1].annotate(f"{detuning:.0f}", (angle, best_vs_angle[index]),
                             textcoords="offset points", xytext=(0, 6),
                             ha="center", fontsize=6)
    axes[1].set_xlabel("pump-probe angle [deg]")
    axes[1].set_ylabel("best gain-derived diagnostic [dB]")
    axes[1].set_title("angle sensitivity; labels show δ [MHz]")
    axes[1].grid(alpha=0.3)
    axes[1].legend(fontsize=7)

    im2 = axes[2].imshow(
        delta_temperature, origin="lower", aspect="auto",
        extent=[config.temperature_axis_c[0], config.temperature_axis_c[-1],
                config.delta_axis_ghz[0], config.delta_axis_ghz[-1]],
        cmap=cmap_dt)
    fig.colorbar(im2, ax=axes[2], label="best gain-derived diagnostic [dB]")
    axes[2].scatter([best["T_C"]], [best["Delta_GHz"]],
                    c="cyan", edgecolor="k", s=70, zorder=5)
    axes[2].scatter([config.reference_temperature_c], [config.reference_delta_ghz],
                    c="white", edgecolor="k", s=45, zorder=5,
                    label=f"coarse {config.version} diagnostic")
    axes[2].set_xlabel("T [C]")
    axes[2].set_ylabel("one-photon detuning Δ [GHz]")
    axes[2].set_title("minimum over δ and angle")
    axes[2].legend(fontsize=7)

    fig.tight_layout()
    config.document_figure.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(config.document_figure, dpi=150)
    print(f"wrote {config.document_figure}")


def run(config: ScanConfig) -> None:
    config.output_dir.mkdir(parents=True, exist_ok=True)
    shape4 = (config.delta_axis_ghz.size, config.temperature_axis_c.size,
              config.angle_axis_deg.size, config.two_photon_axis_mhz.size)
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

    jobs = [(i_delta, i_temp)
            for i_delta in range(config.delta_axis_ghz.size)
            for i_temp in range(config.temperature_axis_c.size)]
    workers = min(os.cpu_count() or 1, 8)
    print(f"Detailed {config.version} gain-derived diagnostic scan")
    print(f"  grid: Δ {shape4[0]} x T {shape4[1]} x angle {shape4[2]} "
          f"x δ {shape4[3]} = {np.prod(shape4)} sampled points")
    print(f"  atomic solves: {len(jobs)} groups on {workers} workers; "
          f"{shape4[2]} angle propagations per group")

    started = time.time()
    with (_reuse_atomic_tables(), blas_single_thread(),
          ThreadPoolExecutor(max_workers=workers) as executor):
        futures = {executor.submit(_compute_delta_temperature, config, *job): job
                   for job in jobs}
        completed = 0
        for future in as_completed(futures):
            group_results = future.result()
            for result in group_results:
                i_delta = result["i_delta"]
                i_temp = result["i_temp"]
                i_angle = result["i_angle"]
                index = (i_delta, i_temp, i_angle)
                xi[index] = result["xi"]
                Gs[index] = result["Gs"]
                Gc[index] = result["Gc"]
                delta_k_z[index] = result["delta_k_z"]
                od_conj[index] = result["od_conj"]
                od_probe[index] = result["od_probe"]
                segment_od[index] = result["segment_od"]
                pump_scatter[index] = result["pump_scatter"]
                od_pump[index] = result["od_pump"]
                spatial_overlap_min[index] = result["spatial_overlap_min"]
                phase_max_change[index] = result["phase_max_change"]
            completed += len(group_results)
            if completed % 50 == 0 or completed == int(np.prod(shape3)):
                print(f"  {completed}/{np.prod(shape3)} angle propagations "
                      f"({time.time() - started:.1f}s)", flush=True)

    gap = Gs - Gc
    detuning_edge = np.zeros(shape4, dtype=bool)
    detuning_edge[..., 0] = True
    detuning_edge[..., -1] = True
    gap_bad = (gap < GAP_MIN) | (gap > GAP_MAX)
    trusted = np.isfinite(xi) & ~gap_bad & ~detuning_edge
    trusted_xi = np.where(trusted, xi, np.nan)
    best_idx = _best(trusted_xi)
    default_i_angle = int(np.argmin(
        np.abs(config.angle_axis_deg - fwm.SEEDED_PHASE_ANGLE_DEG)))
    default_xi = np.full_like(trusted_xi, np.nan)
    default_xi[:, :, default_i_angle, :] = trusted_xi[:, :, default_i_angle, :]
    default_idx = _best(default_xi)
    raw_idx = _best(np.where(np.isfinite(xi) & ~detuning_edge, xi, np.nan))
    rows = [
        _row(config, "trusted 4D best", best_idx, xi, Gs, Gc,
             od_conj, od_probe, delta_k_z),
        _row(config, "trusted default-angle best", default_idx, xi, Gs, Gc,
             od_conj, od_probe, delta_k_z),
        _row(config, "raw no-gap-gate best", raw_idx, xi, Gs, Gc,
             od_conj, od_probe, delta_k_z),
    ]
    for row in rows:
        _print_row(row)

    out_npz = config.output_dir / "squeezing_tpd_angle_detail.npz"
    np.savez(
        out_npz,
        delta_ghz=config.delta_axis_ghz,
        temp_c=config.temperature_axis_c,
        angle_deg=config.angle_axis_deg,
        tpd_mhz=config.two_photon_axis_mhz,
        xi_dB=xi,
        G_s=Gs,
        G_c=Gc,
        gap=gap,
        gap_bad=gap_bad,
        tpd_edge=detuning_edge,
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
        summary_labels=np.array([row["label"] for row in rows]),
        summary_values=np.array([
            [row["Delta_GHz"], row["T_C"], row["angle_deg"], row["TPD_MHz"],
             row["xi_dB"], row["G_s"], row["G_c"], row["gap"], row["od_conj"],
             row["od_probe"], row["delta_k_z"]]
            for row in rows]),
        qe=QE,
        loss=LOSS,
        line_strength=LINE_STRENGTH,
        velocity_step=VELOCITY_STEP,
        velocity_cutoff=VELOCITY_CUTOFF,
        gap_min=GAP_MIN,
        gap_max=GAP_MAX,
    )
    print(f"saved {out_npz}")
    _plot(config, xi, trusted, rows[:2])
