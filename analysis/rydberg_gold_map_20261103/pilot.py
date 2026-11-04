"""Bounded, sampled Gold-neighborhood sensitivity maps using existing GABES APIs.

No production solver is changed. Doppler-on evaluations use explicit finite-IF
probe candidates instead of the UI's static-slope-only candidate shortcut.
"""
from __future__ import annotations

import os

# Each worker solves many small matrices. Avoid nested BLAS oversubscription.
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"

import argparse
import csv
import hashlib
import json
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parents[1]
sys.path.insert(0, str(REPO))

from gabes.schemes.rydberg import MHZ, RydbergEITScheme

TEMPERATURES = [15.0, 20.0, 25.0, 30.0, 33.0, 35.0]
COUPLING_POWERS = [20.0, 30.0, 40.0]
PROBE_POWERS = [2.0, 4.0, 6.0, 8.0, 10.0]
LO_RABIS = [2.0, 3.0, 3.7, 4.5, 6.0]
OPACITIES = [0.0005, 0.001, 0.0015]
MODES = ["off", "on"]
PROBE_WINDOW = 6.0
ENGINE_FILES = [
    "gabes/schemes/rydberg.py", "gabes/rydberg_electrometry.py",
    "gabes/rydberg_opacity.py", "gabes/affine_doppler.py", "gabes/doppler.py",
    "gabes/core.py", "gabes/atoms.py", "gabes/beam.py", "gabes/species.py",
    "gabes/constants.py", "gabes/observables.py",
]


def write_json(path, value):
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def signature():
    files = {name: hashlib.sha256((REPO / name).read_bytes()).hexdigest() for name in ENGINE_FILES}
    files["pilot.py"] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    return hashlib.sha256(json.dumps(files, sort_keys=True).encode()).hexdigest(), files


def task_key(task):
    return json.dumps(task, separators=(",", ":"))


def boundary_flags(row):
    flags = []
    for key, grid in (("probe_power_uw", PROBE_POWERS), ("lo_rabi_mhz", LO_RABIS)):
        if row[key] == grid[0]:
            flags.append(key + "_lower")
        if row[key] == grid[-1]:
            flags.append(key + "_upper")
    if abs(row["probe_detuning_mhz"]) >= PROBE_WINDOW - 0.06:
        flags.append("probe_detuning_boundary")
    return "|".join(flags)


def evaluate(task):
    mode, temperature, coupling, probe, lo, centers = task
    started = time.perf_counter()
    scheme = RydbergEITScheme()
    params = scheme.recommended_defaults(scheme.defaults())["AT electrometry"]
    params.update(doppler=mode, temp_c=temperature, coupling_power_mw=coupling,
                  probe_power_uw=probe, lo_rabi_mhz=lo, mw_detuning_mhz=0.0)
    raw = scheme.compute(params)
    x, nominal_transmission, _ = scheme._transmission(raw["chi_bar"], raw, params)
    if mode == "on":
        targets = list(np.linspace(-PROBE_WINDOW, PROBE_WINDOW, 13))
        for center in centers:
            targets.extend(np.linspace(max(-PROBE_WINDOW, center - 1.0),
                                       min(PROBE_WINDOW, center + 1.0), 25))
        explicit_indices = [int(np.argmin(np.abs(x - value))) for value in targets]
        explicit_indices.extend(scheme._static_slope_indices(x, nominal_transmission).tolist())
        indices = np.unique(explicit_indices)
        indices = indices[np.abs(x[indices]) <= PROBE_WINDOW + 0.03]
        atomic = scheme._superheterodyne_atomic_phasors(
            raw, params, x, nominal_transmission, static_indices=indices)
    else:
        atomic = scheme._superheterodyne_atomic_phasors(raw, params, x, nominal_transmission)
        mask = np.abs(x[atomic[0]]) <= PROBE_WINDOW + 0.03
        atomic = (atomic[0][mask], atomic[1][mask], "bounded full probe scan")

    rows = []
    for scale in OPACITIES:
        scenario = dict(params, effective_opacity_participation=scale)
        _, transmission, _ = scheme._transmission(raw["chi_bar"], raw, scenario)
        result = scheme._superheterodyne_readout(raw, scenario, x, transmission, atomic_response=atomic)
        if result is None:
            raise RuntimeError("No finite-IF response for " + task_key(task))
        sensitivity = result["sensitivity"]
        row = dict(
            doppler=mode, temperature_c=temperature, coupling_power_mw=coupling,
            probe_power_uw=probe, lo_rabi_mhz=lo, opacity=scale,
            probe_detuning_mhz=float(result["optimum_detuning_mhz"]),
            sensitivity_nv_cm_sqrt_hz=float(sensitivity.total_field_asd_nv_cm_per_sqrt_hz),
            psn_nv_cm_sqrt_hz=float(sensitivity.psn_field_asd_nv_cm_per_sqrt_hz),
            transmission=float(result["optimum_transmission"]),
            current_response_a_per_v_m=float(result["current_responsivity_a_per_v_m"]),
            detected_power_w=float(result["signal_detector_power_w"]),
            candidate_count=int(atomic[0].size), refined=bool(centers),
            actual_probe_search="explicit finite-IF candidates" if mode == "on" else "all bounded scan points",
        )
        if not np.isfinite(row["sensitivity_nv_cm_sqrt_hz"]) or row["sensitivity_nv_cm_sqrt_hz"] <= 0:
            raise RuntimeError("Invalid field sensitivity")
        row["boundary_flags"] = boundary_flags(row)
        rows.append(row)
    return dict(task=task, key=task_key(task), rows=rows, elapsed_s=time.perf_counter() - started)


def run_tasks(tasks, workers, checkpoint, records):
    todo = [task for task in tasks if task_key(task) not in records]
    if not todo:
        return
    print(f"Evaluating {len(todo)} control states with {workers} workers", flush=True)
    start = time.perf_counter()
    with checkpoint.open("a", encoding="utf-8") as stream, ProcessPoolExecutor(max_workers=workers) as pool:
        futures = {pool.submit(evaluate, task): task for task in todo}
        for completed, future in enumerate(as_completed(futures), 1):
            result = future.result()
            records[result["key"]] = result
            stream.write(json.dumps(result, ensure_ascii=False, allow_nan=False) + "\n")
            stream.flush()
            if completed == 1 or completed % 25 == 0 or completed == len(todo):
                print(f"{completed}/{len(todo)} states, {time.perf_counter()-start:.1f}s", flush=True)


def export_csv(path, rows):
    with path.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def select_optima(rows):
    selected = []
    for mode in MODES:
        for scale in OPACITIES:
            for temperature in TEMPERATURES:
                for coupling in COUPLING_POWERS:
                    candidates = [r for r in rows if (r["doppler"], r["opacity"], r["temperature_c"],
                                  r["coupling_power_mw"]) == (mode, scale, temperature, coupling)]
                    selected.append(dict(min(candidates, key=lambda r: r["sensitivity_nv_cm_sqrt_hz"])))
    return selected


def make_figures(optima):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm

    def matrix(mode, quantity):
        return np.array([[next(r[quantity] for r in optima if
            (r["doppler"], r["opacity"], r["temperature_c"], r["coupling_power_mw"]) ==
            (mode, 0.001, t, p)) for p in COUPLING_POWERS] for t in TEMPERATURES])

    figure, axes = plt.subplots(1, 2, figsize=(11.6, 5.0), constrained_layout=True)
    for axis, mode, title in zip(axes, MODES, ("Compact model (Doppler off)", "Thermal model (Doppler on)")):
        values = matrix(mode, "sensitivity_nv_cm_sqrt_hz")
        display = axis.imshow(values, origin="lower", aspect="auto", norm=LogNorm(), cmap="viridis_r")
        axis.set_xticks(range(len(COUPLING_POWERS)), [f"{p:g}" for p in COUPLING_POWERS])
        axis.set_yticks(range(len(TEMPERATURES)), [f"{t:g}" for t in TEMPERATURES])
        axis.set(xlabel="Coupling power (mW)", ylabel="Temperature (C)", title=title)
        for i in range(values.shape[0]):
            for j in range(values.shape[1]):
                axis.text(j, i, f"{values[i,j]:.2g}", ha="center", va="center", color="white",
                          bbox=dict(facecolor="black", alpha=0.35, edgecolor="none", pad=2))
        figure.colorbar(display, ax=axis, label="Best sampled sensitivity (nV/cm/sqrt(Hz))")
        axis.plot(1, 1, marker="s", fillstyle="none", color="red", markersize=24, markeredgewidth=1.6)
    figure.suptitle("Gold-neighborhood pilot | 50 mm, 0.15 mm, IF 40 kHz\n"
                    "Optimize sampled probe power, LO Rabi and probe detuning | opacity=0.001 (uncalibrated)")
    figure.savefig(ROOT / "optimized_maps.png", dpi=170)
    plt.close(figure)

    figure, axes = plt.subplots(2, 2, figsize=(10.4, 8.2), constrained_layout=True)
    for i, mode in enumerate(MODES):
        for j, (quantity, label) in enumerate((("probe_power_uw", "Optimal sampled probe power (uW)"),
                                              ("lo_rabi_mhz", "Optimal sampled LO Rabi / 2pi (MHz)"))):
            axis = axes[i, j]
            values = matrix(mode, quantity)
            display = axis.imshow(values, origin="lower", aspect="auto", cmap="cividis")
            axis.set_xticks(range(len(COUPLING_POWERS)), [f"{p:g}" for p in COUPLING_POWERS])
            axis.set_yticks(range(len(TEMPERATURES)), [f"{t:g}" for t in TEMPERATURES])
            axis.set(xlabel="Coupling power (mW)", ylabel="Temperature (C)", title=f"Doppler {mode} | {label}")
            for ii in range(values.shape[0]):
                for jj in range(values.shape[1]):
                    axis.text(jj, ii, f"{values[ii,jj]:g}", ha="center", va="center", color="white",
                              bbox=dict(facecolor="black", alpha=0.35, edgecolor="none", pad=2))
            figure.colorbar(display, ax=axis)
    figure.suptitle("Sampled optimizer settings | edge values mark active search bounds")
    figure.savefig(ROOT / "optimal_controls.png", dpi=170)
    plt.close(figure)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workers", type=int, default=2)
    args = parser.parse_args()
    ROOT.mkdir(parents=True, exist_ok=True)
    run_signature, hashes = signature()
    config = dict(as_of="2026-11-03", signature=run_signature, temperature_c=TEMPERATURES,
        coupling_power_mw=COUPLING_POWERS, probe_power_uw=PROBE_POWERS, lo_rabi_mhz=LO_RABIS,
        opacities=OPACITIES, modes=MODES, probe_window_mhz=PROBE_WINDOW,
        fixed=dict(cell_mm=50, beam_diameter_mm=0.15, if_khz=40, rf_detuning_mhz=0,
                   coupling_detuning_mhz=0, temperature_model="Linked"),
        scope="Sampled local optima, compact homogeneous propagation; not a global optimum or validated forecast",
        provenance="AI-selected bounded pilot around user-adopted Ju2026 Gold; no fit to the sensitivity reference")
    config_path = ROOT / "config.json"
    if config_path.exists() and json.loads(config_path.read_text(encoding="utf-8")) != config:
        raise RuntimeError("Existing pilot inputs or solver hashes differ; use a separate analysis directory")
    write_json(config_path, config)
    checkpoint = ROOT / "checkpoint.jsonl"
    records = {}
    if checkpoint.exists():
        for line in checkpoint.read_text(encoding="utf-8").splitlines():
            result = json.loads(line)
            records[result["key"]] = result
    coarse = [(mode, t, pc, pp, lo, []) for mode in MODES for t in TEMPERATURES
              for pc in COUPLING_POWERS for pp in PROBE_POWERS for lo in LO_RABIS]
    run_tasks(coarse, args.workers, checkpoint, records)
    coarse_rows = [row for key in map(task_key, coarse) for row in records[key]["rows"]]

    # Refine the two best control tuples per cell/scenario and the Gold baseline.
    selected = {}
    for mode in MODES:
        if mode == "off":
            continue  # already searches every bounded point of the 801-point scan
        for scale in OPACITIES:
            for t in TEMPERATURES:
                for pc in COUPLING_POWERS:
                    group = [r for r in coarse_rows if (r["doppler"], r["opacity"], r["temperature_c"],
                             r["coupling_power_mw"]) == (mode, scale, t, pc)]
                    group.sort(key=lambda r: r["sensitivity_nv_cm_sqrt_hz"])
                    for row in group[:2]:
                        key = (mode, t, pc, row["probe_power_uw"], row["lo_rabi_mhz"])
                        selected.setdefault(key, set()).add(row["probe_detuning_mhz"])
    selected.setdefault(("on", 20.0, 30.0, 6.0, 3.7), set()).add(0.0)
    refined = [(*key, sorted(centers)) for key, centers in selected.items()]
    run_tasks(refined, args.workers, checkpoint, records)
    all_rows = coarse_rows + [row for key in map(task_key, refined) for row in records[key]["rows"]]
    optima = select_optima(all_rows)
    baselines = []
    for mode in MODES:
        for scale in OPACITIES:
            candidates = [r for r in all_rows if (r["doppler"], r["opacity"], r["temperature_c"],
                r["coupling_power_mw"], r["probe_power_uw"], r["lo_rabi_mhz"]) ==
                (mode, scale, 20.0, 30.0, 6.0, 3.7)]
            baselines.append(dict(min(candidates, key=lambda r: r["sensitivity_nv_cm_sqrt_hz"])))
    best = [min((r for r in optima if r["doppler"] == mode and r["opacity"] == 0.001),
                key=lambda r: r["sensitivity_nv_cm_sqrt_hz"]) for mode in MODES]
    summary = dict(config=config, baseline=baselines, best_nominal=best, optima=optima,
        coarse_control_states=len(coarse), refined_control_states=len(refined),
        elapsed_worker_s=sum(records[task_key(task)]["elapsed_s"] for task in coarse + refined),
        numerical_status="sampled search complete; convergence checks recorded separately",
        experimental_status="UNCALIBRATED; model-mode discrepancy must be assessed")
    write_json(ROOT / "results.json", summary)
    export_csv(ROOT / "control_points.csv", all_rows)
    export_csv(ROOT / "optima.csv", optima)
    write_json(ROOT / "source_manifest.json", dict(signature=run_signature, code_sha256=hashes,
        python=sys.version, numpy=np.__version__, blas_threads=1))
    make_figures(optima)
    print(json.dumps(dict(baseline=baselines, best_nominal=best), ensure_ascii=False, indent=2), flush=True)


if __name__ == "__main__":
    main()
