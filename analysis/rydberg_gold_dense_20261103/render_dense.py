"""Render complete saved dense results; no interpolation or numerical re-solves.

Run after the dense campaign has written results.json. The renderer refuses
missing, duplicate or non-finite nominal map cells instead of filling them.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm, ListedColormap, LogNorm, Normalize, SymLogNorm
from matplotlib.patches import Patch, Rectangle
import numpy as np

ROOT = Path(__file__).resolve().parent
NOMINAL_OPACITY = .001
FIELD = "sensitivity_nv_cm_sqrt_hz"
GOLD_COLOR = "#fa3459"
BEST_COLOR = "#ff9e1b"


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def edges(values):
    values = np.asarray(values, dtype=float)
    if values.ndim != 1 or values.size < 2 or np.any(np.diff(values) <= 0):
        raise ValueError("Map axes must be strictly increasing with at least two samples")
    mid = (values[:-1] + values[1:]) / 2
    return np.r_[values[0] - (mid[0] - values[0]), mid,
                 values[-1] + values[-1] - mid[-1]]


def finite_number(value):
    return isinstance(value, (float, int)) and np.isfinite(value)


def main():
    source = ROOT / "results.json"
    data = json.loads(source.read_text(encoding="utf-8"))
    config = data["config"]
    temperatures = [float(v) for v in config["temperature_c"]]
    powers = [float(v) for v in config["coupling_power_mw"]]
    modes = list(config["modes"])
    if set(modes) != {"off", "on"} or len(modes) != 2:
        raise ValueError("Expected exactly Doppler off/on modes")
    modes = ["off", "on"]
    if len(temperatures) != 21 or len(powers) != 9:
        raise ValueError("Expected complete 21 by 9 dense map axes")
    if temperatures[0] != 15 or temperatures[-1] != 35 or powers[0] != 20 or powers[-1] != 40:
        raise ValueError("Dense render must retain the original T/Pc domain")
    xedges, yedges = edges(powers), edges(temperatures)
    nominal = {}
    for row in data["optima"]:
        if not np.isclose(float(row["opacity"]), NOMINAL_OPACITY, atol=1e-12, rtol=0):
            continue
        key = (row["doppler"], float(row["temperature_c"]), float(row["coupling_power_mw"]))
        if key in nominal:
            raise ValueError(f"Duplicate nominal map cell: {key}")
        for field in (FIELD, "probe_power_uw", "lo_rabi_mhz", "probe_detuning_mhz"):
            if not finite_number(row.get(field)):
                raise ValueError(f"Invalid {field} at {key}")
        if row[FIELD] <= 0:
            raise ValueError(f"Non-positive sensitivity at {key}")
        nominal[key] = row
    expected = {(mode, t, pc) for mode in modes for t in temperatures for pc in powers}
    if set(nominal) != expected:
        raise ValueError(f"Incomplete nominal map: expected {len(expected)}, found {len(nominal)}")

    # Each branch must be profiled independently across all saved controls.
    # The branch columns attached to the overall winning control row do not
    # provide the competing branch's independent optimum.
    branches = {}
    branch_rows = data.get("branchgaps")
    if not isinstance(branch_rows, list):
        raise ValueError("Expected the separately profiled branchgaps table")
    for row in branch_rows:
        if not np.isclose(float(row["opacity"]), NOMINAL_OPACITY, atol=1e-12, rtol=0):
            continue
        key = (row["doppler"], float(row["temperature_c"]), float(row["coupling_power_mw"]))
        if key in branches:
            raise ValueError(f"Duplicate nominal branch-gap cell: {key}")
        pair = (row.get("central_sensitivity_nv_cm_sqrt_hz"),
                row.get("sideband_sensitivity_nv_cm_sqrt_hz"))
        if not all(finite_number(value) and value > 0 for value in pair):
            raise ValueError(f"Invalid separately profiled branch sensitivities at {key}")
        branches[key] = tuple(float(value) for value in pair)
    if set(branches) != expected:
        raise ValueError(f"Incomplete nominal branch-gap table: expected {len(expected)}, found {len(branches)}")

    def rows(mode):
        return [[nominal[(mode, t, pc)] for pc in powers] for t in temperatures]

    def matrix(mode, field):
        return np.asarray([[row[field] for row in line] for line in rows(mode)], dtype=float)

    best = {mode: min((nominal[(mode, t, pc)] for t in temperatures for pc in powers),
                     key=lambda row: row[FIELD]) for mode in modes}
    gold = {mode: nominal[(mode, 20., 30.)] for mode in modes}
    output_names = []

    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10,
                         "axes.titleweight": "bold", "figure.facecolor": "white"})

    def save(figure, name):
        figure.savefig(ROOT / name, dpi=180, bbox_inches="tight")
        plt.close(figure)
        output_names.append(name)

    def decorate(axis):
        axis.set(xlabel="Coupling power (mW)", ylabel="Temperature (C)")
        axis.set_xticks([20, 25, 30, 35, 40])
        axis.set_yticks([15, 20, 25, 30, 33, 35])
        xx, yy = np.meshgrid(powers, temperatures)
        axis.scatter(xx.ravel(), yy.ravel(), s=3.5, color="white", alpha=.32,
                     linewidths=0, zorder=3)

    def frame(axis, row, color):
        ti = temperatures.index(float(row["temperature_c"]))
        pi = powers.index(float(row["coupling_power_mw"]))
        dx, dy = .06 * (xedges[pi + 1] - xedges[pi]), .08 * (yedges[ti + 1] - yedges[ti])
        axis.add_patch(Rectangle((xedges[pi] + dx, yedges[ti] + dy),
            xedges[pi + 1] - xedges[pi] - 2 * dx,
            yedges[ti + 1] - yedges[ti] - 2 * dy,
            fill=False, edgecolor=color, linewidth=2.1, zorder=5))

    figure, axes = plt.subplots(1, 2, figsize=(12.1, 6.3), layout="constrained")
    for axis, mode in zip(axes, modes):
        values = matrix(mode, FIELD)
        display = axis.pcolormesh(xedges, yedges, values, cmap="viridis_r",
            norm=LogNorm(float(values.min()), float(values.max())),
            edgecolors=(1, 1, 1, .16), linewidth=.25)
        decorate(axis)
        frame(axis, gold[mode], GOLD_COLOR)
        frame(axis, best[mode], BEST_COLOR)
        axis.set_title(f"Doppler {mode} | independent log color scale", fontsize=11)
        figure.colorbar(display, ax=axis, pad=.025, fraction=.047).set_label(
            "PSN sensitivity (nV/cm/sqrt(Hz))")
        b, g = best[mode], gold[mode]
        axis.text(.0, -.17,
            f"Gold (T,Pc) cell after profiling: {g[FIELD]:.4g}\n"
            f"Best sampled cell: {b[FIELD]:.4g} at T={b['temperature_c']:g} C, Pc={b['coupling_power_mw']:g} mW",
            transform=axis.transAxes, fontsize=9, va="top")
    figure.suptitle("Dense Gold-neighborhood sensitivity maps\n"
        "21 x 9 solved coordinates | same domain | assumed opacity participation = 0.001\n"
        "Each cell minimizes sampled probe power, LO Rabi and probe detuning", fontsize=13)
    figure.supxlabel("White dots: solved coordinates. Red frame: Gold (T,Pc). Orange frame: best sampled cell.\n"
        "Discrete cells; no interpolation. Uncalibrated PSN-only model, zero technical noise.",
        fontsize=9)
    save(figure, "dense_sensitivity_maps.png")

    quantities = [("probe_power_uw", "Probe power (uW)", 2., 10.),
                  ("lo_rabi_mhz", "LO Rabi / 2pi (MHz)", 2., 6.),
                  ("probe_detuning_mhz", "Signed probe detuning (MHz)", -6., 6.)]
    figure, axes = plt.subplots(2, 3, figsize=(14.6, 9.0), layout="constrained")
    for mi, mode in enumerate(modes):
        for qi, (field, label, low, high) in enumerate(quantities):
            axis = axes[mi, qi]
            values = matrix(mode, field)
            cmap = "RdBu_r" if field == "probe_detuning_mhz" else "cividis"
            display = axis.pcolormesh(xedges, yedges, values, cmap=cmap,
                norm=Normalize(low, high), edgecolors=(1, 1, 1, .16), linewidth=.25)
            decorate(axis)
            axis.set_title(f"Doppler {mode} | {label}", fontsize=10)
            figure.colorbar(display, ax=axis, pad=.025, fraction=.047).set_label(label)
    figure.suptitle("Controls selected by the dense conditional search\n"
        "Fixed color ranges span the original control bounds", fontsize=13)
    figure.supxlabel("An abrupt change of selected controls can occur between near-degenerate minima.\n"
        "Signed detuning is retained: the finite-IF signal response need not have mirror symmetry.", fontsize=9)
    save(figure, "dense_optimal_controls.png")

    figure, axes = plt.subplots(2, 2, figsize=(12.0, 8.3), layout="constrained")
    colors = ["#335fa6", "#0b9e8b", "#dc8321"]
    for mi, mode in enumerate(modes):
        axis = axes[mi, 0]
        for power, color in zip([20., 30., 40.], colors):
            values = [nominal[(mode, t, power)][FIELD] for t in temperatures]
            axis.semilogy(temperatures, values, "o-", ms=3, lw=1.4,
                          color=color, label=f"Pc = {power:g} mW")
        axis.axvline(33, color=".5", ls=":", lw=1)
        axis.set(xlabel="Temperature (C)", ylabel="PSN sensitivity (nV/cm/sqrt(Hz))",
                 title=f"Doppler {mode} | temperature slices")
        axis.set_xticks([15, 20, 25, 30, 33, 35])
        axis.legend(fontsize=9)
        axis.grid(alpha=.2, which="both")
        axis = axes[mi, 1]
        for temperature, color in zip([20., 30., 35.], colors):
            values = [nominal[(mode, temperature, pc)][FIELD] for pc in powers]
            axis.semilogy(powers, values, "o-", ms=3, lw=1.4,
                          color=color, label=f"T = {temperature:g} C")
        axis.set(xlabel="Coupling power (mW)", ylabel="PSN sensitivity (nV/cm/sqrt(Hz))",
                 title=f"Doppler {mode} | coupling-power slices")
        axis.set_xticks([20, 25, 30, 35, 40])
        axis.legend(fontsize=9)
        axis.grid(alpha=.2, which="both")
    figure.suptitle("Dense solved slices, with hidden controls re-optimized at each point", fontsize=13)
    figure.supxlabel("Markers are actual solved map coordinates; lines only connect them. "
        "Dotted vertical line: paper's 33 C reference.\n"
        "Each panel uses its own logarithmic sensitivity axis; fixed IF = 40 kHz.", fontsize=9)
    save(figure, "dense_slices.png")

    # These masks concern map-coordinate tolerance after profiling, not the
    # separate epsilon-optimal hidden-control bands supplied by the runner.
    figure, axes = plt.subplots(1, 2, figsize=(11.6, 5.5), layout="constrained")
    cmap = ListedColormap(["#eceef3", "#72a5cc", "#f6bd44"])
    norm = BoundaryNorm([-.5, .5, 1.5, 2.5], cmap.N)
    near_best = {}
    for axis, mode in zip(axes, modes):
        values = matrix(mode, FIELD)
        ratio = values / values.min()
        category = np.where(ratio <= 1.01 + 1e-12, 2,
                            np.where(ratio <= 1.03 + 1e-12, 1, 0))
        axis.pcolormesh(xedges, yedges, category, cmap=cmap, norm=norm,
                        edgecolors="white", linewidth=.35)
        decorate(axis)
        frame(axis, best[mode], BEST_COLOR)
        axis.set_title(f"Doppler {mode}", fontsize=11)
        near_best[mode] = {}
        for tolerance in [.01, .03]:
            selected = [nominal[(mode, t, pc)] for ti, t in enumerate(temperatures)
                        for pi, pc in enumerate(powers) if ratio[ti, pi] <= 1 + tolerance + 1e-12]
            near_best[mode][f"{tolerance:.2f}"] = {
                "count": len(selected),
                "temperature_range_c": [min(r["temperature_c"] for r in selected),
                                        max(r["temperature_c"] for r in selected)],
                "coupling_power_range_mw": [min(r["coupling_power_mw"] for r in selected),
                                           max(r["coupling_power_mw"] for r in selected)],
                "coordinates": [[r["temperature_c"], r["coupling_power_mw"]] for r in selected]}
        axis.text(0, -.16, f"Within 1%: {near_best[mode]['0.01']['count']} cells; "
            f"within 3%: {near_best[mode]['0.03']['count']} cells (of {len(temperatures)*len(powers)}).",
            transform=axis.transAxes, fontsize=9, va="top")
    figure.suptitle("Map coordinates close to the best sampled sensitivity\n"
        "Relative tolerance under the same uncalibrated model", fontsize=13)
    figure.legend(handles=[Patch(facecolor=cmap.colors[2], label="Within 1%"),
        Patch(facecolor=cmap.colors[1], label="1-3% above best"),
        Patch(facecolor=cmap.colors[0], label="More than 3% above best")],
        loc="outside lower center", ncol=3, frameon=False, fontsize=9)
    save(figure, "dense_near_best_regions.png")

    signed = {mode: np.array([[100 * (branches[(mode, t, pc)][1] - branches[(mode, t, pc)][0])
                / min(branches[(mode, t, pc)]) for pc in powers] for t in temperatures])
              for mode in modes}
    limit = max(3., max(float(np.abs(v).max()) for v in signed.values()))
    figure, axes = plt.subplots(1, 2, figsize=(11.9, 5.3), layout="constrained")
    for axis, mode in zip(axes, modes):
        display = axis.pcolormesh(xedges, yedges, signed[mode], cmap="RdBu_r",
            norm=SymLogNorm(linthresh=1., vmin=-limit, vmax=limit),
            edgecolors=(1, 1, 1, .16), linewidth=.25)
        decorate(axis)
        axis.set_title(f"Doppler {mode}", fontsize=11)
        figure.colorbar(display, ax=axis, pad=.025, fraction=.047).set_label(
            "Signed branch sensitivity gap (%)")
    figure.suptitle("Separately profiled central and sideband minima\n"
        "100 x (S_sideband - S_central) / min(S_sideband, S_central)", fontsize=13)
    figure.supxlabel("Positive: central branch better. Negative: sideband branch better. "
        "Symmetric-log colors, linear within +/-1%.\n"
        "Each branch independently selects Pp, LO and detuning from all saved candidates. "
        "A small gap means nearly degenerate minima.", fontsize=9)
    save(figure, "dense_branch_gaps.png")
    branch_plot_status = "Rendered from complete branchgaps table; each branch independently profiled over all saved controls"

    diagnostics = {"origin": "AI analysis of complete saved discrete nominal map",
                   "opacity": NOMINAL_OPACITY, "near_best_map_coordinates": near_best,
                   "branch_plot_status": branch_plot_status, "local_slopes": {}}
    for mode in modes:
        log_values = np.log(matrix(mode, FIELD))
        dT, dPc = np.gradient(log_values, temperatures, powers, edge_order=2)
        diagnostics["local_slopes"][mode] = {
            "d_log_sensitivity_per_c_range": [float(dT.min()), float(dT.max())],
            "d_log_sensitivity_per_mw_range": [float(dPc.min()), float(dPc.max())],
            "temperature_adjacent_increase_count": int(np.sum(np.diff(log_values, axis=0) > 1e-10)),
            "coupling_adjacent_increase_count": int(np.sum(np.diff(log_values, axis=1) > 1e-10)),
            "gold_d_log_sensitivity_per_c": float(dT[temperatures.index(20.), powers.index(30.)]),
            "gold_d_log_sensitivity_per_mw": float(dPc[temperatures.index(20.), powers.index(30.)]),
            "scope": "Finite differences of solved profiled map; branch transitions can make derivatives non-smooth"}

    coarse_source = ROOT.parent / "rydberg_gold_map_20261103" / "results.json"
    comparison = []
    if coarse_source.exists():
        coarse = json.loads(coarse_source.read_text(encoding="utf-8"))
        for old in coarse["optima"]:
            if not np.isclose(float(old["opacity"]), NOMINAL_OPACITY, atol=1e-12, rtol=0):
                continue
            key = (old["doppler"], float(old["temperature_c"]), float(old["coupling_power_mw"]))
            current = nominal[key]
            comparison.append({"doppler": key[0], "temperature_c": key[1],
                "coupling_power_mw": key[2], "coarse_sensitivity": old[FIELD],
                "dense_sensitivity": current[FIELD],
                "improvement_percent": 100 * (old[FIELD] - current[FIELD]) / old[FIELD],
                "coarse_controls": {name: old[name] for name in ("probe_power_uw", "lo_rabi_mhz", "probe_detuning_mhz")},
                "dense_controls": {name: current[name] for name in ("probe_power_uw", "lo_rabi_mhz", "probe_detuning_mhz")}})
        diagnostics["coarse_comparison"] = {"rows": comparison,
            "scope": "Same map coordinates and nominal opacity; hidden-control and detuning sampling changed",
            "summary": {mode: {
                "coordinate_count": sum(row["doppler"] == mode for row in comparison),
                "improvement_percent_range": [min(row["improvement_percent"] for row in comparison if row["doppler"] == mode),
                                              max(row["improvement_percent"] for row in comparison if row["doppler"] == mode)]}
                for mode in modes}}
    (ROOT / "render_diagnostics.json").write_text(json.dumps(diagnostics, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    output_names.append("render_diagnostics.json")
    manifest = {"origin": "AI rendering of saved dense campaign results",
        "opacity": NOMINAL_OPACITY, "map_shape": [len(temperatures), len(powers)],
        "scope": "Discrete solved maps and connected solved slices; no interpolated evaluations",
        "input_sha256": {"results.json": sha256(source), "render_dense.py": sha256(Path(__file__))},
        "output_sha256": {name: sha256(ROOT / name) for name in output_names}}
    if coarse_source.exists():
        manifest["input_sha256"]["../rydberg_gold_map_20261103/results.json"] = sha256(coarse_source)
    (ROOT / "render_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"rendered": output_names, "nominal_cells": len(nominal),
                      "branch_plot_status": branch_plot_status}, indent=2))


if __name__ == "__main__":
    main()
