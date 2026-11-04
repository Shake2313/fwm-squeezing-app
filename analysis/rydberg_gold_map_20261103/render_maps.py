"""Render saved discrete pilot results; never re-run or alter numerical inputs."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, Normalize
from matplotlib.patches import Patch, Rectangle
import numpy as np

ROOT = Path(__file__).resolve().parent


def edges(values):
    values = np.asarray(values, dtype=float)
    mid = (values[:-1] + values[1:]) / 2
    return np.r_[values[0] - (mid[0] - values[0]), mid,
                 values[-1] + (values[-1] - mid[-1])]


def main():
    result = json.loads((ROOT / "results.json").read_text(encoding="utf-8"))
    config = result["config"]
    temperatures = config["temperature_c"]
    powers = config["coupling_power_mw"]
    xedges, yedges = edges(powers), edges(temperatures)

    def matrix(mode, quantity):
        return np.array([[next(row[quantity] for row in result["optima"] if
            (row["doppler"], row["opacity"], row["temperature_c"], row["coupling_power_mw"])
            == (mode, .001, temperature, power)) for power in powers]
            for temperature in temperatures])

    def frame(axis, ti, pi, color):
        dx = (xedges[pi + 1] - xedges[pi]) * .025
        dy = (yedges[ti + 1] - yedges[ti]) * .055
        axis.add_patch(Rectangle((xedges[pi] + dx, yedges[ti] + dy),
            xedges[pi + 1] - xedges[pi] - 2 * dx,
            yedges[ti + 1] - yedges[ti] - 2 * dy,
            fill=False, edgecolor=color, linewidth=2.3, zorder=5))

    def draw(axis, values, norm, cmap, format_value):
        display = axis.pcolormesh(xedges, yedges, values, norm=norm, cmap=cmap,
                                 edgecolors="#eeeeee", linewidth=.55)
        axis.set(xlabel="Coupling power (mW)", ylabel="Temperature (C)")
        axis.set_xticks(powers)
        axis.set_yticks(temperatures)
        for ti, temperature in enumerate(temperatures):
            for pi, power in enumerate(powers):
                axis.text(power, temperature, format_value(values[ti, pi]),
                    ha="center", va="center", color="white", fontsize=10,
                    bbox=dict(facecolor="black", alpha=.48, edgecolor="none", pad=2))
        return display

    plt.rcParams.update({"font.size": 10, "axes.titleweight": "bold"})
    figure, axes = plt.subplots(1, 2, figsize=(12.6, 6.4))
    figure.subplots_adjust(left=.065, right=.96, bottom=.21, top=.79, wspace=.39)
    for axis, mode in zip(axes, config["modes"]):
        values = matrix(mode, "sensitivity_nv_cm_sqrt_hz")
        display = draw(axis, values, LogNorm(vmin=values.min(), vmax=values.max()),
                       "viridis_r", lambda value: f"{value:.3g}")
        axis.set_title(f"Doppler {mode} | independent log color scale", fontsize=11)
        figure.colorbar(display, ax=axis, pad=.025, fraction=.047).set_label(
            "PSN sensitivity (nV/cm/sqrt(Hz))")
        frame(axis, temperatures.index(20), powers.index(30), "#ff385c")
        ti, pi = np.unravel_index(np.argmin(values), values.shape)
        frame(axis, ti, pi, "#ffab23")
    figure.suptitle("Gold-neighborhood sampled sensitivity maps\n"
        "50 mm cell | 0.15 mm diameter | IF 40 kHz | assumed opacity = 0.001",
        fontsize=14, y=.96)
    figure.legend(handles=[Patch(facecolor="none", edgecolor="#ff385c", linewidth=2,
        label="Gold (T, Pc) cell, with controls optimized"),
        Patch(facecolor="none", edgecolor="#ffab23", linewidth=2,
        label="Best sampled cell, at map upper bounds")], loc="lower center",
        bbox_to_anchor=(.5, .095), ncol=2, frameon=False, fontsize=10)
    figure.text(.5, .055, "Each tile: minimum over sampled probe power, LO Rabi and probe detuning. "
        "Discrete samples; no interpolation.", ha="center", fontsize=9)
    figure.text(.5, .025, "Uncalibrated compact model; zero technical noise. "
        "Absolute values are conditional model outputs, not measured performance.",
        ha="center", fontsize=9, color="#9b2034")
    figure.savefig(ROOT / "optimized_maps.png", dpi=170, bbox_inches="tight")
    plt.close(figure)

    figure, axes = plt.subplots(2, 2, figsize=(11.5, 9.0))
    figure.subplots_adjust(left=.075, right=.96, bottom=.12, top=.86,
                           hspace=.35, wspace=.35)
    quantities = [("probe_power_uw", "Probe power (uW)", config["probe_power_uw"]),
                  ("lo_rabi_mhz", "LO Rabi / 2pi (MHz)", config["lo_rabi_mhz"])]
    for mi, mode in enumerate(config["modes"]):
        for qi, (quantity, label, grid) in enumerate(quantities):
            axis = axes[mi, qi]
            display = draw(axis, matrix(mode, quantity), Normalize(grid[0], grid[-1]),
                           "cividis", lambda value: f"{value:g}")
            axis.set_title(f"Doppler {mode} | {label}", fontsize=11)
            figure.colorbar(display, ax=axis, pad=.025, fraction=.047,
                            ticks=grid).set_label(label)
    figure.suptitle("Controls selected in nominal-opacity sensitivity maps\n"
                   "Fixed color ranges show full sampled control bounds", fontsize=14, y=.96)
    figure.text(.5, .055, "Pp = 10 uW reaches upper search bound; LO = 2 MHz reaches lower bound. "
        "All values are discrete samples.", ha="center", fontsize=9)
    figure.text(.5, .025, "Fixed coupling and RF detuning = 0; probe detuning searched within +/-6 MHz. "
        "See optima.csv for selected detunings.", ha="center", fontsize=9)
    figure.savefig(ROOT / "optimal_controls.png", dpi=170, bbox_inches="tight")
    plt.close(figure)

    names = ["render_maps.py", "results.json", "optimized_maps.png", "optimal_controls.png"]
    manifest = {"origin": "AI rendering of saved discrete results", "opacity": .001,
                "sha256": {name: hashlib.sha256((ROOT / name).read_bytes()).hexdigest()
                           for name in names}}
    (ROOT / "render_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
