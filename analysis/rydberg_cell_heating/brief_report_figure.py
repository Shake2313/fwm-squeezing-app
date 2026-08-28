"""Build the compact comparison figure used by the Korean brief report."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any, Mapping

import numpy as np


def _metric(point: Mapping[str, Any], name: str) -> float:
    return float(point["metrics"][name]["value"])


def _finite_if_sensitivity(point: Mapping[str, Any]) -> float:
    return float(point["electrometry"]["field_sensitivity_nv_cm_sqrt_hz"])


def build_figure(results_path: Path, output_stem: Path) -> dict[str, float]:
    """Create a two-panel summary and return the plotted comparison values."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    results = json.loads(results_path.read_text(encoding="utf-8"))
    sweeps = results["model"]["sweeps"]
    eit = sweeps["EIT"]
    at = sweeps["AT electrometry"]

    temperature = np.asarray([point["temperature_c"] for point in eit], dtype=float)
    slope = np.asarray(
        [_metric(point, "max_spectral_slope_per_mhz") for point in eit],
        dtype=float,
    )
    contrast = np.asarray(
        [_metric(point, "eit_peak_contrast") for point in eit],
        dtype=float,
    )
    model_sensitivity = np.asarray(
        [_finite_if_sensitivity(point) for point in at],
        dtype=float,
    )

    sensitivity_input = next(
        dataset
        for dataset in results["inputs"]
        if dataset["kind"] == "sensitivity" and dataset["state"] == "loaded"
    )
    measured_temperature = np.asarray(
        [row["temperature_c"] for row in sensitivity_input["rows"]], dtype=float
    )
    measured_sensitivity = np.asarray(
        [
            row["sensitivity_nv_cm_sqrt_hz"]
            for row in sensitivity_input["rows"]
        ],
        dtype=float,
    )
    measured_uncertainty = np.asarray(
        [
            row.get("uncertainty_nv_cm_sqrt_hz") or 0.0
            for row in sensitivity_input["rows"]
        ],
        dtype=float,
    )

    room_index = int(np.argmin(np.abs(temperature - 20.0)))
    measured_room_index = int(np.argmin(np.abs(measured_temperature - 20.0)))
    room_scale = (
        measured_sensitivity[measured_room_index]
        / model_sensitivity[room_index]
    )
    normalized_model = model_sensitivity * room_scale

    colors = {
        "teal": "#0f766e",
        "blue": "#2563eb",
        "orange": "#ea580c",
        "red": "#b91c1c",
        "gray": "#64748b",
    }
    fig, axes = plt.subplots(1, 2, figsize=(11.4, 4.35))

    ax = axes[0]
    ax.plot(
        temperature,
        slope / np.max(slope),
        "o-",
        color=colors["teal"],
        lw=1.8,
        ms=4,
        label="spectral slope / max",
    )
    ax.plot(
        temperature,
        contrast / np.max(contrast),
        "s-",
        color=colors["blue"],
        lw=1.6,
        ms=4,
        label="EIT contrast / max",
    )
    model_optimum = float(
        results["model"]["optima"]["max_spectral_slope"]["temperature_c"]
    )
    ax.axvline(
        model_optimum,
        color=colors["teal"],
        ls="--",
        lw=1.2,
        label=f"GABES slope optimum: {model_optimum:g} C",
    )
    ax.axvline(
        45.0,
        color=colors["orange"],
        ls=":",
        lw=2.0,
        label="PPT clearest EIT: 45 C",
    )
    ax.set_xlabel("Heater/effective temperature parameter [C]")
    ax.set_ylabel("Normalized optical metric [1]")
    ax.set_title("(a) Optical response")
    ax.set_ylim(0.0, 1.10)
    ax.grid(alpha=0.25)
    ax.legend(fontsize=7.6, loc="lower right")

    ax = axes[1]
    ax.errorbar(
        measured_temperature,
        measured_sensitivity,
        yerr=measured_uncertainty,
        fmt="o",
        ms=5,
        capsize=3,
        color=colors["red"],
        label="PPT transcription",
        zorder=4,
    )
    ax.plot(
        temperature,
        model_sensitivity,
        "-s",
        ms=3.5,
        lw=1.7,
        color=colors["teal"],
        label="GABES conditional absolute",
    )
    ax.plot(
        temperature,
        normalized_model,
        "--",
        lw=1.6,
        color=colors["gray"],
        label=f"GABES shape, normalized at 20 C (x{room_scale:.2f})",
    )
    ax.axvspan(40.0, 45.0, color=colors["orange"], alpha=0.08)
    ax.annotate(
        "best PPT point",
        xy=(40.0, 6.68),
        xytext=(42.5, 11.5),
        fontsize=8,
        color=colors["red"],
        arrowprops={"arrowstyle": "->", "color": colors["red"], "lw": 1.0},
    )
    ax.text(
        0.98,
        0.20,
        "lower is better",
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        fontsize=8,
        color="#334155",
    )
    ax.set_xlabel("Heater setpoint [C]")
    ax.set_ylabel(r"RF sensitivity [nV cm$^{-1}$ Hz$^{-1/2}$]")
    ax.set_title("(b) RF sensitivity")
    ax.set_ylim(bottom=0.0)
    ax.grid(alpha=0.25)
    ax.legend(fontsize=7.4, loc="upper left")

    fig.suptitle(
        "Rydberg cell-heating comparison: GABES model vs 2026-07-15 PPT",
        fontsize=12,
        fontweight="bold",
    )
    fig.tight_layout()
    output_stem.parent.mkdir(parents=True, exist_ok=True)
    for suffix in (".pdf", ".png"):
        fig.savefig(
            output_stem.with_suffix(suffix),
            dpi=220 if suffix == ".png" else None,
            bbox_inches="tight",
        )
    plt.close(fig)

    best_measured_index = int(np.argmin(measured_sensitivity))
    best_model_index = int(np.argmin(model_sensitivity))
    return {
        "room_normalization_factor": float(room_scale),
        "best_ppt_temperature_c": float(
            measured_temperature[best_measured_index]
        ),
        "best_ppt_sensitivity_nv_cm_sqrt_hz": float(
            measured_sensitivity[best_measured_index]
        ),
        "best_model_temperature_c": float(temperature[best_model_index]),
        "best_model_sensitivity_nv_cm_sqrt_hz": float(
            model_sensitivity[best_model_index]
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--results", type=Path, required=True)
    parser.add_argument("--output-stem", type=Path, required=True)
    args = parser.parse_args()
    summary = build_figure(args.results.resolve(), args.output_stem.resolve())
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
