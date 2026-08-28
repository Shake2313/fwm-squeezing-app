"""Overlay GABES temperature sweep against the group's paper Fig. 5(a) landmarks.

Uses only real GABES output (from results.json) plus the paper's *stated* landmark
temperatures (PSN optimum ~33 C, Rb vapour-pressure discontinuity 39.30 C) and the
2026-07-15 transcribed PSN points. The paper's full curve is NOT reconstructed; only its
reported optimum position is annotated.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any, Mapping

import numpy as np


def _mv(metrics: Mapping[str, Any], name: str) -> float:
    v = metrics[name]
    return float(v["value"] if isinstance(v, dict) else v)


# Paper landmarks (2606.04354v1, text of Sec. 5 / Fig. 5a) and 2026-07-15 transcription.
PAPER_PSN_OPTIMUM_C = 33.0            # "PSN limit exhibits a minimum around 33 C"
RB_MELTING_C = 39.30                  # vapour-pressure discontinuity (solid->liquid)
EXP_BEST_C = 38.0                     # cell-heating best PSN point (40 C setpoint)
# actual-temperature, PSN limit [nV/cm/sqrt(Hz)] from the 2026-07-15 deck
EXP_T = np.array([20.0, 38.0, 42.5, 44.0])
EXP_PSN = np.array([11.2, 6.72, 17.6, 25.0])
EXP_ERR = np.array([0.4, 0.04, 0.5, 0.1])


def build_figure(results_path: Path, output_stem: Path) -> dict[str, float]:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    results = json.loads(results_path.read_text(encoding="utf-8"))
    at = results["model"]["sweeps"]["AT electrometry"]
    T = np.array([p["temperature_c"] for p in at], dtype=float)
    order = np.argsort(T)
    T = T[order]
    S = np.array([float(p["electrometry"]["field_sensitivity_nv_cm_sqrt_hz"]) for p in at])[order]
    Tp = np.array([_mv(p["metrics"], "transmission_at_resonance") for p in at])[order]
    slope = np.array([_mv(p["metrics"], "max_spectral_slope_per_mhz") for p in at])[order]

    gabes_opt_c = float(results["model"]["optima"]["min_total_field_sensitivity"]["temperature_c"])
    gabes_opt_S = float(results["model"]["optima"]["min_total_field_sensitivity"]["value"])
    Tp_at_gabes_opt = float(np.interp(gabes_opt_c, T, Tp))
    Tp_at_paper_opt = float(np.interp(PAPER_PSN_OPTIMUM_C, T, Tp))

    c = {"teal": "#0f766e", "red": "#b91c1c", "blue": "#2563eb",
         "orange": "#ea580c", "gray": "#64748b", "purple": "#7c3aed"}
    fig, axes = plt.subplots(1, 2, figsize=(11.6, 4.5))

    # ---- Panel (a): PSN-limited sensitivity vs temperature ----
    ax = axes[0]
    ax.plot(T, S, "-", color=c["teal"], lw=2.0, zorder=3,
            label="GABES PSN-limited $S$ (conditional)")
    ax.plot(gabes_opt_c, gabes_opt_S, "o", color=c["teal"], ms=8, zorder=5)
    ax.annotate(f"GABES optimum\n{gabes_opt_c:g} C", xy=(gabes_opt_c, gabes_opt_S),
                xytext=(gabes_opt_c + 1.5, gabes_opt_S + 0.9), fontsize=8.5, color=c["teal"],
                arrowprops={"arrowstyle": "->", "color": c["teal"], "lw": 1.0})
    ax.axvline(PAPER_PSN_OPTIMUM_C, color=c["purple"], ls="--", lw=1.6, zorder=2,
               label=f"paper Fig.5a optimum ~{PAPER_PSN_OPTIMUM_C:g} C")
    ax.axvline(EXP_BEST_C, color=c["red"], ls=":", lw=1.8, zorder=2,
               label=f"experiment best {EXP_BEST_C:g} C")
    ax.axvline(RB_MELTING_C, color=c["gray"], ls="-.", lw=1.2, zorder=1,
               label=f"Rb melting {RB_MELTING_C:g} C")
    ax.set_xlabel("Effective vapour temperature [C]")
    ax.set_ylabel(r"GABES $S$ [nV cm$^{-1}$ Hz$^{-1/2}$] (tech. noise = 0)")
    ax.set_title("(a) PSN-limited sensitivity: where is the optimum?")
    ax.set_ylim(0.0, min(6.5, np.nanmax(S) * 1.05))
    ax.grid(alpha=0.25)

    # secondary axis: transcribed experimental PSN (different absolute scale)
    ax2 = ax.twinx()
    ax2.errorbar(EXP_T, EXP_PSN, yerr=EXP_ERR, fmt="s", color=c["red"], ms=6, capsize=3,
                 zorder=6, label="2026-07-15 PPT PSN limit")
    ax2.set_ylabel("PPT PSN limit [nV cm$^{-1}$ Hz$^{-1/2}$]", color=c["red"])
    ax2.tick_params(axis="y", labelcolor=c["red"])
    ax2.set_ylim(0.0, 28.0)
    h1, l1 = ax.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    ax.legend(h1 + h2, l1 + l2, fontsize=7.3, loc="upper left")

    # ---- Panel (b): probe transmission (the 'why') ----
    ax = axes[1]
    ax.plot(T, Tp * 100.0, "-", color=c["blue"], lw=2.0, label="GABES resonant transmission")
    ax.axvline(gabes_opt_c, color=c["teal"], ls="--", lw=1.4)
    ax.axvline(PAPER_PSN_OPTIMUM_C, color=c["purple"], ls="--", lw=1.4)
    ax.axvline(RB_MELTING_C, color=c["gray"], ls="-.", lw=1.1)
    ax.plot(gabes_opt_c, Tp_at_gabes_opt * 100.0, "o", color=c["teal"], ms=8, zorder=5)
    ax.annotate(f"GABES optimum sits at\nonly {Tp_at_gabes_opt*100:.0f}% transmission",
                xy=(gabes_opt_c, Tp_at_gabes_opt * 100.0),
                xytext=(gabes_opt_c - 21.0, Tp_at_gabes_opt * 100.0 + 22.0),
                fontsize=8.5, color=c["teal"],
                arrowprops={"arrowstyle": "->", "color": c["teal"], "lw": 1.0})
    ax.plot(PAPER_PSN_OPTIMUM_C, Tp_at_paper_opt * 100.0, "o", color=c["purple"], ms=7, zorder=5)
    ax.annotate(f"paper optimum ~{Tp_at_paper_opt*100:.0f}%",
                xy=(PAPER_PSN_OPTIMUM_C, Tp_at_paper_opt * 100.0),
                xytext=(PAPER_PSN_OPTIMUM_C - 8.0, Tp_at_paper_opt * 100.0 + 20.0),
                fontsize=8.5, color=c["purple"],
                arrowprops={"arrowstyle": "->", "color": c["purple"], "lw": 1.0})
    ax.set_xlabel("Effective vapour temperature [C]")
    ax.set_ylabel("Resonant probe transmission [%]")
    ax.set_title("(b) Why GABES runs too hot: photon starvation")
    ax.set_ylim(0.0, 90.0)
    ax.grid(alpha=0.25)
    ax.legend(fontsize=7.6, loc="upper right")

    fig.suptitle("GABES temperature sweep vs paper Fig. 5(a) landmarks (2606.04354v1)",
                 fontsize=12, fontweight="bold")
    fig.tight_layout()
    output_stem.parent.mkdir(parents=True, exist_ok=True)
    for suffix in (".pdf", ".png"):
        fig.savefig(output_stem.with_suffix(suffix), dpi=220 if suffix == ".png" else None,
                    bbox_inches="tight")
    plt.close(fig)

    return {
        "gabes_optimum_c": gabes_opt_c,
        "gabes_optimum_transmission_pct": Tp_at_gabes_opt * 100.0,
        "paper_optimum_c": PAPER_PSN_OPTIMUM_C,
        "paper_optimum_transmission_pct": Tp_at_paper_opt * 100.0,
        "experiment_best_c": EXP_BEST_C,
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
