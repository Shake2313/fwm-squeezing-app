"""
Figure: PPT-internal consistency diagnostic for the 2026-07-15 cell-heating deck.

IMPORTANT: this figure uses ONLY numbers transcribed from the PPT.
No GABES curve is re-plotted here (the model sweep lives in brief_summary.pdf),
so nothing in this figure is model-derived except two annotated optimum markers,
which are labelled GABES and read from generated_brief/results.json.

The PPT transcription is read from ppt_rf_sensitivity_2026-07-15.csv -- the same
file the workflow ingests -- so the deck numbers have a single source. Rows
without an annotated actual temperature (the deck annotates one only for the
heated setpoints) are treated as the room-temperature condition.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from .workflow import (
    _json_dump,
    _serialized_path,
    _sha256,
    _utc_now,
    load_clean_csv,
)

HERE = Path(__file__).resolve().parent
DEFAULT_CSV = HERE / "ppt_rf_sensitivity_2026-07-15.csv"
DEFAULT_RESULTS = HERE / "generated_brief" / "results.json"
DEFAULT_OUTPUT = HERE / "generated_brief" / "ppt_consistency.pdf"


def build_figure(csv_path: Path, results_path: Path, output_path: Path) -> None:
    # ---- PPT transcription (single source: the workflow's status=PPT CSV) ----
    rows = load_clean_csv(csv_path, "sensitivity", "PPT")
    setp: list[str] = []
    tact_list: list[float] = []
    for row in rows:
        actual = row.get("actual_temperature_c")
        if actual is None:
            setp.append("room")
            tact_list.append(float(row["temperature_c"]))
        else:
            setp.append(f"{row['temperature_c']:g} C")
            tact_list.append(float(actual))
    Tact = np.asarray(tact_list)
    sens = np.asarray([row["sensitivity_nv_cm_sqrt_hz"] for row in rows])
    serr = np.asarray(
        [row.get("uncertainty_nv_cm_sqrt_hz") or 0.0 for row in rows])
    psn = np.asarray([row["psn_limit_nv_cm_sqrt_hz"] for row in rows])
    ratio = sens / psn
    # the deck reports EIT as clearest at the 45 C setpoint
    clearest = next(
        (row for row in rows if row["temperature_c"] == 45.0), None)

    # ---- the only model-derived inputs: two GABES optimum markers ----------
    results = json.loads(results_path.read_text(encoding="utf-8"))
    optima = results["model"]["optima"]
    slope_opt = float(optima["max_spectral_slope"]["temperature_c"])
    sens_opt = float(optima["min_total_field_sensitivity"]["temperature_c"])

    fig, axes = plt.subplots(1, 2, figsize=(9.4, 3.5))

    # ---- panel (a) ---------------------------------------------------------
    ax = axes[0]
    if clearest is not None:
        clearest_actual = float(clearest["actual_temperature_c"])
        ax.axvspan(
            min(slope_opt, clearest_actual) - 0.5,
            max(slope_opt, clearest_actual) + 0.5,
            color="#0F766E", alpha=0.10, zorder=0)
    ax.axvline(slope_opt, color="#0F766E", ls=":", lw=1.3, zorder=1,
               label=f"GABES slope opt. {slope_opt:.1f} C")
    ax.axvline(sens_opt, color="#94a3b8", ls="--", lw=1.2, zorder=1,
               label=f"GABES cond. sens. opt. {sens_opt:g} C")

    o = np.argsort(Tact)
    ax.errorbar(Tact[o], sens[o], yerr=serr[o], fmt="o-", color="#C2410C",
                ms=6, lw=1.6, capsize=3, zorder=3, label="measured sensitivity [PPT]")
    ax.plot(Tact[o], psn[o], "s--", color="#17324D", ms=5, lw=1.4, mfc="white",
            zorder=3, label="PSN limit quoted in PPT")
    if "room" in setp:
        i_room = setp.index("room")
        ax.plot([Tact[i_room]], [sens[i_room]], "o", ms=13, mfc="none",
                mec="#C2410C", mew=1.4, zorder=4)
        ax.annotate("room point:\nlikely earlier campaign",
                    xy=(Tact[i_room], sens[i_room]),
                    xytext=(Tact[i_room] + 1.5, 25.5), fontsize=7.5,
                    color="#C2410C",
                    arrowprops=dict(arrowstyle="->", color="#C2410C", lw=0.9))
    if clearest is not None:
        ax.annotate(f"clearest EIT (actual {clearest_actual:g} C)",
                    xy=(clearest_actual, sens[setp.index("45 C")]),
                    xytext=(32.0, 21.0), fontsize=7.5, color="#0F766E",
                    ha="center", va="center",
                    arrowprops=dict(arrowstyle="->", color="#0F766E", lw=0.8))

    setpoint_offsets = {"40 C": (0, -15), "45 C": (-26, 4), "50 C": (24, 2)}
    for lab, x, y in zip(setp, Tact, sens):
        if lab != "room":
            ax.annotate(f"setpoint {lab}", xy=(x, y),
                        xytext=setpoint_offsets.get(lab, (0, -14)),
                        textcoords="offset points", ha="center", fontsize=7,
                        color="#475569",
                        arrowprops=dict(arrowstyle="-", color="#94a3b8", lw=0.6))
    ax.set_xlabel("Actual / effective temperature from PPT PSN note [C]")
    ax.set_ylabel("E-field sens. [nV cm$^{-1}$ Hz$^{-1/2}$]")
    ax.set_title("(a) sensitivity and its own PSN floor track each other", fontsize=9)
    ax.set_xlim(18, 49)
    ax.set_ylim(0, 34)
    ax.legend(fontsize=7.5, loc="upper left", framealpha=0.9)
    ax.grid(alpha=0.25, lw=0.6)

    # ---- panel (b) ---------------------------------------------------------
    ax = axes[1]
    xpos = np.arange(len(rows))
    i_best = int(np.argmin(sens))
    cols = ["#94a3b8" if lab == "room" else
            ("#C2410C" if i == i_best else "#17324D")
            for i, lab in enumerate(setp)]
    ax.bar(xpos, ratio, color=cols, width=0.55, zorder=3)
    ax.axhline(1.0, color="#0F766E", lw=1.4, ls="--", zorder=4)
    # place the reference-line caption above the bar-value labels to avoid overlap
    ax.text(len(rows) - 0.62, 1.20, "measured = own PSN floor", fontsize=7.5,
            color="#0F766E", ha="right", va="bottom")
    for x, r in zip(xpos, ratio):
        ax.text(x, r + 0.03, f"{r:.2f}", ha="center", fontsize=8)
    ax.set_xticks(xpos)
    ax.set_xticklabels([f"{s}\n({t:g} C act.)" if s != "room" else "room\n(n/a)"
                        for s, t in zip(setp, Tact)], fontsize=8)
    ax.set_ylabel("measured / PSN limit  [1]")
    ax.set_ylim(0, 1.45)
    ax.set_title("(b) no excess technical noise at any setpoint", fontsize=9)
    ax.grid(axis="y", alpha=0.25, lw=0.6)

    fig.suptitle("PPT-internal consistency check (2026-07-15 deck numbers only)",
                 fontsize=10, y=1.005)
    fig.tight_layout()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, bbox_inches="tight")
    plt.close(fig)
    print("ratios:", np.round(ratio, 3))
    annotated = [i for i, lab in enumerate(setp) if lab != "room"]
    if annotated:
        i_lo = min(annotated, key=lambda i: Tact[i])
        i_hi = max(annotated, key=lambda i: Tact[i])
        print(f"PSN floor {Tact[i_lo]:g}->{Tact[i_hi]:g} C degradation:",
              round(float(psn[i_hi] / psn[i_lo]), 2))


def register_in_manifest(output_path: Path) -> None:
    """Record this script and its figure (hash + generation time)."""
    manifest_path = output_path.parent / "source_manifest.json"
    if not manifest_path.is_file():
        print(f"manifest not found, skipping registration: {manifest_path}")
        return
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    script_path = Path(__file__).resolve()
    code_entry = {
        "path": _serialized_path(script_path),
        "sha256": _sha256(script_path),
        "bytes": script_path.stat().st_size,
    }
    manifest["code"] = [
        entry for entry in manifest.get("code", [])
        if entry.get("path") != code_entry["path"]
    ] + [code_entry]
    artifact_entry: dict[str, Any] = {
        "path": _serialized_path(output_path),
        "format": "pdf",
        "description": "PPT-internal consistency diagnostic (deck transcription only)",
        "status": "PPT",
        "sha256": _sha256(output_path),
        "generated_at_utc": _utc_now(),
    }
    manifest["artifacts"] = [
        entry for entry in manifest.get("artifacts", [])
        if entry.get("path") != artifact_entry["path"]
    ] + [artifact_entry]
    _json_dump(manifest_path, manifest)
    print(f"registered in {manifest_path}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--csv", type=Path, default=DEFAULT_CSV)
    parser.add_argument("--results", type=Path, default=DEFAULT_RESULTS)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    build_figure(args.csv.resolve(), args.results.resolve(), args.output.resolve())
    register_in_manifest(args.output.resolve())


if __name__ == "__main__":
    main()
