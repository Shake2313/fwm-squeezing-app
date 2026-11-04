"""Derive comparison and local-trend summaries from completed saved samples."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent
PRIOR = ROOT.parent / "rydberg_gold_map_20261103"
S = "sensitivity_nv_cm_sqrt_hz"
CONTROLS = ["temperature_c", "coupling_power_mw", "probe_power_uw", "lo_rabi_mhz", "probe_detuning_mhz"]


def compact(row):
    value = {key: row[key] for key in CONTROLS + [S, "opacity", "doppler", "transmission"]}
    value["model_optical_depth_at_selected_frequency"] = float(-np.log(row["transmission"]))
    return value


def main():
    from dense_study import FixedScanScheme, baseline_parameters

    result = json.loads((ROOT / "results.json").read_text(encoding="utf-8"))
    prior = json.loads((PRIOR / "results.json").read_text(encoding="utf-8"))
    summary = {"origin": "AI interpretation of saved samples; no interpolation or new physical assumptions",
               "mode_summaries": [], "scope": "Sampled conditional PSN model, not experimental calibration"}
    scheme = FixedScanScheme()
    temperature_diagnostics = []
    for temperature in result["config"]["temperature_c"]:
        params = dict(baseline_parameters(), temp_c=temperature)
        state = scheme._temperature_state(params)
        temperature_diagnostics.append({"temperature_c": temperature,
            "model_density_m3": scheme._local_density_from_temperature_state(state),
            "model_transit_dephasing_mhz": scheme._transit_rate_mhz(params)})
    summary["model_temperature_diagnostics"] = temperature_diagnostics
    density_at = {row["temperature_c"]: row["model_density_m3"] for row in temperature_diagnostics}
    summary["model_density_ratio_35c_to_20c"] = density_at[35.] / density_at[20.]
    for mode in result["config"]["modes"]:
        for opacity in result["config"]["opacities"]:
            old = [row for row in prior["optima"] if row["doppler"] == mode and row["opacity"] == opacity]
            new = [row for row in result["optima"] if row["doppler"] == mode and row["opacity"] == opacity]
            lookup = {(row["temperature_c"], row["coupling_power_mw"]): row for row in new}
            comparisons = [{"temperature_c": row["temperature_c"], "coupling_power_mw": row["coupling_power_mw"],
                "coarse": row[S], "dense": lookup[(row["temperature_c"], row["coupling_power_mw"])][S],
                "reduction_percent": 100 * (1 - lookup[(row["temperature_c"], row["coupling_power_mw"])][S] / row[S])}
                for row in old]
            old_best = min(old, key=lambda row: row[S])
            old_coordinates = {(row["temperature_c"], row["coupling_power_mw"]) for row in old}
            new_at_old = min((row for row in new if (row["temperature_c"], row["coupling_power_mw"]) in old_coordinates), key=lambda row: row[S])
            best = min(new, key=lambda row: row[S])
            entry = {"doppler": mode, "opacity": opacity,
                "coarse_best": compact(old_best), "dense_best_at_old_map_coordinates": compact(new_at_old),
                "dense_best": compact(best), "coarse_coordinate_comparisons": comparisons,
                "coarse_to_dense_percent": {
                    "min": min(row["reduction_percent"] for row in comparisons),
                    "median": float(np.median([row["reduction_percent"] for row in comparisons])),
                    "max": max(row["reduction_percent"] for row in comparisons)},
                "hidden_control_and_detuning_refinement_gain_factor": old_best[S] / new_at_old[S],
                "map_coordinate_refinement_gain_factor": new_at_old[S] / best[S],
                "near_best_map_coordinates": [], "temperature_trends": [], "coupling_trends": []}
            for tolerance in [.01, .03]:
                inside = [row for row in new if row[S] <= best[S] * (1 + tolerance)]
                entry["near_best_map_coordinates"].append({"relative_tolerance": tolerance,
                    "count": len(inside),
                    "temperature_c_range": [min(row["temperature_c"] for row in inside), max(row["temperature_c"] for row in inside)],
                    "coupling_power_mw_range": [min(row["coupling_power_mw"] for row in inside), max(row["coupling_power_mw"] for row in inside)],
                    "coordinates": [[row["temperature_c"], row["coupling_power_mw"]] for row in inside]})
            for axis, other, field in [("temperature_c", "coupling_power_mw", "temperature_trends"),
                                       ("coupling_power_mw", "temperature_c", "coupling_trends")]:
                slopes = []
                for fixed in sorted({row[other] for row in new}):
                    line = sorted((row for row in new if row[other] == fixed), key=lambda row: row[axis])
                    values = np.array([row[S] for row in line])
                    points = np.array([row[axis] for row in line])
                    slopes.extend((np.diff(np.log(values)) / np.diff(points)).tolist())
                    minimum = min(line, key=lambda row: row[S])
                    entry[field].append({other: fixed, "best_coordinate": minimum[axis],
                        "minimum_sensitivity": minimum[S],
                        "increasing_segments": int(np.sum(np.diff(values) > 1e-10 * values[:-1])),
                        "segment_count": len(line) - 1})
                entry[field + "_log_derivative_range"] = [min(slopes), max(slopes)]
            summary["mode_summaries"].append(entry)
    nominal = {entry["doppler"]: entry for entry in summary["mode_summaries"] if entry["opacity"] == .001}
    summary["nominal_mode_gap_at_separately_optimized_bests"] = nominal["on"]["dense_best"][S] / nominal["off"]["dense_best"][S]
    summary["sha256"] = {str(path.relative_to(ROOT.parent)): hashlib.sha256(path.read_bytes()).hexdigest()
                         for path in [ROOT / "results.json", PRIOR / "results.json", Path(__file__)]}
    (ROOT / "interpretation.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    for entry in summary["mode_summaries"]:
        if entry["opacity"] == .001:
            print(json.dumps({key: entry[key] for key in ["doppler", "dense_best", "coarse_to_dense_percent", "near_best_map_coordinates"]}, ensure_ascii=False))


if __name__ == "__main__":
    main()
