"""Retrospective structural diagnostics of saved simulations; no new solves."""
from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent
S = "sensitivity_nv_cm_sqrt_hz"


def read_json(name):
    return json.loads((ROOT / name).read_text(encoding="utf-8"))


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def curve_point(row):
    return {key: row[key] for key in ("probe_power_uw", "lo_rabi_mhz",
        "probe_detuning_mhz", S, "transmission")}


def main():
    result = read_json("results.json")
    interpretation = read_json("interpretation.json")
    temperatures = result["config"]["temperature_c"]
    powers = result["config"]["coupling_power_mw"]
    density = {row["temperature_c"]: row["model_density_m3"]
               for row in interpretation["model_temperature_diagnostics"]}
    output = {
        "origin": "⟦AI⟧ Retrospective diagnostics of saved simulation samples",
        "adoption_status": "AI analysis only; not a human-adopted goal or experimental validation",
        "scope": "No new physical solves. Structural reduction applies to this profiled conditional objective and bounded domain only.",
        "separability": [], "nominal_argmin_counts": [], "cost": {},
    }
    for mode in result["config"]["modes"]:
        for opacity in result["config"]["opacities"]:
            rows = [row for row in result["optima"] if row["doppler"] == mode
                    and row["opacity"] == opacity]
            values = {(row["temperature_c"], row["coupling_power_mw"]): row[S] for row in rows}
            matrix = np.array([[values[t, pc] for pc in powers] for t in temperatures])
            logarithm = np.log(matrix)
            additive = logarithm.mean(axis=1)[:, None] + logarithm.mean(axis=0)[None, :] - logarithm.mean()
            additive_errors = np.exp(additive) / matrix - 1
            holdout = []
            for temperature in temperatures:
                for power in powers:
                    if temperature == 20. or power == 30.:
                        continue
                    predicted = values[temperature, 30.] * values[20., power] / values[20., 30.]
                    holdout.append({"temperature_c": temperature, "coupling_power_mw": power,
                        "saved_sensitivity": values[temperature, power], "predicted_sensitivity": predicted,
                        "relative_error": predicted / values[temperature, power] - 1})
            errors = np.array([row["relative_error"] for row in holdout])
            normalized = []
            for power in powers:
                sn = np.array([values[t, power] * density[t] for t in temperatures])
                normalized.append({"coupling_power_mw": power,
                    "max_to_min_SN_ratio_15_to_35c": float(sn.max() / sn.min()),
                    "SN35_to_SN20_ratio": values[35., power] * density[35.] / (values[20., power] * density[20.])})
            output["separability"].append({
                "doppler": mode, "opacity": opacity,
                "in_sample_log_additive": {
                    "formula": "log S_hat(T,Pc) = mean_Pc(log S) + mean_T(log S) - mean_all(log S)",
                    "max_absolute_relative_error": float(np.abs(additive_errors).max()),
                    "rms_relative_error": float(np.sqrt(np.mean(additive_errors ** 2))),
                    "validation_status": "In-sample structural diagnostic; not a holdout test"},
                "gold_cross_holdout": {
                    "formula": "S_hat(T,Pc) = S(T,30 mW)*S(20 C,Pc)/S(20 C,30 mW)",
                    "training_coordinate_count": len(temperatures) + len(powers) - 1,
                    "heldout_coordinate_count": len(holdout),
                    "coordinate_reduction_fraction_before_extra_validation": 1 - (len(temperatures) + len(powers) - 1) / len(rows),
                    "max_absolute_relative_error": float(np.abs(errors).max()),
                    "rms_relative_error": float(np.sqrt(np.mean(errors ** 2))),
                    "validation_status": "Retrospective held-out map coordinates within the same simulation campaign; not independent experimental evidence",
                    "heldout_rows": holdout},
                "known_density_baseline": {
                    "formula": "SN(T,Pc) = S(T,Pc)*N_model(T)",
                    "scope": "Removes the leading inverse-density trend for diagnosis; does not assume exact S proportional to 1/N or linear temperature dependence",
                    "per_coupling_power": normalized}})
            if opacity == .001:
                def counts(key):
                    unique = sorted({row[key] for row in rows})
                    return [{"value": value, "count": sum(row[key] == value for row in rows)} for value in unique]
                output["nominal_argmin_counts"].append({"doppler": mode, "coordinate_count": len(rows),
                    "probe_power_uw": counts("probe_power_uw"), "lo_rabi_mhz": counts("lo_rabi_mhz"),
                    "probe_detuning_mhz": counts("probe_detuning_mhz"),
                    "interpretation_limit": "Counts of joint saved optima do not imply monotonic one-variable response at fixed remaining controls"})

    prior_name = "../rydberg_gold_map_20261103/results.json"
    prior = read_json(prior_name)
    output["prior_pilot_cross_diagnostic"] = []
    for mode in result["config"]["modes"]:
        values = {(row["temperature_c"], row["coupling_power_mw"]): row[S]
                  for row in prior["optima"] if row["doppler"] == mode and row["opacity"] == .001}
        errors = [values[t, 30.] * values[20., pc] / values[20., 30.] / value - 1
                  for (t, pc), value in values.items() if t != 20. and pc != 30.]
        output["prior_pilot_cross_diagnostic"].append({
            "doppler": mode, "opacity": .001, "training_coordinate_count": 8,
            "heldout_coordinate_count": len(errors),
            "max_absolute_relative_error": float(np.max(np.abs(errors))),
            "rms_relative_error": float(np.sqrt(np.mean(np.asarray(errors) ** 2))),
            "scope": "Retrospective check of whether the earlier pilot already contained a separability signal; no prospective experimental validation"})

    with (ROOT / "cells.jsonl").open(encoding="utf-8") as stream:
        for line in stream:
            cell = json.loads(line)
            entry = output["cost"].setdefault(cell["task"][0], {"map_cells": 0, "control_states": 0, "worker_seconds": 0.})
            entry["map_cells"] += 1
            entry["control_states"] += cell["control_evaluations"]
            entry["worker_seconds"] += cell["elapsed_s"]
    total_work = sum(entry["worker_seconds"] for entry in output["cost"].values())
    for entry in output["cost"].values():
        entry["worker_time_fraction"] = entry["worker_seconds"] / total_work
    output["cost_scope"] = "Worker elapsed times from saved cells; fractions are not isolated CPU instruction counts or guaranteed prospective wall-time savings"

    with (ROOT / "control_points.csv").open(encoding="utf-8-sig", newline="") as stream:
        controls = list(csv.DictReader(stream))
    for row in controls:
        for key in ("temperature_c", "coupling_power_mw", "probe_power_uw", "lo_rabi_mhz",
                    "opacity", "probe_detuning_mhz", S, "transmission"):
            row[key] = float(row[key])
    curves = []
    for mode in result["config"]["modes"]:
        cell = [row for row in controls if row["doppler"] == mode and row["opacity"] == .001
                and row["temperature_c"] == 20. and row["coupling_power_mw"] == 30.]
        fixed, profiled = [], []
        for probe in (2., 4., 6., 8., 10.):
            fixed.append(curve_point(next(row for row in cell if row["probe_power_uw"] == probe and row["lo_rabi_mhz"] == 3.7)))
            profiled.append(curve_point(min((row for row in cell if row["probe_power_uw"] == probe), key=lambda row: row[S])))
        curves.append({"doppler": mode, "temperature_c": 20., "coupling_power_mw": 30., "opacity": .001,
            "fixed_gold_LO_3p7MHz": fixed,
            "LO_profiled_over_saved_samples": profiled,
            "detuning_scope": "Each saved row already minimizes its evaluated signed finite-IF detuning candidates",
            "LO_profile_scope": "Minimum over saved LO samples at each fixed Pp; unequal adaptive sampling, not continuous LO global optimization",
            "fixed_LO_probe_minimum": min(fixed, key=lambda row: row[S]),
            "saved_joint_probe_LO_minimum": curve_point(min(cell, key=lambda row: row[S]))})
    output["gold_control_confounding_diagnostic"] = curves
    sources = ["results.json", "interpretation.json", "control_points.csv", "cells.jsonl", "meta_diagnostics.py", prior_name]
    output["source_sha256"] = {name: sha256(ROOT / name) for name in sources}
    path = ROOT / "meta_diagnostics.json"
    path.write_text(json.dumps(output, ensure_ascii=False, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({"output": path.name, "scenario_diagnostics": len(output["separability"]),
        "cost": output["cost"], "gold_control_confounding_diagnostic": curves}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
