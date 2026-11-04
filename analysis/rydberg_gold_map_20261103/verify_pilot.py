"""Verify completed pilot exports and selected numerical reference points.

Run after pilot.py finishes. This verifier leaves the campaign checkpoint and
production solver untouched. Dense thermal checks use existing response APIs.
"""
from __future__ import annotations

import os

os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"

import argparse
import csv
import hashlib
import itertools
import json
import math
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parents[1]
sys.path.insert(0, str(REPO))

from gabes import constants, species
from gabes.schemes.rydberg import RydbergEITScheme


def read_json(path):
    return json.loads(path.read_text(encoding="utf-8"))


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def control_key(row):
    return (row["doppler"], row["temperature_c"], row["coupling_power_mw"],
            row["probe_power_uw"], row["lo_rabi_mhz"])


def cell_key(row):
    return (row["doppler"], row["opacity"], row["temperature_c"], row["coupling_power_mw"])


def read_controls(path):
    with path.open(encoding="utf-8-sig", newline="") as stream:
        rows = list(csv.DictReader(stream))
    numeric = (
        "temperature_c", "coupling_power_mw", "probe_power_uw", "lo_rabi_mhz", "opacity",
        "probe_detuning_mhz", "sensitivity_nv_cm_sqrt_hz", "psn_nv_cm_sqrt_hz",
        "transmission", "current_response_a_per_v_m", "detected_power_w",
    )
    for row in rows:
        for key in numeric:
            row[key] = float(row[key])
        row["candidate_count"] = int(row["candidate_count"])
        row["refined"] = row["refined"] == "True"
    return rows


def thermal_reference(row, window):
    """Check an independent global grid plus a denser neighborhood of its winner."""
    scheme = RydbergEITScheme()
    params = scheme._defaults("AT electrometry")
    params.update(doppler="on", temp_c=row["temperature_c"],
                  coupling_power_mw=row["coupling_power_mw"],
                  probe_power_uw=row["probe_power_uw"], lo_rabi_mhz=row["lo_rabi_mhz"],
                  effective_opacity_participation=row["opacity"])
    raw = scheme.compute(params)
    x, transmission, _ = scheme._transmission(raw["chi_bar"], raw, params)
    center = row["probe_detuning_mhz"]
    targets = np.r_[np.linspace(-window, window, 41),
                    np.linspace(max(-window, center - 0.4), min(window, center + 0.4), 41),
                    center]
    indices = np.unique([int(np.argmin(np.abs(x - target))) for target in targets])
    response = scheme._superheterodyne_atomic_phasors(
        raw, params, x, transmission, static_indices=indices)
    result = scheme._superheterodyne_readout(raw, params, x, transmission, atomic_response=response)
    value = float(result["sensitivity"].total_field_asd_nv_cm_per_sqrt_hz)
    return {
        "controls": control_key(row), "opacity": row["opacity"],
        "saved_sensitivity": row["sensitivity_nv_cm_sqrt_hz"],
        "dense_sensitivity": value,
        "relative_improvement": 1.0 - value / row["sensitivity_nv_cm_sqrt_hz"],
        "saved_probe_detuning_mhz": center,
        "dense_probe_detuning_mhz": float(result["optimum_detuning_mhz"]),
        "candidate_count": len(indices),
        "search_scope": "41 global targets plus 41 local targets, snapped to the existing spectrum grid",
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--thermal-checks", type=int, choices=(0, 1, 2), default=2)
    args = parser.parse_args()
    if not (ROOT / "results.json").is_file():
        raise SystemExit("Completed results.json is required. Wait for pilot.py to finish.")
    config = read_json(ROOT / "config.json")
    summary = read_json(ROOT / "results.json")
    manifest = read_json(ROOT / "source_manifest.json")
    rows = read_controls(ROOT / "control_points.csv")
    checks, warnings = [], []

    def check(name, passed, detail=None):
        checks.append({"id": name, "status": "PASS" if passed else "FAIL", "detail": detail})

    check("completed_config_matches", summary["config"] == config)
    saved_hashes = manifest["code_sha256"]
    current_hashes = {name: digest(ROOT / name if name == "pilot.py" else REPO / name)
                      for name in saved_hashes}
    current_signature = hashlib.sha256(json.dumps(current_hashes, sort_keys=True).encode()).hexdigest()
    check("engine_and_pilot_source_hashes", current_hashes == saved_hashes)
    check("configuration_signature", config["signature"] == manifest["signature"] == current_signature)
    check("current_environment_matches_final_manifest",
          manifest["python"] == sys.version and manifest["numpy"] == np.__version__)

    coarse = [row for row in rows if not row["refined"]]
    expected = set(itertools.product(config["modes"], config["temperature_c"],
        config["coupling_power_mw"], config["probe_power_uw"], config["lo_rabi_mhz"]))
    check("complete_coarse_control_grid", set(map(control_key, coarse)) == expected,
          {"expected": len(expected), "actual": len(set(map(control_key, coarse)))})
    check("complete_unique_opacity_scenarios", len(coarse) == len(expected) * len(config["opacities"])
          and all(sorted(row["opacity"] for row in coarse if control_key(row) == key)
                  == sorted(config["opacities"]) for key in expected))
    numeric_fields = ("sensitivity_nv_cm_sqrt_hz", "psn_nv_cm_sqrt_hz", "transmission",
                      "current_response_a_per_v_m", "detected_power_w")
    check("positive_finite_outputs", all(math.isfinite(row[key]) and row[key] > 0
          for row in rows for key in numeric_fields))
    check("bounded_probe_search", all(abs(row["probe_detuning_mhz"])
          <= config["probe_window_mhz"] + 0.03 for row in rows))
    check("passive_transmission", all(row["transmission"] <= 1.0 + 1e-12 for row in rows))
    check("psn_equals_total_under_zero_technical_noise", all(math.isclose(
          row["psn_nv_cm_sqrt_hz"], row["sensitivity_nv_cm_sqrt_hz"], rel_tol=1e-12) for row in rows))

    # Independent unit/noise reconstruction for the declared equal-arm detector.
    eta = RydbergEITScheme._REF["detector_quantum_efficiency"]
    optical_frequency = species.RB85.line("D2")[1]
    responsivity = eta * constants.ELEMENTARY_CHARGE / (2 * np.pi * constants.HBAR * optical_frequency)
    noise_errors = []
    for row in rows:
        current = responsivity * row["detected_power_w"]
        shot_asd = math.sqrt(4 * constants.ELEMENTARY_CHARGE * current)
        expected_psn = shot_asd / row["current_response_a_per_v_m"] * 1e7
        noise_errors.append(abs(expected_psn / row["psn_nv_cm_sqrt_hz"] - 1.0))
    check("independent_balanced_psn_and_field_unit_reconstruction", max(noise_errors) < 1e-11,
          {"max_relative_error": max(noise_errors)})

    selected = summary["optima"]
    expected_cells = set(itertools.product(config["modes"], config["opacities"],
                        config["temperature_c"], config["coupling_power_mw"]))
    check("complete_unique_map_cells", set(map(cell_key, selected)) == expected_cells
          and len(selected) == len(expected_cells))
    minimum_errors, fixed_control_ratios = [], []
    for winner in selected:
        group = [row for row in rows if cell_key(row) == cell_key(winner)]
        minimum = min(row["sensitivity_nv_cm_sqrt_hz"] for row in group)
        minimum_errors.append(abs(winner["sensitivity_nv_cm_sqrt_hz"] / minimum - 1.0))
        fixed = min(row["sensitivity_nv_cm_sqrt_hz"] for row in group
                    if row["probe_power_uw"] == 6.0 and row["lo_rabi_mhz"] == 3.7)
        fixed_control_ratios.append(winner["sensitivity_nv_cm_sqrt_hz"] / fixed)
    check("each_map_cell_is_saved_control_minimum", max(minimum_errors) < 1e-12)
    check("optimized_map_no_worse_than_fixed_gold_controls", max(fixed_control_ratios) <= 1.0 + 1e-12,
          {"maximum_optimized_to_fixed_ratio": max(fixed_control_ratios)})

    refinement_changes = []
    for row in rows:
        if not row["refined"]:
            continue
        original = next(item for item in coarse if control_key(item) == control_key(row)
                        and item["opacity"] == row["opacity"])
        refinement_changes.append(row["sensitivity_nv_cm_sqrt_hz"]
                                  / original["sensitivity_nv_cm_sqrt_hz"] - 1.0)
    check("refinement_retains_or_improves_coarse_candidates", bool(refinement_changes)
          and max(refinement_changes) < 1e-4,
          {"largest_relative_worsening": max(refinement_changes, default=None)})

    scheme = RydbergEITScheme()
    params = scheme._defaults("AT electrometry")
    raw = scheme.compute(params)
    fresh = scheme._readout(raw, params)["superhet"]["sensitivity"].total_field_asd_nv_cm_per_sqrt_hz
    saved_gold = next(row for row in summary["baseline"] if row["doppler"] == "off" and row["opacity"] == 0.001)
    check("fresh_default_gold_parity", math.isclose(fresh, saved_gold["sensitivity_nv_cm_sqrt_hz"], rel_tol=1e-10)
          and math.isclose(fresh, 6.0823872151659275, rel_tol=1e-10), {"fresh_sensitivity": fresh})

    prior = read_json(ROOT / "sources/prior_work.json")
    archives = []
    import pymupdf
    for source in prior["sources"]:
        if not source.get("archive"):
            continue
        archive = source["archive"]
        path = (ROOT / archive["relative_path"]).resolve()
        with pymupdf.open(path) as document:
            okay = path.is_relative_to(ROOT) and digest(path) == archive["sha256"]
            okay = okay and path.stat().st_size == archive["bytes"] and len(document) == archive["pdf_pages"]
        archives.append({"id": source["id"], "path": archive["relative_path"], "okay": okay})
    check("prior_work_archive_hash_size_pages", bool(archives) and all(item["okay"] for item in archives), archives)
    check("older_apparatus_not_substituted", prior["gold_baseline_policy"]["older_apparatus_substitution"] is False
          and config["fixed"]["cell_mm"] == 50 and config["fixed"]["beam_diameter_mm"] == 0.15)

    dense = []
    thermal_gold = next(row for row in summary["baseline"] if row["doppler"] == "on" and row["opacity"] == 0.001)
    thermal_best = next(row for row in summary["best_nominal"] if row["doppler"] == "on")
    numerical_rows = [thermal_gold, thermal_best][:args.thermal_checks]
    for row in numerical_rows:
        dense.append(thermal_reference(row, config["probe_window_mhz"]))
    if dense:
        check("thermal_dense_candidates_agree_with_saved_local_optima",
              all(abs(item["relative_improvement"]) < 0.02 for item in dense), dense)

    boundaries = []
    for row in summary["best_nominal"]:
        flags = row["boundary_flags"].split("|") if row["boundary_flags"] else []
        for name, grid_key in (("temperature_c", "temperature_c"), ("coupling_power_mw", "coupling_power_mw")):
            grid = config[grid_key]
            if row[name] == min(grid):
                flags.append(name + "_map_lower")
            if row[name] == max(grid):
                flags.append(name + "_map_upper")
        boundaries.append({"mode": row["doppler"], "flags": flags, "controls": control_key(row),
                           "sensitivity": row["sensitivity_nv_cm_sqrt_hz"]})
    if any(item["flags"] for item in boundaries):
        warnings.append("Best sampled settings touch map or control bounds; no unconstrained optimum established.")
    warnings += [
        "Resume parser does not recover a truncated final checkpoint JSONL line.",
        "Campaign signature covers code and inputs, not Python/NumPy; final-manifest environment parity cannot audit earlier resumed processes.",
        "Opacity, anchored Rabi rates, compact propagation and zero technical noise remain model assumptions.",
        "Thermal dense checks refine probe candidates only at selected controls; continuous control-variable convergence is untested.",
    ]
    report = dict(as_of=config["as_of"], origin="⟦AI⟧ Independent pilot verification",
        status="PASS" if all(item["status"] == "PASS" for item in checks) else "FAIL",
        checks=checks, warnings=warnings, best_boundary_flags=boundaries,
        gold_baselines=summary["baseline"], measured_gold_reference_nv_cm_sqrt_hz=12.5,
        scope="Saved sampled map, numerical consistency and source provenance; no experimental predictive validation.",
        artifact_sha256={name: digest(ROOT / name) for name in (
            "results.json", "control_points.csv", "optima.csv", "config.json", "source_manifest.json",
            "sources/prior_work.json", "verify_pilot.py")})
    (ROOT / "validation.json").write_text(json.dumps(report, ensure_ascii=False, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({"status": report["status"], "checks": len(checks),
        "failed": [item["id"] for item in checks if item["status"] == "FAIL"],
        "boundaries": boundaries, "dense_thermal": dense}, ensure_ascii=False, indent=2))
    return 0 if report["status"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
