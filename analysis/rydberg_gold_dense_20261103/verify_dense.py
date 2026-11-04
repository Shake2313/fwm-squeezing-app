"""Independent saved-grid checks and opt-in local/velocity convergence probes.

This verifier never modifies the dense campaign, previous pilot or production
engine. Local checks test the adaptive search; velocity checks independently
refine finite-IF quadrature at a fixed control and detuning.
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
import platform
import sys
import time

import numpy as np
import scipy

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parents[1]
PRIOR = REPO / "analysis/rydberg_gold_map_20261103"
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(ROOT))

from gabes import affine_doppler, constants, core, rydberg_electrometry, species
from gabes.schemes.rydberg import MHZ, RydbergEITScheme


def read_json(path):
    return json.loads(path.read_text(encoding="utf-8"))


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def cell_key(row):
    return (row["doppler"], float(row["opacity"]), float(row["temperature_c"]),
            float(row["coupling_power_mw"]))


def control_key(row):
    return (row["doppler"], float(row["temperature_c"]),
            float(row["coupling_power_mw"]), float(row["probe_power_uw"]),
            float(row["lo_rabi_mhz"]))


def read_controls(path):
    with path.open(encoding="utf-8-sig", newline="") as stream:
        rows = list(csv.DictReader(stream))
    numeric = ("temperature_c", "coupling_power_mw", "probe_power_uw", "lo_rabi_mhz",
               "opacity", "probe_detuning_mhz", "sensitivity_nv_cm_sqrt_hz",
               "psn_nv_cm_sqrt_hz", "transmission", "current_response_a_per_v_m",
               "detected_power_w", "central_sensitivity_nv_cm_sqrt_hz",
               "central_detuning_mhz", "sideband_sensitivity_nv_cm_sqrt_hz",
               "sideband_detuning_mhz", "detuning_step_mhz")
    for row in rows:
        for key in numeric:
            if key in row and row[key] not in ("", None):
                row[key] = float(row[key])
    return rows


def selected_cells(summary, maximum):
    selected = summary["optima"]
    cells = []
    if maximum >= 1:
        cells.extend(row for row in selected if row["opacity"] == .001
                     and row["temperature_c"] == 20 and row["coupling_power_mw"] == 30)
    if maximum >= 2:
        cells.extend(summary["best_nominal"])
    if maximum >= 3:
        # Most nearly tied saved central and sideband branches, away from bounds.
        candidates = [row for row in summary["branchgaps"] if row["doppler"] == "off"
                      and row["opacity"] == .001 and 15 < row["temperature_c"] < 35
                      and 20 < row["coupling_power_mw"] < 40]
        def branch_gap(row):
            center = row.get("central_sensitivity_nv_cm_sqrt_hz", math.inf)
            side = row.get("sideband_sensitivity_nv_cm_sqrt_hz", math.inf)
            return abs(math.log(center / side)) if center > 0 and side > 0 else math.inf
        if candidates:
            key = cell_key(min(candidates, key=branch_gap))
            cells.append(next(row for row in selected if cell_key(row) == key))
    unique = {}
    for row in cells:
        unique[cell_key(row)] = row
    return list(unique.values())


def local_rectangle_reference(row):
    """Exhaustive rectangle around saved controls, on a finer detuning grid."""
    from dense_study import evaluate_control

    pp = float(row["probe_power_uw"])
    lo = float(row["lo_rabi_mhz"])
    probes = sorted({round(float(np.clip(pp + offset, 2, 10)), 10)
                     for offset in (-.5, 0., .5)})
    lattice = sorted(set(np.arange(2, 6.0001, .125).tolist() + [3.7]))
    lower = [value for value in lattice if value < lo - 1e-9]
    upper = [value for value in lattice if value > lo + 1e-9]
    los = ([max(lower)] if lower else []) + [lo] + ([min(upper)] if upper else [])
    records = []
    started = time.perf_counter()
    for power, rabi in itertools.product(probes, los):
        found = evaluate_control(row["doppler"], row["temperature_c"],
            row["coupling_power_mw"], power, rabi, detuning_step=.005,
            thermal_target_step=.25, thermal_refine_step=.025)
        records.extend(found)
    candidates = [candidate for candidate in records if candidate["opacity"] == row["opacity"]]
    winner = min(candidates, key=lambda item: item["sensitivity_nv_cm_sqrt_hz"])
    result = {"saved_control": control_key(row), "opacity": row["opacity"],
        "saved_sensitivity": row["sensitivity_nv_cm_sqrt_hz"],
        "reference_control": control_key(winner),
        "reference_sensitivity": winner["sensitivity_nv_cm_sqrt_hz"],
        "reference_detuning_mhz": winner["probe_detuning_mhz"],
        "relative_improvement": 1 - winner["sensitivity_nv_cm_sqrt_hz"]
                                / row["sensitivity_nv_cm_sqrt_hz"],
        "probe_grid_uw": probes, "lo_grid_mhz": los,
        "control_evaluations": len(probes) * len(los),
        "elapsed_s": time.perf_counter() - started,
        "scope": "Adjacent final-grid 3x3 control rectangle (clipped at bounds, protected LO3.7 if adjacent); 0.005 MHz fixed detuning grid, thermal targets 0.25/0.025 MHz"}
    return result, records


def thermal_quadrature_reference(row):
    """Direct complex velocity sums, independent of the runner's candidate search.

    The initial spacing follows the narrowest pole over the fixed +/-6 MHz scan.
    It is then halved twice. This directly tests the finite-IF response, which
    need not inherit the stationary susceptibility's quadrature error bound.
    """
    started = time.perf_counter()
    class ValidationScan(RydbergEITScheme):
        def _scan(self, params):
            return np.linspace(-6, 6, 1201) * MHZ

    scheme = ValidationScan()
    params = read_json(PRIOR / "resolved_parameters.json")["baseline_parameters"]
    params.update(doppler="on", temp_c=row["temperature_c"],
        coupling_power_mw=row["coupling_power_mw"], probe_power_uw=row["probe_power_uw"],
        lo_rabi_mhz=row["lo_rabi_mhz"], effective_opacity_participation=row["opacity"])
    raw = scheme.compute(params)
    x, transmission, _ = scheme._transmission(raw["chi_bar"], raw, params)
    ground = (raw["rydberg_dephasing_mhz"] + raw["temperature_dephasing_mhz"]
              + raw["density_dephasing_mhz"] + raw["transit_mhz"]) * MHZ
    atom = scheme._atom(ground, raw["rf_dephasing_mhz"] * MHZ)
    detuning = float(row["probe_detuning_mhz"]) * MHZ
    fixed_index = int(np.argmin(abs(raw["scan"] - detuning)))
    if abs(raw["scan"][fixed_index] - detuning) / MHZ > 1e-8:
        raise ValueError("Saved detuning is not on independent fixed grid")
    hamiltonian = lambda value: scheme._hamiltonian_at_probe_detuning(raw, params, value)
    affine = affine_doppler.AffineScan(atom, hamiltonian, (1, 0), raw["probe_rabi_mhz"] * MHZ)
    velocity, weights = affine.quadrature_classes(raw["scan"], raw["T"],
        constants.MASS_85RB, raw["k_vec"])
    dv = float(velocity[1] - velocity[0])
    sigma = math.sqrt(constants.KB * raw["T"] / constants.MASS_85RB)
    base = core.build_liouvillian(hamiltonian(detuning), atom)
    operator = scheme._absorption_quadrature_operator()
    omega_if = float(params["if_khz"]) * 1e3 * 2 * np.pi

    def phasors(v):
        values = np.empty(len(v), dtype=complex)
        for first in range(0, len(v), 256):
            part = slice(first, first + 256)
            liouvillian = base[None] + (raw["k_vec"] * v[part])[:, None, None] * atom.S_v[None]
            response = rydberg_electrometry.weak_signal_response_from_liouvillian(
                liouvillian, 4, omega_if, signal_transition=(2, 3))
            values[part] = response.real_observable_phasor_per_angular_rabi(operator)
        return values

    static = complex(phasors(np.array([0.]))[0])
    variants = []
    for divisor, cutoff in ((1, 6), (2, 6), (4, 6), (4, 8)):
        if divisor == 1 and cutoff == 6:
            grid, weight = velocity, weights
        else:
            step = dv / divisor
            half = math.ceil(cutoff * sigma / step)
            grid = np.arange(-half, half + 1, dtype=float) * step
            weight = np.exp(-.5 * (grid / sigma) ** 2)
            weight /= weight.sum()
        response = phasors(grid)
        value = complex(np.dot(weight, response))
        incoherent_amplitude = float(np.dot(weight, np.abs(response)))
        readout = scheme._superheterodyne_readout(raw, params, x, transmission,
            atomic_response=(np.array([fixed_index]), np.array([value]),
                             "independent direct velocity sum at saved detuning"))
        variants.append({"spacing_divisor": divisor, "cutoff_sigma": cutoff,
            "velocity_classes": len(grid), "dv_m_s": dv / divisor,
            "phasor_real": value.real, "phasor_imag": value.imag,
            "phasor_magnitude": abs(value), "phase_deg": float(np.angle(value, deg=True)),
            "amplitude_participation_A": incoherent_amplitude / abs(static),
            "phase_coherence_C": abs(value) / incoherent_amplitude,
            "coherent_to_static_amplitude_ratio": abs(value) / abs(static),
            "sensitivity_nv_cm_sqrt_hz": float(readout["sensitivity"].total_field_asd_nv_cm_per_sqrt_hz)})
    reference = complex(variants[-1]["phasor_real"], variants[-1]["phasor_imag"])
    for variant in variants:
        value = complex(variant["phasor_real"], variant["phasor_imag"])
        variant["relative_complex_error_vs_dv4_8sigma"] = abs(value - reference) / abs(reference)
        variant["relative_sensitivity_change_vs_reference"] = abs(reference) / abs(value) - 1
    return {"control": control_key(row), "opacity": row["opacity"],
        "detuning_mhz": row["probe_detuning_mhz"], "variants": variants,
        "campaign_sensitivity": row["sensitivity_nv_cm_sqrt_hz"],
        "base_grid_relative_sensitivity_difference": variants[0]["sensitivity_nv_cm_sqrt_hz"]
            / row["sensitivity_nv_cm_sqrt_hz"] - 1,
        "static_phasor_real": static.real, "static_phasor_imag": static.imag,
        "elapsed_s": time.perf_counter() - started,
        "scope": "Fixed winning detuning; initial velocity grid from full fixed +/-6 MHz scan as in campaign. No control or detuning reoptimization."}


def export_csv(path, rows):
    if not rows:
        return
    with path.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--local-checks", type=int, choices=(0, 1, 2, 3), default=0,
                        help="0=none;1=Gold both modes;2=plus best both;3=plus off branch tie")
    parser.add_argument("--quadrature-checks", type=int, choices=(0, 1, 2), default=0,
                        help="Thermal Gold and then thermal best at fixed detuning")
    args = parser.parse_args()
    if not (ROOT / "results.json").is_file():
        raise SystemExit("Completed dense results.json required.")
    summary = read_json(ROOT / "results.json")
    config = read_json(ROOT / "config.json")
    manifest = read_json(ROOT / "source_manifest.json")
    rows = read_controls(ROOT / "control_points.csv")
    checks, warnings = [], []

    def check(identifier, okay, detail=None):
        checks.append({"id": identifier, "status": "PASS" if okay else "FAIL", "detail": detail})

    check("completed_config_matches", summary["config"] == config)
    saved_hashes = manifest["code_sha256"]
    current_hashes = {name: digest(ROOT / name if name == "dense_study.py" else REPO / name)
                      for name in saved_hashes}
    check("engine_inputs_and_dense_source_hashes", current_hashes == saved_hashes
          and config["code_sha256"] == saved_hashes)
    unsigned = {key: value for key, value in config.items() if key != "signature"}
    signature = hashlib.sha256(json.dumps(unsigned, sort_keys=True).encode()).hexdigest()
    check("configuration_signature", signature == config["signature"] == manifest["signature"])
    environment = dict(python=sys.version, numpy=np.__version__, scipy=scipy.__version__,
                       platform=platform.platform(), blas_threads=1)
    check("environment_is_resume_signed_and_matches", config["environment"] == environment
          and manifest["environment"] == environment)
    expected_cells = set(itertools.product(["off", "on"], [.0005, .001, .0015],
        [float(t) for t in range(15, 36)], [float(p) for p in np.arange(20, 40.1, 2.5)]))
    selected = summary["optima"]
    check("complete_unique_map_cells", set(map(cell_key, selected)) == expected_cells
          and len(selected) == len(expected_cells), {"expected": len(expected_cells), "actual": len(selected)})
    check("same_bounded_domain", min(r["temperature_c"] for r in selected) == 15
          and max(r["temperature_c"] for r in selected) == 35
          and min(r["coupling_power_mw"] for r in selected) == 20
          and max(r["coupling_power_mw"] for r in selected) == 40)
    check("gold_geometry_and_fixed_if_unchanged", config["fixed"] == dict(
        cell_mm=50., beam_diameter_mm=.15, if_khz=40., rf_detuning_mhz=0.,
        coupling_detuning_mhz=0., temperature_model="Linked"))
    positive = ("sensitivity_nv_cm_sqrt_hz", "psn_nv_cm_sqrt_hz", "transmission",
                "current_response_a_per_v_m", "detected_power_w")
    check("positive_finite_outputs", all(math.isfinite(row[key]) and row[key] > 0
          for row in rows for key in positive))
    check("bounded_controls_and_detuning", all(2 <= row["probe_power_uw"] <= 10
          and 2 <= row["lo_rabi_mhz"] <= 6 and abs(row["probe_detuning_mhz"]) <= 6 + 1e-9
          for row in rows))
    check("fixed_physical_detuning_grid", config["detuning_step_mhz"] == .01
          and all(abs(row["probe_detuning_mhz"] / .01 - round(row["probe_detuning_mhz"] / .01))
                  < 1e-9 for row in rows))
    check("hidden_controls_on_refined_lattice_or_protected_gold_lo", all(
        abs(row["probe_power_uw"] / .5 - round(row["probe_power_uw"] / .5)) < 1e-9
        and (abs(row["lo_rabi_mhz"] - 3.7) < 1e-9
             or abs((row["lo_rabi_mhz"] - 2) / .125 - round((row["lo_rabi_mhz"] - 2) / .125)) < 1e-9)
        for row in rows))
    check("passive_transmission", all(row["transmission"] <= 1 + 1e-12 for row in rows))
    check("psn_equals_total_under_zero_technical_noise", all(math.isclose(
        row["psn_nv_cm_sqrt_hz"], row["sensitivity_nv_cm_sqrt_hz"], rel_tol=1e-12) for row in rows))
    params = read_json(PRIOR / "resolved_parameters.json")["baseline_parameters"]
    responsivity = params["detector_quantum_efficiency"] * constants.ELEMENTARY_CHARGE / (
        2 * np.pi * constants.HBAR * species.RB85.line("D2")[1])
    noise_errors = [abs(math.sqrt(4 * constants.ELEMENTARY_CHARGE * responsivity * row["detected_power_w"])
        / row["current_response_a_per_v_m"] * 1e7 / row["psn_nv_cm_sqrt_hz"] - 1) for row in rows]
    check("independent_balanced_psn_and_field_unit_reconstruction", max(noise_errors) < 1e-11,
          {"max_relative_error": max(noise_errors)})
    groups = {}
    for row in rows:
        groups.setdefault(cell_key(row), []).append(row)
    all_controls = {}
    for row in rows:
        all_controls.setdefault(control_key(row), []).append(row["opacity"])
    check("complete_unique_opacity_scores_per_control", all(sorted(scales) == [.0005, .001, .0015]
          for scales in all_controls.values()))
    check("reported_state_and_score_counts", summary["control_states"] == len(all_controls)
          and summary["saved_score_rows"] == len(rows) == 3 * len(all_controls)
          and summary["map_coordinates_per_mode"] == 189)
    errors, gold_ratios = [], []
    for winner in selected:
        group = groups[cell_key(winner)]
        minimum = min(item["sensitivity_nv_cm_sqrt_hz"] for item in group)
        errors.append(abs(winner["sensitivity_nv_cm_sqrt_hz"] / minimum - 1))
        gold = [item for item in group if item["probe_power_uw"] == 6 and item["lo_rabi_mhz"] == 3.7]
        if gold:
            gold_ratios.append(winner["sensitivity_nv_cm_sqrt_hz"] / min(
                item["sensitivity_nv_cm_sqrt_hz"] for item in gold))
    check("each_map_cell_is_saved_control_minimum", max(errors) < 1e-12)
    best_errors = []
    for mode in ("off", "on"):
        saved = [row for row in summary["best_nominal"] if row["doppler"] == mode and row["opacity"] == .001]
        expected = min(row["sensitivity_nv_cm_sqrt_hz"] for row in selected
                       if row["doppler"] == mode and row["opacity"] == .001)
        best_errors.append(len(saved) == 1 and math.isclose(saved[0]["sensitivity_nv_cm_sqrt_hz"],
                                                           expected, rel_tol=1e-12))
    check("reported_nominal_best_matches_complete_map", all(best_errors))
    check("optimized_no_worse_than_fixed_gold_controls", len(gold_ratios) == len(expected_cells)
          and max(gold_ratios) <= 1 + 1e-12)
    prior = read_json(PRIOR / "results.json")
    old_controls = set(itertools.product(["off", "on"], prior["config"]["temperature_c"],
        prior["config"]["coupling_power_mw"], prior["config"]["probe_power_uw"], prior["config"]["lo_rabi_mhz"]))
    check("original_25_control_grid_reevaluated_at_original_map_cells",
          old_controls.issubset(set(map(control_key, rows))), {"expected_old_control_states": len(old_controls)})
    comparison = []
    for old in prior["optima"]:
        new = next(item for item in selected if cell_key(item) == cell_key(old))
        comparison.append({"cell": cell_key(old), "old": old["sensitivity_nv_cm_sqrt_hz"],
            "new": new["sensitivity_nv_cm_sqrt_hz"], "relative_change": new["sensitivity_nv_cm_sqrt_hz"]
                / old["sensitivity_nv_cm_sqrt_hz"] - 1})
    check("old_map_values_retained_or_improved_with_fixed_detuning_grid", max(
        item["relative_change"] for item in comparison) <= .001,
        {"maximum_relative_worsening": max(item["relative_change"] for item in comparison)})
    branch_errors = []
    for branch in summary["branchgaps"]:
        group = groups[cell_key(branch)]
        for prefix in ("central", "sideband"):
            quantity = prefix + "_sensitivity_nv_cm_sqrt_hz"
            expected = min(row[quantity] for row in group)
            branch_errors.append(abs(branch[quantity] / expected - 1))
    check("global_saved_central_and_sideband_minima_match", len(summary["branchgaps"])
          == len(expected_cells) and max(branch_errors) < 1e-12)
    epsilon_errors, epsilon_keys = [], []
    for band in summary["epsilon_bands"]:
        group = groups[cell_key(band)]
        best = min(item["sensitivity_nv_cm_sqrt_hz"] for item in group)
        accepted = [item for item in group if item["sensitivity_nv_cm_sqrt_hz"] <= best * (1 + band["epsilon"])]
        epsilon_keys.append((*cell_key(band), band["epsilon"]))
        epsilon_errors.append(band["saved_control_count"] == len(accepted)
            and band["probe_power_min_uw"] == min(item["probe_power_uw"] for item in accepted)
            and band["probe_power_max_uw"] == max(item["probe_power_uw"] for item in accepted)
            and band["lo_rabi_min_mhz"] == min(item["lo_rabi_mhz"] for item in accepted)
            and band["lo_rabi_max_mhz"] == max(item["lo_rabi_mhz"] for item in accepted)
            and band["represented_branches"] == "|".join(sorted(set(item["winning_branch"] for item in accepted))))
    check("epsilon_bands_match_saved_sample_membership", all(epsilon_errors)
          and len(epsilon_keys) == 2 * len(expected_cells)
          and set(epsilon_keys) == {(*key, epsilon) for key in expected_cells for epsilon in (.01, .03)})
    check("adaptive_control_stencils_finished_without_iteration_cap", not summary["failed_convergence"],
          {"failed_paths": len(summary["failed_convergence"]), "max_iterations": config["max_iterations"]})
    check("finite_if_mirror_shortcut_disabled", all(str(row["mirror_symmetry_used"]).lower()
          == "false" for row in rows))
    supplement = None
    if (ROOT / "off_branch_refinement.json").exists():
        supplement = read_json(ROOT / "off_branch_refinement.json")
        source_okay = all(sha == digest(ROOT / name)
                          for name, sha in supplement["source_sha256"].items())
        artifacts_okay = all(sha == digest(ROOT / name)
                             for name, sha in supplement["artifact_sha256"].items())
        from refine_off_branches import complete_off_snapshot
        _, current_off_snapshot = complete_off_snapshot()
        snapshot_okay = (supplement["off_snapshot"]["canonical_sha256"]
                         == current_off_snapshot["canonical_sha256"])
        check("selected_off_finer_control_branch_supplement_integrity", source_okay
              and artifacts_okay and snapshot_okay
              and supplement["off_snapshot"]["complete_off_cells"] == 189)
        check("selected_off_finer_control_branches_stable_within_declared_tolerance",
              supplement["status"] == "PASS" and supplement["maximum_absolute_relative_branch_change"] <= .005,
              {"maximum_absolute_relative_branch_change": supplement["maximum_absolute_relative_branch_change"],
               "minimum_absolute_relative_branch_margin": supplement["minimum_absolute_relative_branch_margin"],
               "branch_selection_changed_cases": supplement["branch_selection_changed_cases"]})
    local, local_rows, quad = [], [], []
    for row in selected_cells(summary, args.local_checks):
        result, records = local_rectangle_reference(row)
        local.append(result)
        local_rows.extend(records)
        print(json.dumps({"local_reference": result}, ensure_ascii=False), flush=True)
    if local:
        check("selected_local_rectangles_and_finer_detuning_agree", all(
            abs(item["relative_improvement"]) <= .005 for item in local), local)
    thermal = [row for row in selected_cells(summary, 2) if row["doppler"] == "on"]
    for row in thermal[:args.quadrature_checks]:
        result = thermal_quadrature_reference(row)
        quad.append(result)
        print(json.dumps({"thermal_quadrature_reference": result}, ensure_ascii=False), flush=True)
    if quad:
        check("independent_velocity_sum_reproduces_saved_fixed_control_sensitivity", all(
            abs(item["base_grid_relative_sensitivity_difference"]) <= 1e-8 for item in quad), quad)
        check("selected_thermal_finite_if_velocity_quadrature_converges", all(
            max(variant["relative_complex_error_vs_dv4_8sigma"] for variant in item["variants"])
            <= .005 for item in quad), quad)
    warnings.extend(["Sampled adaptive hidden-control search; no global continuous optimum established.",
        "Opacity/Rabi anchors and zero technical noise are model assumptions, not calibrated experimental forecasts.",
        "Local reference tests and velocity convergence probes cover selected settings only."])
    if not local:
        warnings.append("Optional exhaustive local/finer-detuning checks were not requested in this receipt.")
    if not quad:
        warnings.append("Optional independent finite-IF velocity quadrature probes were not requested in this receipt.")
    export_csv(ROOT / "verification_local_controls.csv", local_rows)
    report = {"origin": "⟦AI⟧ Independent dense-map verification", "status": "PASS"
        if all(item["status"] == "PASS" for item in checks) else "FAIL", "checks": checks,
        "warnings": warnings, "old_map_comparison": comparison,
        "local_references": local, "thermal_quadrature_references": quad,
        "off_finer_branch_supplement": supplement,
        "scope": "Saved-grid consistency plus explicitly requested selected convergence probes; no experimental predictive validation.",
        "artifact_sha256": {name: digest(ROOT / name) for name in (
            "results.json", "control_points.csv", "optima.csv", "config.json", "source_manifest.json", "verify_dense.py")}}
    for name in ("off_branch_refinement.json", "off_branch_refinement.csv", "refine_off_branches.py",
                 "quadrature_preflight.json", "symmetry_check.json"):
        if (ROOT / name).exists():
            report["artifact_sha256"][name] = digest(ROOT / name)
    (ROOT / "validation.json").write_text(json.dumps(report, ensure_ascii=False, indent=2,
        allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({"status": report["status"], "checks": len(checks),
        "failed": [item["id"] for item in checks if item["status"] == "FAIL"]}, ensure_ascii=False), flush=True)
    return 0 if report["status"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
