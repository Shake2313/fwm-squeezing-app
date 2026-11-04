"""Verify saved Gold evidence and scientific contracts without atomic solves.

Only standard-library arithmetic, CSV/JSON reads, AST inspection and hashes are
used. A numerical-integrity PASS does not establish agreement with experiment.
"""
from __future__ import annotations

import ast
import csv
import hashlib
import json
import math
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parents[1]
MODEL_SHA256 = "e967843a91ee4e39b14153997679e3f8130a3987da316d957b543e45edf6c15b"
POLICIES = {"maximum_static_slope", "bare_two_photon_center"}
CHECKS = []
READ_FILES = set()
HASH_REFERENCES = []


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_json(name):
    path = ROOT / name
    READ_FILES.add(path)
    return json.loads(path.read_text(encoding="utf-8-sig"))


def read_csv(name):
    path = ROOT / name
    READ_FILES.add(path)
    with path.open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def check(name, condition, **details):
    CHECKS.append(dict(name=name, status="PASS" if condition else "FAIL", **details))


def near(first, second, tolerance=2e-11):
    return math.isclose(float(first), float(second), rel_tol=tolerance, abs_tol=2e-12)


def rms(values):
    return math.sqrt(sum(value * value for value in values) / len(values))


def verify_hash_map(name, mapping, base):
    mismatches = []
    for relative, expected in mapping.items():
        path = base / relative
        actual = sha(path) if path.is_file() else None
        HASH_REFERENCES.append(dict(manifest=name, path=str(path.relative_to(REPO)).replace("\\", "/"),
                                    expected_sha256=expected, actual_sha256=actual))
        if actual != expected:
            mismatches.append(relative)
        if path.is_file():
            READ_FILES.add(path)
    check(name, not mismatches, references=len(mapping), mismatches=mismatches)


def coordinates(rows):
    return {(float(row["temperature_c"]), float(row["coupling_power_mw"])) for row in rows}


def value_index(rows):
    return {(row["detuning_policy"], float(row["temperature_c"]), float(row["coupling_power_mw"])):
            float(row["sensitivity_nv_cm_sqrt_hz"]) for row in rows}


def interpolate(value, nodes, ordinates):
    if not nodes[0] <= value <= nodes[-1]:
        raise ValueError("No extrapolation is certified")
    for first in range(len(nodes) - 1):
        if value <= nodes[first + 1]:
            fraction = (value - nodes[first]) / (nodes[first + 1] - nodes[first])
            return ordinates[first] + fraction * (ordinates[first + 1] - ordinates[first])
    return ordinates[-1]


def linked_rb_density(temp_c, boltzmann):
    """Independent closed-form CRC density, used only for saved-data prediction."""
    kelvin = temp_c + 273.15
    exponent = 9.863 - 4215.0 / kelvin if temp_c < 39.30 else 9.318 - 4040.0 / kelvin
    return 10.0 ** exponent / (boltzmann * kelvin)


def literal_assignment(tree, target):
    for node in ast.walk(tree):
        if isinstance(node, ast.Assign) and any(isinstance(item, ast.Name) and item.id == target
                                                for item in node.targets):
            return ast.literal_eval(node.value)
    raise KeyError(target)


def inspect_numeric_anchors():
    paths = [REPO / name for name in ("gabes/constants.py", "gabes/species.py", "gabes/schemes/rydberg.py")]
    trees = [ast.parse(path.read_text(encoding="utf-8-sig")) for path in paths]
    READ_FILES.update(paths)
    charge = literal_assignment(trees[0], "ELEMENTARY_CHARGE")
    hbar = literal_assignment(trees[0], "HBAR")
    boltzmann = literal_assignment(trees[0], "KB")
    rb85 = next(node.value for node in ast.walk(trees[1]) if isinstance(node, ast.Assign)
                and any(isinstance(item, ast.Name) and item.id == "RB85" for item in node.targets))
    frequency = ast.literal_eval(next(item.value for item in rb85.keywords if item.arg == "nu_D2"))
    detector = {}
    for key in ("detector_quantum_efficiency", "detector_path_efficiency", "detector_reference_power_ratio",
                "detector_electronic_noise_pa_sqrt_hz", "detector_rin_per_sqrt_hz"):
        candidates = [ast.literal_eval(node.value) for node in ast.walk(trees[2])
                      if isinstance(node, ast.keyword) and node.arg == key and isinstance(node.value, ast.Constant)]
        detector[key] = candidates[0]
    expected = dict(detector_quantum_efficiency=.95, detector_path_efficiency=.5,
                    detector_reference_power_ratio=1., detector_electronic_noise_pa_sqrt_hz=0.,
                    detector_rin_per_sqrt_hz=0.)
    check("PSN detector assumptions agree with frozen source", detector == expected, assumptions=detector,
          optical_frequency_hz=frequency, constants_read_by="AST only; no model imports")
    return charge, hbar, frequency, boltzmann


def verify_prediction_reports(name, reports, training, direct_rows, tnodes, pnodes, boltzmann, refined):
    truth = value_index(direct_rows)
    results = []
    for report in reports:
        policy = report["detuning_policy"]
        normalized = True if refined else report["density_normalized"]
        thermal = [math.log(training[policy, t, 30.] *
                    (linked_rb_density(t, boltzmann) if normalized else 1.)) for t in tnodes]
        power = [math.log(training[policy, 20., p] / training[policy, 20., 30.]) for p in pnodes]
        errors, consistent = [], True
        for prediction in report["predictions"]:
            temperature, coupling = prediction["temperature_c"], prediction["coupling_power_mw"]
            estimated = math.exp(interpolate(temperature, tnodes, thermal) + interpolate(coupling, pnodes, power))
            if normalized:
                estimated /= linked_rb_density(temperature, boltzmann)
            actual = truth[policy, temperature, coupling]
            stored_estimate = prediction["predicted"] if refined else prediction["predicted_sensitivity"]
            stored_actual = prediction["direct"] if refined else prediction["direct_sensitivity"]
            error = estimated / actual - 1.
            consistent &= near(estimated, stored_estimate) and near(actual, stored_actual) and near(error, prediction["relative_error"])
            errors.append(error)
        maximum, mean_square = max(map(abs, errors)), rms(errors)
        eligible = report["eligible_at_one_percent"] if refined else report["eligible_at_declared_criterion"]
        consistent &= near(maximum, report["maximum_relative_error"]) and near(mean_square, report["rms_relative_error"])
        consistent &= eligible == (maximum <= .01)
        consistent &= len(errors) == 8 and report["training_cross_coordinates"] == (7 if refined else 5)
        results.append(dict(detuning_policy=policy, density_normalized=normalized,
                            maximum_relative_error=maximum, independent_points=len(errors), consistent=bool(consistent)))
    check(name, all(item["consistent"] for item in results), reports=results,
          calculation="Reconstruct log-cross forecast from training CSV only; exact closed-form N(T), no atomic solves")


def main():
    benchmark = read_json("benchmark_results.json")
    current = read_json("current_resolution_audit.json")
    refined = read_json("refined_resolution_audit.json")
    registry = read_json("current_variable_resolution_registry.json")
    figures = read_json("figure_manifest.json")
    source = read_json("sources/gold_observations.json")
    extraction = read_json("sources/source_validation.json")
    historical = read_json("reduction/variable_resolution_registry.json")
    historical_receipt = read_json("reduction/validation.json")
    v1_plan, v2_plan = read_json("current_guard_plan.json"), read_json("refined_guard_plan.json")

    for name, document in (("benchmark", benchmark), ("current audit", current), ("refined audit", refined), ("current registry", registry)):
        verify_hash_map(name + " physical/source hashes", document["source_sha256"], REPO)
        if "artifact_sha256" in document:
            verify_hash_map(name + " artifact hashes", document["artifact_sha256"], ROOT)
    check("frozen CP production model fingerprint", sha(REPO / "gabes/schemes/rydberg.py") == MODEL_SHA256,
          sha256=MODEL_SHA256, dephasing_model="correlated_phase_lindblad")
    verify_hash_map("figure input hashes", figures["input_sha256"], ROOT)
    verify_hash_map("figure renderer hash", {"render_validation.py": figures["code_sha256"]}, ROOT)
    verify_hash_map("figure output hashes", figures["output_sha256"], ROOT)
    verify_hash_map("extracted source artifacts", source["artifact_sha256"], ROOT / "sources")
    verify_hash_map("primary source PDF", {source["source"]["path"]: source["source"]["sha256"]}, REPO)
    check("source extraction receipt and physical scope",
          extraction["extraction_status"] == "PASS" and all(extraction["checks"].values())
          and extraction["observations_sha256"] == sha(ROOT / "sources/gold_observations.json")
          and extraction["physical_protocol_match_status"] == "NOT_ESTABLISHED",
          reported_checks=extraction["checks"], visual_qa="Existing human/agent receipt only; images not reinspected here")
    points = read_csv("sources/fig2b_points.csv")
    psn = sorted(read_csv("sources/fig5b_psn_curve.csv"), key=lambda row: float(row["probe_power_uw"]))
    expected_power = [.5, 1., 2., 4., 6., 8., 10.]
    check("source vector-point coverage", len(points) == 28 and
          {(float(row["probe_power_uw"]), row["condition"], row["quantity"]) for row in points} ==
          {(p, c, q) for p in expected_power for c in ("compensated", "uncompensated") for q in ("eit_fwhm", "eit_peak")},
          points=len(points), experimental_uncertainty="No standard uncertainty inferred from digitization resolution")
    check("source PSN line preserves decreasing power dependence", len(psn) == 7 and all(
          float(first["sensitivity_nv_cm_sqrt_hz"]) > float(second["sensitivity_nv_cm_sqrt_hz"])
          for first, second in zip(psn, psn[1:])), points=len(psn), interpretation="Published theory curve; not measured total noise")
    axis_audit = read_json("sources/scan_axis_audit.json")
    check("Gold unknown scan identity retained", axis_audit["gold_main"]["scan_laser"] == "not reported"
          and axis_audit["gold_main"]["absolute_optical_detunings"] == "not reported"
          and all(not item.get("absolute_transmission_available", False) for item in source["protocols"]),
          prior_protocol="Different apparatus supports hypotheses only")

    spectra = read_csv("spectroscopy_comparison.csv")
    observed = {(float(row["probe_power_uw"]), row["condition"], row["quantity"]): float(row["value"])
                for row in points}
    check("spectroscopy references equal extracted Gold observations", all(
          near(row["source_width_mhz"], observed[float(row["probe_power_uw"]), row["condition"], "eit_fwhm"])
          and near(row["source_peak_arb"], observed[float(row["probe_power_uw"]), row["condition"], "eit_peak"])
          for row in spectra), comparisons=len(spectra), independent_source_points=len(points))
    train_p, heldout_p = set(benchmark["training_powers_uw"]), set(benchmark["heldout_powers_uw"])
    check("power split is disjoint and complete", not train_p & heldout_p and train_p | heldout_p == set(expected_power),
          training=sorted(train_p), heldout=sorted(heldout_p))
    check("spectroscopy coverage and passive outputs", len(spectra) == 56 and all(
          float(row["lo_rabi_mhz"]) == 0. and 0. <= float(row["model_peak_transmission"]) <= 1.
          and 0. <= float(row["model_contrast"]) <= 1. and float(row["eit_width_mhz"]) > 0.
          and (row["partition"] == "training") == (float(row["probe_power_uw"]) in train_p) for row in spectra), rows=len(spectra))
    gains = []
    fit_ok, metrics_ok = True, True
    for mode in ("off", "on"):
        for axis in ("Probe detuning", "Coupling detuning"):
            group = [row for row in spectra if row["doppler"] == mode and row["scan_axis"] == axis]
            training = [row for row in group if row["partition"] == "training"]
            gain = sum(float(row["model_relative_signal"]) * float(row["source_peak_arb"]) for row in training) / sum(
                float(row["model_relative_signal"]) ** 2 for row in training)
            gains.append(dict(doppler=mode, scan_axis=axis, training_only_shared_gain=gain))
            fit_ok &= gain > 0. and len(training) == 8 and len(group) == 14
            for row in group:
                prediction = gain * float(row["model_relative_signal"])
                fit_ok &= near(gain, row["signal_gain_fitted_on_training"]) and near(prediction, row["predicted_signal_arb"])
                fit_ok &= near(float(row["probe_power_uw"]) * float(row["model_contrast"]), row["model_relative_signal"])
                fit_ok &= near(float(row["eit_width_mhz"]) - float(row["source_width_mhz"]), row["width_error_mhz"])
                fit_ok &= near(prediction - float(row["source_peak_arb"]), row["signal_error_arb"])
            for partition in ("training", "heldout"):
                rows = [row for row in group if row["partition"] == partition]
                recorded = next(item for item in benchmark["comparisons"] if item["doppler"] == mode and
                                item["scan_axis"] == axis and item["partition"] == partition)
                metrics_ok &= recorded["count"] == len(rows) and near(recorded["shared_detector_gain"], gain)
                metrics_ok &= near(recorded["width_rmse_mhz"], rms([float(row["width_error_mhz"]) for row in rows]))
                metrics_ok &= near(recorded["signal_rmse_arb"], rms([float(row["signal_error_arb"]) for row in rows]))
    check("one gain per model/axis trained only on training powers", bool(fit_ok), groups=gains,
          sharing="Same gain for compensated/uncompensated conditions and heldout powers")
    check("saved comparison metrics independently reconstructed", bool(metrics_ok), comparisons=len(benchmark["comparisons"]))

    rf_names = ["fixed_lo_rf.csv", "current_model_sparse_pilot.csv", "current_guard_points.csv",
                "refined_axis_points.csv", "refined_guard_points.csv", "independent_detuning_pilot.csv"]
    rf_tables = {name: read_csv(name) for name in rf_names}
    expected_counts = dict(zip(rf_names, [56, 18, 16, 4, 16, 6]))
    check("RF evidence row counts", {name: len(rows) for name, rows in rf_tables.items()} == expected_counts,
          counts={name: len(rows) for name, rows in rf_tables.items()}, total_rows=sum(map(len, rf_tables.values())))
    rows = [row for table in rf_tables.values() for row in table]
    check("fixed LO and IF across RF evidence", all(float(row["lo_rabi_mhz"]) == 3.7 and float(row["if_khz"]) == 40.
          and float(row["lo_carrier_ghz"]) == 37.07104 and float(row["rf_detuning_mhz"]) == 0.
          and row["detuning_policy"] in POLICIES for row in rows), rows=len(rows),
          RF_zero="Analysis assumption; absolute atomic resonance not independently calibrated")
    check("passive finite RF readouts", all(0. < float(row["transmission"]) <= 1. and
          math.isfinite(float(row["sensitivity_nv_cm_sqrt_hz"])) and float(row["sensitivity_nv_cm_sqrt_hz"]) > 0.
          and math.isfinite(float(row["responsivity_a_per_v_m"])) and abs(float(row["responsivity_a_per_v_m"])) > 0.
          and float(row["detected_power_w"]) > 0. for row in rows), rows=len(rows))
    detuning_ok = True
    for row in rows:
        probe, coupling, scan = (float(row[k]) for k in ("chosen_probe_detuning_mhz", "chosen_coupling_detuning_mhz", "chosen_scan_detuning_mhz"))
        if row["scan_axis"] == "Probe detuning":
            detuning_ok &= near(probe, scan) and near(coupling, row["fixed_coupling_detuning_mhz"])
        else:
            detuning_ok &= row["scan_axis"] == "Coupling detuning" and near(coupling, scan) and near(probe, row["fixed_probe_detuning_mhz"])
        if row["detuning_policy"] == "bare_two_photon_center":
            detuning_ok &= near(probe + coupling, 0.)
    check("saved physical optical coordinates agree with scan axis", bool(detuning_ok), rows=len(rows),
          centre_scope="Bare two-photon centre, not a fitted spectral maximum")
    charge, hbar, frequency, boltzmann = inspect_numeric_anchors()
    responsivity = .95 * charge / (2. * math.pi * hbar * frequency)
    psn_errors, power_errors = [], []
    for row in rows:
        power = float(row["detected_power_w"])
        independently_reconstructed = math.sqrt(4. * charge * responsivity * power) / abs(float(row["responsivity_a_per_v_m"])) * 1e7
        psn_errors.append(abs(independently_reconstructed / float(row["sensitivity_nv_cm_sqrt_hz"]) - 1.))
        expected_power = .5 * float(row["probe_power_uw"]) * 1e-6 * float(row["transmission"])
        power_errors.append(abs(power / expected_power - 1.))
    check("equal-arm photon shot noise reconstructed independently", max(psn_errors) < 2e-12 and max(power_errors) < 2e-12,
          rows=len(rows), maximum_relative_PSN_error=max(psn_errors), maximum_relative_power_error=max(power_errors),
          equation="S_E=sqrt(4 e R_PD P_signal)/abs(H_IE)*1e7 nV/cm/sqrt(Hz); P_signal=0.5 Pp T",
          scope="Equal DC arms, QE0.95, zero RIN/electronics; verifies saved readout arithmetic only")

    pilot = rf_tables["current_model_sparse_pilot.csv"]
    old_guards = rf_tables["current_guard_points.csv"]
    new_axis = rf_tables["refined_axis_points.csv"]
    new_guards = rf_tables["refined_guard_points.csv"]
    groups = [coordinates(table) for table in (pilot, old_guards, new_axis, new_guards)]
    check("pilot, V1 guards, added axis and V2 guards disjoint", all(not first & second
          for i, first in enumerate(groups) for second in groups[i + 1:]), counts=list(map(len, groups)),
          unique_atomic_coordinates=len(set.union(*groups)))
    check("27 unique atomic control coordinates honestly counted", len(set.union(*groups)) == 27 ==
          refined["actual_unique_atomic_control_coordinates_across_v1_and_v2"], training_and_guard_costs=[9, 8, 2, 8])
    check("saved guard plans agree with direct CSV coordinates", set(map(tuple, v1_plan["points"])) == groups[1]
          and set(map(tuple, v2_plan["new_training_coordinates"])) == groups[2]
          and set(map(tuple, v2_plan["fresh_guard_coordinates"])) == groups[3]
          and current["guard_plan"] == v1_plan and refined["plan"] == v2_plan,
          chronology="Plans describe predeclared coordinates; hashes cannot independently prove wall-clock chronology")
    check("guard tolerance is explicit AI numerical criterion", v1_plan["criterion_relative"] == .01
          and v2_plan["tolerance_relative"] == .01 and "AI" in v1_plan["criterion_scope"] and
          "AI" in v2_plan["tolerance_status"], tolerance_relative=.01,
          adoption="Not experimental uncertainty or a human-registered goal")
    training = value_index(pilot + new_axis)
    verify_prediction_reports("V1 forecasts reconstructed from five training cross points", current["reports"], training,
                              old_guards, [15., 20., 35.], [20., 30., 40.], boltzmann, False)
    verify_prediction_reports("V2 forecasts reconstructed from seven training cross points", refined["reports"], training,
                              new_guards, [15., 20., 35.], [20., 25., 30., 35., 40.], boltzmann, True)
    failed_v1 = next(report for report in current["reports"] if report["detuning_policy"] == "maximum_static_slope"
                     and report["density_normalized"])
    check("failed V1 retained separately from eligible fresh V2 guards", not failed_v1["eligible_at_declared_criterion"]
          and failed_v1["maximum_relative_error"] > .01 and all(report["eligible_at_one_percent"] for report in refined["reports"]),
          V1_maximum_relative_error=failed_v1["maximum_relative_error"], V2_maximum_relative_errors=[r["maximum_relative_error"] for r in refined["reports"]])

    check("current registry evidence links and scope", registry["evidence_manifest_sha256"] == sha(ROOT / "current_resolution_audit.json")
          and registry["adaptive_refinement_evidence_sha256"] == sha(ROOT / "refined_resolution_audit.json")
          and registry["selected_current_policy"]["reports"] == refined["reports"]
          and registry["previous_failed_guard_preserved"] == "current_resolution_audit.json"
          and registry["selected_current_policy"]["status"] == "ELIGIBLE_ON_TESTED_GUARDS"
          and len(registry["invalidate_when"]) >= 5, roles=len(registry["variables"]))
    registry_csv = read_csv("current_variable_resolution_registry.csv")
    registry_hash = sha(ROOT / "current_variable_resolution_registry.json")
    numeric_roles = {row["variable"] for row in registry_csv if row["maximum_relative_error_static"] or row["maximum_relative_error_center"]}
    report_errors = {row["detuning_policy"]: row["maximum_relative_error"] for row in refined["reports"]}
    csv_numbers_match = all(
        near(row["maximum_relative_error_static"], report_errors["maximum_static_slope"])
        and near(row["maximum_relative_error_center"], report_errors["bare_two_photon_center"])
        for row in registry_csv if row["variable"] in {"temperature_c", "coupling_power_mw"})
    registry_variables = {row["variable"]: row for row in registry["variables"]}
    csv_roles_match = all(row["role"] == registry_variables[row["variable"]]["role"]
                          and row["recommendation"] == registry_variables[row["variable"]]["recommendation"]
                          and (not row["training_nodes"] or json.loads(row["training_nodes"]) ==
                               registry_variables[row["variable"]].get("proposed_training_nodes", []))
                          for row in registry_csv)
    check("human-readable variable registry matches current JSON", len(registry_csv) == len(registry["variables"]) == 7
          and {row["variable"] for row in registry_csv} == {row["variable"] for row in registry["variables"]}
          and all(row["source_registry_sha256"] == registry_hash and row["origin"] == "⟦AI⟧"
                  and row["physical_validation"] == "NOT_ESTABLISHED" for row in registry_csv)
          and numeric_roles == {"temperature_c", "coupling_power_mw"} and csv_numbers_match and csv_roles_match, rows=len(registry_csv),
          numeric_error_roles=sorted(numeric_roles), source_registry_sha256=registry_hash)
    check("historical model eligibility is invalidated", historical["applicability_status"] == "INVALIDATED_FOR_CURRENT_SOURCE"
          and "gabes/schemes/rydberg.py" in historical["current_source_mismatches"]
          and historical["expected_critical_source_sha256"]["gabes/schemes/rydberg.py"] != MODEL_SHA256,
          historical_model=historical["model_tag"], current_applicability=historical["applicability_status"])
    old_hashes = historical_receipt["artifact_sha256"]
    READ_FILES.update(ROOT / "reduction" / name for name in old_hashes)
    old_mismatches = [name for name, digest in old_hashes.items() if sha(ROOT / "reduction" / name) != digest]
    check("historical helper test-code receipt remains applicable", historical_receipt["status"] == "PASS"
          and historical_receipt["tests"]["outcome"].startswith("9 passed") and
          set(old_mismatches) <= {"variable_resolution_registry.json", "variable_resolution_registry.csv"},
          historical_output_hashes_changed=old_mismatches,
          interpretation="Invalidation metadata regenerated after original receipt; current hashes recorded below, original receipt retained")
    check("unresolved experiment states are explicit", benchmark["absolute_transmission_validation"] == "NOT_IDENTIFIABLE_FROM_ARBITRARY_SIGNAL_UNITS"
          and benchmark["physical_contrast_validation"] == "NOT_IDENTIFIABLE_WITHOUT_GAIN_BASELINE"
          and benchmark["protocol_status"].startswith("PARTIAL_MATCH")
          and benchmark["optimization_status"].startswith("EXPERIMENTAL_OPTIMUM_NOT_VALIDATED")
          and registry["physical_model_validation"] == "NOT_ESTABLISHED"
          and refined["experimental_validation_status"] == "NOT_ESTABLISHED",
          absolute_transmission=benchmark["absolute_transmission_validation"], physical_contrast=benchmark["physical_contrast_validation"],
          protocol=benchmark["protocol_status"], physical_model_validation=registry["physical_model_validation"])
    check("no reported-sensitivity or temperature-optimum fit", benchmark["no_sensitivity_fit"] is True
          and benchmark["no_temperature_optimum_fit"] is True, no_sensitivity_fit=True, no_temperature_optimum_fit=True,
          scope="Explicit saved-analysis flags; detector arbitrary-unit gain is fitted separately")

    readme_path = ROOT / "README.md"
    if readme_path.is_file():
        READ_FILES.add(readme_path)
        readme = readme_path.read_text(encoding="utf-8-sig")
        broken = []
        for target in re.findall(r"\]\(([^)]+)\)", readme):
            if "://" in target or target.startswith("#"):
                continue
            target = target.split("#", 1)[0]
            if target != "validation.json" and not (ROOT / target).exists():
                broken.append(target)
        check("README local artifact links", not broken, missing=broken, validation_link="Created by this verifier")
        claims = ["0.588%", "0.304%", "27개", "108 passed", "9 passed", "36 checks PASS", "NOT_IDENTIFIABLE", MODEL_SHA256]
        check("README headline claims and scope match receipts", all(item in readme for item in claims)
              and "Gold 실험의 정량적 재현은 아직 성립하지 않는다" in readme,
              checked_claims=claims, scope="Headline values and linked-artifact presence; not independent line-by-line scientific review")
        numeric_claims = [f"{item['width_rmse_mhz']:.4f}" for item in benchmark["comparisons"] if item["partition"] == "heldout"]
        numeric_claims += [f"{item['signal_rmse_arb']:.5f}" for item in benchmark["comparisons"] if item["partition"] == "heldout"]
        numeric_claims += [f"{item['maximum_relative_error'] * 100.:.3f}%" for item in refined["reports"]]
        numeric_claims += [f"{float(row['sensitivity_nv_cm_sqrt_hz']):.1f}" for row in rf_tables["fixed_lo_rf.csv"]
                           if row["doppler"] == "on" and row["scan_axis"] == "Probe detuning" and float(row["probe_power_uw"]) == 6.]
        check("README current numeric tables match saved results", all(item in readme for item in numeric_claims),
              checked_values=numeric_claims, scope="Current heldout spectroscopy RMSE, refined guard errors and nominal6uW RF readout")

    reported_receipts = [
        dict(origin="⟦AI⟧ Parent-reported actual prior execution", command="python -X utf8 -B -m pytest -q tests/test_rydberg_model_physics.py tests/test_rydberg_eit.py tests/test_rydberg_electrometry.py tests/test_rydberg_opacity.py tests/test_rydberg_experiment.py tests/test_rydberg_experimental_csv.py",
             outcome="108 passed in 10.17s", warnings="2 existing datetime.utcnow DeprecationWarnings",
             subsequent_change="Test function/doc names only, as parent reported", freshly_executed_by_this_verifier=False),
        dict(origin="⟦AI⟧ Persisted historical helper-test receipt", receipt="reduction/validation.json",
             **historical_receipt["tests"], freshly_executed_by_this_verifier=False),
        dict(origin="⟦AI⟧ Parent-reported actual prior execution", command="python -X utf8 -B references/rydberg_experiments/verify_db.py",
             outcome="36 checks PASS", reported_table_counts=dict(papers=26, experiments=39, numbers=171, findings=42),
             freshly_executed_by_this_verifier=False),
    ]
    failures = [item["name"] for item in CHECKS if item["status"] == "FAIL"]
    artifact_hashes = {str(path.relative_to(ROOT)).replace("\\", "/"): sha(path) for path in sorted(READ_FILES)
                       if path.is_relative_to(ROOT)}
    artifact_hashes["verify_analysis.py"] = sha(Path(__file__))
    result = dict(origin="⟦AI⟧ Saved-data integrity and scientific-contract verification", status="FAIL" if failures else "PASS",
        scope="No new atomic solves, fitting, pytest or DB-verifier execution. Standard-library hashes, AST and independent arithmetic only.",
        experimental_validation_status="NOT_ESTABLISHED", unresolved_experiment_is_numerical_failure=False,
        checks_passed=len(CHECKS) - len(failures), checks_total=len(CHECKS), failed_checks=failures, checks=CHECKS,
        hash_reference_count=len(HASH_REFERENCES), unique_hashed_reference_files=len({item["path"] for item in HASH_REFERENCES}),
        hash_references=HASH_REFERENCES, model_sha256=MODEL_SHA256, reported_execution_receipts=reported_receipts,
        intentional_omissions=["No physics campaign or full-suite/GC rerun", "No new source-figure visual QA",
                               "No claim of absolute transmission/physical-contrast identification or Gold optimum validation"],
        artifact_sha256=artifact_hashes)
    (ROOT / "validation.json").write_text(json.dumps(result, ensure_ascii=False, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({key: result[key] for key in ("status", "checks_passed", "checks_total", "failed_checks", "hash_reference_count", "unique_hashed_reference_files", "experimental_validation_status")}, indent=2))
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
