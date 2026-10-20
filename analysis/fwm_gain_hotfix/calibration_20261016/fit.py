"""Offline measured-reference calibration; no production mutation or runtime fit.

The dated fit_spec.md is declared before fitting. Historical calibration files
are inputs to the method review only and are never rewritten by this runner.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
import sys
import time

import numpy as np
from numpy.polynomial import chebyshev
from scipy.optimize import linprog

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
from gabes import constants as c, core
from gabes.schemes import fwm

OUT = ROOT/"docs/grand_challenge/gain_calibration_2026_10_16"
MODEL_ID = "measured_pump_off_ratio_chebyshev_20261016_v1"


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def save(name, value):
    (OUT/name).write_text(json.dumps(value, ensure_ascii=False, indent=2)+"\n", encoding="utf-8")


def measured_data():
    source = json.loads((OUT/"workbook_ooxml.json").read_text(encoding="utf-8"))
    if digest(Path(source["source"])) != source["sha256"]:
        raise ValueError("source workbook changed after OOXML extraction")
    sheet = next(s for s in source["sheets"] if s["name"] == "Sheet1")
    cells = {cell["cell"]: cell["value"] for row in sheet["rows"] for cell in row}
    data = {"eom_MHz": [cells[f"A{r}"] for r in range(5, 28)],
            "probe_voltage_mV": [cells[f"C{r}"] for r in range(5, 28)],
            "conjugate_voltage_mV": [cells[f"D{r}"] for r in range(5, 28)],
            "seed_reference_mV": cells["B1"],
            "seed_reference_power_uW": cells["D1"],
            "pump_power_mW": cells["B2"]}
    for arm in ("probe", "conjugate"):
        data[arm+"_ratio"] = (np.asarray(data[arm+"_voltage_mV"])/data["seed_reference_mV"]).tolist()
    return data, source["sha256"]


def atomic_baselines(data):
    paths = sorted((ROOT/"gabes").rglob("*.py"))
    before = {p.relative_to(ROOT).as_posix(): digest(p) for p in paths}
    power_proxy = data["seed_reference_power_uW"]*1e-6/(.8046*.9709)
    kwargs = dict(T=386.15, P_pump=data["pump_power_mW"]*1e-3,
        P_probe=power_proxy, L=.0125, w_pump=530e-6, w_probe=330e-6,
        line_strength=.74, qe=1., loss_frac=0., detection_efficiency=1.,
        branch=-1, phase_detail=fwm.PHASE_ULTRA, pump_probe_angle_deg=.32,
        model_fidelity=fwm.FIDELITY_BALANCED, gain_closure_enabled=False,
        legacy_gain_calibration_fixture=False, legacy_mean_field_fixture=False,
        floquet_order=3, response_method=fwm.RESPONSE_POLE,
        velocity_step=1., velocity_cutoff=4., fine_points=0)
    start = time.perf_counter()
    with core.blas_single_thread():
        exact = fwm.compute_spectrum(1., coarse_points=23,
            scan_min=1.-3.044, scan_max=1.-3.022, **kwargs)
        center = fwm.branch_center_GHz(1., -1)
        focused = fwm.compute_spectrum(1., coarse_points=401,
            scan_min=center-.060, scan_max=center+.060, **kwargs)
    after = {p.relative_to(ROOT).as_posix(): digest(p) for p in paths}
    changed = [name for name in before if before[name] != after[name]]
    if changed:
        raise RuntimeError("atomic baseline source changed during execution: "+repr(changed))
    expected_axis = 1.-np.asarray(data["eom_MHz"])/1000.
    np.testing.assert_allclose(exact["probe_axis_GHz"], expected_axis, rtol=0., atol=2e-15)
    assert not exact["gain_closure"]["applied"]
    assert not exact["gain_closure"]["output_calibration"]["applied"]
    payload = {
        "source_sha256_before": before, "source_sha256_after": after,
        "source_stable": True, "runtime_seconds": time.perf_counter()-start,
        "assumed_atomic_seed_power_W": power_proxy,
        "seed_proxy_origin": "external pre-window Poff/(fullcell*.aftercell), treated as vapor-input only for conditional comparison",
        "atomic_seed_independently_measured": False,
        "full_cell_transmission_used_as_A0": False,
        "conditions": {key: value for key, value in kwargs.items()
                       if isinstance(value, (float, int, bool, str))},
        "exact": {key: np.asarray(exact[key]).tolist() for key in
                  ("probe_axis_GHz", "G_s", "G_c")},
        "focused": {key: np.asarray(focused[key]).tolist() for key in
                    ("probe_axis_GHz", "G_s", "G_c")},
        "floquet_convergence": exact["floquet_convergence"],
        "scope": "conditional uncorrected atomic comparison; not measured-reference conversion or new independent-input validation"}
    save("uncorrected_baseline.json", payload)
    return payload


def metrics(prediction, target):
    prediction, target = np.asarray(prediction), np.asarray(target)
    relative = prediction/target-1
    return {"max_relative_error": float(np.max(abs(relative))),
            "rms_relative_error": float(np.sqrt(np.mean(relative**2))),
            "rmse_ratio": float(np.sqrt(np.mean((prediction-target)**2)))}


def monotone_order_bound(x, y, eom):
    best = {"lower_bound_max_relative_error": 0.0, "pair": None}
    count = 0
    for i in range(len(x)):
        for j in range(len(x)):
            if x[i] <= x[j] and y[i] > y[j]:
                count += 1
                bound = (y[i]-y[j])/(y[i]+y[j])
                if bound > best["lower_bound_max_relative_error"]:
                    best = {"lower_bound_max_relative_error": float(bound),
                            "pair": {"indices": [i, j], "eom_MHz": [float(eom[i]), float(eom[j])],
                                     "model_gains": [float(x[i]), float(x[j])],
                                     "target_ratios": [float(y[i]), float(y[j])]}}
    return best | {"order_reversal_pair_count": count,
                   "scope": "conditional on a constant positive source-to-measured plane conversion; unknown detuning-dependent A0 or arm throughput can change ordering",
                   "formula": "max over x_i<=x_j,y_i>y_j of (y_i-y_j)/(y_i+y_j)"}


def old_hinge_fit(x, y, knot, conjugate=False):
    basis = np.column_stack((x, np.ones(len(x)), np.maximum(x-knot, 0.)))
    a_ub = np.vstack((np.column_stack((basis, -y)),
                      np.column_stack((-basis, -y)), [-1., 0., -1., 0.]))
    b_ub = np.r_[y, -y, 0.]
    solved = linprog([0., 0., 0., 1.], A_ub=a_ub, b_ub=b_ub,
        bounds=[(0., None), ((0., 0.) if conjugate else (None, None)),
                (None, None), (0., None)], method="highs")
    if not solved.success:
        raise RuntimeError(solved.message)
    return {"fixed_knot": knot, "coefficients_a_b_c": solved.x[:3].tolist(),
            "minimax_relative_error": float(solved.x[3]),
            "prediction": (basis@solved.x[:3]).tolist(),
            "metrics": metrics(basis@solved.x[:3], y),
            "scope": "old monotone gain-only hinge method re-evaluated; no old coefficient reused"}


def log_minimax_coefficients(x, target, degree):
    basis = chebyshev.chebvander(x, degree)
    log_y = np.log(target)
    solved = linprog(np.r_[np.zeros(degree+1), 1.],
        A_ub=np.vstack((np.column_stack((basis, -np.ones(len(x)))),
                        np.column_stack((-basis, -np.ones(len(x)))))),
        b_ub=np.r_[log_y, -log_y],
        bounds=[(None, None)]*(degree+1)+[(0., None)], method="highs")
    if not solved.success:
        raise RuntimeError(solved.message)
    return solved.x[:-1], float(solved.x[-1])


def continuous_peak(coefficients):
    polynomial = chebyshev.Chebyshev(coefficients)
    roots = polynomial.deriv().roots()
    interior = [float(root.real) for root in roots
                if abs(root.imag) < 1e-10 and -1 < root.real < 1]
    points = np.array([-1., 1.]+interior)
    values = np.exp(polynomial(points))
    index = int(np.argmax(values))
    return {"eom_MHz": float(3033.+11.*points[index]), "ratio": float(values[index]),
            "method": "all real Chebyshev derivative roots and domain endpoints"}


def main():
    if not (OUT/"fit_spec.md").is_file():
        raise RuntimeError("declare fit_spec before running calibration")
    data, source_hash = measured_data()
    baseline = atomic_baselines(data)
    eom = np.asarray(data["eom_MHz"])
    coordinate = (eom-3033.)/11.
    train, validation = np.arange(0, 23, 2), np.arange(1, 23, 2)
    candidates, old_methods = [], {}
    for arm, key, knot in (("probe", "G_s", 7.), ("conjugate", "G_c", 6.1)):
        x, y = np.asarray(baseline["exact"][key]), np.asarray(data[arm+"_ratio"])
        old_methods[arm] = {"order_bound": monotone_order_bound(x, y, eom),
                            "old_hinge_method": old_hinge_fit(x, y, knot, arm == "conjugate")}
    for degree in [3, 5, 7]:
        arms = {}
        for arm in ("probe", "conjugate"):
            y = np.asarray(data[arm+"_ratio"])
            coefficients, error = log_minimax_coefficients(coordinate[train], y[train], degree)
            prediction = np.exp(chebyshev.chebval(coordinate, coefficients))
            arms[arm] = {"coefficients": coefficients.tolist(), "train_log_minimax": error,
                         "train_metrics": metrics(prediction[train], y[train]),
                         "validation_metrics": metrics(prediction[validation], y[validation]),
                         "prediction_all_rows": prediction.tolist()}
        score = max(value["validation_metrics"]["max_relative_error"] for value in arms.values())
        candidates.append({"degree": degree, "arms": arms,
                           "joint_validation_max_relative_error": score,
                           "passes_development_threshold": score <= .20})
    eligible = [case for case in candidates if case["passes_development_threshold"]]
    selected = (eligible[0] if eligible else min(candidates, key=lambda case: (
        case["joint_validation_max_relative_error"], case["degree"])))
    result_arms = {}
    for arm in ("probe", "conjugate"):
        y = np.asarray(data[arm+"_ratio"])
        coefficients, error = log_minimax_coefficients(coordinate, y, selected["degree"])
        prediction = np.exp(chebyshev.chebval(coordinate, coefficients))
        measured_peak = np.flatnonzero(y == np.max(y))
        result_arms[arm] = {"coefficients": coefficients.tolist(),
            "all_data_log_minimax": error, "prediction": prediction.tolist(),
            "full_data_metrics": metrics(prediction, y),
            "continuous_peak": continuous_peak(coefficients),
            "measured_peak": {"eom_MHz": eom[measured_peak].tolist(),
                              "ratio": float(np.max(y))}}
    result = {
        "origin": "AI offline measured-ratio development fit under direct user authorization",
        "marker": "⟦AI⟧", "date": "2026-10-16", "model_id": MODEL_ID,
        "source_workbook_sha256": source_hash, "fit_spec_sha256": digest(OUT/"fit_spec.md"),
        "domain_eom_MHz": [3022., 3044.],
        "coordinate": {"center_MHz": 3033., "scale_MHz": 11.},
        "coefficient_convention": "ascending Chebyshev coefficients of natural log ratio",
        "selected_degree": selected["degree"], "arms": result_arms, "data": data,
        "development": {"train_indices": train.tolist(), "validation_indices": validation.tolist(),
            "candidate_degrees": [3, 5, 7], "candidates": candidates,
            "acceptance_threshold_max_relative": .20,
            "acceptance_passed": bool(eligible),
            "selection_rule": "lowest validation-passing degree; otherwise minimum joint validation error, lower degree tie-break",
            "validation_type": "alternating-row development validation; all data already inspected; not untouched holdout"},
        "old_gain_only_method": old_methods,
        "context": {"temperature_C": 113., "opd_GHz": 1., "opd_status": "approximately 1 GHz; exact value not independently calibrated",
            "pump_power_W": data["pump_power_mW"]*1e-3,
            "cell_length_m": .0125, "pump_waist_m": 530e-6, "probe_waist_m": 330e-6,
            "crossing_angle_deg": .32, "seed_reference_mV": 35.,
            "seed_reference_power_W": data["seed_reference_power_uW"]*1e-6,
            "after_cell_transmission": .9709,
            "full_cell_transmission_including_windows": .8046,
            "atomic_pump_off_transmission": None,
            "assumed_atomic_input_proxy_W": baseline["assumed_atomic_seed_power_W"],
            "atomic_input_independently_measured": False},
        "scope": "empirical measured pump-on/pump-off voltage ratio at the supplied apparatus condition; input source unknowns retained",
        "out_of_domain_policy": "unavailable; no extrapolation",
        "atomic_transfer_modified": False, "atomic_gain_modified": False,
        "noise_NF_modified": False, "grand_challenge_no_fit_path_modified": False,
        "statistical_confidence_interval": None, "independent_validation": False,
        "runner_sha256": digest(Path(__file__)),
    }
    save("fit.json", result)
    print(json.dumps({"model_id": MODEL_ID, "selected_degree": result["selected_degree"],
        "acceptance_passed": result["development"]["acceptance_passed"],
        "candidates": [{"degree": case["degree"], "joint_validation_max_relative_error":
                        case["joint_validation_max_relative_error"]} for case in candidates],
        "arms": result_arms, "old_order_bounds": {
            arm: old_methods[arm]["order_bound"] for arm in old_methods}}, indent=2))


if __name__ == "__main__":
    main()
