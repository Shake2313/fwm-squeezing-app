"""Predeclared bounded basis pivot; reuses data/baseline without atomic solves."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
import sys

import numpy as np
from scipy.interpolate import BSpline, PPoly
from scipy.optimize import linprog

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT/"docs/grand_challenge/gain_calibration_2026_10_16"
sys.path.insert(0, str(ROOT))
from analysis.fwm_gain_hotfix.calibration_20261016.fit import metrics


def minimax(basis, target):
    p, n = basis.shape[1], basis.shape[0]
    y = np.log(target)
    fitted = linprog(np.r_[np.zeros(p), 1.],
        A_ub=np.vstack((np.column_stack((basis, -np.ones(n))),
                        np.column_stack((-basis, -np.ones(n))))),
        b_ub=np.r_[y, -y], bounds=[(None, None)]*p+[(0., None)], method="highs")
    if not fitted.success:
        raise RuntimeError(fitted.message)
    return fitted.x[:-1], float(fitted.x[-1])


def peak(knots, coefficients):
    spline = BSpline(knots, coefficients, 3, extrapolate=False)
    stationary = PPoly.from_spline(spline).derivative().roots(extrapolate=False)
    points = np.array([3022., 3044.]+[float(v) for v in stationary
                    if np.isfinite(v) and 3022. < v < 3044.])
    ratios = np.exp(spline(points))
    i = int(np.argmax(ratios))
    return {"eom_MHz": float(points[i]), "ratio": float(ratios[i]),
            "method": "all real piecewise polynomial derivative roots plus domain endpoints"}


def main():
    spec = (OUT/"fit_spec.md").read_text(encoding="utf-8")
    if "bounded basis pivot" not in spec:
        raise RuntimeError("declare pivot before fitting")
    first_path = OUT/"fit.json"
    first_hash = hashlib.sha256(first_path.read_bytes()).hexdigest()
    first = json.loads(first_path.read_text(encoding="utf-8"))
    data = first["data"]
    eom = np.asarray(data["eom_MHz"])
    train = np.asarray(first["development"]["train_indices"])
    validation = np.asarray(first["development"]["validation_indices"])
    candidates = []
    for interior_count in (2, 3, 4):
        interior = np.linspace(3022., 3044., interior_count+2)[1:-1]
        knots = np.r_[np.repeat(3022., 4), interior, np.repeat(3044., 4)]
        design = BSpline.design_matrix(eom, knots, 3, extrapolate=False).toarray()
        arms = {}
        for arm in ("probe", "conjugate"):
            y = np.asarray(data[arm+"_ratio"])
            fitted, error = minimax(design[train], y[train])
            prediction = np.exp(design@fitted)
            full, full_error = minimax(design, y)
            all_prediction = np.exp(design@full)
            measured_peak = np.flatnonzero(y == y.max())
            arms[arm] = {
                "development_train_coefficients": fitted.tolist(),
                "train_log_minimax": error,
                "train_metrics": metrics(prediction[train], y[train]),
                "validation_metrics": metrics(prediction[validation], y[validation]),
                "development_prediction_all_rows": prediction.tolist(),
                "all_data_refit": {"coefficients": full.tolist(),
                    "all_data_log_minimax": full_error,
                    "prediction": all_prediction.tolist(),
                    "full_data_metrics": metrics(all_prediction, y),
                    "continuous_peak": peak(knots, full),
                    "measured_peak": {"eom_MHz": eom[measured_peak].tolist(),
                                      "ratio": float(y.max())}}}
        maximum = max(v["validation_metrics"]["max_relative_error"] for v in arms.values())
        candidates.append({"interior_knots_count": interior_count,
            "coefficient_count_per_arm": int(design.shape[1]),
            "knot_vector_MHz": knots.tolist(), "arms": arms,
            "joint_validation_max_relative_error": maximum,
            "passes_development_threshold": maximum <= .20})
    passed = [case for case in candidates if case["passes_development_threshold"]]
    selected = passed[0] if passed else None
    result = {
        "origin": "AI predeclared bounded spline-basis pivot after failed global Chebyshev development validation",
        "marker": "⟦AI⟧", "date": "2026-10-16",
        "model_id": "measured_pump_off_ratio_bspline_20261016_v1",
        "source_workbook_sha256": first["source_workbook_sha256"],
        "historical_chebyshev_fit_sha256": first_hash,
        "fit_spec_sha256": hashlib.sha256((OUT/"fit_spec.md").read_bytes()).hexdigest(),
        "domain_eom_MHz": [3022., 3044.], "basis_degree": 3,
        "coefficient_convention": "natural-log ratio coefficients in scipy.interpolate.BSpline(knots_MHz,coefficients,3,extrapolate=False)",
        "selected_candidate": (None if selected is None else selected["interior_knots_count"]),
        "selected_interior_knots_count": (None if selected is None else selected["interior_knots_count"]),
        "knot_vector_MHz": (None if selected is None else selected["knot_vector_MHz"]),
        "arms": (None if selected is None else {arm: value["all_data_refit"]
                                                 for arm, value in selected["arms"].items()}),
        "data": data,
        "development": {"train_indices": train.tolist(), "validation_indices": validation.tolist(),
            "candidates": candidates, "candidate_interior_knots_counts": [2, 3, 4],
            "acceptance_threshold_max_relative": .20,
            "acceptance_passed": selected is not None,
            "selection_rule": "fewest uniformly spaced interior knots passing BOTH-arm development validation <=20%; no adoption if all fail",
            "validation_type": "same alternating-row development cohort; previously inspected data, not untouched holdout"},
        "context": first["context"],
        "scope": first["scope"], "out_of_domain_policy": "unavailable; no extrapolation",
        "new_atomic_solves": 0, "atomic_transfer_modified": False,
        "atomic_gain_modified": False, "noise_NF_modified": False,
        "grand_challenge_no_fit_path_modified": False,
        "statistical_confidence_interval": None, "independent_validation": False,
        "runner_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    }
    (OUT/"spline_fit.json").write_text(json.dumps(result, ensure_ascii=False, indent=2)+"\n", encoding="utf-8")
    assert hashlib.sha256(first_path.read_bytes()).hexdigest() == first_hash
    print(json.dumps({"selected_candidate": result["selected_candidate"],
        "acceptance_passed": result["development"]["acceptance_passed"],
        "candidate_validation": [{"interior_knots": v["interior_knots_count"],
            "coefficient_count": v["coefficient_count_per_arm"],
            "arms": {arm: {"train": a["train_metrics"], "validation": a["validation_metrics"],
                      "full_data": a["all_data_refit"]["full_data_metrics"],
                      "peak": a["all_data_refit"]["continuous_peak"]} for arm,a in v["arms"].items()}}
                    for v in candidates]}, indent=2))


if __name__ == "__main__":
    main()
