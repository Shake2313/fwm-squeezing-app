"""Deterministic sampled measured-ratio reference, not a predictive fit."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
import sys

import numpy as np
from scipy.interpolate import PchipInterpolator

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT/"docs/grand_challenge/gain_calibration_2026_10_16"
sys.path.insert(0, str(ROOT))
from analysis.fwm_gain_hotfix.calibration_20261016.fit import metrics

MODEL_ID = "measured_pump_off_ratio_log_pchip_20261016_v1"


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def interpolate(x, y):
    order = np.argsort(x)
    return PchipInterpolator(x[order], np.log(y[order]), extrapolate=False)


def peak_intervals(x, y):
    maximum = float(max(y))
    matches = np.flatnonzero(y == maximum)
    groups = []
    for index in matches:
        if groups and index == groups[-1][-1]+1:
            groups[-1].append(int(index))
        else:
            groups.append([int(index)])
    return {"intervals_eom_MHz": [[float(x[g[0]]), float(x[g[-1]])] for g in groups],
            "ratio": maximum, "method": "shape-preserving interval extrema; retain raw maximum plateau"}


def checks(curve, x, y):
    node = np.exp(curve(x))
    maximum_overshoot, derivative_violation = 0., 0.
    global_extrema = []
    for index in range(len(x)-1):
        lo, hi = x[index:index+2]
        width = hi-lo
        coefficients = curve.c[:, index]
        derivative = np.array([3*coefficients[0], 2*coefficients[1], coefficients[2]])
        nonzero = np.trim_zeros(derivative, "f")
        roots = np.roots(nonzero) if len(nonzero) > 1 else []
        locations = [0., width]+[float(r.real) for r in roots
                                   if abs(r.imag) < 1e-12 and 0 < r.real < width]
        values = np.exp(curve(lo+np.array(locations)))
        expected_low, expected_high = sorted(y[index:index+2])
        maximum_overshoot = max(maximum_overshoot,
            float(max(expected_low-values.min(), values.max()-expected_high, 0.)))
        extrema = [0., width]
        if derivative[0] != 0:
            vertex = -derivative[1]/(2*derivative[0])
            if 0 < vertex < width:
                extrema.append(vertex)
        slopes = np.polyval(derivative, extrema)
        direction = np.sign(y[index+1]-y[index])
        violation = (np.max(abs(slopes)) if direction == 0
                     else max(float(np.max(-direction*slopes)), 0.))
        derivative_violation = max(derivative_violation, float(violation))
        global_extrema.extend(values.tolist())
    outside = np.exp(curve(np.array([x[0]-1e-6, x[-1]+1e-6])))
    node_error = float(np.max(abs(node/y-1)))
    tolerance = 1e-12*max(1., float(max(y)))
    result = {"source_points_positive": bool(np.all(y > 0)),
        "source_RF_strictly_increasing": bool(np.all(np.diff(x) > 0)),
        "exact_node_max_relative_error": node_error,
        "interval_no_overshoot_max_absolute_ratio_error": maximum_overshoot,
        "interval_derivative_sign_max_violation_per_MHz": derivative_violation,
        "positive_at_all_analytic_interval_extrema": bool(min(global_extrema) > 0),
        "outside_support_returns_NaN": bool(np.isnan(outside).all()),
        "all_interval_extrema_peak_ratio": float(max(global_extrema)),
        "passed": bool(node_error < 1e-12 and maximum_overshoot <= tolerance
                       and derivative_violation < 1e-10 and np.isnan(outside).all()
                       and min(global_extrema) > 0)}
    if not result["passed"]:
        raise AssertionError(result)
    return result


def figure(record, first, spline):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from gabes.plot_style import apply_gabes_plot_style
    original = np.asarray(first["data"]["eom_MHz"])
    train = np.asarray(record["development"]["train_indices_original_order"])
    validation = np.asarray(record["development"]["validation_indices_original_order"])
    dense = np.linspace(3022., 3044., 2001)
    fig, axes = plt.subplots(2, 2, figsize=(12., 8.), sharex="col", layout="constrained")
    apply_gabes_plot_style(fig)
    for column, arm in enumerate(("probe", "conjugate")):
        y = np.asarray(first["data"][arm+"_ratio"])
        curve = interpolate(original, y)
        axes[0, column].plot(dense, np.exp(curve(dense)), color="#176b8a", lw=1.8,
                            label="23-point log-PCHIP sampled reference")
        axes[0, column].scatter(original[train], y[train], s=24, color="#222222", zorder=3,
                               label="Original data (development train indices)")
        axes[0, column].scatter(original[validation], y[validation], s=38,
                               facecolors="white", edgecolors="#d2691e", zorder=3,
                               label="Original data (development validation indices)")
        axes[0, column].set(yscale="log", ylabel="Pump-on / pump-off voltage ratio",
                            title=arm.capitalize()+": sampled measured-reference curve")
        for case, label, color, marker in (
                (None, "Log-PCHIP from 12 train points", "#176b8a", "o"),
                (first["development"]["candidates"][-1], "Failed Chebyshev degree 7", "#888888", "s"),
                (spline["development"]["candidates"][1], "Failed uniform cubic spline (3 interior knots)", "#9467bd", "^")):
            if case is None:
                prediction = np.exp(interpolate(original[train], y[train])(original[validation]))
            else:
                item = case["arms"][arm]
                prediction = np.asarray(item.get("prediction_all_rows", item.get(
                    "development_prediction_all_rows")))[validation]
            residual = 100*(prediction/y[validation]-1)
            axes[1, column].plot(original[validation], residual, marker=marker,
                                  color=color, markersize=4, lw=.9, label=label)
        axes[1, column].axhline(0, color="black", lw=.7)
        axes[1, column].axhline(20, color="#b22222", lw=.8, linestyle="--")
        axes[1, column].axhline(-20, color="#b22222", lw=.8, linestyle="--")
        axes[1, column].set(xlabel="Optical pump-seed separation / EOM frequency (MHz)",
                           ylabel="Development validation residual (%)",
                           title="Fixed alternating-row diagnostic; not untouched holdout")
        for ax in axes[:, column]:
            ax.grid(alpha=.25)
            ax.legend(fontsize=7, loc="best")
    fig.suptitle("Measured ratio reference: interpolation checks and failed prediction candidates", fontsize=12)
    path = OUT/"interpolation_comparison.png"
    fig.savefig(path, dpi=180)
    plt.close(fig)
    return path


def main():
    if "sampled reference 전환" not in (OUT/"fit_spec.md").read_text(encoding="utf-8"):
        raise RuntimeError("declare sampled-reference scope before building it")
    first_path, spline_path = OUT/"fit.json", OUT/"spline_fit.json"
    historical = {str(p.name): sha(p) for p in (first_path, spline_path)}
    first = json.loads(first_path.read_text(encoding="utf-8"))
    spline = json.loads(spline_path.read_text(encoding="utf-8"))
    data = first["data"]
    original_x = np.asarray(data["eom_MHz"])
    ascending = np.argsort(original_x)
    x = original_x[ascending]
    train = np.asarray(first["development"]["train_indices"])
    validation = np.asarray(first["development"]["validation_indices"])
    arms, development, source_points = {}, {}, {"eom_MHz": x.tolist()}
    for arm in ("probe", "conjugate"):
        y_original = np.asarray(data[arm+"_ratio"])
        y = y_original[ascending]
        curve = interpolate(x, y)
        arm_checks = checks(curve, x, y)
        arms[arm] = {"continuous_peak": peak_intervals(x, y), "checks": arm_checks,
                     "exact_node_metrics": metrics(np.exp(curve(x)), y)}
        train_curve = interpolate(original_x[train], y_original[train])
        prediction = np.exp(train_curve(original_x[validation]))
        development[arm] = metrics(prediction, y_original[validation]) | {
            "validation_prediction": prediction.tolist(),
            "validation_target": y_original[validation].tolist()}
        source_points[arm+"_ratio"] = y.tolist()
    passed = all(item["checks"]["passed"] for item in arms.values())
    diagnostic_passed = all(item["max_relative_error"] <= .20 for item in development.values())
    result = {"origin": "AI deterministic sampled measured-data reference under direct user calibration authorization",
        "marker": "⟦AI⟧", "date": "2026-10-16", "model_id": MODEL_ID,
        "source_workbook_sha256": first["source_workbook_sha256"],
        "extracted_ooxml_sha256": sha(OUT/"workbook_ooxml.json"),
        "domain_eom_MHz": [3022., 3044.],
        "source_points": source_points,
        "source_points_sha256": hashlib.sha256(json.dumps(source_points,
            sort_keys=True, separators=(",", ":")).encode()).hexdigest(),
        "interpolation": "PCHIP of natural log positive arm ratios, then exp; 23 source points; no fitted/free parameters",
        "source_RF_convention": "lower seed: delta_internal_MHz=3035.732439-EOM_MHz; sheet TPD=EOM_MHz-3036",
        "out_of_domain_policy": "NaN/unavailable; no extrapolation",
        "data_original_order": data, "context": first["context"], "arms": arms,
        "quality": {"sampled_interpolation_validation_passed": passed,
            "development_acceptance_passed": False,
            "status": "sampled reference interpolation; predictive development gate not passed",
            "independent_predictive_validation": False,
            "statistical_confidence_interval": None},
        "development": {"train_indices_original_order": train.tolist(),
            "validation_indices_original_order": validation.tolist(), "arms": development,
            "threshold_max_relative_error": .20,
            "alternating_interpolation_diagnostic_passed": diagnostic_passed,
            "predictive_gate_passed": False,
            "meaning": "same previously inspected cohort; diagnostic interpolation error, not independent prediction acceptance"},
        "historical_failed_fits": historical,
        "scope": "sampled measured pump-on / pump-off voltage ratio reference at supplied apparatus condition; not atomic gain/NF replacement",
        "atomic_gain_modified": False, "canonical_transfer_modified": False,
        "noise_NF_modified": False, "grand_challenge_no_fit_path_modified": False,
        "new_atomic_solves": 0,
        "fit_spec_sha256": sha(OUT/"fit_spec.md"),
        "runner_sha256": sha(Path(__file__)),
    }
    png = figure(result, first, spline)
    result["scientific_comparison_png"] = {"path": png.name, "sha256": sha(png)}
    path = OUT/"interpolation_record.json"
    path.write_text(json.dumps(result, ensure_ascii=False, indent=2)+"\n", encoding="utf-8")
    for name, prior_hash in historical.items():
        assert sha(OUT/name) == prior_hash
    print(json.dumps({"model_id": MODEL_ID, "sampled_checks_passed": passed,
        "development_diagnostic": development,
        "diagnostic_passed": diagnostic_passed,
        "predictive_acceptance": False,
        "peaks": {arm: item["continuous_peak"] for arm, item in arms.items()},
        "record_sha256": sha(path)}, indent=2))


if __name__ == "__main__":
    main()
