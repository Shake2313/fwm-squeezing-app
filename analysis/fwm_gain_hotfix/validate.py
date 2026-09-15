"""Reproduce closure benchmarks, saved-array checks, and conditional comparisons.

Run with an idle CPU after the baseline capture and production edits:
    python -m analysis.fwm_gain_hotfix.validate

This runner records held-out mismatches. It does not certify experimental gain
or physical squeezing, and never changes or refits the production closure.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import platform
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

from gabes.schemes import fwm

from .capture_before import json_default
from .snapshot_convention import snapshot_directory


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
TIERS = {"fast": fwm.FIDELITY_FAST, "balanced": fwm.FIDELITY_BALANCED}
TREND_VALUES = {
    "opd": [0.7, 0.8, 0.9, 1.0, 1.1],
    "pump_mw": [300.0, 450.0, 600.0, 750.0, 900.0],
    "temp_c": [101.0, 111.0, 121.0, 131.0, 141.0],
    "pump_waist_um": [400.0, 465.0, 530.0, 595.0, 660.0],
    "probe_waist_um": [230.0, 280.0, 330.0, 380.0, 430.0],
    "seeded_angle_deg": [0.0, 0.16, 0.32, 0.48, 0.64],
}


def gain_length(gain):
    """Lossless-amplifier diagnostic only; undefined for attenuating gain < 1."""
    return float(np.arccosh(np.sqrt(gain))) if gain >= 1.0 else None


def summarize(raw, tpd):
    op = fwm.operating_point(raw, float(tpd))
    axis = (np.asarray(raw["probe_axis_GHz"])
            - raw["raman_center_minus_GHz"]) * 1e3
    peak = int(np.argmax(raw["G_s"]))
    return {
        "tpd_mhz": float(tpd),
        "operating_point": op,
        "lossless_qL_diagnostic": gain_length(op["G_s"]),
        "scan_peak_G_s": float(raw["G_s"][peak]),
        "scan_peak_tpd_mhz": float(axis[peak]),
        "finite": bool(all(np.all(np.isfinite(raw[key]))
                           for key in ("G_s", "G_c", "S_dB"))),
        "physical_squeezing_unavailable": raw["physical_squeezing_dB"] is None,
        "quantitative_gain_supported": raw["claim_gate"]["quantitative_gain_supported"],
        "floquet_status": raw["floquet_convergence"]["status"],
        "gain_closure": raw.get("gain_closure"),
        "solved_probe_points": raw.get("response_estimator", {}).get("solved_probe_points"),
    }


def source_hashes():
    paths = [ROOT / "gabes/schemes/fwm.py", ROOT / "tests/baseline_focused.npz"]
    paths += list((ROOT / "gabes").glob("*gain*closure*"))
    return {str(path.relative_to(ROOT)).replace("\\", "/"):
            hashlib.sha256(path.read_bytes()).hexdigest() for path in paths}


def compare_saved_arrays(raw, frozen):
    """Record exact identity separately from portable numerical regression."""
    exact, errors = {}, {}
    for key in frozen.files:
        actual, expected = np.asarray(raw[key]), frozen[key]
        exact[key] = bool(np.array_equal(actual, expected))
        rtol, atol = ((0.0, 1e-5) if key in ("S_dB", "gain_referred_noise_dB")
                      else (1e-7, 1e-10))
        errors[key] = {
            "within_tolerance": bool(np.allclose(actual, expected, rtol=rtol, atol=atol)),
            "max_absolute_error": float(np.max(np.abs(actual - expected))),
            "rtol": rtol, "atol": atol,
        }
    return exact, errors


def run(repeats=7):
    before = snapshot_directory()
    capture_path = before / "capture.json"
    baseline = json.loads(capture_path.read_text(encoding="utf-8"))
    timing_path = HERE / "before/capture.json"
    timing_baseline = json.loads(timing_path.read_text(encoding="utf-8"))
    references = json.loads((HERE / "reference_points.json").read_text(encoding="utf-8"))
    scheme = fwm.FWMScheme()
    report = {
        "recorded_utc": datetime.now(timezone.utc).isoformat(),
        "environment": {"python": sys.version, "numpy": np.__version__,
                        "platform": platform.platform()},
        "source_sha256": source_hashes(),
        "snapshot_fixture": {
            "capture": capture_path.relative_to(HERE).as_posix(),
            "capture_sha256": hashlib.sha256(capture_path.read_bytes()).hexdigest(),
            "git_head": baseline["git_head"],
            "source_sha256": baseline["source_sha256"],
            "absorption_population_factors": fwm._reference_population_factors(),
        },
        "historical_timing_fixture": {
            "capture": timing_path.relative_to(HERE).as_posix(),
            "capture_sha256": hashlib.sha256(timing_path.read_bytes()).hexdigest(),
            "note": "Original shared-worktree warmed medians; absorption normalization may differ from the selected array fixture.",
        },
        "benchmarks": {}, "saved_array_regression": {},
        "saved_array_regression_errors": {}, "gold": {},
        "held_out": {}, "trends": {},
        "experimental_gain_validated": False,
        "physical_squeezing_validated": False,
    }
    for tag, tier in TIERS.items():
        params = dict(scheme.defaults(), resolution=tier)
        for enabled in (False, True):
            scheme.compute(dict(params, gain_closure_enabled=enabled))
        samples = {False: [], True: []}
        latest = {}
        for repeat in range(repeats):
            # Alternate order to reduce consistent thermal/order bias.
            for enabled in ((False, True) if repeat % 2 == 0 else (True, False)):
                started = time.perf_counter()
                latest[enabled] = scheme.compute(dict(params, gain_closure_enabled=enabled))
                samples[enabled].append(time.perf_counter() - started)
        old_s = timing_baseline["cases"][f"{tag}_default"]["seconds_median"]
        off_s, on_s = (float(np.median(samples[value])) for value in (False, True))
        report["benchmarks"][tag] = {
            "repeats": repeats, "off_seconds": samples[False], "on_seconds": samples[True],
            "off_median_seconds": off_s, "on_median_seconds": on_s,
            "frozen_before_median_seconds": old_s,
            "on_vs_current_off_ratio": on_s / off_s,
            "on_vs_frozen_before_ratio": on_s / old_s,
            "overhead_milliseconds": (on_s - off_s) * 1e3,
            "below_historical_interactive_budget": on_s <= (0.25 if tag == "fast" else 0.4),
        }
        with np.load(before / f"{tag}_default.npz") as frozen:
            equal, errors = compare_saved_arrays(latest[False], frozen)
        report["saved_array_regression"][f"{tag}_off"] = equal
        report["saved_array_regression_errors"][f"{tag}_off"] = errors
        report["gold"][tag] = {
            "off": summarize(latest[False], -8.0), "on": summarize(latest[True], -8.0)}
        print(f"{tag}: on {on_s:.4f}s / off {off_s:.4f}s / before {old_s:.4f}s; "
              f"on G={fwm.operating_point(latest[True], -8.0)['G_s']:.5g}", flush=True)

    for name in ("ultra_default", "ultra_detuned"):
        params = dict(baseline["cases"][name]["params"], gain_closure_enabled=True)
        raw = scheme.compute(params)
        with np.load(before / f"{name}.npz") as frozen:
            exact, errors = compare_saved_arrays(raw, frozen)
        report["saved_array_regression"][name] = exact
        report["saved_array_regression_errors"][name] = errors
        print(f"{name}: all arrays unchanged="
              f"{all(report['saved_array_regression'][name].values())}", flush=True)

    for point in references["held_out"]:
        comparison = {"reference": point, "tiers": {}}
        target = point["targets"]["probe_gain"]
        for tag, tier in TIERS.items():
            params = dict(scheme.defaults(), resolution=tier, **point["params"])
            row = {}
            for enabled in (False, True):
                raw = scheme.compute(dict(params, gain_closure_enabled=enabled))
                stats = summarize(raw, params["tpd"])
                stats["estimate_over_reference"] = stats["operating_point"]["G_s"] / target
                row["on" if enabled else "off"] = stats
            comparison["tiers"][tag] = row
        report["held_out"][point["id"]] = comparison
        print(f"held-out {point['id']}: " + "; ".join(
            f"{tag} G={row['on']['operating_point']['G_s']:.5g} (reference {target})"
            for tag, row in comparison["tiers"].items()), flush=True)

    # Full deterministic Balanced scan avoids adaptive-node changes obscuring
    # parameter dependence. All other parameters stay at the Gold input values.
    for key, values in TREND_VALUES.items():
        rows = []
        for value in values:
            params = dict(scheme.defaults(), resolution=fwm.FIDELITY_BALANCED)
            params[key] = value
            rows.append({"parameter_value": value, **{
                "on" if enabled else "off": summarize(scheme.compute(dict(
                    params, gain_closure_enabled=enabled)), params["tpd"])
                for enabled in (False, True)}})
        old = [row["off"]["operating_point"]["G_s"] for row in rows]
        new = [row["on"]["operating_point"]["G_s"] for row in rows]
        report["trends"][key] = {
            "rows": rows,
            "adjacent_direction_preserved": (np.sign(np.diff(old))
                                              == np.sign(np.diff(new))).tolist(),
            "note": "Fixed operating detuning; gain peaks may move. Direction agreement is a comparison, not validation.",
        }
        print(f"trend {key}: off={np.round(old, 3)} on={np.round(new, 3)}", flush=True)
    report["all_saved_arrays_unchanged"] = all(
        all(values.values()) for values in report["saved_array_regression"].values())
    report["all_saved_arrays_within_tolerance"] = all(
        result["within_tolerance"]
        for case in report["saved_array_regression_errors"].values()
        for result in case.values())
    report["source_sha256_at_completion"] = source_hashes()
    report["source_unchanged_during_run"] = (
        report["source_sha256"] == report["source_sha256_at_completion"])
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repeats", type=int, default=7)
    parser.add_argument("--output", type=Path, default=HERE / "validation_report.json")
    args = parser.parse_args()
    if args.repeats < 5:
        parser.error("Use at least five warmed repeats")
    report = run(args.repeats)
    args.output.write_text(json.dumps(report, indent=2, ensure_ascii=False,
                                     default=json_default) + "\n", encoding="utf-8")
    print(f"Report: {args.output}", flush=True)
    return 0 if (report["all_saved_arrays_within_tolerance"]
                 and report["source_unchanged_during_run"]) else 1


if __name__ == "__main__":
    raise SystemExit(main())
