"""Completion gate for the FWM Squeezing Fast/Balanced remaster (see GATE.md).

    python -m analysis.fwm_lite.gate --output analysis/fwm_lite/gate_report.json
    python -m analysis.fwm_lite.gate --output NEW.json --skip-pytest   # iteration only

Exit status 0 only when every G1–G5 item passes.  Timings assume nothing else is
running; references are recomputed live from the unchanged Ultra tier.
"""
import argparse
import json
import os
import platform
import re
import subprocess
import sys
import tempfile
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from gabes.schemes import fwm  # noqa: E402

CASES = {
    "default":     dict(),
    "v6_opt":      dict(opd=-1.5, temp_c=110.0, pump_mw=600.0),
    "frontier":    dict(opd=-2.2, temp_c=135.0, pump_mw=600.0),
    "near_res":    dict(opd=-0.2, temp_c=100.0, pump_mw=200.0),
    "blue_res":    dict(opd=0.3,  temp_c=80.0,  pump_mw=300.0),
    "blue_far":    dict(opd=2.5,  temp_c=150.0, pump_mw=1200.0),
    "red_edge":    dict(opd=-3.0, temp_c=130.0, pump_mw=600.0),
    "cold_weak":   dict(opd=0.9,  temp_c=60.0,  pump_mw=50.0),
    "hot":         dict(opd=0.9,  temp_c=150.0, pump_mw=600.0),
    "low_pump":    dict(opd=-1.5, temp_c=120.0, pump_mw=135.0),
    "bright_seed": dict(opd=0.9,  temp_c=121.0, pump_mw=600.0, probe_uw=200.0),
    "long_cell":   dict(opd=0.9,  temp_c=100.0, pump_mw=600.0, cell_mm=50.0,
                        seeded_angle_deg=1.0),
}
LIMITS = {
    "fast": dict(tube=0.10, p95=0.03, logG_tube=0.010, op_dS=0.05, op_G=0.01,
                 min_val=0.05, peak=0.02, pearson=0.998, t_default=0.12, t_any=0.25,
                 speedup=4.5),
    "balanced": dict(tube=0.05, p95=0.015, logG_tube=0.005, op_dS=0.02, op_G=0.005,
                     min_val=0.03, peak=0.01, pearson=0.9995, t_default=0.25, t_any=0.40,
                     speedup=7.0),
}
TIERS = {"fast": fwm.FIDELITY_FAST, "balanced": fwm.FIDELITY_BALANCED}
SCHEME = fwm.FWMScheme()


def params_for(case, resolution):
    # This gate tests the 2026-09-11 common-model solver acceleration. The later
    # gain closure has independent physics/calibration gates under fwm_gain_hotfix.
    params = dict(SCHEME.defaults(), mode=fwm.MODE_SEEDED, resolution=resolution,
                  gain_closure_enabled=False)
    params.update(CASES[case])
    return params


def timed(fn, repeats):
    samples = []
    out = None
    for _ in range(repeats):
        started = time.perf_counter()
        out = fn()
        samples.append(time.perf_counter() - started)
    return float(np.median(samples)), samples, out


def shape(cand, ref, tpd, window_mhz=500.0, tube_mhz=1.0):
    x = (np.asarray(ref["probe_axis_GHz"]) - ref["raman_center_minus_GHz"]) * 1e3
    win = np.abs(x) <= window_mhz
    S_u, S_c = np.asarray(ref["S_dB"]), np.asarray(cand["S_dB"])
    lG_u = np.log10(np.asarray(ref["G_s"]))
    lG_c = np.log10(np.asarray(cand["G_s"]))
    shifts = np.linspace(-tube_mhz, tube_mhz, 9)
    tube_S = np.min([np.abs(np.interp(x + s, x, S_c) - S_u) for s in shifts], axis=0)
    tube_G = np.min([np.abs(np.interp(x + s, x, lG_c) - lG_u) for s in shifts], axis=0)
    op_u, op_c = fwm.operating_point(ref, tpd), fwm.operating_point(cand, tpd)
    i_min = int(np.argmin(np.where(win, S_u, np.inf)))

    def pearson(a, b, min_range):
        a, b = a[win], b[win]
        if np.ptp(a) < min_range:
            return None
        return float(np.corrcoef(a, b)[0, 1])

    return dict(
        tube=float(tube_S[win].max()),
        p95=float(np.percentile(np.abs(S_c - S_u)[win], 95)),
        logG_tube=float(tube_G[win].max()),
        op_dS=float(abs(op_c["S_dB"] - op_u["S_dB"])),
        op_G=float(abs(op_c["G_s"] / op_u["G_s"] - 1.0)),
        min_val=float(max(abs(S_c[win].min() - S_u[win].min()), abs(S_c[i_min] - S_u[i_min]))),
        peak=float(abs(np.max(np.asarray(cand["G_s"])[win]) / np.max(np.asarray(ref["G_s"])[win]) - 1.0)),
        pearson_S=pearson(S_u, S_c, 0.01),
        pearson_logG=pearson(lG_u, lG_c, 1e-3),
    )


def shape_pass(metrics, lim):
    checks = {k: metrics[k] <= lim[k] for k in ("tube", "p95", "logG_tube", "op_dS", "op_G",
                                                 "min_val", "peak")}
    for key in ("pearson_S", "pearson_logG"):
        checks[key] = metrics[key] is None or metrics[key] >= lim["pearson"]
    return checks


def git_unchanged(paths):
    result = subprocess.run(["git", "diff", "--quiet", "HEAD", "--", *paths], cwd=ROOT)
    return result.returncode == 0


def cold_start_seconds():
    code = ("import time, json; t0=time.perf_counter(); from gabes.schemes import fwm; "
            "t1=time.perf_counter(); s=fwm.FWMScheme(); p=s.defaults(); "
            "p['resolution']=fwm.FIDELITY_FAST; p['gain_closure_enabled']=False; "
            "s.compute(p); t2=time.perf_counter(); "
            "print(json.dumps({'import_s': t1-t0, 'first_fast_call_s': t2-t1}))")
    out = subprocess.run([sys.executable, "-c", code], cwd=ROOT, capture_output=True,
                         text=True, encoding="utf-8", check=True)
    return json.loads(out.stdout.strip().splitlines()[-1])


def numba_disabled_arrays():
    code = "\n".join((
        "import sys",
        "import numpy as np",
        "from gabes.schemes import fwm",
        "assert not fwm.kernels.available()",
        "s = fwm.FWMScheme()",
        "out = {}",
        "for tag, tier in (('fast', fwm.FIDELITY_FAST), ('balanced', fwm.FIDELITY_BALANCED)):",
        "    p = s.defaults()",
        "    p['resolution'] = tier",
        "    p['gain_closure_enabled'] = False",
        "    raw = s.compute(p)",
        "    for k in ('G_s', 'G_c', 'S_dB'):",
        "        out[tag + '__' + k] = np.asarray(raw[k])",
        "np.savez(sys.argv[1], **out)",
    ))
    with tempfile.TemporaryDirectory() as tmp:
        target = Path(tmp) / "nonumba.npz"
        env = dict(os.environ, GABES_DISABLE_NUMBA="1")
        subprocess.run([sys.executable, "-c", code, str(target)], cwd=ROOT, env=env, check=True)
        with np.load(target) as data:
            return {k: data[k] for k in data.files}


def run_pytest(extra_env=None, args=("-q", "-p", "no:cacheprovider")):
    env = dict(os.environ, **(extra_env or {}))
    out = subprocess.run([sys.executable, "-m", "pytest", *args], cwd=ROOT, env=env,
                         capture_output=True, text=True, encoding="utf-8", errors="replace")
    tail = out.stdout.strip().splitlines()[-40:]
    summary = next((line for line in reversed(tail) if " passed" in line or " failed" in line), "")
    failed = sorted(set(re.findall(r"FAILED (\S+)", out.stdout)))
    return dict(returncode=out.returncode, summary=summary, failed=failed)


def build(skip_pytest):
    report = {"gate": "GATE.md v1", "environment": {
        "python": platform.python_version(), "numpy": np.__version__,
        "numba_available": fwm.kernels.available(), "cpu_count": os.cpu_count()}}
    items = {}

    # ---- G3.7 Ultra unchanged ---------------------------------------------------------
    ultra_params = params_for("default", fwm.FIDELITY_ULTRA)
    ultra_default = SCHEME.compute(ultra_params)
    with np.load(HERE / "before" / "ultra_default_before.npz") as before:
        identical = all(np.array_equal(np.asarray(ultra_default[k]), before[k])
                        for k in ("probe_axis_GHz", "G_s", "G_c", "S_dB"))
    ultra_dict = fwm.FWM_FIDELITY[fwm.FIDELITY_ULTRA] == dict(
        coarse_points=401, velocity_step=1.0, velocity_cutoff=4.0, phase_detail=fwm.PHASE_ULTRA)
    items["G3.7"] = dict(passed=bool(identical and ultra_dict), bit_identical=identical,
                         tier_dict_unchanged=ultra_dict)

    # ---- references + candidates -----------------------------------------------------
    for tier in TIERS.values():                      # warm JIT/caches for both tiers
        SCHEME.compute(params_for("default", tier))
    cases = {}
    for case in CASES:
        started = time.perf_counter()
        ref = ultra_default if case == "default" else SCHEME.compute(params_for(case, fwm.FIDELITY_ULTRA))
        ultra_s = time.perf_counter() - started
        row = {"ultra_seconds_single": ultra_s,
               "ultra_floquet_status": ref["floquet_convergence"]["status"]}
        for tag, tier in TIERS.items():
            params = params_for(case, tier)
            median, samples, raw = timed(lambda: SCHEME.compute(params), 5)
            est = raw["response_estimator"]
            metrics = shape(raw, ref, params["tpd"])
            checks = shape_pass(metrics, LIMITS[tag])
            integrity = {
                "G3.1_axis_equal": bool(np.array_equal(raw["probe_axis_GHz"], ref["probe_axis_GHz"])),
                "G3.2_finite": bool(all(np.all(np.isfinite(raw[k])) for k in ("G_s", "G_c", "S_dB"))),
                "G3.3_audit": bool(raw["floquet_convergence"]["status"] == "CONVERGED"
                                   and raw["floquet_convergence"]["high_order"] == 3
                                   and raw["floquet_convergence"]["comparison_order"] == 2
                                   and raw["floquet_convergence"].get("solved_scan_points")
                                   == est["solved_probe_points"]
                                   and (tag != "balanced" or est["solved_probe_points"] == 401)),
                "G3.4_guard": bool(est["pole_guard_max_relative"] <= 1e-6),
                "G3.6_provenance": bool(
                    {"method", "velocity_average", "solved_probe_points", "audit_scope",
                     "interpolation"} <= set(est)
                    and (tag != "fast" or est["interpolation"] == "none"
                         or any("interpolated" in r for r in raw["claim_gate"]["reasons"]))),
            }
            row[tag] = dict(seconds_median=median, seconds=samples, metrics=metrics,
                            checks=checks, integrity=integrity,
                            solved_probe_points=est["solved_probe_points"],
                            adaptive_rounds=est.get("adaptive_rounds"),
                            guard_max_relative=est["pole_guard_max_relative"],
                            guard_fallback_rows=est["pole_guard_fallback_rows"],
                            pole_min_imag_MHz=est["pole_min_distance_from_real_axis_MHz"])
        cases[case] = row
        print(f"{case:12s} ultra {ultra_s:5.2f}s | " + " | ".join(
            f"{tag} {row[tag]['seconds_median']*1e3:6.1f} ms n={row[tag]['solved_probe_points']:3d} "
            f"tube {row[tag]['metrics']['tube']:.4f} p95 {row[tag]['metrics']['p95']:.4f} "
            f"{'PASS' if all(row[tag]['checks'].values()) and all(row[tag]['integrity'].values()) else 'FAIL'}"
            for tag in TIERS), flush=True)
    report["cases"] = cases

    before = json.loads((HERE / "before" / "tier_timings_before.json").read_text(encoding="utf-8"))
    for tag in TIERS:
        default_t = cases["default"][tag]["seconds_median"]
        worst_t = max(cases[c][tag]["seconds_median"] for c in CASES)
        before_t = float(np.median([t[0] for t in before[tag]["times"]]))
        lim = LIMITS[tag]
        item = "G1.1" if tag == "fast" else "G1.2"
        items[item] = dict(passed=bool(default_t <= lim["t_default"] and worst_t <= lim["t_any"]),
                           default_s=default_t, worst_case_s=worst_t,
                           limits=(lim["t_default"], lim["t_any"]))
        items.setdefault("G1.3", {"passed": True})
        items["G1.3"][tag] = dict(before_s=before_t, after_s=default_t,
                                  speedup=before_t / default_t, required=lim["speedup"])
        items["G1.3"]["passed"] = bool(items["G1.3"]["passed"] and before_t / default_t >= lim["speedup"])
        items[f"G2.{tag}"] = dict(passed=bool(all(all(cases[c][tag]["checks"].values()) for c in CASES)),
                                  failing=[(c, k) for c in CASES for k, v in cases[c][tag]["checks"].items() if not v])
        for key in ("G3.1_axis_equal", "G3.2_finite", "G3.3_audit", "G3.4_guard", "G3.6_provenance"):
            name = key.split("_")[0]
            ok = all(cases[c][tag]["integrity"][key] for c in CASES)
            items.setdefault(name, {"passed": True})
            items[name][tag] = ok
            items[name]["passed"] = bool(items[name]["passed"] and ok)
    items["ultra_reference_all_converged"] = dict(
        passed=all(cases[c]["ultra_floquet_status"] == "CONVERGED" for c in CASES))

    cold = cold_start_seconds()
    items["G1.4"] = dict(passed=bool(cold["first_fast_call_s"] <= 3.0), **cold)

    nonumba = numba_disabled_arrays()
    deviations = {}
    for tag, tier in TIERS.items():
        raw = SCHEME.compute(params_for("default", tier))
        for key in ("G_s", "G_c", "S_dB"):
            a, b = nonumba[f"{tag}__{key}"], np.asarray(raw[key])
            deviations[f"{tag}__{key}"] = float(np.max(np.abs(a - b) / np.maximum(np.abs(b), 1e-12)))
    items["G3.8"] = dict(passed=bool(max(deviations.values()) <= 1e-9), max_relative=deviations)

    items["G4.3"] = dict(passed=git_unchanged(["tests/baseline_focused.npz",
                                              "tests/baseline_focused_manifest.json"]))
    if skip_pytest:
        items["G4.1"] = dict(passed=False, skipped=True)
        items["G4.2"] = dict(passed=False, skipped=True)
    else:
        before_failures = set(re.findall(r"FAILED (\S+)",
                                         (HERE / "before" / "pytest_before_summary.txt").read_text(encoding="utf-8")))
        full = run_pytest()
        new_failures = sorted(set(full["failed"]) - before_failures)
        items["G4.1"] = dict(passed=not new_failures, summary=full["summary"],
                             new_failures=new_failures, pre_existing=sorted(before_failures))
        targeted = ["tests/test_pole_doppler.py", "tests/test_fwm_fast_tiers.py"]
        on = run_pytest(args=("-q", "-p", "no:cacheprovider", *targeted))
        off = run_pytest({"GABES_DISABLE_NUMBA": "1"}, args=("-q", "-p", "no:cacheprovider", *targeted))
        items["G4.2"] = dict(passed=on["returncode"] == 0 and off["returncode"] == 0,
                             numba=on["summary"], no_numba=off["summary"])

    devlog = HERE / "DEVLOG.md"
    text = devlog.read_text(encoding="utf-8") if devlog.exists() else ""
    required = ("## 기각", "## 재사용", "## 측정", "## 기법")
    items["G5.1"] = dict(passed=all(h in text for h in required), required_headings=required)
    report["items"] = items
    report["passed"] = bool(all(v.get("passed") for k, v in items.items()))
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--skip-pytest", action="store_true")
    args = parser.parse_args()
    report = build(args.skip_pytest)
    report["items"]["G5.2"] = dict(passed=report["passed"], note="exit status mirrors this flag")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, ensure_ascii=False, indent=1, default=float) + "\n",
                           encoding="utf-8")
    for key, value in report["items"].items():
        print(f"{key:32s} {'PASS' if value.get('passed') else 'FAIL'}")
    print("GATE", "PASS" if report["passed"] else "FAIL", "->", args.output)
    return 0 if report["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
