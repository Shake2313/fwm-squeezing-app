"""Reproduce the frozen gain fit and validate the public Fast/Balanced readout.

Run from the repository root with --capture-before before changing production,
--fit to refit frozen baseline data, or no arguments to validate the hotfix.
No fitting runs in the production solver. Historical evidence is read only.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import sys
import time

OUT = Path(__file__).resolve().parent
ROOT = OUT.parents[2]
sys.path.insert(0, str(ROOT))
os.environ.setdefault("NUMBA_CACHE_DIR", str(OUT / "numba_cache"))
os.environ.setdefault("MPLCONFIGDIR", str(OUT / "mpl_cache"))

import numpy as np
from scipy.optimize import linprog
from gabes.schemes import fwm
from analysis.fwm_gain_hotfix.snapshot_convention import snapshot_directory

EOM = np.arange(3039, 3029, -1, dtype=float) / 1000
SEED_OFF_UW = np.array([3.7, 3.7, 3.7, 3.8, 3.7, 3.7, 3.7, 3.8, 3.7, 3.7])
PROBE_UW = np.array([49.7, 56.2, 91, 121, 175, 215, 280, 353, 400, 439.])
CONJ_UW = np.array([47.5, 55, 90, 121, 175, 215, 280, 353, 400, 450.])


def save(name, data):
    (OUT / name).write_text(json.dumps(data, indent=2, ensure_ascii=False) + "\n",
                           encoding="utf-8")


def params(tier, *, gold=False, enabled=True):
    p = dict(fwm.FWMScheme().defaults(), resolution=tier,
             gain_closure_enabled=enabled)
    if not gold:
        p.update(temp_c=118., pump_mw=380., opd=1.06, probe_uw=3.7)
    return p


def exact(tier, *, enabled=True, efficiency=0.8694):
    return fwm.compute_spectrum(
        1.06, T=391.15, P_pump=0.380, P_probe=3.7e-6,
        L=0.0125, w_pump=530e-6, w_probe=330e-6,
        pump_probe_angle_deg=0.32, transit_rate=2*np.pi*100e3,
        branch=-1, line_strength=0.74, floquet_order=3,
        phase_detail=fwm.PHASE_ULTRA, model_fidelity=tier,
        response_method=fwm.RESPONSE_POLE, gain_closure_enabled=enabled,
        detection_efficiency=efficiency, velocity_step=1., velocity_cutoff=4.,
        coarse_points=10, fine_points=0, scan_min=1.06-EOM[0],
        scan_max=1.06-EOM[-1])


def numeric(raw):
    return {k: np.asarray(raw[k]).tolist() for k in
            ("probe_axis_GHz", "G_s", "G_c", "G_s_smallsignal", "G_c_smallsignal",
             "S_dB", "gain_referred_noise_dB")}


def capture():
    cases = {}
    scheme = fwm.FWMScheme()
    for tier in (fwm.FIDELITY_FAST, fwm.FIDELITY_BALANCED):
        name = tier.split()[0]
        cases["exact_" + name] = numeric(exact(tier))
        cases["off_" + name] = numeric(exact(tier, enabled=False))
        cases["gold_" + name] = numeric(scheme.compute(params(tier, gold=True)))
        cases["ui_" + name] = numeric(scheme.compute(params(tier)))
    cases["exact_Ultra"] = numeric(exact(fwm.FIDELITY_ULTRA))
    sources = {name: hashlib.sha256((ROOT / name).read_bytes()).hexdigest()
               for name in ("gabes/fwm_gain_closure.py", "gabes/schemes/fwm.py",
                            "gabes/observables.py", "gabes/constants.py")}
    save("before.json", {"cases": cases, "source_sha256": sources,
                         "devlog_original_bytes": (OUT.parent / "DEVLOG.md").stat().st_size})
    historical = {}
    for directory in (ROOT / "analysis/fwm_gain_hotfix",
                      ROOT / "analysis/tpd_gain_diagnostic_20261001"):
        for path in directory.rglob("*"):
            if (path.is_file() and OUT not in path.parents
                    and not any(part in ("__pycache__", "numba_cache", "mpl_cache")
                                for part in path.parts)):
                historical[path.relative_to(ROOT).as_posix()] = hashlib.sha256(
                    path.read_bytes()).hexdigest()
    save("historical_manifest.json", historical)
    print("Captured baseline and", len(historical), "historical files", flush=True)


def fit():
    baseline = json.loads((OUT / "before.json").read_text(encoding="utf-8"))
    result = {"method": "monotone continuous two-piece affine; minimax relative residual",
              "statistical_confidence_interval": None,
              "independent_validation": False,
              "conditions": params(fwm.FIDELITY_BALANCED),
              "eom_GHz": EOM.tolist(), "pump_off_seed_uW": SEED_OFF_UW.tolist(),
              "probe_uW": PROBE_UW.tolist(), "conjugate_uW": CONJ_UW.tolist(),
              "simulation_cell_input_seed_uW": 3.7}
    for name, knot, powers, key in (("probe", 7., PROBE_UW, "G_s"),
                                    ("conjugate", 6.1, CONJ_UW, "G_c")):
        x = np.array(baseline["cases"]["exact_Balanced"][key])
        y = powers / SEED_OFF_UW
        basis = np.column_stack((x, np.ones(x.size), np.maximum(x-knot, 0.)))
        # Relative error t, nonnegative slopes on both sides of the fixed knee.
        A = np.vstack((np.column_stack((basis, -y)),
                       np.column_stack((-basis, -y)), [-1., 0., -1., 0.]))
        b = np.r_[y, -y, 0.]
        bounds = [(0., None), ((0., 0.) if name == "conjugate" else (None, None)),
                  (None, None), (0., None)]
        solved = linprog([0., 0., 0., 1.], A_ub=A, b_ub=b,
                         bounds=bounds, method="highs")
        if not solved.success:
            raise RuntimeError(solved.message)
        a, intercept, hinge, error = solved.x
        prediction = basis @ solved.x[:3]
        result[name] = {"knot": knot, "a": float(a), "b": float(intercept),
                        "c": float(hinge), "max_relative_residual": float(error),
                        "baseline_gain": x.tolist(), "measured_gain": y.tolist(),
                        "prediction": prediction.tolist(),
                        "rmse": float(np.sqrt(np.mean((prediction-y)**2)))}
    save("fit.json", result)
    print(json.dumps({name: result[name] for name in ("probe", "conjugate")}, indent=2))


def validate():
    # The shared worktree and clean parent use different absorption/noise
    # conventions. Keep both frozen baselines instead of changing tolerances.
    baseline_name = ("before_parent.json" if snapshot_directory().name == "before_parent"
                     else "before.json")
    baseline = json.loads((OUT / baseline_name).read_text(encoding="utf-8"))
    result = {"cases": {}, "bypass_exact_equal": {}, "historical_preservation": {}}
    scheme = fwm.FWMScheme()
    for tier in (fwm.FIDELITY_FAST, fwm.FIDELITY_BALANCED):
        name = tier.split()[0]
        t = time.perf_counter()
        raw = exact(tier)
        ui = scheme.compute(params(tier))
        for label, data in (("exact_", raw), ("ui_", ui)):
            axis = (fwm.constants.NU_HF / 1e9 - EOM)*1000
            gs = np.array([fwm.operating_point(data, delta)["G_s"] for delta in axis])
            gc = np.array([fwm.operating_point(data, delta)["G_c"] for delta in axis])
            result["cases"][label+name] = {
                "G_s": gs.tolist(), "G_c": gc.tolist(),
                "probe_max_relative_error": float(np.max(np.abs(gs/(PROBE_UW/SEED_OFF_UW)-1))),
                "conjugate_max_relative_error": float(np.max(np.abs(gc/(CONJ_UW/SEED_OFF_UW)-1))),
                "gain_closure": data["gain_closure"],
                "physical_squeezing_dB": data["physical_squeezing_dB"],
                "response_estimator": data["response_estimator"],
            }
        for label, raw_bypass in (("off_", exact(tier, enabled=False)),
                                   ("gold_", scheme.compute(params(tier, gold=True)))):
            result["bypass_exact_equal"][label+name] = all(
                np.array_equal(np.asarray(value), baseline["cases"][label+name][key])
                for key, value in numeric(raw_bypass).items())
        low_eta = exact(tier, efficiency=0.5)
        result["bypass_exact_equal"]["efficiency_"+name] = all(
            np.array_equal(raw[key], low_eta[key]) for key in ("G_s", "G_c"))
        ui_low_eta = scheme.compute(dict(params(tier), detection_eff_pct=50.))
        result["bypass_exact_equal"]["ui_efficiency_"+name] = all(
            np.array_equal(ui[key], ui_low_eta[key])
            for key in ("probe_axis_GHz", "G_s", "G_c"))
        result["cases"]["ui_"+name]["display_points"] = len(ui["G_s"])
        print(tier, result["cases"]["exact_"+name]["probe_max_relative_error"],
              "UI", result["cases"]["ui_"+name]["probe_max_relative_error"],
              "seconds", round(time.perf_counter()-t, 2), flush=True)
    result["bypass_exact_equal"]["Ultra"] = all(
        np.array_equal(value, baseline["cases"]["exact_Ultra"][key])
        for key, value in numeric(exact(fwm.FIDELITY_ULTRA)).items())
    manifest = json.loads((OUT / "historical_manifest.json").read_text(encoding="utf-8"))
    unpackaged_cache_entries = []
    for name, digest in manifest.items():
        # The explicitly requested consolidated log gains an appended section.
        if name == "analysis/fwm_gain_hotfix/DEVLOG.md":
            prefix = (ROOT / name).read_bytes()[:baseline["devlog_original_bytes"]]
            assert hashlib.sha256(prefix).hexdigest() == digest, "Original log changed"
            continue
        path = ROOT / name
        if not path.is_file() and name.startswith("analysis/tpd_gain_diagnostic_20261001/"):
            # The user consolidated these records after the original run.
            # Resolve their new location without rewriting the frozen manifest.
            path = ROOT / name.replace("analysis/tpd_gain_diagnostic_20261001/",
                                       "analysis/fwm_gain_hotfix/tpd_gain_diagnostic_20261001/", 1)
        if not path.is_file() and "/pytest_cache/" in name:
            # Generated pytest cache is not distributed with the scientific
            # archive. Preserve its original manifest entry but report the skip.
            unpackaged_cache_entries.append(name)
            continue
        if not path.is_file() or hashlib.sha256(path.read_bytes()).hexdigest() != digest:
            raise AssertionError("Historical evidence changed: " + name)
    result["historical_preservation"] = {
        "checked": len(manifest)-1-len(unpackaged_cache_entries), "unchanged": True,
        "unpackaged_generated_cache_entries": unpackaged_cache_entries}
    assert all(result["bypass_exact_equal"].values()), result["bypass_exact_equal"]
    for label, case in result["cases"].items():
        limit = 0.08
        assert case["probe_max_relative_error"] < limit, label
        assert case["conjugate_max_relative_error"] < limit, label
    save("validation.json", result)
    with (OUT / "comparison.csv").open("w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["EOM_GHz", "measured_probe_gain", "measured_conjugate_gain",
                         "exact_Fast_probe", "exact_Balanced_probe", "UI_Fast_probe",
                         "UI_Balanced_probe", "exact_Balanced_conjugate"])
        for i in range(10):
            writer.writerow([EOM[i], PROBE_UW[i]/SEED_OFF_UW[i], CONJ_UW[i]/SEED_OFF_UW[i]] +
                            [result["cases"][k]["G_s"][i] for k in
                             ("exact_Fast", "exact_Balanced", "ui_Fast", "ui_Balanced")] +
                            [result["cases"]["exact_Balanced"]["G_c"][i]])
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 2, figsize=(10, 4), layout="constrained")
    for ax, name, measured in ((axes[0], "G_s", PROBE_UW/SEED_OFF_UW),
                                (axes[1], "G_c", CONJ_UW/SEED_OFF_UW)):
        ax.scatter(EOM, measured, color="black", label="Measured (fit data)")
        ax.plot(EOM, baseline["cases"]["exact_Balanced"][name], label="Previous closure")
        ax.plot(EOM, result["cases"]["exact_Balanced"][name], label="Hotfix exact")
        ax.plot(EOM, result["cases"]["ui_Fast"][name], "--", label="Fast display")
        ax.set(xlabel="EOM frequency (GHz)", ylabel="Power gain (" + name + ")")
        ax.grid(alpha=0.2)
        ax.legend(fontsize=8)
    fig.savefig(OUT / "comparison.png", dpi=160)
    print("Validation passed; historical evidence unchanged", flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--capture-before", action="store_true")
    parser.add_argument("--fit", action="store_true")
    args = parser.parse_args()
    if args.capture_before:
        capture()
    elif args.fit:
        fit()
    else:
        validate()
