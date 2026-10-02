"""Read-only FWM diagnosis; writes only into this new diagnostic directory.

Run from repository root: python analysis/tpd_gain_diagnostic_20261001/run_diagnostic.py
Existing solver, defaults, calibration and measurements are never changed.
"""
from __future__ import annotations

import csv
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import time

OUT = Path(__file__).resolve().parent
ROOT = OUT.parent.parent
sys.path.insert(0, str(ROOT))
os.environ.setdefault("PYTHONDONTWRITEBYTECODE", "1")
os.environ.setdefault("NUMBA_CACHE_DIR", str(OUT / "numba_cache"))
os.environ.setdefault("MPLCONFIGDIR", str(OUT / "mpl_cache"))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import scipy
from scipy.optimize import minimize_scalar
from gabes import constants, hyperfine, kernels
from gabes.schemes import fwm

EOM = np.arange(3039, 3029, -1, dtype=float) / 1000
SEED = np.array([3.7, 3.7, 3.7, 3.8, 3.7, 3.7, 3.7, 3.8, 3.7, 3.7])
PROBE = np.array([49.7, 56.2, 91, 121, 175, 215, 280, 353, 400, 439.])
CONJ = np.array([47.5, 55, 90, 121, 175, 215, 280, 353, 400, 450.])
MEAS = PROBE / SEED
DELTA = (constants.NU_HF / 1e9 - EOM) * 1000
SCHEME = fwm.FWMScheme()
DEFAULTS = SCHEME.defaults()
PARAMS = dict(DEFAULTS, temp_c=118., pump_mw=380., opd=1.06, probe_uw=3.7)
RESULTS = {}


def plain(value):
    if isinstance(value, np.ndarray):
        return plain(value.tolist())
    if isinstance(value, np.generic):
        return plain(value.item())
    if isinstance(value, complex):
        return {"real": value.real, "imag": value.imag}
    if isinstance(value, dict):
        return {str(k): plain(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [plain(v) for v in value]
    if isinstance(value, float) and not np.isfinite(value):
        return str(value)
    return value


def save_json(name, value):
    (OUT / name).write_text(json.dumps(plain(value), ensure_ascii=False, indent=2), encoding="utf-8")


def source_snapshot():
    files = subprocess.check_output(["git", "ls-files", "-z"], cwd=ROOT).decode("utf-8").split("\0")
    extra = ["gabes/_csv_numeric.py", "sabes/noise.py", "pytest.ini"]
    hashes = {}
    for name in sorted(set(files + extra)):
        path = ROOT / name
        if name and path.is_file():
            hashes[name] = hashlib.sha256(path.read_bytes()).hexdigest()
    return hashes


def fits(g):
    g = np.asarray(g)
    a = float(g @ MEAS / (g @ g))
    aa, b = np.linalg.lstsq(np.column_stack((g, np.ones(g.size))), MEAS, rcond=None)[0]
    out = {}
    for name, pred, coef in (("identity", g, {"a": 1., "b": 0.}),
                             ("scale", a*g, {"a": a, "b": 0.}),
                             ("affine", aa*g+b, {"a": aa, "b": b})):
        residual = MEAS - pred
        out[name] = dict(**coef, rmse=float(np.sqrt(np.mean(residual**2))),
                         mae=float(np.mean(np.abs(residual))),
                         max_abs=float(np.max(np.abs(residual))),
                         max_abs_relative_pct=float(100*np.max(np.abs(residual/MEAS))),
                         r2=float(1-np.sum(residual**2)/np.sum((MEAS-MEAS.mean())**2)),
                         prediction=pred, residual_meas_minus_prediction=residual,
                         relative_meas_minus_prediction_pct=100*residual/MEAS)
    return out


def kwargs(p):
    return dict(T=p["temp_c"]+273.15, P_pump=p["pump_mw"]*1e-3,
                P_probe=p["probe_uw"]*1e-6, L=p["cell_mm"]*1e-3,
                w_pump=p["pump_waist_um"]*1e-6, w_probe=p["probe_waist_um"]*1e-6,
                line_strength=p["line_strength"], mode_overlap_penalty=p["mode_overlap_penalty"],
                polarization_penalty=p["polarization_penalty"],
                zeeman_participation_penalty=p["zeeman_participation_penalty"],
                pump_probe_angle_deg=p["seeded_angle_deg"],
                detection_efficiency=p["detection_eff_pct"]/100,
                loss_frac=p["loss_pct"]/100, qe=p["qe_pct"]/100,
                transit_rate=2*np.pi*p["transit_rate_khz"]*1e3,
                branch=-1, phase_detail=fwm.PHASE_ULTRA,
                model_fidelity=fwm.FIDELITY_BALANCED,
                response_method=fwm.RESPONSE_POLE,
                gain_closure_enabled=True, floquet_order=p["floquet_order"],
                velocity_step=1., velocity_cutoff=4., fine_points=0)


def record(tag, raw, elapsed, p, extra=None):
    meta_keys = ["gain_closure", "floquet_convergence", "response_estimator", "claim_gate",
                 "phase_detail", "model_fidelity", "floquet_order", "n_velocity", "N_atoms",
                 "phase_segments", "segment_absorption_od", "pump_depletion_cap",
                 "ultra_pump_remaining_min", "coupling_norm", "transit_reset_rate_rad_s",
                 "ground_collision_dephasing_rate_rad_s", "ultra_spatial_overlap_min"]
    eom = p["opd"] - raw["probe_axis_GHz"]
    item = dict(seconds=elapsed, params=p, eom_GHz=eom, G_s=raw["G_s"], G_c=raw["G_c"],
                G_s_smallsignal=raw["G_s_smallsignal"],
                metadata={k: raw[k] for k in meta_keys if k in raw})
    if extra:
        item["settings_overrides"] = extra
    RESULTS[tag] = item
    assert np.isfinite(raw["G_s"]).all()
    save_json("results.json", RESULTS)
    print(tag, "sec", round(elapsed, 2), "G endpoints", float(raw["G_s"][0]),
          float(raw["G_s"][-1]), "Floquet", raw["floquet_convergence"]["status"], flush=True)
    return raw


def exact(tag, p=None, eom_min=3.030, eom_max=3.039, points=10, **extra):
    p = dict(PARAMS if p is None else p)
    config = kwargs(p)
    config.update(extra)
    t = time.perf_counter()
    raw = fwm.compute_spectrum(p["opd"], scan_min=p["opd"]-eom_max,
                               scan_max=p["opd"]-eom_min, coarse_points=points, **config)
    return record(tag, raw, time.perf_counter()-t, p, extra)


def displayed(tag, **overrides):
    p = dict(PARAMS, **overrides)
    t = time.perf_counter()
    raw = SCHEME.compute(p)
    record(tag, raw, time.perf_counter()-t, p)
    values = np.array([fwm.operating_point(raw, d)["G_s"] for d in DELTA])
    RESULTS[tag]["measured_eom_readout"] = values
    save_json("results.json", RESULTS)
    return raw


def export_csv(name, headers, rows):
    with (OUT/name).open("w", encoding="utf-8-sig", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(headers)
        writer.writerows(rows)


def main():
    snapshot = source_snapshot()
    save_json("source_manifest_before.json", dict(
        head=subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT).decode().strip(),
        status_before=subprocess.check_output(["git", "status", "--porcelain=v1"], cwd=ROOT).decode("utf-8"),
        tracked_existing_sha256=snapshot, defaults=DEFAULTS, actual_params=PARAMS,
        python=sys.version, numpy=np.__version__, scipy=scipy.__version__, kernels=kernels.available()))
    corrected = exact("exact_corrected")
    raw = exact("exact_uncorrected", gain_closure_enabled=False)
    displayed("ui_fast")
    displayed("ui_balanced", resolution=fwm.FIDELITY_BALANCED)
    exact("grid_uncorrected", response_method=fwm.RESPONSE_GRID,
          model_fidelity=fwm.FIDELITY_ULTRA, gain_closure_enabled=False)
    exact("grid_corrected_control", response_method=fwm.RESPONSE_GRID)
    exact("floquet4", floquet_order=4)
    exact("velocity_half_step", velocity_step=.5)
    exact("velocity_5sigma", velocity_step=.5, velocity_cutoff=5.)
    exact("eta_noise_control", detection_efficiency=.50,
          excess_noise_model={"pump_scatter_kappa": .05})
    exact("seed8_default_control", dict(PARAMS, probe_uw=8.0))
    exact("seed4_7_conditional", dict(PARAMS, probe_uw=4.7))
    exact("opd0_9_default_control", dict(PARAMS, opd=.9))
    temperatures = np.arange(113., 123.0001, .5)
    temperature_tags = []
    for temp in temperatures:
        tag = f"temperature_{temp:.2f}"
        exact(tag, dict(PARAMS, temp_c=float(temp)))
        temperature_tags.append(tag)
    for temp in (113., 115., 117., 118., 119., 121., 123.):
        exact(f"peak_{temp:.2f}", dict(PARAMS, temp_c=temp),
              eom_min=2.990, eom_max=3.050, points=241)
    exact("peak_uncorrected", eom_min=2.990, eom_max=3.050, points=241,
          gain_closure_enabled=False)
    exact("peak_opd0_9", dict(PARAMS, opd=.9),
          eom_min=2.990, eom_max=3.050, points=241)
    summary = {tag: fits(RESULTS[tag]["G_s"]) for tag in
               ["exact_corrected", "exact_uncorrected", "opd0_9_default_control", *temperature_tags]}
    for tag in ("ui_fast", "ui_balanced"):
        summary[tag] = fits(RESULTS[tag]["measured_eom_readout"])
    base = RESULTS["exact_corrected"]["G_s"]
    checks = {}
    for tag, compare in (("grid_uncorrected", RESULTS["exact_uncorrected"]["G_s"]),
                         ("grid_corrected_control", base), ("floquet4", base),
                         ("velocity_half_step", base), ("velocity_5sigma", base),
                         ("eta_noise_control", base), ("seed8_default_control", base),
                         ("seed4_7_conditional", base)):
        values = RESULTS[tag]["G_s"]
        checks[tag] = dict(max_relative_pct=float(100*np.max(np.abs(values/compare-1))),
                           gain=values)
    for tag in ("ui_fast", "ui_balanced"):
        v = RESULTS[tag]["measured_eom_readout"]
        checks[tag] = dict(max_relative_pct=float(100*np.max(np.abs(v/base-1))), gain=v)
    peaks = {}
    for tag, value in RESULTS.items():
        if tag.startswith("peak_"):
            idx = int(np.argmax(value["G_s"]))
            peaks[tag] = dict(eom_GHz=float(value["eom_GHz"][idx]),
                              gain=float(value["G_s"][idx]), edge=idx in (0,len(value["G_s"])-1),
                              grid_step_MHz=.25)
    best = min(temperature_tags, key=lambda t: summary[t]["identity"]["rmse"])
    shape_best = min(temperature_tags, key=lambda t: summary[t]["scale"]["rmse"])
    save_json("summary.json", dict(fits=summary, numerical_checks=checks, peaks=peaks,
                                   common_temperature_min_rmse_tag=best,
                                   common_temperature_with_scale_min_rmse_tag=shape_best,
                                   temperature_density_m3={str(t): hyperfine.number_density(t+273.15)
                                                           for t in (113,115,117,118,119,121,123)}))
    gs, gc = corrected["G_s"], corrected["G_c"]
    export_csv("comparison.csv", ["EOM_GHz", "TPD_GABES_MHz", "seed_uW", "probe_uW", "conjugate_uW",
               "gain_measured", "gain_exact_corrected", "conjugate_gain_exact_corrected",
               "gain_exact_uncorrected", "gain_UI_Fast", "gain_UI_Balanced",
               "residual_measured_minus_sim", "simulation_minus_measured_pct", "measured_over_sim"],
               zip(EOM, DELTA, SEED, PROBE, CONJ, MEAS, gs, gc, raw["G_s"],
                   RESULTS["ui_fast"]["measured_eom_readout"],
                   RESULTS["ui_balanced"]["measured_eom_readout"], MEAS-gs, 100*(gs-MEAS)/MEAS, MEAS/gs))
    export_csv("temperature_fits.csv", ["temperature_C", "raw_RMSE", "scale_a", "scale_RMSE",
               "affine_a", "affine_b", "affine_RMSE", "start_gain", "end_gain"],
               ((RESULTS[t]["params"]["temp_c"], summary[t]["identity"]["rmse"],
                 summary[t]["scale"]["a"], summary[t]["scale"]["rmse"],
                 summary[t]["affine"]["a"], summary[t]["affine"]["b"], summary[t]["affine"]["rmse"],
                 RESULTS[t]["G_s"][0], RESULTS[t]["G_s"][-1]) for t in temperature_tags))
    plt.rcParams.update({"font.family":"DejaVu Sans", "font.size":10, "axes.grid":True,
                         "grid.alpha":.25, "savefig.dpi":180})
    x = EOM
    fig, ax = plt.subplots(2, 2, figsize=(12, 8), layout="constrained")
    ax[0,0].plot(x, MEAS, "ko-", label="Measured Probe / pump-off Seed")
    ax[0,0].plot(x, gs, "o-", label="Corrected, exact EOM rows")
    ax[0,0].plot(x, raw["G_s"], "o-", label="Correction off")
    ax[0,0].set_yscale("log"); ax[0,0].set_ylabel("Probe power gain")
    ax[0,0].legend(fontsize=8)
    ax[0,1].plot(x, MEAS/gs,"o-",label="Measured / corrected simulation")
    ax[0,1].set_ylabel("Gain ratio"); ax[0,1].legend(fontsize=8)
    for y, label in ((MEAS,"Measured"), (gs,"Corrected"), (raw["G_s"],"Uncorrected")):
        ax[1,0].plot(x,y/y[-1],"o-",label=label)
    ax[1,0].set_ylabel("Gain / gain at 3.030 GHz"); ax[1,0].legend(fontsize=8)
    for name in ("identity","scale","affine"):
        ax[1,1].plot(x,summary["exact_corrected"][name]["residual_meas_minus_prediction"],"o-",label=name)
    ax[1,1].axhline(0,color="k",lw=.8); ax[1,1].set_ylabel("Measured - predicted gain")
    ax[1,1].legend(fontsize=8)
    for a in ax.flat:
        a.set_xlabel("EOM frequency (GHz)"); a.set_xlim(3.0395,3.0295)
    fig.suptitle("118 C, 380 mW; OPD 1.06 GHz; input seed 3.7 uW; other defaults")
    fig.savefig(OUT/"gain_comparison.png"); plt.close(fig)
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.6), layout="constrained")
    for temp in (113.,115.,117.,118.,119.,121.,123.):
        tag=f"temperature_{temp:.2f}"
        g=RESULTS[tag]["G_s"]
        ax[0].plot(x,g,label=f"{temp:g} C")
        ax[1].plot(x,g/g[-1],label=f"{temp:g} C")
    ax[0].plot(x,MEAS,"ko--",label="Measured"); ax[1].plot(x,MEAS/MEAS[-1],"ko--",label="Measured")
    for a in ax[:2]:
        a.set_xlabel("EOM frequency (GHz)"); a.set_xlim(3.0395,3.0295); a.legend(fontsize=7)
    ax[0].set_ylabel("Probe gain"); ax[1].set_ylabel("Gain / gain at 3.030 GHz")
    ax[2].plot(temperatures,[summary[t]["identity"]["rmse"] for t in temperature_tags],label="Raw RMSE")
    ax[2].plot(temperatures,[summary[t]["scale"]["rmse"] for t in temperature_tags],label="After common scale")
    ax[2].plot(temperatures,[summary[t]["affine"]["rmse"] for t in temperature_tags],label="After affine")
    ax[2].set_xlabel("Common temperature (C)"); ax[2].set_ylabel("Gain RMSE"); ax[2].legend(fontsize=8)
    fig.savefig(OUT/"temperature_sensitivity.png"); plt.close(fig)
    fig, ax = plt.subplots(figsize=(9,4.8),layout="constrained")
    for tag,label in (("peak_113.00","113 C corrected"),("peak_118.00","118 C corrected"),
                      ("peak_123.00","123 C corrected"),("peak_uncorrected","118 C correction off"),
                      ("peak_opd0_9","118 C OPD 0.9 GHz control")):
        v=RESULTS[tag]; y=np.asarray(v["G_s"])
        ax.plot(v["eom_GHz"],y/y.max(),label=label)
    ax.plot(x,MEAS/MEAS.max(),"ko",label="Measured / measured maximum in table")
    ax.axvline(3.030,color="k",linestyle=":",label="User-reported peak vicinity")
    ax.set_xlim(3.050,2.990); ax.set_xlabel("EOM frequency (GHz)")
    ax.set_ylabel("Gain / maximum in each scanned interval"); ax.legend(fontsize=8)
    fig.savefig(OUT/"peak_sensitivity.png"); plt.close(fig)
    after = source_snapshot()
    changed = [n for n,h in snapshot.items() if after.get(n)!=h]
    save_json("source_preservation_check.json",dict(changed_existing_files=changed,
              existing_source_files_checked=len(snapshot), no_existing_source_changes=not changed))
    if changed:
        raise RuntimeError(f"Existing source changed during run: {changed}")
    print("Finished; source preservation checked", len(snapshot), "files", flush=True)


if __name__ == "__main__":
    main()
