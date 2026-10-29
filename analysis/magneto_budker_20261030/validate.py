"""Historical bounded Budker adapter checks through the unregistered single-pair host.

This reproduces the initial opt-in architecture, not the active unified scheme.

Run: python -m analysis.magneto_budker_20261030.validate --line D1
"""
import argparse
import hashlib
import json
from pathlib import Path
import time

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from gabes.schemes.magneto import SinglePairMagnetoReference
from gabes.schemes import magneto_budker
from gabes.magneto_multilevel import multilevel_response
from gabes.magneto_noise import magnetic_sensitivity
from gabes_ui.export import to_json_primitive


def run(line="D1", output=None):
    out = Path(output) if output else Path(__file__).parent / "results"
    out.mkdir(parents=True, exist_ok=True)
    scheme = SinglePairMagnetoReference()
    preset = next(p for p in scheme.presets() if p.name.startswith(f"Budker {line}"))
    params = dict(scheme.defaults(), **preset.values)
    cases = {
        "reference": params,
        "refined": dict(params, budker_detuning_points=49,
                        budker_velocity_points=10, budker_local_step_ut=5e-6),
        "velocity_preserved": dict(params, budker_velocity_mixing="none"),
    }
    results, raws = {}, {}
    for name, settings in cases.items():
        started = time.perf_counter()
        raw = scheme.compute(settings)
        view = scheme.headless_observables(raw, settings)
        results[name] = dict(seconds=time.perf_counter()-started,
                             params=settings, readout=view["magneto_readout"])
        raws[name] = raw
        print(f"{line} {name}: {results[name]['seconds']:.2f}s, "
              f"minimum {view['magneto_readout']['optimum_detuning_mhz']} MHz", flush=True)
    base, fine = (results[name]["readout"] for name in ("reference", "refined"))
    x, xf = np.asarray(base["detunings_mhz"]), np.asarray(fine["detunings_mhz"])
    indices = np.array([int(np.flatnonzero(xf == value)[0]) for value in x])
    coarse_slope = np.asarray(base["slope_rad_per_t"], float)
    fine_slope = np.asarray(fine["slope_rad_per_t"], float)[indices]
    shared_valid = (np.asarray(base["slope_valid"], bool)
                    & np.asarray(fine["slope_valid"], bool)[indices])
    slope_scale = np.max(np.abs(fine_slope[shared_valid]))
    shared_error = np.max(np.abs(coarse_slope[shared_valid]-fine_slope[shared_valid]))/slope_scale
    # Endpoint intensity scenarios assess attenuation sensitivity only. They
    # do not establish rigorous propagation error bounds or a spatial solution.
    best = fine["optimum_index"]
    selected = np.array([0.0, fine["optimum_detuning_mhz"]])
    attenuation_factor = float(np.min(fine["transmission"]))
    conditions = dict(raws["refined"]["conditions"])
    conditions["intensity_mw_cm2"] *= attenuation_factor
    attenuation = multilevel_response(selected, b_ut=raws["refined"]["b_ut"], **conditions)
    h_t = raws["refined"]["local_step_ut"] * 1e-6
    slope = (attenuation["rotation_rad"][3]-attenuation["rotation_rad"][1])/h_t
    budget = magnetic_sensitivity(
        attenuation["transmission"][2], slope,
        frequency_hz=attenuation["frequency_hz"],
        # Compare the atomic-intensity assumption at the same photon-budget
        # input power; transmission is supplied by each scenario separately.
        intensity_mw_cm2=params["intensity_mw_cm2"],
        beam_diameter_mm=params["beam_diameter_mm"],
        efficiency=params["detection_efficiency"], bandwidth_hz=params["noise_bandwidth_hz"],
        signal_rms_t=params["magnetic_signal_ft"]*1e-15)
    assessment = dict(
        status="conditional model and numerical checks; not quantitative paper reproduction",
        shared_detuning_samples_compared=int(np.count_nonzero(shared_valid)),
        shared_detuning_slope_max_error_over_peak=float(shared_error),
        optimum_shift_mhz=fine["optimum_detuning_mhz"]-base["optimum_detuning_mhz"],
        optimum_asd_relative_change=(fine["optimum_field_asd_t_per_sqrt_hz"]
                                    /base["optimum_field_asd_t_per_sqrt_hz"]-1),
        off_centroid_minimum=bool(best is not None and not fine["optimum_at_scan_edge"]
                                  and abs(fine["optimum_detuning_mhz"]) > 0),
        centroid_snr_improvement=fine["snr_improvement_over_centroid"],
        attenuation_scenario=dict(
            scope="uniform atomic intensity lowered by minimum T; sensitivity check, not a spatial solution or error bound",
            intensity_factor=attenuation_factor, detunings_mhz=selected,
            field_asd_t_per_sqrt_hz=budget["field_asd_t_per_sqrt_hz"],
            nominal_field_asd_t_per_sqrt_hz=[fine["centroid_field_asd_t_per_sqrt_hz"],
                                           fine["optimum_field_asd_t_per_sqrt_hz"]]),
    )
    render = scheme.observables(raws["refined"], cases["refined"])
    render["figure"].savefig(out/f"{line.lower()}_sensitivity.png", dpi=155)
    plt.close(render["figure"])
    root = Path(__file__).resolve().parents[2]
    files = ["gabes/magneto_multilevel.py", "gabes/magneto_noise.py",
             "gabes/schemes/magneto_budker.py", "gabes/schemes/magneto.py"]
    payload = dict(line=line, reference="10.1103/PhysRevA.62.043403",
                   source_sha256={name: hashlib.sha256((root/name).read_bytes()).hexdigest()
                                  for name in files}, cases=results, assessment=assessment)
    (out/f"{line.lower()}_validation.json").write_text(
        json.dumps(to_json_primitive(payload), indent=2, ensure_ascii=False, allow_nan=False),
        encoding="utf-8")
    print(json.dumps(to_json_primitive(assessment), indent=2, ensure_ascii=False), flush=True)
    return payload


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--line", choices=["D1", "D2"], default="D1")
    parser.add_argument("--output")
    args = parser.parse_args()
    run(args.line, args.output)
