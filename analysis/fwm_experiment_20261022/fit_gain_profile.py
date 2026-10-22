"""Reproduce the practical 10/13 gain calibration, independently of GC.

Run from the repository root with ``python analysis/fwm_experiment_20261022/fit_gain_profile.py``.
The JSON is a development fit, not an untouched prediction validation.
"""
from pathlib import Path
import hashlib
import json
import sys

import numpy as np
from scipy.interpolate import PchipInterpolator
from scipy.optimize import least_squares

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from gabes import fwm_measured_gain as measured
from gabes.schemes.fwm import FIDELITY_BALANCED, FWMScheme


def reference_spectrum():
    context = measured.nominal_context()
    scheme = FWMScheme()
    params = dict(
        scheme.defaults(), opd=context["D_GHz"], temp_c=context["T"]-273.15,
        pump_mw=context["P_pump"]*1e3, probe_uw=context["P_probe"]*1e6,
        cell_mm=context["L"]*1e3, pump_waist_um=context["w_pump"]*1e6,
        probe_waist_um=context["w_probe"]*1e6,
        seeded_angle_deg=context["pump_probe_angle_deg"], resolution=FIDELITY_BALANCED,
    )
    return scheme.compute(params)


def main():
    raw = reference_spectrum()
    delta = (raw["probe_axis_GHz"]-raw["raman_center_minus_GHz"])*1e3
    source = measured.observations()
    nodes = measured.internal_delta_mhz_from_eom(source["eom_MHz"])
    result = {
        "model_id": "fwm_20261013_warp_affine_residual_v1",
        "measurement_date": "2026-10-13",
        "source_workbook_sha256": measured.SOURCE_WORKBOOK_SHA256,
        "nominal_context": measured.nominal_context(),
        "source_internal_delta_mhz": nodes.tolist(),
        "reference_model_internal_delta_mhz": delta.tolist(),
        "reference_model_G_s": raw["G_s"].tolist(),
        "reference_model_G_c": raw["G_c"].tolist(),
        "arms": {},
    }
    for arm, key in (("probe", "G_s"), ("conjugate", "G_c")):
        data = source[arm+"_ratio"]
        model = raw[key]
        interp = PchipInterpolator(delta, np.log(np.maximum(model, 1e-30)), extrapolate=False)
        target_peak = float(np.mean(nodes[data >= np.max(data)*(1-1e-12)]))
        model_peak = float(delta[np.argmax(model)])

        def prediction(parameters):
            shift, width, log_scale, offset = parameters
            positions = model_peak+shift+width*(nodes-target_peak)
            return np.exp(log_scale+interp(positions))+offset

        def residual(parameters):
            return np.log(prediction(parameters))-np.log(data)

        fit = least_squares(
            residual, [0., 1., np.log(np.max(data)/np.max(model)), 0.001],
            bounds=([-8., .25, -12., 0.], [8., 4., 3., 1.]),
            xtol=1e-12, ftol=1e-12, gtol=1e-12, max_nfev=1500,
        )
        shift, width, log_scale, offset = fit.x
        base = prediction(fit.x)
        corrections = np.log(data/base)
        positions = model_peak+shift+width*(nodes-target_peak)
        nominal_warped = np.exp(interp(positions))
        result["arms"][arm] = {
            "model_peak_delta_mhz": model_peak,
            "measured_peak_delta_mhz": target_peak,
            "source_center_delta_mhz": float(model_peak+shift),
            "frequency_scale": float(width),
            "gain_scale": float(np.exp(log_scale)),
            "gain_offset": float(offset),
            "log_residual_at_nodes": corrections.tolist(),
            "log_total_coefficient_at_nodes": np.log(data/nominal_warped).tolist(),
            "measured_ratios": data.tolist(),
            "affine_max_relative_error": float(np.max(np.abs(base/data-1))),
            "affine_rms_log_error": float(np.sqrt(np.mean(np.log(base/data)**2))),
            "fit_success": bool(fit.success),
        }
    result["runtime_formula"] = (
        "x_j=c_j+w_j*(delta-d_j); "
        "R_j=G_j(current,x_j)*exp(PCHIP(delta_i,log(R_measured_i/G_j(nominal,x_j_i)))(delta)); "
        "outside measured support, total correction coefficient holds its boundary value while G_j continues. "
        "This carries the fitted affine offset proportionally with current model intensity, so zero source remains zero."
    )
    result["calibration_validation_scope"] = "all 23 observations used; exact-node residual fit; not independent validation"
    result["independent_validation"] = False
    result["runtime_context_gate"] = False
    result["source_file_sha256"] = {
        name: hashlib.sha256((ROOT/name).read_bytes()).hexdigest()
        for name in ("gabes/fwm_measured_gain.py", "gabes/observables.py", "gabes/schemes/fwm.py")
    }
    path = Path(__file__).with_name("gain_profile.json")
    path.write_text(json.dumps(result, ensure_ascii=False, indent=2)+"\n", encoding="utf-8")
    print(json.dumps({arm: {key: fit[key] for key in (
        "source_center_delta_mhz", "frequency_scale", "gain_scale", "gain_offset",
        "affine_max_relative_error", "fit_success")}
        for arm, fit in result["arms"].items()}, indent=2))


if __name__ == "__main__":
    main()
