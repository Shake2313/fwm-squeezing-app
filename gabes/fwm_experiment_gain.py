"""Practical 10/13/2026 FWM gain calibration, separate from atomic/GC solves.

The measured pump-on / pump-off detector ratios calibrate a frequency-warped
model. Frozen positive correction factors act on each CURRENT model spectrum;
they do not replace it with a measured lookup curve. Outside measured support,
only the correction coefficient holds its boundary value: the model continues.
"""
import copy

import numpy as np
from scipy.interpolate import PchipInterpolator

from . import fwm_measured_gain as measured


MODEL_ID = "fwm_20261013_warp_affine_residual_v1"
# Frozen offline fit, reproduced by analysis/fwm_experiment_20261022/fit_gain_profile.py.
# BEGIN FITTED COEFFICIENTS
_ARMS = {'probe': {'model_peak_delta_mhz': 24.749999999999606,
           'measured_peak_delta_mhz': 2.2324390000003405,
           'source_center_delta_mhz': 22.827113312494454,
           'frequency_scale': 1.2392520825296043,
           'gain_scale': 0.09401476007742224,
           'gain_offset': 0.015291576826316149,
           'affine_max_relative_error': 0.25546760317046724,
           'affine_rms_log_error': 0.12136298700494211,
           'log_total_coefficient_at_nodes': [-2.517957851200495,
                                              -2.503995545291295,
                                              -2.5137673070531177,
                                              -2.4485481032565097,
                                              -2.375518819648362,
                                              -2.3266045808010096,
                                              -2.33077689206703,
                                              -2.2870940788820855,
                                              -2.25240072839961,
                                              -2.2386856298749525,
                                              -2.192754219042939,
                                              -2.2193116033664104,
                                              -2.2897526825171437,
                                              -2.3567679242015243,
                                              -2.395190087235892,
                                              -2.4574839135036313,
                                              -2.4650009808774835,
                                              -2.503097335818218,
                                              -2.181096369153177,
                                              -2.323827115233929,
                                              -1.748098562172115,
                                              -1.5595865045394315,
                                              -1.907233438221969]},
 'conjugate': {'model_peak_delta_mhz': 24.749999999999606,
               'measured_peak_delta_mhz': 2.7324390000003405,
               'source_center_delta_mhz': 22.689633036415053,
               'frequency_scale': 1.294620567486792,
               'gain_scale': 0.09409948262794274,
               'gain_offset': 0.17628769405430758,
               'affine_max_relative_error': 0.2876951260476732,
               'affine_rms_log_error': 0.14779889425949747,
               'log_total_coefficient_at_nodes': [-2.566744988284615,
                                                  -2.5837915556673354,
                                                  -2.5077399275123704,
                                                  -2.4262740624082477,
                                                  -2.334979254310135,
                                                  -2.248954102201915,
                                                  -2.2745780939838216,
                                                  -2.227927558690042,
                                                  -2.198546198566624,
                                                  -2.17171466872008,
                                                  -2.198459747510887,
                                                  -2.2051736208963955,
                                                  -2.2670844254470475,
                                                  -2.3737022816810356,
                                                  -2.4278681759505734,
                                                  -2.4995453484526693,
                                                  -2.504238713367706,
                                                  -2.4762389345054436,
                                                  -1.9858426601340071,
                                                  -1.7040635113291327,
                                                  -1.2348914272121547,
                                                  -0.9499447191181953,
                                                  -1.0393488456724986]}}
# END FITTED COEFFICIENTS

_OBSERVATIONS = measured.observations()
_NODES = measured.internal_delta_mhz_from_eom(_OBSERVATIONS["eom_MHz"])
_NOMINAL = measured.nominal_context()
_GAIN_QUANTITY = "pump-on arm output / measured pump-off seed detector output"


def profile():
    """Fresh calibration metadata; empirical fit, not independent validation."""
    return {
        "model_id": MODEL_ID,
        "measurement_date": "2026-10-13",
        "source_workbook_sha256": measured.SOURCE_WORKBOOK_SHA256,
        "calibration_record": "analysis/fwm_experiment_20261022/gain_profile.json",
        "nominal_context": copy.deepcopy(_NOMINAL),
        "eom_support_MHz": list(measured.EOM_SUPPORT_MHZ),
        "internal_delta_support_mhz": [float(_NODES[0]), float(_NODES[-1])],
        "gain_quantity": _GAIN_QUANTITY,
        "gain_reference_plane": "pump-off seed detector output, common measured voltage conversion",
        "method": "per-arm frequency shift/scale, positive affine gain fit and bounded residual correction",
        "formula": "x_j=c_j+w_j*(delta-d_j); R_j=G_j(current,x_j)*C_j(delta)",
        "coefficient_formula": "C_j(delta_i)=measured_ratio_j_i/G_j(nominal,x_j_i)",
        "residual_interpolation": "shape-preserving cubic in log coefficient, within measured support",
        "outside_measurement": "hold boundary correction coefficient; continue current warped model",
        "outside_atomic_scan": "positive continuation to unseeded vacuum/conjugate=0 and probe=1",
        "context_policy": "changed operating conditions retain calibration and are marked as extrapolation",
        "independent_validation": False,
        "node_agreement": "23-point development fit by construction, not a predictive validation result",
        "canonical_transfer_modified": False,
        "noise_covariance_calibrated": False,
        "arms": copy.deepcopy(_ARMS),
    }


def _axis(raw):
    branch = raw.get("branch", -1)
    center_key = "raman_center_plus_GHz" if branch == 1 else "raman_center_minus_GHz"
    center = raw.get(center_key, raw.get("raman_center_minus_GHz"))
    if center is None:
        raise ValueError("calibration requires the Raman scan center")
    axis = (np.asarray(raw["probe_axis_GHz"], dtype=float)-float(center))*1e3
    if axis.ndim != 1 or axis.size < 2 or np.any(~np.isfinite(axis)):
        raise ValueError("calibration requires at least two finite scan points")
    return axis


def _sample_model(axis, gain, positions, *, arm):
    """Positive interpolation with a finite explicit continuation past the scan.

    Log interpolation avoids negative cubic overshoot. Beyond the computed scan,
    a positive exponential continuation approaches the passive asymptote over
    one fifth of the computed span. That scope is reported separately from the
    much narrower experimental calibration interval.
    """
    gain = np.asarray(gain, dtype=float)
    positions = np.asarray(positions, dtype=float)
    if gain.shape != axis.shape or np.any(~np.isfinite(gain)) or np.any(gain < 0):
        raise ValueError("calibration gain arrays must be finite, nonnegative and match scan points")
    order = np.argsort(axis)
    x, g = axis[order], gain[order]
    if np.any(np.diff(x) <= 0):
        raise ValueError("calibration scan points must be distinct")
    outside = (positions < x[0]) | (positions > x[-1])
    if not np.any(g):
        return np.zeros_like(positions), outside
    bounded = np.clip(positions, x[0], x[-1])
    sampled = np.exp(PchipInterpolator(x, np.log(np.maximum(g, 1e-30)), extrapolate=False)(bounded))
    asymptote = 1.0 if arm == "probe" else 0.0
    length = max(float((x[-1]-x[0])/5.0), 1.0)
    distance = np.abs(positions-bounded)
    sampled = np.where(outside, asymptote+(sampled-asymptote)*np.exp(-distance/length), sampled)
    return sampled, outside


def _context_changes(raw):
    keys = {
        "D_GHz": "D_GHz", "T": "temperature_K", "P_pump": "input_pump_power_W",
        "P_probe": "input_seed_power_W", "w_pump": "w_pump_m", "w_probe": "w_probe_m",
        "L": "cell_length_m", "pump_probe_angle_deg": "pump_probe_angle_deg", "branch": "branch",
    }
    changes = []
    for nominal_key, raw_key in keys.items():
        value = raw.get(raw_key)
        try:
            matches = bool(np.isfinite(float(value)) and np.isclose(
                float(value), float(_NOMINAL[nominal_key]), rtol=1e-6, atol=1e-12))
        except (TypeError, ValueError):
            matches = False
        if not matches:
            changes.append(raw_key)
    return changes


def evaluate(raw, delta_mhz, *, enabled=True):
    """Evaluate corrected gains at arbitrary internal detunings, including tails.

    Detector/noise sliders are absent from this API. Its frozen correction acts
    on a newly computed model when physical controls change, without refitting.
    """
    axis = _axis(raw)
    delta = np.asarray(delta_mhz, dtype=float)
    if np.any(~np.isfinite(delta)):
        raise ValueError("calibration detunings must be finite")
    tolerance = 64*np.finfo(float).eps*max(1.0, float(np.max(np.abs(_NODES))))
    within = (delta >= _NODES[0]-tolerance) & (delta <= _NODES[-1]+tolerance)
    changed = _context_changes(raw)
    result = {
        "applied": bool(enabled),
        "internal_delta_mhz": delta.copy(),
        "eom_MHz": measured.eom_mhz_from_internal_delta(delta),
        "within_measured_support": within,
        "context_changed": bool(changed),
        "changed_context_inputs": changed,
        "status": ("DISABLED" if not enabled else
                   "EMPIRICAL_CALIBRATION_CHANGED_CONDITIONS" if changed else "EMPIRICAL_CALIBRATION"),
        "point_status": np.where(within, "CALIBRATED_MEASUREMENT_REGION", "EXTRAPOLATED_FREQUENCY"),
        "source_delta_mhz": {}, "outside_atomic_scan": {},
        "profile": profile(),
        "gain_quantity": _GAIN_QUANTITY,
    }
    for arm, key in (("probe", "G_s"), ("conjugate", "G_c")):
        fit = _ARMS[arm]
        positions = (fit["source_center_delta_mhz"]
                     +fit["frequency_scale"]*(delta-fit["measured_peak_delta_mhz"])) if enabled else delta
        sample, outside = _sample_model(axis, raw[key], positions, arm=arm)
        if enabled:
            bounded = np.clip(delta, _NODES[0], _NODES[-1])
            coefficient = np.exp(PchipInterpolator(
                _NODES, fit["log_total_coefficient_at_nodes"], extrapolate=False)(bounded))
            sample = sample*coefficient
        result[key] = sample
        result["source_delta_mhz"][arm] = positions
        result["outside_atomic_scan"][arm] = outside
    return result


def calibrate_spectrum(raw, *, enabled=False, selected_delta_mhz=None):
    """Return a nonmutating calibrated view on the cached spectrum's own axis."""
    axis = _axis(raw)
    result = evaluate(raw, axis, enabled=enabled)
    if not enabled:
        # Preserve the exact uncalibrated result; log interpolation at its knots
        # would otherwise introduce unnecessary roundoff in a disabled toggle.
        result["G_s"] = np.array(raw["G_s"], copy=True)
        result["G_c"] = np.array(raw["G_c"], copy=True)
    if selected_delta_mhz is not None:
        selected = evaluate(raw, float(selected_delta_mhz), enabled=enabled)
        result["selected_point"] = {
            "G_s": float(selected["G_s"]), "G_c": float(selected["G_c"]),
            "available": True, "applied": bool(enabled),
            "internal_delta_mhz": float(selected_delta_mhz),
            "eom_MHz": float(selected["eom_MHz"]),
            "within_measured_support": bool(selected["within_measured_support"]),
            "context_changed": bool(selected["context_changed"]),
            "status": str(selected["point_status"])+( "; CHANGED_OPERATING_CONDITIONS" if selected["context_changed"] else ""),
        }
    return result
