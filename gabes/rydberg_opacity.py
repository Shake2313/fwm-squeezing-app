"""Declared effective opacity for the compact Rydberg propagation model.

This is a phenomenological participation multiplier, separate from the atomic
dipole/line strength and from the atom-number bookkeeping fraction.  No measured
opacity uncertainty or held-out calibration is available.  Scenario endpoints
are sensitivity checks, not confidence limits or standard uncertainties.
"""

import math


DEFAULT_EFFECTIVE_OPACITY = 0.001
DEFAULT_RELATIVE_UNCERTAINTY = 0.5


def opacity_spec(params):
    """Validate live inputs and return their provenance and scenario bounds."""
    value = float(params.get(
        "effective_opacity_participation", DEFAULT_EFFECTIVE_OPACITY))
    relative = float(params.get(
        "effective_opacity_relative_uncertainty", DEFAULT_RELATIVE_UNCERTAINTY))
    if not math.isfinite(value) or not 0.0 <= value <= 1.0:
        raise ValueError("effective_opacity_participation must be finite in [0, 1]")
    if not math.isfinite(relative) or not 0.0 <= relative <= 1.0:
        raise ValueError(
            "effective_opacity_relative_uncertainty must be finite in [0, 1]")
    lower = max(0.0, value * (1.0 - relative))
    upper = min(1.0, value * (1.0 + relative))
    return {
        "value": value,
        "unit": "1",
        "bounds": [0.0, 1.0],
        "relative_uncertainty": relative,
        "measured_standard_uncertainty": None,
        "scenario_values": [lower, value, upper],
        "uncertainty_kind": "assumed_relative_scenario_half_width_not_standard_uncertainty",
        "status": "ASSUMED",
        "calibration_status": "UNCALIBRATED",
        "prediction_status": "PREDICTED",
        "source": "Legacy compact-model opacity 0.001; user-variable assumption, not measured participation",
        "calibration_dependency": "rydberg-experimental-validation-dataset",
        "formula": "chi_phys = atomic_chi_phys * effective_opacity_participation; T = exp(-k Im(chi_phys) L)",
        "scope": "Optical susceptibility/finite-IF propagation only; effective atom number is independent",
    }


def opacity_samples(params):
    """Unique lower/nominal/upper scenarios, retaining nominal at boundaries."""
    spec = opacity_spec(params)
    samples = [{"label": "nominal", "value": spec["value"]}]
    for label, value in zip(("lower", "upper"),
                            (spec["scenario_values"][0], spec["scenario_values"][2])):
        if not any(sample["value"] == value for sample in samples):
            samples.append({"label": label, "value": value})
    return samples
