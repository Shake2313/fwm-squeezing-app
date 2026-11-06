"""Native PSR UI adapter for the longitudinal pump/noise solver."""
from copy import deepcopy
from functools import lru_cache

import numpy as np

from ..psr_conventions import jones_from_angles, jones_polarization_angles
from ..psr_propagation import MODEL_ID, solve_psr_propagation


_LENGTH_SETTINGS = {
    "Quick": dict(segments=1, max_segments=8, rtol=.02, max_atomic_evaluations=12000),
    "Balanced": dict(segments=2, max_segments=32, rtol=.005, max_atomic_evaluations=30000),
    "Fine": dict(segments=4, max_segments=64, rtol=.001, max_atomic_evaluations=80000),
}


@lru_cache(maxsize=2)
def _propagated(inputs, rf_hz, length_m, accuracy, thermal_settings):
    detuning, intensity, azimuth, ellipticity, velocity_model, velocity, temperature, reset = inputs
    if accuracy not in _LENGTH_SETTINGS:
        raise ValueError("unknown length integration accuracy")
    return solve_psr_propagation(
        detuning, rf_hz, length_m, intensity_mw_cm2=intensity,
        jones_xy=jones_from_angles(np.deg2rad(azimuth), np.deg2rad(ellipticity)),
        velocity_model=velocity_model, velocity_m_s=velocity,
        temp_c=temperature, ground_relaxation_hz=reset,
        thermal_options=dict(thermal_settings), **_LENGTH_SETTINGS[accuracy])


def compute_propagated(params, rf_ratios, thermal_settings):
    center_hz = float(params["rf_frequency_mhz"])*1e6
    length = float(params["cell_mm"])*.001
    if not np.isfinite(center_hz) or center_hz <= 0 or not np.isfinite(length) or length < 0:
        raise ValueError("positive analysis frequency and nonnegative cell length required")
    velocity_model = params["prop_velocity_model"]
    thermal = velocity_model == "thermal"
    raw = dict(kind="psr_rf", model_id=MODEL_ID, center_hz=center_hz, cell_length_m=length,
               unavailable_reason=None, velocity_model=velocity_model,
               propagation_convergence=dict(passed=False, status="propagation incomplete"),
               spectral_observable="Narrowband PSD in the polarization orthogonal to the output carrier",
               local_spatial_frame="Cartesian optical envelope in the local mean common-phase gauge",
               transfer_spatial_frame="Accumulated mean common-phase gauge; output Jones uses this same frame")
    if thermal:
        raw["thermal_convergence"] = dict(passed=False, status="local Maxwell convergence not evaluated")
    inputs = (float(params["obe_detuning_mhz"]), float(params["obe_intensity_mw_cm2"]),
              float(params["obe_azimuth_deg"]), float(params["obe_ellipticity_deg"]), velocity_model,
              float(params["prop_velocity_m_s"]) if not thermal else 0.,
              float(params["temp_c"]), float(params["relaxation_hz"]))
    if inputs[1] == 0:
        raw["unavailable_reason"] = "A nonzero bright input carrier is required for Stokes normalization."
        return raw
    try:
        positive = center_hz*np.asarray(rf_ratios)
        settings = tuple(sorted(thermal_settings[params["prop_thermal_accuracy"]].items())) if thermal else ()
        result = _propagated(inputs, tuple(np.r_[-positive[::-1], positive]), length,
                             params["propagation_accuracy"], settings)
        raw.update(conditions=deepcopy(result["conditions"]),
                   propagation_convergence=deepcopy(result["convergence"]),
                   propagation_performance=deepcopy(result["performance"]),
                   propagation_velocity_convergence=deepcopy(result.get("velocity_convergence")),
                   pump_profile=deepcopy(result["profile"]),
                   segment_diagnostics=deepcopy(result["segment_diagnostics"]),
                   pump_diagnostics=deepcopy(result["pump_diagnostics"]),
                   output_intensity_mw_cm2=result["output_intensity_mw_cm2"],
                   output_jones_xy=result["output_jones_xy"].copy(),
                   carrier_phase_advance_rad=result["carrier_phase_advance_rad"],
                   completed_segments=result["segments"])
        if thermal:
            local_checks = [entry["velocity_convergence"] for entry in result["segment_diagnostics"]]
            passed = bool(local_checks and all(item["converged"] for item in local_checks))
            if length == 0:
                passed = True
            raw["thermal_convergence"] = dict(passed=passed,
                status="Local Maxwell kernels and cutoff converged" if passed else "Local Maxwell checks incomplete",
                scope="Local velocity integrals along the computed pump trajectory",
                local_checks=local_checks)
        if not result["convergence"]["passed"]:
            raw["unavailable_reason"] = "Longitudinal propagation is not converged. "+result["convergence"]["status"]
            return raw
        if result["output_intensity_mw_cm2"] <= np.finfo(float).tiny:
            raw["unavailable_reason"] = "The output carrier is unresolved for Stokes normalization."
            return raw
        transfer = result["transfer"]
        raw.update(omega_rad_s=transfer.frequency_axis.omega_rad_s,
                   frequency_frame=transfer.frequency_axis.frame, mode_labels=transfer.mode_labels,
                   signs=transfer.signs, transfer=transfer.transfer,
                   added_greater=transfer.noise_greater, added_lesser=transfer.noise_lesser,
                   transfer_audit=transfer.audit())
    except (ValueError, np.linalg.LinAlgError, FloatingPointError) as error:
        raw["unavailable_reason"] = f"Longitudinal calculation unavailable. {error}"
    return raw


def profile_readout(raw):
    profile = raw.get("pump_profile")
    if profile is None:
        return None
    angles = np.array([jones_polarization_angles(j) for j in profile["jones_xy"]])
    # Azimuth is pi-periodic. Unwrap only finite contiguous sections, preserving
    # an explicitly undefined azimuth at circular polarization.
    azimuth = angles[:, 0].copy()
    valid = np.flatnonzero(np.isfinite(azimuth))
    for indices in np.split(valid, np.flatnonzero(np.diff(valid) > 1)+1):
        if len(indices):
            azimuth[indices] = .5*np.unwrap(2*azimuth[indices])
    return dict(z_m=profile["z_m"].copy(), intensity_mw_cm2=profile["intensity_mw_cm2"].copy(),
                jones_xy=profile["jones_xy"].copy(), carrier_phase_rad=profile["carrier_phase_rad"].copy(),
                azimuth_deg=np.rad2deg(azimuth), ellipticity_deg=np.rad2deg(angles[:, 1]),
                converged=raw["propagation_convergence"]["passed"],
                status=raw["propagation_convergence"]["status"])


def profile_figure(profile):
    import matplotlib.pyplot as plt
    from ..plot_style import PALETTE

    fig, (intensity_ax, angle_ax) = plt.subplots(2, 1, sharex=True, figsize=(8, 5.4))
    z = profile["z_m"]*1000
    intensity_ax.plot(z, profile["intensity_mw_cm2"], color=PALETTE["cyan"], label="Pump intensity")
    intensity_ax.set_ylabel("Intensity [mW/cm^2]")
    angle_ax.plot(z, profile["azimuth_deg"], color=PALETTE["cyan"], label="Azimuth")
    angle_ax.plot(z, profile["ellipticity_deg"], color=PALETTE["rose"], label="Ellipticity")
    angle_ax.set(xlabel="Position in cell [mm]", ylabel="Polarization angle [deg]")
    angle_ax.legend(fontsize=8)
    if not profile["converged"]:
        intensity_ax.set_title("Diagnostic profile: propagation not converged")
    fig.tight_layout()
    return fig
