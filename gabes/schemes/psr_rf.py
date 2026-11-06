"""Controls and readout for single-velocity and thermal uniform-pump PSR."""
from copy import deepcopy
from functools import lru_cache

import numpy as np

from ..constants import C_LIGHT
from ..psr import UNIFORM_PUMP_OD_LIMIT
from ..psr_conventions import measurement_conventions
from ..psr_full_obe import solve_psr_pump
from ..psr_homodyne import psr_homodyne_readout, psr_output_homodyne_readout
from ..psr_spectrum import psr_local_generator, psr_uniform_transfer
from ..quantum.contracts import GeneratorFrequencyAxis
from ..quantum.traveling import NambuTransfer


RF_MODEL = "Rb87 D1 RF squeezing"
THERMAL_RF_MODEL = "Rb87 D1 thermal RF squeezing"
PROPAGATED_RF_MODEL = "Rb87 D1 propagated RF squeezing"
RF_MODELS = (RF_MODEL, THERMAL_RF_MODEL, PROPAGATED_RF_MODEL)
RF_VIEW = "Analysis frequency"
PHASE_VIEW = "Quadrature phase"
MODEL_ID = "psr_rb87_d1_single_velocity_langevin_rf_v1"
# Two decades around the selected, exactly evaluated analysis frequency.
_RF_RATIOS = np.logspace(-1, 1, 25)


@lru_cache(maxsize=4)
def _pump(detuning_mhz, intensity, velocity, temperature, reset):
    return solve_psr_pump(detuning_mhz, intensity_mw_cm2=intensity, velocity_m_s=velocity,
                          temp_c=temperature, ground_relaxation_hz=reset, jones_xy=(1., 0.))


@lru_cache(maxsize=4)
def _local(pump_inputs, center_hz):
    positive = center_hz * _RF_RATIOS
    return psr_local_generator(_pump(*pump_inputs), np.r_[-positive[::-1], positive])


_THERMAL_SETTINGS = {
    "Quick": dict(velocity_order=2, max_order=6, rtol=.02, max_classes=1800),
    "Balanced": dict(velocity_order=4, max_order=10, rtol=.005, max_classes=3000),
    "Fine": dict(velocity_order=6, max_order=14, rtol=.001, max_classes=6000),
}


@lru_cache(maxsize=2)
def _thermal(pump_inputs, center_hz, accuracy):
    from ..psr_thermal import solve_psr_thermal

    detuning, intensity, temperature, reset = pump_inputs
    positive = center_hz * _RF_RATIOS
    if accuracy not in _THERMAL_SETTINGS:
        raise ValueError("unknown thermal integration accuracy")
    return solve_psr_thermal(
        detuning, np.r_[-positive[::-1], positive], intensity_mw_cm2=intensity,
        temp_c=temperature, ground_relaxation_hz=reset, **_THERMAL_SETTINGS[accuracy])


def compute_rf(params):
    if params.get("calculation") == PROPAGATED_RF_MODEL:
        from .psr_propagated import compute_propagated
        return compute_propagated(params, _RF_RATIOS, _THERMAL_SETTINGS)
    center_hz = float(params["rf_frequency_mhz"]) * 1e6
    length_m = float(params["cell_mm"]) * .001
    if not np.isfinite(center_hz) or center_hz <= 0 or not np.isfinite(length_m) or length_m < 0:
        raise ValueError("positive analysis frequency and nonnegative cell length required")
    inputs = tuple(float(params[k]) for k in ("obe_detuning_mhz", "obe_intensity_mw_cm2",
                                               "obe_velocity_m_s", "temp_c", "relaxation_hz"))
    thermal = params.get("calculation") == THERMAL_RF_MODEL
    model_id = "psr_rb87_d1_maxwell_langevin_rf_v1" if thermal else MODEL_ID
    raw = dict(kind="psr_rf", model_id=model_id, center_hz=center_hz, cell_length_m=length_m,
               uniform_pump_od_limit=UNIFORM_PUMP_OD_LIMIT, unavailable_reason=None,
               spectral_observable="Narrowband PSD at each sampled RF, normalized to carrier shot noise",
               local_spatial_frame="Fixed optical envelope at the evaluation point",
               transfer_spatial_frame="Local output-carrier phase, drift K - i beta diag(signs)")
    if thermal:
        raw["thermal_convergence"] = dict(passed=False, status="thermal calculation incomplete")
    if thermal and inputs[1] == 0:
        raw["unavailable_reason"] = "A nonzero bright carrier is required for Stokes normalization."
        raw["thermal_convergence"] = dict(passed=False, status="not evaluated without a carrier")
        return raw
    try:
        if thermal:
            from ..psr_thermal import thermal_convergence

            result = _thermal((inputs[0], inputs[1], inputs[3], inputs[4]),
                              center_hz, params["thermal_accuracy"])
            state = result["pump_state"]
            quadrature = deepcopy(result["convergence"])
            quadrature["history"] = [{k: entry[k] for k in
                ("interval_order", "cutoff_sigma", "class_count", "weight_sum", "elapsed_seconds",
                 "omitted_maxwell_probability")} for entry in result["history"]]
            raw.update(model_id=result["model_id"],
                       velocity_quadrature=quadrature,
                       thermal_performance=deepcopy(result["performance"]),
                       thermal_convergence=dict(passed=False, status="output convergence not evaluated"))
        else:
            state = _pump(*inputs)
        od = float(state["alpha_per_m"] * length_m)
        beta = np.pi * state["conditions"]["optical_frequency_hz"] / C_LIGHT * state["chi_probe"].real
        raw.update(pump_model_id=state["model_id"], conditions=deepcopy(state["conditions"]),
                   pump_diagnostics=deepcopy(state["diagnostics"]), carrier_phase_rate_rad_per_m=float(beta),
                   excited_fraction=state["excited_fraction"], pump_optical_depth=od,
                   carrier_phase_advance_rad=float(beta * length_m))
        if state["conditions"]["intensity_mw_cm2"] == 0:
            raw["unavailable_reason"] = "A nonzero bright carrier is required for Stokes normalization."
            return raw
        if abs(od) > UNIFORM_PUMP_OD_LIMIT:
            raw["unavailable_reason"] = (f"Uniform-pump range exceeded (|\u03b1L| = {abs(od):.4g}). "
                                         f"Reduce cell length or density to reach |\u03b1L| \u2264 {UNIFORM_PUMP_OD_LIMIT:g}.")
            return raw
        if thermal:
            raw["thermal_convergence"] = thermal_convergence(result, length_m)
            if not raw["thermal_convergence"]["passed"]:
                raw["unavailable_reason"] = ("Thermal velocity integral is not converged. "
                                             + raw["thermal_convergence"]["status"]
                                             + " Increase Thermal integration accuracy.")
                return raw
        local = result["local"] if thermal else _local(inputs, center_hz)
        segment = psr_uniform_transfer(local, state, length_m)
        raw.update(omega_rad_s=segment.frequency_axis.omega_rad_s,
                   frequency_frame=segment.frequency_axis.frame, mode_labels=segment.mode_labels,
                   signs=segment.signs, transfer=segment.transfer,
                   added_greater=segment.noise_greater, added_lesser=segment.noise_lesser,
                   local_drift_per_m=local.drift, reservoir_kinds=local.reservoir_names,
                   local_greater_per_m_by_kind=local.noise_greater_by_reservoir,
                   local_lesser_per_m_by_kind=local.noise_lesser_by_reservoir,
                   local_audit=local.audit(), transfer_audit=segment.audit())
    except (FloatingPointError, np.linalg.LinAlgError, ValueError) as error:
        raw["unavailable_reason"] = f"RF sweep unavailable. Atomic or propagation check failed. {error}"
    return raw


def rf_view(raw, params, include_figures):
    eta = float(params["after_cell_eff_pct"] * params["qe_pct"] / 10000.)
    reason = raw["unavailable_reason"]
    if not 0 <= eta <= 1:
        raise ValueError("external detection efficiency must be in [0,1]")
    if eta == 0 and not reason:
        reason = "A detected carrier is required. External efficiency is zero."
    propagated = raw.get("propagation_convergence") is not None
    noise = None
    if not reason:
        transfer = NambuTransfer(GeneratorFrequencyAxis(raw["omega_rad_s"], raw["frequency_frame"]),
                                 raw["signs"], raw["mode_labels"], raw["transfer"],
                                 raw["added_greater"], raw["added_lesser"])
        try:
            noise = (psr_output_homodyne_readout(transfer, raw["output_jones_xy"], detection_efficiency=eta)
                     if propagated else psr_homodyne_readout(transfer, detection_efficiency=eta))
        except (FloatingPointError, np.linalg.LinAlgError, ValueError) as error:
            reason = f"RF sweep unavailable. Homodyne sideband check failed. {error}"

    point = dict(frequency_hz=raw["center_hz"], available=noise is not None, reason=reason,
                 minimum_db=None, maximum_db=None, min_angle_deg=None)
    phase, spectrum = None, None
    if noise is not None:
        index = len(noise["frequencies_hz"]) // 2
        for key in ("minimum_db", "maximum_db", "min_angle_deg"):
            value = noise[key][index]
            point[key] = float(value) if np.isfinite(value) else None
        phase = dict(phase_degrees=noise["phase_degrees"],
                     **{k: noise[k][index] for k in ("noise_a", "noise_b", "noise_a_db", "noise_b_db")})
        spectrum = {k: noise[k] for k in ("frequencies_hz", "minimum_db", "maximum_db", "min_angle_deg",
                                          "angle_defined", "source_covariance", "covariance",
                                          "source_sideband_covariance", "sideband_covariance", "convention")}
    polarization_readout = None
    if propagated and noise is not None:
        polarization_readout = {key: noise[key] for key in
            ("output_jones_xy", "canonical_output_jones_xy", "orthogonal_jones_xy",
             "output_jones_phase_rad", "lo_offset_rad", "polarization_phase_chart", "output_polarization_analyzer")}
    pump_profile = None
    if propagated:
        from .psr_propagated import profile_readout
        pump_profile = profile_readout(raw)
    readout = dict(model_id=raw["model_id"], selected_point=point, spectrum=spectrum, phase_scan=phase,
                   detection_efficiency=eta, measurement_conventions=measurement_conventions(),
                   thermal_convergence=deepcopy(raw.get("thermal_convergence")),
                   propagation_convergence=deepcopy(raw.get("propagation_convergence")),
                   propagation_velocity_convergence=deepcopy(raw.get("propagation_velocity_convergence")),
                   pump_profile=pump_profile, output_polarization=polarization_readout)
    metrics = [dict(label=label, value=(f"{point[key]:.5g} dB" if point[key] is not None else "unavailable"),
                    tier="hero", help="Orthogonal polarization noise at the selected RF, relative to shot noise.")
               for label, key in (("Minimum RF noise", "minimum_db"), ("Maximum RF noise", "maximum_db"))]
    metrics += [dict(label="Analysis frequency", value=f"{raw['center_hz']/1e6:g} MHz", tier="secondary"),
                dict(label="Minimum-noise phase", tier="secondary",
                     value=f"{point['min_angle_deg']:.2f}°" if point["min_angle_deg"] is not None else "undefined")]
    if propagated and "output_intensity_mw_cm2" in raw:
        metrics.append(dict(label="Output pump intensity", tier="secondary",
                            value=(f"{raw['output_intensity_mw_cm2']:.5g} mW/cm\u00b2"
                                   if raw["propagation_convergence"]["passed"] else "unavailable")))
    if reason:
        metrics.append(dict(label="RF spectrum status", value=reason, kind="status"))
    source = ("Independent Maxwell velocity classes" if raw.get("thermal_convergence") is not None
              else "Single axial velocity")
    polarization = "evolving Jones pump" if propagated else "x-linear pump"
    detail = ("| Quantity | Value |\n|---|---|\n"
              f"| Source | {source}, {polarization}, coherent input with vacuum fluctuations |\n"
              "| Atomic noise | Spontaneous emission and thermal replacement |\n"
              "| Phase reference | Mean output carrier including its dispersive phase |\n"
              "| Detection | Narrowband Stokes PSD, 0 dB at shot noise |\n"
              f"| External transmission × QE | {eta:.5g} |\n")
    if raw.get("thermal_convergence") is not None:
        audit = raw["thermal_convergence"]
        detail += f"| Velocity convergence | {audit['status']} |\n"
        quadrature, performance = raw.get("velocity_quadrature"), raw.get("thermal_performance")
        if quadrature is not None and performance is not None:
            detail += (f"| Final interval order / cutoff | {quadrature['interval_order']} / "
                       f"{quadrature['cutoff_sigma']:g} sigma |\n"
                       f"| Velocity-class evaluations | {performance['class_evaluations']} |\n"
                       f"| Initial integration time | {performance['elapsed_seconds']:.3g} s (reused on cache hits) |\n")
    if propagated:
        detail += f"| Length convergence | {raw['propagation_convergence']['status']} |\n"
        velocity_audit = raw.get("propagation_velocity_convergence")
        if velocity_audit is not None and raw.get("velocity_model") == "thermal":
            detail += f"| Propagated velocity-kernel check | {velocity_audit['passed']} |\n"
        if "output_intensity_mw_cm2" in raw and raw["propagation_convergence"]["passed"]:
            detail += (f"| Input / output intensity | {raw['conditions']['incident_intensity_mw_cm2']:.5g} / "
                       f"{raw['output_intensity_mw_cm2']:.5g} mW/cm\u00b2 |\n"
                       f"| Completed segments | {raw['completed_segments']} |\n"
                       "| Atomic response | Local stationary OBE at each pump midpoint |\n")
    if "pump_optical_depth" in raw:
        detail += (f"| Local pump αL | {raw['pump_optical_depth']:.5g} "
                   f"(uniform-intensity screen ≤ {raw['uniform_pump_od_limit']:g} in magnitude) |\n")
    if reason:
        detail += f"| Status | {reason} |\n"
    figure = _figure(spectrum, phase, point, params["rf_view"]) if include_figures else None
    figure_views = []
    if include_figures and pump_profile is not None:
        from .psr_propagated import profile_figure
        figure_views = [dict(label="RF noise", figure=figure),
                        dict(label="Pump propagation", figure=profile_figure(pump_profile))]
    return dict(metrics=metrics, hero_count=2, psr_readout=readout,
                export_params=dict(raw["active_params"], rf_view=params["rf_view"],
                                   qe_pct=params["qe_pct"], after_cell_eff_pct=params["after_cell_eff_pct"]),
                tables=[dict(title="RF spectrum conditions", markdown=detail)],
                figure=figure, figure_views=figure_views)


def _figure(spectrum, phase, point, view):
    import matplotlib.pyplot as plt
    from ..plot_style import PALETTE

    if spectrum is None:
        fig, ax = plt.subplots(figsize=(8, 3.6))
        ax.set_axis_off()
        ax.text(.5, .5, point["reason"], transform=ax.transAxes, ha="center", va="center", wrap=True)
        return fig
    if view == PHASE_VIEW:
        fig, ax = plt.subplots(figsize=(8, 4.3))
        for key, label, color in (("noise_a_db", "Stokes A", "cyan"), ("noise_b_db", "Stokes B (+90 deg)", "rose")):
            ax.plot(phase["phase_degrees"], phase[key], label=label, color=PALETTE[color])
        ax.set(xlabel="Homodyne phase relative to output carrier [deg]", xlim=(0, 180))
    else:
        fig, (ax, angle_ax) = plt.subplots(2, 1, sharex=True, figsize=(8, 5.4),
                                           gridspec_kw={"height_ratios": [2, 1]})
        frequency = spectrum["frequencies_hz"] / 1e6
        ax.semilogx(frequency, spectrum["minimum_db"], color=PALETTE["cyan"], label="Minimum noise")
        ax.semilogx(frequency, spectrum["maximum_db"], color=PALETTE["rose"], label="Maximum noise")
        angle_ax.semilogx(frequency, spectrum["min_angle_deg"], ".", color=PALETTE["cyan"])
        angle_ax.set(xlabel="Analysis frequency [MHz]", ylabel="Minimum phase [deg]", ylim=(0, 180))
        for panel in (ax, angle_ax):
            panel.axvline(point["frequency_hz"]/1e6, color="#64748B", ls=":", lw=.8)
    ax.axhline(0., color="#64748B", ls="--", lw=1, label="Shot noise")
    ax.set_ylabel("Noise / shot noise [dB]")
    ax.legend(fontsize=8)
    fig.tight_layout()
    return fig
