"""Opt-in multilevel NMOR detuning scan and conditional DC magnetometry."""
import numpy as np

from .. import magneto_multilevel
from ..magneto_noise import magnetic_sensitivity
from ..plot_style import apply_gabes_plot_style
from ..report import derived_table
from .base import ParamSpec


LEGACY = "Magnetic-field scan"
BUDKER = "Budker NMOR detuning scan"
SLOPE_RTOL = 0.05


def param_schema():
    """Physical scan controls; detector controls belong to the host scheme."""
    visible = {"calculation": BUDKER}

    def knob(name, label, default, low=None, high=None, step=None, unit="", **kw):
        return ParamSpec("budker_" + name, label, "Budker NMOR", default,
                         low, high, step, unit, recompute=True,
                         visible_if=visible, **kw)

    return [
        knob("isotope", "Isotope", "Rb85", choices=("Rb85", "Rb87")),
        knob("line", "Optical line", "D1", choices=("D1", "D2")),
        knob("ground_center", "Ground-F reference component", "upper",
             choices=("upper", "lower"),
             help="Detuning zero is the line-strength centroid of this ground-F optical component; all ground and excited hyperfine manifolds remain in the model."),
        knob("density_mode", "Vapor density", "manual", choices=("manual", "vapor"),
             help="Manual density is the selected isotope density. Vapor mode uses the selected isotope's vapor-pressure model."),
        knob("density_cm3", "Manual isotope density", 4.5e9, 0.0, 1.0e12, 1.0e8, "cm⁻³",
             help="Used in manual mode; no isotope abundance is applied again."),
        knob("cell_mm", "NMOR cell length", 100.0, 1.0, 200.0, 1.0, "mm"),
        knob("detuning_max_mhz", "Optical detuning ±", 1200.0, 1.0, 6000.0, 10.0, "MHz"),
        knob("detuning_points", "Optical scan samples", 25, 3, 201, 2,
             help="Zero detuning is included exactly. The optimum is sampled, not a continuous-frequency fit."),
        knob("bias_ut", "Physical longitudinal field", 0.0, -1.0, 1.0, 0.001, "µT"),
        knob("local_step_ut", "Local derivative outer step", 1.0e-5, 1.0e-6, 0.01, 1.0e-6, "µT",
             format="%.6f", advanced=True,
             help="Compare central slopes at ±h and ±h/2; agreement within 5% is required for sensitivity."),
        knob("gamma_out_hz", "Light-region exit rate", 5.0e4, 0.0, 2.0e5, 100.0, "s⁻¹", advanced=True),
        knob("gamma_in_hz", "Dark-region return rate", 20.0, 0.0, 1000.0, 1.0, "s⁻¹", advanced=True),
        knob("wall_relaxation_hz", "Wall population-reset rate", 3.0, 0.01, 1000.0, 0.01, "s⁻¹", advanced=True),
        knob("ground_relaxation_hz", "Light-region population-reset rate", 1.0, 0.01, 1000.0, 0.01, "s⁻¹", advanced=True),
        knob("velocity_mixing", "Dark-region velocity mixing", "complete",
             choices=("complete", "none"), advanced=True,
             help="Complete mixing returns atoms from a common velocity-randomized dark reservoir; none retains velocity classes."),
        knob("velocity_points", "Velocity quadrature order", 6, 2, 20, 1, advanced=True,
             help="Order per resonance-centered interval; change this to assess velocity convergence."),
        knob("doppler", "NMOR Doppler average", "on", choices=("on", "off"), advanced=True),
    ]


def _settings(params):
    return {spec.name: params.get(spec.name, spec.default) for spec in param_schema()}


def _json_value(value):
    if isinstance(value, np.ndarray):
        return _json_value(value.tolist())
    if isinstance(value, np.generic):
        return _json_value(value.item())
    if isinstance(value, dict):
        return {str(key): _json_value(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_value(item) for item in value]
    if isinstance(value, float) and not np.isfinite(value):
        return None
    return value


def compute(params):
    """Cache the five-field atomic stencil; no detection settings enter here."""
    p = _settings(params)
    extent = float(p["budker_detuning_max_mhz"])
    count = int(p["budker_detuning_points"])
    step = float(p["budker_local_step_ut"])
    bias = float(p["budker_bias_ut"])
    if (not np.all(np.isfinite([extent, step, bias])) or extent <= 0 or step <= 0
            or count < 3 or count != p["budker_detuning_points"]):
        raise ValueError("NMOR detuning extent and derivative step must be positive; at least three integer scan samples are required.")
    isotope = str(p["budker_isotope"])
    center = str(p["budker_ground_center"])
    if isotope not in ("Rb85", "Rb87") or center not in ("upper", "lower"):
        raise ValueError("NMOR isotope and ground reference must be Rb85/Rb87 and upper/lower.")
    if p["budker_density_mode"] not in ("manual", "vapor"):
        raise ValueError("NMOR density mode must be manual or vapor.")
    delta = np.unique(np.r_[np.linspace(-extent, extent, count), 0.0])
    fields = bias + step * np.array([-1.0, -0.5, 0.0, 0.5, 1.0])
    if not np.all(np.diff(fields) > 0):
        raise ValueError("Magnetic derivative step is not resolved at the chosen bias.")
    ground_f = (3 if isotope == "Rb85" else 2) - (center == "lower")
    conditions = dict(
        isotope=isotope, line=str(p["budker_line"]), ground_f=int(ground_f),
        intensity_mw_cm2=float(params.get("intensity_mw_cm2", 0.8)),
        temp_c=float(params.get("temp_c", 20.0)),
        cell_length_m=float(p["budker_cell_mm"]) * 1e-3,
        gamma_out_hz=float(p["budker_gamma_out_hz"]),
        gamma_in_hz=float(p["budker_gamma_in_hz"]),
        ground_relaxation_hz=float(p["budker_ground_relaxation_hz"]),
        wall_relaxation_hz=float(p["budker_wall_relaxation_hz"]),
        velocity_mixing=str(p["budker_velocity_mixing"]),
        velocity_points=int(p["budker_velocity_points"]),
        doppler=p["budker_doppler"] == "on", qwp_deg=0.0,
        density_m3=(float(p["budker_density_cm3"]) * 1e6
                    if p["budker_density_mode"] == "manual" else None),
    )
    result = magneto_multilevel.multilevel_response(delta, b_ut=fields, **conditions)
    rotation = np.asarray(result["rotation_rad"], dtype=float)
    transmission = np.asarray(result["transmission"], dtype=float)
    if rotation.shape != (5, delta.size) or transmission.shape != rotation.shape:
        raise ValueError("Multilevel NMOR backend must return (five fields, detunings) arrays.")
    # µT -> T happens once; these slopes therefore have SI units rad/T.
    step_t = step * 1e-6
    coarse = (rotation[4] - rotation[0]) / (2 * step_t)
    fine = (rotation[3] - rotation[1]) / step_t
    scale = np.maximum(np.maximum(np.abs(fine), np.abs(coarse)), np.finfo(float).tiny)
    relative = np.abs(fine - coarse) / scale
    resolved = (np.isfinite(fine) & np.isfinite(coarse) & (relative <= SLOPE_RTOL)
                & (np.signbit(fine) == np.signbit(coarse)))
    if result["diagnostics"].get("solver_status") == "steady-state validation failed":
        resolved[:] = False
    return dict(
        budker_scan=True, valid=True, detunings_mhz=delta, b_ut=fields,
        bias_ut=bias, local_step_ut=step, slope_rad_per_t=fine,
        slope_coarse_rad_per_t=coarse, slope_relative_change=relative,
        slope_valid=resolved, backend=result, conditions=conditions,
        density_mode=p["budker_density_mode"],
    )


def observables(raw, params, include_figures=True):
    """Cheap detector and optimization layer; never invokes an atomic solve."""
    atomic = raw["backend"]
    delta = np.asarray(raw["detunings_mhz"], dtype=float)
    transmission = np.asarray(atomic["transmission"], dtype=float)[2]
    rotation = np.asarray(atomic["rotation_rad"], dtype=float)[2]
    slope = np.asarray(raw["slope_rad_per_t"], dtype=float)
    trusted = np.asarray(raw["slope_valid"], dtype=bool)
    sensitivity = magnetic_sensitivity(
        transmission, np.where(trusted, slope, np.nan),
        frequency_hz=float(atomic["frequency_hz"]),
        intensity_mw_cm2=raw["conditions"]["intensity_mw_cm2"],
        beam_diameter_mm=float(params.get("beam_diameter_mm", 2.0)),
        efficiency=float(params.get("detection_efficiency", 0.8)),
        bandwidth_hz=float(params.get("noise_bandwidth_hz", 1000.0)),
        signal_rms_t=float(params.get("magnetic_signal_ft", 1.0)) * 1e-15,
        signal_peak_limit_t=float(raw["local_step_ut"]) * 0.5e-6,
        linear_input=True, slope_valid=trusted)
    asd = np.asarray(sensitivity["field_asd_t_per_sqrt_hz"], dtype=float)
    available = np.asarray(sensitivity["valid"], dtype=bool)
    zero = int(np.flatnonzero(delta == 0.0)[0])
    best = int(np.argmin(np.where(available, asd, np.inf))) if np.any(available) else None
    edge = best in (0, delta.size - 1) if best is not None else False
    optimum_status = ("no valid sensitivity samples" if best is None else
                      "sampled model minimum at scan edge; expand optical scan" if edge else
                      "sampled model minimum; detuning convergence not assessed")
    improvement = (float(asd[zero] / asd[best])
                   if best is not None and available[zero] else None)
    reference = str(atomic["detuning_reference"])
    scope = ("Conditional weak-optical-excitation, scalar longitudinal, linear-input, "
             "quasi-static NMOR model; photon shot noise only.")
    optical_depth = np.asarray(atomic["optical_depth"], dtype=float)
    max_od = float(np.max(optical_depth))
    propagation_status = (
        "OD > 0.1: homogeneous intensity omits longitudinal depletion; high-OD quantitative use unsupported"
        if max_od > 0.1 else
        "homogeneous intensity; longitudinal depletion omitted")
    readout = _json_value(dict(
        calculation=BUDKER, detunings_mhz=delta, detuning_reference=reference,
        reference_centroid_mhz=atomic["reference_centroid_mhz"],
        isotope=atomic["isotope"], line=atomic["line"], ground_f=atomic["ground_f"],
        density_m3=atomic["density_m3"], density_mode=raw["density_mode"],
        physical_bias_ut=raw["bias_ut"], local_step_ut=raw["local_step_ut"],
        transmission=transmission, rotation_rad=rotation,
        slope_rad_per_t=np.where(trusted, slope, np.nan),
        slope_valid=trusted, slope_relative_change=raw["slope_relative_change"],
        slope_tolerance=SLOPE_RTOL,
        slope_status=["nested stencil consistent" if value else "unresolved local derivative"
                      for value in trusted],
        sensitivity=sensitivity, optimum_index=best,
        optimum_detuning_mhz=None if best is None else float(delta[best]),
        optimum_field_asd_t_per_sqrt_hz=None if best is None else float(asd[best]),
        optimum_snr=None if best is None else sensitivity["signal_to_noise_ratio"][best],
        optimum_at_scan_edge=edge, optimum_status=optimum_status,
        centroid_index=zero,
        centroid_field_asd_t_per_sqrt_hz=float(asd[zero]) if available[zero] else None,
        snr_improvement_over_centroid=improvement,
        numerical_scope="five-field nested derivative check; optical-detuning and velocity convergence must be assessed separately",
        scope=scope, model_scope=atomic["model_scope"],
        max_optical_depth=max_od, propagation_status=propagation_status,
        diagnostics=atomic["diagnostics"], velocity_quadrature=atomic["velocity_quadrature"],
        conditions=raw["conditions"],
    ))
    metrics = [
        dict(label="Best photon-limited sensitivity",
             value="unavailable" if best is None else f"{asd[best] * 1e15:.3g} fT/√Hz",
             tier="hero", **({"kind": "status"} if best is None else {}),
             help="Minimum sampled one-sided field ASD; conditional model, excluding technical and atomic noise."),
        dict(label="Best optical detuning",
             value="unavailable" if best is None else f"{delta[best]:+.4g} MHz",
             tier="hero", **({"kind": "status"} if best is None else {}), help=reference),
        dict(label="Optimum status", value=optimum_status, kind="status"),
        dict(label="Model scope", value=scope, kind="status"),
        dict(label="Velocity convergence", value=atomic["velocity_quadrature"]["convergence_status"], kind="status"),
        dict(label="Optical excitation", value=atomic["diagnostics"]["optical_excitation_status"], kind="status"),
        dict(label="Atomic solver", value=atomic["diagnostics"].get("solver_status", "not assessed"), kind="status"),
        dict(label="Signal linearity", value=sensitivity["signal_linearity_status"], kind="status"),
        dict(label="Propagation scope", value=propagation_status, kind="status"),
    ]
    if improvement is not None:
        metrics.append(dict(label="SNR relative to centroid", value=f"{improvement:.3g} ×",
                            help="Same RMS magnetic signal, incident power, detector and ENBW at the optimum and optical-component centroid."))
    if best is not None and sensitivity["signal_to_noise_ratio"][best] is not None:
        metrics.append(dict(label="Small-field SNR", value=f"{sensitivity['signal_to_noise_ratio'][best]:.3g} ×",
                            help="Declared RMS field divided by photon-limited RMS field noise in the selected ENBW."))
    fig = None
    if include_figures:
        import matplotlib.pyplot as plt

        fig, axes = plt.subplots(3, 1, sharex=True, figsize=(8.5, 8.0))
        plotted = [
            (np.where(trusted, slope * 1e-6, np.nan), "Rotation slope", "rad/uT", "Rotation slope [rad/uT]"),
            (transmission, "Transmission", "1", "Transmission"),
            (asd * 1e15, "Photon-limited field ASD", "fT/sqrt(Hz)", "Field ASD [fT/sqrt(Hz)]"),
        ]
        for ax, (values, label, unit, ylabel) in zip(axes, plotted):
            line, = ax.plot(delta, values, "o-", markersize=3, label=label)
            line._gabes_export_metadata = dict(quantity=label, unit=unit, y_unit=unit,
                x_unit="MHz", x_label="Optical detuning from ground-component centroid",
                detuning_reference=reference, scope=scope,
                physical_bias_ut=float(raw["bias_ut"]), slope_tolerance=SLOPE_RTOL)
            ax.set_ylabel(ylabel)
            ax.axvline(0.0, color="#64748B", linestyle=":", linewidth=1, label="Component centroid")
            if best is not None:
                ax.axvline(delta[best], color="#F97316", linestyle="--", linewidth=1, label="Sampled optimum")
            ax.legend(loc="best", fontsize=8)
        if np.any(available):
            axes[2].set_yscale("log")
        else:
            axes[2].text(0.5, 0.5, "No resolved field sensitivity", ha="center", va="center", transform=axes[2].transAxes)
        axes[-1].set_xlabel("Optical detuning from ground-component centroid [MHz]")
        axes[0].set_title(f"{atomic['isotope']} {atomic['line']} NMOR - conditional DC prediction")
        fig.tight_layout()
        apply_gabes_plot_style(fig)
    rows = [
        ("Isotope / line / reference ground F", f"{atomic['isotope']} {atomic['line']} / F={atomic['ground_f']}"),
        ("Optical detuning zero", reference),
        ("Reference offset from fine-structure centroid", f"{atomic['reference_centroid_mhz']:+.6g} MHz"),
        ("Isotope density", f"{atomic['density_m3']:.6g} m⁻³ ({raw['density_mode']})"),
        ("Physical bias / outer derivative step", f"{raw['bias_ut']:.6g} / {raw['local_step_ut']:.6g} µT"),
        ("Derivative acceptance", "±h and ±h/2 central slopes agree within 5%"),
        ("ENBW / RMS magnetic signal", f"{sensitivity['bandwidth_hz']:.6g} Hz / {sensitivity['signal_rms_t'] * 1e15:.6g} fT"),
        ("Optimum", optimum_status),
        ("Velocity average", atomic["velocity_quadrature"]["convergence_status"]),
        ("Optical elimination", atomic["diagnostics"]["optical_excitation_status"]),
        ("Maximum optical depth", f"{max_od:.6g}"),
        ("Propagation", propagation_status),
    ]
    notes = (scope + " All ground and excited hyperfine manifolds are retained in the "
             "optically eliminated model. The field derivative and transmitted photon flux "
             "refer to the same physical bias and optical detuning. No off-resonance optimum "
             "is imposed. The centroid is a spectral reference, not an individual hyperfine "
             "transition resonance. A sampled minimum needs optical-grid and velocity-order "
             "checks before a quantitative optimum is reported. Linear polarization is fixed "
             "in this mode. Finite-frequency response, technical/electronic noise, spin noise "
             "and longitudinal pump depletion are outside this calculation.\n\n"
             "Reference: [Budker et al., Phys. Rev. A 62, 043403 (2000)]"
             "(https://doi.org/10.1103/PhysRevA.62.043403).")
    return dict(metrics=metrics, hero_count=2, figure=fig, magneto_readout=readout,
                tables=[derived_table(rows), {"title": "Notes", "markdown": notes}])
