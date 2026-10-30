"""One Budker-based multihyperfine engine for every magneto-optics view."""
import math

import numpy as np

from .. import magneto_multilevel, species
from ..magneto_noise import magnetic_sensitivity, shot_noise_reference
from ..plot_style import apply_gabes_plot_style
from ..report import derived_table
from .base import ParamSpec, Preset, Scheme


FIELD_SCAN = "Magnetic field"
OPTICAL_SCAN = "Optical detuning"
TRANSMISSION = "Transmission"
ROTATION = "NMOR rotation"
PARAFFIN = "Paraffin coated cell"
BUFFER = "Buffer gas cell"
LINEARITY_STEP_UT = 1e-6
RESPONSE_DENSITY_M3 = 1e16
RESPONSE_LENGTH_M = .1


def _finite_amplitude_linearity(backend, fields, bias, axis):
    """Qualify the displayed response over the existing +/-1 pT interval."""
    center = int(np.argmin(np.abs(fields-bias)))
    left = int(np.argmin(np.abs(fields-(bias-LINEARITY_STEP_UT))))
    right = int(np.argmin(np.abs(fields-(bias+LINEARITY_STEP_UT))))
    h_t = LINEARITY_STEP_UT*1e-6
    linearity = {}
    for quantity, derivative in (("rotation", "rotation_slope_rad_per_t"),
                                  ("transmission", "transmission_slope_per_t")):
        values = np.asarray(backend["rotation_rad" if quantity == "rotation" else "transmission"], float)
        prediction = np.asarray(backend[derivative], float)[center]*h_t
        error = np.maximum(np.abs(values[right]-values[center]-prediction),
                           np.abs(values[center]-values[left]-prediction))/np.maximum(np.abs(prediction), 1e-30)
        passed = np.isfinite(error) & (error <= .05) & np.asarray(backend["slope_valid"], bool)[center]
        if axis == FIELD_SCAN:
            all_passed, all_error = np.zeros(fields.size, bool), np.full(fields.size, np.nan)
            all_passed[center], all_error[center] = passed[0], error[0]
            passed, error = all_passed, all_error
        linearity[quantity] = dict(valid=passed, relative_error=error)
    return center, linearity


def _display_response(raw, params):
    """Apply density and length after the density-independent atomic solve."""
    mode = params["density_mode"]
    if mode not in ("manual", "vapor"):
        raise ValueError("invalid density mode")
    length = float(params["cell_mm"])*1e-3
    density = (float(params["density_cm3"])*1e6 if mode == "manual" else
               species.number_density(species.ISOTOPES[raw["backend"]["isotope"]],
                                      raw["conditions"]["temp_c"]+273.15))
    if not np.isfinite(length) or length < 0 or not np.isfinite(density) or density < 0:
        raise ValueError("density and cell length must be finite and nonnegative")
    atomic = magneto_multilevel.rescale_readout(raw["backend"], density_m3=density,
                                              cell_length_m=length)
    conditions = dict(raw["conditions"], cell_length_m=length,
                      density_m3=density if mode == "manual" else None)
    _, linearity = _finite_amplitude_linearity(
        atomic, raw["physical_b_ut"], float(params["b_offset_ut"]), raw["scan_axis"])
    return atomic, conditions, linearity


def _limit_snr(budget, passed, errors):
    """Finite-amplitude validation qualifies SNR only, never the analytic ASD."""
    budget = dict(budget)
    passed = np.asarray(passed, bool)
    snr = np.asarray(budget["signal_to_noise_ratio"], float)
    available = np.asarray(budget["valid"], bool)
    in_range = budget["signal_rms_t"] <= LINEARITY_STEP_UT*1e-6/np.sqrt(2)
    valid = np.asarray(budget["snr_valid"], bool) & passed & in_range
    status = np.asarray(budget.get("snr_status", budget["status"]), object)
    status[available & ~passed] = "finite-amplitude linearity not validated at this sample"
    if not in_range:
        status[available] = "signal peak exceeds tested local field interval"
    budget.update(
        signal_to_noise_ratio=np.where(valid, snr, np.nan), snr_valid=valid,
        snr_status=status, signal_peak_limit_t=LINEARITY_STEP_UT*1e-6,
        signal_linearity_valid=passed, signal_linearity_relative_error=errors,
        signal_linearity_status=("signal peak exceeds tested local field interval" if not in_range else
            "finite-amplitude check at declared bias; see per-sample status"),
    )
    return _json(budget)


def _json(value):
    if isinstance(value, np.ndarray):
        return _json(value.tolist())
    if isinstance(value, np.generic):
        return _json(value.item())
    if isinstance(value, dict):
        return {str(k): _json(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json(v) for v in value]
    if isinstance(value, float) and not np.isfinite(value):
        return None
    return value


def _transmission_sensitivity(T, slope, *, frequency_hz, intensity_mw_cm2,
                              beam_diameter_mm, efficiency, bandwidth_hz,
                              signal_rms_t, slope_valid):
    """Fixed noiseless incident reference; PSN of normalized transmission."""
    # Obtain only the common photon budget/availability; the transmission
    # responsivity and its units are kept separate from polarimeter rotation.
    photons = magnetic_sensitivity(
        T, np.ones_like(T), frequency_hz=frequency_hz,
        intensity_mw_cm2=intensity_mw_cm2, beam_diameter_mm=beam_diameter_mm,
        efficiency=efficiency, bandwidth_hz=bandwidth_hz, signal_rms_t=0)
    flux = np.asarray(photons["detected_photon_flux_per_s"], float)
    valid = np.asarray(photons["rotation_valid"], bool) & slope_valid & np.isfinite(slope) & (slope != 0)
    signal_asd = np.full(T.shape, np.nan)
    np.divide(np.sqrt(2)*T, np.sqrt(flux), out=signal_asd, where=flux > 0)
    field_asd = np.full(T.shape, np.nan)
    np.divide(signal_asd, np.abs(slope), out=field_asd, where=valid)
    valid &= np.isfinite(field_asd) & (field_asd > 0)
    field_rms = field_asd*np.sqrt(bandwidth_hz)
    snr = np.full(T.shape, np.nan)
    np.divide(signal_rms_t, field_rms, out=snr, where=valid)
    status = np.asarray(photons["rotation_status"], object)
    status[(flux > 0) & ~slope_valid] = "magnetic response not validated"
    status[(flux > 0) & slope_valid & (slope == 0)] = "zero magnetic slope; unbounded field noise"
    return _json(dict(
        source="coherent-light photon shot noise; normalized transmission",
        formula="S_T^1/2 = sqrt(2 T / incident photon flux); S_B^1/2 = S_T^1/2 / abs(dT/dB)",
        asd_convention="one-sided amplitude spectral density",
        signal_convention="RMS magnetic field; local quasistatic response",
        signal_linearity_status="signal amplitude linearity not assessed",
        transmission=T, transmission_slope_per_t=slope,
        transmission_asd_per_sqrt_hz=signal_asd,
        field_asd_t_per_sqrt_hz=np.where(valid, field_asd, np.nan),
        field_rms_t=np.where(valid, field_rms, np.nan),
        signal_to_noise_ratio=np.where(valid, snr, np.nan),
        valid=valid, status=status, snr_valid=valid,
        bandwidth_hz=bandwidth_hz, signal_rms_t=signal_rms_t,
        incident_power_w=photons["incident_power_w"],
        detected_photon_flux_per_s=flux, efficiency=efficiency,
        beam_diameter_mm=beam_diameter_mm,
    ))


class UnifiedMagnetoScheme(Scheme):
    """B and optical axes select a scan, never a different atomic model."""
    name = "magneto"
    title = "Magneto-optics"
    cluster = "C - Magneto-optics"
    caption = ("Rb D1/D2 Hanle transmission and nonlinear magneto-optical rotation, "
               "using a shared multihyperfine Budker-based model for magnetic-field "
               "and optical-detuning scans.")
    cache_version = "budker-unified-v3-optical-operators"
    defaults_version = "budker-unified-v1"
    supports_headless_observables = True

    def __init__(self, mode=None):
        if mode not in (None, "hanle", "eia", "nmor"):
            raise ValueError("unknown magneto-optics alias")
        self.mode = mode
        self.name = mode or "magneto"

    def param_schema(self):
        field = {"scan_axis": FIELD_SCAN}
        optical = {"scan_axis": OPTICAL_SCAN}
        paraffin = {"cell_type": PARAFFIN}
        signal = TRANSMISSION if self.mode in ("hanle", "eia") else ROTATION
        return [
            ParamSpec("scan_axis", "Scan axis", "Readout", FIELD_SCAN,
                      choices=(FIELD_SCAN, OPTICAL_SCAN), control="segmented",
                      help="Both axes use the same multihyperfine atomic model and conditions."),
            ParamSpec("signal_type", "Signal", "Readout", signal,
                      choices=(TRANSMISSION, ROTATION), control="segmented", recompute=False),
            ParamSpec("isotope", "Isotope", "Atomic", "Rb85", choices=("Rb85", "Rb87")),
            ParamSpec("line", "Optical line", "Atomic", "D1", choices=("D1", "D2")),
            ParamSpec("ground_center", "Ground-F frequency reference", "Atomic", "upper",
                      choices=("upper", "lower"),
                      help="Detuning zero is this ground component's line-strength centroid. All hyperfine manifolds are included."),
            ParamSpec("cell_type", "Cell type", "Cell & beams", PARAFFIN,
                      choices=(PARAFFIN, BUFFER), control="segmented"),
            ParamSpec("intensity_mw_cm2", "Beam intensity", "Fields", 1., 0.01, 20., .01, "mW/cm²"),
            ParamSpec("qwp_deg", "Probe polarization (QWP angle)", "Fields", 45. if self.mode == "eia" else 0.,
                      0., 45., .5, "deg", endpoints=("linear", "circular")),
            ParamSpec("laser_detuning_mhz", "Laser detuning", "Detunings", 0., -6000., 6000., 5., "MHz",
                      visible_if=field, help="From the selected ground-F component's line-strength centroid."),
            ParamSpec("detuning_max_mhz", "Optical detuning ±", "Detunings", 1200., 1., 6000., 10., "MHz", visible_if=optical),
            ParamSpec("detuning_points", "Optical scan samples", "Numerics", 49, 3, 201, 2, visible_if=optical),
            ParamSpec("b_max_ut", "B scan ±", "Fields", .001, .00001, 500., .00001, "µT", visible_if=field, format="%.5f"),
            ParamSpec("scan_points", "B scan samples", "Numerics", 81, 11, 401, 2, visible_if=field),
            ParamSpec("b_offset_ut", "Physical longitudinal bias", "Fields", 0., -50., 50., .00001, "µT", format="%.5f",
                      help="B_physical = B_command + bias; optical scans use this physical field."),
            ParamSpec("residual_transverse_b_ut", "Residual transverse B", "Fields", .0001 if self.mode == "eia" else 0.,
                      0., 5., .00001, "µT", advanced=True, format="%.5f"),
            ParamSpec("transverse_field_angle_deg", "Transverse B angle", "Fields", 0., 0., 180., 1., "deg", advanced=True),
            ParamSpec("temp_c", "Temperature", "Cell & beams", 20., 0., 120., 1., "°C"),
            ParamSpec("cell_mm", "Cell length", "Cell & beams", 100., .5, 200., .5, "mm", recompute=False),
            ParamSpec("density_mode", "Vapor density", "Cell & beams", "manual", choices=("manual", "vapor"), recompute=False),
            ParamSpec("density_cm3", "Manual isotope density", "Cell & beams", 4.5e9, 0., 1e12, 1e8, "cm⁻³",
                      visible_if={"density_mode": "manual"}, recompute=False,
                      help="Selected isotope density; abundance is not applied again."),
            ParamSpec("gamma_out_hz", "Light-region exit rate", "Atomic", 5e4, 0., 2e5, 100., "s⁻¹", advanced=True, visible_if=paraffin),
            ParamSpec("gamma_in_hz", "Dark-region return rate", "Atomic", 20., 0., 1000., 1., "s⁻¹", advanced=True, visible_if=paraffin),
            ParamSpec("wall_relaxation_hz", "Wall population-reset rate", "Atomic", 3., .01, 1000., .01, "s⁻¹", advanced=True, visible_if=paraffin),
            ParamSpec("ground_relaxation_hz", "Ground population-reset rate", "Atomic", 1., .01, 1e6, 1., "s⁻¹", advanced=True,
                      help="Independent isotropic reset rate in the illuminated region; sets effective buffer/transit relaxation when dark exchange is absent."),
            ParamSpec("pure_dephasing_hz", "Ground pure dephasing", "Atomic", 0., 0., 1e6, 1., "s⁻¹", advanced=True),
            ParamSpec("velocity_mixing", "Dark-region velocity mixing", "Atomic", "complete", choices=("complete", "none"),
                      advanced=True, visible_if=paraffin),
            ParamSpec("doppler", "Doppler average", "Numerics", "on", choices=("on", "off"), advanced=True),
            ParamSpec("velocity_points", "Velocity quadrature order", "Numerics", 6, 2, 20, 1, advanced=True,
                      help="Order per resonance-centered interval; velocity convergence is assessed separately."),
            ParamSpec("show_shot_noise", "Show shot noise", "Detection", True, control="checkbox", recompute=False),
            ParamSpec("beam_diameter_mm", "Effective beam diameter", "Detection", 3., .1, 20., .1, "mm", recompute=False, advanced=True),
            ParamSpec("detection_efficiency", "Detection efficiency", "Detection", .8, .01, 1., .01, recompute=False, advanced=True),
            ParamSpec("noise_bandwidth_hz", "Noise bandwidth (ENBW)", "Detection", 1., .1, 1e6, .1, "Hz", recompute=False, advanced=True),
            ParamSpec("magnetic_signal_ft", "Small magnetic signal (RMS)", "Detection", 1., 0., 1e6, 1., "fT", recompute=False, advanced=True,
                      help="Quasistatic RMS field signal. SNR is shown only within the tested ±1 pT interval when its linearity check passes; atomic bandwidth is not validated."),
        ]

    def presets(self):
        # Named starting conditions live in the sidebar Defaults menu. They
        # never select another atomic engine or populate the top-bar selector.
        regimes = dict(
            scan_axis=FIELD_SCAN, isotope="Rb87", line="D1", ground_center="upper",
            signal_type=TRANSMISSION, cell_type=PARAFFIN, temp_c=25., cell_mm=10.,
            density_mode="vapor", beam_diameter_mm=2., laser_detuning_mhz=-407.24,
            b_offset_ut=0., qwp_deg=0., scan_points=121,
            ground_relaxation_hz=1., pure_dephasing_hz=0.,
            gamma_out_hz=8e4, gamma_in_hz=1e3, wall_relaxation_hz=1/.003,
            velocity_mixing="complete", doppler="on", transverse_field_angle_deg=0.,
        )
        # These are new-model starting rates in s^-1. The historical model's
        # kHz pure-dephasing parameters are not identical to isotropic reset.
        help_text = ("Starting conditions in the shared multihyperfine model; "
                     "the resulting feature sign is calculated from the selected conditions.")
        return [
            Preset("Hanle", dict(regimes, cell_type=BUFFER, intensity_mw_cm2=.8, b_max_ut=5.,
                   ground_relaxation_hz=2e4, pure_dephasing_hz=2e3,
                   residual_transverse_b_ut=0.), help=help_text),
            Preset("NMOR", dict(regimes, signal_type=ROTATION, intensity_mw_cm2=.7, b_max_ut=.2,
                   gamma_out_hz=5e4, wall_relaxation_hz=1/.005,
                   residual_transverse_b_ut=.02), help=help_text),
        ]

    def recommended_defaults(self, params):
        return None

    def info(self):
        return ("**Budker-based multihyperfine magneto-optics.** Magnetic-field and optical-detuning "
                "scans solve the same ground-F/Zeeman effective model with coherent light shifts, open "
                "hyperfine pumping and coated-cell velocity redistribution. Ground-F coherences are "
                "secularly omitted; optical excitation is adiabatically eliminated.\n\n"
                "- Budker et al., *PRA* **62**, 043403 (2000), https://doi.org/10.1103/PhysRevA.62.043403.\n"
                "- Budker et al., *RMP* **74**, 1153 (2002), https://arxiv.org/abs/physics/0203077.\n"
                "- Malakyan et al., *PRA* **69**, 013817 (2004), https://arxiv.org/abs/physics/0309061.\n"
                "- Reiter and Sørensen, *PRA* **85**, 032111 (2012), https://arxiv.org/abs/1112.2806.\n\n"
                "Field sensitivity uses an analytic derivative of the steady-state equations. Photon "
                "shot noise, finite-amplitude linearity and experimental total noise remain distinct.")

    def compute(self, params):
        p = dict(self.defaults(), **params)
        axis = p["scan_axis"]
        if axis not in (FIELD_SCAN, OPTICAL_SCAN):
            raise ValueError("scan_axis must be Magnetic field or Optical detuning")
        if p["cell_type"] not in (PARAFFIN, BUFFER):
            raise ValueError("unsupported magneto-optics cell type")
        if p["ground_center"] not in ("upper", "lower") or p["density_mode"] not in ("manual", "vapor"):
            raise ValueError("invalid ground reference or density mode")
        if p["doppler"] not in ("on", "off"):
            raise ValueError("Doppler averaging must be on or off")
        bias = float(p["b_offset_ut"])
        if axis == FIELD_SCAN:
            count, extent = int(p["scan_points"]), float(p["b_max_ut"])
            if count != p["scan_points"] or count < 3 or not np.isfinite(extent) or extent <= 0:
                raise ValueError("positive B scan extent and at least three integer samples required")
            command = np.unique(np.r_[np.linspace(-extent, extent, count), 0.,
                                      [-LINEARITY_STEP_UT, LINEARITY_STEP_UT],
                                      -bias if abs(bias) <= extent else 0.])
            fields = command+bias
            delta = np.array([float(p["laser_detuning_mhz"])])
        else:
            count, extent = int(p["detuning_points"]), float(p["detuning_max_mhz"])
            if count != p["detuning_points"] or count < 3 or not np.isfinite(extent) or extent <= 0:
                raise ValueError("positive optical extent and at least three integer samples required")
            delta = np.unique(np.r_[np.linspace(-extent, extent, count), 0.])
            fields, command = bias+np.array([-LINEARITY_STEP_UT, 0., LINEARITY_STEP_UT]), np.array([0.])
        phi = math.radians(float(p["transverse_field_angle_deg"]))
        transverse = float(p["residual_transverse_b_ut"])
        ground_f = (3 if p["isotope"] == "Rb85" else 2)-(p["ground_center"] == "lower")
        conditions = dict(
            isotope=p["isotope"], line=p["line"], ground_f=int(ground_f),
            intensity_mw_cm2=float(p["intensity_mw_cm2"]), temp_c=float(p["temp_c"]),
            cell_length_m=RESPONSE_LENGTH_M,
            gamma_out_hz=float(p["gamma_out_hz"]) if p["cell_type"] == PARAFFIN else 0.,
            gamma_in_hz=float(p["gamma_in_hz"]) if p["cell_type"] == PARAFFIN else 0.,
            ground_relaxation_hz=float(p["ground_relaxation_hz"]),
            wall_relaxation_hz=float(p["wall_relaxation_hz"]),
            pure_dephasing_hz=float(p["pure_dephasing_hz"]),
            velocity_mixing=p["velocity_mixing"], velocity_points=p["velocity_points"],
            doppler=p["doppler"] == "on", qwp_deg=float(p["qwp_deg"]),
            density_m3=RESPONSE_DENSITY_M3,
            b_x_ut=transverse*math.cos(phi), b_y_ut=transverse*math.sin(phi),
            magnetic_derivative=True,
        )
        backend = dict(magneto_multilevel.multilevel_response(delta, b_ut=fields, **conditions))
        backend.setdefault("cell_length_m", RESPONSE_LENGTH_M)
        expected = (fields.size, delta.size)
        for name in ("transmission", "rotation_rad", "rotation_slope_rad_per_t", "transmission_slope_per_t", "slope_valid"):
            if np.asarray(backend[name]).shape != expected:
                raise ValueError(f"multihyperfine backend {name} must have shape {expected}")
        center, linearity = _finite_amplitude_linearity(backend, fields, bias, axis)
        return dict(valid=True, model="Budker multihyperfine", scan_axis=axis,
                    command_b_ut=command, physical_b_ut=fields, detunings_mhz=delta,
                    conditions=conditions, density_mode="manual", cell_type=p["cell_type"],
                    response_reference=dict(density_m3=RESPONSE_DENSITY_M3,
                                            cell_length_m=RESPONSE_LENGTH_M),
                    operating_field_index=center, linearity=linearity, backend=backend)

    def observables(self, raw, params, include_figures=True):
        p = dict(self.defaults(), **params)
        atomic, conditions, linearity = _display_response(raw, p)
        optical = raw["scan_axis"] == OPTICAL_SCAN
        select = (raw["operating_field_index"], slice(None)) if optical else (slice(None), 0)
        x = np.asarray(raw["detunings_mhz"] if optical else raw["command_b_ut"], float)
        T = np.asarray(atomic["transmission"], float)[select]
        theta = np.asarray(atomic["rotation_rad"], float)[select]
        slope_theta = np.asarray(atomic["rotation_slope_rad_per_t"], float)[select]
        slope_t = np.asarray(atomic["transmission_slope_per_t"], float)[select]
        valid = np.asarray(atomic["slope_valid"], bool)[select]
        if atomic["diagnostics"].get("solver_status") == "steady-state validation failed":
            valid = np.zeros_like(valid)
        linear = atomic.get("rotation_status", "").startswith("linear-input")
        rotation_view = p["signal_type"] == ROTATION
        settings = dict(frequency_hz=atomic["frequency_hz"],
            intensity_mw_cm2=raw["conditions"]["intensity_mw_cm2"],
            beam_diameter_mm=float(p["beam_diameter_mm"]), efficiency=float(p["detection_efficiency"]),
            bandwidth_hz=float(p["noise_bandwidth_hz"]), signal_rms_t=float(p["magnetic_signal_ft"])*1e-15)
        rotation_budget = magnetic_sensitivity(T, slope_theta, **settings, linear_input=linear, slope_valid=valid)
        transmission_budget = _transmission_sensitivity(T, slope_t, **settings, slope_valid=valid)
        rotation_budget = _limit_snr(rotation_budget, linearity["rotation"]["valid"],
                                     linearity["rotation"]["relative_error"])
        transmission_budget = _limit_snr(transmission_budget, linearity["transmission"]["valid"],
                                         linearity["transmission"]["relative_error"])
        sensitivity = rotation_budget if rotation_view else transmission_budget
        asd = np.asarray(sensitivity["field_asd_t_per_sqrt_hz"], float)
        available = np.asarray(sensitivity["valid"], bool)
        op = int(np.flatnonzero(x == 0.)[0])
        best = int(np.argmin(np.where(available, asd, np.inf))) if optical and np.any(available) else None
        chosen = best if optical else op
        edge = best in (0, x.size-1) if best is not None else False
        optimum_status = ("no valid sensitivity samples" if optical and best is None else
            "sampled model minimum at scan edge; expand optical scan" if edge else
            "sampled model minimum; optical-grid convergence not assessed" if optical else
            "local response at declared physical bias")
        improvement = float(asd[op]/asd[best]) if best is not None and available[op] else None
        response_status = atomic["diagnostics"].get("magnetic_response_status", "analytic DC response")
        max_od = float(np.max(atomic["optical_depth"]))
        propagation_status = ("OD > 0.1: homogeneous intensity omits longitudinal depletion" if max_od > .1 else
                              "homogeneous intensity; longitudinal depletion omitted")
        cell_scope = ("coated-cell common dark reservoir" if raw["cell_type"] == PARAFFIN else
                      "effective single-region buffer/transit reset; pressure broadening and diffusion omitted")
        reference = atomic["detuning_reference"]
        noise = shot_noise_reference(T, reference_transmission=float(.5*(T[0]+T[-1])),
            frequency_hz=settings["frequency_hz"], intensity_mw_cm2=settings["intensity_mw_cm2"],
            beam_diameter_mm=settings["beam_diameter_mm"], efficiency=settings["efficiency"],
            bandwidth_hz=settings["bandwidth_hz"], linear_input=linear)
        noise["enabled"] = bool(p["show_shot_noise"])
        readout = _json(dict(
            model=raw["model"], scan_axis=raw["scan_axis"], x=x,
            x_unit="MHz" if optical else "µT", detunings_mhz=raw["detunings_mhz"],
            command_b_ut=raw["command_b_ut"], physical_b_ut=raw["physical_b_ut"],
            physical_bias_ut=float(p["b_offset_ut"]), detuning_reference=reference,
            reference_centroid_mhz=atomic["reference_centroid_mhz"], isotope=atomic["isotope"],
            line=atomic["line"], ground_f=atomic["ground_f"], density_m3=atomic["density_m3"],
            density_mode=p["density_mode"], conditions=conditions,
            transmission=T, rotation_rad=theta, rotation_slope_rad_per_t=slope_theta,
            transmission_slope_per_t=slope_t, slope_valid=valid,
            slope_method="analytic derivative of the coupled steady-state equations",
            slope_relative_residual=np.asarray(atomic.get("slope_relative_residual", np.zeros_like(atomic["slope_valid"])), float)[select],
            sensitivity=sensitivity, rotation_sensitivity=rotation_budget,
            transmission_sensitivity=transmission_budget, shot_noise=noise,
            operating_index=op, optimum_index=best,
            optimum_detuning_mhz=None if best is None else float(x[best]),
            optimum_field_asd_t_per_sqrt_hz=None if best is None else float(asd[best]),
            optimum_snr=None if best is None else sensitivity["signal_to_noise_ratio"][best],
            optimum_at_scan_edge=edge, optimum_status=optimum_status,
            centroid_field_asd_t_per_sqrt_hz=float(asd[op]) if optical and available[op] else None,
            snr_improvement_over_centroid=improvement,
            model_scope=atomic["model_scope"], cell_scope=cell_scope,
            rotation_status=atomic["rotation_status"], diagnostics=atomic["diagnostics"],
            velocity_quadrature=atomic["velocity_quadrature"],
            propagation_status=propagation_status, max_optical_depth=max_od,
        ))
        selected_asd = float(asd[chosen]) if chosen is not None and available[chosen] else None
        if optical:
            heroes = [dict(label="Best photon-limited sensitivity", value="unavailable" if selected_asd is None else f"{selected_asd*1e15:.3g} fT/√Hz", tier="hero"),
                      dict(label="Best optical detuning", value="unavailable" if best is None else f"{x[best]:+.4g} MHz", tier="hero")]
        else:
            signal_value = theta[op]*1e3 if rotation_view else T[op]
            heroes = [dict(label="Rotation at declared bias" if rotation_view else "Transmission at declared bias",
                           value="unavailable" if not np.isfinite(signal_value) else f"{signal_value:.4g}"+(" mrad" if rotation_view else ""), tier="hero"),
                      dict(label="PSN magnetic sensitivity", value="unavailable" if selected_asd is None else f"{selected_asd*1e15:.3g} fT/√Hz", tier="hero")]
        for metric in heroes:
            if metric["value"] == "unavailable":
                metric["kind"] = "status"
        metrics = heroes + [
            dict(label="Model scope", value=atomic["model_scope"], kind="status"),
            dict(label="Atomic solver", value=atomic["diagnostics"]["solver_status"], kind="status"),
            dict(label="Magnetic response", value=response_status, kind="status"),
            dict(label="Optimum status", value=optimum_status, kind="status"),
            dict(label="Optical excitation", value=atomic["diagnostics"]["optical_excitation_status"], kind="status"),
            dict(label="Propagation scope", value=propagation_status, kind="status"),
            dict(label="Velocity convergence", value=atomic["velocity_quadrature"]["convergence_status"], kind="status"),
            dict(label="Cell approximation", value=cell_scope, kind="status"),
            dict(label="Signal linearity", value=sensitivity["signal_linearity_status"], kind="status"),
            dict(label="Rotation readout status", value=atomic["rotation_status"], kind="status"),
        ]
        if chosen is not None and sensitivity["signal_to_noise_ratio"][chosen] is not None:
            metrics.append(dict(label="Small-field SNR", value=f"{sensitivity['signal_to_noise_ratio'][chosen]:.3g} ×",
                                help="Quasistatic estimate for the declared RMS field, shown only when the finite-amplitude check passes at the declared bias; atomic bandwidth is not validated."))
        if improvement is not None:
            metrics.append(dict(label="SNR relative to centroid", value=f"{improvement:.3g} ×"))
        fig = self._figure(raw, p, x, T, theta, slope_theta if rotation_view else slope_t, asd, valid,
                           best, reference, noise) if include_figures else None
        rows = [("Atomic model", "Budker multihyperfine effective-operator model"),
            ("Isotope / line", f"{atomic['isotope']} {atomic['line']}"),
            ("Detuning zero", reference),
            ("Frequency reference from fine centroid", f"{atomic['reference_centroid_mhz']:+.6g} MHz"),
            ("Cell", cell_scope), ("Density", f"{atomic['density_m3']:.6g} m⁻³"),
            ("Physical bias", f"{p['b_offset_ut']:.6g} µT"),
            ("Magnetic slope", "analytic steady-state response; no finite-difference masking"),
            ("Optimum", optimum_status), ("Propagation", propagation_status)]
        notes = ("Both scan axes use the same multihyperfine Hamiltonian, optical pumping and cell dynamics. "
            "The frequency reference selects a ground-component centroid, not a truncated transition. "
            "Sensitivity is photon-noise ASD divided by the local analytic magnetic responsivity. "
            "A true zero responsivity has unbounded field noise and is explicitly unavailable; no interpolation "
            "across a physical singularity is used. The effective model requires weak optical excitation and "
            "omits self-consistent optical propagation, optical Zeeman shifts/linear Faraday background and "
            "technical/electronic/atomic noise. Optical-detuning and velocity convergence are separate checks. "
            "Rotation for elliptical/circular input is unsupported; use the transmission readout.")
        return dict(metrics=metrics, hero_count=2, figure=fig, magneto_readout=readout,
                    tables=[derived_table(rows), {"title": "Notes", "markdown": notes}])

    @staticmethod
    def _figure(raw, p, x, T, theta, slope, asd, valid, best, reference, noise):
        import matplotlib.pyplot as plt
        optical = raw["scan_axis"] == OPTICAL_SCAN
        rotation = p["signal_type"] == ROTATION
        fig, axes = plt.subplots(3 if optical else 1, 1, sharex=True,
                                 figsize=(8.5, 8.) if optical else (8.5, 4.6), squeeze=False)
        axes = axes[:, 0]
        if optical:
            lines = [(np.where(valid, slope*1e-6, np.nan), "Rotation slope" if rotation else "Transmission slope",
                      "rad/uT" if rotation else "1/uT"),
                     (T, "Transmission", "1"), (asd*1e15, "Photon-limited field ASD", "fT/sqrt(Hz)")]
        else:
            lines = [(theta*1e3 if rotation else T, "NMOR rotation" if rotation else "Transmission", "mrad" if rotation else "1")]
        for ax, (values, label, unit) in zip(axes, lines):
            line, = ax.plot(x, values, "-", lw=1.8, label=label)
            line._gabes_export_metadata = dict(quantity=label, y_unit=unit,
                x_unit="MHz" if optical else "uT", detuning_reference=reference,
                model="Budker multihyperfine", slope_method="analytic steady-state derivative")
            ax.set_ylabel(label+f" [{unit}]" if unit != "1" else label)
            ax.axvline(0., color="#64748B", ls=":", lw=1., label="Component centroid" if optical else "Command zero")
            if best is not None:
                ax.axvline(x[best], color="#F97316", ls="--", lw=1., label="Sampled optimum")
            if not np.all(np.isfinite(values)):
                invalid = np.flatnonzero(~np.isfinite(values))
                if invalid.size:
                    ax.plot(x[invalid], np.full(invalid.size, .03), "x", color="#D97706",
                            transform=ax.get_xaxis_transform(), label="Unavailable response")
            ax.legend(loc="best", fontsize=8)
        if optical and np.any(np.isfinite(asd)):
            axes[-1].set_yscale("log")
        if not optical and noise["enabled"]:
            trace = np.asarray(noise["rotation_reference_trace_rad"] if rotation else noise["transmission_reference_trace"], float)
            if np.any(np.isfinite(trace)):
                axes[0].plot(x, trace*(1e3 if rotation else 1.), color="#D97706", alpha=.7, lw=.8, label="Shot noise (simulated)")
                axes[0].legend(loc="best", fontsize=8)
        axes[-1].set_xlabel("Optical detuning from ground-component centroid [MHz]" if optical else
                            "Commanded magnetic field B [uT]")
        axes[0].set_title(f"{raw['backend']['isotope']} {raw['backend']['line']} multihyperfine magneto-optics")
        fig.tight_layout()
        return apply_gabes_plot_style(fig)
