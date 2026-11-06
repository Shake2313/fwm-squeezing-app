"""Polarization self-rotation with an explicit, conditional noise completion."""
from copy import deepcopy
from functools import lru_cache

import numpy as np

from ..psr import compute_psr_response
from ..psr_conventions import jones_from_angles, measurement_conventions
from ..psr_full_obe import solve_psr_pump
from ..psr_noise import psr_gaussian_readout
from .base import ParamSpec, Scheme
from . import psr_rf


PHASE_VIEW = "Quadrature phase scan"
DETUNING_VIEW = "Optical detuning"
CONDITIONAL = "Conditional squeezing"
FULL_OBE = "Rb87 D1 steady state"
SCAN_LIMIT_MHZ = 1500.0
_ATOMIC_INPUTS = (
    "isotope", "line", "intensity_mw_cm2", "temp_c", "cell_length_m",
    "ground_relaxation_hz", "doppler", "velocity_order", "ellipticity_step_deg", "density_m3",
)


def _bracket(axis, value):
    """Locate a selected point in the sampled scan; never extrapolate."""
    if not np.isfinite(value) or value < axis[0] or value > axis[-1]:
        return []
    index = int(np.searchsorted(axis, value))
    if index < len(axis) and axis[index] == value:
        return [index]
    return [index-1, index]


def _reasons(raw, keys, indices):
    return list(dict.fromkeys(reason for key in keys for index in indices
                             for reason in raw[key][index]))


@lru_cache(maxsize=128)
def _cached_selected_response(detuning_mhz, atomic_inputs):
    """Cache one actual atomic solve; callers must not expose mutable entries."""
    return compute_psr_response([detuning_mhz], **dict(atomic_inputs))


def _selected_response(raw, selected, indices):
    if len(indices) == 1:
        return raw, indices[0], "direct atomic solve; existing scan sample reused"
    # The raw scan owns the atomic state. Readout controls, including stale
    # caller copies of atomic parameters, must not change that state.
    inputs = tuple((key, raw["conditions"][key]) for key in _ATOMIC_INPUTS)
    return (deepcopy(_cached_selected_response(selected, inputs)), 0,
            "direct atomic solve at selected detuning; cached independently of detection")


class PSRScheme(Scheme):
    name = "psr"
    cluster = "C — Magneto-optics"
    title = "Polarization squeezing (PSR)"
    caption = ("Explore polarization squeezing, the saturated ⁸⁷Rb D1 atomic state, "
               "and its single-velocity or thermal RF noise spectrum.")
    supports_headless_observables = True
    model_comparison_parameter = "calculation"
    workspace_inputs = ("intensity_mw_cm2", "obe_intensity_mw_cm2")
    cache_version = "psr-longitudinal-langevin-rf-v1"

    def atom_conditions(self, params):
        fixed = ({} if params.get("calculation", CONDITIONAL) == CONDITIONAL else
                 {"Atom": "Rb", "Isotope": "87", "Transition line": "D1"})
        return dict(fields=("isotope", "line", "temp_c", "cell_mm"), line="line", fixed=fixed)

    def param_schema(self):
        weak = {"calculation": CONDITIONAL}
        full = {"calculation": (FULL_OBE, psr_rf.PROPAGATED_RF_MODEL)}
        propagated = {"calculation": psr_rf.PROPAGATED_RF_MODEL}
        single = {"calculation": (FULL_OBE, psr_rf.RF_MODEL)}
        microscopic = {"calculation": (FULL_OBE, *psr_rf.RF_MODELS)}
        squeezing = {"calculation": (CONDITIONAL, *psr_rf.RF_MODELS)}
        rf = {"calculation": psr_rf.RF_MODELS}
        return [
            ParamSpec("calculation", "Calculation", "Model", CONDITIONAL,
                      choices=(CONDITIONAL, FULL_OBE, *psr_rf.RF_MODELS)),
            ParamSpec("isotope", "Atom / Isotope", "Atomic", "Rb87",
                      choices=("Rb87", "Rb85"), visible_if=weak,
                      atom_menu=(("Rb", "87", "Rb87", {}), ("Rb", "85", "Rb85", {}))),
            ParamSpec("line", "Transition line", "Atomic", "D1", choices=("D1", "D2"), visible_if=weak),
            ParamSpec("intensity_mw_cm2", "Carrier intensity", "Light", .1, 0., 10., .01,
                      unit="mW/cm²", visible_if=weak, help="Uniform local intensity; weak optical excitation is assumed."),
            ParamSpec("obe_intensity_mw_cm2", "Pump intensity", "Light", .1, 0., 100000., .001,
                      unit="mW/cm²", visible_if=microscopic, help="Entrance intensity in the propagated model; local intensity in the uniform models."),
            ParamSpec("obe_detuning_mhz", "Pump detuning", "Light", -900., -2000., 2000., .1,
                      unit="MHz", visible_if=microscopic, help="Laser frequency relative to the F = 2 optical centroid."),
            ParamSpec("obe_ellipticity_deg", "Pump ellipticity angle", "Light", 0., -45., 45., .1,
                      unit="°", visible_if=full, help="Zero gives linear polarization. ±45° gives circular polarization."),
            ParamSpec("obe_azimuth_deg", "Pump azimuth", "Light", 0., -90., 90., .1,
                      unit="°", visible_if=full, advanced=True),
            ParamSpec("obe_velocity_m_s", "Axial atomic velocity", "Atomic", 0., -600., 600., .1,
                      unit="m/s", visible_if=single, help="Positive velocity follows the light and lowers the frequency seen by the atom."),
            ParamSpec("prop_velocity_model", "Velocity distribution", "Atomic", "single",
                      choices=("single", "thermal"), choice_labels={"single": "Single velocity", "thermal": "Maxwell thermal"},
                      visible_if=propagated,
                      help="Thermal propagation repeats the converged velocity integral along the cell and can take several minutes."),
            ParamSpec("prop_velocity_m_s", "Axial atomic velocity", "Atomic", 0., -600., 600., .1,
                      unit="m/s", visible_if={**propagated, "prop_velocity_model": "single"}),
            ParamSpec("propagation_accuracy", "Length integration", "Numerics", "Balanced",
                      choices=("Quick", "Balanced", "Fine"), advanced=True, visible_if=propagated,
                      help="Doubles the cell partition and compares the output pump, transfer and noise ellipse."),
            ParamSpec("prop_thermal_accuracy", "Thermal integration", "Numerics", "Balanced",
                      choices=("Quick", "Balanced", "Fine"), advanced=True,
                      visible_if={**propagated, "prop_velocity_model": "thermal"}),
            ParamSpec("temp_c", "Temperature", "Cell", 70., 20., 150., 1., unit="°C"),
            ParamSpec("cell_mm", "Cell length", "Cell", 75., 0., 200., 1., unit="mm", visible_if=squeezing),
            ParamSpec("relaxation_hz", "Ground-state reset rate", "Cell", 1000., 10., 100000., 10.,
                      unit="s⁻¹", help="Isotropic transit replacement; this rate has no 2π conversion."),
            ParamSpec("view", "Scan axis", "Readout", PHASE_VIEW,
                      choices=(PHASE_VIEW, DETUNING_VIEW), control="required_segmented", recompute=False, visible_if=weak),
            ParamSpec("rf_view", "Scan axis", "Readout", psr_rf.RF_VIEW,
                      choices=(psr_rf.RF_VIEW, psr_rf.PHASE_VIEW), control="required_segmented",
                      recompute=False, visible_if=rf),
            ParamSpec("rf_frequency_mhz", "Analysis frequency", "Readout", 1.2, .001, 10., .001,
                      unit="MHz", visible_if=rf,
                      help="The phase scan uses this RF. The spectrum samples one decade on either side."),
            ParamSpec("detuning_mhz", "Optical detuning", "Readout", -900.,
                      -SCAN_LIMIT_MHZ, SCAN_LIMIT_MHZ, 25., unit="MHz", recompute=False, visible_if=weak,
                      help="Relative to the upper ground-F optical line-strength centroid. "
                           "The selected atomic response is solved directly and cached; fixed during the phase scan."),
            ParamSpec("after_cell_eff_pct", "After-cell transmission", "Detection", 95.,
                      0., 100., 1., unit="%", recompute=False, visible_if=squeezing),
            ParamSpec("qe_pct", "Detector efficiency", "Detection", 95.,
                      0., 100., 1., unit="%", recompute=False, visible_if=squeezing),
            ParamSpec("excess_noise_per_od", "Assumed atomic excess noise", "Noise model", 0.,
                      0., 20., .1, recompute=False, visible_if=weak,
                      help="Dimensionless isotropic diffusion per absorption OD. Zero omits atomic "
                           "excess noise; it is not a prediction that the atoms are noiseless."),
            ParamSpec("doppler", "Doppler averaging", "Atomic", True,
                      control="checkbox", advanced=True, visible_if=weak),
            ParamSpec("resolution", "Resolution", "Numerics", "Balanced",
                      choices=("Quick", "Balanced", "Fine"), advanced=True, visible_if=weak,
                      help="41 / 81 / 161 optical points and velocity quadrature order 4 / 6 / 10. "
                           "Compare resolutions before interpreting narrow features."),
            ParamSpec("thermal_accuracy", "Thermal integration", "Numerics", "Balanced",
                      choices=("Quick", "Balanced", "Fine"), advanced=True,
                      visible_if={"calculation": psr_rf.THERMAL_RF_MODEL},
                      help="Refines the velocity integral of both response and atomic noise. "
                           "Convergence and the Maxwell-tail check are reported with the result."),
            ParamSpec("finite_difference_deg", "Ellipticity derivative step", "Numerics", .05,
                      .01, .5, .01, unit="deg", advanced=True, visible_if=weak,
                      help="Symmetric ±h and ±h/2 responses check the local linearization."),
        ]

    def info(self):
        return """**Rb87 D1 propagated RF squeezing**

The pump intensity and Jones polarization evolve along the cell. Each segment
recomputes the saturated atomic response and both ordered Langevin noise kernels,
using a single axial velocity or independent Maxwell classes. The segments use
the shared field propagation and composition engines. Noise is measured in the
polarization orthogonal to the actual output carrier. Length refinement and,
for a thermal medium, local velocity convergence must pass before noise is shown.
The Pump propagation figure shows the intensity and
polarization trajectory. This is a local steady-state medium, without atomic
transport memory along the pump gradient or transverse spatial modes.

**Rb87 D1 RF squeezing**

The saturated atom supplies the frequency-dependent response and Langevin noise
of both polarizations. A bright x-polarized carrier and vacuum in y enter a
uniform thin cell. The analysis-frequency view shows the minimum and maximum
noise and the minimum-noise phase. The phase scan uses the selected RF and a
local oscillator referenced to the output carrier. This mode uses a single
velocity class and narrowband detection. **Rb87 D1 thermal RF squeezing**
averages the same drift and both ordered noise kernels over the one-dimensional
Maxwell velocity distribution before propagation. Thermal integration checks
velocity resolution, the truncated tail and the output noise ellipse.
Both modes use a uniform pump and require
|αL| ≤ 0.3. Optical transmission and detector efficiency act after the cell.

**Rb87 D1 atomic steady state**

At one pump detuning and axial velocity, the full optical Bloch equations
give the populations of the 16 hyperfine and Zeeman sublevels, the optical
coherences, and local absorption. Vary the local pump intensity to explore
optical pumping and saturation. Temperature sets the pure-isotope vapor density.
This zero-field calculation describes the mean atomic state.

**Conditional PSR polarization squeezing**

The existing Rb multihyperfine/Zeeman response supplies the small-ellipticity
self-rotation slope **g = d(dψ/dz)/dε** and absorption **α** at zero magnetic
field. The same dipoles, optical pumping, vapor density and Doppler integration
as the multilevel magneto-optics engine are reused. ε is the ellipticity angle
in radians, not the circular intensity imbalance.

For a bright linear carrier and orthogonal input vacuum, the effective model is
**q → q + κp**, with **κ = gL**. Uniform shear and common vacuum absorption are
integrated together along the cell. External transmission and detector loss
are applied once. The plotted quantity is **Var(Sθ)/|⟨S₁⟩|**,
**Sθ = S₂ cos θ + S₃ sin θ**; shot noise is 0 dB.

**Approximations and range.** Weak optical excitation is eliminated while ground
optical pumping remains; ground hyperfine coherences are secularly omitted.
The pump is spatially uniform, velocity classes are independent, and relaxation
is isotropic transit reset. The noise model assumes quasi-static response, one
effective detected mode, small polarization fluctuations and the same loss for
both polarizations. Pump depletion, finite-RF spectra, wall/buffer-gas dynamics,
multimode propagation and microscopic atomic Langevin noise are omitted.
The absorption of the carrier is only a proxy for dark-polarization loss.

**Atomic excess noise is an assumption.** Its control adds nonnegative isotropic
diffusion per absorption optical depth. Zero is an optimistic omission; a
classical rotation curve does not determine atomic quantum noise. Even inside
the screens, the plotted dB values are conditional estimates, not validated
experimental squeezing predictions. No detected carrier means no Stokes-normalized
noise readout.

The screens are virtual excitation ≤ 0.05, |αL| ≤ 0.3, neglected ellipticity
drift ≤ 0.1, and a 1% (plus 10⁻⁶ m⁻¹ absolute) half-step derivative check.
These are engineering guards, not universal physical boundaries. Outside them,
estimates are retained and marked; failed atomic solves or negative absorption
do not supply a passive-noise curve. Velocity convergence must be checked by
changing Resolution. The default point is an illustrative model setting.
The selected detuning is solved directly. The optical scan joins discrete
samples and can miss narrow structure between them; its grid is not a
convergence certificate for the spectrum.

[Matsko et al., PRA 66, 043815 (2002)](https://arxiv.org/abs/quant-ph/0112072)
develop the PSR squeezing mechanism;
[Lezama et al., PRA 77, 013806 (2008)](https://arxiv.org/abs/0710.0609)
include the atomic noise needed for a microscopic prediction.
"""

    def _active_atomic_params(self, params):
        active = {}
        for spec in self.param_schema():
            visible = all(params.get(key) in (allowed if isinstance(allowed, (tuple, list)) else (allowed,))
                          for key, allowed in (spec.visible_if or {}).items())
            if spec.recompute and visible:
                active[spec.name] = params[spec.name]
        thermal = (params["calculation"] == psr_rf.THERMAL_RF_MODEL or
                   (params["calculation"] == psr_rf.PROPAGATED_RF_MODEL and params["prop_velocity_model"] == "thermal"))
        return dict(active, isotope="Rb87", line="D1", doppler=thermal)

    def model_comparison(self, params, *, resolve_raw, current_raw=None, current_view=None):
        from .polarization_comparison import psr_comparison
        return psr_comparison(self, params, resolve_raw=resolve_raw,
                              current_raw=current_raw, current_view=current_view)

    def compute(self, params):
        p = dict(self.defaults(), **params)
        if p["calculation"] in psr_rf.RF_MODELS:
            result = psr_rf.compute_rf(p)
            result["active_params"] = self._active_atomic_params(p)
            return result
        if p["calculation"] == FULL_OBE:
            state = solve_psr_pump(
                p["obe_detuning_mhz"], intensity_mw_cm2=p["obe_intensity_mw_cm2"],
                velocity_m_s=p["obe_velocity_m_s"], temp_c=p["temp_c"],
                ground_relaxation_hz=p["relaxation_hz"],
                jones_xy=jones_from_angles(np.deg2rad(p["obe_azimuth_deg"]),
                                           np.deg2rad(p["obe_ellipticity_deg"])))
            # Keep the solved density matrix in the result bundle. The large
            # reusable generator and reservoir operators belong to the Python API.
            result = dict(kind="psr_full_obe", active_params=self._active_atomic_params(p), **{key: state[key] for key in (
                "model_id", "conditions", "rho", "ground_populations", "excited_populations",
                "ground_levels", "excited_levels", "polarization_xy_c_m2",
                "excited_fraction", "alpha_per_m",
                "azimuth_rate_per_m", "ellipticity_rate_per_m", "diagnostics")})
            # Streamlit supports complex ndarrays but cannot hash a Python
            # complex scalar in the raw input to its observables cache.
            for key in ("chi_plus", "chi_minus", "chi_probe"):
                value = state[key]
                result[key] = None if value is None else {"real": value.real, "imag": value.imag}
            return result
        if p["calculation"] != CONDITIONAL:
            raise ValueError(f"Unknown PSR calculation {p['calculation']!r}")
        count, order = {"Quick": (41, 4), "Balanced": (81, 6), "Fine": (161, 10)}.get(
            p["resolution"], (41, 4))
        doppler = p["doppler"]
        if isinstance(doppler, str):
            doppler = not doppler.lower().startswith("off")
        return compute_psr_response(
            np.linspace(-SCAN_LIMIT_MHZ, SCAN_LIMIT_MHZ, count),
            isotope=p["isotope"], line=p["line"], intensity_mw_cm2=p["intensity_mw_cm2"],
            temp_c=p["temp_c"], cell_length_m=p["cell_mm"]*1e-3,
            ground_relaxation_hz=p["relaxation_hz"], doppler=bool(doppler),
            velocity_order=order, ellipticity_step_deg=p["finite_difference_deg"])

    def observables(self, raw, params, include_figures=True):
        if raw.get("kind") == "psr_rf":
            return psr_rf.rf_view(raw, dict(self.defaults(), **params), include_figures)
        if raw.get("kind") == "psr_full_obe":
            return self._steady_state_view(raw, include_figures)
        p = dict(self.defaults(), **params)
        axis = np.asarray(raw["detunings_mhz"])
        eta = p["after_cell_eff_pct"]*p["qe_pct"]/10000.
        excess = float(p["excess_noise_per_od"])
        carrier = raw["conditions"]["intensity_mw_cm2"] > 0 and eta > 0
        available = (np.asarray(raw["passive_noise_valid"], dtype=bool)
                     & np.asarray(raw["numerical_valid"], dtype=bool)
                     & (np.asarray(raw["optical_depth"]) >= 0)
                     & (np.asarray(raw["transmission"]) > 0) & carrier)
        spectrum = dict(detunings_mhz=axis.copy(), minimum_db=np.full(axis.shape, np.nan),
                        maximum_db=np.full(axis.shape, np.nan), available=available.copy(),
                        approximation_valid=np.asarray(raw["approximation_valid"]).copy())
        if np.any(available):
            noise = psr_gaussian_readout(
                np.asarray(raw["shear"])[available], np.asarray(raw["optical_depth"])[available],
                excess_noise_per_od=excess, detection_efficiency=eta, phase_degrees=[0.])
            for key in ("minimum_db", "maximum_db"):
                spectrum[key][available] = noise[key]

        selected = float(p["detuning_mhz"])
        indices = _bracket(axis, selected)
        point = dict(detuning_mhz=selected, available=False, approximation_valid=False,
                     minimum_db=None, maximum_db=None, shear=None, optical_depth=None,
                     transmission=None, min_angle_deg=None, approximation_reasons=[],
                     status="unavailable: selected detuning outside computed scan",
                     selection_method="unavailable outside computed scan; no extrapolation",
                     atomic_diagnostics=None,
                     bracket_indices=indices)
        phase = None
        if indices:
            selected_raw, selected_index, method = _selected_response(raw, selected, indices)
            point["selection_method"] = method
            point_indices = [selected_index]
            reasons = _reasons(selected_raw, ("numerical_reasons", "passive_noise_reasons"), point_indices)
            point["approximation_reasons"] = _reasons(selected_raw, ("approximation_reasons",), point_indices)
            point["approximation_valid"] = bool(selected_raw["approximation_valid"][selected_index])
            numerical_valid = bool(selected_raw["numerical_valid"][selected_index])
            point["atomic_diagnostics"] = dict(
                conditions=deepcopy(raw["conditions"]),
                numerical_valid=numerical_valid,
                passive_noise_valid=bool(selected_raw["passive_noise_valid"][selected_index]),
                numerical_reasons=list(selected_raw["numerical_reasons"][selected_index]),
                passive_noise_reasons=list(selected_raw["passive_noise_reasons"][selected_index]),
                health={key: float(values[selected_index]) for key, values in selected_raw["atomic_health"].items()},
                derivative_converged=bool(selected_raw["derivative_converged"][selected_index]),
                derivative_error=float(selected_raw["derivative_error"][selected_index]),
                ellipticity_drift_length=float(selected_raw["ellipticity_drift_length"][selected_index]),
                velocity_class_count=int(selected_raw["velocity_quadrature"]["class_counts"][selected_index]),
            )
            if numerical_valid:
                for key in ("shear", "optical_depth"):
                    value = float(selected_raw[key][selected_index])
                    point[key] = value if np.isfinite(value) else None
                if point["optical_depth"] is not None:
                    with np.errstate(over="ignore"):
                        transmission = float(np.exp(-point["optical_depth"]))
                    point["transmission"] = transmission if np.isfinite(transmission) else None
            if not carrier or point["transmission"] == 0:
                reasons.append("no detected bright carrier for Stokes normalization")
            if point["optical_depth"] is not None and point["optical_depth"] < 0:
                reasons.append("negative absorption unsupported by passive-loss noise model")
            selected_available = (numerical_valid and point["atomic_diagnostics"]["passive_noise_valid"]
                                  and point["shear"] is not None and point["optical_depth"] is not None
                                  and point["optical_depth"] >= 0 and point["transmission"] is not None
                                  and point["transmission"] > 0 and carrier)
            if selected_available:
                phase = psr_gaussian_readout(point["shear"], point["optical_depth"],
                    excess_noise_per_od=excess, detection_efficiency=eta)
                point.update(available=True,
                    **{key: float(phase[key]) for key in ("minimum_db", "maximum_db", "min_angle_deg")})
                point["status"] = ("Conditional estimate · within approximation screens" if point["approximation_valid"]
                                   else "Conditional estimate · outside approximation screens")
            else:
                point["status"] = "unavailable: " + "; ".join(reasons or ["invalid atomic response"])

        noise_assumption = ("Atomic excess noise omitted" if excess == 0 else
                            f"Assumed atomic excess noise: {excess:g} per OD")
        readout = dict(measurement_conventions=measurement_conventions(),
                       selected_point=point, phase_scan=phase, spectrum=spectrum,
                       detection_efficiency=eta, excess_noise_per_od=excess,
                       noise_assumption=noise_assumption, model_scope=raw["model_scope"],
                       approximation_thresholds=raw["approximation_thresholds"],
                       conditions=raw["conditions"], velocity_quadrature=raw["velocity_quadrature"],
                       selection_scope="Selected-point screens use its direct atomic solve. The displayed spectrum joins discrete samples; spectral interpolation accuracy is not certified.")
        metrics = []
        for label, key in (("Minimum polarization noise", "minimum_db"),
                           ("Maximum polarization noise", "maximum_db")):
            metric = dict(label=label, tier="hero", value="unavailable", kind="status",
                          help="Conditional PSR model; 0 dB is shot noise, negative values are squeezed.")
            if point["available"]:
                metric.update(value=f"{point[key]:.2f} dB", kind="numeric", approximate=True,
                              delta="Effective noise model")
            metrics.append(metric)
        for label, value in (
            ("Optical detuning", f"{selected:g} MHz"),
            ("Cell transmission", f"{100*point['transmission']:.1f}%" if point["transmission"] is not None else "unavailable"),
            ("PSR shear κ", f"{point['shear']:.3g}" if point["shear"] is not None else "unavailable"),
            ("Minimum-noise phase", f"{point['min_angle_deg']:.1f}°" if point["available"] else "unavailable")):
            metrics.append(dict(label=label, value=value, tier="secondary"))
        metrics.extend([
            dict(label="Model range", value=point["status"], kind="status",
                 help="; ".join(point["approximation_reasons"]) or readout["selection_scope"]),
            dict(label="Noise assumption", value=noise_assumption, kind="status",
                 help="The classical PSR coefficient does not determine atomic Langevin noise. "
                      "Zero excess diffusion is an optimistic omission."),
        ])
        fig = self._figure(spectrum, phase, point, p["view"]) if include_figures else None
        detail = (
            f"| Quantity | Value |\n|---|---|\n"
            f"| Scope | {point['status']} |\n"
            f"| Approximation warnings | {'; '.join(point['approximation_reasons']) or 'none at selected atomic point'} |\n"
            f"| Noise assumption | {noise_assumption}; not derived from the atom |\n"
            f"| Cell noise model | Uniform shear κ=gL with distributed common attenuation αL |\n"
            f"| External efficiency | {eta:.4g} (after-cell transmission × QE, applied once) |\n"
            f"| Stokes normalization | Var(S₂ cos θ + S₃ sin θ) / abs(mean S₁); detected carrier SQL = 1 |\n"
            f"| Detuning reference | {raw['detuning_reference']} |\n"
            f"| Selected point | {point['selection_method']} |\n"
            f"| Selected-point and scan scope | {readout['selection_scope']} |\n"
            f"| Velocity convergence | {raw['velocity_quadrature']['convergence_status']} |\n"
            f"| Atomic model | {raw['model_scope']} |\n"
        )
        return dict(metrics=metrics, hero_count=2, figure=fig, psr_readout=readout,
                    tables=[dict(title="PSR assumptions and applicability", markdown=detail)])

    @staticmethod
    def _steady_state_view(raw, include_figures):
        conditions = raw["conditions"]
        health = raw["diagnostics"]
        susceptibility = complex(raw["chi_probe"]["real"], raw["chi_probe"]["imag"])
        metrics = [
            dict(label="Excited population", value=f"{100*raw['excited_fraction']:.5g}%",
                 tier="hero", help="Fraction of all atoms in the D1 excited manifold."),
            dict(label="Local absorption α", value=f"{raw['alpha_per_m']:.5g} m⁻¹", tier="hero"),
            dict(label="Atom-frame detuning", value=f"{conditions['effective_detuning_mhz']:.4g} MHz",
                 tier="secondary"),
            dict(label="Steady state", value=health["solver_status"], kind="status"),
        ]
        table = (
            "| Atomic steady state | Value |\n|---|---|\n"
            "| Atom | ⁸⁷Rb D1, 8 ground and 8 excited sublevels |\n"
            "| Observable | Local mean atomic response at one axial velocity and zero magnetic field |\n"
            f"| Complex carrier susceptibility | {susceptibility:.6g} |\n"
            f"| Density matrix trace error | {health['trace_error']:.3g} |\n"
            f"| Hermiticity error | {health['hermiticity_error']:.3g} |\n"
            f"| Minimum density eigenvalue | {health['min_eigenvalue']:.3g} |\n"
            f"| Relative steady-state residual | {health['steady_state_relative_residual']:.3g} |\n"
        )
        for label, key in (("Local rotation dψ/dz", "azimuth_rate_per_m"),
                           ("Local ellipticity change dε/dz", "ellipticity_rate_per_m")):
            rate = raw[key]
            value = f"{rate:.5g} rad/m" if np.isfinite(rate) else "undefined at circular polarization"
            table += f"| {label} | {value} |\n"
        fig = None
        if include_figures:
            import matplotlib.pyplot as plt
            fig, axes = plt.subplots(1, 2, figsize=(8, 4), sharex=True)
            for axis, title, sublevels, populations in (
                (axes[0], "Ground state (5S1/2)", raw["ground_levels"], raw["ground_populations"]),
                (axes[1], "Excited state (5P1/2)", raw["excited_levels"], raw["excited_populations"]),
            ):
                for f in sorted({f for f, m in sublevels}):
                    indices = [i for i, (level_f, _) in enumerate(sublevels) if level_f == f]
                    axis.plot([sublevels[i][1] for i in indices], np.asarray(populations)[indices],
                              marker="o", label=f"F = {f}")
                axis.set(title=title, xlabel="Magnetic sublevel mF", ylabel="Fraction of all atoms",
                         xticks=[-2, -1, 0, 1, 2])
                axis.legend()
            fig.tight_layout()
        return dict(metrics=metrics, hero_count=2, figure=fig, export_params=raw["active_params"],
                    psr_readout=dict(model_id=raw["model_id"], conditions=conditions,
                                     measurement_conventions=measurement_conventions()),
                    tables=[dict(title="Atomic populations and solver checks", markdown=table)])

    @staticmethod
    def _figure(spectrum, phase, point, view):
        import matplotlib.pyplot as plt
        from ..plot_style import PALETTE

        fig, ax = plt.subplots(figsize=(8, 4.3))
        if view == DETUNING_VIEW:
            axis = spectrum["detunings_mhz"]
            plot_axis = axis
            plot_approximation = spectrum["approximation_valid"]
            values = {key: spectrum[key] for key in ("minimum_db", "maximum_db")}
            # Insert the directly solved selected point so curve and hero agree
            # between the discrete scan samples without interpolating the atom.
            if point["bracket_indices"] and not np.any(axis == point["detuning_mhz"]):
                index = np.searchsorted(axis, point["detuning_mhz"])
                plot_axis = np.insert(axis, index, point["detuning_mhz"])
                values = {key: np.insert(value, index, point[key] if point["available"] else np.nan)
                          for key, value in values.items()}
                plot_approximation = np.insert(plot_approximation, index, point["approximation_valid"])
            ax.plot(plot_axis, values["minimum_db"], color=PALETTE["cyan"], label="Minimum noise")
            ax.plot(plot_axis, values["maximum_db"], color=PALETTE["rose"], ls="--", label="Maximum noise")
            outside = ~plot_approximation
            if np.any(outside):
                # Mark intervals adjacent to failed scan samples conservatively.
                # The selected point has its own independently evaluated screen.
                interval_outside = outside[:-1] | outside[1:]
                ax.fill_between(np.repeat(plot_axis, 2)[1:-1], 0., 1.,
                                where=np.repeat(interval_outside, 2), transform=ax.get_xaxis_transform(),
                                color="#94A3B8", alpha=.18, label="Outside approximation screens")
            ax.axvline(point["detuning_mhz"], color="#64748B", lw=.9, ls=":", label="Selected detuning")
            ax.set_xlabel("Optical detuning [MHz]")
            ax.set_xlim(axis[0], axis[-1])
        else:
            if phase is not None:
                ax.plot(phase["phase_degrees"], phase["noise_a_db"],
                        color=PALETTE["cyan"], label="Stokes quadrature A")
                ax.plot(phase["phase_degrees"], phase["noise_b_db"],
                        color=PALETTE["rose"], ls="--", label="Stokes quadrature B (+90 deg)")
            else:
                ax.text(.5, .5, point["status"], ha="center", va="center", wrap=True,
                        transform=ax.transAxes, fontsize=10)
            ax.set_xlabel("Polarization quadrature phase [deg]")
            ax.set_xlim(0., 180.)
            ax.set_xticks(np.arange(0., 181., 30.))
        ax.axhline(0., color="#64748B", ls=":", lw=1., label="Shot noise (0 dB)")
        ax.set_ylabel("Polarization noise / shot noise [dB]")
        ax.legend(loc="best", fontsize=8)
        ax.grid(alpha=.2)
        fig.tight_layout()
        return fig
