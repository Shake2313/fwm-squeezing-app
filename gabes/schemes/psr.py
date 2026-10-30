"""Polarization self-rotation with an explicit, conditional noise completion."""
from copy import deepcopy
from functools import lru_cache

import numpy as np

from ..psr import compute_psr_response
from ..psr_noise import psr_gaussian_readout
from .base import ParamSpec, Scheme


PHASE_VIEW = "Quadrature phase scan"
DETUNING_VIEW = "Optical detuning"
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
    caption = ("A linearly polarized carrier drives atomic self-rotation. "
               "The orthogonal vacuum polarization follows an effective shear-and-loss model.")
    supports_headless_observables = True
    cache_version = "psr-direct-selected-point-v2"

    def param_schema(self):
        return [
            ParamSpec("isotope", "Atom / Isotope", "Atomic", "Rb87",
                      choices=("Rb87", "Rb85"),
                      atom_menu=(("Rb", "87", "Rb87", {}), ("Rb", "85", "Rb85", {}))),
            ParamSpec("line", "Transition line", "Atomic", "D1", choices=("D1", "D2")),
            ParamSpec("intensity_mw_cm2", "Carrier intensity", "Light", .1, 0., 10., .01,
                      unit="mW/cm²", help="Uniform local intensity; weak optical excitation is assumed."),
            ParamSpec("temp_c", "Temperature", "Cell", 70., 20., 150., 1., unit="°C"),
            ParamSpec("cell_mm", "Cell length", "Cell", 75., 0., 200., 1., unit="mm"),
            ParamSpec("relaxation_hz", "Ground-state reset rate", "Cell", 1000., 10., 100000., 10.,
                      unit="s⁻¹", help="Isotropic transit replacement; this rate has no 2π conversion."),
            ParamSpec("view", "Scan axis", "Readout", PHASE_VIEW,
                      choices=(PHASE_VIEW, DETUNING_VIEW), control="required_segmented", recompute=False),
            ParamSpec("detuning_mhz", "Optical detuning", "Readout", -900.,
                      -SCAN_LIMIT_MHZ, SCAN_LIMIT_MHZ, 25., unit="MHz", recompute=False,
                      help="Relative to the upper ground-F optical line-strength centroid. "
                           "The selected atomic response is solved directly and cached; fixed during the phase scan."),
            ParamSpec("after_cell_eff_pct", "After-cell transmission", "Detection", 95.,
                      0., 100., 1., unit="%", recompute=False),
            ParamSpec("qe_pct", "Detector efficiency", "Detection", 95.,
                      0., 100., 1., unit="%", recompute=False),
            ParamSpec("excess_noise_per_od", "Assumed atomic excess noise", "Noise model", 0.,
                      0., 20., .1, recompute=False,
                      help="Dimensionless isotropic diffusion per absorption OD. Zero omits atomic "
                           "excess noise; it is not a prediction that the atoms are noiseless."),
            ParamSpec("doppler", "Doppler averaging", "Atomic", True,
                      control="checkbox", advanced=True),
            ParamSpec("resolution", "Resolution", "Numerics", "Balanced",
                      choices=("Quick", "Balanced", "Fine"), advanced=True,
                      help="41 / 81 / 161 optical points and velocity quadrature order 4 / 6 / 10. "
                           "Compare resolutions before interpreting narrow features."),
            ParamSpec("finite_difference_deg", "Ellipticity derivative step", "Numerics", .05,
                      .01, .5, .01, unit="deg", advanced=True,
                      help="Symmetric ±h and ±h/2 responses check the local linearization."),
        ]

    def info(self):
        return """**Conditional PSR polarization squeezing**

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

    def compute(self, params):
        p = dict(self.defaults(), **params)
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
        readout = dict(selected_point=point, phase_scan=phase, spectrum=spectrum,
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
