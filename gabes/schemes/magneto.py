"""
Unregistered 87Rb D1 single-pair reference and public unified scheme export.

One Hanle-configuration transmission engine now covers EIT-like dips and
EIA/MIA-like peaks.  The zero-field feature sign is a result, not a readout
selector: probe ellipticity, cell relaxation, residual transverse field, and
the addressed Zeeman manifold decide the line shape.
"""
import math
import re
from dataclasses import replace

import numpy as np

from .. import constants, core, kernels, observables, pole_doppler, species, zeeman
from ..magneto_noise import (DEFAULT_BEAM_DIAMETER_MM, DEFAULT_DETECTION_EFFICIENCY,
                             DEFAULT_NOISE_BANDWIDTH_HZ, shot_noise_reference,
                             magnetic_sensitivity)
from ..report import derived_table
from .base import ParamSpec, Preset, Scheme

MAX_LEVELS = 16

RB87 = species.RB87
LINE = "D1"
JE_D1, NU_D1_RB87, GAMMA_D1_MHZ, _, _ = RB87.line(LINE)
GAMMA_D1 = GAMMA_D1_MHZ * species.MHZ
K_D1_RB87 = 2 * np.pi * NU_D1_RB87 / constants.C_LIGHT

# Steck, Rubidium 87 D Line Data, Table 7: D1 far-detuned effective saturation
# intensity for pi-polarized light.  The Hanle geometry uses transverse light,
# so this is an intensity scale rather than a per-subtransition statement.
D1_ISAT_MW_CM2 = 4.4876

CELL_BUFFER = "Buffer gas cell"
CELL_PARAFFIN = "Paraffin coated cell"
SIGNAL_TRANSMISSION = "Transmission"
SIGNAL_NMOR = "NMOR rotation"

LIGHT_STATE_NORMALIZATION = "trace-one conditional light-region state"
DENSITY_CONVENTION = "local vapor density N, applied once"
ABSOLUTE_SCALE_STATUS = "external validation required"
MODEL_SCOPE_STATUS = "closed hyperfine pair; open-channel leakage and repump omitted"
SAMPLED_GRID_RTOL = 0.05

DOPPLER_OFF = "off: v = 0 class only"
DOPPLER_EXACT = "exact untruncated Maxwell average (pole–residue)"
# A pole-form B row is kept only if it reproduces the direct v = 0 solve to this
# relative accuracy; other rows are averaged on FALLBACK velocity grids.
POLE_GUARD_RTOL = 1e-6
POLE_LIMIT_SIGMAS = 1e4          # poles farther than this from the bulk are rounding modes
FALLBACK_SIGMAS = 6.0
FALLBACK_QUADRATURE_ERROR = 1e-6


def _transition_label(Fg, Fe):
    """Human-readable 87Rb D1 transition label for the Atomic dropdown."""
    return f"F={int(round(Fg))} → F'={int(round(Fe))}"


# 87Rb D1: F_g ∈ {1,2}, F'_e ∈ {1,2}; all four satisfy |F'_e − F_g| ≤ 1.
TRANSITION_CHOICES = tuple(_transition_label(g, e) for g in (1, 2) for e in (1, 2))


def _parse_transition(label):
    """(Fg, Fe) from a transition label, or None if it can't be parsed."""
    nums = re.findall(r"\d+", str(label))
    return (int(nums[0]), int(nums[1])) if len(nums) >= 2 else None


def _maxwell_fallback_grid(T, mass, k_vec, gamma_hwhm):
    """Uniform Maxwell velocity grid for B rows the pole form cannot certify.

    The optical response of one velocity class is analytic within the homogeneous
    half-width ``gamma_hwhm`` [rad/s] of the real Δ_eff axis, so the trapezoid
    error is ~exp(-2π·γ/h) for a detuning step h; h is chosen for
    FALLBACK_QUADRATURE_ERROR and the grid spans ±FALLBACK_SIGMAS.
    """
    sigma_v = math.sqrt(constants.KB * T / mass)
    step = 2 * math.pi * gamma_hwhm / math.log(1.0 / FALLBACK_QUADRATURE_ERROR) / k_vec
    half = math.ceil(FALLBACK_SIGMAS * sigma_v / step) * step
    v = np.arange(-half, half + 0.5 * step, step)
    wt = np.exp(-0.5 * (v / sigma_v) ** 2)
    return v, wt / wt.sum()


def _transition_strength(Fg, Fe):
    return species.line_strength(Fg, Fe, RB87.I, RB87.Jg, JE_D1)


def _transition_dipole(Fg, Fe):
    lam = constants.C_LIGHT / NU_D1_RB87
    d2_j = species.reduced_dipole_sq(GAMMA_D1, lam, RB87.Jg, JE_D1)
    return math.sqrt(max(_transition_strength(Fg, Fe), 0.0) * d2_j)


def _qwp_drive_weights(theta_deg):
    """Circular-basis drive weights, normalized to keep old linear scale.

    theta=0 gives weights {+1: 1, -1: 1}; theta=45 gives one circular component.
    """
    th = math.radians(float(theta_deg))
    c, s = math.cos(th), math.sin(th)
    ex = c * c + 1j * s * s
    ey = (1 - 1j) * s * c
    e_plus = (ex - 1j * ey) / math.sqrt(2.0)
    e_minus = (ex + 1j * ey) / math.sqrt(2.0)
    return {+1: math.sqrt(2.0) * e_plus, -1: math.sqrt(2.0) * e_minus}


def _scan_contains_zero(x):
    return (x.ndim == 1 and x.size >= 2 and np.all(np.isfinite(x))
            and np.all(np.diff(x) > 0.0) and x[0] <= 0.0 <= x[-1])


def _at_zero(x, y):
    """Interpolate only inside the sampled physical-field interval."""
    return (float(np.interp(0.0, x, y))
            if _scan_contains_zero(x) and np.all(np.isfinite(y)) else float("nan"))


def _central_feature(x, y):
    """Two half-height crossings about physical zero, with symmetric wings.

    The symmetric window is the largest one about zero contained in the scan.
    Interpolating its endpoints prevents an offset scan from choosing unequal
    physical wing distances. No Lorentzian shape or left/right symmetry is used.
    """
    result = dict(center=float("nan"), background=float("nan"), amplitude=float("nan"),
                  left=float("nan"), right=float("nan"), width=float("nan"),
                  status="physical zero outside scan")
    if not _scan_contains_zero(x) or not np.all(np.isfinite(y)):
        return result
    center = float(np.interp(0.0, x, y))
    radius = min(-float(x[0]), float(x[-1]))
    result.update(center=center, status="unresolved in scan")
    if radius <= 0.0:
        return result
    axis = np.unique(np.r_[x[(x > -radius) & (x < radius)], -radius, 0.0, radius])
    values = np.interp(axis, x, y)
    ic = int(np.searchsorted(axis, 0.0))
    background = float(0.5 * (values[0] + values[-1]))
    amplitude = center - background
    result.update(background=background, amplitude=amplitude)
    scale = max(abs(center), abs(background), np.finfo(float).tiny)
    if abs(amplitude) <= 32 * np.finfo(float).eps * scale:
        return result
    target = background + 0.5 * amplitude

    def crossing(step):
        i = ic
        while 0 <= i + step < axis.size:
            j = i + step
            if (values[j] - target) * np.sign(amplitude) <= 0.0:
                fraction = (target - values[i]) / (values[j] - values[i])
                return float(axis[i] + fraction * (axis[j] - axis[i]))
            i = j
        return float("nan")

    left, right = crossing(-1), crossing(+1)
    result.update(left=left, right=right)
    if np.isfinite(left) and np.isfinite(right) and right > left:
        result.update(width=right - left, status="sampled")
    return result


def _local_derivative(x, y, stride=1):
    """Derivative at physical zero from fourth-degree local interpolation.

    Five actual samples are used, including both sides of zero. The outer
    stencil uses stride=2, so consistency can be checked without a new solve.
    This also handles zero falling between samples after a field offset.
    """
    if not _scan_contains_zero(x) or not np.all(np.isfinite(y)):
        return float("nan")
    ic = int(np.argmin(np.abs(x)))
    indices = ic + stride * np.arange(-2, 3)
    if indices[0] < 0 or indices[-1] >= x.size:
        return float("nan")
    local_x = x[indices]
    if local_x[0] >= 0.0 or local_x[-1] <= 0.0:
        return float("nan")
    scale = float(np.max(np.abs(local_x)))
    z = local_x / scale
    values = y[indices] - np.interp(0.0, x, y)
    # Differentiate the five Lagrange basis polynomials at zero. These small
    # products avoid another linear-system solve in the readout path.
    derivative = 0.0
    for j in range(5):
        denominator = np.prod([z[j] - z[k] for k in range(5) if k != j])
        weight = sum(np.prod([-z[k] for k in range(5) if k not in (j, m)])
                     for m in range(5) if m != j) / denominator
        derivative += weight * values[j]
    return float(derivative / scale)


def _relative_change(a, b):
    if not np.isfinite(a) or not np.isfinite(b):
        return float("nan")
    return abs(a - b) / max(abs(a), abs(b), np.finfo(float).tiny)


def _zero_field_diagnostics(x, alpha, rotation):
    """Readout-only sampled consistency, distinct from external validation."""
    x = np.asarray(x, dtype=float)
    alpha = np.asarray(alpha, dtype=float)
    rotation = np.asarray(rotation, dtype=float)
    feature = _central_feature(x, alpha)
    coarse_indices = np.unique(np.r_[np.arange(0, x.size, 2), x.size - 1])
    coarse = _central_feature(x[coarse_indices], alpha[coarse_indices])
    width_change = _relative_change(feature["width"], coarse["width"])
    if not _scan_contains_zero(x):
        width_status = slope_status = "physical zero outside scan"
    else:
        width_status = ("unresolved in scan" if not np.isfinite(width_change) else
                        "sample-grid consistent" if width_change <= SAMPLED_GRID_RTOL
                        and feature["amplitude"] * coarse["amplitude"] > 0.0 else
                        "nonconverged on sampled grid")
        slope_status = "insufficient local samples"
    slope = _local_derivative(x, rotation)
    coarse_slope = _local_derivative(x, rotation, stride=2)
    slope_change = _relative_change(slope, coarse_slope)
    if np.isfinite(slope_change):
        slope_status = ("sample-grid consistent" if slope_change <= SAMPLED_GRID_RTOL
                        and slope * coarse_slope >= 0.0 else
                        "nonconverged on sampled grid")
    return dict(feature, width_status=width_status, width_relative_change=width_change,
                slope=slope, slope_status=slope_status, slope_relative_change=slope_change)


class SinglePairMagnetoReference(Scheme):
    """Historical single-pair solver for independent regression references.

    Never registered as an application scheme; public MagnetoScheme below
    uses the multihyperfine Budker engine for both scan axes.
    """
    cluster = "C - Magneto-optics"
    presets_group = "Default"
    cache_version = "polarized-two-region-v7-budker-sensitivity"
    defaults_version = "polarized-two-region-v4"
    supports_headless_observables = True

    _REGIME_EMOJI = {
        "EIT dip": "🕳️",
        "EIA peak": "⛰️",
        "Buffer Hanle": "🌫️",
        "Buffer LCA": "✖️",
        "NMOR": "🧭",
    }

    _DEF = {
        "hanle": dict(signal=SIGNAL_TRANSMISSION, cell=CELL_PARAFFIN, Fg=2, Fe=1,
                      intensity=0.8, qwp=0.0, bmax=2.0,
                      title="Hanle transmission",
                      desc="87Rb D1 Hanle transmission with probe polarization."),
        "eia": dict(signal=SIGNAL_TRANSMISSION, cell=CELL_PARAFFIN, Fg=2, Fe=1,
                    intensity=0.8, qwp=45.0, bmax=2.0,
                    title="Polarization-switched EIA/MIA",
                    desc="87Rb D1 Hanle signal switched toward absorption by QWP angle."),
        "nmor": dict(signal=SIGNAL_NMOR, cell=CELL_PARAFFIN, Fg=2, Fe=1,
                     intensity=0.7, qwp=0.0, bmax=2.0,
                     title="Nonlinear magneto-optical rotation (NMOR)",
                     desc="87Rb D1 polarization-rotation readout near zero magnetic field."),
    }

    def __init__(self, mode=None):
        self.mode = mode
        if mode is None:
            self.name = "magneto"
            self.title = "Magneto-optics (Hanle/MOR)"
            self.caption = ("87Rb D1 magnetic-field scans and an optional Rb D1/D2 "
                            "multihyperfine detuning-sensitivity calculation. The field "
                            "scans resolve the ground-state "
                            "Zeeman manifold. Two distinct effects: the Hanle effect "
                            "(zero-field transmission resonance from ground-state "
                            "coherence) and magneto-optical rotation (MOR/NMOR, "
                            "polarization-plane rotation). Probe QWP angle, residual "
                            "transverse field, and cell relaxation decide whether the "
                            "Hanle feature is EIT-like transparency or EIA/MIA-like "
                            "absorption.")
        else:
            self.name = mode
            self.title = self._DEF[mode]["title"]
            self.caption = self._DEF[mode]["desc"]

    def _regime_options(self):
        """Flat, ordered list of emoji-tagged regime labels across both readouts.

        Used both for the segmented `regime` control's choices and as the keys of
        `recommended_defaults`, so a selection maps straight onto a parameter set.
        """
        opts = []
        for regimes in self._REGIMES.values():
            for label in regimes:
                opts.append(f"{self._REGIME_EMOJI.get(label, '•')} {label}")
        return opts

    def param_schema(self):
        d = self._DEF[self.mode or "hanle"]
        specs = []
        if self.mode is None:
            options = self._regime_options()
            specs.append(ParamSpec(
                "regime", "Regime", "Regime", options[0],
                choices=tuple(options), control="segmented",
                applies_defaults=True, recompute=False,
                help="Pick a ready-made magneto-optics regime; selecting one loads "
                     "its full parameter set (including the Signal readout). Same "
                     "in-panel regime selector the other schemes use."))
        specs += [
            ParamSpec("signal_type", "Signal", "Readout", d["signal"],
                      choices=(SIGNAL_TRANSMISSION, SIGNAL_NMOR),
                      control="segmented", recompute=False,
                      help="Transmission shows Hanle absorption/transparency; NMOR shows rotation."),
            ParamSpec("cell_type", "Cell type", "Cell & beams", d["cell"],
                      choices=(CELL_PARAFFIN, CELL_BUFFER), control="segmented",
                      help="Choose the relaxation model and the cell-specific sliders."),
            ParamSpec("transition", "Transition (F → F')", "Atomic",
                      _transition_label(d["Fg"], d["Fe"]), choices=TRANSITION_CHOICES,
                      help="87Rb D1 hyperfine transition F_g → F'_e "
                           "(all four satisfy |F'_e − F_g| ≤ 1)."),
            ParamSpec("intensity_mw_cm2", "Beam intensity", "Fields", d["intensity"],
                      0.01, 20.0, 0.05, "mW/cm²",
                      help="Uniform incident optical intensity before the vapor cell."),
            ParamSpec("qwp_deg", "Probe polarization (QWP angle)", "Fields", d["qwp"],
                      0.0, 45.0, 0.5, "deg", endpoints=("linear", "circular"),
                      help="0 deg is linear; 45 deg is circular. Absolute handedness is hidden."),
            ParamSpec("b_max_ut", "B scan ±", "Fields", d["bmax"],
                      0.02, 500.0, 0.01, "µT",
                      help="Longitudinal magnetic-field scan range around zero."),
            ParamSpec("b_offset_ut", "B zero offset", "Fields", 0.0,
                      -50.0, 50.0, 0.01, "µT", advanced=True,
                      help="Longitudinal bias added to the physical field while "
                           "the displayed scan stays centered at the commanded "
                           "zero. Use it for residual shielding or coil offsets.",
                      advanced_group="Systematics"),
            ParamSpec("laser_detuning_mhz", "Laser detuning", "Detunings", 0.0,
                      -1000.0, 1000.0, 5.0, "MHz", advanced=True,
                      help="Detuning from the selected 87Rb D1 hyperfine transition.",
                      advanced_group="Systematics"),
            ParamSpec("temp_c", "Temperature", "Cell & beams", 25.0,
                      20.0, 120.0, 1.0, "°C",
                      help="Sets Doppler width and the vapor-pressure density when selected."),
            ParamSpec("cell_mm", "Cell length", "Cell & beams", 10.0,
                      0.5, 200.0, 0.5, "mm", recompute=False,
                      advanced=True, advanced_group="Cell geometry",
                      help="Readout path length for the Beer-Lambert conversion."),
            ParamSpec("ne_pressure_torr", "Ne buffer pressure", "Cell & beams", 10.0,
                      0.0, 200.0, 1.0, "Torr",
                      visible_if={"cell_type": CELL_BUFFER},
                      help="Simple Ne pressure broadening for buffer-gas mode."),
            # The headline cell-quality knob (paraffin): the anti-relaxation
            # coating's ground-coherence lifetime. Kept primary, in the cell group.
            ParamSpec("wall_coherence_ms", "Wall coherence lifetime", "Cell & beams", 3.0,
                      0.05, 200.0, 0.05, "ms",
                      visible_if={"cell_type": CELL_PARAFFIN},
                      help="Anti-relaxation coating quality (ground-coherence lifetime). "
                           "The primary paraffin cell knob."),
            # Microscopic relaxation rates below are consequences of the cell /
            # beam (buffer pressure, coating, beam size), not independent dials.
            # Kept as advanced overrides; the cell-property knobs set the regime.
            ParamSpec("buffer_ground_relax_khz", "Buffer ground relaxation", "Atomic", 3.0,
                      0.1, 500.0, 0.5, "kHz", advanced=True,
                      visible_if={"cell_type": CELL_BUFFER},
                      advanced_group="Relaxation overrides",
                      help="Effective transit/reservoir reload into an isotropic "
                           "ground state; exposed here as an override."),
            ParamSpec("collisional_depol_khz", "Ground pure dephasing", "Atomic", 0.5,
                      0.0, 200.0, 0.5, "kHz", advanced=True,
                      visible_if={"cell_type": CELL_BUFFER},
                      advanced_group="Relaxation overrides",
                      help="Ground-state phase randomization without redistribution "
                           "of Zeeman populations; exposed as an override. This rate "
                           "does not describe isotropic spin destruction."),
            ParamSpec("transit_relax_khz", "Transit relaxation", "Atomic", 80.0,
                      0.1, 1000.0, 1.0, "kHz", advanced=True,
                      visible_if={"cell_type": CELL_PARAFFIN},
                      advanced_group="Relaxation overrides",
                      help="Atoms leaving the illuminated region; sets the broad "
                           "pedestal. A beam-size / temperature consequence."),
            ParamSpec("dark_return_khz", "Dark-region return rate", "Atomic", 1.0,
                      0.01, 100.0, 0.01, "kHz", advanced=True,
                      visible_if={"cell_type": CELL_PARAFFIN},
                      advanced_group="Relaxation overrides",
                      help="Rate for wall-preserved atoms returning from the dark region."),
            # Residual transverse field: a systematic (shielding leftover / shim-coil
            # compensation), not a primary scan knob. Advanced — presets carry the
            # values that switch EIT<->EIA and enable circular-light LCA.
            ParamSpec("residual_transverse_b_ut", "Residual transverse B", "Fields", 0.05,
                      0.0, 5.0, 0.005, "µT", advanced=True,
                      help="Weak transverse field. Enables MIA/MIT-like switching and "
                           "the circular-light level-crossing absorption (LCA).",
                      advanced_group="Systematics"),
            ParamSpec("transverse_field_angle_deg", "Transverse B angle", "Fields", 0.0,
                      0.0, 180.0, 1.0, "deg", advanced=True,
                      help="Azimuth of the residual transverse magnetic field.",
                      advanced_group="Systematics"),
            ParamSpec("line_strength", "Line-strength factor", "Detection & scaling", 1.0,
                      0.01, 2.0, 0.01, "", recompute=False, advanced=True,
                      help="Calibration knob for effective optical depth (not a "
                           "physical experimental variable).",
                      advanced_group="Calibration"),
            ParamSpec("show_shot_noise", "Show shot noise", "Detection & scaling", True,
                      control="checkbox", recompute=False,
                      help="Overlay a simulated photon-shot-noise trace and its ±1σ band "
                           "in the plotted y-axis units. Transmission fluctuates around "
                           "the wing background; rotation fluctuates around zero."),
            ParamSpec("beam_diameter_mm", "Effective beam diameter", "Detection & scaling",
                      DEFAULT_BEAM_DIAMETER_MM, 0.1, 20.0, 0.1, "mm",
                      recompute=False, advanced=True, advanced_group="Shot noise",
                      help="Uniform effective illuminated/collected diameter. Incident "
                           "power = beam intensity × π(d/2)². This sets the photon budget "
                           "only; atomic transit rates remain separately specified."),
            ParamSpec("detection_efficiency", "Detection efficiency", "Detection & scaling",
                      DEFAULT_DETECTION_EFFICIENCY, 0.01, 1.0, 0.01, "",
                      recompute=False, advanced=True, advanced_group="Shot noise",
                      help="After-cell optical throughput × detector quantum efficiency. "
                           "Cell transmission is calculated and applied separately once."),
            ParamSpec("noise_bandwidth_hz", "Noise bandwidth (ENBW)", "Detection & scaling",
                      DEFAULT_NOISE_BANDWIDTH_HZ, 0.1, 1e6, 0.1, "Hz",
                      recompute=False, advanced=True, advanced_group="Shot noise",
                      help="One-sided electrical equivalent noise bandwidth. Shot-noise "
                           "RMS scales as √ENBW; 1 kHz is an assumed acquisition setting, "
                           "not a measured detector bandwidth."),
            ParamSpec("doppler", "Doppler averaging", "Numerics", "on",
                      choices=("on", "off"), advanced=True,
                      help="On: exact Maxwell-Boltzmann average of the light-region "
                           "optical Doppler shift (pole–residue form, no velocity "
                           "grid). Off: the v = 0 class only.",
                      advanced_group="Numerics"),
            ParamSpec("scan_points", "B scan points", "Numerics", 121,
                      51, 401, 10, "", advanced=True, advanced_group="Numerics"),
        ]
        from . import magneto_budker
        common = {"intensity_mw_cm2", "temp_c", "beam_diameter_mm",
                  "detection_efficiency", "noise_bandwidth_hz"}
        specs = [replace(s, visible_if={**(s.visible_if or {}),
                                       "calculation": magneto_budker.LEGACY})
                 if s.name not in common else s for s in specs]
        specs.insert(0, ParamSpec(
            "calculation", "Calculation", "Model", magneto_budker.LEGACY,
            choices=(magneto_budker.LEGACY, magneto_budker.BUDKER),
            help="Choose a magnetic-field scan or a multihyperfine optical-detuning "
                 "scan of photon-shot-noise-limited magnetic sensitivity."))
        specs += [ParamSpec(
            "magnetic_signal_ft", "Small magnetic signal (RMS)", "Detection & scaling",
            1.0, 0.0, 1e6, 1.0, "fT", recompute=False, advanced=True,
            advanced_group="Shot noise",
            help="Small-signal amplitude for local, quasistatic magnetic SNR. "
                 "Sinusoidal peak = √2 × RMS must fit inside the tested local "
                 "field interval. The atomic frequency response is not included.")]
        return specs + magneto_budker.param_schema()

    # Regime parameter sets, grouped by readout. The merged entry renders these as
    # readout-contextual one-click default buttons (recommended_defaults).
    _REGIMES = {
        SIGNAL_TRANSMISSION: {
            "EIT dip": dict(
                signal_type=SIGNAL_TRANSMISSION, cell_type=CELL_PARAFFIN,
                transition=_transition_label(2, 1), intensity_mw_cm2=0.8, qwp_deg=0.0,
                b_max_ut=2.0, transit_relax_khz=80.0, dark_return_khz=1.0,
                wall_coherence_ms=3.0, residual_transverse_b_ut=0.05,
                transverse_field_angle_deg=0.0, doppler="on"),
            "EIA peak": dict(
                signal_type=SIGNAL_TRANSMISSION, cell_type=CELL_PARAFFIN,
                transition=_transition_label(2, 1), intensity_mw_cm2=0.8, qwp_deg=45.0,
                b_max_ut=2.0, transit_relax_khz=80.0, dark_return_khz=1.0,
                wall_coherence_ms=3.0, residual_transverse_b_ut=0.08,
                transverse_field_angle_deg=0.0, doppler="on"),
            "Buffer Hanle": dict(
                signal_type=SIGNAL_TRANSMISSION, cell_type=CELL_BUFFER,
                transition=_transition_label(2, 1), intensity_mw_cm2=0.8, qwp_deg=0.0,
                b_max_ut=80.0, ne_pressure_torr=20.0,
                buffer_ground_relax_khz=20.0, collisional_depol_khz=2.0,
                residual_transverse_b_ut=0.0, doppler="on"),
            "Buffer LCA": dict(
                signal_type=SIGNAL_TRANSMISSION, cell_type=CELL_BUFFER,
                transition=_transition_label(2, 2), intensity_mw_cm2=0.2, qwp_deg=45.0,
                b_max_ut=2.0, ne_pressure_torr=20.0,
                buffer_ground_relax_khz=5.0, collisional_depol_khz=0.5,
                residual_transverse_b_ut=0.03, transverse_field_angle_deg=90.0,
                doppler="on"),
        },
        SIGNAL_NMOR: {
            "NMOR": dict(
                signal_type=SIGNAL_NMOR, cell_type=CELL_PARAFFIN,
                transition=_transition_label(2, 1), intensity_mw_cm2=0.7, qwp_deg=0.0,
                b_max_ut=2.0, transit_relax_khz=50.0, dark_return_khz=1.0,
                wall_coherence_ms=5.0, residual_transverse_b_ut=0.02,
                transverse_field_angle_deg=0.0, doppler="on"),
        },
    }

    def presets(self):
        # The merged Hanle/EIA/NMOR entry uses readout-contextual default buttons
        # (recommended_defaults) instead of a flat preset wall; single-mode aliases
        # keep one labelled preset for tests / direct use.
        if self.mode is None:
            from . import magneto_budker
            reference = dict(
                calculation=magneto_budker.BUDKER, budker_isotope="Rb85",
                budker_ground_center="upper", temp_c=20.0,
                budker_density_mode="manual", budker_density_cm3=4.5e9,
                budker_cell_mm=100.0, beam_diameter_mm=3.0,
                budker_bias_ut=0.0, budker_velocity_mixing="complete",
                budker_doppler="on", noise_bandwidth_hz=1.0,
                detection_efficiency=0.8,
            )
            return [
                Preset("Magnetic-field scan", dict(calculation=magneto_budker.LEGACY)),
                Preset("Budker D1 conditions (model)",
                       dict(reference, budker_line="D1", intensity_mw_cm2=1.0),
                       help="Paper optical conditions with assumed effective relaxation and "
                            "detection efficiency; conditional model, not a calibrated reproduction."),
                Preset("Budker D2 conditions (model)",
                       dict(reference, budker_line="D2", intensity_mw_cm2=4.5),
                       help="Paper optical conditions; check weak-excitation validity before "
                            "interpreting this higher-intensity prediction."),
            ]
        d = self._DEF[self.mode]
        return [Preset(f"D1 {self.mode.upper()} default", icon=self.mode.upper(),
                       values=dict(signal_type=d["signal"], cell_type=d["cell"],
                                   transition=_transition_label(d["Fg"], d["Fe"]),
                                   intensity_mw_cm2=d["intensity"], qwp_deg=d["qwp"],
                                   b_max_ut=d["bmax"], doppler="on"))]

    def recommended_defaults(self, params):
        # The merged entry exposes a single segmented "Regime" control
        # (applies_defaults); selecting a regime loads its full parameter set,
        # including the Signal readout. Return every regime (keyed by the same
        # emoji-tagged label as the control's choices) so any selection resolves.
        if self.mode is not None:
            return None
        out = {}
        for regimes in self._REGIMES.values():
            for label, values in regimes.items():
                out[f"{self._REGIME_EMOJI.get(label, '•')} {label}"] = dict(values)
        return out

    def info(self):
        return (
            "**Hanle / magneto-optical rotation and magnetic sensitivity.** "
            "The magnetic-field scan uses the 87Rb D1 single-pair model. "
            "The optional Budker detuning calculation includes multiple Rb "
            "hyperfine manifolds and velocity redistribution in a coated cell, "
            "using a conditional weak-optical-excitation effective model. "
            "Its photon-shot-noise magnetic sensitivity is noise ASD divided by "
            "the local rotation slope, distinct from scan-peak amplitude SNR. "
            "These are two "
            "distinct effects of the ground-state Zeeman manifold near zero field: "
            "the Hanle effect is a transmission resonance from ground-state "
            "coherence, while magneto-optical rotation (MOR/NMOR) is a rotation of "
            "the probe polarization plane. The transmission readout reports "
            "whether the zero-field Hanle feature is EIT-like transparency, EIA/MIA-like "
            "absorption, or a crossover.  Paraffin-coated cells use a two-region "
            "light/dark OBE exchange model so wall-preserved ground coherence can "
            "create a narrow Ramsey feature on a broad transit pedestal.\n\n"
            "**References**\n"
            "- D. Budker et al., sensitive magnetometry based on nonlinear "
            "magneto-optical rotation, *Phys. Rev. A* **62**, 043403 (2000), "
            "https://doi.org/10.1103/PhysRevA.62.043403.\n"
            "- D. Budker et al., resonant nonlinear magneto-optical effects in atoms, "
            "*Rev. Mod. Phys.* **74**, 1153 (2002), https://arxiv.org/abs/physics/0203077.\n"
            "- Yu. P. Malakyan et al., velocity mixing in low-power FM-NMOR, "
            "*Phys. Rev. A* **69**, 013817 (2004), https://arxiv.org/abs/physics/0309061.\n"
            "- S. M. Rochester and D. Budker, NMOR in optically thick media, "
            "*J. Mod. Opt.* **49**, 2543 (2002), https://arxiv.org/abs/physics/0202071.\n"
            "- M. Fleischhauer et al., quantum limit with ac Stark shifts, "
            "*Phys. Rev. A* **62**, 013808 (2000), https://arxiv.org/abs/quant-ph/0001072.\n"
            "- D. A. Steck, *Rubidium 87 D Line Data*, http://steck.us/alkalidata.\n"
            "- H. J. Lee and H. S. Moon, magnetic-field-induced absorption in a "
            "paraffin-coated rubidium vapor cell, *JOSA B* **30**, 2301 (2013).\n"
            "- H. S. Moon and H. J. Kim, Ramsey EIA to MIT transformation in a "
            "paraffin-coated Rb vapor cell, *Optics Express* **22**, 18604 (2014). "
            "Its D2 F=2 → F'=3 transition is a qualitative coating reference, "
            "outside the current D1 transition set.\n"
            "- Y.-J. Yu, H. J. Lee, and H. S. Moon, level-crossing absorption "
            "with narrow spectral width in Rb vapor with buffer gas, "
            "*Phys. Rev. A* **81**, 023416 (2010)."
        )

    def compute(self, params):
        from . import magneto_budker
        if params.get("calculation") == magneto_budker.BUDKER:
            return magneto_budker.compute(params)
        # The UI exposes a single transition dropdown; tests/aliases may still
        # pass explicit Fg/Fe, which take priority.
        if "Fg" in params and "Fe" in params:
            Fg, Fe = int(round(params["Fg"])), int(round(params["Fe"]))
        else:
            parsed = _parse_transition(params.get("transition", _transition_label(2, 1)))
            Fg, Fe = parsed if parsed else (2, 1)
        Fg_allowed = set(int(round(f)) for f in species.f_values(RB87.I, RB87.Jg))
        Fe_allowed = set(int(round(f)) for f in species.f_values(RB87.I, JE_D1))
        strength = _transition_strength(Fg, Fe) if (Fg in Fg_allowed and Fe in Fe_allowed) else 0.0
        valid = (Fg in Fg_allowed and Fe in Fe_allowed and abs(Fe - Fg) <= 1
                 and strength > 1e-12 and (2 * Fg + 1) + (2 * Fe + 1) <= MAX_LEVELS)

        gFg = (species.lande_gf(species.GJ_S12_LEGACY, RB87.g_I, RB87.I, RB87.Jg, Fg)
               if Fg in Fg_allowed else 0.0)
        gFe = (species.lande_gf(RB87.gJ_P12, RB87.g_I, RB87.I, JE_D1, Fe)
               if Fe in Fe_allowed else 0.0)

        cell_type = params.get("cell_type", CELL_PARAFFIN)
        buffer_gamma = (constants.neon_buffer_broadening(
            params.get("ne_pressure_torr", 0.0)) if cell_type == CELL_BUFFER else 0.0)
        gamma_homogeneous = GAMMA_D1 + buffer_gamma
        gamma_elastic = buffer_gamma / 2.0
        gamma_transit = gamma_depol = 0.0

        if cell_type == CELL_BUFFER:
            gamma_transit = 2 * np.pi * params["buffer_ground_relax_khz"] * 1e3
            gamma_depol = 2 * np.pi * params["collisional_depol_khz"] * 1e3
            atom = (zeeman.zeeman_manifold(
                        Fg, Fe, gamma=GAMMA_D1, gamma_gg=gamma_depol,
                        transit_rate=gamma_transit,
                        optical_dephasing=gamma_elastic)
                    if valid else None)
            atom_dark = None
            gamma_out = gamma_in = gamma_wall = 0.0
        else:
            gamma_light = 2 * np.pi * params["transit_relax_khz"] * 1e3
            gamma_out = gamma_light
            gamma_in = 2 * np.pi * params["dark_return_khz"] * 1e3
            gamma_wall = 1.0 / max(params["wall_coherence_ms"] * 1e-3, 1e-9)
            atom = (zeeman.zeeman_manifold(Fg, Fe, gamma=GAMMA_D1, gamma_gg=gamma_light,
                                           transit_rate=0.0)
                    if valid else None)
            atom_dark = (zeeman.zeeman_manifold(Fg, Fe, gamma=GAMMA_D1, gamma_gg=gamma_wall,
                                                transit_rate=0.0)
                         if valid else None)

        b_ut = np.linspace(-params["b_max_ut"], params["b_max_ut"],
                           int(params["scan_points"]))
        b_offset_ut = float(params.get("b_offset_ut", 0.0))
        b_physical_ut = b_ut + b_offset_ut
        b_z = b_physical_ut * 1e-6
        # A real vapor cell always has a small residual transverse field (imperfect
        # shielding); it enters the Zeeman Hamiltonian exactly like any field. It is
        # essential for circular-light level-crossing absorption: σ± light orients
        # the ground state along the beam (z), an eigenstate of the longitudinal B
        # scan, so only a transverse component makes that orientation precess and
        # produce a B=0 resonance. Buffer and paraffin cells share this physics.
        b_perp = params.get("residual_transverse_b_ut", 0.0) * 1e-6
        phi = math.radians(params.get("transverse_field_angle_deg", 0.0))
        b_x = b_perp * math.cos(phi)
        b_y = b_perp * math.sin(phi)

        T = params["temp_c"] + 273.15
        N = species.number_density(RB87, T)
        intensity = max(float(params["intensity_mw_cm2"]), 0.0)
        Om = GAMMA_D1 * math.sqrt(intensity / (2 * D1_ISAT_MW_CM2)) if intensity > 0 else 0.0
        drive = _qwp_drive_weights(params["qwp_deg"])
        dL = 2 * np.pi * params["laser_detuning_mhz"] * 1e6

        chi_probe = np.zeros(b_ut.size, dtype=complex)
        chi_p = np.zeros(b_ut.size, dtype=complex)
        chi_m = np.zeros(b_ut.size, dtype=complex)
        doppler_on = params["doppler"] == "on"
        sigma_v = math.sqrt(constants.KB * T / RB87.mass)
        doppler = dict(method=DOPPLER_EXACT if doppler_on else DOPPLER_OFF,
                       guard_max_relative=None, fallback_rows=0, fallback_classes=0,
                       min_pole_distance_mhz=None)
        if cell_type == CELL_PARAFFIN:
            exchange_rate = gamma_in + gamma_out
            light_region_occupation = (gamma_in / exchange_rate
                                       if exchange_rate > 0 else 1.0)
        else:
            light_region_occupation = 1.0
        light_state_trace = 1.0
        if valid and Om > 0:
            # H₀(bz) = H_xy + bz·H_z is affine in the longitudinal field, and
            # comm_super is linear, so L0(bz) = C_xy + bz·C_z. Build the constant
            # pieces once (no per-B-point kron) and form the whole B-scan stack of
            # Liouvillians by a scalar broadcast. Bit-for-bit the same matrices as
            # the old loop.
            H_xy = self._hamiltonian(atom, (b_x, b_y, 0.0), gFg, gFe, Om, drive)
            H_z = self._hamiltonian(atom, (0.0, 0.0, 1.0), gFg, gFe, 0.0,
                                    {+1: 0.0, -1: 0.0})
            C_xy = core.build_liouvillian(H_xy, atom)
            C_z = core.comm_super(H_z)                       # drive-free → reusable
            L0_all = C_xy[None, :, :] + b_z[:, None, None] * C_z[None, :, :]
            Cd_xy = Ld_all = None
            if cell_type == CELL_PARAFFIN:
                Hd_xy = self._hamiltonian(atom_dark, (b_x, b_y, 0.0), gFg, gFe,
                                          0.0, {+1: 0.0, -1: 0.0})
                Cd_xy = core.build_liouvillian(Hd_xy, atom_dark)
                Ld_all = Cd_xy[None, :, :] + b_z[:, None, None] * C_z[None, :, :]

            def solve(rows, deff):
                """(conditional light-region ρ, occupation) of B rows at detunings deff."""
                if cell_type == CELL_PARAFFIN:
                    return self._steady_state_two_region(
                        L0_all[rows], Ld_all[rows], deff, atom.S_v, atom.n_levels,
                        gamma_out, gamma_in)
                return self._steady_state_buffer(
                    L0_all[rows], deff, atom.S_v, atom.n_levels), None

            if doppler_on:
                def grid_average(rows):
                    v, wt = _maxwell_fallback_grid(
                        T, RB87.mass, K_D1_RB87, 0.5 * gamma_homogeneous)
                    out_p = np.empty(rows.size, dtype=complex)
                    out_m = np.empty(rows.size, dtype=complex)
                    for part in np.array_split(np.arange(rows.size), -(-rows.size // 8)):
                        rho_v, _ = solve(rows[part], dL - K_D1_RB87 * v)
                        cp_v, cm_v, _ = self._coherences(atom, rho_v, Om, drive)
                        out_p[part] = cp_v @ wt
                        out_m[part] = cm_v @ wt
                    return out_p, out_m, int(v.size)

                chi_p, chi_m, occupation_mean = self._exact_doppler_average(
                    atom, C_xy, C_z, Cd_xy, b_z, gamma_out, gamma_in, Om,
                    center=dL, sigma=K_D1_RB87 * sigma_v,
                    grid_average=grid_average, ledger=doppler)
                if Cd_xy is not None:
                    light_region_occupation = occupation_mean
                chi_probe = (np.conj(drive.get(+1, 0.0)) * chi_p
                             + np.conj(drive.get(-1, 0.0)) * chi_m)
            else:
                # The v = 0 class only: the same kernel call and sums as before.
                rho, occupation = solve(slice(None), np.array([dL]))
                if occupation is not None:
                    light_region_occupation = float(np.mean(occupation))
                    light_state_trace = float(np.mean(np.trace(
                        rho, axis1=-2, axis2=-1).real))
                cp, cm, cprobe = self._coherences(atom, rho, Om, drive)   # (nB, 1)
                wt = np.ones(1)
                chi_p = (cp * wt[None, :]).sum(axis=1)
                chi_m = (cm * wt[None, :]).sum(axis=1)
                chi_probe = (cprobe * wt[None, :]).sum(axis=1)

        return dict(
            line=LINE, isotope="87Rb", cell_type=cell_type, b_ut=b_ut,
            b_physical_ut=b_physical_ut, b_offset_ut=b_offset_ut,
            chi_probe=chi_probe, chi_p=chi_p, chi_m=chi_m, N=N, N_eff=N,
            T=T, valid=valid, Fg=Fg, Fe=Fe, gFg=gFg, gFe=gFe,
            gamma=gamma_homogeneous, gamma_natural=GAMMA_D1,
            buffer_gamma=buffer_gamma, gamma_elastic=gamma_elastic,
            gamma_transit=gamma_transit, gamma_depol=gamma_depol,
            gamma_out=gamma_out, gamma_in=gamma_in, gamma_wall=gamma_wall,
            k_vec=K_D1_RB87, dipole=_transition_dipole(Fg, Fe) if valid else 0.0,
            strength=strength, omega_rabi=Om, isat=D1_ISAT_MW_CM2,
            doppler_scale=1.0, doppler_method=doppler["method"],
            doppler_sigma_mhz=K_D1_RB87 * sigma_v / (2 * np.pi) / 1e6,
            doppler_guard_max_relative=doppler["guard_max_relative"],
            doppler_fallback_rows=doppler["fallback_rows"],
            doppler_fallback_classes=doppler["fallback_classes"],
            doppler_min_pole_distance_mhz=doppler["min_pole_distance_mhz"],
            qwp_deg=params["qwp_deg"], b_perp_ut=b_perp * 1e6,
            light_region_occupation=light_region_occupation,
            light_state_trace=light_state_trace,
            light_state_normalization=LIGHT_STATE_NORMALIZATION,
            density_convention=DENSITY_CONVENTION,
            absolute_scale_status=ABSOLUTE_SCALE_STATUS,
        )

    @staticmethod
    def _hamiltonian(atom, b_vec_t, gFg, gFe, Om, drive):
        n = atom.n_levels
        ng = len(atom.ground)
        H = np.zeros((n, n), dtype=complex)
        bx, by, bz = b_vec_t
        Fgx, Fgy, Fgz = zeeman.angular_momentum_matrices((ng - 1) / 2)
        Fex, Fey, Fez = zeeman.angular_momentum_matrices((len(atom.excited) - 1) / 2)
        ground_block = constants.MU_B * gFg / constants.HBAR * (bx * Fgx + by * Fgy + bz * Fgz)
        excited_block = constants.MU_B * gFe / constants.HBAR * (bx * Fex + by * Fey + bz * Fez)
        H[:ng, :ng] = ground_block
        H[ng:, ng:] = excited_block
        for q in (+1, -1):
            Omq = Om * drive.get(q, 0.0)
            if abs(Omq) < 1e-30:
                continue
            for gi, ei, cg in atom.couplings[q]:
                H[ei, gi] += Omq * cg / 2
                H[gi, ei] += np.conj(Omq) * cg / 2
        return H

    @staticmethod
    def _conditional_light_state(rho_light):
        """Separate light-region occupation from its trace-one optical state."""
        occupation = np.trace(rho_light, axis1=-2, axis2=-1)
        if np.any(np.abs(occupation) <= np.finfo(float).eps):
            raise RuntimeError("two-region light-state occupation is zero")
        rho_conditional = rho_light / occupation[..., None, None]
        return rho_conditional, occupation.real

    @staticmethod
    def _steady_state_two_region(L_light0, L_dark0, deff, S_v, n_levels,
                                 gamma_out, gamma_in):
        """Two-region (light/dark) steady state for a stack of B-field points.

        Returns ``(rho_light_conditional, occupation)``. The coupled light/dark
        solve is globally normalized, so occupation is the trace of its raw
        light block. Susceptibility uses only the trace-one conditional state
        because ``N`` is the local vapor density; occupation is diagnostic only.

        L_light0, L_dark0: (nB, M, M) per-B Liouvillians. Only the light region
        carries the optical Δ_eff (Doppler) shift; the dark region is field-
        shifted but unlit. Solves the (B × velocity) batch in memory-bounded
        chunks. Returns the in-beam ρ_light, shape (nB, nv, n_levels, n_levels).
        """
        if kernels.available():
            with core.blas_single_thread():
                rho_l = kernels.magneto_two_region_grid(
                    np.ascontiguousarray(L_light0),
                    np.ascontiguousarray(L_dark0),
                    np.ascontiguousarray(deff, dtype=float),
                    np.ascontiguousarray(S_v),
                    float(gamma_out), float(gamma_in), n_levels)
            return SinglePairMagnetoReference._conditional_light_state(rho_l)
        nB = L_light0.shape[0]
        nv = deff.size
        M = n_levels * n_levels
        eye = core._eye(M)
        rho_l = np.empty((nB, nv, n_levels, n_levels), dtype=complex)
        rows = max(1, int(9.0e6 // max(nv * (2 * M) ** 2, 1)))
        for b0 in range(0, nB, rows):
            sl = slice(b0, b0 + rows)
            Ll = (L_light0[sl][:, None, :, :]
                  - deff[None, :, None, None] * S_v[None, None, :, :])   # (nb,nv,M,M)
            Ld = np.broadcast_to(L_dark0[sl][:, None, :, :], Ll.shape)   # unlit: no Δ_eff
            nb = Ll.shape[0]
            A = np.zeros((nb, nv, 2 * M, 2 * M), dtype=complex)
            A[..., :M, :M] = Ll - gamma_out * eye
            A[..., :M, M:] = gamma_in * eye
            A[..., M:, :M] = gamma_out * eye
            A[..., M:, M:] = Ld - gamma_in * eye
            A[..., 0, :] = 0
            for state in range(n_levels):
                idx = state * n_levels + state
                A[..., 0, idx] = 1
                A[..., 0, M + idx] = 1
            rhs = np.zeros((nb, nv, 2 * M, 1), dtype=complex)
            rhs[..., 0, 0] = 1
            sol = np.linalg.solve(A, rhs)[..., 0]
            rho_l[sl] = sol[..., :M].reshape(nb, nv, n_levels, n_levels)
        return SinglePairMagnetoReference._conditional_light_state(rho_l)

    @staticmethod
    def _steady_state_buffer(L0_all, deff, S_v, n_levels):
        """Single-region steady state for a stack of B-field points.

        L0_all: (nB, M, M). Solves the (B × velocity) batch in memory-bounded
        chunks. Returns ρ, shape (nB, nv, n_levels, n_levels).
        """
        if kernels.available():
            with core.blas_single_thread():
                return kernels.magneto_buffer_grid(
                    np.ascontiguousarray(L0_all),
                    np.ascontiguousarray(deff, dtype=float),
                    np.ascontiguousarray(S_v), n_levels)
        nB = L0_all.shape[0]
        nv = deff.size
        M = n_levels * n_levels
        rho = np.empty((nB, nv, n_levels, n_levels), dtype=complex)
        rows = max(1, int(1.2e7 // max(nv * M * M, 1)))
        for b0 in range(0, nB, rows):
            sl = slice(b0, b0 + rows)
            A = (L0_all[sl][:, None, :, :]
                 - deff[None, :, None, None] * S_v[None, None, :, :])
            rho[sl] = core.steady_state_from_liouvillian(A, n_levels)
        return rho

    @staticmethod
    def _exact_doppler_average(atom, C_xy, C_z, Cd_xy, b_z, gamma_out, gamma_in, Om,
                               center, sigma, grid_average, ledger):
        """Untruncated Maxwell average of the light-region σ± coherences per B row.

        Only the light region carries the optical Doppler shift, so with the trace
        row inserted every B row is ``(A00 + b_z·A01 − Δ_eff·S) x = e₀`` in the real
        Hermitian-generator basis, and each coherence is an exact rational function
        of Δ_eff (`pole_doppler.ShiftInvertSystem`).  Its average over
        Δ_eff ~ Normal(center, sigma²) is closed form (Faddeeva): no velocity grid,
        no cutoff.  Light and dark Liouvillians conserve their own traces, so the
        two-region light occupation does not depend on Δ_eff and dividing by it
        keeps the trace-one conditional state of the grid path.

        A row is kept if its pole form reproduces a direct solve of the v = 0 class
        to POLE_GUARD_RTOL and has no undamped pole inside the Maxwell bulk; other
        rows are replaced by ``grid_average(rows)``.  Returns (chi_p, chi_m, mean
        light occupation) and records the guard in ``ledger``.
        """
        n = atom.n_levels
        M = n * n
        U = core.hermitian_basis(n)
        Uh = U.conj().T

        def real(X):
            return np.ascontiguousarray((Uh @ X @ U).real)

        size = M if Cd_xy is None else 2 * M
        A00 = np.zeros((size, size))
        A01 = np.zeros((size, size))
        S = np.zeros((size, size))
        A00[:M, :M] = real(C_xy)
        A01[:M, :M] = real(C_z)
        S[:M, :M] = real(atom.S_v)
        if Cd_xy is not None:
            eye = np.eye(M)
            A00[:M, :M] -= gamma_out * eye
            A00[:M, M:] = gamma_in * eye
            A00[M:, :M] = gamma_out * eye
            A00[M:, M:] = real(Cd_xy) - gamma_in * eye
            A01[M:, M:] = A01[:M, :M]
        for A in (A00, A01, S):
            A[0] = 0.0
        A00[0, :n] = 1.0                 # trace row: populations are the first n coordinates
        if Cd_xy is not None:
            A00[0, M:M + n] = 1.0
        rhs = np.zeros(size)
        rhs[0] = 1.0
        readouts = np.zeros((3, size), dtype=complex)
        for row, q in ((0, +1), (1, -1)):
            for gi, ei, cg in atom.couplings[q]:
                readouts[row, :M] += cg * U[ei * n + gi]
        readouts[2, :n] = 1.0                       # light-region occupation

        # Expand about the v = +σ class, not by Schur elimination: with no transverse
        # field the paraffin Zeeman populations relax only through the optical
        # coherences, so the Doppler-free block alone is singular.
        b_z = np.ascontiguousarray(b_z, dtype=float)
        system = pole_doppler.ShiftInvertSystem(A00, A01, S, rhs, readouts)
        with core.blas_single_thread():
            poles = system.residues(b_z, center - sigma, pole_limit=POLE_LIMIT_SIGMAS * sigma)
            at_rest = system.direct(b_z, center)
        occupation = at_rest[:, 2].real
        norm = Om if Cd_xy is None else Om * occupation[:, None]
        mean = poles.gaussian_mean(center, sigma)[:, :2] / norm
        rest = poles.value_at(center)[:, :2] / norm
        direct = at_rest[:, :2] / norm
        scale = max(float(np.max(np.abs(direct))), np.finfo(float).tiny)
        relative = np.max(np.abs(rest - direct), axis=1) / scale
        relative = np.where(np.isfinite(relative), relative, np.inf)
        # An undamped pole inside the Maxwell bulk would mean a velocity class
        # without a unique steady state; distant real poles only add smooth terms.
        in_bulk = np.abs(poles.lam.real - center) <= FALLBACK_SIGMAS * sigma
        damping = np.min(np.where(in_bulk, np.abs(poles.lam.imag), np.inf), axis=1)
        failed = ((relative > POLE_GUARD_RTOL) | ~(damping > 1e-9 * sigma)
                  | ~np.all(np.isfinite(mean), axis=1))
        ledger["guard_max_relative"] = float(np.max(relative))
        ledger["min_pole_distance_mhz"] = float(np.min(damping) / (2 * np.pi * 1e6))

        chi_p, chi_m = mean[:, 0].copy(), mean[:, 1].copy()
        rows = np.flatnonzero(failed)
        if rows.size:
            chi_p[rows], chi_m[rows], classes = grid_average(rows)
            ledger["fallback_rows"] = int(rows.size)
            ledger["fallback_classes"] = classes
        return chi_p, chi_m, float(np.mean(occupation))

    @staticmethod
    def _coherences(atom, rho, Om, drive):
        # rho leading dims are arbitrary ((nv,) or (nB, nv)); reduce the last two
        # (Hilbert) axes via the polarization-resolved dipole couplings.
        shape = rho.shape[:-2]
        p = np.zeros(shape, dtype=complex)
        m = np.zeros(shape, dtype=complex)
        for gi, ei, cg in atom.couplings[+1]:
            p += cg * rho[..., ei, gi]
        for gi, ei, cg in atom.couplings[-1]:
            m += cg * rho[..., ei, gi]
        cp = p / Om
        cm = m / Om
        cprobe = (np.conj(drive.get(+1, 0.0)) * p
                  + np.conj(drive.get(-1, 0.0)) * m) / Om
        return cp, cm, cprobe

    def observables(self, raw, params, include_figures=True):
        if raw.get("budker_scan"):
            from . import magneto_budker
            return magneto_budker.observables(raw, params, include_figures=include_figures)

        if not raw["valid"]:
            fig = None
            if include_figures:
                import matplotlib.pyplot as plt

                fig, ax = plt.subplots(figsize=(8.5, 3.0))
                ax.text(0.5, 0.5,
                        f"87Rb D1 F_g={raw['Fg']}, F'_e={raw['Fe']} is not allowed "
                        f"for this compact Zeeman model.",
                        ha="center", va="center", wrap=True)
                ax.axis("off")
            return dict(metrics=[
                            dict(label="Status", value="invalid transition",
                                 kind="status", tier="hero"),
                            dict(label="Transition",
                                 value=f"F={raw['Fg']} → F'={raw['Fe']}",
                                 tier="hero"),
                        ],
                        figure=fig, tables=[])

        x = np.asarray(raw["b_ut"])
        physical_x = np.asarray(raw.get("b_physical_ut", x + raw["b_offset_ut"]))
        ls = params["line_strength"]
        xprobe = observables.chi_phys(
            raw["chi_probe"], raw["N_eff"], dipole=raw["dipole"], line_strength=ls)
        xp = observables.chi_phys(
            raw["chi_p"], raw["N_eff"], dipole=raw["dipole"], line_strength=ls)
        xm = observables.chi_phys(
            raw["chi_m"], raw["N_eff"], dipole=raw["dipole"], line_strength=ls)
        k, L = raw["k_vec"], params["cell_mm"] * 1e-3
        alpha = k * np.imag(xprobe)
        T_trans = np.exp(-alpha * L)
        rotation = 0.25 * k * L * np.real(xp - xm)
        diagnostic = _zero_field_diagnostics(physical_x, alpha, rotation)
        physical_zero_available = _scan_contains_zero(physical_x)
        bg, amp = diagnostic["background"], diagnostic["amplitude"]
        if not np.isfinite(amp):
            feature = "unavailable"
        elif amp < -0.02 * max(abs(bg), abs(diagnostic["center"]), 1e-30):
            feature = "EIT-like dip"
        elif amp > 0.02 * max(abs(bg), abs(diagnostic["center"]), 1e-30):
            feature = "EIA/MIA-like peak"
        else:
            feature = "crossover"
        width_ok = diagnostic["width_status"] == "sample-grid consistent"
        slope_ok = diagnostic["slope_status"] == "sample-grid consistent"
        central_str = f"{diagnostic['width']:.3f} µT" if width_ok else "unavailable"
        physical_T = _at_zero(physical_x, T_trans)
        command_T = _at_zero(x, T_trans)
        physical_rotation = _at_zero(physical_x, rotation)
        command_rotation = _at_zero(x, rotation)
        drive = _qwp_drive_weights(raw["qwp_deg"])
        circular_powers = [abs(drive[q]) ** 2 for q in (+1, -1)]
        linear_input = math.isclose(*circular_powers, rel_tol=0.0, abs_tol=1e-12)
        circular_input = min(circular_powers) <= 1e-12
        rotation_scope = (
            "linear-input approximation" if linear_input else
            "undefined for circular input" if circular_input else
            "unsupported for elliptical input")
        magneto_readout = dict(
            field_unit="µT", default_reference="physical zero",
            physical_zero_in_scan=bool(physical_zero_available),
            physical_zero_command_ut=-float(raw["b_offset_ut"]),
            physical_zero_transmission=physical_T,
            command_zero_transmission=command_T,
            physical_zero_rotation_rad=physical_rotation if linear_input else None,
            command_zero_rotation_rad=command_rotation if linear_input else None,
            rotation_scope=rotation_scope,
            feature=feature,
            central_absorption_per_m=diagnostic["center"],
            background_absorption_per_m=diagnostic["background"],
            feature_amplitude_per_m=diagnostic["amplitude"],
            central_width_ut=diagnostic["width"],
            left_half_height_ut=diagnostic["left"],
            right_half_height_ut=diagnostic["right"],
            width_status=diagnostic["width_status"],
            width_relative_change=diagnostic["width_relative_change"],
            slope_rad_per_ut=diagnostic["slope"] if linear_input else None,
            slope_status=diagnostic["slope_status"] if linear_input else rotation_scope,
            slope_relative_change=diagnostic["slope_relative_change"] if linear_input else None,
            sampled_grid_tolerance=SAMPLED_GRID_RTOL,
            numerical_scope="same-scan coarsening/stencil consistency; no full convergence certificate",
            model_scope=MODEL_SCOPE_STATUS,
            absolute_scale_status=raw["absolute_scale_status"],
        )
        magneto_readout = {
            key: None if isinstance(value, (float, np.floating)) and not np.isfinite(value)
            else value for key, value in magneto_readout.items()}

        # Use symmetric physical-field wings where available, matching the
        # central-feature reference even when command zero has a bias offset.
        wing_radius = min(-float(physical_x[0]), float(physical_x[-1]))
        if physical_zero_available and wing_radius > 0:
            transmission_background = float(np.mean(
                np.interp([-wing_radius, wing_radius], physical_x, T_trans)))
            background_reference = "symmetric physical-field wings"
        else:
            transmission_background = float(0.5 * (T_trans[0] + T_trans[-1]))
            background_reference = "scan-edge mean; physical zero not bracketed"
        noise = shot_noise_reference(
            T_trans, reference_transmission=transmission_background,
            frequency_hz=raw["k_vec"] * constants.C_LIGHT / (2 * np.pi),
            intensity_mw_cm2=params["intensity_mw_cm2"],
            beam_diameter_mm=params.get("beam_diameter_mm", DEFAULT_BEAM_DIAMETER_MM),
            efficiency=params.get("detection_efficiency", DEFAULT_DETECTION_EFFICIENCY),
            bandwidth_hz=params.get("noise_bandwidth_hz", DEFAULT_NOISE_BANDWIDTH_HZ),
            linear_input=linear_input)
        noise.update(enabled=bool(params.get("show_shot_noise", True)),
                     background_reference=background_reference,
                     command_field_ut=x.tolist())
        magneto_readout["shot_noise"] = noise
        local_peak_limit_t = None
        if slope_ok:
            ic = int(np.argmin(np.abs(physical_x)))
            local_peak_limit_t = max(0.0, min(-physical_x[ic-2], physical_x[ic+2])) * 1e-6
        sensitivity = magnetic_sensitivity(
            physical_T, diagnostic["slope"] * 1e6,
            frequency_hz=raw["k_vec"] * constants.C_LIGHT / (2 * np.pi),
            intensity_mw_cm2=params["intensity_mw_cm2"],
            beam_diameter_mm=params.get("beam_diameter_mm", DEFAULT_BEAM_DIAMETER_MM),
            efficiency=params.get("detection_efficiency", DEFAULT_DETECTION_EFFICIENCY),
            bandwidth_hz=params.get("noise_bandwidth_hz", DEFAULT_NOISE_BANDWIDTH_HZ),
            signal_rms_t=params.get("magnetic_signal_ft", 1.0) * 1e-15,
            signal_peak_limit_t=local_peak_limit_t,
            linear_input=linear_input, slope_valid=slope_ok)
        sensitivity.update(physical_bias_ut=0.0, response="quasistatic local rotation slope",
                           model_scope=MODEL_SCOPE_STATUS,
                           absolute_scale_status=raw["absolute_scale_status"])
        magneto_readout["magnetic_sensitivity"] = sensitivity

        def plot_shot_noise(ax, rotation_view=False):
            if not noise["enabled"]:
                return
            trace_key = "rotation_reference_trace_rad" if rotation_view else "transmission_reference_trace"
            trace = np.asarray(noise[trace_key], dtype=float)
            if not np.any(np.isfinite(trace)):
                return
            sigma = (np.asarray(noise["rotation_rms_rad"], dtype=float) * 1e3 if rotation_view
                     else np.full(x.shape, noise["transmission_reference_rms"]))
            center = 0.0 if rotation_view else transmission_background
            line, = ax.plot(x, trace * (1e3 if rotation_view else 1),
                            color="#d97706", lw=1.0, alpha=0.85, zorder=2,
                            label="Shot noise (simulated)")
            line._gabes_export_metadata = dict(
                quantity="simulated_photon_shot_noise_reference",
                y_unit="mrad" if rotation_view else "1", seed=noise["seed"],
                bandwidth_hz=noise["bandwidth_hz"],
                reference_level=center, rms=sigma.tolist())
            ax.fill_between(x, center - sigma, center + sigma,
                            color="#d97706", alpha=0.12, zorder=1,
                            label=f"+/-1 sigma ({noise['bandwidth_hz']:g} Hz ENBW)")

        def mark_field_references(ax):
            ax.axvline(0.0, color="gray", ls=":", lw=0.8,
                       label="Command zero")
            if physical_zero_available and abs(raw["b_offset_ut"]) > 1e-14:
                ax.axvline(-raw["b_offset_ut"], color="black", ls="--", lw=0.8,
                           label="Physical zero")

        def plot_legend(ax):
            # Keep reference artists and their export labels, but avoid
            # repeating the zero guide and noise envelope in the legend.
            handles, labels = ax.get_legend_handles_labels()
            entries = [(h, label) for h, label in zip(handles, labels)
                       if label != "Command zero" and not label.startswith("+/-1 sigma")]
            ax.legend([h for h, _ in entries], [label for _, label in entries], fontsize=8)

        title = (f"87Rb D1 {raw['cell_type']}  F={raw['Fg']} -> F'={raw['Fe']},  "
                 f"QWP={raw['qwp_deg']:.1f} deg, I={params['intensity_mw_cm2']:.2f} mW/cm^2")
        fig = None
        if params.get("signal_type", SIGNAL_TRANSMISSION) == SIGNAL_NMOR:
            metrics = [
                dict(label="Rotation at physical B=0",
                     value=(f"{physical_rotation*1e3:.2f} mrad"
                            if physical_zero_available else "unavailable")),
                dict(label="Slope dtheta/dB",
                     value=(f"{diagnostic['slope']*1e3:.2f} mrad/µT"
                            if slope_ok else "unavailable"),
                     **({"tier": "hero"} if slope_ok else {})),
                dict(label="Peak |rotation|", value=f"{np.max(np.abs(rotation))*1e3:.2f} mrad",
                     **({"tier": "hero"} if slope_ok else {})),
                dict(label="Slope status", value=diagnostic["slope_status"], kind="status",
                     **({} if slope_ok else {"tier": "hero"}),
                     help="A fourth-degree local derivative is compared with a "
                          "twice-spaced stencil. A change above 5% suppresses "
                          "the numeric hero; agreement is sampled consistency, "
                          "not a full convergence certificate."),
                dict(label="Absolute rotation/slope status",
                     value=raw["absolute_scale_status"], kind="status",
                     help="The internal state and density convention is fixed, but "
                          "absolute NMOR rotation and slope still require an "
                          "independent measurement."),
            ]
            note = ("NMOR slope is evaluated at physical zero. Command-zero rotation "
                    "is secondary. Local sampled-grid consistency uses a 5% check "
                    "without additional atomic solves; absolute rotation and slope "
                    "remain externally unvalidated.")
            if not linear_input:
                # Re(chi_+ - chi_-) is a plane-angle approximation for linear
                # input. Elliptical input needs a Jones/Stokes readout, and a
                # circular input has no polarization-plane angle at all.
                metrics = [metric for metric in metrics
                           if metric["label"] == "Absolute rotation/slope status"]
                note = ("The plotted circular-birefringence diagnostic is not a "
                        "polarization-plane rotation for this input. Elliptical "
                        "input requires Jones/Stokes propagation; circular input "
                        "has no defined polarization-plane angle.")
            metrics.append(dict(label="Rotation readout status", value=rotation_scope,
                                kind="status", **({} if linear_input else {"tier": "hero"})))
            if include_figures:
                import matplotlib.pyplot as plt

                fig, ax = plt.subplots(figsize=(8.5, 4.6))
                ax.plot(x, rotation * 1e3, color="#7c3aed", lw=1.8,
                        label="NMOR rotation" if linear_input else "Birefringence diagnostic")
                ax.axhline(0, color="black", lw=0.6)
                plot_shot_noise(ax, rotation_view=True)
                mark_field_references(ax)
                plot_legend(ax)
                ax.set_ylabel("Polarization rotation  [mrad]" if linear_input else
                              "Circular-birefringence diagnostic  [mrad]")
                ax.set_xlabel("Commanded magnetic field B  [uT]")
                ax.set_title(title)
                fig.tight_layout()
        else:
            metrics = [
                dict(label="Transmission at physical B=0",
                     value=f"{physical_T:.3f}" if physical_zero_available else "unavailable",
                     **({"tier": "hero"} if not width_ok else {})),
                dict(label="Zero-field feature", value=feature,
                     **({"tier": "hero"} if width_ok else {})),
                dict(label="Central width", value=central_str,
                     **({"tier": "hero"} if width_ok else {}),
                     help="Distance between both interpolated half-height crossings "
                          "about physical zero; the background uses symmetric wings."),
                dict(label="Width status", value=diagnostic["width_status"], kind="status",
                     **({} if width_ok else {"tier": "hero"}),
                     help="The width and feature sign are checked against every "
                          "second existing B sample. A change above 5% suppresses "
                          "the numeric hero; agreement is sampled consistency."),
                dict(label="Absolute contrast status",
                     value=raw["absolute_scale_status"], kind="status",
                     help="The internal state and density convention is fixed, but "
                          "absolute transmission contrast still requires an "
                          "independent measurement."),
            ]
            note = ("Transmission and feature sign use physical zero and symmetric "
                    "physical-field wings. Command-zero transmission is secondary. "
                    "The width check measures sampled-grid consistency without "
                    "additional atomic solves; absolute contrast remains externally "
                    "unvalidated.")
            if include_figures:
                import matplotlib.pyplot as plt

                fig, axT = plt.subplots(figsize=(8.5, 4.6))
                axT.plot(x, T_trans, color="#1f77b4", lw=1.8, label="Transmission")
                plot_shot_noise(axT)
                mark_field_references(axT)
                plot_legend(axT)
                axT.set_ylabel("Transmission")
                axT.set_xlabel("Commanded magnetic field B  [uT]")
                axT.set_title(title)
                fig.tight_layout()

        rotation_view = params.get("signal_type", SIGNAL_TRANSMISSION) == SIGNAL_NMOR
        signal_amplitude = (np.abs(rotation) if rotation_view else
                            np.abs(T_trans - transmission_background))
        snr_rms = (np.asarray(noise["rotation_rms_rad"], dtype=float) if rotation_view else
                   np.full(x.shape, noise["transmission_reference_rms"]
                           if noise["transmission_status"] == "available" else np.nan))
        snr_values = np.full(x.shape, np.nan)
        np.divide(signal_amplitude, snr_rms, out=snr_values,
                  where=np.isfinite(snr_rms) & (snr_rms > 0))
        peak_index = int(np.nanargmax(snr_values)) if np.any(np.isfinite(snr_values)) else None
        peak_snr = float(snr_values[peak_index]) if peak_index is not None else None
        visibility = ("Unavailable" if peak_snr is None else
                      "Above shot noise" if peak_snr > 1 else
                      "At shot noise" if peak_snr == 1 else "Below shot noise")
        snr_formula = ("max_B |theta(B)| / sigma_theta(B)" if rotation_view else
                       "max_B |T(B) - T_background| / sigma_T_background")
        magneto_readout["shot_noise_snr"] = dict(
            signal_type=params.get("signal_type", SIGNAL_TRANSMISSION),
            convention="sampled signal amplitude / photon-shot-noise RMS; not a power ratio or dB",
            formula=snr_formula, scope="largest ratio in the displayed scan; photon shot noise only",
            reference_level=0.0 if rotation_view else transmission_background,
            signal_unit="rad" if rotation_view else "1",
            amplitude_ratio=[float(v) if np.isfinite(v) else None for v in snr_values],
            peak_amplitude_ratio=peak_snr,
            peak_command_field_ut=float(x[peak_index]) if peak_index is not None else None,
            peak_physical_field_ut=float(physical_x[peak_index]) if peak_index is not None else None,
            peak_signal_amplitude=float(signal_amplitude[peak_index]) if peak_index is not None else None,
            peak_noise_rms=float(snr_rms[peak_index]) if peak_index is not None else None,
            valid_samples=int(np.count_nonzero(np.isfinite(snr_values))),
            bandwidth_hz=noise["bandwidth_hz"], visibility=visibility,
            noise_status=noise["rotation_status"] if rotation_view else noise["transmission_status"],
        )
        snr_help = (f"Largest sampled signal amplitude / shot-noise RMS in this scan: {snr_formula}. "
                    f"ENBW = {noise['bandwidth_hz']:g} Hz. "
                    "1× means signal amplitude equals one noise RMS; it is not a hard detection threshold. "
                    "The ratio uses theoretical signals and RMS, independent of the simulated zigzag. "
                    "Electronics, technical and atomic spin noise are omitted.")
        if peak_index is not None:
            snr_help += (f" Maximum at commanded B = {x[peak_index]:.4g} µT "
                         f"(physical B = {physical_x[peak_index]:.4g} µT).")
        else:
            snr_help += f" Unavailable: {magneto_readout['shot_noise_snr']['noise_status']}."
        # Replace the secondary command-zero value in the strip, preserving
        # that reference in the structured readout and derived table.
        metrics.insert(1, dict(label="Peak shot-noise SNR",
                               value=f"{peak_snr:.4g} ×" if peak_snr is not None else "unavailable",
                               delta=visibility, help=snr_help))
        if rotation_view:
            field_asd = sensitivity["field_asd_t_per_sqrt_hz"]
            magnetic_snr = sensitivity["signal_to_noise_ratio"]
            metrics.extend([
                dict(label="Signal linearity", value=sensitivity["signal_linearity_status"], kind="status"),
                dict(label="PSN magnetic sensitivity",
                     value=f"{field_asd * 1e15:.4g} fT/√Hz" if field_asd is not None else "unavailable",
                     help="Rotation noise ASD divided by the local dtheta/dB at physical B=0. "
                          "Conditional photon-shot-noise estimate; no atomic bandwidth or total-noise claim."),
                dict(label="Small-field SNR",
                     value=f"{magnetic_snr:.4g} ×" if magnetic_snr is not None else "unavailable",
                     help="Declared RMS magnetic signal / (field-noise ASD × √ENBW), "
                          "using the quasistatic slope at physical B=0."),
            ])
        rms_values = np.asarray(noise["rotation_rms_rad"], dtype=float) * 1e3
        finite_rms = rms_values[np.isfinite(rms_values)]
        if rotation_view:
            rms_text = (f"{finite_rms.min():.2e}–{finite_rms.max():.2e} mrad"
                        if finite_rms.size else "unavailable")
            noise_status = noise["rotation_status"]
        else:
            rms_text = (f"{noise['transmission_reference_rms']:.2e} (T)"
                        if noise["transmission_status"] == "available" else "unavailable")
            noise_status = noise["transmission_status"]
        metrics.append(dict(label="Shot noise RMS", value=rms_text,
                            help="Photon shot noise integrated over the selected ENBW. "
                                 "Transmission uses the background; rotation lists the "
                                 "range across B. Detection settings are under Advanced."))
        note += (f"\n\nShot-noise reference: {noise['bandwidth_hz']:g} Hz ENBW, "
                 f"effective uniform beam diameter {noise['beam_diameter_mm']:g} mm, "
                 f"incident power {noise['incident_power_w'] * 1e6:.3g} µW, "
                 f"detection efficiency {noise['efficiency']:.0%}. "
                 "The orange zigzag is a fixed-seed simulated realization; the shaded "
                 "band is ±1σ, not a hard detection threshold. Transmission noise is "
                 f"centered on the {background_reference}; compare signal departures "
                 "from that background. Rotation noise is centered on zero and uses "
                 "the transmitted photon flux at each B. No visual amplification is "
                 "applied. The reference assumes a noiseless incident-power normalization "
                 "and omits technical, electronics and atomic spin noise. "
                 f"Readout status: {noise_status}. Peak shot-noise SNR uses {snr_formula}; "
                 "it describes the largest sampled ratio in this scan, not the value at B=0 "
                 "or a magnetic-field sensitivity.")
        metrics.append(dict(label="Model scope", value=MODEL_SCOPE_STATUS, kind="status",
                            help="Only one ground/excited hyperfine pair is retained. "
                                 "Other ground-F leakage, repumping and interference "
                                 "with the other excited hyperfine manifold are omitted."))
        # The shared strip reserves numeric cells for values; otherwise a
        # status-only hero would fall back to a large unavailable SNR/angle.
        for metric in metrics:
            if metric["value"] == "unavailable":
                metric["kind"] = "status"
        larmor_hz_per_g = constants.MU_B * abs(raw["gFg"]) / (2 * np.pi * constants.HBAR) * 1e-4
        buffer_mhz = raw.get("buffer_gamma", 0.0) / (2 * np.pi) / 1e6
        elastic_mhz = raw.get("gamma_elastic", 0.0) / (2 * np.pi) / 1e6
        derived = derived_table([
            ("Isotope / line", f"{raw['isotope']} {raw['line']}"),
            ("Cell model", f"{raw['cell_type']}"),
            ("Transition", f"F={raw['Fg']} → F'={raw['Fe']}"),
            ("Relative strength S_FF'", f"{raw['strength']:.4f}"),
            ("g_F ground / excited", f"{raw['gFg']:.4f} / {raw['gFe']:.4f}"),
            ("Ground Larmor scale", f"{larmor_hz_per_g/1e6:.3f} MHz/G"),
            ("Ω / Γ", f"{raw['omega_rabi']/raw['gamma_natural']:.3f}"),
            ("Natural Γ / 2π", f"{raw['gamma_natural']/(2*np.pi)/1e6:.4f} MHz"),
            ("Homogeneous optical FWHM / 2π",
             f"{raw['gamma']/(2*np.pi)/1e6:.4f} MHz"),
            ("Ne broadening / 2π", f"{buffer_mhz:.3f} MHz"),
            ("Elastic optical dephasing / 2π", f"{elastic_mhz:.3f} MHz"),
            ("Buffer reservoir reload / 2π",
             f"{raw['gamma_transit']/(2*np.pi)/1e3:.3f} kHz"),
            ("Ground pure dephasing / 2π",
             f"{raw['gamma_depol']/(2*np.pi)/1e3:.3f} kHz"),
            ("I_sat (D1 scale)", f"{raw['isat']:.4f} mW/cm²"),
            ("Longitudinal B0 offset", f"{raw['b_offset_ut']:.4f} µT"),
            ("Physical zero on command axis", f"{-raw['b_offset_ut']:.4f} µT"),
            ("Transmission at command B=0", f"{command_T:.6g}"),
            ("Rotation at command B=0",
             f"{command_rotation * 1e3:.6g} mrad" if linear_input else rotation_scope),
            ("Width sampled-grid change",
             f"{diagnostic['width_relative_change']:.2%}" if np.isfinite(
                 diagnostic['width_relative_change']) else "unavailable"),
            ("Slope sampled-stencil change",
             f"{diagnostic['slope_relative_change']:.2%}" if np.isfinite(
                 diagnostic['slope_relative_change']) else "unavailable"),
            ("Residual transverse B", f"{raw['b_perp_ut']:.4f} µT"),
            ("N(87Rb)", f"{raw['N']:.3e} /m³"),
            ("Optical-state normalization", raw["light_state_normalization"]),
            ("Conditional light-state trace", f"{raw['light_state_trace']:.12f}"),
            ("Light-region occupation",
             f"{raw['light_region_occupation']:.9f} (diagnostic only)"),
            ("Density convention", raw["density_convention"]),
            ("Doppler average", raw.get("doppler_method", DOPPLER_OFF)),
            ("Maxwell Doppler width σ_D / 2π",
             f"{raw.get('doppler_sigma_mhz', float('nan')):.2f} MHz"),
            ("Pole guard (max relative)",
             "n/a" if raw.get("doppler_guard_max_relative") is None
             else f"{raw['doppler_guard_max_relative']:.1e}"),
            ("Grid-fallback B rows", f"{raw.get('doppler_fallback_rows', 0)}"),
            ("Two-region γ_out / 2π", f"{raw['gamma_out']/(2*np.pi)/1e3:.3f} kHz"),
            ("Two-region γ_in / 2π", f"{raw['gamma_in']/(2*np.pi)/1e3:.3f} kHz"),
            ("Wall γ / 2π", f"{raw['gamma_wall']/(2*np.pi):.3f} Hz"),
        ])
        return dict(metrics=metrics, figure=fig, magneto_readout=magneto_readout,
                    hero_count=sum(m.get("tier") == "hero" for m in metrics), tables=[
            {"title": "Notes", "markdown": note},
            derived,
        ])


# One public model for every registered magneto readout, axis and alias.
from .magneto_unified import UnifiedMagnetoScheme as MagnetoScheme
