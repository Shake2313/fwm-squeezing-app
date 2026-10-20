"""
Cluster C - 87Rb D1 magneto-optics.

One Hanle-configuration transmission engine now covers EIT-like dips and
EIA/MIA-like peaks.  The zero-field feature sign is a result, not a readout
selector: probe ellipticity, cell relaxation, residual transverse field, and
the addressed Zeeman manifold decide the line shape.
"""
import math
import re

import numpy as np

from .. import constants, core, kernels, observables, pole_doppler, species, zeeman
from ..lineshape import halfwidth_from_center, lorentz_fwhm
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

GJ_5S12 = 2.002_331_070
GJ_5P12 = 0.666
GI_RB87 = -0.000_995_141_4

CELL_BUFFER = "Buffer gas cell"
CELL_PARAFFIN = "Paraffin coated cell"
SIGNAL_TRANSMISSION = "Transmission"
SIGNAL_NMOR = "NMOR rotation"

LIGHT_STATE_NORMALIZATION = "trace-one conditional light-region state"
DENSITY_CONVENTION = "local vapor density N, applied once"
ABSOLUTE_SCALE_STATUS = "external validation required"

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


def _hyperfine_gf(I, J, F, gJ):
    """Low-field hyperfine Lande g_F in Bohr-magneton units."""
    F = float(F)
    if F == 0:
        return 0.0
    ff = F * (F + 1)
    electronic = gJ * (ff + J * (J + 1) - I * (I + 1)) / (2 * ff)
    nuclear = GI_RB87 * (ff + I * (I + 1) - J * (J + 1)) / (2 * ff)
    return electronic + nuclear


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


class MagnetoScheme(Scheme):
    cluster = "C - Magneto-optics"
    presets_group = "Default"
    cache_version = "polarized-two-region-v6-exact-doppler"
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
            self.caption = ("87Rb D1 near-zero-field readouts of the ground-state "
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
                      help=f"Converted with Steck D1 I_sat = {D1_ISAT_MW_CM2:.4f} mW/cm²."),
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
                      help="Sets 87Rb vapor density and Doppler width."),
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
            ParamSpec("collisional_depol_khz", "Collisional depolarization", "Atomic", 0.5,
                      0.0, 200.0, 0.5, "kHz", advanced=True,
                      visible_if={"cell_type": CELL_BUFFER},
                      advanced_group="Relaxation overrides",
                      help="Spin-destruction ground-coherence loss without population "
                           "reload; a buffer-gas consequence exposed as an override."),
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
            ParamSpec("doppler", "Doppler averaging", "Numerics", "on",
                      choices=("on", "off"), advanced=True,
                      help="On: exact Maxwell-Boltzmann average of the light-region "
                           "optical Doppler shift (pole–residue form, no velocity "
                           "grid). Off: the v = 0 class only.",
                      advanced_group="Numerics"),
            ParamSpec("scan_points", "B scan points", "Numerics", 121,
                      51, 401, 10, "", advanced=True, advanced_group="Numerics"),
        ]
        return specs

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
            return []
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
            "**87Rb D1 Hanle / magneto-optical-rotation model.** These are two "
            "distinct effects of the ground-state Zeeman manifold near zero field: "
            "the Hanle effect is a transmission resonance from ground-state "
            "coherence, while magneto-optical rotation (MOR/NMOR) is a rotation of "
            "the probe polarization plane. The transmission readout reports "
            "whether the zero-field Hanle feature is EIT-like transparency, EIA/MIA-like "
            "absorption, or a crossover.  Paraffin-coated cells use a two-region "
            "light/dark OBE exchange model so wall-preserved ground coherence can "
            "create a narrow Ramsey feature on a broad transit pedestal.\n\n"
            "**References**\n"
            "- D. A. Steck, *Rubidium 87 D Line Data*, http://steck.us/alkalidata.\n"
            "- H. J. Lee and H. S. Moon, magnetic-field-induced absorption in a "
            "paraffin-coated rubidium vapor cell, *JOSA B* **30**, 2301 (2013).\n"
            "- H. S. Moon and H. J. Kim, Ramsey EIA to MIT transformation in a "
            "paraffin-coated Rb vapor cell, *Optics Express* **22**, 18604 (2014).\n"
            "- Y.-J. Yu, H. J. Lee, and H. S. Moon, level-crossing absorption "
            "with narrow spectral width in Rb vapor with buffer gas, "
            "*Phys. Rev. A* **81**, 023416 (2010)."
        )

    def compute(self, params):
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

        gFg = _hyperfine_gf(RB87.I, RB87.Jg, Fg, GJ_5S12) if Fg in Fg_allowed else 0.0
        gFe = _hyperfine_gf(RB87.I, JE_D1, Fe, GJ_5P12) if Fe in Fe_allowed else 0.0

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
            return MagnetoScheme._conditional_light_state(rho_l)
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
        return MagnetoScheme._conditional_light_state(rho_l)

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

        x = raw["b_ut"]
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
        ic = int(np.argmin(np.abs(x)))
        bg = 0.5 * (alpha[0] + alpha[-1])
        amp = alpha[ic] - bg
        if amp < -0.02 * max(abs(bg), abs(alpha[ic]), 1e-30):
            feature = "EIT-like dip"
        elif amp > 0.02 * max(abs(bg), abs(alpha[ic]), 1e-30):
            feature = "EIA/MIA-like peak"
        else:
            feature = "crossover"
        fwhm = lorentz_fwhm(x, alpha)
        central_hw = halfwidth_from_center(x, alpha)
        fwhm_str = f"{fwhm:.3f} µT" if fwhm == fwhm else "n/a"
        central_str = f"{2*central_hw:.3f} µT" if central_hw == central_hw else "n/a"

        title = (f"87Rb D1 {raw['cell_type']}  F={raw['Fg']} -> F'={raw['Fe']},  "
                 f"QWP={raw['qwp_deg']:.1f} deg, I={params['intensity_mw_cm2']:.2f} mW/cm^2")
        fig = None
        if params.get("signal_type", SIGNAL_TRANSMISSION) == SIGNAL_NMOR:
            slope = np.gradient(rotation, x)[ic]
            drive = _qwp_drive_weights(raw["qwp_deg"])
            circular_powers = [abs(drive[q]) ** 2 for q in (+1, -1)]
            linear_input = math.isclose(*circular_powers, rel_tol=0.0, abs_tol=1e-12)
            circular_input = min(circular_powers) <= 1e-12
            rotation_scope = (
                "linear-input approximation" if linear_input else
                "undefined for circular input" if circular_input else
                "unsupported for elliptical input")
            metrics = [
                dict(label="Rotation at B=0", value=f"{rotation[ic]*1e3:.2f} mrad"),
                dict(label="Slope dtheta/dB", value=f"{slope*1e3:.2f} mrad/µT",
                     tier="hero"),
                dict(label="Peak |rotation|", value=f"{np.max(np.abs(rotation))*1e3:.2f} mrad",
                     tier="hero"),
                dict(label="Absolute rotation/slope status",
                     value=raw["absolute_scale_status"], kind="status",
                     help="The internal state and density convention is fixed, but "
                          "absolute NMOR rotation and slope still require an "
                          "independent measurement."),
            ]
            note = ("NMOR readout: zero crossing near B=0; slope is the magnetometer "
                    "signal. Absolute rotation and slope remain externally unvalidated.")
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
                ax.plot(x, rotation * 1e3, color="#7c3aed", lw=1.8)
                ax.axhline(0, color="black", lw=0.6)
                ax.axvline(0, color="gray", ls=":", lw=0.8)
                ax.set_ylabel("Polarization rotation  [mrad]" if linear_input else
                              "Circular-birefringence diagnostic  [mrad]")
                ax.set_xlabel("Magnetic field B  [µT]")
                ax.set_title(title)
                fig.tight_layout()
        else:
            metrics = [
                dict(label="Transmission at B=0", value=f"{T_trans[ic]:.3f}",
                     **({"tier": "hero"} if not np.isfinite(central_hw) else {})),
                dict(label="Zero-field feature", value=feature,
                     **({"tier": "hero"} if np.isfinite(central_hw) else {})),
                dict(label="Feature FWHM", value=fwhm_str),
                dict(label="Central width", value=central_str,
                     **({"tier": "hero"} if np.isfinite(central_hw) else {})),
                dict(label="Absolute contrast status",
                     value=raw["absolute_scale_status"], kind="status",
                     help="The internal state and density convention is fixed, but "
                          "absolute transmission contrast still requires an "
                          "independent measurement."),
            ]
            if not np.isfinite(central_hw):
                metrics.append(dict(
                    label="Width status", value="unresolved in scan", kind="status",
                    help="No finite half-height crossing was found around the "
                         "commanded B=0 sample; widen or recenter the B scan.",
                    tier="hero"))
            note = ("Transmission readout: feature sign is classified from absorption "
                    "at B=0 relative to the scan wings. Absolute contrast remains "
                    "externally unvalidated.")
            if include_figures:
                import matplotlib.pyplot as plt

                fig, axT = plt.subplots(figsize=(8.5, 4.6))
                axT.plot(x, T_trans, color="#1f77b4", lw=1.8)
                axT.axvline(0, color="gray", ls=":", lw=0.8)
                axT.set_ylabel("Transmission")
                axT.set_xlabel("Magnetic field B  [µT]")
                axT.set_title(title)
                fig.tight_layout()

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
            ("Collisional ground dephasing / 2π",
             f"{raw['gamma_depol']/(2*np.pi)/1e3:.3f} kHz"),
            ("I_sat (D1 scale)", f"{raw['isat']:.4f} mW/cm²"),
            ("Longitudinal B0 offset", f"{raw['b_offset_ut']:.4f} µT"),
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
        return dict(metrics=metrics, figure=fig, tables=[
            {"title": "Notes", "markdown": note},
            derived,
        ])
