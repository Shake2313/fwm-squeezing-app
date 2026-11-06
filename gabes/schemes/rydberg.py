"""
Cluster A -- Rydberg-EIT electrometry.

The optical core is a compact 85Rb cascade model:
    5S1/2 F=3 -> 5P3/2 F'=4 -> 40D5/2, with a 37 GHz microwave leg
    40D5/2 -> 39F7/2.

Around the LO-dressed steady state, a first-order finite-IF Liouvillian solve
now provides the weak-SIG superheterodyne response.  It is combined with an
explicit RF dipole convention and a balanced-photodetector noise budget to
report conditional absolute field sensitivity.  This is frequency-domain
linear response, not a full time-domain lock-in simulation.  Full Zeeman /
polarization manifolds and spatial microwave standing waves remain out of
scope; the scalar RF angular factor exposes that missing calibration instead
of hiding it.
"""
import functools
import threading
from collections import OrderedDict

import numpy as np

from .. import affine_doppler, atoms, beam, constants, core, observables, species
from .. import rydberg_electrometry as electrometry
from .. import rydberg_experiment as rydberg_experiment
from .. import rydberg_opacity as opacity
from ..lineshape import level_crossing
from ..report import derived_table
from .base import ExtraView, ParamSpec, Scheme

MHZ = 2 * np.pi * 1e6
PROBE_RABI = 1e-3
BOHR_RADIUS_M = 5.291_772_109_03e-11
# ARC 3.10.2 result for the stretched-state sigma+ matrix element
# 85Rb 40D5/2,mJ=5/2 -> 39F7/2,mJ=7/2.  The actual warm-vapor ensemble can have
# a smaller effective coupling; `rf_angular_factor` is therefore kept explicit.
RF_DIPOLE_STRETCHED_EA0 = 1326.257_243_451_434_3
# 5P3/2 -> 40D5/2 "blue" Rydberg coupling-laser wavelength (Ju et al. Fig. 1a).
# Sets the residual two-photon Doppler ratio (k_probe - k_coupling)/k_probe for
# the counter-propagating geometry the static model assumes.
COUPLING_WAVELENGTH_NM = 481.0

# --- Explicit compact-model optical/transit anchors ---
# Transit-time broadening of the 5S-40D coherence: an atom crosses the beam in
# ~d/v_mp. This phenomenological factor is not an independently validated
# conversion between the reported beam diameter and the measured EIT width.
TRANSIT_FACTOR = 0.6
# Ω_P/2π [MHz] at 6 µW / 0.15 mm. Kept as a live model/calibration input:
# Ju Fig. 2(b)'s ~1.6 MHz point occurs at low probe power, not at 6 µW.
# No claim that this default reproduces the full measured power dependence.
PROBE_RABI_REF_MHZ = 2.0
SUPERHET_ATOMIC_CACHE_SIZE = 16
SUPERHET_CLASS_CHUNK = 16384     # velocity classes per finite-IF batch solve
MIRROR_TIE_RTOL = 1e-5           # mirror slope maxima of a symmetric line are one tie


def _finite_scenario_envelope(curves):
    """Finite scenario envelope and sample counts, retaining all-invalid gaps."""
    values = np.asarray(curves, dtype=float)
    finite = np.isfinite(values)
    counts = np.sum(finite, axis=0)
    lower = np.min(np.where(finite, values, np.inf), axis=0)
    upper = np.max(np.where(finite, values, -np.inf), axis=0)
    return (np.where(counts > 0, lower, np.nan),
            np.where(counts > 0, upper, np.nan), counts)


@functools.lru_cache(maxsize=1)
def _probe_line():
    """Invariant 85Rb D2 probe-line constants (species-only), computed once."""
    rb85 = species.RB85
    Je, nu0, gamma_mhz, _, _ = rb85.line("D2")
    gamma_e = gamma_mhz * MHZ
    lam = constants.C_LIGHT / nu0
    k_vec = 2 * np.pi / lam
    dipole = np.sqrt(species.reduced_dipole_sq(gamma_e, lam, rb85.Jg, Je))
    return dict(nu0=nu0, gamma_mhz=gamma_mhz, gamma_e=gamma_e,
                lam=lam, k_vec=k_vec, dipole=dipole)


@functools.lru_cache(maxsize=1)
def _cascade_skeleton():
    """Constant 4-level cascade topology (decay + per-level Doppler ratios).

    Only the dephasing rates vary per call, so the topology is assembled once and
    `_atom` injects dephasing. The 40D / 39F two-photon levels carry the residual
    Doppler ratio (k_probe - k_coupling)/k_probe; the 5P intermediate level keeps
    the full probe-k ratio of 1.
    """
    line = _probe_line()
    gamma_r = 0.02 * MHZ
    kp = line["k_vec"]
    kc = beam.wavevector_from_wavelength_nm(COUPLING_WAVELENGTH_NM)
    ratio = beam.collinear_residual_k_ratio(kp, kc)
    return dict(
        labels=("5S F=3", "5P F'=4", "40D", "39F"),
        decay=((1, 0, line["gamma_e"]), (2, 1, gamma_r), (3, 2, gamma_r)),
        doppler_levels=(1, 2, 3),
        doppler_ratios=((2, ratio), (3, ratio)),
    )


class RydbergEITScheme(Scheme):
    name = "rydberg_eit"
    cluster = "A — Absorption"
    title = "Rydberg-EIT electrometry"
    caption = ("85Rb cascade EIT / microwave Autler-Townes electrometry. "
               "The static optical spectrum follows the 5S-5P-40D ladder; "
               "the 37 GHz RF leg dresses 40D-39F.")
    cache_version = "rydberg-eit-v8-cp-optical-detunings"
    defaults_version = "rydberg-eit-v7-optical-detunings"
    supports_headless_observables = True
    model_comparison_parameter = "lo_rabi_mhz"
    workspace_inputs = ("probe_power_uw", "coupling_power_mw", "lo_rabi_mhz")
    model_comparison_controls = ("eit_compare", "at_compare", "temperature_compare")
    _superhet_atomic_cache_lock = threading.RLock()

    REFERENCE_SENSITIVITY_NV_CM_SQRT_HZ = 12.5
    REFERENCE_PSN_LIMIT_NV_CM_SQRT_HZ = 11.2
    REFERENCE_UNCERTAINTY_NV_CM_SQRT_HZ = 0.8

    _REF = dict(
        probe_power_uw=6.0,
        coupling_power_mw=30.0,
        beam_diameter_mm=0.15,
        cell_mm=50.0,
        temp_c=20.0,
        temperature_model="Linked",
        heater_setpoint_c=20.0,
        effective_temp_c=20.0,
        cold_spot_temp_c=20.0,
        coupling_rabi_mhz=3.0,
        probe_rabi_ref_mhz=PROBE_RABI_REF_MHZ,
        scan_axis="Probe detuning",
        probe_detuning_mhz=0.0,
        coupling_detuning_mhz=0.0,
        lo_rabi_mhz=3.7,
        mw_detuning_mhz=0.0,
        mw_frequency_ghz=37.0,
        if_khz=40.0,
        rydberg_dephasing_mhz=0.10,
        temp_dephasing_mhz_per_c=0.0,
        density_dephasing_mhz_per_1e16_m3=0.0,
        rf_dephasing_mhz=1.00,
        residual_zeeman_mhz=1.50,
        doppler="off",
        rf_transition_dipole_ea0=RF_DIPOLE_STRETCHED_EA0,
        rf_angular_factor=1.0,
        rf_field_convention="RMS",
        detector_quantum_efficiency=0.95,
        detector_path_efficiency=0.50,
        detector_reference_power_ratio=1.0,
        detector_electronic_noise_pa_sqrt_hz=0.0,
        detector_rin_per_sqrt_hz=0.0,
        detector_rin_correlation=1.0,
        measurement_enbw_hz=1.0,
        effective_opacity_participation=0.001,
        effective_opacity_relative_uncertainty=0.5,
        atom_participation_fraction=1.0,
        beam_overlap_efficiency=1.0,
        sam_enabled="Off",
        sam_source_power_dbm=-40.0,
        sam_antenna_gain_dbi=10.0,
        sam_distance_m=0.30,
        sam_cable_loss_db=0.0,
        sam_additional_loss_db=0.0,
        sam_field_correction=1.0,
        sam_source_power_std_db=0.0,
        sam_antenna_gain_std_db=0.0,
        sam_cable_loss_std_db=0.0,
        sam_additional_loss_std_db=0.0,
        sam_distance_std_m=0.0,
        sam_field_correction_std=0.0,
        sam_antenna_max_dimension_m=0.10,
    )

    def atom_conditions(self, params):
        fields = ("cell_mm", "temperature_model", "temp_c", "heater_setpoint_c",
                  "effective_temp_c", "cold_spot_temp_c")
        return dict(fields=tuple(name for name in fields if name not in self.comparison_input_names(params)),
                    fixed={"Atom": "Rb", "Isotope": "85", "Transition line": "D2"})

    def comparison_input_names(self, params):
        from .rydberg_comparison import comparison_key
        return tuple(dict.fromkeys((*self.workspace_inputs, comparison_key(params))))

    def comparison_control_names(self, params):
        selector = "eit_compare" if params.get("view", "EIT") == "EIT" else "at_compare"
        return (selector,) + (("temperature_compare",) if params.get(selector) == "Temperature"
                             and params.get("temperature_model") == "Separated" else ())

    def param_schema(self):
        r = self._REF
        return [
            ParamSpec("eit_compare", "Compare", "Comparison", "Probe power",
                      choices=("Probe power", "Coupling power", "Temperature", "Residual width"),
                      recompute=False, visible_if={"view": "EIT"}),
            ParamSpec("at_compare", "Compare", "Comparison", "RF Rabi",
                      choices=("RF Rabi", "RF detuning", "Probe power", "Coupling power", "Temperature", "IF frequency"),
                      recompute=False, visible_if={"view": "AT electrometry"}),
            ParamSpec("temperature_compare", "Temperature scan", "Comparison", "Effective vapor temperature",
                      choices=("Effective vapor temperature", "Cold-spot temperature"), recompute=False,
                      visible_if={"temperature_model": "Separated"},
                      help="In Separated mode these are distinct physical temperatures. Heater setpoint stays fixed."),
            ParamSpec("view", "Regime", "Regime", "EIT",
                      choices=("EIT", "AT electrometry"), control="segmented",
                      applies_defaults=True,
                      help="EIT: no microwave dressing (Ju et al. Fig. 2a). "
                           "AT: the Rydberg RF leg is dressed."),
            ParamSpec("probe_power_uw", "Probe power", "Cell & beams",
                      r["probe_power_uw"], 0.1, 20.0, 0.1, "µW",
                      help="780 nm probe power. Drives a weak-probe Ω_P via √(P)/d "
                           "(anchored), so raising it power-broadens the EIT line "
                           "(Ju et al. Fig. 2b). Reference 6 µW."),
            ParamSpec("coupling_power_mw", "Coupling power", "Cell & beams",
                      r["coupling_power_mw"], 1.0, 80.0, 0.5, "mW",
                      help="481 nm 5P→40D coupling-beam power. Drives Ω_c via the "
                           "√(P/P_ref) intensity scaling at fixed waist, anchored so "
                           "the reference power reproduces Coupling Rabi."),
            ParamSpec("beam_diameter_mm", "Beam diameter", "Cell & beams",
                      r["beam_diameter_mm"], 0.05, 1.0, 0.01, "mm",
                      help="Model beam diameter. Sets Ω_c, Ω_P ∝ 1/d (intensity ∝ "
                           "P/d²) AND the transit-time broadening ∝ v_thermal/d — "
                           "a tighter beam is brighter but transit-broadens the EIT. "
                           "The Gold source does not specify a 1/e²/FWHM convention."),
            ParamSpec("cell_mm", "Cell length", "Cell & beams", r["cell_mm"],
                      1.0, 100.0, 0.5, "mm", recompute=False),
            ParamSpec(
                "temperature_model", "Temperature definition", "Cell & beams",
                r["temperature_model"], choices=("Linked", "Separated"),
                control="segmented",
                help="Linked uses one temperature for the heater metadata, atomic "
                     "motion, and vapor pressure. Separated exposes heater setpoint, "
                     "effective vapor temperature, and cold-spot temperature as "
                     "different measured/modelled quantities."),
            ParamSpec(
                "temp_c", "Linked cell temperature", "Cell & beams", r["temp_c"],
                15.0, 80.0, 1.0, "°C",
                visible_if={"temperature_model": "Linked"},
                help="Used for both the atomic velocity/dephasing calculation and "
                     "the saturated vapor pressure in Linked mode."),
            ParamSpec(
                "heater_setpoint_c", "Heater setpoint", "Cell & beams",
                r["heater_setpoint_c"], 15.0, 100.0, 0.5, "°C",
                visible_if={"temperature_model": "Separated"},
                help="Controller setting recorded as metadata; it is not assumed to "
                     "be the vapor temperature."),
            ParamSpec(
                "effective_temp_c", "Effective vapor temperature", "Cell & beams",
                r["effective_temp_c"], 15.0, 100.0, 0.5, "°C",
                visible_if={"temperature_model": "Separated"},
                help="Temperature used for thermal velocity, Doppler averaging, "
                     "transit broadening, and the temperature-dephasing term."),
            ParamSpec(
                "cold_spot_temp_c", "Cold-spot temperature", "Cell & beams",
                r["cold_spot_temp_c"], 15.0, 100.0, 0.5, "°C",
                visible_if={"temperature_model": "Separated"},
                help="Coldest reservoir temperature that fixes the sealed-cell vapor "
                     "pressure. Local density is P_cold/(k_B T_effective)."),
            ParamSpec("coupling_rabi_mhz", "Coupling Rabi (anchor)", "Fields",
                      r["coupling_rabi_mhz"], 0.1, 20.0, 0.1, "MHz", advanced=True,
                      help="Anchor Ω_c/2π at the reference power & waist. The "
                           "Coupling power and Beam diameter sliders scale the "
                           "effective Ω_c around this value."),
            ParamSpec("probe_rabi_ref_mhz", "Probe Rabi (anchor)", "Fields",
                      r["probe_rabi_ref_mhz"], 0.01, 30.0, 0.01, "MHz", advanced=True,
                      help="Model Ω_P/2π at 6 µW and 0.15 mm; scales as sqrt(P)/d. "
                           "This is a calibration parameter, not a measured Gold Rabi rate."),
            ParamSpec("lo_rabi_mhz", "Microwave Rabi Ω_RF", "Fields",
                      r["lo_rabi_mhz"], 0.0, 20.0, 0.1, "MHz",
                      visible_if={"view": "AT electrometry"},
                      help="Rydberg 40D -> 39F dressing Rabi frequency Ω_RF/2π."),
            ParamSpec("mw_detuning_mhz", "Microwave detuning", "Detunings",
                      r["mw_detuning_mhz"], -20.0, 20.0, 0.1, "MHz"),
            ParamSpec("scan_axis", "Optical scan", "Detunings", r["scan_axis"],
                      choices=("Probe detuning", "Coupling detuning"), control="segmented",
                      help="Selects the scanned laser frequency. The other laser's "
                           "detuning is held fixed. The Gold source does not establish "
                           "which optical frequency was scanned."),
            ParamSpec("probe_detuning_mhz", "Fixed probe detuning", "Detunings",
                      r["probe_detuning_mhz"], -100.0, 100.0, 0.1, "MHz",
                      visible_if={"scan_axis": "Coupling detuning"},
                      help="Δp held fixed during a coupling-frequency scan; ignored "
                           "as a fixed value when the probe itself is scanned."),
            ParamSpec("coupling_detuning_mhz", "Fixed coupling detuning", "Detunings",
                      r["coupling_detuning_mhz"], -100.0, 100.0, 0.1, "MHz",
                      visible_if={"scan_axis": "Probe detuning"},
                      help="Δc held fixed during a probe-frequency scan. The bare "
                           "two-photon resonance is Δp + Δc = 0."),
            ParamSpec("mw_frequency_ghz", "Microwave frequency (display only)", "Fields",
                      r["mw_frequency_ghz"], 1.0, 100.0, 0.1, "GHz", recompute=False,
                      help="Display-only metadata; not used in the solve."),
            ParamSpec("rydberg_dephasing_mhz", "Intrinsic dephasing (5S–40D)", "Atomic",
                      r["rydberg_dephasing_mhz"], 0.0, 5.0, 0.01, "MHz", advanced=True,
                      help="Intrinsic (non-transit) 5S–40D coherence broadening "
                           "from laser linewidth etc. The EIT linewidth is the sum "
                           "of this, the transit-time term (from beam diameter), and "
                           "the residual Zeeman term when uncompensated."),
            ParamSpec("temp_dephasing_mhz_per_c", "Temperature dephasing slope", "Atomic",
                      r["temp_dephasing_mhz_per_c"], 0.0, 0.2, 0.001, "MHz/°C",
                      advanced=True,
                      help="Optional phenomenological 5S–40D broadening added above "
                           "the reference temperature. Default 0 preserves the "
                           "Ju et al. reference line shape."),
            ParamSpec(
                "density_dephasing_mhz_per_1e16_m3", "Density dephasing slope",
                "Atomic", r["density_dephasing_mhz_per_1e16_m3"],
                0.0, 2.0, 0.001, "MHz/(10¹⁶ m⁻³)", advanced=True,
                help="Optional phenomenological 5S–40D broadening above the "
                     "20 °C reference density. Temperature and density slopes are "
                     "kept separate; a heating sweep alone may not identify both."),
            ParamSpec("residual_zeeman_mhz", "Residual Zeeman (uncompensated)", "Atomic",
                      r["residual_zeeman_mhz"], 0.0, 3.0, 0.01, "MHz", advanced=True,
                      help="Extra 5S–40D broadening present WITHOUT B-field "
                           "compensation (Ju et al. Fig. 2a blue curve). The EIT "
                           "figure overlays compensated (0) and uncompensated."),
            ParamSpec("rf_dephasing_mhz", "RF dephasing (40D–39F)", "Atomic",
                      r["rf_dephasing_mhz"], 0.0, 5.0, 0.01, "MHz", advanced=True,
                      help="Phenomenological 40D–39F Rydberg–Rydberg coherence "
                           "broadening on the RF-dressed leg."),
            ParamSpec("if_khz", "IF offset", "Detection & scaling", r["if_khz"],
                      1.0, 500.0, 1.0, "kHz", advanced=True, recompute=False,
                      help="LO–SIG beat frequency used by the finite-IF linear-response "
                           "solve. It is not a full time-domain lock-in model."),
            ParamSpec(
                "rf_transition_dipole_ea0", "RF transition dipole", "Detection & scaling",
                r["rf_transition_dipole_ea0"], 1.0, 5000.0, 1.0, "e a₀",
                advanced=True, recompute=False,
                help="40D5/2→39F7/2 matrix element. Default is the ARC 3.10.2 "
                     "stretched-state σ+ value (1326.257 e a0), not an unresolved "
                     "warm-vapor Zeeman average."),
            ParamSpec(
                "rf_angular_factor", "RF angular/polarization factor",
                "Detection & scaling", r["rf_angular_factor"],
                0.01, 1.0, 0.01, "×", advanced=True, recompute=False,
                help="Explicit scalar multiplier for the selected m-state and RF "
                     "polarization. Unity uses the stretched-state dipole; this does "
                     "not replace a full Zeeman/polarization solve."),
            ParamSpec(
                "rf_field_convention", "RF field amplitude", "Detection & scaling",
                r["rf_field_convention"], choices=("RMS", "Peak"),
                advanced=True, recompute=False,
                help="AT splitting is a peak Rabi amplitude. GABES applies the "
                     "sqrt(2) conversion explicitly when reporting an RMS field."),
            ParamSpec(
                "detector_quantum_efficiency", "Detector quantum efficiency",
                "Detection & scaling", r["detector_quantum_efficiency"],
                0.01, 1.0, 0.01, "", advanced=True, recompute=False),
            ParamSpec(
                "detector_path_efficiency", "Signal optical path efficiency",
                "Detection & scaling", r["detector_path_efficiency"],
                0.01, 1.0, 0.01, "×", advanced=True, recompute=False,
                help="Probe fraction delivered to the signal photodiode, including "
                     "the balanced-readout split and post-cell loss. Default 0.5 "
                     "represents an ideal 50:50 split."),
            ParamSpec(
                "detector_reference_power_ratio", "Reference/signal DC power",
                "Detection & scaling", r["detector_reference_power_ratio"],
                0.0, 5.0, 0.05, "×", advanced=True, recompute=False,
                help="Reference-arm DC power relative to the transmitted signal arm. "
                     "The electronic weight is chosen to balance the DC currents."),
            ParamSpec(
                "detector_electronic_noise_pa_sqrt_hz", "Electronics current ASD",
                "Detection & scaling", r["detector_electronic_noise_pa_sqrt_hz"],
                0.0, 100.0, 0.1, "pA/√Hz", advanced=True, recompute=False),
            ParamSpec(
                "detector_rin_per_sqrt_hz", "Optical RIN ASD",
                "Detection & scaling", r["detector_rin_per_sqrt_hz"],
                0.0, 1.0e-3, 1.0e-6, "1/√Hz", advanced=True, recompute=False),
            ParamSpec(
                "detector_rin_correlation", "Signal/reference RIN correlation",
                "Detection & scaling", r["detector_rin_correlation"],
                -1.0, 1.0, 0.01, "", advanced=True, recompute=False),
            ParamSpec(
                "measurement_enbw_hz", "Equivalent noise bandwidth",
                "Detection & scaling", r["measurement_enbw_hz"],
                0.01, 1.0e4, 0.01, "Hz", advanced=True, recompute=False,
                help="Used only to convert the field ASD to an RMS noise floor; "
                     "the displayed sensitivity remains per sqrt(Hz)."),
            ParamSpec(
                "effective_opacity_participation", "Effective opacity participation",
                "Detection & scaling", r["effective_opacity_participation"],
                0.0, 1.0, 0.0001, "×", advanced=True, recompute=False,
                format="%.4f",
                help="Uncalibrated effective multiplier of the atomic optical "
                     "response. The legacy value 0.001 preserves the reference "
                     "spectrum. It is independent of the atom-number fraction."),
            ParamSpec(
                "effective_opacity_relative_uncertainty", "Opacity scenario half-width",
                "Detection & scaling", r["effective_opacity_relative_uncertainty"],
                0.0, 1.0, 0.05, "relative", advanced=True, recompute=False,
                help="Assumed lower/nominal/upper opacity scenarios in the power "
                     "and temperature panels. This is not a measured standard "
                     "uncertainty or a confidence interval."),
            ParamSpec(
                "atom_participation_fraction", "Participating atom fraction",
                "Detection & scaling", r["atom_participation_fraction"],
                0.0, 1.0, 0.01, "×", advanced=True, recompute=False,
                help="Explicit isotope/velocity/Zeeman participation factor used "
                     "only for the reported effective atom number."),
            ParamSpec(
                "beam_overlap_efficiency", "Beam-overlap efficiency",
                "Detection & scaling", r["beam_overlap_efficiency"],
                0.0, 1.0, 0.01, "×", advanced=True, recompute=False,
                help="One-sided geometry factor used only for effective atom number; "
                     "the current lumped OBE already assumes perfect optical overlap."),
            ParamSpec(
                "sam_enabled", "Standard-antenna calibration", "RF source calibration",
                r["sam_enabled"], choices=("Off", "On"), control="segmented",
                recompute=False, advanced=True,
                help="Convert source power to a far-field RF amplitude with an "
                     "explicit uncertainty budget. This does not model horn near fields."),
            ParamSpec(
                "sam_source_power_dbm", "Source power", "RF source calibration",
                r["sam_source_power_dbm"], -120.0, 30.0, 0.1, "dBm",
                recompute=False, advanced=True, visible_if={"sam_enabled": "On"}),
            ParamSpec(
                "sam_antenna_gain_dbi", "Antenna gain", "RF source calibration",
                r["sam_antenna_gain_dbi"], -20.0, 50.0, 0.1, "dBi",
                recompute=False, advanced=True, visible_if={"sam_enabled": "On"}),
            ParamSpec(
                "sam_distance_m", "Antenna distance", "RF source calibration",
                r["sam_distance_m"], 0.001, 20.0, 0.001, "m",
                recompute=False, advanced=True, visible_if={"sam_enabled": "On"}),
            ParamSpec(
                "sam_cable_loss_db", "Cable loss", "RF source calibration",
                r["sam_cable_loss_db"], 0.0, 100.0, 0.1, "dB",
                recompute=False, advanced=True, visible_if={"sam_enabled": "On"}),
            ParamSpec(
                "sam_additional_loss_db", "Additional loss", "RF source calibration",
                r["sam_additional_loss_db"], 0.0, 100.0, 0.1, "dB",
                recompute=False, advanced=True, visible_if={"sam_enabled": "On"}),
            ParamSpec(
                "sam_field_correction", "Field correction", "RF source calibration",
                r["sam_field_correction"], 0.001, 10.0, 0.001, "×",
                recompute=False, advanced=True, visible_if={"sam_enabled": "On"}),
            ParamSpec(
                "sam_source_power_std_db", "Source-power standard uncertainty",
                "RF source calibration", r["sam_source_power_std_db"],
                0.0, 20.0, 0.01, "dB", recompute=False, advanced=True,
                visible_if={"sam_enabled": "On"}),
            ParamSpec(
                "sam_antenna_gain_std_db", "Gain standard uncertainty",
                "RF source calibration", r["sam_antenna_gain_std_db"],
                0.0, 20.0, 0.01, "dB", recompute=False, advanced=True,
                visible_if={"sam_enabled": "On"}),
            ParamSpec(
                "sam_cable_loss_std_db", "Cable-loss standard uncertainty",
                "RF source calibration", r["sam_cable_loss_std_db"],
                0.0, 20.0, 0.01, "dB", recompute=False, advanced=True,
                visible_if={"sam_enabled": "On"}),
            ParamSpec(
                "sam_additional_loss_std_db", "Additional-loss standard uncertainty",
                "RF source calibration", r["sam_additional_loss_std_db"],
                0.0, 20.0, 0.01, "dB", recompute=False, advanced=True,
                visible_if={"sam_enabled": "On"}),
            ParamSpec(
                "sam_distance_std_m", "Distance standard uncertainty",
                "RF source calibration", r["sam_distance_std_m"],
                0.0, 5.0, 0.001, "m", recompute=False, advanced=True,
                visible_if={"sam_enabled": "On"}),
            ParamSpec(
                "sam_field_correction_std", "Field-correction standard uncertainty",
                "RF source calibration", r["sam_field_correction_std"],
                0.0, 5.0, 0.001, "×", recompute=False, advanced=True,
                visible_if={"sam_enabled": "On"}),
            ParamSpec(
                "sam_antenna_max_dimension_m", "Antenna maximum dimension",
                "RF source calibration", r["sam_antenna_max_dimension_m"],
                0.001, 5.0, 0.001, "m", recompute=False, advanced=True,
                visible_if={"sam_enabled": "On"},
                help="Used with the microwave frequency to report r/(2D²/lambda)."),
            ParamSpec("doppler", "Doppler treatment", "Numerics", r["doppler"],
                      choices=("off", "on"), advanced=True,
                      help="off: a single velocity-class approximation. on: exact "
                           "Maxwell average with the counter-propagating optical "
                           "wave vectors. Neither mode is independently calibrated "
                           "to the Gold experiment by selecting this control."),
        ]

    def _defaults(self, view):
        vals = dict(self._REF)
        vals["view"] = view
        if view == "EIT":
            vals["lo_rabi_mhz"] = 0.0
        return vals

    def model_comparison(self, params, *, resolve_raw, current_raw=None, current_view=None):
        from .rydberg_comparison import spectral_comparison

        return spectral_comparison(self, params, resolve_raw=resolve_raw,
                                   current_raw=current_raw, current_view=current_view)

    def recommended_defaults(self, params):
        return {
            "EIT": self._defaults("EIT"),
            "AT electrometry": self._defaults("AT electrometry"),
        }

    def _atom(self, ground_deph, rf_deph):
        """Cascade with completely positive nested diagonal dephasing.

        The optical-Rydberg projector gives gamma_02=ground_deph and the
        additional RF projector gives gamma_23=rf_deph. Consequently the
        0--3 coherence has their sum; independent pairwise damping would not
        in general generate a completely positive map.
        """
        sk = _cascade_skeleton()
        ground_deph, rf_deph = float(ground_deph), float(rf_deph)
        if not all(np.isfinite(rate) and rate >= 0 for rate in (ground_deph, rf_deph)):
            raise ValueError("Dephasing rates must be finite and nonnegative")
        optical = np.sqrt(2.0 * ground_deph) * np.diag([0., 0., 1., 1.])
        rf = np.sqrt(2.0 * rf_deph) * np.diag([0., 0., 0., 1.])
        return atoms.AtomModel(
            name="rydberg_eit_85rb",
            n_levels=4,
            labels=sk["labels"],
            ground=(0,),
            excited=(1, 2, 3),
            decay=sk["decay"],
            dephasing=(),
            collapse_ops=(optical, rf),
            doppler_levels=sk["doppler_levels"],
            doppler_ratios=sk["doppler_ratios"],
        )

    def _coupling_rabi(self, params):
        """Effective coupling Rabi Ω_c/2π [MHz] from the 481 nm beam power and
        waist, FWM-style (fwm.py `_fields_from_params`): Ω ∝ E ∝ √(intensity) =
        √(P / area), area ∝ d². Anchored to the reference operating point, so at
        (P_ref, d_ref) it reproduces the fitted Coupling Rabi anchor exactly."""
        oc_ref = float(params.get("coupling_rabi_mhz", self._REF["coupling_rabi_mhz"]))
        p = float(params.get("coupling_power_mw", self._REF["coupling_power_mw"]))
        d = float(params.get("beam_diameter_mm", self._REF["beam_diameter_mm"]))
        p_ref = self._REF["coupling_power_mw"]
        d_ref = self._REF["beam_diameter_mm"]
        return beam.anchored_rabi_mhz(
            oc_ref, p, p_ref, diameter=d, ref_diameter=d_ref)

    def _probe_rabi(self, params):
        """Weak-probe drive Ω_P/2π [MHz] from probe power & waist, anchored at the
        reference like the coupling: Ω_P ∝ √(P)/d. Finite Ω_P makes the OBE
        saturate, so raising probe power broadens the EIT line (Ju et al. Fig 2b)."""
        p = float(params.get("probe_power_uw", self._REF["probe_power_uw"]))
        d = float(params.get("beam_diameter_mm", self._REF["beam_diameter_mm"]))
        p_ref = self._REF["probe_power_uw"]
        d_ref = self._REF["beam_diameter_mm"]
        return beam.anchored_rabi_mhz(
            float(params.get("probe_rabi_ref_mhz", self._REF["probe_rabi_ref_mhz"])),
            p, p_ref, diameter=d, ref_diameter=d_ref)

    @classmethod
    def _scan_axis(cls, params):
        axis = str(params.get("scan_axis", cls._REF["scan_axis"]))
        if axis not in ("Probe detuning", "Coupling detuning"):
            raise ValueError("scan_axis must select Probe detuning or Coupling detuning")
        return axis

    @classmethod
    def _optical_detunings_at_scan(cls, params, scan_rad_s):
        """Physical (Delta_p, Delta_c) in rad/s for one optical scan coordinate."""
        dp = float(params.get("probe_detuning_mhz", cls._REF["probe_detuning_mhz"])) * MHZ
        dc = float(params.get("coupling_detuning_mhz", cls._REF["coupling_detuning_mhz"])) * MHZ
        if cls._scan_axis(params) == "Probe detuning":
            dp = float(scan_rad_s)
        else:
            dc = float(scan_rad_s)
        if not np.isfinite(dp) or not np.isfinite(dc):
            raise ValueError("Optical detunings must be finite")
        return dp, dc

    @classmethod
    def _resonance_scan_mhz(cls, params):
        fixed_key = ("coupling_detuning_mhz" if cls._scan_axis(params) == "Probe detuning"
                     else "probe_detuning_mhz")
        center = -float(params.get(fixed_key, cls._REF[fixed_key]))
        if not np.isfinite(center):
            raise ValueError("Fixed optical detuning must be finite")
        return center

    @classmethod
    def _scan_label(cls, params):
        return cls._scan_axis(params) + " [MHz]"

    def _temperature_state(self, params):
        """Return heater/effective/cold-spot temperatures without conflating them.

        The legacy/landing configuration remains ``Linked`` so existing callers
        that only pass ``temp_c`` reproduce the old density and linewidth.  The
        separated path is the one intended for heated-cell analysis.
        """
        mode = str(params.get(
            "temperature_model", self._REF["temperature_model"]))
        if mode == "Separated":
            return rydberg_experiment.resolve_cell_temperature(
                float(params.get(
                    "heater_setpoint_c", self._REF["heater_setpoint_c"])),
                effective_vapor_temp_c=float(params.get(
                    "effective_temp_c", self._REF["effective_temp_c"])),
                cold_spot_temp_c=float(params.get(
                    "cold_spot_temp_c", self._REF["cold_spot_temp_c"])),
            )
        linked = float(params.get("temp_c", self._REF["temp_c"]))
        return rydberg_experiment.resolve_cell_temperature(
            linked, effective_vapor_temp_c=linked, cold_spot_temp_c=linked)

    @staticmethod
    def _local_density_from_temperature_state(state):
        """Sealed-cell density at the illuminated effective temperature.

        The cold spot fixes the saturated pressure.  At the beam location the
        ideal-gas density is therefore n = n_cold*T_cold/T_effective.  Linked
        mode reduces exactly to the previous saturated-density calculation.
        """
        cold_k = state.cold_spot_temp_c + 273.15
        effective_k = state.effective_vapor_temp_c + 273.15
        return state.saturated_cold_spot_density_m3() * cold_k / effective_k

    def _transit_rate_mhz(self, params):
        """Transit-time broadening of the 5S–40D coherence /2π [MHz]: an atom
        crosses the beam in ~d/v_mp, so the coherence decays at TRANSIT_FACTOR·
        v_mp/d. Ties the beam diameter (and temperature) to the EIT linewidth —
        the transit-limited regime Ju et al. report."""
        T = self._temperature_state(params).effective_vapor_temp_c + 273.15
        d = float(params.get("beam_diameter_mm", self._REF["beam_diameter_mm"]))
        return beam.transit_broadening_mhz(
            T, d, mass=constants.MASS_85RB, factor=TRANSIT_FACTOR)

    def _scan_chi(self, atom, h_of, scan, probe, vapor):
        """χ̄(scan) for one dephasing configuration: the v = 0 class when `vapor` is
        None, else the exact Maxwell average for vapor = (T, mass, k_vec). The probe
        detuning lives in H. Shared by the compensated / uncompensated solves and the
        Fig. 2(b) probe-power sweep."""
        affine = affine_doppler.AffineScan(atom, h_of, (1, 0), probe)
        return affine.stationary(scan) if vapor is None else affine.maxwell_mean(scan, *vapor)

    @staticmethod
    def _inhomogeneous(chi, scan, fwhm_mhz):
        """Gaussian inhomogeneous broadening of χ̄ over the probe detuning axis:
        a spread of resonance shifts (residual Zeeman) of FWHM `fwhm_mhz` averages
        χ(Δ − δ) over a Gaussian δ-distribution. Returns χ̄ unchanged for ~0 FWHM."""
        if fwhm_mhz <= 1e-6:
            return chi
        dscan = (scan[1] - scan[0]) / MHZ          # MHz per sample (uniform grid)
        sigma = fwhm_mhz / (2 * np.sqrt(2 * np.log(2)))
        half = int(np.ceil(4 * sigma / dscan))
        offs = np.arange(-half, half + 1) * dscan
        g = np.exp(-0.5 * (offs / sigma) ** 2)
        g /= g.sum()
        return np.convolve(chi, g, mode="same")

    def _scan(self, params):
        view = params.get("view", "AT electrometry")
        lo = float(params.get("lo_rabi_mhz", self._REF["lo_rabi_mhz"]))
        # Size the window off the EFFECTIVE Ω_c (the same power/waist-scaled value
        # that enters h_of), not the raw anchor — otherwise raising coupling power
        # or shrinking the beam widens the real feature past a pinned window and
        # the spectrum (and its slope/linewidth metrics) get clipped.
        oc = self._coupling_rabi(params)
        if view == "AT electrometry" and lo > 0:
            half = max(10.0, 3.0 * max(lo, oc))
        else:
            half = max(8.0, 3.0 * oc)
        center = self._resonance_scan_mhz(params)
        return (np.linspace(-half, half, 801) + center) * MHZ

    def compute(self, params):
        line = _probe_line()
        k_vec = line["k_vec"]
        dipole = line["dipole"]
        nu0 = line["nu0"]
        temperature_state = self._temperature_state(params)
        effective_temp_c = temperature_state.effective_vapor_temp_c
        cold_spot_temp_c = temperature_state.cold_spot_temp_c
        T = effective_temp_c + 273.15
        N = self._local_density_from_temperature_state(temperature_state)

        # 5S-40D dephasing budget: intrinsic (laser etc.) + transit-time (beam
        # diameter) is always present; the residual Zeeman term is added only on
        # the uncompensated curve (B-field compensation off).
        intrinsic_base = float(params.get(
            "rydberg_dephasing_mhz", self._REF["rydberg_dephasing_mhz"]))
        temp_slope = float(params.get(
            "temp_dephasing_mhz_per_c", self._REF["temp_dephasing_mhz_per_c"]))
        density_slope = float(params.get(
            "density_dephasing_mhz_per_1e16_m3",
            self._REF["density_dephasing_mhz_per_1e16_m3"]))
        temp_extra = max(effective_temp_c - self._REF["temp_c"], 0.0) * max(
            temp_slope, 0.0)
        reference_density = species.number_density(
            species.RB85, self._REF["temp_c"] + 273.15)
        density_extra = max(N - reference_density, 0.0) / 1.0e16 * max(
            density_slope, 0.0)
        intrinsic = intrinsic_base + temp_extra + density_extra
        transit = self._transit_rate_mhz(params)
        zeeman = float(params.get(
            "residual_zeeman_mhz", self._REF["residual_zeeman_mhz"]))
        rf_deph = float(params.get(
            "rf_dephasing_mhz", self._REF["rf_dephasing_mhz"])) * MHZ
        gd = (intrinsic + transit) * MHZ

        probe = self._probe_rabi(params) * MHZ
        Oc = self._coupling_rabi(params) * MHZ
        Olo = float(params.get("lo_rabi_mhz", self._REF["lo_rabi_mhz"])) * MHZ
        if params.get("view", "EIT") == "EIT":
            Olo = 0.0
        h_fields = dict(probe_rabi_mhz=probe / MHZ, coupling_rabi_mhz=Oc / MHZ,
                        lo_rabi_mhz=Olo / MHZ)

        def h_of(s):
            return self._hamiltonian_at_scan_detuning(h_fields, params, s)

        scan = self._scan(params)
        # Doppler off (default): single static class (kv = 0) — the suppressed-
        # residual, transit-limited regime. Doppler on: Maxwell-average the
        # residual two-photon shift carried by atom.S_v (per-level k ratios).
        vapor = (T, constants.MASS_85RB, k_vec) if params.get("doppler", "off") == "on" else None

        # Compensated EIT: one solve. Uncompensated: residual Zeeman scatters the
        # m-sublevel shifts, so the EIT line is INHOMOGENEOUSLY broadened — model
        # it as a Gaussian (FWHM = residual_zeeman) convolution of the compensated
        # susceptibility over the probe detuning. (A homogeneous dephasing knob
        # only lowers the peak in this Ω_c-limited regime; it does not widen it.)
        chi = self._scan_chi(self._atom(gd, rf_deph), h_of, scan, probe, vapor)
        chi_uncomp = self._inhomogeneous(chi, scan, zeeman)

        return dict(
            scan=scan,
            scan_axis=self._scan_axis(params),
            scan_axis_label=self._scan_label(params),
            fixed_probe_detuning_mhz=float(params.get("probe_detuning_mhz", self._REF["probe_detuning_mhz"])),
            fixed_coupling_detuning_mhz=float(params.get("coupling_detuning_mhz", self._REF["coupling_detuning_mhz"])),
            two_photon_resonance_scan_mhz=self._resonance_scan_mhz(params),
            scan_axis_source="model control; Gold scanned-laser identity not established",
            dephasing_model="correlated_phase_lindblad",
            chi_bar=chi,
            chi_bar_uncomp=chi_uncomp,
            N=N,
            T=T,
            heater_setpoint_c=temperature_state.heater_setpoint_c,
            effective_temp_c=effective_temp_c,
            cold_spot_temp_c=cold_spot_temp_c,
            temperature_model=str(params.get(
                "temperature_model", self._REF["temperature_model"])),
            temperature_effective_source=temperature_state.effective_source,
            temperature_cold_spot_source=temperature_state.cold_spot_source,
            L=float(params.get("cell_mm", self._REF["cell_mm"])) * 1e-3,
            ls=0.001,
            # Compatibility alias only; live propagation uses the named factor.
            ls_status="deprecated legacy nominal combined opacity alias",
            atomic_line_strength=1.0,
            k_vec=k_vec,
            omega0=2 * np.pi * nu0,
            dipole=dipole,
            gamma_mhz=line["gamma_mhz"],
            coupling_rabi_mhz=Oc / MHZ,
            probe_rabi_mhz=probe / MHZ,
            probe_rabi_ref_mhz=float(params.get("probe_rabi_ref_mhz", self._REF["probe_rabi_ref_mhz"])),
            lo_rabi_mhz=Olo / MHZ,
            rydberg_dephasing_mhz=intrinsic_base,
            temperature_dephasing_mhz=temp_extra,
            density_dephasing_mhz=density_extra,
            reference_density_m3=reference_density,
            transit_mhz=transit,
            residual_zeeman_mhz=zeeman,
            rf_dephasing_mhz=rf_deph / MHZ,
            mw_frequency_ghz=float(params.get(
                "mw_frequency_ghz", self._REF["mw_frequency_ghz"])),
            probe_power_uw=float(params.get("probe_power_uw", self._REF["probe_power_uw"])),
            beam_diameter_mm=float(params.get(
                "beam_diameter_mm", self._REF["beam_diameter_mm"])),
        )

    @staticmethod
    def _transparency_maxima(x, y, window_mhz, center_mhz=0.0):
        """Local transmission maxima (pos [MHz], height) within ±window of line
        centre — the dressed AT/EIT transparency peaks."""
        return [(float(x[i]), float(y[i])) for i in range(1, y.size - 1)
                if abs(x[i] - center_mhz) <= window_mhz and y[i] > y[i - 1] and y[i] >= y[i + 1]]

    @staticmethod
    def _opacity_scale(raw, params):
        """Apply the atomic strength and effective participation exactly once."""
        return (float(raw.get("atomic_line_strength", 1.0))
                * opacity.opacity_spec(params)["value"])

    def _transmission(self, chi, raw, params):
        """(x [MHz], Beer-Lambert transmission, χ_phys) for one coherence array."""
        x = raw["scan"] / MHZ
        alpha, xphys = observables.absorption_coefficient(
            chi, raw["k_vec"], raw["N"], dipole=raw["dipole"],
            line_strength=self._opacity_scale(raw, params))
        # Cell length only scales αL here (recompute=False navigate-only knob), so
        # read it from the live params, not the cached raw, or a cell-length change
        # would not update the transmission.
        L = float(params.get("cell_mm", self._REF["cell_mm"])) * 1e-3
        return x, observables.transmission(alpha, L), xphys

    def _rf_coupling(self, params):
        """Explicit 40D--39F field/Rabi calibration used by AT and SIG paths."""
        dipole_ea0 = float(params.get(
            "rf_transition_dipole_ea0", self._REF["rf_transition_dipole_ea0"]))
        angular_factor = float(params.get(
            "rf_angular_factor", self._REF["rf_angular_factor"]))
        convention = str(params.get(
            "rf_field_convention", self._REF["rf_field_convention"])).lower()
        return electrometry.RFDipoleCoupling(
            dipole_ea0 * constants.ELEMENTARY_CHARGE * BOHR_RADIUS_M,
            angular_factor=angular_factor,
            field_amplitude_convention=convention,
        )

    def _sam_calibration(self, params):
        """Optional standard-antenna field and first-order uncertainty budget."""
        if str(params.get("sam_enabled", self._REF["sam_enabled"])) != "On":
            return None
        convention = str(params.get(
            "rf_field_convention", self._REF["rf_field_convention"])).lower()
        return rydberg_experiment.sam_field_calibration(
            float(params.get(
                "sam_source_power_dbm", self._REF["sam_source_power_dbm"])),
            float(params.get(
                "sam_antenna_gain_dbi", self._REF["sam_antenna_gain_dbi"])),
            float(params.get("sam_distance_m", self._REF["sam_distance_m"])),
            cable_loss_db=float(params.get(
                "sam_cable_loss_db", self._REF["sam_cable_loss_db"])),
            additional_loss_db=float(params.get(
                "sam_additional_loss_db", self._REF["sam_additional_loss_db"])),
            field_correction=float(params.get(
                "sam_field_correction", self._REF["sam_field_correction"])),
            source_power_std_db=float(params.get(
                "sam_source_power_std_db", self._REF["sam_source_power_std_db"])),
            antenna_gain_std_db=float(params.get(
                "sam_antenna_gain_std_db", self._REF["sam_antenna_gain_std_db"])),
            cable_loss_std_db=float(params.get(
                "sam_cable_loss_std_db", self._REF["sam_cable_loss_std_db"])),
            additional_loss_std_db=float(params.get(
                "sam_additional_loss_std_db",
                self._REF["sam_additional_loss_std_db"])),
            distance_std_m=float(params.get(
                "sam_distance_std_m", self._REF["sam_distance_std_m"])),
            field_correction_std=float(params.get(
                "sam_field_correction_std", self._REF["sam_field_correction_std"])),
            amplitude_convention=convention,
            frequency_hz=float(params.get(
                "mw_frequency_ghz", self._REF["mw_frequency_ghz"])) * 1.0e9,
            antenna_max_dimension_m=float(params.get(
                "sam_antenna_max_dimension_m",
                self._REF["sam_antenna_max_dimension_m"])),
        )

    def _effective_atom_number(self, raw, params):
        """Gaussian probe/coupling overlap estimate for the illuminated atoms."""
        radius_m = float(raw["beam_diameter_mm"]) * 0.5e-3
        return rydberg_experiment.effective_atom_number(
            raw["N"],
            length_m=float(params.get("cell_mm", self._REF["cell_mm"])) * 1e-3,
            probe_radius_m=radius_m,
            coupling_radius_m=radius_m,
            participation_fraction=float(params.get(
                "atom_participation_fraction",
                self._REF["atom_participation_fraction"])),
            overlap_efficiency=float(params.get(
                "beam_overlap_efficiency",
                self._REF["beam_overlap_efficiency"])),
        )

    @staticmethod
    def _absorption_quadrature_operator():
        """Hermitian O with Tr(O rho) = Im(rho_10) for the probe coherence."""
        operator = np.zeros((4, 4), dtype=complex)
        operator[0, 1] = -0.5j
        operator[1, 0] = +0.5j
        return operator

    @classmethod
    def _hamiltonian_at_scan_detuning(cls, raw, params, detuning_rad_s):
        """One shared optical/RF Hamiltonian for stationary and finite-IF paths."""
        dp, dc = cls._optical_detunings_at_scan(params, detuning_rad_s)
        H = np.zeros((4, 4), dtype=complex)
        H[1, 1] = -dp
        H[2, 2] = -(dp + dc)
        H[3, 3] = -(dp + dc) - float(params.get(
            "mw_detuning_mhz", 0.0)) * MHZ
        H[0, 1] = H[1, 0] = raw["probe_rabi_mhz"] * MHZ / 2.0
        H[1, 2] = H[2, 1] = raw["coupling_rabi_mhz"] * MHZ / 2.0
        H[2, 3] = H[3, 2] = raw["lo_rabi_mhz"] * MHZ / 2.0
        return H

    @classmethod
    def _hamiltonian_at_probe_detuning(cls, raw, params, detuning_rad_s):
        """Compatibility alias: argument follows the selected optical scan axis."""
        return cls._hamiltonian_at_scan_detuning(raw, params, detuning_rad_s)

    def _superheterodyne_atomic_phasors(
            self, raw, params, x, transmission, static_indices=None):
        """Finite-IF atomic response before RF and detector calibration."""
        ground_dephasing = (
            raw["rydberg_dephasing_mhz"]
            + raw["temperature_dephasing_mhz"]
            + raw["density_dephasing_mhz"]
            + raw["transit_mhz"]
        ) * MHZ
        atom = self._atom(ground_dephasing, raw["rf_dephasing_mhz"] * MHZ)
        omega_if = max(float(params.get(
            "if_khz", self._REF["if_khz"])), 1.0e-12) * 1.0e3 * 2.0 * np.pi
        operator = self._absorption_quadrature_operator()

        if params.get("doppler", self._REF["doppler"]) == "on":
            if static_indices is None:
                static_indices = self._static_slope_indices(x, transmission)
            indices = np.asarray(static_indices, dtype=int)
            affine = affine_doppler.AffineScan(
                atom, lambda s: self._hamiltonian_at_scan_detuning(raw, params, s),
                (1, 0), raw["probe_rabi_mhz"] * MHZ)
            velocity, weights = affine.quadrature_classes(
                raw["scan"][indices], raw["T"], constants.MASS_85RB, raw["k_vec"])
            kv = raw["k_vec"] * velocity
            absorption_phasors = np.empty(indices.size, dtype=complex)
            for slot, index in enumerate(indices):
                H = self._hamiltonian_at_scan_detuning(
                    raw, params, raw["scan"][index])
                base = core.build_liouvillian(H, atom)
                class_phasors = np.empty(kv.size, dtype=complex)
                for first in range(0, kv.size, SUPERHET_CLASS_CHUNK):
                    part = slice(first, first + SUPERHET_CLASS_CHUNK)
                    liouvillian = (base[None, :, :]
                                   + kv[part, None, None] * atom.S_v[None, :, :])
                    class_phasors[part] = electrometry.weak_signal_response_from_liouvillian(
                        liouvillian, 4, omega_if, signal_transition=(2, 3)
                    ).real_observable_phasor_per_angular_rabi(operator)
                absorption_phasors[slot] = electrometry.coherent_weighted_average(
                    class_phasors, weights)
            optimization_scope = "static-slope points; coherent Doppler average"
        else:
            indices = np.arange(raw["scan"].size, dtype=int)
            liouvillian = np.asarray([
                core.build_liouvillian(
                    self._hamiltonian_at_scan_detuning(raw, params, detuning), atom)
                for detuning in raw["scan"]
            ])
            response = electrometry.weak_signal_response_from_liouvillian(
                liouvillian, 4, omega_if, signal_transition=(2, 3))
            absorption_phasors = (
                response.real_observable_phasor_per_angular_rabi(operator))
            optimization_scope = "full optical-detuning scan"
        return indices, absorption_phasors, optimization_scope

    @staticmethod
    def _mirror_maxima(values):
        """Indices of the largest `values`, ascending. A symmetric line has two mirror
        maxima that only rounding separates; keeping both stops rounding from picking
        a side of the line."""
        top = int(np.nanargmax(values))
        if values[top] <= 0.0:
            return np.array([top])
        return np.flatnonzero(values >= (1.0 - MIRROR_TIE_RTOL) * values[top])

    @classmethod
    def _static_slope_indices(cls, x, transmission):
        """Scan points of the steepest static slope; the readout keeps the most
        sensitive of them."""
        return cls._mirror_maxima(np.abs(np.gradient(transmission, x)))

    def _cached_superheterodyne_atomic_phasors(
            self, raw, params, x, transmission):
        doppler_mode = str(params.get("doppler", self._REF["doppler"]))
        static_indices = (
            self._static_slope_indices(x, transmission)
            if doppler_mode == "on" else None
        )
        key = (
            doppler_mode,
            None if static_indices is None else tuple(static_indices.tolist()),
            float(params.get("if_khz", self._REF["if_khz"])),
            float(params.get("mw_detuning_mhz", self._REF["mw_detuning_mhz"])),
            self._scan_axis(params),
            float(params.get("probe_detuning_mhz", self._REF["probe_detuning_mhz"])),
            float(params.get("coupling_detuning_mhz", self._REF["coupling_detuning_mhz"])),
            float(raw["probe_rabi_mhz"]), float(raw["coupling_rabi_mhz"]),
            float(raw["lo_rabi_mhz"]), float(raw["rydberg_dephasing_mhz"]),
            float(raw["temperature_dephasing_mhz"]),
            float(raw["density_dephasing_mhz"]), float(raw["transit_mhz"]),
            float(raw["rf_dephasing_mhz"]), int(raw["scan"].size),
            float(raw["scan"][0]), float(raw["scan"][-1]),
            float(raw["T"]) if doppler_mode == "on" else None,
        )
        with self._superhet_atomic_cache_lock:
            cached = getattr(self, "_superhet_atomic_cache", None)
            if not isinstance(cached, OrderedDict):
                cached = OrderedDict()
                self._superhet_atomic_cache = cached
            try:
                result = cached[key]
            except KeyError:
                pass
            else:
                cached.move_to_end(key)
                return result
        result = self._superheterodyne_atomic_phasors(
            raw, params, x, transmission, static_indices=static_indices)
        with self._superhet_atomic_cache_lock:
            cached = self._superhet_atomic_cache
            try:
                existing = cached[key]
            except KeyError:
                cached[key] = result
                if len(cached) > SUPERHET_ATOMIC_CACHE_SIZE:
                    cached.popitem(last=False)
                return result
            cached.move_to_end(key)
            return existing

    def _superheterodyne_readout(
            self, raw, params, x, transmission, atomic_response=None):
        """Finite-IF atomic response with RF and detector calibration."""
        if raw["lo_rabi_mhz"] <= 0.0 or self._opacity_scale(raw, params) == 0.0:
            return None

        response_reused = atomic_response is not None
        if atomic_response is None:
            atomic_response = self._cached_superheterodyne_atomic_phasors(
                raw, params, x, transmission)
        indices, absorption_phasors, optimization_scope = atomic_response
        if response_reused and str(params.get("doppler", self._REF["doppler"])) == "on":
            optimization_scope = "fixed nominal static-gradient sampled point (opacity scenario)"

        # chi_bar = rho_10 / Omega_probe.  Linearizing
        # T=exp[-k Im(chi_phys)L] gives the complex transmission phasor below.
        chi_scale = observables.chi_phys(
            1.0, raw["N"], dipole=raw["dipole"],
            line_strength=self._opacity_scale(raw, params))
        cell_length_m = float(params.get(
            "cell_mm", self._REF["cell_mm"])) * 1.0e-3
        transmission_phasor_per_rabi = (
            -cell_length_m
            * transmission[indices]
            * raw["k_vec"]
            * chi_scale
            * absorption_phasors
            / (raw["probe_rabi_mhz"] * MHZ)
        )
        wavelength_m = 2.0 * np.pi * constants.C_LIGHT / raw["omega0"]
        quantum_efficiency = float(params.get(
            "detector_quantum_efficiency",
            self._REF["detector_quantum_efficiency"]))
        detector_responsivity = electrometry.photodiode_responsivity_a_per_w(
            quantum_efficiency, wavelength_m)
        path_efficiency = float(params.get(
            "detector_path_efficiency",
            self._REF["detector_path_efficiency"]))
        incident_power_w = raw["probe_power_uw"] * 1.0e-6
        current_scale = detector_responsivity * incident_power_w * path_efficiency
        rf_coupling = self._rf_coupling(params)
        reference_ratio = float(params.get(
            "detector_reference_power_ratio",
            self._REF["detector_reference_power_ratio"]))
        reference_weight = 1.0 / reference_ratio if reference_ratio > 0.0 else 0.0
        current_phasors = electrometry.current_responsivity_from_atomic_phasor(
            transmission_phasor_per_rabi, current_scale, rf_coupling)
        candidates = []
        total_field_asd = np.full(indices.size, np.inf, dtype=float)
        for local_index, (global_index, current_phasor) in enumerate(
                zip(indices, current_phasors)):
            responsivity = float(abs(current_phasor))
            if not np.isfinite(responsivity) or responsivity <= 0.0:
                candidates.append(None)
                continue
            signal_power_w = (
                incident_power_w * path_efficiency * transmission[global_index])
            signal_channel = electrometry.PhotodiodeChannel(
                signal_power_w, detector_responsivity)
            reference_channel = (
                electrometry.PhotodiodeChannel(
                    signal_power_w * reference_ratio, detector_responsivity)
                if reference_ratio > 0.0 else None)
            detector = electrometry.BalancedDetector(
                signal=signal_channel,
                reference=reference_channel,
                reference_weight=reference_weight,
                electronic_noise_current_asd_a_per_sqrt_hz=float(params.get(
                    "detector_electronic_noise_pa_sqrt_hz",
                    self._REF["detector_electronic_noise_pa_sqrt_hz"])) * 1.0e-12,
                relative_intensity_noise_per_sqrt_hz=float(params.get(
                    "detector_rin_per_sqrt_hz",
                    self._REF["detector_rin_per_sqrt_hz"])),
                rin_correlation=float(params.get(
                    "detector_rin_correlation",
                    self._REF["detector_rin_correlation"])),
            )
            noise = electrometry.balanced_detector_noise(detector)
            sensitivity = electrometry.electrometry_sensitivity(
                noise, responsivity)
            total_field_asd[local_index] = (
                sensitivity.total_field_asd_v_m_per_sqrt_hz)
            candidates.append((
                current_phasor, signal_power_w, reference_channel,
                noise, sensitivity,
            ))
        local_best = int(np.nanargmin(total_field_asd))
        best_index = int(indices[local_best])
        best = candidates[local_best]
        if best is None:
            return None
        current_phasor, signal_power_w, reference_channel, noise, sensitivity = best
        enbw_hz = float(params.get(
            "measurement_enbw_hz", self._REF["measurement_enbw_hz"]))
        optimum_dp, optimum_dc = self._optical_detunings_at_scan(
            params, float(raw["scan"][best_index]))
        return dict(
            response_detuning_mhz=x[indices],
            current_phasors_a_per_v_m=current_phasors,
            psn_field_asd_v_m_per_sqrt_hz=np.asarray([
                candidate[4].psn_field_asd_v_m_per_sqrt_hz if candidate is not None else np.inf
                for candidate in candidates]),
            technical_field_asd_v_m_per_sqrt_hz=np.asarray([
                candidate[4].technical_field_asd_v_m_per_sqrt_hz if candidate is not None else np.inf
                for candidate in candidates]),
            scan_axis=self._scan_axis(params),
            transmission_phasor_per_angular_rabi=transmission_phasor_per_rabi,
            optimum_index=best_index,
            optimum_detuning_mhz=float(x[best_index]),
            optimum_scan_detuning_mhz=float(x[best_index]),
            optimum_probe_detuning_mhz=optimum_dp / MHZ,
            optimum_coupling_detuning_mhz=optimum_dc / MHZ,
            optimum_two_photon_detuning_mhz=(optimum_dp + optimum_dc) / MHZ,
            optimum_transmission=float(transmission[best_index]),
            response_phase_deg=float(np.angle(current_phasor, deg=True)),
            current_responsivity_a_per_v_m=float(abs(current_phasor)),
            detector_responsivity_a_per_w=float(detector_responsivity),
            signal_detector_power_w=float(signal_power_w),
            reference_detector_power_w=(
                0.0 if reference_channel is None
                else float(reference_channel.optical_power_w)),
            noise=noise,
            sensitivity=sensitivity,
            total_field_asd_v_m_per_sqrt_hz=total_field_asd,
            rms_field_noise_v_m=float(electrometry.asd_to_rms(
                sensitivity.total_field_asd_v_m_per_sqrt_hz, enbw_hz)),
            enbw_hz=enbw_hz,
            rf_coupling=rf_coupling,
            optimization_scope=optimization_scope,
            atomic_response=atomic_response,
        )

    @staticmethod
    def _eit_features(x, T_trans, center_mhz=0.0):
        """(FWHM [MHz], peak contrast) of the central EIT transparency feature.

        Contrast = T(centre) − absorptive floor in a ±5 MHz window. The FWHM is
        measured to the half-contrast level with linear interpolation of the two
        crossings (sub-sample accurate — a bare grid walk stair-steps the
        probe-power sweep). The transmitted-signal amplitude in Fig. 2(b) scales
        as probe power × contrast."""
        ic = int(np.argmin(np.abs(x - center_mhz)))
        peak = float(T_trans[ic])
        win = np.abs(x - center_mhz) <= 5.0
        floor = float(np.min(T_trans[win])) if win.any() else float(np.min(T_trans))
        contrast = peak - floor
        if contrast <= 0:
            return float("nan"), 0.0
        half = floor + 0.5 * contrast
        if ic == 0 or ic == x.size - 1:
            return float("nan"), contrast
        right, _ = level_crossing(x, T_trans, ic, +1, half)
        left, _ = level_crossing(x, T_trans, ic, -1, half)
        return abs(right - left), contrast

    def _readout(self, raw, params, atomic_response=None):
        """Cheap transmission/dispersion arrays + scalar metrics, with no
        matplotlib — the headless path tests and scans reuse without paying the
        figure-build cost. `observables` wraps this and draws the figure."""
        x, T_trans, xphys = self._transmission(raw["chi_bar"], raw, params)
        center_mhz = float(raw.get("two_photon_resonance_scan_mhz", self._resonance_scan_mhz(params)))
        ic = int(np.argmin(np.abs(x - center_mhz)))
        width, _ = self._eit_features(x, T_trans, center_mhz=center_mhz)
        slope = (np.nanmax(np.abs(np.gradient(T_trans, x)))
                 if np.ptp(T_trans) > 0.0 else 0.0)
        if_delta = max(float(params.get("if_khz", self._REF["if_khz"])) / 1000.0,
                       1e-9)
        if_valid = (x >= x[0] + if_delta) & (x <= x[-1] - if_delta)
        if if_valid.any():
            t_hi = np.interp(x[if_valid] + if_delta, x, T_trans)
            t_lo = np.interp(x[if_valid] - if_delta, x, T_trans)
            if_readout = (t_hi - t_lo) / (2.0 * if_delta)
            finite_response = np.isfinite(if_readout)
            if finite_response.any():
                i_if = int(self._mirror_maxima(np.where(
                    finite_response, np.abs(if_readout), -np.inf))[0])
                if_disc = float(abs(if_readout[i_if]))
                if_detuning = (float(x[if_valid][i_if]) if if_disc > 0.0
                               else float("nan"))
            else:
                if_disc = float("nan")
                if_detuning = float("nan")
        else:
            if_disc = float("nan")
            if_detuning = float("nan")

        metrics = []
        at_resolved = None
        at_split_mhz = None
        at_field_v_m = None
        at_split_field_estimate_v_m = None
        at_peak_positions_mhz = []
        superhet = None
        sam = self._sam_calibration(params)
        if raw["lo_rabi_mhz"] > 0:
            dmw = abs(float(params.get("mw_detuning_mhz", self._REF["mw_detuning_mhz"])))
            rf_coupling = self._rf_coupling(params)
            at_field_v_m = rf_coupling.field_from_cyclic_rabi_hz(
                raw["lo_rabi_mhz"] * 1.0e6)
            convention = str(params.get(
                "rf_field_convention", self._REF["rf_field_convention"]))
            metrics.append(dict(
                label="RF field from fitted Ω_RF",
                value=f"{at_field_v_m * 1.0e3:.3f} mV/m {convention}",
                help="Absolute field converted from the RF Rabi parameter used in "
                     "the full optical-spectrum model. For experimental data, fit "
                     "Ω_RF to the spectrum before using this value."))
            window = max(8.0, 2.0 * raw["lo_rabi_mhz"] + dmw)
            peaks = self._transparency_maxima(x, T_trans, window, center_mhz=center_mhz)
            at_resolved = len(peaks) >= 2
            if at_resolved:
                xs = sorted(p[0] for p in sorted(peaks, key=lambda p: p[1])[-2:])
                at_peak_positions_mhz = list(xs)
                at_split_mhz = xs[1] - xs[0]
                metrics.append(dict(label="RF AT splitting", value=f"{at_split_mhz:.2f} MHz",
                                    help="Separation of the two tallest dressed peaks."))
                at_split_field_estimate_v_m = (
                    rf_coupling.field_from_at_splitting_hz(
                        at_split_mhz * 1.0e6,
                        detuning_hz=dmw * 1.0e6)
                    if at_split_mhz >= dmw else None)
                metrics.append(dict(
                    label="AT split field estimate",
                    value=(
                        f"{at_split_field_estimate_v_m * 1.0e3:.3f} mV/m {convention}"
                        if at_split_field_estimate_v_m is not None else "—"),
                    help="Ideal two-level sqrt(split²-detuning²) estimate. Optical "
                         "line pulling and multilevel effects are not corrected; the "
                         "full-spectrum fitted Ω_RF field is preferred."))
            if peaks:
                # Height-weighted centre of the transparency peaks: the symmetric
                # doublet sits near the bare two-photon resonance; a detuned microwave pulls it toward the stronger
                # dressed state, so the shift is readable even when one peak fades.
                wsum = sum(h for _, h in peaks)
                center = sum(px * h for px, h in peaks) / wsum - center_mhz
                metrics.append(dict(label="AT center shift", value=f"{center:+.2f} MHz",
                                    help="Height-weighted centre of the dressed "
                                         "transparency peaks relative to the bare "
                                         "two-photon resonance. Optical line pulling "
                                         "and RF detuning can both shift it."))
            if not at_resolved:
                metrics.append(dict(
                    label="AT status", value="doublet unresolved", kind="status",
                    help="The microwave dressing is on, but the optical spectrum "
                         "does not contain two resolved transparency peaks."))
        else:
            width_valid = np.isfinite(width) and width > 0
            metrics.append(dict(
                label="EIT linewidth",
                value=(f"{width:.2f} MHz" if width_valid else "—"),
                help=("FWHM of the central transparency feature." if width_valid else
                      "Withheld because no finite half-height crossings were found.")))
            if not width_valid:
                metrics.append(dict(
                    label="EIT status", value="window unresolved", kind="status",
                    help="The current opacity/coupling settings do not produce a "
                         "measurable EIT-window FWHM in this scan."))
        metrics.extend([
            dict(label="Max spectral slope", value=f"{slope:.3f} /MHz",
                 help="Largest static dT/dnu; used internally for electrometry tests."),
            dict(label="IF discriminator", value=f"{if_disc:.3f} /MHz",
                 help="Finite-difference transmission discriminator at the selected "
                      "IF offset; a static proxy for a lock-in/superhet readout."),
            dict(label="IF optimum detuning",
                 value=(f"{if_detuning:+.2f} MHz" if np.isfinite(if_detuning) else "—"),
                 help="Selected optical scan coordinate where the finite-IF discriminator is largest."),
            dict(label="Transmission at resonance", value=f"{T_trans[ic]:.3f}"),
        ])
        if sam is not None:
            metrics.append(dict(
                label="SAM RF field",
                value=(f"{sam.field_v_m * 1.0e3:.3f} ± "
                       f"{sam.standard_uncertainty_v_m * 1.0e3:.3f} mV/m "
                       f"{sam.amplitude_convention.upper()}"),
                help="Standard-antenna-method far-field calibration with the "
                     "declared one-standard-deviation uncertainty budget."))
            if sam.far_field_ratio is not None:
                metrics.append(dict(
                    label="SAM far-field ratio",
                    value=f"{sam.far_field_ratio:.3f}",
                    kind="status" if sam.warning else "secondary",
                    help=(sam.warning or
                          "Ratio r/(2D²/lambda); values >= 1 satisfy the declared "
                          "far-field distance criterion.")))
            if at_field_v_m is not None and sam.field_v_m > 0.0:
                metrics.append(dict(
                    label="Fitted-RF/SAM field ratio",
                    value=f"{at_field_v_m / sam.field_v_m:.3f}",
                    help="Optically inferred AT field divided by the independent "
                         "standard-antenna-method field."))
        if raw["lo_rabi_mhz"] > 0:
            superhet = self._superheterodyne_readout(
                raw, params, x, T_trans, atomic_response=atomic_response)
        if superhet is not None:
            sensitivity = superhet["sensitivity"]
            noise = superhet["noise"]
            metrics.extend([
                dict(
                    label="Total sensitivity",
                    value=(f"{sensitivity.total_field_asd_nv_cm_per_sqrt_hz:.2f} "
                           "nV/cm/√Hz"),
                    help="Finite-IF weak-SIG OBE response divided into the total "
                         "balanced-detector ASD (PSN + RIN + electronics). This is "
                         "a conditional model prediction, not the paper's 12.5 value "
                         "injected as an anchor."),
                dict(
                    label="PSN-limited sensitivity",
                    value=(f"{sensitivity.psn_field_asd_nv_cm_per_sqrt_hz:.2f} "
                           "nV/cm/√Hz"),
                    help="Probe photon-shot-noise contribution only, using the "
                         "declared detector QE, optical path, reference arm, RF "
                         "dipole, and RMS/peak convention."),
                dict(
                    label="Technical field noise",
                    value=(f"{sensitivity.technical_field_asd_nv_cm_per_sqrt_hz:.2f} "
                           "nV/cm/√Hz"),
                    help="Quadrature sum of residual RIN and detector electronics, "
                         "referred back to the RF field."),
                dict(
                    label="Superhet optimum detuning",
                    value=f"{superhet['optimum_detuning_mhz']:+.2f} MHz",
                    help="Selected optical scan coordinate minimizing total noise-equivalent field, "
                         "including transmission-dependent shot noise "
                         f"({superhet['optimization_scope']})."),
                dict(
                    label="Superhet response phase",
                    value=f"{superhet['response_phase_deg']:+.1f}°",
                    help="Photocurrent phasor phase relative to exp(-i omega_IF t)."),
                dict(
                    label="RMS noise in ENBW",
                    value=f"{superhet['rms_field_noise_v_m'] * 1.0e7:.2f} nV/cm",
                    help=f"Total white-noise ASD integrated over {superhet['enbw_hz']:.3g} Hz ENBW."),
                dict(
                    label="Detector shot-current ASD",
                    value=(f"{noise.shot_noise_current_asd_a_per_sqrt_hz * 1.0e12:.3f} "
                           "pA/√Hz"),
                    help="One-sided Schottky ASD of the signal and weighted "
                         "reference photodiodes, including dark current if supplied."),
            ])
            mw_detuning = abs(float(params.get(
                "mw_detuning_mhz", self._REF["mw_detuning_mhz"])))
            # A centre shift is only informative once the microwave detuning
            # exceeds both a visible control step and the finite-IF probe span.
            materially_detuned = mw_detuning >= max(0.1, 2.0 * if_delta)
            if at_resolved:
                preferred_heroes = (
                    ("Total sensitivity", "AT center shift", "RF field from fitted Ω_RF",
                     "Transmission at resonance")
                    if materially_detuned else
                    ("Total sensitivity", "PSN-limited sensitivity",
                     "RF field from fitted Ω_RF", "RF AT splitting")
                )
            else:
                preferred_heroes = (
                    "Total sensitivity", "AT status", "PSN-limited sensitivity")
        elif raw["lo_rabi_mhz"] > 0:
            metrics.append(dict(
                label="Sensitivity status", value="no finite optical response",
                kind="status", help="No field sensitivity is available at zero "
                     "effective opacity or without a finite detector response."))
            preferred_heroes = (
                "Sensitivity status", "RF field from fitted Ω_RF", "AT status")
        else:
            preferred_heroes = (
                ("EIT linewidth", "Transmission at resonance", "IF discriminator")
                if np.isfinite(width) and width > 0 else
                ("Transmission at resonance", "EIT status", "IF discriminator")
            )
        by_label = {metric["label"]: metric for metric in metrics}
        hero_labels = [
            label for label in preferred_heroes if label in by_label
        ][:2]
        for label in hero_labels:
            by_label[label]["tier"] = "hero"
        # `partition_metrics` preserves source order, so put the chosen heroes in
        # the scientific priority order above rather than merely tagging them.
        selected = set(hero_labels)
        metrics = ([by_label[label] for label in hero_labels]
                   + [metric for metric in metrics
                      if metric["label"] not in selected])
        metrics.append(dict(
            label="Prediction scope", value="PREDICTED · opacity uncalibrated",
            kind="status", help="The opacity and its scenario range are assumed; "
                 "absolute sensitivity and sweep optima lack held-out validation."))
        opacity_spec = opacity.opacity_spec(params)
        sensitivity_readout = dict(
            status=("NOT_APPLICABLE_NO_LO" if raw["lo_rabi_mhz"] <= 0 else
                    "UNAVAILABLE_ZERO_OPACITY" if opacity_spec["value"] == 0 else
                    "UNAVAILABLE_NO_FINITE_RESPONSE"),
            total_field_asd_nv_cm_per_sqrt_hz=None,
            psn_field_asd_nv_cm_per_sqrt_hz=None,
            optimum_detuning_mhz=None,
            optimum_scan_detuning_mhz=None,
            optimum_probe_detuning_mhz=None,
            optimum_coupling_detuning_mhz=None,
            optimum_two_photon_detuning_mhz=None,
            scan_axis=self._scan_axis(params),
            optimization_scope=None,
        )
        if superhet is not None:
            sensitivity_readout.update(
                status="PREDICTED",
                total_field_asd_nv_cm_per_sqrt_hz=float(
                    superhet["sensitivity"].total_field_asd_nv_cm_per_sqrt_hz),
                psn_field_asd_nv_cm_per_sqrt_hz=float(
                    superhet["sensitivity"].psn_field_asd_nv_cm_per_sqrt_hz),
                optimum_detuning_mhz=superhet["optimum_detuning_mhz"],
                optimum_scan_detuning_mhz=superhet["optimum_scan_detuning_mhz"],
                optimum_probe_detuning_mhz=superhet["optimum_probe_detuning_mhz"],
                optimum_coupling_detuning_mhz=superhet["optimum_coupling_detuning_mhz"],
                optimum_two_photon_detuning_mhz=superhet["optimum_two_photon_detuning_mhz"],
                optimization_scope=superhet["optimization_scope"],
            )
        trust_readout = dict(
            prediction_status="PREDICTED",
            calibration_status="UNCALIBRATED",
            experimental_visibility="NOT_ASSESSED",
            numerical_convergence="NOT_ASSESSED",
            effective_opacity=opacity_spec,
            atomic_line_strength=float(raw.get("atomic_line_strength", 1.0)),
            atom_number_participation_fraction=float(params.get(
                "atom_participation_fraction", self._REF["atom_participation_fraction"])),
            atom_number_scope="reporting only; independent of optical opacity",
            model_scope=("scalar four-level cascade; anchored optical Rabi/transit "
                         "model; pointwise balanced detector; no held-out "
                         "experimental validation"),
            doppler_mode=str(params.get("doppler", self._REF["doppler"])),
            scan_axis=self._scan_axis(params),
            fixed_probe_detuning_mhz=float(params.get("probe_detuning_mhz", self._REF["probe_detuning_mhz"])),
            fixed_coupling_detuning_mhz=float(params.get("coupling_detuning_mhz", self._REF["coupling_detuning_mhz"])),
            two_photon_resonance_scan_mhz=center_mhz,
            dephasing_model=raw.get("dephasing_model", "correlated_phase_lindblad"),
            eit_width_mhz=(float(width) if raw["lo_rabi_mhz"] <= 0
                           and np.isfinite(width) else None),
            at_split_mhz=at_split_mhz,
            max_slope_per_mhz=float(slope),
            if_discriminator_per_mhz=(float(if_disc) if np.isfinite(if_disc) else None),
            if_optimum_detuning_mhz=(float(if_detuning)
                                    if np.isfinite(if_detuning) else None),
            resonance_transmission=float(T_trans[ic]),
            finite_if=sensitivity_readout,
        )
        return dict(
            x=x, T_trans=T_trans, xphys=xphys, width=width, metrics=metrics,
            resonance_scan_mhz=center_mhz,
            max_slope_per_mhz=float(slope), if_discriminator_per_mhz=float(if_disc),
            rydberg_readout=trust_readout, effective_opacity=opacity_spec,
            prediction_status="PREDICTED",
            superhet=superhet, effective_atoms=self._effective_atom_number(raw, params),
            at_split_mhz=at_split_mhz, at_field_v_m=at_field_v_m,
            at_peak_positions_mhz=at_peak_positions_mhz,
            at_split_field_estimate_v_m=at_split_field_estimate_v_m, sam=sam)

    def observables(self, raw, params, include_figures=True):
        ro = self._readout(raw, params)
        x, T_trans, xphys, width = ro["x"], ro["T_trans"], ro["xphys"], ro["width"]
        superhet = ro["superhet"]
        effective_atoms = ro["effective_atoms"]
        sam = ro["sam"]
        view = params.get("view", "EIT")
        center_mhz = ro["resonance_scan_mhz"]

        fig = None
        if include_figures:
            import matplotlib.pyplot as plt

            if view == "EIT":
                # Ju et al. Fig. 2(a): EIT transmission with / without B-field
                # compensation, zoomed onto the transparency feature.
                _, T_uncomp, _ = self._transmission(
                    raw["chi_bar_uncomp"], raw, params)
                fig, axT = plt.subplots(figsize=(8.5, 4.8))
                axT.plot(x, T_uncomp, color="#1f77b4", lw=1.6,
                         label="without B-field compensation")
                axT.plot(x, T_trans, color="#d62728", lw=2.0,
                         label="with compensation")
                axT.set_xlabel(self._scan_label(params))
                axT.set_ylabel("Transmission")
                axT.axvline(center_mhz, color="gray", ls=":", lw=0.8)
                axT.legend(fontsize=9, loc="upper right")
                xlim = max(2.5, 2.5 * width)
                axT.set_xlim(center_mhz - xlim, center_mhz + xlim)
            else:
                fig, (axT, axR, axD) = plt.subplots(
                    3, 1, figsize=(8.5, 8.0), sharex=True)
                axT.plot(x, T_trans, color="#0f766e", lw=1.8)
                axT.set_ylabel("Transmission")
                axT.set_ylim(-0.02, 1.04)
                response_x = (superhet["response_detuning_mhz"]
                              if superhet is not None else np.array([]))
                response_per_khz = (
                    np.abs(superhet["transmission_phasor_per_angular_rabi"])
                    * 2.0 * np.pi * 1.0e3
                    if superhet is not None else np.array([]))
                if params.get("doppler", self._REF["doppler"]) == "on":
                    axR.plot(response_x, response_per_khz, "o",
                             color="#dc2626")
                else:
                    axR.plot(response_x, response_per_khz,
                             color="#dc2626", lw=1.6)
                if superhet is not None:
                    axR.axvline(superhet["optimum_detuning_mhz"], color="gray",
                                ls=":", lw=0.8)
                else:
                    axR.text(0.5, 0.5, "No finite optical response",
                             transform=axR.transAxes, ha="center")
                axR.set_ylabel("|dT| per 1 kHz SIG")
                axR.set_title(
                    f"Finite-IF linear response ({float(params.get('if_khz', self._REF['if_khz'])):.0f} kHz IF)")
                axD.plot(x, np.real(xphys), color="#7c3aed", lw=1.5)
                axD.set_ylabel("Re chi")
                axD.set_xlabel(self._scan_label(params))
                for a in (axT, axR, axD):
                    a.axvline(center_mhz, color="gray", ls=":", lw=0.8)
            fig.tight_layout()

        derived = derived_table([
            ("Ladder", "85Rb 5S1/2 F=3 → 5P3/2 F'=4 → 40D5/2"),
            ("RF leg", "40D5/2 → 39F7/2"),
            ("Microwave frequency",
             f"{float(params.get('mw_frequency_ghz', self._REF['mw_frequency_ghz'])):.1f} GHz"),
            ("Coupling Rabi Ω_c / 2π", f"{raw['coupling_rabi_mhz']:.3f} MHz"),
            ("Probe Rabi Ω_P / 2π", f"{raw['probe_rabi_mhz']:.3f} MHz"),
            ("Probe Rabi anchor / 2π", f"{raw['probe_rabi_ref_mhz']:.3f} MHz at 6 µW, 0.15 mm"),
            ("Optical scan axis", self._scan_axis(params)),
            ("Fixed probe detuning", f"{raw['fixed_probe_detuning_mhz']:+.3f} MHz (used for coupling scan)"),
            ("Fixed coupling detuning", f"{raw['fixed_coupling_detuning_mhz']:+.3f} MHz (used for probe scan)"),
            ("Bare two-photon resonance scan coordinate", f"{center_mhz:+.3f} MHz"),
            ("Dephasing", "CP nested diagonal optical/RF projectors"),
            ("Beam diameter", f"{raw['beam_diameter_mm']:.3f} mm"),
            ("Probe power", f"{raw['probe_power_uw']:.2f} µW"),
            ("Temperature definition", raw["temperature_model"]),
            ("Heater setpoint", f"{raw['heater_setpoint_c']:.2f} °C"),
            ("Effective vapor temperature", f"{raw['effective_temp_c']:.2f} °C"),
            ("Cold-spot temperature", f"{raw['cold_spot_temp_c']:.2f} °C"),
            ("Transit broadening / 2π", f"{raw['transit_mhz']:.3f} MHz"),
            ("Residual Zeeman (uncomp.) / 2π", f"{raw['residual_zeeman_mhz']:.3f} MHz"),
            ("Intrinsic 5S–40D dephasing / 2π", f"{raw['rydberg_dephasing_mhz']:.3f} MHz"),
            ("Temperature dephasing / 2π", f"{raw['temperature_dephasing_mhz']:.3f} MHz"),
            ("Density dephasing / 2π", f"{raw['density_dephasing_mhz']:.3f} MHz"),
            ("N(85Rb)", f"{raw['N']:.3e} /m³"),
            ("Effective Gaussian-overlap atoms", f"{effective_atoms.atoms:.3e}"),
            ("Effective overlap area", f"{effective_atoms.effective_area_m2:.3e} m²"),
            ("Prediction status", "PREDICTED; opacity UNCALIBRATED"),
            ("Atomic line strength", f"{float(raw.get('atomic_line_strength', 1.0)):g}"),
            ("Effective opacity participation", f"{ro['effective_opacity']['value']:.6g}"),
            ("Opacity source", ro["effective_opacity"]["source"]),
            ("Opacity scenario half-width",
             f"{ro['effective_opacity']['relative_uncertainty']:.0%} assumed; "
             "not a standard uncertainty or confidence interval"),
            ("Atom-number participation", "Reporting only; separate from opacity"),
        ])
        tables = [derived]
        if sam is not None:
            tables.append(derived_table([
                ("SAM field", f"{sam.field_v_m:.6e} V/m {sam.amplitude_convention.upper()}"),
                ("SAM standard uncertainty", f"{sam.standard_uncertainty_v_m:.6e} V/m"),
                ("SAM relative standard uncertainty", f"{sam.relative_standard_uncertainty:.3%}"),
                ("Power delivered to antenna", f"{sam.power_at_antenna_w:.6e} W"),
                ("Net antenna gain", f"{sam.net_gain_dbi:.3f} dBi"),
                ("Far-field ratio r/(2D²/λ)", (
                    "not evaluated" if sam.far_field_ratio is None
                    else f"{sam.far_field_ratio:.4f}")),
                ("SAM warning", sam.warning or "none"),
            ]))
        if superhet is not None:
            sensitivity = superhet["sensitivity"]
            noise = superhet["noise"]
            rf_coupling = superhet["rf_coupling"]
            convention = rf_coupling.field_amplitude_convention.upper()
            tables.append(derived_table([
                ("Response model", "first-order finite-IF Liouvillian response"),
                ("Optimization scope", superhet["optimization_scope"]),
                ("RF dipole", f"{float(params.get('rf_transition_dipole_ea0', self._REF['rf_transition_dipole_ea0'])):.3f} e a₀"),
                ("RF angular/polarization factor", f"{rf_coupling.angular_factor:.4f}"),
                ("Field convention", convention),
                ("Optimum scan detuning", f"{superhet['optimum_scan_detuning_mhz']:+.4f} MHz"),
                ("Optimum physical probe detuning", f"{superhet['optimum_probe_detuning_mhz']:+.4f} MHz"),
                ("Optimum physical coupling detuning", f"{superhet['optimum_coupling_detuning_mhz']:+.4f} MHz"),
                ("Detector responsivity", f"{superhet['detector_responsivity_a_per_w']:.6f} A/W"),
                ("Signal detector power", f"{superhet['signal_detector_power_w'] * 1e6:.4f} µW"),
                ("Reference detector power", f"{superhet['reference_detector_power_w'] * 1e6:.4f} µW"),
                ("Current responsivity", f"{superhet['current_responsivity_a_per_v_m']:.6e} A/(V/m)"),
                ("PSN current ASD", f"{noise.shot_noise_current_asd_a_per_sqrt_hz:.6e} A/√Hz"),
                ("RIN current ASD", f"{noise.rin_noise_current_asd_a_per_sqrt_hz:.6e} A/√Hz"),
                ("Electronics current ASD", f"{noise.electronic_noise_current_asd_a_per_sqrt_hz:.6e} A/√Hz"),
                ("PSN field sensitivity", f"{sensitivity.psn_field_asd_nv_cm_per_sqrt_hz:.4f} nV/cm/√Hz"),
                ("Total field sensitivity", f"{sensitivity.total_field_asd_nv_cm_per_sqrt_hz:.4f} nV/cm/√Hz"),
            ]))
        comparison = dict(
            axis_index=0,
            x_unit="MHz",
            raw_x_unit="instrument unit",
            raw_y_unit="detector unit",
            label="Rydberg experimental CSV",
        )
        return dict(
            metrics=ro["metrics"], figure=fig, tables=tables,
            rydberg_readout=ro["rydberg_readout"],
            rydberg_spectrum=dict(x=ro["x"], transmission=ro["T_trans"],
                                  dispersion=np.real(ro["xphys"]),
                                  at_peak_positions_mhz=ro["at_peak_positions_mhz"],
                                  finite_if=ro["superhet"]),
            effective_opacity=ro["effective_opacity"],
            prediction_status="PREDICTED",
            comparison=comparison)

    def extra_views(self):
        """Ju et al. Fig. 2(b): EIT peak amplitude and linewidth vs probe power,
        with / without B-field compensation, plus the cell-heating sweep used
        to diagnose the temperature-dependent RF operating point."""
        def _compute_sweep(params):
            powers = np.linspace(0.5, 10.0, 16)
            spec = opacity.opacity_spec(params)
            scenarios = [dict(sample, comp={"width": [], "amp": []},
                              uncomp={"width": [], "amp": []},
                              prediction_status="PREDICTED")
                         for sample in opacity.opacity_samples(params)]
            for p_uw in powers:
                pr = dict(params, probe_power_uw=float(p_uw), view="EIT")
                raw = self.compute(pr)
                for scenario in scenarios:
                    scenario_params = dict(
                        pr, effective_opacity_participation=scenario["value"])
                    for key, chi_key in (("comp", "chi_bar"),
                                         ("uncomp", "chi_bar_uncomp")):
                        x, transmission, _ = self._transmission(
                            raw[chi_key], raw, scenario_params)
                        width, contrast = self._eit_features(
                            x, transmission, center_mhz=self._resonance_scan_mhz(pr))
                        scenario[key]["width"].append(float(width))
                        scenario[key]["amp"].append(float(p_uw * contrast))
            for scenario in scenarios:
                scenario["optima"] = {}
                for key in ("comp", "uncomp"):
                    amplitudes = np.asarray(scenario[key]["amp"])
                    widths = np.asarray(scenario[key]["width"])
                    valid = np.isfinite(widths)
                    scenario["optima"][key] = dict(
                        peak_amplitude_power_uw=(
                            float(powers[int(np.argmax(amplitudes))])
                            if np.any(amplitudes > 0) else None),
                        minimum_linewidth_power_uw=(
                            float(powers[int(np.argmin(
                                np.where(valid, widths, np.inf)))])
                            if valid.any() else None),
                    )
            nominal = next(s for s in scenarios if s["value"] == spec["value"])
            return dict(
                powers=powers.tolist(), comp=nominal["comp"],
                scan_axis=self._scan_axis(params),
                two_photon_resonance_scan_mhz=self._resonance_scan_mhz(params),
                uncomp=nominal["uncomp"], opacity_scenarios=scenarios,
                effective_opacity=spec, prediction_status="PREDICTED",
                calibration_status="UNCALIBRATED",
                uncertainty_scope="assumed opacity scenarios; not a confidence interval",
                optimum_scope="sampled model optimum; not an experimental recommendation",
            )

        def _render_sweep(s):
            import matplotlib.pyplot as plt
            powers = np.array(s["powers"])
            figF, (aA, aW) = plt.subplots(2, 1, figsize=(8.0, 6.4), sharex=True)
            for ax in (aA, aW):
                ax.grid(alpha=0.3)
            aA.plot(powers, s["uncomp"]["amp"], "s-", color="#1f77b4",
                    label="without compensation")
            aA.plot(powers, s["comp"]["amp"], "o-", color="#d62728",
                    label="with compensation")
            aA.set_ylabel("EIT peak amplitude [arb.]")
            aA.legend(fontsize=9)
            aW.plot(powers, s["uncomp"]["width"], "s-", color="#1f77b4")
            aW.plot(powers, s["comp"]["width"], "o-", color="#d62728")
            labelled_incomplete = set()
            for key, color in (("comp", "#d62728"), ("uncomp", "#1f77b4")):
                for ax, field in ((aA, "amp"), (aW, "width")):
                    curves = np.asarray([v[key][field] for v in s["opacity_scenarios"]])
                    lower, upper, counts = _finite_scenario_envelope(curves)
                    ax.fill_between(powers, lower, upper, color=color, alpha=0.12)
                    if np.any(counts < curves.shape[0]) and ax not in labelled_incomplete:
                        ax.text(0.02, 0.03,
                                "Unavailable scenarios excluded; band covers finite cases only",
                                transform=ax.transAxes, fontsize=7)
                        labelled_incomplete.add(ax)
            aW.set_ylabel("EIT linewidth [MHz]")
            aW.set_xlabel("Probe power [uW]")
            figF.suptitle("PREDICTED: assumed opacity scenarios (not a confidence interval)",
                          fontsize=10)
            figF.tight_layout()
            return figF

        def _compute_temperature_sweep(params):
            # Include the exact nominal values used in the cell-heating deck.
            # This generic panel deliberately treats them as linked simulation
            # temperatures.  The analysis CLI accepts a separate heater ->
            # effective/cold-spot calibration for experimental overlays.
            temperatures = np.array(
                [20.0, 30.0, 35.0, 40.0, 42.0, 45.0, 47.0, 50.0, 53.0, 55.0])
            spec = opacity.opacity_spec(params)
            scenario_fields = (
                "spectra", "widths_mhz", "contrasts", "resonance_transmission",
                "max_slope_per_mhz", "if_discriminator_per_mhz",
                "psn_sensitivity_nv_cm_sqrt_hz", "total_sensitivity_nv_cm_sqrt_hz",
                "superhet_optimum_detuning_mhz", "finite_if_status",
            )
            doppler_on = str(params.get("doppler", self._REF["doppler"])) == "on"
            scenario_scope = (
                "fixed nominal static-gradient sampled point at each temperature"
                if doppler_on else "full optical-detuning scan at each temperature")
            scenarios = [dict(
                sample, **{key: [] for key in scenario_fields},
                prediction_status="PREDICTED", optimization_scope=scenario_scope,
            ) for sample in opacity.opacity_samples(params)]
            densities = []
            transit = []
            effective_atoms = []
            for temp_c in temperatures:
                pr = dict(
                    params, temperature_model="Linked", temp_c=float(temp_c),
                    view="EIT", lo_rabi_mhz=0.0)
                raw = self.compute(pr)
                densities.append(float(raw["N"]))
                transit.append(float(raw["transit_mhz"]))
                effective_atoms.append(float(self._effective_atom_number(raw, pr).atoms))

                at_lo = max(float(params.get(
                    "lo_rabi_mhz", self._REF["lo_rabi_mhz"])),
                    self._REF["lo_rabi_mhz"])
                at_params = dict(
                    params, temperature_model="Linked", temp_c=float(temp_c),
                    view="AT electrometry", lo_rabi_mhz=at_lo)
                at_raw = self.compute(at_params)
                atomic_response = None
                for scenario in scenarios:
                    scenario_params = dict(
                        pr, effective_opacity_participation=scenario["value"])
                    ro = self._readout(raw, scenario_params)
                    x, transmission = ro["x"], ro["T_trans"]
                    center_mhz = self._resonance_scan_mhz(pr)
                    width, contrast = self._eit_features(x, transmission, center_mhz=center_mhz)
                    ic = int(np.argmin(np.abs(x - center_mhz)))
                    scenario["spectra"].append(transmission.tolist())
                    scenario["widths_mhz"].append(float(width))
                    scenario["contrasts"].append(float(contrast))
                    scenario["resonance_transmission"].append(float(transmission[ic]))
                    # Preserve the legacy three-decimal sweep values while
                    # keeping numeric readout fields independent of UI strings.
                    scenario["max_slope_per_mhz"].append(
                        float(f"{ro['max_slope_per_mhz']:.3f}"))
                    scenario["if_discriminator_per_mhz"].append(
                        float(f"{ro['if_discriminator_per_mhz']:.3f}"))
                    at_scenario_params = dict(
                        at_params, effective_opacity_participation=scenario["value"])
                    at_ro = self._readout(
                        at_raw, at_scenario_params, atomic_response=atomic_response)
                    if at_ro["superhet"] is not None and atomic_response is None:
                        atomic_response = at_ro["superhet"]["atomic_response"]
                    sensitivity = at_ro["rydberg_readout"]["finite_if"]
                    for field, source in (
                        ("psn_sensitivity_nv_cm_sqrt_hz", "psn_field_asd_nv_cm_per_sqrt_hz"),
                        ("total_sensitivity_nv_cm_sqrt_hz", "total_field_asd_nv_cm_per_sqrt_hz"),
                    ):
                        value = sensitivity[source]
                        scenario[field].append(float("nan") if value is None else value)
                    scenario["superhet_optimum_detuning_mhz"].append(
                        sensitivity["optimum_detuning_mhz"])
                    scenario["finite_if_status"].append(sensitivity["status"])
            for scenario in scenarios:
                values = np.asarray(scenario["total_sensitivity_nv_cm_sqrt_hz"])
                valid = np.isfinite(values)
                best_index = int(np.argmin(np.where(valid, values, np.inf)))
                scenario["best_sensitivity_temperature_c"] = (
                    float(temperatures[best_index]) if valid.any() else None)
                scenario["best_total_sensitivity_nv_cm_sqrt_hz"] = (
                    float(values[best_index]) if valid.any() else None)
            nominal = next(s for s in scenarios if s["value"] == spec["value"])
            return dict(
                temperatures_c=temperatures.tolist(),
                temperature_semantics="linked effective/cold-spot simulation temperature",
                frequency_mhz=x.tolist(),
                scan_axis=self._scan_axis(params),
                two_photon_resonance_scan_mhz=self._resonance_scan_mhz(params),
                density_m3=densities,
                transit_mhz=transit,
                effective_atoms=effective_atoms,
                **{key: nominal[key] for key in scenario_fields},
                best_sensitivity_temperature_c=nominal["best_sensitivity_temperature_c"],
                best_total_sensitivity_nv_cm_sqrt_hz=nominal["best_total_sensitivity_nv_cm_sqrt_hz"],
                opacity_scenarios=scenarios, effective_opacity=spec,
                prediction_status="PREDICTED", calibration_status="UNCALIBRATED",
                optimization_scope=scenario_scope,
                uncertainty_scope="assumed opacity scenarios; not a confidence interval",
                optimum_scope="sampled model optimum; not an experimental recommendation",
            )

        def _render_temperature_sweep(sweep):
            import matplotlib.pyplot as plt

            temperatures = np.asarray(sweep["temperatures_c"], dtype=float)
            x = np.asarray(sweep["frequency_mhz"], dtype=float)
            spectra = np.asarray(sweep["spectra"], dtype=float)
            fig, axes = plt.subplots(2, 2, figsize=(10.0, 7.4))
            ax_spectrum, ax_feature, ax_sensitivity, ax_density = axes.ravel()
            colors = plt.cm.viridis(np.linspace(0.05, 0.95, temperatures.size))
            for temp_c, transmission, color in zip(
                    temperatures, spectra, colors):
                ax_spectrum.plot(x, transmission, color=color, lw=1.1,
                                 label=f"{temp_c:g} C")
            ax_spectrum.set_xlabel(sweep["scan_axis"] + " [MHz]")
            ax_spectrum.set_ylabel("Transmission")
            ax_spectrum.set_title("Temperature-dependent EIT spectra")
            ax_spectrum.legend(ncol=2, fontsize=7)

            ax_feature.plot(temperatures, sweep["contrasts"],
                            "s-", color="#0f766e", label="EIT contrast")
            ax_feature.set_xlabel("Cell temperature [C]")
            ax_feature.set_ylabel("EIT contrast", color="#0f766e")
            ax_width = ax_feature.twinx()
            ax_width.plot(temperatures, sweep["widths_mhz"],
                          "o--", color="#7c3aed", label="EIT FWHM")
            ax_width.set_ylabel("EIT FWHM [MHz]", color="#7c3aed")
            ax_feature.set_title("EIT contrast and linewidth")

            ax_sensitivity.semilogy(
                temperatures, sweep["psn_sensitivity_nv_cm_sqrt_hz"],
                "o--", label="PSN limit")
            ax_sensitivity.semilogy(
                temperatures, sweep["total_sensitivity_nv_cm_sqrt_hz"],
                "s-", label="total")
            scenarios = sweep["opacity_scenarios"]
            curves = np.asarray([v["total_sensitivity_nv_cm_sqrt_hz"]
                                 for v in scenarios])
            lower, upper, counts = _finite_scenario_envelope(curves)
            incomplete = np.any(counts < curves.shape[0])
            ax_sensitivity.fill_between(
                temperatures, lower, upper, alpha=0.15,
                label=("finite opacity scenarios only" if incomplete else
                       "assumed opacity scenarios"))
            if incomplete:
                ax_sensitivity.text(
                    0.02, 0.03,
                    "Unavailable scenario sensitivity: no finite range bound; "
                    "band covers finite cases only",
                    transform=ax_sensitivity.transAxes, fontsize=7, wrap=True)
            if sweep["best_sensitivity_temperature_c"] is not None:
                ax_sensitivity.axvline(
                    sweep["best_sensitivity_temperature_c"], color="gray",
                    ls=":", lw=1.0)
            ax_sensitivity.set_xlabel("Linked cell temperature [C]")
            ax_sensitivity.set_ylabel("Sensitivity [nV/cm/sqrt(Hz)]")
            ax_sensitivity.set_title("Finite-IF electrometry sensitivity")
            ax_sensitivity.legend(fontsize=8)

            ax_density.semilogy(
                temperatures, sweep["density_m3"], "o-", color="#2563eb",
                label="density")
            ax_density.axvline(39.30, color="gray", ls=":", lw=1.0,
                               label="Rb phase boundary")
            ax_density.set_xlabel("Cell temperature [C]")
            ax_density.set_ylabel("Rb density [m^-3]")
            ax_atoms = ax_density.twinx()
            ax_atoms.semilogy(
                temperatures, sweep["effective_atoms"], "s--",
                color="#dc2626", label="effective atoms")
            ax_atoms.set_ylabel("Effective atoms", color="#dc2626")
            ax_density.set_title("Vapor density and Gaussian-overlap atoms")
            ax_density.legend(fontsize=8)
            for ax in axes.ravel():
                ax.grid(alpha=0.25)
            fig.suptitle("PREDICTED: assumed opacity range (not a confidence interval)",
                         fontsize=10)
            fig.tight_layout()
            return fig

        return [
            ExtraView(
                key="Fig. 2(b): probe-power dependence (peak amplitude & linewidth)",
                description="Sweeps the probe power 0.5–10 µW and extracts the EIT "
                            "peak amplitude and linewidth with and without B-field "
                            "compensation. PREDICTED with assumed opacity scenarios; "
                            "the band is not a confidence interval or external validation.",
                compute=_compute_sweep, render=_render_sweep,
                param_keys=("cell_mm", "effective_opacity_participation",
                            "effective_opacity_relative_uncertainty"),
            ),
            ExtraView(
                key="Cell-heating sweep: EIT and finite-IF sensitivity",
                description="Sweeps the nominal temperatures used in the cell-heating "
                            "study and reports EIT spectra, linewidth/contrast, "
                            "Gaussian-overlap atom number, and the finite-IF PSN/total "
                            "sensitivity. This generic panel uses linked temperatures; "
                            "use the analysis CLI to apply measured heater-to-effective "
                            "and cold-spot calibrations. PREDICTED optima remain "
                            "conditional on the uncalibrated opacity scenario range.",
                compute=_compute_temperature_sweep,
                render=_render_temperature_sweep,
                param_keys=(
                    "cell_mm", "if_khz",
                    "rf_transition_dipole_ea0", "rf_angular_factor",
                    "rf_field_convention", "detector_quantum_efficiency",
                    "detector_path_efficiency", "detector_reference_power_ratio",
                    "detector_electronic_noise_pa_sqrt_hz",
                    "detector_rin_per_sqrt_hz", "detector_rin_correlation",
                    "atom_participation_fraction", "beam_overlap_efficiency",
                    "effective_opacity_participation",
                    "effective_opacity_relative_uncertainty",
                ),
            ),
        ]

    def info(self):
        return (
            "**Rydberg-EIT electrometry.** 85Rb cascade-EIT / microwave Autler-Townes "
            "model with an explicit probe-Rabi calibration anchor and the reported "
            "coupling-Rabi reference of Ju et al. The EIT view models transmission with / without B-field "
            "compensation; the probe-power panel models its power dependence. "
            "The effective optical opacity is an explicit uncalibrated input. "
            "Temperature/power sweeps and their sampled optima are PREDICTED; "
            "their opacity bands are assumed scenarios, not confidence intervals. "
            "Around the LO-dressed steady state, a finite-IF first-order "
            "Liouvillian solve predicts the weak-SIG transmission phasor. Combined "
            "with the declared RF dipole and balanced-detector noise budget, it "
            "reports PSN, technical, and total field sensitivity. This is not a "
            "full time-domain waveform/lock-in simulation.\n\n"
            "Probe (780 nm) and coupling (481 nm) powers drive Ω_P, Ω_c via √(P)/d "
            "intensity scaling (anchored), and the beam diameter also sets the "
            "transit-time broadening. The probe anchor and transit factor have not "
            "been jointly validated against the measured power-dependent line shapes. "
            "Doppler on averages the full residual counter-propagating velocity shifts; "
            "off is a single-class approximation. Completely positive nested diagonal "
            "collapse operators model optical/RF dephasing. Independent probe and "
            "coupling detunings allow either laser to be scanned at a fixed detuning "
            "of the other; the bare two-photon resonance is Δp+Δc=0.\n\n"
            "**Calibration scope.** The default 40D5/2→39F7/2 dipole is the ARC "
            "3.10.2 stretched-state σ+ value; the RF angular/polarization factor, "
            "detector throughput, reference-arm balance, RIN, and electronics noise "
            "remain explicit inputs. Therefore the displayed absolute sensitivity is "
            "conditional on those inputs and is not forced to the paper's 11.2 or "
            "12.5 nV/cm/√Hz reference numbers. The temperature model can also keep "
            "heater setpoint, effective vapor temperature, and cold spot separate. "
            "The CSV comparison panel overlays EIT/AT traces; the batch analysis "
            "workflow additionally accepts RF sweeps and PSD/ASD files with SHA-256 "
            "provenance.\n\n"
            "**References**\n"
            "- [arXiv:2606.04354](https://arxiv.org/abs/2606.04354) — \"Photon "
            "shot-noise-limited Rydberg-EIT electrometry\", Ju et al. (85Rb vapor "
            "cell, 50 mm cell, 6 µW probe, 30 mW coupling, 0.15 mm beam, 37 GHz RF).\n"
            "- [ARC Alkali Rydberg Calculator](https://arc-alkali-rydberg-calculator.readthedocs.io/) "
            "— RF matrix-element provenance for the default stretched-state value."
        )
