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

import numpy as np

from .. import atoms, beam, constants, core, doppler, kernels, observables, species
from .. import rydberg_electrometry as electrometry
from .. import rydberg_experiment as rydberg_experiment
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

# --- Calibration constants fitted to Ju et al. (arXiv:2606.04354) ---
# Transit-time broadening of the 5S-40D coherence: an atom crosses the beam in
# ~d/v_mp, so the coherence decays at ~v_mp/d. The O(1) factor is fit with the
# probe-power anchor so the reference (0.15 mm, 6 µW, 30 mW) lands on the paper's
# 1.6 MHz EIT linewidth (the transit-time term is the zero-probe floor, ~1.4 MHz).
TRANSIT_FACTOR = 0.6
# Weak-probe drive anchor: Ω_P/2π [MHz] at the reference 6 µW / 0.15 mm. A
# first-principles Ω_P from the bare D2 dipole overestimates by ~100x (optical
# pumping + Zeeman/hyperfine sub-structure are not in the lumped 2-level probe
# leg), so the probe drive is anchored like the coupling. Finite Ω_P power-
# broadens the EIT (Ju et al. Fig. 2b); fit so 6 µW sits at ~1.6 MHz. Because
# 6 µW is pinned to the narrow operating point, the broadening over 0–10 µW is
# milder than the paper's wide-range Fig. 2(b) (whose narrow point is at lower
# probe power).
PROBE_RABI_REF_MHZ = 2.0


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
    cache_version = "rydberg-eit-v6"
    defaults_version = "rydberg-eit-v6"
    supports_headless_observables = True

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

    def param_schema(self):
        r = self._REF
        return [
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
                      help="Beam 1/e² diameter. Sets Ω_c, Ω_P ∝ 1/d (intensity ∝ "
                           "P/d²) AND the transit-time broadening ∝ v_thermal/d — "
                           "a tighter beam is brighter but transit-broadens the EIT."),
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
            ParamSpec("lo_rabi_mhz", "Microwave Rabi Ω_RF", "Fields",
                      r["lo_rabi_mhz"], 0.0, 20.0, 0.1, "MHz",
                      help="Rydberg 40D -> 39F dressing Rabi frequency Ω_RF/2π."),
            ParamSpec("mw_detuning_mhz", "Microwave detuning", "Detunings",
                      r["mw_detuning_mhz"], -20.0, 20.0, 0.1, "MHz"),
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
                      help="off (default): the counter-propagating geometry suppresses "
                           "the two-photon Doppler — the transit-limited regime Ju et al. "
                           "operate in (1.6 MHz). on: Maxwell-average the residual "
                           "(k_probe − k_coupling)·v; in this lumped 4-level model that "
                           "over-broadens the EIT (~2.5 MHz floor), so it is a "
                           "what-if, not the calibrated reference."),
        ]

    def _defaults(self, view):
        vals = dict(self._REF)
        vals["view"] = view
        if view == "EIT":
            vals["lo_rabi_mhz"] = 0.0
        return vals

    def recommended_defaults(self, params):
        return {
            "EIT": self._defaults("EIT"),
            "AT electrometry": self._defaults("AT electrometry"),
        }

    def _atom(self, ground_deph, rf_deph):
        """4-level cascade with the constant topology from `_cascade_skeleton`;
        only the two dephasing channels (5S–40D ground–Rydberg and 40D–39F
        Rydberg–Rydberg) are injected per call."""
        sk = _cascade_skeleton()
        deph = (
            (0, 2, ground_deph), (2, 0, ground_deph),
            (0, 3, ground_deph), (3, 0, ground_deph),
            (2, 3, rf_deph), (3, 2, rf_deph),
        )
        return atoms.AtomModel(
            name="rydberg_eit_85rb",
            n_levels=4,
            labels=sk["labels"],
            ground=(0,),
            excited=(1, 2, 3),
            decay=sk["decay"],
            dephasing=deph,
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
            PROBE_RABI_REF_MHZ, p, p_ref, diameter=d, ref_diameter=d_ref)

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

    def _scan_chi(self, atom, h_of, scan, probe, kv, w):
        """χ̄(scan) for one dephasing configuration via the affine kernel
        (Doppler-off kv=[0]; Doppler-on Maxwell grid). Shared by the compensated /
        uncompensated solves and the Fig. 2(b) probe-power sweep."""
        n = atom.n_levels
        if kernels.available():
            base = core.build_liouvillian(h_of(0.0), atom)
            A_coef = core.build_liouvillian(h_of(1.0), atom) - base
            with core.blas_single_thread():
                return kernels.affine_scan_chi(
                    base, A_coef, atom.S_v,
                    np.ascontiguousarray(scan, dtype=float),
                    np.ascontiguousarray(kv, dtype=float),
                    np.ascontiguousarray(w, dtype=float), 1 * n + 0, n) / probe
        chi = np.zeros(scan.size, dtype=complex)
        for i, s in enumerate(scan):
            L0 = core.build_liouvillian(h_of(s), atom)
            # delta = -kv reproduces the kernel's +kv·S_v velocity shift.
            rho = core.steady_state_batched(L0, -kv, atom.S_v, n)
            chi[i] = ((rho[:, 1, 0] / probe) * w).sum()
        return chi

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
        return np.linspace(-half, half, 801) * MHZ

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
        Dmw = float(params.get("mw_detuning_mhz", self._REF["mw_detuning_mhz"])) * MHZ

        def h_of(s):
            H = np.zeros((4, 4), dtype=complex)
            H[1, 1] = -s
            H[2, 2] = -s
            H[3, 3] = -s - Dmw
            H[0, 1] = H[1, 0] = probe / 2
            H[1, 2] = H[2, 1] = Oc / 2
            H[2, 3] = H[3, 2] = Olo / 2
            return H

        scan = self._scan(params)
        # Doppler off (default): single static class (kv = 0) — the suppressed-
        # residual, transit-limited regime. Doppler on: Maxwell-average the
        # residual two-photon shift carried by atom.S_v (per-level k ratios).
        if params.get("doppler", "off") == "on":
            v, w = doppler.velocity_grid(
                T, mass=constants.MASS_85RB, dv=2.0, cutoff_sigma=4.0)
            kv = k_vec * v
        else:
            w, kv = np.ones(1), np.zeros(1)

        # Compensated EIT: one solve. Uncompensated: residual Zeeman scatters the
        # m-sublevel shifts, so the EIT line is INHOMOGENEOUSLY broadened — model
        # it as a Gaussian (FWHM = residual_zeeman) convolution of the compensated
        # susceptibility over the probe detuning. (A homogeneous dephasing knob
        # only lowers the peak in this Ω_c-limited regime; it does not widen it.)
        chi = self._scan_chi(self._atom(gd, rf_deph), h_of, scan, probe, kv, w)
        chi_uncomp = self._inhomogeneous(chi, scan, zeeman)

        return dict(
            scan=scan,
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
            k_vec=k_vec,
            omega0=2 * np.pi * nu0,
            dipole=dipole,
            gamma_mhz=line["gamma_mhz"],
            coupling_rabi_mhz=Oc / MHZ,
            probe_rabi_mhz=probe / MHZ,
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
    def _transparency_maxima(x, y, window_mhz):
        """Local transmission maxima (pos [MHz], height) within ±window of line
        centre — the dressed AT/EIT transparency peaks."""
        return [(float(x[i]), float(y[i])) for i in range(1, y.size - 1)
                if abs(x[i]) <= window_mhz and y[i] > y[i - 1] and y[i] >= y[i + 1]]

    def _transmission(self, chi, raw, params):
        """(x [MHz], Beer-Lambert transmission, χ_phys) for one coherence array."""
        x = raw["scan"] / MHZ
        alpha, xphys = observables.absorption_coefficient(
            chi, raw["k_vec"], raw["N"], dipole=raw["dipole"], line_strength=raw["ls"])
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

    @staticmethod
    def _hamiltonian_at_probe_detuning(raw, params, detuning_rad_s):
        """Rebuild the LO-dressed Hamiltonian at one probe detuning."""
        H = np.zeros((4, 4), dtype=complex)
        H[1, 1] = -detuning_rad_s
        H[2, 2] = -detuning_rad_s
        H[3, 3] = -detuning_rad_s - float(params.get(
            "mw_detuning_mhz", 0.0)) * MHZ
        H[0, 1] = H[1, 0] = raw["probe_rabi_mhz"] * MHZ / 2.0
        H[1, 2] = H[2, 1] = raw["coupling_rabi_mhz"] * MHZ / 2.0
        H[2, 3] = H[3, 2] = raw["lo_rabi_mhz"] * MHZ / 2.0
        return H

    def _superheterodyne_readout(self, raw, params, x, transmission):
        """Finite-IF weak-SIG response and an explicit detector-noise budget.

        The response is solved over the probe scan in the default Doppler-off
        path.  With residual Doppler averaging enabled, solving scan x velocity
        would create an unnecessarily large batch, so the complex velocity-class
        responses are coherently averaged at the static maximum-slope point and
        that limitation is exposed in ``optimization_scope``.
        """
        if raw["lo_rabi_mhz"] <= 0.0:
            return None

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
            static_index = int(np.nanargmax(np.abs(np.gradient(transmission, x))))
            indices = np.array([static_index], dtype=int)
            H = self._hamiltonian_at_probe_detuning(
                raw, params, raw["scan"][static_index])
            base = core.build_liouvillian(H, atom)
            velocity, weights = doppler.velocity_grid(
                raw["T"], mass=constants.MASS_85RB, dv=2.0, cutoff_sigma=4.0)
            kv = raw["k_vec"] * velocity
            liouvillian = base[None, :, :] + kv[:, None, None] * atom.S_v[None, :, :]
            response = electrometry.weak_signal_response_from_liouvillian(
                liouvillian, 4, omega_if, signal_transition=(2, 3))
            class_phasors = response.real_observable_phasor_per_angular_rabi(
                operator)
            absorption_phasor = electrometry.coherent_weighted_average(
                class_phasors, weights)
            absorption_phasors = np.array([absorption_phasor], dtype=complex)
            optimization_scope = "static-slope point; coherent Doppler average"
        else:
            indices = np.arange(raw["scan"].size, dtype=int)
            liouvillian = np.asarray([
                core.build_liouvillian(
                    self._hamiltonian_at_probe_detuning(raw, params, detuning), atom)
                for detuning in raw["scan"]
            ])
            response = electrometry.weak_signal_response_from_liouvillian(
                liouvillian, 4, omega_if, signal_transition=(2, 3))
            absorption_phasors = (
                response.real_observable_phasor_per_angular_rabi(operator))
            optimization_scope = "full probe-detuning scan"

        # chi_bar = rho_10 / Omega_probe.  Linearizing
        # T=exp[-k Im(chi_phys)L] gives the complex transmission phasor below.
        chi_scale = observables.chi_phys(
            1.0, raw["N"], dipole=raw["dipole"], line_strength=raw["ls"])
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
        if best is None:  # pragma: no cover - guarded by finite argmin above
            raise RuntimeError("finite-IF scan produced no calibrated response")
        current_phasor, signal_power_w, reference_channel, noise, sensitivity = best
        enbw_hz = float(params.get(
            "measurement_enbw_hz", self._REF["measurement_enbw_hz"]))
        return dict(
            response_detuning_mhz=x[indices],
            transmission_phasor_per_angular_rabi=transmission_phasor_per_rabi,
            optimum_index=best_index,
            optimum_detuning_mhz=float(x[best_index]),
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
        )

    @staticmethod
    def _eit_features(x, T_trans):
        """(FWHM [MHz], peak contrast) of the central EIT transparency feature.

        Contrast = T(centre) − absorptive floor in a ±5 MHz window. The FWHM is
        measured to the half-contrast level with linear interpolation of the two
        crossings (sub-sample accurate — a bare grid walk stair-steps the
        probe-power sweep). The transmitted-signal amplitude in Fig. 2(b) scales
        as probe power × contrast."""
        ic = int(np.argmin(np.abs(x)))
        peak = float(T_trans[ic])
        win = np.abs(x) <= 5.0
        floor = float(np.min(T_trans[win])) if win.any() else float(np.min(T_trans))
        contrast = peak - floor
        if contrast <= 0:
            return float("nan"), 0.0
        half = floor + 0.5 * contrast

        def crossing(step):
            i = ic
            while 0 < i < x.size - 1 and T_trans[i] >= half:
                i += step
            x0, y0, x1, y1 = x[i - step], T_trans[i - step], x[i], T_trans[i]
            return x1 if y1 == y0 else x0 + (half - y0) * (x1 - x0) / (y1 - y0)

        return abs(crossing(1) - crossing(-1)), contrast

    def _readout(self, raw, params):
        """Cheap transmission/dispersion arrays + scalar metrics, with no
        matplotlib — the headless path tests and scans reuse without paying the
        figure-build cost. `observables` wraps this and draws the figure."""
        x, T_trans, xphys = self._transmission(raw["chi_bar"], raw, params)
        ic = int(np.argmin(np.abs(x)))
        width, _ = self._eit_features(x, T_trans)
        slope = np.nanmax(np.abs(np.gradient(T_trans, x)))
        if_delta = max(float(params.get("if_khz", self._REF["if_khz"])) / 1000.0,
                       1e-9)
        if_valid = (x >= x[0] + if_delta) & (x <= x[-1] - if_delta)
        if if_valid.any():
            t_hi = np.interp(x[if_valid] + if_delta, x, T_trans)
            t_lo = np.interp(x[if_valid] - if_delta, x, T_trans)
            if_readout = (t_hi - t_lo) / (2.0 * if_delta)
            i_if = int(np.nanargmax(np.abs(if_readout)))
            if_disc = float(abs(if_readout[i_if]))
            if_detuning = float(x[if_valid][i_if])
        else:
            if_disc = float("nan")
            if_detuning = float("nan")

        metrics = []
        at_resolved = None
        at_split_mhz = None
        at_field_v_m = None
        at_split_field_estimate_v_m = None
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
            peaks = self._transparency_maxima(x, T_trans, window)
            at_resolved = len(peaks) >= 2
            if at_resolved:
                xs = sorted(p[0] for p in sorted(peaks, key=lambda p: p[1])[-2:])
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
                # doublet sits at 0; a detuned microwave pulls it toward the stronger
                # dressed state, so the shift is readable even when one peak fades.
                wsum = sum(h for _, h in peaks)
                center = sum(px * h for px, h in peaks) / wsum
                metrics.append(dict(label="AT center shift", value=f"{center:+.2f} MHz",
                                    help="Height-weighted centre of the dressed "
                                         "transparency peaks; nonzero when the microwave "
                                         "is detuned off the 40D–39F resonance."))
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
            dict(label="IF optimum detuning", value=f"{if_detuning:+.2f} MHz",
                 help="Probe detuning where the finite-IF discriminator is largest."),
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
            superhet = self._superheterodyne_readout(raw, params, x, T_trans)
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
                    help="Probe detuning minimizing total noise-equivalent field, "
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
        return dict(
            x=x, T_trans=T_trans, xphys=xphys, width=width, metrics=metrics,
            superhet=superhet, effective_atoms=self._effective_atom_number(raw, params),
            at_split_mhz=at_split_mhz, at_field_v_m=at_field_v_m,
            at_split_field_estimate_v_m=at_split_field_estimate_v_m, sam=sam)

    def observables(self, raw, params, include_figures=True):
        ro = self._readout(raw, params)
        x, T_trans, xphys, width = ro["x"], ro["T_trans"], ro["xphys"], ro["width"]
        superhet = ro["superhet"]
        effective_atoms = ro["effective_atoms"]
        sam = ro["sam"]
        view = params.get("view", "EIT")

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
                axT.set_xlabel("Frequency [MHz]")
                axT.set_ylabel("Transmission")
                axT.axvline(0.0, color="gray", ls=":", lw=0.8)
                axT.legend(fontsize=9, loc="upper right")
                axT.set_title(
                    f"85Rb Rydberg-EIT: "
                    f"Omega_c = {raw['coupling_rabi_mhz']:.2f} MHz, "
                    f"probe = {raw['probe_power_uw']:.1f} uW")
                xlim = max(2.5, 2.5 * width)
                axT.set_xlim(-xlim, xlim)
            else:
                fig, (axT, axR, axD) = plt.subplots(
                    3, 1, figsize=(8.5, 8.0), sharex=True)
                axT.plot(x, T_trans, color="#0f766e", lw=1.8)
                axT.set_ylabel("Transmission")
                axT.set_ylim(-0.02, 1.04)
                axT.set_title(
                    f"85Rb Rydberg {view}: "
                    f"Omega_c = {raw['coupling_rabi_mhz']:.2f} MHz, "
                    f"Omega_LO = {raw['lo_rabi_mhz']:.2f} MHz")
                response_x = superhet["response_detuning_mhz"]
                response_per_khz = (
                    np.abs(superhet["transmission_phasor_per_angular_rabi"])
                    * 2.0 * np.pi * 1.0e3)
                if response_x.size > 1:
                    axR.plot(response_x, response_per_khz,
                             color="#dc2626", lw=1.6)
                else:
                    axR.plot(response_x, response_per_khz, "o",
                             color="#dc2626")
                axR.axvline(superhet["optimum_detuning_mhz"], color="gray",
                            ls=":", lw=0.8)
                axR.set_ylabel("|dT| per 1 kHz SIG")
                axR.set_title(
                    f"Finite-IF linear response ({float(params.get('if_khz', self._REF['if_khz'])):.0f} kHz IF)")
                axD.plot(x, np.real(xphys), color="#7c3aed", lw=1.5)
                axD.set_ylabel("Re chi")
                axD.set_xlabel("Probe detuning [MHz]")
                for a in (axT, axR, axD):
                    a.axvline(0.0, color="gray", ls=":", lw=0.8)
            fig.tight_layout()

        derived = derived_table([
            ("Ladder", "85Rb 5S1/2 F=3 → 5P3/2 F'=4 → 40D5/2"),
            ("RF leg", "40D5/2 → 39F7/2"),
            ("Microwave frequency", f"{raw['mw_frequency_ghz']:.1f} GHz"),
            ("Coupling Rabi Ω_c / 2π", f"{raw['coupling_rabi_mhz']:.3f} MHz"),
            ("Probe Rabi Ω_P / 2π", f"{raw['probe_rabi_mhz']:.3f} MHz"),
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
                ("Optimum probe detuning", f"{superhet['optimum_detuning_mhz']:+.4f} MHz"),
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
            comparison=comparison)

    def extra_views(self):
        """Ju et al. Fig. 2(b): EIT peak amplitude and linewidth vs probe power,
        with / without B-field compensation, plus the cell-heating sweep used
        to diagnose the temperature-dependent RF operating point."""
        def _compute_sweep(params):
            powers = np.linspace(0.5, 10.0, 16)
            comp = {"width": [], "amp": []}
            uncomp = {"width": [], "amp": []}
            for p_uw in powers:
                pr = dict(params, probe_power_uw=float(p_uw), view="EIT")
                raw = self.compute(pr)
                x, T_c, _ = self._transmission(raw["chi_bar"], raw, pr)
                x, T_u, _ = self._transmission(raw["chi_bar_uncomp"], raw, pr)
                wc, ac = self._eit_features(x, T_c)
                wu, au = self._eit_features(x, T_u)
                comp["width"].append(wc); comp["amp"].append(p_uw * ac)
                uncomp["width"].append(wu); uncomp["amp"].append(p_uw * au)
            return dict(powers=powers.tolist(), comp=comp, uncomp=uncomp)

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
            aW.set_ylabel("EIT linewidth [MHz]")
            aW.set_xlabel("Probe power [uW]")
            figF.tight_layout()
            return figF

        def _compute_temperature_sweep(params):
            # Include the exact nominal values used in the cell-heating deck.
            # This generic panel deliberately treats them as linked simulation
            # temperatures.  The analysis CLI accepts a separate heater ->
            # effective/cold-spot calibration for experimental overlays.
            temperatures = np.array(
                [20.0, 30.0, 35.0, 40.0, 42.0, 45.0, 47.0, 50.0, 53.0, 55.0])
            spectra = []
            widths = []
            contrasts = []
            peaks = []
            slopes = []
            if_discriminators = []
            densities = []
            transit = []
            effective_atoms = []
            psn_sensitivity = []
            total_sensitivity = []
            superhet_detunings = []
            for temp_c in temperatures:
                pr = dict(
                    params, temperature_model="Linked", temp_c=float(temp_c),
                    view="EIT", lo_rabi_mhz=0.0)
                raw = self.compute(pr)
                ro = self._readout(raw, pr)
                x, transmission, _ = self._transmission(
                    raw["chi_bar"], raw, pr)
                width, contrast = self._eit_features(x, transmission)
                metric_values = {
                    metric["label"]: metric["value"]
                    for metric in ro["metrics"]
                }
                ic = int(np.argmin(np.abs(x)))
                spectra.append(transmission.tolist())
                widths.append(float(width))
                contrasts.append(float(contrast))
                peaks.append(float(transmission[ic]))
                slopes.append(float(metric_values["Max spectral slope"].split()[0]))
                if_discriminators.append(float(
                    metric_values["IF discriminator"].split()[0]))
                densities.append(float(raw["N"]))
                transit.append(float(raw["transit_mhz"]))
                effective_atoms.append(float(ro["effective_atoms"].atoms))

                at_lo = max(float(params.get(
                    "lo_rabi_mhz", self._REF["lo_rabi_mhz"])),
                    self._REF["lo_rabi_mhz"])
                at_params = dict(
                    params, temperature_model="Linked", temp_c=float(temp_c),
                    view="AT electrometry", lo_rabi_mhz=at_lo)
                at_raw = self.compute(at_params)
                at_ro = self._readout(at_raw, at_params)
                sensitivity = at_ro["superhet"]["sensitivity"]
                psn_sensitivity.append(float(
                    sensitivity.psn_field_asd_nv_cm_per_sqrt_hz))
                total_sensitivity.append(float(
                    sensitivity.total_field_asd_nv_cm_per_sqrt_hz))
                superhet_detunings.append(float(
                    at_ro["superhet"]["optimum_detuning_mhz"]))
            finite_sensitivity = np.where(
                np.isfinite(total_sensitivity), total_sensitivity, np.inf)
            best_index = int(np.argmin(finite_sensitivity))
            return dict(
                temperatures_c=temperatures.tolist(),
                temperature_semantics="linked effective/cold-spot simulation temperature",
                frequency_mhz=x.tolist(),
                spectra=spectra,
                widths_mhz=widths,
                contrasts=contrasts,
                resonance_transmission=peaks,
                max_slope_per_mhz=slopes,
                if_discriminator_per_mhz=if_discriminators,
                density_m3=densities,
                transit_mhz=transit,
                effective_atoms=effective_atoms,
                psn_sensitivity_nv_cm_sqrt_hz=psn_sensitivity,
                total_sensitivity_nv_cm_sqrt_hz=total_sensitivity,
                superhet_optimum_detuning_mhz=superhet_detunings,
                best_sensitivity_temperature_c=float(temperatures[best_index]),
                best_total_sensitivity_nv_cm_sqrt_hz=float(
                    total_sensitivity[best_index]),
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
            ax_spectrum.set_xlabel("Probe detuning [MHz]")
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
            fig.tight_layout()
            return fig

        return [
            ExtraView(
                key="Fig. 2(b): probe-power dependence (peak amplitude & linewidth)",
                description="Sweeps the probe power 0.5–10 µW and extracts the EIT "
                            "peak amplitude and linewidth with and without B-field "
                            "compensation — the power broadening of Ju et al. Fig. 2(b).",
                compute=_compute_sweep, render=_render_sweep,
            ),
            ExtraView(
                key="Cell-heating sweep: EIT and finite-IF sensitivity",
                description="Sweeps the nominal temperatures used in the cell-heating "
                            "study and reports EIT spectra, linewidth/contrast, "
                            "Gaussian-overlap atom number, and the finite-IF PSN/total "
                            "sensitivity. This generic panel uses linked temperatures; "
                            "use the analysis CLI to apply measured heater-to-effective "
                            "and cold-spot calibrations.",
                compute=_compute_temperature_sweep,
                render=_render_temperature_sweep,
            ),
        ]

    def info(self):
        return (
            "**Rydberg-EIT electrometry.** 85Rb cascade-EIT / microwave Autler-Townes "
            "model calibrated to the photon-shot-noise-limited vapor-cell experiment "
            "of Ju et al. The EIT view reproduces Fig. 2(a) (transmission with / "
            "without B-field compensation); the probe-power panel reproduces "
            "Fig. 2(b). Around the LO-dressed steady state, a finite-IF first-order "
            "Liouvillian solve predicts the weak-SIG transmission phasor. Combined "
            "with the declared RF dipole and balanced-detector noise budget, it "
            "reports PSN, technical, and total field sensitivity. This is not a "
            "full time-domain waveform/lock-in simulation.\n\n"
            "Probe (780 nm) and coupling (481 nm) powers drive Ω_P, Ω_c via √(P)/d "
            "intensity scaling (anchored), and the beam diameter also sets the "
            "transit-time broadening that limits the EIT linewidth (≈1.6 MHz at the "
            "reference). The counter-propagating geometry suppresses the two-photon "
            "Doppler; an optional Doppler-on mode shows the residual "
            "(k_probe − k_coupling)·v broadening.\n\n"
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
