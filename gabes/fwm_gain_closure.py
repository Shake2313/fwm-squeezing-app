"""Temporary calibrated coupling and power-gain closure for the fast FWM tiers.

This does not change the atomic state or turn a mean-field gain into a quantum
noise prediction. The Gold fit changes two off-diagonal susceptibilities; a
separate OPD-local laboratory fit changes the final power-gain readout only.
See analysis/fwm_gain_hotfix/DEVLOG.md for calibration and comparison evidence.
"""
import math

import numpy as np


MODEL_ID = "sim2025_effective_mixing_participation_v1"
# One fitted number, conditional on the retained 0.74 historical residual.
EFFECTIVE_PARTICIPATION = 0.5594938027

LAB_MODEL_ID = "tpd_20261001_continuous_affine_v1"
# Frozen offline minimax fit: F(G)=a*G+b+c*max(G-k, 0), not per-point lookup.
# In-sample worst relative residual: probe 7.702%, conjugate 7.358%. There are
# no repeat/held-out data, so these are NOT statistical confidence intervals.
# Peak location and other operating conditions remain unvalidated. Replacement
# work is tracked as fwm-tpd-gain-physical-replacement in docs/checklist.json.
LAB_PROBE_AFFINE = (2.64335782913486, -2.320649165812068, 20.159921713415237, 7.0)
LAB_CONJUGATE_AFFINE = (2.4715192372459183, 0.0, 19.933702770307566, 6.1)
LAB_OPD_GHZ = 1.06
# Numerical transition width, NOT a measured detuning tolerance/confidence band.
# Preserve the earlier 0.9-GHz Gold anchor and avoid exporting this fit globally.
LAB_OPD_HALF_WIDTH_GHZ = 0.16


def laboratory_provenance(*, enabled, eligible, D_GHz, branch):
    """Internal ledger for the localized output calibration; no new UI knobs."""
    weight = 0.0
    if enabled and eligible and branch == -1 and D_GHz is not None:
        opd = float(D_GHz)
        distance = abs(opd - LAB_OPD_GHZ) / LAB_OPD_HALF_WIDTH_GHZ
        if (LAB_OPD_GHZ - LAB_OPD_HALF_WIDTH_GHZ < opd
                < LAB_OPD_GHZ + LAB_OPD_HALF_WIDTH_GHZ):
            weight = (1.0 - distance**2)**2
    return {
        "model_id": LAB_MODEL_ID,
        "applied": weight > 0.0,
        "blend_weight": weight,
        "application": "power gains after mean-field propagation and saturation",
        "formula": "F(G)=a*G+b+c*max(G-k,0); output=(1-w)*G+w*F(G)",
        "probe_coefficients_a_b_c_k": list(LAB_PROBE_AFFINE),
        "conjugate_coefficients_a_b_c_k": list(LAB_CONJUGATE_AFFINE),
        "opd_center_GHz": LAB_OPD_GHZ,
        "opd_transition_half_width_GHz": LAB_OPD_HALF_WIDTH_GHZ,
        "transition_is_validation_interval": False,
        "calibration_conditions": {
            "opd_GHz": 1.06, "temperature_C": 118.0, "pump_W": 0.380,
            "cell_input_seed_W": 3.7e-6, "cell_m": 0.0125,
            "pump_waist_m": 530e-6, "probe_waist_m": 330e-6,
            "waist_convention": "1/e^2 intensity radius", "crossing_angle_deg": 0.32,
            "transit_rate_over_2pi_Hz": 100e3, "inherited_residual": 0.74,
            "branch": -1, "eom_interval_GHz": [3.030, 3.039],
        },
        "gain_convention": "measured pump-on output / measured pump-off seed",
        "measured_seed_denominators_uW": [3.7, 3.7, 3.7, 3.8, 3.7, 3.7, 3.7, 3.8, 3.7, 3.7],
        "passive_cell_transmission": None,
        "in_sample_max_relative_residual": {"probe": 0.0770178858941382,
                                             "conjugate": 0.07357959619177322},
        "statistical_confidence_interval": None,
        "independently_validated": False,
        "transfer_matrix_calibrated": False,
        "detector_parameters_used_in_fit": False,
        "extrapolation": "peak location and other EOM/physical conditions unvalidated",
        "development_record": "analysis/fwm_gain_hotfix/tpd_affine_20261001/README_ko.md",
        "replacement_checklist_id": "fwm-tpd-gain-physical-replacement",
    }


def apply_power_gain(G_s, G_c, *, calibration, P_pump, P_seed):
    """Apply a monotone readout fit, preserving no-gain and pump-budget limits.

    This does not calibrate the transfer matrix, microscopic covariance, or a
    physical squeezing spectrum. Downstream algebraic noise remains a diagnostic.
    Inactive contexts return the original arrays without numerical changes.
    """
    weight = calibration["blend_weight"]
    if weight == 0.0:
        return G_s, G_c

    def mapped(gain, coefficients):
        a, b, c, knot = coefficients
        gain = np.asarray(gain, dtype=float)
        fitted = a * gain + b + c * np.maximum(gain-knot, 0.0)
        return gain + weight * (fitted-gain)

    # The upstream soft saturation has already run. Do not saturate twice;
    # enforce only its hard pump-energy envelope if this empirical map hits it.
    cap = max(0.5 * float(P_pump), 0.0) / max(float(P_seed), 1e-30)
    return (np.clip(mapped(G_s, LAB_PROBE_AFFINE), 1.0, 1.0 + cap),
            np.clip(mapped(G_c, LAB_CONJUGATE_AFFINE), 0.0, cap))


def transverse_participation(w_pump, w_probe, w_conjugate=None):
    """Unapplied diagnostic: project I_p/I_p(0) onto Gaussian collected modes.

    For u_j = exp(-r²/w_j²), integral(u_s*u_c*I_p/I_p(0)) /
    integral(u_s*u_c) = a/(a + 2/w_p²), a = 1/w_s² + 1/w_c².
    Matched collected modes give w_p²/(w_p²+w_s²). This is a leading-order
    q proportional to pump-intensity approximation, not a saturated radial solve.
    The existing normalized axial crossing profile remains separate.
    """
    w_conjugate = w_probe if w_conjugate is None else w_conjugate
    waists = tuple(float(w) for w in (w_pump, w_probe, w_conjugate))
    if any(not math.isfinite(w) or w <= 0.0 for w in waists):
        raise ValueError("gain-closure mode waists must be finite and positive")
    wp, ws, wc = waists
    a = 1.0 / ws**2 + 1.0 / wc**2
    return a / (a + 2.0 / wp**2)


def provenance(*, enabled, eligible, w_pump, w_probe, w_conjugate=None,
               D_GHz=None, branch=-1):
    """Fresh, JSON-serializable ledger; no user-adjustable hidden fit parameters."""
    applied = bool(enabled and eligible)
    spatial = transverse_participation(w_pump, w_probe, w_conjugate)
    return {
        "model_id": MODEL_ID,
        "requested": bool(enabled),
        "applied": applied,
        "status": ("semi-empirical estimate; absolute accuracy unknown" if applied
                   else "disabled" if not enabled else "ineligible fidelity"),
        "application": "chi_sc and chi_cs before Maxwell propagation",
        "coupling_multiplier": EFFECTIVE_PARTICIPATION if applied else 1.0,
        "output_calibration": laboratory_provenance(
            enabled=enabled, eligible=eligible, D_GHz=D_GHz, branch=branch),
        "transverse_participation": {
            "value": spatial, "unit": "1",
            "meaning": "Gaussian projected pump intensity / peak pump intensity",
            "source": "analytic Gaussian mode integral; see module docstring",
            "status": "derived under q proportional to local pump intensity",
            "applied": False,
            "uncertainty": None,
            "limitation": "rejected as production law: unsupported pump-waist trend reversal under strong pumping",
        },
        "effective_participation": {
            "value": EFFECTIVE_PARTICIPATION, "unit": "1",
            "meaning": "effective nonlinear mixing susceptibility / reduced-model mixing susceptibility",
            "source": "analysis/fwm_gain_hotfix/exploration_physics/candidate_results.json",
            "status": "fitted to one representative Gold gain",
            "calibration_point": "Sim_2025_gold_final",
            "uncertainty": None,
            "decomposition": "Zeeman, polarization, angular Doppler, Raman and spatial contributions unidentified",
        },
        "calibration": {
            "reference_id": "Sim_2025_gold_final",
            "target_probe_power_gain": 15.5,
            "target_status": "representative literature-compilation midpoint; not a precise measurement",
            "representative_range": [15.0, 16.0],
            "range_is_uncertainty": False,
            "gain_convention": "P_probe_out / P_seed_in; source power gain, no detector loss correction",
            "raw_rounded_probe_power_ratio": 111.0 / 8.0,
            "raw_rounded_conjugate_to_seed_ratio": 109.0 / 8.0,
            "source": "https://www.nature.com/articles/s41598-025-86479-w",
            "repository_provenance": "analysis/fwm_gain_hotfix/reference_points.json",
            "gain_length_diagnostic": {
                "definition": "x = acosh(sqrt(G_s)); ideal-parametric diagnostic only",
                "raw_exact_detuning_gain": 394.29905823951185,
                "x_model": 3.681067366414522,
                "x_target": 2.0470333343160427,
                "x_ratio": 0.5560977647387961,
                "fit_method": "bisection of propagated G_s at exact delta=-8 MHz; not a fit to interpolated display",
            },
            "operating_point": {
                "opd_GHz": 0.9, "tpd_MHz": -8.0, "temperature_C": 121.0,
                "pump_W": 0.6, "seed_W": 8e-6, "cell_m": 0.0125,
                "pump_waist_m": 530e-6, "probe_waist_m": 330e-6,
                "crossing_angle_deg": 0.32,
                "inherited_residual": 0.74, "transit_rate_over_2pi_Hz": 100e3,
            },
        },
        "scope": "Fast/Balanced seeded 85Rb D1; minus-branch Gold calibration only",
        "extrapolation": "other operating points and plus branch: conditional, accuracy unknown",
        "held_out_validation_status": "mixed conditional comparisons; McCormick2008 underpredicted by 4.7x",
        "known_limitations": [
            "Liu2011 probe 6.445 vs 8; McCormick2008 probe 1.926 vs 9",
            "held-out excited-hyperfine detuning references and some beam conventions incomplete",
            "retained local parameter trends checked numerically, not experimentally validated",
        ],
        "physical_squeezing_validated": False,
        "independently_calibrated": False,
        "detector_parameters_used_in_fit": False,
    }


def apply(chi4, multiplier):
    """Scale cross couplings in (ss, cs, sc, cc) order; preserve diagonal drift."""
    if multiplier == 1.0:
        return chi4
    ss, cs, sc, cc = chi4
    return ss, cs * multiplier, sc * multiplier, cc
