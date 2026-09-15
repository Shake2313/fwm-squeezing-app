"""Temporary, explicitly calibrated coupling closure for the two fast FWM tiers.

This does not change the atomic state or turn a mean-field gain into a quantum
noise prediction. Only the two off-diagonal Maxwell susceptibilities change.
See analysis/fwm_gain_hotfix/DEVLOG.md for calibration and comparison evidence.
"""
import math


MODEL_ID = "sim2025_effective_mixing_participation_v1"
# One fitted number, conditional on the retained 0.74 historical residual.
EFFECTIVE_PARTICIPATION = 0.5594938027


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


def provenance(*, enabled, eligible, w_pump, w_probe, w_conjugate=None):
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
