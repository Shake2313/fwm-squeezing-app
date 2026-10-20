"""Separate measured pump-off-reference ratios for one recorded FWM profile.

This empirical readout does not modify atomic power gains, susceptibilities,
canonical transfer, noise factors, or detector power propagation.  The measured
seed denominator is the pump-off detector output, not calibrated vapor input.
"""

import copy

import numpy as np
from scipy.interpolate import PchipInterpolator

from . import constants


MODEL_ID = "measured_pump_off_ratio_log_pchip_20261016_v1"
SOURCE_WORKBOOK_SHA256 = "2653fe73abca75550b2b8665da8c8fff16f5c28f155e029f1e6da155a9409ed0"
REFERENCE_RECORD_SHA256 = "0a52cc123edff9749a9d1bcde0649e5e9cd6961961149304e9c165aab7a50b64"
SEED_INPUT_ASSUMPTION = "full_cell_external_incident_proxy_as_atomic_input"
EOM_SUPPORT_MHZ = (3022.0, 3044.0)
REFERENCE_SEED_VOLTAGE_MV = 35.0
REFERENCE_SEED_POWER_W = 3.063457330415755e-6
FULL_CELL_TRANSMISSION = 0.8046
REFERENCE_AFTER_CELL_TRANSMISSION = 0.9709
NOMINAL_SEED_PROXY_W = (
    REFERENCE_SEED_POWER_W/(FULL_CELL_TRANSMISSION*REFERENCE_AFTER_CELL_TRANSMISSION))

_OBSERVED_PROBE_MV = (
    162., 191., 220., 273., 341., 415., 477., 573., 674., 776., 870., 870.,
    820., 690., 510., 320., 160., 60., 20.5, 4.6, 1.91, 1.17, .86)
_OBSERVED_CONJUGATE_MV = (
    130., 150., 190., 242., 311., 397., 452., 551., 656., 772., 850., 890.,
    860., 760., 590., 400., 220., 100., 46., 18., 11.2, 8.2, 6.5)

_PROBE_LOG_INTERPOLATOR = PchipInterpolator(
    np.arange(3022., 3045.),
    np.log(np.asarray(_OBSERVED_PROBE_MV[::-1])/REFERENCE_SEED_VOLTAGE_MV),
    extrapolate=False)
_CONJUGATE_LOG_INTERPOLATOR = PchipInterpolator(
    np.arange(3022., 3045.),
    np.log(np.asarray(_OBSERVED_CONJUGATE_MV[::-1])/REFERENCE_SEED_VOLTAGE_MV),
    extrapolate=False)

_NOMINAL_CONTEXT = {
    "D_GHz": 1.0,
    "T": 386.15,
    "P_pump": 0.294804,
    "P_probe": NOMINAL_SEED_PROXY_W,
    "w_pump": 530e-6,
    "w_probe": 330e-6,
    "L": 0.0125,
    "pump_probe_angle_deg": 0.32,
    "branch": -1,
    "seed_input_assumption": SEED_INPUT_ASSUMPTION,
}


def nominal_context():
    """Return the explicitly conditional comparison context, not measured input.

    The seed proxy estimates external pre-window incident power by dividing a
    pump-off detector reading by the full-cell and after-cell transmissions.
    Using that number as an atomic-medium input is an additional assumption.
    Neither it nor the full-cell transmission determines atomic pump-off A0.
    """
    return copy.deepcopy(_NOMINAL_CONTEXT)


def observations():
    """Fresh raw voltage ratios; no full-cell/QE/noise correction is applied."""
    probe = np.array(_OBSERVED_PROBE_MV)
    conjugate = np.array(_OBSERVED_CONJUGATE_MV)
    return {
        "eom_MHz": np.arange(3044., 3021., -1.),
        "probe_voltage_mV": probe,
        "conjugate_voltage_mV": conjugate,
        "seed_reference_mV": REFERENCE_SEED_VOLTAGE_MV,
        "probe_ratio": probe/REFERENCE_SEED_VOLTAGE_MV,
        "conjugate_ratio": conjugate/REFERENCE_SEED_VOLTAGE_MV,
    }


def profile():
    """Fresh provenance for the measured-reference ratio, distinct from G_s/G_c."""
    return {
        "model_id": MODEL_ID,
        "source_workbook_sha256": SOURCE_WORKBOOK_SHA256,
        "reference_record_sha256": REFERENCE_RECORD_SHA256,
        "reference_record": "docs/grand_challenge/gain_calibration_2026_10_16/interpolation_record.json",
        "method": "deterministic shape-preserving cubic interpolation of natural log ratio",
        "reference_formula": "R_j(f)=exp(PCHIP(f_i,log(V_j_i/35mV))(f))",
        "interpolation_type": "scipy.interpolate.PchipInterpolator; extrapolate=False",
        "node_agreement": "exact reproduction by construction; not a predictive validation result",
        "tuned_parameter_count": 0,
        "sampled_interpolation_validation_passed": True,
        "predictive_calibration_accepted": False,
        "development_acceptance_passed": False,
        "development_joint_max_relative_error": 0.1992246959912889,
        "development_threshold_max_relative_error": 0.2,
        "development_status": "predictive polynomial/spline trials failed; sampled reference interpolation only",
        "development_metric_scope": "post-selection alternating-row interpolation diagnostic; not independent prediction acceptance",
        "alternating_interpolation_diagnostic_passed": True,
        "alternating_interpolation_max_relative_error": {
            "probe": 0.1992246959912889, "conjugate": 0.15544988918714253},
        "failed_predictive_trials": {
            "chebyshev_record": "docs/grand_challenge/gain_calibration_2026_10_16/fit.json",
            "chebyshev_record_sha256": "535cea48a8c977f2507cd47d7840c5daf771dbb8697d51a731a9125e10dd63c3",
            "spline_record": "docs/grand_challenge/gain_calibration_2026_10_16/spline_fit.json",
            "development_acceptance_passed": False,
            "chebyshev_best_joint_max_relative_error": 0.33156074463810215,
        },
        "validation_type": "sampled data interpolation; all observations inspected; no untouched holdout",
        "independent_validation": False,
        "statistical_confidence_interval": None,
        "all_data_fit_max_relative_error": None,
        "source_observations": "Sheet1 A5:D27; two measured voltages divided by B1=35 mV",
        "sample_count": 23,
        "eom_support_MHz": list(EOM_SUPPORT_MHZ),
        "reference_seed_voltage_mV": REFERENCE_SEED_VOLTAGE_MV,
        "reference_seed_power_W": REFERENCE_SEED_POWER_W,
        "gain_quantity": "pump-on arm output / measured pump-off seed output",
        "reference_plane": "pump-off seed detector output after the cell and after-cell optics",
        "numerator_plane": "pump-on probe/conjugate output at the corresponding common-calibration measurement plane",
        "detector_conversion": "same voltage-to-power conversion for both arms, user confirmed",
        "common_after_cell_transmission": REFERENCE_AFTER_CELL_TRANSMISSION,
        "full_cell_transmission_including_windows": FULL_CELL_TRANSMISSION,
        "full_cell_transmission_is_atomic_A0": False,
        "pump_off_atomic_transmission": None,
        "nominal_context": nominal_context(),
        "one_photon_detuning_status": "approximately +1 GHz; numeric context match is not a measured tolerance",
        "input_seed_proxy_W": NOMINAL_SEED_PROXY_W,
        "input_seed_proxy_formula": "P_off_detector/(full_cell_transmission*after_cell_transmission)",
        "input_seed_proxy_plane": "estimated external incident power before entrance window",
        "input_seed_proxy_adoption": SEED_INPUT_ASSUMPTION,
        "independently_calibrated_vapor_input": False,
        "atomic_gain_modified": False,
        "canonical_transfer_modified": False,
        "NF_modified": False,
        "SABES_detector_power_gain_modified": False,
        "physical_squeezing_validated": False,
        "extrapolation": "not provided outside 3022..3044 MHz",
    }


def eom_mhz_from_internal_delta(delta_mhz):
    """Seed below pump: f_EOM=nu_HF-delta_internal, in MHz."""
    return constants.OMEGA_HF/(2*np.pi*1e6)-np.asarray(delta_mhz, dtype=float)


def internal_delta_mhz_from_eom(eom_mhz):
    """Keep the precise hyperfine origin, rather than workbook's rounded 3036."""
    return constants.OMEGA_HF/(2*np.pi*1e6)-np.asarray(eom_mhz, dtype=float)


def reference_curve(eom_mhz):
    """Evaluate the separate empirical ratio on its measured RF support only.

    This ungated mathematical curve is useful for a labelled reference view.
    Application to a selected apparatus context additionally requires
    :func:`context_eligibility`. Inputs and per-arm arrays retain any shape,
    including scalar shape. Outside support and nonfinite inputs return NaN;
    no endpoint clamping, gain floor, source-power cap or extrapolation exists.
    """
    frequency = np.asarray(eom_mhz, dtype=float)
    finite = np.isfinite(frequency)
    supported = finite & (frequency >= EOM_SUPPORT_MHZ[0]) & (frequency <= EOM_SUPPORT_MHZ[1])
    bounded_frequency = np.where(supported, frequency, EOM_SUPPORT_MHZ[0])
    probe = np.where(supported, np.exp(_PROBE_LOG_INTERPOLATOR(bounded_frequency)), np.nan)
    conjugate = np.where(supported, np.exp(_CONJUGATE_LOG_INTERPOLATOR(bounded_frequency)), np.nan)
    status = np.where(supported, "SAMPLED_EMPIRICAL_REFERENCE_INTERPOLATION",
                      "OUTSIDE_MEASURED_RF_SUPPORT")
    status = np.where(finite, status, "INVALID_RF_FREQUENCY")
    return {
        "eom_MHz": frequency.copy(),
        "measured_gain_probe": probe,
        "measured_gain_conjugate": conjugate,
        "within_support": supported,
        "status": status,
        "profile": profile(),
    }


def context_eligibility(context, *, enabled=False):
    """Fail closed unless a caller explicitly adopts the labelled comparison.

    Numeric tolerances below accommodate floating-point unit conversion only.
    They are not experimental uncertainties, acceptance bands, or evidence that
    the approximate OPD is known to that precision.  Detection/noise controls
    do not identify the source context and do not alter this recorded ratio.
    """
    if not isinstance(context, dict):
        raise TypeError("context must be a dictionary of declared physical inputs")
    reasons, comparisons = [], {}
    for key, target in _NOMINAL_CONTEXT.items():
        if key == "seed_input_assumption":
            continue
        if key not in context:
            comparisons[key] = False
            reasons.append(f"missing context input: {key}")
            continue
        try:
            value = float(context[key])
        except (TypeError, ValueError):
            value = float("nan")
        # Relative roundoff allowance plus a target-scaled absolute term; no
        # absolute SI tolerance large enough to broaden a micowatt seed gate.
        tolerance = 1e-10*max(abs(float(target)), np.finfo(float).tiny)
        match = bool(np.isfinite(value) and abs(value-float(target)) <= tolerance)
        comparisons[key] = match
        if not match:
            reasons.append(f"different context input: {key}")
    assumption_matches = context.get("seed_input_assumption") == SEED_INPUT_ASSUMPTION
    if not assumption_matches:
        reasons.append("explicit external-incident seed proxy adoption required")
    for legacy in ("legacy_gain_calibration_fixture", "legacy_mean_field_fixture"):
        if context.get(legacy, False):
            reasons.append(f"historical fixture cannot compose with current reference profile: {legacy}")
    eligible = not reasons
    requested = bool(enabled)
    return {
        "requested": requested,
        "context_matches": bool(all(comparisons.values())),
        "seed_proxy_assumption_adopted": bool(assumption_matches),
        "eligible": eligible,
        "applied": bool(requested and eligible),
        "status": ("DISABLED" if not requested else
                   "CONDITIONAL_MEASURED_REFERENCE_PROFILE" if eligible else
                   "INELIGIBLE_CONTEXT"),
        "reasons": reasons,
        "numeric_context_matches": comparisons,
        "context_tolerance_scope": "floating-point conversion only; not physical confidence or detuning bandwidth",
        "profile": profile(),
    }
