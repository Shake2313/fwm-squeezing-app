"""Cluster D — 85Rb D1 double-Λ four-wave mixing.

The seeded path is a maintained descendant of ``fwm_obe.py``.  It now differs
materially in normalization, Maxwell-sign/phase conventions, collected-mode
conversion, transit reset, and validation metadata.  Its present atomic response
is a finite-reference-field density-matrix Floquet solve, reported at N_F>=3
only after an adjacent-order full-scan check; it remains a mean-field diagnostic,
not a quantum-Langevin model.  A separate slow reference now solves the exact
static pump state and its infinitesimal trace-zero Nambu response for the
standard minus Raman branch; it is not silently substituted into production.

``compute_spectrum`` / ``operating_point`` retain compatibility keys for older
callers.  ``FWMScheme`` wraps them for the generic Streamlit front-end.
"""
import functools
import math
from dataclasses import dataclass

import numpy as np

from .. import atoms, beam, constants, doppler, hyperfine, kernels, observables, species
from ..constants import K_VEC, OMEGA_HF, OMEGA_EXCITED_HF, rabi_freq
from ..core import (
    blas_single_thread,
    build_liouvillian,
    comm_super,
    floquet_solve_truncated,
    liouvillian_pole_residue_response,
    steady_state_batched,
    trace_zero_liouvillian_response,
)
from ..lineshape import fwhm_interp
from .base import ExtraView, ParamSpec, Preset, Scheme

# =========================================================
# Level indices (this scheme's labelling of the atom model)
# =========================================================
ATOM = atoms.get("double_lambda_rb85")
G1, G2, E2, E3 = 0, 1, 2, 3
GROUND_STATES = (G1, G2)
EXCITED_STATES = (E2, E3)
N_LEVELS = ATOM.n_levels
GROUND_F = {G1: 2, G2: 3}
EXCITED_F = {E2: 2, E3: 3}
TRANSITION_DIPOLE_SCALE = np.zeros((N_LEVELS, N_LEVELS), dtype=float)
for _g in GROUND_STATES:
    for _e in EXCITED_STATES:
        TRANSITION_DIPOLE_SCALE[_g, _e] = np.sqrt(
            3.0 * hyperfine.CF2[(GROUND_F[_g], EXCITED_F[_e])])


def physical_coupling_ledger(branch):
    """One-use ledger for the reduced four-level susceptibility scale.

    ``rho_ss`` is trace normalized over the two representative hyperfine
    manifolds, so its ground-state diagonal already supplies the pump-modified
    manifold population. Multiplying by the equilibrium ``p_F`` again would
    duplicate that population. The remaining structural factor is the
    ``1/[2(2I+1)] = 1/12`` sublevel average used by the AutoOD reference.

    The microscopic drive and polarization readout each contain
    ``sqrt(3 C_F^2)``; their product supplies ``3 C_F^2`` exactly once. Isotope
    abundance is absent here because ``hyperfine.number_density`` already
    returns the pure-85Rb cell density. Residual experimental factors are
    applied separately by ``SeededCouplingFactors``.
    """
    if branch not in BRANCHES:
        raise ValueError(f"branch must be one of {BRANCHES}, got {branch}")
    probe_F = GROUND_F[G2] if branch == -1 else GROUND_F[G1]
    return {
        "probe_ground_F": probe_F,
        "atomic_density_source": "pure-85Rb total density",
        "isotope_factor_applied_here": 1.0,
        "manifold_population_source": "trace-normalized rho_ss",
        "equilibrium_manifold_population_reference": hyperfine.GROUND_POP[probe_F],
        "external_manifold_population_factor": 1.0,
        "ground_sublevel_average": 1.0 / hyperfine.N_GROUND_SUBLEVELS,
        "dipole_strength_source": "sqrt(3 C_F^2) in drive and readout",
        "macroscopic_coupling_norm": 1.0 / hyperfine.N_GROUND_SUBLEVELS,
    }


def physical_coupling_norm(branch):
    """Macroscopic factor from :func:`physical_coupling_ledger`."""
    return physical_coupling_ledger(branch)["macroscopic_coupling_norm"]


@dataclass(frozen=True)
class SeededCouplingFactors:
    """Scalar residuals multiplying the seeded-FWM macroscopic coupling.

    ``reference_residual`` preserves an inherited historical factor.  It has not
    been refitted after the normalization and propagation corrections and is not
    an experimental gain or squeezing calibration.
    The three lab-facing terms are one-sided *additional effective-coupling
    penalties*, not independently calibrated power fractions. Keeping their
    defaults at one avoids inventing an arbitrary decomposition of the reference
    residual already absorbed by the anchor.
    """

    reference_residual: float
    mode_overlap_penalty: float = 1.0
    polarization_penalty: float = 1.0
    zeeman_participation_penalty: float = 1.0

    @property
    def lab_factor(self):
        return (float(self.mode_overlap_penalty)
                * float(self.polarization_penalty)
                * float(self.zeeman_participation_penalty))

    @property
    def combined_residual(self):
        return float(self.reference_residual) * self.lab_factor


@dataclass(frozen=True)
class PumpWeakResponseReference:
    """Slow pump-only, infinitesimal Nambu-response reference.

    ``chi_matrix`` has axes
    ``(analysis_frequency, two_photon_detuning, velocity, output, input)``
    with output/input order ``(probe, conjugate-star)``.  Its four entries are
    reduced susceptibilities per unit angular Rabi frequency and therefore have
    units of seconds.  The named ``chi_cs``/``chi_cc`` properties conjugate the
    second Nambu output row to match the historical positive-frequency tuple
    returned by ``chi_matrix_table``.
    """

    branch: int
    delta_axis_rad_s: np.ndarray
    analysis_frequency_axis_rad_s: np.ndarray
    relative_frequency_rad_s: np.ndarray
    delta_eff_axis_rad_s: np.ndarray
    pump_state: np.ndarray
    chi_matrix: np.ndarray
    diagnostics: dict
    provenance: dict

    @property
    def chi_ss(self):
        return self.chi_matrix[..., 0, 0]

    @property
    def chi_cs(self):
        return np.conj(self.chi_matrix[..., 1, 0])

    @property
    def chi_sc(self):
        return self.chi_matrix[..., 0, 1]

    @property
    def chi_cc(self):
        return np.conj(self.chi_matrix[..., 1, 1])


@dataclass(frozen=True)
class NoncollinearPumpWeakResponseReference:
    """Slow ``(v_z,v_x)`` Maxwell average of the pump-only Nambu response.

    The laboratory optical beat, RF analysis frequency, velocity-shifted Raman
    detuning, and one-photon detuning are stored separately.  ``chi_matrix`` has
    axes ``(analysis_frequency, lab_two_photon_detuning, output, input)`` in the
    raw ``(probe, conjugate-star)`` Nambu convention.
    """

    branch: int
    lab_delta_axis_rad_s: np.ndarray
    lab_optical_beat_axis_rad_s: np.ndarray
    analysis_frequency_axis_rad_s: np.ndarray
    one_photon_detuning_rad_s: float
    pump_k_rad_m: float
    probe_k_axis_rad_m: np.ndarray
    crossing_angle_rad: float
    vx_axis_m_s: np.ndarray
    vz_axis_m_s: np.ndarray
    vx_weights: np.ndarray
    vz_weights: np.ndarray
    raman_shift_grid_rad_s: np.ndarray
    one_photon_effective_axis_rad_s: np.ndarray
    pump_state: np.ndarray
    chi_matrix: np.ndarray
    diagnostics: dict
    provenance: dict

    @property
    def chi_ss(self):
        return self.chi_matrix[..., 0, 0]

    @property
    def chi_cs(self):
        return np.conj(self.chi_matrix[..., 1, 0])

    @property
    def chi_sc(self):
        return self.chi_matrix[..., 0, 1]

    @property
    def chi_cc(self):
        return np.conj(self.chi_matrix[..., 1, 1])

# =========================================================
# FWM experiment configuration (cell, beams, detection, scan)
# =========================================================
L_CELL = 12.5e-3
W_PUMP = 530e-6
W_PROBE = 330e-6
P_PUMP, P_PROBE = 600e-3, 10e-6
T_CELL = 394.15

# Sim et al. reference detection inputs. GABES exposes only their product η;
# SABES retains the device-level optical-loss/QE split.
QE_DETECTOR = 0.92
SEEDED_POST_CELL_LOSS_PCT = 5.5
SEEDED_DETECTION_EFFICIENCY_PCT = (
    QE_DETECTOR * (1.0 - SEEDED_POST_CELL_LOSS_PCT / 100.0) * 100.0)
RESPONSIVITY_AW = 0.59
LOSS_FRAC = 0.0
ETA_TOTAL = QE_DETECTOR * (1.0 - LOSS_FRAC)

# Photon energy in eV·nm, so R[A/W] = QE · λ[nm] / 1239.84. Kept for the
# compatibility API and archived readouts.
_HC_EV_NM = 1239.841984
PROBE_WAVELENGTH_NM = 795.0


def responsivity_AW(qe=QE_DETECTOR, wavelength_nm=PROBE_WAVELENGTH_NM):
    """Photodiode responsivity [A/W] implied by a quantum efficiency."""
    return float(qe) * float(wavelength_nm) / _HC_EV_NM

OMEGA_C_SEED = 0.0
SCAN_MIN_GHZ = -8.0
SCAN_MAX_GHZ = 12.0
SCAN_COARSE_POINTS = 401
RESONANCE_WINDOW_MHZ = 80.0
SCAN_FINE_POINTS = 801
PUMP_OVERLAP_EXCLUSION_MHZ = 1e-3
VELOCITY_STEP_MPS = 1.0
VELOCITY_CUTOFF_SIGMA = 3.0
DELTA_GHZ_LIST = [0.9]
BRANCHES = (-1, +1)
DEFAULT_BRANCH = -1
PHASE_LEGACY = "legacy"
PHASE_BALANCED = "balanced"
PHASE_FINE = "fine"
PHASE_ULTRA = "ultra"
SEEDED_PHASE_ANGLE_DEG = 0.32
ULTRA_PHASE_ITERATIONS = 0   # Option A keeps dispersion in diagonal χ, not in Δk
ULTRA_PROPAGATION_SEGMENTS = 64
SEEDED_REFERENCE_RESIDUAL = 0.74
SEEDED_FLOQUET_ORDER = 3
FLOQUET_COMPLEX_RTOL = 0.01
FLOQUET_TRANSFER_ATOL = 1e-10
FLOQUET_RESPONSE_FLOOR_FRACTION = 1e-8
FLOQUET_GAIN_RTOL = 0.01
FLOQUET_PHASE_TOL_DEG = 0.5
FLOQUET_PHASE_AMPLITUDE_FLOOR = 1e-8

# Optional Ultra diagnostic.  The Maxwell diagonal already carries mean-field arm
# absorption, so this path applies only the phenomenological pump-scatter term.
# Distributed loss vacuum and atomic Langevin diffusion remain unavailable.
# ``kappa`` is an inherited, uncalibrated coefficient; disabling the option removes
# that diagnostic term but does not recover a physical squeezing prediction.
HARDENED_PUMP_SCATTER_KAPPA = 0.1


def assess_floquet_scan_convergence(
        *, high_order, low_order, high_response, low_response,
        high_transfer, low_transfer, scan_axis,
        complex_rtol=FLOQUET_COMPLEX_RTOL,
        transfer_atol=FLOQUET_TRANSFER_ATOL,
        response_floor_fraction=FLOQUET_RESPONSE_FLOOR_FRACTION,
        gain_rtol=FLOQUET_GAIN_RTOL,
        phase_tol_deg=FLOQUET_PHASE_TOL_DEG,
        phase_amplitude_floor=FLOQUET_PHASE_AMPLITUDE_FLOOR):
    """Apply the declared full-scan Floquet truncation criteria.

    ``high_response`` and ``low_response`` are mappings of named complex arrays;
    the transfer arrays have shape ``(scan, 2, 2)`` in the canonical photon-flux
    basis.  A single operating point cannot make this gate pass.
    """
    high_order = int(high_order)
    low_order = int(low_order)
    if high_order < 2 or low_order != high_order - 1:
        raise ValueError("Floquet convergence requires adjacent orders N and N-1")
    if set(high_response) != set(low_response):
        raise ValueError("high/low response components must have identical names")
    high_transfer = np.asarray(high_transfer, dtype=complex)
    low_transfer = np.asarray(low_transfer, dtype=complex)
    scan_axis = np.asarray(scan_axis, dtype=float)
    if scan_axis.ndim != 1 or scan_axis.size < 2:
        raise ValueError("a full-scan Floquet gate requires at least two scan points")
    expected_transfer_shape = (scan_axis.size, 2, 2)
    if (high_transfer.shape != expected_transfer_shape
            or low_transfer.shape != expected_transfer_shape):
        raise ValueError(
            f"high/low transfer arrays must have shape {expected_transfer_shape}")
    nonfinite_entries = int(
        high_transfer.size - np.count_nonzero(np.isfinite(high_transfer))
        + low_transfer.size - np.count_nonzero(np.isfinite(low_transfer))
        + scan_axis.size - np.count_nonzero(np.isfinite(scan_axis)))

    complex_components = {}
    complex_pass = True
    for name in sorted(high_response):
        high = np.asarray(high_response[name], dtype=complex)
        low = np.asarray(low_response[name], dtype=complex)
        if high.shape != (scan_axis.size,) or low.shape != (scan_axis.size,):
            raise ValueError(
                f"response component {name!r} must have full-scan shape "
                f"({scan_axis.size},)")
        component_nonfinite = int(
            high.size - np.count_nonzero(np.isfinite(high))
            + low.size - np.count_nonzero(np.isfinite(low)))
        nonfinite_entries += component_nonfinite
        difference = np.abs(high - low)
        finite_magnitudes = np.concatenate((
            np.abs(high[np.isfinite(high)]),
            np.abs(low[np.isfinite(low)]),
        ))
        component_scale = max(
            float(np.max(finite_magnitudes, initial=0.0)),
            np.finfo(float).tiny)
        # Reduced susceptibilities carry seconds, whereas transfer coefficients
        # are dimensionless.  A fixed 1e-10 absolute term would therefore allow
        # order-unity changes in a 1e-12 s response.  Use a full-scan component
        # scale for the response floor and reserve the fixed SI-free atol for T.
        response_floor = float(response_floor_fraction) * component_scale
        allowed = response_floor + float(complex_rtol) * np.maximum(
            np.abs(high), np.abs(low))
        ratio = np.full(difference.shape, np.inf, dtype=float)
        finite = np.isfinite(difference) & np.isfinite(allowed)
        ratio[finite] = difference[finite] / np.maximum(
            allowed[finite], np.finfo(float).tiny)
        passed = bool(np.all(ratio <= 1.0))
        complex_pass = complex_pass and passed
        complex_components[name] = {
            "max_abs_difference": (
                float("inf") if component_nonfinite else
                float(np.max(difference, initial=0.0))),
            "component_scale": component_scale,
            "absolute_floor": response_floor,
            "max_tolerance_ratio": float(np.max(ratio, initial=0.0)),
            "failed_points": int(np.count_nonzero(ratio > 1.0)),
            "nonfinite_entries": component_nonfinite,
            "passed": passed,
        }

    transfer_difference = np.abs(high_transfer - low_transfer)
    transfer_allowed = float(transfer_atol) + float(complex_rtol) * np.maximum(
        np.abs(high_transfer), np.abs(low_transfer))
    transfer_ratio = np.full(transfer_difference.shape, np.inf, dtype=float)
    finite_transfer = np.isfinite(transfer_difference) & np.isfinite(transfer_allowed)
    transfer_ratio[finite_transfer] = transfer_difference[finite_transfer] / np.maximum(
        transfer_allowed[finite_transfer], np.finfo(float).tiny)
    transfer_pass = bool(np.all(transfer_ratio <= 1.0))

    gain_high = {
        "probe_power": np.abs(high_transfer[:, 0, 0]) ** 2,
        "conjugate_photon_flux": np.abs(high_transfer[:, 1, 0]) ** 2,
    }
    gain_low = {
        "probe_power": np.abs(low_transfer[:, 0, 0]) ** 2,
        "conjugate_photon_flux": np.abs(low_transfer[:, 1, 0]) ** 2,
    }
    gain_components = {}
    gain_pass = True
    for name in gain_high:
        relative = np.full(gain_high[name].shape, np.inf, dtype=float)
        finite_gain = np.isfinite(gain_high[name]) & np.isfinite(gain_low[name])
        relative[finite_gain] = (
            np.abs(gain_high[name][finite_gain] - gain_low[name][finite_gain])
            / np.maximum(1.0, np.abs(gain_high[name][finite_gain])))
        passed = bool(np.all(relative < float(gain_rtol)))
        gain_pass = gain_pass and passed
        gain_components[name] = {
            "max_relative_difference": float(np.max(relative, initial=0.0)),
            "failed_points": int(np.count_nonzero(relative >= float(gain_rtol))),
            "passed": passed,
        }

    transfer_scale = np.maximum(
        1.0,
        np.maximum(
            np.max(np.abs(high_transfer), axis=(-2, -1)),
            np.max(np.abs(low_transfer), axis=(-2, -1))),
    )
    amplitude_threshold = float(phase_amplitude_floor) * transfer_scale[:, None, None]
    phase_mask = ((np.abs(high_transfer) > amplitude_threshold)
                  & (np.abs(low_transfer) > amplitude_threshold))
    wrapped_phase = np.degrees(np.angle(high_transfer * np.conj(low_transfer)))
    reported_phase = np.abs(wrapped_phase[phase_mask])
    max_phase = float(np.max(reported_phase, initial=0.0))
    phase_pass = bool(
        nonfinite_entries == 0 and max_phase < float(phase_tol_deg))

    probe_high = gain_high["probe_power"]
    probe_low = gain_low["probe_power"]
    optimum_finite = (
        probe_high.size > 0
        and np.isfinite(scan_axis).all()
        and np.isfinite(probe_high).all()
        and np.isfinite(probe_low).all())
    if optimum_finite:
        high_index = int(np.nanargmax(probe_high))
        low_index = int(np.nanargmax(probe_low))
        optimum_index_shift = abs(high_index - low_index)
        optimum_axis_shift = float(abs(scan_axis[high_index] - scan_axis[low_index]))
        optimum_pass = optimum_index_shift <= 1
    else:
        high_index = low_index = None
        optimum_index_shift = None
        optimum_axis_shift = None
        optimum_pass = False

    passed = bool(
        nonfinite_entries == 0 and complex_pass and transfer_pass
        and gain_pass and phase_pass and optimum_pass)
    return {
        "status": "CONVERGED" if passed else "UNCONVERGED",
        "passed": passed,
        "high_order": high_order,
        "comparison_order": low_order,
        "full_scan_points": int(scan_axis.size),
        "finite": nonfinite_entries == 0,
        "nonfinite_entries": nonfinite_entries,
        "criteria": {
            "complex_rtol": float(complex_rtol),
            "transfer_atol": float(transfer_atol),
            "response_floor_fraction_of_component_scan_max": float(
                response_floor_fraction),
            "gain_rtol": float(gain_rtol),
            "phase_tolerance_deg": float(phase_tol_deg),
            "phase_amplitude_floor": float(phase_amplitude_floor),
            "optimum_max_grid_intervals": 1,
        },
        "response_components": complex_components,
        "transfer": {
            "max_abs_difference": (
                float("inf") if not np.isfinite(transfer_difference).all()
                else float(np.max(transfer_difference, initial=0.0))),
            "max_tolerance_ratio": float(np.max(transfer_ratio, initial=0.0)),
            "failed_coefficients": int(np.count_nonzero(transfer_ratio > 1.0)),
            "passed": transfer_pass,
        },
        "gains": gain_components,
        "wrapped_phase": {
            "max_difference_deg": max_phase,
            "reported_coefficients": int(np.count_nonzero(phase_mask)),
            "passed": phase_pass,
        },
        "probe_gain_optimum": {
            "high_index": high_index,
            "comparison_index": low_index,
            "index_shift": optimum_index_shift,
            "axis_shift": optimum_axis_shift,
            "passed": bool(optimum_pass),
        },
    }


def pump_weak_response_reference_provenance():
    """Machine-readable boundary of the independent pump-only reference."""
    return {
        "solver_id": "static_pump_trace_zero_nambu_reference",
        "state_equation": "self-consistent trace-one 4-level pump Liouvillian",
        "rotating_frame": (
            "U=exp(+i*Omega_beat*t*|g2><g2|), default minus Raman branch"),
        "frame_equivalence": (
            "unitarily equivalent to the pump-only finite-Floquet state; "
            "executable gauge-parity fixture required"),
        "weak_response": (
            "analytic complex-amplitude derivative on the trace-zero subspace"),
        "dc_solver": "trace-constrained stationary-subspace projection",
        "analysis_frequency": "independent angular-frequency input",
        "reference_fields": "none",
        "supported_branches": (-1,),
        "unsupported_branches": {
            +1: (
                "the inherited plus-branch finite-seed frame is not gauge-"
                "equivalent to the physical static pump frame")},
        "production_default": False,
        "doppler_geometry": "longitudinal one-dimensional reference only",
        "quantum_noise": "not implemented",
    }


def noncollinear_doppler_reference_provenance():
    """Machine-readable boundary of the two-dimensional Maxwell reference."""
    return {
        "solver_id": "pump_trace_zero_nambu_2d_maxwell_reference",
        "atomic_response": pump_weak_response_reference_provenance()["solver_id"],
        "velocity_geometry": "tensor (v_z,v_x) truncated-Maxwell quadrature",
        "one_photon_detuning": "Delta_eff=Delta_lab-k_pump*v_z",
        "two_photon_detuning": (
            "delta_eff=delta_lab+(k_pump-k_probe*cos(theta))*v_z"
            "-k_probe*sin(theta)*v_x"),
        "lab_optical_beat": "independent input; never velocity shifted",
        "analysis_frequency": "independent angular-frequency input",
        "floquet_certificate": (
            "minus-branch static pump frame pinned to the N_F=3 pump-only "
            "Floquet state"),
        "supported_branches": (-1,),
        "production_default": False,
        "quadrature": "Gauss-Legendre nodes with explicit Maxwell weights/cutoff",
        "angular_distribution": "single declared crossing angle; no beam divergence",
        "segmentwise_pump_recomputation": False,
        "quantum_noise": "not implemented",
    }


def seeded_atomic_solver_provenance(*, floquet_order=SEEDED_FLOQUET_ORDER,
                                    convergence=None):
    """Machine-readable provenance for the production seeded atomic response."""
    floquet_order = int(floquet_order)
    convergence_status = (
        "not evaluated" if convergence is None else convergence.get("status"))
    return {
        "solver_id": "finite_seed_finite_floquet_density_matrix",
        "state_equation": "full density-matrix Liouvillian null solve",
        "reference_fields": "finite probe/conjugate reference fields",
        "floquet_order": floquet_order,
        "floquet_modes": tuple(range(-floquet_order, floquet_order + 1)),
        "floquet_convergence_status": convergence_status,
        "floquet_comparison_order": (
            None if convergence is None else convergence.get("comparison_order")),
        "pump_only_self_consistent_nullspace": False,
        "pump_only_weak_response_reference_available": True,
        "noncollinear_doppler_reference_available": True,
        "pump_only_reference_solver_id": (
            pump_weak_response_reference_provenance()["solver_id"]),
        "rate_sylvester_approximation": False,
        "scope": "production seeded-FWM response",
        "target_model_status": (
            "pump-only weak-response and two-dimensional non-collinear Doppler "
            "references are separate from production; quantum-noise response "
            "is not implemented"),
    }


def seeded_parameter_provenance(*, line_strength, transit_rate,
                                pump_scatter_kappa=HARDENED_PUMP_SCATTER_KAPPA):
    """Provenance ledger; sweep ranges are illustrative, not uncertainties."""
    return {
        "ell_s": {
            "value": float(line_strength),
            "source": "inherited historical residual",
            "calibration_status": "not refitted after normalization/sign corrections",
            "illustrative_sweep": (0.666, 0.74, 0.814),
            "range_kind": "illustrative sensitivity sweep; not an uncertainty interval",
        },
        "kappa": {
            "value": float(pump_scatter_kappa),
            "source": "phenomenological pump-scatter diagnostic coefficient",
            "calibration_status": "unvalidated and not fitted",
            "illustrative_sweep": (0.0, 0.1, 0.2),
            "range_kind": "illustrative sensitivity sweep; not an uncertainty interval",
        },
        "gamma_transit": {
            "value_rad_s": float(transit_rate),
            "value_hz": float(transit_rate) / (2.0 * np.pi),
            "source": "inherited transit/residual floor",
            "calibration_status": "independent fit unavailable",
            "illustrative_sweep_hz": (90e3, 100e3, 110e3),
            "range_kind": "illustrative sensitivity sweep; not an uncertainty interval",
        },
    }


def seeded_validation_claim_gate(*, canonical_mode_status,
                                 commutator_defect_max,
                                 floquet_convergence=None,
                                 eom_residual_carrier_power=0.0,
                                 eom_other_sidebands_power=0.0,
                                 eom_spectrum_status="not supplied",
                                 eom_spectrum_application="unapplied"):
    """Return the explicit validation boundary for seeded-FWM outputs."""
    defect = float(np.nanmax(np.asarray(commutator_defect_max, dtype=float)))
    reasons = [
        "atomic response uses finite probe/conjugate reference fields",
        "pump-only weak response exists only as a separate slow reference",
        "two-dimensional non-collinear Doppler averaging exists only as a "
        "separate slow reference; production remains one-dimensional",
        "ell_s, kappa, and gamma_transit lack independent calibration",
        "frequency-dependent atomic Langevin diffusion/covariance is unavailable",
    ]
    floquet_status = (
        "NOT_EVALUATED" if floquet_convergence is None
        else str(floquet_convergence.get("status", "NOT_EVALUATED")))
    if floquet_status != "CONVERGED":
        reasons.append(
            f"full-scan Floquet truncation status is {floquet_status}")
    if "conditional" in str(canonical_mode_status).lower():
        reasons.append("conjugate collected-mode area is assumed rather than measured")
    eom_unapplied = str(eom_spectrum_application).strip().lower() != "applied"
    eom_unsupported = str(eom_spectrum_status).strip().lower() not in {
        "supported", "calibrated"}
    if (float(eom_residual_carrier_power) > 0.0
            or float(eom_other_sidebands_power) > 0.0
            or eom_unapplied or eom_unsupported):
        reasons.append(
            "declared EOM residual components are passed through but their "
            f"spectrum model is {eom_spectrum_status}/{eom_spectrum_application}")
    badges = [
        "MEAN_FIELD_DIAGNOSTIC",
        "QUANTITATIVE_GAIN_UNSUPPORTED",
        "PHYSICAL_SQUEEZING_UNAVAILABLE",
    ]
    if floquet_status != "CONVERGED":
        badges.append("FLOQUET_UNCONVERGED")
    return {
        "level": "MEAN_FIELD_DIAGNOSTIC",
        "badges": tuple(badges),
        "mean_field_gain_diagnostic_available": True,
        "quantitative_gain_supported": False,
        "physical_squeezing_prediction": False,
        "experimental_agreement_claim_allowed": False,
        "physical_claims_allowed": False,
        "commutator_defect_max": defect,
        "floquet_convergence_status": floquet_status,
        "blocked_claims": (
            "absolute quantitative gain",
            "physical squeezing spectrum or bandwidth",
            "agreement with experimental gain/squeezing",
        ),
        "reasons": tuple(reasons),
    }


def seeded_sideband_beat(delta, branch):
    """Positive pump-sideband spacing for the seeded Raman branch.

    The scan coordinate is probe = pump + branch*nu_HF + delta.  Therefore the
    physical optical beat magnitude is nu_HF - delta on the (-) branch and
    nu_HF + delta on the (+) branch.
    """
    if branch not in BRANCHES:
        raise ValueError(f"branch must be one of {BRANCHES}, got {branch}")
    return OMEGA_HF + branch * delta


# =========================================================
# Generic SFWM / biphoton topology layer
# =========================================================
NM = 1e-9
MHZ_ANG = 2 * np.pi * 1e6

MODE_SEEDED = "Gain diagnostic"  # stored value; UI label is "Squeezing"
MODE_SEEDED_LEGACY = "Squeezing"  # accepted implicitly by compute() as seeded
MODE_BIPHOTON = "Biphoton"
MODE_LABELS = {
    MODE_SEEDED: "Squeezing",
    MODE_BIPHOTON: "Biphoton",
}

# Biphoton source model. "Predictive" solves the Doppler-averaged cascade/double-Λ
# biphoton amplitude from first principles (Chen et al. PRR 4, 023132 (2024)
# Eq. (3-5); Kim et al. QST 9, 045006 (2024) Eq. (2); Du, Wen, Rubin JOSAB 25,
# C98 (2008)) — waveform, FWHM, bandwidth, OD reshaping and the rate scaling are
# computed, with only one residual scalar per topology setting the absolute rate
# (the squeezing-mode `line_strength` philosophy). "Calibrated" is the legacy
# reference-injected estimate kept for comparison.
BIPHOTON_PREDICTIVE = "Predictive (first-principles)"
BIPHOTON_CALIBRATED = "Calibrated (reference)"
BIPHOTON_MODELS = (BIPHOTON_PREDICTIVE, BIPHOTON_CALIBRATED)
BIPHOTON_MODEL_LABELS = {
    BIPHOTON_PREDICTIVE: "Reduced model",
    BIPHOTON_CALIBRATED: "Reference model",
}

# Two-photon (ground) coherence dephasing as a fraction of the intermediate Γ;
# sets the EIT/Raman two-photon linewidth. Chen et al. fit γ ≈ 0.02–0.03 Γ.
GROUND_DEPHASING_FRAC = 0.02
# Regularization clip for the complex longitudinal function ρ̄ at high OD
# (a numerical regularizer for the separate biphoton longitudinal model).
PRED_RHO_CLIP = 60.0

# Predictive velocity-grid auto-refinement. The nonlinear source |amp(v)| that the
# velocity-class coherent sum integrates is only ~Γ/k — a few m/s — wide, far
# narrower than the Doppler width σ_v the navigate-only `biphoton_velocity_step`
# is sized for. A step coarser than that resonance aliases the (Fourier) sum, so
# the reported absolute BTW width tracks numerical undersampling, not physics
# (e.g. a factor-~20, non-monotonic swing over vstep=1–12 m/s in the 780/1529 nm
# telecom cascade). The predictive path therefore starts from a step that
# oversamples the *measured* resonance, then halves until the |ψ|² FWHM is stable;
# if the point cap is hit first the width is flagged unconverged.
PRED_V_OVERSAMPLE = 16.0      # starting step = probe-measured resonance FWHM / this
PRED_V_FWHM_TOL = 0.03        # relative |ψ|² FWHM change accepted as converged
PRED_V_MAX_REFINE = 10        # max successive halvings of the velocity step
PRED_V_MAX_POINTS = 40000     # velocity-grid point cap (guards runtime)

TOPOLOGY_RB87_TELECOM = "cascade_rb87_telecom"
TOPOLOGY_CS_BTW = "cascade_cs_btw"
TOPOLOGY_DIAMOND = "diamond_generic"
TOPOLOGY_LABELS = {
    TOPOLOGY_RB87_TELECOM: "⁸⁷Rb telecom cascade",
    TOPOLOGY_CS_BTW: "¹³³Cs cascade",
    TOPOLOGY_DIAMOND: "Generic diamond SFWM",
}
CS_CHANNEL_917 = "6D5/2: 852-917 nm"
CS_CHANNEL_795 = "8S1/2: 852-795 nm"
CS_CHANNEL_LABELS = {
    CS_CHANNEL_917: "852/917 nm (6D₅/₂)",
    CS_CHANNEL_795: "852/795 nm (8S₁/₂)",
}
SIDE_PLUS = "+"
SIDE_MINUS = "-"
SIDE_CHOICES = (SIDE_PLUS, SIDE_MINUS)


@dataclass(frozen=True)
class LevelSpec:
    """Lightweight FWM level metadata used by the generic topology layer."""
    name: str
    energy_hz: float
    gamma_mhz: float = 0.0


@dataclass(frozen=True)
class FieldSpec:
    """A driven or generated optical field in a four-wave-mixing topology."""
    role: str
    lower: int
    upper: int
    wavelength_nm: float
    detuning_mhz: float = 0.0
    rabi_mhz: float = 0.0
    phase_sign: float = 1.0
    direction: float = 1.0
    angle_deg: float = 0.0
    side_sign: float = 0.0

    @property
    def k(self):
        return 2 * np.pi / (self.wavelength_nm * NM)

    @property
    def frequency_hz(self):
        return constants.C_LIGHT / (self.wavelength_nm * NM)


@dataclass(frozen=True)
class TopologySpec:
    """Generic SFWM topology; presets carry the reference-calibrated constants."""
    name: str
    label: str
    family: str
    isotope_name: str
    levels: tuple
    fields: tuple
    signal_role: str
    idler_role: str
    default_temp_c: float
    default_cell_mm: float
    default_pump_uw: float
    default_coupling_mw: float
    pair_rate_cps_per_mw: float
    emission_decay_ns: float
    target_g2_peak: float | None = None
    reference_fwhm_ns: float | None = None
    reference_od: float | None = None
    reference_bandwidth_mhz: float | None = None
    reference_width_ratio: float | None = None
    reference_delta_k: float | None = None
    notes: str = ""

    @property
    def isotope(self):
        return species.ISOTOPES[self.isotope_name]

    @property
    def field_map(self):
        return {f.role: f for f in self.fields}


def _wavevector_nm(wavelength_nm):
    return 2 * np.pi / (float(wavelength_nm) * NM)


def _side_sign(value):
    if isinstance(value, str):
        return 1.0 if value.strip() == SIDE_PLUS else -1.0
    return 1.0 if float(value) >= 0.0 else -1.0


def _side_label(value):
    return SIDE_PLUS if float(value) >= 0.0 else SIDE_MINUS


def transverse_matched_angle_deg(source_wavelength_nm, target_wavelength_nm,
                                 source_angle_deg):
    """Collection angle that cancels transverse k for two generated photons."""
    source_k = _wavevector_nm(source_wavelength_nm)
    target_k = _wavevector_nm(target_wavelength_nm)
    x = source_k / target_k * math.sin(math.radians(float(source_angle_deg)))
    return math.degrees(math.asin(float(np.clip(x, -1.0, 1.0))))


def _field_with(field, *, wavelength_nm=None, detuning_mhz=None, rabi_mhz=None,
                angle_deg=None, side_sign=None):
    return FieldSpec(
        role=field.role,
        lower=field.lower,
        upper=field.upper,
        wavelength_nm=field.wavelength_nm if wavelength_nm is None else wavelength_nm,
        detuning_mhz=field.detuning_mhz if detuning_mhz is None else detuning_mhz,
        rabi_mhz=field.rabi_mhz if rabi_mhz is None else rabi_mhz,
        phase_sign=field.phase_sign,
        direction=field.direction,
        angle_deg=field.angle_deg if angle_deg is None else angle_deg,
        side_sign=field.side_sign if side_sign is None else side_sign,
    )


def phase_mismatch(fields, *, signal_angle_deg=None, idler_angle_deg=None,
                   reference_delta_k=0.0):
    """Longitudinal four-field phase mismatch, with an optional reference offset."""
    total = 0.0
    for field in fields:
        angle = field.angle_deg
        if field.role == "signal" and signal_angle_deg is not None:
            angle = signal_angle_deg
        if field.role == "idler" and idler_angle_deg is not None:
            angle = idler_angle_deg
        total += field.phase_sign * field.k * math.cos(math.radians(angle))
    return total - (reference_delta_k or 0.0)


def phase_mismatch_vector(fields, *, signal_angle_deg=None, idler_angle_deg=None,
                          signal_side_sign=None, idler_side_sign=None,
                          reference_delta_k=0.0):
    """Biphoton vector mismatch: calibrated longitudinal, absolute transverse."""
    delta_k_z_absolute = 0.0
    delta_k_x = 0.0
    for field in fields:
        angle = field.angle_deg
        side = field.side_sign
        if field.role == "signal":
            if signal_angle_deg is not None:
                angle = signal_angle_deg
            if signal_side_sign is not None:
                side = signal_side_sign
        elif field.role == "idler":
            if idler_angle_deg is not None:
                angle = idler_angle_deg
            if idler_side_sign is not None:
                side = idler_side_sign
        angle_rad = math.radians(float(angle))
        delta_k_z_absolute += field.phase_sign * field.k * math.cos(angle_rad)
        delta_k_x += field.phase_sign * field.k * math.sin(angle_rad) * side
    delta_k_z_relative = delta_k_z_absolute - (reference_delta_k or 0.0)
    delta_k_vector = math.hypot(delta_k_z_relative, delta_k_x)
    return {
        "delta_k_z_relative": delta_k_z_relative,
        "delta_k_z_absolute": delta_k_z_absolute,
        "delta_k_x": delta_k_x,
        "delta_k_vector": delta_k_vector,
    }


def phase_matching_weight(delta_k, L):
    """sinc^2(delta_k L / 2), normalized to 1 at perfect phase matching."""
    x = 0.5 * np.asarray(delta_k, dtype=float) * L
    out = np.ones_like(x, dtype=float)
    mask = np.abs(x) > 1e-12
    out[mask] = (np.sin(x[mask]) / x[mask]) ** 2
    return out


def _sinc_complex(x):
    """sinc(x) = sin(x)/x for complex argument (→1 at x→0). The longitudinal
    detuning function Φ = sinc(ρ̄)·e^{iρ̄} (Du et al. Eq. 14, Chen et al. Eq. 3)
    carries a complex ρ̄ when in-cell loss/dispersion (OD) is included."""
    x = np.asarray(x, dtype=complex)
    out = np.ones_like(x)
    mask = np.abs(x) > 1e-9
    out[mask] = np.sin(x[mask]) / x[mask]
    return out


def _bandwidth_from_waveform_mhz(tau_s, psi):
    """Spectral FWHM [MHz] of a biphoton temporal waveform ψ(τ) via its FFT."""
    n = np.asarray(tau_s).size
    if n < 4:
        return float("nan")
    dt = float(tau_s[1] - tau_s[0])
    nfft = 4 * n
    spec = np.fft.fftshift(np.fft.fft(np.asarray(psi), n=nfft))
    freq = np.fft.fftshift(np.fft.fftfreq(nfft, dt))     # Hz
    power = np.abs(spec)**2
    if power.max() <= 0:
        return float("nan")
    above = np.where(power >= 0.5 * power.max())[0]
    if above.size < 2:
        return float("nan")
    return float((freq[above[-1]] - freq[above[0]]) / 1e6)


def phase_mismatch_grid(fields, signal_axis_deg, idler_axis_deg,
                        reference_delta_k=0.0):
    """Vectorized signal/idler longitudinal mismatch grid."""
    signal_axis_deg = np.asarray(signal_axis_deg, dtype=float)
    idler_axis_deg = np.asarray(idler_axis_deg, dtype=float)
    sig, ide = np.meshgrid(signal_axis_deg, idler_axis_deg, indexing="ij")
    total = np.zeros_like(sig, dtype=float)
    for field in fields:
        if field.role == "signal":
            total += field.phase_sign * field.k * np.cos(np.deg2rad(sig))
        elif field.role == "idler":
            total += field.phase_sign * field.k * np.cos(np.deg2rad(ide))
        else:
            total += field.phase_sign * field.k * math.cos(math.radians(field.angle_deg))
    return total - (reference_delta_k or 0.0)


def phase_mismatch_vector_grid(fields, signal_axis_deg, idler_axis_deg,
                               reference_delta_k=0.0):
    """Vectorized signal/idler mismatch magnitude grid for biphoton collection."""
    signal_axis_deg = np.asarray(signal_axis_deg, dtype=float)
    idler_axis_deg = np.asarray(idler_axis_deg, dtype=float)
    sig, ide = np.meshgrid(signal_axis_deg, idler_axis_deg, indexing="ij")
    delta_k_z_absolute = np.zeros_like(sig, dtype=float)
    delta_k_x = np.zeros_like(sig, dtype=float)
    for field in fields:
        if field.role == "signal":
            angle_rad = np.deg2rad(sig)
            delta_k_z_absolute += field.phase_sign * field.k * np.cos(angle_rad)
            delta_k_x += (field.phase_sign * field.k * np.sin(angle_rad)
                          * field.side_sign)
        elif field.role == "idler":
            angle_rad = np.deg2rad(ide)
            delta_k_z_absolute += field.phase_sign * field.k * np.cos(angle_rad)
            delta_k_x += (field.phase_sign * field.k * np.sin(angle_rad)
                          * field.side_sign)
        else:
            angle_rad = math.radians(field.angle_deg)
            delta_k_z_absolute += field.phase_sign * field.k * math.cos(angle_rad)
            delta_k_x += (field.phase_sign * field.k * math.sin(angle_rad)
                          * field.side_sign)
    delta_k_z_relative = delta_k_z_absolute - (reference_delta_k or 0.0)
    return np.hypot(delta_k_z_relative, delta_k_x)


def biphoton_phase_matching_map(fields, L, *, signal_angle_deg,
                                 idler_angle_deg, reference_delta_k=0.0,
                                 span_deg=3.0, points=121):
    """Build the optional 2-D signal/idler collection-acceptance map.

    The map is a presentation diagnostic, not an input to the biphoton solve.
    Keeping it in a standalone helper lets headless and batch callers skip the
    2-D allocation while the figure path can request exactly the same map.
    """
    span_deg = max(float(span_deg), 0.0)
    points = max(int(points), 2)
    signal_axis = np.linspace(max(float(signal_angle_deg) - span_deg, 0.0),
                              float(signal_angle_deg) + span_deg, points)
    idler_axis = np.linspace(max(float(idler_angle_deg) - span_deg, 0.0),
                             float(idler_angle_deg) + span_deg, points)
    delta_k = phase_mismatch_vector_grid(
        fields, signal_axis, idler_axis,
        reference_delta_k=reference_delta_k)
    return signal_axis, idler_axis, phase_matching_weight(delta_k, L)


def energy_mismatch_hz(fields):
    signs = {"pump": 1.0, "coupling": 1.0, "signal": -1.0, "idler": -1.0}
    return sum(signs.get(f.role, 0.0) * f.frequency_hz for f in fields)


def _raw_delta_k(fields):
    return phase_mismatch(fields, reference_delta_k=0.0)


def _with_reference_delta_k(spec):
    return TopologySpec(
        name=spec.name,
        label=spec.label,
        family=spec.family,
        isotope_name=spec.isotope_name,
        levels=spec.levels,
        fields=spec.fields,
        signal_role=spec.signal_role,
        idler_role=spec.idler_role,
        default_temp_c=spec.default_temp_c,
        default_cell_mm=spec.default_cell_mm,
        default_pump_uw=spec.default_pump_uw,
        default_coupling_mw=spec.default_coupling_mw,
        pair_rate_cps_per_mw=spec.pair_rate_cps_per_mw,
        emission_decay_ns=spec.emission_decay_ns,
        target_g2_peak=spec.target_g2_peak,
        reference_fwhm_ns=spec.reference_fwhm_ns,
        reference_od=spec.reference_od,
        reference_bandwidth_mhz=spec.reference_bandwidth_mhz,
        reference_width_ratio=spec.reference_width_ratio,
        reference_delta_k=_raw_delta_k(spec.fields),
        notes=spec.notes,
    )


def _rb87_telecom_spec():
    levels = (
        LevelSpec("5S1/2(F=2)", 0.0, 0.0),
        LevelSpec("5P3/2", constants.C_LIGHT / (780.24 * NM), 6.07),
        LevelSpec("4D5/2", constants.C_LIGHT / (780.24 * NM)
                  + constants.C_LIGHT / (1529.37 * NM), 0.66),
        LevelSpec("5P3/2 collection", constants.C_LIGHT / (780.24 * NM), 6.07),
    )
    fields = (
        FieldSpec("pump", 0, 1, 780.24, phase_sign=+1.0),
        FieldSpec("coupling", 1, 2, 1529.37, phase_sign=-1.0, direction=-1.0),
        FieldSpec("signal", 2, 3, 1529.37, phase_sign=-1.0, angle_deg=1.5,
                  side_sign=+1.0),
        FieldSpec("idler", 3, 0, 780.24, phase_sign=+1.0,
                  angle_deg=transverse_matched_angle_deg(1529.37, 780.24, 1.5),
                  side_sign=+1.0),
    )
    return _with_reference_delta_k(TopologySpec(
        name=TOPOLOGY_RB87_TELECOM,
        label="⁸⁷Rb telecom cascade (5S–5P–4D)",
        family="cascade",
        isotope_name="Rb87",
        levels=levels,
        fields=fields,
        signal_role="signal",
        idler_role="idler",
        default_temp_c=90.0,
        default_cell_mm=12.5,
        default_pump_uw=10.0,
        default_coupling_mw=1.0,
        pair_rate_cps_per_mw=38_000.0,
        emission_decay_ns=0.52,
        target_g2_peak=44.0,
        reference_fwhm_ns=0.56,
        reference_od=112.0,
        reference_bandwidth_mhz=300.0,
        notes=("Reference-calibrated cascade SFWM estimate for the telecom "
               "biphoton source in hot 87Rb."),
    ))


def _cs_btw_spec(channel):
    if channel == CS_CHANNEL_795:
        upper = "8S1/2"
        coupling_nm = 795.0
        upper_gamma_mhz = 1.7
        decay_ns = 1.35
        label = "¹³³Cs cascade (852/795 nm)"
    else:
        upper = "6D5/2"
        coupling_nm = 917.0
        upper_gamma_mhz = 2.6
        decay_ns = 4.1
        label = "¹³³Cs cascade (852/917 nm)"
    pump_nm = 852.35
    levels = (
        LevelSpec("6S1/2(F=4)", 0.0, 0.0),
        LevelSpec("6P3/2(F'=5)", constants.C_LIGHT / (pump_nm * NM), 5.23),
        LevelSpec(upper, constants.C_LIGHT / (pump_nm * NM)
                  + constants.C_LIGHT / (coupling_nm * NM), upper_gamma_mhz),
        LevelSpec("6P3/2 collection", constants.C_LIGHT / (pump_nm * NM), 5.23),
    )
    fields = (
        FieldSpec("pump", 0, 1, pump_nm, phase_sign=+1.0),
        FieldSpec("coupling", 1, 2, coupling_nm, phase_sign=-1.0, direction=-1.0),
        FieldSpec("signal", 2, 3, coupling_nm, phase_sign=-1.0, angle_deg=1.5,
                  side_sign=+1.0),
        FieldSpec("idler", 3, 0, pump_nm, phase_sign=+1.0,
                  angle_deg=transverse_matched_angle_deg(coupling_nm, pump_nm, 1.5),
                  side_sign=+1.0),
    )
    return _with_reference_delta_k(TopologySpec(
        name=TOPOLOGY_CS_BTW,
        label=label,
        family="cascade",
        isotope_name="Cs133",
        levels=levels,
        fields=fields,
        signal_role="signal",
        idler_role="idler",
        default_temp_c=75.0,
        default_cell_mm=12.5,
        default_pump_uw=20.0,
        default_coupling_mw=1.0,
        pair_rate_cps_per_mw=12_000.0,
        emission_decay_ns=decay_ns,
        target_g2_peak=18.0,
        reference_fwhm_ns=decay_ns * 0.42,
        reference_od=10.0,
        reference_width_ratio=3.0,
        notes=("Velocity-class coherent-sum model for the Cs biphoton temporal "
               "waveform comparison."),
    ))


def _diamond_generic_spec(params=None):
    params = params or {}
    pump_nm = float(params.get("diamond_pump_nm", 780.0))
    coupling_nm = float(params.get("diamond_coupling_nm", 776.0))
    signal_nm = float(params.get("diamond_signal_nm", 795.0))
    # Energy conservation: 1/λ_idler = 1/λ_pump + 1/λ_coupling − 1/λ_signal.
    # Guard the reciprocal sum — a zero/negative denominator (e.g. pump=coupling
    # with signal=pump/2) has no physical idler and would otherwise raise
    # ZeroDivisionError. Fall back to NaN; the live UI always supplies an explicit
    # diamond_idler_nm, so this default is only the degenerate-case fallback.
    _inv_idler = 1.0 / pump_nm + 1.0 / coupling_nm - 1.0 / signal_nm
    idler_default = 1.0 / _inv_idler if abs(_inv_idler) > 1e-9 else float("nan")
    idler_nm = float(params.get("diamond_idler_nm", idler_default))
    levels = (
        LevelSpec("g", 0.0, 0.0),
        LevelSpec("e1", constants.C_LIGHT / (pump_nm * NM), 6.0),
        LevelSpec("e2", constants.C_LIGHT / (coupling_nm * NM), 6.0),
        LevelSpec("u", constants.C_LIGHT / (pump_nm * NM)
                  + constants.C_LIGHT / (coupling_nm * NM), 1.0),
    )
    fields = (
        FieldSpec("pump", 0, 1, pump_nm, phase_sign=+1.0),
        FieldSpec("coupling", 0, 2, coupling_nm, phase_sign=+1.0),
        FieldSpec("signal", 3, 1, signal_nm, phase_sign=-1.0, angle_deg=2.0,
                  side_sign=+1.0),
        FieldSpec("idler", 3, 2, idler_nm, phase_sign=-1.0,
                  angle_deg=transverse_matched_angle_deg(signal_nm, idler_nm, 2.0),
                  side_sign=-1.0),
    )
    return _with_reference_delta_k(TopologySpec(
        name=TOPOLOGY_DIAMOND,
        label="Generic diamond SFWM",
        family="diamond",
        isotope_name="Rb87",
        levels=levels,
        fields=fields,
        signal_role="signal",
        idler_role="idler",
        default_temp_c=60.0,
        default_cell_mm=12.5,
        default_pump_uw=20.0,
        default_coupling_mw=1.0,
        pair_rate_cps_per_mw=5_000.0,
        emission_decay_ns=8.0,
        target_g2_peak=None,
        reference_fwhm_ns=None,
        reference_od=None,
        reference_bandwidth_mhz=None,
        reference_width_ratio=None,
        notes=("Generic template only; not tied to a validated diamond reference "
               "preset."),
    ))


def _optional_float(value):
    return None if value is None else float(value)


def _topology_cache_key(params):
    return (
        params.get("topology", TOPOLOGY_RB87_TELECOM),
        params.get("cs_channel", CS_CHANNEL_917),
        float(params.get("diamond_pump_nm", 780.0)),
        float(params.get("diamond_coupling_nm", 776.0)),
        float(params.get("diamond_signal_nm", 795.0)),
        _optional_float(params.get("diamond_idler_nm")),
    )


@functools.lru_cache(maxsize=32)
def _topology_from_key(topo, cs_channel, diamond_pump_nm, diamond_coupling_nm,
                       diamond_signal_nm, diamond_idler_nm):
    if topo == TOPOLOGY_CS_BTW:
        return _cs_btw_spec(cs_channel)
    if topo == TOPOLOGY_DIAMOND:
        params = dict(
            diamond_pump_nm=diamond_pump_nm,
            diamond_coupling_nm=diamond_coupling_nm,
            diamond_signal_nm=diamond_signal_nm,
        )
        if diamond_idler_nm is not None:
            params["diamond_idler_nm"] = diamond_idler_nm
        return _diamond_generic_spec(params)
    return _rb87_telecom_spec()


def topology_from_params(params):
    return _topology_from_key(*_topology_cache_key(params))


def _default_biphoton_geometry(params):
    spec = topology_from_params(params)
    signal = spec.field_map["signal"]
    idler = spec.field_map["idler"]
    return dict(
        signal_angle_deg=signal.angle_deg,
        idler_angle_deg=idler.angle_deg,
        signal_side=_side_label(signal.side_sign),
        idler_side=_side_label(idler.side_sign),
    )


class GenericFWMSolver:
    """Reference-calibrated v3 engine for generic SFWM biphoton estimates."""

    def __init__(self, topology):
        self.topology = topology

    def _fields_from_params(self, params):
        pump_rabi = params.get("pump_biphoton_uw", self.topology.default_pump_uw)
        coupling_rabi = params.get("coupling_mw", self.topology.default_coupling_mw)
        out = []
        for field in self.topology.fields:
            rabi = field.rabi_mhz
            detuning = field.detuning_mhz
            if field.role == "pump":
                rabi = beam.anchored_rabi_mhz(
                    2.0, pump_rabi, self.topology.default_pump_uw)
                detuning = params.get("pump_detuning_mhz", 0.0)
            elif field.role == "coupling":
                rabi = beam.anchored_rabi_mhz(
                    12.0, coupling_rabi, self.topology.default_coupling_mw)
                detuning = params.get("coupling_detuning_mhz", 0.0)
            out.append(_field_with(
                field,
                detuning_mhz=detuning,
                rabi_mhz=rabi,
                angle_deg=params.get(f"{field.role}_angle_deg", field.angle_deg),
                side_sign=_side_sign(params[f"{field.role}_side"])
                if field.role in ("signal", "idler")
                and f"{field.role}_side" in params else field.side_sign,
            ))
        return tuple(out)

    def _leg_optical_depth(self, density, L, gamma_mhz, wavelength_nm):
        """Physical on-resonance optical depth of one cascade leg from its natural
        linewidth, via the AutoOD-validated Γ→|d|² route (`species.reduced_dipole_sq`):
        α₀ = 2 N k |d|² /(ε₀ ℏ Γ), OD = α₀ L."""
        gamma_nat = max(float(gamma_mhz), 1e-3) * MHZ_ANG
        lam = float(wavelength_nm) * NM
        k = 2.0 * np.pi / lam
        d2 = species.reduced_dipole_sq(gamma_nat, lam, 0.5, 0.5)
        alpha0 = 2.0 * density * k * d2 / (constants.EPS_0 * constants.HBAR * gamma_nat)
        return float(alpha0 * L)

    def _apply_longitudinal_response(self, kappa_tau, tau_s, v_grid, weights,
                                     residual_k, two_det, Oc, Gamma_e, gamma_g,
                                     od_phys):
        """Convolve the nonlinear response κ̃(τ) with the linear longitudinal
        function Φ̃(τ) (Du Eq. 15), done as a product in the conjugate domain:
        ψ = FFT⁻¹[FFT(κ̃)·Φ(δ)], with Φ(δ)=sinc(ρ̄)·e^{iρ̄} (Du Eq. 14) and ρ̄(δ)
        the OD-weighted EIT/slow-light phase (Chen Eq. 5). At OD→0, ρ̄→0, Φ→1 and
        ψ→κ̃ exactly (no reshaping)."""
        if od_phys <= 0:
            return np.asarray(kappa_tau)
        n = tau_s.size
        dt = float(tau_s[1] - tau_s[0])
        nfft = 4 * n
        K = np.fft.fft(np.asarray(kappa_tau), n=nfft)
        omega = 2.0 * np.pi * np.fft.fftfreq(nfft, dt)        # δ axis [rad/s]
        d = omega[:, None]
        v = v_grid[None, :]
        w = weights[None, :]
        twoden = Oc**2 - 4.0 * (d + 1j * gamma_g) * (
            d + two_det - residual_k * v + 0.5j * Gamma_e)
        rho = (0.5 * od_phys * Gamma_e) * ((d + 1j * gamma_g) / twoden * w).sum(axis=1)
        rho = (np.clip(rho.real, -PRED_RHO_CLIP, PRED_RHO_CLIP)
               + 1j * np.clip(rho.imag, 0.0, PRED_RHO_CLIP))
        phi = _sinc_complex(rho) * np.exp(1j * rho)
        return np.fft.ifft(K * phi, n=nfft)[:n]

    def _predictive_waveform(self, params, fields, v_grid, weights, pump, coupling,
                             signal, idler, residual_k, density, L, pm_weight,
                             tau_axis_ns, od_value):
        """First-principles Doppler-averaged biphoton amplitude.

        Frequency-domain form of Chen et al. (Phys. Rev. Research 4, 023132 (2024))
        Eq. (3-5), equivalent to Kim et al. (Quantum Sci. Technol. 9, 045006 (2024))
        Eq. (2) and Du, Wen, Rubin (J. Opt. Soc. Am. B 25, C98 (2008)) Eq. (13-18):
        for each signal detuning δ the Maxwell velocity classes are coherently
        summed into the nonlinear coupling κ̃(δ) and the linear longitudinal
        function ρ̄(δ); the joint amplitude A(δ)=κ̃·sinc(ρ̄)·e^{iρ̄} is
        inverse-Fourier-transformed to relative time τ. The two-photon denominator
        carries the Ω_c² Autler-Townes term (no weak-coupling approximation), the
        decay envelope emerges from the transform (no injected lifetime), and the
        optical depth α enters κ̃/ρ̄ (group-delay / Sommerfeld-precursor reshaping).
        """
        Gamma_i = max(self.topology.levels[1].gamma_mhz, 0.1) * MHZ_ANG   # intermediate (idler leg)
        Gamma_e = max(self.topology.levels[2].gamma_mhz, 0.1) * MHZ_ANG   # excited Γ₃
        gamma_g = GROUND_DEPHASING_FRAC * Gamma_i                         # two-photon (ground) dephasing
        Op = pump.rabi_mhz * MHZ_ANG
        Oc = max(coupling.rabi_mhz, 1e-6) * MHZ_ANG
        dp = params.get("pump_detuning_mhz", 0.0) * MHZ_ANG
        two_det = (params.get("pump_detuning_mhz", 0.0)
                   + params.get("coupling_detuning_mhz", 0.0)) * MHZ_ANG

        # optical depth seen by the near-resonant idler leg → Chen's α. OD is a
        # measured quantity in these sources (like cell temperature), so use the
        # reference-anchored, density/L-scaled value where available; the in-cell
        # reshaping (ρ̄) it drives is still computed from first principles.
        od_phys = float(od_value)

        tau_s = tau_axis_ns * 1e-9

        # ---- Nonlinear response κ̃(τ): time-domain Kim et al. Eq. (2) ----
        # Per velocity class amp(v) = Ω_p Ω_c / [4·f₁·f₂ + Ω_c²] with the Ω_c²
        # Autler-Townes term in the two-photon denominator (no weak-coupling
        # approximation). The collective two-photon coherence is the coherent sum
        # over Maxwell velocity classes carrying the single-photon phase
        # e^{i k_P v τ}, ×natural decay e^{−Γ τ/2}, ×H(τ) (τ≥0 implicit). The
        # velocity-sum dephasing — not an injected lifetime — sets the BTW width.
        f1 = 0.5 * Gamma_i + 1j * (dp - pump.direction * pump.k * v_grid)
        f2 = 0.5 * Gamma_e + 1j * (two_det - residual_k * v_grid)
        amp_v = weights * (Op * Oc) / (4.0 * f1 * f2 + Oc**2)
        amp_v = amp_v * od_phys                          # κ ∝ α (OD) → rate ∝ OD²
        coherent = np.exp(1j * pump.k * v_grid[:, None] * tau_s[None, :])
        kappa_tau = (amp_v[:, None] * coherent).sum(axis=0)
        kappa_tau = kappa_tau * np.exp(-0.5 * Gamma_i * tau_s)

        # ---- Linear longitudinal response Φ̃(τ): Du Eq. (15) convolution ----
        # ρ̄(δ) (Chen Eq. 5) is the OD-weighted EIT / slow-light phase; the
        # longitudinal function Φ(δ)=sinc(ρ̄)·e^{iρ̄} (Du Eq. 14) reshapes the
        # waveform (group delay / Sommerfeld precursor) at high OD. Convolved with
        # κ̃(τ); at low OD ρ̄→0, Φ→1, ψ→κ̃ (no reshaping). OFF by default: the
        # lumped 4-level model overestimates the high-OD group-delay broadening
        # (it would smear the validated narrow telecom BTW), so the reshaping is a
        # diagnostic opt-in (`biphoton_od_reshaping`) rather than the default path.
        if params.get("biphoton_od_reshaping", False):
            psi_tau = self._apply_longitudinal_response(
                kappa_tau, tau_s, v_grid, weights, residual_k, two_det,
                Oc, Gamma_e, gamma_g, od_phys)
        else:
            psi_tau = kappa_tau
        psi_tau = psi_tau * math.sqrt(max(pm_weight, 0.0))   # transverse collection

        # Source spectral width from the waveform itself (predictive).
        bandwidth_mhz = _bandwidth_from_waveform_mhz(tau_s, psi_tau)

        # Du regime split: group-delay time τ_g≈(2γ/Ω_c²)·OD·Γ vs Rabi time 2π/Ω_c.
        tau_group = (2.0 * gamma_g / max(Oc**2, 1e-30)) * od_phys * Gamma_e
        tau_rabi = 2.0 * np.pi / max(Oc, 1e-30)
        regime = "group-delay" if tau_group > tau_rabi else "damped-Rabi"

        source_v = amp_v   # velocity-class source for the existing diagnostic plot

        return {
            "psi_tau": psi_tau,
            "source_v": source_v,
            "od_phys": float(od_phys),
            "bandwidth_mhz": float(bandwidth_mhz),
            "regime": regime,
            # Chen et al. ultimate spectral-brightness ceiling ≈ (π/2)·10⁶ pairs/s/MHz
            "brightness_limit_cps_per_mhz": float(0.5 * np.pi * 1e6),
        }

    def _predictive_velocity_step(self, params, pump, coupling, residual_k, T, iso):
        """Velocity-grid step that oversamples the velocity-space resonance of the
        nonlinear source |amp(v)| — the integrand of the velocity-class coherent
        sum in `_predictive_waveform`.

        The resonance is measured directly on a fine probe grid (so it follows the
        Autler-Townes / detuning broadening of the denominator, not just the bare
        linewidth), and the step is set to its FWHM / `PRED_V_OVERSAMPLE`. Falls
        back to the natural-linewidth width Γ/2k if the probe feature is
        ill-defined. This is only the *starting* step; `_converged_predictive_
        waveform` halves it further until the waveform width converges.
        """
        sigma_v = math.sqrt(constants.KB * T / iso.mass)
        Gamma_i = max(self.topology.levels[1].gamma_mhz, 0.1) * MHZ_ANG
        Gamma_e = max(self.topology.levels[2].gamma_mhz, 0.1) * MHZ_ANG
        Oc = max(coupling.rabi_mhz, 1e-6) * MHZ_ANG
        Op = pump.rabi_mhz * MHZ_ANG
        dp = params.get("pump_detuning_mhz", 0.0) * MHZ_ANG
        two_det = (params.get("pump_detuning_mhz", 0.0)
                   + params.get("coupling_detuning_mhz", 0.0)) * MHZ_ANG
        kp = max(abs(pump.k), 1e-30)
        rk = max(abs(residual_k), 1e-30)
        # probe spacing resolves the narrowest bare-linewidth velocity scale
        narrow = min(Gamma_i / (2.0 * kp), Gamma_e / (2.0 * rk))
        dv_probe = max(narrow / 6.0, 1e-3)
        n_probe = int(min(2.0 * 3.2 * sigma_v / dv_probe + 1.0, 200000))
        v = np.linspace(-3.2 * sigma_v, 3.2 * sigma_v, max(n_probe, 64))
        w = np.exp(-v**2 / (2.0 * sigma_v**2))
        f1 = 0.5 * Gamma_i + 1j * (dp - pump.direction * pump.k * v)
        f2 = 0.5 * Gamma_e + 1j * (two_det - residual_k * v)
        integ = np.abs(w * (Op * Oc) / (4.0 * f1 * f2 + Oc**2))
        res_fwhm = fwhm_interp(v, integ)
        if not np.isfinite(res_fwhm) or res_fwhm <= 0:
            res_fwhm = Gamma_i / (2.0 * kp)
        return max(res_fwhm / PRED_V_OVERSAMPLE, 1e-3)

    def _converged_predictive_waveform(self, params, fields, pump, coupling, signal,
                                       idler, residual_k, density, L, pm_weight,
                                       tau_axis_ns, od_value, T, iso, user_step):
        """Predictive biphoton waveform on a velocity grid auto-refined to
        convergence.

        The coherent sum over Maxwell velocity classes is a discretized Fourier
        integral of the narrow source |amp(v)|; a step coarser than that resonance
        aliases it, so the navigate-only `biphoton_velocity_step` (sized for the
        σ_v-wide Doppler profile) is far too coarse and the absolute BTW width
        swings with it. Start from a step that oversamples the measured resonance —
        never coarser than the user step — then halve until the |ψ|² FWHM is stable
        within `PRED_V_FWHM_TOL`. Returns the converged waveform, its velocity grid
        and weights, the step used, and a convergence flag (False if the point cap
        `PRED_V_MAX_POINTS` is reached first → width is qualitative only).
        """
        def solve_on(v_grid, weights):
            wf = self._predictive_waveform(
                params, fields, v_grid, weights, pump, coupling, signal, idler,
                residual_k=residual_k, density=density, L=L, pm_weight=pm_weight,
                tau_axis_ns=tau_axis_ns, od_value=od_value)
            fw = float(fwhm_interp(tau_axis_ns, np.abs(wf["psi_tau"])**2))
            return wf, fw

        dv = min(float(user_step),
                 self._predictive_velocity_step(params, pump, coupling,
                                                residual_k, T, iso))
        v_grid, weights = doppler.velocity_grid(T, dv=dv, cutoff_sigma=3.0,
                                                mass=iso.mass)
        wf, fw = solve_on(v_grid, weights)
        converged = False
        for _ in range(PRED_V_MAX_REFINE):
            step2 = 0.5 * dv
            v2, w2 = doppler.velocity_grid(T, dv=step2, cutoff_sigma=3.0,
                                           mass=iso.mass)
            if v2.size > PRED_V_MAX_POINTS:
                break
            wf2, fw2 = solve_on(v2, w2)
            rel = abs(fw2 - fw) / max(fw2, 1e-12)
            v_grid, weights, wf, dv, fw = v2, w2, wf2, step2, fw2
            if np.isfinite(rel) and rel < PRED_V_FWHM_TOL:
                converged = True
                break
        return wf, v_grid, weights, float(dv), bool(converged)

    def compute_biphoton(self, params):
        T = params.get("biphoton_temp_c", self.topology.default_temp_c) + 273.15
        L = params.get("biphoton_cell_mm", self.topology.default_cell_mm) * 1e-3
        fields = self._fields_from_params(params)
        fmap = {f.role: f for f in fields}
        detail = params.get("phase_detail", "Balanced")
        phase_detail = "Fine" if str(detail).lower() == "fine" else "Balanced"
        iso = self.topology.isotope
        v_step = params.get("biphoton_velocity_step", 2.0)
        pump = fmap["pump"]
        coupling = fmap["coupling"]
        signal = fmap["signal"]
        idler = fmap["idler"]
        delta_k = phase_mismatch(
            fields,
            signal_angle_deg=signal.angle_deg,
            idler_angle_deg=idler.angle_deg,
            reference_delta_k=self.topology.reference_delta_k,
        )
        pm_weight_longitudinal = float(phase_matching_weight(np.array([delta_k]), L)[0])
        delta_k_absolute = phase_mismatch(
            fields,
            signal_angle_deg=signal.angle_deg,
            idler_angle_deg=idler.angle_deg,
            reference_delta_k=0.0,
        )
        pm_weight_absolute = float(phase_matching_weight(
            np.array([delta_k_absolute]), L)[0])
        vector_pm = phase_mismatch_vector(
            fields,
            signal_angle_deg=signal.angle_deg,
            idler_angle_deg=idler.angle_deg,
            reference_delta_k=self.topology.reference_delta_k,
        )
        pm_weight = float(phase_matching_weight(
            np.array([vector_pm["delta_k_vector"]]), L)[0])
        density = species.number_density(iso, T)
        default_density = species.number_density(iso, self.topology.default_temp_c + 273.15)
        residual_k = pump.direction * pump.k + coupling.direction * coupling.k
        if self.topology.reference_od is not None:
            od_estimate = (self.topology.reference_od
                           * density / max(default_density, 1e-30)
                           * L / max(self.topology.default_cell_mm * 1e-3, 1e-30))
        else:
            od_estimate = np.nan

        # Absolute pair rate is reference-anchored. Its current scalar law uses
        # pump power, coupling drive, vector phase matching and sqrt(N/N_ref);
        # OD/cell length shape the waveform path but do not directly scale this
        # rate. The lumped model does not pin the collection coefficient, so the
        # magnitude remains anchored like the squeezing `line_strength` residual.
        pump_mw = params.get("pump_biphoton_uw", self.topology.default_pump_uw) * 1e-3
        coupling_scale = max(params.get("coupling_mw", self.topology.default_coupling_mw),
                             0.0) / max(self.topology.default_coupling_mw, 1e-12)
        pair_rate = (self.topology.pair_rate_cps_per_mw * pump_mw
                     * math.sqrt(max(coupling_scale, 0.0)) * pm_weight
                     * math.sqrt(max(density, 1e-30) / max(default_density, 1e-30)))

        model = params.get("biphoton_model", BIPHOTON_PREDICTIVE)
        predictive = model == BIPHOTON_PREDICTIVE
        tau_axis_ns = np.linspace(0.0, params.get("tau_max_ns", 12.0), 481)

        if predictive:
            od_value = (od_estimate if np.isfinite(od_estimate)
                        else self._leg_optical_depth(
                            density, L, self.topology.levels[1].gamma_mhz,
                            idler.wavelength_nm))
            # The velocity-class coherent sum aliases on a step coarser than the
            # narrow source resonance, so the predictive grid is auto-refined to
            # convergence (the user `biphoton_velocity_step` acts as an upper bound).
            wf, v_grid, weights, v_step, velocity_converged = (
                self._converged_predictive_waveform(
                    params, fields, pump, coupling, signal, idler,
                    residual_k=residual_k, density=density, L=L,
                    pm_weight=pm_weight, tau_axis_ns=tau_axis_ns,
                    od_value=od_value, T=T, iso=iso, user_step=v_step))
            psi_tau = wf["psi_tau"]
            source_v = wf["source_v"]
            od_estimate = wf["od_phys"]
            source_bandwidth_mhz = wf["bandwidth_mhz"]
            regime = wf["regime"]
            brightness_limit = wf["brightness_limit_cps_per_mhz"]
        else:
            v_grid, weights = doppler.velocity_grid(
                T, dv=v_step, cutoff_sigma=3.0, mass=iso.mass)
            lower_gamma = max(self.topology.levels[1].gamma_mhz, 0.1) * MHZ_ANG
            upper_gamma = max(self.topology.levels[2].gamma_mhz, 0.1) * MHZ_ANG
            pump_det = params.get("pump_detuning_mhz", 0.0) * MHZ_ANG
            two_det = (params.get("pump_detuning_mhz", 0.0)
                       + params.get("coupling_detuning_mhz", 0.0)) * MHZ_ANG
            lower = lower_gamma / 2.0 + 1j * (pump_det - pump.direction * pump.k * v_grid)
            upper = upper_gamma / 2.0 + 1j * (two_det - residual_k * v_grid)
            drive = (pump.rabi_mhz * coupling.rabi_mhz) * (MHZ_ANG ** 2)
            source_v = weights * drive / (lower * upper)
            source_v *= math.sqrt(max(pm_weight, 0.0))
            tau_s = tau_axis_ns * 1e-9
            phase_k = abs(idler.k)
            coherent = np.exp(1j * phase_k * v_grid[:, None] * tau_s[None, :])
            psi_tau = (source_v[:, None] * coherent).sum(axis=0)
            psi_tau *= np.exp(-tau_axis_ns / max(self.topology.emission_decay_ns, 1e-12))
            if np.max(np.abs(psi_tau)) > 0:
                psi_tau = psi_tau / np.max(np.abs(psi_tau))
            source_bandwidth_mhz = float(self.topology.reference_bandwidth_mhz or 300.0)
            regime = "calibrated"
            brightness_limit = float("nan")
            # Calibrated width is set by the injected emission lifetime, not the
            # velocity-sum dephasing, so it is stable on the user grid (no refine).
            velocity_converged = True

        return {
            "kind": "biphoton",
            "topology": self.topology,
            "fields": fields,
            "tau_axis_ns": tau_axis_ns,
            "psi_tau": psi_tau,
            "v_grid": v_grid,
            "velocity_weights": weights,
            "source_v": source_v,
            "angle_axis_deg": None,
            "phase_matching": None,
            # The optional 2-D map is built lazily by the figure path. Keep the
            # legacy keys as None so existing raw-result consumers fail soft.
            "signal_angle_axis_deg": None,
            "idler_angle_axis_2d_deg": None,
            "phase_matching_2d": None,
            "delta_k": float(delta_k),
            "delta_k_absolute": float(delta_k_absolute),
            "delta_k_z_relative": float(vector_pm["delta_k_z_relative"]),
            "delta_k_z_absolute": float(vector_pm["delta_k_z_absolute"]),
            "delta_k_x": float(vector_pm["delta_k_x"]),
            "delta_k_vector": float(vector_pm["delta_k_vector"]),
            "phase_match_weight": pm_weight,
            "phase_match_weight_vector": pm_weight,
            "phase_match_weight_longitudinal": pm_weight_longitudinal,
            "phase_match_weight_absolute": pm_weight_absolute,
            "phase_detail": phase_detail,
            "energy_mismatch_hz": float(energy_mismatch_hz(fields)),
            "pair_rate_cps": float(pair_rate),
            "density": float(density),
            "od_estimate": float(od_estimate),
            "source_bandwidth_mhz": float(source_bandwidth_mhz),
            "temperature_K": float(T),
            "cell_length_m": float(L),
            "residual_two_photon_k": float(residual_k),
            "biphoton_model": model,
            "predictive": bool(predictive),
            "regime": regime,
            "brightness_limit_cps_per_mhz": float(brightness_limit),
            "velocity_step_used": float(v_step),
            "velocity_converged": bool(velocity_converged),
            "notes": self.topology.notes,
        }


@functools.lru_cache(maxsize=32)
def _solver_from_topology_key(key):
    return GenericFWMSolver(_topology_from_key(*key))


# =========================================================
# Hamiltonians
# =========================================================
def _add_static_drive(H, ground, omega):
    for excited in EXCITED_STATES:
        omega_ge = omega * TRANSITION_DIPOLE_SCALE[ground, excited]
        H[ground, excited] += omega_ge / 2
        H[excited, ground] += omega_ge / 2


def _add_sideband_drive(H, ground, omega):
    for excited in EXCITED_STATES:
        H[excited, ground] += omega * TRANSITION_DIPOLE_SCALE[ground, excited] / 2


def _polarization_coherence(rho, ground):
    return sum(TRANSITION_DIPOLE_SCALE[ground, e] * rho[:, e, ground]
               for e in EXCITED_STATES)


def pump_hamiltonian_at_Deff_zero(Op_A, Op_B):
    """Static physical pump frame, before the longitudinal Doppler shift.

    For the standard minus Raman branch, applying
    ``U(t)=exp(+i*Omega_beat*t*|g2><g2|)`` to the seeded frame makes the second
    pump leg static and changes its diagonal from ``delta`` to
    ``delta + Omega_beat = OMEGA_HF``.  The resulting pump Hamiltonian is
    independent of probe detuning, as a pump-only state must be.
    """
    H = np.zeros((N_LEVELS, N_LEVELS), dtype=complex)
    H[G2, G2] = OMEGA_HF
    H[E2, E2] = -OMEGA_EXCITED_HF
    H[E3, E3] = 0.0
    _add_static_drive(H, G1, Op_A)
    _add_static_drive(H, G2, Op_B)
    return H


def pump_frame_to_seeded_harmonics(pump_state, n_f=1, branch=-1):
    """Reconstruct seeded-frame pump harmonics from the static pump state.

    This executable gauge map is currently certified only for the standard
    minus Raman branch.  In that frame coherences with ``g2`` in the row occupy
    harmonic ``-1`` and those with ``g2`` in the column occupy ``+1``; every
    ``|n|>=2`` block is identically zero.
    """
    if branch != -1:
        raise NotImplementedError(
            "the inherited plus-branch seeded frame is not gauge-equivalent "
            "to the physical static pump frame")
    n_f = int(n_f)
    if n_f < 1:
        raise ValueError("n_f must be at least 1")
    rho = np.asarray(pump_state, dtype=complex)
    if rho.shape[-2:] != (N_LEVELS, N_LEVELS):
        raise ValueError(
            f"pump_state must end in shape {(N_LEVELS, N_LEVELS)}")
    harmonics = np.zeros(
        rho.shape[:-2] + (2 * n_f + 1, N_LEVELS, N_LEVELS), dtype=complex)
    for row in range(N_LEVELS):
        for column in range(N_LEVELS):
            if row == G2 and column != G2:
                harmonic = -1
            elif row != G2 and column == G2:
                harmonic = +1
            else:
                harmonic = 0
            harmonics[..., n_f + harmonic, row, column] = rho[..., row, column]
    return harmonics


def _pump_nambu_operators(branch):
    """Return weak Nambu drive superoperators and polarization readouts."""
    if branch != -1:
        raise NotImplementedError(
            "pump-only weak response is certified only for the standard minus "
            "Raman branch; the inherited plus-branch frame fails gauge parity")
    probe_ground = G2
    conjugate_ground = G1

    probe_raise = np.zeros((N_LEVELS, N_LEVELS), dtype=complex)
    conjugate_star_lower = np.zeros_like(probe_raise)
    readouts = np.zeros((2, N_LEVELS * N_LEVELS), dtype=complex)
    for excited in EXCITED_STATES:
        probe_scale = TRANSITION_DIPOLE_SCALE[probe_ground, excited]
        conjugate_scale = TRANSITION_DIPOLE_SCALE[conjugate_ground, excited]
        probe_raise[excited, probe_ground] = 0.5 * probe_scale
        conjugate_star_lower[conjugate_ground, excited] = 0.5 * conjugate_scale
        readouts[0, excited * N_LEVELS + probe_ground] = probe_scale
        readouts[1, conjugate_ground * N_LEVELS + excited] = conjugate_scale
    drives = np.stack((comm_super(probe_raise),
                       comm_super(conjugate_star_lower)))
    return drives, readouts


def _pump_reference_system(Op_A, Op_B, Delta_eff_axis, atom, branch):
    """Build the shared state/source/readout pieces for the slow reference."""
    if branch != -1:
        _pump_nambu_operators(branch)  # raises the documented frame error
    deff = np.asarray(Delta_eff_axis, dtype=float)
    if deff.ndim != 1 or deff.size < 1 or not np.isfinite(deff).all():
        raise ValueError("Delta_eff_axis must be a finite non-empty 1-D array")
    H_pump = pump_hamiltonian_at_Deff_zero(Op_A, Op_B)
    L0 = build_liouvillian(H_pump, atom)
    L_batch = L0[None, :, :] - deff[:, None, None] * atom.S_v[None, :, :]
    rho = steady_state_batched(L0, deff, atom.S_v, N_LEVELS)
    rho_vec = rho.reshape(deff.size, N_LEVELS * N_LEVELS)
    drives, readouts = _pump_nambu_operators(branch)
    sources = np.einsum("rij,bj->bir", drives, rho_vec)
    return deff, L0, L_batch, rho, rho_vec, sources, readouts


def pump_only_weak_response_reference(
        Op_A, Op_B, delta_axis, Delta_eff_axis, *, branch=-1, atom=None,
        analysis_frequency_axis_rad_s=(0.0,), T=T_CELL):
    """Solve the self-consistent pump state and infinitesimal 2x2 response.

    ``delta_axis`` is the optical two-photon detuning and
    ``analysis_frequency_axis_rad_s`` is an independent RF/sideband analysis
    frequency.  They are never folded into one input.  In the static physical
    pump frame the Nambu sector frequency is

    ``omega_relative = -OMEGA_HF + delta + Omega_SA``

    for the supported minus Raman branch.  If ``atom`` is omitted, the reference
    uses the same temperature-dependent collisional and thermal-transit-reset
    model as :func:`compute_spectrum`.  The DC singularity is handled by the
    trace-zero stationary-subspace projection in
    :func:`gabes.core.trace_zero_liouvillian_response`.
    """
    if branch != -1:
        _pump_nambu_operators(branch)  # raises with the explicit limitation
    delta = np.atleast_1d(np.asarray(delta_axis, dtype=float))
    analysis = np.atleast_1d(
        np.asarray(analysis_frequency_axis_rad_s, dtype=float))
    if delta.ndim != 1 or delta.size < 1 or not np.isfinite(delta).all():
        raise ValueError("delta_axis must be a finite non-empty 1-D array")
    if analysis.ndim != 1 or analysis.size < 1 or not np.isfinite(analysis).all():
        raise ValueError(
            "analysis_frequency_axis_rad_s must be a finite non-empty 1-D array")
    if atom is None:
        atom = collisional_atom(T)

    (deff, _L0, L_batch, rho, rho_vec, sources, readouts) = (
        _pump_reference_system(
            Op_A, Op_B, Delta_eff_axis, atom, branch))
    relative = -OMEGA_HF + delta[None, :] + analysis[:, None]
    chi = np.empty(
        (analysis.size, delta.size, deff.size, 2, 2), dtype=complex)
    max_response_trace = 0.0
    max_response_residual = 0.0
    eye = np.eye(N_LEVELS * N_LEVELS, dtype=complex)
    diagonal = np.arange(N_LEVELS) * (N_LEVELS + 1)
    for analysis_index in range(analysis.size):
        for delta_index in range(delta.size):
            omega = float(relative[analysis_index, delta_index])
            response = trace_zero_liouvillian_response(
                L_batch, sources, omega, N_LEVELS)
            chi[analysis_index, delta_index] = np.einsum(
                "oi,bik->bok", readouts, response)
            traces = np.sum(response[:, diagonal, :], axis=1)
            max_response_trace = max(
                max_response_trace,
                float(np.max(np.abs(traces), initial=0.0)))
            residual = np.einsum(
                "bij,bjk->bik", L_batch + 1j * omega * eye, response)
            residual += sources
            scale = (
                np.linalg.norm(
                    np.einsum("bij,bjk->bik", L_batch, response), axis=1)
                + abs(omega) * np.linalg.norm(response, axis=1)
                + np.linalg.norm(sources, axis=1))
            normalized = np.linalg.norm(residual, axis=1) / np.maximum(
                scale, np.finfo(float).tiny)
            max_response_residual = max(
                max_response_residual,
                float(np.max(normalized, initial=0.0)))

    pump_residual = np.einsum("bij,bj->bi", L_batch, rho_vec)
    pump_scale = np.linalg.norm(L_batch, axis=(1, 2)) * np.linalg.norm(
        rho_vec, axis=1)
    pump_normalized = np.linalg.norm(pump_residual, axis=1) / np.maximum(
        pump_scale, np.finfo(float).tiny)
    traces = np.trace(rho, axis1=-2, axis2=-1)
    hermiticity = rho - rho.conj().swapaxes(-1, -2)
    minimum_eigenvalue = float(np.min(np.linalg.eigvalsh(rho)))
    diagnostics = {
        "finite": bool(np.isfinite(rho).all() and np.isfinite(chi).all()),
        "max_pump_normalized_residual": float(
            np.max(pump_normalized, initial=0.0)),
        "max_pump_trace_error": float(
            np.max(np.abs(traces - 1.0), initial=0.0)),
        "max_pump_hermiticity_error": float(
            np.max(np.abs(hermiticity), initial=0.0)),
        "minimum_pump_eigenvalue": minimum_eigenvalue,
        "max_response_normalized_residual": max_response_residual,
        "max_response_trace_error": max_response_trace,
        "velocity_classes": int(deff.size),
        "two_photon_points": int(delta.size),
        "analysis_frequency_points": int(analysis.size),
    }
    return PumpWeakResponseReference(
        branch=branch,
        delta_axis_rad_s=delta.copy(),
        analysis_frequency_axis_rad_s=analysis.copy(),
        relative_frequency_rad_s=relative,
        delta_eff_axis_rad_s=deff.copy(),
        pump_state=rho,
        chi_matrix=chi,
        diagnostics=diagnostics,
        provenance=pump_weak_response_reference_provenance(),
    )


def pump_only_weak_response_noncollinear_reference(
        Op_A, Op_B, lab_delta_axis_rad_s, one_photon_detuning_rad_s, *,
        T, pump_k_rad_m, probe_k_axis_rad_m, crossing_angle_rad,
        lab_optical_beat_axis_rad_s=None, branch=-1, atom=None,
        analysis_frequency_axis_rad_s=(0.0,), quadrature_order=24,
        cutoff_sigma=5.0):
    """Two-dimensional Maxwell reference with separated lab/atomic frequencies.

    The standard minus branch is evaluated in the certified static pump frame.
    The laboratory beat is fixed for every velocity class, while the atomic
    weak-field frequency acquires the Raman shift returned by
    :func:`gabes.doppler.noncollinear_raman_shift_rad_s`.  Passing
    ``delta_eff`` through :func:`chi_matrix_table` would be wrong because that
    finite-seed API also reconstructs the laboratory beat from its delta input.

    This is an opt-in slow reference.  It does not replace the production
    one-dimensional finite-seed gain path.
    """
    if branch != -1:
        _pump_nambu_operators(branch)  # raises the certified frame limitation
    T = float(T)
    one_photon = float(one_photon_detuning_rad_s)
    pump_k = float(pump_k_rad_m)
    theta = float(crossing_angle_rad)
    delta_lab = np.atleast_1d(np.asarray(lab_delta_axis_rad_s, dtype=float))
    analysis = np.atleast_1d(
        np.asarray(analysis_frequency_axis_rad_s, dtype=float))
    if delta_lab.ndim != 1 or delta_lab.size < 1 or not np.isfinite(delta_lab).all():
        raise ValueError(
            "lab_delta_axis_rad_s must be a finite non-empty one-dimensional array")
    if analysis.ndim != 1 or analysis.size < 1 or not np.isfinite(analysis).all():
        raise ValueError(
            "analysis_frequency_axis_rad_s must be a finite non-empty 1-D array")
    if not all(np.isfinite(value) for value in (T, one_photon, pump_k, theta)):
        raise ValueError("T, one-photon detuning, pump k, and angle must be finite")
    if T <= 0.0 or pump_k <= 0.0:
        raise ValueError("T and pump_k_rad_m must be positive")

    probe_k = np.asarray(probe_k_axis_rad_m, dtype=float)
    try:
        probe_k = np.broadcast_to(probe_k, delta_lab.shape).astype(float, copy=True)
    except ValueError as exc:
        raise ValueError(
            "probe_k_axis_rad_m must be scalar or match lab_delta_axis_rad_s") from exc
    if not np.isfinite(probe_k).all() or np.any(probe_k <= 0.0):
        raise ValueError("probe wave numbers must be finite and positive")

    if lab_optical_beat_axis_rad_s is None:
        lab_beat = seeded_sideband_beat(delta_lab, branch)
    else:
        try:
            lab_beat = np.broadcast_to(
                np.asarray(lab_optical_beat_axis_rad_s, dtype=float),
                delta_lab.shape,
            ).astype(float, copy=True)
        except ValueError as exc:
            raise ValueError(
                "lab_optical_beat_axis_rad_s must be scalar or match the delta axis") from exc
    if not np.isfinite(lab_beat).all() or np.any(lab_beat <= 0.0):
        raise ValueError("laboratory optical beat values must be finite and positive")
    if atom is None:
        atom = collisional_atom(T)

    velocity, weights = doppler.maxwell_legendre_grid(
        T, order=quadrature_order, cutoff_sigma=cutoff_sigma)
    vx = velocity.copy()
    vz = velocity.copy()
    wx = weights.copy()
    wz = weights.copy()
    deff_z = one_photon - pump_k * vz
    (_deff, _L0, L_z, rho_z, rho_vec_z, sources_z, readouts) = (
        _pump_reference_system(Op_A, Op_B, deff_z, atom, branch))

    nz = vz.size
    nx = vx.size
    M = N_LEVELS * N_LEVELS
    L_pairs = np.broadcast_to(L_z[:, None, :, :], (nz, nx, M, M))
    source_pairs = np.broadcast_to(
        sources_z[:, None, :, :], (nz, nx, M, sources_z.shape[-1]))
    pair_weights = wz[:, None] * wx[None, :]
    shift = np.empty((delta_lab.size, nz, nx), dtype=float)
    chi_average = np.empty(
        (analysis.size, delta_lab.size, 2, 2), dtype=complex)
    diagonal = np.arange(N_LEVELS) * (N_LEVELS + 1)
    max_response_trace = 0.0
    max_response_residual = 0.0
    eye = np.eye(M, dtype=complex)

    for delta_index, probe_wave_number in enumerate(probe_k):
        shift[delta_index] = doppler.noncollinear_raman_shift_rad_s(
            vx[None, :], vz[:, None], pump_k, probe_wave_number, theta)
        for analysis_index, analysis_frequency in enumerate(analysis):
            relative = (
                -lab_beat[delta_index]
                + analysis_frequency
                + shift[delta_index]
            )
            response = trace_zero_liouvillian_response(
                L_pairs, source_pairs, relative, N_LEVELS)
            class_chi = np.einsum("om,zxmi->zxoi", readouts, response)
            chi_average[analysis_index, delta_index] = np.einsum(
                "zx,zxoi->oi", pair_weights, class_chi)

            traces = np.sum(response[..., diagonal, :], axis=-2)
            max_response_trace = max(
                max_response_trace,
                float(np.max(np.abs(traces), initial=0.0)))
            residual = np.einsum("zxij,zxjk->zxik", L_pairs, response)
            residual += 1j * relative[..., None, None] * response
            residual += source_pairs
            response_action = np.einsum(
                "zxij,zxjk->zxik", L_pairs, response)
            scale = (
                np.linalg.norm(response_action, axis=-2)
                + np.abs(relative)[..., None] * np.linalg.norm(response, axis=-2)
                + np.linalg.norm(source_pairs, axis=-2)
            )
            normalized = np.linalg.norm(residual, axis=-2) / np.maximum(
                scale, np.finfo(float).tiny)
            max_response_residual = max(
                max_response_residual,
                float(np.max(normalized, initial=0.0)))

    rho_vector = rho_z.reshape(nz, M)
    pump_residual = np.einsum("zij,zj->zi", L_z, rho_vector)
    pump_scale = np.linalg.norm(L_z, axis=(1, 2)) * np.linalg.norm(
        rho_vector, axis=1)
    pump_normalized = np.linalg.norm(pump_residual, axis=1) / np.maximum(
        pump_scale, np.finfo(float).tiny)
    pump_traces = np.trace(rho_z, axis1=-2, axis2=-1)
    pump_hermiticity = rho_z - rho_z.conj().swapaxes(-1, -2)

    quadrature_rms = []
    analytic_rms = []
    for delta_index, probe_wave_number in enumerate(probe_k):
        mean = float(np.sum(pair_weights * shift[delta_index]))
        variance = float(np.sum(
            pair_weights * (shift[delta_index] - mean)**2))
        quadrature_rms.append(np.sqrt(max(variance, 0.0)))
        analytic_rms.append(doppler.noncollinear_raman_rms_budget(
            T, pump_k, probe_wave_number, theta)["total_rms_rad_s"])

    diagnostics = {
        "finite": bool(np.isfinite(rho_z).all() and np.isfinite(chi_average).all()),
        "max_pump_normalized_residual": float(
            np.max(pump_normalized, initial=0.0)),
        "max_pump_trace_error": float(
            np.max(np.abs(pump_traces - 1.0), initial=0.0)),
        "max_pump_hermiticity_error": float(
            np.max(np.abs(pump_hermiticity), initial=0.0)),
        "minimum_pump_eigenvalue": float(np.min(np.linalg.eigvalsh(rho_z))),
        "max_response_normalized_residual": max_response_residual,
        "max_response_trace_error": max_response_trace,
        "quadrature_order_per_axis": int(vx.size),
        "velocity_pairs": int(nx * nz),
        "cutoff_sigma": float(cutoff_sigma),
        "quadrature_raman_rms_rad_s": np.asarray(quadrature_rms),
        "analytic_raman_rms_rad_s": np.asarray(analytic_rms),
        "lab_beat_velocity_invariant": True,
        "two_photon_points": int(delta_lab.size),
        "analysis_frequency_points": int(analysis.size),
    }
    return NoncollinearPumpWeakResponseReference(
        branch=branch,
        lab_delta_axis_rad_s=delta_lab.copy(),
        lab_optical_beat_axis_rad_s=lab_beat.copy(),
        analysis_frequency_axis_rad_s=analysis.copy(),
        one_photon_detuning_rad_s=one_photon,
        pump_k_rad_m=pump_k,
        probe_k_axis_rad_m=probe_k,
        crossing_angle_rad=theta,
        vx_axis_m_s=vx,
        vz_axis_m_s=vz,
        vx_weights=wx,
        vz_weights=wz,
        raman_shift_grid_rad_s=shift,
        one_photon_effective_axis_rad_s=deff_z,
        pump_state=rho_z,
        chi_matrix=chi_average,
        diagnostics=diagnostics,
        provenance=noncollinear_doppler_reference_provenance(),
    )


def pump_only_pole_residue_reference(
        Op_A, Op_B, delta, Delta_eff, *, branch=-1, atom=None,
        analysis_frequency_rad_s=0.0, T=T_CELL):
    """One-velocity pole/residue audit of the pump-only Nambu response.

    Omitting ``atom`` selects the production collisional/transit-reset model at
    ``T``; an explicit atom keeps archived or adversarial audits reproducible.
    """
    if atom is None:
        atom = collisional_atom(T)
    (deff, _L0, L_batch, rho, rho_vec, sources, readouts) = (
        _pump_reference_system(
            Op_A, Op_B, np.asarray([Delta_eff], dtype=float), atom, branch))
    carrier_relative = -OMEGA_HF + float(delta)
    omega = carrier_relative + float(analysis_frequency_rad_s)
    result = liouvillian_pole_residue_response(
        L_batch[0], sources[0], readouts, omega)
    result.update({
        "branch": branch,
        "delta_rad_s": float(delta),
        "delta_eff_rad_s": float(deff[0]),
        "analysis_frequency_rad_s": float(analysis_frequency_rad_s),
        "carrier_relative_frequency_rad_s": carrier_relative,
        "resolvent_frequency_rad_s": omega,
        "analysis_pole_centers_rad_s": (
            -carrier_relative - np.imag(result["eigenvalues"])),
        "pump_state": rho[0],
    })
    return result


def static_hamiltonian_at_Deff_zero(Op_A, Op_B, Os, delta, branch):
    """H₀ with Δ_eff = 0, so the only v / Δ_eff dependence is added later."""
    H0 = np.zeros((N_LEVELS, N_LEVELS), dtype=complex)
    H0[G2, G2] = delta
    H0[E2, E2] = -OMEGA_EXCITED_HF
    H0[E3, E3] = 0.0
    if branch == -1:
        _add_static_drive(H0, G1, Op_A)
        _add_static_drive(H0, G2, Os)
        return H0
    if branch == +1:
        _add_static_drive(H0, G1, Os)
        _add_static_drive(H0, G2, Op_B)
        return H0
    raise ValueError(f"branch must be one of {BRANCHES}, got {branch}")


def sideband_hamiltonian(Op_A, Op_B, Oc, branch):
    Hp = np.zeros((N_LEVELS, N_LEVELS), dtype=complex)
    if branch == -1:
        _add_sideband_drive(Hp, G1, Oc)
        _add_sideband_drive(Hp, G2, Op_B)
        return Hp
    if branch == +1:
        _add_sideband_drive(Hp, G1, Op_A)
        _add_sideband_drive(Hp, G2, Oc)
        return Hp
    raise ValueError(f"branch must be one of {BRANCHES}, got {branch}")


def sideband_template(Op_A, Op_B, Oc, branch):
    Hp = sideband_hamiltonian(Op_A, Op_B, Oc, branch)
    Cp = comm_super(Hp)
    Cm = comm_super(Hp.conj().T)
    return Cp, Cm


# =========================================================
# χ-matrix table for one fixed atom/dissipator on (δ, Δ_eff)
# =========================================================
def _coherence_weights(ground):
    """w such that Σ_k w[k]·vec(ρ)[k] = Σ_e scale[g,e]·ρ[e,g] (vec row-major)."""
    w = np.zeros(N_LEVELS * N_LEVELS, dtype=complex)
    for e in EXCITED_STATES:
        w[e * N_LEVELS + ground] = TRANSITION_DIPOLE_SCALE[ground, e]
    return w


def chi_matrix_table(Op_A, Op_B, Os_ref, Oc_ref, delta_axis, Delta_eff_axis, branch,
                     atom=ATOM, n_f=1):
    """
    Two finite-Floquet solves per probe-detuning point to extract
    (chi_ss, chi_cs, chi_sc, chi_cc) on a 2-D (delta, Delta_eff) grid.
    ``n_f`` selects the retained harmonics ``-n_f,...,+n_f``.

    `atom` carries the level scheme and its dissipator; pass a temperature-
    dependent model (collisional ground/optical dephasing) to fold density-
    dependent decoherence into the solve. Defaults to the natural-linewidth
    module `ATOM` so the kernel/NumPy regression callers are unchanged.

    With numba available the whole grid runs in one fused compiled continued-
    fraction kernel (`kernels.floquet_chi_grid`); the NumPy delta-loop below is
    the fallback and reference implementation.  Both paths support ``n_f>=1``.
    """
    n_f = int(n_f)
    if n_f < 1:
        raise ValueError("n_f must be at least 1")
    probe_ground = G2 if branch == -1 else G1
    conj_ground = G1 if branch == -1 else G2
    n_d = delta_axis.size
    n_de = Delta_eff_axis.size

    Cp_no_c, Cm_no_c = sideband_template(Op_A, Op_B, 0.0, branch)      # solve 1
    Cp_c, Cm_c = sideband_template(Op_A, Op_B, Oc_ref, branch)        # solve 2

    if kernels.available():
        # H₀(δ) is affine in δ (only H₀[G2,G2] = δ), so L₀(δ) = L₀(0) + δ·C_δ.
        E_g2 = np.zeros((N_LEVELS, N_LEVELS), dtype=complex)
        E_g2[G2, G2] = 1.0
        C_delta = comm_super(E_g2)
        w_probe = _coherence_weights(probe_ground)
        w_conj = _coherence_weights(conj_ground)
        delta_axis = np.ascontiguousarray(delta_axis, dtype=float)
        deff_axis = np.ascontiguousarray(Delta_eff_axis, dtype=float)

        L0_1 = build_liouvillian(
            static_hamiltonian_at_Deff_zero(Op_A, Op_B, Os_ref, 0.0, branch), atom)
        probe_a, conj_a = kernels.floquet_chi_grid(
            L0_1, C_delta, atom.S_v, Cp_no_c, Cm_no_c, delta_axis, deff_axis,
            OMEGA_HF, float(branch), w_probe, w_conj, N_LEVELS, n_f)

        L0_2 = build_liouvillian(
            static_hamiltonian_at_Deff_zero(Op_A, Op_B, 0.0, 0.0, branch), atom)
        probe_b, conj_b = kernels.floquet_chi_grid(
            L0_2, C_delta, atom.S_v, Cp_c, Cm_c, delta_axis, deff_axis,
            OMEGA_HF, float(branch), w_probe, w_conj, N_LEVELS, n_f)

        return probe_a / Os_ref, conj_a / Os_ref, probe_b / Oc_ref, conj_b / Oc_ref

    chi_ss = np.zeros((n_d, n_de), dtype=complex)
    chi_cs = np.zeros((n_d, n_de), dtype=complex)
    chi_sc = np.zeros((n_d, n_de), dtype=complex)
    chi_cc = np.zeros((n_d, n_de), dtype=complex)

    for i, delta in enumerate(delta_axis):
        Omega_beat = seeded_sideband_beat(delta, branch)

        # ---- Solve 1: probe drive only ----
        H0_1 = static_hamiltonian_at_Deff_zero(Op_A, Op_B, Os_ref, delta, branch)
        L0_1 = build_liouvillian(H0_1, atom)
        rho0_a, rhop_a = floquet_solve_truncated(
            L0_1, Cp_no_c, Cm_no_c, Omega_beat, Delta_eff_axis,
            atom.S_v, N_LEVELS, n_f)
        probe_a = _polarization_coherence(rho0_a, probe_ground)
        conj_a = _polarization_coherence(rhop_a, conj_ground)
        chi_ss[i] = probe_a / Os_ref
        chi_cs[i] = conj_a / Os_ref

        # ---- Solve 2: conjugate seed only ----
        H0_2 = static_hamiltonian_at_Deff_zero(Op_A, Op_B, 0.0, delta, branch)
        L0_2 = build_liouvillian(H0_2, atom)
        rho0_b, rhop_b = floquet_solve_truncated(
            L0_2, Cp_c, Cm_c, Omega_beat, Delta_eff_axis,
            atom.S_v, N_LEVELS, n_f)
        probe_b = _polarization_coherence(rho0_b, probe_ground)
        conj_b = _polarization_coherence(rhop_b, conj_ground)
        chi_sc[i] = probe_b / Oc_ref
        chi_cc[i] = conj_b / Oc_ref

    return chi_ss, chi_cs, chi_sc, chi_cc


# =========================================================
# Probe-detuning axis
# =========================================================
def branch_center_GHz(Delta_GHz, branch):
    if branch not in BRANCHES:
        raise ValueError(f"branch must be one of {BRANCHES}, got {branch}")
    return Delta_GHz + branch * constants.NU_HF / 1e9


def probe_scan_axis_GHz(Delta_GHz, coarse_points=None, fine_points=None,
                        window_mhz=None, scan_min=None, scan_max=None,
                        branches=BRANCHES):
    """Probe-frequency scan axis [GHz], referenced to the 85Rb F=2 → F'=3 line."""
    coarse_points = SCAN_COARSE_POINTS if coarse_points is None else coarse_points
    fine_points = SCAN_FINE_POINTS if fine_points is None else fine_points
    window_mhz = RESONANCE_WINDOW_MHZ if window_mhz is None else window_mhz
    scan_min = SCAN_MIN_GHZ if scan_min is None else scan_min
    scan_max = SCAN_MAX_GHZ if scan_max is None else scan_max

    coarse = np.linspace(scan_min, scan_max, coarse_points)
    half_window = window_mhz * 1e-3
    parts = [coarse]
    for branch in branches:
        center = branch_center_GHz(Delta_GHz, branch)
        fmin = max(scan_min, center - half_window)
        fmax = min(scan_max, center + half_window)
        if fmin < fmax:
            parts.append(np.linspace(fmin, fmax, fine_points))
    axis = np.unique(np.concatenate(parts))
    exclusion_GHz = PUMP_OVERLAP_EXCLUSION_MHZ * 1e-3
    return axis[np.abs(axis - Delta_GHz) > exclusion_GHz]


def two_photon_detuning_from_probe_scan(probe_GHz, Delta_GHz, branch):
    delta_Hz = (probe_GHz - branch_center_GHz(Delta_GHz, branch)) * 1e9
    return 2 * np.pi * delta_Hz


def _optical_k_from_offset(offset_rad_s):
    return (constants.OMEGA_D1 + np.asarray(offset_rad_s, dtype=float)) / constants.C_LIGHT


def seeded_option_a_wavenumbers(D_GHz, probe_axis_GHz):
    """Bare frequency-specific pump, probe, and conjugate wave numbers."""
    pump_offset = 2 * np.pi * float(D_GHz) * 1e9
    seed_offset = 2 * np.pi * np.asarray(probe_axis_GHz, dtype=float) * 1e9
    conj_offset = 2.0 * pump_offset - seed_offset
    k_pump = _optical_k_from_offset(pump_offset)
    k_seed = _optical_k_from_offset(seed_offset)
    k_conj = _optical_k_from_offset(conj_offset)
    return k_pump, k_seed, k_conj


def seeded_phase_mismatch_z(D_GHz, probe_axis_GHz,
                            angle_deg=SEEDED_PHASE_ANGLE_DEG):
    """Option-A vacuum/geometric mismatch ``2k_p-(k_s+k_c)cos(theta)``.

    Susceptibility, including its real dispersive part, stays in the diagonal
    Maxwell response. Refractive indices must not be folded into this mismatch.
    """
    theta = math.radians(float(angle_deg))
    k_pump, k_seed, k_conj = seeded_option_a_wavenumbers(
        D_GHz, probe_axis_GHz)
    return 2.0 * k_pump - (k_seed + k_conj) * math.cos(theta)


def _uniform_segment_profile_and_probe_od(
        chi_bar, N_atoms, line_strength, nseg, L=L_CELL):
    """Return a neutral coupling profile plus a probe-OD diagnostic.

    ``chi_bar`` already supplies the probe's local diagonal attenuation in the
    Maxwell matrix.  It cannot also be used as a pump-depletion law for the
    off-diagonal FWM coupling.  Until a pump-frequency propagation solve is
    available, keep that coupling uniform and report the probe OD separately.
    """
    chi = observables.chi_phys(chi_bar, N_atoms, line_strength=line_strength)
    alpha = np.maximum(K_VEC * np.imag(chi), 0.0)
    od = float(np.nanmedian(alpha) * L) if alpha.size else 0.0
    od = max(od, 0.0)
    return np.ones(int(nseg), dtype=float), od


def _arm_linear_od(beam_GHz, T, L=L_CELL, ground_F=2):
    """Independent linear-OD diagnostic at a twin-beam arm's optical frequency.

    It is not applied as a second arm efficiency because the diagonal Maxwell
    response already attenuates the propagated fields.  We sum the requested
    ground-manifold components of the independently validated absorption model.
    Δ=0 is 85Rb F=2→F'=3, hence the reference-line shift.
    """
    from . import absorption
    ref_hz = hyperfine.LINE_SHIFT_HZ[(2, 3)]
    scan = 2 * np.pi * (np.asarray(beam_GHz, dtype=float) * 1e9 + ref_hz)
    params = {"temp_c": float(T) - 273.15, "line_strength": 1.0, "doppler": "on"}
    _alpha_tot, comps, _info = absorption._hyperfine_alpha(scan, params)
    alpha = np.zeros_like(scan, dtype=float)
    for (Fg, _Fe), a in comps.items():
        if Fg == ground_F:
            alpha = alpha + np.asarray(a, dtype=float)
    return np.clip(np.maximum(alpha, 0.0) * L, 0.0, 5.0)


def _pump_scatter_noise(D_GHz, T, L, kappa):
    """Excess intensity-difference noise (normalized to SQL) from pump light
    scattered into the detection mode.

    Anchored to the validated hyperfine absorption path (`schemes.absorption.
    _hyperfine_alpha`, <0.1 % vs lab): the pump's own linear OD at its detuning
    fixes how much pump power is scattered by the vapor; `kappa` is the residual
    fraction reaching the difference photocurrent (PBS extinction + spatial
    filtering leakage). Small far from resonance (Sim's +0.9 GHz window), large
    when the pump sits near a populous manifold (the −2.1 GHz / F=3 region).
    """
    from . import absorption
    # `_hyperfine_alpha`'s scan axis is referenced to the 87Rb F=2→F'=2 marker
    # (hyperfine.LINE_SHIFT_HZ), whereas the FWM Δ is referenced to 85Rb F=2→F'=3.
    # Shift by that line's position so Δ=0 maps to the F=2→F'=3 line center.
    ref_hz = hyperfine.LINE_SHIFT_HZ[(2, 3)]
    scan = np.array([2 * np.pi * (float(D_GHz) * 1e9 + ref_hz)])
    params = {"temp_c": float(T) - 273.15, "line_strength": 1.0, "doppler": "on"}
    alpha, _, _ = absorption._hyperfine_alpha(scan, params)
    od_pump = float(np.clip(float(alpha[0]) * L, 0.0, 50.0))
    return float(kappa) * (1.0 - math.exp(-od_pump)), od_pump


def _gaussian_overlap_profile(nseg, L, w_pump, w_probe, angle_deg):
    """Amplitude overlap of crossed Gaussian beams along the cell."""
    theta = math.radians(float(angle_deg))
    z = ((np.arange(nseg, dtype=float) + 0.5) / nseg - 0.5) * L
    separation = np.abs(z * math.tan(theta))
    waist_sq = max(float(w_pump) ** 2 + float(w_probe) ** 2, 1e-30)
    profile = np.exp(-(separation ** 2) / waist_sq)
    return profile / max(float(np.nanmax(profile)), 1e-30)


def _ultra_segmented_gain(chi_ss_avg, chi_sc_avg, chi_cs_avg, chi_cc_avg,
                          k_probe, k_conj, L, N_atoms, line_strength,
                          delta_k_z, segment_profile, spatial_profile,
                          P_pump, P_seed, conjugate_power_ratio=1.0):
    """Segmented propagation with approximate dynamic pump depletion.

    Diagonal susceptibility already attenuates each field.  ``segment_profile``
    may scale only independently known pump/crossing overlap; a probe-absorption
    estimate must not be recycled as an off-diagonal pump profile or as a second
    post-source loss.
    """
    nseg = int(segment_profile.size)
    dz = L / max(nseg, 1)
    M = observables._gain_matrix_from_chi(
        chi_ss_avg, chi_sc_avg, chi_cs_avg, chi_cc_avg,
        k_probe, k_conj, N_atoms, constants.DIPOLE_D1, line_strength,
        delta_k_z=delta_k_z)
    n = M.shape[0]
    amp = np.zeros((n, 2), dtype=complex)
    amp[:, 0] = math.sqrt(max(float(P_seed), 1e-30))
    T_total = np.broadcast_to(np.eye(2, dtype=complex), (n, 2, 2)).copy()
    pump_remaining = np.full(n, max(float(P_pump), 1e-30), dtype=float)

    for coupling_scale, spatial_scale in zip(segment_profile, spatial_profile):
        pump_scale = np.sqrt(np.clip(pump_remaining / max(float(P_pump), 1e-30),
                                     0.0, 1.0))
        Mz = M.copy()
        scale = float(coupling_scale) * float(spatial_scale) * pump_scale
        Mz[:, 0, 1] *= scale
        Mz[:, 1, 0] *= scale
        Tseg = observables.matrix_exp_2x2(Mz, dz)
        amp = np.einsum("nij,nj->ni", Tseg, amp)
        T_total = Tseg @ T_total
        seed_added = np.maximum(np.abs(amp[:, 0]) ** 2 - float(P_seed), 0.0)
        conj_power = np.maximum(
            np.abs(amp[:, 1]) ** 2 * float(conjugate_power_ratio), 0.0)
        pump_remaining = np.maximum(float(P_pump) - seed_added - conj_power, 0.0)

    G_s = np.abs(amp[:, 0]) ** 2 / max(float(P_seed), 1e-30)
    G_c = np.abs(amp[:, 1]) ** 2 / max(float(P_seed), 1e-30)
    return G_s, G_c, T_total, pump_remaining


def thermal_transit_reset_superoperator(rate, populations=None):
    """Trace-preserving ``rate·(rho_th Tr(rho) - rho)`` superoperator.

    ``rho_th`` defaults to an unpolarized thermal mixture of the representative
    F=2/F=3 ground manifolds with degeneracy weights 5/12 and 7/12.  The reset
    models atoms leaving the optical mode and fresh thermal atoms entering it;
    unlike pure ground-coherence dephasing, it also replenishes populations.
    """
    rate = float(rate)
    if rate < 0.0:
        raise ValueError("transit reset rate must be non-negative")
    if populations is None:
        populations = (
            hyperfine.GROUND_POP[GROUND_F[G1]],
            hyperfine.GROUND_POP[GROUND_F[G2]],
            0.0,
            0.0,
        )
    populations = np.asarray(populations, dtype=float)
    if populations.shape != (N_LEVELS,):
        raise ValueError(f"populations must have shape ({N_LEVELS},)")
    if np.any(populations < 0.0) or not np.isfinite(populations).all():
        raise ValueError("transit target populations must be finite and non-negative")
    total = float(np.sum(populations))
    if total <= 0.0:
        raise ValueError("transit target populations must have positive trace")
    rho_target = np.diag(populations / total).reshape(-1)
    trace_row = np.zeros(N_LEVELS * N_LEVELS, dtype=float)
    trace_row[np.arange(N_LEVELS) * (N_LEVELS + 1)] = 1.0
    return rate * (np.outer(rho_target, trace_row)
                   - np.eye(N_LEVELS * N_LEVELS, dtype=float))


def collisional_atom(T, density=None, *, transit_rate=None):
    """Per-temperature double-Λ atom with collisions and thermal transit reset.

      • optical (g,e) coherence — Rb self-broadening, added collisional
        half-width (Γ_eff − Γ)/2 via the AutoOD-validated
        `hyperfine.self_broadened_gamma`;
      • ground Raman (g₁,g₂) coherence — Rb–Rb pure dephasing plus a
        trace-preserving thermal transit reset.  The reset supplies the inherited
        transit/residual floor exactly once while replenishing populations.

    Because χ̄ now depends on temperature through this dissipator, any cached-χ̄
    scanner must rebuild its table per temperature with this atom (a single
    all-temperature table is no longer valid). Shared by `compute_spectrum` and
    the analysis scanners so the optimum-finder sees the same physics.
    """
    if density is None:
        density = hyperfine.number_density(T)
    transit_rate = constants.GAMMA_GG if transit_rate is None else float(transit_rate)
    if transit_rate < 0.0:
        raise ValueError("transit_rate must be non-negative")
    gamma_collision = constants.ground_coherence_dephasing(
        T, density, floor=0.0)
    gamma_opt = 0.5 * (hyperfine.self_broadened_gamma(density) - constants.GAMMA)
    atom = atoms.double_lambda_rb85(
        gamma_gg=gamma_collision, gamma_opt=max(gamma_opt, 0.0))
    atom.lindblad = atom.lindblad + thermal_transit_reset_superoperator(transit_rate)
    atom.transit_reset_rate = transit_rate
    atom.ground_collision_dephasing_rate = gamma_collision
    atom.thermal_reset_populations = (
        hyperfine.GROUND_POP[GROUND_F[G1]],
        hyperfine.GROUND_POP[GROUND_F[G2]],
        0.0,
        0.0,
    )
    return atom


def _single_branch(branch, branches):
    """Resolve one physical FWM channel; do not merge distinct Raman branches."""
    if branches is None:
        if branch not in BRANCHES:
            raise ValueError(f"branch must be one of {BRANCHES}, got {branch}")
        return branch

    branches = tuple(branches)
    if len(branches) != 1:
        raise ValueError(
            "FWM Raman branches are separate probe/conjugate mode pairs; "
            "compute them one at a time instead of summing susceptibilities."
        )
    only = branches[0]
    if only not in BRANCHES:
        raise ValueError(f"branch must be one of {BRANCHES}, got {only}")
    return only


# =========================================================
# High-level spectrum  (one call → gain + gain-referred diagnostic curves)
# =========================================================
def compute_spectrum(D_GHz, *,
                     T=T_CELL, P_pump=P_PUMP, P_probe=P_PROBE,
                     w_pump=W_PUMP, w_probe=W_PROBE, w_conj=None,
                     line_strength=None,
                     mode_overlap_penalty=1.0, polarization_penalty=1.0,
                     zeeman_participation_penalty=1.0,
                     loss_frac=LOSS_FRAC, qe=QE_DETECTOR,
                     detection_efficiency=None,
                     L=L_CELL,
                     coarse_points=None, fine_points=None, window_mhz=None,
                     scan_min=None, scan_max=None,
                     velocity_step=None, velocity_cutoff=None,
                     branch=DEFAULT_BRANCH, branches=None,
                     phase_detail=PHASE_LEGACY,
                     pump_probe_angle_deg=SEEDED_PHASE_ANGLE_DEG,
                     model_fidelity=None,
                     excess_noise_model=None,
                     floquet_order=SEEDED_FLOQUET_ORDER,
                     enforce_floquet_convergence=True,
                     transit_rate=constants.GAMMA_GG,
                     eom_residual_carrier_power=0.0,
                     eom_other_sidebands_power=0.0,
                     eom_seed_spectrum_provenance="not supplied",
                     eom_seed_spectrum_status="not supplied",
                     eom_seed_spectrum_application="unapplied"):
    """Mean-field seeded-FWM pipeline at Δ = 2π·D_GHz·1e9.

    ``excess_noise_model`` controls only an optional phenomenological pump-scatter
    term in the Ultra gain-referred diagnostic.  Arm attenuation is already in
    the Maxwell drift; distributed loss vacuum and atomic diffusion are not
    supplied.  EOM residual carrier/sideband powers are carried as explicit input
    provenance but remain unapplied until a measured coupling model exists.
    """
    branch = _single_branch(branch, branches)
    floquet_order = int(floquet_order)
    if floquet_order < 1:
        raise ValueError("floquet_order must be at least 1")
    if enforce_floquet_convergence and floquet_order < 3:
        raise ValueError(
            "reported scans require floquet_order >= 3; set "
            "enforce_floquet_convergence=False only for an explicitly historical fixture")
    if line_strength is None:
        line_strength = SEEDED_REFERENCE_RESIDUAL
    factors = SeededCouplingFactors(
        reference_residual=line_strength,
        mode_overlap_penalty=mode_overlap_penalty,
        polarization_penalty=polarization_penalty,
        zeeman_participation_penalty=zeeman_participation_penalty,
    )
    # The reference residual and the three lab-facing relative factors sit on
    # top of the first-principles macroscopic normalization. ``coupling_ls`` is
    # what enters every χ_phys / gain call below.
    coupling_ledger = physical_coupling_ledger(branch)
    coupling_norm = coupling_ledger["macroscopic_coupling_norm"]
    coupling_ls = factors.combined_residual * coupling_norm
    eta = (qe * (1.0 - loss_frac) if detection_efficiency is None
           else float(detection_efficiency))
    if not 0.0 <= eta <= 1.0:
        raise ValueError("detection efficiency must be between 0 and 1")

    eom_residual_carrier_power = float(eom_residual_carrier_power)
    eom_other_sidebands_power = float(eom_other_sidebands_power)
    if eom_residual_carrier_power < 0.0 or eom_other_sidebands_power < 0.0:
        raise ValueError("EOM component powers must be non-negative")
    eom_seed_spectrum_status = str(eom_seed_spectrum_status).strip().lower()
    eom_seed_spectrum_application = str(
        eom_seed_spectrum_application).strip().lower()
    if eom_seed_spectrum_status not in {"not supplied", "unsupported"}:
        raise ValueError(
            "EOM spectrum status cannot be promoted without a calibrated "
            "carrier/sideband transfer model")
    if eom_seed_spectrum_application != "unapplied":
        raise ValueError(
            "EOM residual components must remain unapplied until their atomic "
            "and detector transfer is calibrated")
    wanted_power = max(float(P_probe), 0.0)
    wanted_denominator = max(wanted_power, 1e-30)
    eom_seed_spectrum = {
        "wanted_sideband_power_w": wanted_power,
        "residual_carrier_power_w": eom_residual_carrier_power,
        "other_sidebands_power_w": eom_other_sidebands_power,
        "residual_carrier_to_wanted_ratio": (
            eom_residual_carrier_power / wanted_denominator),
        "other_sidebands_to_wanted_ratio": (
            eom_other_sidebands_power / wanted_denominator),
        "provenance": str(eom_seed_spectrum_provenance),
        "status": eom_seed_spectrum_status,
        "application": eom_seed_spectrum_application,
        "physics_effect": (
            "unapplied: carrier/sideband coupling, phase, polarization, and "
            "technical-noise transfer are not calibrated"),
    }

    Op_A = rabi_freq(P_pump, w_pump)
    Op_B = Op_A
    Os = rabi_freq(P_probe, w_probe)
    Os_ref = Os
    Oc_ref = Os                              # χ̄ is independent of |Ω_ref|

    # Pure-85Rb CRC vapor density, consistent with the AutoOD-validated
    # absorption path (`hyperfine.number_density`). The other (natural-abundance,
    # Steck) `atoms.rb85_density` understated N for the enriched FWM cell.
    N_atoms = hyperfine.number_density(T)
    Delta = 2 * np.pi * D_GHz * 1e9

    # Density-dependent collisions plus trace-preserving transit replenishment.
    atom_T = collisional_atom(T, N_atoms, transit_rate=transit_rate)
    transit_rate = atom_T.transit_reset_rate

    probe_axis_GHz = probe_scan_axis_GHz(
        D_GHz, coarse_points, fine_points, window_mhz, scan_min, scan_max,
        branches=(branch,))
    k_pump_vac, k_probe_vac, k_conj_vac = seeded_option_a_wavenumbers(
        D_GHz, probe_axis_GHz)
    omega_probe = constants.C_LIGHT * k_probe_vac
    omega_conj = constants.C_LIGHT * k_conj_vac
    area_probe = float(observables.gaussian_mode_area(w_probe))
    if w_conj is None:
        w_conj_effective = float(w_probe)
        canonical_mode_status = (
            "conditional: conjugate collection waist not supplied; matched to probe")
    else:
        w_conj_effective = float(w_conj)
        canonical_mode_status = "explicit probe and conjugate collected-mode waists"
    area_conj = float(observables.gaussian_mode_area(w_conj_effective))
    velocity_step = VELOCITY_STEP_MPS if velocity_step is None else velocity_step
    velocity_cutoff = VELOCITY_CUTOFF_SIGMA if velocity_cutoff is None else velocity_cutoff
    v_grid, weights = doppler.velocity_grid(
        T, dv=velocity_step, cutoff_sigma=velocity_cutoff)
    Delta_eff_axis = doppler.build_Delta_eff_axis(Delta, Delta, v_grid)

    delta_axis = two_photon_detuning_from_probe_scan(probe_axis_GHz, D_GHz, branch)
    # All four susceptibility tables share the same Δ_eff interpolation geometry.
    # Build it once for this solve instead of repeating the index/fraction work.
    idx_lo, frac = doppler.interpolation_weights(
        Delta_eff_axis, Delta, v_grid)

    def averaged_response_at_order(order):
        tables = chi_matrix_table(
            Op_A, Op_B, Os_ref, Oc_ref, delta_axis, Delta_eff_axis,
            branch, atom=atom_T, n_f=order)
        return tuple(
            doppler.apply_doppler_average(table, idx_lo, frac, weights)
            for table in tables)

    chi_ss_avg, chi_cs_avg, chi_sc_avg, chi_cc_avg = (
        averaged_response_at_order(floquet_order))
    lower_order_response = None
    if enforce_floquet_convergence:
        lower_order_response = averaged_response_at_order(floquet_order - 1)

    phase_detail = (phase_detail or PHASE_LEGACY).lower()
    delta_k_z = None
    delta_k_z_vacuum = None
    k_probe_prop = np.full_like(probe_axis_GHz, K_VEC, dtype=float)
    k_conj_prop = np.full_like(probe_axis_GHz, K_VEC, dtype=float)
    propagation_segments = 1
    segment_profile = None
    segment_od = 0.0
    spatial_profile = None
    ultra_phase_iterations = 0
    ultra_phase_max_change = 0.0
    ultra_pump_remaining_min = float(P_pump)
    zeeman_status = "inactive"
    zeeman_correction = 1.0
    if phase_detail != PHASE_LEGACY:
        delta_k_z_vacuum = seeded_phase_mismatch_z(
            D_GHz, probe_axis_GHz, angle_deg=pump_probe_angle_deg)
        delta_k_z = delta_k_z_vacuum
        # Option A: frequency-specific bare k in the Maxwell matrix and vacuum /
        # geometric mismatch only. Re(chi) remains in the diagonal response.
        k_probe_prop = k_probe_vac
        k_conj_prop = k_conj_vac
        if phase_detail == PHASE_FINE:
            propagation_segments = 16
            segment_profile, segment_od = _uniform_segment_profile_and_probe_od(
                chi_ss_avg, N_atoms, coupling_ls, propagation_segments, L=L)
        elif phase_detail == PHASE_ULTRA:
            propagation_segments = ULTRA_PROPAGATION_SEGMENTS
            segment_profile, segment_od = _uniform_segment_profile_and_probe_od(
                chi_ss_avg, N_atoms, coupling_ls, propagation_segments, L=L)
            spatial_profile = _gaussian_overlap_profile(
                propagation_segments, L, w_pump, w_probe, pump_probe_angle_deg)
            try:
                from .. import zeeman as _zeeman
                z_atom = _zeeman.rb85_d1_double_lambda_zeeman()
                zeeman_correction = float(
                    getattr(z_atom, "lumped_strength_correction", 1.0))
                zeeman_status = (
                    f"24-level CG-sum consistency check = {zeeman_correction:.4f} "
                    "(lumped 3·C_F² reproduced); diagnostic only — full Floquet "
                    "scan not run in Ultra v1")
            except Exception as exc:                  # pragma: no cover
                zeeman_status = f"unavailable: {exc}"

    if phase_detail == PHASE_ULTRA:
        if segment_profile is None:
            segment_profile = np.ones(propagation_segments, dtype=float)
        if spatial_profile is None:
            spatial_profile = np.ones(propagation_segments, dtype=float)
        # zeeman_correction is a CG-sum consistency diagnostic (≡1.0 by
        # construction), not an active correction, so it is no longer multiplied
        # into the coupling.
        G_s_effective, G_c_effective, _T_effective, pump_remaining = _ultra_segmented_gain(
            chi_ss_avg, chi_sc_avg, chi_cs_avg, chi_cc_avg,
            k_probe_prop, k_conj_prop, L, N_atoms,
            coupling_ls, delta_k_z,
            segment_profile, spatial_profile, P_pump, P_probe,
            conjugate_power_ratio=area_conj / area_probe)
        ultra_pump_remaining_min = float(np.nanmin(pump_remaining))
        # Canonical diagnostics require a linear map. Remove the seed-dependent
        # dynamic-depletion scale while retaining the fixed spatial/pump profile.
        _, _, T_small_signal = observables.gain_from_chi(
            chi_ss_avg, chi_sc_avg, chi_cs_avg, chi_cc_avg,
            k_probe_prop, k_conj_prop, L, N_atoms, line_strength=coupling_ls,
            delta_k_z=delta_k_z, propagation_segments=propagation_segments,
            segment_profile=segment_profile * spatial_profile)
    else:
        G_s_effective, G_c_effective, T_small_signal = observables.gain_from_chi(
            chi_ss_avg, chi_sc_avg, chi_cs_avg, chi_cc_avg,
            k_probe_prop, k_conj_prop, L, N_atoms, line_strength=coupling_ls,
            delta_k_z=delta_k_z, propagation_segments=propagation_segments,
            segment_profile=segment_profile)

    canonical = observables.canonical_transfer_diagnostics(
        T_small_signal, omega_probe, omega_conj, area_probe, area_conj)
    if lower_order_response is not None:
        low_ss, low_cs, low_sc, low_cc = lower_order_response
        _, _, T_lower_order = observables.gain_from_chi(
            low_ss, low_sc, low_cs, low_cc,
            k_probe_prop, k_conj_prop, L, N_atoms, line_strength=coupling_ls,
            delta_k_z=delta_k_z, propagation_segments=propagation_segments,
            segment_profile=(
                segment_profile * spatial_profile
                if phase_detail == PHASE_ULTRA else segment_profile),
        )
        canonical_lower = observables.canonical_transfer_diagnostics(
            T_lower_order, omega_probe, omega_conj, area_probe, area_conj)
        floquet_convergence = assess_floquet_scan_convergence(
            high_order=floquet_order,
            low_order=floquet_order - 1,
            high_response={
                "chi_ss": chi_ss_avg,
                "chi_cs": chi_cs_avg,
                "chi_sc": chi_sc_avg,
                "chi_cc": chi_cc_avg,
            },
            low_response={
                "chi_ss": low_ss,
                "chi_cs": low_cs,
                "chi_sc": low_sc,
                "chi_cc": low_cc,
            },
            high_transfer=canonical["transfer_canonical"],
            low_transfer=canonical_lower["transfer_canonical"],
            scan_axis=probe_axis_GHz,
        )
    else:
        floquet_convergence = {
            "status": "NOT_EVALUATED",
            "passed": False,
            "high_order": floquet_order,
            "comparison_order": None,
            "full_scan_points": int(probe_axis_GHz.size),
            "reason": (
                "adjacent-order full-scan comparison explicitly disabled; "
                "result is a historical/numerical fixture only"),
        }
    G_s_smallsignal = canonical["probe_power_gain"]
    G_c_smallsignal = canonical["conjugate_power_gain"]
    if phase_detail == PHASE_ULTRA:
        # The dynamic-depletion path propagates sqrt(power) amplitudes. Convert
        # its conjugate coefficient for an explicitly unequal collected area.
        G_s = G_s_effective
        G_c = G_c_effective * area_conj / area_probe
    else:
        G_s = G_s_smallsignal
        G_c = G_c_smallsignal

    # The propagation above is linear in the (undepleted) pump. At high density
    # it overshoots the energy the pump can supply, so apply a post-hoc Manley-Rowe
    # cap.  This preserves an energy ledger but is not a self-consistent depleted
    # three-field solve and therefore does not validate the absolute gain.
    G_s, G_c = observables.pump_depletion_saturation(G_s, G_c, P_pump, P_probe)
    hardened_noise = None
    pump_scatter_kappa_used = HARDENED_PUMP_SCATTER_KAPPA
    if phase_detail == PHASE_ULTRA:
        if excess_noise_model is False:
            # The diagonal Maxwell drift already contains the arm attenuation.
            # A second post-source transmission would double count it; explicit
            # distributed Langevin diffusion is a separate, unfinished model.
            gain_referred_noise_db = observables.gain_referred_noise_dB(
                G_s, G_c, eta)
        else:
            cfg = {} if excess_noise_model in (None, True) else dict(excess_noise_model)
            kappa = cfg.get("pump_scatter_kappa", HARDENED_PUMP_SCATTER_KAPPA)
            pump_scatter_kappa_used = float(kappa)
            # The F=2→F'=3 manifold is the linear absorption the lumped double-Λ
            # omits for the generated sidebands. Apply it at BOTH sidebands' own
            # optical frequencies: without the probe term the deep-blue-Δ region
            # parks the probe ON F=2→F'=3 with no penalty (a probe-side mirror of
            # the conjugate artifact).
            conj_ground = G1 if branch == -1 else G2
            manifold_F = GROUND_F[conj_ground]
            conj_GHz = 2.0 * float(D_GHz) - probe_axis_GHz
            od_conj_arr = _arm_linear_od(conj_GHz, T, L, ground_F=manifold_F)
            od_probe_lin = _arm_linear_od(probe_axis_GHz, T, L, ground_F=manifold_F)
            pump_scatter, od_pump = _pump_scatter_noise(D_GHz, T, L, kappa)
            # These OD curves remain useful diagnostics, but applying them as
            # post-source efficiencies would count absorption already present in
            # the diagonal Maxwell drift a second time.  Atomic vacuum/excess
            # noise needs the future distributed covariance solve instead.
            eta_s_arm = eta
            eta_c_arm = eta
            S_lin = observables.balanced_twin_beam_noise(
                G_s, G_c, eta_s_arm, eta_c_arm, reference_weight="dc",
                seed_excess_noise=pump_scatter)
            gain_referred_noise_db = 10.0 * np.log10(np.maximum(S_lin, 1e-30))
            hardened_noise = {
                "kappa": float(kappa),
                "pump_scatter_noise": float(pump_scatter),
                "od_pump": float(od_pump),
                "segment_od_probe": float(segment_od),
                "od_conj_max": float(np.nanmax(od_conj_arr)),
                "od_probe_lin_max": float(np.nanmax(od_probe_lin)),
                "od_conj_arr": od_conj_arr,
                "od_probe_lin_arr": od_probe_lin,
                "atomic_od_application": (
                    "diagnostic only; distributed Langevin covariance unavailable"),
            }
    else:
        gain_referred_noise_db = observables.gain_referred_noise_dB(G_s, G_c, eta)

    atomic_solver_provenance = seeded_atomic_solver_provenance(
        floquet_order=floquet_order, convergence=floquet_convergence)
    parameter_provenance = seeded_parameter_provenance(
        line_strength=line_strength,
        transit_rate=transit_rate,
        pump_scatter_kappa=pump_scatter_kappa_used,
    )
    claim_gate = seeded_validation_claim_gate(
        canonical_mode_status=canonical_mode_status,
        commutator_defect_max=canonical["commutator_defect_max"],
        floquet_convergence=floquet_convergence,
        eom_residual_carrier_power=eom_residual_carrier_power,
        eom_other_sidebands_power=eom_other_sidebands_power,
        eom_spectrum_status=eom_seed_spectrum_status,
        eom_spectrum_application=eom_seed_spectrum_application,
    )

    return {
        "D_GHz": D_GHz,
        "probe_axis_GHz": probe_axis_GHz,
        "G_s": G_s,
        "G_c": G_c,
        "gain_referred_noise_dB": gain_referred_noise_db,
        # Deliberately unavailable until the frequency-dependent microscopic
        # Langevin drift/diffusion and detector transfer are implemented.
        "physical_squeezing_dB": None,
        # Backward-compatible data key.  It is the same algebraic diagnostic,
        # not a physical squeezing prediction.
        "S_dB": gain_referred_noise_db,
        "G_s_smallsignal": G_s_smallsignal,
        "G_c_smallsignal": G_c_smallsignal,
        "G_s_smallsignal_peak": float(np.nanmax(G_s_smallsignal)),
        "G_c_smallsignal_peak": float(np.nanmax(G_c_smallsignal)),
        "T_field_small_signal": T_small_signal,
        "Q_photon_flux": canonical["Q"],
        "T_canonical_small_signal": canonical["transfer_canonical"],
        "conjugate_photon_flux_gain_smallsignal": canonical[
            "conjugate_photon_flux_gain"],
        "photon_flux_gap_smallsignal": canonical["photon_flux_gap"],
        "commutator_defect_smallsignal": canonical["commutator_defect"],
        "commutator_defect_max_smallsignal": canonical["commutator_defect_max"],
        "canonical_mode_status": canonical_mode_status,
        "canonical_mode_area_probe_m2": area_probe,
        "canonical_mode_area_conjugate_m2": area_conj,
        "omega_probe_rad_s": omega_probe,
        "omega_conjugate_rad_s": omega_conj,
        "displayed_gains_posthoc_saturated": True,
        "pump_depletion_cap": 1.0 + 0.5 * P_pump / max(P_probe, 1e-30),
        "phase_detail": phase_detail,
        "model_fidelity": model_fidelity or phase_detail,
        "pump_probe_angle_deg": pump_probe_angle_deg,
        "propagation_convention": (
            "legacy fixed-k" if phase_detail == PHASE_LEGACY
            else "Option A: bare frequency-specific k plus vacuum mismatch"),
        "k_probe_propagation_per_m": k_probe_prop,
        "k_conjugate_propagation_per_m": k_conj_prop,
        "delta_k_z": delta_k_z,
        "delta_k_z_vacuum": delta_k_z_vacuum,
        "phase_segments": propagation_segments,
        "segment_absorption_od": segment_od,
        "ultra_phase_iterations": ultra_phase_iterations,
        "ultra_phase_max_change": ultra_phase_max_change,
        "ultra_dynamic_depletion": phase_detail == PHASE_ULTRA,
        "ultra_in_cell_loss_noise": False,
        "squeezing_status": (
            "unavailable: gain-referred diagnostic only; microscopic distributed "
            "atomic diffusion is not implemented"),
        "validation_level": claim_gate["level"],
        "claim_gate": claim_gate,
        "atomic_solver_provenance": atomic_solver_provenance,
        "pump_weak_response_reference_provenance": (
            pump_weak_response_reference_provenance()),
        "noncollinear_doppler_reference_provenance": (
            noncollinear_doppler_reference_provenance()),
        "floquet_convergence": floquet_convergence,
        "floquet_order": floquet_order,
        "floquet_convergence_enforced": bool(enforce_floquet_convergence),
        "parameter_provenance": parameter_provenance,
        "eom_seed_spectrum": eom_seed_spectrum,
        "transit_reset_rate_rad_s": float(transit_rate),
        "ground_collision_dephasing_rate_rad_s": float(
            atom_T.ground_collision_dephasing_rate),
        "ultra_pump_remaining_min": ultra_pump_remaining_min,
        "ultra_spatial_overlap_min": (float(np.nanmin(spatial_profile))
                                      if spatial_profile is not None else 1.0),
        "zeeman_status": zeeman_status,
        "zeeman_correction": zeeman_correction,
        "eta": eta,
        "qe": qe,
        "w_pump_m": w_pump,
        "w_probe_m": w_probe,
        "w_conjugate_m": w_conj_effective,
        "cell_length_m": L,
        "branch": branch,
        "coupling_norm": coupling_norm,
        "coupling_normalization_ledger": coupling_ledger,
        # Keep the historical key equal to the direct API input.  Consumers can
        # opt into the factorized bookkeeping through the new explicit keys.
        "line_strength_residual": line_strength,
        "mode_overlap_penalty": factors.mode_overlap_penalty,
        "polarization_penalty": factors.polarization_penalty,
        "zeeman_participation_penalty": factors.zeeman_participation_penalty,
        "lab_coupling_factor": factors.lab_factor,
        "combined_line_strength_residual": factors.combined_residual,
        "effective_line_strength": coupling_ls,
        "N_atoms": N_atoms,
        "sigma_v": np.sqrt(constants.KB * T / constants.MASS_85RB),
        "n_velocity": v_grid.size,
        "Op_A_2pi_GHz": Op_A / (2 * np.pi) / 1e9,
        "Os_2pi_MHz": Os / (2 * np.pi) / 1e6,
        "raman_center_minus_GHz": branch_center_GHz(D_GHz, -1),
        "raman_center_plus_GHz": branch_center_GHz(D_GHz, +1),
        "hardened_noise": hardened_noise,
    }


def operating_point(spectrum, delta_mhz, branch=-1):
    """Read gains and the gain-referred diagnostic at a selected δ (MHz)."""
    probe_GHz = (spectrum["raman_center_minus_GHz"] if branch == -1
                 else spectrum["raman_center_plus_GHz"]) + delta_mhz * 1e-3
    x = spectrum["probe_axis_GHz"]
    diagnostic = spectrum.get("gain_referred_noise_dB")
    if diagnostic is None:
        diagnostic = spectrum["S_dB"]
    diagnostic_value = float(np.interp(probe_GHz, x, diagnostic))
    out = {
        "probe_GHz": probe_GHz,
        "G_s": float(np.interp(probe_GHz, x, spectrum["G_s"])),
        "G_c": float(np.interp(probe_GHz, x, spectrum["G_c"])),
        "gain_referred_noise_dB": diagnostic_value,
        "physical_squeezing_dB": None,
        "S_dB": diagnostic_value,  # compatibility alias; not physical squeezing
    }
    for source, target in (
        ("G_s_smallsignal", "G_s_smallsignal"),
        ("G_c_smallsignal", "G_c_smallsignal"),
        ("conjugate_photon_flux_gain_smallsignal", "G_c_flux_smallsignal"),
        ("photon_flux_gap_smallsignal", "photon_flux_gap_smallsignal"),
        ("commutator_defect_max_smallsignal", "commutator_defect_max_smallsignal"),
    ):
        if source in spectrum:
            out[target] = float(np.interp(probe_GHz, x, spectrum[source]))
    if spectrum.get("delta_k_z") is not None:
        out["delta_k_z"] = float(np.interp(probe_GHz, x, spectrum["delta_k_z"]))
    if spectrum.get("delta_k_z_vacuum") is not None:
        out["delta_k_z_vacuum"] = float(np.interp(
            probe_GHz, x, spectrum["delta_k_z_vacuum"]))
    return out


# =========================================================
# Scheme wrapper for the generic front-end
# =========================================================
WINDOW_GHZ = 0.55          # half-width of the focused probe window around (−) Raman
TPD_LIMIT_MHZ = 500.0
# Three tiers. The old coarse "Fast" (121 pts) was dropped and the rest renamed
# down — the real-basis + sideband-symmetry speedups (~1.6×) made it redundant.
# Times are re-estimated for the reference (deployment) environment the old
# labels used, scaled by the measured 1.6×: old ~6 s → ~4 s, old ~20 s → ~12 s.
# The per-tier solver settings are unchanged; only the labels moved.
FIDELITY_FAST = "Fast  (~4 s)"          # was "Balanced  (~6 s)" (181-pt) settings
FIDELITY_BALANCED = "Balanced  (~12 s)"  # was "High fidelity  (~20 s)" (301-pt)
FIDELITY_ULTRA = "Ultra  (slow)"
FIDELITY_LABELS = {
    FIDELITY_FAST: "Fast",
    FIDELITY_BALANCED: "Balanced",
    FIDELITY_ULTRA: "Ultra",
}
FWM_FIDELITY = {
    FIDELITY_FAST:     dict(coarse_points=181, velocity_step=4.0,
                            velocity_cutoff=3.0, phase_detail=PHASE_BALANCED),
    FIDELITY_BALANCED: dict(coarse_points=301, velocity_step=2.0,
                            velocity_cutoff=3.0, phase_detail=PHASE_FINE),
    FIDELITY_ULTRA:    dict(coarse_points=401, velocity_step=1.0,
                            velocity_cutoff=4.0, phase_detail=PHASE_ULTRA),
}
RESOLUTION = FWM_FIDELITY

# Old saved fidelity labels → current labels. Settings-preserving where a tier
# survived (old Balanced/High kept their solver settings under the new names);
# the removed coarse tier falls to the cheapest survivor.
_FIDELITY_LEGACY = {
    "Fast  (~3 s)":           FIDELITY_FAST,      # removed 121-pt tier → cheapest now
    "Balanced  (~6 s)":       FIDELITY_FAST,      # same 181-pt settings, renamed
    "High fidelity  (~20 s)": FIDELITY_BALANCED,  # same 301-pt settings, renamed
    "Fine  (~20 s)":          FIDELITY_BALANCED,  # legacy alias of the old High tier
}


def normalize_fidelity(value):
    """Map old saved labels onto current user-facing fidelity labels."""
    if value in FWM_FIDELITY:
        return value
    return _FIDELITY_LEGACY.get(value, FIDELITY_FAST)


def _detection_efficiency_from_params(params):
    """Return η while preserving legacy loss/QE overrides from SABES and APIs."""
    loss_pct = float(params.get("loss_pct", SEEDED_POST_CELL_LOSS_PCT))
    qe_pct = float(params.get("qe_pct", QE_DETECTOR * 100.0))
    legacy_eta = (qe_pct / 100.0) * (1.0 - loss_pct / 100.0)
    direct_pct = params.get("detection_eff_pct")
    if direct_pct is None:
        return legacy_eta

    legacy_inputs_changed = (
        not np.isclose(loss_pct, SEEDED_POST_CELL_LOSS_PCT, rtol=0.0, atol=1e-12)
        or not np.isclose(qe_pct, QE_DETECTOR * 100.0, rtol=0.0, atol=1e-12)
    )
    return legacy_eta if legacy_inputs_changed else float(direct_pct) / 100.0


class FWMScheme(Scheme):
    name = "fwm"
    cluster = "D — Wave mixing"
    title = "Four-wave mixing (Squeezing / Biphoton)"
    cache_version = "fwm-angular-doppler-reference-v8"
    defaults_version = "fwm-ui-simplification-v1"
    cache_observables = True
    supports_headless_observables = True
    # Dev note: Gain diagnostic is the original 85Rb double-Lambda seeded-gain
    # model (regression-anchored); Biphoton is a newer, less-calibrated
    # spontaneous-FWM source estimate shared across cascade/diamond level
    # schemes rather than fit per atom/transition.
    caption = ("Seeded ⁸⁵Rb D1 mean-field Squeezing indicator and spontaneous "
               "biphoton estimates.")

    def param_schema(self):
        seeded = {"mode": MODE_SEEDED}
        biphoton = {"mode": MODE_BIPHOTON}
        cs_btw = {"mode": MODE_BIPHOTON, "topology": TOPOLOGY_CS_BTW}
        diamond = {"mode": MODE_BIPHOTON, "topology": TOPOLOGY_DIAMOND}
        return [
            ParamSpec("mode", "Mode", "Mode", MODE_SEEDED,
                      choices=(MODE_SEEDED, MODE_BIPHOTON),
                      choice_labels=MODE_LABELS,
                      control="segmented", applies_defaults=True,
                      help="Switches between seeded squeezing and a biphoton "
                           "source model, then loads that mode's defaults."),
            ParamSpec("topology", "Topology", "Model", TOPOLOGY_RB87_TELECOM,
                      choices=(TOPOLOGY_RB87_TELECOM, TOPOLOGY_CS_BTW, TOPOLOGY_DIAMOND),
                      choice_labels=TOPOLOGY_LABELS,
                      visible_if=biphoton, applies_defaults=True,
                      help="Atomic level scheme for the biphoton source."),
            ParamSpec("cs_channel", "Cs wavelength pair", "Model", CS_CHANNEL_917,
                      choices=(CS_CHANNEL_917, CS_CHANNEL_795), visible_if=cs_btw,
                      choice_labels=CS_CHANNEL_LABELS,
                      applies_defaults=True,
                      help="Selects the two wavelengths in the Cs cascade."),
            ParamSpec("biphoton_model", "Source model", "Model", BIPHOTON_PREDICTIVE,
                      choices=BIPHOTON_MODELS, control="segmented",
                      choice_labels=BIPHOTON_MODEL_LABELS,
                      visible_if=biphoton, advanced=True,
                      help="Reduced model calculates the Doppler waveform; its pair-"
                           "rate scale remains reference-anchored. Reference model "
                           "reproduces the stored literature calibration."),
            ParamSpec("opd", "One-photon detuning Δ", "Detunings", 0.9,
                      -3.0, 3.0, 0.1, "GHz",
                      visible_if=seeded,
                      help="Pump detuning from the F=2→F′=3 transition."),
            ParamSpec("tpd", "Two-photon detuning δ", "Detunings", -8.0,
                      -TPD_LIMIT_MHZ, TPD_LIMIT_MHZ, 1.0, "MHz", recompute=False,
                      visible_if=seeded,
                      help="Selects the operating point on the cached curve."),
            ParamSpec("temp_c", "Temperature", "Cell", 121.0,
                      60.0, 150.0, 1.0, "°C", visible_if=seeded),
            ParamSpec("cell_mm", "Cell length", "Cell", 12.5,
                      1.0, 100.0, 0.5, "mm", visible_if=seeded,
                      help="Vapor-cell length L. Enters the Maxwell-Bloch "
                           "propagation exp(M·L), so it recomputes the gain."),
            ParamSpec("transit_rate_khz", "Transit rate", "Cell", 100.0,
                      1.0, 1000.0, 1.0, "kHz", visible_if=seeded, advanced=True,
                      advanced_group="Model provenance",
                      help="Atom replacement rate; 100 kHz is an inherited estimate."),
            ParamSpec("pump_mw", "Pump power", "Beams", 600.0,
                      50.0, 1200.0, 10.0, "mW", visible_if=seeded),
            ParamSpec("probe_uw", "Seed power", "Beams", 8.0,
                      1.0, 200.0, 1.0, "µW", visible_if=seeded,
                      help="Seed power entering the cell."),
            ParamSpec("eom_residual_carrier_uw", "Residual EOM carrier power",
                      "Beams", 0.0, 0.0, 5000.0, 0.1, "µW",
                      visible_if=seeded, advanced=True, hidden=True,
                      advanced_group="EOM spectrum provenance",
                      help="Cell-plane residual carrier power. It is recorded and "
                           "claim-gated but not applied because its coupling, phase, "
                           "polarization, and noise transfer are uncalibrated."),
            ParamSpec("eom_other_sidebands_uw", "Other EOM sidebands power",
                      "Beams", 0.0, 0.0, 5000.0, 0.1, "µW",
                      visible_if=seeded, advanced=True, hidden=True,
                      advanced_group="EOM spectrum provenance",
                      help="Total cell-plane power outside the wanted sideband and "
                           "carrier. Recorded as unapplied provenance only."),
            ParamSpec("pump_waist_um", "Pump waist", "Beams",
                      W_PUMP * 1e6, 50.0, 2000.0, 10.0, "µm",
                      visible_if=seeded,
                      help="Pump 1/e² radius at the cell; sets pump intensity and "
                           "spatial overlap."),
            ParamSpec("probe_waist_um", "Seed waist", "Beams",
                      W_PROBE * 1e6, 50.0, 2000.0, 10.0, "µm",
                      visible_if=seeded,
                      help="Seed 1/e² radius at the cell; used in the seeded drive "
                           "and spatial overlap."),
            ParamSpec("seeded_angle_deg", "Pump–seed angle", "Beams",
                      SEEDED_PHASE_ANGLE_DEG, 0.0, 2.0, 0.05, "°",
                      visible_if=seeded,
                      help="Crossing angle used in the longitudinal phase mismatch."),
            ParamSpec("detection_eff_pct", "Detection efficiency η",
                      "Detection & scaling", SEEDED_DETECTION_EFFICIENCY_PCT,
                      0.0, 100.0, 0.1, "%", visible_if=seeded,
                      help="Total efficiency after the cell, including optical and "
                           "detector loss."),
            ParamSpec("loss_pct", "Loss after cell", "Detection & scaling",
                      SEEDED_POST_CELL_LOSS_PCT,
                      0.0, 50.0, 0.5, "%", visible_if=seeded, hidden=True),
            ParamSpec("qe_pct", "Detector quantum efficiency",
                      "Detection & scaling", QE_DETECTOR * 100.0,
                      50.0, 100.0, 0.01, "%", visible_if=seeded, advanced=True,
                      advanced_group="Detector", hidden=True,
                      help="Photodiode QE. With the post-cell loss knob it forms "
                           "eta = QE·(1−loss) in the gain-referred diagnostic. The "
                           "92% default is a historical model input, not a validation "
                           "of physical squeezing or a measured device calibration."),
            ParamSpec("line_strength", "Inherited residual factor",
                      "Detection & scaling", SEEDED_REFERENCE_RESIDUAL,
                      0.2, 5.0, 0.01, "×",
                      visible_if=seeded, advanced=True, hidden=True,
                      advanced_group="Seeded coupling factorization",
                      help="Backward-compatible dimensionless residual. The physical "
                           "macroscopic normalization — Rb85 D1 hyperfine Clebsch-Gordan "
                           "strengths × 1/[2(2I+1)] — is computed from first principles. "
                           "The trace-normalized rho_ss supplies the manifold population "
                           "exactly once. The inherited 0.74 residual has not been "
                           "refitted after the normalization and Maxwell-convention "
                           "corrections, so it does not anchor measured gain or squeezing. "
                           "It remains separate "
                           "because no measurements identify a unique split of 0.74 among "
                           "overlap, polarization, and Zeeman participation; the three "
                           "unit-default penalties below represent only additional loss "
                           "relative to that anchored setup."),
            ParamSpec("mode_overlap_penalty", "Additional mode-overlap penalty",
                      "Detection & scaling", 1.0,
                      0.0, 1.0, 0.01, "×", visible_if=seeded, advanced=True,
                      hidden=True,
                      advanced_group="Seeded coupling factorization",
                      help="One-sided effective-coupling penalty for unresolved "
                           "transverse pump-seed mode mismatch beyond what the 0.74 "
                           "reference anchor already absorbs. Unity adds no penalty. "
                           "This is separate from Ultra's normalized axial crossing-angle "
                           "profile, so do not use it to count that geometry twice."),
            ParamSpec("polarization_penalty", "Additional polarization penalty",
                      "Detection & scaling", 1.0, 0.0, 1.0, 0.01, "×",
                      visible_if=seeded, advanced=True, hidden=True,
                      advanced_group="Seeded coupling factorization",
                      help="One-sided effective-coupling penalty for extra leakage out "
                           "of the intended double-Lambda polarization channel relative "
                           "to the anchored reference. It is not a raw Stokes or "
                           "optical-power purity."),
            ParamSpec("zeeman_participation_penalty",
                      "Additional Zeeman-participation penalty",
                      "Detection & scaling", 1.0, 0.0, 1.0, 0.01, "×",
                      visible_if=seeded, advanced=True, hidden=True,
                      advanced_group="Seeded coupling factorization",
                      help="One-sided effective-coupling penalty for additional loss of "
                           "addressed m_F pathways relative to the anchored reference. "
                           "Unity keeps the lumped hyperfine model; this scalar does not "
                           "turn on the 24-level Zeeman diagnostic or a full Floquet solve."),
            ParamSpec("biphoton_temp_c", "Temperature", "Cell & beams", 90.0,
                      30.0, 160.0, 1.0, "°C", visible_if=biphoton),
            ParamSpec("biphoton_cell_mm", "Cell length", "Cell & beams", 12.5,
                      1.0, 100.0, 0.5, "mm", visible_if=biphoton,
                      advanced=True),
            ParamSpec("pump_biphoton_uw", "Pump power", "Fields", 10.0,
                      0.1, 200.0, 0.1, "µW", visible_if=biphoton),
            ParamSpec("coupling_mw", "Coupling strength", "Fields", 1.0,
                      0.01, 50.0, 0.01, "×", visible_if=biphoton,
                      help="Relative drive amplitude; this is a scale factor, not "
                           "an absolute power."),
            ParamSpec("pump_detuning_mhz", "Pump detuning", "Detunings", 0.0,
                      -2000.0, 2000.0, 10.0, "MHz", visible_if=biphoton),
            ParamSpec("two_photon_detuning_mhz", "Two-photon detuning", "Detunings", 0.0,
                      -2000.0, 2000.0, 10.0, "MHz", visible_if=biphoton,
                      help="Delta_p + Delta_c. Internally this sets the coupling "
                           "detuning relative to the pump detuning."),
            ParamSpec("coupling_detuning_mhz", "Coupling detuning", "Detunings", 0.0,
                      -2000.0, 2000.0, 10.0, "MHz", visible_if=biphoton,
                      hidden=True, recompute=False),
            ParamSpec("signal_angle_deg", "Signal collection angle", "Phase matching", 1.5,
                      0.0, 10.0, 0.1, "°", visible_if=biphoton),
            ParamSpec("idler_angle_offset_deg", "Idler angle offset", "Phase matching",
                      0.0, -5.0, 5.0, 0.05, "°", visible_if=biphoton,
                      help="Offset from the transverse phase-matched idler angle "
                           "derived from the selected topology and signal angle."),
            ParamSpec("idler_angle_deg", "Idler angle", "Phase matching",
                      transverse_matched_angle_deg(1529.37, 780.24, 1.5),
                      0.0, 10.0, 0.1, "°", visible_if=biphoton,
                      hidden=True, recompute=False),
            ParamSpec("signal_side", "Signal side", "Phase matching", SIDE_PLUS,
                      choices=SIDE_CHOICES, visible_if=biphoton,
                      hidden=True,
                      help="Transverse collection side used in vector phase matching."),
            ParamSpec("idler_side", "Idler side", "Phase matching", SIDE_PLUS,
                      choices=SIDE_CHOICES, visible_if=biphoton,
                      hidden=True,
                      help="Opposite side flips the idler transverse wavevector."),
            ParamSpec("diamond_pump_nm", "Diamond pump wavelength", "Fields", 780.0,
                      300.0, 2000.0, 1.0, "nm", visible_if=diamond,
                      advanced=True),
            ParamSpec("diamond_coupling_nm", "Diamond coupling wavelength", "Fields", 776.0,
                      300.0, 2000.0, 1.0, "nm", visible_if=diamond,
                      advanced=True),
            ParamSpec("diamond_signal_nm", "Diamond signal wavelength", "Fields", 795.0,
                      300.0, 2000.0, 1.0, "nm", visible_if=diamond,
                      advanced=True),
            ParamSpec("diamond_idler_nm", "Diamond idler wavelength", "Fields", 761.702,
                      300.0, 2500.0, 0.001, "nm", visible_if=diamond,
                      advanced=True),
            ParamSpec("signal_eff_pct", "Signal efficiency", "Detection & scaling", 10.0,
                      0.1, 95.0, 0.1, "%", visible_if=biphoton,
                      advanced=True, recompute=False),
            ParamSpec("idler_eff_pct", "Idler efficiency", "Detection & scaling", 10.0,
                      0.1, 95.0, 0.1, "%", visible_if=biphoton,
                      advanced=True, recompute=False),
            ParamSpec("dark_signal_cps", "Signal background", "Detection & scaling", 2000.0,
                      0.0, 100000.0, 100.0, "cps", visible_if=biphoton,
                      advanced=True, recompute=False),
            ParamSpec("dark_idler_cps", "Idler background", "Detection & scaling", 2000.0,
                      0.0, 100000.0, 100.0, "cps", visible_if=biphoton,
                      advanced=True, recompute=False),
            ParamSpec("coincidence_window_ns", "Coincidence window", "Detection & scaling", 1.0,
                      0.01, 100.0, 0.01, "ns", visible_if=biphoton,
                      recompute=False),
            ParamSpec("timing_jitter_ns", "Timing jitter FWHM", "Detection & scaling", 0.55,
                      0.0, 5.0, 0.01, "ns", visible_if=biphoton,
                      advanced=True, recompute=False,
                      help="Net Gaussian signal-idler timing-difference response "
                           "FWHM. It broadens the detected correlation only; the "
                           "intrinsic source FWHM is reported separately."),
            ParamSpec("filter_bandwidth_mhz", "Filter bandwidth", "Detection & scaling", 300.0,
                      1.0, 5000.0, 1.0, "MHz", visible_if=biphoton,
                      recompute=False),
            ParamSpec("tau_max_ns", "Temporal window", "Numerics", 12.0,
                      1.0, 100.0, 1.0, "ns", visible_if=biphoton, advanced=True),
            ParamSpec("biphoton_velocity_step", "Velocity step", "Numerics", 2.0,
                      0.5, 20.0, 0.5, "m/s", visible_if=biphoton, advanced=True,
                      help="Maxwell velocity-grid step. The calibrated source model "
                           "uses it directly; the predictive model treats it as an "
                           "upper bound and auto-refines finer until the biphoton "
                           "width converges (a coarse step aliases the velocity-"
                           "class coherent sum)."),
            ParamSpec("resolution", "Solver detail", "Numerics", FIDELITY_FAST,
                      choices=tuple(FWM_FIDELITY.keys()), advanced=True,
                      choice_labels=FIDELITY_LABELS, visible_if=seeded,
                      help="Controls the scan grid and propagation refinements."),
            ParamSpec("floquet_order", "Floquet order", "Numerics",
                      SEEDED_FLOQUET_ORDER, choices=(3, 4, 5), advanced=True,
                      visible_if=seeded, hidden=True,
                      help="Retained harmonics -N_F,...,+N_F. Every reported scan "
                           "is compared with N_F-1 over all complex response and "
                           "transfer coefficients; a failed gate is returned as "
                           "UNCONVERGED."),
            ParamSpec("phase_detail", "Phase-matching view", "Phase matching", "Balanced",
                      choices=("Balanced", "Fine"), visible_if=biphoton,
                      choice_labels={"Balanced": "1D", "Fine": "1D + 2D"},
                      control="segmented", advanced=True, hidden=True,
                      recompute=False,
                      help="Adds the 2D signal–idler acceptance map without rerunning "
                           "the source model."),
        ]

    def presets(self):
        return []

    def recommended_defaults(self, params):
        return {
            MODE_SEEDED: self._squeezing_defaults(),
            MODE_BIPHOTON: self._biphoton_defaults(params),
        }

    def _squeezing_defaults(self):
        return dict(mode=MODE_SEEDED, opd=0.9, tpd=-8.0, temp_c=121.0,
                    cell_mm=12.5, pump_mw=600.0, probe_uw=8.0,
                    detection_eff_pct=SEEDED_DETECTION_EFFICIENCY_PCT,
                    loss_pct=SEEDED_POST_CELL_LOSS_PCT,
                    transit_rate_khz=100.0,
                    eom_residual_carrier_uw=0.0,
                    eom_other_sidebands_uw=0.0,
                    pump_waist_um=W_PUMP * 1e6, probe_waist_um=W_PROBE * 1e6,
                    qe_pct=QE_DETECTOR * 100.0,
                    line_strength=SEEDED_REFERENCE_RESIDUAL,
                    mode_overlap_penalty=1.0,
                    polarization_penalty=1.0,
                     zeeman_participation_penalty=1.0,
                     resolution=FIDELITY_FAST,
                     floquet_order=SEEDED_FLOQUET_ORDER,
                     seeded_angle_deg=SEEDED_PHASE_ANGLE_DEG)

    def _biphoton_defaults(self, params):
        topology = params.get("topology", TOPOLOGY_RB87_TELECOM)
        base = dict(
            mode=MODE_BIPHOTON,
            topology=topology,
            biphoton_model=params.get("biphoton_model", BIPHOTON_PREDICTIVE),
            biphoton_cell_mm=12.5,
            pump_detuning_mhz=0.0,
            two_photon_detuning_mhz=0.0,
            coupling_detuning_mhz=0.0,
            signal_angle_deg=1.5,
            idler_angle_offset_deg=0.0,
            idler_angle_deg=transverse_matched_angle_deg(1529.37, 780.24, 1.5),
            signal_side=SIDE_PLUS,
            idler_side=SIDE_PLUS,
            signal_eff_pct=10.0,
            idler_eff_pct=10.0,
            dark_signal_cps=2000.0,
            dark_idler_cps=2000.0,
            coincidence_window_ns=1.0,
            filter_bandwidth_mhz=300.0,
            timing_jitter_ns=0.55,
            tau_max_ns=12.0,
            biphoton_velocity_step=2.0,
            phase_detail="Balanced",
        )
        if topology == TOPOLOGY_CS_BTW:
            base.update(cs_channel=params.get("cs_channel", CS_CHANNEL_917),
                        biphoton_temp_c=75.0, pump_biphoton_uw=20.0,
                        coupling_mw=1.0)
        elif topology == TOPOLOGY_DIAMOND:
            base.update(biphoton_temp_c=60.0, pump_biphoton_uw=20.0,
                        coupling_mw=1.0, diamond_pump_nm=780.0,
                        diamond_coupling_nm=776.0, diamond_signal_nm=795.0,
                        diamond_idler_nm=761.702)
        else:
            base.update(topology=TOPOLOGY_RB87_TELECOM, biphoton_temp_c=90.0,
                        pump_biphoton_uw=10.0, coupling_mw=1.0)
        base.update(_default_biphoton_geometry(base))
        return base

    def _biphoton_runtime_params(self, params):
        """Map the compact lab-facing controls onto the backend parameters."""
        out = dict(params)
        if out.get("mode", MODE_SEEDED) != MODE_BIPHOTON:
            return out

        topology = topology_from_params(out)
        signal = topology.field_map["signal"]
        idler = topology.field_map["idler"]
        signal_angle = float(out.get("signal_angle_deg", signal.angle_deg))
        out["signal_angle_deg"] = signal_angle

        if "two_photon_detuning_mhz" in out:
            two_det = float(out.get("two_photon_detuning_mhz", 0.0))
            pump_det = float(out.get("pump_detuning_mhz", 0.0))
            out["coupling_detuning_mhz"] = two_det - pump_det
        else:
            out["two_photon_detuning_mhz"] = (
                float(out.get("pump_detuning_mhz", 0.0))
                + float(out.get("coupling_detuning_mhz", 0.0))
            )

        if "idler_angle_offset_deg" in out:
            matched_idler = transverse_matched_angle_deg(
                signal.wavelength_nm, idler.wavelength_nm, signal_angle)
            out["idler_angle_deg"] = float(np.clip(
                matched_idler + float(out.get("idler_angle_offset_deg", 0.0)),
                0.0, 10.0))
        return out

    def info(self):
        return (
            "**Seeded FWM.** Computes mean-field seed and conjugate gain for the "
            "⁸⁵Rb D1 double-Λ system. The Squeezing indicator is derived from those "
            "gains; without atomic Langevin covariance it shows trends, not a "
            "physical squeezing spectrum. A pump-energy cap prevents unbounded "
            "small-signal gain but is not a depleted three-field solve. Numerical "
            "checks and solver provenance are listed under Model diagnostics.\n\n"
            "**Biphoton.** The Reduced model calculates a Doppler-averaged waveform "
            "and vector phase matching. Absolute widths remain approximate and the "
            "pair-rate scale is anchored to a literature reference. The Reference "
            "model reproduces the stored calibration; agreement by construction is "
            "not an independent validation.\n\n"
            "**References**\n"
            "- G. Sim, H. Kim, H. S. Moon, *Sci. Rep.* **15**, 7727 (2025) "
            "(seeded 85Rb D1 double-Lambda gain & squeezing, regression anchor).\n"
            "- H. Kim, H. Jeong, H. S. Moon, [*Quantum Sci. Technol.* **9**, "
            "045006 (2024)](https://arxiv.org/abs/2402.06872) (Cs cascade BTW, Eq. 2).\n"
            "- H. Jeong, H. Kim, H. S. Moon, *Adv. Quantum Technol.* **7**, "
            "2300108 (2024) (87Rb telecom).\n"
            "- S. Du, J. Wen, M. H. Rubin, [*J. Opt. Soc. Am. B* **25**, C98 "
            "(2008)](https://arxiv.org/abs/0804.3981) (biphoton = nonlinear ⊛ "
            "linear response, Eq. 15).\n"
            "- Chen *et al.*, [*Phys. Rev. Research* **4**, 023132 (2022)]"
            "(https://arxiv.org/abs/2109.09062) (Doppler-averaged hot-vapor SFWM, Eq. 3-5).\n"
            "- J. Park, T. Jeong, H. S. Moon, [*Sci. Rep.* **10**, 16413 "
            "(2020)](https://www.nature.com/articles/s41598-020-73610-2) "
            "(cascade-type warm-atom biphoton waveform)."
        )

    def compute(self, params):
        if params.get("mode", MODE_SEEDED) == MODE_BIPHOTON:
            params = self._biphoton_runtime_params(params)
            return _solver_from_topology_key(
                _topology_cache_key(params)).compute_biphoton(params)
        center = branch_center_GHz(params["opd"], -1)
        fidelity = normalize_fidelity(params["resolution"])
        res = FWM_FIDELITY[fidelity]
        wanted_seed_uw = params["probe_uw"]
        alias_seed_uw = params.get("seed_wanted_sideband_uw", wanted_seed_uw)
        if not np.isclose(alias_seed_uw, wanted_seed_uw, rtol=0.0, atol=1e-12):
            raise ValueError(
                "probe_uw and seed_wanted_sideband_uw disagree; use one "
                "wanted-sideband power at the bridge boundary")
        return compute_spectrum(
            params["opd"],
            T=params["temp_c"] + 273.15,
            P_pump=params["pump_mw"] * 1e-3,
            P_probe=wanted_seed_uw * 1e-6,
            line_strength=params["line_strength"],
            mode_overlap_penalty=params.get("mode_overlap_penalty", 1.0),
            polarization_penalty=params.get("polarization_penalty", 1.0),
            zeeman_participation_penalty=params.get(
                "zeeman_participation_penalty", 1.0),
            L=params["cell_mm"] * 1e-3,
            w_pump=params.get("pump_waist_um", W_PUMP * 1e6) * 1e-6,
            w_probe=params.get("probe_waist_um", W_PROBE * 1e6) * 1e-6,
            qe=params.get("qe_pct", QE_DETECTOR * 100.0) / 100.0,
            loss_frac=params["loss_pct"] / 100.0,
            detection_efficiency=_detection_efficiency_from_params(params),
            coarse_points=res["coarse_points"], fine_points=0,
            scan_min=center - WINDOW_GHZ, scan_max=center + WINDOW_GHZ,
            velocity_step=res["velocity_step"],
            velocity_cutoff=res.get("velocity_cutoff", 3.0),
            phase_detail=res["phase_detail"],
            pump_probe_angle_deg=params.get("seeded_angle_deg", SEEDED_PHASE_ANGLE_DEG),
            model_fidelity=fidelity,
            floquet_order=params.get("floquet_order", SEEDED_FLOQUET_ORDER),
            enforce_floquet_convergence=True,
            transit_rate=(2.0 * np.pi
                          * (params.get("transit_rate_khz", 100.0) * 1e3)),
            eom_residual_carrier_power=(
                params.get("eom_residual_carrier_uw", 0.0) * 1e-6),
            eom_other_sidebands_power=(
                params.get("eom_other_sidebands_uw", 0.0) * 1e-6),
            eom_seed_spectrum_provenance=params.get(
                "eom_seed_spectrum_provenance", "direct FWM input"),
            eom_seed_spectrum_status=params.get(
                "eom_seed_spectrum_status", "not supplied"),
            eom_seed_spectrum_application=params.get(
                "eom_seed_spectrum_application", "unapplied"),
            branch=DEFAULT_BRANCH,
        )

    def observables(self, raw, params, include_figures=True):
        if raw.get("kind") == "biphoton":
            params = self._biphoton_runtime_params(params)
            return self._biphoton_observables(
                raw, params, include_figures=include_figures)
        return self._seeded_observables(
            raw, params, include_figures=include_figures)

    def _seeded_observables(self, raw, params, include_figures=True):
        tpd = params["tpd"]
        op = operating_point(raw, tpd, branch=-1)
        d_axis = (raw["probe_axis_GHz"] - raw["raman_center_minus_GHz"]) * 1e3
        claim_gate = raw.get("claim_gate", {})
        floquet_convergence = raw.get("floquet_convergence", {})

        fig = None
        if include_figures:
            import matplotlib.pyplot as plt

            fig, (axG, axS) = plt.subplots(2, 1, figsize=(8.5, 6.4), sharex=True)
            for ax in (axG, axS):
                ax.grid(alpha=0.3)
            axG.plot(d_axis, raw["G_s"], color="#1f77b4", lw=1.8)
            axG.axvline(tpd, color="crimson", ls="--", lw=1.2)
            axG.axhline(1.0, color="black", lw=0.6)
            axG.scatter([tpd], [op["G_s"]], color="crimson", zorder=5)
            axG.set_ylabel("Seed gain G_s")
            axG.set_title(f"Delta = {params['opd']:.1f} GHz,  "
                          f"T = {params['temp_c']:.0f} C,  eta = {raw['eta']:.3f}")
            if np.nanmax(raw["G_s"]) > 50:
                axG.set_yscale("log")
            axS.plot(d_axis, raw["gain_referred_noise_dB"], color="#2ca02c", lw=1.8)
            axS.axvline(tpd, color="crimson", ls="--", lw=1.2)
            axS.axhline(0.0, color="black", lw=0.6)
            axS.scatter([tpd], [op["gain_referred_noise_dB"]],
                        color="crimson", zorder=5)
            axS.set_ylabel("Squeezing indicator [dB]")
            axS.set_xlabel(
                "Two-photon detuning delta [MHz]   (probe on the - Raman branch)")
            axS.set_xlim(-TPD_LIMIT_MHZ, TPD_LIMIT_MHZ)
            fig.tight_layout()

        metrics = [
            dict(label="Squeezing indicator",
                 value=f"{op['gain_referred_noise_dB']:.2f} dB",
                 delta="physical squeezing unavailable",
                 help="Mean-field gain estimate; atomic Langevin covariance is not "
                      "included.", tier="hero"),
            dict(label="Seed gain G_s", value=f"{op['G_s']:.2f}",
                 help="Mean-field seed power gain at the selected detuning.",
                 tier="hero"),
            dict(label="Conjugate gain G_c", value=f"{op['G_c']:.2f}",
                 help="Generated conjugate power gain at the selected detuning."),
        ]
        cap = raw.get("pump_depletion_cap", float("inf"))
        small_signal = raw.get("G_s_smallsignal_peak", op["G_s"])
        depletion_limited = small_signal > 1.1 * cap
        solver = raw.get("atomic_solver_provenance", {})
        pump_reference = raw.get(
            "pump_weak_response_reference_provenance", {})
        noncollinear_reference = raw.get(
            "noncollinear_doppler_reference_provenance", {})
        solver_detail = FIDELITY_LABELS.get(
            raw.get("model_fidelity"), raw.get("model_fidelity", "unavailable"))
        phase_rows = ""
        if raw.get("phase_detail", PHASE_LEGACY) != PHASE_LEGACY:
            phase_rows = (
                f"| Solver detail | {solver_detail} |\n"
                f"| Phase model | {raw['phase_detail']} |\n"
                f"| Pump–seed angle | {raw['pump_probe_angle_deg']:.3f}° |\n"
                f"| Operating Δk_z | {op.get('delta_k_z', np.nan):.3e} 1/m |\n"
                f"| Vacuum Δk_z | {op.get('delta_k_z_vacuum', np.nan):.3e} 1/m |\n"
                f"| Propagation segments | {raw.get('phase_segments', 1)} |\n"
                f"| Segment absorption OD estimate | {raw.get('segment_absorption_od', 0.0):.3f} |\n"
            )
            if raw.get("phase_detail") == PHASE_ULTRA:
                phase_rows += (
                    f"| Dispersion convention | Option A (diagonal chi only) |\n"
                    f"| Dynamic depletion | {raw.get('ultra_dynamic_depletion', False)} |\n"
                    f"| Min pump remaining | {raw.get('ultra_pump_remaining_min', np.nan):.3e} W |\n"
                    f"| Distributed atomic Langevin covariance | unavailable |\n"
                    f"| Min Gaussian overlap | {raw.get('ultra_spatial_overlap_min', 1.0):.4f} |\n"
                    f"| Zeeman status | {raw.get('zeeman_status', 'inactive')} |\n"
                    f"| Zeeman CG-sum diagnostic | {raw.get('zeeman_correction', 1.0):.4f} |\n"
                )
        depletion_warning = (
            "\n⚠️ **Pump-depletion limited:** the small-signal gain exceeds the "
            "pump-energy budget, so the displayed gain is capped. Lower the "
            "temperature or raise the seed power to return to the linear regime.\n"
            if depletion_limited else ""
        )
        operating_table = (
            f"| Quantity | Value |\n|---|---|\n"
            f"| ⁸⁵Rb density | {raw['N_atoms']:.3e} /m³ |\n"
            f"| Thermal velocity σ_v | {raw['sigma_v']:.1f} m/s |\n"
            f"| Velocity classes | {raw['n_velocity']} |\n"
            f"| Ω_pump / 2π | {raw['Op_A_2pi_GHz']:.3f} GHz |\n"
            f"| Ω_seed / 2π | {raw['Os_2pi_MHz']:.3f} MHz |\n"
            f"| (−) Raman line (probe axis) | {raw['raman_center_minus_GHz']:.3f} GHz |\n"
            f"| Detection efficiency η | {raw['eta']:.4f} |\n"
            f"| Operating probe detuning | {op['probe_GHz']:.4f} GHz |\n"
            f"| Cell length | {raw.get('cell_length_m', L_CELL)*1e3:.1f} mm |\n"
            f"| Pump waist | {raw.get('w_pump_m', W_PUMP)*1e6:.0f} µm |\n"
            f"| Seed waist | {raw.get('w_probe_m', W_PROBE)*1e6:.0f} µm |\n"
            f"| Pump–seed angle | {raw.get('pump_probe_angle_deg', SEEDED_PHASE_ANGLE_DEG):.3f}° |\n"
            + depletion_warning
        )
        claim_reasons = "; ".join(claim_gate.get("reasons", ())) or "none recorded"
        diagnostics_table = (
            f"| Check | Value |\n|---|---|\n"
            f"| Model scope | Gain trends only |\n"
            f"| Validation level | {claim_gate.get('level', 'unavailable')} |\n"
            f"| Claim-gate reasons | {claim_reasons} |\n"
            f"| Quantitative gain claim | "
            f"{'supported' if claim_gate.get('quantitative_gain_supported') else 'UNSUPPORTED'} |\n"
            f"| Physical squeezing claim | "
            f"{'supported' if claim_gate.get('physical_squeezing_prediction') else 'UNAVAILABLE'} |\n"
            f"| Coupling norm 1/[2(2I+1)] (rho_ss supplies population) | "
            f"{raw.get('coupling_norm', float('nan')):.4f} |\n"
            f"| Canonical-mode status | {raw.get('canonical_mode_status', 'unavailable')} |\n"
            f"| Probe / conjugate effective area | "
            f"{raw.get('canonical_mode_area_probe_m2', np.nan):.3e} / "
            f"{raw.get('canonical_mode_area_conjugate_m2', np.nan):.3e} m² |\n"
            f"| Small-signal photon-flux gap | "
            f"{op.get('photon_flux_gap_smallsignal', np.nan):.6f} |\n"
            f"| Bare-map commutator defect max | "
            f"{op.get('commutator_defect_max_smallsignal', np.nan):.3e} |\n"
            f"| Atomic solver | {solver.get('solver_id', 'unavailable')} |\n"
            f"| Pump-only weak-response reference | "
            f"{pump_reference.get('solver_id', 'unavailable')} (separate/slow) |\n"
            f"| Pump-reference supported branch | "
            f"{pump_reference.get('supported_branches', ('unavailable',))} |\n"
            f"| 2-D Raman-Doppler reference | "
            f"{noncollinear_reference.get('solver_id', 'unavailable')} "
            f"(separate/slow; production 1-D) |\n"
            f"| 2-D reference velocity geometry | "
            f"{noncollinear_reference.get('velocity_geometry', 'unavailable')} |\n"
            f"| Floquet order | N_F={solver.get('floquet_order', 'unavailable')} |\n"
            f"| Floquet adjacent-order gate | "
            f"{floquet_convergence.get('status', 'NOT_EVALUATED')} vs "
            f"N_F={floquet_convergence.get('comparison_order', 'unavailable')} |\n"
            f"| Floquet full-scan points | "
            f"{floquet_convergence.get('full_scan_points', 0)} |\n"
            f"| Noise-trace status | {raw.get('squeezing_status', 'unavailable')} |\n"
            f"| Transit reset γ_t / 2π | "
            f"{raw.get('transit_reset_rate_rad_s', np.nan)/(2*np.pi*1e3):.3f} kHz |\n"
            f"| Combined residual coupling | "
            f"{raw.get('combined_line_strength_residual', raw.get('line_strength_residual', float('nan'))):.4f}× |\n"
            f"| Effective coupling scale | {raw.get('effective_line_strength', float('nan')):.4f} |\n"
            f"| Pump-depletion cap on G_s (Manley-Rowe) | {cap:.3e} |\n"
            f"| Small-signal peak G_s (pre-saturation) | {small_signal:.3e} |\n"
            + phase_rows
        )
        return {
            "metrics": metrics,
            "figure": fig,
            "tables": [
                {"title": "Operating point", "markdown": operating_table},
                {"title": "Model diagnostics", "markdown": diagnostics_table},
            ],
        }

    def _biphoton_observables(self, raw, params, include_figures=True):
        topo = raw["topology"]
        predictive = raw.get("predictive", False)
        velocity_converged = raw.get("velocity_converged", True)
        requested_detail = params.get(
            "phase_detail", raw.get("phase_detail", "Balanced"))
        phase_detail = (
            "Fine" if str(requested_detail).lower() == "fine" else "Balanced")
        # Predictive: g²_SI(τ) comes from the computed waveform |ψ|² and the
        # physical accidentals (no target-g² forcing). Calibrated: legacy
        # added-accidental anchoring to the reference g² peak.
        stats = observables.biphoton_stats(
            raw["tau_axis_ns"], raw["psi_tau"], raw["pair_rate_cps"],
            signal_eff=params["signal_eff_pct"] / 100.0,
            idler_eff=params["idler_eff_pct"] / 100.0,
            dark_signal_cps=params["dark_signal_cps"],
            dark_idler_cps=params["dark_idler_cps"],
            coincidence_window_ns=params["coincidence_window_ns"],
            timing_jitter_ns=params["timing_jitter_ns"],
            filter_bandwidth_mhz=params["filter_bandwidth_mhz"],
            source_bandwidth_mhz=raw["source_bandwidth_mhz"],
            target_g2_peak=None if predictive else topo.target_g2_peak,
        )
        width_is_reportable = velocity_converged or not predictive
        source_width_value = (
            f"{stats['source_fwhm_ns']:.2f} ns"
            if width_is_reportable else "unconverged")
        detected_width_value = (
            f"{stats['detected_fwhm_ns']:.2f} ns"
            if width_is_reportable else "unconverged")
        source_width_table = (
            f"{stats['source_fwhm_ns']:.3f} ns"
            if width_is_reportable else "unconverged")
        detected_width_table = (
            f"{stats['detected_fwhm_ns']:.3f} ns"
            if width_is_reportable else "unconverged")

        fig = None
        extra_figures = []
        if include_figures:
            import matplotlib.pyplot as plt

            fig, (axG2, axPM) = plt.subplots(2, 1, figsize=(8.5, 6.4))
            axG2.plot(stats["tau_axis_ns"], stats["g2_SI_tau"],
                      color="#1f77b4", lw=1.8)
            axG2.axhline(2.0, color="black", lw=0.7, ls=":")
            axG2.set_ylabel(r"$g^{(2)}_{SI}(\tau)$")
            title_tag = ("Reduced model · reference-anchored rate" if predictive
                         else "Reference model")
            axG2.set_title(f"{topo.label}: {title_tag}")
            axG2.grid(alpha=0.3)

            angle_axis = raw.get("angle_axis_deg")
            phase_matching = raw.get("phase_matching")
            if angle_axis is None or phase_matching is None:
                angle_axis = np.linspace(max(params["idler_angle_deg"] - 4.0, 0.0),
                                         params["idler_angle_deg"] + 4.0, 181)
                angle_dk = np.array([
                    phase_mismatch_vector(
                        raw["fields"], idler_angle_deg=a,
                        reference_delta_k=topo.reference_delta_k
                    )["delta_k_vector"]
                    for a in angle_axis
                ])
                phase_matching = phase_matching_weight(angle_dk, raw["cell_length_m"])

            axPM.plot(angle_axis, phase_matching,
                      color="#2ca02c", lw=1.8)
            axPM.axvline(params["idler_angle_deg"], color="crimson", lw=1.1, ls="--")
            axPM.set_xlabel("Idler collection angle [deg]")
            axPM.set_ylabel(r"$\mathrm{sinc}^2(|\Delta\mathbf{k}| L / 2)$")
            axPM.grid(alpha=0.3)
            fig.tight_layout()

            amp = np.abs(raw["source_v"])
            amp = amp / np.nanmax(amp) if np.nanmax(amp) > 0 else amp
            phase = np.unwrap(np.angle(raw["source_v"]))
            figV, (axA, axP) = plt.subplots(2, 1, figsize=(8.5, 6.0), sharex=True)
            axA.plot(raw["v_grid"], amp, color="#ff7f0e", lw=1.5)
            axA.set_ylabel("Velocity source amplitude")
            axA.grid(alpha=0.3)
            axP.plot(raw["v_grid"], phase, color="#9467bd", lw=1.3)
            axP.set_xlabel("Atomic velocity [m/s]")
            axP.set_ylabel("Velocity source phase [rad]")
            axP.grid(alpha=0.3)
            figV.tight_layout()

            extra_figures = [("Velocity-class coherent source", figV)]
            if phase_detail == "Fine":
                signal_axis, idler_axis, pm2d = biphoton_phase_matching_map(
                    raw["fields"], raw["cell_length_m"],
                    signal_angle_deg=params["signal_angle_deg"],
                    idler_angle_deg=params["idler_angle_deg"],
                    reference_delta_k=topo.reference_delta_k)
                figPM2, ax2 = plt.subplots(1, 1, figsize=(7.2, 5.2))
                im = ax2.imshow(
                    pm2d.T, origin="lower", aspect="auto",
                    extent=[
                        signal_axis[0], signal_axis[-1],
                        idler_axis[0], idler_axis[-1],
                    ],
                    cmap="viridis", vmin=0.0, vmax=1.0)
                ax2.scatter([params["signal_angle_deg"]], [params["idler_angle_deg"]],
                            color="crimson", s=28, zorder=5)
                ax2.set_xlabel("Signal angle [deg]")
                ax2.set_ylabel("Idler angle [deg]")
                ax2.set_title("2D phase-matching acceptance")
                figPM2.colorbar(
                    im, ax=ax2,
                    label=r"$\mathrm{sinc}^2(|\Delta\mathbf{k}| L / 2)$")
                figPM2.tight_layout()
                extra_figures.append(("2D phase matching", figPM2))

        if predictive:
            g2_help = ("Peak signal-idler cross-correlation from the computed "
                       "waveform |ψ|² and physical accidentals (not forced).")
            rate_help = ("Reference-anchored pair rate after pump, coupling, density, "
                         "and phase-matching scaling.")
            source_fwhm_help = (
                "Intrinsic |ψ(τ)|² FWHM on an auto-refined velocity grid. The "
                "absolute width remains approximate."
                if velocity_converged else
                "Width withheld because velocity-grid refinement hit its point cap.")
            detected_fwhm_help = (
                "Intrinsic waveform convolved with the signal–idler timing response."
                if velocity_converged else
                "Width withheld because the intrinsic waveform did not converge.")
            status_value = (f"Reduced model · {raw.get('regime', '—')}"
                            if velocity_converged
                            else f"Reduced model · {raw.get('regime', '—')} · "
                                 "width unconverged")
            status = dict(
                label="Model status", value=status_value,
                help="Calculates the Doppler waveform and vector phase matching; "
                     "the pair-rate scale remains reference-anchored.")
        else:
            g2_help = ("Peak normalized signal–idler correlation from the stored "
                       "reference calibration.")
            rate_help = "Pair rate from the stored reference calibration."
            source_fwhm_help = (
                "Intrinsic FWHM of the modeled |ψ(τ)|² source waveform before "
                "detector timing response.")
            detected_fwhm_help = (
                "FWHM after convolution with the net Gaussian signal-idler timing-"
                "difference response.")
            status = dict(
                label="Model status", value="Reference model",
                help="Reproduces stored calibration values; it is not an independent "
                     "prediction.")
        metrics = [
            status,
            dict(label="Peak signal–idler correlation",
                 value=f"{stats['g2_peak']:.2f}", help=g2_help,
                 tier="hero"),
            dict(label="CAR", value=f"{stats['CAR']:.1f}",
                 help="True coincidence divided by accidental coincidence."),
            dict(label="Pair rate", value=f"{stats['pair_rate_cps']:.1f} cps",
                 help=rate_help, tier="hero"),
            dict(label="Intrinsic biphoton width",
                 value=source_width_value,
                 help=source_fwhm_help),
            dict(label="Detected biphoton width",
                 value=detected_width_value,
                 help=detected_fwhm_help),
            dict(label="Source bandwidth", value=f"{raw['source_bandwidth_mhz']:.0f} MHz",
                 help=("Spectral FWHM from the waveform (predictive)." if predictive
                       else "Reference source bandwidth.")),
            dict(label="Phase-matching efficiency",
                 value=f"{raw['phase_match_weight']:.3f}",
                 help="Vector sinc^2 phase-matching collection weight."),
        ]

        field_rows = "".join(
            f"| {f.role} | {raw['topology'].levels[f.lower].name} -> "
            f"{raw['topology'].levels[f.upper].name} | {f.wavelength_nm:.2f} nm | "
            f"{f.angle_deg:.2f}° | {_side_label(f.side_sign) if f.side_sign else '0'} |\n"
            for f in raw["fields"]
        )
        pm_warning = (
            "| Warning | Vector phase match < 1e-3; collection geometry is not "
            "physically phase matched. |\n"
            if raw["phase_match_weight"] < 1e-3 else ""
        )
        topology_table = (
            f"| Quantity | Value |\n|---|---|\n"
            f"| Topology | {topo.label} |\n"
            f"| Family | {topo.family} |\n"
            f"| Density | {raw['density']:.3e} /m^3 |\n"
            f"| Cell length | {raw['cell_length_m']*1e3:.2f} mm |\n"
            f"| Delta k z relative | {raw['delta_k_z_relative']:.3e} 1/m |\n"
            f"| Delta k z absolute | {raw['delta_k_z_absolute']:.3e} 1/m |\n"
            f"| Delta k x transverse | {raw['delta_k_x']:.3e} 1/m |\n"
            f"| Delta k vector | {raw['delta_k_vector']:.3e} 1/m |\n"
            f"| Vector phase match | {raw['phase_match_weight']:.3f} |\n"
            f"| Longitudinal phase match | {raw['phase_match_weight_longitudinal']:.3f} |\n"
            f"| Vacuum phase match | {raw['phase_match_weight_absolute']:.3f} |\n"
            + pm_warning +
            f"| Phase-matching view | "
            f"{'1D + 2D' if phase_detail == 'Fine' else '1D'} |\n"
            f"| Energy mismatch | {raw['energy_mismatch_hz']/1e6:.3f} MHz |\n"
            f"| Residual two-photon Doppler k | {raw['residual_two_photon_k']:.3e} 1/m |\n"
            f"| Velocity step{' (auto-refined)' if predictive else ''} | "
            f"{raw.get('velocity_step_used', float('nan')):.3f} m/s"
            f"{'' if velocity_converged else ' · UNCONVERGED'} |\n\n"
            f"| Field | Transition | Wavelength | Angle | Side |\n|---|---|---:|---:|---:|\n"
            + field_rows
        )
        detection_table = (
            "Source estimate with detector/background model:\n\n"
            f"| Quantity | Value |\n|---|---|\n"
            f"| Signal singles | {stats['singles_signal_cps']:.2f} cps |\n"
            f"| Idler singles | {stats['singles_idler_cps']:.2f} cps |\n"
            f"| True coincidence | {stats['coincidence_cps']:.3f} cps |\n"
            f"| Accidental coincidence | {stats['accidental_cps']:.3e} cps |\n"
            f"| Raw accidental before reference calibration | {stats['raw_accidental_cps']:.3e} cps |\n"
            f"| Added unmodelled accidental/background | {stats['added_accidental_cps']:.3e} cps |\n"
            f"| Raw g2 peak before reference calibration | {stats['raw_g2_peak']:.2f} |\n"
            f"| Heralding signal | {stats['heralding_signal']:.3e} |\n"
            f"| Heralding idler | {stats['heralding_idler']:.3e} |\n"
            f"| Cauchy-Schwarz R | {stats['cauchy_schwarz_R']:.2f} |\n"
            f"| Intrinsic biphoton width | {source_width_table} |\n"
            f"| Detected biphoton width | {detected_width_table} |\n"
            f"| Net timing-difference response FWHM | {params['timing_jitter_ns']:.3f} ns |\n"
            f"| Filter transmission estimate | {stats['filter_transmission']:.3f} |\n\n"
            f"{raw['notes']}"
        )
        validation_table = self._reference_validation_table(raw, params, stats)
        return {
            "metrics": metrics,
            "figure": fig,
            "figures": extra_figures,
            "figure_controls": ["phase_detail"],
            "tables": [
                {"title": "Source geometry", "markdown": topology_table},
                {"title": "Detection estimates", "markdown": detection_table},
                {"title": "Literature comparison", "markdown": validation_table},
            ],
        }

    def _reference_validation_table(self, raw, params, stats):
        topo = raw["topology"]
        predictive = raw.get("predictive", False)
        velocity_converged = raw.get("velocity_converged", True)

        def verdict(ok, kind="physical"):
            # Calibration anchors are matched to the reference *by construction*
            # (the reference number is injected into the model), so a "PASS" there
            # would be circular — label it honestly instead. Predictive waveform
            # quantities are computed but absolute-approximate, so they are flagged
            # "predicted" rather than given a pass/fail they would often fail.
            if kind == "calibrated":
                return "by construction"
            if kind == "predicted":
                return "predicted (approx)"
            return "PASS" if ok else "CHECK"

        def row(name, calc, ref, ok, note="", kind="physical"):
            return f"| {name} | {calc} | {ref} | {verdict(ok, kind)} | {note} |\n"

        rows = []
        if topo.name == TOPOLOGY_RB87_TELECOM:
            pump_mw = max(params["pump_biphoton_uw"] * 1e-3, 1e-30)
            rate_per_mw = stats["pair_rate_cps"] / pump_mw
            rows.append(row(
                "Pair rate / pump", f"{rate_per_mw:.0f} cps/mW",
                "38000 cps/mW", abs(rate_per_mw / 38000.0 - 1.0) < 0.15,
                "anchored: pair rate is scaled from this reference number", kind="calibrated"))
            rows.append(row(
                "g2 peak", f"{stats['g2_peak']:.2f} (raw {stats['raw_g2_peak']:.1f})",
                "44(3)", abs(stats["g2_peak"] - 44.0) <= 3.0,
                ("predicted from waveform |ψ|² + physical accidentals"
                 if predictive else
                 "anchored: forced to the target via added-accidental calibration"),
                kind="predicted" if predictive else "calibrated"))
            if predictive and not velocity_converged:
                rows.append(row(
                    "Detected BTW FWHM", "unconverged",
                    "0.56(4) ns", False,
                    "velocity-grid auto-refine hit the point cap — absolute width "
                    "not converged, qualitative only (shape/ordering still indicative)",
                    kind="predicted"))
            else:
                rows.append(row(
                    "Detected BTW FWHM", f"{stats['detected_fwhm_ns']:.3f} ns",
                    "0.56(4) ns",
                    abs(stats["detected_fwhm_ns"] - 0.56) <= 0.04,
                    ("velocity-converged source waveform + net timing response; "
                     "absolute ns carries per-source calibration uncertainty"
                     if predictive else "model waveform + net timing response"),
                    kind="predicted" if predictive else "physical"))
            rows.append(row(
                "OD estimate", f"{raw['od_estimate']:.1f}",
                "112(3)", abs(raw["od_estimate"] - 112.0) <= 3.0,
                "anchored: density/cell scaling of the reference OD", kind="calibrated"))
            rows.append(row(
                "Bandwidth setting", f"{params['filter_bandwidth_mhz']:.0f} MHz",
                "about 300 MHz", abs(params["filter_bandwidth_mhz"] - 300.0) <= 40.0,
                "user filter-bandwidth setting check"))
        elif topo.name == TOPOLOGY_CS_BTW:
            other_channel = CS_CHANNEL_795 if params.get("cs_channel") == CS_CHANNEL_917 else CS_CHANNEL_917
            other_params = dict(params)
            other_params["cs_channel"] = other_channel
            other_raw = _solver_from_topology_key(
                _topology_cache_key(other_params)).compute_biphoton(other_params)
            other_stats = observables.biphoton_stats(
                other_raw["tau_axis_ns"], other_raw["psi_tau"], other_raw["pair_rate_cps"],
                signal_eff=params["signal_eff_pct"] / 100.0,
                idler_eff=params["idler_eff_pct"] / 100.0,
                dark_signal_cps=params["dark_signal_cps"],
                dark_idler_cps=params["dark_idler_cps"],
                coincidence_window_ns=params["coincidence_window_ns"],
                timing_jitter_ns=params["timing_jitter_ns"],
                filter_bandwidth_mhz=params["filter_bandwidth_mhz"],
                source_bandwidth_mhz=other_raw["source_bandwidth_mhz"],
                target_g2_peak=other_raw["topology"].target_g2_peak,
            )
            ratio = max(
                stats["detected_fwhm_ns"], other_stats["detected_fwhm_ns"]
            ) / max(
                min(stats["detected_fwhm_ns"], other_stats["detected_fwhm_ns"]),
                1e-30)
            rows.append(row(
                "Detected BTW width ratio", f"{ratio:.2f}",
                "about 3", abs(ratio - 3.0) <= 0.5,
                ("predicted ordering (917 narrower than 795); absolute ratio approximate"
                 if predictive else
                 "medium model only; full Cs BTW theory is not yet included"),
                kind="predicted" if predictive else "physical"))
            rows.append(row(
                "OD estimate", f"{raw['od_estimate']:.1f}",
                "about 10", abs(raw["od_estimate"] - 10.0) <= 2.0,
                "anchored: density/cell scaling of the reference OD", kind="calibrated"))
        else:
            rows.append(row(
                "Reference validation", "generic diamond template",
                "no paper anchor", False,
                "configure wavelengths manually; no validated default"))

        rows.append(row(
            "Phase matching", f"{raw['phase_match_weight']:.3f}",
            "> 0.90", raw["phase_match_weight"] > 0.90,
            "vector sinc^2(|Delta k| L / 2)"))
        rows.append(row(
            "Energy conservation", f"{raw['energy_mismatch_hz']/1e6:.3f} MHz",
            "near 0 MHz", abs(raw["energy_mismatch_hz"]) < 1e6,
            "wavelength bookkeeping"))

        if predictive:
            intro = (
                "**Reduced model.** Waveform-derived widths are approximate. Pair-"
                "rate rows are reference-anchored; phase matching and energy "
                "conservation are independent geometry checks.\n\n")
        else:
            intro = (
                "**Reference model.** Calibrated rows agree by construction. Phase "
                "matching and energy conservation are independent geometry checks.\n\n")
        return (
            intro
            + "| Check | Calculated | Reference | Verdict | Note |\n|---|---:|---:|---|---|\n"
            + "".join(rows)
        )

    def extra_views(self):
        def _compute_full(params):
            fidelity = normalize_fidelity(params.get("resolution", FIDELITY_FAST))
            fidelity_settings = FWM_FIDELITY[fidelity]
            return full_spectrum(
                params["opd"], params["temp_c"] + 273.15,
                params["pump_mw"], params["probe_uw"], params["line_strength"],
                params["loss_pct"], L=params["cell_mm"] * 1e-3,
                mode_overlap_penalty=params.get("mode_overlap_penalty", 1.0),
                polarization_penalty=params.get("polarization_penalty", 1.0),
                zeeman_participation_penalty=params.get(
                    "zeeman_participation_penalty", 1.0),
                w_pump=params.get("pump_waist_um", W_PUMP * 1e6) * 1e-6,
                w_probe=params.get("probe_waist_um", W_PROBE * 1e6) * 1e-6,
                qe=params.get("qe_pct", QE_DETECTOR * 100.0) / 100.0,
                detection_efficiency=_detection_efficiency_from_params(params),
                floquet_order=params.get("floquet_order", SEEDED_FLOQUET_ORDER),
                phase_detail=fidelity_settings["phase_detail"],
                pump_probe_angle_deg=params.get(
                    "seeded_angle_deg", SEEDED_PHASE_ANGLE_DEG),
                model_fidelity=fidelity,
                velocity_step=fidelity_settings["velocity_step"],
                velocity_cutoff=fidelity_settings.get("velocity_cutoff", 3.0),
                transit_rate=(2.0 * np.pi
                              * params.get("transit_rate_khz", 100.0) * 1e3),
                eom_residual_carrier_power=(
                    params.get("eom_residual_carrier_uw", 0.0) * 1e-6),
                eom_other_sidebands_power=(
                    params.get("eom_other_sidebands_uw", 0.0) * 1e-6),
                eom_seed_spectrum_provenance=params.get(
                    "eom_seed_spectrum_provenance", "direct FWM input"),
                eom_seed_spectrum_status=params.get(
                    "eom_seed_spectrum_status", "not supplied"),
                eom_seed_spectrum_application=params.get(
                    "eom_seed_spectrum_application", "unapplied"))

        def _render_full(full):
            import matplotlib.pyplot as plt
            figF, (aG, aS) = plt.subplots(2, 1, figsize=(8.5, 6.4), sharex=True)
            for ax in (aG, aS):
                ax.grid(alpha=0.3)
            styles = {
                "minus": dict(color="#1f77b4", label="minus Raman branch"),
                "plus": dict(color="#ff7f0e", label="plus Raman branch"),
            }
            for key, style in styles.items():
                spec = full[key]
                aG.plot(spec["probe_axis_GHz"], spec["G_s"], lw=1.4, **style)
                aS.plot(spec["probe_axis_GHz"],
                        spec["gain_referred_noise_dB"], lw=1.4, **style)
            aG.axhline(1.0, color="black", lw=0.6)
            aG.set_ylabel("Seed gain G_s")
            stored_detail = full["minus"].get("model_fidelity", "unknown")
            aG.set_title(
                f"Solver detail: {FIDELITY_LABELS.get(stored_detail, stored_detail)}")
            if max(np.nanmax(full[key]["G_s"]) for key in styles) > 50:
                aG.set_yscale("log")
            aS.axhline(0.0, color="black", lw=0.6)
            aS.set_ylabel("Squeezing indicator [dB]")
            aS.set_xlabel("Probe detuning from F=2 -> F'=3 [GHz]")
            for a in (aG, aS):
                a.axvline(full["D_GHz"], color="gray", ls=":", lw=0.8)
            aG.legend(fontsize=9)
            aS.legend(fontsize=9)
            figF.tight_layout()
            return figF

        return [ExtraView(
            key="Full probe scan (−8 to 12 GHz, slow)",
            description="Runs both Raman branches using the selected geometry and "
                        "solver detail.",
            compute=_compute_full, render=_render_full,
        )]


def full_spectrum(D_GHz, T_K, P_pump_mW, P_probe_uW, line_strength, loss_pct,
                  L=L_CELL, *, mode_overlap_penalty=1.0,
                  polarization_penalty=1.0,
                   zeeman_participation_penalty=1.0,
                   w_pump=W_PUMP, w_probe=W_PROBE, qe=QE_DETECTOR,
                   detection_efficiency=None,
                   floquet_order=SEEDED_FLOQUET_ORDER,
                   phase_detail=PHASE_BALANCED,
                   pump_probe_angle_deg=SEEDED_PHASE_ANGLE_DEG,
                   model_fidelity=FIDELITY_FAST,
                   velocity_step=4.0, velocity_cutoff=3.0,
                   transit_rate=constants.GAMMA_GG,
                  eom_residual_carrier_power=0.0,
                  eom_other_sidebands_power=0.0,
                  eom_seed_spectrum_provenance="direct FWM input",
                  eom_seed_spectrum_status="not supplied",
                  eom_seed_spectrum_application="unapplied"):
    """Wide scan with the two Raman channels calculated independently.

    The minus/plus branches are independent. They are evaluated sequentially
    because each compiled Floquet kernel already occupies the Numba worker pool;
    wrapping both in Python threads oversubscribes that pool and is slower on the
    benchmark hardware. This scheduling choice is not a numerical approximation.
    """
    common = dict(
        T=T_K, P_pump=P_pump_mW * 1e-3, P_probe=P_probe_uW * 1e-6,
        line_strength=line_strength, mode_overlap_penalty=mode_overlap_penalty,
        polarization_penalty=polarization_penalty,
        zeeman_participation_penalty=zeeman_participation_penalty,
        L=L, loss_frac=loss_pct / 100.0,
        w_pump=w_pump, w_probe=w_probe, qe=qe,
        detection_efficiency=detection_efficiency,
        floquet_order=floquet_order, enforce_floquet_convergence=True,
        phase_detail=phase_detail,
        pump_probe_angle_deg=pump_probe_angle_deg,
        model_fidelity=model_fidelity,
        transit_rate=transit_rate,
        eom_residual_carrier_power=eom_residual_carrier_power,
        eom_other_sidebands_power=eom_other_sidebands_power,
        eom_seed_spectrum_provenance=eom_seed_spectrum_provenance,
        eom_seed_spectrum_status=eom_seed_spectrum_status,
        eom_seed_spectrum_application=eom_seed_spectrum_application,
        coarse_points=301, fine_points=401,
        velocity_step=velocity_step, velocity_cutoff=velocity_cutoff)

    def _branch(b):
        with blas_single_thread():
            return compute_spectrum(D_GHz, branch=b, **common)

    minus = _branch(-1)
    plus = _branch(+1)
    return {"D_GHz": D_GHz, "minus": minus, "plus": plus}
