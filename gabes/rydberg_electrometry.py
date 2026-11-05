"""Weak-signal Rydberg electrometry and detector-noise primitives.

This module supplies the pieces between a static, LO-dressed optical Bloch
equation and an absolute RF sensitivity.  It deliberately does *not* implement
full time-domain integration, a lock-in filter, or a Zeeman-resolved Rydberg
manifold.  Those require experiment-specific choices that do not belong in a
generic linear-response primitive.

Conventions
-----------
All Hamiltonians and Rabi frequencies are angular frequencies in rad/s.  In the
frame rotating at the microwave LO, a real weak signal is written

    H_sig(t) = Omega_sig / 2 * [
        exp(-i (omega_if t + phi)) |u><l|
        + exp(+i (omega_if t + phi)) |l><u|],

where ``Omega_sig`` is a *peak* angular Rabi-frequency amplitude.  The returned
density-matrix coefficients obey

    rho(t) = rho_0 + Omega_sig * [
        rho_minus exp(-i omega_if t) + rho_plus exp(+i omega_if t)]
        + O(Omega_sig**2).

Thus the sideband coefficients have units s/(rad), conventionally written s
because the radian is dimensionless.  For a Hermitian observable O, the method
``real_observable_phasor_per_angular_rabi`` returns Q such that

    delta <O>(t) = Re[Omega_sig * Q * exp(-i omega_if t)].

Photocurrent noise amplitudes are one-sided amplitude spectral densities (ASD)
in A/sqrt(Hz).  Electric-field sensitivities are therefore one-sided ASDs in
V/m/sqrt(Hz).  No paper/reference sensitivity is injected: transition dipoles,
optical powers, detector parameters, and readout responsivity must be supplied
by the caller.

The photodiode, balanced-detector and shot-noise model is shared by every scheme
and lives in :mod:`gabes.detection`; it is re-exported here under its original
names.
"""

from dataclasses import dataclass

import numpy as np

from . import constants, core
from .detection import (
    BalancedDetector,
    DetectorNoiseBudget,
    PhotodiodeChannel,
    _require_finite_nonnegative,
    _require_finite_positive,
    asd_to_rms,
    balanced_detector_noise,
    photodiode_responsivity_a_per_w,
)


def _validate_field_convention(convention):
    convention = str(convention).lower()
    if convention not in ("peak", "rms"):
        raise ValueError("field_amplitude_convention must be 'peak' or 'rms'")
    return convention


@dataclass(frozen=True)
class WeakSignalResponse:
    """First-order density-matrix response to a finite-IF RF signal.

    Arrays may carry leading batch dimensions.  This permits callers to solve
    separate velocity classes or spatial segments and average their *complex
    phasors* afterwards.  Averaging magnitudes would discard phase and is not
    equivalent.
    """

    if_angular_frequency_rad_s: float
    signal_transition: tuple
    signal_phase_rad: float
    steady_state: np.ndarray
    rho_minus_per_angular_rabi: np.ndarray
    rho_plus_per_angular_rabi: np.ndarray

    def real_observable_phasor_per_angular_rabi(self, operator):
        """Return the peak phasor of a Hermitian observable per rad/s of SIG.

        The result has the same leading batch shape as the response arrays.  It
        is complex: its magnitude is peak response and its argument is phase
        relative to ``exp(-i omega_if t)``.
        """
        operator = np.asarray(operator, dtype=complex)
        n = self.steady_state.shape[-1]
        if operator.shape != (n, n):
            raise ValueError(f"operator must have shape {(n, n)}")
        if not np.allclose(operator, operator.conj().T, rtol=1e-10, atol=1e-12):
            raise ValueError("operator must be Hermitian for a real readout")
        # Tr(O rho) = sum_ij O_ij rho_ji.  Symmetrising both sidebands makes the
        # real-observable convention robust to floating-point Hermiticity drift.
        y_minus = np.einsum(
            "ij,...ji->...", operator, self.rho_minus_per_angular_rabi)
        y_plus = np.einsum(
            "ij,...ji->...", operator, self.rho_plus_per_angular_rabi)
        return y_minus + np.conj(y_plus)

    def coherence_sidebands_per_angular_rabi(self, row, column):
        """Return ``(minus, plus)`` coefficients for one rho[row, column]."""
        n = self.steady_state.shape[-1]
        if not (0 <= int(row) < n and 0 <= int(column) < n):
            raise IndexError("coherence indices are outside the Hilbert space")
        return (self.rho_minus_per_angular_rabi[..., int(row), int(column)],
                self.rho_plus_per_angular_rabi[..., int(row), int(column)])


def weak_signal_response_from_liouvillian(
        liouvillian, n_levels, if_angular_frequency_rad_s, *,
        signal_transition, signal_phase_rad=0.0, steady_state=None):
    """Solve the finite-IF weak-signal response around a dressed steady state.

    Parameters
    ----------
    liouvillian:
        LO-dressed, time-independent Liouvillian in GABES row-major vec
        convention, shape ``(..., n_levels**2, n_levels**2)`` and units 1/s.
    n_levels:
        Hilbert-space dimension.
    if_angular_frequency_rad_s:
        Positive IF angular frequency.  DC is intentionally rejected because
        the zero-frequency Liouvillian requires a trace-constrained derivative,
        while this API represents a heterodyne sideband.
    signal_transition:
        ``(lower, upper)`` level indices for the RF signal coupling.
    signal_phase_rad:
        SIG phase relative to the LO under the convention in the module docs.
    steady_state:
        Optional precomputed rho_0 with shape ``(..., n_levels, n_levels)``.

    Returns
    -------
    WeakSignalResponse
        Sideband density matrices per unit *peak angular Rabi frequency*.
    """
    n_levels = int(n_levels)
    if n_levels < 2:
        raise ValueError("n_levels must be at least 2")
    omega_if = _require_finite_positive(
        "if_angular_frequency_rad_s", if_angular_frequency_rad_s)
    phase = float(signal_phase_rad)
    if not np.isfinite(phase):
        raise ValueError("signal_phase_rad must be finite")

    try:
        lower, upper = (int(signal_transition[0]), int(signal_transition[1]))
    except (TypeError, ValueError, IndexError) as exc:
        raise ValueError("signal_transition must be a (lower, upper) pair") from exc
    if lower == upper or not (0 <= lower < n_levels and 0 <= upper < n_levels):
        raise ValueError("signal_transition must name two distinct valid levels")

    L0 = np.asarray(liouvillian, dtype=complex)
    m = n_levels * n_levels
    if L0.shape[-2:] != (m, m):
        raise ValueError(f"liouvillian must end in shape {(m, m)}")

    if steady_state is None:
        rho0 = core.steady_state_from_liouvillian(L0, n_levels)
    else:
        rho0 = np.asarray(steady_state, dtype=complex)
        expected = L0.shape[:-2] + (n_levels, n_levels)
        if rho0.shape != expected:
            raise ValueError(f"steady_state must have shape {expected}")

    # Coefficient matrices per unit Omega_sig.  At IF -> 0 their sum is the
    # usual static RWA coupling Omega_sig/2 (|u><l| + |l><u|) for phase zero.
    h_minus = np.zeros((n_levels, n_levels), dtype=complex)
    h_minus[upper, lower] = 0.5 * np.exp(-1j * phase)
    c_minus = core.comm_super(h_minus)

    rho0_vec = rho0.reshape(rho0.shape[:-2] + (m,))
    rhs_minus = -np.einsum("ij,...j->...i", c_minus, rho0_vec)
    eye = np.eye(m, dtype=complex)
    a_minus = L0 + 1j * omega_if * eye
    r_minus = np.linalg.solve(a_minus, rhs_minus[..., None])[..., 0]
    rho_minus = r_minus.reshape(r_minus.shape[:-1] + (n_levels, n_levels))
    rho_plus = np.swapaxes(rho_minus.conj(), -1, -2)

    return WeakSignalResponse(
        if_angular_frequency_rad_s=omega_if,
        signal_transition=(lower, upper),
        signal_phase_rad=phase,
        steady_state=rho0,
        rho_minus_per_angular_rabi=rho_minus,
        rho_plus_per_angular_rabi=rho_plus,
    )


def weak_signal_response(atom, hamiltonian_rad_s, if_angular_frequency_rad_s,
                         *, signal_transition, signal_phase_rad=0.0):
    """Convenience wrapper building the LO-dressed Liouvillian from ``atom``."""
    hamiltonian = np.asarray(hamiltonian_rad_s, dtype=complex)
    expected = (atom.n_levels, atom.n_levels)
    if hamiltonian.shape != expected:
        raise ValueError(f"hamiltonian_rad_s must have shape {expected}")
    return weak_signal_response_from_liouvillian(
        core.build_liouvillian(hamiltonian, atom), atom.n_levels,
        if_angular_frequency_rad_s,
        signal_transition=signal_transition,
        signal_phase_rad=signal_phase_rad)


def coherent_weighted_average(values, weights, axis=0):
    """Average complex response phasors with normalized non-negative weights."""
    values = np.asarray(values)
    weights = np.asarray(weights, dtype=float)
    if weights.ndim != 1:
        raise ValueError("weights must be one-dimensional")
    if not np.all(np.isfinite(weights)) or np.any(weights < 0.0):
        raise ValueError("weights must be finite and non-negative")
    total = weights.sum()
    if total <= 0.0:
        raise ValueError("weights must have a positive sum")
    axis = int(axis)
    if values.shape[axis] != weights.size:
        raise ValueError("weights length must match the selected values axis")
    shape = [1] * values.ndim
    shape[axis] = weights.size
    return np.sum(values * (weights / total).reshape(shape), axis=axis)


@dataclass(frozen=True)
class RFDipoleCoupling:
    """SI conversion between RF electric field and angular Rabi frequency.

    ``transition_dipole_c_m`` is the magnitude of the relevant matrix element in
    C*m.  ``angular_factor`` carries Clebsch-Gordan, polarization projection, or
    other explicitly justified dimensionless participation.  No Rydberg matrix
    element is assumed by this class.

    With ``field_amplitude_convention='peak'``, Omega = d_eff E_peak / hbar.
    With ``'rms'``, the field argument is E_rms and Omega_peak = sqrt(2) d_eff
    E_rms / hbar.  Rabi values returned by this class are always peak amplitudes.
    """

    transition_dipole_c_m: float
    angular_factor: float = 1.0
    field_amplitude_convention: str = "peak"

    def __post_init__(self):
        _require_finite_positive(
            "transition_dipole_c_m", self.transition_dipole_c_m)
        _require_finite_positive("angular_factor", self.angular_factor)
        _validate_field_convention(self.field_amplitude_convention)

    @property
    def effective_dipole_c_m(self):
        return float(self.transition_dipole_c_m) * float(self.angular_factor)

    @property
    def angular_rabi_per_field_rad_s_per_v_m(self):
        peak_per_input = (np.sqrt(2.0)
                          if self.field_amplitude_convention.lower() == "rms"
                          else 1.0)
        return peak_per_input * self.effective_dipole_c_m / constants.HBAR

    def angular_rabi_from_field(self, electric_field_v_m):
        """Peak angular Rabi amplitude [rad/s] from field [V/m]."""
        field = np.asarray(electric_field_v_m, dtype=float)
        if np.any(~np.isfinite(field)) or np.any(field < 0.0):
            raise ValueError("electric_field_v_m must be finite and non-negative")
        out = field * self.angular_rabi_per_field_rad_s_per_v_m
        return float(out) if out.ndim == 0 else out

    def cyclic_rabi_from_field_hz(self, electric_field_v_m):
        """Peak cyclic Rabi amplitude Omega/(2*pi) [Hz] from field [V/m]."""
        return self.angular_rabi_from_field(electric_field_v_m) / (2.0 * np.pi)

    def field_from_angular_rabi(self, angular_rabi_rad_s):
        """Field [V/m] from a peak angular Rabi amplitude [rad/s]."""
        rabi = np.asarray(angular_rabi_rad_s, dtype=float)
        if np.any(~np.isfinite(rabi)) or np.any(rabi < 0.0):
            raise ValueError("angular_rabi_rad_s must be finite and non-negative")
        out = rabi / self.angular_rabi_per_field_rad_s_per_v_m
        return float(out) if out.ndim == 0 else out

    def field_from_cyclic_rabi_hz(self, cyclic_rabi_hz):
        """Field [V/m] from a peak cyclic Rabi amplitude Omega/(2*pi) [Hz]."""
        freq = np.asarray(cyclic_rabi_hz, dtype=float)
        return self.field_from_angular_rabi(2.0 * np.pi * freq)

    def field_from_at_splitting_hz(self, splitting_hz, *, detuning_hz=0.0):
        """Ideal two-level AT splitting estimate of field [V/m].

        The generalized dressed-state splitting is
        ``sqrt(Omega^2 + Delta^2)/(2*pi)``.  This helper removes a known RF
        detuning before applying the dipole calibration.  Optical line pulling,
        unresolved peaks, and multilevel structure are not corrected; those
        require fitting the complete optical spectrum with ``Omega`` as a model
        parameter.
        """
        splitting = np.asarray(splitting_hz, dtype=float)
        detuning = np.asarray(detuning_hz, dtype=float)
        try:
            splitting, detuning = np.broadcast_arrays(splitting, detuning)
        except ValueError as exc:
            raise ValueError("splitting_hz and detuning_hz are not broadcastable") from exc
        if (
            np.any(~np.isfinite(splitting))
            or np.any(~np.isfinite(detuning))
            or np.any(splitting < 0.0)
            or np.any(splitting + np.finfo(float).eps < np.abs(detuning))
        ):
            raise ValueError(
                "splitting_hz must be finite, non-negative, and at least "
                "abs(detuning_hz)")
        rabi_hz = np.sqrt(np.maximum(splitting**2 - detuning**2, 0.0))
        return self.field_from_cyclic_rabi_hz(rabi_hz)


@dataclass(frozen=True)
class ElectrometrySensitivity:
    """Noise-equivalent RF electric-field ASDs for one calibrated readout."""

    current_responsivity_a_per_v_m: float
    psn_field_asd_v_m_per_sqrt_hz: float
    technical_field_asd_v_m_per_sqrt_hz: float
    total_field_asd_v_m_per_sqrt_hz: float

    @property
    def psn_field_asd_nv_cm_per_sqrt_hz(self):
        return self.psn_field_asd_v_m_per_sqrt_hz * 1.0e7

    @property
    def technical_field_asd_nv_cm_per_sqrt_hz(self):
        return self.technical_field_asd_v_m_per_sqrt_hz * 1.0e7

    @property
    def total_field_asd_nv_cm_per_sqrt_hz(self):
        return self.total_field_asd_v_m_per_sqrt_hz * 1.0e7


def electrometry_sensitivity(noise_budget, current_responsivity_a_per_v_m):
    """Convert a detector-noise budget into PSN/technical/total field ASD.

    ``current_responsivity_a_per_v_m`` is the magnitude of the demodulated
    *peak-current phasor* per field amplitude using the same peak/RMS convention
    as the RF calibration.  RBW or lock-in gain must not be folded into this
    quantity; ASD normalisation is per sqrt(Hz).
    """
    if not isinstance(noise_budget, DetectorNoiseBudget):
        raise TypeError("noise_budget must be a DetectorNoiseBudget")
    responsivity = _require_finite_positive(
        "current_responsivity_a_per_v_m", current_responsivity_a_per_v_m)
    return ElectrometrySensitivity(
        current_responsivity_a_per_v_m=responsivity,
        psn_field_asd_v_m_per_sqrt_hz=(
            noise_budget.shot_noise_current_asd_a_per_sqrt_hz / responsivity),
        technical_field_asd_v_m_per_sqrt_hz=(
            noise_budget.technical_noise_current_asd_a_per_sqrt_hz / responsivity),
        total_field_asd_v_m_per_sqrt_hz=(
            noise_budget.total_noise_current_asd_a_per_sqrt_hz / responsivity),
    )


def current_responsivity_from_atomic_phasor(
        atomic_phasor_per_angular_rabi,
        current_per_atomic_readout_a,
        rf_coupling):
    """Compose atomic, optical/detector, and RF-field linear calibrations.

    Returns a complex current responsivity in A/(V/m), with current and field in
    the *same* amplitude convention.  For ``peak`` this is peak current per peak
    field; for ``rms`` the peak atomic phasor is divided by sqrt(2), giving RMS
    current per RMS field.  The two ratios are numerically identical, as they
    must be, so changing a reporting convention cannot change the predicted
    sensitivity.

    The caller supplies the optical propagation/detector slope
    ``current_per_atomic_readout_a``; this prevents a static OBE coherence from
    being mistaken for an absolute detector voltage or current calibration.
    """
    if not isinstance(rf_coupling, RFDipoleCoupling):
        raise TypeError("rf_coupling must be an RFDipoleCoupling")
    atomic = np.asarray(atomic_phasor_per_angular_rabi, dtype=complex)
    current_scale = complex(current_per_atomic_readout_a)
    if not np.isfinite(current_scale.real) or not np.isfinite(current_scale.imag):
        raise ValueError("current_per_atomic_readout_a must be finite")
    output_amplitude_factor = (
        1.0 / np.sqrt(2.0)
        if rf_coupling.field_amplitude_convention.lower() == "rms" else 1.0)
    out = (atomic * current_scale * output_amplitude_factor
           * rf_coupling.angular_rabi_per_field_rad_s_per_v_m)
    return complex(out) if out.ndim == 0 else out
