"""Conditional pump-dressed diagonal transport at the atomic-cell boundaries.

The mean-field Maxwell drift includes coherent mode mixing and diagonal
absorption/gain.  Turning off mixing while keeping that same pumped atomic
state defines a counterfactual diagonal intensity transfer.  A passive net
transfer can be displayed as an effective transmission, but this does not
separate microscopic absorption and stimulated-emission reservoirs.
"""

import numpy as np


def dressed_cell_transmission(maxwell_drift, L):
    """Report in-cell diagonal net transport and available passive efficiencies.

    ``maxwell_drift`` must have trailing shape ``(2, 2)`` in the field basis
    ``(E_probe, E_conjugate.conj())`` and units m^-1.  Its susceptibilities,
    optical frequencies and density normalization must be exactly those used
    for propagation.  ``L`` is the nonnegative scalar atomic-medium length in
    metres.  Windows and detector optics are excluded.

    With mode mixing disabled, ``log(T_j) = 2*Re(A_jj)*L = -alpha_j*L``.
    Raw transmissions retain active gain, including infinity on overflow;
    signed logarithms preserve its direction.  Passive efficiency is NaN
    where the diagonal transfer has net gain.  Only a positive log(T) within
    the reported floating-point tolerance is normalized to unity in that
    *efficiency* field; the raw transmission is unchanged.  Underflow to zero
    is a valid passive transmission.  Off-diagonal entries are not used.

    These are conditional net transmissions, not a microscopic loss-channel
    decomposition.  They must not be applied as a second loss after a full
    transfer matrix that already includes the same diagonal susceptibility.
    """
    drift = np.asarray(maxwell_drift, dtype=complex)
    if drift.ndim < 2 or drift.shape[-2:] != (2, 2):
        raise ValueError("maxwell_drift must end in a 2x2 matrix")
    if np.any(~np.isfinite(drift)):
        raise ValueError("maxwell_drift must be finite")
    length_array = np.asarray(L, dtype=float)
    if length_array.ndim != 0:
        raise ValueError("L must be a finite nonnegative scalar in metres")
    length = float(length_array)
    if not np.isfinite(length) or length < 0.0:
        raise ValueError("L must be a finite nonnegative scalar in metres")

    diagonal = np.real(np.diagonal(drift, axis1=-2, axis2=-1))
    with np.errstate(over="ignore", under="ignore", invalid="ignore"):
        # Multiplication order keeps L=0 exactly transparent for any finite A.
        log_transmission = (diagonal * length) * 2.0
        alpha = -2.0 * diagonal
        transmission = np.exp(log_transmission)
    # Scale the allowance by the computed logarithm.  Infinite active log(T)
    # must never be mistaken for a roundoff-sized deviation.
    finite_magnitude = np.where(
        np.isfinite(log_transmission), np.abs(log_transmission), 1.0)
    tolerance = 64.0 * np.finfo(float).eps * np.maximum(1.0, finite_magnitude)
    passive = log_transmission <= tolerance
    roundoff_unity = passive & (log_transmission > 0.0)
    efficiency = np.where(
        passive, np.where(roundoff_unity, 1.0, transmission), np.nan)
    status = np.where(
        passive, "NET_PASSIVE_DIAGONAL_RESPONSE",
        "ACTIVE_DIAGONAL_RESPONSE_NO_PASSIVE_EFFICIENCY")
    status = np.where(roundoff_unity, "UNITY_WITH_ROUNDOFF_TOLERANCE", status)

    result = {
        "in_cell_transmission_scope": (
            "conditional pump-dressed net diagonal transmission; mode mixing "
            "disabled at the same pumped atomic state; not a microscopic "
            "absorption/emission decomposition"),
        "in_cell_transmission_plane": (
            "atomic medium entrance to exit; excludes windows and detector optics"),
        "in_cell_transmission_applied_to_transfer_again": False,
        "in_cell_transmission_length_m": length,
    }
    for index, arm in enumerate(("probe", "conjugate")):
        prefix = f"in_cell_{arm}"
        result[f"{prefix}_transmission"] = transmission[..., index]
        result[f"{prefix}_log_transmission"] = log_transmission[..., index]
        result[f"{prefix}_signed_alpha_per_m"] = alpha[..., index]
        result[f"{prefix}_passive_efficiency"] = efficiency[..., index]
        result[f"{prefix}_net_passive"] = passive[..., index]
        result[f"{prefix}_log_roundoff_tolerance"] = tolerance[..., index]
        result[f"{prefix}_roundoff_normalized"] = roundoff_unity[..., index]
        result[f"{prefix}_status"] = status[..., index]
    return result
