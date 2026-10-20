"""Conditional passive noise diagnostics for a bright seeded FWM transfer.

This module completes a canonical ``(a, b†)`` transfer with annihilation-only
vacuum reservoirs. It does not derive atomic Langevin diffusion. A positive
completion defect establishes endpoint feasibility, not a microscopic model.
"""

import numpy as np


def passive_l1_readout(transfer_canonical, *, seed_fano=1.0, eta=1.0,
                       eta_probe=None, eta_conjugate=None):
    """Return the bright intensity-difference diagnostic in linear SQL units.

    ``seed_fano`` is the seed amplitude-quadrature variance at the cell input
    (vacuum = 1); the conjugate input has zero coherent mean and is vacuum.
    Seed phase noise cancels in this bright linearization. ``eta`` is a shared
    downstream efficiency for the two photon-flux channels, excluding losses
    already represented by the transfer. All inputs may broadcast over modes.
    The same effective transfer must describe carrier and fluctuations; this
    is a static diagnostic, not an RF spectrum using independent sidebands.
    Single-arm Fano expressions require a bright coherent output in that arm.

    For an annihilation-only vacuum completion, the result is exact. Additional
    phase-insensitive reservoir noise adds a nonnegative term, so the formal
    expression is a conditional lower bound. No such claim applies to arbitrary
    phase-sensitive, correlated atomic reservoirs. Negative defect eigenvalues
    invalidate the passive completion; formal values are retained separately,
    and the display-ready dB is NaN. No noise or eigenvalue is clipped/repaired.
    """
    matrix = np.asarray(transfer_canonical, dtype=complex)
    if matrix.ndim < 2 or matrix.shape[-2:] != (2, 2):
        raise ValueError("transfer_canonical must end in a 2x2 matrix")
    if np.any(~np.isfinite(matrix)):
        raise ValueError("transfer_canonical must be finite")
    fano = np.asarray(seed_fano, dtype=float)
    efficiency = np.asarray(eta, dtype=float)
    if np.any(~np.isfinite(fano)) or np.any(fano <= 0.0):
        raise ValueError("seed_fano must be finite and positive")
    if (np.any(~np.isfinite(efficiency)) or np.any(efficiency < 0.0)
            or np.any(efficiency > 1.0)):
        raise ValueError("eta must be finite and between zero and one")
    shape = np.broadcast_shapes(matrix.shape[:-2], fano.shape, efficiency.shape)
    matrix = np.broadcast_to(matrix, shape + (2, 2))
    fano = np.broadcast_to(fano, shape)
    efficiency = np.broadcast_to(efficiency, shape)
    gs = np.abs(matrix[..., 0, 0]) ** 2
    gc = np.abs(matrix[..., 1, 0]) ** 2
    difference, total = gs - gc, gs + gc
    q1 = (matrix[..., 0, 0].conj() * matrix[..., 0, 1]
          - matrix[..., 1, 0].conj() * matrix[..., 1, 1])
    det_abs2 = np.abs(np.linalg.det(matrix)) ** 2
    adjoint = np.swapaxes(matrix.conj(), -1, -2)
    signature = np.diag([1.0, -1.0])
    defect = signature - matrix @ signature @ adjoint
    # Symmetrization removes roundoff in Hermitian arithmetic only. The raw
    # defect and its eigenvalues are never projected onto the positive cone.
    eigmin = np.linalg.eigvalsh(
        0.5 * (defect + np.swapaxes(defect.conj(), -1, -2)))[..., 0]
    # Backward-error allowance for the small matrix products/eigensolve. A
    # loose relative tolerance would hide real violations at high gain.
    tolerance = (64.0 * np.finfo(float).eps
                 * np.maximum(1.0, np.sum(np.abs(matrix) ** 2, axis=(-2, -1))))
    completion_valid = eigmin >= -tolerance
    numerator = difference + (fano - 1.0) * difference ** 2 + 2.0 * np.abs(q1) ** 2
    with np.errstate(divide="ignore", invalid="ignore"):
        source = np.where(total > 0.0, numerator / total, np.nan)
        detected = 1.0 - efficiency + efficiency * source
        formal_db = np.where(detected > 0.0, 10.0 * np.log10(detected), np.nan)
    physical_variance = np.isfinite(source) & (source >= 0.0)
    readable = completion_valid & physical_variance & (total > 0.0)
    detected_mean = efficiency * total
    status = np.where(completion_valid, "PASSIVE_LOWER_BOUND_ONLY",
                      "L1_COMPLETION_NOT_PASSIVE")
    status = np.where((total > 0.0) & physical_variance, status,
                      "L1_NO_BRIGHT_OUTPUT_OR_INVALID_VARIANCE")
    fp = (fano - 1.0) * gs + 2.0 * np.abs(matrix[..., 0, 1]) ** 2 + 1.0
    fc = (fano - 1.0) * gc + 2.0 * np.abs(matrix[..., 1, 1]) ** 2 - 1.0
    if eta_probe is not None or eta_conjugate is not None:
        ep = np.asarray(eta if eta_probe is None else eta_probe, dtype=float)
        ec = np.asarray(eta if eta_conjugate is None else eta_conjugate, dtype=float)
        for value in (ep, ec):
            if np.any(~np.isfinite(value)) or np.any((value < 0) | (value > 1)):
                raise ValueError("external arm efficiencies must be between zero and one")
        # Jasperse's eta_a, eta_b describe losses AFTER mixing. Cell losses
        # already represented in M must not be applied again here.
        covariance = ((fano - 1.0) * gs * gc
                      + 2.0 * np.real(matrix[..., 0, 0].conj()
                                      * matrix[..., 1, 0]
                                      * matrix[..., 0, 1]
                                      * matrix[..., 1, 1].conj()))
        detected_mean = ep * gs + ec * gc
        detected_variance = (
            ep**2 * gs * fp + ec**2 * gc * fc - 2 * ep * ec * covariance
            + ep * (1 - ep) * gs + ec * (1 - ec) * gc)
        with np.errstate(divide="ignore", invalid="ignore"):
            detected = detected_variance / detected_mean
            formal_db = np.where(detected > 0, 10 * np.log10(detected), np.nan)
        readable = readable & (detected_mean > 0) & (detected >= 0)
    readable = readable & (detected_mean > 0)
    status = np.where((total > 0) & (detected_mean <= 0),
                      "L1_NO_DETECTED_BRIGHT_OUTPUT", status)
    nf = np.where(readable, detected, np.nan)
    nf_db = np.where(readable, formal_db, np.nan)
    return {
        "NF": nf,
        "NF_dB": nf_db,
        "NF_formal": detected,
        "NF_source": np.where(readable, source, np.nan),
        "NF_definition": "Var(N_probe-N_conjugate)/mean(N_probe+N_conjugate)",
        "passive_l1_G_s": gs,
        "passive_l1_G_c": gc,
        "passive_l1_D": difference,
        "passive_l1_T": total,
        "passive_l1_q1": q1,
        "passive_l1_det_abs2": det_abs2,
        "passive_l1_S_src": source,
        "passive_l1_S_det": detected,
        "passive_l1_F_p": fp,
        "passive_l1_F_c": fc,
        "passive_l1_Lambda": source * total - difference * (fp - fc - 1.0),
        "passive_l1_formal_noise_dB": formal_db,
        "passive_l1_noise_dB": nf_db,
        "passive_completion_min_eig": eigmin,
        "passive_completion_tolerance": tolerance,
        "passive_completion_valid": completion_valid,
        "passive_l1_status": status,
        "seed_fano_input": fano,
        "passive_l1_gain_basis": "canonical photon flux; pre-cap/pre-output-calibration",
        "passive_l1_scope": (
            "static bright seed, conjugate vacuum, downstream arm efficiencies; "
            "conditional on phase-insensitive reservoirs; not microscopic atomic noise"),
    }


def passive_l1_from_spectrum(spectrum, *, seed_fano=1.0):
    """Evaluate cached small-signal M without using capped/calibrated gains.

    The transfer includes any applied coupling closure. Its output power fit
    and pump-energy cap do not supply a quantum transfer and are not used here.
    Legacy pump-scatter and user phenomenological noise terms are not included.
    """
    result = passive_l1_readout(
        spectrum["T_canonical_small_signal"], seed_fano=seed_fano,
        eta=spectrum["eta"])
    adaptive = "adaptive" in spectrum.get("response_estimator", {}).get("method", "")
    result["passive_l1_sampling_status"] = (
        "ADAPTIVE_L1_TOLERANCE_NOT_CERTIFIED" if adaptive
        else "CACHED_TRANSFER_SAMPLES; NOT_A_CONTINUOUS_SCAN_CERTIFICATE")
    return result
