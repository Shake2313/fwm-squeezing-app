"""Detection-phase scans for the effective, real Jasperse FWM model.

This readout completes the model's amplitude covariance with shot-noise seed
phase fluctuations. It is an explicitly assumed Gaussian model, independent
of microscopic atomic-noise validation or an RF sideband calculation.
"""

import numpy as np

from .jasperse import jasperse_distributed_noise


_SCAN_ENDPOINT_TOL_MHZ = 1e-9


def experimental_quadrature_phase_scan(
        readout, scan_axis, selected_delta_mhz, *, phase_degrees=None):
    """Return two orthogonal equal-weight collective-quadrature noise traces.

    The mode ordering is (X_probe, P_probe, X_conjugate, P_conjugate), with
    vacuum covariance I and [X, P] = 2i. Define Q=(X_probe-X_conjugate)/sqrt(2)
    and P=(P_probe-P_conjugate)/sqrt(2). Trace A measures Q*cos(theta)+P*sin(theta);
    trace B uses theta+90 degrees. Both local oscillators rotate by theta.
    These homodyne observables differ from bright intensity-difference noise.

    ``readout`` is an experimental_jasperse_readout result. Its amplitude
    covariance retains the selected seed RIN, whereas seed phase variance is
    explicitly assumed to be one and all input X-P correlations are zero.
    The same real mixing and distributed losses supply the phase covariance.
    Downstream loss is applied once. Full covariances, not dB or squeezing
    parameters, are linearly interpolated between adjacent detuning samples.
    A selected point outside ``scan_axis`` returns an unavailable result.
    Endpoint differences within 1e-9 MHz are treated as representation roundoff;
    interpolation uses that endpoint while metadata retains the requested point.
    """
    axis = np.asarray(scan_axis, dtype=float)
    phase = (np.linspace(0.0, 180.0, 361) if phase_degrees is None
             else np.asarray(phase_degrees, dtype=float))
    if (axis.ndim != 1 or not axis.size or not np.all(np.isfinite(axis))
            or np.any(np.diff(axis) <= 0)):
        raise ValueError("scan_axis must be a nonempty, finite, strictly increasing vector")
    if phase.ndim != 1 or not phase.size or not np.all(np.isfinite(phase)):
        raise ValueError("phase_degrees must be a nonempty finite vector")
    selected = float(selected_delta_mhz)
    if not np.isfinite(selected):
        raise ValueError("selected_delta_mhz must be finite")
    result = {
        "available": False,
        "status": "OPERATING_POINT_OUTSIDE_SCAN",
        "phase_degrees": phase.copy(),
        "noise_a": np.full(phase.shape, np.nan),
        "noise_b": np.full(phase.shape, np.nan),
        "noise_a_db": np.full(phase.shape, np.nan),
        "noise_b_db": np.full(phase.shape, np.nan),
        "covariance": np.full((4, 4), np.nan),
        "minimum_variance": np.nan,
        "maximum_variance": np.nan,
        "minimum_db": np.nan,
        "maximum_db": np.nan,
        "crossing_variance": np.nan,
        "crossing_db": np.nan,
        "crossing_phase_degrees": np.array([45.0, 135.0]),
        "selected_delta_mhz": selected,
        "scan_support_mhz": (float(axis[0]), float(axis[-1])),
        "minimum_uncertainty_eigenvalue": np.nan,
        "seed_phase_assumption": "shot-noise seed phase variance = 1; no input X-P correlation",
        "convention": (
            "vacuum = 1; interleaved (X_probe,P_probe,X_conjugate,P_conjugate); "
            "A=(X_probe(theta)-X_conjugate(theta))/sqrt(2); B=A(theta+90 degrees); "
            "theta rotates both local-oscillator phases relative to the model's real coupling"),
        "scope": (
            "Effective real phase-matched Jasperse Gaussian model; amplitude seed RIN retained, "
            "seed phase shot noise assumed; both distributed vacuum-loss reservoirs included; "
            "after-cell efficiency applied once; no microscopic atomic or finite-RF prediction. "
            "Detection-phase rotation leaves source power gains fixed."),
    }
    if (selected < axis[0]-_SCAN_ENDPOINT_TOL_MHZ
            or selected > axis[-1]+_SCAN_ENDPOINT_TOL_MHZ):
        return result

    # GHz subtraction followed by MHz conversion can move a nominal endpoint
    # inward by a few ulps. Snap only differences within the endpoint tolerance.
    interpolation_point = float(np.clip(selected, axis[0], axis[-1]))
    index = int(np.searchsorted(axis, interpolation_point))
    if index < axis.size and axis[index] == interpolation_point:
        indices, weights = np.array([index]), np.array([1.0])
    else:
        fraction = (interpolation_point-axis[index-1])/(axis[index]-axis[index-1])
        indices, weights = np.array([index-1, index]), np.array([1-fraction, fraction])

    def samples(key, trailing=()):
        values = np.asarray(readout[key], dtype=float)
        return np.broadcast_to(values, (axis.size,) + trailing)[indices]

    try:
        amplitude = samples("quadrature_covariance", (2, 2))
        squeeze = samples("effective_squeeze")
        probe_transmission = samples("probe_transmission")
        conjugate_transmission = samples("conjugate_transmission")
        efficiency = samples("eta_after_cell_QE")
    except KeyError:
        result["status"] = "EFFECTIVE_QUADRATURE_INPUT_UNAVAILABLE"
        return result
    if (not all(np.all(np.isfinite(v)) for v in (
            amplitude, squeeze, probe_transmission, conjugate_transmission, efficiency))
            or np.any((efficiency < 0) | (efficiency > 1))):
        result["status"] = "INVALID_EFFECTIVE_QUADRATURE_INPUT"
        return result

    reference = jasperse_distributed_noise(
        squeeze, probe_transmission, conjugate_transmission, seed_fano=1.0)
    # The phase drift has the opposite mixing sign. Conjugation by diag(1,-1)
    # changes that sign while preserving diagonal loss and vacuum injection.
    signs = np.array([1.0, -1.0])
    phase_covariance = (reference["quadrature_covariance"]
                        * signs[:, None] * signs[None, :])
    covariances = np.zeros((indices.size, 4, 4))
    covariances[:, ::2, ::2] = amplitude
    covariances[:, 1::2, 1::2] = phase_covariance
    covariances = (efficiency[:, None, None] * covariances
                   + (1-efficiency[:, None, None]) * np.eye(4))
    commutator = np.kron(np.eye(2), [[0.0, 1.0], [-1.0, 0.0]])
    if (not np.all(np.isfinite(covariances))
            or not np.allclose(covariances, covariances.swapaxes(-1, -2),
                               rtol=1e-12, atol=1e-12)):
        result["status"] = "INVALID_EFFECTIVE_QUADRATURE_COVARIANCE"
        return result
    uncertainty = np.linalg.eigvalsh(covariances + 1j*commutator)[:, 0]
    tolerance = 1e-10 * np.maximum(1.0, np.linalg.norm(covariances, axis=(-2, -1)))
    if np.any(uncertainty < -tolerance):
        result["status"] = "EFFECTIVE_QUADRATURE_UNCERTAINTY_FAILED"
        return result
    covariance = np.einsum("n,nij->ij", weights, covariances)
    q = np.array([1.0, 0.0, -1.0, 0.0])/np.sqrt(2)
    p = np.array([0.0, 1.0, 0.0, -1.0])/np.sqrt(2)
    q_variance, p_variance = float(q@covariance@q), float(p@covariance@p)
    if min(q_variance, p_variance) <= 0:
        result["status"] = "INVALID_EFFECTIVE_QUADRATURE_VARIANCE"
        return result
    theta = np.deg2rad(phase)
    cosine_squared, sine_squared = np.cos(theta)**2, np.sin(theta)**2
    noise_a = q_variance*cosine_squared + p_variance*sine_squared
    noise_b = q_variance*sine_squared + p_variance*cosine_squared
    minimum, maximum = min(q_variance, p_variance), max(q_variance, p_variance)
    crossing = 0.5*(q_variance+p_variance)
    result.update(
        available=True, status="EFFECTIVE_JASPERSE_QUADRATURE_ESTIMATE",
        covariance=covariance, noise_a=noise_a, noise_b=noise_b,
        noise_a_db=10*np.log10(noise_a), noise_b_db=10*np.log10(noise_b),
        minimum_variance=minimum, maximum_variance=maximum,
        minimum_db=float(10*np.log10(minimum)), maximum_db=float(10*np.log10(maximum)),
        crossing_variance=crossing, crossing_db=float(10*np.log10(crossing)),
        minimum_uncertainty_eigenvalue=float(np.linalg.eigvalsh(covariance+1j*commutator)[0]))
    return result
