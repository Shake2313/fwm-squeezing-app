"""Effective Jasperse noise for the experimental FWM tool.

The two gains do not identify three distributed-medium parameters.  We close
that ambiguity with a user-selected common loss, then put the remaining loss
in one arm.  This is an explicit experimental approximation, independent of
the microscopic Grand Challenge noise calculation.
"""

from functools import lru_cache

import numpy as np
from scipy.optimize import brentq

from . import constants
from .jasperse import jasperse_distributed_noise, seed_noise_components


MODEL_ID = "effective_jasperse_common_plus_asymmetric_loss_v1"
MIN_GAIN = 1e-24
MAX_LOG_DETERMINANT = 300.0


def _uniform_gains(squeeze, log_determinant, *, probe_loss):
    """Stable first-column gains for one-arm uniform loss.

    The generator has one diagonal -q and one zero.  Computing its exponential
    through the two eigenvalues avoids overflow of cosh(r) times exp(-q/2).
    """
    half_loss = 0.5 * log_determinant
    radius = np.hypot(squeeze, half_loss)
    if radius == 0.0:
        return 1.0, 0.0
    positive = np.exp(-half_loss + radius)
    negative = np.exp(-half_loss - radius)
    large_weight = 0.5 * (1.0 + half_loss / radius)
    small_weight = 0.5 * squeeze * squeeze / (radius * (radius + half_loss))
    probe = (positive * small_weight + negative * large_weight
             if probe_loss else
             positive * large_weight + negative * small_weight)
    conjugate = (positive - negative) * squeeze / (2.0 * radius)
    return probe * probe, conjugate * conjugate


def _bounded_loss_approximation(gp, gc, *, probe_loss):
    """Preserve the dominant target using the largest representable loss.

    An almost absent conjugate alongside probe amplification can require an
    arbitrarily large conjugate loss coefficient.  Replacing it with a
    lossless amplifier would invent a bright conjugate and optimistic NF.
    Instead retain the distributed reservoir and solve its mixing strength;
    report the remaining gain mismatch as an approximation.
    """
    target = gc if probe_loss else gp
    channel = 1 if probe_loss else 0

    def objective(squeeze):
        return _uniform_gains(squeeze, MAX_LOG_DETERMINANT,
                              probe_loss=probe_loss)[channel] - target

    upper = max(1.0, np.arcsinh(np.sqrt(gc)),
                np.arccosh(np.sqrt(max(gp, 1.0))))
    while objective(upper) < 0.0:
        upper *= 2.0
    return brentq(objective, 0.0, upper, xtol=1e-12)


def _log_transfer(u, v, log_determinant):
    """Symmetric log(M), knowing M00=u, M10=v and det(M)>0.

    The small eigenvalue is obtained from the determinant rather than by
    subtracting nearly equal eigenvalues.  Eigenvector weights likewise avoid
    cancellation for weak conjugate gain.
    """
    determinant = np.exp(-log_determinant)
    w = (determinant + v * v) / u
    radius = np.hypot(u - w, 2.0 * v)
    if radius == 0.0:
        return np.log(u), 0.0, np.log(w)
    large = 0.5 * (u + w + radius)
    log_large = np.log(large)
    log_small = -log_determinant - log_large
    if u >= w:
        plus_probe = 0.5 * (1.0 + (u - w) / radius)
        minus_probe = 2.0 * v * v / (radius * (radius + u - w))
    else:
        plus_probe = 2.0 * v * v / (radius * (radius + w - u))
        minus_probe = 0.5 * (1.0 + (w - u) / radius)
    diagonal_probe = plus_probe * log_large + minus_probe * log_small
    diagonal_conjugate = minus_probe * log_large + plus_probe * log_small
    mixing = (log_large - log_small) * v / radius
    return diagonal_probe, mixing, diagonal_conjugate


def _infer_point(probe_gain, conjugate_gain, common_transmission):
    gp = max(float(probe_gain) / common_transmission, MIN_GAIN)
    gc = max(float(conjugate_gain) / common_transmission, 0.0)
    if gc == 0.0 and gp <= 1.0:
        return 0.0, common_transmission * gp, common_transmission, False
    # Exact lossless limit, including very high gain where subtraction loses
    # small absolute errors.  The residual is reported against original gains.
    if abs(gp - gc - 1.0) <= 8 * np.finfo(float).eps * max(gp, gc, 1.0):
        return np.arcsinh(np.sqrt(gc)), common_transmission, common_transmission, False
    probe_loss = gp - gc <= 1.0
    u, v = np.sqrt(gp), np.sqrt(max(gc, MIN_GAIN))

    def lossless_diagonal(q):
        ap, _, ac = _log_transfer(u, v, q)
        return ac if probe_loss else ap

    # q=-log(det M).  One log(M) diagonal is zero; the other is -q.
    # Factoring the common transmission is exact because scalar damping
    # commutes with both mixing and asymmetric loss.
    low, high = lossless_diagonal(0.0), lossless_diagonal(MAX_LOG_DETERMINANT)
    if low >= 0.0 and high <= 0.0:
        q = brentq(lossless_diagonal, 0.0, MAX_LOG_DETERMINANT, xtol=1e-11)
        _, squeeze, _ = _log_transfer(u, v, q)
        asymmetric = np.exp(-2.0 * q)
        return (squeeze,
                common_transmission * (asymmetric if probe_loss else 1.0),
                common_transmission * (1.0 if probe_loss else asymmetric), False)
    squeeze = _bounded_loss_approximation(gp, gc, probe_loss=probe_loss)
    asymmetric = np.exp(-2.0 * MAX_LOG_DETERMINANT)
    return (squeeze,
            common_transmission * (asymmetric if probe_loss else 1.0),
            common_transmission * (1.0 if probe_loss else asymmetric), True)


@lru_cache(maxsize=12)
def _cached_parameters(probe_values, conjugate_values, common_loss_pct):
    common = 1.0 - common_loss_pct / 100.0
    rows = [_infer_point(gp, gc, common)
            for gp, gc in zip(probe_values, conjugate_values)]
    return np.asarray(rows, dtype=float).T


def effective_parameters(probe_gain, conjugate_gain, *, common_loss_pct=0.0):
    """Infer an effective uniform medium with a declared common-loss closure.

    A fixed common loss supplies the otherwise unidentified third parameter.
    Remaining probe-only or conjugate-only loss is selected by the signed
    gain difference after removal of the common transmission.  Returned
    arrays are fresh; cached inference never depends on seed RIN or detection.
    """
    gp, gc = np.broadcast_arrays(np.asarray(probe_gain, dtype=float),
                                 np.asarray(conjugate_gain, dtype=float))
    if np.any(~np.isfinite(gp)) or np.any(~np.isfinite(gc)):
        raise ValueError("effective Jasperse gains must be finite")
    if np.any(gp < 0.0) or np.any(gc < 0.0):
        raise ValueError("effective Jasperse gains must be nonnegative")
    loss = float(common_loss_pct)
    if not np.isfinite(loss) or not 0.0 <= loss < 100.0:
        raise ValueError("common in-cell loss must be finite and in [0,100)%")
    inferred = _cached_parameters(tuple(gp.ravel()), tuple(gc.ravel()), loss)
    squeeze, ta, tb, approximate = [row.reshape(gp.shape).copy() for row in inferred]
    return {
        "effective_squeeze": squeeze,
        "probe_transmission": ta,
        "conjugate_transmission": tb,
        "inference_approximate": approximate.astype(bool),
        "assumed_common_loss_pct": loss,
        "closure": "fixed common loss plus minimum one-arm asymmetric distributed loss",
    }


def _external_efficiency(raw, params):
    if "after_cell_eff_pct" in params:
        after = float(params["after_cell_eff_pct"]) / 100.0
        qe = float(params.get("qe_pct", 92.0)) / 100.0
        if not 0.0 <= after <= 1.0 or not 0.0 <= qe <= 1.0:
            raise ValueError("after-cell transmission and QE must be in [0,100]%")
        return after * qe
    loss = float(params.get("loss_pct", 5.5))
    qe_pct = float(params.get("qe_pct", 92.0))
    legacy = (1.0 - loss / 100.0) * qe_pct / 100.0
    direct = params.get("detection_eff_pct")
    legacy_changed = (not np.isclose(loss, 5.5, rtol=0.0, atol=1e-12)
                      or not np.isclose(qe_pct, 92.0, rtol=0.0, atol=1e-12))
    if direct is not None:
        return legacy if legacy_changed else float(direct) / 100.0
    if "loss_pct" in params or "qe_pct" in params:
        return legacy
    return float(raw.get("eta", 1.0))


def experimental_jasperse_readout(raw, params, *, gains=None):
    """Evaluate effective two-arm distributed-loss NF over the entire scan.

    ``gains`` may supply the empirical gain-calibrated spectrum.  In that case
    its pump-off-reference units are adopted as effective source-gain units;
    this declared phenomenological identification is not an atomic covariance
    measurement.  Both reservoir contributions are always included.
    """
    source = raw if gains is None else gains
    gp, gc = np.broadcast_arrays(np.asarray(source["G_s"], dtype=float),
                                 np.asarray(source["G_c"], dtype=float))
    effective = effective_parameters(
        gp, gc, common_loss_pct=params.get("jasperse_loss_pct", 0.0))
    eta = _external_efficiency(raw, params)
    seed = seed_noise_components(
        np.broadcast_to(raw.get("input_seed_power_W", float(params.get("probe_uw", 8.0)) * 1e-6),
                        gp.shape),
        np.broadcast_to(raw.get("omega_probe_rad_s", constants.OMEGA_D1), gp.shape),
        classical_rin_per_Hz=float(params.get("seed_classical_rin_scaled", 0.0)) * 1e-14)
    result = jasperse_distributed_noise(
        effective["effective_squeeze"], effective["probe_transmission"],
        effective["conjugate_transmission"], seed_fano=seed["seed_fano"],
        probe_efficiency=eta, conjugate_efficiency=eta)
    # Inference uses a positive floor only to parameterize exactly dark
    # samples; such a floor must not create a detected bright NF.
    bright = (gp + gc > 0.0) & (eta > 0.0)
    nf = np.where(bright, result["NF"], np.nan)
    db = np.where(bright, result["squeezing_degree_dB"], np.nan)
    fit_gs = np.abs(result["transfer"][..., 0, 0]) ** 2
    fit_gc = np.abs(result["transfer"][..., 1, 0]) ** 2
    residual = np.maximum(np.abs(fit_gs - gp) / np.maximum(gp, 1e-12),
                          np.abs(fit_gc - gc) / np.maximum(gc, 1e-12))
    status = np.where(bright, "EFFECTIVE_JASPERSE_ESTIMATE", "NO_DETECTED_BRIGHT_OUTPUT")
    status = np.where(bright & effective["inference_approximate"],
                      "EFFECTIVE_JASPERSE_GAIN_APPROXIMATION", status)
    return dict(result, **effective, NF=nf, NF_dB=db, seed_noise=seed,
                eta_after_cell_QE=eta, model_id=MODEL_ID, status=status,
                fit_relative_residual=residual, G_s=fit_gs, G_c=fit_gc,
                input_G_s=gp.copy(), input_G_c=gc.copy(),
                empirical_gain_calibration_applied=bool(source.get("applied", False)),
                gain_basis=source.get("gain_quantity", "atomic-medium power gain"),
                scope=("Experimental effective uniform real phase-matched Jasperse model; "
                       "fixed common-loss closure and both distributed vacuum reservoirs; "
                       "not a microscopic atomic or finite-RF covariance prediction"))
