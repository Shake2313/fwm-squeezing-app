"""Seed noise units and the explicitly specified Jasperse gain/loss reference.

References: Jasperse (2010), thesis sections 5.2 and A.3; Jasperse et al.,
Opt. Express 19, 3765 (2011), equations (3)--(17). The lowercase UI convention
``a = classical RIN, b = shot RIN`` is declared here; it is not the thesis's
detector-dependent uppercase polynomial calibration ``A P**2 + B P + C``.
"""

import numpy as np

from .constants import HBAR


def seed_noise_components(power_W, omega_rad_s, *, classical_rin_per_Hz=0.0):
    """Return one-sided seed intensity-noise PSDs at the specified input plane.

    ``a`` is classical fractional intensity-noise PSD [1/Hz]. The calculated
    coherent Poisson baseline ``b`` is ``2*hbar*omega/P`` [1/Hz], not an
    adjustable noise floor. Thus the amplitude Fano factor is ``1 + a/b``.
    Absolute optical-power PSDs are also returned [W**2/Hz]. No measurement
    bandwidth, detector responsivity, or electronics calibration is assumed.
    An integrated noise variance requires the effective noise bandwidth.
    All numeric inputs broadcast; P and omega must be positive and a >= 0.
    """
    power, omega, classical = np.broadcast_arrays(
        np.asarray(power_W, dtype=float),
        np.asarray(omega_rad_s, dtype=float),
        np.asarray(classical_rin_per_Hz, dtype=float))
    if np.any(~np.isfinite(power)) or np.any(power <= 0.0):
        raise ValueError("seed power_W must be finite and positive")
    if np.any(~np.isfinite(omega)) or np.any(omega <= 0.0):
        raise ValueError("seed omega_rad_s must be finite and positive")
    if np.any(~np.isfinite(classical)) or np.any(classical < 0.0):
        raise ValueError("classical_rin_per_Hz must be finite and non-negative")
    shot_rin = 2.0 * HBAR * omega / power
    return {
        "seed_power_W": power,
        "seed_omega_rad_s": omega,
        "classical_rin_per_Hz": classical,
        "shot_rin_per_Hz": shot_rin,
        "seed_fano": 1.0 + classical / shot_rin,
        "classical_power_psd_W2_per_Hz": classical * power ** 2,
        "shot_power_psd_W2_per_Hz": 2.0 * HBAR * omega * power,
        "convention": (
            "one-sided optical PSD at the seed input plane; a=classical RIN, "
            "b=2*hbar*omega/P shot RIN [1/Hz]; no integrated bandwidth assumed"),
    }


def jasperse_distributed_noise(
        squeeze_parameter, probe_transmission=1.0, conjugate_transmission=1.0,
        *, seed_fano=1.0, probe_efficiency=1.0, conjugate_efficiency=1.0):
    """Evaluate Jasperse's declared uniform, real, phase-matched reference.

    The constant generator has diagonal ``log(Ta)/2, log(Tb)/2`` and real
    off-diagonals ``S``. Both independent vacuum reservoirs are integrated
    exactly; probe and conjugate losses may differ. ``Ta,Tb`` are intensity
    transmissions with mixing removed. Efficiencies describe *external*
    downstream loss only. The seed is bright, the conjugate input vacuum.

    This is not a bath inferred from an arbitrary susceptibility/transfer M,
    and not a microscopic atomic Langevin or finite-RF prediction. It does
    not establish that a pump-dressed production drift has this generator.
    Negative formal variances are never clipped. Numeric inputs broadcast.
    """
    squeeze, ta, tb, fano, eta_p, eta_c = np.broadcast_arrays(*[
        np.asarray(x, dtype=float) for x in (
            squeeze_parameter, probe_transmission, conjugate_transmission,
            seed_fano, probe_efficiency, conjugate_efficiency)])
    if np.any(~np.isfinite(squeeze)) or np.any(squeeze < 0.0):
        raise ValueError("squeeze_parameter must be finite and non-negative")
    for name, value in (("probe_transmission", ta),
                        ("conjugate_transmission", tb)):
        if (np.any(~np.isfinite(value)) or np.any(value <= 0.0)
                or np.any(value > 1.0)):
            raise ValueError(f"{name} must be finite and in (0, 1]")
    if np.any(~np.isfinite(fano)) or np.any(fano <= 0.0):
        raise ValueError("seed_fano must be finite and positive")
    for name, value in (("probe_efficiency", eta_p),
                        ("conjugate_efficiency", eta_c)):
        if (np.any(~np.isfinite(value)) or np.any(value < 0.0)
                or np.any(value > 1.0)):
            raise ValueError(f"{name} must be finite and in [0, 1]")
    generator = np.zeros(squeeze.shape + (2, 2))
    generator[..., 0, 0] = 0.5 * np.log(ta)
    generator[..., 1, 1] = 0.5 * np.log(tb)
    generator[..., 0, 1] = generator[..., 1, 0] = squeeze
    eigenvalues, vectors = np.linalg.eigh(generator)
    transpose = np.swapaxes(vectors, -1, -2)
    transfer = (vectors * np.exp(eigenvalues)[..., None, :]) @ transpose
    injection = np.zeros_like(generator)
    injection[..., 0, 0] = -np.log(ta)
    injection[..., 1, 1] = -np.log(tb)
    summed = eigenvalues[..., :, None] + eigenvalues[..., None, :]
    # Integral exp(t*(lambda_i+lambda_j)) dt, including its removable zero.
    integral = np.ones_like(summed)
    np.divide(np.expm1(summed), summed, out=integral, where=summed != 0.0)
    reservoir = vectors @ ((transpose @ injection @ vectors) * integral) @ transpose
    input_covariance = np.zeros_like(generator)
    input_covariance[..., 0, 0] = fano
    input_covariance[..., 1, 1] = 1.0
    covariance = (transfer @ input_covariance
                  @ np.swapaxes(transfer, -1, -2) + reservoir)
    attenuation = np.stack((np.sqrt(eta_p), np.sqrt(eta_c)), axis=-1)
    mean = transfer[..., :, 0] * attenuation
    detected_covariance = covariance * attenuation[..., :, None] * attenuation[..., None, :]
    detected_covariance[..., 0, 0] += 1.0 - eta_p
    detected_covariance[..., 1, 1] += 1.0 - eta_c
    difference = mean * np.array([1.0, -1.0])
    numerator = np.einsum("...i,...ij,...j->...", difference,
                          detected_covariance, difference)
    total = np.sum(mean ** 2, axis=-1)
    with np.errstate(invalid="ignore", divide="ignore"):
        nf = np.where(total > 0.0, numerator / total, np.nan)
        db = np.where(nf > 0.0, 10.0 * np.log10(nf), np.nan)
    return {
        "NF": nf,
        "squeezing_degree_dB": db,
        "transfer": transfer,
        "quadrature_covariance": covariance,
        "reservoir_quadrature_covariance": reservoir,
        "detected_probe_gain": mean[..., 0] ** 2,
        "detected_conjugate_gain": mean[..., 1] ** 2,
        "scope": (
            "Jasperse uniform real phase-matched distributed gain/loss reference; "
            "bright seed, vacuum conjugate/reservoirs, static linearization; "
            "external efficiencies exclude in-cell transmissions"),
    }
