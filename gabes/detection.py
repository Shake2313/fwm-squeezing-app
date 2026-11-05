"""Photodetection: photon energy and flux, responsivity, shot noise, detectors.

One implementation, for every scheme, of the conversions between optical power,
photon flux and photocurrent and of coherent-light shot noise:

    E = h nu,   Phi = eta P / (h nu),   R = eta e lambda / (h c),
    S_i = 2 e I  [A^2/Hz],   S_P = 2 hbar omega P  [W^2/Hz],

all one-sided.  It also holds the balanced-detector budget: two photodiode arms,
common-mode RIN with a correlation coefficient, and input-referred electronics
noise of the difference output.  Noise amplitudes are one-sided ASDs per
sqrt(Hz); an RMS needs the electrical equivalent noise bandwidth (asd_to_rms).
"""

from dataclasses import dataclass

import numpy as np

from . import constants


def _require_finite_nonnegative(name, value):
    value = float(value)
    if not np.isfinite(value) or value < 0.0:
        raise ValueError(f"{name} must be finite and non-negative")
    return value


def _require_finite_positive(name, value):
    value = float(value)
    if not np.isfinite(value) or value <= 0.0:
        raise ValueError(f"{name} must be finite and positive")
    return value


def photon_energy(frequency_hz):
    """Photon energy h*nu [J] at optical frequency nu [Hz]."""
    return constants.PLANCK * frequency_hz


def photon_flux(power_w, frequency_hz, efficiency=1.0):
    """Detected photon flux eta*P/(h*nu) [1/s]."""
    return efficiency * power_w / photon_energy(frequency_hz)


def shot_noise_current_psd(current_a):
    """One-sided Schottky current-noise PSD S_i = 2 e I [A^2/Hz]."""
    return 2.0 * constants.ELEMENTARY_CHARGE * current_a


def optical_shot_power_psd(power_w, angular_frequency_rad_s):
    """One-sided optical-power shot-noise PSD 2 hbar omega P [W^2/Hz] of coherent light."""
    return 2.0 * constants.HBAR * angular_frequency_rad_s * power_w


def optical_shot_rin(power_w, angular_frequency_rad_s):
    """One-sided coherent-light shot RIN 2 hbar omega / P [1/Hz]."""
    return 2.0 * constants.HBAR * angular_frequency_rad_s / power_w


def photodiode_responsivity_a_per_w(quantum_efficiency, wavelength_m):
    """Ideal photodiode responsivity eta*e*lambda/(h*c) [A/W]."""
    eta = float(quantum_efficiency)
    if not np.isfinite(eta) or not 0.0 <= eta <= 1.0:
        raise ValueError("quantum_efficiency must be between 0 and 1")
    wavelength = _require_finite_positive("wavelength_m", wavelength_m)
    return eta * constants.ELEMENTARY_CHARGE * wavelength / (
        constants.PLANCK * constants.C_LIGHT)


@dataclass(frozen=True)
class PhotodiodeChannel:
    """One photodiode arm represented by measured optical/detector quantities."""

    optical_power_w: float
    responsivity_a_per_w: float
    dark_current_a: float = 0.0

    def __post_init__(self):
        _require_finite_nonnegative("optical_power_w", self.optical_power_w)
        _require_finite_nonnegative(
            "responsivity_a_per_w", self.responsivity_a_per_w)
        _require_finite_nonnegative("dark_current_a", self.dark_current_a)

    @classmethod
    def from_quantum_efficiency(cls, optical_power_w, quantum_efficiency,
                                wavelength_m, *, dark_current_a=0.0):
        return cls(
            optical_power_w=optical_power_w,
            responsivity_a_per_w=photodiode_responsivity_a_per_w(
                quantum_efficiency, wavelength_m),
            dark_current_a=dark_current_a)

    @property
    def photocurrent_a(self):
        return float(self.optical_power_w) * float(self.responsivity_a_per_w)


@dataclass(frozen=True)
class BalancedDetector:
    """Signal-minus-weighted-reference detector configuration.

    ``electronic_noise_current_asd_a_per_sqrt_hz`` is the input-referred noise of
    the *difference output*.  ``relative_intensity_noise_per_sqrt_hz`` is a
    one-sided fractional optical-power ASD common to the two arms, and
    ``rin_correlation`` is their correlation coefficient (-1 to 1).  Dark-current
    shot noise is included, while dark current does not contribute RIN.
    """

    signal: PhotodiodeChannel
    reference: PhotodiodeChannel | None = None
    reference_weight: float = 1.0
    electronic_noise_current_asd_a_per_sqrt_hz: float = 0.0
    relative_intensity_noise_per_sqrt_hz: float = 0.0
    rin_correlation: float = 1.0

    def __post_init__(self):
        weight = float(self.reference_weight)
        if not np.isfinite(weight):
            raise ValueError("reference_weight must be finite")
        _require_finite_nonnegative(
            "electronic_noise_current_asd_a_per_sqrt_hz",
            self.electronic_noise_current_asd_a_per_sqrt_hz)
        _require_finite_nonnegative(
            "relative_intensity_noise_per_sqrt_hz",
            self.relative_intensity_noise_per_sqrt_hz)
        corr = float(self.rin_correlation)
        if not np.isfinite(corr) or not -1.0 <= corr <= 1.0:
            raise ValueError("rin_correlation must lie between -1 and 1")


@dataclass(frozen=True)
class DetectorNoiseBudget:
    """One-sided current-noise ASD budget for a balanced detector."""

    signal_photocurrent_a: float
    reference_photocurrent_a: float
    shot_noise_current_asd_a_per_sqrt_hz: float
    rin_noise_current_asd_a_per_sqrt_hz: float
    electronic_noise_current_asd_a_per_sqrt_hz: float
    technical_noise_current_asd_a_per_sqrt_hz: float
    total_noise_current_asd_a_per_sqrt_hz: float


def balanced_detector_noise(detector):
    """Calculate PSN, technical, and total one-sided current-noise ASDs."""
    if not isinstance(detector, BalancedDetector):
        raise TypeError("detector must be a BalancedDetector")
    sig = detector.signal
    ref = detector.reference
    i_sig = sig.photocurrent_a
    i_ref = 0.0 if ref is None else ref.photocurrent_a
    dark_sig = float(sig.dark_current_a)
    dark_ref = 0.0 if ref is None else float(ref.dark_current_a)
    weight = float(detector.reference_weight)

    # One-sided Schottky current-noise PSD S_i = 2 e I for independent arms.
    shot_variance_density = shot_noise_current_psd(
        i_sig + dark_sig + weight * weight * (i_ref + dark_ref))
    shot_asd = np.sqrt(max(shot_variance_density, 0.0))

    rin = float(detector.relative_intensity_noise_per_sqrt_hz)
    corr = float(detector.rin_correlation)
    rin_variance_density = rin * rin * (
        i_sig * i_sig + weight * weight * i_ref * i_ref
        - 2.0 * weight * corr * i_sig * i_ref)
    rin_asd = np.sqrt(max(rin_variance_density, 0.0))
    electronic_asd = float(
        detector.electronic_noise_current_asd_a_per_sqrt_hz)
    technical_asd = np.hypot(rin_asd, electronic_asd)
    total_asd = np.hypot(shot_asd, technical_asd)
    return DetectorNoiseBudget(
        signal_photocurrent_a=i_sig,
        reference_photocurrent_a=i_ref,
        shot_noise_current_asd_a_per_sqrt_hz=shot_asd,
        rin_noise_current_asd_a_per_sqrt_hz=rin_asd,
        electronic_noise_current_asd_a_per_sqrt_hz=electronic_asd,
        technical_noise_current_asd_a_per_sqrt_hz=technical_asd,
        total_noise_current_asd_a_per_sqrt_hz=total_asd,
    )


def asd_to_rms(asd_per_sqrt_hz, equivalent_noise_bandwidth_hz):
    """Convert a white-noise ASD to RMS noise over an ENBW."""
    asd = np.asarray(asd_per_sqrt_hz, dtype=float)
    if np.any(~np.isfinite(asd)) or np.any(asd < 0.0):
        raise ValueError("asd_per_sqrt_hz must be finite and non-negative")
    enbw = _require_finite_nonnegative(
        "equivalent_noise_bandwidth_hz", equivalent_noise_bandwidth_hz)
    out = asd * np.sqrt(enbw)
    return float(out) if out.ndim == 0 else out
