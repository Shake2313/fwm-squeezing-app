"""Photon shot noise in transmission and ideal linear-polarization rotation.

One-sided photocurrent PSD S_i = 2 e I, integrated over electrical ENBW.
Transmission uses a fixed, noiseless incident-power normalization. Rotation
uses an ideal balanced polarimeter: d(I_+ - I_-)/dtheta = 2 I_total near balance.
Thus sigma_T = sqrt(2 h nu B T / (eta P_in)) and
sigma_theta = sqrt(h nu B / (2 eta P_in T)) [rad].

The atomic engine uses a uniform intensity. Its effective collected top-hat
area supplies P_in = I pi (d/2)^2; this is not a Gaussian beam integration.
Technical/electronic noise and atomic spin noise are outside this reference.
"""
import numpy as np

from . import constants


NOISE_SEED = 20261023
DEFAULT_BEAM_DIAMETER_MM = 2.0
DEFAULT_DETECTION_EFFICIENCY = 0.8
DEFAULT_NOISE_BANDWIDTH_HZ = 1000.0


def magnetic_sensitivity(transmission, slope_rad_per_t, *, frequency_hz,
                         intensity_mw_cm2,
                         beam_diameter_mm=DEFAULT_BEAM_DIAMETER_MM,
                         efficiency=DEFAULT_DETECTION_EFFICIENCY,
                         bandwidth_hz=DEFAULT_NOISE_BANDWIDTH_HZ,
                         signal_rms_t=0.0, linear_input=True, slope_valid=True,
                         signal_peak_limit_t=None):
    """Refer ideal polarimeter photon noise to a local magnetic-field slope.

    ``transmission`` and ``slope_rad_per_t`` must describe the same operating
    point. Transmission, slope, and the caller's ``slope_valid`` convergence
    mask are broadcast together. The slope is in rad/T, not rad/uT. The returned
    field ASD is one-sided T/sqrt(Hz); multiplying by sqrt(electrical ENBW)
    gives RMS field noise. SNR uses the supplied *RMS* small-signal field.

    This passive-cell, linear-input approximation requires 0 <= T <= 1 and
    positive detected flux. Unlike ``shot_noise_reference``, it does not draw
    Gaussian samples or require a minimum photon count per averaging interval.
    Consequently ENBW never determines ASD availability. Invalid samples use
    JSON null with a status, including zero slopes (unbounded field noise).
    Rotation ASD remains available when only the field slope is invalid.
    An optional symmetric tested field interval limits SNR for a sinusoidal
    signal (peak = sqrt(2) RMS), without changing the infinitesimal field ASD.
    """
    frequency_hz = float(frequency_hz)
    intensity_mw_cm2 = float(intensity_mw_cm2)
    beam_diameter_mm = float(beam_diameter_mm)
    efficiency = float(efficiency)
    bandwidth_hz = float(bandwidth_hz)
    signal_rms_t = float(signal_rms_t)
    if signal_peak_limit_t is not None:
        signal_peak_limit_t = float(signal_peak_limit_t)
        if not np.isfinite(signal_peak_limit_t) or signal_peak_limit_t < 0:
            raise ValueError("Tested signal peak limit must be finite and nonnegative.")
    settings = np.asarray([frequency_hz, intensity_mw_cm2, beam_diameter_mm,
                           efficiency, bandwidth_hz, signal_rms_t])
    if (not np.all(np.isfinite(settings)) or frequency_hz <= 0
            or beam_diameter_mm <= 0 or bandwidth_hz <= 0
            or intensity_mw_cm2 < 0 or not 0 <= efficiency <= 1
            or signal_rms_t < 0):
        raise ValueError("Magnetic sensitivity requires positive finite frequency, "
                         "diameter and ENBW, nonnegative intensity and RMS signal, "
                         "and efficiency in [0, 1].")
    transmission, slope, trusted_slope = np.broadcast_arrays(
        np.asarray(transmission, dtype=float),
        np.asarray(slope_rad_per_t, dtype=float),
        np.asarray(slope_valid, dtype=bool))
    photon_energy = 2 * np.pi * constants.HBAR * frequency_hz
    with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
        area = np.pi * np.square(np.float64(beam_diameter_mm) * 1e-3 / 2)
        incident_power = intensity_mw_cm2 * 10.0 * area
        incident_flux = efficiency * incident_power / photon_energy
    if not np.all(np.isfinite([incident_power, incident_flux])):
        raise ValueError("Magnetic sensitivity requires finite derived power and photon flux.")

    physical_t = np.isfinite(transmission) & (transmission >= 0) & (transmission <= 1)
    flux = np.where(physical_t, incident_flux * np.where(physical_t, transmission, 0), np.nan)
    rotation_ok = physical_t & (flux > 0) & bool(linear_input)
    rotation_status = np.full(transmission.shape, "available", dtype=object)
    rotation_status[~physical_t] = "transmission must be finite and within [0, 1]"
    rotation_status[physical_t & ~(flux > 0)] = "no detected photon flux"
    if not linear_input:
        rotation_status[...] = "unsupported polarization"
    rotation_asd = np.full(transmission.shape, np.nan)
    np.divide(np.sqrt(0.5), np.sqrt(flux), out=rotation_asd, where=rotation_ok)

    status = rotation_status.copy()
    status[rotation_ok & ~trusted_slope] = "slope not validated"
    status[rotation_ok & trusted_slope & ~np.isfinite(slope)] = "nonfinite magnetic slope"
    status[rotation_ok & trusted_slope & np.isfinite(slope) & (slope == 0)] = (
        "zero magnetic slope; unbounded field noise")
    slope_ok = trusted_slope & np.isfinite(slope) & (slope != 0)
    field_asd = np.full(transmission.shape, np.nan)
    with np.errstate(over="ignore", under="ignore", invalid="ignore", divide="ignore"):
        np.divide(rotation_asd, np.abs(slope), out=field_asd, where=rotation_ok & slope_ok)
    valid = rotation_ok & slope_ok & np.isfinite(field_asd) & (field_asd > 0)
    status[rotation_ok & slope_ok & ~valid] = "field ASD outside floating-point range"
    field_asd = np.where(valid, field_asd, np.nan)
    with np.errstate(over="ignore", under="ignore", invalid="ignore", divide="ignore"):
        field_rms = field_asd * np.sqrt(bandwidth_hz)
        rms_valid = valid & np.isfinite(field_rms) & (field_rms > 0)
        snr = np.full(transmission.shape, np.nan)
        np.divide(signal_rms_t, field_rms, out=snr, where=rms_valid)
    snr_valid = rms_valid & np.isfinite(snr)
    rms_status = np.where(valid & ~rms_valid, "field RMS outside floating-point range", status)
    snr_status = np.where(rms_valid & ~snr_valid, "SNR outside floating-point range", rms_status)
    linearity_status = "signal amplitude linearity not assessed"
    if signal_peak_limit_t is not None:
        linearity_status = "sinusoidal peak within tested local field interval"
        if signal_rms_t > signal_peak_limit_t / np.sqrt(2.0):
            linearity_status = "signal peak exceeds tested local field interval"
            snr_status = np.where(rms_valid, linearity_status, snr_status)
            snr_valid = np.zeros_like(snr_valid, dtype=bool)

    def finite_values(values):
        array = np.asarray(values, dtype=float)
        return np.where(np.isfinite(array), array, None).tolist()

    return dict(
        source="coherent-light photon shot noise; ideal linear-input balanced polarimeter",
        scope="local small-signal DC sensitivity; technical, electronics and atomic spin noise omitted",
        asd_convention="one-sided amplitude spectral density",
        bandwidth_convention="one-sided electrical equivalent noise bandwidth",
        signal_convention="RMS magnetic-field amplitude; SNR = signal RMS / field-noise RMS",
        formula="S_theta^1/2 = 1/sqrt(2 photon_flux); S_B^1/2 = S_theta^1/2 / abs(dtheta/dB)",
        beam_model="uniform effective collected area; full beam collection",
        bandwidth_hz=bandwidth_hz, signal_rms_t=signal_rms_t,
        signal_peak_limit_t=signal_peak_limit_t, signal_linearity_status=linearity_status,
        frequency_hz=frequency_hz, intensity_mw_cm2=intensity_mw_cm2,
        beam_diameter_mm=beam_diameter_mm, efficiency=efficiency,
        incident_power_w=float(incident_power), photon_energy_j=float(photon_energy),
        transmission=finite_values(transmission), slope_rad_per_t=finite_values(slope),
        detected_photon_flux_per_s=finite_values(flux),
        rotation_asd_rad_per_sqrt_hz=finite_values(rotation_asd),
        field_asd_t_per_sqrt_hz=finite_values(field_asd),
        field_rms_t=finite_values(np.where(rms_valid, field_rms, np.nan)),
        signal_to_noise_ratio=finite_values(np.where(snr_valid, snr, np.nan)),
        rotation_valid=rotation_ok.tolist(), rotation_status=rotation_status.tolist(),
        valid=valid.tolist(), status=status.tolist(),
        rms_valid=rms_valid.tolist(), rms_status=rms_status.tolist(),
        snr_valid=snr_valid.tolist(), snr_status=snr_status.tolist(),
    )


def shot_noise_reference(transmission, *, reference_transmission, frequency_hz,
                         intensity_mw_cm2, beam_diameter_mm=DEFAULT_BEAM_DIAMETER_MM,
                         efficiency=DEFAULT_DETECTION_EFFICIENCY,
                         bandwidth_hz=DEFAULT_NOISE_BANDWIDTH_HZ,
                         linear_input=True):
    """Return a reproducible Gaussian reference plus its physical RMS budget.

    The transmission trace fluctuates about the wing background, using that
    background's photon flux. The rotation trace fluctuates about zero and uses
    the transmitted flux at each field. Gaussian samples require at least ten
    detected photons per equivalent averaging interval 1/(2 ENBW); unsupported
    samples are gaps, not zero-noise claims. Lists keep headless/export payloads
    directly comparable and JSON serializable.
    """
    settings = np.asarray([frequency_hz, intensity_mw_cm2, beam_diameter_mm,
                           efficiency, bandwidth_hz], dtype=float)
    if (not np.all(np.isfinite(settings)) or frequency_hz <= 0 or beam_diameter_mm <= 0
            or bandwidth_hz <= 0 or intensity_mw_cm2 < 0 or not 0 <= efficiency <= 1):
        raise ValueError("Shot noise requires positive finite frequency, diameter and ENBW, "
                         "nonnegative intensity, and efficiency in [0, 1].")
    transmission = np.asarray(transmission, dtype=float)
    if (not np.all(np.isfinite(transmission)) or np.any(transmission < 0)
            or not np.isfinite(reference_transmission) or reference_transmission < 0):
        raise ValueError("Shot noise requires finite nonnegative transmission.")

    photon_energy = 2 * np.pi * constants.HBAR * frequency_hz
    area = np.pi * (beam_diameter_mm * 1e-3 / 2) ** 2
    incident_power = intensity_mw_cm2 * 10.0 * area
    incident_flux = efficiency * incident_power / photon_energy
    detected_flux = incident_flux * transmission
    reference_flux = incident_flux * reference_transmission
    # With no detected incident photons the transmission normalization is
    # undefined, even though the photocurrent itself has zero shot noise.
    sigma_t = np.full(transmission.shape, np.nan)
    sigma_t_ref = float("nan")
    if incident_flux > 0:
        sigma_t = np.sqrt(2 * bandwidth_hz * transmission / incident_flux)
        sigma_t_ref = float(np.sqrt(2 * bandwidth_hz * reference_transmission / incident_flux))
    transmission_ok = reference_flux / (2 * bandwidth_hz) >= 10
    rotation_ok = (detected_flux / (2 * bandwidth_hz) >= 10) & linear_input
    sigma_theta = np.full(transmission.shape, np.nan)
    np.divide(bandwidth_hz, 2 * detected_flux, out=sigma_theta, where=rotation_ok)
    sigma_theta = np.sqrt(sigma_theta)

    rng = np.random.default_rng(NOISE_SEED)
    t_trace = (reference_transmission + rng.standard_normal(transmission.shape) * sigma_t_ref
               if transmission_ok else np.full(transmission.shape, np.nan))
    theta_trace = rng.standard_normal(transmission.shape) * sigma_theta
    theta_trace = np.where(rotation_ok, theta_trace, np.nan)
    sigma_theta = np.where(rotation_ok, sigma_theta, np.nan)

    def finite_list(values):
        return [float(v) if np.isfinite(v) else None for v in values]

    return dict(
        source="coherent-light photon shot noise; simulated Gaussian reference",
        bandwidth_convention="one-sided electrical equivalent noise bandwidth",
        bandwidth_hz=float(bandwidth_hz), efficiency=float(efficiency),
        beam_diameter_mm=float(beam_diameter_mm),
        beam_model="uniform effective collected area; full beam collection",
        incident_power_w=float(incident_power), photon_energy_j=float(photon_energy),
        detected_photon_flux_per_s=finite_list(detected_flux),
        reference_transmission=float(reference_transmission),
        transmission_normalization="fixed noiseless incident-power reference",
        transmission_rms=finite_list(sigma_t),
        transmission_reference_rms=sigma_t_ref if np.isfinite(sigma_t_ref) else None,
        transmission_reference_trace=finite_list(t_trace),
        transmission_status=("available" if transmission_ok else
                             "insufficient detected photons for Gaussian reference"),
        rotation_rms_rad=finite_list(sigma_theta),
        rotation_reference_trace_rad=finite_list(theta_trace),
        rotation_status=("unsupported polarization" if not linear_input else
                         "insufficient detected photons for Gaussian reference" if not np.any(rotation_ok)
                         else "partial: insufficient detected photons at some fields" if not np.all(rotation_ok)
                         else "linear-input balanced-polarimeter approximation"),
        seed=NOISE_SEED,
        scope="photon shot noise only; technical, electronics and atomic spin noise omitted",
    )
