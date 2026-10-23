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
