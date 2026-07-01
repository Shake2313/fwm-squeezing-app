"""Shared beam-scaling and simple wave-vector geometry helpers."""
import math

from . import constants

MHZ = 2.0 * math.pi * 1e6


def anchored_rabi_mhz(anchor_mhz, power, ref_power, diameter=None,
                      ref_diameter=None, min_diameter=1e-9):
    """Scale an anchored Rabi frequency as sqrt(power) / beam diameter."""
    ref_power = float(ref_power)
    if ref_power <= 0.0:
        raise ValueError("ref_power must be positive")
    scale = math.sqrt(max(float(power), 0.0) / ref_power)
    if diameter is not None or ref_diameter is not None:
        if diameter is None or ref_diameter is None:
            raise ValueError("diameter and ref_diameter must be supplied together")
        scale *= float(ref_diameter) / max(float(diameter), min_diameter)
    return float(anchor_mhz) * scale


def transit_broadening_mhz(temp_k, diameter_mm, *, mass=constants.MASS_85RB,
                           factor=1.0, min_diameter_m=1e-9):
    """Transit broadening / 2pi in MHz from most-probable speed over diameter."""
    v_mp = math.sqrt(2.0 * constants.KB * float(temp_k) / float(mass))
    diameter_m = max(float(diameter_mm) * 1e-3, min_diameter_m)
    return float(factor) * v_mp / diameter_m / MHZ


def wavevector_from_wavelength_nm(wavelength_nm):
    """Optical wave-vector magnitude for a vacuum wavelength in nm."""
    return 2.0 * math.pi / (float(wavelength_nm) * 1e-9)


def collinear_residual_k_ratio(k_probe, k_coupling):
    """Signed collinear residual (k_probe - k_coupling) / k_probe."""
    kp = float(k_probe)
    return (kp - float(k_coupling)) / max(abs(kp), 1e-30)


def residual_wavevector_magnitude(k_probe, k_coupling, angle_rad=0.0):
    """Magnitude of k_probe - k_coupling for an included beam angle."""
    kp = float(k_probe)
    kc = float(k_coupling)
    cos_a = math.cos(float(angle_rad))
    return math.sqrt(max(kp * kp + kc * kc - 2.0 * kp * kc * cos_a, 0.0))


def angled_residual_k_ratio(k_probe, k_coupling, *, angle_mrad=0.0):
    """Unsigned residual |Delta k| / |k_probe| for a small beam-angle mismatch."""
    kp = float(k_probe)
    angle_rad = float(angle_mrad) * 1e-3
    dk = residual_wavevector_magnitude(kp, k_coupling, angle_rad)
    return dk / max(abs(kp), 1e-30)
