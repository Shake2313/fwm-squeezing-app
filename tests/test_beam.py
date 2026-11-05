import math
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from gabes import beam, constants, probe_saturation, species  # noqa: E402


def test_anchored_rabi_scales_with_power_and_diameter():
    base = beam.anchored_rabi_mhz(
        3.0, 10.0, 10.0, diameter=0.5, ref_diameter=0.5)
    high_power = beam.anchored_rabi_mhz(
        3.0, 40.0, 10.0, diameter=0.5, ref_diameter=0.5)
    wider = beam.anchored_rabi_mhz(
        3.0, 10.0, 10.0, diameter=1.0, ref_diameter=0.5)
    assert base == 3.0
    assert high_power == 6.0
    assert wider == 1.5


def test_anchored_rabi_rejects_nonpositive_reference_power():
    with pytest.raises(ValueError):
        beam.anchored_rabi_mhz(3.0, 1.0, 0.0)


def test_transit_broadening_matches_v_over_d_formula():
    temp_k = 293.15
    diameter_mm = 0.15
    factor = 0.6
    v_mp = math.sqrt(2.0 * constants.KB * temp_k / constants.MASS_85RB)
    expected = factor * v_mp / (diameter_mm * 1e-3) / beam.MHZ
    assert abs(beam.transit_broadening_mhz(
        temp_k, diameter_mm, mass=constants.MASS_85RB, factor=factor) - expected) < 1e-12


def test_wavevector_residual_helpers():
    kp = beam.wavevector_from_wavelength_nm(780.0)
    kc = beam.wavevector_from_wavelength_nm(481.0)
    assert abs(beam.collinear_residual_k_ratio(kp, kc) - ((kp - kc) / kp)) < 1e-15

    ratio = beam.angled_residual_k_ratio(kp, kp, angle_mrad=3.0)
    assert abs(ratio - 0.003) < 1e-6


def test_peak_intensity_of_a_gaussian_beam():
    # The on-axis intensity times the effective area pi*w^2/2 is the power.
    P, w = 0.137, 5.3e-4
    assert math.isclose(beam.peak_intensity(P, w) * np.pi * w**2 / 2, P, rel_tol=1e-15)
    arr = beam.peak_intensity(np.array([0.0, 1.0, 2.0]), w)
    assert np.allclose(arr, [0.0, beam.peak_intensity(1.0, w), 2 * beam.peak_intensity(1.0, w)])


def test_two_level_rabi_is_saturation_relation():
    gamma, isat = 2 * np.pi * 6.0666e6, 25.03
    for s in (0.0, 0.3, 1.0, 40.0):
        omega = beam.two_level_rabi(s * isat, gamma, isat)
        assert math.isclose(2 * omega**2 / gamma**2, s, rel_tol=1e-14, abs_tol=1e-300)


@pytest.mark.parametrize("nu_hz, gamma_2pi_hz", [(377.107e12, 5.746e6), (384.230e12, 6.0666e6),
                                                 (335.116e12, 4.575e6)])
def test_rabi_matches_independent_dipole_reference(nu_hz, gamma_2pi_hz):
    # Closed two-level atom: d^2 = 3 pi eps0 hbar c^3 Gamma / omega^3 gives
    # Omega = d E0 / hbar, which must equal Gamma sqrt(I / 2 I_sat) with
    # I_sat = pi h c Gamma / (3 lambda^3).
    gamma = 2 * np.pi * gamma_2pi_hz
    omega0 = 2 * np.pi * nu_hz
    d = math.sqrt(3 * np.pi * constants.EPS_0 * constants.HBAR * constants.C_LIGHT**3
                  * gamma / omega0**3)
    isat = probe_saturation.closed_two_level_saturation_intensity(nu_hz, gamma)
    for intensity in (1.0, 44.84, 3.2e4):
        dipole = d * beam.field_amplitude(intensity) / constants.HBAR
        assert math.isclose(beam.two_level_rabi(intensity, gamma, isat), dipole, rel_tol=1e-12)


def test_callers_use_the_shared_conversion():
    P_w, w_m = 0.6, 5.3e-4
    reference = beam.two_level_rabi(beam.peak_intensity(P_w, w_m), constants.GAMMA, constants.I_SAT)
    assert constants.rabi_freq(P_w, w_m) == reference
    P_mw, w_mm = 137.0, 0.53                 # the mW/mm wrappers convert to SI first
    intensity = beam.peak_intensity(P_mw * 1e-3, w_mm * 1e-3)
    assert probe_saturation.gaussian_peak_intensity(P_mw, w_mm) == intensity
    assert species.pump_rabi_from_power(P_mw, w_mm, constants.GAMMA) == beam.two_level_rabi(
        intensity, constants.GAMMA, constants.I_SAT)
    assert species.pump_rabi_from_power(0.0, 1.0, constants.GAMMA) == 0.0
