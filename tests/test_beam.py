import math
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from gabes import beam, constants  # noqa: E402


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
