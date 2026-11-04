import numpy as np
import pytest

from gabes.lineshape import (fwhm_halfmax, fwhm_interp, halfwidth_from_center,
                             level_crossing, subdoppler_feature)
from gabes.schemes import magneto
from gabes.schemes.rydberg import RydbergEITScheme


def test_fwhm_uses_the_connected_region_around_the_tallest_peak():
    x = np.arange(11.0)
    y = np.array([0.0, 0.0, 0.5, 1.0, 0.5, 0.0,
                  0.0, 0.6, 0.8, 0.6, 0.0])

    assert fwhm_halfmax(x, y) == 2.0
    assert fwhm_interp(x, y) == 2.0


def test_level_crossing_interpolates_from_the_last_inside_sample():
    x = np.arange(7.0)
    y = np.array([0.0, 0.25, 0.75, 1.0, 0.75, 0.25, 0.0])
    assert level_crossing(x, y, 3, +1, 0.5) == (4.5, 5)
    assert level_crossing(x, y, 3, -1, 0.5) == (1.5, 1)
    # A dip is followed with sign=-1 at the same level.
    assert level_crossing(x, 1.0 - y, 3, +1, 0.5, sign=-1.0) == (4.5, 5)


def test_level_crossing_reports_a_walk_off_the_scan_and_stops_at_non_finite():
    x = np.arange(5.0)
    crossing, index = level_crossing(x, np.array([1.0, 1.0, 1.0, 0.9, 0.8]), 2, +1, 0.5)
    assert np.isnan(crossing) and index is None
    crossing, index = level_crossing(x, np.array([1.0, 1.0, 1.0, np.nan, 0.0]), 2, +1, 0.5)
    assert np.isnan(crossing) and index == 3


@pytest.mark.parametrize("hwhm", [0.4, 1.3, 3.0])
def test_scheme_and_generic_widths_match_the_lorentzian_at_their_own_level(hwhm):
    # Each caller picks its own half-height level (peak/2, floor or wing
    # background midpoint); all of them walk the same interpolated edge.
    x = np.linspace(-10.0, 10.0, 401)
    y = 1.0 / (1.0 + (x / hwhm) ** 2)

    def full_width(level):
        return 2.0 * hwhm * np.sqrt(1.0 / level - 1.0)

    ic = int(np.argmax(y))
    assert fwhm_interp(x, y) == pytest.approx(full_width(0.5), rel=2e-3)
    wings = 0.5 * (y[0] + y[-1])
    assert halfwidth_from_center(x, y) == pytest.approx(
        0.5 * full_width(wings + 0.5 * (y[ic] - wings)), rel=2e-3)
    assert magneto._central_feature(x, y)["width"] == pytest.approx(
        full_width(wings + 0.5 * (y[ic] - wings)), rel=2e-3)
    floor = y[np.abs(x) <= 5.0].min()
    width, _ = RydbergEITScheme._eit_features(x, y)
    assert width == pytest.approx(full_width(floor + 0.5 * (y[ic] - floor)), rel=2e-3)


def test_subdoppler_edges_bracket_a_narrow_dip_on_a_doppler_background():
    x = np.linspace(-3.0, 3.0, 2401)
    background = np.exp(-(x / 1.5) ** 2)
    signal = background - 0.05 / (1.0 + (x / 0.02) ** 2)
    feature = subdoppler_feature(x, signal)
    assert feature.resolved
    assert feature.left_half_height < 0.0 < feature.right_half_height
    # The running-median residual is somewhat narrower than the 0.04 dip.
    assert 0.02 < feature.fwhm <= 0.04
