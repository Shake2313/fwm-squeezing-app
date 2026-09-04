import numpy as np

from gabes.lineshape import fwhm_halfmax, fwhm_interp


def test_fwhm_uses_the_connected_region_around_the_tallest_peak():
    x = np.arange(11.0)
    y = np.array([0.0, 0.0, 0.5, 1.0, 0.5, 0.0,
                  0.0, 0.6, 0.8, 0.6, 0.0])

    assert fwhm_halfmax(x, y) == 2.0
    assert fwhm_interp(x, y) == 2.0
