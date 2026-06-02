"""
Phase-2 checks for the FWM twin-beam coincidence observables (ideal parametric).

Analytic: gains from an ideal two-mode squeezer reproduce the textbook photon-
pair statistics. FWM: the spectrum is nonclassical (Cauchy-Schwarz R > 1, cross-
correlation g²_sc > 2) everywhere there is net parametric gain.

    python tests/test_coincidence.py   # or: pytest tests/test_coincidence.py
"""
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from gabes import observables  # noqa: E402
from gabes.schemes import fwm  # noqa: E402


def test_ideal_two_mode_squeezer():
    r = 0.8
    G_s = np.array([np.cosh(r) ** 2])
    G_c = np.array([np.sinh(r) ** 2])
    c = observables.coincidence_stats(G_s, G_c)
    n = np.sinh(r) ** 2                      # mean photons per mode
    assert np.isclose(c["n_s"][0], n)        # G_s − 1 = sinh²r
    assert np.isclose(c["g2_ss"][0], 2.0) and np.isclose(c["g2_cc"][0], 2.0)
    assert np.isclose(c["g2_sc"][0], 2.0 + 1.0 / n)        # textbook TMSV
    assert np.isclose(c["cauchy_schwarz"][0], (2.0 + 1.0 / n) ** 2 / 4.0)
    assert c["cauchy_schwarz"][0] > 1.0                    # nonclassical


def test_fwm_spectrum_nonclassical_in_gain_region():
    center = fwm.branch_center_GHz(0.9, -1)
    spec = fwm.compute_spectrum(
        0.9, T=394.15, P_pump=0.6, P_probe=8e-6, line_strength=0.05, loss_frac=0.055,
        coarse_points=121, fine_points=0, scan_min=center - 0.55, scan_max=center + 0.55,
        velocity_step=5.0, velocity_cutoff=3.0, branches=fwm.BRANCHES)
    c = observables.coincidence_stats(spec["G_s"], spec["G_c"])
    gain = spec["G_s"] > 1.0
    assert gain.any()                                       # there is a gain region
    assert np.all(c["g2_sc"][gain] > 2.0)                   # nonclassical cross-corr
    assert np.all(c["cauchy_schwarz"][gain] > 1.0)          # Cauchy-Schwarz violation
    assert np.all(np.isnan(c["g2_sc"][~gain]))              # undefined without gain
    imax = int(np.argmax(spec["G_s"]))                      # peak gain
    assert np.isfinite(c["g2_sc"][imax]) and c["cauchy_schwarz"][imax] > 1.0


if __name__ == "__main__":
    test_ideal_two_mode_squeezer()
    test_fwm_spectrum_nonclassical_in_gain_region()
    print("Phase-2 coincidence OK (ideal TMSV + FWM gain-region nonclassicality).")
