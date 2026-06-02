"""
Phase-3 physics checks for saturated absorption (SAS).

  no pump  : the background envelope is the Doppler-broadened Voigt (matches OD).
  pump on  : a sub-Doppler Lamb dip appears at line centre (≪ Doppler width).
  two lines: Lamb dips at ±splitting/2 and a crossover dip at the midpoint.

    python tests/test_sas.py    # or: pytest tests/test_sas.py
"""
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from gabes import schemes, constants  # noqa: E402

G = constants.GAMMA
GMHZ = G / (2 * np.pi) / 1e6
SAS = schemes.get("sas")


def _params(**over):
    p = SAS.defaults()
    p.update(over)
    return p


def _xy(raw):
    return raw["scan"] / (2 * np.pi) / 1e6, raw["alpha"]


def _at(x, c):
    return int(np.argmin(np.abs(x - c)))


def test_no_pump_recovers_doppler_voigt():
    raw = SAS.compute(_params(pump_rabi=0.01))
    x, a = _xy(raw)
    above = x[a >= 0.5 * a.max()]
    env_fwhm = above.max() - above.min()
    dopp = raw["dopp_fwhm"] / (2 * np.pi) / 1e6
    assert 0.95 <= env_fwhm / dopp <= 1.05            # background = Voigt
    # no Lamb dip: the absorption peak sits at line centre
    assert a[_at(x, 0.0)] >= a[_at(x, 15 * GMHZ)]


def test_single_line_lamb_dip():
    raw = SAS.compute(_params(pump_rabi=2.0))
    x, a = _xy(raw)
    ic, ish = _at(x, 0.0), _at(x, 15 * GMHZ)
    assert a[ic] < 0.7 * a[ish]                       # sub-Doppler dip at centre
    # dip width ≪ Doppler width
    floor, shoulder = a[ic], a[ish]
    thr = 0.5 * (floor + shoulder)
    lo, hi = ic, ic
    while lo > 0 and a[lo] <= thr:
        lo -= 1
    while hi < a.size - 1 and a[hi] <= thr:
        hi += 1
    dip_fwhm = x[hi] - x[lo]
    assert dip_fwhm < (raw["dopp_fwhm"] / (2 * np.pi) / 1e6) / 5.0


def test_two_line_crossover():
    raw = SAS.compute(_params(transitions="two lines (crossover)",
                              splitting=60.0, pump_rabi=2.0))
    x, a = _xy(raw)
    off = 30 * GMHZ                                   # splitting/2 in MHz
    assert a[_at(x, 0.0)] < a[_at(x, 8 * GMHZ)]       # crossover dip at midpoint
    assert a[_at(x, off)] < a[_at(x, off + 8 * GMHZ)]  # Lamb dip at the line


if __name__ == "__main__":
    test_no_pump_recovers_doppler_voigt()
    test_single_line_lamb_dip()
    test_two_line_crossover()
    print("Phase-3 SAS OK (Voigt background, Lamb dip, crossover).")
