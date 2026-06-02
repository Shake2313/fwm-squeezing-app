"""
Phase-1 physics regression for the absorption cluster (OD / AT / EIT / CPT).

Parameter-free / analytic checks (no frozen baseline needed):
  OD cold  : weak-probe line FWHM = natural Γ.
  OD vapor : Doppler-broadened Voigt FWHM ≈ analytic Gaussian Doppler FWHM.
  AT       : probe doublet splitting = coupling Rabi Ω_c.
  EIT      : strong transparency at two-photon resonance (cold).
  CPT      : sub-natural dark resonance.

    python tests/test_absorption.py    # or: pytest tests/test_absorption.py
"""
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from gabes import schemes, observables, constants  # noqa: E402

G = constants.GAMMA
GMHZ = G / (2 * np.pi) / 1e6
K = constants.K_VEC


def _curve(raw):
    alpha, xp = observables.absorption_coefficient(
        raw["chi_bar"], K, raw["N"], line_strength=raw["ls"])
    x = raw["scan"] / (2 * np.pi) / 1e6
    return x, alpha, observables.transmission(alpha, raw["L"]), xp


def _fwhm(x, y):
    pk = np.nanmax(y)
    above = x[y >= 0.5 * pk]
    return above.max() - above.min()


def _params(scheme, **over):
    p = scheme.defaults()
    p.update(over)
    return p


def test_od_cold_natural_linewidth():
    od = schemes.get("od")
    x, alpha, _, _ = _curve(od.compute(_params(od, temp_c=25, cell_mm=3, doppler="off")))
    assert 0.9 <= _fwhm(x, alpha) / GMHZ <= 1.1


def test_od_doppler_voigt_width():
    od = schemes.get("od")
    raw = od.compute(_params(od, temp_c=50, cell_mm=10, doppler="on"))
    x, alpha, _, _ = _curve(raw)
    analytic = np.sqrt(8 * np.log(2)) * K * raw["sigma_v"] / (2 * np.pi) / 1e6
    ratio = _fwhm(x, alpha) / analytic
    assert 0.98 <= ratio <= 1.06            # Voigt ≥ Gaussian, slightly wider


def test_at_splitting_equals_coupling_rabi():
    at = schemes.get("at")
    for Oc in (6.0, 8.0, 12.0):
        x, alpha, _, _ = _curve(at.compute(_params(at, coupling_rabi=Oc, doppler="off")))
        xl = x[x < 0][np.argmax(alpha[x < 0])]
        xr = x[x > 0][np.argmax(alpha[x > 0])]
        assert abs((xr - xl) / (Oc * GMHZ) - 1.0) < 0.05


def test_eit_transparency():
    eit = schemes.get("eit")
    raw = eit.compute(_params(eit, temp_c=25, cell_mm=3, doppler="off", coupling_rabi=3.0))
    x, _, T, _ = _curve(raw)
    ic = int(np.argmin(np.abs(x)))
    assert T[ic] > 0.8 and T[ic] > 10 * T.min()


def test_cpt_subnatural_dark_resonance():
    cpt = schemes.get("cpt")
    raw = cpt.compute(_params(cpt, doppler="off"))
    x, _, T, _ = _curve(raw)
    ic = int(np.argmin(np.abs(x)))
    peak, floor = T[ic], T.min()
    thr = 0.5 * (peak + floor)
    i, j = ic, ic
    while i > 0 and T[i] >= thr:
        i -= 1
    while j < T.size - 1 and T[j] >= thr:
        j += 1
    dark_fwhm_mhz = x[j] - x[i]
    assert T[ic] > 0.7
    assert dark_fwhm_mhz < GMHZ             # sub-natural


if __name__ == "__main__":
    test_od_cold_natural_linewidth()
    test_od_doppler_voigt_width()
    test_at_splitting_equals_coupling_rabi()
    test_eit_transparency()
    test_cpt_subnatural_dark_resonance()
    print("Phase-1 absorption physics OK (OD/AT/EIT/CPT).")
