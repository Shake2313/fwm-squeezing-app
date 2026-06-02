"""
Phase-4 checks for Zeeman magneto-optics (Hanle / EIA / NMOR).

  CG       : hand-rolled Clebsch-Gordan matches known values + normalisation.
  manifold : CG-branched decay sums to Γ from each excited sublevel.
  Hanle    : zero-field absorption dip (dark resonance, F_e ≤ F_g).
  EIA      : zero-field absorption peak (sign flip, F_e = F_g + 1).
  NMOR     : antisymmetric rotation through zero at B = 0.

    python tests/test_magneto.py    # or: pytest tests/test_magneto.py
"""
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from gabes import schemes, zeeman, observables, constants  # noqa: E402
from gabes.constants import GAMMA, K_VEC  # noqa: E402

cg = zeeman.clebsch_gordan


def test_clebsch_gordan_known_values():
    assert abs(cg(1, 0, 1, 0, 2, 0) - np.sqrt(2 / 3)) < 1e-9
    assert abs(cg(1, 1, 1, -1, 2, 0) - np.sqrt(1 / 6)) < 1e-9
    assert abs(cg(1, 1, 1, -1, 0, 0) - 1 / np.sqrt(3)) < 1e-9
    assert abs(cg(0.5, 0.5, 0.5, -0.5, 1, 0) - 1 / np.sqrt(2)) < 1e-9
    s = sum(cg(2, m1, 1, 1 - m1, 2, 1) ** 2 for m1 in range(-2, 3))
    assert abs(s - 1.0) < 1e-9


def test_manifold_decay_normalised():
    atom = zeeman.zeeman_manifold(2, 3)
    assert atom.n_levels == 12
    for e in atom.excited:
        tot = sum(r for (ei, gi, r) in atom.decay if ei == e)
        assert abs(tot / GAMMA - 1.0) < 1e-9


def _alpha_rot(raw):
    xp = observables.chi_phys(raw["chi_p"], raw["N"], line_strength=raw["ls"])
    xm = observables.chi_phys(raw["chi_m"], raw["N"], line_strength=raw["ls"])
    x = raw["larmor"] / GAMMA
    alpha = K_VEC * np.imag(xp + xm)
    rot = 0.25 * K_VEC * raw["L"] * np.real(xp - xm)
    return x, alpha, rot


def test_hanle_zero_field_dip():
    sc = schemes.get("hanle")
    x, alpha, _ = _alpha_rot(sc.compute(sc.defaults()))
    ic = int(np.argmin(np.abs(x)))
    assert alpha[ic] < alpha[0]                # dip at B=0
    assert np.all(alpha > 0)


def test_eia_zero_field_peak():
    sc = schemes.get("eia")
    x, alpha, _ = _alpha_rot(sc.compute(sc.defaults()))
    ic = int(np.argmin(np.abs(x)))
    assert alpha[ic] > alpha[0]                # peak at B=0 (sign flip vs Hanle)


def test_nmor_zero_crossing():
    sc = schemes.get("nmor")
    x, _, rot = _alpha_rot(sc.compute(sc.defaults()))
    ic = int(np.argmin(np.abs(x)))
    iL = int(np.argmin(np.abs(x + 0.3)))
    iR = int(np.argmin(np.abs(x - 0.3)))
    assert abs(rot[ic]) < 1e-6 * max(np.abs(rot).max(), 1e-30)   # ~0 at B=0
    assert rot[iL] * rot[iR] < 0               # antisymmetric sign flip


def test_invalid_transition_handled():
    sc = schemes.get("hanle")
    p = sc.defaults(); p.update(Fg=1.0, Fe=3.0)      # |Fe−Fg| = 2, not dipole-allowed
    raw = sc.compute(p)
    assert raw["valid"] is False
    view = sc.observables(raw, p)                    # must not crash
    assert view.get("figure") is not None


if __name__ == "__main__":
    test_clebsch_gordan_known_values()
    test_manifold_decay_normalised()
    test_hanle_zero_field_dip()
    test_eia_zero_field_peak()
    test_nmor_zero_crossing()
    test_invalid_transition_handled()
    print("Phase-4 magneto OK (CG, manifold, Hanle dip, EIA peak, NMOR zero-crossing).")
