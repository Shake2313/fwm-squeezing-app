"""
Zeeman manifold support for the magneto-optics schemes (Hanle / EIA / NMOR).

Hand-rolled Clebsch-Gordan (Racah formula — no sympy dependency) and a builder
that turns (F_g, F_e) into an `AtomModel` plus the polarization-resolved dipole
couplings (σ⁺, σ⁻, π) and CG-branched spontaneous emission.
"""
import math

import numpy as np

from . import constants
from .atoms import AtomModel


def clebsch_gordan(j1, m1, j2, m2, J, M):
    """⟨j1 m1 j2 m2 | J M⟩ via the Racah formula. Works for integer/half-integer."""
    if m1 + m2 != M:
        return 0.0
    if J < abs(j1 - j2) or J > j1 + j2:
        return 0.0
    if abs(m1) > j1 or abs(m2) > j2 or abs(M) > J:
        return 0.0

    f = math.factorial
    pref = math.sqrt(
        (2 * J + 1)
        * f(int(j1 + j2 - J)) * f(int(j1 - j2 + J)) * f(int(-j1 + j2 + J))
        / f(int(j1 + j2 + J + 1))
    )
    pref *= math.sqrt(
        f(int(J + M)) * f(int(J - M))
        * f(int(j1 - m1)) * f(int(j1 + m1))
        * f(int(j2 - m2)) * f(int(j2 + m2))
    )

    kmin = int(max(0, j2 - m1 - J, j1 + m2 - J))
    kmax = int(min(j1 + j2 - J, j1 - m1, j2 + m2))
    s = 0.0
    for k in range(kmin, kmax + 1):
        denom = (f(k) * f(int(j1 + j2 - J - k)) * f(int(j1 - m1 - k))
                 * f(int(j2 + m2 - k)) * f(int(J - j2 + m1 + k)) * f(int(J - j1 - m2 + k)))
        s += (-1) ** k / denom
    return pref * s


def zeeman_manifold(Fg, Fe, gamma=None, gamma_gg=None, g_ratio=1.0):
    """
    Build the (F_g ↔ F_e) Zeeman manifold as an AtomModel with extra attributes:
      m_ground, m_excited : magnetic quantum numbers per level index
      g_ratio             : excited/ground Landé g-factor ratio (excited Zeeman)
      couplings[q]        : list of (ground_idx, excited_idx, CG) for q ∈ {−1,0,+1}
    Spontaneous emission is CG-branched (Γ·|CG|² per channel, Σ = Γ from each
    excited); ground-ground Zeeman coherences dephase at γ_gg. All excited levels
    are Doppler-shifted (optical line); ground Zeeman is Doppler-free.
    """
    gamma = constants.GAMMA if gamma is None else gamma
    gamma_gg = constants.GAMMA_GG if gamma_gg is None else gamma_gg

    mg = list(range(-Fg, Fg + 1))
    me = list(range(-Fe, Fe + 1))
    ng, ne = len(mg), len(me)
    n = ng + ne
    ground = tuple(range(ng))
    excited = tuple(range(ng, n))

    decay = []
    couplings = {-1: [], 0: [], +1: []}
    for ie, m_e in enumerate(me):
        eidx = ng + ie
        for ig, m_g in enumerate(mg):
            q = m_e - m_g
            if q not in (-1, 0, 1):
                continue
            cg = clebsch_gordan(Fg, m_g, 1, q, Fe, m_e)
            if abs(cg) < 1e-14:
                continue
            decay.append((eidx, ig, gamma * cg * cg))   # emission |Fe m_e⟩→|Fg m_g⟩
            couplings[q].append((ig, eidx, cg))          # drive ground→excited, pol q

    dephasing = tuple((i, j, gamma_gg) for i in ground for j in ground if i != j)

    atom = AtomModel(
        name=f"zeeman_{Fg}_{Fe}", n_levels=n,
        labels=tuple(f"g{m:+d}" for m in mg) + tuple(f"e{m:+d}" for m in me),
        ground=ground, excited=excited,
        decay=tuple(decay), dephasing=dephasing, doppler_levels=excited,
    )
    atom.m_ground = np.array(mg, dtype=float)
    atom.m_excited = np.array(me, dtype=float)
    atom.g_ratio = g_ratio
    atom.couplings = couplings
    return atom
