"""
Zeeman manifold support for the magneto-optics schemes (Hanle / EIA / NMOR).

Hand-rolled Clebsch-Gordan (Racah formula — no sympy dependency) and a builder
that turns (F_g, F_e) into an `AtomModel` plus the polarization-resolved dipole
couplings (σ⁺, σ⁻, π) and CG-branched spontaneous emission.
"""
import math

import numpy as np

from . import constants, hyperfine
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


def angular_momentum_matrices(F):
    """Return (Fx, Fy, Fz) in the |F,m> basis ordered m=-F..F.

    Matrices are dimensionless: H_B = omega_L * (Bxhat Fx + Byhat Fy + Bzhat Fz).
    """
    mvals = np.arange(-F, F + 1, 1, dtype=float)
    n = mvals.size
    Fp = np.zeros((n, n), dtype=complex)
    for col, m in enumerate(mvals):
        mp = m + 1
        if mp > F:
            continue
        row = int(round(mp + F))
        Fp[row, col] = math.sqrt(F * (F + 1) - m * mp)
    Fm = Fp.conj().T
    Fx = 0.5 * (Fp + Fm)
    Fy = (Fp - Fm) / (2j)
    Fz = np.diag(mvals.astype(complex))
    return Fx, Fy, Fz


def zeeman_manifold(Fg, Fe, gamma=None, gamma_gg=None, g_ratio=1.0,
                    transit_rate=0.0):
    """
    Build the (F_g ↔ F_e) Zeeman manifold as an AtomModel with extra attributes:
      m_ground, m_excited : magnetic quantum numbers per level index
      g_ratio             : excited/ground Landé g-factor ratio (excited Zeeman)
      couplings[q]        : list of (ground_idx, excited_idx, CG) for q ∈ {-1,0,+1}
    Spontaneous emission is CG-branched (Γ·|CG|² per channel, Σ = Γ from each
    excited); ground-ground Zeeman coherences dephase at γ_gg. All excited levels
    are Doppler-shifted (optical line); ground Zeeman is Doppler-free.

    `transit_rate` optionally adds isotropic population reload into the addressed
    ground manifold. It is a compact warm-vapor approximation for atoms entering
    and leaving the beam, useful for practical Hanle/EIA curves where optical
    pumping would otherwise trap all population.
    """
    gamma = constants.GAMMA if gamma is None else gamma
    gamma_gg = constants.GAMMA_GG if gamma_gg is None else gamma_gg

    mg = list(range(-Fg, Fg + 1))
    me = list(range(-Fe, Fe + 1))
    ng, ne = len(mg), len(me)
    n = ng + ne
    ground = tuple(range(ng))
    excited = tuple(range(ng, n))

    # Spontaneous emission as polarization-grouped jump operators Σ_q (q = m_e−m_g):
    #   Σ_q = √Γ Σ ⟨Fg m_g; 1 q | Fe m_e⟩ |Fg m_g⟩⟨Fe m_e|.
    # One operator per q (not one per channel) so that D[Σ_q] carries transfer of
    # coherence: a ground Zeeman coherence |g⟩⟨g'| is refed from an excited
    # coherence ρ_{ee'} sharing the same emitted-photon polarization. This is what
    # turns cycling-type transitions into EIA/LCA. Σ_q^†Σ_q sums to Γ·P_excited, so
    # excited decay and ground-population refilling are byte-identical to the old
    # per-channel form; only the off-diagonal (TOC) gain is new.
    sqrt_gamma = math.sqrt(gamma)
    emission = {q: np.zeros((n, n), dtype=complex) for q in (-1, 0, 1)}
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
            emission[q][ig, eidx] = sqrt_gamma * cg      # |Fg m_g⟩⟨Fe m_e|
            couplings[q].append((ig, eidx, cg))          # drive ground→excited, pol q
    emission_ops = tuple(emission[q] for q in (-1, 0, 1))

    decay = []   # incoherent population reload only (transit through the beam)
    if transit_rate and transit_rate > 0:
        p_ground = 1.0 / ng
        for src in range(n):
            for ig in ground:
                decay.append((src, ig, transit_rate * p_ground))

    dephasing = tuple((i, j, gamma_gg) for i in ground for j in ground if i != j)

    atom = AtomModel(
        name=f"zeeman_{Fg}_{Fe}", n_levels=n,
        labels=tuple(f"g{m:+d}" for m in mg) + tuple(f"e{m:+d}" for m in me),
        ground=ground, excited=excited,
        decay=tuple(decay), dephasing=dephasing, doppler_levels=excited,
        emission_ops=emission_ops,
    )
    atom.m_ground = np.array(mg, dtype=float)
    atom.m_excited = np.array(me, dtype=float)
    atom.g_ratio = g_ratio
    atom.couplings = couplings
    return atom


def rb85_d1_double_lambda_zeeman(gamma=None, gamma_gg=None):
    """Full 85Rb D1 double-lambda Zeeman manifold used by FWM Ultra diagnostics.

    Levels are ordered as all ground states F=2,3 followed by all excited states
    F'=2,3. Coupling amplitudes are normalized so the unpolarized CG sum for each
    Fg->Fe block reproduces the lumped FWM scale 3*C_F^2.
    """
    gamma = constants.GAMMA if gamma is None else gamma
    gamma_gg = constants.GAMMA_GG if gamma_gg is None else gamma_gg

    ground_blocks = [(2, list(range(-2, 3))), (3, list(range(-3, 4)))]
    excited_blocks = [(2, list(range(-2, 3))), (3, list(range(-3, 4)))]
    labels = []
    level_F = []
    level_m = []
    ground = []
    for F, mvals in ground_blocks:
        for m in mvals:
            ground.append(len(labels))
            labels.append(f"gF{F}m{m:+d}")
            level_F.append(F)
            level_m.append(m)
    excited = []
    for F, mvals in excited_blocks:
        for m in mvals:
            excited.append(len(labels))
            labels.append(f"eF{F}m{m:+d}")
            level_F.append(F)
            level_m.append(m)

    n = len(labels)
    couplings = {-1: [], 0: [], +1: []}
    block_sum = {}
    for Fg, _ in ground_blocks:
        for Fe, _ in excited_blocks:
            raw = []
            for gi in ground:
                if level_F[gi] != Fg:
                    continue
                for ei in excited:
                    if level_F[ei] != Fe:
                        continue
                    q = level_m[ei] - level_m[gi]
                    if q not in (-1, 0, 1):
                        continue
                    cg = clebsch_gordan(Fg, level_m[gi], 1, q, Fe, level_m[ei])
                    if abs(cg) > 1e-14:
                        raw.append((q, gi, ei, cg))
            sum_cg2 = sum(abs(cg) ** 2 for _, _, _, cg in raw)
            target = 3.0 * hyperfine.CF2[(Fg, Fe)]
            scale = math.sqrt(target * (2 * Fg + 1) / max(sum_cg2, 1e-30))
            for q, gi, ei, cg in raw:
                couplings[q].append((gi, ei, scale * cg))
            avg = sum(abs(scale * cg) ** 2 for _, _, _, cg in raw) / (2 * Fg + 1)
            block_sum[(Fg, Fe)] = avg

    emission = {q: np.zeros((n, n), dtype=complex) for q in (-1, 0, 1)}
    for Fe, _ in excited_blocks:
        weights = {Fg: hyperfine.CF2[(Fg, Fe)] for Fg, _ in ground_blocks}
        total_w = sum(weights.values())
        for ei in excited:
            if level_F[ei] != Fe:
                continue
            for Fg, _ in ground_blocks:
                branch = weights[Fg] / total_w
                raw = []
                for gi in ground:
                    if level_F[gi] != Fg:
                        continue
                    q = level_m[ei] - level_m[gi]
                    if q not in (-1, 0, 1):
                        continue
                    cg = clebsch_gordan(Fg, level_m[gi], 1, q, Fe, level_m[ei])
                    if abs(cg) > 1e-14:
                        raw.append((q, gi, cg))
                norm = sum(abs(cg) ** 2 for _, _, cg in raw)
                amp_scale = math.sqrt(gamma * branch / max(norm, 1e-30))
                for q, gi, cg in raw:
                    emission[q][gi, ei] = amp_scale * cg

    dephasing = tuple((i, j, gamma_gg) for i in ground for j in ground if i != j)
    atom = AtomModel(
        name="rb85_d1_double_lambda_zeeman",
        n_levels=n,
        labels=tuple(labels),
        ground=tuple(ground),
        excited=tuple(excited),
        decay=(),
        dephasing=dephasing,
        doppler_levels=tuple(excited),
        emission_ops=tuple(emission[q] for q in (-1, 0, 1)),
    )
    atom.level_F = np.array(level_F, dtype=int)
    atom.level_m = np.array(level_m, dtype=float)
    atom.couplings = couplings
    atom.lumped_strengths = block_sum
    ratios = [
        block_sum[key] / (3.0 * hyperfine.CF2[key])
        for key in hyperfine.CF2
    ]
    atom.lumped_strength_correction = float(np.mean(ratios))
    atom.ultra_status = (
        "24-level 85Rb D1 Zeeman manifold available for CG diagnostics")
    return atom
