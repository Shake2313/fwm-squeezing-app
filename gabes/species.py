"""
Alkali species / D-line data and a saturated-absorption manifold builder.

This is the atomic-data layer the strengthened SAS scheme stands on. Everything
here is *data + pure helpers*; the OBE solve lives in `core`, the experiment in
`schemes/sas.py`.

What it provides
  * Per-isotope hyperfine constants (magnetic-dipole A, electric-quadrupole B)
    for the 5S/6S ground and the nP_{1/2} (D1) / nP_{3/2} (D2) excited states,
    optical line-centre (centroid) frequencies, masses, natural linewidths and
    isotopic abundances — for ⁸⁵Rb, ⁸⁷Rb and ¹³³Cs.
  * The Casimir hyperfine-energy formula → every F level placed exactly from A, B
    (no per-line hard-coded tables; self-validating against the known Rb ground /
    excited splittings already in `constants`/`hyperfine`).
  * Relative hyperfine line strengths S(Fg→Fe) from a Wigner-6j (Racah formula),
    which set (i) the weak-probe line weights, (ii) the relative pump Rabi
    frequencies and (iii) the spontaneous-emission branching ratios — the last
    being what drives hyperfine optical pumping (the inverted crossovers).
  * `build_manifold(isotope, line)` → an `atoms.AtomModel` carrying the full
    {Fg}×{Fe} level scheme with CG-branched decay *and* a transit-time relaxation
    toward the thermal ground distribution (regularises the pumping), plus the
    transition table the scheme needs.

References
  [1] Daniel A. Steck, "Rubidium 85 / 87 D Line Data" and "Cesium D Line Data,"
      http://steck.us/alkalidata. (Hyperfine A/B, centroids, masses, linewidths.)
  [2] D. A. Smith and I. G. Hughes, "The role of hyperfine pumping in multilevel
      systems exhibiting saturated absorption," Am. J. Phys. 72, 631 (2004).
      (Multilevel SAS, hyperfine pumping → enhanced/inverted crossovers, transit.)
  [3] K. B. MacAdam, A. Steinbach, C. Wieman, Am. J. Phys. 60, 1098 (1992) and
      D. W. Preston, Am. J. Phys. 64, 1432 (1996). (Rb saturated-absorption,
      crossover bookkeeping.)
"""
import functools
import math
from dataclasses import dataclass, field

import numpy as np

from . import constants, hyperfine
from .atoms import AtomModel
from .zeeman import clebsch_gordan

MHZ = 2 * np.pi * 1e6          # MHz → rad/s
_AMU = 1.660_539_066_60e-27    # kg
_LG, _LE, _SE = 0, 1, 0.5      # alkali D line: nS → nP orbital, electron spin ½

# Self-broadening coefficient β/2π = 0.69e-7 Hz·cm³ (Rb; AutoOD / hyperfine.py).
BETA_SELF = 2 * np.pi * 0.69e-7 * 1e-6        # rad·s⁻¹·m³


# =====================================================================
# Angular-momentum algebra (Wigner 6j via the Racah formula)
# =====================================================================
def _fac(n):
    return math.factorial(int(round(n)))


def _is_int(x):
    return abs(x - round(x)) < 1e-9


def _tri_ok(a, b, c):
    """Triangle condition |a−b| ≤ c ≤ a+b with integer perimeter."""
    return (c >= abs(a - b) - 1e-9 and c <= a + b + 1e-9
            and _is_int(a + b + c))


def _delta(a, b, c):
    return math.sqrt(_fac(a + b - c) * _fac(a - b + c) * _fac(-a + b + c)
                     / _fac(a + b + c + 1))


def wigner6j(j1, j2, j3, j4, j5, j6):
    """{j1 j2 j3 ; j4 j5 j6} via the Racah W-coefficient sum. 0 if non-physical."""
    for (a, b, c) in ((j1, j2, j3), (j1, j5, j6), (j4, j2, j6), (j4, j5, j3)):
        if not _tri_ok(a, b, c):
            return 0.0
    pref = _delta(j1, j2, j3) * _delta(j1, j5, j6) * _delta(j4, j2, j6) * _delta(j4, j5, j3)
    t_lo = max(j1 + j2 + j3, j1 + j5 + j6, j4 + j2 + j6, j4 + j5 + j3)
    t_hi = min(j1 + j2 + j4 + j5, j2 + j3 + j5 + j6, j3 + j1 + j6 + j4)
    s = 0.0
    t = int(math.ceil(t_lo - 1e-9))
    while t <= int(math.floor(t_hi + 1e-9)):
        denom = (_fac(t - j1 - j2 - j3) * _fac(t - j1 - j5 - j6)
                 * _fac(t - j4 - j2 - j6) * _fac(t - j4 - j5 - j3)
                 * _fac(j1 + j2 + j4 + j5 - t) * _fac(j2 + j3 + j5 + j6 - t)
                 * _fac(j3 + j1 + j6 + j4 - t))
        s += (-1) ** t * _fac(t + 1) / denom
        t += 1
    return pref * s


def line_strength(Fg, Fe, I, Jg, Je):
    """
    Relative hyperfine transition strength S(Fg→Fe) (summed over m sublevels):
        S = (2Fe+1)(2Jg+1) { Jg Je 1 ; Fe Fg I }².
    Proportional to |⟨Fe‖er‖Fg⟩|²; sets line weights, Rabi² and decay branching.
    (Validated against the ⁸⁵Rb D1 Clebsch-Gordan factors in `hyperfine.CF2`.)
    """
    sixj = wigner6j(Jg, Je, 1, Fe, Fg, I)
    return (2 * Fe + 1) * (2 * Jg + 1) * sixj * sixj


# =====================================================================
# Hyperfine energies (Casimir formula) and F manifolds
# =====================================================================
def hf_energy_mhz(A, B, I, J, F):
    """Hyperfine shift of level F from the fine-structure centroid [MHz]."""
    K = F * (F + 1) - I * (I + 1) - J * (J + 1)
    E = 0.5 * A * K
    if B and (2 * I - 1) > 0 and (2 * J - 1) > 0:
        E += B * (1.5 * K * (K + 1) - 2 * I * (I + 1) * J * (J + 1)) \
            / (2 * I * (2 * I - 1) * 2 * J * (2 * J - 1))
    return E


def f_values(I, J):
    """Allowed F = |I−J| … I+J (integer steps)."""
    lo = abs(I - J)
    hi = I + J
    n = int(round(hi - lo))
    return [lo + k for k in range(n + 1)]


# =====================================================================
# Isotope data (Steck D-line data sheets)
# =====================================================================
@dataclass(frozen=True)
class Isotope:
    name: str
    label: str            # display label
    I: float              # nuclear spin
    mass: float           # kg
    abundance: float      # natural fractional abundance
    A_S: float            # ground nS_{1/2}  magnetic-dipole A [MHz]
    A_P12: float          # nP_{1/2} (D1 excited) A [MHz]
    A_P32: float          # nP_{3/2} (D2 excited) A [MHz]
    B_P32: float          # nP_{3/2} electric-quadrupole B [MHz]
    nu_D1: float          # D1 line-centre (centroid) frequency [Hz]
    nu_D2: float          # D2 line-centre (centroid) frequency [Hz]
    gamma_D1_mhz: float   # natural linewidth Γ/2π on D1 [MHz]
    gamma_D2_mhz: float   # natural linewidth Γ/2π on D2 [MHz]
    Jg: float = 0.5

    def line(self, which):
        """(Je, centroid ν [Hz], Γ/2π [MHz], A_excited, B_excited) for 'D1'/'D2'."""
        if which == "D1":
            return 0.5, self.nu_D1, self.gamma_D1_mhz, self.A_P12, 0.0
        return 1.5, self.nu_D2, self.gamma_D2_mhz, self.A_P32, self.B_P32


RB85 = Isotope(
    name="Rb85", label="⁸⁵Rb", I=2.5, mass=84.911_789_738 * _AMU, abundance=0.7217,
    A_S=1011.910_813, A_P12=120.527, A_P32=25.0020, B_P32=25.790,
    nu_D1=377.107_385_690e12, nu_D2=384.230_406_373e12,
    gamma_D1_mhz=5.7500, gamma_D2_mhz=6.0666)

RB87 = Isotope(
    name="Rb87", label="⁸⁷Rb", I=1.5, mass=86.909_180_527 * _AMU, abundance=0.2783,
    A_S=3417.341_306, A_P12=407.24, A_P32=84.7185, B_P32=12.4965,
    nu_D1=377.107_463_380e12, nu_D2=384.230_484_468e12,
    gamma_D1_mhz=5.7500, gamma_D2_mhz=6.0666)

CS133 = Isotope(
    name="Cs133", label="¹³³Cs", I=3.5, mass=132.905_451_961 * _AMU, abundance=1.0,
    A_S=2298.157_943, A_P12=291.9201, A_P32=50.275, B_P32=-0.53,
    nu_D1=335.116_048_807e12, nu_D2=351.725_718_501e12,
    gamma_D1_mhz=4.5612, gamma_D2_mhz=5.2227)

ISOTOPES = {iso.name: iso for iso in (RB85, RB87, CS133)}

# Species menu → list of (isotope, weight). Natural Rb overlays both isotopes
# (weighted by abundance); single isotopes use weight 1 (their own density).
SPECIES = {
    "Rb (natural)": [(RB85, RB85.abundance), (RB87, RB87.abundance)],
    "⁸⁵Rb": [(RB85, 1.0)],
    "⁸⁷Rb": [(RB87, 1.0)],
    "¹³³Cs": [(CS133, 1.0)],
}
SPECIES_ORDER = ["Rb (natural)", "⁸⁵Rb", "⁸⁷Rb", "¹³³Cs"]


# =====================================================================
# Number density (Steck vapor-pressure fits; liquid phase, warm cell)
# =====================================================================
def _density_from_logP(log10_P_torr, T):
    P_pa = 10 ** log10_P_torr * 133.322_387_415
    return P_pa / (constants.KB * T)


def number_density(iso, T):
    """
    *Elemental* atomic number density [/m³] at cell temperature T [K]. Total
    (all-isotope) density; the per-isotope share is applied by the caller via the
    species weight (abundance for a natural cell, 1.0 for an isotopically pure one).

    Rb uses the CRC vapor pressure (solid below 39.30 °C, liquid above) — the same
    fit the AutoOD-validated OD scheme uses, so the pump-off limit matches the lab
    tool's absolute scale. Cs uses the Steck liquid-phase fit.
    """
    if iso.name.startswith("Rb"):
        return hyperfine.vapor_pressure_pa(T) / (constants.KB * T)
    logP = 8.22127 - 4006.048 / T - 0.00060194 * T - 0.19623 * np.log10(T)
    return _density_from_logP(logP, T)


def self_broadened_gamma(iso, N):
    """Self-broadened optical linewidth Γ_eff = Γ_nat + β·N [rad/s] (β: Rb value)."""
    return BETA_SELF * N


# =====================================================================
# Absolute line strengths (AutoOD convention) and the pump-power map
# =====================================================================
def _wigner3j(j1, j2, j3, m1, m2, m3):
    if abs(m1 + m2 + m3) > 1e-9:
        return 0.0
    return ((-1) ** int(round(j1 - j2 - m3)) / math.sqrt(2 * j3 + 1)
            * clebsch_gordan(j1, m1, j2, m2, j3, -m3))


def cf2(Fg, Fe, I, Jg, Je):
    """
    Relative hyperfine line strength C_F² in the AutoOD convention (validated:
    reproduces hyperfine.CF2 for ⁸⁵Rb D1). Proportional to (2Fg+1)·line_strength
    but computed the lab tool's way for exact absolute-scale agreement.
    """
    m_max = int(min(Fg, Fe))
    cfsq = sum(_wigner3j(Fe, 1, Fg, m, 0, -m) ** 2 for m in range(-m_max, m_max + 1))
    six1 = abs(wigner6j(Jg, Je, 1, Fe, Fg, I))
    six2 = abs(wigner6j(_LG, _LE, 1, Je, Jg, _SE))
    pref = (2 * Fg + 1) * (2 * Fe + 1) * (2 * Jg + 1) * (2 * Je + 1) * (2 * _LG + 1)
    return pref * cfsq * six1 ** 2 * six2 ** 2


def reduced_dipole_sq(gamma_nat, lam, Jg, Je):
    """|⟨J‖er‖J′⟩|² [C²m²] from the natural linewidth (AutoOD d_func). The six2
    L→J reduction here cancels the one inside C_F², so C_F²·d² is six2-free."""
    six2 = abs(wigner6j(_LG, _LE, 1, Je, Jg, _SE))
    num = 3 * constants.EPS_0 * constants.HBAR * lam ** 3 * gamma_nat
    den = 8 * np.pi ** 2 * six2 ** 2 * (2 * Jg + 1) * (2 * _LG + 1)
    return num / den


def line_integrated_alpha(iso, line, N):
    """
    Per-transition integrated weak-probe absorption ∫α dδ [rad·s⁻¹·m⁻¹] at line
    strength 1 (AutoOD normalisation): ∫α_t = π·k·p_F·C_F²·|d|²·N/(ε₀ℏ)/(2(2I+1)).
    Returns a dict keyed by (Fg, Fe). Density N already carries the species weight.
    """
    Je, nu0, gamma_mhz, _, _ = iso.line(line)
    I, Jg = iso.I, iso.Jg
    lam = constants.C_LIGHT / nu0
    k = 2 * np.pi / lam
    gamma_nat = gamma_mhz * MHZ
    d2 = reduced_dipole_sq(gamma_nat, lam, Jg, Je)
    deg = {F: 2 * F + 1 for F in f_values(I, Jg)}
    ptot = sum(deg.values())
    K = np.pi * k * d2 * N / (constants.HBAR * constants.EPS_0) / (2 * (2 * I + 1))
    out = {}
    for Fg in f_values(I, Jg):
        for Fe in f_values(I, Je):
            S = cf2(Fg, Fe, I, Jg, Je)
            if S > 1e-12:
                out[(Fg, Fe)] = K * (deg[Fg] / ptot) * S
    return out


def pump_rabi_from_power(power_mw, waist_mm, gamma):
    """
    Pump Rabi Ω [rad/s] from beam power and 1/e² waist: I = 2P/(πw²),
    Ω = Γ·√(I/2I_sat) with I_sat = 4.484 mW/cm² (the reference saturation
    intensity; documented in the scheme About). Power 0 → Ω 0 → linear OD limit.
    """
    if power_mw <= 0:
        return 0.0
    w = waist_mm * 1e-3
    I = 2 * (power_mw * 1e-3) / (np.pi * w ** 2)
    return gamma * np.sqrt(I / (2 * constants.I_SAT))


# Per-(species, line) recommended slider values for the "Default" button —
# tuned so the pump-off Doppler dips are clearly visible (peak OD ~0.4–1).
_RECOMMENDED = {
    ("Rb (natural)", "D1"): dict(temp_c=40.0, cell_mm=75.0, pump_power_mw=0.5),
    ("Rb (natural)", "D2"): dict(temp_c=25.0, cell_mm=75.0, pump_power_mw=0.5),
    ("⁸⁵Rb", "D1"): dict(temp_c=40.0, cell_mm=75.0, pump_power_mw=0.5),
    ("⁸⁵Rb", "D2"): dict(temp_c=25.0, cell_mm=75.0, pump_power_mw=0.5),
    ("⁸⁷Rb", "D1"): dict(temp_c=45.0, cell_mm=75.0, pump_power_mw=0.5),
    ("⁸⁷Rb", "D2"): dict(temp_c=30.0, cell_mm=75.0, pump_power_mw=0.5),
    ("¹³³Cs", "D1"): dict(temp_c=30.0, cell_mm=50.0, pump_power_mw=0.5),
    ("¹³³Cs", "D2"): dict(temp_c=22.0, cell_mm=30.0, pump_power_mw=0.5),
}


def recommended(species_key, line):
    """Recommended (temp_c, cell_mm, pump_power_mw) for a species/line, or a default."""
    return _RECOMMENDED.get((species_key, line),
                            dict(temp_c=30.0, cell_mm=50.0, pump_power_mw=0.5))


# =====================================================================
# Manifold builder
# =====================================================================
@dataclass
class Manifold:
    """A built (isotope, line) saturated-absorption level scheme + transition table."""
    iso: Isotope
    line: str
    atom: AtomModel
    Jg: float
    Je: float
    Fg: list                         # ground F values (index order)
    Fe: list                         # excited F values (index order)
    e_ground: np.ndarray             # ground F energies [rad/s], centroid-referenced
    e_excited: np.ndarray            # excited F energies [rad/s]
    p_ground: np.ndarray             # thermal ground populations (2F+1)/Σ
    omega: np.ndarray                # per-transition detuning δe−δg [rad/s]  (n_trans,)
    g_idx: np.ndarray                # ground level index per transition
    e_idx: np.ndarray                # excited level index per transition
    S_abs: np.ndarray                # absorption weight per transition (Σ S·p = 1)
    rabi_rel: np.ndarray             # pump Rabi scale per transition (√(S/Smax))
    gamma: float                     # Γ [rad/s] on this line
    k_vec: float                     # optical wavenumber [1/m]
    nu0: float                       # centroid frequency [Hz]
    H_couplings: list = field(default_factory=list)   # (g_idx, e_idx, rabi_rel)

    @property
    def n_levels(self):
        return self.atom.n_levels


@functools.lru_cache(maxsize=64)
def build_manifold(iso, line, transit_rate=0.0):
    """
    Build the full {Fg}×{Fe} hyperfine manifold for `iso` on `line` ('D1'/'D2').

    Decay is CG-branched: Γ_{Fe→Fg} = Γ · S(Fg,Fe)/Σ_{Fg'} S(Fg',Fe) (Σ_Fg = Γ),
    so an excited state preferentially decays toward the strongest ground link —
    the mechanism that pumps population into the *other* ground hyperfine state.
    A transit-time relaxation moves every state toward the thermal ground
    distribution at `transit_rate` [rad/s] (atoms leaving/entering the beam),
    which regularises the otherwise-runaway optical pumping.
    """
    Je, nu0, gamma_mhz, A_e, B_e = iso.line(line)
    I, Jg = iso.I, iso.Jg
    gamma = gamma_mhz * MHZ

    Fg = f_values(I, Jg)
    Fe = f_values(I, Je)
    ng, ne = len(Fg), len(Fe)
    n = ng + ne
    ground = tuple(range(ng))
    excited = tuple(range(ng, n))

    e_ground = np.array([hf_energy_mhz(iso.A_S, 0.0, I, Jg, F) for F in Fg]) * MHZ
    e_excited = np.array([hf_energy_mhz(A_e, B_e, I, Je, F) for F in Fe]) * MHZ

    deg = np.array([2 * F + 1 for F in Fg], dtype=float)
    p_ground = deg / deg.sum()

    # Relative hyperfine strengths.  The Wigner-6j gives |⟨Fe‖d‖Fg⟩|² ∝
    # line_strength; the *observable* relative strength of a lumped Fg↔Fe line
    # (absorption per ground atom and the spontaneous-emission branching) carries
    # an extra ground-degeneracy factor (2Fg+1):
    #     T(Fg,Fe) = (2Fg+1)·line_strength,   Σ_Fg T = 2Fe+1,  Σ_Fe T = 2Fg+1.
    # T reproduces the AutoOD-validated ⁸⁵Rb D1 factors (T = 9·hyperfine.CF2) and
    # the 49/25 F=3/F=2 manifold absorption ratio. Used for line weight, decay
    # branching and (√T) the relative pump Rabi — one consistent quantity.
    S = np.array([[line_strength(fg, fe, I, Jg, Je) for fe in Fe] for fg in Fg])
    T = deg[:, None] * S
    T_max = T.max()

    # CG-branched spontaneous emission: Γ_{Fe→Fg} = Γ·T/Σ_{Fg'}T = Γ·T/(2Fe+1).
    decay = []
    for ie, fe in enumerate(Fe):
        col = T[:, ie]
        tot = col.sum()
        for ig, fg in enumerate(Fg):
            if col[ig] > 0:
                decay.append((ng + ie, ig, gamma * col[ig] / tot))

    # Transit-time relaxation toward the thermal ground distribution: every state
    # s → ground g at rate γ_t·p(g) (drives ρ → Σ_g p_g|g⟩⟨g| with the laser off).
    if transit_rate > 0:
        for s in range(n):
            for ig in range(ng):
                decay.append((s, ig, transit_rate * p_ground[ig]))

    atom = AtomModel(
        name=f"sas_{iso.name}_{line}", n_levels=n,
        labels=tuple(f"g{F:g}" for F in Fg) + tuple(f"e{F:g}" for F in Fe),
        ground=ground, excited=excited,
        decay=tuple(decay), dephasing=(), doppler_levels=excited,
    )

    # Transition table (every allowed Fg→Fe).
    g_idx, e_idx, omega, S_abs, rabi_rel, H_couplings = [], [], [], [], [], []
    norm = float((T * p_ground[:, None]).sum())      # Σ T·p → 1 (per-atom strength)
    for ig, fg in enumerate(Fg):
        for ie, fe in enumerate(Fe):
            if T[ig, ie] <= 0:
                continue
            g_idx.append(ig)
            e_idx.append(ng + ie)
            omega.append(e_excited[ie] - e_ground[ig])
            S_abs.append(T[ig, ie] / norm)
            rr = math.sqrt(T[ig, ie] / T_max)
            rabi_rel.append(rr)
            H_couplings.append((ig, ng + ie, rr))

    return Manifold(
        iso=iso, line=line, atom=atom, Jg=Jg, Je=Je, Fg=Fg, Fe=Fe,
        e_ground=e_ground, e_excited=e_excited, p_ground=p_ground,
        omega=np.array(omega), g_idx=np.array(g_idx), e_idx=np.array(e_idx),
        S_abs=np.array(S_abs), rabi_rel=np.array(rabi_rel),
        gamma=gamma, k_vec=2 * np.pi * nu0 / constants.C_LIGHT, nu0=nu0,
        H_couplings=H_couplings)


def pump_hamiltonian(man, omega_pump):
    """Rotating-frame pump H₀ at Δ_eff=0: ground/excited HF on the diagonal,
    Ω_pump·√(S/Smax) on each allowed link (the scan enters later via Δ_eff·S_v)."""
    n = man.n_levels
    H = np.zeros((n, n), dtype=complex)
    for ig in range(len(man.Fg)):
        H[ig, ig] = man.e_ground[ig]
    for ie in range(len(man.Fe)):
        e = len(man.Fg) + ie
        H[e, e] = man.e_excited[ie]
    for (g, e, rr) in man.H_couplings:
        H[g, e] = H[e, g] = omega_pump * rr / 2
    return H
