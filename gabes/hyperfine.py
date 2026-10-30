"""
85Rb D1 hyperfine structure for the realistic multi-line OD spectrum.

The atomic data were ported from the lab AutoOD calculator. Absorption
normalization was corrected on 6 September 2026 to count each thermal
magnetic-sublevel population once; software agreement is not experimental
validation. Ported pieces:

  * the four D1 hyperfine transitions  Fg ∈ {2,3} → Fe ∈ {2,3}  with their
    line-center detunings (referenced to the 87Rb F=2→F'=2 marker — the
    AutoOD x-axis origin, so GABES overlays the lab tool and its CSVs);
  * the Clebsch-Gordan line-strength factors C_F²  (Wigner 3j/6j sum);
  * the ground-state population weights  p_F = (2F+1)/Σ(2F+1);
  * self-broadening  Γ_eff = Γ_nat + β·N;
  * a pure-85Rb cell density (f85 = 1) from the CRC vapor pressure.

This module is *data only* and is the ⁸⁵Rb D1 view of `species`, which is
the single source of the vapor pressure, density and β. Kept here: the exact
rational C_F² table, the AutoOD-origin line shifts, and Γ from the legacy
`constants.GAMMA` (5.746 MHz, not the species record's 5.750 MHz).
The Voigt lineshape itself is still produced by the OBE Doppler kernel in
`schemes/absorption.py` (a genuine OBE solve, no analytic wofz / scipy);
hyperfine.py supplies the validated line table, the relative/absolute
strengths, the self-broadened width, and the density.

Reference
  [1] AutoOD lab calculator (rb_od_calculator_optimized_final.py).
  [2] Daniel A. Steck, "Rubidium 85 D Line Data," http://steck.us/alkalidata.
"""
from . import constants, species

# ---- Clebsch-Gordan line-strength factors C_F²  (5S1/2 Fg → 5P1/2 Fe) ----
# Exact rationals from the Wigner 3j/6j sum AutoOD evaluates with sympy
# (verified: (2,2)=10/81, (2,3)=35/81, (3,2)=35/81, (3,3)=28/81).
CF2 = {
    (2, 2): 10.0 / 81.0,
    (2, 3): 35.0 / 81.0,
    (3, 2): 35.0 / 81.0,
    (3, 3): 28.0 / 81.0,
}

# ---- Ground-state population weights  p_F = (2F+1) / Σ(2F+1) ----
# Σ over F ∈ {2,3} = 5 + 7 = 12  (= 2(2I+1), I = 5/2).
GROUND_POP = {2: 5.0 / 12.0, 3: 7.0 / 12.0}
N_GROUND_SUBLEVELS = 12                       # 2(2I+1)

# ---- Line-center detunings [Hz], referenced to the 87Rb F=2→F'=2 marker ----
# AutoOD origin: de8587 = 1065.646 MHz between that marker and 85Rb F=3→F'=3.
# Ground HF = 3035.732 MHz (F=2 above F=3 by that), excited HF (D1) = 361.58 MHz.
_DE8587 = 1065.646e6
_GROUND_HF = constants.NU_GROUND_HF           # 3.035732439e9
_EXC_HF = constants.NU_EXCITED_HF_D1          # 361.58e6
LINE_SHIFT_HZ = {
    (3, 3): _DE8587,
    (3, 2): _DE8587 - _EXC_HF,
    (2, 3): _DE8587 + _GROUND_HF,
    (2, 2): _DE8587 + _GROUND_HF - _EXC_HF,
}

# Iterate in a fixed (low → high detuning) order for stable plotting / metrics.
TRANSITIONS = ((3, 2), (3, 3), (2, 2), (2, 3))

# ---- Effective |d|² in the AutoOD normalization ----
# AutoOD's reduced dipole d² = 3·DIPOLE_D1² (the J=1/2→J=1/2 reduced matrix
# element times the structural Wigner-6j J→L reduction; matches the AutoOD
# value 1.9319e-57 C²m² to < 1e-4).
DIPOLE_SQ = 3.0 * constants.DIPOLE_D1 ** 2

# ---- Self-broadening:  Γ_eff = Γ_nat + β·N,  β/2π = 0.69e-7 Hz·cm³ ----
BETA_SELF = species.BETA_SELF                 # rad·s⁻¹·m³

vapor_pressure_pa = species.rb_vapor_pressure_pa


def number_density(T):
    """
    85Rb number density [/m³] for a *pure-85Rb* cell (f85 = 1), from the CRC
    vapor pressure — matches the AutoOD absolute scale. (Distinct from
    species.rb85_steck_density, the Steck liquid-phase fit at natural
    abundance used by the Λ/Rydberg Generic medium, OD single-line and SAS
    Generic paths.)
    """
    return species.number_density(species.RB85, T)


def self_broadened_gamma(N):
    """Self-broadened optical linewidth  Γ_eff = Γ_nat + β·N  [rad/s]."""
    return constants.GAMMA + species.self_broadening(species.RB85, N)
