"""
85Rb D1 hyperfine structure for the realistic multi-line OD spectrum.

The data and absolute scaling here are validated against the lab AutoOD
calculator (references/AutoOD/rb_od_calculator_optimized_final.py), which
reproduces *measured* 85Rb D1 vapor-cell transmission. Ported pieces:

  * the four D1 hyperfine transitions  Fg ∈ {2,3} → Fe ∈ {2,3}  with their
    line-center detunings (referenced to the 87Rb F=2→F'=2 marker — the
    AutoOD x-axis origin, so GABES overlays the lab tool and its CSVs);
  * the Clebsch-Gordan line-strength factors C_F²  (Wigner 3j/6j sum);
  * the ground-state population weights  p_F = (2F+1)/Σ(2F+1);
  * self-broadening  Γ_eff = Γ_nat + β·N;
  * a pure-85Rb cell density (f85 = 1) from the CRC vapor pressure.

This module is *data only*. The Voigt lineshape itself is still produced by
the OBE Doppler kernel in `schemes/absorption.py` (a genuine OBE solve, no
analytic wofz / scipy); hyperfine.py supplies the validated line table, the
relative/absolute strengths, the self-broadened width, and the density.

Reference
  [1] AutoOD lab calculator (rb_od_calculator_optimized_final.py).
  [2] Daniel A. Steck, "Rubidium 85 D Line Data," http://steck.us/alkalidata.
"""
import numpy as np

from . import constants

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
BETA_SELF = 2 * np.pi * 0.69e-7 * 1e-6        # rad·s⁻¹·m³  (Hz·cm³ → rad·s⁻¹·m³)


def vapor_pressure_pa(T):
    """Rb vapor pressure [Pa] (CRC): solid below 39.30 °C, liquid above (AutoOD)."""
    if T < 273.15 + 39.30:
        return 10 ** (9.863 - 4215.0 / T)
    return 10 ** (9.318 - 4040.0 / T)


def number_density(T):
    """
    85Rb number density [/m³] for a *pure-85Rb* cell (f85 = 1), from the CRC
    vapor pressure — matches the AutoOD absolute scale. (Distinct from
    atoms.rb85_density, which is the Steck liquid-phase fit at natural
    abundance used by the other schemes.)
    """
    return vapor_pressure_pa(T) / (constants.KB * T)


def self_broadened_gamma(N):
    """Self-broadened optical linewidth  Γ_eff = Γ_nat + β·N  [rad/s]."""
    return constants.GAMMA + BETA_SELF * N
