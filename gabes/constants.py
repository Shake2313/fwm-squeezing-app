"""
Physical constants, 85Rb D1 line data, and basic field helpers.

Values ported verbatim from the original fwm_obe.py §1 (no numerical change).

References
  [1] Daniel A. Steck, "Rubidium 85 D Line Data," http://steck.us/alkalidata
      (revision 2.3.4, 8 August 2025).
  [2] G. Sim, H. Kim, and H. S. Moon, Sci. Rep. 15, 7727 (2025),
      doi:10.1038/s41598-025-86479-w.
"""
import numpy as np

# ---- Fundamental constants (SI) ----
HBAR = 1.054571817e-34
KB = 1.380649e-23
C_LIGHT = 299792458.0
EPS_0 = 8.8541878128e-12
ELEMENTARY_CHARGE = 1.602176634e-19

# ---- 85Rb D1 atomic data ----
NU_D1_85RB = 377.107_385_690e12
WAVELENGTH_D1_85RB = C_LIGHT / NU_D1_85RB
GAMMA_2PI = 5.746e6
NU_GROUND_HF = 3.035732439e9
NU_EXCITED_HF_D1 = 361.58e6
NU_HF = NU_GROUND_HF
MASS_85RB = 1.4100e-25
I_SAT = 4.484e-3 * 1e4
DIPOLE_D1 = 2.5377e-29          # ⟨J=1/2‖er‖J=1/2⟩ reduced matrix element [C·m]
RB85_ABUNDANCE = 0.7217

# Effective line-strength factor (Clebsch-Gordan lumping). Calibration knob that
# rescales the *effective* |d|² in Maxwell-Bloch propagation; leaves the OBE
# solve (Γ, Ω via I_sat) unchanged. See README "Traps". Tune to match experiment.
LINE_STRENGTH_FACTOR = 1.0

# ---- Derived angular quantities ----
GAMMA = 2 * np.pi * GAMMA_2PI
OMEGA_HF = 2 * np.pi * NU_HF
OMEGA_EXCITED_HF = 2 * np.pi * NU_EXCITED_HF_D1
K_VEC = 2 * np.pi / WAVELENGTH_D1_85RB
OMEGA_D1 = 2 * np.pi * NU_D1_85RB

# ---- Phenomenological ground-coherence decay ----
GAMMA_GG_2PI = 100e3
GAMMA_GG = 2 * np.pi * GAMMA_GG_2PI


def rabi_freq(power, waist):
    """I = 2P/(π w₀²),  Ω = Γ √(I / 2 I_sat).  rad/s."""
    I = 2 * power / (np.pi * waist**2)
    return GAMMA * np.sqrt(I / (2 * I_SAT))
