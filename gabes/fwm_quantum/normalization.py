"""Reciprocal reduced dipoles and explicitly anchored optical carriers.

The RMS reduction matches unpolarized, pump-off line strengths. It does not
preserve Zeeman coherences, optical pumping or nonlinear interference phases.
"""

from dataclasses import dataclass

import numpy as np

from .. import constants as c, hyperfine, species
from ..quantum.contracts import OpticalDetunings, readonly_array
from ..schemes import fwm


def gaussian_peak_field(power_W, waist_m):
    """Real peak E [V/m], for total power and circular intensity 1/e^2 radius."""
    power, waist = float(power_W), float(waist_m)
    if not np.isfinite([power, waist]).all() or power < 0 or waist <= 0:
        raise ValueError("finite power>=0 and waist>0 required")
    return float(np.sqrt(4*power/(np.pi*waist**2*c.EPS_0*c.C_LIGHT)))


def d1_dipole_from_decay(gamma_s_inverse, centroid_hz):
    """Steck Eq. 38 convention: Jg=Je=1/2, NOT the sqrt(2Jg+1) convention."""
    gamma, nu = float(gamma_s_inverse), float(centroid_hz)
    if not np.isfinite([gamma, nu]).all() or min(gamma, nu) <= 0:
        raise ValueError("positive finite decay and centroid required")
    return float(np.sqrt(3*np.pi*c.EPS_0*c.HBAR*c.C_LIGHT**3*gamma/(2*np.pi*nu)**3))


def validated_transition_scales(scales):
    values = readonly_array(scales, real=True)
    allowed = np.zeros((4, 4), bool)
    allowed[:2, 2:] = True
    if values.shape != (4, 4) or np.any(values < 0) or np.any(values[~allowed] != 0):
        raise ValueError("nonnegative real ground-to-excited 4x4 scales required")
    return values


@dataclass(frozen=True)
class ReducedDipoles:
    base_dipole_C_m: float
    transition_scales: np.ndarray
    convention: str
    scope: str

    def __post_init__(self):
        value = float(self.base_dipole_C_m)
        if not np.isfinite(value) or value <= 0 or not self.convention.strip() or not self.scope.strip():
            raise ValueError("positive dipole and explicit convention/scope required")
        object.__setattr__(self, "base_dipole_C_m", value)
        object.__setattr__(self, "transition_scales", validated_transition_scales(self.transition_scales))

    @property
    def matrix_C_m(self):
        return readonly_array(self.base_dipole_C_m*self.transition_scales, real=True)

    def pump_rabi_rad_s(self, power_W, waist_m):
        return self.base_dipole_C_m*gaussian_peak_field(power_W, waist_m)/c.HBAR


def reduced_dipoles(convention):
    """Named, fixed reductions; no adjustable efficiency or noise scale."""
    if convention == "legacy-reciprocal":
        return ReducedDipoles(c.DIPOLE_D1/np.sqrt(12), fwm.TRANSITION_DIPOLE_SCALE,
            convention, "historical weak-field dipoles reused for pump; not a Zeeman reduction")
    if convention == "uniform-zeeman-rms":
        scales = np.zeros((4, 4))
        for g, Fg in fwm.GROUND_F.items():
            for e, Fe in fwm.EXCITED_F.items():
                scales[g, e] = np.sqrt(species.line_strength(Fg, Fe, 2.5, .5, .5)/3)
        return ReducedDipoles(d1_dipole_from_decay(c.GAMMA, c.NU_D1_85RB), scales,
            convention, "uniform sublevels within each F, weak absorption RMS; positive four-level amplitudes are a conditional nonlinear surrogate")
    raise ValueError("unknown reduced dipole convention")


def optical_carriers(detunings):
    """[pump,probe,conjugate] rad/s; Delta references F=2 -> F'=3.

    Use exactly the splittings in the static Hamiltonian, including its current
    excited splitting rounding, instead of mixing in species.A_P12*3.
    """
    if not isinstance(detunings, OpticalDetunings):
        raise TypeError("OpticalDetunings required")
    # Ground F2=-7/12*split_g, excited F3=+5/12*split_e.
    anchor = c.OMEGA_D1+(7*c.OMEGA_HF+5*c.OMEGA_EXCITED_HF)/12
    pump = anchor+detunings.one_photon_rad_s
    beat = -c.OMEGA_HF+detunings.two_photon_rad_s
    values = readonly_array([pump, pump+beat, pump-beat], real=True)
    if np.any(values <= 0):
        raise ValueError("positive optical carriers required")
    return values


def normalization_audit():
    """Algebraic discrepancy ledger; successful audit does not approve old inputs."""
    legacy_d = c.HBAR*c.GAMMA*np.sqrt(c.EPS_0*c.C_LIGHT/(4*c.I_SAT))
    weak = reduced_dipoles("legacy-reciprocal")
    dipole = d1_dipole_from_decay(c.GAMMA, c.NU_D1_85RB)
    rows = []
    for (Fg, Fe), cf2 in hyperfine.CF2.items():
        S = species.line_strength(Fg, Fe, 2.5, .5, .5)
        rows.append({"Fg": Fg, "Fe": Fe, "CF2": cf2, "S_FFprime": S,
            "single_polarization_sum_over_m_in_dJ2": 3*cf2,
            "uniform_manifold_mean_in_dJ2": S/3,
            "historical_weak_strength_in_dJ2": 3*cf2/12,
            "historical_to_uniform_strength_ratio": (3*cf2/12)/(S/3),
            "required_scalar_after_summed_strength": 1/(2*Fg+1),
            "sum_rule_error": abs(3*cf2-(2*Fg+1)*S/3)})
    return {"legacy_pump_implied_dipole_C_m": legacy_d,
        "historical_weak_base_dipole_C_m": weak.base_dipole_C_m,
        "pump_to_weak_dipole_ratio": legacy_d/weak.base_dipole_C_m,
        "legacy_pump_to_far_detuned_dipole_ratio": legacy_d/(c.DIPOLE_D1/np.sqrt(3)),
        "decay_consistent_dJ_C_m": dipole,
        "stored_dipole_implied_decay_relative_error": (c.DIPOLE_D1/dipole)**2-1,
        "F2_Fprime3_anchor_minus_centroid_hz": (7*c.NU_GROUND_HF+5*c.NU_EXCITED_HF_D1)/12,
        "species_minus_hamiltonian_excited_splitting_hz": 3*species.RB85.A_P12*1e6-c.NU_EXCITED_HF_D1,
        "transition_rows": rows, "single_scalar_matches_both_manifolds": False,
        "legacy_pump_and_weak_dipoles_match": False,
        "full_zeeman_reduction_proven": False}
