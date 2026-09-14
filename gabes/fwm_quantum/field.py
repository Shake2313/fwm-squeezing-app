"""Conditional reduced four-level field bridge in a shared uniform mode area.

The effective weak-field dipole is explicit. The adapter does not infer it from
an empirical squeezing efficiency, reinterpret arbitrary unequal Gaussian modes,
or certify the legacy pump-power/Rabi conversion as an absolute atom-light model.
"""

from dataclasses import dataclass, replace

import numpy as np

from .. import constants, observables
from ..quantum.contracts import AnalysisFrequencyAxis, GeneratorFrequencyAxis, OpticalDetunings, readonly_array
from ..quantum.traveling import LocalNambuGenerator, eliminate_atomic_noise
from ..schemes import fwm
from .model import minus_branch_atomic_frequencies
from .normalization import validated_transition_scales


@dataclass(frozen=True)
class ReducedFieldPair:
    analysis_axis: AnalysisFrequencyAxis
    main: LocalNambuGenerator
    companion: LocalNambuGenerator


def reduced_readout_operators(transition_scales=None):
    """O=(probe lowering on g2, conjugate raising on g1); no private reference import."""
    ops = np.zeros((2, 4, 4), complex)
    scales = (fwm.TRANSITION_DIPOLE_SCALE if transition_scales is None
              else validated_transition_scales(transition_scales))
    for excited in fwm.EXCITED_STATES:
        ops[0, fwm.G2, excited] = scales[fwm.G2, excited]
        ops[1, excited, fwm.G1] = scales[fwm.G1, excited]
    return ops


def reduced_field_pair(atom, detunings, analysis_axis, *, number_density_m3,
                       uniform_area_m2, optical_omega_rad_s, effective_dipole_C_m,
                       phase_mismatch_rad_m=0.0, transition_scales=None):
    """Compute both conjugate sectors independently, retaining both noise orderings.

    optical_omega_rad_s=(probe,conjugate) is explicit, not inferred from RF axes.
    The caller must match the atomic state's one-photon detuning to detunings.
    This same-area phase-selected model has no transverse/Doppler integration.
    """
    if not isinstance(detunings, OpticalDetunings) or not isinstance(analysis_axis, AnalysisFrequencyAxis):
        raise TypeError("separate optical detunings and RF axis required")
    density, area, dipole, mismatch = map(float, (
        number_density_m3, uniform_area_m2, effective_dipole_C_m, phase_mismatch_rad_m))
    if (not np.isfinite([density, area, dipole, mismatch]).all()
            or density < 0 or area <= 0 or dipole <= 0):
        raise ValueError("finite density>=0, area>0, dipole>0 and phase mismatch required")
    omega = readonly_array(optical_omega_rad_s, real=True)
    if omega.shape != (2,) or not np.isfinite(omega).all() or np.any(omega <= 0):
        raise ValueError("two positive optical angular frequencies required")
    q = np.diag(observables.photon_flux_mode_matrix(*omega, area, area))
    coupling = dipole*q/(2*constants.HBAR)
    main_axis = minus_branch_atomic_frequencies(detunings, analysis_axis)
    companion_axis = GeneratorFrequencyAxis(
        fwm.OMEGA_HF-detunings.two_photon_rad_s+analysis_axis.omega_rad_s,
        "rb85-d1-static-pump-companion-sector")
    ops = reduced_readout_operators(transition_scales)
    sectors = []
    for operators, signs, axis in (
            (ops, [1, -1], main_axis),
            (ops.conj().swapaxes(-1, -2), [-1, 1], companion_axis)):
        local = eliminate_atomic_noise(atom, operators, coupling, signs, axis,
                                       linear_density_m_inverse=density*area,
                                       mode_labels=("probe", "conjugate"))
        drift = local.drift-0.5j*mismatch*np.diag(signs)
        local = replace(local, drift=drift)
        if not local.audit()["passed"]:
            raise ValueError("reduced field sector failed local audit")
        sectors.append(local)
    return ReducedFieldPair(analysis_axis, *sectors)
