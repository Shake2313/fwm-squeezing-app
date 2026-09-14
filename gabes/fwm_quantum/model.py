"""Single-velocity reduced D1 pump model with explicitly declared reservoirs."""

import numpy as np

from .. import atoms
from ..quantum.contracts import AnalysisFrequencyAxis, GeneratorFrequencyAxis, OpticalDetunings
from ..quantum.diffusion import stationary_atomic_noise
from ..quantum.reservoirs import ExplicitReservoirs, thermal_reset_channels
from ..schemes import fwm
from .normalization import validated_transition_scales


def reduced_pump_system(pump_rabi_rad_s, effective_one_photon_detuning_rad_s, *,
                        transit_rate_s_inverse=0.0, transition_scales=None):
    """Strong classical pump, no seed/back-action; no implicit collision floor.

    Optional transit is the declared thermal replacement model, not a measured
    transport prediction. No residual coupling or noise-fit coefficient enters.
    """
    pump = float(pump_rabi_rad_s)
    detuning = float(effective_one_photon_detuning_rad_s)
    if not np.isfinite(pump) or pump < 0 or not np.isfinite(detuning):
        raise ValueError("finite nonnegative pump Rabi and finite detuning required")
    radiative = ExplicitReservoirs.from_atom(
        atoms.double_lambda_rb85(gamma_gg=0),
        source="GABES D1 radiative CF2 branching and natural decay")
    reset = thermal_reset_channels(
        transit_rate_s_inverse, [5/12, 7/12, 0, 0],
        source="explicitly supplied Markov thermal replacement rate")
    reservoirs = ExplicitReservoirs(4, radiative.channels+reset)
    h = fwm.pump_hamiltonian_at_Deff_zero(pump, pump)
    if transition_scales is not None:
        scales = validated_transition_scales(transition_scales)
        h[:2, 2:] = pump*scales[:2, 2:]/2
        h[2:, :2] = h[:2, 2:].conj().T
    h = h-detuning*np.diag([0., 0., 1., 1.])
    return h, reservoirs


def reduced_pump_noise(pump_rabi_rad_s, effective_one_photon_detuning_rad_s, *,
                       transit_rate_s_inverse=0.0, transition_scales=None):
    """Stationary microscopic noise from the shared explicit pump system."""
    h, reservoirs = reduced_pump_system(pump_rabi_rad_s, effective_one_photon_detuning_rad_s,
        transit_rate_s_inverse=transit_rate_s_inverse, transition_scales=transition_scales)
    return stationary_atomic_noise(h, reservoirs)


def minus_branch_atomic_frequencies(detunings, analysis_axis):
    """Map lab RF to -OMEGA_HF + delta + Omega in the existing static pump frame.

    This is an atomic sector frequency conversion, not a field-noise/readout
    transformation or a claim that its spectrum is the measured S_minus.
    """
    if not isinstance(detunings, OpticalDetunings) or not isinstance(analysis_axis, AnalysisFrequencyAxis):
        raise TypeError("separate OpticalDetunings and AnalysisFrequencyAxis required")
    return GeneratorFrequencyAxis(
        -fwm.OMEGA_HF+detunings.two_photon_rad_s+analysis_axis.omega_rad_s,
        "rb85-d1-static-pump-minus-sector")
