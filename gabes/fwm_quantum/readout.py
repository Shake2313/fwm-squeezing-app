"""Conditional seeded reduced-cell gain and bright intensity-difference readout.

Uses a single declared atomic model for coherent transfer and sideband noise.
No experimental badge, power/Rabi calibration, Doppler or depletion is supplied.
"""

from dataclasses import dataclass

import numpy as np

from ..constants import HBAR
from ..quantum.contracts import readonly_array
from ..quantum.readout import (
    IntensityDifferenceSpectrum, coherent_carrier_output, intensity_difference_spectrum,
)
from ..quantum.sidebands import SIDEBAND_MODES, SidebandChannels, sideband_channels
from ..quantum.traveling import NambuTransfer, constant_segment
from .field import reduced_field_pair


@dataclass(frozen=True)
class ReducedSeededReadout:
    main_transfer: NambuTransfer
    companion_transfer: NambuTransfer
    sidebands: SidebandChannels
    output_carrier_amplitudes_sqrt_flux: np.ndarray
    output_coherent_powers_W: np.ndarray
    probe_power_gain: float
    conjugate_photon_flux_gain: float
    conjugate_power_gain: float
    spectrum: IntensityDifferenceSpectrum

    def __post_init__(self):
        object.__setattr__(self, "output_carrier_amplitudes_sqrt_flux",
                           readonly_array(self.output_carrier_amplitudes_sqrt_flux))
        object.__setattr__(self, "output_coherent_powers_W", readonly_array(self.output_coherent_powers_W, real=True))


def reduced_seeded_readout(atom, detunings, analysis_axis, detector, *, number_density_m3,
                          uniform_area_m2, optical_omega_rad_s, effective_dipole_C_m,
                          length_m, seed_power_W, seed_phase_rad=0., phase_mismatch_rad_m=0.,
                          transition_scales=None):
    """Joint coherent gains and conditional linearized S_minus, with one seeded input.

    The RF axis must contain zero and both signs of every detector frequency.
    Mean powers omit spontaneous photons; bright-carrier dominance and weak
    atomic back-action must be checked before interpreting this as experiment.
    """
    power, phase = float(seed_power_W), float(seed_phase_rad)
    if not np.isfinite([power, phase]).all() or power <= 0:
        raise ValueError("positive seed power and finite seed phase required")
    omega = readonly_array(optical_omega_rad_s, real=True)
    pair = reduced_field_pair(atom, detunings, analysis_axis, number_density_m3=number_density_m3,
                               uniform_area_m2=uniform_area_m2, optical_omega_rad_s=omega,
                               effective_dipole_C_m=effective_dipole_C_m,
                               phase_mismatch_rad_m=phase_mismatch_rad_m,
                               transition_scales=transition_scales)
    main, companion = constant_segment(pair.main, length_m), constant_segment(pair.companion, length_m)
    seed_flux = power/(HBAR*omega[0])
    beta = coherent_carrier_output(main, analysis_axis, [np.sqrt(seed_flux)*np.exp(1j*phase), 0.])
    state = sideband_channels(main, companion, analysis_axis)
    result = intensity_difference_spectrum(state.vacuum_covariances(), beta, detector,
                                           mode_labels=SIDEBAND_MODES, analysis_axis=state.analysis_axis)
    flux = np.abs(beta)**2
    powers = HBAR*omega*flux
    return ReducedSeededReadout(main, companion, state, beta, powers,
                                float(powers[0]/power), float(flux[1]/seed_flux),
                                float(powers[1]/power), result)
