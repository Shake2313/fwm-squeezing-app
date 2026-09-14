"""Reduced finite-seed periodic atomic noise with shared physical pump inputs."""

from .. import constants as c
from ..quantum.periodic import periodic_atomic_noise
from .inputs import ReducedPowerInputs
from .model import reduced_pump_system
from .normalization import reduced_dipoles, optical_carriers
from .seed_validity import seed_harmonic_hamiltonian


def reduced_periodic_noise(inputs, carrier_amplitudes, *, mean_order=4,
                           phase_samples=96, convention="uniform-zeeman-rms"):
    """Prescribed local coherent carriers; no propagation or Floquet field closure."""
    if not isinstance(inputs, ReducedPowerInputs):
        raise TypeError("ReducedPowerInputs required")
    dipoles = reduced_dipoles(convention)
    h0, reservoirs = reduced_pump_system(dipoles.pump_rabi_rad_s(inputs.pump_power_W, inputs.pump_waist_m),
        inputs.one_photon_rad_s, transit_rate_s_inverse=inputs.transit_rate_s_inverse,
        transition_scales=dipoles.transition_scales)
    v = seed_harmonic_hamiltonian(dipoles, optical_carriers(inputs.detunings)[1:],
                                 inputs.uniform_area_m2, carrier_amplitudes)
    return periodic_atomic_noise(h0, v, reservoirs, -c.OMEGA_HF+inputs.two_photon_rad_s,
                                 mean_order=mean_order, phase_samples=phase_samples)
