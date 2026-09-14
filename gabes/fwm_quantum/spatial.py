"""Two-phase local Rb atom for explicit, possibly nonclosed wavevectors."""

from dataclasses import dataclass

import numpy as np

from .. import constants as c, observables
from ..quantum.contracts import readonly_array
from ..quantum.spatial import rectangle, torus_atomic_noise
from .field import reduced_readout_operators
from .inputs import ReducedPowerInputs
from .kinetic import CarrierGeometry, _pump_system, periodic_velocity_atom
from .normalization import optical_carriers, reduced_dipoles


@dataclass(frozen=True)
class SpatialVelocityAtom:
    atom: object
    wavevectors_rad_m: np.ndarray
    lab_frequencies_rad_s: np.ndarray
    velocity_m_s: np.ndarray
    convective_frequencies_rad_s: np.ndarray
    drive_operators_rad_s: np.ndarray
    optical_lowering_operators: np.ndarray
    coupling_s_inverse_sqrt_flux: np.ndarray
    coordinate_kind: str

    def __post_init__(self):
        for key in ('wavevectors_rad_m', 'lab_frequencies_rad_s', 'velocity_m_s',
                    'convective_frequencies_rad_s', 'drive_operators_rad_s',
                    'optical_lowering_operators', 'coupling_s_inverse_sqrt_flux'):
            object.__setattr__(self, key, readonly_array(getattr(self, key),
                real=key not in ('drive_operators_rad_s', 'optical_lowering_operators')))

    def at_events(self, times_s, positions_m):
        times, positions = readonly_array(times_s, real=True), readonly_array(positions_m, real=True)
        if times.ndim != 1 or positions.shape != (len(times), 3) or not len(times):
            raise ValueError('matched laboratory times and Cartesian positions required')
        phases = times[:, None]*self.lab_frequencies_rad_s-positions@self.wavevectors_rad_m.T
        if self.coordinate_kind == 'closed-single-phase-quotient':
            return self.atom.at_phase(phases[:, 0])
        return self.atom.at_phase(phases)


def spatial_velocity_atom(inputs, geometry, velocity_m_s, carrier_amplitudes, *,
                          mean_orders=(3, 3), convention='uniform-zeeman-rms'):
    """Keep q=kp-k0 and Q=kp+kc-2k0 with drive labels (1,0),(-1,1).

    Fixed plane-wave envelopes, independent ballistic atom and constant pump.
    This does not define spatially collected field noise or thermal propagation.
    """
    if not isinstance(inputs, ReducedPowerInputs) or not isinstance(geometry, CarrierGeometry):
        raise TypeError('declared power inputs and carrier geometry required')
    if inputs.phase_mismatch_rad_m != 0:
        raise ValueError('explicit wavevectors cannot be combined with a scalar mismatch')
    rectangle(mean_orders)  # Validate both cutoffs, including the quotient path.
    velocity, beta = readonly_array(velocity_m_s, real=True), readonly_array(carrier_amplitudes)
    if velocity.shape != (3,) or beta.shape != (2,):
        raise ValueError('one velocity and two canonical laboratory carriers required')
    k0, kp, kc = geometry.wavevectors_rad_m
    phase_vectors = np.array([kp-k0, kp+kc-2*k0])
    lab = np.array([-c.OMEGA_HF+inputs.two_photon_rad_s, 0.])
    frequencies = lab-phase_vectors@velocity
    dipoles = reduced_dipoles(convention)
    omega = optical_carriers(inputs.detunings)[1:]
    cosine = geometry.wavevectors_rad_m[1:, 2]/np.linalg.norm(geometry.wavevectors_rad_m[1:], axis=1)
    flux = np.diag(observables.photon_flux_mode_matrix(*omega, inputs.uniform_area_m2, inputs.uniform_area_m2))
    g = dipoles.base_dipole_C_m*flux/(2*c.HBAR*np.sqrt(cosine))
    ops = reduced_readout_operators(dipoles.transition_scales)
    ops[1] = ops[1].conj().T
    drives = g[:, None, None]*beta[:, None, None]*ops.conj().swapaxes(-1, -2)
    if geometry.closure_audit()['passed']:
        # Q=0 identifies every (h,l) with h at the physical phase theta_2=0.
        # EXACT quotient: combine V_(1,0)+V_(-1,1)^dagger before solving.
        # Averaging independent l copies would erase coherent cross terms.
        atom = periodic_velocity_atom(inputs, geometry, velocity, beta,
            mean_order=mean_orders[0], convention=convention)
        kind = 'closed-single-phase-quotient'
    else:
        h0, reservoirs = _pump_system(inputs, inputs.one_photon_rad_s-k0@velocity, convention)
        atom = torus_atomic_noise(h0, [(1, 0), (-1, 1)], drives, reservoirs, frequencies, mean_orders=mean_orders)
        kind = 'two-phase-convective-lattice'
    return SpatialVelocityAtom(atom, phase_vectors, lab, velocity, frequencies, drives, ops, g, kind)
