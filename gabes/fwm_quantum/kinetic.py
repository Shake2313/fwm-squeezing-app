"""Velocity-resolved reduced fields with an explicit carrier-momentum closure.

The single-phase Floquet construction requires q_probe=-q_conjugate. General
vacuum noncollinear beams fail this condition and need spatial grating dynamics;
their mismatch is not erased, fitted, or substituted into a lab optical carrier.
"""

from dataclasses import dataclass
from math import erf, sqrt

import numpy as np

from .. import constants as c, doppler, observables
from ..quantum.contracts import AnalysisFrequencyAxis, GeneratorFrequencyAxis, readonly_array
from ..quantum.diffusion import stationary_atomic_noise
from ..quantum.periodic import periodic_atomic_noise
from ..quantum.periodic_field import PeriodicFieldPorts, eliminate_periodic_atom
from ..quantum.traveling import LocalNambuGenerator, eliminate_atomic_noise, sum_independent_classes
from .field import reduced_readout_operators
from .inputs import ReducedPowerInputs
from .model import reduced_pump_system
from .normalization import optical_carriers, reduced_dipoles


@dataclass(frozen=True)
class CarrierGeometry:
    """Wavevectors (pump,probe,conjugate), in rad/m in Cartesian (x,y,z).

    Explicit directions use a common z-plane photon flux and overlap area A_z.
    Non-vacuum wavevector magnitudes must be declared; no refractive index is
    inferred, and no mismatch cancellation is performed by this object.
    """
    wavevectors_rad_m: np.ndarray
    source: str

    def __post_init__(self):
        k = readonly_array(self.wavevectors_rad_m, real=True)
        if (k.shape != (3, 3) or np.any(k[:, 2] <= 0)
                or np.linalg.norm(k[0, :2]) > 1e-12*np.linalg.norm(k[0])
                or not isinstance(self.source, str) or not self.source.strip()):
            raise ValueError('three forward wavevectors, pump along +z and provenance required')
        object.__setattr__(self, 'wavevectors_rad_m', k)

    @property
    def mismatch_rad_m(self):
        return readonly_array(2*self.wavevectors_rad_m[0]-self.wavevectors_rad_m[1]-self.wavevectors_rad_m[2], real=True)

    def closure_audit(self, *, temperature_K=None, length_m=None):
        k = self.wavevectors_rad_m
        mismatch = self.mismatch_rad_m
        tolerance = 1e-12*float(np.max(np.linalg.norm(k, axis=1)))
        result = {'passed': bool(np.linalg.norm(mismatch) <= tolerance),
            'mismatch_rad_m': mismatch.tolist(), 'tolerance_rad_m': tolerance,
            'scope': 'single Raman space-time phase only, not Maxwell dispersion or experimental phase matching'}
        if temperature_K is not None:
            temperature = float(temperature_K)
            if not np.isfinite(temperature) or temperature <= 0:
                raise ValueError('positive finite kinetic temperature required')
            result['unresolved_loop_doppler_rms_rad_s'] = float(np.linalg.norm(mismatch)*np.sqrt(c.KB*temperature/c.MASS_85RB))
        if length_m is not None:
            length = float(length_m)
            if not np.isfinite(length) or length < 0:
                raise ValueError('finite nonnegative cell length required')
            result['axial_mismatch_phase_rad'] = float(mismatch[2]*length)
        return result

    def require_closed(self):
        if not self.closure_audit()['passed']:
            raise ValueError(f'single-phase Floquet momentum closure failed: {self.closure_audit()}')

    def atomic_frequencies(self, inputs, velocities_m_s):
        """Return Delta_v, signed nu_probe_v, nu_conjugate_v; preserve lab inputs."""
        velocity = readonly_array(velocities_m_s, real=True)
        if velocity.ndim != 2 or velocity.shape[1] != 3:
            raise ValueError('Cartesian velocity rows required')
        k0, kp, kc = self.wavevectors_rad_m
        beat = -c.OMEGA_HF+inputs.two_photon_rad_s
        return (inputs.one_photon_rad_s-velocity@k0,
                beat-velocity@(kp-k0), -beat-velocity@(kc-k0))

    @classmethod
    def vacuum_beams(cls, inputs, *, probe_angle_rad=0., conjugate_angle_rad=0.):
        theta = readonly_array([0., probe_angle_rad, conjugate_angle_rad], real=True)
        omega = optical_carriers(inputs.detunings)
        k = omega[:, None]/c.C_LIGHT*np.column_stack([np.sin(theta), np.zeros(3), np.cos(theta)])
        return cls(k, 'declared vacuum plane waves; opposite equal angles do not imply vector closure')


@dataclass(frozen=True)
class VelocityQuadrature:
    velocities_m_s: np.ndarray
    probabilities: np.ndarray
    source: str
    omitted_probability: float = 0.
    integrated_velocity_axes: tuple[int, ...] = ()

    def __post_init__(self):
        v = readonly_array(self.velocities_m_s, real=True)
        p = readonly_array(self.probabilities, real=True)
        tail = float(self.omitted_probability)
        axes = tuple(self.integrated_velocity_axes)
        if (v.ndim != 2 or v.shape[1] != 3 or not len(v) or p.shape != (len(v),)
                or np.any(p <= 0) or abs(p.sum()-1) > 1e-12 or not np.isfinite(tail) or not 0 <= tail < 1
                or not isinstance(self.source, str) or not self.source.strip()
                or len(set(axes)) != len(axes) or any(a not in (0, 1, 2) for a in axes)):
            raise ValueError('finite velocities, positive probabilities summing to one, tail and provenance required')
        object.__setattr__(self, 'velocities_m_s', v)
        object.__setattr__(self, 'probabilities', p)
        object.__setattr__(self, 'omitted_probability', tail)
        object.__setattr__(self, 'integrated_velocity_axes', axes)

    @classmethod
    def maxwell_xz(cls, temperature_K, *, longitudinal_order, transverse_order, cutoff_sigma=5.):
        for order in (longitudinal_order, transverse_order):
            if isinstance(order, bool) or int(order) != order or order < 2:
                raise ValueError('integer quadrature orders >=2 required')
        vx, wx = doppler.maxwell_legendre_grid(temperature_K, order=transverse_order, cutoff_sigma=cutoff_sigma)
        vz, wz = doppler.maxwell_legendre_grid(temperature_K, order=longitudinal_order, cutoff_sigma=cutoff_sigma)
        xx, zz = np.meshgrid(vx, vz, indexing='xy')
        return cls(np.column_stack([xx.ravel(), np.zeros(xx.size), zz.ravel()]), (wz[:, None]*wx[None, :]).ravel(),
            'independent Maxwell vx,vz; vy analytically integrated for coplanar x-z waves; retained square normalized to one',
            1-erf(float(cutoff_sigma)/sqrt(2))**2, (1,))


def kinetic_ports(inputs, geometry, *, convention='uniform-zeeman-rms'):
    """Same z-plane canonical photon flux for field driving and readout."""
    if not isinstance(inputs, ReducedPowerInputs) or not isinstance(geometry, CarrierGeometry):
        raise TypeError('physical inputs and declared carrier geometry required')
    if inputs.phase_mismatch_rad_m != 0:
        raise ValueError('do not combine an independent scalar mismatch with explicit wavevectors')
    geometry.require_closed()
    dipoles = reduced_dipoles(convention)
    omega = optical_carriers(inputs.detunings)[1:]
    vectors = geometry.wavevectors_rad_m[1:]
    cosine = vectors[:, 2]/np.linalg.norm(vectors, axis=1)
    q = np.diag(observables.photon_flux_mode_matrix(*omega, inputs.uniform_area_m2, inputs.uniform_area_m2))
    g = dipoles.base_dipole_C_m*q/(2*c.HBAR*np.sqrt(cosine))
    ops = reduced_readout_operators(dipoles.transition_scales)
    ops[1] = ops[1].conj().T
    return PeriodicFieldPorts(('probe', 'conjugate'), [1, -1], ops, g)


def _pump_system(inputs, delta_v, convention):
    dipoles = reduced_dipoles(convention)
    return reduced_pump_system(dipoles.pump_rabi_rad_s(inputs.pump_power_W, inputs.pump_waist_m), delta_v,
        transit_rate_s_inverse=inputs.transit_rate_s_inverse, transition_scales=dipoles.transition_scales)


def periodic_velocity_atom(inputs, geometry, velocity_m_s, carrier_amplitudes, *, mean_order=4,
                            convention='uniform-zeeman-rms'):
    """Local ballistic/slow-envelope approximation; no velocity-changing reservoir."""
    ports = kinetic_ports(inputs, geometry, convention=convention)
    velocity = readonly_array(velocity_m_s, real=True)
    beta = readonly_array(carrier_amplitudes)
    if velocity.shape != (3,) or beta.shape != (2,):
        raise ValueError('one velocity and two laboratory-normalized carrier amplitudes required')
    delta, beat, _ = geometry.atomic_frequencies(inputs, velocity[None])
    h0, reservoirs = _pump_system(inputs, delta[0], convention)
    op, g = ports.lowering_operators, ports.coupling_s_inverse_sqrt_flux
    v = g[0]*beta[0]*op[0].conj().T+g[1]*beta[1].conjugate()*op[1]
    return periodic_atomic_noise(h0, v, reservoirs, beat[0], mean_order=mean_order)


def _stationary_full_field(atom, ports, beat_v, axis, density):
    """Exactly the zero-seed periodic limit, without a redundant harmonic lift."""
    hs, ops, g, signs, labels = ports.nambu_coordinates()
    shape = (len(axis.omega_rad_s), 4, 4)
    drift = np.zeros(shape, complex)
    greater = np.zeros((len(atom.reservoir_names),)+shape, complex)
    lesser = np.zeros_like(greater)
    for harmonic, ids in ((1, [0, 3]), (-1, [2, 1])):
        sector = eliminate_atomic_noise(atom, ops[ids], g[ids], signs[ids],
            GeneratorFrequencyAxis(axis.omega_rad_s+harmonic*beat_v, 'moving-atom-pump-frame-sector'),
            linear_density_m_inverse=density, mode_labels=tuple(labels[j] for j in ids))
        for i, row in enumerate(ids):
            for j, col in enumerate(ids):
                drift[:, row, col] = sector.drift[:, i, j]
                greater[:, :, row, col] = sector.noise_greater_by_reservoir[:, :, i, j]
                lesser[:, :, row, col] = sector.noise_lesser_by_reservoir[:, :, i, j]
    return LocalNambuGenerator(axis, signs, labels, drift, greater, lesser, atom.reservoir_names)


def kinetic_local_field(inputs, geometry, quadrature, analysis_axis, *, carrier_amplitudes=None,
                         mean_order=4, response_order=3, convention='uniform-zeeman-rms', retain_classes=True):
    """Sum local responses and independent covariance sources at the same lab RF.

    None selects the exact pump-state/weak-field limit. Finite carriers select
    velocity-resolved periodic means and full Nambu noise. No power or field
    transmission is averaged across velocities after cell propagation.
    """
    if not isinstance(quadrature, VelocityQuadrature) or not isinstance(analysis_axis, AnalysisFrequencyAxis):
        raise TypeError('velocity quadrature and laboratory RF axis required')
    ports = kinetic_ports(inputs, geometry, convention=convention)
    if quadrature.integrated_velocity_axes and np.any(abs(geometry.wavevectors_rad_m[:, quadrature.integrated_velocity_axes]) > 1e-12):
        raise ValueError('wavevectors depend on a velocity axis that this quadrature integrated out')
    axis = GeneratorFrequencyAxis(analysis_axis.omega_rad_s, 'laboratory-RF-in-common-closed-carrier-frame')
    delta, beats, conjugate_beats = geometry.atomic_frequencies(inputs, quadrature.velocities_m_s)
    classes, states, polarizations = [], [], []
    accumulated = None
    maximum_class_residual = 0.
    pump_cache = {}
    for velocity, probability, dv, beat in zip(quadrature.velocities_m_s, quadrature.probabilities, delta, beats):
        density = inputs.number_density_m3*inputs.uniform_area_m2*probability
        if carrier_amplitudes is None:
            if float(dv) not in pump_cache:
                pump_cache[float(dv)] = stationary_atomic_noise(*_pump_system(inputs, dv, convention))
            atom = pump_cache[float(dv)]
            local = _stationary_full_field(atom, ports, beat, axis, density)
            states.append(atom.stationary_state)
            polarizations.append(np.zeros(2, complex))
        else:
            atom = periodic_velocity_atom(inputs, geometry, velocity, carrier_amplitudes,
                mean_order=mean_order, convention=convention)
            local = eliminate_periodic_atom(atom, ports, axis, response_order=response_order, linear_density_m_inverse=density)
            states.append(atom.state_harmonics[mean_order])
            polarizations.append(np.array([np.trace(op@atom.state_harmonics[mean_order+int(h)])
                for op, h in zip(ports.lowering_operators, ports.carrier_harmonics)]))
        audit = local.audit()
        if not audit['passed']:
            raise ValueError('one velocity class failed canonical commutator/PSD')
        maximum_class_residual = max(maximum_class_residual, audit['maximum_commutator_relative_residual'])
        if retain_classes:
            classes.append(local)
        else:
            if accumulated is None:
                accumulated = [np.zeros_like(x) for x in (local.drift, local.noise_greater_by_reservoir, local.noise_lesser_by_reservoir)]
                reservoir_names = local.reservoir_names
            if reservoir_names != local.reservoir_names:
                raise ValueError('streaming velocity sum requires matched physical reservoir names')
            for total, value in zip(accumulated, (local.drift, local.noise_greater_by_reservoir, local.noise_lesser_by_reservoir)):
                total += value
    combined = (sum_independent_classes(classes) if retain_classes else
        LocalNambuGenerator(axis, local.signs, local.mode_labels, *accumulated,
                            tuple('independent-velocity-sum:'+name for name in reservoir_names)))
    if not combined.audit()['passed']:
        raise ValueError('velocity mixture failed canonical commutator/PSD')
    return {'generator': combined, 'ports': ports, 'class_generators': tuple(classes),
        'maximum_class_commutator_relative_residual': maximum_class_residual,
        'mean_states': readonly_array(states), 'class_polarizations': readonly_array(polarizations),
        'mean_polarization': readonly_array(quadrature.probabilities@np.array(polarizations)),
        'atomic_one_photon_rad_s': readonly_array(delta, real=True),
        'atomic_signed_beats_rad_s': readonly_array(np.column_stack([beats, conjugate_beats]), real=True),
        'lab_optical_omega_rad_s': readonly_array(optical_carriers(inputs.detunings), real=True),
        'lab_signed_beat_rad_s': -c.OMEGA_HF+inputs.two_photon_rad_s,
        'finite_seed': carrier_amplitudes is not None,
        'scope': 'local independent ballistic velocity classes, closed carrier momentum, slow spatial envelope; no class-changing collisions'}


def thermal_pump_cell(inputs, geometry, quadrature, detector, *, convention='uniform-zeeman-rms'):
    """Conditional uniform pump-state thermal gain and bright S_minus.

    The seeded coherent mean is linear in the input because the reference atom
    is pump-only. This does not replace the finite-seed nonlinear kinetic path.
    """
    from ..quantum.periodic_field import paired_frequency_channel
    from ..quantum.readout import intensity_difference_spectrum
    from ..quantum.sidebands import SIDEBAND_MODES
    from ..quantum.traveling import constant_segment
    rf = detector.analysis_axis.omega_rad_s
    axis = AnalysisFrequencyAxis(np.unique(np.r_[-rf, 0., rf]))
    local = kinetic_local_field(inputs, geometry, quadrature, axis, convention=convention, retain_classes=False)
    transfer = constant_segment(local['generator'], inputs.length_m)
    omega = optical_carriers(inputs.detunings)[1:]
    beta_in = np.array([np.sqrt(inputs.seed_power_W/(c.HBAR*omega[0]))*np.exp(1j*inputs.seed_phase_rad), 0.])
    index = np.flatnonzero(axis.omega_rad_s == 0)[0]
    beta = (transfer.transfer[index]@np.r_[beta_in, beta_in.conj()])[:2]
    perm = [0, 1, 6, 7, 4, 5, 2, 3]
    covariances, channel_audits = [], []
    for w in rf:
        index, mirror = [np.flatnonzero(axis.omega_rad_s == sign*w)[0] for sign in (1, -1)]
        channel = paired_frequency_channel(transfer, ('probe', 'conjugate'), index, mirror)
        covariance = channel.apply_covariance(np.eye(8)/2)
        covariances.append(covariance[perm][:, perm])
        channel_audits.append(channel.audit())
    spectrum = intensity_difference_spectrum(covariances, beta, detector,
        mode_labels=SIDEBAND_MODES, analysis_axis=detector.analysis_axis)
    powers = c.HBAR*omega*abs(beta)**2
    return {'local': local, 'transfer': transfer, 'spectrum': spectrum,
        'output_carriers_sqrt_flux': readonly_array(beta), 'output_powers_W': readonly_array(powers, real=True),
        'probe_power_gain': float(powers[0]/inputs.seed_power_W),
        'conjugate_power_gain': float(powers[1]/inputs.seed_power_W),
        'minimum_channel_cp_eigenvalue': min(a.minimum_cp_eigenvalue for a in channel_audits),
        'scope': 'uniform classical pump, thermal pump-state weak-field mean, bright detector, closed carrier wavevectors'}
