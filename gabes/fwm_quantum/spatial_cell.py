"""Finite-aperture two-band field for stationary atomic centers and fixed pump.

Nonclosed carrier geometry enters the same microscopic mean, M and both D
orderings. This is not a moving hot-vapor or complete transverse-mode model.
"""

import numpy as np
from scipy.integrate import solve_ivp

from .. import core, constants as c
from ..quantum.contracts import GeneratorFrequencyAxis, readonly_array
from ..quantum.periodic import periodic_atomic_noise
from ..quantum.periodic_field import PeriodicFieldPorts, eliminate_periodic_atom, paired_frequency_channel, adaptive_field_propagation
from ..quantum.readout import intensity_difference_spectrum
from ..quantum.sidebands import SIDEBAND_MODES
from ..quantum.spatial_modes import TransverseModeGrid, project_local_generator
from ..quantum.traveling import LocalNambuGenerator
from .inputs import ReducedPowerInputs
from .kinetic import CarrierGeometry
from .model import reduced_pump_system
from .normalization import optical_carriers, reduced_dipoles
from .periodic_cell import reduced_periodic_ports


class StationarySpatialMedium:
    """Galerkin projection onto two specified transverse optical profiles.

    Uniform number density and pump over the aperture. No atomic flux through
    transverse/longitudinal boundaries, diffraction or unretained optical ports.
    """

    def __init__(self, inputs, geometry, modes, *, mean_order=4, response_order=3,
                 velocity_m_s=(0., 0., 0.), convention='uniform-zeeman-rms'):
        if not isinstance(inputs, ReducedPowerInputs) or not isinstance(geometry, CarrierGeometry) or not isinstance(modes, TransverseModeGrid):
            raise TypeError('physical inputs, explicit geometry and normalized optical profiles required')
        if inputs.phase_mismatch_rad_m != 0:
            raise ValueError('do not combine explicit wavevectors with a separate scalar mismatch')
        velocity = readonly_array(velocity_m_s, real=True)
        if velocity.shape != (3,) or np.any(velocity != 0):
            raise ValueError('stationary atomic centers required: moving atoms need aperture/inter-slice transport noise')
        for value in (mean_order, response_order):
            if isinstance(value, bool) or int(value) != value or value < 1:
                raise ValueError('positive integer atomic harmonic cutoffs required')
        self.inputs, self.geometry, self.modes = inputs, geometry, modes
        self.mean_order, self.response_order = int(mean_order), int(response_order)
        base = reduced_periodic_ports(inputs, convention=convention)
        vectors = geometry.wavevectors_rad_m[1:]
        cosine = vectors[:, 2]/np.linalg.norm(vectors, axis=1)
        self.ports = PeriodicFieldPorts(base.mode_labels, base.carrier_harmonics,
            base.lowering_operators, base.coupling_s_inverse_sqrt_flux/np.sqrt(cosine))
        dipoles = reduced_dipoles(convention)
        self.h0, self.reservoirs = reduced_pump_system(
            dipoles.pump_rabi_rad_s(inputs.pump_power_W, inputs.pump_waist_m), inputs.one_photon_rad_s,
            transit_rate_s_inverse=inputs.transit_rate_s_inverse, transition_scales=dipoles.transition_scales)
        self.l0 = self.reservoirs.generator(self.h0)
        self.beat_rad_s = -c.OMEGA_HF+inputs.two_photon_rad_s

    def sections(self, z_m):
        z = float(z_m)
        if not np.isfinite(z):
            raise ValueError('finite longitudinal coordinate required')
        q_loop = -self.geometry.mismatch_rad_m
        psi = self.modes.positions_m@q_loop[:2]+z*q_loop[2]
        eta = np.sqrt(self.inputs.uniform_area_m2)*self.modes.mode_values_m_inverse.copy()
        eta[:, 1] *= np.exp(1j*psi)
        # EXACT time-origin quotient removes exp(i*h*q_probe.r) from atom and
        # readout together. Only exp(i*(q_probe+q_conjugate).r) remains here.
        # It eliminates rapid carrier-phase solves without erasing mismatch.
        # Identical eta rows have identical H/B/C/D: sum their atomic counts
        # first. This grouping changes neither the equations nor quadrature.
        unique, inverse = np.unique(eta, axis=0, return_inverse=True)
        areas = np.bincount(inverse, weights=self.modes.area_weights_m2)
        return tuple(zip(unique, areas))

    def harmonic_drive(self, beta, eta):
        amplitudes = readonly_array(beta)
        if amplitudes.shape != (2,):
            raise ValueError('two finite canonical carrier amplitudes required')
        g, ops = self.ports.coupling_s_inverse_sqrt_flux, self.ports.lowering_operators
        local = amplitudes*eta
        return g[0]*local[0]*ops[0].conj().T+g[1]*local[1].conjugate()*ops[1]

    def mean_rate(self, z_m, beta):
        beta = readonly_array(beta)
        if beta.shape != (2,):
            raise ValueError('two finite canonical carrier amplitudes required')
        rate = np.zeros(2, complex)
        # With zero density every atomic source vanishes exactly; no state
        # solve can change the vacuum mean or noise. Keep this exact limit.
        if self.inputs.number_density_m3 == 0:
            return rate
        for eta, area in self.sections(z_m):
            v = self.harmonic_drive(beta, eta)
            states = core.floquet_solve_truncated(self.l0, core.comm_super(v), core.comm_super(v.conj().T),
                self.beat_rad_s, [0.], np.zeros_like(self.l0), 4, n_f=self.mean_order, return_harmonics=True)[0]
            polarization = np.array([np.trace(op@states[self.mean_order+int(h)])
                for op, h in zip(self.ports.lowering_operators, self.ports.carrier_harmonics)])
            rate += -1j*self.inputs.number_density_m3*area*self.ports.coupling_s_inverse_sqrt_flux*eta.conj()*polarization
        return rate

    def local_field(self, z_m, beta, axis):
        if not isinstance(axis, GeneratorFrequencyAxis):
            raise TypeError('common laboratory RF generator axis required')
        beta = readonly_array(beta)
        if beta.shape != (2,):
            raise ValueError('two finite canonical carrier amplitudes required')
        if self.inputs.number_density_m3 == 0:
            names = tuple('independent-area-sum:'+channel.name for channel in self.reservoirs.channels)
            shape = (len(axis.omega_rad_s), 4, 4)
            *_, signs, labels = self.ports.nambu_coordinates()
            return LocalNambuGenerator(axis, signs, labels, np.zeros(shape, complex),
                np.zeros((len(names),)+shape, complex), np.zeros((len(names),)+shape, complex), names)
        accumulated = None
        for eta, area in self.sections(z_m):
            atom = periodic_atomic_noise(self.h0, self.harmonic_drive(beta, eta), self.reservoirs,
                self.beat_rad_s, mean_order=self.mean_order)
            raw = eliminate_periodic_atom(atom, self.ports, axis, response_order=self.response_order,
                linear_density_m_inverse=self.inputs.number_density_m3*area)
            local = project_local_generator(raw, eta)
            if accumulated is None:
                accumulated = [np.zeros_like(x) for x in (local.drift, local.noise_greater_by_reservoir, local.noise_lesser_by_reservoir)]
            for total, value in zip(accumulated, (local.drift, local.noise_greater_by_reservoir, local.noise_lesser_by_reservoir)):
                total += value
        result = LocalNambuGenerator(axis, local.signs, local.mode_labels, *accumulated,
            tuple('independent-area-sum:'+s for s in local.reservoir_names))
        if not result.audit()['passed']:
            raise ValueError('spatially projected generator failed canonical commutator/PSD')
        return result

    def integrate_carriers(self, *, initial_amplitudes=None, rtol=2e-10, atol=2e-12):
        omega = optical_carriers(self.inputs.detunings)[1]
        scale = np.sqrt(self.inputs.seed_power_W/(c.HBAR*omega))
        initial = (np.array([scale*np.exp(1j*self.inputs.seed_phase_rad), 0.]) if initial_amplitudes is None
                   else readonly_array(initial_amplitudes))
        if initial.shape != (2,) or not np.isfinite([rtol, atol]).all() or min(rtol, atol) <= 0:
            raise ValueError('two input amplitudes and positive finite ODE tolerances required')
        length = self.inputs.length_m
        ode = solve_ivp(lambda z, y: self.mean_rate(z, y*scale)/scale, (0., length), initial/scale,
            method='DOP853', dense_output=True, rtol=rtol, atol=atol, max_step=length/4 if length else np.inf)
        if not ode.success:
            raise ValueError('spatial nonlinear mean propagation failed')
        if length == 0:
            def constant_output(z):
                values = np.asarray(z)
                return np.broadcast_to((initial/scale).reshape((2,)+(1,)*values.ndim), (2,)+values.shape).copy()
            ode.sol = constant_output
        return ode, scale


def stationary_spatial_cell(medium, detector, *, propagation_rtol=2e-9):
    """Nonlinear mean, quantum tangent/noise and bright two-mode readout.

    Conditional stationary-center, fixed-pump, fixed-profile Galerkin model.
    Spatial quadrature convergence does not certify discarded optical modes.
    """
    if not isinstance(medium, StationarySpatialMedium):
        raise TypeError('stationary spatial medium required')
    rf = detector.analysis_axis.omega_rad_s
    axis = GeneratorFrequencyAxis(np.unique(np.r_[-rf, 0., rf]), 'laboratory RF; stationary-center spatial projection')
    trajectory, scale = medium.integrate_carriers()
    total, propagation = adaptive_field_propagation(
        lambda z: medium.local_field(z, trajectory.sol(z)*scale, axis), medium.inputs.length_m,
        rtol=propagation_rtol, atol=propagation_rtol/100)
    covariance, cp = [], []
    for w in rf:
        i, j = [int(np.flatnonzero(axis.omega_rad_s == s*w)[0]) for s in (1, -1)]
        channel = paired_frequency_channel(total, medium.ports.mode_labels, i, j)
        value = channel.apply_covariance(np.eye(8)/2)
        perm = [0, 1, 6, 7, 4, 5, 2, 3]
        covariance.append(value[perm][:, perm])
        cp.append(channel.audit().minimum_cp_eigenvalue)
    beta = trajectory.y[:, -1]*scale
    spectrum = intensity_difference_spectrum(covariance, beta, detector,
        mode_labels=SIDEBAND_MODES, analysis_axis=detector.analysis_axis)
    powers = c.HBAR*optical_carriers(medium.inputs.detunings)[1:]*abs(beta)**2
    return {'transfer': total, 'spectrum': spectrum, 'output_carriers_sqrt_flux': readonly_array(beta),
        'output_powers_W': readonly_array(powers, real=True),
        'probe_power_gain': float(powers[0]/medium.inputs.seed_power_W),
        'conjugate_power_gain': float(powers[1]/medium.inputs.seed_power_W),
        'sideband_covariances': readonly_array(covariance, real=True),
        'minimum_channel_cp_eigenvalue': min(cp), 'quantum_propagation': propagation,
        'mean_ode_evaluations': trajectory.nfev,
        'scope': 'conditional stationary atomic centers, fixed pump and two normalized transverse profiles; no thermal transport or mode-completeness claim'}
