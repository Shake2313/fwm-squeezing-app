"""Conditional two-band finite-seed mean and quantum propagation, fixed pump.

Optical bands are the existing phase-selected probe/conjugate carriers. All
atomic harmonics enter elimination; extra physical optical ports, pump quantum
fluctuations, Doppler and transverse collection are not implicit in this model.
"""

import numpy as np
from scipy.integrate import solve_ivp

from .. import constants as c, core, observables
from ..quantum.contracts import AnalysisFrequencyAxis, GeneratorFrequencyAxis, readonly_array
from ..quantum.periodic_field import PeriodicFieldPorts, eliminate_periodic_atom, paired_frequency_channel, adaptive_field_propagation
from ..quantum.readout import intensity_difference_spectrum
from ..quantum.sidebands import SIDEBAND_MODES
from ..quantum.traveling import constant_segment, compose_segments
from .field import reduced_readout_operators
from .inputs import ReducedPowerInputs
from .model import reduced_pump_system
from .normalization import optical_carriers, reduced_dipoles
from .periodic import reduced_periodic_noise
from .seed_validity import seed_harmonic_hamiltonian


def reduced_periodic_ports(inputs, *, convention='uniform-zeeman-rms'):
    if not isinstance(inputs, ReducedPowerInputs):
        raise TypeError('ReducedPowerInputs required')
    if inputs.phase_mismatch_rad_m != 0:
        raise ValueError('this periodic cell requires zero phase mismatch; spatial grating model not implemented')
    dipoles = reduced_dipoles(convention)
    omega = optical_carriers(inputs.detunings)[1:]
    ops = reduced_readout_operators(dipoles.transition_scales)
    ops[1] = ops[1].conj().T
    g = dipoles.base_dipole_C_m*np.diag(observables.photon_flux_mode_matrix(
        *omega, inputs.uniform_area_m2, inputs.uniform_area_m2))/(2*c.HBAR)
    return PeriodicFieldPorts(('probe', 'conjugate'), [1, -1], ops, g)


def integrate_reduced_carriers(inputs, *, mean_order=4, convention='uniform-zeeman-rms',
                               initial_amplitudes=None, rtol=2e-10, atol=2e-12):
    """Nonlinear d beta_j/dz=-i lambda*g_j Tr(O_j rho_(h_j)(beta)).

    Solve_ivp's state is beta/sqrt(input probe flux); dense output has the same
    normalization. Returns (ODE result, normalization). No linearized T*beta.
    """
    if isinstance(mean_order, bool) or int(mean_order) != mean_order or mean_order < 1:
        raise ValueError('positive integer mean order required')
    ports = reduced_periodic_ports(inputs, convention=convention)
    dipoles = reduced_dipoles(convention)
    omega = optical_carriers(inputs.detunings)[1:]
    scale = np.sqrt(inputs.seed_power_W/(c.HBAR*omega[0]))
    initial = (np.array([scale*np.exp(1j*inputs.seed_phase_rad), 0.]) if initial_amplitudes is None
               else readonly_array(initial_amplitudes))
    if initial.shape != (2,):
        raise ValueError('two finite initial coherent amplitudes required')
    h, reservoirs = reduced_pump_system(dipoles.pump_rabi_rad_s(inputs.pump_power_W, inputs.pump_waist_m),
        inputs.one_photon_rad_s, transit_rate_s_inverse=inputs.transit_rate_s_inverse,
        transition_scales=dipoles.transition_scales)
    l0 = reservoirs.generator(h)
    beat = -c.OMEGA_HF+inputs.two_photon_rad_s
    def rhs(z, normalized):
        v = seed_harmonic_hamiltonian(dipoles, omega, inputs.uniform_area_m2, normalized*scale)
        rho = core.floquet_solve_truncated(l0, core.comm_super(v), core.comm_super(v.conj().T),
            beat, [0.], np.zeros_like(l0), 4, n_f=int(mean_order), return_harmonics=True)[0]
        pol = np.array([np.trace(op@rho[int(mean_order+harmonic)])
                        for op, harmonic in zip(ports.lowering_operators, ports.carrier_harmonics)])
        return -1j*inputs.number_density_m3*inputs.uniform_area_m2*ports.coupling_s_inverse_sqrt_flux*pol/scale
    solution = solve_ivp(rhs, (0., inputs.length_m), initial/scale, method='DOP853',
        dense_output=True, rtol=rtol, atol=atol, max_step=inputs.length_m/4 if inputs.length_m else np.inf)
    if not solution.success:
        raise ValueError('nonlinear carrier propagation failed')
    if inputs.length_m == 0:
        # SciPy's constant dense-output branch may cast complex states to real.
        def constant_output(z):
            values = np.asarray(z)
            return np.broadcast_to((initial/scale).reshape((2,)+(1,)*values.ndim), (2,)+values.shape).copy()
        solution.sol = constant_output
    return solution, scale


def reduced_periodic_cell(inputs, detector, *, segments=8, mean_order=4, response_order=3,
                           convention='uniform-zeeman-rms', propagation_rtol=None):
    """Same nonlinear mean trajectory for local M/D and bright detector readout.

    Quantum fluctuations are linearized and each segment freezes its midpoint
    mean. With propagation_rtol, use adaptive continuous-path T/N integration;
    segments then specifies diagnostic midpoint samples only. Otherwise refine
    the midpoint propagation grid independently of the atomic harmonic cutoffs.
    """
    if isinstance(segments, bool) or int(segments) != segments or segments < 1:
        raise ValueError('positive integer segment count required')
    ports = reduced_periodic_ports(inputs, convention=convention)
    rf = detector.analysis_axis.omega_rad_s
    axis = GeneratorFrequencyAxis(np.unique(np.r_[-rf, 0., rf]), 'rb85-finite-seed-Floquet-quasifrequency')
    trajectory, scale = integrate_reduced_carriers(inputs, mean_order=mean_order, convention=convention)
    midpoints = (np.arange(int(segments))+.5)*inputs.length_m/segments
    beta_midpoints = trajectory.sol(midpoints).T*scale
    total, local_generators, atomic_diagnostics = None, [], []
    for beta in beta_midpoints:
        atom = reduced_periodic_noise(inputs, beta, mean_order=mean_order, convention=convention)
        local = eliminate_periodic_atom(atom, ports, axis, response_order=response_order,
            linear_density_m_inverse=inputs.number_density_m3*inputs.uniform_area_m2)
        segment = constant_segment(local, inputs.length_m/segments)
        total = segment if total is None else compose_segments(total, segment)
        local_generators.append(local)
        atomic_diagnostics.append(atom.diagnostics)
    propagation = {'method': 'constant-midpoint', 'segments': int(segments)}
    if propagation_rtol is not None:
        def local_at_z(z):
            atom = reduced_periodic_noise(inputs, trajectory.sol(z)*scale, mean_order=mean_order, convention=convention)
            return eliminate_periodic_atom(atom, ports, axis, response_order=response_order,
                linear_density_m_inverse=inputs.number_density_m3*inputs.uniform_area_m2)
        total, propagation = adaptive_field_propagation(local_at_z, inputs.length_m,
            rtol=propagation_rtol, atol=propagation_rtol/100)
    channels, covariances = [], []
    # General converter orders p+,c+,p-,c-; existing detector uses p+,c-,p-,c+.
    permutation = np.array([0, 1, 6, 7, 4, 5, 2, 3])
    for w in rf:
        idx = int(np.flatnonzero(axis.omega_rad_s == w)[0])
        mirror = int(np.flatnonzero(axis.omega_rad_s == -w)[0])
        channel = paired_frequency_channel(total, ports.mode_labels, idx, mirror)
        covariance = channel.apply_covariance(np.eye(8)/2)
        covariances.append(covariance[permutation][:, permutation])
        channels.append(channel)
    beta_out = trajectory.y[:, -1]*scale
    spectrum = intensity_difference_spectrum(covariances, beta_out, detector,
        mode_labels=SIDEBAND_MODES, analysis_axis=detector.analysis_axis)
    powers = c.HBAR*optical_carriers(inputs.detunings)[1:]*abs(beta_out)**2
    return {'transfer': total, 'channels': tuple(channels), 'spectrum': spectrum,
        'output_carrier_amplitudes_sqrt_flux': readonly_array(beta_out),
        'output_coherent_powers_W': readonly_array(powers, real=True),
        'probe_power_gain': float(powers[0]/inputs.seed_power_W),
        'conjugate_power_gain': float(powers[1]/inputs.seed_power_W),
        'midpoint_z_m': readonly_array(midpoints, real=True),
        'midpoint_carriers_sqrt_flux': readonly_array(beta_midpoints),
        'local_generators': tuple(local_generators), 'atomic_diagnostics': tuple(atomic_diagnostics),
        'sideband_covariances': readonly_array(covariances, real=True),
        'mean_ode_evaluations': trajectory.nfev,
        'quantum_propagation': propagation,
        'scope': 'conditional two retained optical bands, nonlinear seed/conjugate mean, fixed classical pump, linearized quantum fluctuations'}
