"""Moving reduced Rb D1 pump characteristic with reciprocal weak optical ports.

This is the undepleted pump-only linearization. Finite seed saturation, thermal
averaging and self-consistent nonlocal Maxwell propagation are not supplied.
"""

import numpy as np

from .. import constants as c, observables
from ..quantum.contracts import AnalysisFrequencyAxis, readonly_array
from ..quantum.segmented_transport import segmented_wavepacket
from ..quantum.transport import BallisticPath
from .field import reduced_readout_operators
from .inputs import ReducedPowerInputs
from .kinetic import CarrierGeometry
from .model import reduced_pump_system
from .normalization import optical_carriers, reduced_dipoles


def reduced_ballistic_problem(inputs, geometry, path, analysis_axis, *, segments,
                               entry_phase_rad=0., boundary_state=None,
                               pump_center_xy_m=(0., 0.), convention='uniform-zeeman-rms'):
    """Freeze only the transverse Gaussian pump amplitude at segment midpoints.

    Entry phase is the ONE signed lab beat phase. Spatial entry phases follow
    the actual wavevectors, including their nonzero loop mismatch. Pump is
    along +z, constant waist, without diffraction or longitudinal depletion.
    Optical outputs are uniform transverse atomic polarization projections.
    """
    if not isinstance(inputs, ReducedPowerInputs) or not isinstance(geometry, CarrierGeometry):
        raise TypeError('declared reduced power inputs and carrier geometry required')
    if not isinstance(path, BallisticPath) or not isinstance(analysis_axis, AnalysisFrequencyAxis):
        raise TypeError('ballistic path and independent laboratory RF axis required')
    if inputs.transit_rate_s_inverse != 0:
        raise ValueError('explicit path boundary replaces phenomenological transit reset; supply zero transit rate')
    if inputs.phase_mismatch_rad_m != 0:
        raise ValueError('explicit wavevectors cannot be combined with scalar mismatch')
    if isinstance(segments, bool) or int(segments) != segments or segments < 1:
        raise ValueError('positive integer segment count required')
    center = readonly_array(pump_center_xy_m, real=True)
    phase = float(entry_phase_rad)
    if center.shape != (2,) or not np.isfinite(phase):
        raise ValueError('finite pump center and entry beat phase required')
    count = int(segments)
    times = np.linspace(0., path.residence_time_s, count+1)
    midpoints = path.position((times[:-1]+times[1:])/2)
    profile = np.exp(-np.sum((midpoints[:, :2]-center)**2, axis=1)/inputs.pump_waist_m**2)
    dipoles = reduced_dipoles(convention)
    rabi = dipoles.pump_rabi_rad_s(inputs.pump_power_W, inputs.pump_waist_m)
    dv, bp, bc = geometry.atomic_frequencies(inputs, path.velocity_m_s[None])
    h0, reservoirs = reduced_pump_system(0., dv[0], transition_scales=dipoles.transition_scales)
    h1, _ = reduced_pump_system(rabi, dv[0], transition_scales=dipoles.transition_scales)
    hs = h0+profile[:, None, None]*(h1-h0)
    lower = reduced_readout_operators(dipoles.transition_scales)
    lower[1] = lower[1].conj().T
    q = geometry.wavevectors_rad_m[1:]-geometry.wavevectors_rad_m[0]
    phases = np.array([phase, -phase])-q@path.entry_position_m
    lower = lower*np.exp(1j*phases)[:, None, None]
    # Exact demodulation along r=r_in+v*a: nu_j=omega_j-omega_0-q_j.v.
    # In particular nu_p+nu_c=-(kp+kc-2k0).v is retained, not closed by hand.
    ops = np.concatenate([lower, lower.conj().swapaxes(-1, -2)])
    offsets = np.array([bp[0], bc[0], -bp[0], -bc[0]])
    freq = analysis_axis.omega_rad_s[:, None]+offsets
    omega = optical_carriers(inputs.detunings)[1:]
    vectors = geometry.wavevectors_rad_m[1:]
    cosine = vectors[:, 2]/np.linalg.norm(vectors, axis=1)
    flux = np.diag(observables.photon_flux_mode_matrix(*omega, inputs.uniform_area_m2, inputs.uniform_area_m2))
    g = dipoles.base_dipole_C_m*flux/(2*c.HBAR*np.sqrt(cosine))
    drives = ops.conj().swapaxes(-1, -2)*np.tile(g, 2)[:, None, None]
    boundary = np.diag([5/12, 7/12, 0., 0.]) if boundary_state is None else readonly_array(boundary_state)
    return {'hamiltonians': readonly_array(hs), 'reservoirs': reservoirs,
        'boundary_state': readonly_array(boundary), 'durations_s': readonly_array(np.diff(times), real=True),
        'readouts': readonly_array(np.broadcast_to(ops, (count,)+ops.shape)),
        'frequencies_rad_s': readonly_array(freq, real=True),
        'drives': readonly_array(np.broadcast_to(drives, (count,)+drives.shape)),
        'drive_frequencies_rad_s': readonly_array(freq, real=True),
        'metadata': {'mode_labels': ('probe', 'conjugate', 'probe_dagger', 'conjugate_dagger'),
            'entry_phase_rad': phase, 'entry_optical_demodulation_phases_rad': phases.tolist(),
            'carrier_offsets_rad_s': offsets.tolist(), 'pump_amplitudes_rad_s': (rabi*profile).tolist(),
            'coupling_s_inverse_sqrt_flux': g.tolist(), 'one_photon_doppler_rad_s': float(dv[0]),
            'loop_convective_frequency_rad_s': float(bp[0]+bc[0]),
            'boundary_source': 'assumed unpolarized ground F populations' if boundary_state is None else 'caller-supplied state',
            'unused_input_fields': ['number_density_m3', 'length_m', 'seed_power_W', 'seed_phase_rad'],
            'scope': 'one moving pump-only reduced RMS Rb atom; finite seed, full Zeeman, ensemble and Maxwell field absent'}}


def reduced_ballistic_wavepacket(inputs, geometry, path, analysis_axis, **kwargs):
    problem = reduced_ballistic_problem(inputs, geometry, path, analysis_axis, **kwargs)
    packet = segmented_wavepacket(**{k: v for k, v in problem.items() if k != 'metadata'})
    return {**packet, 'metadata': problem['metadata'], 'analysis_axis': analysis_axis}


def average_pump_only_entry_phase(packet):
    """Exact uniform average over the ONE laboratory beat phase, with no solve.

    Pump-only H and rho do not depend on entry phase. The four readout charges
    are s=(1,-1,-1,1), so C_jk(phi)=exp(i*(s_j-s_k)*phi) C_jk(0).
    Uniform integration is exactly delta(s_j,s_k), even for nonclosed spatial
    wavevectors. This proof removes redundant phase solves ONLY in this
    pump-only problem; finite seed saturation invalidates the premise.
    """
    if (packet.get('metadata', {}).get('scope') !=
            'one moving pump-only reduced RMS Rb atom; finite seed, full Zeeman, ensemble and Maxwell field absent'):
        raise ValueError('phase theorem requires the declared pump-only Rb packet')
    charges = np.array([1, -1, -1, 1])
    mask = (charges[:, None] == charges[None, :])
    mean = packet['mean_pulse']
    number = mean[:, :, None]*mean[:, None, :].conj()*mask
    result = {'conditional_'+key: readonly_array(packet[key]*mask) for key in
              ('greater', 'lesser', 'greater_by_source', 'lesser_by_source')}
    result.update({'mean_pulse_phase_average': readonly_array(np.zeros_like(mean)),
        'retarded_response_phase_average': readonly_array(packet['retarded_response']*mask),
        'mean_outer_phase_average': readonly_array(number),
        'raw_greater_phase_average': readonly_array(result['conditional_greater']+number),
        'raw_lesser_phase_average': readonly_array(result['conditional_lesser']+number),
        'residence_time_s': packet['residence_time_s'], 'charges': readonly_array(charges, real=True),
        'scope': 'pump-only uniform common entry-phase average; Poisson stream must use averaged raw moments'})
    return result
