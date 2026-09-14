"""Thermal boundary flux and period-averaged marked-Poisson atom pulses.

All paths are first-exit chords in an explicit rectangular box. Rates come
from n f(v) |v.normal| dA dv, without normalizing the computed occupancy.
This module supplies atomic input ensembles, not an optical quantum channel.
"""

from dataclasses import dataclass
import operator

import numpy as np
from scipy.special import ndtri
from scipy.stats import qmc

from gabes.constants import KB
from .contracts import GeneratorFrequencyAxis, readonly_array
from .transport import BallisticPath


def _positive(value, name):
    result = float(value)
    if not np.isfinite(result) or result <= 0:
        raise ValueError(f'{name} must be finite and positive')
    return result


def _integer(value, name, lower, upper):
    try:
        result = operator.index(value)
    except TypeError as error:
        raise ValueError(f'{name} must be an integer') from error
    if isinstance(value, (bool, np.bool_)) or not lower <= result <= upper:
        raise ValueError(f'{name} must be between {lower} and {upper}')
    return result


@dataclass(frozen=True)
class ThermalInflowQuadrature:
    """Immutable positive boundary-rate nodes; face IDs are x-,x+,y-,y+,z-,z+.

    Use :func:`maxwell_box_inflow` to construct. Arrays have one row per
    incoming chord; ``rate_s_inverse`` already includes density and face area.
    Material walls, collisions, re-entry correlations and a beam aperture are
    not supplied by a finite-box boundary condition.
    """

    lower_corner_m: np.ndarray
    upper_corner_m: np.ndarray
    temperature_K: float
    mass_kg: float
    density_m3: float
    entry_position_m: np.ndarray
    velocity_m_s: np.ndarray
    residence_time_s: np.ndarray
    rate_s_inverse: np.ndarray
    face_index: np.ndarray
    points_per_face_power: int
    seed: int
    source: str

    def __post_init__(self):
        for name in ('lower_corner_m', 'upper_corner_m', 'entry_position_m', 'velocity_m_s',
                     'residence_time_s', 'rate_s_inverse'):
            object.__setattr__(self, name, readonly_array(getattr(self, name), real=True))
        lo, hi = self.lower_corner_m, self.upper_corner_m
        if lo.shape != (3,) or hi.shape != (3,) or np.any(hi <= lo):
            raise ValueError('finite rectangular box needs three strictly increasing bounds')
        for name in ('temperature_K', 'mass_kg', 'density_m3'):
            object.__setattr__(self, name, _positive(getattr(self, name), name))
        power = _integer(self.points_per_face_power, 'points_per_face_power', 0, 20)
        _integer(self.seed, 'seed', 0, 2**32-1)
        if not isinstance(self.source, str) or not self.source.strip():
            raise ValueError('inflow provenance required')
        face = readonly_array(self.face_index, real=True)
        count = 6*2**power
        if (face.shape != (count,) or np.any(face != np.rint(face)) or
                np.any(face < 0) or np.any(face > 5)):
            raise ValueError('one valid boundary face index per quadrature node required')
        face = face.astype(np.int64)
        if np.any(np.bincount(face, minlength=6) != 2**power):
            raise ValueError('all six faces must have the declared quadrature count')
        object.__setattr__(self, 'face_index', np.frombuffer(face.tobytes(), dtype=np.int64))
        entry, velocity, tau, rate = (self.entry_position_m, self.velocity_m_s,
                                     self.residence_time_s, self.rate_s_inverse)
        if (entry.shape != (count, 3) or velocity.shape != (count, 3) or
                tau.shape != (count,) or rate.shape != (count,)):
            raise ValueError('matching entry/velocity (nodes,3) and residence/rate (nodes,) required')
        if np.any(tau <= 0) or np.any(rate <= 0):
            raise ValueError('positive residence times and boundary rates required')
        normal = face//2
        inward = np.where(face % 2 == 0, 1, -1)
        rows = np.arange(count)
        expected_face = np.where(inward == 1, lo[normal], hi[normal])
        if (np.any(entry < lo) or np.any(entry > hi) or
                not np.array_equal(entry[rows, normal], expected_face) or
                np.any(velocity[rows, normal]*inward <= 0)):
            raise ValueError('entry must lie on its declared box face with inward velocity')
        distance = np.where(velocity > 0, hi-entry, entry-lo)
        flight = np.full_like(velocity, np.inf)
        np.divide(distance, np.abs(velocity), out=flight, where=velocity != 0)
        if not np.allclose(tau, flight.min(axis=1), rtol=1e-12, atol=0):
            raise ValueError('residence must be the true first-exit box chord')
        lengths = hi-lo
        sigma = np.sqrt(KB*self.temperature_K/self.mass_kg)
        expected_rate = self.density_m3*(np.prod(lengths)/lengths[normal])*sigma/np.sqrt(2*np.pi)/2**power
        if not np.allclose(rate, expected_rate, rtol=1e-12, atol=0):
            raise ValueError('node rates must equal explicit Maxwell boundary flux, without fitting')

    @property
    def total_arrival_rate_s_inverse(self):
        return float(self.rate_s_inverse.sum())

    @property
    def equilibrium_atom_number(self):
        return float(self.density_m3 * np.prod(self.upper_corner_m-self.lower_corner_m))

    @property
    def mean_occupancy(self):
        return float(self.rate_s_inverse @ self.residence_time_s)

    def path(self, index):
        index = _integer(index, 'path index', 0, len(self.rate_s_inverse)-1)
        return BallisticPath(self.entry_position_m[index], self.velocity_m_s[index],
            self.residence_time_s[index], f'{self.source}; face={self.face_index[index]}, node={index}')

    def occupation_moments(self):
        """Unnormalized phase-space moments divided by analytic nV.

        Uniform-age position moments are integrated analytically on each
        chord. In particular these do not normalize away an occupancy error.
        """
        weight = self.rate_s_inverse*self.residence_time_s/self.equilibrium_atom_number
        velocity = self.velocity_m_s
        midpoint = (self.lower_corner_m+self.upper_corner_m)/2
        entry = self.entry_position_m-midpoint
        displacement = self.residence_time_s[:, None]*velocity
        center = entry+displacement/2
        return {'occupancy_over_nV': float(weight.sum()),
            'velocity_mean_m_s': np.einsum('n,ni->i', weight, velocity),
            'velocity_second_m2_s2': np.einsum('n,ni,nj->ij', weight, velocity, velocity),
            'speed_mean_m_s': float(weight @ np.linalg.norm(velocity, axis=1)),
            'speed_fourth_m4_s4': float(weight @ (np.sum(velocity**2, axis=1)**2)),
            'position_centered_mean_m': np.einsum('n,ni->i', weight, center),
            'position_centered_second_m2': np.einsum('n,ni,nj->ij', weight, center, center)
                + np.einsum('n,ni,nj->ij', weight, displacement, displacement)/12}


def maxwell_box_inflow(lower_corner_m, upper_corner_m, *, temperature_K, mass_kg,
                       density_m3, points_per_face_power=10, seed=0, source):
    """Direct six-face scrambled-Sobol quadrature of Maxwell incoming flux.

    Each face uses 2**points_per_face_power points in five dimensions: two
    uniform positions, one inward Rayleigh normal speed, and two Gaussian
    tangential velocities (sigma=sqrt(kB*T/m)). Each positive node has rate
    n*A*sigma/sqrt(2*pi)/N. Seed and power are reproducible quadrature controls,
    not physical inputs. Refine powers and independent seeds for new geometry
    or wavepacket observables: moment convergence alone cannot certify them.
    """
    lo, hi = (readonly_array(x, real=True) for x in (lower_corner_m, upper_corner_m))
    if lo.shape != (3,) or hi.shape != (3,) or np.any(hi <= lo):
        raise ValueError('finite rectangular box needs three strictly increasing bounds')
    temp = _positive(temperature_K, 'temperature_K')
    mass = _positive(mass_kg, 'mass_kg')
    density = _positive(density_m3, 'density_m3')
    power = _integer(points_per_face_power, 'points_per_face_power', 0, 20)
    seed = _integer(seed, 'seed', 0, 2**32-1)
    if not isinstance(source, str) or not source.strip():
        raise ValueError('independent density, temperature and geometry provenance required')
    lengths = hi-lo
    sigma = np.sqrt(KB*temp/mass)
    count = 2**power
    entries, velocities, rates, times, faces = [], [], [], [], []
    for face, child in enumerate(np.random.SeedSequence(seed).spawn(6)):
        normal = face//2
        tangential = [j for j in range(3) if j != normal]
        inward = 1 if face % 2 == 0 else -1
        # A dyadic midpoint shifts all finite 30-bit Sobol coordinates into
        # the open cube. It avoids corners and infinite Gaussian endpoints;
        # the discarded probability per tail is below 2**-31.
        u = qmc.Sobol(5, scramble=True, bits=30, seed=np.random.default_rng(child)).random_base2(power)
        u += 2.**-31
        entry = np.empty((count, 3))
        entry[:, normal] = lo[normal] if inward == 1 else hi[normal]
        entry[:, tangential] = lo[tangential]+u[:, :2]*lengths[tangential]
        velocity = np.empty_like(entry)
        velocity[:, normal] = inward*sigma*np.sqrt(-2*np.log1p(-u[:, 2]))
        velocity[:, tangential] = sigma*ndtri(u[:, 3:])
        distance = np.where(velocity > 0, hi-entry, entry-lo)
        flight = np.full_like(velocity, np.inf)
        np.divide(distance, np.abs(velocity), out=flight, where=velocity != 0)
        tau = flight.min(axis=1)
        rate = density*np.prod(lengths[tangential])*sigma/np.sqrt(2*np.pi)/count
        if not np.isfinite(rate) or rate <= 0 or np.any(~np.isfinite(tau)) or np.any(tau <= 0):
            raise ValueError('box flux or chord cannot be represented with finite positive SI values')
        entries.append(entry); velocities.append(velocity); times.append(tau)
        rates.append(np.full(count, rate)); faces.append(np.full(count, face, dtype=np.int64))
    arrays = [readonly_array(np.concatenate(a), real=True) for a in (entries, velocities, times, rates)]
    face_ids = np.concatenate(faces)
    face_ids = np.frombuffer(face_ids.tobytes(), dtype=np.int64)
    return ThermalInflowQuadrature(lo, hi, temp, mass, density, *arrays, face_ids, power, seed, source)


@dataclass(frozen=True)
class EntryPhaseQuadrature:
    """One common lab RF beat phase, with normalized positive phase weights."""

    phase_rad: np.ndarray
    probability: np.ndarray
    source: str

    def __post_init__(self):
        phase, probability = (readonly_array(x, real=True) for x in (self.phase_rad, self.probability))
        if phase.ndim != 1 or not phase.size or probability.shape != phase.shape:
            raise ValueError('one scalar common phase per mark and matching probabilities required')
        if np.any(probability <= 0) or not np.isclose(probability.sum(), 1, rtol=1e-13, atol=1e-15):
            raise ValueError('phase probabilities must be positive and sum to one')
        if not isinstance(self.source, str) or not self.source.strip():
            raise ValueError('phase provenance required')
        object.__setattr__(self, 'phase_rad', phase)
        object.__setattr__(self, 'probability', probability)


def uniform_entry_phases(count, *, offset_rad=0.):
    """Periodic trapezoid rule for one uniformly sampled lab RF beat phase."""
    count = _integer(count, 'phase count', 1, 2**20)
    offset = float(offset_rad)
    if not np.isfinite(offset):
        raise ValueError('finite common phase offset required')
    offset = offset % (2*np.pi)
    return EntryPhaseQuadrature(offset+2*np.pi*np.arange(count)/count,
        np.full(count, 1/count), 'uniform common RF beat phase; periodic trapezoid quadrature')


def entry_phase_factors(path, common_phase_rad, wavevectors_rad_m, temporal_harmonics,
                        *, offsets_rad=None):
    """exp[i(k_j.r_entry + h_j*phi + offset_j)] with a SINGLE common phi.

    k_j and offsets are fixed geometry inputs, h_j are signed integer RF
    harmonics. Age evolution k_j.v*a + h_j*omega_RF*a belongs to the caller's
    Hamiltonian. This helper never invents independent random spatial phases.
    """
    if not isinstance(path, BallisticPath):
        raise TypeError('explicit ballistic path required')
    phase = np.asarray(common_phase_rad)
    if phase.ndim != 0 or np.iscomplexobj(phase) or not np.isfinite(phase):
        raise ValueError('one finite real common RF phase required')
    k = readonly_array(wavevectors_rad_m, real=True)
    harmonics = readonly_array(temporal_harmonics, real=True)
    if k.ndim != 2 or k.shape[1] != 3 or harmonics.shape != (len(k),) or not len(k):
        raise ValueError('wavevectors (terms,3) and one harmonic per term required')
    if np.any(harmonics != np.rint(harmonics)):
        raise ValueError('temporal harmonics must be integers of one common RF beat')
    offsets = np.zeros(len(k)) if offsets_rad is None else readonly_array(offsets_rad, real=True)
    if offsets.shape != (len(k),):
        raise ValueError('one fixed spatial phase offset per term required')
    return readonly_array(np.exp(1j*(k@path.entry_position_m+harmonics*float(phase)+offsets)))


def _ordered_matrix(value, shape, label):
    matrix = readonly_array(value)
    if matrix.shape != shape:
        raise ValueError(f'{label} must have shape (frequency, readout, readout)')
    scale = np.maximum(np.max(np.abs(matrix), axis=(-2, -1)), np.finfo(float).tiny)
    scaled = matrix/scale[:, None, None]
    if np.max(np.abs(scaled-scaled.swapaxes(-1, -2).conj())) > 1e-9:
        raise ValueError(f'{label} must be Hermitian')
    if np.linalg.eigvalsh((scaled+scaled.swapaxes(-1, -2).conj())/2).min() < -1e-9:
        raise ValueError(f'{label} must be positive semidefinite')
    return matrix


def average_marked_poisson(inflow, phases, wavepacket_factory, *, source):
    """Stream per-path/per-phase raw wavepacket moments before any subtraction.

    ``wavepacket_factory(BallisticPath, scalar_common_phase_rad)`` must return
    the transport.wavepacket schema: GeneratorFrequencyAxis, connected
    ``greater``/``lesser`` [nf,p,p], ``mean_pulse`` [nf,p], residence_time_s.
    All packets must share frame, grid, readout order/units and coupling
    convention; the latter are the caller's explicitly sourced responsibility.

    Sum_j rate_j Sum_phi p_phi (C_jphi + m_jphi m_jphi^dagger).
    For a periodic lab drive this is the zero-cyclic, period-averaged connected
    stream spectrum. Deterministic periodic mean lines and off-diagonal cyclic
    spectra are excluded. Phase-averaging does not prove full stationarity.
    No optical SQL, detector unit conversion or canonical commutator is implied.
    """
    if not isinstance(inflow, ThermalInflowQuadrature) or not isinstance(phases, EntryPhaseQuadrature):
        raise TypeError('explicit thermal inflow and common entry-phase quadrature required')
    if not callable(wavepacket_factory) or not isinstance(source, str) or not source.strip():
        raise ValueError('wavepacket callback and readout/coupling provenance required')
    axis = None
    internal_greater = internal_lesser = number = mean_sum = None
    total_rate = inflow.total_arrival_rate_s_inverse
    for index, rate in enumerate(inflow.rate_s_inverse):
        path = inflow.path(index)
        for phase, probability in zip(phases.phase_rad, phases.probability):
            packet = wavepacket_factory(path, float(phase))
            packet_axis = packet['frequency_axis']
            if not isinstance(packet_axis, GeneratorFrequencyAxis):
                raise TypeError('wavepacket must declare its GeneratorFrequencyAxis')
            mean = readonly_array(packet['mean_pulse'])
            if mean.ndim != 2 or mean.shape[0] != len(packet_axis.omega_rad_s) or mean.shape[1] == 0:
                raise ValueError('mean_pulse must have shape (frequency, readout)')
            tau = _positive(packet['residence_time_s'], 'wavepacket residence_time_s')
            if not np.isclose(tau, path.residence_time_s, rtol=1e-12, atol=0):
                raise ValueError('wavepacket residence must equal the actual chord residence')
            shape = (len(packet_axis.omega_rad_s), mean.shape[1], mean.shape[1])
            greater = _ordered_matrix(packet['greater'], shape, 'connected greater')
            lesser = _ordered_matrix(packet['lesser'], shape, 'connected lesser')
            if axis is None:
                axis = packet_axis
                internal_greater = np.zeros(shape, complex)
                internal_lesser = np.zeros(shape, complex)
                number = np.zeros(shape, complex)
                mean_sum = np.zeros(mean.shape, complex)
            elif (axis.frame != packet_axis.frame or
                  not np.array_equal(axis.omega_rad_s, packet_axis.omega_rad_s) or
                  shape != number.shape):
                raise ValueError('all marked wavepackets must use the same frame, grid and readout shape')
            weight = rate*probability
            internal_greater += weight*greater
            internal_lesser += weight*lesser
            number += weight*mean[:, :, None]*mean[:, None, :].conj()
            mean_sum += (weight/total_rate)*mean
    return {'frequency_axis': axis, 'greater': readonly_array(internal_greater+number),
        'lesser': readonly_array(internal_lesser+number),
        'internal_greater': readonly_array(internal_greater),
        'internal_lesser': readonly_array(internal_lesser), 'poisson_number': readonly_array(number),
        'arrival_weighted_conditional_mean_pulse': readonly_array(mean_sum),
        'arrival_rate_s_inverse': total_rate, 'mean_occupancy': inflow.mean_occupancy,
        'path_count': len(inflow.rate_s_inverse), 'phase_count': len(phases.phase_rad),
        'source': source, 'inflow_source': inflow.source, 'phase_source': phases.source,
        'scope': 'zero-cyclic period-averaged connected atomic stream spectrum; coherent periodic mean lines excluded; no optical SQL normalization'}
