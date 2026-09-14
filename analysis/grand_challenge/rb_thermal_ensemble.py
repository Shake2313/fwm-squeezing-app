"""Evidence-preserving path reuse for a conditional thermal Rb atomic stream.

Cache physical atomic calculations, never boundary-rate-weighted spectra or
pre-bound convergence certificates. Recompute evidence for each current path.
"""

from collections.abc import Mapping
from dataclasses import asdict, dataclass, field
from functools import cached_property
import hashlib
import itertools
import json
from pathlib import Path
import platform
import sys
import time

import numpy as np
import scipy

from gabes import constants as c
from gabes.fwm_quantum.inputs import ReducedPowerInputs
from gabes.fwm_quantum.kinetic import CarrierGeometry
from gabes.fwm_quantum.smooth_transport import smooth_rb_problem
from gabes.fwm_quantum import transport_ensemble as stream
from gabes.quantum.contracts import AnalysisFrequencyAxis, readonly_array
from gabes.quantum.inflow import maxwell_box_inflow
from gabes.quantum.transport import BallisticPath
from .adjoint_transport_audit import errors as path_errors, METRICS, BUDGETS
from .smooth_transport_audit import fixture


ROOT = Path(__file__).resolve().parents[2]
STREAM_BUDGET = .05
DEFAULT_POWERS = (0, 1, 2)
DEFAULT_SEEDS = (11, 211, 811)
ERROR_DEFINITION = ('Maximum relative Frobenius error per RF AND named source; '
    '128*machine-epsilon times the stated SI dimensional scale is the dark floor; '
    'mean rows use Euclidean norm. No clipping or occupancy renormalization.')


def encode(value):
    if isinstance(value, AnalysisFrequencyAxis):
        return {'analysis_axis_omega_rad_s': value.omega_rad_s.tolist()}
    if isinstance(value, np.ndarray):
        if np.iscomplexobj(value):
            return {'real': value.real.tolist(), 'imag': value.imag.tolist()}
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, Mapping):
        return {str(k): encode(v) for k, v in value.items()}
    if isinstance(value, (tuple, list)):
        return [encode(v) for v in value]
    return value


def decode(value):
    if isinstance(value, dict):
        if set(value) == {'real', 'imag'}:
            return np.asarray(value['real'])+1j*np.asarray(value['imag'])
        if set(value) == {'analysis_axis_omega_rad_s'}:
            return AnalysisFrequencyAxis(value['analysis_axis_omega_rad_s'])
        return {k: decode(v) for k, v in value.items()}
    if isinstance(value, list):
        return [decode(v) for v in value]
    return value


def digest(value):
    return hashlib.sha256(json.dumps(encode(value), sort_keys=True,
        allow_nan=False, separators=(',', ':')).encode()).hexdigest()


def _freeze(value):
    from types import MappingProxyType
    if isinstance(value, Mapping):
        return MappingProxyType({key: _freeze(item) for key, item in value.items()})
    if isinstance(value, list):
        return tuple(_freeze(item) for item in value)
    return value


def consumed_hashes(extra=()):
    files = {Path(__file__).resolve()}
    for name, module in tuple(sys.modules.items()):
        if name == 'gabes' or name.startswith(('gabes.', 'analysis.grand_challenge.')):
            filename = getattr(module, '__file__', None)
            if filename and filename.endswith('.py'):
                files.add(Path(filename).resolve())
    files.update(Path(p).resolve() for p in extra)
    return {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(files)}


@dataclass(frozen=True)
class ThermalRbModel:
    inputs: ReducedPowerInputs
    geometry: CarrierGeometry
    analysis_axis: AnalysisFrequencyAxis
    lower_corner_m: np.ndarray
    upper_corner_m: np.ndarray
    temperature_K: float = 373.
    mass_kg: float = c.MASS_85RB
    boundary_state: np.ndarray = field(default_factory=lambda: np.diag([5/12, 7/12, 0., 0.]))
    pump_center_xy_m: np.ndarray = field(default_factory=lambda: np.zeros(2))
    entry_phase_rad: float = 0.

    def __post_init__(self):
        for key in ('lower_corner_m', 'upper_corner_m', 'pump_center_xy_m'):
            object.__setattr__(self, key, readonly_array(getattr(self, key), real=True))
        object.__setattr__(self, 'boundary_state', readonly_array(self.boundary_state))
        if not isinstance(self.analysis_axis, AnalysisFrequencyAxis):
            raise TypeError('Explicit laboratory RF axis required')
        if not np.isfinite(self.entry_phase_rad):
            raise ValueError('Finite common entry phase required')
        # Validate real box/rates and the native continuous physical problem.
        probe = self.inflow(0, 11).path(0)
        self.problem(probe)
        _physical_exit(self.boundary_state)

    def inflow(self, power, seed):
        return maxwell_box_inflow(self.lower_corner_m, self.upper_corner_m,
            temperature_K=self.temperature_K, mass_kg=self.mass_kg,
            density_m3=self.inputs.number_density_m3, points_per_face_power=power, seed=seed,
            source='Assumed unpolarized open-reservoir Rb column; explicit Maxwell boundary, not a measured cell wall')

    def problem(self, path):
        return smooth_rb_problem(self.inputs, self.geometry, path, self.analysis_axis,
            entry_phase_rad=self.entry_phase_rad, boundary_state=self.boundary_state,
            pump_center_xy_m=self.pump_center_xy_m)

    def identity(self):
        return {'schema': 'conditional-open-column-Rb-model-v1', 'inputs_SI': asdict(self.inputs),
            'wavevectors_rad_m': self.geometry.wavevectors_rad_m,
            'geometry_provenance': self.geometry.source,
            'analysis_axis': self.analysis_axis, 'lower_corner_m': self.lower_corner_m,
            'upper_corner_m': self.upper_corner_m, 'temperature_K': self.temperature_K,
            'mass_kg': self.mass_kg, 'boundary_state': self.boundary_state,
            'pump_center_xy_m': self.pump_center_xy_m, 'entry_phase_rad': self.entry_phase_rad,
            'convention': 'uniform-zeeman-rms; continuous pump-only Gaussian; uniform transverse readouts; unpolarized inflow'}

    @cached_property
    def _contract(self):
        p = self.problem(self.inflow(0, 11).path(0))
        q = self.geometry.wavevectors_rad_m[1:]-self.geometry.wavevectors_rad_m[0]
        _, bp, bc = self.geometry.atomic_frequencies(self.inputs, np.zeros((1, 3)))
        return stream.PumpOnlyConvention(self.analysis_axis, self.inputs.number_density_m3,
            np.tile(p['metadata']['coupling_s_inverse_sqrt_flux'], 2),
            np.linalg.norm(p['readouts'], axis=(-2, -1)),
            ('atomic_inflow',)+tuple('jump:'+channel.name for channel in p['reservoirs'].channels),
            np.array([bp[0], bc[0], -bp[0], -bc[0]]), np.concatenate([q, -q]),
            'conditional-open-column-reduced-Rb-v1', 'smooth',
            'Explicit square open-reservoir column, fixed RMS Rb and continuous Gaussian pump; assumed entry state',
            'One common laboratory beat with actual nonclosed spatial wavevectors and entry phases')

    def convention(self):
        return self._contract


def default_model():
    inputs, geometry, _, axis = fixture()
    side = np.sqrt(inputs.uniform_area_m2)
    lengths = np.array([side, side, inputs.length_m])
    return ThermalRbModel(inputs, geometry, axis, -lengths/2, lengths/2)


@dataclass(frozen=True)
class SolverSpec:
    method: str
    resolution: float
    parameters: Mapping = field(default_factory=dict)

    def __post_init__(self):
        if not isinstance(self.method, str) or not self.method.strip():
            raise ValueError('Named solver method required')
        if not np.isfinite(self.resolution) or self.resolution <= 0:
            raise ValueError('Positive finite numerical resolution required')
        if not isinstance(self.parameters, Mapping):
            raise TypeError('Explicit solver parameters required')
        # Copy through JSON so later caller mutations cannot alter the request.
        parameters = json.loads(json.dumps(encode(self.parameters), allow_nan=False))
        object.__setattr__(self, 'parameters', _freeze(parameters))

    def identity(self):
        return {'method': self.method, 'resolution': self.resolution, 'parameters': dict(self.parameters)}


@dataclass(frozen=True)
class PathPlan:
    primary: tuple[SolverSpec, ...]
    reference: tuple[SolverSpec, ...]

    def __post_init__(self):
        for key, count in (('primary', 3), ('reference', 2)):
            specs = tuple(getattr(self, key))
            if len(specs) != count or not all(isinstance(s, SolverSpec) for s in specs):
                raise ValueError('Three primary and two independent reference calculations required')
            if len({s.method for s in specs}) != 1 or any(a.resolution >= b.resolution for a, b in zip(specs[:-1], specs[1:])):
                raise ValueError('Each solver family must have strictly increasing resolutions')
            object.__setattr__(self, key, specs)
        if self.primary[0].method == self.reference[0].method:
            raise ValueError('Independent reference must name a distinct mathematical solver')

    def identity(self):
        return {key: [s.identity() for s in getattr(self, key)] for key in ('primary', 'reference')}


def _physical_exit(value):
    rho = np.asarray(value, complex)
    if (rho.shape != (4, 4) or not np.isfinite(rho).all()
            or abs(np.trace(rho)-1) > 2e-7 or np.linalg.norm(rho-rho.conj().T) > 2e-7
            or np.linalg.eigvalsh((rho+rho.conj().T)/2).min() < -2e-7):
        raise ValueError('Finite Hermitian trace-one positive Rb density required')


def _validated_packet(model, path, raw):
    if not isinstance(raw, Mapping):
        raise TypeError('Solver must return a complete packet mapping')
    p = model.problem(path)
    packet = dict(raw)
    if 'metadata' in packet and digest(packet['metadata']) != digest(p['metadata']):
        raise ValueError('Cached or computed packet physical metadata differs from the requested model/path')
    packet.setdefault('metadata', p['metadata'])
    packet.setdefault('analysis_axis', model.analysis_axis)
    packet.setdefault('frequencies_rad_s', p['frequencies_rad_s'])
    packet.setdefault('residence_time_s', path.residence_time_s)
    packet.setdefault('source_names', ('atomic_inflow',)+tuple('jump:'+ch.name for ch in p['reservoirs'].channels))
    for key in METRICS:
        if key != 'mean_outer':
            packet[key] = readonly_array(packet[key])
    mean = packet['mean_pulse']
    outer = mean[:, :, None]*mean[:, None, :].conj()
    if 'mean_outer' in packet and not np.allclose(packet['mean_outer'], outer, rtol=1e-12, atol=0):
        raise ValueError('Mean outer must be formed from this same packet mean')
    packet['mean_outer'] = readonly_array(outer)
    _physical_exit(packet['exit_state'])
    # This validates actual shape, phases/frequencies, source closure, PSD and
    # Hermiticity. It does not turn numerical convergence into a boolean claim.
    stream.packet_digest(path, packet, model.convention())
    passed = packet.get('audit', {}).get('passed', True) is True
    packet['audit'] = {'passed': passed,
        'scope': 'Measured local density/ordered covariance/source/phase checks; convergence evidence is separate'}
    return packet


class PathCache:
    """Exclusive immutable JSON records of unweighted physical path solves."""

    def __init__(self, directory, source_hashes):
        self.directory = Path(directory)
        self.source_hashes = _freeze(dict(source_hashes))
        if not self.source_hashes:
            raise ValueError('Nonempty dependency source manifest required')
        self.environment = {'python': platform.python_version(), 'numpy': np.__version__, 'scipy': scipy.__version__}
        self._check_sources()

    def _check_sources(self):
        for name, expected in self.source_hashes.items():
            path = Path(name)
            path = path if path.is_absolute() else ROOT/path
            if hashlib.sha256(path.read_bytes()).hexdigest() != expected:
                raise ValueError('Dependency source changed; create a cache request with the current manifest')

    def identity(self, model, path, spec):
        return {'schema': 'Rb-unweighted-path-calculation-v1', 'model': model.identity(),
            'physical_path': {'entry_position_m': path.entry_position_m, 'velocity_m_s': path.velocity_m_s,
                              'residence_time_s': path.residence_time_s},
            'solver': spec.identity(), 'source_hashes': self.source_hashes, 'environment': self.environment}

    def key(self, model, path, spec):
        return digest(self.identity(model, path, spec))

    def get(self, model, path, spec, provider=None):
        self._check_sources()
        identity = self.identity(model, path, spec)
        key = digest(identity)
        target = self.directory/(key+'.json')
        if target.exists():
            record = json.loads(target.read_text(encoding='utf-8'))
            if (record.get('key') != key or digest(record.get('identity')) != key
                    or record.get('payload_digest') != digest(record.get('packet'))):
                raise ValueError('Corrupt or mismatched immutable path cache record')
            packet = _validated_packet(model, path, decode(record['packet']))
            self._check_sources()
            return packet, {'key': key, 'hit': True, 'payload_digest': record['payload_digest'],
                            'path': str(target), 'original_path_source': record['original_path_source']}
        if provider is None:
            raise FileNotFoundError('No computed path cache record: '+key)
        started = time.monotonic()
        packet = _validated_packet(model, path, provider(model, path, spec))
        self._check_sources()
        record = {'key': key, 'identity': identity, 'original_path_source': path.source,
            'packet': packet, 'payload_digest': digest(packet), 'elapsed_s': time.monotonic()-started}
        payload = json.dumps(encode(record), ensure_ascii=False, indent=2, allow_nan=False)+'\n'
        self.directory.mkdir(parents=True, exist_ok=True)
        # Exclusive final creation: even an interrupted partial record is never
        # silently overwritten or accepted. Corrupt records fail on every read.
        with target.open('x', encoding='utf-8') as handle:
            handle.write(payload)
        return packet, {'key': key, 'hit': False, 'payload_digest': record['payload_digest'],
                        'path': str(target), 'original_path_source': path.source}


def build_path_evidence(model, path, plan, cache, provider=None):
    primary, references, records = [], [], []
    for family, specs in ((primary, plan.primary), (references, plan.reference)):
        for spec in specs:
            packet, record = cache.get(model, path, spec, provider)
            family.append(packet)
            records.append(record)
    p, convention = model.problem(path), model.convention()
    ids = [stream.packet_digest(path, packet, convention) for packet in primary]
    comparisons = []
    for j in range(2):
        comparisons.append(stream.NumericalComparison('path_refinement',
            path_errors(primary[j+1], primary[j], p), dict.fromkeys(METRICS, BUDGETS['primary_refinement']),
            ids[j], ids[j+1], (plan.primary[j].resolution, plan.primary[j+1].resolution),
            ERROR_DEFINITION, 'Actual cached unweighted primary matrices, rebound to current path source and convention'))
    independent_id = plan.reference[-1].method+':'+digest({k: references[-1][k] for k in METRICS})
    comparisons.append(stream.NumericalComparison('independent_reference',
        path_errors(primary[-1], references[-1], p), dict.fromkeys(METRICS, BUDGETS['independent_reference']),
        independent_id, ids[-1], (0., plan.primary[-1].resolution), ERROR_DEFINITION,
        'Actual independently computed full source/response/mean/exit matrices for this exact physical path'))
    independent_refinement = path_errors(references[-1], references[0], p)
    reference_passed = all(np.isfinite(v) and 0 <= v <= BUDGETS['independent_refinement'] for v in independent_refinement.values())
    all_audits = all(packet['audit']['passed'] for packet in primary+references)
    proof = stream.ConvergenceEvidence(ids[-1], 'smooth', tuple(comparisons),
        'Five persisted calculations; actual errors recomputed and all candidate IDs rebound to the current boundary node')
    reasons = stream._evidence_reasons(proof, ids[-1], 'smooth')
    if not reference_passed:
        reasons.append('Independent reference refinement exceeds its unchanged budget')
    if not all_audits:
        reasons.append('A primary or independent local audit failed')
    final = dict(primary[-1])
    final['audit'] = {'passed': not reasons, 'scope': 'Local physical checks plus all five actual numerical calculations'}
    ledger = {'physical_path_key': digest(cache.identity(model, path, plan.primary[-1])),
        'current_path_source': path.source, 'target_digest': ids[-1], 'cache_records': records,
        'primary_refinement': [dict(item.errors) for item in comparisons[:2]],
        'independent_reference': dict(comparisons[-1].errors), 'independent_refinement': independent_refinement,
        'passed': not reasons, 'reasons': reasons}
    return stream.SuppliedPathPacket(final, convention, proof), ledger


def run_grid(model, power, seed, plan, cache, provider=None):
    inflow, convention = model.inflow(power, seed), model.convention()
    ledgers = []
    callbacks = 0
    active_index, factory_completed = None, False

    def factory(path):
        nonlocal callbacks, active_index, factory_completed
        active_index, factory_completed = callbacks, False
        callbacks += 1
        supplied, ledger = build_path_evidence(model, path, plan, cache, provider)
        ledgers.append(ledger)
        factory_completed = True
        return supplied

    try:
        candidate = stream.aggregate_pump_only_stream(inflow, convention, factory)
    except (OSError, ValueError, TypeError, KeyError, RuntimeError, ArithmeticError) as error:
        final_stage = factory_completed and callbacks == len(inflow.rate_s_inverse)
        failed_index = None if final_stage else active_index
        candidate = {'certified': False, 'path_evidence_passed': False, 'spectra': None,
            'candidate_digest': None, 'path_count': len(inflow.rate_s_inverse),
            'consumed_path_count': len(ledgers), 'failed_path_index': failed_index,
            'last_attempted_path_index': active_index,
            'failure_stage': 'last-packet validation or final accumulation' if final_stage else 'path acquisition or validation',
            'failed_path_rate_s_inverse': None if failed_index is None else float(inflow.rate_s_inverse[failed_index]),
            'total_arrival_rate_s_inverse': inflow.total_arrival_rate_s_inverse,
            'reasons': (type(error).__name__+': '+str(error),)}
    row = {'power': power, 'seed': seed, 'model_digest': digest(model.identity()),
        'plan_digest': digest(plan.identity()), 'path_count': len(inflow.rate_s_inverse),
        'callback_count': callbacks, 'cache_hits': sum(r['hit'] for item in ledgers for r in item['cache_records']),
        'cache_misses': sum(not r['hit'] for item in ledgers for r in item['cache_records']),
        'path_ledgers': ledgers, 'candidate': candidate,
        'mean_occupancy': inflow.mean_occupancy, 'equilibrium_nV': inflow.equilibrium_atom_number,
        'occupancy_over_nV': inflow.mean_occupancy/inflow.equilibrium_atom_number,
        'comparison_scales': stream_scales(model, inflow.total_arrival_rate_s_inverse)}
    row['content_digest'] = row_digest(row)
    return candidate, row


def stream_scales(model, rate):
    p = model.problem(model.inflow(0, 11).path(0))
    tau = float(np.max(model.upper_corner_m-model.lower_corner_m)/np.sqrt(c.KB*model.temperature_K/model.mass_kg))
    op = float(np.linalg.norm(p['readouts'], axis=(-2, -1)).max())
    drive = float(np.linalg.norm(p['drives'], axis=(-2, -1)).max())
    values = {key: rate*tau*tau*op*op for key in stream.STREAM_METRICS}
    values['retarded_response'] = rate*tau*tau*op*drive
    return values


def stream_errors(candidate, reference, scales):
    result = {}
    for key in stream.STREAM_METRICS:
        a, b = np.asarray(candidate[key]), np.asarray(reference[key])
        if a.shape != b.shape or not np.isfinite(a).all() or not np.isfinite(b).all():
            raise ValueError('Complete matching finite source/RF matrix arrays required')
        if not np.isfinite(scales[key]) or scales[key] <= 0:
            raise ValueError('Positive finite SI comparison scales required')
        with np.errstate(over='ignore', invalid='ignore', divide='ignore'):
            denominator = np.maximum(np.linalg.norm(b, axis=(-2, -1)), 128*np.finfo(float).eps*scales[key])
            value = float(np.max(np.linalg.norm(a-b, axis=(-2, -1))/denominator))
        if not np.isfinite(value) or value < 0:
            raise ValueError('Nonfinite or negative derived stream comparison error')
        result[key] = value
    return result


def row_digest(row):
    return digest({key: value for key, value in row.items() if key != 'content_digest'})


def convergence_gate(rows, *, powers=DEFAULT_POWERS, seeds=DEFAULT_SEEDS):
    reasons, comparisons = [], []
    powers, seeds = tuple(powers), tuple(seeds)
    if len(powers) < 3 or tuple(sorted(set(powers))) != powers or len(set(seeds)) < 3 or len(set(seeds)) != len(seeds):
        return {'passed': False, 'reasons': ['At least three increasing powers and three distinct seeds required'], 'comparisons': []}
    indexed = {}
    try:
        for row in rows:
            key = (row['power'], row['seed'])
            if key in indexed or key[0] not in powers or key[1] not in seeds:
                reasons.append('Duplicate or unplanned grid/seed')
            indexed[key] = row
            if row['content_digest'] != row_digest(row):
                reasons.append('Grid content digest changed')
            candidate = row['candidate']
            if (not candidate['path_evidence_passed'] or candidate['spectra'] is None
                    or candidate['consumed_path_count'] != 6*2**row['power']
                    or len(row['path_ledgers']) != 6*2**row['power']
                    or not all(item['passed'] for item in row['path_ledgers'])):
                reasons.append('Every physical path requires complete passing numerical evidence')
            else:
                # Existing production check detects altered spectra or metadata.
                stream.certify_pump_only_stream(candidate, None)
            for ledger in row['path_ledgers']:
                if len(ledger['primary_refinement']) != 2:
                    reasons.append('Two measured primary path refinements required')
                measurements = [(values, BUDGETS['primary_refinement']) for values in ledger['primary_refinement']]
                measurements += [(ledger['independent_reference'], BUDGETS['independent_reference']),
                                 (ledger['independent_refinement'], BUDGETS['independent_refinement'])]
                for values, budget in measurements:
                    actual = [values[metric] for metric in METRICS]
                    if not np.isfinite(actual).all() or min(actual) < 0 or max(actual) > budget:
                        reasons.append('Incomplete or failed measured path/reference metrics')
        present = sorted({p for p, seed in indexed})
        if len(present) < 3 or tuple(present) != powers[:len(present)]:
            reasons.append('Three consecutive planned grids required')
        if set(indexed) != set(itertools.product(present, seeds)):
            reasons.append('Every declared seed is required on every evaluated grid')
        if len({row['model_digest'] for row in rows}) != 1 or len({row['plan_digest'] for row in rows}) != 1:
            reasons.append('All grids must share the physical model and numerical plan')
        if reasons:
            return {'passed': False, 'reasons': sorted(set(reasons)), 'comparisons': []}

        def compare(kind, current, previous):
            if kind == 'independent_scramble' and current['candidate']['candidate_digest'] == previous['candidate']['candidate_digest']:
                raise ValueError('Independent scrambles cannot reuse the same candidate computation')
            errors = stream_errors(current['candidate']['spectra'], previous['candidate']['spectra'], current['comparison_scales'])
            return {'kind': kind, 'candidate_grid': [current['power'], current['seed']],
                'reference_grid': [previous['power'], previous['seed']], 'errors': errors,
                'candidate_digest': current['candidate']['candidate_digest'],
                'reference_digest': previous['candidate']['candidate_digest'],
                'budget': STREAM_BUDGET, 'passed': all(np.isfinite(v) and 0 <= v <= STREAM_BUDGET for v in errors.values())}

        for previous, power in zip(present[-3:-1], present[-2:]):
            for seed in seeds:
                comparisons.append(compare('ensemble_refinement', indexed[power, seed], indexed[previous, seed]))
        for power in present[-2:]:
            for current, previous in itertools.permutations(seeds, 2):
                comparisons.append(compare('independent_scramble', indexed[power, current], indexed[power, previous]))
    except (ValueError, KeyError, TypeError, IndexError) as error:
        reasons.append('Invalid or incomplete grid data: '+str(error))
    if not all(row['passed'] for row in comparisons):
        reasons.append('Actual thermal source/noise/number/response error exceeds 5%')
    return {'passed': not reasons and len(comparisons) == 2*len(seeds)**2,
            'reasons': sorted(set(reasons)), 'comparisons': comparisons}
