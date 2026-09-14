"""Validated pump-only atomic stream on a common laboratory RF axis.

Conditional numerical certification is an evidence contract, not a proof of
the caller's physical model or an optical commutator/SQL certification.
"""

from dataclasses import dataclass
import hashlib
import json
from types import MappingProxyType
from collections.abc import Mapping

import numpy as np

from ..quantum.contracts import AnalysisFrequencyAxis, readonly_array
from ..quantum.inflow import ThermalInflowQuadrature
from ..quantum.transport import BallisticPath


MODE_LABELS = ('probe', 'conjugate', 'probe_dagger', 'conjugate_dagger')
CHARGES = (1, -1, -1, 1)
MASK = np.equal.outer(CHARGES, CHARGES)
MASK.setflags(write=False)
PATH_METRICS = ('greater', 'lesser', 'greater_by_source', 'lesser_by_source',
                'mean_pulse', 'mean_outer', 'retarded_response')
STREAM_METRICS = ('greater', 'lesser', 'greater_by_source', 'lesser_by_source',
                  'poisson_number', 'retarded_response')


def _text(value, name):
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f'{name} must name its provenance/convention')
    return value


def _digest(parts):
    h = hashlib.sha256()
    for value in parts:
        if isinstance(value, np.ndarray):
            array = np.ascontiguousarray(value)
            data = json.dumps([array.dtype.str, array.shape]).encode()+array.tobytes()
        else:
            data = json.dumps(value, sort_keys=True, allow_nan=False, separators=(',', ':')).encode()
        h.update(len(data).to_bytes(8, 'big'))
        h.update(data)
    return h.hexdigest()


@dataclass(frozen=True)
class PumpOnlyConvention:
    analysis_axis: AnalysisFrequencyAxis
    number_density_m3: float
    coupling_scales: np.ndarray
    readout_scales: np.ndarray
    source_names: tuple[str, ...]
    carrier_offsets_at_rest_rad_s: np.ndarray
    port_wavevectors_rad_m: np.ndarray
    model_id: str
    model_scope: str
    provenance: str
    phase_provenance: str
    mode_labels: tuple[str, ...] = MODE_LABELS
    readout_units: tuple[str, ...] = ('1',)*4
    drive_units: tuple[str, ...] = ('sqrt(photon/s)',)*4
    phase_charges: tuple[int, ...] = CHARGES
    demodulation: str = 'lab-rf-plus-port-offsets-v1'

    def __post_init__(self):
        if not isinstance(self.analysis_axis, AnalysisFrequencyAxis):
            raise TypeError('common laboratory AnalysisFrequencyAxis required')
        density = float(self.number_density_m3)
        if not np.isfinite(density) or density <= 0:
            raise ValueError('finite positive number_density_m3 required')
        object.__setattr__(self, 'number_density_m3', density)
        for name in ('model_id', 'provenance', 'phase_provenance'):
            _text(getattr(self, name), name)
        if self.model_scope not in ('smooth', 'prescribed_segments', 'analytic_fixture'):
            raise ValueError('declare smooth, prescribed_segments or analytic_fixture model scope')
        for name, expected in (('mode_labels', MODE_LABELS), ('phase_charges', CHARGES),
                               ('readout_units', ('1',)*4), ('drive_units', ('sqrt(photon/s)',)*4)):
            if tuple(getattr(self, name)) != expected:
                raise ValueError(f'incompatible pump-only {name}')
            object.__setattr__(self, name, expected)
        if self.demodulation != 'lab-rf-plus-port-offsets-v1':
            raise ValueError('common lab RF demodulation convention required')
        sources = tuple(self.source_names)
        if not sources or sources[0] != 'atomic_inflow' or len(set(sources)) != len(sources):
            raise ValueError('unique ordered source_names beginning with atomic_inflow required')
        for source in sources:
            _text(source, 'source name')
        object.__setattr__(self, 'source_names', sources)
        for name in ('coupling_scales', 'readout_scales', 'carrier_offsets_at_rest_rad_s', 'port_wavevectors_rad_m'):
            value = readonly_array(getattr(self, name), real=True)
            if value.shape != ((4, 3) if name == 'port_wavevectors_rad_m' else (4,)):
                raise ValueError(f'four Nambu entries required for {name}')
            if name in ('coupling_scales', 'readout_scales'):
                if np.any(value <= 0) or not np.array_equal(value[:2], value[2:]):
                    raise ValueError(f'positive dagger-paired {name} required')
            elif not np.array_equal(value[:2], -value[2:]):
                raise ValueError(f'signed dagger-paired {name} required')
            object.__setattr__(self, name, value)
        if not np.array_equal(self.carrier_offsets_at_rest_rad_s,
                              self.carrier_offsets_at_rest_rad_s[0]*np.array(CHARGES)):
            raise ValueError('rest carrier offsets must belong to one common laboratory beat')

    @property
    def digest(self):
        return _digest([self.analysis_axis.omega_rad_s, self.number_density_m3, self.coupling_scales,
            self.readout_scales, self.source_names, self.carrier_offsets_at_rest_rad_s,
            self.port_wavevectors_rad_m, self.model_id, self.model_scope, self.provenance,
            self.phase_provenance, self.mode_labels, self.readout_units, self.drive_units,
            self.phase_charges, self.demodulation])


@dataclass(frozen=True)
class NumericalComparison:
    """Observed normalized errors and declared budgets, not a passed boolean.

    error_definition must state the norm/scales (including treatment of zero
    signals). IDs bind reference and candidate computations. control_values
    identify actual coarse/fine resolution or distinct independent scrambles.
    The adapter validates this ledger; it cannot verify a fabricated reference.
    """

    kind: str
    errors: Mapping
    tolerances: Mapping
    reference_id: str
    candidate_id: str
    control_values: tuple[float, float]
    error_definition: str
    provenance: str

    def __post_init__(self):
        if self.kind not in ('path_refinement', 'independent_reference', 'ensemble_refinement', 'independent_scramble'):
            raise ValueError('unknown numerical comparison kind')
        for name in ('reference_id', 'candidate_id', 'error_definition', 'provenance'):
            _text(getattr(self, name), name)
        if self.kind in ('independent_reference', 'independent_scramble') and self.reference_id == self.candidate_id:
            raise ValueError('independent controls cannot cite the candidate itself as their reference')
        if not isinstance(self.errors, Mapping) or not isinstance(self.tolerances, Mapping):
            raise TypeError('named numerical error/tolerance mappings required')
        errors = {str(k): float(v) for k, v in self.errors.items()}
        tolerances = {str(k): float(v) for k, v in self.tolerances.items()}
        if not errors or set(errors) != set(tolerances):
            raise ValueError('matching named numerical errors and tolerances required')
        if any(not np.isfinite(v) or v < 0 for v in errors.values()):
            raise ValueError('observed errors must be finite and nonnegative')
        if any(not np.isfinite(v) or v <= 0 for v in tolerances.values()):
            raise ValueError('numerical tolerances must be finite and positive')
        controls = tuple(float(v) for v in self.control_values)
        if len(controls) != 2 or not np.isfinite(controls).all() or controls[0] == controls[1]:
            raise ValueError('two distinct finite numerical controls required')
        if self.kind in ('path_refinement', 'ensemble_refinement') and not 0 < controls[0] < controls[1]:
            raise ValueError('refinement controls must increase from coarse to fine resolution')
        object.__setattr__(self, 'control_values', controls)
        object.__setattr__(self, 'errors', MappingProxyType(errors))
        object.__setattr__(self, 'tolerances', MappingProxyType(tolerances))

    @property
    def passed(self):
        return all(self.errors[name] <= self.tolerances[name] for name in self.errors)


@dataclass(frozen=True)
class ConvergenceEvidence:
    target_digest: str
    model_scope: str
    comparisons: tuple[NumericalComparison, ...]
    provenance: str

    def __post_init__(self):
        if (not isinstance(self.target_digest, str) or len(self.target_digest) != 64 or
                any(c not in '0123456789abcdef' for c in self.target_digest)):
            raise ValueError('target_digest must be a SHA256 computation identifier')
        _text(self.provenance, 'evidence provenance')
        if self.model_scope not in ('smooth', 'prescribed_segments', 'analytic_fixture'):
            raise ValueError('explicit evidence model scope required')
        comparisons = tuple(self.comparisons)
        if not all(isinstance(row, NumericalComparison) for row in comparisons):
            raise TypeError('typed numerical comparisons required')
        object.__setattr__(self, 'comparisons', comparisons)


@dataclass(frozen=True)
class SuppliedPathPacket:
    packet: Mapping
    convention: PumpOnlyConvention
    evidence: ConvergenceEvidence | None = None


def _packet_arrays(path, packet, convention):
    if not isinstance(path, BallisticPath) or not isinstance(convention, PumpOnlyConvention):
        raise TypeError('explicit path and pump-only convention required')
    if not isinstance(packet, Mapping):
        raise TypeError('packet must follow the declared wavepacket mapping schema')
    axis = packet['analysis_axis']
    if not isinstance(axis, AnalysisFrequencyAxis) or not np.array_equal(axis.omega_rad_s, convention.analysis_axis.omega_rad_s):
        raise ValueError('wavepackets must use the same laboratory RF analysis_axis')
    if tuple(packet['source_names']) != convention.source_names:
        raise ValueError('ordered source_names must match the declared reservoirs')
    tau = float(packet['residence_time_s'])
    if not np.isfinite(tau) or not np.isclose(tau, path.residence_time_s, rtol=1e-12, atol=0):
        raise ValueError('packet residence must equal the supplied first-exit path')
    metadata = packet['metadata']
    if 'pump_only' in metadata and metadata['pump_only'] is not True:
        raise ValueError('finite-seed/non-pump-only packets cannot use this phase theorem')
    if 'phase_charges' in metadata and tuple(metadata['phase_charges']) != CHARGES:
        raise ValueError('packet common-phase charges differ from the pump-only theorem')
    if 'demodulation' in metadata and metadata['demodulation'] != convention.demodulation:
        raise ValueError('packet demodulation differs from the common lab RF convention')
    if tuple(metadata['mode_labels']) != convention.mode_labels:
        raise ValueError('packet Nambu mode ordering differs from the convention')
    coupling = readonly_array(metadata['coupling_s_inverse_sqrt_flux'], real=True)
    if coupling.shape == (2,):
        coupling = np.tile(coupling, 2)
    if not np.array_equal(coupling, convention.coupling_scales):
        raise ValueError('packet reciprocal coupling scales differ from the convention')
    _text(metadata['scope'], 'packet model scope provenance')
    if not np.isfinite(float(metadata['entry_phase_rad'])):
        raise ValueError('finite common entry phase required')
    spatial_phase = metadata.get('entry_optical_demodulation_phases_rad')
    if spatial_phase is None and convention.model_scope != 'analytic_fixture':
        raise ValueError('native Rb packet must declare its spatial entry demodulation phases')
    if spatial_phase is not None:
        spatial_phase = readonly_array(spatial_phase, real=True)
        expected_phase = float(metadata['entry_phase_rad'])*np.array(CHARGES[:2])-convention.port_wavevectors_rad_m[:2]@path.entry_position_m
        if spatial_phase.shape != (2,) or np.any(np.abs(np.angle(np.exp(1j*(spatial_phase-expected_phase)))) > 1e-10):
            raise ValueError('packet spatial entry phases differ from the actual path entry position')
    # Some native atomic packets leave density UNUSED. The typed convention
    # supplies it explicitly; a packet that does declare a density must agree.
    if 'number_density_m3' in metadata and float(metadata['number_density_m3']) != convention.number_density_m3:
        raise ValueError('packet density differs from the incoming physical density')
    for key, expected in (('readout_units', convention.readout_units), ('drive_units', convention.drive_units)):
        if key in metadata and tuple(metadata[key]) != expected:
            raise ValueError(f'packet {key} differs from the declared units')
    if 'readout_scales' in metadata and not np.array_equal(metadata['readout_scales'], convention.readout_scales):
        raise ValueError('packet readout scales differ from the convention')
    offsets = readonly_array(metadata['carrier_offsets_rad_s'], real=True)
    expected_offsets = convention.carrier_offsets_at_rest_rad_s-convention.port_wavevectors_rad_m@path.velocity_m_s
    if offsets.shape != (4,) or np.any(np.abs(offsets-expected_offsets) > 5e-14*np.maximum(1, np.abs(expected_offsets))):
        raise ValueError('path Doppler offsets differ from declared signed carrier geometry')
    freq = readonly_array(packet['frequencies_rad_s'], real=True)
    expected = axis.omega_rad_s[:, None]+expected_offsets
    if freq.shape != expected.shape or np.any(np.abs(freq-expected) > 5e-14*np.maximum(1, np.abs(expected))):
        raise ValueError('packet port frequencies must be common lab RF plus the path offsets')
    nf, ns = len(axis.omega_rad_s), len(convention.source_names)
    arrays = {}
    for name in ('greater', 'lesser', 'greater_by_source', 'lesser_by_source', 'mean_pulse', 'retarded_response'):
        value = readonly_array(packet[name])
        shape = (nf, 4) if name == 'mean_pulse' else (ns, nf, 4, 4) if name.endswith('_by_source') else (nf, 4, 4)
        if value.shape != shape:
            raise ValueError(f'{name} has wrong frequency/source/Nambu shape')
        arrays[name] = value
    # A small dark-ordering roundoff floor is stated in s^2 using declared
    # readout scales, not inferred from a missing commutator or added as noise.
    floor = 128*np.finfo(float).eps*tau**2*float(np.max(convention.readout_scales))**2
    for ordering in ('greater', 'lesser'):
        for name in (ordering, ordering+'_by_source'):
            matrix = arrays[name]
            norm = np.max(np.abs(matrix), axis=(-2, -1))
            tolerance = floor+1e-9*norm
            hermitian = matrix.conj().swapaxes(-1, -2)
            if np.any(np.max(np.abs(matrix-hermitian), axis=(-2, -1)) > tolerance):
                raise ValueError(f'{name} is not Hermitian')
            if np.any(np.linalg.eigvalsh((matrix+hermitian)/2)[..., 0] < -tolerance):
                raise ValueError(f'{name} is not positive semidefinite')
        total, summed = arrays[ordering], arrays[ordering+'_by_source'].sum(axis=0)
        closure_scale = np.maximum(np.max(np.abs(total), axis=(-2, -1)), np.max(np.abs(summed), axis=(-2, -1)))
        if np.any(np.max(np.abs(total-summed), axis=(-2, -1)) > floor+1e-10*closure_scale):
            raise ValueError(f'{ordering} does not equal the complete by_source sum')
    return arrays, freq


def _validated_packet_digest(path, packet, convention, arrays, freq):
    return _digest([convention.digest, path.entry_position_m, path.velocity_m_s, path.residence_time_s,
        path.source, freq, float(packet['metadata']['entry_phase_rad']),
        str(packet['metadata']['scope']), *[arrays[name] for name in sorted(arrays)]])


def packet_digest(path, packet, convention):
    """Bind all consumed numerical outputs, frame, path and model declaration."""
    arrays, freq = _packet_arrays(path, packet, convention)
    return _validated_packet_digest(path, packet, convention, arrays, freq)


def _evidence_reasons(evidence, target, scope, *, ensemble=False):
    if evidence is None:
        return ['explicit ensemble convergence evidence missing' if ensemble else 'explicit path convergence evidence missing']
    if not isinstance(evidence, ConvergenceEvidence):
        return ['typed numerical convergence evidence required; audit booleans do not certify convergence']
    if evidence.target_digest != target or evidence.model_scope != scope:
        return ['convergence evidence does not match this computation digest/model scope']
    required = STREAM_METRICS if ensemble else PATH_METRICS
    kinds = ('ensemble_refinement', 'independent_scramble') if ensemble else ('path_refinement', 'independent_reference')
    reasons = []
    for kind in kinds:
        rows = [row for row in evidence.comparisons if row.kind == kind]
        needed = 2 if not ensemble and scope == 'smooth' and kind == 'path_refinement' else 1
        if len(rows) < needed:
            reasons.append(f'{kind}: {needed} explicit comparison(s) required')
            continue
        rows = rows[-needed:]
        if rows[-1].candidate_id != target:
            reasons.append(f'{kind}: final comparison is not bound to the actual candidate')
        if needed == 2 and (rows[0].control_values[1] != rows[1].control_values[0]
                or rows[0].candidate_id != rows[1].reference_id):
            reasons.append('smooth path refinement must link two successive comparisons at three resolutions')
        for row in rows:
            if not set(required) <= set(row.errors):
                reasons.append(f'{kind}: required output metrics are missing')
            if not row.passed:
                reasons.append(f'{kind}: observed numerical error exceeds its declared tolerance')
    return reasons


def _spectrum_digest(result):
    spectra = result['spectra']
    return _digest([result['convention_digest'], result['inflow_digest'], result['packet_chain_digest'],
        result['analysis_axis'].omega_rad_s, result['source_names'], result['model_scope'], result['model_id'],
        result['path_count'], result['consumed_path_count'], result['path_evidence_passed'],
        result['total_arrival_rate_s_inverse'], result['mean_occupancy'], result['equilibrium_nV'],
        result['provenance'], result['phase_provenance'], result['density_ledger'],
        result['scope'], result['numerical_units'],
        *[spectra[name] for name in sorted(spectra)]])


def aggregate_pump_only_stream(inflow, convention, packet_factory):
    """Consume exactly one declared packet per path; never re-solve by phase.

    A failing atomic audit/evidence returns spectra=None with the failing
    path/rate and consumed count. No incomplete stream is renormalized. Valid
    paths produce an explicitly uncertified diagnostic candidate until actual
    ensemble refinement and independent-scramble evidence is attached.
    """
    if not isinstance(inflow, ThermalInflowQuadrature) or not isinstance(convention, PumpOnlyConvention):
        raise TypeError('thermal inflow quadrature and pump-only convention required')
    if inflow.density_m3 != convention.number_density_m3:
        raise ValueError('inflow density must exactly match the consumed model density')
    if not callable(packet_factory):
        raise TypeError('per-path packet factory required')
    inflow_id = _digest([inflow.lower_corner_m, inflow.upper_corner_m, inflow.temperature_K,
        inflow.mass_kg, inflow.density_m3, inflow.entry_position_m, inflow.velocity_m_s,
        inflow.residence_time_s, inflow.rate_s_inverse, inflow.face_index, inflow.source])
    result = {'certified': False, 'path_evidence_passed': False, 'spectra': None,
        'reasons': (), 'candidate_digest': None, 'convention_digest': convention.digest,
        'model_scope': convention.model_scope, 'model_id': convention.model_id,
        'analysis_axis': convention.analysis_axis, 'source_names': convention.source_names,
        'inflow_digest': inflow_id, 'path_count': len(inflow.rate_s_inverse), 'consumed_path_count': 0,
        'total_arrival_rate_s_inverse': inflow.total_arrival_rate_s_inverse,
        'mean_occupancy': inflow.mean_occupancy, 'equilibrium_nV': inflow.equilibrium_atom_number,
        'provenance': convention.provenance, 'phase_provenance': convention.phase_provenance,
        'density_ledger': 'density is included once in boundary rates; no additional density multiplier',
        'scope': 'conditional pump-only zero-cyclic atomic stream; no optical Maxwell channel, SQL or experimental certification'}
    accumulators = None
    chain = hashlib.sha256()
    for index, rate in enumerate(inflow.rate_s_inverse):
        path = inflow.path(index)
        supplied = packet_factory(path)
        result['consumed_path_count'] = index+1
        if not isinstance(supplied, SuppliedPathPacket) or not isinstance(supplied.convention, PumpOnlyConvention):
            raise TypeError('factory must supply a typed SuppliedPathPacket')
        if supplied.convention.digest != convention.digest:
            raise ValueError('packet units/scales/model/provenance/density convention differs from the stream')
        arrays, freq = _packet_arrays(path, supplied.packet, convention)
        token = _validated_packet_digest(path, supplied.packet, convention, arrays, freq)
        reasons = _evidence_reasons(supplied.evidence, token, convention.model_scope)
        if supplied.packet.get('audit', {}).get('passed') is not True:
            reasons.insert(0, 'packet internal physical/numerical audit did not pass')
        if reasons:
            result.update(reasons=tuple(reasons), failed_path_index=index, failed_path_source=path.source,
                failed_path_rate_s_inverse=float(rate), packet_chain_digest=chain.hexdigest())
            return result
        chain.update(bytes.fromhex(token))
        # Pump-only H and boundary rho are entry-phase independent. Readouts
        # carry s=(1,-1,-1,1), reciprocal drives carry -s. Hence C and response
        # transform as D(phi) X D(phi)^dagger; uniform phase integrates to
        # delta(s_j,s_k). Keep mask*(m m†), even though mean phase average=0.
        # This theorem does not apply to finite-seed saturation.
        mean = arrays['mean_pulse']
        number = mean[:, :, None]*mean[:, None, :].conj()*MASK
        terms = {name: arrays[name]*MASK for name in ('greater_by_source', 'lesser_by_source', 'retarded_response')}
        terms['poisson_number'] = number
        if accumulators is None:
            accumulators = {name: np.zeros_like(value) for name, value in terms.items()}
        for name, value in terms.items():
            accumulators[name] += rate*value
    # Totals are assembled from the validated complete source ledger. This
    # keeps source closure visible without counting atomic inflow twice.
    for ordering in ('greater', 'lesser'):
        internal = accumulators[ordering+'_by_source'].sum(axis=0)
        accumulators['internal_'+ordering] = internal
        accumulators[ordering] = internal+accumulators['poisson_number']
    result.update(path_evidence_passed=True, packet_chain_digest=chain.hexdigest(),
        spectra=MappingProxyType({name: readonly_array(value) for name, value in accumulators.items()}),
        reasons=('actual stream quadrature convergence has not been supplied',),
        numerical_units={'covariance': 's', 'retarded_response': 'sqrt(s)'})
    result['candidate_digest'] = _spectrum_digest(result)
    return result


def certify_pump_only_stream(candidate, ensemble_evidence):
    """Attach actual-integrand ensemble evidence without re-evaluating paths.

    The evidence must cover raw ordered spectra, source terms, Poisson number
    and response for this exact candidate. nV or Maxwell moments alone do not
    satisfy it. Certification is limited to the declared model and budgets.
    """
    result = dict(candidate)
    if not result['path_evidence_passed'] or result['spectra'] is None:
        result['certified'] = False
        return result
    current = _spectrum_digest(result)
    if current != result['candidate_digest']:
        raise ValueError('candidate outputs or computation identity changed after accumulation')
    reasons = _evidence_reasons(ensemble_evidence, current, result['model_scope'], ensemble=True)
    result.update(certified=not reasons, reasons=tuple(reasons), ensemble_evidence=ensemble_evidence)
    return result
