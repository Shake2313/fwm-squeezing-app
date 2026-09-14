"""Immutable audit of a conditional pump-only finite-atom stream adapter.

python -m analysis.grand_challenge.transport_ensemble_audit --output NEW.json --plot NEW.png

The numerical fixture is an imposed analytic phase-mark family of a constant
two-level atom. Its zero carrier offsets do not model a nonzero Rb beat. A
separate actual Rb packet checks schema compatibility and must remain rejected
without path and ensemble convergence evidence.
"""

import argparse
from collections.abc import Mapping
from dataclasses import replace
from functools import lru_cache
import hashlib
import json
from pathlib import Path
import platform
import time

import numpy as np
import scipy
from scipy.special import roots_legendre

from gabes import core, constants as c
from gabes.fwm_quantum import transport_ensemble as stream
from gabes.quantum.contracts import AnalysisFrequencyAxis, GeneratorFrequencyAxis
from gabes.quantum.inflow import maxwell_box_inflow
from gabes.quantum.reservoirs import ExplicitReservoirs
from gabes.quantum.transport import integrate_characteristic


ROOT = Path(__file__).resolve().parents[2]
POWERS = (4, 6, 8, 10)
SEEDS = (11, 211)
GAUSS_ORDERS = (24, 48)
AXIS = AnalysisFrequencyAxis.from_hz([0., 1e4, 3e4])
LENGTHS_M = np.array([4e-4, 2e-4, 2e-4])
TEMPERATURE_K, DENSITY_M3 = 373., 3e16
PATH_BUDGET, STREAM_BUDGET = 1e-9, .05
PHASE_BUDGET, TRANSPORT_BUDGET = 2e-12, 2e-8
# Independent literals: this reference never imports the production charge mask.
CHARGES = np.array([1, -1, -1, 1])
PATH_KEYS = ('greater', 'lesser', 'greater_by_source', 'lesser_by_source',
             'mean_pulse', 'mean_outer', 'retarded_response')
STREAM_KEYS = ('greater', 'lesser', 'greater_by_source', 'lesser_by_source',
               'poisson_number', 'retarded_response')
ERROR_DEFINITION = (
    'Maximum across RF rows of the Frobenius difference divided by the larger '
    'of the reference-row norm and 1e-12 times the explicitly stated unit-matched '
    'comparison scale. Source blocks share one Frobenius norm at each RF. '
    'No covariance, rate or numerical discrepancy is clipped or fitted.')


def source_hashes():
    paths = set(ROOT.glob('gabes/quantum/*.py'))
    paths.update(ROOT.glob('gabes/fwm_quantum/*.py'))
    paths.update(ROOT.glob('tests/quantum/test_*transport*.py'))
    paths.update(ROOT.glob('analysis/grand_challenge/reference/*.py'))
    paths.update(ROOT/name for name in (
        'gabes/__init__.py', 'gabes/core.py', 'gabes/constants.py', 'gabes/atoms.py',
        'gabes/observables.py', 'gabes/doppler.py', 'gabes/hyperfine.py',
        'gabes/species.py', 'gabes/schemes/fwm.py',
        'analysis/grand_challenge/normalization_audit.py',
        'analysis/grand_challenge/transport_ensemble_audit.py', 'tests/quantum/test_inflow.py'))
    return {p.relative_to(ROOT).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in sorted(paths)}


def json_value(value):
    if isinstance(value, np.ndarray):
        return ({'real': value.real.tolist(), 'imag': value.imag.tolist()}
                if np.iscomplexobj(value) else value.tolist())
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, Mapping):
        return {str(key): json_value(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [json_value(item) for item in value]
    return value


def toy_atom():
    """H=0, no jumps: all finite-time moments follow from this density matrix."""
    rho = np.array([[.6, .18+.07j], [.18-.07j, .4]])
    lower = np.array([[0., 1.], [0., 0.]])
    z, identity = np.diag([1., -1.]), np.eye(2)
    two = np.array([lower+.25*z+.2*identity,
                    .7*lower.T-.15*z+(.1+.2j)*identity])
    operators = np.concatenate([two, two.conj().swapaxes(-1, -2)])
    coupling = np.array([2., 3., 2., 3.])
    drives = coupling[:, None, None]*operators.conj().swapaxes(-1, -2)
    mean = np.array([np.trace(o@rho) for o in operators])
    greater = np.array([[np.trace(a@b.conj().T@rho) for b in operators]
                        for a in operators])-np.outer(mean, mean.conj())
    lesser = np.array([[np.trace(b.conj().T@a@rho) for b in operators]
                       for a in operators])-np.outer(mean, mean.conj())
    response = np.array([[np.trace(o@(-1j*(v@rho-rho@v))) for v in drives]
                         for o in operators])
    scales = np.tile([np.linalg.norm(o, 2) for o in two], 2)
    return {'rho': rho, 'operators': operators, 'drives': drives, 'coupling': coupling,
            'scales': scales, 'mean': mean, 'greater': greater, 'lesser': lesser,
            'response': response}


def make_inflow(power, seed, *, lengths=LENGTHS_M, density=DENSITY_M3):
    return maxwell_box_inflow(-lengths/2, lengths/2, temperature_K=TEMPERATURE_K,
        mass_kg=c.MASS_85RB, density_m3=density, points_per_face_power=power, seed=seed,
        source='Assumed toy-gas Maxwell boundary reservoir; declared Rb85 mass only, no Rb optical prediction')


def toy_convention(inflow, atom):
    return stream.PumpOnlyConvention(analysis_axis=AXIS, number_density_m3=inflow.density_m3,
        coupling_scales=atom['coupling'], readout_scales=atom['scales'],
        source_names=('atomic_inflow',), carrier_offsets_at_rest_rad_s=np.zeros(4),
        port_wavevectors_rad_m=np.zeros((4, 3)), model_id='constant-two-level-phase-mark-family-v1',
        model_scope='analytic_fixture',
        provenance='H=0, no jumps, explicit coherent mixed state, dimensionless dagger-paired readouts',
        phase_provenance='Imposed uniform common phase mark D(phi)=diag(exp(i*(1,-1,-1,1)*phi)); zero optical offsets, not a nonzero Rb beat')


@lru_cache(maxsize=4)
def gauss_rule(order):
    nodes, weights = roots_legendre(order)
    return (nodes+1)/2, weights/2


def finite_pulse_integrals(tau, *, order=None):
    """F=int exp(i*w*a) da; T=int (tau-u)*exp(i*w*u) du.

    T is the finite-age retarded integral only because every input and output
    port in this fixture has the same w. No unequal-port formula is implied.
    """
    z = 1j*AXIS.omega_rad_s*tau
    if order is not None:
        x, weights = gauss_rule(order)
        oscillation = np.exp(z[:, None]*x)
        return tau*(oscillation@weights), tau**2*(oscillation@(weights*(1-x)))
    f = tau*np.sinc(AXIS.omega_rad_s*tau/(2*np.pi))*np.exp(z/2)
    response = np.empty_like(z)
    small = np.abs(z) < .1
    # Stable analytic Taylor representation of (exp(z)-1-z)/z**2, including DC.
    term = np.full(np.count_nonzero(small), .5, dtype=complex)
    series = term.copy()
    for power in range(1, 16):
        term = term*z[small]/(power+2)
        series += term
    response[small] = series
    response[~small] = (np.expm1(z[~small])-z[~small])/z[~small]**2
    return f, tau**2*response


def toy_packet(path, convention, atom, *, order=None):
    f, response = finite_pulse_integrals(path.residence_time_s, order=order)
    greater = np.abs(f[:, None, None])**2*atom['greater']
    lesser = np.abs(f[:, None, None])**2*atom['lesser']
    return {'analysis_axis': AXIS,
        'frequencies_rad_s': np.broadcast_to(AXIS.omega_rad_s[:, None], (len(f), 4)),
        'source_names': ('atomic_inflow',), 'greater': greater, 'lesser': lesser,
        'greater_by_source': greater[None], 'lesser_by_source': lesser[None],
        'mean_pulse': f[:, None]*atom['mean'],
        'retarded_response': response[:, None, None]*atom['response'],
        'residence_time_s': path.residence_time_s,
        'audit': {'passed': bool(np.linalg.eigvalsh(atom['rho']).min() >= 0
            and min(np.linalg.eigvalsh(atom[key]).min() for key in ('greater', 'lesser')) > -1e-12)},
        'metadata': {'mode_labels': convention.mode_labels,
            'coupling_s_inverse_sqrt_flux': convention.coupling_scales,
            'carrier_offsets_rad_s': np.zeros(4), 'entry_phase_rad': 0.,
            'entry_optical_demodulation_phases_rad': np.zeros(2),
            'scope': 'Imposed analytic phase-mark family of a constant two-level atom; zero port offsets; not actual periodic Rb'}}


def metric_arrays(packet):
    result = {key: packet[key] for key in PATH_KEYS if key != 'mean_outer'}
    mean = packet['mean_pulse']
    result['mean_outer'] = mean[:, :, None]*mean[:, None, :].conj()
    return result


def metric_scales(time_scale, atom, *, rate=1.):
    scale = float(np.max(atom['scales']))
    covariance = 4*rate*time_scale**2*scale**2
    result = {key: covariance for key in set(PATH_KEYS+STREAM_KEYS)}
    result['mean_pulse'] = 2*rate*time_scale*scale
    result['retarded_response'] = covariance*float(np.max(atom['coupling']))
    return result


def observed_errors(candidate, reference, keys, scales):
    errors = {}
    for key in keys:
        a, b = np.asarray(candidate[key]), np.asarray(reference[key])
        if key.endswith('_by_source'):
            a, b = np.moveaxis(a, 1, 0), np.moveaxis(b, 1, 0)
        delta = np.linalg.norm((a-b).reshape(len(AXIS.omega_rad_s), -1), axis=1)
        norm = np.linalg.norm(b.reshape(len(AXIS.omega_rad_s), -1), axis=1)
        errors[key] = float(np.max(delta/np.maximum(norm, 1e-12*scales[key])))
    return errors


def measured_comparison(kind, errors, budget, reference_id, candidate_id, controls, provenance):
    return stream.NumericalComparison(kind=kind, errors=errors,
        tolerances={key: budget for key in errors}, reference_id=reference_id,
        candidate_id=candidate_id, control_values=controls,
        error_definition=ERROR_DEFINITION, provenance=provenance)


def comparison_record(row):
    return {'kind': row.kind, 'errors': dict(row.errors), 'tolerances': dict(row.tolerances),
        'reference_id': row.reference_id, 'candidate_id': row.candidate_id,
        'control_values': list(row.control_values), 'error_definition': row.error_definition,
        'provenance': row.provenance, 'passed': row.passed}


def measured_path_packet(path, convention, atom):
    coarse = toy_packet(path, convention, atom, order=GAUSS_ORDERS[0])
    fine = toy_packet(path, convention, atom, order=GAUSS_ORDERS[1])
    exact = toy_packet(path, convention, atom)
    coarse_id, target, exact_id = [stream.packet_digest(path, item, convention)
                                  for item in (coarse, fine, exact)]
    scales = metric_scales(path.residence_time_s, atom)
    refine = observed_errors(metric_arrays(fine), metric_arrays(coarse), PATH_KEYS, scales)
    analytic = observed_errors(metric_arrays(fine), metric_arrays(exact), PATH_KEYS, scales)
    rows = (measured_comparison('path_refinement', refine, PATH_BUDGET, coarse_id, target,
                GAUSS_ORDERS, 'Actual Gauss-Legendre age/causal-triangle integrals on this same first-exit path'),
            measured_comparison('independent_reference', analytic, PATH_BUDGET,
                'closed-form:'+exact_id, target, (0., GAUSS_ORDERS[1]),
                'Independent finite-pulse sinc/exponential and analytic causal-triangle formula, including DC'))
    evidence = stream.ConvergenceEvidence(target, convention.model_scope, rows,
        'Observed errors from three actual calculations; no zero-error placeholders')
    return stream.SuppliedPathPacket(fine, convention, evidence), exact


def literal_four_phase_average(packet):
    """Independent four-mark sum; preserve average of each mean outer product."""
    names = ('greater', 'lesser', 'greater_by_source', 'lesser_by_source', 'retarded_response')
    result = {name: np.zeros_like(packet[name]) for name in names}
    result['poisson_number'] = np.zeros_like(packet['greater'])
    result['mean_pulse'] = np.zeros_like(packet['mean_pulse'])
    for phase in .173+np.arange(4)*np.pi/2:
        diagonal = np.diag(np.exp(1j*CHARGES*phase))
        for name in names:
            result[name] += (diagonal@packet[name]@diagonal.conj().T)/4
        mean = packet['mean_pulse']@diagonal.T
        result['mean_pulse'] += mean/4
        result['poisson_number'] += mean[:, :, None]*mean[:, None, :].conj()/4
    for name in ('greater', 'lesser'):
        result[name] += result['poisson_number']
    return result


def run_grid(power, seed, atom):
    inflow = make_inflow(power, seed)
    convention = toy_convention(inflow, atom)
    count = 0
    reference = None
    worst = {kind: {key: 0. for key in PATH_KEYS}
             for kind in ('path_refinement', 'independent_reference')}
    worst_record, worst_ratio = None, -1.
    evidence_chain = hashlib.sha256()

    def factory(path):
        nonlocal count, reference, worst_record, worst_ratio
        index = count
        count += 1
        supplied, exact = measured_path_packet(path, convention, atom)
        rows = [comparison_record(row) for row in supplied.evidence.comparisons]
        evidence_chain.update(json.dumps(rows, sort_keys=True, allow_nan=False).encode())
        ratio = max(error/PATH_BUDGET for row in rows for error in row['errors'].values())
        if ratio > worst_ratio:
            worst_ratio, worst_record = ratio, {'path_index': index, 'residence_time_s': path.residence_time_s,
                                               'target_digest': supplied.evidence.target_digest, 'comparisons': rows}
        for row in rows:
            for key, value in row['errors'].items():
                worst[row['kind']][key] = max(worst[row['kind']][key], value)
        phase_average = literal_four_phase_average(exact)
        if reference is None:
            reference = {name: np.zeros_like(value) for name, value in phase_average.items()}
        for name, value in phase_average.items():
            reference[name] += inflow.rate_s_inverse[index]*value
        return supplied

    candidate = stream.aggregate_pump_only_stream(inflow, convention, factory)
    characteristic_time = float(np.max(LENGTHS_M)/np.sqrt(c.KB*TEMPERATURE_K/c.MASS_85RB))
    scales = metric_scales(characteristic_time, atom, rate=inflow.total_arrival_rate_s_inverse)
    phase_errors = (observed_errors(candidate['spectra'], reference, STREAM_KEYS, scales)
                    if candidate['spectra'] is not None else None)
    row = {'power': power, 'seed': seed, 'path_count': len(inflow.rate_s_inverse),
        'callback_count': count, 'consumed_path_count': candidate['consumed_path_count'],
        'callback_once_per_consumed_path': count == candidate['consumed_path_count'],
        'all_paths_consumed': count == len(inflow.rate_s_inverse),
        'path_evidence_passed': candidate['path_evidence_passed'],
        'reasons': list(candidate['reasons']), 'candidate_digest': candidate['candidate_digest'],
        'packet_chain_digest': candidate['packet_chain_digest'],
        'path_evidence_chain_sha256': evidence_chain.hexdigest(),
        'max_path_errors': worst, 'worst_path_evidence': worst_record,
        'literal_four_phase_errors': phase_errors,
        'mean_occupancy': inflow.mean_occupancy, 'equilibrium_nV': inflow.equilibrium_atom_number,
        'occupancy_over_nV': inflow.mean_occupancy/inflow.equilibrium_atom_number,
        'total_arrival_rate_s_inverse': inflow.total_arrival_rate_s_inverse,
        'residence_time_range_s': [float(inflow.residence_time_s.min()), float(inflow.residence_time_s.max())],
        'comparison_scales': scales}
    return candidate, row, reference


def stream_comparisons(candidate, previous, independent, atom, previous_power, power):
    characteristic_time = float(np.max(LENGTHS_M)/np.sqrt(c.KB*TEMPERATURE_K/c.MASS_85RB))
    scales = metric_scales(characteristic_time, atom, rate=candidate['total_arrival_rate_s_inverse'])
    comparisons = []
    for kind, reference, controls, text in (
        ('ensemble_refinement', previous, (2**previous_power, 2**power),
         'Recomputed complete actual toy integrands on two physical boundary Sobol grids, same seed'),
        ('independent_scramble', independent, (SEEDS[1], SEEDS[0]),
         'Recomputed complete actual toy integrands on an independent scrambled boundary Sobol grid')):
        errors = observed_errors(candidate['spectra'], reference['spectra'], STREAM_KEYS, scales)
        comparisons.append(measured_comparison(kind, errors, STREAM_BUDGET,
            reference['candidate_digest'], candidate['candidate_digest'], controls, text))
    return stream.ConvergenceEvidence(candidate['candidate_digest'], candidate['model_scope'],
        tuple(comparisons), 'Actual raw/source/Poisson/response output convergence; phase-space moments are diagnostic only')


def evidence_negative_controls(atom, candidate, reference):
    inflow = make_inflow(0, SEEDS[0])
    convention = toy_convention(inflow, atom)
    audit_only = stream.aggregate_pump_only_stream(inflow, convention,
        lambda path: stream.SuppliedPathPacket(toy_packet(path, convention, atom), convention))
    callbacks = 0

    def one_bad(path):
        nonlocal callbacks
        supplied, _ = measured_path_packet(path, convention, atom)
        index, callbacks = callbacks, callbacks+1
        if index == 2:
            rows = list(supplied.evidence.comparisons)
            errors = dict(rows[0].errors)
            errors['mean_outer'] = 2*PATH_BUDGET
            rows[0] = replace(rows[0], errors=errors,
                provenance='Deliberately failed synthetic evidence gate, not a measured convergence claim')
            supplied = replace(supplied, evidence=replace(supplied.evidence, comparisons=tuple(rows)))
        return supplied

    failed = stream.aggregate_pump_only_stream(inflow, convention, one_bad)
    moment_error = abs(candidate['mean_occupancy']/candidate['equilibrium_nV']-1)
    rows = tuple(measured_comparison(kind, {'occupancy': moment_error}, 1.,
        'analytic nV / moment-only deliberately incomplete control', candidate['candidate_digest'],
        (1., 2.), 'Actual occupancy error, deliberately lacking all actual spectrum/response metrics')
        for kind in ('ensemble_refinement', 'independent_scramble'))
    moment_only = stream.certify_pump_only_stream(candidate, stream.ConvergenceEvidence(
        candidate['candidate_digest'], candidate['model_scope'], rows,
        'Deliberate moment-only control; cannot certify the nonlinear packet integral'))
    absent = stream.certify_pump_only_stream(candidate, None)
    spectra = candidate['spectra']
    norm_number = np.linalg.norm(spectra['poisson_number'])
    erased_error = float(np.linalg.norm(spectra['internal_greater']-spectra['greater'])
                         /np.linalg.norm(spectra['greater']))
    rate = candidate['total_arrival_rate_s_inverse']
    wrong_mean_outer = reference['mean_pulse'][:, :, None]*reference['mean_pulse'][:, None, :].conj()/rate
    mean_first_error = float(np.linalg.norm(wrong_mean_outer-spectra['poisson_number'])/norm_number)
    return {'audit_true_without_path_evidence': {
                'native_audit_boolean': True, 'certified': audit_only['certified'],
                'spectra_is_none': audit_only['spectra'] is None,
                'consumed_path_count': audit_only['consumed_path_count'], 'reasons': list(audit_only['reasons'])},
        'one_failed_path': {'certified': failed['certified'], 'spectra_is_none': failed['spectra'] is None,
            'failed_path_index': failed['failed_path_index'], 'callback_count': callbacks,
            'failed_path_rate_s_inverse': failed['failed_path_rate_s_inverse'],
            'declared_total_rate_unchanged': failed['total_arrival_rate_s_inverse'] == inflow.total_arrival_rate_s_inverse,
            'reasons': list(failed['reasons']), 'failure_injected_for_gate_test': True},
        'path_evidence_without_ensemble_evidence': {'certified': absent['certified'], 'reasons': list(absent['reasons'])},
        'moment_only_evidence': {'actual_occupancy_error': moment_error, 'certified': moment_only['certified'],
            'reasons': list(moment_only['reasons'])},
        'delete_mean_outer_relative_raw_greater_error': erased_error,
        'average_mean_before_outer_relative_number_error': mean_first_error,
        'passed': bool(not audit_only['certified'] and audit_only['spectra'] is None
            and audit_only['consumed_path_count'] == 1 and not failed['certified'] and failed['spectra'] is None
            and failed['failed_path_index'] == 2 and callbacks == 3
            and failed['total_arrival_rate_s_inverse'] == inflow.total_arrival_rate_s_inverse
            and not moment_only['certified'] and not absent['certified']
            and erased_error > .01 and mean_first_error > .99)}


def actual_constant_transport_control(atom):
    """Cross-check all seven packet outputs against the actual ODE transport."""
    inflow = make_inflow(0, SEEDS[0])
    convention = toy_convention(inflow, atom)
    path = inflow.path(int(np.argmax(inflow.residence_time_s)))
    exact = toy_packet(path, convention, atom)
    characteristic = integrate_characteristic(path, lambda age, position: np.zeros((2, 2)),
        ExplicitReservoirs(2, ()), atom['rho'], times_s=[0., path.residence_time_s],
        rtol=2e-11, atol=2e-13)
    axis = GeneratorFrequencyAxis(AXIS.omega_rad_s, 'Equal lab/age frequencies of the zero-offset constant atom control')
    packet = characteristic.wavepacket(axis, lambda age, position: atom['operators'],
        drive_operators=lambda age, position: atom['drives'], rtol=2e-11, atol=2e-13)
    errors = observed_errors(metric_arrays(packet), metric_arrays(exact), PATH_KEYS,
                             metric_scales(path.residence_time_s, atom))
    return {'scope': 'Actual generic two-level transport integration; no optical Maxwell propagation',
        'residence_time_s': path.residence_time_s, 'ode_evaluations': packet['ode_evaluations'],
        'errors': errors, 'tolerance': TRANSPORT_BUDGET,
        'state_audit_passed': characteristic.audit['passed'], 'wavepacket_audit_passed': packet['audit']['passed'],
        'passed': bool(characteristic.audit['passed'] and packet['audit']['passed']
                       and max(errors.values()) < TRANSPORT_BUDGET)}


def actual_rb_missing_evidence_control():
    """Native Rb schema compatibility is deliberately insufficient to certify."""
    from gabes.fwm_quantum.kinetic import CarrierGeometry
    from gabes.fwm_quantum.transport import reduced_ballistic_problem, reduced_ballistic_wavepacket
    from .normalization_audit import conditional_inputs
    inputs = replace(conditional_inputs(), transit_rate_s_inverse=0.)
    geometry = CarrierGeometry.vacuum_beams(inputs, probe_angle_rad=.006, conjugate_angle_rad=-.005)
    inflow = make_inflow(0, SEEDS[0], lengths=np.full(3, 2e-4), density=inputs.number_density_m3)
    axis = AnalysisFrequencyAxis.from_hz([1e6])
    path = inflow.path(0)
    problem = reduced_ballistic_problem(inputs, geometry, path, axis, segments=1)
    q = geometry.wavevectors_rad_m[1:]-geometry.wavevectors_rad_m[0]
    coupling = np.tile(problem['metadata']['coupling_s_inverse_sqrt_flux'], 2)
    convention = stream.PumpOnlyConvention(analysis_axis=axis, number_density_m3=inputs.number_density_m3,
        coupling_scales=coupling, readout_scales=np.ones(4),
        source_names=('atomic_inflow',)+tuple('jump:'+channel.name for channel in problem['reservoirs'].channels),
        carrier_offsets_at_rest_rad_s=(-c.OMEGA_HF+inputs.two_photon_rad_s)*CHARGES,
        port_wavevectors_rad_m=np.concatenate([q, -q]), model_id='native-single-frozen-Rb-packet-diagnostic-v1',
        model_scope='prescribed_segments', provenance='One actual reduced RMS Rb atom at one prescribed constant pump segment',
        phase_provenance='One signed lab beat phase; actual nonclosed vacuum wavevectors and entry spatial phases retained')
    actual = {}

    def factory(item):
        packet = reduced_ballistic_wavepacket(inputs, geometry, item, axis, segments=1)
        actual.update(audit=packet['audit'], metadata=packet['metadata'])
        return stream.SuppliedPathPacket(packet, convention)

    rejected = stream.aggregate_pump_only_stream(inflow, convention, factory)
    return {'scope': 'Actual native schema/phase/Doppler compatibility check at one frozen segment; no Rb ensemble certification',
        'native_audit_passed': actual['audit']['passed'], 'native_metadata': json_value(actual['metadata']),
        'convention_digest': convention.digest, 'certified': rejected['certified'],
        'spectra_is_none': rejected['spectra'] is None, 'consumed_path_count': rejected['consumed_path_count'],
        'path_entry_position_m': path.entry_position_m.tolist(), 'path_velocity_m_s': path.velocity_m_s.tolist(),
        'residence_time_s': path.residence_time_s, 'frequency_hz': axis.frequency_hz.tolist(),
        'wavevectors_rad_m': geometry.wavevectors_rad_m.tolist(),
        'number_density_m3': inputs.number_density_m3, 'reasons': list(rejected['reasons']),
        'passed': bool(actual['audit']['passed'] and not rejected['certified'] and rejected['spectra'] is None
                       and rejected['consumed_path_count'] == 1
                       and 'explicit path convergence evidence missing' in rejected['reasons'])}


def build_report():
    started = time.monotonic()
    before = source_hashes()
    atom = toy_atom()
    rows, convergence = [], []
    previous, previous_power = None, None
    certified, final_reference = None, None
    with core.blas_single_thread():
        for power in POWERS:
            nominal, nominal_row, reference = run_grid(power, SEEDS[0], atom)
            independent, independent_row, _ = run_grid(power, SEEDS[1], atom)
            rows.extend([nominal_row, independent_row])
            print(f'Toy boundary grid: 6 x 2**{power}, both scrambles recomputed; path evidence '
                  f'{nominal["path_evidence_passed"]}/{independent["path_evidence_passed"]}', flush=True)
            if not nominal['path_evidence_passed'] or not independent['path_evidence_passed']:
                certified, final_reference = nominal, reference
                break
            if previous is not None:
                proof = stream_comparisons(nominal, previous, independent, atom, previous_power, power)
                certified = stream.certify_pump_only_stream(nominal, proof)
                convergence.append({'power': power, 'certified': certified['certified'],
                    'comparisons': [comparison_record(row) for row in proof.comparisons],
                    'reasons': list(certified['reasons'])})
                print(f'Actual output comparison: maximum error '
                      f'{max(error for row in proof.comparisons for error in row.errors.values()):.5g}; '
                      f'5% fixture gate {certified["certified"]}', flush=True)
                if power >= 8 and certified['certified']:
                    final_reference = reference
                    break
            previous, previous_power, final_reference = nominal, power, reference
        constant = actual_constant_transport_control(atom)
        native = actual_rb_missing_evidence_control()
        negatives = (evidence_negative_controls(atom, certified, final_reference)
                     if certified['spectra'] is not None else {'passed': False, 'reason': 'No complete candidate spectrum'})
    after = source_hashes()
    phase_passed = all(row['literal_four_phase_errors'] is not None
        and max(row['literal_four_phase_errors'].values()) < PHASE_BUDGET for row in rows)
    callback_passed = all(row['all_paths_consumed'] and row['callback_once_per_consumed_path'] for row in rows)
    implementation_passed = bool(before == after and phase_passed and callback_passed
        and all(row['path_evidence_passed'] for row in rows)
        and negatives['passed'] and constant['passed'] and native['passed'])
    return {'schema': 'gabes-pump-only-transport-ensemble-audit-v1',
        'scope': 'Conditional numerical audit of an imposed phase-mark family of a constant two-level atom; no thermal Rb optical prediction',
        'toy_model': {'hamiltonian': 'zero 2 x 2', 'jump_reservoirs': [],
            'boundary_state': json_value(atom['rho']), 'readout_operators': json_value(atom['operators']),
            'reciprocal_drive_operators': json_value(atom['drives']),
            'coupling_scales_s_inverse_sqrt_flux': atom['coupling'].tolist(),
            'readout_scales': atom['scales'].tolist(), 'mode_labels': list(stream.MODE_LABELS),
            'phase_charges': CHARGES.tolist(), 'carrier_offsets_rad_s': [0.]*4,
            'port_wavevectors_rad_m': np.zeros((4, 3)).tolist(),
            'phase_marks_rad': (.173+np.arange(4)*np.pi/2).tolist(),
            'minimum_boundary_state_eigenvalue': float(np.linalg.eigvalsh(atom['rho']).min()),
            'source_names': ['atomic_inflow']},
        'assumed_boundary_inputs': {'lengths_m': LENGTHS_M.tolist(), 'temperature_K': TEMPERATURE_K,
            'mass_kg': c.MASS_85RB, 'density_m3': DENSITY_M3,
            'input_evidence': 'All geometry/state/gas/readout values assumed for this numerical fixture; no independently measured dataset'},
        'frequency_hz': AXIS.frequency_hz.tolist(),
        'declared_plan': {'available_powers_per_face': list(POWERS), 'seeds': list(SEEDS),
            'gauss_orders': list(GAUSS_ORDERS), 'path_budget': PATH_BUDGET, 'stream_budget': STREAM_BUDGET,
            'per_path_reference_work': 'One factory call supplies the candidate; inside it two numerical age integrals and one closed-form calculation measure evidence. Four phase marks require no repeated transport solves.',
            'literal_phase_budget': PHASE_BUDGET, 'actual_transport_budget': TRANSPORT_BUDGET,
            'stopping_rule': 'Evaluate both seeds at each increasing grid; first passing actual-output gate at power >=8, else finish power10 uncertified',
            'error_definition': ERROR_DEFINITION,
            'comparison_scales': 'For a path: 4*tau^2*max(readout)^2, mean 2*tau*max(readout), response covariance scale*max(coupling). Stream uses tau=max(box length)/thermal sigma and multiplies by total boundary rate. These are unit-matched zero-signal scales, not claimed rigorous error bounds.',
            'density_count': 'Exactly once in physical boundary arrival rates; no occupancy or spectrum renormalization'},
        'grids': rows, 'ensemble_convergence': convergence,
        'final_candidate': {'certified': certified['certified'], 'path_evidence_passed': certified['path_evidence_passed'],
            'candidate_digest': certified['candidate_digest'], 'path_count': certified['path_count'],
            'consumed_path_count': certified['consumed_path_count'], 'model_scope': certified['model_scope'],
            'numerical_units': certified.get('numerical_units'), 'reasons': list(certified['reasons']),
            'spectra': json_value(certified['spectra']),
            'minimum_raw_greater_eigenvalue_s': (float(np.linalg.eigvalsh(certified['spectra']['greater']).min())
                                                if certified['spectra'] is not None else None)},
        'literal_four_phase_control_passed': phase_passed, 'streaming_callback_control_passed': callback_passed,
        'negative_controls': negatives, 'actual_constant_transport_control': constant,
        'actual_native_Rb_missing_evidence_control': native,
        'source_hashes_sha256': before, 'source_hashes_stable': before == after,
        'implementation_controls_passed': implementation_passed,
        'expected_controls_passed': bool(implementation_passed and certified['certified']),
        'environment': {'python': platform.python_version(), 'numpy': np.__version__, 'scipy': scipy.__version__,
            'blas_threads': 1, 'elapsed_s': time.monotonic()-started},
        'limitations': [
            'The zero-offset toy uses an imposed phase-mark family, not the coherence dynamics of a physical nonzero Rb beat.',
            'The 5% stream gate concerns this finite constant-atom integrand only. Independent scrambles are numerical controls, not rigorous statistical confidence bounds.',
            'Only atomic_inflow contributes because the toy has no jumps; nonzero microscopic jump-source ensemble convergence is not established here.',
            'The native reduced Rb calculation is one deliberately uncertified packet; passing its internal audit does not supply path or thermal-ensemble convergence.',
            'Returned retarded response is an integrated finite-age response to a global weak drive; it is not a spatial Maxwell M matrix.',
            'Uniform common phase gives the zero-cyclic time-averaged atomic stream. Other cyclic sectors and coherent mean lines are outside this artifact.',
            'Self-consistent nonlocal optical propagation, field commutators/SQL, finite seed saturation, full Zeeman physics and held-out experimental validation remain downstream.']}


def save_plot(report, path):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 3, figsize=(13, 4.1), layout='constrained')
    for kind, label in (('ensemble_refinement', 'Grid refinement'), ('independent_scramble', 'Independent scramble')):
        rows = report['ensemble_convergence']
        axes[0].loglog([6*2**row['power'] for row in rows],
            [max(item['errors'].values()) for row in rows for item in row['comparisons'] if item['kind'] == kind],
            'o-', label=label)
    axes[0].axhline(STREAM_BUDGET, color='grey', linestyle='--', label='Declared 5% budget')
    axes[0].set(xlabel='Boundary paths', ylabel='Largest output matrix error', title='Actual integrand convergence')
    axes[0].legend(fontsize=8)
    spectra = report['final_candidate']['spectra']
    if spectra is not None:
        for key, label in (('greater', 'Raw greater'), ('internal_greater', 'Connected single atom'), ('poisson_number', 'Poisson mean outer')):
            axes[1].plot(report['frequency_hz'], np.asarray(spectra[key]['real'])[:, 0, 0], 'o-', label=label)
    axes[1].set(xlabel='Laboratory RF (Hz)', ylabel='Toy atomic spectrum (s)', title='One common phase; mean outer retained')
    axes[1].legend(fontsize=8)
    negative = report['negative_controls']
    if 'delete_mean_outer_relative_raw_greater_error' in negative:
        axes[2].bar(['Delete mean\nouter term', 'Average mean\nbefore outer'],
            [negative['delete_mean_outer_relative_raw_greater_error'],
             negative['average_mean_before_outer_relative_number_error']], color=['#c75d6a', '#d99645'])
    axes[2].set(ylabel='Relative matrix error', title='Deliberately wrong Poisson controls')
    for ax in axes:
        ax.grid(axis='y', alpha=.2)
    fig.suptitle('Constant two-level phase-mark fixture | No thermal Rb or optical squeezing certification', fontsize=11)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('xb') as output:
        fig.savefig(output, format='png', dpi=160)
    plt.close(fig)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--plot', type=Path)
    args = parser.parse_args(argv)
    if args.plot is not None and args.output.resolve() == args.plot.resolve():
        raise ValueError('Report and plot require different paths')
    for path in (args.output, args.plot):
        if path is not None and path.exists():
            raise FileExistsError(path)
    report = build_report()
    if not report['source_hashes_stable']:
        raise RuntimeError('Source files changed during the audit; no immutable artifact written')
    payload = json.dumps(report, ensure_ascii=False, indent=2, allow_nan=False)+'\n'
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open('x', encoding='utf-8') as output:
        output.write(payload)
    if args.plot is not None:
        save_plot(report, args.plot)
    print(f'Wrote {args.output}; controls={report["expected_controls_passed"]}', flush=True)
    return 0 if report['expected_controls_passed'] else 1


if __name__ == '__main__':
    raise SystemExit(main())
