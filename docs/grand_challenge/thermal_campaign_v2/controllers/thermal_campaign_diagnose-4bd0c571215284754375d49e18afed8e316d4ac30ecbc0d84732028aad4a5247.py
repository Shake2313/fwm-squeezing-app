"""Attribute a directed frozen-grid difference to complete physical strata.

This read-only cache diagnostic never solves an atom or certifies convergence.
It preserves signed matrix sums and every RF/source comparison. Strata describe
where quadrature differences occur, not independent errors or missing physics.
"""

import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile

import numpy as np

ROOT = Path(__file__).resolve().parent
if ROOT.name == 'tools':
    ROOT = ROOT.parent
sys.path.insert(0, str(ROOT))
from tools import thermal_campaign_compare as comparison

g, c = comparison.g, comparison.c
SCHEMA = 'gabes-frozen-thermal-strata-diagnostic-v1'
FACE_LABELS = ('x_lower', 'x_upper', 'y_lower', 'y_upper', 'z_lower', 'z_upper')
STRATA = (
    ('inlet_face', (), FACE_LABELS, 'integer face index; normal axis = face // 2'),
    ('residence_time_us', (.25, .5, 1., 2., 4.), (), 'residence time in microseconds'),
    ('longitudinal_velocity_sigma', (-2., -1., 0., 1., 2.), (), 'signed v_z / sqrt(k_B T / mass)'),
    ('transverse_speed_sigma', (.5, 1., 2., 3.), (), 'sqrt(v_x^2 + v_y^2) / sqrt(k_B T / mass)'),
)


def captured_controllers():
    if Path(comparison.__file__).resolve().parent != Path(__file__).resolve().parent:
        raise ValueError('Comparison helper must originate beside the diagnostic controller')
    raw = Path(__file__).read_bytes()
    token = hashlib.sha256(raw).hexdigest()
    identity = {'name': Path(__file__).name, 'raw_sha256': token, 'raw_bytes': len(raw),
                'archive': 'controllers/thermal_campaign_diagnose-'+token+'.py'}
    return [(identity, raw), *comparison.captured_controllers()]


def comparison_kind(reference, candidate):
    before, after = tuple(reference['selection']), tuple(candidate['selection'])
    if before == after:
        raise ValueError('Diagnostic requires two distinct grid selections')
    if before[1] == after[1] and after[0] == before[0]+1:
        return 'ensemble_refinement'
    if before[0] == after[0] and before[1] != after[1]:
        return 'independent_scramble'
    raise ValueError('Use adjacent coarse-to-fine powers with one seed, or equal powers with distinct seeds')


def partitions(model, inflow):
    """Fixed, exhaustive disjoint partitions; thresholds belong to upper bins."""
    sigma = float(np.sqrt(c.r.c.KB*model.temperature_K/model.mass_kg))
    velocity = np.asarray(inflow.velocity_m_s)
    values = (np.asarray(inflow.face_index), np.asarray(inflow.residence_time_s)*1e6,
              velocity[:, 2]/sigma, np.linalg.norm(velocity[:, :2], axis=1)/sigma)
    result = {}
    for (name, thresholds, labels, quantity), value in zip(STRATA, values):
        if value.shape != inflow.rate_s_inverse.shape or not np.isfinite(value).all():
            raise ValueError('Finite per-path stratum coordinates required')
        if name == 'inlet_face':
            assignment = value.astype(int)
            if not np.array_equal(assignment, value) or np.any((assignment < 0) | (assignment > 5)):
                raise ValueError('Invalid inlet face')
            groups = [{'label': label, 'face_index': index} for index, label in enumerate(labels)]
        else:
            if name != 'longitudinal_velocity_sigma' and np.any(value < 0):
                raise ValueError('Residence time and transverse speed cannot be negative')
            assignment = np.searchsorted(thresholds, value, side='right')
            bounds = [None, *thresholds, None]
            groups = [{'label': str(index), 'lower_inclusive': bounds[index],
                       'upper_exclusive': bounds[index+1]} for index in range(len(bounds)-1)]
        for index, group in enumerate(groups):
            group['indices'] = np.flatnonzero(assignment == index).tolist()
        # Empty bins are kept. No quantile selection, discarded tails, or
        # renormalization to the occupied subset can enter this partition.
        if sorted(index for group in groups for index in group['indices']) != list(range(len(value))):
            raise ValueError('Strata must partition every physical path exactly once')
        result[name] = {'quantity': quantity, 'thresholds': list(thresholds), 'groups': groups,
                        'sigma_m_s': sigma, 'assignment': assignment.tolist()}
    return result


def _zeros(spectra):
    return {name: np.zeros_like(spectra[name]) for name in c.r.stream.STREAM_METRICS}


def grouped_stream(model, inflow, row, numerical, cache):
    """Literal native-rate contributions, reusing only finest unweighted packets."""
    if tuple(model.convention().phase_charges) != (1, -1, -1, 1):
        raise ValueError('Four-phase pump-only theorem required')
    result = partitions(model, inflow)
    spectra = row['candidate']['spectra']
    for partition in result.values():
        for group in partition['groups']:
            group['spectra'] = _zeros(spectra)
            group['native_arrival_rate_s_inverse'] = float(np.sum(inflow.rate_s_inverse[group['indices']]))
    references = []
    for index, rate in enumerate(inflow.rate_s_inverse):
        packet, record = comparison._checked_read(model, inflow.path(index), numerical.primary[-1],
                                                  cache, row['path_ledgers'][index], 2)
        terms = comparison._phase_terms(packet)
        references.append({'path_index': index, **comparison._reference(record)})
        for partition in result.values():
            group = partition['groups'][partition['assignment'][index]]
            for name in c.r.stream.STREAM_METRICS:
                # Density appears once in the physical boundary arrival rate.
                # Never divide by occupancy, nV, counts, or sums of weights.
                group['spectra'][name] += float(rate)*terms[name]
    for partition in result.values():
        total = _zeros(spectra)
        for group in partition['groups']:
            for name in total:
                total[name] += group['spectra'][name]
        errors = c.r.stream_errors(total, spectra, row['comparison_scales'])
        if any(error > g.AGGREGATION_BUDGET for error in errors.values()):
            raise ValueError('Stratified sums do not reconstruct the validated grid')
        partition['reconstruction_errors'] = errors
    return result, references


def metric_attribution(differences, total_difference, reference, scale):
    """Additive signed projection, with magnitudes exposing cancellation.

    For matrices d_g and D=sum_g d_g, p_g=Re Tr(D† d_g)/||D||/den.
    Sum p_g = ||D||/den (the directed error) at each RF/source. Individual
    magnitudes ||d_g||/den need not add to that error and may exceed it.
    At D=0 the projection is defined as zero while magnitudes remain visible.
    """
    delta = np.asarray(total_difference)
    denominator = np.maximum(np.linalg.norm(reference, axis=(-2, -1)), 128*np.finfo(float).eps*scale)
    norm = np.linalg.norm(delta, axis=(-2, -1))
    unit = np.zeros_like(delta)
    np.divide(delta, norm[..., None, None], out=unit, where=norm[..., None, None] > 0)
    error = norm/denominator
    groups = []
    for difference in differences:
        magnitude = np.linalg.norm(difference, axis=(-2, -1))/denominator
        projection = np.real(np.sum(unit.conj()*difference, axis=(-2, -1)))/denominator
        groups.append({'relative_magnitude': magnitude, 'signed_projection': projection})
    projection_sum = sum((group['signed_projection'] for group in groups), np.zeros_like(error))
    residual = float(np.max(np.abs(projection_sum-error)))
    # Absolute dimensionless residual avoids inflating cancellation near D=0.
    tolerance = g.AGGREGATION_BUDGET*max(1., float(np.max(error)))
    if not np.isfinite(residual) or residual > tolerance:
        raise ValueError('Signed projections do not reconstruct directed errors')
    magnitude_sum = sum((group['relative_magnitude'] for group in groups), np.zeros_like(error))
    if not all(np.isfinite(value).all() for value in (error, magnitude_sum, projection_sum)):
        raise ValueError('Nonfinite stratum attribution')
    return {'relative_error': error, 'denominator': denominator,
            'projection_sum_residual': residual, 'projection_sum_tolerance': tolerance,
            'sum_relative_magnitudes': magnitude_sum, 'groups': groups}


def _location(metric, index, model):
    rf = int(index[-1])
    return {'array_index': list(index), 'rf_index': rf,
            'frequency_hz': float(model.analysis_axis.omega_rad_s[rf]/(2*np.pi)),
            'source': model.convention().source_names[index[0]] if metric.endswith('_by_source') else None}


def difference_diagnostic(reference, candidate, before, after, model, budget):
    left, right = reference['candidate']['spectra'], candidate['candidate']['spectra']
    # This is precisely the directed formal ensemble gate: candidate-reference,
    # reference norm in the denominator, CANDIDATE SI floor. Reverse comparisons
    # have their own denominator; neither directions nor sources are averaged.
    scales = candidate['comparison_scales']
    errors = c.r.stream_errors(right, left, scales)
    delta = {name: right[name]-left[name] for name in c.r.stream.STREAM_METRICS}
    result = {}
    for partition_name in before:
        old, new = before[partition_name], after[partition_name]
        if (old['thresholds'] != new['thresholds'] or old['sigma_m_s'] != new['sigma_m_s']
                or len(old['groups']) != len(new['groups'])):
            raise ValueError('Both grids require identical physical stratum definitions')
        groups = []
        for lo, hi in zip(old['groups'], new['groups']):
            groups.append({'definition': {key: value for key, value in lo.items()
                           if key not in ('indices', 'spectra', 'native_arrival_rate_s_inverse')},
                'reference_indices': lo['indices'], 'candidate_indices': hi['indices'],
                'reference_rate_s_inverse': lo['native_arrival_rate_s_inverse'],
                'candidate_rate_s_inverse': hi['native_arrival_rate_s_inverse'],
                'difference_spectra': {name: hi['spectra'][name]-lo['spectra'][name] for name in delta}})
        summed = {name: sum((group['difference_spectra'][name] for group in groups),
                           np.zeros_like(value)) for name, value in delta.items()}
        # Compare the reconstructed candidate to the actual candidate, using
        # the same full-grid reference denominator even if the difference is 0.
        reconstruction = c.r.stream_errors({name: left[name]+summed[name] for name in delta}, right, scales)
        if any(value > g.AGGREGATION_BUDGET for value in reconstruction.values()):
            raise ValueError('Stratum differences fail full-grid reconstruction')
        metrics = {}
        for name in delta:
            info = metric_attribution([group['difference_spectra'][name] for group in groups],
                                      delta[name], left[name], scales[name])
            error_array = info['relative_error']
            if float(np.max(error_array)) != errors[name]:
                raise ValueError('Directed diagnostic differs from the formal stream error')
            worst = tuple(int(i) for i in np.unravel_index(np.argmax(error_array), error_array.shape))
            rankings = []
            for index, item in enumerate(info.pop('groups')):
                groups[index].setdefault('metrics', {})[name] = item
                rankings.append({'group_index': index,
                    'relative_magnitude': float(item['relative_magnitude'][worst]),
                    'signed_projection': float(item['signed_projection'][worst])})
            info.update(worst_location=_location(name, worst, model),
                        worst_relative_error=errors[name],
                        groups_at_worst=sorted(rankings, key=lambda item: -item['relative_magnitude']))
            metrics[name] = info
        result[partition_name] = {'quantity': old['quantity'], 'thresholds': old['thresholds'],
            'sigma_m_s': old['sigma_m_s'], 'groups': groups, 'metrics': metrics,
            'reference_reconstruction_errors': old['reconstruction_errors'],
            'candidate_reconstruction_errors': new['reconstruction_errors'],
            'difference_reconstruction_errors': reconstruction}
    return {'candidate_grid': [candidate['power'], candidate['seed']],
        'reference_grid': [reference['power'], reference['seed']], 'errors': errors,
        'budget': budget, 'passed': all(value <= budget for value in errors.values()),
        'error_definition': c.r.ERROR_DEFINITION, 'comparison_scales': 'candidate grid, as in convergence_gate',
        'total_difference_spectra': delta, 'partitions': result,
        'attribution_definition': 'Re(sum(conj(D/||D||)*d_group))/denominator; zero when D=0',
        'interpretation': 'Signed projections add to each directed RF/source error. Group magnitudes do not add; cancellation may make their sum larger. These partitions overlap each other and must never be added across partition kinds.',
        'diagnostic_only': True}


def execute(directory, reference_file, candidate_file, output, *, expected_controller_hashes=None):
    directory = Path(directory).resolve()
    output = g._new_output(output, directory)
    plan = c.load_plan(directory)
    c.require_capsule(plan)
    g._loaded_from_capsule()
    originals = {path: path.read_bytes() for path in (directory/'plan.json', directory/'sources.zip')}
    captured = captured_controllers()
    identities = [identity for identity, _ in captured]
    if (expected_controller_hashes is not None
            and list(expected_controller_hashes) != [identity['raw_sha256'] for identity in identities]):
        raise ValueError('Captured diagnostic controllers differ from launch request')
    comparison.archive_controllers(directory, captured)
    saved_reference, ref, model, reference_inflow = comparison.read_grid_report(reference_file, directory, plan)
    saved_candidate, cand, _, candidate_inflow = comparison.read_grid_report(candidate_file, directory, plan)
    kind = comparison_kind(saved_reference, saved_candidate)
    numerical = c.pilot.step_plan(plan['max_steps_s'])
    cache = c.CampaignCache(directory/'cache', plan)
    reference, reference_audit = comparison.reread_grid(saved_reference, model, reference_inflow, numerical, cache, plan)
    candidate, candidate_audit = comparison.reread_grid(saved_candidate, model, candidate_inflow, numerical, cache, plan)
    if (kind == 'independent_scramble'
            and reference['candidate']['candidate_digest'] == candidate['candidate']['candidate_digest']):
        raise ValueError('Independent scrambles cannot reuse one candidate computation')
    before, before_records = grouped_stream(model, reference_inflow, reference, numerical, cache)
    after, after_records = grouped_stream(model, candidate_inflow, candidate, numerical, cache)
    diagnostic = difference_diagnostic(reference, candidate, before, after, model, plan['ensemble_budget'])
    diagnostic['kind'] = kind
    # Reread all five consumed records again. A result changed after its initial
    # path gate must not acquire a new diagnostic seal after attribution.
    for row, inflow in ((reference, reference_inflow), (candidate, candidate_inflow)):
        for index, ledger in enumerate(row['path_ledgers']):
            for level, spec in enumerate(numerical.primary+numerical.reference):
                comparison._checked_read(model, inflow.path(index), spec, cache, ledger, level)
    source_after = cache._check_sources()
    g._loaded_from_capsule()
    if [identity for identity, _ in captured_controllers()] != identities:
        raise ValueError('Diagnostic controller or captured dependency changed during execution')
    for path, raw in originals.items():
        if path.read_bytes() != raw:
            raise ValueError('Campaign plan or source bundle changed during diagnostic')
    for item in (ref, cand):
        if hashlib.sha256(Path(item['path']).read_bytes()).hexdigest() != item['file_sha256']:
            raise ValueError('Input grid report changed during diagnostic')
        comparison._controller_reference(directory, item['controller'])
    return c.write_record(output, {'schema': SCHEMA, 'campaign_sha256': plan['record_sha256'],
        'source_identity': plan['source_identity'], 'source_bundle_sha256': plan['source_bundle_sha256'],
        'source_manifest_after': source_after, 'controllers': identities,
        'input_reports': [ref, cand], 'fresh_grid_aggregation_audits': [reference_audit, candidate_audit],
        'finest_cache_records': [before_records, after_records], 'directed_diagnostic': diagnostic,
        'audit_passed': True, 'new_solve_count': 0, 'thermal_ensemble_converged': False,
        'physical_optical_prediction': False, 'historical_hit_flags_used_as_execution_evidence': False,
        'scope': 'Complete conditional pump-only atomic grid difference, attributed without renormalization. No error bound, cause identification, ensemble certification, optical channel or experimental prediction.'})


def launch(directory, reference_file, candidate_file, output):
    directory = Path(directory).resolve()
    output = g._new_output(output, directory)
    reference_file, candidate_file = Path(reference_file).resolve(), Path(candidate_file).resolve()
    for filename in (reference_file, candidate_file):
        c.records.read_sealed(filename)
    declared = c.records.read_sealed(directory/'plan.json')
    captured = captured_controllers()
    comparison.archive_controllers(directory, captured)
    scratch = ROOT/'.git'/'grand-challenge-runs'
    scratch.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='thermal-diagnose-', dir=scratch) as folder:
        c.extract_bundle(directory, folder)
        controllers = Path(folder)/'tools'
        controllers.mkdir()
        for identity, raw in captured:
            with (controllers/identity['name']).open('xb') as handle:
                handle.write(raw)
        if c.source_paths(folder) != declared['source_paths']:
            raise ValueError('Diagnostic controllers altered the frozen numerical inventory')
        env = os.environ.copy()
        env.pop('PYTHONPATH', None)
        command = [sys.executable, str(controllers/'thermal_campaign_diagnose.py'),
            '--execute', str(directory), '--reference', str(reference_file), '--candidate', str(candidate_file),
            '--output', str(output), '--controllers-sha256', *[identity['raw_sha256'] for identity, _ in captured]]
        subprocess.run(command, cwd=folder, env=env, check=True)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    action = parser.add_mutually_exclusive_group(required=True)
    action.add_argument('--run', type=Path)
    action.add_argument('--execute', type=Path, help=argparse.SUPPRESS)
    parser.add_argument('--reference', type=Path, required=True)
    parser.add_argument('--candidate', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--controllers-sha256', nargs=3, help=argparse.SUPPRESS)
    args = parser.parse_args(argv)
    if args.run:
        launch(args.run, args.reference, args.candidate, args.output)
        return 0
    if args.controllers_sha256 is None:
        parser.error('Captured --execute requires all three --controllers-sha256 values')
    report = execute(args.execute, args.reference, args.candidate, args.output,
                     expected_controller_hashes=args.controllers_sha256)
    print(json.dumps({'audit_passed': report['audit_passed'],
        'directed_diagnostic_passed': report['directed_diagnostic']['passed'],
        'thermal_ensemble_converged': False, 'new_solve_count': 0, 'report': str(args.output)}), flush=True)
    return 0 if report['audit_passed'] else 1


if __name__ == '__main__':
    raise SystemExit(main())
