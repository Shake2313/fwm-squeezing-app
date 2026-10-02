"""Compare two nested frozen thermal grids by rereading their actual path caches.

No atomic provider is invoked. Report seals alone are not numerical evidence:
both complete path gates, spectra and cache references are reconstructed first.
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
from tools import thermal_campaign_grid as g

c = g.c
SCHEMA = 'gabes-frozen-thermal-grid-comparison-v1'
REFERENCE_KEYS = ('key', 'record_sha256', 'payload_digest')


def captured_controllers():
    """Capture each external dependency once; these bytes are copied verbatim."""
    if Path(g.__file__).resolve().parent != Path(__file__).resolve().parent:
        raise ValueError('Grid helper must originate beside the comparison controller')
    result = []
    for filename in (Path(__file__), Path(g.__file__)):
        raw = filename.read_bytes()
        token = hashlib.sha256(raw).hexdigest()
        identity = {'name': filename.name, 'raw_sha256': token, 'raw_bytes': len(raw),
                    'archive': 'controllers/'+filename.stem+'-'+token+'.py'}
        result.append((identity, raw))
    return result


def archive_controllers(directory, captured):
    for identity, raw in captured:
        target = Path(directory)/identity['archive']
        target.parent.mkdir(parents=True, exist_ok=True)
        try:
            with target.open('xb') as handle:
                handle.write(raw)
        except FileExistsError:
            if not target.is_file() or target.read_bytes() != raw:
                raise ValueError('External controller archive changed')


def _controller_reference(directory, identity):
    token = identity['raw_sha256']
    if (identity['name'] != 'thermal_campaign_grid.py'
            or not isinstance(token, str) or len(token) != 64
            or any(char not in '0123456789abcdef' for char in token)
            or type(identity['raw_bytes']) is not int or identity['raw_bytes'] <= 0
            or identity['archive'] != 'controllers/thermal_campaign_grid-'+token+'.py'):
        raise ValueError('Grid report controller reference is invalid')
    raw = (Path(directory)/identity['archive']).read_bytes()
    if len(raw) != identity['raw_bytes'] or hashlib.sha256(raw).hexdigest() != token:
        raise ValueError('Grid report controller archive differs')


def read_grid_report(filename, directory, plan):
    filename = Path(filename).resolve()
    raw = filename.read_bytes()
    report = c.r.decode(c.records.read_sealed(filename))
    if filename.read_bytes() != raw:
        raise ValueError('Grid report changed while being read')
    if (report.get('schema') != g.SCHEMA
            or report.get('campaign_sha256') != plan['record_sha256']
            or report.get('source_identity') != plan['source_identity']
            or report.get('source_bundle_sha256') != plan['source_bundle_sha256']
            or c.provenance.source_identity(report['source_manifest_after']) != plan['source_identity']):
        raise ValueError('Grid report campaign or frozen-source reference differs')
    selection = report['selection']
    if not isinstance(selection, list) or len(selection) != 2:
        raise ValueError('Grid report must select exactly one declared grid')
    model, inflow = g.selected_grid(plan, *selection)
    row = report['row']
    if (row['power'] != selection[0] or row['seed'] != selection[1]
            or row['model_digest'] != c.r.digest(model.identity())
            or row['plan_digest'] != c.r.digest(c.pilot.step_plan(plan['max_steps_s']).identity())
            or row['content_digest'] != c.r.row_digest(row)):
        raise ValueError('Grid report row identity or content digest differs')
    if (report.get('grid_passed') is not True or report.get('path_evidence_passed') is not True
            or report.get('thermal_ensemble_converged') is not False
            or report.get('physical_optical_prediction') is not False
            or report.get('withheld_candidate_digest') is not None
            or row['candidate']['certified'] is not False):
        raise ValueError('Comparison requires passing, explicitly uncertified grid reports')
    _controller_reference(directory, report['controller'])
    reference = {'path': str(filename), 'file_sha256': hashlib.sha256(raw).hexdigest(),
                 'record_sha256': report['record_sha256'], 'selection': selection,
                 'controller': report['controller']}
    return report, reference, model, inflow


def _reference(record):
    # Storage locations and historical hit flags are not the physical identity.
    # They may change upon relocation; neither authorizes nor counts a new solve.
    return {name: record[name] for name in REFERENCE_KEYS}


def _row_semantics(row):
    value = {key: item for key, item in row.items()
             if key not in ('content_digest', 'cache_hits', 'cache_misses')}
    value['path_ledgers'] = [dict(ledger, cache_records=[_reference(record)
        for record in ledger['cache_records']]) for ledger in row['path_ledgers']]
    return value


def _audit_semantics(audit):
    return dict(audit, cache_records=[_reference(record) for record in audit['cache_records']])


def reread_grid(report, model, inflow, numerical, cache, plan):
    power, seed = report['selection']
    candidate, fresh = c.r.run_grid(model, power, seed, numerical, cache)
    if candidate['spectra'] is None or not candidate['path_evidence_passed']:
        raise ValueError('Fresh grid requires complete passing path evidence: '+str(candidate['reasons']))
    if c.r.digest(_row_semantics(report['row'])) != c.r.digest(_row_semantics(fresh)):
        raise ValueError('Saved grid row differs from freshly validated caches and current-path gates')
    direct = g.direct_weighted_audit(model, inflow, numerical, cache, fresh)
    if not direct['passed']:
        raise ValueError('Fresh independent grid aggregation did not pass')
    if c.r.digest(_audit_semantics(report['aggregation_audit'])) != c.r.digest(_audit_semantics(direct)):
        raise ValueError('Saved aggregation audit or cache references differ from actual reread')
    gate = c.r.convergence_gate([fresh], powers=plan['powers'], seeds=plan['seeds'])
    if c.r.digest(report['convergence_gate']) != c.r.digest(gate):
        raise ValueError('Saved grid convergence claim differs from the actual gate')
    return fresh, direct


def nested_indices(coarse_power, fine_power):
    if (type(coarse_power) is not int or type(fine_power) is not int
            or not 0 <= coarse_power < fine_power <= 10 or fine_power != coarse_power+1):
        raise ValueError('Nested comparison requires consecutive bounded powers')
    old_size, new_size = 2**coarse_power, 2**fine_power
    # Exact nested Sobol identity is per face, preserving the within-face index.
    # 2*i happens to work only for p0 -> p1; it is wrong for later refinements.
    return [(index, (index//old_size)*new_size+index % old_size)
            for index in range(6*old_size)]


def _bitwise_equal(left, right):
    a, b = np.asarray(left), np.asarray(right)
    return a.shape == b.shape and a.dtype == b.dtype and a.tobytes() == b.tobytes()


def _checked_read(model, path, spec, cache, ledger, level):
    packet, record = cache.get(model, path, spec)
    if record['hit'] is not True or _reference(record) != _reference(ledger['cache_records'][level]):
        raise ValueError('Cache reference changed after the fresh current-path gate')
    if ledger['current_path_source'] != path.source:
        raise ValueError('Current path source was not rebound')
    if level == 2 and c.r.stream.packet_digest(path, packet, model.convention()) != ledger['target_digest']:
        raise ValueError('Current finest packet digest was not rebound')
    return packet, record


def nested_reuse_audit(model, coarse_inflow, fine_inflow, coarse, fine, numerical, cache):
    pairs = nested_indices(coarse['power'], fine['power'])
    if (coarse['seed'] != fine['seed']
            or len(coarse_inflow.rate_s_inverse) != len(pairs)
            or len(fine_inflow.rate_s_inverse) != 2*len(pairs)):
        raise ValueError('Nested reuse requires the same scramble and complete native grids')
    rows, mapped = [], set()
    for old, new in pairs:
        before, after = coarse_inflow.path(old), fine_inflow.path(new)
        face, offset = divmod(old, 2**coarse['power'])
        if (new in mapped or int(coarse_inflow.face_index[old]) != face
                or int(fine_inflow.face_index[new]) != face
                or any(not _bitwise_equal(getattr(before, name), getattr(after, name))
                       for name in ('entry_position_m', 'velocity_m_s', 'residence_time_s'))
                or fine_inflow.rate_s_inverse[new] != .5*coarse_inflow.rate_s_inverse[old]):
            raise ValueError('Nested map does not preserve exact physical paths and half rates')
        mapped.add(new)
        lo, hi = coarse['path_ledgers'][old], fine['path_ledgers'][new]
        records = []
        for level, spec in enumerate(numerical.primary+numerical.reference):
            one, ref_one = _checked_read(model, before, spec, cache, lo, level)
            two, ref_two = _checked_read(model, after, spec, cache, hi, level)
            if _reference(ref_one) != _reference(ref_two):
                raise ValueError('Nested physical path references different immutable records')
            if any(not _bitwise_equal(one[name], two[name]) for name in c.r.METRICS):
                raise ValueError('Nested raw eight-metric packets are not bitwise identical')
            records.append({'level': level, **_reference(ref_one), 'all_eight_metrics_bitwise_equal': True})
        if (before.source == after.source) != (lo['target_digest'] == hi['target_digest']):
            raise ValueError('Nested target digests do not reflect current source rebinding')
        rows.append({'coarse_index': old, 'fine_index': new, 'face': face, 'offset': offset,
            'physical_path': c.path_identity(before), 'coarse_rate_s_inverse': float(coarse_inflow.rate_s_inverse[old]),
            'fine_rate_s_inverse': float(fine_inflow.rate_s_inverse[new]),
            'coarse_source': before.source, 'fine_source': after.source,
            'coarse_target_digest': lo['target_digest'], 'fine_target_digest': hi['target_digest'],
            'cache_records': records})
    return {'passed': True, 'paired_path_count': len(rows), 'shared_unique_jobs': 5*len(rows),
        'metric_names': c.r.METRICS, 'pairs': rows,
        'new_fine_indices': sorted(set(range(len(fine_inflow.rate_s_inverse)))-mapped),
        'new_solve_count': 0, 'historical_hit_flags_used_as_execution_evidence': False,
        'scope': 'Identical unweighted equations and all five stored calculations; rates and evidence rebound per grid'}


def _phase_terms(packet):
    names = ('greater', 'lesser', 'greater_by_source', 'lesser_by_source', 'retarded_response')
    terms = {name: np.zeros_like(packet[name]) for name in names}
    terms['poisson_number'] = np.zeros_like(packet['greater'])
    # The exact pump-only theorem has only harmonics 0,+/-2. Four roots integrate
    # them without an ODE solve; retain E[m m-dagger], even though E[m] vanishes.
    for root in (1+0j, 1j, -1+0j, -1j):
        diagonal = np.asarray(root)**np.array([1, -1, -1, 1])
        for name in names:
            terms[name] += np.einsum('i,...ij,j->...ij', diagonal, packet[name], diagonal.conj())/4
        mean = packet['mean_pulse']*diagonal
        terms['poisson_number'] += np.einsum('fi,fj->fij', mean, mean.conj())/4
    for name in ('greater', 'lesser'):
        terms[name] += terms['poisson_number']
    return terms


def nested_sum_audit(model, fine_inflow, coarse, fine, numerical, cache, reuse):
    additions = {name: np.zeros_like(fine['candidate']['spectra'][name]) for name in c.r.stream.STREAM_METRICS}
    references = []
    for index in reuse['new_fine_indices']:
        packet, record = _checked_read(model, fine_inflow.path(index), numerical.primary[-1],
                                      cache, fine['path_ledgers'][index], 2)
        terms = _phase_terms(packet)
        for name in additions:
            # Boundary rates already contain density once. No occupancy or
            # sum-of-weights normalization; no duplicate atomic-inflow term.
            additions[name] += float(fine_inflow.rate_s_inverse[index])*terms[name]
        references.append({'fine_index': index, **_reference(record)})
    reconstructed = {name: .5*coarse['candidate']['spectra'][name]+value for name, value in additions.items()}
    errors = c.r.stream_errors(fine['candidate']['spectra'], reconstructed, coarse['comparison_scales'])
    absolute = {name: float(np.max(np.abs(fine['candidate']['spectra'][name]-value)))
                for name, value in reconstructed.items()}
    return {'passed': all(error <= g.AGGREGATION_BUDGET for error in errors.values()),
        'errors': errors, 'budget': g.AGGREGATION_BUDGET, 'max_absolute_entry_residual': absolute,
        'error_definition': c.r.ERROR_DEFINITION, 'formula': 'S_fine = 0.5*S_coarse + sum(new paths at native fine rates)',
        'new_path_sum': additions, 'new_path_records': references,
        'scope': 'Rounding and normalization identity for six atomic stream metrics; no thermal or optical certification'}


def execute(directory, coarse_file, fine_file, output, *, expected_controller_hashes=None):
    directory = Path(directory).resolve()
    output = g._new_output(output, directory)
    plan = c.load_plan(directory)
    c.require_capsule(plan)
    g._loaded_from_capsule()
    captured = captured_controllers()
    identities = [item[0] for item in captured]
    if (expected_controller_hashes is not None
            and list(expected_controller_hashes) != [row['raw_sha256'] for row in identities]):
        raise ValueError('Captured comparison controllers differ from the launch request')
    archive_controllers(directory, captured)
    saved_coarse, coarse_reference, model, coarse_inflow = read_grid_report(coarse_file, directory, plan)
    saved_fine, fine_reference, _, fine_inflow = read_grid_report(fine_file, directory, plan)
    nested_indices(saved_coarse['selection'][0], saved_fine['selection'][0])
    if saved_coarse['selection'][1] != saved_fine['selection'][1]:
        raise ValueError('Nested comparison requires the same scramble seed')
    numerical = c.pilot.step_plan(plan['max_steps_s'])
    cache = c.CampaignCache(directory/'cache', plan)
    coarse, coarse_audit = reread_grid(saved_coarse, model, coarse_inflow, numerical, cache, plan)
    fine, fine_audit = reread_grid(saved_fine, model, fine_inflow, numerical, cache, plan)
    reuse = nested_reuse_audit(model, coarse_inflow, fine_inflow, coarse, fine, numerical, cache)
    split = nested_sum_audit(model, fine_inflow, coarse, fine, numerical, cache, reuse)
    errors = c.r.stream_errors(fine['candidate']['spectra'], coarse['candidate']['spectra'], coarse['comparison_scales'])
    diagnostic = {'errors': errors, 'budget': plan['ensemble_budget'],
        'passed': all(error <= plan['ensemble_budget'] for error in errors.values()),
        'error_definition': c.r.ERROR_DEFINITION, 'reference_scales': 'coarse grid comparison_scales',
        'scope': 'One actual refinement edge only; a passing 5% diagnostic is not ensemble certification'}
    gate = c.r.convergence_gate([coarse, fine], powers=plan['powers'], seeds=plan['seeds'])
    source_after = cache._check_sources()
    g._loaded_from_capsule()
    if [item[0] for item in captured_controllers()] != identities:
        raise ValueError('Comparison controller or grid helper changed during execution')
    for reference in (coarse_reference, fine_reference):
        if hashlib.sha256(Path(reference['path']).read_bytes()).hexdigest() != reference['file_sha256']:
            raise ValueError('Input grid report changed during comparison')
    return c.write_record(output, {'schema': SCHEMA, 'campaign_sha256': plan['record_sha256'],
        'source_identity': plan['source_identity'], 'source_bundle_sha256': plan['source_bundle_sha256'],
        'source_manifest_after': source_after, 'controllers': identities,
        'input_reports': [coarse_reference, fine_reference], 'rows': [coarse, fine],
        'fresh_grid_aggregation_audits': [coarse_audit, fine_audit],
        'nested_reuse': reuse, 'nested_sum': split, 'refinement_diagnostic': diagnostic,
        'convergence_gate': gate, 'audit_passed': bool(reuse['passed'] and split['passed']),
        'new_solve_count': 0, 'thermal_ensemble_converged': bool(gate['passed'] and split['passed']),
        'physical_optical_prediction': False,
        'scope': 'Conditional pump-only atomic comparison and exact nested reuse; no experimental squeezing prediction'})


def launch(directory, coarse_file, fine_file, output):
    directory = Path(directory).resolve()
    output = g._new_output(output, directory)
    coarse_file, fine_file = Path(coarse_file).resolve(), Path(fine_file).resolve()
    # Read seals before creating a capsule; scientific validation remains inside
    # the captured interpreter and is repeated from actual read-only caches.
    c.records.read_sealed(coarse_file)
    c.records.read_sealed(fine_file)
    declared = c.records.read_sealed(directory/'plan.json')
    captured = captured_controllers()
    archive_controllers(directory, captured)
    scratch = ROOT/'.git'/'grand-challenge-runs'
    scratch.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='thermal-compare-', dir=scratch) as folder:
        c.extract_bundle(directory, folder)
        controllers = Path(folder)/'tools'
        controllers.mkdir()
        for identity, raw in captured:
            with (controllers/identity['name']).open('xb') as handle:
                handle.write(raw)
        if c.source_paths(folder) != declared['source_paths']:
            raise ValueError('External comparison controllers altered the frozen numerical inventory')
        env = os.environ.copy()
        env.pop('PYTHONPATH', None)
        command = [sys.executable, str(controllers/'thermal_campaign_compare.py'),
            '--execute', str(directory), '--coarse', str(coarse_file), '--fine', str(fine_file),
            '--output', str(output), '--controllers-sha256', *[item[0]['raw_sha256'] for item in captured]]
        subprocess.run(command, cwd=folder, env=env, check=True)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    action = parser.add_mutually_exclusive_group(required=True)
    action.add_argument('--run', type=Path)
    action.add_argument('--execute', type=Path, help=argparse.SUPPRESS)
    parser.add_argument('--coarse', type=Path, required=True)
    parser.add_argument('--fine', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--controllers-sha256', nargs=2, help=argparse.SUPPRESS)
    args = parser.parse_args(argv)
    if args.run:
        launch(args.run, args.coarse, args.fine, args.output)
        return 0
    if args.controllers_sha256 is None:
        parser.error('Captured --execute requires both --controllers-sha256 values')
    report = execute(args.execute, args.coarse, args.fine, args.output,
                     expected_controller_hashes=args.controllers_sha256)
    print(json.dumps({'audit_passed': report['audit_passed'],
        'refinement_diagnostic_passed': report['refinement_diagnostic']['passed'],
        'thermal_ensemble_converged': report['thermal_ensemble_converged'],
        'new_solve_count': 0, 'report': str(args.output)}), flush=True)
    return 0 if report['audit_passed'] else 1


if __name__ == '__main__':
    raise SystemExit(main())
