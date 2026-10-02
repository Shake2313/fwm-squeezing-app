"""Read-only grid aggregation and independent audit of a frozen thermal campaign.

The controller is captured separately, outside the frozen numerical inventory.
A passing path/grid audit is not thermal ensemble or optical certification.
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
from analysis.grand_challenge import thermal_campaign as c


SCHEMA = 'gabes-frozen-thermal-grid-audit-v1'
AGGREGATION_BUDGET = 1e-10


def controller_identity():
    raw = Path(__file__).read_bytes()
    token = hashlib.sha256(raw).hexdigest()
    return {'name': 'thermal_campaign_grid.py',
            'raw_sha256': token, 'raw_bytes': len(raw),
            'archive': 'controllers/thermal_campaign_grid-'+token+'.py'}


def archive_controller(directory):
    control = controller_identity()
    raw = Path(__file__).read_bytes()
    if hashlib.sha256(raw).hexdigest() != control['raw_sha256']:
        raise ValueError('Controller changed before archival')
    target = Path(directory)/control['archive']
    target.parent.mkdir(parents=True, exist_ok=True)
    try:
        with target.open('xb') as handle:
            handle.write(raw)
    except FileExistsError:
        if not target.is_file() or target.read_bytes() != raw:
            raise ValueError('Controller archive changed')
    return control


def selected_grid(plan, power, seed):
    """Validate every declared path and its native, unrenormalized arrival rate."""
    if (type(power) is not int or type(seed) is not int
            or power not in plan['powers'] or seed not in plan['seeds']):
        raise ValueError('Select one exact declared integer grid')
    rows = [row for row in plan['grids'] if (row['power'], row['seed']) == (power, seed)]
    if len(rows) != 1:
        raise ValueError('Select one exact declared grid')
    model = c.r.default_model()
    inflow = model.inflow(power, seed)
    row = rows[0]
    if (len(row['paths']) != len(inflow.rate_s_inverse)
            or row['mean_occupancy'] != inflow.mean_occupancy
            or row['equilibrium_nV'] != inflow.equilibrium_atom_number):
        raise ValueError('Declared grid count or occupancy differs')
    for index, entry in enumerate(row['paths']):
        if (type(entry['index']) is not int or entry['index'] != index
                or entry['physical_path'] != c.path_identity(inflow.path(index))
                or entry['rate_s_inverse'] != float(inflow.rate_s_inverse[index])):
            raise ValueError('Declared grid path or arrival rate differs')
    return model, inflow


def direct_weighted_audit(model, inflow, numerical, cache, row):
    """Check the stream against literal four-phase transforms of cached packets.

    Mathematically exact pump-only phase integration: charges +/-1 give only
    harmonics 0,+/-2 in D X D-dagger, so four roots of unity integrate them
    exactly. H and inflow state do not depend on the common phase. No additional
    atomic solves are warranted; this theorem does not cover finite seeds.
    """
    candidate = row['candidate']
    count = len(inflow.rate_s_inverse)
    if (not candidate['path_evidence_passed'] or candidate['spectra'] is None
            or candidate['consumed_path_count'] != count
            or len(row['path_ledgers']) != count
            or any(not item['passed'] for item in row['path_ledgers'])):
        return {'performed': False, 'passed': False,
                'reasons': ['Every declared path requires all five passing numerical calculations']}
    convention = model.convention()
    if tuple(convention.phase_charges) != (1, -1, -1, 1):
        raise ValueError('Independent four-phase audit requires the pump-only charge theorem')
    totals, used = None, []
    keys = ('greater', 'lesser', 'greater_by_source', 'lesser_by_source', 'retarded_response')
    roots = np.array([1, 1j, -1, -1j], dtype=complex)
    charges = np.array(convention.phase_charges)
    for index, rate in enumerate(inflow.rate_s_inverse):
        path = inflow.path(index)
        packet, record = cache.get(model, path, numerical.primary[-1])
        previous = row['path_ledgers'][index]['cache_records'][2]
        if (not record['hit'] or any(record.get(key) != previous.get(key)
                for key in ('key', 'record_sha256', 'payload_digest'))):
            raise ValueError('A finest packet changed between the path gate and direct aggregation')
        used.append(record)
        terms = {key: np.zeros_like(packet[key]) for key in keys}
        terms['poisson_number'] = np.zeros_like(packet['greater'])
        for root in roots:
            diagonal = root**charges
            for key in keys:
                terms[key] += np.einsum('i,...ij,j->...ij', diagonal,
                                       packet[key], diagonal.conj())/4
            mean = packet['mean_pulse']*diagonal
            terms['poisson_number'] += np.einsum('fi,fj->fij', mean, mean.conj())/4
        # Number noise is the phase average of m m-dagger, not the outer product
        # of the phase-averaged mean (which vanishes). Atomic-inflow covariance
        # is already in *_by_source and is not added a second time.
        for key in ('greater', 'lesser'):
            terms[key] += terms['poisson_number']
        if totals is None:
            totals = {key: np.zeros_like(value) for key, value in terms.items()}
        for key, value in terms.items():
            # Physical arrival rates contain density once. No nV, weight-sum,
            # residence-time or occupancy normalization is allowed here.
            totals[key] += float(rate)*value
    errors = c.r.stream_errors(candidate['spectra'], totals, row['comparison_scales'])
    passed = all(value <= AGGREGATION_BUDGET for value in errors.values())
    return {'performed': True, 'passed': passed, 'errors': errors,
        'budget': AGGREGATION_BUDGET, 'path_count': count, 'cache_records': used,
        'error_definition': c.r.ERROR_DEFINITION,
        'method': 'Literal four-root phase transforms and direct native-rate sum; each RF and named source checked',
        'limits': 'Conditional pump-only uniform common phase; no finite seed, ensemble convergence, optical channel or SQL',
        'reasons': [] if passed else ['Independent weighted aggregation exceeds its declared budget']}


def _new_output(output, directory):
    output = Path(output)
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    output = output.resolve()
    if output.is_relative_to((Path(directory)/'cache').resolve()):
        raise ValueError('Grid reports must stay outside the immutable packet cache')
    return output


def _loaded_from_capsule():
    for name, module in tuple(sys.modules.items()):
        if name == 'gabes' or name.startswith(('gabes.', 'analysis.grand_challenge.')):
            filename = getattr(module, '__file__', None)
            if filename and not Path(filename).resolve().is_relative_to(c.ROOT.resolve()):
                raise ValueError('Loaded numerical module originates outside the captured source capsule')


def coverage_inventory(directory, plan, numerical, cache, row):
    """Cheap filename inventory only; present files are not numerical evidence."""
    paths, jobs, entry_count = set(), set(), 0
    for declared in plan['grids']:
        model, inflow = selected_grid(plan, declared['power'], declared['seed'])
        entry_count += len(inflow.rate_s_inverse)
        for index in range(len(inflow.rate_s_inverse)):
            path = inflow.path(index)
            paths.add(c.r.digest(c.path_identity(path)))
            for spec in numerical.primary+numerical.reference:
                jobs.add(cache.key(model, path, spec))
    directory = Path(directory)/'cache'
    present = sum((directory/(key+'.json')).is_file() for key in jobs)
    verified = {record['key'] for ledger in row['path_ledgers'] for record in ledger['cache_records']}
    return {'declared_grid_count': len(plan['grids']), 'declared_path_entries': entry_count,
        'unique_physical_paths': len(paths), 'declared_unique_jobs': len(jobs),
        'present_job_files': present, 'missing_job_files': len(jobs)-present,
        'revalidated_current_grid_jobs': len(verified),
        'file_presence_is_validated_evidence': False,
        'scope': 'Read-only file-presence snapshot; only current-grid ledger records were reread and validated'}


def execute(directory, power, seed, output, *, expected_controller_sha256=None):
    directory = Path(directory).resolve()
    output = _new_output(output, directory)
    plan = c.load_plan(directory)
    c.require_capsule(plan)
    _loaded_from_capsule()
    controller = controller_identity()
    if (expected_controller_sha256 is not None
            and controller['raw_sha256'] != expected_controller_sha256):
        raise ValueError('Captured grid controller bytes differ from launch request')
    model, inflow = selected_grid(plan, power, seed)
    if archive_controller(directory) != controller:
        raise ValueError('Captured controller changed before aggregation')
    numerical = c.pilot.step_plan(plan['max_steps_s'])
    cache = c.CampaignCache(directory/'cache', plan)
    # Exact reuse is limited to unweighted packets for identical equations.
    # run_grid/build_path_evidence rebind current path source/digests, recompute
    # every numerical comparison, and apply THIS grid's physical arrival rates.
    candidate, row = c.r.run_grid(model, power, seed, numerical, cache)
    try:
        audit = direct_weighted_audit(model, inflow, numerical, cache, row)
    except (OSError, ValueError, TypeError, KeyError, RuntimeError, ArithmeticError) as error:
        audit = {'performed': True, 'passed': False,
                 'reasons': [type(error).__name__+': '+str(error)]}
    path_passed = candidate['path_evidence_passed']
    withheld_digest = None
    if not audit['passed'] and candidate['spectra'] is not None:
        withheld_digest = candidate['candidate_digest']
        candidate = dict(candidate, spectra=None, candidate_digest=None,
                         path_evidence_passed=False, certified=False,
                         reasons=tuple(audit['reasons']))
        row['candidate'] = candidate
        row['content_digest'] = c.r.row_digest(row)
    # A single grid cannot meet the declared three-grid/three-scramble gate.
    # Keep the real gate result and full row; never promote path success into
    # ensemble success or manufacture missing rows from nested packet reuse.
    gate = c.r.convergence_gate([row], powers=plan['powers'], seeds=plan['seeds'])
    coverage = coverage_inventory(directory, plan, numerical, cache, row)
    after = cache._check_sources()
    _loaded_from_capsule()
    if controller_identity() != controller:
        raise ValueError('Grid execution controller changed during the audit')
    return c.write_record(output, {'schema': SCHEMA, 'selection': [power, seed],
        'campaign_sha256': plan['record_sha256'], 'source_identity': plan['source_identity'],
        'source_bundle_sha256': plan['source_bundle_sha256'], 'controller': controller,
        'source_manifest_after': after, 'row': row, 'aggregation_audit': audit,
        'coverage_inventory': coverage,
        'path_evidence_passed': path_passed, 'grid_passed': bool(path_passed and audit['passed']),
        'withheld_candidate_digest': withheld_digest, 'convergence_gate': gate,
        'thermal_ensemble_converged': gate['passed'], 'physical_optical_prediction': False,
        'scope': 'Conditional pump-only atomic grid; independent aggregation arithmetic only, not experimental squeezing'})


def launch(directory, power, seed, output):
    directory = Path(directory).resolve()
    output = _new_output(output, directory)
    # Validate the requested grid without attempting live-source numerical work.
    # Full model/source/environment validation happens in the captured child.
    declared = c.records.read_sealed(directory/'plan.json')
    if (type(power) is not int or type(seed) is not int
            or sum((g['power'], g['seed']) == (power, seed) for g in declared['grids']) != 1):
        raise ValueError('Select one exact declared integer grid')
    archived = archive_controller(directory)
    scratch = ROOT/'.git'/'grand-challenge-runs'
    scratch.mkdir(parents=True, exist_ok=True)
    raw = Path(__file__).read_bytes()
    controller_hash = hashlib.sha256(raw).hexdigest()
    if controller_hash != archived['raw_sha256']:
        raise ValueError('Controller changed after archival')
    with tempfile.TemporaryDirectory(prefix='thermal-grid-', dir=scratch) as folder:
        c.extract_bundle(directory, folder)
        controller = Path(folder)/'thermal_campaign_grid.py'
        with controller.open('xb') as handle:
            handle.write(raw)
        # A root-level sidecar is outside both captured numerical directories.
        if c.source_paths(folder) != declared['source_paths']:
            raise ValueError('Grid controller altered the frozen numerical inventory')
        env = os.environ.copy()
        env.pop('PYTHONPATH', None)
        command = [sys.executable, str(controller), '--execute', str(directory),
                   '--grid', str(power), str(seed), '--output', str(output),
                   '--controller-sha256', controller_hash]
        subprocess.run(command, cwd=folder, env=env, check=True)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    action = parser.add_mutually_exclusive_group(required=True)
    action.add_argument('--run', type=Path)
    action.add_argument('--execute', type=Path, help=argparse.SUPPRESS)
    parser.add_argument('--grid', nargs=2, type=int, required=True, metavar=('POWER', 'SEED'))
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--controller-sha256', help=argparse.SUPPRESS)
    args = parser.parse_args(argv)
    if args.run:
        launch(args.run, *args.grid, args.output)
        return 0
    if args.controller_sha256 is None:
        parser.error('Captured --execute requires --controller-sha256')
    report = execute(args.execute, *args.grid, args.output,
                     expected_controller_sha256=args.controller_sha256)
    print(json.dumps({'grid_passed': report['grid_passed'],
                      'thermal_ensemble_converged': report['thermal_ensemble_converged'],
                      'report': str(args.output)}), flush=True)
    return 0 if report['grid_passed'] else 1


if __name__ == '__main__':
    raise SystemExit(main())
