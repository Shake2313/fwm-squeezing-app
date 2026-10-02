"""Read-only evidence audit for any nonempty subset of one frozen campaign.

Every supplied grid is reconstructed from its five-calculation path caches.
Subset integrity and complete thermal convergence are separate conclusions.
"""

import argparse
import hashlib
import itertools
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parent
if ROOT.name == 'tools':
    ROOT = ROOT.parent
sys.path.insert(0, str(ROOT))
from tools import thermal_campaign_compare as comparison

g, c = comparison.g, comparison.c
SCHEMA = 'gabes-frozen-thermal-ensemble-audit-v1'


def captured_controllers():
    if Path(comparison.__file__).resolve().parent != Path(__file__).resolve().parent:
        raise ValueError('Comparison helper must originate beside the ensemble controller')
    filename = Path(__file__)
    raw = filename.read_bytes()
    token = hashlib.sha256(raw).hexdigest()
    identity = {'name': filename.name, 'raw_sha256': token, 'raw_bytes': len(raw),
                'archive': 'controllers/'+filename.stem+'-'+token+'.py'}
    return [(identity, raw), *comparison.captured_controllers()]


def input_paths(grid_files):
    if isinstance(grid_files, (str, Path)):
        raise ValueError('Supply a nonempty sequence of grid report paths')
    paths = [Path(filename).resolve() for filename in grid_files]
    if not paths:
        raise ValueError('At least one explicit grid report is required')
    if len(set(paths)) != len(paths):
        raise ValueError('Duplicate grid report paths')
    for path in paths:
        if not path.is_file():
            raise FileNotFoundError(path)
    return paths


def required_grids(plan):
    powers, seeds = tuple(plan['powers']), tuple(plan['seeds'])
    if (len(powers) < 3 or any(type(p) is not int or not 0 <= p <= 10 for p in powers)
            or tuple(sorted(set(powers))) != powers
            or len(seeds) < 3 or len(set(seeds)) != len(seeds)
            or any(type(s) is not int or s < 0 for s in seeds)):
        raise ValueError('Preserve the complete declared powers and independent seeds')
    required = list(itertools.product(powers, seeds))
    declared = [(row['power'], row['seed']) for row in plan['grids']]
    if len(declared) != len(required) or set(declared) != set(required):
        raise ValueError('Campaign grid declarations do not cover the required powers and seeds')
    return required


def coverage(plan, rows):
    required = required_grids(plan)
    supplied = [(row['power'], row['seed']) for row in rows]
    if len(set(supplied)) != len(supplied) or not set(supplied) <= set(required):
        raise ValueError('Duplicate or undeclared grid in ensemble audit')
    present = set(supplied)
    missing = [key for key in required if key not in present]
    return {'required_powers': plan['powers'], 'required_seeds': plan['seeds'],
        'required_grids': required, 'provided_grids': [key for key in required if key in present],
        'missing_grids': missing, 'required_grid_count': len(required),
        'provided_grid_count': len(supplied), 'complete': not missing}


def pairwise_diagnostics(rows, plan):
    """Measured available edges; never fabricate missing controls for the gate."""
    indexed = {(row['power'], row['seed']): row for row in rows}
    refinements, scrambles = [], []

    def measure(kind, current, previous):
        if (kind == 'independent_scramble'
                and current['candidate']['candidate_digest'] == previous['candidate']['candidate_digest']):
            raise ValueError('Independent scrambles cannot reuse one candidate computation')
        # Match the existing formal gate's orientation and dimensional floors:
        # current - reference, with CURRENT grid SI comparison scales.
        errors = c.r.stream_errors(current['candidate']['spectra'],
            previous['candidate']['spectra'], current['comparison_scales'])
        return {'kind': kind, 'candidate_grid': [current['power'], current['seed']],
            'reference_grid': [previous['power'], previous['seed']],
            'candidate_digest': current['candidate']['candidate_digest'],
            'reference_digest': previous['candidate']['candidate_digest'],
            'errors': errors, 'budget': plan['ensemble_budget'],
            'passed': all(value <= plan['ensemble_budget'] for value in errors.values())}

    for seed in plan['seeds']:
        for previous, power in zip(plan['powers'][:-1], plan['powers'][1:]):
            if (previous, seed) not in indexed or (power, seed) not in indexed:
                continue
            item = measure('ensemble_refinement', indexed[power, seed], indexed[previous, seed])
            item['adjacent_planned_powers'] = True
            refinements.append(item)
    for power in plan['powers']:
        available = [seed for seed in plan['seeds'] if (power, seed) in indexed]
        for current, previous in itertools.permutations(available, 2):
            scrambles.append(measure('independent_scramble', indexed[power, current], indexed[power, previous]))
    return {'refinements': refinements, 'independent_scrambles': scrambles,
        'error_definition': c.r.ERROR_DEFINITION, 'comparison_scales': 'candidate grid, as in convergence_gate',
        'diagnostic_only': True,
        'scope': 'Available consecutive declared-power edges and directed seed pairs per power; missing controls stay missing'}


def loaded_source_provenance(manifest):
    """Record actual imported numeric paths and bytes, including helper origins."""
    g._loaded_from_capsule()
    root = c.ROOT.resolve()
    expected = {row['path']: row for row in manifest['files']}
    loaded = []
    for name, module in sorted(tuple(sys.modules.items())):
        if name != 'gabes' and not name.startswith(('gabes.', 'analysis.grand_challenge.')):
            continue
        filename = getattr(module, '__file__', None)
        if filename is None:
            continue
        path = Path(filename).resolve()
        relative = path.relative_to(root).as_posix()
        if relative not in expected:
            raise ValueError('Loaded numerical source is outside the frozen inventory: '+relative)
        raw_hash = hashlib.sha256(path.read_bytes()).hexdigest()
        if raw_hash != expected[relative]['raw_sha256']:
            raise ValueError('Loaded numerical source bytes changed: '+relative)
        loaded.append({'module': name, 'source_path': relative, 'raw_sha256': raw_hash,
                       'normalized_sha256': expected[relative]['normalized_sha256']})
    for filename in (Path(__file__), Path(comparison.__file__), Path(g.__file__)):
        if not filename.resolve().is_relative_to(root):
            raise ValueError('External controller loaded outside the source capsule')
    return {'numerical_modules': loaded,
        'boundary': 'Fresh source-capsule interpreter; controllers separately archived outside the frozen numerical inventory'}


def execute(directory, grid_files, output, *, expected_controller_hashes=None):
    directory = Path(directory).resolve()
    output = g._new_output(output, directory)
    paths = input_paths(grid_files)
    plan = c.load_plan(directory)
    c.require_capsule(plan)
    required = required_grids(plan)
    captured = captured_controllers()
    identities = [item[0] for item in captured]
    if (expected_controller_hashes is not None
            and list(expected_controller_hashes) != [identity['raw_sha256'] for identity in identities]):
        raise ValueError('Captured ensemble controllers differ from the launch request')
    comparison.archive_controllers(directory, captured)
    loaded_source_provenance(plan['source_manifest'])
    supplied = {}
    for path in paths:
        report, reference, model, inflow = comparison.read_grid_report(path, directory, plan)
        key = tuple(report['selection'])
        if key in supplied:
            raise ValueError('Duplicate grid selection in supplied reports')
        supplied[key] = (report, reference, model, inflow)
    numerical = c.pilot.step_plan(plan['max_steps_s'])
    cache = c.CampaignCache(directory/'cache', plan)
    rows, references, audits = [], [], []
    for key in required:
        if key not in supplied:
            continue
        report, reference, model, inflow = supplied[key]
        # Exact packet reuse never carries a convergence label across grids.
        # The existing helper rereads all five calculations, rebinds each path,
        # recomputes every path gate, then independently checks weighted spectra.
        # No provider is supplied; a cache miss is a failed audit, not a solve.
        row, audit = comparison.reread_grid(report, model, inflow, numerical, cache, plan)
        rows.append(row)
        references.append(reference)
        audits.append({'selection': key, 'audit': audit})
    covered = coverage(plan, rows)
    diagnostics = pairwise_diagnostics(rows, plan)
    gate = c.r.convergence_gate(rows, powers=plan['powers'], seeds=plan['seeds'])
    after = cache._check_sources()
    origins = loaded_source_provenance(after)
    if [item[0] for item in captured_controllers()] != identities:
        raise ValueError('Ensemble controller or captured dependency changed during execution')
    for reference in references:
        if hashlib.sha256(Path(reference['path']).read_bytes()).hexdigest() != reference['file_sha256']:
            raise ValueError('Input grid report changed during ensemble audit')
        comparison._controller_reference(directory, reference['controller'])
    # The formal gate can assess a completed prefix of a longer declaration.
    # This campaign-level claim additionally requires ALL declared grids/seeds.
    converged = bool(covered['complete'] and gate['passed'])
    return c.write_record(output, {'schema': SCHEMA, 'campaign_sha256': plan['record_sha256'],
        'source_identity': plan['source_identity'], 'source_bundle_sha256': plan['source_bundle_sha256'],
        'source_manifest_after': after, 'loaded_source_provenance': origins, 'controllers': identities,
        'required_powers': plan['powers'], 'required_seeds': plan['seeds'],
        'input_reports': references, 'coverage': covered, 'rows': rows,
        'fresh_grid_aggregation_audits': audits, 'pairwise_diagnostics': diagnostics,
        'convergence_gate': gate, 'audit_passed': True,
        'audit_scope': 'Every supplied grid is valid against current immutable caches; omitted grids are not audited',
        'thermal_ensemble_converged': converged, 'physical_optical_prediction': False,
        'new_solve_count': 0, 'historical_hit_flags_used_as_execution_evidence': False,
        'scope': 'Conditional open-column pump-only atomic stream; no optical channel, SQL or experimental squeezing certification'})


def launch(directory, grid_files, output):
    directory = Path(directory).resolve()
    output = g._new_output(output, directory)
    paths = input_paths(grid_files)
    for path in paths:
        c.records.read_sealed(path)
    declared = c.records.read_sealed(directory/'plan.json')
    captured = captured_controllers()
    comparison.archive_controllers(directory, captured)
    scratch = ROOT/'.git'/'grand-challenge-runs'
    scratch.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='thermal-ensemble-', dir=scratch) as folder:
        c.extract_bundle(directory, folder)
        controllers = Path(folder)/'tools'
        controllers.mkdir()
        for identity, raw in captured:
            with (controllers/identity['name']).open('xb') as handle:
                handle.write(raw)
        if c.source_paths(folder) != declared['source_paths']:
            raise ValueError('Ensemble controllers altered the frozen numerical inventory')
        env = os.environ.copy()
        env.pop('PYTHONPATH', None)
        command = [sys.executable, str(controllers/'thermal_campaign_ensemble.py'),
            '--execute', str(directory), '--grids', *map(str, paths), '--output', str(output),
            '--controllers-sha256', *[identity['raw_sha256'] for identity, _ in captured]]
        subprocess.run(command, cwd=folder, env=env, check=True)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    action = parser.add_mutually_exclusive_group(required=True)
    action.add_argument('--run', type=Path)
    action.add_argument('--execute', type=Path, help=argparse.SUPPRESS)
    parser.add_argument('--grids', nargs='+', type=Path, required=True, metavar='GRID_JSON')
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--controllers-sha256', nargs=3, help=argparse.SUPPRESS)
    args = parser.parse_args(argv)
    if args.run:
        launch(args.run, args.grids, args.output)
        return 0
    if args.controllers_sha256 is None:
        parser.error('Captured --execute requires all three --controllers-sha256 values')
    report = execute(args.execute, args.grids, args.output,
                     expected_controller_hashes=args.controllers_sha256)
    print(json.dumps({'audit_passed': report['audit_passed'],
        'thermal_ensemble_converged': report['thermal_ensemble_converged'],
        'missing_grid_count': len(report['coverage']['missing_grids']),
        'new_solve_count': 0, 'report': str(args.output)}), flush=True)
    return 0 if report['audit_passed'] else 1


if __name__ == '__main__':
    raise SystemExit(main())
