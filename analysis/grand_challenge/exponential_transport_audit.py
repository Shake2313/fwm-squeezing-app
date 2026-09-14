"""Observed CF4 convergence against immutable, independent Rb source evidence."""

import argparse
import hashlib
import json
import platform
import time
from pathlib import Path

import numpy as np
import scipy
from threadpoolctl import threadpool_info, threadpool_limits

from gabes.fwm_quantum import transport_ensemble as stream
from gabes.quantum.segmented_transport import segmented_wavepacket
from . import adjoint_transport_audit as adj
from .reference.exponential_transport import exponential_wavepacket
from .smooth_transport_audit import ROOT, encode


REFERENCE = ROOT/'docs/grand_challenge/adjoint_transport_report_v1.json'
DEFAULT_SEGMENTS = (4096, 8192, 16384)


def hashes():
    result = adj.hashes()
    for name in ('analysis/grand_challenge/reference/exponential_transport.py',
                 'analysis/grand_challenge/exponential_transport_audit.py',
                 'tests/quantum/test_exponential_transport.py',
                 'tests/quantum/test_exponential_transport_audit.py'):
        result[name] = adj.file_hash(ROOT/name)
    return result


def reference_report():
    report = json.loads(REFERENCE.read_text(encoding='utf-8'))
    parent = adj.parent_report()
    if (report.get('schema') != 'gabes-adjoint-source-response-audit-v1'
            or not report.get('all_declared_controls_passed')
            or not report.get('source_stable_during_run')
            or not report.get('complete_path_evidence_accepted')
            or report['source_sha256'] != adj.hashes()
            or report['parent_report']['sha256'] != adj.file_hash(adj.PARENT)
            or report['path'] != parent['path'] or report['inputs_SI'] != parent['inputs_SI']
            or report['frequency_hz'] != parent['analysis_frequencies_hz']):
        raise ValueError('unchanged, complete selected-path adjoint reference required')
    return report


def validate_segments(segments):
    if (len(segments) != 3 or any(isinstance(n, bool) or not isinstance(n, (int, np.integer))
                                or n < 1 for n in segments)
            or any(b <= a for a, b in zip(segments, segments[1:]))):
        raise ValueError('three strictly increasing positive integer segment counts required')
    return tuple(int(n) for n in segments)


def evidence_for(packets, segments, reference, p, convention, path):
    segments = validate_segments(segments)
    if len(packets) != 3:
        raise ValueError('three independent primary refinements required')
    ids = [stream.packet_digest(path, row, convention) for row in packets]
    comparisons = []
    definition = 'Maximum relative norm per RF AND named source; same 128*eps SI dark floors as immutable adjoint audit'
    for j in range(2):
        comparisons.append(stream.NumericalComparison('path_refinement',
            adj.errors(packets[j+1], packets[j], p), dict.fromkeys(adj.METRICS, adj.BUDGETS['primary_refinement']),
            ids[j], ids[j+1], segments[j:j+2], definition,
            'same affine Lindblad moment lift, positive-duration CF4 substeps; no fitted output adjustment'))
    reference_id = hashlib.sha256(json.dumps(encode({k: reference[k] for k in adj.METRICS}), sort_keys=True).encode()).hexdigest()
    comparisons.append(stream.NumericalComparison('independent_reference',
        adj.errors(packets[-1], reference, p), dict.fromkeys(adj.METRICS, adj.BUDGETS['independent_reference']),
        reference_id, ids[-1], (0., segments[-1]), definition,
        'full-density backward-observable source noise and response; immutable reference hashes verified before reuse'))
    evidence = stream.ConvergenceEvidence(ids[-1], 'smooth', tuple(comparisons),
        'Observed convergence for selected Rb path only; eigenbasis screens are not global error bounds')
    return evidence, ids


def fixed_map_benchmark(p, segments=16):
    """Compare identical midpoint substeps on this process/hardware/BLAS setup."""
    edges = np.linspace(0, p['duration_s'], segments+1)
    h = [p['h0']+p['envelope']((a+b)/2)*p['h1'] for a, b in zip(edges[:-1], edges[1:])]
    begin = time.monotonic()
    original = segmented_wavepacket(h, p['reservoirs'], p['boundary_state'], np.diff(edges),
        np.broadcast_to(p['readouts'], (segments,)+p['readouts'].shape), p['frequencies_rad_s'],
        drives=np.broadcast_to(p['drives'], (segments,)+p['drives'].shape),
        drive_frequencies_rad_s=p['drive_frequencies_rad_s'])
    original_seconds = time.monotonic()-begin
    original['mean_outer'] = original['mean_pulse'][:, :, None]*original['mean_pulse'][:, None, :].conj()
    candidate = exponential_wavepacket(**p, segments=segments, order=2)
    difference = adj.errors(candidate, original, p)
    return {'segments': segments, 'order': 2, 'errors': difference, 'tolerance': 2e-8,
        'passed': max(difference.values()) < 2e-8, 'original_seconds': original_seconds,
        'spectral_seconds': candidate['numerics']['elapsed_seconds'],
        'same_discretization_speed_ratio': original_seconds/candidate['numerics']['elapsed_seconds'],
        'scope': 'Identical coarse midpoint maps only. This does not certify continuum accuracy or an Ultra speedup.'}


def write_audit(output, segments=DEFAULT_SEGMENTS):
    segments = validate_segments(segments)
    output = Path(output)
    if output.exists():
        raise FileExistsError('new immutable output path required')
    before, ref_report = hashes(), reference_report()
    reference_hash = adj.file_hash(REFERENCE)
    p, convention, path = adj.problem()
    metadata = p.pop('metadata')
    reference = adj.decode(ref_report['cases']['adjoint_1']['values'])
    packets = []
    with threadpool_limits(limits=1, user_api='blas'):
        pools = threadpool_info()
        if not pools or any(row['num_threads'] != 1 for row in pools if row['user_api'] == 'blas'):
            raise RuntimeError('verified single-thread BLAS required for benchmark')
        fixed = fixed_map_benchmark(p)
        print('fixed-map comparison: '+json.dumps(fixed), flush=True)
        for n in segments:
            row = exponential_wavepacket(**p, segments=n)
            row['analysis_axis'] = convention.analysis_axis
            row['metadata'] = metadata
            packets.append(row)
            print(json.dumps({'segments': n, 'seconds': row['numerics']['elapsed_seconds'],
                              'independent_errors': adj.errors(row, reference, p)}), flush=True)
    evidence, ids = evidence_for(packets, segments, reference, p, convention, path)
    reasons = stream._evidence_reasons(evidence, ids[-1], 'smooth')
    stable = before == hashes() and reference_hash == adj.file_hash(REFERENCE)
    # Changing just one complex response must invalidate the bound evidence.
    altered = dict(packets[-1]); altered['retarded_response'] = -altered['retarded_response']
    rejects_changed_response = bool(stream._evidence_reasons(evidence, stream.packet_digest(path, altered, convention), 'smooth'))
    passed = stable and fixed['passed'] and not reasons and rejects_changed_response
    report = {'schema': 'gabes-exponential-source-response-audit-v1', 'source_sha256': before,
        'source_stable_during_run': stable, 'reference': {'path': str(REFERENCE.relative_to(ROOT)),
            'sha256': reference_hash, 'reused_cases': ['adjoint_0', 'adjoint_1'],
            'independent_refinement': ref_report['independent_refinement']},
        'segments': segments, 'budgets': adj.BUDGETS, 'fixed_map_benchmark': fixed,
        'environment': {'python': platform.python_version(), 'numpy': np.__version__, 'scipy': scipy.__version__,
                        'threadpools_inside_limit': pools},
        'path': ref_report['path'], 'inputs_SI': ref_report['inputs_SI'], 'frequency_hz': ref_report['frequency_hz'],
        'comparisons': [{'kind': c.kind, 'errors': dict(c.errors), 'tolerances': dict(c.tolerances),
            'reference_id': c.reference_id, 'candidate_id': c.candidate_id, 'control_values': c.control_values,
            'error_definition': c.error_definition, 'provenance': c.provenance, 'passed': c.passed}
            for c in evidence.comparisons], 'packet_digests': ids,
        'complete_path_evidence_accepted': not reasons, 'path_evidence_reasons': reasons,
        'modified_response_digest_rejected': rejects_changed_response,
        'cases': {str(n): {k: v for k, v in row.items() if k != 'analysis_axis'} for n, row in zip(segments, packets)},
        'all_declared_controls_passed': bool(passed),
        'scope': {'complete_selected_path_evidence': bool(passed), 'physical_Rb_ensemble_converged': False,
                  'nonlocal_Maxwell_field': False, 'finite_seed': False, 'absolute_squeezing': False},
        'limitations': ['CF4 convergence is observed, not a rigorous global error bound.',
            'Spectral and block formulas are algebraically identical for the same prescribed substeps.',
            'Old adaptive report wall times have no verified BLAS setup; no comparable continuum speed ratio claimed.',
            'Fixed-map benchmark is neither the public Ultra scheme nor a converged thermal ensemble.']}
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open('x', encoding='utf-8') as handle:
        json.dump(encode(report), handle, ensure_ascii=False, indent=2, allow_nan=False)
        handle.write('\n')
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True)
    parser.add_argument('--segments', nargs=3, type=int, default=DEFAULT_SEGMENTS)
    args = parser.parse_args()
    result = write_audit(args.output, args.segments)
    print('all declared controls passed:', result['all_declared_controls_passed'])
    if not result['all_declared_controls_passed']:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
