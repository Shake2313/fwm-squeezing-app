"""Complete source/response coverage for the hash-verified selected Rb path."""

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import asdict
import hashlib
import json
from pathlib import Path
import platform

import numpy as np
import scipy

from gabes import core
from gabes.fwm_quantum import transport_ensemble as stream
from gabes.fwm_quantum.smooth_transport import smooth_rb_problem
from .smooth_transport_audit import ROOT, fixture, encode
from .reference.adjoint_transport import adjoint_wavepacket


PARENT = ROOT/'docs/grand_challenge/smooth_transport_report_v2.json'
METRICS = stream.PATH_METRICS+('exit_state',)
BUDGETS = {'primary_refinement': 1e-3, 'independent_reference': 5e-6,
           'independent_refinement': 2e-6}
REFERENCE_TOLERANCES = ((1e-9, 1e-12), (3e-10, 3e-13))


def decode(value):
    if isinstance(value, dict):
        if set(value) == {'real', 'imag'}:
            return np.asarray(value['real'])+1j*np.asarray(value['imag'])
        return {k: decode(v) for k, v in value.items()}
    if isinstance(value, list):
        return [decode(v) for v in value]
    return value


def file_hash(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def parent_report():
    report = json.loads(PARENT.read_text(encoding='utf-8'))
    if (report.get('schema') != 'gabes-continuous-rb-transport-v2'
            or not report.get('all_declared_controls_passed')
            or not report.get('source_stable_during_run')):
        raise ValueError('verified continuous selected-path parent required')
    if any(file_hash(ROOT/name) != digest for name, digest in report['source_sha256'].items()):
        raise ValueError('parent source changed: rerun the primary/reference audit')
    inputs, geometry, path, axis = fixture()
    if (report['inputs_SI'] != asdict(inputs)
            or report['path']['entry_position_m'] != path.entry_position_m.tolist()
            or report['path']['velocity_m_s'] != path.velocity_m_s.tolist()
            or report['path']['residence_time_s'] != path.residence_time_s
            or report['analysis_frequencies_hz'] != axis.frequency_hz.tolist()
            or report['wavevectors_rad_m'] != geometry.wavevectors_rad_m.tolist()):
        raise ValueError('parent physical fixture differs')
    return report


def hashes():
    names = list(parent_report()['source_sha256'])+[
        'analysis/grand_challenge/reference/adjoint_transport.py',
        'analysis/grand_challenge/adjoint_transport_audit.py',
        'gabes/fwm_quantum/transport_ensemble.py', 'gabes/quantum/inflow.py',
        'tests/quantum/test_adjoint_transport.py', 'tests/quantum/test_adjoint_transport_audit.py']
    return {name: file_hash(ROOT/name) for name in names}


def problem():
    inputs, geometry, path, axis = fixture()
    p = smooth_rb_problem(inputs, geometry, path, axis)
    q = geometry.wavevectors_rad_m[1:]-geometry.wavevectors_rad_m[0]
    _, bp, bc = geometry.atomic_frequencies(inputs, np.zeros((1, 3)))
    convention = stream.PumpOnlyConvention(
        analysis_axis=axis, number_density_m3=inputs.number_density_m3,
        coupling_scales=np.tile(p['metadata']['coupling_s_inverse_sqrt_flux'], 2),
        readout_scales=np.linalg.norm(p['readouts'], axis=(-2, -1)),
        source_names=('atomic_inflow',)+tuple('jump:'+c.name for c in p['reservoirs'].channels),
        carrier_offsets_at_rest_rad_s=np.array([bp[0], bc[0], -bp[0], -bc[0]]),
        port_wavevectors_rad_m=np.concatenate((q, -q)),
        model_id='reduced-Rb-continuous-Gaussian-selected-path-v1', model_scope='smooth',
        provenance='same prescribed 2 us Rb path as immutable smooth_transport_report_v2; density unused by path',
        phase_provenance='one common lab beat; actual geometric entry offsets retained')
    return p, convention, path


def primary_packets(parent, p, convention, path):
    packets = []
    # All original source hashes/inputs are checked before this reuse. Same
    # primary mathematical solutions need no repeated solve. Their NEW source-
    # and RF-resolved error norms are recomputed below, not copied from old gates.
    for j in range(3):
        row = parent['cases'][f'primary_{j}']
        packet = decode(row['values'])
        if packet['metadata'] != p['metadata']:
            # JSON preserves tuples as lists; compare the canonical encoding.
            if json.dumps(encode(packet['metadata']), sort_keys=True) != json.dumps(encode(p['metadata']), sort_keys=True):
                raise ValueError('parent Rb metadata differs')
        packet.update(analysis_axis=convention.analysis_axis,
                      frequencies_rad_s=p['frequencies_rad_s'], residence_time_s=path.residence_time_s)
        packet['mean_outer'] = packet['mean_pulse'][:, :, None]*packet['mean_pulse'][:, None, :].conj()
        packets.append(packet)
    return packets


def errors(candidate, reference, p):
    """Each RF and each named source has its own norm; dark floors use SI units."""
    t = p['duration_s']
    op = float(np.linalg.norm(p['readouts'], axis=(-2, -1)).max())
    drive = float(np.linalg.norm(p['drives'], axis=(-2, -1)).max())
    eps = 128*np.finfo(float).eps
    result = {}
    for key in METRICS:
        a, b = np.asarray(candidate[key]), np.asarray(reference[key])
        axes = (-1,) if key == 'mean_pulse' else (-2, -1)
        floor = eps if key == 'exit_state' else eps*t*op if key == 'mean_pulse' else eps*t*t*op*drive if key == 'retarded_response' else eps*t*t*op*op
        result[key] = float(np.max(np.linalg.norm(a-b, axis=axes)/np.maximum(np.linalg.norm(b, axis=axes), floor)))
    return result


def evidence_for(packets, reference, p, convention, path):
    ids = [stream.packet_digest(path, row, convention) for row in packets]
    controls = [1/row['numerics']['rtol'] for row in packets]
    comparisons = []
    for j in range(2):
        comparisons.append(stream.NumericalComparison(
            'path_refinement', errors(packets[j+1], packets[j], p),
            dict.fromkeys(METRICS, BUDGETS['primary_refinement']), ids[j], ids[j+1],
            (controls[j], controls[j+1]),
            'max relative Frobenius norm per RF AND named source; 128*eps SI dimensional dark floor; mean row uses Euclidean norm',
            'two successive refinements of hash-verified continuous microscopic-D primary results'))
    reference_id = hashlib.sha256(json.dumps(encode({k: reference[k] for k in METRICS}), sort_keys=True).encode()).hexdigest()
    comparisons.append(stream.NumericalComparison(
        'independent_reference', errors(packets[-1], reference, p),
        dict.fromkeys(METRICS, BUDGETS['independent_reference']), reference_id, ids[-1],
        (0., controls[-1]), comparisons[0].error_definition,
        'independently solved full-density/backward-observable Lindblad-product source noise and causal response'))
    evidence = stream.ConvergenceEvidence(ids[-1], 'smooth', tuple(comparisons),
                                          'selected path only; references and parent SHA256 retained in audit')
    return evidence, ids


def worker(pair):
    before = hashes()
    p, _, _ = problem()
    with core.blas_single_thread():
        result = adjoint_wavepacket(**{k: v for k, v in p.items() if k != 'metadata'}, rtol=pair[0], atol=pair[1])
    return {'rtol': pair[0], 'atol': pair[1], 'source_stable_during_run': before == hashes(), 'values': result}


def write_audit(output, plot):
    output, plot = Path(output), Path(plot)
    if output.resolve() == plot.resolve() or output.exists() or plot.exists():
        raise FileExistsError('new distinct immutable report and plot paths required')
    before, parent, parent_hash = hashes(), parent_report(), file_hash(PARENT)
    p, convention, path = problem()
    packets = primary_packets(parent, p, convention, path)
    cases = {}
    with ProcessPoolExecutor(max_workers=2) as pool:
        futures = {pool.submit(worker, pair): j for j, pair in enumerate(REFERENCE_TOLERANCES)}
        for future in as_completed(futures):
            row = future.result()
            cases[f'adjoint_{futures[future]}'] = row
            print(f"adjoint_{futures[future]} finished: {row['values']['numerics']['elapsed_seconds']:.2f}s", flush=True)
    fine, coarse = cases['adjoint_1']['values'], cases['adjoint_0']['values']
    refinement = errors(fine, coarse, p)
    evidence, ids = evidence_for(packets, fine, p, convention, path)
    reasons = stream._evidence_reasons(evidence, ids[-1], 'smooth')
    # Using total QRT alone must still fail the stronger source/response gate.
    incomplete = stream.NumericalComparison('independent_reference', {'greater': 0., 'lesser': 0.},
        {'greater': 5e-6, 'lesser': 5e-6}, 'total-QRT-only', ids[-1], (0., 1.),
        'deliberately incomplete evidence', 'negative coverage control')
    bad_evidence = stream.ConvergenceEvidence(ids[-1], 'smooth', evidence.comparisons[:-1]+(incomplete,), 'negative control')
    bad_reasons = stream._evidence_reasons(bad_evidence, ids[-1], 'smooth')
    altered = dict(packets[-1])
    altered['retarded_response'] = -altered['retarded_response']
    altered_id = stream.packet_digest(path, altered, convention)
    stale_reasons = stream._evidence_reasons(evidence, altered_id, 'smooth')
    stable = before == hashes() and parent_hash == file_hash(PARENT) and all(row['source_stable_during_run'] for row in cases.values())
    passed = (stable and not reasons and max(refinement.values()) < BUDGETS['independent_refinement']
              and bool(bad_reasons) and bool(stale_reasons))
    report = {'schema': 'gabes-adjoint-source-response-audit-v1', 'source_sha256': before,
        'source_stable_during_run': stable, 'parent_report': {'path': str(PARENT.relative_to(ROOT)), 'sha256': parent_hash,
            'reused_primary_cases': ['primary_0', 'primary_1', 'primary_2'], 'all_parent_source_hashes_match': True},
        'environment': {'python': platform.python_version(), 'numpy': np.__version__, 'scipy': scipy.__version__},
        'budgets': BUDGETS, 'reference_tolerances': REFERENCE_TOLERANCES,
        'path': parent['path'], 'inputs_SI': parent['inputs_SI'], 'frequency_hz': parent['analysis_frequencies_hz'],
        'source_names': convention.source_names, 'primary_packet_digests': ids,
        'comparisons': [{'kind': c.kind, 'errors': dict(c.errors), 'tolerances': dict(c.tolerances),
            'reference_id': c.reference_id, 'candidate_id': c.candidate_id, 'control_values': c.control_values,
            'error_definition': c.error_definition, 'provenance': c.provenance, 'passed': c.passed} for c in evidence.comparisons],
        'independent_refinement': refinement, 'complete_path_evidence_accepted': not reasons,
        'path_evidence_reasons': reasons, 'negative_controls': {'total_QRT_only_rejected': bool(bad_reasons),
            'total_QRT_only_reasons': bad_reasons, 'modified_response_digest_rejected': bool(stale_reasons)},
        'cases': encode(cases), 'all_declared_controls_passed': bool(passed),
        'scope': {'physical_Rb_selected_path_source_and_response_reference': True,
            'complete_selected_path_evidence': bool(passed), 'thermal_Rb_ensemble_converged': False,
            'nonlocal_Maxwell_field': False, 'finite_seed': False, 'absolute_squeezing': False},
        'limitations': ['Observed numerical convergence, not a rigorous global error bound.',
            'All physical jumps remain active in the shared density and adjoint; attribution depends on this fixed jump representation.',
            'The selected path does not certify uncomputed thermal paths or actual thermal integrand convergence.',
            'Retarded response is an integrated atomic response, not the spatially nonlocal Maxwell M.']}
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(11.5, 4.8), layout='constrained')
    x = np.arange(len(METRICS))
    labels = [key.replace('_by_source', '\nby source').replace('_', ' ') for key in METRICS]
    for row in report['comparisons']:
        ax.semilogy(x, [max(row['errors'][k], 1e-16) for k in METRICS], 'o-', label=row['kind']+' '+str(len(ax.lines)+1))
    ax.semilogy(x, [max(refinement[k], 1e-16) for k in METRICS], 's--', label='adjoint refinement')
    for label, limit in BUDGETS.items():
        ax.axhline(limit, ls=':', label=label+' budget')
    ax.set(xticks=x, xticklabels=labels, ylabel='Max relative difference per RF/source',
           title='Selected 2 us Rb source and response: '+('PASS' if passed else 'FAIL')+'; no ensemble/optical claim')
    ax.legend(fontsize=8, ncols=2)
    output.parent.mkdir(parents=True, exist_ok=True)
    plot.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(plot, dpi=160)
    plt.close(fig)
    with output.open('x', encoding='utf-8') as handle:
        json.dump(encode(report), handle, ensure_ascii=False, indent=2, allow_nan=False)
        handle.write('\n')
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True)
    parser.add_argument('--plot', required=True)
    args = parser.parse_args()
    report = write_audit(args.output, args.plot)
    print(json.dumps({k: report[k] for k in ('independent_refinement', 'complete_path_evidence_accepted', 'all_declared_controls_passed')}, indent=2))
    if not report['all_declared_controls_passed']:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
