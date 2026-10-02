"""Plot sealed grid diagnostics without importing or running atomic solvers.

python -m docs.grand_challenge.plot_thermal_grid P0.json P1.json NEW.png
Use --cross-audit AUDIT.json for a verified parent-to-child campaign comparison.
"""

import argparse
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np


METRICS = ('greater', 'lesser', 'greater_by_source', 'lesser_by_source',
           'poisson_number', 'retarded_response')


def read_sealed(path):
    def no_duplicates(items):
        result = {}
        for key, value in items:
            if key in result:
                raise ValueError('Duplicate JSON key')
            result[key] = value
        return result
    report = json.loads(Path(path).read_text(encoding='utf-8'), object_pairs_hook=no_duplicates)
    payload = {k: v for k, v in report.items() if k != 'record_sha256'}
    encoded = json.dumps(payload, sort_keys=True, separators=(',', ':'),
                         ensure_ascii=False, allow_nan=False).encode('utf-8')
    if hashlib.sha256(encoded).hexdigest() != report['record_sha256']:
        raise ValueError('Report content hash differs')
    return report


def read_grid(path):
    report = read_sealed(path)
    if (report['schema'] != 'gabes-frozen-thermal-grid-audit-v1'
            or not report['grid_passed'] or report['row']['candidate']['spectra'] is None):
        raise ValueError('Complete independently audited atomic grid required')
    return report


def array(value):
    if isinstance(value, dict) and set(value) == {'real', 'imag'}:
        return np.asarray(value['real'])+1j*np.asarray(value['imag'])
    return np.asarray(value)


def cross_binding(coarse, fine, audit):
    """Check an existing evidence audit; a plot never grants scientific approval."""
    if (audit.get('schema') != 'gabes-frozen-thermal-cross-grid-comparison-v1'
            or audit.get('audit_passed') is not True
            or audit.get('new_solve_count') != 0
            or audit.get('thermal_ensemble_converged') is not False
            or audit.get('physical_optical_prediction') is not False
            or audit.get('source_identity') != coarse['source_identity']
            or audit.get('source_bundle_sha256') != coarse['source_bundle_sha256']
            or audit.get('nested_reuse', {}).get('passed') is not True
            or audit.get('nested_sum', {}).get('passed') is not True):
        raise ValueError('Passing sealed cross-campaign evidence audit required')
    references = audit.get('input_reports', [])
    if len(references) != 2:
        raise ValueError('Cross-campaign audit must bind both input grids')
    for report, reference in zip((coarse, fine), references):
        if any(reference.get(key) != report[key]
               for key in ('record_sha256', 'campaign_sha256', 'selection')):
            raise ValueError('Cross-campaign audit is bound to different input grids')


def frequency_comparison(coarse, fine, *, cross_audit=None):
    if (coarse['source_identity'] != fine['source_identity']
            or coarse['source_bundle_sha256'] != fine['source_bundle_sha256']
            or coarse['selection'][1] != fine['selection'][1]
            or coarse['selection'][0]+1 != fine['selection'][0]
            or any(coarse['row'][key] != fine['row'][key] for key in ('model_digest', 'plan_digest'))):
        raise ValueError('Consecutive grids of one seed, campaign, model and path plan required')
    if coarse['campaign_sha256'] != fine['campaign_sha256'] and cross_audit is None:
        raise ValueError('Different declarations require a sealed cross-campaign evidence audit')
    if cross_audit is not None:
        cross_binding(coarse, fine, cross_audit)
    one, two = coarse['row']['candidate'], fine['row']['candidate']
    if any(one[key] != two[key] for key in ('analysis_axis', 'source_names')):
        raise ValueError('Matching RF and named-source ordering required')
    omega = np.asarray(one['analysis_axis']['analysis_axis_omega_rad_s'])
    if omega.ndim != 1 or not len(omega) or not np.isfinite(omega).all():
        raise ValueError('Finite one-dimensional analysis-frequency axis required')
    errors = {}
    for key in METRICS:
        a = array(fine['row']['candidate']['spectra'][key])
        b = array(coarse['row']['candidate']['spectra'][key])
        if a.shape != b.shape or not np.isfinite(a).all() or not np.isfinite(b).all():
            raise ValueError('Matching finite grid matrices required')
        scale = coarse['row']['comparison_scales'][key]
        if not np.isfinite(scale) or scale <= 0:
            raise ValueError('Positive dimensional comparison scale required')
        floor = 128*np.finfo(float).eps*scale
        values = np.linalg.norm(a-b, axis=(-2, -1)) / np.maximum(
            np.linalg.norm(b, axis=(-2, -1)), floor)
        expected = (len(one['source_names']), len(omega)) if key.endswith('_by_source') else (len(omega),)
        if values.shape != expected or not np.isfinite(values).all():
            raise ValueError('Complete finite RF/source relative changes required')
        # A source maximum is taken separately at each RF; no averaging can
        # hide a failing source. This is a change diagnostic, not an error bound.
        errors[key] = np.max(values, axis=0) if values.ndim == 2 else values
    return errors


def comparison(coarse, fine, *, cross_audit=None):
    return {key: float(np.max(values)) for key, values in
            frequency_comparison(coarse, fine, cross_audit=cross_audit).items()}


def render(coarse_path, fine_path, output_path, *, cross_audit_path=None):
    output_path = Path(output_path)
    if output_path.exists():
        raise FileExistsError(output_path)
    coarse, fine = read_grid(coarse_path), read_grid(fine_path)
    audit = read_sealed(cross_audit_path) if cross_audit_path is not None else None
    per_frequency = frequency_comparison(coarse, fine, cross_audit=audit)
    errors = {key: float(np.max(values)) for key, values in per_frequency.items()}
    if audit is not None:
        for filename, reference in zip((coarse_path, fine_path), audit['input_reports']):
            if hashlib.sha256(Path(filename).read_bytes()).hexdigest() != reference['file_sha256']:
                raise ValueError('Input grid bytes differ from the cross-campaign audit')
        if errors != audit['refinement_diagnostic']['errors']:
            raise ValueError('Plotted changes differ from the cross-campaign audit')
    frequencies_mhz = np.asarray(coarse['row']['candidate']['analysis_axis']['analysis_axis_omega_rad_s'])/(2*np.pi*1e6)
    rows = fine['row']['path_ledgers']
    labels = ['CF4 coarse to middle', 'CF4 middle to fine',
              'Fine CF4 vs fine adjoint', 'Independent adjoint refinement']
    ratios = np.array([[max(item['primary_refinement'][0].values())/1e-3,
                        max(item['primary_refinement'][1].values())/1e-3,
                        max(item['independent_reference'].values())/5e-6,
                        max(item['independent_refinement'].values())/2e-6] for item in rows])
    fig, axes = plt.subplots(1, 3, figsize=(17.2, 4.8), layout='constrained')
    for j, label in enumerate(labels):
        axes[0].semilogy(np.arange(len(rows)), ratios[:, j], 'o-', markersize=3, label=label)
    axes[0].axhline(1, color='#ad3333', ls='--', label='Path acceptance limit')
    axes[0].set(xlabel=f'Declared p{fine["selection"][0]} boundary path index', ylabel='Error / path budget',
                title='All RF frequencies, sources and eight path metrics')
    # Preserve every plotted path while keeping labels legible on larger grids.
    tick_step = max(1, (len(rows) + 15) // 16)
    ticks = sorted(set(range(0, len(rows), tick_step)) | {len(rows) - 1})
    axes[0].set_xticks(ticks)
    axes[0].legend(fontsize=7)
    axes[0].grid(alpha=.2)
    names = ['Greater', 'Lesser', 'Greater\nby source', 'Lesser\nby source', 'Poisson\nnumber', 'Response']
    axes[1].bar(np.arange(len(METRICS)), [errors[key] for key in METRICS], color='#397c8a')
    axes[1].axhline(.05, color='#ad3333', ls='--', label='Declared ensemble budget (5%)')
    axes[1].set_xticks(np.arange(len(METRICS)), names, fontsize=8)
    axes[1].set(ylabel='Maximum RF/source relative matrix change',
                title=f'p{coarse["selection"][0]} to p{fine["selection"][0]}, seed{fine["selection"][1]}: grid diagnostic')
    axes[1].legend(fontsize=8)
    axes[1].grid(axis='y', alpha=.2)
    for key, name in zip(METRICS, names):
        axes[2].plot(np.arange(len(frequencies_mhz)), per_frequency[key], 'o-',
                     markersize=4, label=name.replace('\n', ' '))
    axes[2].axhline(.05, color='#ad3333', ls='--', label='5% budget')
    axes[2].set_xticks(np.arange(len(frequencies_mhz)), [f'{value:g}' for value in frequencies_mhz])
    axes[2].set(xlabel='Analysis frequency (MHz)', ylabel='Relative matrix change',
                title='Frequency-resolved change; worst source at each RF')
    axes[2].legend(fontsize=7)
    axes[2].grid(alpha=.2)
    fig.supxlabel('Conditional pump-only atomic stream; ensemble convergence and optical squeezing NOT CERTIFIED', fontsize=10)
    with output_path.open('xb') as handle:
        fig.savefig(handle, format='png', dpi=160)
    plt.close(fig)
    print(json.dumps({'coarse_record_sha256': coarse['record_sha256'],
                      'fine_record_sha256': fine['record_sha256'],
                      'cross_audit_record_sha256': audit['record_sha256'] if audit else None,
                      'relative_grid_changes': errors,
                      'analysis_frequency_mhz': frequencies_mhz.tolist(),
                      'relative_changes_by_frequency': {key: value.tolist() for key, value in per_frequency.items()},
                      'thermal_ensemble_converged': False}, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('coarse', type=Path)
    parser.add_argument('fine', type=Path)
    parser.add_argument('output', type=Path)
    parser.add_argument('--cross-audit', type=Path)
    args = parser.parse_args()
    render(args.coarse, args.fine, args.output, cross_audit_path=args.cross_audit)
