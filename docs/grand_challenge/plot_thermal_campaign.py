"""Render immutable path errors without importing the atomic solver.

python -m docs.grand_challenge.plot_thermal_campaign REPORT.json NEW.png
"""

import argparse
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np


def render(report_path, output_path):
    output_path = Path(output_path)
    if output_path.exists():
        raise FileExistsError(output_path)
    report = json.loads(Path(report_path).read_text(encoding='utf-8'))
    payload = {k: v for k, v in report.items() if k != 'record_sha256'}
    encoded = json.dumps(payload, sort_keys=True, separators=(',', ':'),
                         ensure_ascii=False, allow_nan=False).encode('utf-8')
    if hashlib.sha256(encoded).hexdigest() != report['record_sha256']:
        raise ValueError('Path report content hash differs')
    keys = ('greater', 'lesser', 'greater_by_source', 'lesser_by_source',
            'mean_pulse', 'mean_outer', 'retarded_response', 'exit_state')
    labels = ('Greater', 'Lesser', 'Greater\nby source', 'Lesser\nby source',
              'Mean', 'Mean outer', 'Response', 'Exit state')
    ledger = report['ledger']
    comparisons = [
        ('CF4 coarse to middle', ledger['primary_refinement'][0], 1e-3),
        ('CF4 middle to fine', ledger['primary_refinement'][1], 1e-3),
        ('Fine CF4 vs fine adjoint', ledger['independent_reference'], 5e-6),
        ('Independent adjoint refinement', ledger['independent_refinement'], 2e-6)]
    fig, ax = plt.subplots(figsize=(11.5, 4.7), layout='constrained')
    for label, values, budget in comparisons:
        ratios = np.array([values[k]/budget for k in keys])
        if not np.isfinite(ratios).all() or np.any(ratios <= 0):
            raise ValueError('This log plot requires finite positive measured errors')
        ax.semilogy(range(len(keys)), ratios, 'o-', label=label, markersize=4)
    ax.axhline(1, color='#aa3434', ls='--', label='Acceptance limit')
    ax.set_xticks(range(len(keys)), labels)
    ax.set_ylabel('Maximum RF/source error divided by its declared budget')
    ax.set_title('New thermal Rb path: p1 / seed11 / index1')
    ax.grid(alpha=.2, which='both')
    ax.legend(fontsize=8, ncols=2)
    fig.supxlabel('Atomic path validation only; full thermal ensemble and optical squeezing NOT CERTIFIED', fontsize=9)
    with output_path.open('xb') as handle:
        fig.savefig(handle, format='png', dpi=160)
    plt.close(fig)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('report', type=Path)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    render(args.report, args.output)
