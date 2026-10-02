"""Render sealed independent-scramble diagnostics; never run atomic solvers.

python -m docs.grand_challenge.plot_thermal_scrambles ENSEMBLE.json NEW.png
"""

import argparse
import hashlib
import itertools
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from .plot_thermal_grid import METRICS


def read_diagnostics(filename):
    def unique(items):
        result = {}
        for key, value in items:
            if key in result:
                raise ValueError('Duplicate JSON key')
            result[key] = value
        return result

    report = json.loads(Path(filename).read_text(encoding='utf-8'), object_pairs_hook=unique)
    payload = {key: value for key, value in report.items() if key != 'record_sha256'}
    encoded = json.dumps(payload, sort_keys=True, separators=(',', ':'),
                         ensure_ascii=False, allow_nan=False).encode('utf-8')
    if hashlib.sha256(encoded).hexdigest() != report['record_sha256']:
        raise ValueError('Ensemble report content hash differs')
    if (report['schema'] != 'gabes-frozen-thermal-ensemble-audit-v1'
            or report['audit_passed'] is not True):
        raise ValueError('Successful sealed ensemble evidence audit required')
    supplied = [tuple(grid) for grid in report['coverage']['provided_grids']]
    if len(set(supplied)) != len(supplied):
        raise ValueError('Duplicate supplied grid')
    expected = {(a, b) for a, b in itertools.permutations(supplied, 2)
                if a[0] == b[0] and a[1] != b[1]}
    rows = report['pairwise_diagnostics']['independent_scrambles']
    observed = set()
    for row in rows:
        pair = tuple(row['candidate_grid']), tuple(row['reference_grid'])
        if row['kind'] != 'independent_scramble' or pair not in expected or pair in observed:
            raise ValueError('Distinct directed seed pairs at one power required')
        if set(row['errors']) != set(METRICS):
            raise ValueError('All six stream metrics required')
        values = np.asarray([row['errors'][key] for key in METRICS], dtype=float)
        budget = row['budget']
        if (not np.isfinite(values).all() or np.any(values < 0)
                or not np.isfinite(budget) or budget <= 0
                or row['passed'] is not bool(np.all(values <= budget))):
            raise ValueError('Finite changes and consistent acceptance verdict required')
        observed.add(pair)
    if not expected or observed != expected:
        raise ValueError('Every available directed scramble pair required')
    return report, rows


def render(filename, output):
    output = Path(output)
    if output.exists():
        raise FileExistsError(output)
    report, rows = read_diagnostics(filename)
    powers = sorted({row['candidate_grid'][0] for row in rows})
    max_pairs = max(sum(row['candidate_grid'][0] == power for row in rows)
                    for power in powers)
    panel_width = 8.8 if max_pairs > 2 else 6.4
    fig, axes = plt.subplots(1, len(powers), figsize=(panel_width*len(powers), 4.8),
                             squeeze=False, layout='constrained')
    labels = ['Greater', 'Lesser', 'Greater\nby source', 'Lesser\nby source',
              'Poisson\nnumber', 'Response']
    for axis, power in zip(axes[0], powers):
        selected = [row for row in rows if row['candidate_grid'][0] == power]
        width = .8/len(selected)
        for j, row in enumerate(selected):
            # Both directions matter: the denominator is the REFERENCE norm,
            # bounded below by the audit's recorded dimensional dark floor.
            # Taking an average of seed spectra would hide this diagnostic.
            candidate, reference = row['candidate_grid'][1], row['reference_grid'][1]
            axis.bar(np.arange(len(METRICS))+(j-(len(selected)-1)/2)*width,
                     [100*row['errors'][key] for key in METRICS], width,
                     label=f'Seed {candidate}, reference {reference}')
        budgets = {row['budget'] for row in selected}
        if len(budgets) != 1:
            plt.close(fig)
            raise ValueError('One declared budget per comparison panel required')
        budget = budgets.pop()
        axis.axhline(100*budget, color='#ad3333', ls='--', label=f'{100*budget:g}% budget')
        axis.set_xticks(np.arange(len(METRICS)), labels, fontsize=8)
        axis.set(title=f'p{power}: independent scrambles',
                 ylabel='Maximum RF/source relative change (%)')
        axis.grid(axis='y', alpha=.2)
        if len(selected) > 2:
            # Keep every directed-pair label clear of the six metric groups.
            axis.legend(fontsize=8, loc='upper left', bbox_to_anchor=(1.02, 1),
                        borderaxespad=0)
        else:
            axis.legend(fontsize=8)
    covered = report['coverage']
    state = 'PASS' if report['thermal_ensemble_converged'] else 'NOT CERTIFIED'
    fig.supxlabel(f'Atomic stream only | {covered["provided_grid_count"]}/{covered["required_grid_count"]} grids audited'
                  f' | Thermal convergence: {state}\n'
                  'Directed changes, not continuum error bounds or optical squeezing', fontsize=10)
    with output.open('xb') as handle:
        fig.savefig(handle, format='png', dpi=160)
    plt.close(fig)
    print(json.dumps({'ensemble_record_sha256': report['record_sha256'],
                      'directed_comparisons': len(rows), 'powers': powers,
                      'new_solve_count': 0}, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('ensemble', type=Path)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    render(args.ensemble, args.output)
