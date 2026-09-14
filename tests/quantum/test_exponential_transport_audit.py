"""Convergence evidence must remain bound to every source and response."""

import json
import numpy as np
import pytest

from gabes.fwm_quantum import transport_ensemble as stream
from analysis.grand_challenge import exponential_transport_audit as audit


def data():
    report = audit.reference_report()
    p, convention, path = audit.adj.problem()
    reference = audit.adj.decode(report['cases']['adjoint_1']['values'])
    row = {**reference, 'analysis_axis': convention.analysis_axis, 'metadata': p['metadata']}
    return p, convention, path, reference, row


def test_matching_evidence_and_changed_response_digest():
    p, convention, path, reference, row = data()
    evidence, ids = audit.evidence_for([row, row, row], (4, 8, 16), reference, p, convention, path)
    assert not stream._evidence_reasons(evidence, ids[-1], 'smooth')
    changed = {**row, 'retarded_response': -row['retarded_response']}
    digest = stream.packet_digest(path, changed, convention)
    assert stream._evidence_reasons(evidence, digest, 'smooth')
    assert all(set(c.errors) == set(audit.adj.METRICS) for c in evidence.comparisons)


@pytest.mark.parametrize('metric', audit.adj.METRICS)
def test_each_metric_can_prevent_convergence(metric):
    p, convention, path, reference, row = data()
    # Same change in all refinements: refinement passes, independent test fails.
    changed = {**row, metric: np.asarray(row[metric])*1.02}
    if metric in ('greater', 'lesser'):
        changed[metric+'_by_source'] = row[metric+'_by_source']*1.02
    elif metric.endswith('_by_source'):
        changed[metric.removesuffix('_by_source')] = changed[metric].sum(axis=0)
    if metric == 'exit_state':
        # Exit state is an additional audit gate, outside the seven stream inputs.
        changed[metric] = .98*row[metric]+.02*np.eye(4)/4
    evidence, ids = audit.evidence_for([changed]*3, (4, 8, 16), reference, p, convention, path)
    assert not evidence.comparisons[-1].passed
    assert stream._evidence_reasons(evidence, ids[-1], 'smooth')


@pytest.mark.parametrize('segments', [(1, 2), (1, 2, 3, 4), (1, 1, 2), (3, 2, 1),
                                     (0, 2, 4), (True, 2, 4), (1, 2., 4)])
def test_incomplete_or_nonmonotone_refinement_controls_rejected(segments):
    with pytest.raises(ValueError, match='three strictly increasing'):
        audit.validate_segments(segments)


def test_changed_reference_source_cannot_be_reused(tmp_path, monkeypatch):
    report = audit.reference_report()
    report['source_sha256']['analysis/grand_challenge/reference/adjoint_transport.py'] = '0'*64
    path = tmp_path/'changed.json'
    path.write_text(json.dumps(report), encoding='utf-8')
    monkeypatch.setattr(audit, 'REFERENCE', path)
    with pytest.raises(ValueError, match='unchanged, complete'):
        audit.reference_report()


def test_existing_report_is_never_overwritten(tmp_path):
    path = tmp_path/'existing.json'
    path.write_text('preserve', encoding='utf-8')
    with pytest.raises(FileExistsError):
        audit.write_audit(path)
    assert path.read_text(encoding='utf-8') == 'preserve'
