"""Convergence evidence must remain bound to every source and response."""

import json
import numpy as np
import pytest

from gabes.fwm_quantum import transport_ensemble as stream
from analysis.grand_challenge import exponential_transport_audit as audit
from _grand_challenge_audit_fixtures import audit_sources, synthetic_packet, write_json


def data():
    p, convention, path = audit.adj.problem()
    reference = synthetic_packet(p, convention, path)
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


def test_temporary_reference_binds_actual_parent_and_sources(audit_sources):
    report = audit.reference_report()
    assert report['synthetic_test_fixture']
    assert report['source_sha256'] == audit.adj.hashes()
    assert report['parent_report']['sha256'] == audit.adj.file_hash(audit.adj.PARENT)


def test_changed_reference_source_cannot_be_reused(audit_sources):
    report = audit.reference_report()
    report['source_sha256']['analysis/grand_challenge/reference/adjoint_transport.py'] = '0'*64
    write_json(audit.REFERENCE, report)
    with pytest.raises(ValueError, match='unchanged, complete'):
        audit.reference_report()


@pytest.mark.parametrize('damage', ['changed', 'missing'])
def test_reference_rechecks_real_dependency_files(audit_sources, damage):
    source = audit_sources/'analysis/grand_challenge/reference/adjoint_transport.py'
    if damage == 'missing':
        source.unlink()
    else:
        source.write_text('# Changed independent solver dependency\n', encoding='utf-8')
    with pytest.raises(FileNotFoundError if damage == 'missing' else ValueError):
        audit.reference_report()


def test_reference_rejects_parent_artifact_hash_mismatch(audit_sources):
    # The parent still passes its own source/physical checks after this edit.
    parent = audit.adj.parent_report()
    parent['synthetic_note'] = 'Changed artifact, unchanged physical identity'
    write_json(audit.adj.PARENT, parent)
    assert audit.adj.parent_report() == parent
    with pytest.raises(ValueError, match='unchanged, complete'):
        audit.reference_report()


@pytest.mark.parametrize('field', ['schema', 'all_declared_controls_passed',
    'source_stable_during_run', 'complete_path_evidence_accepted', 'path', 'inputs_SI', 'frequency_hz'])
def test_reference_rejects_incomplete_or_mismatched_identity(audit_sources, field):
    report = json.loads(audit.REFERENCE.read_text(encoding='utf-8'))
    if field == 'schema':
        report[field] = 'wrong-schema'
    elif field == 'path':
        report[field]['residence_time_s'] *= 1.01
    elif field == 'inputs_SI':
        report[field]['pump_power_W'] *= 1.01
    elif field == 'frequency_hz':
        report[field][0] += 1.
    else:
        report[field] = False
    write_json(audit.REFERENCE, report)
    with pytest.raises(ValueError, match='unchanged, complete'):
        audit.reference_report()


def test_missing_reference_report_is_not_accepted(audit_sources):
    audit.REFERENCE.unlink()
    with pytest.raises(FileNotFoundError):
        audit.reference_report()


def test_existing_report_is_never_overwritten(tmp_path):
    path = tmp_path/'existing.json'
    path.write_text('preserve', encoding='utf-8')
    with pytest.raises(FileExistsError):
        audit.write_audit(path)
    assert path.read_text(encoding='utf-8') == 'preserve'
