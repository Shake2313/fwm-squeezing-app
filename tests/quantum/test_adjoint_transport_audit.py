"""Provenance, source-wise norms and complete selected-path evidence gates."""

import numpy as np
import pytest

from analysis.grand_challenge import adjoint_transport_audit as audit


def test_parent_reuse_binds_physical_packet_and_all_metrics():
    parent = audit.parent_report()
    p, convention, path = audit.problem()
    packets = audit.primary_packets(parent, p, convention, path)
    evidence, ids = audit.evidence_for(packets, packets[-1], p, convention, path)
    assert len(set(ids)) == 3
    assert set(evidence.comparisons[-1].errors) == set(audit.METRICS)
    assert max(evidence.comparisons[-1].errors.values()) == 0.
    assert evidence.comparisons[0].candidate_id == evidence.comparisons[1].reference_id
    assert not audit.stream._evidence_reasons(evidence, ids[-1], 'smooth')


def test_source_wise_norm_does_not_hide_small_source_behind_large_one():
    p = {'duration_s': 1e-6, 'readouts': np.ones((2, 2, 2)), 'drives': np.ones((2, 2, 2))}
    packet = {key: np.ones((2, 2, 2), complex)*1e-12 for key in audit.METRICS}
    packet['mean_pulse'] = np.ones((2, 2), complex)*1e-6
    packet['exit_state'] = np.eye(2)/2
    for key in ('greater_by_source', 'lesser_by_source'):
        packet[key] = np.ones((2, 2, 2, 2), complex)*1e-12
        packet[key][1] *= 1e-10
    modified = {key: value.copy() for key, value in packet.items()}
    modified['greater_by_source'][1, 1] *= 1.01
    result = audit.errors(modified, packet, p)
    assert result['greater_by_source'] == pytest.approx(.01)
    assert result['greater'] == 0.


def test_reuse_refuses_changed_parent_source(monkeypatch):
    monkeypatch.setattr(audit, 'file_hash', lambda path: 'changed')
    with pytest.raises(ValueError, match='parent source changed'):
        audit.parent_report()


def test_existing_output_refused_before_parent_or_workers(tmp_path, monkeypatch):
    path = tmp_path/'preserved.json'
    path.write_text('preserve')
    monkeypatch.setattr(audit, 'hashes', lambda: pytest.fail('must refuse before reading parent'))
    with pytest.raises(FileExistsError):
        audit.write_audit(path, tmp_path/'unused.png')
    assert path.read_text() == 'preserve'
