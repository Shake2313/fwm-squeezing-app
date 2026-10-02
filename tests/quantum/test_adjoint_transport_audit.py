"""Source-wise gates and real provenance checks on synthetic temporary evidence."""

import json

import numpy as np
import pytest

from analysis.grand_challenge import adjoint_transport_audit as audit
from _grand_challenge_audit_fixtures import audit_sources, synthetic_packet, write_json


def test_parent_reuse_binds_physical_packet_and_all_metrics():
    p, convention, path = audit.problem()
    parent = {'cases': {f'primary_{j}': {'values': audit.encode({
        key: value for key, value in synthetic_packet(
            p, convention, path, correction=1e-6/(j+1), rtol=rtol).items()
        if key != 'analysis_axis'})}
        for j, rtol in enumerate((3e-9, 1e-9, 3e-10))}}
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


def test_parent_validator_accepts_temporary_native_fixture(audit_sources):
    report = audit.parent_report()
    assert report['synthetic_test_fixture']
    assert report['source_sha256']['physics.py'] == audit.file_hash(audit_sources/'physics.py')


@pytest.mark.parametrize('damage', ['changed', 'missing'])
def test_reuse_refuses_changed_or_missing_parent_source(audit_sources, damage):
    source = audit_sources/'physics.py'
    if damage == 'missing':
        source.unlink()
    else:
        source.write_text('# Changed synthetic physics dependency\n', encoding='utf-8')
    with pytest.raises(FileNotFoundError if damage == 'missing' else ValueError):
        audit.parent_report()


@pytest.mark.parametrize('field', ['inputs', 'entry', 'velocity', 'duration', 'rf', 'wavevectors'])
def test_parent_validator_rejects_physical_identity_mismatch(audit_sources, field):
    report = json.loads(audit.PARENT.read_text(encoding='utf-8'))
    if field == 'inputs':
        report['inputs_SI']['pump_power_W'] *= 1.01
    elif field in ('entry', 'velocity'):
        name = {'entry': 'entry_position_m', 'velocity': 'velocity_m_s'}[field]
        report['path'][name][0] += 1e-5
    elif field == 'duration':
        report['path']['residence_time_s'] *= 1.01
    elif field == 'rf':
        report['analysis_frequencies_hz'][0] += 1.
    else:
        report['wavevectors_rad_m'][0][0] += 1.
    write_json(audit.PARENT, report)
    with pytest.raises(ValueError, match='parent physical fixture differs'):
        audit.parent_report()


@pytest.mark.parametrize('field', ['schema', 'all_declared_controls_passed', 'source_stable_during_run'])
def test_parent_validator_rejects_unverified_report(audit_sources, field):
    report = json.loads(audit.PARENT.read_text(encoding='utf-8'))
    report[field] = 'wrong-schema' if field == 'schema' else False
    write_json(audit.PARENT, report)
    with pytest.raises(ValueError, match='verified continuous selected-path parent required'):
        audit.parent_report()


def test_missing_parent_report_is_not_accepted(audit_sources):
    audit.PARENT.unlink()
    with pytest.raises(FileNotFoundError):
        audit.parent_report()


def test_existing_output_refused_before_parent_or_workers(tmp_path, monkeypatch):
    path = tmp_path/'preserved.json'
    path.write_text('preserve')
    monkeypatch.setattr(audit, 'hashes', lambda: pytest.fail('must refuse before reading parent'))
    with pytest.raises(FileExistsError):
        audit.write_audit(path, tmp_path/'unused.png')
    assert path.read_text() == 'preserve'
