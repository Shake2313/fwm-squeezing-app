"""Synthetic gate data and temporary provenance; never historical run evidence."""

from dataclasses import asdict
import json

import numpy as np
import pytest

from analysis.grand_challenge import adjoint_transport_audit as adj
from analysis.grand_challenge import exponential_transport_audit as exponential
from analysis.grand_challenge import rb_thermal_reference_jobs as refs


def write_json(path, value):
    path.write_text(json.dumps(adj.encode(value), allow_nan=False), encoding='utf-8')


def synthetic_packet(problem, convention, path, *, correction=0., rtol=1e-9):
    """Native shapes/metadata with explicit synthetic, nonzero ordered outputs."""
    nf, ports = problem['frequencies_rad_s'].shape
    ns = len(convention.source_names)
    tau = path.residence_time_s
    scale = tau**2*(1+correction)
    sources = np.broadcast_to(np.eye(ports, dtype=complex), (ns, nf, ports, ports)).copy()
    sources *= scale*np.arange(1, ns+1)[:, None, None, None]/(100*ns)
    mean = tau*(1+correction)*np.broadcast_to(
        [.04+.01j, .02-.03j, .04-.01j, .02+.03j], (nf, ports)).copy()
    return {
        'analysis_axis': convention.analysis_axis, 'metadata': problem['metadata'],
        'source_names': convention.source_names,
        'frequencies_rad_s': problem['frequencies_rad_s'], 'residence_time_s': tau,
        'greater_by_source': sources, 'lesser_by_source': .7*sources,
        'greater': sources.sum(axis=0), 'lesser': (.7*sources).sum(axis=0),
        'mean_pulse': mean, 'mean_outer': mean[:, :, None]*mean[:, None, :].conj(),
        'retarded_response': np.broadcast_to(
            (.2+.3j)*scale*np.eye(ports), (nf, ports, ports)).copy(),
        'exit_state': np.asarray(problem['boundary_state'], complex).copy(),
        'numerics': {'rtol': rtol, 'synthetic_test_fixture': True},
    }


@pytest.fixture
def audit_sources(tmp_path, monkeypatch):
    """Run the real historical validators against a wholly temporary source tree.

    These files stand for test dependencies, not copies of certified physics.
    Only path constants change; hashing, fixture matching and source checks run.
    """
    root = tmp_path/'synthetic-audit-sources'
    names = (
        'physics.py',
        'analysis/grand_challenge/reference/adjoint_transport.py',
        'analysis/grand_challenge/adjoint_transport_audit.py',
        'gabes/fwm_quantum/transport_ensemble.py', 'gabes/quantum/inflow.py',
        'tests/quantum/test_adjoint_transport.py', 'tests/quantum/test_adjoint_transport_audit.py',
        refs.DRIVER,
        'analysis/grand_challenge/reference/exponential_transport.py',
        'analysis/grand_challenge/exponential_transport_audit.py',
        'tests/quantum/test_exponential_transport.py', 'tests/quantum/test_exponential_transport_audit.py',
    )
    for name in names:
        source = root/name
        source.parent.mkdir(parents=True, exist_ok=True)
        source.write_text('# Synthetic provenance dependency; no physical calculation.\n', encoding='utf-8')
    monkeypatch.setattr(adj, 'ROOT', root)
    monkeypatch.setattr(adj, 'PARENT', root/'parent.json')
    monkeypatch.setattr(exponential, 'ROOT', root)
    monkeypatch.setattr(exponential, 'REFERENCE', root/'reference.json')
    monkeypatch.setattr(refs, 'ROOT', root)
    inputs, geometry, path, axis = adj.fixture()
    parent = {
        'schema': 'gabes-continuous-rb-transport-v2', 'synthetic_test_fixture': True,
        'all_declared_controls_passed': True, 'source_stable_during_run': True,
        'source_sha256': {'physics.py': adj.file_hash(root/'physics.py')},
        'inputs_SI': asdict(inputs),
        'path': {'entry_position_m': path.entry_position_m.tolist(),
                 'velocity_m_s': path.velocity_m_s.tolist(), 'residence_time_s': path.residence_time_s},
        'analysis_frequencies_hz': axis.frequency_hz.tolist(),
        'wavevectors_rad_m': geometry.wavevectors_rad_m.tolist(),
    }
    write_json(adj.PARENT, parent)
    reference = {
        'schema': 'gabes-adjoint-source-response-audit-v1', 'synthetic_test_fixture': True,
        'all_declared_controls_passed': True, 'source_stable_during_run': True,
        'complete_path_evidence_accepted': True, 'source_sha256': adj.hashes(),
        'parent_report': {'sha256': adj.file_hash(adj.PARENT)},
        'path': parent['path'], 'inputs_SI': parent['inputs_SI'],
        'frequency_hz': parent['analysis_frequencies_hz'],
    }
    write_json(exponential.REFERENCE, reference)
    return root
