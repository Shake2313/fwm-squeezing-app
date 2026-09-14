"""Literal phase integration and strict evidence gates for atomic ensembles."""

from dataclasses import replace
import gc
import weakref

import numpy as np
import pytest

from gabes.constants import MASS_85RB
from gabes.fwm_quantum import transport_ensemble as ensemble
from gabes.quantum.contracts import AnalysisFrequencyAxis
from gabes.quantum.inflow import maxwell_box_inflow


LABELS = ("probe", "conjugate", "probe_dagger", "conjugate_dagger")
SOURCES = ("atomic_inflow", "jump:analytic_reservoir")
CHARGES = np.array([1, -1, -1, 1])
AXIS = AnalysisFrequencyAxis.from_hz([.15e6, 1.25e6])
PATH_METRICS = ("greater", "lesser", "greater_by_source", "lesser_by_source",
                "mean_pulse", "mean_outer", "retarded_response")
ENSEMBLE_METRICS = ("greater", "lesser", "greater_by_source", "lesser_by_source",
                    "poisson_number", "retarded_response")


def inflow(*, density=3e16, power=0):
    return maxwell_box_inflow([-.0002, -.0001, -.0001], [.0002, .0001, .0001],
        temperature_K=373., mass_kg=MASS_85RB, density_m3=density,
        points_per_face_power=power, seed=11, source="explicit analytic transport fixture")


def convention(q, **changes):
    arguments = dict(analysis_axis=AXIS, number_density_m3=q.density_m3,
        coupling_scales=[2., 3., 2., 3.], readout_scales=[1., 1., 1., 1.],
        source_names=SOURCES, carrier_offsets_at_rest_rad_s=[-8e5, 8e5, 8e5, -8e5],
        port_wavevectors_rad_m=[[10., 20., 30.], [4., -6., 9.], [-10., -20., -30.], [-4., 6., -9.]],
        model_id="analytic-shared-phase-fixture-v1", model_scope="analytic_fixture",
        provenance="closed-form synthetic atomic moments, not a Rb prediction",
        phase_provenance="one laboratory phase with exact signed harmonic charges")
    arguments.update(changes)
    return ensemble.PumpOnlyConvention(**arguments)


def raw_packet(path, contract):
    """PSD source blocks with known phase covariance, varying across chords."""
    rng = np.random.default_rng(617)
    a = rng.normal(size=(2, 2, 4, 4))+1j*rng.normal(size=(2, 2, 4, 4))
    b = rng.normal(size=(2, 2, 4, 4))+1j*rng.normal(size=(2, 2, 4, 4))
    scale = path.residence_time_s**2
    greater = scale*(a@a.conj().swapaxes(-1, -2))
    lesser = .3*scale*(b@b.conj().swapaxes(-1, -2))
    mean = path.residence_time_s*np.array([[1+.2j, .6-.3j, .4+.7j, 1.2-.2j],
                                         [.3-.1j, 1+.4j, .7-.6j, .2+.9j]])
    response = scale*(rng.normal(size=(2, 4, 4))+1j*rng.normal(size=(2, 4, 4)))
    offsets = np.asarray(contract.carrier_offsets_at_rest_rad_s)-np.asarray(contract.port_wavevectors_rad_m)@path.velocity_m_s
    return {"analysis_axis": contract.analysis_axis,
        "frequencies_rad_s": contract.analysis_axis.omega_rad_s[:, None]+offsets,
        "source_names": SOURCES, "greater_by_source": greater, "lesser_by_source": lesser,
        "greater": greater.sum(axis=0), "lesser": lesser.sum(axis=0),
        "mean_pulse": mean, "retarded_response": response,
        "residence_time_s": path.residence_time_s, "audit": {"passed": True},
        "metadata": {"mode_labels": LABELS,
            "coupling_s_inverse_sqrt_flux": np.asarray(contract.coupling_scales)[:2].tolist(),
            "carrier_offsets_rad_s": offsets.tolist(), "entry_phase_rad": 0.,
            "entry_optical_demodulation_phases_rad": (-np.asarray(contract.port_wavevectors_rad_m)[:2]@path.entry_position_m).tolist(),
            "scope": "analytic pump-only common-phase fixture; not experimental evidence"}}


def literal_phase_average(packet):
    """Independent integration of D(phi), without constructing a charge mask."""
    result = {key: np.zeros_like(packet[key]) for key in
              ("greater", "lesser", "greater_by_source", "lesser_by_source", "retarded_response")}
    result["poisson_number"] = np.zeros_like(packet["greater"])
    result["mean_pulse"] = np.zeros_like(packet["mean_pulse"])
    for phase in .173+2*np.pi*np.arange(8)/8:
        diagonal = np.diag(np.exp(1j*CHARGES*phase))
        for key in ("greater", "lesser", "greater_by_source", "lesser_by_source", "retarded_response"):
            result[key] += (diagonal@packet[key]@diagonal.conj().T)/8
        mean = packet["mean_pulse"]@diagonal.T
        result["poisson_number"] += mean[:, :, None]*mean[:, None, :].conj()/8
        result["mean_pulse"] += mean/8
    return result


def comparison(kind, target, *, reference="independent reference", controls=(1., 2.), errors=None):
    metrics = ENSEMBLE_METRICS if kind in ("ensemble_refinement", "independent_scramble") else PATH_METRICS
    observed = {name: 1e-12 for name in metrics} if errors is None else errors
    return ensemble.NumericalComparison(kind=kind, errors=observed,
        tolerances={name: 1e-8 for name in observed}, reference_id=reference,
        candidate_id=target, control_values=controls,
        error_definition="relative Frobenius errors with the declared analytic tau-squared floor for zero covariance",
        provenance="synthetic comparison ledger for testing the certification contract")


def evidence(target, *, scope="analytic_fixture", stream=False, two_refinements=False):
    if stream:
        rows = (comparison("ensemble_refinement", target), comparison("independent_scramble", target, controls=(11, 211)))
    else:
        rows = (comparison("path_refinement", target), comparison("independent_reference", target))
        if two_refinements:
            rows = (comparison("path_refinement", "middle-resolution", reference="coarse-resolution"),
                    comparison("path_refinement", target, reference="middle-resolution", controls=(2., 4.)), rows[1])
    return ensemble.ConvergenceEvidence(target, scope, rows, "declared analytic fixture evidence, no apparatus validation")


def supplied(path, contract, *, proof=True, two_refinements=False):
    packet = raw_packet(path, contract)
    token = ensemble.packet_digest(path, packet, contract)
    proof = evidence(token, scope=contract.model_scope, two_refinements=two_refinements) if proof else None
    return ensemble.SuppliedPathPacket(packet, contract, proof)


def candidate(*, q=None, contract=None):
    q = inflow() if q is None else q
    contract = convention(q) if contract is None else contract
    return ensemble.aggregate_pump_only_stream(q, contract,
        lambda path: supplied(path, contract, two_refinements=contract.model_scope == "smooth"))


def test_stream_matches_literal_common_phase_integration_of_every_ordered_source_and_response():
    q = inflow()
    contract = convention(q)
    expected = None
    offsets = []
    for index, rate in enumerate(q.rate_s_inverse):
        packet = raw_packet(q.path(index), contract)
        offsets.append(packet["metadata"]["carrier_offsets_rad_s"])
        phase_average = literal_phase_average(packet)
        if expected is None:
            expected = {name: np.zeros_like(value) for name, value in phase_average.items()}
        for name, value in phase_average.items():
            expected[name] += rate*value
    actual = candidate(q=q, contract=contract)
    spectra = actual["spectra"]
    for name in ("greater_by_source", "lesser_by_source", "retarded_response", "poisson_number"):
        scale = np.linalg.norm(expected[name])
        np.testing.assert_allclose(spectra[name], expected[name], rtol=2e-14, atol=2e-15*scale)
    for ordering in ("greater", "lesser"):
        reference = expected[ordering]+expected["poisson_number"]
        np.testing.assert_allclose(spectra[ordering], reference, rtol=2e-14, atol=2e-15*np.linalg.norm(reference))
        np.testing.assert_allclose(spectra["internal_"+ordering], spectra[ordering+"_by_source"].sum(axis=0), rtol=2e-15)
    assert len(np.unique(np.array(offsets), axis=0)) == len(q.rate_s_inverse)
    assert spectra["poisson_number"][0, 0, 3] != 0
    assert spectra["greater"][0, 1, 2] != 0
    for name in spectra:
        np.testing.assert_array_equal(spectra[name][..., 0, 1], 0.)
    assert actual["path_evidence_passed"] and not actual["certified"]
    assert actual["analysis_axis"] is contract.analysis_axis


def test_mean_term_cannot_be_erased_when_phase_averaged_mean_vanishes():
    q = inflow()
    contract = convention(q)
    def deterministic(path):
        packet = raw_packet(path, contract)
        for name in ("greater", "lesser", "greater_by_source", "lesser_by_source"):
            packet[name] = np.zeros_like(packet[name])
        token = ensemble.packet_digest(path, packet, contract)
        return ensemble.SuppliedPathPacket(packet, contract, evidence(token))
    result = ensemble.aggregate_pump_only_stream(q, contract, deterministic)
    spectra = result["spectra"]
    np.testing.assert_array_equal(spectra["internal_greater"], np.zeros((2, 4, 4)))
    np.testing.assert_array_equal(spectra["greater"], spectra["poisson_number"])
    assert np.linalg.norm(spectra["greater"]) > 0
    average = literal_phase_average(raw_packet(q.path(0), contract))
    assert np.linalg.norm(average["mean_pulse"]) < 1e-14*q.residence_time_s[0]
    assert np.linalg.norm(average["poisson_number"]) > 0


def test_each_path_is_consumed_once_without_retaining_the_entire_packet_history():
    q = inflow(power=2)
    contract = convention(q)
    references, paths = [], []
    def factory(path):
        gc.collect()
        assert sum(item() is not None for item in references) <= 1
        item = supplied(path, contract)
        references.append(weakref.ref(item.packet["greater"]))
        paths.append(path.source)
        return item
    result = ensemble.aggregate_pump_only_stream(q, contract, factory)
    assert len(paths) == len(set(paths)) == len(q.rate_s_inverse)
    assert result["consumed_path_count"] == len(q.rate_s_inverse)
    certified = ensemble.certify_pump_only_stream(result, evidence(result["candidate_digest"], stream=True))
    assert certified["certified"] and len(paths) == len(q.rate_s_inverse)
    assert not result["certified"]
    with pytest.raises(TypeError):
        certified["spectra"]["greater"] = np.zeros((2, 4, 4))
    for array in certified["spectra"].values():
        with pytest.raises(ValueError):
            array.setflags(write=True)


def test_density_is_applied_once_and_occupancy_is_not_renormalized_to_nV():
    a, b = inflow(), inflow(density=9e16)
    left, right = candidate(q=a), candidate(q=b)
    for key in left["spectra"]:
        np.testing.assert_allclose(right["spectra"][key], 3*left["spectra"][key], rtol=3e-15)
    assert left["mean_occupancy"] == a.mean_occupancy
    assert left["equilibrium_nV"] == a.equilibrium_atom_number
    assert left["mean_occupancy"] != left["equilibrium_nV"]
    assert right["mean_occupancy"] == pytest.approx(3*left["mean_occupancy"], rel=3e-15)
    calls = []
    with pytest.raises(ValueError, match="density"):
        ensemble.aggregate_pump_only_stream(a, convention(b), lambda path: calls.append(path))
    assert calls == []


@pytest.mark.parametrize("change", [
    {"mode_labels": tuple(reversed(LABELS))}, {"phase_charges": (1, -1, 1, -1)},
    {"readout_units": ("C m",)*4}, {"drive_units": ("W",)*4},
    {"coupling_scales": [2., 3., 2., 4.]}, {"readout_scales": [1., 0., 1., 0.]},
    {"source_names": tuple(reversed(SOURCES))}, {"source_names": ("atomic_inflow", "atomic_inflow")},
    {"carrier_offsets_at_rest_rad_s": [-8e5, 8e5, 8e5, 8e5]},
    {"carrier_offsets_at_rest_rad_s": [-8e5, 6e5, 8e5, -6e5]},
    {"port_wavevectors_rad_m": np.zeros((4, 2))}, {"provenance": ""},
    {"phase_provenance": ""}, {"demodulation": "atomic rotating frame"},
    {"number_density_m3": np.nan}, {"model_scope": "hot_vapor_experiment"},
])
def test_physical_convention_cannot_hide_order_units_scales_or_provenance_mismatch(change):
    with pytest.raises((ValueError, TypeError)):
        convention(inflow(), **change)


@pytest.mark.parametrize("failure", ["axis", "frequency", "offset", "sources", "nambu", "coupling",
    "units", "density", "readout_scale", "source_closure", "negative_source", "duration", "shape", "scope",
    "spatial_phase", "pump_only", "phase_charges", "demodulation"])
def test_native_packet_contract_is_checked_even_when_old_evidence_claims_success(failure):
    q = inflow()
    contract = convention(q)
    def factory(path):
        item = supplied(path, contract)
        p, m = item.packet, item.packet["metadata"]
        if failure == "axis":
            p["analysis_axis"] = AnalysisFrequencyAxis.from_hz([.16e6, 1.25e6])
        elif failure == "frequency":
            p["frequencies_rad_s"] = p["frequencies_rad_s"]+1.
        elif failure == "offset":
            m["carrier_offsets_rad_s"][0] += 1.
        elif failure == "sources":
            p["source_names"] = tuple(reversed(SOURCES))
        elif failure == "nambu":
            m["mode_labels"] = tuple(reversed(LABELS))
        elif failure == "coupling":
            m["coupling_s_inverse_sqrt_flux"][0] *= 2
        elif failure == "units":
            m["readout_units"] = ("C m",)*4
        elif failure == "density":
            m["number_density_m3"] = 2*q.density_m3
        elif failure == "readout_scale":
            m["readout_scales"] = [2., 2., 2., 2.]
        elif failure == "source_closure":
            p["greater"] = p["greater_by_source"][1]
        elif failure == "negative_source":
            p["greater_by_source"] = -p["greater_by_source"]
            p["greater"] = p["greater_by_source"].sum(axis=0)
        elif failure == "duration":
            p["residence_time_s"] *= 1.1
        elif failure == "shape":
            p["mean_pulse"] = p["mean_pulse"][:, :2]
        elif failure == "scope":
            m["scope"] = ""
        elif failure == "spatial_phase":
            m["entry_optical_demodulation_phases_rad"][0] += .3
        elif failure == "pump_only":
            m["pump_only"] = False
        elif failure == "phase_charges":
            m["phase_charges"] = [1, -1, 1, -1]
        elif failure == "demodulation":
            m["demodulation"] = "rotating frame"
        return item
    with pytest.raises((ValueError, TypeError)):
        ensemble.aggregate_pump_only_stream(q, contract, factory)


def test_strong_rf_response_cannot_hide_a_missing_source_at_a_weak_rf_sample():
    q = inflow()
    contract = convention(q)
    def factory(path):
        packet = raw_packet(path, contract)
        for ordering in ("greater", "lesser"):
            packet[ordering+"_by_source"][:, 1] *= 1e-12
            packet[ordering] = packet[ordering+"_by_source"].sum(axis=0)
        proof = evidence(ensemble.packet_digest(path, packet, contract))
        # The complete weak row is representable and valid. Dropping its first
        # source must not hide under a tolerance set by the strong first RF row.
        packet["greater"][1] = packet["greater_by_source"][1, 1]
        return ensemble.SuppliedPathPacket(packet, contract, proof)
    with pytest.raises(ValueError, match="source sum"):
        ensemble.aggregate_pump_only_stream(q, contract, factory)


@pytest.mark.parametrize("field,value", [("provenance", "changed source"), ("model_id", "other model"),
    ("readout_scales", [2., 2., 2., 2.]), ("number_density_m3", 6e16)])
def test_supplied_packet_signature_must_match_the_same_consumed_stream(field, value):
    q = inflow()
    contract = convention(q)
    changed = replace(contract, **{field: value})
    with pytest.raises(ValueError, match="convention differs"):
        ensemble.aggregate_pump_only_stream(q, contract, lambda path: supplied(path, changed))


@pytest.mark.parametrize("failure", ["missing", "boolean", "physical_audit", "failed_budget", "missing_mean_outer", "wrong_digest", "missing_reference", "tampered_mean"])
def test_unproved_or_failed_path_aborts_with_no_partial_or_renormalized_spectrum(failure):
    q = inflow()
    contract = convention(q)
    calls = []
    def factory(path):
        calls.append(path)
        item = supplied(path, contract)
        if len(calls) == 2:
            proof = item.evidence
            if failure == "missing":
                item = replace(item, evidence=None)
            elif failure == "boolean":
                item = replace(item, evidence={"passed": True})
            elif failure == "physical_audit":
                item.packet["audit"] = {"passed": False}
            elif failure == "failed_budget":
                row = proof.comparisons[0]
                item = replace(item, evidence=replace(proof, comparisons=(replace(row, errors={**row.errors, "greater": 1.}), proof.comparisons[1])))
            elif failure == "missing_mean_outer":
                row = proof.comparisons[0]
                item = replace(item, evidence=replace(proof, comparisons=(replace(row,
                    errors={k: v for k, v in row.errors.items() if k != "mean_outer"},
                    tolerances={k: v for k, v in row.tolerances.items() if k != "mean_outer"}), proof.comparisons[1])))
            elif failure == "wrong_digest":
                item = replace(item, evidence=replace(proof, target_digest="0"*64))
            elif failure == "missing_reference":
                item = replace(item, evidence=replace(proof, comparisons=(proof.comparisons[0],)))
            elif failure == "tampered_mean":
                item.packet["mean_pulse"] *= 2
        return item
    result = ensemble.aggregate_pump_only_stream(q, contract, factory)
    assert len(calls) == result["consumed_path_count"] == 2
    assert result["spectra"] is None and not result["certified"] and not result["path_evidence_passed"]
    assert result["failed_path_index"] == 1
    assert result["failed_path_rate_s_inverse"] == q.rate_s_inverse[1]
    assert result["reasons"] and result["mean_occupancy"] == q.mean_occupancy


def test_smooth_paths_require_three_linked_resolutions_and_never_inherit_segment_only_proof():
    q = inflow()
    contract = convention(q, model_scope="smooth")
    inadequate = ensemble.aggregate_pump_only_stream(q, contract, lambda path: supplied(path, contract))
    assert inadequate["spectra"] is None and any("2 explicit" in reason for reason in inadequate["reasons"])
    good = candidate(q=q, contract=contract)
    assert good["path_evidence_passed"]
    def disconnected(path):
        item = supplied(path, contract, two_refinements=True)
        rows = item.evidence.comparisons
        changed = replace(rows[1], reference_id="unrelated resolution")
        return replace(item, evidence=replace(item.evidence, comparisons=(rows[0], changed, rows[2])))
    bad = ensemble.aggregate_pump_only_stream(q, contract, disconnected)
    assert bad["spectra"] is None and any("link two successive" in reason for reason in bad["reasons"])


@pytest.mark.parametrize("failure", ["missing", "boolean", "wrong_target", "moment_only", "missing_scramble", "missing_number", "failed_budget"])
def test_path_convergence_or_nV_moments_cannot_certify_the_actual_ensemble_integrand(failure):
    current = candidate()
    proof = evidence(current["candidate_digest"], stream=True)
    if failure == "missing":
        proof = None
    elif failure == "boolean":
        proof = {"passed": True}
    elif failure == "wrong_target":
        proof = replace(proof, target_digest="0"*64)
    elif failure == "moment_only":
        proof = replace(proof, comparisons=tuple(replace(row, errors={"nV": 0.}, tolerances={"nV": .01}) for row in proof.comparisons))
    elif failure == "missing_scramble":
        proof = replace(proof, comparisons=(proof.comparisons[0],))
    elif failure == "missing_number":
        rows = tuple(replace(row, errors={k: v for k, v in row.errors.items() if k != "poisson_number"},
            tolerances={k: v for k, v in row.tolerances.items() if k != "poisson_number"}) for row in proof.comparisons)
        proof = replace(proof, comparisons=rows)
    elif failure == "failed_budget":
        row = proof.comparisons[1]
        proof = replace(proof, comparisons=(proof.comparisons[0], replace(row, errors={**row.errors, "retarded_response": 1.})))
    result = ensemble.certify_pump_only_stream(current, proof)
    assert not result["certified"] and result["reasons"]


@pytest.mark.parametrize("field,value", [("analysis_axis", AnalysisFrequencyAxis.from_hz([.16e6, 1.25e6])),
    ("model_scope", "smooth"), ("model_id", "other physical model"),
    ("source_names", tuple(reversed(SOURCES))), ("numerical_units", {"covariance": "A^2/Hz", "retarded_response": "1"})])
def test_candidate_frequency_scope_sources_and_units_cannot_change_after_evidence_binding(field, value):
    current = candidate()
    current[field] = value
    proof = evidence(current["candidate_digest"], scope=current["model_scope"], stream=True)
    with pytest.raises(ValueError, match="changed"):
        ensemble.certify_pump_only_stream(current, proof)


def test_candidate_spectral_data_cannot_be_swapped_after_the_computation_digest():
    current = candidate()
    current["spectra"] = dict(current["spectra"])
    current["spectra"]["poisson_number"] = 2*current["spectra"]["poisson_number"]
    with pytest.raises(ValueError, match="changed"):
        ensemble.certify_pump_only_stream(current, evidence(current["candidate_digest"], stream=True))


def test_numerical_comparisons_own_their_budgets_and_cannot_claim_self_reference():
    errors = {name: 0. for name in PATH_METRICS}
    row = comparison("independent_reference", "candidate", errors=errors)
    errors["greater"] = 1.
    assert row.passed
    with pytest.raises(TypeError):
        row.errors["greater"] = 1.
    with pytest.raises(ValueError):
        comparison("independent_reference", "same-run", reference="same-run")
    for change in ({"errors": {**row.errors, "greater": np.nan}},
                   {"errors": {**row.errors, "greater": -1.}},
                   {"tolerances": {**row.tolerances, "greater": 0.}},
                   {"provenance": ""}, {"control_values": (2., 2.)}):
        with pytest.raises(ValueError):
            replace(row, **change)
