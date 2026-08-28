"""Fast consistency checks for the analytic-reconstruction audit artifacts."""

import json
import sys
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parent.parent
ANALYSIS = ROOT / "analysis" / "squeezing" / "analytic_reconstruction"
GENERATED = ANALYSIS / "generated"
TEX = (ROOT / "docs" / "FWM physics and analytic reconstruction"
       / "squeezing_analytic_reconstruction_v2.tex")


def _audit():
    return json.loads(
        (GENERATED / "convergence_audit.json").read_text(encoding="utf-8"))


def test_generated_markdown_is_exactly_owned_by_the_audit_source():
    sys.path.insert(0, str(ANALYSIS))
    try:
        import convergence_audit
    finally:
        sys.path.pop(0)
    audit = _audit()
    markdown = (GENERATED / "convergence_audit.md").read_text(encoding="utf-8")
    assert markdown == convergence_audit.make_markdown(audit)


def test_parameter_sensitivity_schema_has_explicit_layer_boundaries():
    option = _audit()["option_a_literature_point"]
    provenance = option["parameter_provenance"]
    sensitivity = option["sensitivity"]
    final = option["option_a_N_F_3"]

    assert set(provenance["parameters"]) == {
        "ell_s", "kappa", "gamma_gg_floor_over_2pi_kHz"}
    assert sensitivity["classification"][
        "illustrative_not_parameter_uncertainty"]
    assert not sensitivity["classification"]["physical_squeezing_prediction"]

    ell = sensitivity["ell_s_propagation_only"]
    assert not ell["atomic_response_recomputed"]
    assert [row["ell_s"] for row in ell["rows"]] == [0.666, 0.74, 0.814]
    assert ell["rows"][1]["G_s"] == pytest.approx(final["G_s"], abs=1e-12)

    gamma = sensitivity["gamma_gg_floor_mean_field"]
    assert gamma["atomic_response_recomputed_for_endpoints"]
    assert [row["gamma_gg_floor_over_2pi_kHz"] for row in gamma["rows"]] == [
        90.0, 100.0, 110.0]
    assert gamma["rows"][1]["G_s"] == pytest.approx(final["G_s"], abs=1e-12)

    kappa = sensitivity["kappa_pump_scatter_diagnostic"]
    assert kappa["equation"] == "N_ps = kappa * (1 - exp(-OD_pump))"
    assert "not physical squeezing" in kappa["classification"]


def test_composite_audit_separates_archived_and_current_solver_scope():
    audit = _audit()
    assert audit["scope"].startswith("Composite analysis artifact")
    archived = audit["archived_reduced_model_classification"]
    assert not archived["corrected_option_a"]
    assert "no trace-preserving thermal transit reload" in archived[
        "atomic_dissipator"]
    assert audit["option_a_literature_point"]["classification"][
        "mean_field_propagation"] == "corrected Option A"


def test_option_a_artifact_matches_the_independent_current_anchor():
    option = _audit()["option_a_literature_point"]
    final = option["option_a_N_F_3"]
    assert final["G_s"] == pytest.approx(403.53549729316836, rel=1e-12)
    assert final["G_c"] == pytest.approx(405.8535914600245, rel=1e-12)
    assert final["photon_flux_gap"] == pytest.approx(
        -2.3115427341867303, rel=1e-12)
    assert max(final["independent_reference_parity"].values()) < 2e-11

    fixture = option["algebraic_dilation_fixture"]
    assert "Caves bound" in fixture["classification"]
    assert "not minimum physical noise" in fixture["classification"]
    assert "D_vacuum_fixture_per_m" in fixture
    assert "D_min_per_m" not in fixture


def test_pump_only_reference_artifact_keeps_the_production_boundary_explicit():
    reference = _audit()["option_a_literature_point"][
        "pump_only_weak_response_reference"]
    assert reference["supported_branch"] == -1
    assert not reference["provenance"]["production_default"]
    assert reference["provenance"]["reference_fields"] == "none"
    assert "unsupported" in reference["plus_branch_status"]
    assert reference["frame_equivalence"][
        "max_abs_static_to_floquet_difference"] < 1e-12
    assert reference["pump_state_diagnostics"][
        "max_response_normalized_residual"] < 1e-12
    assert reference["trace_zero_dc_projection"][
        "max_response_normalized_residual"] < 1e-12
    assert max(row["max_abs_direct_pole_difference_seconds"]
               for row in reference["pole_residue_parity"]) < 1e-12
    seed_limit = reference["finite_seed_to_infinitesimal"]
    assert seed_limit["rows"][-1]["worst_normalized_chi_error"] < 1e-5
    assert seed_limit["smallest_fraction_order_errors"]["1"] > 1e-3
    assert seed_limit["smallest_fraction_order_errors"]["3"] < 1e-5
    assert not reference["remaining_limits"][
        "two_dimensional_noncollinear_doppler"]
    assert not reference["remaining_limits"]["microscopic_langevin_diffusion"]


def test_noncollinear_reference_artifact_converges_without_moving_lab_frequencies():
    option = _audit()["option_a_literature_point"]
    reference = option["noncollinear_doppler_reference"]

    assert option["classification"][
        "noncollinear_doppler_reference_implemented"]
    assert not option["classification"][
        "noncollinear_reference_is_production_default"]
    assert not reference["lab_frequency_contract"][
        "lab_optical_beat_velocity_shifted"]
    assert reference["lab_frequency_contract"][
        "atomic_delta_eff_velocity_shifted"]
    assert reference["rms_budget"]["analytic_total_MHz"] == pytest.approx(
        1.38, rel=0.01)
    assert abs(reference["rms_budget"]["quadrature_error_pct"]) < 0.01
    assert all(reference["acceptance"].values())
    assert abs(reference["grid_refinement"][-2][
        "feature_position_error_to_order40_MHz"]) <= 0.1
    assert not reference["remaining_limits"]["production_default"]
    assert not reference["remaining_limits"][
        "finite_seed_2d_production_path"]
    assert not reference["remaining_limits"][
        "microscopic_langevin_diffusion"]


def test_tex_uses_conditional_math_language_and_current_transit_anchor():
    text = TEX.read_text(encoding="utf-8")
    assert "trace-preserving thermal reload channel" in text
    assert "G_p=403.53550" in text
    assert "Raman pole survives velocity average" not in text
    assert "exceptional point is not physics" not in text.lower()
    assert "universal Caves-bound claim" in text
    assert "illustrative sensitivity sweep" in text
