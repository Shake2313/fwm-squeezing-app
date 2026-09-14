"""Fast/Balanced Squeezing tiers: pole–residue responses, adaptive scan and provenance."""
import numpy as np
import pytest

from gabes import kernels
from gabes.schemes import fwm

SMALL = dict(T=394.15, P_pump=0.6, P_probe=8e-6, line_strength=0.74,
             coarse_points=41, fine_points=0, velocity_step=20.0,
             velocity_cutoff=4.0, branch=-1)


def _small(method, phase=fwm.PHASE_ULTRA, D=0.9, **extra):
    center = fwm.branch_center_GHz(D, -1)
    return fwm.compute_spectrum(
        D, scan_min=center - 0.3, scan_max=center + 0.3, phase_detail=phase,
        response_method=method, **SMALL, **extra)


@pytest.mark.parametrize("phase", (fwm.PHASE_BALANCED, fwm.PHASE_FINE, fwm.PHASE_ULTRA))
def test_pole_response_reproduces_the_grid_response(phase):
    grid = _small(fwm.RESPONSE_GRID, phase)
    pole = _small(fwm.RESPONSE_POLE, phase)
    for key in ("G_s", "G_c", "G_s_smallsignal", "G_c_smallsignal"):
        np.testing.assert_allclose(pole[key], grid[key], rtol=1e-7)
    # Ultra's pump-scatter term uses the analytic Voigt in the pole tiers.
    atol = 5e-3 if phase == fwm.PHASE_ULTRA else 1e-9
    np.testing.assert_allclose(pole["S_dB"], grid["S_dB"], rtol=0, atol=atol)
    assert pole["floquet_convergence"]["status"] == grid["floquet_convergence"]["status"]
    estimator = pole["response_estimator"]
    assert estimator["method"] == fwm.RESPONSE_POLE
    assert estimator["solved_probe_points"] == estimator["display_probe_points"] == 41
    assert estimator["interpolation"] == "none"
    assert estimator["pole_guard_max_relative"] <= fwm.POLE_GUARD_RTOL
    assert estimator["pole_guard_fallback_rows"] == 0
    assert grid["response_estimator"]["method"] == fwm.RESPONSE_GRID


def test_numba_disabled_pole_path_matches_the_compiled_path(monkeypatch):
    if not kernels.available():
        pytest.skip("compiled path unavailable")
    compiled = _small(fwm.RESPONSE_POLE, fwm.PHASE_BALANCED)
    monkeypatch.setattr(kernels, "available", lambda: False)
    fallback = _small(fwm.RESPONSE_POLE, fwm.PHASE_BALANCED)
    for key in ("G_s", "G_c", "S_dB"):
        np.testing.assert_allclose(fallback[key], compiled[key], rtol=1e-9, atol=1e-12)


def test_pole_guard_replaces_every_failing_row_with_the_grid_solution(monkeypatch):
    grid = _small(fwm.RESPONSE_GRID, fwm.PHASE_BALANCED)
    monkeypatch.setattr(fwm, "POLE_GUARD_RTOL", -1.0)
    forced = _small(fwm.RESPONSE_POLE, fwm.PHASE_BALANCED)
    assert forced["response_estimator"]["pole_guard_fallback_rows"] == 41
    for key in ("G_s", "G_c", "S_dB"):
        np.testing.assert_allclose(forced[key], grid[key], rtol=1e-12, atol=1e-12)


def test_pole_request_falls_back_to_the_grid_when_identities_fail(monkeypatch):
    monkeypatch.setattr(fwm, "_seeded_pole_systems", lambda *args, **kwargs: None)
    spec = _small(fwm.RESPONSE_POLE_ADAPTIVE, fwm.PHASE_BALANCED)
    assert spec["response_estimator"]["method"] == fwm.RESPONSE_GRID
    assert "fallback_reason" in spec["response_estimator"]


def test_unknown_response_method_is_rejected():
    with pytest.raises(ValueError, match="response_method"):
        _small("velocity_grid_typo", fwm.PHASE_BALANCED)


def test_adaptive_scan_needs_fewer_rows_and_stays_on_the_exact_curve():
    center = fwm.branch_center_GHz(0.9, -1)
    common = dict(T=394.15, P_pump=0.6, P_probe=8e-6, line_strength=0.74,
                  coarse_points=401, fine_points=0, velocity_step=1.0,
                  velocity_cutoff=4.0, branch=-1, phase_detail=fwm.PHASE_ULTRA,
                  scan_min=center - fwm.WINDOW_GHZ, scan_max=center + fwm.WINDOW_GHZ)
    exact = fwm.compute_spectrum(0.9, response_method=fwm.RESPONSE_POLE, **common)
    fast = fwm.compute_spectrum(0.9, response_method=fwm.RESPONSE_POLE_ADAPTIVE, **common)
    estimator = fast["response_estimator"]
    solved = np.asarray(estimator["solved_probe_indices"])
    assert estimator["solved_probe_points"] == solved.size < 200
    assert solved[0] == 0 and solved[-1] == 400
    assert estimator["interpolation"] != "none"
    np.testing.assert_array_equal(fast["probe_axis_GHz"], exact["probe_axis_GHz"])
    # Solved rows are exact; interpolated rows stay near the adaptive tolerance.
    np.testing.assert_allclose(fast["S_dB"][solved], exact["S_dB"][solved], rtol=0, atol=1e-9)
    assert np.max(np.abs(fast["S_dB"] - exact["S_dB"])) < 0.2
    assert np.percentile(np.abs(fast["S_dB"] - exact["S_dB"]), 95) < 0.03
    assert fast["floquet_convergence"]["status"] == "CONVERGED"
    assert fast["floquet_convergence"]["scope"].startswith("solved probe detunings")
    assert any("interpolated" in reason for reason in fast["claim_gate"]["reasons"])
    assert not any("interpolated" in reason for reason in exact["claim_gate"]["reasons"])


def test_app_tiers_use_the_pole_response_on_the_ultra_display_axis():
    scheme = fwm.FWMScheme()
    params = scheme.defaults()
    fast = scheme.compute(dict(params, resolution=fwm.FIDELITY_FAST))
    balanced = scheme.compute(dict(params, resolution=fwm.FIDELITY_BALANCED))
    center = fwm.branch_center_GHz(params["opd"], -1)
    ultra_axis = fwm.probe_scan_axis_GHz(
        params["opd"], 401, 0, None, center - fwm.WINDOW_GHZ, center + fwm.WINDOW_GHZ,
        branches=(-1,))
    assert fast["response_estimator"]["method"] == fwm.RESPONSE_POLE_ADAPTIVE
    assert balanced["response_estimator"]["method"] == fwm.RESPONSE_POLE
    for raw in (fast, balanced):
        np.testing.assert_array_equal(raw["probe_axis_GHz"], ultra_axis)
        for key in ("G_s", "G_c", "S_dB"):
            assert np.all(np.isfinite(raw[key]))
        assert raw["phase_detail"] == fwm.PHASE_ULTRA
        assert raw["hardened_noise"]["absorption_model"].startswith("analytic Voigt")
        view = scheme.observables(raw, params, include_figures=False)
        table = next(t for t in view["tables"] if t["title"] == "Model diagnostics")
        assert "Response estimator" in table["markdown"]


def test_tier_table_keeps_ultra_and_maps_earlier_labels():
    assert fwm.FWM_FIDELITY[fwm.FIDELITY_ULTRA] == dict(
        coarse_points=401, velocity_step=1.0, velocity_cutoff=4.0,
        phase_detail=fwm.PHASE_ULTRA)
    for old in ("Fast  (~4 s)", "Balanced  (~6 s)", "Fast  (~3 s)"):
        assert fwm.normalize_fidelity(old) == fwm.FIDELITY_FAST
    for old in ("Balanced  (~12 s)", "High fidelity  (~20 s)", "Fine  (~20 s)"):
        assert fwm.normalize_fidelity(old) == fwm.FIDELITY_BALANCED
    for tier in (fwm.FIDELITY_FAST, fwm.FIDELITY_BALANCED):
        settings = fwm.FWM_FIDELITY[tier]
        assert settings["coarse_points"] == 401
        assert (settings["velocity_step"], settings["velocity_cutoff"]) == (1.0, 4.0)
        assert "full_scan" in settings


@pytest.mark.parametrize("P_pump, P_seed", ((0.6, 8e-6), (0.05, 2e-4), (1.2, 1e-6)))
def test_compiled_scorer_matches_the_python_segmented_gain(P_pump, P_seed):
    if not kernels.available():
        pytest.skip("compiled path unavailable")
    from gabes import constants, core, observables
    center = fwm.branch_center_GHz(0.9, -1)
    spec = _small(fwm.RESPONSE_POLE, fwm.PHASE_ULTRA)
    rng = np.random.default_rng(5)
    base = 2e-10 * (rng.standard_normal((41, 4)) + 1j * rng.standard_normal((41, 4)))
    chi = tuple(base[:, q] for q in range(4))
    probe = spec["probe_axis_GHz"]
    _, k_probe, k_conj = fwm.seeded_option_a_wavenumbers(0.9, probe)
    dk = fwm.seeded_phase_mismatch_z(0.9, probe)
    N_atoms = spec["N_atoms"]
    profile = fwm._gaussian_overlap_profile(64, fwm.L_CELL, fwm.W_PUMP, fwm.W_PROBE, 0.32)
    G_s, G_c, _, _ = fwm._ultra_segmented_gain(
        chi[0], chi[2], chi[1], chi[3], k_probe, k_conj, fwm.L_CELL, N_atoms, 0.1,
        dk, np.ones(64), profile, P_pump, P_seed, conjugate_power_ratio=1.3)
    M = observables._gain_matrix_from_chi(
        chi[0], chi[2], chi[1], chi[3], k_probe, k_conj, N_atoms,
        constants.DIPOLE_D1, 0.1, delta_k_z=dk)
    got_s, got_c = kernels.segmented_depletion_gain(
        np.ascontiguousarray(M), fwm.L_CELL / 64, np.ascontiguousarray(profile),
        P_pump, P_seed, 1.3, core._EXP_ARG_CLAMP)
    assert center < probe.max()
    # Rows past the Manley-Rowe budget feed the pump debit back into the next
    # step's coupling, which amplifies ~1e-16 round-off to ~1e-10.
    np.testing.assert_allclose(got_s, G_s, rtol=1e-8)
    np.testing.assert_allclose(got_c, G_c, rtol=1e-8, atol=1e-300)


@pytest.mark.parametrize("T, D", ((394.15, 0.9), (408.15, -2.2), (423.15, 2.5), (333.15, 0.3)))
def test_analytic_absorption_diagnostics_track_the_numerical_tables(T, D):
    factors = fwm._reference_population_factors()
    assert factors is not None and set(factors) == {2, 3}
    scatter, od = fwm._pump_scatter_noise(D, T, fwm.L_CELL, 0.1)
    scatter_a, od_a = fwm._pump_scatter_noise_analytic(
        D, T, fwm.L_CELL, 0.1, population_factors=factors)
    assert abs(scatter_a - scatter) < 1e-4
    assert abs(od_a - od) <= 3e-3 * max(od, 1e-3)


def test_reference_population_convention_is_classified_exactly():
    from gabes import hyperfine
    factors = fwm._reference_population_factors()
    for Fg, factor in factors.items():
        assert factor in (1.0, hyperfine.GROUND_POP[Fg])


def test_unrecognized_reference_normalization_keeps_numerical_diagnostics(monkeypatch):
    from gabes.schemes import absorption
    original = absorption._hyperfine_alpha

    def rescaled(scan, params, line_strength=None):
        alpha, components, info = original(scan, params, line_strength)
        return 0.8 * alpha, {key: 0.8 * value for key, value in components.items()}, info

    fwm._reference_population_factors.cache_clear()
    monkeypatch.setattr(absorption, "_hyperfine_alpha", rescaled)
    try:
        assert fwm._reference_population_factors() is None
        spec = _small(fwm.RESPONSE_POLE, fwm.PHASE_ULTRA)
        assert spec["hardened_noise"]["absorption_model"].startswith("numerical")
        assert spec["response_estimator"]["absorption_diagnostics"].startswith("numerical")
    finally:
        monkeypatch.undo()
        fwm._reference_population_factors.cache_clear()
