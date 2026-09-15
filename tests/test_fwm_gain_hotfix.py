"""Gain-closure acceptance, trust boundaries, and pre-hotfix numerical fixtures."""
import json
from pathlib import Path

import numpy as np
import pytest

from gabes import fwm_gain_closure
from gabes.schemes import fwm
from analysis.fwm_gain_hotfix.snapshot_convention import snapshot_directory


HERE = Path(__file__).resolve().parents[1] / "analysis/fwm_gain_hotfix"


@pytest.fixture(scope="module")
def before():
    return snapshot_directory()


@pytest.mark.parametrize("factors,directory", (
    ({2: 1.0, 3: 1.0}, "before"),
    ({2: 5.0 / 12.0, 3: 7.0 / 12.0}, "before_parent")))
def test_snapshot_matches_recorded_absorption_convention(monkeypatch, factors, directory):
    monkeypatch.setattr(fwm, "_reference_population_factors", lambda: factors)
    assert snapshot_directory() == HERE / directory


@pytest.mark.parametrize("factors", (None, {}, {2: 1.0}, {2: 1.0, 3: 7.0 / 12.0}))
def test_snapshot_rejects_unknown_or_mixed_absorption_convention(monkeypatch, factors):
    monkeypatch.setattr(fwm, "_reference_population_factors", lambda: factors)
    with pytest.raises(RuntimeError, match="No pre-hotfix snapshot"):
        snapshot_directory()


def _assert_portable_snapshot(raw, expected):
    # Frozen arrays came from Windows/LAPACK/Numba. Cross-build eigensolver
    # roundoff is not a physics regression. The same-environment validation
    # runner and same-process bypass test below retain exact equality.
    for key in expected.files:
        actual, reference = np.asarray(raw[key]), expected[key]
        assert actual.shape == reference.shape, key
        assert actual.dtype.kind == reference.dtype.kind, key
        if reference.dtype.kind in "biu":
            np.testing.assert_array_equal(actual, reference, err_msg=key)
        elif key in ("S_dB", "gain_referred_noise_dB"):
            np.testing.assert_allclose(actual, reference, rtol=0, atol=1e-5,
                                       err_msg=key)
        else:
            np.testing.assert_allclose(actual, reference, rtol=1e-7, atol=1e-10,
                                       err_msg=key)


@pytest.fixture(scope="module")
def gold_outputs():
    scheme = fwm.FWMScheme()
    outputs = {}
    for tier in (fwm.FIDELITY_FAST, fwm.FIDELITY_BALANCED):
        for enabled in (False, True):
            params = dict(scheme.defaults(), resolution=tier,
                          gain_closure_enabled=enabled)
            outputs[tier, enabled] = scheme.compute(params)
    return outputs


@pytest.mark.parametrize("tier", (fwm.FIDELITY_FAST, fwm.FIDELITY_BALANCED))
def test_gold_gain_estimate_and_physical_squeezing_boundary(gold_outputs, tier):
    raw = gold_outputs[tier, True]
    reference = gold_outputs[tier, False]
    gain = fwm.operating_point(raw, -8.0)["G_s"]
    assert 10.0 <= gain <= 20.0
    assert fwm.operating_point(reference, -8.0)["G_s"] > 300.0
    assert raw["physical_squeezing_dB"] is None
    assert raw["claim_gate"]["quantitative_gain_supported"] is False
    assert raw["eta"] == reference["eta"]
    closure = raw["gain_closure"]
    assert closure["applied"] is True
    assert 0.0 < closure["coupling_multiplier"] < 1.0
    assert closure["application"] == "chi_sc and chi_cs before Maxwell propagation"
    assert closure["calibration"]["target_probe_power_gain"] == 15.5
    assert closure["calibration"]["raw_rounded_probe_power_ratio"] == 111.0 / 8.0
    assert "P_probe_out / P_seed_in" in closure["calibration"]["gain_convention"]
    assert closure["effective_participation"]["uncertainty"] is None
    assert closure["transverse_participation"]["uncertainty"] is None
    assert closure["physical_squeezing_validated"] is False
    assert closure["independently_calibrated"] is False
    assert "mixed conditional comparisons" in closure["held_out_validation_status"]
    assert "underpredicted" in closure["held_out_validation_status"]
    assert reference["gain_closure"]["applied"] is False
    for name in ("G_s", "G_c", "S_dB"):
        assert np.all(np.isfinite(raw[name]))


@pytest.mark.parametrize("tier,name", (
    (fwm.FIDELITY_FAST, "fast_default"),
    (fwm.FIDELITY_BALANCED, "balanced_default")))
def test_switch_off_reproduces_saved_pre_hotfix_arrays(gold_outputs, before, tier, name):
    with np.load(before / f"{name}.npz") as expected:
        _assert_portable_snapshot(gold_outputs[tier, False], expected)


@pytest.mark.parametrize("name", ("ultra_default", "ultra_detuned"))
def test_ultra_matches_all_saved_pre_hotfix_arrays(before, name):
    manifest = json.loads((before / "capture.json").read_text(encoding="utf-8"))
    params = dict(manifest["cases"][name]["params"], gain_closure_enabled=True)
    raw = fwm.FWMScheme().compute(params)
    assert raw["gain_closure"]["applied"] is False
    assert raw["gain_closure"]["coupling_multiplier"] == 1.0
    with np.load(before / f"{name}.npz") as expected:
        _assert_portable_snapshot(raw, expected)


def _tiny_spectrum(**overrides):
    center = fwm.branch_center_GHz(0.9, -1)
    kwargs = dict(T=394.15, P_pump=0.6, P_probe=8e-6,
                  coarse_points=3, fine_points=0,
                  scan_min=center - 0.01, scan_max=center + 0.01,
                  velocity_step=200.0, velocity_cutoff=0.1, branch=-1)
    kwargs.update(overrides)
    return fwm.compute_spectrum(0.9, **kwargs)


@pytest.mark.parametrize("fidelity", (None, fwm.FIDELITY_ULTRA))
def test_untiered_and_ultra_api_ignore_calibration_switch(fidelity):
    kwargs = dict(model_fidelity=fidelity, phase_detail=fwm.PHASE_ULTRA)
    off = _tiny_spectrum(gain_closure_enabled=False, **kwargs)
    on = _tiny_spectrum(gain_closure_enabled=True, **kwargs)
    for key, value in off.items():
        if isinstance(value, np.ndarray) and value.dtype.kind in "biufc":
            np.testing.assert_array_equal(on[key], value, err_msg=key)


@pytest.mark.parametrize("fidelity", ("Fast  (~4 s)", "Balanced  (~12 s)"))
def test_saved_fast_tier_aliases_keep_correction(fidelity):
    raw = _tiny_spectrum(model_fidelity=fidelity, gain_closure_enabled=True)
    assert raw["gain_closure"]["applied"] is True


def test_source_gain_independent_of_detection_efficiency(gold_outputs):
    scheme = fwm.FWMScheme()
    params = dict(scheme.defaults(), resolution=fwm.FIDELITY_BALANCED,
                  gain_closure_enabled=True, detection_eff_pct=40.0)
    changed = scheme.compute(params)
    reference = gold_outputs[fwm.FIDELITY_BALANCED, True]
    assert changed["eta"] != reference["eta"]
    for key in ("G_s", "G_c", "G_s_smallsignal", "G_c_smallsignal"):
        np.testing.assert_array_equal(changed[key], reference[key])
    assert not np.array_equal(changed["S_dB"], reference["S_dB"])


@pytest.mark.parametrize("enabled", (False, True))
@pytest.mark.parametrize("tier", (
    fwm.FIDELITY_FAST, fwm.FIDELITY_BALANCED, fwm.FIDELITY_ULTRA))
def test_full_scan_routes_closure_choice(monkeypatch, enabled, tier):
    captured = {}
    monkeypatch.setattr(fwm, "full_spectrum",
                        lambda *args, **kwargs: captured.update(kwargs) or {"ok": True})
    scheme = fwm.FWMScheme()
    result = scheme.extra_views()[0].compute(dict(
        scheme.defaults(), resolution=tier, gain_closure_enabled=enabled))
    assert result == {"ok": True}
    assert captured["gain_closure_enabled"] is enabled
    assert captured["model_fidelity"] == tier


def test_untiered_full_scan_keeps_both_branches_uncalibrated(monkeypatch):
    calls = []
    monkeypatch.setattr(fwm, "compute_spectrum", lambda *args, **kwargs:
                        calls.append(kwargs) or {"captured": True})
    fwm.full_spectrum(0.9, 394.15, 600.0, 8.0, 0.74, 5.5)
    assert [call["branch"] for call in calls] == [-1, 1]
    assert all(call["model_fidelity"] is None for call in calls)


def test_switch_changes_compute_cache_and_old_params_enable_hotfix(gold_outputs):
    scheme = fwm.FWMScheme()
    assert scheme.defaults()["gain_closure_enabled"] is True
    assert "gain_closure_enabled" in scheme.recompute_keys()
    params = dict(scheme.defaults(), resolution=fwm.FIDELITY_FAST)
    params.pop("gain_closure_enabled")
    raw = scheme.compute(params)
    np.testing.assert_array_equal(raw["G_s"], gold_outputs[fwm.FIDELITY_FAST, True]["G_s"])


def test_runtime_calibration_id_resolves_to_frozen_source(gold_outputs):
    references = json.loads((HERE / "reference_points.json").read_text(encoding="utf-8"))
    ledger = gold_outputs[fwm.FIDELITY_BALANCED, True]["gain_closure"]
    assert ledger["calibration"]["reference_id"] == references["gold"]["id"]
    assert ledger["effective_participation"]["calibration_point"] == references["gold"]["id"]


@pytest.mark.parametrize("waists", ((530e-6, 330e-6, 330e-6),
                                   (450e-6, 230e-6, 400e-6)))
def test_transverse_factor_matches_independent_radial_integral(waists):
    wp, ws, wc = waists
    radius = np.linspace(0.0, 8.0 * max(waists), 20001)
    mode_product = np.exp(-radius**2 / ws**2) * np.exp(-radius**2 / wc**2)
    pump_intensity_fraction = np.exp(-2.0 * radius**2 / wp**2)
    numerical = (np.trapezoid(radius * mode_product * pump_intensity_fraction, radius)
                 / np.trapezoid(radius * mode_product, radius))
    assert fwm_gain_closure.transverse_participation(*waists) == pytest.approx(
        numerical, rel=2e-6)


def test_correction_preserves_diagonal_drift_and_changes_both_cross_couplings():
    # Unequal complex entries detect coefficient-order or conjugation mistakes.
    chi = tuple(np.array([value], dtype=complex)
                for value in (1.0 + 2j, 3.0 - 4j, 5.0 + 6j, 7.0 - 8j))
    corrected = fwm_gain_closure.apply(chi, 0.6)
    assert corrected[0] is chi[0]
    assert corrected[3] is chi[3]
    np.testing.assert_array_equal(corrected[1], 0.6 * chi[1])
    np.testing.assert_array_equal(corrected[2], 0.6 * chi[2])
    assert fwm_gain_closure.apply(chi, 1.0) is chi


def test_fitted_participation_is_fixed_across_unvalidated_waist_extrapolation():
    narrow = fwm_gain_closure.provenance(
        enabled=True, eligible=True, w_pump=400e-6, w_probe=230e-6)
    wide = fwm_gain_closure.provenance(
        enabled=True, eligible=True, w_pump=660e-6, w_probe=430e-6)
    assert narrow["coupling_multiplier"] == wide["coupling_multiplier"]
    assert (narrow["transverse_participation"]["value"]
            != wide["transverse_participation"]["value"])


@pytest.mark.parametrize("bad", (0.0, -1e-3, float("nan"), float("inf")))
def test_mode_projection_rejects_invalid_waist(bad):
    with pytest.raises(ValueError, match="finite and positive"):
        fwm_gain_closure.transverse_participation(bad, 330e-6)
