"""Effective Rydberg opacity is a declared postprocessing assumption."""
import json
from unittest.mock import patch

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pytest

from gabes import observables, rydberg_opacity
from gabes.schemes.rydberg import RydbergEITScheme, _finite_scenario_envelope


@pytest.fixture(scope="module")
def reference():
    scheme = RydbergEITScheme()
    params = {mode: scheme._defaults(mode) for mode in ("EIT", "AT electrometry")}
    raw = {mode: scheme.compute(p) for mode, p in params.items()}
    return scheme, params, raw


@pytest.mark.parametrize("key", ["effective_opacity_participation",
                                 "effective_opacity_relative_uncertainty"])
@pytest.mark.parametrize("value", [-0.001, 1.001, np.nan, np.inf, -np.inf])
def test_opacity_bounds_reject_invalid_inputs(key, value):
    with pytest.raises(ValueError, match=key):
        rydberg_opacity.opacity_spec({key: value})


def test_uncertainty_is_an_assumed_scenario_range_with_boundary_deduplication():
    spec = rydberg_opacity.opacity_spec({})
    assert spec["value"] == 0.001
    assert spec["scenario_values"] == [0.0005, 0.001, 0.0015]
    assert spec["status"] == "ASSUMED"
    assert spec["calibration_status"] == "UNCALIBRATED"
    assert "not_standard_uncertainty" in spec["uncertainty_kind"]
    assert rydberg_opacity.opacity_samples({"effective_opacity_participation": 0}) == [
        {"label": "nominal", "value": 0.0}]
    upper = rydberg_opacity.opacity_samples({"effective_opacity_participation": 1})
    assert upper == [{"label": "nominal", "value": 1.0},
                     {"label": "lower", "value": 0.5}]


def test_finite_envelope_preserves_partial_scenarios_and_all_invalid_gaps():
    curves = [[np.nan, 1.0, np.nan], [2.0, 3.0, np.nan], [4.0, np.nan, np.inf]]
    with np.errstate(all="raise"):
        lower, upper, counts = _finite_scenario_envelope(curves)
    np.testing.assert_array_equal(lower, [2.0, 1.0, np.nan])
    np.testing.assert_array_equal(upper, [4.0, 3.0, np.nan])
    np.testing.assert_array_equal(counts, [2, 2, 0])


def test_opacity_controls_are_navigate_only_and_invalidate_extra_views():
    scheme = RydbergEITScheme()
    names = ("effective_opacity_participation", "effective_opacity_relative_uncertainty")
    schema = {p.name: p for p in scheme.param_schema()}
    for name in names:
        assert not schema[name].recompute
        assert name not in scheme.recompute_keys()
        assert all(name in view.param_keys for view in scheme.extra_views())


def test_default_optical_arrays_and_finite_if_match_legacy_baseline(reference):
    scheme, params, raws = reference
    for mode, raw in raws.items():
        p = params[mode]
        _, current, susceptibility = scheme._transmission(raw["chi_bar"], raw, p)
        alpha, legacy_susceptibility = observables.absorption_coefficient(
            raw["chi_bar"], raw["k_vec"], raw["N"], dipole=raw["dipole"],
            line_strength=0.001)
        legacy = observables.transmission(alpha, p["cell_mm"] * 1e-3)
        np.testing.assert_allclose(current, legacy, rtol=1e-12, atol=0)
        np.testing.assert_allclose(susceptibility, legacy_susceptibility, rtol=1e-12, atol=0)
    eit = scheme._readout(raws["EIT"], params["EIT"])
    at = scheme._readout(raws["AT electrometry"], params["AT electrometry"])
    # Recorded from the current production implementation before this change.
    assert eit["width"] == pytest.approx(1.6130585827288382, rel=1e-12)
    assert eit["T_trans"][400] == pytest.approx(0.9400166150823321, rel=1e-12)
    assert at["superhet"]["sensitivity"].total_field_asd_nv_cm_per_sqrt_hz == pytest.approx(
        6.0823872151659275, rel=1e-12)
    assert at["superhet"]["optimum_detuning_mhz"] == -2.2200000000000006


def test_live_opacity_ignores_legacy_alias_and_atom_count_fraction(reference):
    scheme, params, raws = reference
    p, raw = params["EIT"], raws["EIT"]
    base = scheme._readout(raw, p)
    counted = scheme._readout(raw, dict(p, atom_participation_fraction=0.25))
    np.testing.assert_array_equal(counted["T_trans"], base["T_trans"])
    assert counted["effective_atoms"].atoms == base["effective_atoms"].atoms * 0.25
    changed_params = dict(p, effective_opacity_participation=0.002)
    with patch.object(scheme, "compute", side_effect=AssertionError("unexpected solve")):
        changed = scheme._readout(dict(raw, ls=0.9), changed_params)
    np.testing.assert_allclose(changed["T_trans"], base["T_trans"] ** 2, rtol=1e-12)
    assert changed["effective_atoms"].atoms == base["effective_atoms"].atoms
    # A separately declared atomic strength multiplies the participation once.
    _, half, _ = scheme._transmission(raw["chi_bar"],
                                      dict(raw, atomic_line_strength=0.5), changed_params)
    np.testing.assert_allclose(half, base["T_trans"], rtol=1e-12, atol=0)


def test_opacity_and_scenario_width_reuse_atomic_response(reference):
    _, params, raws = reference
    scheme = RydbergEITScheme()
    p, raw = params["AT electrometry"], raws["AT electrometry"]
    with patch.object(scheme, "_superheterodyne_atomic_phasors",
                      wraps=scheme._superheterodyne_atomic_phasors) as solve:
        nominal = scheme._readout(raw, p)
        varied = scheme._readout(raw, dict(p, effective_opacity_participation=0.0005))
        wider = scheme._readout(raw, dict(p, effective_opacity_relative_uncertainty=1.0))
    assert solve.call_count == 1
    assert varied["superhet"]["sensitivity"].total_field_asd_nv_cm_per_sqrt_hz > (
        nominal["superhet"]["sensitivity"].total_field_asd_nv_cm_per_sqrt_hz)
    np.testing.assert_array_equal(nominal["T_trans"], wider["T_trans"])
    assert wider["superhet"]["sensitivity"] == nominal["superhet"]["sensitivity"]


def test_reused_doppler_scenario_keeps_nominal_sample_and_discloses_scope(reference):
    scheme, params, raws = reference
    p, raw = params["AT electrometry"], raws["AT electrometry"]
    nominal = scheme._readout(raw, p)["superhet"]
    indices, phasors, _ = nominal["atomic_response"]
    fixed = (indices[100:101], phasors[100:101], "nominal sampled point")
    varied = dict(p, doppler="on", effective_opacity_participation=0.002)
    x, transmission, _ = scheme._transmission(raw["chi_bar"], raw, varied)
    with patch.object(scheme, "_cached_superheterodyne_atomic_phasors",
                      side_effect=AssertionError("scenario must reuse phasor")):
        result = scheme._superheterodyne_readout(raw, varied, x, transmission,
                                                atomic_response=fixed)
    assert result["optimum_index"] == indices[100]
    assert result["optimization_scope"] == (
        "fixed nominal static-gradient sampled point (opacity scenario)")


@pytest.mark.parametrize("mode", ["EIT", "AT electrometry"])
def test_zero_opacity_has_no_fabricated_response_or_sensitivity(reference, mode):
    scheme, params, raws = reference
    p = dict(params[mode], effective_opacity_participation=0)
    with patch.object(scheme, "_cached_superheterodyne_atomic_phasors",
                      side_effect=AssertionError("zero opacity must not solve IF")):
        ro = scheme._readout(raws[mode], p)
        view = scheme.observables(raws[mode], p)
    np.testing.assert_array_equal(ro["T_trans"], np.ones_like(ro["T_trans"]))
    assert ro["superhet"] is None
    assert ro["max_slope_per_mhz"] == 0
    typed = view["rydberg_readout"]
    assert typed["if_discriminator_per_mhz"] == 0
    assert typed["if_optimum_detuning_mhz"] is None
    assert next(m["value"] for m in view["metrics"]
                if m["label"] == "IF optimum detuning") == "—"
    assert typed["finite_if"]["total_field_asd_nv_cm_per_sqrt_hz"] is None
    assert typed["finite_if"]["status"] == (
        "UNAVAILABLE_ZERO_OPACITY" if mode == "AT electrometry" else "NOT_APPLICABLE_NO_LO")
    json.dumps(typed, allow_nan=False)
    plt.close(view["figure"])


def test_headless_and_figure_readouts_share_explicit_prediction_boundary(reference):
    scheme, params, raws = reference
    p, raw = params["AT electrometry"], raws["AT electrometry"]
    headless = scheme.headless_observables(raw, p)
    plotted = scheme.observables(raw, p)
    assert headless["rydberg_readout"] == plotted["rydberg_readout"]
    assert headless["rydberg_readout"]["experimental_visibility"] == "NOT_ASSESSED"
    assert headless["prediction_status"] == "PREDICTED"
    assert headless["effective_opacity"]["calibration_status"] == "UNCALIBRATED"
    json.dumps(headless["rydberg_readout"], allow_nan=False)
    plt.close(plotted["figure"])


@pytest.mark.parametrize("value", [0.0, 0.002])
def test_export_preserves_live_opacity_provenance_and_scenario_bounds(reference, value):
    from gabes_ui.export import build_export_payload

    scheme, params, raws = reference
    p = dict(params["AT electrometry"], effective_opacity_participation=value,
             effective_opacity_relative_uncertainty=1.0)
    raw = raws["AT electrometry"]
    view = scheme.headless_observables(raw, p)
    bundle = build_export_payload(
        scheme=scheme, raw=raw, params=p, view=view, plot_payloads=[],
        comparison_payload=None, cache_version=scheme.cache_version,
        readout_cache_version="test")
    decoded = json.loads(json.dumps(bundle, allow_nan=False))
    typed = decoded["view"]["rydberg_readout"]
    assert typed["effective_opacity"]["value"] == value
    assert typed["effective_opacity"]["scenario_values"] == [0.0, value, 2*value]
    assert typed["effective_opacity"]["source"] == rydberg_opacity.opacity_spec(p)["source"]
    assert typed["calibration_status"] == "UNCALIBRATED"
    assert decoded["params"]["effective_opacity_participation"] == value
    if value == 0:
        assert typed["finite_if"]["total_field_asd_nv_cm_per_sqrt_hz"] is None


@pytest.mark.parametrize("view_index,expected_static_calls", [(0, 16), (1, 20)])
def test_scenarios_reuse_each_atomic_sweep_and_report_separate_optima(
        reference, view_index, expected_static_calls):
    _, params, _ = reference
    scheme = RydbergEITScheme()
    # Coarser deterministic grids keep this cache/dependency contract inexpensive;
    # production-grid physics remains covered by the default parity fixture.
    with patch.object(scheme, "_scan", return_value=np.linspace(-12, 12, 101) * 2*np.pi*1e6):
        with patch.object(scheme, "compute", wraps=scheme.compute) as solve:
            with patch.object(scheme, "_superheterodyne_atomic_phasors",
                              wraps=scheme._superheterodyne_atomic_phasors) as response:
                sweep = scheme.extra_views()[view_index].compute(dict(
                    params["AT electrometry"], effective_opacity_relative_uncertainty=1.0))
    assert solve.call_count == expected_static_calls
    assert response.call_count == (0 if view_index == 0 else 10)
    assert sweep["prediction_status"] == "PREDICTED"
    assert sweep["calibration_status"] == "UNCALIBRATED"
    assert "not a confidence interval" in sweep["uncertainty_scope"]
    scenarios = {s["label"]: s for s in sweep["opacity_scenarios"]}
    assert set(scenarios) == {"lower", "nominal", "upper"}
    if view_index == 0:
        assert scenarios["nominal"]["comp"] == sweep["comp"]
        assert scenarios["lower"]["comp"]["amp"] != scenarios["upper"]["comp"]["amp"]
        assert all("optima" in s for s in scenarios.values())
        assert scenarios["lower"]["optima"]["comp"]["peak_amplitude_power_uw"] is None
    else:
        assert scenarios["nominal"]["spectra"] == sweep["spectra"]
        assert scenarios["lower"]["total_sensitivity_nv_cm_sqrt_hz"] != (
            scenarios["upper"]["total_sensitivity_nv_cm_sqrt_hz"])
        for scenario in scenarios.values():
            if scenario["value"] == 0:
                assert scenario["best_sensitivity_temperature_c"] is None
                assert scenario["best_total_sensitivity_nv_cm_sqrt_hz"] is None
                assert set(scenario["finite_if_status"]) == {"UNAVAILABLE_ZERO_OPACITY"}
            else:
                assert scenario["best_sensitivity_temperature_c"] in sweep["temperatures_c"]
                assert scenario["best_total_sensitivity_nv_cm_sqrt_hz"] == min(
                    scenario["total_sensitivity_nv_cm_sqrt_hz"])
    figure = scheme.extra_views()[view_index].render(sweep)
    assert "PREDICTED" in figure._suptitle.get_text()
    if view_index == 0:
        width_axis = figure.axes[1]
        assert len(width_axis.collections) == 2
        assert all(collection.get_paths() for collection in width_axis.collections)
        assert any("finite cases only" in text.get_text() for text in width_axis.texts)
    else:
        sensitivity_axis = figure.axes[2]
        assert sensitivity_axis.collections[0].get_paths()
        assert any("no finite range bound" in text.get_text()
                   for text in sensitivity_axis.texts)
    plt.close(figure)
