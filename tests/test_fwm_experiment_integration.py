"""Practical FWM models stay active across detuning and apparatus controls."""
from copy import deepcopy
import json

import matplotlib
import numpy as np
import pytest

matplotlib.use("Agg")

from gabes.schemes import fwm
from gabes_ui.export import build_export_payload


@pytest.fixture(scope="module")
def experimental_case():
    scheme = fwm.FWMScheme()
    params = dict(scheme.defaults(), opd=-1.4, temp_c=113.0,
                  pump_mw=300.0, probe_uw=3.0, cell_mm=12.5,
                  pump_waist_um=530.0, probe_waist_um=330.0,
                  seeded_angle_deg=0.32, transit_rate_khz=100.0,
                  resolution=fwm.FIDELITY_FAST, tpd=-124.317286437,
                  after_cell_eff_pct=94.5, qe_pct=92.0,
                  seed_classical_rin_scaled=0.0,
                  noise_model=fwm.NOISE_EXPERIMENTAL_JASPERSE,
                  experimental_gain_calibration=True)
    raw = scheme.compute(params)
    return scheme, params, raw


def test_practical_defaults_and_controls_only_change_readout():
    scheme = fwm.FWMScheme()
    specs = {spec.name: spec for spec in scheme.param_schema()}
    assert scheme.defaults()["noise_model"] == fwm.NOISE_EXPERIMENTAL_JASPERSE
    assert scheme.defaults()["experimental_gain_calibration"] is False
    for name in ("noise_model", "experimental_gain_calibration"):
        assert specs[name].recompute is False
        assert name not in scheme.recompute_keys()
    assert fwm.NOISE_PASSIVE_L1 in specs["noise_model"].choices


def test_experimental_nf_fills_the_screenshot_gap_and_preserves_cached_model(experimental_case):
    scheme, params, raw = experimental_case
    original = deepcopy(raw)
    view = scheme.headless_observables(raw, params)
    assert view["noise_model"] == fwm.NOISE_EXPERIMENTAL_JASPERSE
    assert view["experimental_gain_calibration"]["applied"]
    assert view["experimental_jasperse_readout"] is not None
    point = view["noise_operating_point"]
    assert np.isfinite(point["NF"]) and point["NF"] > 0.0
    assert point["NF_dB"] == pytest.approx(10*np.log10(point["NF"]), abs=1e-12)
    assert view["metrics"][0]["value"] == f"{point['NF_dB']:.2f} dB"
    assert not view["passive_l1_operating_point"]["passive_l1_available"]
    for key in ("G_s", "G_c", "T_canonical_small_signal", "gain_referred_noise_dB"):
        np.testing.assert_array_equal(raw[key], original[key])
    assert raw["claim_gate"] == original["claim_gate"]


def test_extrapolated_selected_points_keep_gain_and_noise_available(experimental_case):
    scheme, params, raw = experimental_case
    for delta in (-200.0, -40.0, 40.0, 200.0):
        view = scheme.headless_observables(raw, dict(params, tpd=delta))
        calibration = view["experimental_gain_calibration"]
        assert calibration["applied"]
        point = calibration["selected_point"]
        assert point["available"]
        assert np.isfinite(point["G_s"]) and point["G_s"] >= 0.0
        assert np.isfinite(point["G_c"]) and point["G_c"] >= 0.0
        assert np.isfinite(view["noise_operating_point"]["NF_dB"])


def test_operating_condition_changes_keep_correction_active_and_model_responsive(experimental_case):
    scheme, params, raw = experimental_case
    first = scheme.headless_observables(raw, params)
    changed_params = dict(params, pump_mw=330.0, temp_c=115.0, opd=-1.2)
    changed_raw = scheme.compute(changed_params)
    changed = scheme.headless_observables(changed_raw, changed_params)
    assert changed["experimental_gain_calibration"]["applied"]
    assert changed["experimental_gain_calibration"]["selected_point"]["available"]
    assert np.isfinite(changed["noise_operating_point"]["NF_dB"])
    assert not np.allclose(first["experimental_gain_calibration"]["G_s"],
                           changed["experimental_gain_calibration"]["G_s"])


def test_focused_plot_marker_headless_and_strict_export_agree(experimental_case):
    import matplotlib.pyplot as plt
    scheme, params, raw = experimental_case
    rendered = scheme.observables(raw, params)
    headless = scheme.headless_observables(raw, params)
    try:
        assert rendered["metrics"] == headless["metrics"]
        point = rendered["noise_operating_point"]
        axis = rendered["figure"].axes[1]
        assert float(axis.collections[0].get_offsets()[0, 1]) == pytest.approx(point["NF_dB"])
        assert np.interp(params["tpd"], *axis.lines[0].get_data()) == pytest.approx(point["NF_dB"])
        payload = build_export_payload(
            scheme=scheme, raw=raw, params=params, view=rendered,
            plot_payloads=[], comparison_payload=None,
            cache_version=scheme.cache_version, readout_cache_version="experiment-test")
        decoded = json.loads(json.dumps(payload, allow_nan=False))
        assert decoded["view"]["experimental_gain_calibration"]["applied"]
        assert decoded["view"]["experimental_jasperse_readout"] is not None
        assert decoded["view"]["noise_operating_point"]["NF"] == pytest.approx(point["NF"])
    finally:
        plt.close(rendered["figure"])


def test_detection_and_classical_noise_change_nf_without_refitting_gain(experimental_case):
    scheme, params, raw = experimental_case
    first = scheme.headless_observables(raw, params)
    changed = scheme.headless_observables(raw, dict(
        params, seed_classical_rin_scaled=20.0, qe_pct=60.0, after_cell_eff_pct=80.0))
    for key in ("G_s", "G_c"):
        np.testing.assert_array_equal(first["experimental_gain_calibration"][key],
                                      changed["experimental_gain_calibration"][key])
    assert first["noise_operating_point"]["NF"] != changed["noise_operating_point"]["NF"]
    zero = scheme.headless_observables(raw, dict(params, qe_pct=0.0))
    assert np.all(np.isnan(zero["experimental_jasperse_readout"]["NF"]))
    assert not zero["noise_operating_point"]["available"]


def test_extreme_effective_gain_mismatch_is_visible_in_selected_readout(experimental_case):
    scheme, params, raw = experimental_case
    extreme = dict(raw, G_s=np.full_like(raw["G_s"], 14.6092),
                   G_c=np.full_like(raw["G_c"], .0143289))
    view = scheme.headless_observables(extreme, dict(params, experimental_gain_calibration=False))
    point = view["noise_operating_point"]
    assert point["available"] and point["inference_approximate"]
    assert point["fit_relative_residual"] > .01
    assert "approximate gain match" in view["metrics"][0]["delta"]
    assert "Assumed Jasperse common loss" in view["tables"][0]["markdown"]
    assert "cannot reproduce this pair exactly" in view["tables"][0]["markdown"]


def test_live_cursor_uses_calibrated_gain_and_committed_nf(experimental_case):
    import matplotlib.pyplot as plt
    from gabes_ui.liveplot import preview_data
    scheme, params, raw = experimental_case
    view = scheme.observables(raw, params)
    try:
        figure = view["figure"]
        data, moving = preview_data(figure, "fwm", canonical_noise=True)
        assert data["original"]["tpd"] == pytest.approx(params["tpd"])
        gain_point = view["experimental_gain_calibration"]["selected_point"]
        gain_curve = data["curves"][0]
        assert np.interp(params["tpd"], gain_curve["x"], gain_curve["y"]) == pytest.approx(gain_point["G_s"])
        assert float(figure.axes[0].collections[0].get_offsets()[0, 1]) == pytest.approx(gain_point["G_s"])
        for axis in figure.axes:
            assert axis.lines[1] in moving
            assert axis.collections[0] in moving
        assert data["linear_noise"]
        assert data["curves"][1]["x"] == list(figure.axes[1].lines[0].get_xdata())
    finally:
        plt.close(view["figure"])


def test_full_view_routes_both_practical_models_and_branch_coordinates(experimental_case):
    import matplotlib.pyplot as plt
    from gabes.fwm_experiment_gain import calibrate_spectrum
    from gabes.fwm_experiment_noise import experimental_jasperse_readout
    scheme, params, raw = experimental_case
    plus = dict(raw, branch=1, probe_axis_GHz=(raw["probe_axis_GHz"]
                +raw["raman_center_plus_GHz"]-raw["raman_center_minus_GHz"]))
    full = {"D_GHz": params["opd"], "minus": raw, "plus": plus}
    figure = scheme.extra_views()[0].render(full, params)
    try:
        assert set(figure._gabes_experimental_gain_calibration) == {"minus", "plus"}
        assert set(figure._gabes_experimental_jasperse_readout) == {"minus", "plus"}
        for index, branch in enumerate(("minus", "plus")):
            spec = full[branch]
            gains = calibrate_spectrum(spec, enabled=True)
            noise = experimental_jasperse_readout(spec, params, gains=gains)
            assert figure._gabes_experimental_gain_calibration[branch]["applied"]
            np.testing.assert_array_equal(figure.axes[0].lines[index].get_xdata(), spec["probe_axis_GHz"])
            np.testing.assert_allclose(figure.axes[0].lines[index].get_ydata(), gains["G_s"])
            np.testing.assert_allclose(figure.axes[1].lines[index].get_ydata(), noise["NF_dB"], equal_nan=True)
    finally:
        plt.close(figure)
