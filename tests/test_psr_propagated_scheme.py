"""Native propagated-PSR UI, physical input scope and result consistency."""
import json

import matplotlib
import numpy as np
import pytest

matplotlib.use("Agg")

from gabes import schemes
from gabes.schemes import psr_rf, psr_propagated
from gabes_ui import controls
from gabes_ui.export import build_export_payload, collect_view_plot_payloads


@pytest.fixture(scope="module")
def propagated_case():
    scheme = schemes.get("psr")
    params = dict(scheme.defaults(), calculation=psr_rf.PROPAGATED_RF_MODEL,
                  obe_intensity_mw_cm2=12.5, obe_detuning_mhz=400.,
                  obe_azimuth_deg=18., obe_ellipticity_deg=10., temp_c=40.,
                  cell_mm=25., relaxation_hz=1e5, rf_frequency_mhz=.4,
                  isotope="Rb85", line="D2", obe_velocity_m_s=500., excess_noise_per_od=11.)
    raw = scheme.compute(params)
    assert raw["unavailable_reason"] is None
    return scheme, params, raw


def test_propagated_noise_profile_figures_and_export_agree(propagated_case):
    import matplotlib.pyplot as plt

    scheme, params, raw = propagated_case
    assert raw["propagation_convergence"]["passed"]
    assert raw["output_intensity_mw_cm2"] < .9*params["obe_intensity_mw_cm2"]
    assert len(raw["segment_diagnostics"]) == raw["completed_segments"]
    intensities = [row["intensity_mw_cm2"] for row in raw["segment_diagnostics"]]
    assert np.ptp(intensities) > 1.
    for axis in (psr_rf.RF_VIEW, psr_rf.PHASE_VIEW):
        p = dict(params, rf_view=axis)
        view = scheme.observables(raw, p)
        try:
            readout = view["psr_readout"]
            point = readout["selected_point"]
            assert point["available"]
            eigenvalues = np.linalg.eigvalsh(readout["spectrum"]["covariance"][12])
            np.testing.assert_allclose(10*np.log10(eigenvalues), [point["minimum_db"], point["maximum_db"]], atol=1e-12)
            assert scheme.headless_observables(raw, p)["psr_readout"]["selected_point"] == point
            profile = readout["pump_profile"]
            assert abs(profile["azimuth_deg"][-1]-profile["azimuth_deg"][0]) > .1
            assert abs(profile["ellipticity_deg"][-1]-profile["ellipticity_deg"][0]) > .01
            assert [entry["label"] for entry in view["figure_views"]] == ["RF noise", "Pump propagation"]
            figure = view["figure_views"][1]["figure"]
            np.testing.assert_array_equal(figure.axes[0].lines[0].get_ydata(), profile["intensity_mw_cm2"])
            np.testing.assert_array_equal(figure.axes[1].lines[0].get_ydata(), profile["azimuth_deg"])
            np.testing.assert_array_equal(figure.axes[1].lines[1].get_ydata(), profile["ellipticity_deg"])
            for entry in view["figure_views"]:
                entry["figure"].canvas.draw()
            payload = build_export_payload(scheme=scheme, raw=raw, params=p, view=view,
                plot_payloads=collect_view_plot_payloads(view, view["figure"]), comparison_payload=None,
                cache_version=scheme.cache_version, readout_cache_version="test")
            export = json.loads(json.dumps(payload, allow_nan=False))
            applied = export["params"]
            assert applied["isotope"] == "Rb87" and applied["line"] == "D1" and not applied["doppler"]
            assert {"obe_velocity_m_s", "excess_noise_per_od", "thermal_accuracy", "prop_thermal_accuracy"}.isdisjoint(applied)
            assert applied["prop_velocity_m_s"] == 0.
            assert applied["obe_azimuth_deg"] == 18. and applied["obe_ellipticity_deg"] == 10.
            assert export["saved_control_state"]["obe_velocity_m_s"] == 500.
            assert export["view"]["psr_readout"]["output_polarization"]["polarization_phase_chart"] == "x"
        finally:
            for entry in view["figure_views"]:
                plt.close(entry["figure"])


def test_propagated_controls_exclude_hidden_velocity_and_thermal_inputs(monkeypatch):
    scheme = schemes.get("psr")
    for velocity in ("single", "thermal"):
        params = dict(scheme.defaults(), calculation=psr_rf.PROPAGATED_RF_MODEL, prop_velocity_model=velocity)
        monkeypatch.setattr(controls.st, "session_state", {controls.skey("psr", k): v for k, v in params.items()})
        visible = {sp.name for sp in scheme.param_schema() if controls.param_visible("psr", sp)}
        assert {"prop_velocity_model", "propagation_accuracy", "obe_azimuth_deg", "obe_ellipticity_deg"} <= visible
        assert {"obe_velocity_m_s", "thermal_accuracy", "excess_noise_per_od"}.isdisjoint(visible)
        assert ("prop_velocity_m_s" in visible) == (velocity == "single")
        assert ("prop_thermal_accuracy" in visible) == (velocity == "thermal")
        applied = scheme._active_atomic_params(params)
        assert applied["doppler"] is (velocity == "thermal")
        assert ("prop_velocity_m_s" in applied) == (velocity == "single")
    assert {"rf_view", "after_cell_eff_pct", "qe_pct"}.isdisjoint(scheme.recompute_keys())


def test_detection_reuses_path_but_length_changes_recompute_pump(propagated_case, monkeypatch):
    scheme, params, raw = propagated_case
    original = psr_propagated.solve_psr_propagation
    calls = []
    def spy(*args, **kwargs):
        calls.append(args[2])
        return original(*args, **kwargs)
    monkeypatch.setattr(psr_propagated, "solve_psr_propagation", spy)
    p = dict(params, qe_pct=80., after_cell_eff_pct=40., rf_view=psr_rf.PHASE_VIEW)
    view = scheme.headless_observables(scheme.compute(p), p)
    assert calls == []
    noise = view["psr_readout"]["spectrum"]
    np.testing.assert_allclose(noise["covariance"], .32*noise["source_covariance"]+.68*np.eye(2), atol=1e-13)
    short = scheme.compute(dict(params, cell_mm=12.5))
    assert calls == [.0125]
    assert short["output_intensity_mw_cm2"] > raw["output_intensity_mw_cm2"]


def test_unconverged_path_is_a_diagnostic_profile_without_noise(propagated_case, monkeypatch):
    scheme, params, _ = propagated_case
    monkeypatch.setitem(psr_propagated._LENGTH_SETTINGS, "Quick",
                        dict(segments=1, max_segments=1, rtol=.02, max_atomic_evaluations=20))
    psr_propagated._propagated.cache_clear()
    p = dict(params, propagation_accuracy="Quick")
    raw = scheme.compute(p)
    view = scheme.headless_observables(raw, p)
    assert "not converged" in raw["unavailable_reason"]
    assert not view["psr_readout"]["selected_point"]["available"]
    assert view["psr_readout"]["spectrum"] is None
    assert not view["psr_readout"]["pump_profile"]["converged"]
    assert next(m["value"] for m in view["metrics"] if m["label"] == "Output pump intensity") == "unavailable"


def test_propagated_zero_input_and_zero_length():
    scheme = schemes.get("psr")
    params = dict(scheme.defaults(), calculation=psr_rf.PROPAGATED_RF_MODEL)
    raw = scheme.compute(dict(params, obe_intensity_mw_cm2=0.))
    assert "bright input carrier" in raw["unavailable_reason"]
    assert not scheme.headless_observables(raw, params)["psr_readout"]["selected_point"]["available"]
    raw = scheme.compute(dict(params, cell_mm=0.))
    view = scheme.headless_observables(raw, params)
    np.testing.assert_allclose(view["psr_readout"]["spectrum"]["minimum_db"], 0., atol=1e-12)
    assert raw["output_intensity_mw_cm2"] == pytest.approx(params["obe_intensity_mw_cm2"])
