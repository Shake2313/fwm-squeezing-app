"""Thermal RF integration through native controls, readout, caching and export."""
import json

import matplotlib
import numpy as np
import pytest

matplotlib.use("Agg")

from gabes import schemes, psr_thermal
from gabes.schemes import psr_rf
from gabes_ui import controls
from gabes_ui.export import build_export_payload, collect_view_plot_payloads


@pytest.fixture(scope="module")
def thermal_case():
    scheme = schemes.get("psr")
    params = dict(scheme.defaults(), calculation=psr_rf.THERMAL_RF_MODEL,
                  isotope="Rb85", line="D2", obe_velocity_m_s=593.,
                  obe_azimuth_deg=29., obe_ellipticity_deg=12., excess_noise_per_od=11.)
    # Keep the same RF endpoints and resonance anchors as the 25-point UI,
    # with three directly solved RF samples for bounded integration-test cost.
    with pytest.MonkeyPatch.context() as patch:
        patch.setattr(psr_rf, "_RF_RATIOS", np.array([.1, 1., 10.]))
        psr_rf._thermal.cache_clear()
        raw = scheme.compute(params)
        assert raw["unavailable_reason"] is None
        yield scheme, params, raw
        psr_rf._thermal.cache_clear()


def test_thermal_selected_noise_figures_and_export(thermal_case):
    import matplotlib.pyplot as plt

    scheme, params, raw = thermal_case
    assert raw["conditions"]["doppler"] is True
    assert "velocity_m_s" not in raw["conditions"]
    assert "effective_detuning_mhz" not in raw["conditions"]
    assert raw["pump_model_id"] == psr_thermal.PUMP_MODEL_ID
    assert raw["thermal_convergence"]["passed"]
    for mode in (psr_rf.RF_VIEW, psr_rf.PHASE_VIEW):
        p = dict(params, rf_view=mode)
        view = scheme.observables(raw, p)
        try:
            readout = view["psr_readout"]
            assert readout["model_id"] == psr_thermal.MODEL_ID
            point, spectrum = readout["selected_point"], readout["spectrum"]
            assert point["available"]
            values = np.linalg.eigvalsh(spectrum["covariance"][1])
            np.testing.assert_allclose(10*np.log10(values), [point["minimum_db"], point["maximum_db"]], atol=1e-12)
            assert scheme.headless_observables(raw, p)["psr_readout"]["selected_point"] == point
            view["figure"].canvas.draw()
            export = build_export_payload(
                scheme=scheme, raw=raw, params=p, view=view,
                plot_payloads=collect_view_plot_payloads(view, view["figure"]),
                comparison_payload=None, cache_version=scheme.cache_version, readout_cache_version="test")
            payload = json.loads(json.dumps(export, allow_nan=False))
            assert payload["params"]["doppler"] is True
            assert payload["params"]["isotope"] == "Rb87" and payload["params"]["line"] == "D1"
            assert {"obe_velocity_m_s", "obe_ellipticity_deg", "excess_noise_per_od"}.isdisjoint(payload["params"])
            assert payload["saved_control_state"]["obe_velocity_m_s"] == 593.
            assert payload["raw"]["thermal_performance"]["class_evaluations"] > 1
            assert len(payload["raw"]["velocity_quadrature"]["history"]) >= 3
            assert "Initial integration time" in view["tables"][0]["markdown"]
        finally:
            plt.close(view["figure"])


def test_thermal_control_visibility_and_length_detection_cache(thermal_case, monkeypatch):
    scheme, params, raw = thermal_case
    state = {controls.skey("psr", key): value for key, value in params.items()}
    monkeypatch.setattr(controls.st, "session_state", state)
    visible = {spec.name for spec in scheme.param_schema() if controls.param_visible("psr", spec)}
    assert {"thermal_accuracy", "rf_view", "rf_frequency_mhz", "cell_mm", "obe_intensity_mw_cm2"} <= visible
    assert {"obe_velocity_m_s", "resolution", "doppler", "obe_ellipticity_deg"}.isdisjoint(visible)
    assert "thermal_accuracy" in scheme.recompute_keys()
    assert {"rf_view", "qe_pct", "after_cell_eff_pct"}.isdisjoint(scheme.recompute_keys())

    def forbidden(*args, **kwargs):
        raise AssertionError("length or downstream detection must reuse integrated atomic kernels")
    monkeypatch.setattr(psr_thermal, "solve_psr_thermal", forbidden)
    short = scheme.compute(dict(params, cell_mm=0.))
    assert short["thermal_convergence"]["passed"]
    vacuum = scheme.headless_observables(short, params)["psr_readout"]
    np.testing.assert_allclose(vacuum["spectrum"]["minimum_db"], 0., atol=1e-12)
    changed = dict(params, qe_pct=80., after_cell_eff_pct=40., rf_view=psr_rf.PHASE_VIEW)
    spectrum = scheme.headless_observables(raw, changed)["psr_readout"]["spectrum"]
    np.testing.assert_allclose(spectrum["covariance"], .32*spectrum["source_covariance"]+.68*np.eye(2), atol=1e-13)


def test_unconverged_thermal_integral_never_becomes_a_noise_curve(thermal_case, monkeypatch):
    scheme, params, _ = thermal_case
    monkeypatch.setattr(psr_thermal, "thermal_convergence", lambda *a, **kw:
                        dict(passed=False, status="velocity-class evaluation budget reached"))
    raw = scheme.compute(params)
    view = scheme.headless_observables(raw, params)
    assert "not converged" in raw["unavailable_reason"]
    assert not view["psr_readout"]["selected_point"]["available"]
    assert view["psr_readout"]["spectrum"] is None
    assert view["psr_readout"]["phase_scan"] is None
    assert not view["psr_readout"]["thermal_convergence"]["passed"]


def test_absent_thermal_carrier_skips_the_expensive_integral(monkeypatch):
    scheme = schemes.get("psr")
    def forbidden(*args, **kwargs):
        raise AssertionError("no carrier does not need a thermal noise solve")
    monkeypatch.setattr(psr_thermal, "solve_psr_thermal", forbidden)
    params = dict(scheme.defaults(), calculation=psr_rf.THERMAL_RF_MODEL, obe_intensity_mw_cm2=0.)
    raw = scheme.compute(params)
    assert "nonzero bright carrier" in raw["unavailable_reason"]
    assert not scheme.headless_observables(raw, params)["psr_readout"]["selected_point"]["available"]


def test_failed_thermal_solve_keeps_the_maxwell_source_label(monkeypatch):
    scheme = schemes.get("psr")
    def fail(*args, **kwargs):
        raise ValueError("test atomic convergence failure")
    monkeypatch.setattr(psr_rf, "_thermal", fail)
    params = dict(scheme.defaults(), calculation=psr_rf.THERMAL_RF_MODEL)
    raw = scheme.compute(params)
    view = scheme.headless_observables(raw, params)
    assert "test atomic convergence failure" in raw["unavailable_reason"]
    assert not raw["thermal_convergence"]["passed"]
    assert "Independent Maxwell velocity classes" in view["tables"][0]["markdown"]
    assert view["psr_readout"]["spectrum"] is None
