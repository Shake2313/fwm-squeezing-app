"""End-to-end checks of the PSR scheme's single-velocity RF calculation."""
from collections import Counter
import json

import matplotlib
import numpy as np
import pytest

matplotlib.use("Agg")

from gabes import schemes
from gabes.schemes import psr_rf
from gabes_ui import controls
from gabes_ui.export import build_export_payload, collect_view_plot_payloads


def rf_params(**changes):
    params = schemes.get("psr").defaults()
    # Settings saved in the weak and S1 modes must not become RF inputs.
    params.update(calculation=psr_rf.RF_MODEL, line="D2", isotope="Rb85",
                  obe_azimuth_deg=29., obe_ellipticity_deg=12., excess_noise_per_od=11.)
    return dict(params, **changes)


@pytest.fixture(scope="module")
def rf_case():
    scheme = schemes.get("psr")
    params = rf_params(rf_frequency_mhz=1.234)
    raw = scheme.compute({key: params[key] for key in scheme.recompute_keys()})
    assert raw["unavailable_reason"] is None
    return scheme, params, raw


def test_rf_phase_spectrum_figures_and_export_use_the_same_selected_point(rf_case):
    import matplotlib.pyplot as plt

    scheme, params, raw = rf_case
    views = []
    try:
        for mode in (psr_rf.RF_VIEW, psr_rf.PHASE_VIEW):
            p = dict(params, rf_view=mode)
            headless = scheme.headless_observables(raw, p)
            view = scheme.observables(raw, p)
            views.append(view)
            assert headless["figure"] is None and headless["metrics"] == view["metrics"]
            readout = view["psr_readout"]
            point, spectrum, phase = (readout[k] for k in ("selected_point", "spectrum", "phase_scan"))
            hz = spectrum["frequencies_hz"]
            assert point["available"] and point["frequency_hz"] == pytest.approx(1.234e6)
            selected, = np.flatnonzero(np.isclose(hz, point["frequency_hz"], rtol=1e-12))
            eigenvalues = np.linalg.eigvalsh(spectrum["covariance"][selected])
            assert point["minimum_db"] == pytest.approx(10*np.log10(eigenvalues[0]), abs=1e-12)
            assert point["maximum_db"] == pytest.approx(10*np.log10(eigenvalues[1]), abs=1e-12)
            angles = np.deg2rad(phase["phase_degrees"])
            for key, theta in (("noise_a", angles), ("noise_b", angles+np.pi/2)):
                rows = np.column_stack((np.cos(theta), np.sin(theta)))
                expected = np.einsum("pi,ij,pj->p", rows, spectrum["covariance"][selected], rows)
                np.testing.assert_allclose(phase[key], expected, atol=1e-13)
            axes = view["figure"].axes
            if mode == psr_rf.RF_VIEW:
                np.testing.assert_allclose(axes[0].lines[0].get_ydata(), spectrum["minimum_db"])
                np.testing.assert_allclose(axes[0].lines[1].get_ydata(), spectrum["maximum_db"])
                np.testing.assert_allclose(axes[1].lines[0].get_xdata(), hz/1e6)
                np.testing.assert_allclose(axes[1].lines[0].get_ydata(), spectrum["min_angle_deg"])
            else:
                np.testing.assert_allclose(axes[0].lines[0].get_xdata(), phase["phase_degrees"])
                np.testing.assert_allclose(axes[0].lines[0].get_ydata(), phase["noise_a_db"])
                np.testing.assert_allclose(axes[0].lines[1].get_ydata(), phase["noise_b_db"])
                assert np.all(np.asarray(axes[0].lines[-1].get_ydata()) == 0.)
            view["figure"].canvas.draw()
        payload = build_export_payload(
            scheme=scheme, raw=raw, params=p, view=view,
            plot_payloads=collect_view_plot_payloads(view, view["figure"]),
            comparison_payload=None, cache_version=scheme.cache_version, readout_cache_version="test",
        )
        exported = json.loads(json.dumps(payload, allow_nan=False))
        assert exported["params"]["rf_frequency_mhz"] == 1.234
        assert exported["params"]["rf_view"] == psr_rf.PHASE_VIEW
        assert exported["params"]["isotope"] == "Rb87" and exported["params"]["line"] == "D1"
        assert exported["params"]["doppler"] is False
        assert "excess_noise_per_od" not in exported["params"]
        assert "obe_ellipticity_deg" not in exported["params"] and "obe_azimuth_deg" not in exported["params"]
        assert exported["saved_control_state"] == p
        np.testing.assert_allclose(exported["raw"]["conditions"]["jones_xy"]["real"], [1., 0.])
        assert exported["view"]["psr_readout"]["selected_point"] == point
    finally:
        for view in views:
            plt.close(view["figure"])


def test_rf_controls_and_cache_follow_physical_dependencies(monkeypatch):
    scheme = schemes.get("psr")
    params = rf_params()
    monkeypatch.setattr(controls.st, "session_state", {controls.skey("psr", "calculation"): psr_rf.RF_MODEL})
    visible = {sp.name for sp in scheme.param_schema() if controls.param_visible("psr", sp)}
    assert {"rf_frequency_mhz", "rf_view", "obe_detuning_mhz", "cell_mm", "qe_pct"} <= visible
    assert {"view", "isotope", "line", "doppler", "excess_noise_per_od", "obe_ellipticity_deg"}.isdisjoint(visible)
    assert {"rf_view", "qe_pct", "after_cell_eff_pct"}.isdisjoint(scheme.recompute_keys())

    calls = Counter()
    def spy(name):
        original = getattr(psr_rf, name)
        def counted(*args, **kwargs):
            calls[name] += 1
            return original(*args, **kwargs)
        return counted

    names = ("solve_psr_pump", "psr_local_generator", "psr_uniform_transfer")
    for name in names:
        monkeypatch.setattr(psr_rf, name, spy(name))
    psr_rf._pump.cache_clear()
    psr_rf._local.cache_clear()
    try:
        raw = scheme.compute(params)
        assert tuple(calls[n] for n in names) == (1, 1, 1)
        changed = dict(params, rf_view=psr_rf.PHASE_VIEW, qe_pct=80., after_cell_eff_pct=40.)
        readout = scheme.headless_observables(raw, changed)["psr_readout"]
        spectrum = readout["spectrum"]
        assert readout["detection_efficiency"] == pytest.approx(.32)
        np.testing.assert_allclose(spectrum["covariance"],
                                   .32*spectrum["source_covariance"]+.68*np.eye(2), atol=1e-14)
        assert tuple(calls[n] for n in names) == (1, 1, 1)
        assert scheme.compute(dict(params, cell_mm=35.))["unavailable_reason"] is None
        assert tuple(calls[n] for n in names) == (1, 1, 2)
        assert scheme.compute(dict(params, rf_frequency_mhz=2.4))["unavailable_reason"] is None
        assert tuple(calls[n] for n in names) == (1, 2, 3)
        assert scheme.compute(dict(params, obe_intensity_mw_cm2=.2))["unavailable_reason"] is None
        assert tuple(calls[n] for n in names) == (2, 3, 4)
    finally:
        psr_rf._local.cache_clear()
        psr_rf._pump.cache_clear()


def test_uniform_pump_boundary_and_absent_stokes_carrier(rf_case):
    scheme, params, cached_raw = rf_case
    # Set the boundary length from an actual S1 absorption at the same pump.
    params = dict(params, temp_c=100.)
    reference = scheme.compute(dict(params, cell_mm=1.))
    limit_mm = reference["uniform_pump_od_limit"] / reference["pump_optical_depth"]
    assert 0 < limit_mm < 200
    for factor, available in ((1-1e-6, True), (1+1e-6, False)):
        p = dict(params, cell_mm=limit_mm*factor)
        raw = scheme.compute(p)
        readout = scheme.headless_observables(raw, p)["psr_readout"]
        assert readout["selected_point"]["available"] is available
        if not available:
            assert "Uniform-pump" in readout["selected_point"]["reason"]
            assert readout["spectrum"] is None and readout["phase_scan"] is None
    no_carrier = dict(params, obe_intensity_mw_cm2=0., cell_mm=1.)
    for raw, p in ((scheme.compute(no_carrier), no_carrier), (cached_raw, dict(params, qe_pct=0.))):
        view = scheme.headless_observables(raw, p)
        assert not view["psr_readout"]["selected_point"]["available"]
        assert all(m["value"] == "unavailable" for m in view["metrics"][:2])
        assert view["psr_readout"]["phase_scan"] is None


def test_zero_length_is_measurable_vacuum_with_undefined_squeezing_angle():
    scheme, p = schemes.get("psr"), rf_params(cell_mm=0.)
    readout = scheme.headless_observables(scheme.compute(p), p)["psr_readout"]
    assert readout["selected_point"]["available"]
    assert readout["selected_point"]["min_angle_deg"] is None
    np.testing.assert_allclose(readout["spectrum"]["minimum_db"], 0., atol=1e-13)
    np.testing.assert_allclose(readout["phase_scan"]["noise_b_db"], 0., atol=1e-13)
