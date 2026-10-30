"""Phase plots retain their own axis, operating point and exported evidence."""
import json
from types import SimpleNamespace

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pytest

from gabes_ui import export, liveplot
from gabes_ui.plot_focus import resonance_window


@pytest.fixture
def phase_figure():
    phase = np.linspace(0., 180., 181)
    noise = 2. - 1.5 * np.cos(2 * np.deg2rad(phase))
    figure, axis = plt.subplots()
    axis.plot(phase, 10 * np.log10(noise), label="Joint quadrature")
    axis.plot(phase, 10 * np.log10(4. - noise), label="Orthogonal quadrature")
    axis.axhline(0., label="Shot noise")
    axis.set_xlabel("Analysis phase [deg]")
    axis.set_ylabel("Quadrature noise [dB]")
    yield figure
    plt.close(figure)


def test_phase_plot_has_no_detuning_cursor_or_resonance_zoom(phase_figure, monkeypatch):
    view = dict(figure=phase_figure, squeezing_view="Quadrature phase scan")
    data, moving = liveplot.preview_data(phase_figure, "fwm", canonical_noise=True)
    assert data is None and moving == []
    assert resonance_window(phase_figure, "fwm", view=view) is None
    monkeypatch.setattr(liveplot, "_MOUNTABLE", True)
    monkeypatch.setattr(liveplot, "_COMPONENT", lambda **kw: pytest.fail("Phase plot mounted detuning cursor"))
    assert not liveplot.render_live_plot(phase_figure, "fwm")


def test_phase_view_cannot_be_focused_as_detuning_even_with_two_axes():
    figure, axes = plt.subplots(2, 1)
    try:
        for axis in axes:
            axis.plot(np.linspace(0., 180., 181), np.ones(181))
            axis.axvline(45.)
        assert resonance_window(figure, "fwm") is not None
        assert resonance_window(figure, "fwm", view={
            "squeezing_view": "Quadrature phase scan"}) is None
    finally:
        plt.close(figure)


def test_phase_export_retains_scan_and_fixed_gain_metadata(phase_figure, monkeypatch):
    monkeypatch.setattr(export, "app_revision", lambda: "test-revision")
    scheme = SimpleNamespace(name="fwm", title="FWM", cluster="D", defaults_version="1",
                             defaults=lambda: {"squeezing_view": "Intensity difference"})
    view = dict(
        figure=phase_figure, squeezing_view="Quadrature phase scan",
        quadrature_phase_scan=dict(phase_deg=np.array([0., 90., 180.]),
                                   variance=np.array([.5, 3.5, .5]),
                                   scope="Effective Gaussian model"),
        gain_operating_point=dict(tpd_mhz=np.float64(-8.), probe_gain=np.float64(4.)),
    )
    plots = export.collect_view_plot_payloads(view, phase_figure)
    payload = export.build_export_payload(
        scheme=scheme, raw={"kind": "seeded"},
        params={"tpd": -8., "squeezing_view": "Quadrature phase scan"}, view=view,
        plot_payloads=plots, comparison_payload=None, cache_version="1", readout_cache_version="1")
    saved = json.loads(json.dumps(payload, allow_nan=False))
    assert saved["view"]["squeezing_view"] == "Quadrature phase scan"
    assert saved["view"]["quadrature_phase_scan"]["phase_deg"] == [0., 90., 180.]
    assert saved["view"]["quadrature_phase_scan"]["scope"] == "Effective Gaussian model"
    assert saved["view"]["gain_operating_point"] == {"tpd_mhz": -8., "probe_gain": 4.}
    assert saved["params"]["tpd"] == -8.
    trace = saved["plots"][0]["traces"][0]
    assert trace["x_unit"] == "deg" and trace["y_unit"] == "dB"
    assert export._curve_csv_bytes(trace).decode().splitlines()[0] == (
        "Analysis phase [deg],Quadrature noise [dB]")


def test_late_detuning_gesture_does_not_change_phase_operating_point(monkeypatch):
    state = {
        "fwm__tpd": -8., "fwm__squeezing_view": "Quadrature phase scan",
        "live-plot-fwm": {"gesture": json.dumps(dict(
            n=2, context="old-detuning-plot", values=dict(tpd=45.)))},
    }
    monkeypatch.setattr(liveplot.st, "session_state", state)
    liveplot._commit("live-plot-fwm", "fwm", "old-detuning-plot")
    assert state["fwm__tpd"] == -8.
    assert state["live-plot-fwm-ack"] == {"context": "old-detuning-plot", "n": 2}
