"""Preview math, committed exports, and unavailable intervals stay consistent."""
import base64
import json
from pathlib import Path
import shutil
import subprocess
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pytest

from gabes.experimental_alignment import AxisAlignment
from gabes_ui import liveplot


def test_browser_affine_and_nf_interpolation_match_python():
    node = shutil.which("node")
    if not node:
        pytest.skip("Node.js is needed to exercise browser preview math")
    asset = Path(liveplot.__file__).resolve().parents[1] / "assets/ui/liveplot.js"
    url = "data:text/javascript;base64," + base64.b64encode(asset.read_bytes()).decode()
    source = np.array([1e-6, 1.01e-6, 1.07e-6, 1.2e-6])
    y = np.array([0.2, 0.7, 0.3, 1.1])
    cases, expected = [], []
    for reverse in (False, True):
        for invert in (False, True):
            old = AxisAlignment(5e7, 0.4, reverse, 0.8, -0.2)
            new = AxisAlignment(6e7, -0.1, reverse, 1.3, 0.15)
            old_x = old.transform_x(source)
            old_y = old.transform_y(y, invert=invert)
            cases.append(dict(x=old_x.tolist(), y=old_y.tolist(),
                              pivot=(source.min() + source.max()) / 2,
                              original=dict(shift=old.shift, multiplier=1.,
                                            signal_scale=old.signal_scale, signal_offset=old.signal_offset),
                              values=dict(shift=new.shift, multiplier=1.2,
                                          signal_scale=new.signal_scale, signal_offset=new.signal_offset)))
            expected.append(np.column_stack([new.transform_x(source), new.transform_y(y, invert=invert)]))
    script = f"""
      import {{alignPoint,interpolateCurve}} from {json.dumps(url)};
      const cases={json.dumps(cases)};
      const affine=cases.map(c=>c.x.map((x,i)=>alignPoint(x,c.y[i],c.pivot,c.original,c.values)));
      const c={{x:[-1,1],y:[0,10*Math.log10(9)]}};
      const gap={{x:[-1,0,1],y:[0,null,10]}};
      console.log(JSON.stringify({{affine, nf:interpolateCurve(c,0,true),
        legacy:interpolateCurve(c,0), gap:interpolateCurve(gap,0.5,true),
        outside:interpolateCurve(c,2,true), endpoint:interpolateCurve(gap,-1,true)}}));
    """
    result = subprocess.run([node, "--input-type=module"], input=script, text=True,
                            capture_output=True, check=True, timeout=15)
    output = json.loads(result.stdout)
    np.testing.assert_allclose(output["affine"], expected, rtol=1e-12, atol=1e-12)
    assert output["nf"] == pytest.approx(10 * np.log10(5))
    assert output["legacy"] == pytest.approx(5 * np.log10(9))
    assert output["gap"] is None and output["outside"] is None
    assert output["endpoint"] == 0


@pytest.fixture
def seeded_figure():
    fig, axes = plt.subplots(2, 1, figsize=(8.5, 6.4))
    for axis in axes:
        axis.plot([-1., 0., 1.], [1., 2., 4.])
        axis.axvline(0.2, color="crimson")
        axis.axhline(0)
        axis.scatter([0.2], [2.4], color="crimson")
        axis.set_xlim(-1, 1)
    axes[0].set_yscale("log")
    yield fig
    plt.close(fig)


def test_live_render_restores_all_export_artists_and_uses_fixed_canvas(monkeypatch, seeded_figure):
    captured = {}
    monkeypatch.setattr(liveplot, "_MOUNTABLE", True)
    monkeypatch.setattr(liveplot, "_COMPONENT", lambda **kw: captured.update(kw))
    before = [(line.get_visible(), np.array(line.get_ydata()))
              for axis in seeded_figure.axes for line in axis.lines]
    assert liveplot.render_live_plot(seeded_figure, "fwm", canonical_noise=True)
    data = captured["data"]
    assert data["linear_noise"] is True
    assert data["width"] == 1020 and data["height"] == 768
    assert data["axes"][0]["log"] is True
    assert data["image"].startswith("data:image/png;base64,iVBOR")
    for (visible, values), line in zip(before, [line for axis in seeded_figure.axes for line in axis.lines]):
        assert line.get_visible() == visible
        np.testing.assert_array_equal(line.get_ydata(), values)
    assert all(axis.collections[0].get_visible() for axis in seeded_figure.axes)


def test_snapshot_error_restores_visibility(monkeypatch, seeded_figure):
    monkeypatch.setattr(liveplot, "_MOUNTABLE", True)
    def failed(*args, **kwargs):
        raise OSError("snapshot failed")
    monkeypatch.setattr(seeded_figure, "savefig", failed)
    with pytest.raises(OSError, match="snapshot failed"):
        liveplot.render_live_plot(seeded_figure, "fwm")
    assert all(axis.lines[1].get_visible() and axis.collections[0].get_visible()
               for axis in seeded_figure.axes)


def test_absolute_transmission_preview_has_no_y_adjustments():
    fig, axis = plt.subplots()
    axis.plot([0, 1, 2], [1, .2, 1], label="Experimental CSV — absolute transmission")
    comparison = dict(alignment=dict(shift=0., scale=1., manual_multiplier=1.,
                                    signal_scale=1., signal_offset=0.),
                      descriptor=dict(axis_index=0, axis_x_unit="GHz"),
                      calibration_mode="absolute_dark_reference")
    try:
        data, moving = liveplot.preview_data(fig, "sas", comparison)
        assert len(moving) == 1
        assert [control["name"] for control in data["controls"]] == ["shift", "multiplier"]
    finally:
        plt.close(fig)


def test_commit_updates_alignment_before_widgets_and_invalidates_fit(monkeypatch):
    state = {"live-plot-sas": {"gesture": json.dumps(dict(context="file-hash", values=dict(
        shift=0.4, multiplier=1.2, signal_scale=0.8, signal_offset=0.1)))},
        "sas___csv_alignment_fit": "old fit"}
    monkeypatch.setattr(liveplot.st, "session_state", state)
    liveplot._commit("live-plot-sas", "sas", "file-hash")
    assert state["sas___csv_x_shift"] == 0.4
    assert state["sas___csv_x_multiplier"] == 1.2
    assert state["sas___csv_y_scale"] == 0.8
    assert state["sas___csv_y_offset"] == 0.1
    assert "sas___csv_alignment_fit" not in state


@pytest.mark.parametrize("payload", [dict(context="stale-file", values=dict(shift=1)),
                                      dict(context="current", values=dict(multiplier=0)),
                                      dict(context="current", values=dict(signal_scale=-1)),
                                      dict(context="current", values=dict(shift=float("inf")))])
def test_invalid_or_stale_preview_cannot_commit(monkeypatch, payload):
    state = {"live-plot-sas": {"gesture": json.dumps(payload)}}
    monkeypatch.setattr(liveplot.st, "session_state", state)
    liveplot._commit("live-plot-sas", "sas", "current")
    assert len(state) == 1


def test_sas_shows_only_transmission_but_keeps_other_export_curves():
    from streamlit.testing.v1 import AppTest
    app = AppTest.from_string('''
import streamlit as st
from gabes import schemes
from gabes_ui.plotcard import render_plot_card
scheme = schemes.get('sas')
params = dict(scheme.defaults(), pump_power_mw=0., scan_points=401)
raw = scheme.compute(params)
view = scheme.observables(raw, params)
st.session_state['expected_export_plot_labels'] = ['primary'] + [
    page['label'] for page in view.get('figure_views', [])
    if page['figure'] is not view['figure']]
_, plots, _ = render_plot_card(view, scheme, {}, raw)
st.session_state['export_plot_labels'] = [plot['label'] for plot in plots]
''', default_timeout=30).run()
    assert not app.exception
    assert not app.get("button_group")  # No plot-view selector when there is one view.
    assert app.session_state["export_plot_labels"] == app.session_state["expected_export_plot_labels"]
    assert len(app.session_state["export_plot_labels"]) > 1


def test_csv_cache_uses_current_class_after_streamlit_module_replacement():
    script = '''
import importlib, pickle, sys
from gabes_ui import plotcard
old = sys.modules.pop('gabes.experimental_csv')
current = importlib.import_module('gabes.experimental_csv')
assert old.ExperimentalCSV is not current.ExperimentalCSV
trace = plotcard._cached_experimental_csv(b'0,1\\n1,0.8\\n2,0.2\\n3,0.8\\n4,1\\n', False)
assert type(trace) is current.ExperimentalCSV
pickle.dumps(trace)
'''
    result = subprocess.run([sys.executable, "-c", script], text=True, encoding="utf-8",
                            errors="replace", capture_output=True, timeout=30)
    assert result.returncode == 0, result.stdout + result.stderr


def test_live_component_real_runtime_registry_error_remains_visible(monkeypatch, seeded_figure):
    import importlib
    from types import SimpleNamespace
    from streamlit.runtime import Runtime
    bidi = importlib.import_module("streamlit.components.v2.bidi_component.main")
    monkeypatch.setattr(liveplot, "_MOUNTABLE", True)
    monkeypatch.setattr(bidi, "get_script_run_ctx", lambda: SimpleNamespace(_enqueue=lambda message: None))
    monkeypatch.setattr(Runtime, "instance", lambda: SimpleNamespace(bidi_component_registry={}))
    with pytest.raises(ValueError, match="Component 'gabes_live_plot' is not registered"):
        liveplot.render_live_plot(seeded_figure, "fwm")


@pytest.mark.parametrize("failure_source", ["component", "callback"])
def test_live_component_value_errors_remain_visible_in_apptest(monkeypatch, failure_source):
    from streamlit.testing.v1 import AppTest
    message = "Component 'gabes_live_plot' is not registered"
    def component(**kwargs):
        if failure_source == "component":
            raise ValueError(message)
        kwargs["on_gesture_change"]()
    def fail(*args):
        raise ValueError(message)
    monkeypatch.setattr(liveplot, "_MOUNTABLE", True)
    monkeypatch.setattr(liveplot, "_COMPONENT", component)
    monkeypatch.setattr(liveplot, "_commit", fail)
    app = AppTest.from_string('''
import matplotlib.pyplot as plt
from gabes_ui import liveplot
fig, axes = plt.subplots(2, 1)
for axis in axes:
    axis.plot([-1., 0., 1.], [1., 2., 4.])
    axis.axvline(0.2)
    axis.axhline(0.)
    axis.scatter([0.2], [2.4])
liveplot.render_live_plot(fig, "fwm")
plt.close(fig)
''').run(timeout=10)
    assert len(app.exception) == 1
    assert app.exception[0].message == message
