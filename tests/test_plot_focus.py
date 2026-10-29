"""Zoom changes the view, not the scientific data or cached figure state."""
from copy import deepcopy

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pytest

from gabes_ui.export import collect_view_plot_payloads
from gabes_ui.plot_focus import FocusWindow, focus_signature, focused_axes, resonance_window


@pytest.fixture
def eit_figure():
    x = np.linspace(-70., 70., 601)
    y = .014 / (1 + (x / .23) ** 2)
    fig, axes = plt.subplots(2, 1, sharex=True)
    axes[0].plot(x, y)
    axes[0].set_ylim(-.02, 1.02)
    axes[0].set_ylabel("Transmission")
    axes[1].plot(x, x / (1 + (x / .23) ** 2) * 1e-4)
    for axis in axes:
        axis.axvline(0.)
    yield fig
    plt.close(fig)


def test_eit_zoom_reveals_small_narrow_window_without_changing_samples(eit_figure):
    fig = eit_figure
    before = collect_view_plot_payloads({"figure": fig}, fig)
    limits = [(ax.get_xlim(), ax.get_ylim()) for ax in fig.axes]
    window = resonance_window(fig, "lambda", {"Dc": 0.})
    assert window is not None
    assert window.limits[1] - window.limits[0] < 4.
    with focused_axes(fig, window):
        assert fig.axes[0].get_ylim()[1] < .02
        assert fig.axes[0].get_xlim() == fig.axes[1].get_xlim()
        assert fig.axes[0].lines[0].get_marker() == "."
        assert collect_view_plot_payloads({"figure": fig}, fig) == before
    assert [(ax.get_xlim(), ax.get_ylim()) for ax in fig.axes] == limits
    assert fig.axes[0].lines[0].get_marker() == "None"
    assert collect_view_plot_payloads({"figure": fig}, fig) == before


def test_focus_restores_limits_and_markers_after_failed_renderer(eit_figure):
    fig = eit_figure
    before = [(ax.get_xlim(), ax.get_ylim(), ax.get_autoscalex_on(), ax.get_autoscaley_on())
              for ax in fig.axes]
    with pytest.raises(RuntimeError, match="render failed"):
        with focused_axes(fig, resonance_window(fig, "lambda")):
            raise RuntimeError("render failed")
    assert [(ax.get_xlim(), ax.get_ylim(), ax.get_autoscalex_on(), ax.get_autoscaley_on())
            for ax in fig.axes] == before
    assert all(ax.lines[0].get_marker() == "None" for ax in fig.axes)


def test_gap_at_resonance_is_unavailable_not_a_peak(eit_figure):
    fig = eit_figure
    y = np.array(fig.axes[0].lines[0].get_ydata())
    y[y.size // 2] = np.nan
    fig.axes[0].lines[0].set_ydata(y)
    assert resonance_window(fig, "lambda") is None


def test_fwm_focus_signature_ignores_cursor_but_changes_with_gain():
    x = np.linspace(-500., 500., 501)
    fig, axes = plt.subplots(2, 1, sharex=True)
    try:
        for ax in axes:
            ax.plot(x, np.exp(-((x - 30) / 20) ** 2) + .001)
            ax.axvline(-8.)
        axes[0].set_yscale("log")
        axes[0].set_xlim(-500, 500)
        signature = focus_signature(fig)
        first = resonance_window(fig, "fwm")
        axes[0].lines[1].set_xdata([-20., -20.])
        assert focus_signature(fig) == signature
        assert resonance_window(fig, "fwm") != first
        axes[0].lines[0].set_ydata(axes[0].lines[0].get_ydata() * 2)
        assert focus_signature(fig) != signature
        with focused_axes(fig, first):
            assert axes[0].get_yscale() == "log"
            assert np.isfinite(axes[0].get_ylim()).all()
            assert axes[0].get_ylim()[0] > 0
        assert resonance_window(fig, "fwm", {"kind": "biphoton"}) is None
    finally:
        plt.close(fig)


def test_sas_subdoppler_zoom_keeps_original_resolution_qualification():
    x = np.linspace(-5., 5., 1001)
    y = 1 - .5 * np.exp(-(x / .6) ** 2) + .1 * np.exp(-((x - .1) / .025) ** 2)
    fig, ax = plt.subplots()
    try:
        ax.plot(x, y)
        metrics = [{"label": "SAS resolution", "value": "resolution-limited"}]
        view = {"metrics": metrics}
        before = deepcopy(view)
        window = resonance_window(fig, "sas", view=view)
        assert window is not None
        assert window.limits[0] < .1 < window.limits[1]
        assert window.limits[1] - window.limits[0] < 1
        with focused_axes(fig, window):
            np.testing.assert_array_equal(ax.lines[0].get_xdata(), x)
            np.testing.assert_array_equal(ax.lines[0].get_ydata(), y)
        assert view == before
    finally:
        plt.close(fig)


def test_inverted_axes_and_unavailable_gaps_survive_zoom():
    fig, ax = plt.subplots()
    x = np.linspace(-2., 2., 11)
    y = np.array([0., 1., 2., np.nan, 3., 4., 5., 6., 7., 8., 9.])
    try:
        ax.plot(x, y)
        ax.invert_xaxis()
        ax.invert_yaxis()
        before = (ax.get_xlim(), ax.get_ylim())
        with focused_axes(fig, FocusWindow((-1., 1.), "test")):
            assert ax.get_xlim() == (1., -1.)
            assert ax.get_ylim()[0] > ax.get_ylim()[1]
            np.testing.assert_array_equal(ax.lines[0].get_ydata(), y)
        assert (ax.get_xlim(), ax.get_ylim()) == before
    finally:
        plt.close(fig)


def test_focus_toggle_preserves_full_export_and_figure_limits():
    from streamlit.testing.v1 import AppTest
    app = AppTest.from_string('''
import matplotlib.pyplot as plt
import numpy as np
import streamlit as st
from types import SimpleNamespace
from gabes_ui.plotcard import render_plot_card
x = np.linspace(-70, 70, 601)
fig, ax = plt.subplots()
ax.plot(x, .014 / (1 + (x / .23) ** 2))
ax.set_ylim(-.02, 1.02)
ax.set_title('Test resonance')
view = dict(figure=fig, metrics=[dict(label='Window FWHM', value='0.46 MHz')])
st.session_state['before_limits'] = [ax.get_xlim(), ax.get_ylim()]
_, plots, _ = render_plot_card(view, SimpleNamespace(name='lambda'), {}, dict(Dc=0.))
st.session_state['after_limits'] = [ax.get_xlim(), ax.get_ylim()]
st.session_state['plots'] = plots
''').run(timeout=20)
    assert not app.exception
    before = deepcopy(app.session_state["plots"])
    caption_count = len(app.caption)
    toggle = next(item for item in app.toggle if item.label == "Resonance zoom")
    toggle.set_value(True).run(timeout=20)
    assert not app.exception
    assert app.session_state["plots"] == before
    assert app.session_state["before_limits"] == app.session_state["after_limits"]
    assert len(app.caption) == caption_count
    assert any("Original samples" in caption.value for caption in app.caption)
    toggle.set_value(False).run(timeout=20)
    assert not app.exception
    assert app.session_state["plots"] == before
    assert len(app.caption) == caption_count


def test_fwm_zoom_stays_put_on_cursor_updates_and_resets_when_unavailable():
    from streamlit.testing.v1 import AppTest
    app = AppTest.from_string('''
import matplotlib.pyplot as plt
import numpy as np
import streamlit as st
from types import SimpleNamespace
from gabes_ui.plotcard import render_plot_card
tpd = st.slider('Selected point', -100., 100., -8.)
available = st.checkbox('Available curve', value=True)
x = np.unique(np.r_[np.linspace(-500., 500., 501), tpd])
fig, axes = plt.subplots(2, 1, sharex=True)
for ax in axes:
    ax.plot(x, np.exp(-((x - 30) / 20) ** 2) if available else np.full_like(x, np.nan))
    ax.axvline(tpd)
    ax.set_xlim(-500, 500)
render_plot_card(dict(figure=fig), SimpleNamespace(name='fwm'), {},
                 params=dict(tpd=tpd, gain_calibration=True))
''').run(timeout=20)
    assert not app.exception
    next(t for t in app.toggle if t.label == "Resonance zoom").set_value(True).run(timeout=20)
    first = app.session_state["fwm___resonance_zoom_Plot_window"][1]
    app.slider[0].set_value(-20.).run(timeout=20)
    assert not app.exception
    assert app.session_state["fwm___resonance_zoom_Plot_window"][1] == first
    app.checkbox[0].set_value(False).run(timeout=20)
    assert not app.exception
    toggle = next(t for t in app.toggle if t.label == "Resonance zoom")
    assert not toggle.value and toggle.disabled
    assert any("showing the full range" in caption.value for caption in app.caption)
