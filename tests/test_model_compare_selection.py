"""Saved comparison cursors cannot leak across datasets or point to a different group."""
import pytest

from gabes_ui.model_compare_plot import comparison_selection


@pytest.fixture
def comparison():
    return {
        "scheme": "sas",
        "axis": {"unit": "MHz", "reference": "Rb85 D1 line centroid"},
        "x": [-2., -1., 0., 1., 2.],
        "groups": [{"id": "left", "window_mhz": [-2., 0.]},
                   {"id": "right", "window_mhz": [1., 2.]}],
        "series": [{"power_mw": 0.}, {"power_mw": .5}, {"power_mw": 2.}],
        "default_power_mw": .5,
    }


def test_committed_sample_matches_physical_coordinate_and_context(comparison):
    saved = dict(group_id="right", sample_index=4, power_mw=2., show_powers=True, context="old")
    assert comparison_selection(comparison, saved, "old") == dict(
        group_id="right", sample_index=4, x_mhz=2., power_mw=.5, window_mhz=[1., 2.], show_powers=False)
    assert comparison_selection(comparison, saved, "new") == dict(
        group_id="left", sample_index=1, x_mhz=-1., power_mw=.5, window_mhz=[-2., 0.], show_powers=False)


def test_sample_outside_group_cannot_appear_in_export(comparison):
    saved = dict(group_id="right", sample_index=0, power_mw=1000., context="current")
    point = comparison_selection(comparison, saved, "current")
    assert point["x_mhz"] in (1., 2.)
    assert point["power_mw"] == .5


def test_generic_scan_selects_a_real_sample_without_groups(comparison):
    comparison["groups"] = []
    assert comparison_selection(comparison, None, "generic") == dict(
        group_id=None, sample_index=2, x_mhz=0., power_mw=.5, window_mhz=[-2., 2.], show_powers=False)


def test_custom_window_and_cursor_survive_export_without_a_preset(comparison):
    saved = dict(context="run", group_id="left", window_mhz=[-.2,1.2], sample_index=4)
    point = comparison_selection(comparison, saved, "run")
    assert point["group_id"] is None and point["window_mhz"] == [-.2,1.2]
    assert -.2 <= point["x_mhz"] <= 1.2
    saved["window_mhz"] = [0.1, 0.2]  # no sampled frequency
    assert comparison_selection(comparison, saved, "run")["group_id"] == "left"
    comparison["default_power_mw"] = 0.
    assert comparison_selection(comparison, saved, "run")["power_mw"] == 0.


def test_remount_restores_export_choice_after_streamlit_cleans_widget_state(comparison, monkeypatch):
    from gabes_ui import model_compare_plot as plot

    saved = dict(group_id="right", sample_index=4, power_mw=2., show_powers=True, context="old")
    session = {"plot": {"selection": saved}}
    monkeypatch.setattr(plot.st, "session_state", session)
    monkeypatch.setattr(plot, "_COMPONENT", lambda **kwargs: None)
    first = plot.render_model_compare(comparison, key="plot", context="old")
    del session["plot"]
    assert plot.render_model_compare(comparison, key="plot", context="old") == first
    assert plot.render_model_compare(comparison, key="plot", context="other") == first


def test_free_window_survives_new_power_but_late_gestures_and_new_axes_do_not(comparison,monkeypatch):
    from gabes_ui import model_compare_plot as plot

    session = {"plot": {"selection": dict(context="p1", group_id=None,
                                          window_mhz=[-.25, 1.25], sample_index=3)}}
    monkeypatch.setattr(plot.st, "session_state", session)
    monkeypatch.setattr(plot, "_COMPONENT", lambda **kwargs: None)
    sample = plot.render_model_compare(comparison,key="plot",context="p1")
    comparison["default_power_mw"] = 0.
    session["plot"]["selection"].update(window_mhz=[-2.,0.], sample_index=0)
    carried = plot.render_model_compare(comparison,key="plot",context="p0")
    assert carried == dict(sample,power_mw=0.)
    # Identical point count and endpoints are insufficient when a grid node changes.
    comparison["x"][1] = -1.1
    assert plot.render_model_compare(comparison,key="plot",context="grid2")["group_id"] == "left"
    session["plot"]["selection"].update(context="grid2",window_mhz=[1.,2.],sample_index=4)
    assert plot.render_model_compare(comparison,key="plot",context="grid2")["group_id"] == "right"
    comparison["axis"]["reference"] = "Rb85 D2 line centroid"
    assert plot.render_model_compare(comparison,key="plot",context="D2")["group_id"] == "left"


def test_native_group_command_wins_over_pending_gesture_and_free_drag_clears_group(comparison, monkeypatch):
    from gabes_ui import model_compare_plot as plot

    session = {"plot": {"selection": dict(context="run", group_id="left", sample_index=1)},
               "plot-group-control": "right"}
    monkeypatch.setattr(plot.st, "session_state", session)
    monkeypatch.setattr(plot, "_COMPONENT", lambda **kwargs: None)
    plot._choose_group(comparison, "plot", "run")
    picked = plot.render_model_compare(comparison, key="plot", context="run", show_groups=False)
    assert picked["group_id"] == "right" and picked["window_mhz"] == [1., 2.]
    # A gesture already on the wire before the native click cannot undo it.
    session["plot"]["selection"] = dict(context="run", group_id="left", sample_index=0,
                                        navigation_revision=0)
    plot._receive_selection(comparison, "plot", "run")
    assert session["plot-group-control"] == "right"
    # The acknowledged current browser gesture clears the native preset selection.
    session["plot"]["selection"] = dict(context="run", group_id=None, window_mhz=[-.5, 1.5],
                                        sample_index=3, navigation_revision=1)
    plot._receive_selection(comparison, "plot", "run")
    free = plot.render_model_compare(comparison, key="plot", context="run", show_groups=False)
    assert session["plot-group-control"] is None and free["group_id"] is None
    assert free["window_mhz"] == [-.5, 1.5] and free["x_mhz"] == 1.
