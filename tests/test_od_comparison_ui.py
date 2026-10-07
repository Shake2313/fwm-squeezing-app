"""OD comparison widgets keep plotted curves and exported alignment consistent."""
from types import SimpleNamespace

import numpy as np
import pytest


@pytest.fixture
def comparison_app(monkeypatch):
    from streamlit.delta_generator import DeltaGenerator
    from streamlit.testing.v1 import AppTest
    from gabes_ui import export as ui_export

    captures = []

    def uploaded_csv(_self, *args, **kwargs):
        import streamlit as st
        return SimpleNamespace(
            name="seconds_scope.csv",
            getvalue=lambda: st.session_state["_comparison_source_bytes"],
        )

    def capture_export(_container, **kwargs):
        captures.append(kwargs)

    monkeypatch.setattr(DeltaGenerator, "file_uploader", uploaded_csv)
    monkeypatch.setattr(ui_export, "render_export", capture_export)
    app = AppTest.from_string("""
import numpy as np
import streamlit as st
from gabes import schemes
from gabes_ui.plotcard import render_plot_section

scheme = schemes.get('sas')
if '_comparison_base_view' not in st.session_state:
    params = scheme.defaults()
    params.update(pump_power_mw=0.0, line='D2', temp_c=25.0, scan_points=401)
    raw = scheme.compute(params)
    view = scheme.observables(raw, params)
    # AppTest's segmented-control serializer currently assumes multiselect.
    # Exercise the real default Transmission figure as a single plot page.
    view = dict(view, figure_views=[])
    theory = view['figure'].axes[0].lines[0]
    x = np.asarray(theory.get_xdata())
    y = np.asarray(theory.get_ydata())
    position = np.linspace(0.0, 1.0, 601) ** 1.7
    sample_x = x[0] + (x[-1] - x[0]) * position
    seconds = 1.0e-6 + 2.0e-7 * position
    transmission = np.interp(sample_x, x, y)
    signal = 0.23 + 0.68 * transmission
    st.session_state['_comparison_source_bytes'] = (
        '\\n'.join(f'{a:.15g},{b:.15g}' for a, b in zip(seconds, signal))
    ).encode()
    st.session_state['_comparison_source_transmission'] = transmission
    st.session_state['_comparison_source_seconds'] = seconds
    st.session_state['_comparison_theory_xlim'] = view['figure'].axes[0].get_xlim()
    st.session_state['_comparison_original_line_count'] = len(view['figure'].axes[0].lines)
    st.session_state['_comparison_base_view'] = view
    st.session_state['_comparison_raw'] = raw
    st.session_state['_comparison_params'] = params
    st.session_state['sas___overlay'] = True

render_plot_section(
    st.session_state['_comparison_base_view'], scheme,
    {spec.name: spec for spec in scheme.param_schema()},
    st.session_state['_comparison_raw'], st.session_state['_comparison_params'],
    scheme.cache_version, 'comparison-widget-test',
)
""", default_timeout=30).run()
    assert not app.exception
    assert captures and captures[-1]["comparison_payload"]["enabled"]
    yield app, captures
    import matplotlib.pyplot as plt
    view = app.session_state["_comparison_base_view"]
    plt.close(view["figure"])


def _payload(captures):
    return captures[-1]["comparison_payload"]


def _csv_plot_trace(captures):
    traces = captures[-1]["plot_payloads"][0]["traces"]
    candidates = [trace for trace in traces
                  if trace["figure_curve_label"].startswith("Experimental CSV")]
    assert len(candidates) == 1
    return candidates[0]


def _assert_export_matches_plot(captures):
    plotted = _csv_plot_trace(captures)
    exported = _payload(captures)["comparison_trace"]
    np.testing.assert_array_equal(plotted["x_data"], exported["detuning"])
    np.testing.assert_array_equal(plotted["y_data"], exported["transmission"])


def _apply(app):
    submit = next(button for button in app.button if button.label == "Apply alignment")
    submit.click().run()
    assert not app.exception


def test_seconds_upload_frames_nonuniform_data_and_repeated_apply_keeps_clean_theory(comparison_app):
    app, captures = comparison_app
    payload = _payload(captures)
    xlim = np.asarray(app.session_state["_comparison_theory_xlim"])
    aligned = np.asarray(payload["comparison_trace"]["detuning"])
    assert payload["in_view"]
    assert payload["alignment"]["scale"] > 1e7
    assert aligned.min() > xlim.min() and aligned.max() < xlim.max()
    displayed_ylim = captures[-1]["view"]["figure"].axes[0].get_ylim()
    assert displayed_ylim[0] < 0.0 and displayed_ylim[1] > 1.0
    assert 0.5 * (aligned.min() + aligned.max()) == pytest.approx(xlim.mean())
    # Nonuniform sampling has a different median from its endpoint midpoint.
    assert not np.isclose(np.median(aligned), xlim.mean(), atol=0.01)
    _assert_export_matches_plot(captures)

    fine = app.slider(key="sas___csv_x_multiplier")
    assert fine.min == 0.1 and fine.max == 3.0
    assert fine.proto.form_id
    base_scale = payload["alignment"]["scale"]
    fine.set_value(1.15)
    _apply(app)
    assert _payload(captures)["alignment"]["scale"] == pytest.approx(1.15 * base_scale)
    _assert_export_matches_plot(captures)
    first_applied = _payload(captures)["comparison_trace"]

    _apply(app)
    _assert_export_matches_plot(captures)
    assert _payload(captures)["comparison_trace"] == first_applied
    expected_count = app.session_state["_comparison_original_line_count"] + 1
    assert len(captures[-1]["plot_payloads"][0]["traces"]) == expected_count
    # The persistent cached figure never inherits any experimental curve.
    assert len(app.session_state["_comparison_base_view"]["figure"].axes[0].lines) == expected_count - 1


def test_fit_alignment_improves_seconds_overlay_without_changing_physics(comparison_app):
    app, captures = comparison_app
    app.slider(key="sas___csv_x_multiplier").set_value(0.8)
    app.slider(key="sas___csv_x_shift").set_value(0.4)
    _apply(app)
    original_params = dict(captures[-1]["params"])
    app.button(key="sas___csv_fit").click().run()
    assert not app.exception
    fit = _payload(captures)["alignment_fit"]
    assert fit is not None and fit["success"]
    assert fit["rmse"] < 0.5 * fit["initial_rmse"]
    assert fit["overlap_fraction"] > 0.8
    assert fit["physical_parameters_fitted"] is False
    assert captures[-1]["params"] == original_params
    _assert_export_matches_plot(captures)


def test_axis_reframing_preserves_signal_adjustments(comparison_app):
    app, captures = comparison_app
    app.number_input(key="sas___csv_y_scale").set_value(0.8)
    app.number_input(key="sas___csv_y_offset").set_value(0.1)
    app.checkbox(key="sas___csv_invert").check()
    app.checkbox(key="sas___csv_reverse").check()
    _apply(app)
    transmission = _payload(captures)["comparison_trace"]["transmission"]

    app.button(key="sas___csv_frame").click().run()
    assert not app.exception
    assert _payload(captures)["alignment"]["reverse"]
    assert _payload(captures)["comparison_trace"]["transmission"] == transmission
    _assert_export_matches_plot(captures)

    app.button(key="sas___csv_reset").click().run()
    assert not app.exception
    alignment = _payload(captures)["alignment"]
    assert alignment["scale"] == 1.0 and alignment["shift"] == 0.0
    assert not alignment["reverse"]
    assert _payload(captures)["comparison_trace"]["transmission"] == transmission
    _assert_export_matches_plot(captures)


def test_relative_inversion_and_absolute_calibration_export_the_plotted_transmission(comparison_app):
    app, captures = comparison_app
    relative = np.asarray(_payload(captures)["comparison_trace"]["transmission"])
    app.checkbox(key="sas___csv_invert").check()
    app.number_input(key="sas___csv_y_scale").set_value(0.8)
    app.number_input(key="sas___csv_y_offset").set_value(0.1)
    _apply(app)
    inverted = np.asarray(_payload(captures)["comparison_trace"]["transmission"])
    np.testing.assert_allclose(inverted, 0.8 * (1.0 - relative) + 0.1)
    _assert_export_matches_plot(captures)

    app.selectbox(key="sas___csv_calibration_mode").select(
        "Absolute transmission (gain/offset)").run()
    assert not app.exception
    app.number_input(key="sas___csv_detector_gain").set_value(0.68)
    app.number_input(key="sas___csv_detector_offset").set_value(0.23)
    app.checkbox(key="sas___csv_auto_correct").uncheck().run()
    assert not app.exception
    payload = _payload(captures)
    assert payload["alignment"]["invert"] is False
    assert payload["alignment"]["signal_scale"] == 1.0
    assert payload["alignment"]["signal_offset"] == 0.0
    assert app.checkbox(key="sas___csv_invert").disabled
    assert app.number_input(key="sas___csv_y_scale").disabled
    assert app.number_input(key="sas___csv_y_offset").disabled
    np.testing.assert_allclose(payload["comparison_trace"]["transmission"],
                               app.session_state["_comparison_source_transmission"], atol=1e-12)
    _assert_export_matches_plot(captures)

    # Running the X alignment fitter cannot rewrite measured absolute Y values.
    absolute = payload["comparison_trace"]["transmission"]
    app.button(key="sas___csv_fit").click().run()
    assert not app.exception
    np.testing.assert_array_equal(_payload(captures)["comparison_trace"]["transmission"], absolute)
    _assert_export_matches_plot(captures)


@pytest.mark.parametrize("widget_key", ["sas___csv_x_scale", "sas___csv_y_scale"])
def test_zero_alignment_scale_apply_and_fit_show_readable_errors(comparison_app, widget_key):
    app, captures = comparison_app
    app.number_input(key=widget_key).set_value(0.0)
    _apply(app)
    assert app.error
    assert any("scale" in error.value.lower() and "positive" in error.value.lower()
               for error in app.error)
    assert captures[-1]["comparison_payload"] is None

    app.button(key="sas___csv_fit").click().run()
    assert not app.exception
    assert app.error
    assert any("scale" in error.value.lower() and "positive" in error.value.lower()
               for error in app.error)
    assert captures[-1]["comparison_payload"] is None
    # Invalid transforms leave a clean theoretical curve available to inspect.
    assert len(captures[-1]["plot_payloads"][0]["traces"]) == (
        app.session_state["_comparison_original_line_count"])


@pytest.mark.parametrize("overflow_axis", ["x", "y"])
def test_overflowing_alignment_is_reported_without_crashing_the_plot(comparison_app, overflow_axis):
    app, captures = comparison_app
    if overflow_axis == "x":
        app.number_input(key="sas___csv_x_scale").set_value(1e308)
        app.slider(key="sas___csv_x_multiplier").set_value(3.0)
    else:
        app.number_input(key="sas___csv_y_scale").set_value(1e308)
        app.number_input(key="sas___csv_y_offset").set_value(1e308)
    _apply(app)
    assert app.error
    assert any("finite" in error.value.lower() for error in app.error)
    assert captures[-1]["comparison_payload"] is None

    app.button(key="sas___csv_fit").click().run()
    assert not app.exception
    payload = captures[-1]["comparison_payload"]
    if payload is None:
        assert app.error and any("finite" in error.value.lower() for error in app.error)
    else:
        # Relative Y fitting may replace the unusable starting gain/offset.
        assert payload["alignment_fit"]["success"]
        assert np.all(np.isfinite(payload["comparison_trace"]["transmission"]))
        _assert_export_matches_plot(captures)
