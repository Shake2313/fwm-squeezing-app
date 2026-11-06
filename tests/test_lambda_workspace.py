"""Stage-2 Lambda views keep native scientific arrays and linked control identity."""
import json
from copy import deepcopy

import numpy as np
import pytest

from gabes import schemes
from gabes.core import blas_single_thread
from gabes.schemes.lambda_comparison import _conditions, _readouts, comparison_key
from gabes_ui.spectral_compare import workspace_selection, selected_parameter_change, comparison_plots
from gabes_ui.spectral_sweep import resolve_row


def native_dataset(mode="EIT", comparison="Coupling power", **changes):
    scheme = schemes.get("lambda")
    params = {**scheme.defaults(), **scheme.recommended_defaults({})[mode], "doppler": "off",
              "lambda_compare": comparison, **changes}
    key = comparison_key(params)
    source = {}
    def resolve(p):
        raw = scheme.compute(p)
        source[p[key]] = scheme.headless_observables(raw, p)
        return raw
    with blas_single_thread():
        raw = resolve(params)
        data = scheme.model_comparison(params, resolve_raw=resolve, current_raw=raw, current_view=source[params[key]])
    return scheme, params, data, source


@pytest.mark.parametrize("mode", ["EIT", "AT", "CPT"])
@pytest.mark.parametrize("comparison", ["Coupling power", "Ground relaxation", "Beam angle", "Detuning map"])
def test_all_views_use_native_grids_transmission_depth_and_dispersion(mode, comparison):
    scheme, params, data, source = native_dataset(mode, comparison, coupling_detuning_mhz=2.3)
    key = comparison_key(params)
    assert key == data["varied"]["key"]
    assert params[key] in data["varied"]["values"]
    assert "lambda_compare" not in scheme.recompute_keys()
    for row in data["series"]:
        native = source[row["parameter"]["value"]]["lambda_readout"]
        x = native["x_mhz"]
        if mode == "CPT" and comparison != "Detuning map":
            x = (x-native["center_mhz"])*1000
        np.testing.assert_array_equal(row["x"], x)
        np.testing.assert_array_equal(row["observables"]["transmission"], native["transmission"])
        np.testing.assert_array_equal(row["observables"]["dispersion"], native["dispersion"])
        np.testing.assert_array_equal(row["observables"]["optical_depth"], native["absorption_per_m"]*params["cell_mm"]*1e-3)
        assert selected_parameter_change(scheme, data, {"series_id": row["id"]}, params) in ({}, {key: row["parameter"]["value"]})
    if comparison == "Detuning map":
        assert data["axis"]["unit"] == "MHz"
        assert data["map_guides"][0]["x_per_y"] == 1
        assert data["display"]["window"][0] < min(data["varied"]["values"])
        assert data["display"]["window"][1] > max(data["varied"]["values"])
        for row in data["series"]:
            assert row["x"][300] == pytest.approx(row["parameter"]["value"])
    json.dumps(data, allow_nan=False)


def test_at_peak_positions_are_the_existing_extractor_samples_and_spacing():
    _, _, data, source = native_dataset("AT", coupling_rabi_mhz=15.)
    assert data["display"]["condition_axis"] == "coupling_rabi"
    for row in data["series"]:
        peaks = [m["x"] for m in row["markers"] if m["label"].startswith("Absorption peak")]
        metric = next(m for m in source[row["parameter"]["value"]]["metrics"] if m["label"] == "AT splitting")
        assert peaks == metric.get("sample_positions", [])
        if peaks:
            assert all(p in row["x"] for p in peaks)
            assert peaks[1]-peaks[0] == metric["numeric_value"]
    _, _, weak, _ = native_dataset("AT", coupling_rabi_mhz=.1)
    current = next(row for row in weak["series"] if row["parameter"]["value"] == 1.)
    assert not any(m["label"].startswith("Absorption peak") for m in current["markers"])


def test_cpt_contrast_width_and_quality_share_the_actual_control_axis():
    _, _, data, source = native_dataset("CPT", coupling_rabi_mhz=.8)
    assert data["display"]["default_mode"] == "curves"
    assert data["summary"]["axis"]["unit"] == "mW"
    assert [p["unit"] for p in data["summary"]["panels"]] == ["1", "kHz", "kHz⁻¹"]
    for i, row in enumerate(data["series"]):
        values = {r["key"]: r["value"] for r in row["readouts"]}
        t = source[row["parameter"]["value"]]["lambda_readout"]["transmission"]
        assert values["contrast"] == float(t[300]-np.min(t))
        if values["window_width"] is not None:
            assert values["contrast_per_width"] == values["contrast"]/values["window_width"]
        for series in data["summary"]["series"]:
            assert series["condition_ids"][i] == row["id"]
            assert series["x"][i] == row["parameter"]["value"]


def test_clipped_width_is_withheld_without_mutating_original_metric():
    view = {"metrics": [{"label": "Transmission at resonance", "numeric_value": 1.},
                        {"label": "Window FWHM", "numeric_value": 2., "unit": "kHz"}]}
    readout = dict(coupling_rabi_mhz=1., x_mhz=np.arange(-2,3)/1000,
                   transmission=np.array([.2,.7,1.,.9,.8]), center_mhz=0.)
    before = deepcopy(view)
    rows = {r["key"]: r for r in _readouts(view, readout, "cpt")}
    assert rows["window_width"]["value"] is None
    assert rows["window_width"]["status"] == "scan-edge-limited"
    assert rows["contrast_per_width"]["value"] is None
    assert view == before
    readout["transmission"] = np.array([.2, np.nan, 1., .7, .2])
    gap = {r["key"]: r for r in _readouts(view, readout, "cpt")}
    assert gap["window_width"]["value"] is None


def test_default_modes_and_rabi_coordinate_validation():
    for mode in ("EIT", "AT", "CPT"):
        _, params, data, _ = native_dataset(mode)
        selection = workspace_selection(data, None, "ctx", params)
        assert selection["display_mode"] == ("curves" if mode == "CPT" else "map")
        assert selection["condition_axis"] == ("coupling_rabi" if mode == "AT" else "parameter")
        assert selection["window"][0] <= selection["x"] <= selection["window"][1]
        factor = data["condition_axes"][1]["transform"]["factor"]
        for row in data["series"]:
            assert factor*np.sqrt(row["parameter"]["value"]) == pytest.approx(row["coordinates"]["coupling_rabi"], rel=1e-14)
        forged = workspace_selection(data, {**selection, "condition_axis": "invented"}, "ctx", params)
        assert forged["condition_axis"] == "parameter"


def test_generated_condition_grids_preserve_one_exact_nonround_anchor():
    scheme=schemes.get("lambda")
    for comparison, key, value in [("Coupling power","coupling_power_mw",1.23),
                                   ("Detuning map","coupling_detuning_mhz",2.3)]:
        params={**scheme.defaults(), "lambda_compare":comparison, key:value}
        _, values=_conditions(scheme,params,key)
        assert len(values)==9 and values.count(value)==1


def test_new_comparison_axes_use_the_shared_incremental_row_adapter():
    scheme = schemes.get("lambda")
    for comparison, value in [("Ground relaxation", 2.5), ("Beam angle", 3.), ("Detuning map", 4.2)]:
        params = {**scheme.defaults(), "doppler": "off", "lambda_compare": comparison}
        key = comparison_key(params)
        calls = []
        def resolve(p):
            calls.append(p[key])
            return scheme.compute(p)
        with blas_single_thread():
            row = resolve_row(scheme, params, key, value, resolve)
        assert calls == [value]
        assert row["series"]["parameter"]["value"] == value
        assert len(row["series"]["x"]) == 601


def test_depth_export_stays_finite_when_transmission_underflows():
    _, _, data, _ = native_dataset("EIT", cell_mm=200., temp_c=120.)
    row = data["series"][0]
    assert 0. in row["observables"]["transmission"]
    assert all(np.isfinite(row["observables"]["optical_depth"]))
    payload = next(p for p in comparison_plots(data) if p["label"] == "Optical depth")
    assert payload["traces"][0]["y_data"] == row["observables"]["optical_depth"]


def test_comparison_axis_changes_keep_one_editor_and_canonical_value():
    from streamlit.testing.v1 import AppTest
    app = AppTest.from_string('''
import streamlit as st
from gabes import schemes
from gabes_ui.conditions import render_conditions, render_comparison_controls, render_power_inputs
scheme=schemes.get("lambda")
for name,value in scheme.defaults().items():
    st.session_state.setdefault("lambda__"+name,value)
render_conditions(scheme,scheme.param_schema())
render_comparison_controls(scheme,scheme.param_schema())
render_power_inputs(scheme,scheme.param_schema())
''').run()
    for choice, label, canonical, value in [
        ("Ground relaxation", "Buffer ground relaxation (kHz)", "buffer_ground_relax_khz", "2.5"),
        ("Beam angle", "Beam angle mismatch (mrad)", "beam_angle_mrad", "3.2"),
        ("Detuning map", "Coupling detuning (MHz)", "coupling_detuning_mhz", "-4.1"),
        ("Coupling power", "Coupling power (mW)", "coupling_power_mw", "1.23")]:
        next(w for w in app.selectbox if w.label == "Compare").select(choice).run()
        inputs = [w for w in app.text_input if w.label == label]
        assert len(inputs) == 1
        inputs[0].set_value(value).run()
        assert app.session_state["lambda__"+canonical] == float(value)
        assert not app.exception and not app.slider and not app.error
