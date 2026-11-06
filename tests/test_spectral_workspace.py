"""Shared workspaces preserve native solved samples across navigation and exports."""
import base64
from copy import deepcopy
import csv
from io import BytesIO, StringIO
import json
import math
from pathlib import Path
import shutil
import subprocess
import zipfile

import pytest

from gabes_ui import spectral_compare as workspace
from gabes_ui.export import _build_plot_zip_bytes, _curve_csv_bytes


@pytest.fixture
def comparison():
    return {
        "schema": "gabes-spectral-comparison-v1",
        "scheme": "test",
        "mode": "native grids",
        "axis": {"label": "Probe detuning", "unit": "MHz",
                 "reference": "Unshifted transition"},
        "series": [
            {"id": "wide", "label": "Wide scan",
             "parameter": {"key": "power", "label": "Power", "value": 1., "unit": "mW"},
             "x": [-3., -1., 0., 1., 3.],
             "observables": {"transmission": [.2, .45, .8, .7, .2],
                             "rotation": [-.4, None, .1, .3, .5]},
             "readouts": [{"key": "width", "label": "Width", "value": 1.2, "unit": "MHz"}]},
            {"id": "narrow", "label": "Narrow scan",
             "parameter": {"key": "power", "label": "Power", "value": 2., "unit": "mW"},
             "x": [-2., -.75, .25, 1.75, 2.],
             "observables": {"transmission": [.1, .5, None, .6, .15],
                             "rotation": [-.1, -.04, None, .08, .1]},
             "readouts": [{"key": "width", "label": "Width", "value": None, "unit": "MHz",
                           "status": "unresolved", "reason": "No resolved feature"}]},
        ],
        "panels": [
            {"key": "transmission", "title": "Transmission", "y_label": "T", "unit": "1",
             "curves": [{"series_id": "wide", "observable": "transmission"},
                        {"series_id": "narrow", "observable": "transmission"}]},
            {"key": "rotation", "title": "Rotation", "y_label": "Rotation", "unit": "rad",
             "reference": 0.,
             "curves": [{"series_id": "narrow", "observable": "rotation",
                         "label": "Second condition rotation"}]},
        ],
        "summary": {
            "axis": {"label": "RF drive", "unit": "MHz", "reference": "Input Rabi frequency"},
            "series": [
                {"id": "input", "label": "Applied field", "x": [0., 1., 3.7, 7.4],
                 "observables": {"field": [0., .02, .15, .3]}},
                {"id": "inferred", "label": "Peak-derived field", "x": [0., 1., 3.7, 7.4],
                 "observables": {"field": [None, None, .13, .286]}},
            ],
            "panels": [{"key": "field", "title": "Input and inferred field",
                        "y_label": "E", "unit": "mV/m RMS",
                        "curves": [{"series_id": "input", "observable": "field"},
                                   {"series_id": "inferred", "observable": "field"}]}],
        },
        "fixed": {"temperature": 100.},
        "varied": {"key": "power", "values": [1., 2.]},
        "source": {"model": "test model"},
        "notes": ["No interpolation"],
    }


def test_navigation_union_uses_native_coordinates_and_ignores_summary(comparison):
    before = deepcopy(comparison)
    assert workspace.navigation_grid(comparison) == [
        -3., -2., -1., -.75, 0., .25, 1., 1.75, 2., 3.]
    assert comparison == before
    # Summary points (RF drive here) have a different physical axis.
    assert 3.7 not in workspace.navigation_grid(comparison)


def test_navigation_filters_nonfinite_and_missing_coordinates(comparison):
    comparison["series"][0]["x"] += [float("nan"), float("inf"), None, "bad"]
    assert workspace.navigation_grid(comparison) == [
        -3., -2., -1., -.75, 0., .25, 1., 1.75, 2., 3.]


def test_default_selection_uses_common_scan_domain(comparison):
    assert workspace.selection_for(comparison, None, "run") == {
        "x": 0., "window": [-2., 2.], "context": "run"}


def test_saved_cursor_uses_nearest_native_point_and_stale_context_resets(comparison):
    saved = {"x": .2, "window": [-1., 1.], "context": "old"}
    assert workspace.selection_for(comparison, saved, "old") == {
        "x": .25, "window": [-1., 1.], "context": "old"}
    assert workspace.selection_for(comparison, saved, "new") == {
        "x": 0., "window": [-2., 2.], "context": "new"}
    assert saved["x"] == .2


@pytest.mark.parametrize("window", [[2., -2.], [-4., 1.], [-1., 4.],
                                    [float("nan"), 1.], [0.], None])
def test_invalid_saved_windows_reset_to_common_native_domain(comparison, window):
    selected = workspace.selection_for(
        comparison, {"x": 100., "window": window, "context": "run"}, "run")
    assert selected == {"x": 2., "window": [-2., 2.], "context": "run"}


def test_empty_band_expands_only_to_nearest_native_sample(comparison):
    selected = workspace.selection_for(
        comparison, {"x": .17, "window": [.15, .19], "context": "run"}, "run")
    assert selected == {"x": .25, "window": [.15, .25], "context": "run"}
    assert selected["window"][0] <= selected["x"] <= selected["window"][1]


def test_nearest_cursor_stays_inside_nonempty_custom_band(comparison):
    selected = workspace.selection_for(
        comparison, {"x": .27, "window": [.26, 1.1], "context": "run"}, "run")
    # The globally nearest .25 is outside this zoom; 1 is the available native point.
    assert selected == {"x": 1., "window": [.26, 1.1], "context": "run"}


def test_empty_and_single_sample_navigation_have_defined_selections(comparison):
    comparison["series"] = [{"id": "empty", "x": []}]
    assert workspace.selection_for(comparison, None, "run") is None
    comparison["series"].append({"id": "one", "x": [4.]})
    assert workspace.selection_for(comparison, None, "run") == {
        "x": 4., "window": [3.5, 4.5], "context": "run"}


def test_exports_preserve_exact_arrays_units_missing_values_and_panel_order(comparison):
    before = deepcopy(comparison)
    plots = workspace.comparison_plots(comparison)
    assert [(p["label"], p["source"]) for p in plots] == [
        ("Transmission", "spectral_comparison"),
        ("Rotation", "spectral_comparison"),
        ("Input and inferred field", "derived_summary")]
    wide, narrow = plots[0]["traces"]
    assert wide["x_data"] == [-3., -1., 0., 1., 3.]
    assert narrow["x_data"] == [-2., -.75, .25, 1.75, 2.]
    assert wide["y_data"] == [.2, .45, .8, .7, .2]
    assert narrow["y_data"] == [.1, .5, None, .6, .15]
    assert narrow["quantity_metadata"] == {
        "condition_id": "narrow", "frequency_reference": "Unshifted transition",
        "observable": "transmission"}
    rotation = plots[1]["traces"][0]
    assert rotation["figure_curve_label"] == "Second condition rotation"
    assert rotation["x_data"] == narrow["x_data"]
    assert rotation["y_data"] == [-.1, -.04, None, .08, .1]
    assert rotation["y_unit"] == "rad"
    inferred = plots[2]["traces"][1]
    assert inferred["x_label"] == "RF drive" and inferred["x_unit"] == "MHz"
    assert inferred["y_unit"] == "mV/m RMS"
    assert inferred["y_data"] == [None, None, .13, .286]
    assert comparison == before
    json.dumps(plots, allow_nan=False)


def test_export_csv_zip_retains_native_rows_and_unresolved_gaps(comparison):
    plots = workspace.comparison_plots(comparison)
    narrow = plots[0]["traces"][1]
    rows = list(csv.reader(StringIO(_curve_csv_bytes(narrow).decode("utf-8"))))
    assert rows == [
        ["Probe detuning [MHz]", "T [1]"],
        ["-2.0", "0.1"], ["-0.75", "0.5"], ["0.25", ""],
        ["1.75", "0.6"], ["2.0", "0.15"]]
    with zipfile.ZipFile(BytesIO(_build_plot_zip_bytes(plots))) as archive:
        exported = [list(csv.reader(StringIO(archive.read(name).decode("utf-8"))))
                    for name in archive.namelist()]
    assert len(exported) == 5
    assert rows in exported
    assert [["RF drive [MHz]", "E [mV/m RMS]"],
            ["0.0", ""], ["1.0", ""], ["3.7", "0.13"], ["7.4", "0.286"]] in exported


def test_absent_summary_exports_only_spectral_panels(comparison):
    del comparison["summary"]
    plots = workspace.comparison_plots(comparison)
    assert len(plots) == 2
    assert all(p["source"] == "spectral_comparison" for p in plots)


def test_remount_preserves_accepted_window_and_rejects_late_axis_gestures(comparison, monkeypatch):
    session = {"spectral-test": {"selection": {
        "x": 1., "window": [.2, 1.8], "context": "power1"}}}
    monkeypatch.setattr(workspace.st, "session_state", session)
    calls = []
    monkeypatch.setattr(workspace, "_COMPONENT", lambda **kwargs: calls.append(kwargs))
    first = workspace._render_plot(comparison, "power1")
    # A delayed gesture from power1 cannot replace the last accepted inspection.
    session["spectral-test"]["selection"].update(x=-2., window=[-3., -1.])
    carried = workspace._render_plot(comparison, "power2")
    assert carried == dict(first, context="power2")
    assert calls[-1]["data"]["selection"] == carried
    assert calls[-1]["data"]["grid"] == workspace.navigation_grid(comparison)
    comparison["axis"]["reference"] = "Different transition reference"
    reset = workspace._render_plot(comparison, "other-axis")
    assert reset == workspace.workspace_selection(comparison, None, "other-axis")
    assert reset["series_id"] == "wide" and reset["panel_key"] == "transmission"


def test_component_module_loads_and_accepts_selection_in_apptest(comparison):
    from streamlit.testing.v1 import AppTest

    script = """
import json
import streamlit as st
from gabes_ui.spectral_compare import _render_plot
dataset = json.loads(DATASET)
selected = _render_plot(dataset, "apptest")
st.json(selected)
""".replace("DATASET", repr(json.dumps(comparison)))
    app = AppTest.from_string(script).run()
    assert not app.exception
    selected = app.session_state["spectral-test-accepted"]
    assert selected["context"] == "apptest"
    assert selected["x"] == 0. and selected["window"] == [-2., 2.]
    app.run()
    assert not app.exception


def test_component_serializes_nonfinite_failed_solve_diagnostics_as_null(comparison, monkeypatch):
    comparison["series"][0]["scope"] = {"atomic_diagnostics": {"health": {
        "virtual_excitation_bound": float("nan"), "derivative_error": float("inf"),
        "ellipticity_drift_length": -float("inf")}}}
    monkeypatch.setattr(workspace.st, "session_state", {})
    calls = []
    monkeypatch.setattr(workspace, "_COMPONENT", lambda **kwargs: calls.append(kwargs))
    workspace._render_plot(comparison, "failed-solve")
    json.dumps(calls[-1]["data"], allow_nan=False)
    health = calls[-1]["data"]["dataset"]["series"][0]["scope"]["atomic_diagnostics"]["health"]
    assert health == dict.fromkeys(health)
    assert math.isnan(comparison["series"][0]["scope"]["atomic_diagnostics"]["health"]["virtual_excitation_bound"])


def test_actual_shared_scale_and_native_sampling_in_node():
    node = shutil.which("node")
    if not node:
        pytest.skip("Node.js required for shared spectrum scale tests")
    source = (Path(__file__).resolve().parents[1] /
              "assets/ui/spectrum_plot.js").read_text(encoding="utf-8")
    # Expose the existing internal sampler for this behavioral test only.
    url = "data:text/javascript;base64," + base64.b64encode(
        (source + "\nexport {curveSample};\n").encode()).decode()
    script = r"""
import assert from 'node:assert/strict';
import {spectrumScale, nearestSample, curveSample, sampleInsidePlot, sampleScope, sampleWarningIntervals} from 'MODULE';
const near=(actual,expected)=>assert.ok(
  Math.abs(actual-expected) <= 1e-10 * Math.max(1,Math.abs(expected)),
  String(actual)+' != '+String(expected));
const linear=spectrumScale(-12,18,52,472);
near(linear.forward(-12),52); near(linear.forward(18),472);
near(linear.forward(3),262);
for(const x of [-12,-7.31,0,3,18]) near(linear.inverse(linear.forward(x)),x);
const log=spectrumScale(.001,1000,52,472,true);
near(log.forward(.001),52); near(log.forward(1000),472);
near(log.forward(1),262); near(log.inverse(262),1);
for(const x of [.001,.0042,.1,1,16,1000]) near(log.inverse(log.forward(x)),x);
for(const x of [.001,.01,.1,1,10,100]) near(log.forward(x*10)-log.forward(x),70);
const ordinary=spectrumScale(.001,1000,52,472,false);
near(ordinary.inverse(262),500.0005);
assert.notEqual(log.forward(1),ordinary.forward(1));
const shared=[-3,-1,0,1,2,3,5], native={x:[-1,1,3],y:[10,null,30]};
assert.equal(nearestSample(shared,.5),2); // deterministic lower sample on an exact tie
assert.deepEqual(curveSample(native,shared,-.8),{x:-1,y:10});
assert.deepEqual(curveSample(native,shared,.6),{x:1,y:null}); // do not interpolate across gaps
assert.deepEqual(curveSample(native,shared,2),{x:1,y:null});
assert.deepEqual(curveSample(native,shared,2.1),{x:3,y:30});
assert.deepEqual(curveSample(native,shared,4),{x:null,y:null}); // no invented sampled coordinate
assert.deepEqual(curveSample({x:[],y:[]},shared,0),{x:null,y:null});
assert.deepEqual(curveSample({y:[20,21,22,23,24,25,26]},shared,1.8),{x:2,y:24});
const geometry={left:52,right:587,top:23,bottom:145,
  sx:x=>52+(x-46.65752)/6*535,sy:y=>145-y*100};
assert.equal(sampleInsidePlot({x:46.42768,y:.0000351679},geometry),false);
assert.equal(sampleInsidePlot({x:46.65752,y:.2},geometry),true);
assert.equal(sampleInsidePlot({x:52.65752,y:2},geometry),false);
assert.equal(sampleInsidePlot({x:null,y:null},geometry),false);
const scoped={x:[0,1,2,3,4],sample_status:['available','available outside approximation screens','available','unavailable','available'],
  sample_reasons:[[],['uniform-pump screen'],[],['atomic solve failed'],[]]};
assert.deepEqual(sampleScope(scoped,[],1),{status:'available outside approximation screens',reasons:['uniform-pump screen']});
assert.deepEqual(sampleScope(scoped,[],3),{status:'unavailable',reasons:['atomic solve failed']});
assert.deepEqual(sampleScope(scoped,[],null),{status:'outside scan',reasons:[]});
assert.deepEqual(sampleWarningIntervals([scoped,scoped]),[[.5,1.5],[2.5,3.5]]);
assert.deepEqual(sampleWarningIntervals([{x:[0,1]}]),[]);
""".replace("MODULE", url)
    result = subprocess.run([node, "--input-type=module"], input=script,
                            text=True, encoding="utf-8", capture_output=True, timeout=20)
    assert result.returncode == 0, result.stderr[-4000:]
