"""Shared sampled-curve workspace for the registered scientific schemes."""
from bisect import bisect_left, bisect_right
from pathlib import Path
import hashlib
import json
import math
from threading import current_thread

import streamlit as st

from gabes_ui.component_runtime import missing_apptest_component, spectrum_component_js
from gabes_ui.conditions import render_comparison_controls, render_power_inputs
from gabes_ui import export as ui_export
from gabes_ui.controls import skey
from gabes_ui.spectral_sweep import (
    RowSweep, assemble_dataset, digest, extract_row, resolve_row, row_value,
    scan_values, workspace_key,
)

_ASSETS = Path(__file__).resolve().parents[1] / "assets" / "ui"
_COMPONENT = st.components.v2.component(
    "gabes_spectral_compare",
    html='<section class="mc sc" aria-label="Spectral comparison"></section>',
    css=(_ASSETS / "model_compare.css").read_text(encoding="utf-8") + "\n"
        + (_ASSETS / "spectral_compare.css").read_text(encoding="utf-8"),
    js=spectrum_component_js("spectral_compare.js"),
)


def navigation_grid(dataset):
    return sorted({float(x) for curve in dataset["series"] for x in curve["x"]
                   if isinstance(x, (int, float)) and math.isfinite(x)})


def selection_for(dataset, entry, context):
    grid = navigation_grid(dataset)
    if not grid:
        return None
    lo, hi = grid[0], grid[-1]
    saved = entry if isinstance(entry, dict) and entry.get("context") == context else {}
    window = saved.get("window")
    if not (isinstance(window, (list, tuple)) and len(window) == 2
            and all(isinstance(v, (int, float)) and math.isfinite(v) for v in window)
            and lo <= window[0] < window[1] <= hi):
        common_lo = max(curve["x"][0] for curve in dataset["series"] if curve["x"])
        common_hi = min(curve["x"][-1] for curve in dataset["series"] if curve["x"])
        window = ([common_lo, common_hi] if common_lo < common_hi else
                  [lo, hi] if lo < hi else [lo-.5, hi+.5])
    x = saved.get("x")
    if not isinstance(x, (int, float)) or not math.isfinite(x):
        x = min(max(0., window[0]), window[1])
    x = min(max(float(x), window[0]), window[1])
    index = min(bisect_left(grid, x), len(grid)-1)
    if index and abs(grid[index-1]-x) <= abs(grid[index]-x):
        index -= 1
    first, last = bisect_left(grid, window[0]), bisect_right(grid, window[1])-1
    if first <= last:
        index = min(max(index, first), last)
    else:
        window = [min(window[0], grid[index]), max(window[1], grid[index])]
    return {"x": grid[index], "window": list(window), "context": context}


def comparison_plots(dataset):
    """Export precisely the native arrays selected by each plotted panel."""
    result = []
    for block, source in ((dataset, "spectral_comparison"), (dataset.get("summary"), "derived_summary")):
        if not block:
            continue
        by_id = {s["id"]: s for s in block["series"]}
        for panel in block["panels"]:
            traces = []
            for i, definition in enumerate(panel["curves"]):
                curve = by_id[definition["series_id"]]
                traces.append({
                    "axis_index": 0, "curve_index": i,
                    "figure_curve_label": definition.get("label", curve["label"]),
                    "x_label": block["axis"]["label"], "x_unit": block["axis"]["unit"],
                    "y_label": panel["y_label"], "y_unit": panel.get("unit", ""),
                    "x_data": curve.get("observable_x", {}).get(definition["observable"], curve["x"]), "y_data": curve["observables"][definition["observable"]],
                    "quantity_metadata": {"condition_id": curve["id"],
                                          "frequency_reference": block["axis"].get("reference", ""),
                                          "observable": definition["observable"]},
                })
            result.append({"label": panel["title"], "source": source, "traces": traces})
    return result


def workspace_selection(dataset, entry, context, params=None):
    if not (isinstance(entry, dict) and entry.get("context") == context):
        defaults = dataset.get("display", {})
        entry = dict(context=context, window=defaults.get("window"),
                     display_mode=defaults.get("default_mode"), condition_axis=defaults.get("condition_axis"),
                     profile=defaults.get("profile"))
    selection = selection_for(dataset, entry, context)
    if selection is None:
        return None
    saved = entry if isinstance(entry, dict) and entry.get("context") == context else {}
    rows = dataset["series"]
    picked = next((row for row in rows if row["id"] == saved.get("series_id")), None)
    if picked is None and saved.get("condition"):
        picked = next((row for row in rows if row.get("parameter") == saved["condition"]), None)
    key = dataset.get("varied", {}).get("key")
    if picked is None and params is not None and key in params:
        picked = next((row for row in rows if row_value(row) == params[key]), None)
    picked = picked or rows[0]
    native = [float(x) for x in picked["x"] if isinstance(x, (int, float)) and math.isfinite(x)
              and selection["window"][0] <= x <= selection["window"][1]]
    if native:
        selection["x"] = min(native, key=lambda x: (abs(x-selection["x"]), x))
    panel = next((p for p in dataset["panels"] if p["key"] == saved.get("panel_key")), dataset["panels"][0])
    result = dict(selection, series_id=picked["id"], condition=picked.get("parameter"),
                  panel_key=panel["key"], display_mode="map" if saved.get("display_mode") == "map" else "curves")
    axes = dataset.get("condition_axes", [])
    if axes:
        result["condition_axis"] = next((a["key"] for a in axes if a["key"] == saved.get("condition_axis")), axes[0]["key"])
    profiles = dataset.get("profiles", [])
    if profiles:
        profile = next((p for p in profiles if p["key"] == saved.get("profile")), profiles[0])
        result["profile"] = profile["key"]
        if saved.get("panel_key") not in profile["panels"]:
            result["panel_key"] = profile["default_panel"]
    if dataset.get("optional_overlays"):
        enabled = saved.get("overlays")
        if not isinstance(enabled, list):
            enabled = [o["key"] for o in dataset["optional_overlays"] if o.get("default")]
        result["overlays"] = [o["key"] for o in dataset["optional_overlays"] if o["key"] in enabled]
    return result


def selected_parameter_change(scheme, dataset, selection, params):
    key = dataset.get("varied", {}).get("key")
    row = next((r for r in dataset["series"] if r["id"] == selection["series_id"]), None)
    spec = next((s for s in scheme.param_schema() if s.name == key), None)
    if row is None or spec is None or row.get("parameter", {}).get("key") != key:
        return {}
    value = row_value(row)
    if spec.choices is not None:
        if value not in spec.choices:
            return {}
    elif (isinstance(value, bool) or not isinstance(value, (int, float)) or
          not math.isfinite(value) or (spec.vmin is not None and value < spec.vmin) or
          (spec.vmax is not None and value > spec.vmax)):
        return {}
    return {key: value} if value != params[key] else {}


def _render_plot(dataset, context, *, scheme=None, params=None):
    # Caught native solve failures may leave nonfinite diagnostics even though
    # the spectrum already contains null gaps. The browser requires strict JSON.
    dataset = ui_export.to_json_primitive(dataset)
    key = "spectral-" + dataset["scheme"]
    widget = st.session_state.get(key, {})
    incoming = widget.get("selection") if hasattr(widget, "get") else None
    fresh = isinstance(incoming, dict) and incoming.get("context") == context
    saved_key = key + "-accepted"
    identity = json.dumps([dataset["scheme"], dataset.get("mode"), dataset["axis"], dataset.get("varied", {}).get("key")], sort_keys=True)
    if not isinstance(incoming, dict) or incoming.get("context") != context:
        prior = st.session_state.get(saved_key, {})
        incoming = dict(prior, context=context) if prior.get("axis_key") == identity else None
        if incoming is not None and params is not None:
            key_name = dataset.get("varied", {}).get("key")
            current = next((r for r in dataset["series"] if row_value(r) == params.get(key_name)), None)
            if current is not None:
                incoming.update(series_id=current["id"], condition=current.get("parameter"))
    selection = workspace_selection(dataset, incoming, context, params)
    if selection is None:
        st.info("No finite sampled axis is available for these conditions.")
        return None
    st.session_state[saved_key] = dict(selection, axis_key=identity)
    if fresh and scheme is not None:
        changes = selected_parameter_change(scheme, dataset, selection, params)
        if changes:
            for name, value in changes.items():
                st.session_state[skey(scheme.name, name)] = value
            st.rerun(scope="app")
    summary = dataset.get("summary")
    if summary:
        for row in summary["series"]:
            if "condition_ids" not in row and len(row["x"]) == len(dataset["series"]):
                row["condition_ids"] = [r["id"] for r in dataset["series"]]
    revision = digest([[(r["id"], len(r["x"])) for r in dataset["series"]], dataset.get("sweep")])
    try:
        _COMPONENT(key=key, data={"dataset": dataset, "context": context, "selection": selection,
                                 "grid": navigation_grid(dataset), "revision": revision},
                   default={"selection": None}, on_selection_change=lambda: None, width="stretch")
    except ValueError as error:
        if not missing_apptest_component(error, "gabes_spectral_compare"):
            raise
    except TypeError as error:
        if "bad argument type for built-in operation" not in str(error):
            raise
    return selection


def _workspace(dataset, scheme, params, readout_version, *, resolve_raw=None, raw=None, view=None):
    key = "spectral-workspace-" + scheme.name
    varied = dataset.get("varied", {}).get("key")
    identity = workspace_key(scheme, params, varied, readout_version)
    entry = st.session_state.get(key)
    if entry is None or entry["identity"] != identity or entry["dataset"] is not dataset:
        if entry is not None:
            entry["runner"].cancel()
        runner = RowSweep()
        rows = []
        for i, row in enumerate(dataset["series"]):
            if not row.get("parameter"):
                continue
            payload = extract_row(dataset, i)
            value = row_value(row)
            if varied in params and row["parameter"].get("key") == varied:
                if value == params[varied] and raw is not None:
                    payload.update(_native_raw=raw, _native_view=view)
                elif resolve_raw is not None:
                    payload["_native_raw"] = resolve_raw({**params, varied: value})
            rows.append((value, payload))
        runner.seed(identity, rows)
        entry = dict(identity=identity, dataset=dataset, template=dataset, runner=runner,
                     revision=(entry or {}).get("revision", 0)+1, dirty=False)
        st.session_state[key] = entry
    return entry


def _scan_controls(entry, scheme, params, resolve_raw):
    data = entry["dataset"]
    key = data.get("varied", {}).get("key")
    spec = next((sp for sp in scheme.param_schema() if sp.name == key), None)
    values = [row_value(row) for row in entry["template"]["series"]]
    if (resolve_raw is None or spec is None or spec.choices is not None or len(values) < 2 or
            not all(isinstance(v, (int, float)) and not isinstance(v, bool) and math.isfinite(v) for v in values)):
        return
    prefix = "spectral-scan-" + scheme.name + entry["identity"][:12] + str(entry["revision"])
    def dirty():
        entry["dirty"] = True
    with st.popover("Condition scan"):
        st.caption(data["varied"].get("label", spec.label) + (f" ({spec.unit})" if spec.unit else ""))
        low = st.number_input("From", value=float(min(values)), key=prefix+"-from", on_change=dirty)
        high = st.number_input("To", value=float(max(values)), key=prefix+"-to", on_change=dirty)
        count = st.number_input("Uniform samples", min_value=2, max_value=41,
                                value=min(41, max(2, len(values))), step=1, key=prefix+"-count", on_change=dirty)
        st.caption("The exact condition at scan start is retained alongside the uniform grid. Native frequency resolution is unchanged.")
        if entry["dirty"]:
            entry["dirty"] = False
            try:
                requested = scan_values(low, high, count, params[key], spec)
            except ValueError as error:
                entry["error"] = str(error)
            else:
                from streamlit.runtime.scriptrunner import add_script_run_ctx, get_script_run_ctx
                ctx = get_script_run_ctx()
                committed = dict(params)
                def resolve(value):
                    if ctx is not None:
                        add_script_run_ctx(current_thread(), ctx)
                    return resolve_row(scheme, committed, key, value, resolve_raw)
                entry.pop("error", None)
                entry["runner"].start(entry["identity"], requested, resolve,
                                      request=dict(low=low, high=high, uniform_samples=int(count),
                                                   anchor_value=params[key]))
                st.rerun(scope="app")
        if entry.get("error"):
            st.error(entry["error"])


def _render_spectral_section(dataset, *, scheme, raw, params, view, readout_version,
                             entry, resolve_raw=None, polling=False):
    snapshot = entry["runner"].snapshot()
    if snapshot is not None:
        entry["dataset"] = assemble_dataset(entry["template"], snapshot)
        dataset = entry["dataset"]
        if polling and snapshot["status"] != "running":
            st.rerun(scope="app")
    with st.container(key="gabes_spectrum_toolbar"):
        render_comparison_controls(scheme, scheme.param_schema())
        current = render_power_inputs(scheme, scheme.param_schema())
    if current != params:
        st.rerun(scope="app")
    _scan_controls(entry, scheme, params, resolve_raw)
    if snapshot:
        st.progress(snapshot["completed"]/max(1, snapshot["total"]),
                    text=f'{snapshot["completed"]}/{snapshot["total"]} conditions · {snapshot["status"]} · {snapshot["cached"]} cached')
        if snapshot["status"] == "running" and st.button("Cancel scan", key="cancel-scan-"+scheme.name):
            entry["runner"].cancel()
            st.rerun(scope="app")
        if snapshot["status"] in {"cancelled", "failed"} and st.button("Resume scan", key="resume-scan-"+scheme.name):
            entry["runner"].resume()
            st.rerun(scope="app")
        if snapshot["error"]:
            st.error("Condition scan stopped: " + snapshot["error"])
    context = hashlib.sha256(json.dumps([scheme.cache_version, readout_version, params],
                                       sort_keys=True, allow_nan=False).encode()).hexdigest()
    selection = _render_plot(dataset, context, scheme=scheme, params=params)
    export_box = st.popover("Export", icon=":material/download:", type="tertiary")
    ui_export.render_export(
        export_box, scheme=scheme, raw=raw, params=params, view=view,
        plot_payloads=comparison_plots(dataset), comparison_payload=None,
        cache_version=scheme.cache_version, readout_cache_version=readout_version,
        model_comparison={**dataset, "selection": selection})


def render_spectral_section(dataset, *, scheme, raw, params, view, readout_version, resolve_raw=None):
    entry = _workspace(dataset, scheme, params, readout_version,
                       resolve_raw=resolve_raw, raw=raw, view=view)
    snapshot = entry["runner"].snapshot()
    polling = snapshot is not None and snapshot["status"] == "running"
    fragment = st.fragment(run_every=.75 if polling else None)(_render_spectral_section)
    return fragment(dataset, scheme=scheme, raw=raw, params=params, view=view,
                    readout_version=readout_version, entry=entry, resolve_raw=resolve_raw, polling=polling)
