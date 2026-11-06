"""The model spectra workspace. The component reads prepared samples only."""
import hashlib
import json

import streamlit as st

from gabes_ui import export as ui_export


def comparison_context(scheme, params, schema):
    """Selection events from an older grid must not apply to a new run."""
    identity = [scheme.name, scheme.cache_version, schema, params]
    return hashlib.sha256(json.dumps(identity, sort_keys=True, allow_nan=False).encode()).hexdigest()


def comparison_plot_payloads(dataset):
    """Export the full native T and OD curves shown by the model workspace."""
    plots = []
    for axis_index, (key, title) in enumerate((
            ("transmission", "Transmission T"),
            ("optical_depth", "Optical depth OD"))):
        traces = []
        for curve_index, condition in enumerate(dataset["series"]):
            power = condition["power_mw"]
            label = "Pump off (0 mW)" if power == 0 else f"Pump {power} mW"
            traces.append({
                "axis_index": axis_index,
                "curve_index": curve_index,
                "figure_curve_label": label,
                "x_label": dataset["axis"]["label"],
                "x_unit": dataset["axis"]["unit"],
                "y_label": title,
                "y_unit": "1",
                "x_data": dataset["x"],
                "y_data": condition[key],
                "color": "#4c606c" if power == 0 else "#07655e",
                "quantity_metadata": {
                    "power_mw": power,
                    "frequency_reference": dataset["axis"]["reference"],
                    "sample_scope": "full native sampled range",
                },
            })
        plots.append({"label": title, "source": "model_comparison", "traces": traces})
    return plots


@st.fragment
def render_comparison(dataset, context, *, scheme, raw, params, view, readout_version):
    from gabes_ui.conditions import render_power_inputs
    from gabes_ui.model_compare_plot import render_group_control, render_model_compare

    key = f"model-compare-{scheme.name}"
    with st.container(key="gabes_spectrum_toolbar"):
        group_column, pump_column, probe_column = st.columns([2.4, 1, 1], gap="small")
        with group_column:
            render_group_control(dataset, context=context, key=key)
        current = render_power_inputs(scheme, scheme.param_schema(), columns=(pump_column, probe_column))
    if current != params:
        st.rerun(scope="app")
    selection = render_model_compare(dataset, key=key, context=context, show_groups=False)
    export_box = st.popover("Export", icon=":material/download:", type="tertiary")
    ui_export.render_export(export_box, scheme=scheme, raw=raw, params=params,
                            view=view, plot_payloads=comparison_plot_payloads(dataset), comparison_payload=None,
                            cache_version=scheme.cache_version, readout_cache_version=readout_version,
                            model_comparison={**dataset, "selection": selection})
