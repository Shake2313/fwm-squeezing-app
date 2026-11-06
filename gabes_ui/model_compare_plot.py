"""Interactive rendering of the solved pump comparison. No solver calls live here."""
from array import array
from bisect import bisect_left, bisect_right
from collections.abc import Mapping
from hashlib import sha256
from math import isfinite
from pathlib import Path

import streamlit as st

from .component_runtime import missing_apptest_component, spectrum_component_js

_ASSETS = Path(__file__).resolve().parents[1] / "assets" / "ui"
_COMPONENT = st.components.v2.component(
    "gabes_model_compare",
    html='<section class="mc" aria-label="Model comparison"></section>',
    css=(_ASSETS / "model_compare.css").read_text(encoding="utf-8"),
    js=spectrum_component_js("model_compare.js"),
)


def _sample_nearest(grid, frequency):
    right = min(bisect_left(grid, frequency), len(grid) - 1)
    if right and abs(frequency - grid[right - 1]) <= abs(grid[right] - frequency):
        return right - 1
    return right


def comparison_selection(dataset, entry, context):
    """Validate a client selection against this calculation before export."""
    saved = entry if isinstance(entry, Mapping) and entry.get("context") == context else {}
    groups = dataset["groups"]
    group = next((g for g in groups if g["id"] == saved.get("group_id")),
                 groups[0] if groups else None)
    grid = dataset["x"]
    limits = group["window_mhz"] if group else (grid[0], grid[-1])
    window = saved.get("window_mhz")
    if (isinstance(window, (tuple, list)) and len(window) == 2
            and all(type(v) in (float, int) and isfinite(v) for v in window)
            and grid[0] <= window[0] < window[1] <= grid[-1]
            and bisect_left(grid, window[0]) < bisect_right(grid, window[1])):
        limits = window
        group = next((g for g in groups if list(g["window_mhz"]) == list(limits)), None)
    index = saved.get("sample_index")
    if (type(index) is not int or not 0 <= index < len(grid)
            or not limits[0] <= grid[index] <= limits[1]):
        index = max(bisect_left(grid, limits[0]), min(bisect_right(grid, limits[1])-1,
                    _sample_nearest(grid, sum(limits) / 2)))
    return dict(group_id=group["id"] if group else None, power_mw=dataset["default_power_mw"],
                sample_index=index, x_mhz=grid[index], window_mhz=list(limits), show_powers=False)


def _current_selection(dataset, key, context):
    """Resolve acknowledged gestures and authoritative native toolbar commands."""
    state = st.session_state.get(key, {})
    entry = state.get("selection") if hasattr(state, "get") else getattr(state, "selection", None)
    revision = st.session_state.get(key + "-navigation-revision", 0)
    if isinstance(entry, Mapping) and entry.get("navigation_revision", 0) != revision:
        entry = None
    command = st.session_state.get(key + "-native-selection")
    if command and command.get("context") == context:
        if isinstance(entry, Mapping) and entry.get("context") == context and entry.get("navigation_revision", 0) == revision:
            st.session_state.pop(key + "-native-selection", None)
        else:
            entry = command
    # A power change can produce new spectra on the same physical axis. Only
    # the server's accepted selection can carry over; stale browser gestures
    # still belong to the previous calculation and are discarded.
    axis = dataset["axis"]
    axis_key = (dataset["scheme"], axis["unit"], axis["reference"],
                sha256(array("d", dataset["x"]).tobytes()).digest())
    saved_key = key + "-committed-selection"
    if not isinstance(entry, Mapping) or entry.get("context") != context:
        previous = st.session_state.get(saved_key, {})
        entry = (dict(previous, context=context)
                 if previous.get("axis_key") == axis_key else None)
    selection = comparison_selection(dataset, entry, context)
    st.session_state[saved_key] = dict(selection, context=context, axis_key=axis_key)
    return selection


def _choose_group(dataset, key, context):
    group_key = key + "-group-control"
    group = next((g for g in dataset["groups"] if g["id"] == st.session_state[group_key]), None)
    if group is None:
        return
    window = group["window_mhz"]
    entry = dict(context=context, group_id=group["id"], window_mhz=window)
    selection = comparison_selection(dataset, entry, context)
    revision = st.session_state.get(key + "-navigation-revision", 0) + 1
    st.session_state[key + "-navigation-revision"] = revision
    st.session_state[key + "-native-selection"] = dict(selection, context=context, navigation_revision=revision)


def render_group_control(dataset, *, key="sas-model-compare", context):
    """Native group presets share the lower toolbar with editable beam powers."""
    selection = _current_selection(dataset, key, context)
    group_key = key + "-group-control"
    st.session_state[group_key] = selection["group_id"]
    groups = {g["id"]: g.get("label", g["id"]) for g in dataset["groups"]}
    if groups:
        st.segmented_control("Doppler group", list(groups), key=group_key, format_func=groups.get,
                             on_change=_choose_group, args=(dataset, key, context))
    else:
        st.caption("Doppler group")
        st.caption("Full scan · no grouped transition markers")


def _receive_selection(dataset, key, context):
    selection = _current_selection(dataset, key, context)
    st.session_state[key + "-group-control"] = selection["group_id"]


def render_model_compare(dataset, *, key="sas-model-compare", context, show_groups=True):
    """Mount the browser plot and return its last committed physical selection."""
    selection = _current_selection(dataset, key, context)
    try:
        _COMPONENT(key=key, data=dict(dataset=dataset, context=context, selection=selection,
                                    show_groups=show_groups,
                                    navigation_revision=st.session_state.get(key + "-navigation-revision", 0)),
                   default={"selection": None},
                   on_selection_change=lambda: _receive_selection(dataset, key, context), width="stretch")
    except ValueError as error:
        if not missing_apptest_component(error, "gabes_model_compare"):
            raise
    except TypeError as error:
        # AppTest's v2 registry cannot serialize a component definition.
        if "bad argument type for built-in operation" not in str(error):
            raise
    return selection
