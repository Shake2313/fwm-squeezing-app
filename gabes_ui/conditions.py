"""Shared main-panel editors using canonical parameters and unit-aware entry."""
import math

import streamlit as st

from gabes_ui.controls import current_params, param_visible, selected_atom_row, skey
from gabes_ui.units import display_text, parse_entry


def _field_value(spec, draft):
    if spec.control == "checkbox":
        return bool(draft)
    if spec.choices is not None:
        if draft not in spec.choices:
            raise ValueError(f"Select a valid {spec.label}.")
        return draft
    value, error = parse_entry(str(draft), spec.unit)
    if error or not math.isfinite(value):
        raise ValueError(f"{spec.label}: {error or 'enter a finite number'}.")
    if spec.vmin is not None and value < spec.vmin:
        raise ValueError(f"{spec.label}: minimum {spec.vmin:g} {spec.unit}.".replace(" .", "."))
    if spec.vmax is not None and value > spec.vmax:
        raise ValueError(f"{spec.label}: maximum {spec.vmax:g} {spec.unit}.".replace(" .", "."))
    if isinstance(spec.default, int) and not isinstance(spec.default, bool):
        if not value.is_integer():
            raise ValueError(f"{spec.label}: enter a whole number.")
        return int(value)
    return value


def _write(scheme_name, changes):
    for name, value in changes.items():
        st.session_state[skey(scheme_name, name)] = value


def _commit_draft(key, committed):
    saved = st.session_state.get(key + "-draft", {})
    st.session_state[key + "-draft"] = dict(saved, source=committed)


def _apply_fields(scheme_name, fields, error_key):
    try:
        changes = {sp.name: _field_value(sp, st.session_state[widget_key])
                   for sp, widget_key, key in fields}
    except ValueError as exc:
        st.session_state[error_key] = str(exc)
        return
    _write(scheme_name, changes)
    for spec, widget_key, key in fields:
        _commit_draft(key, changes[spec.name])
    st.session_state.pop(error_key, None)


def _draft(key, committed, display):
    # User edits keep their HTML control and exact draft, including invalid text.
    # External changes (presets or atom/isotope changes) start a fresh control
    # with an explicit default, rather than mixing a default and a state update.
    state_key = key + "-draft"
    saved = st.session_state.get(state_key, {
        "source": st.session_state.get(key + "-source", committed), "revision": 0})
    if saved["source"] != committed:
        saved = dict(source=committed, revision=saved["revision"] + 1)
        st.session_state.pop(key + "-error", None)
    st.session_state[state_key] = saved
    revision = saved["revision"]
    widget_key = f"{key}-refresh-{revision}" if revision else key
    return widget_key, st.session_state.get(widget_key, display)


def _numeric_label(sp):
    return f"{sp.label} ({sp.unit})" if sp.unit else sp.label


def _input(scheme_name, sp, current, key, *, label=None):
    rendered = (display_text(current, None)
                if sp.control != "checkbox" and sp.choices is None else current)
    widget_key, draft = _draft(key, current, rendered)
    options = dict(on_change=_apply_fields,
                   args=(scheme_name, [(sp, widget_key, key)], key+"-error"))
    if sp.control == "checkbox":
        st.checkbox(sp.label, value=bool(draft), key=widget_key, help=sp.help or None, **options)
    elif sp.choices is not None:
        choices = tuple(sp.choices)
        st.selectbox(sp.label, choices, index=choices.index(draft),
                     key=widget_key, help=sp.help or None,
                     format_func=lambda value: (sp.choice_labels or {}).get(value, value), **options)
    else:
        st.text_input(label or _numeric_label(sp), value=draft, key=widget_key,
                      help=sp.help or None, **options)
    error = st.session_state.get(key+"-error")
    if error:
        st.error(error)


def _field_popover(scheme, title, fields, params):
    if not fields:
        return
    with st.popover(title):
        for sp in fields:
            key = f"conditions-{scheme.name}-{sp.name}"
            _input(scheme.name, sp, params[sp.name], key)


def _apply_atom(scheme_name, rows, current, atom_key, field="species", draft_key=None):
    if current and current[0] == st.session_state[atom_key]:
        return
    row = next(r for r in rows if r[0] == st.session_state[atom_key])
    _write(scheme_name, {field: row[2], **row[3]})
    if draft_key:
        _commit_draft(draft_key, row[0])


def _apply_isotope(scheme_name, rows, iso_key, field="species", draft_key=None):
    row = next(r for r in rows if r[1] == st.session_state[iso_key])
    changes = {field: row[2], **row[3]}
    _write(scheme_name, changes)
    if draft_key:
        _commit_draft(draft_key, (row[0], row[1]))


def _atom_editors(scheme, spec, params):
    """Native Atom/Isotope menus inside the shared Atom popover."""
    rows = spec.atom_menu
    # Legacy generic sessions retain their model until an atom is selected.
    value = params[spec.name]
    current = (selected_atom_row(rows, value, params.get)
               if any(row[2] == value for row in rows) else None)
    atom = current[0] if current else value
    atoms = list(dict.fromkeys(row[0] for row in rows))
    key = f"conditions-{scheme.name}-atom"
    widget_key, draft = _draft(key+"-value", atom, atom if atom in atoms else atoms[0])
    st.selectbox("Atom", atoms, index=atoms.index(draft), key=widget_key,
                 disabled=len(atoms) == 1, on_change=_apply_atom,
                 args=(scheme.name, rows, current, widget_key, spec.name, key+"-value"))
    if not current:
        return
    options = [row for row in rows if row[0] == atom]
    if len(options) == 1:
        return
    isotope = current[1]
    key = f"conditions-{scheme.name}-isotope"
    labels = [row[1] for row in options]
    widget_key, draft = _draft(key+"-value", (atom, isotope), isotope)
    st.selectbox("Isotope", labels, index=labels.index(draft), key=widget_key,
                 help=spec.help or None, on_change=_apply_isotope,
                 args=(scheme.name, options, widget_key, spec.name, key+"-value"))


def atom_condition_descriptor(scheme, specs, params=None):
    """Resolve a scheme-declared Atom group without inventing physical controls.

    ``atom_conditions(params)`` declares ordered ``fields`` (including fields
    currently hidden by their ParamSpec), optional read-only ``fixed`` identity
    entries, the ``line`` parameter used in the summary, a compact ``cell``
    label, and optional per-field ``labels``. All widget choices, validation,
    visibility and compute state continue to come from the existing schema.
    """
    params = current_params(scheme.name, scheme) if params is None else params
    metadata = scheme.atom_conditions(params) or {}
    by_name = {sp.name: sp for sp in specs}
    names = tuple(dict.fromkeys(metadata.get("fields", ())))
    fields = [by_name[name] for name in names]
    visible = [sp for sp in fields if sp.control != "derived_slider"
               and param_visible(scheme.name, sp)]
    atom_spec = next((sp for sp in visible if sp.atom_menu), None)
    fixed = dict(metadata.get("fixed", {}))
    atom = fixed.get("Atom", "")
    isotope = fixed.get("Isotope", "")
    if atom_spec is not None:
        value = params[atom_spec.name]
        if any(row[2] == value for row in atom_spec.atom_menu):
            row = selected_atom_row(atom_spec.atom_menu, value, params.get)
            atom, isotope = row[:2]
        else:
            atom, isotope = str(value), ""
    line_name = metadata.get("line")
    line = (params[line_name] if line_name in {sp.name for sp in visible}
            else fixed.get("Transition line", ""))
    parts = [str(value) for value in (isotope if isotope != atom else "", atom,
                                     line, metadata.get("cell", "")) if value]
    return dict(names=frozenset(names), fields=visible, fixed=fixed,
                labels=metadata.get("labels", {}),
                title="Atom" + (" · " + " ".join(parts) if parts else ""))


def render_atom_conditions(scheme, specs, params=None):
    """Render one Atom popover and return every grouped name for caller exclusion."""
    params = current_params(scheme.name, scheme) if params is None else params
    descriptor = atom_condition_descriptor(scheme, specs, params)
    if not descriptor["fields"] and not descriptor["fixed"]:
        return descriptor["names"]
    with st.popover(descriptor["title"]):
        for label, value in descriptor["fixed"].items():
            st.caption(f"{label} · {value}")
        for sp in descriptor["fields"]:
            if sp.atom_menu:
                _atom_editors(scheme, sp, params)
            else:
                _input(scheme.name, sp, params[sp.name],
                       f"conditions-{scheme.name}-{sp.name}",
                       label=descriptor["labels"].get(sp.name))
    return descriptor["names"]


def _value_label(sp, value):
    shown = ("On" if value else "Off") if isinstance(value, bool) else (
        f"{value:.5g}" if isinstance(value, (int, float)) and abs(value) >= 1e6 else
        display_text(value, None) if isinstance(value, (int, float)) else
        str((sp.choice_labels or {}).get(value, value)))
    return f"{shown} {sp.unit}".strip()


def _workspace_inputs(scheme, params):
    declared = getattr(scheme, "comparison_input_names", None)
    return declared(params) if callable(declared) else scheme.workspace_inputs


def render_conditions(scheme, specs):
    """Mount automatic editors before compute; presets and exports share their state."""
    params = current_params(scheme.name, scheme)
    workspace = set(_workspace_inputs(scheme, params)) | set(scheme.model_comparison_controls)
    main = [sp for sp in specs if not sp.advanced and sp.control != "derived_slider"
            and not (sp.applies_defaults and sp.control == "segmented")
            and param_visible(scheme.name, sp)]
    with st.container(key="gabes_conditions"):
        st.html('<h3 class="gabes-conditions-heading">Conditions</h3>')
        with st.container(horizontal=True, gap="small", key="gabes_conditions_options"):
            grouped = render_atom_conditions(scheme, specs, params)
            for sp in main:
                if sp.name not in grouped and sp.name not in workspace:
                    _field_popover(scheme, f"{sp.label} · {_value_label(sp, params[sp.name])}",
                                   [sp], params)
            advanced = [sp for sp in specs if sp.advanced and sp.control != "derived_slider"
                        and sp.name not in grouped and sp.name not in workspace
                        and param_visible(scheme.name, sp)]
            _field_popover(scheme, "Advanced", advanced, params)
    return params



def render_comparison_controls(scheme, specs):
    """Place declared comparison controls beside the comparison workspace."""
    params = current_params(scheme.name, scheme)
    grouped = atom_condition_descriptor(scheme, specs, params)["names"]
    excluded = grouped | set(_workspace_inputs(scheme, params))
    by_name = {sp.name: sp for sp in specs}
    live_controls = getattr(scheme, "comparison_control_names", None)
    for name in (live_controls(params) if callable(live_controls) else scheme.model_comparison_controls):
        sp = by_name[name]
        if name not in excluded and sp.control != "derived_slider" and param_visible(scheme.name, sp):
            _input(scheme.name, sp, params[name], f"comparison-input-{scheme.name}-{name}")
    return params


def render_power_inputs(scheme, specs, *, columns=None):
    """Place a scheme's directly editable physical inputs beside its plots."""
    params = current_params(scheme.name, scheme)
    beams = [sp for name in _workspace_inputs(scheme, params)
             for sp in specs if sp.name == name and param_visible(scheme.name, sp)]
    if not beams:
        return params
    slots = columns if columns is not None else st.columns(len(beams))
    for column, sp in zip(slots, beams):
        with column:
            name = {"pump_power_mw": "Pump power", "probe_power_uw": "Probe power"}.get(sp.name, sp.label)
            _input(scheme.name, sp, params[sp.name], f"power-input-{scheme.name}-{sp.name}",
                   label=f"{name} ({sp.unit})")
    return params
