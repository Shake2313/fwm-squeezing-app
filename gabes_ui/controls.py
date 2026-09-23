"""Control rail — the knobs a scheme declares, rendered generically.

The front-end knows no physics: it reads `param_schema()` and renders one
control per `ParamSpec`. Numeric knobs are ScrubFields (drag, type or step;
`gabes_ui.scrub`), everything else is the matching Streamlit widget.

Session state holds one entry per knob under `skey(scheme, name)`, which is
what presets, default sets and `visible_if` read and write.
"""
from html import escape

import streamlit as st

from gabes_ui.format import slider_format
from gabes_ui.scrub import scrub_field

ADVANCED_FOLD = 6     # an advanced-only group longer than this starts folded


def skey(scheme_name, pname):
    return f"{scheme_name}__{pname}"


def current_params(scheme_name, scheme_obj):
    """Live value of every param knob from session_state, falling back to defaults."""
    return {
        sp.name: st.session_state.get(skey(scheme_name, sp.name), sp.default)
        for sp in scheme_obj.param_schema()
    }


def render_group_header(container, group):
    container.markdown(
        f"<div class='gabes-group-header'>{escape(group)}</div>",
        unsafe_allow_html=True,
    )


def apply_recommended_defaults(scheme_name, scheme_obj, key):
    """on_change for a control whose selection applies the matching
    recommended_defaults set (so choosing it also resets that mode's knobs)."""
    selection = st.session_state.get(key)
    cur = current_params(scheme_name, scheme_obj)
    sets = scheme_obj.recommended_defaults(cur) or {}
    defaults = sets.get(selection) or sets.get(cur.get("mode")) or {}
    for k, v in defaults.items():
        st.session_state[skey(scheme_name, k)] = v


def apply_default_set(scheme_name, scheme_obj, label):
    """Write one named recommended_defaults set into session state."""
    cur = current_params(scheme_name, scheme_obj)
    sets = scheme_obj.recommended_defaults(cur) or {}
    for k, v in (sets.get(label) or {}).items():
        st.session_state[skey(scheme_name, k)] = v


def _hover_label(container, label, help_):
    """A knob label that carries its help on hover, like a ScrubField's.

    Streamlit's own `help=` draws a "?" next to every label; across a rail of
    forty knobs those become the loudest thing on screen (docs/ui_redesign
    decisions.md D4). Widgets that cannot host their own label get this
    instead, and pass `label_visibility="collapsed"`.
    """
    title = f" title='{escape(help_, quote=True)}'" if help_ else ""
    container.markdown(
        f"<div class='gabes-knob-label'{title}>{escape(label)}</div>",
        unsafe_allow_html=True)


def _is_integer_knob(sp):
    numbers = [v for v in (sp.default, sp.vmin, sp.vmax, sp.step) if v is not None]
    return bool(numbers) and all(isinstance(v, int) for v in numbers)


def render_param(container, scheme_name, sp, scheme_obj=None):
    key = skey(scheme_name, sp.name)
    help_ = sp.help or None
    has_state = key in st.session_state
    current = st.session_state.get(key, sp.default)
    choice_labels = getattr(sp, "choice_labels", None) or {}
    choice_kwargs = (
        {"format_func": lambda value: choice_labels.get(value, str(value))}
        if choice_labels else {}
    )
    on_change = None
    if getattr(sp, "applies_defaults", False) and scheme_obj is not None:
        on_change = lambda: apply_recommended_defaults(scheme_name, scheme_obj, key)
    control = getattr(sp, "control", "auto")
    # A ScrubField shows its unit beside the value; the other widgets have
    # nowhere to put it, so they keep it in the label.
    label = sp.label + (f"  [{sp.unit}]" if sp.unit else "")
    if control == "checkbox":
        if has_state:
            return container.checkbox(label, key=key, help=help_,
                                      on_change=on_change)
        return container.checkbox(label, value=bool(current), key=key,
                                  help=help_, on_change=on_change)
    if control == "segmented":
        options = list(sp.choices or ())
        _hover_label(container, label, help_)
        choice_kwargs["label_visibility"] = "collapsed"
        if hasattr(container, "segmented_control"):
            try:
                return container.segmented_control(label, options, key=key,
                                                   on_change=on_change,
                                                   **choice_kwargs)
            except TypeError:
                pass
        if has_state:
            return container.radio(label, options, key=key, horizontal=True,
                                   on_change=on_change, **choice_kwargs)
        idx = options.index(current) if current in options else 0
        return container.radio(label, options, key=key, horizontal=True, index=idx,
                               on_change=on_change, **choice_kwargs)
    if sp.choices is not None:
        options = list(sp.choices)
        _hover_label(container, label, help_)
        choice_kwargs["label_visibility"] = "collapsed"
        if has_state:
            return container.selectbox(label, options, key=key,
                                       on_change=on_change, **choice_kwargs)
        idx = options.index(current) if current in options else 0
        return container.selectbox(label, options, index=idx, key=key,
                                   on_change=on_change, **choice_kwargs)
    return scrub_field(
        container, key=key, label=sp.label, value=current, vmin=sp.vmin,
        vmax=sp.vmax, step=sp.step, unit=sp.unit, default=sp.default,
        fmt=slider_format(sp), help=help_,
        endpoints=getattr(sp, "endpoints", None),
        live=not getattr(sp, "recompute", True),
        integer=_is_integer_knob(sp), on_change=on_change,
    )


def param_visible(scheme_name, sp):
    if getattr(sp, "hidden", False):
        return False
    cond = getattr(sp, "visible_if", None)
    if not cond:
        return True
    for pname, allowed in cond.items():
        cur = st.session_state.get(skey(scheme_name, pname))
        if isinstance(allowed, (set, tuple, list)):
            if cur not in allowed:
                return False
        elif cur != allowed:
            return False
    return True


def render_rail(container, scheme, specs, skip=()):
    """Presets, then grouped knobs. Returns the full parameter dict.

    "Show advanced" reveals advanced knobs inside the group they belong to (a
    pump waist next to the pump power); advanced-only groups follow, folded
    when long.
    """
    scheme_presets = scheme.presets()
    if scheme_presets and getattr(scheme, "presets_group", None):
        render_group_header(container, scheme.presets_group)
    for preset in scheme_presets:
        def _apply(p=preset, sname=scheme.name):
            for k, v in p.values.items():
                st.session_state[skey(sname, k)] = v
        container.button(preset.name, on_click=_apply,
                         use_container_width=True, help=preset.help)

    visible = [sp for sp in specs
               if param_visible(scheme.name, sp)
               and not any(sp is other for other in skip)]
    params = {}
    group_order = []
    for sp in visible:
        if not sp.advanced and sp.group not in group_order:
            group_order.append(sp.group)
    advanced = [sp for sp in visible if sp.advanced]
    advanced_by_group = {}
    for sp in advanced:
        advanced_by_group.setdefault(
            getattr(sp, "advanced_group", "") or sp.group, []).append(sp)
    advanced_key = skey(scheme.name, "_show_advanced")
    show_advanced = bool(st.session_state.get(advanced_key, False))

    for group in group_order:
        render_group_header(container, group)
        for sp in visible:
            if sp.group == group and not sp.advanced:
                params[sp.name] = render_param(container, scheme.name, sp, scheme)
        if show_advanced:
            for sp in advanced_by_group.pop(group, []):
                params[sp.name] = render_param(container, scheme.name, sp, scheme)

    if advanced:
        def _toggle_advanced(k=advanced_key):
            st.session_state[k] = not st.session_state.get(k, False)

        container.button(
            "Hide advanced" if show_advanced else f"Show advanced · {len(advanced)}",
            key=skey(scheme.name, "_advanced_toggle"), on_click=_toggle_advanced,
            type="tertiary", icon=":material/tune:")
        if show_advanced:
            for group, group_specs in advanced_by_group.items():
                if len(group_specs) > ADVANCED_FOLD:
                    box = container.expander(f"{group} · {len(group_specs)}")
                else:
                    render_group_header(container, group)
                    box = container
                for sp in group_specs:
                    params[sp.name] = render_param(box, scheme.name, sp, scheme)

    for sp in specs:
        if sp.name not in params:
            params[sp.name] = st.session_state[skey(scheme.name, sp.name)]
    return params
