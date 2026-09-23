"""Top bar pieces — the brand mark and the two segmented switches beside it.

The bar carries what a scheme is, not how it is tuned: the scheme switcher,
About, the regime/mode switch (or, for schemes without one, the recommended
default sets) and the two links out (docs/ui_redesign/decisions.md D8).
"""
import base64

import streamlit as st

from gabes_ui.controls import (
    apply_default_set,
    apply_recommended_defaults,
    current_params,
    skey,
)
from gabes_ui.readout import clean_choice_label


def brand_html(icon_svg, variant):
    mark = base64.b64encode(icon_svg.encode("utf-8")).decode("ascii")
    return (f"<div class='gabes-brand {variant}'>"
            f"<img src='data:image/svg+xml;base64,{mark}' alt=''>"
            "<b>GABES</b><small>beta</small></div>")


def render_regime_selector(scheme_name, sp, scheme_obj):
    """The scheme's regime/mode switch, shown in the top bar."""
    key = skey(scheme_name, sp.name)
    last_key = skey(scheme_name, f"_{sp.name}_last")
    options = list(sp.choices or ())
    labels = getattr(sp, "choice_labels", None) or {}
    if st.session_state.get(key) in options:
        st.session_state[last_key] = st.session_state[key]

    def _changed():
        # A segmented control can be clicked off; keep a regime selected.
        if st.session_state.get(key) not in options:
            st.session_state[key] = st.session_state.get(last_key, sp.default)
        if getattr(sp, "applies_defaults", False):
            apply_recommended_defaults(scheme_name, scheme_obj, key)

    st.segmented_control(
        sp.label, options, key=key, help=sp.help or None, on_change=_changed,
        format_func=lambda value: clean_choice_label(labels.get(value, value)),
        label_visibility="collapsed",
    )


def render_preset_selector(scheme_name, scheme_obj, sets):
    """Recommended default sets (e.g. OD / SAS) as a segmented control whose
    selection shows which set the current parameters match, if any."""
    key = skey(scheme_name, "_preset")
    labels = list(sets)
    current = current_params(scheme_name, scheme_obj)
    st.session_state[key] = next(
        (label for label in labels
         if all(current.get(k) == v for k, v in (sets[label] or {}).items())),
        None,
    )

    def _changed():
        picked = st.session_state.get(key)
        if picked is not None:
            apply_default_set(scheme_name, scheme_obj, picked)

    st.segmented_control(
        "Defaults", labels, key=key, on_change=_changed,
        format_func=lambda label: clean_choice_label(label).removesuffix(" default"),
        help="Load a ready-made parameter set.", label_visibility="collapsed",
    )
