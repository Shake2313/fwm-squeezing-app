"""
GABES front-end — generic, scheme-driven.

The UI knows nothing about any specific physics. It picks a Scheme from the
registry (gabes.schemes), renders exactly the controls that scheme declares
(param_schema), and draws whatever observables it returns. Adding a scheme adds
controls and plots automatically — no edits here.

Two-tier compute (preserved from the original FWM app): the heavy solve is
cached on the scheme's `recompute` knobs only, so navigate-only knobs (e.g. the
FWM two-photon detuning) update the readout instantly without re-solving.

Run with:
    streamlit run streamlit_app.py
"""
import matplotlib
matplotlib.use("Agg")          # headless server backend (no GUI / Tk)
import streamlit as st
from contextlib import nullcontext
from pathlib import Path

from gabes import schemes
from gabes.core import blas_single_thread

from gabes_ui import plotcard
from gabes_ui import theme as ui_theme
from gabes_ui.update_status import render_update_status
from gabes_ui.controls import (
    current_params, param_visible, render_derived_params, render_rail, skey,
)
from gabes_ui.defaults import render_defaults
from gabes_ui.dropdowns import render_dropdown_policy
from gabes_ui.guide import guide_button
from gabes_ui.plotcard import PLOT_LOCK, render_more, render_plot_section
from gabes_ui.readout import (
    metrics_table_markdown,
    partition_readout,
    present_metrics,
    scheme_label,
    strip_html,
)
from gabes_ui.shell import (
    brand_html,
    apply_scheme_link,
    render_preset_selector,
    render_regime_selector,
)

APP_DIR = Path(__file__).resolve().parent


@st.cache_data(show_spinner=False)
def _asset_text(filename):
    return (APP_DIR / "assets" / filename).read_text(encoding="utf-8")


THEME_BASE = (st.get_option("theme.base") or "light").lower()
TOKENS = ui_theme.tokens(THEME_BASE)
ICON_ASSET = "gabes-mark-v3-dark.svg" if THEME_BASE == "dark" else "gabes-mark-v3.svg"

# User's Guide, served as a static file (config.toml -> server.enableStaticServing).
# Streamlit exposes ./static/<f> at the relative URL "app/static/<f>", which
# resolves correctly both locally and on Streamlit Community Cloud — so the link
# opens from any computer. The file is fully self-contained (images base64-inlined
# by docs/Userguide/build_static_guide.py), so it needs no sibling assets.
_guide_file = APP_DIR / "static" / "GABES_User_Guide.html"
_guide_revision = _guide_file.stat().st_mtime_ns if _guide_file.exists() else 0
GUIDE_URL = f"app/static/GABES_User_Guide.html?v={_guide_revision}"

ICON_SVG = _asset_text(ICON_ASSET)

st.set_page_config(page_title="GABES — Atomic Bloch Equation Solver",
                   page_icon=ICON_SVG, layout="wide")


READOUT_CACHE_VERSION = "v17-rydberg-workspace-20261105"


def _inject_css():
    # Tokens + stylesheet live in gabes_ui/theme.py and assets/ui/gabes.css.
    # The anchor lets the stylesheet hide its own element slot (no layout gap).
    st.markdown(f"<style>{ui_theme.app_css(THEME_BASE)}</style>"
                "<span class='gabes-style-anchor'></span>", unsafe_allow_html=True)


_inject_css()
render_dropdown_policy()


# ----------------------------------------------------------------------
# Cached compute layer (keyed on the scheme + its recompute knobs only)
# ----------------------------------------------------------------------
@st.cache_data(show_spinner=False, max_entries=64)
def _cached_compute(scheme_name, recompute_items, cache_version):
    with blas_single_thread():
        return schemes.get(scheme_name).compute(dict(recompute_items))


@st.cache_data(show_spinner=False, max_entries=64)
def _cached_observables(scheme_name, raw, param_items, cache_version):
    # Matplotlib's font/mathtext/layout caches are process-global. Streamlit can
    # briefly overlap reruns when sliders are moved quickly, so serialize figure
    # construction to avoid layout-time parser crashes.
    with PLOT_LOCK:
        return schemes.get(scheme_name).observables(raw, dict(param_items))


@st.cache_data(show_spinner=False, max_entries=8)
def _cached_model_comparison(scheme_name, param_items, cache_version, readout_version,
                             _current_raw, _current_view):
    """Cache small sweeps by committed inputs; reuse the active single-condition solve."""
    active_scheme = schemes.get(scheme_name)
    solve_keys = active_scheme.recompute_keys()

    def resolve_raw(condition):
        items = tuple(sorted((key, condition[key]) for key in solve_keys))
        return _cached_compute(scheme_name, items, cache_version)

    return active_scheme.model_comparison(dict(param_items), resolve_raw=resolve_raw,
                                          current_raw=_current_raw, current_view=_current_view)


# ----------------------------------------------------------------------
# Router — SABES is a separate site living in the same deployment
# ----------------------------------------------------------------------
# SABES is not a scheme: it is a lab-facing layer over one specific experiment,
# with primary variables (waveplate angles, generator frequency, lens choices)
# instead of processed physics parameters. Routing on a query parameter keeps it
# out of the scheme dropdown and off this page entirely, while still deploying as
# one app. What it needs from the host is handed over explicitly rather than
# imported back, so `sabes_page` never depends on this module's import state.
SABES_QUERY_VALUE = "sabes"


def _sabes_host():
    from types import SimpleNamespace
    return SimpleNamespace(render_fig=plotcard.render_fig, theme_base=THEME_BASE)


def _open_sabes():
    st.query_params["app"] = SABES_QUERY_VALUE


if st.query_params.get("app") == SABES_QUERY_VALUE:
    import sabes_page
    sabes_page.render(host=_sabes_host())
    st.stop()


# ----------------------------------------------------------------------
# Top bar — brand, scheme, about, regime, guide, SABES
# ----------------------------------------------------------------------
# The rail's brand band and the top bar share one 52 px line; the top-bar copy
# of the brand shows only while the sidebar is collapsed (phone, or folded).
all_schemes = schemes.all_schemes()
titles = [s.title for s in all_schemes]
apply_scheme_link(st.session_state, st.query_params, all_schemes)
_scheme_title_migrations = {
    "Four-wave mixing (Gain diagnostic / Biphoton)":
        "Four-wave mixing (Squeezing / Biphoton)",
}
saved_scheme_title = st.session_state.get("_scheme_choice")
if saved_scheme_title in _scheme_title_migrations:
    st.session_state["_scheme_choice"] = _scheme_title_migrations[saved_scheme_title]
elif saved_scheme_title is not None and saved_scheme_title not in titles:
    del st.session_state["_scheme_choice"]

topbar = st.container(horizontal=True, vertical_alignment="center", gap="small",
                      key="gabes_topbar")
with topbar:
    st.markdown(brand_html(ICON_SVG, "gabes-brand--top"), unsafe_allow_html=True, width="content")
    def _scheme_changed():
        selected = all_schemes[titles.index(st.session_state["_scheme_choice"])]
        st.query_params["scheme"] = selected.name
        st.query_params.pop("compare", None)
        st.session_state["_scheme_link_applied"] = selected.name

    choice = st.selectbox("Scheme", titles, key="_scheme_choice", format_func=scheme_label,
                          on_change=_scheme_changed,
                          label_visibility="collapsed", width=240)
scheme = all_schemes[titles.index(choice)]

specs = scheme.param_schema()
defaults_version = getattr(scheme, "defaults_version", "1")
defaults_key = skey(scheme.name, "_defaults_version")
if st.session_state.get(defaults_key) != defaults_version:
    for sp in specs:
        if sp.control == "derived_slider":
            continue
        st.session_state[skey(scheme.name, sp.name)] = sp.default
    st.session_state[defaults_key] = defaults_version
else:
    for sp in specs:
        if sp.control == "derived_slider":
            continue
        key = skey(scheme.name, sp.name)
        # Preserve hidden controls through widget cleanup and send their saved
        # values to the browser when they become visible again.
        st.session_state[key] = st.session_state.get(key, sp.default)

# Reserve the main instrument above saved settings, but mount the defaults
# store before any editable parameter widgets so preset callbacks stay atomic.
is_sas = scheme.name == "sas"
has_spectra = scheme.model_comparison_parameter is not None
if has_spectra:
    workspace = st.container(key="gabes_sas_workspace")
    render_defaults(st.expander("Saved settings", icon=":material/bookmarks:"), scheme)
else:
    workspace = nullcontext()
    st.sidebar.markdown(brand_html(ICON_SVG, "gabes-brand--rail"), unsafe_allow_html=True)
    render_defaults(st.sidebar, scheme)

# The regime/mode switch moves to the top bar: the first segmented control that
# applies a recommended default set (FWM Mode, Λ/Rydberg/magneto Regime).
regime_spec = next(
    (sp for sp in specs
     if getattr(sp, "applies_defaults", False)
     and getattr(sp, "control", "auto") == "segmented"
     and param_visible(scheme.name, sp)),
    None,
)
# Schemes without such a switch offer their recommended sets (e.g. OD / SAS)
# there instead. Probed defensively so a scheme without the hook never breaks.
_rec_fn = getattr(scheme, "recommended_defaults", None)
_rec_sets = None
if callable(_rec_fn):
    try:
        _rec_sets = _rec_fn(current_params(scheme.name, scheme))
    except Exception:
        _rec_sets = None
_mode_driven_defaults = any(getattr(sp, "applies_defaults", False) for sp in specs)

with topbar:
    about = st.popover("About", icon=":material/info:", type="tertiary",
                       help="What this scheme models, and its references.")
    if regime_spec is not None:
        render_regime_selector(scheme.name, regime_spec, scheme)
    elif isinstance(_rec_sets, dict) and _rec_sets and not _mode_driven_defaults:
        render_preset_selector(scheme.name, scheme, _rec_sets)
    with st.container(horizontal=True, horizontal_alignment="right",
                      vertical_alignment="center", gap="small", key="gabes_topbar_end"):
        guide_button(GUIDE_URL)
        st.button("SABES", icon=":material/north_east:", icon_position="right",
                  type="tertiary", on_click=_open_sabes,
                  help="A separate, lab-facing simulator for one experiment: the "
                       "EOM-based 85Rb twin-beam squeezing setup, driven from "
                       "primary settings (waveplate angles, generator frequency, "
                       "lens choices).")
with about:
    st.caption(f"Cluster {scheme.cluster}")
    st.markdown(scheme.caption)
    info = scheme.info()
    if info:
        st.divider()
        st.markdown(info)


def render_spectra():
    comparing = has_spectra and st.query_params.get("compare") != "0"
    derived_slots = {}
    if has_spectra:
        from gabes_ui.conditions import render_conditions, render_power_inputs

        view_key = skey(scheme.name, "_spectra_view")
        model_label = "Model spectra"
        csv_label = "Experimental CSV (legacy)" if is_sas else "Single-condition details"
        st.session_state[view_key] = model_label if comparing else csv_label

        def change_view():
            picked = st.session_state[view_key]
            if picked is None:
                return
            entering = picked == model_label
            if is_sas:
                plotcard.remember_overlay_for_comparison(scheme.name, entering=entering)
            if entering:
                st.query_params.pop("compare", None)
            else:
                st.query_params["compare"] = "0"

        with st.container(horizontal=True, vertical_alignment="center", gap="small",
                          key="gabes_comparison_actions"):
            st.segmented_control("Spectrum view", [model_label, csv_label], key=view_key,
                                 on_change=change_view, label_visibility="collapsed")
            if is_sas:
                guide_button(GUIDE_URL, key="sas_first_spectrum", label="First spectrum",
                             anchor="first-spectrum")
            render_update_status(scheme.name)
        params = render_conditions(scheme, specs)
        if not comparing:
            params = render_power_inputs(scheme, specs)
    else:
        from gabes_ui.conditions import render_atom_conditions
        with st.container(key="gabes_conditions"):
            st.html('<h3 class="gabes-conditions-heading">Conditions</h3>')
            grouped = render_atom_conditions(scheme, specs)
        params = render_rail(st.sidebar, scheme, specs,
                             skip=(regime_spec, *(sp for sp in specs if sp.name in grouped)),
                             derived_slots=derived_slots)
        render_update_status(scheme.name)

    recompute_items = tuple(sorted((k, params[k]) for k in scheme.recompute_keys()))
    cache_version = scheme.cache_version
    from gabes_ui.spectral_sweep import cancel_inactive_sweeps, reusable_solution
    cancel_inactive_sweeps(st.session_state, scheme.name)
    retained = reusable_solution(st.session_state.get("spectral-workspace-"+scheme.name),
                                 scheme, params, READOUT_CACHE_VERSION) if comparing else None
    raw = retained[0] if retained is not None else _cached_compute(scheme.name, recompute_items, cache_version)
    param_items = tuple(sorted(params.items()))
    if comparing:
        view = retained[1] if retained is not None and retained[1] is not None else scheme.headless_observables(raw, params)
    elif scheme.cache_observables:
        view = _cached_observables(scheme.name, raw, param_items, (cache_version, READOUT_CACHE_VERSION))
    else:
        with PLOT_LOCK:
            view = scheme.observables(raw, params)
    render_derived_params(derived_slots, view.get("sidebar_readouts", {}))

    if comparing:
        from gabes_ui.model_compare_section import comparison_context, render_comparison

        from gabes_ui.spectral_sweep import reusable_dataset
        comparison = reusable_dataset(st.session_state.get("spectral-workspace-"+scheme.name),
                                      scheme, params, READOUT_CACHE_VERSION)
        if comparison is None:
            with st.spinner("Preparing spectra…"):
                comparison = _cached_model_comparison(scheme.name, param_items, cache_version,
                                                      READOUT_CACHE_VERSION, raw, view)
        if comparison["schema"] == "gabes-spectral-comparison-v1":
            from gabes_ui.spectral_compare import render_spectral_section
            def resolve_condition(condition):
                items = tuple(sorted((k, condition[k]) for k in scheme.recompute_keys()))
                return _cached_compute(scheme.name, items, cache_version)
            render_spectral_section(comparison, scheme=scheme, raw=raw, params=params,
                                    view=view, readout_version=READOUT_CACHE_VERSION,
                                    resolve_raw=resolve_condition)
        else:
            render_comparison(comparison, comparison_context(scheme, params, comparison["schema"]),
                              scheme=scheme, raw=raw, params=params, view=view,
                              readout_version=READOUT_CACHE_VERSION)
        return

    spec_by_name = {sp.name: sp for sp in specs}
    if is_sas:
        render_plot_section(dict(view, experimental_workspace=True), scheme, spec_by_name, raw, params,
                            cache_version, READOUT_CACHE_VERSION)
        return

    metrics = present_metrics(view.get("metrics", []), scheme_name=scheme.name, raw=raw)
    if metrics:
        heroes, secondary, statuses, overflow = partition_readout(
            metrics, hero_count=view.get("hero_count") or 2)
        with st.container(horizontal=True, vertical_alignment="center", gap="small",
                          key="gabes_readout"):
            st.markdown(strip_html(heroes, secondary, statuses,
                                   compact_statuses=scheme.name == "magneto"), unsafe_allow_html=True)
            if overflow:
                with st.popover(f"+{len(overflow)} more", type="tertiary"):
                    st.markdown(metrics_table_markdown(metrics))

    render_plot_section(view, scheme, spec_by_name, raw, params, cache_version, READOUT_CACHE_VERSION)
    render_more(view, scheme, params, cache_version)


with workspace:
    render_spectra()
