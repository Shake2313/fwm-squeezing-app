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
from pathlib import Path

from gabes import schemes
from gabes.core import blas_single_thread

from gabes_ui import export as ui_export
from gabes_ui import plotcard
from gabes_ui import theme as ui_theme
from gabes_ui.controls import current_params, param_visible, render_rail, skey
from gabes_ui.guide import guide_button
from gabes_ui.plotcard import PLOT_LOCK, render_more, render_plot_card
from gabes_ui.readout import (
    metrics_table_markdown,
    partition_readout,
    strip_html,
)
from gabes_ui.shell import (
    brand_html,
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
GUIDE_URL = "app/static/GABES_User_Guide.html"

ICON_SVG = _asset_text(ICON_ASSET)

st.set_page_config(page_title="GABES — Atomic Bloch Equation Solver",
                   page_icon=ICON_SVG, layout="wide")


READOUT_CACHE_VERSION = "hero-ribbon-v3-single-hero"


def _inject_css():
    # Tokens + stylesheet live in gabes_ui/theme.py and assets/ui/gabes.css.
    # The anchor lets the stylesheet hide its own element slot (no layout gap).
    st.markdown(f"<style>{ui_theme.app_css(THEME_BASE)}</style>"
                "<span class='gabes-style-anchor'></span>", unsafe_allow_html=True)


_inject_css()


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
st.sidebar.markdown(brand_html(ICON_SVG, "gabes-brand--rail"), unsafe_allow_html=True)

all_schemes = schemes.all_schemes()
titles = [s.title for s in all_schemes]
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
    choice = st.selectbox("Scheme", titles, key="_scheme_choice",
                          label_visibility="collapsed", width=310)
scheme = all_schemes[titles.index(choice)]

specs = scheme.param_schema()
defaults_version = getattr(scheme, "defaults_version", "1")
defaults_key = skey(scheme.name, "_defaults_version")
if st.session_state.get(defaults_key) != defaults_version:
    for sp in specs:
        st.session_state[skey(scheme.name, sp.name)] = sp.default
    st.session_state[defaults_key] = defaults_version
else:
    for sp in specs:
        key = skey(scheme.name, sp.name)
        # Preserve hidden controls through widget cleanup and send their saved
        # values to the browser when they become visible again.
        st.session_state[key] = st.session_state.get(key, sp.default)

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


# ----------------------------------------------------------------------
# Control rail
# ----------------------------------------------------------------------
params = render_rail(st.sidebar, scheme, specs, skip=(regime_spec,))


# ----------------------------------------------------------------------
# Compute (cached) + observables
# ----------------------------------------------------------------------
recompute_items = tuple(sorted((k, params[k]) for k in scheme.recompute_keys()))
cache_version = getattr(scheme, "cache_version", "1")
with st.spinner("Solving Bloch equations…"):
    raw = _cached_compute(scheme.name, recompute_items, cache_version)
param_items = tuple(sorted(params.items()))
if getattr(scheme, "cache_observables", False):
    view = _cached_observables(
        scheme.name, raw, param_items,
        (cache_version, READOUT_CACHE_VERSION),
    )
else:
    with PLOT_LOCK:
        view = scheme.observables(raw, params)


# ----------------------------------------------------------------------
# Readout strip → plot card → More
# ----------------------------------------------------------------------
metrics = view.get("metrics", [])
if metrics:
    heroes, secondary, statuses, overflow = partition_readout(
        metrics, hero_count=view.get("hero_count") or 2)
    with st.container(horizontal=True, vertical_alignment="center", gap="small",
                      key="gabes_readout"):
        st.markdown(strip_html(heroes, secondary, statuses), unsafe_allow_html=True)
        if overflow:
            with st.popover(f"+{len(overflow)} more", type="tertiary"):
                st.markdown(metrics_table_markdown(metrics))

spec_by_name = {sp.name: sp for sp in specs}
comparison_payload, plot_payloads, export_box = render_plot_card(view, scheme, spec_by_name)

render_more(view, scheme, params, cache_version)

if export_box is not None:
    ui_export.render_export(
        export_box,
        scheme=scheme,
        raw=raw,
        params=params,
        view=view,
        plot_payloads=plot_payloads,
        comparison_payload=comparison_payload,
        cache_version=cache_version,
        readout_cache_version=READOUT_CACHE_VERSION,
    )
