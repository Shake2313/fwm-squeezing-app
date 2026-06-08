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
from html import escape

from gabes import schemes
from gabes.core import blas_single_thread
from gabes.plot_style import apply_gabes_plot_style

APP_DIR = Path(__file__).resolve().parent


@st.cache_data(show_spinner=False)
def _asset_text(filename):
    return (APP_DIR / "assets" / filename).read_text(encoding="utf-8")


THEME_BASE = (st.get_option("theme.base") or "light").lower()
LOGO_ASSET = "gabes-logo-v3-dark.svg" if THEME_BASE == "dark" else "gabes-logo-v3.svg"
ICON_ASSET = "gabes-mark-v3-dark.svg" if THEME_BASE == "dark" else "gabes-mark-v3.svg"

LOGO_SVG = _asset_text(LOGO_ASSET)
ICON_SVG = _asset_text(ICON_ASSET)
SIDEBAR_LOGO_SVG = LOGO_SVG.replace(
    'width="960" height="232" viewBox="0 0 960 232"',
    'width="738" height="232" viewBox="0 0 738 232"',
    1,
)

st.set_page_config(page_title="GABES — Atomic Bloch Equation Solver",
                   page_icon=ICON_SVG, layout="wide")


GROUP_STYLES = {
    "model": dict(label="model", color="#4F46E5", bg="#EEF2FF"),
    "view": dict(label="view", color="#4F46E5", bg="#EEF2FF"),
    "atomic": dict(label="atomic", color="#475569", bg="#F1F5F9"),
    "fields": dict(label="fields", color="#0284C7", bg="#E6F7FC"),
    "pump": dict(label="field", color="#0284C7", bg="#E6F7FC"),
    "detunings": dict(label="detuning", color="#F43F5E", bg="#FFF1F3"),
    "cell & beams": dict(label="cell", color="#0F766E", bg="#EAF8F5"),
    "detection & scaling": dict(label="readout", color="#7C3AED", bg="#F3E8FF"),
    "numerics": dict(label="numerics", color="#64748B", bg="#F1F5F9"),
    "default": dict(label="preset", color="#0284C7", bg="#E6F7FC"),
}

METRIC_STYLES = [
    (("gain",), dict(kind="gain", color="#F97316", bg="#FFF7ED")),
    (("squeez", "coinc", "nonclass", "pair"), dict(kind="quantum", color="#F43F5E", bg="#FFF1F3")),
    (("od", "trans", "absorp"), dict(kind="optical", color="#0284C7", bg="#E6F7FC")),
    (("rotation", "nmor", "larmor", "mag"), dict(kind="magneto", color="#7C3AED", bg="#F3E8FF")),
    (("density", "temp", "cell", "n("), dict(kind="cell", color="#0F766E", bg="#EAF8F5")),
    (("status",), dict(kind="status", color="#64748B", bg="#F1F5F9")),
]

DEFAULT_STYLE = dict(label="control", color="#0284C7", bg="#E6F7FC")
DEFAULT_METRIC_STYLE = dict(kind="result", color="#2563EB", bg="#EFF6FF")


def _inject_css():
    st.markdown("""
<style>
:root {
  --gabes-bg: #F6F8FB;
  --gabes-surface: #FFFFFF;
  --gabes-ink: #0F172A;
  --gabes-muted: #64748B;
  --gabes-border: #DCE6EF;
  --gabes-grid: #E6EDF5;
  --gabes-primary: #0284C7;
  --gabes-rose: #F43F5E;
}

[data-testid="stAppViewContainer"] {
  background: var(--gabes-bg);
}

[data-testid="stHeader"] {
  height: 0 !important;
  min-height: 0 !important;
  background: transparent !important;
}

[data-testid="stHeader"] > div,
[data-testid="stToolbar"],
[data-testid="stDecoration"],
[data-testid="stStatusWidget"],
.stDeployButton {
  display: none !important;
}

[data-testid="stMainBlockContainer"],
[data-testid="stAppViewBlockContainer"],
section.main > div,
.main .block-container,
.block-container {
  padding-top: 0.85rem !important;
  padding-left: clamp(0.75rem, 1vw, 1rem) !important;
  padding-right: clamp(0.75rem, 1vw, 1rem) !important;
  padding-bottom: 1.35rem !important;
  max-width: none !important;
}

[data-testid="stSidebar"] {
  background: #FBFDFF;
  border-right: 1px solid var(--gabes-border);
}

[data-testid="stSidebarHeader"] {
  display: flex !important;
  height: 0 !important;
  min-height: 0 !important;
  margin: 0 !important;
  padding: 0 !important;
  overflow: visible !important;
  position: relative;
  z-index: 30;
}

[data-testid="stSidebarHeader"] button[data-testid="stBaseButton-headerNoPadding"] {
  position: absolute;
  top: 0.55rem;
  right: -1.25rem;
  width: 2rem !important;
  height: 2rem !important;
  min-width: 2rem !important;
  padding: 0 !important;
  visibility: visible !important;
  opacity: 1 !important;
  pointer-events: auto !important;
  border-radius: 8px !important;
  background: rgba(251, 253, 255, 0.88) !important;
  border: 1px solid rgba(220, 230, 239, 0.9) !important;
  box-shadow: 0 1px 3px rgba(15, 23, 42, 0.08);
}

[data-testid="stSidebar"][aria-expanded="false"] [data-testid="stSidebarContent"] {
  overflow: visible !important;
}

[data-testid="stSidebar"][aria-expanded="false"] [data-testid="stSidebarHeader"] button[data-testid="stBaseButton-headerNoPadding"] {
  transform: translateX(312px);
}

[data-testid="stSidebar"][aria-expanded="false"] [data-testid="stSidebarHeader"] [data-testid="stIconMaterial"] {
  font-size: 0 !important;
}

[data-testid="stSidebar"][aria-expanded="false"] [data-testid="stSidebarHeader"] [data-testid="stIconMaterial"]::before {
  content: "keyboard_double_arrow_right";
  font-family: inherit;
  font-size: 1.25rem;
  line-height: 1;
}

[data-testid="stSidebarContent"],
[data-testid="stSidebarUserContent"] {
  padding-top: 0.85rem !important;
}

[data-testid="stSidebar"] [data-testid="stImage"] {
  margin: 0 0 0.45rem;
}

[data-testid="stSidebar"] [data-testid="stCaptionContainer"] {
  color: var(--gabes-muted);
}

[data-testid="stSidebar"] hr {
  margin: 0.75rem 0 1rem;
  border-color: var(--gabes-border);
}

h1, h2, h3, h4, h5, h6 {
  color: var(--gabes-ink);
  letter-spacing: 0;
}

p, li, label, [data-testid="stMarkdownContainer"] {
  color: var(--gabes-ink);
}

div[data-testid="stButton"] > button {
  border-radius: 8px;
  border: 1px solid var(--gabes-border);
  background: var(--gabes-surface);
  color: var(--gabes-ink);
  box-shadow: 0 1px 2px rgba(15, 23, 42, 0.04);
}

div[data-testid="stButton"] > button:hover {
  border-color: var(--gabes-primary);
  color: var(--gabes-primary);
}

[data-testid="stSelectbox"] div[data-baseweb="select"] > div,
[data-testid="stSlider"] {
  color: var(--gabes-ink);
}

[data-testid="stExpander"] {
  border: 1px solid var(--gabes-border);
  border-radius: 8px;
  background: var(--gabes-surface);
}

.gabes-header {
  margin: 0 0 0.72rem;
  padding: 0.68rem 0.78rem 0.78rem;
  background: var(--gabes-surface);
  border: 1px solid var(--gabes-border);
  border-radius: 8px;
  box-shadow: 0 4px 14px rgba(15, 23, 42, 0.03);
}

.gabes-hairline {
  height: 2px;
  width: 100%;
  border-radius: 999px;
  margin-bottom: 0.55rem;
  background: linear-gradient(90deg, #22D3EE 0%, #2563EB 52%, #F43F5E 100%);
}

.gabes-header-row {
  display: flex;
  flex-wrap: wrap;
  align-items: center;
  gap: 0.34rem;
  margin-bottom: 0.34rem;
}

.gabes-badge {
  display: inline-flex;
  align-items: center;
  min-height: 1.2rem;
  padding: 0.1rem 0.4rem;
  border-radius: 999px;
  border: 1px solid var(--gabes-border);
  background: #F8FAFC;
  color: var(--gabes-muted);
  font-size: 0.68rem;
  font-weight: 650;
  line-height: 1;
}

.gabes-header h1 {
  margin: 0;
  color: var(--gabes-ink);
  font-size: clamp(1.48rem, 2.25vw, 2.12rem);
  line-height: 1.12;
  font-weight: 760;
  letter-spacing: 0;
}

.gabes-header p {
  max-width: 68rem;
  margin: 0.34rem 0 0;
  color: var(--gabes-muted);
  font-size: 0.91rem;
  line-height: 1.42;
}

.gabes-group-header {
  display: flex;
  align-items: center;
  gap: 0.45rem;
  margin: 1rem 0 0.35rem;
  padding: 0.2rem 0 0.25rem;
  border-bottom: 1px solid var(--gabes-border);
  color: var(--gabes-ink);
  font-size: 0.83rem;
  font-weight: 760;
}

.gabes-group-dot {
  width: 0.5rem;
  height: 0.5rem;
  border-radius: 999px;
  flex: 0 0 auto;
}

.gabes-endpoints {
  display: flex;
  justify-content: space-between;
  margin-top: -0.72rem;
  color: var(--gabes-muted);
  font-size: 0.78rem;
}

.gabes-metric-grid {
  display: grid;
  grid-template-columns: repeat(auto-fit, minmax(136px, 1fr));
  gap: 0.32rem;
}

.gabes-metric-card {
  min-height: 3.35rem;
  padding: 0.38rem 0.48rem 0.4rem;
  background: var(--gabes-surface);
  border: 1px solid var(--gabes-border);
  border-left: 3px solid var(--metric-color);
  border-radius: 8px;
  box-shadow: 0 2px 9px rgba(15, 23, 42, 0.025);
}

.gabes-metric-label {
  margin-top: 0;
  color: var(--gabes-muted);
  font-size: 0.7rem;
  line-height: 1.25;
  font-weight: 650;
  white-space: nowrap;
  overflow: hidden;
  text-overflow: ellipsis;
}

.gabes-metric-value {
  margin-top: 0.04rem;
  color: var(--gabes-ink);
  font-size: clamp(0.95rem, 1.22vw, 1.12rem);
  line-height: 1.15;
  font-weight: 760;
}

.gabes-metric-delta {
  margin-top: 0.06rem;
  color: var(--gabes-muted);
  font-size: 0.68rem;
}

.gabes-section-gap {
  height: 0.42rem;
}

.gabes-plot-gap {
  height: 0.56rem;
}

.stDataFrame, [data-testid="stTable"] {
  border-radius: 8px;
}
</style>
""", unsafe_allow_html=True)


_inject_css()


# ----------------------------------------------------------------------
# Cached compute layer (keyed on the scheme + its recompute knobs only)
# ----------------------------------------------------------------------
@st.cache_data(show_spinner=False)
def _cached_compute(scheme_name, recompute_items, cache_version):
    with blas_single_thread():
        return schemes.get(scheme_name).compute(dict(recompute_items))


@st.cache_data(show_spinner=False)
def _cached_extra(scheme_name, view_key, recompute_items, cache_version):
    scheme = schemes.get(scheme_name)
    view = next(v for v in scheme.extra_views() if v.key == view_key)
    with blas_single_thread():
        return view.compute(dict(recompute_items))


@st.cache_data(show_spinner=False, max_entries=64)
def _cached_observables(scheme_name, raw, param_items, cache_version):
    return schemes.get(scheme_name).observables(raw, dict(param_items))


def _close_fig(fig):
    import matplotlib.pyplot as plt
    plt.close(fig)


def _skey(scheme_name, pname):
    return f"{scheme_name}__{pname}"


def _concept_style(group):
    name = group.lower()
    if name in GROUP_STYLES:
        return GROUP_STYLES[name]
    for key, style in GROUP_STYLES.items():
        if key in name:
            return style
    return DEFAULT_STYLE


def _metric_style(label):
    name = label.lower()
    for keywords, style in METRIC_STYLES:
        if any(key in name for key in keywords):
            return style
    return DEFAULT_METRIC_STYLE


def _render_group_header(container, group):
    style = _concept_style(group)
    container.markdown(
        "<div class='gabes-group-header'>"
        f"<span class='gabes-group-dot' style='background:{style['color']}'></span>"
        f"<span>{escape(group)}</span>"
        "</div>",
        unsafe_allow_html=True,
    )


def _apply_recommended_defaults(scheme_name, scheme_obj, key):
    """on_change for a control whose selection applies the matching
    recommended_defaults set (so choosing it also resets that mode's knobs)."""
    selection = st.session_state.get(key)
    cur = {
        sp.name: st.session_state.get(_skey(scheme_name, sp.name), sp.default)
        for sp in scheme_obj.param_schema()
    }
    sets = scheme_obj.recommended_defaults(cur) or {}
    for k, v in (sets.get(selection) or {}).items():
        st.session_state[_skey(scheme_name, k)] = v


def _render_param(container, scheme_name, sp, scheme_obj=None):
    key = _skey(scheme_name, sp.name)
    label = sp.label + (f"  [{sp.unit}]" if sp.unit else "")
    help_ = sp.help or None
    has_state = key in st.session_state
    current = st.session_state.get(key, sp.default)
    on_change = None
    if getattr(sp, "applies_defaults", False) and scheme_obj is not None:
        on_change = lambda: _apply_recommended_defaults(scheme_name, scheme_obj, key)
    if getattr(sp, "control", "auto") == "segmented":
        options = list(sp.choices or ())
        if hasattr(container, "segmented_control"):
            try:
                return container.segmented_control(label, options, key=key,
                                                   help=help_, on_change=on_change)
            except TypeError:
                pass
        if has_state:
            return container.radio(label, options, key=key, help=help_,
                                   horizontal=True, on_change=on_change)
        idx = options.index(current) if current in options else 0
        return container.radio(label, options, key=key, help=help_,
                               horizontal=True, index=idx, on_change=on_change)
    if sp.choices is not None:
        options = list(sp.choices)
        if has_state:
            return container.selectbox(label, options, key=key, help=help_,
                                       on_change=on_change)
        idx = options.index(current) if current in options else 0
        return container.selectbox(label, options, index=idx, key=key, help=help_,
                                   on_change=on_change)
    if has_state:
        val = container.slider(label, sp.vmin, sp.vmax, step=sp.step,
                               key=key, help=help_)
    else:
        val = container.slider(label, sp.vmin, sp.vmax, value=current,
                               step=sp.step, key=key, help=help_)
    endpoints = getattr(sp, "endpoints", None)
    if endpoints:
        left, right = endpoints
        container.markdown(
            "<div class='gabes-endpoints'>"
            f"<span>{escape(str(left))}</span><span>{escape(str(right))}</span></div>",
            unsafe_allow_html=True)
    return val


def _param_visible(scheme_name, sp):
    if getattr(sp, "hidden", False):
        return False
    cond = getattr(sp, "visible_if", None)
    if not cond:
        return True
    for pname, allowed in cond.items():
        cur = st.session_state.get(_skey(scheme_name, pname))
        if isinstance(allowed, (set, tuple, list)):
            if cur not in allowed:
                return False
        elif cur != allowed:
            return False
    return True


def _render_scheme_header(scheme):
    recompute_count = len(scheme.recompute_keys())
    st.markdown(
        "<section class='gabes-header'>"
        "<div class='gabes-hairline'></div>"
        "<div class='gabes-header-row'>"
        f"<span class='gabes-badge'>{escape(scheme.cluster)}</span>"
        f"<span class='gabes-badge'>{recompute_count} solve knobs</span>"
        "<span class='gabes-badge'>scheme-driven</span>"
        "</div>"
        f"<h1>{escape(scheme.title)}</h1>"
        f"<p>{escape(scheme.caption)}</p>"
        "</section>",
        unsafe_allow_html=True,
    )


def _metric_card_html(metric):
    label = str(metric.get("label", ""))
    value = str(metric.get("value", ""))
    delta = metric.get("delta")
    help_text = metric.get("help") or ""
    style = _metric_style(label)
    delta_html = ""
    if delta is not None:
        delta_html = f"<div class='gabes-metric-delta'>{escape(str(delta))}</div>"
    return (
        "<div class='gabes-metric-card' "
        f"title='{escape(str(help_text))}' "
        f"style='--metric-color:{style['color']};--metric-bg:{style['bg']}'>"
        f"<div class='gabes-metric-label'>{escape(label)}</div>"
        f"<div class='gabes-metric-value'>{escape(value)}</div>"
        f"{delta_html}"
        "</div>"
    )


def _render_metrics(metrics):
    cards = "".join(_metric_card_html(metric) for metric in metrics)
    st.markdown(f"<div class='gabes-metric-grid'>{cards}</div>", unsafe_allow_html=True)


# ----------------------------------------------------------------------
# Sidebar — scheme selection
# ----------------------------------------------------------------------
st.sidebar.image(SIDEBAR_LOGO_SVG, width=230)

all_schemes = schemes.all_schemes()
titles = [s.title for s in all_schemes]
choice = st.sidebar.selectbox("Scheme", titles, key="_scheme_choice",
                              help="Pick the experiment / physics to model.")
scheme = all_schemes[titles.index(choice)]
st.sidebar.caption(f"Cluster {scheme.cluster}")
st.sidebar.divider()

specs = scheme.param_schema()
defaults_version = getattr(scheme, "defaults_version", "1")
defaults_key = _skey(scheme.name, "_defaults_version")
if st.session_state.get(defaults_key) != defaults_version:
    for sp in specs:
        st.session_state[_skey(scheme.name, sp.name)] = sp.default
    st.session_state[defaults_key] = defaults_version
else:
    for sp in specs:                               # seed defaults once
        st.session_state.setdefault(_skey(scheme.name, sp.name), sp.default)

# Presets — one click overwrites the relevant sliders.
scheme_presets = scheme.presets()
if scheme_presets and getattr(scheme, "presets_group", None):
    _render_group_header(st.sidebar, scheme.presets_group)
for preset in scheme_presets:
    def _apply(p=preset, sname=scheme.name):
        for k, v in p.values.items():
            st.session_state[_skey(sname, k)] = v
    st.sidebar.button(f"{preset.icon} {preset.name}", on_click=_apply,
                      use_container_width=True, help=preset.help)

# Context-aware default buttons — one per labelled preset the scheme offers for
# the current selection (e.g. "OD default" / "SAS default"). Probed defensively so
# a scheme without the optional hook never breaks the app.
_rec_fn = getattr(scheme, "recommended_defaults", None)
_rec_sets = None
if callable(_rec_fn):
    try:
        # Probe with the live selection (not static defaults) so a scheme can
        # offer readout/mode-dependent default buttons (e.g. magneto shows
        # transmission regimes vs an NMOR default). Keys are stable for SAS/FWM.
        _live_params = {
            sp.name: st.session_state.get(_skey(scheme.name, sp.name), sp.default)
            for sp in scheme.param_schema()
        }
        _rec_sets = _rec_fn(_live_params)
    except Exception:
        _rec_sets = None
# A control flagged applies_defaults (e.g. FWM's Mode) already applies these sets
# on selection, so the standalone "Default" buttons would just duplicate it.
_mode_driven_defaults = any(getattr(sp, "applies_defaults", False) for sp in specs)
if isinstance(_rec_sets, dict) and _rec_sets and not _mode_driven_defaults:
    _render_group_header(st.sidebar, "Default")
    _cols = st.sidebar.columns(len(_rec_sets))
    for _col, _label in zip(_cols, _rec_sets):
        def _apply_default(sname=scheme.name, sc=scheme, lbl=_label):
            cur = {
                sp.name: st.session_state.get(_skey(sname, sp.name), sp.default)
                for sp in sc.param_schema()
            }
            sets = sc.recommended_defaults(cur) or {}
            for k, v in (sets.get(lbl) or {}).items():
                st.session_state[_skey(sname, k)] = v
        _short = _label.replace(" default", "")
        _col.button(_short, on_click=_apply_default, use_container_width=True)

# Controls — grouped sections; advanced/numeric knobs fold into an expander.
visible_specs = [sp for sp in specs if _param_visible(scheme.name, sp)]
params = {}
group_order = []
for sp in visible_specs:
    if not sp.advanced and sp.group not in group_order:
        group_order.append(sp.group)

for g in group_order:
    _render_group_header(st.sidebar, g)
    for sp in visible_specs:
        if sp.group == g and not sp.advanced:
            params[sp.name] = _render_param(st.sidebar, scheme.name, sp, scheme)

advanced = [sp for sp in visible_specs if sp.advanced]
if advanced:
    exp = st.sidebar.expander("Advanced / numerics")
    for sp in advanced:
        params[sp.name] = _render_param(exp, scheme.name, sp, scheme)

for sp in specs:
    if sp.name not in params:
        params[sp.name] = st.session_state[_skey(scheme.name, sp.name)]


# ----------------------------------------------------------------------
# Compute (cached) + observables
# ----------------------------------------------------------------------
recompute_items = tuple(sorted((k, params[k]) for k in scheme.recompute_keys()))
cache_version = getattr(scheme, "cache_version", "1")
with st.spinner("Solving Bloch equations…"):
    raw = _cached_compute(scheme.name, recompute_items, cache_version)
param_items = tuple(sorted(params.items()))
if getattr(scheme, "cache_observables", False):
    view = _cached_observables(scheme.name, raw, param_items, cache_version)
else:
    view = scheme.observables(raw, params)


# ----------------------------------------------------------------------
# Header + readout
# ----------------------------------------------------------------------
_render_scheme_header(scheme)

metrics = view.get("metrics", [])
if metrics:
    _render_metrics(metrics)
    st.markdown("<div class='gabes-section-gap'></div>", unsafe_allow_html=True)

fig = view.get("figure")
if fig is not None:
    apply_gabes_plot_style(fig)
    st.pyplot(fig)
    _close_fig(fig)

for _title, _extra_fig in view.get("figures", []):
    st.markdown("<div class='gabes-plot-gap'></div>", unsafe_allow_html=True)
    st.subheader(_title)
    apply_gabes_plot_style(_extra_fig)
    st.pyplot(_extra_fig)
    _close_fig(_extra_fig)


# ----------------------------------------------------------------------
# Reference / derived tables / optional heavy views
# ----------------------------------------------------------------------
info = scheme.info()
if info:
    with st.expander("Reference / about"):
        st.markdown(info)

for table in view.get("tables", []):
    with st.expander(table["title"]):
        st.markdown(table["markdown"])

for view_def in scheme.extra_views():
    with st.expander(view_def.key):
        st.caption(view_def.description)
        if st.button("Run", key=f"run__{scheme.name}__{view_def.key}"):
            with st.spinner("Running…"):
                data = _cached_extra(scheme.name, view_def.key, recompute_items, cache_version)
            extra_fig = view_def.render(data)
            st.markdown("<div class='gabes-plot-gap'></div>", unsafe_allow_html=True)
            apply_gabes_plot_style(extra_fig)
            st.pyplot(extra_fig)
            _close_fig(extra_fig)
