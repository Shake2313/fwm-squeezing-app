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

APP_DIR = Path(__file__).resolve().parent


@st.cache_data(show_spinner=False)
def _asset_text(filename):
    return (APP_DIR / "assets" / filename).read_text(encoding="utf-8")


LOGO_SVG = _asset_text("gabes-logo-v3.svg")
ICON_SVG = _asset_text("gabes-mark-v3.svg")
SIDEBAR_LOGO_SVG = LOGO_SVG.replace(
    'width="960" height="232" viewBox="0 0 960 232"',
    'width="738" height="232" viewBox="0 0 738 232"',
    1,
)

st.set_page_config(page_title="GABES — Atomic Bloch Equation Solver",
                   page_icon=ICON_SVG, layout="wide")


# ----------------------------------------------------------------------
# Cached compute layer (keyed on the scheme + its recompute knobs only)
# ----------------------------------------------------------------------
@st.cache_data(show_spinner=False)
def _cached_compute(scheme_name, recompute_items, cache_version):
    return schemes.get(scheme_name).compute(dict(recompute_items))


@st.cache_data(show_spinner=False)
def _cached_extra(scheme_name, view_key, recompute_items, cache_version):
    scheme = schemes.get(scheme_name)
    view = next(v for v in scheme.extra_views() if v.key == view_key)
    return view.compute(dict(recompute_items))


@st.cache_data(show_spinner=False, max_entries=64)
def _cached_observables(scheme_name, raw, param_items, cache_version):
    return schemes.get(scheme_name).observables(raw, dict(param_items))


def _close_fig(fig):
    import matplotlib.pyplot as plt
    plt.close(fig)


def _skey(scheme_name, pname):
    return f"{scheme_name}__{pname}"


def _render_param(container, scheme_name, sp):
    key = _skey(scheme_name, sp.name)
    label = sp.label + (f"  [{sp.unit}]" if sp.unit else "")
    help_ = sp.help or None
    if sp.choices is not None:
        return container.selectbox(label, list(sp.choices), key=key, help=help_)
    val = container.slider(label, sp.vmin, sp.vmax, step=sp.step, key=key, help=help_)
    endpoints = getattr(sp, "endpoints", None)
    if endpoints:
        left, right = endpoints
        container.markdown(
            "<div style='display:flex;justify-content:space-between;margin-top:-12px;"
            f"font-size:0.72em;color:#888'><span>{left}</span><span>{right}</span></div>",
            unsafe_allow_html=True)
    return val


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
for preset in scheme.presets():
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
        _rec_sets = _rec_fn(scheme.defaults())
    except Exception:
        _rec_sets = None
if isinstance(_rec_sets, dict) and _rec_sets:
    _cols = st.sidebar.columns(len(_rec_sets))
    for _col, _label in zip(_cols, _rec_sets):
        def _apply_default(sname=scheme.name, sc=scheme, lbl=_label):
            cur = {sp.name: st.session_state[_skey(sname, sp.name)] for sp in sc.param_schema()}
            sets = sc.recommended_defaults(cur) or {}
            for k, v in (sets.get(lbl) or {}).items():
                st.session_state[_skey(sname, k)] = v
        _col.button(_label, on_click=_apply_default, use_container_width=True)

# Controls — grouped sections; advanced/numeric knobs fold into an expander.
params = {}
group_order = []
for sp in specs:
    if not sp.advanced and sp.group not in group_order:
        group_order.append(sp.group)

for g in group_order:
    st.sidebar.subheader(g)
    for sp in specs:
        if sp.group == g and not sp.advanced:
            params[sp.name] = _render_param(st.sidebar, scheme.name, sp)

advanced = [sp for sp in specs if sp.advanced]
if advanced:
    exp = st.sidebar.expander("Advanced / numerics")
    for sp in advanced:
        params[sp.name] = _render_param(exp, scheme.name, sp)


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
st.title(scheme.title)
if scheme.caption:
    st.caption(scheme.caption)

metrics = view.get("metrics", [])
if metrics:
    cols = st.columns(len(metrics))
    for col, m in zip(cols, metrics):
        kwargs = {}
        if m.get("delta") is not None:
            kwargs["delta"] = m["delta"]
            kwargs["delta_color"] = m.get("delta_color", "normal")
        if m.get("help"):
            kwargs["help"] = m["help"]
        col.metric(m["label"], m["value"], **kwargs)

fig = view.get("figure")
if fig is not None:
    st.pyplot(fig)
    _close_fig(fig)

for _title, _extra_fig in view.get("figures", []):
    st.subheader(_title)
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
            st.pyplot(extra_fig)
            _close_fig(extra_fig)
