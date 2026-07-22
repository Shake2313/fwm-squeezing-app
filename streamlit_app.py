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
import base64
import hashlib
import inspect
import json
import matplotlib
matplotlib.use("Agg")          # headless server backend (no GUI / Tk)
import numpy as np
import streamlit as st
import streamlit.components.v1 as components
from io import BytesIO
from pathlib import Path
from html import escape
from threading import RLock

from gabes import schemes
from gabes.core import blas_single_thread
from gabes.experimental_csv import (
    MAX_FILE_BYTES,
    ExperimentalCSVError,
    load_experimental_csv,
)
from gabes.plot_style import PALETTE, apply_gabes_plot_style
from gabes.ui_metrics import partition_metrics, split_metric_value

# Whether this Streamlit exposes the per-widget upload cap (added after the
# project's oldest supported release). A static fact about the installed
# version — probe once, not on every rerun.
_FILE_UPLOADER_HAS_MAX_SIZE = (
    "max_upload_size" in inspect.signature(st.file_uploader).parameters
)

APP_DIR = Path(__file__).resolve().parent
_PLOT_LOCK = RLock()


@st.cache_data(show_spinner=False)
def _asset_text(filename):
    return (APP_DIR / "assets" / filename).read_text(encoding="utf-8")


THEME_BASE = (st.get_option("theme.base") or "light").lower()
LOGO_ASSET = "gabes-logo-v3-dark.svg" if THEME_BASE == "dark" else "gabes-logo-v3.svg"
ICON_ASSET = "gabes-mark-v3-dark.svg" if THEME_BASE == "dark" else "gabes-mark-v3.svg"

# User's Guide, served as a static file (config.toml -> server.enableStaticServing).
# Streamlit exposes ./static/<f> at the relative URL "app/static/<f>", which
# resolves correctly both locally and on Streamlit Community Cloud — so the link
# opens from any computer. The file is fully self-contained (images base64-inlined
# by docs/build_static_guide.py), so it needs no sibling assets.
GUIDE_URL = "app/static/GABES_User_Guide.html"

# BETA badge appended to the wordmark logo (this is "the GABES logo" the app shows).
_BETA_BADGE = (
    '<g transform="translate(582 74)">'
    '<rect width="92" height="38" rx="19" fill="#F43F5E"/>'
    '<text x="46" y="26" text-anchor="middle" '
    'font-family="IBM Plex Sans, Inter, Segoe UI, Arial, sans-serif" '
    'font-size="21" font-weight="800" letter-spacing="2.5" fill="#FFFFFF">BETA</text>'
    '</g></svg>'
)


def _with_beta_badge(svg):
    """Insert the BETA badge just before the closing </svg> tag."""
    return svg.replace('</svg>', _BETA_BADGE, 1)


# Guide launcher: a tiny in-app button that opens the User's Guide in a new tab,
# RENDERED. Streamlit's static handler serves the file as text/plain (a security
# default), so a plain <a> would show source. Instead we fetch the file and open
# it as a text/html Blob — works on any computer with no external hosting. The
# heavy 851 KB file is fetched only on click (then browser-cached); the launcher
# itself is ~2 KB so it costs nothing on rerun.
_GUIDE_LAUNCHER_TMPL = """<!doctype html><html><head><meta charset="utf-8"><style>
html,body{margin:0;padding:0;background:transparent;overflow:hidden;
  font-family:"Pretendard","Apple SD Gothic Neo","Malgun Gothic","Segoe UI",sans-serif;}
.bar{display:flex;justify-content:__JUSTIFY__;align-items:center;}
button.gbtn{display:inline-flex;align-items:center;justify-content:center;gap:.4rem;
  cursor:pointer;__WIDTH__padding:.52rem .95rem;border-radius:9px;border:1px solid #1D4ED8;
  background:linear-gradient(100deg,#0369A1,#2563EB);color:#fff;
  font-weight:720;font-size:.85rem;line-height:1;
  box-shadow:0 2px 10px rgba(37,99,235,.22);transition:filter .12s ease;}
button.gbtn:hover{filter:brightness(1.08);}
button.gbtn .arr{font-size:.95rem;opacity:.9;}
</style></head><body>
<div class="bar">
  <button class="gbtn" id="gbtn" title="사용자 안내서를 새 창에서 엽니다">
    &#128218; User&rsquo;s Guide <span class="arr">&#8599;</span></button>
</div>
<script>
(function(){
  var GURL="__URL__";
  var btn=document.getElementById("gbtn");
  btn.addEventListener("click", function(){
    var old=btn.innerHTML; btn.disabled=true; btn.innerHTML="\\uC5EC\\uB294 \\uC911\\u2026";
    var done=function(){ btn.disabled=false; btn.innerHTML=old; };
    fetch(GURL).then(function(r){ return r.text(); }).then(function(t){
      var u=URL.createObjectURL(new Blob([t],{type:"text/html"}));
      if(!window.open(u,"_blank")){ window.open(GURL,"_blank"); }
      done();
    }).catch(function(){ window.open(GURL,"_blank"); done(); });
  });
})();
</script>
</body></html>"""


def _guide_launcher(container=None, height=56, align="end", full_width=False):
    justify = {"end": "flex-end", "center": "center",
               "start": "flex-start"}.get(align, "flex-end")
    html = (_GUIDE_LAUNCHER_TMPL
            .replace("__JUSTIFY__", justify)
            .replace("__WIDTH__", "width:100%;" if full_width else "")
            .replace("__URL__", GUIDE_URL))
    if container is None:
        components.html(html, height=height)
    else:
        with container:
            components.html(html, height=height)


LOGO_SVG = _with_beta_badge(_asset_text(LOGO_ASSET))
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
    "cell": dict(label="cell", color="#0F766E", bg="#EAF8F5"),
    "beams": dict(label="beams", color="#0891B2", bg="#E0F7FA"),
    "cell geometry": dict(label="geometry", color="#0F766E", bg="#EAF8F5"),
    "systematics": dict(label="systematics", color="#F43F5E", bg="#FFF1F3"),
    "relaxation overrides": dict(label="relaxation", color="#475569", bg="#F1F5F9"),
    "calibration": dict(label="calibration", color="#7C3AED", bg="#F3E8FF"),
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
# Manual cache-buster for _cached_observables: bump when the metrics/figure
# rendering *shape* changes (new tiers, ribbon layout) so stale cached views drop.
READOUT_CACHE_VERSION = "hero-ribbon-v3-single-hero"


def _inject_css():
    st.markdown("""
<style>
:root {
  --gabes-bg: #F6F8FB;
  --gabes-surface: #FFFFFF;
  --gabes-ink: #0F172A;
  --gabes-muted: #64748B;
  --gabes-subtle-ink: #475569;
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
  border: 1px solid rgba(220, 230, 239, 0.78);
  border-radius: 8px;
  background: rgba(255, 255, 255, 0.78);
  box-shadow: none;
}

[data-testid="stSidebar"] [data-testid="stExpander"] {
  margin-top: 0.48rem;
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
  margin: 0.95rem 0 0.28rem;
  padding: 0.12rem 0.04rem 0.1rem;
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

.gabes-advanced-subheader {
  display: flex;
  align-items: center;
  gap: 0.4rem;
  margin: 0.72rem 0 0.18rem;
  padding-top: 0.08rem;
  color: var(--gabes-muted);
  font-size: 0.74rem;
  line-height: 1.2;
  font-weight: 760;
}

.gabes-advanced-subheader::before {
  content: "";
  width: 0.36rem;
  height: 0.36rem;
  border-radius: 999px;
  background: var(--sub-color);
  flex: 0 0 auto;
}

.gabes-endpoints {
  display: flex;
  justify-content: space-between;
  margin-top: -0.72rem;
  color: var(--gabes-muted);
  font-size: 0.78rem;
}

.gabes-readout {
  display: grid;
  width: 100%;
  max-width: 50rem;
  margin-inline: auto;
  gap: 0.46rem;
}

.gabes-hero-grid {
  display: grid;
  grid-template-columns: repeat(2, minmax(0, 1fr));
  gap: 0.46rem;
}

.gabes-hero-grid--single {
  grid-template-columns: minmax(0, 1fr);
}

.gabes-hero-card {
  position: relative;
  min-width: 0;
  min-height: 5.15rem;
  padding: 0.72rem 0.86rem 0.68rem;
  overflow: hidden;
  background: var(--gabes-surface);
  border: 1px solid var(--gabes-border);
  border-radius: 10px;
  box-shadow: 0 3px 12px rgba(15, 23, 42, 0.035);
}

.gabes-hero-card--primary {
  background: linear-gradient(125deg, var(--metric-bg) 0%, var(--gabes-surface) 82%);
  border-color: var(--metric-color);
  box-shadow: 0 4px 16px rgba(15, 23, 42, 0.055);
}

.gabes-hero-card--primary::after {
  content: "";
  position: absolute;
  inset: 0 auto 0 0;
  width: 3px;
  background: var(--metric-color);
}

.gabes-hero-label,
.gabes-ribbon-label {
  min-width: 0;
  color: var(--gabes-muted);
  font-weight: 680;
  line-height: 1.25;
}

.gabes-hero-label {
  font-size: 0.8rem;
  white-space: nowrap;
  overflow: hidden;
  text-overflow: ellipsis;
}

.gabes-hero-value {
  min-width: 0;
  margin-top: 0.18rem;
  color: var(--gabes-ink);
  font-size: clamp(1.45rem, 2.65vw, 2.08rem);
  font-weight: 780;
  line-height: 1.08;
  letter-spacing: -0.02em;
  overflow-wrap: anywhere;
}

.gabes-hero-unit {
  margin-left: 0.28em;
  font-size: 0.58em;
  font-weight: 700;
  letter-spacing: 0;
  white-space: nowrap;
}

.gabes-hero-card--primary .gabes-hero-label,
.gabes-hero-card--primary .gabes-hero-value {
  color: var(--gabes-ink);
}

.gabes-metric-delta {
  margin-top: 0.13rem;
  color: var(--gabes-subtle-ink);
  font-size: 0.68rem;
  line-height: 1.2;
}

.gabes-metric-ribbon {
  display: grid;
  grid-template-columns: repeat(auto-fit, minmax(112px, 1fr));
  gap: 1px;
  overflow: hidden;
  background: var(--gabes-border);
  border: 1px solid var(--gabes-border);
  border-radius: 9px;
}

.gabes-ribbon-item {
  min-width: 0;
  min-height: 3.55rem;
  padding: 0.48rem 0.68rem 0.5rem;
  background: var(--gabes-surface);
}

.gabes-ribbon-label {
  font-size: 0.68rem;
  white-space: nowrap;
  overflow: hidden;
  text-overflow: ellipsis;
}

.gabes-ribbon-value {
  min-width: 0;
  margin-top: 0.12rem;
  color: var(--gabes-ink);
  font-size: clamp(0.88rem, 1.15vw, 1.02rem);
  font-weight: 740;
  line-height: 1.15;
  overflow-wrap: anywhere;
}

@media (max-width: 520px) {
  .gabes-hero-grid {
    grid-template-columns: minmax(0, 1fr);
  }

  .gabes-hero-card {
    min-height: 4.65rem;
  }

  .gabes-metric-ribbon {
    grid-template-columns: repeat(2, minmax(0, 1fr));
  }
}

@media (max-width: 330px) {
  .gabes-metric-ribbon {
    grid-template-columns: minmax(0, 1fr);
  }
}

.gabes-section-gap {
  height: 0.42rem;
}

.gabes-plot-gap {
  height: 0.56rem;
}

/* Cap chart width so a fixed-aspect figure can't overflow viewport height in the no-max-width wide layout. */
[data-testid="stMain"] [data-testid="stImage"] img {
  max-width: 50rem !important;
  height: auto !important;
}

/* stFullScreenFrame stays full-width even once the image inside it is capped; re-center the capped image within it. */
[data-testid="stMain"] [data-testid="stFullScreenFrame"] {
  display: flex;
  justify-content: center;
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
@st.cache_data(show_spinner=False, max_entries=64)
def _cached_compute(scheme_name, recompute_items, cache_version):
    with blas_single_thread():
        return schemes.get(scheme_name).compute(dict(recompute_items))


@st.cache_data(show_spinner=False, max_entries=16)
def _cached_extra(scheme_name, view_key, recompute_items, cache_version):
    scheme = schemes.get(scheme_name)
    view = next(v for v in scheme.extra_views() if v.key == view_key)
    with blas_single_thread():
        return view.compute(dict(recompute_items))


@st.cache_data(show_spinner=False, max_entries=64)
def _cached_observables(scheme_name, recompute_items, param_items, cache_version):
    # `raw` (the full solve output — large numpy arrays) is 1:1 determined by
    # (scheme_name, recompute_items), so key on those small tuples and fetch raw
    # from the compute cache instead of hashing the whole solve on every rerun.
    # This keeps navigate-only knobs cheap (the two-tier design). cache_version is
    # (compute_cache_version, READOUT_CACHE_VERSION); [0] matches _cached_compute.
    raw = _cached_compute(scheme_name, recompute_items, cache_version[0])
    # Matplotlib's font/mathtext/layout caches are process-global. Streamlit can
    # briefly overlap reruns when sliders are moved quickly, so serialize figure
    # construction to avoid layout-time parser crashes.
    with _PLOT_LOCK:
        return schemes.get(scheme_name).observables(raw, dict(param_items))


@st.cache_data(show_spinner=False, max_entries=16)
def _cached_experimental_csv(csv_bytes, denoise):
    """Parse/correct uploaded scope data outside the physics solve cache."""
    return load_experimental_csv(csv_bytes, denoise=denoise)


def _close_fig(fig):
    import matplotlib.pyplot as plt
    plt.close(fig)


def _render_fig(fig):
    """Style, draw, then release a figure (matplotlib figures leak if not closed)."""
    apply_gabes_plot_style(fig)
    st.pyplot(fig)
    _close_fig(fig)


def _figure_data_url(fig):
    """Serialize one styled Matplotlib figure for the client-side carousel."""
    buffer = BytesIO()
    apply_gabes_plot_style(fig)
    fig.savefig(
        buffer,
        format="png",
        dpi=160,
        bbox_inches="tight",
        pad_inches=0.14,
        facecolor=fig.get_facecolor(),
        edgecolor="none",
    )
    encoded = base64.b64encode(buffer.getvalue()).decode("ascii")
    return f"data:image/png;base64,{encoded}"


def _plot_carousel_html(pages):
    """Build an isolated, responsive carousel with arrows and touch dragging."""
    labels = [str(page["label"]) for page in pages]
    labels_json = json.dumps(labels, ensure_ascii=False).replace("<", "\\u003c")
    slides = "".join(
        '<div class="slide" role="group" '
        f'aria-label="{escape(label, quote=True)}" '
        f'aria-hidden="{"false" if index == 0 else "true"}">'
        f'<img src="{url}" alt="{escape(label, quote=True)} graph" '
        'draggable="false"></div>'
        for index, (label, url) in enumerate(
            zip(labels, (page["url"] for page in pages))
        )
    )
    dots = "".join(
        '<button class="dot" type="button" '
        f'aria-label="Show {escape(label, quote=True)}" data-index="{index}"></button>'
        for index, label in enumerate(labels)
    )
    surface = "#FFFFFF" if THEME_BASE != "dark" else "#0E1117"
    ink = "#0F172A" if THEME_BASE != "dark" else "#F1F5F9"
    muted = "#64748B" if THEME_BASE != "dark" else "#94A3B8"
    border = "#DCE6EF" if THEME_BASE != "dark" else "#334155"
    return f"""<!doctype html>
<html><head><meta charset="utf-8"><style>
*{{box-sizing:border-box}}
html,body{{margin:0;padding:0;background:transparent;overflow:hidden;
  font-family:"Pretendard","Apple SD Gothic Neo","Malgun Gothic","Segoe UI",sans-serif;
  color:{ink};}}
.shell{{position:relative;width:100%;padding:.15rem 3.15rem 0;}}
.viewport{{width:100%;overflow:hidden;touch-action:pan-y pinch-zoom;border:1px solid {border};
  border-radius:14px;background:{surface};box-shadow:0 2px 12px rgba(15,23,42,.05);
  cursor:grab;}}
.viewport.dragging{{cursor:grabbing;}}
.track{{display:flex;width:100%;transform:translate3d(0,0,0);
  transition:transform .28s cubic-bezier(.22,.75,.25,1);will-change:transform;}}
.slide{{flex:0 0 100%;min-width:100%;display:grid;place-items:center;background:{surface};}}
.slide img{{display:block;width:100%;height:auto;max-height:650px;object-fit:contain;
  user-select:none;-webkit-user-drag:none;}}
.nav{{position:absolute;z-index:3;top:calc(50% - 1.9rem);width:2.55rem;height:2.55rem;
  border-radius:999px;border:1px solid #BAE6FD;background:#FFFFFF;color:#0369A1;
  box-shadow:0 3px 12px rgba(2,132,199,.18);font-size:1.75rem;line-height:1;
  display:grid;place-items:center;cursor:pointer;transition:.14s ease;}}
.nav:hover:not(:disabled){{background:#E6F7FC;transform:scale(1.04);}}
.nav:focus-visible,.dot:focus-visible,.shell:focus-visible{{outline:2px solid #0284C7;
  outline-offset:2px;}}
.nav:disabled{{opacity:.28;cursor:default;box-shadow:none;}}
.prev{{left:.2rem}} .next{{right:.2rem}}
.meta{{min-height:2.45rem;display:flex;align-items:center;justify-content:center;
  gap:.75rem;color:{muted};font-size:.82rem;padding:.45rem .25rem .15rem;}}
.label{{min-width:8.5rem;color:{ink};font-weight:700;text-align:right;}}
.dots{{display:flex;gap:.38rem;align-items:center;}}
.dot{{width:.5rem;height:.5rem;padding:0;border:0;border-radius:999px;background:{border};
  cursor:pointer;transition:width .18s ease,background .18s ease;}}
.dot.active{{width:1.35rem;background:#0284C7;}}
.hint-mobile{{display:none}}
@media (max-width:700px), (hover:none), (pointer:coarse){{
  .shell{{padding:.1rem 0 0}} .nav{{display:none}}
  .viewport{{border-radius:10px}} .meta{{font-size:.78rem;gap:.55rem}}
  .label{{min-width:0;text-align:center}} .hint-desktop{{display:none}}
  .hint-mobile{{display:inline}}
}}
</style></head><body>
<div class="shell" id="shell" tabindex="0" aria-label="Spectrum graph carousel">
  <button class="nav prev" type="button" aria-label="Previous graph">&#8249;</button>
  <div class="viewport" id="viewport"><div class="track" id="track">{slides}</div></div>
  <button class="nav next" type="button" aria-label="Next graph">&#8250;</button>
  <div class="meta"><span class="label" id="label" aria-live="polite"></span>
    <span class="dots">{dots}</span>
    <span class="hint-desktop">Use arrows</span><span class="hint-mobile">Swipe</span>
  </div>
</div>
<script>
(function(){{
  const labels={labels_json};
  const shell=document.getElementById("shell");
  const viewport=document.getElementById("viewport");
  const track=document.getElementById("track");
  const label=document.getElementById("label");
  const prev=document.querySelector(".prev");
  const next=document.querySelector(".next");
  const dots=Array.from(document.querySelectorAll(".dot"));
  const slideNodes=Array.from(document.querySelectorAll(".slide"));
  let index=0,startX=null,startY=null,lastX=null,horizontal=false,activePointer=null;

  function resizeHost(){{
    const height=Math.ceil(shell.getBoundingClientRect().bottom+3);
    try{{
      const frames=Array.from(window.parent.document.querySelectorAll("iframe"));
      const frame=frames.find(item => item.contentWindow===window);
      if(frame){{ frame.style.height=height+"px"; frame.setAttribute("height",height); }}
    }}catch(_error){{}}
  }}
  function show(nextIndex,animate=true){{
    index=Math.max(0,Math.min(labels.length-1,nextIndex));
    track.style.transition=animate?"transform .28s cubic-bezier(.22,.75,.25,1)":"none";
    track.style.transform=`translate3d(${{-index*100}}%,0,0)`;
    label.textContent=`${{labels[index]}}  ·  ${{index+1}}/${{labels.length}}`;
    prev.disabled=index===0; next.disabled=index===labels.length-1;
    slideNodes.forEach((slide,i)=>{{
      const active=i===index;
      slide.setAttribute("aria-hidden",active?"false":"true");
      if("inert" in slide)slide.inert=!active;
    }});
    dots.forEach((dot,i)=>{{dot.classList.toggle("active",i===index);
      dot.setAttribute("aria-current",i===index?"true":"false");}});
    requestAnimationFrame(resizeHost);
  }}
  prev.addEventListener("click",()=>show(index-1));
  next.addEventListener("click",()=>show(index+1));
  dots.forEach(dot=>dot.addEventListener("click",()=>show(Number(dot.dataset.index))));
  shell.addEventListener("keydown",event=>{{
    if(event.key==="ArrowLeft"){{event.preventDefault();show(index-1);}}
    if(event.key==="ArrowRight"){{event.preventDefault();show(index+1);}}
  }});
  viewport.addEventListener("pointerdown",event=>{{
    if(!event.isPrimary || activePointer!==null)return;
    if(event.pointerType==="mouse" && event.button!==0)return;
    activePointer=event.pointerId;
    startX=lastX=event.clientX;startY=event.clientY;horizontal=false;
    viewport.classList.add("dragging");
    try{{viewport.setPointerCapture(event.pointerId);}}catch(_error){{}}
  }});
  viewport.addEventListener("pointermove",event=>{{
    if(startX===null || event.pointerId!==activePointer)return;
    lastX=event.clientX;
    const dx=lastX-startX,dy=event.clientY-startY;
    if(!horizontal && Math.max(Math.abs(dx),Math.abs(dy))>7){{
      horizontal=Math.abs(dx)>Math.abs(dy);
    }}
    if(!horizontal)return;
    event.preventDefault();
    const edge=(index===0 && dx>0)||(index===labels.length-1 && dx<0);
    const drag=edge?dx*.24:dx;
    track.style.transition="none";
    track.style.transform=`translate3d(calc(${{-index*100}}% + ${{drag}}px),0,0)`;
  }});
  function finish(event){{
    if(startX===null || event.pointerId!==activePointer)return;
    const dx=(lastX===null?startX:lastX)-startX;
    const threshold=Math.min(90,Math.max(42,viewport.clientWidth*.12));
    viewport.classList.remove("dragging");
    if(horizontal && Math.abs(dx)>=threshold)show(index+(dx<0?1:-1));
    else show(index);
    startX=startY=lastX=null;horizontal=false;activePointer=null;
  }}
  viewport.addEventListener("pointerup",finish);
  viewport.addEventListener("pointercancel",finish);
  document.querySelectorAll("img").forEach(img=>{{
    img.addEventListener("dragstart",event=>event.preventDefault());
    if(!img.complete)img.addEventListener("load",resizeHost,{{once:true}});
  }});
  window.addEventListener("resize",resizeHost);
  if(window.ResizeObserver)new ResizeObserver(resizeHost).observe(shell);
  show(0,false);setTimeout(resizeHost,60);
}})();
</script></body></html>"""


def _render_figure_views(figure_views):
    """Render labelled figures as a no-recompute desktop/mobile carousel."""
    pages = [
        {"label": item.get("label", f"View {index + 1}"), "figure": item.get("figure")}
        for index, item in enumerate(figure_views)
        if isinstance(item, dict) and item.get("figure") is not None
    ]
    if not pages:
        return
    if len(pages) == 1:
        _render_fig(pages[0]["figure"])
        return

    figures = [page["figure"] for page in pages]
    try:
        with _PLOT_LOCK:
            for page in pages:
                page["url"] = _figure_data_url(page["figure"])
        components.html(_plot_carousel_html(pages), height=720, scrolling=False)
    finally:
        closed = set()
        for fig in figures:
            if id(fig) not in closed:
                _close_fig(fig)
                closed.add(id(fig))


def _render_experimental_comparison(view, scheme_name):
    """Render a scheme-declared CSV panel and overlay its corrected trace.

    Uploaded bytes and alignment controls intentionally stay outside `params`:
    changing them must reuse both the heavy solve and the cached base figure.
    """
    descriptor = view.get("comparison")
    fig = view.get("figure")
    if not descriptor or fig is None:
        return

    axis_index = int(descriptor.get("axis_index", 0))
    if not 0 <= axis_index < len(fig.axes):
        return
    axis = fig.axes[axis_index]
    xlim = tuple(float(v) for v in axis.get_xlim())
    ylim = tuple(float(v) for v in axis.get_ylim())
    x_unit = descriptor.get("x_unit", "plot unit")
    raw_x_unit = descriptor.get("raw_x_unit", "Arb. unit")
    raw_y_unit = descriptor.get("raw_y_unit", "Arb. unit")

    panel = st.expander("Experimental CSV comparison")
    with panel:
        panel.caption(
            "Column A = detuning, column B = detector signal. Later columns and "
            "non-numeric oscilloscope metadata rows are ignored. GABES processes "
            "up to 10 MiB or 500,000 rows."
        )
        uploader_options = {}
        if _FILE_UPLOADER_HAS_MAX_SIZE:
            # Backend byte check remains the fallback on releases without it.
            uploader_options["max_upload_size"] = MAX_FILE_BYTES // (1024 * 1024)
        uploaded = panel.file_uploader(
            "Oscilloscope CSV",
            type=("csv",),
            key=_skey(scheme_name, "_csv_upload"),
            help="A and B are read as arbitrary units; no header is required.",
            **uploader_options,
        )
        auto_correct = panel.checkbox(
            "Automatic noise correction",
            value=True,
            key=_skey(scheme_name, "_csv_auto_correct"),
            help="Merge repeated A values, reject isolated spikes, and apply "
                 "noise-adaptive local smoothing before 0-1 calibration.",
        )
        show_raw = panel.checkbox(
            "Show unfiltered trace",
            value=False,
            disabled=not auto_correct,
            key=_skey(scheme_name, "_csv_show_raw"),
        )
        if uploaded is None:
            return

        csv_bytes = uploaded.getvalue()
        fingerprint = hashlib.sha256(csv_bytes).hexdigest()
        fingerprint_key = _skey(scheme_name, "_csv_fingerprint")
        scale_key = _skey(scheme_name, "_csv_x_scale")
        shift_key = _skey(scheme_name, "_csv_x_shift")
        reverse_key = _skey(scheme_name, "_csv_reverse")
        invert_key = _skey(scheme_name, "_csv_invert")
        alignment_identity = (
            fingerprint,
            axis_index,
            str(x_unit),
            str(raw_x_unit),
        )
        if st.session_state.get(fingerprint_key) != alignment_identity:
            st.session_state[fingerprint_key] = alignment_identity
            st.session_state[scale_key] = 1.0
            st.session_state[shift_key] = 0.0
            st.session_state[reverse_key] = False
            st.session_state[invert_key] = False
        st.session_state.setdefault(scale_key, 1.0)
        st.session_state.setdefault(shift_key, 0.0)
        st.session_state.setdefault(reverse_key, False)
        st.session_state.setdefault(invert_key, False)

        try:
            trace = _cached_experimental_csv(csv_bytes, auto_correct)
        except ExperimentalCSVError as exc:
            panel.error(str(exc))
            return

        detuning = trace.detuning
        n_detuning = len(detuning)
        pivot = float(
            (detuning[(n_detuning - 1) // 2] + detuning[n_detuning // 2]) / 2.0
        )
        raw_span = float(detuning[-1] - detuning[0])
        plot_span = float(xlim[1] - xlim[0])
        framed_scale = min(max(0.90 * plot_span / raw_span, 1e-9), 1e9)
        framed_shift = 0.5 * (xlim[0] + xlim[1]) - pivot

        action_cols = panel.columns(2)
        if action_cols[0].button(
            "Bring into view",
            key=_skey(scheme_name, "_csv_frame"),
            use_container_width=True,
            help="Map the imported span into the current theoretical plot; "
                 "this does not fit spectral features.",
        ):
            st.session_state[scale_key] = float(framed_scale)
            st.session_state[shift_key] = float(framed_shift)
        if action_cols[1].button(
            "Reset alignment",
            key=_skey(scheme_name, "_csv_reset"),
            use_container_width=True,
        ):
            st.session_state[scale_key] = 1.0
            st.session_state[shift_key] = 0.0
            st.session_state[reverse_key] = False
            st.session_state[invert_key] = False

        align_cols = panel.columns(2)
        scale_now = max(abs(float(st.session_state.get(scale_key, 1.0))), 1e-9)
        scale_step = max(scale_now * 0.01, 1e-6)
        scale_pivot = min(framed_scale, scale_now)
        scale_min = max(scale_pivot / 1000.0, 1e-9)
        scale_max = min(max(framed_scale, scale_now) * 1000.0, 1e9)
        x_scale = align_cols[0].slider(
            f"X scale  [{x_unit}/{raw_x_unit}]",
            scale_min,
            scale_max,
            step=float(scale_step),
            key=scale_key,
            help="Scale around the imported trace centre.",
        )
        shift_now = float(st.session_state.get(shift_key, 0.0))
        shift_step = max(abs(plot_span) / 1000.0, 1e-6)
        shift_margin = 2.0 * max(abs(plot_span), 1e-9)
        shift_min = min(xlim[0] - shift_margin, shift_now)
        shift_max = max(xlim[1] + shift_margin, shift_now)
        x_shift = align_cols[1].slider(
            f"X shift  [{x_unit}]",
            shift_min,
            shift_max,
            step=float(shift_step),
            key=shift_key,
            help="Shift the imported trace along the plotted x axis.",
        )
        option_cols = panel.columns(2)
        reverse = option_cols[0].checkbox(
            "Reverse sweep",
            key=reverse_key,
            help="Reverse the sign of the imported detuning sweep around its centre.",
        )
        invert = option_cols[1].checkbox(
            "Invert transmission",
            key=invert_key,
            help="Use when detector voltage polarity is opposite to transmission.",
        )

        import_diag = trace.import_diagnostics
        correction_diag = trace.correction_diagnostics
        panel.caption(
            f"{uploaded.name}: {import_diag.valid_rows} numeric A/B rows · "
            f"{len(detuning)} unique detuning points · "
            f"{import_diag.duplicate_rows_merged} duplicates merged · "
            f"{import_diag.ignored_rows} rows ignored"
        )
        panel.caption(
            f"Signal calibration [{raw_y_unit}]: "
            f"floor {correction_diag.floor:.6g} · "
            f"ceiling {correction_diag.ceiling:.6g} · "
            f"smoothing window {correction_diag.smoothing_window}"
        )
        if correction_diag.warning:
            panel.warning(str(correction_diag.warning))

        aligned_x = trace.transformed_detuning(
            scale=float(x_scale), shift=float(x_shift), reverse=bool(reverse)
        )
        aligned_y = 1.0 - trace.transmission if invert else trace.transmission
        in_view = ((aligned_x >= xlim[0]) & (aligned_x <= xlim[1])).any()
        if not in_view:
            panel.warning(
                "The imported trace is outside the theoretical x range. "
                "Use Bring into view, then refine X scale and X shift manually."
            )

        raw_overlay = None
        if auto_correct and show_raw:
            raw_x = trace.transformed_detuning(
                scale=float(x_scale), shift=float(x_shift), reverse=bool(reverse)
            )
            raw_y = np.clip(
                (trace.raw_signal - correction_diag.floor)
                / correction_diag.contrast,
                0.0,
                1.0,
            )
            raw_y = 1.0 - raw_y if invert else raw_y
            raw_overlay = (raw_x, raw_y)

    with _PLOT_LOCK:
        if raw_overlay is not None:
            axis.plot(
                raw_overlay[0], raw_overlay[1], color=PALETTE["muted"],
                ls=":", lw=1.0, alpha=0.45, label="CSV · unfiltered", zorder=2,
            )
        axis.plot(
            aligned_x, aligned_y, color=PALETTE["rose"], lw=1.4, alpha=0.88,
            label=descriptor.get("label", "Experimental CSV"), zorder=3,
        )
        # Arbitrary input units must never autoscale the theoretical plot away.
        axis.set_xlim(xlim)
        axis.set_ylim(ylim)
        axis.legend(loc="best")


def _skey(scheme_name, pname):
    return f"{scheme_name}__{pname}"


def _current_params(scheme_name, scheme_obj):
    """Live value of every param knob from session_state, falling back to defaults."""
    return {
        sp.name: st.session_state.get(_skey(scheme_name, sp.name), sp.default)
        for sp in scheme_obj.param_schema()
    }


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


def _render_advanced_subheader(container, group):
    style = _concept_style(group)
    container.markdown(
        "<div class='gabes-advanced-subheader' "
        f"style='--sub-color:{style['color']}'>{escape(group)}</div>",
        unsafe_allow_html=True,
    )


def _apply_recommended_defaults(scheme_name, scheme_obj, key):
    """on_change for a control whose selection applies the matching
    recommended_defaults set (so choosing it also resets that mode's knobs)."""
    selection = st.session_state.get(key)
    cur = _current_params(scheme_name, scheme_obj)
    sets = scheme_obj.recommended_defaults(cur) or {}
    defaults = sets.get(selection) or sets.get(cur.get("mode")) or {}
    for k, v in defaults.items():
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


def _metric_card_html(metric, *, hero=False, primary=False):
    label = str(metric.get("label", ""))
    value = str(metric.get("value", ""))
    delta = metric.get("delta")
    help_text = metric.get("help") or ""
    style = _metric_style(label)
    delta_html = ""
    if delta is not None:
        delta_html = f"<div class='gabes-metric-delta'>{escape(str(delta))}</div>"
    title = escape(str(help_text))
    if hero:
        primary_class = " gabes-hero-card--primary" if primary else ""
        number, unit = split_metric_value(value, kind=metric.get("kind"))
        if unit is None:
            value_html = escape(number)
        else:
            value_html = (
                f"<span class='gabes-hero-number'>{escape(number)}</span> "
                f"<span class='gabes-hero-unit'>{escape(unit)}</span>"
            )
        return (
            f"<article class='gabes-hero-card{primary_class}' "
            f"title='{title}' "
            f"style='--metric-color:{style['color']};--metric-bg:{style['bg']}'>"
            f"<div class='gabes-hero-label'>{escape(label)}</div>"
            f"<div class='gabes-hero-value'>{value_html}</div>"
            f"{delta_html}"
            "</article>"
        )
    return (
        "<div class='gabes-ribbon-item' role='listitem' "
        f"title='{title}'>"
        f"<div class='gabes-ribbon-label'>{escape(label)}</div>"
        f"<div class='gabes-ribbon-value'>{escape(value)}</div>"
        f"{delta_html}"
        "</div>"
    )


def _render_metrics(metrics, *, hero_count=2):
    heroes, ribbon = partition_metrics(metrics, hero_count=hero_count)
    single_class = " gabes-hero-grid--single" if len(heroes) == 1 else ""
    hero_cards = "".join(
        _metric_card_html(metric, hero=True, primary=(index == 0))
        for index, metric in enumerate(heroes)
    )
    ribbon_html = ""
    if ribbon:
        ribbon_cards = "".join(
            _metric_card_html(metric) for metric in ribbon
        )
        ribbon_html = (
            "<div class='gabes-metric-ribbon' role='list'>"
            f"{ribbon_cards}</div>"
        )
    st.markdown(
        "<section class='gabes-readout' aria-label='Key results'>"
        f"<div class='gabes-hero-grid{single_class}'>{hero_cards}</div>"
        f"{ribbon_html}</section>",
        unsafe_allow_html=True,
    )


# ----------------------------------------------------------------------
# Sidebar — scheme selection
# ----------------------------------------------------------------------
st.sidebar.image(SIDEBAR_LOGO_SVG, width=230)
_guide_launcher(st.sidebar, height=46, align="center", full_width=True)

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
        # transmission regimes vs an NMOR default).
        _rec_sets = _rec_fn(_current_params(scheme.name, scheme))
    except Exception:
        _rec_sets = None
# A control flagged applies_defaults (e.g. FWM's Mode) already applies these sets
# on selection, so the standalone "Default" buttons would just duplicate it.
_mode_driven_defaults = any(getattr(sp, "applies_defaults", False) for sp in specs)
def _apply_default_set(sname, sc, label):
    cur = _current_params(sname, sc)
    sets = sc.recommended_defaults(cur) or {}
    for k, v in (sets.get(label) or {}).items():
        st.session_state[_skey(sname, k)] = v

if isinstance(_rec_sets, dict) and _rec_sets and not _mode_driven_defaults:
    if getattr(scheme, "recommended_defaults_as_dropdown", False):
        # Dropdown that keeps the chosen regime visible; picking one loads its
        # full parameter set. Defaults to the first regime offered.
        _options = list(_rec_sets)
        _dkey = _skey(scheme.name, "_default_choice")
        if st.session_state.get(_dkey) not in _options:
            st.session_state[_dkey] = _options[0]

        def _apply_default_dropdown(sname=scheme.name, sc=scheme, dkey=_dkey):
            _apply_default_set(sname, sc, st.session_state.get(dkey))

        st.sidebar.selectbox("Default", _options, key=_dkey,
                             on_change=_apply_default_dropdown,
                             help="Load a ready-made regime's full parameter set.")
    else:
        _render_group_header(st.sidebar, "Default")
        _cols = st.sidebar.columns(len(_rec_sets))
        for _col, _label in zip(_cols, _rec_sets):
            def _apply_default(sname=scheme.name, sc=scheme, lbl=_label):
                _apply_default_set(sname, sc, lbl)
            _short = _label.replace(" default", "")
            _col.button(_short, on_click=_apply_default, use_container_width=True)

# Controls — grouped sections; advanced/numeric knobs fold into an expander.
visible_specs = [sp for sp in specs if _param_visible(scheme.name, sp)]
group_order = []
for sp in visible_specs:
    if not sp.advanced and sp.group not in group_order:
        group_order.append(sp.group)

for g in group_order:
    _render_group_header(st.sidebar, g)
    for sp in visible_specs:
        if sp.group == g and not sp.advanced:
            _render_param(st.sidebar, scheme.name, sp, scheme)

advanced = [sp for sp in visible_specs if sp.advanced]
if advanced:
    exp = st.sidebar.expander("Advanced controls")
    advanced_group_order = []
    for sp in advanced:
        group = getattr(sp, "advanced_group", "") or sp.group
        if group not in advanced_group_order:
            advanced_group_order.append(group)
    show_advanced_subgroups = len(advanced_group_order) > 1
    for group in advanced_group_order:
        group_specs = [
            sp for sp in advanced
            if (getattr(sp, "advanced_group", "") or sp.group) == group
        ]
        if show_advanced_subgroups:
            _render_advanced_subheader(exp, group)
        for sp in group_specs:
            _render_param(exp, scheme.name, sp, scheme)

# Every rendered widget wrote its value to session_state under _skey; hidden and
# currently-invisible specs weren't rendered but keep their last value (or
# default). _current_params gathers both in one pass.
params = _current_params(scheme.name, scheme)


# ----------------------------------------------------------------------
# Compute (cached) + observables
# ----------------------------------------------------------------------
recompute_items = tuple(sorted((k, params[k]) for k in scheme.recompute_keys()))
cache_version = getattr(scheme, "cache_version", "1")
param_items = tuple(sorted(params.items()))
if getattr(scheme, "cache_observables", False):
    with st.spinner("Solving Bloch equations…"):
        view = _cached_observables(
            scheme.name, recompute_items, param_items,
            (cache_version, READOUT_CACHE_VERSION),
        )
else:
    with st.spinner("Solving Bloch equations…"):
        raw = _cached_compute(scheme.name, recompute_items, cache_version)
        with _PLOT_LOCK:
            view = scheme.observables(raw, params)


# ----------------------------------------------------------------------
# Header + readout
# ----------------------------------------------------------------------
_guide_launcher(height=44, align="end")        # top bar: open the guide (new tab)
_render_scheme_header(scheme)

metrics = view.get("metrics", [])
if metrics:
    _render_metrics(metrics, hero_count=view.get("hero_count", 2))
    st.markdown("<div class='gabes-section-gap'></div>", unsafe_allow_html=True)

_render_experimental_comparison(view, scheme.name)

fig = view.get("figure")
if fig is not None:
    figure_views = view.get("figure_views", [])
    if figure_views:
        _render_figure_views(figure_views)
    else:
        _render_fig(fig)

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

for _title, _extra_fig in view.get("figures", []):
    with st.expander(f"Diagnostic plot · {_title}"):
        _render_fig(_extra_fig)

for view_def in scheme.extra_views():
    with st.expander(view_def.key):
        st.caption(view_def.description)
        if st.button("Run", key=f"run__{scheme.name}__{view_def.key}"):
            with st.spinner("Running…"):
                data = _cached_extra(scheme.name, view_def.key, recompute_items, cache_version)
            extra_fig = view_def.render(data)
            st.markdown("<div class='gabes-plot-gap'></div>", unsafe_allow_html=True)
            _render_fig(extra_fig)
