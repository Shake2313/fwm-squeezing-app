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
import importlib
import inspect
import matplotlib
matplotlib.use("Agg")          # headless server backend (no GUI / Tk)
import streamlit as st
from pathlib import Path
from html import escape
from threading import RLock

import gabes as _gabes
import gabes.experimental_csv as _experimental_csv
from gabes import schemes
from gabes.core import blas_single_thread

_EXPERIMENTAL_CSV_IMPORT_LOCK = _gabes.__dict__.setdefault(
    "_streamlit_experimental_csv_import_lock", RLock()
)
_EXPERIMENTAL_CSV_API = (
    "CALIBRATION_ABSOLUTE_DARK_REFERENCE",
    "CALIBRATION_ABSOLUTE_GAIN_OFFSET",
    "CALIBRATION_RELATIVE_EXTREMA",
    "MAX_FILE_BYTES",
    "ExperimentalCSVError",
    "load_experimental_csv",
)
with _EXPERIMENTAL_CSV_IMPORT_LOCK:
    if any(not hasattr(_experimental_csv, name) for name in _EXPERIMENTAL_CSV_API):
        importlib.invalidate_caches()
        _experimental_csv = importlib.reload(_experimental_csv)

    CALIBRATION_ABSOLUTE_DARK_REFERENCE = (
        _experimental_csv.CALIBRATION_ABSOLUTE_DARK_REFERENCE
    )
    CALIBRATION_ABSOLUTE_GAIN_OFFSET = (
        _experimental_csv.CALIBRATION_ABSOLUTE_GAIN_OFFSET
    )
    CALIBRATION_RELATIVE_EXTREMA = _experimental_csv.CALIBRATION_RELATIVE_EXTREMA
    MAX_FILE_BYTES = _experimental_csv.MAX_FILE_BYTES
    ExperimentalCSVError = _experimental_csv.ExperimentalCSVError
    load_experimental_csv = _experimental_csv.load_experimental_csv
from gabes.plot_style import PALETTE, apply_gabes_plot_style
from gabes_ui import export as ui_export
from gabes_ui import theme as ui_theme
from gabes_ui.guide import guide_button
from gabes_ui.readout import (
    caption_html,
    clean_choice_label,
    metrics_table_markdown,
    partition_readout,
    strip_html,
)

APP_DIR = Path(__file__).resolve().parent
_PLOT_LOCK = RLock()


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


@st.cache_data(show_spinner=False, max_entries=16)
def _cached_extra(scheme_name, view_key, param_items, cache_version):
    scheme = schemes.get(scheme_name)
    view = next(v for v in scheme.extra_views() if v.key == view_key)
    with blas_single_thread():
        return view.compute(dict(param_items))


@st.cache_data(show_spinner=False, max_entries=64)
def _cached_observables(scheme_name, raw, param_items, cache_version):
    # Matplotlib's font/mathtext/layout caches are process-global. Streamlit can
    # briefly overlap reruns when sliders are moved quickly, so serialize figure
    # construction to avoid layout-time parser crashes.
    with _PLOT_LOCK:
        return schemes.get(scheme_name).observables(raw, dict(param_items))


@st.cache_data(show_spinner=False, max_entries=16)
def _cached_experimental_csv(
    csv_bytes,
    denoise,
    calibration_mode=CALIBRATION_RELATIVE_EXTREMA,
    dark_signal=None,
    reference_signal=None,
    gain=None,
    offset=None,
):
    """Parse/correct uploaded scope data outside the physics solve cache."""
    return load_experimental_csv(
        csv_bytes,
        denoise=denoise,
        calibration_mode=calibration_mode,
        dark_signal=dark_signal,
        reference_signal=reference_signal,
        gain=gain,
        offset=offset,
    )


def _close_fig(fig):
    import matplotlib.pyplot as plt
    plt.close(fig)


def _render_fig(fig):
    """Style, draw, then release a figure (matplotlib figures leak if not closed)."""
    apply_gabes_plot_style(fig)
    st.pyplot(fig)
    _close_fig(fig)


def _diagnostic_value(obj, *names, default=None):
    for name in names:
        if hasattr(obj, name):
            return getattr(obj, name)
    return default


def _render_experimental_comparison(view, scheme_name, panel):
    """Render a scheme-declared CSV panel into `panel` and overlay its trace.

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

    with panel:
        panel.caption(
            "Column A = detuning, column B = detector signal. Later columns and "
            "non-numeric oscilloscope metadata rows are ignored. GABES processes "
            "up to 10 MiB or 500,000 rows."
        )
        uploader_options = {}
        if "max_upload_size" in inspect.signature(panel.file_uploader).parameters:
            # Streamlit added this per-widget guard after the project's oldest
            # supported release; retain the backend byte check as the fallback.
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
                 "noise-adaptive local smoothing before the selected calibration.",
        )
        show_raw = panel.checkbox(
            "Show unfiltered trace",
            value=False,
            disabled=not auto_correct,
            key=_skey(scheme_name, "_csv_show_raw"),
        )
        calibration_options = {
            "Relative normalized transmission (extrema-scaled)":
                CALIBRATION_RELATIVE_EXTREMA,
            "Absolute transmission (dark/reference)":
                CALIBRATION_ABSOLUTE_DARK_REFERENCE,
            "Absolute transmission (gain/offset)":
                CALIBRATION_ABSOLUTE_GAIN_OFFSET,
        }
        calibration_choice = panel.selectbox(
            "Transmission calibration",
            tuple(calibration_options),
            key=_skey(scheme_name, "_csv_calibration_mode"),
            help="Extrema scaling is relative only. Absolute transmission "
                 "requires measured dark/reference levels or an explicit "
                 "detector law: signal = offset + gain × transmission.",
        )
        calibration_mode = calibration_options[calibration_choice]
        dark_signal = reference_signal = gain = offset = None
        if calibration_mode == CALIBRATION_ABSOLUTE_DARK_REFERENCE:
            calibration_cols = panel.columns(2)
            dark_signal = calibration_cols[0].number_input(
                f"Dark signal  [{raw_y_unit}]",
                value=0.0,
                key=_skey(scheme_name, "_csv_dark_signal"),
                help="Detector signal measured at zero transmitted light.",
            )
            reference_signal = calibration_cols[1].number_input(
                f"Reference signal (T = 1)  [{raw_y_unit}]",
                value=1.0,
                key=_skey(scheme_name, "_csv_reference_signal"),
                help="Detector signal measured for unit-transmission reference light.",
            )
        elif calibration_mode == CALIBRATION_ABSOLUTE_GAIN_OFFSET:
            calibration_cols = panel.columns(2)
            gain = calibration_cols[0].number_input(
                f"Detector gain  [{raw_y_unit}]",
                value=1.0,
                key=_skey(scheme_name, "_csv_detector_gain"),
                help="Gain in signal = offset + gain × transmission; a negative "
                     "gain represents inverted detector polarity.",
            )
            offset = calibration_cols[1].number_input(
                f"Detector offset  [{raw_y_unit}]",
                value=0.0,
                key=_skey(scheme_name, "_csv_detector_offset"),
                help="Signal value corresponding to zero transmission.",
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
            trace = _cached_experimental_csv(
                csv_bytes,
                auto_correct,
                calibration_mode,
                dark_signal,
                reference_signal,
                gain,
                offset,
            )
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
        if trace.is_absolute_transmission:
            # Detector polarity is already part of the signed gain or the
            # dark/reference ordering. A second visual inversion would destroy
            # the declared absolute meaning.
            st.session_state[invert_key] = False
        invert = option_cols[1].checkbox(
            "Invert transmission",
            key=invert_key,
            disabled=trace.is_absolute_transmission,
            help=(
                "Use when detector voltage polarity is opposite to relative "
                "transmission. Absolute calibration takes polarity from its "
                "signed gain and cannot be inverted again."
            ),
        )

        import_diag = trace.import_diagnostics
        correction_diag = trace.correction_diagnostics
        valid_rows = _diagnostic_value(
            import_diag, "valid_rows", "accepted_rows", default="?"
        )
        ignored_rows = _diagnostic_value(
            import_diag, "ignored_rows", "skipped_rows", default="?"
        )
        merged_rows = _diagnostic_value(
            import_diag, "duplicate_rows_merged", "merged_rows",
            "duplicate_rows", default="?"
        )
        panel.caption(
            f"{uploaded.name}: {valid_rows} numeric A/B rows · "
            f"{len(detuning)} unique detuning points · {merged_rows} duplicates "
            f"merged · {ignored_rows} rows ignored"
        )
        panel.caption(
            f"Acquisition order: {import_diag.forward_sample_count} forward · "
            f"{import_diag.reverse_sample_count} reverse · "
            f"{import_diag.sweep_reversal_count} reversal(s) · "
            f"{import_diag.sweep_branch_count} preserved branch(es)"
        )
        if import_diag.directional_merge_warning:
            panel.warning(import_diag.directional_merge_warning)
        floor = _diagnostic_value(correction_diag, "floor", "floor_level")
        ceiling = _diagnostic_value(correction_diag, "ceiling", "ceiling_level")
        window = _diagnostic_value(
            correction_diag, "smoothing_window", "window_size", default=1
        )
        if correction_diag.is_absolute:
            panel.caption(
                f"{trace.transmission_label} [{correction_diag.calibration_source}]: "
                f"offset {correction_diag.calibration_offset:.6g} · "
                f"gain {correction_diag.calibration_gain:.6g} [{raw_y_unit}] · "
                f"smoothing window {window}"
            )
        elif floor is not None and ceiling is not None:
            panel.caption(
                f"{trace.transmission_label} [{raw_y_unit}]: relative floor "
                f"{float(floor):.6g} · relative ceiling {float(ceiling):.6g} · "
                f"smoothing window {window}"
            )
        warnings = _diagnostic_value(
            correction_diag, "warnings", "warning", default=()
        )
        if isinstance(warnings, str):
            warnings = (warnings,)
        for warning in warnings or ():
            panel.warning(str(warning))

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
            raw_y = trace.calibrate_signal(trace.raw_signal)
            raw_y = 1.0 - raw_y if invert else raw_y
            raw_overlay = (raw_x, raw_y)

        overlay_label = (
            f"{descriptor.get('label', 'Experimental CSV')} — "
            f"{trace.transmission_label.lower()}"
        )

    with _PLOT_LOCK:
        if raw_overlay is not None:
            axis.plot(
                raw_overlay[0], raw_overlay[1], color=PALETTE["muted"],
                ls=":", lw=1.0, alpha=0.45, label="CSV · unfiltered", zorder=2,
            )
        axis.plot(
            aligned_x, aligned_y, color=PALETTE["rose"], lw=1.4, alpha=0.88,
            label=overlay_label, zorder=3,
        )
        # Arbitrary input units must never autoscale the theoretical plot away.
        axis.set_xlim(xlim)
        axis.set_ylim(ylim)
        axis.legend(loc="best")
        comparison_warnings = ()
        if isinstance(warnings, str):
            comparison_warnings = (warnings,)
        else:
            comparison_warnings = tuple(warnings or ())
        return {
            "enabled": True,
            "kind": "experimental_csv_overlay",
            "overlay_filename": uploaded.name,
            "source_fingerprint_sha256": fingerprint,
            "descriptor": {
                "label": descriptor.get("label"),
                "axis_index": int(axis_index),
                "axis_x_unit": str(x_unit),
                "raw_x_unit": str(raw_x_unit),
                "raw_y_unit": str(raw_y_unit),
            },
            "calibration_mode": str(calibration_mode),
            "calibration_inputs": {
                "auto_correct": bool(auto_correct),
                "dark_signal": dark_signal,
                "reference_signal": reference_signal,
                "detector_gain": gain,
                "detector_offset": offset,
            },
            "alignment": {
                "scale": float(x_scale),
                "shift": float(x_shift),
                "reverse": bool(reverse),
                "invert": bool(invert),
            },
            "show_raw": bool(show_raw),
            "in_view": bool(in_view),
            "framed_scale": float(framed_scale),
            "framed_shift": float(framed_shift),
            "comparison_trace": {
                "detuning": ui_export.to_json_primitive(aligned_x),
                "transmission": ui_export.to_json_primitive(aligned_y),
                "detuning_unit": str(x_unit),
                "transmission_unit": str(raw_y_unit),
                "transmission_label": str(trace.transmission_label),
            },
            "overlay_statistics": {
                "detuning_points": int(len(aligned_x)),
                "valid_numeric_rows": valid_rows,
                "ignored_rows": ignored_rows,
                "merged_duplicates": merged_rows,
            },
            "warnings": ui_export.to_json_primitive(
                tuple(str(item).strip() for item in comparison_warnings if str(item).strip())
            ),
        }


def _skey(scheme_name, pname):
    return f"{scheme_name}__{pname}"


def _current_params(scheme_name, scheme_obj):
    """Live value of every param knob from session_state, falling back to defaults."""
    return {
        sp.name: st.session_state.get(_skey(scheme_name, sp.name), sp.default)
        for sp in scheme_obj.param_schema()
    }


def _render_group_header(container, group):
    container.markdown(
        f"<div class='gabes-group-header'>{escape(group)}</div>",
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
    choice_labels = getattr(sp, "choice_labels", None) or {}
    choice_kwargs = (
        {"format_func": lambda value: choice_labels.get(value, str(value))}
        if choice_labels else {}
    )
    on_change = None
    if getattr(sp, "applies_defaults", False) and scheme_obj is not None:
        on_change = lambda: _apply_recommended_defaults(scheme_name, scheme_obj, key)
    if getattr(sp, "control", "auto") == "checkbox":
        if has_state:
            return container.checkbox(label, key=key, help=help_,
                                      on_change=on_change)
        return container.checkbox(label, value=bool(current), key=key,
                                  help=help_, on_change=on_change)
    if getattr(sp, "control", "auto") == "segmented":
        options = list(sp.choices or ())
        if hasattr(container, "segmented_control"):
            try:
                return container.segmented_control(label, options, key=key,
                                                   help=help_, on_change=on_change,
                                                   **choice_kwargs)
            except TypeError:
                pass
        if has_state:
            return container.radio(label, options, key=key, help=help_,
                                   horizontal=True, on_change=on_change,
                                   **choice_kwargs)
        idx = options.index(current) if current in options else 0
        return container.radio(label, options, key=key, help=help_,
                               horizontal=True, index=idx, on_change=on_change,
                               **choice_kwargs)
    if sp.choices is not None:
        options = list(sp.choices)
        if has_state:
            return container.selectbox(label, options, key=key, help=help_,
                                       on_change=on_change, **choice_kwargs)
        idx = options.index(current) if current in options else 0
        return container.selectbox(label, options, index=idx, key=key, help=help_,
                                   on_change=on_change, **choice_kwargs)
    # Imported here, not at module level: tests lift this function into a
    # bare AppTest script (tests/test_fwm_excess_noise.py).
    from gabes_ui.format import slider_format as _slider_format
    slider_format = _slider_format(sp)
    slider_kwargs = {"format": slider_format} if slider_format else {}
    if has_state:
        val = container.slider(label, sp.vmin, sp.vmax, step=sp.step,
                               key=key, help=help_, **slider_kwargs)
    else:
        val = container.slider(label, sp.vmin, sp.vmax, value=current,
                               step=sp.step, key=key, help=help_, **slider_kwargs)
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


# ----------------------------------------------------------------------
# Shell pieces (docs/ui_redesign/decisions.md D1, D4)
# ----------------------------------------------------------------------
ADVANCED_FOLD = 6     # an advanced-only group longer than this starts folded


def _brand_html(variant):
    mark = base64.b64encode(ICON_SVG.encode("utf-8")).decode("ascii")
    return (f"<div class='gabes-brand {variant}'>"
            f"<img src='data:image/svg+xml;base64,{mark}' alt=''>"
            "<b>GABES</b><small>beta</small></div>")


def _render_regime_selector(scheme_name, sp, scheme_obj):
    """The scheme's regime/mode switch, shown in the top bar."""
    key = _skey(scheme_name, sp.name)
    last_key = _skey(scheme_name, f"_{sp.name}_last")
    options = list(sp.choices or ())
    labels = getattr(sp, "choice_labels", None) or {}
    if st.session_state.get(key) in options:
        st.session_state[last_key] = st.session_state[key]

    def _changed():
        # A segmented control can be clicked off; keep a regime selected.
        if st.session_state.get(key) not in options:
            st.session_state[key] = st.session_state.get(last_key, sp.default)
        if getattr(sp, "applies_defaults", False):
            _apply_recommended_defaults(scheme_name, scheme_obj, key)

    st.segmented_control(
        sp.label, options, key=key, help=sp.help or None, on_change=_changed,
        format_func=lambda value: clean_choice_label(labels.get(value, value)),
        label_visibility="collapsed",
    )


def _render_preset_selector(scheme_name, scheme_obj, sets):
    """Recommended default sets (e.g. OD / SAS) as a segmented control whose
    selection shows which set the current parameters match, if any."""
    key = _skey(scheme_name, "_preset")
    labels = list(sets)
    current = _current_params(scheme_name, scheme_obj)
    st.session_state[key] = next(
        (label for label in labels
         if all(current.get(k) == v for k, v in (sets[label] or {}).items())),
        None,
    )

    def _changed():
        picked = st.session_state.get(key)
        if picked is not None:
            _apply_default_set(scheme_name, scheme_obj, picked)

    st.segmented_control(
        "Defaults", labels, key=key, on_change=_changed,
        format_func=lambda label: clean_choice_label(label).removesuffix(" default"),
        help="Load a ready-made parameter set.", label_visibility="collapsed",
    )


def _apply_default_set(sname, sc, label):
    cur = _current_params(sname, sc)
    sets = sc.recommended_defaults(cur) or {}
    for k, v in (sets.get(label) or {}).items():
        st.session_state[_skey(sname, k)] = v


def _pop_figure_title(fig):
    """Take the title text off the image; the plot card shows it as a caption."""
    parts = []
    suptitle = getattr(fig, "_suptitle", None)
    if suptitle is not None and suptitle.get_text().strip():
        parts.append(suptitle.get_text())
        suptitle.set_text("")
    for ax in fig.axes:
        for loc in ("left", "center", "right"):
            text = ax.get_title(loc=loc)
            if text.strip():
                parts.append(text)
                ax.set_title("", loc=loc)
    return "  ·  ".join(parts)


def _render_plot_card(view, scheme, spec_by_name):
    """Plot views as tabs, overlay/export tools, figure. Returns (comparison
    payload, plot payloads, export popover) for the export step."""
    fig = view.get("figure")
    pages = [
        (item.get("label") or f"View {index + 1}", item.get("figure"))
        for index, item in enumerate(view.get("figure_views", []) or [])
        if isinstance(item, dict) and item.get("figure") is not None
    ]
    if not pages and fig is not None:
        pages = [("Plot", fig)]
    if not pages:
        return None, [], None

    labels = [label for label, _ in pages]
    selected = labels[0]
    card = st.container(border=True, key="gabes_plotcard")
    with card:
        head = st.container(horizontal=True, vertical_alignment="center",
                            gap="small", key="gabes_plothead")
    with head:
        if len(pages) > 1:
            view_key = _skey(scheme.name, "_plot_view")
            last_key = _skey(scheme.name, "_plot_view_last")
            if st.session_state.get(view_key) not in labels:
                st.session_state[view_key] = (
                    st.session_state.get(last_key)
                    if st.session_state.get(last_key) in labels else labels[0])

            def _keep_view(k=view_key, lk=last_key):
                if st.session_state.get(k) not in labels:
                    st.session_state[k] = st.session_state.get(lk, labels[0])

            selected = st.segmented_control(
                "Plot view", labels, key=view_key, on_change=_keep_view,
                label_visibility="collapsed") or labels[0]
            st.session_state[last_key] = selected
        caption_slot = st.empty()
        tools = st.container(horizontal=True, horizontal_alignment="right",
                             vertical_alignment="center", gap="small",
                             key="gabes_plottools")
    with tools:
        overlay_on = bool(view.get("comparison")) and st.toggle(
            "Overlay data", key=_skey(scheme.name, "_overlay"))
        export_box = st.popover("Export", icon=":material/download:", type="tertiary")

    with card:
        if overlay_on:
            plot_area, panel = st.columns([2.6, 1], gap="medium")
            comparison_payload = _render_experimental_comparison(view, scheme.name, panel)
        else:
            plot_area, comparison_payload = st.container(), None

    plot_payloads = ui_export.collect_view_plot_payloads(view, fig)
    page_fig = dict(pages)[selected]
    caption = _pop_figure_title(page_fig)
    with plot_area:
        for control_name in view.get("figure_controls", []):
            control_spec = spec_by_name.get(control_name)
            if control_spec is not None:
                _render_param(st, scheme.name, control_spec, scheme)
        _render_fig(page_fig)
    if caption:
        caption_slot.markdown(
            f"<div class='g-plot-caption' title='{escape(caption, quote=True)}'>"
            f"{caption_html(caption)}</div>", unsafe_allow_html=True)
    for _label, other in pages:
        if other is not page_fig:
            _close_fig(other)
    return comparison_payload, plot_payloads, export_box


def _render_more(view, scheme, params, cache_version):
    """Tables, diagnostic figures and heavy on-demand views behind pills."""
    items = [(table["title"], "table", table) for table in view.get("tables", [])]
    items += [(f"Diagnostic · {title}", "figure", extra_fig)
              for title, extra_fig in view.get("figures", [])]
    items += [(view_def.key, "extra", view_def) for view_def in scheme.extra_views()]
    if not items:
        return
    labels = [label for label, _, _ in items]
    more_key = _skey(scheme.name, "_more")
    if st.session_state.get(more_key) not in labels:
        st.session_state[more_key] = None
    with st.container(horizontal=True, vertical_alignment="center", gap="small",
                      key="gabes_more"):
        st.markdown("<span class='g-more-label'>More</span>",
                    unsafe_allow_html=True, width="content")
        picked = st.pills("More", labels, key=more_key, label_visibility="collapsed")
    for label, kind, obj in items:
        if kind == "figure" and label != picked:
            _close_fig(obj)
    if picked is None:
        return
    _label, kind, obj = items[labels.index(picked)]
    with st.container(border=True, key="gabes_morepanel"):
        if kind == "table":
            st.markdown(obj["markdown"])
        elif kind == "figure":
            _render_fig(obj)
        else:
            st.caption(obj.description)
            if st.button("Run", key=f"run__{scheme.name}__{obj.key}",
                         icon=":material/play_arrow:"):
                extra_keys = set(scheme.recompute_keys()) | set(obj.param_keys)
                extra_items = tuple(sorted(
                    (key, params[key]) for key in extra_keys if key in params))
                with st.spinner("Running…"):
                    data = _cached_extra(scheme.name, obj.key, extra_items, cache_version)
                _render_fig(obj.render(data, params) if obj.render_with_params
                            else obj.render(data))


# ----------------------------------------------------------------------
# Router — SABES is a separate site living in the same deployment
# ----------------------------------------------------------------------
# SABES is not a scheme: it is a lab-facing layer over one specific experiment,
# with primary variables (waveplate angles, generator frequency, lens choices)
# instead of processed physics parameters. Routing on a query parameter keeps it
# out of the scheme dropdown and off this page entirely, while still deploying as
# one app. The helpers defined above are handed over explicitly rather than
# imported back, so `sabes_page` never depends on this module's import state.
SABES_QUERY_VALUE = "sabes"


def _sabes_host():
    from types import SimpleNamespace
    return SimpleNamespace(render_fig=_render_fig, theme_base=THEME_BASE)


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
st.sidebar.markdown(_brand_html("gabes-brand--rail"), unsafe_allow_html=True)

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
    st.markdown(_brand_html("gabes-brand--top"), unsafe_allow_html=True, width="content")
    choice = st.selectbox("Scheme", titles, key="_scheme_choice",
                          label_visibility="collapsed", width=310)
scheme = all_schemes[titles.index(choice)]

specs = scheme.param_schema()
defaults_version = getattr(scheme, "defaults_version", "1")
defaults_key = _skey(scheme.name, "_defaults_version")
if st.session_state.get(defaults_key) != defaults_version:
    for sp in specs:
        st.session_state[_skey(scheme.name, sp.name)] = sp.default
    st.session_state[defaults_key] = defaults_version
else:
    for sp in specs:
        key = _skey(scheme.name, sp.name)
        # Preserve hidden controls through widget cleanup and send their saved
        # values to the browser when they become visible again.
        st.session_state[key] = st.session_state.get(key, sp.default)

# The regime/mode switch moves to the top bar: the first segmented control that
# applies a recommended default set (FWM Mode, Λ/Rydberg/magneto Regime).
regime_spec = next(
    (sp for sp in specs
     if getattr(sp, "applies_defaults", False)
     and getattr(sp, "control", "auto") == "segmented"
     and _param_visible(scheme.name, sp)),
    None,
)
# Schemes without such a switch offer their recommended sets (e.g. OD / SAS)
# there instead. Probed defensively so a scheme without the hook never breaks.
_rec_fn = getattr(scheme, "recommended_defaults", None)
_rec_sets = None
if callable(_rec_fn):
    try:
        _rec_sets = _rec_fn(_current_params(scheme.name, scheme))
    except Exception:
        _rec_sets = None
_mode_driven_defaults = any(getattr(sp, "applies_defaults", False) for sp in specs)

with topbar:
    about = st.popover("About", icon=":material/info:", type="tertiary",
                       help="What this scheme models, and its references.")
    if regime_spec is not None:
        _render_regime_selector(scheme.name, regime_spec, scheme)
    elif isinstance(_rec_sets, dict) and _rec_sets and not _mode_driven_defaults:
        _render_preset_selector(scheme.name, scheme, _rec_sets)
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
# Presets — one click overwrites the relevant sliders.
scheme_presets = scheme.presets()
if scheme_presets and getattr(scheme, "presets_group", None):
    _render_group_header(st.sidebar, scheme.presets_group)
for preset in scheme_presets:
    def _apply(p=preset, sname=scheme.name):
        for k, v in p.values.items():
            st.session_state[_skey(sname, k)] = v
    st.sidebar.button(preset.name, on_click=_apply,
                      use_container_width=True, help=preset.help)

# Grouped sections. "Show advanced" reveals advanced knobs inside the group they
# belong to (a pump waist next to the pump power); advanced-only groups follow,
# folded when long.
visible_specs = [sp for sp in specs
                 if _param_visible(scheme.name, sp) and sp is not regime_spec]
params = {}
group_order = []
for sp in visible_specs:
    if not sp.advanced and sp.group not in group_order:
        group_order.append(sp.group)
advanced = [sp for sp in visible_specs if sp.advanced]
advanced_by_group = {}
for sp in advanced:
    advanced_by_group.setdefault(getattr(sp, "advanced_group", "") or sp.group, []).append(sp)
advanced_key = _skey(scheme.name, "_show_advanced")
show_advanced = bool(st.session_state.get(advanced_key, False))

for g in group_order:
    _render_group_header(st.sidebar, g)
    for sp in visible_specs:
        if sp.group == g and not sp.advanced:
            params[sp.name] = _render_param(st.sidebar, scheme.name, sp, scheme)
    if show_advanced:
        for sp in advanced_by_group.pop(g, []):
            params[sp.name] = _render_param(st.sidebar, scheme.name, sp, scheme)

if advanced:
    def _toggle_advanced(k=advanced_key):
        st.session_state[k] = not st.session_state.get(k, False)

    st.sidebar.button(
        "Hide advanced" if show_advanced else f"Show advanced · {len(advanced)}",
        key=_skey(scheme.name, "_advanced_toggle"), on_click=_toggle_advanced,
        type="tertiary", icon=":material/tune:")
    if show_advanced:
        for group, group_specs in advanced_by_group.items():
            if len(group_specs) > ADVANCED_FOLD:
                box = st.sidebar.expander(f"{group} · {len(group_specs)}")
            else:
                _render_group_header(st.sidebar, group)
                box = st.sidebar
            for sp in group_specs:
                params[sp.name] = _render_param(box, scheme.name, sp, scheme)

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
    view = _cached_observables(
        scheme.name, raw, param_items,
        (cache_version, READOUT_CACHE_VERSION),
    )
else:
    with _PLOT_LOCK:
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
comparison_payload, plot_payloads, export_box = _render_plot_card(view, scheme, spec_by_name)

_render_more(view, scheme, params, cache_version)

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
