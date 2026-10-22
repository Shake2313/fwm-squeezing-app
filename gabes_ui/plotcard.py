"""Plot card — the figure and everything that hangs off it.

Views as underline tabs (only the selected figure is drawn), the popped figure
title as a caption, experimental CSV overlay in a side panel, and the tables,
diagnostics and heavy extra views behind the More pills
(docs/ui_redesign/decisions.md D8).

The matplotlib lock lives here because this module owns every figure the app
builds; `streamlit_app.py` borrows it for the one call it makes itself.
"""
import copy
import hashlib
import importlib
import inspect
from html import escape
from threading import RLock

import matplotlib.pyplot as plt
import numpy as np
import streamlit as st

import gabes as _gabes
import gabes.experimental_csv as _experimental_csv
from gabes import schemes
from gabes.core import blas_single_thread
from gabes.experimental_alignment import AxisAlignment, frame_alignment, fit_trace_alignment
from gabes.plot_style import PALETTE, apply_gabes_plot_style
from gabes_ui import export as ui_export
from gabes_ui.controls import render_param, skey
from gabes_ui.readout import caption_html

PLOT_LOCK = RLock()

# Temporary user preference; restore the existing tabs by setting this True.
SHOW_SAS_EXTRA_VIEWS = False

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


@st.cache_data(show_spinner=False, max_entries=16)
def _cached_extra(scheme_name, view_key, param_items, cache_version):
    scheme = schemes.get(scheme_name)
    view = next(v for v in scheme.extra_views() if v.key == view_key)
    with blas_single_thread():
        return view.compute(dict(param_items))


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
    # Streamlit can retain a cached wrapper across a source reload while
    # replacing the module. Resolve the current class factory before pickling.
    return importlib.import_module("gabes.experimental_csv").load_experimental_csv(
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


def render_fig(fig):
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
    # Relative extrema data span 0–1 even when the theory's absorption is weak.
    # Include that range so the automatic correction cannot hide the CSV floor.
    ylim = (min(ylim[0], -0.05), max(ylim[1], 1.05))
    x_unit = descriptor.get("x_unit", "plot unit")
    raw_x_unit = descriptor.get("raw_x_unit", "Arb. unit")
    raw_y_unit = descriptor.get("raw_y_unit", "Arb. unit")

    with panel:
        panel.caption(
            "CSV A = X (time/frequency), B = detector signal. Metadata and "
            "later columns are ignored."
        )
        uploader_options = {}
        if "max_upload_size" in inspect.signature(panel.file_uploader).parameters:
            # Streamlit added this per-widget guard after the project's oldest
            # supported release; retain the backend byte check as the fallback.
            uploader_options["max_upload_size"] = MAX_FILE_BYTES // (1024 * 1024)
        uploaded = panel.file_uploader(
            "Oscilloscope CSV",
            type=("csv",),
            key=skey(scheme_name, "_csv_upload"),
            help="A and B are read as arbitrary units; no header is required.",
            **uploader_options,
        )
        auto_correct = panel.checkbox(
            "Automatic noise correction",
            value=True,
            key=skey(scheme_name, "_csv_auto_correct"),
            help="Merge repeated A values, reject isolated spikes, and apply "
                 "noise-adaptive local smoothing before the selected calibration.",
        )
        show_raw = panel.checkbox(
            "Show unfiltered trace",
            value=False,
            disabled=not auto_correct,
            key=skey(scheme_name, "_csv_show_raw"),
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
            key=skey(scheme_name, "_csv_calibration_mode"),
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
                key=skey(scheme_name, "_csv_dark_signal"),
                help="Detector signal measured at zero transmitted light.",
            )
            reference_signal = calibration_cols[1].number_input(
                f"Reference signal (T = 1)  [{raw_y_unit}]",
                value=1.0,
                key=skey(scheme_name, "_csv_reference_signal"),
                help="Detector signal measured for unit-transmission reference light.",
            )
        elif calibration_mode == CALIBRATION_ABSOLUTE_GAIN_OFFSET:
            calibration_cols = panel.columns(2)
            gain = calibration_cols[0].number_input(
                f"Detector gain  [{raw_y_unit}]",
                value=1.0,
                key=skey(scheme_name, "_csv_detector_gain"),
                help="Gain in signal = offset + gain × transmission; a negative "
                     "gain represents inverted detector polarity.",
            )
            offset = calibration_cols[1].number_input(
                f"Detector offset  [{raw_y_unit}]",
                value=0.0,
                key=skey(scheme_name, "_csv_detector_offset"),
                help="Signal value corresponding to zero transmission.",
            )
        if uploaded is None:
            return

        csv_bytes = uploaded.getvalue()
        fingerprint = hashlib.sha256(csv_bytes).hexdigest()
        fingerprint_key = skey(scheme_name, "_csv_fingerprint")
        scale_key = skey(scheme_name, "_csv_x_scale")
        multiplier_key = skey(scheme_name, "_csv_x_multiplier")
        shift_key = skey(scheme_name, "_csv_x_shift")
        reverse_key = skey(scheme_name, "_csv_reverse")
        invert_key = skey(scheme_name, "_csv_invert")
        signal_scale_key = skey(scheme_name, "_csv_y_scale")
        signal_offset_key = skey(scheme_name, "_csv_y_offset")
        fit_key = skey(scheme_name, "_csv_alignment_fit")
        alignment_identity = (
            fingerprint,
            axis_index,
            str(x_unit),
            str(raw_x_unit),
        )
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
        plot_span = float(xlim[1] - xlim[0])
        try:
            framed = frame_alignment(detuning, xlim)
        except ValueError as exc:
            panel.error(str(exc))
            return
        framed_scale, framed_shift = framed.scale, framed.shift

        def set_alignment(alignment, *, invert=False):
            st.session_state[scale_key] = float(alignment.scale)
            st.session_state[multiplier_key] = 1.0
            st.session_state[shift_key] = float(alignment.shift)
            st.session_state[reverse_key] = bool(alignment.reverse)
            st.session_state[invert_key] = bool(invert)
            st.session_state[signal_scale_key] = float(alignment.signal_scale)
            st.session_state[signal_offset_key] = float(alignment.signal_offset)
            st.session_state.pop(fit_key, None)

        if st.session_state.get(fingerprint_key) != alignment_identity:
            st.session_state[fingerprint_key] = alignment_identity
            set_alignment(framed)
        if trace.is_absolute_transmission:
            st.session_state[invert_key] = False
            st.session_state[signal_scale_key] = 1.0
            st.session_state[signal_offset_key] = 0.0

        action_cols = panel.columns(2)
        if action_cols[0].button(
            "Bring into view",
            key=skey(scheme_name, "_csv_frame"),
            use_container_width=True,
            help="Map the imported span into the current theoretical plot; "
                 "this does not fit spectral features.",
        ):
            try:
                reframed = AxisAlignment(
                    scale=framed.scale, shift=framed.shift,
                    reverse=st.session_state[reverse_key],
                    signal_scale=st.session_state[signal_scale_key],
                    signal_offset=st.session_state[signal_offset_key],
                )
                set_alignment(reframed, invert=st.session_state[invert_key])
            except ValueError as exc:
                panel.error(str(exc))
        if action_cols[1].button(
            "Use raw axis",
            key=skey(scheme_name, "_csv_reset"),
            use_container_width=True,
        ):
            try:
                raw_axis = AxisAlignment(
                    signal_scale=st.session_state[signal_scale_key],
                    signal_offset=st.session_state[signal_offset_key],
                )
                set_alignment(raw_axis, invert=st.session_state[invert_key])
            except ValueError as exc:
                panel.error(str(exc))
        # The first curve is the scheme's transmission; later lines are markers.
        theory = axis.lines[0] if axis.lines else None
        fit_context_key = skey(scheme_name, "_csv_fit_context")
        fit_context = (
            fingerprint, auto_correct, calibration_mode,
            dark_signal, reference_signal, gain, offset,
            None if theory is None else hashlib.sha256(
                np.asarray([theory.get_xdata(), theory.get_ydata()]).tobytes()
            ).hexdigest(),
        )
        if st.session_state.get(fit_context_key) != fit_context:
            st.session_state.pop(fit_key, None)
            st.session_state[fit_context_key] = fit_context
        if panel.button(
            "Fit alignment", key=skey(scheme_name, "_csv_fit"),
            use_container_width=True, disabled=theory is None,
            help="Align X scale/shift to the current theoretical transmission. "
                 "For relative data, also fit detector gain/offset. Physical "
                 "parameters are held fixed; choose sweep direction first.",
        ):
            try:
                initial = AxisAlignment(
                    scale=st.session_state[scale_key] * st.session_state[multiplier_key],
                    shift=st.session_state[shift_key],
                    reverse=st.session_state[reverse_key],
                    signal_scale=st.session_state[signal_scale_key],
                    signal_offset=st.session_state[signal_offset_key],
                )
                result = fit_trace_alignment(
                    trace, np.asarray(theory.get_xdata()), np.asarray(theory.get_ydata()),
                    initial=initial, invert=st.session_state[invert_key],
                )
                if result.success:
                    set_alignment(result.alignment, invert=st.session_state[invert_key])
                st.session_state[fit_key] = result
            except ValueError as exc:
                panel.error(str(exc))

        panel.caption(
            "X span is framed automatically; this is not a frequency calibration. "
            "Drag the curve or use the controls below the plot for a live preview. "
            "The precise controls here take effect on Apply alignment."
        )
        # Forms keep slider drags in the browser and send one update on Apply.
        precise = panel.expander("Precise alignment · base scale, direction and signal", expanded=False)
        with precise.form(skey(scheme_name, "_csv_alignment_form")):
            base_scale = st.number_input(
                f"Base X scale  [{x_unit}/{raw_x_unit}]", min_value=0.0,
                key=scale_key, format="%.6g",
                step=float(st.session_state[scale_key] * 0.01),
                help="Automatic coarse scale; scientific notation is accepted.",
            )
            multiplier = st.slider(
                "X fine scale ×", 0.1, 3.0, step=0.005, key=multiplier_key,
                help="Stable relative range around the base scale.",
            )
            shift_now = float(st.session_state[shift_key])
            shift_margin = 2.0 * abs(plot_span)
            x_shift = st.slider(
                f"X shift  [{x_unit}]",
                min(xlim[0] - shift_margin, shift_now),
                max(xlim[1] + shift_margin, shift_now),
                step=abs(plot_span) / 1000.0, key=shift_key,
            )
            reverse = st.checkbox("Reverse sweep", key=reverse_key)
            invert = st.checkbox(
                "Invert transmission", key=invert_key,
                disabled=trace.is_absolute_transmission,
                help="Absolute detector polarity is defined by the calibration gain.",
            )
            signal_scale = st.number_input(
                "Y scale ×", min_value=0.0, step=0.01, key=signal_scale_key,
                disabled=trace.is_absolute_transmission,
            )
            signal_offset = st.number_input(
                "Y offset", step=0.01, key=signal_offset_key,
                disabled=trace.is_absolute_transmission,
            )
            applied = st.form_submit_button("Apply alignment", use_container_width=True)
        if applied:
            st.session_state.pop(fit_key, None)
        x_scale = float(base_scale) * float(multiplier)
        try:
            alignment = AxisAlignment(
                scale=x_scale, shift=float(x_shift), reverse=bool(reverse),
                signal_scale=float(signal_scale), signal_offset=float(signal_offset),
            )
            aligned_x = alignment.transform_x(trace.detuning)
            aligned_y = alignment.transform_y(trace.transmission, invert=bool(invert))
        except ValueError as exc:
            panel.error(str(exc))
            return
        fit_result = st.session_state.get(fit_key)
        if fit_result is not None:
            if fit_result.success:
                panel.caption(
                    f"Alignment RMSE {fit_result.rmse:.4g} · "
                    f"overlap {fit_result.overlap_fraction:.0%}. "
                    "Temperature and isotope composition were held fixed."
                )
            else:
                panel.warning(fit_result.message)

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

        in_view = ((aligned_x >= xlim[0]) & (aligned_x <= xlim[1])).any()
        if not in_view:
            panel.warning(
                "The imported trace is outside the theoretical x range. "
                "Use Bring into view, then refine X scale and X shift manually."
            )

        raw_overlay = None
        if auto_correct and show_raw:
            raw_x = aligned_x
            raw_y = trace.calibrate_signal(trace.raw_signal)
            try:
                raw_y = alignment.transform_y(raw_y, invert=bool(invert))
            except ValueError as exc:
                panel.error(str(exc))
                return
            raw_overlay = (raw_x, raw_y)

        overlay_label = (
            f"{descriptor.get('label', 'Experimental CSV')} — "
            f"{trace.transmission_label.lower()}"
        )

    with PLOT_LOCK:
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
                "base_scale": float(base_scale),
                "manual_multiplier": float(multiplier),
                "signal_scale": float(signal_scale),
                "signal_offset": float(signal_offset),
            },
            "alignment_fit": None if fit_result is None else {
                "success": bool(fit_result.success),
                "rmse": float(fit_result.rmse),
                "initial_rmse": float(fit_result.initial_rmse),
                "overlap_fraction": float(fit_result.overlap_fraction),
                "points_used": int(fit_result.points_used),
                "message": str(fit_result.message),
                "physical_parameters_fitted": False,
            },
            "show_raw": bool(show_raw),
            "in_view": bool(in_view),
            "framed_scale": float(framed_scale),
            "framed_shift": float(framed_shift),
            "comparison_trace": {
                "detuning": ui_export.to_json_primitive(aligned_x),
                "transmission": ui_export.to_json_primitive(aligned_y),
                "detuning_unit": str(x_unit),
                "transmission_unit": "1",
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


def render_plot_card(view, scheme, spec_by_name, raw=None):
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

    # Temporary display policy: keep all physical readouts/export data intact.
    hidden_pages = []
    if scheme.name == "sas" and not SHOW_SAS_EXTRA_VIEWS:
        hidden_pages = [f for label, f in pages if label != "Transmission"]
        pages = [(label, f) for label, f in pages if label == "Transmission"] or [("Transmission", fig)]

    labels = [label for label, _ in pages]
    selected = labels[0]
    card = st.container(border=True, key="gabes_plotcard")
    with card:
        head = st.container(horizontal=True, vertical_alignment="center",
                            gap="small", key="gabes_plothead")
    with head:
        if len(pages) > 1:
            view_key = skey(scheme.name, "_plot_view")
            last_key = skey(scheme.name, "_plot_view_last")
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
            "Overlay data", key=skey(scheme.name, "_overlay"))
        export_box = st.popover("Export", icon=":material/download:", type="tertiary")

    with card:
        plot_area = st.container()
        if overlay_on:
            panel = st.expander("Overlay settings · upload, calibrate and fit", expanded=True)
            comparison_payload = _render_experimental_comparison(view, scheme.name, panel)
        else:
            comparison_payload = None

    plot_payloads = ui_export.collect_view_plot_payloads(view, fig)
    page_fig = dict(pages)[selected]
    caption = _pop_figure_title(page_fig)
    with plot_area:
        for control_name in view.get("figure_controls", []):
            control_spec = spec_by_name.get(control_name)
            if control_spec is not None:
                render_param(st, scheme.name, control_spec, scheme)
        from gabes_ui.liveplot import render_live_plot
        with PLOT_LOCK:
            live = False
            if not (scheme.name == "fwm" and raw and raw.get("kind") == "biphoton"):
                live = render_live_plot(page_fig, scheme.name, comparison_payload,
                                        canonical_noise=bool(raw and "T_canonical_small_signal" in raw),
                                        gain_digits=6 if (raw or {}).get("claim_gate", {}).get("quantitative_gain_supported") else 3)
        if live:
            _close_fig(page_fig)
        else:
            render_fig(page_fig)
    if caption:
        caption_slot.markdown(
            f"<div class='g-plot-caption' title='{escape(caption, quote=True)}'>"
            f"{caption_html(caption)}</div>", unsafe_allow_html=True)
    for _label, other in pages:
        if other is not page_fig:
            _close_fig(other)
    for other in hidden_pages:
        if other is not page_fig:
            _close_fig(other)
    return comparison_payload, plot_payloads, export_box


def _render_plot_and_export(view, scheme, spec_by_name, raw, params,
                            cache_version, readout_cache_version):
    comparison, plots, export_box = render_plot_card(view, scheme, spec_by_name, raw)
    if export_box is not None:
        ui_export.render_export(
            export_box, scheme=scheme, raw=raw, params=params, view=view,
            plot_payloads=plots, comparison_payload=comparison,
            cache_version=cache_version, readout_cache_version=readout_cache_version,
        )


@st.fragment
def _render_comparison_fragment(view, scheme, spec_by_name, raw, params,
                                 cache_version, readout_cache_version):
    # Fragment reruns reuse their input objects. Clone only the plot figures so
    # every rerun starts with clean theory, without accumulated CSV lines/titles.
    with PLOT_LOCK:
        local_view = dict(view)
        local_view.update(copy.deepcopy({
            "figure": view.get("figure"),
            "figure_views": view.get("figure_views", []),
        }))
        _render_plot_and_export(
            local_view, scheme, spec_by_name, raw, params,
            cache_version, readout_cache_version,
        )
        # Closed figures retain their data and can still seed the next clone.
        for figure in [view.get("figure"), *(
                item.get("figure") for item in view.get("figure_views", []))]:
            if figure is not None:
                _close_fig(figure)


def render_plot_section(view, scheme, spec_by_name, raw, params,
                        cache_version, readout_cache_version):
    """Keep CSV and export reruns local; physics controls retain full reruns."""
    renderer = (_render_comparison_fragment
                if view.get("comparison") and not view.get("figure_controls")
                else _render_plot_and_export)
    renderer(view, scheme, spec_by_name, raw, params,
             cache_version, readout_cache_version)


def render_more(view, scheme, params, cache_version):
    """Tables, diagnostic figures and heavy on-demand views behind pills."""
    items = [(table["title"], "table", table) for table in view.get("tables", [])]
    items += [(f"Diagnostic · {title}", "figure", extra_fig)
              for title, extra_fig in view.get("figures", [])]
    items += [(view_def.key, "extra", view_def) for view_def in scheme.extra_views()]
    if not items:
        return
    labels = [label for label, _, _ in items]
    more_key = skey(scheme.name, "_more")
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
            if scheme.name == "fwm" and obj["title"] == "Model diagnostics":
                from gabes_ui.diagnostics import render_diagnostics
                render_diagnostics(obj["markdown"])
            else:
                st.markdown(obj["markdown"])
        elif kind == "figure":
            render_fig(obj)
        else:
            st.caption(obj.description)
            if st.button("Run", key=f"run__{scheme.name}__{obj.key}",
                         icon=":material/play_arrow:"):
                extra_keys = set(scheme.recompute_keys()) | set(obj.param_keys)
                extra_items = tuple(sorted(
                    (key, params[key]) for key in extra_keys if key in params))
                with st.spinner("Running…"):
                    data = _cached_extra(scheme.name, obj.key, extra_items, cache_version)
                render_fig(obj.render(data, params) if obj.render_with_params
                            else obj.render(data))
