"""Result export: arrays, parameters and provenance behind one render.

Moved out of streamlit_app.py (docs/ui_redesign/plan.md P2). The bundle is
built only when a download is clicked (deferred `download_button` data), so
an ordinary rerun no longer serialises every raw array.
"""
import csv
import json
import re
import subprocess
import zipfile
from datetime import datetime
from io import BytesIO, StringIO
from pathlib import Path

import numpy as np
import streamlit as st

REPO_DIR = Path(__file__).resolve().parents[1]


def _safe_token(value, fallback="item"):
    token = str(value).strip() if value is not None else fallback
    token = re.sub(r"[^A-Za-z0-9._-]+", "-", token).strip(".-_")
    return token or fallback


def to_json_primitive(value):
    if value is None:
        return None
    if isinstance(value, (str, bool, int, float)):
        return value
    if isinstance(value, (np.integer, np.floating, np.bool_)):
        return value.item()
    if isinstance(value, (bytes, bytearray)):
        return value.hex()
    if isinstance(value, datetime):
        return value.isoformat()
    if isinstance(value, complex):
        return {
            "__complex__": True,
            "real": float(np.real(value)),
            "imag": float(np.imag(value)),
        }
    if isinstance(value, np.ndarray):
        if np.iscomplexobj(value):
            return {
                "__complex_array__": True,
                "real": to_json_primitive(np.real(value).tolist()),
                "imag": to_json_primitive(np.imag(value).tolist()),
            }
        return to_json_primitive(value.tolist())
    if isinstance(value, (list, tuple, set)):
        return [to_json_primitive(item) for item in value]
    if isinstance(value, dict):
        return {str(k): to_json_primitive(v) for k, v in value.items()}
    return str(value)


def _split_label_unit(label):
    text = str(label or "")
    text = text.strip()
    if text.endswith("]") and "[" in text:
        base, unit = text.rsplit("[", 1)
        return base.strip(), unit[:-1].strip()
    return text, ""


def _collect_figure_curve_payloads(fig):
    payload = {"traces": []}
    if not hasattr(fig, "axes"):
        return payload
    fig_title = getattr(fig, "_suptitle", None)
    payload["figure_title"] = (
        str(fig_title.get_text() if fig_title is not None else "").strip()
        or "Figure"
    )
    for axis_index, axis in enumerate(fig.axes):
        x_axis_label = axis.get_xlabel() or f"Axis {axis_index} X"
        y_axis_label = axis.get_ylabel() or f"Axis {axis_index} Y"
        x_name, x_unit = _split_label_unit(x_axis_label)
        y_name, y_unit = _split_label_unit(y_axis_label)
        for curve_index, curve in enumerate(axis.get_lines()):
            label = str(curve.get_label() or f"curve {curve_index + 1}")
            if label.startswith("_"):
                label = f"line {curve_index + 1}"
            x_data = np.asarray(curve.get_xdata())
            y_data = np.asarray(curve.get_ydata())
            if x_data.size == 0 or y_data.size == 0:
                continue
            n = min(len(x_data), len(y_data))
            payload["traces"].append({
                "axis_index": int(axis_index),
                "curve_index": int(curve_index),
                "figure_curve_label": label,
                "x_label": x_name or "x",
                "y_label": y_name or "y",
                "x_unit": x_unit,
                "y_unit": y_unit,
                "x_data": to_json_primitive(np.asarray(x_data[:n]).tolist()),
                "y_data": to_json_primitive(np.asarray(y_data[:n]).tolist()),
                "color": str(curve.get_color()),
            })
    return payload


def collect_view_plot_payloads(view, main_figure):
    payloads = []
    seen = set()

    def add_payload(label, figure_obj, source):
        if not hasattr(figure_obj, "axes"):
            return
        fid = id(figure_obj)
        if fid in seen:
            return
        seen.add(fid)
        payload = _collect_figure_curve_payloads(figure_obj)
        payload["label"] = str(label)
        payload["source"] = source
        payloads.append(payload)

    add_payload("primary", main_figure, "figure")
    for item in view.get("figure_views", ()):
        if isinstance(item, dict):
            add_payload(item.get("label", "view"), item.get("figure"), "figure_view")
    for figure_title, extra_figure in view.get("figures", ()):
        add_payload(figure_title, extra_figure, "diagnostic")
    return payloads


def _curve_csv_bytes(trace):
    x_data = trace.get("x_data", [])
    y_data = trace.get("y_data", [])
    n = min(len(x_data), len(y_data))
    x_unit = trace.get("x_unit", "")
    y_unit = trace.get("y_unit", "")
    x_name = trace.get("x_label", "x")
    y_name = trace.get("y_label", "y")
    x_hdr = f"{x_name}"
    y_hdr = f"{y_name}"
    if x_unit:
        x_hdr = f"{x_hdr} [{x_unit}]"
    if y_unit:
        y_hdr = f"{y_hdr} [{y_unit}]"
    buf = StringIO()
    writer = csv.writer(buf)
    writer.writerow([x_hdr, y_hdr])
    for idx in range(n):
        writer.writerow([x_data[idx], y_data[idx]])
    return buf.getvalue().encode("utf-8")


def _build_plot_zip_bytes(plot_payloads):
    if not plot_payloads:
        return None
    has_traces = any(
        bool(payload.get("traces")) for payload in plot_payloads if isinstance(payload, dict)
    )
    if not has_traces:
        return None
    zip_buffer = BytesIO()
    with zipfile.ZipFile(zip_buffer, "w", zipfile.ZIP_DEFLATED) as zf:
        for plot_index, plot_payload in enumerate(plot_payloads):
            if not isinstance(plot_payload, dict):
                continue
            plot_label = _safe_token(plot_payload.get("label", f"figure_{plot_index}"), f"figure{plot_index}")
            for trace_index, trace in enumerate(plot_payload.get("traces", ())):
                if not isinstance(trace, dict):
                    continue
                trace_label = _safe_token(
                    trace.get("figure_curve_label", f"trace{trace_index}"),
                    f"trace{trace_index}")
                curve_name = f"{plot_label}/{trace_label}.csv"
                zf.writestr(curve_name, _curve_csv_bytes(trace))
    return zip_buffer.getvalue()


@st.cache_data(show_spinner=False, max_entries=1)
def app_revision():
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "HEAD"],
            cwd=str(REPO_DIR),
            stderr=subprocess.STDOUT,
            text=True,
        ).strip()
    except Exception:
        return "unavailable"


def build_export_payload(*, scheme, raw, params, view, plot_payloads,
                         comparison_payload, cache_version, readout_cache_version):
    return {
        "generated_at_utc": datetime.utcnow().isoformat() + "Z",
        "app_revision": app_revision(),
        "scheme": {
            "name": scheme.name,
            "title": scheme.title,
            "cluster": scheme.cluster,
            "cache_version": cache_version,
            "defaults_version": scheme.defaults_version,
            "readout_cache_version": readout_cache_version,
        },
        "defaults": to_json_primitive(scheme.defaults()),
        "params": to_json_primitive(params),
        "raw": to_json_primitive(raw),
        "view": {
            "metrics": to_json_primitive(view.get("metrics", [])),
            "hero_count": view.get("hero_count", None),
            "tables": to_json_primitive(view.get("tables", [])),
            "figure_views": [
                {"label": item.get("label"), "source": item.get("source")}
                for item in plot_payloads
                if isinstance(item, dict)
            ],
            "comparison": to_json_primitive(comparison_payload),
        },
        "plots": to_json_primitive(plot_payloads),
    }


def render_export(container, *, scheme, raw, params, view, plot_payloads,
                  comparison_payload, cache_version, readout_cache_version):
    """Download buttons for the bundle, all curves, and each curve."""
    def bundle_bytes():
        payload = build_export_payload(
            scheme=scheme, raw=raw, params=params, view=view,
            plot_payloads=plot_payloads, comparison_payload=comparison_payload,
            cache_version=cache_version, readout_cache_version=readout_cache_version,
        )
        return json.dumps(payload, indent=2, ensure_ascii=False).encode("utf-8")

    timestamp = datetime.utcnow().strftime("%Y%m%dT%H%M%SZ")
    scheme_token = _safe_token(scheme.name or "gabes", "gabes")
    with container:
        st.caption("Arrays, parameters and provenance behind this render.")
        st.download_button(
            "Result bundle (JSON)",
            data=bundle_bytes,
            file_name=f"{scheme_token}_result_{timestamp}.json",
            mime="application/json",
            key=f"download_bundle_{scheme_token}",
            on_click="ignore",
            icon=":material/data_object:",
            width="stretch",
        )
        has_traces = any(p.get("traces") for p in plot_payloads if isinstance(p, dict))
        if has_traces:
            st.download_button(
                "All plotted curves (CSV zip)",
                data=lambda: _build_plot_zip_bytes(plot_payloads),
                file_name=f"{scheme_token}_plots_{timestamp}.zip",
                mime="application/zip",
                key=f"download_zip_{scheme_token}",
                on_click="ignore",
                icon=":material/folder_zip:",
                width="stretch",
            )
        export_curves = [
            (plot_index, trace_index, trace)
            for plot_index, plot_payload in enumerate(plot_payloads)
            for trace_index, trace in enumerate(plot_payload.get("traces", ()))
        ]
        if export_curves:
            st.caption("Single curves (CSV)")
            for plot_index, trace_index, curve in export_curves:
                plot_label = plot_payloads[plot_index].get("label", f"Figure {plot_index + 1}")
                trace_label = curve.get("figure_curve_label", f"trace{trace_index}")
                key = _safe_token(f"{scheme_token}_{plot_index}_{trace_index}")
                st.download_button(
                    f"{plot_label}: {trace_label}",
                    data=_curve_csv_bytes(curve),
                    file_name=f"{_safe_token(plot_label)}__{_safe_token(trace_label)}.csv",
                    mime="text/csv",
                    key=f"download_curve_{key}",
                    on_click="ignore",
                    type="tertiary",
                )
