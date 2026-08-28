"""Config-driven Rydberg cell-heating simulation and report artifact workflow.

Run from the repository root:

    python -m analysis.rydberg_cell_heating.workflow \
        --config analysis/rydberg_cell_heating/example_config.json

The workflow distinguishes the implemented first-order finite-IF weak-signal
response from the deliberately deferred full time-domain LO+SIG/lock-in model.
Missing raw data and unavailable optional physics are emitted as ``PENDING``,
never silently substituted.
"""
from __future__ import annotations

import argparse
import csv
import datetime as dt
from dataclasses import asdict
import hashlib
import json
import math
import os
import platform
import subprocess
import sys
from pathlib import Path
from typing import Any, Iterable, Mapping

import numpy as np

try:  # Supports both ``python -m ...`` and direct script execution.
    from .adapters import (
        ALLOWED_STATUSES,
        ElectrometryAdapter,
        OptionalHelperAdapter,
        StaticRydbergAdapter,
        validate_status,
    )
except ImportError:  # pragma: no cover - direct-script convenience path
    ROOT = Path(__file__).resolve().parents[2]
    sys.path.insert(0, str(ROOT))
    from analysis.rydberg_cell_heating.adapters import (  # type: ignore
        ALLOWED_STATUSES,
        ElectrometryAdapter,
        OptionalHelperAdapter,
        StaticRydbergAdapter,
        validate_status,
    )

from gabes.rydberg_experiment import (
    AxialCellProfile,
    axial_cell_profile,
    effective_atom_number,
    integrate_beer_lambert,
    linear_axial_cell_profile,
    resolve_cell_temperature,
    sam_field_calibration,
)
from gabes.rydberg_experimental_csv import load_eit_csv, load_psd_csv, load_rf_sweep_csv


SCHEMA_VERSION = 1
ROOT = Path(__file__).resolve().parents[2]

CSV_SCHEMAS: dict[str, dict[str, set[str]]] = {
    "sensitivity": {
        "required": {"temperature_c", "sensitivity_nv_cm_sqrt_hz"},
        "numeric": {
            "temperature_c", "actual_temperature_c",
            "sensitivity_nv_cm_sqrt_hz",
            "uncertainty_nv_cm_sqrt_hz", "psn_limit_nv_cm_sqrt_hz",
        },
    },
    "temperature_log": {
        "required": {"setpoint_c"},
        "numeric": {
            "setpoint_c", "sensor_left_c", "sensor_right_c", "cold_spot_c",
            "effective_vapor_temp_c", "elapsed_s",
        },
    },
}
TRACE_KINDS = frozenset({"eit_spectrum", "rf_sweep", "psd"})
SUPPORTED_INPUT_KINDS = frozenset(CSV_SCHEMAS) | TRACE_KINDS


def _utc_now() -> str:
    return dt.datetime.now(dt.timezone.utc).replace(microsecond=0).isoformat()


def _json_dump(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(payload, indent=2, ensure_ascii=False, allow_nan=False) + "\n",
        encoding="utf-8",
    )


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _serialized_path(path: str | os.PathLike[str]) -> str:
    """Return a portable, non-identifying path for JSON/report provenance.

    Repository-owned files are written as POSIX paths relative to ``ROOT``.
    External files retain a readable basename plus a short hash of their
    resolved location, which lets repeated references be matched without
    publishing drive letters, user directories, or temporary-directory names.
    The real :class:`Path` is still used for all filesystem operations.
    """

    resolved = Path(path).expanduser().resolve()
    try:
        return resolved.relative_to(ROOT).as_posix()
    except ValueError:
        location_digest = hashlib.sha256(
            str(resolved).encode("utf-8", errors="surrogatepass")
        ).hexdigest()[:16]
        basename = resolved.name or "external-path"
        safe_basename = "".join(
            character if character.isalnum() or character in "._-" else "_"
            for character in basename
        )
        return f"external://{safe_basename}#path-sha256={location_digest}"


def _resolve(base: Path, value: str | os.PathLike[str]) -> Path:
    path = Path(value).expanduser()
    return path.resolve() if path.is_absolute() else (base / path).resolve()


def _numeric(value: str, *, column: str, row_number: int, path: Path) -> float | None:
    text = value.strip()
    if not text:
        return None
    try:
        number = float(text)
    except ValueError as exc:
        raise ValueError(
            f"{path}: row {row_number}, column {column!r} is not numeric: {text!r}"
        ) from exc
    if not math.isfinite(number):
        raise ValueError(
            f"{path}: row {row_number}, column {column!r} must be finite")
    return number


def load_clean_csv(path: Path, kind: str, default_status: str) -> list[dict[str, Any]]:
    """Load one strict, header-based CSV without guessing column meanings."""
    if kind not in CSV_SCHEMAS:
        raise ValueError(f"unsupported input kind {kind!r}")
    schema = CSV_SCHEMAS[kind]
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        headers = set(reader.fieldnames or [])
        missing = schema["required"] - headers
        if missing:
            raise ValueError(f"{path}: missing required columns {sorted(missing)}")
        rows: list[dict[str, Any]] = []
        for row_number, row in enumerate(reader, start=2):
            if not any((value or "").strip() for value in row.values()):
                continue
            parsed: dict[str, Any] = {}
            for column, value in row.items():
                if column is None:
                    continue
                text = (value or "").strip()
                if column in schema["numeric"]:
                    parsed[column] = _numeric(
                        text, column=column, row_number=row_number, path=path)
                elif column == "status":
                    parsed[column] = validate_status(text or default_status)
                else:
                    parsed[column] = text
            parsed.setdefault("status", default_status)
            rows.append(parsed)
    if not rows:
        raise ValueError(f"{path}: no data rows")
    return rows


def _trace_provenance(
    provenance: Any, source_path: Path | None = None
) -> dict[str, Any]:
    """Convert the immutable experimental-loader provenance to JSON values."""
    result = asdict(provenance)
    source_name = str(result.get("source_name", ""))
    if source_name:
        named_path = Path(source_name).expanduser()
        if named_path.is_absolute():
            result["source_name"] = _serialized_path(named_path)
        elif source_path is not None and source_name == str(source_path):
            result["source_name"] = _serialized_path(source_path)
    result["header"] = list(result["header"]) if result["header"] is not None else None
    result["assumptions"] = list(result["assumptions"])
    return result


def _load_trace_input(
    path: Path,
    kind: str,
    spec: Mapping[str, Any],
    status: str,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    """Load a raw trace through the shared GABES experimental CSV layer."""
    options = dict(spec.get("loader", {}))
    options.setdefault("source_name", _serialized_path(path))
    if kind == "eit_spectrum":
        if spec.get("temperature_c") is None:
            raise ValueError(
                f"EIT input {spec.get('id')!r} requires temperature_c metadata")
        temperature = float(spec["temperature_c"])
        trace = load_eit_csv(path, **options)
        rows = [
            {
                "temperature_c": temperature,
                "detuning_mhz": float(x),
                "signal": float(y),
                "status": status,
            }
            for x, y in zip(trace.detuning_mhz, trace.signal)
        ]
        metadata = {
            "temperature_c": temperature,
            "signal_unit": trace.signal_unit,
        }
    elif kind == "rf_sweep":
        trace = load_rf_sweep_csv(path, **options)
        rows = [
            {"drive": float(x), "response": float(y), "status": status}
            for x, y in zip(trace.drive, trace.response)
        ]
        metadata = {
            "temperature_c": spec.get("temperature_c"),
            "drive_quantity": trace.drive_quantity,
            "drive_unit": trace.drive_unit,
            "response_unit": trace.response_unit,
        }
    elif kind == "psd":
        trace = load_psd_csv(path, **options)
        rows = [
            {
                "frequency_hz": float(frequency),
                "asd": float(asd),
                "psd": float(psd),
                "status": status,
            }
            for frequency, asd, psd in zip(trace.frequency_hz, trace.asd, trace.psd)
        ]
        metadata = {
            "temperature_c": spec.get("temperature_c"),
            "quantity": trace.quantity,
            "base_unit": trace.base_unit,
            "input_spectrum_kind": trace.input_spectrum_kind,
        }
    else:  # pragma: no cover - guarded by caller
        raise ValueError(f"unsupported trace kind {kind!r}")
    metadata["parser_provenance"] = _trace_provenance(
        trace.provenance, source_path=path
    )
    return rows, metadata


def load_inputs(
    specs: Iterable[Mapping[str, Any]], config_dir: Path
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    """Load configured datasets and return results plus provenance entries."""
    datasets: list[dict[str, Any]] = []
    provenance: list[dict[str, Any]] = []
    seen_ids: set[str] = set()
    for index, source in enumerate(specs):
        spec = dict(source)
        dataset_id = str(spec.get("id", f"input_{index + 1}")).strip()
        if not dataset_id or dataset_id in seen_ids:
            raise ValueError(f"input ids must be non-empty and unique: {dataset_id!r}")
        seen_ids.add(dataset_id)
        kind = str(spec.get("kind", "")).strip()
        if kind not in SUPPORTED_INPUT_KINDS:
            raise ValueError(f"input {dataset_id!r} has unsupported kind {kind!r}")
        declared_status = validate_status(spec.get("status", "PENDING"))
        required = bool(spec.get("required", False))
        path_value = spec.get("path")
        base = {
            "id": dataset_id,
            "kind": kind,
            "declared_status": declared_status,
            "required": required,
        }
        if not path_value:
            if required:
                raise FileNotFoundError(f"required input {dataset_id!r} has no path")
            datasets.append({
                **base,
                "state": "not_configured",
                "status": "PENDING",
                "rows": [],
                "note": str(spec.get("note", "raw data not configured")),
            })
            provenance.append({**base, "state": "not_configured", "status": "PENDING"})
            continue
        path = _resolve(config_dir, str(path_value))
        serialized_path = _serialized_path(path)
        if not path.is_file():
            if required:
                raise FileNotFoundError(f"required input not found: {path}")
            datasets.append({
                **base,
                "state": "missing_optional",
                "status": "PENDING",
                "rows": [],
                "path": serialized_path,
                "note": str(spec.get("note", "optional raw data file not found")),
            })
            provenance.append({
                **base, "state": "missing_optional", "status": "PENDING",
                "path": serialized_path,
            })
            continue
        if kind in TRACE_KINDS:
            rows, trace_metadata = _load_trace_input(
                path, kind, spec, declared_status)
        else:
            rows = load_clean_csv(path, kind, declared_status)
            trace_metadata = {}
        datasets.append({
            **base,
            "state": "loaded",
            "status": declared_status,
            "path": serialized_path,
            "rows": rows,
            **trace_metadata,
        })
        provenance.append({
            **base,
            "state": "loaded",
            "status": declared_status,
            "path": serialized_path,
            "sha256": _sha256(path),
            "bytes": path.stat().st_size,
            "row_count": len(rows),
            **({"parser": trace_metadata.get("parser_provenance")}
               if trace_metadata else {}),
        })
    return datasets, provenance


def _metric_value(point: Mapping[str, Any], name: str) -> float | None:
    value = point.get("metrics", {}).get(name, {}).get("value")
    if value is None:
        return None
    number = float(value)
    return number if math.isfinite(number) else None


def _best_point(
    points: list[dict[str, Any]], metric: str, *, minimize: bool = False
) -> dict[str, Any] | None:
    candidates = [p for p in points if _metric_value(p, metric) is not None]
    if not candidates:
        return None
    selector = min if minimize else max
    winner = selector(candidates, key=lambda p: float(_metric_value(p, metric)))
    return {
        "temperature_c": winner["temperature_c"],
        "value": _metric_value(winner, metric),
        "status": "PREDICTED",
    }


def _best_helper_point(
    sweeps: Mapping[str, list[dict[str, Any]]], key: str
) -> dict[str, Any] | None:
    candidates: list[tuple[dict[str, Any], float, str]] = []
    for view, points in sweeps.items():
        for point in points:
            found = _helper_metric(point, key)
            if found is not None:
                candidates.append((point, found[0], view))
    if not candidates:
        return None
    point, value, view = min(candidates, key=lambda item: item[1])
    return {
        "temperature_c": point["temperature_c"],
        "value": value,
        "view": view,
        "status": "PREDICTED",
    }


def _helper_metric(point: Mapping[str, Any], key: str) -> tuple[float, str] | None:
    native = point.get("electrometry")
    if isinstance(native, Mapping) and key in native:
        try:
            number = float(native[key])
        except (TypeError, ValueError):
            pass
        else:
            if math.isfinite(number):
                return number, validate_status(native.get("status", "PREDICTED"))
    for result in point.get("optional_helpers", {}).values():
        if not isinstance(result, Mapping) or key not in result:
            continue
        value = result[key]
        if isinstance(value, Mapping):
            status = validate_status(value.get("status", result.get("status", "PREDICTED")))
            value = value.get("value")
        else:
            status = validate_status(result.get("status", "PREDICTED"))
        try:
            number = float(value)
        except (TypeError, ValueError):
            continue
        if math.isfinite(number):
            return number, status
    return None


def _make_figures(
    results: Mapping[str, Any], output_dir: Path, formats: list[str], dpi: int
) -> list[dict[str, Any]]:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    allowed_formats = {"png", "pdf", "svg"}
    unknown = set(formats) - allowed_formats
    if unknown:
        raise ValueError(f"unsupported figure formats: {sorted(unknown)}")
    artifacts: list[dict[str, Any]] = []

    def save(fig: Any, stem: str, description: str) -> None:
        for suffix in formats:
            path = output_dir / f"{stem}.{suffix}"
            fig.savefig(path, dpi=dpi if suffix == "png" else None, bbox_inches="tight")
            artifacts.append({
                "path": _serialized_path(path),
                "format": suffix, "description": description,
                "status": "PREDICTED",
            })
        plt.close(fig)

    sweeps = results["model"]["sweeps"]
    for view, points in sweeps.items():
        if not points or "spectrum" not in points[0]:
            continue
        fig, ax = plt.subplots(figsize=(8.2, 4.8))
        for point in points:
            spectrum = point["spectrum"]
            ax.plot(
                spectrum["detuning_mhz"], spectrum["transmission"], lw=1.3,
                label=f"{point['temperature_c']:g} C",
            )
        if view == "EIT":
            for dataset in results["inputs"]:
                if dataset["kind"] != "eit_spectrum" or dataset["state"] != "loaded":
                    continue
                rows = dataset["rows"]
                detuning = np.asarray([row["detuning_mhz"] for row in rows], dtype=float)
                signal = np.asarray([row["signal"] for row in rows], dtype=float)
                span = float(np.max(signal) - np.min(signal))
                normalized = (
                    (signal - np.min(signal)) / span if span > 0.0 else np.zeros_like(signal))
                ax.plot(
                    detuning, normalized, "--", color="black", lw=1.0, alpha=0.7,
                    label=(f"{dataset['id']} normalized "
                           f"[{dataset['status']}]")
                )
        ax.set_xlabel("Probe detuning [MHz]")
        ax.set_ylabel("Transmission")
        ax.set_title(f"GABES {view}: temperature sweep [PREDICTED]")
        ax.grid(alpha=0.25)
        ax.legend(ncol=2, fontsize=8)
        safe_view = "eit" if view == "EIT" else "at"
        save(fig, f"{safe_view}_temperature_spectra", f"{view} spectrum sweep")

    preferred_view = "EIT" if "EIT" in sweeps else next(iter(sweeps), None)
    if preferred_view and sweeps[preferred_view]:
        points = sweeps[preferred_view]
        temperatures = [p["temperature_c"] for p in points]
        panels = (
            ("eit_peak_contrast", "EIT peak contrast [1]"),
            ("eit_linewidth_mhz", "EIT linewidth [MHz]"),
            ("max_spectral_slope_per_mhz", "Max |dT/dnu| [1/MHz]"),
            ("if_discriminator_per_mhz", "IF discriminator [1/MHz]"),
        )
        fig, axes = plt.subplots(2, 2, figsize=(9.0, 6.6), sharex=True)
        for ax, (metric, label) in zip(axes.flat, panels):
            values = [_metric_value(p, metric) for p in points]
            ax.plot(temperatures, values, "o-", color="#0f766e", lw=1.6)
            ax.set_ylabel(label)
            ax.grid(alpha=0.3)
        for ax in axes[-1, :]:
            ax.set_xlabel("Cell temperature parameter [C]")
        fig.suptitle(f"Static temperature metrics ({preferred_view}) [PREDICTED]")
        fig.tight_layout()
        save(fig, "temperature_metrics", "Static temperature-dependent metrics")

    sensitivity_sets = [
        d for d in results["inputs"]
        if d["kind"] == "sensitivity" and d["state"] == "loaded"
    ]
    helper_series: list[tuple[float, float, str]] = []
    for points in sweeps.values():
        for point in points:
            found = _helper_metric(point, "field_sensitivity_nv_cm_sqrt_hz")
            if found:
                helper_series.append((point["temperature_c"], found[0], found[1]))
        if helper_series:
            break
    if sensitivity_sets or helper_series:
        fig, ax = plt.subplots(figsize=(7.5, 4.8))
        for dataset in sensitivity_sets:
            rows = dataset["rows"]
            x = [row["temperature_c"] for row in rows]
            y = [row["sensitivity_nv_cm_sqrt_hz"] for row in rows]
            yerr = [row.get("uncertainty_nv_cm_sqrt_hz") for row in rows]
            if any(value is not None for value in yerr):
                errors = [0.0 if value is None else value for value in yerr]
                ax.errorbar(x, y, yerr=errors, fmt="o", capsize=3,
                            label=f"{dataset['id']} [{dataset['status']}]")
            else:
                ax.plot(x, y, "o", label=f"{dataset['id']} [{dataset['status']}]")
        if helper_series:
            helper_series.sort()
            ax.plot([p[0] for p in helper_series], [p[1] for p in helper_series],
                    "-s", label=f"finite-IF model [{helper_series[0][2]}]")
        ax.set_xlabel("Temperature [C]")
        ax.set_ylabel("RF sensitivity [nV/cm/sqrt(Hz)]")
        ax.set_title("Sensitivity comparison (status-tagged)")
        ax.grid(alpha=0.3)
        ax.legend(fontsize=8)
        save(fig, "sensitivity_comparison", "Status-tagged sensitivity comparison")

    rf_sets = [
        d for d in results["inputs"]
        if d["kind"] == "rf_sweep" and d["state"] == "loaded"
    ]
    if rf_sets:
        fig, ax = plt.subplots(figsize=(7.5, 4.8))
        for dataset in rf_sets:
            ax.plot(
                [row["drive"] for row in dataset["rows"]],
                [row["response"] for row in dataset["rows"]],
                "o-", ms=3,
                label=f"{dataset['id']} [{dataset['status']}]",
            )
        ax.set_xlabel(f"Drive [{rf_sets[0].get('drive_unit', 'arb.')}]")
        ax.set_ylabel(f"Response [{rf_sets[0].get('response_unit', 'arb.')}]")
        ax.set_title("Imported RF sweeps")
        ax.grid(alpha=0.3)
        ax.legend(fontsize=8)
        save(fig, "rf_sweeps", "Imported RF response sweeps")

    psd_sets = [
        d for d in results["inputs"]
        if d["kind"] == "psd" and d["state"] == "loaded"
    ]
    if psd_sets:
        fig, ax = plt.subplots(figsize=(7.5, 4.8))
        for dataset in psd_sets:
            ax.loglog(
                [row["frequency_hz"] for row in dataset["rows"]],
                [row["asd"] for row in dataset["rows"]],
                label=f"{dataset['id']} [{dataset['status']}]",
            )
        ax.set_xlabel("Frequency [Hz]")
        ax.set_ylabel(f"ASD [{psd_sets[0].get('base_unit', 'arb.')}/sqrt(Hz)]")
        ax.set_title("Imported noise spectra")
        ax.grid(alpha=0.3, which="both")
        ax.legend(fontsize=8)
        save(fig, "noise_spectra", "Imported PSD/ASD traces")
    return artifacts


def _tex_escape(text: Any) -> str:
    value = str(text)
    replacements = {
        "\\": r"\textbackslash{}", "&": r"\&", "%": r"\%", "$": r"\$",
        "#": r"\#", "_": r"\_", "{": r"\{", "}": r"\}",
        "~": r"\textasciitilde{}", "^": r"\textasciicircum{}",
    }
    return "".join(replacements.get(char, char) for char in value)


def _write_tex_macros(results: Mapping[str, Any], path: Path) -> None:
    optima = results["model"]["optima"]
    slope = optima.get("max_spectral_slope") or {}
    if_opt = optima.get("max_if_discriminator") or {}
    loaded = [d for d in results["inputs"] if d["state"] == "loaded"]
    abs_status = results["capabilities"]["absolute_rf_sensitivity"]["status"]
    sensitivity_points: list[tuple[dict[str, Any], float]] = []
    for points in results["model"]["sweeps"].values():
        for point in points:
            found = _helper_metric(point, "field_sensitivity_nv_cm_sqrt_hz")
            if found:
                sensitivity_points.append((point, found[0]))
    best_sensitivity = (
        min(sensitivity_points, key=lambda item: item[1])
        if sensitivity_points else None)
    best_psn = (
        _helper_metric(best_sensitivity[0], "psn_limit_nv_cm_sqrt_hz")
        if best_sensitivity else None)

    sensitivity_rows = [
        row for dataset in loaded if dataset["kind"] == "sensitivity"
        for row in dataset.get("rows", [])
        if row.get("sensitivity_nv_cm_sqrt_hz") is not None
    ]
    room_scale = None
    best_row_psn = None
    optimum_scale = None
    if sensitivity_rows and sensitivity_points:
        room_row = min(
            sensitivity_rows, key=lambda row: abs(row["temperature_c"] - 20.0))
        room_point = min(
            sensitivity_points,
            key=lambda item: abs(item[0]["temperature_c"] - 20.0))
        if room_point[1]:
            room_scale = room_row["sensitivity_nv_cm_sqrt_hz"] / room_point[1]
        best_row = min(
            sensitivity_rows,
            key=lambda row: row["sensitivity_nv_cm_sqrt_hz"])
        best_row_psn = best_row.get("psn_limit_nv_cm_sqrt_hz")
        if best_row_psn and best_sensitivity:
            optimum_scale = best_row_psn / best_sensitivity[1]

    def fmt(value: Any, digits: int = 2) -> str:
        return "--" if value is None else f"{float(value):.{digits}f}"

    macros = {
        "AnalysisIdentifier": results["analysis_id"],
        "AnalysisGeneratedUTC": results["generated_at_utc"],
        "ModelTemperatureCount": len(results["model"]["temperatures_c"]),
        "SlopeOptimumTemperatureC": fmt(slope.get("temperature_c")),
        "SlopeOptimumValue": fmt(slope.get("value"), 4),
        "IFOptimumTemperatureC": fmt(if_opt.get("temperature_c")),
        "IFOptimumValue": fmt(if_opt.get("value"), 4),
        "ExperimentalInputStatus": loaded[0]["status"] if loaded else "PENDING",
        "AbsoluteSensitivityStatus": abs_status,
        "SensitivityOptimumTemperatureC": (
            fmt(best_sensitivity[0]["temperature_c"])
            if best_sensitivity else "--"),
        "TotalSensitivityOptimum": (
            fmt(best_sensitivity[1], 3) if best_sensitivity else "--"),
        "PSNSensitivityAtOptimum": (
            fmt(best_psn[0], 3) if best_psn else "--"),
        "ModelToPPTScaleRoom": fmt(room_scale),
        "PPTBestPSNLimit": fmt(best_row_psn),
        "ModelToPPTScaleOptimum": fmt(optimum_scale, 1),
    }
    lines = [
        "% Generated by analysis/rydberg_cell_heating/workflow.py; do not edit.",
    ]
    for name, value in macros.items():
        lines.append(rf"\providecommand{{\{name}}}{{{_tex_escape(value)}}}")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _git_metadata() -> dict[str, Any]:
    def run(*args: str) -> str | None:
        try:
            proc = subprocess.run(
                ["git", *args], cwd=ROOT, capture_output=True, text=True,
                check=False, timeout=5,
            )
        except (OSError, subprocess.TimeoutExpired):
            return None
        return proc.stdout.strip() if proc.returncode == 0 else None

    return {
        "commit": run("rev-parse", "HEAD"),
        "branch": run("branch", "--show-current"),
        "dirty": bool(run("status", "--porcelain")),
    }


def _source_notes(config: Mapping[str, Any]) -> list[dict[str, Any]]:
    notes = []
    for index, item in enumerate(config.get("source_notes", [])):
        note = dict(item)
        note["id"] = str(note.get("id", f"source_note_{index + 1}"))
        note["status"] = validate_status(note.get("status", "PENDING"))
        notes.append(note)
    return notes


def _temperature_log_rows(datasets: Iterable[Mapping[str, Any]]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for dataset in datasets:
        if dataset.get("kind") == "temperature_log" and dataset.get("state") == "loaded":
            rows.extend(dict(row) for row in dataset.get("rows", []))
    return rows


def _thermal_state(
    heater_setpoint_c: float,
    datasets: Iterable[Mapping[str, Any]],
    thermal_config: Mapping[str, Any],
) -> dict[str, Any]:
    """Resolve one setpoint without conflating setpoint and vapor temperature."""
    tolerance = float(thermal_config.get("setpoint_match_tolerance_c", 0.25))
    matches = [
        row for row in _temperature_log_rows(datasets)
        if row.get("setpoint_c") is not None
        and abs(float(row["setpoint_c"]) - heater_setpoint_c) <= tolerance
    ]
    # A temperature log may contain a warm-up sequence; use the latest elapsed
    # time deterministically and preserve that row in the output.
    row = max(matches, key=lambda item: float(item.get("elapsed_s") or 0.0)) if matches else None
    sensors: list[float] = []
    if row:
        for key in ("sensor_left_c", "sensor_right_c"):
            if row.get(key) is not None:
                sensors.append(float(row[key]))
    offsets = thermal_config.get("sensor_offsets_c")
    if offsets is not None and len(offsets) != len(sensors):
        raise ValueError(
            "model.thermal.sensor_offsets_c must match available temperature sensors")
    state = resolve_cell_temperature(
        heater_setpoint_c,
        sensors,
        sensor_offsets_c=offsets,
        sensor_weights=thermal_config.get("sensor_weights"),
        effective_vapor_temp_c=(
            None if not row else row.get("effective_vapor_temp_c")),
        cold_spot_temp_c=None if not row else row.get("cold_spot_c"),
    )
    row_status = validate_status(row.get("status", "MEASURED")) if row else "PENDING"
    effective_status = (
        row_status if state.effective_source != "heater_setpoint_fallback" else "ASSUMED")
    cold_status = (
        row_status if state.cold_spot_source != "heater_setpoint_fallback" else "ASSUMED")
    return {
        "heater_setpoint_c": state.heater_setpoint_c,
        "sensor_temperatures_c": list(state.sensor_temperatures_c),
        "sensor_offsets_c": list(state.sensor_offsets_c),
        "corrected_sensor_temperatures_c": list(state.corrected_sensor_temperatures_c),
        "effective_vapor_temp_c": state.effective_vapor_temp_c,
        "cold_spot_temp_c": state.cold_spot_temp_c,
        "effective_source": state.effective_source,
        "cold_spot_source": state.cold_spot_source,
        "effective_status": effective_status,
        "cold_spot_status": cold_status,
        "temperature_log_status": row_status,
        "selected_log_row": row,
        "saturated_cold_spot_density_m3": state.saturated_cold_spot_density_m3(),
    }


def _axial_profile_for_point(
    thermal: Mapping[str, Any],
    params: Mapping[str, Any],
    config: Mapping[str, Any],
) -> tuple[AxialCellProfile | None, dict[str, Any]]:
    """Build an optional axial T(z)/n(z) model with explicit provenance."""
    if not bool(config.get("enabled", False)):
        return None, {
            "available": False,
            "status": "PENDING",
            "message": "axial temperature profile is disabled",
        }
    input_status = validate_status(config.get("status", "ASSUMED"))
    length_m = float(params["cell_mm"]) * 1.0e-3
    density_mode = str(config.get("density_mode", "cold_spot_limited"))
    points = int(config.get("points", 21))
    if points < 2:
        raise ValueError("model.axial_profile.points must be at least 2")

    fractions = config.get("z_fraction")
    offsets = config.get("temperature_offsets_c")
    if fractions is not None or offsets is not None:
        if fractions is None or offsets is None:
            raise ValueError(
                "model.axial_profile.z_fraction and temperature_offsets_c "
                "must be supplied together")
        z_fraction = np.asarray(fractions, dtype=float)
        temperature_offsets = np.asarray(offsets, dtype=float)
        if (
            z_fraction.ndim != 1
            or z_fraction.shape != temperature_offsets.shape
            or z_fraction.size < 2
            or not np.isclose(z_fraction[0], 0.0)
            or not np.isclose(z_fraction[-1], 1.0)
        ):
            raise ValueError(
                "axial z_fraction/temperature_offsets_c must be equal-length "
                "1-D arrays spanning exactly 0 to 1")
        z_values = z_fraction * length_m
        temperature_values = (
            float(thermal["effective_vapor_temp_c"]) + temperature_offsets)
        pressure_cold_c = min(
            float(thermal["cold_spot_temp_c"]),
            float(np.min(temperature_values)),
        )
        profile = axial_cell_profile(
            z_values,
            temperature_values,
            density_mode=density_mode,
            pressure_cold_spot_temp_c=pressure_cold_c,
        )
        endpoint_source = "configured piecewise offsets"
    else:
        corrected = list(thermal.get("corrected_sensor_temperatures_c", []))
        use_sensors = (
            str(config.get("endpoint_source", "sensors_or_offsets"))
            == "sensors_or_offsets"
            and len(corrected) >= 2
        )
        if use_sensors:
            left_c, right_c = float(corrected[0]), float(corrected[-1])
            endpoint_source = "corrected end sensors"
        else:
            effective_c = float(thermal["effective_vapor_temp_c"])
            left_c = effective_c + float(config.get("left_offset_c", 0.0))
            right_c = effective_c + float(config.get("right_offset_c", 0.0))
            endpoint_source = "configured offsets from effective vapor temperature"
        if bool(config.get("anchor_cold_spot_at_left", False)):
            left_c = min(left_c, float(thermal["cold_spot_temp_c"]))
            endpoint_source += "; left endpoint includes declared cold spot"
        pressure_cold_c = min(
            float(thermal["cold_spot_temp_c"]), left_c, right_c)
        profile = linear_axial_cell_profile(
            length_m,
            left_c,
            right_c,
            points=points,
            density_mode=density_mode,
            pressure_cold_spot_temp_c=pressure_cold_c,
        )

    result = {
        "available": True,
        "status": "PREDICTED",
        "input_status": input_status,
        "density_mode": profile.density_mode,
        "endpoint_source": endpoint_source,
        "z_m": profile.z_m.tolist(),
        "temperature_c": profile.temperature_c.tolist(),
        "density_m3": profile.density_m3.tolist(),
        "declared_cold_spot_temp_c": float(thermal["cold_spot_temp_c"]),
        "cold_spot_temp_c": profile.cold_spot_temp_c,
        "profile_min_temp_c": profile.profile_min_temp_c,
        "cold_spot_density_m3": profile.cold_spot_density_m3,
        "column_density_m2": profile.column_density_m2,
    }
    return profile, result


def _attach_segmented_beer_lambert(
    point: dict[str, Any],
    raw: Mapping[str, Any],
    profile: AxialCellProfile | None,
) -> None:
    """Attach an axial Beer--Lambert spectrum using a declared approximation.

    The current four-level solve is lumped.  This workflow preserves its local
    line shape and rescales alpha(z) by the pressure-equilibrated density
    profile before segment integration.  The approximation is labelled rather
    than presented as a full spatial OBE solve.
    """
    if profile is None:
        return
    axial = point["axial_profile"]
    spectrum = point.get("spectrum")
    if not isinstance(spectrum, Mapping):
        axial["beer_lambert"] = {
            "available": False,
            "status": "PENDING",
            "message": "include_spectra_in_results is required",
        }
        return
    length_m = float(point["parameters"]["cell_mm"]) * 1.0e-3
    transmission = np.asarray(spectrum["transmission"], dtype=float)
    local_alpha = -np.log(np.clip(transmission, np.finfo(float).tiny, None)) / length_m
    alpha_z = (
        profile.density_m3[:, None]
        / max(float(raw["N"]), np.finfo(float).tiny)
        * local_alpha[None, :]
    )
    integrated = integrate_beer_lambert(profile.z_m, alpha_z)
    axial["beer_lambert"] = {
        "available": True,
        "status": "PREDICTED",
        "method": "segmented density-rescaled Beer-Lambert integration",
        "qualifier": (
            "uses the lumped OBE line shape at the effective vapor temperature; "
            "local density varies with T(z), but local line-shape/dephasing "
            "changes are not re-solved"
        ),
        "used_in_finite_if_sensitivity": False,
        "sensitivity_scope": (
            "finite-IF sensitivity remains the lumped effective-temperature model; "
            "a spatial complex-phasor integration is intentionally not implied"
        ),
        "detuning_mhz": list(spectrum["detuning_mhz"]),
        "transmission": np.asarray(integrated.transmission).tolist(),
        "optical_depth": np.asarray(integrated.optical_depth).tolist(),
    }


def _sam_result(config: Mapping[str, Any]) -> dict[str, Any]:
    """Evaluate an optional standard-antenna-method field calibration."""
    if not bool(config.get("enabled", False)):
        return {
            "available": False,
            "status": "PENDING",
            "message": str(config.get("note", "SAM calibration is disabled")),
        }
    required = ("source_power_dbm", "antenna_gain_dbi", "distance_m")
    missing = [name for name in required if config.get(name) is None]
    if missing:
        return {
            "available": False,
            "status": "PENDING",
            "message": f"missing SAM fields: {', '.join(missing)}",
        }
    optional_names = (
        "cable_loss_db", "additional_loss_db", "field_correction",
        "source_power_std_db", "antenna_gain_std_db", "cable_loss_std_db",
        "additional_loss_std_db", "distance_std_m", "field_correction_std",
        "amplitude_convention", "frequency_hz", "antenna_max_dimension_m",
    )
    kwargs = {
        name: config[name]
        for name in optional_names
        if config.get(name) is not None
    }
    calibration = sam_field_calibration(
        float(config["source_power_dbm"]),
        float(config["antenna_gain_dbi"]),
        float(config["distance_m"]),
        **kwargs,
    )
    return {
        "available": True,
        "status": "PREDICTED",
        "input_status": validate_status(config.get("status", "ASSUMED")),
        "field_v_m": calibration.field_v_m,
        "field_nv_cm": calibration.field_nv_cm,
        "standard_uncertainty_v_m": calibration.standard_uncertainty_v_m,
        "standard_uncertainty_nv_cm": calibration.standard_uncertainty_nv_cm,
        "relative_standard_uncertainty": calibration.relative_standard_uncertainty,
        "power_at_antenna_w": calibration.power_at_antenna_w,
        "net_gain_dbi": calibration.net_gain_dbi,
        "amplitude_convention": calibration.amplitude_convention,
        "far_field_ratio": calibration.far_field_ratio,
        "warning": calibration.warning,
    }


def _sam_at_comparison(
    point: Mapping[str, Any],
    sam: Mapping[str, Any],
    sam_config: Mapping[str, Any],
) -> dict[str, Any]:
    """Compare SAM field with the model/fitted-Rabi field without hiding scope."""
    if not sam.get("available", False):
        return {
            "available": False,
            "status": "PENDING",
            "message": "SAM calibration is unavailable",
        }
    at_field = _metric_value(point, "at_inferred_field_v_m")
    if at_field is None:
        return {
            "available": False,
            "status": "PENDING",
            "message": "AT/full-spectrum RF field is unavailable for this point",
        }
    point_convention = str(
        point["parameters"].get("rf_field_convention", "RMS")).lower()
    sam_convention = str(sam["amplitude_convention"]).lower()
    convention_factor = 1.0
    if sam_convention != point_convention:
        convention_factor = (
            np.sqrt(2.0)
            if sam_convention == "rms" and point_convention == "peak"
            else 1.0 / np.sqrt(2.0)
        )
    sam_field = float(sam["field_v_m"]) * convention_factor
    sam_uncertainty = float(sam["standard_uncertainty_v_m"]) * convention_factor

    sam_frequency = sam_config.get("frequency_hz")
    model_frequency = float(
        point["parameters"].get("mw_frequency_ghz", 37.0)) * 1.0e9
    tolerance = float(sam_config.get("frequency_tolerance_hz", 1.0e6))
    frequency_match = (
        None if sam_frequency is None
        else abs(float(sam_frequency) - model_frequency) <= tolerance
    )
    ratio = at_field / sam_field
    residual = at_field - sam_field
    return {
        "available": frequency_match is not False,
        "status": "PREDICTED" if frequency_match is not False else "PENDING",
        "field_convention": point_convention,
        "sam_field_converted_v_m": sam_field,
        "sam_standard_uncertainty_converted_v_m": sam_uncertainty,
        "model_rf_field_v_m": at_field,
        "model_to_sam_ratio": ratio,
        "model_minus_sam_v_m": residual,
        "ratio_standard_uncertainty_from_sam": (
            ratio * float(sam["relative_standard_uncertainty"])),
        "normalized_residual_using_sam_u": (
            None if sam_uncertainty <= 0.0 else residual / sam_uncertainty),
        "model_frequency_hz": model_frequency,
        "sam_frequency_hz": sam_frequency,
        "frequency_tolerance_hz": tolerance,
        "frequency_match": frequency_match,
        "qualifier": (
            "model_rf_field_v_m is converted from the configured/fitted Ω_RF "
            "parameter and is PREDICTED here; it is not a measured AT fit"
        ),
        "warning": (
            "SAM and model frequencies do not match within tolerance"
            if frequency_match is False else sam.get("warning")
        ),
    }


def _attach_effective_atom_number(
    point: dict[str, Any],
    raw: Mapping[str, Any],
    config: Mapping[str, Any],
    axial_profile: AxialCellProfile | None = None,
) -> None:
    if not bool(config.get("enabled", False)):
        point["effective_atom_number"] = {
            "available": False, "status": "PENDING",
            "message": "effective-atom geometry disabled",
        }
        return
    params = point["parameters"]
    probe_radius_m = float(params["beam_diameter_mm"]) * 0.5e-3
    coupling_diameter = config.get("coupling_beam_diameter_mm")
    coupling_radius_m = (
        None if coupling_diameter is None else float(coupling_diameter) * 0.5e-3)
    density_input: float | AxialCellProfile = (
        float(raw["N"]) if axial_profile is None else axial_profile)
    estimate_kwargs = {
        "probe_radius_m": probe_radius_m,
        "coupling_radius_m": coupling_radius_m,
        "participation_fraction": float(config.get("participation_fraction", 1.0)),
        "overlap_efficiency": float(config.get("overlap_efficiency", 1.0)),
    }
    if axial_profile is None:
        estimate_kwargs["length_m"] = float(params["cell_mm"]) * 1.0e-3
    estimate = effective_atom_number(density_input, **estimate_kwargs)
    input_status = validate_status(config.get("status", "ASSUMED"))
    point["effective_atom_number"] = {
        "available": True,
        "status": "PREDICTED",
        "input_status": input_status,
        "atoms": estimate.atoms,
        "column_density_m2": estimate.column_density_m2,
        "effective_area_m2": estimate.effective_area_m2,
        "participation_fraction": estimate.participation_fraction,
        "overlap_efficiency": estimate.overlap_efficiency,
        "probe_radius_m": probe_radius_m,
        "coupling_radius_m": coupling_radius_m,
        "qualifier": "geometric effective atom estimate, not a fitted atom count",
        "axial_density_integrated": axial_profile is not None,
    }
    point["metrics"]["effective_atom_number"] = {
        "value": estimate.atoms,
        "unit": "atoms",
        "status": "PREDICTED",
        "input_status": input_status,
    }


def run_analysis(
    config_path: str | os.PathLike[str],
    *,
    output_dir: str | os.PathLike[str] | None = None,
    make_figures: bool = True,
) -> dict[str, Path]:
    """Run the configured workflow and return paths to emitted artifacts."""
    config_file = Path(config_path).expanduser().resolve()
    # Capture repository state before creating the output directory or writing
    # generated artifacts, so this run cannot make its own provenance "dirty".
    git_metadata = _git_metadata()
    config = json.loads(config_file.read_text(encoding="utf-8-sig"))
    if int(config.get("schema_version", -1)) != SCHEMA_VERSION:
        raise ValueError(
            f"config schema_version must be {SCHEMA_VERSION}, got "
            f"{config.get('schema_version')!r}")
    config_dir = config_file.parent
    configured_output = output_dir or config.get("output_dir", "generated")
    out = _resolve(config_dir, str(configured_output))
    out.mkdir(parents=True, exist_ok=True)

    datasets, input_provenance = load_inputs(config.get("inputs", []), config_dir)
    model_config = dict(config.get("model", {}))
    temperatures = [float(value) for value in model_config.get("temperatures_c", [])]
    if not temperatures or any(not math.isfinite(value) for value in temperatures):
        raise ValueError("model.temperatures_c must contain finite numeric values")
    if len(set(temperatures)) != len(temperatures):
        raise ValueError("model.temperatures_c must not contain duplicates")
    views = [str(value) for value in model_config.get("views", ["EIT"])]
    if not views:
        raise ValueError("model.views must not be empty")
    overrides = dict(model_config.get("base_params", {}))
    include_spectrum = bool(model_config.get("include_spectra_in_results", True))

    static = StaticRydbergAdapter()
    electrometry = ElectrometryAdapter(static, config.get("electrometry", {}))
    electrometry_summary = electrometry.summary()
    thermal_config = dict(model_config.get("thermal", {}))
    atom_config = dict(model_config.get("effective_atom_number", {}))
    axial_config = dict(model_config.get("axial_profile", {}))
    sam_config = dict(config.get("sam_calibration", {}))
    sam = _sam_result(sam_config)
    helper_adapters = [OptionalHelperAdapter(spec) for spec in config.get("helpers", [])]
    helper_summaries = [helper.discover() for helper in helper_adapters]
    sweeps: dict[str, list[dict[str, Any]]] = {}
    for view in views:
        points: list[dict[str, Any]] = []
        for temperature in temperatures:
            thermal = _thermal_state(temperature, datasets, thermal_config)
            axial_params = static.parameters(view, overrides)
            axial_profile, axial_result = _axial_profile_for_point(
                thermal, axial_params, axial_config)
            solver_thermal = dict(thermal)
            solver_thermal["declared_cold_spot_temp_c"] = thermal["cold_spot_temp_c"]
            if axial_profile is not None:
                solver_thermal["cold_spot_temp_c"] = axial_profile.cold_spot_temp_c
                solver_thermal["cold_spot_source"] = (
                    "minimum of declared cold spot and axial T(z)")
            point = static.temperature_point(
                solver_thermal["effective_vapor_temp_c"], view, overrides,
                include_spectrum=include_spectrum,
                heater_setpoint_c=solver_thermal["heater_setpoint_c"],
                cold_spot_temp_c=solver_thermal["cold_spot_temp_c"],
            )
            # The sweep axis is a heater setpoint.  Preserve the solver's vapor
            # temperature separately instead of silently relabelling it.
            point.result["temperature_c"] = temperature
            point.result["heater_setpoint_c"] = temperature
            point.result["effective_vapor_temp_c"] = solver_thermal[
                "effective_vapor_temp_c"]
            point.result["thermal_state"] = solver_thermal
            point.result["axial_profile"] = axial_result
            _attach_segmented_beer_lambert(
                point.result, point.raw, axial_profile)
            _attach_effective_atom_number(
                point.result, point.raw, atom_config, axial_profile)
            point.result["electrometry"] = electrometry.evaluate(point)
            point.result["sam_comparison"] = _sam_at_comparison(
                point.result, sam, sam_config)
            if helper_adapters:
                point.result["optional_helpers"] = {
                    helper.module_name: helper.evaluate(point) for helper in helper_adapters
                }
            points.append(point.result)
        sweeps[view] = points

    preferred = sweeps.get("EIT") or next(iter(sweeps.values()))
    optima = {
        "max_spectral_slope": _best_point(preferred, "max_spectral_slope_per_mhz"),
        "max_if_discriminator": _best_point(preferred, "if_discriminator_per_mhz"),
        "min_total_field_sensitivity": _best_helper_point(
            sweeps, "field_sensitivity_nv_cm_sqrt_hz"),
    }
    has_absolute = any(
        _helper_metric(point, "field_sensitivity_nv_cm_sqrt_hz") is not None
        for points in sweeps.values() for point in points
    )
    has_superhet = any(
        _helper_metric(point, "superheterodyne_responsivity") is not None
        for points in sweeps.values() for point in points
    )
    has_axial_beer_lambert = any(
        bool(point.get("axial_profile", {}).get(
            "beer_lambert", {}).get("available", False))
        for points in sweeps.values() for point in points
    )
    has_sam_at_comparison = any(
        bool(point.get("sam_comparison", {}).get("available", False))
        for points in sweeps.values() for point in points
    )
    experimental_loaded = any(dataset["state"] == "loaded" for dataset in datasets)
    results: dict[str, Any] = {
        "schema_version": SCHEMA_VERSION,
        "analysis_id": str(config.get("analysis_id", config_file.stem)),
        "generated_at_utc": _utc_now(),
        "evidence_status_vocabulary": sorted(ALLOWED_STATUSES),
        "source_notes": _source_notes(config),
        "inputs": datasets,
        "model": {
            "scheme": static.identity,
            "temperatures_c": temperatures,
            "views": views,
            "base_params": overrides,
            "sweeps": sweeps,
            "optima": optima,
        },
        "optional_helpers": helper_summaries,
        "electrometry_adapter": electrometry_summary,
        "sam_calibration": sam,
        "capabilities": {
            "static_rydberg_spectrum": {
                "status": "PREDICTED", "available": True,
            },
            "experimental_csv_overlay": {
                "status": "MEASURED" if experimental_loaded else "PENDING",
                "available": experimental_loaded,
            },
            "finite_if_superheterodyne": {
                "status": "PREDICTED" if has_superhet else "PENDING",
                "available": has_superhet,
                "note": "first-order frequency-domain weak-SIG response",
            },
            "two_tone_superheterodyne": {
                "status": "PREDICTED" if has_superhet else "PENDING",
                "available": has_superhet,
                "note": (
                    "backward-compatible key: finite-IF linear response only; "
                    "full time-domain LO+SIG/lock-in remains deferred"
                ),
            },
            "absolute_rf_sensitivity": {
                "status": "PREDICTED" if has_absolute else "PENDING",
                "available": has_absolute,
                "note": (
                    "Available only through an explicit RF dipole and detector "
                    "noise/throughput calibration; no literature sensitivity is "
                    "injected as a target."
                ),
            },
            "axial_beer_lambert": {
                "status": "PREDICTED" if has_axial_beer_lambert else "PENDING",
                "available": has_axial_beer_lambert,
                "note": (
                    "diagnostic axial transmission/column-density path; finite-IF "
                    "sensitivity remains the separately labelled lumped model"
                ),
            },
            "sam_field_calibration": {
                "status": sam["status"],
                "available": sam["available"],
            },
            "sam_at_comparison": {
                "status": "PREDICTED" if has_sam_at_comparison else "PENDING",
                "available": has_sam_at_comparison,
            },
        },
    }

    figures_config = dict(config.get("figures", {}))
    formats = [str(value).lower() for value in figures_config.get("formats", ["png", "pdf"])]
    dpi = int(figures_config.get("dpi", 180))
    figures = _make_figures(results, out, formats, dpi) if make_figures else []
    results["figures"] = figures

    results_path = out / "results.json"
    macros_path = out / "results_macros.tex"
    manifest_path = out / "source_manifest.json"
    _json_dump(results_path, results)
    _write_tex_macros(results, macros_path)

    code_paths = [Path(__file__).resolve(), Path(__file__).with_name("adapters.py").resolve()]
    scheme_path = ROOT / "gabes" / "schemes" / "rydberg.py"
    if scheme_path.is_file():
        code_paths.append(scheme_path)
    for name in (
        "rydberg_electrometry.py", "rydberg_experiment.py",
        "rydberg_experimental_csv.py",
    ):
        path = ROOT / "gabes" / name
        if path.is_file():
            code_paths.append(path)
    manifest = {
        "schema_version": SCHEMA_VERSION,
        "analysis_id": results["analysis_id"],
        "generated_at_utc": results["generated_at_utc"],
        "config": {
            "path": _serialized_path(config_file),
            "sha256": _sha256(config_file),
            "bytes": config_file.stat().st_size,
        },
        "inputs": input_provenance,
        "code": [
            {
                "path": _serialized_path(path),
                "sha256": _sha256(path),
                "bytes": path.stat().st_size,
            }
            for path in code_paths
        ],
        "runtime": {
            "python": platform.python_version(),
            "platform": platform.platform(),
            "numpy": np.__version__,
        },
        "git": git_metadata,
        "artifacts": [
            {
                "path": _serialized_path(results_path),
                "role": "machine-readable results",
            },
            {
                "path": _serialized_path(macros_path),
                "role": "generated TeX macros",
            },
            *figures,
        ],
        "scope_note": (
            "The workflow neither embeds nor modifies the source PPTX. A PPT-only "
            "transcription must be supplied as a status=PPT CSV or source note."
        ),
    }
    _json_dump(manifest_path, manifest)
    return {
        "output_dir": out,
        "results": results_path,
        "manifest": manifest_path,
        "tex_macros": macros_path,
    }


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Generate status-tagged Rydberg cell-heating analysis artifacts.")
    parser.add_argument("--config", required=True, help="Path to JSON config")
    parser.add_argument(
        "--output-dir", help="Override output directory (relative to current directory)")
    parser.add_argument(
        "--no-figures", action="store_true", help="Emit JSON/TeX only")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    output_override = None
    if args.output_dir:
        output_override = str(Path(args.output_dir).expanduser().resolve())
    artifacts = run_analysis(
        args.config, output_dir=output_override, make_figures=not args.no_figures)
    for name, path in artifacts.items():
        print(f"{name}: {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
