"""Robust CSV import for Rydberg EIT, RF sweeps, and noise spectra.

The three public loaders standardize only the independent-axis units needed by
analysis.  Detector signals remain in their declared units; no hidden baseline,
normalization, or denoising is applied.  Every returned trace carries a SHA-256
digest and parsing diagnostics so a generated report can distinguish raw data,
unit assumptions, and later fitted quantities.
"""

from __future__ import annotations

import csv
import hashlib
import io
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

import numpy as np


MAX_RYDBERG_CSV_BYTES = 25 * 1024 * 1024
MAX_RYDBERG_CSV_ROWS = 1_000_000
MIN_TRACE_POINTS = 3


class RydbergCSVError(ValueError):
    """Raised when a Rydberg experimental trace cannot be imported safely."""


@dataclass(frozen=True)
class TraceProvenance:
    """Source identity, parser choices, and row-level diagnostics."""

    source_name: str
    sha256: str
    file_size_bytes: int
    encoding: str
    delimiter: str
    header: tuple[str, ...] | None
    x_column: str
    y_column: str
    total_rows: int
    valid_rows: int
    ignored_rows: int
    nonfinite_rows: int
    duplicate_rows_merged: int
    assumptions: tuple[str, ...]


@dataclass(frozen=True)
class EITTrace:
    """Probe detuning in MHz and an unmodified detector signal."""

    detuning_mhz: np.ndarray
    signal: np.ndarray
    signal_unit: str
    provenance: TraceProvenance


@dataclass(frozen=True)
class RFSweepTrace:
    """RF drive and measured response.

    Field drives are standardized to V/m, Rabi-frequency drives to MHz, source
    powers to dBm or W, and generic drives are left in the requested unit.
    """

    drive: np.ndarray
    response: np.ndarray
    drive_quantity: str
    drive_unit: str
    response_unit: str
    provenance: TraceProvenance

    @property
    def field_v_m(self) -> np.ndarray:
        if self.drive_quantity != "field":
            raise AttributeError("This RF sweep is not field-calibrated.")
        return self.drive


@dataclass(frozen=True)
class PSDTrace:
    """Frequency axis plus both ASD and PSD in standardized base units."""

    frequency_hz: np.ndarray
    asd: np.ndarray
    psd: np.ndarray
    quantity: str
    base_unit: str
    input_spectrum_kind: str
    provenance: TraceProvenance


@dataclass(frozen=True)
class _ParsedColumns:
    x: np.ndarray
    y: np.ndarray
    provenance: TraceProvenance
    x_label: str
    y_label: str


def _readonly(values: Sequence[float] | np.ndarray) -> np.ndarray:
    result = np.array(values, dtype=float, copy=True)
    result.setflags(write=False)
    return result


def _read_source(
    source: bytes | bytearray | memoryview | str | Path,
    source_name: str | None,
) -> tuple[bytes, str]:
    if isinstance(source, (bytes, bytearray, memoryview)):
        data = bytes(source)
        name = source_name or "uploaded.csv"
    elif isinstance(source, (str, Path)):
        path = Path(source)
        try:
            data = path.read_bytes()
        except OSError as exc:
            raise RydbergCSVError(f"Could not read CSV source: {exc}") from exc
        name = source_name or str(path)
    else:
        raise TypeError("source must be CSV bytes or a filesystem path.")
    if not data:
        raise RydbergCSVError("CSV source is empty.")
    if len(data) > MAX_RYDBERG_CSV_BYTES:
        raise RydbergCSVError(
            f"CSV exceeds the {MAX_RYDBERG_CSV_BYTES // (1024 * 1024)} MiB limit."
        )
    return data, str(name)


def _decode(data: bytes) -> tuple[str, str]:
    if data.startswith((b"\xff\xfe", b"\xfe\xff")):
        try:
            return data.decode("utf-16"), "utf-16"
        except UnicodeDecodeError as exc:
            raise RydbergCSVError("Malformed UTF-16 CSV text.") from exc
    for encoding in ("utf-8-sig", "utf-8", "cp949"):
        try:
            return data.decode(encoding), encoding
        except UnicodeDecodeError:
            continue
    return data.decode("latin-1"), "latin-1"


def _delimiter(text: str) -> str:
    lines = [line for line in text.splitlines() if line.strip()][:50]
    if not lines:
        raise RydbergCSVError("CSV contains no non-blank lines.")
    candidates = (",", "\t", ";")
    scores = {
        candidate: sum(line.count(candidate) for line in lines)
        for candidate in candidates
    }
    chosen = max(candidates, key=lambda candidate: scores[candidate])
    if scores[chosen] == 0:
        # Whitespace-only numeric files are common spectrum-analyzer exports.
        return "whitespace"
    return chosen


def _rows(text: str, delimiter: str) -> list[list[str]]:
    try:
        if delimiter == "whitespace":
            result = [re.split(r"\s+", line.strip()) for line in text.splitlines()]
        else:
            result = list(csv.reader(io.StringIO(text, newline=""), delimiter=delimiter))
    except csv.Error as exc:
        raise RydbergCSVError(f"Malformed CSV: {exc}") from exc
    if len(result) > MAX_RYDBERG_CSV_ROWS:
        raise RydbergCSVError(
            f"CSV exceeds the {MAX_RYDBERG_CSV_ROWS:,}-row limit."
        )
    return result


def _normal_header(value: str) -> str:
    return re.sub(r"[^a-z0-9]+", "", value.strip().lower())


def _header_matches(value: str, aliases: Sequence[str]) -> bool:
    normalized = _normal_header(value)
    return any(_normal_header(alias) in normalized for alias in aliases)


def _find_named_column(header: Sequence[str], selector: str) -> int:
    target = _normal_header(selector)
    exact = [index for index, value in enumerate(header) if _normal_header(value) == target]
    if len(exact) == 1:
        return exact[0]
    partial = [index for index, value in enumerate(header) if target in _normal_header(value)]
    if len(partial) == 1:
        return partial[0]
    raise RydbergCSVError(f"Column name {selector!r} was not found unambiguously.")


def _as_float(value: str) -> float:
    # Strip ordinary thousands separators only when they cannot be the active
    # CSV delimiter (a quoted "1,234" reaches this function intact).
    cleaned = value.strip().replace("\u2212", "-")
    if re.fullmatch(r"[+-]?\d{1,3}(,\d{3})+(\.\d+)?([eE][+-]?\d+)?", cleaned):
        cleaned = cleaned.replace(",", "")
    return float(cleaned)


def _parse_columns(
    source: bytes | bytearray | memoryview | str | Path,
    *,
    source_name: str | None,
    x_aliases: Sequence[str],
    y_aliases: Sequence[str],
    x_column: int | str | None,
    y_column: int | str | None,
) -> _ParsedColumns:
    data, name = _read_source(source, source_name)
    text, encoding = _decode(data)
    delimiter = _delimiter(text)
    rows = _rows(text, delimiter)

    header_index = None
    inferred_x = None
    inferred_y = None
    if isinstance(x_column, str) or isinstance(y_column, str):
        # An explicitly named column is itself enough to identify a header; the
        # other column may legitimately be selected by a numeric index.
        for index, row in enumerate(rows[:100]):
            if isinstance(x_column, str):
                target = _normal_header(x_column)
                x_candidates = [
                    position for position, value in enumerate(row)
                    if target in _normal_header(value)
                ]
            elif x_column is None:
                x_candidates = [
                    position for position, value in enumerate(row)
                    if _header_matches(value, x_aliases)
                ] or ([0] if row else [])
            else:
                x_candidates = [int(x_column)] if int(x_column) < len(row) else []
            if isinstance(y_column, str):
                target = _normal_header(y_column)
                y_candidates = [
                    position for position, value in enumerate(row)
                    if target in _normal_header(value)
                ]
            elif y_column is None:
                y_candidates = [
                    position for position, value in enumerate(row)
                    if _header_matches(value, y_aliases)
                ] or ([1] if len(row) > 1 else [])
            else:
                y_candidates = [int(y_column)] if int(y_column) < len(row) else []
            distinct = [
                (x, y) for x in x_candidates for y in y_candidates if x != y
            ]
            if distinct:
                header_index = index
                inferred_x, inferred_y = distinct[0]
                break
    else:
        for index, row in enumerate(rows[:100]):
            x_candidates = [
                position for position, value in enumerate(row)
                if _header_matches(value, x_aliases)
            ]
            y_candidates = [
                position for position, value in enumerate(row)
                if _header_matches(value, y_aliases)
            ]
            distinct = [
                (x, y) for x in x_candidates for y in y_candidates if x != y
            ]
            if distinct:
                header_index = index
                inferred_x, inferred_y = distinct[0]
                break

    header = tuple(rows[header_index]) if header_index is not None else None
    if isinstance(x_column, str) or isinstance(y_column, str):
        if header is None:
            raise RydbergCSVError("Named columns require a recognizable header row.")
    if isinstance(x_column, str):
        x_index = _find_named_column(header or (), x_column)
    elif x_column is None:
        x_index = inferred_x if inferred_x is not None else 0
    else:
        x_index = int(x_column)
    if isinstance(y_column, str):
        y_index = _find_named_column(header or (), y_column)
    elif y_column is None:
        y_index = inferred_y if inferred_y is not None else 1
    else:
        y_index = int(y_column)
    if x_index < 0 or y_index < 0 or x_index == y_index:
        raise RydbergCSVError("CSV column indices must be distinct and non-negative.")

    start = header_index + 1 if header_index is not None else 0
    x_values: list[float] = []
    y_values: list[float] = []
    # Rows before/including a recognized header are deliberately not data and
    # therefore count as ignored in the provenance accounting.
    ignored = start
    nonfinite = 0
    for row in rows[start:]:
        if not row or all(not value.strip() for value in row):
            ignored += 1
            continue
        if max(x_index, y_index) >= len(row):
            ignored += 1
            continue
        try:
            x_value = _as_float(row[x_index])
            y_value = _as_float(row[y_index])
        except (ValueError, OverflowError):
            ignored += 1
            continue
        if not np.isfinite(x_value) or not np.isfinite(y_value):
            nonfinite += 1
            continue
        x_values.append(x_value)
        y_values.append(y_value)

    x = np.asarray(x_values, dtype=float)
    y = np.asarray(y_values, dtype=float)
    if x.size == 0:
        raise RydbergCSVError("CSV contains no finite numeric data in selected columns.")
    order = np.argsort(x, kind="mergesort")
    x = x[order]
    y = y[order]
    unique_x, starts, counts = np.unique(x, return_index=True, return_counts=True)
    merged_y = y[starts].copy()
    for group in np.flatnonzero(counts > 1):
        begin = starts[group]
        merged_y[group] = np.median(y[begin : begin + counts[group]])
    duplicates = int(x.size - unique_x.size)
    if unique_x.size < MIN_TRACE_POINTS:
        raise RydbergCSVError(
            f"Trace requires at least {MIN_TRACE_POINTS} unique x values; "
            f"found {unique_x.size}."
        )

    x_label = (
        header[x_index]
        if header is not None and x_index < len(header)
        else f"column {x_index}"
    )
    y_label = (
        header[y_index]
        if header is not None and y_index < len(header)
        else f"column {y_index}"
    )
    assumptions = []
    if header is None:
        assumptions.append("No header recognized; selected columns by index.")
    provenance = TraceProvenance(
        source_name=name,
        sha256=hashlib.sha256(data).hexdigest(),
        file_size_bytes=len(data),
        encoding=encoding,
        delimiter=delimiter,
        header=header,
        x_column=x_label,
        y_column=y_label,
        total_rows=len(rows),
        valid_rows=int(x.size),
        ignored_rows=int(ignored),
        nonfinite_rows=int(nonfinite),
        duplicate_rows_merged=duplicates,
        assumptions=tuple(assumptions),
    )
    return _ParsedColumns(
        x=_readonly(unique_x),
        y=_readonly(merged_y),
        provenance=provenance,
        x_label=x_label,
        y_label=y_label,
    )


def _canonical_unit(unit: str) -> str:
    return (
        str(unit)
        .strip()
        .lower()
        .replace("\u03bc", "u")
        .replace("\u00b5", "u")
        .replace("\u221a", "sqrt")
        .replace("\u00b2", "2")
        .replace("**", "^")
        .replace(" ", "")
    )


def _unit_from_label(label: str, choices: Sequence[str]) -> str | None:
    normalized_label = _canonical_unit(label)
    # Longest first avoids finding V/m inside nV/m.
    for choice in sorted(choices, key=lambda item: len(_canonical_unit(item)), reverse=True):
        if _canonical_unit(choice) in normalized_label:
            return choice
    return None


def _with_assumption(
    provenance: TraceProvenance, assumption: str
) -> TraceProvenance:
    return TraceProvenance(
        **{
            **provenance.__dict__,
            "assumptions": provenance.assumptions + (assumption,),
        }
    )


_FREQUENCY_TO_HZ = {
    "hz": 1.0,
    "khz": 1.0e3,
    "mhz": 1.0e6,
    "ghz": 1.0e9,
}


def _frequency_factor(unit: str, target_unit: str) -> float:
    canonical = _canonical_unit(unit)
    if canonical not in _FREQUENCY_TO_HZ:
        raise RydbergCSVError(f"Unsupported frequency unit: {unit!r}.")
    return _FREQUENCY_TO_HZ[canonical] / _FREQUENCY_TO_HZ[target_unit]


_FIELD_TO_V_M = {
    "v/m": 1.0,
    "mv/m": 1.0e-3,
    "uv/m": 1.0e-6,
    "nv/m": 1.0e-9,
    "v/cm": 1.0e2,
    "mv/cm": 1.0e-1,
    "uv/cm": 1.0e-4,
    "nv/cm": 1.0e-7,
}


def load_eit_csv(
    source: bytes | bytearray | memoryview | str | Path,
    *,
    source_name: str | None = None,
    x_column: int | str | None = None,
    signal_column: int | str | None = None,
    detuning_unit: str | None = None,
    signal_unit: str = "arb.",
) -> EITTrace:
    """Load an EIT spectrum and standardize probe detuning to MHz."""

    parsed = _parse_columns(
        source,
        source_name=source_name,
        x_aliases=("detuning", "probe frequency", "frequency", "freq"),
        y_aliases=("transmission", "eit signal", "photodiode", "voltage", "signal"),
        x_column=x_column,
        y_column=signal_column,
    )
    unit = detuning_unit or _unit_from_label(
        parsed.x_label, ("GHz", "MHz", "kHz", "Hz")
    )
    provenance = parsed.provenance
    if unit is None:
        unit = "MHz"
        provenance = _with_assumption(
            provenance, "Detuning unit was not declared; assumed MHz."
        )
    factor = _frequency_factor(unit, "mhz")
    return EITTrace(
        detuning_mhz=_readonly(parsed.x * factor),
        signal=parsed.y,
        signal_unit=str(signal_unit),
        provenance=provenance,
    )


def load_rf_sweep_csv(
    source: bytes | bytearray | memoryview | str | Path,
    *,
    source_name: str | None = None,
    drive_column: int | str | None = None,
    response_column: int | str | None = None,
    drive_quantity: str = "field",
    drive_unit: str | None = None,
    response_unit: str = "arb.",
) -> RFSweepTrace:
    """Load a weak-SIG/AT RF sweep with an explicit drive quantity."""

    parsed = _parse_columns(
        source,
        source_name=source_name,
        x_aliases=("electric field", "field", "e sig", "sam", "rabi", "power", "drive"),
        y_aliases=("response", "beat", "amplitude", "voltage", "signal", "splitting"),
        x_column=drive_column,
        y_column=response_column,
    )
    quantity = str(drive_quantity).strip().lower().replace(" ", "_")
    provenance = parsed.provenance
    if quantity == "field":
        unit = drive_unit or _unit_from_label(
            parsed.x_label,
            ("nV/cm", "uV/cm", "mV/cm", "V/cm", "nV/m", "uV/m", "mV/m", "V/m"),
        )
        if unit is None:
            unit = "V/m"
            provenance = _with_assumption(
                provenance, "RF field unit was not declared; assumed V/m."
            )
        canonical = _canonical_unit(unit)
        if canonical not in _FIELD_TO_V_M:
            raise RydbergCSVError(f"Unsupported RF field unit: {unit!r}.")
        drive = parsed.x * _FIELD_TO_V_M[canonical]
        output_unit = "V/m"
    elif quantity == "rabi_frequency":
        unit = drive_unit or _unit_from_label(
            parsed.x_label, ("GHz", "MHz", "kHz", "Hz")
        )
        if unit is None:
            unit = "MHz"
            provenance = _with_assumption(
                provenance, "RF Rabi-frequency unit was not declared; assumed MHz."
            )
        drive = parsed.x * _frequency_factor(unit, "mhz")
        output_unit = "MHz"
    elif quantity == "source_power":
        unit = drive_unit or _unit_from_label(parsed.x_label, ("dBm", "mW", "W"))
        if unit is None:
            unit = "dBm"
            provenance = _with_assumption(
                provenance, "RF source-power unit was not declared; assumed dBm."
            )
        canonical = _canonical_unit(unit)
        if canonical == "dbm":
            drive = parsed.x
            output_unit = "dBm"
        elif canonical == "mw":
            drive = parsed.x * 1.0e-3
            output_unit = "W"
        elif canonical == "w":
            drive = parsed.x
            output_unit = "W"
        else:
            raise RydbergCSVError(f"Unsupported source-power unit: {unit!r}.")
    elif quantity == "generic":
        drive = parsed.x
        output_unit = str(drive_unit or "arb.")
        if drive_unit is None:
            provenance = _with_assumption(
                provenance, "Generic RF drive has arbitrary units."
            )
    else:
        raise RydbergCSVError(
            "drive_quantity must be 'field', 'rabi_frequency', "
            "'source_power', or 'generic'."
        )
    return RFSweepTrace(
        drive=_readonly(drive),
        response=parsed.y,
        drive_quantity=quantity,
        drive_unit=output_unit,
        response_unit=str(response_unit),
        provenance=provenance,
    )


def _spectral_scale(unit: str, quantity: str, spectrum_kind: str) -> tuple[float, str]:
    canonical = _canonical_unit(unit).replace("(", "").replace(")", "")
    kind = spectrum_kind
    if kind == "asd":
        suffixes = ("/sqrthz", "/hz^0.5", "/hz0.5")
    else:
        suffixes = ("^2/hz", "2/hz")
    numerator = canonical
    for suffix in suffixes:
        if numerator.endswith(suffix):
            numerator = numerator[: -len(suffix)]
            break

    if quantity == "voltage":
        base_unit = "V"
        amplitude_scales = {"v": 1.0, "mv": 1e-3, "uv": 1e-6, "nv": 1e-9}
    elif quantity == "field":
        base_unit = "V/m"
        amplitude_scales = _FIELD_TO_V_M
    elif quantity == "power":
        base_unit = "W"
        amplitude_scales = {"w": 1.0, "mw": 1e-3, "uw": 1e-6, "nw": 1e-9}
    elif quantity == "generic":
        return 1.0, "arb."
    else:
        raise RydbergCSVError(
            "PSD quantity must be 'voltage', 'field', 'power', or 'generic'."
        )
    if numerator not in amplitude_scales:
        raise RydbergCSVError(f"Unsupported {quantity} spectral unit: {unit!r}.")
    amplitude_scale = amplitude_scales[numerator]
    return (
        amplitude_scale if kind == "asd" else amplitude_scale**2,
        base_unit,
    )


def _infer_spectrum_kind(label: str) -> str | None:
    canonical = _canonical_unit(label)
    if "asd" in canonical or "sqrthz" in canonical or "hz^0.5" in canonical:
        return "asd"
    if "psd" in canonical or "2/hz" in canonical or "^2/hz" in canonical:
        return "psd"
    return None


def _default_spectral_unit(quantity: str, kind: str) -> str:
    numerator = {"voltage": "V", "field": "V/m", "power": "W", "generic": "arb."}[quantity]
    return f"{numerator}/sqrt(Hz)" if kind == "asd" else f"{numerator}^2/Hz"


def load_psd_csv(
    source: bytes | bytearray | memoryview | str | Path,
    *,
    source_name: str | None = None,
    frequency_column: int | str | None = None,
    spectrum_column: int | str | None = None,
    frequency_unit: str | None = None,
    spectrum_kind: str | None = None,
    spectral_unit: str | None = None,
    quantity: str = "voltage",
) -> PSDTrace:
    """Load PSD/ASD data and return both representations in base SI units."""

    parsed = _parse_columns(
        source,
        source_name=source_name,
        x_aliases=("frequency", "freq", "offset"),
        y_aliases=("spectral density", "psd", "asd", "noise", "spectrum"),
        x_column=frequency_column,
        y_column=spectrum_column,
    )
    provenance = parsed.provenance
    freq_unit = frequency_unit or _unit_from_label(
        parsed.x_label, ("GHz", "MHz", "kHz", "Hz")
    )
    if freq_unit is None:
        freq_unit = "Hz"
        provenance = _with_assumption(
            provenance, "Spectrum frequency unit was not declared; assumed Hz."
        )

    kind = (
        str(spectrum_kind).strip().lower()
        if spectrum_kind is not None
        else _infer_spectrum_kind(parsed.y_label)
    )
    if kind is None:
        kind = "asd"
        provenance = _with_assumption(
            provenance, "Spectrum kind was not declared; assumed ASD."
        )
    if kind not in {"asd", "psd"}:
        raise RydbergCSVError("spectrum_kind must be 'asd' or 'psd'.")
    quantity_name = str(quantity).strip().lower()
    if quantity_name not in {"voltage", "field", "power", "generic"}:
        raise RydbergCSVError(
            "quantity must be 'voltage', 'field', 'power', or 'generic'."
        )
    unit = spectral_unit
    if unit is None:
        # Header parsing for arbitrary spectral units is deliberately conservative.
        # Recognize common complete tokens; otherwise keep the assumption explicit.
        candidates = []
        prefixes = ("n", "u", "m", "")
        numerator = {
            "voltage": "V",
            "field": "V/m",
            "power": "W",
            "generic": "arb.",
        }[quantity_name]
        if quantity_name == "field":
            numerator_choices = tuple(prefix + "V/cm" for prefix in prefixes) + tuple(
                prefix + "V/m" for prefix in prefixes
            )
        else:
            numerator_choices = tuple(prefix + numerator for prefix in prefixes)
        for choice in numerator_choices:
            candidates.append(
                f"{choice}/sqrt(Hz)" if kind == "asd" else f"{choice}^2/Hz"
            )
        unit = _unit_from_label(parsed.y_label, tuple(candidates))
    if unit is None:
        unit = _default_spectral_unit(quantity_name, kind)
        provenance = _with_assumption(
            provenance,
            f"Spectral-density unit was not declared; assumed {unit}.",
        )
    scale, base_unit = _spectral_scale(unit, quantity_name, kind)
    values = parsed.y * scale
    if np.any(values < 0.0):
        raise RydbergCSVError("PSD/ASD values must be non-negative.")
    if kind == "asd":
        asd = values
        psd = values**2
    else:
        psd = values
        asd = np.sqrt(values)
    return PSDTrace(
        frequency_hz=_readonly(parsed.x * _frequency_factor(freq_unit, "hz")),
        asd=_readonly(asd),
        psd=_readonly(psd),
        quantity=quantity_name,
        base_unit=base_unit,
        input_spectrum_kind=kind,
        provenance=provenance,
    )


__all__ = [
    "EITTrace",
    "MAX_RYDBERG_CSV_BYTES",
    "MAX_RYDBERG_CSV_ROWS",
    "MIN_TRACE_POINTS",
    "PSDTrace",
    "RFSweepTrace",
    "RydbergCSVError",
    "TraceProvenance",
    "load_eit_csv",
    "load_psd_csv",
    "load_rf_sweep_csv",
]
