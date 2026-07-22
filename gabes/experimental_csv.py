"""Import and robustly correct two-column oscilloscope CSV traces.

Only the first two columns are meaningful: column A is detuning and column B
is the detector signal, both in arbitrary units. The helpers here depend on
the Python standard library and NumPy only.
"""

from __future__ import annotations

import csv
import io
from dataclasses import dataclass

import numpy as np


MAX_FILE_BYTES = 10 * 1024 * 1024
MAX_CSV_ROWS = 500_000
MIN_UNIQUE_POINTS = 5
HAMPEL_THRESHOLD_MAD = 4.5


class ExperimentalCSVError(ValueError):
    """Raised when an experimental CSV cannot be imported safely."""


@dataclass(frozen=True)
class CSVImportDiagnostics:
    """Counts and structural information collected while importing a CSV."""

    encoding: str
    file_size_bytes: int
    total_rows: int
    valid_rows: int
    blank_rows: int
    nonnumeric_rows: int
    nonfinite_rows: int
    extra_column_rows: int
    unique_points: int
    duplicate_rows_merged: int
    segment_count: int
    typical_x_spacing: float
    largest_x_gap: float
    gap_threshold: float

    @property
    def ignored_rows(self) -> int:
        return self.blank_rows + self.nonnumeric_rows + self.nonfinite_rows


@dataclass(frozen=True)
class CorrectionDiagnostics:
    """Numerical choices and quality warnings from signal correction."""

    denoise_enabled: bool
    noise_sigma: float
    hampel_replacements: int
    smoothing_windows: tuple[int, ...]
    floor: float
    ceiling: float
    contrast: float
    contrast_to_noise: float
    relative_contrast: float
    floor_support: int
    ceiling_support: int
    sparse_extrema: bool
    low_contrast: bool
    warning: str | None

    @property
    def smoothing_window(self) -> int:
        return max(self.smoothing_windows, default=1)


@dataclass(frozen=True)
class ExperimentalCSV:
    """A sorted, de-duplicated experimental trace and all correction stages."""

    detuning: np.ndarray
    raw_signal: np.ndarray
    denoised_signal: np.ndarray
    transmission: np.ndarray
    import_diagnostics: CSVImportDiagnostics
    correction_diagnostics: CorrectionDiagnostics

    def transformed_detuning(
        self,
        *,
        scale: float = 1.0,
        shift: float = 0.0,
        reverse: bool = False,
        pivot: float | None = None,
    ) -> np.ndarray:
        return transform_x(
            self.detuning,
            scale=scale,
            shift=shift,
            reverse=reverse,
            pivot=pivot,
        )


def load_experimental_csv(
    data: bytes | bytearray | memoryview,
    *,
    denoise: bool = True,
) -> ExperimentalCSV:
    """Parse and correct an oscilloscope CSV supplied as bytes.

    Text is decoded as BOM-marked UTF-16, UTF-8, or CP949, with a byte-preserving
    fallback so unsupported metadata in ignored columns cannot reject valid A/B
    data. Only columns A and B are inspected. Blank, textual, NaN, and infinite
    rows are counted and ignored. Equal A values are merged by the B median.

    Denoising uses a 4.5-MAD Hampel filter followed by adaptive, x-aware local
    quadratic regression. Large gaps split the trace before correction. The
    chosen signal is robustly mapped to transmission in [0, 1].
    """

    raw_bytes = _coerce_bytes(data)
    text, encoding = _decode_csv(raw_bytes)
    source_x, source_y, counts = _parse_first_two_columns(text)
    detuning, raw_signal = _sort_and_merge(source_x, source_y)
    if detuning.size < MIN_UNIQUE_POINTS:
        raise ExperimentalCSVError(
            f"CSV must contain at least {MIN_UNIQUE_POINTS} unique finite "
            f"detuning values; found {detuning.size}."
        )

    segments, typical_spacing, largest_gap, gap_threshold = _find_segments(
        detuning
    )
    denoised_signal, noise_sigma, replacements, windows = _denoise_signal(
        detuning, raw_signal, segments, enabled=bool(denoise)
    )
    transmission, correction = _normalize_transmission(
        denoised_signal,
        denoise_enabled=bool(denoise),
        noise_sigma=noise_sigma,
        hampel_replacements=replacements,
        smoothing_windows=windows,
    )
    imported = CSVImportDiagnostics(
        encoding=encoding,
        file_size_bytes=len(raw_bytes),
        total_rows=counts["total"],
        valid_rows=counts["valid"],
        blank_rows=counts["blank"],
        nonnumeric_rows=counts["nonnumeric"],
        nonfinite_rows=counts["nonfinite"],
        extra_column_rows=counts["extra"],
        unique_points=int(detuning.size),
        duplicate_rows_merged=int(source_x.size - detuning.size),
        segment_count=len(segments),
        typical_x_spacing=typical_spacing,
        largest_x_gap=largest_gap,
        gap_threshold=gap_threshold,
    )
    return ExperimentalCSV(
        detuning=_readonly(detuning),
        raw_signal=_readonly(raw_signal),
        denoised_signal=_readonly(denoised_signal),
        transmission=_readonly(transmission),
        import_diagnostics=imported,
        correction_diagnostics=correction,
    )


def transform_x(
    x: np.ndarray,
    *,
    scale: float = 1.0,
    shift: float = 0.0,
    reverse: bool = False,
    pivot: float | None = None,
) -> np.ndarray:
    """Apply ``pivot + sign*scale*(x-pivot) + shift`` to an x axis."""

    values = np.asarray(x, dtype=float)
    if values.ndim != 1 or values.size == 0:
        raise ValueError("x must be a non-empty one-dimensional array.")
    if not np.all(np.isfinite(values)):
        raise ValueError("x must contain only finite values.")
    scale = float(scale)
    shift = float(shift)
    if not np.isfinite(scale) or scale <= 0.0:
        raise ValueError("scale must be a finite positive number.")
    if not np.isfinite(shift):
        raise ValueError("shift must be finite.")
    if pivot is None:
        pivot_value = (
            0.5 * float(np.min(values)) + 0.5 * float(np.max(values))
        )
    else:
        pivot_value = float(pivot)
        if not np.isfinite(pivot_value):
            raise ValueError("pivot must be finite.")
    signed_scale = -scale if reverse else scale
    with np.errstate(over="ignore", invalid="ignore"):
        transformed = pivot_value + signed_scale * (values - pivot_value) + shift
    if not np.all(np.isfinite(transformed)):
        raise ValueError("transformed x values must remain finite.")
    return transformed


def _coerce_bytes(data: bytes | bytearray | memoryview) -> bytes:
    if not isinstance(data, (bytes, bytearray, memoryview)):
        raise TypeError("CSV input must be bytes.")
    result = bytes(data)
    if len(result) > MAX_FILE_BYTES:
        raise ExperimentalCSVError(
            f"CSV is larger than the {MAX_FILE_BYTES // (1024 * 1024)} MiB "
            "import limit."
        )
    return result


def _decode_csv(data: bytes) -> tuple[str, str]:
    if data.startswith((b"\xff\xfe", b"\xfe\xff")):
        try:
            return data.decode("utf-16"), "utf-16"
        except UnicodeDecodeError as exc:
            raise ExperimentalCSVError("Malformed UTF-16 CSV text.") from exc

    for encoding in ("utf-8-sig", "utf-8", "cp949"):
        try:
            return data.decode(encoding), encoding
        except UnicodeDecodeError:
            continue
    # A/B numeric fields and CSV delimiters are ASCII. Latin-1's one-to-one
    # byte mapping lets us safely ignore an unsupported byte in oscilloscope
    # metadata or later columns instead of rejecting otherwise valid A/B data.
    return data.decode("latin-1"), "latin-1"


def _parse_first_two_columns(
    text: str,
) -> tuple[np.ndarray, np.ndarray, dict[str, int]]:
    xs: list[float] = []
    ys: list[float] = []
    counts = {
        "total": 0,
        "valid": 0,
        "blank": 0,
        "nonnumeric": 0,
        "nonfinite": 0,
        "extra": 0,
    }
    try:
        reader = csv.reader(io.StringIO(text, newline=""))
        for row in reader:
            counts["total"] += 1
            if counts["total"] > MAX_CSV_ROWS:
                raise ExperimentalCSVError(
                    f"CSV exceeds the {MAX_CSV_ROWS:,}-row import limit."
                )
            if len(row) > 2:
                counts["extra"] += 1
            if not row or all(not cell.strip() for cell in row):
                counts["blank"] += 1
                continue
            if len(row) < 2:
                counts["nonnumeric"] += 1
                continue
            a, b = row[0].strip(), row[1].strip()
            if not a or not b:
                counts["blank"] += 1
                continue
            try:
                x_value = float(a)
                y_value = float(b)
            except ValueError:
                counts["nonnumeric"] += 1
                continue
            if not np.isfinite(x_value) or not np.isfinite(y_value):
                counts["nonfinite"] += 1
                continue
            xs.append(x_value)
            ys.append(y_value)
            counts["valid"] += 1
    except csv.Error as exc:
        raise ExperimentalCSVError(f"Malformed CSV: {exc}") from exc
    return np.asarray(xs, dtype=float), np.asarray(ys, dtype=float), counts


def _sort_and_merge(x: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    if x.size == 0:
        return x.copy(), y.copy()
    order = np.argsort(x, kind="mergesort")
    sorted_x = x[order]
    sorted_y = y[order]
    unique_x, starts, group_sizes = np.unique(
        sorted_x, return_index=True, return_counts=True
    )
    merged_y = sorted_y[starts].copy()
    for group_index in np.flatnonzero(group_sizes > 1):
        start = starts[group_index]
        stop = start + group_sizes[group_index]
        merged_y[group_index] = np.median(sorted_y[start:stop])
    return unique_x, merged_y


def _find_segments(
    x: np.ndarray,
) -> tuple[list[slice], float, float, float]:
    spacing = np.diff(x)
    typical = float(np.median(spacing))
    largest = float(np.max(spacing))
    q1, q3 = np.quantile(spacing, (0.25, 0.75))
    iqr = float(q3 - q1)
    threshold = float(max(8.0 * typical, q3 + 6.0 * iqr))
    boundaries = np.flatnonzero(spacing > threshold) + 1
    starts = np.concatenate(([0], boundaries))
    stops = np.concatenate((boundaries, [x.size]))
    return (
        [slice(int(start), int(stop)) for start, stop in zip(starts, stops)],
        typical,
        largest,
        threshold,
    )


def _denoise_signal(
    x: np.ndarray,
    y: np.ndarray,
    segments: list[slice],
    *,
    enabled: bool,
) -> tuple[np.ndarray, float, int, tuple[int, ...]]:
    segment_noise = [_estimate_noise(y[part]) for part in segments]
    segment_sizes = np.asarray([part.stop - part.start for part in segments])
    noise_sigma = float(
        np.average(np.asarray(segment_noise), weights=segment_sizes)
    )
    if not enabled:
        return y.copy(), noise_sigma, 0, tuple(1 for _ in segments)

    result = y.copy()
    replacements = 0
    windows: list[int] = []
    for part, segment_sigma in zip(segments, segment_noise):
        segment_y, replaced = _hampel_filter(result[part])
        replacements += replaced
        smoothed, window = _adaptive_local_quadratic(
            x[part], segment_y, noise_hint=segment_sigma
        )
        result[part] = smoothed
        windows.append(window)
    return result, noise_sigma, replacements, tuple(windows)


def _estimate_noise(y: np.ndarray) -> float:
    if y.size >= 5:
        differences = np.diff(y, n=2)
        divisor = np.sqrt(6.0)
    elif y.size >= 2:
        differences = np.diff(y)
        divisor = np.sqrt(2.0)
    else:
        return 0.0
    centered = differences - np.median(differences)
    return float(1.4826 * np.median(np.abs(centered)) / divisor)


def _hampel_filter(y: np.ndarray, half_window: int = 5) -> tuple[np.ndarray, int]:
    if y.size < 5:
        return y.copy(), 0
    half_window = min(half_window, max(2, (y.size - 1) // 2))
    width = 2 * half_window + 1
    padded = np.pad(y, (half_window, half_window), mode="reflect")
    windows = np.lib.stride_tricks.sliding_window_view(padded, width)
    local_median = np.median(windows, axis=1)
    local_mad = np.median(np.abs(windows - local_median[:, None]), axis=1)
    local_sigma = 1.4826 * local_mad
    noise = _estimate_noise(y)
    signal_scale = max(
        float(np.max(np.abs(y))), float(np.ptp(y)), np.finfo(float).tiny
    )
    sigma_floor = max(
        0.2 * noise, 32.0 * np.finfo(float).eps * signal_scale
    )
    threshold = HAMPEL_THRESHOLD_MAD * np.maximum(local_sigma, sigma_floor)
    candidate_mask = np.abs(y - local_median) > threshold
    if not np.any(candidate_mask):
        return y.copy(), 0

    # Be deliberately conservative around narrow spectral lines. Preserve a
    # multi-sample candidate run when its depth/height is resolved well above
    # both the measured noise and the full-trace scale. Weak adjacent candidates
    # are still corrected because quantized scope noise often arrives in pairs.
    mask = candidate_mask.copy()
    candidate_indices = np.flatnonzero(candidate_mask)
    candidate_runs = np.split(
        candidate_indices,
        np.flatnonzero(np.diff(candidate_indices) > 1) + 1,
    )
    feature_threshold = max(
        8.0 * noise,
        0.05 * float(np.ptp(y)),
        sigma_floor,
    )
    deviations = np.abs(y - local_median)
    singleton_feature_threshold = max(
        8.0 * noise,
        0.50 * float(np.ptp(y)),
        sigma_floor,
    )
    for run in candidate_runs:
        run_deviation = float(np.max(deviations[run]))
        is_resolved_run = run.size > 1 and run_deviation >= feature_threshold
        is_span_defining_singleton = (
            run.size == 1 and run_deviation >= singleton_feature_threshold
        )
        if is_resolved_run or is_span_defining_singleton:
            mask[run] = False
    if not np.any(mask):
        return y.copy(), 0

    result = y.copy()
    result[mask] = local_median[mask]

    # A single sampled Lamb dip can be mathematically indistinguishable from an
    # impulse.  If filtering would erase essentially all measured contrast,
    # preserve it and let the low-contrast diagnostic guide the user instead.
    original_span = float(np.ptp(y))
    filtered_span = float(np.ptp(result))
    span_tolerance = 64.0 * np.finfo(float).eps * signal_scale
    if (
        original_span > span_tolerance
        and filtered_span < 0.10 * original_span
    ):
        return y.copy(), 0
    return result, int(np.count_nonzero(mask))


def _adaptive_local_quadratic(
    x: np.ndarray,
    y: np.ndarray,
    *,
    noise_hint: float = 0.0,
) -> tuple[np.ndarray, int]:
    if y.size < 5:
        return y.copy(), 1
    # Hampel replacement can make quantized scope plateaus locally noiseless.
    # Keep the pre-Hampel estimate so quantization noise still receives a light,
    # feature-preserving smoothing pass.
    noise = max(_estimate_noise(y), float(noise_hint))
    low, high = np.quantile(y, (0.05, 0.95))
    dynamic_range = max(
        float(high - low), float(np.ptp(y)), np.finfo(float).tiny
    )
    ratio = noise / dynamic_range
    scale = max(
        float(np.max(np.abs(y))), dynamic_range, np.finfo(float).tiny
    )
    if noise <= 32.0 * np.finfo(float).eps * scale or ratio < 2e-5:
        return y.copy(), 1
    half_width = int(np.clip(np.ceil(2.0 + 60.0 * ratio), 2, 12))
    largest_odd_window = y.size if y.size % 2 else y.size - 1
    window = min(2 * half_width + 1, largest_odd_window)
    if window < 5:
        return y.copy(), 1
    return _local_quadratic(x, y, window), int(window)


def _local_quadratic(
    x: np.ndarray,
    y: np.ndarray,
    window: int,
    *,
    batch_size: int = 20_000,
) -> np.ndarray:
    """Unweighted local quadratic regression using actual x coordinates."""

    n_points = x.size
    starts = np.clip(
        np.arange(n_points, dtype=np.int64) - window // 2,
        0,
        n_points - window,
    )
    offsets = np.arange(window, dtype=np.int64)
    result = np.empty_like(y)
    for batch_start in range(0, n_points, batch_size):
        batch_stop = min(batch_start + batch_size, n_points)
        centers = np.arange(batch_start, batch_stop, dtype=np.int64)
        indices = starts[centers, None] + offsets[None, :]
        local_x = x[indices] - x[centers, None]
        radius = np.max(np.abs(local_x), axis=1)
        scaled_x = local_x / radius[:, None]
        local_y = y[indices]

        s0 = np.full(centers.size, float(window))
        s1 = np.sum(scaled_x, axis=1)
        s2 = np.sum(scaled_x**2, axis=1)
        s3 = np.sum(scaled_x**3, axis=1)
        s4 = np.sum(scaled_x**4, axis=1)
        b0 = np.sum(local_y, axis=1)
        b1 = np.sum(local_y * scaled_x, axis=1)
        b2 = np.sum(local_y * scaled_x**2, axis=1)

        normal = np.empty((centers.size, 3, 3), dtype=float)
        normal[:, 0, :] = np.column_stack((s0, s1, s2))
        normal[:, 1, :] = np.column_stack((s1, s2, s3))
        normal[:, 2, :] = np.column_stack((s2, s3, s4))
        rhs = np.column_stack((b0, b1, b2))
        try:
            coefficients = np.linalg.solve(normal, rhs[..., None])[..., 0]
            fitted = coefficients[:, 0]
            # Quadratic edge ringing can invent a new voltage floor/ceiling
            # around a very narrow resonance. Keep every fitted point within
            # the actual signal range of its local window.
            result[centers] = np.clip(
                fitted,
                np.min(local_y, axis=1),
                np.max(local_y, axis=1),
            )
        except np.linalg.LinAlgError:
            for output_index, t_values, signal_values in zip(
                centers, scaled_x, local_y
            ):
                design = np.column_stack(
                    (np.ones(t_values.size), t_values, t_values**2)
                )
                fitted = np.linalg.lstsq(
                    design, signal_values, rcond=None
                )[0][0]
                result[output_index] = np.clip(
                    fitted,
                    np.min(signal_values),
                    np.max(signal_values),
                )
    return result


def _normalize_transmission(
    signal: np.ndarray,
    *,
    denoise_enabled: bool,
    noise_sigma: float,
    hampel_replacements: int,
    smoothing_windows: tuple[int, ...],
) -> tuple[np.ndarray, CorrectionDiagnostics]:
    magnitude = max(float(np.max(np.abs(signal))), np.finfo(float).tiny)
    tolerance = 64.0 * np.finfo(float).eps * magnitude
    minimum = float(np.min(signal))
    maximum = float(np.max(signal))
    full_span = maximum - minimum
    if not np.isfinite(full_span) or full_span <= tolerance:
        raise ExperimentalCSVError(
            "Detector signal has no usable floor-to-ceiling contrast; "
            "check the CSV columns or acquire a wider scan."
        )

    # Estimate plateaus from neighbourhoods of the extrema instead of fixed
    # population quantiles.  A fixed 5th/95th-percentile calibration loses a
    # legitimate Lamb dip when it occupies only a few samples.  The noise-aware
    # band remains robust after Hampel/smoothing while retaining such features.
    extrema_band = min(
        max(2.0 * float(noise_sigma), 0.02 * full_span, tolerance),
        0.25 * full_span,
    )
    floor_values = signal[signal <= minimum + extrema_band]
    ceiling_values = signal[signal >= maximum - extrema_band]
    floor = float(np.median(floor_values))
    ceiling = float(np.median(ceiling_values))
    contrast = ceiling - floor
    if not np.isfinite(contrast) or contrast <= tolerance:
        raise ExperimentalCSVError(
            "Detector signal has no usable floor-to-ceiling contrast; "
            "check the CSV columns or acquire a wider scan."
        )

    normalized = np.clip((signal - floor) / contrast, 0.0, 1.0)
    if noise_sigma > tolerance:
        contrast_to_noise = float(contrast / noise_sigma)
    else:
        contrast_to_noise = float("inf")
    relative_contrast = float(
        contrast
        / max(abs(float(np.median(signal))), contrast, np.finfo(float).tiny)
    )
    floor_support = int(floor_values.size)
    ceiling_support = int(ceiling_values.size)
    minimum_support = max(3, int(np.ceil(0.005 * signal.size)))
    sparse_extrema = (
        floor_support < minimum_support or ceiling_support < minimum_support
    )
    low_contrast = contrast_to_noise < 8.0 or relative_contrast < 1e-4
    warnings = []
    if low_contrast:
        warnings.append(
            "Low detector contrast: automatic 0-1 calibration may amplify "
            "noise; verify the floor and ceiling before fitting."
        )
    if sparse_extrema:
        warnings.append(
            "Sparse calibration level: the estimated floor/ceiling is supported "
            f"by only {floor_support}/{ceiling_support} samples; inspect the "
            "unfiltered trace for an under-resolved line or remaining impulse."
        )
    warning = " ".join(warnings) or None
    diagnostics = CorrectionDiagnostics(
        denoise_enabled=denoise_enabled,
        noise_sigma=float(noise_sigma),
        hampel_replacements=int(hampel_replacements),
        smoothing_windows=smoothing_windows,
        floor=floor,
        ceiling=ceiling,
        contrast=float(contrast),
        contrast_to_noise=contrast_to_noise,
        relative_contrast=relative_contrast,
        floor_support=floor_support,
        ceiling_support=ceiling_support,
        sparse_extrema=bool(sparse_extrema),
        low_contrast=bool(low_contrast),
        warning=warning,
    )
    return normalized, diagnostics


def _readonly(values: np.ndarray) -> np.ndarray:
    result = np.asarray(values, dtype=float)
    result.setflags(write=False)
    return result


__all__ = [
    "CSVImportDiagnostics",
    "CorrectionDiagnostics",
    "ExperimentalCSV",
    "ExperimentalCSVError",
    "HAMPEL_THRESHOLD_MAD",
    "MAX_CSV_ROWS",
    "MAX_FILE_BYTES",
    "MIN_UNIQUE_POINTS",
    "load_experimental_csv",
    "transform_x",
]
