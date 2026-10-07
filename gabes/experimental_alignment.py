"""Display alignment of an imported spectrum, without fitting atomic parameters.

An arbitrary time/voltage trace does not establish a frequency or absolute
transmission calibration. Span framing and optional shape alignment are visual
aids; measured dark/reference calibration remains authoritative.
"""
from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np

from .experimental_csv import ExperimentalCSV, transform_x


@dataclass(frozen=True)
class AxisAlignment:
    """Affine display transforms; ``scale`` is the complete positive X scale."""

    scale: float = 1.0
    shift: float = 0.0
    reverse: bool = False
    signal_scale: float = 1.0
    signal_offset: float = 0.0

    def __post_init__(self):
        if not np.isfinite(self.scale) or self.scale <= 0.0:
            raise ValueError("Alignment scale must be finite and positive.")
        if not np.isfinite(self.shift):
            raise ValueError("Alignment shift must be finite.")
        if not np.isfinite(self.signal_scale) or self.signal_scale <= 0.0:
            raise ValueError("Signal scale must be finite and positive.")
        if not np.isfinite(self.signal_offset):
            raise ValueError("Signal offset must be finite.")

    def transform_x(self, x):
        return transform_x(x, scale=self.scale, shift=self.shift,
                           reverse=self.reverse)

    def transform_y(self, transmission, *, invert=False):
        values = np.asarray(transmission, dtype=float)
        if not np.all(np.isfinite(values)):
            raise ValueError("Transmission must contain only finite values.")
        values = 1.0 - values if invert else values
        with np.errstate(over="ignore", invalid="ignore"):
            transformed = self.signal_scale * values + self.signal_offset
        if not np.all(np.isfinite(transformed)):
            raise ValueError("Transformed transmission must remain finite.")
        return transformed


@dataclass(frozen=True)
class AxisScaleCorrection:
    """Frozen decade base that keeps a separate manual multiplier useful."""

    multiplier: float = 1.0
    target_manual_scale: float = float("nan")
    experimental_span: float = float("nan")
    theoretical_span: float = float("nan")
    valid: bool = False
    reason: str = ""


@dataclass(frozen=True)
class AlignmentFitResult:
    alignment: AxisAlignment
    success: bool
    improved: bool
    rmse: float
    initial_rmse: float
    overlap_fraction: float
    points_used: int
    message: str


def _axis_geometry(x, name):
    axis = np.asarray(x, dtype=float)
    if axis.ndim != 1 or axis.size < 2 or not np.all(np.isfinite(axis)):
        raise ValueError(f"{name} must contain at least two finite axis values.")
    low, high = float(np.min(axis)), float(np.max(axis))
    span = high - low
    if not np.isfinite(span) or span <= 0.0:
        raise ValueError(f"{name} must have a finite positive span.")
    return axis, span, 0.5 * low + 0.5 * high


def frame_alignment(x, xlim, *, fill=0.9, reverse=False):
    """Put the imported span inside a plot using the same pivot as transform_x."""
    _axis, raw_span, raw_midpoint = _axis_geometry(x, "Imported axis")
    _limits, plot_span, plot_midpoint = _axis_geometry(xlim, "Plot limits")
    fill = float(fill)
    if not np.isfinite(fill) or not 0.0 < fill <= 1.0:
        raise ValueError("fill must be finite and in (0, 1].")
    scale = (fill * plot_span) / raw_span
    result = AxisAlignment(scale=scale, shift=plot_midpoint - raw_midpoint,
                           reverse=bool(reverse))
    result.transform_x(x)  # Reject otherwise finite transforms that overflow.
    return result


def estimate_axis_scale_correction(x, xlim, *, manual_min=0.1, manual_max=3.0):
    """Estimate span-only decade preconditioning; never infer spectral features.

    Store the returned correction until a new trace or explicit reframe. Calling
    it on a slider event and compounding it with the prior transform is incorrect.
    """
    try:
        _x, raw_span, _mid = _axis_geometry(x, "Imported axis")
        _limits, plot_span, _mid = _axis_geometry(xlim, "Plot limits")
        lower, upper = float(manual_min), float(manual_max)
        if not (np.isfinite(lower) and np.isfinite(upper)
                and 0.0 < lower < upper):
            raise ValueError("Manual scale bounds must be finite and positive.")
        log_ratio = math.log10(plot_span) - math.log10(raw_span)
        if math.log10(lower) <= log_ratio <= math.log10(upper):
            exponent = 0
        else:
            exponent = math.floor(log_ratio - math.log10(upper)) + 1
        multiplier = float(10.0 ** exponent)
        target = float(10.0 ** (log_ratio - exponent))
        effective = multiplier * target
        if not (np.isfinite(multiplier) and multiplier > 0.0
                and lower <= target <= upper
                and np.isfinite(effective) and effective > 0.0):
            raise ValueError("Axis spans cannot be represented by a usable decade scale.")
        return AxisScaleCorrection(multiplier, target, raw_span, plot_span, True)
    except (TypeError, ValueError, OverflowError) as exc:
        return AxisScaleCorrection(reason=str(exc))


def _robust_affine(source, target):
    """Positive affine signal fit with iteratively reweighted outlier rejection."""
    weight = np.ones(source.size)
    slope, offset = 1.0, 0.0
    for _iteration in range(4):
        total = float(np.sum(weight))
        source_mean = float(np.dot(weight, source) / total)
        target_mean = float(np.dot(weight, target) / total)
        centered = source - source_mean
        variance = float(np.dot(weight, centered * centered))
        if variance <= np.finfo(float).tiny:
            return None
        slope = float(np.dot(weight * centered, target - target_mean) / variance)
        if not np.isfinite(slope) or slope <= 0.0:
            return None
        offset = target_mean - slope * source_mean
        residual = slope * source + offset - target
        sigma = 1.4826 * float(np.median(np.abs(residual - np.median(residual))))
        if sigma <= 32.0 * np.finfo(float).eps:
            break
        weight = np.minimum(1.0, 1.5 * sigma / np.maximum(np.abs(residual), sigma * 1e-12))
    return slope, offset


def _rms(values):
    """Compute RMS without overflowing the square of a finite bad initial gain."""
    magnitude = float(np.max(np.abs(values)))
    if not np.isfinite(magnitude):
        return float("inf")
    if magnitude == 0.0:
        return 0.0
    return magnitude * float(np.sqrt(np.mean((values / magnitude)**2)))


def _nelder_mead(objective, seed, bounds, iterations=110):
    """Small deterministic bounded two-dimensional derivative-free search."""
    lower, upper = np.asarray(bounds[0]), np.asarray(bounds[1])
    simplex = np.asarray([seed, seed + np.array([0.08, 0.0]),
                          seed + np.array([0.0, 0.04])])
    simplex = np.clip(simplex, lower, upper)
    values = np.array([objective(point)[0] for point in simplex])
    for _iteration in range(iterations):
        order = np.argsort(values)
        simplex, values = simplex[order], values[order]
        if np.max(np.abs(simplex[1:] - simplex[0])) < 2e-7:
            break
        center = np.mean(simplex[:2], axis=0)
        reflected = np.clip(2.0 * center - simplex[2], lower, upper)
        reflection = objective(reflected)[0]
        if reflection < values[0]:
            expanded = np.clip(3.0 * center - 2.0 * simplex[2], lower, upper)
            expansion = objective(expanded)[0]
            simplex[2], values[2] = ((expanded, expansion) if expansion < reflection
                                     else (reflected, reflection))
        elif reflection < values[1]:
            simplex[2], values[2] = reflected, reflection
        else:
            contracted = np.clip(0.5 * (center + simplex[2]), lower, upper)
            contraction = objective(contracted)[0]
            if contraction < values[2]:
                simplex[2], values[2] = contracted, contraction
            else:
                simplex[1:] = 0.5 * (simplex[1:] + simplex[0])
                values[1:] = [objective(point)[0] for point in simplex[1:]]
    return simplex[int(np.argmin(values))]


def fit_trace_alignment(trace: ExperimentalCSV, theory_x, theory_transmission, *,
                        initial=None, invert=False, fit_signal=True, reverse=None):
    """Fit only affine display alignment to a fixed theoretical transmission.

    Search is bounded around span framing and the current alignment. Sweep sign
    is explicit: ``reverse=None`` preserves the initial sign. A relative trace
    may receive a positive signal scale/offset; absolute calibration is preserved
    exactly. A successful result is a shape match, not evidence for a physical
    temperature, isotope ratio, or frequency calibration.
    """
    raw_x, raw_span, raw_midpoint = _axis_geometry(trace.detuning, "Imported axis")
    target_x, target_span, target_midpoint = _axis_geometry(theory_x, "Theory axis")
    measured = np.asarray(trace.transmission, dtype=float)
    theory = np.asarray(theory_transmission, dtype=float)
    if (measured.shape != raw_x.shape or theory.shape != target_x.shape
            or not np.all(np.isfinite(measured)) or not np.all(np.isfinite(theory))):
        raise ValueError("Each axis requires a same-length finite transmission array.")
    if raw_x.size < 12 or target_x.size < 12:
        raise ValueError("Spectrum alignment requires at least 12 points per trace.")
    order = np.argsort(target_x)
    target_x, theory = target_x[order], theory[order]
    if np.any(np.diff(target_x) <= 0.0):
        raise ValueError("Theory detuning values must be distinct.")
    measured_contrast = float(np.ptp(measured))
    theory_contrast = float(np.ptp(theory))
    tolerance = 128.0 * np.finfo(float).eps
    if measured_contrast <= tolerance * max(1.0, float(np.max(np.abs(measured)))):
        raise ValueError("A flat imported signal cannot establish spectrum alignment.")
    if theory_contrast <= tolerance * max(1.0, float(np.max(np.abs(theory)))):
        raise ValueError("A flat theoretical spectrum cannot establish alignment.")
    absolute = bool(trace.is_absolute_transmission)
    if absolute and invert:
        raise ValueError("Absolute transmission calibration cannot be inverted.")
    initial = initial or frame_alignment(raw_x, target_x)
    if not isinstance(initial, AxisAlignment):
        raise TypeError("initial must be an AxisAlignment.")
    direction = initial.reverse if reverse is None else bool(reverse)
    if absolute:
        initial = AxisAlignment(initial.scale, initial.shift, direction)
    elif direction != initial.reverse:
        initial = AxisAlignment(initial.scale, initial.shift, direction,
                                initial.signal_scale, initial.signal_offset)
    fit_signal = bool(fit_signal) and not absolute
    # Deterministic decimation bounds button latency on large oscilloscope files.
    indices = np.linspace(0, raw_x.size - 1, min(raw_x.size, 1200), dtype=int)
    normalized_x = (raw_x[indices] - raw_midpoint) / raw_span
    if direction:
        normalized_x = -normalized_x
    source = measured[indices]
    source = 1.0 - source if invert else source
    theory_axis = (target_x - target_midpoint) / target_span
    minimum_overlap = 0.65

    def evaluate(point, *, signal_fit=fit_signal):
        log_width, center = point
        if not np.all(np.isfinite(point)) or log_width > math.log(np.finfo(float).max):
            return float("inf"), None
        with np.errstate(over="ignore", invalid="ignore"):
            mapped = center + math.exp(float(log_width)) * normalized_x
        if not np.all(np.isfinite(mapped)):
            return float("inf"), None
        mask = (mapped >= theory_axis[0]) & (mapped <= theory_axis[-1])
        overlap = float(np.mean(mask))
        count = int(np.count_nonzero(mask))
        if overlap < minimum_overlap or count < 12:
            return float("inf"), None
        desired = np.interp(mapped[mask], theory_axis, theory)
        if float(np.ptp(desired)) < 0.15 * theory_contrast:
            return float("inf"), None
        observed = source[mask]
        signal = (_robust_affine(observed, desired) if signal_fit
                  else (initial.signal_scale, initial.signal_offset))
        if signal is None:
            return float("inf"), None
        signal_scale, signal_offset = signal
        with np.errstate(over="ignore", invalid="ignore"):
            residual = signal_scale * observed + signal_offset - desired
            squared = residual * residual
        if not np.all(np.isfinite(squared)):
            return float("inf"), None
        # Capped residual loss resists isolated spikes. Dividing by local theory
        # variance prevents a free signal scale from preferring a flat baseline.
        cap = max((0.15 * theory_contrast)**2,
                  4.0 * float(np.median(squared)))
        robust_loss = float(np.mean(np.minimum(squared, cap)))
        variance = float(np.var(desired))
        if not np.isfinite(variance) or variance <= np.finfo(float).tiny:
            return float("inf"), None
        score = robust_loss / variance + 0.08 * (1.0 - overlap)**2
        if not np.isfinite(score):
            return float("inf"), None
        return score, (signal_scale, signal_offset,
                       float(np.sqrt(np.mean(squared))), overlap, count)

    initial_log_width = math.log(initial.scale) + math.log(raw_span) - math.log(target_span)
    initial_width = (math.exp(initial_log_width)
                     if initial_log_width < math.log(np.finfo(float).max)
                     else float("inf"))
    initial_center = (raw_midpoint + initial.shift - target_midpoint) / target_span
    initial_point = np.array([initial_log_width, initial_center])
    initial_score, initial_stats = evaluate(initial_point, signal_fit=False)
    initial_rmse = initial_stats[2] if initial_stats is not None else float("inf")
    widths = np.geomspace(0.225, 3.6, 37)
    centers = np.linspace(-0.75, 0.75, 65)
    # A sensible manual starting point contributes its own bounded neighbourhood.
    if 0.025 <= initial_width <= 8.0 and abs(initial_center) <= 2.0:
        widths = np.unique(np.concatenate((widths, initial_width * np.geomspace(0.6, 1.6, 13))))
        centers = np.unique(np.concatenate((centers, initial_center + np.linspace(-0.25, 0.25, 17))))
    lower = np.array([math.log(float(widths[0])), float(centers[0])])
    upper = np.array([math.log(float(widths[-1])), float(centers[-1])])
    candidates = []
    current_score, current_stats = evaluate(initial_point)
    if np.isfinite(current_score):
        candidates.append((current_score, initial_point))
    for width in widths:
        log_width = math.log(float(width))
        for center in centers:
            point = np.array([log_width, center])
            score = evaluate(point)[0]
            if np.isfinite(score):
                candidates.append((score, point))
    if not candidates:
        return AlignmentFitResult(initial, False, False, initial_rmse,
                                  initial_rmse, 0.0, 0,
                                  "No informative overlap; adjust the scan span or sweep direction.")
    candidates.sort(key=lambda item: item[0])
    # Refine separated basins so an accidentally aligned single dip does not
    # prevent finding the multi-feature solution.
    seeds = []
    for _score, point in candidates:
        if all(np.linalg.norm(point - other) > 0.15 for other in seeds):
            seeds.append(point)
        if len(seeds) == 5:
            break
    best_score, best_point = candidates[0]
    for seed in seeds:
        point = _nelder_mead(evaluate, seed, (lower, upper))
        score = evaluate(point)[0]
        if score < best_score:
            best_score, best_point = score, point
    _score, stats = evaluate(best_point)
    signal_scale, signal_offset, rmse, overlap, count = stats
    scale = math.exp(float(best_point[0])) * target_span / raw_span
    shift = target_midpoint + float(best_point[1]) * target_span - raw_midpoint
    alignment = AxisAlignment(scale, shift, direction, signal_scale, signal_offset)
    mapped_full = alignment.transform_x(raw_x)
    full_mask = (mapped_full >= target_x[0]) & (mapped_full <= target_x[-1])
    count = int(np.count_nonzero(full_mask))
    overlap = float(np.mean(full_mask))
    desired_full = np.interp(mapped_full[full_mask], target_x, theory)
    aligned_full = alignment.transform_y(measured, invert=invert)
    rmse = _rms(aligned_full[full_mask] - desired_full)
    # Search is decimated; report statistics on every overlapping imported point.
    try:
        mapped_initial = initial.transform_x(raw_x)
        initial_mask = (mapped_initial >= target_x[0]) & (mapped_initial <= target_x[-1])
        if np.count_nonzero(initial_mask):
            desired_initial = np.interp(mapped_initial[initial_mask], target_x, theory)
            signal_initial = initial.transform_y(measured, invert=invert)
            initial_rmse = _rms(signal_initial[initial_mask] - desired_initial)
    except ValueError:
        initial_rmse = float("inf")
    # Search uses capped residuals to tolerate isolated spikes; acceptance must
    # use the uncapped residual. Otherwise unrelated noise can be flattened by
    # a nearly zero Y slope and score well against a sparse absorption spectrum.
    desired_std = float(np.std(desired_full))
    normalized_rmse = rmse / desired_std if desired_std > 0.0 else float("inf")
    centered_observed = aligned_full[full_mask] - np.mean(aligned_full[full_mask])
    centered_desired = desired_full - np.mean(desired_full)
    correlation_norm = float(np.linalg.norm(centered_observed) * np.linalg.norm(centered_desired))
    correlation = (float(np.dot(centered_observed, centered_desired)) / correlation_norm
                   if correlation_norm > 0.0 else 0.0)
    success = bool(best_score < 0.5 and overlap >= minimum_overlap
                   and np.isfinite(normalized_rmse) and normalized_rmse <= 0.35
                   and np.isfinite(correlation) and correlation >= 0.85)
    improved = bool(success and (not np.isfinite(initial_score)
                    or best_score < initial_score * (1.0 - 1e-5)))
    if not success:
        message = "Spectral shapes did not match reliably; check sweep direction and scan coverage."
    elif improved:
        message = "Affine spectrum alignment improved; inspect the match before interpreting it."
    else:
        message = "Current alignment already matches this spectrum; no material improvement."
    return AlignmentFitResult(alignment, success, improved, rmse, initial_rmse,
                              overlap, count, message)


__all__ = ["AxisAlignment", "AxisScaleCorrection", "AlignmentFitResult",
           "frame_alignment", "estimate_axis_scale_correction", "fit_trace_alignment"]
