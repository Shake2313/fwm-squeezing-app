"""
Numeric line-width extractors shared by the schemes.

These operate on *sampled curves* (transmission, absorption, |ψ|², …). With one
labelled exception they are **shape-agnostic**: they walk a sampled peak to its
half-height and report the width, making no assumption about Lorentzian vs
Gaussian vs Voigt form. Consolidating them here removes the per-scheme copies of
`_window_fwhm` / `_peak_fwhm` / `_fwhm`.

The single model-dependent helper is `lorentz_fwhm`, which *fits* a Lorentzian
(1/signal linear in B²) and is therefore only valid for a Lorentzian-shaped
zero-field feature — kept separate on purpose.
"""
import math
from dataclasses import dataclass

import numpy as np


MIN_SUBDOPPLER_SAMPLES_PER_FWHM = 6.0
MIN_SUBDOPPLER_EDGE_CLEARANCE_FWHM = 1.0


@dataclass(frozen=True)
class SubdopplerFeature:
    """Resolution diagnostics for one sampled sub-Doppler feature.

    All frequency-valued fields use the units of the input axis.  A finite
    width can still be ``resolution-limited`` or ``scan-edge-limited``; callers
    should report its linewidth or lock slope only when :attr:`resolved` is true.
    """

    status: str
    reason: str
    center: float = float("nan")
    left_half_height: float = float("nan")
    right_half_height: float = float("nan")
    fwhm: float = float("nan")
    samples_per_fwhm: float = float("nan")
    scan_edge_distance: float = float("nan")
    amplitude: float = float("nan")

    @property
    def detected(self):
        return np.isfinite(self.center) and np.isfinite(self.fwhm)

    @property
    def resolved(self):
        return self.status == "resolved"


def window_fwhm(x, y, ic):
    """FWHM of a feature peaking at index `ic`, measured to the half-height
    between the peak `y[ic]` and the curve minimum.

    Shape-agnostic: works for a transparency window or a Lamb dip regardless of
    the underlying lineshape. `nan` if the feature is ill-defined.
    """
    peak = y[ic]
    floor = np.nanmin(y)
    if not np.isfinite(peak) or peak <= floor:
        return float("nan")
    thresh = 0.5 * (peak + floor)
    i = ic
    while i > 0 and y[i] >= thresh:
        i -= 1
    j = ic
    while j < y.size - 1 and y[j] >= thresh:
        j += 1
    return float(x[j] - x[i])


def fwhm_halfmax(x, y):
    """FWHM of the tallest peak by the half-maximum samples (no edge
    interpolation); `nan` if ill-defined. Use `fwhm_interp` for sub-sample
    accuracy."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if x.ndim != 1 or y.shape != x.shape or x.size < 2:
        return float("nan")
    finite = np.isfinite(y)
    if not np.any(finite):
        return float("nan")
    ic = int(np.nanargmax(y))
    pk = float(y[ic])
    if not np.isfinite(pk) or pk <= 0:
        return float("nan")
    above = finite & (y >= 0.5 * pk)
    lo = ic
    while lo > 0 and above[lo - 1]:
        lo -= 1
    hi = ic
    while hi < y.size - 1 and above[hi + 1]:
        hi += 1
    return float(x[hi] - x[lo]) if hi > lo else float("nan")


def fwhm_interp(x, y):
    """FWHM of the tallest peak with linear edge interpolation; `nan` if
    ill-defined."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if (x.ndim != 1 or y.shape != x.shape or x.size < 2
            or not np.any(np.isfinite(y))):
        return np.nan
    peak_index = int(np.nanargmax(y))
    peak = float(y[peak_index])
    if peak <= 0:
        return np.nan
    half = 0.5 * peak
    above = np.isfinite(y) & (y >= half)
    lo = peak_index
    while lo > 0 and above[lo - 1]:
        lo -= 1
    hi = peak_index
    while hi < y.size - 1 and above[hi + 1]:
        hi += 1

    def interp_edge(i0, i1):
        if i0 < 0 or i1 >= x.size or y[i1] == y[i0]:
            return x[max(min(i1, x.size - 1), 0)]
        return x[i0] + (half - y[i0]) * (x[i1] - x[i0]) / (y[i1] - y[i0])

    left = interp_edge(lo - 1, lo) if lo > 0 else x[lo]
    right = interp_edge(hi, hi + 1) if hi < x.size - 1 else x[hi]
    return float(max(right - left, 0.0))


def halfwidth_from_center(B, y, frac=0.5):
    """Half-width around zero for a central feature, robust to a broad pedestal.

    Walks right from the B≈0 sample to the `frac` crossing between the central
    amplitude and the wing background, with linear interpolation at the edge.
    Shape-agnostic. `nan` if ill-defined.
    """
    B = np.asarray(B)
    y = np.asarray(y)
    ic = int(np.argmin(np.abs(B)))
    bg = 0.5 * (y[0] + y[-1])
    amp = y[ic] - bg
    if abs(amp) < 1e-30:
        return float("nan")
    target = bg + frac * amp
    right = np.arange(ic, y.size)
    vals = (y[right] - target) * np.sign(amp)
    below = np.where(vals <= 0)[0]
    if below.size == 0 or below[0] == 0:
        return float("nan")
    j = right[below[0]]
    i = j - 1
    y0, y1 = y[i], y[j]
    if y1 == y0:
        return float(abs(B[j] - B[ic]))
    t = (target - y0) / (y1 - y0)
    return float(abs((B[i] + t * (B[j] - B[i])) - B[ic]))


def subdoppler_feature(x, signal, search_window=None, *, min_amplitude=1e-6):
    """Diagnose the sharpest sampled sub-Doppler feature.

    A running-median background is removed before selecting the largest
    residual.  ``search_window=(lo, hi)`` restricts only that peak selection;
    background estimation and half-height crossings still use the complete,
    predeclared scan.  Half-height edges are linearly interpolated.

    Resolution requires at least six sample intervals across the interpolated
    FWHM and one FWHM of clearance between the half-height edges and either scan
    edge.  A detected but undersampled feature retains its provisional numeric
    diagnostics; report its linewidth or lock slope only when its status is
    ``resolved``.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(signal, dtype=float)
    if x.ndim != 1 or y.ndim != 1 or x.size != y.size or x.size < 5:
        return SubdopplerFeature(
            "unresolved", "Need matching one-dimensional scans with at least five samples.")
    if not np.all(np.isfinite(x)) or not np.all(np.isfinite(y)):
        return SubdopplerFeature("unresolved", "The sampled scan contains non-finite values.")
    dx = np.diff(x)
    if not np.all(dx > 0.0):
        return SubdopplerFeature("unresolved", "The frequency axis must be strictly increasing.")

    n = y.size
    win = max(5, (n // 60) | 1)
    win = min(win, n if n % 2 else n - 1)
    pad = win // 2
    ypad = np.pad(y, pad, mode="edge")
    windows = np.lib.stride_tricks.sliding_window_view(ypad, win)
    smooth = np.median(windows, axis=1)
    residual = np.abs(y - smooth)

    candidates = np.ones(n, dtype=bool)
    if search_window is not None:
        lo, hi = map(float, search_window)
        if not np.isfinite(lo) or not np.isfinite(hi) or lo >= hi:
            return SubdopplerFeature("unresolved", "The feature search window is invalid.")
        candidates &= (x >= lo) & (x <= hi)
    if not np.any(candidates):
        return SubdopplerFeature("unresolved", "No samples fall inside the feature search window.")

    indices = np.flatnonzero(candidates)
    ic = int(indices[np.argmax(residual[indices])])
    amplitude = float(residual[ic])
    center = float(x[ic])
    if amplitude <= float(min_amplitude):
        return SubdopplerFeature(
            "unresolved", "No finite sub-Doppler residual exceeds the detection floor.",
            center=center, amplitude=amplitude)

    half = 0.5 * amplitude
    left_below = np.flatnonzero(residual[:ic + 1] <= half)
    right_below = np.flatnonzero(residual[ic:] <= half)
    if left_below.size == 0 or right_below.size == 0:
        return SubdopplerFeature(
            "scan-edge-limited",
            "A half-height crossing is clipped by the displayed scan edge.",
            center=center, scan_edge_distance=0.0, amplitude=amplitude)

    il = int(left_below[-1])
    ir = int(ic + right_below[0])
    if il >= ic or ir <= ic:
        return SubdopplerFeature(
            "unresolved", "The sampled residual has no finite half-height interval.",
            center=center, amplitude=amplitude)

    def _interpolated_crossing(i0, i1):
        y0, y1 = residual[i0], residual[i1]
        if y1 == y0:
            return float(0.5 * (x[i0] + x[i1]))
        fraction = (half - y0) / (y1 - y0)
        return float(x[i0] + fraction * (x[i1] - x[i0]))

    left = _interpolated_crossing(il, il + 1)
    right = _interpolated_crossing(ir - 1, ir)
    width = float(right - left)
    if not np.isfinite(width) or width <= 0.0:
        return SubdopplerFeature(
            "unresolved", "The interpolated half-height width is not positive.",
            center=center, amplitude=amplitude)

    local_steps = dx[(x[:-1] < right) & (x[1:] > left)]
    sample_step = float(np.median(local_steps)) if local_steps.size else float("nan")
    samples_per_fwhm = width / sample_step if sample_step > 0.0 else float("nan")
    edge_distance = float(min(left - x[0], x[-1] - right))
    edge_widths = edge_distance / width

    if edge_widths < MIN_SUBDOPPLER_EDGE_CLEARANCE_FWHM:
        status = "scan-edge-limited"
        reason = "Less than one feature FWHM remains between a half-height edge and the scan edge."
    elif samples_per_fwhm < MIN_SUBDOPPLER_SAMPLES_PER_FWHM:
        status = "resolution-limited"
        reason = "Fewer than six sample intervals span the interpolated feature FWHM."
    else:
        status = "resolved"
        reason = "Half-height edges are clear of the scan boundary and sampled by at least six intervals."

    return SubdopplerFeature(
        status=status,
        reason=reason,
        center=center,
        left_half_height=left,
        right_half_height=right,
        fwhm=width,
        samples_per_fwhm=float(samples_per_fwhm),
        scan_edge_distance=max(edge_distance, 0.0),
        amplitude=amplitude,
    )


def narrowest_subdoppler(x_ghz, T_trans):
    """Backward-compatible ``(width, location)`` sub-Doppler readout.

    New code should use :func:`subdoppler_feature` and check its resolution
    status before presenting the numeric width or a derivative-based lock slope.
    """
    feature = subdoppler_feature(x_ghz, T_trans)
    if not feature.detected:
        return float("nan"), float("nan")
    return feature.fwhm, feature.center


def lorentz_fwhm(B, y):
    """FWHM (in B units) of a zero-field **Lorentzian** feature.

    Model-dependent: assumes a Lorentzian-in-B shape and fits 1/signal linear in
    B² (signal = y − wing baseline). Only valid for a Lorentzian-like central
    feature; do not use as a generic width estimator (use `window_fwhm` /
    `halfwidth_from_center` for that). `nan` if the fit is ill-defined.
    """
    base = 0.5 * (y[0] + y[-1])
    s = y - base
    ic = int(np.argmin(np.abs(B)))
    s0 = s[ic]
    if abs(s0) < 1e-30:
        return float("nan")
    core = (np.sign(s) == np.sign(s0)) & (np.abs(s) >= 0.2 * abs(s0))
    if int(core.sum()) < 3:
        return float("nan")
    a, b = np.polyfit(B[core] ** 2, 1.0 / s[core], 1)
    if a == 0 or b / a <= 0:
        return float("nan")
    return 2.0 * math.sqrt(b / a)
