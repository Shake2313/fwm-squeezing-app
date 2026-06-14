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

import numpy as np


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
    pk = np.nanmax(y)
    if not np.isfinite(pk) or pk <= 0:
        return float("nan")
    above = x[y >= 0.5 * pk]
    return float(above.max() - above.min()) if above.size > 1 else float("nan")


def fwhm_interp(x, y):
    """FWHM of the tallest peak with linear edge interpolation; `nan` if
    ill-defined."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if x.size < 2 or not np.any(np.isfinite(y)):
        return np.nan
    peak = float(np.nanmax(y))
    if peak <= 0:
        return np.nan
    half = 0.5 * peak
    above = np.flatnonzero(y >= half)
    if above.size == 0:
        return np.nan
    lo = int(above[0])
    hi = int(above[-1])

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


def narrowest_subdoppler(x_ghz, T_trans):
    """Width [same units as x] and location of the sharpest Doppler-free feature.

    Subtracts a running-median background and finds the tallest residual, then
    measures its half-width. Shape-agnostic. Returns (nan, nan) if none.
    """
    y = np.asarray(T_trans)
    n = y.size
    win = max(5, (n // 60) | 1)
    pad = win // 2
    ypad = np.pad(y, pad, mode="edge")
    smooth = np.array([np.median(ypad[i:i + win]) for i in range(n)])
    resid = np.abs(y - smooth)
    if resid.max() <= 1e-6:
        return float("nan"), float("nan")
    ic = int(np.argmax(resid))
    half = 0.5 * resid[ic]
    i = ic
    while i > 0 and resid[i] >= half:
        i -= 1
    j = ic
    while j < n - 1 and resid[j] >= half:
        j += 1
    return float(x_ghz[j] - x_ghz[i]), float(x_ghz[ic])


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
