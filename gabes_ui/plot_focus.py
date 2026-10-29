"""Display-only resonance framing of existing samples; never resample a curve."""
from contextlib import contextmanager
from dataclasses import dataclass
import hashlib

import numpy as np

from gabes.lineshape import subdoppler_feature, window_fwhm


@dataclass(frozen=True)
class FocusWindow:
    limits: tuple[float, float]
    description: str


def _sampled_lines(axis):
    """Exclude axes-relative reference lines and non-numerical artists."""
    for line in axis.lines:
        if line.get_transform() != axis.transData:
            continue
        try:
            x = np.asarray(line.get_xdata(), dtype=float)
            y = np.asarray(line.get_ydata(), dtype=float)
        except (TypeError, ValueError):
            continue
        if x.ndim == 1 and x.size >= 5 and y.shape == x.shape:
            yield line, x, y


def focus_signature(fig):
    """Keep a focused window stable when only the selected point changes."""
    curve = next(_sampled_lines(fig.axes[0]), None) if fig.axes else None
    if curve is None:
        return None
    digest = hashlib.blake2b(digest_size=16)
    for values in curve[1:]:
        digest.update(np.ascontiguousarray(values).tobytes())
    return digest.hexdigest()


def resonance_window(fig, scheme_name, raw=None, view=None):
    """Choose a viewing window in the plotted units, without changing readouts.

    FWM focuses the selected operating point. Lambda uses its declared two-photon
    centre, and SAS uses the same sampled-feature diagnostic as its readout.
    Other graphs, including heatmaps and independent-axis biphoton panels, keep
    their original view until a meaningful focus contract is available.
    """
    raw, view = raw or {}, view or {}
    if scheme_name not in {"sas", "lambda", "fwm"} or not fig.axes:
        return None
    if scheme_name == "fwm" and raw.get("kind") == "biphoton":
        return None
    curve = next(_sampled_lines(fig.axes[0]), None)
    if curve is None:
        return None
    _, x, y = curve
    if (not np.all(np.isfinite(x)) or not np.all(np.diff(x) > 0)
            or not np.any(np.isfinite(y))):
        return None
    lo, hi = sorted(fig.axes[0].get_xlim())
    lo, hi = max(lo, float(x[0])), min(hi, float(x[-1]))
    if hi <= lo:
        return None
    labels = {metric.get("label") for metric in view.get("metrics", [])}
    if scheme_name == "fwm":
        if len(fig.axes) != 2 or len(fig.axes[0].lines) < 2:
            return None
        center = float(fig.axes[0].lines[1].get_xdata()[0])
        half = (hi - lo) / 16
        description = "Selected operating point"
    elif scheme_name == "lambda":
        center = float(raw.get("Dc", 0.0)) / (2 * np.pi * 1e6)
        ic = int(np.argmin(np.abs(x - center)))
        # A gap at the resonance is not a zero or a focusable peak.
        if not np.isfinite(y[ic]):
            return None
        if "AT splitting" in labels:
            half = 1.5 * float(raw.get("coupling_rabi_mhz", 0.0))
            description = "Coupling-centred region"
        else:
            # Bound the analysis to the contiguous finite segment at centre.
            left, right = ic, ic
            while left > 0 and np.isfinite(y[left - 1]):
                left -= 1
            while right + 1 < y.size and np.isfinite(y[right + 1]):
                right += 1
            width = window_fwhm(x[left:right + 1], y[left:right + 1], ic - left)
            if not np.isfinite(width) or width <= 0 or ic in (left, right):
                return None
            half = 3 * width
            description = "Two-photon resonance"
    else:
        if "SAS resolution" in labels:
            search = None
            if raw.get("mode") == "generic":
                offsets = np.asarray(raw.get("offsets", []), dtype=float) / (2 * np.pi * 1e6)
                gamma = float(raw.get("gamma_eff", 0.0)) / (2 * np.pi * 1e6)
                if offsets.size and gamma > 0:
                    search = (offsets[0] - 10 * gamma, offsets[0] + 10 * gamma)
            feature = subdoppler_feature(x, y, search_window=search)
            if not feature.detected:
                return None
            center, half = feature.center, 3 * feature.fwhm
            description = "Sampled sub-Doppler feature"
        else:
            if not np.all(np.isfinite(y)):
                return None
            ic = int(np.argmin(y))
            center, width = x[ic], window_fwhm(x, -y, ic)
            if not np.isfinite(width) or width <= 0:
                return None
            half = 1.5 * width
            description = "Absorption resonance"
    if not np.isfinite([center, half]).all() or half <= 0 or not lo <= center <= hi:
        return None
    ic = int(np.argmin(np.abs(x - center)))
    spacing = np.median(np.diff(x[max(0, ic - 4):min(x.size, ic + 5)]))
    # Keep several real samples visible even for a resolution-limited feature.
    half = max(float(half), 4 * float(spacing))
    limits = (max(lo, center - half), min(hi, center + half))
    if limits[1] - limits[0] >= 0.95 * (hi - lo):
        return None
    return FocusWindow(tuple(map(float, limits)), description)


@contextmanager
def focused_axes(fig, window):
    """Temporarily frame axes, restoring cached figures even if render fails."""
    if window is None:
        yield fig
        return
    axis_states = [(ax, ax.get_xlim(), ax.get_ylim(), ax.get_autoscalex_on(),
                    ax.get_autoscaley_on()) for ax in fig.axes]
    line_states = []
    try:
        lo, hi = window.limits
        for ax, xlim, ylim, _, _ in axis_states:
            ax.set_xlim((hi, lo) if xlim[0] > xlim[1] else (lo, hi))
            values = []
            for line, x, y in _sampled_lines(ax):
                inside = (x >= lo) & (x <= hi)
                finite = inside & np.isfinite(y)
                if ax.get_yscale() == "log":
                    finite &= y > 0
                values.extend(y[finite])
                if 0 < np.count_nonzero(finite) <= 40:
                    line_states.append((line, line.get_marker(), line.get_markersize()))
                    line.set_marker(".")
                    line.set_markersize(4)
            if values:
                low, high = float(np.min(values)), float(np.max(values))
                if ax.get_yscale() == "log":
                    low, high = np.log10([low, high])
                    pad = max((high - low) * .08, .04)
                    bounds = (10 ** (low - pad), 10 ** (high + pad))
                else:
                    pad = (max((high - low) * .08, abs(low) * .01)
                           if low != 0 or high != 0 else .05)
                    bounds = (low - pad, high + pad)
                ax.set_ylim(bounds[::-1] if ylim[0] > ylim[1] else bounds)
        yield fig
    finally:
        for line, marker, size in line_states:
            line.set_marker(marker)
            line.set_markersize(size)
        for ax, xlim, ylim, auto_x, auto_y in axis_states:
            ax.set_xlim(xlim)
            ax.set_ylim(ylim)
        # Setting a shared limit changes its siblings' autoscale flags too.
        # Restore flags only after every shared group has its limits restored.
        for ax, _, _, auto_x, auto_y in axis_states:
            ax.set_autoscalex_on(auto_x)
            ax.set_autoscaley_on(auto_y)
