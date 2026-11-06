"""Shared contract for small on-demand sweeps of ordinary scheme conditions."""

import numpy as np
from contextlib import contextmanager
from contextvars import ContextVar


SCHEMA = "gabes-model-comparison-v2"

_condition_values = ContextVar("comparison_condition_values", default=None)


@contextmanager
def comparison_values(key, values):
    """Request explicit conditions from the ordinary comparison adapters.

    Context-local so an incremental sweep cannot change another browser's or
    the main thread's comparison. Atomic solves and readouts stay in schemes.
    """
    token = _condition_values.set((key, tuple(values)))
    try:
        yield
    finally:
        _condition_values.reset(token)


def condition_views(scheme, params, key, values, *, resolve_raw,
                    current_raw=None, current_view=None):
    """Yield each distinct condition in order, reusing the committed solution.

    resolve_raw accepts full parameters; the caller's existing compute cache
    selects recompute keys. A view is only reused at the exact current value.
    No Matplotlib figures are constructed for the other conditions.
    """
    requested = _condition_values.get()
    if requested is not None and requested[0] == key:
        values = requested[1]
    for value in dict.fromkeys(values):
        condition = {**params, key: value}
        if value == params[key] and current_raw is not None:
            raw, view = current_raw, current_view
        else:
            raw, view = resolve_raw(condition), None
        if view is None:
            view = scheme.headless_observables(raw, condition)
        yield condition, raw, view


def sample_values(values):
    """Null denotes an unavailable observation, including a gap in a difference."""
    array = np.asarray(values, dtype=float)
    if np.isfinite(array).all():
        return array.tolist()
    return [float(v) if np.isfinite(v) else None for v in array]


def same_grid(expected, sample_x):
    """Do not silently interpolate a baseline onto a different numerical grid."""
    if not np.array_equal(expected, sample_x):
        raise ValueError("Model comparison requires the same sampled axis for every condition.")
