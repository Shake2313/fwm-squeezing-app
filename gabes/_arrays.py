"""Small array utilities shared across the engine."""

import numpy as np


def readonly_float_array(values, *, copy=True):
    """Return ``values`` as a read-only float64 array.

    ``copy=True`` freezes a new array and leaves the input untouched.
    ``copy=False`` converts with ``np.asarray``: a float64 ndarray input is the
    returned object and becomes read-only in place.
    """
    result = (np.array(values, dtype=float, copy=True) if copy
              else np.asarray(values, dtype=float))
    result.setflags(write=False)
    return result
