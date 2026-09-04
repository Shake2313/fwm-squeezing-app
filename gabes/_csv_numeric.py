"""Numeric helpers shared by the CSV loaders."""

import numpy as np


def sort_and_median_duplicates(x, y):
    """Sort by x and replace duplicate-y groups with their median."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if x.shape != y.shape or x.ndim != 1:
        raise ValueError("x and y must be equal one-dimensional arrays")
    if x.size == 0:
        return x.copy(), y.copy()

    order = np.lexsort((y, x))
    sorted_x = x[order]
    sorted_y = y[order]
    unique_x, starts, counts = np.unique(
        sorted_x, return_index=True, return_counts=True
    )
    middle = starts + counts // 2
    medians = sorted_y[middle].copy()
    even = counts % 2 == 0
    medians[even] = 0.5 * (
        sorted_y[middle[even] - 1] + sorted_y[middle[even]]
    )
    return unique_x, medians
