"""Helpers shared by the CSV loaders: text decoding and duplicate merging."""

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


def decode_csv_text(data, error=ValueError):
    """Decode CSV bytes; return ``(text, encoding)``.

    BOM-marked UTF-16 first (malformed text raises ``error``, the loader's
    exception type), then UTF-8 with or without a BOM, then CP949. Numeric
    fields and delimiters are ASCII, so the final Latin-1 step, which maps every
    byte, lets an unsupported byte in metadata or ignored columns through
    instead of rejecting otherwise valid data. UTF-8 text without a BOM is
    reported as ``utf-8-sig``, the first codec that accepts it.
    """
    if data.startswith((b"\xff\xfe", b"\xfe\xff")):
        try:
            return data.decode("utf-16"), "utf-16"
        except UnicodeDecodeError as exc:
            raise error("Malformed UTF-16 CSV text.") from exc
    for encoding in ("utf-8-sig", "utf-8", "cp949"):
        try:
            return data.decode(encoding), encoding
        except UnicodeDecodeError:
            continue
    return data.decode("latin-1"), "latin-1"
