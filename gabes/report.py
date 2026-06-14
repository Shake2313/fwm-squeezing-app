"""
Shared readout-table formatting for the schemes.

`derived_table` builds the two-column "| Quantity | Value |" markdown block that
every scheme returns in its `observables(...)["tables"]`. Centralising it keeps
the formatting — and the unicode typography (°C, µT, µW, Γ, Ω, 2π) — uniform
across schemes instead of being hand-spelled per scheme.

Note: this is for readout/markdown only. Matplotlib axis/title strings stay ASCII
(the plot lock guards against mathtext layout crashes), so do not route those
through here.
"""
from collections.abc import Mapping


def _iter_rows(rows):
    if isinstance(rows, Mapping):
        return rows.items()
    return rows


def derived_table(rows, title="Derived quantities"):
    """A `{"title", "markdown"}` two-column table from (label, value) rows.

    `rows` is a mapping or an iterable of (label, value) pairs. Values are
    stringified as-is — format numbers before passing them in. Rows whose value
    is None are skipped, so callers can build rows conditionally.
    """
    body = "".join(
        f"| {label} | {value} |\n"
        for label, value in _iter_rows(rows)
        if value is not None
    )
    return {"title": title, "markdown": "| Quantity | Value |\n|---|---|\n" + body}
