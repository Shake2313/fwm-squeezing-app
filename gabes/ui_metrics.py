"""Pure helpers for arranging scheme metrics in the UI."""

from __future__ import annotations

import re
from collections.abc import Iterable, Mapping
from operator import index
from typing import TypeVar


MetricT = TypeVar("MetricT", bound=Mapping[str, object])


_NUMBER_AND_UNIT = re.compile(
    r"^\s*(?P<number>[+\-−]?(?:\d+(?:\.\d*)?|\.\d+)"
    r"(?:[eE][+\-]?\d+)?)\s+(?P<unit>\S+)\s*$"
)


def split_metric_value(
    value: object, *, kind: object = None
) -> tuple[str, str | None]:
    """Separate a displayed number from a compact unit token when unambiguous.

    Scheme metrics are intentionally presentation-ready strings rather than
    numeric records.  Explicit ``kind="status"`` values always remain intact;
    other values split only on the conservative ``number + one unit token``
    shape.  Both returned parts are still plain text and must be escaped by an
    HTML renderer.
    """
    text = str(value).strip()
    if str(kind).lower() == "status":
        return text, None
    match = _NUMBER_AND_UNIT.fullmatch(text)
    if match is None:
        return text, None
    return match.group("number"), match.group("unit")


def partition_metrics(
    metrics: Iterable[MetricT], hero_count: int = 2
) -> tuple[list[MetricT], list[MetricT]]:
    """Split metrics into prominent heroes and an ordered compact ribbon.

    Metrics explicitly tagged with ``tier == "hero"`` take priority.  Any
    unfilled hero slots are then populated from the remaining metrics in their
    original order.  The returned lists contain the original metric objects;
    neither the input iterable nor its metric mappings are modified.
    """
    items = list(metrics)
    limit = max(0, index(hero_count))
    selected = [False] * len(items)
    heroes: list[MetricT] = []

    for position, metric in enumerate(items):
        if len(heroes) == limit:
            break
        if metric.get("tier") == "hero":
            heroes.append(metric)
            selected[position] = True

    if len(heroes) < limit:
        for position, metric in enumerate(items):
            if selected[position]:
                continue
            heroes.append(metric)
            selected[position] = True
            if len(heroes) == limit:
                break

    ribbon = [metric for position, metric in enumerate(items) if not selected[position]]
    return heroes, ribbon
