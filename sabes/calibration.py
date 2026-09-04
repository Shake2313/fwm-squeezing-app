"""
Coefficients, and where each one came from.

SABES predicts lab quantities from lab settings, so it is only as good as the
numbers connecting them. The single most common way a digital twin goes wrong is
that a plausible-looking guess and a measured value are stored the same way and
become indistinguishable a month later. Every coefficient here therefore carries
a `provenance` tag, and the UI is expected to show it:

    fitted     solved from measured data
    datasheet  manufacturer specification for the exact part
    paper      stated in Sim et al., Sci. Rep. 15, 7727 (2025)
    lab        reported by the experimenter, not a formal measurement
    nominal    a placeholder not backed by a measurement
    hand       adjusted interactively during a session

`fitted`, `datasheet`, and `paper` are treated as verified. `lab`, `nominal`,
and `hand` are returned by `unverified()` so the UI does not present them as
calibrated measurements.

Coefficients live in JSON so they can be replaced without touching code, which is
the whole point: when a real etalon scan or P-I curve arrives, only the file
changes.
"""
from dataclasses import dataclass, replace
from pathlib import Path
from typing import Any, Dict

import json

PROVENANCE_ORDER = ("fitted", "datasheet", "paper", "lab", "nominal", "hand")
VERIFIED_PROVENANCE = frozenset({"fitted", "datasheet", "paper"})

DEFAULT_CALIBRATION_PATH = Path(__file__).with_name("calibration_sim2025.json")


@dataclass(frozen=True)
class Coefficient:
    name: str
    value: float
    unit: str = ""
    provenance: str = "nominal"
    note: str = ""

    def __post_init__(self):
        if self.provenance not in PROVENANCE_ORDER:
            raise ValueError(
                f"{self.name}: unknown provenance {self.provenance!r}; "
                f"expected one of {PROVENANCE_ORDER}")

    @property
    def verified(self):
        return self.provenance in VERIFIED_PROVENANCE

    def as_dict(self):
        return {"value": self.value, "unit": self.unit,
                "provenance": self.provenance, "note": self.note}


class Calibration:
    """An immutable-by-convention bag of `Coefficient`s keyed by name."""

    def __init__(self, coefficients):
        self._items: Dict[str, Coefficient] = dict(coefficients)

    # ---- access ----
    def __contains__(self, name):
        return name in self._items

    def __len__(self):
        return len(self._items)

    def __iter__(self):
        return iter(self._items.values())

    def get(self, name) -> Coefficient:
        try:
            return self._items[name]
        except KeyError:
            raise KeyError(
                f"no calibration coefficient named {name!r}. "
                f"Add it to the calibration file rather than hard-coding it."
            ) from None

    def value(self, name, default=None):
        if default is not None and name not in self._items:
            return default
        return self.get(name).value

    def names(self):
        return tuple(sorted(self._items))

    # ---- mutation returns a new object, so a session can branch safely ----
    def with_value(self, name, value, provenance="hand", note=None):
        """Override one coefficient. Defaults to tagging the result `hand`."""
        current = self.get(name)
        updated = replace(
            current, value=float(value), provenance=provenance,
            note=current.note if note is None else note)
        items = dict(self._items)
        items[name] = updated
        return Calibration(items)

    def with_values(self, updates, provenance="hand"):
        """Override changed coefficients with one copy of the calibration map."""
        items = None
        for name, value in updates.items():
            current = self.get(name)
            value = float(value)
            if value == current.value:
                continue
            if items is None:
                items = dict(self._items)
            items[name] = replace(current, value=value, provenance=provenance)
        return self if items is None else Calibration(items)

    def anchored_to(self, name, measured_value, note=""):
        """Promote a coefficient to `fitted` against a measurement.

        This is the intended way to turn interactive tweaking into a recorded
        calibration: adjust until the model matches, then anchor.
        """
        return self.with_value(name, measured_value, provenance="fitted",
                               note=note or "anchored to a measured value")

    # ---- reporting ----
    def unverified(self):
        """Unverified coefficients in declared provenance order."""
        return tuple(sorted(
            (c for c in self._items.values() if not c.verified),
            key=lambda c: (PROVENANCE_ORDER.index(c.provenance), c.name)))

    def table_rows(self):
        return tuple(
            (c.name, c.value, c.unit, c.provenance, c.note)
            for c in sorted(self._items.values(),
                            key=lambda c: (PROVENANCE_ORDER.index(c.provenance),
                                           c.name)))

    # ---- persistence ----
    def to_dict(self) -> Dict[str, Any]:
        return {name: self._items[name].as_dict() for name in sorted(self._items)}

    def save(self, path):
        Path(path).write_text(
            json.dumps(self.to_dict(), indent=2, ensure_ascii=False) + "\n",
            encoding="utf-8")

    @classmethod
    def from_dict(cls, payload):
        return cls({
            name: Coefficient(
                name=name,
                value=float(entry["value"]),
                unit=entry.get("unit", ""),
                provenance=entry.get("provenance", "nominal"),
                note=entry.get("note", ""))
            for name, entry in payload.items()
        })

    @classmethod
    def load(cls, path=None):
        path = Path(path) if path is not None else DEFAULT_CALIBRATION_PATH
        return cls.from_dict(json.loads(path.read_text(encoding="utf-8")))


def default_calibration():
    """The Sim et al. setup as currently known. Reloaded from disk each call."""
    return Calibration.load()
