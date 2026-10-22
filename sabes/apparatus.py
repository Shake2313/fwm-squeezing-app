"""Dated apparatus records for the real SABES table.

This inventory is separate from the optical transfer model and calibration.
Equipment identity, a recorded operating point and an unknown port are useful
even when none supplies a predictive transfer function. Layout node IDs link
the records to the drawing; no display coordinates are treated as measurements.
"""
from copy import deepcopy
from dataclasses import asdict, dataclass, field
from datetime import date
from pathlib import Path
from typing import Any

import json


SCHEMA_VERSION = "sabes_apparatus_v1"
DEFAULT_APPARATUS_PATH = Path(__file__).with_name("apparatus_20261016.json")
RECORD_STATUSES = frozenset({
    "confirmed", "recorded", "tentative", "deferred", "unknown",
    "expectation", "schematic", "historical",
})
SOURCE_ORIGINS = frozenset({
    "human_direct", "lab_record", "schematic", "research",
})


@dataclass(frozen=True)
class SourceRecord:
    id: str
    label: str
    origin: str
    reference: str
    date: str | None = None
    note: str = ""


@dataclass(frozen=True)
class EquipmentRecord:
    id: str
    label: str
    manufacturer: str | None
    model: str | None
    status: str
    source_ids: tuple[str, ...]
    node_ids: tuple[str, ...] = ()
    notes: str = ""
    attributes: dict[str, Any] = field(default_factory=dict)

    @property
    def model_label(self):
        return " ".join(part for part in (self.manufacturer, self.model) if part)


@dataclass(frozen=True)
class ConnectionRecord:
    id: str
    label: str
    category: str
    source_id: str
    target_id: str
    status: str
    source_ids: tuple[str, ...]
    source_port: str | None = None
    target_port: str | None = None
    via: tuple[str, ...] = ()
    notes: str = ""


@dataclass(frozen=True)
class ApparatusFact:
    id: str
    category: str
    label: str
    value: Any
    status: str
    source_ids: tuple[str, ...]
    unit: str = ""
    scope: str = "FWM"
    as_of: str | None = None
    notes: str = ""


@dataclass(frozen=True)
class DeferredQuestion:
    id: str
    equipment_id: str | None
    question: str
    contexts: tuple[str, ...]
    source_ids: tuple[str, ...]
    notes: str = ""


def _records(items, record_type, tuple_fields):
    records = []
    for item in items:
        if not isinstance(item, dict):
            raise ValueError(f"{record_type.__name__} must be an object")
        data = deepcopy(item)
        for key in tuple_fields:
            if key in data:
                if not isinstance(data[key], (list, tuple)):
                    raise ValueError(f"{record_type.__name__}.{key} must be a list")
                data[key] = tuple(data[key])
        try:
            records.append(record_type(**data))
        except TypeError as exc:
            raise ValueError(f"invalid {record_type.__name__}: {exc}") from exc
    return tuple(records)


def _iso_date(value, label):
    try:
        date.fromisoformat(value)
    except (TypeError, ValueError):
        raise ValueError(f"{label}: expected an ISO date, got {value!r}") from None


class ApparatusSnapshot:
    """Validated, JSON-backed inventory, independent of runtime physics.

    Attributes are immutable by convention, like ``Calibration``. ``as_dict``
    returns an independent copy for export or making a revised dated snapshot.
    Source/status fields survive every round trip, including explicit nulls.
    """

    def __init__(self, *, schema_version, title, baseline_date, confirmed_on,
                 update_policy, sources, equipment, connections, facts,
                 deferred=()):
        self.schema_version = schema_version
        self.title = title
        self.baseline_date = baseline_date
        self.confirmed_on = confirmed_on
        self.update_policy = update_policy
        self.sources = tuple(sources)
        self.equipment = tuple(equipment)
        self.connections = tuple(connections)
        self.facts = tuple(facts)
        self.deferred = tuple(deferred)
        self._validate()

    def _validate(self):
        if self.schema_version != SCHEMA_VERSION:
            raise ValueError(f"unsupported apparatus schema {self.schema_version!r}")
        _iso_date(self.baseline_date, "baseline_date")
        _iso_date(self.confirmed_on, "confirmed_on")
        if self.confirmed_on < self.baseline_date:
            raise ValueError("confirmed_on precedes baseline_date")
        if not self.title or not self.update_policy:
            raise ValueError("apparatus title and update_policy must be nonempty")
        if not self.sources or not self.equipment:
            raise ValueError("apparatus requires sources and equipment")
        for group in (self.sources, self.equipment, self.connections,
                      self.facts, self.deferred):
            ids = [item.id for item in group]
            if any(not isinstance(key, str) or not key.strip() for key in ids):
                raise ValueError("record IDs must be nonempty strings")
            if len(ids) != len(set(ids)):
                raise ValueError("duplicate record IDs")
        sources = {item.id for item in self.sources}
        equipment = {item.id for item in self.equipment}
        for item in self.sources:
            if item.origin not in SOURCE_ORIGINS:
                raise ValueError(f"{item.id}: unknown source origin {item.origin!r}")
            if not item.reference:
                raise ValueError(f"{item.id}: missing source reference")
            if item.date is not None:
                _iso_date(item.date, item.id)
        for item in (*self.equipment, *self.connections, *self.facts, *self.deferred):
            if (not item.source_ids
                    or any(not isinstance(key, str) for key in item.source_ids)
                    or set(item.source_ids) - sources):
                raise ValueError(f"{item.id}: missing or unknown source IDs")
            if hasattr(item, "status") and item.status not in RECORD_STATUSES:
                raise ValueError(f"{item.id}: unknown status {item.status!r}")
        for item in self.equipment:
            for name in ("manufacturer", "model"):
                value = getattr(item, name)
                if value is not None and (not isinstance(value, str) or not value.strip()):
                    raise ValueError(f"{item.id}: {name} must be a nonempty string or null")
            if len(item.node_ids) != len(set(item.node_ids)):
                raise ValueError(f"{item.id}: duplicate layout node IDs")
            if any(not isinstance(key, str) or not key for key in item.node_ids):
                raise ValueError(f"{item.id}: invalid layout node ID")
            json.dumps(item.attributes, allow_nan=False)
        for item in self.connections:
            if item.source_id not in equipment or item.target_id not in equipment:
                raise ValueError(f"{item.id}: unknown equipment endpoint")
            if item.source_id == item.target_id:
                raise ValueError(f"{item.id}: connection endpoints must differ")
            for port in (item.source_port, item.target_port):
                if port is not None and (not isinstance(port, str) or not port.strip()):
                    raise ValueError(f"{item.id}: ports must be strings or null")
        for item in self.facts:
            if item.as_of is not None:
                _iso_date(item.as_of, item.id)
            if item.status in {"deferred", "unknown"} and item.value is not None:
                raise ValueError(f"{item.id}: unknown/deferred fact must retain null")
            json.dumps(item.value, allow_nan=False)
        for item in self.deferred:
            if item.equipment_id is not None and item.equipment_id not in equipment:
                raise ValueError(f"{item.id}: unknown deferred equipment ID")
            if not item.contexts:
                raise ValueError(f"{item.id}: deferred question requires a context")

    def get_equipment(self, equipment_id):
        for item in self.equipment:
            if item.id == equipment_id:
                return item
        raise KeyError(f"no apparatus equipment {equipment_id!r}")

    def get_fact(self, fact_id):
        for item in self.facts:
            if item.id == fact_id:
                return item
        raise KeyError(f"no apparatus fact {fact_id!r}")

    def equipment_for_node(self, node_id):
        return tuple(item for item in self.equipment if node_id in item.node_ids)

    def validate_layouts(self, layouts):
        """Reject stale optic mappings; never derive geometry from the drawing."""
        known = {node.id for layout in layouts for node in layout.nodes}
        missing = {key for item in self.equipment for key in item.node_ids} - known
        if missing:
            raise ValueError(f"apparatus references unknown layout nodes {sorted(missing)}")

    def deferred_questions(self, context=None):
        """Return pending questions; callers decide when relevant work starts.

        Exact context keys avoid creating a time-based reminder or making an
        unknown model block unrelated apparatus work.
        """
        return tuple(item for item in self.deferred
                     if context is None or context in item.contexts)

    def equipment_rows(self):
        return tuple({
            "id": item.id, "label": item.label,
            "manufacturer": item.manufacturer, "model": item.model,
            "status": item.status, "node_ids": ", ".join(item.node_ids),
            "source_ids": ", ".join(item.source_ids), "notes": item.notes,
        } for item in self.equipment)

    def connection_rows(self, category=None):
        def endpoint(key, port):
            label = self.get_equipment(key).label
            return f"{label} · {port}" if port else label
        return tuple({
            "id": item.id, "label": item.label, "category": item.category,
            "from": endpoint(item.source_id, item.source_port),
            "to": endpoint(item.target_id, item.target_port),
            "via": " → ".join(item.via), "status": item.status,
            "source_ids": ", ".join(item.source_ids), "notes": item.notes,
        } for item in self.connections if category is None or item.category == category)

    def fact_rows(self, category=None):
        def display(value):
            if value is None:
                return "미확인"
            if isinstance(value, (dict, list, tuple, bool)):
                return json.dumps(value, ensure_ascii=False)
            return str(value)
        return tuple({
            "id": item.id, "category": item.category, "label": item.label,
            "value": display(item.value), "unit": item.unit,
            "status": item.status, "scope": item.scope, "as_of": item.as_of,
            "source_ids": ", ".join(item.source_ids), "notes": item.notes,
        } for item in self.facts if category is None or item.category == category)

    def as_dict(self):
        payload = {key: getattr(self, key) for key in (
            "schema_version", "title", "baseline_date", "confirmed_on", "update_policy")}
        for name in ("sources", "equipment", "connections", "facts", "deferred"):
            payload[name] = [asdict(item) for item in getattr(self, name)]
        # Use JSON's list/null conventions in exported Python dictionaries too.
        return json.loads(json.dumps(payload, ensure_ascii=False, allow_nan=False))

    def to_json(self, **kwargs):
        return json.dumps(self.as_dict(), ensure_ascii=False, allow_nan=False,
                          **kwargs)

    @classmethod
    def from_dict(cls, payload):
        data = deepcopy(payload)
        if not isinstance(data, dict):
            raise ValueError("apparatus snapshot must be a JSON object")
        types = {
            "sources": (SourceRecord, ()),
            "equipment": (EquipmentRecord, ("source_ids", "node_ids")),
            "connections": (ConnectionRecord, ("source_ids", "via")),
            "facts": (ApparatusFact, ("source_ids",)),
            "deferred": (DeferredQuestion, ("source_ids", "contexts")),
        }
        for name, (record_type, tuple_fields) in types.items():
            values = data.get(name, [])
            if not isinstance(values, list):
                raise ValueError(f"{name} must be a list")
            data[name] = _records(values, record_type, tuple_fields)
        try:
            return cls(**data)
        except TypeError as exc:
            raise ValueError(f"invalid apparatus snapshot: {exc}") from exc

    @classmethod
    def load(cls, path):
        with Path(path).open(encoding="utf-8-sig") as handle:
            return cls.from_dict(json.load(handle))


def default_apparatus():
    """Fresh copy of the user-confirmed 2026-10-16 apparatus snapshot."""
    return ApparatusSnapshot.load(DEFAULT_APPARATUS_PATH)
