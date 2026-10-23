"""Portable named parameter snapshots, shared by file and browser storage.

Only input parameters declared by the active scheme are accepted. A missing
parameter uses its current schema default so an older snapshot can still be
loaded when the scheme gains an input. Derived readouts are never persisted.
"""
from __future__ import annotations

import json
import math
import re
from numbers import Real


FORMAT = "gabes-default"
VERSION = 1
MAX_PAYLOAD_BYTES = 256 * 1024
MAX_NAME_LENGTH = 80


class DefaultValidationError(ValueError):
    """A snapshot cannot be safely applied to the active scheme."""


def _name(value):
    if not isinstance(value, str):
        raise DefaultValidationError("Enter a name for this default.")
    value = value.strip()
    if not value:
        raise DefaultValidationError("Enter a name for this default.")
    if len(value) > MAX_NAME_LENGTH:
        raise DefaultValidationError(
            f"Default names must be {MAX_NAME_LENGTH} characters or fewer.")
    if any(ord(char) < 32 or ord(char) == 127 for char in value):
        raise DefaultValidationError("Default names cannot contain control characters.")
    return value


def _integer_input(spec):
    numbers = [v for v in (spec.default, spec.vmin, spec.vmax, spec.step)
               if v is not None]
    return bool(numbers) and all(isinstance(v, int) and not isinstance(v, bool)
                                 for v in numbers)


def _parameter_value(spec, value):
    label = spec.label or spec.name
    if isinstance(spec.default, bool) or spec.control == "checkbox":
        if not isinstance(value, bool):
            raise DefaultValidationError(f"{label} must be true or false.")
    elif isinstance(spec.default, Real):
        if isinstance(value, bool) or not isinstance(value, Real):
            raise DefaultValidationError(f"{label} must be a number.")
        try:
            finite = math.isfinite(value)
        except (OverflowError, ValueError):
            finite = False
        if not finite:
            raise DefaultValidationError(f"{label} must be a finite number.")
        if _integer_input(spec):
            if value != int(value):
                raise DefaultValidationError(f"{label} must be a whole number.")
            value = int(value)
        else:
            value = float(value)
        if spec.vmin is not None and value < spec.vmin:
            raise DefaultValidationError(f"{label} must be at least {spec.vmin}.")
        if spec.vmax is not None and value > spec.vmax:
            raise DefaultValidationError(f"{label} must be at most {spec.vmax}.")
    elif isinstance(spec.default, str):
        if not isinstance(value, str):
            raise DefaultValidationError(f"{label} must be text.")
    elif spec.default is None:
        if value is not None:
            raise DefaultValidationError(f"{label} must be null.")
    else:
        raise DefaultValidationError(f"{label} has an unsupported parameter type.")
    if spec.choices is not None and value not in spec.choices:
        raise DefaultValidationError(f"{label} has an unknown selection.")
    return value


def _params(scheme, values):
    if not isinstance(values, dict):
        raise DefaultValidationError("Default parameters must be a JSON object.")
    specs = [spec for spec in scheme.param_schema()
             if spec.control != "derived_slider"]
    allowed = {spec.name for spec in specs}
    unknown = set(values) - allowed
    if unknown:
        names = ", ".join(sorted(str(name) for name in unknown))
        raise DefaultValidationError(f"Unknown input parameters: {names}.")
    return {spec.name: _parameter_value(spec, values.get(spec.name, spec.default))
            for spec in specs}


def _encoded(record):
    try:
        encoded = json.dumps(record, ensure_ascii=False, allow_nan=False, indent=2)
        size = len(encoded.encode("utf-8"))
    except (TypeError, ValueError, OverflowError, UnicodeError, RecursionError) as exc:
        raise DefaultValidationError("Default must contain valid JSON values.") from exc
    if size > MAX_PAYLOAD_BYTES:
        raise DefaultValidationError("Default files must be 256 KiB or smaller.")
    return encoded


def create_default(scheme, name, params):
    """Capture named input values, validating them against the active schema."""
    record = {
        "format": FORMAT,
        "version": VERSION,
        "scheme": scheme.name,
        "name": _name(name),
        "params": _params(scheme, params),
    }
    _encoded(record)
    return record


def _unique_object(pairs):
    result = {}
    for key, value in pairs:
        if key in result:
            raise DefaultValidationError(f"Duplicate JSON field: {key}.")
        result[key] = value
    return result


def _invalid_constant(value):
    raise DefaultValidationError(f"Non-finite JSON number: {value}.")


def parse_default(scheme, payload):
    """Read UTF-8 JSON or a browser record and return validated input values."""
    if isinstance(payload, (bytes, bytearray)):
        if len(payload) > MAX_PAYLOAD_BYTES:
            raise DefaultValidationError("Default files must be 256 KiB or smaller.")
        try:
            payload = bytes(payload).decode("utf-8-sig")
        except UnicodeError as exc:
            raise DefaultValidationError("Default files must use UTF-8 JSON.") from exc
    if isinstance(payload, str):
        try:
            size = len(payload.encode("utf-8"))
        except UnicodeError as exc:
            raise DefaultValidationError("Default files must use UTF-8 JSON.") from exc
        if size > MAX_PAYLOAD_BYTES:
            raise DefaultValidationError("Default files must be 256 KiB or smaller.")
        try:
            payload = json.loads(payload, object_pairs_hook=_unique_object,
                                 parse_constant=_invalid_constant)
        except DefaultValidationError:
            raise
        except (ValueError, RecursionError) as exc:
            raise DefaultValidationError("This file is not valid default JSON.") from exc
    if not isinstance(payload, dict):
        raise DefaultValidationError("A default must be a JSON object.")
    _encoded(payload)
    expected = {"format", "version", "scheme", "name", "params"}
    if set(payload) != expected:
        raise DefaultValidationError("Default file fields are missing or unsupported.")
    if payload["format"] != FORMAT:
        raise DefaultValidationError("This file is not a GABES default.")
    if type(payload["version"]) is not int or payload["version"] != VERSION:
        raise DefaultValidationError("This default file version is not supported.")
    if payload["scheme"] != scheme.name:
        raise DefaultValidationError(
            f"This default belongs to a different scheme; select {payload['scheme']} first.")
    return create_default(scheme, payload["name"], payload["params"])


def dump_default(record):
    """Serialize a record for downloading or browser storage."""
    return _encoded(record) + "\n"


def default_filename(name):
    """Return a portable download filename, including Windows reserved names."""
    value = _name(name)
    value = re.sub(r'[<>:"/\\|?*\x00-\x1f\x7f]', "-", value).rstrip(" .")
    if value.lower().endswith(".json"):
        value = value[:-5].rstrip(" .")
    value = value or "default"
    if re.match(r"^(CON|PRN|AUX|NUL|COM[1-9]|LPT[1-9])(?:\.|$)", value, re.I):
        value = "default-" + value
    return value + ".json"
