"""Named defaults preserve input settings and reject incompatible snapshots."""
import copy

import pytest

from gabes import schemes
from gabes.schemes.base import ParamSpec
from gabes_ui.defaults_store import (
    MAX_PAYLOAD_BYTES,
    DefaultValidationError,
    create_default,
    default_filename,
    dump_default,
    parse_default,
)


class ExampleScheme:
    name = "example"

    def param_schema(self):
        return [
            ParamSpec("power", "Power", "Fields", 2.5, 0.0, 10.0, 0.1),
            ParamSpec("count", "Count", "Numerics", 3, 1, 9, 1),
            ParamSpec("enabled", "Enabled", "Mode", False, control="checkbox"),
            ParamSpec("mode", "Mode", "Mode", "a", choices=("a", "b")),
            ParamSpec("source", "Source", "Provenance", "unspecified", hidden=True),
            ParamSpec("result", "Result", "Readout", None, control="derived_slider"),
        ]


@pytest.fixture
def example():
    return ExampleScheme()


@pytest.mark.parametrize("scheme", schemes.all_schemes(), ids=lambda item: item.name)
def test_registered_defaults_round_trip(scheme):
    inputs = scheme.defaults()
    record = create_default(scheme, "실험 설정", inputs)
    assert parse_default(scheme, dump_default(record).encode("utf-8")) == record
    assert record["params"] == inputs
    assert not any(spec.name in record["params"] for spec in scheme.param_schema()
                   if spec.control == "derived_slider")


def test_snapshot_preserves_hidden_metadata_and_fills_new_parameters(example):
    record = create_default(example, "  Today  ", {"power": 4.25, "source": "Measured"})
    assert record["name"] == "Today"
    assert record["params"] == {
        "power": 4.25, "count": 3, "enabled": False, "mode": "a", "source": "Measured",
    }
    assert parse_default(example, record) == record
    assert parse_default(example, b"\xef\xbb\xbf" + dump_default(record).encode()) == record


@pytest.mark.parametrize("name", ["", " \t ", "a" * 81, "line\nbreak", "null\x00", 1])
def test_invalid_default_names(example, name):
    with pytest.raises(DefaultValidationError):
        create_default(example, name, {})


@pytest.mark.parametrize("key,value", [
    ("power", True), ("power", "2.5"), ("power", float("nan")),
    ("power", float("inf")), ("power", -0.1), ("power", 10.1),
    ("count", 1.5), ("count", False), ("count", 10**1000),
    ("enabled", 1), ("enabled", "false"),
    ("mode", "unregistered"), ("source", {"code": "no"}),
    ("result", 42), ("unknown", 7),
])
def test_parameter_schema_validation(example, key, value):
    with pytest.raises(DefaultValidationError):
        create_default(example, "Test", {key: value})


def test_integer_numbers_and_precise_float_settings_are_preserved(example):
    record = create_default(example, "Precise", {"count": 4.0, "power": 1.23456789})
    assert type(record["params"]["count"]) is int
    assert record["params"]["power"] == 1.23456789


@pytest.mark.parametrize("field,value", [
    ("format", "arbitrary"), ("version", 2), ("version", True),
    ("version", 1.0), ("scheme", "other"), ("params", []),
])
def test_incompatible_records_fail_before_application(example, field, value):
    record = create_default(example, "Test", {})
    record[field] = value
    with pytest.raises(DefaultValidationError):
        parse_default(example, record)


def test_missing_or_extra_record_fields_are_rejected(example):
    record = create_default(example, "Test", {})
    for modified in [{key: value for key, value in record.items() if key != "params"},
                     dict(record, extra="unused")]:
        with pytest.raises(DefaultValidationError):
            parse_default(example, modified)


@pytest.mark.parametrize("payload", [b"\xff", "{", "[]", "null", "1", '{"a":1,"a":2}',
                                    '{"power":NaN}', '{"power":Infinity}', "1" * 5000])
def test_invalid_json_is_rejected(example, payload):
    with pytest.raises(DefaultValidationError):
        parse_default(example, payload)


def test_oversized_file_and_browser_records_are_rejected(example):
    with pytest.raises(DefaultValidationError):
        parse_default(example, " " * (MAX_PAYLOAD_BYTES + 1))
    record = create_default(example, "Test", {})
    record["params"]["source"] = "x" * MAX_PAYLOAD_BYTES
    with pytest.raises(DefaultValidationError):
        parse_default(example, record)
    with pytest.raises(DefaultValidationError):
        create_default(example, "Test", record["params"])


def test_parser_does_not_mutate_browser_record(example):
    record = create_default(example, "Test", {})
    del record["params"]["enabled"]
    original = copy.deepcopy(record)
    loaded = parse_default(example, record)
    assert record == original
    assert loaded["params"]["enabled"] is False


@pytest.mark.parametrize("name,expected", [
    ("2026-10-13 experiment", "2026-10-13 experiment.json"),
    ("실험 설정", "실험 설정.json"),
    ("folder/name:one?", "folder-name-one-.json"),
    ("CON", "default-CON.json"), ("nul.txt", "default-nul.txt.json"),
    ("Today.json", "Today.json"), ("....", "default.json"),
])
def test_download_filenames_are_portable(name, expected):
    assert default_filename(name) == expected
