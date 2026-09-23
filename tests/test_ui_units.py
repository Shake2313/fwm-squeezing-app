"""Typed numeric entry for the ScrubField (gabes_ui.units)."""
import pytest

from gabes_ui.units import (
    UnitError,
    clamp,
    convert,
    decimals_for,
    display_text,
    parse_entry,
)

# (typed text, knob unit, value in the knob's unit)
ACCEPTED = [
    ("1.5", "MHz", 1.5),
    ("-2e3", "kHz", -2000.0),
    ("1,234.5", "mm", 1234.5),
    ("−7", "°C", -7.0),          # a real minus sign, as a keyboard may give
    ("  12  ", "", 12.0),
    ("0.4 GHz", "MHz", 400.0),
    ("400MHz", "GHz", 0.4),
    ("1 MHz", "kHz", 1000.0),
    ("5 uW", "mW", 0.005),                 # ASCII u for micro
    ("5 μW", "mW", 0.005),            # Greek mu for micro
    ("5 µW", "mW", 0.005),            # micro sign
    ("2 W", "mW", 2000.0),
    ("30 cm", "mm", 300.0),
    ("780 nm", "µm", 0.78),
    ("1 m", "mm", 1000.0),
    ("300 K", "°C", pytest.approx(26.85)),
    ("20 °C", "K", pytest.approx(293.15)),
    ("20 C", "K", pytest.approx(293.15)),
    ("1 G", "µT", pytest.approx(100.0)),
    ("50 µT", "G", pytest.approx(0.5)),
    ("0.5 mT", "µT", 500.0),
    ("3 ms", "ns", 3e6),
    ("1 rad", "mrad", 1000.0),
    ("1 W/cm²", "mW/cm²", 1000.0),
    ("1 Torr", "Torr", 1.0),
    ("2 mTorr", "Torr", pytest.approx(0.002)),
    ("3 dB", "dB", 3.0),
    ("7 %", "%", 7.0),
]

REJECTED = [
    ("", "MHz"),
    ("abc", "MHz"),
    ("MHz", "MHz"),
    ("1.5 mm", "MHz"),
    ("1.5 K", "MHz"),
    ("3 dBm", "dB"),          # a different quantity that merely looks related
    ("1 mK", "°C"),      # offset scale plus prefix: refuse, never mis-shift
    ("5 MHz", ""),            # a plain number knob takes no unit
]


@pytest.mark.parametrize("text,unit,expected", ACCEPTED)
def test_typed_entry_lands_in_the_knobs_own_unit(text, unit, expected):
    value, error = parse_entry(text, unit)
    assert error is None
    assert value == expected


@pytest.mark.parametrize("text,unit", REJECTED)
def test_entry_that_cannot_be_trusted_is_refused_with_a_short_reason(text, unit):
    value, error = parse_entry(text, unit)
    assert value is None
    assert error and len(error) < 40


def test_decibels_do_not_decompose_into_a_deci_prefix():
    assert convert(3.0, "dB", "dB") == 3.0
    with pytest.raises(UnitError):
        convert(3.0, "B", "dB")


def test_a_lone_g_is_gauss_and_a_prefixed_one_is_giga():
    assert convert(1.0, "G", "T") == pytest.approx(1e-4)
    assert convert(1.0, "GHz", "MHz") == pytest.approx(1000.0)


def test_clamp_reports_only_when_it_moved_the_value():
    assert clamp(5.0, 0.0, 10.0) == (5.0, None)
    assert clamp(-1.0, 0.0, 10.0) == (0.0, "clamped to 0")
    assert clamp(12.0, 0.0, 10.0) == (10.0, "clamped to 10")
    assert clamp(12.0, None, None) == (12.0, None)


def test_decimals_come_from_the_slider_format():
    assert decimals_for("%.3f") == 3
    assert decimals_for("%.0f") == 0
    assert decimals_for(None) == 0


def test_display_keeps_the_slider_format_until_precision_would_be_lost():
    assert display_text(40.0, None) == "40"
    assert display_text(40.0, "%.1f") == "40.0"
    assert display_text(86.94, "%.2f") == "86.94"
    # A typed value finer than the format must still be readable in full.
    assert display_text(1.234567, "%.2f") == "1.234567"
    assert display_text(None, "%.2f") == ""
