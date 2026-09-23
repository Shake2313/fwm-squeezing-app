"""ScrubField gesture handling (gabes_ui.scrub)."""
import pytest

from gabes_ui.scrub import decode_gesture, resolve_gesture, state_key

KNOB = dict(unit="mW", vmin=0.0, vmax=2.0, integer=False, default=0.5)


def test_component_keys_drop_the_reserved_double_underscore():
    # Streamlit refuses "__" inside a bidirectional component's id.
    assert state_key("sas__pump_power_mw") == "scrub-sas-pump_power_mw"
    assert "__" not in state_key("fwm__excess_noise_slope")


def test_a_gesture_arrives_as_json_text():
    assert decode_gesture('{"n": 3, "v": 1.25}') == {"n": 3, "v": 1.25}
    assert decode_gesture({"n": 1, "r": True}) == {"n": 1, "r": True}
    assert decode_gesture(None) is None
    assert decode_gesture("[]") is None


def test_a_dragged_value_is_taken_as_it_comes():
    assert resolve_gesture({"n": 1, "v": 1.25}, **KNOB) == (1.25, None)


def test_typed_text_is_parsed_in_the_knobs_unit():
    assert resolve_gesture({"n": 1, "t": "1.25"}, **KNOB) == (1.25, None)
    assert resolve_gesture({"n": 1, "t": "1200 uW"}, **KNOB) == (1.2, None)


def test_out_of_range_entry_lands_on_the_boundary_and_says_so():
    value, note = resolve_gesture({"n": 1, "t": "1.25 W"}, **KNOB)
    assert value == 2.0
    assert note == "clamped to 2"


def test_an_entry_in_the_wrong_quantity_leaves_the_knob_alone():
    value, note = resolve_gesture({"n": 1, "t": "3 MHz"}, **KNOB)
    assert value is None
    assert note == "expected mW"


def test_reset_restores_the_default():
    assert resolve_gesture({"n": 1, "r": True}, **KNOB) == (0.5, None)


def test_integer_knobs_stay_integers():
    knob = dict(unit="", vmin=201, vmax=4001, integer=True, default=1401)
    value, _ = resolve_gesture({"n": 1, "v": 1400.6}, **knob)
    assert value == 1401 and isinstance(value, int)


def test_nothing_to_do_is_not_an_error():
    assert resolve_gesture(None, **KNOB) == (None, None)
    assert resolve_gesture({"n": 1}, **KNOB) == (None, None)


@pytest.mark.parametrize("text,expected", [("2", 2.0), ("0", 0.0), ("-1", 0.0)])
def test_the_range_is_always_respected(text, expected):
    value, _ = resolve_gesture({"n": 1, "t": text}, **KNOB)
    assert value == expected
