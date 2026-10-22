"""Apparatus profiles reach the page's live controls and survive reruns."""
from dataclasses import fields
import math
from pathlib import Path
import sys

import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

import sabes_page as page  # noqa: E402
from sabes.presets import (LAB_PROFILE, PAPER_PROFILE, power_targets,  # noqa: E402
                           profile_calibration, profile_settings)


@pytest.fixture
def state(monkeypatch):
    state = {}
    monkeypatch.setattr(page.st, "session_state", state)
    return state


def _assert_settings_match(actual, expected):
    for field in fields(expected):
        value = getattr(expected, field.name)
        current = getattr(actual, field.name)
        if isinstance(value, float) and math.isnan(value):
            assert math.isnan(current), field.name
        elif isinstance(value, (int, float, tuple)):
            assert current == pytest.approx(value), field.name
        else:
            assert current == value, field.name


def test_first_session_selects_lab_and_passes_its_profile_to_live_settings(state):
    page._seed_session_state()
    assert page._selected_profile() == LAB_PROFILE
    settings, detection = page._current()
    expected_settings, expected_detection = profile_settings(LAB_PROFILE)
    _assert_settings_match(settings, expected_settings)
    _assert_settings_match(detection, expected_detection)
    assert state[page._key("use_recorded_waists")] is True
    assert state[page._key("use_recorded_detection")] is True
    assert state[page._key("_version")] == page.SETTINGS_VERSION


def test_live_calibration_uses_lab_nel_then_keeps_only_allowed_hand_edits(state):
    page._seed_session_state()
    calibration = page._session_calibration()
    assert calibration.value("t_noise_eater") == 0.928
    assert calibration.get("t_noise_eater").provenance == "lab"
    page._seed_transmission_state(calibration)
    state[page._key("t_noise_eater")] = 0.81
    # A retired duplicate EOM transmission must never enter the active model.
    state[page._key("t_eom")] = 0.1
    edited = page._session_calibration()
    assert edited.value("t_noise_eater") == 0.81
    assert edited.get("t_noise_eater").provenance == "hand"
    assert edited.value("t_eom") == calibration.value("t_eom")
    assert calibration.value("t_noise_eater") == 0.928


def test_switching_profile_clears_hand_edits_and_optic_state_but_keeps_selection(state):
    page._seed_session_state()
    state[page._key("eom_rf_dbm")] = 10.0
    state[page._key("t_noise_eater")] = 0.2
    state[page._key("selected")] = "bpd"
    state[page._key("open_optic")] = "bpd"
    state[page._key("_solver_error")] = "old error"
    for part in page.layout.LAYOUTS:
        state[page._probe_key(part.key)] = [20.0, 30.0]
        state[page._key(f"inspect_{part.key}")] = "bpd"
    state[page._key("profile")] = PAPER_PROFILE
    page._profile_changed()
    page._seed_session_state()
    assert page._selected_profile() == PAPER_PROFILE
    assert page._key("t_noise_eater") not in state
    for name in ("selected", "open_optic", "_solver_error"):
        assert page._key(name) not in state
    for part in page.layout.LAYOUTS:
        assert page._probe_key(part.key) not in state
        assert page._key(f"inspect_{part.key}") not in state
    settings, detection = page._current()
    expected_settings, expected_detection = profile_settings(PAPER_PROFILE)
    _assert_settings_match(settings, expected_settings)
    _assert_settings_match(detection, expected_detection)
    assert settings.pump_waist_override_um is None
    assert settings.seed_waist_override_um is None
    assert detection.post_cell_transmission_override is None
    assert detection.qe_pct_override is None
    assert page._session_calibration().to_dict() == profile_calibration(
        PAPER_PROFILE).to_dict()


def test_same_version_rerun_preserves_live_controls_and_overrides(state):
    page._seed_session_state()
    edits = {"eom_rf_dbm": 13.5, "cell_temp_c": 119.0,
             "pump_waist_override_um": 570.0,
             "post_cell_transmission_override": 0.92,
             "qe_pct_override": 87.5, "t_lens_probe": 0.88}
    for name, value in edits.items():
        state[page._key(name)] = value
    page._seed_session_state()
    page._seed_transmission_state(page._session_calibration())
    settings, detection = page._current()
    assert settings.eom_rf_dbm == 13.5
    assert settings.cell_temp_c == 119.0
    assert settings.pump_waist_override_um == 570.0
    assert detection.post_cell_transmission_override == 0.92
    assert detection.qe_pct_override == 87.5
    assert page._session_calibration().value("t_lens_probe") == 0.88
    # Turning off recorded inputs must restore the underlying optical model,
    # while retaining the numeric widget contents for a later re-enable.
    state[page._key("use_recorded_waists")] = False
    state[page._key("use_recorded_detection")] = False
    settings, detection = page._current()
    assert settings.pump_waist_override_um is None
    assert detection.post_cell_transmission_override is None
    assert state[page._key("pump_waist_override_um")] == 570.0
    assert state[page._key("post_cell_transmission_override")] == 0.92


@pytest.mark.parametrize("profile", (LAB_PROFILE, PAPER_PROFILE))
def test_reset_restores_the_selected_profile_without_switching_it(state, profile):
    state[page._key("profile")] = profile
    page._seed_session_state()
    state[page._key("eom_rf_dbm")] = 9.5
    state[page._key("use_recorded_waists")] = profile != LAB_PROFILE
    state[page._key("t_noise_eater")] = 0.4
    page._reset_defaults()
    assert page._selected_profile() == profile
    assert page._key("t_noise_eater") not in state
    actual_settings, actual_detection = page._current()
    expected_settings, expected_detection = profile_settings(profile)
    _assert_settings_match(actual_settings, expected_settings)
    _assert_settings_match(actual_detection, expected_detection)


def test_inverse_controls_use_lab_targets_current_calibration_and_sequential_angles(
        state, monkeypatch):
    page._seed_session_state()
    state[page._key("t_noise_eater")] = 0.84
    calls = []

    def solve(name, answer):
        def run(target, settings, calibration):
            calls.append((name, target, settings, calibration))
            return answer
        return run

    monkeypatch.setattr(page, "solve_split_angle_deg", solve("split", 11.123456))
    monkeypatch.setattr(page, "solve_seed_polarizer_deg", solve("polarizer", 22.234567))
    monkeypatch.setattr(page, "solve_seed_trim_deg", solve("trim", 33.345678))
    page._apply_solved_knobs()
    assert [item[0] for item in calls] == ["split", "polarizer", "trim"]
    assert [item[1] for item in calls] == pytest.approx(power_targets(LAB_PROFILE))
    assert all(item[3].value("t_noise_eater") == 0.84 for item in calls)
    assert calls[1][2].hwp_split_deg == 11.123456
    assert calls[2][2].seed_polarizer_deg == 22.234567
    assert state[page._key("hwp_split_deg")] == 11.1235
    assert state[page._key("seed_polarizer_deg")] == 22.2346
    assert state[page._key("seed_trim_hwp_deg")] == 33.3457


def test_optional_paper_controls_remain_numeric_but_panels_accept_none(
        state, monkeypatch):
    state[page._key("profile")] = PAPER_PROFILE
    page._seed_session_state()
    rendered, controls = [], []
    monkeypatch.setattr(page.st, "markdown", lambda text, **kwargs: rendered.append(text))
    monkeypatch.setattr(page.st, "caption", lambda *args, **kwargs: None)
    monkeypatch.setattr(page.st, "button", lambda *args, **kwargs: False)

    class Container:
        def number_input(self, label, **kwargs):
            assert isinstance(kwargs["value"], (int, float))
            assert kwargs["key"].startswith(page._widget_key(""))
            controls.append((label, kwargs))

    settings, detection = page._current()
    for node in ("seed_l2", "pump_l2", "bpd"):
        page._render_optic_panel(page.layout.PART2, node, settings, detection)
    for name in ("pump_waist_override_um", "seed_waist_override_um",
                 "post_cell_transmission_override", "qe_pct_override"):
        page._render_control(Container(), name)
    assert len(controls) == 4
    assert all(kwargs["disabled"] for _label, kwargs in controls)
    assert sum(text.count("model default") for text in rendered) == 4


def test_widget_callbacks_preserve_model_edits_after_transient_keys_disappear(state):
    page._seed_session_state()
    for name, value in (("post_cell_transmission_override", 0.92),
                        ("use_recorded_detection", False)):
        state[page._widget_key(name)] = value
        page._copy_widget(name)
        del state[page._widget_key(name)]
    page._seed_session_state()
    assert state[page._key("post_cell_transmission_override")] == 0.92
    assert state[page._key("use_recorded_detection")] is False
    assert page._current()[1].post_cell_transmission_override is None
    state[page._widget_key("post_cell_transmission_override")] = 0.5
    page._reset_defaults()
    assert page._widget_key("post_cell_transmission_override") not in state
    assert state[page._key("post_cell_transmission_override")] == 0.9709
    assert state[page._key("use_recorded_detection")] is True


def test_real_streamlit_editor_keeps_edits_when_widgets_are_removed_and_recreated():
    from streamlit.testing.v1 import AppTest

    app = AppTest.from_string("""
import streamlit as st
import sabes_page as page
page._seed_session_state()
st.checkbox("Show editor", key="show_editor", value=True)
st.button("Reset", on_click=page._reset_defaults)
if st.session_state["show_editor"]:
    page._render_bool_control("Use recorded efficiency", "use_recorded_detection")
    page._render_control(st, "post_cell_transmission_override")
settings, detection = page._current()
st.write("Active transmission", detection.post_cell_transmission_override)
""", default_timeout=10)
    app.run()
    assert not app.exception
    field = page._widget_key("post_cell_transmission_override")
    flag = page._widget_key("use_recorded_detection")
    assert app.number_input(key=field).value == 0.9709
    assert app.checkbox(key=flag).value is True
    app.number_input(key=field).set_value(0.92).run()
    assert not app.exception
    assert app.session_state[page._key("post_cell_transmission_override")] == 0.92

    app.checkbox(key="show_editor").uncheck().run()
    assert not app.exception
    assert len(app.number_input) == 0
    assert field not in app.session_state
    assert flag not in app.session_state
    assert app.session_state[page._key("post_cell_transmission_override")] == 0.92
    app.checkbox(key="show_editor").check().run()
    assert not app.exception
    assert app.number_input(key=field).value == 0.92

    app.checkbox(key=flag).uncheck().run()
    assert app.number_input(key=field).disabled is True
    assert app.session_state[page._key("use_recorded_detection")] is False
    app.checkbox(key="show_editor").uncheck().run()
    app.checkbox(key="show_editor").check().run()
    assert not app.exception
    assert app.checkbox(key=flag).value is False
    assert app.number_input(key=field).value == 0.92

    app.checkbox(key="show_editor").uncheck().run()
    app.button[0].click().run()
    assert app.session_state[page._key("profile")] == LAB_PROFILE
    assert app.session_state[page._key("post_cell_transmission_override")] == 0.9709
    app.checkbox(key="show_editor").check().run()
    assert not app.exception
    assert app.number_input(key=field).value == 0.9709
    assert app.checkbox(key=flag).value is True
    assert not app.warning
