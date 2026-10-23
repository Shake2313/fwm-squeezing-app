"""Named defaults apply complete snapshots and acknowledge persistent saves."""
from io import BytesIO

import pytest

from gabes.schemes.base import ParamSpec, Preset, Scheme
from gabes_ui import defaults
from gabes_ui.defaults_store import create_default, dump_default


class SampleScheme(Scheme):
    name = "sample"

    def param_schema(self):
        return [ParamSpec("mode", "Mode", "Setup", "A", choices=("A", "B")),
                ParamSpec("power", "Power", "Setup", 1.0, 0.0, 10.0),
                ParamSpec("flag", "Flag", "Setup", False, control="checkbox"),
                ParamSpec("provenance", "Source", "Setup", "", hidden=True)]

    def presets(self):
        return [Preset("Experiment", {"mode": "B", "power": 4.0, "flag": True}),
                Preset("Reference", {"power": 2.0})]

    def compute(self, params):
        return {}

    def observables(self, raw, params, include_figures=True):
        return {}


@pytest.fixture
def state(monkeypatch):
    values = {}
    monkeypatch.setattr(defaults.st, "session_state", values)
    return values


def test_switching_defaults_resets_hidden_and_previously_enabled_inputs(state):
    scheme = SampleScheme()
    records = defaults.builtin_defaults(scheme)
    defaults._apply(scheme, records["builtin:0"])
    state["sample__provenance"] = "local experiment"
    defaults._apply(scheme, records["builtin:1"])
    assert defaults._params(scheme) == {
        "mode": "A", "power": 2.0, "flag": False, "provenance": ""}


def test_storage_acknowledgement_only_publishes_successful_save(state):
    scheme = SampleScheme()
    name_key = defaults._key(scheme, "default_name")
    state[name_key] = "My experiment"
    defaults._save_browser(scheme, name_key, defaults.builtin_defaults(scheme))
    command = state[defaults._key(scheme, "defaults_command")]
    assert not defaults._saved(scheme)
    defaults._receive(scheme, dict(id=command["id"], action="save",
                                   records=[], error="Browser storage is full."))
    assert not defaults._saved(scheme)
    assert defaults._key(scheme, "defaults_command") not in state
    defaults._save_browser(scheme, name_key, defaults.builtin_defaults(scheme))
    command = state[defaults._key(scheme, "defaults_command")]
    defaults._receive(scheme, dict(id=command["id"], action="save",
                                   records=[command["record"]], error=None))
    assert defaults._saved(scheme)["My experiment"]["params"] == scheme.defaults()


def test_browser_reload_restores_choices_without_overwriting_current_parameters(state):
    scheme = SampleScheme()
    record = create_default(scheme, "Stored", {"power": 7.0})
    state["sample__power"] = 3.0
    defaults._receive(scheme, dict(id="load-fresh", action="load", records=[record], error=None))
    assert defaults._saved(scheme)["Stored"] == record
    assert state["sample__power"] == 3.0


def test_local_import_restores_values_and_offers_named_choice(state):
    scheme = SampleScheme()
    record = create_default(scheme, "Local snapshot", {"mode": "B", "power": 6.25})
    upload_key = defaults._key(scheme, "import_default")
    state[upload_key] = BytesIO(dump_default(record).encode())
    defaults._import(scheme, upload_key, defaults.builtin_defaults(scheme))
    assert defaults._params(scheme) == record["params"]
    assert defaults._saved(scheme)["Local snapshot"] == record


@pytest.mark.parametrize("payload", [b"not json", b'{"scheme":"other"}'])
def test_rejected_import_leaves_current_settings_untouched(state, payload):
    scheme = SampleScheme()
    state["sample__power"] = 8.0
    upload_key = defaults._key(scheme, "import_default")
    state[upload_key] = BytesIO(payload)
    defaults._import(scheme, upload_key, defaults.builtin_defaults(scheme))
    assert state["sample__power"] == 8.0
    assert not defaults._saved(scheme)
    assert state[defaults._key(scheme, "defaults_notice")][1] is True


def test_saved_names_cannot_shadow_builtins(state):
    scheme = SampleScheme()
    state["name"] = "Experiment"
    defaults._save_browser(scheme, "name", defaults.builtin_defaults(scheme))
    assert defaults._key(scheme, "defaults_command") not in state
    assert state[defaults._key(scheme, "defaults_notice")][1] is True


def test_real_dropdown_loads_before_mode_widgets_and_marks_edited_settings():
    from streamlit.testing.v1 import AppTest
    app = AppTest.from_string(
        "import streamlit as st\n"
        "from tests.test_ui_defaults import SampleScheme\n"
        "from gabes_ui.defaults import render_defaults\n"
        "scheme=SampleScheme()\n"
        "render_defaults(st,scheme)\n"
        "st.selectbox('Mode',['A','B'],key='sample__mode')\n"
        "st.number_input('Power',0.,10.,key='sample__power')\n"
        "st.checkbox('Flag',key='sample__flag')\n",
    ).run()
    assert not app.exception
    app.selectbox(key="sample___named_default").select("builtin:0").run()
    assert not app.exception
    assert app.selectbox(key="sample__mode").value == "B"
    assert app.number_input(key="sample__power").value == 4.0
    assert app.checkbox(key="sample__flag").value is True
    app.selectbox(key="sample___named_default").select("builtin:1").run()
    assert not app.exception
    assert app.selectbox(key="sample__mode").value == "A"
    assert app.checkbox(key="sample__flag").value is False
    app.number_input(key="sample__power").set_value(3.25).run()
    assert not app.exception
    assert app.selectbox(key="sample___named_default").value is None


def test_browser_ack_selects_saved_name_when_values_match_a_builtin(monkeypatch):
    from streamlit.testing.v1 import AppTest

    # AppTest's script installs its fake bridge in this shared imported module.
    monkeypatch.setattr(defaults, "render_browser_store", defaults.render_browser_store)
    app = AppTest.from_string(
        "import streamlit as st\n"
        "from tests.test_ui_defaults import SampleScheme\n"
        "from gabes_ui import defaults\n"
        "scheme=SampleScheme()\n"
        "def bridge(container,scheme_name,command=None,on_event=None):\n"
        "    if command:\n"
        "        def ack():\n"
        "            on_event(dict(id=command['id'],action='save',"
        "records=[command['record']],error=None))\n"
        "        container.button('Browser ACK',on_click=ack)\n"
        "defaults.render_browser_store=bridge\n"
        "defaults.render_defaults(st,scheme)\n"
        "st.selectbox('Mode',['A','B'],key='sample__mode')\n"
        "st.number_input('Power',0.,10.,key='sample__power')\n"
        "st.checkbox('Flag',key='sample__flag')\n"
    ).run()
    app.selectbox(key="sample___named_default").select("builtin:1").run()
    app.text_input(key="sample___default_name").set_value("Named reference").run()
    app.button(key="sample___save_browser").click().run()
    assert app.selectbox(key="sample___named_default").value == "builtin:1"
    next(button for button in app.button if button.label == "Browser ACK").click().run()
    assert not app.exception
    assert app.session_state["sample___defaults_revision"] == 1
    assert app.selectbox(key="sample___named_default:1").value == "saved:Named reference"

    # A stale client-side widget value must not replace the logical active name.
    app.session_state["sample___named_default"] = "builtin:1"
    app.run()
    assert not app.exception
    assert app.selectbox(key="sample___named_default:1").value == "saved:Named reference"
