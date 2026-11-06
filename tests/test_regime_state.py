"""Regime widgets preserve canonical parameters when schemes are unmounted."""
from streamlit.proto.WidgetStates_pb2 import WidgetState, WidgetStates
from streamlit.testing.v1 import AppTest
from streamlit.testing.v1.element_tree import get_widget_state


_APP = """
import streamlit as st
from gabes import schemes
from gabes_ui.controls import skey
from gabes_ui.shell import render_regime_selector

name = st.selectbox("Scheme", ["rydberg_eit", "lambda"])
scheme = schemes.get(name)
for param, value in scheme.defaults().items():
    key = skey(name, param)
    st.session_state[key] = st.session_state.get(key, value)
spec = next(sp for sp in scheme.param_schema() if sp.applies_defaults)

def load_preset():
    st.session_state["rydberg_eit__view"] = "EIT"
    st.session_state["rydberg_eit__probe_power_uw"] = 8.5

st.button("Load saved EIT", on_click=load_preset)
st.button("Migrate", on_click=lambda: st.session_state.update(use_legacy=False))
if st.session_state.get("use_legacy", False):
    st.segmented_control(spec.label, list(spec.choices), key=skey(name, spec.name))
else:
    render_regime_selector(name, spec, scheme)
"""


_UNCHANGED = object()


def _run(app, selection=_UNCHANGED):
    """Serialize native segmented controls that AppTest treats as feedback widgets."""
    states = WidgetStates()
    for node in app._tree:
        if node.type == "button_group":
            state = WidgetState(id=node.id)
            value = node.value if selection is _UNCHANGED else selection
            values = [] if value is None else [value]
            labels = [option.content for option in node.options]
            state.int_array_value.data[:] = [labels.index(value) for value in values]
        else:
            state = get_widget_state(node)
        if state is not None:
            states.widgets.append(state)
    return app._run(states)


def test_regime_roundtrip_keeps_selected_mode_and_custom_parameters():
    app = AppTest.from_string(_APP).run()
    app.session_state["rydberg_eit__probe_power_uw"] = 9.
    _run(app, "AT electrometry")
    # Selecting a regime still applies its complete recommended set.
    assert app.session_state["rydberg_eit__view"] == "AT electrometry"
    assert app.session_state["rydberg_eit__probe_power_uw"] == 6.
    app.session_state["rydberg_eit__probe_power_uw"] = 8.75
    app.selectbox[0].select("lambda")
    _run(app)
    assert app.session_state["rydberg_eit__view"] == "AT electrometry"
    app.selectbox[0].select("rydberg_eit")
    _run(app)
    assert not app.exception
    assert app.get("button_group")[0].value == "AT electrometry"
    assert app.session_state["rydberg_eit__probe_power_uw"] == 8.75


def test_clearing_regime_keeps_custom_inputs_without_reapplying_defaults():
    app = AppTest.from_string(_APP).run()
    _run(app, "AT electrometry")
    app.session_state["rydberg_eit__probe_power_uw"] = 8.75
    _run(app, None)
    assert not app.exception
    assert app.session_state["rydberg_eit__view"] == "AT electrometry"
    assert app.get("button_group")[0].value == "AT electrometry"
    assert app.session_state["rydberg_eit__probe_power_uw"] == 8.75


def test_saved_preset_refreshes_regime_widget_without_loading_another_default_set():
    app = AppTest.from_string(_APP).run()
    _run(app, "AT electrometry")
    app.button[0].click()
    _run(app)
    assert not app.exception
    assert app.session_state["rydberg_eit__view"] == "EIT"
    assert app.get("button_group")[0].value == "EIT"
    assert app.session_state["rydberg_eit__probe_power_uw"] == 8.5


def test_existing_canonical_widget_migrates_without_losing_its_regime():
    app = AppTest.from_string(_APP)
    app.session_state["use_legacy"] = True
    app.run()
    _run(app, "AT electrometry")
    app.session_state["rydberg_eit__probe_power_uw"] = 8.75
    app.button[1].click()
    _run(app)
    assert not app.exception
    assert app.get("button_group")[0].key == "rydberg_eit___view_control"
    assert app.get("button_group")[0].value == "AT electrometry"
    app.selectbox[0].select("lambda")
    _run(app)
    app.selectbox[0].select("rydberg_eit")
    _run(app)
    assert app.session_state["rydberg_eit__view"] == "AT electrometry"
    assert app.session_state["rydberg_eit__probe_power_uw"] == 8.75
