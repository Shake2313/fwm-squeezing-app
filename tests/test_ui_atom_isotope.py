"""Element/isotope controls preserve the existing physical parameter contract."""
import json

from streamlit.testing.v1 import AppTest


APP = '''
import streamlit as st
from gabes import schemes
from gabes_ui.controls import render_rail, current_params, param_visible
scheme = schemes.get("sas")
for name, value in scheme.defaults().items():
    st.session_state.setdefault("sas__" + name, value)
def load_old():
    st.session_state["sas__species"] = "⁸⁷Rb"
    st.session_state["sas__rb_composition"] = "Custom"
st.button("Load old settings", on_click=load_old)
specs = {s.name:s for s in scheme.param_schema()}
render_rail(st, scheme, [specs["species"], specs["rb_composition"]])
st.json({**current_params("sas", scheme),
         "fraction_visible":param_visible("sas", specs["rb85_fraction"])})
'''


def test_isotope_selection_keeps_legacy_physics_parameters():
    app = AppTest.from_string(APP).run()
    assert not app.exception
    assert app.selectbox[0].options == ["Rb", "Cs"]
    assert app.selectbox[1].options == ["Nat", "85", "87", "Custom"]
    for label, species, composition in [
        ("85", "⁸⁵Rb", "Natural"), ("87", "⁸⁷Rb", "Natural"),
        ("Custom", "Rb (natural)", "Custom"), ("Nat", "Rb (natural)", "Natural")]:
        app.selectbox[1].select(label).run()
        assert not app.exception
        params = json.loads(app.json[0].value)
        assert (params["species"], params["rb_composition"]) == (species, composition)
        assert params["fraction_visible"] == (label == "Custom")
        assert app.selectbox[0].value == "Rb"
        assert app.selectbox[1].value == label
    app.selectbox[0].select("Cs").run()
    assert len(app.selectbox) == 1
    assert json.loads(app.json[0].value)["species"] == "¹³³Cs"
    app.selectbox[0].select("Rb").run()
    assert app.selectbox[1].value == "Nat"


def test_old_pure_isotope_settings_restore_both_selectors():
    app = AppTest.from_string(APP).run()
    app.selectbox[1].select("Custom").run()
    app.button[0].click().run()
    assert not app.exception
    assert app.selectbox[0].value == "Rb"
    assert app.selectbox[1].value == "87"
    params = json.loads(app.json[0].value)
    assert params["species"] == "⁸⁷Rb"
    assert not params["fraction_visible"]
