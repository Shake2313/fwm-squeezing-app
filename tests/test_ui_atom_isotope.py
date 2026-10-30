"""Atom / Isotope selectors (`ParamSpec.atom_menu`) keep the stored parameter contract."""
import json

import pytest
from streamlit.testing.v1 import AppTest

from gabes import schemes


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

SCHEME_APP = '''
import streamlit as st
from gabes import schemes
from gabes_ui.controls import render_rail, current_params
scheme = schemes.get({name!r})
specs = {{s.name: s for s in scheme.param_schema()}}
render_rail(st, scheme, [specs[{knob!r}], specs["line"]])
st.json(current_params({name!r}, scheme))
'''


def test_isotope_selection_keeps_legacy_physics_parameters():
    app = AppTest.from_string(APP).run()
    assert not app.exception
    assert app.selectbox[0].options == ["Rb", "Cs"]
    assert app.selectbox[1].options == ["Natural", "85", "87", "Custom"]
    for label, species, composition in [
        ("85", "⁸⁵Rb", "Natural"), ("87", "⁸⁷Rb", "Natural"),
        ("Custom", "Rb (natural)", "Custom"), ("Natural", "Rb (natural)", "Natural")]:
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
    assert app.selectbox[1].value == "Natural"


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


def test_lambda_atom_menu_covers_rb_cs_and_generic():
    app = AppTest.from_string(SCHEME_APP.format(name="lambda", knob="species")).run()
    assert not app.exception
    assert app.selectbox[0].options == ["Rb", "Cs", "Generic"]
    assert app.selectbox[1].options == ["Natural", "85", "87"]
    app.selectbox[1].select("87").run()
    assert json.loads(app.json[0].value)["species"] == "⁸⁷Rb"
    app.selectbox[0].select("Generic").run()
    assert [box.label for box in app.selectbox] == ["Atom", "Transition line"]
    assert json.loads(app.json[0].value)["species"] == "Generic"


def test_magneto_shows_a_fixed_atom_and_its_isotopes():
    app = AppTest.from_string(SCHEME_APP.format(name="magneto", knob="isotope")).run()
    assert not app.exception
    assert [box.label for box in app.selectbox] == ["Atom", "Isotope", "Transition line"]
    assert app.selectbox[0].options == ["Rb"] and app.selectbox[0].disabled
    assert app.selectbox[1].options == ["85", "87"]
    app.selectbox[1].select("87").run()
    assert json.loads(app.json[0].value)["isotope"] == "Rb87"


@pytest.mark.parametrize("scheme", schemes.all_schemes(), ids=lambda s: s.name)
def test_atom_menus_reach_every_stored_choice(scheme):
    specs = {sp.name: sp for sp in scheme.param_schema()}
    for sp in specs.values():
        if not sp.atom_menu:
            continue
        assert {row[2] for row in sp.atom_menu} == set(sp.choices)
        for _atom, _isotope, _value, companions in sp.atom_menu:
            for name, value in companions.items():
                assert value in specs[name].choices and specs[name].hidden
