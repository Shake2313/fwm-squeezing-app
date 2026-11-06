"""Automatic main condition editors share physical state with presets and exports."""
import pytest
from streamlit.testing.v1 import AppTest

from gabes.schemes.sas import SASScheme
from gabes_ui.conditions import _field_value


@pytest.mark.parametrize("scheme_name,canonical,selection,expected", [
    ("lambda", "species", "85", "⁸⁵Rb"),
    ("magneto", "isotope", "87", "Rb87"),
    ("psr", "isotope", "85", "Rb85"),
])
def test_shared_atom_editors_write_the_declared_scheme_field(scheme_name, canonical, selection, expected):
    app = AppTest.from_string(f"""
import streamlit as st
from gabes import schemes
from gabes_ui.conditions import render_conditions, render_power_inputs
scheme = schemes.get({scheme_name!r})
for k, v in scheme.defaults().items():
    st.session_state.setdefault(scheme.name+'__'+k, v)
render_conditions(scheme,scheme.param_schema())
render_power_inputs(scheme,scheme.param_schema())
""").run()
    assert not app.exception and not app.slider
    next(w for w in app.selectbox if w.label == "Isotope").select(selection).run()
    assert not app.exception
    assert app.session_state[scheme_name+'__'+canonical] == expected


def test_typed_power_keeps_precision_converts_units_and_enforces_domain():
    spec = next(sp for sp in SASScheme().param_schema() if sp.name == "pump_power_mw")
    assert _field_value(spec, "507.89 µW") == pytest.approx(.50789)
    for value in ("-1", "2.1", "1 GHz", "1e10000", "?", ""):
        with pytest.raises(ValueError):
            _field_value(spec, value)


def test_automatic_laser_inputs_validate_independently_and_presets_refresh_values():
    app = AppTest.from_string("""
import streamlit as st
from gabes.schemes.sas import SASScheme
from gabes_ui.conditions import render_conditions, render_power_inputs
scheme = SASScheme()
for k, v in scheme.defaults().items():
    st.session_state.setdefault('sas__'+k, v)
st.button('Load OD', on_click=lambda: st.session_state.update(sas__pump_power_mw=0., sas__probe_power_uw=0.))
params=render_conditions(scheme,scheme.param_schema())
render_power_inputs(scheme,scheme.param_schema())
""").run()
    def field(label):
        return next(w for w in app.text_input if w.label == label)
    assert not app.exception
    assert not app.slider
    assert not any(b.label == 'Apply' for b in app.button)
    # The browser must hydrate a number, rather than a blank default followed
    # by a separate state update that can trigger automatic validation early.
    assert field("Pump power (mW)").proto.default == "0.5"
    assert field("Probe power (µW)").proto.default == "0"
    field("Pump power (mW)").set_value("0.50789").run()
    field("Probe power (µW)").set_value("70 nW").run()
    field("Probe power (µW)").set_value("999 µW").run()
    assert app.error and app.session_state['sas__probe_power_uw'] == pytest.approx(.07)
    app.button[0].click().run()
    assert not app.error and field("Pump power (mW)").value == "0"
    field("Pump power (mW)").set_value("0.50789").run()
    field("Probe power (µW)").set_value("70 nW").run()
    assert app.session_state["sas__pump_power_mw"] == .50789
    assert app.session_state["sas__probe_power_uw"] == pytest.approx(.07)
    app.button[0].click().run()
    assert not app.exception
    assert field("Pump power (mW)").value == "0"
    assert field("Probe power (µW)").value == "0"
    assert not any("fraction" in w.label for w in app.text_input)
    isotope=next(w for w in app.selectbox if w.label == "Isotope")
    isotope.select("Custom").run()
    field("⁸⁵Rb fraction (Natural 0.7217)").set_value(".4").run()
    assert app.session_state["sas__rb_composition"] == "Custom"
    assert app.session_state["sas__rb85_fraction"] == .4
    next(w for w in app.selectbox if w.label == "Isotope").select("Natural").run()
    assert not any("fraction" in w.label for w in app.text_input)
    assert app.session_state["sas__rb_composition"] == "Natural"
    atom=next(w for w in app.selectbox if w.label == "Atom")
    atom.select("Cs").run()
    assert app.session_state["sas__species"] == "¹³³Cs"
    assert not any(w.label == "Isotope" for w in app.selectbox)



def test_live_native_defaults_survive_rerenders_invalid_drafts_and_saved_settings(monkeypatch):
    # AppTest keeps Python widget values even when the serialized browser default
    # is blank/false/first-option. Check the actual hydration contract as well.
    from streamlit.elements.lib import policies
    monkeypatch.setattr(policies, "_shown_default_value_warning", False)
    app = AppTest.from_string("""
import streamlit as st
from gabes import schemes
from gabes_ui.conditions import render_conditions, render_power_inputs
scheme = schemes.get("psr")
for key, value in scheme.defaults().items():
    st.session_state.setdefault("psr__"+key, value)
st.button("Load saved", on_click=lambda: st.session_state.update(
    psr__intensity_mw_cm2=.25, psr__doppler=True, psr__resolution="Balanced"))
render_conditions(scheme, scheme.param_schema())
render_power_inputs(scheme, scheme.param_schema())
""").run()

    def numeric():
        return next(w for w in app.text_input if w.label.startswith("Carrier intensity"))

    def doppler():
        return next(w for w in app.checkbox if w.label == "Doppler averaging")

    def resolution():
        return next(w for w in app.selectbox if w.label == "Resolution")

    def assert_browser_defaults():
        assert not app.exception and not app.warning
        assert numeric().proto.default == numeric().value
        assert doppler().proto.default == doppler().value
        assert resolution().proto.default == resolution().options.index(resolution().value)

    assert_browser_defaults()
    app.run()
    assert_browser_defaults()
    assert numeric().proto.default == "0.1"
    assert doppler().proto.default is True
    assert resolution().proto.default == 1  # Balanced, after Quick.
    identity = numeric().proto.id
    numeric().set_value("0.50789").run()
    assert_browser_defaults()
    assert numeric().proto.id == identity
    assert app.session_state["psr__intensity_mw_cm2"] == .50789
    numeric().set_value("-0.1").run()
    assert_browser_defaults()
    assert app.error and numeric().proto.default == "-0.1"
    assert numeric().proto.id == identity
    assert app.session_state["psr__intensity_mw_cm2"] == .50789
    doppler().uncheck().run()
    resolution().select("Quick").run()
    assert_browser_defaults()
    assert numeric().value == "-0.1"
    assert app.session_state["psr__doppler"] is False
    assert app.session_state["psr__resolution"] == "Quick"
    app.button[0].click().run()
    assert_browser_defaults()
    assert not app.error and numeric().value == "0.25"
    assert doppler().value is True and resolution().value == "Balanced"
    app.run()
    assert_browser_defaults()
