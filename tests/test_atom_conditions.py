"""Atom/cell condition groups use one popup and the existing physical state."""
from collections import Counter

import pytest
from streamlit.testing.v1 import AppTest


SCHEME_CASES = [
    ("sas", {}),
    ("sas", {"rb_composition": "Custom", "rb85_fraction": .4, "line": "D2", "ne_pressure_torr": 10.5}),
    ("sas", {"species": "\u00b9\u00b3\u00b3Cs"}),
    ("lambda", {}),
    ("lambda", {"species": "Generic"}),
    ("lambda", {"species": "\u00b9\u00b3\u00b3Cs"}),
    ("rydberg_eit", {}),
    ("rydberg_eit", {"temperature_model": "Separated"}),
    ("magneto", {}),
    ("magneto", {"cell_type": "Buffer gas cell"}),
    ("fwm", {}),
    ("fwm", {"mode": "Biphoton"}),
    ("fwm", {"mode": "Biphoton", "topology": "cascade_cs_btw"}),
    ("fwm", {"mode": "Biphoton", "topology": "diamond_generic"}),
    ("psr", {}),
    ("psr", {"calculation": "Rb87 D1 steady state"}),
    ("psr", {"calculation": "Rb87 D1 RF squeezing"}),
]


def _app(scheme_name, updates=None, *, rail=False):
    render = (
        """with st.sidebar:
    grouped = render_atom_conditions(scheme, specs)
params = render_rail(st.sidebar, scheme, specs,
                     skip=[sp for sp in specs if sp.name in grouped])"""
        if rail else "params = render_conditions(scheme, specs)"
    )
    return AppTest.from_string(f"""
import streamlit as st
from gabes import schemes
from gabes_ui.conditions import atom_condition_descriptor, render_atom_conditions, render_comparison_controls, render_conditions
from gabes_ui.controls import render_rail
scheme = schemes.get({scheme_name!r})
specs = scheme.param_schema()
for key, value in scheme.defaults().items():
    st.session_state.setdefault(scheme.name+'__'+key, value)
if not st.session_state.get('_initialized'):
    for key, value in {updates or {}!r}.items():
        st.session_state[scheme.name+'__'+key] = value
    st.session_state['_initialized'] = True
{render}
render_comparison_controls(scheme, specs)
descriptor = atom_condition_descriptor(scheme, specs, params)
st.session_state['_grouped'] = sorted(descriptor['names'])
st.session_state['_visible_grouped'] = [sp.name for sp in descriptor['fields']]
st.session_state['_atom_menus'] = [sp.name for sp in descriptor['fields'] if sp.atom_menu]
""").run()


def _atom_popup(app):
    popups = [popup for popup in app.get("popover")
              if popup.proto.popover.label.startswith("Atom")]
    assert len(popups) == 1
    return popups[0]


def _widgets(app):
    return [widget for kind in ("selectbox", "text_input", "checkbox", "slider")
            for widget in app.get(kind)]


@pytest.mark.parametrize("scheme_name,updates", SCHEME_CASES)
def test_atom_group_contains_supported_controls_once(scheme_name, updates):
    app = _app(scheme_name, updates)
    assert not app.exception and not app.warning and not app.error
    popup = _atom_popup(app)
    # Grouping also owns conditional and hidden fields, so a visibility change
    # cannot expose a second copy in a standalone or Advanced popover.
    visible = set(app.session_state['_visible_grouped'])
    assert visible <= set(app.session_state['_grouped'])
    assert {'temp_c', 'cell_mm'} <= set(app.session_state['_grouped'])
    atom_menus = set(app.session_state['_atom_menus'])
    all_widgets = _widgets(app)
    popup_keys = {widget.key for widget in _widgets(popup)}
    for name in visible - atom_menus:
        key = f"conditions-{scheme_name}-{name}"
        matches = [widget for widget in all_widgets
                   if widget.key == key or widget.key.startswith(key+'-refresh-')]
        assert len(matches) == 1, name
        assert matches[0].key in popup_keys, name
    # Canonical grouping replaces the old separate Atom and Isotope triggers.
    assert not any(p.proto.popover.label.startswith('Isotope') for p in app.get('popover'))
    labels = Counter(widget.label for widget in all_widgets)
    assert labels['Atom'] <= 1 and labels['Isotope'] <= 1


def test_compact_atom_summary_omits_custom_fraction_and_numeric_cell_details():
    app = _app('sas', dict(rb_composition='Custom', rb85_fraction=.4, line='D2',
                           ne_pressure_torr=10.5, temp_c=81., cell_mm=43.5))
    popup = _atom_popup(app)
    assert popup.proto.popover.label == 'Atom \u00b7 Custom Rb D2 Buffer'
    assert not any(value in popup.proto.popover.label for value in ('0.4', '10.5', '81', '43.5', 'Ne'))
    assert any(widget.label == '\u2078\u2075Rb fraction (Natural 0.7217)'
               for widget in popup.text_input)


def test_atom_popup_temperature_variants_remain_conditional_direct_inputs():
    linked = _app('rydberg_eit')
    separated = _app('rydberg_eit', {'temperature_model': 'Separated'})
    assert not linked.exception and not separated.exception
    linked_labels = {widget.label for widget in _atom_popup(linked).text_input}
    separated_labels = {widget.label for widget in _atom_popup(separated).text_input}
    assert 'Linked cell temperature (\u00b0C)' in linked_labels
    assert 'Linked cell temperature (\u00b0C)' not in separated_labels
    assert {'Heater setpoint (\u00b0C)', 'Effective vapor temperature (\u00b0C)',
            'Cold-spot temperature (\u00b0C)'} <= separated_labels
    assert not linked.slider and not separated.slider


@pytest.mark.parametrize('scheme_name,updates', [
    ('rydberg_eit', {}),
    ('fwm', {}),
    ('psr', {'calculation': 'Rb87 D1 steady state'}),
])
def test_fixed_identity_has_no_invented_atom_isotope_or_line_options(scheme_name, updates):
    app = _app(scheme_name, updates)
    assert not app.exception
    popup = _atom_popup(app)
    assert not any(widget.label in ('Atom', 'Isotope', 'Transition line', 'Cell type')
                   for widget in popup.selectbox)


@pytest.mark.parametrize('updates', [{}, {'mode': 'Biphoton'}])
def test_atom_helper_excludes_grouped_fields_from_legacy_fwm_rail(updates):
    app = _app('fwm', updates, rail=True)
    assert not app.exception
    popup = _atom_popup(app)
    for label in ('Temperature (\u00b0C)', 'Cell length (mm)'):
        assert sum(widget.label == label for widget in popup.text_input) == 1
    grouped = set(app.session_state['_grouped'])
    assert not any(widget.key.startswith('fwm__') and widget.key[5:] in grouped
                   for widget in _widgets(app))


def test_sas_cell_type_dropdown_keeps_native_coating_parameter():
    app = _app('sas')
    cell = next(widget for widget in _atom_popup(app).selectbox if widget.label == 'Cell type')
    cell.select('Paraffin coated').run()
    assert not app.exception
    assert app.session_state['sas__paraffin_coated'] is True
    assert 'Paraffin' in _atom_popup(app).proto.popover.label
    next(widget for widget in _atom_popup(app).selectbox if widget.label == 'Cell type').select('Uncoated').run()
    assert app.session_state['sas__paraffin_coated'] is False


def test_atom_number_defaults_invalid_drafts_and_preset_refresh_are_preserved():
    app = AppTest.from_string("""
import streamlit as st
from gabes import schemes
from gabes_ui.conditions import render_conditions
scheme = schemes.get('sas')
for key, value in scheme.defaults().items():
    st.session_state.setdefault('sas__'+key, value)
st.button('Load saved', on_click=lambda: st.session_state.update(sas__temp_c=70., sas__cell_mm=25.))
render_conditions(scheme, scheme.param_schema())
""").run()

    def temperature():
        return next(widget for widget in _atom_popup(app).text_input
                    if widget.label == 'Cell temperature (\u00b0C)')

    def length():
        return next(widget for widget in _atom_popup(app).text_input
                    if widget.label == 'Cell length (mm)')

    assert temperature().proto.default == '40' and length().proto.default == '75'
    identity = temperature().proto.id
    temperature().set_value('81.2345').run()
    temperature().set_value('-2').run()
    assert app.error and temperature().value == '-2'
    assert temperature().proto.default == '-2' and temperature().proto.id == identity
    assert app.session_state['sas__temp_c'] == 81.2345
    length().set_value('4 cm').run()
    assert temperature().value == '-2'
    assert app.session_state['sas__cell_mm'] == 40.
    app.button[0].click().run()
    assert not app.exception and not app.warning and not app.error
    assert temperature().value == temperature().proto.default == '70'
    assert length().value == length().proto.default == '25'


def test_comparison_editors_render_once_with_explicit_native_defaults():
    app = AppTest.from_string("""
import streamlit as st
from gabes.schemes.base import ParamSpec
from gabes_ui.conditions import render_conditions, render_comparison_controls
class Example:
    name = 'example'
    workspace_inputs = ()
    model_comparison_controls = ('axis', 'signal', 'count')
    def normalize_params(self, params):
        return params
    def param_schema(self):
        return [ParamSpec('axis', 'Scan axis', 'Readout', 'Field', choices=('Optical', 'Field')),
                ParamSpec('signal', 'Signal', 'Readout', 'Rotation', choices=('Transmission', 'Rotation'), advanced=True),
                ParamSpec('count', 'Sample count', 'Numerics', 21, 3, 101, 2)]
    def atom_conditions(self, params):
        return {}
scheme = Example()
for sp in scheme.param_schema():
    st.session_state.setdefault('example__'+sp.name, sp.default)
render_conditions(scheme, scheme.param_schema())
with st.container(key='comparison_controls'):
    render_comparison_controls(scheme, scheme.param_schema())
""").run()
    assert not app.exception and not app.warning
    assert not app.get('popover')
    assert len(app.selectbox) == 2 and len(app.text_input) == 1
    for widget in app.selectbox:
        assert widget.proto.default == widget.options.index(widget.value)
    numeric = app.text_input[0]
    assert numeric.proto.default == '21'
    numeric.set_value('33').run()
    numeric = app.text_input[0]
    assert app.session_state['example__count'] == 33 and numeric.proto.default == '33'


def test_descriptor_uses_declared_roles_and_excludes_hidden_details_from_summary(monkeypatch):
    from gabes.schemes.base import ParamSpec
    from gabes_ui import conditions

    class Example:
        name = 'example'

        def atom_conditions(self, params):
            return dict(fields=('medium', 'optical_line', 'length', 'hidden_detail'),
                        line='optical_line', cell='Buffer')

    specs = [ParamSpec('medium', 'Example medium', 'Atomic', 'stored',
                       atom_menu=(('Example', 'Custom', 'stored', {}),)),
             ParamSpec('optical_line', 'Transition line', 'Atomic', 'L2', choices=('L1', 'L2')),
             ParamSpec('length', 'Cell length', 'Cell', 25.),
             ParamSpec('hidden_detail', 'Hidden detail', 'Cell', 13., hidden=True)]
    monkeypatch.setattr(conditions, 'param_visible', lambda scheme, spec: not spec.hidden)
    descriptor = conditions.atom_condition_descriptor(
        Example(), specs, {'medium': 'stored', 'optical_line': 'L2', 'length': 25., 'hidden_detail': 13.})
    assert descriptor['title'] == 'Atom \u00b7 Custom Example L2 Buffer'
    assert descriptor['names'] == {'medium', 'optical_line', 'length', 'hidden_detail'}
    assert [spec.name for spec in descriptor['fields']] == ['medium', 'optical_line', 'length']
