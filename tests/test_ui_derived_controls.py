"""Output-only sidebar controls stay current without entering model inputs."""
import json
import re

from gabes.schemes.base import ParamSpec, Scheme


class DerivedScheme(Scheme):
    name = "derived_test"

    def param_schema(self):
        return [
            ParamSpec("seed", "Seed power", "Noise", 2.0, 1.0, 10.0, 1.0),
            ParamSpec("a", "Classical seed noise a", "Noise", 0.0,
                      0.0, 1.0, 0.01, recompute=False),
            ParamSpec("seed_shot_noise_b", "Seed shot noise b", "Noise", None,
                      0.0, 100.0, 0.01, control="derived_slider"),
            ParamSpec("efficiency", "After-cell efficiency", "Detection", 90.0,
                      0.0, 100.0, 1.0, "%", recompute=False),
            ParamSpec("available", "Calculated transmission available", "Detection",
                      True, control="checkbox", recompute=False),
            ParamSpec("eta_probe_pct", "Total Probe Efficiency", "Detection", None,
                      0.0, 100.0, 0.01, "%", control="derived_slider"),
            ParamSpec("eta_conjugate_pct", "Total Conjugate Efficiency", "Detection",
                      None, 0.0, 100.0, 0.01, "%", control="derived_slider"),
        ]

    def compute(self, params):
        return {}

    def observables(self, raw, params, include_figures=True):
        return {}


def test_derived_rows_never_become_defaults_or_solve_parameters():
    scheme = DerivedScheme()
    assert set(scheme.defaults()) == {"seed", "a", "efficiency", "available"}
    assert scheme.recompute_keys() == ("seed",)


def test_derived_rows_are_read_only_current_and_omitted_from_inputs():
    from streamlit.testing.v1 import AppTest

    app = AppTest.from_string(
        "from tests.test_ui_derived_controls import DerivedScheme\n"
        "import streamlit as st\n"
        "from gabes_ui.controls import (current_params, render_rail, "
        "render_derived_params)\n"
        "scheme = DerivedScheme()\n"
        "slots = {}\n"
        "params = render_rail(st, scheme, scheme.param_schema(), derived_slots=slots)\n"
        "st.session_state['collected_params'] = params\n"
        "st.session_state['current_params'] = current_params(scheme.name, scheme)\n"
        "eta = params['efficiency'] if params['available'] else None\n"
        "render_derived_params(slots, {\n"
        " 'seed_shot_noise_b': {'value': params['seed'] * 2.0},\n"
        " 'eta_probe_pct': {'value': eta, 'status': '' if eta is not None else "
        "'Unavailable: transmission is not passive'},\n"
        " 'eta_conjugate_pct': {'value': params['efficiency'] * 0.5},\n"
        "})\n"
    ).run()
    assert not app.exception
    assert [slider.label for slider in app.slider] == [
        "Seed power", "Classical seed noise a", "After-cell efficiency  [%]",
    ]
    assert all(not slider.proto.disabled for slider in app.slider)

    def row_values():
        return [float(re.search(r"class='g-derived-value'>([^<]+)", m.value)[1])
                for m in app.markdown if "class='g-derived-value'" in m.value
                and ">Unavailable<" not in m.value]

    assert row_values() == [4.0, 90.0, 45.0]
    expected_keys = set(DerivedScheme().defaults())
    assert set(app.session_state["collected_params"]) == expected_keys
    assert set(app.session_state["current_params"]) == expected_keys

    app.slider(key="derived_test__seed").set_value(3.0)
    app.slider(key="derived_test__efficiency").set_value(80.0)
    app.run()
    assert not app.exception
    assert row_values() == [6.0, 80.0, 40.0]
    assert set(app.session_state["collected_params"]) == expected_keys

    app.checkbox(key="derived_test__available").uncheck().run()
    assert not app.exception
    assert len(app.slider) == 3
    assert any("Unavailable: transmission is not passive" in item.value
               for item in app.caption)
    assert row_values() == [6.0, 40.0]
    assert any(">Unavailable<" in m.value for m in app.markdown)


def test_export_preserves_display_values_and_unavailable_status():
    from gabes_ui.export import build_export_payload

    readouts = {
        "seed_shot_noise_b": {"value": 1.4e-14, "status": "calculated"},
        "eta_probe_pct": {"value": None, "status": "unavailable"},
    }
    payload = build_export_payload(
        scheme=DerivedScheme(), raw={}, params={"seed": 2.0},
        view={"sidebar_readouts": readouts,
              "jasperse_readout": {"NF": [0.6, float("nan")]}},
        plot_payloads=[], comparison_payload=None,
        cache_version="test", readout_cache_version="test",
    )
    decoded = json.loads(json.dumps(payload, allow_nan=False))
    assert decoded["view"]["sidebar_readouts"] == readouts
    assert decoded["view"]["jasperse_readout"]["NF"] == [0.6, None]
    assert "seed_shot_noise_b" not in decoded["defaults"]
