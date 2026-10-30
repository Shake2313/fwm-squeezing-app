"""Measured detector-reference interpolation never replaces atomic gain or NF."""
from copy import deepcopy
import json

import matplotlib
import numpy as np
import pytest
from scipy.interpolate import PchipInterpolator

matplotlib.use("Agg")

from gabes import fwm_measured_gain as measured
from gabes.schemes import fwm


def reference_params():
    scheme = fwm.FWMScheme()
    params = dict(scheme.defaults(), **scheme.presets()[0].values)
    params["noise_model"] = fwm.NOISE_PASSIVE_L1
    params["experimental_gain_calibration"] = False
    params["measured_gain_calibration"] = True
    return params


def reference_raw(params=None):
    params = reference_params() if params is None else params
    frequency = measured.observations()["eom_MHz"]
    center = fwm.branch_center_GHz(params["opd"], -1)
    delta = measured.internal_delta_mhz_from_eom(frequency)
    matrix = np.array([[2.0, np.sqrt(3.0)], [np.sqrt(3.0), 2.0]], complex)
    return {
        "D_GHz": params["opd"],
        "probe_axis_GHz": center+delta*1e-3,
        "raman_center_minus_GHz": center,
        "raman_center_plus_GHz": fwm.branch_center_GHz(params["opd"], 1),
        "temperature_K": params["temp_c"]+273.15,
        "input_pump_power_W": params["pump_mw"]*1e-3,
        "input_seed_power_W": params["probe_uw"]*1e-6,
        "cell_length_m": params["cell_mm"]*1e-3,
        "w_pump_m": params["pump_waist_um"]*1e-6,
        "w_probe_m": params["probe_waist_um"]*1e-6,
        "pump_probe_angle_deg": params["seeded_angle_deg"], "branch": -1,
        "legacy_mean_field_fixture": False,
        "gain_closure": {"historical_fixture_requested": False},
        "T_canonical_small_signal": np.tile(matrix, (frequency.size, 1, 1)),
        "G_s": 100.0+np.arange(frequency.size),
        "G_c": 90.0+np.arange(frequency.size),
        "gain_referred_noise_dB": np.full(frequency.size, -5.0),
        "S_dB": np.full(frequency.size, -5.0), "eta": 1.0,
        "N_atoms": 1e18, "sigma_v": 200.0, "n_velocity": 3,
        "Op_A_2pi_GHz": 0.5, "Os_2pi_MHz": 0.1,
        "phase_detail": fwm.PHASE_LEGACY,
        "physical_squeezing_dB": None,
        "claim_gate": {"physical_squeezing_prediction": False},
    }


def test_preset_explicit_proxy_and_default_off_preserve_full_view_identity():
    scheme = fwm.FWMScheme()
    specs = {sp.name: sp for sp in scheme.param_schema()}
    assert scheme.defaults()["measured_gain_calibration"] is False
    assert specs["measured_gain_calibration"].recompute is False
    assert specs["measured_gain_calibration"].hidden is True
    assert specs["seed_input_assumption"].recompute is False
    assert {"measured_gain_calibration", "seed_input_assumption"}.isdisjoint(scheme.recompute_keys())
    preset = scheme.presets()[0]
    assert preset.values["measured_gain_calibration"] is False
    assert preset.values["seed_input_assumption"] == measured.SEED_INPUT_ASSUMPTION
    assert preset.values["probe_uw"]*1e-6 == pytest.approx(
        (35*20/228.5)*1e-6/(0.8046*0.9709), rel=1e-14)
    assert preset.values["tpd"] == pytest.approx(1.732439, abs=1e-9)
    assert scheme.extra_views()[0].key == "Full probe scan (−8 to 12 GHz, slow)"
    assert [p.name for p in scheme.presets()] == [
        "2026-10-13 experiment", "Sim Gold Reference"]
    assert all(view.key != "2026-10-13 measured gain comparison"
               for view in scheme.extra_views())
    assert scheme.internal_measured_gain_comparison().key == "2026-10-13 measured gain comparison"
    assert specs["experimental_gain_calibration"].label == "Empirical gain calibration"
    gold = scheme.presets()[1].values
    assert set(gold) == set(scheme.defaults())
    assert set(preset.values) == set(gold)
    assert gold["mode"] == fwm.MODE_SEEDED
    assert gold["experimental_gain_calibration"] is False
    assert gold["measured_gain_calibration"] is False
    assert gold["seed_input_assumption"] == ""
    assert gold["opd"] == 0.9 and gold["temp_c"] == 121.0
    assert gold["pump_mw"] == 600.0 and gold["probe_uw"] == 8.0


def test_reference_selector_and_proxy_adoption_do_not_enter_cached_atomic_calls(monkeypatch):
    scheme = fwm.FWMScheme()
    params = reference_params()
    baseline = dict(params, measured_gain_calibration=False, seed_input_assumption="")
    for name, compute in (("compute_spectrum", scheme.compute),
                          ("full_spectrum", scheme.extra_views()[0].compute)):
        calls = []
        monkeypatch.setattr(fwm, name, lambda *a, **kw: calls.append((a, kw)) or {})
        compute(baseline)
        compute(params)
        assert calls[0] == calls[1]


def test_focused_headless_both_arms_use_exact_selected_rf_and_keep_atomic_nf():
    import matplotlib.pyplot as plt
    scheme = fwm.FWMScheme()
    params = reference_params()
    params["tpd"] = float(measured.internal_delta_mhz_from_eom(3037.25))
    raw = reference_raw(params)
    original = deepcopy(raw)
    baseline = scheme.headless_observables(raw, dict(params, measured_gain_calibration=False))
    headless = scheme.headless_observables(raw, params)
    view = scheme.observables(raw, params)
    try:
        point = view["measured_gain_reference"]["selected_point"]
        assert point["available"]
        assert point["eom_MHz"] == pytest.approx(3037.25, abs=1e-10)
        observed = measured.observations()
        order = np.argsort(observed["eom_MHz"])
        for arm in ("probe", "conjugate"):
            expected = np.exp(PchipInterpolator(
                observed["eom_MHz"][order], np.log(observed[f"{arm}_ratio"][order]),
                extrapolate=False)(3037.25))
            assert point[f"measured_gain_{arm}"] == pytest.approx(expected, rel=1e-13)
        assert headless["metrics"] == view["metrics"]
        np.testing.assert_array_equal(baseline["jasperse_readout"]["NF"],
                                      view["jasperse_readout"]["NF"])
        assert view["passive_l1_operating_point"]["NF"] == baseline["passive_l1_operating_point"]["NF"]
        labels = [item["label"] for item in view["metrics"]]
        assert "Probe gain, pump-off reference" in labels
        assert "Conjugate gain, pump-off reference" in labels
        assert "Atomic probe gain G_s" in labels
        assert "Atomic conjugate gain G_c" in labels
        assert "does not recalibrate squeezing" in view["metrics"][0]["delta"]
        np.testing.assert_array_equal(view["figure"].axes[0].lines[0].get_ydata(), raw["G_s"])
        table = view["tables"][0]["markdown"]
        assert "pump-off seed detector output" in table
        assert "vapor input not independently measured" in table
        assert "By construction; not independent validation" in table
    finally:
        plt.close(view["figure"])
    for key in ("G_s", "G_c", "T_canonical_small_signal", "S_dB"):
        np.testing.assert_array_equal(raw[key], original[key])


@pytest.mark.parametrize("change", ("temperature", "pump", "seed_proxy", "legacy", "adoption"))
def test_context_mismatches_fail_closed_without_reference_gain_replacement(change):
    params = reference_params()
    raw = reference_raw(params)
    if change == "temperature": raw["temperature_K"] += 1.0
    if change == "pump": raw["input_pump_power_W"] *= 1.01
    if change == "seed_proxy": raw["input_seed_power_W"] *= 1.01
    if change == "legacy": raw["legacy_mean_field_fixture"] = True
    if change == "adoption": params["seed_input_assumption"] = ""
    view = fwm.FWMScheme().headless_observables(raw, params)
    reference = view["measured_gain_reference"]
    assert not reference["applied"]
    assert reference["status"] == "INELIGIBLE_CONTEXT"
    assert not reference["selected_point"]["available"]
    assert np.all(np.isnan(reference["measured_gain_probe"]))
    assert np.all(np.isnan(reference["measured_gain_conjugate"]))
    assert any(item["label"] == "Measured gain reference" and item["value"] == "unavailable"
               for item in view["metrics"])


def test_outside_rf_support_is_unavailable_and_detector_noise_controls_do_not_refit():
    scheme = fwm.FWMScheme()
    params = reference_params()
    raw = reference_raw(params)
    baseline = scheme.headless_observables(raw, params)
    changed = scheme.headless_observables(raw, dict(params,
        after_cell_eff_pct=50.0, qe_pct=60.0, seed_classical_rin_scaled=4.0))
    for arm in ("probe", "conjugate"):
        np.testing.assert_array_equal(baseline["measured_gain_reference"][f"measured_gain_{arm}"],
                                      changed["measured_gain_reference"][f"measured_gain_{arm}"])
    assert not np.array_equal(baseline["jasperse_readout"]["NF"], changed["jasperse_readout"]["NF"])
    outside = scheme.headless_observables(raw, dict(params,
        tpd=float(measured.internal_delta_mhz_from_eom(3044.01))))
    point = outside["measured_gain_reference"]["selected_point"]
    assert not point["available"] and point["status"] == "OUTSIDE_RF_SUPPORT"
    assert np.isnan(point["measured_gain_probe"])


def test_all_23_axis_roundtrips_preserve_endpoint_nodes_and_scalar_array_agreement():
    raw = reference_raw()
    # This construction reproduces the upper endpoint's floating cancellation.
    raw["probe_axis_GHz"] = raw["D_GHz"]-measured.observations()["eom_MHz"]/1000
    reference = fwm._measured_gain_reference_readout(
        raw, enabled=True, seed_input_assumption=measured.SEED_INPUT_ASSUMPTION)
    assert np.all(reference["within_support"])
    for arm in ("probe", "conjugate"):
        np.testing.assert_allclose(reference[f"measured_gain_{arm}"],
                                   measured.observations()[f"{arm}_ratio"], rtol=1e-12)
    for index in (0, -1):
        point = fwm._measured_gain_reference_readout(
            raw, enabled=True, seed_input_assumption=measured.SEED_INPUT_ASSUMPTION,
            selected_delta_mhz=reference["internal_delta_mhz"][index])["selected_point"]
        assert point["available"]
        for arm in ("probe", "conjugate"):
            assert point[f"measured_gain_{arm}"] == reference[f"measured_gain_{arm}"][index]
    # A finite physical offset beyond representational roundoff still fails closed.
    outside = fwm._measured_gain_reference_readout(
        raw, enabled=True, seed_input_assumption=measured.SEED_INPUT_ASSUMPTION,
        selected_delta_mhz=measured.internal_delta_mhz_from_eom(3044+1e-8))
    assert not outside["selected_point"]["available"]


def test_plus_branch_operating_point_cannot_reuse_minus_seed_reference():
    raw = reference_raw()
    raw["measured_gain_reference"] = fwm._measured_gain_reference_readout(
        raw, enabled=True, seed_input_assumption=measured.SEED_INPUT_ASSUMPTION)
    point = fwm.operating_point(raw, measured.internal_delta_mhz_from_eom(3034), branch=1)
    assert not point["measured_gain_reference"]["available"]
    assert point["measured_gain_reference"]["status"] == "INELIGIBLE_CONTEXT"
    assert np.isnan(point["measured_gain_reference"]["measured_gain_probe"])


def test_direct_api_optional_reference_leaves_atomic_solution_and_noise_exact():
    context = measured.nominal_context()
    center = fwm.branch_center_GHz(context["D_GHz"], -1)
    kwargs = dict(context)
    d = kwargs.pop("D_GHz")
    kwargs.update(coarse_points=5, fine_points=0, scan_min=center-0.015,
                  scan_max=center+0.015, velocity_step=200.0, velocity_cutoff=0.1,
                  phase_detail=fwm.PHASE_BALANCED)
    off = fwm.compute_spectrum(d, **kwargs)
    on = fwm.compute_spectrum(d, measured_gain_calibration=True, **kwargs)
    for key in ("G_s", "G_c", "G_s_smallsignal", "G_c_smallsignal",
                "T_canonical_small_signal", "gain_referred_noise_dB", "NF", "NF_dB"):
        np.testing.assert_array_equal(off[key], on[key])
    point = fwm.operating_point(on, float(measured.internal_delta_mhz_from_eom(3034)))
    assert point["measured_gain_reference"]["measured_gain_probe"] == pytest.approx(870/35)
    assert point["measured_gain_reference"]["measured_gain_conjugate"] == pytest.approx(850/35)
    assert not off["measured_gain_reference"]["requested"]
    assert on["measured_gain_reference"]["profile"]["NF_modified"] is False


def test_full_api_forwards_reference_flag_but_plus_branch_fails_context(monkeypatch):
    context = measured.nominal_context()
    params = reference_params()
    calls = []
    def solve(d, **kwargs):
        calls.append(kwargs)
        raw = reference_raw(params)
        raw["branch"] = kwargs["branch"]
        raw["measured_gain_reference"] = fwm._measured_gain_reference_readout(
            raw, enabled=kwargs["measured_gain_calibration"],
            seed_input_assumption=kwargs["seed_input_assumption"])
        return raw
    monkeypatch.setattr(fwm, "compute_spectrum", solve)
    full = fwm.full_spectrum(1.0, 386.15, 294.804, context["P_probe"]*1e6, 0.74, 0.0,
        measured_gain_calibration=True, seed_input_assumption=measured.SEED_INPUT_ASSUMPTION)
    assert [call["branch"] for call in calls] == [-1, 1]
    assert full["minus"]["measured_gain_reference"]["applied"]
    assert not full["plus"]["measured_gain_reference"]["applied"]


def test_comparison_view_has_only_23_measured_ratios_and_reference_interpolation():
    import matplotlib.pyplot as plt
    extra = fwm.FWMScheme().internal_measured_gain_comparison()
    data = extra.compute(reference_params())
    figure = extra.render(data, reference_params())
    try:
        axis = figure.axes[0]
        assert axis.get_xlim() == (3022.0, 3044.0)
        assert len(axis.lines) == 2 and len(axis.collections) == 2
        assert all(len(scatter.get_offsets()) == 23 for scatter in axis.collections)
        assert "by construction" in axis.get_title()
        assert not data["profile"]["independent_validation"]
        assert data["profile"]["tuned_parameter_count"] == 0
        assert "G_s" not in data and "NF" not in data
    finally:
        plt.close(figure)


def test_strict_json_exports_current_reference_profile_status_domain_and_raw_arrays():
    from gabes_ui.export import build_export_payload
    scheme = fwm.FWMScheme()
    params = reference_params()
    raw = reference_raw(params)
    raw["measured_gain_reference"] = fwm._measured_gain_reference_readout(
        raw, enabled=True, seed_input_assumption=measured.SEED_INPUT_ASSUMPTION)
    view = scheme.headless_observables(raw, params)
    payload = build_export_payload(scheme=scheme, raw=raw, params=params, view=view,
        plot_payloads=[], comparison_payload=None, cache_version=scheme.cache_version,
        readout_cache_version="test")
    decoded = json.loads(json.dumps(payload, allow_nan=False))
    for reference in (decoded["raw"]["measured_gain_reference"],
                      decoded["view"]["measured_gain_reference"]):
        assert reference["profile"]["source_workbook_sha256"] == measured.SOURCE_WORKBOOK_SHA256
        assert reference["profile"]["eom_support_MHz"] == [3022.0, 3044.0]
        assert reference["profile"]["independent_validation"] is False
        assert len(reference["within_support"]) == 23
        assert len(reference["measured_gain_probe"]) == 23
        assert len(reference["measured_gain_conjugate"]) == 23
        assert reference["status"]


def test_actual_streamlit_preset_calibration_and_hidden_reference():
    from streamlit.testing.v1 import AppTest
    app = AppTest.from_string(
        "import streamlit as st\n"
        "from gabes.schemes.fwm import FWMScheme\n"
        "from gabes_ui.controls import render_rail, render_derived_params, skey\n"
        "from gabes_ui.defaults import render_defaults\n"
        "from tests.test_fwm_measured_gain_integration import reference_raw\n"
        "scheme=FWMScheme()\n"
        "st.session_state[skey(scheme.name,'_show_advanced')]=True\n"
        "for key,value in scheme.defaults().items():\n"
        "    sk=skey(scheme.name,key)\n"
        "    st.session_state[sk]=st.session_state.get(sk,value)\n"
        "render_defaults(st,scheme)\n"
        "slots={}\n"
        # Streamlit AppTest cannot serialize formatted segmented controls on
        # rerun. Keep mode/readout in state and exercise
        # the real preset, new checkbox and all remaining generic rail controls.
        "specs=scheme.param_schema()\n"
        "params=render_rail(st,scheme,specs,skip=[sp for sp in specs if sp.name in ('mode','squeezing_view')],derived_slots=slots)\n"
        "view=scheme.headless_observables(reference_raw(params),params)\n"
        "render_derived_params(slots,view['sidebar_readouts'])\n"
        "for metric in view['metrics']:\n"
        "    st.metric(metric['label'],metric['value'],metric.get('delta'))\n"
        "st.markdown(view['tables'][0]['markdown'])\n"
        "st.session_state['view']=view\n",
        default_timeout=30,
    ).run()
    assert not app.exception
    assert "fwm__measured_gain_calibration" not in [item.key for item in app.checkbox]
    assert app.checkbox(key="fwm__experimental_gain_calibration").label == "Empirical gain calibration"
    app.selectbox(key="fwm___named_default").select("builtin:0").run()
    assert not app.exception
    assert app.session_state["view"]["experimental_gain_calibration"]["applied"]
    app.checkbox(key="fwm__experimental_gain_calibration").uncheck().run()
    app.selectbox(key="fwm__noise_model").select(fwm.NOISE_PASSIVE_L1).run()
    assert not app.exception
    assert "fwm__measured_gain_calibration" not in [item.key for item in app.checkbox]
    assert not app.session_state["view"]["measured_gain_reference"]["requested"]
    assert app.slider(key="fwm__probe_uw").value == pytest.approx(
        measured.nominal_context()["P_probe"]*1e6, rel=1e-14)
    labels = [metric.label for metric in app.metric]
    assert "Seed gain G_s" in labels
    assert "Conjugate gain G_c" in labels
    assert "Squeezing [dB]" in labels
