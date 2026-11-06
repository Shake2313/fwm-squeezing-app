"""Stage-3 Rydberg views preserve native grids, conventions and IF scope."""
from copy import deepcopy
import json

import numpy as np
import pytest

from gabes import schemes
from gabes.core import blas_single_thread
from gabes.model_comparison import comparison_values
from gabes.schemes.rydberg_comparison import comparison_key
from gabes_ui.spectral_compare import comparison_plots, workspace_selection
from gabes_ui.spectral_sweep import resolve_row


def compare(mode="EIT",choice="Probe power",single=False,**changes):
    scheme=schemes.get("rydberg_eit")
    selector="eit_compare" if mode=="EIT" else "at_compare"
    params={**scheme.defaults(),**scheme.recommended_defaults({})[mode],selector:choice,**changes}
    key=comparison_key(params)
    raws={};views={};calls=[]
    def resolve(p):
        calls.append(p[key]);raw=scheme.compute(p)
        raws[p[key]]=raw;views[p[key]]=scheme.headless_observables(raw,p)
        return raw
    with blas_single_thread():
        raw=resolve(params);calls.clear()
        if single:
            with comparison_values(key,[params[key]]):
                data=scheme.model_comparison(params,resolve_raw=resolve,current_raw=raw,current_view=views[params[key]])
        else:
            data=scheme.model_comparison(params,resolve_raw=resolve,current_raw=raw,current_view=views[params[key]])
    return scheme,params,data,raws,views,calls


@pytest.mark.parametrize("mode,choice",[("EIT","Probe power"),("EIT","Coupling power"),("EIT","Temperature"),
    ("AT electrometry","RF Rabi"),("AT electrometry","RF detuning"),("AT electrometry","Probe power"),
    ("AT electrometry","Coupling power"),("AT electrometry","Temperature"),("AT electrometry","IF frequency")])
def test_comparisons_keep_actual_control_and_optical_grids(mode,choice):
    scheme,params,data,raws,views,calls=compare(mode,choice,scan_axis="Coupling detuning",probe_detuning_mhz=2.)
    key=comparison_key(params)
    assert params[key] not in calls
    assert data["axis"]["reference"].endswith("-2 MHz")
    assert data["varied"]["key"]==key and data["display"]["default_mode"]=="map"
    for row in data["series"]:
        value=row["parameter"]["value"];source=views[value]["rydberg_spectrum"]
        np.testing.assert_array_equal(row["x"],source["x"])
        np.testing.assert_array_equal(row["observables"]["transmission"],source["transmission"])
        assert len(row["x"])==801
        if mode!="EIT" and source["finite_if"] is not None:
            sh=source["finite_if"]
            for quantity in ("if_response","if_phase","total_sensitivity","psn_sensitivity"):
                np.testing.assert_array_equal(row["observable_x"][quantity],sh["response_detuning_mhz"])
            np.testing.assert_array_equal(row["observables"]["if_response"],np.abs(sh["transmission_phasor_per_angular_rabi"])*2*np.pi*1e3)
            np.testing.assert_array_equal(row["observables"]["total_sensitivity"],sh["total_field_asd_v_m_per_sqrt_hz"]*1e7)
            np.testing.assert_array_equal(row["observables"]["psn_sensitivity"],sh["psn_field_asd_v_m_per_sqrt_hz"]*1e7)
            np.testing.assert_array_equal(row["observables"]["if_phase"],np.angle(sh["current_phasors_a_per_v_m"],deg=True))
            assert row["optimum"]["x"]==sh["optimum_detuning_mhz"]
    json.dumps(data,allow_nan=False)
    assert not {"eit_compare","at_compare","temperature_compare"}&set(scheme.recompute_keys())


def test_eit_features_and_optional_residual_are_the_existing_power_analysis():
    scheme,params,data,raws,_,_=compare()
    for row in data["series"]:
        value=row["parameter"]["value"];raw=raws[value]
        for variant,chi_key in (("comp","chi_bar"),("residual","chi_bar_uncomp")):
            x,t,_=scheme._transmission(raw[chi_key],raw,{**params,"probe_power_uw":value})
            width,contrast=scheme._eit_features(x,t)
            assert row["analysis"][variant]["width"]==width
            assert row["analysis"][variant]["contrast"]==contrast
            assert row["analysis"][variant]["amplitude"]==value*contrast
        np.testing.assert_array_equal(row["overlays"]["residual"]["y"],scheme._transmission(raw["chi_bar_uncomp"],raw,params)[1])


@pytest.mark.parametrize("temperature_compare,key",[("Effective vapor temperature","effective_temp_c"),("Cold-spot temperature","cold_spot_temp_c")])
def test_separated_temperature_scan_does_not_replace_heater_or_other_temperature(temperature_compare,key):
    _,params,data,raws,_,_=compare("EIT","Temperature",temperature_model="Separated",
        heater_setpoint_c=60.,effective_temp_c=40.,cold_spot_temp_c=25.,temperature_compare=temperature_compare)
    assert data["varied"]["key"]==key
    for row in data["series"]:
        value=row["parameter"]["value"];raw=raws[value]
        assert raw["heater_setpoint_c"]==60.
        assert raw["effective_temp_c"]==(value if key=="effective_temp_c" else 40.)
        assert raw["cold_spot_temp_c"]==(value if key=="cold_spot_temp_c" else 25.)
        assert row["temperature"]["model"]=="Separated"


def test_finite_if_grid_and_export_are_sparse_with_doppler_not_interpolated():
    _,_,data,_,views,_=compare("AT electrometry","RF Rabi",single=True,doppler="on")
    row=data["series"][0];native=next(iter(views.values()))["rydberg_spectrum"]["finite_if"]
    assert 0<len(row["observable_x"]["total_sensitivity"])<801
    assert "static-slope" in row["finite_if_scope"]
    assert all(p["sampled_only"] for p in data["panels"] if p["key"]!="transmission")
    trace=next(p for p in comparison_plots(data) if p["label"]=="Total noise-equivalent field")["traces"][0]
    np.testing.assert_array_equal(trace["x_data"],native["response_detuning_mhz"])
    assert len(trace["x_data"])==len(trace["y_data"])<801


def test_zero_opacity_rf_off_and_missing_scenarios_are_not_zero_sensitivity():
    _,_,zero,_,_,_=compare("AT electrometry","RF Rabi",single=True,effective_opacity_participation=0.)
    row=zero["series"][0]
    assert row["optimum"] is None and row["observables"]["total_sensitivity"]==[]
    assert next(r["value"] for r in row["readouts"] if r["key"]=="best_total") is None
    _,_,edge,_,_,_=compare("AT electrometry","RF Rabi",single=True,effective_opacity_relative_uncertainty=1.)
    lower=next(s for s in edge["series"][0]["opacity_scenarios"] if s["label"]=="lower")
    assert lower["value"]==0 and lower["analysis"]["total"] is None
    assert next(s for s in edge["summary"]["series"] if s["id"]=="total-lower")["observables"]["value"]==[None]
    _,_,no_residual,_,_,_=compare("EIT","Probe power",single=True,residual_zeeman_mhz=0.)
    assert no_residual["series"][0]["overlays"]=={}
    assert not any(o["key"]=="residual" for o in no_residual["optional_overlays"])


def test_profiles_and_overlay_selection_validate_and_leave_native_arrays_intact():
    _,params,data,_,_,_=compare("AT electrometry","IF frequency",single=True)
    before=deepcopy(data)
    selected=workspace_selection(data,None,"ctx",params)
    assert selected["profile"]=="sensitivity" and selected["panel_key"]=="total_sensitivity"
    selected.update(profile="unknown",overlays=["opacity","injected"])
    valid=workspace_selection(data,selected,"ctx",params)
    assert valid["profile"]=="spectra" and valid["panel_key"]=="transmission"
    assert valid["overlays"]==["opacity"]
    assert data==before


def test_incremental_row_preserves_new_axis_summaries_and_individual_response_grids():
    scheme=schemes.get("rydberg_eit")
    for choice,key,value in [("IF frequency","if_khz",100.),("Temperature","temp_c",32.),("Probe power","probe_power_uw",7.)]:
        params={**scheme.defaults(),**scheme.recommended_defaults({})["AT electrometry"],"at_compare":choice}
        calls=[]
        def resolve(p): calls.append(p[key]);return scheme.compute(p)
        with blas_single_thread(): payload=resolve_row(scheme,params,key,value,resolve)
        assert calls==[value]
        assert payload["series"]["parameter"]["value"]==value
        assert len(payload["summary"]["series"][0]["x"])==1
        assert all(len(y)==len(payload["series"]["observable_x"].get(name,payload["series"]["x"]))
                   for name,y in payload["series"]["observables"].items())


def test_dynamic_temperature_and_if_editors_share_canonical_parameters_once():
    from streamlit.testing.v1 import AppTest
    app=AppTest.from_string('''
import streamlit as st
from gabes import schemes
from gabes_ui.conditions import render_conditions,render_comparison_controls,render_power_inputs
scheme=schemes.get("rydberg_eit")
for name,value in scheme.defaults().items(): st.session_state.setdefault("rydberg_eit__"+name,value)
st.session_state["rydberg_eit__temperature_model"]="Separated"
st.button("AT",on_click=lambda:st.session_state.update(rydberg_eit__view="AT electrometry"))
render_conditions(scheme,scheme.param_schema())
render_comparison_controls(scheme,scheme.param_schema())
render_power_inputs(scheme,scheme.param_schema())
''').run()
    assert not any(w.label=="Temperature scan" for w in app.selectbox)
    next(w for w in app.selectbox if w.label=="Compare").select("Temperature").run()
    next(w for w in app.selectbox if w.label=="Temperature scan").select("Cold-spot temperature").run()
    fields=[w for w in app.text_input if w.label=="Cold-spot temperature (°C)"]
    assert len(fields)==1
    fields[0].set_value("25").run()
    assert app.session_state["rydberg_eit__cold_spot_temp_c"]==25.
    app.button[0].click().run()
    next(w for w in app.selectbox if w.label=="Compare").select("IF frequency").run()
    fields=[w for w in app.text_input if w.label=="IF offset (kHz)"]
    assert len(fields)==1
    fields[0].set_value("0.1 MHz").run()
    assert app.session_state["rydberg_eit__if_khz"]==100.
    assert not app.exception and not app.error and not app.slider
