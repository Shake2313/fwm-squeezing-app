"""Shared scan rows preserve native physics, selection and cancellation scope."""
from copy import deepcopy
from threading import Event
from time import monotonic, sleep
from types import SimpleNamespace

import numpy as np
import pytest

from gabes import schemes
from gabes.model_comparison import comparison_values, condition_views
from gabes.schemes.base import ParamSpec
from gabes_ui.spectral_compare import comparison_plots, selected_parameter_change, workspace_selection
from gabes_ui.spectral_sweep import (
    RowSweep, assemble_dataset, extract_row, resolve_row, reusable_dataset,
    cancel_inactive_sweeps, scan_values, workspace_key,
)


def sample_data():
    return dict(schema="gabes-spectral-comparison-v1",scheme="test",mode="test",
        axis=dict(label="Detuning",unit="MHz",reference="Bare"),
        varied=dict(key="power",label="Power",unit="mW",values=[1.,2.]),
        series=[dict(id=str(v),label=f"Power {v}",parameter=dict(key="power",value=v),
                     x=x,observables=dict(T=y),readouts=[]) for v,x,y in
                [(1.,[-2.,0.,2.],[.1,.8,.1]),(2.,[-3.,.2,3.],[.2,None,.2])]],
        panels=[dict(key="t",title="Transmission",y_label="T",unit="1",
                     curves=[dict(series_id=str(v),observable="T") for v in [1.,2.]])],
        summary=dict(axis=dict(label="Power",unit="mW"),
            series=[dict(id="height",label="Height",x=[1.,2.],observables=dict(y=[.8,None]))],
            panels=[dict(key="height",title="Height",y_label="Height",unit="1",
                         curves=[dict(series_id="height",observable="y")])]))


def test_summary_scope_is_sliced_and_assembled_with_its_native_row():
    data=sample_data()
    summary=data["summary"]["series"][0]
    summary.update(sample_status=["available","unresolved"], sample_reasons=[[],["No crossings"]],
                   observable_status={"y":["available","unresolved"]})
    row=extract_row(data,1)
    scope=row["summary"]["series"][0]
    assert scope["sample_status"]==["unresolved"]
    assert scope["sample_reasons"]==[["No crossings"]]
    snapshot=dict(rows=[(2.,row),(1.,extract_row(data,0))],status="complete",requested_values=[2.,1.])
    rebuilt=assemble_dataset(data,snapshot)["summary"]["series"][0]
    assert rebuilt["x"]==[2.,1.]
    assert rebuilt["sample_status"]==["unresolved","available"]
    assert rebuilt["observable_status"]["y"]==["unresolved","available"]


def finish(runner):
    deadline=monotonic()+3
    while monotonic()<deadline:
        snapshot=runner.snapshot()
        if snapshot["status"]!="running":
            return snapshot
        sleep(.005)
    pytest.fail("Sweep did not complete")


def test_uniform_grid_keeps_exact_current_and_validates_limits():
    spec=ParamSpec("power","Power","Test",1.,0.,10.)
    assert scan_values(0.,4.,3,1.,spec)==[0.,1.,2.,4.]
    assert scan_values(2.,4.,3,1.,spec)==[1.,2.,3.,4.]  # explicit anchor, no hidden input change
    for args in [(4.,0.,3,1.),(-1.,2.,3,1.),(0.,11.,3,1.),(0.,2.,2.5,1.),(0.,2.,42,1.),(0.,float("nan"),3,1.)]:
        with pytest.raises(ValueError):
            scan_values(*args,spec)


def test_condition_override_is_scoped_and_restores_after_errors():
    params={"power":1.,"other":2.}
    scheme=SimpleNamespace(headless_observables=lambda raw,p:p)
    resolve=lambda p:p
    with comparison_values("power",[3.,4.]):
        assert [p["power"] for p,_,_ in condition_views(scheme,params,"power",[1.,2.],resolve_raw=resolve)]==[3.,4.]
        assert [p["other"] for p,_,_ in condition_views(scheme,params,"other",[2.,5.],resolve_raw=resolve)]==[2.,5.]
        with pytest.raises(RuntimeError),comparison_values("power",[9.]):
            raise RuntimeError("abort")
        assert [p["power"] for p,_,_ in condition_views(scheme,params,"power",[1.],resolve_raw=resolve)]==[3.,4.]
    assert [p["power"] for p,_,_ in condition_views(scheme,params,"power",[1.],resolve_raw=resolve)]==[1.]


def test_rows_keep_native_arrays_nulls_and_linked_summary_coordinates():
    data=sample_data();before=deepcopy(data)
    rows=[(s["parameter"]["value"],extract_row(data,i)) for i,s in enumerate(data["series"])]
    assembled=assemble_dataset(data,dict(rows=rows,status="cancelled",requested_values=[1.,2.,3.],completed=2,total=3))
    assert data==before
    assert assembled["series"][1]["x"]==[-3.,.2,3.]
    assert assembled["series"][1]["observables"]["T"]==[.2,None,.2]
    assert assembled["summary"]["series"][0]["condition_ids"]==[r["id"] for r in assembled["series"]]
    assert assembled["summary"]["series"][0]["observables"]["y"]==[.8,None]
    assert comparison_plots(assembled)[0]["traces"][1]["x_data"]==[-3.,.2,3.]
    assert assembled["sweep"]["status"]=="cancelled"
    assert 3. not in assembled["varied"]["values"]
    rows[1][1]["axis"]["reference"]="Different physical zero"
    with pytest.raises(ValueError,match="physical axis"):
        assemble_dataset(data,dict(rows=rows))


def test_cache_reuses_completed_rows_and_fixed_inputs_isolate_entries():
    runner=RowSweep(max_entries=3)
    called=[]
    runner.seed("fixed",[(1.,{"native":1.})])
    def resolve(v):
        called.append(v);return {"native":v}
    runner.start("fixed",[1.,2.,3.],resolve)
    first=finish(runner)
    assert called==[2.,3.] and first["cached"]==1 and first["completed"]==3
    first["rows"][0][1]["native"]=999
    runner.start("fixed",[1.,2.,3.],resolve)
    second=finish(runner)
    assert called==[2.,3.] and second["cached"]==3
    assert second["rows"][0][1]["native"]==1.
    runner.start("other temperature",[1.],resolve)
    assert finish(runner)["cached"]==0 and called==[2.,3.,1.]


def test_cancel_stops_next_rows_and_stale_job_cannot_publish_into_new_scan():
    runner=RowSweep();entered=Event();release=Event();called=[]
    def resolve(v):
        called.append(v)
        if v==1.:
            entered.set();assert release.wait(3)
        return {"native":v}
    runner.start("fixed",[1.,2.,3.],resolve)
    assert entered.wait(3)
    runner.cancel();assert runner.snapshot()["status"]=="cancelled"
    runner.start("fixed",[4.],resolve)
    release.set()
    snapshot=finish(runner)
    assert called==[1.,4.] and snapshot["rows"]==[(4.,{"native":4.})]
    assert snapshot["completed"]==1 and snapshot["status"]=="complete"


def test_failed_row_retains_complete_rows_and_does_not_fill_zero():
    runner=RowSweep();runner.seed("fixed",[(1.,{"native":None})])
    def fail(v):
        raise ValueError("native solve failed")
    runner.start("fixed",[1.,2.,3.],fail)
    snapshot=finish(runner)
    assert snapshot["status"]=="failed" and snapshot["error"]=="native solve failed"
    assert snapshot["rows"]==[(1.,{"native":None})]


def test_retained_solves_survive_global_or_row_cache_eviction_and_stay_out_of_exports():
    runner=RowSweep(max_entries=1)
    def resolve(value):
        return {"value":value,"_native_raw":{"x":[value]},"_native_view":{"T":value}}
    runner.start("fixed",[1.,2.,3.],resolve)
    snapshot=finish(runner)
    assert all(not any(key.startswith("_native_") for key in payload) for _,payload in snapshot["rows"])
    # All rows still plotted by the current job retain their native solution,
    # including the first row that the bounded reusable row cache has evicted.
    raw,view=runner.solution("fixed",1.)
    assert raw=={"x":[1.]} and view=={"T":1.}
    raw["x"][0]=999
    assert runner.solution("fixed",1.)[0]=={"x":[1.]}
    assert runner.solution("different fixed inputs",1.) is None


def test_resume_keeps_completed_rows_and_retries_only_pending_conditions():
    runner=RowSweep();fail_once=[True];called=[]
    def resolve(value):
        called.append(value)
        if value==2. and fail_once[0]:
            fail_once[0]=False
            raise ValueError("temporary failure")
        return {"native":value}
    runner.start("fixed",[1.,2.,3.],resolve)
    assert finish(runner)["status"]=="failed"
    runner.resume()
    snapshot=finish(runner)
    assert called==[1.,2.,2.,3.] and snapshot["completed"]==3
    assert snapshot["requested_values"]==[1.,2.,3.]


def test_frozen_comparison_selection_does_not_move_sweep_range():
    data=sample_data();params=dict(power=1.,temperature=20.)
    scheme=SimpleNamespace(name="test",cache_version="v1")
    runner=RowSweep()
    entry=dict(identity=workspace_key(scheme,params,"power","readout1"),dataset=data,runner=runner)
    assert reusable_dataset(entry,scheme,{**params,"power":2.},"readout1") is data
    assert reusable_dataset(entry,scheme,{**params,"power":4.},"readout1") is None
    assert reusable_dataset(entry,scheme,{**params,"temperature":21.},"readout1") is None
    assert reusable_dataset(entry,scheme,params,"changed readout") is None


def test_selection_and_parameter_commit_use_authoritative_rows():
    data=sample_data();scheme=SimpleNamespace(param_schema=lambda:[ParamSpec("power","Power","Test",1.,0.,10.)])
    selected=workspace_selection(data,None,"run",dict(power=2.))
    assert selected["series_id"]=="2.0"
    assert selected["x"]==.2
    assert selected_parameter_change(scheme,data,{**selected,"condition":{"key":"temperature","value":99}},dict(power=1.))=={"power":2.}
    assert selected_parameter_change(scheme,data,{**selected,"series_id":"forged"},dict(power=1.))=={}
    carried=workspace_selection(data,{**selected,"display_mode":"map","panel_key":"bad"},"run")
    assert carried["display_mode"]=="map" and carried["panel_key"]=="t"
    reset=workspace_selection(data,selected,"other run",dict(power=1.))
    assert reset["series_id"]=="1.0" and reset["display_mode"]=="curves"


def test_real_lambda_row_is_the_native_solve_without_extra_comparison_conditions():
    scheme=schemes.get("lambda")
    params={**scheme.defaults(),"doppler":"off","temp_c":27.,"cell_mm":1.}
    calls=[]
    def resolve(condition):
        calls.append(condition["coupling_power_mw"])
        return scheme.compute(condition)
    payload=resolve_row(scheme,params,"coupling_power_mw",2.,resolve)
    assert calls==[2.]
    raw=scheme.compute({**params,"coupling_power_mw":2.})
    native=scheme.headless_observables(raw,{**params,"coupling_power_mw":2.})["lambda_readout"]
    np.testing.assert_array_equal(payload["series"]["x"],native["x_mhz"])
    np.testing.assert_array_equal(payload["series"]["observables"]["transmission"],native["transmission"])
    np.testing.assert_array_equal(payload["series"]["observables"]["dispersion"],native["dispersion"])


def test_scheme_switch_cancels_only_other_scheme_jobs():
    calls=[]
    session={"other":123,"spectral-workspace-lambda":{"runner":SimpleNamespace(cancel=lambda:calls.append("lambda"))},
             "spectral-workspace-magneto":{"runner":SimpleNamespace(cancel=lambda:calls.append("magneto"))}}
    cancel_inactive_sweeps(session,"magneto")
    assert calls==["lambda"]


@pytest.mark.parametrize("name,key,value",[("rydberg_eit","residual_zeeman_mhz",.2),("magneto","qwp_deg",18.)])
def test_other_schemes_resolve_one_native_condition_without_changing_arrays(name,key,value):
    scheme=schemes.get(name)
    params={**scheme.defaults(),"doppler":"off","scan_points":11,"velocity_points":2,
            "show_shot_noise":False}
    if name == "rydberg_eit": params["eit_compare"] = "Residual width"
    seen=[];raws=[]
    def resolve(condition):
        seen.append(condition[key]);raw=scheme.compute(condition);raws.append(raw);return raw
    payload=resolve_row(scheme,params,key,value,resolve)
    assert seen==[value]
    condition={**params,key:value}
    if name=="magneto":
        native=scheme.headless_observables(raws[0],condition)["magneto_readout"]
        np.testing.assert_array_equal(payload["series"]["x"],native["x"])
        np.testing.assert_allclose(np.array(payload["series"]["observables"]["transmission_pct"]),np.asarray(native["transmission"])*100)
    else:
        x,transmission,_=scheme._transmission(raws[0]["chi_bar_uncomp"],raws[0],condition)
        np.testing.assert_array_equal(payload["series"]["x"],x)
        np.testing.assert_array_equal(payload["series"]["observables"]["transmission"],transmission)
