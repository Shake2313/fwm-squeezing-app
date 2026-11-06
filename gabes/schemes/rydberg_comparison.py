"""Native Rydberg maps, conditional finite-IF readouts and control analyses."""
import numpy as np

from .. import __version__, rydberg_opacity as opacity
from ..model_comparison import condition_views, sample_values

SCHEMA = "gabes-spectral-comparison-v1"
_KEYS = {"Probe power":"probe_power_uw", "Coupling power":"coupling_power_mw",
         "RF Rabi":"lo_rabi_mhz", "RF detuning":"mw_detuning_mhz",
         "IF frequency":"if_khz", "Residual width":"residual_zeeman_mhz"}


def comparison_key(params):
    choice=params.get("eit_compare","Probe power") if params.get("view","EIT")=="EIT" else params.get("at_compare","RF Rabi")
    if choice=="Temperature":
        if params.get("temperature_model","Linked")=="Linked":
            return "temp_c"
        return "cold_spot_temp_c" if params.get("temperature_compare")=="Cold-spot temperature" else "effective_temp_c"
    return _KEYS.get(choice,"probe_power_uw")


def _number(value):
    return float(value) if value is not None and np.isfinite(value) else None


def _reading(key,label,value,unit="",*,status=None,reason=None):
    result=dict(key=key,label=label,value=_number(value),unit=unit)
    if status is not None: result["status"]=status
    if reason: result["reason"]=reason
    return result


def _if_status(status):
    return {"PREDICTED":"predicted · uncalibrated", "NOT_APPLICABLE_NO_LO":"RF off",
            "UNAVAILABLE_ZERO_OPACITY":"zero effective opacity"}.get(status,"no finite optical response")


def _values(params,key,spec):
    current=float(params[key])
    if key=="probe_power_uw": values=[.5,1.,2.,4.,6.,8.,10.]
    elif key=="coupling_power_mw": values=[current/4,current/2,current,current*2,current*4]
    elif key=="lo_rabi_mhz": values=[0.,1.,3.7,7.4]
    elif key=="mw_detuning_mhz": values=np.linspace(current-5.,current+5.,7).tolist()
    elif key=="if_khz": values=[1.,10.,40.,100.,200.,500.]
    elif key=="residual_zeeman_mhz": values=[0.,current]
    else: values=[20.,30.,35.,40.,42.,45.,47.,50.,53.,55.]
    values=sorted({float(np.clip(v,spec.vmin,spec.vmax)) for v in values}|{current})
    if str(params.get("doppler","off"))=="on" and len(values)>3:
        values=sorted({values[0],current,values[-1]})
    return values


def _rf_readouts(scheme,raw,params,readout):
    rf=float(raw["lo_rabi_mhz"])
    split=readout["at_split_mhz"]
    coupling=scheme._rf_coupling(params)
    input_field=coupling.field_from_cyclic_rabi_hz(rf*1e6)
    detuning=abs(float(params["mw_detuning_mhz"]))
    inferred=(coupling.field_from_at_splitting_hz(split*1e6,detuning_hz=detuning*1e6)
              if rf>0 and split is not None and split>=detuning else None)
    status="not applicable" if rf==0 else "unresolved" if split is None else "resolved"
    unit="mV/m "+str(params["rf_field_convention"])
    return [
        _reading("split","Optical transmission-peak separation",split,"MHz",status=status,
                 reason="RF is off." if rf==0 else "The native peak detector did not resolve two transmission maxima." if split is None else None),
        _reading("input_field","Field from input ΩRF",input_field*1e3,unit,status="off" if rf==0 else None),
        _reading("inferred_field","Field estimated from optical peak separation",inferred*1e3 if inferred is not None else None,unit,
                 status=status if inferred is not None or split is None or rf==0 else "unavailable",
                 reason="Ideal two-level conversion; optical line pulling is not corrected."),
        _reading("center_transmission","Bare-resonance transmission",readout["resonance_transmission"])]


def _features(scheme,raw,params,chi):
    x,t,_=scheme._transmission(chi,raw,params)
    width,contrast=scheme._eit_features(x,t,center_mhz=scheme._resonance_scan_mhz(params))
    return dict(width=_number(width) if np.isfinite(width) and width>0 else None,
                contrast=_number(contrast),amplitude=_number(raw["probe_power_uw"]*contrast))


def _analysis(scheme,raw,params,superhet):
    comp=_features(scheme,raw,params,raw["chi_bar"])
    residual=_features(scheme,raw,params,raw["chi_bar_uncomp"])
    return dict(comp=comp,residual=residual,
                total=(_number(superhet["sensitivity"].total_field_asd_nv_cm_per_sqrt_hz) if superhet else None),
                psn=(_number(superhet["sensitivity"].psn_field_asd_nv_cm_per_sqrt_hz) if superhet else None),
                response=(_number(superhet["current_responsivity_a_per_v_m"]) if superhet else None))


def _native_row(scheme,condition,raw,view,key,index):
    native=view["rydberg_spectrum"]
    x=np.asarray(native["x"]); t=np.asarray(native["transmission"]); sh=native["finite_if"]
    value=float(condition[key]); spec=next(sp for sp in scheme.param_schema() if sp.name==key)
    at=condition["view"]!="EIT"
    label=("RF off" if key=="lo_rabi_mhz" and value==0 else
           f"ΩRF/2π = {value:g} MHz" if key=="lo_rabi_mhz" else f"{value:g} {spec.unit}")
    observations=dict(transmission=sample_values(t))
    grids={}
    if sh is not None:
        response_x=sample_values(sh["response_detuning_mhz"])
        current=np.asarray(sh["current_phasors_a_per_v_m"])
        phase=np.where(np.abs(current)>0,np.angle(current,deg=True),np.nan)
        observations.update(
            if_response=sample_values(np.abs(sh["transmission_phasor_per_angular_rabi"])*2*np.pi*1e3),
            if_phase=sample_values(phase),
            total_sensitivity=sample_values(sh["total_field_asd_v_m_per_sqrt_hz"]*1e7),
            psn_sensitivity=sample_values(sh["psn_field_asd_v_m_per_sqrt_hz"]*1e7))
        grids={name:response_x for name in ("if_response","if_phase","total_sensitivity","psn_sensitivity")}
    elif at:
        for name in ("if_response","if_phase","total_sensitivity","psn_sensitivity"):
            observations[name]=[]; grids[name]=[]
    rows=_rf_readouts(scheme,raw,condition,view["rydberg_readout"]) if at else []
    analysis=_analysis(scheme,raw,condition,sh)
    if not at:
        width=analysis["comp"]["width"]
        rows=[_reading("fwhm","EIT FWHM",width,"MHz",status="resolved" if width is not None else "unresolved",
                       reason="Native interpolated crossings of the center-to-local-floor feature within ±5 MHz."),
              _reading("contrast","Center-to-local-floor contrast",analysis["comp"]["contrast"]),
              _reading("amplitude","Transmitted feature amplitude",analysis["comp"]["amplitude"],"µW",
                       reason="Incident probe power × transmission contrast; excludes path/detector efficiency."),
              _reading("center_transmission","Resonance transmission",view["rydberg_readout"]["resonance_transmission"])]
    rows += [_reading("density","Vapor density",raw["N"],"m⁻³"),
             _reading("atoms","Effective atom number",scheme._effective_atom_number(raw,condition).atoms),
             _reading("effective_temperature","Effective vapor temperature",raw["effective_temp_c"],"°C"),
             _reading("cold_spot","Cold-spot temperature",raw["cold_spot_temp_c"],"°C"),
             _reading("heater","Heater setpoint",raw["heater_setpoint_c"],"°C")]
    optimum=None
    if at:
        finite=view["rydberg_readout"]["finite_if"]
        status=_if_status(finite["status"])
        rows += [_reading("if_frequency","IF frequency",condition["if_khz"],"kHz"),
                 _reading("best_total","Best sampled total sensitivity",analysis["total"],"nV/cm/√Hz",status=status),
                 _reading("best_psn","PSN sensitivity at that point",analysis["psn"],"nV/cm/√Hz"),
                 _reading("best_response","Current responsivity at that point",analysis["response"],"A/(V/m)"),
                 _reading("best_x","Best sampled optical coordinate",finite["optimum_detuning_mhz"],"MHz",
                          reason=finite.get("optimization_scope") or status)]
        if sh is not None:
            optimum=dict(x=float(sh["optimum_detuning_mhz"]),label="Inspect sampled optimum",
                         scope=sh["optimization_scope"])
    markers=[dict(x=scheme._resonance_scan_mhz(condition),label="Bare resonance")]
    markers += [dict(x=p,label=f"Transmission peak {i+1}") for i,p in enumerate(native["at_peak_positions_mhz"])]
    if optimum: markers.append(dict(x=optimum["x"],label="Best sampled sensitivity",profiles=["sensitivity"]))
    _,uncomp,_=scheme._transmission(raw["chi_bar_uncomp"],raw,condition)
    row=dict(id=f"condition-{index}",label=label,parameter=dict(key=key,label=spec.label,value=value,unit=spec.unit),
             x=sample_values(x),observables=observations,observable_x=grids,readouts=rows,
             unavailable_observables=({name:_if_status(view["rydberg_readout"]["finite_if"]["status"])
                                       for name in grids} if at and sh is None else {}),
             markers=sorted(markers,key=lambda m:m["x"]),analysis=analysis,optimum=optimum,
             finite_if_scope=(sh["optimization_scope"] if sh else view["rydberg_readout"]["finite_if"]["status"]),
             temperature=dict(model=condition["temperature_model"],heater_setpoint_c=raw["heater_setpoint_c"],
                              effective_temp_c=raw["effective_temp_c"],cold_spot_temp_c=raw["cold_spot_temp_c"]),
             overlays=({"residual":dict(panel="transmission",label=f"Residual width {raw['residual_zeeman_mhz']:g} MHz",
                                        x=sample_values(x),y=sample_values(uncomp))} if raw["residual_zeeman_mhz"]>0 else {}))
    row["opacity_scenarios"]=[dict(label="nominal",value=opacity.opacity_spec(condition)["value"],analysis=analysis)]
    atomic=sh["atomic_response"] if sh is not None else None
    for scenario in opacity.opacity_samples(condition):
        if scenario["label"]=="nominal": continue
        scenario_params={**condition,"effective_opacity_participation":scenario["value"]}
        scenario_ro=scheme._readout(raw,scenario_params,atomic_response=atomic)
        row["opacity_scenarios"].append(dict(**scenario,analysis=_analysis(scheme,raw,scenario_params,scenario_ro["superhet"])))
    return row


def _summary(rows,spec,at):
    xs=[r["parameter"]["value"] for r in rows]
    ids=[r["id"] for r in rows]
    series=[];panels=[]
    def add(sid,label,values,panel,unit,*,optional=None,goal=None,scale=None,note=None):
        series.append(dict(id=sid,label=label,x=xs,condition_ids=ids,observables=dict(value=values)))
        target=next((p for p in panels if p["key"]==panel),None)
        if target is None:
            target=dict(key=panel,title=label,y_label=label,unit=unit,curves=[],note=note)
            if scale: target["scale"]=scale
            panels.append(target)
        definition=dict(series_id=sid,observable="value")
        if optional: definition["optional_overlays"]=optional
        target["curves"].append(definition)
        if goal: target.update(goal=goal,optimum_series_id=sid)
    if at:
        field_unit=next(r["unit"] for r in rows[0]["readouts"] if r["key"]=="input_field")
        for reading,sid,label in [("input_field","input-field","From input ΩRF"),("inferred_field","split-field","Estimated from optical peak separation")]:
            add(sid,label,[next(r["value"] for r in row["readouts"] if r["key"]==reading) for row in rows],"field",field_unit)
        panels[-1].update(title="Input field and peak-based estimate",y_label="Electric field")
        add("split","Optical transmission-peak separation",[next(r["value"] for r in row["readouts"] if r["key"]=="split") for row in rows],"split","MHz")
        add("rf-reference","Input ΩRF/2π reference",[row["rf_rabi_mhz"] for row in rows],"split","MHz")
    for scenario_name in ("nominal","lower","upper"):
        scenarios=[next((s for s in row["opacity_scenarios"] if s["label"]==scenario_name),None) for row in rows]
        opacity_option=[] if scenario_name=="nominal" else ["opacity"]
        suffix="" if scenario_name=="nominal" else " · assumed opacity "+scenario_name
        if at:
            for metric,label,unit,goal in [("response","Current responsivity at sampled optimum","A/(V/m)",None),
                                          ("total","Best sampled total sensitivity","nV/cm/√Hz","min"),
                                          ("psn","PSN at the total-sensitivity optimum","nV/cm/√Hz",None)]:
                add(metric+"-"+scenario_name,label+suffix,[s["analysis"][metric] if s else None for s in scenarios],
                    "response" if metric=="response" else "sensitivity",unit,optional=opacity_option,
                    goal=goal if scenario_name=="nominal" else None,scale="log",
                    note="Computed conditions only. Lower field ASD is better. Opacity scenarios are assumptions, not confidence limits.")
        else:
            for variant in ("comp","residual"):
                optional=opacity_option+([] if variant=="comp" else ["residual"])
                for metric,label,unit,goal in [("contrast","EIT contrast","1",None),("width","EIT FWHM","MHz","min"),
                                              ("amplitude","Transmitted feature amplitude","µW","max")]:
                    add(metric+"-"+variant+"-"+scenario_name,label+(" · residual width" if variant=="residual" else "")+suffix,
                        [s["analysis"][variant][metric] if s else None for s in scenarios],metric,unit,optional=optional,
                        goal=goal if scenario_name=="nominal" and variant=="comp" else None,
                        note="Native ±5 MHz local-floor feature. Opacity scenarios are assumptions; optima are sampled model results.")
    return dict(axis=dict(label=spec.label,unit=spec.unit,reference="Actual varied control; other physical conditions fixed"),
                series=series,panels=panels)


def _residual_comparison(scheme,params,raw):
    curves=[]
    residual=float(raw["residual_zeeman_mhz"])
    for i,(value,chi) in enumerate([(0.,raw["chi_bar"])]+([(residual,raw["chi_bar_uncomp"])] if residual>0 else [])):
        x,t,_=scheme._transmission(chi,raw,params)
        width,contrast=scheme._eit_features(x,t,center_mhz=scheme._resonance_scan_mhz(params))
        curves.append(dict(id=f"residual-{i}",label="Residual broadening removed" if value==0 else f"Residual Zeeman FWHM {value:g} MHz",
                           parameter=dict(key="residual_zeeman_mhz",label="Residual Zeeman FWHM",value=value,unit="MHz"),
                           x=sample_values(x),observables=dict(transmission=sample_values(t)),
                           readouts=[_reading("fwhm","EIT FWHM",width if np.isfinite(width) and width>0 else None,"MHz",
                                              status="resolved" if np.isfinite(width) and width>0 else "unresolved"),
                                     _reading("center_transmission","Resonance transmission",t[np.argmin(np.abs(x-scheme._resonance_scan_mhz(params)))]),
                                     _reading("contrast","Center-to-local-floor contrast",contrast)]))
    return curves,dict(axis=dict(label="Residual Zeeman FWHM",unit="MHz"),
        series=[dict(id="width",label="EIT FWHM",x=[c["parameter"]["value"] for c in curves],
                     observables=dict(fwhm=[c["readouts"][0]["value"] for c in curves]))],
        panels=[dict(key="width",title="Residual broadening and EIT width",y_label="EIT FWHM",unit="MHz",
                     curves=[dict(series_id="width",observable="fwhm")])])


def spectral_comparison(scheme,params,*,resolve_raw,current_raw=None,current_view=None):
    params={**scheme.defaults(),**params}
    at=params["view"]!="EIT"; key=comparison_key(params)
    spec=next(sp for sp in scheme.param_schema() if sp.name==key)
    if current_raw is None: current_raw=resolve_raw(params)
    if not at and key=="residual_zeeman_mhz":
        rows,summary=_residual_comparison(scheme,params,current_raw)
    else:
        rows=[]
        for i,(p,raw,view) in enumerate(condition_views(scheme,params,key,_values(params,key,spec),resolve_raw=resolve_raw,
                                                       current_raw=current_raw,current_view=current_view)):
            if "rydberg_spectrum" not in view: view=scheme.headless_observables(raw,p)
            row=_native_row(scheme,p,raw,view,key,i)
            row["rf_rabi_mhz"]=float(raw["lo_rabi_mhz"])
            rows.append(row)
        summary=_summary(rows,spec,at)
    panels=[dict(key="transmission",title="Transmission spectra",y_label="Transmission T",unit="1",
                 curves=[dict(series_id=row["id"],observable="transmission") for row in rows])]
    convention=str(params["rf_field_convention"])
    if at:
        for observable,title,label,unit,extra in [
            ("if_response","Finite-IF response magnitude","|δT| per 1 kHz SIG","1",dict(cut_group="response")),
            ("if_phase","Finite-IF response phase (wrapped)","Current-response phase","°",dict(cyclic_period=360)),
            ("total_sensitivity","Total noise-equivalent field","Total field ASD · lower is better","nV/cm/√Hz",dict(scale="log",cut_group="sensitivity")),
            ("psn_sensitivity","PSN noise-equivalent field","PSN field ASD · lower is better","nV/cm/√Hz",dict(scale="log",cut_group="sensitivity"))]:
            panels.append(dict(key=observable,title=title,y_label=label,unit=unit,**extra,
                               sampled_only=str(params["doppler"])=="on",
                               curves=[dict(series_id=row["id"],observable=observable) for row in rows]))
    selected=next((r for r in rows if r["parameter"]["value"]==float(params[key])),rows[0])
    center=scheme._resonance_scan_mhz(params)
    half=max(2.5,2.5*(next((r["value"] for r in selected["readouts"] if r["key"]=="fwhm"),None) or 0.))
    if at: half=max(5.,float(params["lo_rabi_mhz"])*1.5+abs(float(params["mw_detuning_mhz"])))
    data=dict(schema=SCHEMA,scheme=scheme.name,mode=params["view"],title=params["view"]+" · "+spec.label+" comparison",
              axis=dict(label=scheme._scan_axis(params),unit="MHz",reference=f"Bare two-photon resonance at {center:g} MHz"),
              series=rows,panels=panels,summary=summary,
              display=dict(default_mode="map",window=[center-half,center+half]),
              fixed={k:v for k,v in params.items() if k!=key},
              varied=dict(key=key,label=spec.label,unit=spec.unit,values=[row["parameter"]["value"] for row in rows]),
              source=dict(model=type(scheme).__name__,package_version=__version__,cache_version=scheme.cache_version,
                          method="existing compute + native readout",interpolation="none",field_convention=convention,
                          prediction_status="PREDICTED",calibration_status="UNCALIBRATED",
                          effective_opacity=opacity.opacity_spec(params),temperature_model=params["temperature_model"]),
              notes=["Scalar ⁸⁵Rb four-level cascade. Conditional model predictions; opacity/dipole/angular factors remain explicit.",
                     "Every condition and observable retains its own solved optical coordinates. Missing or unsolved IF points are blank.",
                     "Heater setpoint, effective vapor temperature and cold spot keep their declared meanings."])
    if key!="residual_zeeman_mhz":
        data["optional_overlays"]=[dict(key="opacity",label="Show assumed opacity scenarios",default=False,**(dict(profiles=["sensitivity"]) if at else {}))]
        if not at and float(params["residual_zeeman_mhz"])>0: data["optional_overlays"].insert(0,dict(key="residual",label="Show residual-width comparison",default=False))
    if at:
        data["profiles"]=[
            dict(key="spectra",label="Spectra",panels=["transmission"],summaries=["field","split"],
                 readouts=["split","input_field","inferred_field","center_transmission","effective_temperature","cold_spot","heater","density"],
                 default_panel="transmission"),
            dict(key="sensitivity",label="Sensitivity",panels=[p["key"] for p in panels],summaries=["response","sensitivity"],
                 readouts=["if_frequency","best_total","best_psn","best_response","best_x","effective_temperature","cold_spot","heater","density","atoms"],
                 default_panel="total_sensitivity",inspect_optimum=True,
                 note=("Finite-IF linear response · "+("IF varies by condition" if key=="if_khz" else f"{params['if_khz']:g} kHz IF")+". Lower field ASD is better. Conditional, uncalibrated model; field ASD is invariant under the RF Peak/RMS reporting choice."),
                 readout_description="Best values refer to solved IF samples; the table reports the cursor's native sample.")]
        data["display"]["profile"]="sensitivity" if key=="if_khz" else "spectra"
        data["notes"] += ["Finite-IF Liouvillian response is separate from the static transmission gradient.",
                          "Doppler on: only the declared static-slope sample points are evaluated; no full-detuning optimum is claimed.",
                          "Peak = √2 × RMS for sinusoidal fields. Optical peak spacing can differ from input ΩRF through detuning and line pulling."]
    return data
