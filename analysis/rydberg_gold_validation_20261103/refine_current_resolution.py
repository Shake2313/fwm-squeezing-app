"""Refine one sparse axis after a failed guard; test on fresh guard points."""
from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

import numpy as np

from benchmark_gold import rf_protocol, source_hashes
from current_resolution_audit import density

ROOT = Path(__file__).resolve().parent
NEW_GUARDS = [(t,p) for t in (18.5,23.,29.,33.) for p in (22.5,37.5)]


def sha(name):
    return hashlib.sha256((ROOT/name).read_bytes()).hexdigest()


def write(name, value):
    (ROOT/name).write_text(json.dumps(value,ensure_ascii=False,indent=2,allow_nan=False)+"\n",encoding="utf-8")


def export(name, rows):
    with (ROOT/name).open("w",encoding="utf-8-sig",newline="") as stream:
        writer=csv.DictWriter(stream,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)


def main():
    before=source_hashes()
    old=json.loads((ROOT/"current_resolution_audit.json").read_text(encoding="utf-8"))
    if before!=old["source_sha256"]:
        raise RuntimeError("Physics changed after first sparse guard")
    plan=dict(origin="⟦AI⟧ Adaptive sparse-axis refinement chosen before fresh guard calculations",
        trigger="V1 maximum-static-slope density-normalized predictor exceeded1% (1.0633%); V1 results preserved",
        refinement="Add Pc25 and35mW on T20 slice; retain T15/20/35 with exact N(T)",
        new_training_coordinates=[[20.,25.],[20.,35.]],fresh_guard_coordinates=NEW_GUARDS,
        guard_scope="None of the fresh guard coordinates was used in V1 or predictor selection",
        tolerance_relative=.01,tolerance_status="AI engineering criterion, not experimental uncertainty or a human-registered goal")
    write("refined_guard_plan.json",plan)
    with (ROOT/"current_model_sparse_pilot.csv").open(encoding="utf-8-sig") as stream:
        old_rows=list(csv.DictReader(stream))
    training={(r["detuning_policy"],float(r["temperature_c"]),float(r["coupling_power_mw"])):
              float(r["sensitivity_nv_cm_sqrt_hz"]) for r in old_rows}
    axis_rows=[]
    for power in (25.,35.):
        found=rf_protocol("on","Probe detuning",temp_c=20.,coupling_power_mw=power)
        axis_rows.extend(found)
        for row in found:
            training[row["detuning_policy"],20.,power]=row["sensitivity_nv_cm_sqrt_hz"]
    export("refined_axis_points.csv",axis_rows)
    guard_rows=[]
    for temperature,power in NEW_GUARDS:
        guard_rows.extend(rf_protocol("on","Probe detuning",temp_c=temperature,coupling_power_mw=power))
    export("refined_guard_points.csv",guard_rows)
    tnodes=np.array([15.,20.,35.]);pnodes=np.array([20.,25.,30.,35.,40.])
    reports=[]
    for policy in ("maximum_static_slope","bare_two_photon_center"):
        thermal=np.array([training[policy,t,30.]*density(t) for t in tnodes])
        power=np.array([training[policy,20.,p]/training[policy,20.,30.] for p in pnodes])
        predictions=[]
        for row in guard_rows:
            if row["detuning_policy"]!=policy:
                continue
            t,p=row["temperature_c"],row["coupling_power_mw"]
            estimate=np.exp(np.interp(t,tnodes,np.log(thermal))+np.interp(p,pnodes,np.log(power)))/density(t)
            actual=row["sensitivity_nv_cm_sqrt_hz"]
            predictions.append(dict(temperature_c=t,coupling_power_mw=p,predicted=float(estimate),
                direct=actual,relative_error=float(estimate/actual-1)))
        error=np.array([row["relative_error"] for row in predictions])
        reports.append(dict(detuning_policy=policy,maximum_relative_error=float(abs(error).max()),
            rms_relative_error=float(np.sqrt(np.mean(error**2))),eligible_at_one_percent=bool(abs(error).max()<=.01),
            training_cross_coordinates=7,fresh_guard_coordinates=8,predictions=predictions))
    if source_hashes()!=before:
        raise RuntimeError("Physical source changed during adaptive check")
    files=["refine_current_resolution.py","current_resolution_audit.json","current_model_sparse_pilot.csv",
           "refined_guard_plan.json","refined_axis_points.csv","refined_guard_points.csv"]
    result=dict(origin="⟦AI⟧ Current CP model adaptive variable-resolution evidence",plan=plan,reports=reports,
        actual_unique_atomic_control_coordinates_across_v1_and_v2=27,
        model_scope="Fixed Pp6/LO3.7/IF40kHz, doppler on, probe scan, fixed delta_c0, opacity.001, L50mm,d.15mm",
        scope_limit="Tests declared guards only; no proof of error everywhere or of experimental model validity",
        experimental_validation_status="NOT_ESTABLISHED",source_sha256=before,
        artifact_sha256={name:sha(name) for name in files})
    write("refined_resolution_audit.json",result)
    registry=json.loads((ROOT/"current_variable_resolution_registry.json").read_text(encoding="utf-8"))
    registry["adaptive_refinement_evidence"]="refined_resolution_audit.json"
    registry["adaptive_refinement_evidence_sha256"]=sha("refined_resolution_audit.json")
    registry["previous_failed_guard_preserved"]="current_resolution_audit.json"
    registry["selected_current_policy"]=dict(temperature_nodes_c=tnodes.tolist(),coupling_power_nodes_mw=pnodes.tolist(),
        exact_density_baseline=True,reports=reports,scope=result["model_scope"],
        status="ELIGIBLE_ON_TESTED_GUARDS" if all(r["eligible_at_one_percent"] for r in reports) else "FURTHER_REFINEMENT_REQUIRED")
    for row in registry["variables"]:
        if row["variable"] in ("temperature_c","coupling_power_mw"):
            row["adaptive_evidence"]="refined_resolution_audit.json"
            row["proposed_training_nodes"]=(tnodes if row["variable"]=="temperature_c" else pnodes).tolist()
        if row["variable"] == "coupling_power_mw":
            row["scope"] = "20–40 mW; 5 training nodes after adaptive refinement; tested 8 fresh cross guards"
        if row["variable"] == "probe/coupling_detuning_mhz":
            row["scope"] = ("Source absolute detunings/scan identity not reported. "
                "Current linewidth comparison is restricted to zero fixed optical detuning; "
                "bare two-photon centre is not a fitted spectral maximum at general detuning.")
    write("current_variable_resolution_registry.json",registry)
    print(json.dumps([{k:r[k] for k in ("detuning_policy","maximum_relative_error","eligible_at_one_percent")} for r in reports],indent=2))


if __name__=="__main__":
    main()
