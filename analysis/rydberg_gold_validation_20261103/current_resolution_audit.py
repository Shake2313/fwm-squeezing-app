"""Predeclared sparse guards for the corrected model; retain resolution evidence."""
from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

import numpy as np

from benchmark_gold import BenchmarkScheme, parameters, rf_protocol, source_hashes

ROOT = Path(__file__).resolve().parent
GUARDS = [(temperature, power) for temperature in (17., 26., 32., 34.) for power in (25., 35.)]
TOLERANCE = .01


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def density(temperature):
    scheme = BenchmarkScheme()
    params = parameters(scheme, "on", "Probe detuning", temp_c=temperature)
    return scheme._local_density_from_temperature_state(scheme._temperature_state(params))


def main():
    before = source_hashes()
    benchmark = json.loads((ROOT/"benchmark_results.json").read_text(encoding="utf-8"))
    if benchmark["source_sha256"] != before:
        raise RuntimeError("Benchmark physical source fingerprint differs")
    pilot = list(csv.DictReader((ROOT/"current_model_sparse_pilot.csv").open(encoding="utf-8-sig")))
    values = {}
    for row in pilot:
        values[row["detuning_policy"], float(row["temperature_c"]), float(row["coupling_power_mw"])] = float(row["sensitivity_nv_cm_sqrt_hz"])
    guard_plan = dict(origin="⟦AI⟧ Guard points fixed before new guard evaluations",
        variables=dict(temperature_c=[17.,26.,32.,34.], coupling_power_mw=[25.,35.]),
        points=GUARDS, criterion_relative=TOLERANCE,
        criterion_scope="AI numerical decision criterion, not experimental uncertainty or user-registered goal",
        training_scope="3x3 fixed-protocol pilot; only its5cross points parameterize the predictor;4corner guards remain separate",
        predictor="log-linear crossed slices; compare raw log S and log[S*N(T)] residuals with exact model density")
    (ROOT/"current_guard_plan.json").write_text(json.dumps(guard_plan, ensure_ascii=False, indent=2)+"\n", encoding="utf-8")
    rows = []
    for temperature, coupling in GUARDS:
        rows.extend(rf_protocol("on", "Probe detuning", temp_c=temperature, coupling_power_mw=coupling))
    with (ROOT/"current_guard_points.csv").open("w", encoding="utf-8-sig", newline="") as stream:
        writer=csv.DictWriter(stream, fieldnames=list(rows[0])); writer.writeheader(); writer.writerows(rows)
    reports = []
    tnodes, pnodes = np.array([15.,20.,35.]), np.array([20.,30.,40.])
    for policy in ("maximum_static_slope", "bare_two_photon_center"):
        for normalize_density in (False, True):
            thermal = np.array([values[policy,t,30.]* (density(t) if normalize_density else 1.) for t in tnodes])
            power = np.array([values[policy,20.,p]/values[policy,20.,30.] for p in pnodes])
            predictions = []
            for row in rows:
                if row["detuning_policy"] != policy:
                    continue
                temperature, coupling = row["temperature_c"], row["coupling_power_mw"]
                predicted = np.exp(np.interp(temperature,tnodes,np.log(thermal))
                                   +np.interp(coupling,pnodes,np.log(power)))
                if normalize_density:
                    predicted /= density(temperature)
                actual = row["sensitivity_nv_cm_sqrt_hz"]
                predictions.append(dict(temperature_c=temperature,coupling_power_mw=coupling,
                    direct_sensitivity=actual,predicted_sensitivity=float(predicted),relative_error=float(predicted/actual-1)))
            errors = np.array([row["relative_error"] for row in predictions])
            reports.append(dict(detuning_policy=policy,density_normalized=normalize_density,
                training_cross_coordinates=5, independent_new_guard_coordinates=len(predictions),
                maximum_relative_error=float(abs(errors).max()),rms_relative_error=float(np.sqrt(np.mean(errors**2))),
                eligible_at_declared_criterion=bool(abs(errors).max() <= TOLERANCE),predictions=predictions))
    after = source_hashes()
    if after != before:
        raise RuntimeError("Physical source changed during new guard evaluations")
    certificate = dict(origin="⟦AI⟧ Corrected-model sparse prediction check",
        model="correlated_phase_lindblad",protocol="Fixed Pp6uW, LO3.7MHz, IF40kHz; two separately declared optical-point policies",
        training_axes=dict(temperature_c=tnodes.tolist(),coupling_power_mw=pnodes.tolist()),
        guard_plan=guard_plan,reports=reports,
        scope="Current CP model, doppler on, delta_c0, L50mm,d.15mm,opacity.001; not a global reoptimized map or an experiment validation",
        physical_validation_status="Spectroscopy not reproduced quantitatively; protocol and absolute gain remain unresolved",
        source_sha256=before,artifact_sha256={name:digest(ROOT/name) for name in [
            "current_resolution_audit.py","current_model_sparse_pilot.csv","current_guard_plan.json","current_guard_points.csv"]})
    (ROOT/"current_resolution_audit.json").write_text(json.dumps(certificate,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    variables = [
        dict(variable="temperature_c",role="Important nonlinear covariate",recommendation="Use exact N(T) baseline and sparse residual-temperature slice with guards",evidence="current_resolution_audit.json",scope="Fixed Pp/LO/probe-scan protocol only; new dephasing or geometry requires recheck"),
        dict(variable="coupling_power_mw",role="Smooth control with weak T interaction in tested CP kernel",recommendation="Sparse1Dslice plus independent cross guards; no removal from physical design",evidence="current_resolution_audit.json",scope="20–40mW;3training nodes; tested8new cross guards"),
        dict(variable="probe_power_uw",role="Control interacting with LO and readout noise",recommendation="Retain explicit Pp–LO interaction; source low-power width and high-power measured sensitivity demand different spacing",evidence="fixed_lo_rf.csv and sources/fig2b_points.csv",scope="Do not infer fixed-LO monotonicity from jointly profiled optimum counts"),
        dict(variable="lo_rabi_mhz",role="Control fixed for Gold reproduction, later joint optimization axis",recommendation="Validate3.7MHz first; then local Pp–LO search after spectroscopy/response gates",evidence="benchmark_results.json",scope="No current calibrated LO optimum"),
        dict(variable="probe/coupling_detuning_mhz",role="Independent physical coordinates",recommendation="Keep scan-axis and fixed counterpart explicit; signed candidates; compare source protocol hypotheses",evidence="independent_detuning_pilot.csv",scope="Source absolute detunings/scan identity not reported"),
        dict(variable="opacity/Rabi/dephasing/beam_parameters",role="Calibration or nuisance parameters",recommendation="Profile against spectra and independent gain/density data; never optimize them freely for lowest S",evidence="diagnostics/model_diagnostic_summary.json",scope="Gaussian candidate improves shape but unvalidated; boundaryGamma0 not adopted"),
        dict(variable="radial_nodes_and_frequency_step",role="Numerical resolution",recommendation="Retain selected12nodes/.0333MHz coarsening evidence; do not spend resolution to fit a much larger model residual",evidence="diagnostics/gaussian_trial.json",scope="Selected Gaussian hypothesis only; joint numerical configuration certificate"),
    ]
    registry=dict(origin="⟦AI⟧ Current-model variable-resolution evidence registry",
        adoption_status="Analysis evidence and recommendations; not a new human-registered goal",
        physical_model_validation="NOT_ESTABLISHED",model="correlated_phase_lindblad",variables=variables,
        current_sparse_evidence=reports,
        historical_registry="reduction/variable_resolution_registry.json",
        historical_evidence_status="Retained historical results; not automatically transferable after the CP/axis change",
        source_sha256=before,
        invalidate_when=["Critical physical source changes","Pp or LO optimized rather than fixed","Optical axis/fixed detuning changes","Opacity/noise/geometry/Rabi calibration changes","Temperature/dephasing policy changes"],
        evidence_manifest_sha256=digest(ROOT/"current_resolution_audit.json"))
    (ROOT/"current_variable_resolution_registry.json").write_text(json.dumps(registry,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    print(json.dumps([{k:r[k] for k in ("detuning_policy","density_normalized","maximum_relative_error","eligible_at_declared_criterion")} for r in reports],indent=2))


if __name__ == "__main__":
    main()
