"""Auditable, scalar-only cross-slice reduction of saved Rydberg maps.

No atomic solves, calibration or physical-model changes. A policy is selected
using the earlier pilot alone; finer saved points then test that frozen rule.
Control labels are never interpolated or averaged.
"""
from __future__ import annotations

import csv
import hashlib
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Mapping

import numpy as np

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parents[2]
PILOT = REPO / "analysis/rydberg_gold_map_20261103"
DENSE = REPO / "analysis/rydberg_gold_dense_20261103"
FIELD = "sensitivity_nv_cm_sqrt_hz"
TOLERANCE = .01  # AI analysis/decision criterion, not a user-registered requirement.
REFERENCE = (20., 30.)
MODEL_TAG = "legacy_compact_RydbergEITScheme_profiled_PSN"
SWEEP_PARAMETER_KEYS = {
    'temp_c','heater_setpoint_c','effective_temp_c','cold_spot_temp_c',
    'coupling_power_mw','probe_power_uw','lo_rabi_mhz','doppler',
    'effective_opacity_participation',
}
CRITICAL_ENGINE = [
    "gabes/schemes/rydberg.py", "gabes/rydberg_electrometry.py",
    "gabes/rydberg_experiment.py", "gabes/rydberg_opacity.py", "gabes/core.py",
    "gabes/affine_doppler.py", "gabes/doppler.py", "gabes/pole_doppler.py",
    "gabes/kernels.py", "gabes/species.py", "gabes/beam.py", "gabes/constants.py",
    "gabes/atoms.py", "gabes/observables.py",
]


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write_json(path, value):
    Path(path).write_text(json.dumps(value, ensure_ascii=False, indent=2,
                                    allow_nan=False) + "\n", encoding="utf-8")


def write_csv(path, rows):
    if not rows:
        return
    with Path(path).open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def surface(rows, field=FIELD):
    """Validate a complete positive scalar surface; reject silent aggregation."""
    values = {}
    for row in rows:
        key = (float(row["temperature_c"]), float(row["coupling_power_mw"]))
        if key in values:
            raise ValueError(f"Duplicate map coordinate: {key}")
        value = float(row[field])
        if not np.isfinite(value) or value <= 0:
            raise ValueError(f"Nonpositive/nonfinite scalar at {key}")
        values[key] = value
    if not values:
        raise ValueError("Empty map")
    temperatures = sorted({key[0] for key in values})
    powers = sorted({key[1] for key in values})
    if len(values) != len(temperatures) * len(powers):
        raise ValueError("Incomplete rectangular map")
    return values, temperatures, powers


@dataclass(frozen=True)
class CrossSlice:
    """A scalar log-additive model fitted ONLY to two crossing 1D slices."""
    temperatures: tuple
    powers: tuple
    log_temperature_slice: tuple
    log_power_slice: tuple
    log_reference: float
    reference: tuple

    @classmethod
    def fit(cls, values: Mapping, reference=REFERENCE):
        temperatures = tuple(sorted({key[0] for key in values}))
        powers = tuple(sorted({key[1] for key in values}))
        if len(temperatures) < 2 or len(powers) < 2:
            raise ValueError("Two coordinates per axis required")
        t0, p0 = reference
        required = {(t, p0) for t in temperatures} | {(t0, p) for p in powers}
        if not required.issubset(values):
            raise ValueError("Reference-cross slices are incomplete")
        for key in required:
            if not np.isfinite(values[key]) or values[key] <= 0:
                raise ValueError("Cross requires positive finite scalars")
        return cls(temperatures, powers,
            tuple(float(np.log(values[t, p0])) for t in temperatures),
            tuple(float(np.log(values[t0, p])) for p in powers),
            float(np.log(values[t0, p0])), (float(t0), float(p0)))

    @property
    def fitted_coordinates(self):
        t0, p0 = self.reference
        return {(t, p0) for t in self.temperatures} | {(t0, p) for p in self.powers}

    def predict(self, temperature, power):
        if not self.temperatures[0] <= temperature <= self.temperatures[-1]:
            raise ValueError("Temperature extrapolation is not allowed")
        if not self.powers[0] <= power <= self.powers[-1]:
            raise ValueError("Coupling-power extrapolation is not allowed")
        log_value = (np.interp(temperature, self.temperatures, self.log_temperature_slice)
                     + np.interp(power, self.powers, self.log_power_slice)
                     - self.log_reference)
        return float(np.exp(log_value))


def evaluate(predictor, values, excluded=()):
    excluded = set(excluded)
    rows = []
    for (temperature, power), truth in sorted(values.items()):
        if (temperature, power) in excluded:
            continue
        prediction = predictor.predict(temperature, power)
        rows.append(dict(temperature_c=temperature, coupling_power_mw=power,
            observed=truth, predicted=prediction, relative_error=prediction / truth - 1))
    if not rows:
        raise ValueError("No evaluation coordinates remain")
    errors = np.array([row["relative_error"] for row in rows])
    summary = dict(count=len(rows), max_absolute_relative_error=float(abs(errors).max()),
        rms_relative_error=float(np.sqrt(np.mean(errors ** 2))),
        worst_coordinate=[rows[int(np.argmax(abs(errors)))][key]
                          for key in ("temperature_c", "coupling_power_mw")])
    return summary, rows


def axis_diagnostics(values):
    """Retrospective smoothness/step diagnostics of the PROFILED scalar surface."""
    temperatures=sorted({key[0] for key in values})
    powers=sorted({key[1] for key in values})
    log_values=np.array([[np.log(values[t,p]) for p in powers] for t in temperatures])
    gradients=np.gradient(log_values,temperatures,powers,edge_order=2)
    result=[]
    for axis,(name,coordinates,unit) in enumerate([
        ('temperature_c',temperatures,'1/degC'),('coupling_power_mw',powers,'1/mW')]):
        gradient=gradients[axis]
        curvature=np.gradient(gradient,coordinates,axis=axis,edge_order=2)
        difference=np.diff(log_values,axis=axis)
        step=max(np.diff(coordinates))
        result.append(dict(variable=name,log_gradient_min=float(gradient.min()),log_gradient_max=float(gradient.max()),
            log_gradient_unit=unit,max_absolute_log_curvature=float(abs(curvature).max()),
            max_step=float(step),max_step_squared_curvature=float(abs(curvature).max()*step**2),
            adjacent_relative_change_min=float(np.expm1(difference).min()),
            adjacent_relative_change_max=float(np.expm1(difference).max()),
            increasing_segment_count=int(np.sum(difference>1e-10)),decreasing_segment_count=int(np.sum(difference<-1e-10)),
            segment_count=int(difference.size),
            scope='Finite differences of saved profiled minima; not fixed-control causal derivatives or a certified interpolation bound'))
    return result


def select_from_pilot(values, tolerance=TOLERANCE):
    """No dense/future truth is accepted by this selection API."""
    predictor = CrossSlice.fit(values)
    error, _ = evaluate(predictor, values, predictor.fitted_coordinates)
    accepted = error["max_absolute_relative_error"] <= tolerance
    return predictor, dict(tolerance_relative=tolerance,
        tolerance_origin="AI-selected illustrative decision tolerance; not a registered user goal",
        training_coordinates=len(predictor.fitted_coordinates), pilot_check=error,
        accepted=accepted,
        selected_action="cross_slices_plus_guards" if accepted else "branch_aware_refinement",
        selection_scope="Earlier pilot scalars only; no fine-map values used to choose this rule")


def guard_plan(temperatures, powers, reference=REFERENCE):
    """Deterministic sparse design, selected without inspecting response values."""
    t0, p0 = reference
    cross = {(float(t), float(p0)) for t in temperatures} | {(float(t0), float(p)) for p in powers}
    def nearest(axis, fraction):
        target = axis[0] + fraction * (axis[-1] - axis[0])
        return min(axis, key=lambda value: (abs(value-target), value))
    guards = {(temperatures[i], powers[j]) for i in (0, -1) for j in (0, -1)}
    guards |= {(nearest(temperatures, fraction), nearest(powers, side))
               for fraction in (.5, .75) for side in (.25, .75)}
    guards -= cross
    return dict(reference=list(reference), fit_coordinates=[list(p) for p in sorted(cross)],
        guard_coordinates=[list(p) for p in sorted(guards)],
        fit_count=len(cross), guard_count=len(guards), initial_count=len(cross | guards),
        full_cartesian_count=len(temperatures)*len(powers),
        guard_origin="Axis endpoints and fixed fractional interior positions; no sensitivity/argmin inspection",
        failure_action="Retain expensive direct solves and refine cells near failing guards; do not average control branches")


def check_fingerprint(expected, root):
    mismatches = []
    for name, checksum in expected.items():
        path = Path(root) / name
        if not path.is_file() or digest(path) != checksum:
            mismatches.append(name)
    return mismatches


def fixed_settings_fingerprint(parameters):
    fixed={key:value for key,value in parameters.items() if key not in SWEEP_PARAMETER_KEYS}
    return hashlib.sha256(json.dumps(fixed,sort_keys=True,separators=(',',':'),allow_nan=False).encode()).hexdigest()


def reuse_rejections(registry, current_parameters, current_sources, model_tag):
    """Parameter/source guard for future callers; absence of rejection is scope eligibility only."""
    rejected=[]
    if model_tag != registry['model_tag']:
        rejected.append('model_tag_changed')
    if fixed_settings_fingerprint(current_parameters) != registry['fixed_model_settings_sha256']:
        rejected.append('fixed_physics_or_noise_settings_changed')
    if current_sources != registry['expected_critical_source_sha256']:
        rejected.append('numerical_sources_changed')
    for parameter,domain_key in [('temp_c','temperature_c'),('probe_power_uw','probe_power_uw'),
                                 ('coupling_power_mw','coupling_power_mw'),('lo_rabi_mhz','lo_rabi_mhz')]:
        low,high=registry['domain'][domain_key]
        if parameter in current_parameters and not low <= current_parameters[parameter] <= high:
            rejected.append(parameter+'_outside_domain')
    if current_parameters.get('doppler') not in ('off','on'):
        rejected.append('untested_doppler_policy')
    if current_parameters.get('effective_opacity_participation') not in registry['domain']['opacity_scenarios']:
        rejected.append('untested_opacity_scenario')
    return rejected


def branch_transitions(rows, branch_rows, tolerance=TOLERANCE):
    """Report discrete label changes without treating them as scalar curvature."""
    indexed = {(row["temperature_c"], row["coupling_power_mw"]): row for row in rows}
    branch = {(row["temperature_c"], row["coupling_power_mw"]): row for row in branch_rows}
    temperatures = sorted({key[0] for key in indexed})
    powers = sorted({key[1] for key in indexed})
    transitions = []
    for axis, first, second in (("temperature_c", temperatures, powers),
                                ("coupling_power_mw", powers, temperatures)):
        for held in second:
            keys = [(value, held) if axis == "temperature_c" else (held, value) for value in first]
            for left, right in zip(keys, keys[1:]):
                a, b = indexed[left], indexed[right]
                if a["winning_branch"] == b["winning_branch"]:
                    continue
                margins = []
                for key in (left, right):
                    pair = [branch[key][prefix + "_sensitivity_nv_cm_sqrt_hz"]
                            for prefix in ("central", "sideband")]
                    margins.append(abs(pair[0]-pair[1])/min(pair))
                transitions.append(dict(axis=axis, left=list(left), right=list(right),
                    left_branch=a["winning_branch"], right_branch=b["winning_branch"],
                    minimum_endpoint_competition_margin=min(margins),
                    at_least_one_endpoint_near_tie=min(margins) <= tolerance,
                    interpretation="Control identity changes; use within-coordinate branch margins, not interpolated argmin coordinates"))
    return transitions


def main():
    pilot = json.loads((PILOT / "results.json").read_text(encoding="utf-8"))
    dense = json.loads((DENSE / "results.json").read_text(encoding="utf-8"))
    interpretation = json.loads((DENSE / "interpretation.json").read_text(encoding="utf-8"))
    densities = {row["temperature_c"]: row["model_density_m3"]
                 for row in interpretation["model_temperature_diagnostics"]}
    audits, prediction_rows, branches, mode_decisions, axis_rows = [], [], [], [], []
    for mode in dense["config"]["modes"]:
        decisions = []
        for opacity in dense["config"]["opacities"]:
            old = [row for row in pilot["optima"] if row["doppler"] == mode and row["opacity"] == opacity]
            new = [row for row in dense["optima"] if row["doppler"] == mode and row["opacity"] == opacity]
            old_values, _, _ = surface(old)
            new_values, temperatures, powers = surface(new)
            axis_rows.extend(dict(doppler=mode,opacity=opacity,**row) for row in axis_diagnostics(new_values))
            # Freeze method selection using the old pilot before accessing the
            # dense truth in either evaluation below.
            old_predictor, decision = select_from_pilot(old_values)
            decisions.append(decision["accepted"])
            transferred, transferred_rows = evaluate(old_predictor, new_values, set(old_values))
            dense_predictor = CrossSlice.fit(new_values)
            cross, cross_rows = evaluate(dense_predictor, new_values, dense_predictor.fitted_coordinates)
            log_matrix = np.array([[np.log(new_values[t, p]) for p in powers] for t in temperatures])
            residual = log_matrix - log_matrix.mean(axis=0)[None, :] - log_matrix.mean(axis=1)[:, None] + log_matrix.mean()
            normalized = []
            for power in powers:
                sn = [new_values[t, power] * densities[t] for t in temperatures]
                normalized.append(dict(coupling_power_mw=power,
                    SN_max_to_min_ratio_15_to_35c=max(sn)/min(sn),
                    SN35_to_SN20_ratio=new_values[35., power]*densities[35.]/(new_values[20., power]*densities[20.])))
            audit = dict(doppler=mode, opacity=opacity, pilot_only_decision=decision,
                old_pilot_8_fit_10_check_scope="Same old sampled map; earlier evidence, not experimental validation",
                pilot_to_new_dense_coordinates=transferred,
                pilot_to_dense_scope="No fine-map fitting: log-linear interpolation of pilot cross; excludes all 18 pilot coordinates; error also contains changed detuning/control discretization",
                retrospective_dense_cross_29_fit_160_check=cross,
                dense_cross_scope="Retrospective structural diagnostic; method was inspected after dense campaign, not blind prospective validation",
                in_sample_log_additive_max_relative_residual=float(np.max(abs(np.exp(residual)-1))),
                known_inverse_density_baseline=normalized)
            audits.append(audit)
            for kind, predicted in (("pilot_to_new_dense", transferred_rows), ("retrospective_dense_cross", cross_rows)):
                prediction_rows.extend(dict(doppler=mode, opacity=opacity, audit_kind=kind, **row) for row in predicted)
            for prefix in ("central", "sideband"):
                subset = [row for row in dense["branchgaps"] if row["doppler"] == mode and row["opacity"] == opacity]
                values, _, _ = surface(subset, prefix + "_sensitivity_nv_cm_sqrt_hz")
                predictor = CrossSlice.fit(values)
                error, _ = evaluate(predictor, values, predictor.fitted_coordinates)
                branches.append(dict(doppler=mode, opacity=opacity, branch=prefix,
                    retrospective_cross_error=error,
                    scope="Independent branch-profiled scalar surfaces; no interpolation of Pp/LO/detuning labels"))
        mode_audits = [a for a in audits if a['doppler'] == mode]
        transfer_max = max(a['pilot_to_new_dense_coordinates']['max_absolute_relative_error'] for a in mode_audits)
        mode_decisions.append(dict(doppler=mode,
            all_opacity_pilot_checks_passed=all(decisions),
            selected_action="cross_slices_plus_guards" if all(decisions) else "branch_aware_refinement",
            decision_rule="All three assumed opacity scenarios must meet 1% pilot error; no fine truth used for selection",
            later_new_coordinate_check_max_relative_error=transfer_max,
            later_check_within_tolerance=transfer_max<=TOLERANCE))

    expected_hashes = {name: dense["config"]["code_sha256"][name] for name in CRITICAL_ENGINE}
    baseline_name = "analysis/rydberg_gold_map_20261103/resolved_parameters.json"
    expected_hashes[baseline_name] = dense["config"]["code_sha256"][baseline_name]
    mismatches = check_fingerprint(expected_hashes, REPO)
    domain = dict(temperature_c=[15.,35.], coupling_power_mw=[20.,40.], probe_power_uw=[2.,10.],
        lo_rabi_mhz=[2.,6.], probe_detuning_mhz=[-6.,6.], opacity_scenarios=[.0005,.001,.0015],
        cell_mm=50., beam_diameter_mm=.15, if_khz=40., delta_c_mhz=0., delta_rf_mhz=0.,
        temperature_policy="Linked", objective="Profiled PSN-only field ASD; zero RIN/electronics")
    baseline_parameters=json.loads((PILOT/'resolved_parameters.json').read_text(encoding='utf-8'))['baseline_parameters']
    common = dict(origin="⟦AI⟧ Saved-model resolution evidence",
        adoption_status="Analysis recommendation; not a human-registered goal",
        model_tag=MODEL_TAG,fixed_model_settings_sha256=fixed_settings_fingerprint(baseline_parameters),
        fixed_settings_hash_excludes=sorted(SWEEP_PARAMETER_KEYS),
        parameter_snapshot=baseline_parameters,
        reuse_check='reuse_rejections(registry,current_parameters,current_source_hashes,model_tag); eligibility is not experimental validation',
        domain=domain, source_campaign_signature=dense["config"]["signature"],
        expected_critical_source_sha256=expected_hashes,
        current_source_mismatches=mismatches,
        applicability_status="INVALIDATED_FOR_CURRENT_SOURCE" if mismatches else "ELIGIBLE_FOR_OLD_MODEL_ONLY",
        physical_validation_status="NOT_ESTABLISHED",
        invalidation_rules=["Any critical source/default/Rabi/opacity/noise-policy change",
            "New density- or temperature-dependent dephasing, separated temperature policy, spatial propagation or manifold",
            "Independent coupling detuning, RF detuning or IF outside fixed scope",
            "Changed geometry, beam convention, scan/control protocol or domain",
            "Guard error above declared decision tolerance"],
        recalibrated_model_policy="Do not inherit eligibility automatically; regenerate cross/guard audit under a new fingerprint")
    nominal = {a["doppler"]: a for a in audits if a["opacity"] == .001}
    entries = [
        dict(variable="temperature_c", role="Important nonlinear covariate; dominant density scale",
            reduction="Thermal mode: separate 1D slice/known N(T) baseline; retain sparse transit/Doppler checks",
            evidence="Density increases 4.7668x from20 to35C; thermal cross factorization passes; do not call T linear",
            policy="CONDITIONAL_REDUCTION", dense_step=1., proposed_seed_nodes=[15.,20.,25.,30.,33.,35.]),
        dict(variable="coupling_power_mw", role="Physical optical control",
            reduction="Thermal mode: 1D slice and guards. Off mode: retain branch-aware interaction refinement",
            evidence="Nominal retrospective cross errors: on %.6g%%; off %.6g%%" %
                (100*nominal['on']['retrospective_dense_cross_29_fit_160_check']['max_absolute_relative_error'],
                 100*nominal['off']['retrospective_dense_cross_29_fit_160_check']['max_absolute_relative_error']),
            policy="MODE_DEPENDENT_REDUCTION", dense_step=2.5, proposed_seed_nodes=[20.,30.,40.]),
        dict(variable="probe_power_uw", role="Physical control interacting with LO",
            reduction="Retain joint local scans and matched fixed-LO validation",
            evidence="All thermal profiled argmins at10uW do not imply monotonic fixed-LO response; fixedLO3.7 has turnover",
            policy="KEEP_INTERACTION", dense_step=.5, proposed_seed_nodes=[2.,4.,6.,8.,10.]),
        dict(variable="lo_rabi_mhz", role="Physical RF control interacting with probe power",
            reduction="Local response/robust plateau refinement; compare with fixedGoldLO3.7 separately",
            evidence="Thermal optima4.625/4.75MHz; near-equal controls exist and do not establish a filled acceptable rectangle",
            policy="KEEP_INTERACTION", dense_step=.125, proposed_seed_nodes=[2.,3.,3.7,4.5,6.]),
        dict(variable="probe_detuning_mhz", role="Signed finite-IF response coordinate",
            reduction="Thermal: center-first plus signed off-center challenges; off: track central and side candidates separately",
            evidence="Thermal189argmins at0; finiteIF sign symmetry rejected; static-slope-only candidates are insufficient",
            policy="ADAPTIVE_SIGNED_PEAK_SEARCH", dense_step=.01, proposed_seed_nodes=[-6.,-3.,-1.,0.,1.,3.,6.]),
        dict(variable="coupling_detuning_mhz", role="Missing independent one/two-photon detuning control in old model",
            reduction="No old-map predictability claim; audit after physical model supports it",
            evidence="Fixed0 in all old maps", policy="UNTESTED_NEW_VARIABLE", dense_step=None, proposed_seed_nodes=[]),
        dict(variable="cold_spot_temperature_c_vs_effective_temperature_c", role="Separate density and kinetic-temperature covariates",
            reduction="Keep known density analytic/postprocessed; audit kinetic/transit residual separately",
            evidence="Old maps link both temperatures; no independent two-temperature evidence",
            policy="UNTESTED_DECONFOUNDING", dense_step=None, proposed_seed_nodes=[]),
        dict(variable="effective_opacity_participation", role="Uncalibrated nuisance parameter",
            reduction="Scenario/calibration axis, never choose its value solely to minimize sensitivity",
            evidence="Three assumed scenarios are not statistical confidence limits; N*eta*L confounding at fixed atomic response",
            policy="CALIBRATE_NOT_OPTIMIZE", dense_step=None, proposed_seed_nodes=[.0005,.001,.0015]),
        dict(variable="beam_diameter_mm_and_overlap", role="Geometry/calibration, then possible design variable",
            reduction="Resolve beam definition and P-to-Rabi conversion before design optimization",
            evidence="Oldfixed0.15mm; 1/e2/FWHM convention unresolved; affects intensity, transit and overlap together",
            policy="UNTESTED_FIXED_GEOMETRY", dense_step=None, proposed_seed_nodes=[.15]),
        dict(variable="if_khz_and_rf_detuning", role="Application/measurement protocol",
            reduction="Include only for specified bandwidth/carrier/linearity question",
            evidence="Only40kHz andzeroRFdetuning evaluated", policy="UNTESTED_OUTSIDE_PROTOCOL", dense_step=None, proposed_seed_nodes=[40.]),
        dict(variable="dephasing_and_technical_noise", role="Measured/model nuisance parameters",
            reduction="Constrain through spectra/noise budget, then propagate uncertainty",
            evidence="ExtraT/densitydephasing0 andzeroRIN/electronics; changing these invalidates old reductions",
            policy="CALIBRATE_NOT_OPTIMIZE", dense_step=None, proposed_seed_nodes=[]),
    ]
    # These are paired-axis prediction errors: neither T nor Pc is physically
    # discarded, and no separate causal error estimate is inferred for one axis.
    for entry in entries[:2]:
        entry['prediction_evidence'] = [dict(doppler=a['doppler'],opacity=a['opacity'],
            pilot_check_max_relative_error=a['pilot_only_decision']['pilot_check']['max_absolute_relative_error'],
            sparse_pilot_to_new_dense_max_relative_error=a['pilot_to_new_dense_coordinates']['max_absolute_relative_error'],
            retrospective_dense_cross_max_relative_error=a['retrospective_dense_cross_29_fit_160_check']['max_absolute_relative_error'],
            evidence_file='predictability_audit.json',
            scope='Paired T/Pc scalar surface; not independent causal importance or a new-model validation') for a in audits]
    registry = dict(**common, decision_tolerance_relative=TOLERANCE,
        tolerance_scope="AI illustrative 1% scalar prediction/decision criterion, not experimental uncertainty",
        mode_policies_selected_from_pilot=mode_decisions, variables=entries,
        proposed_dense_domain_sparse_plan=guard_plan(dense['config']['temperature_c'],dense['config']['coupling_power_mw']),
        future_coarse_domain_seed_plan=guard_plan(pilot['config']['temperature_c'],pilot['config']['coupling_power_mw']),
        actual_historical_policy_fit_and_check=dict(fit_cross_coordinates=8,pilot_check_coordinates=10,
            subsequent_new_dense_check_coordinates=171,retrospective_dense_cross_fit=29,retrospective_dense_cross_check=160),
        calibrated_model_status="PENDING_NEW_PHYSICAL_VALIDATION_AND_RESOLUTION_AUDIT")
    transitions = []
    for mode in dense['config']['modes']:
        r = [x for x in dense['optima'] if x['doppler']==mode and x['opacity']==.001]
        b = [x for x in dense['branchgaps'] if x['doppler']==mode and x['opacity']==.001]
        transitions.extend(dict(doppler=mode,opacity=.001, **row) for row in branch_transitions(r,b))
    with (DENSE/'control_points.csv').open(encoding='utf-8-sig',newline='') as stream:
        saved_controls=list(csv.DictReader(stream))
    fixed_protocol=[]
    for mode in dense['config']['modes']:
        selected=[row for row in saved_controls if row['doppler']==mode and float(row['opacity'])==.001
                  and float(row['coupling_power_mw'])==30 and float(row['probe_power_uw'])==6 and float(row['lo_rabi_mhz'])==3.7]
        normalized={float(row['temperature_c']):float(row[FIELD])*densities[float(row['temperature_c'])]
                    *np.sqrt(float(row['transmission'])) for row in selected}
        fixed_protocol.append(dict(doppler=mode,probe_power_uw=6.,lo_rabi_mhz=3.7,coupling_power_mw=30.,
            residual_definition='S(T)*N(T)*sqrt(transmission at selected frequency), normalized at20C',
            residual_35_to_20_ratio=normalized[35.]/normalized[20.],
            temperature_nodes=[dict(temperature_c=t,residual_normalized_at_20=normalized[t]/normalized[20.]) for t in sorted(normalized)],
            scope='Fixed Pp/LO/Pc but probe detuning retains saved local optimization; not a pure fixed-detuning temperature experiment'))
    report = dict(origin="⟦AI⟧ Resolution/predictability audit; no new physical solves",
        validation_scope="Historical scalar surfaces only; old profile objective, not independent experiment",
        criterion_relative=TOLERANCE, source_campaign_signature=dense['config']['signature'],
        scenarios=audits, axis_resolution_diagnostics=axis_rows,
        separately_profiled_branch_audits=branches,fixed_protocol_density_loss_residual=fixed_protocol,
        mode_decisions_from_pilot=mode_decisions, nominal_argmin_transitions=transitions,
        source_sha256={str(path.relative_to(REPO)).replace('\\','/'):digest(path)
                      for path in [PILOT/'results.json',DENSE/'results.json',DENSE/'interpretation.json',DENSE/'control_points.csv',Path(__file__)]})
    ROOT.mkdir(parents=True,exist_ok=True)
    write_json(ROOT/'predictability_audit.json',report)
    write_json(ROOT/'variable_resolution_registry.json',registry)
    write_csv(ROOT/'heldout_predictions.csv',prediction_rows)
    write_csv(ROOT/'axis_resolution_diagnostics.csv',axis_rows)
    write_csv(ROOT/'argmin_transition_diagnostics.csv',transitions)
    write_csv(ROOT/'variable_resolution_registry.csv',[dict(origin='⟦AI⟧',variable=x['variable'],role=x['role'],policy=x['policy'],
        recommended_reduction=x['reduction'],evidence=x['evidence'],scope_status=registry['applicability_status'],
        physical_validation_status=registry['physical_validation_status'],source_model_signature=dense['config']['signature'],
        evidence_file='predictability_audit.json',tolerance_relative=TOLERANCE) for x in entries])
    write_csv(ROOT/'cross_slice_audit.csv',[dict(doppler=a['doppler'],opacity=a['opacity'],
        pilot_max_error=a['pilot_only_decision']['pilot_check']['max_absolute_relative_error'],
        pilot_policy=a['pilot_only_decision']['selected_action'],
        pilot_to_new_dense_max_error=a['pilot_to_new_dense_coordinates']['max_absolute_relative_error'],
        dense_cross_max_error=a['retrospective_dense_cross_29_fit_160_check']['max_absolute_relative_error'],
        tolerance_relative=TOLERANCE,source_model_signature=dense['config']['signature'],
        scope='Saved old simulation only; not experimental or recalibrated model validation') for a in audits])
    print(json.dumps(dict(mode_decisions=mode_decisions,source_status=registry['applicability_status'],
        transfer_errors=[dict(mode=a['doppler'],opacity=a['opacity'],**a['pilot_to_new_dense_coordinates']) for a in audits],
        registry_variables=len(entries)),indent=2))


if __name__ == '__main__':
    main()
