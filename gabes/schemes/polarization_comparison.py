"""Comparison views of the existing magneto-optics and PSR readouts.

No atomic solve, noise completion or peak fitting is reimplemented here. Each
series comes from the scheme's ordinary headless readout on explicit conditions.
"""
import numpy as np

from .. import __version__
from ..model_comparison import condition_views, sample_values, same_grid


SCHEMA = "gabes-spectral-comparison-v1"


def _number(value):
    return float(value) if value is not None and np.isfinite(value) else None


def _reading(key, label, value, unit="", *, status=None, reason=None):
    result = dict(key=key, label=label, value=_number(value), unit=unit)
    if status:
        result["status"] = status
    if reason:
        result["reason"] = reason
    return result


def _panel(key, title, y_label, unit, series, observable, *, reference=None):
    result = dict(key=key, title=title, y_label=y_label, unit=unit,
                  curves=[dict(series_id=s["id"], observable=observable, label=s["label"]) for s in series])
    if reference is not None:
        result["reference"] = reference
        if unit == "dB":
            result["reference_label"] = "0 dB · shot noise"
    return result


def _source(scheme, method, **extra):
    return dict(model=type(scheme).__name__, method=method, package_version=__version__,
                cache_version=scheme.cache_version, provenance="Numerical model",
                interpolation="none; original scheme samples", missing_value="null", **extra)


def _atomic_conditions(conditions):
    """Keep the native complex Jones vector explicit in a JSON-safe record."""
    if conditions is None:
        return None
    result = dict(conditions)
    if "jones_xy" in result:
        result["jones_xy"] = [dict(real=float(z.real), imag=float(z.imag))
                              for z in np.asarray(result["jones_xy"], complex)]
    return result


def magneto_comparison(scheme, params, *, resolve_raw, current_raw=None, current_view=None):
    from .magneto_unified import (
        FIELD_SCAN, ROTATION, POLARIZATION_COMPARISON, TRANSITION_COMPARISON,
        signal_for_axis,
    )

    params = scheme.normalize_params({**scheme.defaults(), **params})
    mode = params["comparison_mode"]
    if mode == POLARIZATION_COMPARISON:
        varied_key, varied_label, varied_unit = "qwp_deg", "QWP angle", "°"
        values = sorted({0., 9., 18., 27., 36., 45., float(params[varied_key])})
    elif mode == TRANSITION_COMPARISON:
        varied_key, varied_label, varied_unit = "line", "D-line branch", ""
        values = ["D1", "D2"]
    else:
        raise ValueError("unsupported magneto-optics comparison mode")
    expected_signal = signal_for_axis(params["scan_axis"])
    if current_view is not None and current_view["magneto_readout"].get("signal_type") != expected_signal:
        # A legacy cached view may carry an independent Signal selection.
        current_view = None
    show_noise = bool(params["show_shot_noise"])
    result_series, first = [], None
    for condition, raw, view in condition_views(
            scheme, params, varied_key, values, resolve_raw=resolve_raw,
            current_raw=current_raw, current_view=current_view):
        readout = view["magneto_readout"]
        x = np.asarray(readout["x"], float)
        if first is None:
            first = readout
        # Every curve retains its original coordinate array. Native branch
        # grids may differ; comparison never resamples an atomic solution.
        op = readout["operating_index"]
        optical = raw["scan_axis"] != FIELD_SCAN
        picked = readout["optimum_index"] if optical else op
        selected_budget = (readout["rotation_sensitivity"] if
                           signal_for_axis(raw["scan_axis"]) == ROTATION else
                           readout["transmission_sensitivity"])
        sample = op if picked is None else picked
        asd = selected_budget["field_asd_t_per_sqrt_hz"][sample]
        status = selected_budget["status"][sample]
        T = np.asarray(readout["transmission"], float)
        wing = float(.5*(T[0]+T[-1]))
        reference_label = "optical reference" if optical else "command zero"
        readouts = [
            _reading("reference_rotation", f"Rotation at {reference_label}", readout["rotation_rad"][op], "rad",
                     status=readout["rotation_status"]),
            _reading("reference_transmission", f"Transmission at {reference_label}", T[op], "1"),
            _reading("wing_reference", "Mean transmission at scan endpoints", wing, "1"),
            _reading("center_minus_wings", f"Transmission at {reference_label} minus endpoint mean", T[op]-wing, "1",
                     reason="A signed sampled contrast; check the physical field and scan extent."),
        ]
        if show_noise:
            readouts.append(_reading(
                "photon_field_asd", "Best sampled photon shot-noise field limit" if optical else
                "Photon shot-noise field limit at command zero", None if asd is None else asd*1e15,
                "fT/√Hz", status=status,
                reason="Selected signal's photon shot-noise limit; atomic and technical noise are excluded."))
            if optical:
                readouts.extend([
                    _reading("sensitivity_optical_detuning", "Optical detuning at sampled photon-limit minimum",
                             readout["optimum_detuning_mhz"], "MHz",
                             status="available" if picked is not None else "unavailable"),
                    dict(key="sensitivity_optimum_status", label="Optical sensitivity search scope",
                         value=readout["optimum_status"], unit=""),
                ])
        value = condition[varied_key]
        if mode == POLARIZATION_COMPARISON:
            label = f"QWP {value:g}°" + (" (linear)" if value == 0 else " (circular)" if value == 45 else "")
            series_id = f"qwp-{value:g}"
        else:
            label = f"{condition['line']} · F={readout['ground_f']:g} centroid"
            series_id = condition["line"].lower()
        result_series.append(dict(
            id=series_id, label=label,
            parameter=dict(key=varied_key, label=varied_label, value=value, unit=varied_unit),
            x=sample_values(x),
            observables=dict(rotation_mrad=sample_values(np.asarray(readout["rotation_rad"], float)*1000),
                             transmission_pct=sample_values(T*100),
                             photon_field_asd_ft_per_sqrt_hz=sample_values(
                                 np.asarray(selected_budget["field_asd_t_per_sqrt_hz"], float)*1e15)),
            readouts=readouts,
            observable_status=dict(rotation_mrad=[readout["rotation_status"] for _ in x]),
            scope=dict(rotation=readout["rotation_status"], model=readout["model_scope"],
                       cell=readout["cell_scope"], propagation=readout["propagation_status"],
                       velocity_quadrature=readout["velocity_quadrature"],
                       diagnostics=readout["diagnostics"], sensitivity=selected_budget,
                       sensitivity_selection=dict(
                           coordinate=None if picked is None else float(x[picked]),
                           unit=readout["x_unit"], status=readout["optimum_status"],
                           at_scan_edge=readout["optimum_at_scan_edge"])),
            source=_source(scheme, "compute + native magneto_readout",
                           isotope=readout["isotope"], line=readout["line"],
                           ground_f=readout["ground_f"],
                           reference_centroid_mhz=readout["reference_centroid_mhz"],
                           detuning_reference=readout["detuning_reference"],
                           conditions=readout["conditions"]),
        ))
    is_field = first["scan_axis"] == FIELD_SCAN
    bias = float(first["physical_bias_ut"])
    reference = ("Commanded longitudinal field; Bz = Bcommand + longitudinal bias" if is_field else
                 "Each branch's selected ground-F line-strength centroid")
    panels = [
        _panel("rotation", "Polarization rotation", "Rotation", "mrad", result_series, "rotation_mrad", reference=0),
        _panel("transmission", "Transmission at the same coordinate", "Transmission", "%", result_series, "transmission_pct"),
    ]
    panels[0]["note"] = (
        "Rotation is defined only for linear input in this model. Elliptical and circular input curves are unavailable; compare their transmission.")
    if mode == POLARIZATION_COMPARISON:
        panels[1].update(stackable=True, wide=True, stack_label="QWP angle",
                         stack_scale_unit="percentage points T",
                         title="Transmission across incident polarizations")
    if expected_signal != ROTATION:
        panels.reverse()
    if show_noise:
        panel = _panel("photon_field_asd", "Photon shot-noise magnetic sensitivity", "Field ASD", "fT/√Hz",
                       result_series, "photon_field_asd_ft_per_sqrt_hz")
        panel["scale"] = "log"
        panels.append(panel)
    return dict(
        schema=SCHEMA, scheme=scheme.name, mode=first["scan_axis"], comparison_mode=mode,
        primary_panel="rotation" if expected_signal == ROTATION else "transmission",
        axis=dict(label="Commanded longitudinal field" if is_field else "Optical detuning",
                  unit=first["x_unit"], reference=reference),
        series=result_series, panels=panels,
        fixed={k:v for k,v in params.items() if k != varied_key},
        varied=dict(key=varied_key, label=varied_label, unit=varied_unit, values=values),
        markers=([dict(x=0., label="Command B = 0", display_label="Bcmd = 0"),
                  dict(x=-bias, label="Longitudinal component Bz = 0", display_label="Bz = 0")] if is_field else []),
        source=_source(scheme, "compute + magneto_readout", physical_bias_ut=bias,
                       residual_transverse_b_ut=float(params["residual_transverse_b_ut"]),
                       command_b_ut=first["command_b_ut"], physical_b_ut=first["physical_b_ut"],
                       detuning_reference=first["detuning_reference"],
                       comparison_reference=dict(
                           citation="PRA 81, 023416 (2010), Figs. 4–6",
                           doi="10.1103/PhysRevA.81.023416",
                           url="https://doi.org/10.1103/PhysRevA.81.023416",
                           use="Qualitative polarization-dependent transmission comparison; not a numerical reproduction.")),
        notes=[
            "QWP angle varies incident ellipticity from linear (0°) to circular (45°); it is not a linear-polarization azimuth.",
            "Rotation for elliptical/circular input is unsupported in this model and remains unavailable; use the calculated transmission curves for the angle comparison.",
            "D1/D2 comparison uses each complete hyperfine manifold at its ground-F centroid reference. It does not isolate an F→F′ transition or reproduce the paper's conditions.",
            "PRA 81, 023416 (2010), Figs. 4–6 motivates the polarization comparison. Figure-specific transition assignments and quantitative agreement are not established here.",
            "Rotation and transmission use one atomic state at each field or detuning; their vertical units differ.",
            "A transmission enhancement or absorption enhancement is read from the calculated curve, not assumed from a preset name.",
            "At Bcommand = −bias only the longitudinal component is zero. A residual transverse field prevents zero total vector field.",
            "Displayed field sensitivity is a photon shot-noise limit. The quasistatic finite-amplitude screen and atomic bandwidth scope remain in the exported readout.",
        ],
    )


def psr_comparison(scheme, params, *, resolve_raw, current_raw=None, current_view=None):
    from .psr import CONDITIONAL, FULL_OBE
    from .psr_rf import RF_MODELS

    mode = params["calculation"]
    if mode == CONDITIONAL:
        return _conditional_psr(scheme, params, resolve_raw=resolve_raw,
                                current_raw=current_raw, current_view=current_view)
    if mode == FULL_OBE:
        return _full_psr(scheme, params, resolve_raw=resolve_raw,
                        current_raw=current_raw, current_view=current_view)
    if mode in RF_MODELS:
        raw = current_raw if current_raw is not None else resolve_raw(params)
        view = current_view if current_view is not None else scheme.headless_observables(raw, params)
        return _rf_psr(scheme, params, raw, view)
    raise ValueError(f"Unknown PSR calculation {mode!r}")


def _conditional_psr(scheme, params, *, resolve_raw, current_raw=None, current_view=None):
    from .psr import PHASE_VIEW

    varied_key = "excess_noise_per_od"
    values = sorted({0., 1., 2., float(params[varied_key])})
    # This parameter is readout-only. All noise assumptions share the same raw
    # atomic response, including a cached direct solve at the selected detuning.
    raw = current_raw if current_raw is not None else resolve_raw(params)
    phase_view = params["view"] == PHASE_VIEW
    result_series, expected = [], None
    for excess in values:
        condition = {**params, varied_key:excess}
        view = (current_view if excess == params[varied_key] and current_view is not None
                else scheme.headless_observables(raw, condition))
        readout = view["psr_readout"]
        point = readout["selected_point"]
        phase, spectrum = readout["phase_scan"], readout["spectrum"]
        if phase_view:
            x = np.asarray(phase["phase_degrees"] if phase is not None else np.linspace(0.,180.,361), float)
            y1 = phase["noise_a_db"] if phase is not None else np.full(x.shape,np.nan)
            y2 = phase["noise_b_db"] if phase is not None else np.full(x.shape,np.nan)
            available = np.full(x.shape, point["available"], bool)
            applicable = np.full(x.shape, point["approximation_valid"], bool)
            reasons = list(point["approximation_reasons"])
            if not point["available"]:
                reasons.append(point["status"])
            sample_reasons = [list(reasons) for _ in x]
        else:
            x = np.asarray(spectrum["detunings_mhz"], float)
            y1, y2 = spectrum["minimum_db"], spectrum["maximum_db"]
            available = np.asarray(spectrum["available"], bool)
            applicable = np.asarray(spectrum["approximation_valid"], bool)
            sample_reasons = [list(dict.fromkeys(
                reason for key in ("numerical_reasons", "passive_noise_reasons", "approximation_reasons")
                for reason in raw[key][index])) for index in range(x.size)]
        sample_status = [
            "unavailable" if not valid else
            "available" if inside else "available outside approximation screens"
            for valid, inside in zip(available, applicable)
        ]
        if expected is None:
            expected = x
        else:
            same_grid(expected,x)
        result_series.append(dict(
            id=f"excess-{excess:g}", label="nₑₓ = "+f"{excess:g}"+(" (omitted)" if excess == 0 else ""),
            parameter=dict(key=varied_key,label="Assumed excess noise per OD",value=excess,unit="1"),
            x=sample_values(x), observables=dict(noise_a_db=sample_values(y1),noise_b_db=sample_values(y2)),
            sample_status=sample_status, sample_reasons=sample_reasons,
            readouts=[
                _reading("minimum_noise", "Minimum conditional noise", point["minimum_db"], "dB",
                         status=point["status"]),
                _reading("minimum_phase", "Minimum-noise phase", point["min_angle_deg"], "°",
                         status="available" if point["available"] else "unavailable"),
                _reading("transmission", "Transmission at selected optical detuning", point["transmission"],"1"),
            ],
            scope=dict(noise_assumption=readout["noise_assumption"], model=readout["model_scope"],
                       approximation_thresholds=readout["approximation_thresholds"],
                       approximation_reasons=point["approximation_reasons"],
                       selection=readout["selection_scope"], atomic_diagnostics=point["atomic_diagnostics"]),
        ))
    return dict(
        schema=SCHEMA, scheme=scheme.name, mode=params["calculation"],
        axis=dict(label="Analysis phase" if phase_view else "Optical detuning",unit="°" if phase_view else "MHz",
                  reference="Relative to bright carrier amplitude/ellipticity axes" if phase_view else raw["detuning_reference"]),
        series=result_series, panels=[
            _panel("noise_a", "Stokes A" if phase_view else "Minimum noise", "Noise / shot noise","dB",
                   result_series,"noise_a_db",reference=0),
            _panel("noise_b", "Orthogonal Stokes B" if phase_view else "Maximum noise","Noise / shot noise","dB",
                   result_series,"noise_b_db",reference=0),
        ],
        fixed={k:v for k,v in params.items() if k != varied_key},
        varied=dict(key=varied_key,label="Assumed excess noise per OD",unit="1",values=values),
        source=_source(scheme,"shared atomic response + headless conditional PSR readout",
                       optical_detuning_mhz=float(params["detuning_mhz"]),
                       measurement_conventions=readout["measurement_conventions"]),
        notes=[
            "These are conditional Gaussian-noise estimates. The excess diffusion is an assumption, not independently predicted atomic noise.",
            "nₑₓ = 0 omits the excess-noise term. It does not establish noiseless atoms.",
            "Minimum noise and its phase belong to the same assumption. 0 dB is shot noise for this polarization quadrature.",
            "Changing phase, cursor or excess-noise assumption uses the same optical response. An unavailable response is a gap, never a 0 dB trace.",
            "Applicability and reasons are retained for every plotted sample. Estimates outside approximation screens retain their numerical values.",
        ],
    )


def _full_psr(scheme, params, *, resolve_raw, current_raw=None, current_view=None):
    key="obe_intensity_mw_cm2"
    current=float(params[key])
    values=sorted({max(0.,current*.5),current,min(100000.,current*2.)} if current else {0.,.05,.1})
    series=[]
    for condition,raw,view in condition_views(
            scheme,params,key,values,resolve_raw=resolve_raw,current_raw=current_raw,current_view=current_view):
        ground,excited=raw["ground_levels"],raw["excited_levels"]
        labels=[f"5S1/2 F={f:g} mF={m:g}" for f,m in ground]+[f"5P1/2 F={f:g} mF={m:g}" for f,m in excited]
        population=np.r_[raw["ground_populations"],raw["excited_populations"]]
        intensity=float(condition[key])
        series.append(dict(
            id=f"intensity-{intensity:g}",label=f"{intensity:g} mW/cm²",
            parameter=dict(key=key,label="Local pump intensity",value=intensity,unit="mW/cm²"),
            x=list(range(len(population))),observables=dict(population=sample_values(population)),
            readouts=[_reading("absorption","Local absorption coefficient",raw["alpha_per_m"],"m⁻¹"),
                      _reading("excited_fraction","Total excited-state population",raw["excited_fraction"],"1")],
            scope=dict(model=raw["model_id"],diagnostics=raw["diagnostics"],conditions=_atomic_conditions(raw["conditions"])),
        ))
    return dict(
        schema=SCHEMA,scheme=scheme.name,mode=params["calculation"],
        axis=dict(label="Atomic sublevel",unit="index",reference="⁸⁷Rb D1; ordered ground then excited sublevels",
                  tick_labels=labels),
        series=series,panels=[_panel("population","All sixteen atomic populations","Fraction of all atoms","1",
                                     series,"population",reference=0)],
        fixed={k:v for k,v in params.items() if k != key},
        varied=dict(key=key,label="Local pump intensity",unit="mW/cm²",values=values),
        source=_source(scheme,"solve_psr_pump + native ground/excited populations",sublevel_labels=labels),
        notes=["Each curve is the existing sixteen-sublevel steady-state solution at one pump intensity.",
               "All populations use the same detuning, axial velocity and polarization. Local absorption is reported beside each solution.",
               "The x coordinates identify atomic sublevels, not a continuous frequency axis."],
    )


def _rf_psr(scheme, params, raw, view):
    from .psr_rf import PHASE_VIEW
    readout=view["psr_readout"]
    point,phase,spectrum=readout["selected_point"],readout["phase_scan"],readout["spectrum"]
    is_phase=params["rf_view"]==PHASE_VIEW
    if is_phase:
        x=np.asarray(phase["phase_degrees"] if phase is not None else np.linspace(0.,180.,361),float)
        values={key:sample_values(phase[key] if phase is not None else np.full(x.shape,np.nan))
                for key in ("noise_a_db","noise_b_db")}
    else:
        x=np.asarray(spectrum["frequencies_hz"],float)/1e6 if spectrum is not None else raw["center_hz"]/1e6*np.logspace(-1,1,25)
        values={key:sample_values(spectrum[key] if spectrum is not None else np.full(x.shape,np.nan))
                for key in ("minimum_db","maximum_db","min_angle_deg")}
    series=[dict(id="rf",label="Current optical conditions",x=sample_values(x),observables=values,
                 parameter=dict(key="rf_frequency_mhz",label="Selected RF",value=params["rf_frequency_mhz"],unit="MHz"),
                 readouts=[_reading("minimum_noise","Minimum RF noise",point["minimum_db"],"dB",
                                    status="available" if point["available"] else "unavailable",reason=point["reason"]),
                           _reading("maximum_noise","Maximum RF noise",point["maximum_db"],"dB"),
                           _reading("minimum_phase","Minimum-noise phase",point["min_angle_deg"],"°")],
                 scope=dict(model=readout["model_id"],available=point["available"],reason=point["reason"],
                            conditions=_atomic_conditions(raw.get("conditions")),measurement_conventions=readout["measurement_conventions"]))]
    if is_phase:
        panels=[dict(key="noise",title="Polarization noise at selected RF",y_label="Noise / shot noise",unit="dB",reference=0,reference_label="0 dB · shot noise",
                     curves=[dict(series_id="rf",observable=k,label=label) for k,label in
                             (("noise_a_db","Stokes A"),("noise_b_db","Orthogonal Stokes B"))])]
    else:
        panels=[dict(key="noise",title="RF polarization noise",y_label="Noise / shot noise",unit="dB",reference=0,reference_label="0 dB · shot noise",
                     curves=[dict(series_id="rf",observable=k,label=label) for k,label in
                             (("minimum_db","Minimum"),("maximum_db","Maximum"))]),
                _panel("phase","Optimal detection phase","Phase","°",series,"min_angle_deg")]
    return dict(
        schema=SCHEMA,scheme=scheme.name,mode=params["calculation"],
        axis=dict(label="Homodyne phase" if is_phase else "Analysis frequency",unit="°" if is_phase else "MHz",
                  reference="Relative to mean output carrier" if is_phase else "Positive RF sideband frequency",
                  scale="linear" if is_phase else "log"),
        series=series,panels=panels,fixed=dict(params),varied=dict(key=None,label="Existing RF readout",values=[]),
        source=_source(scheme,"native PSR RF headless readout",selected_frequency_hz=raw["center_hz"],
                       optical_scope=("Longitudinal pump intensity and Jones evolution with local stationary atoms"
                                      if raw.get("propagation_convergence") is not None else
                                      "Independent Maxwell velocity classes, uniform undepleted pump"
                                      if raw.get("thermal_convergence") is not None else
                                      "Single axial velocity, uniform undepleted pump"),
                       thermal_convergence=raw.get("thermal_convergence"),
                       propagation_convergence=raw.get("propagation_convergence"),
                       propagation_velocity_convergence=raw.get("propagation_velocity_convergence"),
                       unavailable_reason=point["reason"]),
        notes=["RF analysis frequency differs from optical detuning; the atomic and propagation models remain the existing PSR RF implementation.",
               "The optimum phase is plotted separately from noise because degrees and dB have different units.",
               "Atomic/propagation or homodyne checks can make a result unavailable. Such samples have no numerical noise curve."],
    )
