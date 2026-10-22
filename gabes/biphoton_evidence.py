"""Readout provenance for reduced SFWM models; no source or detector arithmetic."""


def biphoton_output_evidence(raw, stats, params, *, velocity_fwhm_tolerance=0.03):
    """Describe each observable without promoting an anchor to a measurement.

    The velocity refinement checks source FWHM only. It does not certify the
    complete waveform, a Fourier bandwidth, temporal support, or absolute rates.
    A failed check leaves affected numerical arrays available as diagnostics but
    withholds quantitative display values. No inputs or arrays are mutated.
    """
    topology = raw["topology"]
    predictive = bool(raw.get("predictive", False))
    converged = bool(raw.get("velocity_converged", False))
    reference = any(getattr(topology, key, None) is not None for key in
                    ("reference_od", "reference_fwhm_ns", "target_g2_peak"))
    rate_origin = "reference_anchored" if reference else "assumed_template"
    waveform_origin = ("reduced_model" if predictive else "reference_anchored"
                       if getattr(topology, "reference_fwhm_ns", None) is not None
                       else "assumed_template")
    waveform_display = not predictive or converged
    filter_active = (params.get("filter_bandwidth_mhz") is not None
                     and float(params.get("filter_bandwidth_mhz", 0.0)) > 0.0)
    filtered_display = waveform_display or not filter_active
    added_background = (not predictive
                        and getattr(topology, "target_g2_peak", None) is not None)
    background_applied = (
        float(stats["raw_g2_peak"]) > float(stats["g2_peak"])
        if added_background and "raw_g2_peak" in stats and "g2_peak" in stats
        else None if added_background else False)
    calibration = ("reference_added_background" if added_background
                   else "assumed_detector_background")
    numerical_scope = {
        "velocity_check": ("source_fwhm_relative_change" if predictive
                           else "not_assessed_reference_model"),
        "velocity_fwhm_status": ("passed" if converged else "failed_or_capped")
        if predictive else "not_assessed",
        "velocity_fwhm_relative_tolerance": velocity_fwhm_tolerance if predictive else None,
        "whole_waveform_certified": False,
        "time_grid_and_support_certified": False,
        "spectral_bandwidth_certified": False,
        "absolute_rate_certified": False,
    }
    common = ["finite supplied source-delay domain",
              "Gaussian timing-difference response",
              "scalar filter throughput without temporal reshaping"]
    outputs = {}

    def record(names, origin, dependencies=(), assumptions=(), display=True, scope="model_only"):
        for name in names:
            outputs[name] = {
                "origin": origin,
                "status": "qualitative_unconverged" if not display else (
                    "assumed_input" if origin == "assumed_input" else
                    "assumed_template" if not reference or origin == "assumed_template"
                    else "conditional_estimate"),
                "dependencies": list(dependencies),
                "assumptions": list(assumptions),
                "independent_validation": False,
                "display_numeric": bool(display),
                "numerical_scope": scope,
            }

    record(("source_v",), "reduced_model",
           assumptions=("analytic velocity-class source envelope from drives, detunings and linewidths",),
           display=waveform_display, scope=numerical_scope["velocity_check"])
    record(("psi_tau", "source_fwhm_ns"), waveform_origin,
           ("source_v",) if predictive else ("source_v", "emission_decay_ns"),
           ("velocity refinement checks FWHM only; no full waveform certificate",)
           if predictive else ("analytic velocity sum multiplied by the stored emission-lifetime envelope; timing anchor is not independent validation",),
           display=waveform_display, scope=numerical_scope["velocity_check"])
    if predictive:
        record(("source_bandwidth_mhz",), "reduced_model", ("psi_tau",),
               ("no time-grid/support or spectral convergence certificate",),
               display=waveform_display, scope="spectral_convergence_not_assessed")
    else:
        bandwidth_reference = getattr(topology, "reference_bandwidth_mhz", None)
        record(("source_bandwidth_mhz",),
               "reference_anchored" if bandwidth_reference else "assumed_input",
               ("reference_bandwidth_mhz",) if bandwidth_reference else (),
               ("stored bandwidth input, not derived from the waveform",)
               if bandwidth_reference else ("assumed 300 MHz fallback; topology has no reference bandwidth",),
               scope="stored_input_not_spectral_prediction")
    record(("source_pair_rate_cps",), rate_origin,
           assumptions=("pump/coupling/density/phase-match scaling; no derived absolute collection rate",))
    record(("phase_match_weight",), "reference_geometry_constructed",
           assumptions=("reference longitudinal wavevector offset; geometric sinc-squared model",))
    record(("od_estimate",), "reference_anchored" if getattr(
        topology, "reference_od", None) is not None else "reduced_model",
        assumptions=("OD reference scaling or reduced optical-depth estimate",))
    record(("filter_transmission",), "assumed_detector_model", ("source_bandwidth_mhz",),
           (common[2],), display=filtered_display)
    record(("pair_rate_cps", "singles_signal_cps", "singles_idler_cps",
            "all_time_coincidence_cps"), "detector_derived",
           ("source_pair_rate_cps", "filter_transmission"),
           ("detector efficiencies and dark/background rates are supplied inputs",),
           display=filtered_display)
    record(("raw_accidental_cps",), "detector_derived",
           ("singles_signal_cps", "singles_idler_cps", "coincidence_window_ns"),
           ("independent stationary singles multiplied by the coincidence-window width",),
           display=filtered_display)
    record(("detected_fwhm_ns", "fwhm_ns", "detected_delay_density_per_ns",
            "gate_capture_fraction"), "detector_convolution", ("psi_tau",), common,
           display=waveform_display, scope="conditional_on_source_shape")
    record(("coincidence_cps", "gate_captured_coincidence_cps",
            "heralding_signal", "heralding_idler"), "detector_derived",
           ("all_time_coincidence_cps", "gate_capture_fraction"), common,
           display=waveform_display)
    record(("accidental_cps", "added_accidental_cps"), calibration,
           ("raw_accidental_cps", "target_g2_peak", "psi_tau")
           if added_background else ("raw_accidental_cps",),
           ("target g2 bounds the peak by adding background; the target need not be reached" if added_background
            else "independent stationary singles/background model",),
           display=waveform_display if added_background else filtered_display)
    record(("g2_SI_tau", "g2_peak"),
           "reference_anchored" if added_background else "detector_derived",
           ("psi_tau", "all_time_coincidence_cps", "singles_signal_cps", "singles_idler_cps")
           + (("target_g2_peak",) if added_background else ()),
           common + ["histogram correlation is independent of coincidence-window width"],
           display=waveform_display, scope="conditional_on_source_shape")
    record(("CAR",), "reference_anchored" if added_background else "detector_derived",
           ("gate_captured_coincidence_cps", "accidental_cps"), common,
           display=waveform_display, scope="conditional_on_source_shape")
    record(("raw_g2_peak",), "detector_derived",
           ("psi_tau", "all_time_coincidence_cps", "singles_signal_cps", "singles_idler_cps"),
           common + ["before target-g2 added-background calibration; independent of coincidence-window width"],
           display=waveform_display, scope="conditional_on_source_shape")
    record(("raw_CAR",), "detector_derived",
           ("gate_captured_coincidence_cps", "raw_accidental_cps"),
           common + ["before target-g2 added-background calibration"],
           display=waveform_display, scope="conditional_on_source_shape")
    record(("cauchy_schwarz_R", "thermal_arm_comparison_R"), "assumed_thermal_arms",
           ("g2_peak",), ("g2_SS=g2_II=2 is assumed, not computed or measured",),
           display=waveform_display, scope="comparison_only_not_nonclassicality_evidence")
    # The model has no measured/full auto-correlations or quantum covariance.
    record(("physical_cauchy_schwarz_R", "physical_nonclassicality"), "unavailable",
           assumptions=("requires measured/full auto-correlations or an explicit ideal-canonical fixture",),
           display=False, scope="not_available")
    for name in ("physical_cauchy_schwarz_R", "physical_nonclassicality"):
        outputs[name]["status"] = "unavailable"
    record(("energy_mismatch_hz",), "assembly_invariant",
           assumptions=("wavelength bookkeeping is not experimental validation",))
    return {
        "schema_version": "biphoton-output-evidence-v1",
        "source_model": (("reduced_model" if predictive else "reference_model")
                         if reference else "assumed_template_model"),
        "topology": topology.name,
        "independent_validation": False,
        "measured_data_available": False,
        "seeded_fwm_validation_applicable": False,
        "reference_g2_calibration": {
            "requested": added_background,
            "applied": background_applied,
            "target_g2_peak": getattr(topology, "target_g2_peak", None)
            if added_background else None,
        },
        "numerical_scope": numerical_scope,
        "claims": {
            "physical_nonclassicality_available": False,
            "physical_squeezing_available": False,
            "absolute_rate_prediction_available": False,
            "reason": "No measured/full covariance, same-condition auto-correlations or held-out source data.",
        },
        "observables": outputs,
    }
