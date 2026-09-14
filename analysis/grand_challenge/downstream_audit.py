"""Independent RF-band and correlated-input audit, parallel to optical-mode work.

python -m analysis.grand_challenge.downstream_audit --output NEW.json --plot NEW.png
"""

import argparse
from dataclasses import asdict, replace
import hashlib
import json
from pathlib import Path
import platform

import numpy as np
import scipy

from gabes.fwm_quantum.inputs import ReducedPowerInputs, power_normalized_readout
from gabes.quantum.contracts import AnalysisFrequencyAxis, ParameterEvidence
from gabes.quantum.readout import DetectorResponse
from gabes.quantum.spectrum_analysis import SpectralTrace, analyze_squeezing_bands
from gabes.quantum.uncertainty import (
    JointInputUncertainty, propagate_ensemble, propagate_first_order, sample_inputs,
)


ROOT = Path(__file__).resolve().parents[2]


def _jsonable(value):
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, dict):
        return {key: _jsonable(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [_jsonable(item) for item in value]
    return value


def _evidence(key, value, unit, sigma):
    return ParameterEvidence(key, value, unit, "assumed", "downstream analytic fixture",
                             sigma, "declared diagnostic input", (), (),
                             "conditional software check; not measured apparatus uncertainty")


def _trace(frequency, ratio, *, layer="source"):
    return SpectralTrace(AnalysisFrequencyAxis.from_hz(frequency), ratio, layer=layer,
                         provenance="analytic valley fixture; not an atomic spectrum")


def _band_dict(analysis):
    result = _jsonable(asdict(analysis))
    for record, band in zip(result["bands"], analysis.bands):
        record["width_hz"] = band.width_hz
        record["observed_span_hz"] = band.observed_span_hz
    return result


def analytic_band_controls():
    center, scale, floor = 2.1e6, 1.3e6, .25
    root = np.sqrt(1-floor)
    exact_edges = np.array([center-scale*root, center+scale*root])
    cases = []
    for count in (81, 321, 1281):
        frequency = np.linspace(.1e6, 4e6, count)
        ratio = floor+((frequency-center)/scale)**2
        bands = analyze_squeezing_bands(_trace(frequency, ratio))
        band = bands.bands[0]
        edges = np.array([band.lower.frequency_hz, band.upper.frequency_hz])
        cases.append({"samples": count, "maximum_edge_error_hz": float(max(abs(edges-exact_edges))),
                      "analysis": _band_dict(bands)})
    coarse_frequency = np.linspace(.1e6, 4e6, 81)
    refined = analyze_squeezing_bands(_trace(coarse_frequency,
        floor+((coarse_frequency-center)/scale)**2),
        crossing_evaluator=lambda f: floor+((f-center)/scale)**2, crossing_tolerance_hz=.01)
    refined_edges = np.array([refined.bands[0].lower.frequency_hz, refined.bands[0].upper.frequency_hz])
    refined_error = float(max(abs(refined_edges-exact_edges)))
    frequency = np.linspace(.1e6, 4e6, 101)
    flat = analyze_squeezing_bands(_trace(frequency, np.full(frequency.shape, .4)))
    ratio = floor+((frequency-center)/scale)**2
    numerator = ratio*2e-24
    sql = np.full(len(frequency), 2e-24)
    response_power = 1/(1+(frequency/.6e6)**2)
    filtered = SpectralTrace.from_psds(AnalysisFrequencyAxis.from_hz(frequency),
        numerator*response_power, sql*response_power, layer="detected",
        provenance="same electronic low-pass applied to current PSD and SQL")
    filtered_analysis = analyze_squeezing_bands(filtered)
    original = analyze_squeezing_bands(_trace(frequency, ratio))
    lowpass_width_error = abs(filtered_analysis.bands[0].width_hz-original.bands[0].width_hz)
    extreme = analyze_squeezing_bands(_trace([1e5, 1e6, 1e7], [1e308, .5, 1e308]))
    extreme_valid = bool(np.isfinite(extreme.bands[0].width_hz) and extreme.bands[0].width_hz >= 0)
    return {"formula": "R(f)=0.25+((f-center)/scale)^2; linear-ratio threshold interpolation",
            "center_hz": center, "scale_hz": scale, "exact_edges_hz": exact_edges.tolist(),
            "refinements": cases, "locally_refined": _band_dict(refined),
            "locally_refined_maximum_edge_error_hz": refined_error,
            "flat_sub_sql": _band_dict(flat),
            "common_electronic_lowpass_width_change_hz": float(lowpass_width_error),
            "extreme_ratio_crossings_finite_and_ordered": extreme_valid,
            "passed": bool(cases[-1]["maximum_edge_error_hz"] < 3
                and cases[-1]["maximum_edge_error_hz"] < cases[0]["maximum_edge_error_hz"]/50
                and flat.bands[0].width_hz is None and lowpass_width_error < 1e-7
                and refined_error < .01 and extreme_valid)}


def singular_input_controls():
    parameters = tuple(_evidence(f"shared_{i}", 0., "1", 1.) for i in range(5))
    model = JointInputUncertainty(parameters, np.ones((5, 5)))
    samples = sample_inputs(model, 100, seed=0)
    # For exact common errors x0=x1, even a large response cancels identically.
    jacobian = np.array([[1e9, -1e9, 0., 0., 0.]])
    result = propagate_first_order(model, lambda x: jacobian@x,
        jacobian=jacobian, output_ids=("common_error_contrast",), output_units=("1",))
    difference = float(np.max(np.abs(samples.values[:, 0]-samples.values[:, 1])))
    variance = float(result.covariance[0, 0])
    return {"input_count": 5, "correlation": "all entries exactly one",
            "sample_difference_maximum": difference, "amplified_contrast_variance": variance,
            "minimum_correlation_eigenvalue": model.minimum_correlation_eigenvalue,
            "factorization_correlation_max_error": model.factorization_correlation_max_error,
            "passed": difference == 0. and variance == 0.}


def correlated_band_controls():
    model = JointInputUncertainty((
        _evidence("center", 2.1e6, "Hz", 2e4),
        _evidence("scale", 1.3e6, "Hz", 5e4)), [[1., .6], [.6, 1.]])
    root = np.sqrt(.75)
    # Independent closed-form derivative; it never calls the band implementation.
    exact_jacobian = np.array([[1., -root], [1., root], [0., 2*root]])
    exact_covariance = exact_jacobian@model.covariance@exact_jacobian.T
    frequency = np.linspace(.1e6, 4e6, 81)

    def measured_edges(values):
        center, scale = values
        if scale <= 0:
            raise ValueError("positive analytic scale required")
        analysis = analyze_squeezing_bands(_trace(frequency, .25+((frequency-center)/scale)**2),
            crossing_evaluator=lambda f: .25+((f-center)/scale)**2, crossing_tolerance_hz=.0001)
        if len(analysis.bands) != 1 or analysis.bands[0].width_hz is None:
            raise ValueError("finite single band required for this particular edge-uncertainty fixture")
        band = analysis.bands[0]
        return [band.lower.frequency_hz, band.upper.frequency_hz, band.width_hz]

    outputs = dict(output_ids=("lower_edge", "upper_edge", "bandwidth"), output_units=("Hz",)*3)
    linear = propagate_first_order(model, measured_edges, steps=[2e3, 5e3], **outputs)
    ensemble = propagate_ensemble(model, measured_edges, count=512, seed=9112026, **outputs)
    replay = propagate_ensemble(model, measured_edges, count=512, seed=9112026, **outputs)
    relative = np.linalg.norm(linear.covariance-exact_covariance)/np.linalg.norm(exact_covariance)
    sample_covariance_error = np.linalg.norm(ensemble.covariance-exact_covariance)/np.linalg.norm(exact_covariance)
    audit = model.audit_evidence()
    return {"input_values": model.values.tolist(), "input_units": list(model.units),
            "input_standard_uncertainties": model.standard_uncertainties.tolist(),
            "input_correlation": model.correlation.tolist(), "evidence_audit": asdict(audit),
            "output_ids": list(linear.output_ids), "output_units": list(linear.output_units),
            "nominal": linear.nominal.tolist(), "exact_covariance": exact_covariance.tolist(),
            "linear_covariance": linear.covariance.tolist(),
            "linear_standard_uncertainties": linear.standard_uncertainties.tolist(),
            "linear_covariance_relative_error": float(relative),
            "ensemble_count": 512, "ensemble_seed": 9112026,
            "ensemble_covariance": ensemble.covariance.tolist(),
            "ensemble_covariance_relative_error": float(sample_covariance_error),
            "ensemble_quantiles_025_50_975": np.quantile(ensemble.outputs, [.025, .5, .975], axis=0).tolist(),
            "deterministic_replay": bool(np.array_equal(ensemble.outputs, replay.outputs)),
            "passed": bool(relative < .003 and sample_covariance_error < .2
                and np.array_equal(ensemble.outputs, replay.outputs) and not audit.passed)}


def conditional_atomic_controls():
    inputs = ReducedPowerInputs(.6, .00053, 1e18, 1.2e-7, .0125, 8e-6,
        2*np.pi*.9e9, -2*np.pi*8e6, 2*np.pi*1e5)
    positive = AnalysisFrequencyAxis.from_hz([1e5, 1e6, 4e6])
    signed = AnalysisFrequencyAxis(np.r_[-positive.omega_rad_s[::-1], 0., positive.omega_rad_s])
    detector = DetectorResponse(positive, [.85, .85], np.ones((3, 2)), 1., np.zeros(3),
        "assumed 85 percent flat balanced detection; fixed for uncertainty diagnostic")
    model = JointInputUncertainty((
        _evidence("pump_power_W", inputs.pump_power_W, "W", .006),
        _evidence("pump_waist_m", inputs.pump_waist_m, "m", 5.3e-6)), [[1., .6], [.6, 1.]])

    def prediction(values):
        changed = replace(inputs, pump_power_W=values[0], pump_waist_m=values[1])
        result = power_normalized_readout(changed, signed, detector)
        return np.r_[result.probe_power_gain, result.conjugate_power_gain, result.spectrum.quantum_ratio]

    output_ids = ("probe_power_gain", "conjugate_power_gain", "R_0.1MHz", "R_1MHz", "R_4MHz")
    settings = dict(output_ids=output_ids, output_units=("1",)*5)
    coarse = propagate_first_order(model, prediction, steps=model.standard_uncertainties*.1, **settings)
    fine = propagate_first_order(model, prediction, steps=model.standard_uncertainties*.05, **settings)
    diagonal = JointInputUncertainty(model.parameters, np.eye(2))
    uncorrelated = propagate_first_order(diagonal, prediction, jacobian=fine.jacobian, **settings)
    relative = np.linalg.norm(fine.covariance-coarse.covariance)/np.linalg.norm(fine.covariance)
    baseline = power_normalized_readout(inputs, signed, detector)
    trace = SpectralTrace.from_intensity_difference(baseline.spectrum, contribution="quantum")
    analysis = analyze_squeezing_bands(trace)
    audit = model.audit_evidence()
    return {"scope": "existing single-velocity, undepleted, weak-field bright reduced atom; conditional on all fixed inputs",
            "fixed_consumed_input_ledger": _jsonable(inputs.consumed_inputs(detector, convention="uniform-zeeman-rms")),
            "varied_input_evidence": [asdict(item) for item in model.parameters],
            "input_correlation": model.correlation.tolist(), "evidence_audit": asdict(audit),
            "rf_hz": positive.frequency_hz.tolist(), "output_ids": list(output_ids),
            "nominal": fine.nominal.tolist(), "input_jacobian": fine.jacobian.tolist(),
            "output_covariance": fine.covariance.tolist(),
            "output_standard_uncertainties": fine.standard_uncertainties.tolist(),
            "uncorrelated_standard_uncertainties": uncorrelated.standard_uncertainties.tolist(),
            "step_refinement_covariance_relative_change": float(relative),
            "sampled_band_analysis": _band_dict(analysis),
            "passed": bool(relative < 1e-4 and not audit.passed
                and len(analysis.bands) == 1 and analysis.bands[0].width_hz is None)}


def build_report():
    # Record hashes before and after the run, since another task shares this tree.
    # Include the read-only upstream dependency closure, not only the new modules.
    sources = ["gabes/atoms.py", "gabes/core.py", "gabes/constants.py", "gabes/hyperfine.py",
        "gabes/species.py", "gabes/schemes/fwm.py", "gabes/lineshape.py",
        "gabes/observables.py", "gabes/quantum/contracts.py", "gabes/quantum/reservoirs.py",
        "gabes/quantum/diffusion.py", "gabes/quantum/channels.py", "gabes/quantum/traveling.py",
        "gabes/quantum/sidebands.py", "gabes/quantum/readout.py", "gabes/quantum/uncertainty.py",
        "gabes/quantum/spectrum_analysis.py", "gabes/fwm_quantum/model.py",
        "gabes/fwm_quantum/normalization.py", "gabes/fwm_quantum/inputs.py",
        "gabes/fwm_quantum/field.py", "gabes/fwm_quantum/readout.py",
        "analysis/grand_challenge/downstream_audit.py", "gabes/plot_style.py", "requirements-quantum.txt"]
    hashes = {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in sources}
    bands, uncertainty, atomic = analytic_band_controls(), correlated_band_controls(), conditional_atomic_controls()
    singular = singular_input_controls()
    stable = all(hashlib.sha256((ROOT/p).read_bytes()).hexdigest() == h for p, h in hashes.items())
    return {"schema_version": 1, "stage": "parallel RF-band and correlated-input uncertainty tools",
        "expected_controls_passed": bool(bands["passed"] and uncertainty["passed"]
            and atomic["passed"] and singular["passed"] and stable),
        "absolute_hot_vapor_prediction": False, "experimental_validation": False,
        "sources_unchanged_during_run": stable, "source_sha256": hashes,
        "environment": {"python": platform.python_version(), "numpy": np.__version__, "scipy": scipy.__version__},
        "analytic_bands": bands, "correlated_band_uncertainty": uncertainty,
        "singular_input_controls": singular,
        "conditional_atomic_uncertainty": atomic,
        "references": ["https://www.bipm.org/en/doi/10.59161/jcgm102-2011"],
        "limits": [
            "Band interpolation describes the supplied RF samples; it cannot exclude unresolved narrow features.",
            "Scan-limited intervals do not establish intrinsic finite squeezing bandwidth.",
            "Gaussian input PDFs and correlations are declared assumptions here, not measured apparatus distributions.",
            "First-order covariance is a local approximation; step refinement measures derivative numerics, not model validity.",
            "The atomic sensitivity varies only two pump inputs, conditioned on every other ledger value.",
            "No target data, fitted squeezing efficiency, RF filter fit, or out-of-sample experimental test is used.",
            "Current optical spatial-mode work, finite-seed thermal propagation and full-atom physics are separate tasks.",
        ]}


def save_plot(report, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from gabes.plot_style import apply_gabes_plot_style

    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.2), layout="constrained")
    frequency = np.linspace(.1e6, 4e6, 1281)
    ratio = .25+((frequency-2.1e6)/1.3e6)**2
    axes[0].plot(frequency/1e6, 10*np.log10(ratio), label="Analytic fixture")
    axes[0].axhline(0, color="grey", ls="--", label="SQL")
    for edge in report["analytic_bands"]["exact_edges_hz"]:
        axes[0].axvline(edge/1e6, color="grey", ls=":")
    axes[0].set(xlabel="RF frequency [MHz]", ylabel="Noise / SQL [dB]", title="Independent band-edge control")
    axes[0].legend()
    atomic = report["conditional_atomic_uncertainty"]
    values, sigma = np.array(atomic["nominal"])[2:], np.array(atomic["output_standard_uncertainties"])[2:]
    axes[1].errorbar(np.array(atomic["rf_hz"])/1e6, values, yerr=sigma, fmt="o", capsize=5,
                     label="Declared pump input uncertainty (1 sigma)")
    axes[1].set(xlabel="RF frequency [MHz]", ylabel="Linear noise / SQL",
                title="Existing reduced atom: conditional sensitivity")
    axes[1].legend(fontsize=8)
    fig.suptitle("Grand Challenge parallel tools | No experimental validation")
    apply_gabes_plot_style(fig)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("xb") as stream:
        fig.savefig(stream, format="png", dpi=160)
    plt.close(fig)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--plot", type=Path)
    args = parser.parse_args(argv)
    if args.output is not None and args.plot is not None and args.output.resolve() == args.plot.resolve():
        raise ValueError("report and figure must have different paths")
    for path in (args.output, args.plot):
        if path is not None and path.exists():
            raise FileExistsError(path)
    report = build_report()
    body = json.dumps(report, ensure_ascii=False, indent=2, allow_nan=False)+"\n"
    if args.output is None:
        print(body, end="")
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        with args.output.open("x", encoding="utf-8") as stream:
            stream.write(body)
        print(f"Wrote {args.output}")
    if args.plot is not None:
        save_plot(report, args.plot)
        print(f"Wrote {args.plot}")
    return 0 if report["expected_controls_passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
