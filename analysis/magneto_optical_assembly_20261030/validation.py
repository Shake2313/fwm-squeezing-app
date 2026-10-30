"""Compare reusable optical assembly with the committed priorities 1-3 solver.

The frozen backend/scheme use the same live atomic-data dependencies as the new
implementation. Shared comparison helpers and their previously declared numeric
floors are imported without modifying the earlier speed evidence. Every validity
mask, advertised finite SNR and sampled optimum remains a strict comparison.

Run the full campaign once after implementation settles::

    python -m analysis.magneto_optical_assembly_20261030.validation --output validation.json

Direct scheme calls bypass Streamlit's outer result cache. Initial cold-template
compute times include preparation; three alternating warm-template timings do
not cache atomic results. Standalone assembly also measures preparation and cold
versus warm assembly separately. Historical receipts are read-only.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import platform
import statistics
import sys
import time
from pathlib import Path

import numpy as np

from analysis.magneto_speed_20261030.validation import (
    BACKEND_ATOL, _json, array_check, branch_cases, compare_backend,
    compare_readout, exact_check,
)


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
_REFERENCE = None


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    try:
        spec.loader.exec_module(module)
    except BaseException:
        sys.modules.pop(name, None)
        raise
    return module


def load_reference():
    """Isolate frozen scheme -> frozen backend, finishing live imports first."""
    global _REFERENCE
    if _REFERENCE is None:
        import gabes
        from gabes import magneto_multilevel as live_backend
        from gabes.schemes import magneto_unified as live_scheme
        backend = _load("gabes._magneto_optical_before", HERE / "before_magneto_multilevel.py")
        original = gabes.magneto_multilevel
        try:
            gabes.magneto_multilevel = backend
            scheme = _load("gabes.schemes._magneto_optical_before", HERE / "before_magneto_unified.py")
        finally:
            gabes.magneto_multilevel = original
        assert scheme.magneto_multilevel is backend
        assert live_scheme.magneto_multilevel is live_backend
        assert gabes.magneto_multilevel is live_backend
        _REFERENCE = backend, scheme
    return _REFERENCE


def _cache():
    from gabes import magneto_optical
    return magneto_optical


def _prepare(atom, polarization):
    return _cache().prepare_optical_coefficients(atom, polarization)


def _clear():
    _cache().coefficient_cache_clear()


def compare_complete_backend(before, after):
    """Extend the prior physical-output gates to every numeric health diagnostic."""
    report = compare_backend(before, after)
    extra = {}
    for key, value in before["diagnostics"].items():
        if isinstance(value, (int, float, np.number)):
            extra["diagnostics." + key] = array_check(value, after["diagnostics"][key], atol=2e-11)
    for key in ("density_m3", "frequency_hz", "reference_centroid_mhz"):
        extra[key] = exact_check(before[key], after[key])
    extra["velocity_quadrature"] = exact_check(
        json.dumps(_json(before["velocity_quadrature"]), sort_keys=True),
        json.dumps(_json(after["velocity_quadrature"]), sort_keys=True))
    report["additional_checks"] = extra
    report["passed"] = report["passed"] and all(item["passed"] for item in extra.values())
    return report


def _physical_raw(module, raw, params):
    physical, conditions, linearity = module._display_response(raw, params)
    return dict(raw, backend=physical, conditions=conditions, linearity=linearity)


def _evaluate(module, params):
    scheme = module.UnifiedMagnetoScheme()
    start = time.perf_counter()
    raw = scheme.compute(params)
    compute_seconds = time.perf_counter() - start
    start = time.perf_counter()
    readout = scheme.headless_observables(raw, params)["magneto_readout"]
    readout_seconds = time.perf_counter() - start
    return _physical_raw(module, raw, params), readout, compute_seconds, readout_seconds


def compare_scheme_case(params, *, repeats=0):
    from gabes.schemes import magneto_unified as current
    _, reference = load_reference()
    _clear()
    before = _evaluate(reference, params)
    after = _evaluate(current, params)
    cache_after_cold = _cache().coefficient_cache_info()
    backend = compare_complete_backend(before[0]["backend"], after[0]["backend"])
    readout = compare_readout(before[0], after[0], before[1], after[1])
    report = dict(passed=backend["passed"] and readout["passed"], params=params,
                  backend=backend, readout=readout,
                  actual_field_points=len(before[0]["physical_b_ut"]),
                  actual_optical_points=len(before[0]["detunings_mhz"]),
                  cold_preparation_inclusive_compute_seconds=dict(before=before[2], after=after[2]),
                  coefficient_cache_after_cold_compute=cache_after_cold)
    if repeats:
        samples = dict(before=[], after=[])
        reads = dict(before=[], after=[])
        modules = dict(before=reference, after=current)
        for repeat in range(repeats):
            order = ("before", "after") if repeat % 2 == 0 else ("after", "before")
            for label in order:
                # Every call solves the physical ground dynamics again; only
                # velocity-independent optical coefficients may be reused.
                _, _, elapsed, read_elapsed = _evaluate(modules[label], params)
                samples[label].append(elapsed)
                reads[label].append(read_elapsed)
        medians = {label: statistics.median(values) for label, values in samples.items()}
        report["timing"] = dict(repetitions=repeats, warmup_per_side=1,
            order="before/after, after/before, alternating", direct_scheme_compute=True,
            result_cache_used=False, warm_template_compute_samples_seconds=samples,
            warm_template_compute_median_seconds=medians,
            speedup=medians["before"] / medians["after"],
            headless_readout_samples_seconds=reads,
            coefficient_cache_after_warm_computes=_cache().coefficient_cache_info())
    return report


def main_cases():
    _, reference = load_reference()
    base = reference.UnifiedMagnetoScheme().defaults()
    optical = dict(base, scan_axis=reference.OPTICAL_SCAN)
    return dict(optical49x3=optical,
                refined97x3_order10=dict(optical, detuning_points=97, velocity_points=10),
                Rb87D2_optical25x3=dict(optical, isotope="Rb87", line="D2", detuning_points=25),
                baseB83=base)


def small_scheme_cases():
    _, reference = load_reference()
    base = dict(reference.UnifiedMagnetoScheme().defaults(), doppler="off", velocity_points=2,
                scan_points=7, detuning_points=7, detuning_max_mhz=600.)
    return dict(optical_linear=dict(base, scan_axis=reference.OPTICAL_SCAN),
                optical_elliptical=dict(base, scan_axis=reference.OPTICAL_SCAN,
                    signal_type=reference.TRANSMISSION, qwp_deg=17., b_offset_ut=.0002),
                biased_transverse=dict(base, b_offset_ut=.00017, residual_transverse_b_ut=.0003,
                                      transverse_field_angle_deg=31.))


def assembly_cases():
    cases = {}
    for isotope in ("Rb85", "Rb87"):
        for line in ("D1", "D2"):
            for qwp in (0., 17., 45.):
                name = f"{isotope}_{line}_qwp{qwp:g}"
                cases[name] = dict(isotope=isotope, line=line, qwp_deg=qwp,
                                   intensity_mw_cm2=.7, velocity_count=180)
    return cases


def compare_optical_terms(before, after):
    checks = {}
    a, b = before[0], after[0]
    scale = max(1., float(np.max(abs(a), initial=0.)), float(np.max(abs(b), initial=0.)))
    checks["generator"] = array_check(a, b, atol=scale*1e-13, rtol=5e-12)
    for label, index in (("eliminated", 1), ("resolvents", 2)):
        checks[label + ".block_count"] = exact_check(len(before[index]), len(after[index]))
        for bi, (a, b) in enumerate(zip(before[index], after[index])):
            checks[f"{label}.{bi}"] = array_check(a, b, atol=1e-18, rtol=5e-12)
    checks["excitation"] = array_check(before[3], after[3], atol=1e-15, rtol=5e-12)
    return dict(passed=all(item["passed"] for item in checks.values()), checks=checks)


def compare_assembly_case(settings, *, repeats=0):
    from gabes import core, magneto_multilevel as current
    reference, _ = load_reference()
    atom = current.build_multilevel_atom(settings["isotope"], settings["line"])
    old_atom = reference.build_multilevel_atom(settings["isotope"], settings["line"])
    # Nonuniform fixed nodes include resonances and far tails; no random search.
    rng = np.random.default_rng(20261030)
    nodes = np.sort(rng.uniform(-4500., 4500., settings["velocity_count"]))
    polarization = current._polarization(settings["qwp_deg"])
    intensity = settings["intensity_mw_cm2"]
    _clear()
    with core.blas_single_thread():
        before = reference._optical_terms(old_atom, nodes, intensity, polarization)
        after = current._optical_terms(atom, nodes, intensity, polarization)
    report = dict(settings=settings, **compare_optical_terms(before, after),
                  coefficient_cache_after_initial=_cache().coefficient_cache_info())
    if repeats:
        preparation, cold = [], []
        samples = dict(before=[], after=[])
        with core.blas_single_thread():
            for _ in range(repeats):
                _clear()
                start = time.perf_counter()
                _prepare(atom, polarization)
                preparation.append(time.perf_counter() - start)
                _clear()
                start = time.perf_counter()
                current._optical_terms(atom, nodes, intensity, polarization)
                cold.append(time.perf_counter() - start)
            for repeat in range(repeats):
                for label in (("before", "after") if repeat % 2 == 0 else ("after", "before")):
                    module, selected = (reference, old_atom) if label == "before" else (current, atom)
                    start = time.perf_counter()
                    module._optical_terms(selected, nodes, intensity, polarization)
                    samples[label].append(time.perf_counter() - start)
        medians = {label:statistics.median(values) for label, values in samples.items()}
        report["timing"] = dict(repetitions=repeats, velocity_count=len(nodes),
            preparation_only_samples_seconds=preparation, preparation_only_median_seconds=statistics.median(preparation),
            after_cold_preparation_inclusive_samples_seconds=cold,
            after_cold_preparation_inclusive_median_seconds=statistics.median(cold),
            warm_assembly_samples_seconds=samples, warm_assembly_median_seconds=medians,
            warm_speedup=medians["before"] / medians["after"],
            coefficient_cache_info=_cache().coefficient_cache_info())
    return report


def source_hashes():
    paths = [Path(__file__), HERE / "before_magneto_multilevel.py", HERE / "before_magneto_unified.py",
             ROOT / "analysis/magneto_speed_20261030/validation.py"]
    paths.extend(ROOT / "gabes" / name for name in (
        "magneto_multilevel.py", "magneto_optical.py", "magneto_coordinates.py",
        "schemes/magneto_unified.py", "angular.py", "species.py", "zeeman.py",
        "core.py", "constants.py", "magneto_noise.py"))
    return {str(path.relative_to(ROOT)).replace("\\", "/"):hashlib.sha256(path.read_bytes()).hexdigest()
            for path in paths}


def runtime_info():
    import scipy
    from gabes import core
    try:
        from threadpoolctl import threadpool_info
        outside = threadpool_info()
        with core.blas_single_thread():
            inside = threadpool_info()
    except ImportError:
        outside = inside = "threadpoolctl unavailable"
    return dict(python=platform.python_version(), numpy=np.__version__, scipy=scipy.__version__,
                platform=platform.platform(), processor=platform.processor(),
                blas_pools_outside_solve=outside, blas_pools_during_solve=inside)


def run(*, benchmarks=True, repeats=3, include_b=True):
    from gabes import magneto_multilevel as current
    reference, _ = load_reference()
    initial = source_hashes()
    report = dict(provenance="⟦AI⟧ optical-assembly reuse versus committed priorities 1-3; no changed physics or new verification waivers",
        runtime=runtime_info(), source_sha256=initial,
        numeric_comparison_contract=dict(shared_helpers="analysis/magneto_speed_20261030/validation.py",
            relative=2e-8, absolute_by_quantity=BACKEND_ATOL, exact_masks=True,
            exact_sampled_optimum=True, every_advertised_finite_snr_compared=True,
            new_waivers=False, existing_near_zero_asd_conditioning="unchanged helper; raw differences and affected indices retained"),
        timing_contract=dict(result_cache_used=False, before_after_shared_dependencies=True,
            cold_scheme_sample="initial full compute after coefficient cache clear; includes preparation",
            warm_scheme_samples="alternating before/after; only reusable optical coefficients warm",
            preparation_and_cold_assembly="independent three-repeat measurements, clear before each"),
        backend_cases={}, scheme_cases={}, optical_assembly_cases={}, main_cases={})
    for name, kwargs in branch_cases().items():
        before = reference.multilevel_response([-137.], **kwargs)
        after = current.multilevel_response([-137.], **kwargs)
        report["backend_cases"][name] = dict(kwargs=kwargs, **compare_complete_backend(before, after))
        print(f"backend {name}: {report['backend_cases'][name]['passed']}", flush=True)
    for name, params in small_scheme_cases().items():
        report["scheme_cases"][name] = compare_scheme_case(params)
        print(f"scheme {name}: {report['scheme_cases'][name]['passed']}", flush=True)
    for name, settings in assembly_cases().items():
        report["optical_assembly_cases"][name] = compare_assembly_case(settings, repeats=repeats if benchmarks else 0)
        print(f"assembly {name}: {report['optical_assembly_cases'][name]['passed']}", flush=True)
    if benchmarks:
        for name, params in main_cases().items():
            if name == "baseB83" and not include_b:
                continue
            report["main_cases"][name] = compare_scheme_case(params, repeats=1 if name == "baseB83" else repeats)
            entry = report["main_cases"][name]
            print(f"main {name}: {entry['passed']}; speedup={entry['timing']['speedup']:.3f}", flush=True)
    report["source_after_sha256"] = source_hashes()
    report["sources_unchanged_during_run"] = initial == report["source_after_sha256"]
    report["passed"] = report["sources_unchanged_during_run"] and all(
        item["passed"] for section in ("backend_cases", "scheme_cases", "optical_assembly_cases", "main_cases")
        for item in report[section].values())
    return _json(report)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", default="validation.json")
    parser.add_argument("--skip-benchmarks", action="store_true")
    parser.add_argument("--skip-b-regression", action="store_true")
    parser.add_argument("--repeats", type=int, default=3)
    args = parser.parse_args(argv)
    if args.repeats < 1:
        parser.error("--repeats must be positive")
    output = (HERE / args.output).resolve()
    if not output.is_relative_to(HERE) or output.name.startswith("before_"):
        parser.error("reports must stay in the new evidence directory and cannot overwrite frozen baselines")
    report = run(benchmarks=not args.skip_benchmarks, repeats=args.repeats, include_b=not args.skip_b_regression)
    output.write_text(json.dumps(report, ensure_ascii=False, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(f"PASS={report['passed']}; report={output}", flush=True)
    return 0 if report["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
