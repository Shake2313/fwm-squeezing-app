"""Frozen complex-reference comparison for the magneto speed patch.

Run the complete campaign only after implementation settles::

    python -m analysis.magneto_speed_20261030.validation --output validation.json

The four main cases use one warm-up followed by three interleaved repetitions.
Model/constants/JIT caches remain warm. Direct scheme calls bypass Streamlit's
outer compute cache; no additional atomic-result LRU was introduced by this patch.
Headless-readout latency is reported separately.
No historical evidence is overwritten. All reports stay in this new directory.
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


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
BEFORE_BACKEND = HERE / "before_magneto_multilevel.py"
BEFORE_UNIFIED = HERE / "before_magneto_unified.py"
RTOL = 2e-8
BACKEND_ATOL = {
    "chi_plus": 2e-13, "chi_minus": 2e-13, "chi_probe": 2e-13,
    "transmission": 2e-11, "optical_depth": 2e-10, "rotation_rad": 2e-11,
    "chi_plus_slope_per_t": 2e-4, "chi_minus_slope_per_t": 2e-4,
    "chi_probe_slope_per_t": 2e-4, "rotation_slope_rad_per_t": .02,
    "transmission_slope_per_t": .02,
}
_BASELINES = None


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module  # dataclasses resolve their defining module here.
    try:
        spec.loader.exec_module(module)
    except BaseException:
        sys.modules.pop(name, None)
        raise
    return module


def load_baselines():
    """Load frozen modules with relative imports, isolating their backend binding."""
    global _BASELINES
    if _BASELINES is None:
        import gabes
        from gabes import magneto_multilevel as current
        # Finish live scheme/package registration before temporarily changing
        # the package attribute: otherwise its first import could bind the
        # frozen backend into the live UnifiedMagnetoScheme as a side effect.
        from gabes.schemes import magneto_unified as live_unified
        backend = _load("gabes._magneto_speed_before", BEFORE_BACKEND)
        original = gabes.magneto_multilevel
        try:
            # `from .. import magneto_multilevel` binds the attribute during load.
            gabes.magneto_multilevel = backend
            unified = _load("gabes.schemes._magneto_speed_before", BEFORE_UNIFIED)
        finally:
            gabes.magneto_multilevel = original
        assert unified.magneto_multilevel is backend
        assert gabes.magneto_multilevel is current
        assert live_unified.magneto_multilevel is current
        _BASELINES = backend, unified
    return _BASELINES


def array_check(before, after, *, atol, rtol=RTOL):
    """Mixed absolute/relative gate; zero-valued quantities get an explicit floor."""
    a, b = np.asarray(before), np.asarray(after)
    if a.dtype.kind == "O":
        a = np.asarray(before, float)
    if b.dtype.kind == "O":
        b = np.asarray(after, float)
    if a.shape != b.shape:
        return dict(passed=False, before_shape=list(a.shape), after_shape=list(b.shape))
    finite_a, finite_b = np.isfinite(a), np.isfinite(b)
    mask_mismatches = int(np.count_nonzero(finite_a != finite_b))
    same_nonfinite = np.array_equal(np.isnan(a), np.isnan(b)) and np.array_equal(np.isinf(a), np.isinf(b))
    use = finite_a & finite_b
    error = np.abs(a[use] - b[use])
    scale = np.maximum(np.abs(a[use]), np.abs(b[use]))
    allowed = atol + rtol * scale
    ratio = error / np.maximum(allowed, np.finfo(float).tiny)
    return dict(passed=bool(same_nonfinite and not mask_mismatches and np.all(ratio <= 1)),
                max_absolute_error=float(np.max(error, initial=0.)),
                max_tolerance_ratio=float(np.max(ratio, initial=0.)),
                finite_mask_mismatches=mask_mismatches, finite_samples=int(np.count_nonzero(use)),
                absolute_tolerance=atol, relative_tolerance=rtol)


def exact_check(before, after):
    a, b = np.asarray(before), np.asarray(after)
    same_shape = a.shape == b.shape
    return dict(passed=bool(same_shape and np.array_equal(a, b)),
                mismatches=int(np.count_nonzero(a != b)) if same_shape else None)


def compare_backend(before, after):
    checks = {key: array_check(before[key], after[key], atol=atol)
              for key, atol in BACKEND_ATOL.items()}
    for key in ("slope_valid", "transmission_slope_symmetry_zero", "b_ut", "detunings_mhz"):
        checks[key] = exact_check(before[key], after[key])
    for key in ("solver_status", "magnetic_response_status", "optical_excitation_status"):
        checks["diagnostics." + key] = exact_check(before["diagnostics"][key], after["diagnostics"][key])
    return dict(passed=all(value["passed"] for value in checks.values()), checks=checks,
                before_diagnostics=before["diagnostics"], after_diagnostics=after["diagnostics"],
                before_velocity_quadrature=before["velocity_quadrature"],
                after_velocity_quadrature=after["velocity_quadrature"])


def _conditioned_asd_check(a, b, key, before_readout, after_readout):
    """Gate reciprocal ASD only where both slopes are below the declared floor.

    Raw ASD errors and every affected sample remain explicit. Availability and
    all finite SNR values are checked independently, without this conditioning.
    """
    rotation = "slope_rad_per_t" in a
    slope_key = "rotation_slope_rad_per_t" if rotation else "transmission_slope_per_t"
    noise_key = "rotation_asd_rad_per_sqrt_hz" if rotation else "transmission_asd_per_sqrt_hz"
    sa, sb = (np.asarray(readout[slope_key], float) for readout in (before_readout, after_readout))
    va, vb = np.asarray(a[key], float), np.asarray(b[key], float)
    floor = BACKEND_ATOL[slope_key]
    near = (np.isfinite(va) & np.isfinite(vb) & (abs(sa) <= floor) & (abs(sb) <= floor))
    raw = array_check(va, vb, atol=1e-21, rtol=2e-7)
    resolved = array_check(va[~near], vb[~near], atol=1e-21, rtol=2e-7)
    noise = np.minimum(np.asarray(a[noise_key], float), np.asarray(b[noise_key], float))[near]
    if key == "field_rms_t":
        noise *= np.sqrt(float(a["bandwidth_hz"]))
    inv_a, inv_b = 1. / va[near], 1. / vb[near]
    tolerance = floor / noise + 2e-7 * np.maximum(abs(inv_a), abs(inv_b))
    ratio = abs(inv_a - inv_b) / tolerance
    inverse = dict(passed=bool(np.all(np.isfinite(ratio)) and np.all(ratio <= 1.)),
                   max_tolerance_ratio=float(np.max(ratio, initial=0.)),
                   absolute_tolerance_rule="declared absolute slope tolerance / signal noise ASD (or RMS)")
    affected = []
    optical = before_readout["x_unit"] == "MHz"
    for flat_index, ia, ib, tol in zip(np.flatnonzero(near), inv_a, inv_b, tolerance):
        i = int(flat_index)
        affected.append(dict(index=i, before_slope=float(sa[i]), after_slope=float(sb[i]),
            physical_b_ut=float(before_readout["physical_bias_ut"] if optical else before_readout["physical_b_ut"][i]),
            optical_detuning_mhz=float(before_readout["x"][i] if optical else before_readout["detunings_mhz"][0]),
            before_asd_or_rms=float(va[i]), after_asd_or_rms=float(vb[i]),
            absolute_asd_or_rms_difference=float(abs(va[i] - vb[i])),
            before_inverse=float(ia), after_inverse=float(ib), inverse_absolute_tolerance=float(tol)))
    return dict(passed=raw["finite_mask_mismatches"] == 0 and resolved["passed"] and inverse["passed"],
                raw_asd_check=raw, resolved_asd_check=resolved, near_zero_inverse_asd_check=inverse,
                slope_key=slope_key, absolute_slope_floor=floor, affected_samples=affected)


def compare_readout(before_raw, after_raw, before_readout, after_readout):
    checks = {}
    for key in ("transmission", "rotation_rad", "rotation_slope_rad_per_t", "transmission_slope_per_t"):
        checks[key] = array_check(before_readout[key], after_readout[key], atol=BACKEND_ATOL[key])
    for key in ("slope_valid", "operating_index", "optimum_index", "optimum_detuning_mhz", "optimum_at_scan_edge"):
        checks[key] = exact_check(before_readout[key], after_readout[key])
    for budget_name in ("sensitivity", "rotation_sensitivity", "transmission_sensitivity"):
        a, b = before_readout[budget_name], after_readout[budget_name]
        for key in ("valid", "snr_valid", "signal_linearity_valid"):
            checks[budget_name + "." + key] = exact_check(a[key], b[key])
        for key in ("field_asd_t_per_sqrt_hz", "field_rms_t"):
            checks[budget_name + "." + key] = _conditioned_asd_check(a, b, key, before_readout, after_readout)
        checks[budget_name + ".signal_to_noise_ratio"] = array_check(
            a["signal_to_noise_ratio"], b["signal_to_noise_ratio"], atol=1e-12, rtol=2e-7)
    # Keep the historical 5% finite-amplitude rule. Its raw relative error is
    # roundoff-dominated when the predicted signal is zero; report those values
    # without converting a giant relative-on-zero difference into a curve error.
    linearity = {}
    for name, derivative in (("rotation", "rotation_slope_rad_per_t"),
                             ("transmission", "transmission_slope_per_t")):
        a, b = before_raw["linearity"][name], after_raw["linearity"][name]
        checks["linearity." + name + ".valid"] = exact_check(a["valid"], b["valid"])
        remainders = []
        for raw in (before_raw, after_raw):
            center = raw["operating_field_index"]
            fields = np.asarray(raw["physical_b_ut"])
            left = int(np.argmin(abs(fields - (fields[center] - 1e-6))))
            right = int(np.argmin(abs(fields - (fields[center] + 1e-6))))
            values = np.asarray(raw["backend"]["rotation_rad" if name == "rotation" else "transmission"])
            prediction = np.asarray(raw["backend"][derivative])[center] * 1e-12
            remainders.append(np.maximum(abs(values[right] - values[center] - prediction),
                                         abs(values[center] - values[left] - prediction)))
        checks["linearity." + name + ".native_unit_remainder"] = array_check(
            *remainders, atol=5e-15, rtol=2e-7)
        linearity[name] = dict(before_relative_error=a["relative_error"],
                               after_relative_error=b["relative_error"],
                               before_native_unit_remainder=remainders[0], after_native_unit_remainder=remainders[1],
                               unchanged_formula="max(abs(right-center-prediction), abs(center-left-prediction))/max(abs(prediction),1e-30)",
                               derivative=derivative, pass_threshold=.05)
    return dict(passed=all(value["passed"] for value in checks.values()), checks=checks,
                finite_amplitude_linearity=linearity,
                before_optimum_detuning_mhz=before_readout["optimum_detuning_mhz"],
                after_optimum_detuning_mhz=after_readout["optimum_detuning_mhz"],
                before_optimum_snr=before_readout["optimum_snr"],
                after_optimum_snr=after_readout["optimum_snr"])


def main_cases():
    _, baseline = load_baselines()
    scheme = baseline.UnifiedMagnetoScheme()
    base = scheme.defaults()
    presets = {preset.name: dict(base, **preset.values) for preset in scheme.presets()}
    return dict(baseB83=base, optical49x3=dict(base, scan_axis=baseline.OPTICAL_SCAN),
                Hanle123=presets["Hanle"], NMOR123=presets["NMOR"])


def branch_cases():
    """Bounded coverage of model branches, avoiding a large Cartesian campaign."""
    cases = {}
    common = dict(intensity_mw_cm2=.07, density_m3=4.5e15, b_ut=[0., .0002],
                  velocity_points=2, wall_relaxation_hz=3.)
    for isotope in ("Rb85", "Rb87"):
        for line in ("D1", "D2"):
            for doppler in (False, True):
                cases[f"{isotope}_{line}_doppler_{doppler}"] = dict(common, isotope=isotope, line=line,
                                                                      doppler=doppler)
    special = dict(common, isotope="Rb87", line="D1", doppler=False)
    variants = dict(gamma_out_zero=dict(gamma_out_hz=0.), gamma_in_zero=dict(gamma_in_hz=0.),
                    independent_velocities=dict(velocity_mixing="none", doppler=True),
                    transverse_xy=dict(b_x_ut=.0007, b_y_ut=-.0003),
                    elliptical=dict(qwp_deg=17.), circular=dict(qwp_deg=45.),
                    nonzero_bias=dict(b_ut=[.0001, .0003]), lower_ground=dict(ground_f=1),
                    derivative_off=dict(magnetic_derivative=False), zero_intensity=dict(intensity_mw_cm2=0.))
    cases.update({name: dict(special, **changes) for name, changes in variants.items()})
    return cases


def small_scheme_cases():
    _, baseline = load_baselines()
    base = dict(baseline.UnifiedMagnetoScheme().defaults(), doppler="off", velocity_points=2,
                scan_points=7, detuning_points=7, detuning_max_mhz=600.)
    return dict(optical_linear=dict(base, scan_axis=baseline.OPTICAL_SCAN),
                optical_elliptical=dict(base, scan_axis=baseline.OPTICAL_SCAN,
                                       signal_type=baseline.TRANSMISSION, qwp_deg=17., b_offset_ut=.0002),
                biased_transverse=dict(base, b_offset_ut=.00017, residual_transverse_b_ut=.0003,
                                       transverse_field_angle_deg=31.))


def evaluate_scheme(module, params):
    scheme = module.UnifiedMagnetoScheme()
    raw = scheme.compute(params)
    start = time.perf_counter()
    readout = scheme.headless_observables(raw, params)["magneto_readout"]
    readout_seconds = time.perf_counter() - start
    return raw, readout, readout_seconds


def compare_scheme_case(params, *, repeats=0):
    from gabes.schemes import magneto_unified as current
    _, baseline = load_baselines()
    results = {}
    for label, module in (("before", baseline), ("after", current)):
        results[label] = evaluate_scheme(module, params)
    old_raw, old_readout, _ = results["before"]
    new_raw, new_readout, _ = results["after"]
    if hasattr(current, "_display_response"):
        # The new cached raw state uses fixed reference density/path. Compare
        # physical curves after the same rescaling used by the public readout.
        physical, conditions, linearity = current._display_response(new_raw, params)
        new_raw = dict(new_raw, backend=physical, conditions=conditions, linearity=linearity)
    backend = compare_backend(old_raw["backend"], new_raw["backend"])
    readout = compare_readout(old_raw, new_raw, old_readout, new_readout)
    report = dict(passed=backend["passed"] and readout["passed"], params=params,
                  backend=backend, readout=readout,
                  actual_field_points=len(old_raw["physical_b_ut"]),
                  actual_optical_points=len(old_raw["detunings_mhz"]))
    if repeats:
        samples = {"before": [], "after": []}
        reads = {"before": [], "after": []}
        modules = {"before": baseline, "after": current}
        # Initial evaluations above warm code/model/JIT paths for each side.
        for repeat in range(repeats):
            for label in (("before", "after") if repeat % 2 == 0 else ("after", "before")):
                module = modules[label]
                scheme = module.UnifiedMagnetoScheme()
                start = time.perf_counter()
                raw = scheme.compute(params)
                samples[label].append(time.perf_counter() - start)
                start = time.perf_counter()
                scheme.headless_observables(raw, params)
                reads[label].append(time.perf_counter() - start)
        medians = {label: statistics.median(values) for label, values in samples.items()}
        report["timing"] = dict(repetitions=repeats, order="before/after, after/before, alternating",
            warmup_per_side=1, direct_scheme_compute=True, outer_ui_cache_bypassed=True,
            compute_samples_seconds=samples, compute_median_seconds=medians,
            solver_speedup=medians["before"] / medians["after"],
            headless_readout_samples_seconds=reads,
            ui_compute_cache="not timed; external Streamlit cache bypassed by direct scheme calls")
    return report


def preserved_hashes():
    paths = [BEFORE_BACKEND, BEFORE_UNIFIED, ROOT / "tests/test_magneto_full_obe_reference.py",
             ROOT / "gabes/magneto_noise.py", ROOT / "gabes/schemes/magneto.py"]
    for directory in ("docs/magneto_budker_20261030", "analysis/magneto_budker_20261030/results"):
        paths.extend((ROOT / directory).glob("*.json"))
    return {str(path.relative_to(ROOT)).replace("\\", "/"): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in paths if path.exists()}


def _json(value):
    if isinstance(value, np.ndarray):
        return _json(value.tolist())
    if isinstance(value, np.generic):
        return _json(value.item())
    if isinstance(value, dict):
        return {str(k): _json(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json(v) for v in value]
    if isinstance(value, float) and not np.isfinite(value):
        return None
    return value


def run(*, benchmarks=True, repeats=3):
    from gabes import magneto_multilevel as current
    from gabes import core
    import scipy
    try:
        from threadpoolctl import threadpool_info
        pools_outside_solve = threadpool_info()
        with core.blas_single_thread():
            pools_during_solve = threadpool_info()
    except ImportError:
        pools_outside_solve = pools_during_solve = "threadpoolctl unavailable"
    baseline, _ = load_baselines()
    preserved = preserved_hashes()
    report = dict(provenance="⟦AI⟧ frozen-before/new implementation numerical and timing comparison",
                  runtime=dict(python=platform.python_version(), numpy=np.__version__, scipy=scipy.__version__,
                      platform=platform.platform(), processor=platform.processor(),
                      blas_pools_outside_solve=pools_outside_solve, blas_pools_during_solve=pools_during_solve),
                  tolerance_contract=dict(relative=RTOL, absolute_by_quantity=BACKEND_ATOL,
                      exact_masks=True, best_sampled_detuning_exact=True,
                      zero_roundoff="absolute floors; ASD inverse compared only where both slopes are below their declared floor; all affected samples/raw ASD errors retained; finite SNR, masks and optimum never waived"),
                  source_sha256={str(path.relative_to(ROOT)).replace("\\", "/"):hashlib.sha256(path.read_bytes()).hexdigest()
                                 for path in (Path(__file__), ROOT / "gabes/magneto_multilevel.py",
                                              ROOT / "gabes/magneto_coordinates.py",
                                              ROOT / "gabes/schemes/magneto_unified.py",
                                              ROOT / "gabes/species.py", ROOT / "gabes/zeeman.py",
                                              ROOT / "gabes/core.py", ROOT / "gabes/constants.py")},
                  preserved_before_sha256=preserved, backend_cases={}, scheme_cases={}, main_cases={})
    for name, kwargs in branch_cases().items():
        before = baseline.multilevel_response([-137.], **kwargs)
        after = current.multilevel_response([-137.], **kwargs)
        report["backend_cases"][name] = dict(kwargs=kwargs, **compare_backend(before, after))
        print(f"backend {name}: {report['backend_cases'][name]['passed']}", flush=True)
    for name, params in small_scheme_cases().items():
        report["scheme_cases"][name] = compare_scheme_case(params)
        print(f"scheme {name}: {report['scheme_cases'][name]['passed']}", flush=True)
    if benchmarks:
        for name, params in main_cases().items():
            report["main_cases"][name] = compare_scheme_case(params, repeats=repeats)
            entry = report["main_cases"][name]
            print(f"main {name}: {entry['passed']}; speedup={entry['timing']['solver_speedup']:.3f}", flush=True)
    report["preserved_after_sha256"] = preserved_hashes()
    report["preservation_passed"] = preserved == report["preserved_after_sha256"]
    report["source_after_sha256"] = {name:hashlib.sha256((ROOT / name).read_bytes()).hexdigest()
                                     for name in report["source_sha256"]}
    report["sources_unchanged_during_run"] = report["source_sha256"] == report["source_after_sha256"]
    report["passed"] = report["preservation_passed"] and report["sources_unchanged_during_run"] and all(
        item["passed"] for section in ("backend_cases", "scheme_cases", "main_cases") for item in report[section].values())
    return _json(report)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", default="validation.json")
    parser.add_argument("--skip-benchmarks", action="store_true")
    parser.add_argument("--repeats", type=int, default=3)
    args = parser.parse_args(argv)
    if args.repeats < 1:
        parser.error("--repeats must be positive")
    output = (HERE / args.output).resolve()
    if not output.is_relative_to(HERE) or output.name.startswith("before_"):
        parser.error("reports must stay in the new evidence directory and cannot overwrite frozen baselines")
    result = run(benchmarks=not args.skip_benchmarks, repeats=args.repeats)
    output.write_text(json.dumps(result, ensure_ascii=False, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(f"PASS={result['passed']}; report={output}", flush=True)
    return 0 if result["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
