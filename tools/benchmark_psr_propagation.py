"""Reproduce PSR S4 longitudinal pump, noise and thermal-propagation checks.

Example:
    python tools/benchmark_psr_propagation.py --output docs/psr/s4_benchmark.json
The default includes a real 70 C Maxwell calculation. Use
--skip-thermal only for a faster single-velocity development run.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timedelta, timezone
import hashlib
import json
from pathlib import Path
import platform
import subprocess
import sys

import numpy as np
import scipy

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from gabes.psr_conventions import jones_from_angles, jones_polarization_angles
from gabes.psr_full_obe import solve_psr_pump
from gabes.psr_homodyne import psr_homodyne_readout, psr_output_homodyne_readout
from gabes.psr_propagation import solve_psr_propagation
from gabes.psr_spectrum import psr_local_generator, psr_uniform_transfer
from tools.benchmark_psr_thermal import _difference, _measure, _serialize


_SOURCES = (
    "gabes/psr_propagation.py", "gabes/psr_thermal.py", "gabes/psr_full_obe.py",
    "gabes/psr_spectrum.py", "gabes/psr_homodyne.py", "gabes/psr_conventions.py",
    "gabes/magneto_multilevel.py", "gabes/species.py", "gabes/beam.py", "gabes/constants.py",
    "gabes/core.py", "gabes/quantum/traveling.py", "gabes/quantum/contracts.py",
    "gabes/quantum/channels.py", "tools/benchmark_psr_thermal.py",
    "tools/benchmark_psr_propagation.py")


def _hashes():
    return {name: hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in _SOURCES}


def _result(path, efficiency=.9025):
    result = {key: path[key] for key in (
        "model_id", "conditions", "profile", "segments", "output_intensity_mw_cm2",
        "output_jones_xy", "carrier_phase_advance_rad", "convergence", "performance")}
    result["velocity_convergence"] = path.get("velocity_convergence")
    result["output_polarization_angles_rad"] = jones_polarization_angles(path["output_jones_xy"])
    if path["transfer"] is not None:
        result["transfer_audit"] = path["transfer"].audit()
        readout = psr_output_homodyne_readout(path["transfer"], path["output_jones_xy"],
                                             detection_efficiency=efficiency)
        result["readout"] = {key: readout[key] for key in (
            "frequencies_hz", "minimum_db", "maximum_db", "min_angle_deg", "covariance",
            "minimum_variance", "maximum_variance", "detection_efficiency")}
    return result


def _compare_paths(reference, candidate):
    result = {key: _difference(reference[key], candidate[key]) for key in
              ("output_intensity_mw_cm2", "output_jones_xy", "carrier_phase_advance_rad")}
    a, b = reference["transfer"], candidate["transfer"]
    if a is not None and b is not None:
        result["field_transfer"] = {key: _difference(getattr(a, key), getattr(b, key)) for key in
                                    ("transfer", "noise_greater", "noise_lesser")}
        ra, rb = (_result(item)["readout"] for item in (reference, candidate))
        result["readout"] = {key: _difference(ra[key], rb[key]) for key in
                             ("minimum_db", "maximum_db", "covariance", "min_angle_deg")}
    return result


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--skip-thermal", action="store_true")
    parser.add_argument("--max-segments", type=int, default=32)
    parser.add_argument("--reference-segments", type=int, default=128)
    args = parser.parse_args(argv)
    if args.max_segments < 4 or args.reference_segments < 2*args.max_segments:
        parser.error("max-segments >=4 and reference-segments >=2*max-segments required")
    source_before = _hashes()
    rf = np.array([-4e5, 4e5])
    detuning, length = 400., .025
    inputs = dict(intensity_mw_cm2=12.5, temp_c=70., density_m3=1e17,
                  ground_relaxation_hz=1e5, jones_xy=jones_from_angles(.31, .17))
    # Prepare immutable atomic tables before collecting allocation/timing data.
    warm = solve_psr_pump(detuning, **inputs)
    psr_local_generator(warm, rf, method="optimized")
    adaptive, timing = _measure(lambda: solve_psr_propagation(
        detuning, rf, length, **inputs, segments=2, max_segments=args.max_segments))
    reference, reference_timing = _measure(lambda: solve_psr_propagation(
        detuning, rf, length, **inputs, segments=args.reference_segments//2,
        max_segments=args.reference_segments))
    scaled, scaled_timing = _measure(lambda: solve_psr_propagation(
        detuning, rf, length/2, **dict(inputs, density_m3=2*inputs["density_m3"]),
        segments=2, max_segments=args.max_segments))

    linear = dict(inputs, jones_xy=np.array([1., 0.]))
    fixed, fixed_timing = _measure(lambda: solve_psr_propagation(
        detuning, rf, length, **linear, freeze_pump=True))
    state = solve_psr_pump(detuning, **linear)
    direct = psr_uniform_transfer(psr_local_generator(state, rf, method="direct"), state, length)
    frozen_readout = psr_output_homodyne_readout(fixed["transfer"], fixed["output_jones_xy"], detection_efficiency=.9025)
    direct_readout = psr_homodyne_readout(direct, detection_efficiency=.9025)
    uniform = dict(result=_result(fixed), measurement=fixed_timing,
        comparison={key: _difference(getattr(direct, key), getattr(fixed["transfer"], key)) for key in
                    ("transfer", "noise_greater", "noise_lesser")},
        readout_comparison={key: _difference(direct_readout[key], frozen_readout[key]) for key in
                            ("covariance", "minimum_db", "maximum_db")})

    thermal = None
    if not args.skip_thermal:
        thermal_inputs = dict(intensity_mw_cm2=.1, temp_c=70., density_m3=1e15,
                              ground_relaxation_hz=1000., jones_xy=jones_from_angles(.31, .17),
                              velocity_model="thermal", segments=1, max_segments=2,
                              thermal_options=dict(velocity_order=2, max_order=6, rtol=.02),
                              max_atomic_evaluations=10000)
        thermal_path, thermal_timing = _measure(lambda: solve_psr_propagation(
            -900., rf, .001, **thermal_inputs))
        thermal = dict(result=_result(thermal_path), measurement=thermal_timing,
                       scope="Real 70 C Maxwell classes, short uniform-density cell; full pump update and per-segment noise.",
                       requested=dict(detuning_mhz=-900., length_m=.001, rf_hz=rf, **thermal_inputs))

    try:
        head = subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT,
                              capture_output=True, text=True, check=True).stdout.strip()
    except (subprocess.SubprocessError, OSError):
        head = None
    after = _hashes()
    changed = {name: dict(before=source_before[name], after=after[name])
               for name in _SOURCES if source_before[name] != after[name]}
    report = dict(
        provenance="\u27e6AI\u27e7 Numerical S4 implementation benchmark; not experimental validation",
        created_at=datetime.now(timezone(timedelta(hours=9))).isoformat(),
        environment=dict(python=platform.python_version(), numpy=np.__version__, scipy=scipy.__version__,
                         platform=platform.platform(), git_head=head, source_sha256=source_before,
                         source_state="working tree; HEAD alone does not identify uncommitted sources",
                         changed_during_measurement=changed),
        command=dict(max_segments=args.max_segments, reference_segments=args.reference_segments,
                     skip_thermal=args.skip_thermal),
        memory_measure="tracemalloc peak allocations during each measured call, excluding pre-existing inputs; not process RSS",
        requested=dict(detuning_mhz=detuning, rf_hz=rf, length_m=length, detection_efficiency=.9025, **inputs),
        depletion=dict(result=_result(adaptive), measurement=timing),
        finer_grid=dict(result=_result(reference), measurement=reference_timing,
                        comparison=_compare_paths(reference, adaptive)),
        density_length_scaling=dict(result=_result(scaled), measurement=scaled_timing,
                                    comparison=_compare_paths(adaptive, scaled),
                                    scope="At fixed temperature, reset, incident field and constant density, doubling n and halving L preserves the equations in n*z. Geometry-dependent transit/collision changes are not included."),
        uniform_reference=uniform,
        room_temperature_thermal=thermal)
    output = args.output if args.output.is_absolute() else ROOT/args.output
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(_serialize(report), ensure_ascii=False, indent=2, allow_nan=False)+"\n", encoding="utf-8")
    passed = bool(adaptive["convergence"]["passed"] and reference["convergence"]["passed"]
                  and scaled["convergence"]["passed"] and fixed["convergence"]["passed"] and not changed
                  and (thermal is None or thermal["result"]["convergence"]["passed"]))
    print(json.dumps(dict(output=str(output), converged=passed,
                          depletion_seconds=timing["elapsed_seconds"],
                          thermal_seconds=None if thermal is None else thermal["measurement"]["elapsed_seconds"])))
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
