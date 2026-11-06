"""Reproduce the PSR S3 thermal kernel, cutoff and homodyne benchmark.

Run from any directory, for example:
    python tools/benchmark_psr_thermal.py --output docs/psr/s3_benchmark.json
The optional --rf-points 25 reproduces the UI-sized signed RF grid.
--reference-order 20 increases the fixed quadrature reference order.
No timing threshold is an acceptance criterion.
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
from time import perf_counter
import tracemalloc

import numpy as np
import scipy

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from gabes import constants
from gabes.psr_full_obe import solve_psr_pump
from gabes.psr_homodyne import psr_homodyne_readout
from gabes.psr_spectrum import psr_local_generator, psr_uniform_transfer
from gabes.psr_thermal import solve_psr_thermal, thermal_convergence


def _measure(function):
    # Each run starts its own tracker, so pre-existing input arrays are excluded.
    tracemalloc.start()
    start = perf_counter()
    try:
        value = function()
        elapsed = perf_counter() - start
        current, peak = tracemalloc.get_traced_memory()
    finally:
        tracemalloc.stop()
    return value, dict(elapsed_seconds=elapsed, peak_traced_allocation_bytes=peak,
                      retained_traced_allocation_bytes=current)


def _difference(reference, candidate):
    reference, candidate = np.asarray(reference), np.asarray(candidate)
    error = candidate-reference
    norm = float(np.linalg.norm(reference))
    return dict(max_absolute=float(np.max(np.abs(error))),
                relative_frobenius=float(np.linalg.norm(error)/norm) if norm else None)


def _kernels(reference, candidate):
    fields = ("drift", "noise_greater_by_reservoir", "noise_lesser_by_reservoir")
    return {key: _difference(getattr(reference, key), getattr(candidate, key)) for key in fields}


def _readout(result, length, efficiency):
    transfer = psr_uniform_transfer(result["local"], result["pump_state"], length)
    readout = psr_homodyne_readout(transfer, detection_efficiency=efficiency)
    return dict(frequencies_hz=readout["frequencies_hz"],
                minimum_db=readout["minimum_db"], maximum_db=readout["maximum_db"],
                minimum_variance=readout["minimum_variance"],
                maximum_variance=readout["maximum_variance"],
                min_angle_deg=readout["min_angle_deg"], covariance=readout["covariance"],
                local_audit=result["local"].audit(), transfer_audit=transfer.audit())


def _serialize(value):
    if isinstance(value, np.ndarray):
        return _serialize(value.tolist())
    if isinstance(value, np.generic):
        return _serialize(value.item())
    if isinstance(value, complex):
        return dict(real=value.real, imag=value.imag)
    if isinstance(value, dict):
        return {str(key): _serialize(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [_serialize(item) for item in value]
    if isinstance(value, float) and not np.isfinite(value):
        return None
    return value


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--rf-points", type=int, default=1,
                        help="Positive RF samples; 1 is the selected point, 25 is the UI grid.")
    parser.add_argument("--reference-order", type=int, default=12)
    parser.add_argument("--reference-method", choices=("direct", "optimized"), default="optimized",
                        help="Atomic kernel implementation for the fixed quadrature reference.")
    parser.add_argument("--max-order", type=int, default=10)
    args = parser.parse_args(argv)
    if args.rf_points < 1 or args.reference_order < 2 or args.max_order < 4:
        parser.error("rf-points >=1, reference-order >=2 and max-order >=4 required")

    center_hz, length, efficiency = 4.e5, .025, .9025
    positive = (np.array([center_hz]) if args.rf_points == 1 else
                center_hz*np.geomspace(.1, 10., args.rf_points))
    frequencies = np.r_[-positive[::-1], positive]
    inputs = dict(intensity_mw_cm2=12.5, temp_c=70.,
                  ground_relaxation_hz=1.e5, density_m3=1.e16)
    detuning = 400.
    pump = solve_psr_pump(detuning, **inputs)
    # Warm shared basis / reservoir caches before measuring the two algorithms.
    psr_local_generator(pump, [-center_hz, center_hz], method="direct")
    direct, direct_timing = _measure(lambda: psr_local_generator(pump, frequencies, method="direct"))
    optimized, optimized_timing = _measure(
        lambda: psr_local_generator(pump, frequencies, method="optimized"))
    adaptive, adaptive_timing = _measure(lambda: solve_psr_thermal(
        detuning, frequencies, **inputs, velocity_order=4, max_order=args.max_order,
        max_classes=10000, method="optimized"))
    reference, reference_timing = _measure(lambda: solve_psr_thermal(
        detuning, frequencies, **inputs, velocity_order=args.reference_order,
        max_order=args.reference_order, cutoff_sigma=7., tail_sigma=8.,
        max_classes=10000, method=args.reference_method, reference=True))
    adaptive_readout = _readout(adaptive, length, efficiency)
    reference_readout = _readout(reference, length, efficiency)
    convergence = thermal_convergence(adaptive, length, detection_efficiency=efficiency)
    scalar_comparison = {key: _difference(reference["pump_state"][key], adaptive["pump_state"][key])
                         for key in ("chi_probe", "alpha_per_m", "excited_fraction")}
    try:
        head = subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT,
                              capture_output=True, text=True, check=True).stdout.strip()
    except (subprocess.SubprocessError, OSError):
        head = None
    report = dict(
        provenance="\u27e6AI\u27e7 Numerical implementation benchmark; not experimental validation",
        created_at=datetime.now(timezone(timedelta(hours=9))).isoformat(),
        environment=dict(python=platform.python_version(), numpy=np.__version__,
                         scipy=scipy.__version__, platform=platform.platform(), git_head=head,
                         source_state="working tree; HEAD alone does not identify uncommitted sources",
                         source_sha256={name: hashlib.sha256((ROOT/name).read_bytes()).hexdigest()
                                        for name in ("gabes/psr_full_obe.py", "gabes/psr_spectrum.py",
                                                     "gabes/psr_thermal.py", "gabes/psr_homodyne.py",
                                                     "gabes/magneto_multilevel.py", "gabes/species.py",
                                                     "gabes/core.py", "tools/benchmark_psr_thermal.py")}),
        command=dict(rf_points=args.rf_points, reference_order=args.reference_order,
                     reference_method=args.reference_method, max_order=args.max_order),
        conditions=dict(detuning_mhz=detuning, **inputs, rf_hz=frequencies,
                        cell_length_m=length, detection_efficiency=efficiency,
                        adaptive_cutoffs_sigma=[6., 7.], reference_cutoff_sigma=7.,
                        rtol=.005, atol_per_m=1.e-6),
        memory_measure="tracemalloc peak allocations during each measured call, excluding pre-existing inputs; not process RSS",
        atomic=dict(direct=direct_timing, optimized=optimized_timing,
                    comparison=_kernels(direct, optimized),
                    direct_audit=direct.audit(), optimized_audit=optimized.audit()),
        thermal=dict(
            adaptive=dict(measurement=adaptive_timing, performance=adaptive["performance"],
                          convergence=convergence, readout=adaptive_readout,
                          optical_depth=adaptive["pump_state"]["alpha_per_m"]*length,
                          history=[{key: entry[key] for key in
                                    ("interval_order", "cutoff_sigma", "class_count",
                                     "elapsed_seconds", "omitted_maxwell_probability")}
                                   for entry in adaptive["history"]]),
            reference=dict(measurement=reference_timing, performance=reference["performance"],
                           convergence=reference["convergence"], readout=reference_readout),
            comparison=dict(kernels=_kernels(reference["local"], adaptive["local"]),
                            mean=scalar_comparison,
                            readout={key: _difference(reference_readout[key], adaptive_readout[key])
                                     for key in ("minimum_db", "maximum_db", "minimum_variance",
                                                 "maximum_variance", "covariance")})))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(_serialize(report), ensure_ascii=False, indent=2,
                                     allow_nan=False)+"\n", encoding="utf-8")
    print(json.dumps(dict(output=str(args.output), convergence=convergence["status"],
                          class_evaluations=adaptive["performance"]["class_evaluations"],
                          adaptive_seconds=adaptive_timing["elapsed_seconds"],
                          reference_seconds=reference_timing["elapsed_seconds"],
                          min_db=adaptive_readout["minimum_db"].tolist(),
                          max_db=adaptive_readout["maximum_db"].tolist()), indent=2))


if __name__ == "__main__":
    main()
