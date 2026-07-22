"""
A/B benchmark: FWM compute_spectrum with the numba fused kernel vs the NumPy
δ-loop reference path.

    python analysis/squeezing/bench_numba_fwm.py

Three workloads:
  scan-point : the squeezing-scan unit of work (coarse 81, dv=5, ±0.55 GHz)
  app-fast   : the app's "Fast" resolution   (coarse 121, dv=5, full window)
  app-fine   : the app's "Fine" resolution   (coarse 301, dv=2, full window)
"""
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from gabes import kernels  # noqa: E402
from gabes.core import blas_single_thread  # noqa: E402
from gabes.schemes import fwm  # noqa: E402

CENTER = fwm.branch_center_GHz(0.9, -1)
WORKLOADS = {
    "scan-point": dict(coarse_points=81, fine_points=0, velocity_step=5.0,
                       velocity_cutoff=3.0, scan_min=CENTER - 0.55,
                       scan_max=CENTER + 0.55),
    "app-fast": dict(coarse_points=121, fine_points=0, velocity_step=5.0),
    "app-fine": dict(coarse_points=301, fine_points=0, velocity_step=2.0),
}
COMMON = dict(T=394.15, P_pump=0.6, P_probe=8e-6, line_strength=0.74,
              loss_frac=0.055, branch=-1)


def run(name, reps):
    cfg = {**COMMON, **WORKLOADS[name]}
    best = np.inf
    for _ in range(reps):
        t0 = time.perf_counter()
        spec = fwm.compute_spectrum(0.9, **cfg)
        best = min(best, time.perf_counter() - t0)
    return best, spec


def main():
    if not kernels.available():
        print("numba unavailable -- nothing to compare")
        return

    print("warming up JIT...")
    t0 = time.perf_counter()
    run("scan-point", 1)
    print(f"  first call (incl. compile): {time.perf_counter() - t0:.1f}s\n")

    print(f"{'workload':<12} {'numpy ref':>10} {'numpy 1-BLAS':>13} "
          f"{'numba':>9} {'speedup':>8}")
    for name in WORKLOADS:
        reps = 3 if name != "app-fine" else 2
        t_fast, s_fast = run(name, reps)

        orig = fwm.kernels.available
        fwm.kernels.available = lambda: False
        try:
            t_ref, s_ref = run(name, reps)
            with blas_single_thread():
                t_ref1, _ = run(name, reps)
        finally:
            fwm.kernels.available = orig

        worst = np.abs(s_fast["S_dB"] - s_ref["S_dB"]).max()
        print(f"{name:<12} {t_ref:>9.2f}s {t_ref1:>12.2f}s {t_fast:>8.2f}s "
              f"{t_ref / t_fast:>7.1f}x   |dS|max={worst:.2e} dB")


if __name__ == "__main__":
    main()
