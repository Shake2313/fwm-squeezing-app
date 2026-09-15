"""Capture an immutable local pre-hotfix FWM comparison and warm timings.

Run once before production edits: python -m analysis.fwm_gain_hotfix.capture_before
"""
from __future__ import annotations

import hashlib
import json
import os
import platform
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

from gabes.schemes import fwm

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def json_default(value):
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, complex):
        return {"real": value.real, "imag": value.imag}
    raise TypeError(type(value).__name__)


def capture():
    target = HERE / "before"
    target.mkdir(exist_ok=True)
    if (target / "capture.json").exists():
        raise RuntimeError("Pre-hotfix capture exists; do not overwrite it")
    report = {
        "captured_utc": datetime.now(timezone.utc).isoformat(),
        "environment": {
            "python": sys.version,
            "executable": sys.executable,
            "numpy": np.__version__,
            "platform": platform.platform(),
            "cpu_count": os.cpu_count(),
            "numba_available": fwm.kernels.available(),
        },
        "git_head": subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
        "source_sha256": {
            name: digest(ROOT / name)
            for name in (
                "gabes/schemes/fwm.py", "gabes/core.py", "gabes/kernels.py",
                "gabes/observables.py", "gabes/pole_doppler.py",
                "tests/baseline_focused.npz")
        },
        "cases": {},
    }
    if fwm.kernels.available():
        import numba
        report["environment"].update(
            numba=numba.__version__, numba_threads=numba.get_num_threads())
    scheme = fwm.FWMScheme()
    cases = (
        ("fast_default", fwm.FIDELITY_FAST, {}, 7),
        ("balanced_default", fwm.FIDELITY_BALANCED, {}, 7),
        ("ultra_default", fwm.FIDELITY_ULTRA, {}, 5),
        ("ultra_detuned", fwm.FIDELITY_ULTRA,
         dict(opd=1.5, temp_c=110.0, pump_mw=400.0, probe_uw=10.0), 1),
    )
    for name, fidelity, overrides, repeats in cases:
        params = dict(scheme.defaults(), resolution=fidelity)
        params.update(overrides)
        started = time.perf_counter()
        raw = scheme.compute(params)
        warmup_seconds = time.perf_counter() - started
        print(f"{name}: warmup {warmup_seconds:.3f} s", flush=True)
        samples = []
        for _ in range(repeats):
            started = time.perf_counter()
            raw = scheme.compute(params)
            samples.append(time.perf_counter() - started)
        numeric_arrays = {key: value for key, value in raw.items()
                          if isinstance(value, np.ndarray)
                          and value.dtype.kind in "biufc"}
        artifact = target / f"{name}.npz"
        np.savez_compressed(artifact, **numeric_arrays)
        operating = fwm.operating_point(raw, -8.0)
        report["cases"][name] = {
            "params": params, "operating_point_tpd_mhz": -8.0,
            "operating_point": operating,
            "max_G_s": float(np.max(raw["G_s"])),
            "warmup_seconds": warmup_seconds,
            "seconds": samples, "seconds_median": float(np.median(samples)),
            "numeric_array_keys": sorted(numeric_arrays),
            "artifact": artifact.name, "artifact_sha256": digest(artifact),
            "response_estimator": raw.get("response_estimator"),
            "floquet_status": raw["floquet_convergence"]["status"],
        }
        print(f"{name}: median {np.median(samples):.3f} s; "
              f"G_s(-8 MHz)={operating['G_s']:.9g}; "
              f"peak={np.max(raw['G_s']):.9g}", flush=True)
    (target / "capture.json").write_text(
        json.dumps(report, indent=2, ensure_ascii=False, default=json_default)
        + "\n", encoding="utf-8")
    print("CAPTURE COMPLETE", flush=True)


if __name__ == "__main__":
    capture()
