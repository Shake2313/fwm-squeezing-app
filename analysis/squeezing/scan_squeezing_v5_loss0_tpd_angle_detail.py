"""Detailed v5 two-photon-detuning and crossing-angle diagnostic scan."""
from __future__ import annotations

import os

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")

import numpy as np

try:
    from .scan_squeezing_tpd_angle_detail import ScanConfig, run
except ImportError:
    from scan_squeezing_tpd_angle_detail import ScanConfig, run


CONFIG = ScanConfig(
    version="v5",
    delta_axis_ghz=np.round(np.arange(-2.06, -1.94 + 1e-9, 0.02), 3),
    temperature_axis_c=np.round(np.arange(106.0, 124.0 + 1e-9, 1.0), 3),
    angle_axis_deg=np.round(np.arange(0.0, 0.36 + 1e-9, 0.02), 3),
    two_photon_axis_mhz=np.round(np.linspace(-500.0, -180.0, 81), 3),
    reference_delta_ghz=-2.0,
    reference_temperature_c=120.0,
)


if __name__ == "__main__":
    run(CONFIG)
