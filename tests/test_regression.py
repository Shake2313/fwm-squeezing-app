"""
Phase-0 regression: the refactored gabes FWM scheme must reproduce the frozen
known-good baseline captured from the pre-refactor fwm_obe.py.

    python tests/test_regression.py      # or: pytest tests/test_regression.py

Baseline file: tests/baseline_focused.npz  (see capture_baseline.py).
"""
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from gabes.schemes import fwm  # noqa: E402

BASELINE = Path(__file__).resolve().parent / "baseline_focused.npz"
WINDOW_GHZ = 0.55
CONFIGS = {
    "sim_optimum": dict(D_GHz=0.9, T=394.15, P_pump=0.6, P_probe=8e-6,
                        line_strength=0.05, loss_frac=0.055),
    "detuned": dict(D_GHz=1.5, T=383.15, P_pump=0.4, P_probe=10e-6,
                    line_strength=0.05, loss_frac=0.0),
}
KEYS = ("probe_axis_GHz", "G_s", "G_c", "S_dB")


def _spectrum(cfg):
    center = fwm.branch_center_GHz(cfg["D_GHz"], -1)
    return fwm.compute_spectrum(
        cfg["D_GHz"], T=cfg["T"], P_pump=cfg["P_pump"], P_probe=cfg["P_probe"],
        line_strength=cfg["line_strength"], loss_frac=cfg["loss_frac"],
        coarse_points=121, fine_points=0,
        scan_min=center - WINDOW_GHZ, scan_max=center + WINDOW_GHZ,
        velocity_step=5.0, velocity_cutoff=3.0, branches=fwm.BRANCHES,
    )


def test_regression():
    base = np.load(BASELINE)
    for name, cfg in CONFIGS.items():
        spec = _spectrum(cfg)
        for key in KEYS:
            ref = base[f"{name}__{key}"]
            assert np.allclose(spec[key], ref, rtol=1e-9, atol=1e-12), \
                f"{name}/{key} drifted from baseline"


if __name__ == "__main__":
    test_regression()
    print("Phase-0 regression OK - gabes FWM reproduces the baseline.")
