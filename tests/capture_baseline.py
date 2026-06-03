"""
Capture a known-good FWM regression baseline from the CURRENT single-branch
FWM model. Run after intentional physics changes; the resulting .npz is the
frozen anchor that the refactored gabes FWM scheme must reproduce
(see test_regression.py).

    python tests/capture_baseline.py
"""
import os
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")          # no display needed
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from gabes.schemes import fwm

# Two focused configs mirroring the app's "Fast" focused_spectrum settings.
WINDOW_GHZ = 0.55
CONFIGS = {
    "sim_optimum": dict(D_GHz=0.9, T=394.15, P_pump=0.6, P_probe=8e-6,
                        line_strength=0.05, loss_frac=0.055),
    "detuned": dict(D_GHz=1.5, T=383.15, P_pump=0.4, P_probe=10e-6,
                    line_strength=0.05, loss_frac=0.0),
}


def run(cfg):
    center = fwm.branch_center_GHz(cfg["D_GHz"], -1)
    return fwm.compute_spectrum(
        cfg["D_GHz"], T=cfg["T"], P_pump=cfg["P_pump"], P_probe=cfg["P_probe"],
        line_strength=cfg["line_strength"], loss_frac=cfg["loss_frac"],
        coarse_points=121, fine_points=0,
        scan_min=center - WINDOW_GHZ, scan_max=center + WINDOW_GHZ,
        velocity_step=5.0, velocity_cutoff=3.0, branch=-1,
    )


def main():
    saved = {}
    for name, cfg in CONFIGS.items():
        spec = run(cfg)
        saved[f"{name}__probe_axis_GHz"] = spec["probe_axis_GHz"]
        saved[f"{name}__G_s"] = spec["G_s"]
        saved[f"{name}__G_c"] = spec["G_c"]
        saved[f"{name}__S_dB"] = spec["S_dB"]
        print(f"[{name}]  max G_s = {spec['G_s'].max():.4f}"
              f"   min S_dB = {spec['S_dB'].min():.4f}"
              f"   n_points = {spec['probe_axis_GHz'].size}")

    out = Path(__file__).resolve().parent / "baseline_focused.npz"
    np.savez_compressed(out, **saved)
    print(f"\nSaved baseline: {out}")


if __name__ == "__main__":
    main()
