"""
Phase-0 regression: the refactored gabes FWM scheme must reproduce the frozen
single-branch FWM baseline.

    python tests/test_regression.py      # or: pytest tests/test_regression.py

Baseline file: tests/baseline_focused.npz  (see capture_baseline.py).
"""
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from gabes.schemes import fwm  # noqa: E402
from gabes import constants, hyperfine  # noqa: E402

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
        velocity_step=5.0, velocity_cutoff=3.0, branch=-1,
    )


def test_regression():
    base = np.load(BASELINE)
    for name, cfg in CONFIGS.items():
        spec = _spectrum(cfg)
        for key in KEYS:
            ref = base[f"{name}__{key}"]
            assert np.allclose(spec[key], ref, rtol=1e-9, atol=1e-12), \
                f"{name}/{key} drifted from baseline"


def test_branch_summation_rejected():
    center = fwm.branch_center_GHz(0.9, -1)
    try:
        fwm.compute_spectrum(
            0.9, coarse_points=3, fine_points=0,
            scan_min=center - 0.01, scan_max=center + 0.01,
            velocity_step=20.0, velocity_cutoff=1.0, branches=fwm.BRANCHES)
    except ValueError as exc:
        assert "separate probe/conjugate mode pairs" in str(exc)
    else:
        raise AssertionError("multi-branch susceptibility summation must fail")


def test_fwm_uses_rb85_d1_hyperfine_strengths():
    mapping = {
        (fwm.G1, fwm.E2): (2, 2),
        (fwm.G1, fwm.E3): (2, 3),
        (fwm.G2, fwm.E2): (3, 2),
        (fwm.G2, fwm.E3): (3, 3),
    }
    for (g, e), key in mapping.items():
        expected = 3.0 * hyperfine.CF2[key]
        assert np.isclose(fwm.TRANSITION_DIPOLE_SCALE[g, e] ** 2, expected)

    decay = {(e, g): rate for e, g, rate in fwm.ATOM.decay}
    for e, fe in ((fwm.E2, 2), (fwm.E3, 3)):
        total = decay[(e, fwm.G1)] + decay[(e, fwm.G2)]
        assert np.isclose(total, constants.GAMMA)
        w2 = hyperfine.CF2[(2, fe)]
        w3 = hyperfine.CF2[(3, fe)]
        assert np.isclose(decay[(e, fwm.G1)] / total, w2 / (w2 + w3))
        assert np.isclose(decay[(e, fwm.G2)] / total, w3 / (w2 + w3))


def test_fwm_defaults_match_sim_preset():
    scheme = fwm.FWMScheme()
    defaults = scheme.defaults()
    sim = scheme.recommended_defaults(defaults)[fwm.MODE_SEEDED]
    for key in ("opd", "tpd", "temp_c", "pump_mw", "probe_uw", "loss_pct"):
        assert defaults[key] == sim[key]


if __name__ == "__main__":
    test_regression()
    test_branch_summation_rejected()
    test_fwm_uses_rb85_d1_hyperfine_strengths()
    test_fwm_defaults_match_sim_preset()
    print("Phase-0 regression OK - gabes FWM reproduces the baseline.")
