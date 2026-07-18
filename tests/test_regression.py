"""
Phase-0 regression: the refactored gabes FWM scheme must reproduce the frozen
single-branch FWM baseline.

    python tests/test_regression.py      # or: pytest tests/test_regression.py

Baseline file: tests/baseline_focused.npz  (see capture_baseline.py).
"""
import sys
from pathlib import Path
from unittest.mock import patch

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from gabes.schemes import fwm  # noqa: E402
from gabes import constants, hyperfine  # noqa: E402

BASELINE = Path(__file__).resolve().parent / "baseline_focused.npz"
WINDOW_GHZ = 0.55
CONFIGS = {
    "sim_optimum": dict(D_GHz=0.9, T=394.15, P_pump=0.6, P_probe=8e-6,
                        line_strength=0.74, loss_frac=0.055),
    "detuned": dict(D_GHz=1.5, T=383.15, P_pump=0.4, P_probe=10e-6,
                    line_strength=0.74, loss_frac=0.0),
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


def test_seeded_doppler_interpolation_geometry_built_once():
    original = fwm.doppler.interpolation_weights
    center = fwm.branch_center_GHz(0.9, -1)
    with patch.object(fwm.doppler, "interpolation_weights",
                      wraps=original) as counted:
        spec = fwm.compute_spectrum(
            0.9, coarse_points=3, fine_points=0,
            scan_min=center - 0.01, scan_max=center + 0.01,
            velocity_step=200.0, velocity_cutoff=0.1, branch=-1)

    assert counted.call_count == 1
    assert np.all(np.isfinite(spec["G_s"]))


def test_seeded_coupling_factorization_preserves_legacy_api_bitwise():
    center = fwm.branch_center_GHz(0.9, -1)
    common = dict(
        T=394.15, P_pump=0.6, P_probe=8e-6, line_strength=0.74,
        coarse_points=3, fine_points=0,
        scan_min=center - 0.01, scan_max=center + 0.01,
        velocity_step=200.0, velocity_cutoff=0.1, branch=-1,
    )
    legacy = fwm.compute_spectrum(0.9, **common)
    factorized = fwm.compute_spectrum(
        0.9, **common, mode_overlap_penalty=1.0,
        polarization_penalty=1.0, zeeman_participation_penalty=1.0)

    for key in KEYS:
        assert np.array_equal(legacy[key], factorized[key])
    assert factorized["line_strength_residual"] == 0.74
    assert factorized["combined_line_strength_residual"] == 0.74

    scaled = fwm.compute_spectrum(
        0.9, **common, mode_overlap_penalty=0.9,
        polarization_penalty=0.75, zeeman_participation_penalty=0.6)
    expected_residual = 0.74 * 0.9 * 0.75 * 0.6
    assert scaled["lab_coupling_factor"] == 0.9 * 0.75 * 0.6
    assert scaled["combined_line_strength_residual"] == expected_residual
    assert scaled["effective_line_strength"] == (
        expected_residual * scaled["coupling_norm"])


def test_seeded_coupling_factors_multiply_without_hidden_rescaling():
    factors = fwm.SeededCouplingFactors(
        reference_residual=0.8,
        mode_overlap_penalty=0.9,
        polarization_penalty=0.75,
        zeeman_participation_penalty=0.6,
    )
    assert factors.lab_factor == 0.9 * 0.75 * 0.6
    assert factors.combined_residual == 0.8 * factors.lab_factor


def test_seeded_factorized_and_collapsed_coupling_match_in_ultra():
    center = fwm.branch_center_GHz(0.9, -1)
    factors = fwm.SeededCouplingFactors(
        reference_residual=0.74,
        mode_overlap_penalty=0.9,
        polarization_penalty=0.75,
        zeeman_participation_penalty=0.6,
    )
    common = dict(
        T=394.15, P_pump=0.6, P_probe=8e-6,
        coarse_points=3, fine_points=0,
        scan_min=center - 0.01, scan_max=center + 0.01,
        velocity_step=200.0, velocity_cutoff=0.1, branch=-1,
    )
    factor_kwargs = dict(
        line_strength=factors.reference_residual,
        mode_overlap_penalty=factors.mode_overlap_penalty,
        polarization_penalty=factors.polarization_penalty,
        zeeman_participation_penalty=factors.zeeman_participation_penalty,
    )

    for phase_detail in (fwm.PHASE_BALANCED, fwm.PHASE_ULTRA):
        factorized = fwm.compute_spectrum(
            0.9, **common, **factor_kwargs, phase_detail=phase_detail)
        collapsed = fwm.compute_spectrum(
            0.9, **common, line_strength=factors.combined_residual,
            phase_detail=phase_detail)
        for key in KEYS:
            assert np.array_equal(factorized[key], collapsed[key])
        assert factorized["ultra_spatial_overlap_min"] == (
            collapsed["ultra_spatial_overlap_min"])
        assert factorized["effective_line_strength"] == (
            collapsed["effective_line_strength"])


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


def test_seeded_sideband_beat_matches_tpd_convention():
    delta = 2 * np.pi * 12.0e6
    assert np.isclose(
        fwm.seeded_sideband_beat(delta, -1),
        constants.OMEGA_HF - delta)
    assert np.isclose(
        fwm.seeded_sideband_beat(delta, +1),
        constants.OMEGA_HF + delta)


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
    for key in ("opd", "tpd", "temp_c", "pump_mw", "probe_uw", "loss_pct",
                "line_strength", "mode_overlap_penalty", "polarization_penalty",
                "zeeman_participation_penalty"):
        assert defaults[key] == sim[key]
    assert defaults["line_strength"] == 0.74
    assert defaults["mode_overlap_penalty"] == 1.0
    assert defaults["polarization_penalty"] == 1.0
    assert defaults["zeeman_participation_penalty"] == 1.0


if __name__ == "__main__":
    test_regression()
    test_seeded_doppler_interpolation_geometry_built_once()
    test_seeded_coupling_factorization_preserves_legacy_api_bitwise()
    test_seeded_coupling_factors_multiply_without_hidden_rescaling()
    test_seeded_factorized_and_collapsed_coupling_match_in_ultra()
    test_branch_summation_rejected()
    test_seeded_sideband_beat_matches_tpd_convention()
    test_fwm_uses_rb85_d1_hyperfine_strengths()
    test_fwm_defaults_match_sim_preset()
    print("Phase-0 regression OK - gabes FWM reproduces the baseline.")
