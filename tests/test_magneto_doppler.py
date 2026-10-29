"""
Unregistered single-pair reference: exact Maxwell light-region average.

  Exactness : the closed-form average matches a converged velocity grid solved
              by the unchanged steady-state kernels (forced guard fallback).
  Off path  : doppler="off" never builds a pole system.
  Fallback  : the grid step resolves the homogeneous optical Lorentzian.
  Provenance: raw keys and the Derived table name the Doppler method.

    pytest tests/test_magneto_doppler.py
"""
import sys
from pathlib import Path
from unittest.mock import patch

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from gabes import constants, doppler, pole_doppler  # noqa: E402
from gabes.schemes import magneto  # noqa: E402
from gabes.schemes.magneto import SinglePairMagnetoReference  # noqa: E402

FIELDS = ("chi_probe", "chi_p", "chi_m")


def _params(regime, **over):
    params = SinglePairMagnetoReference().defaults()
    for regimes in SinglePairMagnetoReference._REGIMES.values():
        if regime in regimes:
            params.update(regimes[regime])
    params.update(over)
    return params


CASES = {
    "paraffin EIT dip": _params("EIT dip", scan_points=5),
    "paraffin detuned": _params("EIT dip", scan_points=5, laser_detuning_mhz=200.0),
    # No transverse field: the Doppler-free Zeeman populations have no relaxation
    # of their own, so only the shift-invert pole form applies.
    "paraffin F1-F2 no transverse field": _params(
        "EIT dip", scan_points=5, transition=magneto._transition_label(1, 2),
        residual_transverse_b_ut=0.0),
    "buffer LCA F2-F2": _params("Buffer LCA", scan_points=5),
}


def _relative(a, b):
    return float(np.max(np.abs(a - b)) / max(float(np.max(np.abs(b))), 1e-300))


@pytest.mark.parametrize("name", list(CASES))
def test_exact_average_matches_a_converged_velocity_grid(name):
    params = CASES[name]
    exact = SinglePairMagnetoReference().compute(params)
    with patch.object(magneto, "POLE_GUARD_RTOL", -1.0):
        grid = SinglePairMagnetoReference().compute(params)

    assert exact["doppler_method"] == magneto.DOPPLER_EXACT
    assert exact["doppler_fallback_rows"] == 0
    assert exact["doppler_guard_max_relative"] <= magneto.POLE_GUARD_RTOL
    assert exact["doppler_min_pole_distance_mhz"] > 0.0
    assert grid["doppler_fallback_rows"] == params["scan_points"]
    assert grid["doppler_fallback_classes"] > 50
    for key in FIELDS:
        assert _relative(exact[key], grid[key]) <= 1e-5, key


def test_doppler_off_never_builds_a_pole_system():
    params = _params("EIT dip", scan_points=5, doppler="off")
    with patch.object(pole_doppler, "ShiftInvertSystem",
                      side_effect=AssertionError("pole path used")), \
            patch.object(pole_doppler, "AffineShiftSystem",
                         side_effect=AssertionError("pole path used")):
        raw = SinglePairMagnetoReference().compute(params)
    assert raw["doppler_method"] == magneto.DOPPLER_OFF
    assert raw["doppler_guard_max_relative"] is None
    assert raw["doppler_fallback_rows"] == 0


def test_fallback_grid_resolves_the_homogeneous_lorentzian():
    T = 298.15
    gamma_hwhm = 0.5 * magneto.GAMMA_D1
    v, wt = magneto._maxwell_fallback_grid(T, magneto.RB87.mass, magneto.K_D1_RB87,
                                           gamma_hwhm)
    step = magneto.K_D1_RB87 * (v[1] - v[0])
    assert np.exp(-2 * np.pi * gamma_hwhm / step) <= (
        magneto.FALLBACK_QUADRATURE_ERROR * (1 + 1e-9))
    sigma_v = np.sqrt(constants.KB * T / magneto.RB87.mass)
    assert v[-1] >= magneto.FALLBACK_SIGMAS * sigma_v
    assert np.isclose(wt.sum(), 1.0, rtol=0.0, atol=1e-12)

    # The grid average of the unit-area Lorentzian reproduces the Voigt profile.
    sigma_d = magneto.K_D1_RB87 * sigma_v
    for detuning in (0.0, 0.7 * sigma_d, 3.0 * sigma_d):
        arg = detuning - magneto.K_D1_RB87 * v
        grid = np.sum(wt * (gamma_hwhm / np.pi) / (arg ** 2 + gamma_hwhm ** 2))
        voigt = doppler.voigt_profile(detuning, 2 * gamma_hwhm, sigma_d)
        assert abs(grid / voigt - 1.0) <= 1e-5


def test_exact_average_is_reported_and_replaces_the_class_knob():
    params = _params("EIT dip", scan_points=5)
    scheme = SinglePairMagnetoReference()
    raw = scheme.compute(params)
    view = scheme.headless_observables(raw, params)
    derived = next(table["markdown"] for table in view["tables"]
                   if table["title"] == "Derived quantities")
    assert magneto.DOPPLER_EXACT in derived
    assert "Grid-fallback B rows" in derived
    assert raw["doppler_scale"] == 1.0 and raw["N_eff"] == raw["N"]
    assert "velocity_classes" not in {spec.name for spec in scheme.param_schema()}


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-q"]))
