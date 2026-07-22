"""
Pin the numba fast path (kernels.floquet_chi_grid) to the NumPy reference
implementation of the FWM χ̄ table.

    python tests/test_kernels.py      # or: pytest tests/test_kernels.py

The compiled LU uses LAPACK's pivot rule, so agreement is expected at the
~1e-13 relative level — far inside the 1e-9 regression tolerance.
"""
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from gabes import constants, doppler, kernels  # noqa: E402
from gabes.constants import rabi_freq  # noqa: E402
from gabes.schemes import fwm  # noqa: E402

import pytest  # noqa: E402

pytestmark = pytest.mark.skipif(
    not kernels.available(), reason="numba not available")


def _chi_tables(branch):
    Op = rabi_freq(0.6, fwm.W_PUMP)
    Os = rabi_freq(8e-6, fwm.W_PROBE)
    D_GHz = 0.9
    Delta = 2 * np.pi * D_GHz * 1e9

    center = fwm.branch_center_GHz(D_GHz, branch)
    probe_GHz = np.linspace(center - 0.4, center + 0.4, 41)
    delta_axis = fwm.two_photon_detuning_from_probe_scan(probe_GHz, D_GHz, branch)

    v_grid, _ = doppler.velocity_grid(394.15, dv=12.0, cutoff_sigma=3.0)
    deff_axis = doppler.build_Delta_eff_axis(Delta, Delta, v_grid)

    fast = fwm.chi_matrix_table(Op, Os, Os, delta_axis, deff_axis, branch)

    orig = fwm.kernels.available
    fwm.kernels.available = lambda: False
    try:
        ref = fwm.chi_matrix_table(Op, Os, Os, delta_axis, deff_axis, branch)
    finally:
        fwm.kernels.available = orig
    return fast, ref


@pytest.mark.parametrize("branch", [-1, +1])
def test_kernel_matches_reference(branch):
    fast, ref = _chi_tables(branch)
    names = ("chi_ss", "chi_cs", "chi_sc", "chi_cc")
    for name, a, b in zip(names, fast, ref):
        scale = np.abs(b).max()
        worst = np.abs(a - b).max() / scale
        assert worst < 1e-10, f"branch {branch} {name}: rel diff {worst:.3e}"


# --- magneto (P2): the kernel path must match the NumPy solver path ----------
from gabes.schemes.magneto import MagnetoScheme  # noqa: E402

_MAG_BASE = dict(
    signal_type="Transmission", cell_type="Paraffin coated cell",
    transition="F=2 → F'=1", intensity_mw_cm2=0.8, qwp_deg=0.0, b_max_ut=2.0,
    laser_detuning_mhz=0.0, temp_c=25.0, cell_mm=10.0, ne_pressure_torr=10.0,
    wall_coherence_ms=3.0, transit_relax_khz=80.0, dark_return_khz=1.0,
    buffer_ground_relax_khz=3.0, collisional_depol_khz=0.5,
    residual_transverse_b_ut=0.05, transverse_field_angle_deg=0.0,
    line_strength=1.0, doppler="on", velocity_classes=9, scan_points=61)

_MAG_CASES = {
    "paraffin Fe1": dict(_MAG_BASE),
    "paraffin Fe2": {**_MAG_BASE, "transition": "F=2 → F'=2"},
    "buffer Fe1": {**_MAG_BASE, "cell_type": "Buffer gas cell"},
    "nmor": {**_MAG_BASE, "signal_type": "NMOR rotation"},
}


def _magneto_with(numba_on, params):
    orig = kernels.NUMBA_AVAILABLE
    kernels.NUMBA_AVAILABLE = numba_on
    try:
        return MagnetoScheme().compute(params)
    finally:
        kernels.NUMBA_AVAILABLE = orig


@pytest.mark.parametrize("name", list(_MAG_CASES))
def test_magneto_kernel_matches_reference(name):
    params = _MAG_CASES[name]
    fast = _magneto_with(True, params)
    ref = _magneto_with(False, params)
    for key in ("chi_probe", "chi_p", "chi_m"):
        scale = max(np.abs(ref[key]).max(), 1e-300)
        worst = np.abs(fast[key] - ref[key]).max() / scale
        assert worst < 1e-9, f"{name}/{key}: rel diff {worst:.3e}"


# --- Λ / Rydberg (P3): affine-scan kernel must match the NumPy loop ----------
from gabes.schemes.absorption import LambdaScheme  # noqa: E402
from gabes.schemes.rydberg import RydbergEITScheme  # noqa: E402


def _scheme_with(numba_on, scheme, params):
    orig = kernels.NUMBA_AVAILABLE
    kernels.NUMBA_AVAILABLE = numba_on
    try:
        return scheme.compute(params)
    finally:
        kernels.NUMBA_AVAILABLE = orig


def _affine_cases():
    cases = []
    for mode in ("eit", "at", "cpt"):
        s = LambdaScheme(mode)
        d = s._default_values(mode)
        cases.append((f"{mode}-off", s, {**d, "doppler": "off"}))
        cases.append((f"{mode}-on", s, {**d, "doppler": "on", "temp_c": 50.0}))
    ryd = RydbergEITScheme()
    cases.append(("ryd-at", ryd, ryd._defaults("AT electrometry")))
    cases.append(("ryd-eit", ryd, ryd._defaults("EIT")))
    return cases


@pytest.mark.parametrize("case", _affine_cases(), ids=lambda c: c[0])
def test_affine_scan_kernel_matches_reference(case):
    _name, scheme, params = case
    fast = _scheme_with(True, scheme, params)["chi_bar"]
    ref = _scheme_with(False, scheme, params)["chi_bar"]
    scale = max(np.abs(ref).max(), 1e-300)
    worst = np.abs(fast - ref).max() / scale
    assert worst < 1e-11, f"rel diff {worst:.3e}"


if __name__ == "__main__":
    for br in (-1, +1):
        test_kernel_matches_reference(br)
    print("FWM kernel == reference (rel < 1e-10) on both branches")
    for nm in _MAG_CASES:
        test_magneto_kernel_matches_reference(nm)
    print("magneto kernel == reference (rel < 1e-9) on all cases")
    for c in _affine_cases():
        test_affine_scan_kernel_matches_reference(c)
    print("affine-scan kernel == reference (rel < 1e-11) on Λ/Rydberg")
