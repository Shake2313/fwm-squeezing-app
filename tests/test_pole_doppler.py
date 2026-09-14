"""Exact pole–residue Doppler reduction, Faddeeva/Voigt helpers, adaptive scan sampling."""
import math

import numpy as np
import pytest

from gabes import adaptive_scan, core, doppler, hyperfine, kernels, pole_doppler
from gabes.core import build_liouvillian, comm_super
from gabes.schemes import fwm


def test_faddeeva_matches_scipy_in_the_upper_half_plane():
    special = pytest.importorskip("scipy.special")
    rng = np.random.default_rng(7)
    z = np.concatenate([
        rng.uniform(-60, 60, 3000) + 1j * 10 ** rng.uniform(-8, 2, 3000),
        rng.uniform(-3000, 3000, 500) + 1j * 10 ** rng.uniform(-8, 0, 500),
        np.linspace(-6, 6, 601) + 0j,
    ])
    reference = special.wofz(z)
    relative = np.abs(doppler.faddeeva_upper(z) - reference) / np.abs(reference)
    assert np.max(relative) < 1e-11


def test_maxwell_resolvent_mean_matches_quadrature_in_both_half_planes():
    lam = np.array([0.4 + 0.02j, -1.7 - 0.03j, 6.0 + 0.001j, 0.1 - 0.5j])
    x = np.linspace(-14.0, 14.0, 1_400_001)
    density = np.exp(-0.5 * x**2) / math.sqrt(2.0 * math.pi)
    dx = x[1] - x[0]
    brute = np.array([np.sum(density / (value - x)) * dx for value in lam])
    # Shift/scale the Gaussian to exercise center and sigma.
    got = doppler.maxwell_resolvent_mean(2.0 * lam + 3.0, 3.0, 2.0) * 2.0
    np.testing.assert_allclose(got, brute, rtol=1e-9)


def test_voigt_profile_has_unit_area_and_gaussian_limit():
    sigma = 1.3
    x = np.linspace(-400.0, 400.0, 2_000_001)
    area = np.sum(doppler.voigt_profile(x, 0.02, sigma)) * (x[1] - x[0])
    assert abs(area - 1.0) < 1e-4          # Lorentzian wings beyond ±400 are ~3e-5
    # Gaussian limit; keep the Lorentzian wing γ/(2πx²) far below the Gaussian at 6σ.
    core_x = np.linspace(-6.0 * sigma, 6.0 * sigma, 241)
    gauss = np.exp(-0.5 * (core_x / sigma) ** 2) / (sigma * math.sqrt(2 * math.pi))
    # Weideman's expansion has a ~1e-14 absolute floor on Re w in the far wing.
    np.testing.assert_allclose(doppler.voigt_profile(core_x, 1e-14, sigma), gauss,
                               rtol=1e-6, atol=1e-13)


def test_voigt_profile_matches_scipy():
    special = pytest.importorskip("scipy.special")
    x = np.linspace(-80.0, 80.0, 4001)
    for gamma_fwhm, sigma in ((0.02, 1.3), (1.0, 0.4), (5.0, 2.0)):
        np.testing.assert_allclose(
            doppler.voigt_profile(x, gamma_fwhm, sigma),
            special.voigt_profile(x, sigma, 0.5 * gamma_fwhm), rtol=1e-8, atol=1e-15)


def test_cubic_spline_matches_scipy_not_a_knot():
    interpolate = pytest.importorskip("scipy.interpolate")
    rng = np.random.default_rng(11)
    x = np.sort(rng.uniform(0.0, 10.0, 29))
    x[0], x[-1] = 0.0, 10.0
    y = rng.standard_normal((29, 3)) + 1j * rng.standard_normal((29, 3))
    q = np.linspace(0.0, 10.0, 777)
    reference = (interpolate.CubicSpline(x, y.real)(q)
                 + 1j * interpolate.CubicSpline(x, y.imag)(q))
    np.testing.assert_allclose(adaptive_scan.cubic_spline(x, y, q), reference,
                               rtol=0, atol=1e-11)


def test_cubic_spline_reproduces_cubics_and_its_nodes():
    x = np.array([-2.0, -1.1, 0.0, 0.3, 1.7, 2.5, 4.0])
    f = lambda t: (1.0 - 2.0 * t + 0.5 * t**2 - 0.1 * t**3) * (1.0 + 0.5j)
    q = np.linspace(-2.0, 4.0, 101)
    np.testing.assert_allclose(adaptive_scan.cubic_spline(x, f(x), q), f(q),
                               rtol=1e-12, atol=1e-12)
    np.testing.assert_allclose(adaptive_scan.cubic_spline(x, f(x), x), f(x), atol=1e-13)


def test_refine_scan_nodes_concentrates_on_a_sharp_feature():
    x = np.linspace(-1.0, 1.0, 401)
    f = 1.0 / (x - 0.3137 + 0.004j)
    known = {}

    def solve(rows):
        for r in rows:
            known[int(r)] = f[int(r)]

    def disagrees(nodes, candidates):
        predicted = adaptive_scan.cubic_spline(
            x[nodes], np.array([known[int(n)] for n in nodes]), x[candidates])
        return np.abs(predicted - f[candidates]) > 1e-2 * np.abs(f[candidates]) + 1.0

    nodes, rounds = adaptive_scan.refine_scan_nodes(401, solve, disagrees, start_stride=32)
    assert nodes[0] == 0 and nodes[-1] == 400
    assert set(known) == set(int(n) for n in nodes)
    assert nodes.size < 120 and rounds >= 3
    spacing = np.diff(x[nodes])
    near = np.abs(x[nodes][:-1] - 0.3137) < 0.02
    assert spacing[near].max() < spacing.max()           # refinement went to the pole
    full = adaptive_scan.cubic_spline(x[nodes], f[nodes], x)
    assert np.max(np.abs(full - f)) < 0.05 * np.max(np.abs(f))


def _seeded_systems(T, D_GHz, branch, order):
    Op = fwm.rabi_freq(0.6, fwm.W_PUMP)
    Os = fwm.rabi_freq(8e-6, fwm.W_PROBE)
    atom = fwm.collisional_atom(T, hyperfine.number_density(T))
    systems = fwm._seeded_pole_systems(Op, Op, Os, Os, branch, atom, (order,))
    assert systems is not None
    return Op, Os, atom, systems


@pytest.mark.parametrize("branch", (-1, 1))
@pytest.mark.parametrize("order", (2, 3))
def test_floquet_real_form_matches_the_compiled_kernel(branch, order):
    T, D = 408.15, -2.2
    Op, Os, atom, systems = _seeded_systems(T, D, branch, order)
    center = fwm.branch_center_GHz(D, branch)
    delta = fwm.two_photon_detuning_from_probe_scan(
        np.linspace(center - 0.5, center + 0.5, 5), D, branch)
    Delta = 2 * np.pi * D * 1e9
    D_eff = np.array([Delta - 1.3e9, Delta, Delta + 0.7e9])
    tables = fwm.chi_matrix_table(Op, Op, Os, Os, delta, D_eff, branch, atom=atom, n_f=order)
    for seed, columns in ((1, (0, 1)), (2, (2, 3))):
        system, reference = systems[(seed, order)]
        poles = system.residues(delta)
        for j, value in enumerate(D_eff):
            got = poles.value_at(value) / reference
            for readout, column in enumerate(columns):
                np.testing.assert_allclose(got[:, readout], tables[column][:, j], rtol=1e-9)


def test_pole_maxwell_average_is_exact_against_a_fine_velocity_grid():
    T, D, branch, order = 394.15, 0.9, -1, 3
    Op, Os, atom, systems = _seeded_systems(T, D, branch, order)
    center = fwm.branch_center_GHz(D, branch)
    delta = fwm.two_photon_detuning_from_probe_scan(
        np.array([center - 0.3, center + 0.04, center + 0.4]), D, branch)
    Delta = 2 * np.pi * D * 1e9
    sigma = fwm.K_VEC * math.sqrt(fwm.constants.KB * T / fwm.constants.MASS_85RB)
    v, w = doppler.velocity_grid(T, dv=0.5, cutoff_sigma=8.0)
    tables = fwm.chi_matrix_table(Op, Op, Os, Os, delta, Delta - fwm.K_VEC * v, branch,
                                  atom=atom, n_f=order)
    for seed, columns in ((1, (0, 1)), (2, (2, 3))):
        system, reference = systems[(seed, order)]
        exact = system.residues(delta).gaussian_mean(Delta, sigma) / reference
        for readout, column in enumerate(columns):
            np.testing.assert_allclose(exact[:, readout], tables[column] @ w, rtol=1e-8)


def test_discrete_mean_reproduces_the_interpolated_grid_average():
    T, D, branch, order = 394.15, 0.9, -1, 2
    Op, Os, atom, systems = _seeded_systems(T, D, branch, order)
    center = fwm.branch_center_GHz(D, branch)
    delta = fwm.two_photon_detuning_from_probe_scan(
        np.linspace(center - 0.5, center + 0.5, 4), D, branch)
    Delta = 2 * np.pi * D * 1e9
    v, w = doppler.velocity_grid(T, dv=25.0, cutoff_sigma=4.0)
    axis = doppler.build_Delta_eff_axis(Delta, Delta, v)
    idx_lo, frac = doppler.interpolation_weights(axis, Delta, v)
    nodes = np.zeros(axis.size)
    np.add.at(nodes, idx_lo, w * (1.0 - frac))
    np.add.at(nodes, idx_lo + 1, w * frac)
    tables = fwm.chi_matrix_table(Op, Op, Os, Os, delta, axis, branch, atom=atom, n_f=order)
    for seed, columns in ((1, (0, 1)), (2, (2, 3))):
        system, reference = systems[(seed, order)]
        mean = system.residues(delta).discrete_mean(axis, nodes) / reference
        for readout, column in enumerate(columns):
            expected = doppler.apply_doppler_average(tables[column], idx_lo, frac, w)
            np.testing.assert_allclose(mean[:, readout], expected, rtol=1e-9)


def test_numpy_fallback_matches_the_compiled_stages(monkeypatch):
    if not kernels.available():
        pytest.skip("compiled path unavailable")
    T, D, branch, order = 394.15, 0.9, 1, 3
    _, _, _, systems = _seeded_systems(T, D, branch, order)
    system, _ = systems[(2, order)]
    delta = np.linspace(-2.0e9, 2.0e9, 7)
    nodes = np.linspace(-3e9, 3e9, 301)
    weights = np.exp(-0.5 * (nodes / 1.5e9) ** 2)
    weights /= weights.sum()
    compiled = system.residues(delta).discrete_mean(nodes, weights)
    monkeypatch.setattr(kernels, "available", lambda: False)
    fallback = system.residues(delta).discrete_mean(nodes, weights)
    np.testing.assert_allclose(fallback, compiled, rtol=1e-9)


def test_affine_shift_system_validates_its_inputs():
    eye = np.eye(3)
    with pytest.raises(TypeError):
        pole_doppler.AffineShiftSystem(eye + 0j, eye, eye, np.ones(3), np.ones((1, 3)))
    with pytest.raises(ValueError):
        pole_doppler.AffineShiftSystem(eye, eye, np.zeros((3, 3)), np.ones(3), np.ones((1, 3)))
