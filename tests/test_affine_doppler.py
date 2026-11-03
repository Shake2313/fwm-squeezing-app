"""Affine-scan Doppler averages: exact Maxwell mean, resolving grids and per-class solves agree."""
import numpy as np
import pytest

from gabes import affine_doppler, atoms, constants, core, doppler, kernels, pole_doppler, schemes
from gabes.affine_doppler import AffineScan

T, MASS, K = 330.0, constants.MASS_85RB, constants.K_VEC
SIGMA = K * np.sqrt(constants.KB * T / MASS)


def _lambda_scan():
    atom = atoms.lambda3(two_photon_doppler_ratio=0.02)
    probe, coupling = 1e-3 * constants.GAMMA, 0.4 * constants.GAMMA

    def build_H0(s):
        H = np.zeros((3, 3), dtype=complex)
        H[1, 1] = -s
        H[0, 2] = H[2, 0] = probe / 2
        H[1, 2] = H[2, 1] = coupling / 2
        return H

    return atom, build_H0, probe


def _maxwell_scan():
    atom, build_H0, probe = _lambda_scan()
    affine = AffineScan(atom, build_H0, (2, 0), probe, atom.S_v.copy())
    return affine, np.linspace(-3, 3, 41) * constants.GAMMA


def test_maxwell_strip_grid_covers_the_bulk_at_the_quadrature_step():
    strip = 0.2 * constants.GAMMA
    v, w = doppler.maxwell_strip_grid(T, MASS, K, strip, error=1e-6, cutoff_sigma=6.0)
    step = v[1] - v[0]
    np.testing.assert_allclose(K * step, 2 * np.pi * strip / np.log(1e6), rtol=1e-12)
    assert np.isclose(w.sum(), 1.0)
    assert v[-1] >= 6.0 * SIGMA / K
    np.testing.assert_allclose(np.sum(w * (K * v) ** 2), SIGMA ** 2, rtol=1e-6)


def test_class_sum_is_the_weighted_per_class_steady_state():
    atom, build_H0, probe = _lambda_scan()
    scan = np.linspace(-2, 2, 9) * constants.GAMMA
    kv = np.array([-3.0, 0.0, 1.5]) * constants.GAMMA
    w = np.array([0.2, 0.5, 0.3])
    shift = atom.S_v.copy()
    got = AffineScan(atom, build_H0, (2, 0), probe, shift)._class_sum(scan, kv, w)
    expected = np.array([
        sum(wj * core.steady_state_batched(
            core.build_liouvillian(build_H0(s), atom) - s * shift + kj * atom.S_v,
            np.zeros(1), atom.S_v, 3)[0, 2, 0] for kj, wj in zip(kv, w)) / probe
        for s in scan])
    np.testing.assert_allclose(got, expected, rtol=1e-9, atol=1e-12)


def test_stationary_matches_the_single_zero_class():
    atom, build_H0, probe = _lambda_scan()
    scan = np.linspace(-1, 1, 7) * constants.GAMMA
    affine = AffineScan(atom, build_H0, (2, 0), probe, atom.S_v)
    np.testing.assert_allclose(affine.stationary(scan),
                               affine._class_sum(scan, np.zeros(1), np.ones(1)),
                               rtol=1e-10, atol=1e-13)
    if kernels.available():   # the kernel path is the same call for both
        np.testing.assert_array_equal(affine.stationary(scan),
                                      affine._class_sum(scan, np.zeros(1), np.ones(1)))


def test_stationary_with_the_scan_in_the_hamiltonian_is_the_single_zero_class():
    atom, build_H0, probe = _lambda_scan()
    scan = np.linspace(-1, 1, 7) * constants.GAMMA
    affine = AffineScan(atom, build_H0, (2, 0), probe)
    np.testing.assert_array_equal(affine.stationary(scan),
                                  affine._class_sum(scan, np.zeros(1), np.ones(1)))


def test_maxwell_mean_is_the_limit_of_fine_velocity_grids():
    affine, scan = _maxwell_scan()
    v, w = doppler.maxwell_strip_grid(T, MASS, K, 0.1 * constants.GAMMA,
                                      error=1e-10, cutoff_sigma=8.0)
    np.testing.assert_allclose(affine.maxwell_mean(scan, T, MASS, K),
                               affine._class_sum(scan, K * v, w), rtol=1e-8, atol=0.0)


def test_a_two_metre_per_second_grid_aliases_the_optical_resonance():
    affine, scan = _maxwell_scan()
    exact = affine.maxwell_mean(scan, T, MASS, K)
    v, w = doppler.velocity_grid(T, mass=MASS, dv=2.0, cutoff_sigma=8.0)
    assert np.max(np.abs(affine._class_sum(scan, K * v, w) / exact - 1.0)) > 1e-4


def test_quadrature_classes_resolve_the_narrowest_pole():
    affine, scan = _maxwell_scan()
    v, w = affine.quadrature_classes(scan, T, MASS, K)
    np.testing.assert_allclose(affine._class_sum(scan, K * v, w),
                               affine.maxwell_mean(scan, T, MASS, K), rtol=1e-6, atol=0.0)


def test_strip_classes_stay_inside_their_floor_and_ceiling():
    floored, _ = affine_doppler._strip_classes(np.zeros(1), T, MASS, K)
    capped, _ = affine_doppler._strip_classes(np.full(1, np.inf), T, MASS, K)
    step = 2 * np.pi * affine_doppler.STRIP_FLOOR_SIGMAS * SIGMA
    np.testing.assert_allclose(K * (floored[1] - floored[0]),
                               step / np.log(1.0 / affine_doppler.QUADRATURE_ERROR), rtol=1e-9)
    assert 20 < capped.size < 100


def test_rows_missing_the_guard_are_summed_on_a_resolving_grid(monkeypatch):
    affine, scan = _maxwell_scan()
    exact = affine.maxwell_mean(scan, T, MASS, K)
    pole_rows = AffineScan._pole_rows

    def one_bad_row(self, *args):
        poles, damping, failed = pole_rows(self, *args)
        poles.res[7] = 1e30
        failed[7] = True
        return poles, damping, failed

    monkeypatch.setattr(AffineScan, "_pole_rows", one_bad_row)
    got = affine.maxwell_mean(scan, T, MASS, K)
    assert got[7] != exact[7]
    np.testing.assert_allclose(got, exact, rtol=1e-6, atol=0.0)


def test_guard_rejects_a_row_off_the_direct_solve(monkeypatch):
    affine, scan = _maxwell_scan()
    assert not affine._pole_rows(scan, SIGMA)[2].any()
    monkeypatch.setattr(affine_doppler, "POLE_GUARD_RTOL", -1.0)
    assert affine._pole_rows(scan, SIGMA)[2].all()


def test_guard_rejects_an_undamped_pole_inside_the_bulk(monkeypatch):
    affine, scan = _maxwell_scan()
    residues = pole_doppler.ShiftInvertSystem.residues

    def real_pole(self, *args, **kwargs):
        poles = residues(self, *args, **kwargs)
        poles.lam[3, 0] = 0.5 * SIGMA
        return poles

    monkeypatch.setattr(pole_doppler.ShiftInvertSystem, "residues", real_pole)
    assert affine._pole_rows(scan, SIGMA)[2][3]


def test_singular_expansion_falls_back_to_a_resolving_grid(monkeypatch):
    affine, scan = _maxwell_scan()
    scan = scan[:3]
    exact = affine.maxwell_mean(scan, T, MASS, K)

    def singular(*args, **kwargs):
        raise np.linalg.LinAlgError("Singular matrix")

    monkeypatch.setattr(pole_doppler.ShiftInvertSystem, "residues", singular)
    np.testing.assert_allclose(affine.maxwell_mean(scan, T, MASS, K), exact, rtol=1e-6, atol=0.0)


@pytest.mark.parametrize("name, params", [
    ("lambda", {"view": "EIT"}),
    ("lambda", {"view": "CPT", "species": "⁸⁷Rb", "line": "D2"}),
    ("rydberg_eit", {"view": "EIT"}),
    ("rydberg_eit", {"view": "AT electrometry"}),
])
def test_doppler_on_scheme_matches_a_resolving_velocity_grid(monkeypatch, name, params):
    scheme = schemes.get(name)
    p = scheme.defaults()
    p.update(params, doppler="on")
    class_sum, maxwell_mean = AffineScan._class_sum, AffineScan.maxwell_mean
    fallback, checked = [], []

    def checked_mean(self, scan, T, mass, k_vec):
        out = maxwell_mean(self, scan, T, mass, k_vec)
        rows = np.arange(0, scan.size, 40)
        v, w = self.quadrature_classes(scan[rows], T, mass, k_vec)
        np.testing.assert_allclose(out[rows], class_sum(self, scan[rows], k_vec * v, w),
                                   rtol=1e-6, atol=0.0)
        checked.append(rows.size)
        return out

    def counted_sum(self, scan, kv, w):
        fallback.append(scan.size)
        return class_sum(self, scan, kv, w)

    monkeypatch.setattr(AffineScan, "maxwell_mean", checked_mean)
    monkeypatch.setattr(AffineScan, "_class_sum", counted_sum)
    scheme.compute(p)
    assert checked and not fallback


@pytest.mark.parametrize("name", ["lambda", "rydberg_eit"])
def test_doppler_off_never_reaches_the_maxwell_mean(monkeypatch, name):
    scheme = schemes.get(name)
    p = scheme.defaults()
    p.update(doppler="off")

    def refuse(*args):
        raise AssertionError("Doppler off is the single v = 0 class")

    monkeypatch.setattr(AffineScan, "maxwell_mean", refuse)
    assert np.all(np.isfinite(scheme.compute(p)["chi_bar"]))
