"""Affine-scan Doppler averages: pole form, class sum and explicit per-class solves agree."""
import numpy as np
import pytest

from gabes import affine_doppler, atoms, constants, core, kernels, pole_doppler, schemes
from gabes.affine_doppler import AffineScan


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
    kv, w = affine_doppler.maxwell_shift_classes(330.0, constants.MASS_85RB, constants.K_VEC, True)
    return affine, np.linspace(-3, 3, 41) * constants.GAMMA, kv, w


def test_maxwell_shift_classes_off_is_one_stationary_class():
    kv, w = affine_doppler.maxwell_shift_classes(330.0, constants.MASS_85RB,
                                                 constants.K_VEC, False)
    np.testing.assert_array_equal(kv, [0.0])
    np.testing.assert_array_equal(w, [1.0])
    kv, w = affine_doppler.maxwell_shift_classes(330.0, constants.MASS_85RB,
                                                 constants.K_VEC, True)
    assert kv.size == w.size > 100 and np.isclose(w.sum(), 1.0)


def test_grid_mean_is_the_weighted_per_class_steady_state():
    atom, build_H0, probe = _lambda_scan()
    scan = np.linspace(-2, 2, 9) * constants.GAMMA
    kv = np.array([-3.0, 0.0, 1.5]) * constants.GAMMA
    w = np.array([0.2, 0.5, 0.3])
    shift = atom.S_v.copy()
    got = affine_doppler.AffineScan(atom, build_H0, (2, 0), probe, shift).grid_mean(scan, kv, w)
    expected = np.array([
        sum(wj * core.steady_state_batched(
            core.build_liouvillian(build_H0(s), atom) - s * shift + kj * atom.S_v,
            np.zeros(1), atom.S_v, 3)[0, 2, 0] for kj, wj in zip(kv, w)) / probe
        for s in scan])
    np.testing.assert_allclose(got, expected, rtol=1e-9, atol=1e-12)


def test_stationary_matches_the_single_zero_class():
    atom, build_H0, probe = _lambda_scan()
    scan = np.linspace(-1, 1, 7) * constants.GAMMA
    affine = affine_doppler.AffineScan(atom, build_H0, (2, 0), probe, atom.S_v)
    np.testing.assert_allclose(affine.stationary(scan),
                               affine.grid_mean(scan, np.zeros(1), np.ones(1)),
                               rtol=1e-10, atol=1e-13)
    if kernels.available():   # the kernel path is the same call for both
        np.testing.assert_array_equal(affine.stationary(scan),
                                      affine.grid_mean(scan, np.zeros(1), np.ones(1)))


def test_pole_mean_matches_the_class_sum_on_a_maxwell_grid():
    affine, scan, kv, w = _maxwell_scan()
    sigma = float(np.sqrt(np.sum(w * kv ** 2)))
    mean, failed = affine._pole_mean(scan, kv, w, sigma)
    assert not failed.any()
    np.testing.assert_allclose(mean, affine._class_sum(scan, kv, w), rtol=1e-8, atol=0.0)


def test_stationary_class_never_reaches_the_pole_form(monkeypatch):
    affine, scan, _, _ = _maxwell_scan()
    expected = affine._class_sum(scan, np.zeros(1), np.ones(1))

    def refuse(*args):
        raise AssertionError("kv = 0 has no velocity average to reduce")

    monkeypatch.setattr(AffineScan, "_pole_mean", refuse)
    np.testing.assert_array_equal(affine.grid_mean(scan, np.zeros(1), np.ones(1)), expected)


def test_rows_missing_the_guard_are_summed_class_by_class(monkeypatch):
    affine, scan, kv, w = _maxwell_scan()
    expected = affine._class_sum(scan, kv, w)
    pole_mean = AffineScan._pole_mean

    def one_bad_row(self, *args):
        mean, failed = pole_mean(self, *args)
        mean[7], failed[7] = 1e30, True
        return mean, failed

    monkeypatch.setattr(AffineScan, "_pole_mean", one_bad_row)
    got = affine.grid_mean(scan, kv, w)
    assert got[7] == expected[7]
    np.testing.assert_allclose(got, expected, rtol=1e-8, atol=0.0)


def test_guard_rejects_a_row_off_the_direct_solve(monkeypatch):
    affine, scan, kv, w = _maxwell_scan()
    monkeypatch.setattr(affine_doppler, "POLE_GUARD_RTOL", -1.0)
    sigma = float(np.sqrt(np.sum(w * kv ** 2)))
    assert affine._pole_mean(scan, kv, w, sigma)[1].all()
    np.testing.assert_array_equal(affine.grid_mean(scan, kv, w), affine._class_sum(scan, kv, w))


def test_singular_expansion_falls_back_to_the_class_sum(monkeypatch):
    affine, scan, kv, w = _maxwell_scan()

    def singular(*args, **kwargs):
        raise np.linalg.LinAlgError("Singular matrix")

    monkeypatch.setattr(pole_doppler.ShiftInvertSystem, "residues", singular)
    np.testing.assert_array_equal(affine.grid_mean(scan, kv, w), affine._class_sum(scan, kv, w))


@pytest.mark.parametrize("name, params", [
    ("lambda", {"view": "EIT"}),
    ("lambda", {"view": "CPT", "species": "⁸⁷Rb", "line": "D2"}),
    ("rydberg_eit", {"view": "EIT"}),
    ("rydberg_eit", {"view": "AT electrometry"}),
])
def test_doppler_on_scheme_matches_the_velocity_grid(monkeypatch, name, params):
    scheme = schemes.get(name)
    p = scheme.defaults()
    p.update(params, doppler="on")
    calls = []
    class_sum = AffineScan._class_sum

    def spy(self, scan, kv, w):
        calls.append(kv.size)
        return class_sum(self, scan, kv, w)

    monkeypatch.setattr(AffineScan, "_class_sum", spy)
    got = scheme.compute(p)["chi_bar"]
    assert not calls
    monkeypatch.setattr(AffineScan, "grid_mean", class_sum)
    np.testing.assert_allclose(got, scheme.compute(p)["chi_bar"], rtol=1e-6, atol=0.0)
