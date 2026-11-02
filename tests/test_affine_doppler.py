"""Affine-scan Doppler averages: kernel, fallback and explicit per-class solves agree."""
import numpy as np

from gabes import affine_doppler, atoms, constants, core, kernels


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
