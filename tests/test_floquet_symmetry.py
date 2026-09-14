"""Exact finite-cutoff symmetry: independent dense oracle and guarded fallback."""

import numpy as np
import pytest

from gabes import core, kernels


@pytest.fixture(autouse=True)
def single_blas_thread():
    with core.blas_single_thread():
        yield


def physical_grid(n=3, n_f=3, branch=-1):
    rng = np.random.default_rng(42+n)
    z = rng.normal(size=(n, n))+1j*rng.normal(size=(n, n))
    h0 = (z+z.conj().T)/2
    drive = .4*(rng.normal(size=(n, n))+1j*rng.normal(size=(n, n)))
    trace = np.eye(n).reshape(-1)
    populations = np.arange(1., n+1)
    target = np.diag(populations/populations.sum()).reshape(-1)
    reset = .7*(np.outer(target, trace)-np.eye(n*n))
    return [core.comm_super(h0)+reset, core.comm_super(np.diag(np.arange(n))),
        core.comm_super(np.diag(np.linspace(0., .9, n))),
        core.comm_super(drive), core.comm_super(drive.conj().T),
        np.array([-.8, .2, 1.7]), np.array([-3.1, .7]), 9.3, branch,
        rng.normal(size=n*n)+1j*rng.normal(size=n*n),
        rng.normal(size=n*n)+1j*rng.normal(size=n*n), n, n_f]


@pytest.mark.parametrize('n,n_f,branch', [(2, 1, -1), (2, 4, 1), (3, 2, -1),
    (3, 3, 1), (4, 1, 1), (4, 3, -1), (4, 5, 1)])
def test_symmetry_kernel_matches_full_density_block_at_finite_cutoff(n, n_f, branch):
    args = physical_grid(n, n_f, branch)
    actual = kernels.floquet_chi_grid(*args)
    l0, cd, sv, cp, cm, deltas, deff, whf, br, wp, wc, _, _ = args
    # A conjugate without the raw vec transpose permutation is NOT the identity.
    assert np.max(abs(cm-cp.conj())) > .1
    for i, delta in enumerate(deltas):
        modes = core.floquet_solve_direct(l0+delta*cd, cp, cm, whf+br*delta,
            deff, sv, n, n_f, return_harmonics=True)
        expected = (modes[:, n_f].reshape(len(deff), -1)@wp,
                    modes[:, n_f+1].reshape(len(deff), -1)@wc)
        for value, reference in zip(actual, expected):
            np.testing.assert_allclose(value[i], reference, rtol=1e-11, atol=1e-12)
        diagnostics = core.floquet_solution_diagnostics(l0+delta*cd, cp, cm,
            whf+br*delta, deff, sv, modes)
        assert diagnostics['max_normalized_residual'] < 1e-11


@pytest.mark.parametrize('changed', [0, 1, 2, 4, 7])
def test_broken_symmetry_or_complex_frequency_uses_general_solver(monkeypatch, changed):
    args = physical_grid()
    if changed == 7:
        args[changed] += .03j
    else:
        # Even a one-ULP defect must not be called mathematically zero.
        old = args[changed][0, 1]
        args[changed][0, 1] = np.nextafter(old.real, np.inf)+1j*old.imag
    expected = kernels._floquet_chi_grid_general(*args)
    def forbidden(*args):
        raise AssertionError('non-symmetric input reached the exact reduction')
    monkeypatch.setattr(kernels, '_floquet_chi_grid_hermitian', forbidden)
    actual = kernels.floquet_chi_grid(*args)
    for a, b in zip(actual, expected):
        np.testing.assert_array_equal(a, b)


@pytest.mark.skipif(not kernels.available(), reason='Numba call-count instrumentation')
def test_exact_reduction_requires_only_nf_complex_and_one_real_factorization(monkeypatch):
    args = physical_grid(n=2, n_f=4)
    expected = kernels._floquet_chi_grid_general(*args)
    counts = {'complex': 0, 'real': 0}
    complex_factor, real_factor = kernels._lu_factor, kernels._lu_factor_real
    def count_complex(a, p):
        counts['complex'] += 1
        return complex_factor(a, p)
    def count_real(a, p):
        counts['real'] += 1
        return real_factor(a, p)
    monkeypatch.setattr(kernels, '_lu_factor', count_complex)
    monkeypatch.setattr(kernels, '_lu_factor_real', count_real)
    monkeypatch.setattr(kernels, '_floquet_chi_grid_hermitian', kernels._floquet_chi_grid_hermitian.py_func)
    actual = kernels.floquet_chi_grid(*args)
    points = len(args[5])*len(args[6])
    assert counts == {'complex': points*args[-1], 'real': points}
    for a, b in zip(actual, expected):
        np.testing.assert_allclose(a, b, rtol=1e-11, atol=1e-12)


def test_wrapper_preserves_readonly_inputs_and_invalid_order_error():
    args = physical_grid()
    saved = []
    for value in args:
        if isinstance(value, np.ndarray):
            saved.append(value.copy())
            value.setflags(write=False)
    kernels.floquet_chi_grid(*args)
    for before, after in zip(saved, (a for a in args if isinstance(a, np.ndarray))):
        np.testing.assert_array_equal(before, after)
    args[-1] = 0
    with pytest.raises(ValueError, match='n_f'):
        kernels.floquet_chi_grid(*args)
