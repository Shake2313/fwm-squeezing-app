"""The one trace-bordered steady-state solve in gabes.core and its callers."""
import numpy as np
import pytest

from gabes import core


def _lindblad(n, seed):
    rng = np.random.default_rng(seed)
    a = rng.standard_normal((n, n)) + 1j * rng.standard_normal((n, n))
    L = core.comm_super(0.5 * (a + a.conj().T))
    for _ in range(3):
        c = rng.standard_normal((n, n)) + 1j * rng.standard_normal((n, n))
        L = L + core.dissipator(0.6 * c)
    return L


def test_population_trace_rows():
    assert np.array_equal(core.population_trace(3), [1, 0, 0, 0, 1, 0, 0, 0, 1])
    assert np.array_equal(core.population_trace(3, "hermitian"), [1, 1, 1, 0, 0, 0, 0, 0, 0])
    assert not core.population_trace(3).flags.writeable
    with pytest.raises(ValueError):
        core.population_trace(3, "spherical")


@pytest.mark.parametrize("n", [2, 4, 6])
def test_stationary_solve_is_trace_one_null_vector(n):
    L = np.stack([_lindblad(n, seed) for seed in range(3)])
    kept = L.copy()
    x = core.stationary_solve(L, core.population_trace(n))
    assert np.array_equal(L, kept)                     # input untouched by default
    assert np.allclose(x @ core.population_trace(n), 1.0, atol=1e-13)
    residual = np.einsum("bij,bj->bi", L, x)
    assert np.max(np.abs(residual)) < 1e-12 * np.max(np.abs(L))
    rho = core.steady_state_from_liouvillian(L, n)
    assert np.array_equal(rho.reshape(3, n * n), x)


@pytest.mark.parametrize("n", [2, 3, 5])
def test_real_hermitian_frame_matches_complex_vec_frame(n):
    L0 = _lindblad(n, 7)
    S_v = core.comm_super(np.diag(np.linspace(-1.0, 1.0, n)).astype(complex))
    deff = np.linspace(-2.0, 2.0, 5)
    real = core.steady_state_batched(L0, deff, S_v, n)
    vec = core.steady_state_from_liouvillian(L0[None] - deff[:, None, None] * S_v[None], n)
    assert np.max(np.abs(real - vec)) < 1e-12


def test_coupled_regions_share_one_global_trace():
    n, M = 3, 9
    light, dark = _lindblad(n, 1), _lindblad(n, 2)
    g_out, g_in = 0.4, 0.15
    eye = np.eye(M)
    A = np.block([[light - g_out * eye, g_in * eye], [g_out * eye, dark - g_in * eye]])
    trace = np.concatenate((core.population_trace(n),) * 2)
    x = core.stationary_solve(A, trace)
    occupation = x[:M] @ core.population_trace(n)
    assert abs(x @ trace - 1.0) < 1e-13
    assert abs(occupation - g_in / (g_in + g_out)) < 1e-12   # exchange balance


@pytest.mark.parametrize("dtype", [complex, float])
def test_factored_states_and_trace_zero_response(dtype):
    n = 3
    L = np.stack([_lindblad(n, seed) for seed in range(4)])
    if dtype is float:                                 # real Hermitian frame
        L = core.to_real_liouvillian(L, n)
        trace = core.population_trace(n, "hermitian")
    else:
        trace = core.population_trace(n)
    states, factors = core.stationary_factored(L, trace)
    assert np.max(np.abs(states - core.stationary_solve(L, trace))) < 1e-12
    forcing = np.random.default_rng(3).standard_normal(states.shape).astype(states.dtype)
    response = core.stationary_response(factors, forcing)
    assert np.max(np.abs(response @ trace)) < 1e-12     # row 0 carries the trace
    residual = np.einsum("bij,bj->bi", L, response) + forcing
    assert np.max(np.abs(residual[:, 1:])) < 1e-11 * np.max(np.abs(L))


def test_floquet_solve_is_truncated_order_one():
    n = 4
    L0 = _lindblad(n, 11)
    V = np.tril(np.random.default_rng(5).standard_normal((n, n)), -1).astype(complex)
    Cp, Cm = core.comm_super(V), core.comm_super(V.conj().T)
    S_v = core.comm_super(np.diag([0.0, 0.0, 1.0, 1.0]).astype(complex))
    deff = np.linspace(-1.0, 1.0, 4)
    old = core.floquet_solve(L0, Cp, Cm, 0.8, deff, S_v, n)
    new = core.floquet_solve_truncated(L0, Cp, Cm, 0.8, deff, S_v, n, n_f=1)
    for a, b in zip(old, new):
        assert np.array_equal(a, b)
