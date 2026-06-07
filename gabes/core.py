"""
Physics-agnostic OBE engine.

Nothing here knows about a specific atom or experiment: every routine takes the
Hamiltonian / super-operators / dimension as arguments. The level scheme lives
in atoms.AtomModel; the experiment lives in schemes/.

Ported from fwm_obe.py with the only change that the Hilbert-space dimension is
derived from array shapes / passed explicitly instead of a module global.
"""
from contextlib import nullcontext
from functools import lru_cache

import numpy as np

try:
    from threadpoolctl import threadpool_limits as _threadpool_limits
except Exception:                                       # pragma: no cover
    _threadpool_limits = None


def blas_single_thread():
    """Context manager limiting BLAS/LAPACK to one thread.

    Every GABES solve is a large batch of *tiny* dense systems (M = n² ≤ 64).
    For matrices this small the per-call LAPACK threading overhead outweighs any
    parallelism, so a single thread is uniformly faster on the whole workload
    (measured ~25-30% across FWM / Λ / magneto). Scope it around the heavy solve
    so unrelated work is unaffected. No-op if threadpoolctl is unavailable.
    """
    if _threadpool_limits is None:
        return nullcontext()
    return _threadpool_limits(limits=1, user_api="blas")


@lru_cache(maxsize=None)
def _eye(n):
    """Cached n×n complex identity. Read-only: callers must not mutate it."""
    e = np.eye(n, dtype=complex)
    e.flags.writeable = False
    return e


def comm_super(H):
    """vec rule (row-major, vec[N·i+j] = ρ_{ij}):  −i [H, ·].

    Dimension is taken from H.shape (was a fixed module global before).
    """
    n = H.shape[0]
    eye = _eye(n)
    return -1j * (np.kron(H, eye) - np.kron(eye, H.T))


def build_liouvillian(H, atom):
    """L₀ = −i[H,·] + (Lindblad dissipator + dephasing) for the given AtomModel.

    Equivalent to the original `L0_of(H0) = comm_super(H0) + L_LINDBLAD`.
    `atom.lindblad` already bundles spontaneous emission + ground dephasing.
    """
    return comm_super(H) + atom.lindblad


def floquet_solve(L0_at_Deff_zero, Cp, Cm, Omega_beat, Delta_eff_axis,
                  S_v, n_levels):
    """
    3-mode sideband (Floquet) steady state, batched over every Δ_eff.

    Generalisation of the original `solve_sidebands_batched`: the velocity-shift
    super-operator `S_v` and `n_levels` are passed in (were 4-level globals).
    The trace-normalisation row is row 0 (= ρ_index(ground₀, ground₀)).

    Returns ρ_0 and ρ_+1, each reshaped to (N_deff, n_levels, n_levels).
    """
    n = Delta_eff_axis.size
    M = n_levels * n_levels
    eye_rho = _eye(M)
    trace_row = 0

    L0_batch = (L0_at_Deff_zero[None, :, :]
                - Delta_eff_axis[:, None, None] * S_v[None, :, :])

    iO = 1j * Omega_beat * eye_rho
    A_minus = L0_batch - iO[None, :, :]
    A_plus = L0_batch + iO[None, :, :]

    Cm_batch = np.broadcast_to(Cm, (n, M, M))
    Cp_batch = np.broadcast_to(Cp, (n, M, M))

    Am_inv_Cm = np.linalg.solve(A_minus, Cm_batch)
    Ap_inv_Cp = np.linalg.solve(A_plus, Cp_batch)

    minus_feedback = Cp_batch @ Am_inv_Cm
    plus_feedback = Cm_batch @ Ap_inv_Cp

    A_eff = L0_batch - minus_feedback - plus_feedback
    A_eff[:, trace_row, :] = 0
    for state in range(n_levels):
        A_eff[:, trace_row, state * n_levels + state] = 1

    rhs = np.zeros((n, M, 1), dtype=complex)
    rhs[:, trace_row, 0] = 1
    rho_0_vec = np.linalg.solve(A_eff, rhs)
    # ρ₊₁ = −A₊⁻¹ (Cp ρ₀) = −(A₊⁻¹ Cp) ρ₀  (associativity): reuse the factor
    # already solved for the +sideband feedback instead of a second linear solve.
    rho_p1_vec = -(Ap_inv_Cp @ rho_0_vec)

    rho_0 = rho_0_vec[:, :, 0].reshape(n, n_levels, n_levels)
    rho_p1 = rho_p1_vec[:, :, 0].reshape(n, n_levels, n_levels)
    return rho_0, rho_p1


def steady_state_batched(L0_at_Deff_zero, Delta_eff_axis, S_v, n_levels):
    """
    Single-mode steady state ρ (Lρ = 0, trace = 1), batched over every Δ_eff.

    Same velocity-batching idea as `floquet_solve` but without the sideband
    coupling — this is the engine for the absorption-cluster schemes (OD / AT /
    EIT / CPT), where a (weak) probe sits inside H₀ and we want the steady-state
    coherence. Trace-normalisation replaces row 0 (= ρ_index(level₀, level₀)).

    Returns ρ reshaped to (N_deff, n_levels, n_levels).
    """
    n = Delta_eff_axis.size
    M = n_levels * n_levels
    A = (L0_at_Deff_zero[None, :, :]
         - Delta_eff_axis[:, None, None] * S_v[None, :, :])
    A[:, 0, :] = 0
    for state in range(n_levels):
        A[:, 0, state * n_levels + state] = 1
    rhs = np.zeros((n, M, 1), dtype=complex)
    rhs[:, 0, 0] = 1
    rho_vec = np.linalg.solve(A, rhs)
    return rho_vec[:, :, 0].reshape(n, n_levels, n_levels)


def steady_state_from_liouvillian(L_batch, n_levels, trace_row=0):
    """Steady state ρ (Lρ = 0, trace = 1) for a stack of Liouvillians.

    `L_batch` has shape (..., M, M) with arbitrary leading batch dims (e.g.
    scan × velocity); `np.linalg.solve` batches over all of them. This is the
    generic engine the schemes use to collapse outer Python scan/B-field loops
    into a single batched solve (each Liouvillian already carries whatever
    Hamiltonian / velocity shift the caller folded in). `L_batch` is copied, not
    mutated. Returns ρ reshaped to (..., n_levels, n_levels).
    """
    A = np.array(L_batch, dtype=complex)              # own copy (trace row edited)
    M = n_levels * n_levels
    A[..., trace_row, :] = 0
    for state in range(n_levels):
        A[..., trace_row, state * n_levels + state] = 1
    rhs = np.zeros(A.shape[:-1] + (1,), dtype=complex)
    rhs[..., trace_row, 0] = 1
    rho_vec = np.linalg.solve(A, rhs)
    return rho_vec[..., 0].reshape(A.shape[:-2] + (n_levels, n_levels))


def matrix_exp_2x2(M, L):
    """Closed-form exp(M·L) for a batched stack of complex 2×2 matrices."""
    s = 0.5 * (M[..., 0, 0] + M[..., 1, 1])
    q00 = M[..., 0, 0] - s
    q01 = M[..., 0, 1]
    q10 = M[..., 1, 0]
    q11 = M[..., 1, 1] - s
    c2 = q00 * q11 - q01 * q10           # = det(q),  with trace(q)=0
    c = np.sqrt(-c2 + 0j)                # exp(q L) = cosh(cL)I + sinh(cL)/c · q
    big = np.abs(c) > 1e-30
    safe_c = np.where(big, c, 1.0)
    sinh_over_c = np.where(big, np.sinh(c * L) / safe_c, L * np.ones_like(c))
    cosh_cL = np.cosh(c * L)
    exp_sL = np.exp(s * L)

    out = np.empty_like(M)
    out[..., 0, 0] = exp_sL * (cosh_cL + sinh_over_c * q00)
    out[..., 0, 1] = exp_sL * (sinh_over_c * q01)
    out[..., 1, 0] = exp_sL * (sinh_over_c * q10)
    out[..., 1, 1] = exp_sL * (cosh_cL + sinh_over_c * q11)
    return out
