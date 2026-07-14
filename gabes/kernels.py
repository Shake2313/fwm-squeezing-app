"""
Numba-compiled solve kernels (optional fast path).

The hot FLOPs of every GABES scheme flow through a handful of batched
tiny-matrix linear solves (M = n_levels² ≤ ~36 for the dense-batch schemes).
At that size the per-call LAPACK overhead dominates `np.linalg.solve`, the
NumPy expression tree allocates large (batch, M, M) temporaries, and the GIL
serialises everything around the solve. These kernels fuse the whole
per-grid-point pipeline (assemble → factor → solve → contract) into one
GIL-free parallel loop with hand-rolled complex LU.

Numba is an optional dependency: if it is missing (or GABES_DISABLE_NUMBA is
set) callers fall back to the NumPy path, so results never depend on it being
installed — only the speed does.

Numerical compatibility: the LU uses LAPACK's zgetf2 pivot rule
(max |Re| + |Im|) and the same right-looking update order, so pivot choices
match `np.linalg.solve` and results agree to ~1e-13 relative (well inside the
1e-9 regression tolerance).
"""
import os

import numpy as np

from . import core                                      # hermitian_basis (real frame)

try:                                                    # pragma: no cover
    if os.environ.get("GABES_DISABLE_NUMBA"):
        raise ImportError("disabled via GABES_DISABLE_NUMBA")
    from numba import njit, prange
    NUMBA_AVAILABLE = True
except ImportError:                                     # pragma: no cover
    NUMBA_AVAILABLE = False

    def njit(*args, **kwargs):                          # no-op decorator stub
        def wrap(fn):
            return fn
        if args and callable(args[0]):
            return args[0]
        return wrap

    prange = range


def available():
    """True if the compiled fast path can be used."""
    return NUMBA_AVAILABLE


@njit(cache=True)
def _lu_factor(A, piv):
    """In-place LU with partial pivoting (LAPACK zgetf2 conventions).

    Pivot magnitude is |Re|+|Im| (LAPACK's cabs1) so pivot choices — and hence
    rounding — track `np.linalg.solve` on the same matrix.
    """
    n = A.shape[0]
    for k in range(n):
        p = k
        amax = abs(A[k, k].real) + abs(A[k, k].imag)
        for i in range(k + 1, n):
            a = abs(A[i, k].real) + abs(A[i, k].imag)
            if a > amax:
                amax = a
                p = i
        piv[k] = p
        if p != k:
            for j in range(n):
                tmp = A[k, j]
                A[k, j] = A[p, j]
                A[p, j] = tmp
        inv = 1.0 / A[k, k]
        for i in range(k + 1, n):
            lik = A[i, k] * inv
            A[i, k] = lik
            for j in range(k + 1, n):
                A[i, j] -= lik * A[k, j]


@njit(cache=True)
def _lu_solve(A, piv, B):
    """Solve A X = B in place for B (M, K) using a factored (A, piv)."""
    n = A.shape[0]
    K = B.shape[1]
    for k in range(n):
        p = piv[k]
        if p != k:
            for j in range(K):
                tmp = B[k, j]
                B[k, j] = B[p, j]
                B[p, j] = tmp
    for i in range(1, n):                       # L y = P b (unit lower)
        for k in range(i):
            lik = A[i, k]
            for j in range(K):
                B[i, j] -= lik * B[k, j]
    for i in range(n - 1, -1, -1):              # U x = y
        for k in range(i + 1, n):
            uik = A[i, k]
            for j in range(K):
                B[i, j] -= uik * B[k, j]
        inv = 1.0 / A[i, i]
        for j in range(K):
            B[i, j] *= inv


@njit(cache=True)
def _lu_factor_real(A, piv):
    """Real in-place LU with partial pivoting (LAPACK dgetf2 conventions).

    The real-arithmetic twin of `_lu_factor`, used by the steady-state kernels
    that run in the Hermitian-generator basis (`core.hermitian_basis`), where the
    Liouvillian is real. Pivot magnitude is |A| (LAPACK's idamax) so pivots — and
    rounding — track `np.linalg.solve` on the same real matrix.
    """
    n = A.shape[0]
    for k in range(n):
        p = k
        amax = abs(A[k, k])
        for i in range(k + 1, n):
            a = abs(A[i, k])
            if a > amax:
                amax = a
                p = i
        piv[k] = p
        if p != k:
            for j in range(n):
                tmp = A[k, j]
                A[k, j] = A[p, j]
                A[p, j] = tmp
        inv = 1.0 / A[k, k]
        for i in range(k + 1, n):
            lik = A[i, k] * inv
            A[i, k] = lik
            for j in range(k + 1, n):
                A[i, j] -= lik * A[k, j]


@njit(cache=True)
def _lu_solve_real(A, piv, B):
    """Solve A X = B in place for real B (M, K) using a factored (A, piv)."""
    n = A.shape[0]
    K = B.shape[1]
    for k in range(n):
        p = piv[k]
        if p != k:
            for j in range(K):
                tmp = B[k, j]
                B[k, j] = B[p, j]
                B[p, j] = tmp
    for i in range(1, n):                       # L y = P b (unit lower)
        for k in range(i):
            lik = A[i, k]
            for j in range(K):
                B[i, j] -= lik * B[k, j]
    for i in range(n - 1, -1, -1):              # U x = y
        for k in range(i + 1, n):
            uik = A[i, k]
            for j in range(K):
                B[i, j] -= uik * B[k, j]
        inv = 1.0 / A[i, i]
        for j in range(K):
            B[i, j] *= inv


@njit(cache=True, parallel=True)
def floquet_chi_grid(L0_base, C_delta, S_v, Cp, Cm, delta_axis, deff_axis,
                     omega_hf, branch, w_probe, w_conj, n_levels):
    """
    Fused FWM χ̄ grid: the 3-mode sideband (Floquet) steady state contracted
    with coherence weights, over the full (δ, Δ_eff) grid in one pass.

    Equivalent to the `chi_matrix_table` δ-loop over `core.floquet_solve`,
    using the affine structure of the Hamiltonian:

        L0(δ, Δ_eff) = L0_base + δ·C_delta − Δ_eff·S_v
        Ω_beat(δ)    = omega_hf + branch·δ

    Per grid point: factor A₋ = L0 − iΩ_beat, solve Xm = A₋⁻¹Cm, and get the
    +sideband for free from the ±conjugate symmetry (P = vec transpose swap):

        A₊ = P·A₋*·P,  Cp = P·Cm*·P  ⇒  Xp = A₊⁻¹Cp = P·Xm*·P,
        Cm·A₊⁻¹Cp = P·(Cp·A₋⁻¹Cm)*·P,

    so the second factor + M-column solve + one M×M product are replaced by an
    O(M²) conjugate-transpose-permute. Then form the feedback
    A_eff = L0 − Cp·Xm − P·(Cp·Xm)*·P, replace the trace row, solve for ρ₀, back
    out ρ₊₁ = −Xp·ρ₀, and contract with w_probe / w_conj:

        chi_probe[i, j] = Σ_k w_probe[k]·ρ₀[k]
        chi_conj[i, j]  = Σ_k w_conj[k]·ρ₊₁[k]

    Workspaces are hoisted per δ row (one prange task), so the inner Δ_eff loop
    runs allocation-free.
    """
    n_d = delta_axis.size
    n_de = deff_axis.size
    M = n_levels * n_levels
    chi_probe = np.empty((n_d, n_de), np.complex128)
    chi_conj = np.empty((n_d, n_de), np.complex128)
    perm = np.empty(M, np.int64)                    # vec transpose swap n·i+j→n·j+i
    for a in range(n_levels):
        for b in range(n_levels):
            perm[a * n_levels + b] = b * n_levels + a

    for i in prange(n_d):
        delta = delta_axis[i]
        ob = omega_hf + branch * delta
        Am = np.empty((M, M), np.complex128)
        Aeff = np.empty((M, M), np.complex128)
        Xm = np.empty((M, M), np.complex128)
        Xp = np.empty((M, M), np.complex128)
        Mfb = np.empty((M, M), np.complex128)       # minus-sideband feedback Cp·Xm
        rhs = np.empty((M, 1), np.complex128)
        rho1 = np.empty(M, np.complex128)
        piv = np.empty(M, np.int64)

        for j in range(n_de):
            deff = deff_axis[j]
            for r in range(M):
                for c in range(M):
                    l0 = L0_base[r, c] + delta * C_delta[r, c] - deff * S_v[r, c]
                    Am[r, c] = l0
                    Aeff[r, c] = l0
                    Xm[r, c] = Cm[r, c]
                Am[r, r] -= 1j * ob

            _lu_factor(Am, piv)
            _lu_solve(Am, piv, Xm)                  # Xm = A₋⁻¹ Cm

            for r in range(M):                      # Xp = A₊⁻¹Cp = P·Xm*·P
                pr = perm[r]
                for c in range(M):
                    Xp[r, c] = np.conj(Xm[pr, perm[c]])

            for r in range(M):                      # minus feedback Cp·Xm
                for c in range(M):
                    acc = 0.0 + 0.0j
                    for k in range(M):
                        acc += Cp[r, k] * Xm[k, c]
                    Mfb[r, c] = acc

            for r in range(M):                      # A_eff −= Cp·Xm + P·(Cp·Xm)*·P
                pr = perm[r]
                for c in range(M):
                    Aeff[r, c] -= Mfb[r, c] + np.conj(Mfb[pr, perm[c]])

            for c in range(M):                      # trace-normalisation row 0
                Aeff[0, c] = 0.0
            for s in range(n_levels):
                Aeff[0, s * n_levels + s] = 1.0
            for r in range(M):
                rhs[r, 0] = 0.0
            rhs[0, 0] = 1.0

            _lu_factor(Aeff, piv)
            _lu_solve(Aeff, piv, rhs)               # rhs = ρ₀

            for r in range(M):                      # ρ₊₁ = −Xp ρ₀
                acc = 0.0 + 0.0j
                for k in range(M):
                    acc += Xp[r, k] * rhs[k, 0]
                rho1[r] = -acc

            cp_acc = 0.0 + 0.0j
            cc_acc = 0.0 + 0.0j
            for k in range(M):
                cp_acc += w_probe[k] * rhs[k, 0]
                cc_acc += w_conj[k] * rho1[k]
            chi_probe[i, j] = cp_acc
            chi_conj[i, j] = cc_acc

    return chi_probe, chi_conj


@njit(cache=True, parallel=True)
def _magneto_buffer_grid_real(L0_all, deff, S_v, n_levels):
    """Real-frame single-region magneto (B, v) grid; see `magneto_buffer_grid`.

    Inputs are real (Hermitian-generator basis); the trace row is Σ_{s<n} r_s = 1.
    Returns the real coordinate vectors r, shape (nB, nv, M).
    """
    nB = L0_all.shape[0]
    nv = deff.size
    M = n_levels * n_levels
    out = np.empty((nB, nv, M), np.float64)
    for b in prange(nB):
        A = np.empty((M, M), np.float64)
        rhs = np.empty((M, 1), np.float64)
        piv = np.empty(M, np.int64)
        for j in range(nv):
            d = deff[j]
            for r in range(M):
                for c in range(M):
                    A[r, c] = L0_all[b, r, c] - d * S_v[r, c]
            for c in range(M):
                A[0, c] = 0.0
            for s in range(n_levels):
                A[0, s] = 1.0
            for r in range(M):
                rhs[r, 0] = 0.0
            rhs[0, 0] = 1.0
            _lu_factor_real(A, piv)
            _lu_solve_real(A, piv, rhs)
            for r in range(M):
                out[b, j, r] = rhs[r, 0]
    return out


def magneto_buffer_grid(L0_all, deff, S_v, n_levels):
    """
    Single-region magneto steady state over the (B, v) grid (the buffer-cell
    path). Mirrors `core.steady_state_from_liouvillian` on
    A = L0_all[b] − deff[j]·S_v with the trace row replaced, parallel over B.

    L0_all and S_v are Hermiticity-preserving, so they are changed once into the
    real Hermitian-generator basis (`core.hermitian_basis`) and the (B × v) grid
    is solved in real arithmetic (~2–4× fewer flops); ρ is reconstructed as
    vec(ρ) = U·r. Returns ρ as (nB, nv, n_levels, n_levels).
    """
    nB = L0_all.shape[0]
    nv = int(np.asarray(deff).size)
    U = core.hermitian_basis(n_levels)
    Uh = U.conj().T
    L0_r = np.ascontiguousarray(np.real(Uh @ L0_all @ U))   # batched BLAS matmul
    S_vr = np.ascontiguousarray(np.real(Uh @ S_v @ U))
    r = _magneto_buffer_grid_real(
        L0_r, np.ascontiguousarray(deff, dtype=np.float64), S_vr, n_levels)
    rho_vec = r @ U.T                                       # real r → complex vec(ρ)
    return rho_vec.reshape(nB, nv, n_levels, n_levels)


@njit(cache=True, parallel=True)
def _magneto_two_region_grid_real(L_light0, L_dark0, deff, S_v,
                                  gamma_out, gamma_in, n_levels):
    """Real-frame two-region magneto (B, v) grid; see `magneto_two_region_grid`.

    Inputs are real (each n-level block in the Hermitian-generator basis); the γ
    exchange is ∝ I (basis-invariant) and the trace row is the two-region
    population sum Σ_{s<n} (r_light[s] + r_dark[s]) = 1. Returns the real
    light-block coordinate vectors, shape (nB, nv, M).
    """
    nB = L_light0.shape[0]
    nv = deff.size
    M = n_levels * n_levels
    M2 = 2 * M
    out = np.empty((nB, nv, M), np.float64)
    for b in prange(nB):
        A = np.zeros((M2, M2), np.float64)           # reused across the v loop
        rhs = np.zeros((M2, 1), np.float64)
        piv = np.empty(M2, np.int64)
        for j in range(nv):
            d = deff[j]
            for r in range(M):
                for c in range(M):
                    A[r, c] = L_light0[b, r, c] - d * S_v[r, c]
                    A[r, M + c] = 0.0
                    A[M + r, c] = 0.0
                    A[M + r, M + c] = L_dark0[b, r, c]
                A[r, r] -= gamma_out
                A[r, M + r] += gamma_in
                A[M + r, r] += gamma_out
                A[M + r, M + r] -= gamma_in
            for c in range(M2):                      # global trace row
                A[0, c] = 0.0
            for s in range(n_levels):                # populations = first n coords
                A[0, s] = 1.0
                A[0, M + s] = 1.0
            for r in range(M2):                      # reset RHS (solved in place)
                rhs[r, 0] = 0.0
            rhs[0, 0] = 1.0
            _lu_factor_real(A, piv)
            _lu_solve_real(A, piv, rhs)              # rhs now holds the solution
            for r in range(M):                       # light block only
                out[b, j, r] = rhs[r, 0]
    return out


def magneto_two_region_grid(L_light0, L_dark0, deff, S_v,
                            gamma_out, gamma_in, n_levels):
    """
    Two-region (light/dark) magneto steady state over the (B, v) grid (the
    paraffin-cell path). Builds the 2M×2M exchange system

        [ L_light − Δ_eff·S_v − γ_out·I        γ_in·I            ]
        [        γ_out·I                  L_dark − γ_in·I         ]

    (only the light block carries the optical Δ_eff shift), replaces the global
    trace row with the two-region population sum, solves, and returns the in-beam
    light block ρ_light as (nB, nv, n_levels, n_levels).

    The hot path: at 2M up to ~200 the (nB,nv,2M,2M) solve dominates. Each block
    is Hermiticity-preserving and the γ exchange is ∝ I (basis-invariant), so the
    system is changed once into the real Hermitian-generator basis
    (`core.hermitian_basis`) and solved in real arithmetic (~2–4× fewer flops,
    half the memory); ρ_light is reconstructed as vec(ρ) = U·r_light. Pin BLAS to
    one thread around the call so the prange threads don't oversubscribe.
    """
    nB = L_light0.shape[0]
    nv = int(np.asarray(deff).size)
    U = core.hermitian_basis(n_levels)
    Uh = U.conj().T
    L_light_r = np.ascontiguousarray(np.real(Uh @ L_light0 @ U))   # batched BLAS
    L_dark_r = np.ascontiguousarray(np.real(Uh @ L_dark0 @ U))
    S_vr = np.ascontiguousarray(np.real(Uh @ S_v @ U))
    r = _magneto_two_region_grid_real(
        L_light_r, L_dark_r, np.ascontiguousarray(deff, dtype=np.float64),
        S_vr, float(gamma_out), float(gamma_in), n_levels)
    rho_vec = r @ U.T                                              # real → vec(ρ)
    return rho_vec.reshape(nB, nv, n_levels, n_levels)


@njit(cache=True, parallel=True)
def _affine_scan_chi_real(base, A_coef, B_coef, scan, kv, weights, w_coh, n_levels):
    """Real-arithmetic core of `affine_scan_chi` (see it for the physics).

    `base, A_coef, B_coef` are the REAL Hermitian-generator-basis coefficient
    matrices; `w_coh` (complex, length M) reconstructs the target coherence from
    the real coordinate vector r as ρ[coh] = Σ_a w_coh[a]·r_a. The trace row is
    Σ_{s<n} r_s = 1 (populations are the first n real coordinates).
    """
    ns = scan.size
    nv = kv.size
    M = n_levels * n_levels
    out = np.empty(ns, np.complex128)
    for i in prange(ns):
        s = scan[i]
        A = np.empty((M, M), np.float64)
        rhs = np.empty((M, 1), np.float64)
        piv = np.empty(M, np.int64)
        acc = 0.0 + 0.0j
        for j in range(nv):
            kvj = kv[j]
            for r in range(M):
                for c in range(M):
                    A[r, c] = base[r, c] + s * A_coef[r, c] + kvj * B_coef[r, c]
            for c in range(M):
                A[0, c] = 0.0
            for st in range(n_levels):
                A[0, st] = 1.0                  # trace: populations = first n coords
            for r in range(M):
                rhs[r, 0] = 0.0
            rhs[0, 0] = 1.0
            _lu_factor_real(A, piv)
            _lu_solve_real(A, piv, rhs)
            coh = 0.0 + 0.0j
            for a in range(M):
                coh += w_coh[a] * rhs[a, 0]
            acc += weights[j] * coh
        out[i] = acc
    return out


def affine_scan_chi(base, A_coef, B_coef, scan, kv, weights, coh_idx, n_levels):
    """
    One probe coherence χ̄(scan) = Σ_v weights[v]·ρ[coh] / (caller's Ω), for a
    Hamiltonian that is affine in the scan variable s (Λ EIT/AT/CPT, Rydberg).

    The per-(s, v) Liouvillian is assembled from three constant matrices:

        L(s, kv) = base + s·A_coef + kv·B_coef

    where (computed once by the caller from two `build_liouvillian` evaluations)

        base    = L₀(s=0)                          [Lindblad already folded in]
        A_coef  = dL/ds − S_v                      [scan shift + two-photon term]
        B_coef  = S_v                              [optical Doppler shift, k·v]

    This collapses the scheme's Python `for s in scan` loop (rebuilding a
    Liouvillian and solving a velocity batch each step) into one parallel pass:
    prange over the scan, an allocation-free inner velocity loop, trace-row
    normalisation, hand-rolled LU, and the velocity-weighted contraction of the
    coherence at flat index `coh_idx = e·n_levels + g`.

    The three constant matrices are Hermiticity-preserving, so they are changed
    once into the real Hermitian-generator basis (`core.hermitian_basis`) and the
    whole scan is solved in real arithmetic (~2–4× fewer flops, half the memory);
    ρ[coh] is reconstructed from the row `U[coh_idx]`. Bit-for-bit the same χ̄ as
    the complex form to machine precision. Covers all three cases identically:
      Λ Doppler-on   : kv = k·v (Maxwell grid),  weights = Maxwell weights
      Λ Doppler-off  : kv = [0],                 weights = [1]
      Rydberg        : kv = [0], B_coef = S_v = 0 (no Doppler), weights = [1]
    """
    U = core.hermitian_basis(n_levels)
    Uh = U.conj().T
    base_r = np.ascontiguousarray(np.real(Uh @ base @ U))
    A_coef_r = np.ascontiguousarray(np.real(Uh @ A_coef @ U))
    B_coef_r = np.ascontiguousarray(np.real(Uh @ B_coef @ U))
    w_coh = np.ascontiguousarray(U[coh_idx])
    return _affine_scan_chi_real(
        base_r, A_coef_r, B_coef_r,
        np.ascontiguousarray(scan, dtype=np.float64),
        np.ascontiguousarray(kv, dtype=np.float64),
        np.ascontiguousarray(weights, dtype=np.float64),
        w_coh, n_levels)
