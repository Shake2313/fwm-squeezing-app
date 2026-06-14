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


@njit(cache=True, parallel=True)
def floquet_chi_grid(L0_base, C_delta, S_v, Cp, Cm, delta_axis, deff_axis,
                     omega_hf, branch, w_probe, w_conj, n_levels):
    """
    Fused FWM χ̄ grid: the 3-mode sideband (Floquet) steady state contracted
    with coherence weights, over the full (δ, Δ_eff) grid in one pass.

    Equivalent to the `chi_matrix_table` δ-loop over `core.floquet_solve`,
    using the affine structure of the Hamiltonian:

        L0(δ, Δ_eff) = L0_base + δ·C_delta − Δ_eff·S_v
        Ω_beat(δ)    = omega_hf − branch·δ

    Per grid point: factor A∓ = L0 ∓ iΩ_beat, form the sideband feedback
    A_eff = L0 − Cp·A₋⁻¹Cm − Cm·A₊⁻¹Cp, replace the trace row, solve for ρ₀,
    back out ρ₊₁ = −(A₊⁻¹Cp)ρ₀, and contract with w_probe / w_conj:

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

    for i in prange(n_d):
        delta = delta_axis[i]
        ob = omega_hf - branch * delta
        Am = np.empty((M, M), np.complex128)
        Ap = np.empty((M, M), np.complex128)
        Aeff = np.empty((M, M), np.complex128)
        Xm = np.empty((M, M), np.complex128)
        Xp = np.empty((M, M), np.complex128)
        rhs = np.empty((M, 1), np.complex128)
        rho1 = np.empty(M, np.complex128)
        piv = np.empty(M, np.int64)

        for j in range(n_de):
            deff = deff_axis[j]
            for r in range(M):
                for c in range(M):
                    l0 = L0_base[r, c] + delta * C_delta[r, c] - deff * S_v[r, c]
                    Am[r, c] = l0
                    Ap[r, c] = l0
                    Aeff[r, c] = l0
                    Xm[r, c] = Cm[r, c]
                    Xp[r, c] = Cp[r, c]
                Am[r, r] -= 1j * ob
                Ap[r, r] += 1j * ob

            _lu_factor(Am, piv)
            _lu_solve(Am, piv, Xm)                  # Xm = A₋⁻¹ Cm
            _lu_factor(Ap, piv)
            _lu_solve(Ap, piv, Xp)                  # Xp = A₊⁻¹ Cp

            for r in range(M):                      # A_eff −= Cp·Xm + Cm·Xp
                for c in range(M):
                    acc = 0.0 + 0.0j
                    for k in range(M):
                        acc += Cp[r, k] * Xm[k, c] + Cm[r, k] * Xp[k, c]
                    Aeff[r, c] -= acc

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
def magneto_buffer_grid(L0_all, deff, S_v, n_levels):
    """
    Single-region magneto steady state over the (B, v) grid (the buffer-cell
    path). Mirrors `core.steady_state_from_liouvillian` on
    A = L0_all[b] − deff[j]·S_v with the trace row replaced, but builds each
    M×M system in place and runs the grid in parallel (prange over B).

    Returns ρ as (nB, nv, n_levels, n_levels).
    """
    nB = L0_all.shape[0]
    nv = deff.size
    M = n_levels * n_levels
    out = np.empty((nB, nv, n_levels, n_levels), np.complex128)
    for b in prange(nB):
        A = np.empty((M, M), np.complex128)
        rhs = np.empty((M, 1), np.complex128)
        piv = np.empty(M, np.int64)
        for j in range(nv):
            d = deff[j]
            for r in range(M):
                for c in range(M):
                    A[r, c] = L0_all[b, r, c] - d * S_v[r, c]
            for c in range(M):
                A[0, c] = 0.0
            for s in range(n_levels):
                A[0, s * n_levels + s] = 1.0
            for r in range(M):
                rhs[r, 0] = 0.0
            rhs[0, 0] = 1.0
            _lu_factor(A, piv)
            _lu_solve(A, piv, rhs)
            for r in range(n_levels):
                for c in range(n_levels):
                    out[b, j, r, c] = rhs[r * n_levels + c, 0]
    return out


@njit(cache=True, parallel=True)
def magneto_two_region_grid(L_light0, L_dark0, deff, S_v,
                            gamma_out, gamma_in, n_levels):
    """
    Two-region (light/dark) magneto steady state over the (B, v) grid (the
    paraffin-cell path). Builds the 2M×2M exchange system

        [ L_light − Δ_eff·S_v − γ_out·I        γ_in·I            ]
        [        γ_out·I                  L_dark − γ_in·I         ]

    in place (only the light block carries the optical Δ_eff shift), replaces
    the global trace row 0 with the two-region population sum, solves, and
    returns the in-beam light block ρ_light as (nB, nv, n_levels, n_levels).

    This is the hot path: at 2M up to ~200 the giant (nB,nv,2M,2M) NumPy
    temporary and its serial batched solve dominate; the kernel fuses assembly
    and solve per grid point and parallelises over B. Like the other kernels it
    uses the hand-rolled `_lu_factor`/`_lu_solve` (LAPACK zgetf2 conventions, so
    results match `np.linalg.solve` to ~1e-13) rather than numba's
    `np.linalg.solve`, which needs a SciPy/LAPACK runtime that is not guaranteed
    to be present. Pin BLAS to one thread around the call so the prange threads
    don't oversubscribe.
    """
    nB = L_light0.shape[0]
    nv = deff.size
    M = n_levels * n_levels
    M2 = 2 * M
    out = np.empty((nB, nv, n_levels, n_levels), np.complex128)
    for b in prange(nB):
        A = np.zeros((M2, M2), np.complex128)        # reused across the v loop
        rhs = np.zeros((M2, 1), np.complex128)
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
            for s in range(n_levels):
                idx = s * n_levels + s
                A[0, idx] = 1.0
                A[0, M + idx] = 1.0
            for r in range(M2):                      # reset RHS (solved in place)
                rhs[r, 0] = 0.0
            rhs[0, 0] = 1.0
            _lu_factor(A, piv)
            _lu_solve(A, piv, rhs)                   # rhs now holds the solution
            for r in range(n_levels):                # light block only
                for c in range(n_levels):
                    out[b, j, r, c] = rhs[r * n_levels + c, 0]
    return out


@njit(cache=True, parallel=True)
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

    Covers all three cases with the same code:
      Λ Doppler-on   : kv = k·v (Maxwell grid),  weights = Maxwell weights
      Λ Doppler-off  : kv = [0],                 weights = [1]
      Rydberg        : kv = [0], B_coef = S_v = 0 (no Doppler), weights = [1]
    """
    ns = scan.size
    nv = kv.size
    M = n_levels * n_levels
    out = np.empty(ns, np.complex128)
    for i in prange(ns):
        s = scan[i]
        A = np.empty((M, M), np.complex128)
        rhs = np.empty((M, 1), np.complex128)
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
                A[0, st * n_levels + st] = 1.0
            for r in range(M):
                rhs[r, 0] = 0.0
            rhs[0, 0] = 1.0
            _lu_factor(A, piv)
            _lu_solve(A, piv, rhs)
            acc += weights[j] * rhs[coh_idx, 0]
        out[i] = acc
    return out
