"""
Exact velocity-class reduction by poles and residues.

Many GABES steady states depend on atomic velocity only through an affine
optical detuning shift.  After the trace row is inserted, one row of a scan is

    (A00 + t*A01 - D*S) x = rhs,        D = Delta - k*v,

and the observables are linear readouts ``r_i . x``.  Split the coordinates into
the range of ``S`` (optical coherences, subscript R) and its null space (N).
``D`` appears only in the R-R block, so Gaussian elimination of N is exact and
independent of ``D``:

    K = A_RR - A_RN A_NN^-1 A_NR,      Z = S_RR^-1 K,
    r_i(D) = c0_i - sum_k res_ik / (lam_k - D),      lam = eig(Z).

Every row is therefore an exactly known rational function of ``D``.  A velocity
average costs one eigendecomposition per row instead of one linear solve per
velocity class:

- :meth:`PoleResidues.discrete_mean` evaluates any velocity measure, e.g. the
  truncated uniform Maxwell grid that a brute-force tier uses, exactly;
- :meth:`PoleResidues.gaussian_mean` evaluates the untruncated Maxwell average
  in closed form with the Faddeeva function.

The eigenproblem is delegated to LAPACK through NumPy on a small thread pool;
the Schur complements, residues and pole sums run in compiled kernels when Numba
is available and in NumPy otherwise.  :func:`residues_many` pools the
eigenproblems of several systems, so small scan batches still run in parallel.
:func:`floquet_real_form` builds the real system for a finite-Floquet steady
state whose couplings are adjoint partners.
"""
import math
import os
from concurrent.futures import ThreadPoolExecutor

import numpy as np

from . import doppler, kernels

_POOL = None
_POOL_WORKERS = max(1, int(os.cpu_count() or 1))
_MIN_EIG_CHUNK = 4          # below this a task costs more to dispatch than to run


def _pool():
    global _POOL
    if _POOL is None:
        _POOL = ThreadPoolExecutor(max_workers=_POOL_WORKERS,
                                   thread_name_prefix="gabes-poles")
    return _POOL


def _eig_stacks(stacks):
    """Eigen-decompose several real (rows, n, n) stacks on one pooled task list."""
    total = sum(stack.shape[0] for stack in stacks)
    chunk = max(_MIN_EIG_CHUNK, math.ceil(total / _POOL_WORKERS))
    tasks = [(k, i) for k, stack in enumerate(stacks)
             for i in range(0, stack.shape[0], chunk)]

    def one(task):
        k, i = task
        lam, V = np.linalg.eig(stacks[k][i:i + chunk])
        return lam.astype(np.complex128, copy=False), V.astype(np.complex128, copy=False)

    if len(tasks) <= 1 or _POOL_WORKERS == 1:
        parts = [one(task) for task in tasks]
    else:
        parts = list(_pool().map(one, tasks))
    out, cursor = [], 0
    for stack in stacks:
        count = math.ceil(stack.shape[0] / chunk)
        pieces = parts[cursor:cursor + count]
        cursor += count
        out.append((np.concatenate([p[0] for p in pieces]),
                    np.concatenate([p[1] for p in pieces])))
    return out


class PoleResidues:
    """Rows of ``r(D) = c0 - sum_k res[..., k] / (lam[k] - D)``.

    ``lam`` has shape (rows, n_poles), ``res`` (rows, n_readouts, n_poles) and
    ``c0`` (rows, n_readouts).
    """

    def __init__(self, lam, res, c0):
        self.lam = np.ascontiguousarray(lam, dtype=np.complex128)
        self.res = np.ascontiguousarray(res, dtype=np.complex128)
        self.c0 = np.ascontiguousarray(c0, dtype=np.complex128)

    @property
    def min_abs_imag(self):
        """Smallest distance of a pole from the real detuning axis [rad/s]."""
        return float(np.min(np.abs(self.lam.imag))) if self.lam.size else math.inf

    def discrete_mean(self, nodes, weights):
        """Exact ``sum_j weights_j * r(nodes_j)`` for every row (no interpolation)."""
        nodes = np.ascontiguousarray(nodes, dtype=np.float64)
        weights = np.ascontiguousarray(weights, dtype=np.float64)
        if kernels.available():
            return kernels.pole_discrete_mean(self.lam, self.res, self.c0, nodes, weights)
        out = np.empty(self.c0.shape, dtype=np.complex128)
        wsum = float(np.sum(weights))
        step = max(1, int(2_000_000 // max(1, self.lam.shape[1] * nodes.size)))
        for i in range(0, self.lam.shape[0], step):
            lam = self.lam[i:i + step, :, None]
            s = np.sum(weights[None, None, :] / (lam - nodes[None, None, :]), axis=2)
            out[i:i + step] = (self.c0[i:i + step] * wsum
                               - np.einsum("rok,rk->ro", self.res[i:i + step], s))
        return out

    def value_at(self, D):
        """Readouts at one real detuning ``D`` for every row."""
        return self.discrete_mean(np.array([float(D)]), np.array([1.0]))

    def gaussian_mean(self, center, sigma):
        """Untruncated average over ``D ~ Normal(center, sigma**2)``."""
        mean = doppler.maxwell_resolvent_mean(self.lam, center, sigma)
        return self.c0 - np.einsum("rok,rk->ro", self.res, mean)


class AffineShiftSystem:
    """Pole–residue factory for ``(A00 + t*A01 - D*S) x = rhs`` with readouts ``W x``.

    ``A00``, ``A01``, ``S`` and ``rhs`` must be real; ``readouts`` may be complex.
    ``S`` restricted to its row/column support must be invertible.
    """

    def __init__(self, A00, A01, S, rhs, readouts):
        A00 = np.asarray(A00)
        A01 = np.asarray(A01)
        S = np.asarray(S)
        rhs = np.asarray(rhs)
        for name, arr in (("A00", A00), ("A01", A01), ("S", S), ("rhs", rhs)):
            if np.iscomplexobj(arr):
                raise TypeError(f"{name} must be real")
        readouts = np.atleast_2d(np.asarray(readouts, dtype=np.complex128))
        n = A00.shape[0]
        if A00.shape != (n, n) or A01.shape != (n, n) or S.shape != (n, n):
            raise ValueError("A00, A01 and S must be square with equal size")
        if rhs.shape != (n,) or readouts.shape[1] != n:
            raise ValueError("rhs/readout length must match the system size")
        support = np.any(S != 0.0, axis=0) | np.any(S != 0.0, axis=1)
        R = np.flatnonzero(support)
        N = np.flatnonzero(~support)
        if R.size == 0:
            raise ValueError("S has no support; the response does not depend on D")
        S_RR = S[np.ix_(R, R)]
        if np.linalg.cond(S_RR) > 1e12:
            raise ValueError("S restricted to its support must be invertible")
        take = lambda A, a, b: np.ascontiguousarray(A[np.ix_(a, b)], dtype=np.float64)
        self.size = n
        self.R, self.N = R, N
        self.NN0, self.NN1 = take(A00, N, N), take(A01, N, N)
        self.NR0, self.NR1 = take(A00, N, R), take(A01, N, R)
        self.RN0, self.RN1 = take(A00, R, N), take(A01, R, N)
        self.RR0, self.RR1 = take(A00, R, R), take(A01, R, R)
        self.S_RR_inv = np.ascontiguousarray(np.linalg.inv(S_RR))
        self.eN = np.ascontiguousarray(rhs[N], dtype=np.float64)
        self.eR = np.ascontiguousarray(rhs[R], dtype=np.float64)
        self.WR = np.ascontiguousarray(readouts[:, R])
        self.WN = np.ascontiguousarray(readouts[:, N])

    def _schur(self, t):
        """Stage 1: ``Z``, ``g``, ``c0`` and ``weff`` for every scan parameter."""
        if kernels.available():
            return kernels.pole_schur_rows(
                self.NN0, self.NN1, self.NR0, self.NR1, self.RN0, self.RN1,
                self.RR0, self.RR1, self.S_RR_inv, self.eN, self.eR,
                self.WR, self.WN, t)
        tt = t[:, None, None]
        NN = self.NN0[None] + tt * self.NN1[None]
        NR = self.NR0[None] + tt * self.NR1[None]
        RN = self.RN0[None] + tt * self.RN1[None]
        RR = self.RR0[None] + tt * self.RR1[None]
        rhs = np.concatenate(
            [NR, np.broadcast_to(self.eN[None, :, None], (t.size, self.eN.size, 1))],
            axis=2)
        sol = np.linalg.solve(NN, rhs)
        F, y = sol[:, :, :-1], sol[:, :, -1]
        Z = self.S_RR_inv[None] @ (RR - RN @ F)
        u = (RN @ y[:, :, None])[:, :, 0] - self.eR[None]
        g = u @ self.S_RR_inv.T
        c0 = y @ self.WN.T
        weff = self.WR[None] - np.einsum("rn,bnm->brm", self.WN, F)
        return Z, g, c0, weff

    def _residues(self, V, g, weff):
        """Stage 3: ``res = (weff V) * (V^-1 g)``."""
        if kernels.available():
            return kernels.pole_residue_rows(V, np.ascontiguousarray(g), weff)
        b = np.linalg.solve(V, g[:, :, None].astype(np.complex128))[:, :, 0]
        return (weff @ V) * b[:, None, :]

    def _empty(self):
        nR, n_ro = self.RR0.shape[0], self.WR.shape[0]
        return PoleResidues(np.empty((0, nR)), np.empty((0, n_ro, nR)),
                            np.empty((0, n_ro)))

    def residues(self, t):
        """Pole–residue rows for the scan parameters ``t``."""
        return residues_many([(self, t)])[0]


def residues_many(requests):
    """Pole–residue rows for several ``(AffineShiftSystem, t)`` requests.

    Stages 1 and 3 run per system in compiled kernels; the eigenproblems of all
    requests are pooled by matrix size so a small refinement batch still spreads
    over the worker threads.
    """
    staged = []
    for system, t in requests:
        t = np.ascontiguousarray(np.atleast_1d(t), dtype=np.float64)
        staged.append((system, t, system._schur(t) if t.size else None))
    groups = {}
    for index, (_, _, stage) in enumerate(staged):
        if stage is not None:
            groups.setdefault(stage[0].shape[1], []).append(index)
    sizes = list(groups)
    stacks = [np.ascontiguousarray(np.concatenate([staged[i][2][0] for i in groups[n]]))
              for n in sizes]
    eig_of = {}
    for n, (lam, V) in zip(sizes, _eig_stacks(stacks) if stacks else ()):
        cursor = 0
        for i in groups[n]:
            rows = staged[i][1].size
            eig_of[i] = (lam[cursor:cursor + rows], V[cursor:cursor + rows])
            cursor += rows
    out = []
    for index, (system, _, stage) in enumerate(staged):
        if stage is None:
            out.append(system._empty())
            continue
        _, g, c0, weff = stage
        lam, V = eig_of[index]
        out.append(PoleResidues(lam, system._residues(np.ascontiguousarray(V), g, weff), c0))
    return out


def floquet_real_form(L, L_t, Cp, S, n_f, omega, omega_t, n_levels, readouts):
    """Real finite-Floquet system for :class:`AffineShiftSystem`.

    Harmonic equations in Hermitian operator coordinates (``h = -n_f..n_f``)::

        (L + t*L_t + i*h*(omega + t*omega_t)) x_h + Cp x_(h-1) + conj(Cp) x_(h+1) = 0

    with row 0 of the zero harmonic replaced by the trace (sum of the first
    ``n_levels`` population coordinates) equal to one.  For real ``L``, ``L_t``
    and ``S`` and adjoint-partner couplings, complex conjugation combined with
    harmonic reversal maps solutions to solutions, so ``x_(-h) = conj(x_h)``.
    Writing ``x_h = a_h + i*b_h`` gives an equivalent real system of the same
    size with unknowns ``(a_0, a_1, b_1, ..., a_n_f, b_n_f)``.

    ``readouts`` is a sequence of ``(harmonic, complex weight vector)`` with
    ``harmonic >= 0``; the readout is ``weight . x_harmonic``.
    Returns ``(A00, A01, S_ext, rhs, W)``.
    """
    L = np.asarray(L, dtype=np.float64)
    L_t = np.asarray(L_t, dtype=np.float64)
    S = np.asarray(S, dtype=np.float64)
    Cp = np.asarray(Cp, dtype=np.complex128)
    n_f = int(n_f)
    if n_f < 1:
        raise ValueError("n_f must be at least 1")
    M = L.shape[0]
    D = M * (2 * n_f + 1)
    eye = np.eye(M)
    P, Q = Cp.real, Cp.imag
    A00 = np.zeros((D, D))
    A01 = np.zeros((D, D))

    def a(h):
        return 0 if h == 0 else (2 * h - 1) * M

    def b(h):
        return 2 * h * M

    def put(A, row, col, X):
        A[row:row + M, col:col + M] += X

    put(A00, 0, a(0), L)
    put(A01, 0, a(0), L_t)
    put(A00, 0, a(1), 2.0 * P)
    put(A00, 0, b(1), 2.0 * Q)
    for h in range(1, n_f + 1):
        re, im = a(h), b(h)
        put(A00, re, a(h), L)
        put(A01, re, a(h), L_t)
        put(A00, re, b(h), -h * omega * eye)
        put(A01, re, b(h), -h * omega_t * eye)
        put(A00, im, b(h), L)
        put(A01, im, b(h), L_t)
        put(A00, im, a(h), h * omega * eye)
        put(A01, im, a(h), h * omega_t * eye)
        if h == 1:
            put(A00, re, a(0), P)
            put(A00, im, a(0), Q)
        else:
            put(A00, re, a(h - 1), P)
            put(A00, re, b(h - 1), -Q)
            put(A00, im, a(h - 1), Q)
            put(A00, im, b(h - 1), P)
        if h < n_f:
            put(A00, re, a(h + 1), P)
            put(A00, re, b(h + 1), Q)
            put(A00, im, a(h + 1), -Q)
            put(A00, im, b(h + 1), P)
    A00[0, :] = 0.0
    A01[0, :] = 0.0
    A00[0, :n_levels] = 1.0
    S_ext = np.kron(np.eye(2 * n_f + 1), S)
    S_ext[0, :] = 0.0
    rhs = np.zeros(D)
    rhs[0] = 1.0
    W = np.zeros((len(readouts), D), dtype=np.complex128)
    for r, (harmonic, weight) in enumerate(readouts):
        harmonic = int(harmonic)
        if not 0 <= harmonic <= n_f:
            raise ValueError("readout harmonics must lie in 0..n_f")
        weight = np.asarray(weight, dtype=np.complex128)
        W[r, a(harmonic):a(harmonic) + M] = weight
        if harmonic > 0:
            W[r, b(harmonic):b(harmonic) + M] = 1j * weight
    return A00, A01, S_ext, rhs, W
