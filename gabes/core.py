"""
Physics-agnostic OBE engine.

Nothing here knows about a specific atom or experiment: every routine takes the
Hamiltonian / super-operators / dimension as arguments. The level scheme lives
in atoms.AtomModel; the experiment lives in schemes/.

Originally factored from ``fwm_obe.py`` and since generalized to arbitrary
Hilbert-space dimension and finite Floquet order, with independent dense-block
references for numerical validation.
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


@lru_cache(maxsize=None)
def _transpose_perm(n):
    """Cached vec transpose-swap permutation: index n·i+j → n·j+i (M = n²).

    P = perm as a gather: (P·X·P)[r,c] = X[perm[r], perm[c]]. It is an involution
    (perm[perm[k]] = k). Used to exploit the FWM ±sideband conjugate symmetry.
    Read-only.
    """
    perm = np.empty(n * n, dtype=np.intp)
    for i in range(n):
        for j in range(n):
            perm[i * n + j] = j * n + i
    perm.flags.writeable = False
    return perm


@lru_cache(maxsize=None)
def hermitian_basis(n):
    """Orthonormal (Hilbert–Schmidt) Hermitian operator basis for an n-level ρ.

    Returns U (M×M complex, M = n²) whose columns are vec(Gₐ) in the row-major
    vec convention (vec[n·i+j] = ρ_ij) used throughout core. The generators are:
      · n real diagonal projectors |a⟩⟨a|          → coefficient rₐ = ρ_aa
      · symmetric  (|i⟩⟨j|+|j⟩⟨i|)/√2   for i<j
      · antisymmetric  i(|i⟩⟨j|−|j⟩⟨i|)/√2  for i<j

    U is unitary (U†U = I). For any Hermitian ρ the coordinate vector r = U†·vec(ρ)
    is REAL, and any Hermiticity-preserving (physical) Liouvillian L becomes the
    REAL matrix Re(U†·L·U) — so a steady-state solve on it runs in real arithmetic
    (~2–4× fewer flops, half the memory) at the same dimension. Because the first
    n coordinates are exactly the populations (rₛ = ρ_ss), the trace condition
    Σ_s ρ_ss = 1 is just Σ_{s<n} r_s = 1. Read-only: callers must not mutate it.
    """
    M = n * n
    U = np.zeros((M, M), dtype=complex)
    col = 0
    for a in range(n):                              # diagonal projectors
        U[a * n + a, col] = 1.0
        col += 1
    inv_sqrt2 = 1.0 / np.sqrt(2.0)
    for i in range(n):
        for j in range(i + 1, n):
            U[i * n + j, col] = inv_sqrt2           # symmetric
            U[j * n + i, col] = inv_sqrt2
            col += 1
            U[i * n + j, col] = 1j * inv_sqrt2      # antisymmetric
            U[j * n + i, col] = -1j * inv_sqrt2
            col += 1
    U.flags.writeable = False
    return U


def to_real_liouvillian(L, n):
    """Change a (…, M, M) complex Liouvillian into the real Hermitian-generator
    frame: Lʳ = Re(U†·L·U), with U = hermitian_basis(n). Exact for any
    Hermiticity-preserving L (the discarded imaginary part is rounding noise)."""
    U = hermitian_basis(n)
    Uh = U.conj().T
    return np.real(Uh @ L @ U)


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
    `atom.lindblad` already bundles spontaneous emission, generic collapse
    operators, and any explicitly declared coherence dephasing.
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

    Cm_batch = np.broadcast_to(Cm, (n, M, M))
    Cp_batch = np.broadcast_to(Cp, (n, M, M))

    Am_inv_Cm = np.linalg.solve(A_minus, Cm_batch)
    minus_feedback = Cp_batch @ Am_inv_Cm

    # ±sideband conjugate symmetry. For the vec transpose-swap permutation P, a
    # Hermiticity-preserving L₀ obeys P·L₀*·P = L₀ and the couplings obey
    # Cp = P·Cm*·P, so A₊ = P·A₋*·P. Hence A₊⁻¹Cp = P·(A₋⁻¹Cm)*·P and the +side
    # feedback is the conjugate-transpose-permute of the − side — the second
    # factorisation + M-column solve (and one M×M product) drop out entirely.
    perm = _transpose_perm(n_levels)
    Ap_inv_Cp = Am_inv_Cm.conj()[:, perm][:, :, perm]
    plus_feedback = minus_feedback.conj()[:, perm][:, :, perm]

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


def floquet_solve_truncated(L0_at_Deff_zero, Cp, Cm, Omega_beat,
                            Delta_eff_axis, S_v, n_levels, n_f=1, *,
                            return_harmonics=False):
    """Finite ``[-n_f, +n_f]`` Floquet steady state by block elimination.

    The harmonic convention is

    ``(L0 + i*n*Omega) rho_n + Cp rho_(n-1) + Cm rho_(n+1) = 0``.

    Positive and negative harmonic chains are eliminated by exact block
    continued fractions at the selected finite boundary.  Only ``M x M``
    systems (``M=n_levels**2``) are solved, rather than one dense
    ``(2*n_f+1)M`` system per detuning.  The zero harmonic has trace one and all
    nonzero harmonics have trace zero as a consequence of the trace-preserving
    block equations; :func:`floquet_solution_diagnostics` verifies both.

    With ``return_harmonics=False`` the compatibility return is ``(rho_0,
    rho_+1)``.  Otherwise a single array with shape
    ``(n_detuning, 2*n_f+1, n_levels, n_levels)`` is returned in ascending
    harmonic order ``-n_f,...,+n_f``.
    """
    n_f = int(n_f)
    if n_f < 1:
        raise ValueError("n_f must be at least 1")
    if int(n_levels) != n_levels or n_levels < 1:
        raise ValueError("n_levels must be a positive integer")
    n_levels = int(n_levels)
    M = n_levels * n_levels
    L0 = np.asarray(L0_at_Deff_zero, dtype=complex)
    Cp = np.asarray(Cp, dtype=complex)
    Cm = np.asarray(Cm, dtype=complex)
    S_v = np.asarray(S_v, dtype=complex)
    expected = (M, M)
    for name, value in (("L0", L0), ("Cp", Cp), ("Cm", Cm), ("S_v", S_v)):
        if value.shape != expected:
            raise ValueError(f"{name} must have shape {expected}, got {value.shape}")

    deff = np.asarray(Delta_eff_axis, dtype=float)
    if deff.ndim != 1:
        raise ValueError("Delta_eff_axis must be one-dimensional")
    n_deff = deff.size
    eye = _eye(M)
    L_batch = L0[None, :, :] - deff[:, None, None] * S_v[None, :, :]
    Cp_batch = np.broadcast_to(Cp, (n_deff, M, M))
    Cm_batch = np.broadcast_to(Cm, (n_deff, M, M))

    # R[h-1] maps rho_(h-1) -> rho_h for h=1,...,n_f.
    R = [None] * n_f
    next_R = None
    for harmonic in range(n_f, 0, -1):
        A = L_batch + 1j * harmonic * float(Omega_beat) * eye[None, :, :]
        if next_R is not None:
            A = A + Cm_batch @ next_R
        current = -np.linalg.solve(A, Cp_batch)
        R[harmonic - 1] = current
        next_R = current

    # Q[h+n_f] maps rho_(h+1) -> rho_h for h=-n_f,...,-1.
    Q = [None] * n_f
    previous_Q = None
    for harmonic in range(-n_f, 0):
        A = L_batch + 1j * harmonic * float(Omega_beat) * eye[None, :, :]
        if previous_Q is not None:
            A = A + Cp_batch @ previous_Q
        current = -np.linalg.solve(A, Cm_batch)
        Q[harmonic + n_f] = current
        previous_Q = current

    A_eff = L_batch + Cp_batch @ Q[-1] + Cm_batch @ R[0]
    A_eff = np.array(A_eff, copy=True)
    A_eff[:, 0, :] = 0.0
    diagonal = np.arange(n_levels) * (n_levels + 1)
    A_eff[:, 0, diagonal] = 1.0
    rhs = np.zeros((n_deff, M, 1), dtype=complex)
    rhs[:, 0, 0] = 1.0
    rho0 = np.linalg.solve(A_eff, rhs)

    if not return_harmonics:
        rho1 = R[0] @ rho0
        shape = (n_deff, n_levels, n_levels)
        return (rho0[:, :, 0].reshape(shape),
                rho1[:, :, 0].reshape(shape))

    harmonic_vectors = np.empty(
        (n_deff, 2 * n_f + 1, M, 1), dtype=complex)
    harmonic_vectors[:, n_f] = rho0
    current = rho0
    for harmonic in range(1, n_f + 1):
        current = R[harmonic - 1] @ current
        harmonic_vectors[:, n_f + harmonic] = current
    current = rho0
    for harmonic in range(-1, -n_f - 1, -1):
        current = Q[harmonic + n_f] @ current
        harmonic_vectors[:, n_f + harmonic] = current
    return harmonic_vectors[..., 0].reshape(
        n_deff, 2 * n_f + 1, n_levels, n_levels)


def floquet_solve_direct(L0_at_Deff_zero, Cp, Cm, Omega_beat,
                         Delta_eff_axis, S_v, n_levels, n_f=1, *,
                         return_harmonics=False):
    """Independent dense finite-block Floquet reference.

    This intentionally does not share the continued-fraction assembly.  It is
    slower and intended for parity/golden tests and small held-out grids.
    """
    n_f = int(n_f)
    n_levels = int(n_levels)
    if n_f < 1:
        raise ValueError("n_f must be at least 1")
    M = n_levels * n_levels
    L0 = np.asarray(L0_at_Deff_zero, dtype=complex)
    Cp = np.asarray(Cp, dtype=complex)
    Cm = np.asarray(Cm, dtype=complex)
    S_v = np.asarray(S_v, dtype=complex)
    expected = (M, M)
    if any(value.shape != expected for value in (L0, Cp, Cm, S_v)):
        raise ValueError(f"Floquet blocks must all have shape {expected}")
    deff = np.asarray(Delta_eff_axis, dtype=float)
    if deff.ndim != 1:
        raise ValueError("Delta_eff_axis must be one-dimensional")

    sectors = np.arange(-n_f, n_f + 1)
    n_sector = sectors.size
    total = n_sector * M
    block = np.zeros((deff.size, total, total), dtype=complex)
    L_batch = L0[None, :, :] - deff[:, None, None] * S_v[None, :, :]
    eye = _eye(M)
    for q, harmonic in enumerate(sectors):
        sl = slice(q * M, (q + 1) * M)
        block[:, sl, sl] = (
            L_batch + 1j * harmonic * float(Omega_beat) * eye[None, :, :])
        if q > 0:
            block[:, sl, slice((q - 1) * M, q * M)] = Cp
        if q + 1 < n_sector:
            block[:, sl, slice((q + 1) * M, (q + 2) * M)] = Cm

    rhs = np.zeros((deff.size, total, 1), dtype=complex)
    # A literal unit trace row is many orders of magnitude smaller than the
    # Liouvillian-rate rows and makes the deliberately large dense reference
    # needlessly ill-conditioned.  Row scaling leaves the solution unchanged.
    trace_scale = np.maximum(
        np.max(np.abs(block), axis=(-2, -1)), 1.0)
    diagonal = np.arange(n_levels) * (n_levels + 1)
    for q, harmonic in enumerate(sectors):
        row = q * M
        block[:, row, :] = 0.0
        block[:, row, q * M + diagonal] = trace_scale[:, None]
        if harmonic == 0:
            rhs[:, row, 0] = trace_scale
    solution = np.linalg.solve(block, rhs)[..., 0]
    harmonics = solution.reshape(
        deff.size, n_sector, n_levels, n_levels)
    if return_harmonics:
        return harmonics
    return harmonics[:, n_f], harmonics[:, n_f + 1]


def _solve_block_tridiagonal(diagonal, lower, upper, rhs):
    """Solve a batched block-tridiagonal system by block Thomas elimination.

    ``diagonal/lower/upper`` have shape ``(batch, sectors, M, M)`` and ``rhs``
    has shape ``(batch, sectors, M, K)``.  ``lower[:, 0]`` and
    ``upper[:, -1]`` are ignored.  The routine is intentionally private: the
    public Floquet response API below validates the physical block structure
    and trace constraint before calling it.
    """
    diagonal = np.asarray(diagonal, dtype=complex)
    lower = np.asarray(lower, dtype=complex)
    upper = np.asarray(upper, dtype=complex)
    rhs = np.asarray(rhs, dtype=complex)
    if diagonal.ndim != 4 or diagonal.shape[-1] != diagonal.shape[-2]:
        raise ValueError("diagonal must have shape (batch, sectors, M, M)")
    if lower.shape != diagonal.shape or upper.shape != diagonal.shape:
        raise ValueError("lower and upper must match diagonal")
    if rhs.shape[:3] != diagonal.shape[:3]:
        raise ValueError("rhs must have shape (batch, sectors, M, K)")

    batch, sectors, M, _ = diagonal.shape
    K = rhs.shape[-1]
    forward = np.empty((batch, sectors, M, K), dtype=complex)
    upper_solved = np.empty(
        (batch, max(sectors - 1, 0), M, M), dtype=complex)

    previous_upper = None
    previous_forward = None
    for q in range(sectors):
        block = np.array(diagonal[:, q], copy=True)
        effective_rhs = np.array(rhs[:, q], copy=True)
        if q:
            block -= lower[:, q] @ previous_upper
            effective_rhs -= lower[:, q] @ previous_forward
        if q + 1 < sectors:
            combined = np.concatenate((upper[:, q], effective_rhs), axis=-1)
            solved = np.linalg.solve(block, combined)
            previous_upper = solved[:, :, :M]
            previous_forward = solved[:, :, M:]
            upper_solved[:, q] = previous_upper
        else:
            previous_forward = np.linalg.solve(block, effective_rhs)
        forward[:, q] = previous_forward

    solution = np.empty_like(forward)
    solution[:, -1] = forward[:, -1]
    for q in range(sectors - 2, -1, -1):
        solution[:, q] = (
            forward[:, q] - upper_solved[:, q] @ solution[:, q + 1])
    return solution


def _coerce_floquet_perturbations(dL0, dCp, dCm, M):
    """Validate/stack first-order Floquet perturbation blocks."""
    stacked = []
    n_response = None
    for name, value in (("dL0", dL0), ("dCp", dCp), ("dCm", dCm)):
        array = np.asarray(value, dtype=complex)
        if array.ndim == 2:
            array = array[None, :, :]
        if array.ndim != 3 or array.shape[1:] != (M, M):
            raise ValueError(
                f"{name} must have shape (responses, {M}, {M}) or ({M}, {M})")
        if n_response is None:
            n_response = array.shape[0]
        elif array.shape[0] != n_response:
            raise ValueError("dL0, dCp, and dCm must have the same response count")
        stacked.append(array)
    if n_response is None or n_response < 1:
        raise ValueError("at least one perturbation is required")
    return (*stacked, n_response)


def floquet_linear_response_truncated(
        L0_at_Deff_zero, Cp, Cm, Omega_beat, Delta_eff_axis, S_v, n_levels,
        dL0, dCp, dCm, n_f=1):
    """Pump-only periodic state and exact first-order Floquet response.

    First, the trace-one pump-only state ``rho_pump[n]`` is solved on the finite
    harmonic interval ``[-n_f,+n_f]``.  For every declared perturbation ``a``,
    this routine then solves

    ``A_pump delta_rho_a = -(dA_a) rho_pump``

    on the *same* finite Floquet operator, with ``Tr(delta_rho_a[0])=0``.  The
    perturbation blocks ``dL0[a]``, ``dCp[a]``, and ``dCm[a]`` are derivatives
    per unit complex field amplitude; they need not themselves be Hermiticity
    preserving.  Consequently this is a Wirtinger/complex-amplitude response,
    not a finite real seed divided by its amplitude.

    Returns ``(pump_harmonics, response_harmonics)`` with shapes
    ``(batch, 2*n_f+1, n_levels, n_levels)`` and
    ``(batch, responses, 2*n_f+1, n_levels, n_levels)``.
    """
    n_f = int(n_f)
    n_levels = int(n_levels)
    if n_f < 1:
        raise ValueError("n_f must be at least 1")
    M = n_levels * n_levels
    L0 = np.asarray(L0_at_Deff_zero, dtype=complex)
    Cp = np.asarray(Cp, dtype=complex)
    Cm = np.asarray(Cm, dtype=complex)
    S_v = np.asarray(S_v, dtype=complex)
    expected = (M, M)
    if any(value.shape != expected for value in (L0, Cp, Cm, S_v)):
        raise ValueError(f"Floquet blocks must all have shape {expected}")
    dL0, dCp, dCm, n_response = _coerce_floquet_perturbations(
        dL0, dCp, dCm, M)
    deff = np.asarray(Delta_eff_axis, dtype=float)
    if deff.ndim != 1:
        raise ValueError("Delta_eff_axis must be one-dimensional")

    pump = floquet_solve_truncated(
        L0, Cp, Cm, Omega_beat, deff, S_v, n_levels, n_f,
        return_harmonics=True)
    sectors = 2 * n_f + 1
    harmonics = np.arange(-n_f, n_f + 1)
    eye = _eye(M)
    L_batch = L0[None, :, :] - deff[:, None, None] * S_v[None, :, :]
    diagonal = (
        L_batch[:, None, :, :]
        + 1j * float(Omega_beat) * harmonics[None, :, None, None]
        * eye[None, None, :, :])
    diagonal = np.array(diagonal, copy=True)
    lower = np.broadcast_to(
        Cp, (deff.size, sectors, M, M)).copy()
    upper = np.broadcast_to(
        Cm, (deff.size, sectors, M, M)).copy()
    lower[:, 0] = 0.0
    upper[:, -1] = 0.0

    pump_vec = pump.reshape(deff.size, sectors, M)
    source = np.empty(
        (deff.size, sectors, M, n_response), dtype=complex)
    for response in range(n_response):
        source[..., response] = -np.einsum(
            "mk,bqk->bqm", dL0[response], pump_vec)
        source[:, 1:, :, response] -= np.einsum(
            "mk,bqk->bqm", dCp[response], pump_vec[:, :-1])
        source[:, :-1, :, response] -= np.einsum(
            "mk,bqk->bqm", dCm[response], pump_vec[:, 1:])

    # Replace one redundant zero-harmonic equation by Tr(delta rho_0)=0.
    # The complete extended row includes both neighbouring blocks, so their
    # corresponding rows must be removed as well.
    trace_scale = np.maximum(
        np.max(np.abs(diagonal), axis=(1, 2, 3)), 1.0)
    trace_indices = np.arange(n_levels) * (n_levels + 1)
    diagonal[:, n_f, 0, :] = 0.0
    diagonal[:, n_f, 0, trace_indices] = trace_scale[:, None]
    lower[:, n_f, 0, :] = 0.0
    upper[:, n_f, 0, :] = 0.0
    source[:, n_f, 0, :] = 0.0

    response = _solve_block_tridiagonal(
        diagonal, lower, upper, source)
    response = response.transpose(0, 3, 1, 2).reshape(
        deff.size, n_response, sectors, n_levels, n_levels)
    return pump, response


def floquet_linear_response_direct(
        L0_at_Deff_zero, Cp, Cm, Omega_beat, Delta_eff_axis, S_v, n_levels,
        dL0, dCp, dCm, n_f=1):
    """Independent dense-block reference for pump-only linear response.

    This assembles and solves the full extended Floquet matrix.  It deliberately
    does not call the block-Thomas response solver and is intended for held-out
    parity tests and small diagnostic grids.
    """
    n_f = int(n_f)
    n_levels = int(n_levels)
    if n_f < 1:
        raise ValueError("n_f must be at least 1")
    M = n_levels * n_levels
    L0 = np.asarray(L0_at_Deff_zero, dtype=complex)
    Cp = np.asarray(Cp, dtype=complex)
    Cm = np.asarray(Cm, dtype=complex)
    S_v = np.asarray(S_v, dtype=complex)
    expected = (M, M)
    if any(value.shape != expected for value in (L0, Cp, Cm, S_v)):
        raise ValueError(f"Floquet blocks must all have shape {expected}")
    dL0, dCp, dCm, n_response = _coerce_floquet_perturbations(
        dL0, dCp, dCm, M)
    deff = np.asarray(Delta_eff_axis, dtype=float)
    if deff.ndim != 1:
        raise ValueError("Delta_eff_axis must be one-dimensional")

    pump = floquet_solve_direct(
        L0, Cp, Cm, Omega_beat, deff, S_v, n_levels, n_f,
        return_harmonics=True)
    harmonics = np.arange(-n_f, n_f + 1)
    sectors = harmonics.size
    total = sectors * M
    block = np.zeros((deff.size, total, total), dtype=complex)
    L_batch = L0[None, :, :] - deff[:, None, None] * S_v[None, :, :]
    eye = _eye(M)
    for q, harmonic in enumerate(harmonics):
        sl = slice(q * M, (q + 1) * M)
        block[:, sl, sl] = (
            L_batch + 1j * harmonic * float(Omega_beat) * eye[None, :, :])
        if q:
            block[:, sl, slice((q - 1) * M, q * M)] = Cp
        if q + 1 < sectors:
            block[:, sl, slice((q + 1) * M, (q + 2) * M)] = Cm

    pump_vec = pump.reshape(deff.size, sectors, M)
    source = np.empty(
        (deff.size, sectors, M, n_response), dtype=complex)
    for response in range(n_response):
        source[..., response] = -np.einsum(
            "mk,bqk->bqm", dL0[response], pump_vec)
        source[:, 1:, :, response] -= np.einsum(
            "mk,bqk->bqm", dCp[response], pump_vec[:, :-1])
        source[:, :-1, :, response] -= np.einsum(
            "mk,bqk->bqm", dCm[response], pump_vec[:, 1:])
    rhs = source.reshape(deff.size, total, n_response)

    trace_scale = np.maximum(np.max(np.abs(block), axis=(1, 2)), 1.0)
    trace_indices = np.arange(n_levels) * (n_levels + 1)
    trace_row = n_f * M
    block[:, trace_row, :] = 0.0
    block[:, trace_row, n_f * M + trace_indices] = trace_scale[:, None]
    rhs[:, trace_row, :] = 0.0
    solution = np.linalg.solve(block, rhs)
    response = solution.reshape(
        deff.size, sectors, M, n_response).transpose(0, 3, 1, 2)
    response = response.reshape(
        deff.size, n_response, sectors, n_levels, n_levels)
    return pump, response


def trace_zero_liouvillian_response(
        liouvillian, sources, angular_frequency_rad_s, n_levels, *,
        trace_tolerance=1e-10):
    """Solve a driven Liouvillian resolvent on the trace-zero subspace.

    The convention is

    ``delta_rho = -(L + i*omega*I)^(-1) source``.

    ``liouvillian`` has shape ``(..., M, M)`` and ``sources`` has shape
    ``(..., M, K)`` (or ``(..., M)`` for one source), where ``M=n_levels**2``.
    ``angular_frequency_rad_s`` is either one scalar or an array broadcastable
    to the Liouvillian batch shape; it never occupies a response-column axis.
    At nonzero frequency the ordinary resolvent is nonsingular.  At exactly DC,
    an SVD finds the complete left stationary subspace and an augmented solve
    imposes orthogonality to every conserved mode.  This is a full projected
    Drazin action even when the Liouvillian has more than one stationary state;
    a singular inverse is never taken at DC, and a source outside the range of
    ``L`` fails closed.
    """
    n_levels = int(n_levels)
    if n_levels < 1:
        raise ValueError("n_levels must be positive")
    M = n_levels * n_levels
    L = np.asarray(liouvillian, dtype=complex)
    if L.shape[-2:] != (M, M):
        raise ValueError(f"liouvillian must end in shape {(M, M)}")
    source = np.asarray(sources, dtype=complex)
    squeeze = False
    if source.ndim == L.ndim - 1 and source.shape[-1] == M:
        source = source[..., None]
        squeeze = True
    if source.shape[:-2] != L.shape[:-2] or source.shape[-2] != M:
        raise ValueError(
            "sources must have shape liouvillian_batch + (M, K) or (M,)")
    omega = np.asarray(angular_frequency_rad_s, dtype=float)
    try:
        omega_batch = np.broadcast_to(omega, L.shape[:-2])
    except ValueError as exc:
        raise ValueError(
            "angular_frequency_rad_s must be scalar or broadcastable to the "
            "Liouvillian batch shape") from exc
    if not np.isfinite(omega_batch).all():
        raise ValueError("angular_frequency_rad_s must be finite")
    if not (np.isfinite(L).all() and np.isfinite(source).all()):
        raise ValueError("liouvillian and sources must be finite")

    diagonal = np.arange(n_levels) * (n_levels + 1)
    source_trace = np.sum(source[..., diagonal, :], axis=-2)
    # Normalize every response column independently.  Sharing one scale across
    # K would let a large trace-free source hide a smaller traceful source and
    # make the contract depend on which unrelated RHS columns are bundled.
    source_scale = np.maximum(
        np.max(np.abs(source), axis=-2), np.finfo(float).tiny)
    normalized_trace = np.max(
        np.abs(source_trace) / source_scale, initial=0.0)
    if normalized_trace > float(trace_tolerance):
        raise ValueError(
            "sources must be trace-free for a Liouvillian response; "
            f"normalized trace defect is {normalized_trace:.3e}")

    leading = L.shape[:-2]
    K = source.shape[-1]
    A_flat = np.asarray(L, dtype=complex).reshape((-1, M, M))
    rhs_flat = -np.asarray(source, dtype=complex).reshape((-1, M, K))
    omega_flat = np.asarray(omega_batch, dtype=float).reshape(-1)
    response_flat = np.empty_like(rhs_flat)
    dc_mask = omega_flat == 0.0
    nonzero_indices = np.flatnonzero(~dc_mask)
    if nonzero_indices.size:
        shifted = (
            A_flat[nonzero_indices]
            + 1j * omega_flat[nonzero_indices, None, None] * _eye(M)
        )
        response_flat[nonzero_indices] = np.linalg.solve(
            shifted, rhs_flat[nonzero_indices])

    for index in np.flatnonzero(dc_mask):
        matrix = A_flat[index]
        U, singular, _Vh = np.linalg.svd(matrix, full_matrices=True)
        largest = max(float(singular[0]) if singular.size else 0.0, 1.0)
        null_tolerance = max(matrix.shape) * np.finfo(float).eps * largest
        nullity = int(np.count_nonzero(singular <= null_tolerance))
        if nullity < 1:
            raise np.linalg.LinAlgError(
                "DC Liouvillian has no resolved stationary subspace")
        left_null = U[:, -nullity:]
        augmented = np.concatenate(
            (matrix, left_null.conj().T), axis=0)
        augmented_rhs = np.concatenate(
            (rhs_flat[index], np.zeros((nullity, K), dtype=complex)),
            axis=0)
        solution, *_ = np.linalg.lstsq(
            augmented, augmented_rhs,
            rcond=null_tolerance / largest)
        residual = matrix @ solution - rhs_flat[index]
        residual_scale = (
            np.linalg.norm(matrix @ solution, axis=0)
            + np.linalg.norm(rhs_flat[index], axis=0))
        normalized_residual = np.linalg.norm(residual, axis=0) / np.maximum(
            residual_scale, np.finfo(float).tiny)
        constraint_error = np.linalg.norm(
            left_null.conj().T @ solution, axis=0)
        solution_scale = np.maximum(
            np.linalg.norm(solution, axis=0), np.finfo(float).tiny)
        normalized_constraint = constraint_error / solution_scale
        if (np.max(normalized_residual, initial=0.0) > trace_tolerance
                or np.max(normalized_constraint, initial=0.0)
                > trace_tolerance):
            raise np.linalg.LinAlgError(
                "DC source is incompatible with the Liouvillian range or "
                "stationary-subspace projection")
        response_flat[index] = solution
    response = response_flat.reshape(leading + (M, K))
    return response[..., 0] if squeeze else response


def liouvillian_pole_residue_response(
        liouvillian, sources, readouts, angular_frequency_rad_s, *,
        zero_mode_absolute_tolerance=None, residue_relative_tolerance=1e-10):
    """Evaluate a diagonalizable Liouvillian resolvent as an exact pole sum.

    This is a deliberately small, non-batched validation primitive.  For
    ``L=V diag(lambda) V^-1`` and the convention used by
    :func:`trace_zero_liouvillian_response`, each readout/input element is

    ``sum_j residue[j] / (lambda[j] + i*omega)``.

    The returned pole centers are ``-Im(lambda)`` and half-widths are
    ``-Re(lambda)``.  ``visibility`` is ``abs(residue)/half_width`` for stable
    nonstationary modes, not the bare residue magnitude.  At DC the stationary
    pole is projected out only when its
    residue is numerically zero; a nonzero stationary residue is rejected.
    By default the zero-mode threshold is the matrix-rank scale
    ``M*eps*||L||_2``, so a stiff Liouvillian does not relabel a genuine slow
    pole as stationary.  An explicit absolute threshold may be supplied for a
    fixture with independently known spectral accuracy.
    """
    L = np.asarray(liouvillian, dtype=complex)
    if L.ndim != 2 or L.shape[0] != L.shape[1]:
        raise ValueError("liouvillian must be one square matrix")
    M = L.shape[0]
    source = np.asarray(sources, dtype=complex)
    if source.ndim == 1:
        source = source[:, None]
    if source.ndim != 2:
        raise ValueError("sources must be a vector or a two-dimensional matrix")
    readout = np.asarray(readouts, dtype=complex)
    if readout.ndim == 1:
        readout = readout[None, :]
    if readout.ndim != 2:
        raise ValueError("readouts must be a vector or a two-dimensional matrix")
    if source.shape[0] != M or readout.shape[1] != M:
        raise ValueError("sources/readouts must match the Liouvillian dimension")
    if not (np.isfinite(L).all() and np.isfinite(source).all()
            and np.isfinite(readout).all()):
        raise ValueError("liouvillian, sources, and readouts must be finite")
    omega = float(angular_frequency_rad_s)
    if not np.isfinite(omega):
        raise ValueError("angular_frequency_rad_s must be finite")

    eigenvalues, right = np.linalg.eig(L)
    condition = float(np.linalg.cond(right))
    if not np.isfinite(condition):
        raise np.linalg.LinAlgError("Liouvillian eigenvector matrix is singular")
    left_coordinates = np.linalg.solve(right, source)
    read_coordinates = readout @ right
    residues = -np.einsum(
        "rj,jk->jrk", read_coordinates, left_coordinates)
    denominator = eigenvalues + 1j * omega
    matrix_scale = max(float(np.linalg.norm(L, ord=2)), 1.0)
    if zero_mode_absolute_tolerance is None:
        zero_threshold = max(L.shape) * np.finfo(float).eps * matrix_scale
    else:
        zero_threshold = float(zero_mode_absolute_tolerance)
        if not np.isfinite(zero_threshold) or zero_threshold < 0.0:
            raise ValueError(
                "zero_mode_absolute_tolerance must be finite and nonnegative")
    residue_tolerance = float(residue_relative_tolerance)
    if not np.isfinite(residue_tolerance) or residue_tolerance < 0.0:
        raise ValueError(
            "residue_relative_tolerance must be finite and nonnegative")
    stationary = np.abs(eigenvalues) <= zero_threshold
    stationary_residue_max = float(
        np.max(np.abs(residues[stationary]), initial=0.0))
    zero = np.abs(denominator) <= zero_threshold
    terms = np.zeros_like(residues)
    if np.any(zero):
        # This must be relative to the actual response scale.  A unit floor
        # would make the fail-closed zero-pole check depend on arbitrary source
        # or readout units and could silently discard a small but entirely
        # stationary response.
        residue_scale = max(
            float(np.max(np.abs(residues), initial=0.0)),
            np.finfo(float).tiny,
        )
        if np.max(np.abs(residues[zero]), initial=0.0) > (
                residue_tolerance * residue_scale):
            raise np.linalg.LinAlgError(
                "stationary pole has nonzero residue; trace-zero projection failed")
    terms[~zero] = residues[~zero] / denominator[~zero, None, None]
    response = np.sum(terms, axis=0)
    half_widths = -np.real(eigenvalues)
    visible_mode = (~stationary) & (half_widths > zero_threshold)
    visibility = np.zeros_like(np.abs(residues), dtype=float)
    visibility[visible_mode] = (
        np.abs(residues[visible_mode])
        / half_widths[visible_mode, None, None]
    )
    return {
        "response": response,
        "eigenvalues": eigenvalues,
        "residues": residues,
        "pole_centers_rad_s": -np.imag(eigenvalues),
        "half_widths_rad_s": half_widths,
        "visibility": visibility,
        "stationary_mask": stationary,
        "stationary_residue_max": stationary_residue_max,
        "nonstationary_mask": ~stationary,
        "eigenvector_condition": condition,
        "stationary_modes_projected": int(np.count_nonzero(zero)),
    }


def floquet_solution_diagnostics(L0_at_Deff_zero, Cp, Cm, Omega_beat,
                                 Delta_eff_axis, S_v, harmonics):
    """Trace and normalized equation residuals for a Floquet solution stack."""
    harmonics = np.asarray(harmonics, dtype=complex)
    if harmonics.ndim != 4 or harmonics.shape[1] % 2 != 1:
        raise ValueError(
            "harmonics must have shape (batch, 2*n_f+1, n_levels, n_levels)")
    n_deff, n_sector, n_levels, n_levels_2 = harmonics.shape
    if n_levels != n_levels_2:
        raise ValueError("harmonic density matrices must be square")
    deff = np.asarray(Delta_eff_axis, dtype=float)
    if deff.shape != (n_deff,):
        raise ValueError("Delta_eff_axis length must match the harmonic batch")
    M = n_levels * n_levels
    n_f = (n_sector - 1) // 2
    L0 = np.asarray(L0_at_Deff_zero, dtype=complex)
    Cp = np.asarray(Cp, dtype=complex)
    Cm = np.asarray(Cm, dtype=complex)
    S_v = np.asarray(S_v, dtype=complex)
    omega = np.asarray(Omega_beat, dtype=float)
    finite_inputs = all(np.isfinite(value).all()
                        for value in (L0, Cp, Cm, S_v, omega, deff, harmonics))
    nonfinite_entries = sum(
        int(np.size(value) - np.count_nonzero(np.isfinite(value)))
        for value in (L0, Cp, Cm, S_v, omega, deff, harmonics)
    )
    if not finite_inputs:
        return {
            "n_f": n_f,
            "finite": False,
            "nonfinite_entries": nonfinite_entries,
            "max_normalized_residual": float("inf"),
            "max_trace_error": float("inf"),
            "max_nonzero_harmonic_trace": float("inf"),
            "max_rho0_hermiticity_error": float("inf"),
        }
    L_batch = L0[None, :, :] - deff[:, None, None] * S_v[None, :, :]
    vectors = harmonics.reshape(n_deff, n_sector, M, 1)
    eye = _eye(M)
    worst = 0.0
    for q, harmonic in enumerate(range(-n_f, n_f + 1)):
        central = ((L_batch + 1j * harmonic * float(Omega_beat)
                    * eye[None, :, :]) @ vectors[:, q])
        scale = np.linalg.norm(central[..., 0], axis=-1)
        residual = central
        if q > 0:
            term = Cp[None, :, :] @ vectors[:, q - 1]
            residual = residual + term
            scale = scale + np.linalg.norm(term[..., 0], axis=-1)
        if q + 1 < n_sector:
            term = Cm[None, :, :] @ vectors[:, q + 1]
            residual = residual + term
            scale = scale + np.linalg.norm(term[..., 0], axis=-1)
        normalized = np.linalg.norm(residual[..., 0], axis=-1) / np.maximum(
            scale, np.finfo(float).tiny)
        worst = max(worst, float(np.max(normalized, initial=0.0)))

    traces = np.trace(harmonics, axis1=-2, axis2=-1)
    target = np.zeros_like(traces)
    target[:, n_f] = 1.0
    return {
        "n_f": n_f,
        "finite": True,
        "nonfinite_entries": 0,
        "max_normalized_residual": worst,
        "max_trace_error": float(np.max(np.abs(traces - target), initial=0.0)),
        "max_nonzero_harmonic_trace": float(np.max(
            np.abs(np.delete(traces, n_f, axis=1)), initial=0.0)),
        "max_rho0_hermiticity_error": float(np.max(
            np.abs(harmonics[:, n_f]
                   - harmonics[:, n_f].conj().swapaxes(-1, -2)), initial=0.0)),
    }


def steady_state_batched(L0_at_Deff_zero, Delta_eff_axis, S_v, n_levels):
    """
    Single-mode steady state ρ (Lρ = 0, trace = 1), batched over every Δ_eff.

    Same velocity-batching idea as `floquet_solve` but without the sideband
    coupling — this is the engine for the absorption-cluster schemes (OD / AT /
    EIT / CPT), where a (weak) probe sits inside H₀ and we want the steady-state
    coherence.

    Solved in the **real Hermitian-generator basis** (see `hermitian_basis`):
    L₀ and S_v are Hermiticity-preserving, so one change of basis (done once, not
    per Δ_eff) turns the whole batch into a REAL linear solve of the same
    dimension — ~2–4× fewer flops and half the memory of the complex |i⟩⟨j| form,
    bit-for-bit the same ρ (the real solve matches the old complex solve to
    machine precision). The trace row is Σ_{s<n} r_s = 1 (populations are the
    first n real coordinates); ρ is reconstructed as vec(ρ) = U·r at the end.

    Returns ρ reshaped to (N_deff, n_levels, n_levels).
    """
    n = Delta_eff_axis.size
    nl = n_levels
    M = nl * nl
    U = hermitian_basis(nl)
    Uh = U.conj().T
    L0r = np.real(Uh @ L0_at_Deff_zero @ U)             # one-time → real frame
    S_vr = np.real(Uh @ S_v @ U)
    A = (L0r[None, :, :] - Delta_eff_axis[:, None, None] * S_vr[None, :, :])
    A[:, 0, :] = 0.0
    A[:, 0, :nl] = 1.0                                  # trace: Σ populations = 1
    rhs = np.zeros((n, M, 1))
    rhs[:, 0, 0] = 1.0
    r = np.linalg.solve(A, rhs)[:, :, 0]                # REAL batched solve
    rho_vec = r @ U.T                                   # back to complex vec(ρ)
    return rho_vec.reshape(n, nl, nl)


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


# float64 overflows exp/cosh past ~709; cap the per-exponent real part well
# below that so a runaway *linear* gain stays finite (the physical bound is then
# enforced by the Manley-Rowe pump-depletion saturation downstream). Two terms
# at the cap multiply, so 2·cap must stay under ~709 → 350.
_EXP_ARG_CLAMP = 350.0


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
    # Clamp the real part of the exponent arguments so cosh/sinh/exp cannot
    # overflow to +inf. At extreme density/length/coupling the linear gain this
    # represents is unphysical (it overshoots the pump's energy budget); without
    # the clamp it returns inf, and the downstream pump-depletion saturation then
    # evaluates inf/(1+inf) = NaN, silently poisoning the whole gain/squeezing
    # curve. A huge-but-finite value is instead capped to the energy bound by
    # that saturation. The clamp is a no-op in the validated regime (|Re·L| ≪ cap).
    cL = c * L
    sL = s * L
    cL = np.clip(cL.real, -_EXP_ARG_CLAMP, _EXP_ARG_CLAMP) + 1j * cL.imag
    sL = np.clip(sL.real, -_EXP_ARG_CLAMP, _EXP_ARG_CLAMP) + 1j * sL.imag
    sinh_over_c = np.where(big, np.sinh(cL) / safe_c, L * np.ones_like(c))
    cosh_cL = np.cosh(cL)
    exp_sL = np.exp(sL)

    out = np.empty_like(M)
    out[..., 0, 0] = exp_sL * (cosh_cL + sinh_over_c * q00)
    out[..., 0, 1] = exp_sL * (sinh_over_c * q01)
    out[..., 1, 0] = exp_sL * (sinh_over_c * q10)
    out[..., 1, 1] = exp_sL * (cosh_cL + sinh_over_c * q11)
    return out
