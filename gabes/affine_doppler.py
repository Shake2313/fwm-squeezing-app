"""
Doppler averages of a scan whose Liouvillian is affine in scan and velocity:

    L(s, kv) = base + s·A + kv·S_v

The Λ (two-photon detuning) and Rydberg ladder scans have this form whenever
`build_H0(s)` is affine in the scan variable s. `AffineScan` builds the three
constant operators once; `grid_mean` returns the probe coherence ρ_eg / Ω_p
averaged over a discrete set of velocity classes. Every scan point is a rational
function of kv (`pole_doppler`), so a Maxwell grid costs one small
eigendecomposition per point; a point that misses the guard is summed class by
class (the compiled affine kernel when Numba is available, a batched solve
otherwise), as is the single stationary class. `maxwell_shift_classes` gives the
classes for a Maxwell vapor, or kv = 0 with Doppler off.
"""
import numpy as np

from . import core, doppler, kernels, pole_doppler

POLE_GUARD_RTOL = 1e-6      # pole form vs a direct solve of the v = 0 class, every row
POLE_LIMIT_SIGMAS = 1e8     # farther poles are rounding modes; 1e4 already cuts a real one


def maxwell_shift_classes(T, mass, k_vec, doppler_on, dv=2.0, cutoff_sigma=4.0):
    """(kv [rad/s], weights) of the Maxwell velocity classes; Doppler off → kv = 0."""
    if not doppler_on:
        return np.zeros(1), np.ones(1)
    v, w = doppler.velocity_grid(T, mass=mass, dv=dv, cutoff_sigma=cutoff_sigma)
    return k_vec * v, w


class AffineScan:
    """Constant pieces of L(s, kv) = base + s·A + kv·S_v for one probe readout.

    `scan_shift` is the optical scan operator, entering as −s·scan_shift: the Λ
    scan moves the optical detuning through it while a residual Raman Doppler
    term belongs only in S_v and must not change a stationary atom's scan
    Hamiltonian. Pass None when the scan already lives in `build_H0(s)` (the
    Rydberg probe detuning). The s-independent drive entries are bit-identical
    between build_H0(0) and build_H0(1), so their difference isolates dL/ds
    exactly with no cancellation.
    """

    def __init__(self, atom, build_H0, coherence, probe_omega, scan_shift=None):
        self.atom = atom
        self.build_H0 = build_H0
        self.coherence = coherence
        self.probe_omega = probe_omega
        self.scan_shift = scan_shift
        self.base = core.build_liouvillian(build_H0(0.0), atom)
        dL_ds = core.build_liouvillian(build_H0(1.0), atom) - self.base
        self.A = dL_ds if scan_shift is None else dL_ds - scan_shift

    def _kernel(self, scan, kv, w):
        e, g = self.coherence
        n = self.atom.n_levels
        with core.blas_single_thread():
            out = kernels.affine_scan_chi(
                self.base, self.A, self.atom.S_v,
                np.ascontiguousarray(scan, dtype=float),
                np.ascontiguousarray(kv, dtype=float),
                np.ascontiguousarray(w, dtype=float), e * n + g, n)
        return out / self.probe_omega

    def grid_mean(self, scan, kv, w):
        """χ̄(scan) = Σ_j w_j ρ_eg(s, kv_j) / Ω_p over the given velocity classes."""
        sigma = float(np.sqrt(np.sum(w * kv ** 2) / np.sum(w)))
        if sigma == 0.0:
            return self._class_sum(scan, kv, w)
        try:
            out, failed = self._pole_mean(scan, kv, w, sigma)
        except np.linalg.LinAlgError:
            return self._class_sum(scan, kv, w)
        if failed.any():
            out[failed] = self._class_sum(scan[failed], kv, w)
        return out

    def _pole_mean(self, scan, kv, w, sigma):
        """Pole-form `grid_mean` and the mask of rows its guard rejects.

        With the trace row inserted a row is (A00 + s·A01 − D·S) x = e₀ in the real
        Hermitian-generator basis, with D = −kv. A row is kept if its pole form
        reproduces a direct solve of the v = 0 class to POLE_GUARD_RTOL and has no
        undamped pole inside the velocity grid.
        """
        e, g = self.coherence
        n = self.atom.n_levels
        A00, A01, S = (np.ascontiguousarray(core.to_real_liouvillian(X, n))
                       for X in (self.base, self.A, self.atom.S_v))
        for X in (A00, A01, S):
            X[0] = 0.0
        A00[0, :n] = 1.0                 # trace row: populations are the first n coordinates
        rhs = np.zeros(n * n)
        rhs[0] = 1.0
        # Expand about kv = +σ, not by Schur elimination: the ground populations
        # relax only through the optical coherences, so the Doppler-free block alone
        # is singular.
        system = pole_doppler.ShiftInvertSystem(A00, A01, S, rhs, core.hermitian_basis(n)[e * n + g])
        with core.blas_single_thread():
            poles = system.residues(scan, -sigma, pole_limit=POLE_LIMIT_SIGMAS * sigma)
            at_rest = system.direct(scan, 0.0)[:, 0]
        rest = poles.value_at(0.0)[:, 0]
        mean = poles.discrete_mean(-kv, w)[:, 0]
        in_grid = np.abs(poles.lam.real) <= np.max(np.abs(kv))
        damping = np.min(np.where(in_grid, np.abs(poles.lam.imag), np.inf), axis=1)
        failed = ~((np.abs(rest - at_rest) <= POLE_GUARD_RTOL * np.abs(at_rest))
                   & (damping > 1e-9 * sigma) & np.isfinite(mean))
        return mean / self.probe_omega, failed

    def _class_sum(self, scan, kv, w):
        """`grid_mean` solved class by class."""
        if kernels.available():
            return self._kernel(scan, kv, w)
        e, g = self.coherence
        n = self.atom.n_levels
        out = np.zeros(scan.size, dtype=complex)
        for i, s in enumerate(scan):
            L0 = core.build_liouvillian(self.build_H0(s), self.atom)
            if self.scan_shift is not None:
                L0 = L0 - s * self.scan_shift
            # Δ = −kv reproduces the kernel's +kv·S_v velocity shift.
            rho = core.steady_state_batched(L0, -kv, self.atom.S_v, n)
            out[i] = ((rho[:, e, g] / self.probe_omega) * w).sum()
        return out

    def stationary(self, scan):
        """χ̄(scan) of the v = 0 class with the scan carried by `scan_shift`."""
        if kernels.available():
            return self._kernel(scan, np.zeros(1), np.ones(1))
        e, g = self.coherence
        n = self.atom.n_levels
        out = np.zeros(scan.size, dtype=complex)
        for i, s in enumerate(scan):
            L0 = core.build_liouvillian(self.build_H0(s), self.atom)
            rho = core.steady_state_batched(L0, np.array([s]), self.scan_shift, n)
            out[i] = rho[0, e, g] / self.probe_omega
        return out
