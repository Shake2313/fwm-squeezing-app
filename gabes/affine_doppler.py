"""
Doppler averages of a scan whose Liouvillian is affine in scan and velocity:

    L(s, kv) = base + s·A + kv·S_v

The Λ (two-photon detuning) and Rydberg ladder scans have this form whenever
`build_H0(s)` is affine in the scan variable s. `AffineScan` builds the three
constant operators once. Every scan point is a rational function of kv
(`pole_doppler`), so `maxwell_mean` returns the probe coherence ρ_eg / Ω_p averaged
over the Maxwell distribution in closed form: one small eigendecomposition per
point, no velocity grid and no cutoff. A point that misses the guard is summed on
a grid fine enough for its poles (the compiled affine kernel when Numba is
available, a batched solve otherwise). `stationary` is the single v = 0 class.
"""
import math

import numpy as np

from . import constants, core, doppler, kernels, pole_doppler

POLE_GUARD_RTOL = 1e-6      # pole form vs a direct solve of the v = 0 class, every row
POLE_LIMIT_SIGMAS = 1e8     # farther poles are rounding modes; 1e4 already cuts a real one
BULK_SIGMAS = 6.0           # velocity span where a pole must be damped and a grid reaches
STRIP_FLOOR_SIGMAS = 1e-4   # narrowest pole a grid resolves; real ones stay above 9e-4 σ
QUADRATURE_ERROR = 1e-8     # nominal trapezoid error of a resolving grid


def _shift_sigma(T, mass, k_vec):
    """Maxwell standard deviation of kv [rad/s]."""
    return k_vec * math.sqrt(constants.KB * T / mass)


def _strip_classes(damping, T, mass, k_vec):
    """Maxwell (v, weights) resolving the narrowest of the row `damping`s."""
    sigma = _shift_sigma(T, mass, k_vec)
    strip = float(np.clip(np.min(damping), STRIP_FLOOR_SIGMAS * sigma, sigma))
    return doppler.maxwell_strip_grid(T, mass, k_vec, strip, QUADRATURE_ERROR, BULK_SIGMAS)


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

    def maxwell_mean(self, scan, T, mass, k_vec):
        """χ̄(scan) = E[ρ_eg(s, kv)] / Ω_p for kv = k_vec·v of a Maxwell vapor."""
        sigma = _shift_sigma(T, mass, k_vec)
        try:
            poles, damping, failed = self._pole_rows(scan, sigma)
        except np.linalg.LinAlgError:
            out = np.empty(scan.size, dtype=complex)
            damping = np.zeros(scan.size)
            failed = np.ones(scan.size, dtype=bool)
        else:
            out = poles.gaussian_mean(0.0, sigma)[:, 0] / self.probe_omega
            failed |= ~np.isfinite(out)
        if failed.any():
            v, w = _strip_classes(damping[failed], T, mass, k_vec)
            out[failed] = self._class_sum(scan[failed], k_vec * v, w)
        return out

    def quadrature_classes(self, scan, T, mass, k_vec):
        """Maxwell (v, weights) resolving the narrowest pole of these scan points, for
        a response there with no closed form (the finite-IF phasor)."""
        return _strip_classes(self._pole_rows(scan, _shift_sigma(T, mass, k_vec))[1],
                              T, mass, k_vec)

    def _pole_rows(self, scan, sigma):
        """Pole rows of the coherence, their damping and the rows the guard rejects.

        With the trace row inserted a row is (A00 + s·A01 − D·S) x = e₀ in the real
        Hermitian-generator basis, with D = −kv. A row is kept if its pole form
        reproduces a direct solve of the v = 0 class to POLE_GUARD_RTOL and has no
        undamped pole inside the Maxwell bulk.
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
        in_bulk = np.abs(poles.lam.real) <= BULK_SIGMAS * sigma
        damping = np.min(np.where(in_bulk, np.abs(poles.lam.imag), np.inf), axis=1)
        failed = ~((np.abs(rest - at_rest) <= POLE_GUARD_RTOL * np.abs(at_rest))
                   & (damping > 1e-9 * sigma))
        return poles, damping, failed

    def _class_sum(self, scan, kv, w):
        """χ̄(scan) = Σ_j w_j ρ_eg(s, kv_j) / Ω_p solved class by class."""
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
        """χ̄(scan) of the v = 0 class."""
        if kernels.available() or self.scan_shift is None:
            return self._class_sum(scan, np.zeros(1), np.ones(1))
        e, g = self.coherence
        n = self.atom.n_levels
        out = np.zeros(scan.size, dtype=complex)
        for i, s in enumerate(scan):
            L0 = core.build_liouvillian(self.build_H0(s), self.atom)
            rho = core.steady_state_batched(L0, np.array([s]), self.scan_shift, n)
            out[i] = rho[0, e, g] / self.probe_omega
        return out
