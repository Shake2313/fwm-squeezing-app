"""
Doppler averages of a scan whose Liouvillian is affine in scan and velocity:

    L(s, kv) = base + s·A + kv·S_v

The Λ (two-photon detuning) and Rydberg ladder scans have this form whenever
`build_H0(s)` is affine in the scan variable s. `AffineScan` builds the three
constant operators once; `grid_mean` returns the probe coherence ρ_eg / Ω_p
averaged over a discrete set of velocity classes (the compiled affine kernel
when Numba is available, a per-scan-point batched solve otherwise).
`maxwell_shift_classes` gives those classes for a Maxwell vapor, or the single
stationary class kv = 0 with Doppler off.
"""
import numpy as np

from . import core, doppler, kernels


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
