"""Local, prescribed-carrier Floquet test of the pump-only weak-seed state.

This diagnoses back-action. It does not replace stationary diffusion by an
averaged Floquet state or claim to propagate finite-seed quantum noise.
"""

import numpy as np

from .. import constants as c, core, observables
from ..quantum.contracts import readonly_array
from .field import reduced_readout_operators
from .normalization import ReducedDipoles


def seed_harmonic_hamiltonian(dipoles, optical_omega_rad_s, uniform_area_m2, carrier_amplitudes):
    """H/hbar=H0+V exp(-i*beat*t)+V† exp(+i*beat*t), with beat=-hf+delta."""
    if not isinstance(dipoles, ReducedDipoles):
        raise TypeError("explicit ReducedDipoles required")
    omega = readonly_array(optical_omega_rad_s, real=True)
    beta = readonly_array(carrier_amplitudes)
    area = float(uniform_area_m2)
    if (omega.shape != (2,) or np.any(omega <= 0) or beta.shape != (2,)
            or not np.isfinite(area) or area <= 0):
        raise ValueError("two carriers, two optical frequencies and positive uniform area required")
    q = np.diag(observables.photon_flux_mode_matrix(*omega, area, area))
    g = dipoles.base_dipole_C_m*q/(2*c.HBAR)
    ops = reduced_readout_operators(dipoles.transition_scales)
    return readonly_array(g[0]*beta[0]*ops[0].conj().T+g[1]*beta[1].conjugate()*ops[1].conj().T)


def local_seed_backaction(atom, harmonic_hamiltonian_rad_s, beat_rad_s, *, n_f, phase_samples=96):
    """Solve periodic mean state and compare its driven harmonic to weak response.

    Positivity is sampled in phase, not proven at all times. Residual outside
    the truncated harmonic ladder is included, along with the boundary solve.
    """
    v = readonly_array(harmonic_hamiltonian_rad_s)
    beat = float(beat_rad_s)
    if v.shape != (4, 4) or not np.isfinite(beat) or beat == 0:
        raise ValueError("4-level finite harmonic drive and nonzero beat required")
    if isinstance(n_f, bool) or int(n_f) != n_f or n_f < 1:
        raise ValueError("positive integer Floquet order required")
    if isinstance(phase_samples, bool) or int(phase_samples) != phase_samples or phase_samples < 2*n_f+1:
        raise ValueError("enough integer phase samples required")
    cp, cm = core.comm_super(v), core.comm_super(v.conj().T)
    zero = np.zeros((16, 16), complex)
    harmonics = core.floquet_solve_truncated(atom.generator, cp, cm, beat, [0.], zero, 4,
                                            n_f=n_f, return_harmonics=True)[0]
    diag = core.floquet_solution_diagnostics(atom.generator, cp, cm, beat, [0.], zero, harmonics[None])
    # Trace-bordered independent first-order equation; no inverse of the drift.
    equation = np.zeros((17, 17), complex)
    equation[:16, :16] = atom.generator+1j*beat*np.eye(16)
    trace = np.eye(4).reshape(-1)
    equation[:16, 16], equation[16, :16] = trace, trace
    rhs = np.r_[-cp@atom.stationary_state.reshape(-1), 0.]
    weak = np.linalg.solve(equation, rhs)[:16].reshape(4, 4)
    difference = harmonics[n_f+1]-weak
    weak_norm = float(np.linalg.norm(weak))
    mean_change = harmonics[n_f]-atom.stationary_state
    phases = 2*np.pi*np.arange(phase_samples)/phase_samples
    rho = np.einsum("tn,nij->tij", np.exp(-1j*np.outer(phases, np.arange(-n_f, n_f+1))), harmonics)
    minimum = float(np.linalg.eigvalsh((rho+rho.conj().swapaxes(-1, -2))/2).min())
    omitted = max(np.linalg.norm(cp@harmonics[-1].reshape(-1)),
                  np.linalg.norm(cm@harmonics[0].reshape(-1)))
    scale = max(np.linalg.norm(atom.generator)*np.linalg.norm(harmonics), np.finfo(float).tiny)
    # At tiny seed powers the zero-harmonic forcing is O(P); dividing its
    # cancellation residual by that forcing alone spuriously fails as P->0.
    ratios = []
    for k, harmonic in enumerate(range(-n_f, n_f+1)):
        terms = [(atom.generator+1j*harmonic*beat*np.eye(16))@harmonics[k].reshape(-1)]
        cancellation_scale = (np.linalg.norm(atom.generator)+abs(harmonic*beat))*np.linalg.norm(harmonics[k])
        for neighbor, op in ((k-1, cp), (k+1, cm)):
            if 0 <= neighbor < len(harmonics):
                terms.append(op@harmonics[neighbor].reshape(-1))
                cancellation_scale += np.linalg.norm(op)*np.linalg.norm(harmonics[neighbor])
        tolerance = 1e-9*sum(np.linalg.norm(t) for t in terms)+128*np.finfo(float).eps*cancellation_scale
        ratios.append(np.linalg.norm(sum(terms))/max(tolerance, np.finfo(float).tiny))
    diag.update({"phase_sample_count": phase_samples, "phase_sampled_minimum_state_eigenvalue": minimum,
                 "phase_sampled_maximum_hermiticity_error": float(np.linalg.norm(rho-rho.conj().swapaxes(-1, -2), axis=(-2, -1)).max()),
                 "omitted_boundary_relative_residual": float(omitted/scale),
                 "maximum_cancellation_aware_residual_ratio": float(max(ratios)),
                 "cancellation_aware_residual_passed": bool(max(ratios) <= 1.)})
    return {"harmonics": readonly_array(harmonics), "weak_harmonic": readonly_array(weak),
            "mean_state_trace_distance": float(np.abs(np.linalg.eigvalsh((mean_change+mean_change.conj().T)/2)).sum()/2),
            "first_harmonic_relative_error": float(np.linalg.norm(difference)/weak_norm) if weak_norm else 0.,
            "first_harmonic_absolute_error": float(np.linalg.norm(difference)),
            "mean_excited_population": float(np.trace(harmonics[n_f, 2:, 2:]).real),
            "diagnostics": diag}
