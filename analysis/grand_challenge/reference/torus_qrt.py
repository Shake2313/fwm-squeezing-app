"""Independent full-density trajectory relaxation and two-time regression.

No mean lattice, atomic basis drift or diffusion imports. Integrate actual
moving phases; finite history/delay cutoffs require separate refinement.
"""

from itertools import product

import numpy as np
from scipy.integrate import solve_ivp

from .periodic_qrt import time_domain_qrt


def trajectory_qrt(generators, generator_labels, frequencies_rad_s, operators, base_frequencies_rad_s,
                   output_labels, *, phase_samples=(12, 12), history=40., delay=40., rtol=2e-10, atol=2e-12):
    generators, labels, frequencies, ops, omega, outputs = [np.asarray(x) for x in
        (generators, generator_labels, frequencies_rad_s, operators, base_frequencies_rad_s, output_labels)]
    n, m = ops.shape[-1], len(ops)
    phases = np.array(list(product(*(2*np.pi*np.arange(count)/count for count in phase_samples))))
    actual = (omega[:, None]+outputs@frequencies).reshape(-1)
    nf, nh, size = len(omega), len(outputs), n*n
    readout = ops.swapaxes(-1, -2).reshape(m, size)
    positive = np.zeros((nf, nh, m, nh, m), complex)
    states = []
    maximum_tail = 0.
    def generator(t, theta):
        return np.einsum('q,qab->ab', np.exp(-1j*labels@(theta+frequencies*t)), generators)
    for theta in phases:
        mean = solve_ivp(lambda t, y: generator(t, theta)@y, (-history, 0.),
            (np.eye(n)/n).astype(complex).reshape(-1), method='DOP853', rtol=rtol, atol=atol)
        if not mean.success:
            raise ValueError('independent history integration failed')
        rho = mean.y[:, -1].reshape(n, n)
        states.append(rho)
        mu = np.einsum('iab,ba->i', ops, rho)
        source = ((ops-mu[:, None, None]*np.eye(n))@rho).reshape(m, size).T
        initial = np.zeros((1+len(actual), size, m), complex)
        initial[0] = source
        def ode(t, raw):
            values = raw.reshape(initial.shape)
            out = np.empty_like(values)
            out[0] = generator(t, theta)@values[0]
            out[1:] = np.exp(1j*actual*t)[:, None, None]*values[0]
            return out.reshape(-1)
        correlation = solve_ivp(ode, (0., delay), initial.reshape(-1),
            method='DOP853', rtol=rtol, atol=atol)
        if not correlation.success:
            raise ValueError('independent regression integration failed')
        result = correlation.y[:, -1].reshape(initial.shape)
        maximum_tail = max(maximum_tail, float(np.linalg.norm(result[0])))
        integrals = (readout@result[1:]).reshape(nf, nh, m, m)
        for i, a in enumerate(outputs):
            for j, b in enumerate(outputs):
                positive[:, i, :, j, :] += integrals[:, i]*np.exp(1j*(a-b)@theta)/len(phases)
    positive = positive.reshape(nf, nh*m, nh*m)
    return {'greater': positive+positive.conj().swapaxes(-1, -2), 'phase_states': np.array(states),
        'phases_rad': phases, 'maximum_regression_endpoint_norm': maximum_tail,
        'phase_samples': list(phase_samples), 'history': history, 'delay': delay, 'rtol': rtol, 'atol': atol,
        'scope': 'direct trajectories with finite history/delay; refine all three controls before comparison'}


def static_grating_qrt(l0, probe_pair, conjugate_pair, beat_rad_s, operators,
                        base_frequencies_rad_s, output_labels, *, grating_samples=8, phase_samples=16):
    """At v.Q=0 integrate each fixed grating phase over one exact beat period.

    Resum the time-domain QRT tail geometrically, then Fourier-project the
    independent spatial phase. It is not averaged away before solving the atom.
    """
    outputs = np.asarray(output_labels, int)
    harmonics = np.unique(outputs[:, 0])
    m = len(operators)
    out = np.zeros((len(base_frequencies_rad_s), len(outputs)*m, len(outputs)*m), complex)
    phase_rows, states = [], []
    for psi in 2*np.pi*np.arange(grating_samples)/grating_samples:
        # Each pair is (L(V), L(V^dagger)); L(V)^dagger is not L(V^dagger).
        lp = probe_pair[0]+conjugate_pair[1]*np.exp(1j*psi)
        lm = probe_pair[1]+conjugate_pair[0]*np.exp(-1j*psi)
        reference = time_domain_qrt(l0, lp, lm, beat_rad_s, operators, base_frequencies_rad_s,
            harmonics, phase_samples=phase_samples)
        states.extend(reference['phase_states'])
        phase_rows.extend(np.column_stack([reference['physical_phases_rad'], np.full(phase_samples, psi)]))
        for i, a in enumerate(outputs):
            hi = np.flatnonzero(harmonics == a[0])[0]
            for j, b in enumerate(outputs):
                hj = np.flatnonzero(harmonics == b[0])[0]
                out[:, i*m:(i+1)*m, j*m:(j+1)*m] += (
                    reference['ordered_spectrum'][:, hi*m:(hi+1)*m, hj*m:(hj+1)*m]
                    *np.exp(1j*(a[1]-b[1])*psi)/grating_samples)
    return {'greater': out, 'phases_rad': np.array(phase_rows), 'phase_states': np.array(states),
        'grating_samples': grating_samples, 'phase_samples': phase_samples,
        'scope': 'zero loop convection only; independent period QRT at each spatial phase with exact geometric tail'}
