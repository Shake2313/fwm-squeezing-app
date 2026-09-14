"""Full density-operator response and QRT for a moving pump-state atom.

No atomic drift, jump diffusion, production field elimination or velocity
aggregation imports. The readouts are decomposed into Hermitian components.
"""

import numpy as np

from .atomic_qrt import qrt_ordered_spectrum


def pump_state_field_reference(generator, readouts, couplings, signs, carrier_harmonics,
                                beat_rad_s, lab_rf_rad_s, linear_density):
    l = np.asarray(generator, complex)
    ops = np.asarray(readouts, complex)
    n, count = ops.shape[-1], len(ops)
    trace = np.eye(n).reshape(-1)
    equation = l.copy(); equation[0] = trace
    rhs = np.zeros(n*n, complex); rhs[0] = 1.
    rho = np.linalg.solve(equation, rhs).reshape(n, n)
    result = [np.zeros((len(lab_rf_rad_s), count, count), complex) for _ in range(3)]
    for harmonic in np.unique(carrier_harmonics):
        indices = np.flatnonzero(carrier_harmonics == harmonic)
        selected = ops[indices]
        hermitian = np.array([o for op in selected for o in ((op+op.conj().T)/2, (op-op.conj().T)/(2j))])
        weight = np.zeros((len(indices), len(hermitian)), complex)
        for j, idx in enumerate(indices):
            weight[j, 2*j:2*j+2] = -1j*signs[idx]*couplings[idx]*np.array([1., 1j])
        frequencies = np.asarray(lab_rf_rad_s)+harmonic*beat_rad_s
        ordered = qrt_ordered_spectrum(l, rho, hermitian, frequencies)
        reversed_order = qrt_ordered_spectrum(l, rho, hermitian, -frequencies).swapaxes(-1, -2)
        dp = linear_density*(weight@ordered@weight.conj().T)
        dm = linear_density*(weight@reversed_order@weight.conj().T)
        for fi, w in enumerate(frequencies):
            bordered = np.zeros((n*n+1, n*n+1), complex)
            bordered[:-1, :-1] = -l-1j*w*np.eye(n*n)
            bordered[:-1, -1], bordered[-1, :-1] = rho.reshape(-1), trace
            driving = np.zeros((n*n+1, len(indices)), complex)
            for j, idx in enumerate(indices):
                op = ops[idx].conj().T
                driving[:-1, j] = (-1j*couplings[idx]*(op@rho-rho@op)).reshape(-1)
            response = np.linalg.solve(bordered, driving)[:-1]
            for i, idx in enumerate(indices):
                values = -1j*signs[idx]*couplings[idx]*linear_density*(ops[idx].T.reshape(-1)@response)
                result[0][fi, idx, indices] = values
                result[1][fi, idx, indices] = dp[fi, i]
                result[2][fi, idx, indices] = dm[fi, i]
    return {'drift': result[0], 'greater': result[1], 'lesser': result[2], 'state': rho}
