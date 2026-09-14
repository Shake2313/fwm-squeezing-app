"""Independent full-Liouville forced harmonic response and ordered current.

Does not import atomic drift, field elimination, quadrature conversion or
photocurrent production formulas. Trace-zero constraints are explicit per block.
"""

import numpy as np


def forced_liouville_response(l0, lp, lm, beat, states, operators, harmonics, couplings,
                              signs, frequencies, *, response_order, linear_density):
    n, count = states.shape[-1], len(operators)
    size, block = n*n, n*n+1
    hs = np.arange(-response_order, response_order+1)
    order = len(states)//2
    trace = np.eye(n).reshape(-1)
    source = np.zeros((len(hs)*block, count), complex)
    for hi, harmonic in enumerate(hs):
        for j, (op, field_h, g) in enumerate(zip(operators, harmonics, couplings)):
            q = int(harmonic-field_h)
            if abs(q) <= order:
                rho = states[q+order]
                source[hi*block:hi*block+size, j] = (-1j*g*(op.conj().T@rho-rho@op.conj().T)).reshape(-1)
    results = []
    for w in frequencies:
        matrix = np.zeros((len(hs)*block,)*2, complex)
        for hi, h in enumerate(hs):
            i = hi*block
            matrix[i:i+size, i:i+size] = -1j*(w+h*beat)*np.eye(size)-l0
            matrix[i:i+size, i+size] = trace/n
            matrix[i+size, i:i+size] = trace
            if hi:
                matrix[i:i+size, i-block:i-block+size] = -lp
            if hi+1 < len(hs):
                matrix[i:i+size, i+block:i+block+size] = -lm
        response = np.linalg.solve(matrix, source)
        results.append(np.array([-1j*s*g*linear_density*(op.T.reshape(-1)@response[int(h+response_order)*block:int(h+response_order)*block+size])
                                for op, h, g, s in zip(operators, harmonics, couplings, signs)]))
    return np.array(results)


def ordered_current_psd(greater, beta, eta, response, balance, elementary_charge):
    """Direct y=beta* a(+RF)+beta a†(-RF), with loss vacuum, one-sided PSD."""
    transmission = np.tile(eta, 2)
    output = np.sqrt(transmission)[:, None]*greater*np.sqrt(transmission)[None, :]
    output = output+np.diag(np.r_[1-np.asarray(eta), [0., 0.]])
    bp, bc = np.sqrt(eta)*beta
    hp, hc = response
    row = elementary_charge*np.array([hp*bp.conjugate(), -balance*hc*bc.conjugate(), hp*bp, -balance*hc*bc])
    return float(2*(row@output@row.conj()).real)
