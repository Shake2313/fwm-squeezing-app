"""Independent I x J=1/2 dipole reference, with no CF2/6j/GABES Zeeman import.

Closed-form I +/- 1/2 spinors transform Pauli spherical tensors. Overall
electronic phase is arbitrary; strengths and closure identities are not.
"""

import numpy as np


def d1_dipole_reference():
    # Uncoupled basis |mI, mJ>, with mJ=+1/2,-1/2 in each nuclear block.
    nuclear = np.arange(-2.5, 3., 1.)
    labels = [(F, m) for F in (2, 3) for m in range(-F, F+1)]
    U = np.zeros((12, 12))
    for k, (F, m) in enumerate(labels):
        plus, minus = np.sqrt((3+m)/6), np.sqrt((3-m)/6)
        coeff = (plus, minus) if F == 3 else (-minus, plus)
        for j, mj in enumerate((.5, -.5)):
            matches = np.flatnonzero(nuclear == m-mj)
            if len(matches):
                U[2*matches[0]+j, k] = coeff[j]
    tensors = [np.array([[0., 0.], [np.sqrt(2/3), 0.]]),
               np.diag([1., -1.])/np.sqrt(3),
               np.array([[0., -np.sqrt(2/3)], [0., 0.]])]
    dipoles = np.array([U.T@np.kron(np.eye(6), t)@U for t in tensors])
    # Rows excited, columns ground; in units of Steck dJ.
    rows = []
    for Fg in (2, 3):
        gi = [k for k, x in enumerate(labels) if x[0] == Fg]
        for Fe in (2, 3):
            ei = [k for k, x in enumerate(labels) if x[0] == Fe]
            strength = np.array([np.sum(abs(d[np.ix_(ei, gi)])**2)/(2*Fg+1) for d in dipoles])
            rows.append({"Fg": Fg, "Fe": Fe, "mean_strength_by_q": strength.tolist()})
    emission = sum(d@d.T for d in dipoles)
    absorption = sum(d.T@d for d in dipoles)
    return {"labels": labels, "dipoles": dipoles,
            "basis_unitarity_error": float(np.max(abs(U.T@U-np.eye(12)))),
            "emission_closure_error": float(np.max(abs(emission-np.eye(12)))),
            "absorption_closure_error": float(np.max(abs(absorption-np.eye(12)))),
            "transition_rows": rows}
