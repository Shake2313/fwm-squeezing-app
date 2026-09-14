"""Four-mode Fock-space current reference, independent of Gaussian contractions."""

import numpy as np


def two_pair_current_reference(occupation_a, occupation_b, phase_a, phase_b, *, cutoff=14):
    """Modes p_low,p_high,c_low,c_high; pair p_low/c_high and p_high/c_low.

    Compare the one-bin RF offset. Explicit bilinear operators in a truncated
    Fock state supply the fourth moment. Additional vacuum frequency bins supply
    the missing spontaneous shot term at the two optical band edges.
    """
    n_a, n_b = float(occupation_a), float(occupation_b)
    if min(n_a, n_b) < 0 or cutoff < 3:
        raise ValueError("nonnegative populations and sufficient Fock cutoff required")
    la = np.sqrt(n_a/(n_a+1))*np.exp(1j*phase_a)
    lb = np.sqrt(n_b/(n_b+1))*np.exp(1j*phase_b)
    state = np.zeros((cutoff,)*4, complex)
    for k in range(cutoff):
        for ell in range(cutoff):
            state[k, ell, ell, k] = la**k*lb**ell
    state /= np.linalg.norm(state)

    def apply(value, mode, dagger):
        moved = np.moveaxis(value, mode, 0)
        out = np.zeros_like(moved)
        factors = np.sqrt(np.arange(1, cutoff)).reshape((-1,)+(1,)*3)
        if dagger:
            out[1:] = factors*moved[:-1]
        else:
            out[:-1] = factors*moved[1:]
        return np.moveaxis(out, 0, mode)

    # J(+df)† = a_high† a_low in each arm.
    vectors = [apply(apply(state, lo, False), hi, True) for lo, hi in ((0, 1), (2, 3))]
    covariance = np.array([[np.vdot(a, b) for b in vectors] for a in vectors])
    populations = np.array([np.linalg.norm(apply(state, j, False))**2 for j in range(4)])
    # The missing high->outside-vacuum beat contributes n_high to each auto PSD.
    covariance += np.diag(populations[[1, 3]])
    return {"count_covariance_per_hz": covariance,
            "occupations": populations, "omitted_probability_bound": float((n_a/(n_a+1))**cutoff+(n_b/(n_b+1))**cutoff)}
