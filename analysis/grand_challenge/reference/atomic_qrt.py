"""Independent density-operator QRT and adjoint-Einstein references.

No import from quantum.diffusion and no use of its atomic drift or noise.
The physical inputs H/L, stationary rho and measured operators are shared.
Frequency convention here is integral exp(+i*omega*tau) C(tau) d tau.
"""

import numpy as np


def _inputs(generator, rho, operators):
    rho, generator, operators = (np.asarray(v, complex) for v in (rho, generator, operators))
    if rho.ndim != 2 or rho.shape[0] != rho.shape[1]:
        raise ValueError("square state required")
    n = rho.shape[0]
    if (generator.shape != (n*n, n*n) or operators.ndim != 3
            or operators.shape[1:] != (n, n) or not len(operators)):
        raise ValueError("reference dimensions do not match")
    if not all(np.isfinite(v).all() for v in (rho, generator, operators)):
        raise ValueError("reference data must be finite")
    if (abs(np.trace(rho)-1) > 1e-10 or np.linalg.norm(rho-rho.conj().T) > 1e-10
            or np.linalg.eigvalsh((rho+rho.conj().T)/2)[0] < -1e-10):
        raise ValueError("physical trace-one reference state required")
    if np.linalg.norm(operators-operators.conj().swapaxes(-1, -2)) > 1e-10:
        raise ValueError("Hermitian reference operators required")
    return generator, rho, operators


def adjoint_einstein_diffusion(generator, rho, operators):
    """D_ij=<L*(Fi Fj)-(L*Fi)Fj-Fi(L*Fj)>; no jump formula."""
    generator, rho, operators = _inputs(generator, rho, operators)
    n = rho.shape[0]

    def adjoint(operator):
        return (generator.conj().T@operator.reshape(-1)).reshape(n, n)

    action = [adjoint(operator) for operator in operators]
    result = np.zeros((len(operators), len(operators)), complex)
    for i, left in enumerate(operators):
        for j, right in enumerate(operators):
            result[i, j] = np.trace(rho@(
                adjoint(left@right)-action[i]@right-left@action[j]))
    return result


def qrt_ordered_spectrum(generator, rho, operators, generator_omega_rad_s):
    """Full Liouville-space, trace-bordered QRT integral including both times.

    For tau>=0, C_ij(tau)=Tr[delta Fi exp(L*tau)(delta Fj rho)].
    The negative-time part is supplied by Hermitian conjugation of the
    one-sided integral, not by discarding imaginary cross correlations.
    """
    generator, rho, operators = _inputs(generator, rho, operators)
    raw_omega = np.asarray(generator_omega_rad_s)
    if np.iscomplexobj(raw_omega) and np.any(raw_omega.imag != 0):
        raise ValueError("real generator frequencies required")
    omega = np.asarray(raw_omega.real, float)
    if omega.ndim != 1 or not omega.size or not np.isfinite(omega).all():
        raise ValueError("finite generator-frequency vector required")
    n, m = rho.shape[0], len(operators)
    scale = max(float(np.linalg.norm(generator, 2)), np.finfo(float).tiny)
    if np.linalg.norm(generator@rho.reshape(-1)) > 1e-10*scale*np.linalg.norm(rho):
        raise ValueError("QRT reference state must be stationary")
    poles = np.linalg.eigvals(generator)
    resolution = 64*np.finfo(float).eps*(n*n)*scale
    stationary = np.abs(poles) <= resolution
    if (np.count_nonzero(stationary) != 1
            or np.any(poles[~stationary].real >= -resolution)):
        raise ValueError("nondecaying reference correlations require explicit elastic terms")
    means = np.array([np.trace(operator@rho) for operator in operators])
    centered = operators-means[:, None, None]*np.eye(n)
    source = np.column_stack([(operator@rho).reshape(-1) for operator in centered])
    readout = np.vstack([operator.T.reshape(-1) for operator in centered])
    trace = np.eye(n).reshape(-1)
    # Remove only roundoff in the known zero-trace sources, along rho_ss.
    source -= np.outer(rho.reshape(-1), trace@source)
    rhs = np.zeros((n*n+1, m), complex)
    rhs[:n*n] = source/scale
    result = []
    for frequency in omega:
        bordered = np.zeros((n*n+1, n*n+1), complex)
        bordered[:n*n, :n*n] = (-generator-1j*frequency*np.eye(n*n))/scale
        bordered[:n*n, n*n] = rho.reshape(-1)
        bordered[n*n, :n*n] = trace
        solution = np.linalg.solve(bordered, rhs)[:n*n]
        one_sided = readout@solution
        result.append(one_sided+one_sided.conj().T)
    return np.array(result)
