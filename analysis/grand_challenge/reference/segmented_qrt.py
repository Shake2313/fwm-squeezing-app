"""Independent full-density QRT for a finite piecewise-constant protocol.

Only H, explicit reservoirs, the entry density and readout operators are used.
No atomic drift A, microscopic diffusion D, or propagated covariance enters.
Raw positive-time moments are integrated before global mean subtraction.
"""

import numpy as np
from scipy.linalg import expm

from gabes.quantum.reservoirs import ExplicitReservoirs


def _real_array(value, name):
    raw = np.asarray(value)
    if np.iscomplexobj(raw) and np.any(raw.imag != 0):
        raise ValueError(f"{name} must be real")
    result = np.asarray(raw.real, float)
    if not np.isfinite(result).all():
        raise ValueError(f"{name} must be finite")
    return result


def segmented_qrt(hamiltonians, reservoirs, boundary_state, durations_s,
                  readouts, frequencies_rad_s):
    """Integrate Y_j = integral_0^T exp(i*w_j*t) O_j(t) dt by exponentials.

    Parameters are H[S,n,n] in rad/s, one shared ExplicitReservoirs, a physical
    trace-one boundary_state[n,n], strictly positive durations_s[S], arbitrary
    complex readouts[S,p,n,n], and real frequencies_rad_s[nf,p]. Neither the
    entry state nor any segment needs to be stationary. Frequencies are fixed
    across segments; H and O may jump at their boundaries.

    The returned dict contains connected greater[f,j,k] = <dY_j dY_k^dagger>
    and lesser[f,j,k] = <dY_k^dagger dY_j>, mean_pulse[nf,p], exit_state[n,n],
    and raw_greater/raw_lesser before mean subtraction. Moment units are
    readout units squared times s**2; mean units are readout units times s.
    residence_time_s permits the existing Poisson-pulse conversion. No
    clipping, trace renormalization, or covariance repair is performed.
    """
    if not isinstance(reservoirs, ExplicitReservoirs):
        raise TypeError("one shared ExplicitReservoirs instance required")
    h = np.asarray(hamiltonians, complex)
    ops = np.asarray(readouts, complex)
    rho = np.asarray(boundary_state, complex)
    dt = _real_array(durations_s, "durations_s")
    frequencies = _real_array(frequencies_rad_s, "frequencies_rad_s")
    n = reservoirs.n_levels
    if h.ndim != 3 or h.shape[1:] != (n, n) or not len(h):
        raise ValueError("hamiltonians must have shape (S,n,n) with S > 0")
    segments = len(h)
    if dt.shape != (segments,) or np.any(dt <= 0):
        raise ValueError("durations_s must have shape (S,) and be strictly positive")
    duration = float(dt.sum())
    if not np.isfinite(duration):
        raise ValueError("total duration must be finite")
    if (ops.ndim != 4 or ops.shape[0] != segments
            or ops.shape[2:] != (n, n) or not ops.shape[1]):
        raise ValueError("readouts must have shape (S,p,n,n) with p > 0")
    ports = ops.shape[1]
    if frequencies.ndim != 2 or frequencies.shape[1] != ports or not len(frequencies):
        raise ValueError("frequencies_rad_s must have shape (nf,p) with nf > 0")
    if rho.shape != (n, n):
        raise ValueError("boundary_state must have shape (n,n)")
    if not all(np.isfinite(value).all() for value in (h, ops, rho)):
        raise ValueError("Hamiltonians, readouts and boundary state must be finite")
    if (abs(np.trace(rho)-1) > 1e-10
            or np.linalg.norm(rho-rho.conj().T) > 1e-10
            or np.linalg.eigvalsh((rho+rho.conj().T)/2)[0] < -1e-10):
        raise ValueError("physical trace-one boundary state required")
    # Validate every H before doing any propagation. These are full n**2
    # density generators, not complete-basis atomic drift/covariance maps.
    generators = [reservoirs.generator(value) for value in h]

    n2 = n*n
    x_end = n2 + 2*ports*n2
    z_end = x_end + 2*ports*ports
    dimension = z_end + ports
    x_indices = np.arange(n2, x_end).reshape(2, ports, n2)
    z_indices = np.arange(x_end, z_end).reshape(2, ports, ports)
    mean_indices = np.arange(z_end, dimension)
    states = np.zeros((len(frequencies), dimension), complex)
    states[:, :n2] = rho.ravel()
    exit_density = rho.ravel().copy()

    # Let X_k(t) = integral_0^t U(t,s)[O_k(s)^dagger rho(s)]
    #                         * exp(-i*w_k*s) ds/T (reverse product for lesser).
    # x_k = exp(i*w_k*t) X_k gives x'_k = (L+i*w_k)x_k + source/T.
    # If B_jk is the physical positive-time triangle divided by T**2,
    # Z_jk = exp(-i*(w_j-w_k)*t) B_jk obeys
    # Z'_jk = i*(w_k-w_j) Z_jk + Tr(O_j x_k)/T.
    # m'_j = -i*w_j*m_j + Tr(O_j rho)/T accumulates the rotating mean.
    # All sources are RAW and linear in rho. In u=t/T these equations form
    # one constant homogeneous block system per segment. Thus exp(G*du) is
    # the exact semigroup step (up to floating-point exponential error): no
    # quadrature or repeated fast-phase ODE solve is needed, even at GHz.
    # x, Z and m carry continuously across boundaries because their phase
    # origins are global t=0; jumps in H/O change derivatives, not states.
    for generator, operators, step in zip(generators, ops, dt):
        block = np.zeros((dimension, dimension), complex)
        block[:n2, :n2] = duration*generator
        trace_rows = operators.swapaxes(-1, -2).reshape(ports, n2)
        block[z_end:, :n2] = trace_rows
        for ordering in range(2):
            for k in range(ports):
                xs = slice(x_indices[ordering, k, 0], x_indices[ordering, k, -1]+1)
                adjoint = operators[k].conj().T
                # Row-major vec(A rho B) = (A tensor B.T) vec(rho).
                source = (np.kron(adjoint, np.eye(n)) if ordering == 0
                          else np.kron(np.eye(n), adjoint.T))
                block[xs, :n2] = source
                block[xs, xs] = duration*generator
                block[z_indices[ordering, :, k], xs] = trace_rows
        for f, omega in enumerate(frequencies):
            rotating = block.copy()
            rotating[x_indices, x_indices] += 1j*duration*omega[None, :, None]
            rotating[z_indices, z_indices] += (
                1j*duration*(omega[None, :]-omega[:, None])[None])
            rotating[mean_indices, mean_indices] = -1j*duration*omega
            states[f] = expm(rotating*(step/duration)) @ states[f]
        # A small separate full-density step makes the reported exit state
        # independent of the requested frequency grid/readout amplitudes.
        exit_density = expm(generator*step) @ exit_density

    mean = np.exp(1j*frequencies*duration)*states[:, z_end:]*duration
    triangles = states[:, x_end:z_end].reshape(len(frequencies), 2, ports, ports)
    phase = np.exp(1j*duration*(frequencies[:, :, None]-frequencies[:, None, :]))
    physical_triangles = triangles*phase[:, None]*duration**2
    # The other half-plane is the Hermitian partner, not a numerical PSD
    # symmetrization. Equal times have zero measure for these ordinary pulses.
    raw = physical_triangles + physical_triangles.conj().swapaxes(-1, -2)
    mean_outer = mean[:, :, None]*mean[:, None, :].conj()
    return {
        "greater": raw[:, 0]-mean_outer,
        "lesser": raw[:, 1]-mean_outer,
        "raw_greater": raw[:, 0],
        "raw_lesser": raw[:, 1],
        "mean_pulse": mean,
        "exit_state": exit_density.reshape(n, n),
        "residence_time_s": duration,
        "qrt_block_dimension": dimension,
        "qrt_exponentials": segments*len(frequencies),
        "density_exponentials": segments,
        "scope": "independent full-density segmented QRT; raw moments then global mean subtraction",
    }
