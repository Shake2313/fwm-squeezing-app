"""Independent full-density QRT for H(t) = h0 + envelope(t) * h1.

Raw operator sources, positive-time triangles and means are evolved together.
No atomic drift A, microscopic diffusion D or covariance propagator is used.
"""

import numpy as np
from scipy.integrate import solve_ivp
from scipy.sparse import coo_matrix

from gabes.quantum.reservoirs import ExplicitReservoirs


def _real_array(value, name):
    raw = np.asarray(value)
    if np.iscomplexobj(raw) and np.any(raw.imag != 0):
        raise ValueError(f"{name} must be real")
    result = np.asarray(raw.real, float)
    if not np.isfinite(result).all():
        raise ValueError(f"{name} must be finite")
    return result


def _positive_scalar(value, name):
    result = _real_array(value, name)
    if result.ndim != 0 or result <= 0:
        raise ValueError(f"{name} must be a strictly positive scalar")
    return float(result)


def _sparse_builder(dimension):
    """Accumulate small dense blocks without allocating a dense lifted matrix."""
    rows, columns, entries = [], [], []

    def add(row_indices, column_indices, block):
        row, column = np.nonzero(block)
        rows.append(np.asarray(row_indices)[row])
        columns.append(np.asarray(column_indices)[column])
        entries.append(block[row, column])

    def finish():
        matrix = coo_matrix((np.concatenate(entries),
                             (np.concatenate(rows), np.concatenate(columns))),
                            shape=(dimension, dimension)).tocsr()
        matrix.eliminate_zeros()
        return matrix

    return add, finish


def smooth_qrt(h0, h1, envelope, reservoirs, boundary_state, duration_s,
               readouts, frequencies_rad_s, *, rtol=2e-11, atol=2e-14,
               max_step_s=None):
    """Integrate Y_j = integral_0^T exp(i*w_j*t) O_j(t) dt by raw QRT.

    h0/h1[n,n] are Hermitian in rad/s. The real scalar envelope callback
    receives physical seconds. Reservoirs are fixed ExplicitReservoirs; the
    entry density[n,n] is physical and trace one, but need not be stationary.
    Fixed complex readouts[p,n,n] and real frequencies[nf,p] allow unequal
    signed frequencies within each row. Rows do not correlate with each other.

    Returns connected greater/lesser and raw_greater/raw_lesser[nf,p,p],
    mean_pulse[nf,p], exit_state[n,n], residence_time_s and solver diagnostics.
    Means have operator units times s, moments operator-product units times
    s**2. Nothing is clipped, trace-normalized or repaired at the output.

    DOP853 evolves dimensionless u=t/T and normalized pulse coordinates;
    rtol/atol apply to these coordinates, not directly to SI pulse moments.
    max_step_s defaults to T/64. Set it to resolve narrow envelope features
    and independently refine tolerances for tiny or rapidly rotating outputs.
    Only the endpoint is retained. All frequency rows share one density.
    """
    if not isinstance(reservoirs, ExplicitReservoirs):
        raise TypeError("one shared ExplicitReservoirs instance required")
    if not callable(envelope):
        raise TypeError("envelope must be callable with time in seconds")
    n = reservoirs.n_levels
    h0, h1 = np.asarray(h0, complex), np.asarray(h1, complex)
    rho, ops = np.asarray(boundary_state, complex), np.asarray(readouts, complex)
    for name, value in (("h0", h0), ("h1", h1)):
        if value.shape != (n, n):
            raise ValueError(f"{name} must have shape (n,n)")
        if not np.isfinite(value).all():
            raise ValueError(f"{name} must be finite")
        scale = max(float(np.linalg.norm(value)), np.finfo(float).tiny)
        if np.linalg.norm(value-value.conj().T) > 1e-12*scale:
            raise ValueError(f"{name} must be Hermitian")
    if ops.ndim != 3 or ops.shape[1:] != (n, n) or not len(ops):
        raise ValueError("readouts must have shape (p,n,n) with p > 0")
    if rho.shape != (n, n):
        raise ValueError("boundary_state must have shape (n,n)")
    if not np.isfinite(ops).all() or not np.isfinite(rho).all():
        raise ValueError("readouts and boundary state must be finite")
    if (abs(np.trace(rho)-1) > 1e-10
            or np.linalg.norm(rho-rho.conj().T) > 1e-10
            or np.linalg.eigvalsh((rho+rho.conj().T)/2)[0] < -1e-10):
        raise ValueError("physical trace-one boundary state required")
    duration = _positive_scalar(duration_s, "duration_s")
    rtol = _positive_scalar(rtol, "rtol")
    atol = _positive_scalar(atol, "atol")
    if rtol >= 1:
        raise ValueError("rtol must be less than one")
    max_step = (1/64 if max_step_s is None
                else _positive_scalar(max_step_s, "max_step_s")/duration)
    frequencies = _real_array(frequencies_rad_s, "frequencies_rad_s")
    ports = len(ops)
    if frequencies.ndim != 2 or frequencies.shape[1] != ports or not len(frequencies):
        raise ValueError("frequencies_rad_s must have shape (nf,p) with nf > 0")

    # Row-major vec(A rho B) = (A tensor B.T) vec(rho). Construct L1
    # directly: subtracting two dissipative generators can lose a small drive.
    eye = np.eye(n)
    l0 = duration*reservoirs.generator(h0)
    l1 = -1j*duration*(np.kron(h1, eye)-np.kron(eye, h1.T))
    omega = duration*frequencies
    if not all(np.isfinite(value).all() for value in (l0, l1, omega)):
        raise ValueError("dimensionless generators and frequencies must be finite")
    n2, nf = n*n, len(frequencies)
    x_count, z_count = 2*ports*n2, 2*ports*ports
    auxiliary_size = x_count+z_count+ports
    dimension = n2+nf*auxiliary_size
    density_indices = np.arange(n2)
    offsets = n2+np.arange(nf)*auxiliary_size
    x_indices = (offsets[:, None]+np.arange(x_count)).reshape(nf, 2, ports, n2)
    z_indices = (offsets[:, None]+x_count+np.arange(z_count)).reshape(nf, 2, ports, ports)
    mean_indices = offsets[:, None]+x_count+z_count+np.arange(ports)
    trace_rows = ops.swapaxes(-1, -2).reshape(ports, n2)
    sources = [(np.kron(op.conj().T, eye), np.kron(eye, op.conj())) for op in ops]
    add0, finish0 = _sparse_builder(dimension)
    add1, finish1 = _sparse_builder(dimension)
    add0(density_indices, density_indices, l0)
    add1(density_indices, density_indices, l1)
    for f in range(nf):
        add0(mean_indices[f], density_indices, trace_rows)
        add0(mean_indices[f], mean_indices[f], np.diag(-1j*omega[f]))
        for ordering in range(2):
            for k in range(ports):
                xs = x_indices[f, ordering, k]
                zs = z_indices[f, ordering, :, k]
                add0(xs, density_indices, sources[k][ordering])
                add0(xs, xs, l0+1j*omega[f, k]*np.eye(n2))
                add1(xs, xs, l1)
                add0(zs, xs, trace_rows)
                add0(zs, zs, np.diag(1j*(omega[f, k]-omega[f])))
    g0, g1 = finish0(), finish1()
    initial = np.zeros(dimension, complex)
    initial[:n2] = rho.ravel()

    def rhs(u, state):
        value = _real_array(envelope(u*duration), "envelope(time_s)")
        if value.ndim != 0:
            raise ValueError("envelope(time_s) must be a real scalar")
        return g0@state+float(value)*(g1@state)

    ode = solve_ivp(rhs, (0., 1.), initial, method="DOP853", t_eval=[1.],
                    rtol=rtol, atol=atol, max_step=max_step, dense_output=False)
    if not ode.success:
        raise ValueError(f"independent smooth full-density QRT failed: {ode.message}")
    final = ode.y[:, -1]
    if not np.isfinite(final).all():
        raise ValueError("independent smooth full-density QRT returned nonfinite state")
    mean = duration*np.exp(1j*omega)*final[mean_indices]
    phase = np.exp(1j*(omega[:, :, None]-omega[:, None, :]))
    triangles = duration**2*phase[:, None]*final[z_indices]
    # This adjoint is the other integration half-plane, not a PSD repair.
    raw = triangles+triangles.conj().swapaxes(-1, -2)
    mean_outer = mean[:, :, None]*mean[:, None, :].conj()
    return {
        "greater": raw[:, 0]-mean_outer,
        "lesser": raw[:, 1]-mean_outer,
        "raw_greater": raw[:, 0],
        "raw_lesser": raw[:, 1],
        "mean_pulse": mean,
        "exit_state": final[:n2].reshape(n, n),
        "residence_time_s": duration,
        "evaluations": ode.nfev,
        "qrt_block_dimension": dimension,
        "affine_generator_nnz": (g0.nnz, g1.nnz),
        "stored_time_points": len(ode.t),
        "solver_method": "DOP853",
        "rtol": rtol,
        "atol": atol,
        "max_step_s": max_step*duration,
        "scope": "independent full-density smooth QRT; raw moments then global mean subtraction",
    }
