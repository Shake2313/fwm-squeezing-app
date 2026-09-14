"""Finite-path source noise from backward observables and Lindblad products.

Independent of forward atomic drift/diffusion/covariance lifts. The same physical
Hamiltonian and explicitly named jump operators define the reference problem.
"""

import time

import numpy as np
from scipy.integrate import solve_ivp

from gabes.quantum.reservoirs import ExplicitReservoirs
from .smooth_qrt import _positive_scalar, _real_array


def adjoint_wavepacket(h0, h1, envelope, reservoirs, boundary_state, duration_s,
                      readouts, frequencies_rad_s, *, drives=None,
                      drive_frequencies_rad_s=None, rtol=1e-10, atol=1e-12,
                      max_step_s=None):
    """Backward-observable source decomposition of Y_j=int exp(i*w_j*t)O_j dt.

    Return both connected orderings, each named source, mean, response, exit
    density. Physical jumps are fixed; H(t)=h0+envelope(t)*h1 may vary smoothly.
    The response drive is V_k exp(-i*v_k*t), including complex Nambu columns.
    All RF rows share ONE independently solved forward density. Backward
    terminal value is zero, not a final-state stationarity assumption.
    """
    if not isinstance(reservoirs, ExplicitReservoirs) or not callable(envelope):
        raise TypeError('explicit reservoirs and real envelope callback required')
    n = reservoirs.n_levels
    hs = np.asarray([h0, h1], complex)
    if (hs.shape != (2, n, n) or not np.isfinite(hs).all()
            or any(np.linalg.norm(h-h.conj().T) > 1e-12*max(np.linalg.norm(h), 1e-300) for h in hs)):
        raise ValueError('two finite Hermitian Hamiltonians required')
    rho0 = np.asarray(boundary_state, complex)
    if (rho0.shape != (n, n) or not np.isfinite(rho0).all()
            or abs(np.trace(rho0)-1) > 1e-10 or np.linalg.norm(rho0-rho0.conj().T) > 1e-10
            or np.linalg.eigvalsh((rho0+rho0.conj().T)/2).min() < -1e-10):
        raise ValueError('physical trace-one entry density required')
    ops = np.asarray(readouts, complex)
    freq = _real_array(frequencies_rad_s, 'readout frequencies')
    if (ops.ndim != 3 or ops.shape[1:] != (n, n) or not len(ops) or not np.isfinite(ops).all()
            or freq.ndim != 2 or freq.shape[1] != len(ops) or not len(freq)):
        raise ValueError('fixed nonempty readouts and matched frequency rows required')
    duration = _positive_scalar(duration_s, 'duration_s')
    rtol, atol = _positive_scalar(rtol, 'rtol'), _positive_scalar(atol, 'atol')
    if max(rtol, atol) >= 1:
        raise ValueError('tolerances must be below one')
    step = duration/64 if max_step_s is None else _positive_scalar(max_step_s, 'max_step_s')
    nf, ports = freq.shape
    if drives is None:
        if drive_frequencies_rad_s is not None:
            raise ValueError('drive frequencies require drives')
        v, vf = np.zeros((0, n, n), complex), np.zeros((nf, 0))
    else:
        v = np.asarray(drives, complex)
        vf = _real_array(drive_frequencies_rad_s, 'drive frequencies')
        if (v.ndim != 3 or v.shape[1:] != (n, n) or not len(v)
                or not np.isfinite(v).all() or vf.shape != (nf, len(v))):
            raise ValueError('fixed nonempty drives and matched frequency rows required')
    nr, inputs, nj = n*n, len(v), len(reservoirs.channels)
    eye = np.eye(n)
    l0 = duration*reservoirs.generator(hs[0])
    l1 = -1j*duration*(np.kron(hs[1], eye)-np.kron(eye, hs[1].T))
    omega, nu = duration*freq, duration*vf
    eta = 1+abs(omega)
    if not all(np.isfinite(x).all() for x in (l0, l1, omega, nu)):
        raise ValueError('dimensionless generators and frequencies must be finite')

    def amplitude(u):
        value = _real_array(envelope(float(u)*duration), 'envelope')
        if value.ndim != 0:
            raise ValueError('envelope must return one real scalar')
        return float(value)

    started = time.monotonic()
    # rho(t) is RF/source independent. Solve once; interpolation is separately
    # refined with the same ODE budget. Never invert a dissipative propagator.
    density = solve_ivp(lambda u, y: l0@y+amplitude(u)*(l1@y), (0., 1.),
                        rho0.ravel(), method='DOP853', rtol=rtol, atol=atol,
                        max_step=step/duration, dense_output=True)
    if not density.success:
        raise ValueError('forward density failed: '+density.message)
    ksize = nf*ports*nr
    qshape = (2, nj, nf, ports, ports)
    qend = ksize+int(np.prod(qshape))
    size = qend+nf*ports*inputs
    jumps = np.asarray([c.operator for c in reservoirs.channels]).reshape(nj, n, n)
    jump_dag = jumps.conj().swapaxes(-1, -2)
    forcing = eta[:, :, None]*ops.reshape(1, ports, nr)
    noise_frequency = omega[:, :, None]-omega[:, None, :]
    response_frequency = omega[:, :, None]-nu[:, None, :]

    def rhs(u, state):
        k = state[:ksize].reshape(nf, ports, n, n)
        q = state[ksize:qend].reshape(qshape)
        response = state[qend:].reshape(nf, ports, inputs)
        rho = density.sol(u).reshape(n, n)
        flat = k.reshape(nf*ports, nr)
        # Hilbert-Schmidt adjoint acts on column vec(K); rows multiply L.conj().
        kd = -(flat@l0.conj()+amplitude(u)*(flat@l1.conj())).reshape(nf, ports, nr)
        kd -= 1j*omega[:, :, None]*flat.reshape(nf, ports, nr)+forcing
        # Lindblad product identity, evaluated directly on backward operators:
        # Gamma_r(X,Y)=[L_r†,X][Y,L_r]. No forward D, A or C is used.
        # Greater: Gamma(K_j,K_k†); lesser: Gamma(K_k†,K_j).
        left = jump_dag[:, None, None]@k[None]-k[None]@jump_dag[:, None, None]
        right = k[None]@jumps[:, None, None]-jumps[:, None, None]@k[None]
        greater = np.einsum('rfjab,rfkcb,ca->rfjk', left, left.conj(), rho, optimize=False)
        lesser = np.einsum('rfkba,rfjbc,ca->rfjk', right.conj(), right, rho, optimize=False)
        qd = -duration*np.stack([greater, lesser])-1j*noise_frequency[None, None]*q
        # Exact phase/units changes: K=eta*exp(-i*w*t)*A/T. Q and response
        # likewise carry their remaining phase. Undo ALL factors at u=0.
        comm = k[:, :, None]@v[None, None]-v[None, None]@k[:, :, None]
        susceptibility = -1j*np.einsum('fjkab,ba->fjk', comm, rho)
        rd = -susceptibility-1j*response_frequency*response
        return np.concatenate((kd.ravel(), qd.ravel(), rd.ravel()))

    backward = solve_ivp(rhs, (1., 0.), np.zeros(size, complex), method='DOP853',
                         rtol=rtol, atol=atol, max_step=step/duration, t_eval=[0.])
    if not backward.success:
        raise ValueError('backward observables failed: '+backward.message)
    final = backward.y[:, -1]
    a = final[:ksize].reshape(nf, ports, n, n)*duration/eta[:, :, None, None]
    mean = np.einsum('fjab,ba->fj', a, rho0)
    mean_outer = mean[:, :, None]*mean[:, None, :].conj()
    greater_in = np.einsum('fjab,fkcb,ca->fjk', a, a.conj(), rho0)-mean_outer
    lesser_in = np.einsum('fkba,fjbc,ca->fjk', a.conj(), a, rho0)-mean_outer
    bath = final[ksize:qend].reshape(qshape)*duration**2/eta[None, None, :, :, None]/eta[None, None, :, None, :]
    greater = np.concatenate((greater_in[None], bath[0]))
    lesser = np.concatenate((lesser_in[None], bath[1]))
    retarded = final[qend:].reshape(nf, ports, inputs)*duration**2/eta[:, :, None]
    if not all(np.isfinite(x).all() for x in (greater, lesser, mean, retarded)):
        raise ValueError('nonfinite adjoint wavepacket')
    return {'greater': greater.sum(axis=0), 'lesser': lesser.sum(axis=0),
            'greater_by_source': greater, 'lesser_by_source': lesser,
            'mean_pulse': mean, 'mean_outer': mean_outer,
            'retarded_response': retarded, 'exit_state': density.y[:, -1].reshape(n, n),
            'source_names': ('atomic_inflow',)+tuple('jump:'+c.name for c in reservoirs.channels),
            'frequencies_rad_s': freq, 'residence_time_s': duration,
            'numerics': {'rtol': rtol, 'atol': atol, 'max_step_s': step,
                         'forward_density_evaluations': density.nfev,
                         'backward_evaluations': backward.nfev,
                         'density_mesh_points': len(density.t),
                         'backward_complex_variables': size,
                         'elapsed_seconds': time.monotonic()-started},
            'scope': 'independent backward-observable atomic source/response reference; no Maxwell field'}
