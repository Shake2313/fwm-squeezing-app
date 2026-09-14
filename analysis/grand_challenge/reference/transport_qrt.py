"""Full-density driven Liouville and two-time QRT along a finite path.

No complete-basis drift, jump-product D, or propagated covariance is consumed.
Integrate the positive-time QRT triangle and its adjoint; include the actual
boundary state in the independently evolved density matrix.
"""

import numpy as np
from scipy.integrate import solve_ivp


def mean_pulse(path, hamiltonian, reservoirs, boundary_state, readout, *, rtol=2e-12, atol=2e-14):
    """Direct full-density mean integral for independent finite differences."""
    n, duration = reservoirs.n_levels, path.residence_time_s
    ports = len(readout(0., path.position(0.)))
    def rhs(age, y):
        r = path.position(age)
        rho = y[:n*n].reshape(n, n)
        return np.r_[reservoirs.generator(hamiltonian(age, r))@rho.ravel(),
            np.einsum('jab,ba->j', readout(age, r), rho)/duration]
    ode = solve_ivp(rhs, (0.,duration), np.r_[np.asarray(boundary_state).ravel(),np.zeros(ports)],
        method='DOP853',rtol=rtol,atol=atol,max_step=duration/64)
    if not ode.success:
        raise ValueError('independent mean pulse failed')
    return ode.y[n*n:,-1]*duration


def characteristic_qrt(path, hamiltonian, reservoirs, boundary_state, readout, frequencies_rad_s,
                       *, drives=None, rtol=2e-12, atol=2e-14, intervals=32):
    n, duration = reservoirs.n_levels, path.residence_time_s
    w = np.asarray(frequencies_rad_s)
    ports = len(readout(0., path.position(0.)))
    inputs = 0 if drives is None else len(drives(0., path.position(0.)))
    # X_k(t)=integral_0^t U(t,s)[delta O_k(s)† rho(s)] exp(-i*w*s) ds/T.
    # Reversed operator ordering uses rho(s) delta O_k(s)† instead.
    shapes = [(n, n), (2, len(w), ports, n, n), (2, len(w), ports, ports),
              (len(w), ports), (len(w), inputs, n, n), (len(w), ports, inputs)]
    lengths = [int(np.prod(s)) for s in shapes]
    edges = np.r_[0, np.cumsum(lengths)]
    def unpack(y):
        return [y[a:b].reshape(s) for a, b, s in zip(edges[:-1], edges[1:], shapes)]
    y = np.r_[np.asarray(boundary_state).ravel(), np.zeros(sum(lengths[1:]), complex)]
    def rhs(age, y):
        rho, x, _triangle, _mean, response, _output = unpack(y)
        r = path.position(age)
        generator = reservoirs.generator(hamiltonian(age, r))
        ops = np.asarray(readout(age, r))
        means = np.einsum('jab,ba->j', ops, rho)
        centered_adjoint = ops.conj().swapaxes(-1, -2)-means.conj()[:, None, None]*np.eye(n)
        source = np.stack([centered_adjoint@rho, rho@centered_adjoint])
        xd = np.einsum('ab,ofkb->ofka', generator, x.reshape(2, len(w), ports, n*n)).reshape(x.shape)
        xd += np.exp(-1j*w*age)[None, :, None, None, None]*source[:, None]/duration
        triangle = np.exp(1j*w*age)[None, :, None, None]*np.einsum('jab,ofkba->ofjk', ops, x)/duration
        mean = np.exp(1j*w*age)[:, None]*means/duration
        v = np.zeros((0, n, n), complex) if drives is None else np.asarray(drives(age, r))
        rd = np.einsum('ab,fkb->fka', generator, response.reshape(len(w), inputs, n*n)).reshape(response.shape)
        rd += 1j*w[:, None, None, None]*response-1j*(v@rho-rho@v)
        output = np.einsum('jab,fkba->fjk', ops, response)/duration
        return np.concatenate([value.ravel() for value in
            ((generator@rho.ravel()).reshape(n,n), xd, triangle, mean, rd, output)])
    ode = solve_ivp(rhs, (0., duration), y, method='DOP853', rtol=rtol, atol=atol,
        max_step=duration/intervals, dense_output=True)
    if not ode.success:
        raise ValueError('independent finite-time QRT failed')
    rho, _x, triangles, mean, _response, response = unpack(ode.y[:, -1])
    full = (triangles+triangles.conj().swapaxes(-1, -2))*duration**2
    return {'greater': full[0], 'lesser': full[1], 'mean_pulse': mean*duration,
        'retarded_response': response*duration, 'exit_state': rho, 'evaluations': ode.nfev,
        'scope': 'independent full-density positive-time QRT plus Hermitian partner; no microscopic D used'}
