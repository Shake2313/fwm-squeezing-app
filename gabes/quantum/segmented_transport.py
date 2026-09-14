"""Exact finite-atom moments for a prescribed piecewise constant characteristic.

Fast internal/carrier phases are integrated by matrix exponentials, not sampled
away. Only freezing a smooth envelope into segments is an approximation. The
finite atom is not a canonical optical port or a Maxwell propagation channel.
"""

import numpy as np
from scipy.linalg import block_diag, expm

from .contracts import readonly_array
from .diffusion import _state, traceless_hermitian_basis
from .reservoirs import ExplicitReservoirs
from .transport import _moments


def _ordered_audit(value, absolute_floor):
    adj = value.conj().swapaxes(-1, -2)
    scale = np.linalg.norm(value, axis=(-2, -1))
    tolerance = absolute_floor+2e-7*scale
    minimum = np.linalg.eigvalsh((value+adj)/2)[..., 0]
    herm = np.linalg.norm(value-adj, axis=(-2, -1))
    return {'passed': bool(np.all(minimum >= -tolerance) and np.all(herm <= tolerance)),
        'minimum_eigenvalue': float(minimum.min()), 'absolute_roundoff_floor_s2': float(absolute_floor),
        'maximum_hermiticity_tolerance_ratio': float(np.max(herm/tolerance)),
        'minimum_eigenvalue_tolerance_ratio': float(np.min(minimum/tolerance))}


def segmented_wavepacket(hamiltonians, reservoirs, boundary_state, durations_s,
                         readouts, frequencies_rad_s, *, drives=None,
                         drive_frequencies_rad_s=None):
    """Y_j = integral exp(+i*w_j*t) O_j(t) dt, with ordered source covariance.

    H: (segment,n,n); O: (segment,port,n,n); w: (frequency,port).
    Optional V: (segment,input,n,n), w_in: (frequency,input), describes a
    perturbation V_k(t)*epsilon_k*exp(-i*w_in,k*t). Complex V columns are
    independent; supply their physical conjugates explicitly when needed.
    Readout and drive phases refer to one common age origin, never each segment.
    """
    if not isinstance(reservoirs, ExplicitReservoirs):
        raise TypeError('explicit reservoirs required')
    n = reservoirs.n_levels
    h, durations, ops, freq = (readonly_array(hamiltonians), readonly_array(durations_s, real=True),
                             readonly_array(readouts), readonly_array(frequencies_rad_s, real=True))
    if (durations.ndim != 1 or not len(durations) or np.any(durations <= 0)
            or h.shape != (len(durations), n, n)
            or ops.ndim != 4 or ops.shape[0] != len(durations) or ops.shape[2:] != (n, n)
            or not ops.shape[1] or freq.ndim != 2 or freq.shape[1] != ops.shape[1] or not len(freq)):
        raise ValueError('matched finite H, positive segment durations, readouts and port frequencies required')
    if np.linalg.norm(h-h.conj().swapaxes(-1, -2)) > 1e-12*max(np.linalg.norm(h), 1e-300):
        raise ValueError('Hermitian Hamiltonians required')
    rho0 = _state(boundary_state, n)
    duration = float(durations.sum())
    if not np.isfinite(duration):
        raise ValueError('finite total duration required')
    nf, ports = freq.shape
    if drives is None:
        if drive_frequencies_rad_s is not None:
            raise ValueError('drive frequencies require drive columns')
        v, wf = np.zeros((len(durations), 0, n, n), complex), np.zeros((nf, 0))
    else:
        v, wf = readonly_array(drives), readonly_array(drive_frequencies_rad_s, real=True)
        if (v.ndim != 4 or v.shape[0] != len(durations) or v.shape[2:] != (n, n)
                or not v.shape[1] or wf.shape != (nf, v.shape[1])):
            raise ValueError('matched drive columns and frequencies required')
    inputs = v.shape[1]
    f = traceless_hermitian_basis(n)
    basis, m = f.reshape(-1, n*n).T, len(f)
    dim, nr = m+ports, n*n
    ns = 1+len(reservoirs.channels)
    # Independently derived microscopic D_r(rho) is LINEAR in rho. Cache
    # [Lr†,Fi][Fj,Lr], not a commutator-restoring residual or fitted coefficient.
    products = []
    for channel in reservoirs.channels:
        comm = f@channel.operator-channel.operator@f
        products.append(np.einsum('iab,jbc->ijac', comm.conj().swapaxes(-1, -2), comm))
    products = np.asarray(products).reshape(ns-1, m, m, n, n)
    _, cin = _moments(rho0, f)
    cov = np.zeros((2, ns, nf, dim, dim), complex)
    cov[0, 0, :, :m, :m], cov[1, 0, :, :m, :m] = cin, cin.T
    # Each independent source/order is an affine forcing by the SAME actual
    # rho. A block exponential computes all convolutions, even if L or the
    # Kronecker sum is singular. No inverse of a decaying propagator is used.
    forcing = np.zeros((2*(ns-1), dim, dim, nr), complex)
    for order in range(2):
        p = products if order == 0 else products.swapaxes(1, 2)
        forcing[order*(ns-1):(order+1)*(ns-1), :m, :m] = p.swapaxes(-1, -2).reshape(ns-1, m, m, nr)
    dmap = forcing.reshape(2*(ns-1), dim*dim, nr)
    state = np.array(rho0)
    states, atomic_sources = [state.copy()], [cov[:, :, 0, :m, :m].copy()]
    # Mean and retarded response use an independent smaller linear lift.
    response_size = nr+ports+inputs*(m+ports)
    small = np.zeros((nf, response_size), complex)
    small[:, :nr] = rho0.ravel()
    exponentials, reused = 0, 0
    previous = [None]*nf
    for segment, dt in enumerate(durations):
        l = reservoirs.generator(h[segment])
        a = basis.conj().T@l@basis
        if np.linalg.norm(a.imag) > 1e-11*max(np.linalg.norm(l), 1e-300):
            raise ValueError('Hermitian atomic drift must be real')
        a = a.real
        c = np.einsum('jab,iba->ji', ops[segment], f)
        for fi, w in enumerate(freq):
            key = (float(dt), h[segment].tobytes(), ops[segment].tobytes(), v[segment].tobytes())
            cached = previous[fi] if previous[fi] is not None and previous[fi][0] == key else None
            aug = np.zeros((dim, dim), complex)
            aug[:m, :m], aug[m:, :m], aug[m:, m:] = a, c/duration, -1j*np.diag(w)
            k = np.kron(aug, np.eye(dim))+np.kron(np.eye(dim), aug.conj())
            full = np.zeros((dim*dim+len(dmap)*nr,)*2, complex)
            full[:dim*dim, :dim*dim] = k
            if len(dmap):
                full[:dim*dim, dim*dim:] = np.concatenate(dmap, axis=1)
                full[dim*dim:, dim*dim:] = block_diag(*([l]*len(dmap)))
            # Mathematical equality, not a stationary-state assumption: with
            # identical H/O/V/dt the augmented linear map is identical for ANY
            # entering rho and source covariance. Reuse its exponential, but
            # apply it to the current state. Do not reinstate repeated solves.
            e = expm(full*dt) if cached is None else cached[1]
            cov[:, :, fi] = (e[:dim*dim, :dim*dim]@cov[:, :, fi].reshape(2*ns, -1).T).T.reshape(2, ns, dim, dim)
            for order in range(2):
                for source in range(ns-1):
                    start = dim*dim+(order*(ns-1)+source)*nr
                    cov[order, source+1, fi] += (e[:dim*dim, start:start+nr]@state.ravel()).reshape(dim, dim)
            b = np.zeros((response_size, response_size), complex)
            b[:nr, :nr] = l
            b[nr:nr+ports, :nr] = ops[segment].swapaxes(-1, -2).reshape(ports, nr)/duration
            b[nr:nr+ports, nr:nr+ports] = -1j*np.diag(w)
            for j in range(inputs):
                start = nr+ports+j*(m+ports)
                # Row-major vec: vec(V rho-rho V)=(V x I-I x V.T) vec(rho).
                bl = basis.conj().T@(-1j*(np.kron(v[segment,j], np.eye(n))-np.kron(np.eye(n), v[segment,j].T)))
                b[start:start+m, :nr] = bl
                b[start:start+m, start:start+m] = a+1j*wf[fi,j]*np.eye(m)
                b[start+m:start+m+ports, start:start+m] = c/duration
                b[start+m:start+m+ports, start+m:start+m+ports] = -1j*np.diag(w-wf[fi,j])
            eb = expm(b*dt) if cached is None else cached[2]
            small[fi] = eb@small[fi]
            exponentials += 2 if cached is None else 0
            reused += 0 if cached is None else 2
            previous[fi] = (key, e, eb)
        # All RF rows must be forced by the SAME entry density on the next
        # segment. A view into small[0] would change when the first RF advances.
        state = small[0, :nr].reshape(n, n).copy()
        _state(state, n)
        states.append(state.copy())
        atomic_sources.append(cov[:, :, 0, :m, :m].copy())
    phase = np.exp(1j*freq*duration)
    ordered = cov[:, :, :, m:, m:]*phase[None,None,:,:,None]*phase.conj()[None,None,:,None,:]*duration**2
    mean = small[:, nr:nr+ports]*phase*duration
    response = np.empty((nf, ports, inputs), complex)
    for j in range(inputs):
        start = nr+ports+j*(m+ports)+m
        response[:, :, j] = small[:, start:start+ports]*np.exp(1j*(freq-wf[:,j,None])*duration)*duration
    states, atomic_sources = np.array(states), np.array(atomic_sources)
    exact = np.array([_moments(rho, f)[1] for rho in states])
    evolved = atomic_sources[:, 0].sum(axis=1)
    scale = max(np.linalg.norm(exact), 1e-300)
    # A dark ordering can be exactly zero. Retain its signed roundoff; testing
    # it against its own ~epsilon norm is meaningless. The absolute floor is
    # in s^2, from the finite-pulse operator bound, not an added noise term.
    floor = 64*np.finfo(float).eps*duration**2*max(float(np.linalg.norm(ops, axis=(-2,-1)).max())**2, 1e-300)
    audit = {
        'atomic_covariance_relative_residual': float(np.linalg.norm(exact-evolved)/scale),
        'atomic_commutator_relative_residual': float(np.linalg.norm(exact-exact.swapaxes(-1,-2)-evolved+evolved.swapaxes(-1,-2))/scale),
        'minimum_state_eigenvalue': float(np.linalg.eigvalsh(states).min()),
        'greater': _ordered_audit(ordered[0], floor), 'lesser': _ordered_audit(ordered[1], floor)}
    audit['passed'] = bool(audit['atomic_covariance_relative_residual'] < 2e-7
        and audit['atomic_commutator_relative_residual'] < 2e-7
        and audit['minimum_state_eigenvalue'] >= -1e-10
        and audit['greater']['passed'] and audit['lesser']['passed'])
    if not audit['passed']:
        raise ValueError(f'segmented atomic covariance audit failed: {audit}')
    return {'greater_by_source': readonly_array(ordered[0]), 'lesser_by_source': readonly_array(ordered[1]),
        'greater': readonly_array(ordered[0].sum(axis=0)), 'lesser': readonly_array(ordered[1].sum(axis=0)),
        'mean_pulse': readonly_array(mean), 'retarded_response': readonly_array(response),
        'exit_state': readonly_array(state), 'states': readonly_array(states),
        'atomic_covariance_by_source': readonly_array(atomic_sources),
        'source_names': ('atomic_inflow',)+tuple('jump:'+r.name for r in reservoirs.channels),
        'frequencies_rad_s': freq, 'residence_time_s': duration, 'audit': audit,
        'matrix_exponentials': exponentials, 'identical_exponentials_reused': reused,
        'scope': 'exact prescribed piecewise constant finite atom; no optical SQL, Maxwell feedback or Gaussian-atom claim'}
