"""Exact spectral source convolutions and positive-substep CF4 propagation.

Frequency oscillations remain in matrix exponentials. Only the time-dependent
Hamiltonian integration is approximate. Source/order/response refinements are
required before numerical certification of a physical path.
"""

import time
import numpy as np
from scipy.linalg import block_diag, eig, expm

from gabes.quantum.contracts import readonly_array
from gabes.quantum.diffusion import _state, traceless_hermitian_basis
from gabes.quantum.reservoirs import ExplicitReservoirs
from gabes.quantum.transport import _moments
from gabes.quantum.segmented_transport import _ordered_audit


def _divided_exponential(left, right, dt):
    """Integral exp(left*(dt-s))*exp(right*s) ds, including equal eigenvalues."""
    left, right = np.broadcast_arrays(left, right)
    delta = (right-left)*dt
    small = abs(delta) < 1e-3
    result = np.empty_like(delta, dtype=complex)
    z = delta[small]
    # Removable singularity: dt*exp(left*dt)*phi_1(delta). No inverse at zero.
    phi = 1+z*(.5+z*(1/6+z*(1/24+z*(1/120+z/720))))
    result[small] = dt*np.exp(left[small]*dt)*phi
    # Difference of bounded exponentials avoids 0*inf for decaying modes.
    result[~small] = (np.exp(right[~small]*dt)-np.exp(left[~small]*dt))/(right[~small]-left[~small])
    return result


def _density_modes(l, dmap):
    lam, r = eig(l)
    condition = float(np.linalg.cond(r))
    residual = float(np.linalg.norm(l@r-r*lam[None])/max(np.linalg.norm(l)*np.linalg.norm(r), 1e-300))
    if not np.isfinite(condition) or condition > 1e6:
        return lam, r, None, None, condition, residual
    return lam, r, np.linalg.inv(r), dmap@r, condition, residual


def _drift_modes(aug):
    beta, s = eig(aug)
    condition = float(np.linalg.cond(s))
    residual = float(np.linalg.norm(aug@s-s*beta[None])/max(np.linalg.norm(aug)*np.linalg.norm(s), 1e-300))
    return beta, s, condition, residual


def _source_map(aug, l, dmap, dt, max_condition, density_modes=None, drift_modes=None):
    """Build a source map once, independent of the entering density/covariance."""
    dim, nr = len(aug), len(l)
    beta, s, cond_s, residual = _drift_modes(aug) if drift_modes is None else drift_modes
    lam, r, invr, d_modes, cond_r, rho_residual = _density_modes(l, dmap) if density_modes is None else density_modes
    condition = max(cond_s, cond_r)
    # Conditioning amplification involves BOTH covariance basis factors and
    # the density basis. This is a conservative numerical screen, not a global
    # error proof; fixed-map parity and path refinement remain mandatory.
    roundoff_scale = np.finfo(float).eps*cond_s*cond_s*cond_r
    if (invr is not None and np.isfinite(condition) and condition <= max_condition
            and roundoff_scale <= 1e-8 and max(residual, rho_residual) <= 1e-12):
        invs = np.linalg.inv(s)
        forcing = d_modes.reshape(len(dmap), dim, dim, nr)
        transformed = np.einsum('ai,sijb,cj->sacb', invs, forcing, invs.conj(), optimize=True)
        summed = beta[:, None]+beta.conj()[None, :]
        kernel = _divided_exponential(summed[:, :, None], lam[None, None, :], dt)
        if np.isfinite(kernel).all():
            homogeneous = (s*np.exp(beta*dt)[None])@invs
            return ('spectral', homogeneous, s, invr, transformed*kernel[None], condition)
    # A defective/ill-conditioned eigenbasis is not an approximation license.
    # Use the original inverse-free block exponential for THIS step instead.
    full = np.zeros((dim*dim+len(dmap)*nr,)*2, complex)
    full[:dim*dim, :dim*dim] = np.kron(aug, np.eye(dim))+np.kron(np.eye(dim), aug.conj())
    if len(dmap):
        full[:dim*dim, dim*dim:] = np.concatenate(dmap, axis=1)
        full[dim*dim:, dim*dim:] = block_diag(*([l]*len(dmap)))
    return ('block', expm(full*dt), None, None, None, condition)


def _apply_source_map(mapping, cov, rho):
    """Apply the same map to all named sources; no source-dependent atom solves."""
    kind, e, s, invr, weighted, _ = mapping
    ns, dim, nr = cov.shape[1], cov.shape[-1], rho.size
    if kind == 'spectral':
        result = e@cov@e.conj().T
        # Exact identity: vec(D_r(rho(s))) is linear in the SAME rho(s).
        # In the two eigenbases its convolution is scalar divided exponentials.
        # Preserve every r and ordering. Never replace with total noise or fit.
        addition = np.einsum('sijk,k->sij', weighted, invr@rho.ravel(), optimize=True)
        addition = s@addition@s.conj().T
        if ns > 1:
            result[:, 1:] += addition.reshape(2, ns-1, dim, dim)
        return result
    result = (e[:dim*dim, :dim*dim]@cov.reshape(2*ns, -1).T).T.reshape(cov.shape)
    for order in range(2):
        for source in range(ns-1):
            start = dim*dim+(order*(ns-1)+source)*nr
            result[order, source+1] += (e[:dim*dim, start:start+nr]@rho.ravel()).reshape(dim, dim)
    return result


def _first_moment_map(mapping, density_modes, drift_modes, trace_rows, drive_columns,
                      w, wf, dt):
    """Exact triangular forcing in the SAME two bases as the source map.

    x_k'=(B+i*nu_k*I)x_k+G_k*rho; rho'=L*rho. Its convolution is
    S [S^-1 G_k R * F(beta+i*nu_k, lambda)] R^-1, including resonances.
    Mean rows use B=-i*w. These are algebraic identities, not new solves,
    an adiabatic approximation or a fitted recentering of source noise.
    """
    lam, r, invr, _, _, _ = density_modes
    beta, s, _, _ = drift_modes
    if mapping[0] != 'spectral':
        return None
    invs = np.linalg.inv(s)
    mean_forcing = (trace_rows@r)*_divided_exponential(-1j*w[:, None], lam[None], dt)
    projected = invs@drive_columns@r
    forcing = projected*_divided_exponential(
        beta[None, :, None]+1j*wf[:, None, None], lam[None, None], dt)
    return (np.exp(-1j*w*dt), mean_forcing@invr,
            mapping[1][None]*np.exp(1j*wf*dt)[:, None, None], s@forcing@invr)


def _apply_first_moment_map(mapping, small, rho, nr, ports, dim):
    phase, mean_forcing, propagation, forcing = mapping
    result = small.copy()
    result[nr:nr+ports] = phase*small[nr:nr+ports]+mean_forcing@rho.ravel()
    x = small[nr+ports:].reshape(-1, dim)
    result[nr+ports:] = (np.einsum('kij,kj->ki', propagation, x)
                         + forcing@rho.ravel()).ravel()
    return result


def exponential_wavepacket(h0, h1, envelope, reservoirs, boundary_state, duration_s,
                          readouts, frequencies_rad_s, *, drives=None,
                          drive_frequencies_rad_s=None, segments=64, order=4,
                          max_condition=1e6, sample_count=17):
    """Same smooth-wavepacket input convention; segments controls macro steps.

    order=2: midpoint. order=4: two commutator-free exponentials per macro step.
    Both use exact fixed-generator source convolutions and complex responses.
    No Gaussian shortcut for density, source noise or carrier phases is made.
    """
    started = time.monotonic()
    if (isinstance(segments, bool) or not isinstance(segments, (int, np.integer)) or segments < 1
            or isinstance(order, bool) or order not in (2, 4)):
        raise ValueError('positive integer segments and order 2 or 4 required')
    duration = float(duration_s)
    if not np.isfinite(duration) or duration <= 0 or not callable(envelope):
        raise ValueError('positive finite duration and envelope callback required')
    if not np.isfinite(max_condition) or max_condition < 1:
        raise ValueError('finite eigenbasis condition limit >=1 required')
    if sample_count is not None and (isinstance(sample_count, bool)
            or not isinstance(sample_count, (int, np.integer)) or sample_count < 2):
        raise ValueError('integer sample_count >=2 or None required')
    hs = readonly_array([h0, h1])
    if hs.ndim != 3 or hs.shape[1] != hs.shape[2] or any(
            np.linalg.norm(h-h.conj().T) > 1e-12*max(np.linalg.norm(h), 1e-300) for h in hs):
        raise ValueError('two finite Hermitian Hamiltonians required')
    edges = np.linspace(0., duration, int(segments)+1)
    def value(t):
        f = np.asarray(envelope(float(t)))
        if f.ndim != 0 or np.iscomplexobj(f) or not np.isfinite(f):
            raise ValueError('envelope must return one finite real value')
        return float(f)
    values, durations = [], []
    for begin, end in zip(edges[:-1], edges[1:]):
        dt = end-begin
        if order == 2:
            values.append(value((begin+end)/2)); durations.append(dt)
        else:
            f1, f2 = value(begin+dt*(.5-np.sqrt(3)/6)), value(begin+dt*(.5+np.sqrt(3)/6))
            a1, a2 = (3-2*np.sqrt(3))/12, (3+2*np.sqrt(3))/12
            # CF4 for the ENTIRE affine moment lift G0+fG1. a1+a2=1/2:
            # each exponential has positive dt/2 and the unchanged jumps,
            # readout, drive and carrier terms. Real extrapolated H is legal;
            # never clip f_eff or rescale noise. Later map multiplies on left.
            values.extend((2*(a2*f1+a1*f2), 2*(a1*f1+a2*f2)))
            durations.extend((dt/2, dt/2))
    hamiltonians = hs[0]+np.asarray(values)[:, None, None]*hs[1]
    count = len(durations)
    ops = np.broadcast_to(readouts, (count,)+np.shape(readouts))
    v = None if drives is None else np.broadcast_to(drives, (count,)+np.shape(drives))
    stride = 2 if order == 4 else 1
    samples = np.unique(np.rint(np.linspace(0, segments, segments+1 if sample_count is None
                                             else min(segments+1, sample_count))).astype(int))
    result = _exponential_segments(hamiltonians, reservoirs, boundary_state, durations, ops,
                                  frequencies_rad_s, drives=v, drive_frequencies_rad_s=drive_frequencies_rad_s,
                                  max_condition=max_condition, history_steps=set(samples*stride))
    # CF4's first auxiliary half-step is NOT the physical midpoint's fourth-
    # order approximation. Expose histories only at completed macro steps.
    result['sample_ages_s'] = readonly_array(edges[samples], real=True)
    mean = result['mean_pulse']
    result['mean_outer'] = readonly_array(mean[:, :, None]*mean[:, None, :].conj())
    result['numerics'].update(segments=int(segments), order=int(order), substeps=count,
                             max_condition=max_condition, elapsed_seconds=time.monotonic()-started)
    result['numerics']['stored_time_points'] = len(samples)
    result['scope'] = 'CF time integration of prescribed smooth atomic path; convergence separate from quantum consistency'
    return result


def _exponential_segments(hamiltonians, reservoirs, boundary_state, durations_s,
                         readouts, frequencies_rad_s, *, drives=None,
                         drive_frequencies_rad_s=None, max_condition=1e6, history_steps=None):
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
    # Audit EVERY auxiliary step, retain only requested physical macro times.
    # Sum squared Frobenius norms equals the former stacked-history norm;
    # bounded memory changes storage only, never weakens a quantum check.
    norm_squared, error_squared, comm_squared = float(np.linalg.norm(cin)**2), 0., 0.
    minimum_state_eigenvalue = float(np.linalg.eigvalsh(state).min())
    # Mean and retarded response use an independent smaller linear lift.
    response_size = nr+ports+inputs*(m+ports)
    small = np.zeros((nf, response_size), complex)
    small[:, :nr] = rho0.ravel()
    response_exponentials, response_spectral, reused = 0, 0, 0
    spectral_maps, block_maps, max_observed_condition = 0, 0, 0.
    eta = 1+abs(freq)*duration
    previous = [None]*nf
    previous_density = None
    for segment, dt in enumerate(durations):
        l = reservoirs.generator(h[segment])
        a = basis.conj().T@l@basis
        if np.linalg.norm(a.imag) > 1e-11*max(np.linalg.norm(l), 1e-300):
            raise ValueError('Hermitian atomic drift must be real')
        a = a.real
        density_key = l.tobytes()
        # The density eigenproblem and all D_r density-mode coefficients are
        # mathematically RF-independent. Reuse exactly; do not repeat by RF.
        if previous_density is None or previous_density[0] != density_key:
            previous_density = (density_key, _density_modes(l, dmap))
        # rho propagation is mathematically independent of RF and sources.
        # Compute it once per substep; never restore duplicate density solves.
        rho_step = expm(l*dt)
        c = np.einsum('jab,iba->ji', ops[segment], f)
        for fi, w in enumerate(freq):
            key = (float(dt), h[segment].tobytes(), ops[segment].tobytes(), v[segment].tobytes())
            cached = previous[fi] if previous[fi] is not None and previous[fi][0] == key else None
            aug = np.zeros((dim, dim), complex)
            aug[:m, :m], aug[m:, :m], aug[m:, m:] = a, c*eta[fi,:,None]/duration, -1j*np.diag(w)
            drift_modes = _drift_modes(aug) if cached is None else None
            mapping = _source_map(aug, l, dmap, dt, max_condition, previous_density[1], drift_modes) if cached is None else cached[1]
            cov[:, :, fi] = _apply_source_map(mapping, cov[:, :, fi], state)
            if cached is None:
                spectral_maps += mapping[0] == 'spectral'
                block_maps += mapping[0] == 'block'
                max_observed_condition = max(max_observed_condition, mapping[-1])
            trace_rows = ops[segment].swapaxes(-1, -2).reshape(ports, nr)*eta[fi,:,None]/duration
            drive_columns = np.zeros((inputs, dim, nr), complex)
            for j in range(inputs):
                drive_columns[j, :m] = basis.conj().T@(-1j*(np.kron(v[segment,j], np.eye(n))-np.kron(np.eye(n), v[segment,j].T)))
            moments = (_first_moment_map(mapping, previous_density[1], drift_modes, trace_rows,
                       drive_columns, w, wf[fi], dt) if cached is None else cached[2])
            if moments is not None:
                small[fi] = _apply_first_moment_map(moments, small[fi], state, nr, ports, dim)
                response_spectral += cached is None
                reused += 2 if cached is not None else 0
                previous[fi] = (key, mapping, moments, None)
                continue
            b = np.zeros((response_size, response_size), complex)
            b[:nr, :nr] = l
            b[nr:nr+ports, :nr] = trace_rows
            b[nr:nr+ports, nr:nr+ports] = -1j*np.diag(w)
            for j in range(inputs):
                start = nr+ports+j*(m+ports)
                # Row-major vec: vec(V rho-rho V)=(V x I-I x V.T) vec(rho).
                b[start:start+m, :nr] = drive_columns[j, :m]
                b[start:start+m, start:start+m] = a+1j*wf[fi,j]*np.eye(m)
                b[start+m:start+m+ports, start:start+m] = c*eta[fi,:,None]/duration
                b[start+m:start+m+ports, start+m:start+m+ports] = -1j*np.diag(w-wf[fi,j])
            eb = expm(b*dt) if cached is None else cached[3]
            small[fi] = eb@small[fi]
            response_exponentials += cached is None
            reused += 0 if cached is None else 2
            previous[fi] = (key, mapping, None, eb)
        # All RF rows must be forced by the SAME entry density on the next
        # segment. A view into small[0] would change when the first RF advances.
        state = (rho_step@state.ravel()).reshape(n, n)
        small[:, :nr] = state.ravel()
        _state(state, n)
        exact = _moments(state, f)[1]
        evolved = cov[0, :, 0, :m, :m].sum(axis=0)
        norm_squared += float(np.linalg.norm(exact)**2)
        error_squared += float(np.linalg.norm(exact-evolved)**2)
        comm_squared += float(np.linalg.norm(exact-exact.T-evolved+evolved.T)**2)
        minimum_state_eigenvalue = min(minimum_state_eigenvalue, float(np.linalg.eigvalsh(state).min()))
        if history_steps is None or segment+1 in history_steps:
            states.append(state.copy())
            atomic_sources.append(cov[:, :, 0, :m, :m].copy())
    phase = np.exp(1j*freq*duration)/eta
    ordered = cov[:, :, :, m:, m:]*phase[None,None,:,:,None]*phase.conj()[None,None,:,None,:]*duration**2
    mean = small[:, nr:nr+ports]*phase*duration
    response = np.empty((nf, ports, inputs), complex)
    for j in range(inputs):
        start = nr+ports+j*(m+ports)+m
        response[:, :, j] = small[:, start:start+ports]*np.exp(1j*(freq-wf[:,j,None])*duration)*duration/eta
    states, atomic_sources = np.array(states), np.array(atomic_sources)
    scale = max(np.sqrt(norm_squared), 1e-300)
    # A dark ordering can be exactly zero. Retain its signed roundoff; testing
    # it against its own ~epsilon norm is meaningless. The absolute floor is
    # in s^2, from the finite-pulse operator bound, not an added noise term.
    floor = 64*np.finfo(float).eps*duration**2*max(float(np.linalg.norm(ops, axis=(-2,-1)).max())**2, 1e-300)
    audit = {
        'atomic_covariance_relative_residual': float(np.sqrt(error_squared)/scale),
        'atomic_commutator_relative_residual': float(np.sqrt(comm_squared)/scale),
        'minimum_state_eigenvalue': minimum_state_eigenvalue,
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
        'numerics': {'spectral_maps': spectral_maps, 'block_fallback_maps': block_maps,
            'maximum_eigenvector_condition': max_observed_condition if np.isfinite(max_observed_condition) else None,
            'cached_exponential_reuses': reused, 'small_response_exponentials': int(response_exponentials),
            'spectral_response_maps': int(response_spectral)},
        'scope': 'exact spectral/block maps of prescribed substeps; no optical field claim'}
