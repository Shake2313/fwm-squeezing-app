"""Two-phase convective atomic mean and microscopic ordered noise.

This is an extended phase-coordinate calculation, not a physical optical-mode
channel. No Floquet Brillouin zone exists for a general irrational frequency
ratio; different lattice labels must not be counted as independent photons.
"""

from dataclasses import dataclass
from itertools import product

import numpy as np

from .. import core
from .contracts import GeneratorFrequencyAxis, readonly_array
from .diffusion import traceless_hermitian_basis
from .periodic import _order, _psd
from .reservoirs import ExplicitReservoirs


def rectangle(orders):
    if len(orders) != 2:
        raise ValueError('two harmonic cutoffs required')
    a, b = [_order(x, 'harmonic cutoff') for x in orders]
    return np.array(list(product(range(-a, a+1), range(-b, b+1))), dtype=int)


def phase_series(labels, coefficients, phases):
    return np.einsum('pq,q...->p...', np.exp(-1j*np.asarray(phases)@np.asarray(labels).T), coefficients)


def lattice_matrix(labels, coefficients, output_labels):
    """Block (a,b) is coefficient_(a-b); neither conjugate nor transpose."""
    lookup = {tuple(q): c for q, c in zip(labels, coefficients)}
    size = coefficients.shape[-1]
    out = np.zeros((len(output_labels)*size,)*2, complex)
    for i, a in enumerate(output_labels):
        for j, b in enumerate(output_labels):
            value = lookup.get(tuple(a-b))
            if value is not None:
                out[i*size:(i+1)*size, j*size:(j+1)*size] = value
    return out


@dataclass(frozen=True)
class TorusAtomicNoise:
    frequencies_rad_s: np.ndarray
    operators: np.ndarray
    state_labels: np.ndarray
    state_coefficients: np.ndarray
    generator_labels: np.ndarray
    generator_coefficients: np.ndarray
    drift_coefficients: np.ndarray
    diffusion_by_reservoir: np.ndarray
    reservoirs: ExplicitReservoirs
    diagnostics: dict

    def __post_init__(self):
        for key in ('frequencies_rad_s', 'operators', 'state_labels', 'state_coefficients',
                    'generator_labels', 'generator_coefficients', 'drift_coefficients', 'diffusion_by_reservoir'):
            object.__setattr__(self, key, readonly_array(getattr(self, key),
                real=key in ('frequencies_rad_s', 'state_labels', 'generator_labels')))

    def at_phase(self, phases_rad):
        phases = readonly_array(phases_rad, real=True)
        if phases.ndim != 2 or phases.shape[1] != 2 or not len(phases):
            raise ValueError('nonempty rows of two phases required')
        derivative = (-1j*(self.state_labels@self.frequencies_rad_s))[:, None, None]*self.state_coefficients
        return {'state': phase_series(self.state_labels, self.state_coefficients, phases),
            'state_derivative': phase_series(self.state_labels, derivative, phases),
            'drift': phase_series(self.generator_labels, self.drift_coefficients, phases),
            'diffusion_by_reservoir': np.array([phase_series(self.state_labels, d, phases)
                for d in self.diffusion_by_reservoir])}

    def spectrum(self, axis, *, response_orders, output_labels):
        """Selected rows of R D R†, retaining each physical reservoir.

        The carrier/grating labels are auxiliary coordinates. A physical field
        projection and spatial integration are still required before readout.
        No frequency is folded, rounded to a rational ratio, or deduplicated.
        """
        if not isinstance(axis, GeneratorFrequencyAxis):
            raise TypeError('explicit generator-frequency axis required')
        labels = rectangle(response_orders)
        selected = readonly_array(output_labels, real=True)
        if (selected.ndim != 2 or selected.shape[1] != 2 or not len(selected)
                or len({tuple(q) for q in selected}) != len(selected)):
            raise ValueError('distinct output lattice labels required')
        lookup = {tuple(q): i for i, q in enumerate(labels)}
        if any(tuple(q) not in lookup for q in selected):
            raise ValueError('output label outside response lattice')
        m = len(self.operators)
        drift = lattice_matrix(self.generator_labels, self.drift_coefficients, labels)
        drift += np.kron(np.diag(1j*(labels@self.frequencies_rad_s)), np.eye(m))
        maximum_real = float(np.linalg.eigvals(drift).real.max())
        if maximum_real >= 0:
            raise ValueError('nondecaying finite response lattice')
        ids = np.concatenate([np.arange(lookup[tuple(q)]*m, (lookup[tuple(q)]+1)*m) for q in selected])
        eye = np.eye(len(drift))
        # EXACT selected resolvent rows: A.T X=E.T -> X.T=E A^-1.
        # This avoids solving unused right-hand sides, without dropping modes.
        resolvents = [np.linalg.solve((-1j*w*eye-drift).T, eye[ids].T).T for w in axis.omega_rad_s]
        greater, lesser = [], []
        for coefficients in self.diffusion_by_reservoir:
            pair = [lattice_matrix(self.state_labels, value, labels)
                    for value in (coefficients, coefficients.swapaxes(-1, -2))]
            for out, diffusion in zip((greater, lesser), pair):
                if not _psd(diffusion)['passed']:
                    raise ValueError('nonpositive finite lattice diffusion; refine mean, never clip')
                out.append([r@diffusion@r.conj().T for r in resolvents])
        greater, lesser = np.array(greater), np.array(lesser)
        audit = {'greater': _psd(greater), 'lesser': _psd(lesser),
            'maximum_finite_lift_real_eigenvalue_s_inverse': maximum_real}
        audit['passed'] = audit['greater']['passed'] and audit['lesser']['passed']
        if not audit['passed']:
            raise ValueError('ordered lattice spectrum failed PSD')
        return {'frequency_axis': axis, 'output_labels': readonly_array(selected, real=True),
            'greater_by_reservoir': readonly_array(greater), 'lesser_by_reservoir': readonly_array(lesser),
            'greater': readonly_array(greater.sum(axis=0)), 'lesser': readonly_array(lesser.sum(axis=0)),
            'reservoir_names': tuple(c.name for c in self.reservoirs.channels), 'audit': audit,
            'scope': 'phase-coordinate atomic kernel; not independent optical modes or a squeezing spectrum'}


def torus_atomic_noise(h0, drive_labels, drives, reservoirs, frequencies_rad_s, *, mean_orders,
                       phase_samples=None, consistency_rtol=1e-7):
    """H=H0+sum_a(V_a exp(-i a.theta)+h.c.), with theta_dot=frequencies.

    The two phase coordinates may include a static spatial grating. A zero
    convective frequency alone does not authorize quotienting a spatial phase.
    """
    if not isinstance(reservoirs, ExplicitReservoirs) or not reservoirs.channels:
        raise TypeError('nonempty explicit reservoirs required')
    h0, drives = readonly_array(h0), readonly_array(drives)
    frequencies, drive_labels = readonly_array(frequencies_rad_s, real=True), readonly_array(drive_labels, real=True)
    n = reservoirs.n_levels
    if (h0.shape != (n, n) or drives.ndim != 3 or drives.shape[1:] != (n, n)
            or frequencies.shape != (2,) or drive_labels.shape != (len(drives), 2)
            or np.any(drive_labels != np.round(drive_labels)) or np.any(np.all(drive_labels == 0, axis=1))
            or not np.isfinite(consistency_rtol) or not 0 < consistency_rtol <= 1e-5):
        raise ValueError('matched Hamiltonians, nonzero integer drive labels and finite two-phase frequencies required')
    labels = rectangle(mean_orders)
    if not all(np.any(np.all(labels == q, axis=1)) for q in drive_labels):
        raise ValueError('mean lattice must contain every drive label')
    generators = {(0, 0): reservoirs.generator(h0)}
    for q, drive in zip(drive_labels.astype(int), drives):
        for key, value in ((tuple(q), drive), (tuple(-q), drive.conj().T)):
            generators[key] = generators.get(key, np.zeros((n*n, n*n), complex))+core.comm_super(value)
    generator_labels = np.array(list(generators))
    generator_coefficients = np.array(list(generators.values()))
    operators = traceless_hermitian_basis(n)
    m = len(operators)
    basis = operators.reshape(m, n*n).T
    drift = np.array([basis.conj().T@g@basis for g in generator_coefficients])
    lifted = lattice_matrix(generator_labels, drift, labels)
    lifted += np.kron(np.diag(1j*(labels@frequencies)), np.eye(m))
    forcing = np.zeros((len(labels), m), complex)
    for i, q in enumerate(labels):
        generator = generators.get(tuple(q))
        if generator is not None:
            forcing[i] = basis.conj().T@generator@(np.eye(n).reshape(-1)/n)
    coordinates = np.linalg.solve(lifted, -forcing.reshape(-1)).reshape(len(labels), m)
    states = (coordinates@basis.T).reshape(len(labels), n, n)
    states[np.flatnonzero(np.all(labels == 0, axis=1))[0]] += np.eye(n)/n
    diffusion = []
    for channel in reservoirs.channels:
        comm = operators@channel.operator-channel.operator@operators
        products = comm.conj().swapaxes(-1, -2)[:, None]@comm[None, :]
        diffusion.append(np.einsum('ijab,qba->qij', products, states))
    diffusion = np.array(diffusion)
    provisional = TorusAtomicNoise(frequencies, operators, labels, states, generator_labels,
        generator_coefficients, drift, diffusion, reservoirs, {})
    counts = tuple(4*int(o)+3 for o in mean_orders) if phase_samples is None else tuple(phase_samples)
    if len(counts) != 2 or any(_order(c, 'phase samples') < 4*int(o)+3 for c, o in zip(counts, mean_orders)):
        raise ValueError('phase grid must resolve products of retained mean harmonics')
    phases = np.array(list(product(*(2*np.pi*np.arange(c)/c for c in counts))))
    sampled = provisional.at_phase(phases)
    rho, rho_dot, a, d = (sampled[k] for k in ('state', 'state_derivative', 'drift', 'diffusion_by_reservoir'))
    mu = np.einsum('iab,pba->pi', operators, rho)
    mu_dot = np.einsum('iab,pba->pi', operators, rho_dot)
    products = operators[:, None]@operators[None, :]
    covariance = np.einsum('ijab,pba->pij', products, rho)-mu[:, :, None]*mu[:, None, :]
    cdot = np.einsum('ijab,pba->pij', products, rho_dot)-mu_dot[:, :, None]*mu[:, None, :]-mu[:, :, None]*mu_dot[:, None, :]
    total_d = d.sum(axis=0)
    terms = [cdot, -a@covariance, -covariance@a.swapaxes(-1, -2), -total_d]
    scale = np.maximum(sum(np.linalg.norm(t, axis=(-2, -1)) for t in terms), np.finfo(float).tiny)
    k, kdot = covariance-covariance.swapaxes(-1, -2), cdot-cdot.swapaxes(-1, -2)
    defect_k = kdot-a@k-k@a.swapaxes(-1, -2)-total_d+total_d.swapaxes(-1, -2)
    derivative = (phase_series(generator_labels, generator_coefficients, phases)@rho.reshape(len(phases), n*n, 1))[..., 0]
    generator_scale = sum(np.linalg.norm(g) for g in generator_coefficients)
    diagnostics = {'phase_grid': list(counts),
        'maximum_state_equation_relative_residual': float(np.linalg.norm(derivative-rho_dot.reshape(len(phases), n*n), axis=-1).max()/generator_scale),
        'maximum_state_trace_error': float(abs(np.trace(rho, axis1=-2, axis2=-1)-1).max()),
        'maximum_state_hermiticity_error': float(np.linalg.norm(rho-rho.conj().swapaxes(-1, -2), axis=(-2, -1)).max()),
        'minimum_sampled_state_eigenvalue': float(np.linalg.eigvalsh((rho+rho.conj().swapaxes(-1, -2))/2).min()),
        'maximum_dynamic_einstein_relative_residual': float((np.linalg.norm(sum(terms), axis=(-2, -1))/scale).max()),
        'maximum_dynamic_commutator_relative_residual': float((np.linalg.norm(defect_k, axis=(-2, -1))/scale).max()),
        'phase_diffusion': _psd(d), 'phase_ordered_covariance': _psd(covariance),
        'scope': 'sampled local convective atom; separate cutoff/stability/reference checks required'}
    diagnostics['passed'] = bool(diagnostics['maximum_state_equation_relative_residual'] < consistency_rtol
        and diagnostics['maximum_state_trace_error'] < 1e-9 and diagnostics['maximum_state_hermiticity_error'] < 1e-9
        and diagnostics['minimum_sampled_state_eigenvalue'] >= -1e-10
        and diagnostics['maximum_dynamic_einstein_relative_residual'] < consistency_rtol
        and diagnostics['maximum_dynamic_commutator_relative_residual'] < consistency_rtol
        and diagnostics['phase_diffusion']['passed'] and diagnostics['phase_ordered_covariance']['passed'])
    if not diagnostics['passed']:
        raise ValueError(f'two-phase atomic consistency failed: {diagnostics}')
    return TorusAtomicNoise(frequencies, operators, labels, states, generator_labels,
        generator_coefficients, drift, diffusion, reservoirs, diagnostics)
