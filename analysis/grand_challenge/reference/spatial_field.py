"""Independent finite-aperture full-Liouville response and time-domain QRT.

Use original probe/conjugate plane-wave phases, never the production quotient
or profile-projection helper. The optical profiles enter the operators directly.
"""

import numpy as np

from gabes import core
from gabes.quantum.diffusion import traceless_hermitian_basis
from .periodic_field import forced_liouville_response
from .periodic_qrt import time_domain_qrt


def spatial_field_reference(h0, reservoirs, beat, lowering_operators, couplings, beta,
                            positions_m, area_weights_m2, mode_values_m_inverse,
                            reference_area_m2, relative_wavevectors_rad_m, z_m,
                            number_density_m3, frequencies_rad_s, *, mean_order=4,
                            response_order=3, qrt=True, phase_samples=16):
    n = reservoirs.n_levels
    atomic_ops = traceless_hermitian_basis(n)
    m = len(atomic_ops)
    frequencies = np.asarray(frequencies_rad_s)
    axis = np.unique(np.r_[frequencies, -frequencies])
    results = {name: np.zeros((len(frequencies), 4, 4), complex) for name in ('drift', 'greater', 'lesser')}
    mean_rate = np.zeros(2, complex)
    l0 = reservoirs.generator(h0)
    hs, signs, gs = np.array([1, -1, -1, 1]), np.array([1., 1., -1., -1.]), np.tile(couplings, 2)
    for xy, area, profiles in zip(positions_m, area_weights_m2, mode_values_m_inverse):
        position = np.r_[xy, z_m]
        chi = np.sqrt(reference_area_m2)*profiles*np.exp(1j*np.asarray(relative_wavevectors_rad_m)@position)
        physical_ops = chi.conj()[:, None, None]*lowering_operators
        ops = np.concatenate([physical_ops, physical_ops.conj().swapaxes(-1, -2)])
        v = (couplings[0]*beta[0]*physical_ops[0].conj().T
             +couplings[1]*np.conjugate(beta[1])*physical_ops[1])
        lp, lm = core.comm_super(v), core.comm_super(v.conj().T)
        states = core.floquet_solve_direct(l0, lp, lm, beat, [0.], np.zeros_like(l0), n,
            n_f=mean_order, return_harmonics=True)[0]
        density = number_density_m3*area
        results['drift'] += forced_liouville_response(l0, lp, lm, beat, states, ops, hs, gs,
            signs, frequencies, response_order=response_order, linear_density=density)
        for j in range(2):
            mean_rate[j] += -1j*density*couplings[j]*np.trace(physical_ops[j]@states[mean_order+int(hs[j])])
        if qrt:
            ref = time_domain_qrt(l0, lp, lm, beat, atomic_ops, axis, [-1, 1],
                phase_samples=phase_samples, rtol=2e-12, atol=2e-14)
            readout = np.zeros((4, 2*m), complex)
            for j, (h, op, g, sign) in enumerate(zip(hs, ops, gs, signs)):
                block = 0 if h == -1 else m
                readout[j, block:block+m] = -1j*sign*g*np.einsum('ab,kba->k', op, atomic_ops)
            flip = np.r_[np.arange(m, 2*m), np.arange(m)]
            for i, w in enumerate(frequencies):
                index, mirror = [int(np.flatnonzero(axis == s*w)[0]) for s in (1, -1)]
                dg = ref['ordered_spectrum'][index]
                dl = ref['ordered_spectrum'][mirror].T[flip][:, flip]
                results['greater'][i] += density*readout@dg@readout.conj().T
                results['lesser'][i] += density*readout@dl@readout.conj().T
    return {**results, 'mean_rate': mean_rate,
        'scope': 'stationary atomic centers; original carrier phases, full-density forcing and independent time QRT'}
