"""Microscopic Floquet elimination into explicitly retained optical bands.

Each physical band has one annihilator, its conjugate, a signed carrier
harmonic and an explicit atomic lowering/readout operator. Atomic harmonics
are internal response coordinates, not additional propagating optical ports.
"""

from dataclasses import dataclass

import numpy as np
from scipy.integrate import solve_ivp

from .channels import GaussianChannel
from .contracts import GeneratorFrequencyAxis, readonly_array
from .periodic import PeriodicAtomicNoise
from .traveling import LocalNambuGenerator, NambuTransfer, _same_coordinates


@dataclass(frozen=True)
class PeriodicFieldPorts:
    mode_labels: tuple[str, ...]
    carrier_harmonics: np.ndarray
    lowering_operators: np.ndarray
    coupling_s_inverse_sqrt_flux: np.ndarray

    def __post_init__(self):
        labels = tuple(self.mode_labels)
        h = readonly_array(self.carrier_harmonics, real=True)
        ops = readonly_array(self.lowering_operators)
        g = readonly_array(self.coupling_s_inverse_sqrt_flux, real=True)
        count = len(labels)
        if (not count or len(set(labels)) != count or any(not isinstance(s, str) or not s.strip() for s in labels)
                or h.shape != (count,) or np.any(h != np.round(h)) or g.shape != h.shape or np.any(g <= 0)
                or ops.ndim != 3 or ops.shape[0] != count or ops.shape[1] != ops.shape[2]):
            raise ValueError("unique physical ports, integer carrier harmonics, operators and positive couplings required")
        object.__setattr__(self, "mode_labels", labels)
        for key, value in (("carrier_harmonics", h), ("lowering_operators", ops), ("coupling_s_inverse_sqrt_flux", g)):
            object.__setattr__(self, key, value)

    def nambu_coordinates(self):
        return (np.r_[self.carrier_harmonics, -self.carrier_harmonics],
                np.concatenate([self.lowering_operators, self.lowering_operators.conj().swapaxes(-1, -2)]),
                np.tile(self.coupling_s_inverse_sqrt_flux, 2),
                np.repeat([1., -1.], len(self.mode_labels)),
                tuple(s+suffix for suffix in (":a", ":adag") for s in self.mode_labels))


def periodic_field_couplings(atom, ports, *, response_order):
    """Independent density-operator driving B and Maxwell readout C0.

    B_(h,i)=-i*g_i*Tr(F [O_i†,rho_(h-h_i)]); C0 reads atomic harmonic h_i.
    Neither B nor C0 is inferred from a commutator defect or noise covariance.
    """
    if not isinstance(atom, PeriodicAtomicNoise) or not isinstance(ports, PeriodicFieldPorts):
        raise TypeError("validated periodic atom and explicit optical ports required")
    if ports.lowering_operators.shape[1:] != atom.state_harmonics.shape[1:]:
        raise ValueError("atomic dimensions must match")
    lift = atom.lift(response_order)
    hs, ops, g, signs, labels = ports.nambu_coordinates()
    if np.any(abs(hs) > response_order):
        raise ValueError("response ladder must contain every optical carrier harmonic")
    m, count = len(atom.operators), len(signs)
    b = np.zeros((len(lift.drift), count), complex)
    c0 = np.zeros((count, len(lift.drift)), complex)
    mean_order = len(atom.state_harmonics)//2
    for i, (hi, op, coupling, sign) in enumerate(zip(hs, ops, g, signs)):
        row = int(hi+response_order)*m
        c0[i, row:row+m] = -1j*sign*coupling*np.einsum('ab,kba->k', op, atom.operators)
        for j, h in enumerate(lift.harmonics):
            q = int(h-hi)
            if abs(q) <= mean_order:
                rho = atom.state_harmonics[q+mean_order]
                drive = -1j*coupling*(op.conj().T@rho-rho@op.conj().T)
                b[j*m:(j+1)*m, i] = np.einsum('kab,ba->k', atom.operators, drive)
    return lift, b, c0, signs, labels


def eliminate_periodic_atom(atom, ports, frequency_axis, *, response_order, linear_density_m_inverse):
    """Return photon-flux local M and both ordered D, all in m^-1."""
    if not isinstance(atom, PeriodicAtomicNoise) or not isinstance(ports, PeriodicFieldPorts):
        raise TypeError('periodic atom and physical ports required')
    if not isinstance(frequency_axis, GeneratorFrequencyAxis):
        raise TypeError("explicit Floquet quasifrequency axis required")
    density = float(linear_density_m_inverse)
    if not np.isfinite(density) or density < 0:
        raise ValueError("finite nonnegative linear density required")
    if np.any(abs(frequency_axis.omega_rad_s) >= abs(atom.beat_rad_s)/2):
        raise ValueError("quasifrequencies must be strictly inside one Floquet zone")
    lift, b, c0, signs, labels = periodic_field_couplings(atom, ports, response_order=response_order)
    drift, greater, lesser = [], [], []
    eye = np.eye(len(lift.drift))
    for w in frequency_axis.omega_rad_s:
        cr = np.linalg.solve((-1j*w*eye-lift.drift).T, c0.T).T
        drift.append(density*cr@b)
        greater.append(density*(cr@lift.greater_by_reservoir@cr.conj().T))
        lesser.append(density*(cr@lift.lesser_by_reservoir@cr.conj().T))
    result = LocalNambuGenerator(frequency_axis, signs, labels, drift,
        np.swapaxes(greater, 0, 1), np.swapaxes(lesser, 0, 1), lift.reservoir_names)
    if not result.audit()['passed']:
        raise ValueError(f'periodic field commutator/PSD failed; no noise repair: {result.audit()}')
    return result


def paired_frequency_channel(transfer, mode_labels, index, mirror):
    """Full Nambu at +/-RF -> 2*n physical sideband Gaussian channel.

    Ordering is (all bands:+, all bands:-), quadratures (x,p) per mode.
    Retains inter-sector normal and anomalous correlations. Point samples mean
    the narrow normalized spectral-mode limit; DC is not two independent modes.
    """
    if not isinstance(transfer, NambuTransfer) or not transfer.audit()['passed']:
        raise ValueError('valid Nambu transfer required')
    labels = tuple(mode_labels)
    n = len(labels)
    if (not np.array_equal(transfer.signs, np.repeat([1., -1.], n))
            or transfer.mode_labels != tuple(s+suffix for suffix in (':a', ':adag') for s in labels)):
        raise ValueError('complete annihilator/creator ordering and matched labels required')
    w = transfer.frequency_axis.omega_rad_s
    if w[index] <= 0 or not np.isclose(w[index], -w[mirror], rtol=1e-13, atol=1e-12):
        raise ValueError('strictly positive RF and its unique negative partner required')
    flip = np.r_[np.arange(n, 2*n), np.arange(n)]
    for a, b in ((transfer.transfer, transfer.transfer), (transfer.noise_greater, transfer.noise_lesser)):
        left, right = a[index], b[mirror].conj()[flip][:, flip]
        if np.linalg.norm(left-right) > 1e-8*max(np.linalg.norm(left), np.linalg.norm(right), np.finfo(float).tiny):
            raise ValueError('full Nambu frequency reflection failed')
    tp, tm = transfer.transfer[index], transfer.transfer[mirror]
    u, v = np.zeros((2*n, 2*n), complex), np.zeros((2*n, 2*n), complex)
    u[:n, :n], u[n:, n:] = tp[:n, :n], tm[:n, :n]
    v[:n, n:], v[n:, :n] = tp[:n, n:], tm[:n, n:]
    x = np.zeros((4*n, 4*n))
    x[::2, ::2], x[::2, 1::2] = (u+v).real, -(u-v).imag
    x[1::2, ::2], x[1::2, 1::2] = (u+v).imag, (u-v).real
    np_, nm = [(transfer.noise_greater[i]+transfer.noise_lesser[i])/2 for i in (index, mirror)]
    k, anomalous = np.zeros_like(u), np.zeros_like(u)
    k[:n, :n], k[n:, n:] = np_[:n, :n], nm[:n, :n]
    anomalous[:n, n:], anomalous[n:, :n] = np_[:n, n:], nm[:n, n:]
    y = np.zeros_like(x)
    y[::2, ::2], y[1::2, 1::2] = (k+anomalous).real, (k-anomalous).real
    y[::2, 1::2], y[1::2, ::2] = (anomalous-k).imag, (anomalous+k).imag
    modes = tuple(s+suffix for suffix in (':+', ':-') for s in labels)
    result = GaussianChannel(x, y, modes, modes, 'full Floquet Nambu paired narrowband modes')
    if not result.audit().passed:
        raise ValueError('full paired-frequency quantum channel failed CP')
    return result


def adaptive_field_propagation(local_at_z, length_m, *, rtol=2e-10, atol=2e-12):
    """Integrate dT=M T and dN=M N+N M†+D along a prescribed mean path.

    Full complex ordered covariances are integrated separately; positivity and
    the global field commutator are checked, never projected or repaired.
    """
    length = float(length_m)
    if not np.isfinite([length, rtol, atol]).all() or length < 0 or rtol <= 0 or atol <= 0:
        raise ValueError('finite length>=0 and positive solver tolerances required')
    first = local_at_z(0.)
    if not isinstance(first, LocalNambuGenerator) or not first.audit()['passed']:
        raise ValueError('validated local Nambu generators required')
    nf, n = len(first.frequency_axis.omega_rad_s), len(first.signs)
    initial = np.zeros((3, nf, n, n), complex)
    initial[0] = np.eye(n)
    max_local_residual = first.audit()['maximum_commutator_relative_residual']
    def rhs(z, raw):
        nonlocal max_local_residual
        local = local_at_z(z)
        _same_coordinates(first, local)
        audit = local.audit()
        if not audit['passed']:
            raise ValueError('invalid local generator along propagation')
        max_local_residual = max(max_local_residual, audit['maximum_commutator_relative_residual'])
        t, ng, nl = raw.reshape(initial.shape)
        m = local.drift
        md = m.conj().swapaxes(-1, -2)
        return np.array([m@t, m@ng+ng@md+local.noise_greater, m@nl+nl@md+local.noise_lesser]).reshape(-1)
    ode = solve_ivp(rhs, (0., length), initial.reshape(-1), method='DOP853', rtol=rtol, atol=atol,
                    max_step=length/4 if length else np.inf)
    if not ode.success:
        raise ValueError('adaptive field propagation failed')
    t, ng, nl = ode.y[:, -1].reshape(initial.shape)
    result = NambuTransfer(first.frequency_axis, first.signs, first.mode_labels, t, ng, nl)
    if not result.audit()['passed']:
        raise ValueError('adaptive propagation failed global commutator/PSD; no covariance repair')
    return result, {'method': 'DOP853', 'rtol': rtol, 'atol': atol, 'evaluations': ode.nfev,
                    'maximum_local_commutator_relative_residual': max_local_residual}
