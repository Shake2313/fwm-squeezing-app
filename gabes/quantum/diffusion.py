"""Single-atom stationary ordered diffusion and two-point spectra.

This is exact second-order Markov regression on a complete traceless Hermitian
operator basis, not a claim that a finite atom is a Gaussian bosonic system.
Atomic reservoirs have no field/ensemble/mode-area normalization in this module.
"""

from dataclasses import dataclass

import numpy as np

from .. import core
from .contracts import GeneratorFrequencyAxis, readonly_array
from .reservoirs import ExplicitReservoirs, audit_generator


def traceless_hermitian_basis(n_levels):
    """Hilbert-Schmidt orthonormal F_i, with Tr F_i=0, i=1..n**2-1.

    Diagonal generators precede core's symmetric/antisymmetric off-diagonals.
    The identity is excluded; means are handled through the affine drift.
    """
    if isinstance(n_levels, bool) or int(n_levels) != n_levels or n_levels < 2:
        raise ValueError("at least two integer atomic levels required")
    n = int(n_levels)
    operators = []
    for k in range(1, n):
        diagonal = np.zeros(n)
        diagonal[:k], diagonal[k] = 1, -k
        operators.append(np.diag(diagonal)/np.sqrt(k*(k+1)))
    operators.extend(core.hermitian_basis(n)[:, n:].T.reshape(-1, n, n))
    return readonly_array(operators)


def _state(value, n):
    rho = readonly_array(value)
    if rho.shape != (n, n):
        raise ValueError("atomic state dimension mismatch")
    if abs(np.trace(rho)-1) > 1e-10 or np.linalg.norm(rho-rho.conj().T) > 1e-10:
        raise ValueError("state must be trace-one and Hermitian")
    if np.linalg.eigvalsh((rho+rho.conj().T)/2)[0] < -1e-10:
        raise ValueError("state must be positive semidefinite")
    return rho


def _operators(value, n):
    operators = readonly_array(value)
    if operators.ndim != 3 or operators.shape[1:] != (n, n) or not len(operators):
        raise ValueError("operators must have shape (count,n,n)")
    for operator in operators:
        if np.linalg.norm(operator-operator.conj().T) > 1e-12*max(
                np.linalg.norm(operator), np.finfo(float).tiny):
            raise ValueError("ordered diffusion requires Hermitian operators")
    return operators


def ordered_jump_diffusion(reservoirs, rho, operators):
    """Return D[k,i,j]=Tr rho [L_k^dagger,F_i][F_j,L_k], in s**-1.

    This jump-commutator formula is evaluated without the drift or covariance.
    No ordering/symmetrization or factor of two is discarded. Each explicit
    independent reservoir has a Hermitian PSD Gram matrix for a physical rho.
    """
    if not isinstance(reservoirs, ExplicitReservoirs):
        raise TypeError("explicit reservoirs are required")
    n = reservoirs.n_levels
    rho = _state(rho, n)
    operators = _operators(operators, n)
    result = np.zeros((len(reservoirs.channels), len(operators), len(operators)), complex)
    for k, channel in enumerate(reservoirs.channels):
        commutator = operators@channel.operator-channel.operator@operators
        result[k] = np.einsum(
            "iab,jbc,ca->ij", commutator.conj().swapaxes(-1, -2), commutator, rho,
            optimize=True)
    return readonly_array(result)


@dataclass(frozen=True)
class AtomicNoiseDiagnostics:
    stationary_relative_residual: float
    stationary_minimum_eigenvalue: float
    decay_gap_s_inverse: float
    stability_resolution_s_inverse: float
    lyapunov_relative_residual: float
    commutator_relative_residual: float
    diffusion_hermiticity_relative_residual: float
    minimum_ordered_diffusion_eigenvalue_s_inverse: float
    minimum_ordered_covariance_eigenvalue: float


@dataclass(frozen=True)
class AtomicSpectrum:
    frequency_axis: GeneratorFrequencyAxis
    ordered: np.ndarray       # seconds; integral C(tau)*exp(+i*Omega*tau) d tau
    symmetrized: np.ndarray   # S_sym(Omega)=(S_ord(Omega)+S_ord(-Omega).T)/2

    def __post_init__(self):
        if not isinstance(self.frequency_axis, GeneratorFrequencyAxis):
            raise TypeError("generator-frame frequency axis required")
        ordered = readonly_array(self.ordered)
        sym = readonly_array(self.symmetrized)
        if (ordered.ndim != 3 or ordered.shape != sym.shape
                or ordered.shape[0] != len(self.frequency_axis.omega_rad_s)
                or ordered.shape[1] != ordered.shape[2]):
            raise ValueError("atomic spectra must have shape (frequency,operator,operator)")
        object.__setattr__(self, "ordered", ordered)
        object.__setattr__(self, "symmetrized", sym)

    def audit(self, *, rtol=1e-10):
        """Signed spectral PSD/Hermiticity checks, relative to atomic PSD units."""
        if not np.isfinite(rtol) or not 0 < rtol < 1:
            raise ValueError("rtol must be finite and between zero and one")
        result = {}
        for name in ("ordered", "symmetrized"):
            matrices = getattr(self, name)
            scales = np.maximum(np.linalg.norm(matrices, axis=(1, 2)), np.finfo(float).tiny)
            adjoints = matrices.conj().transpose(0, 2, 1)
            hp = np.linalg.norm(matrices-adjoints, axis=(1, 2))/scales
            minimum = np.linalg.eigvalsh((matrices+adjoints)/2)[:, 0]
            result[name] = {
                "passed": bool(np.all(hp <= rtol) and np.all(minimum >= -rtol*scales)),
                "max_hermiticity_relative_residual": float(np.max(hp)),
                "minimum_eigenvalue_seconds": float(np.min(minimum)),
                "minimum_relative_eigenvalue": float(np.min(minimum/scales)),
            }
        result["passed"] = all(result[name]["passed"] for name in ("ordered", "symmetrized"))
        return result


@dataclass(frozen=True)
class AtomicNoiseModel:
    """Build with stationary_atomic_noise; coordinates are atomic, not canonical."""

    generator: np.ndarray
    stationary_state: np.ndarray
    operators: np.ndarray
    means: np.ndarray
    drift: np.ndarray
    ordered_covariance: np.ndarray
    channel_diffusion: np.ndarray
    reservoir_names: tuple[str, ...]
    diagnostics: AtomicNoiseDiagnostics

    def __post_init__(self):
        for name in ("generator", "stationary_state", "operators", "means", "drift",
                     "ordered_covariance", "channel_diffusion"):
            object.__setattr__(self, name, readonly_array(
                getattr(self, name), real=name in {"means", "drift"}))

    @property
    def ordered_diffusion(self):
        return self.channel_diffusion.sum(axis=0)

    @property
    def symmetrized_diffusion(self):
        d = self.ordered_diffusion
        return (d+d.T)/2

    def spectrum(self, frequency_axis):
        if not isinstance(frequency_axis, GeneratorFrequencyAxis):
            raise TypeError("explicit generator-frame frequencies required")
        d = self.ordered_diffusion
        sym = (d+d.T)/2
        eye = np.eye(len(self.operators))
        ordered, symmetrized = [], []
        for omega in frequency_axis.omega_rad_s:
            response = np.linalg.solve(-self.drift-1j*omega*eye, eye)
            ordered.append(response@d@response.conj().T)
            symmetrized.append(response@sym@response.conj().T)
        result = AtomicSpectrum(frequency_axis, ordered, symmetrized)
        if not result.audit().get("passed"):
            raise ValueError("atomic spectrum failed numerical PSD/Hermiticity audit")
        return result


def stationary_atomic_noise(hamiltonian_rad_s, reservoirs, *, stationary_state=None):
    """Construct complete-basis atomic drift/noise for a mixing GKSL generator.

    Nondecaying traceless modes are rejected: their elastic delta peaks or
    initial-state dependence cannot be replaced with an ordinary continuous PSD.
    An explicit stationary state is accepted only after validation, never fitted.
    """
    if not isinstance(reservoirs, ExplicitReservoirs):
        raise TypeError("explicit reservoirs required")
    n = reservoirs.n_levels
    generator = reservoirs.generator(hamiltonian_rad_s)
    if not audit_generator(generator).passed:
        raise ValueError("generator failed numerical GKSL audit")
    operators = traceless_hermitian_basis(n)
    basis = operators.reshape(n*n-1, n*n).T
    complex_drift = basis.conj().T@generator@basis
    scale = max(float(np.linalg.norm(generator, 2)), np.finfo(float).tiny)
    if np.linalg.norm(complex_drift.imag) > 1e-12*scale:
        raise ValueError("Hermitian atomic coordinates did not give a real drift")
    drift = complex_drift.real
    gap = -float(np.max(np.linalg.eigvals(drift).real))
    resolution = 64*np.finfo(float).eps*(n*n)*scale
    if gap <= resolution:
        raise ValueError("nondecaying or unresolved atomic modes: continuous stationary spectrum unavailable")
    if stationary_state is None:
        stationary_state = core.steady_state_from_liouvillian(generator, n)
    rho = _state(stationary_state, n)
    stationary_residual = float(np.linalg.norm(generator@rho.reshape(-1))
                                / (scale*np.linalg.norm(rho)))
    if stationary_residual > 1e-10:
        raise ValueError("supplied state is not stationary for this generator")
    means = np.einsum("iab,ba->i", operators, rho).real
    centered = operators-means[:, None, None]*np.eye(n)
    covariance = np.einsum("iab,jbc,ca->ij", centered, centered, rho, optimize=True)
    by_channel = ordered_jump_diffusion(reservoirs, rho, operators)
    diffusion = by_channel.sum(axis=0)
    d_scale = max(float(np.linalg.norm(diffusion)), np.finfo(float).tiny)
    defect = drift@covariance+covariance@drift.T+diffusion
    moment_residual = float(np.linalg.norm(defect)/d_scale)
    commutator = covariance-covariance.T
    comm_defect = drift@commutator+commutator@drift.T+diffusion-diffusion.T
    # Scale by full ordered diffusion also for commuting observables at equilibrium.
    comm_residual = float(np.linalg.norm(comm_defect)/d_scale)
    hp_residual = float(np.linalg.norm(diffusion-diffusion.conj().T)/d_scale)
    d_min = float(np.linalg.eigvalsh((diffusion+diffusion.conj().T)/2)[0])
    c_min = float(np.linalg.eigvalsh((covariance+covariance.conj().T)/2)[0])
    if (moment_residual > 1e-8 or comm_residual > 1e-8 or hp_residual > 1e-10
            or d_min < -1e-10*d_scale or c_min < -1e-10):
        raise ValueError("atomic second-moment consistency failed; no covariance repair applied")
    diagnostics = AtomicNoiseDiagnostics(
        stationary_residual, float(np.linalg.eigvalsh((rho+rho.conj().T)/2)[0]),
        gap, resolution, moment_residual, comm_residual, hp_residual, d_min, c_min)
    return AtomicNoiseModel(generator, rho, operators, means, drift, covariance,
                            by_channel, tuple(c.name for c in reservoirs.channels), diagnostics)
