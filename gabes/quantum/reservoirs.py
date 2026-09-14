"""Explicit atomic reservoirs and complete-positivity checks.

Matrices use core's row-major density vectorization. A successful generator
audit proves numerical GKSL consistency at the supplied setting, not reservoir
provenance or the physical accuracy of a collision/transport model.
"""

from dataclasses import dataclass
from math import isqrt

import numpy as np

from ..core import comm_super
from .contracts import readonly_array


def _operator(value):
    result = readonly_array(value)
    if result.ndim != 2 or result.shape[0] != result.shape[1] or not result.size:
        raise ValueError("operator must be a nonempty square matrix")
    return result


@dataclass(frozen=True)
class CollapseChannel:
    """A jump operator in s**(-1/2), including sqrt(rate)."""

    name: str
    operator: np.ndarray
    kind: str
    source: str

    def __post_init__(self):
        if not all(isinstance(x, str) and x.strip()
                   for x in (self.name, self.kind, self.source)):
            raise ValueError("channel name, kind and source are required")
        object.__setattr__(self, "operator", _operator(self.operator))


@dataclass(frozen=True)
class ExplicitReservoirs:
    n_levels: int
    channels: tuple[CollapseChannel, ...]

    def __post_init__(self):
        if isinstance(self.n_levels, bool) or int(self.n_levels) != self.n_levels or self.n_levels < 1:
            raise ValueError("n_levels must be a positive integer")
        object.__setattr__(self, "n_levels", int(self.n_levels))
        channels = tuple(self.channels)
        if len({c.name for c in channels}) != len(channels):
            raise ValueError("channel names must be unique")
        if any(c.operator.shape != (self.n_levels, self.n_levels) for c in channels):
            raise ValueError("collapse operator dimension mismatch")
        object.__setattr__(self, "channels", channels)

    def dissipator(self):
        n = self.n_levels
        eye = np.eye(n)
        result = np.zeros((n*n, n*n), complex)
        for channel in self.channels:
            jump = channel.operator
            product = jump.conj().T @ jump
            result += (np.kron(jump, jump.conj())
                       - 0.5 * np.kron(product, eye)
                       - 0.5 * np.kron(eye, product.T))
        return result

    def generator(self, hamiltonian_rad_s):
        h = _operator(hamiltonian_rad_s)
        if h.shape != (self.n_levels, self.n_levels):
            raise ValueError("Hamiltonian dimension mismatch")
        scale = max(float(np.linalg.norm(h)), np.finfo(float).tiny)
        if np.linalg.norm(h-h.conj().T) > 1e-12 * scale:
            raise ValueError("Hamiltonian must be Hermitian")
        return comm_super(h) + self.dissipator()

    @classmethod
    def from_atom(cls, atom, *, source):
        """Import only explicit channels, with exact assembly parity.

        Nonzero legacy coherence dephasing and untracked superoperator edits
        are rejected even if their *combined* generator passes a CP audit.
        """
        if any(rate != 0 for _i, _j, rate in atom.dephasing):
            raise ValueError("nonzero legacy dephasing has no explicit reservoir construction")
        channels = []
        n = atom.n_levels
        for index, (excited, ground, rate) in enumerate(atom.decay):
            if not np.isfinite(rate) or rate < 0:
                raise ValueError("decay rates must be finite and non-negative")
            jump = np.zeros((n, n), complex)
            jump[ground, excited] = np.sqrt(rate)
            channels.append(CollapseChannel(
                f"decay:{index}:{excited}->{ground}", jump,
                "spontaneous_emission", source))
        for attr, kind in (("emission_ops", "spontaneous_emission"),
                           ("collapse_ops", "other")):
            for index, jump in enumerate(getattr(atom, attr)):
                channels.append(CollapseChannel(f"{attr}:{index}", jump, kind, source))
        result = cls(n, tuple(channels))
        supplied = _operator(atom.lindblad)
        assembled = result.dissipator()
        if supplied.shape != assembled.shape:
            raise ValueError("atomic generator dimension mismatch")
        scale = max(float(np.linalg.norm(supplied)),
                    float(np.linalg.norm(assembled)), np.finfo(float).tiny)
        if np.linalg.norm(supplied-assembled) > 1e-12 * scale:
            raise ValueError("atomic generator contains untracked dissipator terms")
        return result


def thermal_reset_channels(rate_s_inverse, populations, *, source, name="transit"):
    """J_ij=sqrt(rate*p_i)|i><j| gives rate*(rho_th*Tr(rho)-rho).

    Populations specify a diagonal replacement state and must already sum to
    one. This is a declared Markov replacement model, not a beam-transport law.
    """
    rate = float(rate_s_inverse)
    p = readonly_array(populations, real=True)
    if not np.isfinite(rate) or rate < 0:
        raise ValueError("reset rate must be finite and non-negative")
    if p.ndim != 1 or not p.size or np.any(p < 0):
        raise ValueError("reset populations must be a nonnegative vector")
    if abs(float(p.sum())-1.0) > 1e-12:
        raise ValueError("reset populations must sum to one")
    channels = []
    for i, population in enumerate(p):
        if population == 0 or rate == 0:
            continue
        for j in range(p.size):
            jump = np.zeros((p.size, p.size), complex)
            jump[i, j] = np.sqrt(rate*population)
            channels.append(CollapseChannel(
                f"{name}:{i}<-{j}", jump, "transit_reset", source))
    return tuple(channels)


def choi_matrix(superoperator):
    """Unnormalized J(S)=sum_ij |i><j| tensor S(|i><j|); no symmetrization."""
    matrix = _operator(superoperator)
    n = isqrt(matrix.shape[0])
    if n*n != matrix.shape[0]:
        raise ValueError("superoperator dimension must be a perfect square")
    return matrix.reshape(n, n, n, n).transpose(2, 0, 3, 1).reshape(n*n, n*n)


@dataclass(frozen=True)
class GeneratorAudit:
    passed: bool
    trace_relative_residual: float
    hermiticity_relative_residual: float
    minimum_conditional_choi_eigenvalue_s_inverse: float
    conditional_choi_tolerance_s_inverse: float
    generator_norm_s_inverse: float
    conditional_choi_norm_s_inverse: float
    scope: str = "numerical GKSL consistency of the supplied generator only"


def audit_generator(generator, *, rtol=1e-10):
    """Check TP, HP and conditional CP, retaining the signed raw eigenvalue.

    CCP tolerance uses the projected dissipative scale, plus a floating-point
    floor for projection cancellation. A large Hamiltonian is not the rtol
    scale for dissipative negativity. No eigenvalue is clipped or repaired.
    """
    if not np.isfinite(rtol) or not 0 < rtol < 1:
        raise ValueError("rtol must be finite and between zero and one")
    matrix = _operator(generator)
    choi = choi_matrix(matrix)
    n = isqrt(matrix.shape[0])
    norm = float(np.linalg.norm(matrix, 2))
    scale = max(norm, np.finfo(float).tiny)
    trace = np.eye(n).reshape(-1)
    trace_residual = float(np.linalg.norm(trace@matrix) / (np.sqrt(n)*scale))
    hp_residual = float(np.linalg.norm(choi-choi.conj().T, 2) / scale)
    omega = trace/np.sqrt(n)
    p = np.eye(n*n)-np.outer(omega, omega)
    projected = p @ ((choi+choi.conj().T)/2) @ p
    eigenvalues = np.linalg.eigvalsh((projected+projected.conj().T)/2)
    minimum = float(eigenvalues[0])
    dissipative_norm = float(np.max(np.abs(eigenvalues)))
    tolerance = rtol*dissipative_norm + 64*np.finfo(float).eps*norm
    return GeneratorAudit(
        bool(trace_residual <= rtol and hp_residual <= rtol and minimum >= -tolerance),
        trace_residual, hp_residual, minimum, float(tolerance), norm, dissipative_norm)
