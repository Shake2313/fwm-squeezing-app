"""Finite temporal-mode Gaussian channels in interleaved (x,p) coordinates.

These are algebraic/device channels. In particular a vacuum-loss fixture does
not determine the microscopic diffusion of a pumped atomic vapor.
"""

from dataclasses import dataclass

import numpy as np

from .contracts import readonly_array


def canonical_commutator(n_modes):
    if isinstance(n_modes, bool) or int(n_modes) != n_modes or n_modes < 1:
        raise ValueError("n_modes must be a positive integer")
    return np.kron(np.eye(int(n_modes)), [[0., 1.], [-1., 0.]])


def _modes(names):
    if isinstance(names, str):
        raise ValueError("mode names must be a sequence")
    names = tuple(names)
    if not names or any(not isinstance(x, str) or not x.strip() for x in names):
        raise ValueError("nonempty mode names required")
    if len(set(names)) != len(names):
        raise ValueError("mode names must be unique")
    return names


def _covariance(value, dimension):
    matrix = readonly_array(value, real=True)
    if matrix.shape != (dimension, dimension):
        raise ValueError("covariance dimension mismatch")
    if np.linalg.norm(matrix-matrix.T, 2) > 1e-12*max(1., np.linalg.norm(matrix, 2)):
        raise ValueError("covariance must be symmetric")
    return matrix


def covariance_uncertainty_minimum(covariance):
    matrix = readonly_array(covariance, real=True)
    if matrix.ndim != 2 or not matrix.shape[0] or matrix.shape[0] % 2:
        raise ValueError("covariance must have positive even dimension")
    matrix = _covariance(matrix, matrix.shape[0])
    j = canonical_commutator(matrix.shape[0]//2)
    return float(np.linalg.eigvalsh((matrix+matrix.T)/2 + 0.5j*j)[0])


@dataclass(frozen=True)
class ChannelAudit:
    passed: bool
    minimum_noise_eigenvalue: float
    minimum_cp_eigenvalue: float
    tolerance: float
    scope: str = "finite temporal-mode channel consistency only"


@dataclass(frozen=True)
class GaussianChannel:
    """V_out = X V_in X.T + Y, vacuum I/2; mean fields are separate.

    Input/output mode labels are ordered and checked during composition.
    Rectangular channels (e.g. discarding uncollected modes) are supported.
    """

    transfer: np.ndarray
    added_covariance: np.ndarray
    input_modes: tuple[str, ...]
    output_modes: tuple[str, ...]
    source: str

    def __post_init__(self):
        inputs, outputs = _modes(self.input_modes), _modes(self.output_modes)
        x = readonly_array(self.transfer, real=True)
        if x.shape != (2*len(outputs), 2*len(inputs)):
            raise ValueError("transfer dimension does not match named modes")
        y = _covariance(self.added_covariance, 2*len(outputs))
        if not isinstance(self.source, str) or not self.source.strip():
            raise ValueError("channel source is required")
        object.__setattr__(self, "input_modes", inputs)
        object.__setattr__(self, "output_modes", outputs)
        object.__setattr__(self, "transfer", x)
        object.__setattr__(self, "added_covariance", y)

    def audit(self, *, rtol=1e-10):
        if not np.isfinite(rtol) or not 0 < rtol < 1:
            raise ValueError("rtol must be finite and between zero and one")
        x, y = self.transfer, self.added_covariance
        jin = canonical_commutator(len(self.input_modes))
        jout = canonical_commutator(len(self.output_modes))
        defect = jout - x@jin@x.T
        cp = (y+y.T)/2 + 0.5j*defect
        minimum_y = float(np.linalg.eigvalsh((y+y.T)/2)[0])
        minimum_cp = float(np.linalg.eigvalsh((cp+cp.conj().T)/2)[0])
        tolerance = rtol*max(1., float(np.linalg.norm(y, 2)),
                             float(np.linalg.norm(defect, 2)))
        return ChannelAudit(minimum_y >= -tolerance and minimum_cp >= -tolerance,
                            minimum_y, minimum_cp, float(tolerance))

    def apply_covariance(self, covariance):
        if not self.audit().passed:
            raise ValueError("channel violates complete positivity")
        v = _covariance(covariance, 2*len(self.input_modes))
        tolerance = 1e-10*max(1., float(np.linalg.norm(v, 2)))
        if covariance_uncertainty_minimum(v) < -tolerance:
            raise ValueError("input covariance violates quantum uncertainty")
        result = self.transfer@v@self.transfer.T + self.added_covariance
        if not np.isfinite(result).all():
            raise ValueError("output covariance overflowed")
        output_tolerance = 1e-10*max(1., float(np.linalg.norm(result, 2)))
        if covariance_uncertainty_minimum(result) < -output_tolerance:
            raise ValueError("output covariance violates quantum uncertainty")
        return result

    @classmethod
    def identity(cls, modes):
        modes = _modes(modes)
        return cls(np.eye(2*len(modes)), np.zeros((2*len(modes), 2*len(modes))),
                   modes, modes, "analytic identity channel")

    @classmethod
    def vacuum_attenuator(cls, modes, transmission, *, source):
        modes = _modes(modes)
        eta = readonly_array(transmission, real=True)
        if eta.shape != (len(modes),) or np.any((eta < 0) | (eta > 1)):
            raise ValueError("one intensity transmission in [0,1] per mode required")
        eta = np.repeat(eta, 2)
        return cls(np.diag(np.sqrt(eta)), np.diag((1-eta)/2), modes, modes, source)


def compose_channels(first, second):
    """Propagate through first, then second; do not reorder named quadratures."""
    if first.output_modes != second.input_modes:
        raise ValueError("intermediate mode order mismatch")
    if not first.audit().passed or not second.audit().passed:
        raise ValueError("cannot compose a channel violating complete positivity")
    x = second.transfer@first.transfer
    y = (second.transfer@first.added_covariance@second.transfer.T
         + second.added_covariance)
    result = GaussianChannel(x, y, first.input_modes, second.output_modes,
                             first.source + " -> " + second.source)
    if not result.audit().passed:
        raise ValueError("composed channel failed complete positivity")
    return result
