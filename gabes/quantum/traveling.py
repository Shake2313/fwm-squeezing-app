"""Local photon-flux Nambu generators from independent atomic reservoirs.

Assumes a declared, phase-selected traveling-mode sector, fixed classical pump,
independent atoms in a common uniform transverse area, and weak quantum fields.
No real-quadrature or experimental squeezing interpretation is implicit.
"""

from dataclasses import dataclass

import numpy as np
from scipy.linalg import expm

from .contracts import GeneratorFrequencyAxis, readonly_array
from .diffusion import AtomicNoiseModel


def _signs(value):
    signs = readonly_array(value, real=True)
    if signs.ndim != 1 or not signs.size or not np.isin(signs, [-1, 1]).all():
        raise ValueError("nonempty Nambu signs (+1 annihilator, -1 creator) required")
    return signs


def _matrix_audit(matrices, rtol, ordering_scale):
    own_scale = np.linalg.norm(matrices, axis=(-2, -1))
    scales = np.maximum(np.maximum(own_scale, ordering_scale), np.finfo(float).tiny)
    adjoint = matrices.conj().swapaxes(-1, -2)
    hp_absolute = np.linalg.norm(matrices-adjoint, axis=(-2, -1))
    eig = np.linalg.eigvalsh((matrices+adjoint)/2)[..., 0]
    # Vacuum can annihilate one ordering exactly. Its tiny projected roundoff
    # must be compared to the paired reservoir scale, not divided by itself.
    tolerance = rtol*own_scale+64*np.finfo(float).eps*ordering_scale
    return {
        "passed": bool(np.all(hp_absolute <= tolerance) and np.all(eig >= -tolerance)),
        "maximum_hermiticity_relative_residual": float(np.max(hp_absolute/scales)),
        "minimum_eigenvalue": float(np.min(eig)),
        "minimum_relative_eigenvalue": float(np.min(eig/scales)),
    }


def _arrays(instance, names, *, per_channel=False):
    if not isinstance(instance.frequency_axis, GeneratorFrequencyAxis):
        raise TypeError("explicit generator-frequency axis required")
    signs = _signs(instance.signs)
    object.__setattr__(instance, "signs", signs)
    labels = tuple(instance.mode_labels)
    if (len(labels) != len(signs) or len(set(labels)) != len(labels)
            or any(not isinstance(n, str) or not n.strip() for n in labels)):
        raise ValueError("unique nonempty physical mode labels required")
    object.__setattr__(instance, "mode_labels", labels)
    shape = (len(instance.frequency_axis.omega_rad_s), len(signs), len(signs))
    for name in names:
        value = readonly_array(getattr(instance, name))
        expected = ((len(instance.reservoir_names),)+shape
                    if per_channel and name != "drift" else shape)
        if value.shape != expected:
            raise ValueError(f"{name} has incompatible frequency/mode/channel dimensions")
        object.__setattr__(instance, name, value)


@dataclass(frozen=True)
class LocalNambuGenerator:
    frequency_axis: GeneratorFrequencyAxis
    signs: np.ndarray
    mode_labels: tuple[str, ...]
    drift: np.ndarray                     # m^-1
    noise_greater_by_reservoir: np.ndarray # m^-1; <f_i f_j^dagger>
    noise_lesser_by_reservoir: np.ndarray  # m^-1; <f_j^dagger f_i>
    reservoir_names: tuple[str, ...]

    def __post_init__(self):
        names = tuple(self.reservoir_names)
        if not names or len(set(names)) != len(names) or any(not n for n in names):
            raise ValueError("unique nonempty reservoir labels required")
        object.__setattr__(self, "reservoir_names", names)
        _arrays(self, ("drift", "noise_greater_by_reservoir", "noise_lesser_by_reservoir"),
                per_channel=True)

    @property
    def noise_greater(self):
        return self.noise_greater_by_reservoir.sum(axis=0)

    @property
    def noise_lesser(self):
        return self.noise_lesser_by_reservoir.sum(axis=0)

    def audit(self, *, rtol=1e-8):
        if not np.isfinite(rtol) or not 0 < rtol < 1:
            raise ValueError("finite relative tolerance between zero and one required")
        j = np.diag(self.signs)
        terms = (self.drift@j, j@self.drift.conj().swapaxes(-1, -2),
                 self.noise_greater, -self.noise_lesser)
        scale = np.maximum(sum(np.linalg.norm(t, axis=(-2, -1)) for t in terms),
                           np.finfo(float).tiny)
        residual = np.linalg.norm(sum(terms), axis=(-2, -1))/scale
        result = {"maximum_commutator_relative_residual": float(np.max(residual))}
        ordering_scale = (np.linalg.norm(self.noise_greater_by_reservoir, axis=(-2, -1))
                          +np.linalg.norm(self.noise_lesser_by_reservoir, axis=(-2, -1)))
        for name in ("noise_greater_by_reservoir", "noise_lesser_by_reservoir"):
            result[name] = _matrix_audit(getattr(self, name), rtol, ordering_scale)
        result["passed"] = bool(np.all(residual <= rtol) and all(
            result[k]["passed"] for k in ("noise_greater_by_reservoir", "noise_lesser_by_reservoir")))
        return result


def eliminate_atomic_noise(atom, readout_operators, coupling_s_inverse_sqrt_flux,
                           signs, frequency_axis, *, linear_density_m_inverse, mode_labels):
    """Eliminate a complete atomic basis, deriving M and both noise orderings.

    H_int/hbar = sum_j g_j b_j O_j^dagger + h.c. in the selected sector.
    For number density n and shared uniform area A, lambda=n*A. The slice-average
    atomic diffusion is D_atom/(lambda*dz), hence D_field=lambda*C0 R D R† C0†.
    For velocity classes pass lambda_v=n*A*w_v separately and sum the results.
    g=d*Q/(2*hbar) has units s^-1/2; b has photon-flux units s^-1/2.
    """
    if not isinstance(atom, AtomicNoiseModel):
        raise TypeError("validated AtomicNoiseModel required")
    if not isinstance(frequency_axis, GeneratorFrequencyAxis):
        raise TypeError("explicit generator-frequency axis required")
    signs = _signs(signs)
    ops = readonly_array(readout_operators)
    n = atom.stationary_state.shape[0]
    if ops.shape != (len(signs), n, n):
        raise ValueError("one atomic readout operator per Nambu coordinate required")
    g = readonly_array(coupling_s_inverse_sqrt_flux, real=True)
    if g.shape != signs.shape or np.any(g <= 0):
        raise ValueError("one positive coupling per Nambu coordinate required")
    density = float(linear_density_m_inverse)
    if not np.isfinite(density) or density < 0:
        raise ValueError("finite nonnegative linear density required")
    rho, f = atom.stationary_state, atom.operators
    w = np.einsum("jab,iba->ji", ops, f)
    # Density-operator driving commutators; neither noise nor atomic K is used.
    sources = np.column_stack([
        np.einsum("iab,ba->i", f, -1j*g[k]*(o.conj().T@rho-rho@o.conj().T))
        for k, o in enumerate(ops)])
    c0 = -1j*np.diag(signs*g)@w
    m = len(signs)
    shape = (len(frequency_axis.omega_rad_s), m, m)
    drift = np.zeros(shape, complex)
    greater = np.zeros((len(atom.reservoir_names),)+shape, complex)
    lesser = np.zeros_like(greater)
    eye = np.eye(len(f))
    for idx, omega in enumerate(frequency_axis.omega_rad_s):
        cr = c0@np.linalg.solve(-atom.drift-1j*omega*eye, eye)
        drift[idx] = density*cr@sources
        greater[:, idx] = density*(cr@atom.channel_diffusion@cr.conj().T)
        lesser[:, idx] = density*(cr@atom.channel_diffusion.swapaxes(-1, -2)@cr.conj().T)
    result = LocalNambuGenerator(frequency_axis, signs, mode_labels, drift, greater, lesser, atom.reservoir_names)
    if not result.audit()["passed"]:
        raise ValueError("microscopic field generator failed commutator/PSD audit")
    return result


def _same_coordinates(first, second):
    if (first.frequency_axis.frame != second.frequency_axis.frame
            or not np.array_equal(first.frequency_axis.omega_rad_s, second.frequency_axis.omega_rad_s)
            or first.mode_labels != second.mode_labels
            or not np.array_equal(first.signs, second.signs)):
        raise ValueError("identical generator frequency frames, axes and Nambu signs required")


def sum_independent_classes(classes):
    """Each class already contains its number density: sum covariances, not amplitudes."""
    classes = tuple(classes)
    if not classes or not all(isinstance(c, LocalNambuGenerator) for c in classes):
        raise TypeError("one or more local atomic class generators required")
    first = classes[0]
    for item in classes:
        _same_coordinates(first, item)
        if not item.audit()["passed"]:
            raise ValueError("invalid atomic class")
    return LocalNambuGenerator(
        first.frequency_axis, first.signs, first.mode_labels, sum(c.drift for c in classes),
        np.concatenate([c.noise_greater_by_reservoir for c in classes]),
        np.concatenate([c.noise_lesser_by_reservoir for c in classes]),
        tuple(f"class{idx}:{name}" for idx, c in enumerate(classes) for name in c.reservoir_names))


@dataclass(frozen=True)
class NambuTransfer:
    frequency_axis: GeneratorFrequencyAxis
    signs: np.ndarray
    mode_labels: tuple[str, ...]
    transfer: np.ndarray
    noise_greater: np.ndarray
    noise_lesser: np.ndarray

    def __post_init__(self):
        _arrays(self, ("transfer", "noise_greater", "noise_lesser"))

    def audit(self, *, rtol=1e-8):
        if not np.isfinite(rtol) or not 0 < rtol < 1:
            raise ValueError("finite relative tolerance between zero and one required")
        j = np.diag(self.signs)
        transported = self.transfer@j@self.transfer.conj().swapaxes(-1, -2)
        terms = (transported, self.noise_greater, -self.noise_lesser)
        scale = np.maximum(sum(np.linalg.norm(t, axis=(-2, -1)) for t in terms)
                           +np.linalg.norm(j), np.finfo(float).tiny)
        residual = np.linalg.norm(sum(terms)-j, axis=(-2, -1))/scale
        result = {"maximum_commutator_relative_residual": float(np.max(residual))}
        ordering_scale = (np.linalg.norm(self.noise_greater, axis=(-2, -1))
                          +np.linalg.norm(self.noise_lesser, axis=(-2, -1)))
        for name in ("noise_greater", "noise_lesser"):
            result[name] = _matrix_audit(getattr(self, name), rtol, ordering_scale)
        result["passed"] = bool(np.all(residual <= rtol) and all(
            result[k]["passed"] for k in ("noise_greater", "noise_lesser")))
        return result

    def vacuum_output(self):
        """Ordered spectral covariances for vacuum inputs; not a detector PSD."""
        if not self.audit()["passed"]:
            raise ValueError("invalid field transfer")
        t, td = self.transfer, self.transfer.conj().swapaxes(-1, -2)
        return (t@np.diag((1+self.signs)/2)@td+self.noise_greater,
                t@np.diag((1-self.signs)/2)@td+self.noise_lesser)


def constant_segment(local, length_m):
    """Exact constant-M covariance integral using an augmented Kronecker exponential.

    Pump/state/density are prescribed within this segment. No depletion is solved.
    Retarded time removes vacuum propagation; declared phase mismatch stays in M.
    """
    if not isinstance(local, LocalNambuGenerator) or not local.audit()["passed"]:
        raise ValueError("valid microscopic local generator required")
    length = float(length_m)
    if not np.isfinite(length) or length < 0:
        raise ValueError("finite nonnegative segment length required")
    m = len(local.signs)
    identity = np.eye(m)
    transfers, greater, lesser = [], [], []
    for drift, dp, dm in zip(local.drift, local.noise_greater, local.noise_lesser):
        augmented = np.zeros((m*m+2, m*m+2), complex)
        # Row-major vec(N), dN/dz=M N+N M†+D.
        augmented[:m*m, :m*m] = np.kron(drift, identity)+np.kron(identity, drift.conj())
        augmented[:m*m, -2] = dp.reshape(-1)
        augmented[:m*m, -1] = dm.reshape(-1)
        solution = expm(augmented*length)
        transfers.append(expm(drift*length))
        greater.append(solution[:m*m, -2].reshape(m, m))
        lesser.append(solution[:m*m, -1].reshape(m, m))
    result = NambuTransfer(local.frequency_axis, local.signs, local.mode_labels, transfers, greater, lesser)
    if not result.audit()["passed"]:
        raise ValueError("propagation failed commutator/PSD audit; no covariance repair applied")
    return result


def compose_segments(first, second):
    """Spatial order first then second, transporting first-segment noise through second."""
    _same_coordinates(first, second)
    if not first.audit()["passed"] or not second.audit()["passed"]:
        raise ValueError("valid segment transfers required")
    t, td = second.transfer, second.transfer.conj().swapaxes(-1, -2)
    result = NambuTransfer(first.frequency_axis, first.signs, first.mode_labels, t@first.transfer,
                           t@first.noise_greater@td+second.noise_greater,
                           t@first.noise_lesser@td+second.noise_lesser)
    if not result.audit()["passed"]:
        raise ValueError("composed propagation failed commutator/PSD audit")
    return result
