"""Versioned frequency and independent-input contracts for quantum research."""

from dataclasses import dataclass

import numpy as np


CONVENTION_ID = "gabes-quantum-v1"


def readonly_array(value, *, real=False):
    """Own a finite array; reject complex data at real-quadrature boundaries."""
    raw = np.asarray(value)
    if real and np.iscomplexobj(raw) and np.any(raw.imag != 0):
        raise ValueError("real coordinates required; convert the complex basis explicitly")
    result = np.array(raw.real if real else raw,
                      dtype=float if real else complex, copy=True)
    if not np.isfinite(result).all():
        raise ValueError("array must be finite")
    # An immutable byte buffer also prevents callers from re-enabling writes.
    return np.frombuffer(result.tobytes(), dtype=result.dtype).reshape(result.shape)


@dataclass(frozen=True)
class AnalysisFrequencyAxis:
    """Independent RF angular frequencies, not optical/TPD scan coordinates.

    Signed frequencies and DC are allowed; inverse Fourier convention is
    f(t) = integral f(Omega) exp(-i Omega t) dOmega/(2*pi).
    """

    omega_rad_s: np.ndarray

    def __post_init__(self):
        values = readonly_array(self.omega_rad_s, real=True)
        if values.ndim != 1 or not values.size:
            raise ValueError("analysis frequency axis must be a nonempty vector")
        if values.size > 1 and np.any(np.diff(values) <= 0):
            raise ValueError("analysis frequencies must be strictly increasing")
        object.__setattr__(self, "omega_rad_s", values)

    @classmethod
    def from_hz(cls, frequencies_hz):
        return cls(2 * np.pi * readonly_array(frequencies_hz, real=True))

    @property
    def frequency_hz(self):
        return self.omega_rad_s / (2 * np.pi)


@dataclass(frozen=True)
class OpticalDetunings:
    """One operating condition; neither value is an RF analysis frequency."""

    one_photon_rad_s: float
    two_photon_rad_s: float

    def __post_init__(self):
        for name in ("one_photon_rad_s", "two_photon_rad_s"):
            value = float(getattr(self, name))
            if not np.isfinite(value):
                raise ValueError("optical detunings must be finite")
            object.__setattr__(self, name, value)


@dataclass(frozen=True)
class GeneratorFrequencyAxis:
    """Frequencies relative to a declared atomic generator frame, in rad/s.

    These can include optical beat offsets and must not be labelled as the
    laboratory RF axis. The adapter performing that shift owns the conversion.
    """

    omega_rad_s: np.ndarray
    frame: str

    def __post_init__(self):
        values = AnalysisFrequencyAxis(self.omega_rad_s).omega_rad_s
        if not isinstance(self.frame, str) or not self.frame.strip():
            raise ValueError("generator frequency frame must be named")
        object.__setattr__(self, "omega_rad_s", values)


@dataclass(frozen=True)
class ParameterEvidence:
    """Evidence for one scalar input; independence is declared, not inferred.

    An independent beam-profile fit is allowed. A fit to the target FWM gain
    or squeezing is not. Dataset IDs retain that distinction across runs.
    """

    parameter_id: str
    value: float | None
    unit: str
    status: str  # independent | assumed | target_fitted | unknown
    source_id: str = ""
    uncertainty: float | None = None  # standard uncertainty, same units as value
    estimation_method: str = ""
    dataset_ids: tuple[str, ...] = ()
    observables_used: tuple[str, ...] = ()
    applicability: str = ""
    covariance_group: str = ""

    def __post_init__(self):
        if not self.parameter_id.strip() or not self.unit.strip():
            raise ValueError("parameter ID and unit are required (use '1' if dimensionless)")
        if self.status not in {"independent", "assumed", "target_fitted", "unknown"}:
            raise ValueError("unknown input-evidence status")
        for name in ("value", "uncertainty"):
            value = getattr(self, name)
            if value is not None:
                value = float(value)
                if not np.isfinite(value) or (name == "uncertainty" and value < 0):
                    raise ValueError(f"invalid {name}")
                object.__setattr__(self, name, value)
        for name in ("dataset_ids", "observables_used"):
            values = getattr(self, name)
            if isinstance(values, str):
                raise ValueError(f"{name} must contain nonempty strings")
            values = tuple(values)
            if any(not isinstance(v, str) or not v.strip() for v in values):
                raise ValueError(f"{name} must contain nonempty strings")
            object.__setattr__(self, name, values)


@dataclass(frozen=True)
class InputAudit:
    passed: bool
    reasons: tuple[str, ...]
    # This is a metadata gate, never an experimental-validation claim.
    scope: str = "declared input evidence only"


def audit_independent_inputs(parameters, *, required_ids, target_dataset_ids=()):
    """Check required inputs and data leakage without silently accepting omissions."""
    parameters = tuple(parameters)
    def id_set(values):
        if isinstance(values, str):
            raise ValueError("input/dataset IDs must be sequences, not strings")
        values = tuple(values)
        if any(not isinstance(v, str) or not v.strip() for v in values):
            raise ValueError("input/dataset IDs must be nonempty strings")
        return set(values)

    required = id_set(required_ids)
    target = id_set(target_dataset_ids)
    reasons = []
    if not required:
        reasons.append("required input IDs must be declared")
    ids = [p.parameter_id for p in parameters]
    if len(ids) != len(set(ids)):
        reasons.append("duplicate parameter IDs")
    for missing in sorted(required - set(ids)):
        reasons.append(f"{missing}: missing input")
    for parameter in parameters:
        prefix = parameter.parameter_id + ": "
        if parameter.status != "independent":
            reasons.append(prefix + parameter.status)
        for name in ("value", "uncertainty"):
            if getattr(parameter, name) is None:
                reasons.append(prefix + "missing " + name)
        for name in ("source_id", "estimation_method", "applicability"):
            if not getattr(parameter, name).strip():
                reasons.append(prefix + "missing " + name)
        if target.intersection(parameter.dataset_ids):
            reasons.append(prefix + "target dataset used to estimate input")
    return InputAudit(not reasons, tuple(reasons))
