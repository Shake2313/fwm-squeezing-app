"""Correlated scalar-input uncertainty, separate from quantum field noise.

All coordinates retain their declared physical units. Correlations and positive
semidefiniteness are checked in dimensionless coordinates, so a density and an
atomic dipole do not share an absolute covariance tolerance. Numerical results
are conditional on the supplied joint distribution and callback; evidence
audits are metadata checks and never experimental validation.
"""

from dataclasses import dataclass, field

import numpy as np

from .contracts import InputAudit, ParameterEvidence, audit_independent_inputs, readonly_array


def _names(values, label, *, unique=False):
    if isinstance(values, str):
        raise ValueError(f"{label} must be a sequence of nonempty strings")
    values = tuple(values)
    if any(not isinstance(v, str) or not v.strip() for v in values):
        raise ValueError(f"{label} must contain nonempty strings")
    if unique and len(set(values)) != len(values):
        raise ValueError(f"{label} must be unique")
    return values


def _exact_correlation_groups(correlation):
    """Consolidate exact common errors, without thresholding positive modes.

    In a correlation matrix, R[i,j]=+/-1 requires the complete corresponding
    rows and columns to agree up to that sign. Verify this exact identity before
    replacing them with one latent coordinate. Nearly-perfect correlations do
    not qualify, however small their remaining positive variance is.
    """
    representatives, groups, signs = [], [], []
    for index in range(len(correlation)):
        for group, representative in enumerate(representatives):
            sign = correlation[index, representative]
            if abs(sign) != 1:
                continue
            if not (np.array_equal(correlation[index], sign * correlation[representative]) and
                    np.array_equal(correlation[:, index], sign * correlation[:, representative])):
                raise ValueError("exact +/-1 correlations require identical signed rows and columns")
            groups.append(group)
            signs.append(float(sign))
            break
        else:
            groups.append(len(representatives))
            signs.append(1.)
            representatives.append(index)
    return representatives, tuple(groups), tuple(signs)


@dataclass(frozen=True)
class CorrelationEvidence:
    """Evidence for the whole supplied correlation matrix, including its zeros.

    Independence of an input's calibration from the target experiment does not
    imply that its errors are independent of the other calibrated inputs.
    """

    status: str
    source_id: str = ""
    estimation_method: str = ""
    dataset_ids: tuple[str, ...] = ()
    applicability: str = ""

    def __post_init__(self):
        if self.status not in {"independent", "assumed", "target_fitted", "unknown"}:
            raise ValueError("unknown correlation-evidence status")
        for name in ("source_id", "estimation_method", "applicability"):
            if not isinstance(getattr(self, name), str):
                raise ValueError(f"{name} must be a string")
        object.__setattr__(self, "dataset_ids", _names(self.dataset_ids, "dataset IDs"))


@dataclass(frozen=True)
class JointInputUncertainty:
    """Ordered marginal evidence and an explicit joint correlation matrix.

    Every value and standard uncertainty must be supplied, including explicit
    zero for an exact parameter. Correlation has unit diagonal; exact parameters
    have zero off-diagonal entries (their diagonal one is a latent convention).
    Missing/assumed provenance can be used for sensitivity diagnostics, but
    cannot pass :meth:`audit_evidence`. No unit conversion or marginal
    distribution shape is inferred from the evidence.
    """

    parameters: tuple[ParameterEvidence, ...]
    correlation: np.ndarray
    correlation_evidence: CorrelationEvidence | None = None
    minimum_correlation_eigenvalue: float = field(init=False)
    correlation_psd_tolerance: float = field(init=False)
    clipped_negative_eigenvalue_magnitude: float = field(init=False)
    factorization_correlation_max_error: float = field(init=False)
    _latent_factor: np.ndarray = field(init=False, repr=False)
    _group_indices: tuple[int, ...] = field(init=False, repr=False)
    _group_signs: tuple[float, ...] = field(init=False, repr=False)

    def __post_init__(self):
        parameters = tuple(self.parameters)
        if not parameters or any(not isinstance(p, ParameterEvidence) for p in parameters):
            raise ValueError("a nonempty sequence of ParameterEvidence is required")
        _names((p.parameter_id for p in parameters), "parameter IDs", unique=True)
        if any(p.value is None or p.uncertainty is None for p in parameters):
            raise ValueError("every input needs a value and explicit standard uncertainty")
        if self.correlation_evidence is not None and not isinstance(
                self.correlation_evidence, CorrelationEvidence):
            raise TypeError("correlation_evidence must be CorrelationEvidence or None")
        corr = readonly_array(self.correlation, real=True)
        n = len(parameters)
        if corr.shape != (n, n):
            raise ValueError("correlation shape must match the ordered inputs")
        tolerance = 128 * np.finfo(float).eps * n
        if not np.allclose(corr, corr.T, atol=tolerance, rtol=0):
            raise ValueError("correlation must be symmetric in dimensionless coordinates")
        if not np.allclose(corr.diagonal(), 1, atol=tolerance, rtol=0):
            raise ValueError("correlation diagonal must match unit marginal variances")
        exact = np.array([p.uncertainty == 0 for p in parameters])
        off_diagonal = np.array(corr, copy=True)
        np.fill_diagonal(off_diagonal, 0.)
        if np.any(off_diagonal[exact] != 0) or np.any(off_diagonal[:, exact] != 0):
            raise ValueError("exact inputs require zero off-diagonal correlations")
        eigenvalues, eigenvectors = np.linalg.eigh((corr + corr.T) / 2)
        psd_tolerance = tolerance * max(1., eigenvalues[-1])
        if eigenvalues[0] < -psd_tolerance:
            raise ValueError("correlation must be positive semidefinite")
        representatives, groups, signs = _exact_correlation_groups(corr)
        if len(representatives) < n:
            reduced = corr[np.ix_(representatives, representatives)]
            factor_eigenvalues, factor_eigenvectors = np.linalg.eigh((reduced + reduced.T) / 2)
        else:
            factor_eigenvalues, factor_eigenvectors = eigenvalues, eigenvectors
        # Only eigenvalues within the documented roundoff tolerance are clipped.
        latent_factor = factor_eigenvectors * np.sqrt(np.maximum(factor_eigenvalues, 0))[None, :]
        factor = latent_factor[list(groups)] * np.array(signs)[:, None]
        object.__setattr__(self, "parameters", parameters)
        object.__setattr__(self, "correlation", readonly_array(corr, real=True))
        object.__setattr__(self, "_latent_factor", readonly_array(latent_factor, real=True))
        object.__setattr__(self, "_group_indices", groups)
        object.__setattr__(self, "_group_signs", signs)
        object.__setattr__(self, "minimum_correlation_eigenvalue", float(eigenvalues[0]))
        object.__setattr__(self, "correlation_psd_tolerance", float(psd_tolerance))
        object.__setattr__(self, "clipped_negative_eigenvalue_magnitude", float(max(0., -factor_eigenvalues[0])))
        object.__setattr__(self, "factorization_correlation_max_error",
                           float(np.max(np.abs(factor @ factor.T - corr))))

    @classmethod
    def from_covariance(cls, parameters, covariance, *, correlation_evidence=None):
        """Accept C[i,j] in unit[i]*unit[j], checking each marginal's own scale.

        The supplied diagonal must agree with the evidence uncertainties; it
        never overwrites them. Exact parameters require identically zero rows
        and columns. A scalar absolute tolerance on SI covariance is not used.
        """
        parameters = tuple(parameters)
        if not parameters or any(not isinstance(p, ParameterEvidence) or
                p.value is None or p.uncertainty is None for p in parameters):
            raise ValueError("every input needs ParameterEvidence with value and uncertainty")
        covariance = readonly_array(covariance, real=True)
        n = len(parameters)
        if covariance.shape != (n, n):
            raise ValueError("covariance shape must match the ordered inputs")
        sigma = np.array([p.uncertainty for p in parameters])
        exact = sigma == 0
        if np.any(covariance[exact] != 0) or np.any(covariance[:, exact] != 0):
            raise ValueError("exact inputs require zero covariance rows and columns")
        scale = np.where(exact, 1., sigma)
        with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
            correlation = (covariance / scale[:, None]) / scale[None, :]
        correlation[exact, exact] = 1.
        return cls(parameters, correlation, correlation_evidence)

    @property
    def parameter_ids(self):
        return tuple(p.parameter_id for p in self.parameters)

    @property
    def units(self):
        return tuple(p.unit for p in self.parameters)

    @property
    def values(self):
        return readonly_array([p.value for p in self.parameters], real=True)

    @property
    def standard_uncertainties(self):
        return readonly_array([p.uncertainty for p in self.parameters], real=True)

    @property
    def covariance(self):
        sigma = self.standard_uncertainties
        with np.errstate(over="ignore", invalid="ignore"):
            covariance = self.correlation * np.outer(sigma, sigma)
        return readonly_array(covariance, real=True)

    def audit_evidence(self, *, required_ids=None, target_dataset_ids=()):
        """Audit declared marginal and correlation provenance, not the callback.

        Supply all consumed input IDs to detect omitted evidence. By default
        only this model's declared IDs are checked; fixed inputs hidden inside a
        callback are outside this scope. Correlation datasets also undergo the
        target-data leakage check.
        """
        targets = _names(target_dataset_ids, "target dataset IDs")
        required = self.parameter_ids if required_ids is None else required_ids
        base = audit_independent_inputs(self.parameters, required_ids=required,
                                        target_dataset_ids=targets)
        reasons = list(base.reasons)
        record = self.correlation_evidence
        if record is None:
            reasons.append("correlation: missing joint evidence")
        else:
            if record.status != "independent":
                reasons.append("correlation: " + record.status)
            for name in ("source_id", "estimation_method", "applicability"):
                if not getattr(record, name).strip():
                    reasons.append("correlation: missing " + name)
            if set(targets).intersection(record.dataset_ids):
                reasons.append("correlation: target dataset used to estimate input correlations")
        return InputAudit(not reasons, tuple(reasons),
                          "declared scalar and correlation evidence only; no experimental validation")


@dataclass(frozen=True)
class InputSamples:
    model: JointInputUncertainty
    values: np.ndarray
    seed: int
    distribution: str = field(default="untruncated joint Gaussian", init=False)
    generator: str = field(default="numpy.random.PCG64", init=False)

    def __post_init__(self):
        object.__setattr__(self, "values", readonly_array(self.values, real=True))


def _integer(value, label, minimum):
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, (int, np.integer)) or value < minimum:
        raise ValueError(f"{label} must be an integer >= {minimum}")
    return int(value)


def sample_inputs(model, count, *, seed):
    """Draw every requested sample once using an explicitly seeded PCG64.

    The Gaussian law is an explicit modeling choice, not a consequence of
    standard uncertainties. No bounds, truncation, rejection or target-based
    resampling are applied. Callback/domain validation belongs to propagation.
    """
    count, seed = _integer(count, "count", 1), _integer(seed, "seed", 0)
    rng = np.random.Generator(np.random.PCG64(seed))
    latent = rng.standard_normal((count, len(model._latent_factor))) @ model._latent_factor.T
    # Indexing one latent draw preserves an exact common error; a full singular
    # eigendecomposition could introduce small spurious independent errors.
    standardized = latent[:, model._group_indices] * np.array(model._group_signs)
    with np.errstate(over="ignore", invalid="ignore"):
        values = model.values + standardized * model.standard_uncertainties
    # Preserve exact coordinates, including the sign bit of an exact zero.
    exact = model.standard_uncertainties == 0
    values[:, exact] = model.values[exact]
    return InputSamples(model, values, seed)


class InputEvaluationError(ValueError):
    """A callback failed at an identified point; no samples have been discarded."""

    def __init__(self, context, input_values, cause):
        self.context = context
        self.input_values = readonly_array(input_values, real=True)
        super().__init__(f"{context}: callback failed; no samples discarded or resampled: {cause}")


def _evaluate(callback, values, output_count, context):
    values = readonly_array(values, real=True)
    try:
        output = readonly_array(callback(values), real=True)
        if output.ndim == 0:
            output = output.reshape(1)
        if output.shape != (output_count,):
            raise ValueError("callback must return one real scalar per output ID")
        return output
    except Exception as exc:
        raise InputEvaluationError(context, values, exc) from exc


def _output_schema(output_ids, output_units):
    ids = _names(output_ids, "output IDs", unique=True)
    units = _names(output_units, "output units")
    if not ids or len(units) != len(ids):
        raise ValueError("each output needs an ID and explicit unit")
    return ids, units


@dataclass(frozen=True)
class FirstOrderResult:
    model: JointInputUncertainty
    output_ids: tuple[str, ...]
    output_units: tuple[str, ...]
    nominal: np.ndarray
    jacobian: np.ndarray
    covariance: np.ndarray
    steps: np.ndarray | None
    experimental_validation: bool = field(default=False, init=False)

    def __post_init__(self):
        for name in ("nominal", "jacobian", "covariance", "steps"):
            if getattr(self, name) is not None:
                object.__setattr__(self, name, readonly_array(getattr(self, name), real=True))

    @property
    def standard_uncertainties(self):
        return readonly_array(np.sqrt(self.covariance.diagonal()), real=True)


def propagate_first_order(model, callback, *, output_ids, output_units, jacobian=None, steps=None):
    """Return f(mean) and J C J.T, using physical-unit derivatives.

    Supply a real (outputs, inputs) Jacobian, or an explicit per-input vector of
    central-difference steps in the input units. Step selection/convergence is
    the caller's responsibility. Exact inputs are not perturbed and their
    numerical Jacobian columns are zero (not claims about their derivatives).
    Callbacks receive an immutable input vector in ``model.parameter_ids`` order.
    """
    ids, units = _output_schema(output_ids, output_units)
    n = len(model.parameters)
    if (jacobian is None) == (steps is None):
        raise ValueError("supply exactly one of jacobian or explicit central-difference steps")
    nominal = _evaluate(callback, model.values, len(ids), "nominal input")
    sigma = model.standard_uncertainties
    if jacobian is not None:
        jacobian = readonly_array(jacobian, real=True)
        if jacobian.shape != (len(ids), n):
            raise ValueError("Jacobian shape must be (outputs, inputs)")
    else:
        steps = readonly_array(steps, real=True)
        if steps.shape != (n,) or np.any(steps < 0) or np.any(steps[sigma > 0] <= 0):
            raise ValueError("positive physical-unit steps required for every uncertain input")
        jacobian = np.zeros((len(ids), n))
        steps = np.where(sigma > 0, steps, 0.)
        for index in np.flatnonzero(sigma > 0):
            plus, minus = np.array(model.values), np.array(model.values)
            plus[index] += steps[index]
            minus[index] -= steps[index]
            if not np.isfinite([plus[index], minus[index]]).all() or (
                    plus[index] == model.values[index] or minus[index] == model.values[index]):
                raise ValueError(f"{model.parameter_ids[index]}: step is not finite and representable")
            name = model.parameter_ids[index]
            high = _evaluate(callback, plus, len(ids), f"central difference + {name}")
            low = _evaluate(callback, minus, len(ids), f"central difference - {name}")
            # Form the quotient in two halves so a finite +/- step spanning
            # almost the floating-point range does not overflow its denominator.
            denominator = plus[index]/2 - minus[index]/2
            jacobian[:, index] = (high/2 - low/2) / denominator
    # Scaling before multiplication avoids a badly conditioned raw SI covariance
    # and gives a positive-semidefinite output covariance even at singular rho.
    with np.errstate(over="ignore", invalid="ignore"):
        scaled_jacobian = jacobian * sigma[None, :]
        common_response = np.zeros((len(ids), len(model._latent_factor)))
        # Sum derivatives of exact common errors before multiplying their
        # correlation factor, preserving cancellation of equal opposite terms.
        for index, (group, sign) in enumerate(zip(model._group_indices, model._group_signs)):
            common_response[:, group] += sign * scaled_jacobian[:, index]
        response_factor = common_response @ model._latent_factor
        covariance = response_factor @ response_factor.T
    return FirstOrderResult(model, ids, units, nominal, jacobian, covariance, steps)


@dataclass(frozen=True)
class EnsembleResult:
    samples: InputSamples
    output_ids: tuple[str, ...]
    output_units: tuple[str, ...]
    outputs: np.ndarray
    mean: np.ndarray
    covariance: np.ndarray
    experimental_validation: bool = field(default=False, init=False)

    def __post_init__(self):
        for name in ("outputs", "mean", "covariance"):
            object.__setattr__(self, name, readonly_array(getattr(self, name), real=True))

    @property
    def standard_uncertainties(self):
        return readonly_array(np.sqrt(self.covariance.diagonal()), real=True)


def propagate_ensemble(model, callback, *, count, seed, output_ids, output_units):
    """Evaluate all Gaussian samples; fail at the first invalid callback result.

    Outputs and inputs are retained for inspection. Covariance is the unbiased
    finite-sample estimate (count-1 denominator), not a confidence interval or a
    numerical-convergence certificate. Domain violations must raise in the
    callback; finite unphysical numbers cannot be inferred by this generic layer.
    """
    count = _integer(count, "count", 2)
    ids, units = _output_schema(output_ids, output_units)
    samples = sample_inputs(model, count, seed=seed)
    outputs = np.array([_evaluate(callback, point, len(ids), f"sample {index}, seed {samples.seed}")
                        for index, point in enumerate(samples.values)])
    # Referencing a sample prevents avoidable overflow from summing a large
    # common offset; constant outputs have exactly zero covariance.
    offsets = outputs - outputs[0]
    centered_mean = offsets.mean(axis=0)
    mean = outputs[0] + centered_mean
    deviations = offsets - centered_mean
    covariance = deviations.T @ deviations / (count - 1)
    return EnsembleResult(samples, ids, units, outputs, mean, covariance)
