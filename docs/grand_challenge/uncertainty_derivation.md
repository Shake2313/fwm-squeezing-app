# Correlated uncertainty in independently specified inputs

`gabes.quantum.uncertainty` propagates uncertainty in the scalar inputs of a
research calculation. This is separate from the quantum Langevin diffusion and
photocurrent spectra inside that calculation. It neither adds fitted field noise
nor changes the production FWM model.

The backend accepts an ordered tuple of `ParameterEvidence` records. Every
record must supply its numerical value, unit, and standard uncertainty; exact
inputs require an explicit zero uncertainty. Input order and units remain
attached to the model. There are no implicit unit conversions, default
uncertainties, or correlations inferred from `covariance_group` labels.

## Joint covariance and numerical scaling

For input means $\mu_i$, standard uncertainties $s_i$, and an explicitly supplied
correlation matrix $R$, the input covariance is

\[
C_{ij}=s_i R_{ij}s_j.
\]

`JointInputUncertainty(parameters, correlation, correlation_evidence=None)`
requires a symmetric positive-semidefinite correlation matrix with unit diagonal.
An exact input has zero off-diagonal correlations; its diagonal one is only a
latent convention and its covariance row and column are exactly zero.
`JointInputUncertainty.from_covariance(parameters, covariance, ...)` accepts
covariance elements in the products of the declared units and checks their
diagonal against the supplied marginal uncertainties. It does not replace those
uncertainties with values extracted from the matrix.

Validation and eigendecomposition use the dimensionless $R$. Applying one
absolute covariance tolerance to, for example, number density and atomic dipole
would hide errors in the smaller coordinate. Symmetry and diagonal tolerance is
$128\epsilon n$; the negative-eigenvalue tolerance additionally multiplies this
by $\max(1,\lambda_{\max})$. Symmetrization and clipping of negative eigenvalues
are limited to these floating-point tolerances. There is no diagonal jitter,
regularization, or repair of a materially nonpositive covariance.

The supplied correlation is retained verbatim. Public diagnostics expose the
signed `minimum_correlation_eigenvalue`, `correlation_psd_tolerance`,
`clipped_negative_eigenvalue_magnitude`, and the maximum elementwise difference
`factorization_correlation_max_error` between the factorized covariance and the
supplied correlation. These characterize numerical factorization, not physical
noise corrections. `minimum_correlation_eigenvalue` describes the full supplied
matrix; `clipped_negative_eigenvalue_magnitude` describes the matrix actually
factorized after the following exact consolidation.

Inputs with exactly $R_{ij}=+1$ or $-1$ share a single standardized error, with
the declared sign. Their complete rows and columns must agree exactly up to
that sign; inconsistent exact-correlation declarations are rejected. These
inputs are consolidated before factorization and reconstructed from the same
latent draw. Their signed Jacobian columns are summed before propagation. This
preserves exact shared-error cancellation: directly diagonalizing a large
singular matrix can produce tiny spurious positive eigenvalues that become
visible after a large derivative amplifies them. Nearly-perfect correlations
are not consolidated, and positive eigenvalues are never dropped by a numerical
rank threshold.

Other singular or nearly singular matrices still use ordinary floating-point
eigendecomposition. Small factorization error can be amplified by a large or
strongly cancelling Jacobian; the reported correlation-factorization defect is
not an automatic bound on output accuracy. General exact linear constraints
beyond identical or opposite input errors need a dedicated parameterization or
an independent higher-precision check when such amplification matters.

## First-order propagation

For a callback $f$ with real scalar or vector output and Jacobian
$J_{ai}=\partial f_a/\partial x_i$, the first-order result is

\[
y_0=f(\mu),\qquad U_y=JCJ^T.
\]

The implementation factors $R=BB^T$ numerically and forms
$A=(J\,\mathrm{diag}(s))B$, then $U_y=AA^T$. Scaling the Jacobian columns before
multiplication avoids constructing a poorly scaled SI covariance in this path.
The nominal output is $f(\mu)$, which need not equal the nonlinear ensemble mean.

```python
result = propagate_first_order(
    model, callback,
    output_ids=("probe_gain", "conjugate_gain"),
    output_units=("1", "1"),
    steps=[0.0001, 1e-7],  # one physical-unit step per ordered input
)
```

Supply exactly one of a real `jacobian` with shape `(outputs, inputs)` or an
explicit vector of physical-unit central-difference `steps`. Callbacks receive
an immutable one-dimensional array in `model.parameter_ids` order. Inputs with
exact zero uncertainty are not perturbed; their numerical Jacobian columns are
zero placeholders, not estimates of their physical derivatives. Callers must
check numerical step convergence and model convergence separately. A domain
failure on either side of a central difference is reported with the parameter
name and the failed input vector; no one-sided fallback is substituted.

`FirstOrderResult` retains `model`, `output_ids`, `output_units`, `nominal`,
`jacobian`, `covariance`, and the actual `steps` (or `None` for a supplied
Jacobian), and exposes `standard_uncertainties`.

## Reproducible nonlinear ensembles

`sample_inputs(model, count, seed=...)` explicitly assumes an **untruncated joint
Gaussian** input law. Means and standard uncertainties alone do not establish
that law. The sampler uses a private NumPy PCG64 generator and records its seed;
it leaves the global random state unchanged. Samples use $x=\mu+\mathrm{diag}(s)Bz$
with independent standard-normal $z$. Exact parameters are copied bit for bit.
Replay is deterministic for the same model, seed, and numerical software
environment; eigensolver and normal-generator implementations should be recorded
when comparing environments.

```python
ensemble = propagate_ensemble(
    model, callback, count=1000, seed=20260911,
    output_ids=("probe_gain", "conjugate_gain"),
    output_units=("1", "1"),
)
```

`EnsembleResult.samples` retains the joint model, full input array, seed,
distribution label, and generator label. Its `outputs`, `mean`, `covariance`,
`output_ids`, `output_units`, and `standard_uncertainties` describe the evaluated
outputs. The output covariance uses the unbiased $N-1$ denominator and requires
at least two samples. Individual input/output samples remain available; a
finite ensemble is not a Monte Carlo convergence certificate or a confidence
interval.

Callbacks must reject unphysical inputs or outputs themselves, for example a
negative sampled optical power or an unstable solver result. The generic layer
also rejects nonfinite, complex, and incorrectly shaped outputs. At the first
failure it raises `InputEvaluationError`, retaining `context`, `input_values`,
and the chained original error. No partial covariance is returned, no sample is
dropped, and no clipping, truncation, retry, or target-based resampling occurs.
If substantial probability lies beyond the physical domain, the declared joint
Gaussian model must be reconsidered using independent input information.

## Evidence and limitations

`model.audit_evidence(required_ids=..., target_dataset_ids=...)` reuses the
existing scalar `audit_independent_inputs` contract and adds separate
`CorrelationEvidence`. That record declares the status, source, estimation
method, datasets, and applicability of the whole correlation matrix, including
any assumed zero correlations. The status `independent` means independently
estimated from the target experiment, not statistically uncorrelated inputs.
Missing, assumed, unknown, target-fitted, and target-leaking evidence cannot
pass. Numerical calculations with incomplete evidence remain conditional
diagnostics and preserve the original evidence status.

By default the metadata audit checks only the model's declared parameter IDs.
Callers must supply all required consumed IDs to detect omissions and must bind
fixed callback inputs to their own consumed-input audit. A passed metadata audit
does not establish measurement quality, the Gaussian law, applicability of the
physical model, or experimental agreement. Propagation results explicitly retain
`experimental_validation=False`.

The backend does not infer uncertain reservoir models, bound omitted physics,
fit target gain or squeezing, decide whether a finite output is physical, or
automatically identify a squeezing valley. A caller can propagate edge locations
only after defining a well-posed observable with an explicit failure policy for
missing or discontinuous edges. First-order error bars become unreliable for
strong nonlinearities and discontinuities; Gaussian input draws are not a
replacement for independently characterized bounded or asymmetric input laws.

The multivariate covariance law and distribution-propagation distinction follow
the primary metrology treatments in JCGM 101 and 102, available from the
[BIPM JCGM publications](https://www.bipm.org/en/committees/jc/jcgm/publications).
This bounded implementation does not claim full GUM compliance or validation.
