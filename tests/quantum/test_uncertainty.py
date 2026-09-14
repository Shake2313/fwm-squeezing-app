"""Independent analytic and domain-failure checks for input uncertainty."""

from dataclasses import FrozenInstanceError, replace

import numpy as np
import pytest

from gabes.quantum.contracts import ParameterEvidence
from gabes.quantum.uncertainty import (
    CorrelationEvidence, InputEvaluationError, JointInputUncertainty,
    propagate_ensemble, propagate_first_order, sample_inputs,
)


def parameter(name, value, uncertainty, unit="1", **kwargs):
    defaults = dict(status="independent", source_id="separate calibration",
                    estimation_method="independent calibration fit",
                    dataset_ids=("calibration",), applicability="fixture operating point")
    defaults.update(kwargs)
    return ParameterEvidence(name, value, unit, uncertainty=uncertainty, **defaults)


def correlation_evidence(**kwargs):
    defaults = dict(status="independent", source_id="joint calibration record",
                    estimation_method="covariance of paired calibration observations",
                    dataset_ids=("paired-calibration",), applicability="fixture operating point")
    defaults.update(kwargs)
    return CorrelationEvidence(**defaults)


def model():
    return JointInputUncertainty((parameter("a", 10., 2.), parameter("b", 20., 3.)),
                                [[1., .5], [.5, 1.]], correlation_evidence())


def test_affine_correlated_output_covariance_matches_independent_hand_calculation():
    joint = model()
    # Var(a+2b)=4+4*9+4*3=52; Var(3a-b)=9*4+9-6*3=27.
    # Cov(a+2b,3a-b)=3*4-2*9+(6-1)*3=9.
    jacobian = np.array([[1., 2.], [3., -1.]])
    callback = lambda x: jacobian @ x + [7., -4.]
    result = propagate_first_order(joint, callback, output_ids=("sum", "difference"),
                                   output_units=("1", "1"), jacobian=jacobian)
    np.testing.assert_allclose(result.nominal, [57., 6.], atol=1e-14)
    np.testing.assert_allclose(result.covariance, [[52., 9.], [9., 27.]], atol=5e-14)
    np.testing.assert_allclose(result.standard_uncertainties, np.sqrt([52., 27.]))
    assert result.steps is None and not result.experimental_validation
    numerical = propagate_first_order(joint, callback, output_ids=result.output_ids,
                                       output_units=result.output_units, steps=[.01, .02])
    np.testing.assert_allclose(numerical.covariance, result.covariance, rtol=1e-12)


def test_singular_anticorrelation_gives_exact_cancellation_without_diagonal_jitter():
    joint = JointInputUncertainty((parameter("a", 0., 2.), parameter("b", 0., 2.)),
                                  [[1., -1.], [-1., 1.]])
    result = propagate_first_order(joint, lambda x: [x.sum(), x[0]-x[1]],
        output_ids=("sum", "difference"), output_units=("1", "1"),
        jacobian=[[1., 1.], [1., -1.]])
    np.testing.assert_allclose(result.covariance, [[0., 0.], [0., 16.]], atol=1e-14)
    values = sample_inputs(joint, 100, seed=9).values
    np.testing.assert_array_equal(values[:, 0], -values[:, 1])


@pytest.mark.parametrize("signs", [[1., 1., 1., 1., 1.], [1., -1., 1., -1., -1.]])
def test_five_exact_common_errors_preserve_draw_identity_and_amplified_cancellation(signs):
    signs = np.array(signs)
    parameters = tuple(parameter(f"x{i}", 0., 1.) for i in range(5))
    joint = JointInputUncertainty(parameters, np.outer(signs, signs))
    values = sample_inputs(joint, 100, seed=0).values
    for index in range(5):
        np.testing.assert_array_equal(values[:, index], signs[index] * values[:, 0])
    jacobian = np.array([[1e9, -signs[1]*1e9, 0., 0., 0.]])
    result = propagate_first_order(joint, lambda x: jacobian @ x,
        output_ids=("amplified_difference",), output_units=("1",), jacobian=jacobian)
    np.testing.assert_array_equal(result.covariance, [[0.]])
    assert joint.factorization_correlation_max_error == 0.


def test_exact_common_errors_remain_exact_when_correlated_to_another_latent_error():
    parameters = tuple(parameter(f"x{i}", 0., 1.) for i in range(4))
    joint = JointInputUncertainty(parameters,
        [[1., -1., .3, -.3], [-1., 1., -.3, .3], [.3, -.3, 1., -1.], [-.3, .3, -1., 1.]])
    values = sample_inputs(joint, 100, seed=0).values
    np.testing.assert_array_equal(values[:, 0], -values[:, 1])
    np.testing.assert_array_equal(values[:, 2], -values[:, 3])
    result = propagate_first_order(joint, lambda x: x[:2].sum()*1e9,
        output_ids=("amplified_sum",), output_units=("1",), jacobian=[[1e9, 1e9, 0., 0.]])
    np.testing.assert_array_equal(result.covariance, [[0.]])


@pytest.mark.parametrize("rho", [1.-1e-12, 1.-8*np.finfo(float).eps])
def test_nearly_perfect_correlation_retains_small_positive_difference_variance(rho):
    # This positive mode is below the negative-roundoff tolerance; it must not
    # be dropped by a blanket rank cutoff or treated as an exact common error.
    parameters = (parameter("a", 0., 1.), parameter("b", 0., 1.))
    joint = JointInputUncertainty(parameters, [[1., rho], [rho, 1.]])
    jacobian = np.array([[1e9, -1e9]])
    result = propagate_first_order(joint, lambda x: jacobian @ x,
        output_ids=("amplified_difference",), output_units=("1",), jacobian=jacobian)
    expected = 2*(1-rho)*1e18
    assert result.covariance[0, 0] == pytest.approx(expected, rel=1e-12)
    assert result.covariance[0, 0] > 0
    samples = sample_inputs(joint, 100, seed=0).values
    assert np.any(samples[:, 0] != samples[:, 1])


def test_exact_pair_with_inconsistent_signed_rows_is_rejected_before_consolidation():
    parameters = tuple(parameter(f"x{i}", 0., 1.) for i in range(3))
    # An exact pair implies identical full rows, even if its violation is small
    # enough to hide inside the ordinary eigenvalue roundoff tolerance.
    with pytest.raises(ValueError, match="identical signed rows"):
        JointInputUncertainty(parameters, [[1., 1., .2], [1., 1., .2+1e-14], [.2, .2+1e-14, 1.]])


def test_disparate_si_scales_have_no_shared_absolute_covariance_tolerance():
    parameters = (parameter("density", 1e20, 2e18, "m^-3"),
                  parameter("dipole", 1e-29, 3e-31, "C m"))
    covariance = np.array([[4e36, 3e-13], [3e-13, 9e-62]])
    joint = JointInputUncertainty.from_covariance(parameters, covariance)
    np.testing.assert_allclose(joint.correlation, [[1., .5], [.5, 1.]], rtol=1e-14, atol=0)
    np.testing.assert_allclose(joint.covariance, covariance, rtol=1e-14, atol=0)
    callback = lambda x: x[0]/2e18 + x[1]/3e-31
    result = propagate_first_order(joint, callback, output_ids=("scaled_sum",),
                                   output_units=("1",), steps=[2e16, 3e-33])
    np.testing.assert_allclose(result.covariance, [[3.]], rtol=1e-10)
    wrong = covariance.copy()
    wrong[1, 1] *= 2
    with pytest.raises(ValueError, match="diagonal"):
        JointInputUncertainty.from_covariance(parameters, wrong)
    wrong = covariance.copy()
    wrong[0, 1] += 1e-14
    with pytest.raises(ValueError, match="symmetric"):
        JointInputUncertainty.from_covariance(parameters, wrong)


def test_extreme_but_representable_si_variances_are_validated_dimensionlessly():
    parameters = (parameter("large", 0., 1e150), parameter("small", 0., 1e-150))
    joint = JointInputUncertainty.from_covariance(parameters, [[1e300, -.4], [-.4, 1e-300]])
    np.testing.assert_allclose(joint.correlation, [[1., -.4], [-.4, 1.]], rtol=1e-14)
    result = propagate_first_order(joint, lambda x: x[0]*1e-150+x[1]*1e150,
        output_ids=("normalized",), output_units=("1",), jacobian=[[1e-150, 1e150]])
    np.testing.assert_allclose(result.covariance, [[1.2]], atol=1e-14)


@pytest.mark.parametrize("bad,match", [
    ([[1., .5], [.4, 1.]], "symmetric"),
    ([[1., 2.], [2., 1.]], "positive semidefinite"),
    ([[2., 0.], [0., 1.]], "diagonal"),
    ([[1.]], "shape"),
    ([[1., np.nan], [np.nan, 1.]], "finite"),
    ([[1., 1j], [-1j, 1.]], "real coordinates"),
])
def test_invalid_correlations_are_rejected(bad, match):
    with pytest.raises(ValueError, match=match):
        JointInputUncertainty(model().parameters, bad)


def test_covariance_requires_full_psd_not_only_valid_pairwise_correlations():
    parameters = tuple(parameter(name, 1., 1.) for name in ("a", "b", "c"))
    with pytest.raises(ValueError, match="positive semidefinite"):
        JointInputUncertainty(parameters, [[1., -.9, -.9], [-.9, 1., -.9], [-.9, -.9, 1.]])


def test_exact_parameters_are_bitwise_fixed_and_not_finite_differenced():
    parameters = (parameter("varying", 1., .2), parameter("exact", 2., 0.),
                  parameter("signed_zero", -0., 0.))
    joint = JointInputUncertainty.from_covariance(parameters, np.diag([.04, 0., 0.]))
    points = sample_inputs(joint, 200, seed=85)
    assert np.all(points.values[:, 1] == 2.)
    assert np.all(np.signbit(points.values[:, 2]))
    calls = []
    def callback(x):
        assert x[1] == 2. and np.signbit(x[2])
        calls.append(x)
        return x[0] * x[1]
    result = propagate_first_order(joint, callback, output_ids=("product",),
                                   output_units=("1",), steps=[.01, 0., 1.])
    assert len(calls) == 3
    np.testing.assert_array_equal(result.steps, [.01, 0., 0.])
    np.testing.assert_allclose(result.jacobian, [[2., 0., 0.]])
    np.testing.assert_allclose(result.covariance, [[.16]])
    for covariance in ([[.04, 1e-40, 0.], [1e-40, 0., 0.], [0., 0., 0.]],
                       [[.04, 0., 0.], [0., 1e-40, 0.], [0., 0., 0.]]):
        with pytest.raises(ValueError, match="exact inputs"):
            JointInputUncertainty.from_covariance(parameters, covariance)
    with pytest.raises(ValueError, match="exact inputs"):
        JointInputUncertainty(parameters, [[1., .1, 0.], [.1, 1., 0.], [0., 0., 1.]])


def test_sampling_replays_seed_and_matches_analytic_affine_statistics():
    joint = model()
    before = np.random.get_state()
    points = sample_inputs(joint, 100000, seed=3035)
    after = np.random.get_state()
    np.testing.assert_array_equal(before[1], after[1])
    assert before[2:] == after[2:]
    np.testing.assert_array_equal(points.values, sample_inputs(joint, 100000, seed=3035).values)
    assert not np.array_equal(points.values[:10], sample_inputs(joint, 10, seed=3036).values)
    np.testing.assert_allclose(np.mean(points.values, axis=0), [10., 20.], atol=.025, rtol=0)
    expected_covariance = np.array([[4., 3.], [3., 9.]])
    sampling_error = np.sqrt((np.outer([4., 9.], [4., 9.]) + expected_covariance**2) / 99999)
    assert np.all(np.abs(np.cov(points.values.T)-expected_covariance) < 6*sampling_error)
    assert points.seed == 3035 and points.generator == "numpy.random.PCG64"
    assert points.distribution == "untruncated joint Gaussian"
    result = propagate_ensemble(joint, lambda x: [x[0]+2*x[1], 3*x[0]-x[1]],
        count=10000, seed=20, output_ids=("sum", "difference"), output_units=("1", "1"))
    np.testing.assert_allclose(result.mean, [50., 10.], atol=.25, rtol=0)
    np.testing.assert_allclose(result.covariance, [[52., 9.], [9., 27.]], atol=.8, rtol=0)
    np.testing.assert_allclose(result.covariance, np.cov(result.outputs.T), atol=1e-13)
    assert not result.experimental_validation


def test_ensemble_retains_nonlinearity_that_first_order_omits():
    joint = JointInputUncertainty((parameter("x", 0., 2.),), [[1.]])
    linear = propagate_first_order(joint, lambda x: x[0]**2, output_ids=("square",),
                                  output_units=("1",), jacobian=[[0.]])
    ensemble = propagate_ensemble(joint, lambda x: x[0]**2, count=12000, seed=912,
                                  output_ids=("square",), output_units=("1",))
    assert linear.nominal[0] == 0 and linear.covariance[0, 0] == 0
    assert ensemble.mean[0] == pytest.approx(4., abs=.12)
    assert ensemble.covariance[0, 0] == pytest.approx(32., abs=2.)


def test_constant_exact_ensemble_has_zero_covariance_without_roundoff_variance():
    joint = JointInputUncertainty((parameter("x", 1e20, 0., "m^-3"),), [[1.]])
    result = propagate_ensemble(joint, lambda x: [x[0], 1e-29], count=31, seed=0,
                                output_ids=("density", "dipole"), output_units=("m^-3", "C m"))
    np.testing.assert_array_equal(result.mean, [1e20, 1e-29])
    np.testing.assert_array_equal(result.covariance, np.zeros((2, 2)))


def test_input_and_correlation_evidence_are_separate_and_target_leakage_cannot_pass():
    joint = model()
    assert joint.audit_evidence(target_dataset_ids=("FWM-target",)).passed
    assert not joint.audit_evidence(required_ids=("a", "b", "fixed_detector")).passed
    for changed in (replace(joint, correlation_evidence=None),
                    replace(joint, correlation_evidence=correlation_evidence(status="assumed")),
                    replace(joint, correlation_evidence=correlation_evidence(source_id="")),
                    replace(joint, correlation_evidence=correlation_evidence(dataset_ids=("FWM-target",))),
                    replace(joint, correlation_evidence=correlation_evidence(status="target_fitted")),
                    replace(joint, parameters=(replace(joint.parameters[0], status="assumed"), joint.parameters[1])),
                    replace(joint, parameters=(replace(joint.parameters[0], status="target_fitted"), joint.parameters[1])),
                    replace(joint, parameters=(replace(joint.parameters[0], dataset_ids=("FWM-target",)), joint.parameters[1]))):
        audit = changed.audit_evidence(target_dataset_ids=("FWM-target",))
        assert not audit.passed and audit.reasons
        assert "no experimental validation" in audit.scope
        # A conditional diagnostic does not relabel its incomplete metadata.
        samples = sample_inputs(changed, 2, seed=0)
        assert not samples.model.audit_evidence(target_dataset_ids=("FWM-target",)).passed


def test_missing_numeric_uncertainty_and_duplicate_ids_never_become_exact_defaults():
    p = model().parameters[0]
    for parameters in ((), (replace(p, uncertainty=None),), (replace(p, value=None),), (p, p)):
        with pytest.raises(ValueError):
            JointInputUncertainty(parameters, np.eye(len(parameters)))


def test_roundoff_eigenvalue_clipping_retains_supplied_correlation_and_diagnostics():
    supplied = np.array([[1., 1.+1e-14], [1.+1e-14, 1.]])
    joint = JointInputUncertainty(model().parameters, supplied)
    np.testing.assert_array_equal(joint.correlation, supplied)
    assert joint.minimum_correlation_eigenvalue < 0
    assert 0 < joint.clipped_negative_eigenvalue_magnitude < joint.correlation_psd_tolerance
    assert joint.factorization_correlation_max_error > 0
    assert joint.factorization_correlation_max_error < joint.correlation_psd_tolerance
    supplied[0, 1] = supplied[1, 0] = 1. + 1e-10
    with pytest.raises(ValueError, match="positive semidefinite"):
        JointInputUncertainty(model().parameters, supplied)


@pytest.mark.parametrize("bad_output", [np.nan, np.inf, 1j, [1., 2.], [[1.]]])
def test_invalid_nonlinear_output_fails_with_original_point_and_no_resampling(bad_output):
    joint = JointInputUncertainty((parameter("x", 0., 1.),), [[1.]])
    expected = sample_inputs(joint, 5, seed=7)
    seen = []
    def callback(x):
        seen.append(x)
        return x[0]**2 if len(seen) == 1 else bad_output
    with pytest.raises(InputEvaluationError) as caught:
        propagate_ensemble(joint, callback, count=5, seed=7,
                           output_ids=("square",), output_units=("1",))
    assert len(seen) == 2
    assert "sample 1, seed 7" in caught.value.context
    np.testing.assert_array_equal(caught.value.input_values, expected.values[1])
    assert caught.value.__cause__ is not None


def test_physical_boundary_failure_is_reported_and_not_clipped_or_rejected():
    joint = JointInputUncertainty((parameter("positive_power", .01, 1., "W"),), [[1.]])
    expected = sample_inputs(joint, 10, seed=3)
    failure_index = int(np.flatnonzero(expected.values[:, 0] <= 0)[0])
    seen = []
    def positive_power(x):
        seen.append(x)
        if x[0] <= 0:
            raise ValueError("power must be positive")
        return np.log(x[0])
    with pytest.raises(InputEvaluationError, match="power must be positive") as caught:
        propagate_ensemble(joint, positive_power, count=10, seed=3,
                           output_ids=("log_power",), output_units=("1",))
    assert len(seen) == failure_index+1
    assert f"sample {failure_index}" in caught.value.context
    with pytest.raises(InputEvaluationError, match="central difference - positive_power"):
        propagate_first_order(joint, positive_power, output_ids=("log_power",),
                              output_units=("1",), steps=[.1])


def test_arrays_are_owned_immutable_and_callback_cannot_change_joint_values():
    corr = np.eye(2)
    joint = JointInputUncertainty(model().parameters, corr)
    corr[:] = 5
    np.testing.assert_array_equal(joint.correlation, np.eye(2))
    for values in (joint.correlation, joint.values, joint.covariance,
                   sample_inputs(joint, 2, seed=1).values):
        with pytest.raises(ValueError):
            values.setflags(write=True)
    with pytest.raises(FrozenInstanceError):
        joint.parameters = ()
    def mutating(x):
        x[0] = 0
        return 1.
    with pytest.raises(InputEvaluationError, match="read-only"):
        propagate_first_order(joint, mutating, output_ids=("value",), output_units=("1",),
                              jacobian=[[0., 0.]])


def test_explicit_seed_count_steps_and_output_units_are_required():
    joint = model()
    for seed in (None, True, -1, 1.5):
        with pytest.raises(ValueError, match="seed"):
            sample_inputs(joint, 2, seed=seed)
    for count in (0, True, 2.5):
        with pytest.raises(ValueError, match="count"):
            sample_inputs(joint, count, seed=1)
    with pytest.raises(ValueError, match="count"):
        propagate_ensemble(joint, lambda x: x[0], count=1, seed=1,
                           output_ids=("a",), output_units=("1",))
    for steps in ([0., 1.], [-1., 1.], [1.], [1e-300, 1e-300]):
        with pytest.raises(ValueError):
            propagate_first_order(joint, lambda x: x[0], output_ids=("a",),
                                  output_units=("1",), steps=steps)
    with pytest.raises(ValueError, match="exactly one"):
        propagate_first_order(joint, lambda x: x[0], output_ids=("a",), output_units=("1",))
    with pytest.raises(ValueError, match="explicit unit"):
        propagate_first_order(joint, lambda x: x[0], output_ids=("a",),
                              output_units=(), jacobian=[[1., 0.]])
