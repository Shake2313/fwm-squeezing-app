"""No-fit provenance and physically distinguishable channel checks."""

from dataclasses import replace

import numpy as np
import pytest

from gabes.quantum.channels import (
    GaussianChannel, canonical_commutator, compose_channels,
    covariance_uncertainty_minimum,
)
from gabes.quantum.contracts import (
    AnalysisFrequencyAxis, OpticalDetunings, ParameterEvidence,
    audit_independent_inputs,
)


MODES = ("probe", "conjugate")


def _two_mode_squeezer(r):
    c, s = np.cosh(r), np.sinh(r)
    x = np.array([[c, 0, s, 0], [0, c, 0, -s],
                  [s, 0, c, 0], [0, -s, 0, c]])
    return GaussianChannel(x, np.zeros((4, 4)), MODES, MODES,
                           "ideal two-mode symplectic fixture")


def test_signed_rf_axis_converts_once_and_is_distinct_from_optical_detunings():
    hz = np.array([-1e6, 0., 1e6])
    axis = AnalysisFrequencyAxis.from_hz(hz)
    optical = OpticalDetunings(2*np.pi*0.9e9, -2*np.pi*8e6)
    np.testing.assert_allclose(axis.omega_rad_s, [-2*np.pi*1e6, 0, 2*np.pi*1e6])
    np.testing.assert_allclose(axis.frequency_hz, hz)
    hz[:] = 42
    assert axis.omega_rad_s[1] == 0
    assert optical.two_photon_rad_s != axis.omega_rad_s[0]
    with pytest.raises((ValueError, TypeError)):
        AnalysisFrequencyAxis(optical)
    for values in ([1., 1.], [np.nan], [[1., 2.]], [1j]):
        with pytest.raises(ValueError):
            AnalysisFrequencyAxis(values)


def _independent_input():
    return ParameterEvidence(
        "pump_waist", 530e-6, "m", "independent", "beam camera run 17",
        3e-6, "Gaussian fit to independently measured beam profile",
        ("beam-profile-17",), ("beam intensity profile",), "cell center, 600 mW")


def test_independent_calibration_fit_is_allowed_but_target_fit_and_leakage_are_not():
    evidence = _independent_input()
    kwargs = dict(required_ids=("pump_waist",), target_dataset_ids=("FWM-target",))
    assert audit_independent_inputs([evidence], **kwargs).passed
    for changed in (
        replace(evidence, status="target_fitted"),
        replace(evidence, status="assumed"),
        replace(evidence, dataset_ids=("FWM-target",)),
        replace(evidence, uncertainty=None),
        replace(evidence, source_id=""),
    ):
        result = audit_independent_inputs([changed], **kwargs)
        assert not result.passed
        assert result.reasons


def test_omitted_and_duplicate_input_records_cannot_pass_the_no_fit_gate():
    evidence = _independent_input()
    assert not audit_independent_inputs([], required_ids=("pump_waist",)).passed
    assert not audit_independent_inputs([evidence], required_ids=()).passed
    assert not audit_independent_inputs([evidence, evidence],
                                         required_ids=("pump_waist",)).passed
    assert not audit_independent_inputs([evidence],
                                         required_ids=("pump_waist", "density")).passed


def test_unequal_passive_losses_preserve_coherent_vacuum_and_compose_as_transmissions():
    first = GaussianChannel.vacuum_attenuator(MODES, [0.4, 0.8], source="fixture A")
    second = GaussianChannel.vacuum_attenuator(MODES, [0.6, 0.3], source="fixture B")
    composed = compose_channels(first, second)
    direct = GaussianChannel.vacuum_attenuator(MODES, [0.24, 0.24], source="direct")
    assert first.audit().passed and composed.audit().passed
    # Independent vacuum-reservoir dilation for this analytic loss fixture.
    j = canonical_commutator(2)
    b = np.diag(np.sqrt(np.repeat([0.6, 0.2], 2)))
    np.testing.assert_allclose(first.transfer@j@first.transfer.T+b@j@b.T, j, atol=1e-15)
    np.testing.assert_allclose(composed.transfer, direct.transfer, atol=1e-15)
    np.testing.assert_allclose(composed.added_covariance, direct.added_covariance, atol=1e-15)
    np.testing.assert_allclose(composed.apply_covariance(np.eye(4)/2), np.eye(4)/2, atol=1e-15)


def test_ideal_squeezing_and_vacuum_loss_have_analytic_epr_variance():
    r, eta = 0.9, 0.72
    squeezer = _two_mode_squeezer(r)
    j = canonical_commutator(2)
    np.testing.assert_allclose(squeezer.transfer@j@squeezer.transfer.T, j, atol=1e-15)
    assert squeezer.audit().passed
    loss = GaussianChannel.vacuum_attenuator(MODES, [eta, eta], source="fixture")
    channel = compose_channels(squeezer, loss)
    v = channel.apply_covariance(np.eye(4)/2)
    difference = np.array([1., 0, -1., 0])
    assert difference@v@difference == pytest.approx(1-eta+eta*np.exp(-2*r), abs=1e-14)
    assert covariance_uncertainty_minimum(v) >= -1e-14
    np.testing.assert_allclose(v, loss.apply_covariance(
        squeezer.apply_covariance(np.eye(4)/2)), atol=1e-15)
    # Reversing gain and loss changes the added noise even with commuting X.
    reverse = compose_channels(loss, squeezer).apply_covariance(np.eye(4)/2)
    assert not np.allclose(v, reverse)


def test_positive_noise_alone_does_not_make_a_physical_channel():
    bad = GaussianChannel(np.sqrt(0.5)*np.eye(2), 0.01*np.eye(2),
                          ("probe",), ("probe",), "missing loss vacuum fixture")
    audit = bad.audit()
    assert audit.minimum_noise_eigenvalue > 0
    assert audit.minimum_cp_eigenvalue < -0.2
    assert not audit.passed
    with pytest.raises(ValueError, match="complete positivity"):
        bad.apply_covariance(np.eye(2)/2)
    # Positive covariance also needs the uncertainty relation.
    with pytest.raises(ValueError, match="uncertainty"):
        GaussianChannel.identity(("probe",)).apply_covariance(np.eye(2)*0.1)


def test_rectangular_collection_preserves_subsystem_covariance():
    collect = GaussianChannel(np.eye(4)[:2], np.zeros((2, 2)), MODES,
                              ("probe",), "trace out uncollected conjugate")
    assert collect.audit().passed
    full = _two_mode_squeezer(0.7).apply_covariance(np.eye(4)/2)
    np.testing.assert_allclose(collect.apply_covariance(full), full[:2, :2], atol=1e-15)


def test_mode_order_and_complex_rf_basis_cannot_be_silently_mixed():
    with pytest.raises(ValueError, match="mode order"):
        compose_channels(GaussianChannel.identity(MODES),
                         GaussianChannel.identity(tuple(reversed(MODES))))
    with pytest.raises(ValueError, match="real coordinates"):
        GaussianChannel(np.eye(2)*1j, np.eye(2)/2, ("a",), ("a",), "Nambu fixture")
    with pytest.raises(ValueError, match="symmetric"):
        GaussianChannel(np.eye(2), [[1, 0.2], [0, 1]], ("a",), ("a",), "bad Y")


def test_channel_owns_immutable_matrices_and_has_no_implicit_mean_field_map():
    x = np.eye(2)
    channel = GaussianChannel(x, np.zeros((2, 2)), ("a",), ("a",), "identity")
    x[:] = 0
    np.testing.assert_array_equal(channel.transfer, np.eye(2))
    with pytest.raises(ValueError):
        channel.transfer.setflags(write=True)
