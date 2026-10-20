"""Independent diagonal transport and efficiency-domain controls."""

import numpy as np
import pytest

from gabes.fwm_efficiency import dressed_cell_transmission


def test_independent_beer_lambert_intensity_and_phase_control():
    # A field attenuation amplitude is exp(-alpha*L/2), not exp(-alpha*L).
    alpha_probe, alpha_conjugate, length = 23.0, 71.0, 0.0125
    diagonal = np.array([
        -alpha_probe / 2 + 937j, -alpha_conjugate / 2 - 743j])
    drift = np.diag(diagonal)
    actual = dressed_cell_transmission(drift, length)
    expected = np.abs(np.exp(diagonal * length)) ** 2
    for index, arm in enumerate(("probe", "conjugate")):
        assert actual[f"in_cell_{arm}_transmission"] == pytest.approx(expected[index])
        assert actual[f"in_cell_{arm}_passive_efficiency"] == pytest.approx(expected[index])
        assert actual[f"in_cell_{arm}_net_passive"]
    assert actual["in_cell_probe_signed_alpha_per_m"] == alpha_probe
    assert actual["in_cell_conjugate_signed_alpha_per_m"] == alpha_conjugate


def test_full_mixing_gain_cannot_be_confused_with_cell_transmission():
    # An ideal parametric amplifier has |T00|^2=cosh(qL)^2>1 while each
    # counterfactual diagonal transfer is exactly unity.
    length, q = 0.2, 10.0
    drift = np.array([[0.0, q], [q, 0.0]], dtype=complex)
    actual = dressed_cell_transmission(drift, length)
    assert np.cosh(q * length) ** 2 > 10.0
    for arm in ("probe", "conjugate"):
        assert actual[f"in_cell_{arm}_transmission"] == 1.0
        assert actual[f"in_cell_{arm}_passive_efficiency"] == 1.0


def test_signed_active_diagonal_is_retained_without_false_passive_efficiency():
    drift = np.array([[-12.0, 9.0j], [17.0, 3.0]], dtype=complex)
    actual = dressed_cell_transmission(drift, 0.1)
    assert actual["in_cell_probe_transmission"] == pytest.approx(np.exp(-2.4))
    assert actual["in_cell_conjugate_transmission"] == pytest.approx(np.exp(0.6))
    assert actual["in_cell_conjugate_signed_alpha_per_m"] == -6.0
    assert not actual["in_cell_conjugate_net_passive"]
    assert np.isnan(actual["in_cell_conjugate_passive_efficiency"])
    assert actual["in_cell_conjugate_status"] == (
        "ACTIVE_DIAGONAL_RESPONSE_NO_PASSIVE_EFFICIENCY")


def test_stacked_detunings_keep_arm_domains_independent():
    drift = np.zeros((2, 3, 2, 2), dtype=complex)
    drift[..., 0, 0] = np.array([[-1.0, 2.0, 0.0], [1.0, -2.0, -3.0]])
    drift[..., 1, 1] = -4.0 + 3.0j
    actual = dressed_cell_transmission(drift, 0.2)
    assert actual["in_cell_probe_transmission"].shape == (2, 3)
    np.testing.assert_array_equal(
        actual["in_cell_probe_net_passive"], np.real(drift[..., 0, 0]) <= 0)
    assert np.all(actual["in_cell_conjugate_net_passive"])
    assert np.all(np.isfinite(actual["in_cell_conjugate_passive_efficiency"]))


def test_zero_length_is_transparent_even_with_gain_and_mixing():
    drift = np.array([[1e308, 1e300j], [-1e300j, -1e308]])
    actual = dressed_cell_transmission(drift, 0.0)
    for arm in ("probe", "conjugate"):
        assert actual[f"in_cell_{arm}_transmission"] == 1.0
        assert actual[f"in_cell_{arm}_passive_efficiency"] == 1.0


def test_underflow_and_overflow_preserve_signed_log_transport():
    drift = np.diag([-1e6, 1e6]).astype(complex)
    actual = dressed_cell_transmission(drift, 0.1)
    assert actual["in_cell_probe_transmission"] == 0.0
    assert actual["in_cell_probe_passive_efficiency"] == 0.0
    assert actual["in_cell_probe_log_transmission"] == -200000.0
    assert np.isposinf(actual["in_cell_conjugate_transmission"])
    assert actual["in_cell_conjugate_log_transmission"] == 200000.0
    assert np.isnan(actual["in_cell_conjugate_passive_efficiency"])


def test_overflowed_log_cannot_pass_the_roundoff_guard():
    drift = np.diag([-1e308, 1e308]).astype(complex)
    actual = dressed_cell_transmission(drift, 1e308)
    assert np.isneginf(actual["in_cell_probe_log_transmission"])
    assert actual["in_cell_probe_passive_efficiency"] == 0.0
    assert np.isposinf(actual["in_cell_conjugate_log_transmission"])
    assert not actual["in_cell_conjugate_net_passive"]
    assert np.isnan(actual["in_cell_conjugate_passive_efficiency"])


def test_roundoff_normalization_is_explicit_and_raw_gain_is_retained():
    epsilon = np.finfo(float).eps
    drift = np.diag([epsilon, 1e-10]).astype(complex)
    actual = dressed_cell_transmission(drift, 1.0)
    assert actual["in_cell_probe_transmission"] > 1.0
    assert actual["in_cell_probe_passive_efficiency"] == 1.0
    assert actual["in_cell_probe_roundoff_normalized"]
    assert actual["in_cell_probe_status"] == "UNITY_WITH_ROUNDOFF_TOLERANCE"
    assert not actual["in_cell_conjugate_net_passive"]


@pytest.mark.parametrize("drift", (np.zeros(4), np.zeros((2, 3)),
                                  np.full((2, 2), np.nan),
                                  np.full((2, 2), np.inf)))
def test_invalid_drift_rejected(drift):
    with pytest.raises(ValueError, match="maxwell_drift"):
        dressed_cell_transmission(drift, 0.01)


@pytest.mark.parametrize("length", (-1.0, np.nan, np.inf, [0.1]))
def test_invalid_length_rejected(length):
    with pytest.raises(ValueError, match="L must be"):
        dressed_cell_transmission(np.zeros((2, 2)), length)
