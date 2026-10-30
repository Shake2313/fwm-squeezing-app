"""Physical limits and interpolation of the effective FWM quadrature readout."""

from copy import deepcopy

import numpy as np
import pytest

from gabes.fwm_experiment_noise import experimental_jasperse_readout
from gabes.fwm_quadrature import experimental_quadrature_phase_scan
from gabes.jasperse import jasperse_distributed_noise


def _readout(squeeze, *, ta=1.0, tb=1.0, seed_fano=1.0, eta=1.0):
    squeeze = np.atleast_1d(squeeze)
    result = jasperse_distributed_noise(squeeze, ta, tb, seed_fano=seed_fano)
    return dict(result, effective_squeeze=squeeze, probe_transmission=ta,
                conjugate_transmission=tb, eta_after_cell_QE=eta)


@pytest.mark.parametrize("gain", [1.0, 1.01, 2.0, 10.0, 100.0])
def test_lossless_gain_matches_analytic_joint_quadrature_extrema(gain):
    readout = experimental_jasperse_readout(
        {"G_s": np.array([gain]), "G_c": np.array([gain-1])}, {})
    result = experimental_quadrature_phase_scan(readout, [0.0], 0.0)
    r = np.arccosh(np.sqrt(gain))
    assert result["available"]
    assert result["minimum_variance"] == pytest.approx(np.exp(-2*r), rel=1e-9)
    assert result["maximum_variance"] == pytest.approx(np.exp(2*r), rel=1e-12)
    assert result["noise_a"][0] == pytest.approx(result["minimum_variance"])
    assert result["noise_b"][0] == pytest.approx(result["maximum_variance"])
    # Direct bright photon-number difference uses unequal mean-field weights.
    if gain > 1:
        assert result["minimum_variance"] != pytest.approx(readout["NF"][0])


@pytest.mark.parametrize("eta", [0.0, 0.3, 0.85, 1.0])
def test_external_loss_is_applied_once_and_vacuum_limit_is_sql(eta):
    source = experimental_quadrature_phase_scan(_readout(0.9, ta=.7, tb=.4), [0.0], 0.0)
    detected = experimental_quadrature_phase_scan(
        _readout(0.9, ta=.7, tb=.4, eta=eta), [0.0], 0.0)
    for name in ("noise_a", "noise_b"):
        np.testing.assert_allclose(detected[name], 1-eta+eta*source[name], rtol=2e-13)
    vacuum = experimental_quadrature_phase_scan(
        _readout(0.0, ta=.3, tb=.8, eta=eta), [0.0], 0.0)
    np.testing.assert_allclose(vacuum["noise_a"], 1.0, atol=1e-14)
    np.testing.assert_allclose(vacuum["noise_b_db"], 0.0, atol=1e-13)


def test_amplitude_excess_noise_preserved_without_inventing_seed_phase_noise():
    clean = experimental_quadrature_phase_scan(_readout(.6), [0.0], 0.0)
    noisy = experimental_quadrature_phase_scan(_readout(.6, seed_fano=100), [0.0], 0.0)
    assert noisy["noise_a"][0] > clean["noise_a"][0]
    assert noisy["noise_b"][0] == pytest.approx(clean["noise_b"][0])
    assert noisy["minimum_variance"] > 1.0  # No forced squeezing or reciprocal anti-squeezing.
    np.testing.assert_allclose(noisy["covariance"][1::2, 1::2],
                               clean["covariance"][1::2, 1::2])
    assert "shot-noise" in noisy["seed_phase_assumption"]


def test_mixed_loss_noisy_seed_and_interpolated_covariance_obey_uncertainty():
    readout = _readout([.2, 1.3], ta=[.9, .2], tb=[.4, .8], seed_fano=[2., 7.], eta=.77)
    result = experimental_quadrature_phase_scan(readout, [-1., 2.], .5)
    j = np.kron(np.eye(2), [[0., 1.], [-1., 0.]])
    assert result["available"]
    assert np.linalg.eigvalsh(result["covariance"]+1j*j)[0] >= -1e-12
    assert result["minimum_variance"]*result["maximum_variance"] >= 1-1e-12
    endpoints = [experimental_quadrature_phase_scan(readout, [-1., 2.], x)
                 for x in (-1., 2.)]
    np.testing.assert_allclose(result["covariance"],
                               .5*(endpoints[0]["covariance"]+endpoints[1]["covariance"]))
    assert result["minimum_variance"] == pytest.approx(
        .5*(endpoints[0]["noise_a"][0]+endpoints[1]["noise_a"][0]))


@pytest.mark.parametrize("selected", [-1.001, 2.001])
def test_outside_scan_is_unavailable_without_endpoint_clamping(selected):
    result = experimental_quadrature_phase_scan(_readout([.2, .9]), [-1., 2.], selected)
    assert not result["available"]
    assert result["status"] == "OPERATING_POINT_OUTSIDE_SCAN"
    assert np.isnan(result["noise_a"]).all()
    assert np.isnan(result["covariance"]).all()


@pytest.mark.parametrize("selected,index", [(-30.0, 0), (30.0, 1)])
def test_nominal_endpoint_survives_ghz_to_mhz_roundoff(selected, index):
    axis = np.array([-30.0+2e-13, 30.0-2e-13])
    readout = _readout([.2, .9])
    result = experimental_quadrature_phase_scan(readout, axis, selected)
    reference = experimental_quadrature_phase_scan(readout, axis, axis[index])
    assert result["available"]
    assert result["selected_delta_mhz"] == selected
    np.testing.assert_array_equal(result["covariance"], reference["covariance"])


@pytest.mark.parametrize("selected", [-30.0-2e-9, 30.0+2e-9])
def test_true_out_of_support_point_is_not_snapped(selected):
    result = experimental_quadrature_phase_scan(_readout([.2, .9]), [-30., 30.], selected)
    assert not result["available"]
    assert result["status"] == "OPERATING_POINT_OUTSIDE_SCAN"


def test_pi_periodicity_orthogonal_swap_and_crossings_in_linear_units():
    phase = np.array([0., 12., 45., 90., 102., 135., 180., 192., 225.])
    result = experimental_quadrature_phase_scan(_readout(.8), [0.0], 0.0, phase_degrees=phase)
    np.testing.assert_allclose(result["noise_a"][:3], result["noise_a"][6:])
    np.testing.assert_allclose(result["noise_a"][:3], result["noise_b"][3:6])
    for index in (2, 5, 8):
        assert result["noise_a"][index] == pytest.approx(result["noise_b"][index])
        assert result["noise_a"][index] == pytest.approx(result["crossing_variance"])
    assert result["crossing_variance"] > 1  # Crossings are not forced onto SQL.
    np.testing.assert_allclose(result["noise_a_db"], 10*np.log10(result["noise_a"]))


def test_input_is_unchanged_and_default_phase_range_is_closed_half_turn():
    readout = _readout([.2, .7])
    original = deepcopy(readout)
    result = experimental_quadrature_phase_scan(readout, [0., 1.], .4)
    for key in readout:
        if isinstance(readout[key], np.ndarray):
            np.testing.assert_array_equal(readout[key], original[key])
    np.testing.assert_array_equal(result["phase_degrees"][[0, -1]], [0., 180.])


def test_unphysical_covariance_is_rejected_without_repair():
    readout = _readout(.5)
    readout["quadrature_covariance"] = np.eye(2)[None]*.01
    result = experimental_quadrature_phase_scan(readout, [0.], 0.)
    assert not result["available"]
    assert result["status"] == "EFFECTIVE_QUADRATURE_UNCERTAINTY_FAILED"


@pytest.mark.parametrize("axis", [[], [1., 0.], [0., 0.], [0., np.nan]])
def test_invalid_scan_axis_rejected(axis):
    with pytest.raises(ValueError, match="scan_axis"):
        experimental_quadrature_phase_scan(_readout(.5), axis, 0.)
