"""Analytic and independently integrated checks of the effective PSR channel."""

import numpy as np
import pytest
from scipy.integrate import solve_ivp

from gabes.psr_noise import psr_gaussian_readout


@pytest.mark.parametrize("shear", [-7.0, -1.0, 0.0, .03, 2.0, 20.0])
def test_lossless_shear_matches_analytic_squeezing_eigenvalues(shear):
    result = psr_gaussian_readout(shear, 0.0)
    maximum = 1+shear**2/2+abs(shear)*np.sqrt(1+shear**2/4)
    assert result["maximum_variance"] == pytest.approx(maximum)
    assert result["minimum_variance"] == pytest.approx(1/maximum, rel=1e-10)
    assert np.linalg.det(result["covariance"]) == pytest.approx(1.0)
    np.testing.assert_allclose(result["source_added_covariance"], 0.0)


def test_shear_sign_mirrors_phase_without_changing_noise_extrema():
    phase = np.array([0., 13., 35., 90., 141., 180.])
    positive = psr_gaussian_readout(1.7, .4, phase_degrees=phase)
    negative = psr_gaussian_readout(-1.7, .4, phase_degrees=-phase)
    np.testing.assert_allclose(positive["noise_a"], negative["noise_a"])
    assert positive["minimum_variance"] == pytest.approx(negative["minimum_variance"])
    assert (positive["min_angle_deg"]+negative["min_angle_deg"]) % 180 == pytest.approx(0.)


@pytest.mark.parametrize("depth", [0., 1e-14, 1e-5, .001, .1, 5., 100., 1e6])
def test_no_shear_pure_absorption_preserves_vacuum(depth):
    result = psr_gaussian_readout(0., depth)
    np.testing.assert_allclose(result["covariance"], np.eye(2), atol=3e-15)
    np.testing.assert_allclose(result["noise_a_db"], 0., atol=1e-13)
    assert result["channel_audit_passed"]


@pytest.mark.parametrize("shear,depth,excess", [
    (1.2, 0., 0.), (1.2, 1e-7, .3), (-2.1, .7, 0.), (3.0, 4., .8),
])
def test_distributed_channel_matches_independent_covariance_ode(shear, depth, excess):
    drift = np.array([[-depth/2, shear], [0., -depth/2]])
    diffusion = depth*(1+excess)*np.eye(2)/2

    def derivative(_z, flattened):
        covariance = flattened.reshape(2, 2)
        return (drift@covariance+covariance@drift.T+diffusion).ravel()

    solved = solve_ivp(derivative, (0., 1.), (np.eye(2)/2).ravel(),
                       rtol=2e-11, atol=2e-13)
    result = psr_gaussian_readout(shear, depth, excess_noise_per_od=excess)
    assert solved.success
    np.testing.assert_allclose(result["source_covariance"], 2*solved.y[:, -1].reshape(2, 2),
                               rtol=5e-11, atol=5e-12)


def test_distributed_loss_is_not_replaced_by_endpoint_loss():
    shear, depth = 2., 1.1
    actual = psr_gaussian_readout(shear, depth)
    shear_only = np.array([[1., shear], [0., 1.]])
    endpoint = np.exp(-depth)*(shear_only@shear_only.T)+(1-np.exp(-depth))*np.eye(2)
    assert not np.allclose(actual["covariance"], endpoint)
    assert actual["covariance"][0, 1] > endpoint[0, 1]


def test_excess_diffusion_is_positive_and_can_remove_squeezing():
    result = psr_gaussian_readout(1.1, .9, excess_noise_per_od=[0., 1., 10.])
    assert np.all(np.diff(result["minimum_variance"]) > 0)
    assert np.all(np.diff(result["maximum_variance"]) > 0)
    assert result["minimum_variance"][-1] > 1
    for difference in np.diff(result["covariance"], axis=0):
        assert np.linalg.eigvalsh(difference).min() > 0
    assert "optimistically omits" in " ".join(result["assumptions"])


def test_excess_noise_parameter_is_per_absorption_od():
    result = psr_gaussian_readout(0., [.0, .4, 2.], excess_noise_per_od=.7)
    expected = 1+.7*(1-np.exp(-np.array([0., .4, 2.])))
    np.testing.assert_allclose(result["minimum_variance"], expected)
    np.testing.assert_allclose(result["maximum_variance"], expected)


@pytest.mark.parametrize("depth", [1e-16, 1e-4, .001, 800., 1e8, 1e300])
def test_optical_depth_limits_stay_finite(depth):
    result = psr_gaussian_readout(1.7, depth)
    assert np.isfinite(result["covariance"]).all()
    assert np.isfinite(result["noise_a_db"]).all()
    if depth < 1e-10:
        expected = np.array([[1+1.7**2, 1.7], [1.7, 1.]])
        np.testing.assert_allclose(result["covariance"], expected, atol=1e-13)
    if depth > 1e7:
        np.testing.assert_allclose(result["covariance"], np.eye(2), atol=2e-8)


def test_large_depth_and_proportional_shear_retain_finite_local_noise():
    result = psr_gaussian_readout(2e300, 1e300)
    np.testing.assert_allclose(result["covariance"], [[9., 2.], [2., 1.]])
    assert result["channel_audit_passed"]


@pytest.mark.parametrize("eta", [0., .35, 1.])
def test_detection_loss_applied_once(eta):
    result = psr_gaussian_readout(.8, .7, excess_noise_per_od=.2,
                                 detection_efficiency=eta)
    expected = eta*result["source_covariance"]+(1-eta)*np.eye(2)
    np.testing.assert_allclose(result["covariance"], expected, atol=1e-14)
    if eta == 0:
        np.testing.assert_allclose(result["covariance"], np.eye(2))
        np.testing.assert_allclose(result["noise_a_db"], 0., atol=1e-13)


def test_channel_cp_uncertainty_and_bright_stokes_commutator():
    result = psr_gaussian_readout([-.9, 0., 2.], [.1, 3., .8],
                                 excess_noise_per_od=.4, detection_efficiency=.72)
    j = np.array([[0., 1.], [-1., 0.]])
    carrier_photons = 1e6
    stokes_map = np.sqrt(2*carrier_photons)*np.eye(2)
    np.testing.assert_allclose(stokes_map@j@stokes_map.T, 2*carrier_photons*j)
    for covariance, transfer, noise in zip(
            result["covariance"], result["detected_transfer"], result["detected_added_covariance"]):
        cp_matrix = noise+.5j*(j-transfer@j@transfer.T)
        assert np.linalg.eigvalsh(cp_matrix).min() >= -1e-12
        stokes_covariance = carrier_photons*covariance
        assert np.linalg.eigvalsh(stokes_covariance+1j*carrier_photons*j).min() >= -1e-6
        assert np.linalg.det(covariance) >= 1-1e-12
    assert np.all(result["channel_audit_passed"])
    assert np.all(result["uncertainty_minimum"] >= -1e-12)


def test_phase_scan_periodicity_orthogonal_swap_and_minimum_angle():
    base = psr_gaussian_readout(1.8, .4)
    angle = float(base["min_angle_deg"])
    phases = angle+np.array([0., 45., 90., 180., 225., 270.])
    result = psr_gaussian_readout(1.8, .4, phase_degrees=phases)
    np.testing.assert_allclose(result["noise_a"][:3], result["noise_a"][3:])
    assert result["noise_a"][0] == pytest.approx(result["minimum_variance"])
    assert result["noise_a"][2] == pytest.approx(result["maximum_variance"])
    assert result["noise_a"][0] == pytest.approx(result["noise_b"][2])
    assert result["noise_a"][1] == pytest.approx(result["noise_b"][1])
    assert result["noise_a"][1] == pytest.approx(result["crossing_variance"])


def test_broadcast_shapes_and_inputs_preserved():
    shear = np.array([[.0], [.5]])
    original = shear.copy()
    result = psr_gaussian_readout(shear, [0., .2, .8], phase_degrees=[0., 30., 60., 90.])
    assert result["covariance"].shape == (2, 3, 2, 2)
    assert result["minimum_db"].shape == (2, 3)
    assert result["noise_a"].shape == (2, 3, 4)
    np.testing.assert_array_equal(shear, original)
    result["shear"][0, 0] = 9
    np.testing.assert_array_equal(shear, original)


@pytest.mark.parametrize("kwargs", [
    {"shear": np.nan}, {"optical_depth": -1.}, {"excess_noise_per_od": -.1},
    {"detection_efficiency": 1.01}, {"phase_degrees": []}, {"phase_degrees": [np.inf]},
])
def test_invalid_model_inputs_rejected(kwargs):
    values = {"shear": 1., "optical_depth": .1}
    values.update(kwargs)
    with pytest.raises(ValueError):
        psr_gaussian_readout(**values)
