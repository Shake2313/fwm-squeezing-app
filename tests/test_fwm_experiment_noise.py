"""Experimental effective-medium closure and observable noise behavior."""

from copy import deepcopy

import numpy as np
import pytest

from gabes.fwm_experiment_noise import (
    _cached_parameters, effective_parameters, experimental_jasperse_readout,
)
from gabes.jasperse import jasperse_distributed_noise
from gabes.passive_readout import passive_l1_readout


def raw_from_transfer(matrix):
    return {
        "G_s": np.abs(matrix[..., 0, 0]) ** 2,
        "G_c": np.abs(matrix[..., 1, 0]) ** 2,
        "input_seed_power_W": 3e-6,
        "eta": 1.0,
    }


@pytest.mark.parametrize("gain", [1., 1.01, 2., 10., 1000.])
def test_lossless_limit_matches_analytic_bright_twin_beam_noise(gain):
    result = experimental_jasperse_readout(
        {"G_s": np.array([gain]), "G_c": np.array([gain - 1.])}, {})
    assert result["NF"][0] == pytest.approx(1. / (2. * gain - 1.), rel=1e-9)
    assert result["effective_squeeze"][0] == pytest.approx(np.arccosh(np.sqrt(gain)))
    assert result["probe_transmission"][0] == pytest.approx(1.)
    assert result["conjugate_transmission"][0] == pytest.approx(1.)


def test_conjugate_reservoir_produces_finite_nf_where_passive_l1_is_invalid():
    reference = jasperse_distributed_noise(1., 1., .8)
    raw = raw_from_transfer(reference["transfer"])
    assert not passive_l1_readout(reference["transfer"])["passive_completion_valid"]
    result = experimental_jasperse_readout(
        raw, {"after_cell_eff_pct": 71., "qe_pct": 100.})
    assert result["NF"] == pytest.approx(.542363992958049, abs=2e-12)
    assert result["effective_squeeze"] == pytest.approx(1., abs=1e-11)
    assert result["conjugate_transmission"] == pytest.approx(.8, abs=1e-11)
    assert result["fit_relative_residual"] < 1e-11


def test_both_arm_losses_are_recovered_with_the_explicit_common_loss_closure():
    # Independent forward medium fixes the otherwise unidentified parameter.
    reference = jasperse_distributed_noise(1.2, .8, .35, seed_fano=1.)
    raw = raw_from_transfer(reference["transfer"])
    result = experimental_jasperse_readout(raw, {"jasperse_loss_pct": 20.})
    assert result["effective_squeeze"] == pytest.approx(1.2, abs=2e-11)
    assert result["probe_transmission"] == pytest.approx(.8, abs=2e-11)
    assert result["conjugate_transmission"] == pytest.approx(.35, abs=2e-11)
    assert result["NF"] == pytest.approx(reference["NF"], abs=2e-11)


def test_pure_attenuation_does_not_invent_squeezing():
    result = experimental_jasperse_readout({"G_s": np.array([.5]), "G_c": np.array([0.])}, {})
    assert result["NF"][0] == pytest.approx(1.)
    assert result["effective_squeeze"][0] == 0.


def test_conjugate_brighter_than_probe_has_finite_distributed_noise():
    reference = jasperse_distributed_noise(3., .1, 1.)
    raw = raw_from_transfer(reference["transfer"])
    assert raw["G_c"] > raw["G_s"]
    result = experimental_jasperse_readout(raw, {})
    assert np.isfinite(result["NF"]) and result["NF"] > 0.
    assert result["NF"] == pytest.approx(reference["NF"], rel=1e-10)
    assert result["fit_relative_residual"] < 1e-10


def test_readout_only_controls_reuse_inference_and_preserve_input():
    raw = {"G_s": np.array([2.3, 8., 1.2]), "G_c": np.array([1.2, 7.1, .1]),
           "input_seed_power_W": 3e-6}
    saved = deepcopy(raw)
    _cached_parameters.cache_clear()
    baseline = experimental_jasperse_readout(raw, {})
    before = _cached_parameters.cache_info()
    changed = experimental_jasperse_readout(
        raw, {"seed_classical_rin_scaled": 100., "after_cell_eff_pct": 80., "qe_pct": 90.})
    after = _cached_parameters.cache_info()
    assert after.hits == before.hits + 1
    np.testing.assert_array_equal(changed["effective_squeeze"], baseline["effective_squeeze"])
    assert not np.array_equal(changed["NF"], baseline["NF"])
    assert changed["eta_after_cell_QE"] == pytest.approx(.72)
    for key in ("G_s", "G_c"):
        np.testing.assert_array_equal(raw[key], saved[key])


def test_seed_excess_and_external_loss_have_expected_monotonic_effects():
    raw = {"G_s": np.array([10.]), "G_c": np.array([9.]), "input_seed_power_W": 3e-6}
    baseline = experimental_jasperse_readout(raw, {})
    excess = experimental_jasperse_readout(raw, {"seed_classical_rin_scaled": 10.})
    lossy = experimental_jasperse_readout(raw, {"after_cell_eff_pct": 50., "qe_pct": 100.})
    assert excess["NF"][0] > baseline["NF"][0]
    assert lossy["NF"][0] == pytest.approx(.5 + .5 * baseline["NF"][0])


def test_empirical_gains_are_used_as_a_declared_effective_approximation():
    raw = {"G_s": np.array([2.]), "G_c": np.array([1.])}
    gains = {"G_s": np.array([20.]), "G_c": np.array([19.]), "applied": True,
             "gain_quantity": "pump-off reference"}
    result = experimental_jasperse_readout(raw, {}, gains=gains)
    assert result["NF"][0] == pytest.approx(1 / 39.)
    assert result["empirical_gain_calibration_applied"]
    assert result["gain_basis"] == "pump-off reference"
    assert "effective" in result["scope"]
    assert raw["G_s"][0] == 2.


def test_scalar_seed_inputs_expand_to_the_spectrum_for_sidebar_interpolation():
    raw = {"G_s": np.array([2., 3., 4.]), "G_c": np.array([1., 2., 3.]),
           "input_seed_power_W": 3e-6, "omega_probe_rad_s": 2.37e15}
    result = experimental_jasperse_readout(raw, {})
    for name in ("shot_rin_per_Hz", "seed_fano", "seed_power_W"):
        assert result["seed_noise"][name].shape == (3,)


@pytest.mark.parametrize("settings,eta", [
    ({"detection_eff_pct": 80., "loss_pct": 5.5, "qe_pct": 92.}, .8),
    ({"detection_eff_pct": 80., "loss_pct": 10., "qe_pct": 92.}, .9*.92),
    ({"detection_eff_pct": 80., "qe_pct": 85.}, .945*.85),
])
def test_direct_efficiency_alias_preserves_explicit_legacy_overrides(settings, eta):
    result = experimental_jasperse_readout({"G_s": np.array([2.]), "G_c": np.array([1.])}, settings)
    assert result["eta_after_cell_QE"] == pytest.approx(eta)


def test_near_unity_probe_and_extreme_conjugate_tails_remain_finite():
    result = experimental_jasperse_readout(
        {"G_s": np.array([1., .999999, 1.000001]),
         "G_c": np.array([0., 1e-100, 1e-12])}, {})
    assert np.all(np.isfinite(result["NF"]))
    assert np.all(result["NF"] > 0.)


@pytest.mark.parametrize("raw,params", [
    ({"G_s": np.array([0.]), "G_c": np.array([0.])}, {}),
    ({"G_s": np.array([2.]), "G_c": np.array([1.])}, {"after_cell_eff_pct": 0.}),
    ({"G_s": np.array([2.]), "G_c": np.array([1.])}, {"qe_pct": 0.}),
])
def test_dark_detection_has_no_normalized_nf(raw, params):
    result = experimental_jasperse_readout(raw, params)
    assert np.all(np.isnan(result["NF"]))
    assert np.all(result["status"] == "NO_DETECTED_BRIGHT_OUTPUT")


def test_unrepresentable_zero_conjugate_amplifier_is_marked_as_approximate():
    result = experimental_jasperse_readout({"G_s": np.array([2.]), "G_c": np.array([0.])}, {})
    assert np.isfinite(result["NF"][0])
    assert result["inference_approximate"][0]
    assert result["status"][0] == "EFFECTIVE_JASPERSE_GAIN_APPROXIMATION"
    assert result["fit_relative_residual"][0] > 0.


def test_large_loss_fallback_keeps_the_reservoir_instead_of_inventing_a_twin_beam():
    # Actual changed-condition spectrum: the former lossless fallback replaced
    # conjugate .0143 with13.609 and incorrectly predictedNF .1614.
    probe, conjugate = 14.609208115154768, .014328900148894413
    result = experimental_jasperse_readout(
        {"G_s": np.array([probe]), "G_c": np.array([conjugate])}, {})
    assert result["inference_approximate"][0]
    assert result["G_s"][0] == pytest.approx(probe, rel=1e-11)
    assert result["G_c"][0] < .1
    assert result["fit_relative_residual"][0] < 4.
    assert result["conjugate_transmission"][0] < 1e-250
    assert result["NF"][0] > 1.


@pytest.mark.parametrize("probe,conjugate", [(.2, 1e-5), (.1, 1e-7)])
def test_strong_common_loss_fallback_preserves_dominant_gain_and_reservoir_noise(probe, conjugate):
    result = experimental_jasperse_readout(
        {"G_s": np.array([probe]), "G_c": np.array([conjugate])},
        {"jasperse_loss_pct": 95.})
    assert result["inference_approximate"][0]
    assert result["G_s"][0] == pytest.approx(probe, rel=1e-11)
    # A lossless fallback after common loss generated probe-.05 conjugate;
    # bounded distributed loss substantially improves that discrepancy.
    old_conjugate = probe - .05
    old_residual = abs(old_conjugate - conjugate) / conjugate
    assert result["fit_relative_residual"][0] < old_residual / 100.
    assert result["NF"][0] > 1.


def test_bounded_loss_fallback_matches_the_same_forward_distributed_medium():
    result = experimental_jasperse_readout({"G_s": np.array([2.]), "G_c": np.array([0.])}, {})
    expected = jasperse_distributed_noise(
        result["effective_squeeze"], result["probe_transmission"], result["conjugate_transmission"])
    np.testing.assert_allclose(result["NF"], expected["NF"])
    assert np.any(expected["reservoir_quadrature_covariance"] > 0.)


@pytest.mark.parametrize("probe,conjugate,loss", [
    (np.nan, 1., 0.), (-1., 1., 0.), (1., -1., 0.),
    (1., 1., 100.), (1., 1., -.1), (1., 1., np.nan),
])
def test_invalid_effective_inputs_raise(probe, conjugate, loss):
    with pytest.raises(ValueError):
        effective_parameters(probe, conjugate, common_loss_pct=loss)
