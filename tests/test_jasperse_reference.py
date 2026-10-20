"""Independent Jasperse controls, including the printed 2010/2011 differences."""

import numpy as np
import pytest

from gabes.constants import HBAR
from gabes.jasperse import jasperse_distributed_noise, seed_noise_components
from gabes.passive_readout import passive_l1_readout


def test_seed_shot_noise_units_and_power_scaling():
    power = np.array([1e-6, 2e-6, 4e-6])
    omega = 2.0 * np.pi * 377.1e12
    result = seed_noise_components(power, omega, classical_rin_per_Hz=1e-13)
    energy = HBAR * omega
    np.testing.assert_allclose(result["shot_rin_per_Hz"], 2 * energy / power)
    np.testing.assert_allclose(result["shot_power_psd_W2_per_Hz"], 2 * energy * power)
    np.testing.assert_allclose(result["classical_power_psd_W2_per_Hz"], 1e-13 * power ** 2)
    np.testing.assert_allclose(result["seed_fano"], 1 + 1e-13 * power / (2 * energy))
    # Integrating the one-sided shot PSD over 1/(2*tau) recovers counting variance.
    tau = 3.7e-5
    photons = power * tau / energy
    integrated = result["shot_power_psd_W2_per_Hz"] / (2 * tau)
    np.testing.assert_allclose(integrated, (energy / tau) ** 2 * photons)


def test_seed_zero_classical_and_broadcast():
    result = seed_noise_components(np.ones((2, 1)) * 1e-6,
                                   np.ones(3) * 2e15)
    assert result["seed_fano"].shape == (2, 3)
    np.testing.assert_array_equal(result["seed_fano"], 1.0)
    np.testing.assert_array_equal(result["classical_power_psd_W2_per_Hz"], 0.0)


@pytest.mark.parametrize("power,omega,classical", [
    (0, 1, 0), (-1, 1, 0), (np.nan, 1, 0), (np.inf, 1, 0),
    (1, 0, 0), (1, -1, 0), (1, np.nan, 0), (1, np.inf, 0),
    (1, 1, -1), (1, 1, np.nan), (1, 1, np.inf),
])
def test_seed_invalid_physical_inputs(power, omega, classical):
    with pytest.raises(ValueError):
        seed_noise_components(power, omega, classical_rin_per_Hz=classical)


def test_oe_eq3_and_classical_seed_limit():
    gain = np.array([1., 2., 10., 100.])
    squeeze = np.arccosh(np.sqrt(gain))
    np.testing.assert_allclose(jasperse_distributed_noise(squeeze)["NF"],
                               1 / (2 * gain - 1), rtol=1e-10)
    np.testing.assert_allclose(jasperse_distributed_noise(squeeze, seed_fano=3)["NF"],
                               3 / (2 * gain - 1), rtol=1e-10)


def test_oe_eq5_unequal_external_losses():
    gain = np.array([1.1, 2., 10., 100.])[:, None]
    eta_p = np.array([.8, .9, .7])
    eta_c = np.array([.8, .7, .9])
    expected = 1 + 2 * (gain - 1) * (
        gain * (eta_p - eta_c) ** 2 - eta_c ** 2) / (
            gain * eta_p + (gain - 1) * eta_c)
    result = jasperse_distributed_noise(
        np.arccosh(np.sqrt(gain)), probe_efficiency=eta_p,
        conjugate_efficiency=eta_c)
    np.testing.assert_allclose(result["NF"], expected, rtol=1e-10, atol=2e-12)


def _paper_special_case(squeeze, transmission, conjugate_loss=False):
    """2011 OE (16)/(17), independently transcribed after rendering PDF pages."""
    ta, tb = (1., transmission) if conjugate_loss else (transmission, 1.)
    xi = .25 * np.sqrt(16 * squeeze ** 2 + np.log(ta / tb) ** 2)
    chi = np.arctanh(np.log(ta / tb) / (4 * xi))
    den = np.cosh(2 * xi + chi)
    if conjugate_loss:
        return (1 - .71 * 2 * squeeze * np.cosh(xi + chi) ** 2 / (xi * den)
                + .71 * np.sqrt(tb) * squeeze
                * (4 * squeeze - np.log(tb) * np.sinh(2 * xi + chi)) ** 2
                / (8 * xi ** 3 * den))
    return (1 - .71 * 2 * squeeze * np.sinh(xi) ** 2 / (xi * den)
            + .71 * np.sqrt(ta) * squeeze * np.log(ta) ** 2
            * np.sinh(xi) ** 4 / (2 * xi ** 3 * den))


@pytest.mark.parametrize("conjugate_loss", [False, True])
def test_oe_eq16_eq17_against_exact_reservoir_integral(conjugate_loss):
    squeeze = np.array([.1, .3, 1., 2., 3.])[:, None]
    transmission = np.array([.01, .1, .5, .9, 1.])
    expected = _paper_special_case(squeeze, transmission, conjugate_loss)
    ta, tb = (1., transmission) if conjugate_loss else (transmission, 1.)
    result = jasperse_distributed_noise(
        squeeze, ta, tb, probe_efficiency=.71, conjugate_efficiency=.71)
    np.testing.assert_allclose(result["NF"], expected, rtol=1e-10, atol=2e-12)


def test_pure_probe_absorption_reproduces_current_passive_scope():
    squeeze = np.array([0., .3, 1., 2.])[:, None, None]
    ta = np.array([.1, .6, 1.])[None, :, None]
    fano = np.array([1., 3.])
    reference = jasperse_distributed_noise(
        squeeze, ta, seed_fano=fano, probe_efficiency=.71,
        conjugate_efficiency=.71)
    passive = passive_l1_readout(reference["transfer"], seed_fano=fano, eta=.71)
    assert np.all(passive["passive_completion_valid"])
    np.testing.assert_allclose(reference["NF"], passive["passive_l1_S_det"], atol=1e-12)


def test_conjugate_loss_is_not_annihilation_only_passive_completion():
    result = jasperse_distributed_noise(
        1., 1., .8, probe_efficiency=.71, conjugate_efficiency=.71)
    passive = passive_l1_readout(result["transfer"], eta=.71)
    assert float(result["NF"]) == pytest.approx(.542363992958049, abs=1e-12)
    assert not passive["passive_completion_valid"]
    assert np.isnan(passive["passive_l1_noise_dB"])
    assert result["NF"] > passive["passive_l1_S_det"]


def test_no_mixing_coherent_seed_and_no_bright_detection():
    result = jasperse_distributed_noise(0., .3, .7,
                                       probe_efficiency=.6, conjugate_efficiency=.2)
    assert result["NF"] == pytest.approx(1.)
    absent = jasperse_distributed_noise(1., probe_efficiency=0.,
                                       conjugate_efficiency=0.)
    assert np.isnan(absent["NF"])


@pytest.mark.parametrize("kwargs", [
    {"squeeze_parameter": -1}, {"squeeze_parameter": np.nan},
    {"probe_transmission": 0}, {"probe_transmission": 1.01},
    {"conjugate_transmission": -.1}, {"conjugate_transmission": np.inf},
    {"seed_fano": 0}, {"seed_fano": np.nan},
    {"probe_efficiency": -.1}, {"probe_efficiency": np.inf},
    {"conjugate_efficiency": 1.01}, {"conjugate_efficiency": np.nan},
])
def test_reference_rejects_undefined_or_nonpassive_inputs(kwargs):
    with pytest.raises(ValueError):
        jasperse_distributed_noise(**({"squeeze_parameter": 1.} | kwargs))
