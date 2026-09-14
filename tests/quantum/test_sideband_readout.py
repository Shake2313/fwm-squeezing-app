"""Optical sideband, finite wavepacket and physical current/SQL normalization checks."""

from dataclasses import replace

import numpy as np
import pytest

from analysis.grand_challenge.reference.direct_readout import direct_current_psd
from gabes import constants as c
from gabes.fwm_quantum.field import reduced_field_pair
from gabes.fwm_quantum.model import reduced_pump_noise
from gabes.quantum.channels import GaussianChannel, covariance_uncertainty_minimum
from gabes.quantum.contracts import AnalysisFrequencyAxis, GeneratorFrequencyAxis, OpticalDetunings
from gabes.quantum.readout import DetectorResponse, coherent_carrier_output, intensity_difference_spectrum
from gabes.quantum.sidebands import SIDEBAND_MODES, TopHatBand, average_spectral_channels, sideband_channels
from gabes.quantum.traveling import NambuTransfer, constant_segment
from gabes.schemes import fwm


def _ideal(gain, *, phase=.73):
    axis = AnalysisFrequencyAxis.from_hz([-4e6, -1e6, 0., 1e6, 4e6])
    t = np.array([[np.sqrt(gain), np.sqrt(gain-1)*np.exp(1j*phase)],
                  [np.sqrt(gain-1)*np.exp(-1j*phase), np.sqrt(gain)]])
    values = np.broadcast_to(t, (5, 2, 2))
    zero = np.zeros_like(values)
    main = NambuTransfer(GeneratorFrequencyAxis(-2e10+axis.omega_rad_s, "main"), [1, -1],
                          ("probe", "conjugate"), values, zero, zero)
    companion = NambuTransfer(GeneratorFrequencyAxis(2e10+axis.omega_rad_s, "companion"), [-1, 1],
                               main.mode_labels, values[::-1].conj(), zero, zero)
    return main, companion, axis


def _detector(axis, eta=(1., 1.), balance=1., h=None, electronics=None):
    nf = len(axis.omega_rad_s)
    return DetectorResponse(axis, eta, np.ones((nf, 2), complex) if h is None else h, balance,
                             np.zeros(nf) if electronics is None else electronics, "declared test detector")


def _read(state, beta, detector):
    return intensity_difference_spectrum(state.vacuum_covariances(), beta, detector,
                                         mode_labels=SIDEBAND_MODES, analysis_axis=state.analysis_axis)


@pytest.mark.parametrize("gain", [1., 1.7, 4., 15.])
def test_bright_seed_ideal_intensity_difference_and_equal_detection_loss(gain):
    main, companion, axis = _ideal(gain)
    state = sideband_channels(main, companion, axis)
    flux = 7e11
    beta = coherent_carrier_output(main, axis, [np.sqrt(flux)*np.exp(.41j), 0.])
    np.testing.assert_allclose(np.abs(beta)**2, flux*np.array([gain, gain-1]), atol=.01, rtol=3e-15)
    for eta in (1., .83, .2):
        result = _read(state, beta, _detector(state.analysis_axis, (eta, eta)))
        expected = 1-eta+eta/(2*gain-1)
        np.testing.assert_allclose(result.quantum_ratio, expected, rtol=2e-12, atol=3e-14)
        np.testing.assert_allclose(result.sql_psd_A2_Hz,
                                   2*c.ELEMENTARY_CHARGE**2*eta*flux*(2*gain-1), rtol=3e-15)
        assert result.minimum_covariance_uncertainty_eigenvalue > -1e-10


def test_unequal_loss_balance_and_complex_detector_phase_have_analytic_noise():
    gain, flux, balance = 3.6, 8e12, .81
    main, companion, axis = _ideal(gain)
    state = sideband_channels(main, companion, axis)
    beta = coherent_carrier_output(main, axis, [np.sqrt(flux), 0.])
    eta = np.array([.73, .91])
    h = np.array([[.7+.2j, .4-.3j], [1.4-.1j, .8+.7j]])
    result = _read(state, beta, _detector(state.analysis_axis, eta, balance, h=h))
    p, q = eta
    hp, hc = h.T
    denominator = p*gain*abs(hp)**2+balance**2*q*(gain-1)*abs(hc)**2
    numerator = (p*gain*(1+2*p*(gain-1))*abs(hp)**2
                 +balance**2*q*(gain-1)*(1+2*q*(gain-1))*abs(hc)**2
                 -4*balance*p*q*gain*(gain-1)*(hp*hc.conj()).real)
    np.testing.assert_allclose(result.quantum_ratio, numerator/denominator, rtol=4e-14)


def test_coherent_sql_uses_actual_means_balance_response_and_electronics_units():
    main, companion, axis = _ideal(1.)
    state = sideband_channels(main, companion, axis)
    beta = np.array([np.sqrt(7e10)*np.exp(.8j), np.sqrt(3e11)*np.exp(-.2j)])
    h = np.array([[2.-3j, .7+2j], [.3-.1j, 4.+.2j]])
    eta, balance = np.array([.65, .81]), .4
    electronics = np.array([2e-26, 5e-27])
    result = _read(state, beta, _detector(state.analysis_axis, eta, balance, h, electronics))
    sql = 2*c.ELEMENTARY_CHARGE**2*(eta[0]*abs(beta[0])**2*abs(h[:, 0])**2
          +balance**2*eta[1]*abs(beta[1])**2*abs(h[:, 1])**2)
    np.testing.assert_allclose(result.quantum_ratio, 1., atol=5e-16)
    np.testing.assert_allclose(result.sql_psd_A2_Hz, sql, rtol=3e-15)
    np.testing.assert_allclose(result.total_ratio, 1+electronics/sql, rtol=3e-15)


def test_mode_projection_retains_orthogonal_vacuum_when_transfer_varies_in_band():
    theta = .67
    rotations = [np.array([[np.cos(t), -np.sin(t)], [np.sin(t), np.cos(t)]])
                 for t in (-theta, theta)]
    channels = [GaussianChannel(x, np.zeros((2, 2)), ("a",), ("a",), "dispersive phase") for x in rotations]
    filtered = average_spectral_channels(channels, [.5, .5], source="two equal spectral weights")
    np.testing.assert_allclose(filtered.transfer, np.eye(2)*np.cos(theta), atol=1e-16)
    np.testing.assert_allclose(filtered.added_covariance, np.eye(2)*np.sin(theta)**2/2, rtol=3e-15)
    np.testing.assert_allclose(filtered.apply_covariance(np.eye(2)/2), np.eye(2)/2, atol=1e-16)
    assert not GaussianChannel(filtered.transfer, np.zeros((2, 2)), ("a",), ("a",), "missing leakage").audit().passed
    squeezed = np.diag([np.exp(-1.1), np.exp(1.1)])/2
    np.testing.assert_allclose(filtered.apply_covariance(squeezed),
                               np.cos(theta)**2*squeezed+np.sin(theta)**2*np.eye(2)/2, atol=1e-16)


def test_top_hat_normalization_and_frequency_dependent_channel_convergence():
    results = []
    for order in (4, 8, 16):
        band = TopHatBand(1e6, 6e5, order)
        axis, weights = band.quadrature()
        assert weights.sum() == pytest.approx(1., abs=3e-16)
        assert (weights*axis.frequency_hz).sum() == pytest.approx(1e6)
        # Equal average of a variable phase has analytic sinc attenuation.
        rotations = []
        for freq in axis.frequency_hz:
            phase = 2*(freq-band.center_hz)/band.bandwidth_hz
            x = np.array([[np.cos(phase), -np.sin(phase)], [np.sin(phase), np.cos(phase)]])
            rotations.append(GaussianChannel(x, np.zeros((2, 2)), ("a",), ("a",), "spectral phase"))
        ch = average_spectral_channels(rotations, weights, source="normalized top-hat Gauss integral")
        results.append(ch.transfer[0, 0])
        np.testing.assert_allclose(ch.apply_covariance(np.eye(2)/2), np.eye(2)/2, atol=5e-16)
    assert abs(results[1]-np.sin(1)) < abs(results[0]-np.sin(1))*1e-5
    assert results[-1] == pytest.approx(np.sin(1), abs=3e-16)


@pytest.fixture
def reduced():
    optical = OpticalDetunings(2*np.pi*.9e9, -2*np.pi*8e6)
    atom = reduced_pump_noise(fwm.rabi_freq(.6, fwm.W_PUMP), optical.one_photon_rad_s,
                               transit_rate_s_inverse=c.GAMMA_GG)
    axis = AnalysisFrequencyAxis.from_hz([-4e6, -1e6, 0., 1e6, 4e6])
    pair = reduced_field_pair(atom, optical, axis, number_density_m3=1e18, uniform_area_m2=1.2e-7,
                              optical_omega_rad_s=[c.OMEGA_D1-c.OMEGA_HF, c.OMEGA_D1+c.OMEGA_HF],
                              effective_dipole_C_m=c.DIPOLE_D1/np.sqrt(12))
    return constant_segment(pair.main, .0125), constant_segment(pair.companion, .0125), axis


def test_microscopic_readout_agrees_with_direct_ordered_nambu_contraction(reduced):
    main, companion, axis = reduced
    state = sideband_channels(main, companion, axis)
    beta = coherent_carrier_output(main, axis, [np.sqrt(3e12)*np.exp(.57j), 0.])
    h = np.array([[.8+.17j, 1.3-.41j], [1.5-.9j, .3+.01j]])
    detector = _detector(state.analysis_axis, (.73, .92), .88, h)
    result = _read(state, beta, detector)
    greater, lesser = main.vacuum_output()
    direct = direct_current_psd(greater, lesser, [3, 4], [1, 0], beta,
                                 detector.transmissions, h, detector.balance)
    np.testing.assert_allclose(result.quantum_psd_A2_Hz, direct, rtol=3e-12)
    wrong_side = direct_current_psd(greater, lesser, [3, 4], [3, 4], beta,
                                    detector.transmissions, h, detector.balance)
    assert np.max(abs(wrong_side-direct)/direct) > 1e-5
    assert all(c.audit().passed for c in state.channels)
    assert min(covariance_uncertainty_minimum(v) for v in state.vacuum_covariances()) > -1e-10


def test_carrier_uses_rf_zero_and_outputs_creation_coordinate_conjugate(reduced):
    main, _, axis = reduced
    alpha = np.array([3.+4j, .1+.2j])
    output = coherent_carrier_output(main, axis, alpha)
    expected = main.transfer[2]@np.array([alpha[0], alpha[1].conjugate()])
    np.testing.assert_allclose(output, [expected[0], expected[1].conjugate()], rtol=1e-15)
    assert np.linalg.norm(output-main.transfer[0]@alpha) > .01
    with pytest.raises(ValueError, match="RF=0"):
        coherent_carrier_output(main, AnalysisFrequencyAxis(axis.omega_rad_s+123.), alpha)


def test_bad_companion_missing_mirror_uncertainty_and_zero_sql_are_rejected(reduced):
    main, companion, axis = reduced
    bad_t = companion.transfer.copy(); bad_t[:, 0] *= np.exp(.1j)
    # Retain valid channel algebra by rotating its reservoir noises as well.
    phase = np.diag([np.exp(.1j), 1.])
    bad = replace(companion, transfer=bad_t,
                  noise_greater=phase@companion.noise_greater@phase.conj().T,
                  noise_lesser=phase@companion.noise_lesser@phase.conj().T)
    with pytest.raises(ValueError, match="companion conjugation"):
        sideband_channels(main, bad, axis)
    with pytest.raises(ValueError, match="reflected RF"):
        shifted_axis = AnalysisFrequencyAxis(axis.omega_rad_s+20.)
        shifted_main = replace(main, frequency_axis=GeneratorFrequencyAxis(main.frequency_axis.omega_rad_s+20., "main"))
        shifted_comp = replace(companion, frequency_axis=GeneratorFrequencyAxis(companion.frequency_axis.omega_rad_s+20., "companion"))
        sideband_channels(shifted_main, shifted_comp, shifted_axis)
    state = sideband_channels(main, companion, axis)
    detector = _detector(state.analysis_axis)
    with pytest.raises(ValueError, match="SQL"):
        _read(state, [0., 0.], detector)
    with pytest.raises(ValueError, match="uncertainty"):
        intensity_difference_spectrum(np.tile(np.eye(8)*.1, (2, 1, 1)), [1., 1.], detector,
                                       mode_labels=SIDEBAND_MODES, analysis_axis=state.analysis_axis)
    with pytest.raises(ValueError, match="axes"):
        intensity_difference_spectrum(state.vacuum_covariances(), [1., 1.], detector,
                                       mode_labels=SIDEBAND_MODES, analysis_axis=AnalysisFrequencyAxis.from_hz([2e6, 3e6]))
    with pytest.raises(ValueError):
        TopHatBand(1e6, 2e6, 8)
    with pytest.raises(ValueError):
        detector.transmissions.setflags(write=True)
