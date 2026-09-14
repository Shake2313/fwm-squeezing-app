"""Thermal/paired analytic limits and explicit Fock-space fourth moments."""

import numpy as np
import pytest

from analysis.grand_challenge.reference.fock_photocurrent import two_pair_current_reference
from gabes import constants as c
from gabes.quantum.contracts import AnalysisFrequencyAxis
from gabes.quantum.photocounting import (
    PairedOpticalBins, box_autocorrelation, gaussian_count_covariance,
    quadratic_readout_correction, corrected_gaussian_readout,
)
from gabes.quantum.readout import DetectorResponse, IntensityDifferenceSpectrum


def test_rectangular_thermal_band_has_white_shot_noise_and_triangular_bunching():
    bins = PairedOpticalBins([-1.5, -.5, .5, 1.5], 1., np.tile([.3, .7], (4, 1)), np.zeros(4), "thermal")
    freq = np.array([-5, -1.25, 0., .5, 4., 10.])
    cov = gaussian_count_covariance(bins, AnalysisFrequencyAxis.from_hz(freq))
    for j, n in enumerate((.3, .7)):
        np.testing.assert_allclose(cov[:, j, j], 4*n+np.maximum(4-abs(freq), 0)*n*n)
    np.testing.assert_array_equal(cov[:, 0, 1], 0.)


def test_fractional_complex_bin_overlap_has_no_wraparound_and_keeps_phase():
    x = np.array([1., 2j])
    r = box_autocorrelation(x, 2., [-4., -3., -2., -1., 0., 1., 2., 3., 4.])
    np.testing.assert_allclose(r, [0., -2j, -4j, 5-2j, 10., 5+2j, 4j, 2j, 0.])


@pytest.mark.parametrize("eta", [1., .85, .3])
def test_unseeded_pure_twin_beams_and_detection_loss(eta):
    n = .4
    bins = PairedOpticalBins([-1.5, -.5, .5, 1.5], 1., np.full((4, 2), n), np.full(4, np.sqrt(n*(n+1))), "paired")
    axis = AnalysisFrequencyAxis.from_hz([0., .3, 2., 5.])
    cov = gaussian_count_covariance(bins.attenuate([eta, eta]), axis)
    difference = np.einsum('i,fij,j->f', [1, -1], cov, [1, -1]).real
    np.testing.assert_allclose(difference/(8*n*eta), 1-eta*np.maximum(4-axis.frequency_hz, 0)/4, atol=1e-14)


def test_quadratic_psd_matches_independent_fock_current_with_complex_pair_phase():
    n = np.array([.1, .2])
    phase = np.array([.2, -.7])
    bins = PairedOpticalBins([-.5, .5], 1., np.column_stack([n, n[::-1]]), np.sqrt(n*(n+1))*np.exp(1j*phase), "two pairs")
    cov = gaussian_count_covariance(bins, AnalysisFrequencyAxis.from_hz([1.]))[0]
    ref = two_pair_current_reference(*n, *phase, cutoff=20)
    np.testing.assert_allclose(cov, ref["count_covariance_per_hz"], atol=2e-13, rtol=0)


def test_asymmetric_detector_weights_electronics_and_total_mean_sql():
    axis = AnalysisFrequencyAxis.from_hz([.25, .5])
    n = np.array([.1, .2])
    bins = PairedOpticalBins([-.5, .5], 1., np.column_stack([n, n[::-1]]),
        np.sqrt(n*(n+1))*np.exp(1j*np.array([.2, -.7])), "asymmetric")
    detector = DetectorResponse(axis, [.7, .4], [[1, 1j], [.8+.3j, -.4j]], .6, [1e-39, 2e-39], "known detector")
    correction = quadratic_readout_correction(bins, detector)
    # Independently weight the already validated two-beam count matrix.
    mean = n.sum()*np.array([.7, .4])
    np.testing.assert_allclose(correction.spontaneous_detected_flux_s_inverse, mean)
    weights = detector.current_response*np.array([1, -.6])
    expected_sql = 2*c.ELEMENTARY_CHARGE**2*np.sum(abs(weights)**2*mean, axis=1)
    np.testing.assert_allclose(correction.added_sql_psd_A2_Hz, expected_sql, atol=0, rtol=1e-14)
    bright = IntensityDifferenceSpectrum(axis, expected_sql*10, expected_sql*20,
        detector.electronics_difference_psd_A2_Hz, np.full(2, .5), np.full(2, .5), mean*20, 0., 0., "same fixture")
    full = corrected_gaussian_readout(bright, correction)
    np.testing.assert_allclose(full['quantum_ratio'], (bright.quantum_psd_A2_Hz+correction.quadratic_psd_A2_Hz)/(expected_sql*21))
    np.testing.assert_allclose(full['total_ratio']-full['quantum_ratio'], detector.electronics_difference_psd_A2_Hz/(expected_sql*21))


def test_vacuum_has_no_spontaneous_flux_or_quadratic_counts_and_invalid_moments_fail():
    bins = PairedOpticalBins([-.5, .5], 1., np.zeros((2, 2)), np.zeros(2), "vacuum")
    np.testing.assert_array_equal(gaussian_count_covariance(bins, AnalysisFrequencyAxis.from_hz([0., 10.])), 0.)
    with pytest.raises(ValueError, match="unphysical"):
        PairedOpticalBins([-.5, .5], 1., np.zeros((2, 2)), np.ones(2), "bad")
    with pytest.raises(ValueError, match="uniform"):
        PairedOpticalBins([-1., 1.], 1., np.zeros((2, 2)), np.zeros(2), "bad")
    with pytest.raises(ValueError):
        bins.centers_hz.setflags(write=True)
