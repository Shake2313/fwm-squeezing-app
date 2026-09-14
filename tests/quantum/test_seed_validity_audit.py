"""Moment extraction parity and immutable artifact preflight."""

import numpy as np
import pytest

from analysis.grand_challenge.seed_validity_audit import main, _field
from analysis.grand_challenge.normalization_audit import conditional_inputs
from gabes.fwm_quantum.model import reduced_pump_noise
from gabes.fwm_quantum.normalization import reduced_dipoles
from gabes.quantum.contracts import AnalysisFrequencyAxis
from gabes.quantum.photocounting import paired_optical_bins
from gabes.quantum.sidebands import sideband_channels
from gabes.quantum.traveling import constant_segment


def test_broadband_ordered_occupations_agree_with_four_sideband_quadratures():
    inputs = conditional_inputs()
    d = reduced_dipoles('uniform-zeeman-rms')
    atom = reduced_pump_noise(d.pump_rabi_rad_s(inputs.pump_power_W, inputs.pump_waist_m), inputs.one_photon_rad_s,
        transit_rate_s_inverse=inputs.transit_rate_s_inverse, transition_scales=d.transition_scales)
    axis = AnalysisFrequencyAxis.from_hz([-3e5, -1e5, 1e5, 3e5])
    pair = _field(atom, inputs, d, axis)
    main, companion = constant_segment(pair.main, inputs.length_m), constant_segment(pair.companion, inputs.length_m)
    bins = paired_optical_bins(main, companion, axis, bin_width_hz=2e5, source='test band')
    cov = sideband_channels(main, companion, axis).vacuum_covariances()
    occupations = (cov[:, ::2, ::2].diagonal(axis1=-2, axis2=-1)+cov[:, 1::2, 1::2].diagonal(axis1=-2, axis2=-1)-1)/2
    np.testing.assert_allclose(occupations[:, 0], bins.occupations[2:, 0], atol=1e-13)
    np.testing.assert_allclose(occupations[:, 1], bins.occupations[1::-1, 1], atol=1e-13)
    np.testing.assert_allclose(occupations[:, 2], bins.occupations[1::-1, 0], atol=1e-13)
    np.testing.assert_allclose(occupations[:, 3], bins.occupations[2:, 1], atol=1e-13)


def test_existing_report_is_not_overwritten(tmp_path):
    path = tmp_path/'result.json'
    path.write_bytes(b'previous run')
    with pytest.raises(FileExistsError):
        main(['--output', str(path)])
    assert path.read_bytes() == b'previous run'
