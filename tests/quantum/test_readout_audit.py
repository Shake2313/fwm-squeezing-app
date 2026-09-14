"""Joint gain/noise artifacts keep conditional results separate from validation."""

import json

import numpy as np
import pytest

from analysis.grand_challenge.readout_audit import main


def test_readout_report_has_joint_results_normalized_filters_and_no_experimental_badge(tmp_path):
    path = tmp_path/"readout.json"
    assert main(["--output", str(path)]) == 0
    original = path.read_bytes()
    report = json.loads(original)
    assert report["expected_controls_passed"]
    assert report["conditional_linearized_s_minus_implemented"]
    assert not report["absolute_hot_vapor_prediction"] and not report["experimental_validation"]
    assert len(report["configurations"]) == 3
    for case in report["configurations"]:
        assert case["passed"]
        assert len(case["source_spectrum"]["quantum_ratio"]) == 40
        assert case["direct_nambu_psd_maximum_relative_error"] < 1e-9
        assert case["probe_coherent_power_gain"] == pytest.approx(
            case["output_coherent_powers_W"][0]/report["operating_inputs"]["seed_power_W"])
    np.testing.assert_allclose(report["configurations"][0]["source_spectrum"]["quantum_ratio"], 1., atol=1e-12)
    for band in report["normalized_band_convergence"]:
        assert band["channel_audit"]["passed"]
        assert band["unit_gain_real_bandpass_variance_A2"] == pytest.approx(
            band["band_average_quantum_psd_A2_Hz"]*band["bandwidth_hz"], rel=1e-14, abs=0.)
    with pytest.raises(FileExistsError):
        main(["--output", str(path)])
    assert path.read_bytes() == original
