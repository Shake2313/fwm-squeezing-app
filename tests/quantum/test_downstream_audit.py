"""The independent downstream audit reaches the existing atomic readout."""

import json

import numpy as np
import pytest

from analysis.grand_challenge.downstream_audit import main


def test_parallel_report_preserves_censoring_correlations_and_conditional_claims(tmp_path):
    path = tmp_path/"downstream.json"
    assert main(["--output", str(path)]) == 0
    original = path.read_bytes()
    report = json.loads(original)
    assert report["expected_controls_passed"]
    assert report["sources_unchanged_during_run"]
    assert not report["absolute_hot_vapor_prediction"]
    assert not report["experimental_validation"]
    assert report["singular_input_controls"]["amplified_contrast_variance"] == 0.
    bands = report["analytic_bands"]
    assert bands["locally_refined_maximum_edge_error_hz"] < .01
    assert bands["flat_sub_sql"]["bands"][0]["width_hz"] is None
    assert bands["flat_sub_sql"]["bands"][0]["observed_span_hz"] == pytest.approx(3.9e6)
    uncertainty = report["correlated_band_uncertainty"]
    assert uncertainty["deterministic_replay"]
    assert not uncertainty["evidence_audit"]["passed"]
    np.testing.assert_allclose(uncertainty["linear_covariance"], uncertainty["exact_covariance"], rtol=.003)
    atomic = report["conditional_atomic_uncertainty"]
    assert len(atomic["output_standard_uncertainties"]) == 5
    assert all(x > 0 for x in atomic["output_standard_uncertainties"])
    assert atomic["sampled_band_analysis"]["bands"][0]["width_hz"] is None
    assert not atomic["evidence_audit"]["passed"]
    # All outputs depend on the same two input draws; don't reduce to five
    # independent error bars or discard pump power/waist correlations.
    covariance = np.array(atomic["output_covariance"])
    assert covariance[0, 1] != 0
    assert np.linalg.matrix_rank(covariance, tol=1e-10*np.linalg.norm(covariance)) <= 2
    assert not np.allclose(atomic["output_standard_uncertainties"],
                          atomic["uncorrelated_standard_uncertainties"], rtol=.01, atol=0)
    with pytest.raises(FileExistsError):
        main(["--output", str(path)])
    assert path.read_bytes() == original


def test_report_and_plot_cannot_share_an_output_path(tmp_path):
    path = tmp_path/"output"
    with pytest.raises(ValueError, match="different paths"):
        main(["--output", str(path), "--plot", str(path)])
    assert not path.exists()
