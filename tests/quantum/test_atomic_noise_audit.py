"""Reproducible atomic spectrum output with explicit field-noise limits."""

import json

import pytest

from analysis.grand_challenge.atomic_noise_audit import main


def test_atomic_noise_report_keeps_input_frames_and_physical_claims_separate(tmp_path):
    path = tmp_path/"atomic.json"
    assert main(["--output", str(path)]) == 0
    original = path.read_bytes()
    report = json.loads(original)
    assert report["expected_controls_passed"]
    assert report["atomic_ordered_diffusion_implemented"]
    assert not report["field_diffusion_implemented"]
    assert not report["physical_squeezing_prediction"]
    assert not report["experimental_validation"]
    assert len(report["configurations"]) == 3
    assert min(report["operating_inputs"]["rf_analysis_hz"]) > 0
    for case in report["configurations"]:
        assert case["passed"]
        shifted = case["scans"][1]
        assert max(shifted["generator_omega_rad_s"]) < 0
        assert shifted["spectral_audit"]["passed"]
    with pytest.raises(FileExistsError):
        main(["--output", str(path)])
    assert path.read_bytes() == original
