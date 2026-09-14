"""Persisted field audits retain numerical and physical scope boundaries."""

import json

import pytest

from analysis.grand_challenge.field_noise_audit import main


def test_field_audit_evidence_and_non_overwriting_report(tmp_path):
    path = tmp_path/"field.json"
    assert main(["--output", str(path)]) == 0
    original = path.read_bytes()
    report = json.loads(original)
    assert report["expected_controls_passed"]
    assert report["reduced_field_diffusion_implemented"]
    assert not report["physical_squeezing_prediction"]
    assert not report["experimental_validation"]
    assert len(report["configurations"]) == 3
    for case in report["configurations"]:
        assert case["noise_only_factor_1_over_12_rejected"]
        assert len(case["sectors"]) == 2
        assert case["field_greater_qrt_relative_error"] < 2e-8
        assert case["field_lesser_qrt_relative_error"] < 2e-8
        assert case["passed"]
    with pytest.raises(FileExistsError):
        main(["--output", str(path)])
    assert path.read_bytes() == original
