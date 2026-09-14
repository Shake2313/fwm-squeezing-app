"""Reproducible research entry point and non-overwriting report contract."""

import json

import pytest

from analysis.grand_challenge.foundation_audit import main


def test_foundation_command_emits_scoped_report_and_preserves_existing_artifacts(tmp_path):
    output = tmp_path / "run" / "report.json"
    assert main(["--output", str(output)]) == 0
    original = output.read_bytes()
    report = json.loads(original)
    assert report["expected_controls_passed"]
    assert not report["physical_squeezing_prediction"]
    assert not report["microscopic_diffusion_implemented"]
    assert not report["experimental_validation"]
    assert not report["generator_audits"]["legacy_default_counterexample"]["passed"]
    assert report["generator_audits"]["legacy_combined_121C"]["passed"]
    assert len(report["source_sha256"]) >= 8
    with pytest.raises(FileExistsError):
        main(["--output", str(output)])
    assert output.read_bytes() == original
