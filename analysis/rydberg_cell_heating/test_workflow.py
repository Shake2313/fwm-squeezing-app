from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from analysis.rydberg_cell_heating import workflow as workflow_module
from analysis.rydberg_cell_heating.adapters import validate_status
from analysis.rydberg_cell_heating.workflow import load_clean_csv, run_analysis
from gabes import species


def _config(tmp_path: Path) -> Path:
    payload = {
        "schema_version": 1,
        "analysis_id": "test-cell-heating",
        "output_dir": "generated",
        "model": {
            "temperatures_c": [20.0, 40.0],
            "views": ["EIT"],
            "base_params": {"doppler": "off", "if_khz": 40.0},
            "include_spectra_in_results": True,
        },
        "inputs": [
            {
                "id": "missing_raw",
                "kind": "eit_spectrum",
                "path": None,
                "required": False,
                "status": "PENDING",
            }
        ],
        "source_notes": [{"id": "assumption", "status": "ASSUMED"}],
        "helpers": [],
        "figures": {"formats": ["png"], "dpi": 80},
    }
    path = tmp_path / "config.json"
    path.write_text(json.dumps(payload), encoding="utf-8")
    return path


def test_status_vocabulary_is_strict():
    assert validate_status(" measured ") == "MEASURED"
    with pytest.raises(ValueError, match="unknown evidence status"):
        validate_status("probably")


def test_serialized_paths_are_repo_relative_or_safe_external_ids(tmp_path: Path):
    internal = workflow_module._serialized_path(Path(workflow_module.__file__))
    assert internal == "analysis/rydberg_cell_heating/workflow.py"
    assert "\\" not in internal

    external = tmp_path / "operator trace.csv"
    external_id = workflow_module._serialized_path(external)
    assert external_id.startswith(
        "external://operator_trace.csv#path-sha256="
    )
    assert str(tmp_path) not in external_id


def test_clean_sensitivity_csv(tmp_path: Path):
    path = tmp_path / "sensitivity.csv"
    path.write_text(
        "temperature_c,sensitivity_nv_cm_sqrt_hz,uncertainty_nv_cm_sqrt_hz,status\n"
        "40,6.68,2.55,PPT\n",
        encoding="utf-8",
    )
    rows = load_clean_csv(path, "sensitivity", "PENDING")
    assert rows == [{
        "temperature_c": 40.0,
        "sensitivity_nv_cm_sqrt_hz": 6.68,
        "uncertainty_nv_cm_sqrt_hz": 2.55,
        "status": "PPT",
    }]


def test_smoke_run_emits_status_tagged_artifacts(tmp_path: Path):
    artifacts = run_analysis(_config(tmp_path))
    assert artifacts["results"].is_file()
    assert artifacts["manifest"].is_file()
    assert artifacts["tex_macros"].is_file()
    assert (artifacts["output_dir"] / "eit_temperature_spectra.png").is_file()
    assert (artifacts["output_dir"] / "temperature_metrics.png").is_file()

    results = json.loads(artifacts["results"].read_text(encoding="utf-8"))
    assert results["inputs"][0]["status"] == "PENDING"
    assert results["capabilities"]["static_rydberg_spectrum"]["status"] == "PREDICTED"
    assert results["capabilities"]["absolute_rf_sensitivity"]["status"] == "PENDING"
    assert len(results["model"]["sweeps"]["EIT"]) == 2
    assert results["model"]["optima"]["max_spectral_slope"] is not None


def test_shared_eit_loader_and_temperature_metadata_are_preserved(tmp_path: Path):
    trace = tmp_path / "eit_40c.csv"
    trace.write_text(
        "Probe detuning [MHz],Photodiode voltage [V]\n"
        "-1,0.2\n0,0.8\n1,0.2\n",
        encoding="utf-8",
    )
    config = json.loads(_config(tmp_path).read_text(encoding="utf-8"))
    config["inputs"] = [{
        "id": "eit_40c",
        "kind": "eit_spectrum",
        "path": trace.name,
        "temperature_c": 40.0,
        "required": True,
        "status": "MEASURED",
        "loader": {"source_name": str(trace)},
    }]
    path = tmp_path / "trace_config.json"
    path.write_text(json.dumps(config), encoding="utf-8")
    artifacts = run_analysis(path, make_figures=False)
    results = json.loads(artifacts["results"].read_text(encoding="utf-8"))
    dataset = results["inputs"][0]
    assert dataset["temperature_c"] == 40.0
    assert dataset["parser_provenance"]["sha256"]
    assert [row["detuning_mhz"] for row in dataset["rows"]] == [-1.0, 0.0, 1.0]
    assert dataset["path"].startswith(
        "external://eit_40c.csv#path-sha256="
    )
    assert dataset["parser_provenance"]["source_name"] == dataset["path"]

    manifest = json.loads(artifacts["manifest"].read_text(encoding="utf-8"))
    assert manifest["config"]["path"].startswith(
        "external://trace_config.json#path-sha256="
    )
    assert manifest["inputs"][0]["path"] == dataset["path"]
    code_paths = {entry["path"] for entry in manifest["code"]}
    assert "analysis/rydberg_cell_heating/workflow.py" in code_paths
    assert "analysis/rydberg_cell_heating/adapters.py" in code_paths
    assert all("\\" not in path for path in code_paths)
    assert all(
        artifact["path"].startswith("external://")
        for artifact in manifest["artifacts"]
    )


def test_git_metadata_is_captured_before_generated_writes(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    config_path = _config(tmp_path)
    generated = tmp_path / "generated"
    expected = {"commit": "test-head", "branch": "test", "dirty": False}
    calls = []

    def capture_before_write():
        calls.append(True)
        assert not generated.exists()
        return expected

    monkeypatch.setattr(workflow_module, "_git_metadata", capture_before_write)
    artifacts = run_analysis(config_path, make_figures=False)
    manifest = json.loads(artifacts["manifest"].read_text(encoding="utf-8"))
    assert calls == [True]
    assert manifest["git"] == expected


def test_calibrated_finite_if_chain_can_emit_absolute_sensitivity(tmp_path: Path):
    config = json.loads(_config(tmp_path).read_text(encoding="utf-8"))
    config["model"]["temperatures_c"] = [40.0]
    config["model"]["views"] = ["AT electrometry"]
    config["electrometry"] = {
        "enabled": True,
        "status": "ASSUMED",
        "transition_dipole_c_m": 2.0e-27,
        "angular_factor": 0.5,
        "field_amplitude_convention": "peak",
        "quantum_efficiency": 0.8,
        "signal_probe_power_scale": 1.0,
        "electronic_noise_current_asd_a_per_sqrt_hz": 1.0e-12,
        "probe_detuning_search_points": 11,
    }
    path = tmp_path / "electrometry_config.json"
    path.write_text(json.dumps(config), encoding="utf-8")
    artifacts = run_analysis(path, make_figures=False)
    results = json.loads(artifacts["results"].read_text(encoding="utf-8"))
    point = results["model"]["sweeps"]["AT electrometry"][0]
    assert point["electrometry"]["status"] == "PREDICTED"
    assert point["electrometry"]["field_sensitivity_nv_cm_sqrt_hz"] > 0.0
    assert results["capabilities"]["two_tone_superheterodyne"]["available"]
    assert results["capabilities"]["absolute_rf_sensitivity"]["available"]


def test_cold_spot_temperature_is_passed_into_solver(tmp_path: Path):
    temperature_log = tmp_path / "temperature.csv"
    temperature_log.write_text(
        "setpoint_c,effective_vapor_temp_c,cold_spot_c,status\n"
        "50,50,40,MEASURED\n",
        encoding="utf-8",
    )
    config = json.loads(_config(tmp_path).read_text(encoding="utf-8"))
    config["model"]["temperatures_c"] = [50.0]
    config["model"]["axial_profile"] = {
        "enabled": True,
        "status": "ASSUMED",
        "left_offset_c": 0.0,
        "right_offset_c": 0.0,
        "points": 3,
        "density_mode": "cold_spot_limited",
    }
    config["inputs"] = [{
        "id": "temperature",
        "kind": "temperature_log",
        "path": temperature_log.name,
        "required": True,
        "status": "MEASURED",
    }]
    path = tmp_path / "cold_spot_config.json"
    path.write_text(json.dumps(config), encoding="utf-8")
    artifacts = run_analysis(path, make_figures=False)
    results = json.loads(artifacts["results"].read_text(encoding="utf-8"))
    point = results["model"]["sweeps"]["EIT"][0]
    expected = (
        species.number_density(species.RB85, 40.0 + 273.15)
        * (40.0 + 273.15) / (50.0 + 273.15)
    )
    assert point["parameters"]["temperature_model"] == "Separated"
    assert point["parameters"]["cold_spot_temp_c"] == 40.0
    assert np.isclose(point["metrics"]["number_density_m3"]["value"], expected)
    assert point["axial_profile"]["cold_spot_temp_c"] == 40.0
    assert np.allclose(point["axial_profile"]["density_m3"], expected)


def test_native_electrometry_axial_and_sam_paths_are_reported(tmp_path: Path):
    config = json.loads(_config(tmp_path).read_text(encoding="utf-8"))
    config["model"]["temperatures_c"] = [40.0]
    config["model"]["views"] = ["AT electrometry"]
    config["model"]["effective_atom_number"] = {
        "enabled": True,
        "status": "ASSUMED",
        "participation_fraction": 0.5,
        "overlap_efficiency": 0.8,
    }
    config["model"]["axial_profile"] = {
        "enabled": True,
        "status": "ASSUMED",
        "z_fraction": [0.0, 0.5, 1.0],
        "temperature_offsets_c": [-2.0, 1.0, 2.0],
        "density_mode": "cold_spot_limited",
    }
    config["electrometry"] = {
        "enabled": True,
        "mode": "scheme",
        "status": "ASSUMED",
    }
    config["sam_calibration"] = {
        "enabled": True,
        "status": "ASSUMED",
        "source_power_dbm": -40.0,
        "antenna_gain_dbi": 10.0,
        "distance_m": 0.3,
        "source_power_std_db": 0.2,
        "distance_std_m": 0.002,
        "amplitude_convention": "rms",
        "frequency_hz": 37.0e9,
        "antenna_max_dimension_m": 0.1,
    }
    path = tmp_path / "integrated_config.json"
    path.write_text(json.dumps(config), encoding="utf-8")
    artifacts = run_analysis(path, make_figures=False)
    results = json.loads(artifacts["results"].read_text(encoding="utf-8"))
    point = results["model"]["sweeps"]["AT electrometry"][0]
    assert point["electrometry"]["field_sensitivity_nv_cm_sqrt_hz"] > 0.0
    assert point["axial_profile"]["beer_lambert"]["available"]
    assert point["parameters"]["cold_spot_temp_c"] == 38.0
    assert point["axial_profile"]["cold_spot_temp_c"] == 38.0
    assert not point["axial_profile"]["beer_lambert"][
        "used_in_finite_if_sensitivity"]
    assert point["effective_atom_number"]["axial_density_integrated"]
    assert results["sam_calibration"]["standard_uncertainty_v_m"] > 0.0
    assert point["sam_comparison"]["available"]
    assert point["sam_comparison"]["frequency_match"]
    assert point["sam_comparison"]["model_to_sam_ratio"] > 0.0
    assert results["capabilities"]["finite_if_superheterodyne"]["available"]
    assert results["capabilities"]["axial_beer_lambert"]["available"]
    assert results["capabilities"]["sam_field_calibration"]["available"]
    assert results["capabilities"]["sam_at_comparison"]["available"]
    assert results["model"]["optima"]["min_total_field_sensitivity"] is not None


def test_axial_capability_requires_emitted_segmented_spectrum(tmp_path: Path):
    config = json.loads(_config(tmp_path).read_text(encoding="utf-8"))
    config["model"]["temperatures_c"] = [40.0]
    config["model"]["include_spectra_in_results"] = False
    config["model"]["axial_profile"] = {
        "enabled": True,
        "status": "ASSUMED",
        "left_offset_c": 0.0,
        "right_offset_c": 0.0,
    }
    path = tmp_path / "no_axial_spectrum_config.json"
    path.write_text(json.dumps(config), encoding="utf-8")
    artifacts = run_analysis(path, make_figures=False)
    results = json.loads(artifacts["results"].read_text(encoding="utf-8"))
    point = results["model"]["sweeps"]["EIT"][0]
    assert not point["axial_profile"]["beer_lambert"]["available"]
    assert not results["capabilities"]["axial_beer_lambert"]["available"]


@pytest.mark.parametrize("mode,include_spectrum", [
    ("scheme", True), ("scheme", False), ("configured", False),
])
def test_opacity_temperature_scenarios_reuse_nominal_states_and_report_optima(
        tmp_path, monkeypatch, mode, include_spectrum):
    from gabes import schemes

    config = json.loads(_config(tmp_path).read_text(encoding="utf-8"))
    config["model"]["views"] = ["EIT", "AT electrometry"]
    config["model"]["include_spectra_in_results"] = include_spectrum
    config["model"]["base_params"].update(
        effective_opacity_participation=.001,
        effective_opacity_relative_uncertainty=1.0)
    config["model"]["axial_profile"] = {
        "enabled": True, "status": "ASSUMED", "points": 3,
        "left_offset_c": 0.0, "right_offset_c": 0.0,
    }
    config["electrometry"] = {
        "enabled": True, "mode": mode, "status": "ASSUMED",
        "transition_dipole_c_m": 2e-27, "angular_factor": .5,
        "field_amplitude_convention": "peak", "quantum_efficiency": .8,
        "signal_probe_power_scale": 1.0, "probe_detuning_search_points": 11,
    }
    path = tmp_path / "opacity.json"
    path.write_text(json.dumps(config), encoding="utf-8")
    scheme = schemes.get("rydberg_eit")
    original_compute = scheme.compute
    solves = []

    def counted(params):
        solves.append(dict(params))
        return original_compute(params)

    monkeypatch.setattr(scheme, "compute", counted)
    artifacts = run_analysis(path, make_figures=False)
    assert len(solves) == 4  # two temperatures x two views; none for scenarios
    result = json.loads(artifacts["results"].read_text(encoding="utf-8"))
    json.dumps(result, allow_nan=False)
    model = result["model"]
    opacity = model["effective_opacity"]
    assert opacity["status"] == "ASSUMED"
    assert opacity["calibration_status"] == "UNCALIBRATED"
    scenarios = model["opacity_scenarios"]
    assert scenarios["status"] == "PREDICTED"
    assert "not measured confidence" in scenarios["scope"]
    assert scenarios["samples"] == [
        {"label": "nominal", "value": .001},
        {"label": "lower", "value": 0.0},
        {"label": "upper", "value": .002},
    ]
    assert scenarios["sweeps"]["nominal"] == model["sweeps"]
    assert scenarios["optima"]["nominal"] == model["optima"]
    assert scenarios["optima"]["lower"]["min_total_field_sensitivity"] is None
    assert scenarios["optima"]["lower"]["max_spectral_slope"] is None
    assert scenarios["optima"]["lower"]["max_if_discriminator"] is None
    assert scenarios["optima"]["upper"]["min_total_field_sensitivity"]["value"] > 0.0
    for view in model["views"]:
        for index, nominal in enumerate(model["sweeps"][view]):
            for label in ("lower", "upper"):
                point = scenarios["sweeps"][label][view][index]
                assert ("spectrum" in point) is include_spectrum
                assert point["effective_opacity"]["value"] == (0.0 if label == "lower" else .002)
            lower = scenarios["sweeps"]["lower"][view][index]
            upper = scenarios["sweeps"]["upper"][view][index]
            assert lower["metrics"]["transmission_at_resonance"]["value"] == 1.0
            assert upper["metrics"]["transmission_at_resonance"]["value"] < nominal[
                "metrics"]["transmission_at_resonance"]["value"]
            if include_spectrum:
                # Uniform axial density preserves the already scaled lumped T;
                # multiplying opacity again here would fail this parity check.
                np.testing.assert_allclose(upper["axial_profile"]["beer_lambert"]["transmission"],
                                           upper["spectrum"]["transmission"], rtol=1e-14)
    manifest = json.loads(artifacts["manifest"].read_text(encoding="utf-8"))
    assert "gabes/rydberg_opacity.py" in {row["path"] for row in manifest["code"]}


def test_zero_opacity_deduplicates_workflow_scenarios(tmp_path):
    config = json.loads(_config(tmp_path).read_text(encoding="utf-8"))
    config["model"]["temperatures_c"] = [40.0]
    config["model"]["base_params"]["effective_opacity_participation"] = 0.0
    path = tmp_path / "zero-opacity.json"
    path.write_text(json.dumps(config), encoding="utf-8")
    artifact = run_analysis(path, make_figures=False)
    result = json.loads(artifact["results"].read_text(encoding="utf-8"))
    scenarios = result["model"]["opacity_scenarios"]
    assert scenarios["samples"] == [{"label": "nominal", "value": 0.0}]
    assert set(scenarios["sweeps"]) == {"nominal"}
    assert all(value is None for value in scenarios["optima"]["nominal"].values())
    assert all(value is None for value in result["model"]["optima"].values())
