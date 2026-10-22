"""Live apparatus profiles and portable sessions, without a Bloch solve."""
from dataclasses import fields, replace
import json
from pathlib import Path
import sys

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from sabes import bridge  # noqa: E402
from sabes.apparatus import default_apparatus  # noqa: E402
from sabes.beamline import SetupSettings  # noqa: E402
from sabes.calibration import default_calibration  # noqa: E402
from sabes.detection import DetectionSettings  # noqa: E402
from sabes.presets import (LAB_PROFILE, PAPER_PROFILE, PROFILES,  # noqa: E402
                           LAB_AFTER_CELL_TRANSMISSION, LAB_PUMP_W,
                           LAB_SEED_INPUT_PROXY_W, power_targets,
                           profile_calibration, profile_settings)
from sabes.session_export import json_value, session_json, session_payload  # noqa: E402


def _profile(profile):
    calibration = profile_calibration(profile)
    settings, detection = profile_settings(profile, calibration)
    result = bridge.run(settings, calibration, detection, solve=False)
    assert result.raw is None
    return result, detection, calibration


def test_paper_profile_preserves_default_settings_calibration_and_bridge_inputs():
    result, detection, calibration = _profile(PAPER_PROFILE)
    original = SetupSettings()
    for field in fields(original):
        expected = getattr(original, field.name)
        actual = getattr(result.settings, field.name)
        if isinstance(expected, float) and np.isnan(expected):
            assert np.isnan(actual), field.name
        else:
            assert actual == expected, field.name
    assert detection == DetectionSettings()
    assert calibration.to_dict() == default_calibration().to_dict()
    assert result.params == bridge.run(solve=False).params
    assert power_targets(PAPER_PROFILE) == (0.600, 7e-3, 8e-6)


def test_lab_profile_drives_the_actual_source_and_detection_inputs():
    result, detection, calibration = _profile(LAB_PROFILE)
    assert result.settings.cell_temp_c == 113.0
    assert result.settings.opd_ghz == 1.0
    assert result.settings.eom_frequency_hz == 3.039e9
    assert result.chain.pump_power_w == pytest.approx(LAB_PUMP_W, rel=1e-7)
    assert result.chain.wanted_seed_sideband_power_w == pytest.approx(
        LAB_SEED_INPUT_PROXY_W, rel=1e-7)
    assert result.params["pump_mw"] == pytest.approx(294.8)
    assert result.params["probe_uw"] == pytest.approx(3.922)
    assert result.params["pump_waist_um"] == 530.0
    assert result.params["probe_waist_um"] == 330.0
    assert result.params["after_cell_eff_pct"] == pytest.approx(
        LAB_AFTER_CELL_TRANSMISSION * 100.0)
    assert result.params["qe_pct"] == pytest.approx(90.5)
    assert result.eta == pytest.approx(LAB_AFTER_CELL_TRANSMISSION * 0.905)
    assert detection.post_cell_transmission_override == LAB_AFTER_CELL_TRANSMISSION
    eater = calibration.get("t_noise_eater")
    assert eater.value == 0.928
    assert eater.provenance == "lab"
    assert "2026-01-07" in eater.note
    assert "Not a new 10/16 measurement" in eater.note


def test_profile_calibrations_are_independent_and_do_not_keep_session_edits():
    lab = profile_calibration(LAB_PROFILE)
    edited = lab.with_value("t_noise_eater", 0.5)
    assert edited.get("t_noise_eater").provenance == "hand"
    assert profile_calibration(LAB_PROFILE).value("t_noise_eater") == 0.928
    paper = profile_calibration(PAPER_PROFILE)
    assert paper.to_dict() == default_calibration().to_dict()
    assert lab.value("t_noise_eater") != paper.value("t_noise_eater")


@pytest.mark.parametrize("profile", PROFILES)
def test_export_preserves_separate_rf_power_planes_and_unknown_termination(profile):
    result, detection, calibration = _profile(profile)
    result = bridge.run(replace(result.settings, eom_rf_dbm=14.5),
                        calibration, detection, solve=False)
    payload = session_payload(result, detection, calibration,
                              default_apparatus(), profile)
    scope = payload["scope"]
    assert scope["rf_generator_display_dbm"] == 11.8
    assert scope["rf_eom_port_measured_dbm"] is None
    assert scope["rf_model_input_dbm"] == 14.5
    assert payload["settings"]["eom_rf_dbm"] == 14.5
    assert scope["detector_gain_switch_v_per_a"] == 1e5
    assert scope["detector_load_measured_ohm"] is None
    facts = {item["id"]: item for item in payload["apparatus"]["facts"]}
    assert facts["eom_port_power"]["value"] is None
    assert facts["detector_load"]["value"] is None
    assert facts["table_coordinates"]["value"] is None
    assert facts["rf_display_setting"]["status"] == "recorded"
    assert "assumptions" in scope["instrument_models"]
    assert "no suppression" in scope["noise_eater"]


def test_export_records_live_knobs_calibration_and_their_derived_result():
    base, detection, calibration = _profile(LAB_PROFILE)
    settings = replace(base.settings, eom_rf_dbm=14.5,
                       seed_trim_hwp_deg=base.settings.seed_trim_hwp_deg + 1.0)
    detection = replace(detection, pd_defocus_mm=5.0,
                        post_cell_transmission_override=0.91,
                        qe_pct_override=88.0)
    calibration = calibration.with_value(
        "t_m_seed_1", 0.8, note="Live session throughput edit")
    result = bridge.run(settings, calibration, detection, solve=False)
    payload = session_payload(result, detection, calibration,
                              default_apparatus(), LAB_PROFILE)
    assert payload["settings"]["seed_trim_hwp_deg"] == settings.seed_trim_hwp_deg
    assert payload["detection_settings"]["pd_defocus_mm"] == 5.0
    assert payload["calibration"]["t_m_seed_1"] == {
        "value": 0.8, "unit": calibration.get("t_m_seed_1").unit,
        "provenance": "hand", "note": "Live session throughput edit"}
    assert payload["derived_parameters"] == json_value(result.params)
    assert payload["derived_parameters"]["probe_uw"] != base.params["probe_uw"]
    assert payload["derived_parameters"]["after_cell_eff_pct"] == 91.0
    assert payload["derived_parameters"]["qe_pct"] == 88.0
    at_cell = next(row for row in payload["source_power_budget"]
                   if row["path"] == "seed" and row["stage"] == "at cell")
    assert at_cell["power_w"] == result.chain.seed.total_power_w
    # The payload is independent; editing a downloaded record cannot alter inputs.
    payload["calibration"]["t_m_seed_1"]["value"] = 0.2
    payload["apparatus"]["facts"][0]["value"] = "changed in export"
    assert calibration.value("t_m_seed_1") == 0.8
    assert default_apparatus().get_fact("table_coordinates").value is None


def test_session_json_replaces_unavailable_numbers_with_null_and_is_strict_json():
    result, detection, calibration = _profile(LAB_PROFILE)
    params = dict(result.params, unavailable=np.array([np.nan, np.inf, -np.inf]),
                  nested={"value": np.float32(np.nan), "count": np.int64(2)})
    result = replace(result, params=params)
    serialized = session_json(result, detection, calibration,
                              default_apparatus(), LAB_PROFILE)
    payload = json.loads(serialized, parse_constant=lambda value:
                         pytest.fail(f"non-JSON numeric constant: {value}"))
    json.dumps(payload, allow_nan=False)
    assert payload["schema_version"] == "sabes-session-v1"
    assert payload["calculation_profile"] == LAB_PROFILE
    assert payload["apparatus"]["baseline_date"] == "2026-10-16"
    assert payload["settings"]["ecdl_temp_c"] is None
    assert payload["derived_parameters"]["unavailable"] == [None, None, None]
    assert payload["derived_parameters"]["nested"] == {"value": None, "count": 2}
    assert payload["operating_point"] == {
        "gains": None, "NF_dB": None, "claim_gate": None}
    assert np.isnan(result.params["unavailable"][0])


@pytest.mark.parametrize("function", [profile_calibration, power_targets,
                                      profile_settings])
def test_unknown_profile_fails_instead_of_silently_using_a_reference(function):
    with pytest.raises(ValueError, match="unknown SABES profile"):
        function("unregistered apparatus profile")
