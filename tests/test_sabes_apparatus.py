"""The current apparatus inventory must preserve identity and evidence scope."""
from copy import deepcopy
import json
from pathlib import Path
import sys

import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from sabes.apparatus import ApparatusSnapshot, default_apparatus  # noqa: E402
from sabes.layout import parts  # noqa: E402


def test_default_snapshot_is_dated_and_layout_mappings_resolve():
    snapshot = default_apparatus()
    assert snapshot.baseline_date == "2026-10-16"
    assert snapshot.confirmed_on == "2026-10-21"
    snapshot.validate_layouts(parts.LAYOUTS)
    assert {item.id for item in snapshot.equipment_for_node("eom")} == {
        "fiber_eom", "rf_generator"}
    assert snapshot.equipment_for_node("not_an_optic") == ()
    assert all(snapshot.equipment_for_node(node_id) for node_id in (
        "m_sql_bypass", "hwp_sql", "bs_sql", "flip_sql_probe", "flip_sql_conj"))


def test_current_equipment_overrides_historical_model_candidates():
    snapshot = default_apparatus()
    assert snapshot.get_equipment("ecdl_ta_controller").model == "LDC501"
    assert snapshot.get_equipment("spectrum_analyzer").model == "DSA875"
    assert snapshot.get_equipment("fiber_eom").model == "PM-0S5-10-PFA-PFA-780-UL"
    assert snapshot.get_equipment("etalon_controller").model == "TC200"
    assert snapshot.get_equipment("pbs").model == "CCM1-PBS252/M"
    sql_pbs = snapshot.get_equipment("sql_beamsplitter")
    assert sql_pbs.model == "CCM1-PBS252/M"
    assert sql_pbs.manufacturer == "Thorlabs"
    assert sql_pbs.status == "confirmed"
    assert "sql_confirmed_procedure" in sql_pbs.source_ids
    assert snapshot.get_equipment("cylindrical_lens").status == "tentative"


def test_rf_and_monitor_connections_preserve_actual_route_and_unknown_ports():
    snapshot = default_apparatus()
    connections = {item.id: item for item in snapshot.connections}
    rf = connections["rf_generator_to_eom"]
    assert (rf.source_id, rf.target_id) == ("rf_generator", "fiber_eom")
    assert rf.via == ("N-to-SMA adapter", "SMA cable")
    assert rf.target_port is None
    assert connections["detector_rf_to_sa"].target_id == "spectrum_analyzer"
    for key, port in (("detector_monitor_plus", "Monitor +"),
                      ("detector_monitor_minus", "Monitor −")):
        monitor = connections[key]
        assert monitor.source_port == port
        assert monitor.target_id == "oscilloscope"
        assert monitor.target_port is None
    assert snapshot.get_fact("rf_display_setting").as_of == "2026-09-03"
    assert snapshot.get_fact("eom_port_power").value is None


def test_physical_state_is_not_promoted_to_feedback_or_noise_calibration():
    snapshot = default_apparatus()
    eater = snapshot.get_equipment("noise_eater").attributes
    assert eater["inserted"] and eater["powered_on"]
    assert not eater["optimized"]
    assert eater["after_node_ids"] == ["etalon_3", "gt_seed"]
    assert eater["model_scope"] == "transmission_only"
    assert eater["claimed_noise_suppression_db"] is None
    assert snapshot.get_fact("sas_lock").value is False
    assert snapshot.get_fact("table_grounded").value is False


def test_fwm_heater_readouts_and_hanle_expectation_keep_distinct_scope():
    snapshot = default_apparatus()
    assert snapshot.get_fact("fwm_vt_setpoint").value == 2.262
    assert snapshot.get_fact("fwm_vt_tp1").value == 2.227
    assert snapshot.get_fact("fwm_od_fit_temperature").value == 113.0
    assert snapshot.get_fact("heater_schematic_feedback").status == "schematic"
    similarity = snapshot.get_fact("hanle_vt_similarity")
    assert similarity.scope == "Hanle"
    assert similarity.status == "expectation"
    assert snapshot.get_fact("hanle_vt_calibration").value is None
    assert snapshot.get_fact("cold_spot_temperature").value is None
    assert snapshot.get_fact("hanle_before_cell_power").value == 27.17
    assert snapshot.get_fact("hanle_after_cell_power").value == 26.53


def test_unknown_fiber_specs_and_sql_ratio_remain_null():
    snapshot = default_apparatus()
    fiber = snapshot.get_equipment("pump_fiber").attributes
    assert fiber["connector_type"] == "APC-to-APC"
    assert not fiber["oz_high_power_in_use"]
    assert all(fiber[key] is None for key in ("mfd_um", "length_m", "pm", "ar_coating"))
    assert snapshot.get_equipment("sql_beamsplitter").attributes[
        "measured_splitting_ratio"] is None
    assert snapshot.get_fact("sql_total_detected_power_matched").value is True
    assert snapshot.get_fact("sql_total_detected_power_matched").status == "confirmed"


def test_contextual_deferral_does_not_raise_unrelated_questions():
    snapshot = default_apparatus()
    assert snapshot.get_equipment("waveform_generator").model is None
    assert snapshot.get_equipment("waveform_generator").status == "deferred"
    assert snapshot.deferred_questions("apparatus_inventory") == ()
    question, = snapshot.deferred_questions("coil_control")
    assert question.id == "waveform_generator_model"
    assert question.equipment_id == "waveform_generator"


def test_export_roundtrip_keeps_sources_statuses_nulls_and_independent_data(tmp_path):
    snapshot = default_apparatus()
    path = tmp_path / "apparatus.json"
    path.write_text(snapshot.to_json(indent=2), encoding="utf-8")
    restored = ApparatusSnapshot.load(path)
    assert restored.as_dict() == snapshot.as_dict()
    assert json.loads(restored.to_json())["equipment"][0]["model"] is None
    exported = restored.as_dict()
    exported["equipment"][0]["attributes"]["sas_locked"] = True
    assert restored.get_equipment("ecdl").attributes["sas_locked"] is False
    rows = restored.connection_rows("Monitor")
    assert len(rows) == 2 and all(row["to"] == "Oscilloscope" for row in rows)
    assert all(row["scope"] == "Hanle" for row in restored.fact_rows("hanle"))


@pytest.mark.parametrize("field,value,match", [
    ("schema_version", "future_v99", "unsupported apparatus schema"),
    ("baseline_date", "yesterday", "expected an ISO date"),
    ("confirmed_on", "2026-10-15", "precedes baseline_date"),
])
def test_invalid_snapshot_metadata_is_rejected(field, value, match):
    data = default_apparatus().as_dict()
    data[field] = value
    with pytest.raises(ValueError, match=match):
        ApparatusSnapshot.from_dict(data)


@pytest.mark.parametrize("mutation,match", [
    (lambda data: data["equipment"].append(deepcopy(data["equipment"][0])),
     "duplicate record IDs"),
    (lambda data: data["equipment"][0].update(status="guessed"), "unknown status"),
    (lambda data: data["equipment"][0].update(source_ids=["missing"]),
     "missing or unknown source IDs"),
    (lambda data: data["connections"][0].update(target_id="missing"),
     "unknown equipment endpoint"),
    (lambda data: data["equipment"][0].update(model=""),
     "nonempty string or null"),
    (lambda data: data["facts"].append({
        "id": "unknown_with_number", "category": "geometry", "label": "unknown",
        "value": 0.0, "status": "unknown", "source_ids": ["user_20261021"]}),
     "must retain null"),
])
def test_invalid_record_evidence_and_links_are_rejected(mutation, match):
    data = default_apparatus().as_dict()
    mutation(data)
    with pytest.raises(ValueError, match=match):
        ApparatusSnapshot.from_dict(data)


def test_stale_layout_node_mapping_is_rejected_separately_from_loading():
    data = default_apparatus().as_dict()
    data["equipment"][0]["node_ids"] = ["obsolete_node"]
    snapshot = ApparatusSnapshot.from_dict(data)
    with pytest.raises(ValueError, match="unknown layout nodes.*obsolete_node"):
        snapshot.validate_layouts(parts.LAYOUTS)


def test_missing_lookup_raises_instead_of_supplying_a_default():
    snapshot = default_apparatus()
    with pytest.raises(KeyError, match="no apparatus equipment"):
        snapshot.get_equipment("missing")
    with pytest.raises(KeyError, match="no apparatus fact"):
        snapshot.get_fact("missing")
