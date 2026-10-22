"""Portable SABES session record with measurement planes and unknowns intact."""
from dataclasses import asdict, is_dataclass
from datetime import datetime, timezone
import json
import math


def json_value(value):
    """Strict JSON: unavailable numeric values become null, not NaN/Infinity."""
    if is_dataclass(value):
        return json_value(asdict(value))
    if isinstance(value, dict):
        return {str(k): json_value(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_value(v) for v in value]
    if hasattr(value, "tolist"):
        return json_value(value.tolist())
    if isinstance(value, float):
        return value if math.isfinite(value) else None
    if value is None or isinstance(value, (str, int, bool)):
        return value
    raise TypeError(f"unsupported SABES export value: {type(value).__name__}")


def session_payload(result, detection, calibration, apparatus, profile):
    """Export actual active inputs, not merely the shipped apparatus JSON."""
    raw = result.raw or {}
    rows = [{"path": path, "stage": stage, "power_w": power,
             "waist_radius_m": waist}
            for path, stage, power, waist in result.chain.budget_rows()]
    return json_value({
        "schema_version": "sabes-session-v1",
        "exported_at_utc": datetime.now(timezone.utc).isoformat(),
        "calculation_profile": profile,
        "apparatus": apparatus.as_dict(),
        "settings": result.settings,
        "detection_settings": detection,
        "calibration": {c.name: c.as_dict() for c in calibration},
        "derived_parameters": result.params,
        "detection_geometry": result.geometry,
        "source_power_budget": rows,
        "operating_point": {
            "gains": result.gains,
            "NF_dB": result.gain_referred_noise_db,
            "claim_gate": raw.get("claim_gate"),
        },
        "scope": {
            "display_coordinates": "schematic units; physical x/y/z unknown",
            "distance_plane": "cell output end to D-mirror",
            "rf_generator_display_dbm": apparatus.get_fact("rf_display_setting").value,
            "rf_eom_port_measured_dbm": apparatus.get_fact("eom_port_power").value,
            "rf_model_input_dbm": result.settings.eom_rf_dbm,
            "detector_gain_switch_v_per_a": apparatus.get_fact("detector_gain_switch").value,
            "detector_load_measured_ohm": apparatus.get_fact("detector_load").value,
            "instrument_models": "Virtual traces; 50-ohm SA conversion and "
                                 "stock noise floor are assumptions. Actual "
                                 "Monitor +/- volts are not calibrated here.",
            "noise_eater": "Transmission only; no suppression model.",
            "sql_reference": "Ideal coherent bypass power before unmeasured "
                             "SQL attenuation; not a measured matched SQL.",
        },
    })


def session_json(result, detection, calibration, apparatus, profile):
    return json.dumps(session_payload(result, detection, calibration,
                                     apparatus, profile),
                      ensure_ascii=False, indent=2, allow_nan=False)
