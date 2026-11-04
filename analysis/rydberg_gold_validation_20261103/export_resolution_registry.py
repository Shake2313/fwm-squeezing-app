"""Export the current, scoped variable-resolution evidence as a flat table."""
from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent


def main():
    path = ROOT / "current_variable_resolution_registry.json"
    registry = json.loads(path.read_text(encoding="utf-8"))
    selected = registry["selected_current_policy"]
    errors = {r["detuning_policy"]: r["maximum_relative_error"] for r in selected["reports"]}
    rows = []
    for item in registry["variables"]:
        sparse = item["variable"] in ("temperature_c", "coupling_power_mw")
        rows.append(dict(
            origin="⟦AI⟧",
            variable=item["variable"], role=item["role"],
            recommendation=item["recommendation"],
            training_nodes=json.dumps(item.get("proposed_training_nodes", [])),
            evidence=item.get("adaptive_evidence", item["evidence"]),
            scope=item["scope"],
            numerical_status=selected["status"] if sparse else "NO_COARSENING_CERTIFICATE",
            maximum_relative_error_static=errors["maximum_static_slope"] if sparse else "",
            maximum_relative_error_center=errors["bare_two_photon_center"] if sparse else "",
            physical_validation=registry["physical_model_validation"],
            source_registry_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
            invalidate_when="; ".join(registry["invalidate_when"]),
        ))
    with (ROOT / "current_variable_resolution_registry.csv").open(
            "w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(f"Exported {len(rows)} current variable roles")


if __name__ == "__main__":
    main()
