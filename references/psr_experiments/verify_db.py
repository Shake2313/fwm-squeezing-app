"""Focused PSR database checks; no numerical engine or campaign execution."""
from __future__ import annotations

import copy
import csv
import io
import json
import subprocess
import sys

from build_db import (AS_OF, ROOT, UNIT_ALIASES, build_outputs, evidence_verification, export_tables,
                      json_bytes, quantities, read_records, stale_outputs, validate)


def fixture():
    """Small validation fixture, never emitted into the literature database."""
    return [{
        "schema_version": 1, "paper_id": "validation_fixture", "title": "Validation fixture",
        "authors": ["Fixture"], "year": 2000, "publication": {"status": "fixture"},
        "study_type": "theory", "topics": ["validator_fixture"],
        "scope": {"species": "Rb", "line": "D1", "configuration": "PSR", "measurement": "quadrature",
                  "inclusion": "context_only", "reason": "Validation fixture only"},
        "sources": [{"source_id": "main", "url": "https://example.invalid/fixture", "access": "full_text",
                     "version": "fixture", "path": None, "access_method": "fixture", "access_note": "Fixture only"}],
        "experiments": [{"experiment_id": "case", "label": "Fixture case", "kind": "theory_prediction",
            "conditions": [{"quantity": "pump_power", "value": 1, "unit": "mW", "value_type": "reported",
                            "source_id": "main", "locator": "Fixture p. 1"}], "results": []}],
        "findings": [{"finding_id": "result", "claim": "Fixture claim", "source_id": "main",
                      "locator": "Fixture p. 1", "category": "fixture", "model_implication": "Fixture only"}],
    }]


def cell(value):
    if value is None:
        return ""
    if isinstance(value, (dict, list)):
        return json.dumps(value, ensure_ascii=False, allow_nan=False)
    return str(value)


def main():
    papers = read_records()
    checks = []

    def check(name, okay, detail=None):
        checks.append({"id": name, "status": "PASS" if okay else "FAIL", "detail": detail})

    receipt = validate(papers)
    check("schema_evidence_source_hashes_pdf_pages", receipt["status"] == "PASS", receipt)
    if receipt["status"] != "PASS":
        result = {"date": AS_OF, "origin": "⟦AI⟧ Focused PSR DB verification", "status": "FAIL", "checks": checks}
        (ROOT / "validation_detail.json").write_bytes(json_bytes(result))
        print(json.dumps({"status": "FAIL", "failed": [c["id"] for c in checks]}))
        return 1

    payload = json.loads((ROOT / "database.json").read_text(encoding="utf-8"))
    check("combined_json_lossless", payload["papers"] == papers)
    outputs = build_outputs(papers, receipt)
    stale = stale_outputs(outputs)
    check("all_generated_views_and_manifest_fresh", not stale, stale)
    tables = {}
    expected_counts = {"papers.csv": len(papers), "experiments.csv": receipt["experiment_count"],
                       "observations.csv": receipt["quantity_count"], "findings.csv": receipt["finding_count"]}
    for name, expected in expected_counts.items():
        content = (ROOT / name).read_bytes()
        tables[name] = list(csv.DictReader(io.StringIO(content.decode("utf-8-sig"), newline="")))
        check(name + "_encoding_and_row_count", content.startswith(b"\xef\xbb\xbf") and len(tables[name]) == expected,
              {"expected": expected, "actual": len(tables[name])})

    errors = []
    rows = {row["observation_id"]: row for row in tables["observations.csv"]}
    for p in papers:
        sources = {s["source_id"]: s for s in p["sources"]}
        for exp in p["experiments"]:
            for role, index, q in quantities(exp):
                key = f"{p['paper_id']}/{exp['experiment_id']}/{role}/{index}"
                row = rows.get(key, {})
                for field, value in q.items():
                    if field == "verification" and q["value_type"] == "not_reported":
                        continue
                    if row.get(field) != cell(value):
                        errors.append(f"{key}: {field}")
                source = sources.get(q.get("source_id"), {})
                verification = evidence_verification(q, source) if q["value_type"] != "not_reported" else "not_reported"
                if row.get("normalized_unit") != UNIT_ALIASES.get(q["unit"], q["unit"]) or row.get("verification") != verification:
                    errors.append(f"{key}: normalized unit or evidence status")
                if row.get("source_sha256") != cell(source.get("sha256")):
                    errors.append(f"{key}: source checksum")
    check("observations_all_fields_roundtrip", not errors, errors)

    errors = []
    rows = {(r["paper_id"], r["finding_id"]): r for r in tables["findings.csv"]}
    for p in papers:
        sources = {s["source_id"]: s for s in p["sources"]}
        for finding in p["findings"]:
            row = rows.get((p["paper_id"], finding["finding_id"]), {})
            if any(row.get(key) != cell(value) for key, value in finding.items()):
                errors.append(f"{p['paper_id']}/{finding['finding_id']}: fields")
            if row.get("verification") != evidence_verification(finding, sources[finding["source_id"]]):
                errors.append(f"{p['paper_id']}/{finding['finding_id']}: evidence status")
    check("findings_all_fields_and_verification_roundtrip", not errors, errors)

    # A changed expected view must be detected without touching the actual file.
    changed = dict(outputs)
    changed["index.html"] += b"\nSTALE PROBE"
    check("stale_html_detected", "index.html" in stale_outputs(changed))
    before = {name: (ROOT / name).stat().st_mtime_ns for name in outputs}
    command = subprocess.run([sys.executable, str(ROOT / "build_db.py"), "--check"],
                             cwd=ROOT, capture_output=True, text=True, encoding="utf-8")
    after = {name: (ROOT / name).stat().st_mtime_ns for name in outputs}
    check("check_command_read_only", command.returncode == 0 and before == after,
          {"returncode": command.returncode, "output_files_untouched": before == after})

    nominal = fixture()
    check("minimal_fixture_valid", validate(nominal, check_files=False)["status"] == "PASS")

    configurations = copy.deepcopy(nominal)
    parent_scope = {"species": "87Rb", "line": "D1", "configuration": "PSR",
                    "measurement": "vacuum_quadrature", "seed_mode": "vacuum_input"}
    configurations[0]["scope"].update(parent_scope)
    configurations[0]["apparatus_group"] = "source_lab"
    override = copy.deepcopy(configurations[0]["experiments"][0])
    override.update(experiment_id="seeded_case", species="85Rb", line="D2",
                    configuration="seeded_sensor", measurement="polarization_noise",
                    seed_mode="squeezed_input", apparatus_group="sensor_lab")
    configurations[0]["experiments"].append(override)
    exported = export_tables(configurations)["experiments.csv"].decode("utf-8-sig")
    rows = {r["experiment_id"]: r for r in csv.DictReader(io.StringIO(exported, newline=""))}
    fields = [*parent_scope, "apparatus_group"]
    expected_parent = {**parent_scope, "apparatus_group": "source_lab"}
    check("experiment_metadata_inheritance_and_overrides",
          all(rows["case"][key] == expected_parent[key] and rows["seeded_case"][key] == override[key] for key in fields))

    def q(records):
        return records[0]["experiments"][0]["conditions"][0]

    mutations = [
        ("duplicate_papers_rejected", lambda r: r.append(copy.deepcopy(r[0]))),
        ("duplicate_experiments_rejected", lambda r: r[0]["experiments"].append(copy.deepcopy(r[0]["experiments"][0]))),
        ("duplicate_findings_rejected", lambda r: r[0]["findings"].append(copy.deepcopy(r[0]["findings"][0]))),
        ("duplicate_sources_rejected", lambda r: r[0]["sources"].append(copy.deepcopy(r[0]["sources"][0]))),
        ("unknown_quantity_source_rejected", lambda r: q(r).update(source_id="missing")),
        ("unknown_finding_source_rejected", lambda r: r[0]["findings"][0].update(source_id="missing")),
        ("absent_value_rejected", lambda r: q(r).update(value=None)),
        ("missing_is_not_zero", lambda r: q(r).update(value_type="not_reported", value=0)),
        ("derived_formula_required", lambda r: q(r).update(value_type="derived", derivation={})),
        ("reversed_interval_rejected", lambda r: q(r).update(value_type="range", value=None, lower=2, upper=1)),
        ("nonfinite_value_rejected", lambda r: q(r).update(value=float("nan"))),
        ("empty_locator_rejected", lambda r: q(r).update(locator="")),
        ("negative_squeezing_reduction_rejected", lambda r: q(r).update(quantity="squeezing_measured", value=-3, unit="dB")),
        ("nonnumeric_squeezing_reduction_rejected", lambda r: q(r).update(quantity="squeezing_measured", value="three", unit="dB")),
        ("repository_escape_rejected", lambda r: r[0]["sources"][0].update(path="../../../../outside.pdf", sha256="0" * 64)),
        ("abstract_not_fulltext_rejected", lambda r: r[0]["sources"][0].update(access="abstract_only")),
        ("online_fulltext_audit_required", lambda r: r[0]["sources"][0].pop("access_note")),
    ]
    for name, mutate in mutations:
        bad = copy.deepcopy(nominal)
        mutate(bad)
        rejected = validate(bad, check_files=False)
        check(name, rejected["status"] == "FAIL", rejected["errors"])

    metadata = copy.deepcopy(nominal)
    metadata[0]["sources"][0]["access"] = "abstract_only"
    q(metadata)["verification"] = "metadata_only"
    metadata[0]["findings"][0]["verification"] = "metadata_only"
    check("explicit_metadata_only_supported", validate(metadata, check_files=False)["status"] == "PASS")
    metadata_exports = build_outputs(metadata, validate(metadata, check_files=False))
    check("metadata_only_visible_in_experiment_summary", b"[metadata_only]" in metadata_exports["experiments.csv"])
    theory = copy.deepcopy(nominal)
    theory[0]["experiments"] = []
    check("theory_without_experiments_supported", validate(theory, check_files=False)["status"] == "PASS")

    local = next(((pi, si) for pi, p in enumerate(papers) for si, s in enumerate(p["sources"]) if s.get("path")), None)
    if local:
        bad = copy.deepcopy(papers)
        bad[local[0]]["sources"][local[1]]["sha256"] = "0" * 64
        invalid = validate(bad)
        check("source_checksum_corruption_rejected", invalid["status"] == "FAIL" and any("checksum mismatch" in e for e in invalid["errors"]))
    else:
        check("source_checksum_corruption_rejected", False, "No local source available to exercise checksum validation")

    result = {"date": AS_OF, "origin": "⟦AI⟧ Focused PSR DB verification",
              "status": "PASS" if all(c["status"] == "PASS" for c in checks) else "FAIL", "checks": checks,
              "scope": "Data structure, evidence provenance and exports only; no solver, fit or scientific acceptance test."}
    (ROOT / "validation_detail.json").write_bytes(json_bytes(result))
    print(json.dumps({"status": result["status"], "checks": len(checks), "failed": [c["id"] for c in checks if c["status"] != "PASS"]}, ensure_ascii=False))
    return 0 if result["status"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
