"""Build the PSR literature evidence views without importing the solver.

Canonical input: papers/*.json. Run --check to validate and detect stale views
without writing files. Paths resolve from this file, independent of cwd.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import math
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parents[1]
AS_OF = "2026-10-30"
UNIT_ALIASES = {"µW": "uW", "μW": "uW", "µm": "um", "μm": "um", "°C": "degC", "%": "percent"}
REDUCTIONS = {"squeezing_measured", "squeezing_electronic_noise_subtracted", "squeezing_loss_corrected"}


def json_bytes(value):
    return (json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False) + "\n").encode("utf-8")


def read_records():
    records = []
    for path in sorted((ROOT / "papers").glob("*.json")):
        paper = json.loads(path.read_text(encoding="utf-8-sig"))
        if path.stem != paper.get("paper_id"):
            raise ValueError(f"{path.name}: filename must match paper_id")
        records.append(paper)
    return records


def quantities(exp):
    for role in ("conditions", "results"):
        for index, item in enumerate(exp.get(role, [])):
            yield role, index, item


def source_path(path):
    """Allow repository-local reuse, including ../, but no external paths."""
    if not isinstance(path, str) or not path or "\\" in path or ":" in path:
        raise ValueError("source path must be a relative POSIX path")
    relative = Path(path)
    target = (ROOT / relative).resolve()
    if relative.is_absolute() or not target.is_relative_to(REPO):
        raise ValueError("source path outside repository")
    return target


def finite_numbers(value):
    if isinstance(value, float):
        return math.isfinite(value)
    if isinstance(value, dict):
        return all(finite_numbers(v) for v in value.values())
    if isinstance(value, list):
        return all(finite_numbers(v) for v in value)
    return True


def evidence_verification(item, source):
    return item.get("verification", "full_text_checked" if source.get("access") == "full_text" else "metadata_only")


def validate(papers, *, check_files=True):
    import jsonschema

    schema = json.loads((ROOT / "schema.json").read_text(encoding="utf-8"))
    validator = jsonschema.Draft202012Validator(schema)
    errors, warnings = [], []
    ids, exp_ids, finding_ids = set(), set(), set()
    experiment_count = quantity_count = source_count = finding_count = 0
    if not papers:
        errors.append("no canonical paper records")
    for paper in papers:
        pid = paper.get("paper_id", "MISSING") if isinstance(paper, dict) else "MISSING"
        if pid in ids:
            errors.append(f"duplicate paper_id {pid}")
        ids.add(pid)
        schema_errors = list(validator.iter_errors(paper))
        errors.extend(f"{pid}/{'.'.join(map(str, e.absolute_path))}: {e.message}" for e in schema_errors)
        if schema_errors:
            continue
        if not finite_numbers(paper):
            errors.append(f"{pid}: nonfinite numeric value")
        source_map = {s["source_id"]: s for s in paper["sources"]}
        if len(source_map) != len(paper["sources"]):
            errors.append(f"{pid}: duplicate source_id")
        source_count += len(paper["sources"])
        for source in paper["sources"]:
            label = f"{pid}/{source['source_id']}"
            path = source.get("path")
            if source["access"] != "full_text":
                warnings.append(f"{label}: {source['access']}; evidence must be marked metadata_only")
            if not path:
                if source["access"] == "full_text":
                    warnings.append(f"{label}: full text inspected online; no local snapshot")
                    if not source.get("access_method") or not source.get("access_note"):
                        errors.append(f"{label}: online full text needs access_method and access_note")
                continue
            try:
                target = source_path(path)
            except ValueError as exc:
                errors.append(f"{label}: {exc}")
                continue
            if not source.get("sha256"):
                errors.append(f"{label}: local source lacks checksum")
            if not check_files:
                continue
            if not target.is_file():
                errors.append(f"{label}: source missing {path}")
                continue
            content = target.read_bytes()
            if hashlib.sha256(content).hexdigest() != source.get("sha256"):
                errors.append(f"{label}: checksum mismatch {path}")
            if target.suffix.lower() == ".pdf":
                if not content[:1024].lstrip().startswith(b"%PDF"):
                    errors.append(f"{label}: invalid PDF header")
                    continue
                try:
                    import pymupdf
                    with pymupdf.open(target) as document:
                        if source.get("pages") != len(document):
                            errors.append(f"{label}: PDF page count mismatch ({source.get('pages')} != {len(document)})")
                except Exception as exc:
                    errors.append(f"{label}: PDF unreadable: {exc}")

        def check_evidence(item, label):
            source = source_map.get(item.get("source_id"))
            if source is None:
                errors.append(f"{label}: unknown source")
                return
            if not item.get("locator", "").strip():
                errors.append(f"{label}: missing locator")
            if source["access"] != "full_text" and item.get("verification") != "metadata_only":
                errors.append(f"{label}: non-full-text evidence requires explicit metadata_only verification")

        for exp in paper["experiments"]:
            experiment_count += 1
            eid = f"{pid}/{exp['experiment_id']}"
            if eid in exp_ids:
                errors.append(f"duplicate experiment_id {eid}")
            exp_ids.add(eid)
            for role, index, q in quantities(exp):
                quantity_count += 1
                label = f"{eid}/{role}/{index} {q['quantity']}"
                if q["value_type"] == "not_reported":
                    if q.get("source_id"):
                        check_evidence(q, label)
                    continue
                check_evidence(q, label)
                if q["value_type"] == "range" and q["lower"] > q["upper"]:
                    errors.append(f"{label}: reversed interval")
                if q["quantity"] in REDUCTIONS:
                    values = (q.get("lower"), q.get("upper")) if q["value_type"] == "range" else (q["value"],)
                    if q["unit"] != "dB" or any(not isinstance(v, (int, float)) or isinstance(v, bool) or v < 0 for v in values):
                        errors.append(f"{label}: squeezing reduction must be nonnegative dB below SQL")
        for finding in paper["findings"]:
            finding_count += 1
            fid = f"{pid}/{finding['finding_id']}"
            if fid in finding_ids:
                errors.append(f"duplicate finding_id {fid}")
            finding_ids.add(fid)
            check_evidence(finding, fid)
    valid_objects = [p for p in papers if isinstance(p, dict)]
    return {
        "schema_version": 1, "validation_date": AS_OF,
        "status": "PASS" if not errors else "FAIL", "paper_count": len(papers),
        "experiment_count": experiment_count, "quantity_count": quantity_count,
        "finding_count": finding_count, "source_count": source_count,
        "study_type_counts": dict(Counter(p.get("study_type", "missing") for p in valid_objects)),
        "full_text_paper_count": sum(any(s.get("access") == "full_text" for s in p.get("sources", []) if isinstance(s, dict)) for p in valid_objects),
        "errors": errors, "warnings": warnings,
        "scope": "Structure, evidence links, source bytes and page counts; scientific interpretation requires source review.",
    }


def display_value(q):
    if q["value_type"] == "not_reported":
        return "미보고"
    value = f"{q['lower']} … {q['upper']}" if q["value_type"] == "range" else str(q["value"])
    prefix = {"approximate": "≈", "lower_bound": "≥", "upper_bound": "≤"}.get(q["value_type"], "")
    return f"{prefix}{value} {q['unit']}".strip()


def csv_bytes(rows, preferred):
    fields = preferred + sorted(set().union(*(r.keys() for r in rows)) - set(preferred)) if rows else preferred
    stream = io.StringIO(newline="")
    writer = csv.DictWriter(stream, fieldnames=fields, lineterminator="\n")
    writer.writeheader()
    for row in rows:
        writer.writerow({k: json.dumps(v, ensure_ascii=False, allow_nan=False) if isinstance(v, (dict, list)) else v for k, v in row.items()})
    return stream.getvalue().encode("utf-8-sig")


def source_columns(source):
    return {"source_url": source.get("url"), "source_path": source.get("path"),
            "source_sha256": source.get("sha256"), "source_access": source.get("access"),
            "source_version": source.get("version")}


def export_tables(papers):
    paper_rows, experiment_rows, observation_rows, finding_rows = [], [], [], []
    for p in papers:
        pid = p["paper_id"]
        paper_rows.append(dict(paper_id=pid, title=p["title"], authors=p["authors"], year=p["year"],
            study_type=p["study_type"], topics=p["topics"], **p["publication"], **p["scope"],
            apparatus_group=p.get("apparatus_group"), experiment_count=len(p["experiments"]),
            finding_count=len(p["findings"]), source_urls=[s["url"] for s in p["sources"]],
            full_text_available=any(s["access"] == "full_text" for s in p["sources"])))
        sm = {s["source_id"]: s for s in p["sources"]}
        for exp in p["experiments"]:
            row = dict(paper_id=pid, experiment_id=exp["experiment_id"], title=p["title"], study_type=p["study_type"],
                       label=exp["label"], kind=exp["kind"], inclusion=p["scope"]["inclusion"], notes=exp.get("notes", []))
            for key in ("species", "line", "configuration", "measurement", "seed_mode"):
                row[key] = exp.get(key, p["scope"].get(key))
            row["apparatus_group"] = exp.get("apparatus_group", p.get("apparatus_group"))
            verifications = set()
            for role, index, q in quantities(exp):
                name = q["quantity"]
                source = sm.get(q.get("source_id"), {})
                verification = evidence_verification(q, source) if q["value_type"] != "not_reported" else "not_reported"
                verifications.add(verification)
                value = display_value(q) + (" [metadata_only]" if verification == "metadata_only" else "")
                row[name] = " | ".join(filter(None, (row.get(name, ""), value)))
                observation_rows.append({"paper_id": pid, "experiment_id": exp["experiment_id"],
                    "observation_id": f"{pid}/{exp['experiment_id']}/{role}/{index}", "role": role, **q,
                    "normalized_unit": UNIT_ALIASES.get(q["unit"], q["unit"]),
                    "verification": verification,
                    **source_columns(source)})
            row["evidence_verification"] = sorted(verifications)
            experiment_rows.append(row)
        for finding in p["findings"]:
            source = sm[finding["source_id"]]
            finding_rows.append({"paper_id": pid, "title": p["title"], "study_type": p["study_type"], **finding,
                "verification": evidence_verification(finding, source), **source_columns(source)})
    return {
        "papers.csv": csv_bytes(paper_rows, ["paper_id", "title", "authors", "year", "study_type", "topics", "journal", "doi", "arxiv_id", "status", "inclusion", "full_text_available", "experiment_count", "finding_count"]),
        "experiments.csv": csv_bytes(experiment_rows, ["paper_id", "experiment_id", "title", "study_type", "label", "kind", "inclusion", "species", "line", "configuration", "measurement", "seed_mode", "apparatus_group", "evidence_verification"]),
        "observations.csv": csv_bytes(observation_rows, ["paper_id", "experiment_id", "observation_id", "role", "quantity", "value", "lower", "upper", "unit", "normalized_unit", "value_type", "verification", "basis", "reference", "plane", "locator", "evidence", "notes", "source_id", "source_url", "source_path", "source_sha256", "source_access"]),
        "findings.csv": csv_bytes(finding_rows, ["paper_id", "finding_id", "title", "study_type", "category", "claim", "model_implication", "verification", "source_id", "locator", "source_url", "source_path", "source_sha256", "source_access"]),
    }


def build_outputs(papers, receipt):
    outputs = export_tables(papers)
    outputs["database.json"] = json_bytes({"schema_version": 1, "as_of": AS_OF,
        "purpose": "PSR literature evidence and model-development context", "unit_aliases": UNIT_ALIASES, "papers": papers})
    outputs["validation.json"] = json_bytes(receipt)
    template = (ROOT / "viewer.html").read_text(encoding="utf-8")
    data = json.dumps(papers, ensure_ascii=False, allow_nan=False).replace("<", "\\u003c").replace("&", "\\u0026")
    receipt_data = json.dumps(receipt, ensure_ascii=False).replace("<", "\\u003c").replace("&", "\\u0026")
    outputs["index.html"] = template.replace("__DB_DATA__", data).replace("__DB_RECEIPT__", receipt_data).encode("utf-8")
    members = sorted((ROOT / "papers").glob("*.json")) + [ROOT / name for name in (
        "SCHEMA.md", "schema.json", "build_db.py", "verify_db.py", "viewer.html", "README.md", "SEARCH.md", "REVIEW.md")]
    hashes = {p.relative_to(ROOT).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest() for p in members if p.is_file()}
    hashes.update({name: hashlib.sha256(content).hexdigest() for name, content in outputs.items()})
    outputs["manifest.json"] = json_bytes({"as_of": AS_OF, "files": dict(sorted(hashes.items())),
        "sources": {s["path"]: s.get("sha256") for p in papers for s in p["sources"] if s.get("path")}})
    return outputs


def stale_outputs(outputs):
    return [name for name, content in outputs.items() if not (ROOT / name).is_file() or (ROOT / name).read_bytes() != content]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true", help="Validate and compare generated bytes; write nothing")
    args = parser.parse_args()
    try:
        papers = read_records()
        receipt = validate(papers)
        if receipt["status"] == "PASS":
            outputs = build_outputs(papers, receipt)
            if args.check:
                stale = stale_outputs(outputs)
                receipt = {**receipt, "stale_outputs": stale}
                if stale:
                    receipt["status"] = "FAIL"
                    receipt["errors"] = ["stale or missing generated views: " + ", ".join(stale)]
            else:
                for name, content in outputs.items():
                    (ROOT / name).write_bytes(content)
        print(json.dumps(receipt, ensure_ascii=False, indent=2))
        return 0 if receipt["status"] == "PASS" else 1
    except (ValueError, OSError) as exc:
        print(json.dumps({"status": "FAIL", "errors": [str(exc)]}, ensure_ascii=False))
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
