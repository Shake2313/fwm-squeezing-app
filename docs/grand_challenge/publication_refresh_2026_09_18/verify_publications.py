"""Publication preservation checks; complements, never replaces, scientific review.

Run from any directory. Does not rebuild documents or alter historical evidence.
Writes the requested new JSON report exclusively.
"""
import argparse
import hashlib
from html.parser import HTMLParser
import json
from pathlib import Path
import re
from urllib.parse import unquote, urlparse

ROOT = Path(__file__).resolve().parents[3]
HERE = Path(__file__).resolve().parent


def read(path):
    return path.read_text(encoding="utf-8")


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def atlas(path):
    html = read(path)
    match = re.search(r'<script type="application/json" id="atlas-data">(.*?)</script>', html, re.S)
    if not match:
        raise ValueError(f"Missing embedded graph in {path}")
    return html, json.loads(match.group(1))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    manifest = json.loads(read(ROOT / "docs/grand_challenge/current_publications.json"))
    checks = []

    def check(name, passed, **evidence):
        checks.append({"name": name, "passed": bool(passed), **evidence})

    previous = json.loads(read(HERE / "historical_artifacts_before.json"))["sha256"]
    changed = [p for p, digest in previous.items() if not (ROOT / p).exists() or sha(ROOT / p) != digest]
    check("historical_documents_unchanged", not changed, file_count=len(previous), changed=changed)

    artifacts = {}
    for name, item in manifest["publications"].items():
        for role in ("source", "artifact"):
            path = ROOT / item[role]
            check(f"{name}_{role}_exists", path.is_file(), path=item[role])
            if path.is_file():
                artifacts[item[role]] = {"bytes": path.stat().st_size, "sha256": sha(path)}
    if not all(row["passed"] for row in checks):
        raise ValueError("Publication files missing or historical files changed: " + str(checks))

    theory = ROOT / manifest["publications"]["analytic_reconstruction"]["source"]
    old_theory = theory.parent / "squeezing_analytic_reconstruction_v2.tex"
    pattern = r"\\label\{([^}]+)\}"
    old_labels = set(re.findall(pattern, read(old_theory)))
    current_labels = re.findall(pattern, read(theory))
    missing = sorted(old_labels - set(current_labels))
    duplicate = sorted(x for x in set(current_labels) if current_labels.count(x) > 1)
    check("analytic_labels_preserved", not missing, previous=len(old_labels), current=len(current_labels), missing=missing)
    check("analytic_labels_unique", not duplicate, duplicates=duplicate)

    report = ROOT / manifest["publications"]["squeezing_report"]["source"]
    figure_pattern = r"\\includegraphics(?:\[[^]]*\])?\{([^}]+)\}"
    old_figs = set(re.findall(figure_pattern, read(report.parent / "squeezing_report_v7.tex")))
    new_figs = set(re.findall(figure_pattern, read(report)))
    # Basenames permit relocation into a dedicated new asset directory.
    missing_figs = sorted({Path(p).name for p in old_figs} - {Path(p).name for p in new_figs})
    check("squeezing_scientific_figures_preserved", not missing_figs, previous=sorted(old_figs), missing=missing_figs)

    html_path = ROOT / manifest["publications"]["quotient_structure"]["artifact"]
    html, graph = atlas(html_path)
    _, old_graph = atlas(html_path.parent / "fwm_quotient_structure_v3.html")
    old_nodes = {n["id"] for g in old_graph["groups"] for n in g["nodes"]}
    new_nodes = [n["id"] for g in graph["groups"] for n in g["nodes"]]
    old_pairs = {(e["source"], e["target"]) for e in old_graph["edges"]}
    new_pairs = {(e["source"], e["target"]) for e in graph["edges"]}
    check("quotient_node_identity_preserved", old_nodes <= set(new_nodes), previous=len(old_nodes), current=len(new_nodes), missing=sorted(old_nodes-set(new_nodes)))
    check("quotient_relations_preserved", old_pairs <= new_pairs, previous=len(old_pairs), current=len(new_pairs), missing=sorted(old_pairs-new_pairs))
    check("quotient_node_ids_unique", len(new_nodes) == len(set(new_nodes)))
    bad_ends = [e for e in graph["edges"] if e["source"] not in new_nodes or e["target"] not in new_nodes]
    check("quotient_edge_endpoints_exist", not bad_ends, invalid=bad_ends)
    deps = re.findall(r'<(?:script|link)[^>]+(?:src|href)="(https?://[^"]+)"', html)
    check("quotient_no_network_runtime_dependency", not deps, dependencies=deps)
    class Links(HTMLParser):
        def __init__(self):
            super().__init__()
            self.ids, self.anchors = set(), set()
            self.local_files = set()

        def handle_starttag(self, tag, attrs):
            values = dict(attrs)
            if values.get("id"):
                self.ids.add(values["id"])
            if values.get("href", "").startswith("#"):
                self.anchors.add(values["href"][1:])
            href = values.get("href", "")
            parsed = urlparse(href)
            if tag == "a" and not parsed.scheme and parsed.path:
                self.local_files.add(unquote(parsed.path))

    links = Links()
    links.feed(html)
    ids, anchors = links.ids, links.anchors
    missing_anchors = sorted(anchors-ids)
    check("quotient_internal_links_resolve", not missing_anchors, unresolved=missing_anchors)
    missing_files = sorted(path for path in links.local_files
                           if not (html_path.parent / path).is_file())
    check("quotient_local_source_and_companion_links_resolve", not missing_files,
          local_file_count=len(links.local_files), unresolved=missing_files)

    try:
        import fitz
    except ImportError:
        check("pdf_text_inspection_available", False, reason="PyMuPDF unavailable; visual inspection remains separate")
    else:
        for name in ("squeezing_report", "analytic_reconstruction"):
            path = ROOT / manifest["publications"][name]["artifact"]
            with fitz.open(path) as pdf:
                texts = [page.get_text() for page in pdf]
                blank = [i+1 for i, text in enumerate(texts) if len(text.strip()) < 4]
                unknown_refs = [i+1 for i, text in enumerate(texts) if "??" in text]
                check(name+"_pdf_text", len(pdf) > 0 and not blank and not unknown_refs,
                      pages=len(pdf), near_empty_pages=blank, unresolved_reference_pages=unknown_refs)

    result = {"scope": "Structural preservation, build/link integrity and archival byte preservation; semantic coverage and rendered-page review are recorded separately.",
              "all_checks_passed": all(row["passed"] for row in checks),
              "checks": checks, "artifact_hashes": artifacts}
    with args.output.open("x", encoding="utf-8") as stream:
        json.dump(result, stream, ensure_ascii=False, indent=2)
        stream.write("\n")
    print(json.dumps(result, ensure_ascii=True, indent=2))
    if not result["all_checks_passed"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
