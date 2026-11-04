"""Acquire declared public sources and bibliographic metadata for manual review.

Downloads never imply full-text verification. Canonical evidence is curated in
papers/*.json after reading the archived source. Existing PDFs are reused.
"""
from __future__ import annotations

import hashlib
import json
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import pymupdf
import requests
from bs4 import BeautifulSoup

ROOT = Path(__file__).resolve().parent
AS_OF = "2026-11-03"

# Source URLs were discovered through publisher, author, NIST and arXiv pages.
SEEDS = [
    ("mohapatra_2007_prl", "10.1103/PhysRevLett.98.113003", "quant-ph/0612200",
     ["https://harvest.aps.org/v2/journals/articles/10.1103/PhysRevLett.98.113003/fulltext"]),
    ("sedlacek_2012_natphys", "10.1038/nphys2423", "1205.4461", []),
    ("sedlacek_2013_prl", "10.1103/PhysRevLett.111.063001", "1304.4299",
     ["https://harvest.aps.org/v2/journals/articles/10.1103/PhysRevLett.111.063001/fulltext"]),
    ("holloway_2014_ieee", "10.1109/TAP.2014.2360208", "1405.7066",
     ["https://tsapps.nist.gov/publication/get_pdf.cfm?pub_id=915858"]),
    ("simons_2016_apl", "10.1063/1.4947231", None,
     ["https://tsapps.nist.gov/publication/get_pdf.cfm?pub_id=920398"]),
    ("kumar_2017_scirep", "10.1038/srep42981", "1610.09550",
     ["https://www.nature.com/articles/srep42981.pdf"]),
    ("kumar_2017_oe", "10.1364/OE.25.008625", "1702.00494",
     ["https://opg.optica.org/oe/viewmedia.cfm?uri=oe-25-8-8625&seq=0"]),
    ("simons_2019_apl", "10.1063/1.5088821", None,
     ["https://tsapps.nist.gov/publication/get_pdf.cfm?pub_id=927353"]),
    ("gordon_2019_aipadv", "10.1063/1.5095633", "1903.09712",
     ["https://tsapps.nist.gov/publication/get_pdf.cfm?pub_id=927624"]),
    ("jing_2020_natphys", "10.1038/s41567-020-0918-5", None,
     ["https://laserspec.sxu.edu.cn/docs/2020-09/6f9a83e33f82405caa24290a90e1a35f.pdf"]),
    ("prajapati_2021_apl", "10.1063/5.0069195", "2108.13226",
     ["https://tsapps.nist.gov/publication/get_pdf.cfm?pub_id=933094"]),
    ("meyer_2021_pra", "10.1103/PhysRevA.104.043103", "2105.10494",
     ["https://harvest.aps.org/v2/journals/articles/10.1103/PhysRevA.104.043103/fulltext"]),
    ("cai_2022_photonics", "10.3390/photonics9040250", None,
     ["https://mdpi-res.com/d_attachment/photonics/photonics-09-00250/article_deploy/photonics-09-00250.pdf"]),
    ("su_2022_oe", "10.1364/OE.444894", "2111.13408",
     ["https://opg.optica.org/oe/viewmedia.cfm?uri=oe-30-2-1499&seq=0"]),
    ("berweger_2022_pra", "10.1103/PhysRevApplied.19.034078", "2209.02583", []),
    ("wang_2023_oe", None, "2303.06421",
     ["https://laserspec.sxu.edu.cn/docs/2023-06/b3890c5ab31c4230abd9adce20c42c64.pdf"]),
    ("zhang_2024_quantum_scaling", None, "2307.15492", []),
    ("yang_2024_pra", "10.1103/PhysRevA.109.032609", None,
     ["https://harvest.aps.org/v2/journals/articles/10.1103/PhysRevA.109.032609/fulltext"]),
    ("tu_2024_sciadv", "10.1126/sciadv.ads0683", "2307.15617", []),
    ("prajapati_2024_bandwidth", None, "2402.17942",
     ["https://tsapps.nist.gov/publication/get_pdf.cfm?pub_id=957473"]),
    ("ju_2026_preprint", None, "2606.04354v1", []),
    ("continuous_2026_comphys", "10.1038/s42005-026-02529-3", None,
     ["https://www.nature.com/articles/s42005-026-02529-3.pdf"]),
    ("schlossberger_2026_jscheme", None, "2603.08678v1", []),
    ("mems_2026_micronano", "10.1038/s41378-026-01216-1", "2509.01911",
     ["https://www.nature.com/articles/s41378-026-01216-1.pdf"]),
    ("offresonant_2026_sensors", "10.3390/s26185887", "2606.10151",
     ["https://www.mdpi.com/1424-8220/26/18/5887/pdf"]),
    ("multichannel_2025_natcom", "10.1038/s41467-025-64439-2", None,
     ["https://www.nature.com/articles/s41467-025-64439-2.pdf"]),
]


def collect(seed):
    pid, doi, arxiv, urls = seed
    record = {"paper_id": pid, "doi": doi, "arxiv_id": arxiv,
              "accessed_on": AS_OF, "attempts": []}
    session = requests.Session()
    session.headers["User-Agent"] = "GABES-Literature-Reference/1.0 (research archive)"
    if arxiv:
        try:
            response = session.get("https://arxiv.org/abs/" + arxiv, timeout=30)
            response.raise_for_status()
            soup = BeautifulSoup(response.text, "html.parser")
            meta = {}
            for item in soup.select('meta[name^="citation_"]'):
                meta.setdefault(item["name"], []).append(item.get("content", ""))
            record["arxiv_metadata"] = meta
            if not doi:
                doi = next(iter(meta.get("citation_doi", [])), None)
                record["doi"] = doi
        except Exception as exc:
            record["attempts"].append({"url": "https://arxiv.org/abs/" + arxiv, "error": str(exc)})
    if doi:
        try:
            response = session.get("https://api.crossref.org/works/" + doi, timeout=30)
            response.raise_for_status()
            record["crossref"] = response.json()["message"]
        except Exception as exc:
            record["attempts"].append({"url": "https://api.crossref.org/works/" + doi, "error": str(exc)})
    candidates = list(urls)
    if arxiv:
        candidates.append("https://arxiv.org/pdf/" + arxiv)
    output = ROOT / "sources" / (pid + ".pdf")
    previous = ROOT / "sources" / (pid + ".access.json")
    if output.exists() and previous.exists():
        old = json.loads(previous.read_text(encoding="utf-8"))
        for key in ("source_url", "final_url", "path", "sha256", "pages", "review_status", "access_method"):
            if key in old:
                record[key] = old[key]
        record["attempts"] += old.get("attempts", [])
    else:
        for url in candidates:
            try:
                response = session.get(url, timeout=40)
                response.raise_for_status()
                if not response.content[:1024].lstrip().startswith(b"%PDF"):
                    raise ValueError("Response is not PDF")
                document = pymupdf.open(stream=response.content, filetype="pdf")
                output.write_bytes(response.content)
                record.update(source_url=url, final_url=response.url, path="sources/" + output.name,
                              sha256=hashlib.sha256(response.content).hexdigest(), pages=len(document),
                              review_status="downloaded_pending_manual_review")
                text = "\n".join(f"\n=== PDF PAGE {i + 1} ===\n" + page.get_text(sort=True)
                                 for i, page in enumerate(document))
                output.with_suffix(".txt").write_text(text, encoding="utf-8")
                document.close()
                break
            except Exception as exc:
                record["attempts"].append({"url": url, "error": str(exc)})
    previous.write_text(json.dumps(record, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    return record


def main():
    seeds = [s for s in SEEDS if len(sys.argv) == 1 or s[0] in sys.argv[1:]]
    report_path = ROOT / "research_sources.json"
    records = json.loads(report_path.read_text(encoding="utf-8")) if report_path.exists() else {}
    with ThreadPoolExecutor(max_workers=4) as pool:
        futures = [pool.submit(collect, seed) for seed in seeds]
        for future in as_completed(futures):
            record = future.result()
            records[record["paper_id"]] = record
            title = record.get("crossref", {}).get("title", record.get("arxiv_metadata", {}).get("citation_title", []))
            print(record["paper_id"], record.get("pages", "NO_PDF"), str(title)[:160], flush=True)
    report_path.write_text(json.dumps(records, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
