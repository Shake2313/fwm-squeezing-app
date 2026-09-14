"""Refine the independent QRT only; reuse hash-verified primary audit results.

The original failed report is retained. Reuse is valid because every source and
test in its manifest must still match byte for byte; no physical solve changes.
"""

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import copy
import hashlib
import json
from pathlib import Path

import numpy as np

from .smooth_transport_audit import (
    ROOT, KEYS, QKEYS, REFINEMENTS, TOLERANCES, encode, hashes, relative, worker,
)


def validated_parent(path):
    path = Path(path)
    report = json.loads(path.read_text(encoding="utf-8"))
    if report.get("schema") != "gabes-continuous-rb-transport-v1":
        raise ValueError("expected the original continuous-path audit")
    if not report.get("source_stable_during_run") or report.get("source_sha256") != hashes():
        raise ValueError("parent source manifest no longer matches: rerun the primary audit")
    if report.get("tolerances") != TOLERANCES:
        raise ValueError("parent acceptance criteria differ")
    for j, (rtol, atol) in enumerate(REFINEMENTS):
        case = report["cases"][f"primary_{j}"]
        if (case["rtol"], case["atol"]) != (rtol, atol) or not case["source_stable_during_run"]:
            raise ValueError("parent primary refinement provenance differs")
    return report


def decoded(values, keys):
    return {key: np.asarray(values[key]["real"]) + 1j*np.asarray(values[key]["imag"])
            for key in keys}


def manifest():
    result = hashes()
    name = Path(__file__).resolve().relative_to(ROOT).as_posix()
    result[name] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    return result


def write_audit(parent_path, output, plot):
    parent_path, output, plot = map(Path, (parent_path, output, plot))
    if output.resolve() == plot.resolve() or output.exists() or plot.exists():
        raise FileExistsError("new distinct immutable report and plot paths required")
    parent = validated_parent(parent_path)
    parent_hash = hashlib.sha256(parent_path.read_bytes()).hexdigest()
    before = manifest()
    jobs = [("QRT_coarse", "reference", 2e-12, 2e-16),
            ("QRT_fine", "reference", 2e-13, 2e-17)]
    results = {}
    with ProcessPoolExecutor(max_workers=2) as pool:
        for future in as_completed([pool.submit(worker, job) for job in jobs]):
            row = future.result()
            results[row["label"]] = row
            print(f"{row['label']} finished: {row['elapsed_seconds']:.2f}s", flush=True)
    # These primary values are the same previously solved mathematical problem,
    # with all source bytes and tolerances verified above. Repeating that solve
    # adds no independent evidence; refine only the under-resolved QRT reference.
    primary = [parent["cases"][f"primary_{j}"] for j in range(len(REFINEMENTS))]
    pvalues = [decoded(case["values"], KEYS) for case in primary]
    refinements = [{key: relative(new[key], old[key]) for key in KEYS}
                   for old, new in zip(pvalues[:-1], pvalues[1:])]
    reference = results["QRT_fine"]["values"]
    independent = {key: relative(pvalues[-1][key], reference[key]) for key in QKEYS}
    qrefinement = {key: relative(reference[key], results["QRT_coarse"]["values"][key])
                   for key in QKEYS}
    stable = (before == manifest()
              and parent_hash == hashlib.sha256(parent_path.read_bytes()).hexdigest()
              and all(row["source_stable_during_run"] for row in results.values()))
    quantum = all(row["values"]["audit"]["passed"] for row in primary)
    converged = (len(refinements) == TOLERANCES["required_successive_primary_comparisons"]
                 and all(max(row.values()) < TOLERANCES["successive_primary_relative"] for row in refinements)
                 and max(independent.values()) < TOLERANCES["independent_QRT_relative"]
                 and max(qrefinement.values()) < TOLERANCES["independent_QRT_refinement_relative"])
    report = copy.deepcopy(parent)
    report.update(schema="gabes-continuous-rb-transport-v2",
                  source_sha256=before, source_stable_during_run=stable,
                  primary_refinements=refinements, independent_QRT_errors=independent,
                  independent_QRT_refinement=qrefinement, quantum_controls_passed=quantum,
                  continuous_Gaussian_selected_path_converged=converged,
                  all_declared_controls_passed=bool(converged and quantum and stable),
                  parent_report={"path": str(parent_path), "sha256": parent_hash,
                                 "reused_cases": [f"primary_{j}" for j in range(len(primary))],
                                 "all_parent_source_hashes_match": True,
                                 "reason": "tighten independent QRT tolerances; retain failed v1"},
                  previous_QRT_cases={key: parent["cases"][key] for key in results})
    report["cases"].update(encode(results))
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5), layout="constrained")
    for key in KEYS:
        rows = report["frozen_envelope_error_vs_continuous"]
        axes[0].loglog([row["segments"] for row in rows],
                       [row["relative_error_vs_continuous"][key] for row in rows], "o-", label=key)
    x = np.arange(len(QKEYS))
    axes[1].semilogy(x, [parent["independent_QRT_refinement"][key] for key in QKEYS], "x--", label="original QRT refinement")
    axes[1].semilogy(x, [qrefinement[key] for key in QKEYS], "o-", label="tighter QRT refinement")
    axes[1].semilogy(x, [independent[key] for key in QKEYS], "s-", label="primary vs tighter QRT")
    axes[1].axhline(TOLERANCES["independent_QRT_refinement_relative"], color="black", ls=":", label="QRT refinement limit")
    axes[1].axhline(TOLERANCES["independent_QRT_relative"], color="gray", ls="--", label="agreement limit")
    axes[0].set(xlabel="Frozen Gaussian subdivisions", ylabel="Relative difference", title="Same physics: envelope approximation error")
    axes[1].set(xticks=x, xticklabels=QKEYS, ylabel="Relative difference", title="Independent raw QRT: unchanged criteria")
    for ax in axes:
        ax.legend(fontsize=8)
    fig.suptitle("Selected moving Rb path: " + ("CONVERGED" if report["all_declared_controls_passed"] else "UNCONVERGED") + "; no optical squeezing claim")
    output.parent.mkdir(parents=True, exist_ok=True)
    plot.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(plot, dpi=160)
    plt.close(fig)
    with output.open("x", encoding="utf-8") as handle:
        json.dump(report, handle, ensure_ascii=False, indent=2, allow_nan=False)
        handle.write("\n")
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--parent", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--plot", required=True)
    args = parser.parse_args()
    report = write_audit(args.parent, args.output, args.plot)
    print(json.dumps({key: report[key] for key in ("independent_QRT_errors", "independent_QRT_refinement", "all_declared_controls_passed")}, indent=2))
    if not report["all_declared_controls_passed"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
