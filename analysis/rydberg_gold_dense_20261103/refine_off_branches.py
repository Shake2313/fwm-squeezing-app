"""Independent, bounded finer-control checks of three completed OFF map cells.

No campaign inputs or checkpoints are modified. Both profiled optical-detuning
branches receive their own local rectangle, and each branch is rescored over
the union of actually evaluated controls at that cell.
"""
from __future__ import annotations

import os

os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"

import csv
import hashlib
import itertools
import json
from pathlib import Path
import sys
import time

import numpy as np

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parents[1]
sys.path.insert(0, str(ROOT))

from dense_study import evaluate_control

CASES = [(20., 30.), (35., 37.5), (23., 30.)]
OPACITY = .001


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def complete_off_snapshot():
    """Read only complete lines from an append-only campaign checkpoint."""
    records = {}
    lines = (ROOT / "cells.jsonl").read_bytes().splitlines(keepends=True)
    ignored_partial_tail = False
    for index, line in enumerate(lines):
        try:
            record = json.loads(line)
        except json.JSONDecodeError:
            if index == len(lines) - 1 and not line.endswith(b"\n"):
                ignored_partial_tail = True
                break
            raise RuntimeError("Invalid complete/non-final checkpoint line")
        if record["task"][0] == "off":
            records[tuple(record["task"])] = record
    expected = {("off", float(t), float(pc)) for t in range(15, 36)
                for pc in np.arange(20, 40.0001, 2.5)}
    if set(records) != expected:
        raise RuntimeError(f"Completed OFF map required: {len(records)}/189 unique cells")
    snapshot = [records[key] for key in sorted(records)]
    canonical = json.dumps(snapshot, sort_keys=True, ensure_ascii=False,
                           separators=(",", ":"), allow_nan=False).encode("utf-8")
    manifest = dict(complete_off_cells=len(snapshot),
        canonicalization="Sorted by (mode,T,Pc); json.dumps sort_keys=True,separators=(',',':'),ensure_ascii=False",
        canonical_sha256=hashlib.sha256(canonical).hexdigest(),
        canonical_bytes=len(canonical), ignored_incomplete_appended_tail=ignored_partial_tail,
        scope="Only the completed 189 OFF cell records; later thermal appends do not alter this snapshot hash")
    return records, manifest


def local_grid(center, radius, step, bounds, protected=()):
    values = {round(float(value), 10) for value in np.arange(center - radius,
              center + radius + step / 3, step) if bounds[0] - 1e-10 <= value <= bounds[1] + 1e-10}
    values.add(float(center))
    for value in protected:
        if abs(value - center) <= radius + 1e-10 and bounds[0] <= value <= bounds[1]:
            values.add(float(value))
    return sorted(values)


def control(row):
    return dict(probe_power_uw=row["probe_power_uw"], lo_rabi_mhz=row["lo_rabi_mhz"])


def main():
    started = time.perf_counter()
    cells, snapshot = complete_off_snapshot()
    all_rows, cases = [], []
    for temperature, coupling in CASES:
        record = cells[("off", temperature, coupling)]
        original = next(row for row in record["branchgaps"] if row["opacity"] == OPACITY)
        cache, grids = {}, {}
        for branch in ("central", "sideband"):
            pp = original[branch + "_probe_power_uw"]
            lo = original[branch + "_lo_rabi_mhz"]
            probes = local_grid(pp, .5, .25, (2., 10.))
            los = local_grid(lo, .125, .0625, (2., 6.), protected=(3.7,))
            points = list(itertools.product(probes, los))
            grids[branch] = dict(probe_power_uw=probes, lo_rabi_mhz=los,
                                 control_points=len(points))
            for probe, rabi in points:
                key = (probe, rabi)
                if key not in cache:
                    cache[key] = evaluate_control("off", temperature, coupling,
                        probe, rabi, detuning_step=.005, off_target_step=.05)
        rows = [row for key in sorted(cache) for row in cache[key]]
        all_rows.extend(rows)
        nominal = [row for row in rows if row["opacity"] == OPACITY]
        comparisons = []
        for branch in ("central", "sideband"):
            quantity = branch + "_sensitivity_nv_cm_sqrt_hz"
            winner = min(nominal, key=lambda row: row[quantity])
            old = original[quantity]
            own_points = set(itertools.product(grids[branch]["probe_power_uw"],
                                               grids[branch]["lo_rabi_mhz"]))
            own = min((row for row in nominal if (row["probe_power_uw"], row["lo_rabi_mhz"])
                       in own_points), key=lambda row: row[quantity])
            comparisons.append(dict(branch=branch, old_sensitivity=old,
                finer_sensitivity=winner[quantity], relative_change=winner[quantity] / old - 1,
                relative_improvement=1 - winner[quantity] / old,
                old_control=dict(probe_power_uw=original[branch + "_probe_power_uw"],
                    lo_rabi_mhz=original[branch + "_lo_rabi_mhz"],
                    detuning_mhz=original[branch + "_detuning_mhz"]),
                finer_control=dict(**control(winner), detuning_mhz=winner[branch + "_detuning_mhz"]),
                own_rectangle_best_sensitivity=own[quantity],
                minimum_scope="Branch-constrained minimum over union of both local control rectangles"))
        branch_values = {item["branch"]: item["finer_sensitivity"] for item in comparisons}
        old_gap = original["sideband_over_central_minus_one"]
        fine_gap = branch_values["sideband"] / branch_values["central"] - 1
        old_winner = "central" if old_gap >= 0 else "sideband"
        fine_winner = "central" if fine_gap >= 0 else "sideband"
        case = dict(temperature_c=temperature, coupling_power_mw=coupling,
            opacity=OPACITY, grids=grids, actual_unique_control_evaluations=len(cache),
            branch_comparisons=comparisons, old_sideband_over_central_minus_one=old_gap,
            finer_sideband_over_central_minus_one=fine_gap,
            absolute_relative_branch_margin=abs(fine_gap),
            old_winning_branch=old_winner, finer_winning_branch=fine_winner,
            branch_selection_changed=old_winner != fine_winner,
            old_best_sensitivity=min(original[prefix + "_sensitivity_nv_cm_sqrt_hz"]
                                     for prefix in ("central", "sideband")),
            finer_best_sensitivity=min(branch_values.values()))
        cases.append(case)
        print(json.dumps(case, ensure_ascii=False), flush=True)
    maximum = max(abs(item["relative_change"]) for case in cases for item in case["branch_comparisons"])
    report = dict(origin="⟦AI⟧ Independent OFF local branch refinement", status="PASS" if maximum <= .005 else "REQUIRES_REVIEW",
        scope="Selected local control rectangles only; no new global optimum claim; original campaign untouched",
        settings=dict(cases=CASES, opacity=OPACITY, probe_radius_uw=.5,
            probe_step_uw=.25, lo_radius_mhz=.125, lo_step_mhz=.0625,
            detuning_step_mhz=.005, off_global_detuning_target_step_mhz=.05),
        off_snapshot=snapshot, cases=cases,
        maximum_absolute_relative_branch_change=maximum,
        minimum_absolute_relative_branch_margin=min(case["absolute_relative_branch_margin"] for case in cases),
        branch_selection_changed_cases=sum(case["branch_selection_changed"] for case in cases),
        actual_control_evaluations=sum(case["actual_unique_control_evaluations"] for case in cases),
        scored_rows=len(all_rows), elapsed_s=time.perf_counter() - started,
        interpretation="PASS means selected branches stable within0.5%; a near-tie branch identity may still change. Extents describe evaluated samples, not a filled acceptable region.",
        source_sha256={name: sha(ROOT / name) for name in ("refine_off_branches.py", "dense_study.py", "config.json")})
    with (ROOT / "off_branch_refinement.csv").open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(all_rows[0]))
        writer.writeheader()
        writer.writerows(all_rows)
    report["artifact_sha256"] = {"off_branch_refinement.csv": sha(ROOT / "off_branch_refinement.csv")}
    (ROOT / "off_branch_refinement.json").write_text(json.dumps(report, ensure_ascii=False,
        indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({key: value for key, value in report.items() if key not in ("cases",)},
                     ensure_ascii=False), flush=True)
    return 0 if report["status"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
