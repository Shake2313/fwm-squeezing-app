"""Cheap frozen-source inflow diagnostics; no atomic solve or certification.

Run from the repository root with Python -B. A rerun needs --output pointing to
a new JSON in this audit directory because historical output is never replaced.
The source ZIP is imported directly, without extracting or changing its files.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import platform
import sys
import time
import zipfile


def sha256(raw):
    return hashlib.sha256(raw).hexdigest()


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"),
                      ensure_ascii=False, allow_nan=False).encode("utf-8")


def unique_object(pairs):
    result = {}
    for name, value in pairs:
        if name in result:
            raise ValueError(f"Duplicate JSON key: {name}")
        result[name] = value
    return result


def main():
    started = time.perf_counter()
    sys.dont_write_bytecode = True
    probe = Path(__file__).resolve()
    root = probe.parents[3]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=probe.with_suffix(".json"))
    args = parser.parse_args()
    output = args.output.resolve()
    if output.parent != probe.parent or output.suffix != ".json":
        raise ValueError("Output must be a new JSON inside this audit directory")
    if output.exists():
        raise FileExistsError(output)

    campaign = root / "docs/grand_challenge/thermal_campaign_v2"
    plan_path, bundle = campaign / "plan.json", campaign / "sources.zip"
    captured = {path: path.read_bytes() for path in (probe, plan_path, bundle)}
    plan = json.loads(captured[plan_path].decode("utf-8-sig"),
                      object_pairs_hook=unique_object)
    plan_payload = {key: value for key, value in plan.items()
                    if key != "record_sha256"}
    if sha256(canonical(plan_payload)) != plan["record_sha256"]:
        raise ValueError("Campaign plan seal does not match its payload")
    if sha256(captured[bundle]) != plan["source_bundle_sha256"]:
        raise ValueError("Source ZIP differs from the sealed campaign plan")
    with zipfile.ZipFile(bundle) as archive:
        entries = plan["source_manifest"]["files"]
        if sorted(archive.namelist()) != sorted(row["path"] for row in entries):
            raise ValueError("ZIP source inventory differs from plan manifest")
        for row in entries:
            raw = archive.read(row["path"])
            if len(raw) != row["raw_bytes"] or sha256(raw) != row["raw_sha256"]:
                raise ValueError(f"ZIP source differs from manifest: {row['path']}")

    # A fresh interpreter and first-position ZIP ensure that repository modules
    # cannot silently replace the frozen geometry implementation.
    if any(name == "gabes" or name.startswith("gabes.") for name in sys.modules):
        raise ValueError("Run in a fresh interpreter; gabes already imported")
    sys.path.insert(0, str(bundle))
    import numpy as np
    import scipy
    from gabes.quantum.inflow import maxwell_box_inflow

    frozen_prefix = str(bundle) + "\\" if sys.platform == "win32" else str(bundle) + "/"
    loaded_sources = []
    for name, module in list(sys.modules.items()):
        if name == "gabes" or name.startswith("gabes."):
            location = getattr(module, "__file__", "") or ""
            if not location.startswith(frozen_prefix):
                raise ValueError(f"Loaded gabes source outside frozen ZIP: {name}")
            relative = location[len(frozen_prefix):].replace("\\", "/")
            identity = next(row for row in entries if row["path"] == relative)
            loaded_sources.append({"module": name, "zip_member": relative,
                                   "raw_sha256": identity["raw_sha256"]})

    model = plan["model"]
    powers = [2, 3, 4, 5, 6, 8, 10, 14]
    seeds = plan["seeds"]
    cutoff_s = 4e-6
    sigma2 = 1.380649e-23 * model["temperature_K"] / model["mass_kg"]
    declared = {(grid["power"], grid["seed"]): grid for grid in plan["grids"]}
    rows, declaration_checks = [], []
    geometry_started = time.perf_counter()
    for power in powers:
        for seed in seeds:
            inflow = maxwell_box_inflow(
                model["lower_corner_m"], model["upper_corner_m"],
                temperature_K=model["temperature_K"], mass_kg=model["mass_kg"],
                density_m3=model["inputs_SI"]["number_density_m3"],
                points_per_face_power=power, seed=seed,
                source="Audit geometry only: sealed conditional fixture; no atomic solve")
            if (power, seed) in declared:
                expected = declared[power, seed]["paths"]
                actual = [{"index": index, "physical_path": {
                    "entry_position_m": inflow.entry_position_m[index].tolist(),
                    "velocity_m_s": inflow.velocity_m_s[index].tolist(),
                    "residence_time_s": float(inflow.residence_time_s[index])},
                    "rate_s_inverse": float(inflow.rate_s_inverse[index])}
                    for index in range(len(inflow.rate_s_inverse))]
                if actual != expected:
                    raise ValueError(f"Generated paths differ from declaration: {power}, {seed}")
                declaration_checks.append({"grid": [power, seed],
                                           "path_count": len(actual), "exact_match": True})
            tail = inflow.residence_time_s >= cutoff_s
            moments = inflow.occupation_moments()
            rows.append({
                "power": power, "seed": seed, "points_per_face": 2**power,
                "path_count": len(inflow.rate_s_inverse),
                "minimum_residence_us": float(inflow.residence_time_s.min() * 1e6),
                "maximum_residence_us": float(inflow.residence_time_s.max() * 1e6),
                "occupancy_over_nV": inflow.mean_occupancy / inflow.equilibrium_atom_number,
                "occupation_velocity_second_over_sigma2":
                    (np.diag(moments["velocity_second_m2_s2"]) / sigma2).tolist(),
                "tail_path_count": int(tail.sum()),
                "tail_path_count_by_face": np.bincount(
                    inflow.face_index[tail], minlength=6).tolist(),
                "tail_arrival_fraction": float(inflow.rate_s_inverse[tail].sum()
                                                / inflow.total_arrival_rate_s_inverse),
                "tail_occupancy_over_nV": float(
                    inflow.rate_s_inverse[tail] @ inflow.residence_time_s[tail]
                    / inflow.equilibrium_atom_number),
            })
    geometry_elapsed = time.perf_counter() - geometry_started
    actual_environment = {"python": platform.python_version(), "numpy": np.__version__,
                          "scipy": scipy.__version__, "platform": platform.platform()}
    environment_matches = {key: actual_environment[key] == plan["environment"][key]
                           for key in actual_environment}
    for path, raw in captured.items():
        if path.read_bytes() != raw:
            raise ValueError(f"Consumed input changed during the probe: {path}")
    result = {
        "schema": "gabes-grand-challenge-geometry-audit-v1",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "command_from_repository_root":
            "python -B docs/grand_challenge/audit_2026_10_02/geometry_probe.py",
        "rerun": "Use --output NEW.json within the same audit directory; existing outputs are exclusive.",
        "campaign_plan": "docs/grand_challenge/thermal_campaign_v2/plan.json",
        "campaign_plan_raw_sha256": sha256(captured[plan_path]),
        "campaign_plan_record_sha256": plan["record_sha256"],
        "campaign_source_identity": plan["source_identity"],
        "source_bundle_raw_sha256": sha256(captured[bundle]),
        "source_manifest_raw_entries_verified": len(entries),
        "probe_raw_sha256": sha256(captured[probe]),
        "loaded_frozen_sources": sorted(loaded_sources, key=lambda row: row["module"]),
        "geometry_configuration": {key: model[key] for key in (
            "lower_corner_m", "upper_corner_m", "temperature_K", "mass_kg")},
        "number_density_m3": model["inputs_SI"]["number_density_m3"],
        "powers": powers, "seeds": seeds, "tail_residence_cutoff_s": cutoff_s,
        "declared_grid_exact_path_checks": declaration_checks,
        "selected_grid_estimates": rows,
        "environment": actual_environment,
        "environment_matches_sealed_plan": environment_matches,
        "geometry_seconds": geometry_elapsed,
        "total_seconds_before_output_write": time.perf_counter() - started,
        "new_atomic_solves": 0, "atomic_cache_packets_read": 0,
        "source_extracted": False, "historical_inputs_modified": False,
        "thermal_ensemble_converged": False, "physical_optical_prediction": False,
        "scope": "Finite frozen-source Maxwell inflow geometry diagnostics only.",
        "caveats": [
            "Arrival fractions and occupancy/nV contributions are numerical geometry estimates, not atomic noise fractions or tail error bounds.",
            "High-power geometry moments do not certify atomic spectra, the ensemble gate, boundary physics, an optical channel, or experimental predictions.",
            "Uncomputed declared paths are inspected geometrically only; this does not complete their five-resolution atomic evidence or grid audit.",
            "The p14 estimates are an inexpensive tail reconnaissance, not a rigorous continuous-measure reference or confidence interval.",
            "Maximum sampled residence is a runtime warning, not a physical residence cutoff; the thermal measure retains unbounded slow-atom tails.",
            "Input hashes detect byte changes and bind provenance; they do not authenticate the original scientific assumptions.",
        ],
    }
    result["record_sha256"] = sha256(canonical(result))
    with output.open("x", encoding="utf-8", newline="\n") as handle:
        json.dump(result, handle, ensure_ascii=False, allow_nan=False, indent=2)
        handle.write("\n")
    print(json.dumps({"output": str(output), "record_sha256": result["record_sha256"],
                      "geometry_seconds": geometry_elapsed, "new_atomic_solves": 0}))


if __name__ == "__main__":
    main()
