"""Actual fixed-domain dense maps with adaptive sampled control refinement.

Analysis-only scan subclass; shared production numerical code remains unchanged.
Every saved map coordinate is solved. Hidden controls use two branch-specific
local lattice descents, not a full five-dimensional Cartesian exhaustive search.
"""
from __future__ import annotations

import os

os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"

import argparse
import csv
import hashlib
import json
import platform
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from functools import lru_cache
from pathlib import Path

import numpy as np
import scipy

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parents[1]
PRIOR = REPO / "analysis/rydberg_gold_map_20261103"
sys.path.insert(0, str(REPO))

from gabes import affine_doppler, constants, core
from gabes import rydberg_electrometry as electrometry
from gabes.schemes.rydberg import MHZ, SUPERHET_CLASS_CHUNK, RydbergEITScheme

TEMPERATURES = [float(t) for t in range(15, 36)]
COUPLING_POWERS = np.arange(20, 40.0001, 2.5).tolist()
OPACITIES = [.0005, .001, .0015]
MODES = ["off", "on"]
OLD_TEMPERATURES = [15., 20., 25., 30., 33., 35.]
OLD_COUPLING_POWERS = [20., 30., 40.]
OLD_PROBE_POWERS = [2., 4., 6., 8., 10.]
OLD_LO_RABIS = [2., 3., 3.7, 4.5, 6.]
PROBE_BOUNDS = (2., 10.)
LO_BOUNDS = (2., 6.)
PROBE_WINDOW = 6.
BRANCH_CUT = .75
CONTROL_STAGES = [(1., .25), (.5, .125)]
MAX_ITERATIONS = 12
IMPROVEMENT_RTOL = 1e-9


def write_json(path, value):
    path.write_text(json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False)
                    + "\n", encoding="utf-8")


@lru_cache(maxsize=1)
def baseline_parameters():
    return json.loads((PRIOR / "resolved_parameters.json").read_text(
        encoding="utf-8"))["baseline_parameters"]


@lru_cache(maxsize=1)
def prior_results():
    return json.loads((PRIOR / "results.json").read_text(encoding="utf-8"))


class FixedScanScheme(RydbergEITScheme):
    def __init__(self, detuning_step=.01):
        self.detuning_step = float(detuning_step)
        self.fixed_velocity = None
        self.fixed_weights = None
        super().__init__()

    def _scan(self, params):
        intervals = round(2 * PROBE_WINDOW / self.detuning_step)
        if not np.isclose(intervals * self.detuning_step, 2 * PROBE_WINDOW):
            raise ValueError("Detuning step must divide the 12 MHz domain")
        return np.linspace(-PROBE_WINDOW, PROBE_WINDOW, intervals + 1) * MHZ

    def _superheterodyne_atomic_phasors(self, raw, params, x, transmission, static_indices=None):
        if params["doppler"] != "on":
            return super()._superheterodyne_atomic_phasors(
                raw, params, x, transmission, static_indices=static_indices)
        indices = np.asarray(static_indices, dtype=int)
        gd = (raw["rydberg_dephasing_mhz"] + raw["temperature_dephasing_mhz"]
              + raw["density_dephasing_mhz"] + raw["transit_mhz"]) * MHZ
        atom = self._atom(gd, raw["rf_dephasing_mhz"] * MHZ)
        omega_if = float(params["if_khz"]) * 1e3 * 2 * np.pi
        operator = self._absorption_quadrature_operator()
        if self.fixed_velocity is None:
            affine = affine_doppler.AffineScan(atom,
                lambda s: self._hamiltonian_at_probe_detuning(raw, params, s),
                (1, 0), raw["probe_rabi_mhz"] * MHZ)
            # Same pole-resolving quadrature as production, selected from the
            # entire fixed scan once per control, independent of response subset.
            self.fixed_velocity, self.fixed_weights = affine.quadrature_classes(
                raw["scan"], raw["T"], constants.MASS_85RB, raw["k_vec"])
        kv = raw["k_vec"] * self.fixed_velocity
        phasors = np.empty(indices.size, dtype=complex)
        for slot, index in enumerate(indices):
            H = self._hamiltonian_at_probe_detuning(raw, params, raw["scan"][index])
            base = core.build_liouvillian(H, atom)
            classes = np.empty(kv.size, dtype=complex)
            for first in range(0, kv.size, SUPERHET_CLASS_CHUNK):
                part = slice(first, first + SUPERHET_CLASS_CHUNK)
                liouvillian = base[None, :, :] + kv[part, None, None] * atom.S_v[None, :, :]
                classes[part] = electrometry.weak_signal_response_from_liouvillian(
                    liouvillian, 4, omega_if, signal_transition=(2, 3)
                ).real_observable_phasor_per_angular_rabi(operator)
            phasors[slot] = electrometry.coherent_weighted_average(classes, self.fixed_weights)
        return indices, phasors, "fixed full-domain velocity grid; coherent Doppler average"


def target_indices(x, step, centers=(), radius=0.):
    targets = list(np.arange(-PROBE_WINDOW, PROBE_WINDOW + step / 3, step))
    for center in centers:
        left, right = max(-PROBE_WINDOW, center - radius), min(PROBE_WINDOW, center + radius)
        targets.extend(np.arange(left, right + step / 3, step))
        targets.append(center)
    return np.unique(np.array([int(np.argmin(np.abs(x - value))) for value in targets], dtype=int))


def local_indices(x, centers, radius, step):
    targets = []
    for center in centers:
        left, right = max(-PROBE_WINDOW, center - radius), min(PROBE_WINDOW, center + radius)
        targets.extend(np.arange(left, right + step / 3, step))
        targets.append(center)
    return np.unique(np.array([int(np.argmin(np.abs(x - value))) for value in targets], dtype=int))


def atomic_at_indices(scheme, raw, params, x, transmission, indices):
    if params["doppler"] == "on":
        return scheme._superheterodyne_atomic_phasors(
            raw, params, x, transmission, static_indices=indices)
    # Production Doppler-off API intentionally ignores static_indices. Subset
    # its scan explicitly, then map its local indices back to the fixed grid.
    subset = dict(raw, scan=raw["scan"][indices])
    atomic = scheme._superheterodyne_atomic_phasors(
        subset, params, x[indices], transmission[indices])
    return indices[atomic[0]], atomic[1], "explicit physical detuning targets"


def boundary_flags(row):
    flags = []
    for key, bounds in (("probe_power_uw", PROBE_BOUNDS), ("lo_rabi_mhz", LO_BOUNDS),
                        ("temperature_c", (15., 35.)), ("coupling_power_mw", (20., 40.))):
        if np.isclose(row[key], bounds[0], rtol=0, atol=1e-10):
            flags.append(key + "_lower")
        if np.isclose(row[key], bounds[1], rtol=0, atol=1e-10):
            flags.append(key + "_upper")
    if abs(row["probe_detuning_mhz"]) >= PROBE_WINDOW - row["detuning_step_mhz"] * 1.01:
        flags.append("probe_detuning_boundary")
    return "|".join(flags)


def score_atomic(scheme, raw, params, x, atomic):
    scores = []
    for scale in OPACITIES:
        scenario = dict(params, effective_opacity_participation=scale)
        _, transmission, _ = scheme._transmission(raw["chi_bar"], raw, scenario)
        result = scheme._superheterodyne_readout(
            raw, scenario, x, transmission, atomic_response=atomic)
        if result is None:
            raise RuntimeError("No finite-IF sensitivity")
        scores.append((scale, transmission, result))
    return scores


def evaluate_control(mode, t, pc, pp, lo, detuning_step=.01,
                     thermal_target_step=.5, thermal_refine_step=.05,
                     off_target_step=.1, off_refine_radius=.12,
                     thermal_refine_radius=.3):
    """Return three opacity-score rows for one actual atomic calculation.

    A fixed full-domain velocity grid permits response reuse across refinement.
    Both detuning signs are retained: finite-IF response is not mirror symmetric.
    """
    started = time.perf_counter()
    scheme = FixedScanScheme(detuning_step)
    params = dict(baseline_parameters(), doppler=mode, temp_c=float(t),
                  coupling_power_mw=float(pc), probe_power_uw=float(pp),
                  lo_rabi_mhz=float(lo), mw_detuning_mhz=0.)
    raw = scheme.compute(params)
    x, transmission, _ = scheme._transmission(raw["chi_bar"], raw, params)
    target_step = thermal_target_step if mode == "on" else off_target_step
    indices = target_indices(x, target_step)
    indices = np.union1d(indices, scheme._static_slope_indices(x, transmission)).astype(int)
    atomic = atomic_at_indices(scheme, raw, params, x, transmission, indices)
    initial = score_atomic(scheme, raw, params, x, atomic)
    centers = []
    for _, _, result in initial:
        values = result["total_field_asd_v_m_per_sqrt_hz"]
        detunings = x[atomic[0]]
        for mask in (abs(detunings) < BRANCH_CUT, abs(detunings) >= BRANCH_CUT):
            positions = np.flatnonzero(mask & np.isfinite(values))
            if positions.size:
                centers.append(float(detunings[positions[np.argmin(values[positions])]]))
    radius = thermal_refine_radius if mode == "on" else off_refine_radius
    refine_step = thermal_refine_step if mode == "on" else detuning_step
    refined_indices = np.union1d(indices, local_indices(x, centers, radius, refine_step)).astype(int)
    extra_indices = np.setdiff1d(refined_indices, indices)
    if extra_indices.size:
        extra = atomic_at_indices(scheme, raw, params, x, transmission, extra_indices)
        combined_indices = np.r_[atomic[0], extra[0]]
        order = np.argsort(combined_indices)
        refined_atomic = (combined_indices[order], np.r_[atomic[1], extra[1]][order],
                          "global physical targets plus branch-local refinement")
    else:
        refined_atomic = atomic
    final = score_atomic(scheme, raw, params, x, refined_atomic)
    rows = []
    for scale, transmission, result in final:
        sensitivity = result["sensitivity"]
        row = dict(doppler=mode, temperature_c=float(t), coupling_power_mw=float(pc),
            probe_power_uw=float(pp), lo_rabi_mhz=float(lo), opacity=scale,
            probe_detuning_mhz=float(result["optimum_detuning_mhz"]),
            sensitivity_nv_cm_sqrt_hz=float(sensitivity.total_field_asd_nv_cm_per_sqrt_hz),
            psn_nv_cm_sqrt_hz=float(sensitivity.psn_field_asd_nv_cm_per_sqrt_hz),
            transmission=float(result["optimum_transmission"]),
            current_response_a_per_v_m=float(result["current_responsivity_a_per_v_m"]),
            detected_power_w=float(result["signal_detector_power_w"]),
            coarse_candidate_count=int(indices.size), candidate_count=int(refined_atomic[0].size),
            detuning_step_mhz=float(detuning_step),
            actual_probe_search="fixed physical targets plus branch-local finite-IF refinement",
            control_key=f"{float(pp):g}|{float(lo):g}",
            mirror_symmetry_used=False,
            mirror_degenerate=False,
            quadrature_class_count=0 if scheme.fixed_velocity is None else int(scheme.fixed_velocity.size),
            quadrature_dv_m_s=0. if scheme.fixed_velocity is None else float(np.diff(scheme.fixed_velocity)[0]),
            quadrature_velocity_limit_m_s=0. if scheme.fixed_velocity is None else float(scheme.fixed_velocity[-1]))
        values = result["total_field_asd_v_m_per_sqrt_hz"] * 1e7
        detunings = x[refined_atomic[0]]
        for branch, mask in (("central", abs(detunings) < BRANCH_CUT),
                             ("sideband", abs(detunings) >= BRANCH_CUT)):
            positions = np.flatnonzero(mask & np.isfinite(values))
            if not positions.size:
                raise RuntimeError("A detuning branch has no finite candidates")
            position = positions[np.argmin(values[positions])]
            row[branch + "_sensitivity_nv_cm_sqrt_hz"] = float(values[position])
            row[branch + "_detuning_mhz"] = float(detunings[position])
        row["winning_branch"] = "central" if abs(row["probe_detuning_mhz"]) < BRANCH_CUT else "sideband"
        row["boundary_flags"] = boundary_flags(row)
        row["elapsed_s"] = time.perf_counter() - started
        if not np.isfinite(row["sensitivity_nv_cm_sqrt_hz"]) or row["sensitivity_nv_cm_sqrt_hz"] <= 0:
            raise RuntimeError("Nonpositive or nonfinite sensitivity")
        if not 0 < row["transmission"] <= 1 + 1e-10:
            raise RuntimeError("Nonpassive transmission")
        rows.append(row)
    return rows


def control_key(pp, lo):
    return (round(float(pp), 8), round(float(lo), 8))


def cell_key(task):
    return json.dumps(task, separators=(",", ":"))


def axis_neighbors(value, step, bounds, protected=()):
    # Consecutive canonical lattice points on each side; exact Gold LO is
    # retained as an extra physical point rather than shifting the lattice.
    lattice = np.arange(bounds[0], bounds[1] + step / 3, step).tolist() + list(protected)
    lower = [point for point in lattice if point < value - 1e-9]
    upper = [point for point in lattice if point > value + 1e-9]
    return ([max(lower)] if lower else []) + ([min(upper)] if upper else [])


def nearby_seeds(mode, t, pc):
    temps = sorted(OLD_TEMPERATURES, key=lambda v: (abs(v-t), v))[:2]
    powers = sorted(OLD_COUPLING_POWERS, key=lambda v: (abs(v-pc), v))[:2]
    return {(row["probe_power_uw"], row["lo_rabi_mhz"])
            for row in prior_results()["optima"] if row["doppler"] == mode
            and row["temperature_c"] in temps and row["coupling_power_mw"] in powers}


def search_cell(task):
    """Run actual local lattice searches for both optical-detuning branches."""
    mode, t, pc = task
    started = time.perf_counter()
    cache = {}

    def evaluate(pp, lo):
        key = control_key(pp, lo)
        if key not in cache:
            cache[key] = evaluate_control(mode, t, pc, *key)
        return cache[key]

    seeds = ({(10., 2.), (4., 3.7), (4., 4.5), (6., 3.7)} if mode == "off"
             else {(10., 4.5), (10., 6.), (6., 3.7), (2., 2.)})
    seeds.update(nearby_seeds(mode, t, pc))
    if t in OLD_TEMPERATURES and pc in OLD_COUPLING_POWERS:
        seeds.update((pp, lo) for pp in OLD_PROBE_POWERS for lo in OLD_LO_RABIS)
    for seed in sorted(seeds):
        evaluate(*seed)

    paths = {}
    for oi, opacity in enumerate(OPACITIES):
        for branch in ("central", "sideband"):
            quantity = branch + "_sensitivity_nv_cm_sqrt_hz"
            paths[(oi, branch)] = min(cache, key=lambda key: cache[key][oi][quantity])
    convergence = []
    for pp_step, lo_step in CONTROL_STAGES:
        for (oi, branch), position in list(paths.items()):
            quantity = branch + "_sensitivity_nv_cm_sqrt_hz"
            converged = False
            for iteration in range(1, MAX_ITERATIONS + 1):
                pp, lo = position
                neighbors = [control_key(value, lo) for value in axis_neighbors(pp, pp_step, PROBE_BOUNDS)]
                neighbors += [control_key(pp, value) for value in axis_neighbors(
                    lo, lo_step, LO_BOUNDS, protected=(3.7,))]
                for neighbor in neighbors:
                    evaluate(*neighbor)
                best = min([position] + neighbors, key=lambda key: cache[key][oi][quantity])
                if cache[best][oi][quantity] >= cache[position][oi][quantity] * (1 - IMPROVEMENT_RTOL):
                    converged = True
                    break
                position = best
            paths[(oi, branch)] = position
            convergence.append(dict(opacity=OPACITIES[oi], branch=branch,
                probe_step_uw=pp_step, lo_step_mhz=lo_step, iterations=iteration,
                converged=converged, final_control=list(position),
                sensitivity=cache[position][oi][quantity]))

    # Later opacity/branch descents may add a better point for an earlier path.
    # Close every final global branch winner under its final cross stencil.
    for iteration in range(1, MAX_ITERATIONS + 1):
        changed = False
        for oi, opacity in enumerate(OPACITIES):
            for branch in ("central", "sideband"):
                quantity = branch + "_sensitivity_nv_cm_sqrt_hz"
                position = min(cache, key=lambda key: cache[key][oi][quantity])
                pp, lo = position
                neighbors = [control_key(value, lo) for value in axis_neighbors(pp, .5, PROBE_BOUNDS)]
                neighbors += [control_key(pp, value) for value in axis_neighbors(lo, .125, LO_BOUNDS, (3.7,))]
                value_before = cache[position][oi][quantity]
                for neighbor in neighbors:
                    evaluate(*neighbor)
                winner = min(cache, key=lambda key: cache[key][oi][quantity])
                changed |= cache[winner][oi][quantity] < value_before * (1-IMPROVEMENT_RTOL)
        if not changed:
            break
    convergence.append(dict(opacity="all", branch="final-global-branch-closure",
        probe_step_uw=.5, lo_step_mhz=.125, iterations=iteration,
        converged=not changed, final_control=[], sensitivity=0.))

    rows = [row for key in sorted(cache) for row in cache[key]]
    optima, branchgaps, epsilon_bands = [], [], []
    for oi, opacity in enumerate(OPACITIES):
        candidates = [cache[key][oi] for key in sorted(cache)]
        best = dict(min(candidates, key=lambda row: row["sensitivity_nv_cm_sqrt_hz"]))
        optima.append(best)
        central = min(candidates, key=lambda row: row["central_sensitivity_nv_cm_sqrt_hz"])
        sideband = min(candidates, key=lambda row: row["sideband_sensitivity_nv_cm_sqrt_hz"])
        c, s = central["central_sensitivity_nv_cm_sqrt_hz"], sideband["sideband_sensitivity_nv_cm_sqrt_hz"]
        branchgaps.append(dict(doppler=mode, temperature_c=t, coupling_power_mw=pc,
            opacity=opacity, central_sensitivity_nv_cm_sqrt_hz=c,
            sideband_sensitivity_nv_cm_sqrt_hz=s, sideband_over_central_minus_one=s/c-1,
            central_probe_power_uw=central["probe_power_uw"], central_lo_rabi_mhz=central["lo_rabi_mhz"],
            central_detuning_mhz=central["central_detuning_mhz"],
            sideband_probe_power_uw=sideband["probe_power_uw"], sideband_lo_rabi_mhz=sideband["lo_rabi_mhz"],
            sideband_detuning_mhz=sideband["sideband_detuning_mhz"]))
        for epsilon in (.01, .03):
            accepted = [row for row in candidates if row["sensitivity_nv_cm_sqrt_hz"] <=
                        best["sensitivity_nv_cm_sqrt_hz"] * (1 + epsilon)]
            epsilon_bands.append(dict(doppler=mode, temperature_c=t, coupling_power_mw=pc,
                opacity=opacity, epsilon=epsilon, saved_control_count=len(accepted),
                probe_power_min_uw=min(row["probe_power_uw"] for row in accepted),
                probe_power_max_uw=max(row["probe_power_uw"] for row in accepted),
                lo_rabi_min_mhz=min(row["lo_rabi_mhz"] for row in accepted),
                lo_rabi_max_mhz=max(row["lo_rabi_mhz"] for row in accepted),
                represented_branches="|".join(sorted(set(row["winning_branch"] for row in accepted))),
                scope="Bounds over evaluated samples; not a filled rectangular acceptable region"))
    return dict(task=list(task), key=cell_key(task), rows=rows, optima=optima,
        branchgaps=branchgaps, epsilon_bands=epsilon_bands, convergence=convergence,
        control_evaluations=len(cache), initial_seed_count=len(seeds), elapsed_s=time.perf_counter()-started)


def campaign_config():
    engine_files = sorted((REPO / "gabes").glob("*.py")) + [REPO / "gabes/schemes/rydberg.py"]
    inputs = [PRIOR / name for name in ("resolved_parameters.json", "results.json", "pilot.py")]
    source_hashes = {str(path.relative_to(REPO)).replace("\\", "/"):
                     hashlib.sha256(path.read_bytes()).hexdigest() for path in engine_files + inputs}
    source_hashes["dense_study.py"] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    environment = dict(python=sys.version, numpy=np.__version__, scipy=scipy.__version__,
                       platform=platform.platform(), blas_threads=1)
    config = dict(as_of="2026-11-03", temperature_c=TEMPERATURES, coupling_power_mw=COUPLING_POWERS,
        probe_bounds_uw=list(PROBE_BOUNDS), lo_bounds_mhz=list(LO_BOUNDS), opacities=OPACITIES,
        modes=MODES, probe_window_mhz=PROBE_WINDOW, detuning_step_mhz=.01,
        control_stages=[dict(probe_step_uw=pp, lo_step_mhz=lo) for pp, lo in CONTROL_STAGES],
        max_iterations=MAX_ITERATIONS, branch_cut_mhz=BRANCH_CUT,
        fixed=dict(cell_mm=50., beam_diameter_mm=.15, if_khz=40., rf_detuning_mhz=0.,
                   coupling_detuning_mhz=0., temperature_model="Linked"),
        search_scope="Actual dense T/Pc map; adaptive two-branch local control lattice, not full 5D Cartesian",
        detuning_symmetry="Both signs retained; fresh API check rejected finite-IF mirror shortcut before campaign freeze",
        thermal_quadrature="Once per control: production pole-resolving grid on full fixed scan, reused for new response candidates",
        provenance="AI-selected refinement of user-adopted Gold; same bounded domain; no sensitivity fit",
        code_sha256=source_hashes, environment=environment)
    config["signature"] = hashlib.sha256(json.dumps(config, sort_keys=True).encode()).hexdigest()
    return config


def read_checkpoint(path):
    records = {}
    if not path.exists():
        return records
    content = path.read_bytes()
    lines = content.splitlines(keepends=True)
    valid_end = 0
    for index, line in enumerate(lines):
        try:
            record = json.loads(line)
        except json.JSONDecodeError:
            if index != len(lines) - 1:
                raise RuntimeError("Corrupt non-final checkpoint line; original retained")
            # Preserve exact corrupt bytes before retaining only verified lines.
            backup = path.with_suffix(path.suffix + ".truncated_tail")
            if backup.exists() and backup.read_bytes() != line:
                raise RuntimeError("Existing checkpoint recovery backup differs")
            backup.write_bytes(line)
            path.write_bytes(content[:valid_end])
            print(f"Recovered truncated final checkpoint line; original tail: {backup.name}", flush=True)
            break
        records[record["key"]] = record
        valid_end += len(line)
    if records and not path.read_bytes().endswith(b"\n"):
        with path.open("ab") as stream:
            stream.write(b"\n")
    return records


def export_csv(path, rows):
    with path.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def export_results(config, records, tasks):
    selected = [records[cell_key(task)] for task in tasks]
    rows = [row for record in selected for row in record["rows"]]
    optima = [row for record in selected for row in record["optima"]]
    baselines = [row for row in rows if (row["temperature_c"], row["coupling_power_mw"],
        row["probe_power_uw"], row["lo_rabi_mhz"]) == (20., 30., 6., 3.7)]
    best = [min((row for row in optima if row["doppler"] == mode and row["opacity"] == .001),
                key=lambda row: row["sensitivity_nv_cm_sqrt_hz"]) for mode in MODES]
    failed = [dict(task=record["task"], **entry) for record in selected
              for entry in record["convergence"] if not entry["converged"]]
    summary = dict(config=config, baseline=baselines, best_nominal=best, optima=optima,
        branchgaps=[row for record in selected for row in record["branchgaps"]],
        epsilon_bands=[row for record in selected for row in record["epsilon_bands"]],
        map_coordinates_per_mode=len(TEMPERATURES)*len(COUPLING_POWERS),
        control_states=sum(record["control_evaluations"] for record in selected),
        saved_score_rows=len(rows), failed_convergence=failed,
        elapsed_worker_s=sum(record["elapsed_s"] for record in selected),
        numerical_status="Dense sampled map complete; independent checks separate",
        experimental_status="UNCALIBRATED conditional PSN model")
    write_json(ROOT / "results.json", summary)
    export_csv(ROOT / "control_points.csv", rows)
    export_csv(ROOT / "optima.csv", optima)
    export_csv(ROOT / "branchgaps.csv", summary["branchgaps"])
    export_csv(ROOT / "epsilon_bands.csv", summary["epsilon_bands"])
    write_json(ROOT / "source_manifest.json", dict(signature=config["signature"],
        code_sha256=config["code_sha256"], environment=config["environment"]))
    print(json.dumps(dict(control_states=summary["control_states"],
        best_nominal=best, failed_convergence=len(failed)), ensure_ascii=False, indent=2), flush=True)


def symmetry_diagnostic():
    """Fresh production API check of a proposed, rejected runtime shortcut."""
    cases = []
    for mode, pp, lo, t, pc in [("off", 6., 3.7, 20., 30.),
            ("off", 10., 2., 35., 40.), ("on", 6., 3.7, 20., 30.),
            ("on", 10., 4.5, 35., 40.)]:
        scheme = FixedScanScheme()
        params = dict(baseline_parameters(), doppler=mode, temp_c=t,
                      coupling_power_mw=pc, probe_power_uw=pp, lo_rabi_mhz=lo)
        raw = scheme.compute(params)
        x, transmission, _ = scheme._transmission(raw["chi_bar"], raw, params)
        targets = [-3., -1.5, -.25, .25, 1.5, 3.]
        indices = np.array([int(np.argmin(abs(x-value))) for value in targets])
        if mode == "on":
            atomic = RydbergEITScheme._superheterodyne_atomic_phasors(
                scheme, raw, params, x, transmission, static_indices=indices)
        else:
            atomic = atomic_at_indices(scheme, raw, params, x, transmission, indices)
        pairs = []
        for left in range(3):
            right = len(targets)-1-left
            negative, positive = atomic[1][left], atomic[1][right]
            pairs.append(dict(abs_detuning_mhz=abs(targets[left]),
                negative_transmission=float(transmission[indices[left]]),
                positive_transmission=float(transmission[indices[right]]),
                negative_response_phasor=[float(negative.real), float(negative.imag)],
                positive_response_phasor=[float(positive.real), float(positive.imag)],
                response_magnitude_relative_difference=float(abs(abs(negative)-abs(positive))/abs(negative))))
        cases.append(dict(doppler=mode, temperature_c=t, coupling_power_mw=pc,
                          probe_power_uw=pp, lo_rabi_mhz=lo, pairs=pairs))
    result = dict(origin="AI pre-freeze numerical diagnostic", status="SHORTCUT_REJECTED",
        candidate_assumption="Finite-IF |response| mirror symmetric under probe detuning sign reversal",
        conclusion="Static transmission is mirror symmetric, but finite-IF weak-SIG response is not; retain both signs",
        campaign_impact="No campaign values use the rejected shortcut", cases=cases)
    write_json(ROOT / "symmetry_check.json", result)
    print(json.dumps(result, ensure_ascii=False, indent=2), flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--benchmark", action="store_true", help="Two controls only; no map artifacts written")
    parser.add_argument("--benchmark-cell", action="store_true", help="Two Gold-coordinate searches only")
    parser.add_argument("--symmetry-check", action="store_true", help="Check rejected detuning mirror shortcut with fresh APIs")
    args = parser.parse_args()
    if args.symmetry_check:
        symmetry_diagnostic()
        return
    if args.benchmark:
        for mode in MODES:
            started = time.perf_counter()
            rows = evaluate_control(mode, 20., 30., 6., 3.7)
            print(json.dumps(dict(mode=mode, elapsed_s=time.perf_counter()-started, rows=rows),
                             ensure_ascii=False, indent=2), flush=True)
        return
    if args.benchmark_cell:
        for mode in MODES:
            record = search_cell([mode, 20., 30.])
            print(json.dumps({key: value for key, value in record.items() if key != "rows"},
                             ensure_ascii=False, indent=2), flush=True)
        return
    ROOT.mkdir(parents=True, exist_ok=True)
    config = campaign_config()
    config_path = ROOT / "config.json"
    if config_path.exists() and json.loads(config_path.read_text(encoding="utf-8")) != config:
        raise RuntimeError("Existing campaign code/inputs/environment differ; preserve frozen campaign")
    write_json(config_path, config)
    checkpoint = ROOT / "cells.jsonl"
    records = read_checkpoint(checkpoint)
    tasks = [[mode, t, pc] for mode in MODES for t in TEMPERATURES for pc in COUPLING_POWERS]
    todo = [task for task in tasks if cell_key(task) not in records]
    print(f"Actual cells: {len(tasks)}, completed: {len(records)}, pending: {len(todo)}; workers={args.workers}", flush=True)
    started = time.perf_counter()
    if todo:
        with checkpoint.open("a", encoding="utf-8") as stream, ProcessPoolExecutor(max_workers=args.workers) as pool:
            futures = {pool.submit(search_cell, task): task for task in todo}
            for index, future in enumerate(as_completed(futures), 1):
                record = future.result()
                records[record["key"]] = record
                stream.write(json.dumps(record, ensure_ascii=False, allow_nan=False) + "\n")
                stream.flush()
                if index == 1 or index % 5 == 0 or index == len(todo):
                    states = sum(r["control_evaluations"] for r in records.values())
                    print(f"{len(records)}/{len(tasks)} cells, {states} control states, "
                          f"{time.perf_counter()-started:.1f}s elapsed; last={record['task']}", flush=True)
    export_results(config, records, tasks)


if __name__ == "__main__":
    main()
