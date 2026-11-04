"""Protocol-separated Gold checks and a small current-model resolution pilot.

No absolute-transmission calibration is inferred from arbitrary signal units.
No RF sensitivity or temperature optimum is fitted to the published values.
"""
from __future__ import annotations

import os

for variable in ("OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "OMP_NUM_THREADS"):
    os.environ[variable] = "1"

import csv
import hashlib
import json
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parents[1]
sys.path.insert(0, str(REPO))

from gabes import affine_doppler, constants, core, rydberg_electrometry
from gabes.schemes.rydberg import MHZ, RydbergEITScheme, SUPERHET_CLASS_CHUNK

POWERS = [.5, 1., 2., 4., 6., 8., 10.]
TRAIN = [.5, 2., 6., 10.]
HELDOUT = [1., 4., 8.]
AXES = ["Probe detuning", "Coupling detuning"]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_json(name, value):
    (ROOT / name).write_text(json.dumps(value, ensure_ascii=False, indent=2,
        allow_nan=False) + "\n", encoding="utf-8")


def csv_rows(name, rows):
    with (ROOT / name).open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def source_hashes():
    paths = sorted((REPO / "gabes").glob("*.py")) + [
        REPO / "gabes/schemes/rydberg.py", Path(__file__),
        ROOT / "sources/gold_observations.json", ROOT / "sources/fig2b_points.csv"]
    return {str(p.relative_to(REPO)).replace("\\", "/"): sha(p) for p in paths}


class BenchmarkScheme(RydbergEITScheme):
    """Fixed physical scan and one velocity quadrature per control state."""

    def __init__(self):
        super().__init__()
        self.velocity = self.weights = None

    def _scan(self, params):
        return (np.linspace(-12., 12., 1601) + self._resonance_scan_mhz(params)) * MHZ

    def _superheterodyne_atomic_phasors(self, raw, params, x, transmission, static_indices=None):
        if params["doppler"] != "on":
            return super()._superheterodyne_atomic_phasors(raw, params, x, transmission, static_indices)
        indices = np.asarray(static_indices, dtype=int)
        rate = sum(raw[k] for k in ("rydberg_dephasing_mhz", "temperature_dephasing_mhz",
                                    "density_dephasing_mhz", "transit_mhz")) * MHZ
        atom = self._atom(rate, raw["rf_dephasing_mhz"] * MHZ)
        if self.velocity is None:
            affine = affine_doppler.AffineScan(atom,
                lambda delta: self._hamiltonian_at_scan_detuning(raw, params, delta),
                (1, 0), raw["probe_rabi_mhz"] * MHZ)
            self.velocity, self.weights = affine.quadrature_classes(
                raw["scan"], raw["T"], constants.MASS_85RB, raw["k_vec"])
        kv = raw["k_vec"] * self.velocity
        operator = self._absorption_quadrature_operator()
        omega_if = params["if_khz"] * 1e3 * 2 * np.pi
        values = np.empty(indices.size, dtype=complex)
        for slot, index in enumerate(indices):
            base = core.build_liouvillian(self._hamiltonian_at_scan_detuning(
                raw, params, raw["scan"][index]), atom)
            classes = np.empty(kv.size, dtype=complex)
            for first in range(0, kv.size, SUPERHET_CLASS_CHUNK):
                part = slice(first, first + SUPERHET_CLASS_CHUNK)
                matrix = base[None] + kv[part, None, None] * atom.S_v[None]
                response = rydberg_electrometry.weak_signal_response_from_liouvillian(
                    matrix, 4, omega_if, signal_transition=(2, 3))
                classes[part] = response.real_observable_phasor_per_angular_rabi(operator)
            values[slot] = np.dot(self.weights, classes)
        return indices, values, "Explicit protocol points; fixed full-scan velocity quadrature"


def parameters(scheme, mode, axis, power=6., **extra):
    result = scheme.recommended_defaults(scheme.defaults())["AT electrometry"]
    result.update(doppler=mode, scan_axis=axis, probe_power_uw=power,
        coupling_power_mw=30., coupling_rabi_mhz=3., probe_rabi_ref_mhz=2.,
        cell_mm=50., beam_diameter_mm=.15, temp_c=20.,
        temperature_model="Linked", probe_detuning_mhz=0., coupling_detuning_mhz=0.,
        lo_rabi_mhz=3.7, mw_frequency_ghz=37.07104, mw_detuning_mhz=0., if_khz=40.)
    result.update(extra)
    return result


def rf_protocol(mode, axis, power=6., **extra):
    scheme = BenchmarkScheme()
    params = parameters(scheme, mode, axis, power, **extra)
    raw = scheme.compute(params)
    x, transmission, _ = scheme._transmission(raw["chi_bar"], raw, params)
    center = scheme._resonance_scan_mhz(params)
    center_index = int(np.argmin(abs(x-center)))
    slopes = scheme._static_slope_indices(x, transmission)
    indices = np.union1d(slopes, [center_index]).astype(int)
    atomic = scheme._superheterodyne_atomic_phasors(raw, params, x, transmission,
                                                   static_indices=indices)
    result = []
    for policy, selected in (("maximum_static_slope", slopes),
                             ("bare_two_photon_center", np.array([center_index]))):
        mask = np.isin(atomic[0], selected)
        reused = atomic[0][mask], atomic[1][mask], policy
        readout = scheme._superheterodyne_readout(raw, params, x, transmission, atomic_response=reused)
        if readout is None:
            raise RuntimeError("No finite response at declared protocol points")
        result.append(dict(doppler=mode, scan_axis=axis, probe_power_uw=power,
            temperature_c=params["temp_c"], coupling_power_mw=params["coupling_power_mw"],
            fixed_probe_detuning_mhz=params["probe_detuning_mhz"],
            fixed_coupling_detuning_mhz=params["coupling_detuning_mhz"],
            lo_rabi_mhz=3.7, lo_carrier_ghz=37.07104, rf_detuning_mhz=0., if_khz=40., detuning_policy=policy,
            chosen_scan_detuning_mhz=float(readout["optimum_scan_detuning_mhz"]),
            chosen_probe_detuning_mhz=float(readout["optimum_probe_detuning_mhz"]),
            chosen_coupling_detuning_mhz=float(readout["optimum_coupling_detuning_mhz"]),
            sensitivity_nv_cm_sqrt_hz=float(readout["sensitivity"].total_field_asd_nv_cm_per_sqrt_hz),
            transmission=float(readout["optimum_transmission"]),
            responsivity_a_per_v_m=float(readout["current_responsivity_a_per_v_m"]),
            detected_power_w=float(readout["signal_detector_power_w"]),
            velocity_classes=0 if scheme.velocity is None else len(scheme.velocity),
            candidate_count=int(mask.sum()),
            scope="Fixed LO and declared optical-point policy; no claim of exact unpublished Gold detunings"))
    return result


def main():
    started = time.perf_counter()
    before = source_hashes()
    references = list(csv.DictReader((ROOT / "sources/fig2b_points.csv").open(encoding="utf-8")))
    observed = {(float(row["probe_power_uw"]), row["condition"], row["quantity"]): float(row["value"])
                for row in references}
    rows, arrays = [], {}
    for mode in ("off", "on"):
        for axis in AXES:
            for power in POWERS:
                scheme = BenchmarkScheme()
                params = parameters(scheme, mode, axis, power, view="EIT", lo_rabi_mhz=0.)
                raw = scheme.compute(params)
                for condition, key in (("compensated", "chi_bar"), ("uncompensated", "chi_bar_uncomp")):
                    x, transmission, _ = scheme._transmission(raw[key], raw, params)
                    width, contrast = scheme._eit_features(x, transmission,
                        center_mhz=scheme._resonance_scan_mhz(params))
                    rows.append(dict(doppler=mode, scan_axis=axis, probe_power_uw=power,
                        condition=condition, lo_rabi_mhz=0., eit_width_mhz=float(width),
                        model_contrast=float(contrast), model_peak_transmission=float(transmission[len(x)//2]),
                        model_relative_signal=float(power*contrast),
                        source_width_mhz=observed[power, condition, "eit_fwhm"],
                        source_peak_arb=observed[power, condition, "eit_peak"],
                        partition="training" if power in TRAIN else "heldout"))
                    prefix=f"{mode}_{AXES.index(axis)}_{power:g}_{condition}"
                    arrays[prefix+"_x"] = x
                    arrays[prefix+"_transmission"] = transmission
            print(f"Static spectroscopy done: Doppler {mode}, {axis}", flush=True)
    comparisons = []
    for mode in ("off", "on"):
        for axis in AXES:
            group = [row for row in rows if row["doppler"] == mode and row["scan_axis"] == axis]
            train = [row for row in group if row["partition"] == "training"]
            model = np.array([row["model_relative_signal"] for row in train])
            target = np.array([row["source_peak_arb"] for row in train])
            gain = float(np.dot(model, target)/np.dot(model, model))
            for row in group:
                row["signal_gain_fitted_on_training"] = gain
                row["predicted_signal_arb"] = gain*row["model_relative_signal"]
                row["width_error_mhz"] = row["eit_width_mhz"]-row["source_width_mhz"]
                row["signal_error_arb"] = row["predicted_signal_arb"]-row["source_peak_arb"]
            for split in ("training", "heldout"):
                part = [row for row in group if row["partition"] == split]
                comparisons.append(dict(doppler=mode, scan_axis=axis, partition=split,
                    count=len(part), shared_detector_gain=gain,
                    width_rmse_mhz=float(np.sqrt(np.mean([row["width_error_mhz"]**2 for row in part]))),
                    maximum_absolute_width_error_mhz=max(abs(row["width_error_mhz"]) for row in part),
                    signal_rmse_arb=float(np.sqrt(np.mean([row["signal_error_arb"]**2 for row in part]))),
                    comparison_scope="One detector gain shared across powers and compensation conditions; absolute transmission/physical contrast not identified"))
    csv_rows("spectroscopy_comparison.csv", rows)
    np.savez_compressed(ROOT/"matched_spectra.npz", **arrays)

    rf_rows = []
    for mode in ("off", "on"):
        for axis in AXES:
            for power in POWERS:
                rf_rows.extend(rf_protocol(mode, axis, power))
            print(f"Fixed-LO RF checks done: Doppler {mode}, {axis}", flush=True)
    csv_rows("fixed_lo_rf.csv", rf_rows)

    pilot = []
    for temperature in (15., 20., 35.):
        for coupling in (20., 30., 40.):
            pilot.extend(rf_protocol("on", "Probe detuning", 6., temp_c=temperature,
                                     coupling_power_mw=coupling))
    csv_rows("current_model_sparse_pilot.csv", pilot)
    pilot_errors = []
    for policy in ("maximum_static_slope", "bare_two_photon_center"):
        values = {(row["temperature_c"], row["coupling_power_mw"]): row["sensitivity_nv_cm_sqrt_hz"]
                  for row in pilot if row["detuning_policy"] == policy}
        errors = []
        for (temperature, coupling), target in values.items():
            if temperature == 20 or coupling == 30:
                continue
            predicted = values[temperature,30.]*values[20.,coupling]/values[20.,30.]
            errors.append(dict(temperature_c=temperature, coupling_power_mw=coupling,
                predicted=predicted, direct=target, relative_error=predicted/target-1))
        maximum = max(abs(row["relative_error"]) for row in errors)
        pilot_errors.append(dict(detuning_policy=policy, training_coordinates=5,
            guard_coordinates=4, maximum_relative_error=maximum,
            candidate_sparse_rule_eligible_at_one_percent=maximum<=.01,
            evidence_scope="New CP model, fixed Pp6/LO3.7, 3x3 pilot only; not a certificate for a reoptimized 21x9 map",
            errors=errors))
    detuning_rows = []
    for delta in (-2., 0., 2.):
        detuning_rows.extend(rf_protocol("on", "Probe detuning", 6., coupling_detuning_mhz=delta))
    csv_rows("independent_detuning_pilot.csv", detuning_rows)
    after = source_hashes()
    if before != after:
        raise RuntimeError("Source changed during benchmark; preserve files but do not accept this receipt")
    reference = json.loads((ROOT/"sources/gold_observations.json").read_text(encoding="utf-8"))
    result = dict(origin="⟦AI⟧ Gold protocol-separated spectroscopy and fixed-LO analysis",
        model="CP nested optical/RF phase channels; explicit physical optical detunings",
        protocol_status="PARTIAL_MATCH: LO/powers/geometry declared; source optical scan identity and absolute detunings not reported",
        absolute_transmission_validation="NOT_IDENTIFIABLE_FROM_ARBITRARY_SIGNAL_UNITS",
        physical_contrast_validation="NOT_IDENTIFIABLE_WITHOUT_GAIN_BASELINE",
        amplitude_calibration_scope="One positive detector gain fitted only on training powers, shared with heldout powers and both B-field conditions",
        training_powers_uw=TRAIN, heldout_powers_uw=HELDOUT,
        comparisons=comparisons, source_protocols=reference["protocols"],
        source_consistency_issues=reference["source_issues"],
        current_model_predictability=pilot_errors,
        optimization_status="EXPERIMENTAL_OPTIMUM_NOT_VALIDATED: physical protocol and absolute readout calibration remain unresolved",
        reference_lo_carrier_ghz=37.07104,
        rf_detuning_scope="Nominal delta_RF=0 is an analysis assumption; source carrier is reported but independent atomic RF zero calibration is not supplied",
        next_grid_policy="Use sparse T/Pc cross plus guards only where current-model error criterion passes; retain Pp/LO interaction and signed optical detuning candidates",
        no_sensitivity_fit=True, no_temperature_optimum_fit=True,
        elapsed_seconds=time.perf_counter()-started, source_sha256=before)
    write_json("benchmark_results.json", result)
    print(json.dumps(dict(comparisons=comparisons, current_model_predictability=pilot_errors,
                         elapsed_seconds=result["elapsed_seconds"]), indent=2), flush=True)


if __name__ == "__main__":
    main()
