"""Source-pinned old/current/explicit-historical mean-gain comparison; no fit."""
import hashlib
import json
import sys
import types
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT))

from gabes.core import blas_single_thread
from gabes.schemes import fwm


def loaded_snapshot(filename, module_name, package):
    path = HERE / filename
    source = path.read_bytes()
    snapshot = types.ModuleType(module_name)
    snapshot.__package__ = package
    snapshot.__file__ = str(path)
    sys.modules[module_name] = snapshot
    exec(compile(source, str(path), "exec"), snapshot.__dict__)
    return snapshot, hashlib.sha256(source).hexdigest()


def main():
    old, old_sha = loaded_snapshot(
        "fwm_prepatch.py.txt", "gabes.schemes._mean_gain_before", "gabes.schemes")
    old_obs, old_obs_sha = loaded_snapshot(
        "observables_prepatch.py.txt", "gabes._mean_gain_obs_before", "gabes")
    assert old_sha == "e3d422163426445a3eb4e82b96a8056d7fd88d5b82a88ba696ba51bffd3b2e04"
    assert old_obs_sha == "0b708458aefb55abbb76ab4d7e62e8288e8aa322c40161e0f1b733067612c8e8"
    old.observables = old_obs
    center = fwm.branch_center_GHz(1.0, -1)
    optical_eta = 20.3063457 / 26.0
    kwargs = dict(
        T=386.15, P_pump=.2948, P_probe=3.06e-6 / optical_eta,
        line_strength=.74, detection_efficiency=.905 * optical_eta,
        coarse_points=401, fine_points=0, scan_min=center - .060,
        scan_max=center + .060, velocity_step=1., velocity_cutoff=4.,
        branch=-1, phase_detail=fwm.PHASE_ULTRA,
        model_fidelity=fwm.FIDELITY_BALANCED, response_method=fwm.RESPONSE_POLE,
        gain_closure_enabled=False)
    with blas_single_thread():
        before = old.compute_spectrum(1.0, **kwargs)
        historical = fwm.compute_spectrum(1.0, **kwargs, legacy_mean_field_fixture=True)
        current = fwm.compute_spectrum(1.0, **kwargs)
    parity = {key: bool(np.array_equal(
                  value, historical[key], equal_nan=value.dtype.kind in "fc"))
              for key, value in before.items() if isinstance(value, np.ndarray)}
    assert all(parity.values()), parity
    lab_tpd = -(current["probe_axis_GHz"] - center) * 1e3
    summaries = {}
    for name, spectrum in (("prepatch", before), ("current", current)):
        peak = int(np.argmax(spectrum["G_s"]))
        summaries[name] = {
            "peak_probe_gain": float(spectrum["G_s"][peak]),
            "peak_lab_tpd_MHz": float(lab_tpd[peak]),
            "peak_small_signal_gain": float(spectrum["G_s_smallsignal"][peak]),
            "absorbing_displayed_probe_points": int(np.sum(spectrum["G_s"] < 1.0)),
            "invalid_passive_points": int(np.sum(~spectrum["passive_completion_valid"])),
            "physical_squeezing_dB": spectrum["physical_squeezing_dB"],
        }
    changes = {}
    for key in ("G_s", "G_c", "G_s_smallsignal", "G_c_smallsignal",
                "omega_probe_rad_s", "delta_k_z_vacuum"):
        a, b = np.asarray(before[key]), np.asarray(current[key])
        changes[key] = {"max_abs": float(np.max(np.abs(b-a))),
                        "max_relative": float(np.max(np.abs(b-a) / np.maximum(np.abs(a), 1e-300)))}
    np.savez_compressed(HERE / "current_production_spectrum.npz", lab_tpd_MHz=lab_tpd,
                        **{k: v for k, v in current.items() if isinstance(v, np.ndarray)})
    result = {
        "record": "FWM-MEAN-GAIN-2026-10-15-PRODUCTION",
        "origin": "source_evidence", "no_fitting_performed": True,
        "prepatch_fwm_sha256": old_sha, "prepatch_observables_sha256": old_obs_sha,
        "current_fwm_sha256": hashlib.sha256(Path(fwm.__file__).read_bytes()).hexdigest(),
        "points": 401, "new_meanfield_scans": 3, "new_thermal_transport_solves": 0,
        "session_prior_comparison_attempt_scans": 3,
        "prior_attempt": "array comparator rejected Unicode status arrays after solves; dtype handling corrected, no physics change",
        "input_seed_W": kwargs["P_probe"],
        "input_seed_status": "AI handoff conditional vapor-entry assumption; measured pump-off output does not identify it",
        "detection_efficiency": kwargs["detection_efficiency"],
        "detection_efficiency_status": "AI handoff value; atomic/downstream split and RF calibration unresolved",
        "human_confirmed": "blue+1GHz from85Rb F2-Fprime3; seed below pump; pure85 with possible~0.1%87",
        "full_spectrum_summary": summaries,
        "historical_fixture_arrays_bitwise_equal": parity,
        "all_historical_fixture_arrays_bitwise_equal": all(parity.values()),
        "current_vs_prepatch": changes,
        "current_gain_reference": current["gain_reference"],
        "experimental_peak_from_ai_handoff": {
            "gain": 870.0 / 35.0, "lab_tpd_MHz": -2.2,
            "raw_data_independently_acquired_here": False},
        "minimum_same_plane_peak_gain_ratio_A0_1": summaries["current"]["peak_probe_gain"] / (870.0 / 35.0),
        "lab_peak_location_gap_MHz": summaries["current"]["peak_lab_tpd_MHz"] - (-2.2),
        "scope": "Conditional current mean-field diagnostic; static canonical L1 validity not microscopic atomic noise or experimental agreement.",
    }
    (HERE / "production_comparison.json").write_text(
        json.dumps(result, ensure_ascii=False, allow_nan=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"summaries": summaries, "all_historical_arrays_equal": all(parity.values()),
                      "array_count": len(parity), "minimum_gain_ratio": result["minimum_same_plane_peak_gain_ratio_A0_1"]}))


if __name__ == "__main__":
    main()
