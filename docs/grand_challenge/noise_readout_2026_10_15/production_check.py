"""Reproduce the AI handoff's explicitly conditional production example."""
import hashlib
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
from gabes.core import blas_single_thread
from gabes.passive_readout import passive_l1_from_spectrum
from gabes.schemes import fwm

OUT = Path(__file__).resolve().parent


def main():
    optical_eta = 20.3063457 / 26.0
    eta = 0.905 * optical_eta
    center = fwm.branch_center_GHz(1.0, -1)
    with blas_single_thread():
        raw = fwm.compute_spectrum(
            1.0, T=386.15, P_pump=294.8e-3,
            P_probe=(3.06 / optical_eta) * 1e-6, line_strength=0.74,
            detection_efficiency=eta, coarse_points=401, fine_points=0,
            scan_min=center - 0.060, scan_max=center + 0.060,
            velocity_step=1.0, velocity_cutoff=4.0, branch=-1,
            phase_detail=fwm.PHASE_ULTRA, model_fidelity=fwm.FIDELITY_BALANCED,
            response_method=fwm.RESPONSE_POLE, gain_closure_enabled=False)
    coherent = passive_l1_from_spectrum(raw)
    noisy = passive_l1_from_spectrum(raw, seed_fano=4.05)
    lab = -(raw["probe_axis_GHz"] - center) * 1e3
    rows = []
    for target in (-36.0, -24.0, -10.8, 0.0, 20.1, 39.9):
        i = int(np.argmin(np.abs(lab - target)))
        rows.append({
            "lab_tpd_mhz": float(lab[i]),
            "D_canonical": float(coherent["passive_l1_D"][i]),
            "legacy_dB": float(raw["gain_referred_noise_dB"][i]),
            "formal_l1_F1_dB": float(coherent["passive_l1_formal_noise_dB"][i]),
            "formal_l1_F4_05_dB": float(noisy["passive_l1_formal_noise_dB"][i]),
            "completion_valid": bool(coherent["passive_completion_valid"][i]),
            "display_l1_F1_dB": (float(coherent["passive_l1_noise_dB"][i])
                                  if np.isfinite(coherent["passive_l1_noise_dB"][i]) else None),
        })
    invalid = ~coherent["passive_completion_valid"]
    peak = int(np.argmax(raw["G_s"]))
    np.savez_compressed(OUT / "production_spectrum.npz", lab_tpd_mhz=lab,
                        **{k: v for k, v in raw.items() if isinstance(v, np.ndarray)},
                        l1_f1_formal_dB=coherent["passive_l1_formal_noise_dB"],
                        l1_f4_05_formal_dB=noisy["passive_l1_formal_noise_dB"])
    baseline = ROOT / "tests/baseline_focused.npz"
    summary = {
        "record": "FWM-NOISE-2026-10-15-PRODUCTION",
        "origin": "source_evidence",
        "seed_fano_4_05_status": "handoff assumption; not independently measured here",
        "operating_point_human_confirmed": {
            "opd": "blue +1 GHz from85Rb F2-Fprime3", "branch": -1,
            "cell": "pure85Rb, possible~0.1%87Rb contamination"},
        "eta_handoff_value_unchanged": eta,
        "eta_value_source": "AI handoff; not independently re-confirmed in this chat",
        "eta_scope": "handoff total value used conditionally; downstream/cell split unresolved",
        "points": len(lab), "invalid_passive_points": int(np.sum(invalid)),
        "invalid_lab_interval_mhz": [float(np.min(lab[invalid])), float(np.max(lab[invalid]))],
        "minimum_completion_eigenvalue": float(np.min(coherent["passive_completion_min_eig"])),
        "peak_lab_tpd_mhz": float(lab[peak]), "peak_displayed_gain": float(raw["G_s"][peak]),
        "baseline_sha256": hashlib.sha256(baseline.read_bytes()).hexdigest(),
        "physical_squeezing_dB": raw["physical_squeezing_dB"],
        "rows": rows,
    }
    (OUT / "production_summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
