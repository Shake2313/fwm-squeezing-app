"""Small conditional Gaussian-beam spectroscopy trial, no production edits."""
from __future__ import annotations

import csv
import hashlib
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parents[2]
sys.path.insert(0, str(ROOT))
from linewidth_profiles import Candidate
from gabes.schemes.rydberg import MHZ

TRAIN = [.5, 2., 6., 10.]
HOLDOUT = [1., 4., 8.]
OPACITIES = [.001, .003, .01, .03, .1, .3, 1.]
ANCHORS = [3., 4., 5., 6.]
GROUND_RATES = [0., .05, .1, .15, .25]


def evaluate(gamma, anchor, powers, nodes=12, scan_step=1/30):
    candidate = Candidate("nested_cp", anchor, gamma)
    candidate._scan = lambda params: np.linspace(-10., 10., round(20./scan_step)+1)*MHZ
    q, weights = np.polynomial.laguerre.laggauss(nodes)
    params = candidate.recommended_defaults(candidate.defaults())["EIT"]
    rows = []
    for power in powers:
        curves = {opacity: [] for opacity in OPACITIES}
        for u, weight in zip(q, weights):
            local = dict(params, doppler="on", probe_power_uw=power*np.exp(-u),
                coupling_power_mw=30.*np.exp(-u), coupling_rabi_mhz=3.,
                rydberg_dephasing_mhz=gamma, lo_rabi_mhz=0., view="EIT")
            raw = candidate.compute(local)
            for opacity in OPACITIES:
                x, curve, _ = candidate._transmission(raw["chi_bar"], raw,
                    dict(local, effective_opacity_participation=opacity))
                curves[opacity].append(weight*curve)
        for opacity in OPACITIES:
            transmission = np.sum(curves[opacity], axis=0)
            width, contrast = candidate._eit_features(x, transmission)
            rows.append(dict(total_ground_coherence_rate_mhz=gamma,
                probe_rabi_peak_anchor_mhz=anchor, coupling_rabi_peak_mhz=3.,
                opacity_participation=opacity, probe_power_uw=power,
                eit_fwhm_mhz=float(width), eit_contrast=float(contrast),
                relative_voltage_amplitude=power*float(contrast),
                transmission_at_zero=float(transmission[len(x)//2]),
                radial_nodes=nodes))
    return rows


def main():
    data_path = ROOT.parent / "sources/fig2b_points.csv"
    sources = [Path(__file__), ROOT/'linewidth_profiles.py', data_path, REPO/'gabes/schemes/rydberg.py']
    before = {str(path.relative_to(REPO)).replace('\\','/'):hashlib.sha256(path.read_bytes()).hexdigest() for path in sources}
    with data_path.open(encoding="utf-8-sig", newline="") as stream:
        source = list(csv.DictReader(stream))
    observations = {(float(row["probe_power_uw"]), row["quantity"]):float(row["value"])
                    for row in source if row["condition"] == "compensated"}

    train_rows = []
    for gamma in GROUND_RATES:
        for anchor in ANCHORS:
            train_rows.extend(evaluate(gamma, anchor, TRAIN))
            print("Trial", gamma, anchor, flush=True)
    scores = []
    for gamma in GROUND_RATES:
        for anchor in ANCHORS:
            for opacity in OPACITIES:
                rows = [row for row in train_rows if row["total_ground_coherence_rate_mhz"] == gamma
                        and row["probe_rabi_peak_anchor_mhz"] == anchor
                        and row["opacity_participation"] == opacity]
                model = np.array([row["relative_voltage_amplitude"] for row in rows])
                measured = np.array([observations[row["probe_power_uw"], "eit_peak"] for row in rows])
                gain = float(model@measured/(model@model))
                width_error = np.array([row["eit_fwhm_mhz"]-observations[row["probe_power_uw"],"eit_fwhm"] for row in rows])
                amp_error = gain*model-measured
                wr = float(np.sqrt(np.mean(width_error**2)))
                ar = float(np.sqrt(np.mean(amp_error**2)))
                scores.append(dict(total_ground_coherence_rate_mhz=gamma,
                    probe_rabi_peak_anchor_mhz=anchor, opacity_participation=opacity,
                    one_global_amplitude_gain=gain, train_width_rmse_mhz=wr,
                    train_amplitude_rmse_arb=ar, score=(wr/3.)**2+(ar/.735)**2))
    winner = min(scores, key=lambda row:row["score"])
    holdout = [row for row in evaluate(winner["total_ground_coherence_rate_mhz"],
        winner["probe_rabi_peak_anchor_mhz"], HOLDOUT)
        if row["opacity_participation"] == winner["opacity_participation"]]
    for row in holdout:
        row["measured_eit_fwhm_mhz"] = observations[row["probe_power_uw"], "eit_fwhm"]
        row["predicted_amplitude_arb"] = winner["one_global_amplitude_gain"] * row["relative_voltage_amplitude"]
        row["measured_amplitude_arb"] = observations[row["probe_power_uw"], "eit_peak"]
    converged = [row for row in evaluate(winner["total_ground_coherence_rate_mhz"],
        winner["probe_rabi_peak_anchor_mhz"], TRAIN+HOLDOUT, nodes=24, scan_step=.01)
        if row["opacity_participation"] == winner["opacity_participation"]]
    coarse_by_power = {row["probe_power_uw"]:row for row in train_rows+holdout
        if row["total_ground_coherence_rate_mhz"] == winner["total_ground_coherence_rate_mhz"]
        and row["probe_rabi_peak_anchor_mhz"] == winner["probe_rabi_peak_anchor_mhz"]
        and row["opacity_participation"] == winner["opacity_participation"]}
    for row in converged:
        coarse = coarse_by_power[row["probe_power_uw"]]
        row["width_difference_from_nodes12_step_1over30_mhz"] = row["eit_fwhm_mhz"]-coarse["eit_fwhm_mhz"]
        row["relative_amplitude_change_from_nodes12_step_1over30_mhz"] = row["relative_voltage_amplitude"]/coarse["relative_voltage_amplitude"]-1
    winner["touches_tested_parameter_bound"] = (
        winner["total_ground_coherence_rate_mhz"] in (GROUND_RATES[0], GROUND_RATES[-1])
        or winner["probe_rabi_peak_anchor_mhz"] in (ANCHORS[0], ANCHORS[-1])
        or winner["opacity_participation"] in (OPACITIES[0], OPACITIES[-1]))
    after = {str(path.relative_to(REPO)).replace('\\','/'):hashlib.sha256(path.read_bytes()).hexdigest() for path in sources}
    result = dict(origin="⟦AI⟧ Bounded Gaussian-beam hypothesis test, not a calibrated production model",
        physics_hypothesis="Reported coupling Rabi3MHz is interpreted as Gaussian beam peak; paper does not establish peak-vs-mean. Equal probe/coupling radii; local Omega exp(-r^2/w^2).",
        quadrature="u=2r^2/w^2; detector integrates exp(-u)*T(u), T(u)=exp[-alpha(u)L]; Gauss-Laguerre12 then selected-hypothesis24. No average-chi-then-exponent shortcut.",
        tested_parameters=dict(ground_rate_mhz=GROUND_RATES, probe_peak_anchor_mhz=ANCHORS, opacity_participation=OPACITIES),
        training_powers_uw=TRAIN, heldout_powers_uw=HOLDOUT,
        amplitude_gain_scope="One gain fitted to all four training powers, applied unchanged to holdout. No absolute-transmission calibration inferred.",
        score_scope="Unweighted relative-scale engineering comparison: (widthRMSE/3MHz)^2+(amplitudeRMSE/0.735arb)^2. Not a statistical chi-square or confidence interval.",
        train_rows=train_rows, train_scores=scores, selected_training_hypothesis=winner,
        heldout_rows=holdout, selected_hypothesis_convergence_rows=converged,
        source_sha256=before, source_sha256_after=after,
        changed_sources_during_run=[name for name in before if before[name]!=after[name]])
    (ROOT/'gaussian_trial.json').write_text(json.dumps(result, ensure_ascii=False, indent=2, allow_nan=False)+"\n", encoding="utf-8")
    print(json.dumps(dict(selected=winner,holdout=holdout), indent=2))


if __name__ == "__main__":
    main()
