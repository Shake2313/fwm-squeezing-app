"""Hypothetical, immutable single-Cmix sensitivity; no production correction edit.

Run from repository root: python analysis/tpd_gain_diagnostic_20261001/mixing_sensitivity.py
Only mixing_sensitivity.json and caches in this diagnostic directory are written.
"""
from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import sys

OUT = Path(__file__).resolve().parent
ROOT = OUT.parent.parent
sys.path.insert(0, str(ROOT))
sys.dont_write_bytecode = True
os.environ["PYTHONDONTWRITEBYTECODE"] = "1"
os.environ["NUMBA_CACHE_DIR"] = str(OUT / "numba_cache")
os.environ["MPLCONFIGDIR"] = str(OUT / "mpl_cache")
os.environ.setdefault("GABES_DISABLE_NUMBA", "1")

import numpy as np
from scipy.optimize import minimize_scalar
from gabes import doppler, fwm_gain_closure, hyperfine, observables
from gabes.schemes import fwm
import run_diagnostic as d

HYPOTHETICAL_C_FIXED = 0.8790717343
SOURCE_FILES = (
    "gabes/schemes/fwm.py", "gabes/fwm_gain_closure.py", "gabes/constants.py",
    "gabes/doppler.py", "gabes/observables.py", "gabes/hyperfine.py",
)


def source_hashes():
    return {name: hashlib.sha256((ROOT / name).read_bytes()).hexdigest()
            for name in SOURCE_FILES}


def predictor(params, eom_GHz):
    """Freeze atomic responses, vary only a diagnostic off-diagonal multiplier.

    Spreading the Maxwell weights onto the same Delta_eff nodes reproduces
    compute_spectrum's linear velocity interpolation exactly. No function,
    production constant, default, or source file is replaced or changed.
    """
    p = dict(params)
    eom = np.asarray(eom_GHz, dtype=float)
    T = p["temp_c"] + 273.15
    D = p["opd"]
    Delta = 2 * np.pi * D * 1e9
    probe = D - eom
    delta = fwm.two_photon_detuning_from_probe_scan(probe, D, -1)
    P_pump, P_seed = p["pump_mw"] * 1e-3, p["probe_uw"] * 1e-6
    L = p["cell_mm"] * 1e-3
    wp, ws = p["pump_waist_um"] * 1e-6, p["probe_waist_um"] * 1e-6
    angle = p["seeded_angle_deg"]
    density = hyperfine.number_density(T)
    atom = fwm.collisional_atom(
        T, density, transit_rate=2 * np.pi * p["transit_rate_khz"] * 1e3)
    Op, Os = fwm.rabi_freq(P_pump, wp), fwm.rabi_freq(P_seed, ws)
    v, weights = doppler.velocity_grid(T, dv=1.0, cutoff_sigma=4.0)
    nodes = doppler.build_Delta_eff_axis(Delta, Delta, v)
    lo, fraction = doppler.interpolation_weights(nodes, Delta, v)
    node_weights = np.zeros(nodes.size)
    np.add.at(node_weights, lo, weights * (1 - fraction))
    np.add.at(node_weights, lo + 1, weights * fraction)
    pole_rows = fwm.chi_pole_rows(
        Op, Op, Os, Os, delta, -1, atom=atom, n_f=p["floquet_order"],
        guard_delta_eff=Delta)
    if pole_rows is None:
        raise RuntimeError("Diagnostic pole guard failed; do not report this fit.")
    chi = pole_rows.mean(nodes, node_weights)
    _, k_probe, k_conj = fwm.seeded_option_a_wavenumbers(D, probe)
    delta_k = fwm.seeded_phase_mismatch_z(D, probe, angle_deg=angle)
    coupling = (p["line_strength"] * p["mode_overlap_penalty"]
                * p["polarization_penalty"] * p["zeeman_participation_penalty"]
                * fwm.physical_coupling_norm(-1))
    nseg = fwm.ULTRA_PROPAGATION_SEGMENTS
    segment_profile = np.ones(nseg)
    spatial = fwm._gaussian_overlap_profile(nseg, L, wp, ws, angle)

    def gain(hypothetical_C):
        ss, cs, sc, cc = fwm_gain_closure.apply(chi, float(hypothetical_C))
        gs, gc, _, _ = fwm._ultra_segmented_gain(
            ss, sc, cs, cc, k_probe, k_conj, L, density, coupling, delta_k,
            segment_profile, spatial, P_pump, P_seed)
        return observables.pump_depletion_saturation(gs, gc, P_pump, P_seed)

    return gain


def parity(params, eom, gain):
    """Assert helper parity against the unmodified public production pipeline."""
    eom = np.asarray(eom, dtype=float)
    if eom.size == 1:
        # The existing full-scan Floquet gate requires at least two points.
        # Include the requested Gold point and one 0.1-MHz neighbor for parity.
        eom = np.array([eom[0], eom[0] - 0.0001])
        gain = predictor(params, eom)
    out = {}
    for enabled, C in ((True, fwm_gain_closure.EFFECTIVE_PARTICIPATION),
                       (False, 1.0)):
        kwargs = d.kwargs(params)
        kwargs["gain_closure_enabled"] = enabled
        raw = fwm.compute_spectrum(
            params["opd"], scan_min=params["opd"] - float(np.max(eom)),
            scan_max=params["opd"] - float(np.min(eom)),
            coarse_points=len(eom), **kwargs)
        gs, gc = gain(C)
        np.testing.assert_allclose(gs, raw["G_s"], rtol=1e-9, atol=1e-10)
        np.testing.assert_allclose(gc, raw["G_c"], rtol=1e-9, atol=1e-10)
        out["current" if enabled else "off"] = {
            "C": C, "eom_GHz": eom, "G_s": gs, "G_c": gc,
            "max_relative_G_s": float(np.max(np.abs(gs / raw["G_s"] - 1))),
            "max_relative_G_c": float(np.max(np.abs(gc / raw["G_c"] - 1))),
            "floquet_status": raw["floquet_convergence"]["status"],
        }
    return out


def fit_single_C(gain, upper):
    objective = lambda C: float(np.mean((gain(C)[0] - d.MEAS) ** 2))
    grid = np.linspace(0.0, upper, 101)
    scores = np.array([objective(C) for C in grid])
    j = int(np.argmin(scores))
    lo, hi = grid[max(j - 1, 0)], grid[min(j + 1, grid.size - 1)]
    fit = minimize_scalar(objective, bounds=(lo, hi), method="bounded",
                          options={"xatol": 1e-12})
    candidates = [float(fit.x), float(grid[j]), 0.0, upper]
    C = min(candidates, key=objective)
    gs, gc = gain(C)
    return dict(C=C, search_bounds=[0.0, upper], mse=objective(C),
                fits=d.fits(gs), G_s=gs, G_c=gc,
                endpoint_ratio=float(gs[-1] / gs[0]),
                attenuation_bound_exceeded=bool(C > 1.0),
                optimizer_success=bool(fit.success))


def main():
    before = source_hashes()
    sensitivity = []
    dense_eom = np.linspace(3.070, 2.970, 201)
    for temp in (113.0, 118.0, 123.0):
        p = dict(d.PARAMS, temp_c=temp)
        gain = predictor(p, d.EOM)
        reconstruction = parity(p, d.EOM, gain)
        unconstrained = fit_single_C(gain, 2.0)
        attenuation_constrained = fit_single_C(gain, 1.0)
        peak_gain = predictor(p, dense_eom)(unconstrained["C"])[0]
        j = int(np.argmax(peak_gain))
        unconstrained["peak"] = dict(
            eom_GHz=float(dense_eom[j]), G_s=float(peak_gain[j]),
            eom_scan_bounds_GHz=[2.970, 3.070], grid_spacing_MHz=0.5,
            endpoint_maximum=bool(j in (0, dense_eom.size - 1)))
        sensitivity.append(dict(temperature_C=temp, params=p,
                                parity=reconstruction,
                                hypothetical_single_C=unconstrained,
                                hypothetical_attenuation_only_C=attenuation_constrained))

    gold = dict(d.DEFAULTS, temp_c=121.0, pump_mw=600.0, opd=0.9, probe_uw=8.0)
    gold_effect = []
    # -8 MHz is the code's frozen calibration point. +8 is a sign audit only:
    # the primary paper reports -8 MHz without an explicit algebraic mapping.
    for tpd in (-8.0, 8.0):
        eom = np.array([fwm.constants.NU_HF / 1e9 - tpd * 1e-3])
        gain = predictor(gold, eom)
        reconstruction = parity(gold, eom, gain)
        current_gs, current_gc = gain(fwm_gain_closure.EFFECTIVE_PARTICIPATION)
        alternative_gs, alternative_gc = gain(HYPOTHETICAL_C_FIXED)
        gold_effect.append(dict(
            GABES_tpd_MHz=tpd, eom_GHz=eom,
            role=("frozen production calibration point" if tpd == -8.0
                  else "opposite sign diagnostic; not an adopted reference"),
            parity=reconstruction,
            current_C=fwm_gain_closure.EFFECTIVE_PARTICIPATION,
            current_G_s=current_gs, current_G_c=current_gc,
            hypothetical_C=HYPOTHETICAL_C_FIXED,
            hypothetical_G_s=alternative_gs, hypothetical_G_c=alternative_gc,
            relative_G_s_change_pct=100 * (alternative_gs / current_gs - 1)))

    after = source_hashes()
    assert before == after, "Production sources changed during the diagnostic."
    result = dict(
        kind="hypothetical parameter sensitivity only; no correction implemented or applied",
        method="same atomic Floquet response and Maxwell measure; diagnostic Cmix before propagation",
        residual_convention="measurement minus prediction, as in run_diagnostic.py",
        measured_eom_GHz=d.EOM, measured_gain=d.MEAS, actual_params=d.PARAMS,
        measured_endpoint_ratio=float(d.MEAS[-1] / d.MEAS[0]),
        current_production_C=fwm_gain_closure.EFFECTIVE_PARTICIPATION,
        hypothetical_C_fixed_118C=HYPOTHETICAL_C_FIXED,
        temperature_sensitivity=sensitivity,
        frozen_gold_effect=gold_effect,
        source_preservation=dict(before_sha256=before, after_sha256=after,
                                 all_listed_sources_unchanged=True),
        limitations=[
            "One common multiplier per whole ten-point curve; no pointwise fit.",
            "C above one is a mathematical diagnostic, not attenuation participation.",
            "Temperature bounds are sensitivity scenarios, not confidence intervals.",
            "No measurement uncertainty or statistical significance is inferred.",
            "A trial C is used only in this helper; production constants remain frozen.",
            "Primary-paper TPD sign must be independently resolved before recalibration.",
        ])
    (OUT / "mixing_sensitivity.json").write_text(
        json.dumps(d.plain(result), ensure_ascii=False, indent=2), encoding="utf-8")
    for item in sensitivity:
        fit = item["hypothetical_single_C"]
        print("Hypothetical Cmix", item["temperature_C"], fit["C"],
              "RMSE", fit["fits"]["identity"]["rmse"], flush=True)
    print("Frozen Gold effect", json.dumps(d.plain(gold_effect)), flush=True)


if __name__ == "__main__":
    main()
