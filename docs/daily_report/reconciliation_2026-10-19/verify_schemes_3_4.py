"""⟦AI⟧ One-pass report checks; diagnostic model changes are never applied."""
import json
import sys
from pathlib import Path
from unittest.mock import patch

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from gabes import schemes, core, observables, zeeman
from gabes.schemes.rydberg import MHZ


def main():
    out = {}
    sc = schemes.get("rydberg_eit")
    p = dict(sc.defaults(), **sc.recommended_defaults(sc.defaults())["AT electrometry"])
    raw = sc.compute(p)
    out["opacity_sensitivity_nv_cm_sqrt_hz"] = {
        str(ls): sc._readout(dict(raw, ls=ls), p)["superhet"]["sensitivity"]
        .total_field_asd_nv_cm_per_sqrt_hz for ls in [.0005, .001, .002]
    }
    hot = dict(p, view="EIT", lo_rabi_mhz=0., temp_c=80.)
    hot_raw = sc.compute(hot)
    x, transmission, _ = sc._transmission(hot_raw["chi_bar"], hot_raw, hot)
    out["80_c"] = {
        "transmission_zero": float(np.interp(0, x, transmission)),
        "heroes": [m for m in sc.headless_observables(hot_raw, hot)["metrics"]
                   if m.get("tier") == "hero"],
    }
    weak = dict(p, lo_rabi_mhz=.1)
    out["weak_lo_heroes"] = [
        m for m in sc.headless_observables(sc.compute(weak), weak)["metrics"]
        if m.get("tier") == "hero"]
    alt = dict(p, if_khz=200., detector_quantum_efficiency=.2,
               detector_path_efficiency=.2, cell_mm=20.,
               rf_transition_dipole_ea0=700., rf_angular_factor=.5,
               detector_reference_power_ratio=.7,
               detector_electronic_noise_pa_sqrt_hz=5.,
               detector_rin_per_sqrt_hz=1e-6, detector_rin_correlation=.6,
               atom_participation_fraction=.3, beam_overlap_efficiency=.4)
    rkeys = {s.name for s in sc.param_schema() if s.recompute}
    out["extra_views"] = {}
    for view in sc.extra_views():
        full = view.compute(alt)
        filtered = view.compute({k: v for k, v in alt.items()
                                 if k in rkeys or k in view.param_keys})
        parity = all(
            np.array_equal(np.asarray(full[k]), np.asarray(filtered[k]), equal_nan=True)
            if np.asarray(full[k]).dtype.kind in "fc" else full[k] == filtered[k]
            for k in full)
        with patch.object(sc, "compute", wraps=sc.compute) as counter:
            view.compute(dict(alt, detector_electronic_noise_pa_sqrt_hz=6.))
        out["extra_views"][view.key] = {
            "filtered_full_parity": parity,
            "detector_change_atomic_compute_count": counter.call_count}

    gamma = sum(raw[k] for k in ["rydberg_dephasing_mhz", "temperature_dephasing_mhz",
                                "density_dephasing_mhz", "transit_mhz"]) * MHZ
    atom = sc._atom(gamma, raw["rf_dephasing_mhz"] * MHZ)
    a = core.build_liouvillian(sc._hamiltonian_at_probe_detuning(raw, p, 0.), atom)
    b = core.build_liouvillian(sc._hamiltonian_at_probe_detuning(raw, p, 1.), atom) - a
    old = np.asarray([core.build_liouvillian(
        sc._hamiltonian_at_probe_detuning(raw, p, d), atom) for d in raw["scan"]])
    new = a[None] + raw["scan"][:, None, None] * b[None]
    out["affine_assembly"] = {
        "array_equal": bool(np.array_equal(old, new)),
        "max_abs_error": float(np.max(abs(old - new)))}

    sc = schemes.get("magneto")
    regs = {k: dict(sc.defaults(), **v)
            for k, v in sc.recommended_defaults(sc.defaults()).items()}
    out["magneto_presets"] = {}
    for name, p in regs.items():
        raw = sc.compute(p)
        out["magneto_presets"][name] = {
            "feature": next((m["value"] for m in sc.headless_observables(raw, p)["metrics"]
                             if m["label"] == "Zero-field feature"), None),
            "fallback_rows": raw["doppler_fallback_rows"]}
    p = next(dict(v) for k, v in regs.items() if "NMOR" in k)
    out["nmor_slope_mrad_per_ut"] = {}
    for points in [121, 401, 801]:
        p["scan_points"] = points
        raw = sc.compute(p)
        xp = observables.chi_phys(raw["chi_p"], raw["N_eff"],
                                 dipole=raw["dipole"], line_strength=p["line_strength"])
        xm = observables.chi_phys(raw["chi_m"], raw["N_eff"],
                                 dipole=raw["dipole"], line_strength=p["line_strength"])
        rotation = .25 * raw["k_vec"] * p["cell_mm"] * 1e-3 * np.real(xp - xm)
        x = raw["b_ut"]
        i = int(np.argmin(abs(x)))
        step = x[i + 1] - x[i]
        out["nmor_slope_mrad_per_ut"][points] = {
            "three_point": float(np.gradient(rotation, x)[i] * 1e3),
            "five_point": float((rotation[i - 2] - 8 * rotation[i - 1]
                                 + 8 * rotation[i + 1] - rotation[i + 2])
                                / (12 * step) * 1e3)}
    p = dict(sc.defaults(), b_offset_ut=.25, scan_points=401)
    raw = sc.compute(p)
    alpha = raw["k_vec"] * np.imag(observables.chi_phys(
        raw["chi_probe"], raw["N_eff"], dipole=raw["dipole"],
        line_strength=p["line_strength"]))
    transmission = np.exp(-alpha * p["cell_mm"] * 1e-3)
    out["zero_reference"] = {
        "command_zero_transmission": float(np.interp(0, raw["b_ut"], transmission)),
        "physical_zero_transmission": float(np.interp(0, raw["b_physical_ut"], transmission)),
        "displayed_metric": next(m for m in sc.headless_observables(raw, p)["metrics"]
                                 if m["label"] == "Transmission at B=0")}
    p = dict(sc.defaults(), scan_points=201)
    before = sc.compute(p)
    original = zeeman.zeeman_manifold

    def exchange_only(*args, **kwargs):
        if kwargs.get("gamma_gg") == 2 * np.pi * p["transit_relax_khz"] * 1e3:
            kwargs = dict(kwargs, gamma_gg=0.)
        return original(*args, **kwargs)

    with patch("gabes.schemes.magneto.zeeman.zeeman_manifold", side_effect=exchange_only):
        after = sc.compute(p)
    out["paraffin_exchange_only_max_chi_relative_change"] = float(
        np.max(abs(after["chi_probe"] - before["chi_probe"]))
        / np.max(abs(before["chi_probe"])))
    result = {"provenance": "⟦AI⟧ Independent current-code checks, diagnostic model changes not applied",
              "results": out}
    Path(__file__).with_name("schemes_3_4_probes.json").write_text(
        json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(result, ensure_ascii=True, indent=2))


if __name__ == "__main__":
    main()
