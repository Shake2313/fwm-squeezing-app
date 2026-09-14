"""Dipole/degeneracy/carrier audit and conditional joint gain/noise comparison.

python -m analysis.grand_challenge.normalization_audit --output NEW.json --plot NEW.png
"""

import argparse
from dataclasses import asdict
import hashlib
import json
from pathlib import Path
import platform

import numpy as np
import scipy

from gabes import constants as c
from gabes.fwm_quantum.inputs import ReducedPowerInputs, audit_consumed_inputs, power_normalized_readout
from gabes.fwm_quantum.model import reduced_pump_noise
from gabes.fwm_quantum.normalization import normalization_audit, optical_carriers, reduced_dipoles
from gabes.fwm_quantum.readout import reduced_seeded_readout
from gabes.quantum.contracts import AnalysisFrequencyAxis, ParameterEvidence
from gabes.quantum.readout import DetectorResponse, intensity_difference_spectrum
from gabes.quantum.sidebands import SIDEBAND_MODES
from .reference.direct_readout import direct_current_psd
from .reference.uncoupled_d1 import d1_dipole_reference


ROOT = Path(__file__).resolve().parents[2]


def conditional_inputs():
    return ReducedPowerInputs(.6, 530e-6, 1e18, 1.2e-7, .0125, 8e-6,
                              2*np.pi*.9e9, -2*np.pi*8e6, c.GAMMA_GG)


def _complex(value):
    value = np.asarray(value)
    return {"real": value.real.tolist(), "imag": value.imag.tolist()}


def build_report():
    inputs = conditional_inputs()
    axis = AnalysisFrequencyAxis.from_hz(np.arange(1, 41)*1e5)
    rf = AnalysisFrequencyAxis(np.r_[-axis.omega_rad_s[::-1], 0., axis.omega_rad_s])
    detector = DetectorResponse(axis, [.85, .85], np.ones((40, 2)), 1., np.zeros(40),
                                "declared conditional detector; no measurement")
    source_detector = DetectorResponse(axis, [1., 1.], np.ones((40, 2)), 1., np.zeros(40),
                                      "cell output")
    algebra = normalization_audit()
    reference = d1_dipole_reference()
    reference_error = max(abs(v-row["uniform_manifold_mean_in_dJ2"])
        for row, ref in zip(algebra["transition_rows"], reference["transition_rows"], strict=True)
        for v in ref["mean_strength_by_q"])
    cases = []
    for name in ("historical-mixed", "historical-with-carrier-anchor", "legacy-reciprocal", "uniform-zeeman-rms"):
        if name.startswith("historical"):
            coupling = reduced_dipoles("legacy-reciprocal")
            pump = c.rabi_freq(inputs.pump_power_W, inputs.pump_waist_m)
            atom = reduced_pump_noise(pump, inputs.one_photon_rad_s,
                                     transit_rate_s_inverse=inputs.transit_rate_s_inverse)
            carriers = optical_carriers(inputs.detunings)
            if name == "historical-mixed":
                carriers = carriers-2*np.pi*algebra["F2_Fprime3_anchor_minus_centroid_hz"]
            result = reduced_seeded_readout(atom, inputs.detunings, rf, detector,
                number_density_m3=inputs.number_density_m3, uniform_area_m2=inputs.uniform_area_m2,
                optical_omega_rad_s=carriers[1:], effective_dipole_C_m=coupling.base_dipole_C_m,
                length_m=inputs.length_m, seed_power_W=inputs.seed_power_W)
            pump_matrix = pump*coupling.transition_scales[:2, 2:]
        else:
            coupling = reduced_dipoles(name)
            pump = coupling.pump_rabi_rad_s(inputs.pump_power_W, inputs.pump_waist_m)
            carriers = optical_carriers(inputs.detunings)
            result = power_normalized_readout(inputs, rf, detector, convention=name)
            pump_matrix = pump*coupling.transition_scales[:2, 2:]
        cov = result.sidebands.vacuum_covariances()
        source = intensity_difference_spectrum(cov, result.output_carrier_amplitudes_sqrt_flux,
                    source_detector, mode_labels=SIDEBAND_MODES, analysis_axis=axis)
        greater, lesser = result.main_transfer.vacuum_output()
        direct = direct_current_psd(greater, lesser, np.arange(41, 81), np.arange(39, -1, -1),
            result.output_carrier_amplitudes_sqrt_flux, detector.transmissions,
            detector.current_response, detector.balance)
        direct_error = float(np.max(abs(direct-result.spectrum.quantum_psd_A2_Hz)/direct))
        cp = [asdict(ch.audit()) for ch in result.sidebands.channels]
        transfer_audits = [result.main_transfer.audit(), result.companion_transfer.audit()]
        cases.append({"name": name, "pump_base_rabi_rad_s": pump,
            "pump_transition_rabi_rad_s": pump_matrix.tolist(),
            "weak_transition_dipoles_C_m": coupling.matrix_C_m[:2, 2:].tolist(),
            "optical_carriers_rad_s": carriers.tolist(),
            "probe_power_gain": result.probe_power_gain,
            "conjugate_power_gain": result.conjugate_power_gain,
            "conjugate_photon_flux_gain": result.conjugate_photon_flux_gain,
            "source_S_minus_db": source.quantum_db.tolist(),
            "detected_S_minus_db": result.spectrum.quantum_db.tolist(),
            "detected_quantum_psd_A2_Hz": result.spectrum.quantum_psd_A2_Hz.tolist(),
            "detected_sql_psd_A2_Hz": result.spectrum.sql_psd_A2_Hz.tolist(),
            "source_covariance": cov.tolist(), "channel_audits": cp,
            "transfer_audits": transfer_audits,
            "direct_nambu_maximum_relative_error": direct_error,
            "source_minimum_uncertainty_eigenvalue": source.minimum_covariance_uncertainty_eigenvalue,
            "main_transfer": _complex(result.main_transfer.transfer),
            "passed": bool(direct_error < 1e-9 and all(item["passed"] for item in cp+transfer_audits))})
    ledger = inputs.consumed_inputs(detector, convention="uniform-zeeman-rms")
    evidence = [ParameterEvidence(key, value, unit, "assumed", source_id="conditional normalization fixture",
                estimation_method="declared inherited value; not an independent measurement",
                applicability="single-velocity RMS surrogate")
                for key, (value, unit) in ledger["scalars"].items()]
    input_audit = audit_consumed_inputs(ledger, evidence)
    source_paths = sorted(set(
        list(ROOT.glob("gabes/quantum/*.py"))+list(ROOT.glob("gabes/fwm_quantum/*.py"))+
        [ROOT/p for p in ("gabes/constants.py", "gabes/core.py", "gabes/atoms.py", "gabes/hyperfine.py",
            "gabes/species.py", "gabes/zeeman.py", "gabes/observables.py", "gabes/schemes/fwm.py",
            "gabes/plot_style.py", "analysis/grand_challenge/normalization_audit.py",
            "analysis/grand_challenge/reference/uncoupled_d1.py",
            "analysis/grand_challenge/reference/direct_readout.py")]))
    return {"stage": "S1 physical normalization audit", "absolute_hot_vapor_prediction": False,
        "experimental_validation": False, "algebra": algebra,
        "independent_uncoupled_reference": {k: v for k, v in reference.items() if k != "dipoles"},
        "reference_strength_maximum_error": reference_error,
        "analysis_rf_hz": axis.frequency_hz.tolist(), "conditional_inputs": asdict(inputs),
        "consumed_input_ledger": ledger, "input_evidence": [asdict(e) for e in evidence],
        "independent_input_audit": asdict(input_audit), "configurations": cases,
        "expected_controls_passed": bool(reference_error < 1e-14
            and all(reference[k] < 1e-14 for k in ("basis_unitarity_error", "emission_closure_error", "absorption_closure_error"))
            and all(case["passed"] for case in cases) and not input_audit.passed),
        "environment": {"python": platform.python_version(), "numpy": np.__version__, "scipy": scipy.__version__},
        "source_sha256": {str(p.relative_to(ROOT)).replace('\\', '/'): hashlib.sha256(p.read_bytes()).hexdigest() for p in source_paths},
        "sources": [{"url": "https://steck.us/alkalidata/rubidium85numbers.pdf",
                     "sections": "Eqs. 34-45, 49-50; Tables 4, 7, 8; revision 2.3.4"}],
        "limits": ledger["assumptions"]+[
            "A scalar 1/12 after sublevel-summed strength cannot reproduce both unpolarized manifold absorption weights.",
            "RMS strengths fix this weak-absorption condition only; they do not prove full-atom nonlinear squeezing.",
            "Dipole uses the existing natural decay 2pi*5.746 MHz; no silent migration to the species table 5.750 MHz.",
            "Carrier anchor is corrected without changing Delta in the atomic Hamiltonian; it affects optical power normalization only.",
            "Input evidence is bound to consumed values/units. Fixture assumptions deliberately fail the independent-input gate."]}


def save_plot(report, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from gabes.plot_style import PALETTE, apply_gabes_plot_style
    fig, axes = plt.subplots(1, 2, figsize=(11.8, 4.9), layout="constrained")
    selected = [report["configurations"][k] for k in (0, 2, 3)]
    colors = [PALETTE["cyan"], PALETTE["rose"], "#D97706"]
    labels = ["Historical mixed", "Shared historical dipoles", "Shared Zeeman RMS"]
    for case, label, color in zip(selected, labels, colors, strict=True):
        axes[0].plot(np.array(report["analysis_rf_hz"])/1e6, case["detected_S_minus_db"],
                     label=label, color=color)
    axes[0].set(xlabel="Analysis frequency [MHz]", ylabel="S_minus / SQL [dB]", title="Declared 85% detection")
    axes[0].legend(fontsize=8)
    axes[1].bar(np.arange(3), [x["probe_power_gain"]-1 for x in selected], color=colors)
    axes[1].set_xticks(np.arange(3), ["Historical\nmixed", "Shared\nhistorical", "Shared\nZeeman RMS"])
    axes[1].set(ylabel="Coherent probe gain minus 1", title="Same power, density and detunings")
    fig.suptitle("Conditional normalization comparison\n600 mW, 530 um peak pump; single velocity; no experimental validation")
    apply_gabes_plot_style(fig)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("xb") as stream:
        fig.savefig(stream, format="png", dpi=170)
    plt.close(fig)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--plot", type=Path)
    args = parser.parse_args(argv)
    for path in (args.output, args.plot):
        if path is not None and path.exists():
            raise FileExistsError(path)
    report = build_report()
    body = json.dumps(report, ensure_ascii=False, indent=2, allow_nan=False)+"\n"
    if args.output is None:
        print(body, end="")
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        with args.output.open("x", encoding="utf-8") as stream:
            stream.write(body)
        print(f"Wrote {args.output}")
    if args.plot is not None:
        save_plot(report, args.plot)
        print(f"Wrote {args.plot}")
    return 0 if report["expected_controls_passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
