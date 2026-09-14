"""Gaussian quadratic current, optical-band convergence and local seed back-action.

python -m analysis.grand_challenge.seed_validity_audit --output NEW.json --plot NEW.png
"""

import argparse
from dataclasses import asdict
import hashlib
import json
from pathlib import Path
import platform

import numpy as np
import scipy

from gabes import constants as c, core
from gabes.fwm_quantum.field import reduced_field_pair, reduced_readout_operators
from gabes.fwm_quantum.inputs import power_normalized_readout
from gabes.fwm_quantum.model import reduced_pump_noise
from gabes.fwm_quantum.normalization import reduced_dipoles, optical_carriers
from gabes.fwm_quantum.seed_validity import seed_harmonic_hamiltonian, local_seed_backaction
from gabes.quantum.contracts import AnalysisFrequencyAxis
from gabes.quantum.photocounting import (PairedOpticalBins, paired_optical_bins,
    gaussian_count_covariance, quadratic_readout_correction, corrected_gaussian_readout)
from gabes.quantum.readout import DetectorResponse, intensity_difference_spectrum, coherent_carrier_output
from gabes.quantum.sidebands import SIDEBAND_MODES
from gabes.quantum.traveling import constant_segment
from .normalization_audit import conditional_inputs
from .reference.fock_photocurrent import two_pair_current_reference


ROOT = Path(__file__).resolve().parents[2]


def _complex(a):
    a = np.asarray(a)
    return {"real": a.real.tolist(), "imag": a.imag.tolist()}


def _field(atom, inputs, dipoles, axis):
    return reduced_field_pair(atom, inputs.detunings, axis, number_density_m3=inputs.number_density_m3,
        uniform_area_m2=inputs.uniform_area_m2, optical_omega_rad_s=optical_carriers(inputs.detunings)[1:],
        effective_dipole_C_m=dipoles.base_dipole_C_m, transition_scales=dipoles.transition_scales,
        phase_mismatch_rad_m=inputs.phase_mismatch_rad_m)


def build_report():
    with core.blas_single_thread():
        return _build_report()


def _build_report():
    inputs = conditional_inputs()
    dipoles = reduced_dipoles("uniform-zeeman-rms")
    omega = optical_carriers(inputs.detunings)[1:]
    atom = reduced_pump_noise(dipoles.pump_rabi_rad_s(inputs.pump_power_W, inputs.pump_waist_m),
        inputs.one_photon_rad_s, transit_rate_s_inverse=inputs.transit_rate_s_inverse,
        transition_scales=dipoles.transition_scales)
    positive = AnalysisFrequencyAxis.from_hz(np.arange(1, 41)*1e5)
    rf = AnalysisFrequencyAxis(np.r_[-positive.omega_rad_s[::-1], 0., positive.omega_rad_s])
    detector = DetectorResponse(positive, [.85, .85], np.ones((40, 2)), 1., np.zeros(40), "declared 85% flat detector")
    result = power_normalized_readout(inputs, rf, detector)
    baseline_covariance = result.sidebands.vacuum_covariances()

    def collect(half_width, width):
        count = int(round(2*half_width/width))
        centers = (np.arange(count)-(count-1)/2)*width
        axis = AnalysisFrequencyAxis.from_hz(centers)
        field = _field(atom, inputs, dipoles, axis)
        main, companion = constant_segment(field.main, inputs.length_m), constant_segment(field.companion, inputs.length_m)
        bins = paired_optical_bins(main, companion, axis, bin_width_hz=width,
            source=f"declared ideal optical top-hat in each arm, offsets +/-{half_width:g} Hz; vacuum outside")
        correction = quadratic_readout_correction(bins, detector)
        full = corrected_gaussian_readout(result.spectrum, correction)
        summary = {"optical_half_width_hz": half_width, "bin_width_hz": width, "bin_count": count,
            "source_spontaneous_flux_s_inverse": bins.spontaneous_flux_s_inverse.tolist(),
            "detected_spontaneous_flux_s_inverse": correction.spontaneous_detected_flux_s_inverse.tolist(),
            "quadratic_psd_A2_Hz": correction.quadratic_psd_A2_Hz.tolist(),
            "added_sql_psd_A2_Hz": correction.added_sql_psd_A2_Hz.tolist(),
            "corrected_S_minus_db": full["quantum_db"].tolist(),
            "maximum_correction_db_at_8uW": float(np.max(abs(full["quantum_db"]-result.spectrum.quantum_db))),
            "bin_moment_audit": bins.audit(),
            "maximum_propagated_commutator_residual": max(main.audit()["maximum_commutator_relative_residual"],
                                                          companion.audit()["maximum_commutator_relative_residual"])}
        return bins, correction, summary

    bandwidth = []
    for half_width in (8e6, 32e6, 128e6, 512e6, 1024e6):
        _, _, summary = collect(half_width, .5e6)
        bandwidth.append(summary)
    convergence = []
    for width in (1e6, .5e6, .25e6, .125e6, .0625e6, .03125e6):
        bins, correction, summary = collect(1024e6, width)
        convergence.append(summary)
    relative_grid_change = max(
        float(np.max(abs(np.array(convergence[-1][key])-convergence[-2][key])/
                         np.maximum(abs(np.array(convergence[-1][key])), np.finfo(float).tiny)))
        for key in ("source_spontaneous_flux_s_inverse", "quadratic_psd_A2_Hz"))
    bandwidth_flux_change = float(np.max(abs(np.array(bandwidth[-1]["source_spontaneous_flux_s_inverse"])
        /np.array(bandwidth[-2]["source_spontaneous_flux_s_inverse"])-1)))

    # Mean carriers are prescribed from the weak-field propagation at z=0,L/2,L.
    dc = AnalysisFrequencyAxis.from_hz([0.])
    local_field = _field(atom, inputs, dipoles, dc)
    seed_input = [np.sqrt(inputs.seed_power_W/(c.HBAR*omega[0])), 0.]
    local_carriers = [coherent_carrier_output(constant_segment(local_field.main, z), dc, seed_input)
                      for z in (0., inputs.length_m/2, inputs.length_m)]
    beat = -c.OMEGA_HF+inputs.two_photon_rad_s
    operators = reduced_readout_operators(dipoles.transition_scales)
    sweep = []
    powers = np.unique(np.r_[np.logspace(-14, -2, 13), inputs.seed_power_W])
    baseline_details = []
    for power in powers:
        ratio = power/inputs.seed_power_W
        bright = intensity_difference_spectrum(baseline_covariance,
            result.output_carrier_amplitudes_sqrt_flux*np.sqrt(ratio), detector,
            mode_labels=SIDEBAND_MODES, analysis_axis=positive)
        full = corrected_gaussian_readout(bright, correction)
        checks = []
        for fraction, beta in zip((0., .5, 1.), local_carriers):
            v = seed_harmonic_hamiltonian(dipoles, omega, inputs.uniform_area_m2, beta*np.sqrt(ratio))
            ordered = [local_seed_backaction(atom, v, beat, n_f=nf) for nf in (2, 3, 4)]
            last = ordered[-1]
            weak_pol = np.einsum("iab,ba->i", operators, last["weak_harmonic"])
            finite_pol = np.einsum("iab,ba->i", operators, last["harmonics"][5])
            pol_error = float(np.linalg.norm(finite_pol-weak_pol)/np.linalg.norm(weak_pol))
            order_change = float(np.linalg.norm(last["harmonics"][1:-1]-ordered[-2]["harmonics"]))
            check = {"z_fraction": fraction, "mean_state_trace_distance": last["mean_state_trace_distance"],
                "first_harmonic_relative_error": last["first_harmonic_relative_error"],
                "polarization_relative_error": pol_error,
                "order_3_to_4_harmonics_absolute_change": order_change,
                "diagnostics": last["diagnostics"]}
            checks.append(check)
            if power == inputs.seed_power_W:
                dense = core.floquet_solve_direct(atom.generator, core.comm_super(v), core.comm_super(v.conj().T),
                    beat, [0.], np.zeros((16, 16)), 4, n_f=4, return_harmonics=True)[0]
                baseline_details.append({**check, "prescribed_carriers_sqrt_flux": _complex(beta),
                    "harmonics": _complex(last["harmonics"]),
                    "dense_floquet_maximum_absolute_difference": float(np.max(abs(dense-last["harmonics"])))})
        db_change = float(np.max(abs(full["quantum_db"]-bright.quantum_db)))
        fraction_sql = float(np.max(correction.added_sql_psd_A2_Hz/bright.sql_psd_A2_Hz))
        pol_max = max(check["polarization_relative_error"] for check in checks)
        sweep.append({"seed_power_W": float(power), "bright_S_minus_db": bright.quantum_db.tolist(),
            "corrected_S_minus_db": full["quantum_db"].tolist(),
            "maximum_change_db": db_change, "spontaneous_to_coherent_sql_ratio": fraction_sql,
            "quadratic_to_bright_psd_maximum_ratio": float(np.max(correction.quadratic_psd_A2_Hz/bright.quantum_psd_A2_Hz)),
            "detected_total_flux_s_inverse": full["detected_total_flux_s_inverse"].tolist(),
            "local_backaction": checks, "maximum_local_polarization_relative_error": pol_max,
            "passes_declared_diagnostic_thresholds": bool(db_change < .01 and fraction_sql < .01 and pol_max < .01)})

    n = np.array([.1, .2])
    phase = np.array([.2, -.7])
    fixture = PairedOpticalBins([-.5, .5], 1., np.column_stack([n, n[::-1]]), np.sqrt(n*(n+1))*np.exp(1j*phase), "Fock control")
    predicted = gaussian_count_covariance(fixture, AnalysisFrequencyAxis.from_hz([1.]))[0]
    fock = two_pair_current_reference(*n, *phase, cutoff=20)
    fock_error = float(np.max(abs(predicted-fock["count_covariance_per_hz"])))
    all_checks = [check for point in sweep for check in point["local_backaction"]]
    controls = bool(relative_grid_change < 1e-4 and fock_error < 1e-11
        and all(x["bin_moment_audit"]["passed"] for x in bandwidth+convergence)
        and all(x["order_3_to_4_harmonics_absolute_change"] < 1e-8
                and x["diagnostics"]["cancellation_aware_residual_passed"]
                and x["diagnostics"]["phase_sampled_minimum_state_eigenvalue"] > -1e-10 for x in all_checks)
        and all(x["dense_floquet_maximum_absolute_difference"] < 1e-11 for x in baseline_details))
    sources = sorted(set(list(ROOT.glob("gabes/quantum/*.py"))+list(ROOT.glob("gabes/fwm_quantum/*.py"))+
        [ROOT/p for p in ("gabes/core.py", "gabes/constants.py", "gabes/atoms.py", "gabes/hyperfine.py",
            "gabes/species.py", "gabes/zeeman.py", "gabes/observables.py", "gabes/schemes/fwm.py", "gabes/plot_style.py",
            "analysis/grand_challenge/normalization_audit.py", "analysis/grand_challenge/seed_validity_audit.py",
            "analysis/grand_challenge/reference/fock_photocurrent.py")]))
    return {"stage": "S1 Gaussian photocount correction and local finite-seed validity",
        "expected_controls_passed": controls, "absolute_hot_vapor_prediction": False, "experimental_validation": False,
        "full_unfiltered_fluorescence_bounded": False, "finite_seed_quantum_noise_implemented": False,
        "conditional_inputs": asdict(inputs), "dipole_convention": dipoles.convention,
        "analysis_rf_hz": positive.frequency_hz.tolist(), "detector_transmissions": detector.transmissions.tolist(),
        "optical_bandwidth_scan": bandwidth, "fixed_band_grid_convergence": convergence,
        "maximum_relative_grid_change": relative_grid_change,
        "512_to_1024MHz_relative_collected_flux_change": bandwidth_flux_change,
        "selected_optical_bins": {"centers_hz": bins.centers_hz.tolist(), "bin_width_hz": bins.width_hz,
            "occupations": bins.occupations.tolist(), "anomalous_pc": _complex(bins.anomalous_pc)},
        "seed_power_sweep": sweep, "baseline_local_floquet_details": baseline_details,
        "fock_reference": {"cutoff": 20, "maximum_absolute_covariance_error": fock_error,
                           "omitted_probability_bound": fock["omitted_probability_bound"]},
        "diagnostic_thresholds": {"maximum_readout_change_db": .01, "spontaneous_to_coherent_sql": .01,
            "maximum_local_polarization_relative_error": .01,
            "scope": "declared diagnostic tolerances, not rigorous full-model error bounds or fitted parameters"},
        "environment": {"python": platform.python_version(), "numpy": np.__version__, "scipy": scipy.__version__},
        "source_sha256": {str(p.relative_to(ROOT)).replace('\\', '/'): hashlib.sha256(p.read_bytes()).hexdigest() for p in sources},
        "sources": ["https://pmc.ncbi.nlm.nih.gov/articles/PMC11667585/", "https://arxiv.org/abs/1007.1610"],
        "limits": [
            "Gaussian Wick factorization closes fourth moments; connected non-Gaussian atomic/field cumulants are not bounded.",
            "Finite optical top-hat collection is declared, not measured. Vacuum outside restores the full spontaneous shot term.",
            "Optical band changes are physical collection changes. Fixed-band numerical convergence does not bound unfiltered tails.",
            "Only the selected two-mode sectors are used. Far-away resonances, other collected modes/polarizations and technical pump noise are not included.",
            "Floquet means use prescribed carriers at three positions. This is a back-action diagnostic, not self-consistent finite-seed propagation or Floquet Langevin diffusion.",
            "Uniform Zeeman RMS, single velocity, Gaussian pump peak, declared reset and undepleted pump retain their previous limitations."]}


def save_plot(report, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from gabes.plot_style import PALETTE, apply_gabes_plot_style
    fig, axes = plt.subplots(2, 2, figsize=(11.4, 8.2), layout="constrained")
    sweep = report["seed_power_sweep"]
    powers = np.array([r["seed_power_W"] for r in sweep])*1e6
    axes[0, 0].semilogx(powers, [r["bright_S_minus_db"][9] for r in sweep], color=PALETTE["cyan"], label="Bright carrier")
    axes[0, 0].semilogx(powers, [r["corrected_S_minus_db"][9] for r in sweep], color=PALETTE["rose"], label="Gaussian quadratic + full SQL")
    axes[0, 0].set(xlabel="Seed power [uW]", ylabel="S_minus at 1 MHz [dB]", title="Declared +/-1.024 GHz optical bands")
    axes[0, 0].legend(fontsize=8)
    axes[0, 1].loglog(powers, [r["spontaneous_to_coherent_sql_ratio"] for r in sweep], color=PALETTE["cyan"], label="Spontaneous / coherent SQL")
    axes[0, 1].loglog(powers, [r["maximum_local_polarization_relative_error"] for r in sweep], color=PALETTE["rose"], label="Local finite-seed polarization error")
    axes[0, 1].axhline(.01, color="#64748B", ls="--", label="Declared 1% diagnostic")
    axes[0, 1].set(xlabel="Seed power [uW]", ylabel="Ratio", title="Distinct lower and upper seed limits")
    axes[0, 1].legend(fontsize=8)
    bands = report["optical_bandwidth_scan"]
    for j, (label, color) in enumerate(zip(("Probe", "Conjugate"), (PALETTE["cyan"], PALETTE["rose"]))):
        axes[1, 0].semilogx([r["optical_half_width_hz"]/1e6 for r in bands],
            [r["source_spontaneous_flux_s_inverse"][j]/1e6 for r in bands], "o-", color=color, label=label)
    axes[1, 0].set(xlabel="Optical collection half-width [MHz]", ylabel="Spontaneous flux [million / s]", title="Collection-band dependence; no tail bound")
    axes[1, 0].legend(fontsize=8)
    optical = report["selected_optical_bins"]
    for j, (label, color) in enumerate(zip(("Probe", "Conjugate"), (PALETTE["cyan"], PALETTE["rose"]))):
        axes[1, 1].semilogy(np.array(optical["centers_hz"])/1e6, np.array(optical["occupations"])[:, j], color=color, label=label)
    axes[1, 1].set(xlabel="Optical envelope offset [MHz]", ylabel="Spectral occupation", title="Beyond the 0.1-4 MHz analysis band")
    axes[1, 1].legend(fontsize=8)
    fig.suptitle("Conditional Gaussian readout and seed-validity audit\nUniform Zeeman RMS; prescribed pump; no experimental validation")
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
