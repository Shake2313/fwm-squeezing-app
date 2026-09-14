"""Joint conditional reduced gain/S_minus and normalized-sideband readout audit.

python -m analysis.grand_challenge.readout_audit --output NEW.json --plot NEW.png
"""

import argparse
from dataclasses import asdict
import hashlib
import json
from pathlib import Path
import platform

import numpy as np
import scipy

from gabes import constants as c, hyperfine
from gabes.fwm_quantum.model import reduced_pump_noise
from gabes.fwm_quantum.readout import reduced_seeded_readout
from gabes.quantum.channels import covariance_uncertainty_minimum
from gabes.quantum.contracts import CONVENTION_ID, AnalysisFrequencyAxis, OpticalDetunings
from gabes.quantum.readout import DetectorResponse, intensity_difference_spectrum
from gabes.quantum.sidebands import SIDEBAND_MODES, TopHatBand, average_spectral_channels
from gabes.schemes import fwm
from .reference.direct_readout import direct_current_psd


ROOT = Path(__file__).resolve().parents[2]


def _complex(value):
    value = np.asarray(value)
    return {"real": value.real.tolist(), "imag": value.imag.tolist()}


def _detector(axis, eta):
    nf = len(axis.omega_rad_s)
    return DetectorResponse(axis, [eta, eta], np.ones((nf, 2)), 1., np.zeros(nf),
                             "declared equal-efficiency flat-response detector; no independent measurement claimed")


def _spectrum(result):
    return {"quantum_psd_A2_Hz": result.quantum_psd_A2_Hz.tolist(),
            "sql_psd_A2_Hz": result.sql_psd_A2_Hz.tolist(),
            "electronics_psd_A2_Hz": result.electronics_psd_A2_Hz.tolist(),
            "quantum_ratio": result.quantum_ratio.tolist(), "total_ratio": result.total_ratio.tolist(),
            "quantum_db": result.quantum_db.tolist(), "total_db": result.total_db.tolist(),
            "detected_carrier_flux_s_inverse": result.detected_carrier_flux_s_inverse.tolist(),
            "minimum_covariance_uncertainty_eigenvalue": result.minimum_covariance_uncertainty_eigenvalue,
            "maximum_readout_commutator_absolute": result.maximum_readout_commutator_absolute}


def build_report():
    optical = OpticalDetunings(2*np.pi*.9e9, -2*np.pi*8e6)
    positive = AnalysisFrequencyAxis.from_hz(np.arange(1, 41)*1e5)
    rf = AnalysisFrequencyAxis(np.r_[-positive.omega_rad_s[::-1], 0., positive.omega_rad_s])
    pump_carrier = c.OMEGA_D1+optical.one_photon_rad_s
    beat = -c.OMEGA_HF+optical.two_photon_rad_s
    settings = dict(number_density_m3=1e18, uniform_area_m2=1.2e-7,
                    optical_omega_rad_s=[pump_carrier+beat, pump_carrier-beat],
                    effective_dipole_C_m=c.DIPOLE_D1/np.sqrt(hyperfine.N_GROUND_SUBLEVELS),
                    length_m=.0125, seed_power_W=8e-6, phase_mismatch_rad_m=0.)
    configurations = []
    detector = _detector(positive, .85)
    for name, pump, reset in (
            ("pump_off_with_declared_reset", 0., c.GAMMA_GG),
            ("pumped_radiative_only", fwm.rabi_freq(.6, fwm.W_PUMP), 0.),
            ("pumped_with_declared_reset", fwm.rabi_freq(.6, fwm.W_PUMP), c.GAMMA_GG)):
        atom = reduced_pump_noise(pump, optical.one_photon_rad_s, transit_rate_s_inverse=reset)
        result = reduced_seeded_readout(atom, optical, rf, detector, **settings)
        covariance = result.sidebands.vacuum_covariances()
        source = intensity_difference_spectrum(covariance, result.output_carrier_amplitudes_sqrt_flux,
                  _detector(positive, 1.), mode_labels=SIDEBAND_MODES, analysis_axis=positive)
        greater, lesser = result.main_transfer.vacuum_output()
        direct = direct_current_psd(greater, lesser, np.arange(41, 81), np.arange(39, -1, -1),
                  result.output_carrier_amplitudes_sqrt_flux, detector.transmissions,
                  detector.current_response, detector.balance)
        reference_error = float(np.max(abs(direct-result.spectrum.quantum_psd_A2_Hz)/direct))
        equal_loss_error = float(np.max(abs(result.spectrum.quantum_ratio-(.15+.85*source.quantum_ratio))))
        cp = [asdict(channel.audit()) for channel in result.sidebands.channels]
        # Finite sampled band only: not a bound on all spontaneous fluorescence.
        mode_occupations = (np.diagonal(covariance, axis1=1, axis2=2)[:, ::2]
                            +np.diagonal(covariance, axis1=1, axis2=2)[:, 1::2]-1)/2
        band_flux = np.trapezoid(mode_occupations[:, [0, 1]]+mode_occupations[:, [2, 3]],
                                 positive.frequency_hz, axis=0)
        configurations.append({
            "name": name, "pump_rabi_rad_s": pump, "reset_rate_s_inverse": reset,
            "probe_coherent_power_gain": result.probe_power_gain,
            "conjugate_coherent_power_gain": result.conjugate_power_gain,
            "conjugate_coherent_photon_flux_gain": result.conjugate_photon_flux_gain,
            "output_coherent_powers_W": result.output_coherent_powers_W.tolist(),
            "output_carrier_amplitudes_sqrt_flux": _complex(result.output_carrier_amplitudes_sqrt_flux),
            "source_spectrum": _spectrum(source), "detected_spectrum": _spectrum(result.spectrum),
            "sideband_transfer": [channel.transfer.tolist() for channel in result.sidebands.channels],
            "sideband_added_covariance": [channel.added_covariance.tolist() for channel in result.sidebands.channels],
            "source_sideband_covariance": covariance.tolist(), "gaussian_channel_audits": cp,
            "sampled_0p1_to_4MHz_spontaneous_flux_s_inverse": band_flux.tolist(),
            "direct_nambu_psd_maximum_relative_error": reference_error,
            "equal_detection_loss_ratio_maximum_error": equal_loss_error,
            "passed": bool(reference_error < 1e-9 and equal_loss_error < 1e-10
                           and all(item["passed"] for item in cp)
                           and source.minimum_covariance_uncertainty_eigenvalue > -1e-10),
        })
    # The final atom is the unchanged pumped/reset baseline. Integrate a finite
    # 0.8--1.2 MHz wavepacket, recomputing its spectral nodes at each order.
    band_results = []
    for order in (4, 8, 16):
        band = TopHatBand(1e6, 4e5, order)
        axis, weights = band.quadrature()
        all_rf = AnalysisFrequencyAxis(np.r_[-axis.omega_rad_s[::-1], 0., axis.omega_rad_s])
        result = reduced_seeded_readout(atom, optical, all_rf, _detector(axis, .85), **settings)
        filtered = average_spectral_channels(result.sidebands.channels, weights,
                                             source="normalized flat positive band and reflected partner filters")
        center_axis = AnalysisFrequencyAxis.from_hz([band.center_hz])
        filtered_covariance = filtered.apply_covariance(np.eye(8)/2)
        measurement = intensity_difference_spectrum(filtered_covariance[None],
            result.output_carrier_amplitudes_sqrt_flux, _detector(center_axis, .85),
            mode_labels=SIDEBAND_MODES, analysis_axis=center_axis)
        expected = float(weights@result.spectrum.quantum_psd_A2_Hz)
        computed = float(measurement.quantum_psd_A2_Hz[0])
        band_results.append({
            "center_hz": band.center_hz, "bandwidth_hz": band.bandwidth_hz, "quadrature_order": order,
            "nodes_hz": axis.frequency_hz.tolist(), "normalized_weights": weights.tolist(),
            "effective_transfer": filtered.transfer.tolist(), "effective_added_covariance": filtered.added_covariance.tolist(),
            "covariance": filtered_covariance.tolist(), "channel_audit": asdict(filtered.audit()),
            "uncertainty_minimum": covariance_uncertainty_minimum(filtered_covariance),
            "band_average_quantum_psd_A2_Hz": computed,
            "unit_gain_real_bandpass_variance_A2": computed*band.bandwidth_hz,
            "quantum_ratio": float(measurement.quantum_ratio[0]),
            "direct_weighted_psd_relative_error": abs(computed-expected)/expected,
        })
    sources = [
        "gabes/atoms.py", "gabes/core.py", "gabes/constants.py", "gabes/hyperfine.py",
        "gabes/observables.py", "gabes/schemes/fwm.py", "gabes/quantum/contracts.py",
        "gabes/quantum/reservoirs.py", "gabes/quantum/diffusion.py", "gabes/quantum/traveling.py",
        "gabes/quantum/channels.py", "gabes/quantum/sidebands.py", "gabes/quantum/readout.py",
        "gabes/fwm_quantum/model.py", "gabes/fwm_quantum/field.py", "gabes/fwm_quantum/readout.py",
        "analysis/grand_challenge/reference/direct_readout.py",
        "analysis/grand_challenge/readout_audit.py", "gabes/plot_style.py", "requirements-quantum.txt",
    ]
    convergence = abs(band_results[-1]["quantum_ratio"]-band_results[-2]["quantum_ratio"])
    return {
        "schema_version": 1, "convention_id": CONVENTION_ID,
        "stage": "S1 conditional bright-carrier gain and intensity-difference readout",
        "conditional_linearized_s_minus_implemented": True,
        "absolute_hot_vapor_prediction": False, "experimental_validation": False,
        "expected_controls_passed": bool(all(case["passed"] for case in configurations)
            and all(b["channel_audit"]["passed"] and b["uncertainty_minimum"] > -1e-10
                    and b["direct_weighted_psd_relative_error"] < 1e-10 for b in band_results)
            and convergence < 1e-8),
        "environment": {"python": platform.python_version(), "numpy": np.__version__, "scipy": scipy.__version__},
        "source_sha256": {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in sources},
        "operating_inputs": {**settings, "one_photon_rad_s": optical.one_photon_rad_s,
            "two_photon_rad_s": optical.two_photon_rad_s, "analysis_rf_hz": positive.frequency_hz.tolist(),
            "detector_transmissions": detector.transmissions.tolist(), "detector_balance": detector.balance,
            "detector_response_A_per_A": _complex(detector.current_response),
            "electronics_difference_psd_A2_Hz": detector.electronics_difference_psd_A2_Hz.tolist(),
            "provenance": "unchanged conditional field-audit inputs plus declared 8 uW seed and 85% flat detection; not independent apparatus measurements"},
        "mode_order": SIDEBAND_MODES, "quadrature_order": "interleaved x,p; vacuum I/2",
        "configurations": configurations, "normalized_band_convergence": band_results,
        "band_order_8_to_16_ratio_change": convergence,
        "limits": [
            "Single velocity, selected transitions/sectors, one uniform area and prescribed pump; no angular-Doppler, full Zeeman, collection or depletion prediction.",
            "The 1/12 weak-field normalization and legacy pump-power/Rabi conversion remain conditional; no measured-input/no-fit validation badge is issued.",
            "Bright-carrier linearization drops delta-a-dagger delta-a photocurrent and spontaneous mean powers. Sampled-band photon flux is only a partial diagnostic; omitted fluorescence is not bounded.",
            "Positive RF uses four sidebands; missing reflected frequencies or companion symmetry fail explicitly. RF=0 is used only for coherent carrier transfer.",
            "Finite-band projection uses flat normalized filters, vacuum orthogonal input modes and a declared quadrature; arbitrary filter phases or correlated external spectral inputs need a broader adapter.",
            "Electronic response is calibrated A/A; separately declared electronics PSD is added, never silently subtracted. No measured RBW/VBW transfer, imbalance fit or target squeezing coefficient enters.",
        ],
    }


def save_plot(report, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from gabes.plot_style import PALETTE, apply_gabes_plot_style

    case = report["configurations"][-1]
    freq = np.array(report["operating_inputs"]["analysis_rf_hz"])/1e6
    source, detected = case["source_spectrum"], case["detected_spectrum"]
    fig, axes = plt.subplots(2, 1, figsize=(8.8, 6.8), sharex=True, layout="constrained")
    axes[0].plot(freq, source["quantum_db"], color=PALETTE["cyan"], label="Cell output (conditional)")
    axes[0].plot(freq, detected["quantum_db"], color=PALETTE["rose"], label="Declared 85% detection")
    axes[0].axhline(0., color="#64748B", ls="--", label="SQL")
    axes[0].set_ylabel("Intensity difference / SQL [dB]")
    lo = min(0., *source["quantum_db"], *detected["quantum_db"])
    hi = max(0., *source["quantum_db"], *detected["quantum_db"])
    margin = max(.03, .12*(hi-lo))
    axes[0].set_ylim(lo-margin, hi+margin)
    axes[0].legend(loc="center right")
    axes[1].plot(freq, np.array(detected["quantum_psd_A2_Hz"])/1e-24,
                 color=PALETTE["cyan"], label="Detected quantum PSD")
    axes[1].plot(freq, np.array(detected["sql_psd_A2_Hz"])/1e-24,
                 color=PALETTE["rose"], ls="--", label="Same-current SQL")
    axes[1].set_xlabel("Analysis frequency [MHz]")
    axes[1].set_ylabel("One-sided PSD [1e-24 A^2/Hz]")
    axes[1].legend(loc="center right")
    axes[1].set_xlim(.1, 4.)
    fig.suptitle("Conditional reduced FWM readout\nSingle velocity, prescribed pump, 8 uW coherent seed")
    axes[0].set_title(f"Coherent gains: probe {case['probe_coherent_power_gain']:.5f}; "
                       f"conjugate {case['conjugate_coherent_power_gain']:.5f}", fontsize=10)
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
