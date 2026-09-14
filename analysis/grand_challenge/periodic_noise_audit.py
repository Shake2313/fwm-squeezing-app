"""Reduced finite-seed atomic Floquet drift/diffusion and independent cyclic QRT.

python -m analysis.grand_challenge.periodic_noise_audit --output NEW.json --plot NEW.png
"""

import argparse
from dataclasses import replace
import hashlib
import json
from pathlib import Path
import platform

import numpy as np
import scipy

from gabes import core, constants as c
from gabes.fwm_quantum.field import reduced_readout_operators
from gabes.fwm_quantum.inputs import INPUT_UNITS, power_normalized_readout
from gabes.fwm_quantum.model import reduced_pump_system
from gabes.fwm_quantum.normalization import reduced_dipoles, optical_carriers
from gabes.fwm_quantum.periodic import reduced_periodic_noise
from gabes.quantum.contracts import AnalysisFrequencyAxis, GeneratorFrequencyAxis
from gabes.quantum.diffusion import stationary_atomic_noise
from gabes.quantum.readout import DetectorResponse
from .normalization_audit import conditional_inputs
from .reference.atomic_qrt import adjoint_einstein_diffusion
from .reference.periodic_qrt import time_domain_qrt


ROOT = Path(__file__).resolve().parents[2]


def _complex(value):
    value = np.asarray(value)
    return {"real": value.real.tolist(), "imag": value.imag.tolist()}


def _crop(spectrum, output=(-1, 0, 1), ordering="greater"):
    return np.concatenate([np.concatenate([spectrum.block(i, j, ordering=ordering) for j in output], axis=2) for i in output], axis=1)


def _relative(a, b):
    return float(np.linalg.norm(a-b)/max(np.linalg.norm(b), np.finfo(float).tiny))


def build_report():
    with core.blas_single_thread():
        return _build_report()


def _build_report():
    inputs = conditional_inputs()
    dipoles = reduced_dipoles("uniform-zeeman-rms")
    omega = optical_carriers(inputs.detunings)[1:]
    carrier = np.array([np.sqrt(inputs.seed_power_W/(c.HBAR*omega[0])), 0.])
    rf = AnalysisFrequencyAxis.from_hz([-1e6, 0., 1e6])
    detector = DetectorResponse(AnalysisFrequencyAxis.from_hz([1e6]), [.85, .85], [[1., 1.]], 1., [0.], "declared detector for prescribed mean only")
    mean = power_normalized_readout(inputs, rf, detector)
    exit_carrier = mean.output_carrier_amplitudes_sqrt_flux
    axis = GeneratorFrequencyAxis(2*np.pi*1e6*np.array([-4., -1., 0., .1, .5, 1., 2., 3., 4.]),
                                  "rb85-static-pump-Floquet-quasifrequency")
    reference_axis = GeneratorFrequencyAxis([0., 2*np.pi*1e6], axis.frame)
    cases = []
    models = []
    for name, beta in (("zero-seed", carrier*0), ("8uW-cell-entrance", carrier),
                       ("8uW-prescribed-cell-exit", exit_carrier),
                       ("1mW-prescribed-cell-exit", exit_carrier*np.sqrt(1e-3/inputs.seed_power_W))):
        model_orders = [reduced_periodic_noise(inputs, beta, mean_order=order) for order in (3, 4, 5)]
        model = model_orders[-1]
        models.append(model)
        mean_change = float(np.max(abs(model.state_harmonics[1:-1]-model_orders[-2].state_harmonics)))
        spectra = [model.lift(order).spectrum(reference_axis) for order in (1, 2, 3, 4, 5)]
        order_changes = [_relative(_crop(b), _crop(a)) for a, b in zip(spectra[:-1], spectra[1:])]
        lift = model.lift(5)
        spectrum = lift.spectrum(axis)
        blocks = _crop(spectrum)
        reference = time_domain_qrt(model.reservoirs.generator(model.hamiltonian_zero_rad_s),
            core.comm_super(model.hamiltonian_plus_rad_s), core.comm_super(model.hamiltonian_plus_rad_s.conj().T),
            model.beat_rad_s, model.operators, reference_axis.omega_rad_s, [-1, 0, 1],
            phase_samples=16, rtol=2e-12, atol=2e-14)
        ref_error = _relative(_crop(spectra[-1]), reference["ordered_spectrum"])
        state_error = float(np.max(abs(model.at_phase(reference["physical_phases_rad"])["state"]-reference["phase_states"])))
        # Independent adjoint product-rule Einstein check at an off-grid phase.
        phi = .371
        sample = model.at_phase([phi])
        H = model.hamiltonian_zero_rad_s+model.hamiltonian_plus_rad_s*np.exp(-1j*phi)+model.hamiltonian_plus_rad_s.conj().T*np.exp(1j*phi)
        independent_d = adjoint_einstein_diffusion(model.reservoirs.generator(H), sample["state"][0], model.operators)
        d_error = _relative(sample["diffusion_by_reservoir"][:, 0].sum(axis=0), independent_d)
        # Diagnostic ablation: same periodic A, discard D_q for q!=0.
        center = len(model.state_harmonics)//2
        averaged = replace(lift,
            greater_by_reservoir=np.array([np.kron(np.eye(len(lift.harmonics)), d[center]) for d in model.diffusion_harmonics_by_reservoir]),
            lesser_by_reservoir=np.array([np.kron(np.eye(len(lift.harmonics)), d[center].T) for d in model.diffusion_harmonics_by_reservoir]))
        averaged_spectrum = averaged.spectrum(axis)
        average_error = _relative(_crop(averaged_spectrum), blocks)
        harmonics_d = model.diffusion_harmonics_by_reservoir.sum(axis=0)
        time_dependent_noise = float(np.linalg.norm(np.delete(harmonics_d, center, axis=0))/np.linalg.norm(harmonics_d))
        m = len(model.operators)
        diagonal = np.zeros_like(blocks)
        for k in range(3):
            diagonal[:, k*m:(k+1)*m, k*m:(k+1)*m] = blocks[:, k*m:(k+1)*m, k*m:(k+1)*m]
        cross_ratio = float(np.linalg.norm(blocks-diagonal)/np.linalg.norm(blocks))
        op = reduced_readout_operators(dipoles.transition_scales)[0]
        weight = np.einsum("ab,iba->i", op, model.operators)
        projected = np.einsum("i,fij,j->f", weight, spectrum.block(1, 1, ordering="lesser"), weight.conj()).real
        projected_avg = np.einsum("i,fij,j->f", weight, averaged_spectrum.block(1, 1, ordering="lesser"), weight.conj()).real
        cases.append({"name": name, "prescribed_carriers_sqrt_flux": _complex(beta),
            "mean_orders": [3, 4, 5], "mean_order_4_to_5_maximum_change": mean_change,
            "mean_diagnostics": model.diagnostics,
            "response_orders": [1, 2, 3, 4, 5], "adjacent_response_interior_relative_changes": order_changes,
            "response_lift_audit": lift.audit(), "spectrum_audit": spectrum.audit(),
            "independent_time_qrt_relative_error": ref_error, "independent_periodic_state_maximum_error": state_error,
            "independent_product_rule_diffusion_relative_error": d_error,
            "time_dependent_diffusion_harmonics_relative_norm": time_dependent_noise,
            "cross_harmonic_spectrum_relative_norm": cross_ratio,
            "averaged_noise_spectrum_relative_change": average_error,
            "probe_dipole_lesser_spectrum_s": projected.tolist(),
            "probe_dipole_lesser_spectrum_with_averaged_noise_s": projected_avg.tolist(),
            "state_harmonics": _complex(model.state_harmonics), "atomic_drift_harmonics": _complex(model.drift_harmonics),
            "diffusion_harmonics_by_reservoir": _complex(model.diffusion_harmonics_by_reservoir),
            "reservoir_names": list(lift.reservoir_names),
            "ordered_spectrum_interior": _complex(blocks), "lesser_spectrum_interior": _complex(_crop(spectrum, ordering="lesser")),
            "time_domain_qrt_interior": _complex(reference["ordered_spectrum"]),
            "passed": bool(mean_change < 1e-10 and order_changes[-1] < 1e-8 and ref_error < 1e-7
                            and state_error < 1e-9 and d_error < 1e-10 and spectrum.audit()["passed"])})
    # Phase integration and ODE tolerance convergence on the 8uW exit fixture.
    model = models[2]
    references = []
    for samples, tol in ((16, 2e-10), (16, 2e-12), (32, 2e-12)):
        references.append(time_domain_qrt(model.reservoirs.generator(model.hamiltonian_zero_rad_s),
            core.comm_super(model.hamiltonian_plus_rad_s), core.comm_super(model.hamiltonian_plus_rad_s.conj().T),
            model.beat_rad_s, model.operators, reference_axis.omega_rad_s, [-1, 0, 1],
            phase_samples=samples, rtol=tol, atol=tol/100))
    phase_change = _relative(references[-1]["ordered_spectrum"], references[-2]["ordered_spectrum"])
    tolerance_change = _relative(references[1]["ordered_spectrum"], references[0]["ordered_spectrum"])
    h, reservoirs = reduced_pump_system(dipoles.pump_rabi_rad_s(inputs.pump_power_W, inputs.pump_waist_m),
        inputs.one_photon_rad_s, transit_rate_s_inverse=inputs.transit_rate_s_inverse, transition_scales=dipoles.transition_scales)
    stationary = stationary_atomic_noise(h, reservoirs)
    zero = models[0].lift(2).spectrum(axis)
    zero_errors = []
    for harmonic in (-1, 0, 1):
        reference = stationary.spectrum(GeneratorFrequencyAxis(axis.omega_rad_s+harmonic*model.beat_rad_s, "stationary shifted reference"))
        zero_errors.append(_relative(zero.block(harmonic, harmonic), reference.ordered))
    sources = sorted(set(list(ROOT.glob("gabes/quantum/*.py"))+list(ROOT.glob("gabes/fwm_quantum/*.py"))+
        [ROOT/p for p in ("gabes/core.py", "gabes/atoms.py", "gabes/constants.py", "gabes/hyperfine.py", "gabes/species.py",
            "gabes/zeeman.py", "gabes/observables.py", "gabes/schemes/fwm.py", "gabes/plot_style.py",
            "analysis/grand_challenge/normalization_audit.py", "analysis/grand_challenge/periodic_noise_audit.py",
            "analysis/grand_challenge/reference/atomic_qrt.py", "analysis/grand_challenge/reference/periodic_qrt.py")]))
    return {"stage": "S1 periodic microscopic atomic drift and harmonic-correlated diffusion",
        "periodic_atomic_quantum_noise_implemented": True, "finite_seed_field_quantum_noise_implemented": False,
        "absolute_hot_vapor_prediction": False, "experimental_validation": False,
        "expected_controls_passed": bool(all(c["passed"] for c in cases) and phase_change < 1e-8
            and tolerance_change < 1e-7 and max(zero_errors) < 1e-10),
        "quasifrequencies_rad_s": axis.omega_rad_s.tolist(), "reference_quasifrequencies_rad_s": reference_axis.omega_rad_s.tolist(),
        "beat_rad_s": model.beat_rad_s, "output_harmonics": [-1, 0, 1], "atomic_operator_order": "complete traceless Hermitian basis",
        "prescribed_inputs": {k: {"value": getattr(inputs, k), "unit": unit, "status": "assumed"}
                              for k, unit in INPUT_UNITS.items()},
        "carrier_scope": "8uW fixture defines the entrance; exit carriers use the prior weak-field mean; 1mW scales that prescribed exit by sqrt(125).",
        "dipole_convention": dipoles.convention, "configurations": cases,
        "zero_seed_maximum_relative_spectral_error": max(zero_errors),
        "time_domain_reference_phase_16_to_32_relative_change": phase_change,
        "time_domain_reference_rtol_2e10_to_2e12_relative_change": tolerance_change,
        "environment": {"python": platform.python_version(), "numpy": np.__version__, "scipy": scipy.__version__},
        "source_sha256": {str(p.relative_to(ROOT)).replace('\\', '/'): hashlib.sha256(p.read_bytes()).hexdigest() for p in sources},
        "sources": ["https://arxiv.org/abs/2005.08249", "https://arxiv.org/abs/1512.05562"],
        "limits": [
            "Atomic K(t)=<commutator> is state-dependent, not the canonical field commutator J.",
            "Mean and response harmonic orders are independent cutoffs; compare retained interior blocks, not moving outer boundaries.",
            "Time-domain QRT retains a one-period propagator and geometric infinite-period tail; coherent periodic means are subtracted.",
            "Only explicit constant jump operators and one sinusoidal Hamiltonian drive are currently supported.",
            "Full finite-seed field elimination must include both Nambu sectors and additional harmonics/mode labels; the previous two-sector field/readout is not reused as a finite-seed result.",
            "Prescribed local carriers, reduced RMS dipoles, single velocity and declared reservoirs remain conditional; no depletion, full atom or experimental validation."]}


def save_plot(report, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from gabes.plot_style import PALETTE, apply_gabes_plot_style
    fig, axes = plt.subplots(1, 2, figsize=(11.2, 4.8), layout="constrained")
    f = np.array(report["quasifrequencies_rad_s"])/(2*np.pi*1e6)
    mask = f > 0
    baseline = np.array(report["configurations"][0]["probe_dipole_lesser_spectrum_s"])
    for case, label, color in zip(report["configurations"][1:], ("8 uW entrance", "8 uW prescribed exit", "1 mW prescribed exit"),
                                   (PALETTE["cyan"], PALETTE["rose"], "#D97706")):
        axes[0].semilogy(f[mask], np.array(case["probe_dipole_lesser_spectrum_s"])[mask]/baseline[mask],
                     "o-", color=color, label=label)
    axes[0].axhline(1., color=PALETTE["muted"], linestyle="--", linewidth=1)
    axes[0].set(xlabel="Base frequency [MHz]", ylabel="Atomic probe spectrum / zero-seed spectrum",
                title="h=1 atomic lesser spectrum vs zero seed")
    axes[0].legend(fontsize=8)
    case = report["configurations"][2]
    m = 15
    block = case["ordered_spectrum_interior"]
    matrix = np.array(block["real"])[5]+1j*np.array(block["imag"])[5]
    norms = np.array([[np.linalg.norm(matrix[i*m:(i+1)*m, j*m:(j+1)*m]) for j in range(3)] for i in range(3)])
    normalized = norms/np.sqrt(np.diag(norms)[:, None]*np.diag(norms)[None, :])
    plot = axes[1].imshow(normalized, vmin=0, vmax=1, cmap="Blues")
    for i in range(3):
        for j in range(3):
            axes[1].text(j, i, f"{normalized[i,j]:.3f}", ha="center", va="center",
                         color="white" if normalized[i,j] > .6 else "#172554", fontsize=10)
    axes[1].set_xticks(range(3), [-1, 0, 1]); axes[1].set_yticks(range(3), [-1, 0, 1])
    axes[1].set(xlabel="Harmonic k", ylabel="Harmonic h", title="8 uW exit: cross-harmonic block norms at 1 MHz")
    fig.colorbar(plot, ax=axes[1], label="Normalized block norm")
    fig.suptitle("Periodic microscopic atomic noise\nConditional reduced atom; field squeezing is not computed here")
    apply_gabes_plot_style(fig)
    axes[1].grid(False)
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
    text = json.dumps(report, ensure_ascii=False, indent=2, allow_nan=False)+"\n"
    if args.output is None:
        print(text, end="")
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        with args.output.open("x", encoding="utf-8") as stream:
            stream.write(text)
        print(f"Wrote {args.output}")
    if args.plot is not None:
        save_plot(report, args.plot)
        print(f"Wrote {args.plot}")
    return 0 if report["expected_controls_passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
