"""S1 atomic diffusion/QRT audit; no optical field squeezing is computed.

python -m analysis.grand_challenge.atomic_noise_audit --output NEW_REPORT.json
"""

import argparse
from dataclasses import asdict, replace
import hashlib
import json
from pathlib import Path
import platform

import numpy as np

from gabes import atoms, constants
from gabes.fwm_quantum.model import reduced_pump_noise, minus_branch_atomic_frequencies
from gabes.quantum.contracts import (
    CONVENTION_ID, AnalysisFrequencyAxis, GeneratorFrequencyAxis, OpticalDetunings,
)
from gabes.quantum.reservoirs import thermal_reset_channels
from gabes.schemes import fwm
from .reference.atomic_qrt import adjoint_einstein_diffusion, qrt_ordered_spectrum


ROOT = Path(__file__).resolve().parents[2]


def _complex_array(value):
    value = np.asarray(value)
    return {"real": value.real.tolist(), "imag": value.imag.tolist()}


def build_report():
    pump = fwm.rabi_freq(0.6, fwm.W_PUMP)
    optical = OpticalDetunings(2*np.pi*0.9e9, -2*np.pi*8e6)
    # Positive RF coordinates remain separate from the shifted atomic frequency.
    rf = AnalysisFrequencyAxis.from_hz(np.linspace(0.1e6, 4e6, 40))
    axes = (
        GeneratorFrequencyAxis(2*np.pi*np.array([-4e6, -1e6, 0, 1e6, 4e6]),
                               "rb85-d1-static-pump-near-zero"),
        minus_branch_atomic_frequencies(optical, rf),
    )
    configurations = []
    for name, pump_rabi, reset_rate in (
        ("pump_off_with_declared_reset", 0., constants.GAMMA_GG),
        ("pumped_radiative_only", pump, 0.),
        ("pumped_with_declared_reset", pump, constants.GAMMA_GG),
    ):
        model = reduced_pump_noise(pump_rabi, optical.one_photon_rad_s,
                                   transit_rate_s_inverse=reset_rate)
        einstein = adjoint_einstein_diffusion(
            model.generator, model.stationary_state, model.operators)
        d_error = float(np.linalg.norm(einstein-model.ordered_diffusion)
                        / np.linalg.norm(model.ordered_diffusion))
        ref_atom = replace(atoms.double_lambda_rb85(gamma_gg=0), collapse_ops=tuple(
            c.operator for c in thermal_reset_channels(
                reset_rate, [5/12, 7/12, 0, 0], source="explicit reference reset fixture")))
        old = fwm.pump_only_weak_response_reference(
            pump_rabi, pump_rabi, [optical.two_photon_rad_s],
            [optical.one_photon_rad_s], atom=ref_atom)
        state_error = float(np.linalg.norm(model.stationary_state-old.pump_state[0]))
        scans = []
        for axis in axes:
            calculated = model.spectrum(axis)
            expected = qrt_ordered_spectrum(model.generator, model.stationary_state,
                                            model.operators, axis.omega_rad_s)
            errors = np.linalg.norm(calculated.ordered-expected, axis=(1, 2))/np.maximum(
                np.linalg.norm(expected, axis=(1, 2)), np.finfo(float).tiny)
            scans.append({
                "frame": axis.frame,
                "generator_omega_rad_s": axis.omega_rad_s.tolist(),
                "ordered_diffusion_vs_qrt_max_relative_error": float(np.max(errors)),
                "per_frequency_relative_error": errors.tolist(),
                "spectral_audit": calculated.audit(),
                "ordered_spectrum_seconds": _complex_array(calculated.ordered),
                "symmetrized_spectrum_seconds": _complex_array(calculated.symmetrized),
            })
        configurations.append({
            "name": name,
            "pump_rabi_rad_s": pump_rabi,
            "transit_rate_s_inverse": reset_rate,
            "stationary_state": _complex_array(model.stationary_state),
            "existing_pump_reference_state_error": state_error,
            "operators": _complex_array(model.operators),
            "atomic_drift_s_inverse": model.drift.tolist(),
            "ordered_covariance": _complex_array(model.ordered_covariance),
            "reservoir_names": list(model.reservoir_names),
            "per_reservoir_ordered_diffusion_s_inverse": _complex_array(model.channel_diffusion),
            "einstein_vs_jump_relative_error": d_error,
            "diagnostics": asdict(model.diagnostics),
            "scans": scans,
            "passed": bool(d_error < 1e-8 and state_error < 1e-10 and all(
                scan["spectral_audit"]["passed"]
                and scan["ordered_diffusion_vs_qrt_max_relative_error"] < 2e-8
                for scan in scans)),
        })
    sources = [
        "gabes/atoms.py", "gabes/core.py", "gabes/constants.py", "gabes/hyperfine.py",
        "gabes/schemes/fwm.py", "gabes/quantum/contracts.py", "gabes/quantum/reservoirs.py",
        "gabes/quantum/diffusion.py", "gabes/fwm_quantum/model.py",
        "analysis/grand_challenge/reference/atomic_qrt.py",
        "analysis/grand_challenge/atomic_noise_audit.py",
    ]
    return {
        "schema_version": 1,
        "stage": "S1 atomic diffusion and independent QRT reference",
        "convention_id": CONVENTION_ID,
        "expected_controls_passed": all(c["passed"] for c in configurations),
        "atomic_ordered_diffusion_implemented": True,
        "field_diffusion_implemented": False,
        "physical_squeezing_prediction": False,
        "experimental_validation": False,
        "environment": {"python": platform.python_version(), "numpy": np.__version__},
        "source_sha256": {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in sources},
        "operating_inputs": {
            "pump_power_W_when_on": 0.6, "pump_waist_1e2_radius_m": fwm.W_PUMP,
            "one_photon_rad_s": optical.one_photon_rad_s,
            "two_photon_rad_s": optical.two_photon_rad_s,
            "rf_analysis_hz": rf.frequency_hz.tolist(),
            "input_provenance": "repository nominal/reference inputs; no new independent apparatus measurements",
        },
        "configurations": configurations,
        "limits": [
            "Single atom, one velocity, static classical pump; exact second-order regression of the declared Markov model.",
            "The optional reset rate is a conditional input, not an experimentally validated collision/transport law.",
            "The two methods share the physical generator/state/operators but not diffusion or response implementations.",
            "Atomic operator commutators are state dependent; no canonical bosonic commutator is assigned to the atom.",
            "Density/mode normalization, field M and D, propagation, collected-mode covariance and detected S_minus remain pending.",
        ],
    }


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, help="create a new JSON report without overwriting")
    args = parser.parse_args(argv)
    report = build_report()
    body = json.dumps(report, ensure_ascii=False, indent=2, allow_nan=False)+"\n"
    if args.output is None:
        print(body, end="")
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        with args.output.open("x", encoding="utf-8") as stream:
            stream.write(body)
        print(f"Wrote {args.output}")
    return 0 if report["expected_controls_passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
