"""Conditional reduced microscopic field M/D and distributed-noise audit.

python -m analysis.grand_challenge.field_noise_audit --output NEW_REPORT.json
"""

import argparse
from dataclasses import replace
import hashlib
import json
from pathlib import Path
import platform

import numpy as np
import scipy

from gabes import atoms, constants as c, hyperfine, observables
from gabes.fwm_quantum.field import reduced_field_pair, reduced_readout_operators
from gabes.fwm_quantum.model import reduced_pump_noise
from gabes.quantum.contracts import CONVENTION_ID, AnalysisFrequencyAxis, OpticalDetunings
from gabes.quantum.reservoirs import thermal_reset_channels
from gabes.quantum.traveling import constant_segment, compose_segments
from gabes.schemes import fwm
from .reference.atomic_qrt import qrt_ordered_spectrum


ROOT = Path(__file__).resolve().parents[2]


def _complex(value):
    value = np.asarray(value)
    return {"real": value.real.tolist(), "imag": value.imag.tolist()}


def _error(actual, expected):
    return float(np.max(np.linalg.norm(actual-expected, axis=(-2, -1))/np.maximum(
        np.linalg.norm(expected, axis=(-2, -1)), np.finfo(float).tiny)))


def build_report():
    optical = OpticalDetunings(2*np.pi*.9e9, -2*np.pi*8e6)
    axis = AnalysisFrequencyAxis.from_hz(np.arange(-40, 41)*1e5)
    # Use the repository's declared optical carrier convention for parity.
    pump_carrier = c.OMEGA_D1+optical.one_photon_rad_s
    beat = -c.OMEGA_HF+optical.two_photon_rad_s
    carriers = [pump_carrier+beat, pump_carrier-beat]
    settings = dict(number_density_m3=1e18, uniform_area_m2=1.2e-7,
                    optical_omega_rad_s=carriers,
                    effective_dipole_C_m=c.DIPOLE_D1/np.sqrt(hyperfine.N_GROUND_SUBLEVELS))
    length = .0125
    configurations = []
    for name, pump, reset in (
            ("pump_off_with_declared_reset", 0., c.GAMMA_GG),
            ("pumped_radiative_only", fwm.rabi_freq(.6, fwm.W_PUMP), 0.),
            ("pumped_with_declared_reset", fwm.rabi_freq(.6, fwm.W_PUMP), c.GAMMA_GG)):
        atom = reduced_pump_noise(pump, optical.one_photon_rad_s, transit_rate_s_inverse=reset)
        pair = reduced_field_pair(atom, optical, axis, **settings)
        explicit = replace(atoms.double_lambda_rb85(gamma_gg=0), collapse_ops=tuple(
            channel.operator for channel in thermal_reset_channels(
                reset, [5/12, 7/12, 0, 0], source="identical explicit reference fixture")))
        old = fwm.pump_only_weak_response_reference(
            pump, pump, [optical.two_photon_rad_s], [optical.one_photon_rad_s], atom=explicit,
            analysis_frequency_axis_rad_s=axis.omega_rad_s)
        chi = [v[:, 0, 0] for v in (old.chi_ss, old.chi_sc, old.chi_cs, old.chi_cc)]
        _, _, t_field = observables.gain_from_chi(
            *chi, carriers[0]/c.C_LIGHT, carriers[1]/c.C_LIGHT, length,
            settings["number_density_m3"], dipole=c.DIPOLE_D1,
            line_strength=1/hyperfine.N_GROUND_SUBLEVELS)
        expected_t, _ = observables.canonical_transfer_from_field(
            t_field, *carriers, settings["uniform_area_m2"], settings["uniform_area_m2"])
        q = np.diag(observables.photon_flux_mode_matrix(
            *carriers, settings["uniform_area_m2"], settings["uniform_area_m2"]))
        g = settings["effective_dipole_C_m"]*q/(2*c.HBAR)
        w = np.einsum("jab,iba->ji", reduced_readout_operators(), atom.operators)
        radiation = -1j*np.diag([g[0], -g[1]])@w
        omega = pair.main.frequency_axis.omega_rad_s
        references = (
            qrt_ordered_spectrum(atom.generator, atom.stationary_state, atom.operators, omega),
            qrt_ordered_spectrum(atom.generator, atom.stationary_state, atom.operators,
                                 -omega).transpose(0, 2, 1))
        qrt_errors = []
        for source, actual in zip(references, (pair.main.noise_greater, pair.main.noise_lesser)):
            expected = settings["number_density_m3"]*settings["uniform_area_m2"]*(
                radiation@source@radiation.conj().T)
            qrt_errors.append(_error(actual, expected))
        main_transfer = constant_segment(pair.main, length)
        split = compose_segments(constant_segment(pair.main, length*.37),
                                  constant_segment(pair.main, length*.63))
        parity_error = _error(main_transfer.transfer, expected_t)
        split_error = max(_error(getattr(main_transfer, n), getattr(split, n))
                          for n in ("transfer", "noise_greater", "noise_lesser"))
        companion_error = max(
            _error(pair.companion.drift, pair.main.drift[::-1].conj()),
            _error(pair.companion.noise_greater, pair.main.noise_lesser[::-1].conj()),
            _error(pair.companion.noise_lesser, pair.main.noise_greater[::-1].conj()))
        bad = replace(pair.main,
                      noise_greater_by_reservoir=pair.main.noise_greater_by_reservoir/12,
                      noise_lesser_by_reservoir=pair.main.noise_lesser_by_reservoir/12)
        sectors = []
        for local in (pair.main, pair.companion):
            propagated = constant_segment(local, length)
            greater, lesser = propagated.vacuum_output()
            sectors.append({
                "frame": local.frequency_axis.frame,
                "generator_omega_rad_s": local.frequency_axis.omega_rad_s.tolist(),
                "mode_labels": local.mode_labels, "nambu_signs": local.signs.tolist(),
                "drift_m_inverse": _complex(local.drift),
                "reservoir_names": local.reservoir_names,
                "noise_greater_by_reservoir_m_inverse": _complex(local.noise_greater_by_reservoir),
                "noise_lesser_by_reservoir_m_inverse": _complex(local.noise_lesser_by_reservoir),
                "local_audit": local.audit(),
                "transfer": _complex(propagated.transfer),
                "added_noise_greater": _complex(propagated.noise_greater),
                "added_noise_lesser": _complex(propagated.noise_lesser),
                "propagation_audit": propagated.audit(),
                "vacuum_output_greater": _complex(greater),
                "vacuum_output_lesser": _complex(lesser),
                "minimum_vacuum_output_greater_eigenvalue": float(np.linalg.eigvalsh(greater).min()),
                "minimum_vacuum_output_lesser_eigenvalue": float(np.linalg.eigvalsh(lesser).min()),
            })
        configurations.append({
            "name": name, "pump_rabi_rad_s": pump, "reset_rate_s_inverse": reset,
            "atomic_stationary_state": _complex(atom.stationary_state),
            "existing_maxwell_transfer_relative_error": parity_error,
            "field_greater_qrt_relative_error": qrt_errors[0],
            "field_lesser_qrt_relative_error": qrt_errors[1],
            "companion_relative_error": companion_error,
            "constant_segment_split_relative_error": split_error,
            "noise_only_factor_1_over_12_rejected": not bad.audit()["passed"],
            "sectors": sectors,
            "passed": bool(max(parity_error, *qrt_errors, companion_error, split_error) < 2e-8
                           and not bad.audit()["passed"]
                           and all(s["local_audit"]["passed"] and s["propagation_audit"]["passed"]
                                   and s["minimum_vacuum_output_greater_eigenvalue"] > -1e-10
                                   and s["minimum_vacuum_output_lesser_eigenvalue"] > -1e-10
                                   for s in sectors)),
        })
    sources = [
        "gabes/atoms.py", "gabes/core.py", "gabes/constants.py", "gabes/hyperfine.py",
        "gabes/observables.py", "gabes/schemes/fwm.py", "gabes/quantum/contracts.py",
        "gabes/quantum/reservoirs.py", "gabes/quantum/diffusion.py", "gabes/quantum/traveling.py",
        "gabes/fwm_quantum/model.py", "gabes/fwm_quantum/field.py",
        "analysis/grand_challenge/reference/atomic_qrt.py",
        "analysis/grand_challenge/field_noise_audit.py", "requirements-quantum.txt",
    ]
    return {
        "schema_version": 1, "convention_id": CONVENTION_ID,
        "stage": "S1 conditional local field M/D and prescribed-segment propagation",
        "expected_controls_passed": all(v["passed"] for v in configurations),
        "reduced_field_diffusion_implemented": True,
        "physical_squeezing_prediction": False, "experimental_validation": False,
        "environment": {"python": platform.python_version(), "numpy": np.__version__, "scipy": scipy.__version__},
        "source_sha256": {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in sources},
        "operating_inputs": {**settings, "one_photon_rad_s": optical.one_photon_rad_s,
            "two_photon_rad_s": optical.two_photon_rad_s, "analysis_rf_hz": axis.frequency_hz.tolist(),
            "length_m": length, "phase_mismatch_rad_m": 0.,
            "input_provenance": "conditional numerical fixtures, not independent measurements",
            "effective_dipole_source": "DIPOLE_D1/sqrt(N_GROUND_SUBLEVELS); declared reduced 1/12 Maxwell normalization, applied reciprocally to weak-field drive and emission",
            "pump_rabi_source": "legacy rabi_freq(0.6 W, W_PUMP), held as explicit conditional Rabi; absolute pump-power/dipole normalization remains unaudited"},
        "configurations": configurations,
        "limits": [
            "Selected double-Lambda transitions and phase-selected two-mode sectors; additional optical modes/cross-sector processes are outside this closure.",
            "Single velocity, independent atoms, one common uniform transverse area and constant classical pump; no collision, angular-Doppler or transverse-mode convergence is claimed.",
            "Slice normalization cancels the common area at fixed density and pump Rabi; this is not invariance under changing a real Gaussian beam waist at fixed power.",
            "Greater and lesser noise are independently projected from atomic jump diffusion, not inferred from M or repaired to pass a commutator check.",
            "Segment propagation integrates distributed source noise; segment states/pumps are prescribed and depletion is not solved.",
            "Finite-temporal-mode/quadrature readout, SQL, collected intensity-difference spectrum and independent apparatus input ledger remain pending.",
        ],
    }


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, help="create a new report without overwriting")
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
