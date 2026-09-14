"""Run the S0 foundations: python -m analysis.grand_challenge.foundation_audit.

Use --output PATH to create a new report. Existing reports are never replaced.
This command audits foundations; it does not predict atomic squeezing.
"""

import argparse
from dataclasses import asdict
import hashlib
import json
from pathlib import Path
import platform

import numpy as np

from gabes import atoms, constants
from gabes.quantum.channels import GaussianChannel, compose_channels
from gabes.quantum.contracts import CONVENTION_ID, AnalysisFrequencyAxis
from gabes.quantum.reservoirs import (
    ExplicitReservoirs, audit_generator, thermal_reset_channels,
)
from gabes.schemes.fwm import collisional_atom, thermal_transit_reset_superoperator


ROOT = Path(__file__).resolve().parents[2]


def build_report():
    atom = atoms.double_lambda_rb85(gamma_gg=0)
    radiative = ExplicitReservoirs.from_atom(
        atom, source="gabes.atoms.double_lambda_rb85: CF2 branching and natural GAMMA")
    reset = ExplicitReservoirs(4, thermal_reset_channels(
        constants.GAMMA_GG, [5/12, 7/12, 0, 0],
        source="declared Markov thermal replacement; nominal rate, not measured transport"))
    reset_reference = thermal_transit_reset_superoperator(constants.GAMMA_GG)
    reset_error = float(np.linalg.norm(reset.dissipator()-reset_reference)
                        / np.linalg.norm(reset_reference))
    combined = ExplicitReservoirs(4, radiative.channels+reset.channels)
    generator_results = {
        "explicit_radiative": asdict(audit_generator(radiative.dissipator())),
        "explicit_thermal_reset": asdict(audit_generator(reset.dissipator())),
        "explicit_radiative_plus_reset": asdict(audit_generator(combined.dissipator())),
        "legacy_default_counterexample": asdict(audit_generator(
            atoms.double_lambda_rb85().lindblad)),
        "legacy_combined_121C": asdict(audit_generator(collisional_atom(394.15).lindblad)),
    }
    modes = ("probe", "conjugate")
    loss = GaussianChannel.vacuum_attenuator(modes, [0.71, 0.86], source="analytic fixture")
    channel = compose_channels(GaussianChannel.identity(modes), loss)
    vacuum_error = float(np.max(np.abs(channel.apply_covariance(np.eye(4)/2)-np.eye(4)/2)))
    expected = all(generator_results[k]["passed"] for k in (
        "explicit_radiative", "explicit_thermal_reset", "explicit_radiative_plus_reset",
        "legacy_combined_121C"))
    expected = (expected and not generator_results["legacy_default_counterexample"]["passed"]
                and reset_error < 1e-12 and vacuum_error < 1e-12 and channel.audit().passed)
    paths = [
        "gabes/atoms.py", "gabes/core.py", "gabes/constants.py", "gabes/hyperfine.py",
        "gabes/schemes/fwm.py", "gabes/quantum/contracts.py",
        "gabes/quantum/reservoirs.py", "gabes/quantum/channels.py",
        "analysis/grand_challenge/foundation_audit.py",
    ]
    return {
        "schema_version": 1,
        "stage": "S0 foundations, initial implementation",
        "convention_id": CONVENTION_ID,
        "expected_controls_passed": bool(expected),
        "physical_squeezing_prediction": False,
        "microscopic_diffusion_implemented": False,
        "experimental_validation": False,
        "environment": {"python": platform.python_version(), "numpy": np.__version__},
        "source_sha256": {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in paths},
        "generator_audits": generator_results,
        "reset_assembly_relative_error": reset_error,
        "passive_channel_audit": asdict(channel.audit()),
        "passive_vacuum_max_error": vacuum_error,
        "analysis_axis_example_rad_s": AnalysisFrequencyAxis.from_hz(
            [-1e6, 0, 1e6]).omega_rad_s.tolist(),
        "limits": [
            "A failed legacy-default audit is the expected negative control, not repaired physics.",
            "The combined 121 C generator passing does not identify its independent reservoirs.",
            "The explicit thermal reset uses a nominal rate solely to check algebraic equivalence.",
            "No diffusion, atomic field propagation, RF spectrum or measured no-fit prediction is supplied.",
        ],
    }


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, help="create a new JSON report (no overwrite)")
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
