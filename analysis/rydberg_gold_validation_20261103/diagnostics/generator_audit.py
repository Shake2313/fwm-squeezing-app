"""Check stationary positivity and conditional Choi positivity of legacy model."""
from __future__ import annotations

import os
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"

import hashlib
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parents[2]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(ROOT))

from gabes import core
from gabes.schemes.rydberg import MHZ
from baseline_spectroscopy import FixedScan


def main():
    scheme = FixedScan()
    params = json.loads((REPO/"analysis/rydberg_gold_map_20261103/resolved_parameters.json").read_text(encoding="utf-8"))["baseline_parameters"]
    rows = []
    for power in (.5, 6., 10.):
        for lo in (0., 3.7):
            pr = dict(params, probe_power_uw=power, lo_rabi_mhz=lo)
            raw = scheme.compute(pr)
            gd = (raw["rydberg_dephasing_mhz"] + raw["transit_mhz"]) * MHZ
            atom = scheme._atom(gd, raw["rf_dephasing_mhz"] * MHZ)
            matrices, labels = [], []
            for delta in (-6., -3., 0., 3., 6.):
                base = core.build_liouvillian(scheme._hamiltonian_at_probe_detuning(raw, pr, delta * MHZ), atom)
                for velocity in (-600., -100., 0., 100., 600.):
                    matrices.append(base + raw["k_vec"] * velocity * atom.S_v)
                    labels.append((delta, velocity))
            rho = core.steady_state_from_liouvillian(np.array(matrices), 4)
            eig = np.linalg.eigvalsh((rho + rho.conj().transpose(0, 2, 1)) / 2)
            k = np.unravel_index(np.argmin(eig), eig.shape)[0]
            rows.append(dict(probe_power_uw=power, lo_rabi_mhz=lo,
                min_stationary_eigenvalue=float(eig.min()),
                min_case_detuning_mhz=labels[k][0], min_case_velocity_m_s=labels[k][1],
                sampled_state_count=len(labels),
                nonpositive_state_count=int(np.sum(eig.min(axis=1) < -1e-10))))

    n = 4
    choi = atom.lindblad.reshape(n, n, n, n).transpose(2, 0, 3, 1).reshape(n*n, n*n)
    identity = np.eye(n).reshape(n*n) / np.sqrt(n)
    q = np.eye(n*n) - np.outer(identity, identity)
    conditional = q @ choi @ q
    eigenvalues = np.linalg.eigvalsh((conditional + conditional.conj().T) / 2) / MHZ
    result = dict(origin="⟦AI⟧ Analysis-only generator audit; no production edits",
        criterion="A trace-preserving Hermiticity-preserving generator is GKSL iff its Choi matrix projected orthogonal to vectorized identity is positive semidefinite",
        conditional_choi_eigenvalues_per_2pi_mhz=eigenvalues.tolist(),
        conditional_choi_minimum_per_2pi_mhz=float(eigenvalues.min()),
        stationary_samples=rows,
        interpretation="Negative conditional Choi eigenvalue diagnoses phenomenological pairwise dephasing as non-GKSL. A finite stationary-state sample is not a universal positivity proof; no change to legacy numerical results is implied.",
        source_sha256={str(p.relative_to(REPO)).replace('\\','/'):hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__), ROOT/'baseline_spectroscopy.py', REPO/'analysis/rydberg_gold_map_20261103/resolved_parameters.json', REPO/'gabes/schemes/rydberg.py', REPO/'gabes/atoms.py', REPO/'gabes/core.py']})
    ROOT.mkdir(parents=True, exist_ok=True)
    (ROOT/'generator_audit.json').write_text(json.dumps(result, ensure_ascii=False, indent=2)+"\n", encoding="utf-8")
    print(json.dumps(result, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
