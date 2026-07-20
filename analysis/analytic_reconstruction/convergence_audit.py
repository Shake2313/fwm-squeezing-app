"""Numerical convergence audit for the reduced FWM reconstruction.

This analysis-only script deliberately leaves the production solver untouched.
It extends the exact equations in ``ref_solver.py`` from the hard-coded
``n=-1,0,+1`` Floquet truncation to ``n=-N_F,...,+N_F`` using block continued
fractions.  The recurrence is algebraically identical to the existing solver at
``N_F=1`` and is checked against a direct finite-block solve for higher orders.

The script also refines the *existing one-dimensional* Maxwell velocity
quadrature.  That test does not add the non-collinear two-photon Doppler shift;
the generated report calls out this limitation explicitly.  Gain values inherit
``ref_solver.py``'s archived dressed-k/refractive-mismatch propagation and are
therefore isolation tests, not corrected Option-A predictions.

Separately, the script evaluates the repository's literature benchmark with the
corrected no-double-count Option-A propagation convention (bare optical wave
numbers plus vacuum/geometric mismatch).  It also constructs one minimum-vacuum
Gaussian commutator completion.  That completion is deliberately labelled a
mathematical dilation, not microscopic atomic diffusion or a squeezing spectrum.

Run from the repository root::

    python analysis/analytic_reconstruction/convergence_audit.py

Outputs are written below ``analysis/analytic_reconstruction/generated/``.
"""

from __future__ import annotations

import argparse
from contextlib import contextmanager
import json
import math
import os
from pathlib import Path
import sys
import tempfile
import time

import numpy as np


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import ref_solver as reference  # noqa: E402
from gabes import constants as gabes_constants, doppler as gabes_doppler  # noqa: E402
from gabes.core import (  # noqa: E402
    blas_single_thread,
    build_liouvillian,
    matrix_exp_2x2,
)
from gabes.schemes import fwm as fwm_scheme  # noqa: E402


MODEL_ARGS = {
    "D_GHz": -1.5,
    "T": 383.15,
    "P_pump": 0.6,
    "P_probe": 8e-6,
    "ls": 0.74,
    "loss": 0.055,
    "qe": 0.92,
    "L": 12.5e-3,
    "wp": 530e-6,
    "ws": 330e-6,
    "theta_deg": 0.32,
    "kappa": 0.1,
}
FIXED_DELTA_MHZ = -280.0

LITERATURE_ARGS = {
    "D_GHz": 0.9,
    "delta_MHz": -8.0,
    "T": 394.15,
    "P_pump": 0.6,
    "P_probe": 8e-6,
    "line_strength_residual": 0.74,
    "L": 12.5e-3,
    "w_pump": 530e-6,
    "w_probe": 330e-6,
    "theta_deg": 0.32,
    "branch": -1,
    "velocity_step_m_per_s": 1.0,
    "velocity_cutoff_sigma": 5.0,
    "external_qe": 0.92,
    "external_path_loss": 0.055,
}


def floquet_solve_truncated(
    L0: np.ndarray,
    Cp: np.ndarray,
    Cm: np.ndarray,
    omega_beat: float,
    deff_axis: np.ndarray,
    n_f: int,
    *,
    S_v: np.ndarray = reference.S_V,
    n_levels: int = reference.NL,
) -> tuple[np.ndarray, np.ndarray]:
    """Return rho_0 and rho_+1 for a finite ``[-n_f,+n_f]`` Floquet block.

    The harmonic equations are

        (L0 + i*n*omega_beat) rho_n + Cp rho_{n-1} + Cm rho_{n+1} = 0.

    Eliminating the positive and negative chains gives an exact continued
    fraction at the chosen finite boundary.  Only 16x16 batched solves are
    required, rather than a dense ``(2*n_f+1)*16`` solve at every velocity.
    """
    if n_f < 1:
        raise ValueError("n_f must be at least 1")

    deff_axis = np.asarray(deff_axis, dtype=float)
    n_deff = deff_axis.size
    M = L0.shape[0]
    eye = np.eye(M, dtype=complex)
    L_batch = L0[None, :, :] - deff_axis[:, None, None] * S_v[None, :, :]
    Cp_batch = np.broadcast_to(Cp, (n_deff, M, M))
    Cm_batch = np.broadcast_to(Cm, (n_deff, M, M))

    # rho_n = R_n rho_{n-1}, n > 0; eliminate from +N_F toward +1.
    R = None
    for harmonic in range(n_f, 0, -1):
        A = L_batch + 1j * harmonic * omega_beat * eye[None, :, :]
        if R is not None:
            A = A + Cm_batch @ R
        R = -np.linalg.solve(A, Cp_batch)

    # rho_n = Q_n rho_{n+1}, n < 0; eliminate from -N_F toward -1.
    Q = None
    for harmonic in range(-n_f, 0):
        A = L_batch + 1j * harmonic * omega_beat * eye[None, :, :]
        if Q is not None:
            A = A + Cp_batch @ Q
        Q = -np.linalg.solve(A, Cm_batch)

    A_eff = L_batch + Cp_batch @ Q + Cm_batch @ R
    A_eff[:, 0, :] = 0.0
    for state in range(n_levels):
        A_eff[:, 0, state * n_levels + state] = 1.0
    rhs = np.zeros((n_deff, M, 1), dtype=complex)
    rhs[:, 0, 0] = 1.0
    rho0_vec = np.linalg.solve(A_eff, rhs)
    rho1_vec = R @ rho0_vec
    shape = (n_deff, n_levels, n_levels)
    return rho0_vec[:, :, 0].reshape(shape), rho1_vec[:, :, 0].reshape(shape)


def _direct_floquet_one(
    L0: np.ndarray,
    Cp: np.ndarray,
    Cm: np.ndarray,
    omega_beat: float,
    deff: float,
    n_f: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Dense finite-block reference used only to validate the recurrence."""
    M = L0.shape[0]
    sectors = np.arange(-n_f, n_f + 1)
    K = sectors.size
    block = np.zeros((K * M, K * M), dtype=complex)
    L = L0 - deff * reference.S_V
    eye = np.eye(M, dtype=complex)
    for q, harmonic in enumerate(sectors):
        sl = slice(q * M, (q + 1) * M)
        block[sl, sl] = L + 1j * harmonic * omega_beat * eye
        if q > 0:
            block[sl, slice((q - 1) * M, q * M)] = Cp
        if q + 1 < K:
            block[sl, slice((q + 1) * M, (q + 2) * M)] = Cm

    rhs = np.zeros(K * M, dtype=complex)
    trace_indices = np.arange(reference.NL) * reference.NL + np.arange(reference.NL)
    for q, harmonic in enumerate(sectors):
        row = q * M
        block[row, :] = 0.0
        block[row, q * M + trace_indices] = 1.0
        rhs[row] = 1.0 if harmonic == 0 else 0.0
    solution = np.linalg.solve(block, rhs).reshape(K, M)
    q0 = n_f
    return solution[q0], solution[q0 + 1]


def solver_self_checks() -> dict[str, float]:
    """Check N_F=1 parity and recurrence/direct-block parity."""
    delta = 2 * np.pi * FIXED_DELTA_MHZ * 1e6
    Op = reference.rabi(MODEL_ARGS["P_pump"], MODEL_ARGS["wp"])
    Os = reference.rabi(MODEL_ARGS["P_probe"], MODEL_ARGS["ws"])
    density = reference.density(MODEL_ARGS["T"])
    vbar = np.sqrt(
        8 * reference.KB * MODEL_ARGS["T"] / (np.pi * reference.MASS / 2)
    )
    gamma_gg = reference.GAMMA_GG0 + density * reference.XSEC * vbar
    gamma_opt = 0.5 * reference.BETA_SELF * density
    L0 = (
        reference.comm_super(reference.H0_build(Op, Os, delta))
        + reference.lindblad(gamma_gg, gamma_opt)
    )
    Hp = reference.Hp_build(Op, 0.0)
    Cp = reference.comm_super(Hp)
    Cm = reference.comm_super(Hp.conj().T)
    omega_beat = reference.OMEGA_HF - delta
    deff = np.linspace(2 * np.pi * -3.0e9, 2 * np.pi * -0.1e9, 7)

    rho0_old, rho1_old = reference.floquet(L0, Cp, Cm, omega_beat, deff)
    rho0_new, rho1_new = floquet_solve_truncated(
        L0, Cp, Cm, omega_beat, deff, 1
    )
    checks = {
        "n_f_1_rho0_max_abs_difference": float(np.max(np.abs(rho0_old - rho0_new))),
        "n_f_1_rho1_max_abs_difference": float(np.max(np.abs(rho1_old - rho1_new))),
    }

    for n_f in (2, 3):
        r0, r1 = floquet_solve_truncated(
            L0, Cp, Cm, omega_beat, np.array([deff[3]]), n_f
        )
        d0, d1 = _direct_floquet_one(L0, Cp, Cm, omega_beat, deff[3], n_f)
        checks[f"n_f_{n_f}_direct_block_rho0_max_abs_difference"] = float(
            np.max(np.abs(r0[0].reshape(-1) - d0))
        )
        checks[f"n_f_{n_f}_direct_block_rho1_max_abs_difference"] = float(
            np.max(np.abs(r1[0].reshape(-1) - d1))
        )
    return checks


def _complex_payload(value: complex) -> dict[str, float]:
    value = complex(value)
    return {"real": float(value.real), "imag": float(value.imag)}


def _complex_matrix_payload(matrix: np.ndarray) -> list[list[dict[str, float]]]:
    matrix = np.asarray(matrix, dtype=complex)
    return [[_complex_payload(value) for value in row] for row in matrix]


def _literature_context(probe_power: float | None = None) -> dict:
    """Build shared atomic/velocity data for the corrected literature point."""
    p = LITERATURE_ARGS
    probe_power = p["P_probe"] if probe_power is None else float(probe_power)
    branch = int(p["branch"])
    delta = 2 * np.pi * p["delta_MHz"] * 1e6
    Delta = 2 * np.pi * p["D_GHz"] * 1e9
    Op = gabes_constants.rabi_freq(p["P_pump"], p["w_pump"])
    Os = gabes_constants.rabi_freq(probe_power, p["w_probe"])
    density = fwm_scheme.hyperfine.number_density(p["T"])
    atom = fwm_scheme.collisional_atom(p["T"], density)
    velocity, weights = gabes_doppler.velocity_grid(
        p["T"],
        dv=p["velocity_step_m_per_s"],
        cutoff_sigma=p["velocity_cutoff_sigma"],
    )
    deff = gabes_doppler.build_Delta_eff_axis(Delta, Delta, velocity)
    return {
        "branch": branch,
        "delta": delta,
        "Delta": Delta,
        "Op": Op,
        "Os": Os,
        "P_probe": probe_power,
        "density": density,
        "atom": atom,
        "velocity": velocity,
        "weights": weights,
        "deff": deff,
        "omega_beat": fwm_scheme.seeded_sideband_beat(delta, branch),
    }


def _literature_reduced_chi(context: dict, n_f: int) -> tuple[complex, ...]:
    """Doppler-averaged (ss, sc, cs, cc) reduced responses at one N_F."""
    branch = context["branch"]
    Op = context["Op"]
    Os = context["Os"]
    delta = context["delta"]
    atom = context["atom"]
    deff = context["deff"]
    omega_beat = context["omega_beat"]
    probe_ground = fwm_scheme.G2 if branch == -1 else fwm_scheme.G1
    conj_ground = fwm_scheme.G1 if branch == -1 else fwm_scheme.G2

    # Probe-seeded response: rho_0 -> chi_ss, rho_+1 -> chi_cs.
    Cp, Cm = fwm_scheme.sideband_template(Op, Op, 0.0, branch)
    L0 = build_liouvillian(
        fwm_scheme.static_hamiltonian_at_Deff_zero(Op, Op, Os, delta, branch),
        atom,
    )
    rho0, rho1 = floquet_solve_truncated(
        L0,
        Cp,
        Cm,
        omega_beat,
        deff,
        n_f,
        S_v=atom.S_v,
        n_levels=fwm_scheme.N_LEVELS,
    )
    chi_ss = fwm_scheme._polarization_coherence(rho0, probe_ground) / Os
    chi_cs = fwm_scheme._polarization_coherence(rho1, conj_ground) / Os

    # Conjugate-seeded response: rho_0 -> chi_sc, rho_+1 -> chi_cc.
    Cp, Cm = fwm_scheme.sideband_template(Op, Op, Os, branch)
    L0 = build_liouvillian(
        fwm_scheme.static_hamiltonian_at_Deff_zero(Op, Op, 0.0, delta, branch),
        atom,
    )
    rho0, rho1 = floquet_solve_truncated(
        L0,
        Cp,
        Cm,
        omega_beat,
        deff,
        n_f,
        S_v=atom.S_v,
        n_levels=fwm_scheme.N_LEVELS,
    )
    chi_sc = fwm_scheme._polarization_coherence(rho0, probe_ground) / Os
    chi_cc = fwm_scheme._polarization_coherence(rho1, conj_ground) / Os

    def average(table: np.ndarray) -> complex:
        return gabes_doppler.doppler_average(
            table[None, :],
            deff,
            context["Delta"],
            context["velocity"],
            context["weights"],
        )[0]

    return average(chi_ss), average(chi_sc), average(chi_cs), average(chi_cc)


def _option_a_transfer(context: dict, n_f: int) -> dict:
    """Corrected bare-k/geometric-mismatch drift and transfer for one N_F."""
    p = LITERATURE_ARGS
    chi_ss, chi_sc, chi_cs, chi_cc = _literature_reduced_chi(context, n_f)
    effective_line_strength = (
        p["line_strength_residual"]
        * fwm_scheme.physical_coupling_norm(context["branch"])
    )
    coupling = (
        -2.0
        * context["density"]
        * effective_line_strength
        * gabes_constants.DIPOLE_D1**2
        / (gabes_constants.EPS_0 * gabes_constants.HBAR)
    )

    probe_GHz = (
        fwm_scheme.branch_center_GHz(p["D_GHz"], context["branch"])
        + p["delta_MHz"] * 1e-3
    )
    pump_offset = 2 * np.pi * p["D_GHz"] * 1e9
    probe_offset = 2 * np.pi * probe_GHz * 1e9
    conjugate_offset = 2 * pump_offset - probe_offset
    k_pump = (gabes_constants.OMEGA_D1 + pump_offset) / gabes_constants.C_LIGHT
    k_probe = (gabes_constants.OMEGA_D1 + probe_offset) / gabes_constants.C_LIGHT
    k_conj = (gabes_constants.OMEGA_D1 + conjugate_offset) / gabes_constants.C_LIGHT
    delta_k_vac = 2 * k_pump - (k_probe + k_conj) * np.cos(
        np.radians(p["theta_deg"])
    )

    # Symmetric rotating frame. Re(chi) occurs only in the diagonal response;
    # delta_k_vac contains no refractive-index contribution.
    M = np.array(
        [
            [
                0.5j * k_probe * coupling * chi_ss + 0.5j * delta_k_vac,
                0.5j * k_probe * coupling * chi_sc,
            ],
            [
                -0.5j * k_conj * coupling * np.conj(chi_cs),
                -0.5j * k_conj * coupling * np.conj(chi_cc)
                - 0.5j * delta_k_vac,
            ],
        ],
        dtype=complex,
    )
    transfer = matrix_exp_2x2(M[None, :, :], p["L"])[0]
    G_s = float(abs(transfer[0, 0]) ** 2)
    G_c = float(abs(transfer[1, 0]) ** 2)
    return {
        "N_F": int(n_f),
        "reduced_chi_seconds": {
            "ss": _complex_payload(chi_ss),
            "sc": _complex_payload(chi_sc),
            "cs": _complex_payload(chi_cs),
            "cc": _complex_payload(chi_cc),
        },
        "effective_line_strength": float(effective_line_strength),
        "probe_offset_GHz": float(probe_GHz),
        "conjugate_offset_GHz": float(2 * p["D_GHz"] - probe_GHz),
        "k_pump_per_m": float(k_pump),
        "k_probe_per_m": float(k_probe),
        "k_conjugate_per_m": float(k_conj),
        "delta_k_vac_per_m": float(delta_k_vac),
        "M_per_m": _complex_matrix_payload(M),
        "T": _complex_matrix_payload(transfer),
        "G_s": G_s,
        "G_c": G_c,
        "gain_gap": G_s - G_c,
        "probe_power_uW": float(1e6 * context["P_probe"]),
        "probe_rabi_over_2pi_MHz": float(context["Os"] / (2 * np.pi * 1e6)),
        "arg_chi_sc_rad": float(np.angle(chi_sc)),
        "arg_chi_cs_rad": float(np.angle(chi_cs)),
        "arg_chi_sc_deg": float(np.degrees(np.angle(chi_sc))),
        "arg_chi_cs_deg": float(np.degrees(np.angle(chi_cs))),
        "evaluation_type": "fixed literature operating point",
        "delta_mhz": float(p["delta_MHz"]),
        "_M_array": M,
        "_T_array": transfer,
    }


def _constant_matrix_integral(
    M: np.ndarray, source: np.ndarray, length: float, points: int
) -> np.ndarray:
    """Gauss-Legendre integral of exp(Mu) source exp(M^dagger u)."""
    nodes, weights = np.polynomial.legendre.leggauss(points)
    u = 0.5 * length * (nodes + 1.0)
    weights = 0.5 * length * weights
    propagators = matrix_exp_2x2(
        np.broadcast_to(M, (points, 2, 2)), u
    )
    return np.einsum(
        "n,nij,jk,nlk->il",
        weights,
        propagators,
        source,
        propagators.conj(),
    )


def _bright_seed_noise(
    covariance: np.ndarray, transfer: np.ndarray, weight: float
) -> dict[str, float]:
    """Linearized bright-seed intensity difference for z=(a_s,a_c^dagger)."""
    beta_s = transfer[0, 0]
    beta_c = np.conj(transfer[1, 0])
    G_s = float(abs(beta_s) ** 2)
    G_c = float(abs(beta_c) ** 2)
    vector = np.array([np.conj(beta_s), -weight * beta_c], dtype=complex)
    variance = float(2 * np.real(vector @ covariance @ np.conj(vector)))
    shot_noise = float(G_s + weight**2 * G_c)
    S = variance / shot_noise
    return {
        "electronic_weight": float(weight),
        "variance": variance,
        "shot_noise_reference": shot_noise,
        "S_linear": float(S),
        "S_dB": float(10 * np.log10(S)),
    }


def option_a_literature_diagnostic() -> dict:
    """Reproducible Option-A transfer plus a minimum-vacuum dilation diagnostic.

    The dilation enforces the canonical commutator for the supplied static M but
    does not derive atomic Langevin diffusion. It must not be interpreted as a
    microscopic or frequency-dependent squeezing calculation.
    """
    p = LITERATURE_ARGS
    context = _literature_context()
    with blas_single_thread():
        convergence = [_option_a_transfer(context, n_f) for n_f in (1, 2, 3)]
        seed_linearity = []
        for probe_uW in (2.0, 4.0, 8.0):
            if probe_uW == 8.0:
                response = convergence[-1]
            else:
                response = _option_a_transfer(
                    _literature_context(probe_uW * 1e-6), 3
                )
            seed_linearity.append(
                {
                    key: response[key]
                    for key in (
                        "probe_power_uW",
                        "probe_rabi_over_2pi_MHz",
                        "G_s",
                        "G_c",
                        "gain_gap",
                        "arg_chi_sc_deg",
                        "arg_chi_cs_deg",
                    )
                }
            )
    weak_reference = seed_linearity[0]
    for row in seed_linearity:
        row["G_s_change_from_2uW_pct"] = float(
            100 * (row["G_s"] - weak_reference["G_s"]) / weak_reference["G_s"]
        )
        row["G_c_change_from_2uW_pct"] = float(
            100 * (row["G_c"] - weak_reference["G_c"]) / weak_reference["G_c"]
        )
        row["gain_gap_change_from_2uW_pct"] = float(
            100
            * (row["gain_gap"] - weak_reference["gain_gap"])
            / weak_reference["gain_gap"]
        )
    final = convergence[-1]
    M = final.pop("_M_array")
    transfer = final.pop("_T_array")
    for row in convergence[:-1]:
        row.pop("_M_array")
        row.pop("_T_array")

    J = np.diag([1.0, -1.0]).astype(complex)
    K = -(M @ J + J @ M.conj().T)
    eigenvalues, eigenvectors = np.linalg.eigh(K)
    B = eigenvectors @ np.diag(np.sqrt(np.abs(eigenvalues)))
    J_f = np.diag(np.sign(eigenvalues)).astype(complex)
    factorized_K = B @ J_f @ B.conj().T
    D_min = 0.5 * B @ B.conj().T

    quadrature_points = 400
    commutator_integral = _constant_matrix_integral(
        M, K, p["L"], quadrature_points
    )
    diffusion_integral = _constant_matrix_integral(
        M, D_min, p["L"], quadrature_points
    )
    # Independent lower-order quadrature check.
    commutator_integral_200 = _constant_matrix_integral(M, K, p["L"], 200)
    diffusion_integral_200 = _constant_matrix_integral(M, D_min, p["L"], 200)

    bare_commutator = transfer @ J @ transfer.conj().T
    completed_commutator = bare_commutator + commutator_integral
    V_in = 0.5 * np.eye(2, dtype=complex)
    V_transfer_only = transfer @ V_in @ transfer.conj().T
    V_out = V_transfer_only + diffusion_integral

    G_s = final["G_s"]
    G_c = final["G_c"]
    unweighted = _bright_seed_noise(V_out, transfer, 1.0)
    dc_balanced = _bright_seed_noise(V_out, transfer, G_s / G_c)
    transfer_only = _bright_seed_noise(V_transfer_only, transfer, 1.0)
    ideal = 1.0 / (2.0 * G_s - 1.0)
    eta_ext = p["external_qe"] * (1.0 - p["external_path_loss"])
    detected = 1.0 - eta_ext + eta_ext * unweighted["S_linear"]
    ideal_detected = 1.0 - eta_ext + eta_ext * ideal

    literature_gain = 15.5
    literature_squeezing_db = -7.8
    literature_bandwidth_mhz = 3.5
    return {
        "classification": {
            "mean_field_propagation": "corrected Option A",
            "dispersion_counting": "bare k in M; vacuum/geometric mismatch only",
            "dilation": "minimum-vacuum mathematical Gaussian completion",
            "microscopic_atomic_diffusion": False,
            "frequency_dependent": False,
            "physical_squeezing_prediction": False,
            "weak_field_reference_linearity_tested": True,
            "weak_field_reference_test_range_uW": [2.0, 4.0, 8.0],
            "canonical_mode_normalization": (
                "assumes the classical two-component amplitudes are identified "
                "with photon-flux-normalized canonical modes; production code "
                "does not independently derive this scaling"
            ),
        },
        "parameters": {
            **p,
            "temperature_C": p["T"] - 273.15,
            "n_velocity": int(context["velocity"].size),
            "n_delta_eff": int(context["deff"].size),
            "atomic_density_per_m3": float(context["density"]),
            "sigma_v_m_per_s": float(
                np.sqrt(gabes_constants.KB * p["T"] / gabes_constants.MASS_85RB)
            ),
            "angular_two_photon_doppler_included": False,
        },
        "floquet_convergence": convergence,
        "seed_reference_linearity": seed_linearity,
        "option_a_N_F_3": final,
        "minimum_vacuum_mathematical_dilation": {
            "units": {
                "J": "dimensionless",
                "K": "m^-1",
                "B": "m^-1/2",
                "J_f": "dimensionless",
                "D_min": "m^-1",
                "commutator_integral": "dimensionless",
                "V_out": "dimensionless",
            },
            "J": _complex_matrix_payload(J),
            "K_per_m": _complex_matrix_payload(K),
            "K_eigenvalues_per_m": [float(value) for value in eigenvalues],
            "B_per_sqrt_m": _complex_matrix_payload(B),
            "J_f": _complex_matrix_payload(J_f),
            "D_min_per_m": _complex_matrix_payload(D_min),
            "commutator_integral": _complex_matrix_payload(commutator_integral),
            "V_out": _complex_matrix_payload(V_out),
            "quadrature_points": quadrature_points,
            "factorization_residual_max": float(np.max(np.abs(factorized_K - K))),
            "bare_commutator_residual_max": float(
                np.max(np.abs(bare_commutator - J))
            ),
            "completed_commutator_residual_max": float(
                np.max(np.abs(completed_commutator - J))
            ),
            "completed_commutator_residual_relative_frobenius": float(
                np.linalg.norm(completed_commutator - J) / np.linalg.norm(J)
            ),
            "quadrature_200_to_400_commutator_change_max": float(
                np.max(np.abs(commutator_integral - commutator_integral_200))
            ),
            "quadrature_200_to_400_covariance_change_max": float(
                np.max(np.abs(diffusion_integral - diffusion_integral_200))
            ),
            "transfer_only_invalid_unweighted": transfer_only,
            "completed_unweighted": unweighted,
            "completed_dc_balanced": dc_balanced,
            "ideal_bogoliubov_matched_to_G_s": {
                "S_linear": float(ideal),
                "S_dB": float(10 * np.log10(ideal)),
            },
            "symmetric_external_efficiency": float(eta_ext),
            "completed_unweighted_after_external_loss": {
                "S_linear": float(detected),
                "S_dB": float(10 * np.log10(detected)),
            },
            "ideal_bogoliubov_after_external_loss": {
                "S_linear": float(ideal_detected),
                "S_dB": float(10 * np.log10(ideal_detected)),
            },
        },
        "repository_literature_benchmark_comparison": {
            "benchmark_source": "README.md, Sim et al. 85Rb optimum repository benchmark",
            "benchmark_G_s": literature_gain,
            "model_G_s": G_s,
            "G_s_discrepancy": G_s - literature_gain,
            "G_s_relative_discrepancy_pct": 100 * (G_s - literature_gain) / literature_gain,
            "benchmark_squeezing_dB": literature_squeezing_db,
            "mathematical_dilation_detected_squeezing_dB": float(
                10 * np.log10(detected)
            ),
            "squeezing_discrepancy_dB": float(
                10 * np.log10(detected) - literature_squeezing_db
            ),
            "benchmark_bandwidth_MHz": literature_bandwidth_mhz,
            "model_bandwidth_MHz": None,
            "comparison_status": (
                "diagnostic only: mathematical dilation is not a microscopic "
                "or frequency-dependent experimental prediction"
            ),
        },
    }


@contextmanager
def _patched_reference_floquet(n_f: int):
    original = reference.floquet

    def solve(L0, Cp, Cm, omega_beat, deff_axis):
        return floquet_solve_truncated(L0, Cp, Cm, omega_beat, deff_axis, n_f)

    reference.floquet = solve
    try:
        yield
    finally:
        reference.floquet = original


def _phase_deg(value: complex) -> float:
    return float(np.degrees(np.angle(value)))


def _wrap_deg(value: float) -> float:
    return float((value + 180.0) % 360.0 - 180.0)


def run_model_case(
    *,
    n_f: int,
    velocity_step: float,
    velocity_cutoff: float,
    coarse: int,
    window: float,
) -> dict[str, float]:
    """Run one exact reduced-model case and extract fixed/legacy-optimum metrics."""
    handle, tmp_name = tempfile.mkstemp(suffix=".npz")
    os.close(handle)
    try:
        started = time.perf_counter()
        with _patched_reference_floquet(n_f), blas_single_thread():
            result = reference.run(
                **MODEL_ARGS,
                coarse=coarse,
                window=window,
                vstep=velocity_step,
                vcut=velocity_cutoff,
                save=tmp_name,
            )
        elapsed = time.perf_counter() - started
        with np.load(tmp_name) as data:
            delta_axis_mhz = data["delta_ax"] / (2 * np.pi * 1e6)
            i_fixed = int(np.argmin(np.abs(delta_axis_mhz - FIXED_DELTA_MHZ)))
            i_opt = int(result["i"])
            Gs = data["Gs"]
            Gc = data["Gc"]
            avg_sc = data["avg_sc"]
            avg_cs = data["avg_cs"]
            fixed = {
                "delta_mhz": float(delta_axis_mhz[i_fixed]),
                "G_s": float(Gs[i_fixed]),
                "G_c": float(Gc[i_fixed]),
                "gain_gap": float(Gs[i_fixed] - Gc[i_fixed]),
                "arg_chi_sc_deg": _phase_deg(avg_sc[i_fixed]),
                "arg_chi_cs_deg": _phase_deg(avg_cs[i_fixed]),
            }
            optimum = {
                "delta_star_mhz": float(delta_axis_mhz[i_opt]),
                "G_s": float(Gs[i_opt]),
                "G_c": float(Gc[i_opt]),
                "gain_gap": float(Gs[i_opt] - Gc[i_opt]),
                "arg_chi_sc_deg": _phase_deg(avg_sc[i_opt]),
                "arg_chi_cs_deg": _phase_deg(avg_cs[i_opt]),
                "legacy_xi_finite_db": float(data["Sf"][i_opt]),
                "legacy_xi_ideal_db": float(data["Si"][i_opt]),
            }
            scan_step = float(np.median(np.diff(delta_axis_mhz)))
            n_velocity = int(data["v"].size)
            n_deff = int(data["deff"].size)
        return {
            "N_F": n_f,
            "velocity_step_m_per_s": velocity_step,
            "velocity_cutoff_sigma": velocity_cutoff,
            "n_velocity": n_velocity,
            "n_delta_eff": n_deff,
            "scan_points": coarse,
            "scan_half_window_ghz": window,
            "scan_step_mhz": scan_step,
            "runtime_s": elapsed,
            "segment_od": float(result["seg_od"]),
            "fixed_operating_point": fixed,
            "legacy_reduced_objective_optimum": optimum,
        }
    finally:
        try:
            os.remove(tmp_name)
        except FileNotFoundError:
            pass


def _relative_percent(new: float, old: float) -> float:
    return 100.0 * abs(new - old) / max(abs(new), 1e-300)


def floquet_successive_changes(rows: list[dict]) -> list[dict[str, float]]:
    changes = []
    for old, new in zip(rows, rows[1:]):
        f0 = old["fixed_operating_point"]
        f1 = new["fixed_operating_point"]
        o0 = old["legacy_reduced_objective_optimum"]
        o1 = new["legacy_reduced_objective_optimum"]
        changes.append(
            {
                "from_N_F": old["N_F"],
                "to_N_F": new["N_F"],
                "fixed_G_s_relative_change_pct": _relative_percent(f1["G_s"], f0["G_s"]),
                "fixed_G_c_relative_change_pct": _relative_percent(f1["G_c"], f0["G_c"]),
                "fixed_gain_gap_relative_change_pct": _relative_percent(
                    f1["gain_gap"], f0["gain_gap"]
                ),
                "fixed_arg_chi_sc_change_deg": abs(
                    _wrap_deg(f1["arg_chi_sc_deg"] - f0["arg_chi_sc_deg"])
                ),
                "fixed_arg_chi_cs_change_deg": abs(
                    _wrap_deg(f1["arg_chi_cs_deg"] - f0["arg_chi_cs_deg"])
                ),
                "delta_star_change_mhz": abs(o1["delta_star_mhz"] - o0["delta_star_mhz"]),
            }
        )
    return changes


def _errors_to_reference(rows: list[dict]) -> list[dict[str, float]]:
    ref = rows[-1]["fixed_operating_point"]
    out = []
    for row in rows:
        fixed = row["fixed_operating_point"]
        out.append(
            {
                "velocity_step_m_per_s": row["velocity_step_m_per_s"],
                "velocity_cutoff_sigma": row["velocity_cutoff_sigma"],
                "n_velocity": row["n_velocity"],
                "G_s": fixed["G_s"],
                "G_c": fixed["G_c"],
                "gain_gap": fixed["gain_gap"],
                "arg_chi_sc_deg": fixed["arg_chi_sc_deg"],
                "arg_chi_cs_deg": fixed["arg_chi_cs_deg"],
                "G_s_error_to_last_pct": _relative_percent(ref["G_s"], fixed["G_s"]),
                "G_c_error_to_last_pct": _relative_percent(ref["G_c"], fixed["G_c"]),
                "gain_gap_error_to_last_pct": _relative_percent(
                    ref["gain_gap"], fixed["gain_gap"]
                ),
                "arg_chi_sc_error_to_last_deg": abs(
                    _wrap_deg(fixed["arg_chi_sc_deg"] - ref["arg_chi_sc_deg"])
                ),
                "arg_chi_cs_error_to_last_deg": abs(
                    _wrap_deg(fixed["arg_chi_cs_deg"] - ref["arg_chi_cs_deg"])
                ),
            }
        )
    return out


def _fmt(value: float, digits: int = 6) -> str:
    return f"{value:.{digits}g}"


def _fmt_complex_payload(value: dict[str, float], digits: int = 9) -> str:
    real = value["real"]
    imag = value["imag"]
    if abs(imag) < 5e-15:
        return f"{real:.{digits}g}"
    return f"{real:.{digits}g}{imag:+.{digits}g}i"


def _append_complex_matrix_table(
    lines: list[str], label: str, matrix: list[list[dict[str, float]]], unit: str
) -> None:
    lines += [
        "",
        f"**{label}** {unit}",
        "",
        "| row | column 1 | column 2 |",
        "|---:|---:|---:|",
    ]
    for index, row in enumerate(matrix, start=1):
        lines.append(
            f"| {index} | `{_fmt_complex_payload(row[0])}` | "
            f"`{_fmt_complex_payload(row[1])}` |"
        )


def make_markdown(audit: dict) -> str:
    lines = [
        "# Reduced-model Floquet and velocity convergence audit",
        "",
        "This is an analysis-only audit of `analysis/analytic_reconstruction/ref_solver.py`; "
        "it does not modify the LaTeX report or production solver.",
        "",
        "**Classification warning:** the initial Floquet/velocity gain tables inherit the archived "
        "`ref_solver.py` propagation, which uses dressed optical wave numbers and "
        "a refractive phase mismatch together. They isolate Floquet/velocity "
        "numerics inside that shared implementation; they are **not** predictions "
        "of the corrected no-double-count Option-A propagation. A separate corrected "
        "Option-A literature-point diagnostic appears later.",
        "",
        "## Solver self-checks",
        "",
        "| check | max absolute difference |",
        "|---|---:|",
    ]
    for key, value in audit["solver_self_checks"].items():
        lines.append(f"| `{key}` | {value:.3e} |")

    lines += [
        "",
        "## Floquet truncation at the common operating point",
        "",
        f"Common point: $\\Delta/2\\pi=-1.50$ GHz, $T=110$ C, "
        f"$\\delta/2\\pi={FIXED_DELTA_MHZ:.0f}$ MHz. The velocity grid is "
        "5 m/s to 3 sigma. Phases are modulo 360 degrees.",
        "",
        "| N_F | G_s | G_c | G_s-G_c | arg chi_sc (deg) | arg chi_cs (deg) |",
        "|---:|---:|---:|---:|---:|---:|",
    ]
    for row in audit["floquet_rows"]:
        f = row["fixed_operating_point"]
        lines.append(
            f"| {row['N_F']} | {_fmt(f['G_s'])} | {_fmt(f['G_c'])} | "
            f"{_fmt(f['gain_gap'])} | {_fmt(f['arg_chi_sc_deg'])} | "
            f"{_fmt(f['arg_chi_cs_deg'])} |"
        )

    lines += [
        "",
        "| change | rel. G_s | rel. G_c | rel. gap | phase sc | phase cs | delta-star shift |",
        "|---|---:|---:|---:|---:|---:|---:|",
    ]
    for row in audit["floquet_successive_changes"]:
        lines.append(
            f"| {row['from_N_F']} to {row['to_N_F']} | "
            f"{_fmt(row['fixed_G_s_relative_change_pct'])}% | "
            f"{_fmt(row['fixed_G_c_relative_change_pct'])}% | "
            f"{_fmt(row['fixed_gain_gap_relative_change_pct'])}% | "
            f"{_fmt(row['fixed_arg_chi_sc_change_deg'])} deg | "
            f"{_fmt(row['fixed_arg_chi_cs_change_deg'])} deg | "
            f"{_fmt(row['delta_star_change_mhz'])} MHz |"
        )

    lines += [
        "",
        "## Legacy reduced-objective minimizer",
        "",
        "`delta_star` below is only the minimizer of the reconstruction's legacy, "
        "gain-only squeezing objective (including its gap gate); it is not a "
        "commutator-preserving quantum prediction. Scan spacing is 5 MHz.",
        "",
        "| N_F | delta_star (MHz) | G_s | G_c | gap | arg chi_sc | arg chi_cs | legacy xi finite |",
        "|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for row in audit["floquet_rows"]:
        o = row["legacy_reduced_objective_optimum"]
        lines.append(
            f"| {row['N_F']} | {_fmt(o['delta_star_mhz'])} | {_fmt(o['G_s'])} | "
            f"{_fmt(o['G_c'])} | {_fmt(o['gain_gap'])} | "
            f"{_fmt(o['arg_chi_sc_deg'])} | {_fmt(o['arg_chi_cs_deg'])} | "
            f"{_fmt(o['legacy_xi_finite_db'])} dB |"
        )

    lines += [
        "",
        "## One-dimensional velocity-step refinement",
        "",
        f"All rows use N_F=3, cutoff 5 sigma, and the fixed "
        f"$\\delta/2\\pi={FIXED_DELTA_MHZ:.0f}$ MHz point. Errors are relative "
        "to the last row.",
        "",
        "| dv (m/s) | points | G_s | G_c | gap | G_s err. | gap err. | phase-sc err. |",
        "|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for row in audit["velocity_step_errors"]:
        lines.append(
            f"| {_fmt(row['velocity_step_m_per_s'])} | {row['n_velocity']} | "
            f"{_fmt(row['G_s'])} | {_fmt(row['G_c'])} | {_fmt(row['gain_gap'])} | "
            f"{_fmt(row['G_s_error_to_last_pct'])}% | "
            f"{_fmt(row['gain_gap_error_to_last_pct'])}% | "
            f"{_fmt(row['arg_chi_sc_error_to_last_deg'])} deg |"
        )

    lines += [
        "",
        "## One-dimensional velocity-cutoff refinement",
        "",
        f"All rows use N_F=3, dv=2.5 m/s, and the fixed "
        f"$\\delta/2\\pi={FIXED_DELTA_MHZ:.0f}$ MHz point. Errors are relative "
        "to the last row.",
        "",
        "| cutoff (sigma) | points | G_s | G_c | gap | G_s err. | gap err. | phase-sc err. |",
        "|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for row in audit["velocity_cutoff_errors"]:
        lines.append(
            f"| {_fmt(row['velocity_cutoff_sigma'])} | {row['n_velocity']} | "
            f"{_fmt(row['G_s'])} | {_fmt(row['G_c'])} | {_fmt(row['gain_gap'])} | "
            f"{_fmt(row['G_s_error_to_last_pct'])}% | "
            f"{_fmt(row['gain_gap_error_to_last_pct'])}% | "
            f"{_fmt(row['arg_chi_sc_error_to_last_deg'])} deg |"
        )

    option = audit["option_a_literature_point"]
    op = option["parameters"]
    final = option["option_a_N_F_3"]
    dilation = option["minimum_vacuum_mathematical_dilation"]
    benchmark = option["repository_literature_benchmark_comparison"]
    lines += [
        "",
        "## Corrected Option-A literature-point diagnostic",
        "",
        "This section is separate from the archived tables above. It uses bare, "
        "frequency-specific optical wave numbers in the susceptibility terms and "
        "only vacuum/geometric phase mismatch. No refractive-index contribution is "
        "inserted into the mismatch.",
        "",
        f"Operating point: $\\Delta/2\\pi={op['D_GHz']:+.3f}$ GHz, "
        f"$\\delta/2\\pi={op['delta_MHz']:+.3f}$ MHz, "
        f"$T={op['temperature_C']:.1f}$ C, pump={1e3*op['P_pump']:.0f} mW, "
        f"seed={1e6*op['P_probe']:.0f} uW, $L={1e3*op['L']:.1f}$ mm, "
        f"$\\theta={op['theta_deg']:.2f}$ deg. The one-dimensional velocity grid "
        f"has {op['n_velocity']} points ($dv={op['velocity_step_m_per_s']:.1f}$ m/s, "
        f"cutoff {op['velocity_cutoff_sigma']:.0f} sigma). Angular two-photon Doppler "
        "broadening is not included.",
        "",
        "All rows below are evaluated at the same fixed "
        "$\\delta/2\\pi=-8$ MHz literature point; no detuning optimization is "
        "mixed into this table.",
        "",
        "| N_F | G_s | G_c | G_s-G_c | arg chi_sc (deg) | arg chi_cs (deg) |",
        "|---:|---:|---:|---:|---:|---:|",
    ]
    for row in option["floquet_convergence"]:
        lines.append(
            f"| {row['N_F']} | {_fmt(row['G_s'], 9)} | {_fmt(row['G_c'], 9)} | "
            f"{_fmt(row['gain_gap'], 9)} | {_fmt(row['arg_chi_sc_deg'], 9)} | "
            f"{_fmt(row['arg_chi_cs_deg'], 9)} |"
        )
    lines += [
        "",
        "### Weak-field reference-amplitude check",
        "",
        "The atomic response is still the inherited approximate four-level model, "
        "and its finite seed/reference field enters the steady solve. The following "
        "N_F=3 check changes only that reference from 2 to 8 uW; it tests numerical "
        "weak-field linearity, not microscopic or experimental validity.",
        "",
        "| seed (uW) | Omega_s/2pi (MHz) | G_s | G_c | gap | dG_s vs 2uW | dgap vs 2uW |",
        "|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for row in option["seed_reference_linearity"]:
        lines.append(
            f"| {row['probe_power_uW']:.0f} | "
            f"{row['probe_rabi_over_2pi_MHz']:.6f} | "
            f"{row['G_s']:.9f} | {row['G_c']:.9f} | "
            f"{row['gain_gap']:.9f} | "
            f"{row['G_s_change_from_2uW_pct']:+.6f}% | "
            f"{row['gain_gap_change_from_2uW_pct']:+.6f}% |"
        )
    lines += [
        "",
        f"Vacuum/geometric mismatch: $\\Delta k_{{\\rm vac}}="
        f"{final['delta_k_vac_per_m']:.9f}\\,\\mathrm{{m^{{-1}}}}$.",
    ]
    _append_complex_matrix_table(lines, "M", final["M_per_m"], "[m^-1]")
    _append_complex_matrix_table(lines, "T = exp(M L)", final["T"], "[dimensionless]")

    lines += [
        "",
        "## Minimum-vacuum mathematical commutator completion",
        "",
        "**This is a mathematical dilation only.** It is not microscopic atomic "
        "diffusion, is not frequency dependent, and must not be presented as a "
        "Langevin-corrected squeezing-spectrum prediction. It additionally assumes "
        "the classical two-mode amplitudes have canonical photon-flux normalization.",
        "",
        f"For $J=\\mathrm{{diag}}(1,-1)$, the eigenvalues of "
        f"$K=-(MJ+JM^\\dagger)$ are "
        f"{dilation['K_eigenvalues_per_m'][0]:.9f} and "
        f"{dilation['K_eigenvalues_per_m'][1]:.9f} m^-1. Both are positive, so "
        "$J_f=I_2$ for the displayed eigenfactor.",
    ]
    _append_complex_matrix_table(lines, "K", dilation["K_per_m"], "[m^-1]")
    _append_complex_matrix_table(lines, "B", dilation["B_per_sqrt_m"], "[m^-1/2]")
    _append_complex_matrix_table(lines, "D_min = B B^dagger / 2", dilation["D_min_per_m"], "[m^-1]")
    _append_complex_matrix_table(lines, "V_out for vacuum input", dilation["V_out"], "[dimensionless]")

    lines += [
        "",
        "| commutator/diffusion check | max residual |",
        "|---|---:|",
        f"| Bare transfer $TJT^\\dagger-J$ | "
        f"{dilation['bare_commutator_residual_max']:.3e} |",
        f"| Factorization $BJ_fB^\\dagger-K$ | "
        f"{dilation['factorization_residual_max']:.3e} m^-1 |",
        f"| Completed output commutator | "
        f"{dilation['completed_commutator_residual_max']:.3e} |",
        f"| 200-to-400 point covariance-integral change | "
        f"{dilation['quadrature_200_to_400_covariance_change_max']:.3e} |",
        "",
        "Here `max residual` means the entrywise norm "
        "$\\max_{ij}|R_{ij}|$.",
        "",
        "| bright-seed diagnostic | S_- | dB | classification |",
        "|---|---:|---:|---|",
        f"| Bare T only | {dilation['transfer_only_invalid_unweighted']['S_linear']:.7f} | "
        f"{dilation['transfer_only_invalid_unweighted']['S_dB']:.3f} | invalid; commutator not restored |",
        f"| Minimum-vacuum dilation, unweighted | "
        f"{dilation['completed_unweighted']['S_linear']:.7f} | "
        f"{dilation['completed_unweighted']['S_dB']:.3f} | mathematical diagnostic |",
        f"| Minimum-vacuum dilation, DC-balanced | "
        f"{dilation['completed_dc_balanced']['S_linear']:.7f} | "
        f"{dilation['completed_dc_balanced']['S_dB']:.3f} | mathematical diagnostic |",
        f"| Ideal Bogoliubov matched to G_s | "
        f"{dilation['ideal_bogoliubov_matched_to_G_s']['S_linear']:.7f} | "
        f"{dilation['ideal_bogoliubov_matched_to_G_s']['S_dB']:.3f} | counterfactual benchmark |",
        f"| Minimum dilation after external eta={dilation['symmetric_external_efficiency']:.4f} | "
        f"{dilation['completed_unweighted_after_external_loss']['S_linear']:.7f} | "
        f"{dilation['completed_unweighted_after_external_loss']['S_dB']:.3f} | external-loss diagnostic |",
        "",
        f"Against the repository literature benchmark (`README.md`, Sim et al. "
        f"85Rb optimum), the corrected Option-A "
        f"mean-field result is $G_s={benchmark['model_G_s']:.3f}$ versus "
        f"approximately {benchmark['benchmark_G_s']:.1f} "
        f"({benchmark['G_s_relative_discrepancy_pct']:.1f}%). The mathematical "
        f"dilation after external loss is "
        f"{benchmark['mathematical_dilation_detected_squeezing_dB']:.3f} dB versus "
        f"the reported scale near {benchmark['benchmark_squeezing_dB']:.1f} dB. "
        "No bandwidth comparison is available because the static model has no "
        "spectrum-analyzer frequency.",
    ]

    angular = audit["geometry_diagnostic"]
    lines += [
        "",
        "## Interpretation and API limitations",
        "",
        "- The production APIs `gabes.core.floquet_solve(...)` and "
        "`gabes.kernels.floquet_chi_grid(...)` are fixed to N_F=1. "
        "`gabes.schemes.fwm.chi_matrix_table(...)` selects that fused kernel when "
        "Numba is available and exposes no truncation-order argument.",
        "- `gabes.schemes.fwm.compute_spectrum(...)` exposes `velocity_step` and "
        "`velocity_cutoff`, but not N_F. Its one-dimensional path calls "
        "`gabes.doppler.velocity_grid(...)`, `build_Delta_eff_axis(...)`, and "
        "`doppler_average(...)`.",
        "- The N_F=1 result is not Floquet-converged: N_F=2 changes the common-point "
        "gains substantially. N_F=2 and N_F=3 agree to the precision reported here.",
        "- Velocity refinement here converges only the existing collinear integral "
        "$\\Delta_{eff}=\\Delta-kv$. The current susceptibility API keeps delta "
        "independent of velocity and therefore cannot represent the crossing-angle "
        "two-photon Doppler distribution.",
        f"- At theta={angular['theta_deg']:.2f} deg, sigma_v="
        f"{angular['sigma_v_m_per_s']:.3f} m/s, and lambda="
        f"{angular['lambda_nm']:.3f} nm, the omitted one-sigma angular width is "
        f"{angular['angular_two_photon_sigma_mhz']:.3f} MHz.",
        "- A correct geometry extension must separate the lab beat frequency from "
        "the velocity-shifted atomic two-photon detuning. Reusing the current "
        "`floquet_chi_grid` with `delta_eff` would incorrectly shift both because "
        "that kernel computes `omega_beat = omega_hf + branch*delta` internally.",
        "- The continued-fraction extension assumes the same periodic Hamiltonian "
        "with only +/-1 Fourier couplings. It is exact for that finite truncation, "
        "but it does not repair the pump steady-state, quantum-Langevin, or "
        "four-level-model limitations.",
        "- The initial/common-point convergence tables retain `ref_solver.py`'s dressed-k plus "
        "refractive-mismatch convention. Use them only to diagnose convergence of "
        "the archived calculation. The separate Option-A section above supplies "
        "the corrected bare-k/vacuum-mismatch literature-point calculation.",
    ]
    return "\n".join(lines) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=HERE / "generated",
        help="Generated JSON/Markdown directory (default: analysis-local generated/)",
    )
    args = parser.parse_args()

    checks = solver_self_checks()

    # Five-MHz scan spacing over the report's full +/-700 MHz window.
    floquet_rows = [
        run_model_case(
            n_f=n_f,
            velocity_step=5.0,
            velocity_cutoff=3.0,
            coarse=281,
            window=0.7,
        )
        for n_f in (1, 2, 3)
    ]

    # Keep the report's 81-point scan for the quadrature refinements. This keeps
    # ref_solver's scan-global segment-OD prescription identical between rows and
    # includes delta=-280 MHz exactly.
    cache: dict[tuple[float, float], dict] = {}

    def velocity_case(dv: float, cutoff: float) -> dict:
        key = (dv, cutoff)
        if key not in cache:
            cache[key] = run_model_case(
                n_f=3,
                velocity_step=dv,
                velocity_cutoff=cutoff,
                coarse=81,
                window=0.7,
            )
        return cache[key]

    # Five sigma is required here: the nominal pump detuning is so large that the
    # 3-sigma truncation has not yet reached a stable tail contribution.  At the
    # converged cutoff the 5 m/s grid is already extremely stable, but retain the
    # finer rows to demonstrate it rather than assume it.
    velocity_step_rows = [velocity_case(dv, 5.0) for dv in (10.0, 5.0, 2.5, 1.25)]
    velocity_cutoff_rows = [
        velocity_case(2.5, cut) for cut in (2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0)
    ]

    # Separate corrected propagation at the experimentally established literature
    # point. This does not reuse ref_solver.run's dressed-k mismatch convention.
    option_a = option_a_literature_diagnostic()

    sigma_v = math.sqrt(reference.KB * MODEL_ARGS["T"] / reference.MASS)
    angular_sigma_mhz = (
        math.radians(MODEL_ARGS["theta_deg"]) * sigma_v / reference.LAM / 1e6
    )
    audit = {
        "scope": "Archived reduced-model numerical isolation test; no production/report edits",
        "propagation_classification": {
            "implementation": "analysis/analytic_reconstruction/ref_solver.py",
            "inherits_dressed_k_plus_refractive_mismatch": True,
            "corrected_option_a": False,
            "permitted_use": "Floquet and velocity numerical-isolation diagnostics only",
        },
        "parameters": {
            **MODEL_ARGS,
            "branch": -1,
            "fixed_delta_mhz": FIXED_DELTA_MHZ,
        },
        "solver_self_checks": checks,
        "floquet_rows": floquet_rows,
        "floquet_successive_changes": floquet_successive_changes(floquet_rows),
        "velocity_step_rows": velocity_step_rows,
        "velocity_step_errors": _errors_to_reference(velocity_step_rows),
        "velocity_cutoff_rows": velocity_cutoff_rows,
        "velocity_cutoff_errors": _errors_to_reference(velocity_cutoff_rows),
        "geometry_diagnostic": {
            "theta_deg": MODEL_ARGS["theta_deg"],
            "sigma_v_m_per_s": sigma_v,
            "lambda_nm": reference.LAM * 1e9,
            "angular_two_photon_sigma_mhz": angular_sigma_mhz,
            "implemented_in_audited_solver": False,
        },
        "option_a_literature_point": option_a,
    }

    args.output_dir.mkdir(parents=True, exist_ok=True)
    json_path = args.output_dir / "convergence_audit.json"
    md_path = args.output_dir / "convergence_audit.md"
    json_path.write_text(json.dumps(audit, indent=2), encoding="utf-8")
    md_path.write_text(make_markdown(audit), encoding="utf-8")
    print(md_path)
    print(json_path)


if __name__ == "__main__":
    main()
