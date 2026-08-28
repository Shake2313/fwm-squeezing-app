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
numbers plus vacuum/geometric mismatch).  It also constructs one algebraic
commutator/diffusion fixture with a chosen vacuum-reservoir covariance.  That
fixture is not a Caves bound, microscopic atomic diffusion, or a squeezing
spectrum.

Run from the repository root::

    python analysis/squeezing/analytic_reconstruction/convergence_audit.py

Outputs are written below ``analysis/squeezing/analytic_reconstruction/generated/``.
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
from scipy.linalg import expm as scipy_expm


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import ref_solver as reference  # noqa: E402
from gabes import (  # noqa: E402
    constants as gabes_constants,
    doppler as gabes_doppler,
    observables as gabes_observables,
)
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


def _literature_context(
    probe_power: float | None = None,
    gamma_gg_floor_khz: float | None = None,
) -> dict:
    """Build shared atomic/velocity data for the corrected literature point."""
    p = LITERATURE_ARGS
    probe_power = p["P_probe"] if probe_power is None else float(probe_power)
    branch = int(p["branch"])
    delta = 2 * np.pi * p["delta_MHz"] * 1e6
    Delta = 2 * np.pi * p["D_GHz"] * 1e9
    Op = gabes_constants.rabi_freq(p["P_pump"], p["w_pump"])
    Os = gabes_constants.rabi_freq(probe_power, p["w_probe"])
    density = fwm_scheme.hyperfine.number_density(p["T"])
    current_floor_khz = gabes_constants.GAMMA_GG_2PI / 1e3
    floor_khz = (
        current_floor_khz
        if gamma_gg_floor_khz is None
        else float(gamma_gg_floor_khz)
    )
    floor_rad_s = 2 * np.pi * floor_khz * 1e3
    # The current production path represents the inherited gamma_gg floor as a
    # trace-preserving thermal transit reset, while Rb-Rb ground-coherence
    # collisions remain a separate pure-dephasing term.  Always call that path;
    # the sensitivity override changes only its transit-rate argument.
    atom = fwm_scheme.collisional_atom(
        p["T"],
        density,
        transit_rate=(None if gamma_gg_floor_khz is None else floor_rad_s),
    )
    transit_rate = float(atom.transit_reset_rate)
    gamma_collision = float(atom.ground_collision_dephasing_rate)
    gamma_opt = 0.5 * (
        fwm_scheme.hyperfine.self_broadened_gamma(density)
        - gabes_constants.GAMMA
    )
    atomic_model_path = (
        "gabes.schemes.fwm.collisional_atom -> "
        "gabes.schemes.fwm.thermal_transit_reset_superoperator + "
        "gabes.atoms.double_lambda_rb85"
    )
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
        "atomic_model_path": atomic_model_path,
        "gamma_gg_floor_over_2pi_kHz": float(transit_rate / (2 * np.pi * 1e3)),
        "gamma_gg_collision_over_2pi_kHz": float(
            gamma_collision / (2 * np.pi * 1e3)),
        "gamma_gg_total_over_2pi_kHz": float(
            (transit_rate + gamma_collision) / (2 * np.pi * 1e3)),
        "gamma_opt_over_2pi_kHz": float(max(gamma_opt, 0.0) / (2 * np.pi * 1e3)),
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


def _doppler_average_table(table: np.ndarray, context: dict) -> np.ndarray:
    """Average a ``(scan, Delta_eff)`` response on the audit velocity grid."""
    return gabes_doppler.doppler_average(
        np.asarray(table, dtype=complex),
        context["deff"],
        context["Delta"],
        context["velocity"],
        context["weights"],
    )


def _noncollinear_doppler_audit(context: dict) -> dict:
    """Grid/cutoff convergence of the separate two-dimensional reference."""
    p = LITERATURE_ARGS
    branch = int(p["branch"])
    effective_line_strength = (
        p["line_strength_residual"]
        * fwm_scheme.physical_coupling_norm(branch)
    )
    theta = np.deg2rad(p["theta_deg"])

    def evaluate(delta_mhz, order, cutoff_sigma, angle_rad=theta):
        delta_mhz = np.asarray(delta_mhz, dtype=float)
        delta = 2.0 * np.pi * delta_mhz * 1e6
        probe_GHz = (
            fwm_scheme.branch_center_GHz(p["D_GHz"], branch)
            + delta_mhz * 1e-3
        )
        k_pump, k_probe, k_conj = fwm_scheme.seeded_option_a_wavenumbers(
            p["D_GHz"], probe_GHz)
        started = time.perf_counter()
        response = fwm_scheme.pump_only_weak_response_noncollinear_reference(
            context["Op"],
            context["Op"],
            delta,
            context["Delta"],
            T=p["T"],
            pump_k_rad_m=k_pump,
            probe_k_axis_rad_m=k_probe,
            crossing_angle_rad=angle_rad,
            branch=branch,
            atom=context["atom"],
            analysis_frequency_axis_rad_s=(0.0,),
            quadrature_order=order,
            cutoff_sigma=cutoff_sigma,
        )
        mismatch = fwm_scheme.seeded_phase_mismatch_z(
            p["D_GHz"], probe_GHz, angle_deg=np.rad2deg(angle_rad))
        G_s, G_c, transfer = gabes_observables.gain_from_chi(
            response.chi_ss[0],
            response.chi_sc[0],
            response.chi_cs[0],
            response.chi_cc[0],
            k_probe,
            k_conj,
            p["L"],
            context["density"],
            dipole=gabes_constants.DIPOLE_D1,
            line_strength=effective_line_strength,
            delta_k_z=mismatch,
        )
        peak = int(np.nanargmax(G_s))
        return {
            "order_per_axis": int(order),
            "velocity_pairs": int(order * order),
            "cutoff_sigma": float(cutoff_sigma),
            "delta_scan_MHz": delta_mhz,
            "G_s": np.asarray(G_s, dtype=float),
            "G_c": np.asarray(G_c, dtype=float),
            "T_field": np.asarray(transfer, dtype=complex),
            "feature_delta_MHz": float(delta_mhz[peak]),
            "feature_G_s": float(G_s[peak]),
            "feature_G_c": float(G_c[peak]),
            "response": response,
            "runtime_seconds": float(time.perf_counter() - started),
        }

    coarse_delta = np.arange(-100.0, 100.0 + 1e-12, 2.0)
    coarse = evaluate(coarse_delta, 24, 5.0)
    feature_center = coarse["feature_delta_MHz"]
    feature_delta = np.arange(
        feature_center - 2.0, feature_center + 2.0 + 1e-12, 0.05)
    refinement = [evaluate(feature_delta, order, 5.0)
                  for order in (24, 32, 40)]
    cutoff_rows = [evaluate(feature_delta, 40, cutoff)
                   for cutoff in (4.5, 5.0)]
    final = refinement[-1]
    for row in refinement:
        row["feature_G_s_error_to_order40_pct"] = float(
            100.0 * (row["feature_G_s"] - final["feature_G_s"])
            / final["feature_G_s"])
        row["feature_position_error_to_order40_MHz"] = float(
            row["feature_delta_MHz"] - final["feature_delta_MHz"])
    cutoff_final = cutoff_rows[-1]
    for row in cutoff_rows:
        row["feature_G_s_error_to_5sigma_pct"] = float(
            100.0 * (row["feature_G_s"] - cutoff_final["feature_G_s"])
            / cutoff_final["feature_G_s"])
        row["feature_position_error_to_5sigma_MHz"] = float(
            row["feature_delta_MHz"] - cutoff_final["feature_delta_MHz"])

    nominal = evaluate([p["delta_MHz"]], 40, 5.0)
    nominal_response = nominal["response"]
    velocity = nominal_response.vz_axis_m_s
    weights = nominal_response.vz_weights
    one_dimensional = fwm_scheme.pump_only_weak_response_reference(
        context["Op"],
        context["Op"],
        [context["delta"]],
        context["Delta"] - nominal_response.pump_k_rad_m * velocity,
        branch=branch,
        atom=context["atom"],
        analysis_frequency_axis_rad_s=(0.0,),
    )
    one_dimensional_chi = np.einsum(
        "v,voi->oi", weights, one_dimensional.chi_matrix[0, 0])
    probe_GHz = np.asarray([
        fwm_scheme.branch_center_GHz(p["D_GHz"], branch)
        + p["delta_MHz"] * 1e-3
    ])
    _kp, k_probe, k_conj = fwm_scheme.seeded_option_a_wavenumbers(
        p["D_GHz"], probe_GHz)
    mismatch = fwm_scheme.seeded_phase_mismatch_z(
        p["D_GHz"], probe_GHz, angle_deg=p["theta_deg"])
    one_d_named = (
        np.asarray([one_dimensional_chi[0, 0]]),
        np.asarray([one_dimensional_chi[0, 1]]),
        np.asarray([np.conj(one_dimensional_chi[1, 0])]),
        np.asarray([np.conj(one_dimensional_chi[1, 1])]),
    )
    one_d_Gs, one_d_Gc, _ = gabes_observables.gain_from_chi(
        *one_d_named,
        k_probe,
        k_conj,
        p["L"],
        context["density"],
        dipole=gabes_constants.DIPOLE_D1,
        line_strength=effective_line_strength,
        delta_k_z=mismatch,
    )
    actual_budget = gabes_doppler.noncollinear_raman_rms_budget(
        p["T"], nominal_response.pump_k_rad_m,
        nominal_response.probe_k_axis_rad_m[0], theta)
    collinear_budget = gabes_doppler.noncollinear_raman_rms_budget(
        p["T"], nominal_response.pump_k_rad_m,
        nominal_response.probe_k_axis_rad_m[0], 0.0)
    numeric_rms = float(
        nominal_response.diagnostics["quadrature_raman_rms_rad_s"][0])

    def public_row(row, error_keys):
        payload = {
            "order_per_axis": row["order_per_axis"],
            "velocity_pairs": row["velocity_pairs"],
            "cutoff_sigma": row["cutoff_sigma"],
            "feature_delta_MHz": row["feature_delta_MHz"],
            "feature_G_s": row["feature_G_s"],
            "feature_G_c": row["feature_G_c"],
            "runtime_seconds": row["runtime_seconds"],
        }
        payload.update({key: row[key] for key in error_keys})
        return payload

    return {
        "classification": (
            "slow two-dimensional pump-only/trace-zero Maxwell reference; "
            "not the one-dimensional finite-seed production path"),
        "provenance": nominal_response.provenance,
        "lab_frequency_contract": {
            "lab_optical_beat_velocity_shifted": False,
            "rf_analysis_frequency_independent": True,
            "atomic_delta_eff_velocity_shifted": True,
            "angular_frequency_internal_units": True,
        },
        "rms_budget": {
            "temperature_K": float(p["T"]),
            "angle_deg": float(p["theta_deg"]),
            "analytic_total_MHz": float(actual_budget["total_rms_hz"] * 1e-6),
            "analytic_transverse_MHz": float(
                actual_budget["transverse_rms_hz"] * 1e-6),
            "analytic_axial_at_angle_kHz": float(
                actual_budget["axial_rms_hz"] * 1e-3),
            "analytic_collinear_residual_kHz": float(
                collinear_budget["total_rms_hz"] * 1e-3),
            "quadrature_total_MHz": float(numeric_rms / (2.0 * np.pi * 1e6)),
            "quadrature_error_pct": float(
                100.0 * (numeric_rms - actual_budget["total_rms_rad_s"])
                / actual_budget["total_rms_rad_s"]),
        },
        "coarse_feature_scan": {
            "delta_min_MHz": float(coarse_delta[0]),
            "delta_max_MHz": float(coarse_delta[-1]),
            "step_MHz": float(coarse_delta[1] - coarse_delta[0]),
            "feature_delta_MHz": coarse["feature_delta_MHz"],
        },
        "refined_feature_scan": {
            "delta_min_MHz": float(feature_delta[0]),
            "delta_max_MHz": float(feature_delta[-1]),
            "step_MHz": float(feature_delta[1] - feature_delta[0]),
        },
        "grid_refinement": [public_row(row, (
            "feature_G_s_error_to_order40_pct",
            "feature_position_error_to_order40_MHz",
        )) for row in refinement],
        "cutoff_refinement": [public_row(row, (
            "feature_G_s_error_to_5sigma_pct",
            "feature_position_error_to_5sigma_MHz",
        )) for row in cutoff_rows],
        "acceptance": {
            "successive_grid_gain_change_below_1pct": bool(
                abs(refinement[-2]["feature_G_s_error_to_order40_pct"]) < 1.0),
            "successive_grid_feature_shift_at_most_0p1MHz": bool(
                abs(refinement[-2][
                    "feature_position_error_to_order40_MHz"]) <= 0.1),
            "cutoff_gain_change_below_0p5pct": bool(
                abs(cutoff_rows[0]["feature_G_s_error_to_5sigma_pct"]) < 0.5),
            "cutoff_feature_shift_at_most_0p1MHz": bool(
                abs(cutoff_rows[0][
                    "feature_position_error_to_5sigma_MHz"]) <= 0.1),
        },
        "nominal_literature_point": {
            "delta_MHz": float(p["delta_MHz"]),
            "one_dimensional_G_s": float(one_d_Gs[0]),
            "one_dimensional_G_c": float(one_d_Gc[0]),
            "two_dimensional_G_s": nominal["feature_G_s"],
            "two_dimensional_G_c": nominal["feature_G_c"],
            "G_s_change_pct": float(
                100.0 * (nominal["feature_G_s"] - one_d_Gs[0]) / one_d_Gs[0]),
            "max_response_normalized_residual": float(
                nominal_response.diagnostics["max_response_normalized_residual"]),
            "max_response_trace_error": float(
                nominal_response.diagnostics["max_response_trace_error"]),
        },
        "remaining_limits": {
            "production_default": False,
            "finite_seed_2d_production_path": False,
            "beam_angular_distribution": False,
            "segmentwise_pump_recomputation": False,
            "microscopic_langevin_diffusion": False,
        },
    }
def _pump_reference_audit(context: dict) -> dict:
    """Full-scan pump-state/weak-response and pole-residue validation."""
    delta_mhz = np.asarray([-30.0, -8.0, 0.0, 12.0, 40.0])
    delta = 2 * np.pi * delta_mhz * 1e6
    analysis_mhz = np.asarray([0.0, 0.1, 1.0, 4.0])
    analysis = 2 * np.pi * analysis_mhz * 1e6
    started = time.perf_counter()
    reference_response = fwm_scheme.pump_only_weak_response_reference(
        context["Op"], context["Op"], delta, context["deff"],
        branch=-1, atom=context["atom"],
        analysis_frequency_axis_rad_s=analysis,
    )
    reference_runtime = time.perf_counter() - started
    reference_components = {
        "ss": reference_response.chi_ss[0],
        "sc": reference_response.chi_sc[0],
        "cs": reference_response.chi_cs[0],
        "cc": reference_response.chi_cc[0],
    }
    reference_average = {
        name: _doppler_average_table(value, context)
        for name, value in reference_components.items()
    }
    reference_scale = max(
        float(np.max(np.abs(value), initial=0.0))
        for value in reference_average.values())

    finite_rows = []
    seed_started = time.perf_counter()
    for fraction in (1.0, 0.3, 0.1, 0.03):
        seed_rabi = context["Os"] * fraction
        tables = fwm_scheme.chi_matrix_table(
            context["Op"], context["Op"], seed_rabi, seed_rabi,
            delta, context["deff"], -1, atom=context["atom"], n_f=3)
        finite_average = {
            name: _doppler_average_table(table, context)
            for name, table in zip(("ss", "cs", "sc", "cc"), tables)
        }
        component_errors = {
            name: float(np.max(np.abs(
                finite_average[name] - reference_average[name]), initial=0.0))
            for name in reference_average
        }
        worst = max(component_errors.values()) / max(
            reference_scale, np.finfo(float).tiny)
        finite_rows.append({
            "rabi_fraction": float(fraction),
            "seed_power_uW": float(1e6 * context["P_probe"] * fraction**2),
            "seed_rabi_over_2pi_MHz": float(
                seed_rabi / (2 * np.pi * 1e6)),
            "worst_normalized_chi_error": float(worst),
            "component_max_abs_errors_seconds": component_errors,
        })
    finite_seed_runtime = time.perf_counter() - seed_started
    slope_rows = finite_rows[1:]
    convergence_slope = float(np.polyfit(
        np.log([row["rabi_fraction"] for row in slope_rows]),
        np.log([row["worst_normalized_chi_error"] for row in slope_rows]),
        1,
    )[0])
    pairwise_slopes = np.diff(np.log([
        row["worst_normalized_chi_error"] for row in slope_rows
    ])) / np.diff(np.log([
        row["rabi_fraction"] for row in slope_rows
    ]))

    smallest_fraction = finite_rows[-1]["rabi_fraction"]
    seed_rabi = context["Os"] * smallest_fraction
    order_limit_errors = {}
    for n_f in (1, 2, 3):
        tables = fwm_scheme.chi_matrix_table(
            context["Op"], context["Op"], seed_rabi, seed_rabi,
            delta, context["deff"], -1, atom=context["atom"], n_f=n_f)
        finite_average = {
            name: _doppler_average_table(table, context)
            for name, table in zip(("ss", "cs", "sc", "cc"), tables)
        }
        order_limit_errors[str(n_f)] = float(max(
            np.max(np.abs(finite_average[name] - reference_average[name]),
                   initial=0.0)
            for name in reference_average
        ) / max(reference_scale, np.finfo(float).tiny))

    sample_indices = np.unique(np.asarray(
        [0, context["deff"].size // 2, context["deff"].size - 1], dtype=int))
    sample_deff = context["deff"][sample_indices]
    static_sample = reference_response.pump_state[sample_indices]
    mapped = fwm_scheme.pump_frame_to_seeded_harmonics(
        static_sample, n_f=3, branch=-1)
    Cp, Cm = fwm_scheme.sideband_template(
        context["Op"], context["Op"], 0.0, branch=-1)
    frame_error = 0.0
    for detuning in delta:
        L0 = build_liouvillian(
            fwm_scheme.static_hamiltonian_at_Deff_zero(
                context["Op"], context["Op"], 0.0, detuning, branch=-1),
            context["atom"],
        )
        periodic = fwm_scheme.floquet_solve_truncated(
            L0, Cp, Cm, fwm_scheme.seeded_sideband_beat(detuning, -1),
            sample_deff, context["atom"].S_v, fwm_scheme.N_LEVELS, 3,
            return_harmonics=True,
        )
        frame_error = max(frame_error, float(np.max(np.abs(periodic - mapped))))

    center_index = int(np.argmin(np.abs(context["deff"] - context["Delta"])))
    delta_index = int(np.flatnonzero(delta_mhz == -8.0)[0])
    pole_rows = []
    first_poles = None
    for analysis_index, analysis_value in enumerate(analysis):
        poles = fwm_scheme.pump_only_pole_residue_reference(
            context["Op"], context["Op"], delta[delta_index],
            context["deff"][center_index], branch=-1, atom=context["atom"],
            analysis_frequency_rad_s=analysis_value,
        )
        if first_poles is None:
            first_poles = poles
        direct = reference_response.chi_matrix[
            analysis_index, delta_index, center_index]
        pole_rows.append({
            "analysis_frequency_MHz": float(analysis_mhz[analysis_index]),
            "max_abs_direct_pole_difference_seconds": float(
                np.max(np.abs(poles["response"] - direct), initial=0.0)),
            "eigenvector_condition": float(poles["eigenvector_condition"]),
            "stationary_residue_max": float(poles["stationary_residue_max"]),
        })
    nonstationary = np.asarray(first_poles["nonstationary_mask"], dtype=bool)
    visibility = np.max(first_poles["visibility"], axis=(1, 2))
    ranked = np.flatnonzero(nonstationary)[
        np.argsort(visibility[nonstationary])[::-1]]
    visible_poles = []
    for index in ranked[:8]:
        eigenvalue = first_poles["eigenvalues"][index]
        visible_poles.append({
            "eigenvalue_per_s": _complex_payload(eigenvalue),
            "analysis_center_MHz": float(
                first_poles["analysis_pole_centers_rad_s"][index]
                / (2 * np.pi * 1e6)),
            "half_width_MHz": float(
                first_poles["half_widths_rad_s"][index]
                / (2 * np.pi * 1e6)),
            "max_residue_over_half_width": float(visibility[index]),
        })

    dc_analysis = fwm_scheme.OMEGA_HF - context["delta"]
    dc_reference = fwm_scheme.pump_only_weak_response_reference(
        context["Op"], context["Op"], [context["delta"]],
        [context["deff"][center_index]], branch=-1, atom=context["atom"],
        analysis_frequency_axis_rad_s=[dc_analysis],
    )
    return {
        "classification": (
            "slow independent pump-only/trace-zero weak-response reference; "
            "not the production finite-seed gain path"),
        "provenance": reference_response.provenance,
        "supported_branch": -1,
        "plus_branch_status": (
            "unsupported: inherited plus-branch finite-seed frame fails physical "
            "static-pump gauge parity"),
        "delta_scan_MHz": [float(value) for value in delta_mhz],
        "analysis_frequency_scan_MHz": [float(value) for value in analysis_mhz],
        "pump_state_diagnostics": reference_response.diagnostics,
        "frame_equivalence": {
            "sampled_delta_points": int(delta.size),
            "sampled_velocity_classes": int(sample_deff.size),
            "N_F": 3,
            "max_abs_static_to_floquet_difference": float(frame_error),
            "outer_harmonics_exactly_zero": True,
        },
        "finite_seed_to_infinitesimal": {
            "N_F": 3,
            "rows": finite_rows,
            "rabi_error_power_law_slope": convergence_slope,
            "pairwise_rabi_error_slopes": [
                float(value) for value in pairwise_slopes],
            "smallest_fraction_order_errors": order_limit_errors,
            "reference_full_scan_max_abs_seconds": reference_scale,
        },
        "trace_zero_dc_projection": {
            "analysis_frequency_MHz": float(dc_analysis / (2 * np.pi * 1e6)),
            "relative_resolvent_frequency_rad_s": float(
                dc_reference.relative_frequency_rad_s[0, 0]),
            "max_response_normalized_residual": float(
                dc_reference.diagnostics["max_response_normalized_residual"]),
            "max_response_trace_error": float(
                dc_reference.diagnostics["max_response_trace_error"]),
        },
        "pole_residue_parity": pole_rows,
        "most_visible_nonstationary_poles": visible_poles,
        "runtime_seconds": {
            "pump_reference_4x5_scan": float(reference_runtime),
            "finite_seed_N_F3_four_amplitudes": float(finite_seed_runtime),
        },
        "remaining_limits": {
            "two_dimensional_noncollinear_doppler": False,
            "segmentwise_pump_recomputation": False,
            "microscopic_langevin_diffusion": False,
        },
    }


def _independent_option_a_reference(
        context: dict, reduced_chi: tuple[complex, ...],
        effective_line_strength: float, k_probe: float, k_conj: float,
        delta_k_vac: float) -> dict[str, np.ndarray]:
    """Literal SI/Scipy reference, independent of production propagation code."""
    chi_ss, chi_sc, chi_cs, chi_cc = reduced_chi
    prefactor = (
        -2.0 * context["density"] * gabes_constants.DIPOLE_D1**2
        * effective_line_strength
        / (gabes_constants.EPS_0 * gabes_constants.HBAR)
    )
    chi_ss_phys = prefactor * chi_ss
    chi_sc_phys = prefactor * chi_sc
    chi_cs_phys = prefactor * chi_cs
    chi_cc_phys = prefactor * chi_cc
    M_field = np.array([
        [0.5j * k_probe * chi_ss_phys - 0.5j * delta_k_vac,
         0.5j * k_probe * chi_sc_phys],
        [-0.5j * k_conj * np.conj(chi_cs_phys),
         -0.5j * k_conj * np.conj(chi_cc_phys) + 0.5j * delta_k_vac],
    ], dtype=complex)
    transfer_field = scipy_expm(M_field * LITERATURE_ARGS["L"])
    area_probe = 0.5 * np.pi * LITERATURE_ARGS["w_probe"]**2
    area_conj = area_probe
    omega_probe = gabes_constants.C_LIGHT * k_probe
    omega_conj = gabes_constants.C_LIGHT * k_conj
    q_probe = np.sqrt(
        2.0 * gabes_constants.HBAR * omega_probe
        / (gabes_constants.EPS_0 * gabes_constants.C_LIGHT * area_probe))
    q_conj = np.sqrt(
        2.0 * gabes_constants.HBAR * omega_conj
        / (gabes_constants.EPS_0 * gabes_constants.C_LIGHT * area_conj))
    Q = np.diag([q_probe, q_conj])
    Q_inv = np.diag([1.0 / q_probe, 1.0 / q_conj])
    return {
        "M_field": M_field,
        "T_field": transfer_field,
        "Q": Q,
        "M_canonical": Q_inv @ M_field @ Q,
        "T_canonical": Q_inv @ transfer_field @ Q,
    }


def _option_a_transfer_from_reduced_chi(
    context: dict,
    n_f: int,
    reduced_chi: tuple[complex, ...],
    *,
    line_strength_residual: float | None = None,
) -> dict:
    """Option-A propagation of one fixed reduced-susceptibility tuple."""
    p = LITERATURE_ARGS
    chi_ss, chi_sc, chi_cs, chi_cc = reduced_chi
    ell_s = (
        p["line_strength_residual"]
        if line_strength_residual is None
        else float(line_strength_residual)
    )
    effective_line_strength = (
        ell_s
        * fwm_scheme.physical_coupling_norm(context["branch"])
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

    independent = _independent_option_a_reference(
        context, (chi_ss, chi_sc, chi_cs, chi_cc), effective_line_strength,
        k_probe, k_conj, delta_k_vac)

    # Production path is checked against the literal SI/Scipy reference above;
    # sharing the atomic chi is intentional, sharing Maxwell/Q code is not.
    chi_args = tuple(
        np.asarray([value], dtype=complex)
        for value in (chi_ss, chi_sc, chi_cs, chi_cc)
    )
    M_field = gabes_observables._gain_matrix_from_chi(
        *chi_args,
        np.asarray([k_probe]), np.asarray([k_conj]), context["density"],
        gabes_constants.DIPOLE_D1, effective_line_strength,
        delta_k_z=np.asarray([delta_k_vac]),
    )[0]
    _, _, transfer_stack = gabes_observables.gain_from_chi(
        *chi_args,
        np.asarray([k_probe]), np.asarray([k_conj]), p["L"], context["density"],
        dipole=gabes_constants.DIPOLE_D1,
        line_strength=effective_line_strength,
        delta_k_z=np.asarray([delta_k_vac]),
    )
    transfer_field = transfer_stack[0]
    area_probe = float(gabes_observables.gaussian_mode_area(p["w_probe"]))
    area_conj = area_probe
    canonical = gabes_observables.canonical_transfer_diagnostics(
        transfer_field,
        gabes_constants.C_LIGHT * k_probe,
        gabes_constants.C_LIGHT * k_conj,
        area_probe,
        area_conj,
    )
    M_canonical, Q = gabes_observables.canonical_transfer_from_field(
        M_field,
        gabes_constants.C_LIGHT * k_probe,
        gabes_constants.C_LIGHT * k_conj,
        area_probe,
        area_conj,
    )
    transfer = canonical["transfer_canonical"]
    parity = {
        "M_field_max_abs_difference": float(
            np.max(np.abs(M_field - independent["M_field"]))),
        "T_field_max_abs_difference": float(
            np.max(np.abs(transfer_field - independent["T_field"]))),
        "Q_max_abs_difference": float(
            np.max(np.abs(Q - independent["Q"]))),
        "M_canonical_max_abs_difference": float(
            np.max(np.abs(M_canonical - independent["M_canonical"]))),
        "T_canonical_max_abs_difference": float(
            np.max(np.abs(transfer - independent["T_canonical"]))),
    }
    for key, difference in parity.items():
        if difference > 2e-11:
            raise AssertionError(
                f"production Option-A differs from independent reference: "
                f"{key}={difference:.3e}")
    G_s = float(canonical["probe_power_gain"])
    G_c = float(canonical["conjugate_power_gain"])
    G_c_flux = float(canonical["conjugate_photon_flux_gain"])
    return {
        "N_F": int(n_f),
        "line_strength_residual": float(ell_s),
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
        "mode_area_probe_m2": area_probe,
        "mode_area_conjugate_m2": area_conj,
        "Q": _complex_matrix_payload(Q),
        "M_field_per_m": _complex_matrix_payload(M_field),
        "M_canonical_per_m": _complex_matrix_payload(M_canonical),
        "T_field": _complex_matrix_payload(transfer_field),
        "T_canonical": _complex_matrix_payload(transfer),
        "independent_reference_parity": parity,
        # Backward-compatible aliases now explicitly mean canonical matrices.
        "M_per_m": _complex_matrix_payload(M_canonical),
        "T": _complex_matrix_payload(transfer),
        "G_s": G_s,
        "G_c": G_c,
        "gain_gap": G_s - G_c,
        "G_c_photon_flux": G_c_flux,
        "photon_flux_gap": G_s - G_c_flux,
        "probe_power_uW": float(1e6 * context["P_probe"]),
        "probe_rabi_over_2pi_MHz": float(context["Os"] / (2 * np.pi * 1e6)),
        "arg_chi_sc_rad": float(np.angle(chi_sc)),
        "arg_chi_cs_rad": float(np.angle(chi_cs)),
        "arg_chi_sc_deg": float(np.degrees(np.angle(chi_sc))),
        "arg_chi_cs_deg": float(np.degrees(np.angle(chi_cs))),
        "evaluation_type": "fixed literature operating point",
        "delta_mhz": float(p["delta_MHz"]),
        "_M_array": M_canonical,
        "_T_array": transfer,
    }


def _option_a_transfer(context: dict, n_f: int) -> dict:
    """Recompute the atomic response, then apply corrected Option-A propagation."""
    return _option_a_transfer_from_reduced_chi(
        context,
        n_f,
        _literature_reduced_chi(context, n_f),
    )


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


def _response_gain_row(response: dict) -> dict[str, float]:
    """Select JSON-stable mean-field observables from an Option-A response."""
    return {
        "G_s": float(response["G_s"]),
        "G_c_power": float(response["G_c"]),
        "G_c_photon_flux": float(response["G_c_photon_flux"]),
        "photon_flux_gap": float(response["photon_flux_gap"]),
    }


def _reduced_chi_from_response(response: dict) -> tuple[complex, ...]:
    """Recover the stored (ss, sc, cs, cc) tuple without another atomic solve."""
    payload = response["reduced_chi_seconds"]
    return tuple(
        complex(payload[key]["real"], payload[key]["imag"])
        for key in ("ss", "sc", "cs", "cc")
    )


def _parameter_sensitivity(
    context: dict,
    current_response: dict,
) -> dict:
    """Illustrative parameter sensitivities with explicit layer boundaries."""
    p = LITERATURE_ARGS
    reduced_chi = _reduced_chi_from_response(current_response)

    # ell_s changes propagation only.  The N_F=3 atomic response is deliberately
    # frozen so this is not misread as a refit or as an uncertainty interval.
    ell_rows = []
    for ell_s in (0.666, 0.74, 0.814):
        response = (
            current_response
            if np.isclose(ell_s, p["line_strength_residual"], rtol=0.0, atol=1e-15)
            else _option_a_transfer_from_reduced_chi(
                context,
                3,
                reduced_chi,
                line_strength_residual=ell_s,
            )
        )
        ell_rows.append({
            "ell_s": float(ell_s),
            "effective_line_strength": float(response["effective_line_strength"]),
            **_response_gain_row(response),
        })
    ell_span = ell_rows[-1]["ell_s"] - ell_rows[0]["ell_s"]
    ell_derivatives = {
        f"d{key}_dell_s": float(
            (ell_rows[-1][key] - ell_rows[0][key]) / ell_span)
        for key in ("G_s", "G_c_power", "G_c_photon_flux", "photon_flux_gap")
    }

    # Recompute the N_F=3 atomic response when the ground-coherence floor moves.
    # The central 100-kHz row reuses the already-computed current response.
    gamma_rows = []
    current_floor = float(context["gamma_gg_floor_over_2pi_kHz"])
    for floor_khz in (90.0, 100.0, 110.0):
        if np.isclose(floor_khz, current_floor, rtol=0.0, atol=1e-12):
            gamma_context = context
            response = current_response
            response_source = "reused current N_F=3 atomic solve"
        else:
            gamma_context = _literature_context(gamma_gg_floor_khz=floor_khz)
            response = _option_a_transfer(gamma_context, 3)
            response_source = "recomputed N_F=3 atomic solve"
        gamma_rows.append({
            "gamma_gg_floor_over_2pi_kHz": float(floor_khz),
            "gamma_gg_collision_over_2pi_kHz": float(
                gamma_context["gamma_gg_collision_over_2pi_kHz"]),
            "gamma_gg_total_over_2pi_kHz": float(
                gamma_context["gamma_gg_total_over_2pi_kHz"]),
            "gamma_opt_over_2pi_kHz": float(
                gamma_context["gamma_opt_over_2pi_kHz"]),
            "atomic_response_source": response_source,
            **_response_gain_row(response),
        })
    gamma_span = (
        gamma_rows[-1]["gamma_gg_floor_over_2pi_kHz"]
        - gamma_rows[0]["gamma_gg_floor_over_2pi_kHz"]
    )
    gamma_derivatives = {
        f"d{key}_d_gamma_gg_floor_kHz": float(
            (gamma_rows[-1][key] - gamma_rows[0][key]) / gamma_span)
        for key in ("G_s", "G_c_power", "G_c_photon_flux", "photon_flux_gap")
    }

    # kappa is only the coefficient of the existing normalized pump-scatter
    # diagnostic.  It is not an atomic diffusion coefficient or a squeezing fit.
    current_kappa = float(fwm_scheme.HARDENED_PUMP_SCATTER_KAPPA)
    pump_scatter_value, pump_od = fwm_scheme._pump_scatter_noise(
        p["D_GHz"], p["T"], p["L"], current_kappa)
    pump_scatter_slope = float(1.0 - np.exp(-pump_od))
    if not np.isclose(
        pump_scatter_value,
        current_kappa * pump_scatter_slope,
        rtol=2e-15,
        atol=1e-15,
    ):
        raise AssertionError("pump-scatter diagnostic is inconsistent with kappa law")

    return {
        "classification": {
            "illustrative_not_parameter_uncertainty": True,
            "physical_squeezing_prediction": False,
            "microscopic_atomic_diffusion": False,
        },
        "ell_s_propagation_only": {
            "classification": (
                "illustrative +/-10% sweep; fixed stored/current N_F=3 reduced chi; "
                "not a fit or uncertainty interval"
            ),
            "atomic_response_recomputed": False,
            "reduced_chi_source": "option_a_N_F_3.reduced_chi_seconds",
            "rows": ell_rows,
            "central_finite_difference": ell_derivatives,
        },
        "gamma_gg_floor_mean_field": {
            "classification": (
                "actual N_F=3 atomic-response recomputation through the current "
                "thermal-transit-reset path with all non-gamma inputs fixed; "
                "central row reuses the identical current solve"
            ),
            "atomic_response_recomputed_for_endpoints": True,
            "rows": gamma_rows,
            "central_finite_difference_per_kHz": gamma_derivatives,
        },
        "kappa_pump_scatter_diagnostic": {
            "equation": "N_ps = kappa * (1 - exp(-OD_pump))",
            "classification": (
                "SQL-normalized technical pump-scatter diagnostic only; not "
                "physical squeezing or microscopic diffusion"
            ),
            "kappa_current": current_kappa,
            "OD_pump": float(pump_od),
            "dN_ps_dkappa": pump_scatter_slope,
            "N_ps_current": float(pump_scatter_value),
        },
    }


def _parameter_provenance(context: dict) -> dict:
    """Machine-readable sources and solver paths for every added table."""
    p = LITERATURE_ARGS
    return {
        "schema_version": 1,
        "current_atomic_path": context["atomic_model_path"],
        "parameters": {
            "ell_s": {
                "current_value": float(p["line_strength_residual"]),
                "units": "dimensionless",
                "source": "LITERATURE_ARGS.line_strength_residual",
                "status": "inherited unrefitted residual; not first-principles",
            },
            "kappa": {
                "current_value": float(fwm_scheme.HARDENED_PUMP_SCATTER_KAPPA),
                "units": "dimensionless",
                "source": "gabes.schemes.fwm.HARDENED_PUMP_SCATTER_KAPPA",
                "status": "uncalibrated technical-noise diagnostic coefficient",
            },
            "gamma_gg_floor_over_2pi_kHz": {
                "current_value": float(gabes_constants.GAMMA_GG_2PI / 1e3),
                "units": "kHz",
                "source": "gabes.constants.GAMMA_GG_2PI",
                "status": (
                    "thermal transit-reset rate; density collision dephasing is "
                    "added separately in the current atomic path"
                ),
            },
        },
        "table_solver_provenance": {
            "floquet_convergence": {
                "atomic_response": (
                    "analysis.squeezing.analytic_reconstruction.convergence_audit."
                    "floquet_solve_truncated"
                ),
                "atomic_model": context["atomic_model_path"],
                "propagation": (
                    "gabes.observables Option-A path, checked against literal "
                    "SI scipy.linalg.expm reference"
                ),
            },
            "seed_reference_linearity": {
                "atomic_response": "recomputed N_F=3 for each finite seed reference",
                "atomic_model": context["atomic_model_path"],
                "propagation": "same independently checked Option-A path",
            },
            "pump_only_weak_response_reference": {
                "atomic_response": (
                    "gabes.schemes.fwm.pump_only_weak_response_reference"
                ),
                "state": "static physical pump frame, trace-one 16x16 null solve",
                "response": (
                    "two-column complex Nambu derivative; projected trace-zero "
                    "DC and ordinary finite-frequency resolvent"
                ),
                "production_default": False,
            },
            "noncollinear_doppler_reference": {
                "atomic_response": (
                    "gabes.schemes.fwm."
                    "pump_only_weak_response_noncollinear_reference"
                ),
                "velocity_average": (
                    "tensor Gauss-Legendre quadrature of the truncated Maxwell "
                    "distribution over independent (v_z, v_x)"
                ),
                "frequency_contract": (
                    "Omega_beat,lab and Omega_SA fixed; only Delta_eff and "
                    "delta_eff depend on velocity"
                ),
                "production_default": False,
            },
            "ell_s_propagation_only": {
                "atomic_response": "stored current N_F=3 reduced susceptibility",
                "atomic_model": context["atomic_model_path"],
                "propagation": "Option-A rerun only; atomic response held fixed",
            },
            "gamma_gg_floor_mean_field": {
                "atomic_response": "N_F=3 recomputed at 90/100/110 kHz floor",
                "atomic_model": (
                    "production collisional_atom thermal-reset path with "
                    "sensitivity-only transit-rate override"
                ),
                "propagation": "same independently checked Option-A path",
            },
            "kappa_pump_scatter_diagnostic": {
                "solver": (
                    "gabes.schemes.fwm._pump_scatter_noise -> "
                    "gabes.schemes.absorption._hyperfine_alpha"
                ),
                "atomic_covariance": "not computed",
            },
            "algebraic_dilation_fixture": {
                "solver": "static eigenfactorization plus Gauss-Legendre integral",
                "atomic_diffusion": "not computed",
                "physical_bound": False,
            },
        },
    }


def option_a_literature_diagnostic() -> dict:
    """Reproducible Option-A transfer plus an algebraic dilation fixture.

    The fixture enforces the canonical commutator for the supplied static M with
    one chosen vacuum-reservoir covariance. It is neither a Caves bound nor an
    atomic Langevin derivation and must not be interpreted as a microscopic or
    frequency-dependent squeezing calculation.
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
                        "G_c_photon_flux",
                        "photon_flux_gap",
                        "arg_chi_sc_deg",
                        "arg_chi_cs_deg",
                    )
                }
            )
        pump_reference = _pump_reference_audit(context)
        noncollinear_reference = _noncollinear_doppler_audit(context)
        sensitivity = _parameter_sensitivity(context, convergence[-1])
        parameter_provenance = _parameter_provenance(context)
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
        row["photon_flux_gap_change_from_2uW_pct"] = float(
            100
            * (row["photon_flux_gap"] - weak_reference["photon_flux_gap"])
            / weak_reference["photon_flux_gap"]
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
    D_vacuum_fixture = 0.5 * B @ B.conj().T

    quadrature_points = 400
    commutator_integral = _constant_matrix_integral(
        M, K, p["L"], quadrature_points
    )
    diffusion_integral = _constant_matrix_integral(
        M, D_vacuum_fixture, p["L"], quadrature_points
    )
    # Independent lower-order quadrature check.
    commutator_integral_200 = _constant_matrix_integral(M, K, p["L"], 200)
    diffusion_integral_200 = _constant_matrix_integral(
        M, D_vacuum_fixture, p["L"], 200)

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
            "population_normalization": (
                "trace-normalized rho_ss supplies p_F once; external structural "
                "factor is 1/[2(2I+1)]"
            ),
            "dilation": (
                "algebraic commutator/diffusion fixture with one chosen vacuum "
                "reservoir covariance; not a Caves bound"
            ),
            "microscopic_atomic_diffusion": False,
            "frequency_dependent": False,
            "physical_squeezing_prediction": False,
            "weak_field_reference_linearity_tested": True,
            "weak_field_reference_test_range_uW": [2.0, 4.0, 8.0],
            "pump_only_trace_zero_reference_implemented": True,
            "pump_only_reference_is_production_default": False,
            "noncollinear_doppler_reference_implemented": True,
            "noncollinear_reference_is_production_default": False,
            "canonical_mode_normalization": (
                "explicit Q^-1 T_field Q using Gaussian A=pi*w^2/2; the "
                "literature audit uses an explicitly declared equal collected "
                "probe/conjugate area"
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
            "separate_two_dimensional_reference_available": True,
        },
        "floquet_convergence": convergence,
        "seed_reference_linearity": seed_linearity,
        "pump_only_weak_response_reference": pump_reference,
        "noncollinear_doppler_reference": noncollinear_reference,
        "parameter_provenance": parameter_provenance,
        "sensitivity": sensitivity,
        "option_a_N_F_3": final,
        "algebraic_dilation_fixture": {
            "classification": (
                "commutator-preserving algebraic fixture only; chosen vacuum "
                "reservoir covariance, not minimum physical noise or Caves bound"
            ),
            "units": {
                "J": "dimensionless",
                "K": "m^-1",
                "B": "m^-1/2",
                "J_f": "dimensionless",
                "D_vacuum_fixture": "m^-1",
                "commutator_integral": "dimensionless",
                "V_out": "dimensionless",
            },
            "J": _complex_matrix_payload(J),
            "K_per_m": _complex_matrix_payload(K),
            "K_eigenvalues_per_m": [float(value) for value in eigenvalues],
            "B_per_sqrt_m": _complex_matrix_payload(B),
            "J_f": _complex_matrix_payload(J_f),
            "D_vacuum_fixture_per_m": _complex_matrix_payload(D_vacuum_fixture),
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
            "fixture_unweighted": unweighted,
            "fixture_dc_balanced": dc_balanced,
            "ideal_bogoliubov_matched_to_G_s": {
                "S_linear": float(ideal),
                "S_dB": float(10 * np.log10(ideal)),
            },
            "symmetric_external_efficiency": float(eta_ext),
            "fixture_unweighted_after_external_loss": {
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
            "algebraic_fixture_detected_diagnostic_dB": float(
                10 * np.log10(detected)
            ),
            "fixture_diagnostic_difference_to_benchmark_dB": float(
                10 * np.log10(detected) - literature_squeezing_db
            ),
            "benchmark_bandwidth_MHz": literature_bandwidth_mhz,
            "model_bandwidth_MHz": None,
            "comparison_status": (
                "diagnostic only: algebraic dilation fixture is not a microscopic "
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
        "This is an analysis-only audit of `analysis/squeezing/analytic_reconstruction/ref_solver.py`; "
        "it does not modify the LaTeX report or production solver.",
        "",
        "**Classification warning:** the initial Floquet/velocity gain tables inherit the archived "
        "`ref_solver.py` propagation, which uses dressed optical wave numbers and "
        "a refractive phase mismatch together. They isolate Floquet/velocity "
        "numerics inside that shared implementation. Its dissipator also applies "
        "the inherited ground-coherence rate as coherence-only damping, without "
        "the production thermal transit reload channel. They are **not** predictions "
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
    provenance = option["parameter_provenance"]
    sensitivity = option["sensitivity"]
    dilation = option["algebraic_dilation_fixture"]
    benchmark = option["repository_literature_benchmark_comparison"]
    pump_reference = option["pump_only_weak_response_reference"]
    noncollinear = option["noncollinear_doppler_reference"]
    lines += [
        "",
        "## Corrected Option-A literature-point diagnostic",
        "",
        "This section is separate from the archived tables above. It uses bare, "
        "frequency-specific optical wave numbers in the susceptibility terms and "
        "only vacuum/geometric phase mismatch. No refractive-index contribution is "
        "inserted into the mismatch. The trace-normalized rho_ss supplies the "
        "manifold population once; the external structural factor is 1/[2(2I+1)].",
        "",
        f"Operating point: $\\Delta/2\\pi={op['D_GHz']:+.3f}$ GHz, "
        f"$\\delta/2\\pi={op['delta_MHz']:+.3f}$ MHz, "
        f"$T={op['temperature_C']:.1f}$ C, pump={1e3*op['P_pump']:.0f} mW, "
        f"seed={1e6*op['P_probe']:.0f} uW, $L={1e3*op['L']:.1f}$ mm, "
        f"$\\theta={op['theta_deg']:.2f}$ deg. The one-dimensional velocity grid "
        f"has {op['n_velocity']} points ($dv={op['velocity_step_m_per_s']:.1f}$ m/s, "
        f"cutoff {op['velocity_cutoff_sigma']:.0f} sigma). This base finite-seed "
        "table is one-dimensional; the separate slow two-dimensional reference "
        "below includes angular two-photon Doppler broadening.",
        "",
        "All rows below are evaluated at the same fixed "
        "$\\delta/2\\pi=-8$ MHz literature point; no detuning optimization is "
        "mixed into this table.",
        "",
        "| N_F | G_s | G_c power | G_c flux | flux gap | arg chi_sc (deg) | arg chi_cs (deg) |",
        "|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for row in option["floquet_convergence"]:
        lines.append(
            f"| {row['N_F']} | {_fmt(row['G_s'], 9)} | {_fmt(row['G_c'], 9)} | "
            f"{_fmt(row['G_c_photon_flux'], 9)} | "
            f"{_fmt(row['photon_flux_gap'], 9)} | {_fmt(row['arg_chi_sc_deg'], 9)} | "
            f"{_fmt(row['arg_chi_cs_deg'], 9)} |"
        )
    lines += [
        "",
        "### Weak-field reference-amplitude check",
        "",
        "The production atomic response is the current four-level finite-Floquet "
        "density-matrix model, and its finite seed/reference field enters the "
        "steady solve. The following "
        "N_F=3 check changes only that reference from 2 to 8 uW; it tests numerical "
        "weak-field linearity, not microscopic or experimental validity.",
        "",
        "| seed (uW) | Omega_s/2pi (MHz) | G_s | G_c power | flux gap | dG_s vs 2uW | dflux gap vs 2uW |",
        "|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for row in option["seed_reference_linearity"]:
        lines.append(
            f"| {row['probe_power_uW']:.0f} | "
            f"{row['probe_rabi_over_2pi_MHz']:.6f} | "
            f"{row['G_s']:.9f} | {row['G_c']:.9f} | "
            f"{row['photon_flux_gap']:.9f} | "
            f"{row['G_s_change_from_2uW_pct']:+.6f}% | "
            f"{row['photon_flux_gap_change_from_2uW_pct']:+.6f}% |"
        )
    pump_diag = pump_reference["pump_state_diagnostics"]
    frame = pump_reference["frame_equivalence"]
    seed_limit = pump_reference["finite_seed_to_infinitesimal"]
    dc_projection = pump_reference["trace_zero_dc_projection"]
    lines += [
        "",
        "### Pump-only stationary state and infinitesimal response reference",
        "",
        "This is a separate slow reference, not the production finite-seed gain "
        "path. For the standard minus Raman branch, the explicit rotating frame "
        "$U=\\exp(+i\\Omega_{\\rm beat}t|g_2\\rangle\\langle g_2|)$ makes both "
        "pump couplings static and gives $\\delta+\\Omega_{\\rm beat}=\\omega_{\\rm HF}$. "
        "The inherited plus-branch frame is not gauge-equivalent and is therefore "
        "reported as unsupported rather than silently mapped.",
        "",
        "| pump/reference check | value |",
        "|---|---:|",
        f"| max normalized pump residual | "
        f"{pump_diag['max_pump_normalized_residual']:.3e} |",
        f"| max pump trace error | {pump_diag['max_pump_trace_error']:.3e} |",
        f"| minimum pump-state eigenvalue | "
        f"{pump_diag['minimum_pump_eigenvalue']:.9e} |",
        f"| max trace-zero response residual | "
        f"{pump_diag['max_response_normalized_residual']:.3e} |",
        f"| max response trace error | "
        f"{pump_diag['max_response_trace_error']:.3e} |",
        f"| static-frame vs N_F=3 pump Floquet max difference | "
        f"{frame['max_abs_static_to_floquet_difference']:.3e} |",
        f"| projected DC relative frequency | "
        f"{dc_projection['relative_resolvent_frequency_rad_s']:.3e} rad/s |",
        f"| projected DC residual | "
        f"{dc_projection['max_response_normalized_residual']:.3e} |",
        "",
        "The finite-reference susceptibility approaches the analytic complex-field "
        "derivative as the reference amplitude is reduced. The table reports the "
        "empirical full-scan error; it does not impose the nominal quadratic law:",
        "",
        "| Rabi fraction | seed power (uW) | worst normalized chi error |",
        "|---:|---:|---:|",
    ]
    for row in seed_limit["rows"]:
        lines.append(
            f"| {row['rabi_fraction']:.3f} | {row['seed_power_uW']:.6f} | "
            f"{row['worst_normalized_chi_error']:.9e} |"
        )
    lines += [
        "",
        f"Fitted Rabi-amplitude error slope: "
        f"{seed_limit['rabi_error_power_law_slope']:.6f}; pairwise slopes "
        f"{', '.join(f'{value:.6f}' for value in seed_limit['pairwise_rabi_error_slopes'])}. "
        "The asymptotic perturbative expectation is 2, but the moving full-scan "
        "maximum is reported without relabelling it as exact quadratic convergence.",
        "",
        "| N_F at smallest reference amplitude | worst normalized chi error |",
        "|---:|---:|",
    ]
    for order, error in seed_limit["smallest_fraction_order_errors"].items():
        lines.append(f"| {order} | {error:.9e} |")
    lines += [
        "",
        "Direct resolvent and diagonalizable Liouvillian pole/residue expansion:",
        "",
        "| analysis frequency (MHz) | max direct/pole difference (s) | eigenvector condition | stationary residue max |",
        "|---:|---:|---:|---:|",
    ]
    for row in pump_reference["pole_residue_parity"]:
        lines.append(
            f"| {row['analysis_frequency_MHz']:.3f} | "
            f"{row['max_abs_direct_pole_difference_seconds']:.3e} | "
            f"{row['eigenvector_condition']:.6f} | "
            f"{row['stationary_residue_max']:.3e} |"
        )
    rms = noncollinear["rms_budget"]
    nominal_2d = noncollinear["nominal_literature_point"]
    lines += [
        "",
        "### Two-dimensional non-collinear Raman-Doppler reference",
        "",
        "This is another separate slow reference, not the production finite-seed "
        "scan. It tensor-averages the pump-only trace-zero response over "
        "independent Maxwellian $v_z$ and $v_x$. The laboratory optical beat and "
        "spectrum-analyzer frequency remain fixed while only the atomic detunings "
        "are velocity shifted:",
        "",
        "$\\Delta_{\\rm eff}=\\Delta-k_pv_z$,  "
        "$\\delta_{\\rm eff}=\\delta+(k_p-k_s\\cos\\theta)v_z"
        "$-k_s\\sin\\theta\\,v_x$.",
        "",
        "| Raman-Doppler width budget | value |",
        "|---|---:|",
        f"| angle | {rms['angle_deg']:.3f} deg |",
        f"| analytic total rms | {rms['analytic_total_MHz']:.9f} MHz |",
        f"| analytic transverse rms | {rms['analytic_transverse_MHz']:.9f} MHz |",
        f"| analytic axial rms at angle | {rms['analytic_axial_at_angle_kHz']:.6f} kHz |",
        f"| collinear residual rms | {rms['analytic_collinear_residual_kHz']:.6f} kHz |",
        f"| order-40 quadrature rms | {rms['quadrature_total_MHz']:.9f} MHz |",
        f"| quadrature error | {rms['quadrature_error_pct']:+.6f}% |",
        "",
        "Grid refinement uses the same 0.05 MHz refined feature scan and a "
        "five-sigma velocity cutoff:",
        "",
        "| order/axis | velocity pairs | feature delta (MHz) | feature G_s | feature G_c | G_s err. vs order 40 | shift vs order 40 | runtime (s) |",
        "|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for row in noncollinear["grid_refinement"]:
        lines.append(
            f"| {row['order_per_axis']} | {row['velocity_pairs']} | "
            f"{row['feature_delta_MHz']:.6f} | {row['feature_G_s']:.9f} | "
            f"{row['feature_G_c']:.9f} | "
            f"{row['feature_G_s_error_to_order40_pct']:+.6f}% | "
            f"{row['feature_position_error_to_order40_MHz']:+.6f} MHz | "
            f"{row['runtime_seconds']:.3f} |"
        )
    lines += [
        "",
        "Cutoff refinement uses order 40 per velocity axis:",
        "",
        "| cutoff (sigma) | feature delta (MHz) | feature G_s | G_s err. vs 5 sigma | shift vs 5 sigma |",
        "|---:|---:|---:|---:|---:|",
    ]
    for row in noncollinear["cutoff_refinement"]:
        lines.append(
            f"| {row['cutoff_sigma']:.1f} | {row['feature_delta_MHz']:.6f} | "
            f"{row['feature_G_s']:.9f} | "
            f"{row['feature_G_s_error_to_5sigma_pct']:+.6f}% | "
            f"{row['feature_position_error_to_5sigma_MHz']:+.6f} MHz |"
        )
    lines += [
        "",
        "| nominal -8 MHz comparison | value |",
        "|---|---:|",
        f"| one-dimensional G_s | {nominal_2d['one_dimensional_G_s']:.9f} |",
        f"| two-dimensional G_s | {nominal_2d['two_dimensional_G_s']:.9f} |",
        f"| G_s change | {nominal_2d['G_s_change_pct']:+.6f}% |",
        f"| max response normalized residual | {nominal_2d['max_response_normalized_residual']:.3e} |",
        f"| max response trace error | {nominal_2d['max_response_trace_error']:.3e} |",
        "",
        "All grid/cutoff acceptance flags pass. Remaining exclusions are the "
        "finite-seed two-dimensional production path, a beam angular "
        "distribution, segmentwise pump-state recomputation, and microscopic "
        "Langevin diffusion.",
    ]
    lines += [
        "",
        f"Vacuum/geometric mismatch: $\\Delta k_{{\\rm vac}}="
        f"{final['delta_k_vac_per_m']:.9f}\\,\\mathrm{{m^{{-1}}}}$.",
        f"Equal declared Gaussian collected-mode areas: "
        f"$A_p=A_c={final['mode_area_probe_m2']:.9e}\\,\\mathrm{{m^2}}$.",
    ]
    _append_complex_matrix_table(
        lines, "Q (E = Q a)", final["Q"], "[field per sqrt(photon/s)]")
    _append_complex_matrix_table(lines, "M_field", final["M_field_per_m"], "[m^-1]")
    _append_complex_matrix_table(lines, "T_field", final["T_field"], "[dimensionless]")
    _append_complex_matrix_table(
        lines, "M_canonical = Q^-1 M_field Q",
        final["M_canonical_per_m"], "[m^-1]")
    _append_complex_matrix_table(
        lines, "T_canonical = Q^-1 T_field Q",
        final["T_canonical"], "[dimensionless]")
    lines += [
        "",
        "The production Maxwell/Q path is independently reconstructed from the "
        "stored reduced susceptibilities using literal SI prefactors and "
        "`scipy.linalg.expm`:",
        "",
        "| independent reference comparison | max absolute difference |",
        "|---|---:|",
    ]
    for key, value in final["independent_reference_parity"].items():
        lines.append(f"| `{key}` | {value:.3e} |")

    lines += [
        "",
        "## Parameter provenance and illustrative sensitivity",
        "",
        "These sweeps are machine-readable under `parameter_provenance` and "
        "`sensitivity` in the JSON artifact. They are illustrative local tests, "
        "not parameter uncertainties or fits.",
        "",
        f"Current atomic construction: `{provenance['current_atomic_path']}`.",
        "",
        "| parameter | current value | source | status |",
        "|---|---:|---|---|",
    ]
    for name, entry in provenance["parameters"].items():
        lines.append(
            f"| `{name}` | {entry['current_value']:.9g} {entry['units']} | "
            f"`{entry['source']}` | {entry['status']} |"
        )
    lines += [
        "",
        "| table/block | atomic or numerical solver provenance |",
        "|---|---|",
    ]
    for name, entry in provenance["table_solver_provenance"].items():
        detail = "; ".join(f"{key}: {value}" for key, value in entry.items())
        lines.append(f"| `{name}` | {detail} |")

    ell = sensitivity["ell_s_propagation_only"]
    lines += [
        "",
        "### Propagation-only ell_s sweep",
        "",
        ell["classification"].capitalize() + ".",
        "",
        "| ell_s | effective line strength | G_s | G_c power | G_c flux | flux gap |",
        "|---:|---:|---:|---:|---:|---:|",
    ]
    for row in ell["rows"]:
        lines.append(
            f"| {row['ell_s']:.3f} | {row['effective_line_strength']:.9f} | "
            f"{row['G_s']:.9f} | {row['G_c_power']:.9f} | "
            f"{row['G_c_photon_flux']:.9f} | {row['photon_flux_gap']:.9f} |"
        )

    gamma = sensitivity["gamma_gg_floor_mean_field"]
    lines += [
        "",
        "### Ground-coherence-floor mean-field sensitivity",
        "",
        gamma["classification"].capitalize() + ".",
        "",
        "| gamma_gg transit floor / 2pi (kHz) | collision dephasing / 2pi (kHz) | total coherence rate / 2pi (kHz) | G_s | G_c power | G_c flux | flux gap | source |",
        "|---:|---:|---:|---:|---:|---:|---:|---|",
    ]
    for row in gamma["rows"]:
        lines.append(
            f"| {row['gamma_gg_floor_over_2pi_kHz']:.1f} | "
            f"{row['gamma_gg_collision_over_2pi_kHz']:.6f} | "
            f"{row['gamma_gg_total_over_2pi_kHz']:.6f} | "
            f"{row['G_s']:.9f} | {row['G_c_power']:.9f} | "
            f"{row['G_c_photon_flux']:.9f} | {row['photon_flux_gap']:.9f} | "
            f"{row['atomic_response_source']} |"
        )

    kappa = sensitivity["kappa_pump_scatter_diagnostic"]
    lines += [
        "",
        "### Pump-scatter kappa diagnostic",
        "",
        f"`{kappa['equation']}`. {kappa['classification']}. At "
        f"kappa={kappa['kappa_current']:.6g} and "
        f"OD_pump={kappa['OD_pump']:.9g}, "
        f"dN_ps/dkappa={kappa['dN_ps_dkappa']:.9g} and "
        f"N_ps={kappa['N_ps_current']:.9g}.",
    ]

    lines += [
        "",
        "## Algebraic commutator/diffusion fixture",
        "",
        "**This is an algebraic dilation fixture with one chosen vacuum-reservoir "
        "covariance.** It is not a Caves bound, microscopic atomic diffusion, or "
        "frequency dependent, and must not be presented as a Langevin-corrected "
        "squeezing-spectrum prediction. The canonical basis is constructed "
        "explicitly above from the declared frequencies and mode areas.",
        "",
        f"For $J=\\mathrm{{diag}}(1,-1)$, the eigenvalues of "
        f"$K=-(MJ+JM^\\dagger)$ are "
        f"{dilation['K_eigenvalues_per_m'][0]:.9f} and "
        f"{dilation['K_eigenvalues_per_m'][1]:.9f} m^-1. Both are positive, so "
        "$J_f=I_2$ for the displayed eigenfactor.",
    ]
    _append_complex_matrix_table(lines, "K", dilation["K_per_m"], "[m^-1]")
    _append_complex_matrix_table(lines, "B", dilation["B_per_sqrt_m"], "[m^-1/2]")
    _append_complex_matrix_table(
        lines,
        "D_vacuum_fixture = B B^dagger / 2",
        dilation["D_vacuum_fixture_per_m"],
        "[m^-1]",
    )
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
        f"| Algebraic vacuum fixture, unweighted | "
        f"{dilation['fixture_unweighted']['S_linear']:.7f} | "
        f"{dilation['fixture_unweighted']['S_dB']:.3f} | algebraic diagnostic |",
        f"| Algebraic vacuum fixture, DC-balanced | "
        f"{dilation['fixture_dc_balanced']['S_linear']:.7f} | "
        f"{dilation['fixture_dc_balanced']['S_dB']:.3f} | algebraic diagnostic |",
        f"| Ideal Bogoliubov matched to G_s | "
        f"{dilation['ideal_bogoliubov_matched_to_G_s']['S_linear']:.7f} | "
        f"{dilation['ideal_bogoliubov_matched_to_G_s']['S_dB']:.3f} | counterfactual benchmark |",
        f"| Algebraic fixture after external eta={dilation['symmetric_external_efficiency']:.4f} | "
        f"{dilation['fixture_unweighted_after_external_loss']['S_linear']:.7f} | "
        f"{dilation['fixture_unweighted_after_external_loss']['S_dB']:.3f} | external-loss diagnostic |",
        "",
        f"Against the repository literature benchmark (`README.md`, Sim et al. "
        f"85Rb optimum), the corrected Option-A "
        f"mean-field result is $G_s={benchmark['model_G_s']:.3f}$ versus "
        f"approximately {benchmark['benchmark_G_s']:.1f} "
        f"({benchmark['G_s_relative_discrepancy_pct']:.1f}%). The algebraic "
        f"fixture diagnostic after external loss is "
        f"{benchmark['algebraic_fixture_detected_diagnostic_dB']:.3f} dB versus "
        f"the reported scale near {benchmark['benchmark_squeezing_dB']:.1f} dB. "
        "No bandwidth comparison is available because the static model has no "
        "spectrum-analyzer frequency.",
    ]

    angular = audit["geometry_diagnostic"]
    lines += [
        "",
        "## Interpretation and API limitations",
        "",
        "- Production exposes `gabes.core.floquet_solve_truncated(...)` and an "
        "independently assembled dense-block reference for arbitrary finite N_F. "
        "`gabes.kernels.floquet_chi_grid(...)` and "
        "`gabes.schemes.fwm.chi_matrix_table(...)` carry the same order argument; "
        "tests pin compiled/reference parity at N_F=1,2,3.",
        "- `gabes.schemes.fwm.compute_spectrum(...)` defaults to N_F=3 and compares "
        "the entire reported scan with N_F=2 using complex-response/transfer, gain, "
        "wrapped-phase, and optimum-shift criteria. It also exposes `velocity_step` "
        "and `velocity_cutoff`. Its one-dimensional path calls "
        "`gabes.doppler.velocity_grid(...)`, `build_Delta_eff_axis(...)`, and "
        "`doppler_average(...)`.",
        "- The archived N_F=1 result is not Floquet-converged: N_F=2 changes the "
        "common-point gains substantially. N_F=2 and N_F=3 agree to the precision "
        "reported here; production no longer promotes a single-point comparison.",
        "- The production velocity refinement still converges only the collinear "
        "integral $\\Delta_{eff}=\\Delta-kv$. The opt-in two-dimensional reference "
        "above represents the crossing-angle Raman-Doppler distribution without "
        "changing the laboratory beat frequency.",
        f"- At theta={angular['theta_deg']:.2f} deg, sigma_v="
        f"{angular['sigma_v_m_per_s']:.3f} m/s, and lambda="
        f"{angular['lambda_nm']:.3f} nm, the one-sigma angular width is "
        f"{angular['angular_two_photon_sigma_mhz']:.3f} MHz; the separate reference "
        "resolves this width explicitly.",
        "- The two-dimensional reference separates the lab beat frequency from "
        "the velocity-shifted atomic two-photon detuning. Reusing the production "
        "`floquet_chi_grid` with `delta_eff` remains invalid because that kernel "
        "also derives `omega_beat = omega_hf + branch*delta` from the same value.",
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
        "scope": (
            "Composite analysis artifact: archived reduced-model Floquet/velocity "
            "isolation plus a separate current Option-A literature-point diagnostic"
        ),
        "archived_reduced_model_classification": {
            "implementation": "analysis/squeezing/analytic_reconstruction/ref_solver.py",
            "inherits_dressed_k_plus_refractive_mismatch": True,
            "corrected_option_a": False,
            "atomic_dissipator": (
                "coherence-only inherited gamma_gg; no trace-preserving thermal "
                "transit reload"
            ),
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
