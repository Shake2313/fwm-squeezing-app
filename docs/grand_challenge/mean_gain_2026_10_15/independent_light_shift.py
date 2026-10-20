"""Independent four-level pump-Hamiltonian / Raman-response audit (AI supplement).

No production Hamiltonian, commutator, stationary solver, response solver or
Floquet implementation is used to form the reference. Production imports are
deferred to the comparison. This is a cheap deterministic local-atom control,
not a Doppler/gain prediction or a full-Zeeman reduction.
"""
from pathlib import Path
import json
import numpy as np

from gabes import constants as c, hyperfine

UNIT = 2 * np.pi * 1e6
POWER_W = 0.2948
WAIST_M = 530e-6
TEMPERATURE_K = 386.15
CF = np.array([[10, 35], [35, 28]], dtype=float) / 81
SCALE = np.sqrt(3 * CF)
OMEGA = c.GAMMA * np.sqrt(POWER_W / (np.pi * WAIST_M**2 * c.I_SAT))
HF = c.OMEGA_HF
EHF = c.OMEGA_EXCITED_HF


def hamiltonian(delta_eff, a=OMEGA, b=OMEGA):
    h = np.diag([0, HF, -delta_eff - EHF, -delta_eff]).astype(complex)
    for g, drive in enumerate((a, b)):
        for j, e in enumerate((2, 3)):
            h[g, e] = h[e, g] = drive * SCALE[g, j] / 2
    return h


def generator(h):
    density = hyperfine.number_density(TEMPERATURE_K)
    optical_extra = hyperfine.BETA_SELF * density / 2
    ground_extra = density * c.RB_GROUND_COHERENCE_XSECTION * c.mean_relative_speed(TEMPERATURE_K)
    reload_rate = c.GAMMA_GG
    thermal = np.diag([5 / 12, 7 / 12, 0, 0])
    jumps = []
    for j, e in enumerate((2, 3)):
        for g in (0, 1):
            jump = np.zeros((4, 4), complex)
            jump[g, e] = np.sqrt(c.GAMMA * CF[g, j] / CF[:, j].sum())
            jumps.append(jump)

    def action(rho):
        out = -1j * (h @ rho - rho @ h)
        for jump in jumps:
            square = jump.conj().T @ jump
            out += jump @ rho @ jump.conj().T - (square @ rho + rho @ square) / 2
        out[0, 1] -= ground_extra * rho[0, 1]
        out[1, 0] -= ground_extra * rho[1, 0]
        for g in (0, 1):
            for e in (2, 3):
                out[g, e] -= optical_extra * rho[g, e]
                out[e, g] -= optical_extra * rho[e, g]
        out += reload_rate * (thermal * np.trace(rho) - rho)
        return out

    basis = np.eye(16).reshape(16, 4, 4)
    return np.stack([action(rho).reshape(16) for rho in basis], axis=1)


def state(l):
    matrix, rhs = l.copy(), np.zeros(16, complex)
    matrix[0] = np.eye(4).reshape(16)
    rhs[0] = 1
    return np.linalg.solve(matrix, rhs).reshape(4, 4)


def trace_zero_basis():
    columns = []
    for j in range(16):
        if j in (0, 5, 10, 15):
            continue
        column = np.zeros(16)
        column[j] = 1
        columns.append(column)
    for j in (0, 5, 10):
        column = np.zeros(16)
        column[j], column[15] = 1, -1
        columns.append(column)
    return np.linalg.qr(np.stack(columns, axis=1))[0]


PROJECT = trace_zero_basis()


def weak_response(l, rho, delta):
    probe, conjugate_star = np.zeros((4, 4), complex), np.zeros((4, 4), complex)
    for j, e in enumerate((2, 3)):
        probe[e, 1] = SCALE[1, j] / 2
        conjugate_star[0, e] = SCALE[0, j] / 2
    sources = np.stack([-1j * (q @ rho - rho @ q).reshape(16)
                        for q in (probe, conjugate_star)], axis=1)
    shifted = l + 1j * (-HF + delta) * np.eye(16)
    response = PROJECT @ np.linalg.solve(PROJECT.T @ shifted @ PROJECT, -PROJECT.T @ sources)
    matrices = response.T.reshape(2, 4, 4)
    outputs = np.stack([sum(SCALE[1, j] * matrices[:, e, 1] for j, e in enumerate((2, 3))),
                        sum(SCALE[0, j] * matrices[:, 0, e] for j, e in enumerate((2, 3)))])
    residual = np.linalg.norm(shifted @ response + sources) / (
        np.linalg.norm(shifted @ response) + np.linalg.norm(sources))
    return outputs, float(residual)


def dressed_and_pole(delta_eff):
    h = hamiltonian(delta_eff)
    eigenvalues, eigenvectors = np.linalg.eigh(h)
    ground_indices = np.argmax(abs(eigenvectors[:2])**2, axis=1)
    gap = eigenvalues[ground_indices[1]] - eigenvalues[ground_indices[0]]
    dressed = (HF - gap) / UNIT
    l = generator(h)
    ev = np.linalg.eigvals(l)
    poles = (HF - ev.imag) / UNIT
    closest = int(np.argmin(abs(poles - dressed)))
    second = [OMEGA**2 / 4 * sum(SCALE[g, j]**2 / (
        delta_eff + joffset + g * HF) for j, joffset in enumerate((EHF, 0))) / UNIT
        for g in (0, 1)]
    return {
        "delta_eff_MHz": float(delta_eff / UNIT),
        "second_order_ground_shifts_MHz": second,
        "second_order_differential_shift_MHz": second[0] - second[1],
        "exact_ground_energy_shifts_MHz": [float(eigenvalues[ground_indices[0]] / UNIT),
                                            float((eigenvalues[ground_indices[1]] - HF) / UNIT)],
        "exact_dressed_Raman_delta_MHz": float(dressed),
        "dissipative_Raman_pole_delta_MHz": float(poles[closest]),
        "dissipative_Raman_HWHM_MHz": float(-ev[closest].real / UNIT),
    }


def signed_zeeman_pump_control():
    # Established independent I x J reference, no production CG normalization.
    # Illustrative linear pump / zero B; this is not a full-atom FWM response.
    from analysis.grand_challenge.reference.uncoupled_d1 import d1_dipole_reference
    reference = d1_dipole_reference()
    labels = reference["labels"]
    peak_field = np.sqrt(4 * POWER_W / (np.pi * WAIST_M**2 * c.EPS_0 * c.C_LIGHT))
    omega_j = c.DIPOLE_D1 * peak_field / c.HBAR
    ground = [0 if f == 2 else HF for f, m in labels]
    excited = [-1000 * UNIT - (EHF if f == 2 else 0) for f, m in labels]
    h = np.diag(ground + excited).astype(complex)
    h[12:, :12] = omega_j * reference["dipoles"][1] / 2
    h[:12, 12:] = h[12:, :12].conj().T
    values, vectors = np.linalg.eigh(h)
    indices = np.argmax(abs(vectors[:12])**2, axis=1)
    rows = []
    for m in range(-2, 3):
        a, b = labels.index((2, m)), labels.index((3, m))
        rows.append({"m": m, "same_m_dressed_Raman_delta_MHz": float(
            (HF - values[indices[b]] + values[indices[a]]) / UNIT)})
    return {
        "scope": "Illustrative exact 24-level Hamiltonian, linear q=0 pump, B=0, v=0; no dissipator or probe polarization/gain solve.",
        "adoption_status": "diagnostic_only_not_adopted",
        "origin": "analysis/grand_challenge/reference/uncoupled_d1.py",
        "rows": rows,
    }


def main():
    # Reference assembled before loading the production scheme.
    deltas_MHz = np.array([-10, 0, 2.2, 23.7, 52.269382386541935, 70.])
    deff_MHz = np.array([750., 1000., 1250.])
    reference = np.empty((len(deltas_MHz), len(deff_MHz), 2, 2), complex)
    states, generators, residuals = [], [], []
    for k, deff in enumerate(deff_MHz):
        l = generator(hamiltonian(deff * UNIT))
        rho = state(l)
        states.append(rho)
        generators.append(l)
        for j, delta in enumerate(deltas_MHz):
            reference[j, k], residual = weak_response(l, rho, delta * UNIT)
            residuals.append(residual)

    from gabes.schemes import fwm
    production_reference = fwm.pump_only_weak_response_reference(
        OMEGA, OMEGA, deltas_MHz * UNIT, deff_MHz * UNIT, T=TEMPERATURE_K)
    actual_atom = fwm.collisional_atom(TEMPERATURE_K)
    expected_tuple = (reference[..., 0, 0], np.conj(reference[..., 1, 0]),
                      reference[..., 0, 1], np.conj(reference[..., 1, 1]))
    all_norm = np.linalg.norm(reference)
    finite = []
    for power_uw in (3.06 / (20.3063457 / 26), 0.1, 0.01, 0.001, 0.0001):
        seed = fwm.rabi_freq(power_uw * 1e-6, 330e-6)
        got = fwm.chi_matrix_table(OMEGA, OMEGA, seed, seed,
                                  deltas_MHz * UNIT, deff_MHz * UNIT,
                                  -1, atom=actual_atom, n_f=3)
        error = np.stack([a - b for a, b in zip(got, expected_tuple)])
        finite.append({"reference_seed_power_uW": power_uw,
                       "relative_Frobenius_error": float(np.linalg.norm(error) / all_norm),
                       "max_absolute_error_seconds": float(abs(error).max()),
                       "per_detuning_relative_Frobenius_error": [
                           float(np.linalg.norm(error[:, j]) / np.linalg.norm(reference[j]))
                           for j in range(len(deltas_MHz))]})
    manual_L_error = max(np.linalg.norm(manual - (
        fwm.build_liouvillian(fwm.pump_hamiltonian_at_Deff_zero(OMEGA, OMEGA), actual_atom)
        - deff * UNIT * actual_atom.S_v)) / np.linalg.norm(manual)
        for manual, deff in zip(generators, deff_MHz))
    report = {
        "provenance": "⟦AI⟧ Independent local-atom reference, not a measured/full-atom gain model.",
        "inputs": {"pump_power_W": POWER_W, "pump_waist_m": WAIST_M,
                   "temperature_K": TEMPERATURE_K, "pump_Rabi_over_2pi_MHz": OMEGA / UNIT,
                   "delta_axis_MHz": deltas_MHz.tolist(), "delta_eff_axis_MHz": deff_MHz.tolist(),
                   "branch": -1, "transit_rate_rad_s": c.GAMMA_GG},
        "stark_controls": [dict(velocity_mps=v, **dressed_and_pole(1000 * UNIT - c.K_VEC * v))
                           for v in (-400., -200., 0., 200., 400.)],
        "manual_L_relative_error_vs_production": float(manual_L_error),
        "pump_state_max_error_vs_production": float(abs(np.array(states) - production_reference.pump_state).max()),
        "manual_response_max_normalized_residual": max(residuals),
        "response_relative_error_vs_production_reference": float(np.linalg.norm(
            reference - production_reference.chi_matrix[0]) / all_norm),
        "finite_seed_comparisons": finite,
        "signed_zeeman_pump_control": signed_zeeman_pump_control(),
    }
    target = Path(__file__).with_name("independent_light_shift.json")
    target.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2, ensure_ascii=True))


if __name__ == "__main__":
    main()
