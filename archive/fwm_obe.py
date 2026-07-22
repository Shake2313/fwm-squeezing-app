"""
Backward-compatibility shim for the FWM backend.

The physics moved into the `gabes/` package (Phase-0 refactor). This module
re-exports the FWM scheme's public API under the old `fwm_obe` names so existing
imports and the `python fwm_obe.py` CLI keep working unchanged. New code should
import from `gabes` / `gabes.schemes.fwm` directly.

CLI:  python fwm_obe.py   →  prints a summary and saves the gain/squeezing figure.
"""
from pathlib import Path
import os

os.environ.setdefault("OMP_NUM_THREADS", str(os.cpu_count() or 1))
os.environ.setdefault("OPENBLAS_NUM_THREADS", str(os.cpu_count() or 1))
os.environ.setdefault("MKL_NUM_THREADS", str(os.cpu_count() or 1))

import numpy as np

from gabes import atoms, constants
from gabes.constants import (  # noqa: F401  (re-exported for compatibility)
    HBAR, KB, C_LIGHT, EPS_0, ELEMENTARY_CHARGE,
    NU_D1_85RB, WAVELENGTH_D1_85RB, GAMMA_2PI, NU_GROUND_HF, NU_EXCITED_HF_D1,
    NU_HF, MASS_85RB, I_SAT, DIPOLE_D1, RB85_ABUNDANCE, LINE_STRENGTH_FACTOR,
    GAMMA, OMEGA_HF, OMEGA_EXCITED_HF, K_VEC, OMEGA_D1,
    GAMMA_GG_2PI, GAMMA_GG, rabi_freq,
)
from gabes.atoms import rb85_density
from gabes.core import comm_super, build_liouvillian, floquet_solve, matrix_exp_2x2
from gabes.doppler import velocity_grid, build_Delta_eff_axis, doppler_average
from gabes.observables import gain_from_chi, intensity_difference_squeezing_dB
from gabes.schemes.fwm import (  # noqa: F401
    G1, G2, E2, E3, GROUND_STATES, EXCITED_STATES, N_LEVELS,
    L_CELL, W_PUMP, W_PROBE, P_PUMP, P_PROBE, T_CELL,
    QE_DETECTOR, RESPONSIVITY_AW, LOSS_FRAC, ETA_TOTAL,
    OMEGA_C_SEED, SCAN_MIN_GHZ, SCAN_MAX_GHZ, SCAN_COARSE_POINTS,
    RESONANCE_WINDOW_MHZ, SCAN_FINE_POINTS, PUMP_OVERLAP_EXCLUSION_MHZ,
    VELOCITY_STEP_MPS, VELOCITY_CUTOFF_SIGMA, DELTA_GHZ_LIST, BRANCHES,
    static_hamiltonian_at_Deff_zero, sideband_hamiltonian, sideband_template,
    chi_matrix_table, branch_center_GHz, probe_scan_axis_GHz,
    two_photon_detuning_from_probe_scan, compute_spectrum, operating_point,
    ATOM,
)

SCRIPT_DIR = Path(__file__).resolve().parent

# Legacy aliases for the 4-level globals the old module exposed at top level.
L_LINDBLAD = ATOM.lindblad
S_V_SUPER = ATOM.S_v
SIGMA_V = np.sqrt(KB * T_CELL / MASS_85RB)


def L0_of(H0):
    """L₀ = comm_super(H0) + Lindblad (old top-level helper, 4-level model)."""
    return build_liouvillian(H0, ATOM)


def rho_index(row, col):
    return row * N_LEVELS + col


def main():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    Op_A = rabi_freq(P_PUMP, W_PUMP)
    Os = rabi_freq(P_PROBE, W_PROBE)
    N_atoms = rb85_density(T_CELL)

    print(f"T                = {T_CELL:7.2f} K")
    print(f"N(85Rb)          = {N_atoms:.3e} /m^3")
    print(f"Omega_pA / 2pi   = {Op_A /(2*np.pi)/1e9:7.3f} GHz")
    print(f"Omega_s  / 2pi   = {Os  /(2*np.pi)/1e6:7.3f} MHz")
    print(f"sigma_v          = {SIGMA_V:7.1f} m/s")
    print(f"L_cell           = {L_CELL*1e3:7.3f} mm")
    print(f"eta_total        = {ETA_TOTAL:7.4f} (QE {QE_DETECTOR:.4f} × (1-loss {LOSS_FRAC:.3f}))")

    v_grid, _ = velocity_grid(T_CELL, dv=VELOCITY_STEP_MPS, cutoff_sigma=VELOCITY_CUTOFF_SIGMA)
    dv = v_grid[1] - v_grid[0]
    print(f"dv               = {dv:7.1f} m/s  (N_v = {len(v_grid)})")

    results = {}
    for D_GHz in DELTA_GHZ_LIST:
        spec = compute_spectrum(D_GHz, T=T_CELL)
        results[D_GHz] = (spec["probe_axis_GHz"], spec["G_s"],
                          spec["G_c"], spec["S_dB"])
        print(f"Delta / 2pi = {D_GHz:.1f} GHz done"
              f"  max G_s = {spec['G_s'].max():.3f}"
              f"  best squeezing = {spec['S_dB'].min():.2f} dB")

    fig, (axG, axS) = plt.subplots(2, 1, figsize=(8.0, 7.2), sharex=True)
    use_log = any(np.nanmax(results[d][1]) > 50 for d in DELTA_GHZ_LIST)
    for D_GHz in DELTA_GHZ_LIST:
        probe_axis_GHz, G_s, _G_c, S_dB = results[D_GHz]
        label = f"Δ/2π = {D_GHz:.1f} GHz"
        axG.plot(probe_axis_GHz, G_s, label=label)
        axS.plot(probe_axis_GHz, S_dB, label=label)

    for ax in (axG, axS):
        ax.axvline(0, color="gray", ls=":", lw=0.8)
        ax.grid(alpha=0.3)
        ax.set_xlim(SCAN_MIN_GHZ, SCAN_MAX_GHZ)
    if use_log:
        axG.set_yscale("log")
    axG.set_ylabel("Probe (seed) power gain $G_s$")
    axG.set_title("Probe / seed gain")
    axG.axhline(1.0, color="black", lw=0.6)
    axG.legend(fontsize=9)
    axS.set_ylabel("Intensity-difference squeezing [dB]")
    axS.axhline(0.0, color="black", lw=0.6)
    axS.set_title(f"Squeezing (η = {ETA_TOTAL:.4f}, "
                  f"line-strength factor = {LINE_STRENGTH_FACTOR})")
    axS.set_xlabel(r"Probe detuning from $F=2 \rightarrow F'=3$   [GHz]")
    axS.legend(fontsize=9)
    fig.tight_layout()

    out = SCRIPT_DIR / "fwm_double_lambda_gain_squeezing.png"
    fig.savefig(out, dpi=130)
    print(f"Saved figure: {out}")


if __name__ == "__main__":
    main()
