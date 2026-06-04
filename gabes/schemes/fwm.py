"""
Cluster D — 85Rb D1 double-Λ four-wave mixing.

The physics is ported verbatim from the original fwm_obe.py (Hamiltonians,
χ̄-matrix table, probe scan, Doppler average, Maxwell-Bloch propagation). The
only change is that shared engine pieces now come from gabes.core / .doppler /
.observables and the 4-level structure from gabes.atoms.

`compute_spectrum` / `operating_point` keep their original signatures so the CLI
shim and the regression test can call them directly. `FWMScheme` wraps them for
the generic Streamlit front-end.
"""
import numpy as np

from .. import atoms, constants, doppler, hyperfine, observables
from ..constants import K_VEC, OMEGA_HF, OMEGA_EXCITED_HF, rabi_freq
from ..core import build_liouvillian, comm_super, floquet_solve
from .base import ExtraView, ParamSpec, Preset, Scheme

# =========================================================
# Level indices (this scheme's labelling of the atom model)
# =========================================================
ATOM = atoms.get("double_lambda_rb85")
G1, G2, E2, E3 = 0, 1, 2, 3
GROUND_STATES = (G1, G2)
EXCITED_STATES = (E2, E3)
N_LEVELS = ATOM.n_levels
GROUND_F = {G1: 2, G2: 3}
EXCITED_F = {E2: 2, E3: 3}
TRANSITION_DIPOLE_SCALE = np.zeros((N_LEVELS, N_LEVELS), dtype=float)
for _g in GROUND_STATES:
    for _e in EXCITED_STATES:
        TRANSITION_DIPOLE_SCALE[_g, _e] = np.sqrt(
            3.0 * hyperfine.CF2[(GROUND_F[_g], EXCITED_F[_e])])

# =========================================================
# FWM experiment configuration (cell, beams, detection, scan)
# =========================================================
L_CELL = 12.5e-3
W_PUMP = 530e-6
W_PROBE = 330e-6
P_PUMP, P_PROBE = 600e-3, 10e-6
T_CELL = 394.15

QE_DETECTOR = 0.9047
RESPONSIVITY_AW = 0.58
LOSS_FRAC = 0.0
ETA_TOTAL = QE_DETECTOR * (1.0 - LOSS_FRAC)

OMEGA_C_SEED = 0.0
SCAN_MIN_GHZ = -8.0
SCAN_MAX_GHZ = 12.0
SCAN_COARSE_POINTS = 401
RESONANCE_WINDOW_MHZ = 80.0
SCAN_FINE_POINTS = 801
PUMP_OVERLAP_EXCLUSION_MHZ = 1e-3
VELOCITY_STEP_MPS = 1.0
VELOCITY_CUTOFF_SIGMA = 3.0
DELTA_GHZ_LIST = [0.9]
BRANCHES = (-1, +1)
DEFAULT_BRANCH = -1


# =========================================================
# Hamiltonians
# =========================================================
def _add_static_drive(H, ground, omega):
    for excited in EXCITED_STATES:
        omega_ge = omega * TRANSITION_DIPOLE_SCALE[ground, excited]
        H[ground, excited] += omega_ge / 2
        H[excited, ground] += omega_ge / 2


def _add_sideband_drive(H, ground, omega):
    for excited in EXCITED_STATES:
        H[excited, ground] += omega * TRANSITION_DIPOLE_SCALE[ground, excited] / 2


def _polarization_coherence(rho, ground):
    return sum(TRANSITION_DIPOLE_SCALE[ground, e] * rho[:, e, ground]
               for e in EXCITED_STATES)


def static_hamiltonian_at_Deff_zero(Op_A, Op_B, Os, delta, branch):
    """H₀ with Δ_eff = 0, so the only v / Δ_eff dependence is added later."""
    H0 = np.zeros((N_LEVELS, N_LEVELS), dtype=complex)
    H0[G2, G2] = delta
    H0[E2, E2] = -OMEGA_EXCITED_HF
    H0[E3, E3] = 0.0
    if branch == -1:
        _add_static_drive(H0, G1, Op_A)
        _add_static_drive(H0, G2, Os)
        return H0
    if branch == +1:
        _add_static_drive(H0, G1, Os)
        _add_static_drive(H0, G2, Op_B)
        return H0
    raise ValueError(f"branch must be one of {BRANCHES}, got {branch}")


def sideband_hamiltonian(Op_A, Op_B, Oc, branch):
    Hp = np.zeros((N_LEVELS, N_LEVELS), dtype=complex)
    if branch == -1:
        _add_sideband_drive(Hp, G1, Oc)
        _add_sideband_drive(Hp, G2, Op_B)
        return Hp
    if branch == +1:
        _add_sideband_drive(Hp, G1, Op_A)
        _add_sideband_drive(Hp, G2, Oc)
        return Hp
    raise ValueError(f"branch must be one of {BRANCHES}, got {branch}")


def sideband_template(Op_A, Op_B, Oc, branch):
    Hp = sideband_hamiltonian(Op_A, Op_B, Oc, branch)
    Cp = comm_super(Hp)
    Cm = comm_super(Hp.conj().T)
    return Cp, Cm


# =========================================================
# χ-matrix table   (T- and Δ-independent)
# =========================================================
def chi_matrix_table(Op_A, Op_B, Os_ref, Oc_ref, delta_axis, Delta_eff_axis, branch):
    """
    Two solves per probe-detuning point to extract (χ̄_ss, χ̄_cs, χ̄_sc, χ̄_cc)
    on a 2-D (δ, Δ_eff) grid. Returns each array (n_delta, n_deff), complex.
    """
    probe_ground = G2 if branch == -1 else G1
    conj_ground = G1 if branch == -1 else G2
    n_d = delta_axis.size
    n_de = Delta_eff_axis.size

    chi_ss = np.zeros((n_d, n_de), dtype=complex)
    chi_cs = np.zeros((n_d, n_de), dtype=complex)
    chi_sc = np.zeros((n_d, n_de), dtype=complex)
    chi_cc = np.zeros((n_d, n_de), dtype=complex)

    Cp_no_c, Cm_no_c = sideband_template(Op_A, Op_B, 0.0, branch)      # solve 1
    Cp_c, Cm_c = sideband_template(Op_A, Op_B, Oc_ref, branch)        # solve 2

    for i, delta in enumerate(delta_axis):
        Omega_beat = OMEGA_HF - branch * delta

        # ---- Solve 1: probe drive only ----
        H0_1 = static_hamiltonian_at_Deff_zero(Op_A, Op_B, Os_ref, delta, branch)
        L0_1 = build_liouvillian(H0_1, ATOM)
        rho0_a, rhop_a = floquet_solve(
            L0_1, Cp_no_c, Cm_no_c, Omega_beat, Delta_eff_axis, ATOM.S_v, N_LEVELS)
        probe_a = _polarization_coherence(rho0_a, probe_ground)
        conj_a = _polarization_coherence(rhop_a, conj_ground)
        chi_ss[i] = probe_a / Os_ref
        chi_cs[i] = conj_a / Os_ref

        # ---- Solve 2: conjugate seed only ----
        H0_2 = static_hamiltonian_at_Deff_zero(Op_A, Op_B, 0.0, delta, branch)
        L0_2 = build_liouvillian(H0_2, ATOM)
        rho0_b, rhop_b = floquet_solve(
            L0_2, Cp_c, Cm_c, Omega_beat, Delta_eff_axis, ATOM.S_v, N_LEVELS)
        probe_b = _polarization_coherence(rho0_b, probe_ground)
        conj_b = _polarization_coherence(rhop_b, conj_ground)
        chi_sc[i] = probe_b / Oc_ref
        chi_cc[i] = conj_b / Oc_ref

    return chi_ss, chi_cs, chi_sc, chi_cc


# =========================================================
# Probe-detuning axis
# =========================================================
def branch_center_GHz(Delta_GHz, branch):
    if branch not in BRANCHES:
        raise ValueError(f"branch must be one of {BRANCHES}, got {branch}")
    return Delta_GHz + branch * constants.NU_HF / 1e9


def probe_scan_axis_GHz(Delta_GHz, coarse_points=None, fine_points=None,
                        window_mhz=None, scan_min=None, scan_max=None,
                        branches=BRANCHES):
    """Probe-frequency scan axis [GHz], referenced to the 85Rb F=2 → F'=3 line."""
    coarse_points = SCAN_COARSE_POINTS if coarse_points is None else coarse_points
    fine_points = SCAN_FINE_POINTS if fine_points is None else fine_points
    window_mhz = RESONANCE_WINDOW_MHZ if window_mhz is None else window_mhz
    scan_min = SCAN_MIN_GHZ if scan_min is None else scan_min
    scan_max = SCAN_MAX_GHZ if scan_max is None else scan_max

    coarse = np.linspace(scan_min, scan_max, coarse_points)
    half_window = window_mhz * 1e-3
    parts = [coarse]
    for branch in branches:
        center = branch_center_GHz(Delta_GHz, branch)
        fmin = max(scan_min, center - half_window)
        fmax = min(scan_max, center + half_window)
        if fmin < fmax:
            parts.append(np.linspace(fmin, fmax, fine_points))
    axis = np.unique(np.concatenate(parts))
    exclusion_GHz = PUMP_OVERLAP_EXCLUSION_MHZ * 1e-3
    return axis[np.abs(axis - Delta_GHz) > exclusion_GHz]


def two_photon_detuning_from_probe_scan(probe_GHz, Delta_GHz, branch):
    delta_Hz = (probe_GHz - branch_center_GHz(Delta_GHz, branch)) * 1e9
    return 2 * np.pi * delta_Hz


def _single_branch(branch, branches):
    """Resolve one physical FWM channel; do not merge distinct Raman branches."""
    if branches is None:
        if branch not in BRANCHES:
            raise ValueError(f"branch must be one of {BRANCHES}, got {branch}")
        return branch

    branches = tuple(branches)
    if len(branches) != 1:
        raise ValueError(
            "FWM Raman branches are separate probe/conjugate mode pairs; "
            "compute them one at a time instead of summing susceptibilities."
        )
    only = branches[0]
    if only not in BRANCHES:
        raise ValueError(f"branch must be one of {BRANCHES}, got {only}")
    return only


# =========================================================
# High-level spectrum  (one call → gain + squeezing curves)
# =========================================================
def compute_spectrum(D_GHz, *,
                     T=T_CELL, P_pump=P_PUMP, P_probe=P_PROBE,
                     w_pump=W_PUMP, w_probe=W_PROBE,
                     line_strength=None, loss_frac=LOSS_FRAC, qe=QE_DETECTOR,
                     coarse_points=None, fine_points=None, window_mhz=None,
                     scan_min=None, scan_max=None,
                     velocity_step=None, velocity_cutoff=None,
                     branch=DEFAULT_BRANCH, branches=None):
    """Full pipeline for one one-photon detuning Δ = 2π·D_GHz·1e9 (see README)."""
    branch = _single_branch(branch, branches)
    if line_strength is None:
        line_strength = constants.LINE_STRENGTH_FACTOR
    eta = qe * (1.0 - loss_frac)

    Op_A = rabi_freq(P_pump, w_pump)
    Op_B = Op_A
    Os = rabi_freq(P_probe, w_probe)
    Os_ref = Os
    Oc_ref = Os                              # χ̄ is independent of |Ω_ref|

    N_atoms = atoms.rb85_density(T)
    Delta = 2 * np.pi * D_GHz * 1e9

    probe_axis_GHz = probe_scan_axis_GHz(
        D_GHz, coarse_points, fine_points, window_mhz, scan_min, scan_max,
        branches=(branch,))
    velocity_step = VELOCITY_STEP_MPS if velocity_step is None else velocity_step
    velocity_cutoff = VELOCITY_CUTOFF_SIGMA if velocity_cutoff is None else velocity_cutoff
    v_grid, weights = doppler.velocity_grid(
        T, dv=velocity_step, cutoff_sigma=velocity_cutoff)
    Delta_eff_axis = doppler.build_Delta_eff_axis(Delta, Delta, v_grid)

    chi_ss_avg = np.zeros(probe_axis_GHz.size, dtype=complex)
    chi_cs_avg = np.zeros(probe_axis_GHz.size, dtype=complex)
    chi_sc_avg = np.zeros(probe_axis_GHz.size, dtype=complex)
    chi_cc_avg = np.zeros(probe_axis_GHz.size, dtype=complex)

    delta_axis = two_photon_detuning_from_probe_scan(probe_axis_GHz, D_GHz, branch)
    ch_ss, ch_cs, ch_sc, ch_cc = chi_matrix_table(
        Op_A, Op_B, Os_ref, Oc_ref, delta_axis, Delta_eff_axis, branch)
    chi_ss_avg += doppler.doppler_average(ch_ss, Delta_eff_axis, Delta, v_grid, weights)
    chi_cs_avg += doppler.doppler_average(ch_cs, Delta_eff_axis, Delta, v_grid, weights)
    chi_sc_avg += doppler.doppler_average(ch_sc, Delta_eff_axis, Delta, v_grid, weights)
    chi_cc_avg += doppler.doppler_average(ch_cc, Delta_eff_axis, Delta, v_grid, weights)

    G_s, G_c, _ = observables.gain_from_chi(
        chi_ss_avg, chi_sc_avg, chi_cs_avg, chi_cc_avg,
        K_VEC, K_VEC, L_CELL, N_atoms, line_strength=line_strength)
    S_dB = observables.intensity_difference_squeezing_dB(G_s, G_c, eta)

    return {
        "D_GHz": D_GHz,
        "probe_axis_GHz": probe_axis_GHz,
        "G_s": G_s,
        "G_c": G_c,
        "S_dB": S_dB,
        "eta": eta,
        "branch": branch,
        "N_atoms": N_atoms,
        "sigma_v": np.sqrt(constants.KB * T / constants.MASS_85RB),
        "n_velocity": v_grid.size,
        "Op_A_2pi_GHz": Op_A / (2 * np.pi) / 1e9,
        "Os_2pi_MHz": Os / (2 * np.pi) / 1e6,
        "raman_center_minus_GHz": branch_center_GHz(D_GHz, -1),
        "raman_center_plus_GHz": branch_center_GHz(D_GHz, +1),
    }


def operating_point(spectrum, delta_mhz, branch=-1):
    """Read G_s / G_c / squeezing at a chosen δ (MHz) on the selected branch."""
    probe_GHz = (spectrum["raman_center_minus_GHz"] if branch == -1
                 else spectrum["raman_center_plus_GHz"]) + delta_mhz * 1e-3
    x = spectrum["probe_axis_GHz"]
    return {
        "probe_GHz": probe_GHz,
        "G_s": float(np.interp(probe_GHz, x, spectrum["G_s"])),
        "G_c": float(np.interp(probe_GHz, x, spectrum["G_c"])),
        "S_dB": float(np.interp(probe_GHz, x, spectrum["S_dB"])),
    }


# =========================================================
# Scheme wrapper for the generic front-end
# =========================================================
WINDOW_GHZ = 0.55          # half-width of the focused probe window around (−) Raman
TPD_LIMIT_MHZ = 500.0
RESOLUTION = {
    "Fast  (~3 s)":     dict(coarse_points=121, velocity_step=5.0),
    "Balanced  (~6 s)": dict(coarse_points=181, velocity_step=4.0),
    "Fine  (~20 s)":    dict(coarse_points=301, velocity_step=2.0),
}


class FWMScheme(Scheme):
    name = "fwm"
    cluster = "D — Wave mixing"
    title = "85Rb D1 double-Λ four-wave mixing"
    cache_version = "cg-weighted-single-branch-fwm-v1"
    cache_observables = True
    caption = ("Seed/probe gain and intensity-difference squeezing vs two-photon "
               "detuning. OPD and cell parameters recompute (cached); TPD navigates "
               "instantly.")

    def param_schema(self):
        return [
            ParamSpec("opd", "OPD — one-photon detuning Δ", "Detunings", 0.9,
                      -1.0, 3.0, 0.1, "GHz",
                      help="ω_pump = ω(F=2→F'=3) + Δ. Sets where the pump sits; recomputes."),
            ParamSpec("tpd", "TPD — two-photon detuning δ", "Detunings", 0.0,
                      -TPD_LIMIT_MHZ, TPD_LIMIT_MHZ, 1.0, "MHz", recompute=False,
                      help="ω_seed = ω_pump − ν_HF + δ. Navigates the curve instantly (no recompute)."),
            ParamSpec("temp_c", "Temperature", "Cell & beams", 121.0,
                      60.0, 150.0, 1.0, "°C"),
            ParamSpec("pump_mw", "Pump power", "Cell & beams", 600.0,
                      50.0, 1200.0, 10.0, "mW"),
            ParamSpec("probe_uw", "Seed / probe power", "Cell & beams", 10.0,
                      1.0, 200.0, 1.0, "µW"),
            ParamSpec("loss_pct", "Loss after cell", "Detection & scaling", 0.0,
                      0.0, 50.0, 0.5, "%", help="Folds into η = QE × (1 − loss)."),
            ParamSpec("ls", "FWM coupling scale", "Detection & scaling", 0.05,
                      0.01, 1.0, 0.01, "",
                      help="Residual propagation-coupling scale after applying "
                           "Rb85 D1 hyperfine Clebsch-Gordan strengths."),
            ParamSpec("resolution", "Resolution", "Numerics", "Balanced  (~6 s)",
                      choices=tuple(RESOLUTION.keys()), advanced=True),
        ]

    def presets(self):
        return [Preset(
            "Sim et al. 85Rb optimum",
            values=dict(opd=0.9, tpd=-8.0, temp_c=121.0, pump_mw=600.0,
                        probe_uw=8.0, loss_pct=5.5),
            help="One click → squeezing-optimised 85Rb conditions from Sim, Kim & "
                 "Moon, Sci. Rep. 15, 7727 (2025): Δ = 0.9 GHz, δ = −8 MHz, "
                 "T = 121 °C, pump 600 mW, seed 8 µW, loss 5.5 %.",
        )]

    def info(self):
        return (
            "**Source:** G. Sim, H. Kim, H. S. Moon, *Sci. Rep.* **15**, 7727 (2025).\n\n"
            "| Parameter | Paper value | Slider |\n|---|---|---|\n"
            "| One-photon detuning Δ | ≈ 0.9 GHz | ✅ |\n"
            "| Two-photon detuning δ | −8 MHz | ✅ |\n"
            "| Cell temperature | 121 °C | ✅ |\n"
            "| Pump power | 600 mW | ✅ |\n"
            "| Probe-seed power | 8 µW (squeezing run) | ✅ |\n"
            "| Optical loss after cell | 5.5 % | ✅ (loss) |\n"
            "| Cell length | 12.5 mm | fixed |\n"
            "| Pump / seed waist w₀ | 530 / 330 µm (1/e² radius) | fixed |\n"
            "| Measured result | gain ≈ 15, IDS −7.8 dB | — |\n\n"
            "ℹ️ **Line-strength factor** is a model calibration knob, not a paper "
            "value. Tune it until on-resonance gain matches the paper's ≈ 15."
        )

    def compute(self, params):
        center = branch_center_GHz(params["opd"], -1)
        res = RESOLUTION[params["resolution"]]
        return compute_spectrum(
            params["opd"],
            T=params["temp_c"] + 273.15,
            P_pump=params["pump_mw"] * 1e-3,
            P_probe=params["probe_uw"] * 1e-6,
            line_strength=params["ls"],
            loss_frac=params["loss_pct"] / 100.0,
            coarse_points=res["coarse_points"], fine_points=0,
            scan_min=center - WINDOW_GHZ, scan_max=center + WINDOW_GHZ,
            velocity_step=res["velocity_step"], velocity_cutoff=3.0,
            branch=DEFAULT_BRANCH,
        )

    def observables(self, raw, params):
        import matplotlib.pyplot as plt

        tpd = params["tpd"]
        op = operating_point(raw, tpd, branch=-1)
        d_axis = (raw["probe_axis_GHz"] - raw["raman_center_minus_GHz"]) * 1e3

        fig, (axG, axS) = plt.subplots(2, 1, figsize=(8.5, 6.4), sharex=True)
        for ax in (axG, axS):
            ax.grid(alpha=0.3)
        axG.plot(d_axis, raw["G_s"], color="#1f77b4", lw=1.8)
        axG.axvline(tpd, color="crimson", ls="--", lw=1.2)
        axG.axhline(1.0, color="black", lw=0.6)
        axG.scatter([tpd], [op["G_s"]], color="crimson", zorder=5)
        axG.set_ylabel("Seed / probe gain  $G_s$")
        axG.set_title(f"Δ = {params['opd']:.1f} GHz,  T = {params['temp_c']:.0f} °C,  "
                      f"η = {raw['eta']:.3f}")
        if np.nanmax(raw["G_s"]) > 50:
            axG.set_yscale("log")
        axS.plot(d_axis, raw["S_dB"], color="#2ca02c", lw=1.8)
        axS.axvline(tpd, color="crimson", ls="--", lw=1.2)
        axS.axhline(0.0, color="black", lw=0.6)
        axS.scatter([tpd], [op["S_dB"]], color="crimson", zorder=5)
        axS.set_ylabel("Intensity-difference\nsqueezing  [dB]")
        axS.set_xlabel("Two-photon detuning δ  [MHz]   (probe on the − Raman branch)")
        axS.set_xlim(-TPD_LIMIT_MHZ, TPD_LIMIT_MHZ)
        fig.tight_layout()

        metrics = [
            dict(label="Seed / probe gain  G_s", value=f"{op['G_s']:.2f}",
                 help="Power gain of the seeded probe through the cell."),
            dict(label="Squeezing", value=f"{op['S_dB']:.2f} dB",
                 delta="below shot noise" if op["S_dB"] < 0 else "above shot noise",
                 delta_color="inverse"),
            dict(label="Conjugate gain  G_c", value=f"{op['G_c']:.2f}",
                 help="Generated conjugate power gain (drives the twin-beam squeezing)."),
        ]
        derived = (
            f"| Quantity | Value |\n|---|---|\n"
            f"| N(85Rb) | {raw['N_atoms']:.3e} /m³ |\n"
            f"| σ_v (1-D thermal) | {raw['sigma_v']:.1f} m/s |\n"
            f"| Velocity classes | {raw['n_velocity']} |\n"
            f"| Ω_pump / 2π | {raw['Op_A_2pi_GHz']:.3f} GHz |\n"
            f"| Ω_seed / 2π | {raw['Os_2pi_MHz']:.3f} MHz |\n"
            f"| (−) Raman line (probe axis) | {raw['raman_center_minus_GHz']:.3f} GHz |\n"
            f"| Detection η = QE·(1−loss) | {raw['eta']:.4f} |\n"
            f"| Operating probe detuning | {op['probe_GHz']:.4f} GHz |\n\n"
            f"Fixed: cell L = {L_CELL*1e3:.1f} mm · pump w₀ {W_PUMP*1e6:.0f} µm · "
            f"seed w₀ {W_PROBE*1e6:.0f} µm · QE {QE_DETECTOR*100:.2f}% · "
            f"responsivity {RESPONSIVITY_AW} A/W @ 795 nm · pump⊥probe at PBS."
        )
        # ---- Twin-beam coincidence / correlations (ideal parametric) ----
        coinc = observables.coincidence_stats(raw["G_s"], raw["G_c"])
        ci = int(np.argmin(np.abs(d_axis - tpd)))      # nearest operating-point index
        g2sc_op = coinc["g2_sc"][ci]
        R_op = coinc["cauchy_schwarz"][ci]
        ns_op = coinc["n_s"][ci]
        nc_op = coinc["n_c"][ci]
        has_gain = bool(np.isfinite(g2sc_op))

        figC, (axN, axG2) = plt.subplots(2, 1, figsize=(8.5, 6.4), sharex=True)
        for ax in (axN, axG2):
            ax.grid(alpha=0.3)
        axN.plot(d_axis, coinc["n_s"], color="#1f77b4", lw=1.6, label="signal $n_s = G_s-1$")
        axN.plot(d_axis, coinc["n_c"], color="#ff7f0e", lw=1.2, ls="--", label="conjugate $n_c = G_c$")
        axN.axvline(tpd, color="crimson", ls="--", lw=1.2)
        axN.set_ylabel("Generated photons / mode")
        axN.set_title("Twin-beam photon pairs (ideal parametric, gain region)")
        axN.legend(fontsize=9)
        axG2.plot(d_axis, coinc["g2_sc"], color="#2ca02c", lw=1.6)
        axG2.axhline(2.0, color="black", lw=0.7, ls=":")
        axG2.text(d_axis.min(), 2.05, "classical bound g²=2", fontsize=8, va="bottom")
        axG2.axvline(tpd, color="crimson", ls="--", lw=1.2)
        axG2.set_ylabel("Cross-correlation $g^{(2)}_{sc}(0)$")
        axG2.set_xlabel("Two-photon detuning δ  [MHz]   (probe on the − Raman branch)")
        axG2.set_xlim(-TPD_LIMIT_MHZ, TPD_LIMIT_MHZ)
        axG2.set_ylim(1.8, 12)
        figC.tight_layout()

        if has_gain:
            coinc_rows = (
                f"| Generated photons/mode n_s = G_s−1 | {ns_op:.3f} |\n"
                f"| Conjugate photons n_c = G_c | {nc_op:.3f} |\n"
                f"| Cross-correlation g²_sc(0) = 2 + 1/n_s | {g2sc_op:.3f} |\n"
                f"| Auto-correlation g²_ss = g²_cc | 2.000 (thermal each arm) |\n"
                f"| Cauchy-Schwarz R = g²_sc²/(g²_ss g²_cc) | {R_op:.3f} |\n"
                f"| Nonclassical (R > 1) | {'yes' if R_op > 1 else 'no'} |\n"
            )
        else:
            coinc_rows = "| (no net gain at this δ — no photon pairs) | — |\n"
        coinc_table = (
            f"Ideal (lossless) twin-beam photon-pair statistics at the operating "
            f"point (δ = {tpd:.0f} MHz), from the parametric gain — the spontaneous "
            f"/ coincidence-counting regime:\n\n"
            f"| Quantity | Value |\n|---|---|\n"
            + coinc_rows +
            "\ng²_sc > 2 (and R > 1) violates the classical Cauchy-Schwarz "
            "inequality — the hallmark of FWM photon-pair correlations. g²_sc → 2 "
            "at high gain, → large near threshold. Propagation loss is not modelled "
            "here (as in the squeezing panel)."
        )

        return {
            "metrics": metrics,
            "figure": fig,
            "figures": [("Twin-beam coincidence / correlations", figC)],
            "tables": [
                {"title": "Derived quantities", "markdown": derived},
                {"title": "Twin-beam coincidence (spontaneous)", "markdown": coinc_table},
            ],
        }

    def extra_views(self):
        def _compute_full(params):
            return full_spectrum(
                params["opd"], params["temp_c"] + 273.15,
                params["pump_mw"], params["probe_uw"], params["ls"], params["loss_pct"])

        def _render_full(full):
            import matplotlib.pyplot as plt
            figF, (aG, aS) = plt.subplots(2, 1, figsize=(8.5, 6.4), sharex=True)
            for ax in (aG, aS):
                ax.grid(alpha=0.3)
            styles = {
                "minus": dict(color="#1f77b4", label="minus Raman branch"),
                "plus": dict(color="#ff7f0e", label="plus Raman branch"),
            }
            for key, style in styles.items():
                spec = full[key]
                aG.plot(spec["probe_axis_GHz"], spec["G_s"], lw=1.4, **style)
                aS.plot(spec["probe_axis_GHz"], spec["S_dB"], lw=1.4, **style)
            aG.axhline(1.0, color="black", lw=0.6)
            aG.set_ylabel("Seed / probe gain  $G_s$")
            if max(np.nanmax(full[key]["G_s"]) for key in styles) > 50:
                aG.set_yscale("log")
            aS.axhline(0.0, color="black", lw=0.6)
            aS.set_ylabel("Squeezing [dB]")
            aS.set_xlabel(r"Probe detuning from $F=2\to F'=3$  [GHz]")
            for a in (aG, aS):
                a.axvline(full["D_GHz"], color="gray", ls=":", lw=0.8)
            aG.legend(fontsize=9)
            aS.legend(fontsize=9)
            figF.tight_layout()
            return figF

        return [ExtraView(
            key="Full −8…12 GHz probe scan (slow, both Raman branches)",
            description="The focused view zooms on the (−) Raman line. This runs the "
                        "original wide scan showing both branches.",
            compute=_compute_full, render=_render_full,
        )]


def full_spectrum(D_GHz, T_K, P_pump_mW, P_probe_uW, line_strength, loss_pct):
    """Wide scan with the two Raman channels calculated independently."""
    common = dict(
        T=T_K, P_pump=P_pump_mW * 1e-3, P_probe=P_probe_uW * 1e-6,
        line_strength=line_strength, loss_frac=loss_pct / 100.0,
        coarse_points=301, fine_points=401, velocity_step=2.0)
    return {
        "D_GHz": D_GHz,
        "minus": compute_spectrum(D_GHz, branch=-1, **common),
        "plus": compute_spectrum(D_GHz, branch=+1, **common),
    }
