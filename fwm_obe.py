"""
85Rb D1 double-Λ FWM minimal 4-level density-matrix model.

Atomic level labels (this convention)
-------------------------------------
  g₁ = |F=2⟩ (lower ground)        e₂ = |F'=2⟩ (lower excited)
  g₂ = |F=3⟩ (upper ground)        e₃ = |F'=3⟩ (upper excited)

  Δ (one-photon)  : ω_pump = ω_{F=2→F'=3} + Δ
  δ (two-photon)  : ω_seed = ω_pump − ν_HF + δ          (standard Stokes / branch −)

Field layout
------------
  Pump A : |g₁⟩↔|e₂,e₃⟩            (ω_p, one-photon detuning Δ from g₁↔e₃)
  Pump B : |g₂⟩↔|e₂,e₃⟩            (same ω_p; one-photon detuning Δ + ν_HF from g₂↔e₃)
  Probe (−) : |g₂⟩↔|e₂,e₃⟩, ω_s = ω_p − ν_HF + δ₋       (standard FWM seed)
  Probe (+) : |g₁⟩↔|e₂,e₃⟩, ω_s = ω_p + ν_HF + δ₊       (other Raman branch)
  Conjugate : energy conservation, 2 ω_p = ω_s + ω_c

Rotating frame + RWA + 3-mode sideband expansion (modes −1, 0, +1).
Each velocity gives a 16×16 Bloch system (3-mode block-eliminated to ρ_0).

Optimisation
------------
The only v-dependence of L0 is the excited-state diagonal shift
        H₀(v) = H₀(v=0) − Δ_eff · diag(0,0,1,1),   Δ_eff = Δ − k·v
so      L₀(Δ_eff) = L₀(Δ_eff=0) − Δ_eff · S_v,    S_v = comm_super(diag(0,0,1,1)).
That means every velocity in a Doppler average shares the same L₀ template and
only differs by a single diagonal super-operator scaled by Δ_eff. The Bloch
solve is therefore vectorised over Δ_eff with one batched np.linalg.solve call
instead of one Python loop per velocity.

The single-velocity response ρ(δ, Δ_eff) is built once as a 2-D table:
    R(δ, Δ_eff)  is independent of temperature   (T only sets Maxwell weights)
                 and independent of Δ            (Δ is an offset on the Δ_eff axis)
so multiple (T, Δ) configurations reuse the same table.

Observable pipeline
-------------------
Two solves at unit (Ω_s, Ω_c) extract the linear-response χ̄ matrix
        ρ_probe  =  χ̄_ss · Ω_s + χ̄_sc · Ω_c*
        ρ_conj   =  χ̄_cs · Ω_s + χ̄_cc · Ω_c*
which is converted to a physical susceptibility
        χ_xy(ω) = −2 N |d|² / (ε₀ ℏ) · χ̄_xy
and propagated through L_cell with a 2×2 Maxwell-Bloch matrix exponential to
give amplitude gains G_s = |T₀₀|², G_c = |T₁₀|² for an input seed in the probe
channel and vacuum in the conjugate. Intensity-difference squeezing follows
from the standard FWM-twin-beam formula with detection efficiency η.
"""

from pathlib import Path
import os

os.environ.setdefault("OMP_NUM_THREADS", str(os.cpu_count() or 1))
os.environ.setdefault("OPENBLAS_NUM_THREADS", str(os.cpu_count() or 1))
os.environ.setdefault("MKL_NUM_THREADS", str(os.cpu_count() or 1))

import numpy as np
import matplotlib.pyplot as plt

SCRIPT_DIR = Path(__file__).resolve().parent


# =========================================================
# 1. Constants
# =========================================================
HBAR = 1.054571817e-34
KB = 1.380649e-23
C_LIGHT = 299792458.0
EPS_0 = 8.8541878128e-12
ELEMENTARY_CHARGE = 1.602176634e-19

# References
# [1] Daniel A. Steck, "Rubidium 85 D Line Data," http://steck.us/alkalidata
#     (revision 2.3.4, 8 August 2025).
# [2] G. Sim, H. Kim, and H. S. Moon, Sci. Rep. 15, 7727 (2025),
#     doi:10.1038/s41598-025-86479-w.
NU_D1_85RB = 377.107_385_690e12
WAVELENGTH_D1_85RB = C_LIGHT / NU_D1_85RB
GAMMA_2PI = 5.746e6
NU_GROUND_HF = 3.035732439e9
NU_EXCITED_HF_D1 = 361.58e6
NU_HF = NU_GROUND_HF
MASS_85RB = 1.4100e-25
I_SAT = 4.484e-3 * 1e4
DIPOLE_D1 = 2.5377e-29          # ⟨J=1/2‖er‖J=1/2⟩ reduced matrix element [C·m]
RB85_ABUNDANCE = 0.7217

# Effective line-strength factor.
#
# The reduced J-basis matrix element 2.5377e-29 C·m above is what gives the
# full natural width Γ in a one-channel 2-level radiator. Real D1 atoms split
# this oscillator strength over multiple F,m_F → F',m_F' Zeeman channels
# according to Clebsch-Gordan factors, so a single Λ-pair in a 4-level
# reduction couples weaker than the J-basis element suggests. The 4-level
# model lumps all m-substates together and otherwise has no way to express
# that selection; LINE_STRENGTH_FACTOR scales the *effective* |d|² used in
# the Maxwell-Bloch propagation, leaving the OBE solve (which uses Γ and Ω
# via I_sat, both already empirical) unchanged. Tune to match experiment.
LINE_STRENGTH_FACTOR = 1.0

GAMMA = 2 * np.pi * GAMMA_2PI
OMEGA_HF = 2 * np.pi * NU_HF
OMEGA_EXCITED_HF = 2 * np.pi * NU_EXCITED_HF_D1
K_VEC = 2 * np.pi / WAVELENGTH_D1_85RB
OMEGA_D1 = 2 * np.pi * NU_D1_85RB

# Cell + beam geometry (Sim et al. 2025).
# Pump diameter 530 μm → 1/e² waist 265 μm; probe diameter 330 μm → 165 μm.
# Polarisations: pump and probe orthogonal at the PBS. The 4-level effective
# model averages over Zeeman m-substates so it does not resolve the σ⁺/σ⁻ split
# explicitly; the orthogonal-polarisation geometry is reflected only through
# the independent powers/waists of the two channels.
L_CELL  = 12.5e-3
W_PUMP  = 530e-6 / 2
W_PROBE = 330e-6 / 2
P_PUMP, P_PROBE = 600e-3, 10e-6
T_CELL = 394.15

# Detection.
QE_DETECTOR     = 0.9047
RESPONSIVITY_AW = 0.58           # consistent with QE × eλ/(hc) at 795 nm
LOSS_FRAC       = 0.0            # additional fractional intensity loss after the cell (0..1)
ETA_TOTAL       = QE_DETECTOR * (1.0 - LOSS_FRAC)

# Phenomenological ground-coherence decay.
GAMMA_GG_2PI = 100e3
GAMMA_GG     = 2 * np.pi * GAMMA_GG_2PI

SIGMA_V = np.sqrt(KB * T_CELL / MASS_85RB)

# FWM scan parameters.
OMEGA_C_SEED         = 0.0
SCAN_MIN_GHZ         = -8.0
SCAN_MAX_GHZ         = 12.0
SCAN_COARSE_POINTS   = 401
RESONANCE_WINDOW_MHZ = 80.0
SCAN_FINE_POINTS     = 801
PUMP_OVERLAP_EXCLUSION_MHZ = 1e-3
VELOCITY_STEP_MPS    = 1.0
VELOCITY_CUTOFF_SIGMA = 3.0
DELTA_GHZ_LIST       = [0.9]
BRANCHES             = (-1, +1)

G1, G2, E2, E3 = 0, 1, 2, 3
GROUND_STATES  = (G1, G2)
EXCITED_STATES = (E2, E3)
N_LEVELS = 4
RHO_DIM  = N_LEVELS * N_LEVELS
TRACE_RHO_ROW = G1 * N_LEVELS + G1

TRACE_RHS = np.zeros(RHO_DIM, dtype=complex)
TRACE_RHS[TRACE_RHO_ROW] = 1
EYE_RHO = np.eye(RHO_DIM, dtype=complex)


# =========================================================
# 2. Basic helpers
# =========================================================
def rabi_freq(power, waist):
    """I = 2P/(π w₀²),  Ω = Γ √(I / 2 I_sat).  rad/s."""
    I = 2 * power / (np.pi * waist**2)
    return GAMMA * np.sqrt(I / (2 * I_SAT))


def rb85_density(T):
    """
    Rb-85 atomic number density [/m³] from Steck Rb vapor pressure
    (liquid phase, T > 312.46 K). 85-Rb isotopic abundance applied.
    """
    log10_P_torr = (15.88253 - 4529.635 / T
                    + 0.00058663 * T - 2.99138 * np.log10(T))
    P_pa = 10 ** log10_P_torr * 133.322387415
    N_total = P_pa / (KB * T)
    return RB85_ABUNDANCE * N_total


def rho_index(row, col):
    return row * N_LEVELS + col


# =========================================================
# 3. Hamiltonians
# =========================================================
def _add_static_drive(H, ground, omega):
    for excited in EXCITED_STATES:
        H[ground, excited] += omega / 2
        H[excited, ground] += omega / 2


def _add_sideband_drive(H, ground, omega):
    for excited in EXCITED_STATES:
        H[excited, ground] += omega / 2


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


# =========================================================
# 4. Super-operators
# =========================================================
def comm_super(H):
    """vec rule (row-major, vec[Ni+j] = ρ_{ij}):  −i [H, ·]."""
    eye = np.eye(N_LEVELS, dtype=complex)
    return -1j * (np.kron(H, eye) - np.kron(eye, H.T))


def _build_lindblad_fixed():
    eye = np.eye(N_LEVELS, dtype=complex)

    def dissipator(Lop):
        LdL = Lop.conj().T @ Lop
        return (np.kron(Lop, Lop.conj())
                - 0.5 * np.kron(LdL, eye)
                - 0.5 * np.kron(eye, LdL.T))

    kets = np.eye(N_LEVELS, dtype=complex)
    rate = GAMMA / 2
    D = np.zeros((RHO_DIM, RHO_DIM), dtype=complex)
    for excited in EXCITED_STATES:
        for ground in GROUND_STATES:
            L = np.sqrt(rate) * np.outer(kets[ground], kets[excited])
            D += dissipator(L)
    D[rho_index(G1, G2), rho_index(G1, G2)] -= GAMMA_GG
    D[rho_index(G2, G1), rho_index(G2, G1)] -= GAMMA_GG
    return D


L_LINDBLAD = _build_lindblad_fixed()


def _velocity_shift_super():
    """S_v = comm_super(diag(0,0,1,1)) so that L₀(Δ_eff) = L₀(0) − Δ_eff · S_v."""
    Dee = np.zeros((N_LEVELS, N_LEVELS), dtype=complex)
    Dee[E2, E2] = 1.0
    Dee[E3, E3] = 1.0
    return comm_super(Dee)


S_V_SUPER = _velocity_shift_super()


def L0_of(H0):
    return comm_super(H0) + L_LINDBLAD


def sideband_template(Op_A, Op_B, Oc, branch):
    Hp = sideband_hamiltonian(Op_A, Op_B, Oc, branch)
    Cp = comm_super(Hp)
    Cm = comm_super(Hp.conj().T)
    return Cp, Cm


# =========================================================
# 5. Batched sideband solver over Δ_eff
# =========================================================
def solve_sidebands_batched(L0_at_Deff_zero, Cp, Cm, Omega_beat, Delta_eff_axis):
    """
    Solve the 3-mode sideband steady state for every Δ_eff in the array.

    Returns ρ_0 and ρ_+1 reshaped to (N_deff, 4, 4) each.
    """
    n = Delta_eff_axis.size
    M = RHO_DIM
    L0_batch = (L0_at_Deff_zero[None, :, :]
                - Delta_eff_axis[:, None, None] * S_V_SUPER[None, :, :])

    iO = 1j * Omega_beat * EYE_RHO
    A_minus = L0_batch - iO[None, :, :]
    A_plus  = L0_batch + iO[None, :, :]

    Cm_batch = np.broadcast_to(Cm, (n, M, M))
    Cp_batch = np.broadcast_to(Cp, (n, M, M))

    Am_inv_Cm = np.linalg.solve(A_minus, Cm_batch)
    Ap_inv_Cp = np.linalg.solve(A_plus,  Cp_batch)

    minus_feedback = Cp_batch @ Am_inv_Cm
    plus_feedback  = Cm_batch @ Ap_inv_Cp

    A_eff = L0_batch - minus_feedback - plus_feedback
    A_eff[:, TRACE_RHO_ROW, :] = 0
    for state in range(N_LEVELS):
        A_eff[:, TRACE_RHO_ROW, rho_index(state, state)] = 1

    rhs = np.broadcast_to(TRACE_RHS[:, None], (n, M, 1)).copy()
    rho_0_vec  = np.linalg.solve(A_eff, rhs)
    cp_rho0 = np.einsum("ij,njk->nik", Cp, rho_0_vec)
    rho_p1_vec = -np.linalg.solve(A_plus, cp_rho0)

    rho_0  = rho_0_vec[:, :, 0].reshape(n, N_LEVELS, N_LEVELS)
    rho_p1 = rho_p1_vec[:, :, 0].reshape(n, N_LEVELS, N_LEVELS)
    return rho_0, rho_p1


# =========================================================
# 6. χ-matrix table   (T- and Δ-independent)
# =========================================================
def chi_matrix_table(Op_A, Op_B, Os_ref, Oc_ref, delta_axis, Delta_eff_axis, branch):
    """
    Two solves per probe-detuning point to extract the full linear-response
    matrix (χ̄_ss, χ̄_cs, χ̄_sc, χ̄_cc) on a 2-D (δ, Δ_eff) grid.

    Solve 1: (Ω_s = Os_ref, Ω_c = 0)        → χ̄_ss, χ̄_cs
    Solve 2: (Ω_s = 0,      Ω_c = Oc_ref)   → χ̄_sc, χ̄_cc

    Returns each array with shape (n_delta, n_deff), complex.
    """
    probe_ground = G2 if branch == -1 else G1
    conj_ground  = G1 if branch == -1 else G2
    n_d  = delta_axis.size
    n_de = Delta_eff_axis.size

    chi_ss = np.zeros((n_d, n_de), dtype=complex)
    chi_cs = np.zeros((n_d, n_de), dtype=complex)
    chi_sc = np.zeros((n_d, n_de), dtype=complex)
    chi_cc = np.zeros((n_d, n_de), dtype=complex)

    Cp_no_c, Cm_no_c = sideband_template(Op_A, Op_B, 0.0,     branch)   # solve 1
    Cp_c,    Cm_c    = sideband_template(Op_A, Op_B, Oc_ref,  branch)   # solve 2

    for i, delta in enumerate(delta_axis):
        Omega_beat = OMEGA_HF - branch * delta

        # ---- Solve 1: probe drive only ----
        H0_1 = static_hamiltonian_at_Deff_zero(Op_A, Op_B, Os_ref, delta, branch)
        L0_1 = L0_of(H0_1)
        rho0_a, rhop_a = solve_sidebands_batched(
            L0_1, Cp_no_c, Cm_no_c, Omega_beat, Delta_eff_axis)
        probe_a = sum(rho0_a[:, e, probe_ground] for e in EXCITED_STATES)
        conj_a  = sum(rhop_a[:, e, conj_ground]  for e in EXCITED_STATES)
        chi_ss[i] = probe_a / Os_ref
        chi_cs[i] = conj_a  / Os_ref

        # ---- Solve 2: conjugate seed only ----
        H0_2 = static_hamiltonian_at_Deff_zero(Op_A, Op_B, 0.0, delta, branch)
        L0_2 = L0_of(H0_2)
        rho0_b, rhop_b = solve_sidebands_batched(
            L0_2, Cp_c, Cm_c, Omega_beat, Delta_eff_axis)
        probe_b = sum(rho0_b[:, e, probe_ground] for e in EXCITED_STATES)
        conj_b  = sum(rhop_b[:, e, conj_ground]  for e in EXCITED_STATES)
        chi_sc[i] = probe_b / Oc_ref
        chi_cc[i] = conj_b  / Oc_ref

    return chi_ss, chi_cs, chi_sc, chi_cc


# =========================================================
# 7. Doppler average  (T-dependence enters here, not in R)
# =========================================================
def velocity_grid(T=T_CELL, dv=None, cutoff_sigma=None):
    sigma = np.sqrt(KB * T / MASS_85RB)
    dv = VELOCITY_STEP_MPS if dv is None else dv
    cutoff_sigma = VELOCITY_CUTOFF_SIGMA if cutoff_sigma is None else cutoff_sigma
    v_limit = np.ceil(cutoff_sigma * sigma / dv) * dv
    v   = np.arange(-v_limit, v_limit + 0.5 * dv, dv)
    pdf = np.exp(-v**2 / (2 * sigma**2)) / (np.sqrt(2 * np.pi) * sigma)
    w   = pdf * dv
    return v, w / w.sum()


def build_Delta_eff_axis(Delta_min, Delta_max, v_grid):
    dv = v_grid[1] - v_grid[0]
    step = K_VEC * dv
    lo = Delta_min - K_VEC * v_grid.max()
    hi = Delta_max - K_VEC * v_grid.min()
    n = int(np.ceil((hi - lo) / step)) + 1
    return np.linspace(lo, hi, n)


def doppler_average(chi_table, Delta_eff_axis, Delta, v_grid, weights):
    """
    Sum   Σ_v weights(v) · χ_table[δ, Δ_eff = Δ − k·v]   for every δ row.
    Linear interpolation along the Δ_eff axis.
    """
    deff_v = Delta - K_VEC * v_grid
    idx_float = np.interp(deff_v, Delta_eff_axis, np.arange(Delta_eff_axis.size))
    n_de = Delta_eff_axis.size
    idx_lo = np.clip(np.floor(idx_float).astype(int), 0, n_de - 2)
    frac = (idx_float - idx_lo).astype(chi_table.dtype)
    lo_part = chi_table[:, idx_lo]
    hi_part = chi_table[:, idx_lo + 1]
    interp = lo_part * (1 - frac)[None, :] + hi_part * frac[None, :]
    return interp @ weights


# =========================================================
# 8. Propagation  →  Gain  →  Squeezing
# =========================================================
def matrix_exp_2x2(M, L):
    """Closed-form exp(M·L) for a batched stack of complex 2×2 matrices."""
    s   = 0.5 * (M[..., 0, 0] + M[..., 1, 1])
    q00 =  M[..., 0, 0] - s
    q01 =  M[..., 0, 1]
    q10 =  M[..., 1, 0]
    q11 =  M[..., 1, 1] - s
    c2  = q00 * q11 - q01 * q10           # = det(q),  with trace(q)=0
    c   = np.sqrt(-c2 + 0j)               # exp(q L) = cosh(cL)I + sinh(cL)/c · q
    big = np.abs(c) > 1e-30
    safe_c = np.where(big, c, 1.0)
    sinh_over_c = np.where(big, np.sinh(c * L) / safe_c, L * np.ones_like(c))
    cosh_cL = np.cosh(c * L)
    exp_sL  = np.exp(s * L)

    out = np.empty_like(M)
    out[..., 0, 0] = exp_sL * (cosh_cL + sinh_over_c * q00)
    out[..., 0, 1] = exp_sL * (         sinh_over_c * q01)
    out[..., 1, 0] = exp_sL * (         sinh_over_c * q10)
    out[..., 1, 1] = exp_sL * (cosh_cL + sinh_over_c * q11)
    return out


def gain_from_chi(chi_ss_avg, chi_sc_avg, chi_cs_avg, chi_cc_avg,
                  k_probe, k_conj, L, N_atoms,
                  dipole=None, line_strength=None):
    """
    Linearised Maxwell-Bloch propagation:

        d/dz [Ω_s, Ω_c*]ᵀ = M · [Ω_s, Ω_c*]ᵀ
        M[0,0] =  i k_s / 2 · χ_phys_ss
        M[0,1] =  i k_s / 2 · χ_phys_sc
        M[1,0] = −i k_c / 2 · χ_phys_cs*
        M[1,1] = −i k_c / 2 · χ_phys_cc*
        χ_phys_xy = −2 N |d_eff|² / (ε₀ ℏ) · χ̄_xy
        |d_eff|²  = LINE_STRENGTH_FACTOR · |d|²    (Clebsch-Gordan rescaling)

    Returns probe amplitude gain G_s = |T₀₀|² and conjugate intensity
    G_c = |T₁₀|² for input (Ω_s, Ω_c* = 0).
    """
    if dipole is None:
        dipole = DIPOLE_D1
    if line_strength is None:
        line_strength = LINE_STRENGTH_FACTOR
    coupling = -2.0 * N_atoms * line_strength * dipole**2 / (EPS_0 * HBAR)
    n = chi_ss_avg.size
    M = np.zeros((n, 2, 2), dtype=complex)
    M[:, 0, 0] =  0.5j * k_probe * coupling * chi_ss_avg
    M[:, 0, 1] =  0.5j * k_probe * coupling * chi_sc_avg
    M[:, 1, 0] = -0.5j * k_conj  * coupling * chi_cs_avg.conj()
    M[:, 1, 1] = -0.5j * k_conj  * coupling * chi_cc_avg.conj()

    T = matrix_exp_2x2(M, L)
    G_s = np.abs(T[:, 0, 0]) ** 2
    G_c = np.abs(T[:, 1, 0]) ** 2
    return G_s, G_c, T


def intensity_difference_squeezing_dB(G_s, G_c, eta):
    """
    Ideal twin-beam intensity-difference noise (lossless):
        S_ideal = (G_s + G_c − 2√(G_s·G_c − 1)) / (G_s + G_c)
    Symmetric detection efficiency η on both arms:
        S(η) = η · S_ideal + (1 − η)
    Returns 10·log10(S) in dB (negative ↔ squeezed).
    """
    G_s = np.asarray(G_s, dtype=float)
    G_c = np.asarray(G_c, dtype=float)
    cross = np.sqrt(np.maximum(G_s * G_c - 1.0, 0.0))
    S_ideal = (G_s + G_c - 2.0 * cross) / np.maximum(G_s + G_c, 1e-30)
    S_ideal = np.clip(S_ideal, 0.0, None)
    S = eta * S_ideal + (1.0 - eta)
    return 10.0 * np.log10(np.maximum(S, 1e-30))


# =========================================================
# 9. Probe-detuning axis
# =========================================================
def branch_center_GHz(Delta_GHz, branch):
    if branch not in BRANCHES:
        raise ValueError(f"branch must be one of {BRANCHES}, got {branch}")
    return Delta_GHz + branch * NU_HF / 1e9


def probe_scan_axis_GHz(Delta_GHz, coarse_points=None, fine_points=None,
                        window_mhz=None, scan_min=None, scan_max=None):
    """
    Probe-frequency scan axis [GHz], referenced to the 85Rb F=2 → F'=3 line.
    The pump (at ω_pump = ω_{F=2→F'=3} + Δ) lands at probe_axis = Delta_GHz;
    the standard Stokes / branch (−) Raman peak sits at Delta_GHz − ν_HF/1e9.
    """
    coarse_points = SCAN_COARSE_POINTS if coarse_points is None else coarse_points
    fine_points   = SCAN_FINE_POINTS   if fine_points   is None else fine_points
    window_mhz    = RESONANCE_WINDOW_MHZ if window_mhz  is None else window_mhz
    scan_min = SCAN_MIN_GHZ if scan_min is None else scan_min
    scan_max = SCAN_MAX_GHZ if scan_max is None else scan_max

    coarse = np.linspace(scan_min, scan_max, coarse_points)
    half_window = window_mhz * 1e-3
    parts = [coarse]
    for branch in BRANCHES:
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


# =========================================================
# 9b. High-level spectrum  (one call → gain + squeezing curves)
# =========================================================
def compute_spectrum(D_GHz, *,
                     T=T_CELL, P_pump=P_PUMP, P_probe=P_PROBE,
                     w_pump=W_PUMP, w_probe=W_PROBE,
                     line_strength=None, loss_frac=LOSS_FRAC, qe=QE_DETECTOR,
                     coarse_points=None, fine_points=None, window_mhz=None,
                     scan_min=None, scan_max=None,
                     velocity_step=None, velocity_cutoff=None,
                     branches=BRANCHES):
    """
    Full pipeline for one one-photon detuning Δ = 2π·D_GHz·1e9.

    Sweeps the probe over the scan axis, Doppler-averages the χ̄ matrix, and
    propagates to seed gain G_s, conjugate gain G_c, and intensity-difference
    squeezing S_dB. Returns a dict; every physical knob is an argument so the
    function is side-effect free and cache-friendly.
    """
    if line_strength is None:
        line_strength = LINE_STRENGTH_FACTOR
    eta = qe * (1.0 - loss_frac)

    Op_A = rabi_freq(P_pump, w_pump)
    Op_B = Op_A
    Os   = rabi_freq(P_probe, w_probe)
    Os_ref = Os
    Oc_ref = Os                              # χ̄ is independent of |Ω_ref|; any nonzero value works

    N_atoms = rb85_density(T)
    Delta = 2 * np.pi * D_GHz * 1e9

    probe_axis_GHz = probe_scan_axis_GHz(
        D_GHz, coarse_points, fine_points, window_mhz, scan_min, scan_max)
    v_grid, weights = velocity_grid(T, velocity_step, velocity_cutoff)
    Delta_eff_axis = build_Delta_eff_axis(Delta, Delta, v_grid)

    chi_ss_avg = np.zeros(probe_axis_GHz.size, dtype=complex)
    chi_cs_avg = np.zeros(probe_axis_GHz.size, dtype=complex)
    chi_sc_avg = np.zeros(probe_axis_GHz.size, dtype=complex)
    chi_cc_avg = np.zeros(probe_axis_GHz.size, dtype=complex)

    for branch in branches:
        delta_axis = two_photon_detuning_from_probe_scan(probe_axis_GHz, D_GHz, branch)
        ch_ss, ch_cs, ch_sc, ch_cc = chi_matrix_table(
            Op_A, Op_B, Os_ref, Oc_ref, delta_axis, Delta_eff_axis, branch)
        chi_ss_avg += doppler_average(ch_ss, Delta_eff_axis, Delta, v_grid, weights)
        chi_cs_avg += doppler_average(ch_cs, Delta_eff_axis, Delta, v_grid, weights)
        chi_sc_avg += doppler_average(ch_sc, Delta_eff_axis, Delta, v_grid, weights)
        chi_cc_avg += doppler_average(ch_cc, Delta_eff_axis, Delta, v_grid, weights)

    G_s, G_c, _ = gain_from_chi(
        chi_ss_avg, chi_sc_avg, chi_cs_avg, chi_cc_avg,
        K_VEC, K_VEC, L_CELL, N_atoms, line_strength=line_strength)
    S_dB = intensity_difference_squeezing_dB(G_s, G_c, eta)

    return {
        "D_GHz": D_GHz,
        "probe_axis_GHz": probe_axis_GHz,
        "G_s": G_s,
        "G_c": G_c,
        "S_dB": S_dB,
        "eta": eta,
        "N_atoms": N_atoms,
        "sigma_v": np.sqrt(KB * T / MASS_85RB),
        "n_velocity": v_grid.size,
        "Op_A_2pi_GHz": Op_A / (2 * np.pi) / 1e9,
        "Os_2pi_MHz": Os / (2 * np.pi) / 1e6,
        "raman_center_minus_GHz": branch_center_GHz(D_GHz, -1),
        "raman_center_plus_GHz": branch_center_GHz(D_GHz, +1),
    }


def operating_point(spectrum, delta_mhz, branch=-1):
    """
    Read G_s and squeezing at a chosen two-photon detuning δ (in MHz) on the
    selected Raman branch, by interpolating an already-computed spectrum.
    Standard FWM seed is branch = −1.
    """
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
# 10. Main
# =========================================================
def main():
    Op_A = rabi_freq(P_PUMP,  W_PUMP)
    Os   = rabi_freq(P_PROBE, W_PROBE)
    N_atoms = rb85_density(T_CELL)

    print(f"T                = {T_CELL:7.2f} K")
    print(f"N(85Rb)          = {N_atoms:.3e} /m^3")
    print(f"Omega_pA / 2pi   = {Op_A /(2*np.pi)/1e9:7.3f} GHz")
    print(f"Omega_s  / 2pi   = {Os  /(2*np.pi)/1e6:7.3f} MHz")
    print(f"sigma_v          = {SIGMA_V:7.1f} m/s")
    print(f"L_cell           = {L_CELL*1e3:7.3f} mm")
    print(f"eta_total        = {ETA_TOTAL:7.4f} (QE {QE_DETECTOR:.4f} × (1-loss {LOSS_FRAC:.3f}))")

    v_grid, _ = velocity_grid(T_CELL)
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
