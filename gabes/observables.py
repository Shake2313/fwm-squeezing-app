"""
χ → physical observables.

Today: linearised Maxwell-Bloch propagation → seed/conjugate gain and
intensity-difference squeezing (FWM). Ported verbatim from fwm_obe.py.
Absorption / OD (Im χ → Beer-Lambert) for the absorption-cluster schemes lands
in Phase 1; the FWM transfer matrix T returned by `gain_from_chi` is also what
the Phase-2 coincidence panel will reuse.
"""
import numpy as np

from . import constants
from .core import matrix_exp_2x2


def gain_from_chi(chi_ss_avg, chi_sc_avg, chi_cs_avg, chi_cc_avg,
                  k_probe, k_conj, L, N_atoms,
                  dipole=None, line_strength=None):
    """
    Linearised Maxwell-Bloch propagation:

        d/dz [Ω_s, Ω_c*]ᵀ = M · [Ω_s, Ω_c*]ᵀ
        χ_phys_xy = −2 N |d_eff|² / (ε₀ ℏ) · χ̄_xy,
        |d_eff|²  = LINE_STRENGTH_FACTOR · |d|²

    Returns probe amplitude gain G_s = |T₀₀|², conjugate G_c = |T₁₀|², and the
    full 2×2 transfer matrix stack T (reused later by the coincidence panel).
    """
    if dipole is None:
        dipole = constants.DIPOLE_D1
    if line_strength is None:
        line_strength = constants.LINE_STRENGTH_FACTOR
    coupling = -2.0 * N_atoms * line_strength * dipole**2 / (constants.EPS_0 * constants.HBAR)
    n = chi_ss_avg.size
    M = np.zeros((n, 2, 2), dtype=complex)
    M[:, 0, 0] = 0.5j * k_probe * coupling * chi_ss_avg
    M[:, 0, 1] = 0.5j * k_probe * coupling * chi_sc_avg
    M[:, 1, 0] = -0.5j * k_conj * coupling * chi_cs_avg.conj()
    M[:, 1, 1] = -0.5j * k_conj * coupling * chi_cc_avg.conj()

    T = matrix_exp_2x2(M, L)
    G_s = np.abs(T[:, 0, 0]) ** 2
    G_c = np.abs(T[:, 1, 0]) ** 2
    return G_s, G_c, T


def coincidence_stats(G_s, G_c):
    """
    Equal-time twin-beam (signal/conjugate) photon statistics for the FWM
    parametric process, in the **ideal (lossless) parametric** limit set by the
    gains G_s, G_c — consistent with how `intensity_difference_squeezing_dB`
    idealises the twin beams (propagation loss is not modelled with quantum
    Langevin noise; folding it in via the bare transfer matrix would corrupt the
    photon statistics, so we use the gain directly).

    A two-mode squeezed vacuum with mean photon number n per mode obeys:
        n_pairs = G_s − 1            (signal photons generated from vacuum)
        g²_ss = g²_cc = 2            (each arm thermal)
        g²_sc(0) = 2 + 1/n_pairs     (cross-correlation, > 2)
        Cauchy-Schwarz  R = [g²_sc]² / (g²_ss g²_cc) = (2 + 1/n)²/4  > 1.
    The cross-correlation g²_sc → 2 at high gain and diverges at low pair flux;
    R > 1 everywhere in the gain region is the nonclassical photon-pair signature.
    Only meaningful where there is net gain (G_s > 1); elsewhere set to NaN.

    Returns a dict of per-point arrays (plus a `gain_mask`).
    """
    G_s = np.asarray(G_s, dtype=float)
    G_c = np.asarray(G_c, dtype=float)
    gain = G_s > 1.0
    n_pairs = np.where(gain, G_s - 1.0, np.nan)
    with np.errstate(divide="ignore", invalid="ignore"):
        g2_sc = 2.0 + 1.0 / n_pairs
    g2_auto = np.where(gain, 2.0, np.nan)
    R = g2_sc ** 2 / (g2_auto * g2_auto)
    return {
        "n_s": n_pairs,
        "n_c": np.where(gain, G_c, np.nan),
        "g2_ss": g2_auto, "g2_cc": g2_auto, "g2_sc": g2_sc,
        "cauchy_schwarz": R, "gain_mask": gain,
    }


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
# Absorption-cluster observables (OD / AT / EIT / CPT)
# =========================================================
def chi_phys(chi_bar, N_atoms, dipole=None, line_strength=None):
    """
    Physical linear susceptibility from the dimensionless χ̄ = ρ_probe / Ω_probe.

    Same coupling convention as `gain_from_chi`:
        χ_phys = −2 N · LINE_STRENGTH_FACTOR · |d|² / (ε₀ ℏ) · χ̄
    so a passive transition is absorptive (Im χ_phys > 0 on resonance) and the
    line-strength factor is the same calibration knob used by the FWM path.
    """
    if dipole is None:
        dipole = constants.DIPOLE_D1
    if line_strength is None:
        line_strength = constants.LINE_STRENGTH_FACTOR
    coupling = -2.0 * N_atoms * line_strength * dipole**2 / (constants.EPS_0 * constants.HBAR)
    return coupling * chi_bar


def absorption_coefficient(chi_bar, k, N_atoms, dipole=None, line_strength=None):
    """α = k · Im(χ_phys)  [1/m].  Returns (α, χ_phys)."""
    xp = chi_phys(chi_bar, N_atoms, dipole=dipole, line_strength=line_strength)
    return k * np.imag(xp), xp


def transmission(alpha, L):
    """Beer-Lambert intensity transmission T = exp(−αL)."""
    return np.exp(-alpha * L)


def optical_density(alpha, L):
    """Base-10 optical density OD = −log10(T) = αL / ln10."""
    return alpha * L / np.log(10.0)


def group_index(chi_phys_axis, detuning_axis, omega0):
    """
    Group index n_g = n + ω dn/dω with n ≈ 1 + Re(χ)/2 (dilute vapor).
    `detuning_axis` is the probe detuning (rad/s); ω ≈ omega0 near resonance.
    Returned per point via a centred gradient of Re(χ).
    """
    n_re = 1.0 + 0.5 * np.real(chi_phys_axis)
    dn_dw = np.gradient(n_re, detuning_axis)
    return n_re + omega0 * dn_dw
