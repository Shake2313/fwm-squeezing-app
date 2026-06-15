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
from .lineshape import fwhm_interp


def _gain_matrix_from_chi(chi_ss_avg, chi_sc_avg, chi_cs_avg, chi_cc_avg,
                          k_probe, k_conj, N_atoms, dipole, line_strength,
                          delta_k_z=None):
    """Build the 2x2 Maxwell-Bloch propagation matrix stack."""
    chi_ss_avg = np.asarray(chi_ss_avg, dtype=complex)
    chi_sc_avg = np.asarray(chi_sc_avg, dtype=complex)
    chi_cs_avg = np.asarray(chi_cs_avg, dtype=complex)
    chi_cc_avg = np.asarray(chi_cc_avg, dtype=complex)
    k_probe = np.asarray(k_probe, dtype=complex)
    k_conj = np.asarray(k_conj, dtype=complex)

    coupling = -2.0 * N_atoms * line_strength * dipole**2 / (constants.EPS_0 * constants.HBAR)
    n = chi_ss_avg.size
    M = np.zeros((n, 2, 2), dtype=complex)
    M[:, 0, 0] = 0.5j * k_probe * coupling * chi_ss_avg
    M[:, 0, 1] = 0.5j * k_probe * coupling * chi_sc_avg
    M[:, 1, 0] = -0.5j * k_conj * coupling * chi_cs_avg.conj()
    M[:, 1, 1] = -0.5j * k_conj * coupling * chi_cc_avg.conj()
    if delta_k_z is not None:
        dk = np.asarray(delta_k_z, dtype=float)
        M[:, 0, 0] += 0.5j * dk
        M[:, 1, 1] -= 0.5j * dk
    return M


def gain_from_chi(chi_ss_avg, chi_sc_avg, chi_cs_avg, chi_cc_avg,
                  k_probe, k_conj, L, N_atoms,
                  dipole=None, line_strength=None, delta_k_z=None,
                  propagation_segments=1, segment_profile=None):
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
    M = _gain_matrix_from_chi(
        chi_ss_avg, chi_sc_avg, chi_cs_avg, chi_cc_avg,
        k_probe, k_conj, N_atoms, dipole, line_strength, delta_k_z=delta_k_z)

    n = M.shape[0]
    nseg = max(int(propagation_segments or 1), 1)
    if nseg <= 1 and segment_profile is None:
        T = matrix_exp_2x2(M, L)
    else:
        if segment_profile is None:
            profile = np.ones(nseg, dtype=float)
        else:
            profile = np.asarray(segment_profile, dtype=float)
            if profile.size != nseg:
                raise ValueError("segment_profile length must match propagation_segments")
        dz = L / nseg
        T = np.broadcast_to(np.eye(2, dtype=complex), (n, 2, 2)).copy()
        for scale in profile:
            Mz = M.copy()
            Mz[:, 0, 1] *= scale
            Mz[:, 1, 0] *= scale
            T = matrix_exp_2x2(Mz, dz) @ T
    G_s = np.abs(T[:, 0, 0]) ** 2
    G_c = np.abs(T[:, 1, 0]) ** 2
    return G_s, G_c, T


def pump_depletion_saturation(G_s, G_c, P_pump, P_seed):
    """
    Energy-conservation (pump-depletion) saturation of the small-signal FWM gains.

    The undepleted-pump linear propagation in `gain_from_chi` returns a
    *small-signal* gain that, at high density, would extract more power than the
    pump can physically supply (e.g. G_s·P_seed ≫ P_pump). Non-degenerate FWM is
    a Manley-Rowe process — two pump photons create one signal + one conjugate
    photon — so at full conversion the seeded signal adds at most half the pump
    power and the generated conjugate the other half:

        (G_s − 1)·P_seed → P_pump/2,    G_c·P_seed → P_pump/2   (high gain).

    A smooth homogeneous-saturation form leaves the small-signal gain untouched
    where (G−1)·P_seed ≪ P_pump and enforces the energy bound at high gain:

        G_s_sat = 1 + (G_s−1) / (1 + (G_s−1)·P_seed / P_cap),   P_cap = P_pump/2
        G_c_sat =      G_c    / (1 +  G_c   ·P_seed / P_cap)

    so (G_s−1) and G_c saturate identically (preserving the twin-beam relation
    G_c ≈ G_s − 1). Returns the saturated (G_s, G_c).
    """
    P_seed = max(float(P_seed), 1e-30)
    P_cap = max(0.5 * float(P_pump), 1e-30)
    gain_part = np.maximum(np.asarray(G_s, dtype=float) - 1.0, 0.0)
    conj = np.maximum(np.asarray(G_c, dtype=float), 0.0)
    G_s_sat = 1.0 + gain_part / (1.0 + gain_part * P_seed / P_cap)
    G_c_sat = conj / (1.0 + conj * P_seed / P_cap)
    return G_s_sat, G_c_sat


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


def _smooth_same(y, x, fwhm):
    if fwhm is None or fwhm <= 0:
        return y
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if x.size < 3:
        return y
    dx = float(np.median(np.diff(x)))
    if dx <= 0:
        return y
    sigma = fwhm / 2.354820045
    half = max(int(np.ceil(4.0 * sigma / dx)), 1)
    grid = np.arange(-half, half + 1) * dx
    kernel = np.exp(-0.5 * (grid / sigma) ** 2)
    kernel /= np.sum(kernel)
    full = np.convolve(y, kernel, mode="full")
    start = (kernel.size - 1) // 2
    return full[start:start + y.size]


def biphoton_stats(tau_axis_ns, waveform, pair_rate_cps, *,
                   signal_eff=0.1, idler_eff=0.1,
                   dark_signal_cps=0.0, dark_idler_cps=0.0,
                   coincidence_window_ns=1.0, timing_jitter_ns=0.0,
                   filter_bandwidth_mhz=None, source_bandwidth_mhz=300.0,
                   target_g2_peak=None):
    """
    Reference-calibrated spontaneous-SFWM biphoton readout.

    `waveform` is the complex velocity-class coherent sum. Rates are a calibrated
    source estimate; detector efficiency, background/dark counts, timing jitter,
    coincidence window and finite filter bandwidth are folded into the returned
    count-rate and correlation observables.
    """
    tau_axis_ns = np.asarray(tau_axis_ns, dtype=float)
    waveform = np.asarray(waveform, dtype=complex)
    intensity = np.abs(waveform) ** 2
    if np.nanmax(intensity) > 0:
        intensity = intensity / np.nanmax(intensity)

    if filter_bandwidth_mhz and filter_bandwidth_mhz > 0:
        filter_transmission = min(1.0, float(filter_bandwidth_mhz)
                                  / max(float(source_bandwidth_mhz), 1e-12))
    else:
        filter_transmission = 1.0
    pair_rate = max(float(pair_rate_cps), 0.0) * filter_transmission

    intensity = _smooth_same(intensity, tau_axis_ns, timing_jitter_ns)
    if np.nanmax(intensity) > 0:
        intensity = intensity / np.nanmax(intensity)

    eta_s = np.clip(float(signal_eff), 0.0, 1.0)
    eta_i = np.clip(float(idler_eff), 0.0, 1.0)
    singles_signal = pair_rate * eta_s + max(float(dark_signal_cps), 0.0)
    singles_idler = pair_rate * eta_i + max(float(dark_idler_cps), 0.0)
    coincidence = pair_rate * eta_s * eta_i
    accidental = singles_signal * singles_idler * max(float(coincidence_window_ns), 0.0) * 1e-9
    raw_accidental = accidental
    raw_car = coincidence / max(raw_accidental, 1e-30)
    added_accidental = 0.0
    if target_g2_peak is not None:
        target_car = max(float(target_g2_peak) - 1.0, 0.0)
        if target_car > 0 and target_car < raw_car:
            accidental = coincidence / target_car
            added_accidental = max(accidental - raw_accidental, 0.0)
            car = target_car
        else:
            car = raw_car
    else:
        car = raw_car
    g2_tau = 1.0 + car * intensity
    g2_peak = float(np.nanmax(g2_tau)) if g2_tau.size else np.nan
    return {
        "g2_SI_tau": g2_tau,
        "tau_axis_ns": tau_axis_ns,
        "fwhm_ns": fwhm_interp(tau_axis_ns, intensity),
        "pair_rate_cps": pair_rate,
        "singles_signal_cps": singles_signal,
        "singles_idler_cps": singles_idler,
        "coincidence_cps": coincidence,
        "accidental_cps": accidental,
        "raw_accidental_cps": raw_accidental,
        "added_accidental_cps": added_accidental,
        "CAR": car,
        "raw_CAR": raw_car,
        "heralding_signal": coincidence / max(singles_idler, 1e-30),
        "heralding_idler": coincidence / max(singles_signal, 1e-30),
        "cauchy_schwarz_R": g2_peak ** 2 / 4.0,
        "g2_peak": g2_peak,
        "raw_g2_peak": 1.0 + raw_car,
        "filter_transmission": filter_transmission,
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


def segmented_loss_noise_squeezing_dB(
        G_s, G_c, eta, *, in_cell_loss_frac=0.0,
        seed_excess_noise=0.0, pump_scatter_noise=0.0,
        eom_residual_noise=0.0):
    """Ultra FWM squeezing with in-cell loss and additive technical noise."""
    G_s = np.asarray(G_s, dtype=float)
    G_c = np.asarray(G_c, dtype=float)
    cross = np.sqrt(np.maximum(G_s * G_c - 1.0, 0.0))
    S_ideal = (G_s + G_c - 2.0 * cross) / np.maximum(G_s + G_c, 1e-30)
    S_ideal = np.clip(S_ideal, 0.0, None)
    tau_cell = 1.0 - np.clip(float(in_cell_loss_frac), 0.0, 1.0)
    S_cell = tau_cell * S_ideal + (1.0 - tau_cell)
    tech = (max(float(seed_excess_noise), 0.0)
            + max(float(pump_scatter_noise), 0.0)
            + max(float(eom_residual_noise), 0.0))
    S = eta * S_cell + (1.0 - eta) + tech
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
