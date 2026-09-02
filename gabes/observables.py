"""Convert reduced susceptibilities into propagation diagnostics.

The module supplies Beer--Lambert absorption, the classical seed/conjugate
transfer matrix, and explicitly labelled gain-referred noise fixtures.  It does
not derive the microscopic Langevin diffusion needed for a physical FWM
squeezing spectrum.
"""
import numpy as np

from . import constants
from .core import matrix_exp_2x2
from .lineshape import fwhm_interp


def _gain_matrix_from_chi(chi_ss_avg, chi_sc_avg, chi_cs_avg, chi_cc_avg,
                          k_probe, k_conj, N_atoms, dipole, line_strength,
                          delta_k_z=None):
    """Build the electric-field Maxwell matrix in the Option-A convention.

    GABES uses positive-frequency fields proportional to ``exp(-i omega t)``
    and defines a passive susceptibility by ``Im(chi_phys) > 0``.  Therefore
    ``dE/dz = +i k chi_phys E / 2`` and a passive diagonal response attenuates
    intensity as ``exp[-k Im(chi_phys) z]``.  Real susceptibility is present
    only in these diagonal terms; ``delta_k_z = 2*k_pump-k_probe-k_conj``
    must be the vacuum/geometric mismatch and must not contain a refractive-
    index correction.  With forward carriers ``exp[i(k*z-omega*t)]``, the
    symmetric amplitudes are ``(E_probe*exp[-i*delta_k*z/2],
    E_conj^* * exp[+i*delta_k*z/2])`` and their mismatch terms are therefore
    ``(-i*delta_k/2, +i*delta_k/2)``.
    """
    chi_ss_avg = np.asarray(chi_ss_avg, dtype=complex)
    chi_sc_avg = np.asarray(chi_sc_avg, dtype=complex)
    chi_cs_avg = np.asarray(chi_cs_avg, dtype=complex)
    chi_cc_avg = np.asarray(chi_cc_avg, dtype=complex)
    k_probe = np.asarray(k_probe, dtype=complex)
    k_conj = np.asarray(k_conj, dtype=complex)

    chi_ss_phys = chi_phys(
        chi_ss_avg, N_atoms, dipole=dipole, line_strength=line_strength)
    chi_sc_phys = chi_phys(
        chi_sc_avg, N_atoms, dipole=dipole, line_strength=line_strength)
    chi_cs_phys = chi_phys(
        chi_cs_avg, N_atoms, dipole=dipole, line_strength=line_strength)
    chi_cc_phys = chi_phys(
        chi_cc_avg, N_atoms, dipole=dipole, line_strength=line_strength)
    n = chi_ss_avg.size
    M = np.zeros((n, 2, 2), dtype=complex)
    M[:, 0, 0] = 0.5j * k_probe * chi_ss_phys
    M[:, 0, 1] = 0.5j * k_probe * chi_sc_phys
    M[:, 1, 0] = -0.5j * k_conj * chi_cs_phys.conj()
    M[:, 1, 1] = -0.5j * k_conj * chi_cc_phys.conj()
    if delta_k_z is not None:
        dk = np.asarray(delta_k_z, dtype=float)
        M[:, 0, 0] -= 0.5j * dk
        M[:, 1, 1] += 0.5j * dk
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

    Returns the raw electric-field-basis power coefficients
    ``G_s = |T₀₀|²`` and ``G_c = |T₁₀|²``, plus the full 2×2 transfer
    matrix stack.  Photon-flux and collected-power interpretations require the
    explicit mode conversion in :func:`canonical_transfer_from_field`.
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


def gaussian_mode_area(waist):
    """Effective area ``pi*w^2/2`` for a Gaussian 1/e^2 intensity radius."""
    waist = np.asarray(waist, dtype=float)
    if np.any(~np.isfinite(waist)) or np.any(waist <= 0.0):
        raise ValueError("Gaussian mode waist must be finite and positive")
    return 0.5 * np.pi * waist**2


def photon_flux_mode_matrix(omega_probe, omega_conj, area_probe, area_conj):
    """Return ``Q`` in ``E = Q a`` for two traveling-wave collected modes.

    ``a = (a_probe, a_conj^dagger)^T`` is photon-flux normalized and

        Q_j = sqrt[2 hbar omega_j / (epsilon_0 c A_j)].

    Frequencies and effective areas may be scalars or broadcast-compatible
    arrays.  Explicit areas prevent an electric-field coefficient from being
    silently interpreted as a photon-flux or power coefficient.
    """
    omega_probe, omega_conj, area_probe, area_conj = np.broadcast_arrays(
        np.asarray(omega_probe, dtype=float),
        np.asarray(omega_conj, dtype=float),
        np.asarray(area_probe, dtype=float),
        np.asarray(area_conj, dtype=float),
    )
    values = (omega_probe, omega_conj, area_probe, area_conj)
    if any(np.any(~np.isfinite(value)) or np.any(value <= 0.0)
           for value in values):
        raise ValueError("Mode frequencies and areas must be finite and positive")
    q_probe = np.sqrt(
        2.0 * constants.HBAR * omega_probe
        / (constants.EPS_0 * constants.C_LIGHT * area_probe))
    q_conj = np.sqrt(
        2.0 * constants.HBAR * omega_conj
        / (constants.EPS_0 * constants.C_LIGHT * area_conj))
    Q = np.zeros(q_probe.shape + (2, 2), dtype=float)
    Q[..., 0, 0] = q_probe
    Q[..., 1, 1] = q_conj
    return Q


def canonical_transfer_from_field(
        transfer_field, omega_probe, omega_conj, area_probe, area_conj):
    """Convert an electric-field transfer to the photon-flux canonical basis."""
    transfer_field = np.asarray(transfer_field, dtype=complex)
    if transfer_field.shape[-2:] != (2, 2):
        raise ValueError("transfer_field must end in a 2x2 matrix")
    Q = photon_flux_mode_matrix(
        omega_probe, omega_conj, area_probe, area_conj)
    leading = np.broadcast_shapes(transfer_field.shape[:-2], Q.shape[:-2])
    transfer_field = np.broadcast_to(transfer_field, leading + (2, 2))
    Q = np.broadcast_to(Q, leading + (2, 2))
    q_probe = Q[..., 0, 0]
    q_conj = Q[..., 1, 1]
    transfer_canonical = np.array(transfer_field, copy=True)
    transfer_canonical[..., 0, 1] *= q_conj / q_probe
    transfer_canonical[..., 1, 0] *= q_probe / q_conj
    return transfer_canonical, Q


def canonical_transfer_diagnostics(
        transfer_field, omega_probe, omega_conj, area_probe, area_conj):
    """Canonical gains, photon-flux gap, and bare commutator defect.

    The defect is diagnostic only: a dissipative cell generally needs explicit
    reservoir channels before the complete input-output map is canonical.
    """
    transfer_canonical, Q = canonical_transfer_from_field(
        transfer_field, omega_probe, omega_conj, area_probe, area_conj)
    J = np.diag([1.0, -1.0]).astype(complex)
    defect = (transfer_canonical @ J
              @ np.swapaxes(transfer_canonical.conj(), -1, -2) - J)
    probe_gain = np.abs(transfer_canonical[..., 0, 0])**2
    conjugate_flux = np.abs(transfer_canonical[..., 1, 0])**2
    omega_probe, omega_conj = np.broadcast_arrays(
        np.asarray(omega_probe, dtype=float),
        np.asarray(omega_conj, dtype=float))
    omega_ratio = np.broadcast_to(
        omega_conj / omega_probe, probe_gain.shape)
    conjugate_power = omega_ratio * conjugate_flux
    return {
        "Q": Q,
        "transfer_canonical": transfer_canonical,
        "probe_power_gain": probe_gain,
        "conjugate_photon_flux_gain": conjugate_flux,
        "conjugate_power_gain": conjugate_power,
        "photon_flux_gap": probe_gain - conjugate_flux,
        "commutator_defect": defect,
        "commutator_defect_max": np.max(np.abs(defect), axis=(-2, -1)),
    }


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
    source estimate; detector efficiency, background/dark counts, the net
    signal-idler timing-difference response, coincidence window and finite filter
    bandwidth are folded into the returned count-rate and correlation observables.

    ``source_fwhm_ns`` is measured from ``|waveform|^2`` before timing smoothing;
    ``detected_fwhm_ns`` is measured after it. The legacy ``fwhm_ns`` key remains
    an alias for the detected width.
    """
    tau_axis_ns = np.asarray(tau_axis_ns, dtype=float)
    waveform = np.asarray(waveform, dtype=complex)
    intensity = np.abs(waveform) ** 2
    if np.nanmax(intensity) > 0:
        intensity = intensity / np.nanmax(intensity)
    source_fwhm_ns = fwhm_interp(tau_axis_ns, intensity)

    if filter_bandwidth_mhz and filter_bandwidth_mhz > 0:
        filter_transmission = min(1.0, float(filter_bandwidth_mhz)
                                  / max(float(source_bandwidth_mhz), 1e-12))
    else:
        filter_transmission = 1.0
    pair_rate = max(float(pair_rate_cps), 0.0) * filter_transmission

    intensity = _smooth_same(intensity, tau_axis_ns, timing_jitter_ns)
    if np.nanmax(intensity) > 0:
        intensity = intensity / np.nanmax(intensity)
    detected_fwhm_ns = fwhm_interp(tau_axis_ns, intensity)

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
        "source_fwhm_ns": source_fwhm_ns,
        "detected_fwhm_ns": detected_fwhm_ns,
        "fwhm_ns": detected_fwhm_ns,
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


def ideal_twin_beam_noise(G_s, G_c):
    """
    Ideal (lossless) twin-beam intensity-difference noise normalized to the SQL.

    For a lossless seeded FWM amplifier the photon-number difference N_s − N_c is
    a conserved quantity of the Bogoliubov transformation (|u|²−|v|²=1): one finds
    the operator identity (N_s − N_c)|_out = (G_s − G_c)·(N_s − N_c)|_in, so with a
    bright coherent seed (conjugate from vacuum) and SQL = ⟨N_s⟩+⟨N_c⟩,

        S_ideal = (G_s − G_c)² / (G_s + G_c).

    In the lossless domain G_s − G_c = 1 and this reduces to the standard
    1/(G_s + G_c) = 1/(2G − 1) (G_s=G, G_c=G−1; Sim et al., Sci. Rep. 15, 7727
    (2025); McCormick et al., Opt. Lett. 32, 178 (2007)). S_ideal ≤ 1 in that
    domain; it is clipped to [0, 1] so passive-loss regions of a probe scan (no
    parametric gain, G_s−G_c arbitrary) cannot spuriously report sub-SQL noise.
    """
    G_s = np.asarray(G_s, dtype=float)
    G_c = np.asarray(G_c, dtype=float)
    total = np.maximum(G_s + G_c, 1e-30)
    return np.clip((G_s - G_c) ** 2 / total, 0.0, 1.0)


def gain_referred_noise_dB(G_s, G_c, eta):
    """Algebraic gain-referred intensity-difference diagnostic in dB.

    This completes a mean-field gain pair with the ideal lossless twin-beam
    identity and symmetric detection efficiency η:
        S_ideal = (G_s − G_c)² / (G_s + G_c)      [see `ideal_twin_beam_noise`]
        S(η)    = η · S_ideal + (1 − η)

    A negative result is *not* a physical squeezing prediction unless the
    frequency-dependent atomic Langevin diffusion and collected-mode covariance
    are supplied independently.  The seeded-FWM path does not yet supply them.
    """
    S_ideal = ideal_twin_beam_noise(G_s, G_c)
    S = eta * S_ideal + (1.0 - eta)
    return 10.0 * np.log10(np.maximum(S, 1e-30))


def intensity_difference_squeezing_dB(G_s, G_c, eta):
    """Backward-compatible alias for :func:`gain_referred_noise_dB`.

    The historical name is retained for callers, but it must not be interpreted
    as a physical squeezing spectrum without a microscopic noise covariance.
    """
    return gain_referred_noise_dB(G_s, G_c, eta)


def balanced_twin_beam_noise(
        G_s, G_c, eta_s=1.0, eta_c=1.0, *, reference_weight=1.0,
        source_noise=None,
        seed_excess_noise=0.0, reference_excess_noise=0.0):
    """
    Twin-beam intensity-difference noise with unequal arm efficiencies.

    The returned value is linear noise power normalized to the coherent-state
    shot-noise level of the same weighted photocurrent difference,

        D = I_s - w I_c,      S = Var(D) / (I_s + w^2 I_c).

    `eta_s` and `eta_c` are total intensity efficiencies after the FWM source
    (cell transmission, path loss and detector QE).  `reference_weight` may be a
    scalar electronic gain, or one of:

    - "raw": w = 1.
    - "dc": w = <I_s>/<I_c>, which cancels the mean photocurrents.
    - "shot": w = sqrt(<I_s>/<I_c>), which minimizes the normalized quantum
      noise for fixed detected powers.

    `source_noise`, when provided, replaces the ideal GABES source covariance
    with a measured/inferred source intensity-difference noise in linear units.
    This is useful for paper-anchored resonant-FWM operating points where the
    compact four-level gain model is not the final source metrology reference.

    For eta_s == eta_c == eta, w == 1 and source_noise is None this reduces
    exactly to `intensity_difference_squeezing_dB(..., eta)` in linear units.
    """
    G_s = np.asarray(G_s, dtype=float)
    G_c = np.asarray(G_c, dtype=float)
    eta_s = np.clip(np.asarray(eta_s, dtype=float), 0.0, 1.0)
    eta_c = np.clip(np.asarray(eta_c, dtype=float), 0.0, 1.0)

    mean_s = eta_s * G_s
    mean_c = eta_c * G_c
    if source_noise is None:
        # Source intensity-difference covariance consistent with the conserved
        # (N_s − N_c) result: at w=1, eta_s=eta_c this gives exactly
        # Var/SQL = (G_s − G_c)²/(G_s + G_c) = `ideal_twin_beam_noise` (lossless
        # → 1/(2G−1)). cov0 = [(G_s + G_c) − (G_s − G_c)²]/2 (= G_c at G_s−G_c=1),
        # clipped to the Cauchy-Schwarz bound √(G_s G_c) and ≥ 0.
        cov0 = 0.5 * ((G_s + G_c) - (G_s - G_c) ** 2)
        cov0 = np.clip(cov0, 0.0, np.sqrt(np.maximum(G_s * G_c, 0.0)))
    else:
        S0 = np.clip(np.asarray(source_noise, dtype=float), 0.0, None)
        cov0 = 0.5 * (G_s + G_c) * (1.0 - S0)
        cov0 = np.clip(cov0, 0.0, np.sqrt(np.maximum(G_s * G_c, 0.0)))
    cov = eta_s * eta_c * cov0

    if isinstance(reference_weight, str):
        mode = reference_weight.lower()
        if mode == "raw":
            w = np.asarray(1.0)
        elif mode == "dc":
            w = mean_s / np.maximum(mean_c, 1e-30)
        elif mode in ("shot", "optimal", "opt"):
            w = np.sqrt(mean_s / np.maximum(mean_c, 1e-30))
        else:
            raise ValueError("reference_weight must be scalar, raw, dc, or shot")
    else:
        w = np.asarray(reference_weight, dtype=float)

    var = mean_s + w * w * mean_c - 2.0 * w * cov
    shot = mean_s + w * w * mean_c
    tech = (max(float(seed_excess_noise), 0.0)
            + max(float(reference_excess_noise), 0.0))
    return np.maximum(var / np.maximum(shot, 1e-30) + tech, 1e-30)


def balanced_twin_beam_squeezing_dB(
        G_s, G_c, eta_s=1.0, eta_c=1.0, *, reference_weight=1.0,
        source_noise=None,
        seed_excess_noise=0.0, reference_excess_noise=0.0):
    """dB wrapper for `balanced_twin_beam_noise`."""
    S = balanced_twin_beam_noise(
        G_s, G_c, eta_s, eta_c, reference_weight=reference_weight,
        source_noise=source_noise,
        seed_excess_noise=seed_excess_noise,
        reference_excess_noise=reference_excess_noise)
    return 10.0 * np.log10(S)


def segmented_loss_noise_squeezing_dB(
        G_s, G_c, eta, *, in_cell_loss_frac=0.0,
        seed_excess_noise=0.0, pump_scatter_noise=0.0,
        eom_residual_noise=0.0):
    """Legacy gain-referred loss/noise fixture, not a microscopic prediction."""
    S_ideal = ideal_twin_beam_noise(G_s, G_c)
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


def refractive_index(chi_phys_axis):
    """*Phase* refractive index n ≈ 1 + Re(χ)/2 (dilute vapor, |χ| ≪ 1).

    The single place the dilute-vapor convention lives: `group_index` and every
    scheme-level dispersion readout call this instead of re-deriving n.
    """
    return 1.0 + 0.5 * np.real(chi_phys_axis)


def single_pass_phase(chi_phys_axis, k, L):
    """Single-pass phase shift relative to vacuum, φ = k·(n − 1)·L  [rad].

    A *phase* quantity (not a group delay); pair it with `group_index` when the
    envelope velocity is what matters.
    """
    return k * (refractive_index(chi_phys_axis) - 1.0) * L


def group_index(chi_phys_axis, detuning_axis, omega0):
    """
    Group index n_g = n + ω dn/dω with n ≈ 1 + Re(χ)/2 (dilute vapor).
    `detuning_axis` is the probe detuning (rad/s); ω ≈ omega0 near resonance.
    Returned per point via a centred gradient of Re(χ).
    """
    n_re = refractive_index(chi_phys_axis)
    dn_dw = np.gradient(n_re, detuning_axis)
    return n_re + omega0 * dn_dw
