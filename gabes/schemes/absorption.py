"""
Cluster A — weak-probe absorption / dispersion schemes.

OD  : 2-level Doppler-broadened absorption (the validation backbone).
EIT : Λ system, weak probe + strong coupling → transparency window + dispersion.
AT  : same Λ in the strong-coupling regime → Autler-Townes doublet (split ≈ Ω_c).
CPT : same Λ, narrow two-photon scan → sub-natural dark resonance.

All share one engine: build a rotating-frame H₀ with a weak probe (and optional
strong coupling), solve the steady state per velocity class (core.steady_state_
batched), and read the probe coherence ρ_eg/Ω_p = χ̄. Co-propagating geometry, so
the single excited-state Doppler shift (atom.S_v) is exact: the optical line is
Doppler-broadened while the Λ two-photon resonance stays Doppler-free.

Spectroscopic knobs are in units of Γ (natural linewidth) so the physics is
atom-agnostic and validation is direct (AT splitting reads Ω_c). Temperature sets
the Doppler width and the 85Rb number density (hence the absorption scale).
"""
import numpy as np

from .. import atoms, constants, doppler, hyperfine, observables, species
from ..constants import GAMMA, K_VEC, OMEGA_D1
from .. import core, kernels
from ..lineshape import fwhm_halfmax, window_fwhm
from ..report import derived_table
from .base import ParamSpec, Preset, Scheme

PROBE_RABI = 1e-3              # weak probe, in units of Γ
GAMMA_MHZ = GAMMA / (2 * np.pi) / 1e6
MHZ = 2 * np.pi * 1e6
KHZ = 2 * np.pi * 1e3

# numpy ≥ 2.0 renamed np.trapz → np.trapezoid; support both.
_trapz = getattr(np, "trapezoid", getattr(np, "trapz", None))


_TABLE_STEP = GAMMA / 30.0          # Δ_eff sampling for the Doppler χ̄ table


def _ne_buffer_gamma(params):
    return constants.neon_buffer_broadening(params.get("ne_pressure_torr", 0.0))


def _affine_scan_coeffs(atom, build_H0):
    """Constant pieces of L(s, kv) = base + s·A_coef + kv·B_coef.

    Valid only when `build_H0(s)` is affine in the scan variable s — true for
    the Λ (two-photon detuning) and Rydberg ladders. The s-independent drive
    entries are bit-identical between build_H0(0) and build_H0(1), so the
    difference isolates dL/ds exactly with no cancellation. The optical Doppler
    shift enters separately as −Δ_eff·S_v with Δ_eff = s − kv, hence
    A_coef = dL/ds − S_v (the s part) and B_coef = S_v (the +k·v part).
    """
    base = core.build_liouvillian(build_H0(0.0), atom)
    dL_ds = core.build_liouvillian(build_H0(1.0), atom) - base
    return base, dL_ds - atom.S_v, atom.S_v


def _solve_chi_avg(atom, build_H0, probe_coh, scan, params, h_dep,
                   probe_omega=None, k_vec=None, mass=None):
    """
    χ̄(scan) = probe coherence ρ_eg / Ω_p, Doppler-averaged over the Maxwell
    velocity distribution (single excited-state shift, co-propagating geometry).

    Doppler off: a single v = 0 solve per scan point.
    Doppler on, scan-independent H (OD): one fine χ̄(Δ_eff) table + interpolated
        Doppler average — decouples velocity sampling from the solve (accurate
        Voigt without millions of solves), mirroring the FWM Δ_eff trick.
    Doppler on, scan-dependent H (Λ): per scan point, solve directly at the
        Maxwell velocity classes. The Λ two-photon feature is Doppler-free
        (exact); only the broad optical background is averaged.
    """
    e, g = probe_coh
    om_p = PROBE_RABI * GAMMA if probe_omega is None else probe_omega
    k_vec = K_VEC if k_vec is None else k_vec
    mass = constants.MASS_85RB if mass is None else mass
    n = atom.n_levels
    doppler_on = params.get("doppler", "off") == "on"

    if not doppler_on:
        if not h_dep:
            L0 = core.build_liouvillian(build_H0(0.0), atom)
            rho = core.steady_state_batched(L0, scan, atom.S_v, n)
            return rho[:, e, g] / om_p
        # Scan-dependent H, single (v = 0) class: Δ_eff = s. Fold the scan loop
        # into the affine kernel (kv = 0, unit weight) when available.
        if kernels.available():
            base, A_coef, B_coef = _affine_scan_coeffs(atom, build_H0)
            with core.blas_single_thread():
                out = kernels.affine_scan_chi(
                    base, A_coef, B_coef,
                    np.ascontiguousarray(scan, dtype=float),
                    np.zeros(1), np.ones(1), e * n + g, n)
            return out / om_p
        out = np.zeros(scan.size, dtype=complex)
        for i, s in enumerate(scan):
            L0 = core.build_liouvillian(build_H0(s), atom)
            rho = core.steady_state_batched(L0, np.array([s]), atom.S_v, n)
            out[i] = rho[0, e, g] / om_p
        return out

    T = params["temp_c"] + 273.15

    if not h_dep:
        # OD: scan-independent H solved once on a fine Δ_eff table, then averaged.
        # The solve cost is independent of the velocity grid, so keep it fine
        # (dv = 1) — this path feeds the AutoOD-validated absolute scale.
        v, w = doppler.velocity_grid(T, mass=mass, dv=1.0, cutoff_sigma=4.0)
        kvmax = k_vec * np.abs(v).max()
        L0 = core.build_liouvillian(build_H0(0.0), atom)
        lo, hi = scan.min() - kvmax, scan.max() + kvmax
        deff = np.linspace(lo, hi, int((hi - lo) / _TABLE_STEP) + 2)
        rho = core.steady_state_batched(L0, deff, atom.S_v, n)
        table = rho[:, e, g] / om_p
        return doppler.doppler_average_1d(table, deff, scan, v, w, k_vec=k_vec)

    # Scan-dependent H (Λ): the two-photon feature is Doppler-free, so the optical
    # background must be solved at the Maxwell classes per scan point — this path
    # is solve-bound (scan_points × velocity classes linear solves). The Doppler-
    # free feature is insensitive to velocity sampling and the broad optical
    # background is smooth, so a coarser Maxwell quadrature (dv = 2 m/s vs 1)
    # ~halves the solves while keeping transmission < 0.02% and the Re χ / group-
    # index readout < 1% against the fine grid.
    v, w = doppler.velocity_grid(T, mass=mass, dv=2.0, cutoff_sigma=4.0)
    kv = k_vec * v
    if kernels.available():
        base, A_coef, B_coef = _affine_scan_coeffs(atom, build_H0)
        with core.blas_single_thread():
            out = kernels.affine_scan_chi(
                base, A_coef, B_coef,
                np.ascontiguousarray(scan, dtype=float),
                np.ascontiguousarray(kv, dtype=float),
                np.ascontiguousarray(w, dtype=float), e * n + g, n)
        return out / om_p
    out = np.zeros(scan.size, dtype=complex)
    for i, s in enumerate(scan):
        L0 = core.build_liouvillian(build_H0(s), atom)
        rho = core.steady_state_batched(L0, s - kv, atom.S_v, n)
        out[i] = ((rho[:, e, g] / om_p) * w).sum()
    return out


def _hyperfine_alpha(scan, params):
    """
    Realistic 85Rb D1 absorption coefficient α(scan) [1/m]: an incoherent sum of
    the four hyperfine transitions (Fg∈{2,3}→Fe∈{2,3}), each a Doppler-broadened
    Voigt placed at its line-center detuning, weighted by the validated
    Clebsch-Gordan strength C_F² and ground population p_F, with self-broadening.

    The Voigt *shape* is produced once by the OBE Doppler kernel (a genuine OBE
    solve — no analytic wofz) on a fine grid, normalised to unit area, then
    shifted/scaled per line. The per-line integrated absorption reproduces the
    lab-validated AutoOD scaling:
        ∫α_F dδ = π · k · p_F · C_F² · |d|² · N / (ε₀ℏ) / (2(2I+1)).

    Returns (alpha, components, info) where `components` maps (Fg,Fe)→α_F(scan)
    and `info` carries N, Γ_eff, σ_v, Doppler FWHM for the readout.
    """
    T = params["temp_c"] + 273.15
    doppler_on = params.get("doppler", "on") == "on"
    ls = params.get("line_strength", 1.0)

    N = hyperfine.number_density(T)
    buffer_gamma = _ne_buffer_gamma(params)
    gamma_eff = hyperfine.self_broadened_gamma(N) + buffer_gamma
    atom = atoms.two_level(gamma=gamma_eff)

    def build_H0(_s):
        H = np.zeros((2, 2), dtype=complex)
        H[0, 1] = H[1, 0] = PROBE_RABI * GAMMA / 2
        return H

    sigma_v = np.sqrt(constants.KB * T / constants.MASS_85RB)
    dopp_fwhm = np.sqrt(8 * np.log(2)) * K_VEC * sigma_v
    width = dopp_fwhm if doppler_on else gamma_eff

    shifts = {k: 2 * np.pi * v for k, v in hyperfine.LINE_SHIFT_HZ.items()}
    smin, smax = min(shifts.values()), max(shifts.values())

    # One unit-area Voigt shape on a fine grid covering (scan − shift) ± wings.
    lo = scan.min() - smax - 12 * width
    hi = scan.max() - smin + 12 * width
    grid = np.linspace(lo, hi, int((hi - lo) / (width / 40.0)) + 2)
    chi = _solve_chi_avg(atom, build_H0, (1, 0), grid, params, h_dep=False)
    shape, _ = observables.absorption_coefficient(chi, K_VEC, 1.0, line_strength=1.0)
    S_unit = shape / _trapz(shape, grid)          # ∫ S_unit dδ = 1

    # Absolute per-line integrated absorption (AutoOD normalisation).
    K = (np.pi * K_VEC * hyperfine.DIPOLE_SQ * N
         / (constants.HBAR * constants.EPS_0) / hyperfine.N_GROUND_SUBLEVELS)

    alpha = np.zeros_like(scan)
    components = {}
    for (Fg, Fe) in hyperfine.TRANSITIONS:
        I_line = ls * K * hyperfine.GROUND_POP[Fg] * hyperfine.CF2[(Fg, Fe)]
        line = I_line * np.interp(scan - shifts[(Fg, Fe)], grid, S_unit,
                                  left=0.0, right=0.0)
        alpha += line
        components[(Fg, Fe)] = line

    info = dict(N=N, gamma_eff=gamma_eff, buffer_gamma=buffer_gamma,
                sigma_v=sigma_v, dopp_fwhm=dopp_fwhm)
    return alpha, components, info


# =========================================================
# OD — 2-level Doppler-broadened absorption
# =========================================================
class ODScheme(Scheme):
    name = "od"
    cluster = "A — Absorption"
    title = "Optical absorption (OD)"
    caption = ("Weak-probe absorption. The full 85Rb D1 four-line hyperfine "
               "spectrum (validated against the lab AutoOD calculator), or a single "
               "bare 2-level Voigt. Doppler-broadened in a warm vapor; natural-"
               "linewidth Lorentzian when Doppler is off.")

    def param_schema(self):
        return [
            ParamSpec("model", "Absorption model", "Model", "85Rb D1 hyperfine",
                      choices=("85Rb D1 hyperfine", "single 2-level"),
                      help="Hyperfine: the full 85Rb D1 four-line spectrum (CG "
                      "strengths, ground populations, self-broadening) validated "
                      "against the lab AutoOD calculator. Single: one bare 2-level "
                      "Voigt — the validation backbone the Λ schemes reduce to."),
            ParamSpec("temp_c", "Temperature", "Cell & beams", 50.0, 20.0, 200.0, 1.0, "°C",
                      help="Sets the 85Rb density (absorption scale) and Doppler width."),
            ParamSpec("cell_mm", "Cell length", "Cell & beams", 10.0, 0.5, 200.0, 0.5, "mm",
                      recompute=False),
            ParamSpec("line_strength", "Line-strength factor", "Detection & scaling", 1.0,
                      0.01, 2.0, 0.01, "", help="Effective |d|² calibration knob. Single: "
                      "=1.0 reproduces the textbook 3λ²/2π cross-section. Hyperfine: "
                      "=1.0 reproduces the validated AutoOD absolute scale."),
            ParamSpec("doppler", "Doppler (vapor motion)", "Numerics", "on",
                      choices=("on", "off"), advanced=True),
            ParamSpec("ne_pressure_torr", "Ne buffer pressure", "Cell & beams", 0.0,
                      0.0, 200.0, 1.0, "Torr", advanced=True,
                      help="Fixed-neon pressure broadening only; pressure shift and "
                      "Dicke narrowing are not included."),
        ]

    def presets(self):
        return [
            Preset("85Rb D1 cell (90 °C, 12.5 mm)",
                   values=dict(model="85Rb D1 hyperfine", temp_c=90.0, cell_mm=12.5),
                   icon="🧪", help="AutoOD-validated warm-cell spectrum."),
            Preset("Single line (cold)",
                   values=dict(model="single 2-level", temp_c=25.0, cell_mm=3.0,
                               doppler="off"), icon="📏"),
        ]

    # ---- hyperfine (validated full-D1) ----
    def _compute_hyperfine(self, params):
        T = params["temp_c"] + 273.15
        sigma_v = np.sqrt(constants.KB * T / constants.MASS_85RB)
        dopp_fwhm = np.sqrt(8 * np.log(2)) * K_VEC * sigma_v
        N = hyperfine.number_density(T)
        gamma_eff = hyperfine.self_broadened_gamma(N) + _ne_buffer_gamma(params)
        width = dopp_fwhm if params["doppler"] == "on" else gamma_eff

        shifts = np.array([2 * np.pi * v for v in hyperfine.LINE_SHIFT_HZ.values()])
        margin = max(2 * np.pi * 1.5e9, 6 * width)
        lo, hi = shifts.min() - margin, shifts.max() + margin
        n = int(np.clip((hi - lo) / (width / 8.0), 1201, 8000))
        scan = np.linspace(lo, hi, n)

        alpha, components, info = _hyperfine_alpha(scan, params)
        return dict(model="hyperfine", scan=scan, alpha=alpha,
                    components={f"{Fg}-{Fe}": v for (Fg, Fe), v in components.items()},
                    L=params["cell_mm"] * 1e-3, T=T, **info)

    def _observables_hyperfine(self, raw, params):
        import matplotlib.pyplot as plt
        x = raw["scan"] / (2 * np.pi) / 1e9                      # GHz
        alpha = raw["alpha"]
        L = params["cell_mm"] * 1e-3            # navigate-only knob: read live
        OD = observables.optical_density(alpha, L)
        T_trans = observables.transmission(alpha, L)

        fig, (axT, axOD) = plt.subplots(2, 1, figsize=(8.5, 6.4), sharex=True)
        axT.plot(x, T_trans, color="#1f77b4", lw=1.6)
        axT.set_ylabel("Transmission")
        axT.set_ylim(-0.02, 1.04)
        axT.set_title(f"85Rb D1 hyperfine:  T = {params['temp_c']:.0f} °C,  "
                      f"L = {params['cell_mm']:.1f} mm,  Doppler {params['doppler']}")
        axOD.plot(x, OD, color="#d62728", lw=1.6)
        axOD.set_ylabel("Optical density  (−log₁₀T)")
        axOD.set_xlabel("Detuning  [GHz]  (ref: 87Rb F=2→F′=2)")
        for (Fg, Fe), sh in hyperfine.LINE_SHIFT_HZ.items():
            for a in (axT, axOD):
                a.axvline(sh / 1e9, color="gray", ls=":", lw=0.7)
            axOD.annotate(f"{Fg}→{Fe}′", (sh / 1e9, 0), xytext=(0, 2),
                          textcoords="offset points", ha="center", va="bottom",
                          fontsize=7, color="gray")
        fig.tight_layout()

        buffer_mhz = raw.get("buffer_gamma", 0.0) / (2 * np.pi) / 1e6
        self_gamma = raw["gamma_eff"] - GAMMA - raw.get("buffer_gamma", 0.0)
        sb_mhz = self_gamma / (2 * np.pi) / 1e6
        metrics = [
            dict(label="Peak OD", value=f"{np.nanmax(OD):.3f}",
                 help="Largest −log₁₀(T) across the four-line spectrum."),
            dict(label="Min transmission", value=f"{np.nanmin(T_trans):.3e}"),
            dict(label="Doppler FWHM", value=f"{raw['dopp_fwhm']/(2*np.pi)/1e6:.0f} MHz",
                 help="Gaussian (Doppler) width of each hyperfine line."),
        ]
        rows = "".join(
            f"| {Fg}→{Fe}′ | {hyperfine.LINE_SHIFT_HZ[(Fg, Fe)]/1e9:.3f} | "
            f"{hyperfine.GROUND_POP[Fg]*hyperfine.CF2[(Fg, Fe)]:.4f} |\n"
            for (Fg, Fe) in hyperfine.TRANSITIONS)
        derived = derived_table([
            ("N(85Rb), pure cell", f"{raw['N']:.3e} /m³"),
            ("σ_v (1-D)", f"{raw['sigma_v']:.1f} m/s"),
            ("Doppler FWHM", f"{raw['dopp_fwhm']/(2*np.pi)/1e6:.1f} MHz"),
            ("Γ_eff / 2π (self + Ne broadened)",
             f"{raw['gamma_eff']/(2*np.pi)/1e6:.3f} MHz (+{sb_mhz:.3f})"),
            ("Ne buffer pressure", f"{params.get('ne_pressure_torr', 0.0):.0f} Torr"),
            ("Ne broadening / 2π", f"{buffer_mhz:.3f} MHz"),
            ("Natural Γ / 2π", f"{GAMMA_MHZ:.3f} MHz"),
        ])
        lines = ("| Line | Center [GHz] | p_F·C_F² |\n|---|---|---|\n" + rows)
        return dict(metrics=metrics, figure=fig, tables=[
            derived,
            {"title": "Hyperfine lines (relative strength)", "markdown": lines}])

    # ---- single bare 2-level line (validation backbone) ----
    def _compute_single(self, params):
        T = params["temp_c"] + 273.15
        sigma_v = np.sqrt(constants.KB * T / constants.MASS_85RB)
        dopp_fwhm = np.sqrt(8 * np.log(2)) * K_VEC * sigma_v
        buffer_gamma = _ne_buffer_gamma(params)
        gamma_eff = GAMMA + buffer_gamma
        half = (max(10 * gamma_eff, 3.5 * dopp_fwhm)
                if params["doppler"] == "on" else 12 * gamma_eff)
        scan = np.linspace(-half, half, 601)

        atom = atoms.two_level(gamma=gamma_eff)

        def build_H0(_s):
            H = np.zeros((2, 2), dtype=complex)
            H[0, 1] = H[1, 0] = PROBE_RABI * GAMMA / 2
            return H

        chi_bar = _solve_chi_avg(atom, build_H0, (1, 0), scan, params, h_dep=False)
        N = atoms.rb85_density(T)
        return dict(model="single", scan=scan, chi_bar=chi_bar, N=N, T=T,
                    L=params["cell_mm"] * 1e-3, ls=params["line_strength"],
                    sigma_v=sigma_v, dopp_fwhm=dopp_fwhm,
                    gamma_eff=gamma_eff, buffer_gamma=buffer_gamma)

    def _observables_single(self, raw, params):
        import matplotlib.pyplot as plt
        x = raw["scan"] / (2 * np.pi) / 1e6                      # MHz
        alpha, _ = observables.absorption_coefficient(
            raw["chi_bar"], K_VEC, raw["N"], line_strength=raw["ls"])
        L = params["cell_mm"] * 1e-3           # navigate-only knob: read live
        OD = observables.optical_density(alpha, L)
        T_trans = observables.transmission(alpha, L)

        fig, (axT, axOD) = plt.subplots(2, 1, figsize=(8.5, 6.4), sharex=True)
        axT.plot(x, T_trans, color="#1f77b4", lw=1.8)
        axT.set_ylabel("Transmission")
        axT.set_ylim(-0.02, 1.02)
        axT.set_title(f"Single 2-level:  T = {params['temp_c']:.0f} °C,  "
                      f"L = {params['cell_mm']:.0f} mm,  Doppler {params['doppler']}")
        axOD.plot(x, OD, color="#d62728", lw=1.8)
        axOD.set_ylabel("Optical density  (−log₁₀T)")
        axOD.set_xlabel("Probe detuning  [MHz]")
        for a in (axT, axOD):
            a.axvline(0, color="gray", ls=":", lw=0.8)
        fig.tight_layout()

        fwhm_mhz = fwhm_halfmax(x, OD)
        buffer_mhz = raw.get("buffer_gamma", 0.0) / (2 * np.pi) / 1e6
        metrics = [
            dict(label="Peak OD", value=f"{np.nanmax(OD):.3f}",
                 help="On-line optical density −log₁₀(T)."),
            dict(label="Min transmission", value=f"{np.nanmin(T_trans):.3f}"),
            dict(label="Line FWHM", value=f"{fwhm_mhz:.1f} MHz",
                 help="Voigt width (Doppler on) or natural Γ (Doppler off)."),
        ]
        derived = derived_table([
            ("Ne buffer pressure", f"{params.get('ne_pressure_torr', 0.0):.0f} Torr"),
            ("Ne broadening / 2π", f"{buffer_mhz:.3f} MHz"),
            ("Γ_eff / 2π", f"{raw.get('gamma_eff', GAMMA)/(2*np.pi)/1e6:.3f} MHz"),
            ("N(85Rb)", f"{raw['N']:.3e} /m³"),
            ("σ_v (1-D)", f"{raw['sigma_v']:.1f} m/s"),
            ("Doppler FWHM", f"{raw['dopp_fwhm']/(2*np.pi)/1e6:.1f} MHz"),
            ("Natural Γ / 2π", f"{GAMMA_MHZ:.3f} MHz"),
        ])
        return dict(metrics=metrics, figure=fig, tables=[derived])

    def compute(self, params):
        if params.get("model", "85Rb D1 hyperfine") == "single 2-level":
            return self._compute_single(params)
        return self._compute_hyperfine(params)

    def observables(self, raw, params):
        if raw.get("model") == "single":
            return self._observables_single(raw, params)
        return self._observables_hyperfine(raw, params)


LAMBDA_GENERIC = "Generic"


def _medium_from_params(params, default_temp_c=50.0):
    """Scalar D-line medium used by the Lambda/Rydberg linear readout models."""
    T = float(params.get("temp_c", default_temp_c)) + 273.15
    line = params.get("line", "D1")
    key = params.get("species", "Rb (natural)")
    if key == LAMBDA_GENERIC:
        return dict(
            key=key, label=key, line=line, T=T, N=atoms.rb85_density(T),
            gamma=GAMMA, gamma_mhz=GAMMA_MHZ, k_vec=K_VEC, omega0=OMEGA_D1,
            dipole=constants.DIPOLE_D1, mass=constants.MASS_85RB,
        )

    comps = species.SPECIES.get(key, species.SPECIES["Rb (natural)"])
    iso_ref = max(comps, key=lambda item: item[1])[0]
    Je, nu0, gamma_mhz, _, _ = iso_ref.line(line)
    gamma = gamma_mhz * MHZ
    lam = constants.C_LIGHT / nu0
    total_weight = sum(weight for _, weight in comps)
    N = sum(species.number_density(iso, T) * weight for iso, weight in comps)
    mass = sum(iso.mass * weight for iso, weight in comps) / max(total_weight, 1e-30)
    dipole = np.sqrt(species.reduced_dipole_sq(gamma, lam, iso_ref.Jg, Je))
    return dict(
        key=key, label=key, line=line, T=T, N=N, gamma=gamma,
        gamma_mhz=gamma_mhz, k_vec=2 * np.pi / lam, omega0=2 * np.pi * nu0,
        dipole=dipole, mass=mass,
    )


def _mkey(params, new_key, old_key, scale, default):
    """Read a physical-unit key, with old Gamma-unit Lambda params as fallback."""
    if new_key in params:
        return float(params[new_key])
    if old_key in params:
        return float(params[old_key]) * scale
    return default


# =========================================================
# Λ system — EIT / AT / CPT (one engine, three presentations)
# =========================================================
class LambdaScheme(Scheme):
    cluster = "A — Absorption"
    cache_version = "physical-units-v2"
    defaults_version = "physical-units-v2"

    _TITLE = {
        "eit": "Electromagnetically induced transparency (EIT)",
        "at":  "Autler-Townes splitting (AT)",
        "cpt": "Coherent population trapping (CPT)",
    }
    _CAPTION = {
        "eit": ("Λ system: a strong coupling field opens a transparency "
                "window for a weak probe at two-photon resonance."),
        "at":  ("Strong coupling dresses the excited state; the weak-probe "
                "absorption line splits into a doublet separated by ≈ Ω_c."),
        "cpt": ("Narrow two-photon (Raman) scan: atoms are pumped into a "
                "dark state, giving a sub-natural absorption dip."),
    }

    # Per-regime defaults converted from the old Gamma-unit controls into lab units.
    _DEF = {
        "eit": dict(oc_mhz=3.0 * GAMMA_MHZ, gg_khz=0.01 * GAMMA_MHZ * 1e3,
                    temp=50.0, cell=15.0, dopp="on"),
        "at":  dict(oc_mhz=8.0 * GAMMA_MHZ, gg_khz=0.01 * GAMMA_MHZ * 1e3,
                    temp=25.0, cell=3.0, dopp="off"),
        "cpt": dict(oc_mhz=1.0 * GAMMA_MHZ, gg_khz=0.005 * GAMMA_MHZ * 1e3,
                    temp=25.0, cell=3.0, dopp="off"),
    }

    def __init__(self, mode=None):
        self.mode = mode                       # None → merged EIT/AT/CPT dropdown entry
        if mode is None:
            self.name = "lambda"
            self.title = "Λ coherence (EIT / AT / CPT)"
            self.caption = ("One Λ system, three regimes of the same coupling field Ω_c: a "
                            "weak coupling opens an EIT transparency window, a strong one "
                            "splits the line into an Autler-Townes doublet, and a narrow "
                            "two-photon scan reveals the CPT dark resonance.")
        else:
            self.name = mode
            self.title = self._TITLE[mode]
            self.caption = self._CAPTION[mode]

    def _mode(self, params):
        """Effective regime: pinned (alias instance) or the merged `view` knob."""
        return self.mode or str(params.get("view", "eit")).lower()

    def param_schema(self):
        d = self._DEF[self.mode or "eit"]
        specs = []
        if self.mode is None:
            specs.append(ParamSpec(
                "view", "Regime", "Regime", "EIT", choices=("EIT", "AT", "CPT"),
                control="segmented", applies_defaults=True,
                help="EIT: weak-coupling transparency window. AT: strong-coupling "
                     "doublet (split ≈ Ω_c). CPT: narrow two-photon dark resonance. "
                     "One Λ engine — the presets set textbook values for each."))
        specs += [
            ParamSpec("species", "Atom / isotope", "Atomic", "Rb (natural)",
                      choices=tuple(species.SPECIES_ORDER) + (LAMBDA_GENERIC,),
                      help="Scalar D-line medium. Natural Rb uses abundance-weighted density."),
            ParamSpec("line", "Transition line", "Atomic", "D1",
                      choices=("D1", "D2"),
                      help="Sets linewidth, wavelength, Doppler scale and dipole moment."),
            ParamSpec("coupling_detuning_mhz", "Coupling detuning", "Detunings", 0.0,
                      -80.0, 80.0, 0.1, "MHz"),
            ParamSpec("coupling_rabi_mhz", "Coupling Rabi", "Fields", d["oc_mhz"],
                      0.1, 120.0, 0.1, "MHz",
                      help="Strong dressing-field Rabi frequency on the second leg.",
                      endpoints=("weak EIT", "strong AT")),
            ParamSpec("buffer_ground_relax_khz", "Buffer ground relaxation", "Atomic",
                      d["gg_khz"], 0.0, 2000.0, 1.0, "kHz",
                      help="Ground-state coherence relaxation (buffer-gas collisions); "
                           "sets the EIT/CPT linewidth floor. Same quantity as the "
                           "magneto buffer ground relaxation."),
            ParamSpec("temp_c", "Temperature", "Cell & beams", d["temp"], 20.0, 200.0, 1.0, "°C"),
            ParamSpec("cell_mm", "Cell length", "Cell & beams", d["cell"], 0.5, 200.0, 0.5, "mm",
                      recompute=False),
            ParamSpec("line_strength", "Line-strength factor", "Detection & scaling", 1.0,
                      0.01, 2.0, 0.01, "", advanced=True),
            ParamSpec("doppler", "Doppler (vapor motion)", "Numerics", d["dopp"],
                      choices=("on", "off"), advanced=True),
        ]
        return specs

    def presets(self):
        if self.mode is None:
            return []
        label = self._TITLE[self.mode].split("(")[-1].rstrip(")")
        return [Preset(f"{label} default", values=self._default_values(self.mode))]

    def _default_values(self, mode):
        d = self._DEF[mode]
        view = {"eit": "EIT", "at": "AT", "cpt": "CPT"}[mode]
        return dict(
            view=view, coupling_detuning_mhz=0.0,
            coupling_rabi_mhz=d["oc_mhz"], buffer_ground_relax_khz=d["gg_khz"],
            temp_c=d["temp"], cell_mm=d["cell"], doppler=d["dopp"],
            line_strength=1.0,
        )

    def recommended_defaults(self, params):
        if self.mode is not None:
            return {self._TITLE[self.mode].split("(")[-1].rstrip(")"):
                    self._default_values(self.mode)}
        return {
            "EIT": self._default_values("eit"),
            "AT": self._default_values("at"),
            "CPT": self._default_values("cpt"),
        }

    def _scan(self, params):
        m = self._mode(params)
        medium = _medium_from_params(params)
        gamma_mhz = medium["gamma_mhz"]
        Oc = _mkey(params, "coupling_rabi_mhz", "coupling_rabi", gamma_mhz,
                   self._DEF[m]["oc_mhz"])
        Dc = _mkey(params, "coupling_detuning_mhz", "coupling_detuning", gamma_mhz, 0.0)
        gg = _mkey(params, "buffer_ground_relax_khz", "gamma_gg", gamma_mhz * 1e3,
                   self._DEF[m]["gg_khz"]) / (gamma_mhz * 1e3)
        Oc_g = Oc / gamma_mhz
        Dc_g = Dc / gamma_mhz
        if m == "cpt":
            dark = gg + Oc_g**2          # power-broadened dark width (Gamma units)
            half = max(0.05, 12.0 * dark)
        elif m == "at":
            half = max(10.0, 3.0 * Oc_g)
        else:
            half = max(8.0, 4.0 * Oc_g)
        return (Dc_g + np.linspace(-half, half, 601)) * medium["gamma"]

    def compute(self, params):
        medium = _medium_from_params(params)
        gamma = medium["gamma"]
        gamma_mhz = medium["gamma_mhz"]
        Oc_mhz = _mkey(params, "coupling_rabi_mhz", "coupling_rabi", gamma_mhz,
                       self._DEF[self._mode(params)]["oc_mhz"])
        Dc_mhz = _mkey(params, "coupling_detuning_mhz", "coupling_detuning", gamma_mhz, 0.0)
        gg_khz = _mkey(params, "buffer_ground_relax_khz", "gamma_gg", gamma_mhz * 1e3,
                       self._DEF[self._mode(params)]["gg_khz"])
        atom = atoms.lambda3(gamma_gg=gg_khz * KHZ, gamma=gamma)
        Oc = Oc_mhz * MHZ
        Dc = Dc_mhz * MHZ
        probe = PROBE_RABI * gamma

        def build_H0(s):
            H = np.zeros((3, 3), dtype=complex)
            H[1, 1] = Dc - s                      # two-photon detuning (Doppler-free)
            H[0, 2] = H[2, 0] = probe / 2
            H[1, 2] = H[2, 1] = Oc / 2
            return H

        scan = self._scan(params)
        chi_bar = _solve_chi_avg(
            atom, build_H0, (2, 0), scan, params, h_dep=True,
            probe_omega=probe, k_vec=medium["k_vec"], mass=medium["mass"])
        return dict(scan=scan, chi_bar=chi_bar, N=medium["N"], T=medium["T"], Dc=Dc,
                    L=params["cell_mm"] * 1e-3, ls=params.get("line_strength", 1.0),
                    gamma=gamma, gamma_mhz=gamma_mhz, k_vec=medium["k_vec"],
                    omega0=medium["omega0"], dipole=medium["dipole"],
                    medium_label=medium["label"], line=medium["line"],
                    coupling_rabi_mhz=Oc_mhz, buffer_ground_relax_khz=gg_khz)

    def observables(self, raw, params):
        import matplotlib.pyplot as plt
        m = self._mode(params)
        x = raw["scan"] / (2 * np.pi) / 1e6                      # MHz
        alpha, xphys = observables.absorption_coefficient(
            raw["chi_bar"], raw.get("k_vec", K_VEC), raw["N"],
            dipole=raw.get("dipole"), line_strength=raw["ls"])
        L = params["cell_mm"] * 1e-3           # navigate-only knob: read live
        T_trans = observables.transmission(alpha, L)
        center = raw["Dc"] / (2 * np.pi) / 1e6                   # two-photon resonance, MHz

        fig, (axT, axD) = plt.subplots(2, 1, figsize=(8.5, 6.4), sharex=True)
        axT.plot(x, T_trans, color="#1f77b4", lw=1.8)
        axT.set_ylabel("Transmission")
        axT.set_ylim(-0.02, 1.02)
        axT.set_title(
            f"{self._TITLE[m].split('(')[0].strip()}: "
            f"{raw.get('medium_label', '85Rb')} {raw.get('line', 'D1')}, "
            f"Omega_c = {raw['coupling_rabi_mhz']:.2f} MHz, "
            f"buffer relax = {raw['buffer_ground_relax_khz']:.1f} kHz, "
            f"Doppler {params.get('doppler', 'off')}")
        axD.plot(x, np.real(xphys), color="#9467bd", lw=1.6)
        axD.set_ylabel("Re chi  (dispersion)")
        axD.set_xlabel("Probe detuning  [MHz]")
        for a in (axT, axD):
            a.axvline(center, color="gray", ls=":", lw=0.8)
        fig.tight_layout()

        ic = int(np.argmin(np.abs(x - center)))
        metrics = self._metrics(m, x, alpha, T_trans, xphys, ic, center, raw, params)
        derived = derived_table([
            ("Medium", f"{raw.get('medium_label', '85Rb')} {raw.get('line', 'D1')}"),
            ("Γ / 2π", f"{raw.get('gamma_mhz', GAMMA_MHZ):.4f} MHz"),
            ("Ω_c / 2π", f"{raw['coupling_rabi_mhz']:.4f} MHz"),
            ("Buffer ground relax / 2π", f"{raw['buffer_ground_relax_khz']:.3f} kHz"),
            ("N", f"{raw['N']:.3e} /m³"),
        ])
        return dict(metrics=metrics, figure=fig, tables=[derived])

    def _metrics(self, mode, x, alpha, T_trans, xphys, ic, center, raw, params):
        if mode == "at":
            # Two absorption maxima straddling the (transparent) two-photon resonance.
            left = x < center
            right = x > center
            xl = x[left][np.argmax(alpha[left])] if left.any() else np.nan
            xr = x[right][np.argmax(alpha[right])] if right.any() else np.nan
            split = abs(xr - xl)
            expected = raw.get("coupling_rabi_mhz",
                               params.get("coupling_rabi", 0.0) * GAMMA_MHZ)
            return [
                dict(label="AT splitting", value=f"{split:.1f} MHz",
                     help="Separation of the two absorption maxima."),
                dict(label="Expected ≈ Ω_c", value=f"{expected:.1f} MHz",
                     help="Coupling Rabi frequency Ω_c/2π."),
                dict(label="Transmission at center", value=f"{T_trans[ic]:.3f}"),
            ]
        # EIT / CPT: transparency at two-photon resonance + dark-resonance width.
        win_fwhm = window_fwhm(x, T_trans, ic)
        ng = observables.group_index(xphys, raw["scan"], raw.get("omega0", OMEGA_D1))[ic]
        unit = "kHz" if mode == "cpt" else "MHz"
        scale = 1e3 if mode == "cpt" else 1.0
        return [
            dict(label="Transmission at resonance", value=f"{T_trans[ic]:.3f}",
                 help="Probe transmission at two-photon resonance (transparency)."),
            dict(label="Window FWHM", value=f"{win_fwhm*scale:.2f} {unit}",
                 help="Width of the transparency feature."),
            dict(label="Group index n_g", value=f"{ng:.3e}",
                 help="n_g = n + ω dn/dω at resonance (slow light when ≫ 1)."),
        ]
