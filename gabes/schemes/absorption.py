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

from .. import atoms, constants, doppler, hyperfine, observables
from ..constants import GAMMA, K_VEC, OMEGA_D1
from .. import core
from .base import ParamSpec, Preset, Scheme

PROBE_RABI = 1e-3              # weak probe, in units of Γ
GAMMA_MHZ = GAMMA / (2 * np.pi) / 1e6

# numpy ≥ 2.0 renamed np.trapz → np.trapezoid; support both.
_trapz = getattr(np, "trapezoid", getattr(np, "trapz", None))


_TABLE_STEP = GAMMA / 30.0          # Δ_eff sampling for the Doppler χ̄ table


def _solve_chi_avg(atom, build_H0, probe_coh, scan, params, h_dep):
    """
    χ̄(scan) = probe coherence ρ_eg / Ω_p, Doppler-averaged over the Maxwell
    velocity distribution (single excited-state shift, co-propagating geometry).

    Doppler off: a single v = 0 solve per scan point.
    Doppler on, scan-independent H (OD): one fine χ̄(Δ_eff) table + interpolated
        Doppler average — decouples velocity sampling from the solve (accurate
        Voigt without millions of solves), mirroring the FWM Δ_eff trick.
    Doppler on, scan-dependent H (Λ): per scan point, a local fine Δ_eff table
        around s, then interpolated average. The Λ two-photon feature is
        Doppler-free (exact); only the broad optical background is averaged.
    """
    e, g = probe_coh
    om_p = PROBE_RABI * GAMMA
    n = atom.n_levels
    doppler_on = params.get("doppler", "off") == "on"

    if not doppler_on:
        if not h_dep:
            L0 = core.build_liouvillian(build_H0(0.0), atom)
            rho = core.steady_state_batched(L0, scan, atom.S_v, n)
            return rho[:, e, g] / om_p
        out = np.zeros(scan.size, dtype=complex)
        for i, s in enumerate(scan):
            L0 = core.build_liouvillian(build_H0(s), atom)
            rho = core.steady_state_batched(L0, np.array([s]), atom.S_v, n)
            out[i] = rho[0, e, g] / om_p
        return out

    T = params["temp_c"] + 273.15
    v, w = doppler.velocity_grid(T, dv=1.0, cutoff_sigma=4.0)
    kvmax = K_VEC * np.abs(v).max()

    if not h_dep:
        L0 = core.build_liouvillian(build_H0(0.0), atom)
        lo, hi = scan.min() - kvmax, scan.max() + kvmax
        deff = np.linspace(lo, hi, int((hi - lo) / _TABLE_STEP) + 2)
        rho = core.steady_state_batched(L0, deff, atom.S_v, n)
        table = (rho[:, e, g] / om_p)[None, :]
        return np.array([doppler.doppler_average(table, deff, s, v, w)[0] for s in scan])

    out = np.zeros(scan.size, dtype=complex)
    for i, s in enumerate(scan):
        L0 = core.build_liouvillian(build_H0(s), atom)
        lo, hi = s - kvmax, s + kvmax
        deff = np.linspace(lo, hi, int((hi - lo) / _TABLE_STEP) + 2)
        rho = core.steady_state_batched(L0, deff, atom.S_v, n)
        table = (rho[:, e, g] / om_p)[None, :]
        out[i] = doppler.doppler_average(table, deff, s, v, w)[0]
    return out


def _fwhm(x, y):
    """FWHM of a single peak in y over x (same units as x); nan if ill-defined."""
    pk = np.nanmax(y)
    if not np.isfinite(pk) or pk <= 0:
        return float("nan")
    above = x[y >= 0.5 * pk]
    return float(above.max() - above.min()) if above.size > 1 else float("nan")


def _window_fwhm(x, y, ic):
    """
    Width of a transparency feature: a peak in y (transmission) at index ic
    sitting above an adjacent absorption floor. Walks out from ic to the
    half-height between the peak and the deepest neighbouring value.
    """
    peak = y[ic]
    floor = np.nanmin(y)
    if not np.isfinite(peak) or peak <= floor:
        return float("nan")
    thresh = 0.5 * (peak + floor)
    i = ic
    while i > 0 and y[i] >= thresh:
        i -= 1
    j = ic
    while j < y.size - 1 and y[j] >= thresh:
        j += 1
    return float(x[j] - x[i])


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
    gamma_eff = hyperfine.self_broadened_gamma(N)
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

    info = dict(N=N, gamma_eff=gamma_eff, sigma_v=sigma_v, dopp_fwhm=dopp_fwhm)
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
            ParamSpec("cell_mm", "Cell length", "Cell & beams", 10.0, 0.5, 200.0, 0.5, "mm"),
            ParamSpec("line_strength", "Line-strength factor", "Detection & scaling", 1.0,
                      0.01, 2.0, 0.01, "", help="Effective |d|² calibration knob. Single: "
                      "=1.0 reproduces the textbook 3λ²/2π cross-section. Hyperfine: "
                      "=1.0 reproduces the validated AutoOD absolute scale."),
            ParamSpec("doppler", "Doppler (vapor motion)", "Numerics", "on",
                      choices=("on", "off"), advanced=True),
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
        gamma_eff = hyperfine.self_broadened_gamma(N)
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
        OD = observables.optical_density(alpha, raw["L"])
        T_trans = observables.transmission(alpha, raw["L"])

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

        sb_mhz = (raw["gamma_eff"] - GAMMA) / (2 * np.pi) / 1e6
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
        derived = (
            f"| Quantity | Value |\n|---|---|\n"
            f"| N(85Rb), pure cell | {raw['N']:.3e} /m³ |\n"
            f"| σ_v (1-D) | {raw['sigma_v']:.1f} m/s |\n"
            f"| Doppler FWHM | {raw['dopp_fwhm']/(2*np.pi)/1e6:.1f} MHz |\n"
            f"| Γ_eff/2π (self-broadened) | {raw['gamma_eff']/(2*np.pi)/1e6:.3f} MHz "
            f"(+{sb_mhz:.3f}) |\n"
            f"| Natural Γ/2π | {GAMMA_MHZ:.3f} MHz |\n"
        )
        lines = ("| Line | Center [GHz] | p_F·C_F² |\n|---|---|---|\n" + rows)
        return dict(metrics=metrics, figure=fig, tables=[
            {"title": "Derived quantities", "markdown": derived},
            {"title": "Hyperfine lines (relative strength)", "markdown": lines}])

    # ---- single bare 2-level line (validation backbone) ----
    def _compute_single(self, params):
        T = params["temp_c"] + 273.15
        sigma_v = np.sqrt(constants.KB * T / constants.MASS_85RB)
        dopp_fwhm = np.sqrt(8 * np.log(2)) * K_VEC * sigma_v
        half = max(10 * GAMMA, 3.5 * dopp_fwhm) if params["doppler"] == "on" else 12 * GAMMA
        scan = np.linspace(-half, half, 601)

        atom = atoms.two_level()

        def build_H0(_s):
            H = np.zeros((2, 2), dtype=complex)
            H[0, 1] = H[1, 0] = PROBE_RABI * GAMMA / 2
            return H

        chi_bar = _solve_chi_avg(atom, build_H0, (1, 0), scan, params, h_dep=False)
        N = atoms.rb85_density(T)
        return dict(model="single", scan=scan, chi_bar=chi_bar, N=N, T=T,
                    L=params["cell_mm"] * 1e-3, ls=params["line_strength"],
                    sigma_v=sigma_v, dopp_fwhm=dopp_fwhm)

    def _observables_single(self, raw, params):
        import matplotlib.pyplot as plt
        x = raw["scan"] / (2 * np.pi) / 1e6                      # MHz
        alpha, _ = observables.absorption_coefficient(
            raw["chi_bar"], K_VEC, raw["N"], line_strength=raw["ls"])
        OD = observables.optical_density(alpha, raw["L"])
        T_trans = observables.transmission(alpha, raw["L"])

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

        fwhm_mhz = _fwhm(x, OD)
        metrics = [
            dict(label="Peak OD", value=f"{np.nanmax(OD):.3f}",
                 help="On-line optical density −log₁₀(T)."),
            dict(label="Min transmission", value=f"{np.nanmin(T_trans):.3f}"),
            dict(label="Line FWHM", value=f"{fwhm_mhz:.1f} MHz",
                 help="Voigt width (Doppler on) or natural Γ (Doppler off)."),
        ]
        derived = (
            f"| Quantity | Value |\n|---|---|\n"
            f"| N(85Rb) | {raw['N']:.3e} /m³ |\n"
            f"| σ_v (1-D) | {raw['sigma_v']:.1f} m/s |\n"
            f"| Doppler FWHM | {raw['dopp_fwhm']/(2*np.pi)/1e6:.1f} MHz |\n"
            f"| Natural Γ/2π | {GAMMA_MHZ:.3f} MHz |\n"
        )
        return dict(metrics=metrics, figure=fig,
                    tables=[{"title": "Derived quantities", "markdown": derived}])

    def compute(self, params):
        if params.get("model", "85Rb D1 hyperfine") == "single 2-level":
            return self._compute_single(params)
        return self._compute_hyperfine(params)

    def observables(self, raw, params):
        if raw.get("model") == "single":
            return self._observables_single(raw, params)
        return self._observables_hyperfine(raw, params)


# =========================================================
# Λ system — EIT / AT / CPT (one engine, three presentations)
# =========================================================
class LambdaScheme(Scheme):
    cluster = "A — Absorption"

    def __init__(self, mode):
        self.mode = mode
        self.name = mode
        if mode == "eit":
            self.title = "Electromagnetically induced transparency (EIT)"
            self.caption = ("Λ system: a strong coupling field opens a transparency "
                            "window for a weak probe at two-photon resonance.")
        elif mode == "at":
            self.title = "Autler-Townes splitting (AT)"
            self.caption = ("Strong coupling dresses the excited state; the weak-probe "
                            "absorption line splits into a doublet separated by ≈ Ω_c.")
        else:
            self.title = "Coherent population trapping (CPT)"
            self.caption = ("Narrow two-photon (Raman) scan: atoms are pumped into a "
                            "dark state, giving a sub-natural absorption dip.")

    # Per-mode defaults chosen so the feature is unsaturated and clearly visible
    # (Rb is very absorbing — ls=1.0 is the true cross-section, so cells are short).
    _DEF = {
        "eit": dict(oc=3.0, gg=0.01, temp=50.0, cell=15.0, dopp="on"),
        "at":  dict(oc=8.0, gg=0.01, temp=25.0, cell=3.0, dopp="off"),
        "cpt": dict(oc=1.0, gg=0.005, temp=25.0, cell=3.0, dopp="off"),
    }

    def param_schema(self):
        d = self._DEF[self.mode]
        return [
            ParamSpec("coupling_detuning", "Coupling detuning Δ_c", "Detunings", 0.0,
                      -10.0, 10.0, 0.1, "Γ"),
            ParamSpec("coupling_rabi", "Coupling Rabi Ω_c", "Fields", d["oc"],
                      0.1, 20.0, 0.1, "Γ", help="Strong dressing field on g₂↔e."),
            ParamSpec("gamma_gg", "Ground dephasing γ_gg", "Atomic", d["gg"],
                      0.0, 0.5, 0.001, "Γ", help="Sets the EIT/CPT dark-resonance floor."),
            ParamSpec("temp_c", "Temperature", "Cell & beams", d["temp"], 20.0, 200.0, 1.0, "°C"),
            ParamSpec("cell_mm", "Cell length", "Cell & beams", d["cell"], 0.5, 200.0, 0.5, "mm"),
            ParamSpec("line_strength", "Line-strength factor", "Detection & scaling", 1.0,
                      0.01, 2.0, 0.01, ""),
            ParamSpec("doppler", "Doppler (vapor motion)", "Numerics", d["dopp"],
                      choices=("on", "off"), advanced=True),
        ]

    def presets(self):
        if self.mode == "eit":
            return [Preset("Vapor EIT", values=dict(coupling_rabi=3.0, gamma_gg=0.01,
                                                    temp_c=50.0, doppler="on"))]
        if self.mode == "at":
            return [Preset("Strong-coupling doublet",
                           values=dict(coupling_rabi=8.0, gamma_gg=0.01, doppler="off"))]
        return [Preset("Dark resonance",
                       values=dict(coupling_rabi=1.0, gamma_gg=0.005, doppler="off"))]

    def _scan(self, params):
        Oc = params["coupling_rabi"]
        Dc = params["coupling_detuning"]
        gg = params["gamma_gg"]
        if self.mode == "cpt":
            dark = gg + Oc**2          # power-broadened dark width (Γ units)
            half = max(0.05, 12.0 * dark)
        elif self.mode == "at":
            half = max(10.0, 3.0 * Oc)
        else:
            half = max(8.0, 4.0 * Oc)
        return (Dc + np.linspace(-half, half, 601)) * GAMMA

    def compute(self, params):
        atom = atoms.lambda3(gamma_gg=params["gamma_gg"] * GAMMA)
        Oc = params["coupling_rabi"] * GAMMA
        Dc = params["coupling_detuning"] * GAMMA

        def build_H0(s):
            H = np.zeros((3, 3), dtype=complex)
            H[1, 1] = Dc - s                      # two-photon detuning (Doppler-free)
            H[0, 2] = H[2, 0] = PROBE_RABI * GAMMA / 2
            H[1, 2] = H[2, 1] = Oc / 2
            return H

        scan = self._scan(params)
        chi_bar = _solve_chi_avg(atom, build_H0, (2, 0), scan, params, h_dep=True)
        T = params["temp_c"] + 273.15
        N = atoms.rb85_density(T)
        return dict(scan=scan, chi_bar=chi_bar, N=N, T=T, Dc=Dc,
                    L=params["cell_mm"] * 1e-3, ls=params["line_strength"])

    def observables(self, raw, params):
        import matplotlib.pyplot as plt
        x = raw["scan"] / (2 * np.pi) / 1e6                      # MHz
        alpha, xphys = observables.absorption_coefficient(
            raw["chi_bar"], K_VEC, raw["N"], line_strength=raw["ls"])
        T_trans = observables.transmission(alpha, raw["L"])
        center = raw["Dc"] / (2 * np.pi) / 1e6                   # two-photon resonance, MHz

        fig, (axT, axD) = plt.subplots(2, 1, figsize=(8.5, 6.4), sharex=True)
        axT.plot(x, T_trans, color="#1f77b4", lw=1.8)
        axT.set_ylabel("Transmission")
        axT.set_ylim(-0.02, 1.02)
        axT.set_title(f"{self.title.split('(')[0].strip()}:  Ω_c = {params['coupling_rabi']:.1f} Γ,  "
                      f"γ_gg = {params['gamma_gg']:.3f} Γ,  Doppler {params['doppler']}")
        axD.plot(x, np.real(xphys), color="#9467bd", lw=1.6)
        axD.set_ylabel("Re χ  (dispersion)")
        axD.set_xlabel("Probe detuning  [MHz]")
        for a in (axT, axD):
            a.axvline(center, color="gray", ls=":", lw=0.8)
        fig.tight_layout()

        ic = int(np.argmin(np.abs(x - center)))
        metrics = self._metrics(x, alpha, T_trans, xphys, ic, center, raw, params)
        return dict(metrics=metrics, figure=fig, tables=[])

    def _metrics(self, x, alpha, T_trans, xphys, ic, center, raw, params):
        if self.mode == "at":
            # Two absorption maxima straddling the (transparent) two-photon resonance.
            left = x < center
            right = x > center
            xl = x[left][np.argmax(alpha[left])] if left.any() else np.nan
            xr = x[right][np.argmax(alpha[right])] if right.any() else np.nan
            split = abs(xr - xl)
            expected = params["coupling_rabi"] * GAMMA_MHZ      # Ω_c in MHz
            return [
                dict(label="AT splitting", value=f"{split:.1f} MHz",
                     help="Separation of the two absorption maxima."),
                dict(label="Expected ≈ Ω_c", value=f"{expected:.1f} MHz",
                     help="Coupling Rabi frequency Ω_c/2π."),
                dict(label="Transmission at center", value=f"{T_trans[ic]:.3f}"),
            ]
        # EIT / CPT: transparency at two-photon resonance + dark-resonance width.
        win_fwhm = _window_fwhm(x, T_trans, ic)
        ng = observables.group_index(xphys, raw["scan"], OMEGA_D1)[ic]
        unit = "kHz" if self.mode == "cpt" else "MHz"
        scale = 1e3 if self.mode == "cpt" else 1.0
        return [
            dict(label="Transmission at resonance", value=f"{T_trans[ic]:.3f}",
                 help="Probe transmission at two-photon resonance (transparency)."),
            dict(label="Window FWHM", value=f"{win_fwhm*scale:.2f} {unit}",
                 help="Width of the transparency feature."),
            dict(label="Group index n_g", value=f"{ng:.3e}",
                 help="n_g = n + ω dn/dω at resonance (slow light when ≫ 1)."),
        ]
