"""
Cluster C - 87Rb D1 magneto-optics (Hanle / EIA / NMOR).

The practical default is the 87Rb D1 F_g=1 -> F_e=2 open transition, a
laboratory Hanle-configuration EIA peak.  The solver keeps the compact Zeeman
OBE engine, but the axes and atomic constants are now experimental:

  * magnetic-field scan in microtesla, converted with Steck low-field g_F values
  * 87Rb D1 frequency, linewidth, mass, dipole strength, and vapor density
  * beam intensity in mW/cm^2, mapped through Steck's D1 saturation intensity
  * optional Doppler averaging and transit/repopulation relaxation

This is still a single-addressed-F manifold model, not a full two-ground-state
open-transition propagation code.  It is designed to give a realistic, useful
zero-field peak/dip shape with physically named knobs.
"""
import math

import numpy as np

from .. import constants, core, observables, species, zeeman
from .base import ParamSpec, Preset, Scheme

MAX_LEVELS = 16

RB87 = species.RB87
LINE = "D1"
JE_D1, NU_D1_RB87, GAMMA_D1_MHZ, _, _ = RB87.line(LINE)
GAMMA_D1 = GAMMA_D1_MHZ * species.MHZ
K_D1_RB87 = 2 * np.pi * NU_D1_RB87 / constants.C_LIGHT

# Steck, Rubidium 87 D Line Data, Table 7: D1 far-detuned effective saturation
# intensity for pi-polarized light.  The Hanle geometry uses transverse linear
# light (sigma+ + sigma-), so this is an intensity scale rather than a claim that
# every Zeeman subtransition has the same saturation intensity.
D1_ISAT_MW_CM2 = 4.4876

GJ_5S12 = 2.002_331_070
GJ_5P12 = 0.666
GI_RB87 = -0.000_995_141_4


def _hyperfine_gf(I, J, F, gJ):
    """Low-field hyperfine Lande g_F in Bohr-magneton units."""
    F = float(F)
    if F == 0:
        return 0.0
    ff = F * (F + 1)
    electronic = gJ * (ff + J * (J + 1) - I * (I + 1)) / (2 * ff)
    nuclear = GI_RB87 * (ff + I * (I + 1) - J * (J + 1)) / (2 * ff)
    return electronic + nuclear


def _thermal_velocity_grid(T, mass, n_classes, cutoff_sigma=3.0):
    """Small Maxwell-Boltzmann quadrature grid for interactive Doppler averaging."""
    n_classes = int(max(1, n_classes))
    if n_classes <= 1:
        return np.array([0.0]), np.array([1.0])
    sigma_v = math.sqrt(constants.KB * T / mass)
    v = np.linspace(-cutoff_sigma * sigma_v, cutoff_sigma * sigma_v, n_classes)
    wt = np.exp(-0.5 * (v / sigma_v) ** 2)
    wt /= wt.sum()
    return v, wt


def _doppler_dilution(T, mass, k_vec, gamma, detuning, v, wt):
    """
    Correct the optical-depth scale for an intentionally coarse velocity grid.

    The Zeeman solve can use ~20 velocity classes interactively, but the optical
    Lorentzian is only a few m/s wide in velocity units.  Without this scalar
    correction the v=0 class over-represents resonant atoms.  We match the coarse
    scalar Lorentzian average to a fine-grid Voigt average at the selected laser
    detuning.
    """
    if v.size <= 1:
        return 1.0
    hwhm = gamma / 2.0

    def profile(delta):
        return hwhm * hwhm / (delta * delta + hwhm * hwhm)

    coarse = float(np.sum(wt * profile(detuning - k_vec * v)))
    if coarse <= 0:
        return 1.0

    sigma_v = math.sqrt(constants.KB * T / mass)
    dv = max(0.5, hwhm / k_vec / 4.0)
    n_fine = int(np.clip(8 * sigma_v / dv, 801, 4001))
    vf = np.linspace(-4 * sigma_v, 4 * sigma_v, n_fine)
    wf = np.exp(-0.5 * (vf / sigma_v) ** 2)
    wf /= wf.sum()
    exact = float(np.sum(wf * profile(detuning - k_vec * vf)))
    return exact / coarse if exact > 0 else 1.0


def _transition_strength(Fg, Fe):
    return species.line_strength(Fg, Fe, RB87.I, RB87.Jg, JE_D1)


def _transition_dipole(Fg, Fe):
    lam = constants.C_LIGHT / NU_D1_RB87
    d2_j = species.reduced_dipole_sq(GAMMA_D1, lam, RB87.Jg, JE_D1)
    return math.sqrt(max(_transition_strength(Fg, Fe), 0.0) * d2_j)


class MagnetoScheme(Scheme):
    cluster = "C - Magneto-optics"
    presets_group = "Default"
    cache_version = "2"
    defaults_version = "2"

    _DEF = {
        "hanle": dict(view="Hanle", Fg=2, Fe=1, intensity=0.8, gg_khz=20.0,
                      bmax=80.0, title="Hanle effect",
                      desc="87Rb D1 dark Hanle resonance on F=2 -> F'=1."),
        "eia": dict(view="EIA", Fg=1, Fe=2, intensity=2.0, gg_khz=30.0,
                    bmax=120.0, title="Electromagnetically induced absorption (EIA)",
                    desc="87Rb D1 bright Hanle/EIA peak on the open F=1 -> F'=2 transition."),
        "nmor": dict(view="NMOR", Fg=2, Fe=1, intensity=0.7, gg_khz=10.0,
                     bmax=80.0, title="Nonlinear magneto-optical rotation (NMOR)",
                     desc="87Rb D1 polarization-rotation readout near zero magnetic field."),
    }

    def __init__(self, mode=None):
        self.mode = mode
        if mode is None:
            self.name = "magneto"
            self.title = "Magneto-optics (Hanle/EIA/NMOR)"
            self.caption = ("Practical 87Rb D1 Hanle-configuration magneto-optics. "
                            "The default button simulates the D1 F=1 -> F'=2 EIA "
                            "absorption peak with a real B-field axis, Steck atomic "
                            "data, Doppler averaging, and transit relaxation.")
        else:
            self.name = mode
            self.title = self._DEF[mode]["title"]
            self.caption = self._DEF[mode]["desc"]

    def _mode(self, params):
        """Effective readout: pinned alias instance or the merged hidden `view`."""
        return self.mode or str(params.get("view", "EIA")).lower()

    def param_schema(self):
        d = self._DEF[self.mode or "eia"]
        specs = []
        if self.mode is None:
            specs.append(ParamSpec(
                "view", "Readout", "Readout", d["view"],
                choices=("Hanle", "EIA", "NMOR"), recompute=False,
                help="Internal readout selected by the default buttons.",
                hidden=True))
        specs += [
            ParamSpec("Fg", "Ground F_g", "Atomic", float(d["Fg"]), 1.0, 2.0, 1.0, "",
                      help="87Rb D1 ground hyperfine level."),
            ParamSpec("Fe", "Excited F'_e", "Atomic", float(d["Fe"]), 1.0, 2.0, 1.0, "",
                      help="87Rb D1 excited hyperfine level. |F'_e-F_g| <= 1."),
            ParamSpec("intensity_mw_cm2", "Beam intensity", "Fields", d["intensity"],
                      0.01, 20.0, 0.05, "mW/cm²",
                      help="Single linearly-polarized beam intensity. Converted with "
                           f"Steck D1 I_sat = {D1_ISAT_MW_CM2:.4f} mW/cm²."),
            ParamSpec("b_max_ut", "B scan ±", "Fields", d["bmax"],
                      1.0, 500.0, 1.0, "µT",
                      help="Longitudinal magnetic-field scan range around zero."),
            ParamSpec("laser_detuning_mhz", "Laser detuning", "Detunings", 0.0,
                      -1000.0, 1000.0, 5.0, "MHz",
                      help="Detuning from the selected 87Rb D1 hyperfine transition."),
            ParamSpec("ground_relax_khz", "Ground relaxation γg/2π", "Atomic", d["gg_khz"],
                      0.1, 1000.0, 1.0, "kHz",
                      help="Transit / wall / spin relaxation. Sets the zero-field width "
                           "and repopulates the addressed ground manifold."),
            ParamSpec("temp_c", "Temperature", "Cell & beams", 25.0,
                      20.0, 120.0, 1.0, "°C",
                      help="Sets 87Rb vapor density and Doppler width."),
            ParamSpec("cell_mm", "Cell length", "Cell & beams", 10.0,
                      0.5, 200.0, 0.5, "mm", recompute=False),
            ParamSpec("line_strength", "Line-strength factor", "Detection & scaling", 1.0,
                      0.01, 2.0, 0.01, "", recompute=False,
                      help="Experimental calibration knob for effective optical depth."),
            ParamSpec("doppler", "Doppler averaging", "Numerics", "on",
                      choices=("on", "off"), advanced=True,
                      help="Coarse Maxwell-Boltzmann average. Keep on for the practical "
                           "87Rb D1 peak; off is faster and closer to the textbook toy."),
            ParamSpec("velocity_classes", "Velocity classes", "Numerics", 21,
                      1, 81, 2, "", advanced=True),
            ParamSpec("scan_points", "B scan points", "Numerics", 161,
                      51, 401, 10, "", advanced=True),
        ]
        return specs

    def presets(self):
        if self.mode is None:
            return [
                Preset("D1 Hanle dip", icon="🔻", values=dict(
                    view="Hanle", Fg=2.0, Fe=1.0, intensity_mw_cm2=0.8,
                    ground_relax_khz=20.0, b_max_ut=80.0, temp_c=25.0,
                    cell_mm=10.0, doppler="on")),
                Preset("D1 EIA peak", icon="🔺", values=dict(
                    view="EIA", Fg=1.0, Fe=2.0, intensity_mw_cm2=2.0,
                    ground_relax_khz=30.0, b_max_ut=120.0, temp_c=25.0,
                    cell_mm=10.0, doppler="on")),
                Preset("D1 NMOR rotation", icon="🧭", values=dict(
                    view="NMOR", Fg=2.0, Fe=1.0, intensity_mw_cm2=0.7,
                    ground_relax_khz=10.0, b_max_ut=80.0, temp_c=25.0,
                    cell_mm=10.0, doppler="on")),
            ]
        d = self._DEF[self.mode]
        icon = {"hanle": "🔻", "eia": "🔺", "nmor": "🧭"}[self.mode]
        return [Preset(f"D1 {self.mode.upper()} default", icon=icon,
                       values=dict(Fg=float(d["Fg"]), Fe=float(d["Fe"]),
                                   intensity_mw_cm2=d["intensity"],
                                   ground_relax_khz=d["gg_khz"],
                                   b_max_ut=d["bmax"], temp_c=25.0,
                                   cell_mm=10.0, doppler="on"))]

    def info(self):
        return (
            "**87Rb D1 magneto-optics.** The practical EIA preset is the "
            "Hanle-configuration D1 `F_g=1 -> F'_e=2` bright resonance.  The "
            "Hamiltonian is still a compact single-addressed-F Zeeman manifold, "
            "with transit/repopulation relaxation standing in for atoms entering "
            "and leaving the laser beam.\n\n"
            "Atomic constants are taken from Steck's *Rubidium 87 D Line Data*: "
            "D1 centroid, natural linewidth, hyperfine strengths, dipole moment, "
            "D1 saturation-intensity scale, and low-field Lande g-factors.  The "
            "87Rb vapor density uses the same Steck/CRC vapor-pressure fit as "
            "the OD/SAS data layer.\n\n"
            "**References**\n"
            "- D. A. Steck, *Rubidium 87 D Line Data*, http://steck.us/alkalidata.\n"
            "- A. S. Zibrov and A. B. Matsko, induced absorption on the open "
            "87Rb D1 `F_g=1 -> F_e=2` transition, arXiv:physics/0512199.\n"
            "- D. V. Brazhnikov et al., high-contrast 87Rb D1 Hanle EIA, "
            "*Laser Physics Letters* 11, 125702 (2014)."
        )

    def compute(self, params):
        Fg, Fe = int(round(params["Fg"])), int(round(params["Fe"]))
        Fg_allowed = set(int(round(f)) for f in species.f_values(RB87.I, RB87.Jg))
        Fe_allowed = set(int(round(f)) for f in species.f_values(RB87.I, JE_D1))
        strength = _transition_strength(Fg, Fe) if (Fg in Fg_allowed and Fe in Fe_allowed) else 0.0
        valid = (Fg in Fg_allowed and Fe in Fe_allowed and abs(Fe - Fg) <= 1
                 and strength > 1e-12 and (2 * Fg + 1) + (2 * Fe + 1) <= MAX_LEVELS)

        gFg = _hyperfine_gf(RB87.I, RB87.Jg, Fg, GJ_5S12) if Fg in Fg_allowed else 0.0
        gFe = _hyperfine_gf(RB87.I, JE_D1, Fe, GJ_5P12) if Fe in Fe_allowed else 0.0
        gamma_g = 2 * np.pi * params["ground_relax_khz"] * 1e3
        g_ratio = gFe / gFg if abs(gFg) > 1e-12 else 0.0
        atom = (zeeman.zeeman_manifold(Fg, Fe, gamma=GAMMA_D1, gamma_gg=gamma_g,
                                       g_ratio=g_ratio, transit_rate=gamma_g)
                if valid else None)

        b_ut = np.linspace(-params["b_max_ut"], params["b_max_ut"],
                           int(params["scan_points"]))
        b_t = b_ut * 1e-6
        larmor = constants.MU_B * gFg * b_t / constants.HBAR

        T = params["temp_c"] + 273.15
        N = species.number_density(RB87, T)
        intensity = max(float(params["intensity_mw_cm2"]), 0.0)
        Om = GAMMA_D1 * math.sqrt(intensity / (2 * D1_ISAT_MW_CM2)) if intensity > 0 else 0.0
        dL = 2 * np.pi * params["laser_detuning_mhz"] * 1e6

        chi_p = np.zeros(b_ut.size, dtype=complex)
        chi_m = np.zeros(b_ut.size, dtype=complex)
        velocity_count = 1
        doppler_scale = 1.0
        if valid and Om > 0:
            if params["doppler"] == "on":
                v, wt = _thermal_velocity_grid(T, RB87.mass, params["velocity_classes"])
                doppler_scale = _doppler_dilution(T, RB87.mass, K_D1_RB87, GAMMA_D1, dL, v, wt)
            else:
                v, wt = np.array([0.0]), np.array([1.0])
            velocity_count = int(v.size)
            deff = dL - K_D1_RB87 * v
            cpl_p, cpl_m = atom.couplings[+1], atom.couplings[-1]
            for j, OmL in enumerate(larmor):
                H = self._hamiltonian(atom, OmL, Om)
                L0 = core.build_liouvillian(H, atom)
                rho = core.steady_state_batched(L0, deff, atom.S_v, atom.n_levels)
                cp = sum(cg * rho[:, ei, gi] for gi, ei, cg in cpl_p) / Om
                cm = sum(cg * rho[:, ei, gi] for gi, ei, cg in cpl_m) / Om
                chi_p[j] = (cp * wt).sum()
                chi_m[j] = (cm * wt).sum()

        return dict(
            line=LINE, isotope="87Rb", b_ut=b_ut, larmor=larmor,
            chi_p=chi_p, chi_m=chi_m, N=N, N_eff=N * doppler_scale, T=T,
            L=params["cell_mm"] * 1e-3, ls=params["line_strength"],
            valid=valid, Fg=Fg, Fe=Fe, gFg=gFg, gFe=gFe,
            gamma=GAMMA_D1, k_vec=K_D1_RB87, dipole=_transition_dipole(Fg, Fe) if valid else 0.0,
            strength=strength, omega_rabi=Om, isat=D1_ISAT_MW_CM2,
            velocity_count=velocity_count, doppler_scale=doppler_scale,
        )

    @staticmethod
    def _hamiltonian(atom, OmL, Om):
        n = atom.n_levels
        H = np.zeros((n, n), dtype=complex)
        for i in atom.ground:
            H[i, i] = OmL * atom.m_ground[i]
        for k, e in enumerate(atom.excited):
            H[e, e] = atom.g_ratio * OmL * atom.m_excited[k]
        for q in (+1, -1):
            for gi, ei, cg in atom.couplings[q]:
                H[gi, ei] += Om * cg / 2
                H[ei, gi] += Om * cg / 2
        return H

    def observables(self, raw, params):
        import matplotlib.pyplot as plt

        m = self._mode(params)
        if not raw["valid"]:
            fig, ax = plt.subplots(figsize=(8.5, 3.0))
            ax.text(0.5, 0.5,
                    f"87Rb D1 F_g={raw['Fg']}, F'_e={raw['Fe']} is not an allowed "
                    f"addressed transition for this compact Zeeman model.",
                    ha="center", va="center", wrap=True)
            ax.axis("off")
            return dict(metrics=[dict(label="Status", value="invalid transition")],
                        figure=fig, tables=[])

        x = raw["b_ut"]
        xphys_p = observables.chi_phys(
            raw["chi_p"], raw["N_eff"], dipole=raw["dipole"], line_strength=raw["ls"])
        xphys_m = observables.chi_phys(
            raw["chi_m"], raw["N_eff"], dipole=raw["dipole"], line_strength=raw["ls"])
        k, L = raw["k_vec"], raw["L"]
        alpha = k * np.imag(xphys_p + xphys_m)
        T_trans = np.exp(-alpha * L)
        OD = alpha * L / np.log(10.0)
        rotation = 0.25 * k * L * np.real(xphys_p - xphys_m)
        ic = int(np.argmin(np.abs(x)))

        title = (f"87Rb D1 {m.upper()}  F={raw['Fg']}→F'={raw['Fe']},  "
                 f"I={params['intensity_mw_cm2']:.2f} mW/cm²,  "
                 f"γg/2π={params['ground_relax_khz']:.0f} kHz")
        if m == "nmor":
            fig, ax = plt.subplots(figsize=(8.5, 4.6))
            ax.plot(x, rotation * 1e3, color="#9467bd", lw=1.8)
            ax.axhline(0, color="black", lw=0.6)
            ax.axvline(0, color="gray", ls=":", lw=0.8)
            ax.set_ylabel("Polarization rotation  [mrad]")
            ax.set_xlabel("Magnetic field B  [µT]")
            ax.set_title(title)
            fig.tight_layout()
            slope = np.gradient(rotation, x)[ic]
            metrics = [
                dict(label="Rotation at B=0", value=f"{rotation[ic]*1e3:.2f} mrad"),
                dict(label="Slope dθ/dB", value=f"{slope*1e3:.2f} mrad/µT"),
                dict(label="Peak |rotation|", value=f"{np.max(np.abs(rotation))*1e3:.2f} mrad"),
            ]
            note = "NMOR readout: zero crossing near B=0; slope is the magnetometer signal."
        else:
            fig, (axT, axA) = plt.subplots(2, 1, figsize=(8.5, 6.4), sharex=True)
            axT.plot(x, T_trans, color="#1f77b4", lw=1.8)
            axT.axvline(0, color="gray", ls=":", lw=0.8)
            axT.set_ylabel("Transmission")
            axT.set_title(title)
            axA.plot(x, OD, color="#d62728", lw=1.8)
            axA.axvline(0, color="gray", ls=":", lw=0.8)
            axA.set_ylabel("Optical density")
            axA.set_xlabel("Magnetic field B  [µT]")
            fig.tight_layout()
            bg = 0.5 * (alpha[0] + alpha[-1])
            contrast = (alpha[ic] - bg) / abs(bg) if bg != 0 else 0.0
            kind = "dip (transparency)" if alpha[ic] < bg else "peak (enhanced)"
            metrics = [
                dict(label="OD at B=0", value=f"{OD[ic]:.3f}"),
                dict(label="Zero-field feature", value=kind),
                dict(label="Contrast vs edge", value=f"{contrast*100:+.1f} %"),
            ]
            note = ("Hanle readout: a dark resonance appears as reduced absorption at B=0."
                    if m == "hanle"
                    else "EIA readout: the practical D1 preset targets an enhanced absorption peak at B=0.")

        larmor_hz_per_g = constants.MU_B * abs(raw["gFg"]) / (2 * np.pi * constants.HBAR) * 1e-4
        derived = (
            "| Quantity | Value |\n|---|---|\n"
            f"| Isotope / line | {raw['isotope']} {raw['line']} |\n"
            f"| Transition | F={raw['Fg']} -> F'={raw['Fe']} |\n"
            f"| Relative strength S_FF' | {raw['strength']:.4f} |\n"
            f"| g_F ground / excited | {raw['gFg']:.4f} / {raw['gFe']:.4f} |\n"
            f"| Ground Larmor scale | {larmor_hz_per_g/1e6:.3f} MHz/G |\n"
            f"| Ω/Γ | {raw['omega_rabi']/raw['gamma']:.3f} |\n"
            f"| Γ/2π | {raw['gamma']/(2*np.pi)/1e6:.4f} MHz |\n"
            f"| I_sat(D1 scale) | {raw['isat']:.4f} mW/cm² |\n"
            f"| N(87Rb) | {raw['N']:.3e} /m³ |\n"
            f"| Doppler OD scale | {raw['doppler_scale']:.3e} |\n"
            f"| Doppler classes | {raw['velocity_count']} |\n"
        )
        return dict(metrics=metrics, figure=fig, tables=[
            {"title": "Notes", "markdown": note},
            {"title": "Derived quantities", "markdown": derived},
        ])
