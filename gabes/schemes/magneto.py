"""
Cluster C — Zeeman magneto-optics (Hanle / EIA / NMOR).

One linearly-polarised beam (= σ⁺ + σ⁻ about the B axis) on an F_g ↔ F_e
transition; the magnetic field (Larmor shift Ω_L ∝ B) is the scan axis.

  Hanle : F_e ≤ F_g → ground dark state at B = 0 → reduced absorption
          (transmission peak), a narrow zero-field resonance.
  EIA   : F_e = F_g + 1 → coherence transfer flips the sign → enhanced
          absorption at B = 0.
  NMOR  : aligned ground state makes the medium circularly birefringent →
          the linear polarization rotates; φ ∝ Re(χ₊ − χ₋), a steep zero
          crossing at B = 0 (atomic-magnetometer signal).

Engine: build H(Ω_L) = Zeeman (ground + excited) + σ⁺/σ⁻ drive, solve the
self-consistent steady state per velocity class (optical detuning / Doppler via
the excited shift), and read the σ⁺ and σ⁻ coherences →
    χ̄_± = Σ_{q=±1} CG · ρ_eg / Ω.
Absorption ∝ Im(χ₊ + χ₋); rotation ∝ Re(χ₊ − χ₋). The B-resonance lives in the
ground Zeeman coherence (Doppler-free); its width is set by γ_gg + optical pumping.
"""
import numpy as np

from .. import atoms, constants, doppler, observables, zeeman
from ..constants import GAMMA, K_VEC
from .. import core
from .base import ParamSpec, Preset, Scheme

MAX_LEVELS = 16          # dimension guard (≤ 256-dim ρ)


class MagnetoScheme(Scheme):
    cluster = "C — Magneto-optics"

    # Transitions with all ground sublevels coupled by σ± (linear ⊥ B), so the
    # dark state is B-sensitive. F=2→1 (Fe=Fg−1) → Hanle/NMOR dip; F=1→2
    # (Fe=Fg+1) → EIA peak. (F→0 traps population in the uncoupled m=0 ground
    # state — no signal — so it is avoided as a default.)
    _DEF = {
        "hanle": dict(Fg=2, Fe=1, rabi=0.5, gg=0.005, title="Hanle effect",
                      desc="Ground-state dark resonance: a narrow transmission peak at "
                           "zero field (F_e ≤ F_g)."),
        "eia":   dict(Fg=1, Fe=2, rabi=0.5, gg=0.005, title="Electromagnetically induced absorption (EIA)",
                      desc="On F_e = F_g + 1 the zero-field resonance flips sign — enhanced "
                           "absorption at B = 0."),
        "nmor":  dict(Fg=2, Fe=1, rabi=0.5, gg=0.002, title="Nonlinear magneto-optical rotation (NMOR)",
                      desc="Polarization rotation φ ∝ Re(χ₊ − χ₋): a steep zero crossing at "
                           "B = 0 — the atomic-magnetometer signal."),
    }

    def __init__(self, mode):
        self.mode = mode
        self.name = mode
        self.title = self._DEF[mode]["title"]
        self.caption = self._DEF[mode]["desc"]

    def param_schema(self):
        d = self._DEF[self.mode]
        return [
            ParamSpec("Fg", "Ground F_g", "Atomic", float(d["Fg"]), 1.0, 3.0, 1.0, ""),
            ParamSpec("Fe", "Excited F_e", "Atomic", float(d["Fe"]), 0.0, 4.0, 1.0, "",
                      help="Dipole-allowed needs |F_e − F_g| ≤ 1."),
            ParamSpec("rabi", "Beam Rabi Ω", "Fields", d["rabi"], 0.05, 5.0, 0.05, "Γ",
                      help="Optical pumping strength (power-broadens the B resonance)."),
            ParamSpec("gamma_gg", "Ground relaxation γ_gg", "Atomic", d["gg"],
                      0.0, 0.2, 0.001, "Γ", help="Transit / collisional ground-coherence decay; sets the B-resonance width."),
            ParamSpec("laser_detuning", "Laser detuning δ_L", "Detunings", 0.0,
                      -10.0, 10.0, 0.1, "Γ"),
            ParamSpec("larmor_max", "Larmor scan ±", "Detunings", 1.0, 0.05, 10.0, 0.05, "Γ",
                      help="B-field range as the ground Larmor shift Ω_L (∝ B)."),
            ParamSpec("temp_c", "Temperature", "Cell & beams", 40.0, 20.0, 200.0, 1.0, "°C"),
            ParamSpec("cell_mm", "Cell length", "Cell & beams", 2.0, 0.5, 200.0, 0.5, "mm"),
            ParamSpec("line_strength", "Line-strength factor", "Detection & scaling", 1.0,
                      0.01, 2.0, 0.01, ""),
            ParamSpec("g_ratio", "g_Fe / g_Fg", "Atomic", 1.0, -2.0, 2.0, 0.1, "", advanced=True),
            ParamSpec("doppler", "Doppler (vapor motion)", "Numerics", "off",
                      choices=("on", "off"), advanced=True),
            ParamSpec("scan_points", "Scan points", "Numerics", 201, 51, 801, 50, "", advanced=True),
        ]

    def presets(self):
        d = self._DEF[self.mode]
        return [Preset(f"{self.mode.upper()} default",
                       values=dict(Fg=float(d["Fg"]), Fe=float(d["Fe"]),
                                   rabi=d["rabi"], gamma_gg=d["gg"]), icon="🧲")]

    def compute(self, params):
        Fg, Fe = int(round(params["Fg"])), int(round(params["Fe"]))
        valid = abs(Fe - Fg) <= 1 and (2 * Fg + 1) + (2 * Fe + 1) <= MAX_LEVELS
        atom = zeeman.zeeman_manifold(Fg, Fe, gamma_gg=params["gamma_gg"] * GAMMA,
                                      g_ratio=params["g_ratio"]) if valid else None

        Lmax = params["larmor_max"] * GAMMA
        larmor = np.linspace(-Lmax, Lmax, int(params["scan_points"]))
        T = params["temp_c"] + 273.15
        N = atoms.rb85_density(T)

        chi_p = np.zeros(larmor.size, dtype=complex)
        chi_m = np.zeros(larmor.size, dtype=complex)
        if valid:
            Om = params["rabi"] * GAMMA
            dL = params["laser_detuning"] * GAMMA
            if params["doppler"] == "on":
                v, wt = doppler.velocity_grid(T, dv=4.0, cutoff_sigma=3.0)
            else:
                v, wt = np.array([0.0]), np.array([1.0])
            deff = dL - K_VEC * v
            ng = len(atom.ground)
            cpl_p, cpl_m = atom.couplings[+1], atom.couplings[-1]
            for j, OmL in enumerate(larmor):
                H = self._hamiltonian(atom, OmL, Om)
                L0 = core.build_liouvillian(H, atom)
                rho = core.steady_state_batched(L0, deff, atom.S_v, atom.n_levels)
                cp = sum(cg * rho[:, ei, gi] for gi, ei, cg in cpl_p) / Om
                cm = sum(cg * rho[:, ei, gi] for gi, ei, cg in cpl_m) / Om
                chi_p[j] = (cp * wt).sum()
                chi_m[j] = (cm * wt).sum()

        return dict(larmor=larmor, chi_p=chi_p, chi_m=chi_m, N=N,
                    L=params["cell_mm"] * 1e-3, ls=params["line_strength"],
                    valid=valid, Fg=Fg, Fe=Fe)

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
        if not raw["valid"]:
            fig, ax = plt.subplots(figsize=(8.5, 3.0))
            ax.text(0.5, 0.5, f"F_g={raw['Fg']}, F_e={raw['Fe']} not dipole-allowed "
                    f"(need |F_e−F_g| ≤ 1) or too large (≤ {MAX_LEVELS} levels).",
                    ha="center", va="center", wrap=True)
            ax.axis("off")
            return dict(metrics=[dict(label="Status", value="invalid transition")],
                        figure=fig, tables=[])

        x = raw["larmor"] / GAMMA                                  # Ω_L / Γ  (∝ B)
        xphys_p = observables.chi_phys(raw["chi_p"], raw["N"], line_strength=raw["ls"])
        xphys_m = observables.chi_phys(raw["chi_m"], raw["N"], line_strength=raw["ls"])
        k, L = K_VEC, raw["L"]
        alpha = k * np.imag(xphys_p + xphys_m)
        T_trans = np.exp(-alpha * L)
        rotation = 0.25 * k * L * np.real(xphys_p - xphys_m)       # Faraday angle [rad]
        ic = int(np.argmin(np.abs(x)))

        if self.mode == "nmor":
            fig, ax = plt.subplots(figsize=(8.5, 4.6))
            ax.plot(x, rotation * 1e3, color="#9467bd", lw=1.8)
            ax.axhline(0, color="black", lw=0.6)
            ax.axvline(0, color="gray", ls=":", lw=0.8)
            ax.set_ylabel("Polarization rotation  [mrad]")
            ax.set_xlabel("Larmor shift  Ω_L / Γ   (∝ B)")
            ax.set_title(f"NMOR  F={raw['Fg']}→{raw['Fe']},  Ω={params['rabi']:.2f} Γ,  "
                         f"γ_gg={params['gamma_gg']:.3f} Γ")
            fig.tight_layout()
            slope = np.gradient(rotation, x)[ic]                   # dφ/d(Ω_L/Γ)
            metrics = [
                dict(label="Rotation at B=0", value=f"{rotation[ic]*1e3:.2f} mrad"),
                dict(label="Sensitivity dφ/dB", value=f"{slope*1e3:.1f} mrad/Γ",
                     help="Slope of the zero crossing — magnetometer responsivity."),
                dict(label="Peak |rotation|", value=f"{np.max(np.abs(rotation))*1e3:.2f} mrad"),
            ]
            note = "Zero crossing at B=0; the steep slope is the magnetometer signal."
        else:
            fig, (axT, axA) = plt.subplots(2, 1, figsize=(8.5, 6.4), sharex=True)
            axT.plot(x, T_trans, color="#1f77b4", lw=1.8)
            axT.axvline(0, color="gray", ls=":", lw=0.8)
            axT.set_ylabel("Transmission")
            axT.set_title(f"{self.mode.upper()}  F={raw['Fg']}→{raw['Fe']},  Ω={params['rabi']:.2f} Γ,  "
                          f"γ_gg={params['gamma_gg']:.3f} Γ")
            axA.plot(x, alpha * L / np.log(10), color="#d62728", lw=1.8)
            axA.axvline(0, color="gray", ls=":", lw=0.8)
            axA.set_ylabel("Optical density")
            axA.set_xlabel("Larmor shift  Ω_L / Γ   (∝ B)")
            fig.tight_layout()
            bg = alpha[0]                                          # far-field baseline
            contrast = (alpha[ic] - bg) / abs(bg) if bg != 0 else 0.0
            kind = "dip (transparency)" if alpha[ic] < bg else "peak (enhanced)"
            metrics = [
                dict(label="Absorption at B=0", value=f"{alpha[ic]*L/np.log(10):.3f} OD"),
                dict(label="Zero-field feature", value=kind),
                dict(label="Contrast vs baseline", value=f"{contrast*100:+.0f} %"),
            ]
            note = ("Hanle: dark-state transparency dip at B=0." if self.mode == "hanle"
                    else "EIA: enhanced-absorption peak at B=0 (sign flip vs Hanle).")
        return dict(metrics=metrics, figure=fig, tables=[{"title": "Notes", "markdown": note}])
