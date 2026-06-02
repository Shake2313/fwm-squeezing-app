"""
Cluster B — saturated absorption spectroscopy (SAS).

A strong pump and a weak probe from the same laser counter-propagate through the
cell, so a moving atom sees them with opposite Doppler shifts:
    pump atom-frame detuning  = Δ + k·v
    probe atom-frame detuning = Δ − k·v        (Δ = ω − ω₀, the scan)
The pump burns a hole in the velocity distribution (saturates the population) at
its resonant class; the weak probe samples that population. At line centre both
address the v = 0 atoms → a sub-Doppler "Lamb dip". With two excited transitions
a "crossover" dip appears at the midpoint.

Model (standard hole-burning / rate picture): the pump-saturated populations come
from a full steady-state OBE solve per velocity class (core.steady_state_batched,
Δ_eff = Δ + k·v); the probe contributes the bare weak-probe absorption lineshape
(same scale as the OD scheme) weighted by the pump-saturated population difference:

    α_SAS(Δ) = Σ_v f(v) Σ_i (ρ_gg − ρ_ee_i)(v) · α_probe((Δ − k·v) − offset_i)

The no-pump limit (ρ_gg − ρ_ee → 1) reproduces the Doppler-broadened Voigt of OD.
"""
import numpy as np

from .. import atoms, constants, doppler, observables
from ..constants import GAMMA, K_VEC
from .. import core
from .base import ParamSpec, Preset, Scheme

PROBE_RABI = 1e-3                       # weak probe, in units of Γ
GAMMA_MHZ = GAMMA / (2 * np.pi) / 1e6


class SASScheme(Scheme):
    name = "sas"
    cluster = "B — Sub-Doppler"
    title = "Saturated absorption (SAS)"
    caption = ("Counter-propagating pump + weak probe. The pump burns velocity-class "
               "holes, giving sub-Doppler Lamb dips (and a crossover for two lines) "
               "on the Doppler-broadened line.")

    def param_schema(self):
        return [
            ParamSpec("pump_rabi", "Pump Rabi Ω_pump", "Fields", 2.0,
                      0.1, 10.0, 0.1, "Γ", help="Saturating counter-propagating beam."),
            ParamSpec("transitions", "Transitions", "Atomic", "single line",
                      choices=("single line", "two lines (crossover)")),
            ParamSpec("splitting", "Excited splitting", "Atomic", 60.0,
                      5.0, 200.0, 1.0, "Γ", help="Only used for two lines."),
            ParamSpec("temp_c", "Temperature", "Cell & beams", 50.0, 20.0, 200.0, 1.0, "°C"),
            ParamSpec("cell_mm", "Cell length", "Cell & beams", 10.0, 0.5, 200.0, 0.5, "mm"),
            ParamSpec("line_strength", "Line-strength factor", "Detection & scaling", 1.0,
                      0.01, 2.0, 0.01, ""),
            ParamSpec("scan_points", "Scan points", "Numerics", 1501,
                      401, 3001, 100, "", advanced=True),
        ]

    def presets(self):
        return [
            Preset("Single-line Lamb dip",
                   values=dict(transitions="single line", pump_rabi=2.0,
                               temp_c=50.0, cell_mm=10.0), icon="🎯"),
            Preset("Crossover (two lines)",
                   values=dict(transitions="two lines (crossover)", splitting=60.0,
                               pump_rabi=2.0, temp_c=50.0, cell_mm=10.0), icon="✛"),
        ]

    def _config(self, params):
        two = params["transitions"].startswith("two")
        if two:
            split = params["splitting"] * GAMMA
            return 2, np.array([-split / 2, split / 2])
        return 1, np.array([0.0])

    def compute(self, params):
        n_exc, offsets = self._config(params)
        atom = atoms.sas_atom(n_exc)
        Op = params["pump_rabi"] * GAMMA

        # Pump Hamiltonian (fixed; velocity enters only via Δ_eff = Δ + k·v).
        Hp = np.zeros((atom.n_levels, atom.n_levels), dtype=complex)
        for i, e in enumerate(atom.excited):
            Hp[e, e] = offsets[i]
            Hp[0, e] = Hp[e, 0] = Op / 2
        L0_pump = core.build_liouvillian(Hp, atom)

        T = params["temp_c"] + 273.15
        sigma_v = np.sqrt(constants.KB * T / constants.MASS_85RB)
        dopp_fwhm = np.sqrt(8 * np.log(2)) * K_VEC * sigma_v
        v, wt = doppler.velocity_grid(T, dv=3.0, cutoff_sigma=3.5)
        kv = K_VEC * v
        N = atoms.rb85_density(T)

        off_span = float(np.abs(offsets).max())
        half = max(3.5 * dopp_fwhm, off_span + 0.4 * dopp_fwhm)
        scan = np.linspace(-half, half, int(params["scan_points"]))

        # Bare weak-probe absorption lineshape (single 2-level transition), same
        # scale as the OD scheme; sampled by interpolation at (Δ − k·v) − offset_i.
        two_lvl = atoms.two_level()
        Hpr = np.zeros((2, 2), dtype=complex)
        Hpr[0, 1] = Hpr[1, 0] = PROBE_RABI * GAMMA / 2
        L0_probe = core.build_liouvillian(Hpr, two_lvl)
        kvmax = float(np.abs(kv).max())
        flo, fhi = scan.min() - kvmax - off_span, scan.max() + kvmax + off_span
        fine = np.linspace(flo, fhi, int((fhi - flo) / (GAMMA / 20)) + 2)
        rho_pr = core.steady_state_batched(L0_probe, fine, two_lvl.S_v, 2)
        chi_pr = rho_pr[:, 1, 0] / (PROBE_RABI * GAMMA)
        alpha_L, _ = observables.absorption_coefficient(
            chi_pr, K_VEC, N, line_strength=params["line_strength"])

        # α_SAS(Δ) = Σ_v f(v) Σ_i (ρ_gg − ρ_ee_i)(v) · α_probe((Δ − k·v) − offset_i)
        alpha = np.zeros(scan.size)
        excited = atom.excited
        for j, D in enumerate(scan):
            rho = core.steady_state_batched(L0_pump, D + kv, atom.S_v, atom.n_levels)
            rho_gg = rho[:, 0, 0].real
            contrib = np.zeros(v.size)
            for i, e in enumerate(excited):
                w_i = rho_gg - rho[:, e, e].real
                contrib += w_i * np.interp((D - kv) - offsets[i], fine, alpha_L)
            alpha[j] = float((wt * contrib).sum())

        return dict(scan=scan, alpha=alpha, L=params["cell_mm"] * 1e-3,
                    dopp_fwhm=dopp_fwhm, N=N, sigma_v=sigma_v,
                    offsets=offsets, two=(n_exc == 2))

    def observables(self, raw, params):
        import matplotlib.pyplot as plt
        x = raw["scan"] / (2 * np.pi) / 1e6                    # MHz
        alpha = raw["alpha"]
        T_trans = observables.transmission(alpha, raw["L"])
        OD = observables.optical_density(alpha, raw["L"])
        offs_mhz = raw["offsets"] / (2 * np.pi) / 1e6

        fig, (axT, axA) = plt.subplots(2, 1, figsize=(8.5, 6.4), sharex=True)
        axT.plot(x, T_trans, color="#1f77b4", lw=1.6)
        axT.set_ylabel("Transmission")
        axT.set_title(f"SAS: Ω_pump = {params['pump_rabi']:.1f} Γ,  T = {params['temp_c']:.0f} °C")
        axA.plot(x, OD, color="#d62728", lw=1.6)
        axA.set_ylabel("Optical density")
        axA.set_xlabel("Probe detuning  [MHz]")
        for ax in (axT, axA):
            for off in offs_mhz:
                ax.axvline(off, color="gray", ls=":", lw=0.7)
            if raw["two"]:
                ax.axvline(0.0, color="green", ls=":", lw=0.7)   # crossover midpoint
        fig.tight_layout()

        # Sub-Doppler dip at the first transition centre (transmission peak).
        ic = int(np.argmin(np.abs(x - offs_mhz[0])))
        sub_fwhm = _peak_fwhm(x, T_trans, ic)
        dopp_mhz = raw["dopp_fwhm"] / (2 * np.pi) / 1e6
        metrics = [
            dict(label="Doppler FWHM", value=f"{dopp_mhz:.0f} MHz",
                 help="Width of the Doppler-broadened background."),
            dict(label="Sub-Doppler dip FWHM", value=f"{sub_fwhm:.1f} MHz",
                 help="Lamb-dip width (≪ Doppler — natural / power-broadened)."),
            dict(label="Peak OD", value=f"{np.nanmax(OD):.2f}"),
        ]
        note = ("Two transitions: Lamb dips at ±splitting/2 and a **crossover** dip "
                "at the midpoint (green line)." if raw["two"]
                else "Single transition: one Lamb dip at line centre.")
        return dict(metrics=metrics, figure=fig,
                    tables=[{"title": "Notes", "markdown": note}])


def _peak_fwhm(x, y, ic):
    """FWHM of a narrow transmission peak at index ic above its adjacent floor."""
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
