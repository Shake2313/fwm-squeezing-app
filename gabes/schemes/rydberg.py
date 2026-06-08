"""
Cluster A -- Rydberg-EIT electrometry.

v1 is a static 85Rb cascade model:
    5S1/2 F=3 -> 5P3/2 F'=4 -> 40D5/2, with a 37 GHz microwave leg
    40D5/2 -> 39F7/2.

The model intentionally stops at the optical spectrum and microwave
Autler-Townes splitting. Time-domain superheterodyne demodulation and public
experiment-match overlays are left out; reference sensitivity numbers are kept
as internal constants for tests.
"""
import numpy as np

from .. import atoms, constants, core, observables, species
from .base import ParamSpec, Scheme

MHZ = 2 * np.pi * 1e6
PROBE_RABI = 1e-3


def _window_fwhm(x, y, ic):
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


class RydbergEITScheme(Scheme):
    name = "rydberg_eit"
    cluster = "A — Absorption"
    title = "Rydberg-EIT electrometry"
    caption = ("85Rb cascade EIT / microwave Autler-Townes electrometry. "
               "The static optical spectrum follows the 5S-5P-40D ladder; "
               "the 37 GHz RF leg dresses 40D-39F.")
    cache_version = "rydberg-eit-v1"
    defaults_version = "rydberg-eit-v1"

    REFERENCE_SENSITIVITY_NV_CM_SQRT_HZ = 12.5
    REFERENCE_PSN_LIMIT_NV_CM_SQRT_HZ = 11.2
    REFERENCE_UNCERTAINTY_NV_CM_SQRT_HZ = 0.8

    _REF = dict(
        probe_power_uw=6.0,
        coupling_power_mw=30.0,
        beam_diameter_mm=0.15,
        cell_mm=50.0,
        temp_c=20.0,
        coupling_rabi_mhz=3.0,
        lo_rabi_mhz=3.7,
        mw_detuning_mhz=0.0,
        mw_frequency_ghz=37.0,
        if_khz=40.0,
        rydberg_dephasing_mhz=1.00,
        doppler="off",
    )

    def param_schema(self):
        r = self._REF
        return [
            ParamSpec("view", "Regime", "Regime", "AT electrometry",
                      choices=("EIT", "AT electrometry"), control="segmented",
                      applies_defaults=True,
                      help="EIT: no microwave dressing. AT: the Rydberg RF leg is dressed."),
            ParamSpec("probe_power_uw", "Probe power", "Cell & beams",
                      r["probe_power_uw"], 0.1, 20.0, 0.1, "uW", recompute=False,
                      help="Experimental optical readout power; the OBE remains weak-probe linear."),
            ParamSpec("coupling_power_mw", "Coupling power", "Cell & beams",
                      r["coupling_power_mw"], 1.0, 80.0, 0.5, "mW", recompute=False,
                      help="Experimental 480 nm beam power metadata; Omega_c is the fitted OBE knob."),
            ParamSpec("beam_diameter_mm", "Beam diameter", "Cell & beams",
                      r["beam_diameter_mm"], 0.05, 1.0, 0.01, "mm", recompute=False),
            ParamSpec("cell_mm", "Cell length", "Cell & beams", r["cell_mm"],
                      1.0, 100.0, 0.5, "mm", recompute=False),
            ParamSpec("temp_c", "Temperature", "Cell & beams", r["temp_c"],
                      15.0, 80.0, 1.0, "deg C"),
            ParamSpec("coupling_rabi_mhz", "Coupling Rabi", "Fields",
                      r["coupling_rabi_mhz"], 0.1, 20.0, 0.1, "MHz",
                      help="Optical 5P -> 40D coupling Rabi frequency Omega_c/2pi."),
            ParamSpec("lo_rabi_mhz", "Microwave LO Rabi", "Fields",
                      r["lo_rabi_mhz"], 0.0, 20.0, 0.1, "MHz",
                      help="Rydberg 40D -> 39F dressing Rabi frequency."),
            ParamSpec("mw_detuning_mhz", "Microwave detuning", "Detunings",
                      r["mw_detuning_mhz"], -20.0, 20.0, 0.1, "MHz"),
            ParamSpec("mw_frequency_ghz", "Microwave frequency", "Fields",
                      r["mw_frequency_ghz"], 1.0, 100.0, 0.1, "GHz", recompute=False),
            ParamSpec("rydberg_dephasing_mhz", "Rydberg dephasing", "Atomic",
                      r["rydberg_dephasing_mhz"], 0.0, 5.0, 0.01, "MHz",
                      help="Phenomenological 5S-40D coherence broadening."),
            ParamSpec("if_khz", "IF offset", "Detection & scaling", r["if_khz"],
                      1.0, 500.0, 1.0, "kHz", advanced=True, recompute=False),
            ParamSpec("doppler", "Doppler treatment", "Numerics", r["doppler"],
                      choices=("off",), advanced=True,
                      help="v1 uses the counter-propagating, Doppler-suppressed static model."),
        ]

    def _defaults(self, view):
        vals = dict(self._REF)
        vals["view"] = view
        if view == "EIT":
            vals["lo_rabi_mhz"] = 0.0
        return vals

    def recommended_defaults(self, params):
        return {
            "EIT": self._defaults("EIT"),
            "AT electrometry": self._defaults("AT electrometry"),
        }

    def _atom(self, gamma_e, rydberg_deph):
        gamma_r = 0.02 * MHZ
        deph = (
            (0, 2, rydberg_deph), (2, 0, rydberg_deph),
            (0, 3, rydberg_deph), (3, 0, rydberg_deph),
            (2, 3, rydberg_deph), (3, 2, rydberg_deph),
        )
        return atoms.AtomModel(
            name="rydberg_eit_85rb",
            n_levels=4,
            labels=("5S F=3", "5P F'=4", "40D", "39F"),
            ground=(0,),
            excited=(1, 2, 3),
            decay=((1, 0, gamma_e), (2, 1, gamma_r), (3, 2, gamma_r)),
            dephasing=deph,
            doppler_levels=(),
        )

    def _scan(self, params):
        view = params.get("view", "AT electrometry")
        lo = float(params.get("lo_rabi_mhz", self._REF["lo_rabi_mhz"]))
        oc = float(params.get("coupling_rabi_mhz", self._REF["coupling_rabi_mhz"]))
        if view == "AT electrometry" and lo > 0:
            half = max(10.0, 3.0 * max(lo, oc))
        else:
            half = max(8.0, 3.0 * oc)
        return np.linspace(-half, half, 801) * MHZ

    def compute(self, params):
        rb85 = species.RB85
        Je, nu0, gamma_mhz, _, _ = rb85.line("D2")
        gamma_e = gamma_mhz * MHZ
        lam = constants.C_LIGHT / nu0
        k_vec = 2 * np.pi / lam
        dipole = np.sqrt(species.reduced_dipole_sq(gamma_e, lam, rb85.Jg, Je))
        T = float(params.get("temp_c", self._REF["temp_c"])) + 273.15
        N = species.number_density(rb85, T)
        rydberg_deph = float(params.get(
            "rydberg_dephasing_mhz", self._REF["rydberg_dephasing_mhz"])) * MHZ
        atom = self._atom(gamma_e, rydberg_deph)
        probe = PROBE_RABI * gamma_e
        Oc = float(params.get("coupling_rabi_mhz", self._REF["coupling_rabi_mhz"])) * MHZ
        Olo = float(params.get("lo_rabi_mhz", self._REF["lo_rabi_mhz"])) * MHZ
        if params.get("view", "AT electrometry") == "EIT":
            Olo = 0.0
        Dmw = float(params.get("mw_detuning_mhz", self._REF["mw_detuning_mhz"])) * MHZ

        def h_of(s):
            H = np.zeros((4, 4), dtype=complex)
            H[1, 1] = -s
            H[2, 2] = -s
            H[3, 3] = -s - Dmw
            H[0, 1] = H[1, 0] = probe / 2
            H[1, 2] = H[2, 1] = Oc / 2
            H[2, 3] = H[3, 2] = Olo / 2
            return H

        scan = self._scan(params)
        chi = np.zeros(scan.size, dtype=complex)
        for i, s in enumerate(scan):
            L0 = core.build_liouvillian(h_of(s), atom)
            rho = core.steady_state_batched(L0, np.array([0.0]), atom.S_v, atom.n_levels)
            chi[i] = rho[0, 1, 0] / probe

        return dict(
            scan=scan,
            chi_bar=chi,
            N=N,
            T=T,
            L=float(params.get("cell_mm", self._REF["cell_mm"])) * 1e-3,
            ls=0.001,
            k_vec=k_vec,
            omega0=2 * np.pi * nu0,
            dipole=dipole,
            gamma_mhz=gamma_mhz,
            coupling_rabi_mhz=Oc / MHZ,
            lo_rabi_mhz=Olo / MHZ,
            rydberg_dephasing_mhz=rydberg_deph / MHZ,
            mw_frequency_ghz=float(params.get(
                "mw_frequency_ghz", self._REF["mw_frequency_ghz"])),
            probe_power_uw=float(params.get("probe_power_uw", self._REF["probe_power_uw"])),
            beam_diameter_mm=float(params.get(
                "beam_diameter_mm", self._REF["beam_diameter_mm"])),
        )

    @staticmethod
    def _at_split(x, y, center=0.0, window_mhz=8.0):
        local = []
        for i in range(1, y.size - 1):
            if abs(x[i] - center) > window_mhz:
                continue
            if y[i] > y[i - 1] and y[i] >= y[i + 1]:
                local.append((x[i], y[i]))
        left = [(xx, yy) for xx, yy in local if xx < center]
        right = [(xx, yy) for xx, yy in local if xx > center]
        if not left or not right:
            return float("nan")
        xl = max(left, key=lambda item: item[1])[0]
        xr = max(right, key=lambda item: item[1])[0]
        return float(abs(xr - xl))

    def observables(self, raw, params):
        import matplotlib.pyplot as plt
        x = raw["scan"] / MHZ
        alpha, xphys = observables.absorption_coefficient(
            raw["chi_bar"], raw["k_vec"], raw["N"], dipole=raw["dipole"],
            line_strength=raw["ls"])
        T_trans = observables.transmission(alpha, raw["L"])
        ic = int(np.argmin(np.abs(x)))
        width = _window_fwhm(x, T_trans, ic)
        slope = np.nanmax(np.abs(np.gradient(T_trans, x)))

        fig, (axT, axD) = plt.subplots(2, 1, figsize=(8.5, 6.4), sharex=True)
        axT.plot(x, T_trans, color="#0f766e", lw=1.8)
        axT.set_ylabel("Transmission")
        axT.set_ylim(-0.02, 1.04)
        axT.set_title(
            f"85Rb Rydberg {params.get('view', 'AT electrometry')}: "
            f"Omega_c = {raw['coupling_rabi_mhz']:.2f} MHz, "
            f"Omega_LO = {raw['lo_rabi_mhz']:.2f} MHz")
        axD.plot(x, np.real(xphys), color="#7c3aed", lw=1.5)
        axD.set_ylabel("Re chi")
        axD.set_xlabel("Probe detuning [MHz]")
        for a in (axT, axD):
            a.axvline(0.0, color="gray", ls=":", lw=0.8)
        fig.tight_layout()

        metrics = []
        if raw["lo_rabi_mhz"] > 0:
            split = self._at_split(x, T_trans)
            metrics.append(dict(label="RF AT splitting", value=f"{split:.2f} MHz",
                                help="Separation of the dressed Rydberg EIT peaks."))
        else:
            metrics.append(dict(label="EIT linewidth", value=f"{width:.2f} MHz",
                                help="FWHM of the central transparency feature."))
        metrics.extend([
            dict(label="Max spectral slope", value=f"{slope:.3f} /MHz",
                 help="Largest static dT/dnu; used internally for electrometry tests."),
            dict(label="Transmission at resonance", value=f"{T_trans[ic]:.3f}"),
        ])

        derived = (
            "| Quantity | Value |\n|---|---|\n"
            "| Ladder | 85Rb 5S1/2 F=3 -> 5P3/2 F'=4 -> 40D5/2 |\n"
            "| RF leg | 40D5/2 -> 39F7/2 |\n"
            f"| Microwave frequency | {raw['mw_frequency_ghz']:.1f} GHz |\n"
            f"| Probe power | {raw['probe_power_uw']:.2f} uW |\n"
            f"| Beam diameter | {raw['beam_diameter_mm']:.3f} mm |\n"
            f"| Gamma_5P/2pi | {raw['gamma_mhz']:.4f} MHz |\n"
            f"| Rydberg dephasing/2pi | {raw['rydberg_dephasing_mhz']:.3f} MHz |\n"
            f"| N(85Rb) | {raw['N']:.3e} /m^3 |\n"
        )
        return dict(metrics=metrics, figure=fig,
                    tables=[{"title": "Derived quantities", "markdown": derived}])

    def info(self):
        return (
            "**Rydberg-EIT electrometry.** Static 85Rb cascade-EIT model based on "
            "the photon-shot-noise-limited vapor-cell experiment. The displayed "
            "spectrum shows EIT and microwave Autler-Townes dressing only; "
            "time-domain superheterodyne demodulation is not simulated in v1.\n\n"
            "Reference: arXiv:2606.04354, 85Rb vapor cell, 50 mm cell, 6 uW probe, "
            "30 mW coupling beam, 0.15 mm beam diameter, 37 GHz RF transition."
        )
