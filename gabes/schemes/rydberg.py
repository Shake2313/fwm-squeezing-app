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
import functools

import numpy as np

from .. import atoms, constants, core, doppler, kernels, observables, species
from ..lineshape import window_fwhm
from ..report import derived_table
from .base import ParamSpec, Scheme

MHZ = 2 * np.pi * 1e6
PROBE_RABI = 1e-3
# 5P3/2 -> 40D5/2 "blue" Rydberg coupling-laser wavelength. Sets the residual
# two-photon Doppler ratio (k_probe - k_coupling)/k_probe for the counter-
# propagating geometry the static model assumes.
COUPLING_WAVELENGTH_NM = 480.0


@functools.lru_cache(maxsize=1)
def _probe_line():
    """Invariant 85Rb D2 probe-line constants (species-only), computed once."""
    rb85 = species.RB85
    Je, nu0, gamma_mhz, _, _ = rb85.line("D2")
    gamma_e = gamma_mhz * MHZ
    lam = constants.C_LIGHT / nu0
    k_vec = 2 * np.pi / lam
    dipole = np.sqrt(species.reduced_dipole_sq(gamma_e, lam, rb85.Jg, Je))
    return dict(nu0=nu0, gamma_mhz=gamma_mhz, gamma_e=gamma_e,
                lam=lam, k_vec=k_vec, dipole=dipole)


@functools.lru_cache(maxsize=1)
def _cascade_skeleton():
    """Constant 4-level cascade topology (decay + per-level Doppler ratios).

    Only the dephasing rates vary per call, so the topology is assembled once and
    `_atom` injects dephasing. The 40D / 39F two-photon levels carry the residual
    Doppler ratio (k_probe - k_coupling)/k_probe; the 5P intermediate level keeps
    the full probe-k ratio of 1.
    """
    line = _probe_line()
    gamma_r = 0.02 * MHZ
    kp = line["k_vec"]
    kc = 2 * np.pi / (COUPLING_WAVELENGTH_NM * 1e-9)
    ratio = (kp - kc) / kp
    return dict(
        labels=("5S F=3", "5P F'=4", "40D", "39F"),
        decay=((1, 0, line["gamma_e"]), (2, 1, gamma_r), (3, 2, gamma_r)),
        doppler_levels=(1, 2, 3),
        doppler_ratios=((2, ratio), (3, ratio)),
    )


class RydbergEITScheme(Scheme):
    name = "rydberg_eit"
    cluster = "A — Absorption"
    title = "Rydberg-EIT electrometry"
    caption = ("85Rb cascade EIT / microwave Autler-Townes electrometry. "
               "The static optical spectrum follows the 5S-5P-40D ladder; "
               "the 37 GHz RF leg dresses 40D-39F.")
    cache_version = "rydberg-eit-v2"
    defaults_version = "rydberg-eit-v2"

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
        rf_dephasing_mhz=1.00,
        doppler="off",
    )

    def param_schema(self):
        r = self._REF
        return [
            ParamSpec("view", "Regime", "Regime", "AT electrometry",
                      choices=("EIT", "AT electrometry"), control="segmented",
                      applies_defaults=True,
                      help="EIT: no microwave dressing. AT: the Rydberg RF leg is dressed."),
            ParamSpec("probe_power_uw", "Probe power (display only)", "Cell & beams",
                      r["probe_power_uw"], 0.1, 20.0, 0.1, "µW", recompute=False,
                      help="Display-only: the probe stays weak-probe linear "
                           "(Ω_probe = 1e-3 Γ), so its power only scales the readout, "
                           "not the solve."),
            ParamSpec("coupling_power_mw", "Coupling power", "Cell & beams",
                      r["coupling_power_mw"], 1.0, 80.0, 0.5, "mW",
                      help="480 nm 5P→40D coupling-beam power. Drives Ω_c via the "
                           "√(P/P_ref) intensity scaling at fixed waist, anchored so "
                           "the reference power reproduces Coupling Rabi."),
            ParamSpec("beam_diameter_mm", "Beam diameter", "Cell & beams",
                      r["beam_diameter_mm"], 0.05, 1.0, 0.01, "mm",
                      help="Coupling-beam 1/e² diameter. Ω_c ∝ 1/d at fixed power "
                           "(intensity ∝ P/d²)."),
            ParamSpec("cell_mm", "Cell length", "Cell & beams", r["cell_mm"],
                      1.0, 100.0, 0.5, "mm", recompute=False),
            ParamSpec("temp_c", "Temperature", "Cell & beams", r["temp_c"],
                      15.0, 80.0, 1.0, "°C"),
            ParamSpec("coupling_rabi_mhz", "Coupling Rabi (anchor)", "Fields",
                      r["coupling_rabi_mhz"], 0.1, 20.0, 0.1, "MHz", advanced=True,
                      help="Anchor Ω_c/2π at the reference power & waist. The "
                           "Coupling power and Beam diameter sliders scale the "
                           "effective Ω_c around this value."),
            ParamSpec("lo_rabi_mhz", "Microwave Rabi Ω_RF", "Fields",
                      r["lo_rabi_mhz"], 0.0, 20.0, 0.1, "MHz",
                      help="Rydberg 40D -> 39F dressing Rabi frequency Ω_RF/2π."),
            ParamSpec("mw_detuning_mhz", "Microwave detuning", "Detunings",
                      r["mw_detuning_mhz"], -20.0, 20.0, 0.1, "MHz"),
            ParamSpec("mw_frequency_ghz", "Microwave frequency (display only)", "Fields",
                      r["mw_frequency_ghz"], 1.0, 100.0, 0.1, "GHz", recompute=False,
                      help="Display-only metadata; not used in the solve."),
            ParamSpec("rydberg_dephasing_mhz", "Rydberg dephasing (5S–40D)", "Atomic",
                      r["rydberg_dephasing_mhz"], 0.0, 5.0, 0.01, "MHz",
                      help="Phenomenological 5S–40D ground–Rydberg coherence "
                           "broadening; sets the EIT linewidth."),
            ParamSpec("rf_dephasing_mhz", "RF dephasing (40D–39F)", "Atomic",
                      r["rf_dephasing_mhz"], 0.0, 5.0, 0.01, "MHz", advanced=True,
                      help="Phenomenological 40D–39F Rydberg–Rydberg coherence "
                           "broadening on the RF-dressed leg."),
            ParamSpec("if_khz", "IF offset", "Detection & scaling", r["if_khz"],
                      1.0, 500.0, 1.0, "kHz", advanced=True, recompute=False),
            ParamSpec("doppler", "Doppler treatment", "Numerics", r["doppler"],
                      choices=("off", "on"), advanced=True,
                      help="off: counter-propagating, Doppler-suppressed static model. "
                           "on: Maxwell-average the residual two-photon Doppler "
                           "(k_probe − k_coupling)·v, so cell temperature also broadens "
                           "the coherence, not just the vapor density."),
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

    def _atom(self, ground_deph, rf_deph):
        """4-level cascade with the constant topology from `_cascade_skeleton`;
        only the two dephasing channels (5S–40D ground–Rydberg and 40D–39F
        Rydberg–Rydberg) are injected per call."""
        sk = _cascade_skeleton()
        deph = (
            (0, 2, ground_deph), (2, 0, ground_deph),
            (0, 3, ground_deph), (3, 0, ground_deph),
            (2, 3, rf_deph), (3, 2, rf_deph),
        )
        return atoms.AtomModel(
            name="rydberg_eit_85rb",
            n_levels=4,
            labels=sk["labels"],
            ground=(0,),
            excited=(1, 2, 3),
            decay=sk["decay"],
            dephasing=deph,
            doppler_levels=sk["doppler_levels"],
            doppler_ratios=sk["doppler_ratios"],
        )

    def _coupling_rabi(self, params):
        """Effective coupling Rabi Ω_c/2π [MHz] from the 480 nm beam power and
        waist, FWM-style (fwm.py `_fields_from_params`): Ω ∝ E ∝ √(intensity) =
        √(P / area), area ∝ d². Anchored to the reference operating point, so at
        (P_ref, d_ref) it reproduces the fitted Coupling Rabi anchor exactly."""
        oc_ref = float(params.get("coupling_rabi_mhz", self._REF["coupling_rabi_mhz"]))
        p = float(params.get("coupling_power_mw", self._REF["coupling_power_mw"]))
        d = float(params.get("beam_diameter_mm", self._REF["beam_diameter_mm"]))
        p_ref = self._REF["coupling_power_mw"]
        d_ref = self._REF["beam_diameter_mm"]
        scale = np.sqrt(max(p, 0.0) / p_ref) * (d_ref / max(d, 1e-9))
        return oc_ref * scale

    def _scan(self, params):
        view = params.get("view", "AT electrometry")
        lo = float(params.get("lo_rabi_mhz", self._REF["lo_rabi_mhz"]))
        # Size the window off the EFFECTIVE Ω_c (the same power/waist-scaled value
        # that enters h_of), not the raw anchor — otherwise raising coupling power
        # or shrinking the beam widens the real feature past a pinned window and
        # the spectrum (and its slope/linewidth metrics) get clipped.
        oc = self._coupling_rabi(params)
        if view == "AT electrometry" and lo > 0:
            half = max(10.0, 3.0 * max(lo, oc))
        else:
            half = max(8.0, 3.0 * oc)
        return np.linspace(-half, half, 801) * MHZ

    def compute(self, params):
        line = _probe_line()
        gamma_e = line["gamma_e"]
        k_vec = line["k_vec"]
        dipole = line["dipole"]
        nu0 = line["nu0"]
        T = float(params.get("temp_c", self._REF["temp_c"])) + 273.15
        N = species.number_density(species.RB85, T)
        ground_deph = float(params.get(
            "rydberg_dephasing_mhz", self._REF["rydberg_dephasing_mhz"])) * MHZ
        rf_deph = float(params.get(
            "rf_dephasing_mhz", self._REF["rf_dephasing_mhz"])) * MHZ
        atom = self._atom(ground_deph, rf_deph)
        probe = PROBE_RABI * gamma_e
        Oc = self._coupling_rabi(params) * MHZ
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
        n = atom.n_levels
        # Doppler off: single static class (kv = 0). Doppler on: Maxwell-average
        # the residual two-photon shift carried by atom.S_v (per-level k ratios),
        # kv = k_probe·v — the same affine-kernel path the Λ scheme uses.
        if params.get("doppler", "off") == "on":
            v, w = doppler.velocity_grid(
                T, mass=constants.MASS_85RB, dv=2.0, cutoff_sigma=4.0)
            kv = k_vec * v
        else:
            w, kv = np.ones(1), np.zeros(1)
        if kernels.available():
            base = core.build_liouvillian(h_of(0.0), atom)
            A_coef = core.build_liouvillian(h_of(1.0), atom) - base
            with core.blas_single_thread():
                chi = kernels.affine_scan_chi(
                    base, A_coef, atom.S_v,
                    np.ascontiguousarray(scan, dtype=float),
                    np.ascontiguousarray(kv, dtype=float),
                    np.ascontiguousarray(w, dtype=float), 1 * n + 0, n) / probe
        else:
            chi = np.zeros(scan.size, dtype=complex)
            for i, s in enumerate(scan):
                L0 = core.build_liouvillian(h_of(s), atom)
                # delta = -kv reproduces the kernel's +kv·S_v velocity shift.
                rho = core.steady_state_batched(L0, -kv, atom.S_v, atom.n_levels)
                chi[i] = ((rho[:, 1, 0] / probe) * w).sum()

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
            gamma_mhz=line["gamma_mhz"],
            coupling_rabi_mhz=Oc / MHZ,
            lo_rabi_mhz=Olo / MHZ,
            rydberg_dephasing_mhz=ground_deph / MHZ,
            rf_dephasing_mhz=rf_deph / MHZ,
            mw_frequency_ghz=float(params.get(
                "mw_frequency_ghz", self._REF["mw_frequency_ghz"])),
            probe_power_uw=float(params.get("probe_power_uw", self._REF["probe_power_uw"])),
            beam_diameter_mm=float(params.get(
                "beam_diameter_mm", self._REF["beam_diameter_mm"])),
        )

    @staticmethod
    def _transparency_maxima(x, y, window_mhz):
        """Local transmission maxima (pos [MHz], height) within ±window of line
        centre — the dressed AT/EIT transparency peaks."""
        return [(float(x[i]), float(y[i])) for i in range(1, y.size - 1)
                if abs(x[i]) <= window_mhz and y[i] > y[i - 1] and y[i] >= y[i + 1]]

    def _readout(self, raw, params):
        """Cheap transmission/dispersion arrays + scalar metrics, with no
        matplotlib — the headless path tests and scans reuse without paying the
        figure-build cost. `observables` wraps this and draws the figure."""
        x = raw["scan"] / MHZ
        alpha, xphys = observables.absorption_coefficient(
            raw["chi_bar"], raw["k_vec"], raw["N"], dipole=raw["dipole"],
            line_strength=raw["ls"])
        # Cell length only scales αL here (recompute=False navigate-only knob), so
        # read it from the live params, not the cached raw, or a cell-length change
        # would not update the transmission.
        L = float(params.get("cell_mm", self._REF["cell_mm"])) * 1e-3
        T_trans = observables.transmission(alpha, L)
        ic = int(np.argmin(np.abs(x)))
        width = window_fwhm(x, T_trans, ic)
        slope = np.nanmax(np.abs(np.gradient(T_trans, x)))

        metrics = []
        if raw["lo_rabi_mhz"] > 0:
            dmw = abs(float(params.get("mw_detuning_mhz", self._REF["mw_detuning_mhz"])))
            window = max(8.0, 2.0 * raw["lo_rabi_mhz"] + dmw)
            peaks = self._transparency_maxima(x, T_trans, window)
            if len(peaks) >= 2:
                xs = sorted(p[0] for p in sorted(peaks, key=lambda p: p[1])[-2:])
                metrics.append(dict(label="RF AT splitting", value=f"{xs[1] - xs[0]:.2f} MHz",
                                    help="Separation of the two tallest dressed peaks."))
            if peaks:
                # Height-weighted centre of the transparency peaks: the symmetric
                # doublet sits at 0; a detuned microwave pulls it toward the stronger
                # dressed state, so the shift is readable even when one peak fades.
                wsum = sum(h for _, h in peaks)
                center = sum(px * h for px, h in peaks) / wsum
                metrics.append(dict(label="AT center shift", value=f"{center:+.2f} MHz",
                                    help="Height-weighted centre of the dressed "
                                         "transparency peaks; nonzero when the microwave "
                                         "is detuned off the 40D–39F resonance."))
        else:
            metrics.append(dict(label="EIT linewidth", value=f"{width:.2f} MHz",
                                help="FWHM of the central transparency feature."))
        metrics.extend([
            dict(label="Max spectral slope", value=f"{slope:.3f} /MHz",
                 help="Largest static dT/dnu; used internally for electrometry tests."),
            dict(label="Transmission at resonance", value=f"{T_trans[ic]:.3f}"),
        ])
        return dict(x=x, T_trans=T_trans, xphys=xphys, metrics=metrics)

    def observables(self, raw, params):
        import matplotlib.pyplot as plt
        ro = self._readout(raw, params)
        x, T_trans, xphys = ro["x"], ro["T_trans"], ro["xphys"]

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

        derived = derived_table([
            ("Ladder", "85Rb 5S1/2 F=3 → 5P3/2 F'=4 → 40D5/2"),
            ("RF leg", "40D5/2 → 39F7/2"),
            ("Microwave frequency", f"{raw['mw_frequency_ghz']:.1f} GHz"),
            ("Coupling Rabi Ω_c / 2π", f"{raw['coupling_rabi_mhz']:.3f} MHz"),
            ("Beam diameter", f"{raw['beam_diameter_mm']:.3f} mm"),
            ("Probe power (display only)", f"{raw['probe_power_uw']:.2f} µW"),
            ("Γ_5P / 2π", f"{raw['gamma_mhz']:.4f} MHz"),
            ("5S–40D dephasing / 2π", f"{raw['rydberg_dephasing_mhz']:.3f} MHz"),
            ("40D–39F dephasing / 2π", f"{raw['rf_dephasing_mhz']:.3f} MHz"),
            ("N(85Rb)", f"{raw['N']:.3e} /m³"),
        ])
        return dict(metrics=ro["metrics"], figure=fig, tables=[derived])

    def info(self):
        return (
            "**Rydberg-EIT electrometry.** Static 85Rb cascade-EIT model based on "
            "the photon-shot-noise-limited vapor-cell experiment. The displayed "
            "spectrum shows EIT and microwave Autler-Townes dressing only; "
            "time-domain superheterodyne demodulation is not simulated.\n\n"
            "The 480 nm coupling power and beam diameter drive Ω_c via the "
            "√(P/d²) intensity scaling, anchored to the reference operating point. "
            "An optional Doppler-on mode Maxwell-averages the residual two-photon "
            "shift (k_probe − k_coupling)·v for the counter-propagating geometry.\n\n"
            "Reference: arXiv:2606.04354, 85Rb vapor cell, 50 mm cell, 6 uW probe, "
            "30 mW coupling beam, 0.15 mm beam diameter, 37 GHz RF transition."
        )
