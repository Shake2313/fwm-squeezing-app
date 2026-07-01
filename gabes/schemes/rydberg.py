"""
Cluster A -- Rydberg-EIT electrometry.

v1 is a static 85Rb cascade model:
    5S1/2 F=3 -> 5P3/2 F'=4 -> 40D5/2, with a 37 GHz microwave leg
    40D5/2 -> 39F7/2.

The model intentionally stops at the optical spectrum and microwave
Autler-Townes splitting. Full time-domain superheterodyne demodulation and public
experiment-match overlays are left out; a finite-IF discriminator proxy is
reported from the static spectrum. Reference sensitivity numbers are kept as
internal constants for tests.
"""
import functools

import numpy as np

from .. import atoms, beam, constants, core, doppler, kernels, observables, species
from ..report import derived_table
from .base import ExtraView, ParamSpec, Scheme

MHZ = 2 * np.pi * 1e6
PROBE_RABI = 1e-3
# 5P3/2 -> 40D5/2 "blue" Rydberg coupling-laser wavelength (Ju et al. Fig. 1a).
# Sets the residual two-photon Doppler ratio (k_probe - k_coupling)/k_probe for
# the counter-propagating geometry the static model assumes.
COUPLING_WAVELENGTH_NM = 481.0

# --- Calibration constants fitted to Ju et al. (arXiv:2606.04354) ---
# Transit-time broadening of the 5S-40D coherence: an atom crosses the beam in
# ~d/v_mp, so the coherence decays at ~v_mp/d. The O(1) factor is fit with the
# probe-power anchor so the reference (0.15 mm, 6 µW, 30 mW) lands on the paper's
# 1.6 MHz EIT linewidth (the transit-time term is the zero-probe floor, ~1.4 MHz).
TRANSIT_FACTOR = 0.6
# Weak-probe drive anchor: Ω_P/2π [MHz] at the reference 6 µW / 0.15 mm. A
# first-principles Ω_P from the bare D2 dipole overestimates by ~100x (optical
# pumping + Zeeman/hyperfine sub-structure are not in the lumped 2-level probe
# leg), so the probe drive is anchored like the coupling. Finite Ω_P power-
# broadens the EIT (Ju et al. Fig. 2b); fit so 6 µW sits at ~1.6 MHz. Because
# 6 µW is pinned to the narrow operating point, the broadening over 0–10 µW is
# milder than the paper's wide-range Fig. 2(b) (whose narrow point is at lower
# probe power).
PROBE_RABI_REF_MHZ = 2.0


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
    kc = beam.wavevector_from_wavelength_nm(COUPLING_WAVELENGTH_NM)
    ratio = beam.collinear_residual_k_ratio(kp, kc)
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
    cache_version = "rydberg-eit-v4"
    defaults_version = "rydberg-eit-v4"
    supports_headless_observables = True

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
        rydberg_dephasing_mhz=0.10,
        temp_dephasing_mhz_per_c=0.0,
        rf_dephasing_mhz=1.00,
        residual_zeeman_mhz=1.50,
        doppler="off",
    )

    def param_schema(self):
        r = self._REF
        return [
            ParamSpec("view", "Regime", "Regime", "EIT",
                      choices=("EIT", "AT electrometry"), control="segmented",
                      applies_defaults=True,
                      help="EIT: no microwave dressing (Ju et al. Fig. 2a). "
                           "AT: the Rydberg RF leg is dressed."),
            ParamSpec("probe_power_uw", "Probe power", "Cell & beams",
                      r["probe_power_uw"], 0.1, 20.0, 0.1, "µW",
                      help="780 nm probe power. Drives a weak-probe Ω_P via √(P)/d "
                           "(anchored), so raising it power-broadens the EIT line "
                           "(Ju et al. Fig. 2b). Reference 6 µW."),
            ParamSpec("coupling_power_mw", "Coupling power", "Cell & beams",
                      r["coupling_power_mw"], 1.0, 80.0, 0.5, "mW",
                      help="481 nm 5P→40D coupling-beam power. Drives Ω_c via the "
                           "√(P/P_ref) intensity scaling at fixed waist, anchored so "
                           "the reference power reproduces Coupling Rabi."),
            ParamSpec("beam_diameter_mm", "Beam diameter", "Cell & beams",
                      r["beam_diameter_mm"], 0.05, 1.0, 0.01, "mm",
                      help="Beam 1/e² diameter. Sets Ω_c, Ω_P ∝ 1/d (intensity ∝ "
                           "P/d²) AND the transit-time broadening ∝ v_thermal/d — "
                           "a tighter beam is brighter but transit-broadens the EIT."),
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
            ParamSpec("rydberg_dephasing_mhz", "Intrinsic dephasing (5S–40D)", "Atomic",
                      r["rydberg_dephasing_mhz"], 0.0, 5.0, 0.01, "MHz", advanced=True,
                      help="Intrinsic (non-transit) 5S–40D coherence broadening "
                           "from laser linewidth etc. The EIT linewidth is the sum "
                           "of this, the transit-time term (from beam diameter), and "
                           "the residual Zeeman term when uncompensated."),
            ParamSpec("temp_dephasing_mhz_per_c", "Temperature dephasing slope", "Atomic",
                      r["temp_dephasing_mhz_per_c"], 0.0, 0.2, 0.001, "MHz/°C",
                      advanced=True,
                      help="Optional phenomenological 5S–40D broadening added above "
                           "the reference temperature. Default 0 preserves the "
                           "Ju et al. reference line shape."),
            ParamSpec("residual_zeeman_mhz", "Residual Zeeman (uncompensated)", "Atomic",
                      r["residual_zeeman_mhz"], 0.0, 3.0, 0.01, "MHz", advanced=True,
                      help="Extra 5S–40D broadening present WITHOUT B-field "
                           "compensation (Ju et al. Fig. 2a blue curve). The EIT "
                           "figure overlays compensated (0) and uncompensated."),
            ParamSpec("rf_dephasing_mhz", "RF dephasing (40D–39F)", "Atomic",
                      r["rf_dephasing_mhz"], 0.0, 5.0, 0.01, "MHz", advanced=True,
                      help="Phenomenological 40D–39F Rydberg–Rydberg coherence "
                           "broadening on the RF-dressed leg."),
            ParamSpec("if_khz", "IF offset", "Detection & scaling", r["if_khz"],
                      1.0, 500.0, 1.0, "kHz", advanced=True, recompute=False),
            ParamSpec("doppler", "Doppler treatment", "Numerics", r["doppler"],
                      choices=("off", "on"), advanced=True,
                      help="off (default): the counter-propagating geometry suppresses "
                           "the two-photon Doppler — the transit-limited regime Ju et al. "
                           "operate in (1.6 MHz). on: Maxwell-average the residual "
                           "(k_probe − k_coupling)·v; in this lumped 4-level model that "
                           "over-broadens the EIT (~2.5 MHz floor), so it is a "
                           "what-if, not the calibrated reference."),
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
        """Effective coupling Rabi Ω_c/2π [MHz] from the 481 nm beam power and
        waist, FWM-style (fwm.py `_fields_from_params`): Ω ∝ E ∝ √(intensity) =
        √(P / area), area ∝ d². Anchored to the reference operating point, so at
        (P_ref, d_ref) it reproduces the fitted Coupling Rabi anchor exactly."""
        oc_ref = float(params.get("coupling_rabi_mhz", self._REF["coupling_rabi_mhz"]))
        p = float(params.get("coupling_power_mw", self._REF["coupling_power_mw"]))
        d = float(params.get("beam_diameter_mm", self._REF["beam_diameter_mm"]))
        p_ref = self._REF["coupling_power_mw"]
        d_ref = self._REF["beam_diameter_mm"]
        return beam.anchored_rabi_mhz(
            oc_ref, p, p_ref, diameter=d, ref_diameter=d_ref)

    def _probe_rabi(self, params):
        """Weak-probe drive Ω_P/2π [MHz] from probe power & waist, anchored at the
        reference like the coupling: Ω_P ∝ √(P)/d. Finite Ω_P makes the OBE
        saturate, so raising probe power broadens the EIT line (Ju et al. Fig 2b)."""
        p = float(params.get("probe_power_uw", self._REF["probe_power_uw"]))
        d = float(params.get("beam_diameter_mm", self._REF["beam_diameter_mm"]))
        p_ref = self._REF["probe_power_uw"]
        d_ref = self._REF["beam_diameter_mm"]
        return beam.anchored_rabi_mhz(
            PROBE_RABI_REF_MHZ, p, p_ref, diameter=d, ref_diameter=d_ref)

    def _transit_rate_mhz(self, params):
        """Transit-time broadening of the 5S–40D coherence /2π [MHz]: an atom
        crosses the beam in ~d/v_mp, so the coherence decays at TRANSIT_FACTOR·
        v_mp/d. Ties the beam diameter (and temperature) to the EIT linewidth —
        the transit-limited regime Ju et al. report."""
        T = float(params.get("temp_c", self._REF["temp_c"])) + 273.15
        d = float(params.get("beam_diameter_mm", self._REF["beam_diameter_mm"]))
        return beam.transit_broadening_mhz(
            T, d, mass=constants.MASS_85RB, factor=TRANSIT_FACTOR)

    def _scan_chi(self, atom, h_of, scan, probe, kv, w):
        """χ̄(scan) for one dephasing configuration via the affine kernel
        (Doppler-off kv=[0]; Doppler-on Maxwell grid). Shared by the compensated /
        uncompensated solves and the Fig. 2(b) probe-power sweep."""
        n = atom.n_levels
        if kernels.available():
            base = core.build_liouvillian(h_of(0.0), atom)
            A_coef = core.build_liouvillian(h_of(1.0), atom) - base
            with core.blas_single_thread():
                return kernels.affine_scan_chi(
                    base, A_coef, atom.S_v,
                    np.ascontiguousarray(scan, dtype=float),
                    np.ascontiguousarray(kv, dtype=float),
                    np.ascontiguousarray(w, dtype=float), 1 * n + 0, n) / probe
        chi = np.zeros(scan.size, dtype=complex)
        for i, s in enumerate(scan):
            L0 = core.build_liouvillian(h_of(s), atom)
            # delta = -kv reproduces the kernel's +kv·S_v velocity shift.
            rho = core.steady_state_batched(L0, -kv, atom.S_v, n)
            chi[i] = ((rho[:, 1, 0] / probe) * w).sum()
        return chi

    @staticmethod
    def _inhomogeneous(chi, scan, fwhm_mhz):
        """Gaussian inhomogeneous broadening of χ̄ over the probe detuning axis:
        a spread of resonance shifts (residual Zeeman) of FWHM `fwhm_mhz` averages
        χ(Δ − δ) over a Gaussian δ-distribution. Returns χ̄ unchanged for ~0 FWHM."""
        if fwhm_mhz <= 1e-6:
            return chi
        dscan = (scan[1] - scan[0]) / MHZ          # MHz per sample (uniform grid)
        sigma = fwhm_mhz / (2 * np.sqrt(2 * np.log(2)))
        half = int(np.ceil(4 * sigma / dscan))
        offs = np.arange(-half, half + 1) * dscan
        g = np.exp(-0.5 * (offs / sigma) ** 2)
        g /= g.sum()
        return np.convolve(chi, g, mode="same")

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
        k_vec = line["k_vec"]
        dipole = line["dipole"]
        nu0 = line["nu0"]
        T = float(params.get("temp_c", self._REF["temp_c"])) + 273.15
        N = species.number_density(species.RB85, T)

        # 5S-40D dephasing budget: intrinsic (laser etc.) + transit-time (beam
        # diameter) is always present; the residual Zeeman term is added only on
        # the uncompensated curve (B-field compensation off).
        intrinsic_base = float(params.get(
            "rydberg_dephasing_mhz", self._REF["rydberg_dephasing_mhz"]))
        temp_slope = float(params.get(
            "temp_dephasing_mhz_per_c", self._REF["temp_dephasing_mhz_per_c"]))
        temp_extra = max(float(params.get("temp_c", self._REF["temp_c"]))
                         - self._REF["temp_c"], 0.0) * max(temp_slope, 0.0)
        intrinsic = intrinsic_base + temp_extra
        transit = self._transit_rate_mhz(params)
        zeeman = float(params.get(
            "residual_zeeman_mhz", self._REF["residual_zeeman_mhz"]))
        rf_deph = float(params.get(
            "rf_dephasing_mhz", self._REF["rf_dephasing_mhz"])) * MHZ
        gd = (intrinsic + transit) * MHZ

        probe = self._probe_rabi(params) * MHZ
        Oc = self._coupling_rabi(params) * MHZ
        Olo = float(params.get("lo_rabi_mhz", self._REF["lo_rabi_mhz"])) * MHZ
        if params.get("view", "EIT") == "EIT":
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
        # Doppler off (default): single static class (kv = 0) — the suppressed-
        # residual, transit-limited regime. Doppler on: Maxwell-average the
        # residual two-photon shift carried by atom.S_v (per-level k ratios).
        if params.get("doppler", "off") == "on":
            v, w = doppler.velocity_grid(
                T, mass=constants.MASS_85RB, dv=2.0, cutoff_sigma=4.0)
            kv = k_vec * v
        else:
            w, kv = np.ones(1), np.zeros(1)

        # Compensated EIT: one solve. Uncompensated: residual Zeeman scatters the
        # m-sublevel shifts, so the EIT line is INHOMOGENEOUSLY broadened — model
        # it as a Gaussian (FWHM = residual_zeeman) convolution of the compensated
        # susceptibility over the probe detuning. (A homogeneous dephasing knob
        # only lowers the peak in this Ω_c-limited regime; it does not widen it.)
        chi = self._scan_chi(self._atom(gd, rf_deph), h_of, scan, probe, kv, w)
        chi_uncomp = self._inhomogeneous(chi, scan, zeeman)

        return dict(
            scan=scan,
            chi_bar=chi,
            chi_bar_uncomp=chi_uncomp,
            N=N,
            T=T,
            L=float(params.get("cell_mm", self._REF["cell_mm"])) * 1e-3,
            ls=0.001,
            k_vec=k_vec,
            omega0=2 * np.pi * nu0,
            dipole=dipole,
            gamma_mhz=line["gamma_mhz"],
            coupling_rabi_mhz=Oc / MHZ,
            probe_rabi_mhz=probe / MHZ,
            lo_rabi_mhz=Olo / MHZ,
            rydberg_dephasing_mhz=intrinsic_base,
            temperature_dephasing_mhz=temp_extra,
            transit_mhz=transit,
            residual_zeeman_mhz=zeeman,
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

    def _transmission(self, chi, raw, params):
        """(x [MHz], Beer-Lambert transmission, χ_phys) for one coherence array."""
        x = raw["scan"] / MHZ
        alpha, xphys = observables.absorption_coefficient(
            chi, raw["k_vec"], raw["N"], dipole=raw["dipole"], line_strength=raw["ls"])
        # Cell length only scales αL here (recompute=False navigate-only knob), so
        # read it from the live params, not the cached raw, or a cell-length change
        # would not update the transmission.
        L = float(params.get("cell_mm", self._REF["cell_mm"])) * 1e-3
        return x, observables.transmission(alpha, L), xphys

    @staticmethod
    def _eit_features(x, T_trans):
        """(FWHM [MHz], peak contrast) of the central EIT transparency feature.

        Contrast = T(centre) − absorptive floor in a ±5 MHz window. The FWHM is
        measured to the half-contrast level with linear interpolation of the two
        crossings (sub-sample accurate — a bare grid walk stair-steps the
        probe-power sweep). The transmitted-signal amplitude in Fig. 2(b) scales
        as probe power × contrast."""
        ic = int(np.argmin(np.abs(x)))
        peak = float(T_trans[ic])
        win = np.abs(x) <= 5.0
        floor = float(np.min(T_trans[win])) if win.any() else float(np.min(T_trans))
        contrast = peak - floor
        if contrast <= 0:
            return float("nan"), 0.0
        half = floor + 0.5 * contrast

        def crossing(step):
            i = ic
            while 0 < i < x.size - 1 and T_trans[i] >= half:
                i += step
            x0, y0, x1, y1 = x[i - step], T_trans[i - step], x[i], T_trans[i]
            return x1 if y1 == y0 else x0 + (half - y0) * (x1 - x0) / (y1 - y0)

        return abs(crossing(1) - crossing(-1)), contrast

    def _readout(self, raw, params):
        """Cheap transmission/dispersion arrays + scalar metrics, with no
        matplotlib — the headless path tests and scans reuse without paying the
        figure-build cost. `observables` wraps this and draws the figure."""
        x, T_trans, xphys = self._transmission(raw["chi_bar"], raw, params)
        ic = int(np.argmin(np.abs(x)))
        width, _ = self._eit_features(x, T_trans)
        slope = np.nanmax(np.abs(np.gradient(T_trans, x)))
        if_delta = max(float(params.get("if_khz", self._REF["if_khz"])) / 1000.0,
                       1e-9)
        if_valid = (x >= x[0] + if_delta) & (x <= x[-1] - if_delta)
        if if_valid.any():
            t_hi = np.interp(x[if_valid] + if_delta, x, T_trans)
            t_lo = np.interp(x[if_valid] - if_delta, x, T_trans)
            if_readout = (t_hi - t_lo) / (2.0 * if_delta)
            i_if = int(np.nanargmax(np.abs(if_readout)))
            if_disc = float(abs(if_readout[i_if]))
            if_detuning = float(x[if_valid][i_if])
        else:
            if_disc = float("nan")
            if_detuning = float("nan")

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
            dict(label="IF discriminator", value=f"{if_disc:.3f} /MHz",
                 help="Finite-difference transmission discriminator at the selected "
                      "IF offset; a static proxy for a lock-in/superhet readout."),
            dict(label="IF optimum detuning", value=f"{if_detuning:+.2f} MHz",
                 help="Probe detuning where the finite-IF discriminator is largest."),
            dict(label="Transmission at resonance", value=f"{T_trans[ic]:.3f}"),
        ])
        return dict(x=x, T_trans=T_trans, xphys=xphys, width=width, metrics=metrics)

    def observables(self, raw, params, include_figures=True):
        ro = self._readout(raw, params)
        x, T_trans, xphys, width = ro["x"], ro["T_trans"], ro["xphys"], ro["width"]
        view = params.get("view", "EIT")

        fig = None
        if include_figures:
            import matplotlib.pyplot as plt

            if view == "EIT":
                # Ju et al. Fig. 2(a): EIT transmission with / without B-field
                # compensation, zoomed onto the transparency feature.
                _, T_uncomp, _ = self._transmission(
                    raw["chi_bar_uncomp"], raw, params)
                fig, axT = plt.subplots(figsize=(8.5, 4.8))
                axT.plot(x, T_uncomp, color="#1f77b4", lw=1.6,
                         label="without B-field compensation")
                axT.plot(x, T_trans, color="#d62728", lw=2.0,
                         label="with compensation")
                axT.set_xlabel("Frequency [MHz]")
                axT.set_ylabel("Transmission")
                axT.axvline(0.0, color="gray", ls=":", lw=0.8)
                axT.legend(fontsize=9, loc="upper right")
                axT.set_title(
                    f"85Rb Rydberg-EIT: "
                    f"Omega_c = {raw['coupling_rabi_mhz']:.2f} MHz, "
                    f"probe = {raw['probe_power_uw']:.1f} uW")
                xlim = max(2.5, 2.5 * width)
                axT.set_xlim(-xlim, xlim)
            else:
                fig, (axT, axD) = plt.subplots(2, 1, figsize=(8.5, 6.4),
                                               sharex=True)
                axT.plot(x, T_trans, color="#0f766e", lw=1.8)
                axT.set_ylabel("Transmission")
                axT.set_ylim(-0.02, 1.04)
                axT.set_title(
                    f"85Rb Rydberg {view}: "
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
            ("Probe Rabi Ω_P / 2π", f"{raw['probe_rabi_mhz']:.3f} MHz"),
            ("Beam diameter", f"{raw['beam_diameter_mm']:.3f} mm"),
            ("Probe power", f"{raw['probe_power_uw']:.2f} µW"),
            ("Transit broadening / 2π", f"{raw['transit_mhz']:.3f} MHz"),
            ("Residual Zeeman (uncomp.) / 2π", f"{raw['residual_zeeman_mhz']:.3f} MHz"),
            ("Intrinsic 5S–40D dephasing / 2π", f"{raw['rydberg_dephasing_mhz']:.3f} MHz"),
            ("Temperature dephasing / 2π", f"{raw['temperature_dephasing_mhz']:.3f} MHz"),
            ("N(85Rb)", f"{raw['N']:.3e} /m³"),
        ])
        return dict(metrics=ro["metrics"], figure=fig, tables=[derived])

    def extra_views(self):
        """Ju et al. Fig. 2(b): EIT peak amplitude and linewidth vs probe power,
        with / without B-field compensation."""
        def _compute_sweep(params):
            powers = np.linspace(0.5, 10.0, 16)
            comp = {"width": [], "amp": []}
            uncomp = {"width": [], "amp": []}
            for p_uw in powers:
                pr = dict(params, probe_power_uw=float(p_uw), view="EIT")
                raw = self.compute(pr)
                x, T_c, _ = self._transmission(raw["chi_bar"], raw, pr)
                x, T_u, _ = self._transmission(raw["chi_bar_uncomp"], raw, pr)
                wc, ac = self._eit_features(x, T_c)
                wu, au = self._eit_features(x, T_u)
                comp["width"].append(wc); comp["amp"].append(p_uw * ac)
                uncomp["width"].append(wu); uncomp["amp"].append(p_uw * au)
            return dict(powers=powers.tolist(), comp=comp, uncomp=uncomp)

        def _render_sweep(s):
            import matplotlib.pyplot as plt
            powers = np.array(s["powers"])
            figF, (aA, aW) = plt.subplots(2, 1, figsize=(8.0, 6.4), sharex=True)
            for ax in (aA, aW):
                ax.grid(alpha=0.3)
            aA.plot(powers, s["uncomp"]["amp"], "s-", color="#1f77b4",
                    label="without compensation")
            aA.plot(powers, s["comp"]["amp"], "o-", color="#d62728",
                    label="with compensation")
            aA.set_ylabel("EIT peak amplitude [arb.]")
            aA.legend(fontsize=9)
            aW.plot(powers, s["uncomp"]["width"], "s-", color="#1f77b4")
            aW.plot(powers, s["comp"]["width"], "o-", color="#d62728")
            aW.set_ylabel("EIT linewidth [MHz]")
            aW.set_xlabel("Probe power [uW]")
            figF.tight_layout()
            return figF

        return [ExtraView(
            key="Fig. 2(b): probe-power dependence (peak amplitude & linewidth)",
            description="Sweeps the probe power 0.5–10 µW and extracts the EIT peak "
                        "amplitude and linewidth with and without B-field "
                        "compensation — the power broadening of Ju et al. Fig. 2(b).",
            compute=_compute_sweep, render=_render_sweep,
        )]

    def info(self):
        return (
            "**Rydberg-EIT electrometry.** 85Rb cascade-EIT / microwave Autler-Townes "
            "model calibrated to the photon-shot-noise-limited vapor-cell experiment "
            "of Ju et al. The EIT view reproduces Fig. 2(a) (transmission with / "
            "without B-field compensation); the probe-power panel reproduces "
            "Fig. 2(b). Time-domain superheterodyne demodulation is not simulated.\n\n"
            "Probe (780 nm) and coupling (481 nm) powers drive Ω_P, Ω_c via √(P)/d "
            "intensity scaling (anchored), and the beam diameter also sets the "
            "transit-time broadening that limits the EIT linewidth (≈1.6 MHz at the "
            "reference). The counter-propagating geometry suppresses the two-photon "
            "Doppler; an optional Doppler-on mode shows the residual "
            "(k_probe − k_coupling)·v broadening.\n\n"
            "**References**\n"
            "- [arXiv:2606.04354](https://arxiv.org/abs/2606.04354) — \"Photon "
            "shot-noise-limited Rydberg-EIT electrometry\", Ju et al. (85Rb vapor "
            "cell, 50 mm cell, 6 µW probe, 30 mW coupling, 0.15 mm beam, 37 GHz RF)."
        )
