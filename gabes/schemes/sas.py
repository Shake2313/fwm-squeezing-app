"""
Cluster A — absorption spectroscopy (OD / SAS), one scheme.

A weak probe measures the vapor's absorption. A counter-propagating *pump* of
power P saturates and optically pumps the velocity class it is resonant with;
the probe then samples those prepared populations:

    pump  atom-frame detuning = Δ + k·v
    probe atom-frame detuning = Δ − k·v        (Δ = laser − line centre, the scan)

  * **Pump off (P = 0) → OD.** The probe sees the linear, Doppler-broadened
    multi-line absorption — exactly the validated 85Rb D1 hyperfine spectrum
    (AutoOD calibration) at pump = 0, generalised to any isotope/line.
  * **Pump on → SAS.** Velocity-selective hole burning + hyperfine optical
    pumping give Doppler-free Lamb dips and (enhanced / inverted) crossovers.

Model
-----
For a chosen isotope/line `gabes.species` builds the full {Fg}×{Fe} hyperfine
manifold (CG-branched decay + transit-time relaxation). The absorption is

    α(δ) = Σ_components Σ_(Fg→Fe)  A_(Fg,Fe) · ĝ_(Fg,Fe)(δ)

  A_t = ∫α_t dδ = ls·π·k·p_Fg·C_F²·|d|²·N/(ε₀ℏ)/(2(2I+1))   (AutoOD absolute area)
  ĝ_t(δ) = Σ_v f(v)·[(ρ_Fg − ρ_Fe)(Δ+k·v)/p_Fg]·L̂((δ − k·v) − ω_t)   (unit area at P=0)

ρ comes from one OBE steady state per velocity class; the pump Hamiltonian is
scan-independent (the scan enters only via Δ_eff = Δ + k·v), so it is solved once
on a fine Δ_eff grid and interpolated. At P=0 the bracket is 1, ĝ_t is a unit
Voigt, and α reduces to the validated OD spectrum (∫α_t = A_t).

A generic Γ-unit hole-burning toy (one ground + one/two excited states) is kept
under Advanced for the bare picture without atomic data.

References:  Smith & Hughes, Am. J. Phys. 72, 631 (2004) (hyperfine pumping);
Preston, Am. J. Phys. 64, 1432 (1996); Steck alkali D-line data; the lab AutoOD
calculator (pump-off 85Rb D1 absolute scale).
"""
import numpy as np

from .. import atoms, constants, doppler, observables, species
from ..constants import GAMMA, K_VEC
from .. import core
from ..lineshape import narrowest_subdoppler, window_fwhm
from .base import ParamSpec, Scheme

PROBE_RABI = 1e-3                       # weak probe, in units of Γ
GAMMA_MHZ = GAMMA / (2 * np.pi) / 1e6
GENERIC = "Generic (Γ units)"
_trapz = getattr(np, "trapezoid", getattr(np, "trapz", None))   # numpy ≥2.0 rename


class SASScheme(Scheme):
    name = "sas"
    cluster = "A — Absorption"
    title = "Absorption spectroscopy (OD / SAS)"
    cache_version = "2"            # merged model: bust the old SAS compute cache
    defaults_version = "2"         # new param schema: reseed sidebar defaults
    supports_headless_observables = True
    caption = ("Weak-probe absorption with a counter-propagating pump. Pump off → "
               "linear Doppler-broadened OD (validated 85Rb D1 hyperfine scale); "
               "pump on → Doppler-free saturated-absorption Lamb dips and crossovers "
               "with hyperfine optical pumping. ⁸⁵Rb / ⁸⁷Rb / ¹³³Cs · D1/D2 or natural Rb.")

    def param_schema(self):
        return [
            ParamSpec("pump_power_mw", "Pump beam power", "Pump", 0.5, 0.0, 2.0, 0.01,
                      "mW", endpoints=("◀ OD", "SAS ▶"),
                      help="Counter-propagating saturating beam. Pull to 0 → linear "
                      "absorption (OD); raise → Doppler-free SAS features. Converted "
                      "to a Rabi frequency via the beam waist and I_sat (see About)."),
            ParamSpec("species", "Atom / isotope", "Atomic", "Rb (natural)",
                      choices=tuple(species.SPECIES_ORDER) + (GENERIC,),
                      help="Natural Rb overlays ⁸⁵Rb+⁸⁷Rb by abundance. Generic = bare "
                      "Γ-unit hole-burning model (no atomic data)."),
            ParamSpec("line", "Transition line", "Atomic", "D1",
                      choices=("D1", "D2"),
                      help="D1 (nP₁/₂) or D2 (nP₃/₂). Sets the excited hyperfine "
                      "manifold; ignored in Generic mode."),
            ParamSpec("temp_c", "Temperature", "Cell & beams", 40.0, 20.0, 200.0, 1.0, "°C",
                      help="Sets the vapor density (absorption scale) and Doppler width."),
            ParamSpec("cell_mm", "Cell length", "Cell & beams", 75.0, 0.5, 200.0, 0.5, "mm",
                      recompute=False),
            ParamSpec("ne_pressure_torr", "Ne buffer pressure", "Cell & beams", 0.0,
                      0.0, 200.0, 1.0, "Torr", advanced=True,
                      help="Fixed-neon pressure broadening only; pressure shift and "
                      "Dicke narrowing are not included."),
            ParamSpec("waist_mm", "Pump beam waist (1/e²)", "Pump", 1.0, 0.1, 5.0, 0.05,
                      "mm", advanced=True, help="Sets the pump intensity I = 2P/(πw²) "
                      "for the power→Rabi conversion."),
            ParamSpec("transit_khz", "Transit relaxation γ_t/2π", "Atomic", 100.0,
                      5.0, 2000.0, 5.0, "kHz", advanced=True,
                      help="Atoms leaving/entering the beam relax toward the thermal "
                      "ground state — regularises hyperfine pumping. Smaller γ_t → "
                      "stronger inverted crossovers."),
            ParamSpec("line_strength", "Line-strength factor", "Detection & scaling", 1.0,
                      0.01, 2.0, 0.01, "", advanced=True, recompute=False,
                      help="Effective |d|² calibration knob. =1.0 reproduces the "
                      "AutoOD-validated 85Rb D1 absolute scale at pump = 0."),
            ParamSpec("transitions", "Generic: transitions", "Atomic", "single line",
                      choices=("single line", "two lines (crossover)"), advanced=True,
                      help="Generic mode only."),
            ParamSpec("splitting", "Generic: excited splitting", "Atomic", 60.0,
                      5.0, 200.0, 1.0, "Γ", advanced=True, help="Generic mode only."),
            ParamSpec("scan_points", "Scan points", "Numerics", 1401,
                      401, 4001, 100, "", advanced=True),
        ]

    def recommended_defaults(self, params):
        """Two labelled presets — 'OD default' (pump off) and 'SAS default' (pump
        on) — for the current atom/line. Same temperature & cell; they differ only
        in the pump power."""
        if params.get("species", "Rb (natural)") == GENERIC:
            base, sas_power = dict(temp_c=50.0, cell_mm=10.0), 0.5
        else:
            rec = species.recommended(params.get("species", "Rb (natural)"),
                                      params.get("line", "D1"))
            base = dict(temp_c=rec["temp_c"], cell_mm=rec["cell_mm"])
            sas_power = rec["pump_power_mw"]
        return {"OD default": {**base, "pump_power_mw": 0.0},
                "SAS default": {**base, "pump_power_mw": sas_power}}

    def info(self):
        return (
            "**Absorption spectroscopy (OD / SAS).** A weak probe measures the "
            "vapor absorption; a counter-propagating pump of power *P* saturates "
            "the resonant velocity class.\n\n"
            "- **Pump off (P = 0):** linear Doppler-broadened absorption (OD). For "
            "85Rb D1 the absolute scale reproduces the lab AutoOD calculator "
            "(`references/AutoOD/`) to <0.1 %.\n"
            "- **Pump on:** Doppler-free Lamb dips + crossovers; CG-branched decay "
            "pumps population into the other ground hyperfine state, enhancing / "
            "inverting the crossovers (the dominant feature of real alkali SAS).\n\n"
            "The **probe** is fixed weak (Ω_probe = 1e-3 Γ) — its only role is to "
            "read the populations linearly; it is not a user knob. The **pump power** "
            "maps to a Rabi frequency via I = 2P/(πw²), Ω = Γ·√(I/2I_sat) with "
            "I_sat = 4.484 mW/cm² and the beam waist *w* (Advanced).\n\n"
            "Atomic data (hyperfine A/B, line centres, masses, linewidths) from the "
            "Steck D-line data sheets; Wigner-6j/3j line strengths in the AutoOD "
            "convention. Rb densities use the CRC vapor pressure (AutoOD), Cs the "
            "Steck fit.\n\n"
            "**References**\n"
            "- D. A. Smith & I. G. Hughes, *Am. J. Phys.* **72**, 631 (2004).\n"
            "- D. W. Preston, *Am. J. Phys.* **64**, 1432 (1996).\n"
            "- D. A. Steck, *Rubidium 85 / 87 & Cesium D Line Data*, http://steck.us/alkalidata."
        )

    # =================================================================
    # compute  (dispatch)
    # =================================================================
    def compute(self, params):
        if params.get("species", "Rb (natural)") == GENERIC:
            return self._compute_generic(params)
        return self._compute_species(params)

    # ---- realistic multilevel OD/SAS model ----
    def _compute_species(self, params):
        line = params["line"]
        comps = species.SPECIES[params["species"]]
        T = params["temp_c"] + 273.15
        gt = 2 * np.pi * params["transit_khz"] * 1e3
        power, waist = params["pump_power_mw"], params["waist_mm"]
        buffer_gamma = constants.neon_buffer_broadening(params.get("ne_pressure_torr", 0.0))

        iso_ref = max(comps, key=lambda c: c[1])[0]
        nu_ref = iso_ref.line(line)[1]

        built, omega_all, dopp_fwhm, gamma_max = [], [], 0.0, 0.0
        for iso, weight in comps:
            man = species.build_manifold(iso, line, transit_rate=gt)
            offset = 2 * np.pi * (man.nu0 - nu_ref)
            sigma_v = np.sqrt(constants.KB * T / iso.mass)
            dopp_fwhm = max(dopp_fwhm, np.sqrt(8 * np.log(2)) * man.k_vec * sigma_v)
            N = species.number_density(iso, T) * weight
            Op = species.pump_rabi_from_power(power, waist, man.gamma)
            gamma_eff = man.gamma + species.self_broadened_gamma(iso, N) + buffer_gamma
            gamma_max = max(gamma_max, gamma_eff)
            built.append(dict(man=man, offset=offset, N=N, Op=Op, iso=iso,
                              gamma_eff=gamma_eff))
            omega_all.append(man.omega + offset)
        omega_all = np.concatenate(omega_all)

        margin = max(3.5 * dopp_fwhm, 6.0 * gamma_max)
        scan = np.linspace(omega_all.min() - margin, omega_all.max() + margin,
                           int(params["scan_points"]))

        alpha = np.zeros(scan.size)
        markers = []
        ng_of = {id(b["man"]): len(b["man"].Fg) for b in built}
        for b in built:
            alpha += self._component_alpha(b, scan, T)
            man, off = b["man"], b["offset"]
            ng = ng_of[id(man)]
            for t in range(man.omega.size):
                fg = man.Fg[man.g_idx[t]]
                fe = man.Fe[man.e_idx[t] - ng]
                markers.append((float((man.omega[t] + off) / (2 * np.pi) / 1e9),
                                f"{man.iso.label} {fg:g}→{fe:g}′"))

        return dict(mode="species", scan=scan, alpha_unit=alpha,
                    dopp_fwhm=dopp_fwhm, buffer_gamma=buffer_gamma,
                    gamma_eff_max=gamma_max, markers=markers,
                    species=params["species"], line=line)

    def _component_alpha(self, b, scan, T):
        """Σ_(Fg→Fe) A_t · ĝ_t(δ) for one isotope (1/m, line strength 1)."""
        man, offset, N, Op = b["man"], b["offset"], b["N"], b["Op"]
        iso, k = b["iso"], man.k_vec
        ng = len(man.Fg)
        gamma_eff = b["gamma_eff"]

        v, wt = doppler.velocity_grid(T, mass=iso.mass, dv=1.0, cutoff_sigma=4.0)
        kv = k * v

        # Pump steady state ρ(Δ_eff): scan-independent H, one fine table.
        Hp = species.pump_hamiltonian(man, Op)
        L0 = core.build_liouvillian(Hp, man.atom)
        om = man.omega
        de_lo, de_hi = om.min() - 14 * gamma_eff, om.max() + 14 * gamma_eff
        deff = np.linspace(de_lo, de_hi, int((de_hi - de_lo) / (gamma_eff / 8.0)) + 2)
        pops = _pump_pops(L0, deff, man.atom.S_v, man.n_levels)

        # Homogeneous unit-area lineshape: the weak-probe 2-level absorption is a
        # Lorentzian of FWHM Γ_eff, ∫L̂ dδ = 1. Evaluated analytically at the probe
        # detuning (exact, no truncation); its Doppler sum is the Voigt.
        hwhm = gamma_eff / 2.0

        # Absolute per-line integrated absorption A_t (AutoOD normalisation, ls=1).
        Aline = species.line_integrated_alpha(iso, line=man.line, N=N)

        DS = scan - offset
        deff_grid = DS[:, None] + kv[None, :]                # pump Δ_eff (ns, nv)
        probe_base = DS[:, None] - kv[None, :]               # probe arg base
        levels = set(man.g_idx.tolist()) | set(man.e_idx.tolist())
        pop_at = {lvl: np.interp(deff_grid.ravel(), deff, pops[:, lvl]).reshape(deff_grid.shape)
                  for lvl in levels}

        alpha = np.zeros(scan.size)
        for t in range(om.size):
            g, e = man.g_idx[t], man.e_idx[t]
            fg, fe = man.Fg[g], man.Fe[e - ng]
            A_t = Aline[(fg, fe)]
            w = (pop_at[g] - pop_at[e]) / man.p_ground[g]    # 1 at pump off
            arg = probe_base - om[t]
            Lp = (hwhm / np.pi) / (arg ** 2 + hwhm ** 2)     # unit-area Lorentzian
            alpha += A_t * ((w * Lp) @ wt)
        return alpha

    # ---- generic Γ-unit hole-burning toy (pedagogical) ----
    def _compute_generic(self, params):
        two_lines = params["transitions"].startswith("two")
        if two_lines:
            split = params["splitting"] * GAMMA
            n_exc, offsets = 2, np.array([-split / 2, split / 2])
        else:
            n_exc, offsets = 1, np.array([0.0])

        buffer_gamma = constants.neon_buffer_broadening(params.get("ne_pressure_torr", 0.0))
        gamma_eff = GAMMA + buffer_gamma
        atom = atoms.sas_atom(n_exc, gamma=gamma_eff)
        Op = species.pump_rabi_from_power(params["pump_power_mw"], params["waist_mm"], GAMMA)
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
        half = max(3.5 * dopp_fwhm, off_span + 0.4 * dopp_fwhm, 10 * gamma_eff)
        scan = np.linspace(-half, half, int(params["scan_points"]))

        two_lvl = atoms.two_level(gamma=gamma_eff)
        Hpr = np.zeros((2, 2), dtype=complex)
        Hpr[0, 1] = Hpr[1, 0] = PROBE_RABI * GAMMA / 2
        L0_probe = core.build_liouvillian(Hpr, two_lvl)
        kvmax = float(np.abs(kv).max())
        flo, fhi = scan.min() - kvmax - off_span, scan.max() + kvmax + off_span
        fine = np.linspace(flo, fhi, int((fhi - flo) / (gamma_eff / 20)) + 2)
        rho_pr = core.steady_state_batched(L0_probe, fine, two_lvl.S_v, 2)
        chi_pr = rho_pr[:, 1, 0] / (PROBE_RABI * GAMMA)
        alpha_L, _ = observables.absorption_coefficient(chi_pr, K_VEC, N)

        # The pump H₀ is scan-independent — the scan enters only through the pump
        # Δ_eff = D + k·v — so solve the pump populations ONCE on a fine Δ_eff
        # table and interpolate, instead of re-solving the OBE at every scan
        # point. Same Δ_eff-table trick the realistic species OD/SAS path uses;
        # cuts the pump solves from (scan_points × velocity) down to one table.
        ns, nv = scan.size, v.size
        dlo, dhi = scan.min() + kv.min(), scan.max() + kv.max()
        n_grid = int((dhi - dlo) / (gamma_eff / 8.0)) + 2
        deff_p = np.linspace(dlo, dhi, n_grid)
        pops = _pump_pops(L0_pump, deff_p, atom.S_v, atom.n_levels)   # (n_grid, n)

        dgrid = (scan[:, None] + kv[None, :]).ravel()                 # pump Δ_eff
        rho_gg = np.interp(dgrid, deff_p, pops[:, 0]).reshape(ns, nv)
        probe_base = scan[:, None] - kv[None, :]
        alpha = np.zeros(ns)
        for i, e in enumerate(atom.excited):
            pe = np.interp(dgrid, deff_p, pops[:, e]).reshape(ns, nv)
            aL = np.interp((probe_base - offsets[i]).ravel(),
                           fine, alpha_L).reshape(ns, nv)
            alpha += ((rho_gg - pe) * aL * wt[None, :]).sum(axis=1)

        return dict(mode="generic", scan=scan, alpha_unit=alpha,
                    dopp_fwhm=dopp_fwhm, buffer_gamma=buffer_gamma,
                    gamma_eff=gamma_eff, offsets=offsets, two=(n_exc == 2))

    # =================================================================
    # observables  (dispatch)
    # =================================================================
    def observables(self, raw, params, include_figures=True):
        L = params["cell_mm"] * 1e-3
        alpha = raw["alpha_unit"] * params["line_strength"]
        if raw["mode"] == "species":
            return self._obs_species(
                raw, params, alpha, L, include_figures=include_figures)
        return self._obs_generic(
            raw, params, alpha, L, include_figures=include_figures)

    def _obs_species(self, raw, params, alpha, L, include_figures=True):
        x = raw["scan"] / (2 * np.pi) / 1e9                  # GHz (relative)
        T_trans = observables.transmission(alpha, L)
        OD = observables.optical_density(alpha, L)
        buffer_mhz = raw.get("buffer_gamma", 0.0) / (2 * np.pi) / 1e6
        pump = params["pump_power_mw"]
        pump_on = pump > 0
        regime = "OD (pump off)" if pump <= 0 else f"SAS, P = {pump:.2f} mW"

        fig = None
        figure_views = []
        if include_figures:
            import matplotlib.pyplot as plt

            title = (f"{raw['species']} {raw['line']} — {regime}:  "
                     f"T = {params['temp_c']:.0f} °C, "
                     f"L = {params['cell_mm']:.0f} mm")
            xlabel = "Relative frequency  [GHz]  (ref: line centroid)"
            fig, axT = plt.subplots(1, 1, figsize=(8.5, 4.35))
            axT.plot(x, T_trans, color="#0284C7", lw=1.3)
            axT.set_ylabel("Transmission")
            axT.set_xlabel(xlabel)
            axT.set_title(title)
            fig_od, axA = plt.subplots(1, 1, figsize=(8.5, 4.35))
            axA.plot(x, OD, color="#F43F5E", lw=1.3)
            axA.set_ylabel("Optical density")
            axA.set_xlabel(xlabel)
            axA.set_title(title)
            for gx, _lbl in raw["markers"]:
                for ax in (axT, axA):
                    ax.axvline(gx, color="gray", ls=":", lw=0.5, alpha=0.6)
            fig.tight_layout()
            fig_od.tight_layout()
            figure_views = [
                {"label": "Transmission", "figure": fig},
                {"label": "Optical density", "figure": fig_od},
            ]

        sub_fwhm, sub_at = narrowest_subdoppler(x, T_trans)
        if pump <= 0:                                        # OD limit: no holes burned
            sub_fwhm = float("nan")
        feature_metric = None
        if np.isfinite(sub_fwhm):
            feature_metric = dict(
                label="Narrowest sub-Doppler", value=f"{sub_fwhm*1e3:.1f} MHz",
                help=f"Sharpest Doppler-free feature (near {sub_at:.2f} GHz).")
        elif pump_on:
            feature_metric = dict(
                label="SAS status", value="sub-Doppler unresolved", kind="status",
                help="The pump is on, but no finite sub-Doppler feature width was "
                     "resolved in the displayed spectrum.")
        peak_metric = dict(label="Peak OD", value=f"{np.nanmax(OD):.2f}")
        if pump_on:
            metrics = _lock_readout_metrics(
                x, T_trans, mhz_per_x=1000.0,
                feature_at=sub_at, feature_fwhm=sub_fwhm)
            metrics[0]["tier"] = "hero"
            metrics.append(peak_metric)
            if feature_metric is not None:
                metrics.append(feature_metric)
        else:
            peak_metric["tier"] = "hero"
            metrics = [peak_metric]
        if params.get("ne_pressure_torr", 0.0) != 0.0:
            metrics.append(dict(
                label="Buffer Gas Broadening", value=f"{buffer_mhz:.1f} MHz"))
        rows = "".join(f"| {lbl} | {gx*1e3:.1f} |\n" for gx, lbl in raw["markers"])
        table = ("Hyperfine transitions (Lamb-dip centres); crossovers appear at the "
                 "midpoint of any two sharing a ground state, enhanced/inverted by "
                 "hyperfine pumping.\n\n| Transition | Center [MHz] |\n|---|---|\n" + rows)
        return dict(
            metrics=metrics,
            hero_count=1,
            figure=fig,
            figure_views=figure_views,
            tables=[{"title": "Hyperfine lines", "markdown": table}],
            comparison={
                "axis_index": 0,
                "x_unit": "GHz",
                "raw_x_unit": "Arb. unit",
                "raw_y_unit": "Arb. unit",
                "label": "Experimental CSV",
            },
        )

    def _obs_generic(self, raw, params, alpha, L, include_figures=True):
        x = raw["scan"] / (2 * np.pi) / 1e6                  # MHz
        T_trans = observables.transmission(alpha, L)
        OD = observables.optical_density(alpha, L)
        offs_mhz = raw["offsets"] / (2 * np.pi) / 1e6
        buffer_mhz = raw.get("buffer_gamma", 0.0) / (2 * np.pi) / 1e6

        fig = None
        figure_views = []
        if include_figures:
            import matplotlib.pyplot as plt

            title = (f"Generic SAS: P = {params['pump_power_mw']:.2f} mW, "
                     f"T = {params['temp_c']:.0f} °C")
            xlabel = "Probe detuning  [MHz]"
            fig, axT = plt.subplots(1, 1, figsize=(8.5, 4.35))
            axT.plot(x, T_trans, color="#0284C7", lw=1.6)
            axT.set_ylabel("Transmission")
            axT.set_xlabel(xlabel)
            axT.set_title(title)
            fig_od, axA = plt.subplots(1, 1, figsize=(8.5, 4.35))
            axA.plot(x, OD, color="#F43F5E", lw=1.6)
            axA.set_ylabel("Optical density")
            axA.set_xlabel(xlabel)
            axA.set_title(title)
            for ax in (axT, axA):
                for off in offs_mhz:
                    ax.axvline(off, color="gray", ls=":", lw=0.7)
                if raw["two"]:
                    ax.axvline(0.0, color="green", ls=":", lw=0.7)
            fig.tight_layout()
            fig_od.tight_layout()
            figure_views = [
                {"label": "Transmission", "figure": fig},
                {"label": "Optical density", "figure": fig_od},
            ]

        ic = int(np.argmin(np.abs(x - offs_mhz[0])))
        sub_fwhm = window_fwhm(x, T_trans, ic)
        pump_on = params["pump_power_mw"] > 0
        feature_metric = None
        if pump_on and np.isfinite(sub_fwhm):
            feature_metric = dict(
                label="Sub-Doppler dip FWHM", value=f"{sub_fwhm:.1f} MHz")
        elif pump_on:
            feature_metric = dict(
                label="SAS status", value="sub-Doppler unresolved", kind="status",
                help="The pump is on, but no finite Lamb-dip width was resolved "
                     "around the selected transition.")
        peak_metric = dict(label="Peak OD", value=f"{np.nanmax(OD):.2f}")
        if pump_on:
            metrics = _lock_readout_metrics(
                x, T_trans, mhz_per_x=1.0,
                feature_at=float(x[ic]), feature_fwhm=sub_fwhm)
            metrics[0]["tier"] = "hero"
            metrics.append(peak_metric)
            if feature_metric is not None:
                metrics.append(feature_metric)
        else:
            peak_metric["tier"] = "hero"
            metrics = [peak_metric]
        if params.get("ne_pressure_torr", 0.0) != 0.0:
            metrics.append(dict(
                label="Buffer Gas Broadening", value=f"{buffer_mhz:.1f} MHz"))
        note = ("Two transitions: Lamb dips at ±splitting/2 and a **crossover** dip "
                "at the midpoint (green)." if raw["two"]
                else "Single transition: one Lamb dip at line centre.")
        return dict(
            metrics=metrics,
            hero_count=1,
            figure=fig,
            figure_views=figure_views,
            tables=[{"title": "Notes", "markdown": note}],
            comparison={
                "axis_index": 0,
                "x_unit": "MHz",
                "raw_x_unit": "Arb. unit",
                "raw_y_unit": "Arb. unit",
                "label": "Experimental CSV",
            },
        )


# =====================================================================
# helpers
# =====================================================================
def _pump_pops(L0, deff_axis, S_v, n, chunk=1500):
    """Diagonal populations ρ_ii(Δ_eff) on a fine axis, in memory-safe chunks."""
    pops = np.empty((deff_axis.size, n))
    for s in range(0, deff_axis.size, chunk):
        sl = slice(s, s + chunk)
        rho = core.steady_state_batched(L0, deff_axis[sl], S_v, n)
        pops[sl] = np.einsum("vii->vi", rho).real
    return pops


def _lock_readout_metrics(x, T_trans, mhz_per_x,
                          feature_at=None, feature_fwhm=None):
    """Finite-slope lock-point proxy from the displayed transmission curve.

    When a sub-Doppler feature has been detected, search only within one FWHM
    of its centre.  This prevents a steeper Doppler-envelope flank (or another
    hyperfine manifold) from being reported as the pump-on laser-lock point.
    Failed-feature-detection paths retain the full-spectrum search.
    """
    if len(x) < 2:
        return []
    slope_per_x = np.gradient(T_trans, x)
    candidates = np.isfinite(slope_per_x)
    local_search = (
        feature_at is not None
        and feature_fwhm is not None
        and np.isfinite(feature_at)
        and np.isfinite(feature_fwhm)
        and float(feature_fwhm) > 0.0
    )
    if local_search:
        local = np.abs(np.asarray(x) - float(feature_at)) <= float(feature_fwhm)
        if np.any(candidates & local):
            candidates &= local
        else:
            local_search = False
    if not candidates.any():
        return []
    indices = np.flatnonzero(candidates)
    i = int(indices[np.argmax(np.abs(slope_per_x[indices]))])
    slope_per_mhz = abs(float(slope_per_x[i])) / float(mhz_per_x)
    detuning_mhz = float(x[i]) * float(mhz_per_x)
    slope_help = (
        "Largest |dT/dΔ| within one FWHM of the detected sub-Doppler "
        "feature; proxy for a laser-lock discriminator."
        if local_search else
        "Largest |dT/dΔ| on the displayed spectrum; proxy for a laser-lock "
        "discriminator."
    )
    return [
        dict(label="Lock Slope", value=f"{slope_per_mhz:.4f} /MHz",
             help=slope_help),
        dict(label="Lock Detuning", value=f"{detuning_mhz:+.1f} MHz",
             help="Detuning where the lock-slope proxy is largest."),
    ]


