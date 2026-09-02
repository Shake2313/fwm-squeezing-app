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
from ..lineshape import subdoppler_feature
from .base import ParamSpec, Scheme

PROBE_RABI = 1e-3                       # weak probe, in units of Γ
GAMMA_MHZ = GAMMA / (2 * np.pi) / 1e6
GENERIC = "Generic (Γ units)"
PARAFFIN_REFERENCE_T1_S = 25.1e-3   # 87Rb population T1 at 300 K (Bandi et al.)
_trapz = getattr(np, "trapezoid", getattr(np, "trapz", None))   # numpy ≥2.0 rename


def _temperatures(params):
    """(density T, Doppler T) in kelvin, plus whether the two are tied.

    ``temp_c`` is the **vapor-density** temperature (cold spot / reservoir): it
    fixes the number density, hence the absorption scale and the density-
    dependent self-broadening. The **Doppler-width** temperature sets the
    thermal velocity distribution (Doppler width and velocity classes) and
    follows the cell temperature unless ``constrain_doppler_temp`` is released.
    Buffer-gas broadening is pressure-only here (no temperature enters), and the
    transit rate stays a user knob; see `collisional-coefficient-provenance-and-
    pressure-shift` in docs/checklist.json.

    A params dict from before these knobs existed stays tied, so legacy calls
    reproduce their previous spectra exactly.
    """
    density_c = float(params["temp_c"])
    tied = bool(params.get("constrain_doppler_temp", True))
    doppler_c = density_c if tied else float(
        params.get("doppler_temp_c", density_c))
    return density_c + 273.15, doppler_c + 273.15, tied


class SASScheme(Scheme):
    name = "sas"
    cluster = "A — Absorption"
    title = "Absorption spectroscopy (OD / SAS)"
    cache_version = "4"            # + dispersive susceptibility, split temperatures
    defaults_version = "3"         # new param schema: reseed sidebar defaults
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
            ParamSpec("temp_c", "Cell temperature", "Cell & beams", 40.0, 20.0, 200.0,
                      1.0, "°C",
                      help="Vapor-density temperature: it fixes the number density "
                      "(absorption scale) and the density-dependent self-broadening. "
                      "It also sets the Doppler width unless the Advanced "
                      "Doppler-width temperature is released."),
            ParamSpec("cell_mm", "Cell length", "Cell & beams", 75.0, 0.5, 200.0, 0.5, "mm",
                      recompute=False),
            ParamSpec(
                "constrain_doppler_temp", "Tie Doppler width to cell temperature",
                "Cell & beams", True, advanced=True, control="checkbox",
                advanced_group="Cell temperatures",
                help="On (default): one temperature sets both the vapor density and "
                     "the Doppler width. Off: the cell temperature keeps setting the "
                     "density (cold spot / reservoir) while the Doppler-width "
                     "temperature below sets the thermal velocity distribution.",
            ),
            ParamSpec("doppler_temp_c", "Doppler-width temperature", "Cell & beams",
                      40.0, 20.0, 300.0, 1.0, "°C", advanced=True,
                      advanced_group="Cell temperatures",
                      visible_if={"constrain_doppler_temp": False},
                      help="Temperature of the atoms crossing the beam: it sets the "
                      "thermal velocity distribution (Doppler width, velocity "
                      "classes) only. The vapor density still follows the cell "
                      "temperature. Ignored while the tie above is on."),
            ParamSpec(
                "paraffin_coated", "Paraffin-coated cell", "Cell & beams", False,
                advanced=True, control="checkbox",
                visible_if={"species": tuple(species.SPECIES_ORDER)},
                help="Preserve ground-hyperfine population between velocity-randomized "
                     "beam passages using a nominal coated-cell reservoir. This changes "
                     "atomic pumping only; it does not apply a coating transmission loss.",
            ),
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
            "The Advanced **Paraffin-coated cell** switch retains the two ground-"
            "hyperfine populations between beam passages while rethermalizing velocity. "
            "It uses the existing transit rate for beam exchange and a fixed reference "
            "population lifetime T1 = 25.1 ms, measured for a paraffin-coated 87Rb cell "
            "at 300 K. This is a quasi-static population-memory model, not a calibrated "
            "cell-geometry model; it preserves no Zeeman coherence and adds no window "
            "or coating throughput loss.\n\n"
            "**Two temperatures.** The **cell temperature** is the vapor-density "
            "temperature: it fixes the number density (absorption scale) and the "
            "density-dependent self-broadening. The **Doppler-width temperature** "
            "sets the thermal velocity distribution only. They are tied by default; "
            "release the Advanced tie to model a cold spot below the beam-path "
            "temperature. Buffer-gas broadening here is pressure-only (no "
            "temperature dependence) and the transit rate stays a separate knob.\n\n"
            "**Dispersion.** The absorption and its Kramers-Kronig partner come from "
            "one complex line profile, so the **Dispersion** view shows the refractive "
            "index n = 1 + Re χ/2 (dilute vapor) built from the *same* velocity-"
            "resolved population difference as the absorption — under pumping it is "
            "the saturated medium's dispersion, not a weak-probe curve pasted next to "
            "a saturated spectrum. The reported **Peak phase shift** is a single-pass "
            "phase φ = k(n−1)L, not a group delay.\n\n"
            "The reported **Gaussian Doppler FWHM** is calculated from the thermal "
            "velocity distribution for one optical line; it is not a measured "
            "Voigt or multi-line envelope FWHM. Sub-Doppler widths use interpolated "
            "half-height edges and report samples per FWHM plus scan-edge clearance. "
            "A local lock slope is shown only when that resolution status is "
            "resolved.\n\n"
            "Atomic data (hyperfine A/B, line centres, masses, linewidths) from the "
            "Steck D-line data sheets; Wigner-6j/3j line strengths in the AutoOD "
            "convention. Rb densities use the CRC vapor pressure (AutoOD), Cs the "
            "Steck fit.\n\n"
            "**References**\n"
            "- D. A. Smith & I. G. Hughes, *Am. J. Phys.* **72**, 631 (2004).\n"
            "- D. W. Preston, *Am. J. Phys.* **64**, 1432 (1996).\n"
            "- T. Bandi, C. Affolderbach & G. Mileti, *J. Appl. Phys.* **111**, "
            "124906 (2012), https://doi.org/10.1063/1.4729925.\n"
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
        T_density, T_doppler, tied = _temperatures(params)
        gt = 2 * np.pi * params["transit_khz"] * 1e3
        power, waist = params["pump_power_mw"], params["waist_mm"]
        paraffin_coated = bool(params.get("paraffin_coated", False))
        # With no pump, coating memory cannot change a thermal OD spectrum. Keep
        # that limit on the legacy path exactly (and avoid an unnecessary solve).
        paraffin_memory = paraffin_coated and power > 0.0
        buffer_gamma = constants.neon_buffer_broadening(params.get("ne_pressure_torr", 0.0))

        iso_ref = max(comps, key=lambda c: c[1])[0]
        nu_ref = iso_ref.line(line)[1]

        built, omega_all, dopp_fwhm, gamma_max = [], [], 0.0, 0.0
        for iso, weight in comps:
            man = species.build_manifold(
                iso, line, transit_rate=0.0 if paraffin_memory else gt)
            offset = 2 * np.pi * (man.nu0 - nu_ref)
            sigma_v = np.sqrt(constants.KB * T_doppler / iso.mass)
            dopp_fwhm = max(dopp_fwhm, np.sqrt(8 * np.log(2)) * man.k_vec * sigma_v)
            N = species.number_density(iso, T_density) * weight
            Op = species.pump_rabi_from_power(power, waist, man.gamma)
            gamma_eff = man.gamma + species.self_broadened_gamma(iso, N) + buffer_gamma
            gamma_max = max(gamma_max, gamma_eff)
            built.append(dict(man=man, offset=offset, N=N, Op=Op, iso=iso,
                              gamma_eff=gamma_eff, transit_rate=gt,
                              paraffin_memory=paraffin_memory))
            omega_all.append(man.omega + offset)
        omega_all = np.concatenate(omega_all)

        margin = max(3.5 * dopp_fwhm, 6.0 * gamma_max)
        scan = np.linspace(omega_all.min() - margin, omega_all.max() + margin,
                           int(params["scan_points"]))

        alpha = np.zeros(scan.size)
        chi_real = np.zeros(scan.size)
        markers = []
        ng_of = {id(b["man"]): len(b["man"].Fg) for b in built}
        for b in built:
            component_alpha, component_chi_real = self._component_response(
                b, scan, T_doppler)
            alpha += component_alpha
            chi_real += component_chi_real
            man, off = b["man"], b["offset"]
            ng = ng_of[id(man)]
            for t in range(man.omega.size):
                fg = man.Fg[man.g_idx[t]]
                fe = man.Fe[man.e_idx[t] - ng]
                markers.append((float((man.omega[t] + off) / (2 * np.pi) / 1e9),
                                f"{man.iso.label} {fg:g}→{fe:g}′"))

        # Each component's own k converts its A_t into a susceptibility; the
        # single k reported here only turns that susceptibility into a phase, so
        # the reference isotope's value is enough (isotope shifts are ~1e-7 of k).
        k_ref = float(next(b["man"].k_vec for b in built if b["iso"] is iso_ref))
        return dict(mode="species", scan=scan, alpha_unit=alpha,
                    chi_real_unit=chi_real, k_vec=k_ref,
                    dopp_fwhm=dopp_fwhm, buffer_gamma=buffer_gamma,
                    gamma_eff_max=gamma_max, markers=markers,
                    species=params["species"], line=line,
                    density_temp_c=T_density - 273.15,
                    doppler_temp_c=T_doppler - 273.15,
                    doppler_temp_tied=tied,
                    paraffin_coated=paraffin_coated,
                    paraffin_t1_s=(PARAFFIN_REFERENCE_T1_S
                                   if paraffin_coated else None))

    def _component_response(self, b, scan, T):
        """(α, Re χ) for one isotope at line strength 1.

        Both quadratures come from the same complex line profile

            χ_t(δ) = (A_t / k) · i / (π (Γ_eff/2 − i·δ_t)),

        whose imaginary part is the unit-area absorption profile ĝ_t used for α
        (∫ĝ dδ = 1) and whose real part is its Kramers-Kronig partner. Sign and
        normalisation follow `observables.chi_phys` (α = k·Im χ, absorptive for a
        passive transition), so n = 1 + Re χ/2 rises below resonance. The same
        velocity-resolved population difference and the same Doppler quadrature
        weight both parts — a saturated absorption spectrum is never paired with
        a weak-probe dispersion.
        """
        man, offset, N, Op = b["man"], b["offset"], b["N"], b["Op"]
        iso, k = b["iso"], man.k_vec
        ng = len(man.Fg)
        gamma_eff = b["gamma_eff"]

        v, wt = doppler.velocity_grid(T, mass=iso.mass, dv=1.0, cutoff_sigma=4.0)
        kv = k * v

        # Pump steady state ρ(Δ_eff): scan-independent H, one fine table.
        Hp = species.pump_hamiltonian(man, Op)
        om = man.omega
        DS = scan - offset
        covered_detunings = None
        if b.get("paraffin_memory", False):
            covered_detunings = (
                float(DS.min() + kv.min()), float(DS.max() + kv.max()))
        deff = _pump_detuning_axis(
            om, gamma_eff, covered_detunings=covered_detunings)
        if b.get("paraffin_memory", False):
            basis_pops = _basis_reset_pump_pops(
                man, Hp, deff, b["transit_rate"])
        else:
            L0 = core.build_liouvillian(Hp, man.atom)
            pops = _pump_pops(L0, deff, man.atom.S_v, man.n_levels)

        # Homogeneous unit-area lineshape: the weak-probe 2-level absorption is a
        # Lorentzian of FWHM Γ_eff, ∫L̂ dδ = 1. Evaluated analytically at the probe
        # detuning (exact, no truncation); its Doppler sum is the Voigt.
        hwhm = gamma_eff / 2.0

        # Absolute per-line integrated absorption A_t (AutoOD normalisation, ls=1).
        Aline = species.line_integrated_alpha(iso, line=man.line, N=N)

        deff_grid = DS[:, None] + kv[None, :]                # pump Δ_eff (ns, nv)
        probe_base = DS[:, None] - kv[None, :]               # probe arg base
        levels = set(man.g_idx.tolist()) | set(man.e_idx.tolist())
        if b.get("paraffin_memory", False):
            reservoir = _coated_ground_populations(
                man, basis_pops, deff, deff_grid, wt,
                cycle_rate=b["transit_rate"])
            pop_at = {}
            flat_deff = deff_grid.ravel()
            for lvl in levels:
                local = np.zeros_like(deff_grid)
                for source in range(ng):
                    conditioned = np.interp(
                        flat_deff, deff, basis_pops[source, :, lvl]
                    ).reshape(deff_grid.shape)
                    local += reservoir[:, source, None] * conditioned
                pop_at[lvl] = local
        else:
            pop_at = {
                lvl: np.interp(deff_grid.ravel(), deff, pops[:, lvl]).reshape(
                    deff_grid.shape)
                for lvl in levels
            }

        alpha = np.zeros(scan.size)
        chi_real = np.zeros(scan.size)
        for t in range(om.size):
            g, e = man.g_idx[t], man.e_idx[t]
            fg, fe = man.Fg[g], man.Fe[e - ng]
            A_t = Aline[(fg, fe)]
            w = (pop_at[g] - pop_at[e]) / man.p_ground[g]    # 1 at pump off
            arg = probe_base - om[t]
            denom = arg ** 2 + hwhm ** 2
            Lp = (hwhm / np.pi) / denom                      # unit-area Lorentzian
            Dp = (arg / np.pi) / denom                       # KK dispersive partner
            alpha += A_t * _velocity_correlated_average(w, Lp, wt)
            chi_real -= (A_t / k) * _velocity_correlated_average(w, Dp, wt)
        return alpha, chi_real

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

        T_density, T_doppler, tied = _temperatures(params)
        sigma_v = np.sqrt(constants.KB * T_doppler / constants.MASS_85RB)
        dopp_fwhm = np.sqrt(8 * np.log(2)) * K_VEC * sigma_v
        v, wt = doppler.velocity_grid(T_doppler, dv=3.0, cutoff_sigma=3.5)
        kv = K_VEC * v
        N = atoms.rb85_density(T_density)

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
        # Both quadratures of one weak-probe susceptibility: α = k·Im χ is the
        # absorption table, Re χ its dispersive partner on the same fine axis.
        alpha_L, chi_L = observables.absorption_coefficient(chi_pr, K_VEC, N)
        chi_real_L = np.real(chi_L)

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
        chi_real = np.zeros(ns)
        for i, e in enumerate(atom.excited):
            pe = np.interp(dgrid, deff_p, pops[:, e]).reshape(ns, nv)
            probe_arg = (probe_base - offsets[i]).ravel()
            aL = np.interp(probe_arg, fine, alpha_L).reshape(ns, nv)
            xL = np.interp(probe_arg, fine, chi_real_L).reshape(ns, nv)
            alpha += ((rho_gg - pe) * aL * wt[None, :]).sum(axis=1)
            chi_real += ((rho_gg - pe) * xL * wt[None, :]).sum(axis=1)

        return dict(mode="generic", scan=scan, alpha_unit=alpha,
                    chi_real_unit=chi_real, k_vec=float(K_VEC),
                    dopp_fwhm=dopp_fwhm, buffer_gamma=buffer_gamma,
                    density_temp_c=T_density - 273.15,
                    doppler_temp_c=T_doppler - 273.15,
                    doppler_temp_tied=tied,
                    gamma_eff=gamma_eff, offsets=offsets, two=(n_exc == 2))

    # =================================================================
    # observables  (dispatch)
    # =================================================================
    def observables(self, raw, params, include_figures=True):
        L = params["cell_mm"] * 1e-3
        # α and Re χ are the two quadratures of one susceptibility, so the same
        # line-strength calibration scales both.
        alpha = raw["alpha_unit"] * params["line_strength"]
        chi_real_unit = raw.get("chi_real_unit")
        chi_real = (None if chi_real_unit is None
                    else chi_real_unit * params["line_strength"])
        if raw["mode"] == "species":
            return self._obs_species(
                raw, params, alpha, chi_real, L, include_figures=include_figures)
        return self._obs_generic(
            raw, params, alpha, chi_real, L, include_figures=include_figures)

    def _obs_species(self, raw, params, alpha, chi_real, L, include_figures=True):
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
                     f"{_temperature_title(raw, params)}, "
                     f"L = {params['cell_mm']:.0f} mm")
            xlabel = "Relative frequency  [GHz]  (ref: line centroid)"
            marker_x = [gx for gx, _lbl in raw["markers"]]
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
            for gx in marker_x:
                for ax in (axT, axA):
                    ax.axvline(gx, color="gray", ls=":", lw=0.5, alpha=0.6)
            fig.tight_layout()
            fig_od.tight_layout()
            figure_views = [
                {"label": "Transmission", "figure": fig},
                {"label": "Optical density", "figure": fig_od},
            ]
            if chi_real is not None:
                figure_views.append({
                    "label": "Dispersion",
                    "figure": _dispersion_figure(
                        x, chi_real, xlabel, title, markers=marker_x),
                })

        feature = subdoppler_feature(x, T_trans) if pump_on else None
        peak_metric = dict(label="Peak OD", value=f"{np.nanmax(OD):.2f}")
        broad_metric = _gaussian_doppler_metric(raw["dopp_fwhm"], raw)
        if pump_on:
            status_metric, feature_metrics = _subdoppler_metrics(
                feature, mhz_per_x=1000.0)
            if feature.resolved:
                metrics = _lock_readout_metrics(
                    x, T_trans, mhz_per_x=1000.0,
                    feature_at=feature.center, feature_fwhm=feature.fwhm)
                metrics[0]["tier"] = "hero"
                metrics.extend(feature_metrics)
                metrics.append(status_metric)
            else:
                status_metric["tier"] = "hero"
                metrics = [status_metric, *feature_metrics]
            metrics.extend([peak_metric, broad_metric])
        else:
            peak_metric["tier"] = "hero"
            metrics = [peak_metric, broad_metric]
        metrics.extend(_dispersion_metrics(chi_real, raw.get("k_vec"), L))
        metrics.extend(_temperature_metrics(raw))
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

    def _obs_generic(self, raw, params, alpha, chi_real, L, include_figures=True):
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
                     f"{_temperature_title(raw, params)}")
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
            if chi_real is not None:
                figure_views.append({
                    "label": "Dispersion",
                    "figure": _dispersion_figure(
                        x, chi_real, xlabel, title,
                        markers=list(offs_mhz) + ([0.0] if raw["two"] else [])),
                })

        pump_on = params["pump_power_mw"] > 0
        gamma_mhz = raw["gamma_eff"] / (2 * np.pi) / 1e6
        selected_line = float(offs_mhz[0])
        feature_window = (
            selected_line - 10.0 * gamma_mhz,
            selected_line + 10.0 * gamma_mhz,
        )
        feature = (
            subdoppler_feature(x, T_trans, search_window=feature_window)
            if pump_on else None
        )
        peak_metric = dict(label="Peak OD", value=f"{np.nanmax(OD):.2f}")
        broad_metric = _gaussian_doppler_metric(raw["dopp_fwhm"], raw)
        if pump_on:
            status_metric, feature_metrics = _subdoppler_metrics(
                feature, mhz_per_x=1.0)
            if feature.resolved:
                metrics = _lock_readout_metrics(
                    x, T_trans, mhz_per_x=1.0,
                    feature_at=feature.center, feature_fwhm=feature.fwhm)
                metrics[0]["tier"] = "hero"
                metrics.extend(feature_metrics)
                metrics.append(status_metric)
            else:
                status_metric["tier"] = "hero"
                metrics = [status_metric, *feature_metrics]
            metrics.extend([peak_metric, broad_metric])
        else:
            peak_metric["tier"] = "hero"
            metrics = [peak_metric, broad_metric]
        metrics.extend(_dispersion_metrics(chi_real, raw.get("k_vec"), L))
        metrics.extend(_temperature_metrics(raw))
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
def _pump_detuning_axis(omega, gamma_eff, covered_detunings=None):
    """Fine pump-state axis, optionally covering every requested detuning.

    The legacy path uses the historical ±14Γ extent. A long-lived coated
    reservoir must instead solve (not edge-clamp) the full scan+velocity domain,
    because tiny one-pass errors accumulate over many returns.
    """
    omega = np.asarray(omega, dtype=float)
    gamma_eff = float(gamma_eff)
    de_lo = float(omega.min() - 14.0 * gamma_eff)
    de_hi = float(omega.max() + 14.0 * gamma_eff)
    if covered_detunings is not None:
        covered_lo, covered_hi = map(float, covered_detunings)
        de_lo = min(de_lo, covered_lo)
        de_hi = max(de_hi, covered_hi)
    count = int((de_hi - de_lo) / (gamma_eff / 8.0)) + 2
    return np.linspace(de_lo, de_hi, count)


def _velocity_correlated_average(prepared_population, probe_profile, weights):
    """Average population×probe response in the same velocity class."""
    return (prepared_population * probe_profile) @ weights


def _basis_reset_atom(atom, transit_rate, source_ground):
    """Clone a species manifold atom with transit reload into one ground F.

    The basis-reset steady states span the response to any incoherent incoming
    ground-hyperfine distribution. The uncoated path does not call this helper.
    """
    reset = tuple(
        (state, source_ground, float(transit_rate))
        for state in range(atom.n_levels)
    )
    atom_kwargs = dict(
        name=f"{atom.name}_reset_g{source_ground}",
        n_levels=atom.n_levels,
        labels=atom.labels,
        ground=atom.ground,
        excited=atom.excited,
        decay=atom.decay + reset,
        dephasing=atom.dephasing,
        doppler_levels=atom.doppler_levels,
        doppler_ratios=atom.doppler_ratios,
        emission_ops=atom.emission_ops,
    )
    if hasattr(atom, "collapse_ops"):
        atom_kwargs["collapse_ops"] = atom.collapse_ops
    return atoms.AtomModel(**atom_kwargs)


def _basis_reset_pump_pops(man, pump_hamiltonian, deff_axis, transit_rate):
    """Pump populations conditioned on each incoming ground hyperfine state."""
    conditioned = []
    for source_ground in man.atom.ground:
        atom = _basis_reset_atom(man.atom, transit_rate, source_ground)
        L0 = core.build_liouvillian(pump_hamiltonian, atom)
        conditioned.append(
            _pump_pops(L0, deff_axis, atom.S_v, atom.n_levels)
        )
    return np.asarray(conditioned)


def _ground_decay_branching(man):
    """Return B[g,e], including only natural excited-state decay channels."""
    grounds = tuple(man.atom.ground)
    excited = tuple(man.atom.excited)
    g_pos = {level: i for i, level in enumerate(grounds)}
    e_pos = {level: i for i, level in enumerate(excited)}
    branching = np.zeros((len(grounds), len(excited)))
    for source, target, rate in man.atom.decay:
        if source in e_pos and target in g_pos:
            branching[g_pos[target], e_pos[source]] += float(rate)
    totals = branching.sum(axis=0)
    if np.any(totals <= 0.0):
        raise ValueError("every excited hyperfine state needs a ground decay path")
    return branching / totals[None, :]


def _coated_ground_transfer(man, basis_pops, deff_axis, deff_grid, weights):
    """Velocity-averaged ground-F transfer matrix for one pump scan.

    ``transfer[scan, destination, source]`` includes spontaneous decay of any
    excited population after the atom exits the illuminated region. Velocity is
    averaged only for this dark reservoir map; the absorption calculation keeps
    each pump/probe velocity correlation intact.
    """
    grounds = np.asarray(man.atom.ground, dtype=int)
    excited = np.asarray(man.atom.excited, dtype=int)
    branching = _ground_decay_branching(man)
    exit_pops = basis_pops[:, :, grounds].copy()
    exit_pops += np.einsum(
        "sxe,ge->sxg", basis_pops[:, :, excited], branching)

    n_scan = deff_grid.shape[0]
    n_ground = grounds.size
    transfer = np.empty((n_scan, n_ground, n_ground))
    flat_deff = deff_grid.ravel()
    for source in range(n_ground):
        for destination in range(n_ground):
            local = np.interp(
                flat_deff, deff_axis, exit_pops[source, :, destination]
            ).reshape(deff_grid.shape)
            transfer[:, destination, source] = local @ weights

    # Interpolation/linear-solve roundoff can make a probability about -1e-16.
    transfer = np.maximum(transfer, 0.0)
    column_sum = transfer.sum(axis=1, keepdims=True)
    if np.any(column_sum <= 0.0):
        raise ValueError("coated-cell ground transfer lost probability")
    return transfer / column_sum


def _stationary_ground_populations(
        transfer, thermal_populations, cycle_rate,
        wall_t1=PARAFFIN_REFERENCE_T1_S):
    """Solve the quasi-static coated-cell ground-population reservoir.

    The existing transit rate closes the one-zone model as both bright-region
    exit and return cadence. Wall relaxation drives the reservoir toward the
    thermal hyperfine populations at 1/T1 (not 2π/T1).
    """
    transfer = np.asarray(transfer, dtype=float)
    thermal = np.asarray(thermal_populations, dtype=float)
    n_ground = thermal.size
    identity = np.eye(n_ground)
    wall_rate = 1.0 / float(wall_t1)
    A = (float(cycle_rate) * (identity[None, :, :] - transfer)
         + wall_rate * identity[None, :, :])
    rhs = np.broadcast_to(wall_rate * thermal, transfer.shape[:1] + (n_ground,))
    populations = np.linalg.solve(A, rhs[..., None])[..., 0]
    populations = np.maximum(populations, 0.0)
    totals = populations.sum(axis=1, keepdims=True)
    if np.any(totals <= 0.0):
        raise ValueError("coated-cell reservoir has no ground population")
    return populations / totals


def _coated_ground_populations(
        man, basis_pops, deff_axis, deff_grid, weights, cycle_rate):
    transfer = _coated_ground_transfer(
        man, basis_pops, deff_axis, deff_grid, weights)
    return _stationary_ground_populations(
        transfer, man.p_ground, cycle_rate=cycle_rate)


def _pump_pops(L0, deff_axis, S_v, n, chunk=1500):
    """Diagonal populations ρ_ii(Δ_eff) on a fine axis, in memory-safe chunks."""
    pops = np.empty((deff_axis.size, n))
    for s in range(0, deff_axis.size, chunk):
        sl = slice(s, s + chunk)
        rho = core.steady_state_batched(L0, deff_axis[sl], S_v, n)
        pops[sl] = np.einsum("vii->vi", rho).real
    return pops


def _dispersion_figure(x, chi_real, xlabel, title, markers=()):
    """Refractive index over the scan (ASCII axis strings — mathtext lock)."""
    import matplotlib.pyplot as plt

    index_ppm = (observables.refractive_index(chi_real) - 1.0) * 1e6
    fig, ax = plt.subplots(1, 1, figsize=(8.5, 4.35))
    ax.plot(x, index_ppm, color="#7C3AED", lw=1.3)
    ax.axhline(0.0, color="gray", ls="--", lw=0.6, alpha=0.7)
    ax.set_ylabel("Refractive index  n - 1  [ppm]")
    ax.set_xlabel(xlabel)
    ax.set_title(title)
    for gx in markers:
        ax.axvline(gx, color="gray", ls=":", lw=0.5, alpha=0.6)
    fig.tight_layout()
    return fig


def _dispersion_metrics(chi_real, k_vec, L):
    """Single-pass phase readout from the dispersive quadrature (or nothing)."""
    if chi_real is None or k_vec is None:
        return []
    phase = observables.single_pass_phase(chi_real, float(k_vec), float(L))
    peak = float(np.nanmax(np.abs(phase)))
    value = f"{peak:.2f} rad" if peak >= 1.0 else f"{peak * 1e3:.1f} mrad"
    return [dict(
        label="Peak phase shift",
        value=value,
        help="Largest single-pass |φ| = k·(n−1)·L from the dispersive (real) "
             "part of the same susceptibility that gives the absorption, with "
             "n = 1 + Re χ/2. A phase quantity — not a group delay. It is the "
             "medium's phase alone: a deeply absorbed spectral region can carry "
             "a large phase with almost no light left to measure it.",
    )]


def _temperature_metrics(raw):
    """Report the Doppler temperature only when it is released from the cell."""
    if raw.get("doppler_temp_tied", True):
        return []
    return [dict(
        label="Doppler-width temperature",
        value=f"{float(raw['doppler_temp_c']):.0f} °C",
        help="Released from the cell temperature: it sets the thermal velocity "
             "distribution only. The vapor density (absorption scale) and the "
             f"self-broadening still follow the cell temperature "
             f"{float(raw['density_temp_c']):.0f} °C.",
    )]


def _temperature_title(raw, params):
    """Plot-title fragment naming whichever temperatures are in play."""
    density_c = float(raw.get("density_temp_c", params["temp_c"]))
    if raw.get("doppler_temp_tied", True):
        return f"T = {density_c:.0f} °C"
    return (f"T_density = {density_c:.0f} °C, "
            f"T_Doppler = {float(raw['doppler_temp_c']):.0f} °C")


def _gaussian_doppler_metric(dopp_fwhm, raw=None):
    width_mhz = float(dopp_fwhm) / (2 * np.pi) / 1e6
    source = "the cell temperature"
    if raw is not None and not raw.get("doppler_temp_tied", True):
        source = ("the released Doppler-width temperature "
                  f"({float(raw['doppler_temp_c']):.0f} °C)")
    return dict(
        label="Gaussian Doppler FWHM",
        value=f"{width_mhz:.1f} MHz",
        help=f"Calculated thermal Gaussian FWHM for one optical line, from "
             f"{source}. It is not a measured Voigt or multi-line envelope FWHM.",
    )


def _subdoppler_metrics(feature, mhz_per_x):
    """User-facing resolution ledger for a :class:`SubdopplerFeature`."""
    status = dict(
        label="SAS resolution",
        value=feature.status,
        kind="status",
        help=feature.reason,
    )
    if not feature.detected:
        return status, []

    scale = float(mhz_per_x)
    metrics = [
        dict(
            label="Sub-Doppler FWHM",
            value=f"{feature.fwhm * scale:.2f} MHz",
            help=f"Interpolated residual half-height width near "
                 f"{feature.center * scale:+.1f} MHz; trust as a linewidth only "
                 "when SAS resolution is resolved.",
        ),
        dict(
            label="Half-height edges",
            value=(f"{feature.left_half_height * scale:+.2f} to "
                   f"{feature.right_half_height * scale:+.2f} MHz"),
            help="Linear interpolation between the samples bracketing each "
                 "half-height crossing.",
        ),
        dict(
            label="Samples / FWHM",
            value=f"{feature.samples_per_fwhm:.1f}",
            help="Interpolated FWHM divided by the local median sample spacing; "
                 "at least 6 is required for resolved status.",
        ),
        dict(
            label="Scan-edge distance",
            value=f"{feature.scan_edge_distance * scale:.1f} MHz",
            help="Nearest distance from either interpolated half-height edge to "
                 "the displayed scan boundary; at least one FWHM is required.",
        ),
    ]
    return status, metrics


def _lock_readout_metrics(x, T_trans, mhz_per_x,
                          feature_at=None, feature_fwhm=None,
                          search_window=None):
    """Local finite-slope lock proxy for a *resolved* sub-Doppler feature.

    There is deliberately no full-spectrum fallback: a Doppler-envelope flank
    is not a sub-Doppler laser-lock discriminator.
    """
    if len(x) < 2:
        return []
    slope_per_x = np.gradient(T_trans, x)
    candidates = np.isfinite(slope_per_x)
    if search_window is None:
        valid_feature = (
            feature_at is not None
            and feature_fwhm is not None
            and np.isfinite(feature_at)
            and np.isfinite(feature_fwhm)
            and float(feature_fwhm) > 0.0
        )
        if not valid_feature:
            return []
        lo = float(feature_at) - float(feature_fwhm)
        hi = float(feature_at) + float(feature_fwhm)
    else:
        lo, hi = map(float, search_window)
    if not np.isfinite(lo) or not np.isfinite(hi) or lo >= hi:
        return []
    local = (np.asarray(x) >= lo) & (np.asarray(x) <= hi)
    if np.any(candidates & local):
        candidates &= local
    else:
        return []
    if not candidates.any():
        return []
    indices = np.flatnonzero(candidates)
    i = int(indices[np.argmax(np.abs(slope_per_x[indices]))])
    slope_per_mhz = abs(float(slope_per_x[i])) / float(mhz_per_x)
    detuning_mhz = float(x[i]) * float(mhz_per_x)
    return [
        dict(label="Lock Slope", value=f"{slope_per_mhz:.4f} /MHz",
             help="Largest |dT/dΔ| within the resolved sub-Doppler feature "
                  "window; proxy for a laser-lock discriminator."),
        dict(label="Lock Detuning", value=f"{detuning_mhz:+.1f} MHz",
             help="Detuning where the lock-slope proxy is largest."),
    ]


