"""
Cluster B — saturated absorption spectroscopy (SAS).

A strong pump and a weak probe from the same laser counter-propagate through the
cell, so a moving atom sees them with opposite Doppler shifts:
    pump  atom-frame detuning = Δ + k·v
    probe atom-frame detuning = Δ − k·v        (Δ = laser − line centre, the scan)
The pump saturates / optically pumps the velocity class it is resonant with; the
weak probe samples those prepared populations. Where pump and probe address the
*same* atoms (v such that they hit two transitions sharing a level) a sub-Doppler
feature appears: a Lamb dip at each transition and a crossover at every midpoint.

Realistic multilevel model (the default)
----------------------------------------
`gabes.species` builds the full hyperfine manifold for a chosen isotope/line
(⁸⁵Rb, ⁸⁷Rb, ¹³³Cs · D1/D2): all ground F and excited F′ as lumped levels, with
CG-branched spontaneous emission *and* a transit-time relaxation toward the
thermal ground distribution. A single OBE steady state per velocity class gives
the pump-prepared populations; the weak probe reads them:

    α(Δ) = Σ_components Σ_v f(v) Σ_(Fg→Fe) S·(ρ_Fg − ρ_Fe)(Δ+k·v)
                                          · α_probe((Δ − k·v) − ω_(Fg,Fe))

Because the pump Hamiltonian is scan-independent (the scan enters only through
Δ_eff = Δ + k·v on the excited diagonal), the pump steady state is a 1-D function
of Δ_eff — solved once on a fine grid and interpolated (exact, fast; the same
Δ_eff trick the OD scheme uses). The CG decay branching transfers population into
the *other* ground hyperfine state (hyperfine optical pumping), which turns the
relevant crossovers into enhanced/inverted transmission peaks — the dominant
feature of real alkali SAS that a single-ground model cannot produce.

No pump → ρ_Fg = thermal, ρ_Fe = 0, and the expression reduces to the
Doppler-broadened multi-line absorption of the OD scheme.

Generic (Γ-units) mode
----------------------
A pedagogical fallback (one ground + one/two excited states, splitting set by
hand in natural-linewidth units) is kept under Advanced for teaching the bare
hole-burning picture without atomic data.

References:  Smith & Hughes, Am. J. Phys. 72, 631 (2004) (hyperfine pumping);
Preston, Am. J. Phys. 64, 1432 (1996); Steck alkali D-line data.
"""
import numpy as np

from .. import atoms, constants, doppler, observables, species
from ..constants import GAMMA, K_VEC
from .. import core
from .base import ParamSpec, Scheme

PROBE_RABI = 1e-3                       # weak probe, in units of Γ
GAMMA_MHZ = GAMMA / (2 * np.pi) / 1e6
GENERIC = "Generic (Γ units)"


class SASScheme(Scheme):
    name = "sas"
    cluster = "B — Sub-Doppler"
    title = "Saturated absorption (SAS)"
    caption = ("Counter-propagating pump + weak probe. Doppler-free Lamb dips and "
               "crossovers on the Doppler-broadened line — for a real alkali "
               "isotope/line (⁸⁵Rb, ⁸⁷Rb, ¹³³Cs · D1/D2) with hyperfine optical "
               "pumping, or a generic Γ-unit hole-burning model.")

    def param_schema(self):
        return [
            ParamSpec("species", "Atom / isotope", "Atomic", "Rb (natural)",
                      choices=tuple(species.SPECIES_ORDER) + (GENERIC,),
                      help="Real alkali species overlay the full hyperfine "
                      "spectrum (Natural Rb overlays ⁸⁵Rb+⁸⁷Rb by abundance). "
                      "Generic = bare Γ-unit hole-burning model."),
            ParamSpec("line", "Transition line", "Atomic", "D1",
                      choices=("D1", "D2"),
                      help="D1 (nP₁/₂) or D2 (nP₃/₂). Sets the excited hyperfine "
                      "manifold; ignored in Generic mode."),
            ParamSpec("temp_c", "Temperature", "Cell & beams", 25.0, 20.0, 200.0, 1.0, "°C",
                      help="Sets the vapor density (absorption scale) and Doppler width."),
            ParamSpec("cell_mm", "Cell length", "Cell & beams", 50.0, 0.5, 200.0, 0.5, "mm",
                      recompute=False),
            ParamSpec("pump_rabi", "Pump Rabi Ω_pump", "Fields", 2.0,
                      0.1, 10.0, 0.1, "Γ", advanced=True,
                      help="Saturating counter-propagating beam. Fixed at a sensible "
                      "value by default; tune here for power-broadening studies."),
            ParamSpec("transit_khz", "Transit relaxation γ_t/2π", "Atomic", 100.0,
                      5.0, 2000.0, 5.0, "kHz", advanced=True,
                      help="Atoms leaving/entering the beam relax toward the thermal "
                      "ground state — regularises hyperfine optical pumping. "
                      "Smaller γ_t → stronger inverted crossovers."),
            ParamSpec("line_strength", "Line-strength factor", "Detection & scaling", 1.0,
                      0.01, 2.0, 0.01, "", advanced=True, recompute=False,
                      help="Effective |d|² calibration knob (absolute OD scale)."),
            ParamSpec("transitions", "Generic: transitions", "Atomic", "single line",
                      choices=("single line", "two lines (crossover)"), advanced=True,
                      help="Generic mode only."),
            ParamSpec("splitting", "Generic: excited splitting", "Atomic", 60.0,
                      5.0, 200.0, 1.0, "Γ", advanced=True,
                      help="Generic mode only (two-line splitting)."),
            ParamSpec("scan_points", "Scan points", "Numerics", 1401,
                      401, 4001, 100, "", advanced=True),
        ]

    # No presets: the Atomic-section species/line dropdowns (and the Advanced
    # knobs for the generic mode) already cover every configuration directly.

    def info(self):
        return (
            "**Saturated absorption spectroscopy.** Counter-propagating pump and "
            "weak probe burn velocity-selective holes, revealing Doppler-free Lamb "
            "dips and crossovers. The realistic mode solves the full hyperfine "
            "manifold (CG-branched decay + transit relaxation), so **hyperfine "
            "optical pumping** turns the appropriate crossovers into enhanced / "
            "inverted transmission peaks.\n\n"
            "Atomic data (hyperfine A/B constants, line-centre frequencies, masses, "
            "linewidths) from the Steck alkali D-line data sheets; Wigner-6j line "
            "strengths validated against the ⁸⁵Rb D1 Clebsch-Gordan factors.\n\n"
            "**References**\n"
            "- D. A. Smith & I. G. Hughes, *Am. J. Phys.* **72**, 631 (2004) — "
            "hyperfine pumping in multilevel saturated absorption.\n"
            "- D. W. Preston, *Am. J. Phys.* **64**, 1432 (1996) — Rb SAS / crossovers.\n"
            "- D. A. Steck, *Rubidium 85 / 87 & Cesium D Line Data*, http://steck.us/alkalidata."
        )

    # =================================================================
    # compute  (dispatch)
    # =================================================================
    def compute(self, params):
        if params.get("species", "Rb (natural)") == GENERIC:
            return self._compute_generic(params)
        return self._compute_species(params)

    # ---- realistic multilevel species/line model ----
    def _compute_species(self, params):
        line = params["line"]
        comps = species.SPECIES[params["species"]]
        T = params["temp_c"] + 273.15
        Op_scale = params["pump_rabi"]
        gt = 2 * np.pi * params["transit_khz"] * 1e3        # transit rate [rad/s]

        # Common relative-frequency axis: reference to the highest-abundance
        # component's line centre (0 GHz at that centroid).
        iso_ref = max(comps, key=lambda c: c[1])[0]
        nu_ref = iso_ref.line(line)[1]

        built = []
        omega_all = []
        dopp_fwhm = 0.0
        for iso, weight in comps:
            man = species.build_manifold(iso, line, transit_rate=gt)
            offset = 2 * np.pi * (man.nu0 - nu_ref)         # axis shift [rad/s]
            sigma_v = np.sqrt(constants.KB * T / iso.mass)
            dfw = np.sqrt(8 * np.log(2)) * man.k_vec * sigma_v
            dopp_fwhm = max(dopp_fwhm, dfw)
            N = species.number_density(iso, T) * weight
            built.append(dict(man=man, offset=offset, sigma_v=sigma_v,
                              dfw=dfw, N=N))
            omega_all.append(man.omega + offset)
        omega_all = np.concatenate(omega_all)

        margin = 3.5 * dopp_fwhm
        lo, hi = omega_all.min() - margin, omega_all.max() + margin
        scan = np.linspace(lo, hi, int(params["scan_points"]))

        # Probe Lorentzian (bare weak-probe absorption of one unit transition,
        # density applied per component) on a fine grid near resonance.
        alpha = np.zeros(scan.size)
        markers = []
        for b in built:
            man = b["man"]
            alpha += self._component_alpha(man, scan, b["offset"], T,
                                           b["sigma_v"], b["N"], Op_scale)
            for k, om in enumerate(man.omega):
                markers.append((float((om + b["offset"]) / (2 * np.pi) / 1e9),
                                f"{man.iso.label} {man.Fg[man.g_idx[k]]:g}→{man.Fe[man.e_idx[k]-len(man.Fg)]:g}′"))

        return dict(mode="species", scan=scan, alpha_unit=alpha,
                    dopp_fwhm=dopp_fwhm, markers=markers,
                    species=params["species"], line=line, nu_ref=nu_ref)

    def _component_alpha(self, man, scan, offset, T, sigma_v, N, Op_scale):
        """Doppler-summed weak-probe absorption (1/m, ls=1) for one isotope."""
        n = man.n_levels
        gamma = man.gamma
        k_vec = man.k_vec

        v, wt = doppler.velocity_grid(T, mass=man.iso.mass, dv=3.0, cutoff_sigma=3.5)
        kv = k_vec * v

        # Pump steady state ρ(Δ_eff): scan-independent H, one fine table.
        Hp = species.pump_hamiltonian(man, Op_scale * gamma)
        L0 = core.build_liouvillian(Hp, man.atom)
        om = man.omega
        de_lo = om.min() - 12 * gamma
        de_hi = om.max() + 12 * gamma
        deff = np.linspace(de_lo, de_hi, int((de_hi - de_lo) / (gamma / 8.0)) + 2)
        pops = _pump_pops(L0, deff, man.atom.S_v, n)        # (n_deff, n_levels)

        # Bare weak-probe Lorentzian (unit transition, this density).
        two = atoms.two_level(gamma=gamma)
        Hpr = np.zeros((2, 2), dtype=complex)
        Hpr[0, 1] = Hpr[1, 0] = PROBE_RABI * gamma / 2
        L0pr = core.build_liouvillian(Hpr, two)
        pgrid = np.linspace(-18 * gamma, 18 * gamma,
                            int(36 * gamma / (gamma / 20.0)) + 2)
        rho_pr = core.steady_state_batched(L0pr, pgrid, two.S_v, 2)
        chi_pr = rho_pr[:, 1, 0] / (PROBE_RABI * gamma)
        alpha_L, _ = observables.absorption_coefficient(chi_pr, k_vec, N)

        # α(Δ) = Σ_v f(v) Σ_t S·(ρ_Fg − ρ_Fe)(Δ+kv) · α_L((Δ − offset − kv) − ω_t).
        DS = (scan - offset)                                 # laser-from-centroid (ns,)
        deff_grid = DS[:, None] + kv[None, :]                # (ns, nv)
        probe_base = DS[:, None] - kv[None, :]               # (ns, nv)
        pop_at = {lvl: np.interp(deff_grid.ravel(), deff, pops[:, lvl]).reshape(deff_grid.shape)
                  for lvl in set(man.g_idx.tolist()) | set(man.e_idx.tolist())}

        alpha = np.zeros(scan.size)
        for t in range(om.size):
            w = pop_at[man.g_idx[t]] - pop_at[man.e_idx[t]]      # (ns, nv)
            Lp = np.interp((probe_base - om[t]).ravel(), pgrid, alpha_L,
                           left=0.0, right=0.0).reshape(w.shape)
            alpha += man.S_abs[t] * ((w * Lp) @ wt)
        return alpha

    # ---- generic Γ-unit hole-burning toy (pedagogical) ----
    def _compute_generic(self, params):
        two = params["transitions"].startswith("two")
        if two:
            split = params["splitting"] * GAMMA
            n_exc, offsets = 2, np.array([-split / 2, split / 2])
        else:
            n_exc, offsets = 1, np.array([0.0])

        atom = atoms.sas_atom(n_exc)
        Op = params["pump_rabi"] * GAMMA
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

        two_lvl = atoms.two_level()
        Hpr = np.zeros((2, 2), dtype=complex)
        Hpr[0, 1] = Hpr[1, 0] = PROBE_RABI * GAMMA / 2
        L0_probe = core.build_liouvillian(Hpr, two_lvl)
        kvmax = float(np.abs(kv).max())
        flo, fhi = scan.min() - kvmax - off_span, scan.max() + kvmax + off_span
        fine = np.linspace(flo, fhi, int((fhi - flo) / (GAMMA / 20)) + 2)
        rho_pr = core.steady_state_batched(L0_probe, fine, two_lvl.S_v, 2)
        chi_pr = rho_pr[:, 1, 0] / (PROBE_RABI * GAMMA)
        alpha_L, _ = observables.absorption_coefficient(chi_pr, K_VEC, N)

        alpha = np.zeros(scan.size)
        for j, D in enumerate(scan):
            rho = core.steady_state_batched(L0_pump, D + kv, atom.S_v, atom.n_levels)
            rho_gg = rho[:, 0, 0].real
            contrib = np.zeros(v.size)
            for i, e in enumerate(atom.excited):
                w_i = rho_gg - rho[:, e, e].real
                contrib += w_i * np.interp((D - kv) - offsets[i], fine, alpha_L)
            alpha[j] = float((wt * contrib).sum())

        return dict(mode="generic", scan=scan, alpha_unit=alpha,
                    dopp_fwhm=dopp_fwhm, offsets=offsets, two=(n_exc == 2))

    # =================================================================
    # observables  (dispatch)
    # =================================================================
    def observables(self, raw, params):
        L = params["cell_mm"] * 1e-3
        ls = params["line_strength"]
        alpha = raw["alpha_unit"] * ls
        if raw["mode"] == "species":
            return self._obs_species(raw, params, alpha, L)
        return self._obs_generic(raw, params, alpha, L)

    def _obs_species(self, raw, params, alpha, L):
        import matplotlib.pyplot as plt
        x = raw["scan"] / (2 * np.pi) / 1e9                  # GHz (relative)
        T_trans = observables.transmission(alpha, L)
        OD = observables.optical_density(alpha, L)
        dopp_mhz = raw["dopp_fwhm"] / (2 * np.pi) / 1e6

        fig, (axT, axA) = plt.subplots(2, 1, figsize=(8.5, 6.4), sharex=True)
        axT.plot(x, T_trans, color="#1f77b4", lw=1.3)
        axT.set_ylabel("Transmission")
        axT.set_title(f"SAS — {raw['species']} {raw['line']}:  "
                      f"T = {params['temp_c']:.0f} °C, L = {params['cell_mm']:.0f} mm, "
                      f"Ω_pump = {params['pump_rabi']:.1f} Γ")
        axA.plot(x, OD, color="#d62728", lw=1.3)
        axA.set_ylabel("Optical density")
        axA.set_xlabel("Relative frequency  [GHz]  (ref: line centroid)")
        for gx, _lbl in raw["markers"]:
            for ax in (axT, axA):
                ax.axvline(gx, color="gray", ls=":", lw=0.5, alpha=0.6)
        fig.tight_layout()

        sub_fwhm, sub_at = _narrowest_subdoppler(x, T_trans)
        metrics = [
            dict(label="Doppler FWHM", value=f"{dopp_mhz:.0f} MHz",
                 help="Width of the Doppler-broadened background lines."),
            dict(label="Narrowest sub-Doppler", value=f"{sub_fwhm*1e3:.1f} MHz"
                 if np.isfinite(sub_fwhm) else "—",
                 help="FWHM of the sharpest Doppler-free feature "
                 f"(near {sub_at:.2f} GHz)." if np.isfinite(sub_fwhm) else ""),
            dict(label="Peak OD", value=f"{np.nanmax(OD):.2f}"),
        ]
        rows = "".join(f"| {lbl} | {gx*1e3:.1f} |\n" for gx, lbl in raw["markers"])
        table = ("Hyperfine transitions (Lamb-dip centres); crossovers appear at "
                 "the midpoint of any two sharing a ground state, enhanced/inverted "
                 "by hyperfine pumping.\n\n"
                 "| Transition | Center [MHz] |\n|---|---|\n" + rows)
        return dict(metrics=metrics, figure=fig,
                    tables=[{"title": "Hyperfine lines", "markdown": table}])

    def _obs_generic(self, raw, params, alpha, L):
        import matplotlib.pyplot as plt
        x = raw["scan"] / (2 * np.pi) / 1e6                  # MHz
        T_trans = observables.transmission(alpha, L)
        OD = observables.optical_density(alpha, L)
        offs_mhz = raw["offsets"] / (2 * np.pi) / 1e6

        fig, (axT, axA) = plt.subplots(2, 1, figsize=(8.5, 6.4), sharex=True)
        axT.plot(x, T_trans, color="#1f77b4", lw=1.6)
        axT.set_ylabel("Transmission")
        axT.set_title(f"SAS (generic): Ω_pump = {params['pump_rabi']:.1f} Γ, "
                      f"T = {params['temp_c']:.0f} °C")
        axA.plot(x, OD, color="#d62728", lw=1.6)
        axA.set_ylabel("Optical density")
        axA.set_xlabel("Probe detuning  [MHz]")
        for ax in (axT, axA):
            for off in offs_mhz:
                ax.axvline(off, color="gray", ls=":", lw=0.7)
            if raw["two"]:
                ax.axvline(0.0, color="green", ls=":", lw=0.7)
        fig.tight_layout()

        ic = int(np.argmin(np.abs(x - offs_mhz[0])))
        sub_fwhm = _peak_fwhm(x, T_trans, ic)
        dopp_mhz = raw["dopp_fwhm"] / (2 * np.pi) / 1e6
        metrics = [
            dict(label="Doppler FWHM", value=f"{dopp_mhz:.0f} MHz"),
            dict(label="Sub-Doppler dip FWHM", value=f"{sub_fwhm:.1f} MHz"),
            dict(label="Peak OD", value=f"{np.nanmax(OD):.2f}"),
        ]
        note = ("Two transitions: Lamb dips at ±splitting/2 and a **crossover** dip "
                "at the midpoint (green)." if raw["two"]
                else "Single transition: one Lamb dip at line centre.")
        return dict(metrics=metrics, figure=fig,
                    tables=[{"title": "Notes", "markdown": note}])


# =====================================================================
# helpers
# =====================================================================
def _pump_pops(L0, deff_axis, S_v, n, chunk=1500):
    """Diagonal populations ρ_ii(Δ_eff) on a fine axis, solved in memory-safe chunks."""
    pops = np.empty((deff_axis.size, n))
    for s in range(0, deff_axis.size, chunk):
        sl = slice(s, s + chunk)
        rho = core.steady_state_batched(L0, deff_axis[sl], S_v, n)
        pops[sl] = np.einsum("vii->vi", rho).real
    return pops


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


def _narrowest_subdoppler(x_ghz, T_trans):
    """
    Width [GHz] and location of the sharpest Doppler-free feature: detrend the
    transmission with a broad median, then find the narrowest prominent peak in
    |residual|. Returns (fwhm_ghz, center_ghz); nan if none stands out.
    """
    y = np.asarray(T_trans)
    n = y.size
    win = max(5, (n // 60) | 1)                              # broad smoothing window
    pad = win // 2
    ypad = np.pad(y, pad, mode="edge")
    smooth = np.array([np.median(ypad[i:i + win]) for i in range(n)])
    resid = np.abs(y - smooth)
    if resid.max() <= 1e-9:
        return float("nan"), float("nan")
    ic = int(np.argmax(resid))
    half = 0.5 * resid[ic]
    i = ic
    while i > 0 and resid[i] >= half:
        i -= 1
    j = ic
    while j < n - 1 and resid[j] >= half:
        j += 1
    return float(x_ghz[j] - x_ghz[i]), float(x_ghz[ic])
