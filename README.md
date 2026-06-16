# GABES — Generic Atomic Bloch Equation Solver

Ensemble optical-Bloch-equation solver for warm-vapor / cold-atom spectroscopy,
with a **scheme-driven** Streamlit front-end: pick a physics scheme in the
sidebar and it declares its own controls and plots. Started as a single 85Rb D1
double-Λ four-wave-mixing model (now one scheme among several).

## Schemes (current)

| Cluster | Scheme | Output |
|---|---|---|
| A — Absorption | **OD / SAS** | weak-probe absorption with a counter-propagating pump. Pump off → Doppler-broadened OD (validated ⁸⁵Rb D1 hyperfine scale); pump on → Doppler-free Lamb dips + crossovers with **hyperfine optical pumping**. ⁸⁵Rb / ⁸⁷Rb / ¹³³Cs · D1/D2 or natural Rb; generic Γ-unit fallback |
| A | **Lambda coherence (EIT / AT / CPT)** | one 3-level Lambda engine with regime-driven defaults, physical MHz/kHz controls, Rb/Cs D-line media, EIT transparency, AT splitting, and CPT dark resonance |
| A | **Rydberg-EIT electrometry** | 85Rb cascade EIT / microwave AT static spectrum for the 5S-5P-40D ladder and 37 GHz 40D-39F RF leg; reference sensitivity numbers stay in internal tests |
| C — Magneto-optics | **Hanle / EIA / NMOR** | two distinct effects vs B: the **Hanle** effect (zero-field transmission dip/peak, EIA variant) from ground-state coherence, and **magneto-optical rotation** (MOR/NMOR, polarization-plane rotation) — both over the Zeeman manifold |
| D — Wave mixing | **FWM** | legacy seeded 85Rb D1 double-Λ gain/squeezing, plus generic SFWM biphoton source estimates (`g²_SI(τ)`, CAR, rates, phase matching, velocity-class BTW) |

Roadmap (parking lot): slow-light / group-index readout, Raman gain, higher-order
wave mixing, Bell-Bloom magnetometry, Na D-lines (SAS species data); time-domain
(STIRAP, Ramsey) and two-time correlations (Mollow, g²(τ)) would need new engine layers.

## Layout

- `gabes/` — the package:
  - `constants.py` — physical constants + 85Rb D1 line data + `rabi_freq`.
  - `core.py` — physics-agnostic engine: super-operators, Liouvillian, the 3-mode
    Floquet solver, the single-mode steady-state solver, 2×2 matrix exp.
  - `doppler.py` — Maxwell velocity grid, Δ_eff axis, Doppler average.
  - `atoms.py` — `AtomModel` (level scheme as data) + factories
    (`two_level`, `lambda3`, `sas_atom`, `double_lambda_rb85`) + Rb85 vapor density.
  - `hyperfine.py` — 85Rb D1 hyperfine line table (4 transitions + shifts), CG
    line strengths C_F², ground populations p_F, self-broadening Γ(N), and the
    pure-85Rb CRC density — the data/scaling for the AutoOD-validated full-D1 OD.
  - `zeeman.py` — hand-rolled Clebsch-Gordan + `zeeman_manifold(F_g, F_e)` builder
    (σ±/π couplings, CG-branched decay) for the magneto-optics schemes.
  - `species.py` — alkali D-line data (⁸⁵Rb/⁸⁷Rb/¹³³Cs hyperfine A/B constants,
    line-centre frequencies, masses, linewidths, abundances; Steck), Wigner-6j
    (Racah) line strengths, Casimir hyperfine energies, Steck vapor density, and
    the SAS hyperfine-manifold builder `build_manifold(iso, line)` (CG-branched
    decay + transit-time relaxation).
  - `observables.py` — gain, squeezing, legacy twin-beam coincidence, calibrated biphoton statistics, absorption / OD / dispersion.
  - `schemes/` — experiment plugins: `base.py` (`Scheme`/`ParamSpec`/`Preset`/`ExtraView`),
    `absorption.py` (Lambda EIT/AT/CPT + the unregistered `ODScheme` validation
    primitive), `rydberg.py` (Rydberg-EIT electrometry), `sas.py`
    (the merged **Absorption OD/SAS** scheme on `species.py`),
    `magneto.py` (Hanle/EIA/NMOR),
    `fwm.py`, `__init__.py` (registry).
- `streamlit_app.py` — generic UI. Renders only the selected scheme's
  `param_schema()`; caches the heavy solve on `recompute` knobs only (so
  navigate-only knobs like the FWM two-photon detuning update instantly).
- `fwm_obe.py` — backward-compat shim re-exporting the FWM API and the
  `python fwm_obe.py` CLI (physics now lives in `gabes/`).
- `tests/` — regression + physics validation; `baseline_focused.npz` is the
  frozen pre-refactor FWM anchor.
- `requirements.txt` — streamlit, numpy, matplotlib.
- `references/` — Sim et al. 2025 paper PDF + Steck Rb85 data + OE stabilization paper. The Rydberg-EIT scheme cites arXiv:2606.04354 for its 85Rb reference defaults; the Biphoton defaults cite the Cs biphoton-temporal-waveform paper and the 87Rb telecom biphoton source paper from the app's Reference panel.

## Run

```
streamlit run streamlit_app.py
```

**Trap**: `python streamlit_app.py` does NOT work (no server, prints "run it with
streamlit run"). FWM CLI backend: `python fwm_obe.py`.

## Tests

```
python tests/test_regression.py      # FWM bit-identical to the pre-refactor baseline
python tests/test_absorption.py      # OD width, full-D1 AutoOD scale + line ratio, Lambda AT/EIT/CPT
python tests/test_rydberg_eit.py     # 85Rb Rydberg-EIT reference defaults, linewidth, RF AT split
python tests/test_sas.py             # 6j↔CF2, HF splittings, no-pump→OD (49/25), hyperfine-pumping crossovers, generic

python tests/test_magneto.py         # CG values, Hanle dip, EIA peak, NMOR zero-crossing
python tests/test_coincidence.py     # twin-beam photon-pair statistics
python tests/test_fwm_generic.py     # generic SFWM topology + biphoton detector model
python tests/test_schemes_render.py  # every registered scheme computes + renders
```
(or `pytest tests/`)

## Adding a scheme

Subclass `gabes.schemes.base.Scheme` — declare `param_schema()`, `compute(params)`,
`observables(raw, params)` (and optionally `presets`, `info`, `extra_views`) — then
add an instance to the list in `gabes/schemes/__init__.py`. No UI edits: the
sidebar controls and the plots follow `param_schema()` and the observables dict.

## FWM conventions (must-know, not obvious from numbers)

### Legacy seeded gain / squeezing

- Levels: g₁=F=2, g₂=F=3, e₂=F'=2, e₃=F'=3.
- OPD Δ (one-photon): ω_pump = ω(F=2→F'=3) + Δ.
- TPD δ (two-photon): ω_seed = ω_pump − ν_HF + δ.   ν_HF = 3.0357 GHz.
- Plot x-axis ref = **F=2→F'=3** line. (−) Raman branch = standard FWM seed, at Δ − ν_HF.
- Beam waists W_PUMP=530 µm, W_PROBE=330 µm = **1/e² radius** (paper convention). Not diameter.

### Generic SFWM / biphoton mode

- Mode selector: **Squeezing** keeps the regression-anchored 85Rb
  double-Λ model; **Biphoton** switches to the generic SFWM source
  estimate.
- Topologies: `cascade_rb87_telecom` (87Rb 5S1/2-5P3/2-4D5/2, 780/1529 nm),
  `cascade_cs_btw` (133Cs 852-917 nm or 852-795 nm BTW comparison), and
  `diamond_generic` (four-level user-wavelength template; not a validated paper
  preset).
- Biphoton readout: `g²_SI(τ)`, FWHM, pair-rate estimate, singles, true and
  accidental coincidences, CAR, heralding estimates, and Cauchy-Schwarz R.
- The waveform is a coherent sum over Doppler velocity classes. Biphoton v3 phase
  matching uses calibrated longitudinal Δk plus absolute transverse Δk, with a
  strict vector `sinc²(|Δk| L / 2)` collection weight.
- Reference anchors: the 87Rb telecom source is calibrated to order
  `g²_SI≈44`, OD≈112, bandwidth≈300 MHz, and coincidence rate≈38,000 cps/mW;
  the Cs BTW preset exposes the wavelength-dependent temporal-width change
  reported for the 852-917 nm and 852-795 nm cascade channels.
- v3 limitation: absolute pair rate is a **calibrated source estimate**, not a
  full quantum-Langevin propagation/noise calculation.

## Traps

Current FWM note: the 4-level FWM model applies the real 85Rb D1 hyperfine
`C_F^2` values `(10,35,35,28)/81` to Rabi couplings, polarization readout, and
spontaneous-emission branching. The macroscopic-coupling scale is now
**first-principles**: `fwm.physical_coupling_norm` supplies `p_F/[2(2I+1)]` (the
ground-population fraction × ground-sublevel degeneracy, ≈0.0486 on the (−)
branch) — the factor the validated absorption path applies explicitly and the
lumped 4-level + total-density convention otherwise omits. Density also uses the
pure-85Rb CRC fit (`hyperfine.number_density`), consistent with that path. The
**Line-strength factor** is now only a dimensionless **residual** calibration
(default `1.0`), parked under *Advanced* because it is not an experimentally
tunable variable. Ultra fidelity adds the slow propagation refinements, but the
full 24-level Zeeman Floquet scan is still reported as a diagnostic rather than
used as the default full-scan solver.

1. **FWM gain is exponentially sensitive at high density.** At paper optimum
   T=121 °C linear Maxwell-Bloch still over-amplifies; the residual Line-strength
   factor multiplies the physical coupling on top of `p_F/[2(2I+1)]`, so leave it
   at `1.0` unless an absolute-gain measurement says otherwise. Pump-depletion
   (Manley-Rowe) saturation caps the runaway; squeezing depth is η-limited.
2. **FWM Raman branches are separate mode pairs, not one summed susceptibility.**
   The Sim et al. 85Rb operating point uses the standard red-detuned seed on the
   (−) Raman branch. Compute `branch=-1` and `branch=+1` independently; do not add
   their χ matrices into one 2×2 propagation matrix. The old branch-summed model
   created artificial high-gain extrema (for example near +70 MHz TPD).
3. **Rb is very absorbing — absorption schemes use short cells.** ls=1.0 is the
   true cross-section, so on-resonance OD saturates in a cm-scale cell; the
   OD/EIT/AT/CPT defaults use mm-scale cells / moderate T to keep features visible.
4. **Twin-beam coincidence is the ideal (lossless) parametric estimate** from the
   gain (n=G_s−1, g²_sc=2+1/n, R=g²_sc²/4>1), valid in the gain region — like the
   squeezing panel, propagation loss is not modelled with quantum Langevin noise.
5. **OD is the pump-off limit of the Absorption (OD/SAS) scheme.** Pump power = 0
   → linear Doppler-broadened absorption; raising it burns the SAS sub-Doppler
   features on the *same* spectrum. For ⁸⁵Rb D1 the pump-off limit reproduces the
   lab AutoOD calculator (`references/AutoOD/`) to <0.1 %: the absolute scale uses
   the CRC vapor-pressure density (Rb) and the AutoOD C_F²·|d|² normalisation
   (`species.line_integrated_alpha`, `species.cf2`/`reduced_dipole_sq`). The old
   **single 2-level** OD model is kept as an internal validation primitive
   (`schemes.absorption.ODScheme`, *no longer registered*) that the Λ schemes
   reduce to and the analytic FWHM=Γ tests use. The probe is fixed weak; only the
   **pump power [mW]** is a knob (→ Rabi via I=2P/πw² and I_sat).
6. **SAS line weight is `(2Fg+1)·line_strength`, not `line_strength`.** The
   observable strength of a lumped Fg↔Fe hyperfine line — absorption per ground
   atom *and* the spontaneous-emission branching — carries the ground degeneracy:
   `T = (2Fg+1)(2Fe+1)(2Jg+1){6j}²`. One quantity drives the line weight, the
   CG decay branching (→ hyperfine pumping) and (√T) the relative pump Rabi. It
   reproduces the validated ⁸⁵Rb D1 `CF2` (`T = 9·CF2`) and the 49/25 F=3/F=2
   manifold ratio. **Hyperfine pumping** (decay into the *other* ground state)
   is what turns crossovers into enhanced/inverted transmission peaks — the
   dominant feature of real alkali SAS; a single-ground model cannot make them.
   A transit-time relaxation `γ_t` (atoms leaving the beam, an Advanced knob)
   regularises the pumping: without it the dark ground state saturates, and a
   smaller `γ_t` gives stronger inverted crossovers.

## Speed (why the architecture)

- L₀(Δ_eff) = L₀_base − Δ_eff·S_v (only the excited diagonal shifts with velocity).
  → all velocities stacked, one batched `np.linalg.solve`.
- The χ̄(δ, Δ_eff) table is **T- and Δ-independent** (T → Maxwell weights only; Δ → axis offset).
- Full FWM −8…12 GHz scan ~227 s. App focused window ~6 s/recompute, `st.cache_data`.
  TPD slider = instant (navigates the cached curve, no recompute).

## FWM future physics work

- Ultra fidelity includes fixed-iteration refractive phase matching, 64-segment
  propagation, Gaussian overlap, segment-wise pump-budget depletion, and in-cell
  loss/noise before detector efficiency.
- Technical noise channels (seed excess noise, pump scattering, EOM residual
  carrier/sideband noise) are present as zero-default internal terms until
  calibrated measurements are available.
- The 24-level 85Rb D1 Zeeman manifold is built for CG diagnostics and correction
  bookkeeping; the full Zeeman Floquet solve remains future work because the
  density matrix jumps from 4-level `M=16` to 24-level `M=576`.

## Sim et al. 85Rb optimum (Squeezing default, FWM scheme)

G. Sim, H. Kim, H. S. Moon, Sci. Rep. **15**, 7727 (2025). 85Rb squeezing-optimal:

| Δ | δ | T | pump | seed | loss | → result |
|---|---|---|---|---|---|---|
| 0.9 GHz | −8 MHz | 121 °C | 600 mW | 8 µW | 5.5 % | gain ≈ 15, IDS −7.8 dB |

Fixed geom: cell L=12.5 mm, QE 90.47 %, responsivity 0.58 A/W @ 795 nm, pump⊥probe.

## Deploy

- **Streamlit Community Cloud** (easiest), or Render/Railway/Fly/Cloud Run via Docker.
- Free tier = weak CPU → 6 s recompute slower. Real bottleneck = CPU, not host.
- Git remote: `github.com/Shake2313/fwm-squeezing-app` (private). gh acct `Shake2313`.
  Repo name ≠ folder name.
