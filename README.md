# GABES — Generic Atomic Bloch Equation Solver

Ensemble optical-Bloch-equation solver for warm-vapor / cold-atom spectroscopy,
with a **scheme-driven** Streamlit front-end: pick a physics scheme in the
sidebar and it declares its own controls and plots. Started as a single 85Rb D1
double-Λ four-wave-mixing model (now one scheme among several).

## Schemes (current)

| Cluster | Scheme | Output |
|---|---|---|
| A — Absorption | **OD / SAS** | weak-probe absorption with a counter-propagating pump. Pump off → Doppler-broadened OD (validated ⁸⁵Rb D1 hyperfine scale); pump on → Doppler-free Lamb dips + crossovers with **hyperfine optical pumping**. ⁸⁵Rb / ⁸⁷Rb / ¹³³Cs · D1/D2 or natural Rb; generic Γ-unit fallback |
| A | **EIT** | transparency window + dispersion (slow light) |
| A | **AT** | Autler-Townes doublet (splitting = Ω_c) |
| A | **CPT** | sub-natural dark resonance |
| C — Magneto-optics | **Hanle / EIA / NMOR** | zero-field dip / peak / polarization rotation vs B (Zeeman manifold) |
| D — Wave mixing | **FWM** | seed/probe gain G_s, intensity-difference squeezing, twin-beam coincidence |

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
  - `observables.py` — gain, squeezing, twin-beam coincidence, absorption / OD / dispersion.
  - `schemes/` — experiment plugins: `base.py` (`Scheme`/`ParamSpec`/`Preset`/`ExtraView`),
    `absorption.py` (EIT/AT/CPT + the unregistered `ODScheme` validation
    primitive), `sas.py` (the merged **Absorption OD/SAS** scheme on `species.py`),
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
- `references/` — Sim et al. 2025 paper PDF + Steck Rb85 data + OE stabilization paper.

## Run

```
streamlit run streamlit_app.py
```

**Trap**: `python streamlit_app.py` does NOT work (no server, prints "run it with
streamlit run"). FWM CLI backend: `python fwm_obe.py`.

## Tests

```
python tests/test_regression.py      # FWM bit-identical to the pre-refactor baseline
python tests/test_absorption.py      # OD width, full-D1 AutoOD scale + line ratio, AT split = Ω_c, EIT/CPT
python tests/test_sas.py             # 6j↔CF2, HF splittings, no-pump→OD (49/25), hyperfine-pumping crossovers, generic

python tests/test_magneto.py         # CG values, Hanle dip, EIA peak, NMOR zero-crossing
python tests/test_coincidence.py     # twin-beam photon-pair statistics
python tests/test_schemes_render.py  # every registered scheme computes + renders
```
(or `pytest tests/`)

## Adding a scheme

Subclass `gabes.schemes.base.Scheme` — declare `param_schema()`, `compute(params)`,
`observables(raw, params)` (and optionally `presets`, `info`, `extra_views`) — then
add an instance to the list in `gabes/schemes/__init__.py`. No UI edits: the
sidebar controls and the plots follow `param_schema()` and the observables dict.

## FWM conventions (must-know, not obvious from numbers)

- Levels: g₁=F=2, g₂=F=3, e₂=F'=2, e₃=F'=3.
- OPD Δ (one-photon): ω_pump = ω(F=2→F'=3) + Δ.
- TPD δ (two-photon): ω_seed = ω_pump − ν_HF + δ.   ν_HF = 3.0357 GHz.
- Plot x-axis ref = **F=2→F'=3** line. (−) Raman branch = standard FWM seed, at Δ − ν_HF.
- Beam waists W_PUMP=530 µm, W_PROBE=330 µm = **1/e² radius** (paper convention). Not diameter.

## Traps

Current FWM note: the 4-level FWM model now applies the real 85Rb D1 hyperfine
`C_F^2` values `(10,35,35,28)/81` to Rabi couplings, polarization readout, and
spontaneous-emission branching. The remaining app knob is a residual **FWM
coupling scale**, not a literal line-strength constant; keep it near the
calibrated regime (~0.05) until phase mismatch, loss/noise, pump depletion, and
full Zeeman structure are included.

1. **FWM gain is exponentially sensitive at high density.** At paper optimum
   T=121 °C linear Maxwell-Bloch still over-amplifies if the residual coupling
   scale is pushed to 1.0. Keep it near the calibrated regime (~0.05) unless you
   are deliberately stress-testing the ideal model.
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

- Add longitudinal phase mismatch `Δk_z` to the probe/conjugate propagation
  matrix. Sim et al. use a finite pump-probe angle, and Turnbull et al. show that
  phase mismatch controls the gain/absorption tradeoff for squeezed twin beams.
- Add explicit loss/noise physics rather than folding everything into a scalar
  detection efficiency: propagation loss with Langevin noise, seed excess noise,
  pump scattering, and EOM residual carrier/sideband multimode noise.
- Benchmark the cost/benefit of those additions before making them default. Record
  runtime slowdown, shift in optimum TPD, gain match to Sim et al. (~15), and IDS
  accuracy improvement versus the current single-branch model.

## Sim et al. 85Rb optimum (the ⚡ button, FWM scheme)

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
