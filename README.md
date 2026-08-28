# GABES — Generic Atomic Bloch Equation Solver

Ensemble optical-Bloch-equation solver for warm-vapor / cold-atom spectroscopy,
with a **scheme-driven** Streamlit front-end: pick a physics scheme in the
sidebar and it declares its own controls and plots. Started as a single 85Rb D1
double-Λ four-wave-mixing model (now one scheme among several).

## Schemes (current)

| Cluster | Scheme | Output |
|---|---|---|
| A — Absorption | **OD / SAS** | weak-probe absorption with a counter-propagating pump. Pump off → Doppler-broadened OD (validated ⁸⁵Rb D1 hyperfine scale); pump on → Doppler-free Lamb dips + crossovers with **hyperfine optical pumping**. An Advanced paraffin-cell switch adds ground-hyperfine population memory between velocity-randomized beam passages without adding a coating-throughput loss. Imports oscilloscope A/B CSV data for robust correction and manual transmission-overlay alignment. ⁸⁵Rb / ⁸⁷Rb / ¹³³Cs · D1/D2 or natural Rb; generic Γ-unit fallback |
| A | **Lambda coherence (EIT / AT / CPT)** | one 3-level Lambda engine with regime-driven defaults, physical MHz/kHz controls, Rb/Cs D-line media, EIT transparency, AT splitting, and CPT dark resonance |
| A | **Rydberg-EIT electrometry** | 85Rb cascade EIT / microwave AT spectrum for the 5S-5P-40D ladder and 37 GHz 40D-39F RF leg; first-order finite-IF weak-SIG response, AT→field calibration, balanced-detector PSN/RIN/electronics budget, temperature/cold-spot sweep, and conditional absolute sensitivity |
| C — Magneto-optics | **Hanle / EIA / NMOR** | two distinct effects vs B: the **Hanle** effect (zero-field transmission dip/peak, EIA variant) from ground-state coherence, and **magneto-optical rotation** (MOR/NMOR, polarization-plane rotation) — both over the Zeeman manifold |
| D — Wave mixing | **FWM** | seeded 85Rb D1 double-Λ mean-field gain diagnostic (physical squeezing claim-gated), plus generic SFWM biphoton source estimates (`g²_SI(τ)`, CAR, rates, phase matching, velocity-class BTW) |

Roadmap (parking lot): slow-light / group-index readout, Raman gain, higher-order
wave mixing, Bell-Bloom magnetometry, Na D-lines (SAS species data); time-domain
(STIRAP, Ramsey) and two-time correlations (Mollow, g²(τ)) would need new engine layers.

## Layout

- `gabes/` — the package:
  - `constants.py` — physical constants + 85Rb D1 line data + `rabi_freq`.
  - `core.py` — physics-agnostic engine: super-operators, Liouvillian, general
    finite-order Floquet continued fractions plus an independent dense-block
    reference, the single-mode steady-state solver, projected trace-zero
    Liouvillian response and pole-residue checks, and 2×2 matrix exp.
  - `doppler.py` — Maxwell velocity grids, Δ_eff axis, one-dimensional Doppler
    average, and the opt-in two-dimensional non-collinear Raman-Doppler reference
    quadrature.
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
  - `experimental_csv.py` — strict A/B oscilloscope CSV import, robust denoising,
    explicitly relative extrema normalization or evidence-backed absolute
    transmission calibration, acquisition-order sweep diagnostics, and manual
    x-axis transforms.
  - `rydberg_electrometry.py` — finite-IF Liouvillian response, explicit RF
    dipole/peak-RMS conventions, balanced-detector noise, and absolute field ASD.
  - `rydberg_experiment.py` — heater/effective/cold-spot temperature states,
    axial density/Beer–Lambert integration, effective atom number, SAM uncertainty,
    and temperature/density dephasing fits.
  - `rydberg_experimental_csv.py` — unit-aware EIT/RF/PSD trace import with
    duplicate handling and SHA-256 provenance.
  - `schemes/` — experiment plugins: `base.py` (`Scheme`/`ParamSpec`/`Preset`/`ExtraView`),
    `absorption.py` (Lambda EIT/AT/CPT + the unregistered `ODScheme` validation
    primitive), `rydberg.py` (Rydberg-EIT electrometry), `sas.py`
    (the merged **Absorption OD/SAS** scheme on `species.py`),
    `magneto.py` (Hanle/EIA/NMOR),
    `fwm.py`, `__init__.py` (registry).
- `streamlit_app.py` — generic UI. Renders only the selected scheme's
  `param_schema()`; caches the heavy solve on `recompute` knobs only (so
  navigate-only knobs like the FWM two-photon detuning update instantly).
- `archive/` — retired legacy entry points and artifacts. It preserves the former
  `fwm_obe.py` compatibility shim and its generated
  `fwm_double_lambda_gain_squeezing.png`; the active physics lives in `gabes/`.
- `tests/` — regression + physics validation; `baseline_focused.npz` is an
  intentional compatibility snapshot documented by `baseline_focused_manifest.json`.
- `analysis/rydberg_cell_heating/` — config-driven temperature sweep, clean CSV
  overlays, provenance manifest, figures, and XeLaTeX report inputs.
- `requirements.txt` — streamlit, numpy, matplotlib.
- `references/` — Sim et al. 2025 paper PDF + Steck Rb85 data + OE stabilization paper. The Rydberg-EIT scheme cites arXiv:2606.04354 for its 85Rb reference defaults; the Biphoton defaults cite the Cs biphoton-temporal-waveform paper and the 87Rb telecom biphoton source paper from the app's Reference panel.

## Run

```
streamlit run streamlit_app.py
```

**Trap**: `python streamlit_app.py` does NOT work (no server, prints "run it with
streamlit run"). The supported application entry point is the Streamlit command
above. Root-level `fwm_obe` CLI/import compatibility is retired; its shim is kept
under `archive/` for historical reference only.

### Experimental CSV comparison (OD / SAS)

Open **Experimental CSV comparison** above the OD/SAS plot and upload a CSV.
Column A is always detuning and column B is always the detector signal; later
columns and non-numeric oscilloscope metadata rows are ignored. Raw A/B values
are treated as arbitrary units. GABES merges repeated detuning samples, removes
isolated spikes, and smooths according to the measured noise. Extrema scaling
maps the robust signal floor/ceiling to **relative normalized transmission**
0/1; it is not an absolute transmission calibration. Absolute transmission
requires measured dark/reference levels or an explicit detector gain and
offset. Sweep reversals are counted before sorting, and the original forward
and reverse branches are retained; the compatibility overlay still merges
equal detunings and warns that this loses hysteresis information. Use x scale,
x shift, and the sweep-direction control to align the corrected trace with the
simulated frequency axis. The graph opens on **Transmission**; use desktop arrow
buttons or a mobile left/right swipe to switch to **Optical density**. These
display-only controls do not rerun the Bloch-equation solve.

### Rydberg finite-IF electrometry and cell heating

The Rydberg scheme separates heater setpoint, effective vapor temperature, and
cold-spot temperature; the cold spot fixes vapor pressure while the effective
temperature controls motion and temperature-dependent dephasing. In the AT view,
GABES solves the first-order weak-SIG response around the LO-dressed steady state
at the declared IF, then refers balanced-detector PSN, RIN, and electronics noise
back to RF field through an explicit 40D–39F dipole and peak/RMS convention.
Displayed absolute sensitivities are conditional predictions: the literature
11.2/12.5 nV/cm/√Hz values are not injected as anchors.

The advanced SAM controls convert source power/gain/loss/distance to field with a
standard-uncertainty budget and flag a failed `2D²/λ` far-field criterion. The
batch workflow additionally supports axial `T(z),n(z)`, segmented Beer–Lambert,
effective atom number, and unit-aware EIT/RF/PSD CSV overlays:

```powershell
python -m analysis.rydberg_cell_heating.workflow --config analysis/rydberg_cell_heating/example_config.json
```

Full time-domain LO+SIG/lock-in processing, full Zeeman/polarization, spatial
horn/standing-wave RF fields, and three-photon imaging remain deliberately
deferred.

## Tests

```
python tests/test_regression.py      # FWM matches the versioned compatibility baseline
python tests/test_absorption.py      # OD width, full-D1 AutoOD scale + line ratio, Lambda AT/EIT/CPT
python tests/test_rydberg_eit.py     # 85Rb Rydberg-EIT reference defaults, linewidth, RF AT split
python tests/test_rydberg_electrometry.py # finite-IF response and detector/noise conventions
python tests/test_rydberg_experiment.py # thermal/axial/SAM/dephasing experiment helpers
python tests/test_rydberg_experimental_csv.py # unit-aware EIT/RF/PSD import
python -m pytest -q analysis/rydberg_cell_heating/test_workflow.py # report workflow
python tests/test_sas.py             # 6j↔CF2, HF splittings, no-pump→OD (49/25), hyperfine-pumping crossovers, generic
python tests/test_experimental_csv.py # A/B parsing, correction, ReferenceOD, manual alignment

python tests/test_magneto.py         # CG values, Hanle dip, EIA peak, NMOR zero-crossing
python tests/test_coincidence.py     # twin-beam photon-pair statistics
python tests/test_fwm_generic.py     # generic SFWM topology + biphoton detector model
python tests/test_fwm_floquet.py     # N_F=1,2,3 direct/kernel parity + full-scan gate
python tests/test_fwm_pump_response.py # pump-only frame, trace-zero Nambu response + poles
python tests/test_fwm_angular_doppler.py # 2-D Maxwell geometry + grid/cutoff convergence
python tests/test_schemes_render.py  # every registered scheme computes + renders
```
(or `pytest tests/`)

## Adding a scheme

Subclass `gabes.schemes.base.Scheme` — declare `param_schema()`, `compute(params)`,
`observables(raw, params)` (and optionally `presets`, `info`, `extra_views`) — then
add an instance to the list in `gabes/schemes/__init__.py`. No UI edits: the
sidebar controls and the plots follow `param_schema()` and the observables dict.

## FWM conventions (must-know, not obvious from numbers)

### Seeded mean-field gain diagnostic

- Levels: g₁=F=2, g₂=F=3, e₂=F'=2, e₃=F'=3.
- OPD Δ (one-photon): ω_pump = ω(F=2→F'=3) + Δ.
- TPD δ (two-photon): ω_seed = ω_pump − ν_HF + δ.   ν_HF = 3.0357 GHz.
- Plot x-axis ref = **F=2→F'=3** line. (−) Raman branch = standard FWM seed, at Δ − ν_HF.
- Beam waists = **1/e² radius** (paper convention), not diameter. They are sidebar
  knobs (`pump_waist_um`, `probe_waist_um`) whose defaults are the paper geometry
  530 / 330 µm; `W_PUMP` / `W_PROBE` remain as those defaults. The pump waist sets
  the pump Rabi through I=2P/πw² **and** the Gaussian crossing overlap; the seed
  waist is nearly inert (its own Rabi is weak and the overlap is already ≈1), so it
  matters mainly for the downstream divergence λ/πw₀, which this scheme does not model.

### Generic SFWM / biphoton mode

- Mode selector: **Gain diagnostic** selects the 85Rb
  double-Λ model; **Biphoton** switches to the generic SFWM source
  estimate.
- Topologies: `cascade_rb87_telecom` (87Rb 5S1/2-5P3/2-4D5/2, 780/1529 nm),
  `cascade_cs_btw` (133Cs 852-917 nm or 852-795 nm BTW comparison), and
  `diamond_generic` (four-level user-wavelength template; not a validated paper
  preset).
- The default Biphoton sidebar is intentionally lab-facing: temperature, pump and
  coupling drive, pump detuning, two-photon detuning, collection geometry, filter
  bandwidth, and coincidence window. Detector calibration, manual wavelengths,
  source-model comparison, and numerical diagnostics live under Advanced.
- Biphoton readout: `g²_SI(τ)`, intrinsic source FWHM, timing-response-broadened
  detected FWHM, pair-rate estimate, singles, true and accidental coincidences,
  CAR, heralding estimates, and Cauchy-Schwarz R. The legacy `fwhm_ns` API key is
  retained as an alias for the detected width.
- The waveform is a coherent sum over Doppler velocity classes. Biphoton v3 phase
  matching uses calibrated longitudinal Δk plus absolute transverse Δk, with a
  strict vector `sinc²(|Δk| L / 2)` collection weight.
- Reference anchors: the 87Rb telecom source is calibrated to order
  `g²_SI≈44`, OD≈112, bandwidth≈300 MHz, and coincidence rate≈38,000 cps/mW;
  the Cs BTW preset exposes the wavelength-dependent temporal-width change
  reported for the 852-917 nm and 852-795 nm cascade channels.
- **Advanced source model toggle** (`biphoton_model`, default **Predictive**):
  - **Predictive** solves the Doppler-averaged cascade biphoton amplitude from
    first principles (Kim *et al.* QST 9, 045006 (2024) Eq. 2; Du, Wen, Rubin
    JOSAB 25, C98 (2008); Chen *et al.* PRR 4, 023132 (2022) Eq. 3-5): the
    two-photon denominator carries the **Ω_c² Autler-Townes term** (not a
    weak-coupling drive), the BTW is the collective velocity-class coherent sum
    with **natural-linewidth decay** (no injected lifetime), the source bandwidth
    comes from the waveform, and `g²_SI(τ)` is computed from `|ψ|²` with physical
    accidentals (no target-g² forcing). The wavelength-dependent BTW width
    **ordering** (917 narrower than 795) emerges.
  - **Calibrated** is the legacy reference-injected estimate (decay, bandwidth,
    g² target forced) kept for comparison.
- **Honest limits of Predictive** (documented in `info()` / the in-app validation
  table): absolute ns-widths and the Cs channel ratio are **approximate**; the
  wavelength ordering emerges, while exact per-source widths still need the
  deferred Rabi/dephasing calibration. At the default telecom point the modeled
  source FWHM is about 0.17 ns; the 0.55 ns net signal-idler timing response
  broadens the detected FWHM to about 0.50 ns, which is the quantity compared
  with the measured 0.56(4) ns. That agreement is not an independent validation.
  The absolute pair rate stays **reference-anchored** and currently
  scales with pump power, coupling drive, vector phase matching, and the square
  root of the density ratio; OD/cell length do not directly scale the rate. An OD
  waveform-reshaping path (Du/Chen ρ̄, group-delay/precursor) is implemented but
  **off by default** (`biphoton_od_reshaping`) since the lumped model overestimates
  it at high OD. Full quantum-Langevin noise is still future work
  (`docs/checklist.json`).

## Traps

Current FWM note: the 4-level FWM model applies the real 85Rb D1 hyperfine
`C_F^2` values `(10,35,35,28)/81` to Rabi couplings, polarization readout, and
spontaneous-emission branching. The macroscopic structural factor is
`fwm.physical_coupling_norm = 1/[2(2I+1)] = 1/12`; the trace-normalized density
matrix supplies the actual pump-modified manifold population exactly once.
Multiplying by an external equilibrium `p_F` would duplicate it. Density uses the
pure-85Rb CRC fit (`hyperfine.number_density`), consistent with that path. The
remaining seeded coupling is explicitly factorized under *Advanced* as
`reference residual × additional mode-overlap penalty × additional polarization
penalty × additional Zeeman-participation penalty`. The backward-compatible
reference residual is `0.74`; it has not been refitted after the normalization,
Maxwell-sign, mismatch-sign, and transit-reset corrections. The three lab-facing
penalties default to `1.0`; no unsupported numerical split of `0.74` is implied.
They represent one-sided extra coupling loss relative to the inherited setup,
not independently predicted overlap, Stokes purity, or population fractions.
The finite-reference-field atomic solve defaults to `N_F=3` and compares every
reported complex susceptibility/transfer coefficient, gain, wrapped phase, and
probe-gain optimum with `N_F=2` across the full displayed scan. A failed check is
returned as `UNCONVERGED`; trace/positivity alone never certifies truncation.
Separately, a slow reference for the standard (−) branch solves the self-consistent
pump-only state in an explicitly gauge-equivalent static frame and evaluates the
full infinitesimal 2×2 Nambu response on the trace-zero subspace at an independent
analysis frequency. It is an audit path, not the production default: production
still uses the finite-seed Floquet response. The inherited (+) branch fails the
  physical static-pump gauge-parity check and is rejected by this reference instead
  of being silently reinterpreted. An opt-in slow wrapper tensor-averages that
  branch over independent `(v_z,v_x)` Maxwell velocities with the lab beat and RF
  analysis frequency held fixed. At 121 °C and 0.32° it reproduces the analytic
  1.380 MHz Raman rms width and passes the documented grid/cutoff gates. Production
  remains one-dimensional, and neither reference supplies microscopic Langevin
  diffusion.
Ultra fidelity adds the slow propagation refinements, but the full
24-level Zeeman Floquet scan is still reported as a diagnostic rather than used
as the default full-scan solver.

1. **FWM gain is exponentially sensitive at high density.** At paper optimum
   T=121 °C linear Maxwell-Bloch still over-amplifies. The `0.74` reference
   residual and the three unit-default lab factors multiply the physical coupling
   on top of `1/[2(2I+1)]`; change a lab factor only when a corresponding
   effective-coupling measurement is available. `mode_overlap_penalty` is an
   additional unresolved transverse mode-matching penalty and is separate from
   Ultra's normalized axial crossing-angle profile. A post-hoc Manley–Rowe cap
   limits runaway power but is not a self-consistent depleted three-field solve.
   Absolute gain is therefore marked unsupported, and the dB curve is only a
   gain-referred algebraic diagnostic because atomic Langevin diffusion is absent.
2. **FWM Raman branches are separate mode pairs, not one summed susceptibility.**
   The Sim et al. 85Rb operating point uses the standard red-detuned seed on the
   (−) Raman branch. The production finite-seed path computes `branch=-1` and
   `branch=+1` independently; do not add their χ matrices into one 2×2 propagation
   matrix. The pump-only weak-response reference is certified only for
   `branch=-1`; `branch=+1` is explicitly unsupported until its inherited frame is
   corrected. The old branch-summed model created artificial high-gain extrema
   (for example near +70 MHz TPD).
3. **Rb is very absorbing — absorption schemes use short cells.** ls=1.0 is the
   true cross-section, so on-resonance OD saturates in a cm-scale cell; the
   OD/EIT/AT/CPT defaults use mm-scale cells / moderate T to keep features visible.
4. **Twin-beam coincidence is the ideal (lossless) parametric estimate** from the
   gain (n=G_s−1, g²_sc=2+1/n, R=g²_sc²/4>1), valid only as an ideal comparison;
   distributed atomic Langevin noise is not modelled.
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
   The Advanced **Paraffin-coated cell** checkbox replaces the immediate thermal
   reset with a two-ground-F reservoir: velocity is rethermalized between beam
   passages, while hyperfine population relaxes with the nominal `T1 = 25.1 ms`
   measured for one paraffin-coated ⁸⁷Rb cell at 300 K
   ([Bandi *et al.*, J. Appl. Phys. 111, 124906 (2012)](https://doi.org/10.1063/1.4729925)).
   The existing `γ_t` closes the MVP as both beam-exit and return cadence, so the
   coated result is semi-quantitative until cell/beam geometry and a cell-specific
   T1 are supplied. The switch changes atomic populations only: direct coating or
   window transmission remains exactly unity.

## Speed (why the architecture)

- L₀(Δ_eff) = L₀_base − Δ_eff·S_v (only the excited diagonal shifts with velocity).
  → all velocities stacked, one batched `np.linalg.solve`.
- The reduced χ̄ table is rebuilt for each temperature because the collisional
  broadening and thermal transit-reset dissipator depend on T and density.
  Within one such table, Δ shifts the effective-detuning axis and T sets the
  Maxwell velocity weights.
- The seeded app evaluates the compiled `N_F=3` response and a full-scan `N_F=2`
  comparison on every recompute. Runtime therefore depends on fidelity/grid size;
  the result is cached, and the TPD slider only navigates that cached curve.

## FWM future physics work

- Every non-legacy seeded tier uses Option-A bare frequency-specific wavevectors
  and vacuum/geometric mismatch. Ultra adds 64-segment propagation, Gaussian
  overlap, and approximate segment-wise pump-budget depletion.
- The separate pump-only reference includes a converged two-dimensional
  non-collinear Maxwell average. Promotion into the finite-seed production scan,
  beam-divergence averaging, and segmentwise pump-state recomputation remain future
  work.
- The diagonal Maxwell drift carries mean-field attenuation. Distributed in-cell
  vacuum/atomic covariance is unavailable and is not applied a second time.
  Pump scatter is an optional phenomenological diagnostic. SABES passes wanted
  sideband, residual carrier, and other-sideband powers end-to-end, but the latter
  components remain explicitly `unapplied/unsupported` until calibrated.
- The 24-level 85Rb D1 Zeeman manifold is built for CG diagnostics and correction
  bookkeeping; the full Zeeman Floquet solve remains future work because the
  density matrix jumps from 4-level `M=16` to 24-level `M=576`.

## Sim et al. 85Rb experimental reference (not a reproduced prediction)

G. Sim, H. Kim, H. S. Moon, Sci. Rep. **15**, 7727 (2025). 85Rb squeezing-optimal:

| Δ | δ | T | pump | seed | loss | → result |
|---|---|---|---|---|---|---|
| 0.9 GHz | −8 MHz | 121 °C | 600 mW | 8 µW | 5.5 % | gain ≈ 15, IDS −7.8 dB |

Geometry: cell L=12.5 mm, pump⊥probe. Pump/seed waist, crossing angle and detector
QE are knobs whose defaults describe this apparatus point. The current reduced
model does not reproduce its absolute gain or physical squeezing.

**QE caveat.** The scheme default is the inherited `qe_pct = 92 %`; it is not a
current squeezing calibration. The paper's detector — a PDB450A whose photodiodes are
swapped for Hamamatsu S3883 (0.58 A/W @ 795 nm, φ1.5 mm) — is **90.45 %**
(=1240·0.58/795), matching the 90.47 % quoted above. Using the true device QE lowers
the lossless floor `10·log10(1−η)` from −8.84 dB to −8.38 dB. Both are kept available
rather than silently re-anchored: 92 % is baseline-compatible, 90.45 % is the measured
device. These floors belong only to the ideal-vacuum algebraic completion and are
not physical bounds for the current atomic model.

## Deploy

- **Streamlit Community Cloud** (easiest), or Render/Railway/Fly/Cloud Run via Docker.
- Free tier = weak CPU → 6 s recompute slower. Real bottleneck = CPU, not host.
- Git remote: `github.com/Shake2313/fwm-squeezing-app` (public). gh acct `Shake2313`.
  Repo name ≠ folder name.
