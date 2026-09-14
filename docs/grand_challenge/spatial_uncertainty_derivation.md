# Conditional input uncertainty for the stationary spatial field

2026-09-11. `analysis/grand_challenge/spatial_uncertainty_audit.py` connects the
stationary spatial mean/quantum-field solver to the correlated scalar-input
uncertainty machinery. It computes joint uncertainty in two gains and three
detected intensity-difference noise ratios. This is a conditional numerical
audit, with no absolute hot-vapor or experimental-validation claim.

The underlying spatial model and its nominal controls are recorded in
[the stationary spatial field report](s1_spatial_field_report_v2.json).
The uncertainty contract is derived in
[uncertainty_derivation.md](uncertainty_derivation.md); RF band status follows
[spectrum_analysis_derivation.md](spectrum_analysis_derivation.md).
The new outputs are
[spatial_uncertainty_report_v1.json](spatial_uncertainty_report_v1.json) and
[spatial_uncertainty_v1.png](spatial_uncertainty_v1.png).

## Fixed model and two varied inputs

The audit reads the parent report's physical inputs and detector ledger. Its
nominal fixture is:

| Quantity | Value or convention |
|---|---|
| Pump power and Gaussian waist | 0.6 W; 530 µm, 1/e² intensity radius |
| Seed power | 8 µW |
| Density; cell length | 10¹⁸ m⁻³; 12.5 mm |
| One-/two-photon detuning | 2π × 0.9 GHz; −2π × 8 MHz |
| Transit/reset rate | 2π × 0.1 MHz |
| Pump, probe, conjugate angles | 0, +5, −4 mrad; common x-z plane |
| Wavevectors | `CarrierGeometry.vacuum_beams`, no fitted dispersion |
| Transverse profiles | Coextensive normalized top-hats, fixed 400 × 300 µm rectangle |
| Area; velocity | 1.2 × 10⁻⁷ m²; stationary centers, (0, 0, 0) m/s |
| Dipole convention | `uniform-zeeman-rms` |
| Detector RF samples | 0.1, 1 and 4 MHz |
| Detector | Both transmissions × QE 0.85; current responses 1; balance 1; electronics PSD 0 |

The explicitly varied input vector and assumed standard uncertainties are

```text
x = (pump_power_W, conjugate_angle_rad)
mu = (0.6 W, -0.004 rad)
sigma = (0.006 W, 0.00001 rad) = (6 mW, 10 microrad).
```

These values are declared diagnostic assumptions, not apparatus measurements.
The base correlation is zero. The same derivative is also evaluated with
correlations −0.6 and +0.6; these are predeclared sensitivity controls, not fits
to gain or squeezing. Both marginal `ParameterEvidence` records and the
`CorrelationEvidence` record have status `assumed`.

All other consumed inputs, detector settings and transverse profiles are held
fixed. Conditioning on them does not assign zero uncertainty to their unknown
experimental errors. The report retains the full consumed-input ledger and
separate evidence audits for the two selected inputs and all consumed inputs.
Neither audit can grant a complete independent-input/no-fit claim here.

## Rebuilding the physical calculation

`SpatialUncertaintyEvaluator.medium(values, order_x)` creates a new
`ReducedPowerInputs` with the requested pump power, then rebuilds
`CarrierGeometry.vacuum_beams` with the requested conjugate angle and fixed
probe angle. It also constructs the declared transverse quadrature and a new
`StationarySpatialMedium`.

Consequently, an angle perturbation changes the actual wavevectors, transverse
and longitudinal loop phase, and the common-z-plane coupling normalization.
It is not represented by an output penalty or by a separate scalar
`phase_mismatch_rad_m`; that scalar must remain zero. A pump perturbation
rebuilds the pump Hamiltonian and the subsequent atomic mean/response/noise
calculation. The pump remains classical and undepleted along the cell.

`evaluate(values, order_x=..., propagation_rtol=...)` runs
`stationary_spatial_cell` and returns the five dimensionless outputs

```text
y = (probe_power_gain, conjugate_power_gain,
     R_0.1MHz, R_1MHz, R_4MHz).
```

The ratios come from the quantum contribution of the detected
`IntensityDifferenceSpectrum`, with its own same-current SQL. Electronics is
zero in this fixture. Gains retain their separately computed probe and
conjugate output-power definitions; no ideal gain identity is imposed.

Every newly calculated point must pass the transfer audit, have finite channel
CP and output covariance-uncertainty eigenvalues no lower than −10⁻¹⁰, and
provide finite nonnegative outputs on the declared RF axis. Each point records
the wavevectors, loop wavevector, transfer audit, eigenvalues, physical PSD and
SQL, ODE diagnostics, RF band status and elapsed time. Nominal parent checks do
not substitute for these perturbed-point checks.

The evaluator caches duplicate requests only within its current run, keyed by
the two varied inputs, spatial order and quantum propagation tolerance. It does
not load cached nominal outputs from the parent artifact.

## Central differences and spatial refinement

The fixed numerical choices are mean harmonic order 4, response harmonic order
3, mean ODE `rtol=2e-10, atol=2e-12`, quantum propagation `rtol=2e-9`, and
transverse `Ny=1`. The audit runs under `core.blas_single_thread()`.

For each input, the physical-unit central difference is

```text
J_ai(h_i) = [y_a(mu + h_i e_i) - y_a(mu - h_i e_i)] / (2 h_i).

coarse h = (0.0006 W, 0.000001 rad)
fine   h = (0.0003 W, 0.0000005 rad).
```

The three cases are `Nx8_coarse_steps`, `Nx8_fine_steps` and
`Nx12_fine_steps`. A two-input Jacobian needs the nominal point and four
perturbed points. Sharing the Nx8 nominal gives nine Nx8 solves; the Nx12 case
adds five solves, for **14 unique spatial-cell evaluations**.

Two comparisons have distinct meanings:

1. `step_refinement_at_Nx8` compares coarse and fine differences at a fixed
   spatial quadrature. It tests local derivative-step sensitivity.
2. `quadrature_refinement_at_fine_steps` compares Nx8 and Nx12 with the same fine
   differences. It tests spatial discretization sensitivity of the derivatives
   and their propagated uncertainty.

The comparison coordinates are the uncertainty-weighted response
`S = J diag(sigma)` and the output standard uncertainties `u`. Both use the
predeclared componentwise acceptance rule

```text
|S_test - S_reference| <= 1e-8 + 1e-3 |S_reference|
|u_test - u_reference| <= 1e-8 + 1e-3 |u_reference|.
```

The report also stores the relative covariance change and each comparison's
maximum fraction of its tolerance. Separately, current and historical nominal
outputs are compared with absolute tolerances 10⁻⁷ for gains and 10⁻⁶ dB for the
three noise ratios. These deterministic changes are not added to input
covariance as independent random variances.

The parent report's harmonic-order and ODE refinements concern its nominal
fixtures. This audit records that inherited evidence, including its mean-ODE
trajectory and quantum-ODE comparisons, but does not repeat harmonic/ODE
refinement at every newly perturbed input. Passing the two comparisons above
therefore does not establish joint convergence in every numerical parameter.

## Joint covariance and signed correlation contributions

For the declared two-input correlation matrix `R_x`, the local first-order
output covariance is

```text
C_x = diag(sigma) R_x diag(sigma)
U_y = J C_x J^T = S R_x S^T
u_a = sqrt((U_y)_aa).
```

`run_linear_case` uses `propagate_first_order` with the actual cell callback and
explicit steps. `linear_summary` retains the nominal output, Jacobian in each
input's units, weighted response, full output covariance, standard
uncertainties, output correlations and zero-variance output indices. Thus the
three RF errors and the two gain errors are not treated as independent merely
because they have different output names.

Writing `s_P` and `s_theta` for the two columns of `S` gives

```text
U_y(rho) = s_P s_P^T + s_theta s_theta^T
           + rho (s_P s_theta^T + s_theta s_P^T).
```

The first two terms are stored as
`uncorrelated_covariance_contributions`. The signed cross term for `rho=+0.6`
is stored as `cross_covariance_at_rho_plus_0_6`; the −0.6 case reverses its sign.
Cross terms are not separate positive quantum-noise sources. They can increase
or reduce an output's variance depending on the signs of its sensitivities.

`correlation_cases_using_same_J` recomputes this algebra for correlations
−0.6, 0 and +0.6 using the same Nx12 fine-step Jacobian and nominal outputs.
No additional physical solve or parameter refit is required for that controlled
comparison. Correlation changes the assumed input errors, not the nominal
atomic state or microscopic Langevin diffusion.

These standard uncertainties are local first-order results. Derivative-step
convergence does not prove the accuracy of a nonlinear uncertainty distribution,
the equality of the nominal output and an ensemble mean, or a confidence or
coverage interval. No Gaussian ensemble or large-excursion bound is calculated.
The figure's error bars denote one propagated standard uncertainty only.

## Historical evidence and current source snapshot

`load_parent` requires the passed stationary spatial v2 schema, stationary
velocity, declared RF samples and dipole convention, and passing nominal Nx8
and Nx12 rows with the +5/−4 mrad geometry. The parent file remains immutable;
its path and SHA-256 are retained in `parent_artifact`.

Preparation found three parent-source mismatches: `gabes/kernels.py`,
`gabes/doppler.py` and `gabes/schemes/fwm.py`. The new report records the actual
changed-source list with old and current hashes rather than presenting the
parent as an identical-code result. Current nominal outputs are recalculated
even if a future execution finds no source mismatch.

The audit captures its own explicit dependency hashes before and after the
calculation. `source_stable_during_run` must pass before
`expected_controls_passed` can pass. Nominal agreement with a historical report
is a regression comparison, not a substitute for current-source provenance or
experimental validation.

## RF and physical limits

Only three positive RF points are evaluated. `SpectralTrace` and
`analyze_squeezing_bands` preserve their detected-plane provenance and domain
censoring. A below-SQL interval reaching the sampled endpoints has
`width_hz=None`; its `observed_span_hz` is not an intrinsic bandwidth. Three
below-SQL samples also cannot rule out an unresolved peak between them. The
audit computes no refined RF minimum or crossing-to-crossing bandwidth.

The model retains stationary centers, a fixed classical pump and two fixed
top-hat optical profiles. Moving-atom transport, aperture and inter-slice noise,
walkoff, diffraction, additional collection modes, pump depletion, full-atom
extensions and measured detector/input uncertainties remain outside this audit.
Internal quantum consistency and conditional input propagation do not establish
an absolute hot-vapor prediction or held-out experimental agreement.

## Execution and results

Use new output names for another execution; the command refuses to overwrite
an existing report or figure:

```powershell
python -m analysis.grand_challenge.spatial_uncertainty_audit --output NEW_SPATIAL_UNCERTAINTY.json --plot NEW_SPATIAL_UNCERTAINTY.png
python -m pytest -q tests/quantum/test_spatial_uncertainty_audit.py
```

The immutable v1 report passes all declared controls. All **14 unique cell
evaluations** completed, with summed cell runtime **747.33 s** on the recorded
environment; this is elapsed audit time under concurrent project work, not a
benchmark. Source stability passed and all **34 source hashes** were checked
again after completion. The generated figure was visually inspected.

The selected Nx12 fine-step result, with independent input errors (`rho=0`), is:

| Output | Nominal | Conditional standard uncertainty |
|---|---:|---:|
| Probe power gain | 1.027772739898 | 0.000608066903 |
| Conjugate power gain | 0.030651495471 | 0.000608699722 |
| R, 0.1 MHz | 0.954088169635 | 0.000927417583 |
| R, 1 MHz | 0.954099308513 | 0.000927269957 |
| R, 4 MHz | 0.954270324115 | 0.000925028868 |

At 1 MHz the two uncertainty-weighted response terms are **−0.000261303415**
from pump power and **+0.000889691013** from conjugate angle. Their opposite signs
explain why a positive correlation reduces the propagated error. The 1 MHz
standard uncertainty is **0.001067147944**, **0.000927269957** and
**0.000762138054** for correlations −0.6, 0 and +0.6 respectively. This sensitivity
comparison retains the same nominal output and Jacobian; no parameter is fitted
to a measured spectrum. It does not establish the actual apparatus correlation.

The derivative step comparison has maximum weighted-response change
**3.88998×10⁻¹⁰**, maximum standard-uncertainty change **3.47895×10⁻¹⁰**, and
relative covariance change **7.46086×10⁻⁷**. The largest component uses
**0.00043244** of its predeclared weighted-response tolerance.

The Nx8→Nx12 comparison at the fine steps has maximum weighted-response change
**9.12603×10⁻¹³**, maximum standard-uncertainty change **1.08900×10⁻¹²**, and
relative covariance change **9.59467×10⁻¹⁰**. The largest component uses
**4.38702×10⁻⁶** of its weighted-response tolerance. Both nominal grid outputs
equal the corresponding saved parent gain/noise values at the recorded floating
precision, despite the explicitly recorded three source mismatches.

Every new point passes its own transfer audit. The maximum recorded local
commutator relative residual is **6.47310×10⁻¹⁴**; the smallest channel CP
eigenvalue is **2.43190×10⁻⁴** and the smallest output covariance-uncertainty
eigenvalue is **4.26535×10⁻³**. All 14 detected three-point band analyses retain
`width_hz=None`. These controls do not validate an intrinsic bandwidth.

The dedicated adapter tests pass: **33 passed (2.36 s)** in the final targeted
run. They include independent analytic wavevector derivatives and joint
five-output covariance, an actual zero-density spatial solve, rejected solver
results and cache isolation, parent/source provenance, and immutable output
prechecks. The coordinating main task completed the final repository-wide
`python -m pytest -q`: **961 passed, 1 failed (216.03 s)**, including all 33
spatial-uncertainty tests. The sole failure is the pre-existing missing
`FWM_physics.tex` in the repository-visibility documentation test. The main
task also rechecked this report's controls and 34 source hashes and owns the
shared checklist/research-log integration. This task did not modify those
shared files or any existing solver.

Report SHA-256:
`58a57a572674a78834d391f399cde40080d69c72b0f3bb08b3948f3860285346`.
