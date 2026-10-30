# Polarization squeezing through self-rotation

Request, 2026-10-30: implement PSR polarization squeezing as a separate scheme,
reuse existing atomic/polarization calculations, and state the approximations
and range of validity.

Choose **Polarization squeezing (PSR)** in GABES, or open the app with
`?scheme=psr`. **Scan axis** switches between a polarization-quadrature phase
scan at fixed optical detuning and minimum/maximum noise versus optical
detuning. Both plots show the 0 dB shot-noise level. The two phase traces are
orthogonal Stokes quadratures, with a 90° phase offset and 180° period; the
linear variance is sinusoidal, and its dB representation need not be.

The illustrative default is pure Rb87 D1, 70 °C, 75 mm, 0.1 mW/cm²,
1000 s⁻¹ isotropic ground-state reset, and −900 MHz detuning. Optical detuning
is relative to the line-strength centroid of the upper ground-F component.
After-cell transmission and detector efficiency default to 95% each. The
atomic excess-noise assumption defaults to zero and is labelled as omitted
on screen. This default is not an experimental calibration.

## ⟦AI⟧ Model implemented

The atomic part in [`gabes/psr.py`](../../gabes/psr.py) calls the existing
[`magneto_multilevel.multilevel_response`](../../gabes/magneto_multilevel.py).
It reuses Rb85/Rb87 D1/D2 level structures, signed dipoles, ground-state optical
pumping, vapor density, Jones polarization convention, and resonance-centered
Maxwell averaging. The new scheme does not modify the shared atomic solver.

Five input polarizations, 0 and ±h, ±h/2, give a local finite-difference slope
at the linearly polarized state. With circular amplitudes
`E± = (Ex ∓ i Ey)/sqrt(2)`, propagation `dE±/dz = i k χ± E±/2`, and
`S3 = 2 Im(Ex* Ey)`, the local Jones derivative determines the polarization
azimuth ψ and ellipticity angle ε. In particular,

\[
g=\left.\frac{\partial(d\psi/dz)}{\partial\epsilon}\right|_{\epsilon=0},
\qquad \alpha=-\frac{1}{S_0}\frac{dS_0}{dz},
\qquad \kappa=gL,\quad a=\alpha L.
\]

Both angles are radians; g and α are m⁻¹. The shared QWP convention has
ε = −β near linear polarization. Differentiating the Stokes angles avoids
treating QWP angle, circular intensity imbalance and ellipticity as identical.

For a bright linear carrier of real amplitude A and weak orthogonal mode b,
define `q=(b+b†)/sqrt(2)`, `p=(b−b†)/(i sqrt(2))`. Then
`ε ≈ p/(sqrt(2)A)` and `ψ ≈ q/(sqrt(2)A)`. The effective PSR mapping is
`q → q + κp`, `p → p`, with vacuum covariance I/2. Thus κ is gL, not 2gL.

[`gabes/psr_noise.py`](../../gabes/psr_noise.py) reuses the existing
[`GaussianChannel`](../../gabes/quantum/channels.py) for exact integration of
constant shear and common absorption. At normalized distance t ∈ [0,1],

\[
\frac{dV}{dt}=MV+VM^T+D,\qquad
M=\begin{pmatrix}-a/2&\kappa\\0&-a/2\end{pmatrix},\qquad
D=\frac{a(1+n_{\rm ex})}{2}I.
\]

`n_ex ≥ 0` is the **assumed atomic excess noise** per absorption optical depth.
It is not obtained from the susceptibility. The default `n_ex=0` keeps vacuum
loss but omits atomic excess noise. This is an optimistic conditional model.
The detector channel then applies `η = transmission_after_cell × QE` once:
`Vdet = η Vcell + (1−η) I/2`. The bright carrier transmission is `exp(−a)`.

The plot reports `10 log10(2 uθᵀ Vdet uθ)`, with `uθ=(cos θ,sin θ)`.
For the linear bright carrier this equals the linearized polarization variance
`Var(S2 cos θ + S3 sin θ)/|mean S1|`. Negative dB indicates squeezing in this
conditional model. Zero input intensity, zero detected throughput or zero
transmitted carrier makes the measured Stokes readout unavailable; it does
not display an invented measured vacuum level.

## ⟦AI⟧ Applicability and omitted physics

- Optical excitation is adiabatically eliminated in the weak-excitation
  approximation; ground optical pumping and Zeeman coherences remain. Ground
  hyperfine coherences are secularly omitted. Full optical saturation is absent.
- The magnetic field is zero. Isotropic transit reset and independent velocity
  classes describe an uncoated cell without explicit buffer-gas/wall dynamics.
  The reset control is s⁻¹, with no 2π conversion.
- The local pump and coefficients are uniform along the cell. Longitudinal pump
  depletion and multimode spatial evolution are absent. Pump absorption is used
  as a common-loss proxy for the orthogonal vacuum; its true polarization-dependent
  absorption is not derived. Significant ellipticity drift violates the pure-shear
  approximation.
- The static atomic slope is used for one effective, quasi-static detected mode.
  No finite-RF transfer, noise bandwidth, mode duration or MHz squeezing spectrum
  is predicted. Linearized Stokes normalization assumes a sufficiently bright
  carrier relative to the orthogonal-mode fluctuations.
- Atomic Langevin diffusion, fluorescence, spin noise, technical laser noise
  and their correlations are not inferred from classical rotation. The optional
  isotropic excess diffusion is a sensitivity parameter, not a microscopic or
  calibrated noise spectrum. A positive channel/uncertainty audit checks the
  assumed model's consistency, not the experiment's squeezing level.

Per-point engineering screens are virtual excitation bound ≤ 0.05,
`|αL| ≤ 0.3`, `|L ∂(dε/dz)/∂ε| ≤ 0.1`, and agreement between half/full derivative
steps within `10⁻⁶ m⁻¹ + 1% × |fine derivative|` for both rotation and drift.
These thresholds are not universal physical limits or experimental validation.
Estimates outside them remain visible, with a status warning and shaded detuning
regions. Numerical failures and negative absorption cannot feed the passive-loss
noise model and leave gaps instead.

Quick / Balanced / Fine use 41 / 81 / 161 optical points and velocity quadrature
order 4 / 6 / 10 **per integration interval**. Velocity convergence is not
certified automatically. Compare resolutions for each operating regime.
The selected detuning uses a direct atomic solve, reusing an exact scan sample
when available and otherwise caching one additional point. The cache follows
the scan's frozen atomic conditions, including isotope density and numerical
settings; changing detector/noise/view controls reuses that atomic point.
Its own health diagnostics and applicability screens determine its status,
independently of neighboring scan points. No extrapolation outside the scan is
performed. The detuning plot inserts the directly evaluated point so the curve
and displayed minimum/maximum agree; a failed direct point creates a gap.
The background spectrum still joins discrete samples and can miss narrow
optical features between them. Resolution controls the sampled spectrum and
velocity quadrature; it does not certify continuous spectral interpolation.

Result JSON includes atomic conditions, susceptibility stencil, health diagnostics,
scope and thresholds, `psr_readout.selected_point` with its atomic diagnostics,
both phase traces, the
conditional covariance and channel matrices, and detuning extrema. Curves can
also be exported as CSV through the existing export panel.

## ⟦AI⟧ Implementation validation — 2026-10-30

Review correction: the former selected-point interpolation failed for Rb85 D1,
Doppler off, 20 °C, 100 mm, 0.1 mW/cm², reset 1000 s⁻¹, +200 MHz. Direct selected
solves now give shear −0.350903, minimum −1.237066 dB and minimum-noise phase
49.852° at Quick, Balanced and Fine (external efficiency 0.9025). Focused
regressions compare the direct response, verify independent point validity and
failure gaps, and check atomic-cache reuse across detection changes.

The final isolated commit candidate was checked with:

```text
python -m pytest -q tests/test_fwm_quadrature.py tests/test_fwm_phase_integration.py tests/test_fwm_phase_ui.py tests/test_psr.py tests/test_psr_noise.py tests/test_psr_scheme.py tests/test_squeezing_mode_selection.py tests/quantum/test_contracts_and_channels.py tests/test_fwm_measured_gain_integration.py tests/test_fwm_passive_l1_integration.py tests/test_fwm_experiment_integration.py tests/test_jasperse_integration.py tests/test_plot_focus.py tests/test_liveplot.py tests/test_ui_defaults.py tests/test_defaults_store.py tests/test_ui_atom_isotope.py tests/test_ui_format.py --durations=10
```

**313 passed** in 50.06 s. Coverage includes atomic symmetry/scaling and derivative
refinement, four isotope/line combinations, an independent covariance ODE,
channel complete positivity and uncertainty, external/excess-noise limits,
absent carrier, invalid-point handling, both figures, direct selected-point
calculations/cache isolation, mode deselection, saved defaults and JSON.
After the final defaults-before-controls ordering fix, the affected suites
`test_squeezing_mode_selection.py`, `test_ui_defaults.py` and
`test_defaults_store.py` were rerun: **70 passed** in 4.27 s. Independent reviews
of the numerical/cache change, UI lifecycle and isolated commit scope found no
remaining actionable defects. Existing export `datetime.utcnow()` deprecation
warnings remain.
Shared atom/quantum solvers were not changed; full regression and unrelated
Grand Challenge campaigns were intentionally omitted.

Edge browser checks: phase/detuning switching, zero-dB references, shaded
approximation range, default noise omission label, and an excess-noise change
from 0 to 2 per OD (minimum moved from about −1.38 to +0.41 dB). Moving the
detuning to −500 MHz exposed the outside-screen warning. No browser errors;
existing Streamlit theme warnings remained. Screenshots:
[phase](../../output/playwright/psr_phase_20261030.png),
[detuning](../../output/playwright/psr_detuning_20261030.png).

Independent bounded velocity-order spot checks found order 4→8 changes of at
most 3.1×10⁻⁵ relative in g in the default Rb87 D1 wing; a wider four-isotope/line
spot check reached 0.40%. This supports the illustrative default, without
certifying velocity convergence throughout the controls' ranges.

## References

- [Matsko et al., *Vacuum squeezing in atomic media via self-rotation*, PRA 66,
  043815 (2002)](https://arxiv.org/abs/quant-ph/0112072): PSR mechanism and its
  connection to vacuum squeezing.
- [Lezama et al., *Numerical investigation of the quantum fluctuations of optical
  fields transmitted through an atomic medium*, PRA 77, 013806
  (2008)](https://arxiv.org/abs/0710.0609): multilevel quantum-noise calculation,
  including atomic noise missing from a mean-field rotation curve.
- [Mikhailov et al., *Vacuum squeezing via polarization self-rotation and excess
  noise in hot Rb vapors* (2009)](https://arxiv.org/abs/0903.3156): experimental
  and theoretical role of excess noise and the atomic level structure.
