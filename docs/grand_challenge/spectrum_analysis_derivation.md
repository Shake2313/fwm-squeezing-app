# RF squeezing band extraction

2026-09-11. `gabes/quantum/spectrum_analysis.py` adds a downstream numerical
readout that can run independently of the spatial atom-to-field elimination.
It consumes the positive laboratory RF axis and the PSD conventions already
defined in [conventions.md](conventions.md) and
[readout_derivation.md](readout_derivation.md). It does not derive new atomic
physics or validate experimental squeezing.

## Observable and boundaries

The input is the linear ratio

```text
R(f) = S_difference(f) / S_SQL(f),       f > 0
t = 10^(threshold_db / 10),             threshold_db <= 0
B_t = { f in the declared RF domain : R(f) < t }.
```

`SpectralTrace` carries `AnalysisFrequencyAxis`, source/detected layer,
provenance, explicit sample validity and numerical convergence masks, and a
measurement-filter scope. The default scope is pointwise PSD with no RBW/VBW
convolution. Source and detected traces remain separate objects; the code does
not invert a detector transfer or infer an atomic bandwidth from detector roll-off.

`from_psds` requires matching one-sided A²/Hz numerator and SQL samples. Finite
numerators must be nonnegative and finite SQL values must be positive. Nonfinite
data remain missing. Matching units, current planes and SQL conventions are
declared by the caller; this adapter cannot verify an external measurement.
`from_intensity_difference` consumes the established quantum readout and
explicitly selects quantum-only or total (quantum plus electronics) PSD. It
defaults to the detected plane and never silently subtracts electronics.

Adjacent accepted samples straddling the threshold define the linear-ratio
crossing estimate

```text
f_cross = f0 + (t - R0) (f1 - f0) / (R1 - R0).
```

Interpolation is in linear PSD/SQL ratio, not dB. Exact equality samples and
plateaus are excluded from the strict subthreshold set. They therefore separate
components, including when below-threshold samples lie on both sides of an
isolated equality point. A technical-noise peak likewise separates two bands;
the extractor never spans it with one combined width.

Every edge is tagged `crossing`, `domain`, `invalid` or `unconverged`. Missing,
nonfinite, explicitly invalid and numerically rejected samples are never bridged
by interpolation. At a missing/invalid or unconverged neighbor, the reported
edge frequency is the last accepted subthreshold sample, **not an inferred
threshold crossing**. `width_hz` exists only when both edges are crossings;
`observed_span_hz` records the accepted span even when an edge is censored.
Thus a flat ideal spectrum below SQL has a domain-limited observed span and no
intrinsic finite bandwidth. Extending the scan changes the observed span.

## Numerical status and refinement

`analyze_squeezing_bands` reports every component, sampled minima and one of:

- `bands_found`: at least one sampled subthreshold component exists; consult
  `complete_domain` and invalid/unconverged counts for missing regions.
- `no_squeezing`: every supplied sample is accepted and none is below SQL when
  the selected threshold is 0 dB.
- `no_target_band`: every supplied sample is accepted but none reaches the
  selected negative-dB target; weaker squeezing may still exist.
- `incomplete`: no accepted sample is subthreshold but gaps prevent that
  conclusion for the full sampled domain.
- `unavailable`: no sample is accepted.

`no_squeezing` is always limited to the supplied resolution and domain. It does
not rule out an unresolved narrow band between samples. `complete_domain` means
that every *supplied sample* is accepted; it is not a claim of infinite RF
coverage or sufficient RF sampling.

An optional `crossing_evaluator(frequency_hz)` refines existing crossing
brackets by bisection. The provider must return a finite nonnegative linear
ratio and reproduce both supplied endpoint values (relative tolerance 1e-10,
absolute tolerance 1e-14) while preserving each endpoint's strict threshold side
or equality. This rejects a different source/provider rather than returning a
precise artificial crossing. Refinement requires an explicit tolerance in Hz.
The midpoint estimate, final bracket, bracket width and iteration count are
returned. For a continuous provider with a crossing in that bracket, the midpoint
location error is at most half the bracket width. Exhausting the iteration limit
reports `unresolved`; it does not erase the observed existence of the band.

Local refinement does not improve sampled minima, discover missed bands,
validate provider covariance, propagate calibration/model uncertainty, or assess
the convergence of Floquet, velocity or spatial calculations. Those statuses
remain separate. Independent RF-grid refinement is needed for global coverage.
Input uncertainty must be propagated by reevaluating the physical/measurement
provider and retaining changing band topology and censoring. Root-location
precision alone is not a prediction interval.

## Independent acceptance fixtures

The analytic valley

```text
R(f) = 0.25 + ((f - 2.1 MHz) / 1.3 MHz)^2
```

has SQL crossings `2.1 MHz ± 1.3 MHz sqrt(0.75)` and a 0.5-ratio target band
from 1.45 to 2.75 MHz. Tests compare coarse/fine interpolation and local brackets
against these closed-form values. Other fixtures cover split bands, equality
plateaus, missing/unconverged holes, domain censoring, unavailable data, and
provider inconsistency.

A shared scalar electronic transfer gives

```text
R_detected(f) = |H(f)|² S_difference(f) / (|H(f)|² S_SQL(f)) = R(f)
```

where the denominator is nonzero and there is no separately added electronics
term. The low-pass fixture checks this algebra at the same current plane and
preserves the source band edges. It makes no claim that unequal arm responses,
optical loss or additive detector electronics leave squeezing unchanged.

Reproduction:

```powershell
python -m pytest -q tests/quantum/test_spectrum_analysis.py
```

This completes a bounded extraction primitive for
`fwm-rf-squeezing-spectrum-and-bandwidth`. Absolute hot-vapor gain/noise
prediction, calibrated measurement filtering, uncertainty-aware experimental
band edges and held-out validation remain separate Grand Challenge requirements.
