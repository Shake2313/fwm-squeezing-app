# FWM thermal and excess-noise report

This folder contains the Korean report on coherent-mode contamination, hot-vapor excess noise, laser/electronic noise, and balanced SQL calibration in seeded atomic four-wave mixing.

## Reproduce the numerical artifacts

From the repository root:

```powershell
python docs\fwm_thermal_excess_noise_report\scripts\generate_analysis.py
```

The script regenerates four CSV files, `data/analysis_summary.json`, and the four PDF figures. It calls the current GABES/SABES code for vapor-density and damping scales, the deterministic EOM/filter power ledger, and the two-dimensional non-collinear Maxwell-averaged weak response. The JSON records strict finite/null values, runtime versions, hashes for the analysis script and every local Python source under `gabes/` and `sabes/` (a conservative import-closure superset), and an 18/24/32-order velocity-quadrature convergence check.

The generated frequency response is **not** a physical squeezing spectrum. GABES does not yet supply the microscopic frequency-dependent atomic Langevin diffusion matrix or the measured detector transfer functions required for that claim.

## Build the report

```powershell
xelatex report.tex
bibtex report
xelatex report.tex
xelatex report.tex
```

The document uses XeLaTeX, `kotex`, Malgun Gothic, and standard BibTeX. The final deliverable is `report.pdf`.

## Evidence conventions

- **Theorem:** exact counting statistics and symmetric covariance-matched post-loss algebra under stated assumptions; unequal loss is treated separately.
- **GABES/SABES calculation:** deterministic model output from the repository.
- **Phenomenology:** a sensitivity or explicitly uncalibrated input-ratio proxy, not a microscopic identification or Fano-factor measurement.
- **Measurement:** a value reported by a cited experiment.
- **Open:** a physical claim that still requires covariance, transfer-function, or calibration data.

## Internal experimental provenance for the 1.8–2 MHz feature

The report treats the narrow feature as a recurring, setup-specific technical-noise problem. It can locally erase squeezing or rise above the plotted SQL, but its coexistence with strong low-frequency squeezing does not identify it as the cause of a broadband excess-noise floor.

- `230630_심기성.pptx`, slides 34–36: seed/EOM-path noise and transfer into the FWM difference signal.
- `230728_심기성.pptx`, slides 5 and 22: EOM-added noise and residual-carrier dependence.
- `231020_심기성.pptx`, slides 14–15: 1.8–1.9 MHz structure and reduction after removing an unused laser/driver.
- `231110_심기성.pptx`, slides 14, 16–17, and 20: electrical-noise transfer and 100 kHz–4 MHz passive-filter design.
- `231208_심기성.pptx`, slide 8: explicit 1 MHz and 1.5–2 MHz noise features.
- `231222_심기성.pptx`, slides 8–11: laser/PZT/MOPA-side attribution and passive-filter comparisons.
- `240213_심기성.pptx`, slides 4 and 6–10: direct 6.9 ± 0.2 dB at 100 kHz, a persistent near-2 MHz line, and separate slow-light/PD-imbalance bandwidth hypotheses. (`240212_심기성.pptx` was not present.)
- [2026-05-28 experiment log](https://app.notion.com/p/36e6cba14fee81449db8ead16b5dbc72): approximately 2.4 MHz peak and an open noise-eater test; squeezing was not measured that day because alignment failed.
- [2026-07-13 experiment log](https://app.notion.com/p/39c6cba14fee81d293bfc138689db30e): current near-2 MHz narrow peak and pedestal, evidence audit, and proposed separation tests. The `20261713` label in the request was interpreted as this `20260713` entry.
