# Rydberg cell-heating validation workflow

This directory contains a reproducible, config-driven bridge between the GABES
Rydberg-EIT/finite-IF electrometry model, clean experimental CSV files, and a Korean
XeLaTeX validation report. It does **not** parse, embed, or modify the source
presentation.

Run from the repository root:

```powershell
python -m analysis.rydberg_cell_heating.workflow --config analysis/rydberg_cell_heating/example_config.json
```

The default output directory is `analysis/rydberg_cell_heating/generated/` and
contains:

- `results.json`: model sweeps, imported rows, capability states, and optima.
- `source_manifest.json`: SHA-256 provenance for config, inputs, and code.
- `results_macros.tex`: generated values consumed by the report source.
- temperature-spectrum and metric figures in PNG/PDF.

## Evidence statuses

Every external or derived quantity must carry one of these tags:

- `MEASURED`: derived directly from an available raw measurement.
- `PPT`: manually transcribed from a slide; not raw-data verified.
- `REFERENCE`: taken from a cited paper or calibrated reference.
- `FITTED`: inferred by fitting measurements.
- `PREDICTED`: independently emitted by the configured model.
- `ASSUMED`: unmeasured analysis input.
- `PENDING`: unavailable or not yet wired.

The legacy static finite-difference discriminator has units of inverse MHz and is
never relabelled as a field sensitivity. The separate finite-IF Liouvillian path
reports `nV/cm/sqrt(Hz)` only after an RF dipole convention and detector/noise
chain are declared. The example enables that conditional path with a
`REFERENCE` ARC dipole and explicitly `ASSUMED` detector/geometry inputs; it does
not fit or inject the paper's sensitivity value.

## Clean CSV contracts

CSV files are UTF-8, comma-separated, and header based. Empty optional numeric
cells are allowed; non-finite values and guessed columns are not.

EIT spectrum (one trace file per configured temperature; `temperature_c` is
metadata in the JSON input entry):

```text
Probe detuning [MHz],Photodiode voltage [V]
-2.0,0.154
-1.9,0.159
```

RF sensitivity summary:

```text
temperature_c,sensitivity_nv_cm_sqrt_hz,uncertainty_nv_cm_sqrt_hz,psn_limit_nv_cm_sqrt_hz,status,note
40,6.68,2.55,,PPT,manual transcription pending raw-data verification
```

Temperature log:

```text
setpoint_c,sensor_left_c,sensor_right_c,effective_vapor_temp_c,cold_spot_c,elapsed_s,status,note
40,39.8,40.2,,,1800,MEASURED,cold spot not independently measured
```

Paths in the config are resolved relative to the config file. Set `required` to
`true` only when a missing file should abort the run.

Raw EIT, RF sweep, and PSD traces are parsed by
`gabes.rydberg_experimental_csv`, which preserves SHA-256 provenance, unit
assumptions, ignored-row diagnostics, and unmodified detector signals. Configure
RF inputs with `kind: "rf_sweep"` and PSD inputs with `kind: "psd"`; loader
keyword arguments go in the input entry's `loader` object.

The model also uses `gabes.rydberg_experiment` to keep heater setpoint, sensor,
effective vapor, and cold-spot temperatures separate. The resolved cold spot is
passed into the solver, so it controls vapor pressure rather than merely being
reported as metadata. When no temperature log is supplied, the identity
setpoint-to-vapor mapping is explicitly tagged `ASSUMED`. The optional effective
atom number is a geometric estimate whose participation and overlap inputs
retain their own `ASSUMED` status.

Set `model.axial_profile.enabled` to `true` to add a configured or sensor-derived
`T(z)`/`n(z)` profile, axial column density, axial `N_eff`, and segmented
Beer--Lambert integration. The current report path rescales the lumped OBE line
shape by the local density; its result is labelled as such and is not claimed to
be a full spatial OBE solve or a spatial finite-IF sensitivity. The pressure cold
spot is the colder of the declared reservoir and sampled `T(z)` minimum, and the
same value is passed to both the lumped solver and axial density model. The
underlying `integrate_beer_lambert` API also accepts independently solved local
absorption arrays.

## Absolute electrometry chain

The `electrometry` config block can reuse the scheme-native finite-IF response
(`mode: "scheme"`) or run a separately configured detector chain
(`mode: "configured"`). Both paths retain the RF transition dipole/angular
factor, photodiode shot noise, RIN/electronics, and PSN/technical/total field ASD
as explicit quantities. Neither uses a paper sensitivity as an anchor. This is
a first-order frequency-domain weak-SIG calculation; full time-domain LO+SIG
waveforms and lock-in filters remain deliberately deferred.

## SAM calibration and uncertainty

Set `sam_calibration.enabled` to `true` and supply source power [dBm], antenna
gain [dBi], distance [m], losses/correction, and their one-standard-deviation
uncertainties. The workflow returns RMS or peak field, propagated standard
uncertainty, delivered power, and `r/(2D^2/lambda)` with a warning when the
declared geometry is outside the usual far-field criterion. It is a standard
antenna far-field conversion, not a horn near-field or standing-wave model. For
AT points the workflow also converts peak/RMS conventions, checks an optional
frequency tolerance, and reports the SAM/model field ratio, residual, and
SAM-derived uncertainty. The model field remains tagged `PREDICTED` until an
experimental full-spectrum fit supplies the RF Rabi parameter.

## Optional helper protocol

The adapter discovers report-aware helper modules only when they declare:

```python
REPORT_ADAPTER_API = 1

def report_temperature_point(*, params, raw, readout, helper_config):
    return {
        "status": "PREDICTED",
        "field_sensitivity_nv_cm_sqrt_hz": 10.0,
        "psn_limit_nv_cm_sqrt_hz": 9.0,
        "superheterodyne_responsivity": 1.0,
    }
```

This protocol keeps future `gabes.rydberg_readout` and
`gabes.rydberg_experiment` APIs out of the CSV/report layer. Modules with a
different or absent API marker are recorded as `PENDING` rather than called.
