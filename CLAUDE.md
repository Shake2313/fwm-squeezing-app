# CLAUDE.md

Orientation for AI assistants working in this repo.

## Where things live
- **Basics, architecture, physics traps, references** → `README.md`. Read it first.
- **Deferred / planned work (not done yet on purpose)** → `docs/checklist.json`.
  Anything intentionally left for later (e.g. wiring Rydberg beam power to the OBE
  drive, a predictive biphoton model) is recorded there, not in code TODOs.

## Quick mental model
- `gabes/` is a generic atomic Bloch-equation engine; `streamlit_app.py` is a
  scheme-driven front-end that knows **no** specific physics — it renders whatever
  controls a `Scheme` declares (`param_schema`) and draws whatever it returns
  (`observables`). Adding a scheme = one entry in `gabes/schemes/__init__.py`.
- **Two-tier compute**: the heavy solve is cached on a scheme's `recompute=True`
  knobs only; `recompute=False` knobs are navigate-only and update the readout
  instantly (no re-solve). Put a knob in the right tier based on whether it enters
  the solve or only the post-processing.
- Schemes share engine pieces in `core.py` / `doppler.py` / `observables.py`
  (physical observables), `lineshape.py` (FWHM extractors), `report.py`
  (derived-table markdown). Reuse these instead of re-implementing per scheme.

## Conventions
- Gold references: **FWM Squeezing** for physics realism, **Absorption
  Spectroscopy (OD/SAS)** for cross-scheme uniformity. New schemes should match
  their naming/units/structure.
- UI labels, units and markdown tables use unicode (°C, µT, µW, Γ, Ω, 2π).
  Matplotlib axis strings stay ASCII (mathtext layout-lock guard — see
  `streamlit_app.py` plot lock).
- Run `python -m pytest -q` before declaring done; `test_regression.py` pins the
  FWM baseline (`tests/baseline_focused.npz`).
