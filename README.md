# GABES — General Atomic Bloch Equation Solver

Name aspirational. Current impl = **85Rb D1 double-Λ four-wave-mixing** density-matrix
model. Output: seed/probe gain G_s + intensity-difference squeezing dB. Streamlit front.

## What

- `fwm_obe.py` — physics backend. Steady-state OBE, 4-level (|F=2⟩,|F=3⟩,|F'=2⟩,|F'=3⟩),
  3-mode sideband expansion, Doppler average → χ̄ matrix → Maxwell-Bloch 2×2 propagation
  → G_s, G_c, squeezing. CLI `main()` plots a spectrum.
- `streamlit_app.py` — interactive UI. Sliders + "⚡ Sim et al. 85Rb optimum" preset button.
- `requirements.txt` — streamlit, numpy, matplotlib.
- `references/` — Sim et al. 2025 paper PDF + Steck Rb85 data + OE stabilization paper.

## Run

```
streamlit run streamlit_app.py
```

**Trap**: `python streamlit_app.py` does NOT work (no server, prints "run it with streamlit run").
CLI backend: `python fwm_obe.py`.

## Conventions (must-know, not obvious from numbers)

- Levels: g₁=F=2, g₂=F=3, e₂=F'=2, e₃=F'=3.
- OPD Δ (one-photon): ω_pump = ω(F=2→F'=3) + Δ.
- TPD δ (two-photon): ω_seed = ω_pump − ν_HF + δ.   ν_HF = 3.0357 GHz.
- Plot x-axis ref = **F=2→F'=3** line. (−) Raman branch = standard FWM seed, at Δ − ν_HF.
- Beam waists W_PUMP=530 µm, W_PROBE=330 µm = **1/e² radius** (paper convention). Not diameter.

## Traps

1. **LINE_STRENGTH_FACTOR = calibration knob, NOT physical.** Effective |d|² rescaling
   (Clebsch-Gordan lumping). `fwm_obe.py` default = 1.0. App default = 0.05.
2. **Gain is exponentially sensitive at high density.** At paper optimum T=121°C
   (N≈1.3e19/m³) linear Maxwell-Bloch gives nonsense (G~1e25) at ls=1.0. Drop ls (~0.05)
   or lower T. Model has no pump depletion / saturation → knife-edge near resonance.
3. **Must sum BOTH Raman branches.** Single (−) branch changes result materially
   (off-resonant (+) tail sets background absorption the exponential gain feels).
   `compute_spectrum(..., branches=BRANCHES)`.

## Speed (why the architecture)

- L₀(Δ_eff) = L₀_base − Δ_eff·S_v (only excited diagonal shifts with velocity).
  → all velocities stacked, one batched `np.linalg.solve`. 1180 Python loop → 1 call.
- R(δ, Δ_eff) table is **T- and Δ-independent** (T → Maxwell weights only; Δ → axis offset).
- Full −8…12 GHz scan ~227 s. App focused window ~6 s/recompute, `st.cache_data`.
  TPD slider = instant (navigates cached curve, no recompute).

## Sim et al. 85Rb optimum (the ⚡ button)

G. Sim, H. Kim, H. S. Moon, Sci. Rep. **15**, 7727 (2025). 85Rb squeezing-optimal:

| Δ | δ | T | pump | seed | loss | → result |
|---|---|---|---|---|---|---|
| 0.9 GHz | −8 MHz | 121 °C | 600 mW | 8 µW | 5.5 % | gain ≈ 15, IDS −7.8 dB |

Fixed geom: cell L=12.5 mm, QE 90.47 %, responsivity 0.58 A/W @ 795 nm, pump⊥probe.

## Deploy

- **Vercel ✗** — serverless, no persistent server, Streamlit can't run.
- **Streamlit Community Cloud ✓** (easiest), or Render/Railway/Fly/Cloud Run via Docker.
- Free tier = weak CPU → 6 s recompute slower. Real bottleneck = CPU, not host.
- Git remote: `github.com/Shake2313/fwm-squeezing-app` (private). gh acct `Shake2313`.
  **Repo name ≠ folder name** — repo stays `fwm-squeezing-app` after folder → GABES.

## Rename safe

No absolute paths. `fwm_obe.py` uses `__file__`, app uses relative import. Folder rename
GABES does not break code or git (git tracks by remote URL, not folder name).
