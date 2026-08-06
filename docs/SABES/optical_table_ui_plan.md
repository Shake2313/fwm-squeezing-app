# SABES — optical-table UI plan

Plan for the next stretch of SABES work: an interactive optical-table view,
click-to-edit optics, and virtual lab instruments. Written after Phases 0–4
landed (commit `f3898a5`), when the sidebar had grown to 23 controls in 7 groups
and the parameter count itself became the usability problem.

Tracked in `docs/checklist.json` under the `sabes-*` ids.

## Where SABES is now

Working, headless-testable, and deployed behind `?app=sabes`:

- `sabes/beamstate.py` — spectrum + Gaussian mode (M²) + Jones polarization
- `sabes/devices.py` — closed-form device laws (EOM Bessel ladder, etalon Airy +
  leak floor, saturable amplifier, fiber coupling, telescopes, waveplates)
- `sabes/calibration.py` + `calibration_sim2025.json` — 42 coefficients, each
  carrying a provenance tag; 17 verified, 25 placeholders
- `sabes/beamline.py` — the Sim et al. source chain as an **ordered list**
- `sabes/detection.py` — pre-solve geometry, post-solve readout
- `sabes/bridge.py` — settings → GABES params → squeezing → detector
- `sabes_page.py` — three compute tiers (optics / solve / readout)

Validated: fed the paper's own geometry the pipeline reproduces a hand-driven
GABES run to 1e-6 dB at both Fast and Ultra fidelity.

## The problem this plan solves

Two things, and they are not the same thing.

1. **Legibility.** A wall of numeric inputs does not tell anyone which knob is
   which optic, or what is upstream of what. A picture of the table does.
2. **Parameter count.** 23 controls is too many to scan. This is *not* fixed by
   drawing the table — it is fixed by moving parameters into per-optic popups,
   which is Stage D. Expect no reduction from Stage C alone.

## Decisions taken

| Question | Decision |
|---|---|
| Does the layout drive physics, or only the drawing? | **Drives physics** |
| How are the two setup drawings presented? | **Two tabs** (part 1 / part 2) |
| What about optics the model does not have? | **Draw everything, only modelled optics are clickable** |
| Oscilloscope / SA trace realism | **Both, as switchable modes** |

## Corrections owed to the model

Reading `references/IDS setup part 1.png` turned up three ordering errors in
`sabes/beamline.py`, plus several elements the model has no representation for.
These are real discrepancies, not cosmetic ones, and should be fixed as part of
Stage A rather than carried into the drawing.

| Item | Current model | Drawing |
|---|---|---|
| Glan-Taylor | one, in the seed trim | **several** — one in part 1 after the three etalons, plus a GTP on each arm in part 2 ahead of the cell PBS. The modelled one is *incomplete*, not misplaced |
| HWP between etalons | absent | **present between etalon 1 and etalon 2** |
| Collimators | treated identically | pump **F230**, seed **F220** |

Not modelled at all: optical isolators (×2–3), the cylindrical lens at the MOPA
input (it shapes the amplifier beam and therefore feeds `ta_m2_out`), the **noise
eater** in the seed path (an intensity-noise servo — directly relevant to a
squeezing measurement), and the **Rb-natural reference arm** (cell + PBS + PD +
f=300 lens) used to lock the ECDL.

The reference arm is worth a note: it is exactly the thing SABES deliberately
declined to model from first principles (ECDL frequency from temperature /
current / PZT is not predictable through mode hops). Seeing it in the drawing
*justifies* that decision — the detuning is set by looking at a saturated
absorption signal, not by reading the laser's dials. It also suggests the first
instrument demo: put a photodiode and an oscilloscope on that arm and show the
SAS trace. GABES already has an SAS scheme, so the physics is free.

---

## Stage A — layout model

The foundation all three requested features share. `sabes/beamline.py` currently
has no notion of position; the drawing, the hit-testing and "what is the beam at
this point" all need one.

### The central design point: two coordinate systems

The setup drawings are **schematics, not scale drawings**. In part 2 the pump
telescope lenses must sit `f1 + f2 = 450 mm` apart, but nothing makes the pixel
distance agree. Using drawing coordinates as optical path lengths would silently
corrupt the physics.

So the layout carries two independent things:

```
display   x, y, angle for rendering.  Source: the setup drawings. Never enters physics.
path      per-segment optical path length [m], with a provenance tag.
          Source: design values, measurements, or nominal.
```

They are allowed to disagree — schematics always do. What matters is that the
disagreement is visible rather than hidden.

This also means **the table does not need to be surveyed with a tape measure.**
Only segments where the beam meaningfully diverges enter the physics, and there
are few of them:

| Segment | Value | Provenance |
|---|---|---|
| Collimator → L1 | `f1` | design (afocal arrangement) |
| Telescope L1 ↔ L2 | `f1 + f2` (550 / 450 mm) | design |
| Cell → D-shaped mirrors | 584 mm | measured (23 inch) |
| D-mirrors → iris → focusing lens | currently 0.9 m | **nominal — worth measuring** |
| Everything else | — | nominal |

### Shape

```
sabes/layout/
  table.py    Placement(device_id, x, y, angle), Segment(a, b, path_m, beam_id)
  parts.py    the part 1 and part 2 layouts as data
  trace.py    walk the graph, propagate BeamState along segments
  probe.py    beam_at(x, y) -> BeamState | None
```

`beam_at()` is the whole point: **one function feeds all five instruments.**

**Done when:** power, waist and squeezing are unchanged, propagation distances
come from the layout rather than from loose coefficients, and the three ordering
corrections above are applied.

## Stage B — canvas component spike

The only unvalidated piece in this plan, so it is cut out and proven first.

Streamlit cannot return click coordinates natively. The chosen route is a
**custom bidirectional component with hand-written vanilla JS and no npm build**:
`declare_component(path=...)` pointed at a static folder whose `index.html`
implements the postMessage protocol directly (`streamlit:componentReady`,
`streamlit:setFrameHeight`, `streamlit:setComponentValue`). Roughly 80 lines, no
dependency, and it deploys as ordinary package data.

**Spike scope:** three SVG rectangles plus empty space; clicking returns either
`{"kind": "optic", "id": ...}` or `{"kind": "point", "x": ..., "y": ...}`.
Must be verified **both locally and on the deployed app**.

**Fallback if it fails:** the `streamlit-image-coordinates` package. Costs a
dependency, loses hover feedback and SVG interactivity, and forces server-side
PNG rendering per interaction — Stage E survives, Stage D gets worse.

## Stage C — the optical table

Renderer driven by the Stage A layout. Two tabs, matching the drawings.

### Rendering conventions

**Beam thickness cannot be to scale.** A 530 µm beam on a 1 m table is 0.05 % of
the width — invisible. Thickness is therefore a deliberate visual encoding
(log-scaled, roughly µW→3 px through W→12 px) and the legend must say so.

Path identity keeps the drawings' own colour convention: orange for the
high-power pump, red dashed for the seed, blue for the conjugate.

Modelled optics are visually distinguished from decorative ones (isolators,
steering mirrors, beam blocks), so "clickable" and "in the physics" mean the same
thing and nobody has to guess.

## Stage D — click to edit

`@st.dialog` is native in the pinned Streamlit, so the popup itself is free. The
work is hit-testing and migrating parameters out of the sidebar.

| Optic popup | Parameters it takes over |
|---|---|
| ECDL | output power; temperature / current / PZT as recorded annotations |
| MOPA | amplifier current, input HWP |
| Split PBS HWP | pump / seed division |
| Polarizer | EOM input power |
| Fiber EOM | generator offset, RF drive, input HWP |
| **Etalon ×3** | **one detuning each** — currently the UI forces a shared value |
| Glan-Taylor + trim HWP/QWP | seed power trim |
| Telescope lenses ×4 | focal lengths (they are swappable, so they are settings) |
| Rb cell | temperature |
| D-shaped mirrors | separation [mm] → crossing angle |
| Iris, focusing lenses, photodiode | radius, focal length, defocus |

**Sidebar after the migration:** fidelity, Solve knobs, Reset, instrument
selector, part tabs. Roughly 23 → 5.

Side benefit: per-etalon detuning becomes possible. The model already accepts a
3-tuple; only the UI was collapsing it, and in the lab each stage has its own
knob and its own temperature.

## Stage E — virtual instruments

Kept free of any layout or UI dependency so the package is reusable elsewhere.
Input is a `BeamState` (or a photocurrent signal); output is a `Reading`.

```
sabes/instruments/
  base.py               Instrument protocol, Reading, HeadKind(optical | photocurrent)
  power_meter.py        total power, plus a per-spectral-line breakdown
  beam_profiler.py      w_x / w_y, M², 2-D intensity map
  wavemeter.py          spectral line list
  photodiode.py         BeamState -> PhotocurrentSignal (the head for the two below)
  oscilloscope.py       scan mode / time-series mode
  spectrum_analyzer.py  RF noise power vs frequency: SNL, squeezed, electronic floor
```

Optical heads are placed by clicking a point on the table; the oscilloscope and
spectrum analyser need a photodiode as their head, exactly as in the lab.

**Priority: the spectrum analyser and the wavemeter.** The SA is the instrument
that actually measures squeezing, so it reproduces the paper's Fig. 3/4. The
wavemeter, dropped either side of the etalon chain, shows the filtering in one
picture — which is the clearest single demonstration of why SABES exists.

### Honest limits

The model is **steady-state**. There is no stochastic time evolution, and no
quantum Langevin noise (that remains deferred — see `fwm-quantum-langevin-noise`
in the checklist). Consequently:

- **Oscilloscope, scan mode** — a swept-laser trace. This is genuinely what a
  scope shows in this experiment, GABES already computes it, and it carries no
  extra assumptions. Default.
- **Oscilloscope, time-series mode** — samples drawn from the modelled noise PSD.
  Statistically correct, visually convincing, but it re-expresses what is already
  known rather than deriving anything new.
- **Spectrum analyser** — constructed from the noise budget (shot noise,
  squeezed level, electronic floor). Reproduces the measured shape; it is **not**
  a Langevin simulation.

Every synthesised trace must carry that label on screen. A plot that looks like
instrument output invites being read as evidence.

---

## Sizing and order

| Stage | Size | Risk | State |
|---|---|---|---|
| A layout | medium | low — headless and testable | **done** |
| B spike | small | was the only unvalidated piece | **done — fallback not needed** |
| C table | large | medium — ~60 optics of placement data | next |
| D popups | medium | low — `st.dialog` is native | |
| E instruments | large | medium — representing the SA honestly | |

## What A and B settled

**B: the hand-written component works.** No npm, no lockfile, no build step —
`sabes/components/` declares the component over a static `frontend/index.html`
implementing the postMessage protocol in vanilla JS. Verified in the browser at
`?app=sabes&dev=canvas`: optic clicks return an id, empty clicks return a
coordinate **in spec units** (checked at three render scales including 0.42×,
where reading screen pixels would have been wrong by more than 2×), and a repeat
click on the same optic still re-fires thanks to a monotonic `seq` — Streamlit
reruns only when a component value *changes*, so without it the second click
would vanish. **Still outstanding: confirm on Streamlit Community Cloud.**

**A: the layout drives the physics.** `sabes/layout/` holds the two coordinate
systems, and `detection` now asks it for `cell → D-mirror` and `cell → lens`
rather than carrying `cell_to_dmirror_m` and `pd_lens_distance_m` as
coefficients — both deleted. `beam_at(layout, x, y, chain)` returns the beam at a
point with its mode propagated *to* that point. Results are unchanged; the bridge
now reproduces a hand-driven GABES run to exactly zero at both fidelities.

Two bugs surfaced while building it, both of the kind that would have been very
hard to see later:

- Chain stages were keyed by name alone. Both arms have a stage called
  `"at cell"`, so a probe on the pump silently reported the **seed's** power.
- `Telescope.apply` collapsed the afocal pair into one step, so a probe between
  L1 and L2 got the state from *after* the telescope. The beamline now records
  the intermediate.

Still nominal and worth a tape measure: **D-mirror → iris (100 mm)** and
**iris → lens (216 mm)**. Together they set the beam size at the focusing lens
and therefore the photodiode spot, which is what the 4 W/cm² linearity finding
rests on.

## Open questions

1. **D-mirrors → focusing lens distance.** Currently nominal 0.9 m. It sets the
   beam size at the lens and therefore the focal spot, so it drives the
   photodiode power-density conclusion (the as-built f=75/100 mm lenses come out
   at 9.6 W/cm² against a 4 W/cm² linearity limit).
2. **Glan-Taylor position.** The drawing puts it after the etalons and before the
   delivery fiber; the model has it in part 2 after the telescope.
3. **Noise eater** — where in the seed path, and is it in use? Determines whether
   it is decorative or clickable.

## Related deferred work

- `fwm-quantum-langevin-noise` — the proper home for real noise propagation; the
  SA's synthesised traces are a stand-in until it exists.
- Residual pump *noise* (as opposed to power) is not modelled. GABES has a
  pump-scatter coefficient (`kappa`) at Ultra fidelity that is where it belongs.
- 25 of 42 calibration coefficients are placeholders. The three that matter most:
  the EOM `V_π`, the amplifier P-I curve, and the detector NEP after the S3883
  photodiode swap.
