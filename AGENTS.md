# AGENTS.md

Instructions for every coding agent here (Codex, Claude Code, others); the single
source. `CLAUDE.md` only imports it; `docs/grand_challenge/AGENTS.md` adds GC rules.

## Where things live

- **Basics, architecture, physics traps, references** → `README.md`. Read it first.
- **Front-end design refactor** → `docs/ui_redesign/README.md`: plan, decisions,
  dev log, before/after metrics from `tools/ui_audit.py`. Read it before touching
  layout or styling in `streamlit_app.py` or `gabes_ui/`.
- **Deferred or planned work and user decisions** → `docs/checklist.json`.
  Anything intentionally left for later, and every pending or recorded value
  decision, lives there, not in code TODOs.
- **Grand Challenge research turns** → `docs/grand_challenge/AGENTS.md`, then
  `.agents/skills/grand-challenge-head-manager/SKILL.md` and its durable
  hypermemory, then `blueprint.md` for the living master design/current card and
  `NEXT_SESSION.md`. See **Grand Challenge Head Manager** below.
- **Lean in Grand Challenge** → `docs/grand_challenge/lean_workflow.md`: an evolving
  working default for selected formalization tasks. Propose better equivalent
  methods to the user, revise the guide within the authorized scope, and append its
  bottom revision history. Lean is not a blanket requirement for GC work.
- **Scientific publication policy** → `docs/grand_challenge/publication_policy.md`;
  current artifacts → `docs/grand_challenge/current_publications.json`.
  Squeezing report owns experiment agreement/frontier exploration; analytic
  reconstruction owns theory/precision; quotient structure owns the explicit
  logical and order-theoretic organization. Each latest document must stand
  alone. Preserve substantive content when revising; document every
  omission/correction in Grand Challenge coverage records. Development logs and
  version differences belong only in `docs/grand_challenge/`, not the three
  scientific documents. This policy continues after Grand Challenge completion.

## Quick mental model

- `gabes/` is a generic atomic Bloch-equation engine; `streamlit_app.py` with
  `gabes_ui/` is a scheme-driven front-end that knows **no** specific physics. It
  renders whatever controls a `Scheme` declares (`param_schema`) and draws
  whatever it returns (`observables`). Adding a scheme = one entry in
  `gabes/schemes/__init__.py`.
- **Two-tier compute**: the heavy solve is cached on a scheme's `recompute=True`
  knobs only; `recompute=False` knobs are navigate-only and update the readout
  instantly. Put a knob in the tier matching whether it enters the solve or only
  the post-processing.

## Unified engine (development principle)

User policy, 2026-11-04: GABES's original purpose is to develop the many
calculations of atomic optics in one integrated system, so that the code is
consistent. This is the main development direction. See
[README → Development principle](README.md#development-principle).

### ⟦AI⟧ Execution rules implementing that request

- One implementation per piece of physics or numerics, in a shared `gabes/`
  module that schemes call: atomic data (`species`, `angular`); level structures
  and dissipators (`species.level_structure`, `core.dissipator`,
  `AtomModel.with_reset`); steady-state and Floquet solves (`core`); beam power
  to intensity, field and Rabi frequency (`beam`); cell temperatures, cold-spot
  density and axial Beer–Lambert propagation (`vapor_cell`); velocity averages
  (`doppler`, `affine_doppler`, `pole_doppler`); physical observables
  (`observables`); interpolated widths (`lineshape.level_crossing`); derived
  tables (`report`); photodetection and shot noise (`detection`, with
  `constants.PLANCK`).
- Before adding a helper, look for the shared one. If a scheme needs a variant,
  give the shared code an explicit parameter instead of copying it. A
  scheme-local copy of shared physics is a defect to consolidate.
- Consolidate in two steps.
  - Refactor first, keeping results. Snapshot every caller's raw outputs and
    displayed readouts before and after: recursively, including extra views,
    at the defaults and at the recompute slider ends. The refactor lands when
    the outputs are bit-identical. User decision 2026-11-04: it may also land
    when every caller differs by at most 1e-12 relative and every displayed
    readout is identical; record the difference in the commit and the checklist.
  - Value changes come second. Choosing between duplicated values, conventions
    or measures changes results. It needs a user decision recorded in
    `docs/checklist.json`, lands in its own commit, and recaptures the affected
    baselines (for example `tests/baseline_focused.npz`).
- When copies disagree, do not pick a winner silently. Pin each caller's
  current value, record the split and its size, and ask the user.
- Sealed records stay intact. Do not edit sealed scripts or receipts to route
  them through shared code; new scripts use the shared code from the start.
- Consolidate files that another session is editing only after that work lands.
- A new scheme uses the shared modules from the start and matches the gold
  references: FWM squeezing for physics realism, OD/SAS for cross-scheme
  uniformity (naming, units, structure).

## Grand Challenge and practical FWM

User policy, 2026-10-22: Grand Challenge and GABES FWM have distinct purposes.
Source: [practical FWM implementation request](docs/fwm_experiment_2026_10_22/README.md).

- Grand Challenge pursues fundamental analytic squeezing and retains its own
  registered acceptance criteria.
- Practical FWM serves the user's experiment. Approximate models, empirical
  coefficients, affine gain corrections and extrapolation are permitted tools.
  Grand Challenge no-fit or microscopic-completeness requirements must not block
  their implementation or use.
- A requested gain-calibration toggle must continue the calibrated model outside
  the measured interval and across operating controls; identify extrapolation
  and changed calibration conditions without silently disabling the correction.
  A bounded measured-reference interpolation alone does not satisfy this request.
- Keep measurement source, correction formula and calibration/extrapolation
  scope explicit. Empirical agreement is useful FWM evidence; it does not by
  itself establish Grand Challenge milestones or a noise covariance.

## Change-scoped validation

User policy, 2026-10-23: reduce costly regression runs for minor patches,
omit unrelated GC checks, and consolidate validation across development turns.
This replaces the unconditional full-suite requirement in earlier instructions.

### ⟦AI⟧ Execution rules implementing that request

- Choose checks from the changed behavior and its callers/shared dependencies,
  not the chat title or the folder containing the work log. Briefly state scope.
- Documentation, comments, labels and styling without behavioral changes need
  relevant link/content/render checks, not the numerical regression suite.
- Localized behavior changes need focused tests for that behavior and affected
  integrations. A minor patch may be completed and committed with those checks;
  a full-suite run is not required merely because code changed.
- Run `python -m pytest -q --durations=20` for broad shared-engine changes,
  cross-scheme numerical changes, test-discovery/dependency changes with broad
  impact, an unresolved cross-cutting regression, or an explicit full-validation
  request/release requirement. A commit alone does not require a full run.
- `tests/test_regression.py` pins the FWM compatibility baseline
  (`tests/baseline_focused.npz`); it is not required for every patch.
- Across a multi-turn work package, use focused checks during development. When
  full validation is warranted, run it once after edits and review fixes settle,
  before the final commit or completion claim. Intermediate progress need not
  trigger it. Reuse a passed run while relevant code, tests, data and environment
  remain unchanged; later localized edits need only affected checks unless they
  invalidate the broader result. State what each result covers.
- GC suites/audits and expensive research campaigns are separate scopes. Omit
  them when their model, shared engine, inputs, validators and evidence contracts
  are unaffected. Changes to shared solvers may require GC checks even outside
  GC folders. Do not treat all `tests/quantum` cases as unrelated by default.
  Existing registered scientific acceptance criteria remain intact when claiming
  their milestones; routine maintenance does not rerun sealed campaigns.
- Report commands, outcomes and intentional omissions briefly. Focused success
  is not a full-suite PASS. Record deferred required validation before handoff.
  Do not restart a long run solely for unrelated prose edits or duplicate a run
  already covering the same state. Review cost/progress before another long run.

## Conventions

- Instructions and reviews between agents: English. Final user reports and summary
  Markdown: Korean, `caveman` full style. Drop filler; keep equations,
  conditions, errors and validation scope. Existing records need no bulk compression.
- UI labels, units and markdown tables use unicode (°C, µT, µW, Γ, Ω, 2π).
  Matplotlib axis strings stay ASCII (mathtext layout lock; see the plot lock in
  `streamlit_app.py`).

## Internal Text Provenance

User policy, 2026-10-13: keep supplementary verification text brief. Separate
direct user requests and registered project goals from AI-added interpretations,
recommendations and validation notes. Prefix those additions with `⟦AI⟧`, or
place them in one section whose heading has that marker. Preserve the marker,
source and adoption status when summarizing or transferring them into memory.
An unmarked historical sentence does not establish human authorship.

Keep essential mathematical assumptions, units and result scope with the claim.
Consolidate repeated limitations once. This is an internal development convention;
frontend, user guides and public Docs use ordinary reader-facing prose.
Details: [internal provenance policy](docs/grand_challenge/internal_provenance.md).

## Text encoding on Korean Windows

PowerShell's default code page here is CP949, not UTF-8, so non-ASCII text breaks.

- Edit text with Python, not PowerShell, and pass `encoding="utf-8"` to every
  read and write.
- Do not pass Korean or other non-ASCII text as a shell argument or an inline
  string. Write it to a file first and read that file from Python.
- Edit Office files (PPTX/DOCX/XLSX) only with python-pptx, python-docx or
  openpyxl, never through shell text pipes.
- Prefer bash to PowerShell. If PowerShell is unavoidable, first run
  `[Console]::OutputEncoding = [System.Text.Encoding]::UTF8; $OutputEncoding = [System.Text.Encoding]::UTF8; chcp 65001`.

## Branch Policy

- Work on the current branch by default.
- Do not create, switch to, or open new git/Codex branches unless the user
  explicitly asks for it, or there is a clear technical blocker that requires an
  isolated branch.
- If a temporary branch is truly necessary, explain why, keep it narrowly
  scoped, commit or stash the work before leaving it, and close/delete the
  branch once the work is merged or no longer needed.

## Grand Challenge Head Manager

For every Grand Challenge session, including work outside its documentation
folder, read and apply
[Grand Challenge Head Manager](.agents/skills/grand-challenge-head-manager/SKILL.md)
at start, relevant checkpoints and closure. Use its persistent linked memory
to check long-term direction, closed work, failed approaches and repetition.
Apply the latest user request; old AI proposals are working guidance, not new
approval requirements. The same blueprint remains the master design. End every
research turn with a principal-investigator review; update the same blueprint,
close bounded work packages, and state the single next question. Sealed failures
and official milestone criteria remain intact.
