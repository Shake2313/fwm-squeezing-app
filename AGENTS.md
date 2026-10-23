# AGENTS.md

Guidance for Codex and other coding agents working in this repository.

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


## Branch Policy
- Work on the current branch by default.
- Do not create, switch to, or open new git/Codex branches unless the user
  explicitly asks for it, or there is a clear technical blocker that requires an
  isolated branch.
- If a temporary branch is truly necessary, explain why, keep it narrowly
  scoped, commit or stash the work before leaving it, and close/delete the
  branch once the work is merged or no longer needed.

## Project Orientation
- Read `CLAUDE.md` and `README.md` before making non-trivial changes.

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
