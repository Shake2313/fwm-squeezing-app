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
- Run `python -m pytest -q` before declaring code changes complete.
