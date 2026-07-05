# AGENTS.md

Guidance for Codex and other coding agents working in this repository.

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
