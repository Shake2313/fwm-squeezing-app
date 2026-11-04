# GABES Lean pilot (not formally adopted)

One-session feasibility pilot requested 2026-10-28. Independent nested Lake project;
not a production dependency, CI gate, or GC milestone.

- Pinned Lean 4.35.0, core + Std only; no Mathlib or external packages.
- Exact rational scalar covariance maps: g is the power/covariance multiplier.
- Vacuum variance 1/2 (hbar=1); independent input/reservoirs and stage reservoirs.
- Physical pure loss additionally requires 0 <= eta <= 1.
- Six universally quantified algebraic identities, seven exact witnesses/counterexamples.
- No real-matrix, quantum CP, microscopic FWM, or Python equivalence claim.

Run from the repository root in PowerShell:

```powershell
.\analysis\grand_challenge\lean_pilot\setup.ps1
.\analysis\grand_challenge\lean_pilot\run_pilot.ps1
```

The portable toolchain lives under the user's temporary directory by default.
setup.ps1 accepts -InstallRoot; run_pilot.ps1 accepts the corresponding -LeanBin.
No global PATH change is made. Dependencies/build output are excluded locally.

run_pilot.ps1 requires a successful build, rechecks Main.lean, audits all 13 theorem
axiom lists, and requires semantic failure diagnostics for two deliberately false
theorems under negative/. These negative files are outside the default build.
Only propext, Classical.choice and Quot.sound are allowed. The positive cases use
grind or decide +kernel, never native evaluation.

The default result directory is docs/grand_challenge/lean_pilot_2026_10_28/run.
Use a fresh -OutputDir for a later run to preserve this pilot's evidence.
See the [pilot report](../../../docs/grand_challenge/lean_pilot_2026_10_28/README.md).

Official sources:
- https://github.com/leanprover/lean4/releases/tag/v4.35.0
- https://lean-lang.org/doc/reference/latest/Tactic-Proofs/Tactic-Reference
- https://lean-lang.org/doc/reference/latest/ValidatingProofs/
