# Grand Challenge code review — 2026-10-02

Scope: pending Grand Challenge campaign/provenance code, external controllers,
regressions, retained evidence and direct publication dependencies. Reviewed on
the existing `main` branch. Independent reviews covered campaign storage,
controller provenance, tests and publication dependencies.

## Findings and repairs

| Finding | Failure | Repair and regression |
|---|---|---|
| CR-1 / P2 / interrupted record publication | Direct final-path JSON writes could leave an incomplete immutable record that blocked resume. | Seal/serialize first, write and fsync a same-directory temporary, publish without overwrite. Tests cover interruption, I/O failure and existing/concurrent destinations. |
| CR-2 / P2 / consumed cache mutation | Cached numerical validation could precede later byte changes; separate native reads could validate different bytes. | Bind each native record seal to the captured raw snapshot; recheck consumed files before publication. Tests cover before/during/after-read mutation and A→B→A substitution. |
| CR-3 / P2 / controller first-capture gap | A just-created controller archive could change before its initial snapshot and become the accepted baseline. | Initial archive SHA-256 and length must match the captured controller descriptor; final byte checks remain. Six dependency-specific negative controls. |
| CR-4 / P2 / sealed report/request substitution | Hashing file A around a separate accessor could consume transient file B while preserving A's trusted file hash. | Parse and validate the exact captured bytes, including duplicate-key and seal checks. Grid reports, cross-comparison requests, parent proofs and transfers share this rule. |

The extension destination-race test now intercepts the extension's own
publication function, so the test cannot accidentally intercept unrelated
transactional JSON renames. Finalized source changes stopped before the full
suite. Earlier overlapping reviews correctly triggered frozen-source guards;
those transient failures are not the final candidate test result.

## Validation and commit scope

Final results: [test summary](test_summary.json), [scope and preservation](verification.json).
The full suite runs with `python -m pytest -q` in a clean export of the candidate
index. Unrelated working-tree physics/app changes cannot supply missing code.
After documentation receives the measured results, its focused consistency
checks run again; tested Python bytes must still match the final candidate.

Commit scope includes campaign caches/declarations/source ZIPs/controller
archives and the existing scientific-publication policy/manifest. Only direct
external document dependencies are included. Three historical website files are
ordinary byte snapshots; the embedded website repository and its metadata are
not changed or included. Narrow Git attributes preserve raw evidence across
Windows checkout. Historical whitespace is retained where hashes bind it.
The user's unrelated staged entries and working files remain separate.

## Scientific and historical limits

No new thermal pilot, source migration or physical certification is performed.
Frozen v2 remains 96 unique paths, 480 records and 4/9 audited grids; its mandatory
p2→p3/seed11 edge still fails. Official milestones remain **0/4**. The next
science question stays the bounded p3/seed811 pilot in the same
[blueprint](../blueprint.md).

The transactional writer repair changes current source and future capsules.
Existing immutable source ZIPs retain their original writer. Their records must
still pass native validation; an interrupted legacy write may require an
explicit recovery/source-migration decision. Never rewrite a frozen ZIP,
overwrite a corrupt final cache record, or silently relabel old evidence as
current-source output. Current external audit guards retain the original
numerical cache API and work with those historical capsules.

The earlier [governance audit](../audit_2026_10_02/README.md) and its 2091-pass
working-tree regression remain dated evidence from before these code repairs.
They are not the new clean-candidate regression or a renewed thermal audit.
