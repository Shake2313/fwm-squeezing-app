# Next session

Governance snapshot: **2026-10-02**. Atomic science snapshot: **2026-09-28**.
Read the SAME [living blueprint](blueprint.md) first: dashboard, closed scopes,
open dependencies, current experiment card, mandatory end-turn PI review.
Goal: independent inputs -> atomic mean/noise/response -> optical channel ->
validated absolute hot-vapor FWM gain and S_minus. Official milestones **0/4**.

## Verified state and closures

- Conditional frozen v2 `[2,3,4]`, seeds `[11,211,811]`, 108-source ZIP:
  **96 unique paths / 480 records / 4-of-9 grids audited**, p2/all seeds+p3/11.
  Submitted path/native/arithmetic checks pass; ensemble convergence fails.
- P2->p3/11: greater **10.5141%**, lesser **7.8477%**, source-resolved
  **13.7886% / 12.3309%**, Poisson **3.6776%**, response **36.1694%**.
  Five of six exceed 5%; all six directed p2 scramble comparisons also fail.
- Existing `[2,3,4]` cannot pass: both mandatory edges are required and the first
  already fails. Preserve the negative decision. **Do not automatically run the
  remaining 960 solves to fill the declaration.**
- Blueprint accepts six bounded historical evidence packages and closes the
  v2 PASS-feasibility question as FAIL. These are not seven GC milestones or a
  project completion percentage; original four milestone criteria stay unchanged.
- Thermal ensemble, boundary history, nonlocal Maxwell, independent-input
  absolute squeezing and untouched experimental holdout remain uncertified.

## Next primary package: WP-THERM-PILOT

Preregistered question/card lives in blueprint; edit it in place when executing.
Choose **p3/seed811 only**: [frozen-source geometry probe](audit_2026_10_02/geometry_probe.json)
finds two >=4 us chords (indices20,45; max6.01058 us), versus none in p3/11.
Geometry estimates do not certify atomic-noise tails or true integration error.

1. Inspect current dirty tree. Use frozen numerical source/model/environment and
   exact-spec reuse; current unrelated changes must not enter old capsule.
2. [Execution contract](thermal_grid_execution.md): confirm120 reused+120 new
   requests and distinct outputs. Forecast wall/memory allowing for longer chords;
   record operational allocation before launch. Historical p3/11 wall10,730 s
   is a reference only; seed811 cost is unknown. No automatic queue.
3. Revalidate/rebind120 p2/811 records and perform at most **120 new solves**
   (24 new paths x5 specs). Preserve old evidence.
4. Audit batch/grid; compare p2->p3/811 and BOTH directed p3 seed11<->811 pairs.
   Diagnose all RF/source matrices and signed strata. Keep original path gates
   and all six 5% ensemble metrics. Candidate-SI-floor ensemble/diagnosis controls
   pilot routing; preserve standalone coarse-floor comparison separately.
5. All metrics pass -> candidate evidence for separate finer-window design;
   seed211 and final edge/scramble conditions remain missing. Any >5% ->
   PIVOT_REVIEW of quadrature strategy versus boundary pilot. Path failure or
   incomplete allocation -> INCONCLUSIVE; preserve checkpoint, no certification.
6. End with PI review. **No automatic p4/p5 extension**. A harder path needing extra
   refinement gets a new bounded card; never loosen gate to force completion.

Parallel small work allowed: missing-input/provenance inventory and a small
nonlocal optical-closure design. Boundary-history [plan](progress_2026_09_28/boundary_history_plan.md)
remains unexecuted; keep collection Q fixed and distinguish state-history totals
from noise-source ancestry. Independent data preparation need not wait for thermal PASS.

## Integrity and interpretation

- Native arrival weights, one density factor, mean outer Poisson term, all complex
  source/RF matrices and signed cancellation remain mandatory. Geometry occupancy
  or missing tails never justify renormalizing rates or deleting paths.
- Candidate SI floor in ensemble/diagnosis differs from coarse SI floor in
  standalone comparison. Do not silently merge conventions; retain both reports.
- Two consecutive turns without new discriminating evidence/dependency closure/
  counterexample require strategy review before a third similar expansion.
- [Document map](DOCUMENT_MAP.md) supplies relevant details. Avoid full log/cache/
  matrix loading at startup.
- [October2 audit](audit_2026_10_02/README.md): pre-code-review working-tree regression **2091 passed /
  3 skipped / 0 failed**,803.24 s. Prior2077-pass result remains dated evidence.
  Source-bound historical reports and current working-tree regression are separate.

At turn end: update blueprint's latest PI review/current card, append research
log, synchronize only relevant checklist state and this handoff. Do not fork a
parallel master design or alter official acceptance criteria.

## Code integrity follow-up — 2026-10-02

[Code review](code_review_2026_10_02/README.md): clean commit candidate
**2044 passed / 4 skipped / 0 failed**. Current writer publishes complete JSON without overwrite;
external audits bind consumed cache seals, archived controllers, reports and
trusted requests to captured bytes. Historical ZIPs retain their old writer:
an interrupted legacy cache write requires explicit recovery/source-migration;
never overwrite a bad final record or silently change the frozen source.
Scientific state and the p3/811 pilot card above stay unchanged.
