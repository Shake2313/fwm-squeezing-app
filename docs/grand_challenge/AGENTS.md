# Grand Challenge context rules

Follow repository `AGENTS.md` and `CLAUDE.md`. After repository orientation:

- Read [blueprint.md](blueprint.md) for the living master design, scoped work
  packages, active experiment card, and end-turn principal-investigator review.
  Then read [NEXT_SESSION.md](NEXT_SESSION.md) for the compact execution handoff.
- Human entry: [README.md](README.md). Topic lookup: [DOCUMENT_MAP.md](DOCUMENT_MAP.md).
- Load only task-relevant derivations, contract sections, or evidence fields. Long research log, cache trees, source ZIPs, and full JSON matrices are not default context.
- Check snapshot date. Old "latest" / "next" text describes its own checkpoint. Resolve conflicts against current checklist and sealed evidence; update handoff after verified progress.
- Preserve historical evidence bytes and paths: reports, cache records, declarations, source ZIPs, archived controllers, figures, logs, validators. New results need distinct outputs; failed gates stay recorded.
- Audit success, path accuracy, ensemble convergence, optical validity, and experimental validation are separate claims.
- Before publication work, read [publication_policy.md](publication_policy.md); current editions live in [current_publications.json](current_publications.json).
- Agent notes: compact English; user summaries: concise Korean. Preserve equations, assumptions, units, tolerances, counterexamples, and provenance. Compression must not change claims.
- Code changes: run repository-required `python -m pytest -q`. Documentation-only navigation: verify links, status sources, and preservation.

## Research-turn contract

- Advance one primary work-package question per research turn. Freeze hypothesis,
  scope, gate, output, new/reused solve forecast, and stopping/pivot rule before
  observing results. Long derivations/batches may span turns; record verified
  partial progress and unresolved conditions without claiming completion.
- Before ending each research turn, perform the principal-investigator review in
  the blueprint: new claim/evidence, PASS/FAIL/INCONCLUSIVE, critical-path impact,
  value of the next computation, and one next question. Prefer a different reviewer
  or independent mathematical/implementation route where available.
- Edit the SAME blueprint's current dashboard, card, work-package status, and
  latest review. Append historical decisions to research_log or dated evidence;
  synchronize current checklist fields and NEXT_SESSION. Do not fork master plans.
- Accepted scoped packages remain closed until a relevant model/code change
  invalidates their evidence. Broader physics gets a new dependent package.
  Official GC milestones retain their original acceptance criteria.
- A failed immutable declaration cannot be rescued by filling its other grids.
  Further runs require a specific diagnostic/reuse decision. No automatic higher
  power/extra seeds or retrospective gate changes after a failed pilot.
- Two consecutive turns without new discriminating evidence, a closed dependency,
  or a counterexample trigger strategy review before a third similar expansion.
  Token use and file/figure/test counts alone are not progress measures.
- Repeat successful full audits/tests only for relevant changes, new failures,
  unresolved concerns, or release requirements. Preserve necessary native checks
  for new paths and all exact-source/reuse contracts.
