---
name: grand-challenge-head-manager
description: Review GABES Grand Challenge sessions from the project head's long-term perspective, using durable linked memory to keep work aligned with user intent, choose reusable research methods, prevent repeated work, and separate AI commentary from project goals. Use at session start, meaningful direction changes, and closure; not for unrelated UI or engine work.
---

# Grand Challenge Head Manager

Policy source: the user's 2026-10-13 request. Project root is three directories
above this skill. Apply the repository instructions and
[internal provenance policy](../../../docs/grand_challenge/internal_provenance.md).

## Responsibility

Keep the full path in view: independent apparatus inputs -> one atomic
mean/response/noise model -> optical propagation -> detected absolute gain and
S_minus -> experimental validation. Judge a session by the dependency it closes,
counterexample it finds, or decision it resolves. Support the user's requested
work, including deliberate exploratory work and changed priorities.

Current human instructions take precedence over old plans, memory summaries and
AI working defaults. Record a changed goal or scope and adapt the next work;
retain the identity and result of historical experiments. Do not ask for extra
approval solely because an old AI suggestion set a candidate count or budget.

## Durable hypermemory

At start read [memory.json](../../../docs/grand_challenge/head_manager/memory.json),
then the current dashboard/card in [blueprint](../../../docs/grand_challenge/blueprint.md)
and [handoff](../../../docs/grand_challenge/NEXT_SESSION.md). Read only relevant
registry fields, ledger events and linked evidence. Use
[memory format](references/hypermemory.md) when writing or reconstructing memory.

The memory connects requests, decisions, closed scopes, failures, methods and
dependencies through stable IDs and source links. Blueprint remains the master
design; checklist owns registered milestones; sealed results own calculation
facts. Reconcile stale memory against those sources before using its snapshot.
Unknown historical authorship stays unknown; a label such as "official" does
not establish that a human wrote every clause.

## Session lifecycle

**Start.** Identify the latest user request and its relation to the long-term
goal. Match it to the active card or create a scoped card in the same blueprint.
Check closed work and failed attempts before choosing computations. Select a
useful method from [research methods](references/research-methods.md); ordinary
maintenance can proceed without a new research method. Preserve multiple user
deliverables when the request spans more than one question.

**Checkpoint.** Review before an expensive batch, a model/source/gate change,
publication, a substantial pivot, or resuming after a context break. Save a
verified resume point after meaningful partial progress. Check whether the work
still resolves the user's question, whether a different approach would be more
discriminating, and whether repeated text has displaced substantive work.

⟦AI⟧ Long numerical runs: a saved `RUNNING` label is not process evidence.
Check the recorded PID/controller and completion receipt before resuming. Prefer
an independently running, hidden batch with distinct logs and a completion receipt
when a chat interruption would otherwise lose work. Preserve finished cache records
and keep explicit user stop/pause requests authoritative.

**Close.** Review new evidence, scope, dependency impact and the value of the
next computation. Update changed current fields in the same blueprint,
registry and handoff; append a concise event to the durable ledger and detailed
history only when needed. Refresh the memory snapshot from those sources and
give the user the outcome, verification and next action.

Run these reviews in the current agent. Delegate an independent mathematical,
implementation or adversarial route when it adds useful evidence. Separate
contexts before exchanging conclusions for blind re-derivation. A lightweight
review is enough for a small documentation task; do not manufacture extra logs.

## Long-horizon review

- Trace the selected work through open dependencies to the original objective.
  Keep a valuable local success closed and attach broader physics as new work.
- Before repeating a computation or audit, name the changed input/model/code,
  new hypothesis, failed check or release requirement that justifies it.
- Apply the project's two-stagnant-turn strategy review to repeated expansions.
  Continue toward the user's goal through a better question, method or scope.
- Preserve model, numerics, apparatus inputs and measurement definition as
  separate axes. Keep current card costs distinct from global resource policy.
- For a new experiment declare its question, distinguishing result, reuse/new
  solve forecast, evaluation split and stop/pivot rule before seeing results.
- Review internal writing for repeated disclaimers. Keep necessary mathematical
  assumptions and result typing close to the claim; put supplementary AI
  interpretations once in a short `⟦AI⟧` section. Preserve that origin in memory
  and summaries. Public/frontend prose follows its own reader-facing policy.

## Review output

Use a compact internal record:

`⟦AI⟧ <question> -> <new evidence or decision> -> <effect on long-term path> -> <next action>`

Link the source and distinguish an observation, working default and human
decision. Do not turn a reviewer recommendation into a new user restriction.
Store only changes worth recovering in a later session.
