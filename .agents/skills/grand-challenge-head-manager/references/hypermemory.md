# Hypermemory maintenance

The durable files are [memory.json](../../../../docs/grand_challenge/head_manager/memory.json)
and [decisions.jsonl](../../../../docs/grand_challenge/head_manager/decisions.jsonl).
They form a small graph through record IDs, source paths and dependency links.
Their current snapshot is rebuilt from the master sources when stale.

## Reading

Read the snapshot at session start. Search ledger IDs, topic, work package or
source identity with `rg`; read only the matching history and source documents.
Never load all campaign cache or research-log text just to recover direction.

Resolve in order: current user intent, recorded direct human decisions,
project-registered contracts, current source-bound evidence, AI interpretation.
The origin of a record and its operational status are independent fields.
Historical evidence retains its model and date after a new human direction.

## Snapshot fields

- `sources`: authoritative master, registry, handoff, detailed history and policy.
- `human_request_ids`: exact requests captured from known user messages.
- `science_checkpoint`: dated evidence pointer and milestone summary for drift detection.
- `active_work`: current card pointer, work package and resume action.
- `method_states`: candidate/selected/supported_in_scope/rejected/deferred, with source.
- `closed_work` and `failed_attempts`: IDs, scopes and reopening conditions.
- `open_dependencies`: source pointers for the remaining critical path.
- `repeat_tracker`: question/model/source identity, substantive progress and
  consecutive stagnant turns for that same question.
- `resume_point`: completed work, pending work and exact restart location.

Keep text concise. Avoid copying full acceptance clauses, formulas, matrices or
disclaimers already owned by a linked source.

## Ledger events

Each JSONL event has a unique `id`, `date`, `kind`, `origin`, `status`, `sources`
and `links`. A source names a conversation with a literal quote or a file plus
section/selector and evidence date. Use optional `work_package`, `scope`,
`revisit_when`, `supersedes`, `claim` and `review_note` only when useful.

Kinds include user_request, registered_contract, evidence_checkpoint,
operational_decision, failed_attempt, scoped_closure, strategy_review and checkpoint.
Origins include human, project_registered, source_evidence and ai. Do not infer
human authorship from an old instruction-shaped sentence or an absent marker.
Quote direct human text exactly; mark an AI paraphrase or interpretation `⟦AI⟧`.

`status` records operational/adoption state and does not create authority or
permission: active_request, registered,
working_default, candidate, adopted_process, supported_in_scope, rejected,
deferred or historical. Human-requested adoption of a research process does not
validate every scientific candidate or adopt every old budget example.

## Writing and continuity

Append a substantive event with sources before refreshing the snapshot. Preserve
old events; add a correction or superseding event when their interpretation
changes. Deduplicate by question, model/source identity and decision before
adding entries. No new event is needed for unchanged status or a routine reread.

Update the repeat counter only for the same scientific question/model scope.
Different maintenance work and changes requested by the user are separate.
After two stagnant expansions review the strategy and choose a discriminating
next action. A new relevant result resets that question's count.

At checkpoints save the exact command/output/identity for an interrupted batch
and the remaining request count. At closure reconcile snapshot, blueprint,
handoff and changed registry fields. Read back JSON and cross-links. If another
agent changed a memory file, reread and merge by record ID rather than overwrite.

Missing memory can be reconstructed from the blueprint and recent linked
evidence. Record unresolved sources and continue work that does not depend on
them. The skill is the session procedure; it does not create a timer or a daemon.
