# Grand Challenge Head Manager 운영 체계 적용

사용자 요청일: 2026-10-13. 요청 기록: [결정 ledger](../head_manager/decisions.jsonl)의
REQ-2026-10-13-METHODS, PROVENANCE, MANAGER.

## 적용한 운영 절차

- 프로젝트 스킬: [Grand Challenge Head Manager](../../../.agents/skills/grand-challenge-head-manager/SKILL.md).
- 매 GC 세션의 시작·중간 checkpoint·종료에서 사용자 요청, 장기 목적, 종결 작업,
  실패와 현재 의존성을 대조. Root/GC AGENTS와 CLAUDE에 진입점 등록.
- 장기 메모리: [snapshot](../head_manager/memory.json)과 append-only ledger를
  요청·근거·결정·방법·실패·재개 지점의 ID/출처 링크로 연결.
- 방법표: bridge fixture, ensemble variance reduction, paired history,
  blind 재유도/negative control, candidate evaluator, input 식별성, mechanism 축약.
- 내부 기호 `⟦AI⟧`: 사용자 직접 요청·등록된 목적·증거와 AI 검증 해석/제안을 구분.
  [규약](../internal_provenance.md). 공개 본문·가이드·프론트엔드에는 적용하지 않음.

## ⟦AI⟧ 편집과 검토

10-07 제안서의 문서 수준 한계를 한 검증 절로 통합하고 AI 적용안의 출처를 표시.
수식·물리 조건·평가 표본 독립성·잡음 source 계약·수치·출처는 유지.
MERGE10·제목/출처 표현 REWRITE; 관련 없는 난이도 순위 논평 삭제.

[편집 전 원본](ai_research_strategy_2026_10_07.original.md)은 바이트 그대로 보관.
원본은 당시 조사 기록으로 보존하며 현행 규칙·다음 카드의 진입점으로 사용하지 않음.
원본 SHA-256: `2230a3a89c993adf52e7bd56c3375368c17bbb6f2e41f16deb077f96f77ae8c8`.

Skill 형식 검사 PASS. 기존 목표/no-fit/milestone, 기술 blueprint 본문,
연구 일지 원래 prefix와 선택 봉인 증거 보존 확인.
검증 수치·독립 시나리오 결과: [verification.json](verification.json).

독립 검토에서 plan과 사용자 권한을 분리하고 evidence heading을 교정.
면책 문구를 합치는 과정에서 미수행 Lean 빌드를 미래 의무로 만들지 않도록 교정.
두 가상 세션으로 사용자 요청의 exploratory 우선·네 후보 비교와 CV 평가 표본
중복 처리를 확인. 가상 요청은 실제 사용자 결정 ledger에 추가하지 않음.

운영 절차를 도입한 이번 작업의 science solve0. 다음 과학 카드의 범위와
검증 조건은 같은 blueprint가 소유. 새 요청이 오면 그 요청에 맞춰 작업과 메모리 갱신.
