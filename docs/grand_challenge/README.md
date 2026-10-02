# Grand Challenge — 시작 안내

**관리 감사: 2026-10-02. 과학 근거: 2026-09-28. 공식 milestone0/4 유지.**
실제 진전 있음. 검증된 범위의 소과제6개 종결 인정, v2 통과 가능성은 실패로 종결.
이7개는 전체 과제의 완료율 아님. [총괄 설계·현재 카드](blueprint.md)를 같은 파일에서 계속 갱신.

| 목적 | 먼저 읽을 문서 |
|---|---|
| 연구 방향·종결 단위·의존성·턴 종료 총책임자 검토 | [살아 있는 총괄 설계](blueprint.md) |
| 다음 한 턴 실행 | [NEXT_SESSION](NEXT_SESSION.md) — 짧은 영어 인계 |
| 이번 감사 판단·재현·검증 | [2026-10-02 감사](audit_2026_10_02/README.md), [후속 코드 리뷰](code_review_2026_10_02/README.md) |
| 특정 이론·실행법·과거 근거 | [DOCUMENT_MAP](DOCUMENT_MAP.md) |
| 공식 네 완료 관문·상태 | [체크리스트](../checklist.json) |
| 과학 문서 개정 | [발행·보존 원칙](publication_policy.md), [현행 판본](current_publications.json) |

Frozen v2는 **96고유 경로·480계산 기록·4/9격자 감사**. P2→P3/seed11은5/6metric이5% 초과,
response36.1694%. 필수 첫 edge 실패 때문에 기존 `[2,3,4]`는 남은5격자를 채워도 PASS 불가.
[봉인된 과학 진행 보고](progress_2026_09_28/README.md)와 실패 원본 보존.

다음 주 공정: **p3/seed811 독립 scramble·tail pilot, 새120계산 상한**.
결과 뒤 새 적분 window 설계 또는 quadrature/경계 전략 전환 판단. 자동 p4/p5 확대 금지.
이번 감사는 새 원자 계산 없이 geometry만 확인. Thermal·optical·실험 squeezing 미인증.

매 연구 턴 한 주 질문 전진 → 검증 → 총책임자 검토 → 같은 blueprint 현황·카드 갱신.
독립 입력 ledger·작은 nonlocal 광장 설계는 병행 가능. 원본 유도·실패·cache·source ZIP은 보존.
날짜가 지난 문서의 “최신/다음”은 작성 당시 기준; 현행 판단과 분리.
