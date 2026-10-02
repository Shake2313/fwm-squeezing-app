# 2026-09-18 과학 문서 반영 기록

사용자 지시: squeezing report는 실험 정합성·frontier, analytic reconstruction은 이론 정밀도, quotient structure는 순서론적 논리 구조를 소유. 각 최신 문서는 이전 판 없이 완결. 실질 내용 보존; 개발·버전 기록은 Grand Challenge에만 둠. 완료 이후에도 같은 방향 유지.

## 산출물과 지속 적용

현행 경로는 [current_publications.json](../current_publications.json), 지속 정책은 [publication_policy.md](../publication_policy.md). 정책을 `CLAUDE.md`, `docs/checklist.json`, Grand Challenge blueprint에서 연결. 루트와 문서별 README는 현행 문서로 연결.

| 산출물 | 반영 범위 |
|---|---|
| squeezing report v8 | 실험값·적합값·조건부 계산 구분; 분포 잡음·유한 seed·공차; surrogate frontier의 조건과 반례; 열적 앙상블의 미수렴과 강한 pump 경계 조건 |
| analytic reconstruction v3 | 기존 유도·증명·수치 사례; GKSL·ordered diffusion·QRT; reciprocal normalization; sideband·periodic·spatial·moving-atom 이론; 경계·Poisson noise; nonlocal 광장과 수렴의 범위 |
| quotient structure v4 | 기존 계산 도식; 고정 비교 계약 아래 정보 준순서와 동치류 부분순서; 가정 함의·수치 세분·증거 관계 구별; 조건부 microscopic/transport 연결 |

물리 엔진이나 기존 원시 계산은 이 문서 갱신을 위해 변경하지 않음. 기존의 광범위한 작업 트리 변경·삭제·staging은 그대로 둠. 이전 판 파일을 옮기거나 삭제하지 않음. 기존 `fwm_quotient_site`는 보존한 옛 출판물이며 이번의 현행 standalone HTML과 구별.

## 보존과 교정

- [historical_artifacts_before.json](historical_artifacts_before.json): 기존 문서·생성기·site 31개 파일의 시작 SHA-256.
- [squeezing_coverage.md](squeezing_coverage.md): v1–v7·paper의 주제별 매핑, 포함·통합·교정 이유. v7의 23 equation·3 align 블록, 기존 label과 그림 5개 보존. 오래된 대리 목적함수의 값은 모델과 실패 조건을 붙여 보존.
- [analytic_coverage.md](analytic_coverage.md): v2의 133개 label 보존, 최종 188개. v1에만 있던 반례/수치도 독립 본문으로 복원. 일반 Lyapunov spectrum의 네 고유값 합, manifold별 정규화, 강구동 population·noise 구분 교정.
- [quotient_coverage.md](quotient_coverage.md): 기존 61개 node ID와 98개 source–target 관계 보존. 최종 12개 group·87개 node·142개 relation. 개발 일정/완료율 대신 실제 모델·표본·수렴 조건만 과학 본문에 둠.
- 증거·한계의 범위를 세 문서에서 교차 검토. 미시적 계산이 존재한다는 사실과 실제 열원자 장치의 절대 예측이 성립한다는 주장을 구별.

## 확인 기록

`tests/test_fwm_docs_consistency.py`의 이론 소스 경로를 고정 v2에서 현행 manifest로 연결. 원래 검사 내용은 유지. 새 v3 대상으로 문서 일관성 검사 7개 통과.

Quotient는 로컬 서버와 Codex Browser에서 실제 렌더·동작 확인: 전체 지도, L/T 그룹 선택, quotient/thermal convergence 상세, 검증 보기, 상세→전체 문서 앵커 이동, 수식과 조건 표시. 브라우저 오류 로그 없음. 내부 링크와 독립 HTML의 네트워크 자산 의존성은 별도 구조 검사로 확인. 인쇄 대화상자와 다운로드는 실행하지 않음.

이론 본문에서 추가 사용한 문헌의 서지·주제는 원저자 기록으로 확인: [Gaussian Quantum Information](https://arxiv.org/abs/1110.3234), [Numerical time propagation of quantum systems in radiation fields](https://arxiv.org/abs/1205.1379). 프로젝트의 수치 결과는 이 문헌을 출처로 돌리지 않고 해당 계산 artifact로 연결.

최종 PDF는 squeezing 23쪽, analytic 53쪽. 전 페이지 렌더·시각 검토 완료. 참조 누락·빠진 glyph·페이지 밖 본문 없음. 이론 PDF의 분수 렌더 오류를 발견해 교정하고 최종 전체 페이지를 다시 렌더. Squeezing 문헌 비교의 계산 설정도 원시 JSON과 대조해 pole response/Ultra phase 설정을 명시.

최종 구조 검사 **19/19 통과**: 이전 31개 파일의 SHA-256 동일, analytic label 133개·report 그림 5개·quotient node 61개/관계 98개 보존, 내부 앵커 및 로컬 근거/동반 파일 26개 연결 정상. 구조 검사는 의미 보존이나 물리 검증을 대신하지 않는다.

전체 `python -m pytest -q`: **1942 passed, 3 skipped, 1 failed**, 822.78초. 유일한 실패는 이전 점검과 동일한 `FWM_physics.tex` 부재. 현행 이론 문서 검사 7개 포함해 새 실패 없음. 현재 문서 갱신이 열적 수렴·독립 실험 검증 또는 Grand Challenge milestone을 완료한 것은 아님.

- 재실행 가능한 구조 점검: [verify_publications.py](verify_publications.py)
- 최종 구조 검사 및 산출물 해시: [validation.json](validation.json)
- 브라우저 실측 확인: [browser_qa.json](browser_qa.json)
- 이론 PDF 페이지/경계 확인: [analytic_pdf_validation.json](analytic_pdf_validation.json)
- Squeezing 전체 보존·빌드·시각 확인: [squeezing_build/validation.json](squeezing_build/validation.json)
- 전체 테스트 출력: [pytest_full.log](pytest_full.log)
- 이전 점검에서 이미 있던 실패: `test_repository_visibility_wording_is_consistently_public`가 작업 트리에서 삭제된 `FWM_physics.tex`를 읽음. 실패를 숨기려고 파일을 복원하거나 과학적 assertion을 제거하지 않음.
