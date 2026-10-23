# SABES 상단 요약 리뷰 · 2026-10-23

직접 사용자 요청: SABES 한정 코드 리뷰, 이슈 즉시 수정, 추가 검수 1회 후 이번 상단 변경 커밋·푸시.

## ⟦AI⟧ 리뷰·추가 검수

- 1차: 상단 세 카드·상태 배지·Detection/Power budget 이동, NF unavailable 처리·Gain 순서·HTML escaping·계산 cache 영향 확인. 수정이 필요한 결함 없음.
- 2차: 같은 변경 재검토, 실제 Chromium 1440/720/390/320px, 상세 탭 값 보존·Lab↔Sim 왕복·제목/배지/사이드바 버튼 확인. 추가 결함 없음.
- `python -m pytest -q tests/test_sabes_page.py tests/test_sabes_current_nf.py tests/test_docs_consistency.py`: **31 통과**, 3.78s. 현재 작업 트리의 페이지·NF 연동 검사.
- 커밋 대상 페이지를 분리해 같은 검사 실행: **30 통과·1 deselect**, 1.46s. 미커밋 canvas 이벤트 검사만 제외. [분리 소스 검증](candidate_validation.json).
- [화면 수치](metrics.json): 세 카드, 카드 overflow·문서 가로 overflow 0. 프로필 변경 후 계산 완료를 기다려 NF 상세 값 일치 확인. 브라우저 console 오류 0.
- PNG: `output/playwright/sabes_readout_review_20261023/`. 브라우저 검수 범위는 상단·상세 탭이며, 화면에 함께 찍힌 미커밋 테이블 배치는 이번 커밋에서 제외.
- 전체 수치 회귀·GC 감사 생략: 계산식·공유 solver·입력·검증 계약 무수정. 10/22 전체 회귀는 [이전 audit](../sabes_readout_20261022/summary.md)의 당시 결과이며, 이번 검증은 위 focused 검사.
- 커밋 전 형식 검사에서 10/22 검증 로그의 중복 CR 줄바꿈 제거. 텍스트·통과 수치 보존, 분리된 diff의 whitespace 검사 통과.

커밋은 `sabes_page.py`의 `render`·`_inject_css` 변경, 디자인 규칙의 상단 요약 절, 해당 로그·audit만 포함. 광학 테이블·iframe 이벤트 수정 및 다른 작업의 staged 항목은 보존.
