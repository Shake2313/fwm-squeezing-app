# SABES 별표·범례 검증 · 2026-10-22

사용자 요청: 기존 디자인 규약 공유, 광학 테이블 표현 기준 강화, 반복 `lumped`→`*`, 선 굵기 범례 설명 제거.

## ⟦AI⟧ 검증 범위·결과

| 항목 | 결과 |
|---|---|
| 규약 연결 | D13·SABES design_rules·기존 optical_table_ui_plan 연결 |
| 전체 회귀 | `python -m pytest -q`: 2570개 통과·3개 skip, 1025.09s |
| 최종 문구 회귀 | 전체 회귀 후 중복 상세창 캡션만 삭제. 관련 73개 재통과, 6.00s |
| 문서·차이 검사 | 문서 일관성4개 통과, scoped `git diff --check` 통과 |
| 실제 화면 | Chromium, 1440×900·390×844, Part1/Part2 screenshot 직접 확인 |
| 별표 | Noise eater·PBS·D-mirror·SQL 등 상세 모델 미구성 부품에 표시. EOM·cell은 기존 상세 모델 상태 유지 |
| 상세 표시 | Noise eater 선택·popup 제목의 별표가 Markdown 강조로 사라지지 않음. 닫기 정상 |
| 범례 | `* Detailed component model not yet implemented` 표시, 기존 선 굵기 문장 제거 |
| 좁은 화면 | 범례 오른쪽291.875px≤viewport390px. 기존 Inspect hardware 상세 접근 유지 |
| 계산 범위 | SABES 물리·기하 소스9개가 HEAD와 동일. 표시 flag·문구만 변경 |

[화면·경계 기록](presentation_validation.json), [전체 회귀 receipt](full_pytest_receipt.json), [최종 문구 후속 검증](final_copy_validation.json).
전체 회귀 중 소스239개 변경 없음. 그 뒤 상세창의 중복 모델 설명 캡션을 삭제하고 관련 테스트와 실제 상세창을 재확인. 캡션의 투과율 적용 주장이 aggregate 효율 override에서는 성립하지 않을 수 있어 삭제했으며, 범례·tooltip의 별표 설명은 유지. 앞서 캡션 표시 수정으로 중단한 실행 기록은 `superseded_before_caption_fix_*`에 보존.
전체 회귀 종료 이후 다른 작업의 비관련 테스트 변경은 후속 receipt에 구분해 기록·보존. 전체 회귀 결과를 그 이후의 비관련 변경 검증으로 확장하지 않음.

화면 자료: `output/playwright/sabes_design_20261022/`의 `part1_1440.png`, `part2_1440.png`, `part1_390.png`, `part2_390.png`, `legend_390.png`, 최종 `final_dialog_390.png` 및 대응 snapshot. 최종 상세창 제목의 `Noise eater*`, 투과율 입력, 닫기를 확인.

새 디자인 규약의 전체 마감 검수는 별도 범위. 현재 Part2의 `Balanced detector` 라벨 잘림·인접 라벨 겹침6쌍은 측정해 다음 배치 검수 입력으로 보존. 이번 문구 치환으로 전체 재배치·zoom·font 전달까지 완료했다고 판정하지 않음.
