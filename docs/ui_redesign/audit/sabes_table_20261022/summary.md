# SABES 광학 테이블 마감 · 2026-10-22

직접 사용자 요청: 설정한 디자인 규약을 실제 Part1/Part2에 적용. 라벨 겹침·잘림 해결, 옵틱 기호·광로 배치 정리.

## ⟦AI⟧ 반영·검증

| 항목 | 결과 |
|---|---|
| 배치 | Part1 1200×680, Part2 1200×640의 schematic 좌표. 수평·수직 광로 정렬, 입력 PBS 되짚기 제거, CL을 MOPA 입력에 정렬 |
| 라벨 | Plex Sans 12px, 불투명 글자색, 회전 없는 별도 레이어·개별 anchor. 두 Part 모두 경계 잘림·라벨 간·glyph/beam 겹침0건 |
| 기호 | 슬림한 lens ellipse, waveplate hatch, PBS diagonal, cell windows, iris aperture, BPD 두 입력. 90° steering/SQL mirror의 도식 반사 방향 일치. 기존 console 토큰·광선 식별색 유지 |
| 검출기 | 서로 다른 두 표시 포트에 probe/conjugate 연결. 렌더링과 빔 측정이 같은 `Layout.points` 사용 |
| 보기 | Fit·100%·확대/축소, 스크롤·마우스/터치 이동, Part별 보기 상태 보존 |
| 조작 | Enter/Space 선택, 반복 클릭, Part 재진입 첫 클릭, 확대 후 빔 측정 통과. 드래그가 측정점 클릭으로 처리되지 않음 |
| 물리 | 68 node·69 ordered link 의미 필드, 28 source stage 출력, 경로 손실·거리·검출 geometry·물리 소스8개 보존 |
| 회귀 | 전체2627개 통과·3 skip(1110.05s). 마지막11개 미러 기울기와 반사 방향 회귀 검사 추가 후 SABES·문서250개 통과(48.30s) |

[측정값](metrics.json), [변경 전 물리 계약](physics_before.json), [물리 불변 검증](physics_after_validation.json), [전체 회귀](full_pytest_receipt.json), [최종 기호 검증](final_presentation_receipt.json).

전체 검사 실행 중 이번 SABES 소스는 불변. 다른 세션이 수정한 FWM·공통 UI 소스4개는 receipt에 별도 기록·보존. 전체 검사 후 이번 소스 변경은 11개 display angle과 반사 방향 검사뿐임을 hash 재구성으로 확인하고 관련250개를 재실행. 최종 화면에서도 라벨 충돌0건 재확인. D-mirror의 도식 분리 각도로 실제 기계 각도를 역산하지 않음.

화면 검증: Chromium1440×900·390×844, 기본 Fit·100%·테이블200%. 브라우저200% 배율의 레이아웃은 720×450 CSS viewport/DPR2로 등가 검증; 네이티브 브라우저 메뉴 조작과 구분. 실제 폰트 로드 후 SVG bounds에 2 display-unit 여유, 광선 stroke 반폭을 포함해 교차 판정. 화면을 직접 확인한 결과와 함께 기록.

로컬 이미지·snapshot: `output/playwright/sabes_table_20261022/`의 `part1_1440.png`, `part2_1440.png`, `part1_390_fit.png`, `part2_390_fit.png`, `part1_390_detail.png`, `part2_390_detail.png`, `part2_zoom200_1440.png`, `part*_browser200_equivalent.png`. PNG는 기존 저장소 정책에 따라 비추적.

최초 기본 `python -m pytest -q`는 다른 작업의 `output/commit_fwm_experiment_20261022/release/tests`까지 수집해 모듈 충돌120건으로 중단. 생성물 폴더만 제외한 `python -m pytest -q --ignore=output`으로 정식 저장소 검사를 실행. 최초 기록은 `collection_collision_*`에 로컬 보존.
