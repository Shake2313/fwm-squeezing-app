# 결정 기록

번복도 삭제 말고 새 항목으로 추가 (`D<n>`, 날짜, 대체 대상 명시).

## D1 · 레이아웃: 시안 B Instrument console — 2026-09-21
- **결정**: 상단 바(브랜드 · scheme 전환 · ⓘ About · regime/mode · Guide · SABES) + 좌측 control rail(sidebar) + plot 우선 메인
  (readout strip → plot 카드(탭·툴바) → More pills).
- **통과 조건**: Streamlit 유지 · 계산 경로 무수정 · 300 KB 넘는 JS 라이브러리 금지.
- **평가** (1–5점, 가중치 합 100 → 100점 환산):

| 후보 | 절제 20 | 조작 루프 20 | 도달 15 | 완성도 15 | 경량 10 | 위험 10 | 계약 5 | 반응형 5 | 점수 |
|---|---|---|---|---|---|---|---|---|---|
| **B Instrument console** | 5 | 5 | 5 | 5 | 4 | 3 | 4 | 4 | **92** |
| A Refined sidebar | 3 | 4 | 5 | 3 | 5 | 5 | 5 | 3 | 80 |
| C 3-pane workbench | 2 | 4 | 5 | 4 | 4 | 3 | 4 | 1 | 70 |
| D 노트북식 흐름 | 4 | 1 | 4 | 4 | 5 | 4 | 4 | 4 | 70 |
| E 전면 plot + floating HUD | 5 | 4 | 3 | 4 | 2 | 1 | 3 | 2 | 68 |
| F SPA 재작성 | 통과 조건 탈락 (Streamlit 이탈, scheme 계약 파기, 수개월) | | | | | | | | — |

- 기준 정의: 절제 = 첫 화면 정보량·중복; 조작 루프 = control과 plot 동시 가시; 도달 = 모든 정보 ≤1 조작;
  위험 = Streamlit 1.54 구현 난이도·testid CSS 의존; 계약 = scheme 계약 변경량.
- **B 대가**: 상단 바 CSS 배치, metric 표시 계약 소폭 추가(D5). **A**는 fallback (위험 최저, 절제 약함).

## D2 · 시각 언어: 정돈된 계측기, 밝은 테마 — 2026-09-21
- 글꼴: **IBM Plex Sans**(UI) + **IBM Plex Mono**(숫자, tabular). self-host woff2 (~150 KB, 캐시). 로고 BETA 배지가 이미 Plex 지정.
- 색: 중립 쿨 그레이 바탕 + 강조 1(#0284C7 채움 / #0369A1 텍스트) + 상태 2(주의 #92400E on #FEF3C7, 정상 #047857).
  muted 텍스트 #5B6B80 (흰 바탕 5.3:1, AA).
- 금지: 그라데이션 워시, 좌측 테두리 강조 카드, 그룹별 색 점, 크롬 이모지.
- 대안: 다크 오실로스코프(→ P5 테마 옵션, matplotlib·export 흰 배경과 충돌), 에디토리얼 serif(숫자 판독성 낮음).

## D3 · 정밀 입력: ScrubField — 2026-09-21
- **결정**: 값 = spinbutton. 좌우 드래그 scrub(Shift ×0.1), 클릭/Enter 직접 입력(단위 변환), ↑↓·PgUp/PgDn·Home/End,
  범위 밖 → 경계 clamp + 알림, 기본값과 다르면 점(클릭 = 복원). 3 px track = 거친 조정. 제스처 1회 = rerun 1회;
  navigate-only(`recompute=False`) knob은 드래그 중 125 ms 간격 스트리밍.
- **구현 수단**: `st.components.v2` (1.54 확인: 인라인 HTML/CSS/JS, iframe 없음, 양방향 state, npm 빌드 불필요). 실패 시 `st.slider` fallback.
- **평가** (6기준 동일 가중, 30점): 슬라이더+number_input 병렬 18 · 정밀 모드 토글 21 · 값 클릭 popover 22 ·
  slider 키입력 hack 17 · **ScrubField 27**.
- 프로토타입: canvas `ScrubField 프로토타입` 보드.

## D4 · 정보 감량 규칙 — 2026-09-21
판별 질문: "knob 돌리는 동안 매 순간 봐야 하나?"
- **유지**: control, 핵심 metric, plot.
- **중복 병합**: 제목 = scheme 전환기 · Cluster = 전환 메뉴 구획명 · Default 버튼 = regime/mode segmented · plot 제목(파라미터 반복) = rail 값.
- **판정에 근거 부착**: status metric → 해당 값 옆 칩, 근거 metric(`evidence`)은 칩 hover 카드.
  예: SAS `resolution-limited` 칩 ← Samples/FWHM(≥6 필요), Scan-edge distance(≥1 FWHM), Half-height edges.
- **1 조작 뒤로**: caption·Reference → ⓘ · CSV 비교 → plot 툴바 "Overlay data" · export → 툴바 · 표/진단/무거운 스캔 → More pills ·
  Advanced → 각 그룹 안 인라인 공개 (많으면 서브그룹 접힘+개수).
- **라벨 대신 동작**: navigate-only = 드래그 중 실시간 반영 (배지 없음) · min/max = clamp 알림 · help = 라벨 hover.
- **삭제**: 장식(헤어라인·left-border·색 점·이모지) · "Use arrows" · iframe 고정 높이 공백 · "40.00" 꼬리 0
  (유효숫자는 보존: 86.94 유지).
- **부풀림 방지**: 새 요소는 기존 요소 ≥1 대체 또는 조건부만. 허용 예: 변경점(기본값이면 0개), 전체 metric 버튼(숨긴 metric 있을 때만).
- plot 제목: 앱 화면에선 숨김, PNG/JSON export엔 유지. 파생값 포함 제목(Rydberg Ω_c 등)은 metric으로 승격(P4).

## D5 · metric 표시 계약 추가 (예정, P4) — 2026-09-21
- `tier`: `hero` | `key` | `detail` (기존 `hero`·ribbon 호환: 미지정 = `key`).
- status metric(`kind="status"`): `attach_to` = 대상 metric label, `evidence` = 근거 label 목록.
  미지정 fallback: status 칩을 strip 끝에 표시.
- `delta`: 값 아래 muted 보조줄 (칩 아님).
- 순수 표시 메타데이터. 계산 무관. `gabes/schemes/*` 수정이라 계산 세션 커밋 후 진행.

## D6 · 개발 인프라 — 2026-09-21
- 문서 폴더 `docs/ui_redesign/` 단일 출처 (README 허브 · plan · decisions · log · audit).
  `docs/checklist.json`(과학 작업 레지스트리)엔 넣지 않음 — 구조·소유 세션 상이.
- 측정 `tools/ui_audit.py`: Edge/Chrome headless + DevTools(tornado websocket, 추가 패키지 0).
  5 scheme × (1440×900, 390×844) 지표 + 첫 화면 PNG. 단계마다 tag 추가, `--compare`로 비교.

## D7 · 숫자 표시 규칙 · 폰트 적용 범위 — 2026-09-21
- slider 자릿수 = max(step 자릿수, min(기본값 자릿수, step 자릿수 + 2)). 명시 `format` 우선, 정수 slider는 Streamlit 기본.
  예: step 1 → `40`; step 0.1·기본 86.94 → `86.94`; 계산된 기본 1326.2572434514343·step 1 → `1326.26`.
  이유: 손으로 정한 off-grid 기본값은 그대로 보이되, 계산값의 부동소수 꼬리는 막음.
- Plex Mono = hero 숫자·slider 값·눈금. ribbon 값은 Plex Sans tabular (좁은 칸에서 Mono가 줄바꿈 늘림).
- status 문자열(`resolution-limited` 등)은 Sans (숫자로 시작할 때만 Mono, `looks_numeric`).
- matplotlib 글꼴은 DejaVu Sans 유지 (woff2 불가). 재검토 P5.

## D8 · P2 셸 구현 규칙 — 2026-09-22
- **상단 바 regime**: `applies_defaults`인 첫 segmented param(FWM Mode, Λ·Rydberg·magneto Regime)을 상단 바로 이동.
  없으면(SAS) `recommended_defaults` 세트를 segmented로 — 현재 param이 한 세트와 정확히 같을 때만 그 세트 선택 표시
  (값 하나라도 바꾸면 선택 해제 = "지금 프리셋 상태 아님"을 정직하게 표시). 선택지 앞 이모지는 표시할 때 제거.
- **브랜드**: rail 상단 52 px 띠 + 상단 바 같은 높이·선 → 한 줄로 보임. 상단 바 복사본은 sidebar 접힘일 때만, 폰에선 숨김.
- **plot 제목**: 이미지에서 떼어 카드 머리의 muted 캡션으로 (ASCII→기호 복원: Ω<sub>c</sub>, µW, °C, η, →).
  D4의 "숨김" 대신 캡션 유지 — 스크린샷 출처 정보 + 파생값(Rydberg Ω_c, Λ buffer relax)이 아직 metric에 없음. P4 승격 후 재검토.
- **plot view**: carousel(3장 모두 PNG 인코딩) → 밑줄 탭 segmented, 선택된 그림만 렌더 (rerun당 savefig 3 → 1).
- **Overlay data**: 토글 → 우측 패널(2.6:1 열), plot 가림 없음. 토글 끄면 업로드 위젯이 렌더되지 않아 파일 재업로드 필요 — 수용.
- **Export**: popover. JSON 번들·CSV zip은 클릭 시 생성(`download_button` callable data, `on_click="ignore"`) → 일반 rerun에서 직렬화 0.
- **More**: `st.pills`, 기본 선택 없음. 표·진단 그림·무거운 스캔.
- **Advanced**: 버튼 "Show advanced · N". 켜면 각 그룹 안에 인라인, advanced 전용 그룹은 버튼 아래, 7개 이상이면 접힌 expander + 개수.
- **readout strip**: status → 칩 (문구로 톤: limited/required/unavailable… → 주의색), 나머지 hero(`partition_metrics`) +
  보조 4칸, 초과분 "+N more" popover = 전 metric + delta + help 표. `delta` → 보조줄.
- **CSS 주입**: `st.markdown(<style>…<span class=gabes-style-anchor>)` + `:has(.gabes-style-anchor){display:none}`.
  `st.html` style-only는 1.54에서 event container에 들어가지만 실제로 적용 안 됨(빈 element) — 사용 금지.
- **Streamlit 1.54 배치 사실**: key 준 container는 `stLayoutWrapper`로 감싸짐 → 폭·bleed·sticky는 wrapper에 (`max-width` 해제 필요).
  segmented 버튼 testid `stBaseButton-segmented_control(Active)`, pills `stBaseButton-pills(Active)`.
