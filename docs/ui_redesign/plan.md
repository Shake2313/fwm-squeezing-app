# 단계 계획

상태 표기: `[x]` 완료 · `[ ]` 대기 · `[~]` 진행 중 · `[-]` 보류(사유 명시).
단계 종료 = 체크리스트 완료 + 완료 기준 충족 + `audit/<단계>` 측정 + `log.md` 기록.

## 목표 수치 (baseline → 목표)
baseline = `audit/baseline/summary.md` (commit `2998dc1`, 1440×900 / 390×844). 측정 도구 `tools/ui_audit.py`.

| 지표 | baseline (sas · lambda · rydberg · magneto · fwm) | 목표 | 담당 단계 |
|---|---|---|---|
| plot 첫 화면 완전 노출 (desktop) | F · T · T · T · F | 5/5 T | P2 |
| plot 상단 y, desktop | 478 · 385 · 423 · 406 · 377 | ≤ 260 | P2 |
| plot 상단 y, 폰 | 732 · 511 · 586 · 638 · 494 | ≤ 360 | P2·P5 |
| rail 높이 (Advanced 닫힘) | 1402 · 1558 · 1610 · 1666 · 1922 | ≤ 900 (1화면) | P2·P3 |
| rail px / main control | 234 · 173 · 161 · 185 · 148 | ≤ 60 | P3 |
| rail 높이, Advanced 펼침 | 2236 · 2017 · 3657 · 2672 · 2261 | ≤ 1300 | P2·P3 |
| help `?` glyph | 35 · 22 · 62 · 42 · 34 | 0 (라벨 hover) | P3 |
| iframe | 2 · 1 · 1 · 1 · 1 | 0 | P2 |
| main expander | 4 · 2 · 6 · 4 · 5 | 0–1 (More pills 대체) | P2 |
| 초기 전송량 | 1680 KB / 67 요청 | ≤ 1680 + 폰트 150 KB | 전 단계 |
| switch 후 ready (ms) | 5267 · 3146 · 1725 · 1726 · 1728 | 악화 금지 | 전 단계 |

`text leaves, 1st screen`은 기록만 (control이 압축되면 오히려 늘 수 있어 목표 아님).
ready ms 잡음: 같은 commit 두 번 측정에서 lambda 2424 → 3146 (+30 %). 단일 측정 비교 금지 — 악화 의심 시 3회 중앙값.
SAS ready = 첫 페이지 로드 포함 (다른 scheme은 전환 후).

## P0 · 준비 — 2026-09-21
- [x] 현재 화면 진단 + 후보 평가 + 시안 (canvas) — `decisions.md` D1–D4
- [x] 작업 폴더 `docs/ui_redesign/` (README · plan · decisions · log)
- [x] 측정 도구 `tools/ui_audit.py` (Edge/Chrome headless, 추가 패키지 0)
- [x] baseline 측정 + 스크린샷 10장 `audit/baseline/`
- [x] CLAUDE.md "Where things live"에 포인터 1줄
- [x] git index 정리: stale staged 40개 unstage(`git restore --staged`, worktree 무변경) + 5일 된 `.git/index.lock`·`AUTO_MERGE.lock` 삭제
- [x] 커밋: P0 파일만 (`docs/ui_redesign/**`, `tools/ui_audit.py`, CLAUDE.md 해당 3줄) — 임시 index + commit-tree

## P1 · 토큰 · 타이포 · 포맷 (계산 코드 무접촉) — 2026-09-21
- [x] IBM Plex woff2 12종 self-host → `static/fonts/` + `OFL.txt` (jsDelivr `@fontsource/*@5`, 207 KB, 사용자 승인)
      Sans latin 400–700 · latin-ext/greek 400–600 · Mono latin 400–500
- [x] `.streamlit/config.toml`: Plex `font`/`headingFont`/`codeFont` + `[[theme.fontFaces]]` 12개(unicodeRange),
      `borderColor`, `showWidgetBorder`, `showSidebarBorder`, `baseRadius`/`buttonRadius` 0.5rem, `linkColor`, sidebar 배경
- [x] 토큰 단일 출처 `gabes_ui/theme.py` (LIGHT + iframe용 DARK 최소) → `--g-*` CSS 변수.
      `_inject_css` 458줄 → `assets/ui/gabes.css` (Streamlit testid 목록 머리 주석)
- [x] 장식 제거: 그라데이션 헤어라인 · hero left-border+gradient · 그룹 색 점(`GROUP_STYLES`/`METRIC_STYLES`/
      `_concept_style`/`_metric_style` 삭제) · BETA 로즈 pill → 외곽선 · Guide 그라데이션 버튼 → ghost + SVG 아이콘(이모지 제거)
- [x] 숫자 포맷 `gabes_ui/format.py::slider_format` (`decisions.md` D7) + `tests/test_ui_format.py` 23개.
      slider 표시 `40.00 → 40`, `75.00 → 75.0`, `0.00 → 0.0`; hero 숫자 Plex Mono, ribbon Sans tabular
- [-] matplotlib 폰트: 보류 — woff2는 matplotlib 불가, TTF 별도 다운로드 필요. P2에서 앱 plot 제목 숨기면
      남는 건 축 라벨·눈금뿐 → P5 재검토
- [x] 완료 기준: `streamlit_app.py` 하드코딩 hex 0 · audit `p1` · 스크린샷 확인 · SABES 페이지 오류 0 · pytest (`log.md`)

## P2 · 셸 레이아웃 — 2026-09-22
- [x] 상단 바: 브랜드(접힘 때만) · scheme 전환기(`role=combobox`, aria-label "Scheme") · About popover(cluster + caption + `info()`) ·
      regime/preset segmented (`decisions.md` D8) · Guide(`st.components.v2`, iframe 제거) · SABES
- [x] rail: 로고·Guide·Scheme·Cluster·SABES·Default 버튼 제거, 52 px 브랜드 띠
- [x] Advanced: 그룹 안 인라인 공개, advanced 전용 그룹 하단, 7개 이상 접힘+개수
- [x] readout strip: status 칩 · delta 보조줄 · "+N more" popover(전 metric + help)
- [x] plot 카드: 밑줄 탭(선택 그림만 렌더) · 제목 → 캡션 · Overlay data 토글 → 우측 패널 · Export popover(지연 생성)
- [x] More pills: 표 · 진단 그림 · extra view(Run)
- [~] 모듈 분할: `gabes_ui/{export,guide,readout}.py` 분리 → `streamlit_app.py` 2007 → ~1070줄.
      상단 바·rail·plot 카드 조립은 아직 `streamlit_app.py` (테스트가 `_render_param` 등을 AST로 떼어 씀) → P3에서 controls 모듈과 함께
- [x] 테스트 계약 갱신: `tests/test_sabes_page.py` 라우터 순서 검사 → `st.sidebar.` / `key="gabes_topbar"` 기준
- [x] `tests/test_ui_readout.py` 10개 (partition·strip·칩·캡션·이모지 제거)
- [x] `tools/ui_audit.py`: 버튼형 Advanced 토글 대기, 전환 재시도, 폰은 sidebar 닫고 전환, 실패 시 `_failure.png`
- [x] 완료 기준: plot 5/5 첫 화면 · iframe 0 · main expander 0 · audit `p2` · `design:design-critique`(`log.md`) · pytest

## P3 · ScrubField
- [ ] `st.components.v2` 인라인 컴포넌트 (html/css/js). props: label, unit, min, max, step, value, default, endpoints,
      format, live(`recompute=False`), help
- [ ] 동작: `decisions.md` D3 전부. 단위 표: power(W…nW), length(m…µm), temperature(°C/K), frequency(Hz…GHz),
      field(T/mT/µT/G), time(s…ns)
- [ ] 상태 동기: `session_state[key]`가 진실 원천. preset/default 적용 → data로 전달 → 컴포넌트 반영. 제스처당 rerun 1
- [ ] fallback: v2 없음/예외 → `st.slider`
- [ ] 단위 파서 Python 미러 + 표 기반 테스트, AppTest로 fallback 경로
- [ ] 접근성: `role=spinbutton`, aria-value*, focus-visible, 라벨 hover help (`aria-describedby`) → `design:accessibility-review`
- [ ] 모듈 분할 마무리: controls(rail 조립·ScrubField)·shell(상단 바)·plotcard → `gabes_ui/`, AST 테스트 대상 함수 이동 시 테스트도 갱신
- [ ] slider 끝점 캡션 행(◀ OD / SAS ▶)을 track 양끝 인라인으로 (P2 critique)
- [ ] 완료 기준: rail px/control ≤ 60 · help glyph 0 · rail ≤ 900 · ready ms 악화 없음 · audit `p3`

## P4 · readout 계약 (`gabes/schemes/*` 수정 — 계산 세션 커밋 후)
- [ ] `base.py` 계약 문서화: `tier`, `attach_to`, `evidence` (`decisions.md` D5)
- [ ] 5 scheme 태깅. SAS: FWHM + resolution 칩(근거 3) · FWM: 3칸 + 보조줄 · Rydberg/magneto/lambda 동일 규칙
- [ ] magneto regime 이모지 라벨 정리, 긴 choice label 단축 (`design:ux-copy`)
- [ ] figure 제목 속 파생값(Rydberg Ω_c 등) → `detail` metric 승격 → 캡션 유지 여부 재결정 (D8)
- [ ] status 칩을 대상 값 옆으로 (`attach_to`) — P2에선 strip 끝 (critique)
- [ ] SAS strip: 근거 진단값(edges·samples·scan-edge) 대신 Doppler FWHM·phase가 보이도록 tier 태깅
- [ ] 완료 기준: strip ≤ 5칸 · 모든 metric ≤ 1 조작 도달 · pytest · audit `p4`

## P5 · 마감
- [ ] 폰 390: plot 상단 ≤ 360 (P2: 349–530), 보조 metric 폰에선 2칸, plot 머리 Export 행 축소, 터치 타깃 44 px
- [ ] 다크: `theme.dark.*` + CSS 토큰 + matplotlib 팔레트(`dataviz`) — 선택
- [ ] a11y 감사(WCAG AA), ux-copy 점검 (Export 단일 곡선 이름 "primary: line N", 긴 extra view 이름)
- [ ] SABES 페이지 토큰 정렬
- [ ] User Guide 스크린샷 9장 재캡처 → `docs/Userguide/build_static_guide.py`
- [ ] 최종 audit `final` + `--compare baseline final` 표 → `log.md`

## 추가 후보 (미확정, 착수 전 사용자 확인)
- FWM plot의 δ 마커 직접 드래그 (navigate-only knob, 픽셀↔데이터 매핑 필요)
- scheme deep link (`?scheme=`) — audit·공유 편의
