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

## P1 · 토큰 · 타이포 · 포맷 (계산 코드 무접촉)
- [ ] IBM Plex Sans(400/500/600/700) + Plex Mono(400/500) woff2 self-host → `static/fonts/` + `OFL.txt`.
      출처·파일명·용량 명시 후 사용자 승인 받고 다운로드
- [ ] `.streamlit/config.toml`: `font`/`codeFont`/`headingFont` + `[[theme.fontFaces]]`(`app/static/fonts/…`),
      `baseRadius`, `buttonRadius`, `borderColor`, `showWidgetBorder`, `showSidebarBorder`, `primaryColor` 토큰화.
      1.54 옵션 확인 완료(`theme.light.*`/`theme.dark.*` 분리 지원 → P5 다크에 재사용)
- [ ] `_inject_css` → `assets/ui/gabes.css` 한 파일 + CSS 변수 토큰 (`decisions.md` D2 색표). testid 선택자 목록 주석화
- [ ] 장식 제거: 그라데이션 헤어라인, hero left-border+gradient, group 색 점(`GROUP_STYLES`/`METRIC_STYLES` 색 → 중립),
      BETA 로즈 pill → 외곽선 태그, Guide 그라데이션 버튼 → ghost
- [ ] 숫자 포맷: `format` 미지정 slider에 step 기반 자릿수 + 실제 값 유효숫자 보존(86.94 유지). 순수 함수 + 단위 테스트
- [ ] matplotlib 폰트: `gabes/plot_style.py` FONT_STACK 앞에 Plex 추가 검토 (없으면 fallback; 축 문자열 ASCII 유지)
- [ ] 완료 기준: pytest 통과 · audit `p1` · 하드코딩 hex 색 `streamlit_app.py`에서 토큰 참조로 대체 · 스크린샷 비교

## P2 · 셸 레이아웃
- [ ] 상단 바: mark+wordmark · scheme 전환기(`role=combobox`, aria-label "Scheme" 유지) · ⓘ About popover
      (caption + `info()` + cluster) · regime/mode segmented (`applies_defaults` param 또는 `recommended_defaults`) ·
      Guide(components v2 런처, iframe 제거) · SABES 링크
- [ ] rail: 로고·Guide·Scheme·Cluster·SABES 제거, 그룹 헤더 중립 small caps, Default 버튼 제거(상단 segmented로 통합)
- [ ] Advanced: 각 그룹 안 인라인 공개, advanced 전용 그룹은 하단, 항목 > 6 그룹은 접힘 + 개수
- [ ] readout strip: hero/ribbon → 한 줄 strip. status → 칩, `delta` → 보조줄, 전체 metric popover(숨긴 metric 있을 때만)
- [ ] plot 카드: carousel iframe → 탭 + `st.image`/`st.pyplot`, figure 제목 앱에선 숨김(export 유지),
      툴바: Overlay data(comparison 있을 때만, 우측 패널) · Export popover
- [ ] More pills: tables · diagnostic figures · extra views(Run) — `segmented_control` 선택 없음 기본
- [ ] 모듈 분할: `streamlit_app.py`(2000줄) → `gabes_ui/{theme,layout,controls,readout,plot,overlay,export,guide}.py`,
      `streamlit_app.py`는 라우팅·캐시만. SABES 라우터 유지
- [ ] 완료 기준: plot 5/5 첫 화면 · iframe 0 · main expander ≤ 1 · pytest · audit `p2` · `design:design-critique`

## P3 · ScrubField
- [ ] `st.components.v2` 인라인 컴포넌트 (html/css/js). props: label, unit, min, max, step, value, default, endpoints,
      format, live(`recompute=False`), help
- [ ] 동작: `decisions.md` D3 전부. 단위 표: power(W…nW), length(m…µm), temperature(°C/K), frequency(Hz…GHz),
      field(T/mT/µT/G), time(s…ns)
- [ ] 상태 동기: `session_state[key]`가 진실 원천. preset/default 적용 → data로 전달 → 컴포넌트 반영. 제스처당 rerun 1
- [ ] fallback: v2 없음/예외 → `st.slider`
- [ ] 단위 파서 Python 미러 + 표 기반 테스트, AppTest로 fallback 경로
- [ ] 접근성: `role=spinbutton`, aria-value*, focus-visible, 라벨 hover help (`aria-describedby`) → `design:accessibility-review`
- [ ] 완료 기준: rail px/control ≤ 60 · help glyph 0 · rail ≤ 900 · ready ms 악화 없음 · audit `p3`

## P4 · readout 계약 (`gabes/schemes/*` 수정 — 계산 세션 커밋 후)
- [ ] `base.py` 계약 문서화: `tier`, `attach_to`, `evidence` (`decisions.md` D5)
- [ ] 5 scheme 태깅. SAS: FWHM + resolution 칩(근거 3) · FWM: 3칸 + 보조줄 · Rydberg/magneto/lambda 동일 규칙
- [ ] magneto regime 이모지 라벨 정리, 긴 choice label 단축 (`design:ux-copy`)
- [ ] figure 제목 속 파생값(Rydberg Ω_c 등) → `detail` metric 승격
- [ ] 완료 기준: strip ≤ 5칸 · 모든 metric ≤ 1 조작 도달 · pytest · audit `p4`

## P5 · 마감
- [ ] 폰 390: plot 상단 ≤ 360, Controls 버튼 → sidebar
- [ ] 다크: `theme.dark.*` + CSS 토큰 + matplotlib 팔레트(`dataviz`) — 선택
- [ ] a11y 감사(WCAG AA), ux-copy 점검
- [ ] SABES 페이지 토큰 정렬
- [ ] User Guide 스크린샷 9장 재캡처 → `docs/Userguide/build_static_guide.py`
- [ ] 최종 audit `final` + `--compare baseline final` 표 → `log.md`

## 추가 후보 (미확정, 착수 전 사용자 확인)
- FWM plot의 δ 마커 직접 드래그 (navigate-only knob, 픽셀↔데이터 매핑 필요)
- scheme deep link (`?scheme=`) — audit·공유 편의
