# GABES UI 리디자인 — 작업 허브

프론트엔드(`streamlit_app.py` 계열) 디자인 리팩터링의 단일 출처. 계산 코드와 독립 진행.
새 세션은 이 파일 → `plan.md` 현재 단계 → `log.md` 최신 항목 순서로 읽고 시작.

## 현재 상태
- 방향 확정: **시안 B · Instrument console** + IBM Plex Sans/Mono + ScrubField 정밀 입력 (`decisions.md` D1–D3).
- 진행 단계: P0 · P1(`4cb011f`) · P2 셸(`7a95e9e`) · P3 ScrubField 완료. 다음: **P4 readout 계약**.
- 푸시: 보이는 뷰가 어느 정도 완성되면(사용자 판단) — 그 전엔 로컬 커밋만.
- 시안 canvas (비공개): https://claude.ai/artifact/C8d3dLBYDjH3ZJm84FPEQw
  — 현재 화면 진단, B(SAS·FWM·폰), A(보수안), ScrubField 동작 프로토타입.

## 파일 지도
| 파일 | 용도 | 갱신 시점 |
|---|---|---|
| `README.md` | 허브: 상태·규칙·파일 지도 | 단계 전환 시 |
| `plan.md` | 단계별 작업 체크리스트 + 완료 기준 + 목표 수치 | 작업 완료·범위 변경 시 |
| `decisions.md` | 결정 기록(평가표·대안·근거) | 결정 추가·번복 시 (번복도 새 항목) |
| `log.md` | 날짜별 개발 로그, 최신이 위 | 작업 세션마다 |
| `audit/<tag>/` | `tools/ui_audit.py` 측정값(`metrics.json`, `summary.md`) + 첫 화면 스크린샷 | 단계 종료마다 새 tag |

스크린샷 PNG는 repo `.gitignore`(`*.png`) 정책으로 git 제외 — Google Drive 동기 폴더에만 보존. 수치는 git에 남음.

## 작업 규칙
- **범위**: `streamlit_app.py`(라우터·상단 바·흐름), UI 패키지 `gabes_ui/`, `assets/`, `static/`,
  `.streamlit/config.toml`, `sabes_page.py`의 공유 스타일. 계산 경로(`compute`, `observables`의 물리값) 무수정.
  - `theme` 토큰 · `format` 숫자 포맷 · `units` 단위 파서 · `scrub` ScrubField 컴포넌트 ·
    `controls` rail 조립 · `shell` 상단 바 조각 · `plotcard` 그림·오버레이·More · `readout` strip · `export` · `guide`
- **예외**: P4 metric 표시 계약(`tier`/`attach_to`/`evidence`)만 `gabes/schemes/*` 수정. 계산 세션 커밋 후 진행.
- **감량 원칙** (`decisions.md` D4): 새 요소는 기존 요소 ≥1 대체 또는 조건부 표시만 허용. "혹시 쓸모" 추가 금지.
- **matplotlib 축 문자열 ASCII 유지** (CLAUDE.md mathtext layout-lock). UI 라벨·단위는 unicode.
- **Streamlit 버전**: 1.54 기준 설계 (`st.components.v2`, `st.popover`, `segmented_control`). 버전 올리면 `ui_audit` 재측정.
- **git**: main 직접 작업. 작업트리에 다른 세션 staged/미커밋 변경 섞임 → 커밋은 내 파일만
  (임시 index + commit-tree). 다른 세션 staged 파일이 남아 있으므로 일반 `git commit` 금지.
  `MM` 표시는 `git diff HEAD -- <f>`로 실제 변경 여부 확인 (stale index 사례: `log.md` 2026-09-21).

## 검증 절차 (단계마다)
1. `python -m pytest -q` (FWM baseline 포함).
2. `python tools/ui_audit.py --tag <단계>` → `python tools/ui_audit.py --compare baseline <단계>`.
3. 비교표를 `log.md`에 붙이고 `plan.md` 목표 수치 대조.
4. 시각 변화 → 스크린샷 확인 (`audit/<tag>/*.png`), 필요 시 `design:design-critique`·`design:accessibility-review`.

## ui_audit DOM 계약 (리디자인 중 유지, 바꾸면 도구도 같은 커밋에서 수정)
- `stApp`의 `data-test-script-state` 속성 (Streamlit 제공).
- scheme 전환기 = `role=combobox`, `aria-label`에 "Scheme" 포함; 옵션 = `role=option`, 텍스트 = scheme title.
- 메인 plot = `stMain` 안 최대 `img`/`iframe`/`canvas`.
- control rail = Streamlit sidebar.
- Advanced 토글 = sidebar 안 summary/button, 텍스트 "Advanced…" 또는 "Show advanced".
- 숫자 knob(ScrubField)은 **shadow DOM 안**(`.stBidiComponent` 호스트 → `.sf`). 일반 `document.querySelector`로 안 잡힘 —
  DOM을 훑는 도구는 `el.shadowRoot`를 재귀로 들어가야 함. 값 = `[role=spinbutton]`.

## 유용 스킬
`design:design-critique`, `design:accessibility-review`, `design:design-system`, `design:ux-copy`,
`design:design-handoff`, `dataviz`, `run`, `simplify`, `code-review`.
