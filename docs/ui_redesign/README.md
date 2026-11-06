# GABES UI 리디자인 — 작업 허브

프론트엔드(`streamlit_app.py` 계열) 디자인 리팩터링의 단일 출처. 계산 코드와 독립 진행.
새 세션은 이 파일 → `plan.md` 현재 단계 → `log.md` 최신 항목 순서로 읽고 시작.

SABES에도 같은 시각 규약 적용. 광학 테이블의 glyph·광로·라벨·모델 표시·화면 검수 기준은 [SABES 디자인 규칙](../SABES/design_rules.md), 사용자 결정은 D13 참조.

## 현재 상태
- ⟦AI⟧ 2026-11-05 세 뷰 개선안의 **새 1단계 공용 지도·선택 체계 완료**: Map/Curves, 조건·단면·입력·판독 연동, 네이티브 행 스캔과 캐시·중단/재개. [단계 계획](plan.md), [검증 기록](audit/common_map_stage1_2026_11_05.md). 아래의 기존 G05–G08 완료 이력과 구별한다.
- ⟦AI⟧ 같은 계획의 **2단계 Λ 뷰 완료**: EIT/AT 지도·T/OD·분산, AT 네이티브 피크, CPT 대비/폭/대비÷폭, Δp×Δc 이조 지도와 추가 비교 축. [검증 기록](audit/lambda_stage2_2026_11_05.md).
- ⟦AI⟧ **3단계 리드버그 뷰 완료**: EIT 출력 지도·잔류 폭 비교, AT RF 지도·피크/장 환산, 출력·온도 분석, Spectra/Sensitivity와 네이티브 유한 IF 크기·위상·PSN/총 ASD 연결. [검증 기록](audit/rydberg_stage3_2026_11_05.md). 새 4–5단계는 대기한다.
- 다변수 그래프 3단계 **G05–G08 적용**: Λ·Rydberg·Magneto·PSR의 조건 비교와 공용 탐색, 가이드 연결 완료. 실행 범위·검증은 최신 개발 로그 참조.
- 공통 Atom 팝오버를 6개 스킴에 적용. Magneto는 자동 Signal과 편광각·D1/D2 비교로 개편. 입력·표시·모델 범위는 [D18](decisions.md#d18) 참조.
- 조작 방향: 메인 화면의 **Conditions** 팝오버와 단위가 붙은 직접 입력으로 옮기고, 사이드바와 변수 슬라이더를 점진적으로 줄인다. OD/SAS 기본 화면과 이어지는 스킴 작업은 [D17](decisions.md#d17)을 따른다. 기존 셸·타이포그래피의 설계 이력은 D1–D3에 있다.
- 이름 붙인 설정: **DEFAULTS** 드롭다운·브라우저/로컬 저장·가져오기 반영 (`decisions.md` D14, 2026-10-22).
- 진행 단계: P0 · P1(`4cb011f`) · P2 셸(`7a95e9e`) · P3 ScrubField 완료. 다음: **P4 readout 계약**.
- 푸시: P0–P3 네 커밋 `origin/main`에 반영 (2026-09-23, 사용자 승인). 이후에도 단계마다 승인 받고 커밋·푸시.
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
- 메인 plot은 `stMain`의 `img`, `iframe`, `canvas`, OD/SAS overview의 `.mc-overview svg` 중 가장 큰 그림으로 측정. 인라인 컴포넌트의 shadow root 포함.
- OD/SAS는 메인 화면의 `.st-key-gabes_conditions`에서 조건과 Advanced 팝오버를 제공. 사이드바 왕복 측정은 `N/A`.
- 사이드바가 있는 스킴은 데스크톱·모바일에서 실제 열기·닫기를 검사. Advanced 토글은 sidebar 안 "Advanced controls", "Show advanced" 또는 "Hide advanced"를 포함하는 summary/button.
- 숫자 knob(ScrubField)은 **shadow DOM 안**(`.stBidiComponent` 호스트 → `.sf`). 일반 `document.querySelector`로 안 잡힘 —
  DOM을 훑는 도구는 `el.shadowRoot`를 재귀로 들어가야 함. 값 = `[role=spinbutton]`.

## 유용 스킬
`design:design-critique`, `design:accessibility-review`, `design:design-system`, `design:ux-copy`,
`design:design-handoff`, `dataviz`, `run`, `simplify`, `code-review`.
