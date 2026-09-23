# 개발 로그

최신이 위. 항목 = 날짜 · 단계 · 한 일 · 측정/검증 · 남은 일·주의.

## 2026-09-23 · P3 ScrubField
**한 일**
- P2 커밋 `7a95e9e`. 사전 검증: 후보 트리를 별도로 export해 전체 pytest → 실패 36개 전부 `tests/quantum`(기존),
  UI·앱 테스트 실패 0.
- **ScrubField** `gabes_ui/scrub.py` (`st.components.v2`, 인라인 HTML/CSS/JS, shadow DOM): 숫자 knob 전부 교체.
  값 = `role=spinbutton`. 드래그(Shift ×0.1) · 클릭/Enter 입력 · ↑↓·PgUp/PgDn·Home/End · 3 px track 거친 조정 ·
  범위 밖 clamp + 알림 한 줄 · 기본값과 다르면 점(클릭 = 복원) · 끝점 캡션은 track 양끝 인라인.
- **단위 파서** `gabes_ui/units.py`: 브라우저는 사람이 친 문자열만 보내고 해석은 Python에서만 (D9).
  `1.25 W`→mW 환산, `300 K`→°C, `1 G`→µT, `1,234.5`, `−7`(진짜 빼기표), `1.5e3` 허용. 다른 물리량이면 거부.
- **help `?` 제거**: ScrubField는 라벨 hover, select·segmented는 직접 그린 라벨(hover) + 위젯 라벨 숨김.
  checkbox만 `?` 유지(라벨이 클릭 영역).
- **모듈 분할**: `gabes_ui/controls.py`(rail 조립·`render_param`) · `shell.py`(브랜드·regime/preset) ·
  `plotcard.py`(그림·CSV 오버레이·More·matplotlib 락). `streamlit_app.py` 1072 → 265줄.
  AST로 함수를 떼어 쓰던 `tests/test_fwm_excess_noise.py`는 모듈 import로 교체,
  `tests/test_experimental_csv.py`의 stale 모듈 복구 계약은 `gabes_ui.plotcard` 기준으로 갱신.

**함정 (재발 방지)**
- bidi component id에 `__` 금지 → knob key를 `scrub-scheme-param`으로 접음. 안 그러면 mount가 통째로 예외.
- `setStateValue("value", {객체})`는 Python까지 오지 않음(콜백 미발생). JSON **문자열**로 보내야 도착.
- `AppTest`는 컴포넌트 매니저를 MagicMock으로 바꿔 mount가 `TypeError` → slider fallback 경로가 실제로 돌아감.
  그래서 AppTest로는 ScrubField 자체를 검증할 수 없음 (브라우저 CDP 스크립트로 확인).
- `stMarkdownContainer`에 `margin-bottom: -1rem`이 기본으로 걸려 있어 sidebar의 16 px element gap을 상쇄함.
  직접 그린 라벨에 음수 margin을 더 주면 위젯이 라벨 위로 올라와 글자를 덮음.
- `scrollHeight`는 뷰포트(900) 밑으로 안 내려감 → rail이 1화면에 들어온 뒤로는 비교 불가. `rail_content_height` 신설.
- `stTooltipIcon`은 help 붙은 **버튼**(About·SABES)도 감쌈 → 물음표 개수와 다름. `help_icons`(버튼 제외) 신설.
- 모듈을 옮기면 **소스를 직접 들여다보는 테스트**가 조용히 깨짐. 후보 트리 전체 pytest가
  `test_streamlit_import_repairs_stale_experimental_csv_module`(streamlit_app 네임스페이스에서 CSV API 바인딩 확인)을 잡아냄
  → `gabes_ui.plotcard` 기준으로 갱신. 덤으로 그 테스트의 subprocess를 `encoding="utf-8", errors="replace"`로 바꿈:
  자식 traceback에 한글 경로가 들어가면 cp949 디코딩이 터져 실패 메시지 자체가 못 만들어졌음.
- bash heredoc으로 Python을 쓰면 `\n`·`\w` 이스케이프가 한 겹 벗겨짐 → 이스케이프 많은 파일은 Write/Edit 도구로.

**측정** `--compare baseline p3` (1440×900 / 390×844):

| 지표 | sas | lambda | rydberg | magneto | fwm |
|---|---|---|---|---|---|
| rail 내용 높이 | (1402) → 599 | (1558) → 688 | (1610) → 650 | (1666) → 776 | (1922) → 862 |
| rail, Advanced 펼침 | 2236 → 1192 | 2017 → 1078 | 3657 → 1350 | 2672 → 1559 | 2261 → 1269 |
| help `?` icon | 35* → 0 | 22* → 0 | 62* → 0 | 42* → 0 | 34* → 0 |
| text leaves, 1st screen | 56 → 44 | 38 → 34 | 55 → 38 | 45 → 44 | 39 → 32 |

\* baseline은 glyph 지표(버튼 hover 포함). P3의 같은 지표는 전 scheme 7 = About·SABES·preset 버튼 hover 영역.
plot 상단·iframe·expander는 P2에서 바뀐 뒤 그대로 (238·229·229·229·251 px, 0, 0).

- **knob 1행 = 34 px + 10 px 간격.** rail 전체 평균은 65–100 px/control (baseline 148–234) — 분모에 안 들어가는
  브랜드 띠(44)·그룹 머리글(26×4)·Advanced 버튼(40)이 SAS처럼 knob 적은 scheme에서 평균을 끌어올림. 목표 60은 knob 행만 충족.
- **첫 로드 시간**: 같은 기계에서 P2 트리와 P3 트리를 각각 띄워 cold 3회 — P2 중앙값 2871 ms, P3 3545 ms (+0.7 s).
  scheme 전환 시간은 변화 없음(rydberg 1725→1706, magneto 1726→1720, fwm 1728→1714). 컴포넌트 런타임 1회 비용으로 판단.
  (audit의 sas 5267→8813은 streamlit 서버 cold start까지 포함한 값이라 과장됨.)

**검증**
- 브라우저 CDP로 실제 입력: 드래그 +60 px → 0.50→1.10, `1.25 W` 입력 → 2.00 + "clamped to 2",
  `3 MHz` 입력 → 값 유지 + "expected mW", track 80 % 클릭 → 1.60, 복원 점 → 0.50.
  키보드 ↑↑·PgUp·Shift+↓·End·Home 모두 동작하고 포커스 유지, `aria-valuenow/min/max/text`·`aria-describedby` 채워짐.
- navigate-only knob(FWM δ) 0.7 s 드래그(125 ms 간격 12회 전송) → 실제 rerun 1회. Streamlit이 합침.
- About·Export·More·탭·Overlay·Show advanced·프리셋 클릭 smoke → `stException` 0. FWM Biphoton 뷰(그림 위 knob) 정상.
- pytest: 작업트리 전체 2027 통과 · 2 실패 → 하나는 기존(`test_docs_consistency`, 다른 세션이 지운 파일),
  하나는 내 탓(모듈 이동으로 깨진 CSV 복구 계약) → 고침. 커밋 후보 트리 재검증: 1542 통과, 실패 36개 전부 기존 `tests/quantum`.

**남은 일·주의**
- P3 커밋 대기 (사용자 승인).
- clamp·거부 알림은 한 번 보이고 다음 rerun에 사라짐(pop). 오래 남겨 낡은 정보가 되지 않게 한 선택.
- checkbox `?` 3개(sas 2·fwm 1)는 유지. 라벨 자체가 클릭 영역이라 hover 대체가 어려움.
- 폰(390)은 P5: rail 접힘 상태, plot 상단 349–529 px.

## 2026-09-22 · P2 셸 레이아웃
**한 일**
- P1 커밋 `4cb011f`. 사전 검증: HEAD 단독 export와 HEAD+P1 export 각각 전체 pytest → 실패 집합 동일
  (둘 다 `tests/quantum` 감사 36개 — 커밋 안 된 다른 세션 캐시·테스트에 의존), P1 신규 실패 0.
- 셸 (`decisions.md` D8): 상단 바(브랜드·scheme·About·regime/preset·Guide·SABES) + rail 브랜드 띠 · readout strip(칩·+N more) ·
  plot 카드(밑줄 탭·캡션·Overlay 토글·Export popover) · More pills · Advanced 인라인 공개.
- 분리: `gabes_ui/export.py`(지연 생성 다운로드), `gabes_ui/guide.py`(v2 컴포넌트, iframe 제거), `gabes_ui/readout.py`(순수 함수).
  `streamlit_app.py` 2007(P0) → 1072줄. carousel·hero 카드·헤더 카드 코드 삭제.
- 테스트: `tests/test_ui_readout.py` 10개, `tests/test_sabes_page.py` 라우터 순서 검사를 새 구조 기준으로 갱신.
- `tools/ui_audit.py`: 버튼형 Advanced 토글 후 ready 대기, scheme 전환 4회 재시도, 폰은 sidebar 닫고 상단 바에서 전환, 실패 시 `_failure.png`.

**함정 (재발 방지)**
- `st.html("<style>…</style>")`(style-only) → event container에 빈 element만, CSS 미적용. `st.markdown` + 숨김 anchor로 복귀.
- key 준 container는 `stLayoutWrapper`가 감쌈: 폭·bleed·sticky는 wrapper에, `max-width: none` 없으면 폭 확장 무시됨.
- 버튼형 토글은 rerun 유발 → 자동화는 고정 sleep 대신 ready 대기 필요.
- 폰에서 sidebar 열면 상단 바(scheme 전환기)를 덮음.

**측정** `--compare baseline p2` (1440×900 / 390×844):

| 지표 | sas | lambda | rydberg | magneto | fwm |
|---|---|---|---|---|---|
| plot 상단 (desktop) | 478 → 238 | 385 → 229 | 423 → 229 | 406 → 229 | 377 → 251 |
| plot 첫 화면 완전 노출 | F → T | T → T | T → T | T → T | F → T |
| rail 높이 | 1402 → 900 | 1558 → 1036 | 1610 → 1091 | 1666 → 1128 | 1922 → 1400 |
| rail, Advanced 펼침 | 2236 → 1633 | 2017 → 1471 | 3657 → 2029 | 2672 → 2178 | 2261 → 1757 |
| help `?` | 35 → 17 | 22 → 17 | 62 → 21 | 42 → 23 | 34 → 27 |
| iframe | 2 → 0 | 1 → 0 | 1 → 0 | 1 → 0 | 1 → 0 |
| main expander | 4 → 0 | 2 → 0 | 6 → 0 | 4 → 0 | 5 → 0 |
| plot 상단 (폰) | 732 → 529 | 511 → 349 | 586 → 416 | 638 → 480 | 494 → 390 |

- 초기 전송 1680 → 1709 KB (P1 1752 — iframe 2개 제거분 상쇄). ready ms: lambda 3146 → 3846 외 불변 — lambda는 측정마다
  2424–3846 요동(첫 전환에 solve 포함), 악화 판단 보류 → P3에서 3회 중앙값.
- 상호작용 스모크(headless, 예외 0): About · Export · More pill · plot 탭 · Overlay 패널 · Show advanced · OD 프리셋.
- pytest: 1975 passed · 3 skipped · 1 failed (기존 `test_docs_consistency`, 다른 세션 삭제 파일).

**design-critique 요약** (p2 스크린샷)
- 첫인상: 시선 순서 plot → hero 값 → 주의 칩. 목적 명확, 크롬 조용함. 목표 달성.
- 지금 고침: 상단 바가 오른쪽 30 px 못 미침(wrapper `max-width`) · "+2"가 strip 밖에 떠 있음 → 마지막 셀로 ·
  폰에서 "+2" 세로 기둥 → strip 아래 한 줄. magneto에서 상단 바 줄바꿈 → scheme 칸 310 px, 끝 묶음 nowrap.
- 이관: `?` glyph 17–27 · slider 끝점 캡션 별도 행 · rail > 900 px (P3) / 상태 칩이 대상 값과 멂 · SAS strip에 근거 진단값이
  Doppler·phase보다 앞섬 · 캡션이 rail 값 반복 (P4) / 폰 plot 상단 349–529 · 터치 타깃 36–40 px (P5) /
  Export 단일 곡선 이름 "primary: line N" (ux-copy, P5).
- 대비: muted #5B6B80 5.3:1, 주의 칩 6.4:1, 캡션 12 px muted — AA 통과. 탭·popover·토글 모두 실제 버튼 요소.

**남은 일·주의**
- P2 커밋 대기 (사용자 승인). Overlay 토글을 끄면 업로드한 CSV 재업로드 필요 (D8, 수용).

## 2026-09-21 · P1 토큰 · 타이포 · 포맷
**한 일**
- git 정리 (P0 승인분): stale staged 40개 unstage, lock 2개 삭제. P0 커밋 `ae456ab` (임시 index, 내 파일만).
- 폰트: IBM Plex woff2 12종 + `OFL.txt` → `static/fonts/` (207 KB). `config.toml` `[[theme.fontFaces]]` 12개, unicodeRange로
  필요한 부분집합만 로드. Streamlit static handler는 `.woff2`를 안전 확장자로 서빙 (text/plain 강제 없음) 확인.
- 토큰: `gabes_ui/theme.py` (LIGHT, iframe용 DARK 최소) → `--g-*`. `streamlit_app.py`의 인라인 CSS 458줄 →
  `assets/ui/gabes.css`. `streamlit_app.py` 하드코딩 hex 색 0개 (carousel·Guide 런처·BETA 배지도 토큰).
- 장식 제거: 헤어라인, hero 좌측 테두리+그라데이션, 그룹 색 점 (`GROUP_STYLES`·`METRIC_STYLES`·`_concept_style`·`_metric_style` 삭제),
  로즈 BETA pill → 외곽선, Guide 그라데이션+이모지 → ghost+SVG.
- 숫자: `gabes_ui/format.py` (`decisions.md` D7). slider `40.00→40`, `75.00→75.0`, `0.00→0.0`, `0.50` 유지.
  첫 규칙(기본값 자릿수 무제한)이 Rydberg `RF transition dipole` 계산값 1326.2572434514343을 6자리로 보여 → step+2 상한으로 수정.
- `_render_param`의 format import는 함수 안 (tests/test_fwm_excess_noise.py가 이 함수를 AST로 떼어 실행).

**측정** `--compare baseline p1` — 레이아웃 지표는 거의 불변 (P1은 레이아웃 아님):

| 지표 | sas | lambda | rydberg | magneto | fwm |
|---|---|---|---|---|---|
| plot 상단 (desktop) | 478 → 489 | 385 → 375 | 423 → 412 | 406 → 396 | 377 → 366 |
| rail 높이 | 1402 → 1389 | 1558 → 1545 | 1610 → 1600 | 1666 → 1710 | 1922 → 1909 |
| rail, Advanced 펼침 | 2236 → 2272 | 2017 → 1998 | 3657 → 3783 | 2672 → 2707 | 2261 → 2242 |
| plot 상단 (폰) | 732 → 732 | 511 → 504 | 586 → 599 | 638 → 633 | 494 → 490 |

- ±10–45 px 변동 = 글꼴 metric 차이(Source Sans → Plex)와 그룹 헤더 small caps 여백. 목표 달성은 P2·P3 몫.
- 초기 전송 1680 → 1752 KB (+72 KB, 폰트 부분집합 5개), 요청 67 → 73. 예산(+150 KB) 안.
- help glyph·iframe·expander 불변 (P2·P3 대상).
- SABES 페이지(`?app=sabes`) Plex 적용, 예외 0.

**pytest**: 1965 passed · 3 skipped · 1 failed (P0와 같은 `test_docs_consistency` — 다른 세션이 지운 `FWM_physics.tex`).
신규 `tests/test_ui_format.py` 23개 포함.

**남은 일·주의**
- matplotlib 글꼴은 DejaVu 유지 (P5 재검토). plot 제목 굵은 DejaVu가 Plex UI와 어긋남 → P2에서 앱 화면 제목 숨김으로 해소.
- ribbon "Half-height edges" 3줄 줄바꿈은 baseline과 동일 — P2 strip·P4 evidence 이동으로 해소 예정.
- P1 커밋 대기 (사용자 승인 후, 분리 worktree pytest → commit-tree).

## 2026-09-21 · P0 준비
**한 일**
- 현재 UI 진단(앱 실행·DOM 측정), 레이아웃 후보 6 · 입력 후보 5 · 시각 방향 3 평가 → B · Plex · ScrubField 확정 (`decisions.md`).
- 시안 canvas: https://claude.ai/artifact/C8d3dLBYDjH3ZJm84FPEQw (비공개). plot 이미지는 실제 GABES 출력(SAS·FWM 기본값), 제목만 제거.
  Before 보드는 markup 재현 — 실제 스크린샷은 이제 `audit/baseline/*.png`.
- `docs/ui_redesign/` 생성, CLAUDE.md "Where things live"에 포인터 1줄.
- `tools/ui_audit.py` 작성: Streamlit 자동 기동(빈 포트) → Edge headless + DevTools(tornado websocket) →
  5 scheme × (1440×900, 390×844) 측정 + 첫 화면 PNG → `audit/<tag>/{metrics.json,summary.md}`.

**측정** — `audit/baseline/summary.md`. 핵심:
- plot 첫 화면 완전 노출: SAS·FWM 실패. plot 상단 y 377–478 (desktop), 494–732 (폰).
- rail 1402–1922 px (1.6–2.1 화면), rail px/main control 148–234.
- help `?` glyph 22–62, iframe 1–2 (가이드 런처 + carousel), main expander 2–6.
- rail, Advanced 펼침: 2017–3657 px (Rydberg 최대).
- 초기 전송 1680 KB / 67 요청. scheme 전환 후 ready 1.7–5.3 s (SAS는 첫 로드 포함).
- 도구 검증: 2회 실행. 1회차 Advanced 토글 미검출(expander 텍스트 앞 아이콘 ligature) → 매칭 수정 후 재측정.
  레이아웃 지표 2회 동일, ready ms만 잡음 (lambda 2424 → 3146).

**정정** — 계획 단계 보고 오류:
- "`streamlit_app.py`·`sabes_page.py`에 UI 미커밋분(export +371줄)" → 틀림. 두 파일 worktree == HEAD.
  원인: git index에 옛 버전이 staged(`MM` 표시). 실제 UI 파일 미커밋 변경 0 → 계산 작업과 충돌 없음.

**git 정리 (사용자 승인, 같은 날 처리)** — stale 40개 `git restore --staged`, lock 2개 삭제(git 프로세스 없음 확인). 이전 상태 기록:
- staged 57개 중 40개가 stale: worktree == HEAD인데 index는 옛 내용. 이 상태로 `git commit` → 해당 40개 파일이
  옛 버전으로 되돌아감 (예: `streamlit_app.py` export 기능 삭제).
  복구안: `git restore --staged <40개>` (worktree 무변경, index만 HEAD로). 목록 재현:
  `for f in $(git diff --cached --name-only); do git diff --quiet HEAD -- "$f" && echo "$f"; done`
- `.git/index.lock`(2026-09-16 20:30, 0 B)·`.git/AUTO_MERGE.lock` 잔존 → index 쓰는 git 명령 실패.
  git 프로세스 없음 확인 후 삭제 필요 (Google Drive 동기 잔재).
- 나머지 staged 17개 + 미커밋 38개 파일 = 계산·문서 세션 작업. UI 리디자인은 손대지 않음.

**pytest** (공유 작업트리, 다른 세션 미커밋 포함): 1942 passed · 3 skipped · 1 failed.
실패 = `test_docs_consistency::test_repository_visibility_wording_is_consistently_public` —
다른 세션이 삭제한 `docs/FWM physics and analytic reconstruction/FWM_physics.tex`를 읽음. UI 작업 무관 (P0는 Python 모듈 무수정).

**P1 폰트 후보** (HEAD 요청으로 용량 확인, 미다운로드): jsDelivr `@fontsource/ibm-plex-sans@5` · `@fontsource/ibm-plex-mono@5`, OFL-1.1.
Sans latin 400/500/600/700 = 22.6/24.2/24.3/22.8 KB, latin-ext 400 = 16.0 KB, greek 400 = 9.9 KB; Mono latin 400/500 = 14.7/14.9 KB.
Plex Mono엔 greek 부분집합 없음 → 숫자 외 그리스 문자는 Sans/시스템 fallback.

**다음**: P1 — Plex woff2 다운로드 승인 → config.toml 테마 → CSS 토큰 파일 → 장식 제거 → 숫자 포맷.
