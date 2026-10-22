# 개발 로그

최신이 위. 항목 = 날짜 · 단계 · 한 일 · 측정/검증 · 남은 일·주의.

## 2026-10-22 · 접힌 사이드바 다시 열기

직접 사용자 요청: 사이드바를 접으면 다시 열 수 없는 문제 확인·수정.

### ⟦AI⟧ 재현·수정·검증

- 폭 448px로 확장 후 접으면 기존 화살표가 화면 밖(x≈−123px)으로 이동. 실제 열기 버튼도
  상단 toolbar 숨김에 포함. 고정 312px 이동 우회 제거, Streamlit의 실제 열기 버튼 복원.
- 열기 버튼 44×44px, 상단 바와 겹치지 않도록 여백 확보. 창·사이드바 폭에 독립적으로 표시.
- [UI audit](audit/sidebarfix_20261022/summary.md): 5개 스킴 모두 desktop 접기/다시 열기·phone 열기/접기 통과.
  포인터 hit-test·실제 클릭·상태 전환을 검사하도록 감사 도구 보완. 폭 448px·Enter 재열기도 확인.
- 기존 audit 대비 플롯 위치·크기 동일. 800/390/320px에서도 버튼 접근·가로 overflow 없음 확인.
- `python -m pytest -q --ignore=output`: **2570 통과·3 skip**, 1034.36s.
  임시 검증 사본 제외; 기본 명령의 수집 대상도 2573개로 일치 확인.

## 2026-10-21 · 외부 웹 디자인 피드백 재검증·패치

직접 사용자 요청: 외부 AI 피드백을 1회 독립 재검증하고 확인된 사항 패치.
안내서 꾸밈(첨부 D11)은 의도적으로 유지. 기존 결정 번호와 충돌하여 D12에 기록.

### ⟦AI⟧ 독립 재검증·반영

출처: [외부 AI 피드백](<C:/Users/박수인/Downloads/gabes_web_design_feedback.md>).
문서의 지시·규약 초안 전체를 사용자 정책으로 채택하지 않고 코드·실제 Edge 화면과 대조.

| 항목 | 확인·반영 |
|---|---|
| FWM 경고 위계 | 11.2px 보조줄·칩 0개 확인. 값에 조건부 NF 경고 칩 부착. |
| 지표명 | 기존 사용자 지정 `Squeezing Degree` 유지. 값과 중복되는 `[dB]`만 HTML 라벨에서 생략. Guide 용어 통일. |
| 이득 정밀도 | 정량 gain 미지원인 모델값은 ≈·3 유효숫자 표시. 지원/실측참조 값과 원본 계산·내보내기 보존. |
| Guide 중복 | 배지 줄 제거, 카드 4개는 서로 다른 탐색 링크로 변경. FWM 해석은 07장에 모으고 다른 곳은 링크. |
| Guide 정확성 | STEP 1을 상단 Scheme 선택으로 수정. classical RIN a·계산 b·after-cell/QE 설명을 현재 코드에 맞춤. |
| 스킴 이름·800px | 표시 이름만 축약, 내부 식별자·저장 상태 유지. 800px FWM 상단 1행, 버튼 접근성 이름 유지. |
| HTML 기호 | `eta_ext`, `G_s/G_c`, `η_probe` 하첨자 복원. Matplotlib ASCII 규칙 유지. |
| 읽기 전용 값 | disabled slider 3개 → 값·단위 행. 상태·계산 근거는 hover, 실제 unavailable 사유는 표시. 입력/solve key 제외 유지. |
| 사이드바 기본 접힘 | 새 1440px 세션에서 `aria-expanded=true`. 재현 안 됨. 사용자 접기 상태·기본 auto 유지. |
| 글자·대비 | 앱 보조줄·배지·그룹 제목, Guide 작은 라벨 12px 이상. 대비 부족한 Guide 텍스트 3곳 수정. |
| 진단 표 | Claims / Solver / Geometry & coupling / Passivity로 묶고 claim 이유를 목록화. 행·수치·이유 전문 보존. |
| 딥링크 | `?scheme=<registry name>` 지원, 수동 변경 시 URL 동기화. 유효하지 않은 ID 무시. 기존 임의 해시를 공식 링크로 간주하지 않음. |
| Guide 같은 탭 | 기존 버튼도 Blob 새 탭으로 열림 확인. 동작 유지. 개발 서버의 오래된 static ETag 문제는 문서 revision URL·재검증 fetch로 보완. |
| 안내서 꾸밈 | 그라데이션 6개·카드 4개·색 callout 유지. 요약은 강조색, 실제 한계는 경고색 사용. |

- 문구 블록 판정: MERGE 7 · REWRITE 5 · DELETE 1. UTF-8·원본 mixed newline 유지.
  FWM 그림의 기존 gain-noise 진단 표기도 확인하여 해당 1장만 현재 기본 NF로 재생성
  (`−8.19 dB`, gain `394.88`). 다른 8장 보존, static guide 재빌드.
- `python -m pytest -q`: **2542 통과·3 skip**, 955.77 s. 최신 UI/CSV/FWM 관련 추가 **120 통과**,
  Guide 재빌드 후 docs consistency **4 통과**. 계산 경로 변경 없음.
- 커밋 분리 검증(2026-10-22): 다른 미커밋 작업을 제외한 예정본에서 관련 **115 통과·1 skip**.
  기본 FWM `−8.19 dB`·gain `394.88`, Guide 소스/그림과 배포 HTML의 빌드 일치 확인.
- [전](audit/feedback_before_20261021/summary.md) / [후](audit/feedback_after_20261021/summary.md):
  desktop plot 5/5 첫 화면 유지. FWM rail content 1372→1099px, plot top 251→260px.
  전환 시간은 단일 표본이므로 속도 개선 근거로 사용하지 않음.
- [Guide 검사](../../output/playwright/design_feedback_20261021/guide_after.json): 일반 HTML 텍스트 381곳
  대비 부족 0·12px 미만 0, 내부 anchor 끊김 0. 수식 첨자·로고·그림 문자는 해당 자동검사에서 제외;
  그라데이션 텍스트는 별도 시각 확인. 대비 기준: [W3C SC 1.4.3](https://www.w3.org/WAI/WCAG22/Understanding/contrast-minimum.html).
- 실제 화면: [FWM 1440](../../output/playwright/design_feedback_20261021/fwm_1440.png),
  [800](../../output/playwright/design_feedback_20261021/fwm_800.png),
  [Guide](../../output/playwright/design_feedback_20261021/guide_1440.png).
  `?scheme=fwm`, 수동 `?scheme=sas` 전환·Biphoton 전환·Guide 새 탭 검증.

## 2026-10-21 · OD 표시·Overlay·드래그 미리보기

직접 사용자 요청: Transmission만 임시 표시, 팝업 없는 Overlay 개선,
OD 정렬·FWM TPD 드래그 중 즉시 미리보기. 반영 완료.

### ⟦AI⟧ 구현·검증

- 표시 flag로 OD 추가 뷰만 숨김. 물리 계산·기존 내보내기 유지.
  Overlay 설정은 아래 접이식 영역. 그래프 아래 X/Y 미리보기 조정·직접 곡선 드래그 지원.
- 고정 Matplotlib 배경 + SVG 레이어. FWM rail·그래프 양쪽에서 TPD 선택.
  드래그 중 서버 전송 0회, 종료 시 1회 확인. 확정 readout·정렬·내보내기 동기화.
- 비균일 축·반전·relative Y 변환과 Python 결과 일치. canonical NF 보간 후 dB 변환,
  invalid gap 유지. 절대 교정 데이터의 live Y 조정 차단.
- 개발 중 Streamlit 모듈 교체 뒤 CSV dataclass pickle 실패 재현. 캐시 함수가 현재 모듈의
  loader를 조회하도록 수정·별도 프로세스 회귀 검증.
- `python -m pytest -q`: 2494 통과·3 skip (922.64 s).
  전체 수집 시작 후 추가한 신규 11개 포함 UI/CSV 관련 최종 101 통과 (16.21 s).
- [UI audit](audit/liveplot_20261021/summary.md): desktop 기본 플롯 5/5 첫 화면 노출, iframe 0.
  390×844에서 live 컨트롤 1열·가로 overflow 없음. 실제 Edge 검사·수학 검증 통과.
- 브라우저 자료: [OD 드래그](../../output/playwright/live_patch_20261021/od_gesture.json),
  [FWM 드래그](../../output/playwright/live_patch_20261021/fwm_gesture.json),
  [OD desktop](../../output/playwright/live_patch_20261021/od_desktop.png),
  [OD mobile](../../output/playwright/live_patch_20261021/od_mobile.png).

## 2026-10-07 · AutoOD 기반 OD 비교·조정 개선

- AutoOD 기능 대조: 임의 Rb 혼합비, 큰 X 배율 자동 보정, Y 배율·이동,
  갱신 지연. 기존 GABES의 상대/절대 투과율 구분·잡음 보정 유지.
- 업로드 즉시 X span 자동 프레이밍. Base X scale + 고정 0.1–3 미세 배율.
  비균일 샘플도 transform_x와 같은 양끝 중점 사용. 초 단위 CSV 지원.
- 조정 위젯 form 적용: 드래그 중 서버 갱신 없음; Apply 한 번만 그림 갱신.
  CSV·내보내기 fragment 분리. 원본 그림 복제해 반복 적용 시 곡선 누적 방지.
- Fit alignment: 고정 이론 곡선에 X affine 정렬·상대 Y gain/offset 추정.
  절대 투과율 Y 변환 차단. 주파수 교정·물리 파라미터 추정과 구분.
- Natural/Custom Rb 조성·고정 ⁸⁵Rb 주파수 원점. 기존 기본·순수 메뉴 유지.
- 독립 리뷰: 0 배율 후 Fit 예외 발견·수정; 입력/변환 유효성 검사 추가.
  합성 초 단위·300배 신호·역방향·비균일 축·절대 교정·출력 일치 회귀 추가.
- 무관한 잡음·램프·사인 신호 성공 오판 수정: 최종 uncapped normalized RMSE
  ≤0.35·상관 ≥0.85 필요. 최적화 중 robust loss와 최종 성공 판정 분리.
- 관련 회귀 99통과·1skip; 다른 화면 회귀 89통과. 별도 리뷰 추가 지적 없음.
- 실제 브라우저: 초 단위 업로드 자동 프레이밍, 드래그 중 이미지 무변경,
  Apply 후 갱신, 합성 OD 정렬 RMSE 0.00126·100% 겹침 확인.
  Custom 0.35 입력 → ⁸⁵Rb/⁸⁷Rb 35/65% 표시·계산 확인.
- 최종 전체 `python -m pytest -q --durations=10`: 2204통과·4skip·실패0,
  981.17초. 실행 프로세스의 `OPENBLAS_NUM_THREADS=1`만 지정.
  검증 로그·화면: `.playwright-cli/od_final_pytest.txt`,
  `.playwright-cli/od_alignment_verified.png`.


## 2026-09-23 · 세로 데스크톱 대응 (P3 후속)
**한 일** (`decisions.md` D10, CSS만)
- P0–P3 네 커밋 `origin/main` 푸시 (`2998dc1..808c5a5`, 사용자 승인). 푸시 전 `git fetch`로 원격 선행 없음 확인.
- 1080×1920에서 발견한 진짜 버그 2개:
  - plot 캡션이 Overlay·Export 위로 겹쳐 그려짐. Streamlit markdown 블록이 내용 폭을 가져서 줄지 않음 → `min-width: 0` 고정.
  - readout strip이 두 줄로 접히면 행 사이에 구분선이 없고 행 끝에 빈 칸이 남음 → 구분선을 `gap: 1px`로 바꾸고 칩 칸도 늘어나게.
- ≤1439: scheme 이름 칸이 먼저 줄어들어 상단 바 한 줄 유지(FWM 88 → 52 px). ≤1200: 캡션이 자기 줄로.
- ≥1100 높이: 블록 간격·카드 여백만 소폭 확대.

**측정** (빈 세로 공간 = 마지막 내용 아래)

| 뷰포트 | scheme | 상단 바 | strip | plot | 빈 공간 | 캡션 겹침 |
|---|---|---|---|---|---|---|
| 1440×900 | SAS | 52 | 72 | 1079×516 | 69 | 없음 |
| 1080×1920 | SAS | 52 | 140 | 724×346 | 1149 | 없음 |
| 1080×1920 | magneto | 107 | 110 | 724×368 | 1101 | 없음 |
| 1080×1920 | FWM | 52 | 86 | 724×523 | 1026 | 없음 |
| 1200×1600 | SAS | 52 | 140 | 844×404 | 778 | 없음 |
| 1280×800 | SAS | 52 | 120 | 922×441 | -3 | 없음 |

- 1440×900 화면은 변화 없음(픽셀 동일), 폰 390은 오히려 개선(행 구분선 생김). 가로 스크롤·예외 0.
- 1080×1920에서 Overlay 패널·More 표를 열면 빈 공간 1140 → 462 px.

**남은 일·주의**
- 세로 화면의 빈 공간은 더 줄이지 않음. 이유는 D10 — 스펙트럼 비율 고정 + 그림이 이미 열 폭을 꽉 채움.
  카드 늘리기·세로 가운데 정렬 둘 다 시안을 만들어 보고 기각(각각 빈 흰 상자 / 상단 바와 560 px 분리).
- magneto 상단 바는 ≤1439에서 두 줄. regime 선택지 5개라 한 줄에 안 들어감.
- 더 채우려면 창 높이를 Python이 알아야 함 → 보고용 컴포넌트 + rerun 1회. 사용자 결정 대기.

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
