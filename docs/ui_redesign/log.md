# 개발 로그

최신이 위. 항목 = 날짜 · 단계 · 한 일 · 측정/검증 · 남은 일·주의.

## 2026-10-30 · Atomic 표시 통일 (D16)

직접 사용자 요청: 스킴별 Atomic 표시를 OD/SAS Atom / Isotope / Transition line 구조로 통일, `Nat` → `Natural`.

### ⟦AI⟧ 구현·검증

- `ParamSpec.atom_menu` + `controls._render_atom_menu`; SAS·Λ·Magneto가 선언. 저장 파라미터·계산 경로 무변경. 결정 [D16](decisions.md).
- `tests/test_ui_atom_isotope.py` 9 passed(SAS 기존 2 + Λ·Magneto·메뉴 도달성). UI·스킴 focused 18파일 347 passed.
- 내장 브라우저: SAS Rb/Natural/D1, Λ Rb/Natural/D1, Magneto Rb(비활성)/85/D1. Magneto 87 전환 시 재계산(PSN 감도 8.78→14.9 fT/√Hz). 콘솔·서버 오류 0.
- `ui_audit` 미실행: rail 행 수만 Λ·Magneto 각 +1, 측정 지표 변화 범위 밖.

## 2026-10-30 · 상단 관측값 조사·이름 가독성 개선

직접 사용자 요청: OD는 Peak OD 중심, SAS는 선폭보다 주파수 간격을 중시한다는 의견과 함께 모든 스킴의 과밀한 표시줄 조사. 작고 회색인 관측값 이름 개선.

### ⟦AI⟧ 구현·검증

- [조사·모드별 권고](audit/readout_20261030/summary.md), [23개 조합 표시 목록](audit/readout_20261030/inventory.json). 공통 자동 배치가 보조 4개를 노출하며, 상태 제외 후 SAS의 resolution-limited FWHM까지 주요값으로 승격함을 확인. 표시 항목 감량·SAS 두 공명 간격 구현은 권고이며 미적용.
- 공통 이름 스타일: 12→14.4px, 굵기600, 본문색, 줄바꿈 허용. 계산·결과·표시 항목 목록 무변경.
- readout/design 관련 pytest **29 통과**, 1.46s. 실제 FWM 데스크톱/모바일 및 23개 표시줄·102개 이름 모음의 두 화면 폭 검수: 잘림0·가로 overflow0. 전체 수치 회귀·GC 생략.
- 커밋 직전 이번 변경만 포함한 후보 디렉터리 검증: 기존 93개 표시·상호작용 검사 + Magneto/Rydberg readout-render 2개, **95 통과**, 17.13s. 다른 작업의 원자 모델·UI 확장 없이 통과.

## 2026-10-30 · 다른 스킴 상단 중복 문구 제거

직접 사용자 요청: 앞서 확인한 다른 스킴의 중복 문구도 제거.

### ⟦AI⟧ 구현·검증

- OD/SAS(원자종·Generic·보조 광학 뷰), Lambda(EIT/AT/CPT), Rydberg(EIT/AT), Magneto(현재 통합 모델·이전 단일 전이 참조)의 설정 반복 제목 제거. 축·단위·범례·결과값 유지. finite-IF 설명·PREDICTED 가정·Biphoton 출처/수렴 상태 무변경.
- `python -m pytest -q tests/test_ui_readout.py tests/test_ui_design_feedback.py tests/test_magneto_readout.py::test_readout_never_solves_or_mutates_raw_and_render_matches tests/test_rydberg_opacity.py::test_headless_and_figure_readouts_share_explicit_prediction_boundary --durations=3`: **31 통과**, 7.19s.
- 12개 선택 뷰의 제목·렌더 확인. 최초 확인에서 현재 Magneto가 별도 통합 renderer를 사용함을 확인하고 해당 제목도 제거; 앞서 통과한 8개 렌더 재사용, Magneto 4개 확인 완료. [검수 그림](../../output/playwright/scheme_headers_20261030/overview.png)의 바깥 제목은 검수용 구분명. 앱 캡션에는 없음.
- 표시 문구만 변경. 전체 수치 회귀·GC campaign 생략. 필수 검증 이월 없음.

## 2026-10-30 · FWM 라벨·그래프 상단 중복 문구 감량

직접 사용자 요청: Probe / conjugate gain → Gain, Squeezing Degree → Squeezing, Two-photon detuning delta → Two-photon Detuning. 그래프 위 설정·선택점 요약 제거. 다른 스킴 상단 과다 문구 점검.

### ⟦AI⟧ 구현·검증

- 기본 FWM 축·detuning 입력·Squeezing readout/표/설명 통일. 단위·Probe/Conjugate 범례·보정 모델별 gain 표기 유지. 새/이전 Squeezing 라벨 모두 기존 상태 안내 처리.
- 주 그림의 Δ/T/η/RIN 제목 제거. 빈 캡션 미생성. 확정 상태의 δ/Gain/Noise 반복문 제거; 드래그·갱신 중 미리보기 유지, Noise 표기는 Squeezing으로 통일.
- 다른 스킴 점검만 수행. OD/SAS의 원자·파워·온도·셀 길이, Lambda의 모드·Ω·완충기체·Doppler, Rydberg의 Ω·probe/LO, legacy Magneto의 셀·전이·QWP·광세기 제목은 설정과 중복. 감량 후보. Biphoton의 rate 출처/수렴 실패, Rydberg의 finite-IF·PREDICTED/가정 범위, Budker의 conditional DC 문구는 해석에 필요. 일괄 제거 부적절.
- 표시·상호작용 관련 pytest **93 통과**, 12.49s: `test_fwm_excess_noise`, `test_ui_design_feedback`, `test_ui_readout`, `test_plotcard_alignment`, `test_liveplot`, `test_update_status` 및 `test_fwm_generic`의 short-label/info/full-scan-render 3개. 최초 검사에서 새 라벨의 상태 안내 누락 발견·수정 후 통과.
- Chrome 1440×1000: 새 축/입력/readout 라벨·상단 캡션0·드래그 preview→ready 확인. 확정 후 요약문 공백, 앱 예외0·브라우저 오류0. 기존 Streamlit theme 경고 존재. [화면](../../output/playwright/fwm_labels_20261030.png).
- 계산식·엔진 무변경. 전체 수치 회귀·GC campaign 생략. 다른 스킴 감량 후보는 미수정.


## 2026-10-30 · UI 배포 후보 검수

직접 사용자 요청: 내장 브라우저 재검수 후 UI 변경 커밋·푸시·공개 사이트 반영.

### ⟦AI⟧ 검증

- 기존 HEAD에서 UI 변경 19개 파일만 구성한 독립 후보 검증. 다른 staged/물리 모델 변경 제외.
- 관련 pytest 161개 통과(22.60s): liveplot, scrub, focus, alignment, status, atom/isotope, isotope ratio, UI readout/design, CSV, Magneto noise. 전체 suite·GC는 범위 밖으로 생략.
- 내장 브라우저: Rb/87 선택·SAS 공명 확대, FWM 연속 드래그 마지막 33.6165 MHz 확정, 완료 문구 숨김, Magneto 범례 검수.

## 2026-10-30 · Magneto 범례·완료 상태 감량

직접 사용자 의견: ±1σ·command zero는 그래프 표현으로 충분하며, 동기화 완료 문구는 불필요.

### ⟦AI⟧ 구현·검증

- Magneto 범례에서 ±1σ·Command zero 제외. 음영·기준선·원래 artist/export 라벨 유지. offset이 있는 Physical zero는 구별을 위해 범례 유지. Ready 문구 비움; 갱신·미리보기·오류와 고정 상태 영역 높이 유지.
- `python -m pytest -q tests/test_update_status.py tests/test_magneto_noise.py -k 'browser_update or removed_plot or overlay_uses_axis_units or readout_controls_reuse_atomic_scan' --durations=3`: **5 통과**, 3.61s. 현재 범례 함수의 독립 렌더 확인 통과(두 항목 제외, 원래 artist 보존). 초기 import 누락은 해당 진행 중 모듈 추가 후 해소되어 재검사 완료. 계산 엔진 무변경, 전체 suite·GC 생략.

## 2026-10-30 · OD/SAS 원자·동위원소 선택 분리

직접 사용자 요청: Atom의 Rb 표기를 통일하고 아래 Isotope에서 Nat / 85 / 87 / Custom 선택.

### ⟦AI⟧ 구현·검증

- Atom: Rb / Cs. Rb Isotope: Nat / 85 / 87 / Custom. Custom에서만 혼합비 표시. 기존 저장 설정·순수 동위원소의 계산 입력·주파수 기준 유지하는 UI 매핑 적용.
- `python -m pytest -q tests/test_ui_atom_isotope.py tests/test_od_isotope_ratio.py --durations=3`: **21 통과**, 7.55s. AppTest 선택 전환·기존 87Rb 설정 복원 및 기존 isotope 수치 검사 확인. 계산 엔진 무변경으로 전체 suite·GC 생략.

## 2026-10-29 · 연속 입력·갱신 상태·그래프 배치·공명 확대

직접 사용자 요청: 네 가지 UX 개선을 구현하되 물리 계산을 유지하고, 적용 후 브라우저에서 다시 조작해 확인.

### ⟦AI⟧ 구현·검증

- 그래프 요청은 하나씩 전송하고 대기 중에는 마지막 입력을 보존. FWM의 선택 위치 삽입 샘플을 새 물리 조건으로 오인하던 문제 수정. 숫자 입력 초안·커서·드래그도 컴포넌트 재생성에 보존하며, 늦은 그래프 응답이 더 최근의 직접 입력을 덮지 않도록 순서 검사.
- Preview / Updating / Ready 상태를 구분하고 이전 결과임을 명시. 모드 전환으로 그래프가 제거돼도 Updating에 남지 않음. 캡션·확대 안내 공간을 예약하고 그래프 높이를 화면에 맞춤.
- SAS·Λ coherence·seeded FWM에 기존 샘플 기반 Resonance zoom. 축 범위·표식만 임시 변경하고 복원하며, 전체 범위 export·원본 배열·물리 계산 경로 보존. 성긴 샘플은 점으로 표시하고 확대가 해상도를 높이지 않음을 안내.
- `python -m pytest -q tests/test_liveplot.py tests/test_plot_focus.py tests/test_plotcard_alignment.py tests/test_ui_scrub.py tests/test_update_status.py tests/test_ui_readout.py tests/test_ui_design_feedback.py tests/test_experimental_csv.py --durations=10`: **112 통과**, 13.82s. 마지막 그래프/입력창 상태 이벤트·교체된 DOM 정리 순서 보완 후 `tests/test_liveplot.py tests/test_ui_scrub.py tests/test_plot_focus.py` 재검사 **49 통과**, 8.04s(회귀 검사 1개 추가).
- 내장 브라우저의 로컬 앱: 응답 전 두 번 드래그 후 마지막 값 확정, 갱신 중 Biphoton 전환 후 Ready, 직접 숫자 편집, EIT/SAS 공명 확대 확인. EIT 투과율 0.014·FWHM 0.46 MHz, SAS FWHM 18.87 MHz 유지. SAS 확대 전후 plot stage의 Y=387.55 px·높이=366.20 px 동일(1122×884 화면).
- 계산 시간 단축 배율은 측정하지 않음. 이번 검증은 입력 연속성·표시·export 보존 범위이며, 전체 수치 회귀·GC campaign은 변경 범위 밖이라 생략. 공개 서버 배포·커밋·푸시는 수행하지 않음.

## 2026-10-23 · 자기광학 상태 배지 감량

직접 사용자 요청: `sample-grid consistent`·`external validation required`·긴 hyperfine 제한 세 배지 과밀 해소. 삭제 또는 디자인으로 정리.

### ⟦AI⟧ 구현·검증

- 자기광학 상단에 중립 `Model estimate` disclosure1개. 기존 상태·설명은 클릭/키보드로 펼침. 수치 카드와 같은 줄로 배치해 기본 상태의 별도 경고 행 제거.
- 실제 조치 필요 시 `Check resolution`·`Check polarization`·`Check field range`·`Check detection` 등 한 개의 경고 배지. 색·기호·라벨 함께 사용. 모델 추정치를 실험 검증 완료로 표시하지 않음.
- 표시 계층만 변경. 계산·SNR·기존 상태 descriptor/export 보존. 다른 스킴은 기존 chip 배치 유지.
- `python -m pytest -q tests/test_ui_readout.py tests/test_ui_metrics.py tests/test_ui_design_feedback.py`: **31 통과**, 1.29s. 상세 보존·경고 우선순위·HTML escape·다른 숫자 셀 보존 확인.
- 실제 Chrome1440×900/390×844: 닫힘·펼침·Space 조작·rotation의 `Check resolution` 확인. 모바일 닫힘/펼침 가로 overflow0, 앱 예외0·브라우저 error0. 기존 theme warning 존재.
- [기본 화면](../../output/playwright/magneto_status_20261023_closed.png) · [상세 펼침](../../output/playwright/magneto_status_20261023_open.png) · [모바일](../../output/playwright/magneto_status_20261023_mobile.png).
- 수치 경로 무변경으로 전체 numerical regression·GC campaign 생략. 필수 검증 이월 없음.

## 2026-10-23 · SABES 상단 요약 리뷰·커밋

직접 사용자 요청: SABES 리뷰, 이슈 수정, 추가 검수1회 후 상단 변경 커밋·푸시.

### ⟦AI⟧ 검수

- 코드 리뷰·추가 검수에서 수정할 결함 없음. 페이지·NF 연동·문서 검사31개 통과(3.78s), 실제 브라우저 카드/탭/프로필 전환 및 1440/720/390/320px 검수 통과.
- 분리된 커밋 소스로 30개 통과·미커밋 canvas 검사1개 제외(1.46s). 이전 검증 로그의 중복 CR 줄바꿈만 정리.
- 상단 변경·해당 문서/audit만 분리. 이전 광학 테이블·이벤트 변경과 다른 staged 작업 보존. 계산 경로 무변경으로 전체 수치 회귀·GC 생략.
- [검수 기록](audit/sabes_readout_review_20261023/summary.md).

## 2026-10-22 · SABES 상단 요약 감량

직접 사용자 요청: 잡음 수준·Gain G_s/G_c·Detector clearance 세 카드. 나머지는 상세 탭·상태 배지로 이동, 잘못된 라벨 수정.

### ⟦AI⟧ 구현·검증

- `Noise level [dB]`는 계산 NF의 dB 값. Gain 두 값은 probe/seed G_s·conjugate G_c 순서로 묶고 `Pump gain` 제거.
- 제목 옆 모델 상태 배지. Detection에 선형 NF·계산 NF·algebraic floor·기존 RIN 진단, Power budget에 잔류 carrier 비율과 적용 범위 보존.
- 계산 경로·입력·수치 변경 없음. 1440px 한 줄 세 카드, 390/320px 세로 배치. 720px 화면 포함 카드 잘림·문서 가로 overflow 없음. 접힌 사이드바 열기 버튼과 제목 겹침 방지.
- [UI audit](audit/sabes_readout_20261022/summary.md): 전체 `python -m pytest -q` **2707 통과·3 skip**, 1191.07s. 최종 제목 여백 CSS 적용 후 관련31개 통과. 화면·수치·로그 수집 범위 기록.

## 2026-10-23 · 자기광학 코드 리뷰·수정

직접 사용자 요청: 자기광학만 코드 리뷰, 발견 사항1회 수정 후 문제 없으면 현재 브랜치 커밋·푸시.

### ⟦AI⟧ 검토·조치

- 검토 범위: Hanle/EIA/NMOR의 field reference·광학 readout·샷 노이즈/SNR·캐시·export·관련 검사. 공유 solver·타 스킴·GC 변경 제외.
- [P2] `magneto_noise.py`의 양수 광량 검사 때문에 기존0-drive 계산의 readout이 ValueError. 광량/효율0을 검출 photon budget0으로 처리하고 noise/SNR을 unavailable로 반환하도록 수정. 음수·비유한 입력은 계속 거절.
- [P2] SNR·각도 등을 정의할 수 없을 때 상태 hero의 대체값으로 `unavailable` 숫자 카드가 표시됨. 자기광학의 unavailable metric을 status로 지정해 chip으로 표시. 타원/원편광·physical zero 밖·광량0 경계 검사 추가.
- 한 차례 수정 후 `python -m pytest -q tests/test_magneto_noise.py tests/test_magneto_readout.py tests/test_magneto.py tests/test_magneto_doppler.py tests/test_ui_readout.py --durations=10`: **82 통과**, 9.55s. 기존 export UTC 경고2건. 추가 미해결 결함 없음.
- 다른 미커밋 변경을 제외한 커밋 후보에서도 위82개 **14.44s 통과**, `python -m pytest -q tests/test_schemes_render.py -k 'Magneto or magneto' --durations=5` **7개·2.63s 통과**. 최종 관련89개. 이후 문서에 결과만 보충, 검증된 코드 보존. [검증 receipt](../../output/magneto_review_20261023/validated.json).
- 전체 suite·GC campaign 생략: 국소 자기광학 readout 수정, 공유 수치 엔진·GC 계약 무변경. 커밋에는 자기광학 코드/검사·readout cache 버전·관련 문서 부분만 포함.

## 2026-10-23 · 자기광학 상단 SNR

직접 사용자 요청: Transmission/rotation에서 신호가 샷 노이즈보다 몇 배 큰지 SNR 표시. 상단 관측량 중 덜 유의미한 값 하나 교체.

### ⟦AI⟧ 구현·검증

- 중복도가 높은 `command B=0` 값을 상단에서 빼고 `Peak shot-noise SNR`로 교체. 기존 command-zero 값은 Derived quantities·JSON 보존. 관측량 셀 수 유지.
- 표시 기준으로 **스캔 내 최대 진폭/RMS 비율** 채택: Transmission `max |T−T_bg|/σ_T,bg`, rotation `max |θ|/σ_θ(B)`. 같은 B의 신호·잡음 대응. `×` 값·Above/At/Below shot noise 표시. 1×는 1 RMS 비교이며 검출 보장 임계값 아님.
- 모의 지그재그·표시 토글과 독립. ENBW·검출효율·광량 변경에 갱신. 타원/원편광 rotation·검출광량0은 unavailable. 식·최대 지점·각 B의 비율·범위·조건 JSON 보존. Readout cache만 갱신; 원자 solve 경로 무변경.
- `python -m pytest -q tests/test_magneto_noise.py tests/test_magneto_readout.py tests/test_magneto.py tests/test_ui_readout.py --durations=10`: **69 통과**, 7.09s. 기존 export UTC 경고2건. 배경 DC 배제·B별 잡음·0/저/고 SNR 표시·ENBW 역제곱근·토글 독립·export·기존 field-reference 확인.
- `python -m pytest -q tests/test_schemes_render.py -k 'Magneto or magneto' --durations=5`: **7 통과**, 2.66s. 공용 readout·자기광학5 regime/offset 계약 유지. `git diff --check` 대상 파일 오류0.
- 검출광량0 SNR 보류 추가 확인: `python -m pytest -q tests/test_magneto_noise.py -k no_detected_light_has_unavailable_snr` **2 통과**, 2.99s. 총 관련78개 통과; 소스 무변경 상태에 검증 추가.
- 실제 Chrome: 기본 Transmission **1172×**, rotation **2245×** 확인. 이전 답변의 636×는 Transmission의 **B=0** 값; 새 표시는 스캔 최대값. 1440×900·390×844 SNR/상태 표시 정상, 모바일 가로 overflow0·앱 예외0·브라우저 error0. 기존 theme warning 존재.
- [Transmission](../../output/playwright/magneto_snr_20261023_transmission.png) · [rotation](../../output/playwright/magneto_snr_20261023_rotation.png) · [모바일](../../output/playwright/magneto_snr_20261023_mobile.png).
- 전체 suite·GC campaign 생략: 자기광학 readout·캐시 표시 변경에 한정. 공유 solver·GC 모델·입력·검증 계약 무변경. 필수 검증 이월 없음.

## 2026-10-23 · 자기광학 샷 노이즈 오버레이

직접 사용자 요청: Transmission/rotation의 y축 단위로 noisy한 샷 노이즈 기준선 동시 표시. 신호 관측 가능 수준 판단 지원.

### ⟦AI⟧ 구현·검증

- `Show shot noise` 기본 ON. 주황색 고정 seed 모의선·±1σ 음영·RMS 값 표시. Transmission 기준은 대칭 physical-field wing 배경, rotation은 0 mrad·각 B의 투과 광량. 시각적 진폭 확대 없음.
- Advanced → Shot noise: 유효 균일 빔 지름2 mm·검출효율0.8·ENBW1 kHz를 가정 기본값으로 채택. 광세기×면적으로 입사파워 산출. 모두 `recompute=False`; 원자 계산·기존 신호 보존.
- 단측 전기 ENBW의 coherent-light 광자 샷 노이즈. 고정 무잡음 입사 기준 정규화·균형 편광계 가정. 타원/원편광 rotation과 저광자수 Gaussian 샘플은 표시 보류. 전자·기술·원자 spin 잡음 제외. 조건·seed·곡선·RMS export 보존.
- `python -m pytest -q tests/test_magneto_noise.py tests/test_magneto_readout.py tests/test_magneto.py tests/test_magneto_doppler.py --durations=10`: **59 통과**, 9.96s. 기존 export UTC 경고2건. 독립 Poisson 계수 유도·전력/효율/ENBW 스케일·mrad 변환·토글·export·불투명 극한·UI 입력 확인.
- `python -m pytest -q tests/test_schemes_render.py -k 'Magneto or magneto' --durations=5`: 자기광학5 regime·offset 화면 계약 **7 통과**, 2.93s. `git diff --check` 대상 파일 오류0.
- 실제 Chrome: ON/OFF·rotation 전환·ENBW1→4 kHz에서 RMS4.43e-6→8.87e-6 확인. 1440×900·390×844 표시 확인, 모바일 가로 overflow0·앱 예외0·브라우저 error0. 기존 theme 빈 색상 warning 존재.
- [Transmission 화면](../../output/playwright/magneto_snl_20261023_transmission.png) · [rotation 4 kHz 화면](../../output/playwright/magneto_snl_20261023_rotation.png) · [모바일](../../output/playwright/magneto_snl_20261023_mobile.png).
- 전체 suite·GC campaign 생략: 자기광학 readout/표시 추가에 한정, 공유 solver·GC 모델/입력/검증 계약 무변경. 필수 검증 이월 없음.

## 2026-10-22 · 기본 conjugate gain·이전 토글 제거

직접 사용자 요청: `Measured-data reference only` 토글 제거. Gain calibration OFF에서도 conjugate gain 기본 표시.

### ⟦AI⟧ 구현·검증

- Advanced에서도 이전 토글 숨김. 내부 분석·저장 파일 호환용 입력은 유지.
- Probe·conjugate와 범례 상시 표시. OFF conjugate는 주황색 점선. 기존 커서·선택점 순서 유지, gain·NF 계산식 보존.
- 관련47개 통과. 전체 `python -m pytest -q`: **2707 통과·3 skip**, 1174.02s. 기존 UTC 경고13개.
  작업 소스3개 전후SHA 동일. [검증 receipt](../../output/playwright/fwm_conjugate_20261022/validation.json).
- [UI audit](audit/fwm_conjugate_20261022/summary.md): desktop/phone 정상, FWM Advanced 높이1652→1588px.
  [Calibration OFF 화면](../../output/playwright/fwm_conjugate_20261022/calibration_off.png)에서 두 곡선·범례 확인, 실제 브라우저 오류0.

## 2026-10-22 · DEFAULTS·설정 저장

직접 사용자 요청: 중복 lab reference 제거, 두 기본값 드롭다운 전환, 이름 붙인 브라우저·로컬 저장/가져오기.
보정명 `Empirical gain calibration`, 이전 측정비교는 화면에서 제거.

### ⟦AI⟧ 구현·검증

- `DEFAULTS`에 `2026-10-13 experiment`·`Sim Gold Reference`. 모드·숨김 입력까지 전체 복원.
- `Save as default`: 이름 지정 → 브라우저 localStorage 또는 로컬 JSON 다운로드. `Import`: JSON 검증 후 복원.
- 저장소 오류·잘못된 파일은 현재 설정 보존. 저장 응답 유실은 재전송, 이전 선택 응답은 새 드롭다운에 적용되지 않도록 분리.
- 이전 비교는 [내부 분석 CLI](../../analysis/fwm_experiment_20261022/measured_reference_comparison.py)로 보존. 물리 계산식·GC 변경 없음.
- [UI audit](audit/fwm_defaults_20261022/summary.md): 5개 스킴 desktop/phone 사이드바 왕복 정상, desktop 플롯 5/5 첫 화면 유지.
  FWM 플롯 상단 desktop 260→253px·phone 407→368px. 저장 메뉴 조건부 표시, 좁은 rail 줄바꿈 확인.
- 전체 `python -m pytest -q`: **2704 통과·3 skip**, 1188.18s. 수집 후 최종 이름 선택 동기화 수정은 UI 두 파일에 한정;
  최종 관련 **97 통과**, 22.63s. 기존 UTC 경고. 다른 세션의 SABES 두 파일 변경은 receipt에 별도 기록.
  임시 검증 복사본 중복 수집은 `pytest.ini`의 output 제외로 해결.
- 실제 브라우저: 저장·새로고침 후 목록 복원·선택·로컬 다운로드·가져오기·같은 수치의 사용자 이름 유지 확인.
  [검증 receipt](../../output/playwright/fwm_defaults_20261022/validation.json), [저장 메뉴](../../output/playwright/fwm_defaults_20261022/save_options.png).

## 2026-10-22 · SABES 디자인 규약·모델 별표

직접 사용자 요청: 기존 디자인 규약 동일 적용, 광학 테이블 표현 기준 강화.
반복 `lumped`를 `*`로 교체. 선 굵기 설명 대신 별표 뜻을 범례에 표시.

### ⟦AI⟧ 반영·검증

- [디자인 규칙](../SABES/design_rules.md), D13, 기존 광학 테이블 계획 연결. glyph·광로 anchor·라벨 겹침/잘림·선택·좁은 화면 기준 정리.
- `*`는 상세 부품 모델 미구성. 투과율·기하 일부 적용과 양립하며, 기종 미상·클릭 가능 여부와 독립. 이름 없는 mirror의 단독 별표 라벨은 생략.
- 그림·소자 선택 목록·상세 제목·tooltip의 반복 `(lumped)` 제거. 범례 설명 교체, 보조 문구 대비 보완. 광로·계수·계산·광선 식별색은 그대로.
- [UI audit](audit/sabes_asterisk_20261022/summary.md): Chromium 1440×900·390×844에서 Part1/Part2, 범례·선택·상세 제목 확인. 별표 표시와 상세창 닫기 정상.
- 전체 `python -m pytest -q`: **2570 통과·3 skip**, 1025.09s. 이후 중복 상세창 캡션만 제거하고 관련 테스트 **73개 재통과**, 6.00s. 최종 상세창도 다시 확인.
- 문서 일관성4개 통과, scoped `git diff --check` 통과. 전체 회귀 종료 이후 다른 작업에서 수정된 비관련 테스트는 후속 검증 receipt에 별도 기록·보존.
- 새 품질 규약은 이후 테이블 설계·검수의 기준. 이번 패치는 표기 정리이며, 기존 모든 배치의 전면 재설계 완료를 뜻하지 않음.

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
