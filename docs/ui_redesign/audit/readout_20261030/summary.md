# 관측값 표시줄 전수 조사 · 2026-10-30

## 직접 사용자 요청

- 모든 스킴 상단 관측값의 과밀 여부 조사.
- OD는 Peak OD 하나면 충분하다는 의견. SAS는 FWHM보다 공명 간 주파수 간격을 중시.
- 작고 회색인 관측값 이름의 가독성 개선.

## ⟦AI⟧ 조사 범위·공통 원인

조사 시점의 작업트리 기준. FWM quadrature는 당시 별도 진행 중인 구현을 포함하며, 이번 표시 패치의 신규 기능이 아님.

등록된 5개 스킴의 주요 모드·분기 23개를 현재 renderer로 확인. [전체 표시 목록](inventory.json). OD/SAS의 약한/유한 probe, Lambda EIT/AT/CPT, Rydberg EIT/공명·비공명 AT, Magneto 두 스캔축×두 신호, FWM intensity/보정/위상, Biphoton 세 topology×두 모델 포함. 일부 스캔은 빠른 격자 사용. 수치는 UI 분기 확인용이며 물리 정확도·수렴 검증 아님. 추가 조건부 분기는 소스 확인.

- [공통 분류](../../../../gabes_ui/readout.py): 주요값 1–2개 뒤에 목록상 다음 수치 최대 4개를 자동 노출. 사용 목적·선택된 그래프와 무관. 상태 항목은 먼저 제외되므로, 상태가 hero여도 나머지 수치가 빈 주요값 자리를 채움.
- [앱 조립](../../../../streamlit_app.py): Magneto만 상태를 하나의 설명 패널로 접음. 나머지는 상태 칩을 그대로 표시. 수치가 넘칠 때에만 `+N more` 생성.
- 상세 이관 시 값·단위·출처·가용성은 유지하고, 상세 열기 경로도 항상 제공해야 함. 단순히 수치 목록을 삭제하는 방식은 부적절.

## ⟦AI⟧ 모드별 판정·권고 — 표시 항목 변경은 아직 미적용

| 모드 | 현재 노출·문제 | 기본 표시 권고 | 상세/조건부로 이동 |
|---|---|---|---|
| OD | A10, τ, Gaussian Doppler FWHM, Peak phase shift: 4개. A10·τ는 같은 투과율의 상수배 환산. Transmission만 보이는 화면에 phase도 노출 | **Peak OD 1개**. 현재 주값 기준이면 `Peak OD (A10)` 또는 이름+도움말로 −log10(T) 명시 | τ, Doppler 폭, phase. 유한 probe의 포화·모델 범위는 상세와 짧은 상태 표시 |
| SAS | 기본 설정에서 `resolution-limited`인데 Sub-Doppler FWHM이 주요값으로 자동 승격. Half-height edges, Samples/FWHM, Scan-edge distance, A10까지 5개. resolved이면 Lock Slope·Lock Detuning이 앞에 옴 | **선택한 두 공명의 Δν 1개**. 필요 시 두 공명 위치를 선택 UI에 표시 | FWHM·경계·샘플 수·스캔 끝 거리·lock slope. 간격 미선택/미분해 상태는 명시 |
| Lambda EIT/CPT | 공명 투과율, Window FWHM, group index: 3개 | **공명 투과율·window 폭 2개**. SAS와 달리 이 모드에서는 폭 자체가 유용한 관측량 | Group index는 slow-light 분석 때 표시 |
| Lambda AT | AT splitting, 중앙 투과율, Expected ≈ Ωc: 3개 | **AT splitting**, 필요 시 중앙 투과율 | Expected ≈ Ωc는 입력에서 환산한 비교값이므로 상세 |
| Rydberg EIT | EIT linewidth·공명 투과율 뒤에 Max spectral slope·IF discriminator·IF optimum detuning: 5개 | **EIT linewidth·공명 투과율 2개** | IF 및 slope 진단은 전기장 감도 분석에 한정 |
| Rydberg AT | Total/PSN sensitivity 2개와 RF parameter 환산 전기장·AT splitting·근사 전기장·center shift 등 총 6개. 기본 noise 조건에서 Total과 PSN 값 동일 | **RF AT splitting·Total sensitivity 2개** | PSN 분해, 두 전기장 환산, 중심 이동, 검출기 잡음. 비공명 설정에서만 중심 이동 보조 표시 |
| Magneto 자기장 스캔 | bias의 신호·PSN 감도, 조건부 Small-field SNR. 상태는 이미 접힘 | **선택 신호·PSN 감도 2개** | SNR은 선언된 입력 신호의 관측 가능성을 볼 때 조건부 표시. 감도 불가 상태는 숨기지 않음 |
| Magneto 광주파수 스캔 | 최적 감도·최적 detuning, SNR·centroid 대비 SNR 증가율 | **최적 감도·최적 detuning 2개** | SNR 비교율 및 모델 세부. 광주파수 최적값과 현재 bias 신호를 혼용하지 않음 |
| FWM intensity | Squeezing·probe gain·conjugate gain: 기본 3개는 각각 다른 결과. 보정 시 원래 atomic gain 2개 추가 | **Squeezing·probe/conjugate gain 3개 유지** | 보정 전 atomic gain은 비교 상세. 보정/외삽·noise 모델 조건은 결과와 연결 |
| FWM quadrature | 최저/최고 noise, 두 gain, 고정 detuning: 5개 | **최저/최고 quadrature noise 2개** | 위상에 따라 변하지 않는 gain·detuning |
| Biphoton | g²·Pair rate 뒤에 Model status·CAR·intrinsic/detected width: 6칸. Model status에 kind=status가 없어 일반 수치 칸을 점유 | **g²·Pair rate·CAR 3개**, 시간폭 분석 시 detected width 추가 | intrinsic width·bandwidth·phase matching. 모델 출처·수렴 상태는 수치 칸 대신 짧은 상태 표시 |

OD의 표기 주의: 현재 주값 `Peak absorbance A10 = −log10(T)`와 보조값 `τ = −ln(T)`는 다름. `Peak OD`로 줄일 때 수치 정의를 바꾸면 안 됨. 유한 probe에서는 effective collected-transmission 값이라는 범위도 보존.

SAS의 Δν는 현재 상단 metric에 없음. [Hyperfine lines 표](../../../../gabes/schemes/sas.py)는 전이 중심만 제공. 사용자 선택 두 선/두 커서의 간격을 별도로 구현해야 함. 원자종·전이·crossover 선택 없이 임의의 인접 피크를 자동 선택하면 사용자가 의도한 간격과 달라질 수 있음. 기존 FWHM 숫자를 Δν로 재명명할 수 없음.

Rydberg의 `RF field from fitted Ω_RF`는 현재 시뮬레이션 Rabi 파라미터의 전기장 환산. 실험 데이터를 실제로 fitting한 측정 결과로 표현하면 안 됨. Total과 PSN 감도는 기술 잡음이 생기면 달라지므로 값이 같다는 이유로 계산을 삭제하지 않고 PSN을 상세로 이동하는 권고.

## ⟦AI⟧ 이번 적용·검증

- [공통 이름 스타일](../../../../assets/ui/gabes.css): 0.75rem→0.9rem, 일반 두께→600, muted 회색→본문 ink 색. 말줄임 제거·줄바꿈 허용. 숫자·설명과의 계층 유지. 표시 항목 선정·계산 무변경.
- `python -m pytest -q tests/test_ui_readout.py tests/test_ui_design_feedback.py --durations=3`: **29 통과**, 1.46s.
- 전체 수치 회귀·GC 생략: CSS 변경과 표시 조사만 수행.
- 실제 FWM 화면 1440×1000·390×844: 이름 14.4px/600/`rgb(14,23,38)` 확인, 잘림0. [데스크톱](../../../../output/playwright/readout_audit_20261030/fwm_desktop.png) · [모바일](../../../../output/playwright/readout_audit_20261030/fwm_mobile.png).
- 동일 HTML/CSS로 23개 표시줄·102개 이름을 렌더한 모음도 두 폭에서 이름 잘림0·가로 overflow0. 실제 앱 23개 설정의 개별 브라우저 재실행을 뜻하지 않음.
