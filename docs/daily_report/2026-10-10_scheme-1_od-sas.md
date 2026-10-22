# GABES 일일 Scheme 물리 검토 — 2026-10-10

## 선택과 한 줄 판정

- 현지 날짜: 2026-10-10, Asia/Seoul.
- `n = (10 mod 5) + 1 = 1`.
- 대상: **Scheme 1 — Absorption spectroscopy (OD / SAS)**, registry key `sas`.
- 기준: `main` HEAD `39be4ae` + 검토 시작 시점의 기존 dirty working tree. 기존 변경은 보존했고 production 코드는 수정하지 않았다.
- 직전 Scheme 1 검토(2026-10-05)의 HEAD `fcb312d` 뒤 `39be4ae`가 Rb 조성 조절과 CSV affine 정렬을 추가했다. 핵심 OBE·선형흡수·SAS 응답식은 그대로다.
- **판정:** 실제 원자물리를 구현한다. 약한-probe OD, hyperfine 위치·세기, Doppler/Voigt 폭, Beer–Lambert 전파, 같은 속도군의 counter-propagating pump/probe 상관, hyperfine optical pumping은 실험 설계와 line assignment에 유용하다. 새 Rb isotope fraction과 trace alignment도 실험 workflow에 실질적 도움이 된다. 그러나 절대 pump Rabi, pressure shift, pump/probe 공통 collisional dephasing, polarization/mF, pump depletion, 검증된 measured holdout이 없으므로 **조건 탐색용 semi-quantitative reference**이지 절대 SAS 선폭·contrast·lock slope·cell parameter metrology 모델은 아니다.

## 현재 다섯 scheme의 정의와 순서

실제 dropdown 순서는 [`gabes/schemes/__init__.py:19`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/schemes/__init__.py:19>)의 `_SCHEMES` 순서다. `ODScheme`은 등록되지 않은 2-level 검증 primitive다.

| n | registry key | 현재 정의 | 코드 |
|---:|---|---|---|
| 1 | `sas` | Absorption spectroscopy: pump off OD, pump on SAS | [`sas.py:216`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/schemes/sas.py:216>) |
| 2 | `lambda` | Λ coherence: EIT / AT / CPT | [`absorption.py:466`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/schemes/absorption.py:466>) |
| 3 | `rydberg_eit` | ⁸⁵Rb cascade Rydberg-EIT / microwave AT electrometry | [`rydberg.py:94`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/schemes/rydberg.py:94>) |
| 4 | `magneto` | ⁸⁷Rb D1 magneto-optics: Hanle / EIA / MOR/NMOR | [`magneto.py:121`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/schemes/magneto.py:121>) |
| 5 | `fwm` | Four-wave mixing: mean-field Squeezing indicator / Biphoton | [`fwm.py:3626`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/schemes/fwm.py:3626>) |

## 선행 개선안·TODO·issue 검색

먼저 Scheme 1 이전 일일 보고서 14개, [`docs/checklist.json`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/docs/checklist.json>), README, User Guide, 코드의 `TODO/FIXME`, 로컬 issue/proposal/DEVLOG 파일명을 검색했다. 원격 공개 GitHub Issues API도 2026-10-10에 확인했고 issue는 0건이었다.

- **기존 개선안이 충분히 있다.** 새 physics roadmap이 없는 상태가 아니다.
- Scheme 1 전용 별도 local issue note는 없다. 실질적인 issue ledger는 checklist와 이전 일일 보고서다.
- 직접 코드 TODO는 Ne broadening scalar를 gas/species/line coefficient table과 pressure shift로 승격하라는 한 건이다: [`constants.py:79`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/constants.py:79>).
- 새로 드러난 ledger 문제: `od-sas-isotope-fraction-and-element-coverage`는 여전히 “fixed menu”와 `needs_decision`을 말하지만, Rb custom fraction은 `39be4ae`에서 구현·테스트됐다. Rb fraction 부분은 `done`으로 분리하고 Na/K 등 새 원소 결정만 별도 `needs_decision`으로 남기는 편이 정확하다: [`checklist:4179`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/docs/checklist.json:4179>), [`sas.py:94`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/schemes/sas.py:94>).

### 기존 제안의 계산비용과 물리 보존성

| 기존 항목 | 현재 상태 | 실행비용 | 물리 보존·실험 가치 |
|---|---|---:|---|
| collision coefficient provenance + pressure shift | P1 `ready` ([`checklist:2841`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/docs/checklist.json:2841>)) | coefficient lookup·detuning/width 보정 O(1), 같은 Liouvillian 차원·solve 수 | **negligible runtime.** zero-pressure anchor를 exact 보존하면서 잘못된 Cs fallback, pressure shift 누락, pump/probe width 불일치를 교정할 수 있다. 문헌 조사·검증 effort는 크다. |
| species/line/polarization별 pump Rabi·saturation | P1 `ready` ([`checklist:2881`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/docs/checklist.json:2881>)) | O(1) lookup/폐형식 | **negligible runtime.** pump-off OD는 exact 불변, pump-on power axis는 의도적으로 교정된다. 편광과 open-transition convention을 함께 고정해야 한다. |
| custom Rb isotope fraction | checklist는 P2 `needs_decision`, 코드는 구현 완료 | 내부 조성은 natural Rb와 같은 두 isotope solve 후 weight만 변경; 추가 차수 없음. 0/1 endpoint는 한 component를 건너뛰어 더 싸다 | total elemental density와 고정 Rb-85 centroid를 명시하고 legacy natural spectrum을 array-exact 보존한다. 새 원소 지원과 ledger 상태는 별도 과제다. |
| affine CSV alignment | 코드·README 구현 완료 | 고정 이론곡선 후처리. `min(Ncsv,1200)`점, 약 2,405개 coarse 후보 + 최대 5개 local basin; OBE solve 0회 | sweep 방향과 absolute transmission을 보존하는 유용한 display aid. **물리 파라미터 추정·주파수 교정·모델 검증은 아니다.** |
| physical spectrum fitting | P2 `ready` ([`checklist:4143`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/docs/checklist.json:4143>)) | 최소 `N_eval × compute`; pumped SAS는 매 평가마다 pump table solve | negligible 아님. uncertainty, identifiability, fit/holdout 분리, 평가 횟수 제한이 필요하다. 현재 affine alignment와 구분 유지. |
| dense scan×velocity 제거 | P3 `ready` ([`checklist:3531`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/docs/checklist.json:3531>)) | fused loop는 같은 O(`Nscan Nv`), pole reduction은 더 작음 | fused loop는 현 discrete measure·interpolation을 보존하면 roundoff 범위의 pure optimization. pole형은 검증된 Γ/8 table을 바꾸므로 독립 gate 필요. |
| weak-probe Zeeman/polarization reference | P2 research `parked` ([`checklist:4209`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/docs/checklist.json:4209>)) | hyperfine-Zeeman diagonalization + mode/Jones propagation | 실험 가치는 크지만 무부하가 아니다. 측정 geometry·polarization·holdout이 먼저다. |
| segmented cell / pump depletion | P2 `ready` ([`checklist:4243`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/docs/checklist.json:4243>)) | 대략 `Nz × atomic response` | high-OD pumped SAS에는 필요하지만 수렴 parameter가 필요한 별도 reference mode다. |
| coated-cell transport / full VCC | P2 research `parked` | dark-time·geometry·velocity-coupling 상태 추가, 현 separability 상실 | cell별 데이터 없이는 식별 불가. population-memory MVP를 기본 이상으로 과장하지 않는 현재 경계가 맞다. |
| result export / provenance | P2 `in_progress` ([`checklist:4270`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/docs/checklist.json:4270>)) | O(`Nscan`) serialization, 추가 solve 0회 | physics exact. 새 alignment/export provenance는 좋아졌지만 named arrays·round trip·lazy extra-view 범위가 남았다. |

우선순위는 기존 P1 두 건이 맞다. 런타임을 거의 늘리지 않으면서 현재 가장 위험한 “정밀해 보이는 잘못된 power/pressure 축”을 직접 고친다.

## 이번 변경: 실험적 가치와 신뢰 경계

### 1. 사용자 정의 Rb isotope fraction

- UI는 natural Rb에 대해 `Natural`/`Custom`을 선택하고 ⁸⁵Rb fraction을 0–1로 받는다: [`sas.py:249`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/schemes/sas.py:249>).
- `_isotope_composition`은 총 elemental density를 complementary ⁸⁵Rb/⁸⁷Rb share로 나누고 zero-share component를 생략한다. mixture의 frequency origin은 dominant isotope에 따라 튀지 않도록 항상 Rb-85 centroid다: [`sas.py:94`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/schemes/sas.py:94>), [`sas.py:423`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/schemes/sas.py:423>).
- raw output은 fraction, per-isotope share, frequency reference, reference-matched status를 남긴다: [`sas.py:565`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/schemes/sas.py:565>).
- legacy natural result array exact, endpoint/pure-isotope equivalence(선언된 centroid offset 포함), total density 보존, dominant-isotope crossing의 축 안정성, finite-probe passivity를 테스트한다: [`test_od_isotope_ratio.py:21`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/tests/test_od_isotope_ratio.py:21>), [`test_od_isotope_ratio.py:65`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/tests/test_od_isotope_ratio.py:65>), [`test_od_isotope_ratio.py:77`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/tests/test_od_isotope_ratio.py:77>), [`test_od_isotope_ratio.py:109`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/tests/test_od_isotope_ratio.py:109>).

**연구자 판정:** enriched/depleted Rb cell의 line-strength 비율과 isotope overlap을 보는 데 실제로 유용하다. 그러나 fraction을 데이터에서 추정하거나 uncertainty를 주는 fit은 아니며, isotope-dependent vapor pressure/activity나 cell 제조 편차까지 측정하는 모델도 아니다. 수동 입력을 알고 있는 cell composition으로 취급해야 한다.

### 2. CSV `Fit alignment`

- 모듈 자체가 “without fitting atomic parameters”라고 선언한다: [`experimental_alignment.py:1`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/experimental_alignment.py:1>).
- 최대 1,200점으로 decimate하고, 축 width/center를 bounded search하며 relative trace에만 robust positive affine Y fit을 허용한다. absolute calibration은 변형하지 않는다: [`experimental_alignment.py:199`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/experimental_alignment.py:199>), [`experimental_alignment.py:229`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/experimental_alignment.py:229>), [`experimental_alignment.py:242`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/experimental_alignment.py:242>).
- false positive를 막기 위해 overlap, normalized RMSE, correlation gate를 함께 쓴다: [`experimental_alignment.py:361`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/experimental_alignment.py:361>).
- random/ramp/sine mismatch 거부, absolute transmission 보존, outlier 내성, UI export와 plotted trace의 array equality가 테스트된다: [`test_experimental_alignment.py:109`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/tests/test_experimental_alignment.py:109>), [`test_experimental_alignment.py:175`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/tests/test_experimental_alignment.py:175>), [`test_experimental_alignment.py:189`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/tests/test_experimental_alignment.py:189>), [`test_od_comparison_ui.py:137`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/tests/test_od_comparison_ui.py:137>).

오늘 synthetic 1,201-point multipeak trace의 median은 **0.954 s**였고 `scale=399.9983`, `shift=0.699998`, RMSE `1.75×10⁻⁴`를 복구했다. 고정 이론곡선 뒤의 UI 후처리라 OBE overhead는 0이다. 이 결과가 실제 oscilloscope time axis의 Hz calibration을 만들어 주는 것은 아니다. etalon/known sideband/scan calibration 없이 fitted X scale을 원자 주파수 교정값으로 인용하면 안 된다.

## 구현 물리 감사

### 실제로 구현한 유용한 물리

1. **원자·hyperfine 구조.** ⁸⁵Rb/⁸⁷Rb/¹³³Cs D1/D2의 hyperfine A/B, line centre, mass, natural linewidth, allowed F→F′와 Wigner-6j/3j 세기, CG-branched decay를 갖는다: [`species.py:130`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/species.py:130>), [`species.py:337`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/species.py:337>).
2. **약한-probe 절대 면적.** `C_F²`가 magnetic sublevels를 이미 합하므로 `1/[2(2I+1)]`을 한 번만 적용한다. transition 면적은 `A_t = πk C_F²|d|²N/(ε₀ℏ)/[2(2I+1)]`: [`species.py:248`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/species.py:248>).
3. **SAS kinematics.** pump detuning `Δ+kv`, counter-propagating probe detuning `Δ−kv`를 같은 velocity class에서 곱해 평균한다. factorized average가 아니다: [`sas.py:600`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/schemes/sas.py:600>), [`sas.py:651`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/schemes/sas.py:651>), [`test_sas.py:249`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/tests/test_sas.py:249>).
4. **Hyperfine optical pumping.** pump steady state, CG-branched spontaneous decay와 thermal transit reload가 population transfer와 enhanced/inverted crossover를 낸다: [`species.py:365`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/species.py:365>), [`sas.py:624`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/schemes/sas.py:624>), [`test_sas.py:339`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/tests/test_sas.py:339>).
5. **흡수·분산의 한 complex profile.** 같은 population difference와 Doppler quadrature로 Im/Re χ를 만든다. 독립 Voigt와 OBE sign fixture가 있다: [`sas.py:600`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/schemes/sas.py:600>), [`test_sas.py:594`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/tests/test_sas.py:594>), [`test_sas.py:607`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/tests/test_sas.py:607>).
6. **유한 probe.** transition별 two-level saturation/power broadening, longitudinal depletion과 Gaussian radial average가 있다: [`probe_saturation.py:46`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/probe_saturation.py:46>), [`probe_saturation.py:95`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/probe_saturation.py:95>), [`probe_saturation.py:166`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/probe_saturation.py:166>). Rb D2는 software reference scope이며 D1/Cs는 generalized qualitative estimate라고 표시한다.
7. **실험 workflow.** density/Doppler temperature 분리, opt-in paraffin population memory, absolute/relative CSV calibration과 bidirectional sweep 보존, custom isotope ratio, fixed-theory alignment/export가 있다.

### 실험 reference로서의 핵심 한계

| 질문 | 현재 사용 가능성 |
|---|---|
| isotope/hyperfine line assignment, scan span, 온도·길이·조성 경향 | 좋음 |
| pump-off weak-probe Rb OD, Doppler/Voigt scale | 유용. 원자상수와 software-reference closure 중심 |
| Rb D2 finite-probe trend와 Gaussian depletion | 계산 reference로 유용. 독립 실험 validation 아님 |
| pumped SAS contrast/FWHM/lock slope, 절대 pump-power calibration | 정성/반정량 |
| high-pressure buffer cell, pressure shift, coated-cell absolute linewidth | 부적합 |
| Cs absolute saturation/self-broadening | 부적합. self-broadening은 Rb scalar를 상속 |
| polarization, mF optical pumping, magnetic field, Faraday/dichroism | 모델 밖 |
| 데이터로부터 temperature/cell length/isotope fraction 추정 | 현재 alignment 기능 밖; uncertainty-aware physical fit 없음 |

특히 probe homogeneous width는 `natural + self + buffer`를 사용하지만 pump Liouvillian에는 같은 elastic optical dephasing이 들어가지 않는다: [`sas.py:442`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/schemes/sas.py:442>), [`sas.py:456`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/schemes/sas.py:456>), [`sas.py:629`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/schemes/sas.py:629>). 또한 `self_broadened_gamma(iso, N)`은 `iso`를 받지만 실제로는 단일 Rb `BETA_SELF`만 쓴다: [`species.py:210`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/species.py:210>). pressure가 올라갈수록 precise-looking SAS hole/crossover를 정량 reference로 쓰기 어려운 이유다.

Pump Rabi도 모든 종·선에 고정 `I_sat=4.484 mW/cm²`를 사용한다: [`constants.py:30`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/constants.py:30>), [`species.py:272`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/species.py:272>). line별 closed-two-level finite-probe helper와 비교한 현재/helper Ω 비는 Rb D1/D2, Cs D1/D2에서 각각 약 **0.5776 / 0.6102 / 0.4309 / 0.4958**이다. helper식도 open hyperfine·polarization-resolved pump의 정답은 아니며, 이 차이는 P1 provenance가 필요한 정량화다.

## 오늘 재현한 수치

기본값: natural Rb D1, 40 °C, 75 mm, pump 0.5 mW, probe 0 µW, pump/probe waist 1 mm, transit `γ_t/2π=100 kHz`, 1,401 points. warm compute 3회 median.

| 경로 | 결과 |
|---|---:|
| 기본 SAS compute | 0.7825 s median |
| SAS `Tmin` | 0.1180495075 |
| SAS peak `αL` | 2.136651188 |
| 표시 sub-Doppler FWHM | 18.86836 MHz |
| samples/FWHM | 2.3418 |
| resolution | `resolution-limited` |
| pump-off OD compute | 0.3707 s median |
| pump-off `Tmin` | 0.1050956546 |
| pump-off peak `αL` | 2.252884348 |
| custom fraction ⁸⁵Rb=0.35 compute | 0.7025 s median, fractions `0.35/0.65` |

기본 수치는 2026-10-05와 동일하다. 표시 FWHM은 2.34 samples밖에 걸치지 않으므로 실험 hero 값으로 인용하면 안 된다. 코드도 six-sample gate 전에 lock slope를 숨기고, 4×/8× fixed-window reference에서만 수렴을 검사한다: [`test_sas.py:406`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/tests/test_sas.py:406>), [`test_sas.py:442`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/tests/test_sas.py:442>).

## 동작을 바꾸지 않는 순수 코드 최적화

오늘 기본 SAS cProfile 1회는 0.838 s였다.

- `_component_response`: 0.803 s cumulative / 0.323 s own.
- `np.interp`: 8회, 0.237 s.
- `_velocity_correlated_average`: 16회, 0.150 s.
- `_pump_pops`: 0.090 s.

즉 선형 solve보다 dense scan×velocity interpolation·temporary·matvec가 우세하다: [`sas.py:600`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/schemes/sas.py:600>).

1. **Pump-off scalar fast path — 즉시 후보.** `_velocity_correlated_average(1.0, profile, weights)`는 불필요한 full dense multiply/copy 뒤 matvec를 한다: [`sas.py:1070`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/schemes/sas.py:1070>). synthetic `1401×1403`, 9회 median에서 **6.4255 → 0.8587 ms, 7.48×**, `np.array_equal=True`, max difference 0. `np.isscalar(prepared_population)`이면 `probe_profile @ weights`로 분기하면 물리·배열 exact다.
2. **Pump-off finite-probe 중복 matvec 제거.** `transition_alpha`를 만든 직후 pump-off에서 같은 `Lp @ wt`를 `thermal_alpha`로 다시 계산한다: [`sas.py:690`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/schemes/sas.py:690>), [`sas.py:694`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/schemes/sas.py:694>). `Op<=0`이면 `thermal_alpha = transition_alpha` 재사용이 exact하다.
3. **기존 P3 fused compiled loop.** pump-on의 `deff_grid`, `pop_at`, `arg`, `denom`, `Lp`, `Dp` 대형 temporary와 `np.interp`를 velocity-order 고정 compiled loop로 합친다. 현 1 m/s·4σ discrete measure, Γ/8 pump table, interpolation 규칙, transition sum order를 regression gate로 고정하면 behavior를 허용 roundoff 내 보존할 수 있다. 3× 목표는 profile상 여전히 현실적인 후보다.
4. **Legendre node cache.** finite-probe observable마다 `leggauss(32)`를 다시 만든다: [`probe_saturation.py:195`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/probe_saturation.py:195>). order-key immutable cache는 수식·node를 바꾸지 않는 작은 최적화다.
5. **Finite-probe frequency inversion loop compile.** 각 frequency에서 `np.interp` 두 번을 호출하는 loop를 같은 monotone table과 연산 순서로 compile한다: [`probe_saturation.py:236`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/probe_saturation.py:236>). finite-probe batch에서만 우선순위가 높다.

`build_manifold`는 이미 bounded LRU cache다: [`species.py:336`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/species.py:336>). 여기의 추가 cache는 우선순위가 낮다. affine alignment는 약 0.95 s지만 coarse grid·robust gate를 줄이면 acceptance behavior가 바뀔 수 있으므로, profile 없이 pure optimization이라고 부르지 않는다.

## 문서와 예제 감사

### 잘 된 부분

- README는 새 Fit alignment를 affine display match로 제한하고 physical parameter estimation이 아니라고 명시한다: [`README.md:121`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/README.md:121>).
- custom isotope fraction, fixed Rb-85 centroid, cached theory reuse와 export provenance를 설명한다: [`README.md:129`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/README.md:129>), [`README.md:135`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/README.md:135>).
- 새 unit/UI tests가 실질적인 실행 예제 역할을 하며, 특히 seconds-axis scope trace, fixed theory, export 일치를 재현한다: [`test_od_comparison_ui.py:28`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/tests/test_od_comparison_ui.py:28>).

### 남은 문서·예제 문제

- User Guide는 여전히 ⁸⁵Rb D1을 “AutoOD와 0.1% 이내”라고 표현해 독립 experimental validation처럼 읽힌다: [`GABES_User_Guide_v2.html:409`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/docs/Userguide/GABES_User_Guide_v2.html:409>), [`:922`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/docs/Userguide/GABES_User_Guide_v2.html:922>). 실제 의미는 normalization-corrected software-reference closure다.
- Guide는 finite-probe knob가 natural Rb D2에서만 나타난다고 쓰고 기본을 `0 mW`로 적는다: [`Guide:556`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/docs/Userguide/GABES_User_Guide_v2.html:556>), [`Guide:573`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/docs/Userguide/GABES_User_Guide_v2.html:573>). 현재 UI는 built-in D1/D2와 Rb/Cs에 0–100 µW를 노출하고 scope qualifier를 붙인다: [`sas.py:236`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/schemes/sas.py:236>).
- Guide에는 custom isotope fraction과 Fit alignment workflow가 없다. 이것은 새 기능의 가장 직접적인 사용 예이므로 업데이트 우선순위가 높다.
- Guide의 `od.png`, `sas.png`는 2026-06-08 생성본이다. 9월 normalization correction, 현재 isotope controls와 readout을 반영하지 않는다.
- standalone `examples/` 디렉터리와 재현 가능한 Scheme 1 config/script는 없다. README와 tests가 de facto examples다. `SAS.defaults()` → OD/SAS preset → `headless_observables` → CSV alignment/export를 잇는 짧은 example은 추가 solve를 만들지 않는 문서성 개선이다.

## 테스트

### Scheme 1 집중 회귀

비-UI 물리·CSV·alignment·isotope-ratio 경로:

```text
python -m pytest -q tests/test_sas.py tests/test_absorption.py \
  tests/test_probe_saturation.py tests/test_experimental_csv.py \
  tests/test_experimental_alignment.py tests/test_od_isotope_ratio.py \
  tests/test_lineshape.py --tb=short
186 passed, 1 skipped in 90.45s
```

Streamlit comparison UI:

```text
python -m pytest -q tests/test_od_comparison_ui.py --tb=short
8 passed in 20.58s
```

처음 두 묶음을 한 프로세스로 실행했을 때 약 73% 뒤 pytest teardown에서 Windows native `0x80000003`으로 프로세스가 summary 없이 종료됐다. 두 묶음을 분리하면 총 **194 passed, 1 skipped**로 통과했으므로 Scheme 1 계산 assertion 실패로 재현되지는 않았다. 환경/teardown 결함은 실패한 결합 실행과 분리 통과를 함께 기록한다.

### 저장소 필수 전체 회귀

```text
python -m pytest -q
2 failed, 2202 passed, 4 skipped in 2844.96s (0:47:24)
```

두 실패는 직전 2026-10-08/09 전체 회귀와 동일하며 Scheme 1 집중 테스트 실패가 아니다.

1. [`tests/quantum/test_thermal_campaign_extend.py:367`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/tests/quantum/test_thermal_campaign_extend.py:367>): Windows Grand Challenge fresh capsule이 중첩 `.git/grand-challenge-runs/...` 아래 `ballistic_linear_channel.py`를 만들 때 destination parent가 없어 `FileNotFoundError`.
2. [`tests/test_fwm_excess_noise.py:541`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/tests/test_fwm_excess_noise.py:541>): Streamlit ScrubField fallback에서 gain-proportional slope `0.02`의 기대 metric `0.04` 대신 이전 constant noise `2.0`이 유지됨.

Scheme 1 관련 집중 회귀는 위의 분리 실행에서 전부 통과했다.

## 최종 우선순위

1. **P1, runtime 거의 0:** `(gas, species, line)` 충돌계수·pressure shift·elastic optical dephasing을 sourced table로 만들고 pump/probe가 한 homogeneous-width convention을 쓰게 한다.
2. **P1, O(1):** species/line/polarization별 pump Rabi/saturation 정의를 독립 근거와 함께 고정한다. pump-off OD anchor는 exact 보존한다.
3. **P1 trust, O(1):** checklist의 isotope item을 “Rb fraction done / 새 원소 needs decision”으로 분리하고 User Guide의 scope·단위·새 workflow·이미지를 현재 코드와 동기화한다.
4. **P2 사용 규칙:** 기본 SAS FWHM·lock slope는 `resolution-limited`이면 인용하지 않는다. fixed-window resolved refinement 또는 실험 scan calibration을 사용한다.
5. **성능:** exact scalar fast path와 pump-off 중복 matvec부터 적용하고, 그 다음 discrete-measure preserving fused loop를 benchmark/regression gate와 함께 검토한다.
6. **비-negligible physics:** uncertainty-aware physical fitting, segmented pump depletion, full VCC, Zeeman/polarization은 측정 target·holdout·수렴 budget을 갖춘 opt-in reference로 진행한다. 현재 affine alignment를 그 대체물로 취급하지 않는다.

## 결론

Scheme 1은 “그럴듯한 Lorentzian 그림”보다 훨씬 낫다. 원자상수, hyperfine branching, thermal velocity, same-class pump/probe 상관과 optical pumping이 연결된 실제 semi-quantitative SAS/OD 모델이고, 조성·온도·길이·line assignment를 탐색하는 실험 reference로 유용하다. 이번 custom Rb fraction과 trace alignment는 계산 차원을 늘리지 않으면서 실제 lab workflow를 개선했다. 다만 정량 SAS 계측의 가장 중요한 power/pressure/polarization 계약과 독립 measured validation은 아직 비어 있다. 따라서 현재 가장 과학적으로 경제적인 다음 단계는 새 복잡한 solver가 아니라 **기존 P1 Rabi/collision provenance를 같은 solve 차원에서 닫고, 문서·ledger를 현재 구현과 맞춘 뒤, exact dense arithmetic 최적화를 이식하는 것**이다.
