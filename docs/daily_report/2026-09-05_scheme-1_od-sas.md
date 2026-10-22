# GABES 일일 Scheme 물리 검토 — 2026-09-05

## 1. 오늘의 선택과 검토 기준

- 현지 날짜(Asia/Seoul)의 일(day)은 5이므로
  \(n=(5\bmod5)+1=1\)이다.
- 따라서 오늘 검토한 대상은 **Scheme 1 — Absorption spectroscopy (OD / SAS)**
  (`sas`)이다.
- 검토 기준 커밋은 `d98bdf4` (`perf(core): streamline solves and readouts`)이다.
  작업 트리에는 이 검토 전부터 다른 영역의 변경이 있었으며, 본 작업은 코드나 기존
  변경을 수정하지 않고 이 보고서만 추가했다.

현재 드롭다운의 다섯 scheme은
[`gabes/schemes/__init__.py:19`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/schemes/__init__.py:19>)의
등록 순서로 결정된다.

| 순번 | registry key | 현재 표시명 | 정의 위치 |
|---:|---|---|---|
| 1 | `sas` | Absorption spectroscopy (OD / SAS) | [`sas.py:175`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/schemes/sas.py:175>) |
| 2 | `lambda` | Λ coherence (EIT / AT / CPT) | [`absorption.py:466`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/schemes/absorption.py:466>) |
| 3 | `rydberg_eit` | Rydberg-EIT electrometry | [`rydberg.py:95`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/schemes/rydberg.py:95>) |
| 4 | `magneto` | Magneto-optics (Hanle/MOR) | [`magneto.py:131`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/schemes/magneto.py:131>) |
| 5 | `fwm` | Four-wave mixing (Squeezing / Biphoton) | [`fwm.py:2987`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/schemes/fwm.py:2987>) |

`ODScheme`는 단일 2-level 검증 primitive로 남아 있지만 드롭다운에는 없고, OD는
Scheme 1에서 pump power를 0으로 둔 극한이다
([`gabes/schemes/__init__.py:12`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/schemes/__init__.py:12>)).

## 2. 먼저 찾은 기존 제안과 이슈

새 제안을 만들기 전에 다음 자료를 확인했다.

- Scheme 1 선행 보고서 8건: 2026-06-25, 07-10, 07-20, 07-30, 08-05,
  08-10, 08-25, 08-30. 최신 선행 검토는
  [`2026-08-30_scheme-1_od-sas.md`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/docs/daily_report/2026-08-30_scheme-1_od-sas.md:1>)이다.
- 공식 개선 항목은
  [`docs/checklist.json`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/docs/checklist.json:219>)에
  이미 충분히 기록되어 있다. 따라서 “기존 제안 없음”에 해당하지 않는다.
- 코드 내 Scheme 1 관련 TODO는 Ne를 gas/species/line별 계수표와 압력 이동으로
  승격하자는 한 건이다
  ([`gabes/constants.py:79`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/constants.py:79>)).
- 원격 저장소의 GitHub issue API를 2026-09-05에 확인한 결과 open 0건,
  closed 0건이었다. 별도 issue note는 없었다.

### 기존 제안의 계산비용과 기본 물리 보존성

| 기존 항목 | 상태/규모 | 계산비용 판단 | 물리를 거의 무비용으로 보존할 수 있는가 |
|---|---|---|---|
| 충돌계수 provenance + pressure shift | P1, ready, large | 계수 조회와 detuning/linewidth 보정 자체는 transition당 O(1)로 사실상 무시 가능. 큰 비용은 데이터 수집·검증이다. | **예.** Rb-85 D1 기준 경로를 고정하면서 지원 조합만 올바른 계수로 바꿀 수 있다. 알려지지 않은 조합은 타 종 계수를 재사용하지 말고 `unsupported/qualitative`로 표시해야 한다. 기준은 [`checklist:219`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/docs/checklist.json:219>)에 있다. |
| species/line별 pump Rabi·saturation | P1, ready, medium | lookup/폐형식 계산뿐이라 runtime 증가는 무시 가능하다. | **예.** pump-off OD는 그대로이고 pumped 결과만 물리적으로 교정된다. 현 pump 경로의 단일 `I_SAT` 재사용은 [`species.py:272`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/species.py:272>)에 남아 있다. |
| 실험 spectrum fitting | P2, ready, large | 약 \(O(N_\mathrm{eval}N_\mathrm{scan})\). pump-off도 반복 계산이고 pumped fit은 매 평가마다 OBE table을 다시 풀어 비싸다. | **기본 경로에는 아니오.** 명시적 실행, 평가 횟수 제한, cache, fit/holdout 분리가 필요하다. [`checklist:1245`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/docs/checklist.json:1245>)의 제한이 타당하다. |
| isotope fraction | P2, needs decision, medium | natural Rb가 이미 두 isotope를 합산하므로 비율만 바꾸는 비용은 거의 없다. | **예.** 기존 natural/pure preset을 exact default로 유지할 수 있다. 새 원소 추가는 별도 원자·증기압·충돌 데이터가 필요하다. [`checklist:1281`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/docs/checklist.json:1281>) |
| segmented cell / pump depletion | P2, ready, large | 대략 segment 수에 선형 비례하므로 무시할 수 없다. | **조건부.** 1 segment를 현재 경로와 exact하게 하고 다중 segment는 opt-in reference mode로 두어야 한다. [`checklist:1345`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/docs/checklist.json:1345>) |
| result export + provenance | P2, ready, medium | 계산 뒤 O(`scan_points`) 직렬화뿐이라 OBE 대비 작다. | **예.** 물리 배열을 재계산하지 말고 현재 named arrays를 그대로 내보내면 된다. [`checklist:1372`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/docs/checklist.json:1372>) |
| coated-cell 확장, full VCC, weak-probe magneto-optics | P2, research/parked | geometry/velocity/mF 차원이 늘어 계산량과 검증면이 모두 커진다. | **아니오.** 데이터·측정 target·runtime budget이 생기기 전에는 기본 모델에 넣지 않는 현재 판단이 맞다. [`coating:315`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/docs/checklist.json:315>), [`VCC:614`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/docs/checklist.json:614>), [`magneto-optics:1311`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/docs/checklist.json:1311>) |

Sub-Doppler 신뢰도 판정, opt-in paraffin population memory, CSV 절대보정과 sweep
보존, 분산 quadrature, density/Doppler 온도 분리는 이미 `done`이다
([`checklist:278`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/docs/checklist.json:278>),
[`299`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/docs/checklist.json:299>),
[`339`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/docs/checklist.json:339>),
[`1189`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/docs/checklist.json:1189>),
[`1217`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/docs/checklist.json:1217>)).

## 3. 실험 원자물리 관점의 판정

### 결론

**Scheme 1은 실제 물리를 구현한 유용한 실험 계획·해석 보조 도구다.** 특히
약한 probe의 Rb OD, hyperfine line 위치와 넓은 Doppler envelope, pump-induced
Lamb dip/crossover의 존재와 상대적 경향, weak-probe dispersion의 부호와 규모를
빠르게 확인하는 데 가치가 있다.

다만 **실험 정밀 기준(reference model)이라기보다 계산 가능한 semi-quantitative
reference**로 보는 것이 맞다. 절대 OD 기준은 ⁸⁵Rb D1 pump-off에 강하고,
finite probe는 Rb D2 독립 프로그램에 대한 software regression일 뿐 실험 검증이
아니다. pumped finite-probe, Cs, buffer gas, coated cell의 정량 추론은 현재 표시되는
qualifier를 반드시 결과와 함께 남겨야 한다.

### 실제로 잘 구현된 물리

1. **원자·hyperfine 구조와 절대 OD**

   - Rb natural/⁸⁵Rb/⁸⁷Rb/¹³³Cs와 D1/D2 manifold를 구성하고, Wigner 계수에서
     transition strength와 절대 integrated absorption을 만든다
     ([`species.py:225`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/species.py:225>),
     [`species.py:248`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/species.py:248>)).
   - spontaneous decay branching, thermal ground population, transit reload를 같은
     line-strength 구조에서 구성한다
     ([`species.py:337`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/species.py:337>)).
   - pump-off는 thermal population factor `w=1`을 정확히 사용하고, unit-area
     Lorentzian을 Maxwell velocity class와 합쳐 Voigt형 OD를 만든다
     ([`sas.py:560`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/schemes/sas.py:560>),
     [`sas.py:619`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/schemes/sas.py:619>)).

2. **SAS와 optical pumping**

   - pump detuning은 `Δ+kv`, 역행 probe는 `Δ-kv`로 두고 동일한 velocity class의
     pump-prepared population 차이를 probe와 곱한다
     ([`sas.py:590`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/schemes/sas.py:590>)).
     이는 velocity average를 서로 독립 평균으로 잘못 factorize하지 않는 중요한 구현이다.
   - CG-branched decay와 transit relaxation이 ground-hyperfine optical pumping을
     만들기 때문에 단순 Lorentzian hole보다 실제 alkali SAS의 crossover 증강·반전을
     더 잘 반영한다.

3. **Weak-probe dispersion과 두 온도**

   - absorption과 Re χ는 동일한 complex line profile과 같은 population/velocity
     weight를 사용한다
     ([`sas.py:539`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/schemes/sas.py:539>)).
     따라서 pump-on dispersion도 별도 곡선을 붙인 것이 아니라 같은 준비 상태의
     Kramers–Kronig partner다.
   - density temperature와 Doppler-width temperature를 기본적으로 묶되, cold spot과
     illuminated column이 다른 셀을 위해 풀 수 있다
     ([`sas.py:214`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/schemes/sas.py:214>)).

4. **Finite probe와 공간 전파**

   - Gaussian `I0=2P/(πw²)`, line별 `I_sat`, transition별 saturation, power broadening,
     `dI/dz=-α(I)I`, 32-point Gaussian power average를 구현한다
     ([`probe_saturation.py:34`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/probe_saturation.py:34>),
     [`probe_saturation.py:95`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/probe_saturation.py:95>),
     [`probe_saturation.py:166`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/probe_saturation.py:166>)).
   - Rb D2 reference-matched, Rb D2 extrapolation, D1/Cs generalized estimate를
     구분하고 pump+probe의 frozen-population 한계를 출력에 드러낸다
     ([`sas.py:92`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/schemes/sas.py:92>),
     [`sas.py:113`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/schemes/sas.py:113>)).
   - 분리 근사가 음의 흡수를 만들면 residual만 축소해 passivity를 강제하고 그 사실을
     metric으로 보고한다
     ([`sas.py:466`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/schemes/sas.py:466>),
     [`sas.py:1188`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/schemes/sas.py:1188>)).

5. **실험 데이터 취급**

   - CSV는 relative extrema와 증거가 있는 absolute dark/reference 또는 gain/offset을
     구분하고, acquisition-order sweep branch를 보존한다
     ([`experimental_csv.py:214`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/experimental_csv.py:214>),
     [`experimental_csv.py:443`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/experimental_csv.py:443>)).
     단순히 측정 trace의 최솟값·최댓값을 0·1로 바꾼 뒤 absolute transmission이라
     부르는 오류를 피한 점은 실험 reference로서 좋다.

### 정량 reference로 사용할 때의 핵심 한계

1. **Pump Rabi scale이 아직 종·선별이 아니다.** `pump_rabi_from_power()`는 모든
   isotope/D1/D2에 고정 `I_SAT=4.484 mW/cm²`를 쓴다
   ([`constants.py:30`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/constants.py:30>),
   [`species.py:272`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/species.py:272>)).
   0.5 mW, 1 mm waist에서 확인한 Ω/Γ는 ⁸⁵Rb D1/D2, ⁸⁷Rb D1/D2,
   ¹³³Cs D1/D2 모두 `1.883984233`으로 동일했다. Pumped spectrum의 power axis를
   정밀 calibration으로 쓰면 안 된다.

2. **Self broadening이 Cs까지 Rb 계수를 재사용한다.** 함수가 `iso`를 받지만 실제로는
   단일 `BETA_SELF`만 곱한다
   ([`species.py:49`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/species.py:49>),
   [`species.py:210`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/species.py:210>)).
   Cs에 그 값을 적용하면 Γself/2π가 100/150/200 °C에서 각각
   `1.0186/15.1751/124.0098 MHz`까지 커진다. 이는 작은 correction이 아니라 고온
   spectrum을 지배할 수 있는 잘못된 fallback이다.

3. **Ne는 broadening-only scalar다.** pressure shift, 온도/선 의존성, Dicke narrowing,
   speed-dependent collision과 VCC가 없다. 현재 UI도 이를 명시한다
   ([`sas.py:247`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/schemes/sas.py:247>)).

4. **Pump는 셀 입구 세기로 한 번만 계산된다.** probe는 finite-power일 때 z 방향으로
   전파하지만 pump depletion과 셀 gradient는 풀지 않는다. 고 OD에서 SAS contrast와
   hole depth는 낙관적일 수 있다.

5. **Hyperfine-lumped 모델이라 mF, polarization, B field가 없다.** Faraday/Voigt
   rotation, dichroism, optical-pumping orientation/alignment을 Scheme 1에서 추론하면
   안 된다. 이 문제는 별도 weak-probe magneto-optical reference mode로 이미 parked다.

6. **동시 finite pump + finite probe는 joint OBE가 아니다.** Pump population을 먼저
   고정한 뒤 probe broadening을 적용하므로 probe-induced hyperfine repumping과 coherent
   cross-line coupling이 빠진다. Passivity limiter는 비물리적 gain을 막는 안전장치이지
   누락된 동역학을 복원하는 모델은 아니다.

7. **Paraffin은 population-memory MVP다.** 고정 `T1=25.1 ms`와 기존 transit rate를
   재사용하며 Zeeman coherence, geometry별 dark-time distribution, wall shift를 포함하지
   않는다
   ([`sas.py:313`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/schemes/sas.py:313>)).

## 4. 오늘의 수치 spot check

같은 머신의 warm median 또는 직접 계산값이며 절대 runtime은 환경에 따라 달라진다.

| 확인 | 결과 | 해석 |
|---|---:|---|
| 기본 SAS compute | `0.807 s` | natural Rb D1, pump 0.5 mW, probe 0, 40 °C, 75 mm, 1401 points |
| 기본 sub-Doppler FWHM | `21.56 MHz` | `2.7 samples/FWHM`, 따라서 상태는 `resolution-limited` |
| 기본 peak OD / Gaussian Doppler FWHM | `0.54 / 518.7 MHz` | broad envelope 규모는 유용하지만 sub-Doppler width 정량 인용에는 local refinement가 필요 |
| 기본 weak-probe peak phase | `421.4 mrad` | phase이며 group delay가 아님 |
| ⁸⁵Rb D1 pump-off vs 기존 ODScheme | area `-0.0393169%`, peak `-0.0125734%` | 문서의 `<0.1%` absolute-scale 주장을 만족 |
| natural Rb D2, pump off, probe 100 µW | effective peak OD `0.23`, peak s₀ `3.81` | `Rb D2 reference-matched`; software reference 범위 |
| natural Rb D1, pump off, probe 6 µW | effective peak OD `0.55`, peak s₀ `0.26` | `Generalized estimate`; 정량 검증으로 승격 불가 |

기본 SAS에서 resolution 상태가 hero metric이 되고, resolved일 때만 lock slope를 내는
로직은 [`sas.py:812`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/schemes/sas.py:812>)와
그 회귀 테스트
[`tests/test_sas.py:382`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/tests/test_sas.py:382>)에 있다.

## 5. 문서·테스트·예제 평가

### 문서

- 현재 [`README.md:289`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/README.md:289>)는
  pump-off exact limit, Rb D2 reference tier, D1/Cs generalized tier, paraffin,
  분리 온도와 weak-probe dispersion을 비교적 정확하게 설명한다.
- 반면 HTML User Guide는 최신 코드와 어긋난다.
  [`GABES_User_Guide_v2.html:556`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/docs/Userguide/GABES_User_Guide_v2.html:556>)은
  finite probe를 natural Rb D2에서만 쓸 수 있다고 하고,
  [`line 572`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/docs/Userguide/GABES_User_Guide_v2.html:572>)는
  단위를 `mW`로 적는다. 현재 코드는 모든 built-in D line에 0–100 µW control을 보이고
  D1/Cs를 generalized로 표시한다
  ([`sas.py:195`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/schemes/sas.py:195>)).
- User Guide 표에는 Dispersion, phase, paraffin, density/Doppler 두 온도도 빠져 있다
  ([`GABES_User_Guide_v2.html:527`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/docs/Userguide/GABES_User_Guide_v2.html:527>)).
- OD/SAS 그림
  [`od.png`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/docs/Userguide/userguide_assets/od.png>)와
  [`sas.png`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/docs/Userguide/userguide_assets/sas.png>)는
  2026-06-08 생성본이다. 현재의 분리된 Transmission/Optical density/Dispersion view와
  provenance metric을 보여 주지 못한다.

### 테스트

- 제품 테스트: `python -m pytest -q tests` → **577 passed in 334.20 s**.
- Scheme 1 + CSV + 별도 Rb-D2 reference 테스트:
  **135 passed, 1 failed in 63.22 s**.
- 실패한 것은
  [`test_complete_spectra_match_preoptimization_golden_hashes`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/references/AutoOD-Rb/test_beam_saturation.py:180>)의
  raw array SHA-256 고정값이다. 같은 실행에서 weak/finite numerical anchors는 reference
  문서의 값과 일치했다: weak `Tmin=0.468823515134`, peak OD `0.757528881708`;
  finite `Tmin=1.17288769e-11`, peak OD `25.1689672059`
  ([`model.tex:648`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/references/AutoOD-Rb/docs/autood_natural_rb_d2_model.tex:648>)).
  따라서 이 실패는 현재 확인 범위에서는 physics mismatch 증거가 아니라 BLAS/FFT/플랫폼에
  취약한 byte-exact oracle 문제다. 값 기반 tolerance anchor를 주 oracle로 삼는 편이 낫다.
- repository root의 의무 실행 `python -m pytest -q`는 제품 테스트 전에
  `references/AutoOD/test_autood_gui_adjustments.py`를 수집하다 `sympy` 미설치로 중단됐다.
  `references/*`는 원래 gitignore 대상이고
  ([`.gitignore:19`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/.gitignore:19>)),
  기본 [`requirements.txt`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/requirements.txt:1>)에도
  `sympy`가 없다. reference test를 공식 수트에서 제외하거나 별도 optional test
  dependency/profile로 분리해야 재현 가능한 검증 계약이 된다.

### 예제와 데이터

- 추적되는 `examples/` 디렉터리나 Scheme 1 standalone example은 없다. 실제 사용 예는
  README, Streamlit UI, tests가 대신한다. 실험 연구자에게는 “pump-off OD → pump-on SAS →
  finite probe → CSV overlay”를 고정 조건으로 재현하는 작은 scripted example이 유용하다.
- 로컬의 독립 reference
  [`AutoOD-NatRbD2`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/references/AutoOD-Rb/README.txt:1>)는
  전이별 two-level saturation과 Gaussian propagation을 잘 문서화하지만, 스스로도
  mathematical/software verification이며 finite-power experimental validation은 아니라고
  명시한다
  ([`model.tex:681`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/references/AutoOD-Rb/docs/autood_natural_rb_d2_model.tex:681>)).
- 로컬 CSV
  [`OD refence_NatRbD2_50C6uW150um^2 .csv`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/references/OD refence_NatRbD2_50C6uW150um^2 .csv>)는
  7,674 valid rows, 1,350 unique x, 6,324 duplicate merges, 102 gap segments로 읽혔다.
  기본 보정은 `relative_extrema`이고 98개 sample이 보정 범위 밖이었다. 헤더, x 단위,
  dark/reference, cell length, 주파수 calibration, pump 상태가 없어 파일명만으로
  `50C/6uW/150um²` 조건을 확정할 수 없다. 따라서 현재는 parser stress/example data이지
  물리 validation dataset이 아니다. SHA-256은
  `9118e088366bac6fd6a6832f7d1d00e40084ea2d24ca98ee62cedcdb7c666de9`이다.

## 6. 동작을 바꾸지 않는 순수 코드 최적화

### 최신 커밋에서 이미 반영된 것

1. **Pump-off OBE solve 생략**

   Pump가 0이면 Hamiltonian/Liouvillian/`_pump_pops()`를 만들지 않고 `w=1`을 쓴다
   ([`sas.py:563`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/schemes/sas.py:563>),
   [`test_sas.py:58`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/tests/test_sas.py:58>)).
   parent `d98bdf4^`와 비교한 ⁸⁵Rb D1, 90 °C, 12.5 mm, 1401-point warm median은
   `0.35625 → 0.18391 s`, **1.94×**였다. scan은 exact, alpha 최대 상대차는
   `5.63e-16`이었다. 선행 보고서가 제안한 가장 큰 Scheme 1 compute 개선이 반영됐다.

2. **Zero residual FFT 생략**

   finite probe의 pump-off residual처럼 profile이 정확히 0이면 FFT table을 즉시 0으로
   반환한다
   ([`probe_saturation.py:123`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/probe_saturation.py:123>),
   [`test_probe_saturation.py:92`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/tests/test_probe_saturation.py:92>)).

3. **Duplicate-heavy CSV median vectorization**

   `(x,y)` lexicographic sort와 중앙 index로 duplicate median을 벡터화했다
   ([`_csv_numeric.py:6`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/_csv_numeric.py:6>),
   [`experimental_csv.py:508`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/experimental_csv.py:508>)).
   500,000 rows/각 x 5회 반복 fixture에서 기존 group loop `2.15537 s`, 현재
   `0.18347 s`, **11.75×**였고 x/y 결과는 exact equal이었다.

### 오늘 확인한 남은 안전 후보

1. **Pump-off finite-probe에서 thermal profile 중복 matvec 제거 — P3**

   `transition_alpha`는 pump-off에서 이미 `A_t*(Lp@wt)`인데 바로 아래에서
   `thermal_alpha`를 같은 식으로 다시 계산한다
   ([`sas.py:629`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/schemes/sas.py:629>)).
   pump-off일 때 `thermal_alpha=transition_alpha`를 재사용한 in-memory 후보는 scan,
   alpha, Re χ, intensity grid, finite-probe table이 모두 bitwise exact였고 natural Rb D2,
   100 µW, 1401 points에서 `0.75527 → 0.74578 s`, **1.013×**였다.
   안전하지만 이득이 작으므로 다른 코드 작업과 함께 처리할 항목이다.

2. **Gauss–Legendre node/weight 캐시 — P3**

   light observables가 호출될 때마다 `leggauss(32)`를 다시 만든다
   ([`probe_saturation.py:195`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/probe_saturation.py:195>)).
   immutable `(nodes, weights)` 캐시 후보는 transmission/effective alpha가 bitwise exact였고
   propagation median이 `0.02685 → 0.02588 s`, **1.038×**였다. Cell length나
   line-strength slider처럼 heavy solve를 재사용하는 조작에서만 보이는 작은 개선이다.

3. **Coated basis-conditioned multi-RHS — 보류**

   선행 보고서는 해당 stage `0.12159 → 0.09901 s`(1.23×), population 최대차
   `3.65e-12`를 측정했다
   ([`2026-08-30 report:327`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/docs/daily_report/2026-08-30_scheme-1_od-sas.md:327>)).
   이는 매우 작지만 bitwise exact는 아니며 전체 compute 이득도 작다. Linear-system
   conditioning과 trace-row identity 회귀가 갖춰지기 전에는 “동작 불변” 최적화로
   분류하지 않는 편이 안전하다.

## 7. 우선순위 권고

1. **먼저 P1 두 건을 닫는다.** species/line별 pump Rabi·saturation과 collision
   broadening/shift provenance는 lookup 비용만 추가하면서 가장 큰 false precision을
   없앤다. 특히 Cs의 Rb self-broadening fallback은 정량 사용을 막는 현재 최우선 물리 결함이다.
2. **User Guide와 이미지를 현재 UI로 갱신한다.** µW control, all-line generalized tier,
   Dispersion, 두 온도, paraffin, passivity qualifier를 반영해야 실험자가 올바른 범위를 읽는다.
3. **측정 validation은 metadata 계약부터 시작한다.** dark/reference, detuning 단위와
   calibration, cell/beam geometry, pump/probe polarization와 불확도를 갖춘 pump-off Rb D2
   dataset을 먼저 확보한다. 현재 로컬 CSV로 finite-power accuracy를 주장하지 않는다.
4. **기본 SAS scan을 linewidth reference로 자동 승격하지 않는다.** 현재 2.7
   samples/FWHM 진단을 유지하고, 정량 FWHM/lock slope는 명시적 local refinement 결과만
   export한다.
5. **순수 성능은 이미 큰 항목이 반영됐다.** 남은 exact 후보는 1–4% 수준이므로 물리
   provenance와 문서 동기화보다 뒤에 둔다.

## 최종 평가

Scheme 1은 단순 toy가 아니다. 실제 hyperfine structure, CG branching, velocity-correlated
pump/probe, thermal averaging, absolute OD normalization, weak-probe dispersion, Gaussian
finite-power propagation을 한 인터랙티브 경로에 결합한다는 점에서 실험 설계에 유용하다.
하지만 현재 가장 믿을 수 있는 정량 기준은 **pump-off ⁸⁵Rb D1 OD와 software-matched
pump-off Rb D2 finite-probe 범위**다. Pumped power calibration, Cs 고온 linewidth,
buffer-gas spectrum, joint pump/probe, coated-cell transport는 qualifier가 붙은 경향 예측으로
제한해야 한다. 다음 품질 상승은 더 큰 solver가 아니라 **종·선별 상수 provenance,
재현 가능한 실험 metadata, 문서와 UI의 일치**에서 얻는 편이 계산비용 대비 가장 크다.
