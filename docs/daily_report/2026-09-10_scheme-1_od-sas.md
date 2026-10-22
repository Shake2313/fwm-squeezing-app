# GABES 일일 Scheme 물리 검토 — 2026-09-10

## 1. 오늘의 선택과 한 줄 판정

- 현지 날짜는 2026-09-10, Asia/Seoul이고 day = 10이다.
- 따라서 \(n=(10\bmod 5)+1=1\)이며, 오늘 대상은 **Scheme 1 —
  Absorption spectroscopy (OD / SAS)** (`sas`)이다.
- 검토 기준 HEAD는 `c0c46f0`이다. 작업 트리에는 이 실행 전부터 Scheme 1을
  포함한 광범위한 staged/unstaged 변경이 있었다. 이번 실행은 그 변경을 수정하지
  않고 이 보고서만 추가했다.

**판정:** Scheme 1은 실제 원자·열·전파 물리를 담은 유용한 **반정량적 실험
계획/해석 도구**다. 특히 pump-off Rb 약한-probe 흡수, hyperfine 위치와 상대
세기, Doppler envelope, velocity-correlated hyperfine pumping, 약한-probe 분산,
유한 probe의 Gaussian 전파는 연구실에서 쓸 가치가 있다. 그러나 현 상태를
그대로 **절대 pumped-SAS, buffer/coated-cell, Cs 고온, 편광/자기장 또는 유한-power
실험 기준 모델**로 인용하면 안 된다.

이번 주기의 가장 큰 변화는 2026-09-06 작업 트리에 들어온 **thermal sublevel
population 중복 계수 제거**다. 이 교정은 oscillator-strength sum rule과 맞고
계산비용도 늘리지 않는다. 반면 다음 두 신뢰 계약은 아직 동기화되지 않았다.

1. 공통 `Optical density`는 \(-\log_{10}T\)인데 AutoOD와 원자물리 문헌의 OD는
   \(-\ln T=\alpha L\)다. 같은 스펙트럼을 `OD`라고 부르면서 값이 정확히
   \(\ln 10\)배 달라진다.
2. 절대 흡수 배열이 크게 바뀌었지만 `SASScheme.cache_version`은 여전히 `7`이다.
   앱 cache key와 export provenance가 이 값을 소비하므로 physics revision을
   올려야 한다.

## 2. 현재 다섯 scheme의 정의와 순서

드롭다운 순서는
[`gabes/schemes/__init__.py:19`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/schemes/__init__.py:19>)의
`_SCHEMES` 순서다.

| 순번 | registry key | 현재 표시명 | 정의 위치 |
|---:|---|---|---|
| 1 | `sas` | Absorption spectroscopy (OD / SAS) | [`sas.py:175`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/schemes/sas.py:175>) |
| 2 | `lambda` | Λ coherence (EIT / AT / CPT) | [`absorption.py:466`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/schemes/absorption.py:466>) |
| 3 | `rydberg_eit` | Rydberg-EIT electrometry | [`rydberg.py:94`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/schemes/rydberg.py:94>) |
| 4 | `magneto` | Magneto-optics (Hanle/MOR) | [`magneto.py:131`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/schemes/magneto.py:131>) |
| 5 | `fwm` | Four-wave mixing (Squeezing / Biphoton) | [`fwm.py:2986`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/schemes/fwm.py:2986>) |

`ODScheme`는 단일/85Rb-D1 검증 primitive로 남지만 등록 목록에는 없다. 실제 UI의
OD는 Scheme 1에서 `pump_power_mw=0`인 극한이다
([`__init__.py:12-24`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/schemes/__init__.py:12>)).

## 3. 기존 문서·일일 보고서·TODO·issue를 먼저 검색한 결과

새 기능을 상상하기 전에 다음을 확인했다.

- Scheme 1 선행 보고서 9건: 2026-06-25, 07-10, 07-20, 07-30, 08-05,
  08-10, 08-25, 08-30, 09-05. 최신 선행 보고서는
  [`2026-09-05_scheme-1_od-sas.md`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/docs/daily_report/2026-09-05_scheme-1_od-sas.md:1>)다.
- 통합 등록부
  [`docs/checklist.json`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/docs/checklist.json:794>)에는
  충돌계수, Rabi/saturation, fitting, isotope fraction, segmented cell, Zeeman,
  coating/VCC, export 제안이 이미 구체적으로 있다. 따라서 **기존 개선안이 없는
  경우가 아니다.**
- 코드의 직접 Scheme 1 TODO는 Ne broadening scalar를 gas/species/line별 표와
  pressure shift로 승격하자는 항목이다
  ([`gabes/constants.py:79-89`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/constants.py:79>)).
- 별도의 로컬 `issue`, `TODO`, `proposal` note 파일은 발견되지 않았다.
- `gh` CLI는 이 환경에 없었지만, 2026-09-10에 공개
  [GitHub Issues 페이지](https://github.com/Shake2313/fwm-squeezing-app/issues?q=is%3Aissue)를
  직접 확인했으며 issue 수는 0, 검색 결과도 없었다.

현재 Scheme 1 관련 checklist의 핵심 상태는 다음과 같다.

| 항목 | 상태 | 오늘의 판단 |
|---|---|---|
| collision provenance + pressure shift | P1 / large / `ready` | 여전히 최우선 물리 교정 |
| species/line Rabi + saturation provenance | P1 / medium / `ready` | 여전히 최우선 power-axis 교정 |
| sub-Doppler readout trust | P1 / small / `done` | 상태 gate는 동작하나 기본 grid는 아직 미수렴 |
| paraffin population memory | P2 / medium / `done` | opt-in 반정량 모델로 적절 |
| CSV calibration + sweep direction | P1 / medium / `done` | 절대/상대 의미와 branch 보존이 좋음 |
| dispersion, split temperature | P2 / `done` | 기본 solver를 무겁게 하지 않고 유용한 물리 추가 |
| experimental fitting | P2 / large / `ready` | 명시적 실행과 fit/holdout 분리가 필요 |
| isotope fraction | P2 / medium / `needs_decision` | Rb fraction 자체는 싸지만 새 원소는 데이터가 먼저 |
| segmented cell | P2 / large / `ready` | pump depletion에는 필요하나 segment 수만큼 비쌈 |
| weak-probe Zeeman, full VCC, coating 확장 | P2 / research / `parked` | 측정 target과 runtime budget 전에는 보류가 맞음 |
| result export/provenance | P2 / medium / `in_progress` | physics/default/readout version 정확성이 특히 중요 |

근거 위치는 각각
[`checklist:794`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/docs/checklist.json:794>),
[`832`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/docs/checklist.json:832>),
[`853`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/docs/checklist.json:853>),
[`874`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/docs/checklist.json:874>),
[`914`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/docs/checklist.json:914>),
[`1924`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/docs/checklist.json:1924>),
[`1960`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/docs/checklist.json:1960>),
[`1990`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/docs/checklist.json:1990>),
[`2024`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/docs/checklist.json:2024>),
[`2051`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/docs/checklist.json:2051>)이다.

## 4. 9월 5일 이후의 핵심 physics 변경: 열적 population을 한 번만 셈

### 4.1 교정 내용과 물리 판단

현재 전이별 약한-probe 적분 흡수는

\[
\int\alpha_t\,d\delta=
\frac{\pi k\,C_{F,t}^2|d|^2N}{\epsilon_0\hbar\,2(2I+1)}
\]

로 계산된다
([`species.py:248-269`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/species.py:248>)).
`C_F²`는 이미 초기 자기부분준위의 합을 포함한다. 따라서 각 sublevel의 열적
population인 `1/[2(2I+1)]`을 적용한 뒤 ground-manifold population `p_F`를 다시
곱하면 중복이다. 현재 코드는 그 추가 `p_F`를 제거했다.

이 판단은 단순히 기존 GABES 숫자에 맞춘 것이 아니다.

- 현재 테스트는 Rb-85, Rb-87, Cs-133의 D1/D2 전체에서 lifetime oscillator sum과
  두 ground manifold의 축퇴도 비를 독립 항등식으로 검사한다
  ([`tests/test_sas.py:91-109`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/tests/test_sas.py:91>)).
- Rb-85 D1은 F=3/F=2 적분비가 종전 49/25가 아니라 7/5가 된다
  ([`tests/test_sas.py:308-315`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/tests/test_sas.py:308>)).
- 수정된 별도 AutoOD-D2 reference도 population을 한 번만 쓰고 D2 lifetime sum으로
  spectrum 전체를 재구성한다
  ([`test_beam_saturation.py:186-212`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/references/AutoOD-Rb/test_beam_saturation.py:186>),
  [`AutoOD_Natural_Rb_D2_Model.tex:351-383`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/references/AutoOD-Rb/docs/AutoOD_Natural_Rb_D2_Model.tex:351>)).
- 공식 Steck Rb-85 표의 D1 line strengths는 각 Fg row에서 합이 1이며, 현재
  `line_strength()` 정의와 일치한다
  ([Steck Rb-85 D-line data](https://steck.us/alkalidata/rubidium85numbers.pdf)).

따라서 이번 정규화 변경은 **물리적으로 정당하며 계산량은 동일한 O(전이 수)**다.
배열의 scalar weight만 달라지고 solver 차원·velocity grid·solve count는 늘지 않는다.

### 4.2 현재 수치에 미친 영향

이 교정은 작은 문서 수정이 아니라 절대 흡수 baseline 변경이다. 최신 선행 보고서의
기본점은 peak OD 0.54, sub-Doppler FWHM 21.56 MHz, peak phase 421.4 mrad였다
([`2026-09-05 report:188-194`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/docs/daily_report/2026-09-05_scheme-1_od-sas.md:188>)).
현재 같은 기본점은 base-10 peak OD 0.928, FWHM 18.87 MHz, peak phase
688.3 mrad다. 즉 과거 report/그림의 절대 숫자는 현재 기준값으로 재사용하면 안 된다.

수정된 independent reference의 baseline 이유·source hash·환경은
[`thermal_normalization_baseline.json`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/references/AutoOD-Rb/thermal_normalization_baseline.json>)에
기록되어 있고 reference 테스트 29개가 통과했다. 다만 `.gitignore:19`가
`references/*`를 제외하므로 이 source와 NPZ는 일반 clone에서 자동 보존되는 공식
회귀 자산이 아니다
([`.gitignore:19`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/.gitignore:19>)).
제품 테스트에 복사된 scalar anchor와 independent sum-rule 테스트는 남지만, 장기
재현성을 위해 reference source/hash를 versioned artifact로 승격할지 결정해야 한다.

### 4.3 physics cache version이 갱신되지 않음

`SASScheme.cache_version`은 여전히 `7`이다
([`sas.py:175-180`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/schemes/sas.py:175>)).
앱의 heavy cache는 이 값을 key로 사용하고
([`streamlit_app.py:637-643`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/streamlit_app.py:637>)),
export에도 그대로 기록한다
([`streamlit_app.py:864-873`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/streamlit_app.py:864>)).
이번처럼 절대 spectrum이 바뀐 교정 뒤에는 `cache_version`을 새 physics revision으로
올리는 것이 필요하다. 비용은 O(1)이고 physics는 전혀 바꾸지 않는다.

## 5. 실제로 잘 구현된 실험 원자물리

### 5.1 hyperfine/열/절대 흡수

- `build_manifold()`는 isotope/line의 Fg/Fe 구조, hyperfine energy, CG-branched
  decay와 transit reset을 구성한다
  ([`species.py:337-416`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/species.py:337>)).
- pump-off에서는 별도 OBE population solve를 생략하고 `w=1`인 Voigt 합으로 정확히
  축약한다
  ([`sas.py:563-588`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/schemes/sas.py:563>),
  [`tests/test_sas.py:59-67`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/tests/test_sas.py:59>)).
- number density는 elemental density와 isotope abundance를 구분하고, 자연 Rb는
  두 isotope를 실제로 합산한다
  ([`species.py:192-207`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/species.py:192>),
  [`sas.py:365-415`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/schemes/sas.py:365>)).

### 5.2 SAS의 중요한 velocity correlation과 hyperfine pumping

Pump는 `scan + kv`, 역행 probe는 `scan - kv`를 보며, pump-prepared population과
probe Lorentzian을 같은 velocity class에서 곱한 뒤 평균한다
([`sas.py:590-631`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/schemes/sas.py:590>)).
이 상관을 두 개의 독립 Doppler 평균으로 factorize하지 않은 것은 실제 Lamb dip과
crossover를 만들기 위한 핵심이다. CG-branched spontaneous decay와 transit reload가
ground-hyperfine optical pumping을 만들므로, 단일 two-level hole-burning toy보다
실험적인 crossover 증강/반전을 더 잘 설명한다. 이 물리의 중요성은 Smith와 Hughes의
multilevel SAS 연구와도 맞는다
([DOI 10.1119/1.1652039](https://doi.org/10.1119/1.1652039)).

### 5.3 약한-probe dispersion과 두 온도

- absorption과 dispersion은 같은 complex profile, 같은 population difference,
  같은 Doppler quadrature에서 나온다
  ([`sas.py:539-553`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/schemes/sas.py:539>)).
- density temperature와 Doppler-width temperature를 기본적으로 묶되, cold spot과
  beam path가 다른 셀을 위해 분리할 수 있다
  ([`sas.py:56-75`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/schemes/sas.py:56>)).
- finite probe에서는 공간적으로 달라지는 nonlinear phase를 억지로 정의하지 않고
  Dispersion view를 숨긴다
  ([`sas.py:743-758`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/schemes/sas.py:743>)).

### 5.4 유한 probe와 실험 CSV

- Gaussian peak intensity, line-specific closed-two-level `I_sat`, power
  broadening, `dI/dz=-alpha(I)I`, 32-point radial quadrature가 구현되어 있다
  ([`probe_saturation.py:34-82`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/probe_saturation.py:34>),
  [`166-260`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/probe_saturation.py:166>)).
- Rb D2 reference-scope와 D1/Cs generalized estimate를 결과에 구분하고, pump+probe
  sequential approximation과 passivity limiter를 명시한다
  ([`sas.py:92-168`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/schemes/sas.py:92>),
  [`466-501`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/schemes/sas.py:466>)).
- CSV는 relative extrema normalization과 증거가 있는 absolute calibration을
  구분하고 forward/reverse acquisition branch를 보존한다
  ([`experimental_csv.py:214`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/experimental_csv.py:214>),
  [`443`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/experimental_csv.py:443>),
  [`tests/test_experimental_csv.py:167-224`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/tests/test_experimental_csv.py:167>),
  [`290-315`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/tests/test_experimental_csv.py:290>)).

## 6. 실험 reference로서 남은 핵심 한계

### 6.1 `OD`의 log convention이 코드와 reference에서 다름 — 새 P1 trust issue

공통 함수는

```text
Optical density = -log10(T) = alpha*L/ln(10)
```

를 반환한다
([`observables.py:540-547`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/observables.py:540>)).
반면 수정된 AutoOD 기술 문서는 `OD=-ln(T)=alpha*L`을 쓴다
([`AutoOD_Natural_Rb_D2_Model.tex:399-406`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/references/AutoOD-Rb/docs/AutoOD_Natural_Rb_D2_Model.tex:399>)).

직접 계산한 85Rb D1, 50 °C, 12.5 mm 값은 다음과 같다.

| 동일 spectrum | peak |
|---|---:|
| natural optical depth `alpha*L = -ln(T)` | `1.2556791912` |
| 현재 GABES `Optical density = -log10(T)` | `0.5453345438` |
| 비 | `2.302585093 = ln(10)` |

이는 raw absorption physics 오류는 아니지만 실험 비교·CSV export·논문 표에서 factor
2.3026 오류를 만들 수 있다. 최소 무비용 조치는 기존 값을 **Decadic absorbance
`A10=-log10(T)`**로 명확히 바꾸고, 원자물리 관례의 **Optical depth
`tau=-ln(T)=alpha L`**를 함께 내보내는 것이다. 호환성을 위해 기존 배열을 유지해야
한다면 label/version을 바꾸고 둘 사이의 exact conversion test를 추가하면 된다.

### 6.2 pump power axis는 아직 species/line/polarization별이 아님

`pump_rabi_from_power()`는 모든 선택에 Rb-85 D1 상수 `I_SAT=4.484 mW/cm²`를 쓴다
([`constants.py:22-31`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/constants.py:22>),
[`species.py:272-282`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/species.py:272>)).
그 결과 0.5 mW, 1 mm waist에서 natural/pure Rb, Cs의 D1/D2 모두 정확히
`Omega/Gamma = 1.8839842329`다. Steck의 Rb-85 표만 해도 D1 far-detuned effective,
D2 isotropic, D2 cycling `I_sat`가 서로 다르다. 따라서 현재 pump power 축은 trend
knob이지 정밀 Rabi calibration이 아니다. checklist의 P1 항목은 lookup/폐형식
계산만 추가하므로 runtime O(1)로 교정 가능하다. 단, 편광과 open transition의
effective convention을 먼저 선언해야 한다.

또한 현재 `rabi_rel`은 `T=(2Fg+1)S`의 제곱근을 쓴다
([`species.py:365-407`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/species.py:365>)).
Rb-85 D1의 Fg=3 가지는 per-ground-atom `sqrt(S/Smax)`와 비교해 각각
`1.000 vs 0.845`, `0.894 vs 0.756`이다. full mF/polarization average가 없는 lumped
Hamiltonian에서 어느 effective Rabi를 써야 하는지는 별도 provenance/검증 문제다.
새 thermal absorption normalization이 맞다는 사실이 이 pump coupling까지 자동으로
검증하지는 않는다.

### 6.3 충돌 넓힘이 pump와 probe에 일관되게 들어가지 않음

- `self_broadened_gamma(iso,N)`는 `iso`를 받지만 모든 종/선에 하나의 Rb 계수를 쓴다
  ([`species.py:210-212`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/species.py:210>)).
- 이 fallback을 Cs에 적용하면 Gamma_self/2pi가 100/150/200 °C에서 각각
  `1.0186 / 15.1751 / 124.0098 MHz`가 된다. 고온에서는 spectrum을 지배하는
  값이므로 단순 작은 오차가 아니다.
- Ne는 `3.91 MHz/Torr` broadening-only scalar이며 pressure shift, 온도 의존성,
  Dicke narrowing/VCC가 없다
  ([`constants.py:79-89`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/constants.py:79>)).
- `gamma_eff`는 probe Lorentzian과 pump detuning-grid spacing에는 쓰이지만,
  pump Liouvillian의 optical coherence dephasing에는 들어가지 않는다
  ([`sas.py:393-413`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/schemes/sas.py:393>),
  [`568-585`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/schemes/sas.py:568>)).

Rb D1 self-broadening 계수 자체는 Weller 등의 절대흡수 실험에 근거가 있지만
([Weller et al.](https://arxiv.org/abs/1107.3092)), Cs/D2와 buffer-gas 경로로 그대로
확장할 근거는 아니다. species/line/gas별 table, sign/unit/validity/uncertainty, pressure
shift와 elastic optical dephasing을 분리하면 solver 차원과 solve 수를 그대로 둔 채
물리를 개선할 수 있다.

### 6.4 finite pump + finite probe, uniform cell, mF와 coating

- finite probe는 pump-prepared population을 고정한 후 별도 two-level broadening을
  적용한다. probe-induced hyperfine repumping과 coherent joint solve가 없다.
- pump Rabi는 셀 입구 power에서 한 번 계산되고 z 방향 pump depletion/gradient가
  없다. 고 OD에서 SAS contrast는 낙관적일 수 있다.
- scalar hyperfine manifold라 mF, polarization, magnetic field, Faraday/Voigt
  dichroism·rotation을 예측하지 않는다.
- paraffin은 fixed `T1=25.1 ms`와 기존 transit rate를 사용하는 one-zone population
  memory다
  ([`sas.py:313-319`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/schemes/sas.py:313>)).
  이 값은 한 87Rb coated cell 측정치이며 universal coating constant가 아니다
  ([Bandi et al. primary PDF](https://libra.unine.ch/server/api/core/bitstreams/00e19b8f-092e-4821-8a8d-d21291de1f50/content)).

## 7. 오늘의 수치 spot check

동일 머신의 직접 실행값이다. `OD10`은 현재 UI의 base-10 값, `tau`는 자연로그
optical depth다.

| 확인 | 결과 | 해석 |
|---|---:|---|
| 기본 SAS compute, 1401점 | `0.88–1.04 s` | natural Rb D1, pump 0.5 mW, 40 °C, 75 mm |
| 기본 Tmin / peak tau / peak OD10 | `0.11805 / 2.13665 / 0.92794` | OD convention을 반드시 함께 기록해야 함 |
| 기본 sub-Doppler FWHM | `18.868 MHz`, `2.342 samples/FWHM` | `resolution-limited`; 정량 선폭 아님 |
| UI 최대 4001점 | `17.656 MHz`, `6.261 samples/FWHM`, `2.45 s` | status는 resolved이나 8x 기준보다 2.0% 큼 |
| 고정-window 5601 / 11201점 | `17.3596 / 17.3108 MHz` | 마지막 두 FWHM 차 0.282%; slope 차 약 2.64% |
| Gaussian Doppler FWHM | `518.694 MHz` | 한 optical line의 계산된 Gaussian 폭 |
| 기본 peak phase | `688.263 mrad` | single-pass phase, group delay가 아님 |
| pump-off 1401점 | `0.57 s`, Tmin `0.10510`, peak OD10 `0.97842` | pump OBE를 생략하는 exact fast path |
| natural Rb D2, OD default, probe 100 uW | effective peak OD10 `0.42`, peak s0 `3.813` | Rb-D2 software-reference scope |
| natural Rb D1, OD default, probe 6 uW | effective peak OD10 `0.96`, peak s0 `0.255` | generalized qualitative estimate |
| paraffin, 401점, pump on | `0.393 -> 1.066 s` (`2.71x`) | opt-in 연구 경로로 적절 |
| paraffin, 401점, pump off | 배열 exact equal | coating throughput을 잘못 곱하지 않음 |

기본 resolution status와 lock-slope gate는
[`sas.py:812-835`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/schemes/sas.py:812>)에,
4x/8x 고정-window 수렴 테스트는
[`tests/test_sas.py:442-494`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/tests/test_sas.py:442>)에 있다.
Gate가 잘못된 숫자를 숨기는 것은 좋지만, `resolved`는 physics convergence가 아니라
최소 sampling gate임을 기억해야 한다.

## 8. 문서·테스트·예제 평가

### 8.1 문서

- README는 thermal population once, 7/5 ratio, finite-probe tier, split
  temperature, dispersion과 주요 한계를 현재 코드에 맞게 많이 갱신했다
  ([`README.md:289-336`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/README.md:289>)).
- 그러나 첫 scheme 표와 in-app info는 여전히 `validated 85Rb D1 hyperfine scale`,
  `lab AutoOD ... <0.1%`라고 쓴다
  ([`README.md:12`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/README.md:12>),
  [`sas.py:296-306`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/schemes/sas.py:296>)).
  현재 증거는 **normalization-audited/software-reference matched**이지, 별도 held-out
  experiment로 절대 정확도를 검증한 것은 아니다.
- HTML User Guide는 현 코드보다 뒤처져 있다. natural Rb D2에서만 finite probe가
  보인다고 하고 power 기본값/단위를 0 mW로 쓴다
  ([`GABES_User_Guide_v2.html:551-584`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/docs/Userguide/GABES_User_Guide_v2.html:551>)).
  실제 UI는 모든 built-in D line에 0–100 uW를 제공하고 D1/Cs를 generalized로
  표시한다
  ([`sas.py:188-213`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/schemes/sas.py:188>)).
- User Guide 표는 Dispersion/phase, paraffin, 두 온도, base-10 OD convention도
  충분히 드러내지 않는다
  ([`GABES_User_Guide_v2.html:524-530`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/docs/Userguide/GABES_User_Guide_v2.html:524>)).

### 8.2 예제

Scheme 1의 독립 실행 예제/config는 없고 문서 예제는
[`od.png`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/docs/Userguide/userguide_assets/od.png>)와
[`sas.png`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/docs/Userguide/userguide_assets/sas.png>)가 전부다.
둘 다 2026-06-08 생성본이며:

- 종전 weaker absorption baseline을 보여 현재 thermal normalization과 다르고,
- Transmission과 OD를 한 figure에 쌓아 현재 carousel 구조와 다르며,
- Dispersion, finite-probe qualifier, resolution status와 provenance를 보여 주지 않는다.

따라서 시각적 개념 예제로는 쓸 수 있지만 **현재 numerical example** 또는 figure
regression으로는 부적합하다. current params, physics/cache version, git revision을
동반한 재생성 예제 하나가 필요하다.

### 8.3 테스트

실행 결과:

- `python -m pytest -q tests/test_sas.py tests/test_absorption.py tests/test_probe_saturation.py tests/test_experimental_csv.py`
  -> **132 passed in 103.95 s**.
- `python -m pytest -q references/AutoOD-NatRbD2/test_beam_saturation.py`
  -> **29 passed in 7.98 s**.
- render/headless/UI/docs 묶음은 **33 passed, 1 failed**. 실패는 Scheme 1이 아니라
  작업 트리에서 삭제 상태인
  `docs/FWM physics and analytic reconstruction/FWM_physics.tex`를
  `test_docs_consistency.py`가 읽으려 한 `FileNotFoundError`다.
- 저장소 지침의 `python -m pytest -q`도 실행했지만 제품 테스트 전에
  `references/AutoOD/test_autood_gui_adjustments.py`가 기본 requirements에 없는
  optional `sympy`를 import하여 collection에서 중단됐다: **1 collection error in
  44.71 s**.

좋은 점은 oscillator sum, 7/5, pump-off skip, coating invariance, dispersion sign,
split temperature, finite-probe reference anchors, passivity, CSV absolute calibration과
sweep branch를 폭넓게 검사한다는 것이다. 다만
[`test_pump_off_reproduces_autood_85rb_d1`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/tests/test_sas.py:260>)은
같은 저장소의 `ODScheme`와 비교하므로 독립 검증이 아니다. 별도 AutoOD test와
analytic sum rule이 이를 보완하지만, held-out 실험 spectrum은 여전히 없다.

## 9. 기존 개선안의 계산비용과 physics 보존성

| 기존 개선안 | 계산비용 | 기존 물리를 거의 무비용으로 보존 가능한가 | 권고 |
|---|---|---|---|
| thermal population once 교정 | solver/배열 크기 동일, weight O(1) | **예.** 이미 물리적으로 올바르게 반영됨 | cache/doc/example/version 동기화 후 확정 |
| collision provenance + pressure shift | coefficient lookup/offset O(1); optical dephasing도 같은 Liouvillian 차원 | **예.** zero-pressure 및 Rb-85 D1 anchor를 exact 고정 가능 | P1 |
| species/line/polarization Rabi | lookup/폐형식 O(1) | **예.** pump-off는 정확히 불변; pump-on power axis만 교정 | P1 |
| OD/log convention 명시 | array 하나의 O(N) scale 또는 label만 변경 | **예.** alpha/T는 불변, exact ln(10) 변환 | 새 P1 trust item으로 등록 권고 |
| export/provenance completion | O(Nscan) serialization, 추가 solve 없음 | **예.** 단 physics/cache version을 먼저 올바르게 해야 함 | P1에 가깝게 선행 |
| isotope fraction | 기존 두 Rb component의 weight만 변경 | **예.** natural/pure defaults exact 유지 가능 | 실제 enriched cell 수요 결정 후 |
| spectrum fitting | 대략 `N_eval x compute`; 현 50회면 pump-off 약 25–30 s, pump-on 약 45–55 s | 기본 interactive 경로에는 **아니오** | explicit bounded action + cache + holdout |
| segmented cell/pump depletion | 대략 segment 수에 선형; self-consistent pump state면 더 큼 | 1 segment exact는 가능, multi-segment는 비쌈 | opt-in convergence mode |
| calibrated coating/VCC/mF-Zeeman | 상태/velocity/geometry 차원 증가 | **아니오** | 측정 target·held-out data·runtime budget 전에는 parked 유지 |

## 10. 동작을 바꾸지 않는 순수 코드 최적화

### 10.1 pump-off scalar population fast reduction — 우선순위 높음

현재 helper는

```python
return (prepared_population * probe_profile) @ weights
```

를 항상 실행한다
([`sas.py:996-998`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/schemes/sas.py:996>)).
Pump-off에서는 `prepared_population`이 정확히 scalar `1.0`이므로 큰
`scan x velocity` 임시 배열을 만들 필요가 없다. `probe_profile @ weights`로 바로
가면 된다.

동일 process warm median에서 runtime은 다음과 같았다.

| 경로 | 현재 | scalar fast | speedup | 결과 |
|---|---:|---:|---:|---|
| pump off, probe 0 uW | `0.4671 s` | `0.3655 s` | `1.278x` | `alpha_unit` array exact |
| pump off, probe 100 uW | `0.5249 s` | `0.4343 s` | `1.209x` | alpha와 finite table array exact |

대표 1401x1403 transition reduction 자체는 `6.48 -> 0.64 ms`였다. 이 분기는
수학·부동소수 결과가 exact하고 pump-on에는 닿지 않으므로 가장 안전한 순수 최적화다.

### 10.2 pump-off finite-probe thermal matvec 재사용

Pump-off에서 `transition_alpha`와 `thermal_alpha`는 같은 `A_t*(Lp@wt)`인데 현재
연속으로 두 번 계산한다
([`sas.py:629-638`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/schemes/sas.py:629>)).
`thermal_alpha = transition_alpha`를 쓰면 exact하며 finite-probe OD의 transition당
matvec 하나를 제거한다. 10.1보다 작지만 함께 적용하기 좋다.

### 10.3 32-point Gauss-Legendre nodes/weights cache

`propagate_gaussian_spectrum()`은 navigate-only `cell_mm` 또는 `line_strength`가
바뀔 때마다 `leggauss(32)`를 다시 만든다
([`probe_saturation.py:166-199`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/probe_saturation.py:166>)).
32점 nodes/weights를 read-only module constant나 `lru_cache`로 재사용하면 된다.
측정상 headless finite-probe observables median은 `22.67 -> 21.50 ms` (`1.055x`)였고
metric은 exact equal이었다. 전체 heavy solve에 비하면 작지만 navigate interaction에는
반복적으로 이득이다.

### 10.4 낮은 우선순위

`line_integrated_alpha()`의 Wigner algebra는 기본 natural Rb compute에서 isotope당
약 217 us에 불과하다. `(iso,line)` 단위 per-density coefficient cache는 exact하게
가능하지만 전체 약 0.9 s solve의 0.05% 미만이므로 우선순위가 낮다. Pump-on의
population interpolation을 한 번에 vectorize하는 실험은 같은 배열에서 오히려
`0.17 -> 0.51 s`로 느려졌고 80 MB급 임시 배열을 만들어 권고하지 않는다.

## 11. 권고 순서

1. **신뢰 계약을 먼저 고친다.** `Optical depth=-ln(T)`와
   `Decadic absorbance=-log10(T)`를 구분하고 UI/export/test에 exact conversion을
   고정한다. thermal normalization 변경에 맞춰 Scheme 1 `cache_version`을 올린다.
2. **정규화 교정의 증거를 정식 보존한다.** ignored reference에만 있는 source/hash와
   baseline을 versioned artifact 또는 재현 가능한 fixture로 승격하고, README/info의
   `validated`를 `normalization-audited/software-reference matched`로 낮춘다.
3. **기존 P1 두 건을 실행한다.** species/line/gas별 collision/shift/dephasing과
   pump Rabi convention을 출처·유효범위·불확도와 함께 넣는다. 둘 다 solver 차원을
   늘리지 않으며 pump-off reference를 그대로 보존할 수 있다.
4. **기본 해상도는 status로만 끝내지 않는다.** 현재 1401점 FWHM은 11201점 대비
   약 9.0% 크다. 정량 lock/linewidth export에는 4001점 또는 bounded local refinement와
   convergence 값까지 저장한다.
5. **순수 성능 개선**은 scalar pump-off reduction, 중복 thermal matvec 제거,
   Legendre cache 순으로 적용한다. 모두 physics 결과를 exact 보존할 수 있다.
6. fitting, segmented pump depletion, VCC, full Zeeman/coating 확장은 현재처럼 opt-in
   또는 parked로 유지하고 먼저 실제 calibration/holdout data와 runtime budget을 정한다.

## 12. 최종 결론

Scheme 1은 “가짜 모양 생성기”가 아니다. current thermal normalization은
oscillator-strength sum과 맞고, pump/probe velocity correlation, hyperfine pumping,
Voigt absorption/dispersion, finite-probe Gaussian propagation, CSV calibration이 실제
실험에 유용하다. 특히 **pump-off Rb 약한-probe 흡수와 qualitative SAS 구조**는
좋은 참고점이다.

그러나 실험자가 지금 숫자를 가져갈 때는 세 가지를 반드시 함께 적어야 한다.

- `OD`가 natural optical depth인지 base-10 absorbance인지,
- pump power/Rabi와 collision coefficient가 어떤 species/line/polarization convention인지,
- sub-Doppler feature가 grid-converged인지 단지 resolution gate를 통과했는지.

이 세 계약을 명확히 하고 physics/cache version과 예제를 갱신하면, 큰 solver를 새로
만들지 않고도 실험 reference로서의 신뢰도를 가장 크게 높일 수 있다.
