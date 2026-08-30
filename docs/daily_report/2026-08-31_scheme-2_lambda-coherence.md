# 2026-08-31 Scheme 2 물리 검토 — Lambda coherence (EIT / AT / CPT)

## 1. 오늘의 선택과 현재 다섯 scheme

서울 현지 날짜의 day-of-month는 `31`이므로

```text
n = (day mod 5) + 1 = (31 mod 5) + 1 = 2
```

이다. UI 드롭다운의 실제 순서는 `_SCHEMES`가 정한다
(`gabes/schemes/__init__.py:12-25`). 현재 정의와 순서는 다음과 같다.

| 순번 | 등록 인스턴스 | UI scheme | 주된 물리 |
|---:|---|---|---|
| 1 | `SASScheme()` | OD / SAS | Doppler 흡수, 포화흡수, hyperfine optical pumping |
| 2 | `LambdaScheme()` | Λ coherence (EIT / AT / CPT) | 축약 3준위 Λ 결맞음 |
| 3 | `RydbergEITScheme()` | Rydberg-EIT electrometry | cascade EIT, microwave AT, electrometry |
| 4 | `MagnetoScheme()` | Hanle / EIA / NMOR | Zeeman OBE의 투과·회전 |
| 5 | `FWMScheme()` | FWM | seeded mean-field FWM, SFWM biphoton |

`eit`, `at`, `cpt`는 테스트와 직접 호출을 위한 alias이고 드롭다운 항목은 아니다
(`gabes/schemes/__init__.py:27-39`). README의 현재 scheme 표도 같은 다섯 엔진을
기술한다 (`README.md:8-16`). 따라서 오늘 검토 대상은 2번 `LambdaScheme()`다.

검토 시점은 branch `main`, HEAD
`ce60b5c029495400e014c0f1ce7d281b17046d60`이다. 작업 트리에는 기존 untracked
8월 29·30일 보고서가 있었고 보존했다. Lambda 본체
`gabes/schemes/absorption.py`의 마지막 변경은 2026-07-20 commit `84aba45`다.
8월 26일 보고서가 검사한 working-tree Lambda 동작과 오늘의 수치 결과 사이에
의미 있는 변화는 발견하지 못했다.

## 2. 먼저 검색한 기존 제안, TODO, issue note

먼저 다음을 검색·확인했다.

- Scheme 2 선행 보고서 9개: `2026-06-26`, `07-01`, `07-06`, `07-16`,
  `07-31`, `08-01`, `08-06`, `08-21`, `08-26`
  (`docs/daily_report/*scheme-2_lambda-coherence.md`)
- 통합 작업 등록부 `docs/checklist.json`과 코드 TODO
- 코드: `gabes/schemes/absorption.py`, `gabes/atoms.py`, `gabes/beam.py`,
  `gabes/kernels.py`, `gabes/lineshape.py`, `gabes/observables.py`
- 테스트: `tests/test_absorption.py`, `tests/test_kernels.py`,
  `tests/test_headless_observables.py`, `tests/test_schemes_render.py`,
  `tests/test_docs_consistency.py`
- 문서·예제: `README.md`, `CLAUDE.md`,
  `docs/Userguide/GABES_User_Guide_v2.html`, `eit.png`, `at.png`, `analysis/`,
  `examples/`
- 로컬 파일명과 본문의 `TODO`, `issue`, `proposal`, `deferred`, `parked` 검색
- 공개 GitHub의 [open issue](https://github.com/Shake2313/fwm-squeezing-app/issues)와
  [closed issue](https://github.com/Shake2313/fwm-squeezing-app/issues?q=is%3Aissue+state%3Aclosed)

별도의 Lambda 전용 TODO/issue/proposal 파일과 실행 가능한 Lambda 예제는 없다.
공개 GitHub issue도 오늘 확인 시 open 0건, closed 0건이다. 그러나 기존 개선안은
통합 checklist와 선행 보고서에 충분히 구체적으로 남아 있다.

| 기존 항목 | 상태 | 오늘의 판단 |
|---|---|---|
| `lambda-spectrum-validity` | P1, `ready` | cold-limit residual-Doppler 오류와 미수렴 linewidth/group-index가 그대로 재현됨. 최우선 (`docs/checklist.json:479-501`) |
| `collisional-coefficient-provenance-and-pressure-shift` | P1, `ready` | gas/species/line별 출처·부호·단위·불확도 없이 scalar 폭을 일반화하면 false precision. 계산은 싸지만 조사·검증 effort가 큼 (`docs/checklist.json:209-243`) |
| `paraffin-coated-cell-extensions` | P2, `parked` | OD/SAS population reservoir를 Lambda에 복사하면 안 되고 ground coherence와 bright/dark Ramsey dynamics가 필요. geometry·dark-time·holdout data 전에는 유지 보류 (`docs/checklist.json:304-325`) |
| `full-velocity-changing-collision-kernel` | P2, `parked` | 선택한 scheme/kernel/dataset/runtime budget이 있을 때만 opt-in reference solver로 시작 (`docs/checklist.json:603-622`) |
| `lambda-hyperfine-resolved-manifold` | P2, `parked` | isotope/transition/polarization과 holdout spectrum을 특정한 뒤 opt-in 구현 (`docs/checklist.json:625-643`) |

코드의 직접 TODO는 `gabes/constants.py:79-89`의 Ne broadening-only scalar를
gas/species/line table과 pressure shift로 승격하라는 항목이다. 오늘은 새 대형 기능을
더 제안하지 않고, 위 제안의 물리 가치와 계산비용을 현재 코드에서 재검증했다.

## 3. 현재 구현이 실제로 푸는 물리

### 3.1 정상상태 3준위 Lindblad OBE

`atoms.lambda3()`는 `g1`, `g2`, `e`의 세 상태를 만들고 들뜬 상태가 두
바닥상태로 각각 `Γ/2`씩 붕괴하도록 한다. `gamma_gg`는 population exchange가
아니라 `rho_g1,g2`와 그 켤레만 감쇠시키는 Raman-coherence `T2` rate다
(`gabes/atoms.py:168-190`; 실제 dissipator 적용은 `gabes/atoms.py:90-92`).

Hamiltonian은 약한 probe `Ω_p = 10^-3 Γ`, coupling field, coupling detuning과
two-photon scan을 직접 포함한다 (`gabes/schemes/absorption.py:29-31, 633-656`).
따라서 다음은 임의의 선모양 합성이 아니라 정상상태 density matrix에서 나온다.

- destructive interference가 만드는 EIT transparency와 정상분산
- 강한 coupling의 dressed states가 만드는 Autler–Townes doublet
- Raman-coherence dephasing이 제한하는 dark resonance
- coupling power와 beam diameter에 따른 `Ω_c ∝ sqrt(P)/d` 경향
  (`gabes/schemes/absorption.py:584-591`; `gabes/beam.py:9-20`)

이는 실험적으로 유용한 실제 원자물리다. 단, power-to-Rabi 관계는 transition
dipole과 polarization에서 절대적으로 계산하지 않고 사용자가 정한
`coupling_rabi_mhz` anchor를 `sqrt(P)/d`로 스케일한다. 따라서 power trend는
유용하지만 서로 다른 isotope/line/polarization의 절대 Rabi calibration은 아니다.

또한 CPT는 균형 잡힌 두 optical field를 쓰는 별도 clock-CPT 엔진이 아니라 같은
weak-probe Lambda 모델의 좁은 scan preset이다
(`gabes/schemes/absorption.py:482-500, 624-631`). Clock light shift, magnetic
sensitivity, optical pumping을 정량 예측하는 reference로 읽으면 안 된다.

### 3.2 scalar D-line medium

Rb/Cs와 D1/D2 선택은 온도별 밀도, 자연선폭, 질량, 파장과 reduced dipole을
바꾼다 (`gabes/schemes/absorption.py:418-443`). 그러나 구체적인 `Fg -> Fe`,
Clebsch–Gordan 계수, Zeeman population, polarization selection, hyperfine optical
pumping은 없다.

특히 `Rb (natural)`은 85Rb와 87Rb susceptibility를 각각 계산해 비결맞음 합산하지
않는다. abundance가 가장 큰 isotope의 line constant와 abundance-weighted
density/mass를 하나의 scalar Lambda에 넣는다
(`gabes/schemes/absorption.py:430-442`; mixture 정의는 `gabes/species.py:176-182`).
따라서 특정 isotope/hyperfine transition의 절대 contrast·linewidth reference가
아니라 **유효 scalar medium**이다. 사용자 가이드는 이 점을 비교적 정직하게
밝힌다 (`docs/Userguide/GABES_User_Guide_v2.html:828-830`).

### 3.3 Maxwell 평균, Beer–Lambert, group index

Doppler-on Lambda는 4-sigma Maxwell grid를 `dv = 2 m/s`로 만들고 scan × velocity의
9×9 real-Hermitian Liouvillian을 푼 뒤 probe coherence만 가중 수축한다
(`gabes/schemes/absorption.py:118-140`; `gabes/kernels.py:437-518`). 기본 50 °C
EIT는 601 scan point와 711 velocity class를 사용한다.

그 뒤 `chi_bar`는 밀도와 dipole을 포함한 `chi_phys`,
`alpha = k Im(chi_phys)`, `T = exp(-alpha L)`로 변환된다
(`gabes/observables.py:517-542`). Group index는 dilute-medium 근사의
`n_g = n + omega dn/domega`를 국소 수치 미분한 값이다
(`gabes/observables.py:550-558`). CW susceptibility와 local dispersion 계산으로는
타당하지만 finite-pulse delay, distortion, bandwidth, detector SNR을 푸는
slow-light propagation model은 아니다.

## 4. 실험 원자물리 reference로서의 적합성

현재 Scheme 2는 **실제 물리를 구현한 반정량적 실험 planning/reference**다.
교과서 모양을 임의로 그리는 도구보다 훨씬 유용하지만, 특정 alkali Lambda 실험의
절대 spectroscopy·clock·buffer-cell·slow-light reference는 아니다.

| 실험 질문 | 적합성 | 이유 |
|---|---|---|
| coupling power/diameter 변화와 EIT→AT 경향 | 좋음 | 정상상태 OBE와 anchor-preserving Rabi scaling |
| AT split 대 측정 `Ω_c` sanity check | 좋음 | resolved local absorption maxima와 `split ≈ Ω_c`를 직접 검사 |
| Raman `T2` 증가가 EIT/CPT를 흐리는 방향 | 유용 | 단일 coherence-dephasing rate이며 transit, diffusion, population exchange와 분리되지 않음 |
| warm-vapor beam alignment sensitivity | 오류 수정 후 유용 | residual-Doppler 방향은 유용하지만 현재 cold-limit 위반 |
| 특정 isotope/hyperfine/polarization의 절대 contrast·width | 낮음 | scalar 3-level, symmetric branching, Zeeman/optical pumping 부재 |
| clock-CPT light shift·자기장 민감도 | 낮음 | balanced bichromatic CPT가 아닌 weak-probe preset |
| slow-light pulse delay·왜곡·효율 | 낮음 | local `n_g`만 있고 pulse propagation·검출 모델이 없음 |
| coated-cell EIT/Ramsey narrowing | 현재 부적합 | ground coherence 보존, dark-time 분포, cell/beam geometry가 없음 |

온증기 EIT에서 velocity-changing collisions, diffusion/transit, Raman dephasing은
서로 대체할 수 없는 채널이다. 따라서 `buffer_ground_relax_khz` 하나를 모든
buffer-cell 이동·충돌 물리로 해석해서는 안 된다. Coated-cell EIT도 bright/dark
interval과 Ramsey phase를 요구하므로 OD/SAS의 population-only reservoir를 복사하는
방식은 물리적으로 부적절하다.

## 5. 2026-08-31 직접 재계산

### 5.1 기본 출력과 실행시간

Numba warm-up 뒤 각 기본값을 7회 계산한 median이다. 절대 시간은 현재 환경에만
해당한다.

| regime | compute median | 현재 hero/readout |
|---|---:|---|
| EIT | `34.258 ms` | `T_res = 0.0135583`, FWHM `0.45968 MHz`, `n_g = 100,305` |
| AT | `3.204 ms` | split `46.0 MHz`, expected `46.0 MHz`, `T_center = 0.989090` |
| CPT | `3.140 ms` | `T_res = 0.706328`, FWHM `923.32 kHz`, `n_g = 95,255` |

AT의 strong-coupling scale은 유용하다. EIT/CPT width와 EIT `n_g`는 아래 수치
해상도 문제 때문에 정밀 숫자로 사용하면 안 된다.

### 5.2 P1 — 기본 EIT linewidth와 group index 미수렴

모든 regime은 601점 uniform scan을 쓴다
(`gabes/schemes/absorption.py:609-631`). `window_fwhm()`은 half-height edge를
보간하지 않고 바깥 sample 사이의 거리를 반환한다
(`gabes/lineshape.py:53-71`). 기본 EIT 표시 FWHM은 중앙 grid의 정확히 두 칸이다.

| scan | 중앙 step | 표시/sample FWHM | edge-interpolated FWHM | center `n_g` | warm solve |
|---|---:|---:|---:|---:|---:|
| uniform 601 | `0.229840 MHz` | `0.459680 MHz` | `0.229899 MHz` | `100,305` | `32.921 ms` |
| `sinh(5u)` center-clustered 601 | `0.015488 MHz` | `0.186156 MHz` | `0.159141 MHz` | `112,708` | `33.653 ms` |
| uniform 9,601 | `0.014365 MHz` | `0.172380 MHz` | `0.159102 MHz` | `112,693` | `479.810 ms` |

현재 hero FWHM은 수렴 보간값보다 약 `2.89×` 크고 `n_g`는 약 `11%` 작다.
같은 601점을 중심에 모으면 9,601점 기준에 대해 FWHM `0.025%`, group index
`0.013%` 안에서 일치한다. 오늘 runtime은 uniform 601보다 약 `2.2%` 늘어
checklist의 `<=20%` budget을 충분히 만족한다.

공용 `subdoppler_feature()`는 보간 edge, samples/FWHM, scan-edge clearance와
status를 이미 계산한다 (`gabes/lineshape.py:141-249`). 현재 기본 EIT에 적용하면
`resolution-limited`, `1.0003 samples/FWHM`이다. Lambda의 transparency-floor
정의는 SAS running-median background와 다르므로 구조와 status contract를
재사용하되 feature definition은 Lambda용으로 명시하는 편이 안전하다. Checklist는
numeric hero에 최소 8 samples/FWHM을 요구한다 (`docs/checklist.json:492-499`).

### 5.3 P1 — residual Doppler의 cold-limit 위반

현재 affine 조립은

```text
A_coef = dL/ds - S_v
B_coef = S_v
L(s, kv) = base + s A_coef + kv B_coef
```

이다 (`gabes/schemes/absorption.py:45-57`; `gabes/kernels.py:477-490`).
Angle-dependent ground-level Doppler coefficient가 optical scan shift에도 들어가므로,
`v = 0`인 Doppler-off 계산에서도 angle이 complex susceptibility를 바꾼다. 실제
residual `Delta k · v`는 `v = 0`에서 정확히 0이어야 한다.

오늘 cold EIT의 0 대 10 mrad 비교는 다음을 재현했다.

- `max |Delta chi_bar| / max |chi_bar| = 2.942038%`
- `|Delta k| / k = 0.009999958`
- center transmission은 우연히 양쪽 모두 `0.925093`이지만 전체 complex curve는 다름

반면 warm 50 °C, 3 mm에서는 10 mrad가 공명 투과를
`0.423098 -> 0.00267946`로 낮춘다. Intended alignment physics는 강하고 유용하지만,
현재 operator bookkeeping이 cold invariant를 깨뜨린다.

Optical scan operator `S_opt`와 velocity-only residual operator `S_v`를 분리해
`A = dL/ds - S_opt`, `B = S_v`로 만들면 matrix 차원, scan 수, velocity class 수가
그대로다. 추가 solve 없이 정확성만 회복할 수 있다. Doppler-off 0/10 mrad 전체
complex array invariance를 checklist tolerance로 회귀해야 한다.

### 5.4 collinear optical-frequency difference 누락

`_two_photon_k_ratio()`는 probe와 coupling에 같은 `medium["k_vec"]`를 넘겨
angle 0에서 residual을 정확히 0으로 만든다
(`gabes/schemes/absorption.py:594-597`; `gabes/beam.py:36-55`). 실제 hyperfine
Lambda의 두 광장은 ground splitting만큼 주파수가 달라 co-propagating이어도 작은
종방향 `Delta k`가 남는다. 선행 계산의 25 °C Gaussian FWHM은 85Rb 약 `4.07 kHz`,
87Rb `9.07 kHz`, 133Cs `9.86 kHz`였다
(`docs/daily_report/2026-08-06_scheme-2_lambda-coherence.md:221-240`).

기본 EIT에는 작지만 협폭 CPT에는 무시하기 어렵다. Pure-isotope mode에서 실제 두
wavevector를 사용하면 coefficient scalar만 바뀌며 추가 solve가 없다. 앞의
cold-limit operator 분리가 먼저다.

### 5.5 큰 group index를 usable slow light로 읽으면 안 된다

기본 EIT의 local `n_g = 100,305`와 15 mm cell을 그대로 delay로 환산하면
`(n_g-1)L/c ≈ 5.019 us`다. 그러나 같은 점의 transmission은 `1.356%`뿐이고,
linewidth도 미수렴이다. 코드에는 pulse bandwidth, pulse reshaping, delay-efficiency,
detector SNR이 없다. 따라서 이 값은 **local CW dispersion diagnostic**이지 usable
slow-light delay 예측이 아니다. User Guide가 가파른 분산을 조건 없이 “느린 빛의
서명”이라고 부르는 문구 (`docs/Userguide/GABES_User_Guide_v2.html:595-597`)에는
최소한 transmission과 finite-pulse 한계를 함께 붙여야 한다.

## 6. 기존 물리 개선안의 비용과 물리 보존성

| 개선 | 물리 보존 여부 | 계산비용 | 권고 |
|---|---|---:|---|
| scan/velocity operator 분리 | 잘못 섞인 항만 분리, 모델 차원 동일 | 사실상 0 | 즉시 P1 |
| center-clustered 601 + 보간/status | 동일 OBE를 더 적절한 grid에서 평가 | 오늘 `+2.2%`; 후처리 `O(N)` | broad-wing/off-center 회귀 후 P1 |
| transmission·group-delay·visibility status | 기존 배열의 해석만 개선 | `O(N)`, 추가 solve 없음 | P1과 함께 |
| pure-isotope collinear `Delta k` | 실제 wavevector bookkeeping | scalar coefficient, 추가 solve 없음 | P2 |
| sourced pressure shift·elastic optical width | line centre·homogeneous width 개선 | table lookup, solve 차원 동일 | provenance 확보 후 P1 |
| coated-cell Lambda ground-coherence/bright-dark model | 현재 population-only OD/SAS 모델을 복사하지 않음 | dark-time/trajectory quadrature로 중간~높음 | geometry·T2·holdout trace가 있을 때 opt-in research |
| full VCC/diffusion | velocity redistribution·Dicke/Ramsey 물리 추가 | velocity classes 결합, separability 상실 | one-scheme/one-kernel/one-dataset opt-in 연구 solver |
| hyperfine/Zeeman Lambda | line assignment·편광·optical pumping 추가 | 상태 수와 Liouvillian 급증 | 특정 transition과 holdout spectrum 이후 opt-in |

저차 pressure shift와 elastic optical width는 runtime 관점에서는 싸다. Checklist가
effort를 `large`로 둔 이유는 연산량이 아니라 gas/species/line별 출처·부호·단위·
불확도·유효범위를 정하는 연구 작업 때문이다. 검증 데이터 없는 phenomenological
Dicke knob는 빠르더라도 reference 신뢰도를 낮추므로 추가하지 않는 편이 맞다.

Coated-cell Lambda도 단순한 `T2` exponential 하나는 계산상 싸지만, 그것만으로는
bright/dark return과 Ramsey phase를 검증하지 못한다. Cell/beam geometry와 dark-time
distribution을 넣으면 scan × dark-time 또는 trajectory quadrature가 필요해진다.
현재 데이터 없이 “negligible overhead” 구현을 서두르기보다 parked 상태를 유지하는
것이 물리 보존에 유리하다.

## 7. 동작을 바꾸지 않는 순수 코드 최적화

### 7.1 `line_strength`를 navigate-only로 이동 — 가장 큰 확실한 이득

`ParamSpec.recompute` 기본은 `True`다 (`gabes/schemes/base.py:44-56`). Lambda의
`line_strength`는 이를 그대로 상속한다 (`gabes/schemes/absorption.py:557-558`).
그러나 이 값은 Hamiltonian이나 `chi_bar`에 들어가지 않고 physical susceptibility를
만들 때만 곱해진다 (`gabes/schemes/absorption.py:664-682`;
`gabes/observables.py:517-537`).

기본 EIT에서 `line_strength = 1 -> 1.5`를 오늘 재검증한 결과:

- `chi_bar`는 bitwise identical
- 전체 재 solve median `33.399 ms`
- 같은 raw를 새 strength로 headless remap한 median `0.1162 ms`
- 최종 metric 동일, 불필요한 경로 차이 약 `287×`

`line_strength`에 `recompute=False`를 주고 `observables()`가 live
`params["line_strength"]`를 읽도록 하면 public output과 물리를 보존하면서
불필요한 solve를 없앨 수 있다. Cache equivalence, `chi_bar` 불변, `alpha` 선형
scaling을 회귀로 고정해야 한다. 이것이 가장 우선할 순수 성능 개선이다.

### 7.2 immutable medium/atom template cache — 안전하지만 낮은 우선순위

오늘 `_medium_from_params()` median은 `19.1 us`, `atoms.lambda3()`는 `196.3 us`였다.
Scalar key의 immutable medium record와 3-level dissipator/template를 작은 LRU로
캐시하면 반복 sweep의 Python/array construction을 줄일 수 있다. 합계는 기본 EIT
solve의 약 `0.63%`이므로 우선순위는 낮다. Mutable ndarray 공유와 round-key cache는
정확성 위험에 비해 가치가 작으므로 피해야 한다.

### 7.3 affine inner-loop hoist는 근거 부족

8월 26일의 bitwise-equivalent prototype은 current `29.299 ms`, hoisted
`29.732 ms`로 오히려 `1.5%` 느렸다
(`docs/daily_report/2026-08-26_scheme-2_lambda-coherence.md:297-309`). LU solve가
지배하고 temporary matrix가 이득을 상쇄한다. 여러 CPU의 프로파일에서 일관된
이득이 나오기 전에는 적용하지 않는 편이 맞다.

## 8. 문서, 테스트, 예제의 reference 품질

### 8.1 문서

- README의 Lambda 한 줄은 기능을 정확히 요약하지만 scalar/weak-probe/
  anchor-calibrated Rabi/natural-Rb effective-medium 한계를 말하지 않는다
  (`README.md:13`).
- README roadmap은 slow-light/group-index readout을 미래 항목으로 쓰지만 local
  group index는 이미 출력한다 (`README.md:18-20`;
  `gabes/schemes/absorption.py:771-787`). 미래 항목을 finite-pulse propagation과
  delay-efficiency validation으로 좁혀야 한다.
- Module docstring은 knob가 `Gamma` 단위라고 하지만 UI는 이미 MHz/kHz 기반
  `physical-units-v3`다 (`gabes/schemes/absorption.py:16-18, 462-465`).
- `buffer_ground_relax_khz`의 help는 buffer-gas collisions처럼 보이지만 구현은
  단일 Raman-coherence dephasing이다 (`gabes/schemes/absorption.py:549-553`;
  `gabes/atoms.py:168-190`). `Ground-coherence dephasing`처럼 좁혀야 한다.
- Checklist의 coated-cell 항목은 Lambda 경로를 존재하지 않는
  `schemes/lambda_system.py`로 적는다 (`docs/checklist.json:304-307`). 실제 위치
  `schemes/absorption.py`로 고치는 것은 runtime 0의 provenance 정리다.

### 8.2 테스트

현재 테스트는 다음에 강하다.

- AT split `≈ Ω_c`, power/diameter Rabi scaling
  (`tests/test_absorption.py:110-129`)
- warm angle에서 EIT broadening, cold transparency, sub-natural CPT
  (`tests/test_absorption.py:132-175`)
- unresolved hero/status contract (`tests/test_absorption.py:195-216`)
- Lambda/Rydberg affine kernel과 NumPy reference parity
  (`tests/test_kernels.py:130-165`)
- headless/normal render contract (`tests/test_headless_observables.py:38-89`;
  `tests/test_schemes_render.py:140-177`)

그러나 다음 핵심 회귀는 아직 없다.

- Doppler-off angle invariance
- default FWHM/group-index high-resolution convergence와 interpolated edge
- collinear hyperfine `Delta k`
- natural-Rb scalar/effective-medium status
- transmission, group delay, loss, experimental-visibility의 동시 표시
- `line_strength` navigate-only cache equivalence
- Lambda example의 parameter/code/grid provenance

현재 warm-angle test는 intended effect만 검사하므로 cold operator 오류를 잡지 못한다.
또 CPT test는 sample-count 기반 폭이 자연선폭보다 작다는 넓은 조건만 검사한다
(`tests/test_absorption.py:161-175`). 내부 회귀 통과는 외부 실험 검증을 뜻하지 않는다.

### 8.3 예제

실행 가능한 Lambda 전용 script/config는 없다. 가장 가까운 코드 예제는
`tests/test_absorption.py:110-216`의 direct API 호출이다. User Guide의 `eit.png`와
`at.png`는 현상을 시각적으로 잘 보여 주지만 둘 다 2026-06-08 생성 파일이며 생성
script와 parameter/commit/grid manifest가 없다.

- EIT PNG: `Ω_c = 6.00 MHz`, `gamma_gg = 10.0 kHz`; 현재 기본값 약
  `17.24 MHz`, `57.46 kHz`와 다름
- AT PNG: `Ω_c = 45.97 MHz`, `gamma_gg = 57.5 kHz`; 현재 AT 기본 scale과 가까움
- CPT PNG 또는 실행 가능한 CPT 예제 없음
- 어느 그림에도 resolution status, samples/FWHM, model-scope badge가 없음

Power/diameter, beam angle, dephasing, cell length sweep을 headless CSV/Markdown으로
저장하고 parameter manifest, commit, grid/status, scalar-model badge를 포함하는
Lambda 예제를 추가하면 solver를 바꾸지 않고도 실험 reference 재현성을 크게 높일
수 있다.

## 9. 검증

현재 working tree에서 다음을 실행했다.

```powershell
python -m pytest tests/test_absorption.py tests/test_kernels.py `
  tests/test_headless_observables.py tests/test_schemes_render.py `
  tests/test_docs_consistency.py -q
# 53 passed in 15.32s

python -m pytest -q
# 487 passed in 251.85s
```

추가로 기본 EIT/AT/CPT, uniform 601/9,601점과 center-clustered 601점 수렴,
Doppler-off/on 0/10 mrad, 공용 resolution helper, local group delay,
`line_strength` remap, medium/factory construction benchmark를 직접 계산했다.

## 10. 최종 판단과 권고 순서

Scheme 2는 정상상태 3준위 Lindblad OBE, Maxwell velocity average, lab-facing
Rabi scaling과 Beer–Lambert propagation을 결합한다. 따라서 **EIT/AT/CPT의 기본
scale과 knob sensitivity를 빠르게 보는 반정량적 실험 reference**로는 유용하다.
특히 AT splitting과 coupling power/diameter sweep은 실험 전 sanity check에 적합하다.

그러나 현재 기본 EIT의 hero linewidth와 group index는 grid 미수렴이고,
beam-angle residual Doppler는 cold-limit를 위반한다. 특정 hyperfine/Zeeman,
clock-CPT, buffer/coated-cell, finite-pulse slow-light의 절대 reference로 사용하면 안
된다. `n_g`가 커도 default transmission이 `1.36%`라는 사실을 함께 보아야 한다.

권고 순서는 다음과 같다.

1. `lambda-spectrum-validity` P1을 완료한다: optical scan/velocity operator 분리,
   center-clustered grid, interpolated width, samples/edge/visibility status,
   transmission과 group-delay/`n_g` 동시 표시.
2. `line_strength`를 navigate-only로 옮겨 약 `287×`의 불필요한 재계산을 없앤다.
3. Pure-isotope collinear `Delta k`, natural-Rb scalar 경고와 Raman-dephasing
   명칭을 추가한다.
4. 출처와 uncertainty가 있는 pressure shift/elastic width만 저차 모델에 넣는다.
5. 재현 가능한 Lambda sweep example과 manifest를 만들고 오래된 EIT/AT PNG를
   current status-aware 출력으로 교체한다.
6. Coated-cell ground-coherence transport, full VCC, hyperfine/Zeeman Lambda는
   cell geometry, sourced rates, calibration/holdout spectrum과 runtime budget이
   정해진 뒤 opt-in research/reference mode로 시작한다.

이 순서는 compact 3-level 모델의 유용한 실제 물리를 유지하면서 거의 같은
interactive 비용으로 정밀해 보이는 숫자를 과신할 위험을 먼저 줄인다.
