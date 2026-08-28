# 2026-07-31 Scheme 2 물리 검토 — Lambda coherence

## 1. 오늘의 선택과 현재 다섯 scheme

서울 현지 날짜의 day-of-month는 `31`이므로

```text
n = (31 mod 5) + 1 = 2
```

이다. 현재 dropdown 레지스트리의 실제 순서는 다음과 같다
(`gabes/schemes/__init__.py:12-24`, `README.md:10-16`).

| 순번 | 등록 인스턴스 / 이름 | 표시 scheme | 핵심 물리 |
|---:|---|---|---|
| 1 | `SASScheme()` / `sas` | Absorption spectroscopy (OD / SAS) | pump-off OD, pump-on SAS와 hyperfine pumping |
| 2 | `LambdaScheme()` / `lambda` | Lambda coherence (EIT / AT / CPT) | 축약 3준위 Lambda 결맞음 |
| 3 | `RydbergEITScheme()` / `rydberg_eit` | Rydberg-EIT electrometry | cascade EIT, microwave AT, electrometry |
| 4 | `MagnetoScheme()` / `magneto` | Hanle / EIA / NMOR | Zeeman manifold의 투과·회전 |
| 5 | `FWMScheme()` / `fwm` | FWM (Squeezing / Biphoton) | seeded FWM과 SFWM source estimate |

`eit`, `at`, `cpt` alias는 직접 테스트·호출용이며 dropdown 순서를 바꾸지
않는다 (`gabes/schemes/__init__.py:27-39`). 따라서 오늘 대상은 2번
`LambdaScheme`이다. 검토 기준 HEAD는 `a82bbf7`, 현재 branch는 `main`이다.

## 2. 먼저 찾은 기존 개선안·보고서·TODO·issue note

확인한 범위는 다음과 같다.

- Scheme 2 선행 보고서:
  `docs/daily_report/2026-06-26_scheme-2_lambda-coherence.md`,
  `2026-07-01_scheme-2_lambda-coherence.md`,
  `2026-07-06_scheme-2_lambda-coherence.md`,
  `2026-07-16_scheme-2_lambda-coherence.md`
- 계획/TODO: `docs/checklist.json`, `gabes/constants.py:79-89`
- 구현: `gabes/schemes/absorption.py`, `gabes/atoms.py`,
  `gabes/beam.py`, `gabes/kernels.py`, `gabes/observables.py`,
  `gabes/lineshape.py`
- 테스트: `tests/test_absorption.py:110-216`,
  `tests/test_kernels.py:102-136`, headless/render smoke tests
- 문서·예제: `README.md`, `docs/Userguide/GABES_User_Guide_v2.html`,
  `docs/Userguide/userguide_assets/eit.png`, `at.png`

별도의 로컬 issue/proposal/roadmap 파일은 없고, 저장소 안의 해당 역할은
`docs/checklist.json`과 일일 보고서가 맡고 있다. 기존 개선 제안은 **있다**.

| 기존 항목 | 현재 상태 | 물리 보존과 계산 비용 평가 |
|---|---|---|
| lab-facing `sqrt(P)/d` Rabi | 구현 완료 | solve 전 scalar 변환이다. 차원·solve 수 증가가 없고 실험 노브 대응성이 높다 (`absorption.py:529-548, 583-591`). |
| residual two-photon Doppler / beam angle | 구현 완료 | 기존 velocity average의 coefficient만 바꾸므로 추가 solve가 없다. 다만 아래의 cold-limit 오류는 아직 남았다 (`checklist.json:109-113`). |
| figureless/headless readout | 구현 완료 | 물리를 바꾸지 않고 Matplotlib을 생략한다 (`checklist.json:57-61`; `tests/test_headless_observables.py:38-67`). |
| gas/species/line별 pressure shift·broadening·Dicke proxy | `deferred`, GROUP B | scalar shift와 같은 차원의 dissipator/유효 Doppler 계수로 구현하면 solve 차원·횟수는 그대로다. 검증된 계수표가 전제다 (`checklist.json:50-54`). |
| full velocity-changing collision(VCC) | `deferred`, GROUP C | velocity class를 서로 결합하므로 현재의 분리 가능한 Maxwell 평균을 깨며 무겁다 (`checklist.json:130-134`). |
| hyperfine/Zeeman-resolved Lambda + optical pumping | `deferred`, GROUP C | 실험 충실도는 크게 오르지만 density-matrix 차원과 모든 preset의 의미가 바뀐다 (`checklist.json:137-141`). |

7월 16일 보고서 뒤 Lambda 물리 코어 변경은 없다. 7월 20일
`84aba45`에서 수치적으로 존재하지 않는 EIT/AT feature를 hero 숫자로
올리지 않도록 unresolved status가 추가됐다
(`gabes/schemes/absorption.py:723-795`,
`tests/test_absorption.py:195-216`). 이는 좋은 무부하 readout 개선이지만,
오늘 검토에서 해상도와 실험적 검출 가능성에 관한 두 빈틈을 추가로 확인했다.

## 3. 현재 구현하는 실제 물리

### 3.1 축약 3준위 정상상태 OBE

`atoms.lambda3()`는 `g1, g2, e` 세 상태를 만들고, 들뜬 상태가 두
바닥상태로 `Gamma/2`씩 붕괴하도록 한다. `buffer_ground_relax_khz`는
population reload가 아니라 `rho_g1g2`와 그 켤레의 Raman-coherence
dephasing이다 (`gabes/atoms.py:157-179`, 특히 `175-176`).

Hamiltonian은 고정 약한 probe `Omega_p = 10^-3 Gamma`, coupling
`Omega_c`, 그리고 two-photon detuning을 직접 포함한다
(`gabes/schemes/absorption.py:29-31, 633-656`). 따라서 다음은 실제
정상상태 OBE에서 나온다.

- two-photon resonance의 EIT transparency와 정상분산
- 강한 coupling에서의 Autler-Townes doublet
- 협폭 dark resonance
- coupling power·beam diameter·Raman dephasing에 대한 방향성

다만 CPT는 별도 bichromatic clock-CPT 엔진이 아니다. probe/coupling 비가
기본적으로 약 `10^-3`인 같은 weak-probe Lambda의 좁은 scan preset이다.
“CPT dark resonance”의 개념 확인에는 유효하지만, 균형 잡힌 두 광장,
광펌핑, light shift, clock contrast를 예측하는 reference는 아니다.

### 3.2 scalar D-line medium과 실험 노브

species/line 선택은 온도 밀도, 자연선폭, 파장, 질량, reduced dipole을
바꾼다 (`gabes/schemes/absorption.py:418-443, 521-526`). coupling은
측정된 Rabi anchor를

```text
Omega_c proportional to sqrt(power) / (1/e^2 diameter)
```

로 환산한다 (`gabes/beam.py:9-20`). 이는 절대 polarization·CG·mode overlap을
예측하지는 않지만, 실험에서 보정한 Rabi를 power/diameter sweep으로 옮기는
방식으로는 정직하고 실용적이다.

매질은 특정 `Fg -> Fe` leg, Clebsch-Gordan 계수, Zeeman population을
고르지 않는 scalar model이다. 특히 기본 `Rb (natural)`은 abundance로
분리한 두 isotope의 susceptibility 합이 아니다. `85Rb`를 reference
isotope로 골라 natural-Rb 전체 밀도를 한 동일 Lambda에 넣는다
(`absorption.py:430-442`). 오늘 직접 비교한 기본 EIT에서
`Rb (natural)`과 `85Rb`는 `N`, `omega0`, `chi_bar`가 **bitwise 동일**했다.
실제 natural-Rb 셀의 두 isotope가 같은 Raman resonance에 결맞게 참여한다는
뜻으로 읽으면 안 된다.

### 3.3 Doppler 평균과 Beer-Lambert readout

Doppler-on Lambda는 4-sigma Maxwell grid, `dv = 2 m/s`를 쓰고 각 scan
point에서 velocity-weighted coherence만 수축한다
(`gabes/schemes/absorption.py:118-140`;
`gabes/kernels.py:413-495`). 비영 beam angle은 `|Delta k|/k`를
`AtomModel.S_v`에 넣어 warm-vapor EIT/CPT를 넓히거나 씻어낸다
(`gabes/beam.py:42-55`; `gabes/atoms.py:165-178`).

`chi_bar`는 density와 dipole을 거쳐 physical susceptibility,
`alpha = k Im(chi)`, `T = exp(-alpha L)`로 바뀐다
(`gabes/observables.py:393-418`). 분산의 수치 미분으로 group index도
계산한다 (`gabes/observables.py:426-434`).

## 4. 실험 원자물리 reference로서의 판단

### 유용한 용도

현재 코드는 다음 질문에 답하는 **유용한 real-physics,
semi-quantitative 실험 planning reference**다.

- coupling power 또는 beam diameter를 바꿀 때 EIT에서 AT로 넘어가는 scale
- AT split이 보정된 `Omega_c`를 따라가는지
- ground-coherence `T2`가 EIT/CPT window를 얼마나 흐리는지
- mrad 정렬 오차가 warm-vapor transparency를 얼마나 씻어내는지
- 같은 susceptibility에서 cell length만 바꾼 Beer-Lambert transmission

기본 현재 출력은 다음과 같다. warm-up 후 이 환경의 median compute 시간도
함께 적었다.

| regime | compute median | 현재 hero/readout |
|---|---:|---|
| EIT | `0.108 s` | `T_res = 0.014`, FWHM `0.46 MHz`, `n_g = 1.003e5` |
| AT | `0.024 s` | split `46.0 MHz`, expected `46.0 MHz`, `T_center = 0.989` |
| CPT | `0.019 s` | `T_res = 0.706`, FWHM `923.32 kHz`, `n_g = 9.525e4` |

테스트는 AT split, `sqrt(P)/d` scaling, warm angle broadening, cold EIT,
sub-natural CPT, Numba/NumPy parity, unresolved status를 검증한다
(`tests/test_absorption.py:110-216`; `tests/test_kernels.py:116-136`).

### 믿지 말아야 할 절대값

다음에는 아직 1차 실험 reference가 아니다.

- natural-Rb isotope/hyperfine line assignment와 absolute participating density
- polarization·Zeeman·CG branching·optical pumping에 따른 contrast
- 실제 buffer-cell의 pressure-shifted line centre, homogeneous width,
  diffusion/velocity-changing collision
- clock-CPT의 light shift·균형 광장·절대 linewidth
- 유한 pulse의 delay, distortion, bandwidth와 검출 SNR

따라서 출력 linewidth·contrast·group index를 논문 피팅값으로 직접 쓰기보다,
실험 조건의 **방향성과 order of magnitude**를 잡는 데 쓰는 것이 안전하다.

## 5. 오늘 확인한 핵심 문제

### P0 — 기본 EIT FWHM과 group index가 scan 해상도에 미수렴

`_scan()`은 모든 regime에 601점을 주고
(`gabes/schemes/absorption.py:609-631`), `window_fwhm()`은 half-height
경계를 sample 단위로 걷되 보간하지 않는다
(`gabes/lineshape.py:19-37`). 기본 EIT의 현재 표시 FWHM `0.45968 MHz`는
격자의 정확히 **2칸**이다.

같은 물리 범위에서 점수만 바꾼 진단은 다음과 같다.

| uniform points | step | 현재 sample FWHM | edge-interpolated FWHM | center `n_g` |
|---:|---:|---:|---:|---:|
| 601 | 0.22984 MHz | 0.45968 MHz | 0.22990 MHz | 100,305 |
| 1,201 | 0.11492 MHz | 0.22984 MHz | 0.14169 MHz | 116,233 |
| 2,401 | 0.05746 MHz | 0.22984 MHz | 0.16174 MHz | 113,989 |
| 4,801 | 0.02873 MHz | 0.17238 MHz | 0.15884 MHz | 112,966 |
| 9,601 | 0.01436 MHz | 0.17238 MHz | 0.15910 MHz | 112,693 |

즉 현재 hero FWHM은 수렴값보다 약 `2.9x` 크고, group index는 약 `11 %`
작다. CPT도 현재 `923.32 kHz`가 edge-interpolated 수렴값
약 `803.54 kHz`보다 `15 %` 크다. 현재 테스트는 “sub-natural”이라는
넓은 조건만 보므로 이 문제를 잡지 못한다
(`tests/test_absorption.py:161-175`).

**거의 무부하 개선:**

1. 먼저 `samples_per_width`, `resolution-limited` status와 half-height edge
   interpolation을 추가한다. 추가 OBE solve 없이 `O(N)` 또는 상수 시간이다.
2. EIT/CPT는 같은 601 solve를 중앙에 집중하는 비균일 grid를 쓸 수 있다.
   오늘 `sinh(5u)` 601점 진단은 center step `0.01549 MHz`,
   FWHM `0.15914 MHz`, `n_g = 112,708`을 주어 9,601점 uniform 결과와
   각각 약 `0.03 %`, `0.01 %` 안에서 일치했고 compute 시간은 같은
   약 `0.1-0.2 s` 범위였다. 넓은 wing sampling이 성겨지는 trade-off가
   있으므로 regime별 회귀가 필요하다.
3. publication/reference용 uniform local refinement는 opt-in으로 둔다.
   9,601점 full scan은 약 `1.84 s`라 기본 UI에 그대로 쓰기에는 불필요하게
   무겁다.

### P0 — residual Doppler가 cold spectrum도 바꾸는 기존 오류

기존 7월 16일 보고서의 가장 중요한 지적은 아직 유효하다.
`_affine_scan_coeffs()`가

```text
A_coef = dL/ds - S_v
B_coef = S_v
```

를 만든다 (`gabes/schemes/absorption.py:45-57`). residual ratio가
ground level까지 `S_v`에 들어가므로 (`gabes/atoms.py:165-178`),
velocity가 0이어도 scan coefficient에 `-r*s`가 섞인다. 물리적인
residual Doppler는 `Delta k dot v`이므로 `Doppler=off`에서 angle은
스펙트럼을 바꾸지 않아야 한다.

오늘 current HEAD에서 cold EIT, 0 대 10 mrad 비교 결과:

- `max |Delta chi| / max |chi| = 2.942e-2`
- 최대 transmission 차이 `7.885e-3`

warm-vapor 0/1/5/10 mrad 결과는 이전 보고서와 그대로 재현됐다.

| angle | `|Delta k|/k` | `T_res` | sample window FWHM |
|---:|---:|---:|---:|
| 0 mrad | 0 | 0.750724 | 0.45968 MHz |
| 1 mrad | 0.001000 | 0.326469 | 0.91936 MHz |
| 5 mrad | 0.005000 | 0.159116 | 4.13712 MHz |
| 10 mrad | 0.010000 | 0.138894 | 7.81456 MHz |

**권고:** scan operator와 velocity operator를 분리해 residual ratio는
`kv`에만 곱한다. solver 차원, velocity 수, solve 횟수가 같으므로 runtime
증가는 사실상 없다. `doppler="off"`에서 전체 `chi_bar`의 angle invariance를
새 회귀로 고정해야 한다. 그 전까지 warm angle sweep은 정렬 경향성으로만
사용하고 정량 linewidth에는 해상도 문제까지 함께 표시해야 한다.

### P1 — “수치적으로 존재”와 “실험에서 검출 가능”을 구분하지 않음

7월 20일 readout은 half-height crossing이나 두 local maximum이 없을 때
숫자를 숨긴다 (`absorption.py:723-795`). 그러나 신호의 절대 contrast나
검출 noise floor는 보지 않는다. 허용 범위의
`20 C`, `0.5 mm`, `line_strength=0.01` EIT에서

```text
T_min = 0.9994266
T_max = 0.9999261
contrast = 4.995e-4
```

인데도 hero는 `T_res = 1.000`, `FWHM = 0.92 MHz`라고 표시한다.
이 status는 “모델 배열에 feature가 있음”만 뜻하며 실험적으로 resolved라는
뜻은 아니다.

**권고:** modeled transmission contrast, samples-per-width, 선택 가능한
detector noise/contrast threshold를 readout status에 넣는다. 모두 기존
배열의 `O(N)` 후처리라 추가 solve가 없다. noise floor를 모를 때는
`numerically resolved; experimental visibility not assessed`라고 명시하면
된다.

### P1 — group index만으로 slow-light reference처럼 보임

`group_index()`는 dilute-medium local derivative이며 finite pulse 전파나
distortion을 풀지 않는다 (`gabes/observables.py:426-434`). 수렴 grid로
환산한 기본 EIT는 대략

- `n_g = 1.127e5`
- `L = 15 mm`에서 local group delay 약 `5.64 us`
- 수렴 FWHM과의 delay-bandwidth product 약 `0.90`
- 그러나 resonance transmission은 겨우 `1.36 %`

이다. 큰 `n_g`만 보면 좋은 slow-light 조건처럼 보이지만 대부분의 광이
사라진다. 사용자 가이드는 이를 “느린 빛의 서명”으로 설명한다
(`docs/Userguide/GABES_User_Guide_v2.html:584-586`).

**권고:** 추가 solve 없이 group delay, delay-bandwidth product,
`T_res`와 함께 loss/status를 표시한다. finite-pulse 예측이 아니라 local
linear-response estimate임을 help와 README roadmap에 명시한다.

### P1 — `buffer_ground_relax_khz`의 rate 의미가 너무 넓음

UI help는 이를 “buffer-gas collisions”라고 하지만 실제로는 ground
coherence dephasing 하나다 (`absorption.py:549-553`;
`atoms.py:175-176`). population exchange, transit escape/reload, diffusion,
spin destruction이 분리되지 않는다. 이미 존재하는
`beam.transit_broadening_mhz()`로 50 C, 1 mm Rb의 단순 ballistic scale을
계산하면 약 `40.0 kHz`로, 기본 `57.46 kHz`와 같은 규모다
(`gabes/beam.py:23-28`). 현재 이름만으로 이를 모두 buffer collision로
해석하면 잘못된 실험 추론을 할 수 있다.

**권고:** 우선 이름을 `ground-coherence dephasing / 2pi`로 바꾸고,
Derived table에 optional transit-scale diagnostic을 둔다. 이후 검증된
pressure/diffusion model이 있을 때 residual dephasing, transit/diffusion,
spin-destruction을 scalar rate budget으로 분리한다. 같은 3준위 Liouvillian에
합산되는 rate라 solve 차원과 횟수는 늘지 않는다.

## 6. 기존 물리 개선안의 우선순위와 비용

| 개선 | 실험적 가치 | 예상 비용 | 판단 |
|---|---|---:|---|
| resolution status + edge interpolation | 숫자 과신 방지 | 추가 solve 없음 | 즉시 적용할 P0 |
| EIT/CPT center-clustered scan | 같은 비용에서 수렴 readout | 같은 601 solve | wing 회귀 후 높은 우선순위 |
| cold residual-Doppler coefficient 수정 | 정확성 복원 | 사실상 0 | 즉시 적용할 P0 |
| contrast/detectability status | 실험 검출성과 수치 feature 분리 | `O(N)` 후처리 | 높은 우선순위 |
| natural-Rb scalar semantics 경고 / isotope 강제 | participating density 오해 방지 | 0 | 즉시 가능한 문서·UI 개선 |
| pressure shift + homogeneous width table | buffer-cell line centre·폭 | scalar/table, 같은 차원 | 계수 검증 후 낮은 부하로 적용 |
| phenomenological Dicke proxy | residual-Doppler/width 경향 | 같은 차원 | 적용 범위를 명시한 opt-in |
| full VCC | velocity redistribution | velocity-class coupled solve | heavy mode 필요 |
| full hyperfine/Zeeman Lambda | line assignment·편광·광펌핑 | 훨씬 큰 Liouvillian | 별도 opt-in reference solver |

## 7. 동작을 바꾸지 않는 순수 코드 최적화

### 7.1 `line_strength`를 navigate-only로 이동

`line_strength`는 현재 `recompute=True` 기본값을 상속한다
(`gabes/schemes/absorption.py:557-558`). 그러나 Hamiltonian이나
`chi_bar`에는 들어가지 않고, raw에 저장된 뒤 physical susceptibility를
곱할 때만 사용된다 (`absorption.py:664-668, 680-684`;
`observables.py:393-413`).

기본 EIT의 불필요한 재 solve median은 `0.108 s`, 같은 raw에서
line-strength-only headless remap은 `0.171 ms`였다. 약 `630x` 차이다.

**안전한 변경:** `ParamSpec(..., recompute=False)`로 바꾸고
`observables()`가 live `params["line_strength"]`를 읽게 한다. 여러
line strength에 대해 `chi_bar` bitwise 동일, `alpha` 선형 scaling,
transmission equivalence를 회귀하면 물리 동작을 보존하면서 slider의
heavy solve를 없앨 수 있다.

### 7.2 affine kernel에서 scan-constant matrix를 velocity loop 밖으로 hoist

현재 kernel은 모든 `(scan, velocity)` 쌍에서

```text
base + s*A_coef + kv*B_coef
```

의 9x9 전체 matrix를 다시 조립한다
(`gabes/kernels.py:426-449`). `base + s*A_coef`를 scan당 한 번 만들고,
inner loop에서는 copy와 `kv*B_coef`만 더한 진단은 default EIT
711 velocity classes에서:

| 구현 | median |
|---|---:|
| 현재 | 0.1344 s |
| hoist | 0.1193 s |
| speedup | `1.13x` |

`chi`의 max absolute difference는 정확히 `0`이었다. 기존
Numba/NumPy `<1e-11` parity test와 Doppler on/off 회귀를 유지하면
behavior-preserving 최적화다.

### 7.3 작은 object cache는 후순위

기존 보고서의 `_medium_from_params()`와 `atoms.lambda3()` immutable cache
후보를 다시 측정했다. median construction은 각각 약 `43.5 us`,
`332 us`였다. 기본 EIT solve 약 `0.1 s`에 비해 작으므로 단일 UI
recompute의 우선순위는 낮다. 반복 batch에서는 bounded cache가 유효하지만,
먼저 line-strength cache 분리와 kernel hoist를 적용하는 편이 효과가 크다.

headless readout은 이미 구현되어 있으므로 자동 보고서·parameter sweep은
계속 `headless_observables()`를 써야 한다.

## 8. 문서·테스트·예제의 reference 품질

- README의 한 줄 설명은 현재 기능을 잘 요약하지만 scalar/weak-probe,
  natural-Rb effective-density 한계를 말하지 않는다 (`README.md:13`).
- 사용자 가이드의 limits 표는 full polarization/Zeeman 모델이 아니라고
  정직하게 적는다 (`GABES_User_Guide_v2.html:794-803`). 그러나 CPT가
  weak-probe preset이고 `buffer_ground_relax`가 순수 Raman `T2`이며,
  group index가 loss를 포함한 pulse-delay validation이 아니라는 설명은
  부족하다.
- `eit.png`와 `at.png`는 유용한 시각 예제다. 다만 해당 PNG를 만드는
  Lambda 전용 script/config가 저장소에 없어 현재 코드에서 exact
  reproduction하기 어렵다. 실행 가능한 예제는 사실상
  `tests/test_absorption.py:110-216`의 direct API 호출뿐이다.
- 테스트는 textbook invariant와 smoke contract에는 좋지만,
  cold-angle invariance, FWHM/grid convergence, natural-Rb 대 isotope
  semantics, feature contrast threshold, group-delay/loss status는 고정하지
  않는다.

낮은 비용으로 `analysis/lambda_lab_sweep.py` 같은 headless 예제를 추가해
power/diameter, angle, ground dephasing sweep과 파라미터 provenance를
CSV/Markdown으로 내보내면 실험 reference 재현성이 크게 좋아진다.

## 9. 검증

보고서 작성 과정에서 current HEAD로 다음을 직접 실행했다.

- 현재 EIT/AT/CPT 기본 solve와 headless metrics
- cold 0/10 mrad invariance 진단
- warm 0/1/5/10 mrad sweep
- 301-9,601점 uniform resolution sweep
- 601점 center-clustered grid 진단
- low-contrast allowed-range readout
- natural-Rb 대 `85Rb` bitwise 비교
- line-strength remap 및 object-construction benchmark
- affine inner-loop hoist parity/benchmark

최종 테스트 결과:

```powershell
python -m pytest tests/test_absorption.py tests/test_kernels.py tests/test_headless_observables.py tests/test_schemes_render.py -q
# 52 passed in 50.38 s

python -m pytest -q
# 232 passed in 163.63 s
```

## 10. 최종 결론

Scheme 2는 교과서 그림만 그리는 장난감은 아니다. 실제 3준위 정상상태
OBE, Maxwell velocity average, lab-facing Rabi scaling, residual alignment,
Beer-Lambert propagation을 결합해 **EIT/AT/CPT 실험의 scale과 knob
sensitivity를 빠르게 보는 반정량 reference**로 충분히 유용하다.

그러나 현재 가장 큰 위험은 더 큰 원자 manifold의 부재보다 먼저,
**기본 EIT hero FWHM과 group index가 scan grid에 미수렴**이라는 점이다.
같은 601 solve를 중앙에 배치하고 resolution/status와 보간을 추가하면
runtime을 거의 늘리지 않고 이 문제를 크게 줄일 수 있다. 동시에 기존
cold residual-Doppler coefficient 오류를 같은 차원에서 고쳐야 한다.

그 다음은 modeled contrast와 experimental visibility의 구분, natural-Rb
scalar semantics, ground-dephasing rate budget을 명확히 하는 일이다. full
VCC와 hyperfine/Zeeman solver는 가치가 크지만 interactive 기본 경로에
넣을 저비용 수정이 아니며, 검증 데이터와 별도 heavy mode 설계가 필요하다.
