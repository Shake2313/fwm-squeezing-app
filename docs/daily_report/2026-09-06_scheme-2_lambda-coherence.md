# 2026-09-06 Scheme 2 — Lambda coherence (EIT / AT / CPT) 물리 검토

## 1. 선택 규칙과 현재 다섯 scheme

- 로컬 날짜: `2026-09-06` (Asia/Seoul)
- 계산: `n = (6 mod 5) + 1 = 2`
- 오늘 대상: **Scheme 2, `LambdaScheme` — Λ coherence (EIT / AT / CPT)**
- 검토 시점: `main`, `HEAD d98bdf4`

실행 중인 registry의 순서는 다음과 같다
(`gabes/schemes/__init__.py:12-25`, `README.md:8-16`).

| 번호 | registry key | 표시 이름 | 구현 위치 |
|---:|---|---|---|
| 1 | `sas` | Absorption spectroscopy (OD / SAS) | `gabes/schemes/sas.py` |
| 2 | `lambda` | Λ coherence (EIT / AT / CPT) | `gabes/schemes/absorption.py:466-800` |
| 3 | `rydberg_eit` | Rydberg-EIT electrometry | `gabes/schemes/rydberg.py` |
| 4 | `magneto` | Magneto-optics (Hanle/MOR) | `gabes/schemes/magneto.py` |
| 5 | `fwm` | Four-wave mixing (Squeezing / Biphoton) | `gabes/schemes/fwm.py` |

`eit`, `at`, `cpt` 단일-regime 인스턴스는 테스트와 direct API용 alias이며 dropdown의
별도 scheme이 아니다 (`gabes/schemes/__init__.py:27-34`). `ODScheme`도 Lambda의
2-level 극한을 검사하는 내부 primitive일 뿐 현재 dropdown에는 등록되지 않는다.

작업 트리는 검토 시작 전부터 다수의 staged/unstaged 변경을 포함했다. 공유 파일
`gabes/schemes/absorption.py`와 `gabes/species.py`의 현재 미커밋 변경은 OD/SAS
thermal-sublevel normalization 영역이며 Lambda class의 계산 경로는 현재 `HEAD`와
동일하다. 기존 사용자 변경은 수정하지 않았다.

## 2. 먼저 검색한 기존 제안, TODO, issue note

검색 범위는 `CLAUDE.md`, `README.md`, `docs/checklist.json`, 기존 Scheme 2 일일 보고서
11개, 저장소의 `TODO/FIXME/XXX/issue/proposal` 문자열, issue/TODO 이름의 파일,
그리고 공개 GitHub issue 목록이다. Scheme 2 선행 보고서는 `2026-06-26`, `07-01`,
`07-06`, `07-16`, `07-31`, `08-01`, `08-06`, `08-21`, `08-26`, `08-31`, `09-01`이다.

따라서 **기존 개선안이 없는 상태가 아니다**. 현재 관련 checklist는 다음처럼 잘
정리되어 있다.

| 기존 항목 | 상태 | 이번 판단 |
|---|---|---|
| `lambda-spectrum-validity` | P1, `ready` | residual-Doppler cold-limit 오류와 EIT/CPT grid 신뢰성 문제를 함께 고치는 최우선 항목 (`docs/checklist.json:490-512`) |
| `collisional-coefficient-provenance-and-pressure-shift` | P1, `ready` | 연산은 싸지만 gas/species/line별 계수의 출처·부호·단위·불확도 조사가 큰 일 (`docs/checklist.json:219-253`) |
| `paraffin-coated-cell-extensions` | P2, `parked` | OD/SAS population reservoir를 복사하지 말고 Lambda ground coherence와 bright/dark Ramsey dynamics가 필요 (`docs/checklist.json:315-334`) |
| `full-velocity-changing-collision-kernel` | P2, `parked` | scheme, collision kernel, 검증 spectrum, runtime budget이 정해진 opt-in reference만 타당 (`docs/checklist.json:614-633`) |
| `lambda-hyperfine-resolved-manifold` | P2, `parked` | isotope/transition/polarization과 held-out spectrum을 특정한 뒤 opt-in으로 시작 (`docs/checklist.json:636-653`) |
| `figureless-observables-paths` | `done` | Lambda도 metric/table 계산에서 Matplotlib을 생략한다 (`docs/checklist.json:359-363`; `gabes/schemes/base.py:131-140`) |

코드 TODO는 `gabes/constants.py:81`의 collisional coefficient table 승격 메모 외에
Lambda 전용 신규 제안을 찾지 못했다. 별도의 issue/TODO/proposal 이름 파일도 없다.
공개 [GitHub issue 목록](https://github.com/Shake2313/fwm-squeezing-app/issues?q=is%3Aissue)은
검토 시점에 open/closed 모두 합쳐 0건이었다. 즉 이번 권고는 새 아이디어를 늘리는
것보다 이미 준비된 P1의 실효성과 비용을 다시 확인하는 데 초점을 둔다.

## 3. 9월 1일 검토 이후 달라진 점

이전 Scheme 2 보고서의 기준 `c922ca2`에서 현재 `d98bdf4`로 이동했다. Lambda에
직접 관련된 변화는 다음과 같다.

1. **`line_strength` navigate-only 최적화가 완료됐다.** 파라미터가
   `recompute=False`가 되었고 (`gabes/schemes/absorption.py:561-562`), 저장된
   `chi_bar`를 관측량으로 바꿀 때 live 값을 적용한다
   (`gabes/schemes/absorption.py:681-688`). 회귀도 추가됐다
   (`tests/test_absorption.py:110-132`). 이는 이전 보고서의 가장 큰 순수 코드
   최적화 제안을 그대로 해결한다.
2. 공용 steady-state/affine 경로가 Hermitian-generator real basis에서 계산된다
   (`gabes/core.py:868-901`; `gabes/kernels.py:437-518`). 물리 모델은 바꾸지 않고
   complex density matrix를 real coordinates로 푼다.
3. 공용 `fwhm_halfmax`/`fwhm_interp`의 peak 연결성은 개선됐지만, Lambda가 실제로
   쓰는 것은 여전히 보간하지 않는 `window_fwhm()`이다
   (`gabes/lineshape.py:53-70`; `gabes/schemes/absorption.py:776-790`). 따라서
   `lambda-spectrum-validity`의 linewidth/resolution 요구는 완료되지 않았다.
4. EIT/AT/CPT Hamiltonian, 601점 scan, residual-k 전달, scalar medium, 문서와
   예제 자산에는 물리적 변경이 없다. 아래의 기존 P1 문제도 그대로 재현된다.

## 4. 현재 구현이 실제로 푸는 물리

### 4.1 정상상태 3준위 OBE

`atoms.lambda3()`는 `g1`, `g2`, `e`의 세 준위를 만들고 excited state가 두 ground
state로 각각 `Gamma/2`씩 붕괴하도록 한다. `gamma_gg`는 population exchange가 아니라
`rho_g1,g2`와 그 켤레에 직접 들어가는 phenomenological ground-coherence `T2` rate다
(`gabes/atoms.py:168-190`; 직접 dephasing 적용은 `gabes/atoms.py:90-92`). 실제
hyperfine branching과 polarization-dependent optical pumping은 이 모델에 없다.

Hamiltonian에는 weak probe `Omega_p = 10^-3 Gamma`, coupling field, coupling detuning,
probe/two-photon scan이 직접 들어간다 (`gabes/schemes/absorption.py:29-31, 637-667`).
따라서 EIT transparency, 정상분산, strong-coupling dressed-state AT doublet,
ground-coherence-limited dark resonance는 임의 선모양 합성이 아니라 steady-state
density-matrix 해의 결과다. 이 부분은 유용한 실제 원자물리다.

CPT는 별도의 balanced bichromatic clock-CPT 엔진이 아니라 같은 weak-probe Lambda
엔진의 좁은 scan preset이다 (`gabes/schemes/absorption.py:486-493, 628-635`). 따라서
clock light shift, magnetic sensitivity, pumping transient를 정량 예측하는 모델은 아니다.

### 4.2 매질과 실험 노브

Rb/Cs, D1/D2 선택은 vapor density, mass, wavelength, natural linewidth와
lifetime-derived reduced dipole을 바꾼다 (`gabes/schemes/absorption.py:422-447`;
`gabes/species.py:153-207,239-247`). Coupling power와 beam diameter는 사용자가 정한
Rabi anchor를

```text
Omega_c proportional to sqrt(P) / diameter
```

로 스케일한다 (`gabes/beam.py:9-20`; `gabes/schemes/absorption.py:587-595`). 따라서
power/diameter trend에는 유용하지만 transition dipole, polarization, mode overlap에서
절대 `Omega_c`를 예측하지는 않는다.

`Rb (natural)`은 85Rb와 87Rb susceptibility를 따로 계산해 incoherent sum하지 않는다.
가장 abundance가 큰 isotope의 line constants와 abundance-weighted density/mass를
하나의 scalar Lambda에 넣는다 (`gabes/schemes/absorption.py:434-446`;
`gabes/species.py:173-181`). Pure-isotope 선택도 특정 `Fg -> Fe`, `mF`, polarization을
선택하지 않는다. 따라서 매질 이름은 실제 hyperfine-resolved 실험선보다 **유효 scalar
D-line medium**을 뜻한다.

### 4.3 Doppler 평균과 관측량

Doppler-on 경로는 `dv = 2 m/s`, `+-4 sigma` Maxwell grid에서 scan x velocity마다
9x9 density-matrix Liouvillian을 푼다 (`gabes/schemes/absorption.py:118-140`). 기본
50 C EIT는 601 scan point와 711 velocity classes다. Affine Numba kernel은 scan을
병렬화하고 각 velocity의 probe coherence만 누적한다 (`gabes/kernels.py:437-518`).

이후

```text
chi_phys = -2 N |d|^2 line_strength chi_bar / (epsilon_0 hbar)
alpha = k Im(chi_phys)
T = exp(-alpha L)
```

를 적용한다 (`gabes/observables.py:517-542`). `n_g`는 dilute-medium 근사의
`n = 1 + Re(chi)/2`와 local finite difference로 계산한다
(`gabes/observables.py:550-576`). 이는 CW local dispersion 진단으로는 의미가 있지만,
pulse bandwidth, distortion, delay-efficiency, detector SNR을 포함한 slow-light
propagation 결과는 아니다.

## 5. 실험 원자물리 reference로서의 적합성

현재 Scheme 2는 **실제 물리를 구현한 반정량적 experiment-planning/reference**다.
교과서 곡선을 임의로 그리는 도구보다 낫지만 특정 alkali Lambda 실험의 absolute
spectroscopy나 clock/slow-light reference는 아니다.

| 실험 질문 | 적합성 | 이유 |
|---|---|---|
| `Omega_c` 증가에 따른 EIT -> AT, split sanity check | 좋음 | 3-level OBE와 resolved absorption maxima를 직접 사용 |
| coupling power/diameter 변화의 방향과 scale | 유용 | anchor-preserving `sqrt(P)/d`; 절대 Rabi calibration은 아님 |
| Raman `T2`가 EIT/CPT를 흐리는 경향 | 유용 | 단일 coherence-dephasing rate; transit/diffusion/population exchange와 분리되지 않음 |
| warm-vapor beam misalignment sensitivity | **수정 후 유용** | residual-Doppler 방향은 맞지만 현재 cold-limit invariant를 위반 |
| 특정 isotope/hyperfine/polarization의 contrast/width | 낮음 | scalar 3-level, symmetric branching, Zeeman/pumping 부재 |
| clock-CPT light shift와 자기장 민감도 | 낮음 | balanced two-field clock model과 time dynamics가 없음 |
| slow-light pulse delay/왜곡/효율 | 낮음 | local `n_g`만 있고 finite pulse propagation과 visibility gate가 없음 |
| buffer/coated-cell EIT | 현재 부적합 | collision, diffusion, bright/dark return, Ramsey phase를 단일 `T2`가 대신함 |

## 6. 현재 작업 트리의 직접 재계산

### 6.1 기본 출력과 실행시간

Numba warm-up 뒤 11회 median이다. 시간은 현재 Windows/Python 환경의 상대 비교용이다.

| regime | compute median | 현재 표시값 |
|---|---:|---|
| EIT | `33.882 ms` | `T_res = 0.0135583`, FWHM `0.45968 MHz`, `n_g = 100,305` |
| AT | `3.523 ms` | split `46.0 MHz`, expected `46.0 MHz`, `T_center = 0.989090` |
| CPT | `3.174 ms` | `T_res = 0.706328`, FWHM `923.32 kHz`, `n_g = 95,255` |

AT의 strong-coupling scale은 유용하다. EIT/CPT의 linewidth와 EIT `n_g`는 아래
grid 진단을 통과하지 못한다.

### 6.2 P1 — EIT grid와 hero linewidth/group index

모든 regime은 601점 uniform scan을 쓴다
(`gabes/schemes/absorption.py:613-635`). 현재 `window_fwhm()`은 half-height를 넘은
바깥 sample 사이의 거리를 그대로 반환하며 interpolation이나 samples/FWHM status가
없다 (`gabes/lineshape.py:53-70`).

| EIT scan | 중앙 step | 현재 sample width | 보간 FWHM | samples/FWHM | center `n_g` | solve median |
|---|---:|---:|---:|---:|---:|---:|
| uniform 601 | `0.229840 MHz` | `0.459680 MHz` | `0.229899 MHz` | `1.0003` | `100,305` | `33.182 ms` |
| `sinh(5u)` clustered 601 | `0.015488 MHz` | `0.186156 MHz` | `0.159141 MHz` | `10.27` | `112,708` | `33.420 ms` |
| uniform 9,601 | `0.014365 MHz` | `0.172380 MHz` | `0.159102 MHz` | `11.08` | `112,693` | `476.933 ms` |

현재 hero FWHM `0.45968 MHz`는 converged/interpolated `0.15910 MHz`보다 약 `2.89x`
크고, `n_g`는 약 `11.0%` 작다. 같은 601점을 중앙에 모으면 9,601점 기준과 FWHM
`0.025%`, `n_g` `0.013%` 이내로 일치했고 실행시간 증가는 이번 median에서 `0.7%`였다.
이는 checklist의 `<=20%` budget 안에서 동일 OBE와 동일 점 수를 유지하는 개선이다.

### 6.3 CPT도 현재 resolution gate를 통과하지 못함

기본 CPT의 중앙 step은 `0.230829 MHz`다.

| CPT scan | 현재 sample width | 보간 FWHM | center `n_g` |
|---|---:|---:|---:|
| uniform 601 | `0.923317 MHz` | `0.805733 MHz` | `95,255.0` |
| `sinh(5u)` clustered 601 | `0.834352 MHz` | `0.803540 MHz` | `95,272.0` |
| uniform 9,601 | `0.807903 MHz` | `0.803539 MHz` | `95,272.0` |

표시 width는 converged/interpolated 값보다 약 `14.9%` 크며 uniform 601의
samples/FWHM은 약 `3.49`다. EIT보다는 덜 심하지만 checklist가 요구하는 최소 8을
충족하지 않는다. 현재 `width_valid`는 finite positive width만 확인하므로 이 값을
여전히 hero로 승격한다 (`gabes/schemes/absorption.py:781-799`).

### 6.4 P1 — residual Doppler가 cold spectrum도 바꿈

현재 affine 조립은

```text
A_coef = dL/ds - S_v
B_coef = S_v
L(s, kv) = base + s A_coef + kv B_coef
```

이다 (`gabes/schemes/absorption.py:45-57`; `gabes/kernels.py:477-505`). Angle-dependent
ground-level residual coefficient가 optical scan shift에도 섞여 `v = 0`인 Doppler-off
계산까지 바꾼다. 실제 residual `Delta k dot v`는 `v = 0`에서 정확히 0이어야 한다.

오늘 0 대 10 mrad 비교를 재현하면 다음과 같다.

- Doppler off: `max |Delta chi_bar| / max |chi_bar| = 2.942038%`
- 10 mrad에서 `|Delta k|/k = 0.009999958`
- Doppler on, 50 C, 3 mm: resonance transmission `0.423098 -> 0.00267946`

따라서 warm-vapor alignment effect 자체는 크고 유용하지만 현재 operator bookkeeping은
cold invariant를 깬다. Optical scan operator와 velocity-only residual operator를
분리하면 matrix dimension, scan points, velocity classes, solve 횟수가 모두 그대로다.

### 6.5 collinear frequency difference와 slow-light 해석

`_two_photon_k_ratio()`는 probe와 coupling에 같은 `medium["k_vec"]`를 넘겨 beam
angle 0에서 residual을 정확히 0으로 만든다
(`gabes/schemes/absorption.py:598-601`; `gabes/beam.py:42-55`). 실제 ground-state
splitting 때문에 co-propagating optical fields도 작은 longitudinal `Delta k`를 갖는다.
이전 검토의 analytic thermal Gaussian estimate는 25 C에서 85Rb `4.07 kHz`,
87Rb `9.07 kHz`, 133Cs `9.86 kHz`였다
(`docs/daily_report/2026-08-06_scheme-2_lambda-coherence.md:210-233`). Pure-isotope
mode부터 실제 두 wavevector를 쓰는 것은 scalar coefficient 변경일 뿐 추가 solve가 없다.

기본 EIT의 표시 `n_g = 100,305`를 15 mm cell delay로 단순 환산하면 약 `5.02 us`지만
같은 점의 transmission은 `1.36%`뿐이고 linewidth도 미수렴이다. 따라서 이 값은 local
CW dispersion diagnostic이다. usable pulse delay 예측으로 해석할 수 없다.

## 7. 기존 물리 개선안의 계산비용과 물리 보존성

| 개선 | 물리 보존성 | 계산비용 | 권고 |
|---|---|---:|---|
| scan/velocity operator 분리 | 잘못 섞인 Doppler 항만 분리, 3-level 모델 동일 | 사실상 0; 추가 solve 없음 | 즉시 P1 |
| center-clustered 601 scan | 동일 OBE를 feature에 적합한 grid에서 평가 | 이번 EIT `+0.7%` | off-center detuning/AT wing 회귀 후 P1 |
| interpolated width + samples/edge/contrast status | 계산된 배열의 해석만 개선 | `O(N)`, 추가 solve 없음 | P1과 함께 |
| pure-isotope collinear `Delta k` | 실제 probe-coupling wavevector bookkeeping | scalar coefficient, 추가 solve 없음 | operator 분리 후 P2 |
| `n_g`, transmission, group delay, bandwidth gate 동시 표시 | 새 원자모델 없이 과신만 방지 | `O(N)`, 추가 solve 없음 | P1과 함께 |
| sourced pressure shift/elastic optical width | line center와 homogeneous width 개선 | `O(1)` lookup/arithmetic, solve dimension 동일 | provenance 확보 후 P1 |
| natural-Rb isotope-resolved incoherent sum | scalar mixture보다 실제 isotope 구조에 가까움 | isotope별 solve로 대략 2배 | opt-in 또는 pure isotope 우선 |
| coated-cell bright/dark coherence | population-only 복사를 피하고 Ramsey physics 보존 | dark-time/trajectory quadrature로 중간~높음 | geometry/T2/holdout data 전에는 parked |
| full VCC/diffusion | velocity redistribution과 Dicke/Ramsey physics | velocity classes 결합, separability 상실 | 한 scheme/한 dataset의 research solver |
| hyperfine/Zeeman Lambda | line assignment, polarization, pumping 추가 | Hilbert/Liouville dimension 급증 | 특정 transition의 opt-in reference |

Pressure shift와 elastic width는 runtime이 아니라 provenance가 병목이다. 출처 없는
빠른 phenomenological knob는 experimental reference 신뢰도를 오히려 낮춘다. Coated
cell도 단일 exponential `T2`만 추가하면 싸지만, bright/dark return과 Ramsey phase를
재현하지 못하므로 물리 보존 개선이라고 보기 어렵다.

## 8. 동작을 바꾸지 않는 순수 코드 최적화

### 8.1 완료 — `line_strength` navigate-only

현재 `line_strength`는 heavy cache key에서 빠졌고 같은 raw `chi_bar`를 live remap한다.
오늘 확인 결과 `line_strength = 1.0`과 `1.5`로 다시 `compute()`한 `chi_bar`는
array-exact했다.

| 작업 | median |
|---|---:|
| 불필요한 EIT full compute를 가정한 비용 | `33.409 ms` |
| 같은 raw의 `headless_observables` remap | `0.1179 ms` |
| 비용비 | 약 `283x` |

따라서 이전 보고서의 가장 큰 behavior-preserving 최적화는 완료됐다. 이 항목을 다시
제안할 필요는 없다.

### 8.2 남은 후보 — immutable construction cache는 낮은 우선순위

기본 EIT에서 현재 median은 다음과 같다.

| construction | median |
|---|---:|
| `_medium_from_params()` | `18.8 us` |
| `atoms.lambda3()` | `182.0 us` |
| `velocity_grid()` | `17.9 us` |
| 합계 | `218.7 us` |

이는 EIT full compute의 약 `0.65%`다. Exact parameter tuple로 작은 LRU를 만들고
returned arrays를 read-only로 유지하면 behavior-preserving이지만 우선순위는 낮다.
Doppler-off AT/CPT에서는 상대 비중이 더 크므로 cold scan 반복 workload가 실제로
병목일 때만 재측정할 가치가 있다.

### 8.3 현재 권고하지 않는 후보

- 이전 affine inner-loop hoist prototype은 현재보다 약 `1.5%` 느렸다
  (`docs/daily_report/2026-08-26_scheme-2_lambda-coherence.md:297-309`). LU solve가
  지배하므로 여러 CPU에서 일관된 이득이 나오기 전에는 적용하지 않는다.
- Figureless path와 real-basis solve는 이미 구현됐다. 동일 기능을 다른 캐시나 별도
  solver로 중복시키지 않는다.
- 9,601점 uniform scan은 정확하지만 이번 EIT에서 `476.9 ms`, 즉 601점의 약 `14.4x`다.
  같은 정확도를 내는 clustered 601이 있으므로 interactive default 최적화가 아니다.

## 9. 문서, 테스트, 예제의 reference 품질

### 9.1 문서

- README의 Lambda 설명은 기능을 맞게 요약하지만 scalar/weak-probe/Rabi-anchor/
  natural-Rb effective-medium 한계를 밝히지 않는다 (`README.md:13`).
- 모듈 docstring은 knob가 `Gamma` 단위라고 쓰지만 UI는 이미 MHz/kHz의
  `physical-units-v3`다 (`gabes/schemes/absorption.py:16-18, 468-469`).
- `buffer_ground_relax_khz`는 이름과 help가 buffer-gas collision처럼 보이지만
  구현은 직접 지정한 ground-coherence dephasing rate다
  (`gabes/schemes/absorption.py:553-557`; `gabes/atoms.py:168-190`).
- User Guide는 EIT dispersion을 “느린 빛의 서명”으로 소개하면서 현재의 낮은
  transmission, unresolved grid, finite-pulse 한계를 붙이지 않는다
  (`docs/Userguide/GABES_User_Guide_v2.html:595-606`). Lambda 한계 표도 full
  polarization/Zeeman이 아니라는 한 문장뿐이다 (`...:849-851`).
- User Guide는 각 scheme의 Reference/about panel에 1차 문헌이 있다고 말하지만
  (`...:915-925`), `LambdaScheme`은 `info()`를 override하지 않는다
  (`gabes/schemes/base.py:116-117`; Lambda method 목록 `466-800`). Lambda 전용
  source/provenance panel은 사실상 없다.
- Checklist의 coated-cell 항목은 존재하지 않는 `schemes/lambda_system.py`를 가리킨다
  (`docs/checklist.json:315-317`). 실제 파일은 `schemes/absorption.py`다.

### 9.2 테스트

현재 테스트가 잘 고정하는 것은 다음과 같다.

- AT split `~ Omega_c`, power/diameter Rabi scaling
  (`tests/test_absorption.py:135-154`)
- warm 10 mrad EIT broadening, cold transparency, sub-natural CPT
  (`tests/test_absorption.py:157-200`)
- line-strength live remap과 unresolved extreme/weak-AT status
  (`tests/test_absorption.py:110-132,220-241`)
- Lambda/Rydberg affine kernel과 NumPy reference parity
  (`tests/test_kernels.py:102-136`)
- headless와 정상 figure render contract
  (`tests/test_headless_observables.py:38-89`; `tests/test_schemes_render.py:97-107`)

그러나 다음 핵심 회귀는 없다.

- Doppler-off angle invariance
- EIT/CPT default FWHM와 `n_g`의 high-resolution convergence
- Lambda-specific interpolated edge, samples/FWHM, scan-edge, contrast gate
- co-propagating hyperfine `Delta k`
- natural-Rb effective-medium provenance/status
- transmission/loss/pulse bandwidth와 group-delay 해석의 동시 gate
- Lambda parameter/grid/commit manifest를 가진 standalone example

Warm-angle test는 intended warm effect만 확인하므로 cold operator 오류를 잡지 못한다.
Extreme-opacity status test도 half-height crossing 자체가 없을 때만 막으며, 기본 EIT의
1 sample/FWHM이나 CPT의 3.49 samples/FWHM을 막지 못한다.

### 9.3 예제와 그림

Lambda 전용 실행 script/config/CSV example은 없다. 가장 가까운 실행 예제는
`tests/test_absorption.py:135-241`의 direct API 호출이다. User Guide에는 EIT/AT PNG만
있고 CPT 그림은 없다.

- `docs/Userguide/userguide_assets/eit.png`: 2026-06-08, `Omega_c = 6.00 MHz`,
  `gamma_gg = 10.0 kHz`; 현재 EIT default와 다름
- `docs/Userguide/userguide_assets/at.png`: 2026-06-08, `Omega_c = 45.97 MHz`,
  `gamma_gg = 57.5 kHz`; 현재 AT scale과는 가깝지만 생성 manifest가 없음
- 두 그림 모두 source script, full parameter set, commit, grid, resolution status,
  scalar-model badge가 없다.

`coupling power/diameter`, beam angle, dephasing, cell length sweep을 headless CSV/Markdown로
내보내고 parameter JSON, commit, grid/status, scalar-model badge를 묶는 예제가 있으면
solver를 바꾸지 않고 experimental reference 재현성을 크게 높일 수 있다.

## 10. 검증

현재 작업 트리에서 실행한 결과다.

```powershell
python -m pytest -q tests/test_absorption.py tests/test_kernels.py `
  tests/test_headless_observables.py tests/test_schemes_render.py `
  tests/test_docs_consistency.py
# 54 passed in 15.07s

python -m pytest -q tests analysis
# 595 passed in 276.85s

python -m pytest -q
# collection error: references/AutoOD가 undeclared optional dependency sympy를 요구
```

루트 전체 수집 중단은 제품 Lambda 실패가 아니다.
`references/AutoOD/test_autood_gui_adjustments.py`가
`sympy.physics.wigner`를 import하지만 현재 `requirements.txt` 환경에 `sympy`가 없어서
collection이 멈췄다. 제품/분석 범위 595개는 모두 통과했다.

추가로 EIT/AT/CPT defaults, EIT/CPT uniform 601/9,601 및 `sinh(5u)` clustered 601,
Doppler-off/on 0/10 mrad, `line_strength` cache equivalence, construction overhead를
현재 작업 트리에서 직접 계산했다.

## 11. 최종 판단과 권고 순서

Scheme 2는 정상상태 3-level OBE, Maxwell average, anchored Rabi scaling,
Beer-Lambert propagation을 결합하므로 **EIT/AT/CPT의 line-shape와 knob sensitivity를
빠르게 보는 반정량적 실험 reference**로 유용하다. 특히 resolved AT splitting과
coupling power/diameter sweep은 실험 전 sanity check에 적합하다.

그러나 기본 EIT hero linewidth는 수렴값보다 `2.89x` 크고, CPT width도 약 `14.9%`
크다. Beam-angle residual Doppler는 cold-limit를 위반한다. Specific hyperfine/Zeeman,
clock-CPT, buffer/coated-cell, finite-pulse slow-light의 absolute reference로 사용하면 안
된다. 큰 `n_g`는 낮은 transmission과 grid status를 반드시 함께 보아야 한다.

권고 순서는 다음과 같다.

1. `lambda-spectrum-validity` P1을 완료한다: scan/velocity operator 분리,
   Lambda-specific clustered scan, interpolated width, samples/edge/contrast status,
   transmission과 group-delay/`n_g` 동시 표시.
2. Pure-isotope co-propagating `Delta k`, natural-Rb scalar/effective-medium 경고,
   Raman-dephasing 명칭을 정리한다.
3. 출처와 uncertainty가 있는 pressure shift/elastic width만 low-order model에 넣는다.
4. Reproducible Lambda sweep example과 manifest를 만들고 6월 EIT/AT 그림을 current
   status-aware 출력으로 교체하며 CPT 예제를 추가한다.
5. Construction cache는 실제 cold repeated-scan profile에서 이득이 확인될 때만 한다.
6. Coated-cell coherence transport, full VCC, hyperfine/Zeeman Lambda는 cell geometry,
   sourced rates, calibration/held-out spectrum, runtime budget이 확보된 뒤 opt-in
   research/reference mode로 시작한다.

이 순서는 compact 3-level 모델의 실제 물리를 보존하면서 거의 같은 interactive
비용으로 정밀해 보이는 숫자를 과신할 위험부터 줄인다.
