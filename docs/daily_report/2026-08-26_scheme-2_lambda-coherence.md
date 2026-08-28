# 2026-08-26 Scheme 2 물리 검토 — Lambda coherence (EIT / AT / CPT)

## 1. 오늘의 선택과 현재 다섯 scheme

서울 현지 날짜의 day-of-month는 `26`이므로

```text
n = (day mod 5) + 1 = (26 mod 5) + 1 = 2
```

이다. UI 드롭다운 순서는 `_SCHEMES`가 정한다
(`gabes/schemes/__init__.py:12-24`). 현재 정의와 순서는 다음과 같다.

| 순번 | 등록 인스턴스 | UI scheme | 주된 물리 |
|---:|---|---|---|
| 1 | `SASScheme()` | OD / SAS | Doppler 흡수, 포화흡수, hyperfine pumping |
| 2 | `LambdaScheme()` | Λ coherence (EIT / AT / CPT) | 축약 3준위 Λ 결맞음 |
| 3 | `RydbergEITScheme()` | Rydberg-EIT electrometry | cascade EIT, microwave AT, electrometry |
| 4 | `MagnetoScheme()` | Hanle / EIA / NMOR | Zeeman OBE의 투과·회전 |
| 5 | `FWMScheme()` | FWM | seeded mean-field FWM, SFWM biphoton |

`eit`, `at`, `cpt`는 테스트와 직접 호출용 alias이며 드롭다운 항목을 늘리지
않는다 (`gabes/schemes/__init__.py:27-39`). README의 현재 scheme 표도 같은 다섯
엔진을 기술한다 (`README.md:8-16`). 따라서 오늘 대상은 2번 `LambdaScheme()`다.

검토 시점은 branch `main`, HEAD
`c964a724b1ad370e53af1a89c78f8df04fd37983`이다. 작업 트리에는 기존 사용자
변경이 많이 있어 보존했다. `gabes/schemes/absorption.py`에는 HEAD 대비 diff가
없다. `gabes/atoms.py`의 일반 `collapse_ops` 지원과 `gabes/kernels.py`의 FWM
Floquet 변경은 현재 Lambda factory와 affine-scan 구간을 바꾸지 않는다. 즉
2026-08-21 Scheme 2 검토 뒤 Lambda 동작 자체의 변경은 발견하지 못했다.

## 2. 먼저 검색한 기존 제안, TODO, issue note

다음 자료를 먼저 확인했다.

- Scheme 2 선행 보고서 8개: `2026-06-26`, `07-01`, `07-06`, `07-16`,
  `07-31`, `08-01`, `08-06`, `08-21`
  (`docs/daily_report/*scheme-2_lambda-coherence.md`)
- 통합 계획과 TODO: `docs/checklist.json`, `gabes/constants.py:79-89`
- 코드: `gabes/schemes/absorption.py`, `gabes/atoms.py`, `gabes/beam.py`,
  `gabes/kernels.py`, `gabes/observables.py`, `gabes/lineshape.py`
- 테스트: `tests/test_absorption.py`, `tests/test_kernels.py`,
  `tests/test_headless_observables.py`, `tests/test_schemes_render.py`
- 문서·예제: `README.md`, `CLAUDE.md`,
  `docs/Userguide/GABES_User_Guide_v2.html`, `eit.png`, `at.png`, `analysis/`
- 로컬 파일명과 본문의 `TODO`, `issue`, `proposal` 검색 및 공개 GitHub Issues

별도의 Lambda TODO/issue/proposal 파일과 실행 가능한 Lambda 전용 예제는 없다.
공개 GitHub issue도 현재 0건이다. 그러나 기존 개선 제안은 통합 checklist에
명확히 정리돼 있다.

| 기존 항목 | 상태 | 현재 판단 |
|---|---|---|
| `lambda-spectrum-validity` | P1, `ready` | residual Doppler의 cold-limit 오류, grid 미수렴, visibility/group-delay 해석을 함께 고치는 핵심 항목 (`docs/checklist.json:439-461`) |
| `collisional-coefficient-provenance-and-pressure-shift` | P1, `ready` | gas/species/line별 출처·단위·부호·불확도를 먼저 정하고 scalar shift/width를 적용 (`docs/checklist.json:209-243`) |
| `full-velocity-changing-collision-kernel` | P2, `parked` | 대상 scheme·kernel·검증 spectrum·runtime budget이 정해진 뒤 opt-in reference solver로 시작 (`docs/checklist.json:563-582`) |
| `lambda-hyperfine-resolved-manifold` | P2, `parked` | 특정 isotope/transition/polarization과 holdout spectrum을 고른 뒤 opt-in으로 구현 (`docs/checklist.json:585-603`) |

따라서 이번에는 새 물리 기능을 임의로 제안하기보다 위 제안이 현재 코드에서
여전히 필요한지, 물리를 보존하면서 계산비용을 거의 늘리지 않고 구현할 수
있는지를 재검증했다.

## 3. 현재 구현이 실제로 푸는 물리

### 3.1 정상상태 3준위 Lindblad OBE

`atoms.lambda3()`는 `g1`, `g2`, `e`의 세 상태를 만들고 들뜬 상태가 두
바닥상태로 각각 `Gamma/2`씩 붕괴하도록 한다. `gamma_gg`는 population reload가
아니라 `rho_g1,g2`와 그 켤레를 직접 감쇠시키는 Raman-coherence `T2` rate다
(`gabes/atoms.py:168-190`).

Hamiltonian은 약한 probe `Omega_p = 10^-3 Gamma`, coupling field, coupling
detuning과 two-photon scan을 직접 포함한다
(`gabes/schemes/absorption.py:29-31, 633-656`). 따라서 다음 현상은 임의의
선모양 합성이 아니라 정상상태 density matrix에서 나온다.

- two-photon destructive interference의 EIT transparency와 정상분산
- 강한 coupling이 만드는 Autler-Townes dressed doublet
- Raman coherence dephasing에 의해 제한되는 dark resonance
- coupling power와 beam diameter에 따른 `Omega_c proportional to sqrt(P)/d`
  scaling (`gabes/schemes/absorption.py:583-591`; `gabes/beam.py:9-20`)

이는 실제 유용한 원자물리다. 상온 Rb D1의 weak-probe/strong-coupling Lambda
EIT에서 단순 Doppler 포함 이론이 실험과 정성적으로 맞는다는 고전적 근거도 있다
([Li & Xiao, PRA 51, R2703 (1995)](https://doi.org/10.1103/PhysRevA.51.R2703)).

단, `CPT`는 균형 잡힌 두 optical field를 쓰는 별도 clock-CPT 엔진이 아니라 같은
weak-probe Lambda 모델의 좁은 scan preset이다
(`gabes/schemes/absorption.py:482-500, 624-631`). light shift, magnetic
sensitivity, 실제 optical pumping을 예측하는 clock reference로 읽으면 안 된다.

### 3.2 scalar D-line medium

Rb/Cs와 D1/D2 선택은 온도별 밀도, 자연선폭, 질량, 파장과 reduced dipole을
바꾼다 (`gabes/schemes/absorption.py:418-443, 521-526`). 그러나 구체적인
`Fg -> Fe`, CG coefficient, Zeeman population, polarization selection과 hyperfine
optical pumping은 없다. 특히 `Rb (natural)`은 85Rb와 87Rb susceptibility를
비결맞음 합산하지 않고, abundance가 가장 큰 isotope의 line constant와
abundance-weighted density/mass를 한 scalar Lambda에 넣는다
(`gabes/schemes/absorption.py:430-442`).

따라서 이 경로는 특정 isotope/hyperfine transition의 절대 contrast와 linewidth
reference가 아니라 **유효 scalar medium**이다. 사용자 가이드는 이 한계를
정확히 한 줄로 밝힌다 (`docs/Userguide/GABES_User_Guide_v2.html:828-830`).

### 3.3 Maxwell 평균, Beer-Lambert, group index

Doppler-on Lambda는 4-sigma Maxwell grid를 `dv=2 m/s`로 만들고 각
scan x velocity의 9x9 real-Hermitian Liouvillian을 푼 뒤 probe coherence만
속도 가중 수축한다 (`gabes/schemes/absorption.py:118-140`;
`gabes/kernels.py:436-518`). 기본 50 degC EIT는 601 scan point와 711 velocity
class를 쓴다.

그 뒤 `chi_bar`는 밀도와 dipole을 포함한 `chi_phys`,
`alpha = k Im(chi_phys)`, `T = exp(-alpha L)`로 바뀐다
(`gabes/observables.py:520-542`). group index는 dilute-medium 근사의
`Re(chi)` 국소 수치 미분이다 (`gabes/observables.py:550-558`). CW susceptibility와
local dispersion으로서는 올바르지만 finite-pulse delay, distortion, bandwidth,
검출 SNR을 푸는 slow-light propagation model은 아니다. EIT의 CW 응답과 pulse
propagation이 별도 단계라는 점은 표준 review에도 명확하다
([Fleischhauer et al., RMP 77, 633 (2005)](https://doi.org/10.1103/RevModPhys.77.633)).

## 4. 실험 원자물리 reference로서의 적합성

현재 Scheme 2는 **실제 물리를 구현한 반정량적 실험 planning/reference**다.
교과서용 임의 곡선보다 유용하지만, 특정 alkali Lambda 실험의 절대
spectroscopy·clock·slow-light reference는 아니다.

| 실험 질문 | 적합성 | 이유 |
|---|---|---|
| coupling power/diameter 변화와 EIT-to-AT scale | 좋음 | 정상상태 OBE와 보정된 Rabi scaling |
| AT split 대 측정 `Omega_c` sanity check | 좋음 | 두 local absorption maximum과 `split approximately Omega_c`를 직접 검사 |
| Raman `T2` 증가가 EIT/CPT를 흐리는 방향 | 유용 | 실제 coherence dephasing이나 population exchange·transit·diffusion과 분리되지 않음 |
| warm-vapor beam alignment sensitivity | 오류 수정 후 유용 | residual Doppler 방향은 맞지만 현재 cold-limit 위반 |
| 특정 isotope/hyperfine/polarization의 절대 contrast·width | 낮음 | scalar 3-level, 대칭 branching, Zeeman/optical pumping 부재 |
| clock-CPT light shift·자기장 민감도 | 낮음 | balanced bichromatic CPT가 아니라 weak-probe preset |
| slow-light pulse delay·왜곡·효율 | 낮음 | local `n_g`만 있고 pulse propagation과 detector model이 없음 |

온증기 EIT의 VCC, transit/influx, Raman dephasing은 서로 대체할 수 없는 채널이다.
이를 분리해 실험과 비교한 근거는
[Ghosh et al., arXiv:0901.3790](https://arxiv.org/abs/0901.3790)에 있다. 그러므로
현재 `buffer_ground_relax_khz` 하나를 모든 충돌·이동 물리로 해석해서는 안 된다.

## 5. 2026-08-26 현재 작업 트리 직접 재계산

### 5.1 기본 출력

Numba warm-up 뒤 각 기본값을 5회 계산한 median이다. 절대 시간은 현재 환경에
한정된다.

| regime | compute median | 현재 표시 metric |
|---|---:|---|
| EIT | `33.666 ms` | `T_res=0.0135583`, FWHM `0.45968 MHz`, `n_g=100,305` |
| AT | `3.281 ms` | split `46.0 MHz`, expected `46.0 MHz`, `T_center=0.989090` |
| CPT | `3.218 ms` | `T_res=0.706328`, FWHM `923.32 kHz`, `n_g=95,255` |

AT의 strong-coupling scale은 유용하다. EIT/CPT width와 EIT `n_g`는 다음 grid
문제 때문에 정밀 숫자로 쓰면 안 된다.

### 5.2 P1 — 기본 EIT linewidth와 group index 미수렴

모든 regime은 601점 uniform scan을 쓴다
(`gabes/schemes/absorption.py:609-631`). `window_fwhm()`은 half-height edge를
보간하지 않고 바깥 sample 차이를 반환한다 (`gabes/lineshape.py:53-71`). 기본
EIT 표시 FWHM은 중앙 grid 두 칸뿐이다.

| scan | 중앙 step | 표시/sample FWHM | edge-interpolated FWHM | center `n_g` | warm solve |
|---|---:|---:|---:|---:|---:|
| uniform 601 | `0.229840 MHz` | `0.459680 MHz` | `0.229899 MHz` | `100,305` | `33.666 ms` |
| center-clustered 601 | `0.015488 MHz` | `0.186156 MHz` | `0.159141 MHz` | `112,708` | `33.267 ms` |
| uniform 9,601 | `0.014365 MHz` | `0.172380 MHz` | `0.159102 MHz` | `112,693` | `470.504 ms` |

현재 표시 FWHM은 수렴 보간값보다 약 `2.89x` 크고, `n_g`는 약 `11%` 작다.
동일한 601점을 중심에 집중하면 FWHM은 9,601점 기준과 `0.025%`, `n_g`는
`0.014%` 안에서 맞고 이번 측정에서는 오히려 `1.2%` 빨랐다. 즉 **동일한 OBE와
동일한 solve 수를 유지하며 사실상 무비용으로 수치 물리를 보존할 수 있다**.

공용 `subdoppler_feature()`는 보간 edge, samples/FWHM, scan-edge clearance와
status를 이미 계산한다 (`gabes/lineshape.py:141-153`). 현재 기본 EIT에 적용하면
`resolution-limited`, `1.0003 samples/FWHM`이며 후처리 median은 `0.226 ms`다.
Lambda의 전체-window floor 정의는 SAS의 running-median background와 다르므로
edge interpolation/status 구조만 공유하는 것이 안전하다. checklist의 요구치는
numeric hero에 최소 8 samples/FWHM이다 (`docs/checklist.json:452-459`).

### 5.3 P1 — residual Doppler의 cold-limit 위반

현재 `_affine_scan_coeffs()`는

```text
A_coef = dL/ds - S_v
B_coef = S_v
L(s, kv) = base + s A_coef + kv B_coef
```

로 조립한다 (`gabes/schemes/absorption.py:45-57`; `gabes/kernels.py:477-490`).
angle-dependent ground coefficient가 scan 항에도 들어가므로 velocity가 0인
Doppler-off 계산에서도 angle이 spectrum을 바꾼다. 실제 residual Doppler
`Delta k dot v`는 `v=0`에서 정확히 사라져야 한다.

현재 cold EIT에서 0 대 10 mrad를 비교하면

- `max |Delta chi_bar| / max |chi_bar| = 2.94204%`
- `|Delta k|/k = 0.00999996`

로, 2026-07-16 이후 보고된 오류가 그대로다. 반면 warm 50 degC, 3 mm에서는
10 mrad가 공명 투과를 `0.423098 -> 0.002679`로 낮추므로 intended alignment
physics 자체는 강하고 유용하다.

optical scan-shift operator와 velocity-only two-photon residual operator를
분리하면 matrix 차원, scan 수, velocity class 수가 그대로다. 따라서 추가 solve
없이 정확성만 회복할 수 있다. `doppler="off"`에서 angle 0/10 mrad의 전체 complex
`chi_bar` invariance를 `rtol=1e-12`, `atol=1e-14`로 고정해야 한다
(`docs/checklist.json:453-455`).

### 5.4 P2 — collinear optical frequency difference 누락

`_two_photon_k_ratio()`는 probe와 coupling에 같은 `medium["k_vec"]`를 넘겨
angle 0에서 residual을 정확히 0으로 만든다
(`gabes/schemes/absorption.py:593-597`; `gabes/beam.py:36-55`). 실제 hyperfine
Lambda의 두 광장은 ground splitting만큼 주파수가 다르므로 co-propagating이어도
작은 종방향 `Delta k`가 남는다. 선행 2026-08-06 보고서가 current species
constant로 계산한 25 degC Gaussian FWHM은 85Rb `4.07 kHz`, 87Rb `9.07 kHz`,
133Cs `9.86 kHz`였다.

기본 EIT에는 작지만 협폭 CPT에는 무시하기 어렵다. pure-isotope mode에서 실제
두 wavevector를 만들면 기존 velocity coefficient의 scalar만 바뀌어 추가 solve가
없다. 다만 앞의 cold-limit operator 분리가 먼저다.

### 5.5 수치 feature와 실험 visibility, group index와 usable delay

EIT/CPT는 finite half-height crossing만 있으면 linewidth를 hero로 올리고
(`gabes/schemes/absorption.py:771-795`), detector noise와 modeled contrast,
samples/FWHM은 보지 않는다. AT는 두 local maximum과 1% central contrast를
요구해 더 보수적이다 (`gabes/schemes/absorption.py:723-769`).

interpolated width, resolution/edge status, `T_res`, local group-delay estimate
`(n_g-1)L/c`와 visibility-not-assessed 문구는 기존 배열의 `O(N)` 후처리다.
추가 OBE solve가 없으며, 현재 `n_g`를 finite-pulse 예측으로 과신하는 위험을
줄인다.

## 6. 기존 물리 개선안의 비용과 우선순위

| 개선 | 물리 보존 여부 | 계산비용 | 권고 |
|---|---|---:|---|
| scan/velocity operator 분리 | 잘못 섞인 항만 분리, 모델 차원 동일 | 사실상 0 | 즉시 P1 |
| center-clustered 601 + 보간/status | 동일 OBE를 적절한 grid에 평가 | 이번에는 `-1.2%`; 후처리 `0.226 ms` | broad-wing 회귀 후 P1 |
| pure-isotope collinear `Delta k` | analytic wavevector 보정 | scalar coefficient, 추가 solve 없음 | P2 |
| visibility/group-delay trust status | 계산 배열 해석만 개선 | `O(N)`, 추가 solve 없음 | P1과 함께 |
| sourced pressure shift·elastic width | line center·homogeneous width 개선 | scalar/table, solve 차원 동일 | provenance 확보 후 P1 |
| full VCC/diffusion | velocity redistribution·Dicke/Ramsey 물리 추가 | velocity class 결합, separability 상실 | 검증 data가 있는 opt-in 연구 solver |
| hyperfine/Zeeman Lambda | line assignment·편광·optical pumping 추가 | 상태 수와 Liouvillian 급증 | 특정 transition의 opt-in reference |

저차 pressure shift와 elastic width는 runtime 관점에서는 싸다. 그러나 현재 Ne
상수는 fixed broadening-only이며 계수 table과 pressure shift가 TODO다
(`gabes/constants.py:79-89`). checklist가 effort를 `large`로 둔 이유는 연산량이
아니라 gas/species/line별 출처·부호·단위·불확도·유효범위를 정하는 연구 작업
때문이다 (`docs/checklist.json:209-243`). 검증 데이터 없는 phenomenological
Dicke knob는 빠르더라도 reference 신뢰도를 낮추므로 추가하지 않는 편이 맞다.

full VCC는 기본 EIT의 711 velocity class가 독립이라는 구조를 깨뜨린다. 실제로
VCC·diffusion·Dicke·Ramsey 효과가 별도 운동 물리라는 근거도 있다
([Firstenberg et al., arXiv:0801.3660](https://arxiv.org/abs/0801.3660)). 따라서
one-scheme/one-kernel/one-dataset opt-in solver가 타당하다. hyperfine/Zeeman
mode도 compact 3-level 기본값을 바꾸지 말고 held-out spectrum이 있는 별도
reference mode로 두어야 한다.

## 7. 동작을 바꾸지 않는 순수 코드 최적화

### 7.1 `line_strength`를 navigate-only로 이동 — 가장 큰 확실한 이득

`ParamSpec.recompute` 기본은 `True`다 (`gabes/schemes/base.py:44-56`). Lambda의
`line_strength`는 이를 그대로 상속한다
(`gabes/schemes/absorption.py:557-558`). 그러나 이 값은 Hamiltonian이나
`chi_bar`에 들어가지 않고 physical susceptibility를 만들 때만 곱해진다
(`gabes/schemes/absorption.py:664-682`; `gabes/observables.py:520-537`).

기본 EIT에서 `line_strength=1 -> 1.5`를 재검증한 결과:

- `chi_bar`는 bitwise identical
- 전체 재 solve median `33.293 ms`
- 같은 raw를 새 strength로 headless remap하는 median `0.1126 ms`
- 최종 metric은 동일, 불필요한 경로 차이는 약 `296x`

`line_strength`에 `recompute=False`를 주고 `observables()`가 live
`params["line_strength"]`를 읽게 하면 물리와 사용자-visible 결과를 보존하면서
재계산을 없앨 수 있다. cache equivalence, `chi_bar` 불변, `alpha` 선형 scaling을
회귀로 고정해야 한다.

### 7.2 작은 immutable template cache — 낮은 우선순위

현재 `_medium_from_params()` median은 `16.7 us`, `atoms.lambda3()`는 `184.5 us`다.
scalar key의 immutable medium tuple과 dissipator/template를 캐시하면 긴 parameter
sweep에서 누적 비용을 줄일 수 있다. 그러나 합쳐도 기본 EIT solve의 약 `0.60%`라
우선순위는 낮고, mutable ndarray를 공유하지 않도록 해야 한다.

### 7.3 affine scan-constant hoist는 이번에는 권고하지 않음

과거 보고서는 velocity loop 밖으로 `base + s*A_coef`를 옮기는 후보를 약
`1.13x`로 기록했다. 현재 601x711 기본 EIT에서 같은 산술 순서를 보존한 독립
prototype을 다시 재자 결과:

- current kernel `29.299 ms`
- hoisted prototype `29.732 ms`
- output bitwise identical, 그러나 `0.985x`로 약 `1.5%` 느림

LU solve가 지배하고 추가 temporary matrix가 이득을 상쇄한 것으로 보인다.
현재 환경의 증거로는 이 최적화를 구현하지 않는 편이 맞다. 프로파일과 여러 CPU
환경에서 일관된 이득이 확인되기 전에는 checklist에 올릴 가치가 낮다.

## 8. 문서, 테스트, 예제의 reference 품질

### 문서

- README의 Lambda 한 줄은 기능을 정확히 요약하지만 scalar/weak-probe/
  natural-Rb effective-medium 한계를 말하지 않는다 (`README.md:13`).
- README roadmap은 slow-light/group-index readout을 미래 항목으로 쓰지만
  local group index는 이미 출력한다 (`README.md:18-20`;
  `gabes/schemes/absorption.py:771-787`). pulse propagation이 미래라는 식으로
  문구를 좁혀야 한다.
- module docstring은 knob가 `Gamma` 단위라고 하지만 UI는 이미 MHz/kHz 기반
  `physical-units-v3`다 (`gabes/schemes/absorption.py:16-18, 462-465`).
- 사용자 가이드는 full polarization/Zeeman model이 아님을 정직하게 밝힌다
  (`docs/Userguide/GABES_User_Guide_v2.html:828-830`). 반면 EIT 그림의 분산을
  조건 없이 “느린 빛의 서명”이라고 부르며 transmission loss와 finite-pulse
  한계를 같이 말하지 않는다 (`docs/Userguide/GABES_User_Guide_v2.html:586-602`).
- `buffer_ground_relax_khz`의 help는 buffer collision처럼 보이지만 구현은 단일
  Raman-coherence dephasing이다 (`gabes/schemes/absorption.py:549-553`;
  `gabes/atoms.py:168-190`). 명칭과 help를 더 좁혀야 한다.

### 테스트

현재 테스트는 다음에 강하다.

- AT split `approximately Omega_c`, power/diameter Rabi scaling
  (`tests/test_absorption.py:110-129`)
- warm angle에서 EIT broadening, cold transparency, sub-natural CPT
  (`tests/test_absorption.py:132-175`)
- unresolved hero/status contract (`tests/test_absorption.py:195-216`)
- Lambda/Rydberg Numba kernel과 NumPy reference parity
  (`tests/test_kernels.py:102-136`)
- headless와 normal render contract (`tests/test_headless_observables.py:38-89`;
  `tests/test_schemes_render.py:97-107, 140-177`)

그러나 다음 회귀는 아직 없다.

- Doppler-off angle invariance
- default FWHM/group-index high-resolution convergence와 interpolated edge
- collinear hyperfine `Delta k`
- natural-Rb scalar status
- visibility와 group delay/loss의 동시 표시
- `line_strength` navigate-only cache equivalence
- Lambda example artifact의 parameter/code/grid provenance

### 예제

실행 가능한 Lambda 전용 script/config는 없다. 가장 가까운 코드 예제는
`tests/test_absorption.py:110-216`의 direct API 호출이다. 가이드의 `eit.png`와
`at.png`는 현상을 잘 보여 주지만 생성 script와 manifest가 없다. EIT PNG의
`Omega_c=6.00 MHz`, `gamma_gg=10.0 kHz`는 현재 EIT 기본값 약 `17.24 MHz`,
`57.46 kHz`와 다르고, CPT 그림은 없다.

headless power/diameter, angle, dephasing, cell-length sweep을 CSV/Markdown으로
저장하고 parameter manifest, commit, grid/status, model-scope badge를 포함하는
Lambda 예제를 추가하면 solver를 바꾸지 않고도 실험 reference 재현성을 크게
높일 수 있다.

## 9. 검증

현재 working tree에서 다음을 실행했다.

```powershell
python -m pytest tests/test_absorption.py tests/test_kernels.py `
  tests/test_headless_observables.py tests/test_schemes_render.py `
  tests/test_sas.py tests/test_docs_consistency.py -q
# 71 passed in 26.56s

python -m pytest -q
# 476 passed in 245.27s
```

추가로 기본 EIT/AT/CPT, uniform 601/9,601점과 center-clustered 601점 수렴,
Doppler-off/on 0/10 mrad, 공용 resolution helper, `line_strength` remap,
medium/factory microbenchmark와 affine-hoist prototype을 직접 계산했다.

## 10. 최종 판단과 권고 순서

Scheme 2는 정상상태 3준위 Lindblad OBE, Maxwell velocity average, lab-facing
Rabi scaling과 Beer-Lambert propagation을 결합한다. 따라서 **EIT/AT/CPT의 기본
scale과 노브 민감도를 빠르게 보는 반정량적 실험 reference**로는 유용하다.
특히 AT splitting과 coupling power/diameter sweep은 실험 전 sanity check에
적합하다.

그러나 현재 기본 EIT의 hero linewidth와 group index는 grid 미수렴이고,
beam-angle residual Doppler는 cold-limit를 위반한다. 특정 hyperfine/Zeeman,
clock-CPT, buffer-cell, finite-pulse slow-light의 절대 reference로 사용하면 안 된다.

권고 순서는 다음과 같다.

1. `lambda-spectrum-validity` P1을 완료한다: scan/velocity operator 분리,
   center-clustered grid, interpolated width, samples/edge/visibility status,
   transmission과 group-delay/`n_g` 동시 표시.
2. `line_strength`를 navigate-only로 옮겨 약 `296x`의 불필요한 재계산을 없앤다.
3. pure-isotope collinear `Delta k`, natural-Rb scalar 경고와 Raman-dephasing
   명칭을 추가한다.
4. 출처와 uncertainty가 있는 pressure shift/elastic width만 저차 모델에 넣는다.
5. full VCC와 hyperfine/Zeeman Lambda는 검증 spectrum과 runtime budget이 정해진
   뒤 opt-in research/reference mode로 시작한다.

이 순서는 compact 3-level 모델의 유용한 실제 물리를 유지하면서, 거의 같은
interactive 비용으로 실험가가 정밀해 보이는 숫자를 과신할 위험을 가장 크게
줄인다.
