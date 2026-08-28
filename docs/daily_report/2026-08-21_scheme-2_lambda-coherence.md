# 2026-08-21 Scheme 2 물리 검토 — Lambda coherence

## 1. 오늘의 선택과 현재 다섯 scheme

서울 현지 날짜의 day-of-month는 `21`이므로

```text
n = (day mod 5) + 1 = (21 mod 5) + 1 = 2
```

이다. 드롭다운 순서는 `gabes/schemes/__init__.py:12-24`의 `_SCHEMES`가
정하며, 현재 정의는 다음과 같다.

| 순번 | 등록 인스턴스 | UI scheme | 핵심 범위 |
|---:|---|---|---|
| 1 | `SASScheme()` | OD / SAS | Doppler 흡수, 포화흡수, hyperfine pumping |
| 2 | `LambdaScheme()` | Λ coherence (EIT / AT / CPT) | 축약 3준위 Λ 결맞음 |
| 3 | `RydbergEITScheme()` | Rydberg-EIT electrometry | cascade EIT, microwave AT, electrometry |
| 4 | `MagnetoScheme()` | Hanle / EIA / NMOR | Zeeman OBE의 투과·회전 |
| 5 | `FWMScheme()` | FWM | seeded mean-field FWM, SFWM biphoton |

`eit`, `at`, `cpt`는 직접 호출과 테스트를 위한 alias일 뿐 드롭다운 항목을
늘리지 않는다 (`gabes/schemes/__init__.py:27-39`). README의 표도 같은 다섯
엔진을 설명한다 (`README.md:8-16`). 따라서 오늘 대상은 두 번째
`LambdaScheme()`이다.

검토 시점은 branch `main`, HEAD
`c964a724b1ad370e53af1a89c78f8df04fd37983`이다. 작업 트리에는 기존 사용자
변경이 많이 있으므로 모두 보존했다. `gabes/schemes/absorption.py`,
`gabes/beam.py`, `gabes/kernels.py`, `tests/test_absorption.py`에는 로컬 diff가
없다. `gabes/atoms.py`의 로컬 변경은 일반 `collapse_ops` 지원 추가이고
`lambda3()`의 대칭 붕괴·Raman dephasing 의미는 그대로다. 최근 공용
`gabes/lineshape.py`에는 SAS용 해상도 진단이 추가됐지만 Lambda 출력은 아직
이를 사용하지 않는다.

## 2. 먼저 검색한 기존 제안, TODO, issue note

다음 범위를 먼저 확인했다.

- Scheme 2 선행 보고서 7개:
  `docs/daily_report/2026-06-26_scheme-2_lambda-coherence.md`,
  `2026-07-01`, `2026-07-06`, `2026-07-16`, `2026-07-31`, `2026-08-01`,
  `2026-08-06`
- 통합 계획과 TODO: `docs/checklist.json`, `gabes/constants.py:79-89`
- 코드: `gabes/schemes/absorption.py`, `gabes/atoms.py`, `gabes/beam.py`,
  `gabes/kernels.py`, `gabes/observables.py`, `gabes/lineshape.py`
- 테스트: `tests/test_absorption.py`, `tests/test_kernels.py`,
  `tests/test_headless_observables.py`, `tests/test_schemes_render.py`,
  공용 line-shape 테스트가 있는 `tests/test_sas.py`
- 문서·예제: `README.md`, `CLAUDE.md`,
  `docs/Userguide/GABES_User_Guide_v2.html`,
  `docs/Userguide/userguide_assets/eit.png`, `at.png`, `analysis/`

별도의 로컬 issue/proposal/roadmap 파일과 Lambda 전용 실행 예제는 발견하지
못했다. 그러나 **기존 개선안은 명확히 존재한다**. 현재 통합 checklist의 관련
항목은 다음과 같다.

| 기존 항목 | 상태 | 현재 판단 |
|---|---|---|
| `lambda-spectrum-validity` | P1, `ready` | residual Doppler 오류, grid 미수렴, visibility/group-delay 표시를 함께 고치는 핵심 항목 (`docs/checklist.json:439-461`) |
| `collisional-coefficient-provenance-and-pressure-shift` | P1, `ready` | gas/species/line별 출처·단위·부호·불확도를 먼저 정한 뒤 scalar shift/width를 적용 (`docs/checklist.json:209-243`) |
| `full-velocity-changing-collision-kernel` | P2, `parked` | 대상 scheme·충돌 kernel·검증 spectrum·runtime budget이 정해질 때만 opt-in 연구 solver로 시작 (`docs/checklist.json:563-582`) |
| `lambda-hyperfine-resolved-manifold` | P2, `parked` | 특정 isotope/transition/polarization과 holdout spectrum을 고른 뒤 opt-in reference mode로 구현 (`docs/checklist.json:585-603`) |

이미 반영된 과거 제안도 있다. coupling power와 1/e² beam diameter가 보정된
Rabi anchor를 `sqrt(P)/d`로 바꾸며 (`gabes/schemes/absorption.py:529-548,
583-591`; `gabes/beam.py:9-20`), headless observables는 figure 없이 metric/table만
계산한다 (`gabes/schemes/absorption.py:462-466`;
`tests/test_headless_observables.py:38-67`). 두 항목 모두 물리 차원이나 solve 수를
늘리지 않는다.

## 3. 현재 구현이 실제로 푸는 물리

### 3.1 3준위 정상상태 Lindblad OBE

`atoms.lambda3()`는 `g1`, `g2`, `e`의 세 상태를 만들고 들뜬 상태가 두
바닥상태로 각각 `Γ/2`씩 붕괴하게 한다. `gamma_gg`는 population reload가 아니라
`rho_g1g2`와 그 켤레를 직접 감쇠시키는 Raman-coherence `T2` rate다
(`gabes/atoms.py:168-190`).

Hamiltonian은 고정 약한 probe `Omega_p = 10^-3 Γ`, coupling `Omega_c`, coupling
detuning, two-photon scan을 직접 포함한다
(`gabes/schemes/absorption.py:29-31, 633-656`). 그러므로 다음은 단순 선모양
합성이 아니라 정상상태 density-matrix 해에서 나온다.

- two-photon destructive interference에 의한 EIT transparency와 정상분산
- 강한 coupling이 만드는 Autler–Townes dressed doublet
- Raman coherence의 좁은 dark resonance
- coupling strength와 ground-coherence dephasing 변화의 올바른 정성적 경향

다만 `CPT`는 독립적인 bichromatic clock-CPT 엔진이 아니다. 같은 weak-probe
Lambda 모델의 좁은 scan preset이다 (`absorption.py:468-500, 624-631`). 균형 잡힌
두 광장, 실제 optical pumping, light shift, magnetic sensitivity를 포함하는
clock reference로 읽으면 안 된다.

### 3.2 scalar D-line medium과 lab-facing knob

Rb/Cs와 D1/D2 선택은 온도별 원자수밀도, 자연선폭, 질량, 파장, reduced dipole을
바꾼다 (`gabes/schemes/absorption.py:419-443, 521-526`). coupling은 절대 CG,
편광, 공간 overlap을 처음부터 계산하지 않고 실험적으로 보정한 Rabi anchor를
power와 diameter에 따라 옮긴다. 한 operating point의 측정 `Omega_c`를 power/waist
sweep으로 확장하는 용도에는 정직하고 유용하다.

그러나 이 매질에는 구체적인 `Fg -> Fe`, CG coefficient, Zeeman population,
hyperfine optical pumping이 없다. 특히 `Rb (natural)`은 85Rb와 87Rb의 서로 다른
광학·Raman susceptibility를 합산하지 않는다. 가장 abundance가 큰 isotope의
line constant를 쓰고, abundance-weighted density와 mass를 한 scalar Lambda에 넣는다
(`absorption.py:430-442`; `gabes/species.py:127-180`). 현재 50 °C D1에서 natural
Rb와 pure 85Rb는 `N`과 optical centroid가 같고 mass만 약 `0.655%` 다르며,
`chi_bar` 최대 상대 차이는 약 `0.323%`였다. 이는 natural-Rb의 두 isotope가 한
dark state에 결맞게 참여한다는 뜻이 아니라 **유효 scalar medium**이라는 뜻이다.

### 3.3 Doppler 평균, Beer–Lambert, group index

Doppler-on Lambda는 4-sigma Maxwell grid를 `dv=2 m/s`로 만들고 각
scan×velocity 정상상태를 푼 뒤 probe coherence만 속도 가중 수축한다
(`gabes/schemes/absorption.py:118-140`; `gabes/kernels.py:413-495`). 기본 50 °C
EIT에서는 711 velocity class와 601 scan point를 사용한다. beam angle은
`|Delta k|/k`를 ground level의 velocity coefficient에 추가한다
(`gabes/beam.py:42-55`; `gabes/atoms.py:176-189`).

이후 `chi_bar`는 density와 dipole을 포함한 `chi_phys`,
`alpha = k Im(chi)`, `T = exp(-alpha L)`로 바뀐다
(`gabes/observables.py:524-542`). group index는 dilute-medium 근사에서
`Re(chi)`의 국소 수치 미분이다 (`gabes/observables.py:550-558`). 이는 실제
정상상태 선형응답 물리지만 finite-pulse propagation, distortion, detector SNR을
계산한 slow-light experiment model은 아니다.

외부 물리 기준과도 대조했다.

- Li와 Xiao의 room-temperature Rb D1 실험은 weak probe/strong pump의
  Doppler-free Lambda EIT와 Doppler를 포함한 단순 이론의 정성적 유효성을 직접
  보인다: <https://doi.org/10.1103/PhysRevA.51.R2703>
- Ghosh 등의 hot-vapor 모델은 VCC, transit/influx, Raman coherence decay가 서로
  다른 역할을 하므로 하나의 `gamma_gg`로 치환할 수 없음을 보여 준다:
  <https://arxiv.org/abs/0901.3790>
- EIT review는 continuous-wave susceptibility와 pulse propagation을 구분한다.
  따라서 현재 local `n_g`를 pulse delay·fidelity 예측으로 승격하면 안 된다:
  <https://doi.org/10.1103/RevModPhys.77.633>

## 4. 실험 원자물리 reference로서의 적합성

현재 Scheme 2는 **실제 물리를 구현한 반정량적 실험 planning/reference**다.
교과서 곡선을 그리는 장난감보다 훨씬 낫지만, 특정 alkali Lambda 실험의 절대
spectroscopy 또는 clock reference는 아니다.

| 연구 질문 | 적합성 | 이유 |
|---|---|---|
| coupling power/diameter 변화로 EIT에서 AT로 넘어가는 scale | 좋음 | 정상상태 OBE와 보정된 `sqrt(P)/d` Rabi를 사용 |
| AT split 대 측정 `Omega_c` sanity check | 좋음 | 두 local maximum과 `split ≈ Omega_c`가 직접 테스트됨 |
| Raman `T2`가 EIT/CPT를 흐리는 방향 | 유용 | 실제 coherence dephasing이지만 population exchange·transit·diffusion과 분리되지 않음 |
| warm-vapor beam-alignment sensitivity | 수정 후 유용 | residual Doppler 방향은 맞지만 현재 cold-limit 오류가 있음 |
| 특정 isotope/hyperfine/polarization의 절대 contrast·linewidth | 낮음 | scalar 3-level, 대칭 branching, Zeeman/optical pumping 부재 |
| clock-CPT light shift·자기장 민감도 | 낮음 | weak-probe preset이며 balanced bichromatic pumping이 아님 |
| slow-light pulse delay·왜곡·효율 | 낮음 | local group index만 있고 pulse propagation과 detection이 없음 |

실험 전에는 `Omega_c`, ground-coherence dephasing, alignment, temperature,
cell length의 경향과 대략적인 scale을 찾는 데 쓸 수 있다. 논문 피팅이나 장비
specification에는 실제 transition assignment, beam profile, pressure/collision
budget, detector threshold와 held-out spectrum이 추가로 필요하다.

## 5. 2026-08-21 직접 재계산과 핵심 문제

### 5.1 기본 출력

Numba warm-up 뒤 각 기본값을 5회 계산한 median이다. 절대 시간은 이 환경에
한정된다.

| regime | compute median | 현재 metric |
|---|---:|---|
| EIT | `34.32 ms` | `T_res=0.013558`, 표시 FWHM `0.45968 MHz`, `n_g=100,305` |
| AT | `3.40 ms` | split `46.0 MHz`, expected `46.0 MHz`, `T_center=0.98909` |
| CPT | `3.27 ms` | `T_res=0.70633`, 표시 FWHM `923.32 kHz`, `n_g=95,255` |

AT는 강한 coupling scale을 잘 읽는다. 그러나 EIT와 CPT의 표시 linewidth는 아래
grid 문제 때문에 신뢰 가능한 숫자가 아니다.

### 5.2 P1 — 기본 EIT FWHM과 group index 미수렴

`LambdaScheme._scan()`은 모든 regime에 601점 uniform grid를 사용한다
(`gabes/schemes/absorption.py:609-631`). `window_fwhm()`은 half-height edge를
보간하지 않고 바깥 sample의 차를 반환한다 (`gabes/lineshape.py:53-71`). 기본
EIT 표시 FWHM은 grid의 정확히 두 칸이다.

| scan | 중앙 step | sample FWHM | edge-interpolated FWHM | center `n_g` | warm solve |
|---|---:|---:|---:|---:|---:|
| uniform 601 | `0.22984 MHz` | `0.45968 MHz` | `0.22990 MHz` | `100,305` | `33.8 ms` |
| `sinh(5u)` center-clustered 601 | `0.01549 MHz` | `0.18616 MHz` | `0.15914 MHz` | `112,708` | `34.7 ms` |
| uniform 9,601 | `0.014365 MHz` | `0.17238 MHz` | `0.15910 MHz` | `112,693` | `486 ms` |

현재 hero FWHM은 수렴 보간값보다 약 `2.89x` 크고, `n_g`는 약 `11%` 작다.
같은 601점을 중앙에 집중하면 FWHM은 9,601점 결과와 `0.03%`, group index는
`0.02%` 안에서 일치했다. 15회 median의 clustered/uniform runtime 비는
`1.026`으로, solve 수를 늘리지 않고 checklist의 20% budget을 충분히 만족한다.

공용 `subdoppler_feature()`는 이미 interpolated edges, samples per FWHM,
scan-edge clearance와 resolution status를 구현한다
(`gabes/lineshape.py:20-50, 141-153`). 현 Lambda 기본 배열에 적용하면 EIT는
약 `1.00`, CPT는 약 `3.19` sample/FWHM으로 둘 다 `resolution-limited`가 된다.
601점에서 이 진단의 median은 `0.214 ms`였다. 다만 이 helper의 running-median
background 정의는 Lambda의 전체-window floor와 완전히 같지 않으므로, 코드를
그대로 호출하기보다 **edge interpolation과 status 구조를 공유**하는 편이 안전하다.

### 5.3 P1 — residual Doppler가 cold spectrum도 바꿈

현재 affine coefficient는

```text
A_coef = dL/ds - S_v
B_coef = S_v
L(s, kv) = base + s A_coef + kv B_coef
```

로 조립된다 (`gabes/schemes/absorption.py:45-57`;
`gabes/kernels.py:454-482`). angle-dependent ground coefficient가 scan 항에도
들어가므로, velocity가 0인 Doppler-off 계산에서도 angle이 spectrum을 바꾼다.
물리적인 residual Doppler `Delta k · v`는 `v=0`에서 사라져야 한다.

현재 cold EIT, angle 0 대 10 mrad 재현값:

- `max |Delta chi_bar| / max |chi_bar| = 2.94204%`
- 3 mm cell에서 `max |Delta T| = 7.88490e-3`
- `|Delta k|/k`: `0`에서 `0.00999996`

scan optical-shift operator와 velocity-only residual operator를 분리하면 matrix
차원, scan 수, velocity class 수가 그대로이므로 runtime 증가 없이 정확성을
회복할 수 있다. `doppler="off"`의 0/10 mrad 전체 complex-`chi` invariance와
Numba/NumPy parity를 회귀로 고정해야 한다.

### 5.4 P2 — collinear optical frequency difference가 빠짐

`_two_photon_k_ratio()`는 probe와 coupling에 같은 `medium["k_vec"]`를 넘겨
angle 0에서 residual을 정확히 0으로 만든다
(`gabes/schemes/absorption.py:593-597`; `gabes/beam.py:36-55`). 실제 hyperfine
Lambda의 두 광장은 ground splitting만큼 주파수가 다르므로 co-propagating이어도
작은 종방향 `Delta k`가 남는다. 2026-08-06 보고서의 current species constants
기반 25 °C Gaussian FWHM은 85Rb `4.07 kHz`, 87Rb `9.07 kHz`, 133Cs
`9.86 kHz`였다 (`docs/daily_report/2026-08-06_scheme-2_lambda-coherence.md:210-232`).

기본 EIT에는 작은 보정이지만 협폭 clock-CPT에서는 무시하기 어렵다. pure-isotope
mode에서 두 실제 wavevector를 만들면 기존 velocity coefficient의 scalar만
바뀌므로 추가 solve는 없다. 앞의 cold-limit operator 분리를 먼저 해야 한다.

### 5.5 P2 — 수치 feature와 실험 visibility, group index와 usable delay

현재 EIT/CPT는 finite half-height crossing만 있으면 linewidth를 hero로 올린다
(`absorption.py:771-795`). modeled contrast, samples/FWHM, detector noise floor는
보지 않는다. 반대로 AT는 두 local maximum과 1% central contrast를 요구해 더
보수적이다 (`absorption.py:723-769`).

다음은 모두 기존 배열의 `O(N)` 후처리이므로 추가 OBE solve가 없다.

1. interpolated width, samples/FWHM, edge clearance를 status에 포함한다.
2. noise 입력이 없으면 `numerically resolved; experimental visibility not assessed`
   라고 표시한다.
3. `n_g`와 함께 `T_res`, local group-delay estimate `(n_g-1)L/c`, 제한적
   delay-bandwidth proxy를 보이고 finite-pulse prediction이 아님을 명시한다.

## 6. 기존 물리 개선안의 비용과 우선순위

| 개선 | 물리 보존 | 계산비용 | 권고 |
|---|---|---:|---|
| scan/velocity operator 분리 | 오류만 제거, 모델 차원 동일 | 사실상 0 | 즉시 P1 |
| center-clustered 601 + edge interpolation/status | 동일 OBE를 더 적절한 grid에서 평가 | 이번 median `+2.6%`; 후처리 `~0.2 ms` | broad wing 회귀 후 즉시 P1 |
| collinear hyperfine `Delta k` | analytic wavevector correction | scalar coefficient, 추가 solve 없음 | pure-isotope부터 P2 |
| contrast/visibility/group-delay status | 계산된 배열의 해석만 개선 | `O(N)`, 추가 solve 없음 | P1과 함께 적용 |
| Raman-dephasing 명칭과 natural-Rb scalar 경고 | 오해만 줄임 | runtime 0 | 즉시 문서/UI 수정 |
| sourced pressure shift·elastic width | line centre·homogeneous width 개선 | scalar/table, solve 차원 동일 | 계수 provenance 확보 후 P1 |
| full VCC/diffusion | velocity redistribution과 Dicke/Ramsey physics | velocity class 결합, 현재 fast separability 상실 | 검증 data가 있는 opt-in 연구 solver |
| hyperfine/Zeeman Lambda | line assignment·편광·optical pumping | 더 큰 state/Liouvillian, preset 재정의 | 특정 transition의 opt-in reference mode |

`gabes/constants.py:79-89`의 Ne 값은 fixed broadening-only scalar이며 TODO도
gas/species/line coefficient table과 pressure shift를 요구한다. pressure shift와
elastic width는 runtime 측면에서는 저렴하지만, checklist가 effort `large`로 둔
이유는 **계수의 출처·단위·부호·불확도와 지원 범위**를 정하는 연구 작업 때문이다
(`docs/checklist.json:209-243`). 검증 dataset 없는 phenomenological Dicke knob는
빠르더라도 실험 reference 신뢰도를 낮추므로 넣지 않는 현재 정책이 타당하다.

full VCC는 기본 50 °C EIT의 711 velocity class가 서로 독립이라는 현재 구조를
깨뜨린다. 따라서 단순 scalar overhead가 아니며, one-scheme/one-kernel/one-dataset
reference solver로 제한해야 한다. hyperfine/Zeeman mode도 compact 3-level의
설명 가능성과 속도를 기본값에서 희생하지 말고 opt-in으로 분리하는 것이 맞다.

## 7. 동작을 바꾸지 않는 순수 코드 최적화

### 7.1 `line_strength`를 navigate-only로 이동 — 가장 큰 즉시 효과

`ParamSpec.recompute`의 기본은 `True`이고
(`gabes/schemes/base.py:44-56`), Lambda `line_strength`는 이를 그대로 상속한다
(`gabes/schemes/absorption.py:557-558`). 하지만 `line_strength`는 Hamiltonian이나
`chi_bar`에 들어가지 않고 solve 결과의 physical susceptibility를 곱할 때만 쓰인다
(`absorption.py:664-668, 677-684`; `gabes/observables.py:524-537`).

오늘 기본 EIT에서:

- line-strength-only 재 solve median: `33.70 ms`
- 같은 raw의 headless remap median: `0.109 ms`
- 불필요한 경로 차이: 약 `308x`
- `chi_bar`는 bitwise 동일했고 최종 metric도 동일했다.

`line_strength`에 `recompute=False`를 지정하고 `observables()`가
`raw["ls"]` 대신 live `params["line_strength"]`를 읽게 하면 사용자-visible
동작과 물리를 보존하면서 이 재계산을 없앨 수 있다. cache path, 직접 호출,
여러 line strength에서 `chi_bar` 불변·`alpha` 선형 scaling·metric equivalence를
테스트해야 한다.

### 7.2 affine kernel의 scan-constant matrix hoist

`_affine_scan_chi_real()`은 각 scan point의 모든 velocity에서 9x9 matrix 전체를

```text
base + s*A_coef + kv*B_coef
```

로 다시 조립한다 (`gabes/kernels.py:426-449`). `base + s*A_coef`와 고정 trace
row를 scan당 한 번 만들고 velocity loop에는 copy와 `kv*B_coef`만 남길 수 있다.
2026-08-01/08-06 진단은 bit-identical `chi`와 약 `1.13x` 개선을 기록했고, 해당
kernel은 이후 바뀌지 않았다
(`docs/daily_report/2026-08-06_scheme-2_lambda-coherence.md:292-299`).
Numba/NumPy parity와 Doppler on/off 회귀를 유지하면 동작 보존형 최적화다.

### 7.3 낮은 우선순위의 object/template cache

현재 median은 `_medium_from_params()` 약 `23.2 us`, `atoms.lambda3()` 약
`189 us`로 합쳐도 기본 EIT solve `34 ms`의 1% 미만이다. immutable template 또는
작은 scalar-key cache는 대규모 batch에는 도움이 되지만, mutable array 공유와
cache-key 복잡성을 감수할 만큼 현재의 우선 병목은 아니다. 먼저
`line_strength` cache boundary와 kernel hoist를 처리하는 편이 낫다.

이미 구현된 `headless_observables()`는 Matplotlib을 건너뛰므로 자동 sweep과
보고서 경로에서 계속 사용해야 한다 (`tests/test_headless_observables.py:46-67`).

## 8. 문서, 테스트, 예제의 reference 품질

### 문서

- README의 Lambda 한 줄은 기능은 정확히 요약하지만 scalar/weak-probe/
  natural-Rb effective-medium 한계를 말하지 않는다 (`README.md:13`).
- README roadmap은 group-index readout을 아직 미래 항목으로 둔다
  (`README.md:18-20`). 실제 출력은 이미 `absorption.py:771-787`에 있다.
- module docstring은 knob가 `Gamma` 단위라고 적지만 실제 UI는 MHz/kHz의
  `physical-units-v3`이다 (`absorption.py:16-18, 462-465`).
- 사용자 가이드는 Lambda가 full polarization/Zeeman model이 아님을 정직하게
  밝힌다 (`docs/Userguide/GABES_User_Guide_v2.html:813-815`). 반면 EIT 그림의
  가파른 분산을 조건 없이 “느린 빛의 서명”이라고 부르며 loss와 finite-pulse
  한계를 함께 표시하지 않는다 (`...:586-602`).
- `buffer_ground_relax_khz`의 UI help는 buffer collision처럼 보이지만 실제
  구현은 단일 Raman-coherence dephasing이다
  (`absorption.py:549-553`; `atoms.py:168-190`). 명칭을
  `ground-coherence dephasing / 2pi`로 바꾸는 것이 정확하다.

### 테스트

현재 테스트는 다음 내부 불변량과 UI contract에 강하다.

- AT split `≈ Omega_c`, power/diameter Rabi scaling
  (`tests/test_absorption.py:110-129`)
- warm angle에서 EIT broadening, cold EIT transparency, sub-natural CPT
  (`tests/test_absorption.py:132-175`)
- unresolved hero status (`tests/test_absorption.py:195-216`)
- Lambda/Rydberg Numba affine kernel과 NumPy reference의 `1e-11` parity
  (`tests/test_kernels.py:102-136`)
- headless와 figure-render contract (`tests/test_headless_observables.py:38-89`;
  `tests/test_schemes_render.py`)

그러나 다음은 아직 테스트하지 않는다.

- Doppler-off angle invariance
- default FWHM/group-index grid convergence와 interpolated edge
- collinear hyperfine `Delta k`
- natural-Rb scalar semantics/status
- experimental visibility와 group delay/loss 동시 표시
- `line_strength` navigate-only cache equivalence
- Lambda example artifact의 parameter/code/grid provenance

### 예제

실행 가능한 Lambda 전용 script/config는 없다. 가장 가까운 코드 예제는
`tests/test_absorption.py:110-216`의 direct API 호출이다. 사용자 가이드의
`eit.png`와 `at.png`는 현상을 시각적으로 잘 보여 주지만 생성 script와 manifest가
없다. 특히 EIT PNG 제목의 `Omega_c=6.00 MHz`, `gamma_gg=10.0 kHz`는 현재 EIT
기본값 약 `17.24 MHz`, `57.46 kHz`와 다르고, CPT 그림은 없다.

낮은 비용으로 headless power/diameter, angle, dephasing, cell-length sweep을
CSV/Markdown으로 저장하는 Lambda 예제를 추가하면 실험 reference 재현성이 크게
좋아진다. output에는 parameter manifest, commit, grid/status, 모델-scope badge를
포함해야 한다. 이는 solver 물리나 runtime을 바꾸지 않는 문서성 개선이다.

## 9. 검증

현재 working tree에서 다음을 실행했다.

```powershell
python -m pytest tests/test_absorption.py tests/test_kernels.py `
  tests/test_headless_observables.py tests/test_schemes_render.py `
  tests/test_sas.py -q
# 70 passed in 19.76s

python -m pytest -q
# 439 passed in 82.42s
```

추가로 기본 EIT/AT/CPT, uniform 601/9,601점과 center-clustered 601점 수렴,
Doppler-off 0/10 mrad 비교, current `subdoppler_feature` 상태/비용,
natural-Rb 대 pure-isotope medium, line-strength remap equivalence와 microbenchmark를
실행했다.

## 10. 최종 판단과 권고 순서

Scheme 2는 실제 3준위 Lindblad OBE, Maxwell velocity average, lab-facing Rabi
scaling, Beer–Lambert propagation을 결합한다. 따라서 **EIT/AT/CPT의 기본 scale과
노브 민감도를 빠르게 보는 반정량적 실험 reference**로는 유용하다. 특히 AT
splitting과 coupling-power/diameter sweep은 실험 전 sanity check에 적합하다.

가장 먼저 할 일은 더 큰 atomic manifold가 아니다.

1. `lambda-spectrum-validity` P1을 완료한다: scan/velocity operator 분리,
   center-clustered grid, interpolated width, samples/edge/visibility status,
   group index와 transmission/group-delay의 동시 표시.
2. `line_strength`를 navigate-only로 옮겨 불필요한 solve를 제거한다.
3. pure-isotope의 collinear hyperfine `Delta k`, natural-Rb scalar 경고,
   Raman-dephasing 명칭을 추가한다.
4. 출처와 uncertainty가 있는 pressure shift/elastic width만 저차 모델에 넣는다.
5. full VCC와 hyperfine/Zeeman Lambda는 검증 spectrum과 runtime budget이 정해진
   뒤 opt-in research/reference mode로 시작한다.

이 순서는 현재 compact model의 실제 물리를 보존하면서도, 거의 같은 interactive
비용으로 실험가가 숫자를 과신할 위험을 가장 크게 줄인다.
