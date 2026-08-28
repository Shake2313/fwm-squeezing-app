# 2026-08-08 Scheme 4 물리 검토: Hanle / EIA / NMOR

## 1. 오늘의 선택과 현재 다섯 scheme 순서

- 현지 날짜: 2026-08-08 (Asia/Seoul)
- 계산: `n = (8 mod 5) + 1 = 4`
- 선택: **Scheme 4, `magneto` — Magneto-optics (Hanle/MOR)**

현재 드롭다운 순서는 실제 registry인 `gabes/schemes/__init__.py:19-25`에서 다음과 같다.

| 순서 | registry 이름 | 사용자-facing scheme | 주 출력 |
|---:|---|---|---|
| 1 | `sas` | OD / SAS | Doppler OD와 hyperfine-pumping SAS |
| 2 | `lambda` | Lambda coherence (EIT / AT / CPT) | 3준위 투과·분산·군지수 |
| 3 | `rydberg_eit` | Rydberg-EIT electrometry | cascade EIT, microwave AT, finite-IF electrometry |
| 4 | `magneto` | Magneto-optics (Hanle/MOR) | Hanle/EIA transmission, NMOR rotation |
| 5 | `fwm` | Four-wave mixing (Squeezing / Biphoton) | seeded gain/squeezing, SFWM biphoton |

README의 사용자 설명도 같은 다섯 항목을 같은 순서로 둔다 (`README.md:8-16`).

## 2. 선행 문서·제안·TODO·issue 검색

새 제안을 만들기 전에 다음을 먼저 확인했다.

- 구현: `gabes/schemes/magneto.py:1-746`
- Zeeman/CG/TOC: `gabes/zeeman.py:17-140`
- 실수기저 Numba kernel: `gabes/kernels.py:276-410`
- 물리 감수율: `gabes/observables.py:393-413`
- 사용자 문서: `README.md:8-16`, `34-56`, `118-135`와
  `docs/Userguide/GABES_User_Guide_v2.html:520-545`, `608-619`
- 핵심 테스트: `tests/test_magneto.py:26-256`, `tests/test_kernels.py:62-99`,
  `tests/test_resonant_hanle_reference.py:14-131`
- 문헌 대조 실행: `tests/verify_hanle_eit_eia.py:1-277`
- 실험 예제: `analysis/squeezing/resonant_hanle_squeezing_reference.py:352-412`,
  `469-626`, `657-727`, `829-930`, `1006-1158`
- 예제 설정: `analysis/squeezing/resonant_hanle_experiment_config.example.json:1-47`
- 현재 계획: `docs/checklist.json:103-138`, `231-319`, `460-479`
- 이전 Scheme 4 보고서: 2026-06-23, 06-28, 07-03, 07-13, 07-23, 07-28,
  08-03의 `docs/daily_report/*scheme-4_hanle-eia-nmor.md`

별도 로컬 issue 노트는 찾지 못했다. 현재 작업 목록은 `docs/checklist.json`과 일일 보고서에
모여 있다. 공개 [GitHub Issues](https://github.com/Shake2313/fwm-squeezing-app/issues)도
2026-08-08 확인 시 `Issues 0`이다.

`magneto.py`와 `tests/test_magneto.py`의 마지막 commit은 모두
`84aba45a`(2026-07-20, `Refine scheme metric hierarchy`)이며, 2026-08-03 검토 뒤 Magneto
소스·테스트·예제에 별도 commit이나 working-tree diff가 없다. 따라서 이번 수치 변화는 새 코드가
아니라 더 엄격한 격자와 읽기 의미를 점검한 결과다.

## 3. 현재 구현이 담은 실제 원자물리

이 scheme은 임의 Lorentzian을 그리는 toy model은 아니다.

1. 선택한 87Rb D1 `Fg -> Fe`의 모든 `mF` 상태를 만들고, `q=-1,0,+1`의
   Clebsch–Gordan 결합을 계산한다 (`gabes/zeeman.py:89-118`).
2. 자발방출을 편광별 `Sigma_q` jump operator로 묶어 excited-state coherence가 ground-state
   coherence로 전달되는 transfer of coherence(TOC)를 보존한다
   (`gabes/zeeman.py:96-118`). 테스트는 자연방출 합과 TOC matrix element를 직접 고정한다
   (`tests/test_magneto.py:35-61`).
3. QWP 각도를 복소 `sigma+/-` drive amplitude로 변환하고
   (`gabes/schemes/magneto.py:113-124`), longitudinal scan과 residual transverse field를
   ground/excited Zeeman Hamiltonian에 넣는다 (`gabes/schemes/magneto.py:431-479`,
   `512-530`).
4. 파라핀 셀은 illuminated/dark density matrix의 교환 OBE로 wall-preserved coherence와
   Ramsey narrowing을 표현한다 (`gabes/schemes/magneto.py:419-429`, `481-488`,
   `533-576`). Buffer cell은 단일 영역 모델이다 (`gabes/schemes/magneto.py:411-418`,
   `490-491`, `578-601`).
5. transmission은 `T=exp(-alpha L)`, NMOR은
   `theta = kL Re(chi_+ - chi_-)/4`로 계산한다
   (`gabes/schemes/magneto.py:642-653`). 따라서 NMOR은 transmission에 이름만 바꾼 값이 아니다.
6. 별도 Hanle/FWM 예제는 측정 CSV가 있으면
   `signal = offset + scale * model(B_scale*B_measured + B_offset)`을 피팅한다
   (`analysis/squeezing/resonant_hanle_squeezing_reference.py:536-626`). 이는 compact model을
   실측 trace의 초기 forward model로 쓰는 좋은 실험 workflow다.

87Rb D1 F=2 -> F'=1 파라핀 셀에서 선형편광 CPT와 원편광 MIA가 잔류 횡장·벽면 결맞음에
따라 전환된다는 정성적 목표는 실제 문헌과 맞는다. Lee와 Moon은 0.12 mG CPT와 0.20 mG
흡수선을 보고했다
([JOSA B 30, 2301 (2013)](https://opg.optica.org/josab/abstract.cfm?uri=josab-30-8-2301)).
Buffer-cell circular LCA의 기준은 2.4 mG이다
([Phys. Rev. A 81, 023416 (2010)](https://doi.org/10.1103/PhysRevA.81.023416)).

## 4. 실험물리 연구자 관점의 적합성

### 유용한 용도

- 전이, 편광, 잔류 횡장, wall lifetime, ground relaxation을 바꿀 때 dip/peak 부호와 폭의
  방향성을 확인하는 실험 설계
- TOC 기반 intrinsic EIA, circular-light LCA, broad transit pedestal와 narrow Ramsey core,
  NMOR 영점 교차의 sanity check
- shielding/coil offset과 effective relaxation의 초기값 탐색
- 측정 trace를 넣은 뒤 B축·signal scale을 보정하는 반정량 forward model
- 실수기저 kernel과 headless observables를 이용한 parameter sweep

### 부적합한 용도

- 실측 보정과 수렴 검증 없이 절대 linewidth, contrast, NMOR slope를 인용하는 것
- 현재 scheme 단독으로 절대 `pT/sqrt(Hz)` 감도를 예측하는 것
- Ne 압력 하나로 pressure shift, diffusion, spin destruction, Dicke narrowing을 정량 예측하는 것
- 다른 ground hyperfine manifold로의 leakage와 repump가 중요한 open transition의 절대 population
- 임의 타원·순환편광에서 현재 scalar를 그대로 “편광면 회전”으로 해석하는 것

따라서 현재 평가는 **실제 물리가 유용한 정성적·반정량적 실험 reference**이다. 다만
**absolute magnetometry 또는 absolute line-shape reference는 아니다**. 내부 일관성 테스트는 강하지만,
held-out measured Hanle trace로 absolute scale과 linewidth를 검증한 상태는 아니다.

## 5. 기존 개선안의 계산비용과 물리 보존 평가

현재 체크리스트는 과거의 넓은 제안을 더 엄격한 항목으로 재정리했다.

| 기존 항목 | 현재 상태 | 계산비용 | 물리를 거의 그대로 보존할 수 있는가 | 판단 |
|---|---|---|---|---|
| `magneto-dynamics-semantics` (`docs/checklist.json:231-249`) | P1, ready, medium | Liouville 차원·`B x v` solve 수는 동일. 자연방출, elastic optical dephasing, transit reload, depolarization의 고정 matrix 항만 분리하므로 runtime overhead는 거의 0 | **예.** zero-rate/zero-pressure limit, trace, Hermiticity, positivity를 고정하면 intended physics를 더 정확히 보존한다. 기존 숫자는 바뀌어야 하지만 이는 오류 의미의 수정이다. | 가장 먼저 구현할 물리 수정 |
| `magneto-light-region-normalization` (`docs/checklist.json:251-270`) | P1, ready, medium | solve 뒤 trace와 나눗셈: `O(n_B n_v n_level^2)`. 새 solve 없음 | **예.** `N`을 local density로 쓰기로 한 decision과 conditional trace-one state를 일치시킨다. 부호·폭 trend는 회귀로 보존하고 absolute scale은 의도적으로 수정한다. | absolute contrast/rotation 전에 필수 |
| `magneto-observable-convergence-and-trust` (`docs/checklist.json:272-291`) | P1, ready, medium | physical-zero 재중심, 상태표시, local fit은 `O(n_B)`로 사실상 무료. 그러나 641-class reference solve는 무겁다 | **분리해야 한다.** 저비용 trust/readout 수정은 즉시 가능하고, 고비용 reference mode는 opt-in이어야 한다. | UI 기본은 상태표시, 배치 검증만 고해상도 |
| `collisional-coefficient-provenance-and-pressure-shift` (`docs/checklist.json:103-138`) | P1, ready, large | pressure shift는 detuning scalar, elastic dephasing도 같은 matrix 차원. runtime overhead는 거의 0이나 자료 조사·API 작업은 큼 | **예.** gas/species/line, 단위·부호·불확도·유효범위를 명시하고 unsupported 조합을 막으면 된다. | 계산비용보다 provenance가 어려움 |
| `geometry-aware-buffer-relaxation-budget` (`docs/checklist.json:294-319`) | P2, needs_decision, large | 선택한 analytic geometry에서 scalar rate budget을 선계산하면 solve overhead는 작음 | **조건부.** pressure만으로 자동값을 만들면 물리를 훼손한다. beam radius, cell geometry, wall 조건, diffusion/spin-destruction source가 먼저 필요하다. | 현재처럼 decision-blocked가 타당 |
| `full-velocity-changing-collision-kernel` (`docs/checklist.json:460-479`) | P2, parked, research | velocity classes가 결합되어 독립 batch solve를 깨고 메모리·시간이 크게 증가 | **아니오.** negligible-overhead 범주가 아니다. kernel·dataset·runtime budget이 정해진 opt-in reference에만 적합하다. | 기본 interactive 경로에 넣지 말 것 |

### 현재 코드에서 재확인한 rate 의미

Buffer mode는
`gamma_g = 2*pi*(buffer_ground_relax_khz + collisional_depol_khz)`를 만든 뒤 같은 값을
`gamma_gg`와 `transit_rate`에 모두 넘긴다 (`gabes/schemes/magneto.py:411-416`). 이번 실행에서
`20/2 kHz`와 `2/20 kHz`를 바꾸어도 `chi_probe`, `chi_p`, `chi_m`의 최대 차이는 모두 정확히
`0.0`이었다. 즉 UI의 두 knob는 현재 물리적으로 구별되지 않는다.

또한 20 Torr Ne의 `78.2 MHz` FWHM broadening이 자연폭 `5.75 MHz`와 더해져
`gamma_opt/2pi = 83.95 MHz`, 자연값의 14.6배가 되고, 이 전체 값이 `Sigma_q` 자발방출
operator의 `gamma`로 들어간다 (`gabes/schemes/magneto.py:407-415`,
`gabes/zeeman.py:96-118`; coefficient는 `gabes/constants.py:79-89`). Elastic pressure
broadening이 population lifetime과 TOC 공급률까지 14.6배로 바꾸는 의미이므로 P1 분리는 정당하다.

## 6. 현재 수치 신뢰 경계

### 6.1 기본 다섯 regime

현재 기본 grid(`scan_points=121`, `velocity_classes=9`, Doppler on)의 대표값이다.

| regime | absorption feature [m^-1] | central width [µT] | T(command B=0) | NMOR slope [mrad/µT] |
|---|---:|---:|---:|---:|
| EIT dip | -0.020047 | 0.16661 | 0.998251 | -1.0191 |
| EIA peak | +0.072676 | 0.22316 | 0.999270 | -0.1976 |
| Buffer Hanle | -4.04604 | 7.58078 | 0.756499 | -2.3356 |
| Buffer LCA | +0.005224 | 2.27758 | 0.655267 | -0.00374 |
| NMOR | -0.081280 | 0.13212 | 0.998137 | -2.55246 |

마지막 NMOR 행 외의 slope는 주 readout이 아니며 진단값일 뿐이다. 이 표의 절대값은 아래 격자와
정규화 문제 때문에 검증된 reference 값으로 해석하면 안 된다.

### 6.2 Doppler/B-grid 수렴과 실제 비용

`velocity_classes` UI 범위는 1-41이고 `scan_points`는 최대 401이다
(`gabes/schemes/magneto.py:291-298`). 체크리스트가 요구하는 641-class 기준은 headless에서만
직접 넣어 실행했다.

| regime | grid `(B,v)` | feature [m^-1] | central width [µT] | NMOR slope [mrad/µT] | compute |
|---|---:|---:|---:|---:|---:|
| EIT | 121 x 9 | -0.020047 | 0.16661 | -1.0191 | warm median 약 0.142 s |
| EIT | 401 x 321 | +0.003545 | 1.08273 | -1.4197 | 13.81 s |
| EIT | 401 x 641 | +0.004031 | 1.02444 | -1.4512 | 21.40 s |
| NMOR | 121 x 9 | -0.081280 | 0.13212 | -2.5525 | 약 0.14 s |
| NMOR | 401 x 321 | -0.054691 | 0.11057 | -6.0421 | 14.02 s |
| NMOR | 401 x 641 | -0.055384 | 0.11012 | -6.1743 | 22.01 s |

핵심은 단순한 몇 퍼센트 오차가 아니다.

- 기본 EIT는 641 classes에서 **dip에서 peak로 부호가 뒤집힌다**.
- 321 -> 641에서도 EIT 중앙폭은 약 5.7% 변하여 체크리스트의 5% 기준을 아직 넘는다.
- 기본 NMOR slope는 641-class 값보다 약 59% 작다. 321 -> 641은 약 2.2%로 수렴하지만
  기본 9 classes는 그 기준을 대표하지 못한다.
- 401 x 641의 반환 `rho(B,v,8,8)`만 약 **251.0 MiB**이며, 각 실행은 약 22초였다.
  기본 warmed solve 대비 약 150배이므로 default UI에 그대로 넣는 것은 부적절하다.

따라서 `coarse/nonconverged` 상태표시는 solve overhead 없이 즉시 추가할 수 있지만,
641-class schedule은 opt-in reference 또는 CI의 느린 검증 경로여야 한다.

### 6.3 파라핀 light-region normalization

현재 two-region 식은 `Tr(rho_light + rho_dark)=1`을 강제한다
(`gabes/schemes/magneto.py:562-575`, `gabes/kernels.py:364-375`). 기본
`gamma_out/2pi=80 kHz`, `gamma_in/2pi=1 kHz`에서 이번 capture는 모든 `(B,v)` 점에 대해

`Tr(rho_light) = gamma_in/(gamma_in+gamma_out) = 1/81`

을 재확인했다. 이 unnormalised light block을 local vapor density `N_eff`와 곱한다
(`gabes/schemes/magneto.py:493-505`, `604-618`, `642-653`;
`gabes/observables.py:393-407`). 결과는 다음과 같다.

| convention | alpha(B=0) [m^-1] | T(B=0) | transmission feature |
|---|---:|---:|---:|
| 현재 unnormalised light block | 0.175039 | 0.998251 | 0.0002001 |
| conditional `rho_light/Tr(rho_light)` 진단 | 14.17818 | 0.867811 | 0.0139780 |

정확히 81배의 absolute absorption 차이다. 체크리스트의 “`N`은 local density, susceptibility는
trace-one conditional state” 결정 (`docs/checklist.json:263-268`)이 현재 구현보다 명확하다.
추가 matrix solve가 없으므로 계산비용은 무시할 수 있다.

### 6.4 `b_offset_ut`의 solve와 readout 기준이 다르다

`compute()`는 `b_physical_ut=b_ut+b_offset_ut`를 Hamiltonian에 올바르게 넣는다
(`gabes/schemes/magneto.py:431-479`). 그러나 `observables()`는 commanded axis인 `b_ut`의 0을
중심으로 feature, width, NMOR slope를 추출한다 (`gabes/schemes/magneto.py:642-680`).
`tests/test_magneto.py:170-176`은 두 축 차이만 확인하며 readout 재중심은 고정하지 않는다.

이전 2026-08-03 진단의 0.25 µT offset에서는 command-zero 분류가 `crossover`, 폭이
0.7247 µT였지만 physical-zero feature는 여전히 -0.01927 m^-1, 폭 약 0.2040 µT였다.
`raw["b_physical_ut"]`를 기준으로 interpolation/local fit하는 것은 `O(n_B)` postprocess이며
새 solve가 없다.

### 6.5 문헌 검증 출력의 과장

`python tests/verify_hanle_eit_eia.py`의 현재 출력은 다음과 같다.

- fixed residual field에서 linear CPT dip / circular MIA peak 부호 전환: **PASS**
- paraffin linear/circular FWHM: `1.195/1.739 mG`
- 문헌: `0.12/0.20 mG`, 즉 각각 약 `9.96x/8.70x` 차이
- buffer LCA의 표시된 low-power/low-relaxation 점: `4.722 mG`
- Yu 등의 기준: `2.4 mG`, 약 `1.97x` 차이

그런데 script는 paraffin 결과를 `MATCH`, `same sub-mG order`로 쓰고
(`tests/verify_hanle_eit_eia.py:260-271`), buffer LCA도 현재 표에 2.4 mG가 없는데
“reaches the reference ~2.4 mG”라고 출력한다 (`tests/verify_hanle_eit_eia.py:174-206`).
부호/trend PASS와 절대폭 CHECK/FAIL을 분리하는 기존 P1은 그대로 유효하다.

## 7. 이번에 새로 확인한 NMOR 편광 readout 범위 문제

이 특정 항목은 기존 Scheme 4 보고서와 현재 checklist에서 찾지 못했다.

UI는 NMOR readout에서도 QWP를 0-45 deg 전 범위로 허용한다
(`gabes/schemes/magneto.py:200-216`). 그러나 `_coherences()`는 `chi_p`와 `chi_m`을 각각의
circular input amplitude가 아니라 공통 total `Omega`로 나누고
(`gabes/schemes/magneto.py:604-618`), readout은 모든 QWP에 선형입사광의 단순
`Re(chi_+-chi_-)` rotation 식을 그대로 적용한다 (`gabes/schemes/magneto.py:644-653`).
출력 Stokes vector나 ellipse major-axis angle은 계산하지 않는다.

NMOR preset에서 QWP=45 deg로 바꾼 이번 실행은 `|E_+|=0`, `|E_-|=sqrt(2)`인 순환편광인데도
다음 hero를 냈다.

- `Rotation at B=0 = 0.01 mrad`
- `Slope dtheta/dB = -0.10 mrad/µT` (hero)
- `Peak |rotation| = 0.01 mrad` (hero)
- note: “zero crossing near B=0”

순환편광에는 정의된 선형 편광면 major axis가 없으므로 이 값을 그대로 “polarization-plane
rotation”이라 부를 수 없다. 타원편광 NMOR 자체는 실제 현상이지만 ellipticity-dependent
polarimetry가 필요하다
([Matsko et al., *Nonlinear Magneto-Optical Rotation of Elliptically Polarized Light*](https://arxiv.org/abs/physics/0210107)).

조심스러운 저비용 순서는 다음과 같다.

1. 먼저 QWP가 선형편광 기준을 벗어나면 현재 scalar rotation을 `unsupported/diagnostic`으로
   내리고 circular endpoint에서는 hero를 숨긴다. 비용은 `O(1)`, solve 변화 없음.
2. thin/undepleted propagation 근사를 명시할 수 있을 때만 circular-channel normalization과
   Jones/Stokes postprocess를 추가해 ellipse orientation과 ellipticity를 함께 보고한다.
   이는 보통 `O(n_B)`이고 새 atomic solve가 필요 없지만, 강흡수·self-consistent propagation까지
   요구하면 negligible-overhead라고 단정하면 안 된다.
3. `tests/test_magneto.py:214-223`의 QWP=0 antisymmetry 테스트에 더해 circular endpoint의
   “plane angle undefined” status와 타원편광 Stokes fixture를 고정한다.

## 8. behavior를 바꾸지 않는 순수 코드 최적화

1. **Kernel-side Doppler-weighted coherence contraction**

   현재 kernel은 모든 `(B,v)`의 full `rho`를 복원해 Python으로 돌려주고
   (`gabes/kernels.py:399-410`), Python이 필요한 세 coherence를 다시 축약한다
   (`gabes/schemes/magneto.py:493-496`, `604-618`). 401 x 641 기준 반환 배열만 251.0 MiB다.
   kernel이 velocity loop 안에서 동일 순서로 `chi_+`, `chi_-`, `chi_probe` 가중합만 누적하면
   solver와 관측량 정의를 유지하면서 이 allocation과 memory traffic을 거의 없앨 수 있다.
   `tests/test_kernels.py:91-99`와 같은 상대오차 `<1e-9` 동등성 검증을 유지해야 한다.

2. **Affine real-basis assembly를 앞당기기**

   Python은 `C_xy+B*C_z`의 complex light/dark stack을 먼저 만들고
   (`gabes/schemes/magneto.py:469-485`), wrapper가 다시 `U^dagger L U`로 real basis 변환한다
   (`gabes/kernels.py:399-405`). `C_xy`, `C_z`, dark coefficient만 한 번 real basis로 옮기고
   kernel 쪽에서 B-affine 조립하면 complex stack과 batched transform을 피할 수 있다.
   1번보다 구현·동등성 검증 부담은 크지만 물리는 바뀌지 않는다.

3. **Angular-momentum matrix 소형 cache**

   `_hamiltonian()` 호출 때 같은 `Fx,Fy,Fz`를 다시 만든다
   (`gabes/schemes/magneto.py:517-520`, `gabes/zeeman.py:48-66`). `(Fg,Fe)`별 작은 LRU cache는
   안전하지만 heavy solve에 비해 이득은 작으므로 우선순위는 낮다.

현재 Numba 실수기저 path는 warmed default EIT에서 0.186 s, NumPy fallback은 0.463 s로
약 2.49배 빨랐고, 세 susceptibility의 상대차는 `1.4e-14`-`6.4e-14`였다. 이미 구현된
실수기저 kernel, BLAS thread 제한, affine B 조립은 유지할 가치가 크다
(`gabes/kernels.py:310-410`, `gabes/schemes/magneto.py:469-491`).
`zeeman_manifold()`와 Hermitian basis는 이미 cache되어 있으므로 다시 제안하지 않는다
(`gabes/zeeman.py:69-70`, `gabes/core.py:60-93`). Headless observables도 완료 상태다
(`gabes/schemes/magneto.py:132`, `620-746`; `docs/checklist.json:200-204`).

## 9. 권고 우선순위와 최종 판단

1. **P1 absolute-scale correctness:** light-region conditional normalization을 구현하고 occupation을
   별도 진단값으로 노출한다. 새 solve가 없다.
2. **P1 dissipator semantics:** 자연방출/TOC와 Ne elastic dephasing을 분리하고, buffer transit
   reload와 collisional depolarization을 독립 dissipator로 만든다. matrix 크기와 solve 수는 같다.
3. **P1 trust/readout:** physical zero와 command zero를 분리하고, coarse Doppler·open-hyperfine·
   문헌폭 불일치 상태를 hero/verification에 명시한다. NMOR slope는 local odd fit 또는 고차 미분을
   쓴다. 대부분 postprocess라 저비용이다.
4. **P1 NMOR polarization scope:** QWP가 선형편광 기준을 벗어나면 현재 rotation hero를 제한하고,
   이후 Jones/Stokes readout을 검증한다.
5. **P1/P2 reference path:** 321/641-class schedule은 batch/slow test로 유지하되, 기본 UI에서는
   status와 저비용 readout 수정으로 신뢰 경계를 솔직하게 드러낸다.
6. **P2 low-overhead physics:** sourced pressure shift와 elastic coefficient table을 넣고, geometry와
   diffusion source가 정해진 뒤에만 buffer relaxation budget을 자동화한다.
7. **P2 performance:** kernel-side coherence contraction을 먼저, affine real-basis assembly를 다음에 한다.
8. **Parked:** full VCC는 target dataset과 runtime budget이 생길 때까지 interactive 기본 경로에서 제외한다.

결론적으로 `MagnetoScheme`은 full addressed `mF` manifold, CG, TOC, 실제 편광 drive, vector B,
two-region Ramsey exchange와 dispersive NMOR를 갖춘 **실질적인 원자물리 모델**이다. 현상 부호,
parameter sensitivity, 실측 trace fitting의 초기 모델로는 충분히 가치가 있다.

그러나 기본 Doppler grid가 EIT 부호와 NMOR slope를 틀릴 수 있고, light-region normalization은
absolute scale을 81배 바꾸며, buffer rate와 lifetime 의미가 겹친다. `b_offset`과 임의 QWP의
NMOR hero도 측정량 정의와 어긋난다. 따라서 현 버전은 measured calibration과 convergence status를
곁들인 compact semi-quantitative reference로 사용하고, 절대 linewidth·contrast·rotation·자력계
감도의 독립 reference로 인용하지 않는 것이 안전하다.

## 10. 검증 기록

- registry 직접 확인: `sas -> lambda -> rydberg_eit -> magneto -> fwm`
- 공개 GitHub Issues: 0
- read-only 수치 진단:
  - 다섯 built-in regime 대표 출력
  - buffer rate-swap exact equality
  - paraffin light-block trace와 conditional normalization
  - 401 x 321 및 401 x 641 Doppler/B-grid schedule
  - QWP 0/20/40/45 deg NMOR readout
- `python tests/verify_hanle_eit_eia.py`: 부호 전환 PASS, paraffin 폭
  `1.195/1.739 mG` 대 문헌 `0.12/0.20 mG`, buffer 표시점 `4.722 mG` 대 문헌 `2.4 mG`
- targeted:
  `python -m pytest -q tests/test_magneto.py tests/test_kernels.py tests/test_headless_observables.py tests/test_schemes_render.py tests/test_resonant_hanle_reference.py`
  -> **63 passed in 11.56 s**
- 전체 `python -m pytest -q` -> **380 passed in 61.22 s**
- 기존 dirty README/checklist/Rydberg/SABES/분석/보고서 작업은 보존했으며 production code는 수정하지 않았다.
