# 2026-07-28 Scheme 4 물리 검토: Hanle / EIA / NMOR

## 1. 오늘의 선택과 현재 다섯 scheme 순서

Asia/Seoul 현지 날짜의 day는 28이므로

`n = (28 mod 5) + 1 = 4`

이다. 드롭다운 순서는 `gabes/schemes/__init__.py:19-24`의 `_SCHEMES`가 결정한다.

| 순번 | 등록 인스턴스 / 이름 | 현재 표시 제목 | 주 출력 |
|---:|---|---|---|
| 1 | `SASScheme()` / `sas` | Absorption spectroscopy (OD / SAS) | pump-off OD, pump-on SAS |
| 2 | `LambdaScheme()` / `lambda` | Lambda coherence (EIT / AT / CPT) | 3준위 흡수·분산, EIT/AT/CPT |
| 3 | `RydbergEITScheme()` / `rydberg_eit` | Rydberg-EIT electrometry | cascade EIT, microwave AT, finite-IF electrometry |
| 4 | `MagnetoScheme()` / `magneto` | Magneto-optics (Hanle/MOR) | Hanle/EIA transmission과 NMOR rotation |
| 5 | `FWMScheme()` / `fwm` | Four-wave mixing (Squeezing / Biphoton) | seeded gain/squeezing과 SFWM biphoton |

README의 사용자용 목록도 같은 순서다 (`README.md:8-17`). 따라서 오늘 검토 대상은
`MagnetoScheme`이다 (`gabes/schemes/magneto.py:127-173`).

## 2. 선행 제안·문서·테스트·예제·issue 노트 선검색

새 제안을 만들기 전에 다음을 먼저 조사했다.

- 구현: `gabes/schemes/magneto.py:1-746`
- Zeeman/CG/TOC 원자 모델: `gabes/zeeman.py:17-140`, `gabes/atoms.py:21-98`
- 수치 커널과 공통 관측량: `gabes/kernels.py:276-410`,
  `gabes/core.py:96-121`, `gabes/observables.py:393-410`
- 사용자 문서: `README.md:8-17`, `README.md:34-35`, `README.md:131`,
  `docs/Userguide/GABES_User_Guide_v2.html:590-602`
- scheme 자체 reference 패널: `gabes/schemes/magneto.py:368-388`
- 물리/회귀 테스트: `tests/test_magneto.py:26-256`,
  `tests/test_kernels.py:62-96`, `tests/test_headless_observables.py:38-87`,
  `tests/test_schemes_render.py:118-125`
- 문헌 대조 예제: `tests/verify_hanle_eit_eia.py:1-277`
- 실제 분석 예제:
  `analysis/squeezing/resonant_hanle_squeezing_reference.py:352-412`,
  `536-626`, `900-987`, `1706-1783` 및
  `analysis/squeezing/resonant_hanle_experiment_config.example.json:1-47`
- 계획/TODO: `docs/checklist.json:50-60`, `88-92`, `130-134`,
  `gabes/constants.py:79-89`
- 이전 Scheme 4 보고서: 2026-06-23, 06-28, 07-03, 07-13, 07-23의
  `docs/daily_report/*scheme-4_hanle-eia-nmor.md`

`examples/` 디렉터리는 없지만, resonant-Hanle 분석 스크립트는 FWM probe를 87Rb D1
Hanle cell에 통과시키고 balanced detection, 절대 감도, 손실 sweep, 측정 CSV 보정까지
실행하는 실험 지향 예제다. 이 예제 자체도 모델을 “compact and semi-quantitative”라고
명시한다 (`analysis/squeezing/resonant_hanle_squeezing_reference.py:1142-1150`).

로컬 별도 issue-note 파일은 없었다. GitHub connector로
`Shake2313/fwm-squeezing-app`의 open/closed issue를 2026-07-28에 다시 검색했으며 둘 다
0건이었다. 2026-07-23 검토 이후 magneto 구현·테스트에는 새 commit이나 working-tree
변경이 없고, 분석 문서 위치만 `analysis/squeezing/` 아래로 정리됐다.

### 이미 있던 개선안과 현재 비용 판단

| 기존 항목 | 상태 | 물리 보존과 계산 비용 | 이번 판단 |
|---|---|---|---|
| `buffer-gas-pressure-shift` | deferred | pressure shift는 detuning scalar, 저차 Dicke 보정도 같은 Liouville 차원이라 solve overhead가 거의 없다 | 여전히 비용 대비 가치가 좋다. gas·species·line별 계수 근거가 필요하다 (`docs/checklist.json:50-54`). |
| `magneto-buffer-relaxation-map` | deferred | solve 전 pressure·temperature·beam geometry에서 두 scalar rate를 정하므로 사실상 0 | 아래에서 확인한 rate 의미 중복을 먼저 고친 뒤 구현해야 한다 (`docs/checklist.json:88-92`). |
| `full-velocity-changing-collision-kernel` | deferred, GROUP C | velocity class가 서로 결합되어 현재 separable `(B,v)` solve를 깨뜨린다 | 정밀 buffer-cell fit에는 유용하지만 interactive 기본 경로에는 무겁다 (`docs/checklist.json:130-134`). |
| `figureless-observables-paths` | done | figure 생성을 건너뛰며 물리는 동일하다 | 이미 `supports_headless_observables=True`와 `include_figures=False`로 완료됐다 (`gabes/schemes/magneto.py:132`, `620-746`). |
| `b_offset_ut` | done | B축 scalar 이동만 추가한다 | shielding/coil zero 보정에 유용하고 비용은 0에 가깝다 (`gabes/schemes/magneto.py:220-225`, `431-435`). |

2026-07-23 보고서가 추가로 제안한 Ne optical pure-dephasing 분리, Doppler 수렴
상태 표시, NMOR slope 고차 추출, kernel-side coherence contraction도 아직 구현되지 않았다
(`docs/daily_report/2026-07-23_scheme-4_hanle-eia-nmor.md:121-208`).

## 3. 현재 구현이 담고 있는 실제 원자물리

이 scheme은 단순 Lorentzian toy curve보다 훨씬 유용하다.

1. 선택한 87Rb D1 `Fg -> Fe`의 모든 `mF` 상태를 만들고 `q=-1,0,+1` 전이의
   Clebsch-Gordan 계수를 계산한다 (`gabes/zeeman.py:89-118`).
2. 자발방출을 채널별 독립 jump가 아니라 편광별 `Sigma_q`로 묶어 excited coherence가
   ground coherence로 전달되는 transfer of coherence(TOC)를 보존한다
   (`gabes/zeeman.py:96-118`). 이 때문에 `Fe=Fg+1`의 intrinsic EIA와
   `Fe<=Fg`의 EIT 부호 규칙을 테스트할 수 있다 (`tests/test_magneto.py:179-192`).
3. QWP 각도는 실제 `sigma+/-` 복소 drive amplitude를 바꾸고
   (`gabes/schemes/magneto.py:113-124`), longitudinal scan과 residual transverse field는
   ground/excited Zeeman Hamiltonian에 들어간다 (`gabes/schemes/magneto.py:431-479`,
   `511-530`). 그래서 linear Hanle dip, circular MIA/EIA, transverse-field LCA가
   같은 OBE에서 나온다.
4. paraffin cell은 illuminated/dark 두 영역의 density matrix를 교환하여 wall-preserved
   coherence와 Ramsey narrowing을 표현한다 (`gabes/schemes/magneto.py:419-429`,
   `481-488`, `532-576`). Buffer cell은 단일 영역에 optical broadening과 phenomenological
   ground rates를 쓴다 (`gabes/schemes/magneto.py:407-418`, `578-601`).
5. transmission은 Beer-Lambert `exp(-alpha L)`, NMOR는
   `theta = kL Re(chi_+ - chi_-)/4`로 계산된다 (`gabes/schemes/magneto.py:642-653`).
   NMOR는 transmission curve의 이름만 바꾼 것이 아니라 circular susceptibility의
   dispersive 차이다.
6. 분석 예제는 측정 Hanle CSV가 있으면
   `signal = offset + scale * model(B_scale*B_measured + B_offset)`을 fit한다
   (`analysis/squeezing/resonant_hanle_squeezing_reference.py:527-626`). 이처럼 measured
   trace로 effective parameter를 보정하는 사용법은 현재 모델의 적절한 용도다.

[Lee와 Moon의 paraffin-cell 실험](https://doi.org/10.1364/JOSAB.30.002301)은
87Rb D1 `F=2 -> F'=1`에서 linear-polarization 0.12 mG CPT가 circular-polarization
0.20 mG absorption으로 바뀌는 것을 보고한다. 현재 코드가 polarization sign switch와
two-region 방향성을 재현하는 것은 실제 물리와 맞는다. 그러나 오늘 재실행한
`tests/verify_hanle_eit_eia.py`는 1.195/1.739 mG를 냈다. 각각 문헌값보다 약
10.0/8.7배 넓은데도 마지막에 `MATCH`와 `same sub-mG order`를 출력한다
(`tests/verify_hanle_eit_eia.py:243-271`). 이 문구는 정량 판정으로 부정확하다.

### 현재 대표 출력

kernel warm-up 뒤 현재 기본 regime을 실행한 값이다. 시간은 이 실행 환경의 wall-clock이고
절대 benchmark가 아니라 상대 비용 지표다.

| regime | compute | absorption feature amplitude [m^-1] | central width [uT] | T(B=0) |
|---|---:|---:|---:|---:|
| EIT dip | 0.457 s | -0.02005 | 0.1666 | 0.99825 |
| EIA peak | 0.428 s | +0.07268 | 0.2232 | 0.99927 |
| Buffer Hanle | 0.349 s | -4.0460 | 7.5808 | 0.75650 |
| Buffer LCA | 0.418 s | +0.005224 | 2.2776 | 0.65527 |
| NMOR | 0.479 s | -0.08128 | 0.1321 | 0.99814 |

부호·편광·residual-field·wall-lifetime의 방향성 연구에는 충분히 의미가 있다. 반면 아래
신뢰 경계 때문에 absolute linewidth, absolute contrast, absolute NMOR slope나
magnetometer sensitivity의 독립 예측 reference로는 아직 적합하지 않다.

## 4. 이번 실행에서 새로 확인한 신뢰 경계

### 4.1 두 buffer-ground rate knob가 수학적으로 완전히 같은 knob다

UI 설명은 `buffer_ground_relax_khz`를 ground-state relaxation, `collisional_depol_khz`를
spin-destruction ground-coherence loss로 구분한다
(`gabes/schemes/magneto.py:251-262`). 그러나 compute는

```python
gamma_g = 2*pi*(buffer_ground_relax_khz + collisional_depol_khz)
zeeman_manifold(..., gamma_gg=gamma_g, transit_rate=gamma_g)
```

처럼 합계 하나만 만든다 (`gabes/schemes/magneto.py:411-415`). 따라서
`(20,2) kHz`와 `(2,20) kHz`를 서로 바꾼 두 실행의 `chi_probe`, `chi_p`, `chi_m`
최대 차이는 모두 정확히 `0.0`이었다. 두 control은 현재 독립적인 물리 control이 아니다.

더구나 `transit_rate`는 모든 source state에서 ground manifold로 population reload
Lindblad channel을 만들고 (`gabes/zeeman.py:120-126`), `gamma_gg`는 ground off-diagonal에
별도 dephasing을 다시 더한다 (`gabes/zeeman.py:127-134`). ground coherence 한 개에 대한
현재 Lindblad diagonal을 계산하면 `-2*gamma_g`였다. 즉 합산 rate는 population reload로
한 번, 명시적 dephasing으로 한 번 들어간다.

이 문제는 solver 차원이나 grid를 늘리지 않고 고칠 수 있다.

- `gamma_transit = 2*pi*buffer_ground_relax_khz`는 `transit_rate`에만 넣는다.
  이 channel 자체가 ground coherence를 `gamma_transit`만큼 감쇠시킨다.
- `gamma_depol = 2*pi*collisional_depol_khz`는 추가 `gamma_gg` 또는, 더 물리적으로,
  같은 `mF` manifold 안의 spin-randomization Lindblad map으로 분리한다.
- `(20,2)`와 `(2,20)`이 다른 population/contrast를 내는 테스트, laser-off coherence
  decay가 선언한 총 rate와 같은 테스트를 추가한다.

행렬 차원, `(B,v)` 점 수, solve 횟수는 모두 그대로이므로 overhead는 사실상 없다.
계수의 실험적 의미를 바로잡는 것이 pressure-to-relaxation 자동 매핑보다 먼저다.

### 4.2 다른 ground hyperfine manifold로의 자발방출 leakage가 사라진다

현재 `zeeman_manifold(Fg,Fe)`는 선택한 ground `Fg` 하나만 만들고
(`gabes/zeeman.py:89-94`), 각 excited `mF`의 총 자연방출 `Gamma`를 모두 그 same-`Fg`
manifold로 되돌린다 (`gabes/zeeman.py:96-118`). 테스트도 각 excited state의
`sum Sigma_q^dagger Sigma_q = Gamma`가 현재 manifold 안에서 성립하도록 고정한다
(`tests/test_magneto.py:35-45`).

그러나 87Rb D1의 각 `F'`는 `F=1`과 `F=2` 양쪽으로 decay할 수 있다.
[Steck의 87Rb D-line data, Table 8](https://steck.us/alkalidata/rubidium87numbers.pdf)의
`S_FF'`와 저장소 자체의 full-hyperfine branching 규약
`T=(2Fg+1)S_FF'` (`gabes/species.py:364-382`)을 쓰면 귀환율은 다음과 같다.

| addressed transition | addressed `Fg`로 귀환 | 다른 ground `F`로 leakage | repump 없이 5회 산란 뒤 bright-manifold 잔류 확률 |
|---|---:|---:|---:|
| `F=2 -> F'=1` (기본 Hanle) | 5/6 | 1/6 | 0.402 |
| `F=1 -> F'=2` (intrinsic EIA test) | 1/2 | 1/2 | 0.0313 |
| `F=2 -> F'=2` (LCA/분석 예제) | 1/2 | 1/2 | 0.0313 |
| `F=1 -> F'=1` | 1/6 | 5/6 | 1.29e-4 |

실제 잔류 population은 transit, wall return, repump, detuning, saturation과 경쟁하므로 이
마지막 열을 그대로 실험 contrast로 읽으면 안 된다. 하지만 코드가 모든 경우 귀환율 1을
강제한다는 사실은 명확하다. 따라서 현재 intrinsic-EIA test는 TOC 부호를 검증하지만,
open-transition optical pumping과 steady-state contrast까지 검증하지는 않는다.

가장 싼 개선은 zero-solve-cost warning이다. 각 transition에 physical branching ratio와
“다른 hyperfine ground state/repump 미포함” 상태를 Derived table과 hero metric에 표시하면
결과의 용도를 오해하지 않는다.

정확한 수정은 negligible-overhead가 아니다.

- 기본 `F=2 -> F'=1`에 다른 ground manifold까지 모두 넣으면 level 수가
  `8 -> 11`, paraffin two-region dense system이 `128 -> 242`가 된다. 단순 cubic
  factorization 기준 비용은 약 `(242/128)^3 = 6.8`배까지 늘 수 있다.
- dark reservoir 한 level과 measured repump/transit rate를 쓰는 compact model은
  default system을 `128 -> 162`로 늘려 대략 2배 dense-solve 비용이 예상되지만,
  full Zeeman ground manifold보다 싸다.
- 동일 차원의 사후 bright-fraction rate proxy는 매우 싸지만 coherence와 optical pumping의
  feedback을 정확히 보존하지 못하므로 “approximate” reference mode로만 허용해야 한다.

### 4.3 이전 보고서의 P0 경계도 그대로 남아 있다

1. **Ne optical broadening 의미:** `gamma_opt = GAMMA_D1 + buffer_gamma`를
   `zeeman_manifold(..., gamma=gamma_opt)`에 넘기므로
   (`gabes/schemes/magneto.py:407-415`), pressure broadening이 optical pure dephasing이
   아니라 excited population decay와 TOC refeeding까지 가속한다
   (`gabes/zeeman.py:96-118`). 자연방출 `Gamma`와 elastic optical dephasing을 분리하는
   수정은 같은 차원이라 overhead가 거의 없다.
2. **Doppler 수렴:** 기본 EIT preset에서 velocity class를 9에서 641로 늘리면
   `doppler_scale 0.0560 -> 0.9972`, feature amplitude
   `-0.02005 -> +0.004031 m^-1`, width `0.1666 -> 1.0253 uT`로 변했다. 부호까지
   바뀐다. 이 실행의 compute time은 `0.398 -> 22.692 s`였다. 현재 scalar
   `_doppler_dilution()`은 linear Voigt 크기만 맞추며 nonlinear optical pumping/Zeeman
   line shape를 수렴시키지 못한다 (`gabes/schemes/magneto.py:68-100`, `460-496`).
3. **NMOR B-grid 수렴:** 같은 기본 NMOR에서 121/201/401점 중앙미분 slope는 각각
   `-2.5525/-2.9511/-3.1503 mrad/uT`였다. hero slope는 여전히
   `np.gradient()` 한 번으로 계산된다 (`gabes/schemes/magneto.py:671-680`).

Doppler는 adaptive/composite quadrature의 reference mode가 필요해 계산이 비싸지만, 먼저
`coarse/qualitative only` 상태를 표시하는 것은 비용 0이다. NMOR는 local odd-polynomial fit,
5-point derivative, samples-per-width 상태를 같은 grid에 적용할 수 있어 거의 무료다.

## 5. 순수 코드 최적화

물리와 출력을 바꾸지 않는 우선순위는 다음과 같다.

1. **kernel 안에서 필요한 coherence만 Doppler-weighted contraction한다.** 현재 real-basis
   kernel이 모든 `(B,v)`의 full `rho`를 complex 행렬로 복원한 뒤
   (`gabes/kernels.py:399-410`), Python이 `chi_+`, `chi_-`, `chi_probe`를 다시 순회·평균한다
   (`gabes/schemes/magneto.py:493-496`, `604-618`). default `F=2 -> F'=1`,
   `121 x 641` reference grid의 complex `rho`만 약 75.7 MiB다. real solution에서 필요한
   두 polarization coherence를 직접 누적하면 결과는 같고 이 materialization과 메모리
   traffic을 제거할 수 있다. Numba/NumPy parity test는 유지해야 한다.
2. **affine real-basis coefficient를 kernel에 전달한다.** 현재 Python이
   `C_xy + B*C_z` complex stack을 만든 뒤 (`gabes/schemes/magneto.py:469-485`),
   kernel wrapper가 다시 `U^dagger L U`를 batched 수행한다 (`gabes/kernels.py:399-405`).
   `C_xy`, `C_z`, dark coefficient를 한 번만 real basis로 바꾸고 kernel에서 B-affine
   조립하면 큰 complex stack과 transform을 피할 수 있다. 물리는 동일하지만 1번보다
   구현·동등성 검증 비용이 크다.
3. **angular-momentum matrix의 작은 read-only cache는 안전하지만 미세 개선이다.**
   `_hamiltonian()` 호출마다 `angular_momentum_matrices()`를 다시 만든다
   (`gabes/schemes/magneto.py:511-520`, `gabes/zeeman.py:48-66`). 이전 측정은 호출당
   수십 us 수준이므로 heavy solve보다 우선순위가 낮다.
4. `zeeman_manifold()`와 Hermitian basis는 이미 LRU cache다
   (`gabes/zeeman.py:69-70`, `gabes/core.py:60-93`). 같은 template cache를 다시
   제안할 이유가 없다. Headless observables도 이미 완료됐다.

## 6. 권고 우선순위와 최종 판단

1. **P0 rate 의미 수정:** buffer ground relaxation과 collisional depolarization을 분리하고
   현행 `-2*(rate sum)` coherence decay를 명시적 rate budget으로 고친다. 차원·solve 수
   변화가 없어 가장 싼 correctness 개선이다.
2. **P0 optical dissipator 수정:** 자연방출/TOC `Gamma`와 Ne optical pure dephasing을
   분리한다. 이 역시 거의 zero-overhead다.
3. **P0 신뢰 경계 표시:** Doppler quadrature 상태와 open-hyperfine branching/repump
   미포함 상태를 hero/Derived readout에 표시한다. 추가 solve는 필요 없다.
4. **P1 reference 정확도:** adaptive Doppler reference mode와 수렴 검사를 추가하고,
   NMOR slope는 local fit 또는 고차 중앙미분으로 추출한다. 전자는 고비용, 후자는 거의
   무료다.
5. **P1 기존 저부하 물리:** pressure shift, 저차 Dicke correction, 그리고 의미가 분리된
   pressure-to-relaxation mapping을 문헌계수와 override를 보존한 채 구현한다.
6. **P1/P2 open-transition model:** measured repump가 있는 compact reservoir 또는
   두 ground hyperfine Zeeman manifold를 opt-in reference fidelity로 제공한다. full model은
   interactive 기본값으로 쓰기에는 수 배 비싸다.
7. **P2 성능:** kernel-side coherence contraction을 먼저 하고 affine real-basis assembly를
   뒤따르게 한다. 이는 reference Doppler refinement의 메모리·시간 비용을 줄인다.

결론적으로 `MagnetoScheme`은 Zeeman manifold, CG, TOC, 실제 편광 drive, vector B,
two-region Ramsey exchange, dispersive NMOR까지 갖춘 **실제 물리가 유용한 정성적·반정량적
실험 reference**다. dip/peak/zero-crossing의 방향, QWP·residual-field·wall-lifetime 감도,
측정 trace fit의 초기 forward model로는 충분히 가치가 있다.

그러나 현재 두 buffer rate control은 사실상 하나이고 coherence rate가 중복 적용되며, 모든
자연방출을 addressed hyperfine ground manifold로 되돌린다. 기본 Doppler grid는 feature
부호도 수렴시키지 못하고, 문헌 검증 script의 linewidth MATCH 문구도 실제 숫자와 맞지 않는다.
따라서 현 버전을 absolute linewidth, absolute contrast, absolute NMOR slope 또는
`pT/sqrt(Hz)` 자력계 성능의 독립 예측 reference로 사용하면 안 된다. measured trace로
B축·signal scale·relaxation/repump를 보정하고 수렴 상태를 확인하는 compact model로 쓰는
것이 적절하다.

## 7. 검증 기록

- 현재 registry를 직접 출력해 `sas -> lambda -> rydberg_eit -> magneto -> fwm` 순서를 확인했다.
- `tests/verify_hanle_eit_eia.py`를 재실행했다. polarization sign switch와 TOC 부호는
  재현됐지만 1.195/1.739 mG 대 문헌 0.12/0.20 mG의 불일치도 재확인했다.
- buffer-rate swap, Lindblad coherence diagonal, hyperfine branching, Doppler-class sweep,
  NMOR B-grid sweep는 production 파일을 수정하지 않은 read-only 진단으로 수행했다.
- `python -m pytest tests/test_magneto.py tests/test_kernels.py
  tests/test_headless_observables.py tests/test_schemes_render.py
  tests/test_resonant_hanle_reference.py -q`: **63 passed in 32.80 s**.
- AGENTS.md가 요구하는 전체 `python -m pytest -q`:
  **232 passed in 100.37 s**.
