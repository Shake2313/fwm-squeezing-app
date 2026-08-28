# 2026-08-03 Scheme 4 물리 검토: Hanle / EIA / NMOR

## 1. 오늘의 선택과 현재 다섯 scheme 순서

Asia/Seoul 현지 날짜의 day는 3이므로

`n = (3 mod 5) + 1 = 4`

이다. 실제 dropdown 순서는 `gabes/schemes/__init__.py:19-24`의 `_SCHEMES`가 정한다.

| 순번 | 등록 이름 | 표시 제목 | 핵심 출력 |
|---:|---|---|---|
| 1 | `sas` | Absorption spectroscopy (OD / SAS) | pump-off OD, pump-on SAS |
| 2 | `lambda` | Λ coherence (EIT / AT / CPT) | EIT/AT/CPT 투과·분산 |
| 3 | `rydberg_eit` | Rydberg-EIT electrometry | cascade EIT, microwave AT, finite-IF electrometry |
| 4 | `magneto` | Magneto-optics (Hanle/MOR) | Hanle/EIA transmission, NMOR rotation |
| 5 | `fwm` | Four-wave mixing (Squeezing / Biphoton) | seeded gain/squeezing, SFWM biphoton |

따라서 오늘의 대상은 `MagnetoScheme`이다 (`gabes/schemes/magneto.py:127-173`).

## 2. 선행 문서·제안·테스트·예제·issue 선검색

새 개선안을 만들기 전에 다음을 확인했다.

- 현재 구현: `gabes/schemes/magneto.py:1-746`
- Zeeman/CG/TOC 원자 모델: `gabes/zeeman.py:17-140`, `gabes/atoms.py:21-98`
- 실수기저 고속 커널: `gabes/kernels.py:276-410`
- 사용자 문서: `README.md:8-17`, `README.md:34-35`, `README.md:131`,
  `docs/Userguide/GABES_User_Guide_v2.html:538-541`, `608-619`, `791-805`
- 물리·회귀 테스트: `tests/test_magneto.py:26-256`, `tests/test_kernels.py:62-99`,
  `tests/test_resonant_hanle_reference.py:14-131`
- 문헌 대조: `tests/verify_hanle_eit_eia.py:1-277`
- 실험 분석 예제: `analysis/squeezing/resonant_hanle_squeezing_reference.py:352-412`,
  `527-626`, `900-987`, `1138-1150` 및
  `analysis/squeezing/resonant_hanle_experiment_config.example.json:1-47`
- 계획/TODO: `docs/checklist.json:50-61`, `88-92`, `130-134`,
  `gabes/constants.py:79-89`
- 이전 Scheme 4 보고서: 2026-06-23, 06-28, 07-03, 07-13, 07-23, 07-28의
  `docs/daily_report/*scheme-4_hanle-eia-nmor.md`

별도 로컬 issue/TODO 파일은 없고, 계획은 `docs/checklist.json`과 이전 일일 보고서에 모여 있다.
공개 GitHub Issues도 2026-08-03 현재 open/closed 모두 0건이다
([repository issues](https://github.com/Shake2313/fwm-squeezing-app/issues)).
2026-07-28 검토 뒤 Magneto 구현·테스트의 commit 또는 working-tree 변경은 없었다.

### 이미 있던 개선안과 계산비용 판단

| 기존 제안 | 상태 | 계산비용과 물리 보존 평가 | 오늘의 판단 |
|---|---|---|---|
| `buffer-gas-pressure-shift` | deferred | pressure shift는 detuning scalar, 저차 Dicke 보정도 같은 Liouville 차원에서 가능하므로 solve overhead가 거의 없다 (`docs/checklist.json:50-54`) | 여전히 비용 대비 가치가 높다. gas/species/line별 실험계수와 적용 범위를 명시해야 한다. |
| `magneto-buffer-relaxation-map` | deferred | pressure·temperature·beam geometry에서 기본 scalar rate 두 개를 solve 전에 정하므로 사실상 무료다 (`docs/checklist.json:88-92`) | 아래의 rate 의미 중복을 먼저 고친 뒤 measured override를 보존해 구현해야 한다. |
| `full-velocity-changing-collision-kernel` | deferred, heavy | velocity class가 서로 결합되어 현재 독립 `(B,v)` solve를 깨므로 negligible-overhead가 아니다 (`docs/checklist.json:130-134`) | 정밀 buffer-cell fit에는 유용하지만 interactive 기본 경로에는 부적합하다. |
| headless observables | done | figure만 생략하고 물리는 그대로다 (`gabes/schemes/magneto.py:132`, `620-746`) | 이미 완료됐으므로 다시 제안할 필요가 없다. |
| `b_offset_ut` | done | Hamiltonian의 B scalar를 이동시키는 비용만 든다 (`gabes/schemes/magneto.py:220-225`, `431-435`) | 물리 solve는 구현됐지만, 아래 5.1의 readout 중심 버그가 남았다. |

이전 보고서가 제안한 자연방출/Ne pure-dephasing 분리, buffer rate 의미 분리,
Doppler 수렴 상태, NMOR 고차 기울기, open-hyperfine 경고·reservoir, kernel-side
coherence contraction도 아직 구현되지 않았다
(`docs/daily_report/2026-07-23_scheme-4_hanle-eia-nmor.md:121-247`,
`docs/daily_report/2026-07-28_scheme-4_hanle-eia-nmor.md:123-263`).

## 3. 현재 구현이 담은 실제 원자물리

이 scheme은 단순 Lorentzian 생성기가 아니다.

1. 선택한 87Rb D1 `Fg -> Fe`의 모든 `mF` 상태와 `q=-1,0,+1` Clebsch–Gordan
   결합을 만든다 (`gabes/zeeman.py:89-118`).
2. 자발방출을 개별 population channel이 아니라 편광별 `Σ_q` jump operator로 묶어
   transfer of coherence(TOC)를 보존한다 (`gabes/zeeman.py:96-118`). 그 결과
   `Fe=Fg+1` intrinsic EIA와 `Fe<=Fg` EIT 부호 규칙을 재현한다
   (`tests/test_magneto.py:179-192`).
3. QWP 각도를 실제 복소 `σ+/-` drive amplitude로 변환하고
   (`gabes/schemes/magneto.py:113-124`), longitudinal scan과 residual transverse field를
   ground/excited Zeeman Hamiltonian에 넣는다 (`gabes/schemes/magneto.py:431-479`,
   `511-530`).
4. 파라핀 셀은 illuminated/dark density matrix의 교환 OBE로 wall-preserved coherence와
   Ramsey narrowing을 표현한다 (`gabes/schemes/magneto.py:419-429`, `481-488`,
   `532-576`). Buffer cell은 단일 영역에서 optical broadening과 phenomenological
   ground rates를 쓴다 (`gabes/schemes/magneto.py:407-418`, `578-601`).
5. transmission은 Beer–Lambert `T=exp(-alpha L)`, NMOR은
   `theta=kL Re(chi_+-chi_-)/4`로 계산한다 (`gabes/schemes/magneto.py:642-653`).
   즉 NMOR은 transmission의 이름만 바꾼 신호가 아니다.
6. 별도 분석 예제는 측정 CSV가 있으면
   `signal = offset + scale*model(B_scale*B_measured+B_offset)`을 피팅한다
   (`analysis/squeezing/resonant_hanle_squeezing_reference.py:527-626`). 현재 모델을
   측정 trace에 보정해 쓰는 실험 workflow와 잘 맞는다.

### 현재 대표 출력

아래는 현재 기본 grid(`scan_points=121`, `velocity_classes=9`, Doppler on)를 kernel warm-up을
포함해 한 번씩 실행한 값이다. 시간은 이 실행 환경의 상대 지표이며 절대 benchmark가 아니다.

| regime | compute | absorption feature [m^-1] | central width [µT] | T(B command=0) |
|---|---:|---:|---:|---:|
| EIT dip | 1.043 s | -0.020047 | 0.16661 | 0.998251 |
| EIA peak | 0.295 s | +0.072676 | 0.22316 | 0.999270 |
| Buffer Hanle | 0.180 s | -4.04604 | 7.58078 | 0.756499 |
| Buffer LCA | 0.270 s | +0.005224 | 2.27758 | 0.655267 |
| NMOR | 0.317 s | -0.081280 | 0.13212 | 0.998137 |

NMOR 기본 hero 중앙미분 기울기는 `-2.55246 mrad/µT`이다. 아래 신뢰 경계 때문에
이 표의 절대 폭·대비·회전 기울기를 검증된 절대값으로 해석해서는 안 된다.

## 4. 실험물리 연구자 관점의 등급

### 유용한 용도

- 전이·편광·잔류 횡장·wall lifetime·ground relaxation을 바꿀 때 dip/peak와 폭의 방향성이
  어떻게 변하는지 조사하는 실험 설계
- 파라핀 셀의 broad transit pedestal와 narrow Ramsey core, TOC 기반 intrinsic EIA,
  circular-light LCA, NMOR zero crossing의 sanity check
- 측정 trace를 넣기 전 coil offset·B scale·signal scale·effective relaxation의 초기값 탐색
- 고속 kernel과 headless 경로를 이용한 parameter sensitivity sweep

### 부적합한 용도

- 측정 보정과 수렴 확인 없이 절대 linewidth, 절대 contrast, 절대 NMOR slope를 인용하는 것
- 현재 scheme 자체만으로 `pT/sqrt(Hz)` magnetometer 성능을 독립 예측하는 것
- Ne 압력만으로 diffusion·spin destruction·pressure shift·Dicke narrowing까지 정량 예측하는 것
- repump가 없는 open hyperfine transition의 steady-state population을 절대값으로 믿는 것

따라서 현재 코드는 **실제 물리가 유용한 정성적·반정량적 실험 reference**이지만,
**absolute metrology reference는 아니다**. 사용자 가이드도 전체 도구를 최종 fitting 엔진보다
parameter map으로 쓰라고 설명한다 (`docs/Userguide/GABES_User_Guide_v2.html:377-403`).

## 5. 이번 검토에서 새로 확인한 신뢰 경계

### 5.1 `b_offset_ut`가 물리 공진은 옮기지만 hero readout은 명령축 0에 고정된다

`compute()`는 `b_physical_ut = b_ut + b_offset_ut`를 Hamiltonian에 올바르게 사용한다
(`gabes/schemes/magneto.py:431-479`). 그러나 `observables()`는 다시 `x=raw["b_ut"]`를
잡고 `argmin(abs(x))`, 즉 **commanded B=0**에서 feature 부호·폭·NMOR slope를 계산한다
(`gabes/schemes/magneto.py:642-672`). 실제 zero field는 `b_physical_ut=0`, 곧
`b_ut=-b_offset_ut`에 있다.

기본 EIT 조건의 read-only 진단은 다음과 같았다.

| offset [µT] | 현재 command-zero amplitude [m^-1] | physical-zero amplitude [m^-1] | 현재 hero width [µT] | physical-zero 중심 width [µT] | 현재 분류 |
|---:|---:|---:|---:|---:|---|
| 0.00 | -0.020047 | -0.020047 | 0.1666 | 0.1666 | EIT-like dip |
| 0.25 | -0.003097 | -0.019266 | 0.7247 | 0.2040 | crossover |
| 0.50 | -0.001815 | -0.020121 | 1.5568 | 0.1673 | crossover |

즉 0.25 µT offset만으로 실제 공진은 거의 그대로인데 hero가 `crossover`, 폭 `0.725 µT`로
보고된다. `tests/test_magneto.py:170-176`은 두 축의 차이만 검사하고, readout 재중심은
검사하지 않는다. `tests/test_magneto.py:247-256`도 50 µT의 극단적 unresolved 상태만 다룬다.

개선은 새 solve 없이 가능하다. physical-zero 지표는 `raw["b_physical_ut"]`로 추출하고,
plot에는 command zero와 physical zero를 둘 다 표시한다. command-zero readout도 필요하면
라벨을 명시적으로 분리한다. 선형 보간/local extremum fit까지 포함해도 `O(n_B)` postprocess라
solve overhead는 무시할 수 있고, 기존 Hamiltonian 물리는 보존된다.

### 5.2 파라핀 2영역의 `rho_light` trace가 절대 감수율에 그대로 들어간다

현재 2영역 방정식은 light/dark density matrix를 합친 trace를 1로 둔다
(`gabes/schemes/magneto.py:562-575`). Trace-preserving Liouvillian에서 정상상태 light trace는

`Tr(rho_light) = gamma_in / (gamma_in + gamma_out)`

이다. 기본 `gamma_out/2pi=80 kHz`, `gamma_in/2pi=1 kHz`에서는 read-only capture로 모든
`(B,v)` 점에서 `Tr(rho_light)=1/81=0.012345679...`를 확인했다. 그 뒤
`_coherences()`는 이 unnormalised `rho_light`를 그대로 사용하고
(`gabes/schemes/magneto.py:493-496`, `604-618`), `chi_phys()`는 다시 local vapor density
`N_eff`를 곱한다 (`gabes/schemes/magneto.py:644-653`, `gabes/observables.py:393-407`).

이 convention은 명시되지 않았다. `N`이 셀의 국소 증기 밀도라면 beam 안의 감수율에는
conditional light-state `rho_light/Tr(rho_light)`가 들어가야 하므로, 현재 절대 흡수와 회전은
light/dark 체류확률을 한 번 더 곱해 희석할 가능성이 크다. 기본 EIT에서 비교하면 다음과 같다.

| convention | alpha(B=0) [m^-1] | T(B=0) | transmission feature |
|---|---:|---:|---:|
| 현재 unnormalised light block | 0.17504 | 0.998251 | 0.000200 |
| conditional `rho_light/Tr(rho_light)` 진단 | 14.1782 | 0.867811 | 0.013978 |

흡수계수는 정확히 81배 차이가 난다. 이는 폭·부호보다 **절대 contrast, rotation, slope,
감도 예측**에 직접 영향을 준다. 다만 원래 의도가 whole-cell occupancy density를 쓰는 것이라면
beam/cell volume fraction과 Maxwell coupling normalization을 명시해야 한다. 현재 코드에는 그
geometry가 없으므로 먼저 convention test와 문서화가 필요하다.

권고 구현은 (a) susceptibility용 conditional light state와 (b) 진단용 light occupation fraction을
분리하는 것이다. trace 계산과 나눗셈은 `O(n_B n_v n_level)`이고 새 matrix solve가 없어
overhead는 negligible하다. `gamma_in/out`이 달라도 local unpumped density가 변하지 않는지,
현재 부호·폭 trend는 유지되는지를 회귀 테스트로 고정해야 한다.

## 6. 기존 P0/P1 신뢰 경계 재평가

| 문제 | 현재 근거 | 비용 | 판단 |
|---|---|---|---|
| Ne broadening이 자연방출/TOC rate로 들어감 | `gamma_opt=Gamma+buffer_gamma`를 emission `gamma`로 전달한다 (`gabes/schemes/magneto.py:407-415`, `gabes/zeeman.py:96-118`) | 자연방출 `Gamma`와 optical pure dephasing을 같은 차원에서 분리하므로 거의 0 | **P0**. population lifetime/TOC와 elastic optical FWHM을 분리해야 한다. |
| buffer rate 두 knob가 사실상 합 하나 | 두 rate를 합쳐 `gamma_gg`와 `transit_rate`에 동일 전달한다 (`gabes/schemes/magneto.py:411-415`) | 같은 matrix 크기와 solve 수 | **P0**. transit reload와 추가 depolarization을 별도 dissipator 의미로 분리한다. |
| 다른 ground hyperfine manifold 누락 | 선택한 `Fg` 하나만 만들고 모든 `Gamma`를 그 manifold로 되돌린다 (`gabes/zeeman.py:89-118`, `tests/test_magneto.py:35-45`) | 경고는 0; compact reservoir는 대략 2배, full 두-ground Zeeman 2-region은 dense solve 기준 최대 약 6.8배 | 먼저 경고/branching 표, 이후 opt-in reservoir/reference fidelity가 적절하다. |
| 기본 Doppler quadrature 미수렴 | coarse grid와 scalar Voigt dilution (`gabes/schemes/magneto.py:68-100`, `460-496`) | 상태 표시는 0; adaptive reference는 비쌈 | **P0 status, P1 reference mode**. 이전 진단에서 9→641 class가 EIT 부호와 폭까지 바꿨다. |
| NMOR 중앙미분 미수렴 | `np.gradient(...)[ic]` (`gabes/schemes/magneto.py:671-680`) | local odd fit/5-point derivative는 거의 0 | **P1**. 이전 121/201/401점 값은 -2.5525/-2.9511/-3.1503 mrad/µT였다. |
| 문헌 검증 문구 과장 | 현재 실행은 1.195/1.739 mG, 문헌 0.12/0.20 mG인데 `MATCH`, `same sub-mG order` 출력 (`tests/verify_hanle_eit_eia.py:243-271`) | 출력 판정 수정만 필요 | **P0 trust**. 부호 PASS와 절대폭 FAIL/CHECK를 분리한다. |

저부하 물리 개선의 순서는 `rate/optical dissipator 의미 수정 -> 신뢰상태 표시 -> pressure shift,
저차 Dicke, pressure-to-relaxation mapping`이 적절하다. Full VCC는 이 범주에 포함시키면 안 된다.

## 7. 순수 코드 최적화 후보

물리와 출력 behavior를 바꾸지 않는 후보는 다음과 같다.

1. **kernel-side Doppler-weighted coherence contraction**: 현재 real-basis kernel은 모든
   `(B,v)`의 full `rho`를 complex 배열로 복원하고 (`gabes/kernels.py:399-410`), Python이 필요한
   세 coherence를 다시 축약한다 (`gabes/schemes/magneto.py:493-496`, `604-618`). 기본 8-level,
   `121 x 641` reference grid의 complex `rho`만 약 75.7 MiB다. kernel이 `chi_+`, `chi_-`,
   `chi_probe`의 weighted sum만 반환하면 결과를 유지하면서 allocation과 memory traffic을 줄인다.
2. **affine real-basis assembly를 한 번만 수행**: Python이 `C_xy+B*C_z` complex stack을 만든 뒤
   (`gabes/schemes/magneto.py:469-485`) wrapper가 다시 `U^dagger L U`를 batch transform한다
   (`gabes/kernels.py:399-405`). `C_xy`, `C_z`, dark coefficient를 한 번 real basis로 옮겨 kernel에서
   B-affine 조립하면 대형 complex stack과 반복 transform을 피할 수 있다. 1번보다 구현·동등성
   검증 비용은 크다.
3. **angular-momentum matrix 소형 cache**: `_hamiltonian()` 호출마다 같은 `Fx,Fy,Fz`를 만든다
   (`gabes/schemes/magneto.py:511-520`, `gabes/zeeman.py:48-66`). 안전하지만 heavy solve에 비해
   미세 최적화다.

`zeeman_manifold()`는 이미 LRU cache다 (`gabes/zeeman.py:69-70`). Headless observables와
`cell_mm`/`line_strength`의 navigate-only postprocess도 이미 구현돼 있어 중복 제안하지 않는다
(`gabes/schemes/magneto.py:233-236`, `286-290`, `620-746`).

## 8. 권고 우선순위와 최종 판단

1. **P0 absolute-scale semantics:** 파라핀 `rho_light`의 conditional/global convention을 정하고,
   local vapor density와 일관되게 정규화한다. 추가 solve 없이 가능하다.
2. **P0 readout correctness:** `b_offset_ut`가 있을 때 physical-zero feature와 command-zero
   readout을 분리한다. `O(n_B)` postprocess다.
3. **P0 dissipator semantics:** Ne optical pure dephasing과 자연방출/TOC를 분리하고,
   buffer transit reload와 collisional depolarization도 서로 다른 rate로 구현한다.
4. **P0 trust labels:** coarse Doppler, open-hyperfine/repump 미포함, 문헌 linewidth mismatch를
   hero/Derived/verification 출력에 명시한다. solve 비용은 없다.
5. **P1 reference accuracy:** adaptive Doppler reference와 수렴 검사를 추가하고, NMOR slope는
   local odd fit 또는 고차 중앙미분으로 추출한다. 전자는 고비용, 후자는 거의 무료다.
6. **P1 기존 저부하 물리:** pressure shift, 저차 Dicke correction, 의미가 분리된
   pressure-to-relaxation mapping을 문헌계수와 override를 보존해 구현한다.
7. **P2 performance:** kernel-side coherence contraction을 먼저, affine real-basis assembly를
   그 다음에 한다.

결론적으로 `MagnetoScheme`은 full `mF` manifold, CG, TOC, 실제 편광 drive, vector B,
two-region Ramsey exchange와 dispersive NMOR를 갖춘 **실질적인 원자물리 모델**이다. 현상 부호,
parameter sensitivity, measured trace의 초기 forward model로는 충분히 가치가 있다.

그러나 현재 absolute scale은 light-region normalization convention이 불명확하고, `b_offset` hero는
물리 공진 중심을 놓치며, buffer dissipator·open transition·Doppler/B-grid 수렴에도 중요한 경계가
남아 있다. 따라서 현 버전을 절대 linewidth·contrast·NMOR slope·자력계 감도의 독립 reference로
사용해서는 안 된다. 측정 자료로 B축·signal scale·relaxation/repump를 보정하고 수렴 상태를 확인하는
compact semi-quantitative model로 쓰는 것이 적절하다.

## 9. 검증 기록

- registry 직접 출력: `sas -> lambda -> rydberg_eit -> magneto -> fwm`
- `python tests/verify_hanle_eit_eia.py`: polarization 부호 전환 재현;
  linear/circular FWHM `1.195/1.739 mG` 대 문헌 `0.12/0.20 mG` 불일치 확인
- read-only 진단: `b_offset` command/physical-zero 지표, 파라핀 light-block trace 및 conditional
  normalization, 다섯 regime 대표 출력
- targeted:
  `python -m pytest tests/test_magneto.py tests/test_kernels.py tests/test_headless_observables.py tests/test_schemes_render.py tests/test_resonant_hanle_reference.py -q`
  -> **63 passed in 22.46 s**
- 전체: `python -m pytest -q` -> **232 passed in 64.15 s**
- 기존 dirty Rydberg/README/checklist/분석/보고서 작업은 보존했고 production code는 수정하지 않았다.
