# 2026-07-13 Scheme 4 물리 검토: Hanle / EIA / NMOR

## 선택 결과와 다섯 스킴의 현재 순서

현지 날짜의 일(day)은 13이므로

`n = (13 mod 5) + 1 = 4`

이다. 현재 드롭다운 등록 순서는 다음과 같다.

1. OD / SAS — `SASScheme()`
2. Lambda coherence (EIT / AT / CPT) — `LambdaScheme()`
3. Rydberg-EIT electrometry — `RydbergEITScheme()`
4. Hanle / EIA / NMOR — `MagnetoScheme()`
5. FWM — `FWMScheme()`

이 순서는 `_SCHEMES` 리스트가 결정한다 (`gabes/schemes/__init__.py:14-25`). 따라서 오늘의 검토 대상은 4번 `MagnetoScheme`이다 (`gabes/schemes/magneto.py:127-173`). README의 다섯 스킴 표도 같은 구성을 설명한다 (`README.md:10-16`).

## 조사 범위

- 스킴 설명·파라미터·프리셋: `gabes/schemes/magneto.py:127-388`
- 계산 본체와 관측량: `gabes/schemes/magneto.py:390-509`, `620-729`
- Zeeman manifold와 CG/TOC 방출: `gabes/zeeman.py:17-140`
- 실수기저 Numba 커널: `gabes/kernels.py:276-410`
- 단위·버퍼가스 상수: `gabes/constants.py:79-89`
- 물리 테스트: `tests/test_magneto.py:26-241`
- 커널 동등성·headless·렌더 테스트: `tests/test_kernels.py:62-99`, `tests/test_headless_observables.py:38-89`, `tests/test_schemes_render.py:31-42`
- 문헌 대조 실행 예: `tests/verify_hanle_eit_eia.py:1-277`
- 실제 데이터 보정 예: `analysis/resonant_hanle_squeezing_reference.py:352-412`, `527-625`
- 실험 설정 예: `analysis/resonant_hanle_experiment_config.example.json`
- 기존 개선안: `docs/checklist.json:22-33`, `60-64`, `102-106`
- 이전 Scheme 4 검토: `docs/daily_report/2026-06-23_scheme-4_hanle-eia-nmor.md`, `docs/daily_report/2026-06-28_scheme-4_hanle-eia-nmor.md`, `docs/daily_report/2026-07-03_scheme-4_hanle-eia-nmor.md`

별도 최상위 `examples/` 디렉터리는 없다. 대신 위의 검증 스크립트와 `analysis/` 설정·보정 경로가 실행 가능한 예제 역할을 한다.

## 구현된 물리의 평가

### 실제로 유용한 물리

이 스킴은 임의 Lorentzian을 그리는 toy model보다 훨씬 낫다.

1. `87Rb D1`의 선택된 `Fg -> Fe` 전이에 대해 모든 `mF` 준위를 만들고, `sigma+`, `sigma-`, `pi` 결합의 Clebsch-Gordan 계수를 계산한다 (`gabes/zeeman.py:89-118`).
2. 자발방출을 전이별 독립 jump가 아니라 동일 방출 편광별 `Sigma_q`로 묶어 excited coherence가 ground coherence로 전달되는 transfer of coherence(TOC)를 보존한다 (`gabes/zeeman.py:96-118`). 이 때문에 선형 편광·무잔류장 조건에서 `Fe = Fg + 1`은 intrinsic EIA, `Fe <= Fg`는 EIT가 되는 부호 규칙이 실제 계산에서 나온다 (`tests/test_magneto.py:179-192`).
3. QWP 각도를 원형기저 구동 진폭으로 바꾸고 (`gabes/schemes/magneto.py:113-124`), 종방향 스캔장과 횡방향 잔류장을 같은 Zeeman Hamiltonian에 넣는다 (`gabes/schemes/magneto.py:431-478`, `511-530`). 따라서 편광과 잔류장이 Hanle dip, MIA/EIA peak, circular-light LCA를 바꾸는 경로가 명시적이다 (`tests/test_magneto.py:123-138`, `195-211`).
4. 파라핀 셀은 illuminated/dark 두 영역의 교환 OBE로 wall-preserved coherence와 Ramsey narrowing을 표현하고 (`gabes/schemes/magneto.py:419-429`, `481-488`, `532-576`), 버퍼 셀은 단일 영역에서 ground relaxation과 collisional depolarization을 포함한다 (`gabes/schemes/magneto.py:407-418`, `489-491`, `578-601`).
5. NMOR는 transmission을 재사용한 가짜 readout이 아니라 `Re(chi_+ - chi_-)`에서 회전각을 만들고 영점 기울기를 magnetometer metric으로 보고한다 (`gabes/schemes/magneto.py:638-647`, `665-672`).
6. `87Rb` 밀도, D1 dipole strength, Beer-Lambert transmission, Doppler 평균을 사용한다 (`gabes/schemes/magneto.py:447-467`, `498-507`, `636-647`). 실험 데이터에는 별도 분석 경로가 B축 scale/offset과 신호 affine scale을 피팅한다 (`analysis/resonant_hanle_squeezing_reference.py:536-625`).

따라서 이 코드는 전이 선택, 편광, 잔류장, wall lifetime, ground relaxation을 바꿀 때 신호의 부호·상대 폭·민감도가 어떻게 변하는지를 조사하는 **실험 설계 및 sanity-check용 반정량 레퍼런스**로 유용하다.

### 정량 레퍼런스로서의 한계

절대 linewidth와 절대 magnetometer calibration에는 주의가 필요하다.

- 한 번에 하나의 `Fg -> Fe` manifold만 풀므로 다른 ground/excited hyperfine manifold를 통한 광펌핑, repump, 인접선의 동시 기여는 없다 (`gabes/schemes/magneto.py:390-416`).
- 파라핀의 `wall_coherence_ms`, `transit_relax_khz`, `dark_return_khz`는 coating·셀 형상·빔 크기에서 예측되는 값이 아니라 effective knob이다 (`gabes/schemes/magneto.py:241-273`, `419-429`).
- 버퍼 압력은 현재 optical homogeneous broadening만 직접 바꾸며 (`gabes/constants.py:79-89`, `gabes/schemes/magneto.py:407-415`), pressure shift, diffusion/Dicke narrowing, pressure-linked ground relaxation, velocity-changing collision은 없다. 저장소 자체도 이 모델을 compact·semi-quantitative로 명시한다 (`analysis/resonant_hanle_squeezing_reference.py:1147-1150`).
- Doppler 계산은 작은 등간격 velocity grid에 fine-grid scalar dilution을 곱하는 방식이다 (`gabes/schemes/magneto.py:68-100`, `460-467`). line-centre OD 크기 보정에는 싸고 유용하지만, 충돌로 결합된 velocity-class dynamics의 대체물은 아니다.
- `line_strength`는 절대 OD용 calibration factor이고 실험 조절량이 아니다 (`gabes/schemes/magneto.py:286-290`, `636-644`).

이번 실행에서 `python tests/verify_hanle_eit_eia.py`를 다시 돌렸다. 편광 전환과 전이별 부호는 재현됐지만 절대 폭은 다음과 같았다.

| 비교 | 계산 | 스크립트의 문헌 기준 | 판단 |
|---|---:|---:|---|
| 파라핀, linear CPT | 1.195 mG | 약 0.12 mG | 약 10배 넓음 |
| 파라핀, circular MIA | 1.739 mG | 약 0.20 mG | 약 8.7배 넓음 |
| 저출력 buffer circular LCA 예 | 4.722 mG | 약 2.4 mG | 약 2배 넓음 |

즉 `tests/verify_hanle_eit_eia.py:243-271`의 “same sub-mG order”라는 출력 문구는 현재 수치와 맞지 않는다. 같은 스크립트의 LCA 설명 (`tests/verify_hanle_eit_eia.py:174-206`)도 실행된 세 점 중 어느 것도 2.4 mG에 도달하지 않았다. 현재 자동 테스트는 부호와 narrowing 방향은 고정하지만 문헌 절대 폭을 assert하지 않는다 (`tests/test_magneto.py:123-167`, `179-223`). 이 때문에 스킴은 현상 분류에는 신뢰할 수 있어도 문헌 linewidth 재현 모델이라고 소개하면 과장이다.

## 기존 개선안과 계산 비용

기존 제안이 있으므로 새 물리를 무리하게 덧붙이기보다 우선순위를 재평가한다.

### 1. 버퍼가스 pressure shift + 저차 Dicke narrowing

`buffer-gas-pressure-shift`는 Ne broadening 상수를 gas/species/line table로 일반화하고 pressure shift와 phenomenological Dicke narrowing을 OD/SAS, Magneto, Lambda에 공통 적용하자는 GROUP B 항목이다 (`docs/checklist.json:22-26`).

- pressure shift는 `dL`에 압력 의존 scalar offset을 더하는 전처리라 solver 차원과 `(B, v)` 점 수를 바꾸지 않는다 (`gabes/schemes/magneto.py:452-479`). 계산 오버헤드는 사실상 무시할 수 있다.
- 저차 Dicke correction을 effective Doppler width 또는 homogeneous-width 보정으로 제한하면 마찬가지로 배열 크기와 solve 횟수는 그대로다. 대화형 속도를 거의 보존할 수 있다.
- 반면 full VCC는 velocity classes를 서로 결합해 현재 독립 `(B, v)` solve를 깨므로 GROUP C로 분리된 것이 타당하다 (`docs/checklist.json:102-106`). 이는 negligible-overhead 개선이 아니다.

### 2. Ne 압력에서 ground relaxation/depolarization 기본값 산출

`magneto-buffer-relaxation-map`은 `ne_pressure_torr`로부터 `buffer_ground_relax_khz`와 `collisional_depol_khz` 기본값을 정하고, 고급 override는 유지하자는 GROUP B 항목이다 (`docs/checklist.json:60-64`). 현재는 20 Torr preset조차 두 rate를 따로 지정한다 (`gabes/schemes/magneto.py:318-329`).

이 매핑은 solve 전에 scalar 두 개를 고르는 작업이라 계산 비용은 사실상 0이다. 경험계수와 온도·셀 형상 의존성을 문헌 또는 측정으로 정당화할 수 있다면, 물리를 훼손하지 않고 buffer-cell preset의 실험적 해석력을 크게 높일 수 있다. 단일 보편식으로 과신하지 않도록 source/cell metadata와 override를 남겨야 한다.

### 3. 이미 완료된 headless 관측량 경로

`figureless-observables-paths`는 완료 상태다 (`docs/checklist.json:29-33`). Magneto는 `observables(..., include_figures=False)`를 제공하고 (`gabes/schemes/magneto.py:620-729`), 테스트가 Matplotlib을 호출하지 않음을 확인한다 (`tests/test_headless_observables.py:46-67`). 아래 측정에서 headless readout은 0.3-0.8 ms로 heavy solve보다 두세 자릿수 작았다.

### 이번 근거로 추가할 수 있는 조심스러운 개선

가장 먼저 할 일은 새 고차 물리를 넣는 것이 아니라 **문헌 대조 스크립트의 주장과 검증 기준을 바로잡는 것**이다. `verify_hanle_eit_eia.py`가 계산/문헌 폭의 비를 출력하고, “정성 부호 PASS / 절대 폭 CHECK 또는 FAIL”을 분리하도록 하면 runtime 증가는 거의 없다. 이후 파라핀 effective relaxation과 LCA 조건을 실제 셀·빔 자료로 보정하되, 단지 목표 폭에 맞춘 hidden factor는 추가하지 않는 편이 좋다.

## 현재 성능과 순수 코드 최적화

현재 작업 트리에서 Numba가 활성화된 상태로 `scan_points=121`, `velocity_classes=9`, Doppler on을 사용해 각 조건을 한 번 예열한 뒤 3회 측정했다.

| 조건 | warm compute 3회 | headless observables | 대표 결과 |
|---|---:|---:|---|
| Paraffin EIT dip | 0.188 / 0.149 / 0.160 s | 0.790 ms | feature amp `-2.00e-2` |
| Paraffin EIA peak | 0.175 / 0.155 / 0.169 s | 0.435 ms | feature amp `+7.27e-2` |
| Buffer Hanle | 0.046 / 0.045 / 0.049 s | 0.315 ms | `T(0)=0.7565` |
| NMOR | 0.141 / 0.173 / 0.165 s | 0.525 ms | slope `-2.55 mrad/uT` |

동일한 대표 paraffin 조건은 실수기저 Numba kernel에서 0.143 s, NumPy fallback에서 0.505 s로 측정되어 현재 fast path가 약 3.5배 빨랐다. 커널은 Hermitian generator basis로 바꿔 복소 solve를 실수 solve로 만들고 B축을 병렬화한다 (`gabes/kernels.py:276-330`, `333-410`). B 의존 Liouvillian도 affine broadcast로 한 번에 조립한다 (`gabes/schemes/magneto.py:469-479`). 이들은 이미 반영된 주요 순수 최적화다.

남은 behavior-preserving 후보는 다음 순서가 합리적이다.

1. **커널에서 coherence contraction까지 융합**: 현재 커널은 모든 `(B,v)`의 full density matrix를 복원해 반환하고 (`gabes/kernels.py:329-330`, `409-410`), Python에서 `chi_+`, `chi_-`, probe coherence를 다시 순회·평균한다 (`gabes/schemes/magneto.py:493-496`, `604-618`). 필요한 세 coherence와 Doppler 합만 커널에서 직접 내보내면 full `rho` materialization과 메모리 트래픽을 줄일 수 있다. 내부 API만 바꾸므로 물리는 동일하며, 기존 fast/fallback parity test를 그대로 확장할 수 있다 (`tests/test_kernels.py:91-99`).
2. **실수기저 affine 계수를 직접 전달**: 현재는 복소 `L0_all`/`Ld_all` stack을 만든 뒤 커널 wrapper가 매 호출 `U^dagger L U`를 수행한다 (`gabes/schemes/magneto.py:477-485`, `gabes/kernels.py:323-329`, `401-409`). `C_xy`, `C_z`, B축을 실수기저로 한 번 변환해 커널에서 affine 조립하면 대형 복소 stack과 batched basis transform을 피할 수 있다. solver 결과는 바뀌지 않지만 구현·동등성 테스트 비용은 1번보다 크다.
3. **작은 immutable helper cache**: `_hamiltonian()`은 매번 ground/excited angular-momentum matrix를 다시 만든다 (`gabes/schemes/magneto.py:511-520`; `gabes/zeeman.py:48-66`). `angular_momentum_matrices(F)`를 read-only LRU cache로 만들 수 있다. `_doppler_dilution()`도 동일 `(T, gamma, detuning, grid)` 반복 sweep에서 재사용 가능하다 (`gabes/schemes/magneto.py:80-100`). 둘 다 안전하지만 전체 runtime 절감은 커널 융합보다 작을 가능성이 높다.

이전 보고서의 `zeeman_manifold` template cache 제안은 현재 이미 구현돼 있다 (`gabes/zeeman.py:69-70`). `core.hermitian_basis`와 identity도 캐시되어 있으므로 (`gabes/core.py:36-41`, `60-93`) 같은 제안을 다시 할 필요가 없다.

## 검증

- 문헌 대조 실행: `python tests/verify_hanle_eit_eia.py` — 완료, 37.1 s
- 관련 테스트: `python -m pytest tests/test_magneto.py tests/test_kernels.py tests/test_headless_observables.py tests/test_schemes_render.py -q`
- 결과: `33 passed in 11.39s`

- 전체 저장소: `python -m pytest -q` — `118 passed in 25.44s`

## 최종 판단

Scheme 4는 GABES의 다섯 스킴 중에서도 실제 원자물리 구조가 비교적 충실하다. Zeeman manifold, CG 결합, TOC 기반 intrinsic EIA, two-region Ramsey model, transverse-field LCA, `chi_+ - chi_-` NMOR readout은 실험가가 부호와 경향을 이해하는 데 실질적으로 쓸 수 있다.

그러나 현재 문헌 대조 실행은 파라핀 linewidth가 약 한 자릿수 크게 나오므로 **절대 linewidth 또는 절대 magnetometer 성능의 독립 예측 레퍼런스**로 사용하면 안 된다. 측정 trace에 B축·신호 scale과 relaxation parameters를 보정하는 전단 모델로 사용하는 것이 적절하다. 물리 가치 대비 비용이 가장 좋은 다음 단계는 기존의 pressure shift/저차 Dicke correction과 pressure-to-relaxation mapping이며, 우선 문헌 대조 스크립트의 과도한 MATCH 문구를 정량 기준으로 교정해야 한다.
