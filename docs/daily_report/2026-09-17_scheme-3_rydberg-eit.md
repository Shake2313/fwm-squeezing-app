# GABES 일일 Scheme 물리 검토 — 2026-09-17

## 0. 선택 결과와 검토 범위

- 로컬 날짜의 일(day) `17`에 대해 `n = (17 mod 5) + 1 = 3`이므로 오늘의
  대상은 **Scheme 3 — Rydberg-EIT electrometry**다.
- 현재 dropdown registry는 `gabes/schemes/__init__.py:19-25`에 다음 순서로
  정의된다.

  1. `SASScheme` — OD / SAS
  2. `LambdaScheme` — EIT / AT / CPT
  3. `RydbergEITScheme` — Rydberg-EIT electrometry
  4. `MagnetoScheme` — Hanle / EIA / NMOR
  5. `FWMScheme` — seeded FWM / biphoton

- `gabes/schemes/__init__.py:27-39`의 Lambda/Magneto alias는 API 직접 조회용이며
  위 다섯 dropdown 순서를 바꾸지 않는다.
- `CLAUDE.md`, `README.md`, `docs/checklist.json`, Scheme 3 선행 보고서 13개
  (2026-06-22…2026-09-12), 코드의 `TODO`/`FIXME`, 테스트, 분석 workflow,
  예제 config, 문헌 요약을 먼저 검색했다.
- GitHub REST API로 `Shake2313/fwm-squeezing-app`의 open/closed issue를 재확인했으며
  issue는 0건이다. 저장소 내 Rydberg 전용 TODO/issue 파일은 없고, 실제
  개선안은 `docs/checklist.json`과 일일 보고서에 축적돼 있다.

## 1. 직전 검토 이후 변화

현재 `HEAD` 는 `2998dc1`이다. 9월 12일 검토의 `c0c46f0`에서 현재
HEAD까지 Rydberg 원자계 코드, electrometry/experiment/CSV helper, 테스트,
분석 workflow에 관련된 commit diff는 없다. 현재 working-tree 파일의 blob도
각각 HEAD와 일치한다. `git status`에는 다른 사용자의 staged/unstaged 변경이
많으므로 이를 건드리지 않았다.

직전 보고서 이후 Rydberg에 직접 추가된 내용은 코드가 아니라 새
checklist 제안 `rydberg-doppler-on-pole-residue`다
(`docs/checklist.json:2716-2730`). 9월 14일 FWM/Magneto에서 검증한 pole–residue
엔진을 Rydberg Doppler-on 경로에 이식하자는 P3/small/ready 제안이다.
`analysis/pole_reuse/DEVLOG.md:160-173`은 이를 남은 이식 대상으로 명시한다.

공통 결과 export는 commit `8f97f05`에서 추가되었지만, checklist는 아직
`in_progress`다. full parameter/version/revision, JSON, CSV/ZIP, overlay provenance는
구현됐으나 named arrays, round-trip 재현 테스트, lazy extra-view export가
남아 있다(`streamlit_app.py:816-947`; `docs/checklist.json:3418-3439`).

## 2. 결론

**판정: 실제 Rydberg cascade EIT, microwave AT, 유한-IF 약신호 검출, 균형
검출기 noise budget을 연결한 유용한 반정량적 실험 설계 도구다. 그러나
절대 감도, 온도 optimum, velocity-resolved opacity, Zeeman/polarization, RF field map의
독립 기준으로 사용하기에는 아직 부족하다.**

적합한 용도는 다음과 같다.

- ⁸⁵Rb `5S–5P–40D` cascade EIT와 약 37 GHz `40D–39F` RF dressing의 line shape,
  linewidth, AT splitting, parameter trend 탐색.
- LO-dressed 정상상태 주위의 finite-IF weak-SIG response와 PSN/RIN/electronics
  제한의 조건부 비교.
- heater setpoint, effective vapor temperature, cold spot을 분리한 heated-cell
  시나리오, axial Beer–Lambert, effective atom number, SAM 불확도, CSV provenance
  workflow의 초기 설계.
- 각 입력의 provenance와 현장 calibration을 별도로 제공한 후의 조건부
  sensitivity budget 비교.

부적합한 용도는 다음과 같다.

- 기본 `6.08 nV/cm/√Hz`나 내부 온도 sweep optimum을 실장치의 절대 예측으로
  사용하는 일.
- Doppler-off empirical reference와 Doppler-on velocity-resolved mode를 동일한 opacity
  calibration으로 교체 사용하는 일.
- 낮은 contrast, 거의 0인 투과 power, scan-edge optimum, under-resolved AT에서도
  나오는 정밀한 숫자를 관측 가능한 feature로 해석하는 일.
- full Zeeman/polarization, cell standing wave, 비선형/time-domain lock-in, 가열에 따른
  technical drift를 현재 compact model이 대표한다고 가정하는 일.

## 3. 현재 구현이 담고 있는 실제 물리

### 3.1 원자계와 optical spectrum

- `RydbergEITScheme`은 ⁸⁵Rb
  `5S1/2 F=3 → 5P3/2 F'=4 → 40D5/2 → 39F7/2` 4준위 cascade를 명시한다
  (`gabes/schemes/rydberg.py:1-15`, `72-91`, `417-437`).
- probe와 481 nm coupling power는 beam diameter와 함께 `Ω ∝ √P/d`로 연결된다
  (`gabes/schemes/rydberg.py:439-461`). 실험 knob trend에는 유용하지만 절대 optical
  dipole에서 나온 값이 아니라 reference operating point에 맞춘 effective anchor다.
  코드도 bare D2 dipole이 optical pumping/Zeeman/hyperfine 누락으로 약 100배
  과대평가한다고 명시한다(`gabes/schemes/rydberg.py:41-55`).
- cold spot은 vapor pressure를, effective temperature는 local density, thermal motion,
  transit/dephasing을 정한다(`gabes/schemes/rydberg.py:463-505`). 이 분리는 heated-cell
  실험에 실용적이고, heater metadata만 바꿔서는 원자 solve가 달라지지
  않는다(`tests/test_rydberg_eit.py:456-475`).
- 801점 spectrum은 detuning-affine Liouvillian kernel을 사용하고, compiled kernel이
  없을 때도 같은 steady state를 풀도록 구성돼 있다
  (`gabes/schemes/rydberg.py:507-527`, `544-627`). 단순 선형/Voigt 그림이
  아니라 실제 포화 4-level OBE다.
- Doppler-on은 5P에 probe `k`, 40D/39F에 `(k_probe-k_coupling)` residual
  ratio를 넣어 Maxwell velocity class를 평균한다
  (`gabes/schemes/rydberg.py:72-91`, `610-619`; `gabes/atoms.py:95-109`). 기본
  Doppler-off는 residual 2-photon Doppler뿐 아니라 5P 1-photon Doppler도 함께
  끄는 calibrated stationary reference다. 따라서 두 mode는 단순한 fast/accurate tier가
  아니라 **다른 calibration regime**으로 봐야 한다.
- uncompensated residual Zeeman은 full magnetic-sublevel solve가 아니라 susceptibility의
  Gaussian convolution proxy다(`gabes/schemes/rydberg.py:530-542`, `621-627`).
  B-compensation trend에는 유용하지만 polarization/B-dependent absolute spectrum을
  대체하지는 못한다.
- Beer–Lambert 변환은 표준적이지만 `ls=0.001`이 숨은 effective opacity
  factor로 고정돼 있다(`gabes/schemes/rydberg.py:629-644`, `672-681`). 이는
  알려진 ⁸⁵Rb D2 cycling-line strength가 아니라 unresolved velocity fraction,
  ground/mF participation, optical pumping, mode/polarization, fit amplitude를 혼합한
  empirical participation anchor다.

### 3.2 finite-IF electrometry와 detector chain

- 약한 SIG는 LO-dressed Liouvillian 주위에서 `(L₀+iω_IF I) x = b`를 풀어
  sideband density matrix를 구한다(`gabes/rydberg_electrometry.py:116-194`). 이는
  정적 slope에 단위만 붙인 것이 아니라 실제 frequency-domain linear response다.
- Hermitian RF perturbation의 minus sideband만 풀고 plus sideband을 conjugate transpose로
  복원한다(`gabes/rydberg_electrometry.py:173-185`). 물리와 Hermiticity를
  보존하는 올바른 최적화이며 테스트가 solve 횟수까지 고정한다
  (`tests/test_rydberg_electrometry.py:39-94`).
- Doppler-on finite-IF는 complex phasor를 velocity weight로 먼저 평균한다
  (`gabes/schemes/rydberg.py:787-806`; `gabes/rydberg_electrometry.py:211-227`).
  magnitude를 먼저 평균하는 오류를 피한 점은 좋지만, 이 모드는 full
  detuning scan이 아니라 static maximum-slope 한 점만 해석한다.
- RF dipole, angular/polarization factor, RMS/peak convention이 분리돼 있고 AT split의
  detuning 보정과 field 변환도 명시적이다
  (`gabes/rydberg_electrometry.py:231-316`). ARC의 dipole API가 `e a₀` 단위의
  state-resolved matrix element를 반환한다는 정의와도 단위체계가 맞는다
  ([ARC documentation](https://arc-alkali-rydberg-calculator.readthedocs.io/en/latest/alkali_atom_functions.html)).
- balanced detector는 one-sided Schottky ASD, correlated RIN, electronics noise를 합성하고
  field responsivity로 나눈다(`gabes/rydberg_electrometry.py:403-483`). 문헌 감도
  숫자를 결과에 주입하지 않는 설계는 올바르다.
- 현재 reference arm power는 각 detuning에서의 signal transmission에 비례해
  점별로 다시 설정된다(`gabes/schemes/rydberg.py:902-927`). 실제 fixed arm
  또는 한 operating detuning에서만 맞춘 balance가 아니라 ideal pointwise-balanced
  envelope이므로 sensitivity scan을 낙관적으로 만들 수 있다.
- low transmission/contrast, edge optimum, samples-per-feature, successive local refinement을
  통합한 validity gate가 아직 없다. contrast가 양수이기만 하면 linewidth를
  내고, unresolved weak-LO 조건에서도 `Total sensitivity`를 hero로 유지한다
  (`gabes/schemes/rydberg.py:979-1003`, `1078-1094`, `1127-1207`;
  `tests/test_rydberg_eit.py:134-141`).

### 3.3 실험 분석 helper와 예제

- `gabes/rydberg_experiment.py:53-95`, `171-338`, `368-448`은 온도 상태,
  axial profile/segmented Beer–Lambert, Gaussian-overlap atom number를 분리한다.
- `gabes/rydberg_experiment.py:452-594`는 SAM field, first-order standard uncertainty,
  `2D²/λ` far-field 경고를 제공하고 3-D horn/cell solver라고 과장하지 않는다.
- `gabes/rydberg_experimental_csv.py:35-100`, `441-691`은 EIT/RF/PSD trace의
  단위, 중복, non-finite row, SHA-256 provenance와 unit assumption을 보존한다.
- workflow는 `MEASURED/PPT/REFERENCE/FITTED/PREDICTED/ASSUMED/PENDING`을
  분리한다(`analysis/rydberg_cell_heating/README.md:22-39`). 이는 실험 참조
  코드로서 매우 좋은 설계다.
- 다만 예제 config의 한 `base_params`에 문헌 조건과 미측정 detector/path/noise
  가정이 혼재한다(`analysis/rydberg_cell_heating/example_config.json:8-27`). raw
  EIT, sensitivity, temperature log는 모두 `PENDING`이고(`:53-78`), RF dipole만
  `REFERENCE`, 전체 electrometry는 `ASSUMED`다(`:88-102`). 따라서 이는 pipeline
  demo이지 held-out validation example이 아니다.

## 4. 문헌 기준과 오늘의 수치 재검증

Ju *et al.*은 1.6 MHz EIT linewidth와 37 GHz에서
`12.5(8) nV/cm/√Hz` 실측 감도를 보고한다
([arXiv:2606.04354](https://arxiv.org/abs/2606.04354)). 저장소의 apparatus 요약은
paper PSN을 `11.2(4) nV/cm/√Hz`로 기록하고, paper 온도 모델 optimum을 약
33 °C로 요약한다(`references/rydberg_ju_experiment_context.md:20-22`,
`65-71`, `94-100`, `140-155`).

현재 worktree에서 직접 계산한 결과는 다음과 같다.

| 점검 | 결과 | 해석 |
|---|---:|---|
| 기본 EIT `T(0)` | 0.940016615 | 기본점에서 관측 가능 |
| 기본 EIT contrast | 0.135391492 | 정상적인 calibrated-reference 형태 |
| 기본 EIT FWHM | 1.613058583 MHz | 문헌 1.6 MHz scale 재현 |
| 기본 AT `T(0)` | 0.806029283 | detector power가 남는 조건 |
| 기본 AT split | 3.496500 MHz | LO Rabi 3.7 MHz와 같은 scale |
| finite-IF optimum | -2.220000 MHz | 기본 scan 내부 |
| 기본 PSN = total | 6.082387215 nV/cm/√Hz | paper PSN 11.2보다 1.84배 낙관적 |
| 기본 heating sweep optimum | 47 °C, 0.971731 nV/cm/√Hz | paper 요약 33 °C와 불일치; 절대 optimum으로 사용 불가 |

기본 spectrum이 빠르고 문헌 scale을 재현하는 점은 유용하지만, 그 일치는
Rabi, transit factor, opacity가 reference point에 맞춰졌기 때문에 독립 예측 검증이
아니다. 특히 heating sweep은 순수 model 내부에서 47 °C와 paper 요약의
33 °C가 다르므로, 현장 운전점 권고가 아니라 `PREDICTED` scenario로만
다루어야 한다.

### 4.1 숨은 opacity 민감도

AT raw의 `ls`만 진단적으로 바꿔 원자 solve를 다시 하지 않은 결과다.

| `ls` | total sensitivity [nV/cm/√Hz] |
|---:|---:|
| 0.0005 | 11.781990083 |
| 0.0010 | 6.082387215 |
| 0.0020 | 3.238775723 |

단일 hidden scalar의 4배 범위에서 감도가 3.64배 변한다. 이것이
`rydberg-effective-opacity-and-sweep-trust`가 full Zeeman 보다 먼저인 이유다.

### 4.2 관측 불가 조건의 precise-looking result

80 °C EIT에서 `T(0)=contrast=1.798923×10⁻⁶`인데도 FWHM은
`0.307015 MHz`로 보고되고 static IF optimum은 -9 MHz scan 경계 바로 안쪽인
`-8.96 MHz`를 고른다. 즉 숫자 계산은 정상적이어도 실험적 특징 검출은
실패한 조건이다. contrast/power/edge gate는 선택적 표시 기능이 아니라
오해 방지에 필요한 물리 계약이다.

### 4.3 Doppler-off와 Doppler-on의 큰 간격

| 기본 조건 | Doppler off | Doppler on |
|---|---:|---:|
| EIT `T(0)` | 0.940016615 | 0.996638726 |
| EIT contrast | 0.135391492 | 0.001885647 |
| EIT FWHM | 1.613059 MHz | 2.921208 MHz |
| AT split | 3.496500 MHz | 2.386500 MHz |
| AT total sensitivity | 6.082387 | 1007.280836 nV/cm/√Hz |
| finite-IF optimization | full 801-point scan | static-slope 1-point coherent average |

이 격차는 Doppler-on을 빠르게 만드는 것이 곧 물리 정확도 개선은 아님을
보여 준다. 먼저 stationary-reference와 velocity-resolved mode의 opacity/participation을
별도로 검증해야 한다. pole–residue는 이 trust boundary를 해결하지 않는다.

## 5. 기존 개선 제안의 계산 비용과 물리 보존성

기존 제안이 충분하므로 새로운 대형 physics proposal은 필요하지 않다.

| 기존 항목 | 비용·물리 보존성 | 판단 |
|---|---|---|
| extra-view dependency correctness (`docs/checklist.json:2295-2312`, P1 ready) | cache key/계층 분리는 O(1) metadata와 작은 memory 비용. IF만 phasor, cell/opacity만 Beer–Lambert, detector/SAM만 후처리를 무효화하면 원자 물리는 그대로다. | 현재 `ExtraView.param_keys`와 phasor LRU는 일부만 반영. **P1** |
| effective opacity / sweep trust (`docs/checklist.json:2315-2335`, P1 ready) | `ls=0.001`을 named input으로 옮기는 기본 비용 0. default fixture로 현재 spectrum을 `≤1e-12` 보존 가능. uncertainty sweep만 선택적 선형 비용. | 절대값 신뢰도의 가장 싼 개선. **P1** |
| readout validity guardrails (`docs/checklist.json:2338-2360`, P1 needs_decision) | contrast, power, samples/FWHM, edge distance는 이미 있는 801점 배열의 O(N) 후처리. fixed reference arm은 원자 physics를 바꾸지 않고 receiver model을 실제 조건에 맞춤. | threshold 정책만 결정하면 됨. **최우선 P1** |
| held-out validation (`docs/checklist.json:2363-2385`, P2 blocked_external) | runtime보다 측정, metadata, uncertainty, hash/provenance 비용이 지배. 내부 계산으로 대체 불가. | 데이터 전에 absolute/optimum은 `PREDICTED`. |
| Doppler-on pole–residue (`docs/checklist.json:2716-2730`, P3 ready) | velocity class당 solve 대신 scan row당 작은 eigendecomposition 1회. `PoleResidues.discrete_mean`으로 현재 2 m/s·4σ 절단 Maxwell 측도를 그대로 합산하면 물리/수치 계약 보존. untruncated `gaussian_mean`은 비용은 작지만 결과를 바꾸는 별도 물리 결정. | 오늘 EIT compute off 3.91 ms, on 147.98 ms. default가 off이고 heating sweep도 현재 303 ms이므로 절대 이득이 작아 **P3가 맞음**. |
| Zeeman/polarization reference (`docs/checklist.json:2899-2923`, P2 parked) | Hilbert/Liouville 차원과 velocity/B grid가 함께 커짐. 특정 state/polarization/B dataset 없이는 비식별 자유도만 증가. | opt-in research mode로 parked가 타당. |
| measured transfer function (`docs/checklist.json:2926-2950`, P2 blocked_external) | 측정 complex transfer를 현재 phasor에 곱하는 것은 O(Nf), 실질적으로 negligible. linear response 물리를 보존. | 측정 transfer 확보 후 full time-domain보다 먼저. |
| RF field-map import (`docs/checklist.json:2953-2977`, P2 blocked_external) | O(Nspace×Nscan), mesh convergence 비용. 내부 3-D EM 구현보다 provenance-tagged 외부 map import가 타당. | named geometry/map이 있을 때 opt-in. |

### pole–residue 제안의 구체적 gate

`gabes/pole_doppler.py:1-31`, `82-123`은 affine detuning response를 정확한 유리함수로
표현하고 discrete 측도와 continuous Gaussian 평균을 별도 API로 제공한다.
Rydberg 이식은 다음 조건을 만족할 때만 behavior-preserving 최적화다.

1. 현재 `velocity_grid(T, dv=2, cutoff_sigma=4)` 노드와 weight를
   `discrete_mean` 순서까지 그대로 사용한다.
2. 각 row에서 하나 이상의 direct solve guard를 두고 max relative error `≤1e-6`를
   통과한다.
3. Doppler-off 경로는 bitwise 불변으로 둔다.
4. finite-IF Doppler-on의 coherent phasor average는 pole row를 재사용하거나, 단일
   static-slope path를 유지하는 이유를 명시한다.
5. untruncated Faddeeva 평균을 채택하려면 속도 patch와 분리해 opacity/
   participation을 재검증한다.

## 6. 동작을 바꾸지 않는 순수 코드 최적화

### 6.1 finite-IF Doppler-off Liouvillian affine batch 조립

`_superheterodyne_atomic_phasors`는 801개 detuning에서 Python list comprehension으로
`core.build_liouvillian`을 호출한다(`gabes/schemes/rydberg.py:808-815`). 이
Hamiltonian/Liouvillian은 detuning에 정확히 affine이므로 기준 행렬과 1 rad/s
coefficient를 한 번만 만들고 broadcasting하면 된다.

- 현재 list 조립 중앙값: `36.281 ms`
- affine batch 중앙값: `3.093 ms`
- 조립 단계 speedup: `11.73×`
- 검증: `np.array_equal=True`, 최대 절대차 `0`

grid 축소나 근사가 아니라 **똑같은 행렬을 더 싸게 조립**한다. exact-array
fixture를 추가하면 가장 안전한 첫 성능 patch다.

### 6.2 detector scan vectorization

atomic phasor cache hit 후에도 801점 각각에서 `PhotodiodeChannel`,
`BalancedDetector`, noise, sensitivity dataclass를 만든다
(`gabes/schemes/rydberg.py:908-946`). 오늘 AT readout은 cold `73.48 ms`, phasor
cache hit 후에도 `20.32 ms`였다. signal/reference current, Schottky variance,
correlated RIN, electronics, field ASD를 NumPy 배열식으로 계산하고 최적점 하나에만
dataclass를 만들면 식·단위·결과를 그대로 보존할 수 있다. scalar helper
결과와 전체 배열의 `allclose` fixture를 함께 두어야 한다.

### 6.3 heating extra view의 네 계층 cache

temperature sweep은 온도 10개에서 EIT와 AT `compute/_readout`을 반복한다
(`gabes/schemes/rydberg.py:1384-1465`). 오늘 중앙 runtime은 `303.14 ms`였다.
다음 경계를 별도 bounded cache로 나누면 물리를 바꾸지 않고 detector-only
조정에서 원자 solve를 없앨 수 있다.

1. temperature/probe/coupling/dephasing → atomic raw spectrum
2. IF/LO/RF dephasing → finite-IF atomic phasor
3. cell/opacity → Beer–Lambert transmission
4. QE/path/reference/noise/RF dipole/SAM → detector/field conversion

### 6.4 formatted metric 재파싱 제거

temperature sweep은 `_readout`의 표시용 문자열에서 slope/discriminator를
`split()[0]`으로 다시 숫자로 읽는다(`gabes/schemes/rydberg.py:1412-1423`). 큰
병목은 아니지만 format/locale 변경에 취약하다. numeric field를 readout 구조에
반환하고 renderer만 format하면 동일 물리를 더 안전하게 재사용한다.

### 6.5 pole–residue의 위치

pole–residue 이식도 현재 discrete measure를 유지하면 pure optimization에 포함된다.
다만 default path는 이미 4 ms 내외이고, 전체 AT readout에서는 affine 조립과
detector loop의 절대 시간이 더 직접적이다. 따라서 순서는 **affine batch →
detector vectorization/cache layering → Doppler-on pole 이식**이 타당하다.

## 7. 문서·테스트·예제의 참조 적합성

### 잘 된 부분

- README는 finite-IF, detector budget, cold spot, SAM, conditional sensitivity, deferred
  scope를 현재 코드와 맞게 설명한다(`README.md:115-137`).
- Scheme/helper/workflow 테스트는 defaults, linewidth, AT split, finite-IF equation,
  Hermiticity, RMS/peak, detector noise, cache, temperature/cold spot, axial/SAM, CSV
  provenance를 넓게 다룬다. 특히 low-IF response를 독립 static finite difference와
  비교하는 `tests/test_rydberg_eit.py:309-337`은 단순 golden value보다 강한
  물리 테스트다.
- 분석 workflow는 clean CSV, SHA-256 source manifest, evidence status, generated JSON/
  figure/TeX macro를 제공하고 spatial approximation의 한계를 명시한다
  (`analysis/rydberg_cell_heating/README.md:14-39`, `72-117`).
- 공통 app export는 같은 render에 사용한 parameter/raw/plot/provenance를 외부로
  가져갈 수 있게 해 재현성을 높였다.

### 수정이 필요한 부분

1. HTML User Guide는 Rydberg를 여전히 정적 EIT/AT와 dispersion 위주로만
   설명한다(`docs/Userguide/GABES_User_Guide_v2.html:535-538`, `617-625`).
   finite-IF, detector budget, temperature/cold spot, SAM, conditional status와 현재
   3-panel AT figure가 빠져 있다.
2. `tests/test_rydberg_eit.py:1-5`의 module docstring은 “public UI shows the static
   spectrum only”라고 해 현재 sensitivity UI와 모순된다.
3. weak-LO 테스트는 unresolved AT에서도 `Total sensitivity` hero가 남아야 한다고
   고정한다(`tests/test_rydberg_eit.py:134-141`). guardrail 구현 시 numeric hero
   suppression/status 계약으로 바꿔야 한다.
4. high-temperature tiny-contrast, edge optimum, fixed reference arm, successive refinement,
   Doppler-mode별 opacity provenance를 acceptance test가 직접 고정하지 않는다.
5. 공통 export는 round-trip, complex/named arrays, lazy extra view까지 테스트해야
   experiment record로서 완결된다.
6. example config의 문헌 reference 값과 미측정 assumption에 parameter-group별 status를
   붙여야 “paper apparatus reference”와 “current detector assumption”이 분리된다.

## 8. 검증 결과

- Scheme 3 집중 6개 테스트 파일: **69 passed in 22.22 s**
  - `tests/test_rydberg_eit.py`
  - `tests/test_rydberg_electrometry.py`
  - `tests/test_rydberg_experiment.py`
  - `tests/test_rydberg_experimental_csv.py`
  - `analysis/rydberg_cell_heating/test_workflow.py`
  - `analysis/rydberg_cell_heating/test_adapter_integration_regressions.py`
- AGENTS.md가 요구한 root `python -m pytest -q`는 **collection 실패**했다.
  23.74 s 후 `references/AutoOD/test_autood_gui_adjustments.py`가 선택 의존성
  `sympy` 없이 `rb_od_calculator_optimized_final.py`를 import해
  `ModuleNotFoundError: No module named 'sympy'`가 났다. Scheme 3 테스트 실패가
  아니라 수집 환경 문제다.
- 기존 staged/unstaged 사용자 변경은 보존했고 production source, checklist,
  example은 수정하지 않았다. 이 보고서만 추가했다.

## 9. 권장 실행 순서

1. **P1 — observability/fixed-reference gate:** low contrast/power, edge, unresolved,
   unconverged 결과의 numeric hero를 억제하고 reference-arm balance detuning을
   하나로 고정한다.
2. **P1 — effective opacity 공개:** `ls=0.001`을 provenance, 범위, uncertainty,
   Doppler-mode calibration scope가 있는 `effective_opacity_participation`으로 분리한다.
3. **P1 — exact 성능 patch:** affine Liouvillian batch, detector vectorization, 4-layer
   cache를 이 순서로 구현하고 exact-array/scalar-parity/cache-invalidation 테스트를
   둔다.
4. **P1 — 문서·테스트 동기화:** User Guide, stale test docstring, weak-LO hero
   contract, export round-trip을 현재 기능과 guardrail 정책에 맞춘다.
5. **P2 — held-out data:** measured EIT/RF/PSD/temperature, fixed reference arm, complex
   readout transfer가 있을 때 절대 감도와 온도 optimum을 재검증한다.
6. **P3 — Doppler-on pole–residue:** Doppler mode가 실제 기본/sweep 병목이 된 뒤,
   현재 discrete measure와 row guard를 보존하며 이식한다.

현재 가장 중요한 일은 full Zeeman이나 time-domain 엔진을 늘리는 것이
아니라, **이미 있는 compact physics의 관측 가능성, opacity provenance, fixed receiver
조건을 명시적으로 gate하는 것**이다. 이 세 가지는 추가 solve가 거의 없으며,
실험 참조 코드로서의 신뢰도를 가장 크게 높인다.
