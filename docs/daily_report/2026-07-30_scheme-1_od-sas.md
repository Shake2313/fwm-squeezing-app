# 2026-07-30 Scheme 1 물리 리뷰 — OD / SAS

## 선택 규칙과 현재 다섯 scheme

- 로컬 날짜는 `2026-07-30`, day-of-month는 `30`이다.
- `n = (day mod 5) + 1 = (30 mod 5) + 1 = 1`이므로 오늘은 **Scheme 1**을 검토한다.
- 실제 UI 순서는 registry의 `_SCHEMES`가 정한다. 현재 코드와 README의 순서는 서로 일치한다
  (`gabes/schemes/__init__.py:12-24`, `README.md:8-16`).

| 번호 | 현재 scheme | registry 구현 | 주된 물리/출력 |
|---:|---|---|---|
| 1 | OD / SAS | `SASScheme()` | pump-off Doppler OD, pump-on Lamb dip/crossover |
| 2 | Lambda coherence (EIT / AT / CPT) | `LambdaScheme()` | 3준위 coherence, EIT/AT/CPT |
| 3 | Rydberg-EIT electrometry | `RydbergEITScheme()` | cascade EIT, microwave AT, electrometry |
| 4 | Hanle / EIA / NMOR | `MagnetoScheme()` | 영자기장 투과와 자기광학 회전 |
| 5 | FWM | `FWMScheme()` | seeded gain/squeezing, SFWM biphoton |

오늘의 대상은 내부 이름 `sas`, UI 제목 `Absorption spectroscopy (OD / SAS)`인
`SASScheme`이다 (`gabes/schemes/sas.py:53-63`).

## 조사 범위와 기존 제안 검색

현재 구현을 판단하기 전에 다음 자료를 먼저 검색했다.

- 저장소 지침과 개요: `CLAUDE.md`, `README.md`
- 계획/TODO: `docs/checklist.json`, `gabes/constants.py`
- 이전 Scheme 1 보고서:
  `docs/daily_report/2026-06-25_scheme-1_od-sas.md`,
  `docs/daily_report/2026-07-10_scheme-1_od-sas.md`,
  `docs/daily_report/2026-07-20_scheme-1_od-sas.md`
- 구현: `gabes/schemes/sas.py`, `gabes/species.py`, `gabes/core.py`,
  `gabes/doppler.py`, `gabes/lineshape.py`, `gabes/experimental_csv.py`,
  `streamlit_app.py`
- 테스트: `tests/test_sas.py`, `tests/test_experimental_csv.py`,
  `tests/test_absorption.py`, `tests/test_headless_observables.py`,
  `tests/test_schemes_render.py`
- 예제/사용 설명: `docs/Userguide/GABES_User_Guide_v2.html:518-573`,
  `README.md:70-92`. 별도 `examples/` 디렉터리는 없으며, 사용자 가이드의
  `od.png`/`sas.png`가 현재의 계산 예제다.
- 로컬에서 별도 issue/proposal/note 파일은 발견되지 않았다. 공개 GitHub
  [Issues](https://github.com/Shake2313/fwm-squeezing-app/issues?q=is%3Aissue)는
  2026-07-30 확인 시 open/closed 모두 합쳐 0건이었다.

기존 개선 제안은 **있다**. 체크리스트에는 buffer-gas pressure
shift/Dicke narrowing, low-order polarization/Zeeman proxy, full
velocity-changing collision(VCC)이 남아 있다
(`docs/checklist.json:50-54`, `docs/checklist.json:116-120`,
`docs/checklist.json:130-134`). 2026-07-20 보고서는 여기에 더해
species-mode pump OBE의 collisional optical dephasing 누락, sub-Doppler
readout의 해상도 의존성, 전 species/line에 고정된 `I_sat`, pump-off fast
path와 running-median vector화를 구체적으로 제안했다
(`docs/daily_report/2026-07-20_scheme-1_od-sas.md:74-144`).

이전 Scheme 1 검토 뒤에는 두 가지 사용자-facing 변경이 들어왔다.

- `1b3b276`: 실험 CSV import/correction/overlay 추가
- `bd44273`: Transmission/OD 플롯 탐색과 metric hierarchy 변경

따라서 이번 보고서는 기존 제안을 재평가하면서, 새 CSV와 lock readout이 실험
reference로서 갖는 신뢰 경계를 새로 검토한다.

## 결론 요약

**OD/SAS는 실제 원자물리를 구현하며, pure warm-vapor cell의 선 배정,
절대 OD scale sanity check, pump-power/transit-rate 경향, SAS lock 후보 탐색에
유용한 semi-quantitative 실험 레퍼런스다.** 단순 Lorentzian 합이 아니라 같은
hyperfine manifold에서 pump-off OD와 pump-on SAS가 연속적으로 이어지고,
CG-branched spontaneous decay와 hyperfine optical pumping이 crossover를
만든다 (`gabes/schemes/sas.py:4-30`, `gabes/species.py:337-422`).

다만 현재 상태를 다음 용도의 절대 reference로 사용하면 안 된다.

- buffer-gas cell의 pressure-shifted line centre와 Dicke/VCC linewidth
- 편광/Zeeman-resolved contrast와 line strength
- 변조·복조·검출기 전달함수를 포함한 실제 servo error-signal slope
- CSV overlay로부터 얻은 절대 transmission/OD
- feature 검출에 실패한 조건에서의 hero `Lock Slope`

오늘 새로 확인한 가장 중요한 신뢰 문제는 두 가지다.

1. sub-Doppler feature가 검출되지 않아도 full-spectrum 기울기를 `Lock Slope`
   hero로 표시한다. 진단 조건에서는 이 값의 위치가 가장 가까운 hyperfine
   marker에서 **2.23 GHz** 떨어진 스캔 외곽이었다.
2. CSV의 자동 floor/ceiling 보정은 물리적인 `T=0.70–0.92` trace도 정확히
   `0–1`로 다시 매핑한다. 따라서 현재 overlay는 **상대 line-shape 비교**이지
   absolute transmission/OD calibration이 아니다.

둘 다 추가 OBE solve 없이 고칠 수 있다. feature 미검출 상태를 hero로 되돌리고,
CSV를 “relative normalized transmission”으로 명시하거나 dark/reference level을
받는 calibration mode를 더하면 된다.

## 현재 구현이 담는 실제 물리

### 1. 같은 모델 안의 pump-off OD와 pump-on SAS

펌프와 프로브는 반대 방향이므로 원자 좌표계 detuning이 각각
`Δ + kv`, `Δ - kv`가 된다. 펌프를 끄면 population factor가 1인 선형 Voigt
흡수로 축약되고, 펌프를 켜면 velocity-selective saturation과 optical pumping이
같은 spectrum에 Lamb dip/crossover를 만든다
(`gabes/schemes/sas.py:4-30`, `gabes/schemes/sas.py:215-246`).

이 구조는 OD와 SAS를 서로 무관한 toy curve 두 개로 붙인 것보다 실험적으로 훨씬
유용하다. 온도, cell length, isotope, D1/D2, pump power, waist, transit
relaxation을 한 모델에서 바꾸며 경향을 볼 수 있다
(`gabes/schemes/sas.py:65-106`).

### 2. 실제 alkali hyperfine data와 optical pumping

Rb-85, Rb-87, Cs-133의 hyperfine A/B, centroid frequency, mass, natural
linewidth, abundance가 명시되어 있다 (`gabes/species.py:130-182`). 각
`Fg→Fe` 전이의 상대 세기는 Wigner-6j/3j에서 계산되며, ground degeneracy를
포함한 동일한 `T(Fg,Fe)`가 line weight, decay branching, pump Rabi의 상대
크기에 일관되게 들어간다 (`gabes/species.py:366-414`).

`build_manifold()`는 모든 허용 hyperfine 전이와 CG-branched decay를 만들고,
모든 상태를 thermal ground distribution으로 되돌리는 transit relaxation으로
무한 optical pumping을 regularize한다 (`gabes/species.py:337-399`). 이 decay
redistribution이 실제 alkali SAS에서 중요한 enhanced/inverted crossover를
만든다.

### 3. 절대 OD 기준과 테스트

약한 probe의 line-integrated absorption은 AutoOD convention의 dipole과
`C_F²` normalization을 사용한다 (`gabes/species.py:226-270`). 테스트는 다음을
보호한다.

- 85Rb D1 `C_F²`, reduced dipole, CRC density
  (`tests/test_sas.py:58-68`)
- known hyperfine splitting과 manifold line count
  (`tests/test_sas.py:71-84`)
- decay branching과 합계
  (`tests/test_sas.py:87-98`)
- pump-off AutoOD absolute scale와 49/25 manifold ratio
  (`tests/test_sas.py:101-127`)
- pump-on sharp feature와 hyperfine-pumping crossover
  (`tests/test_sas.py:130-154`)
- natural-Rb isotope overlay와 species/line별 preset
  (`tests/test_sas.py:157-171`)

사용자 가이드도 같은 natural-Rb D1 조건의 pump-off OD와 pump-on SAS 그림을
나란히 제공한다 (`docs/Userguide/GABES_User_Guide_v2.html:549-572`).

오늘 기본값의 현재 출력은 다음과 같았다.

| 항목 | 현재 값 |
|---|---:|
| Peak OD | `0.54` |
| Narrowest sub-Doppler | `32.2 MHz` |
| Lock Slope | `0.0241 /MHz` |
| Lock Detuning | `-1304.9 MHz` |
| 첫 headless compute | `1.205 s` |

이 수치 중 Peak OD와 큰 선 구조는 유용하지만, 아래 해상도 진단 때문에
sub-Doppler FWHM과 lock slope의 표시 자릿수를 정량적으로 신뢰해서는 안 된다.

## 새 CSV 비교 경로의 가치와 한계

### 구현상 강점

CSV importer는 첫 두 열만 읽고 metadata/invalid row를 진단하며, 같은 x를
median으로 합치고 큰 gap을 분리한다
(`gabes/experimental_csv.py:127-190`, `gabes/experimental_csv.py:260-344`).
자동 보정은 4.5-MAD Hampel filter와 실제 x 좌표를 쓰는 local quadratic
regression을 적용한다 (`gabes/experimental_csv.py:347-486`). 좁은 line을
impulse로 오인하지 않도록 multi-sample feature와 span-defining singleton을
보존하는 방어도 있다 (`gabes/experimental_csv.py:389-456`).

테스트는 CP949/UTF-16, invalid row, duplicate/gap, noisy Gaussian feature,
1–5 sample Lamb dip, sparse extrema, 실제 `ReferenceOD.csv`를 다룬다
(`tests/test_experimental_csv.py:25-263`). UI는 원본 trace 선택 표시, 보정
floor/ceiling/window, warning, x scale/shift/reverse, detector polarity inversion을
노출한다 (`streamlit_app.py:890-1072`).

이 경로는 “scope trace를 빠르게 깨끗하게 만들어 계산 곡선과 모양을 맞춰 보는”
실험실 workflow에는 유용하다.

### 새 발견 1 — 자동 0–1 mapping은 절대 transmission을 보존하지 않는다

`_normalize_transmission()`은 trace의 extrema 부근에서 floor/ceiling을 잡은 뒤
`(signal-floor)/(ceiling-floor)`를 `[0,1]`로 clip한다
(`gabes/experimental_csv.py:558-641`). 실제 dark level 또는 별도 reference
beam을 받는 calibration은 아니다. 현재 warning은 low contrast와 sparse
extrema만 검사한다 (`gabes/experimental_csv.py:605-639`).

진단용으로 noise 없는 실제 transmission

```text
baseline = 0.92
minimum  = 0.70
```

을 CSV로 넣었을 때 importer는 이를

```text
baseline = 1.00
minimum  = 0.00
```

으로 바꾸었다. 알고리즘의 명시된 동작과 일치하지만, 이 trace를 이론의 absolute
transmission 또는 `OD=-ln(T)`와 비교하면 contrast를 과대평가한다.

**권고:** 즉시 UI/README label을 `Relative normalized transmission`으로 명시하고,
후속으로 `(dark, reference)` 또는 `(gain, offset/I0)`를 받는 calibration mode를
추가한다. 이는 O(N) affine post-processing뿐이므로 solve overhead가 없고,
물리적 absolute contrast를 보존할 수 있다. 자동 extrema mode는 빠른
shape-overlay용으로 유지하면 된다.

### 새 발견 2 — non-monotonic/bidirectional sweep가 한 곡선으로 붕괴한다

`_sort_and_merge()`는 원래 취득 순서를 버리고 x로 정렬한 뒤 정확히 같은 x의 y를
median으로 합친다 (`gabes/experimental_csv.py:310-324`). 이는 단일 monotonic
sweep에는 적절하지만, 왕복 scan·여러 sweep·piezo hysteresis가 섞인 파일은 서로
다른 branch를 하나로 평균한다. `Reverse sweep`은 이미 합쳐진 x축의 부호만 바꾸며
원래 sweep branch를 복원하지 않는다 (`streamlit_app.py:1049-1052`).

진단용 왕복 9-row trace는 5 unique point로 줄었고, forward/return의 서로 다른
값 네 쌍이 median으로 합쳐졌다. 현재 `CSVImportDiagnostics`에는 원본 x 방향전환
횟수나 monotonic run 수가 없다 (`gabes/experimental_csv.py:27-44`).

**권고:** 가장 싼 1단계는 sort 전에 `sign(diff(x))`의 유효 방향전환 수를 O(N)으로
세어 non-monotonic trace warning을 내는 것이다. 물리/성능 영향은 사실상 0이다.
그다음 필요할 때만 monotonic run별 overlay 또는 사용자가 고른 한 sweep을 표시한다.
자동 bin-average는 hysteresis를 숨길 수 있으므로 기본값으로 두지 않는 편이 안전하다.

## 새 lock/readout 문제

### feature 미검출인데 full-spectrum slope가 hero가 됨

species mode는 `narrowest_subdoppler()`가 `nan`을 반환해도
`_lock_readout_metrics()`를 호출하고 첫 metric을 hero로 만든다
(`gabes/schemes/sas.py:365-386`). helper는 feature가 유효하지 않으면 의도적으로
full-spectrum gradient의 최대값을 찾는다
(`gabes/schemes/sas.py:508-552`). 현재 테스트도 미검출 조건에서 `Lock Slope`
hero를 요구한다 (`tests/test_sas.py:215-220`).

테스트와 같은 조건

```text
pump = 0.01 mW
temperature = 200 °C
cell = 200 mm
```

에서 오늘 얻은 결과는 다음과 같다.

| 항목 | 결과 |
|---|---:|
| sub-Doppler feature | `unresolved` |
| Peak OD | `19881.44` |
| hero Lock Slope | `0.0000 /MHz` |
| Lock Detuning | `+6886.3 MHz` |
| 가장 가까운 hyperfine marker까지 거리 | `2231.5 MHz` |

이는 “feature가 없으므로 lock 후보를 특정할 수 없음”이라는 상황에서 broad
envelope/scan edge의 수치가 가장 중요한 실험 readout으로 승격된 경우다. 특히
`0.0000 /MHz`와 2.23 GHz marker distance는 실제 lock discriminator 기준으로
유용하지 않다.

**권고:** feature가 finite하지 않으면 `SAS status: sub-Doppler unresolved`를 hero로
되돌리고, full-spectrum slope는 숨기거나 `Broad-envelope slope (not a SAS lock)`
같은 secondary diagnostic으로 명확히 낮춘다. 추가 solve는 없으며 이미 계산한
O(N) gradient의 label/metric hierarchy만 바뀐다. feature가 검출된 정상 경로의
한-FWHM local search는 그대로 유지해야 한다.

### 문서가 약속하는 Doppler FWHM과 현재 readout이 불일치

사용자 가이드는 OD/SAS 주요 출력에 `Doppler FWHM`을 포함한다
(`docs/Userguide/GABES_User_Guide_v2.html:523-529`). compute는 여전히 thermal
Doppler FWHM을 계산해 raw에 보관한다 (`gabes/schemes/sas.py:168-203`). 그러나
현재 observables는 이를 metric으로 내지 않으며, 테스트가 부재를 명시적으로
고정한다 (`tests/test_sas.py:189-197`, `tests/test_sas.py:240-256`).

이 값이 단일 isotope thermal width의 최대값이라는 의미를 명확히 할 수 있다면
secondary metric으로 되살리는 편이 실험가에게 유용하다. 아니면 guide를 현재
readout에 맞춰 고쳐야 한다. 계산값은 이미 있으므로 runtime overhead는 O(1)이다.

## 기존 개선 제안과 계산비용 재평가

| 기존 제안 | 물리적 가치 | 계산비용 | 오늘의 판단 |
|---|---|---|---|
| gas/species/line별 pressure broadening + shift | buffer-cell line centre와 폭 | table lookup, scalar offset | **거의 0**, coefficient 검증 후 우선 적용 |
| phenomenological Dicke narrowing | 단순 pressure broadening 편향 완화 | 동일 solver dimension | **낮음**, 적용 범위를 명시 |
| added homogeneous width를 pump optical coherence에도 반영 | saturation/hole-burning 일관성 | 같은 Liouvillian 크기·solve 수 | **낮음, correctness 우선** |
| sub-Doppler resolution status + edge interpolation | 과도한 숫자 정밀도 방지 | 추가 solve 없이 O(N) | **거의 0, 즉시 가능** |
| species/line/polarization별 power-to-Rabi calibration | pump mW의 종간/선간 비교 | table lookup | **거의 0**, 검증 자료가 필요 |
| low-order polarization/Zeeman proxy | 비대칭·contrast·branching 설명 | scalar/저차 보정이면 낮음 | 낮지만 임의 fit knob가 되지 않게 주의 |
| full VCC | buffer-rich cell의 velocity redistribution | velocity class 결합, block/iterative solve | **높음**, 별도 heavy mode/설계 필요 |

### pump OBE의 homogeneous broadening 불일치

species path는
`gamma_eff = gamma_nat + self broadening + Ne broadening`을 계산하지만
(`gabes/schemes/sas.py:163-179`), 이 값은 pump table 간격과 probe Lorentzian에
사용될 뿐이다 (`gabes/schemes/sas.py:219-246`). pump Liouvillian의 atom은
natural `gamma`로 먼저 만들어지고 추가 optical dephasing이 없다
(`gabes/species.py:349-399`).

2026-07-20 진단에서는 같은 density-matrix 차원에 added FWHM/2의 optical
pure-dephasing surrogate를 넣었을 때 Rb-85 D2, 1.5 mW, 30 °C의 최대 alpha
변화가 Ne 5 Torr에서 27.9%, 20 Torr에서 38.8%였고 runtime은 사실상 같았다
(`docs/daily_report/2026-07-20_scheme-1_od-sas.md:74-82`). 이 수치는 정답 계수가
아니라 현재 누락에 대한 민감도다. spontaneous population decay는 natural Γ로
유지하고, 검증된 collisional optical-dephasing convention을 별도로 넣어야 한다.

### sub-Doppler readout의 해상도 의존성

`narrowest_subdoppler()`는 Python running median을 뺀 residual의 half-height
경계를 sample 단위로 찾고 edge interpolation을 하지 않는다
(`gabes/lineshape.py:107-130`). 오늘 동일 기본 조건의 현재 코드에서 재측정했다.

| scan points | sample 간격 | FWHM | Lock Slope | Lock Detuning | compute |
|---:|---:|---:|---:|---:|---:|
| 1401 | `8.057 MHz` | `32.23 MHz` | `0.0241 /MHz` | `-1304.9 MHz` | `1.156 s` |
| 4001 | `2.820 MHz` | `22.56 MHz` | `0.0316 /MHz` | `-1301.7 MHz` | `3.134 s` |

1401→4001에서 FWHM은 약 30%, slope는 약 31% 바뀌며 compute는 2.71배
느려졌다. 따라서 기본값을 무조건 4001로 올리기보다, 먼저
`samples_per_width`, `resolution-limited`, half-height linear interpolation을
추가하는 것이 좋다. 이 상태 표시와 보간은 추가 OBE solve 없이 가능하다.
publication용일 때만 local refinement를 opt-in으로 돌리는 편이 낫다.

### 고정 `I_sat`

모든 isotope와 D1/D2의 pump Rabi가 동일한
`I_sat = 4.484 mW/cm²`를 사용한다
(`gabes/species.py:273-283`, `gabes/schemes/sas.py:134-137`). 그래서 pump-power
trend는 유용하지만 서로 다른 species/line의 mW를 절대 saturation 기준으로
직접 비교하면 안 된다. validated line/polarization convention별 lookup을
도입하는 runtime 비용은 사실상 0이지만, 근거 없는 보정 knob를 추가해서는 안 된다.

## 동작을 바꾸지 않는 순수 코드 최적화

### 1. pump-off analytic population fast path

`pump_power_mw=0`이면 펌프 population factor는 모든 detuning에서 정확히
`w=1`이다 (`gabes/schemes/sas.py:27-30`). 현재도 Hamiltonian/Liouvillian과
`_pump_pops()` table을 만들고 보간한다
(`gabes/schemes/sas.py:215-246`, `gabes/schemes/sas.py:498-505`).

현재 코드에서 85Rb D1, 90 °C, 12.5 mm, 1401점으로 다시 측정한 결과:

| 경로 | median runtime |
|---|---:|
| 현재 pump-off | `0.5840 s` |
| analytic `w=1` | `0.1804 s` |
| speedup | **3.24×** |

두 scan axis는 정확히 같았고, 최대 상대 alpha 차이는 `5.63×10⁻¹⁶`이었다.
AutoOD, natural-Rb overlay, 여러 pressure/temperature 회귀로 보호하면 가장
안전하고 효과가 큰 pure optimization이다.

### 2. running median vectorization

`narrowest_subdoppler()`는 각 sample에서 Python list comprehension으로 median을
계산한다 (`gabes/lineshape.py:113-119`). 동일 padding/window를
`sliding_window_view`와 axis median으로 바꾼 1401점 진단은:

| 구현 | median runtime |
|---|---:|
| 현재 loop | `70.67 ms` |
| vectorized | `0.134 ms` |
| speedup | **529×** |
| 최대 차이 | `0` |

전체 OBE solve가 아니라 observables 단계의 최적화지만, headless batch와
interactive readout을 그대로 보존하며 줄일 수 있는 비용이다.

### 3. duplicate-heavy CSV median merge vector화

현재 `_sort_and_merge()`는 duplicate group마다 Python loop에서 `np.median`을
호출한다 (`gabes/experimental_csv.py:310-324`). 500,000행, 각 x가 5번 반복되는
synthetic scope trace에서, 한 번의 `(x,y)` lexsort 뒤 각 group의 중앙 index를
벡터로 선택하는 동등 구현을 비교했다.

| 구현 | runtime |
|---|---:|
| 현재 group loop | `5.097 s` |
| vectorized median indices | `0.569 s` |
| speedup | **8.96×** |
| 최대 merged-y 차이 | `0` |

CSV row limit 자체가 500,000이므로 실제 상한 workload에 의미가 있다
(`gabes/experimental_csv.py:17-19`, `gabes/experimental_csv.py:277-280`).
odd/even duplicate 수, unsorted input, signed zero, large finite values를 회귀한 뒤
적용할 수 있는 pure optimization이다.

### 4. 낮은 우선순위

- `Show unfiltered trace`일 때 `transformed_detuning()`을 corrected/raw용으로 두 번
  같은 인자로 호출한다 (`streamlit_app.py:1049-1063`). `aligned_x`를 재사용하면
  O(N) 배열 계산 하나를 없앨 수 있으나 위 세 항목보다 효과는 작다.
- pump-on의 interpolation/transition accumulation fused kernel은 여전히 주된
  compute 병목 후보지만, 단순 NumPy index arithmetic이 오히려 느려질 수 있다는
  이전 진단이 있다. 먼저 위의 exact fast paths를 적용하고, 그다음 별도 benchmark와
  strict tolerance test 아래에서만 검토하는 편이 안전하다.

## 검증 결과

Scheme 1 관련:

```text
python -m pytest tests/test_sas.py tests/test_experimental_csv.py \
  tests/test_absorption.py tests/test_headless_observables.py \
  tests/test_schemes_render.py -q

80 passed in 41.28s
```

전체 저장소:

```text
MPLBACKEND=Agg python -m pytest -q

232 passed in 109.92s
```

테스트는 현재 구현 동작을 잘 보호하지만, 일부는 오늘 지적한 trust boundary도
의도적으로 고정한다. 특히 feature 미검출인데 `Lock Slope`를 hero로 유지하는
테스트와 Doppler FWHM 부재를 요구하는 테스트는 물리 reference 관점에서
재검토 대상이다 (`tests/test_sas.py:189-197`, `tests/test_sas.py:215-256`).

## 최종 우선순위

1. **P1·무부하:** feature 미검출 시 numeric lock hero를 내리지 말고
   `sub-Doppler unresolved` status를 hero로 표시한다.
2. **P1/P2·무부하:** CSV 자동 0–1 결과를 `relative normalized`로 명시하고,
   absolute 비교에는 dark/reference calibration mode를 제공한다.
3. **P2·거의 무부하:** non-monotonic/bidirectional CSV scan을 O(N)으로 감지해
   branch collapse warning을 낸다.
4. **P2·무부하:** sub-Doppler `samples_per_width`/resolution status와 edge
   interpolation을 추가하고, guide와 Doppler FWHM readout 불일치를 해소한다.
5. **P2·저부하 correctness:** buffer-gas table/pressure shift와 함께 added
   homogeneous width를 pump optical dephasing에도 일관되게 반영한다.
6. **순수 성능:** pump-off analytic path, running-median vector화,
   duplicate-heavy CSV merge vector화를 적용한다.
7. **필요할 때만 heavy:** full VCC와 full Zeeman/polarization solver는 대상 셀,
   편광, 허용 runtime, 검증 데이터가 정해진 뒤 별도 mode로 설계한다.

종합하면 Scheme 1의 원자물리 코어는 실제 실험에 유용하다. 특히 pump-off
absolute OD anchor와 hyperfine-pumping SAS 구조는 강점이다. 현재 가장 값싼
신뢰도 향상은 더 큰 원자모델을 넣는 것이 아니라, **검출되지 않은 feature에는
lock 숫자를 주지 않고, 상대 정규화 CSV를 절대 transmission처럼 보이지 않게
하며, 이미 알려진 수치 해상도를 상태로 노출하는 것**이다. 이 세 가지는
interactive 성능과 현재 물리를 보존하면서 실험 reference의 정직성을 크게 높인다.
