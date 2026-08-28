# 2026-08-05 Scheme 1 물리 리뷰 — OD / SAS

## 선택 규칙과 현재 다섯 scheme

- 로컬 날짜는 `2026-08-05`, day-of-month는 `5`이다.
- `n = (day mod 5) + 1 = (5 mod 5) + 1 = 1`이므로 오늘의 대상은 **Scheme 1**이다.
- UI 순서는 `gabes/schemes/__init__.py`의 `_SCHEMES`가 정하며, README의 표와
  일치한다 (`gabes/schemes/__init__.py:12-25`, `README.md:8-16`).

| 번호 | 현재 scheme | registry 구현 | 주된 물리/출력 |
|---:|---|---|---|
| 1 | OD / SAS | `SASScheme()` | pump-off Doppler OD, pump-on Lamb dip/crossover |
| 2 | Lambda coherence (EIT / AT / CPT) | `LambdaScheme()` | 3준위 결맞음, EIT/AT/CPT |
| 3 | Rydberg-EIT electrometry | `RydbergEITScheme()` | cascade EIT, microwave AT, finite-IF electrometry |
| 4 | Hanle / EIA / NMOR | `MagnetoScheme()` | 영자기장 투과와 자기광학 회전 |
| 5 | FWM | `FWMScheme()` | seeded gain/squeezing, SFWM biphoton |

Scheme 1의 내부 이름은 `sas`, UI 제목은 `Absorption spectroscopy (OD / SAS)`이다
(`gabes/schemes/sas.py:53-63`).

## 먼저 검색한 기존 제안과 변경 이력

물리 판단 전에 다음 자료를 검색했다.

- 저장소 지침/개요: `CLAUDE.md`, `README.md`
- 계획/TODO: `docs/checklist.json`, `gabes/constants.py:79-89`
- 이전 Scheme 1 보고서:
  `docs/daily_report/2026-06-25_scheme-1_od-sas.md`,
  `2026-07-10_scheme-1_od-sas.md`,
  `2026-07-20_scheme-1_od-sas.md`,
  `2026-07-30_scheme-1_od-sas.md`
- 구현: `gabes/schemes/sas.py`, `gabes/species.py`, `gabes/core.py`,
  `gabes/doppler.py`, `gabes/lineshape.py`, `gabes/experimental_csv.py`,
  `streamlit_app.py`
- 테스트: `tests/test_sas.py`, `tests/test_experimental_csv.py`,
  `tests/test_absorption.py`, `tests/test_headless_observables.py`,
  `tests/test_schemes_render.py`
- 예제/설명: `docs/Userguide/GABES_User_Guide_v2.html:518-573`,
  `docs/Userguide/userguide_assets/od.png`,
  `docs/Userguide/userguide_assets/sas.png`, `README.md:81-100`
- 공개 GitHub Issues는 2026-08-05 GitHub API 확인 기준 open/closed 합계 `0`건이다.
  로컬에도 별도 issue-note 파일은 없었다.

기존 개선 제안은 **있다**.

- gas/species/line별 buffer-gas broadening/shift와 phenomenological Dicke
  narrowing (`docs/checklist.json:50-54`)
- OD/SAS용 저차 polarization/Zeeman proxy
  (`docs/checklist.json:116-120`)
- velocity class를 결합하는 full VCC kernel
  (`docs/checklist.json:130-134`)
- 이전 일일 보고서의 pump optical-dephasing 일관성, sub-Doppler resolution
  status/interpolation, species/line별 power-to-Rabi calibration, CSV 절대 보정,
  sweep-branch 경고, pump-off fast path, running-median/duplicate-merge vectorization
  (`docs/daily_report/2026-07-30_scheme-1_od-sas.md:277-432`)

체크리스트의 `_daily_report_sync`는 아직 “2026-07-09까지 반영”이라고 되어 있어
(`docs/checklist.json:3`), 7월 20일/30일의 Scheme 1 신뢰성 제안은 checklist item으로
승격되지 않은 상태다. `sas.py`, `species.py`, `experimental_csv.py`와 관련 테스트에는
2026-07-30 검토 이후 commit이 없으며, 현재 uncommitted 변경도 Scheme 1 코어를
건드리지 않는다.

## 결론 요약

**OD/SAS는 실제 원자물리를 구현하며, pure warm-vapor cell의 선 배정, 절대 OD
sanity check, 온도/펌프/통과시간 경향, SAS lock 후보 탐색에 유용한
semi-quantitative 실험 레퍼런스다.** 같은 hyperfine manifold에서 pump-off OD와
pump-on SAS가 연속적으로 이어지고, Wigner-6j/3j 선세기, CG-branched decay,
hyperfine optical pumping, Maxwell 속도 평균을 포함한다
(`gabes/schemes/sas.py:4-30`, `gabes/species.py:337-422`).

다만 다음 용도의 절대 reference로는 아직 부적합하다.

- buffer-cell의 pressure-shifted line centre와 Dicke/VCC linewidth
- 편광/Zeeman-resolved SAS contrast와 절대 saturation power
- 고온 Cs의 self-broadened linewidth
- 변조·복조·검출기 전달함수를 포함한 실제 servo error signal
- 자동 0–1 CSV overlay에서의 absolute transmission/OD
- feature가 미검출되었거나 격자에 충분히 분해되지 않은 조건의 lock slope/FWHM

오늘의 새 핵심 발견은 두 가지다.

1. `species.py`가 **Rb용이라고 명시한 self-broadening 계수**를 Cs에도 그대로
   적용한다. 권장 저온 Cs 프리셋에서는 무시할 수 있지만, 허용 slider의 고온부에서는
   추가폭이 자연선폭보다 훨씬 커져 Cs 정량 reference를 훼손한다.
2. 공개 문서는 85Rb D1 AutoOD 일치를 `<0.1 %`라고 말하지만, 해당 회귀 테스트는
   `<1 %`만 요구한다. 현재 값은 실제로 0.1% 안에 있지만, 문서 약속을 테스트가
   충분히 보호하지 않는다.

두 항목 모두 solver 차원이나 solve 수를 늘리지 않고 해결할 수 있다.

## 현재 구현이 담는 실제 물리

### 1. pump-off OD와 pump-on SAS가 한 모델이다

이상적인 역행 pump/probe에 대해 원자 좌표계 detuning은 각각 `Δ + kv`, `Δ - kv`다.
펌프를 끄면 population factor가 1인 Doppler-Voigt 흡수로 축약되고, 펌프를 켜면
velocity-selective saturation/hole burning과 optical pumping이 같은 스펙트럼에
Lamb dip/crossover를 만든다 (`gabes/schemes/sas.py:4-30`,
`gabes/schemes/sas.py:205-247`).

온도, cell length, isotope/natural abundance, D1/D2, pump power, waist, transit
relaxation, Ne pressure를 실험 단위로 바꿀 수 있다
(`gabes/schemes/sas.py:65-107`). `cell_mm`와 `line_strength`는 post-process knob로
분리되어 불필요한 OBE 재계산도 피한다 (`gabes/schemes/sas.py:82-99`,
`gabes/schemes/sas.py:318-325`).

### 2. 실제 alkali hyperfine data와 optical pumping

Rb-85, Rb-87, Cs-133의 hyperfine A/B, centroid frequency, mass, natural
linewidth, abundance가 명시되어 있다 (`gabes/species.py:127-182`). 허용
`Fg→Fe` 선세기는 Wigner-6j/3j로 계산된다 (`gabes/species.py:75-103`,
`gabes/species.py:226-270`). Ground degeneracy를 포함한 동일한
`T(Fg,Fe)`가 흡수 weight, spontaneous branching, 상대 pump Rabi에 들어간다
(`gabes/species.py:366-414`).

`build_manifold()`는 모든 허용 hyperfine 전이와 CG-branched decay를 만들고,
모든 상태를 thermal ground distribution으로 되돌리는 transit relaxation으로
무한 optical pumping을 regularize한다 (`gabes/species.py:337-399`). 실제 alkali
SAS에서 중요한 enhanced/inverted crossover를 단일 Lorentzian 합보다 훨씬
현실적으로 재현하는 부분이다.

단, 한 `F` manifold를 하나의 lumped state로 취급하므로 mF, laser polarization,
잔류 자기장, beam-profile diffusion은 풀지 않는다. 따라서 line centre/큰 구조와
경향은 유용하지만 absolute contrast/branching calibration에는 한계가 있다.

### 3. 절대 OD 기준과 테스트

약한 probe의 line-integrated absorption은 AutoOD convention의 reduced dipole,
`C_F²`, CRC Rb vapor density를 사용한다 (`gabes/species.py:193-213`,
`gabes/species.py:226-270`). 테스트는 다음을 보호한다.

- 85Rb D1 `C_F²`, dipole, CRC density (`tests/test_sas.py:58-68`)
- known hyperfine splitting과 manifold line count (`tests/test_sas.py:71-84`)
- spontaneous branching의 값과 합계 (`tests/test_sas.py:87-98`)
- pump-off AutoOD scale와 49/25 manifold ratio (`tests/test_sas.py:101-127`)
- pump-on sharp feature와 transit-rate에 따른 hyperfine-pumping crossover
  (`tests/test_sas.py:129-154`)
- natural-Rb isotope overlay, species/line preset, generic fallback, figure/headless
  contract (`tests/test_sas.py:157-294`, `tests/test_headless_observables.py:33-64`)

오늘 85Rb D1, 90 °C, 12.5 mm, pump off 조건에서 내부 AutoOD validation primitive와
비교한 차이는 integrated area `0.0393 %`, peak `0.0126 %`였다. 따라서 현재
`<0.1 %` 문구 자체는 재현된다. 그러나 테스트 허용치는 area/peak 모두 `<1 %`다
(`tests/test_sas.py:101-111`). 반면 app/README는 `<0.1 %`를 약속한다
(`gabes/schemes/sas.py:128-130`, `README.md:252-257`). 같은 계산의 threshold를
0.1%에 맞추면 runtime 추가 없이 문서 약속을 직접 보호할 수 있다.

### 4. 문서와 예제의 참고 가치

사용자 가이드는 동일한 natural-Rb D1 조건의 pump-off `od.png`와 pump-on
`sas.png`를 나란히 제시해 물리적 차이를 직관적으로 보여 준다
(`docs/Userguide/GABES_User_Guide_v2.html:549-572`). 다만 두 PNG는 2026-06-08에
생성된 stacked Transmission/OD 레이아웃이고, 현재 app은 두 figure view를
분리한다 (`gabes/schemes/sas.py:336-363`). SAS 그림의 매우 좁은 spike도 아래의
격자 제한을 함께 명시하지 않는다. 따라서 예제는 정성적 교육 자료로는 좋지만,
현재 metric/layout의 재현 예제로 보려면 새 provenance와 resolution status로
재생성하는 편이 안전하다.

## 오늘 새로 확인한 trust 문제

### 1. Rb self-broadening 계수가 Cs에도 무조건 적용됨

`BETA_SELF`의 주석은 명시적으로 “Rb value”라고 되어 있다
(`gabes/species.py:49-50`). 그런데 `self_broadened_gamma(iso, N)`는 `iso`를 전혀
사용하지 않고 `βN`만 반환하며 (`gabes/species.py:211-213`), Scheme 1은 모든
isotope component에 이를 적용한다 (`gabes/schemes/sas.py:169-179`). Cs도 예외가
아니다.

현재 코드가 Cs에 더하는 `β_Rb N_Cs / 2π`는 다음과 같다.

| Cs 온도 | 현재 추가 self width | 비교 |
|---:|---:|---|
| 22 °C | `0.00219 MHz` | D2 권장점, 자연폭 5.2227 MHz에 비해 무시 가능 |
| 30 °C | `0.00478 MHz` | D1 권장점, 자연폭 4.5612 MHz에 비해 무시 가능 |
| 100 °C | `1.0186 MHz` | 더는 정밀 linewidth에서 무시하기 어려움 |
| 150 °C | `15.175 MHz` | 자연폭보다 큼 |
| 200 °C | `124.010 MHz` | 허용 slider 끝, 자연폭을 압도 |

온도 slider는 모든 species에 20–200 °C를 허용한다
(`gabes/schemes/sas.py:80-81`). 따라서 기본 Cs 그림은 거의 영향을 받지 않지만,
고온 Cs 스윕의 linewidth/contrast를 현재 값으로 정량 해석해서는 안 된다.

**권고:** self-broadening도 `(species, line)` coefficient table로 만들고, 검증된 Cs
계수가 없으면 0으로 조용히 가정하지 말고 `unsupported/qualitative` status를 낸다.
Rb natural mixture에서는 `number_density()`가 elemental total density를 반환한 뒤
caller가 isotope abundance를 곱하므로 (`gabes/species.py:193-207`,
`gabes/schemes/sas.py:169-176`), `βN`에 total Rb와 isotope-partial density 중 무엇을
쓸지도 문헌 convention으로 명시해야 한다. 둘 다 table lookup/scalar multiply라
계산 오버헤드는 사실상 0이다.

### 2. 문서의 AutoOD 허용오차보다 회귀가 10배 느슨함

현재 구현은 오늘 측정에서 `<0.1 %`를 만족하지만, `test_pump_off_reproduces_autood_85rb_d1`
은 `<1 %`만 요구한다 (`tests/test_sas.py:101-111`). 향후 0.5% drift가 생겨도 테스트는
통과하면서 공개 설명과 모순될 수 있다.

**권고:** 현재 검증 조건에서 area/peak tolerance를 `1e-3` 이하로 맞추고, 가능하면
`references/AutoOD/ReferenceOD.csv`의 조건/축 convention을 별도 provenance test로
고정한다. 기존 solve를 그대로 쓰므로 추가 runtime은 없거나 매우 작다.

## 재확인된 기존 신뢰 경계와 계산비용

### 1. sub-Doppler readout은 아직 grid-limited다

`narrowest_subdoppler()`는 running-median residual의 half-height 경계를 sample 단위로
찾으며 edge interpolation을 하지 않는다 (`gabes/lineshape.py:107-130`). 현재 기본
조건을 재측정한 결과다.

| scan points | sample 간격 | FWHM | Lock Slope | Lock Detuning | warm compute |
|---:|---:|---:|---:|---:|---:|
| 1401 | `8.057 MHz` | `32.23 MHz` | `0.0241 /MHz` | `-1304.9 MHz` | `0.533 s` |
| 4001 | `2.820 MHz` | `22.56 MHz` | `0.0316 /MHz` | `-1301.7 MHz` | `1.289 s` |

점 수를 2.86배로 늘리자 FWHM은 약 30%, slope는 약 31% 바뀌었다. 기본점을 무조건
늘리기보다 `samples_per_width`, `resolution-limited` status, half-height linear
interpolation을 먼저 넣는 것이 낫다. 이는 O(N) post-process이고 추가 OBE solve가
없다. publication용 local refinement만 opt-in으로 두면 된다.

사용자 가이드는 주요 출력에 `Doppler FWHM`을 적지만
(`docs/Userguide/GABES_User_Guide_v2.html:523-529`), compute가 보관한 값
(`gabes/schemes/sas.py:168-203`)을 현재 observables는 내지 않으며 테스트가 그 부재를
고정한다 (`tests/test_sas.py:189-197`, `tests/test_sas.py:240-256`). 계산값은 이미
있으므로 metric을 복구하거나 guide를 고치는 비용은 O(1)이다.

### 2. feature 미검출인데 full-spectrum slope가 hero다

feature가 `nan`이어도 `_lock_readout_metrics()`는 full-spectrum gradient 최대값을
찾고, 첫 metric을 hero로 만든다 (`gabes/schemes/sas.py:365-386`,
`gabes/schemes/sas.py:508-552`). 테스트도 이 동작을 요구한다
(`tests/test_sas.py:215-220`).

재현 조건 `pump=0.01 mW`, `T=200 °C`, `L=200 mm`에서 결과는 다음과 같다.

- `SAS status = sub-Doppler unresolved`
- `Peak OD = 19881.44`
- hero `Lock Slope = 0.0000 /MHz`
- `Lock Detuning = +6886.3 MHz`
- 가장 가까운 hyperfine marker와 거리 `2231.5 MHz`

feature가 없다는 status를 hero로 올리고, full-envelope slope는 숨기거나 secondary로
내리는 것이 물리적으로 정직하다. 이미 계산한 O(N) gradient의 hierarchy/label만
바꾸므로 solve overhead는 0이다.

### 3. CSV overlay는 상대 line-shape 비교다

CSV importer는 robust parsing, duplicate merge, Hampel filter, x-aware local
quadratic smoothing을 제공하고 좁은 feature를 지우지 않도록 방어한다
(`gabes/experimental_csv.py:127-190`, `gabes/experimental_csv.py:347-555`). 실험실
scope trace를 빠르게 정리하는 도구로는 유용하다.

그러나 `_normalize_transmission()`은 trace 자체의 extrema를 floor/ceiling으로 잡아
항상 `[0,1]`로 affine map한다 (`gabes/experimental_csv.py:558-641`;
`tests/test_experimental_csv.py:83-95`). 실제 dark/reference calibration이 아니므로
물리적 `T=0.70–0.92`도 `0–1`이 된다. UI도 이를 단순 “0-1 calibration”과 “Signal
calibration”으로 부른다 (`streamlit_app.py:899-905`, `streamlit_app.py:1031-1040`).

또 `_sort_and_merge()`는 취득 순서를 버리고 x로 정렬한 뒤 동일 x의 y를 median으로
합친다 (`gabes/experimental_csv.py:310-324`). 왕복 scan/hysteresis branch가 한 곡선으로
붕괴할 수 있고, UI의 `Reverse sweep`은 이미 합친 x축의 부호만 바꾼다
(`streamlit_app.py:1002-1011`, `streamlit_app.py:1049-1052`).

가장 싼 개선은 다음 두 가지다.

- label을 `relative normalized transmission`으로 바꾸고, absolute 모드에는
  `(dark, reference)` 또는 `(gain, offset/I0)`를 받는다: O(N), OBE overhead 0.
- sort 전에 `sign(diff(x))` 방향전환 수를 세어 non-monotonic warning을 낸다:
  O(N), 메모리/시간 미미. 실제 branch overlay는 필요할 때만 추가한다.

### 4. added homogeneous width가 pump OBE에는 들어가지 않는다

species path는 `gamma_eff = gamma_nat + self + Ne broadening`을 계산하지만
(`gabes/schemes/sas.py:163-179`), pump Liouvillian은 natural-Γ manifold로 먼저
만들어진다 (`gabes/species.py:349-399`, `gabes/schemes/sas.py:215-221`).
`gamma_eff`는 pump table grid와 probe Lorentzian에만 쓰인다
(`gabes/schemes/sas.py:219-246`).

이전 진단의 optical pure-dephasing surrogate는 Rb-85 D2에서 max alpha를 Ne 5 Torr
`27.9 %`, 20 Torr `38.8 %` 바꿨지만 solver 차원과 runtime은 거의 같았다
(`docs/daily_report/2026-07-20_scheme-1_od-sas.md:74-82`). 이 민감도는 정답 계수가
아니며, natural population decay와 collisional optical dephasing을 분리한 문헌
convention이 필요하다. 구현 자체는 같은 Liouvillian 크기/solve 수라 저비용이다.

### 5. 모든 species/line이 같은 Rb 기준 `I_sat`를 쓴다

`pump_rabi_from_power()`는 모든 Rb/Cs D1/D2에 `I_sat = 4.484 mW/cm²`를 쓴다
(`gabes/species.py:273-283`, `gabes/schemes/sas.py:134-137`). 따라서 각 조건 안의
pump-power trend는 유용하지만 종/선간 절대 saturation power 비교는 아직 정량
기준이 아니다. polarization convention이 붙은 검증된 lookup table은 O(1)이지만,
검증 데이터 없는 자유 fit knob를 늘려서는 안 된다.

## 기존 개선안의 비용/물리 보존성 평가

| 개선안 | 물리적 가치 | 계산비용 | 판단 |
|---|---|---|---|
| species/line별 self-broadening table/status | 고온 Cs/자연 Rb linewidth 신뢰 | O(1) lookup | **즉시 가능한 새 P1/P2 trust fix** |
| AutoOD `<0.1 %` 회귀 강화 | 공개 검증 주장 보호 | 같은 solve | **즉시 가능** |
| gas/species/line buffer broadening + shift | buffer-cell line centre/폭 | table lookup + scalar offset | **거의 0**, coefficient 검증 필요 |
| phenomenological Dicke narrowing | 단순 pressure broadening 편향 완화 | 동일 solver dimension | **낮음**, 유효범위 표기 필요 |
| added width의 pump optical dephasing | saturation/hole-burning 일관성 | 동일 차원·동일 solve 수 | **낮음, correctness 우선** |
| resolution status/interpolation | 숫자 과신 방지 | O(N) post-process | **거의 0** |
| dark/reference CSV calibration | absolute contrast 보존 | O(N) affine transform | **OBE overhead 0** |
| 저차 polarization/Zeeman proxy | contrast/비대칭 설명 | scalar/transition weight면 낮음 | calibration 없으면 임의 knob 위험 |
| full VCC | collision-rich velocity redistribution | velocity classes 결합 | **높음**, 별도 heavy mode/설계 필요 |
| full mF/polarization manifold | absolute contrast/branching | density-matrix 차원 급증 | **높음**, 검증 데이터와 범위 합의 필요 |

저부하 제안은 현재 핵심 물리와 AutoOD anchor를 보존할 수 있다. 특히 status/label,
table lookup, post-process calibration은 기존 OBE 결과를 바꾸지 않는다. Pump optical
dephasing과 Dicke proxy는 결과를 의도적으로 바꾸므로 feature flag, 기존 pure-cell
회귀, 문헌 계수 provenance가 필요하다.

## 동작을 바꾸지 않는 순수 코드 최적화

### 1. pump-off analytic population fast path

펌프가 0이면 population factor는 정확히 `w=1`인데도 현재 Hamiltonian/Liouvillian,
`_pump_pops()` table과 보간을 수행한다 (`gabes/schemes/sas.py:215-246`,
`gabes/schemes/sas.py:498-505`). 85Rb D1, 90 °C, 12.5 mm, 1401점의 warm median:

| 경로 | runtime |
|---|---:|
| 현재 | `0.2420 s` |
| analytic `w=1` | `0.0810 s` |
| speedup | **2.99×** |

scan axis는 bitwise 동일했고 max relative alpha difference는 `5.63×10⁻¹⁶`이었다.
AutoOD/natural-Rb/pressure/temperature 회귀로 보호하면 가장 가치가 큰 exact fast
path다.

### 2. running median vectorization

`narrowest_subdoppler()`는 각 sample마다 Python loop에서 median을 구한다
(`gabes/lineshape.py:113-119`). 동일 padding/window를
`sliding_window_view(...); median(axis=1)`로 바꾼 1401점 benchmark:

| 구현 | runtime |
|---|---:|
| 현재 loop | `21.64 ms` |
| vectorized | `0.160 ms` |
| speedup | **약 135×** |
| 최대 차이 | `0` |

전체 solve보다 작지만 headless sweep/readout에는 의미 있고 동작을 그대로 보존한다.

### 3. duplicate-heavy CSV median merge vectorization

`_sort_and_merge()`는 duplicate group마다 Python loop에서 `np.median`을 호출한다
(`gabes/experimental_csv.py:310-324`). 500,000행, x당 5개 중복 synthetic trace에서
`lexsort((y,x))` 후 group 중앙 index를 벡터로 선택한 동등 구현의 median은 다음과
같았다.

| 구현 | runtime |
|---|---:|
| 현재 group loop | `1.669 s` |
| vectorized median index | `0.1459 s` |
| speedup | **11.44×** |
| 최대 merged-y 차이 | `0` |

CSV의 실제 row limit가 500,000이므로 (`gabes/experimental_csv.py:17-19`) 상한
workload에 직접 의미가 있다. odd/even group, signed zero, large finite values를
회귀해야 한다.

### 4. 낮은 우선순위

- raw overlay에서 같은 인자의 `transformed_detuning()`을 두 번 호출한다
  (`streamlit_app.py:1049-1064`). `aligned_x` 재사용은 O(N) 배열 계산 하나를
  없애지만 위 세 항목보다 영향이 작다.
- pump-on의 uniform `deff` grid에 대해 `np.interp`의 index/weight를 level마다 다시
  찾는다 (`gabes/schemes/sas.py:231-246`). 공통 interpolation geometry를 재사용할
  여지는 있으나, strict equivalence/memory benchmark 뒤에 진행해야 한다.
- `build_manifold()`는 이미 `lru_cache(maxsize=64)`이므로
  (`gabes/species.py:337-338`) 추가 angular-data cache는 우선순위가 낮다.

## 검증 결과

Scheme 1 관련 테스트:

```text
python -m pytest tests/test_sas.py tests/test_experimental_csv.py \
  tests/test_absorption.py tests/test_headless_observables.py \
  tests/test_schemes_render.py -q

80 passed in 22.76s
```

전체 저장소 테스트:

```text
MPLBACKEND=Agg python -m pytest -q

299 passed in 127.96s
```

테스트는 현재 코어를 잘 보호하지만, 일부는 오늘 지적한 trust boundary를 그대로
고정한다. 특히 unresolved feature의 numeric lock hero
(`tests/test_sas.py:215-220`), Doppler FWHM 부재
(`tests/test_sas.py:189-197`, `tests/test_sas.py:240-256`), AutoOD 1% tolerance
(`tests/test_sas.py:101-111`)는 실험 reference 관점에서 재검토 대상이다.

## 최종 우선순위

1. **P1/P2·무부하:** Rb `β`를 Cs에 적용하지 않도록 species/line별
   self-broadening provenance를 만들고, 미검증 범위에는 status를 표시한다.
2. **P2·무부하:** AutoOD `<0.1 %` 공개 주장과 회귀 tolerance를 일치시킨다.
3. **P1·무부하:** feature 미검출 시 numeric lock hero 대신
   `sub-Doppler unresolved`를 hero로 표시한다.
4. **P1/P2·무부하:** CSV 자동 결과를 `relative normalized`로 명시하고,
   absolute 모드에는 dark/reference calibration을 제공한다.
5. **P2·거의 무부하:** bidirectional sweep 경고와 sub-Doppler
   `samples_per_width`/interpolation/status를 추가한다.
6. **P2·저부하 correctness:** 검증된 buffer/self-broadening table과 함께 added
   homogeneous width를 pump optical dephasing에도 일관되게 반영한다.
7. **순수 성능:** pump-off analytic path, running-median vector화,
   duplicate-heavy CSV merge vector화를 적용한다.
8. **필요할 때만 heavy:** full VCC와 full mF/Zeeman/polarization solve는 대상 셀,
   편광, 허용 runtime, 검증 데이터가 정해진 뒤 별도 mode로 설계한다.

종합하면 Scheme 1의 강점은 **실제 hyperfine line data, AutoOD-anchored pump-off OD,
hyperfine-pumping SAS를 하나의 빠른 모델로 연결한 것**이다. 현재 가장 좋은 다음
단계는 solver를 무겁게 만드는 것이 아니라, 종별 linewidth coefficient의 provenance,
문서-테스트 허용오차, unresolved/resolution/CSV calibration status를 바로잡는 것이다.
이들은 interactive 성능을 사실상 유지하면서 실험 reference의 신뢰도를 크게 높인다.
