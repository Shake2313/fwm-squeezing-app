# 2026-09-02 Scheme 3 물리 검토 — Rydberg-EIT electrometry

## 1. 선택 규칙과 현재 다섯 scheme

- 현지 날짜: 2026-09-02 (Asia/Seoul), 일자 `2`
- `n = (2 mod 5) + 1 = 3`
- 오늘의 대상: **Scheme 3 — Rydberg-EIT electrometry**

드롭다운 순서는 `gabes/schemes/__init__.py:19-25`의 `_SCHEMES`가 결정한다. alias는
직접 호출·테스트용일 뿐 드롭다운 순서를 바꾸지 않는다(`gabes/schemes/__init__.py:27-39`).
README의 현재 기능 표도 같은 다섯 엔진을 같은 순서로 설명한다(`README.md:7-16`).

| 번호 | 등록 클래스 / registry name | 현재 정의 |
|---:|---|---|
| 1 | `SASScheme()` / `sas` | Absorption spectroscopy — OD / SAS |
| 2 | `LambdaScheme()` / `lambda` | Lambda coherence — EIT / AT / CPT |
| **3** | `RydbergEITScheme()` / `rydberg_eit` | **⁸⁵Rb cascade EIT / microwave AT electrometry** |
| 4 | `MagnetoScheme()` / `magneto` | Hanle / EIA / NMOR |
| 5 | `FWMScheme()` / `fwm` | Seeded FWM gain diagnostic / generic SFWM biphoton |

Scheme 3의 공개 이름, 제목과 cache/default version은 `gabes/schemes/rydberg.py:91-100`에
정의돼 있다.

## 2. 먼저 검색한 기존 개선 제안, TODO와 issue note

새 제안을 처음부터 만들 필요는 없었다. 다음을 먼저 확인했다.

- 2026-06-22부터 2026-08-27까지의 Scheme 3 일일 보고서 10개
- 통합 작업 등록부 `docs/checklist.json`
- Scheme 3 코드·테스트·분석 경로와 저장소 파일명의 `TODO`, `FIXME`, `issue`,
  `proposal`, `deferred`, `parked` 검색
- `references/rydberg_ju_experiment_context.md`의 실험·논문 문맥
- 공개 GitHub 저장소의 open/closed issue 목록

별도의 Rydberg 전용 TODO/issue/proposal 파일과 코드 inline `TODO`/`FIXME`는 찾지 못했다.
2026-09-02에 GitHub API로 `Shake2313/fwm-squeezing-app`의 open/closed issue 전체를
조회한 결과도 0건이다. 그러나 기존 개선안은 아래처럼 `docs/checklist.json`에 이미 충분히
구체적으로 존재한다.

| 기존 항목 | 현재 상태 | 근거 위치 |
|---|---|---|
| ExtraView dependency와 cache correctness | P1 `ready`, medium | `docs/checklist.json:504-521` |
| effective opacity 공개와 sweep trust | P1 `ready`, medium | `docs/checklist.json:524-544` |
| 저투과·저대비·경계·해상도 readout guardrail | P1 `needs_decision`, medium | `docs/checklist.json:547-568` |
| held-out EIT/RF/PSD/temperature dataset | P2 `blocked_external`, large | `docs/checklist.json:572-593` |
| limited Zeeman/polarization reference | P2 `parked`, research | `docs/checklist.json:898-920` |
| measured readout transfer function | P2 `blocked_external`, medium | `docs/checklist.json:923-945` |
| provenance-tagged RF field-map import | P2 `blocked_external`, large | `docs/checklist.json:948-971` |

직전 Scheme 3 보고서의 검토 기준 HEAD는 `ce60b5c`였고, 현재 HEAD는 `c922ca2`다.
두 commit 사이 및 현재 working tree에서 Scheme 3의 source, test, cell-heating workflow,
checklist, User Guide에는 변경이 없다. 현재 관련 물리와 결함 판정은 따라서 그대로다.
사용자가 진행 중인 SABES/FWM 파일과 2026-09-01 보고서는 건드리지 않았다.

## 3. 현재 구현이 담은 실제 물리

### 3.1 축약 4준위 cascade OBE

광학 코어는 ⁸⁵Rb
`5S1/2 F=3 -> 5P3/2 F'=4 -> 40D5/2 -> 39F7/2`의 4준위 정상상태 OBE다.
780 nm probe, 481 nm counter-propagating coupling, 약 37 GHz RF LO를 Hamiltonian에 넣고,
5P spontaneous decay, 작은 Rydberg decay, 5S–40D와 40D–39F dephasing을 분리한다
(`gabes/schemes/rydberg.py:1-15`, `69-88`, `413-433`, `554-623`).

- probe와 coupling Rabi는 기준 실험점에 고정한 `Omega proportional to sqrt(P)/d` 스케일이다
  (`gabes/schemes/rydberg.py:39-53`, `435-457`). 이것은 실험 노브의 방향을 담지만 bare-state
  dipole만으로 얻은 first-principles 절대 Rabi는 아니다.
- heater setpoint, effective vapor temperature, cold spot을 분리한다. effective temperature는
  motion/transit/dephasing에, cold spot은 sealed-cell vapor pressure에 들어간다
  (`gabes/schemes/rydberg.py:459-501`, `559-580`).
- residual two-photon Doppler는 level별 wave-vector ratio를 쓰고, Doppler-on에서 Maxwell
  velocity class를 평균한다(`gabes/schemes/rydberg.py:69-88`, `503-523`, `606-615`).
- residual Zeeman은 full magnetic manifold가 아니라 susceptibility의 Gaussian convolution이다
  (`gabes/schemes/rydberg.py:617-623`).
- coherence는 D2 dipole, density, line factor와 Beer–Lambert 법칙으로 transmission에 변환된다
  (`gabes/schemes/rydberg.py:625-677`).

따라서 기준점 주변 EIT 폭, power/waist/transit trade-off, RF AT splitting과 microwave-detuning
방향을 보는 데 실제로 유용하다. 반면 특정 polarization, m-state population, stray B-field와
절대 optical depth를 예측하는 full atomic reference는 아니다.

### 3.2 finite-IF atomic superheterodyne와 detector chain

정적 AT 자나 slope proxy만 있는 것이 아니라, LO-dressed steady state 주변의 frequency-domain
small-signal chain이 구현돼 있다.

- `L0 +/- i omega_IF` sideband를 풀며 SIG phase와 peak angular-Rabi convention을 명시한다
  (`gabes/rydberg_electrometry.py:116-199`).
- Hermitian absorption quadrature와 velocity-class complex phasor의 coherent average를 쓴다
  (`gabes/rydberg_electrometry.py:69-113`, `216-232`;
  `gabes/schemes/rydberg.py:748-821`).
- RF dipole, angular/polarization factor와 peak/RMS field convention을 분리한다
  (`gabes/rydberg_electrometry.py:235-321`).
- 두 photodiode의 one-sided Schottky ASD, correlated RIN, electronics ASD를 합쳐 field ASD와
  ENBW noise로 환산한다(`gabes/rydberg_electrometry.py:362-530`).
- Doppler-off 기본 경로는 801개 probe detuning에서 total noise-equivalent field를 최소화한다
  (`gabes/schemes/rydberg.py:769-923`).

이는 atomic mixer의 선형·유한 IF 한계를 계산하는 유용한 실제 물리다. 다만 full LO+SIG
time trace, mixer/lock-in transfer function, sampling, colored technical noise는 없다. 코드와
분석 문서도 이 경계를 명시한다(`gabes/rydberg_electrometry.py:1-35`;
`analysis/rydberg_cell_heating/README.md:96-104`).

### 3.3 heated-cell, axial, SAM, CSV/report 계층

실험 reference 코드로서 특히 좋은 부분이다.

- setpoint/sensor/effective/cold-spot temperature와 fallback source를 보존한다
  (`gabes/rydberg_experiment.py:47-162`).
- cold-spot-limited `n(z)`, column density와 segmented Beer–Lambert를 독립 primitive로 둔다
  (`gabes/rydberg_experiment.py:165-338`).
- Gaussian overlap atom 수에서 geometry, participation과 overlap factor를 분리한다
  (`gabes/rydberg_experiment.py:341-443`).
- SAM은 far-field field, 표준불확도와 `r/(2D^2/lambda)` 경고를 계산한다
  (`gabes/rydberg_experiment.py:446-589`).
- EIT/RF/PSD CSV는 raw detector signal을 임의 보정하지 않고 단위 가정, 무시/중복 row,
  SHA-256 provenance를 보존한다(`gabes/rydberg_experimental_csv.py:1-8`, `211-379`,
  `446-718`).
- workflow는 `MEASURED/PPT/REFERENCE/FITTED/PREDICTED/ASSUMED/PENDING`을 구분하고,
  axial density-rescaled spectrum이 spatial finite-IF OBE가 아님을 밝힌다
  (`analysis/rydberg_cell_heating/README.md:22-39`, `78-117`;
  `analysis/rydberg_cell_heating/workflow.py:1187-1251`).

## 4. 실험 원자물리 reference 적합성 판정

**판정: 장난감 line shape가 아니라 실제 실험 설계와 sanity check에 유용한
semi-quantitative Rydberg-EIT/electrometry reference다. 그러나 absolute RF sensitivity,
cell-heating optimum, 편광·자기장 의존성을 확정하는 계측 reference로는 아직 부적합하다.**

로컬 실험 문맥은 20 °C 부근 EIT linewidth 1.6 MHz, PSN limit 11.2 nV/cm/sqrt(Hz),
측정 감도 12.5(8) nV/cm/sqrt(Hz), 온도 optimum 약 33 °C와 실험 probe optimum 6 µW를
기록한다(`references/rydberg_ju_experiment_context.md:10-22`, `65-71`, `140-159`).
GABES는 level topology, power/beam/IF knob, cold spot, finite-IF와 detector 식을 잘 반영하지만,
다음이 absolute 인용을 막는다.

1. `raw["ls"] = 0.001`이 opacity와 sensitivity를 지배하지만 사용자 parameter도 provenance가
   있는 calibration도 아니다(`gabes/schemes/rydberg.py:625-643`).
2. extra view는 실제로 IF, QE, optical path, reference arm, cell과 RF dipole을 사용하지만 앱은
   recompute knob만 전달한다(`streamlit_app.py:631-644`, `1583-1586`, `1634-1640`).
3. low-contrast/near-zero-power/scan-edge 결과에도 정밀한 linewidth와 sensitivity hero가 남을 수
   있다(`gabes/schemes/rydberg.py:926-1051`, `1121-1154`).
4. reference arm power를 각 detuning의 signal power에서 다시 만들어 실제 fixed arm이 아니라
   이상적인 pointwise-balanced envelope를 계산한다(`gabes/schemes/rydberg.py:849-885`).
5. full Zeeman/polarization, spatial RF field, measured readout transfer와 held-out dataset이 없다.

## 5. 2026-09-02 수치 재검증

현재 working tree, warm JIT 상태에서 측정했다. timing은 이 Windows 환경의 median이므로 절대
benchmark보다 병목과 speedup 비율을 본다.

| 조건 | 현재 결과 |
|---|---|
| EIT default | `T(0)=0.940016615`, contrast `0.135391492`, FWHM `1.613058583 MHz` |
| AT default | `T(0)=0.806029283`, split `3.4965 MHz`, optimum `-2.2200 MHz` |
| AT sensitivity | PSN = total = `6.082387215 nV/cm/sqrt(Hz)` |
| warm timing | EIT compute/readout `8.34/0.77 ms`, AT compute/readout `9.98/169.89 ms` |

### 5.1 숨은 effective opacity

동일 atomic state에서 `raw["ls"]`만 바꿔 optical/detector chain을 다시 계산했다.

| effective `ls` | AT `T(0)` | total sensitivity |
|---:|---:|---:|
| 0.0005 | 0.897791336 | 11.781990 nV/cm/sqrt(Hz) |
| **0.0010** | **0.806029283** | **6.082387 nV/cm/sqrt(Hz)** |
| 0.0020 | 0.649683205 | 3.238776 nV/cm/sqrt(Hz) |

2배 opacity 선택이 감도를 약 1.9배 움직인다. 이 값은 atomic line strength와 분리된 effective
participation/calibration으로 공개되기 전에는 absolute atomic prediction으로 인용하면 안 된다.

### 5.2 extra-view cache dependency 오류

`IF=200 kHz`, QE `0.2`, path efficiency `0.2`, cell `20 mm`, microwave metadata `55 GHz`에서
temperature extra view를 앱과 같은 recompute-only dictionary와 full parameter로 호출했다.

| 경로 | 20 °C sensitivity | 선택된 best T |
|---|---:|---:|
| 앱형 recompute-only 입력 | 6.082387 | 47 °C |
| direct full-parameter 입력 | 52.184629 | 55 °C |

`ExtraView` 자체에는 dependency contract가 없고(`gabes/schemes/base.py:78-84`), 현재 test는
full params를 직접 넘길 뿐 app cache path를 검증하지 않는다(`tests/test_rydberg_eit.py:419-431`).

### 5.3 관측 불가능한 고온 hero

80 °C EIT는 `T(0)=contrast=1.798923e-6`인데도 FWHM `0.307015 MHz`를 내고,
IF optimum은 `-8.96 MHz`로 `-9 MHz` scan edge에 붙는다. 수치상 half-height가 존재하는 것과
실험에서 관측 가능한 것, 중앙 EIT feature를 계속 추적한 것은 서로 다른 조건이다.

### 5.4 grid 좌표 의존성

| scan points | step | AT split | sensitivity optimum | total sensitivity |
|---:|---:|---:|---:|---:|
| 401 | 0.055500 MHz | 3.4410 MHz | -2.220000 MHz | 6.082387 |
| 801 | 0.027750 MHz | 3.4965 MHz | -2.220000 MHz | 6.082387 |
| 1601 | 0.013875 MHz | 3.4965 MHz | -2.206125 MHz | 6.078733 |

기본값은 대체로 안정적이지만 peak와 optimum은 grid point다
(`gabes/schemes/rydberg.py:540-552`, `662-666`, `894-906`, `994-1005`). 보간·경계 status 없이
두 자리 소수로 표시하면 수치 해상도보다 정밀해 보인다.

## 6. 기존 개선안의 계산비용과 물리 보존성

| 기존 항목 | 비용 | 물리 보존성과 실험 의미 |
|---|---:|---|
| ExtraView dependency/caches 분리 | key O(k), 구현 medium | 같은 방정식에 실제 사용자 입력을 전달한다. detector-only 변경은 atomic re-solve를 줄인다. |
| effective opacity 공개 + uncertainty | 기본 O(1), uncertainty O(Nscan) | 현재 default `0.001`에서 curve를 `rtol <= 1e-12`로 보존하면서 empirical participation과 atomic strength를 분리할 수 있다. |
| contrast/power/edge/resolution status | O(Nscan) post-process | spectrum을 바꾸지 않고 과신만 막는다. calibration이 없으면 `visibility not assessed`가 정직하다. |
| fixed reference arm | detector O(Nscan) | atomic response는 그대로고 sensitivity만 더 현실적인 readout 가정으로 의도적으로 바뀐다. |
| AT/optimum local interpolation | 3-point O(1) | grid stair-step만 줄인다. curvature와 edge status를 함께 내야 한다. |
| opt-in 401/801/1601 refinement | 기본 비용 0, 실행 시 O(Ngrid) finite-IF | 기본 UI에는 spacing/status만 싸게 표시하고 reference refinement만 요청 시 실행하면 된다. |
| parameter별 evidence status | O(k) metadata | 수치 불변. published apparatus는 `REFERENCE`, 미측정 throughput/angular/noise는 `ASSUMED`로 나눈다. |
| held-out validation data | 외부 실험 비용이 지배 | absolute reference에 필수이며 fit/holdout를 사전 분리해야 한다. |
| measured readout transfer | complex filter O(Nf) | atomic solve 불변. full time-domain보다 먼저 검증할 저비용 현실화다. |
| RF field-map import | O(Nspace*Nscan), cached reduction 가능 | uniform-map limit을 scalar SAM과 비교하는 opt-in reference가 적절하다. |
| full Zeeman/polarization | Liouville 차원 급증, research | absolute polarization/B dependence에는 중요하지만 negligible-overhead가 아니다. held-out residual이 요구할 때만 활성화한다. |

앞의 일곱 항목은 solve 차원을 늘리지 않거나 O(Nscan) 이하 후처리다. 현재 spectrum physics와
default curve를 보존하면서 신뢰성을 크게 높일 수 있다. full Zeeman, 3-D EM solver와 generic
time-domain lock-in부터 시작하는 것은 현재 증거와 비용에 맞지 않는다.

## 7. 동작을 바꾸지 않는 순수 코드 최적화

### 7.1 finite-IF Liouvillian affine assembly — 최우선

`_superheterodyne_readout()`은 801 detuning마다 Python에서 Hamiltonian과 Liouvillian을 조립한다
(`gabes/schemes/rydberg.py:811-818`). detuning dependence는 정확히 affine이므로
`L(s)=L(0)+s[L(1)-L(0)]` broadcast로 바꿀 수 있다.

- current list assembly median: `88.633 ms`
- affine broadcast median: `3.760 ms`
- **24.40x faster**, `np.array_equal=True`, max absolute difference `0`

방정식, 배열 순서와 부동소수 결과까지 그대로 유지한다.

### 7.2 이미 구한 steady state `rho0` 재사용

finite-IF helper는 optional `steady_state`를 지원하지만(`gabes/rydberg_electrometry.py:116-171`),
Scheme path는 같은 Liouvillian batch의 `rho0`를 helper 안에서 다시 푼다.

- fresh response median: `36.765 ms`
- supplied `rho0` median: `25.358 ms`
- **1.45x faster**, `rho_minus`와 `rho_plus` bit-identical

full density matrix를 장기 global cache하지 않고 한 atomic-phasor 계산 lifetime에서만 재사용하면
mutable cache 위험 없이 동작을 보존한다.

### 7.3 atomic phasor / optical propagation / detector-SAM cache 분리

QE, path, reference balance, RIN, electronics, ENBW, RF dipole와 SAM은 LO-dressed atomic phasor를
바꾸지 않는다(`gabes/schemes/rydberg.py:823-923`). phasor까지만 cache하고 cell/opacity optical
propagation과 detector/SAM을 분리하면 calibration slider가 즉시 반응한다. 이 구조는 extra-view
correctness와 성능을 같은 dependency contract로 해결한다.

## 8. 문서, 테스트와 example 평가

### 잘된 점

- README는 finite-IF, conditional sensitivity, cold spot, SAM과 deferred scope를 현재 코드에
  맞게 설명한다(`README.md:105-127`).
- Scheme tests는 reference knobs, EIT linewidth, AT split/field conversion, low-IF limit,
  detector noise, separated temperatures, SAM과 extra view를 검사한다
  (`tests/test_rydberg_eit.py:29-439`).
- primitive tests는 sideband residual/Hermiticity, peak/RMS, coherent averaging,
  PSN/RIN/electronics, axial integration, dephasing identifiability, SAM uncertainty와 CSV
  provenance를 폭넓게 검사한다(`tests/test_rydberg_electrometry.py:38-203`;
  `tests/test_rydberg_experiment.py:19-173`; `tests/test_rydberg_experimental_csv.py:14-109`).
- workflow tests는 strict status, portable path/hash manifest, cold-spot 전달, scheme/configured
  finite-IF, axial/SAM capability를 검사한다
  (`analysis/rydberg_cell_heating/test_workflow.py:44-299`).
- 실행 가능한 `analysis/rydberg_cell_heating/example_config.json`은 raw EIT/RF/temperature를
  `PENDING`, model output을 `PREDICTED`, RF dipole을 `REFERENCE`로 구분하고 report artifact와
  source manifest를 생성한다(`analysis/rydberg_cell_heating/README.md:8-39`;
  `analysis/rydberg_cell_heating/example_config.json:53-107`).
- adapter는 private Scheme API 의존성을 한 파일에 격리하고 출력에 `PREDICTED` 상태를 붙인다
  (`analysis/rydberg_cell_heating/adapters.py:90-228`).

### 보완할 점

- User Guide는 Rydberg를 여전히 정적 EIT linewidth, AT split, dispersion으로만 설명한다
  (`docs/Userguide/GABES_User_Guide_v2.html:536-539`, `606-616`). finite-IF/noise/cold-spot/SAM과
  conditional-absolute 경계가 없다.
- `docs/Userguide/userguide_assets/rydberg.png`는 2026-06-08의 transmission/dispersion 2-panel
  그림이다. 현재 AT output은 transmission/finite-IF response/dispersion 3-panel이다
  (`gabes/schemes/rydberg.py:1193-1223`).
- test module docstring의 “public UI shows the static spectrum only”는 현재 공개 sensitivity와
  모순된다(`tests/test_rydberg_eit.py:1-5`).
- tests는 hidden opacity와 uncertainty, 33 °C/6 µW external optimum, 80 °C hero 억제,
  fixed reference arm, sub-grid convergence와 app-path extra-view dependency를 검증하지 않는다.
  오히려 현재 AT sensitivity hero를 회귀값으로 고정한다(`tests/test_rydberg_eit.py:99-125`).
- example은 좋은 pipeline demonstration이지만 held-out experiment-model agreement 예제는 아니다.
  `electrometry.status=ASSUMED` 하나가 published apparatus와 실제 미측정 calibration을 함께 묶는다
  (`analysis/rydberg_cell_heating/example_config.json:98-102`). 로컬 실험 문맥은 같은 apparatus의
  공개값을 `REFERENCE`로 분류해야 한다고 명시한다
  (`references/rydberg_ju_experiment_context.md:26-46`).

## 9. 검증 상태

- Scheme 3, finite-IF, heated-cell, CSV, workflow, kernel/headless/render/docs 관련:
  **102 passed in 51.89 s**
- 프로젝트 본체 범위 `python -m pytest -q tests analysis`:
  **509 passed in 465.03 s**
- 정책 명령 `python -m pytest -q`는 production test failure가 아니라 collection 단계에서
  `references/AutoOD/test_autood_gui_adjustments.py`와
  `references/AutoOD-NatRbD2/test_gui_adjustments.py`가 optional `sympy`를 import하지 못해 중단됐다.
  `requirements.txt:1-6`에는 `sympy`가 없다. 본체 test path를 명시하거나 reference-test extra를
  별도로 선언하기 전에는 root-wide collection을 재현 가능한 green gate로 볼 수 없다.
- production code는 수정하지 않았다. 기존 dirty working tree를 보존하고 이 보고서만 추가했다.

## 10. 결론과 권장 순서

Scheme 3은 compact cascade OBE, residual Doppler, power/beam/transit scaling, finite-IF complex
response, balanced-detector noise, cold-spot/axial density, SAM과 provenance-aware CSV를 갖춘
유용한 실험 설계 도구다. 단위/convention/evidence status를 primitive 수준에서 분리한 점은
실험물리 reference 코드로 좋은 방향이다.

그러나 현재 `6.082 nV/cm/sqrt(Hz)`, 47 °C optimum과 내부 power optimum은 **조건부 모델 출력**이다.
숨은 opacity, extra-view cache dependency 오류, pointwise-balanced reference arm,
observability/convergence guardrail과 held-out data 부재 때문에 absolute 결과로 인용할 수 없다.

권장 순서는 다음과 같다.

1. **P1 correctness + speed:** ExtraView dependency contract와 app-path regression을 추가하고,
   atomic state / finite-IF phasor / optical / detector-SAM cache를 분리한다.
2. **P1 trust:** `ls=0.001`을 provenance와 uncertainty가 있는 effective-opacity participation으로
   공개하고 temperature/power sweep을 `PREDICTED / EXTERNALLY UNVALIDATED`로 표시한다.
3. **P1 readout:** fixed reference arm, contrast/detected-power/edge/samples-per-width status,
   AT/optimum interpolation과 opt-in refinement gate를 함께 구현한다.
4. **P1 docs/tests:** parameter별 `REFERENCE`/`ASSUMED`를 나누고 User Guide, static image와 stale
   test docstring을 현재 finite-IF 기능에 맞춘다. app-path cache와 high-temperature failure fixture를
   추가한다.
5. **순수 최적화:** bit-identical affine Liouvillian assembly와 local `rho0` 재사용을 적용한다.
6. 그 뒤 held-out EIT/RF/PSD/temperature와 detector/SAM calibration으로 compact model을 검증하고,
   잔차가 요구할 때만 Zeeman/polarization, measured transfer와 RF field-map reference를 활성화한다.

1–5는 solve 차원을 늘리지 않거나 O(Nscan) 이하 후처리이며, 현재 spectrum physics와 default
curve를 보존하면서 correctness, 신뢰성과 속도를 함께 높일 수 있다.
