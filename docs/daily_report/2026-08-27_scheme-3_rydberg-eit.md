# 2026-08-27 Scheme 3 물리 검토 — Rydberg-EIT electrometry

## 1. 선택 규칙과 현재 다섯 스킴

- 현지 날짜: 2026-08-27 (Asia/Seoul), 일자 `27`
- `n = (27 mod 5) + 1 = 3`
- 오늘의 대상: **Scheme 3 — Rydberg-EIT electrometry**

드롭다운 순서는 `gabes/schemes/__init__.py:19-24`의 `_SCHEMES`가 결정한다.
README의 현재 기능 표(`README.md:8-16`)와도 순서와 정의가 일치한다.

| 번호 | 등록 클래스 / registry name | 현재 정의 |
|---:|---|---|
| 1 | `SASScheme()` / `sas` | Absorption spectroscopy — OD / SAS |
| 2 | `LambdaScheme()` / `lambda` | Λ coherence — EIT / AT / CPT |
| **3** | `RydbergEITScheme()` / `rydberg_eit` | **⁸⁵Rb cascade EIT / microwave AT electrometry** |
| 4 | `MagnetoScheme()` / `magneto` | Hanle / EIA / NMOR |
| 5 | `FWMScheme()` / `fwm` | seeded FWM gain diagnostic / generic SFWM biphoton |

Scheme 3의 공개 이름·제목·cache/default version은
`gabes/schemes/rydberg.py:91-100`에 정의돼 있다.

## 2. 먼저 확인한 기존 제안과 변경 여부

새 개선안을 처음부터 만들 필요는 없었다. Scheme 3 선행 일일 보고서는 2026-06-22부터
2026-08-22까지 9개이며, 07-22, 08-02, 08-07, 08-12, 08-22 보고서가 현재 문제를
재현 가능한 제안으로 정리했다. 로컬 파일명 검색과 Scheme 3 경로의 inline 검색에서는 별도
`TODO`/`ISSUE`/proposal 파일 및 `TODO`/`FIXME`가 없었다. 활성 제안은
`docs/checklist.json`에 통합돼 있다. 연결된 공개 GitHub 저장소
`Shake2313/fwm-squeezing-app`도 2026-08-27 현재 열린 이슈 0건, 닫힌 이슈 0건이다.

현재 핵심 항목은 다음과 같다.

- **P1 ready — extra-view dependency correctness:** extra view가 IF, QE, optical path,
  reference arm, cell length, RF dipole을 쓰지만 앱 cache key에는 recompute knob만 들어간다
  (`docs/checklist.json:464-480`; `streamlit_app.py:609-622`, `1555-1559`, `1606-1612`).
- **P1 ready — effective opacity / sweep trust:** 고정된 `ls=0.001`을 provenance와 불확도를
  갖는 effective participating opacity로 공개하고, 온도·파워 sweep을 검증 전
  `PREDICTED`로 제한한다 (`docs/checklist.json:484-503`).
- **P1 needs_decision — readout validity:** 저투과·저대비·scan-edge·저해상도에서
  linewidth/AT/sensitivity hero를 억제하고, fixed reference arm과 local refinement를 요구한다
  (`docs/checklist.json:507-528`).
- **P2 blocked_external:** held-out EIT/RF/PSD/temperature dataset, measured readout transfer
  function, provenance-tagged RF field map은 외부 측정이 필요하다
  (`docs/checklist.json:532-552`, `883-931`).
- **P2 parked:** full Zeeman/polarization reference는 상태쌍·편광·B-field·검증 spectrum·runtime
  budget이 정해질 때까지 보류돼 있다 (`docs/checklist.json:858-880`).

2026-08-22 검토 뒤 새 커밋은 없다. 현재 working tree의 기본 결과, cache reproducer,
고온 관측성 실패, grid 의존성도 당시와 동일했다. 따라서 기존 항목을 완료로 바꾸거나 새로운
대규모 물리 확장을 제안할 근거는 없다.

## 3. 현재 구현이 담은 실제 물리

### 3.1 축약 4준위 cascade OBE

코어는 ⁸⁵Rb
`5S1/2 F=3 → 5P3/2 F'=4 → 40D5/2 → 39F7/2`의 4준위 정상상태 OBE다.
780 nm probe, 481 nm counter-propagating coupling, 약 37 GHz RF LO를 Hamiltonian에 넣고,
5P spontaneous decay, 작은 Rydberg decay, 5S–40D와 40D–39F dephasing을 분리한다
(`gabes/schemes/rydberg.py:1-15`, `69-88`, `413-433`, `554-623`).

- probe/coupling Rabi는 `Ω ∝ √P/d`로 기준 실험점에 anchor된다
  (`gabes/schemes/rydberg.py:39-53`, `435-457`).
- heater setpoint, effective vapor temperature, cold spot을 분리한다. effective temperature는
  thermal motion/transit에, cold spot은 sealed-cell pressure에 들어간다
  (`gabes/schemes/rydberg.py:459-501`, `559-580`).
- residual two-photon Doppler는 `(k_probe-k_coupling)/k_probe`를 Rydberg level에 적용하고,
  Doppler-on에서 Maxwell velocity class를 평균한다
  (`gabes/schemes/rydberg.py:69-88`, `503-523`, `606-615`).
- residual Zeeman은 full magnetic manifold가 아니라 compensated susceptibility의 Gaussian
  convolution이다 (`gabes/schemes/rydberg.py:617-623`).
- coherence를 D2 dipole, density, effective line factor와 Beer–Lambert 법칙으로 transmission에
  바꾼다 (`gabes/schemes/rydberg.py:625-677`).

따라서 기준점 주변 EIT 폭, power/waist/transit trade-off, RF AT splitting과 microwave detuning
방향을 보는 데 실제로 유용하다. 반면 probe/coupling Rabi와 transit factor는 bare-state
first-principles 예측이 아니라 같은 논문 조건에 맞춘 anchor다. 기준 선폭을 다시 얻는 것만으로는
독립 물리 검증이 아니다.

### 3.2 finite-IF atomic superheterodyne와 detector chain

정적 AT 자나 slope proxy에 그치지 않고 frequency-domain small-signal chain이 구현돼 있다.

- LO-dressed 정상상태 주위에서 `L₀ ± iω_IF` sideband를 풀고 SIG phase 및 peak angular-Rabi
  convention을 명시한다 (`gabes/rydberg_electrometry.py:116-199`).
- Hermitian absorption quadrature와 velocity-class complex phasor의 coherent average를 쓴다
  (`gabes/rydberg_electrometry.py:70-113`, `216-232`;
  `gabes/schemes/rydberg.py:748-821`).
- RF dipole, angular/polarization factor, peak/RMS field convention을 분리한다
  (`gabes/rydberg_electrometry.py:235-321`).
- 두 photodiode의 one-sided Schottky ASD, correlated RIN, electronics ASD를 합쳐 field ASD와
  ENBW noise로 환산한다 (`gabes/rydberg_electrometry.py:395-530`).
- 기본 Doppler-off path는 801개 probe detuning에서 total noise-equivalent field를 최소화한다
  (`gabes/schemes/rydberg.py:769-923`).

이는 atomic mixer의 선형·유한 IF 한계로서 유용한 실제 물리다. 다만 full LO+SIG time trace,
mixer/lock-in transfer function, sampling, colored technical noise는 없다. 코드와 분석 문서도 이
경계를 명시한다 (`gabes/rydberg_electrometry.py:1-7`;
`analysis/rydberg_cell_heating/README.md:96-104`).

### 3.3 heated-cell, axial, SAM, CSV/report 계층

실험 reference로서 가장 잘 설계된 부분이다.

- setpoint/sensor/effective/cold-spot temperature와 source를 보존한다
  (`gabes/rydberg_experiment.py:48-162`).
- cold-spot-limited `n(z)`, column density, segmented Beer–Lambert를 별도 primitive로 둔다
  (`gabes/rydberg_experiment.py:166-338`).
- Gaussian overlap atom 수에서 geometry, participation, overlap factor를 분리한다
  (`gabes/rydberg_experiment.py:341-443`).
- SAM은 far-field field, 표준불확도와 `r/(2D²/λ)` 경고를 계산한다
  (`gabes/rydberg_experiment.py:447-589`).
- EIT/RF/PSD CSV는 단위 가정, 무시/중복 row, SHA-256 provenance를 보존한다
  (`gabes/rydberg_experimental_csv.py:211-379`, `446-718`).
- workflow는 `MEASURED/PPT/REFERENCE/FITTED/PREDICTED/ASSUMED/PENDING`을 구분하고,
  axial density-rescaled spectrum이 spatial finite-IF OBE가 아님을 명시한다
  (`analysis/rydberg_cell_heating/README.md:22-39`, `78-117`;
  `analysis/rydberg_cell_heating/workflow.py:1037-1265`).

## 4. 실험물리 reference 적합성 판정

**판정: 장난감 line shape가 아니라 실제 실험 설계와 sanity check에 쓸 수 있는
semi-quantitative Rydberg-EIT/electrometry reference다. 그러나 absolute RF sensitivity,
cell-heating optimum, 편광·자기장 의존성을 확정하는 계측 reference로는 아직 부적합하다.**

로컬 실험 문맥은 20 °C 부근 EIT linewidth 1.6 MHz, 계산 PSN limit
11.2 nV/cm/√Hz, 측정 감도 12.5(8) nV/cm/√Hz, 온도 optimum 약 33 °C를 기록한다
(`references/rydberg_ju_experiment_context.md:10-22`, `65-71`, `140-159`). GABES는 상태
topology, power/beam/IF knob, cold spot, finite-IF/detector 식을 잘 반영한다. 하지만 다음이
absolute 인용을 막는다.

1. `raw["ls"]=0.001`이 opacity와 sensitivity를 지배하지만 사용자 parameter/provenance가 아니다
   (`gabes/schemes/rydberg.py:625-643`).
2. extra view가 일부 실제 detector/IF/cell 입력을 받지 못한다
   (`streamlit_app.py:609-622`, `1555-1559`, `1606-1612`).
3. low-contrast/edge 결과와 Doppler-on 단일점 최적화도 정밀한 hero가 될 수 있다
   (`gabes/schemes/rydberg.py:772-821`, `926-1051`, `1137-1154`).
4. reference arm은 detuning마다 `signal power × ratio`로 다시 만들어져 실제 고정 arm이 아니라
   이상적인 pointwise balance envelope다 (`gabes/schemes/rydberg.py:849-885`).
5. full Zeeman/polarization, spatial RF field, measured readout transfer 및 held-out dataset이 없다.

## 5. 오늘의 수치 재검증

현재 working tree, warm JIT 상태에서 측정했다. timing은 이 Windows 환경의 median이며 절대
benchmark보다 병목 비율을 본다.

| 조건 | 현재 결과 |
|---|---|
| EIT default | `T(0)=0.940016615`, contrast `0.135391492`, FWHM `1.613058583 MHz` |
| AT default | `T(0)=0.806029283`, split `3.4965 MHz`, optimum `-2.2200 MHz` |
| AT sensitivity | PSN = total = `6.082387215 nV/cm/√Hz` |
| warm timing | EIT compute/readout `3.84/0.23 ms`, AT compute/readout `3.64/69.70 ms` |

### 5.1 숨은 effective opacity

동일 atomic state에서 `raw["ls"]`만 바꿔 optical/detector chain을 다시 계산했다.

| effective `ls` | AT `T(0)` | total sensitivity |
|---:|---:|---:|
| 0.0005 | 0.897791336 | 11.781990 nV/cm/√Hz |
| **0.0010** | **0.806029283** | **6.082387 nV/cm/√Hz** |
| 0.0020 | 0.649683205 | 3.238776 nV/cm/√Hz |

2배 opacity 선택이 감도를 약 1.9배 움직인다. 이 값은 atomic line strength와 분리된 effective
participation/calibration으로 공개되기 전에는 absolute atomic prediction으로 인용하면 안 된다.

### 5.2 외부 optimum과 불일치

- 내장 10점 temperature extra view: best **47 °C**, `0.971731 nV/cm/√Hz`.
- 선행 검토의 0.5 µW 간격 probe-power sensitivity sweep: best **4.0 µW**, `5.911989`;
  6 µW에서는 `6.082387`.
- 로컬 reference context: 온도 optimum 약 **33 °C**, 실험 probe optimum **6 µW**
  (`references/rydberg_ju_experiment_context.md:20-22`, `35-43`, `65-71`).

논문 값을 강제로 맞춰서는 안 된다. 현재 sweep을 `PREDICTED/EXTERNALLY UNVALIDATED`로 표시하고
opacity, dephasing, detector/calibration uncertainty를 함께 탐색해야 한다.

### 5.3 저투과·scan-edge hero

80 °C EIT는 `T(0)=contrast=1.798923×10⁻⁶`인데 FWHM `0.307015 MHz`와 transmission이
계속 hero다. IF optimum은 `-8.96 MHz`, scan edge는 `-9 MHz`다. 수치적 half-height가 존재하는
것, 실험에서 관측 가능한 것, 중앙 EIT를 추적한 것은 서로 다른 판정이다.

### 5.4 extra-view cache 불일치

`IF=200 kHz`, `QE=0.2`, path efficiency `0.2`, cell `20 mm`, microwave metadata `55 GHz`에서
temperature extra view를 앱형 recompute-only dictionary와 전체 parameter로 호출했다.

| 경로 | 20 °C sensitivity | 선택된 best T |
|---|---:|---:|
| 앱과 같은 recompute-only dictionary | `6.082387` | 47 °C |
| direct full-parameter | `52.184629` | 55 °C |

`ExtraView`에는 dependency contract가 없고(`gabes/schemes/base.py:78-84`), 앱은
recompute-only dictionary를 전달한다. 현재 extra-view test는 full params를 직접 넘기므로 UI 경로를
검증하지 못한다 (`tests/test_rydberg_eit.py:419-431`).

### 5.5 grid 좌표 의존성

| scan points | Δx | AT split | sensitivity optimum | total sensitivity |
|---:|---:|---:|---:|---:|
| 401 | 0.055500 MHz | 3.4410 MHz | -2.220000 MHz | 6.082387 |
| 801 | 0.027750 MHz | 3.4965 MHz | -2.220000 MHz | 6.082387 |
| 1601 | 0.013875 MHz | 3.4965 MHz | -2.206125 MHz | 6.078733 |

기본 결과는 대체로 안정적이지만 AT peak와 optimum은 grid point다
(`gabes/schemes/rydberg.py:540-552`, `662-666`, `894-906`, `994-1005`). 보간 없이 split과
optimum을 0.01 MHz로 표시하는 것은 수치 해상도보다 정밀하다.

## 6. 기존 개선안의 계산비용과 물리 보존성

| 기존 항목 | 상태 | 추가 비용 | 물리 보존 / 실험 의미 |
|---|---|---:|---|
| ExtraView dependency와 cache 계층 분리 | P1 `ready`, medium | key는 O(k); detector-only 변경은 atomic re-solve를 제거 | 같은 방정식에 실제 사용자 입력을 전달한다. steady state / finite-IF phasor / Beer–Lambert / detector-SAM 경계를 분리해야 한다. |
| effective opacity 공개 + uncertainty | P1 `ready`, medium | 기본 표시는 O(1), uncertainty sweep은 O(Nscan); OBE 차원 증가 없음 | 현재 default를 그대로 보존하면서 empirical participation과 atomic strength를 구분한다. |
| contrast/power/edge/resolution status | P1 `needs_decision`, medium | O(Nscan) 후처리, solve 없음 | 이미 계산한 배열을 진단하므로 spectrum physics는 그대로다. calibration이 없으면 `visibility not assessed`가 정직하다. |
| fixed reference arm | readout P1에 통합 | detector O(Nscan), atomic solve 없음 | atomic response는 보존하고 실제 한 arm의 DC power와 명시된 balance detuning을 복원한다. sensitivity는 의도적으로 더 현실적으로 바뀐다. |
| AT/optimum interpolation | readout P1에 통합 | 3점 local fit O(1) | grid stair-step만 줄인다. curvature와 edge status를 함께 내야 한다. |
| deterministic 401/801/1601 reference | opt-in 권장 | 기본 path 0; 요청 시 O(Ngrid) finite-IF | 기본 UI는 spacing/status만 싸게 표시하고 reference refinement만 opt-in으로 두면 된다. |
| parameter별 evidence status | 기존 issue note | O(k) metadata, 수치 변화 없음 | published apparatus 값은 `REFERENCE`, path/angular/noise/linked-temperature 가정은 `ASSUMED`로 분리한다. |
| held-out EIT/RF/PSD/temperature dataset | P2 `blocked_external`, large | 계산보다 측정·provenance 비용이 큼 | absolute reference에 필수이며 fit/holdout를 사전 분리해야 한다. |
| measured readout transfer function | P2 `blocked_external`, medium | complex filter O(Nf), atomic solve 없음 | full time-domain보다 훨씬 싸며 측정된 선형 chain을 우선 검증한다. |
| RF field-map import | P2 `blocked_external`, large | O(Nspace·Nscan) 또는 cached reduction | uniform map이 scalar SAM으로 환원되는 별도 reference path가 적절하다. |
| full Zeeman/polarization | P2 `parked`, research | Hilbert/Liouville 차원 급증 | absolute polarization/B 의존성에는 중요하지만 negligible-overhead가 아니다. held-out residual이 요구할 때만 활성화한다. |

앞의 여섯 항목과 evidence-status 수정은 solve 차원을 늘리지 않으며, 잘 구현하면 현재 default curve를
`≤1e-12` 상대오차로 유지하거나 readout 가정만 더 현실화할 수 있다. full Zeeman, 3-D field,
time-domain부터 시작하는 것은 현재 증거와 비용에 맞지 않는다.

## 7. 동작을 바꾸지 않는 순수 코드 최적화

### 7.1 finite-IF Liouvillian affine assembly — 최우선

현재 `_superheterodyne_readout()`은 801 detuning마다 Python에서 Hamiltonian과 Liouvillian을
조립한다 (`gabes/schemes/rydberg.py:811-818`). detuning dependence는 정확히 affine이므로
`L(s)=L(0)+s[L(1)-L(0)]` broadcast로 바꿀 수 있다.

- current list assembly median: `34.519 ms`
- affine broadcast median: `1.819 ms`
- **18.98× faster**, `np.array_equal=True`, max absolute difference `0`

방정식, 배열 순서, 부동소수 결과를 그대로 유지하는 가장 확실한 최적화다.

### 7.2 이미 구한 정상상태 `rho₀` 재사용

finite-IF helper는 `steady_state` optional 인자를 지원하지만
(`gabes/rydberg_electrometry.py:116-171`), Scheme path는 동일 Liouvillian batch의 `rho₀`를 다시 푼다.

- fresh response median: `13.811 ms`
- supplied `rho₀` median: `9.429 ms`
- **1.46× faster**, `rho_minus/rho_plus` bit-identical

full density matrix를 장기 global cache하지 말고 atomic-phasor 계산 lifetime에만 보존하면 mutable
cache와 메모리 위험을 피할 수 있다.

### 7.3 atomic phasor / optical propagation / detector-SAM cache 분리

QE, path, reference balance, RIN, electronics, ENBW, RF dipole, SAM은 LO-dressed atomic phasor를
바꾸지 않는다 (`gabes/schemes/rydberg.py:837-923`). phasor까지만 cache하고 cell/opacity optical
propagation과 detector/SAM을 별도로 재계산하면 calibration slider가 즉시 반응하며 계산 동작은
동일하다. 이는 성능과 extra-view correctness를 한 dependency contract로 함께 해결한다.

## 8. 문서·테스트·예제 평가

### 잘된 점

- README는 finite-IF, conditional sensitivity, cold spot, SAM, deferred scope를 현재 코드와
  대체로 일치하게 설명한다 (`README.md:105-127`).
- Scheme test는 reference knobs, EIT linewidth, AT split/field conversion, low-IF limit,
  finite-IF optimum, detector noise, separated temperatures, SAM, extra view를 검사한다
  (`tests/test_rydberg_eit.py:29-439`).
- primitive tests는 sideband residual/Hermiticity, peak/RMS, coherent average, PSN/RIN/electronics,
  axial integration, dephasing identifiability, SAM uncertainty, CSV units/provenance를 폭넓게 검증한다
  (`tests/test_rydberg_electrometry.py:38-203`; `tests/test_rydberg_experiment.py:19-173`;
  `tests/test_rydberg_experimental_csv.py:14-109`).
- cell-heating workflow test는 status vocabulary, artifact emission, cold spot 전달, native finite-IF,
  axial/SAM capability를 검증한다 (`analysis/rydberg_cell_heating/test_workflow.py:43-246`).
- example config는 raw EIT/RF/temperature input을 `PENDING`, model을 `PREDICTED`, ARC dipole을
  `REFERENCE`로 두어 raw data 없이 validation 완료처럼 보이지 않게 한다
  (`analysis/rydberg_cell_heating/example_config.json:53-107`).
- adapter는 report가 private scheme API를 사용한다는 사실을 격리·기록하고, 모든 출력에
  `PREDICTED` 상태를 붙인다 (`analysis/rydberg_cell_heating/adapters.py:90-228`, `231-354`).

### 보완할 점

- 사용자 가이드는 Rydberg를 여전히 정적 spectrum, EIT linewidth, AT split, dispersion으로만
  설명하며 finite-IF/noise/cold-spot/SAM을 누락한다
  (`docs/Userguide/GABES_User_Guide_v2.html:536-539`, `606-616`). 2026-06-08 생성
  `docs/Userguide/userguide_assets/rydberg.png`도 현재 AT의 transmission/finite-IF response/dispersion
  3패널이 아니라 예전 transmission/dispersion 2패널이다.
- test module docstring의 “public UI shows the static spectrum only”도 현재 공개 sensitivity와
  모순된다 (`tests/test_rydberg_eit.py:1-5`).
- 가이드는 numerics에서 scan resolution/velocity step을 확인하라고 하지만
  (`docs/Userguide/GABES_User_Guide_v2.html:789-798`), Rydberg UI에는 Doppler on/off만 있고 두 grid는
  고정이다. 임의 knob 증설보다 deterministic convergence/status panel이 안전하다.
- tests는 hidden opacity와 uncertainty, 33 °C/6 µW external optimum, high-temperature hero 억제,
  fixed reference detector, sub-grid convergence, app-path extra-view dependency를 검증하지 않는다.
  오히려 기본 AT sensitivity hero를 현재 동작으로 고정한다 (`tests/test_rydberg_eit.py:99-125`).
- example은 pipeline demonstration으로는 좋지만 held-out experiment-model agreement 예제는 아니다.
  published apparatus parameter 전체를 `electrometry.status=ASSUMED`로 묶은 부분도
  (`analysis/rydberg_cell_heating/example_config.json:98-102`) parameter별 evidence status로 나눠야 한다.

## 9. 검증

- Scheme 3 + finite-IF + heated-cell + CSV + workflow 관련:
  `60 passed in 16.25 s`
- 저장소 정책의 전체 `python -m pytest -q`:
  `476 passed in 241.64 s`
- production code는 수정하지 않았다. 기존 dirty working tree를 보존하고 오늘 보고서만 추가했다.

## 10. 결론과 권장 순서

Scheme 3은 compact cascade OBE, residual Doppler, power/beam/transit scaling, finite-IF complex
response, balanced-detector noise, cold-spot/axial density, SAM 및 provenance-aware CSV를 갖춘 유용한
실험 설계 도구다. 특히 단위/convention과 evidence status를 primitive 수준에서 분리한 점은 연구실
reference 코드로 좋은 방향이다.

그러나 현재 `6.082 nV/cm/√Hz`, 47 °C optimum, 4 µW optimum은 **조건부 모델 출력**이다.
숨은 opacity, extra-view cache dependency 오류, pointwise-balanced reference arm,
observability/convergence guardrail 부재와 held-out data 부재 때문에 absolute 결과로 인용할 수 없다.

권장 순서는 다음과 같다.

1. **P1:** ExtraView별 dependency contract와 app-path 회귀 test를 추가하고 atomic state /
   finite-IF phasor / optical / detector-SAM cache를 분리한다.
2. **P1:** `ls=0.001`을 provenance와 uncertainty가 있는 effective opacity/participation으로 공개하고
   temperature/probe-power sweep을 `PREDICTED/EXTERNALLY UNVALIDATED`로 표시한다.
3. **P1:** fixed reference arm, contrast/detected-power/edge/samples-per-width status,
   AT/optimum interpolation과 opt-in refinement gate를 함께 구현한다.
4. **P1 문서:** published apparatus parameter는 `REFERENCE`, 실제 미측정 calibration은
   `ASSUMED`로 parameter별 분리하고 사용자 가이드·정적 이미지·test docstring을 현재 기능에 맞춘다.
5. **속도:** bit-identical affine Liouvillian assembly와 `rho₀` 재사용을 적용한다.
6. 그 뒤 held-out EIT/RF/PSD/temperature와 detector/SAM calibration으로 compact model을 검증하고,
   잔차가 요구할 때만 full Zeeman, measured transfer function, RF field-map reference를 활성화한다.

1–5는 solve 차원을 늘리지 않거나 O(Nscan) 이하 후처리이며, 현재 spectrum physics와 default curve를
보존하면서 신뢰성과 속도를 함께 높일 수 있다.
