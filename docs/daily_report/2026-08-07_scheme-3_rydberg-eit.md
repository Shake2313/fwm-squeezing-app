# 2026-08-07 Scheme 3 물리 검토: Rydberg-EIT electrometry

## 1. 오늘의 선택과 현재 5개 scheme 순서

- 로컬 날짜는 2026-08-07이고, `n = (7 mod 5) + 1 = 3`이다.
- 드롭다운 순서는 `gabes/schemes/__init__.py:19-25`의 `_SCHEMES`가 결정한다.
- 따라서 오늘의 대상은 **Scheme 3, `rydberg_eit` — Rydberg-EIT electrometry**이다
  (`gabes/schemes/rydberg.py:91-100`).

| 번호 | registry 객체 | 현재 표시 이름 | 코드 근거 |
|---:|---|---|---|
| 1 | `SASScheme()` | Absorption spectroscopy (OD / SAS) | `gabes/schemes/sas.py:53-56` |
| 2 | `LambdaScheme()` | Λ coherence (EIT / AT / CPT) | `gabes/schemes/absorption.py:492-500` |
| 3 | `RydbergEITScheme()` | Rydberg-EIT electrometry | `gabes/schemes/rydberg.py:91-100` |
| 4 | `MagnetoScheme()` | Magneto-optics (Hanle/MOR) | `gabes/schemes/magneto.py:157-169` |
| 5 | `FWMScheme()` | Four-wave mixing (Squeezing / Biphoton) | `gabes/schemes/fwm.py:1816-1819` |

평가 대상은 현재 `main`의 `e1e3cdc` 위 working tree다. Scheme 3 핵심 소스와 테스트는
2026-07-23 이후 수정되지 않았고, 2026-08-02 Scheme 3 검토 뒤의 커밋은 SABES와 공통
체크리스트 작업이다. 기존 사용자 변경은 건드리지 않고 이 보고서만 추가했다.

## 2. 먼저 찾은 기존 제안, TODO, daily report, issue note

기존 개선안은 **충분히 존재한다**.

- `docs/checklist.json:8-40`에는 power-to-Rabi, finite-IF electrometry, heated-cell,
  experimental report, axial/SAM 계층이 `done`으로 기록돼 있다.
- `docs/checklist.json:116-120`의 저차 polarization/Zeeman proxy는 `deferred`다.
- `docs/checklist.json:158-183`의 full Zeeman/polarization/ionization, time-domain
  superheterodyne, 3-D RF field, three-photon imaging은 명시적인 GROUP C다.
- Scheme 3의 이전 보고서 `2026-06-22`, `06-27`, `07-02`, `07-12`, `07-22`,
  `08-02`를 대조했다. 특히 `2026-07-22_scheme-3_rydberg-eit.md:107-176`과
  `2026-08-02_scheme-3_rydberg-eit.md:132-227`의 숨은 opacity, 저대비 hero,
  grid 정밀도, ideal reference-arm 제안은 아직 유효하다.
- Rydberg 경로에는 별도 `TODO`/`FIXME` 또는 issue-note 파일이 없다. 공개
  [GitHub Issues](https://github.com/Shake2313/fwm-squeezing-app/issues?q=is%3Aissue)도
  2026-08-07 확인 시 open/closed 합계 0건이다.

따라서 새 기능 목록을 임의로 늘리기보다 기존 저비용 신뢰성 제안의 현재 비용을 다시
평가하고, 이번에는 기준 논문의 **probe-power 및 온도 최적점**과 직접 대조했다.

## 3. 현재 구현이 담는 실제 물리

### 3.1 정적 4준위 cascade와 광학 전파

축약 원자계는 ⁸⁵Rb
`5S1/2 F=3 → 5P3/2 F'=4 → 40D5/2 → 39F7/2`이다. 5P spontaneous decay,
Rydberg decay, 5S–40D 및 40D–39F coherence dephasing, 780 nm probe, 481 nm
counter-propagating coupling, 37 GHz RF LO를 steady-state OBE에 넣는다
(`gabes/schemes/rydberg.py:4-15`, `69-88`, `413-433`, `554-623`).

- probe/coupling Rabi는 기준점에 anchor된 `Ω ∝ √P/d`를 따른다
  (`gabes/schemes/rydberg.py:435-457`).
- effective vapor temperature와 공통 beam diameter가 transit broadening을 정한다
  (`gabes/schemes/rydberg.py:493-501`).
- `doppler=on`은 residual `(k_probe-k_coupling)v`를 Maxwell 평균한다
  (`gabes/schemes/rydberg.py:503-523`, `606-615`).
- optical coherence는 physical susceptibility와 Beer–Lambert transmission으로 연결된다
  (`gabes/schemes/rydberg.py:625-643`, `668-677`; `gabes/observables.py:393-418`).
- 중앙 EIT FWHM은 half-height crossing을 선형 보간하지만, AT peak와 sensitivity optimum은
  여전히 scan grid 위치를 사용한다 (`gabes/schemes/rydberg.py:926-950`, `994-1005`,
  `894-906`).

기준 논문은 같은 4준위 cascade, 6 µW probe, 30 mW coupling, 0.15 mm beam,
50 mm cell, 3.0 MHz coupling Rabi, 3.7 MHz LO Rabi, 40 kHz IF를 사용하고 1.6 MHz
EIT linewidth를 보고한다. 이 실험 조건은 코드와 테스트에 정확히 고정돼 있다
(`gabes/schemes/rydberg.py:39-53`, `106-136`; `tests/test_rydberg_eit.py:29-61`;
[Ju et al., arXiv:2606.04354](https://arxiv.org/pdf/2606.04354)).

### 3.2 finite-IF weak-SIG와 detector chain

- LO-dressed Liouvillian 주위에서 `L₀ ± iω_IF` 두 sideband의 1차 응답을 푼다
  (`gabes/rydberg_electrometry.py:116-199`).
- Doppler class의 복소 phasor를 coherent average한다
  (`gabes/rydberg_electrometry.py:69-113`, `216-232`).
- RF dipole, angular/polarization factor, peak/RMS convention을 SI field와 Rabi 사이에
  명시한다 (`gabes/rydberg_electrometry.py:235-321`).
- signal/reference photodiode의 one-sided Schottky ASD, correlated RIN, electronics ASD를
  합친다 (`gabes/rydberg_electrometry.py:324-488`).
- probe detuning마다 noise-equivalent field를 계산해 최소점을 고른다
  (`gabes/schemes/rydberg.py:769-923`).

sideband 방정식, Hermiticity, low-IF static derivative limit, high-IF roll-off,
peak/RMS round trip, balanced shot/RIN/electronics noise는 전용 테스트로 검증된다
(`tests/test_rydberg_electrometry.py:38-203`; `tests/test_rydberg_eit.py:265-353`).
즉 수학적 frequency-domain linear response는 실제 유용한 물리다. 다만 논문의 time trace,
lock-in phase/filter, sampling, colored noise를 직접 재현하는 엔진은 아니다
(`gabes/rydberg_electrometry.py:1-35`; `docs/checklist.json:165-169`).

### 3.3 온도·cold spot·SAM·실험 데이터 계층

- heater setpoint, sensor/effective vapor temperature, cold spot을 분리한다
  (`gabes/rydberg_experiment.py:48-162`).
- cold spot vapor pressure로 illuminated hot region의 local density를 계산한다
  (`gabes/schemes/rydberg.py:459-491`).
- axial `T(z),n(z)`, segmented Beer–Lambert, Gaussian-overlap `N_eff` primitive가 있다
  (`gabes/rydberg_experiment.py:165-443`).
- SAM은 far-field scalar field, 표준불확도, `r/(2D²/λ)` 유효성 경고를 제공한다
  (`gabes/rydberg_experiment.py:447-589`).
- EIT/RF/PSD CSV loader는 단위, duplicate, SHA-256 provenance를 보존한다
  (`gabes/rydberg_experimental_csv.py:340-373`, `446-718`).

`analysis/rydberg_cell_heating`은 `MEASURED/PPT/REFERENCE/FITTED/PREDICTED/ASSUMED/PENDING`
상태를 분리하고, axial path가 local OBE 재해가 아닌 density-rescaled diagnostic임을 명시한다
(`analysis/rydberg_cell_heating/README.md:22-39`, `78-117`). 이는 실험 notebook/report
기반으로 좋은 설계다.

## 4. 현재 수치와 기준 논문 대조

JIT warm 상태에서 현재 코드를 직접 실행했다. 시간은 이 Windows 환경의 median이다.

| 조건 | 현재 결과 | 시간 |
|---|---|---:|
| EIT default, Doppler off | `T(0)=0.94002`, contrast `0.13539`, FWHM `1.61306 MHz` | compute `10.19 ms`, readout `0.62 ms` |
| AT default, Doppler off | `T(0)=0.80603`, split `3.4965 MHz`, optimum `-2.220 MHz` | compute `10.20 ms`, finite-IF readout `203.20 ms` |
| AT default sensitivity | PSN = total = `6.08239 nV/cm/√Hz` | RIN/electronics 기본값이 0이므로 동일 |

논문 값을 코드에 강제로 주입하지 않는 선택 자체는 옳다
(`tests/test_rydberg_eit.py:64-81`). 그러나 현재 테스트는 “논문 값과 다르다”만 확인하고,
그 차이가 실험 reference로 허용 가능한지는 검증하지 않는다.

### 4.1 새 P1: 내장 온도 최적점은 논문의 33 °C를 재현하지 않는다

논문은 계산된 PSN sensitivity가 약 33 °C에서 최소이고, 20 °C·6 µW에서
`11.2 nV/cm/√Hz`라고 보고한다. 현재 내장 10-point cell-heating extra view는 다음을 낸다.

| 온도 | 현재 PSN/total sensitivity |
|---:|---:|
| 20 °C | `6.0824 nV/cm/√Hz` |
| 30 °C | `2.3905 nV/cm/√Hz` |
| 35 °C | `1.6312 nV/cm/√Hz` |
| 40 °C | `1.1896 nV/cm/√Hz` |
| **47 °C** | **`0.9717 nV/cm/√Hz` (내장 최적점)** |
| 55 °C | `1.1517 nV/cm/√Hz` |

온도 grid와 최적점 선택은 `gabes/schemes/rydberg.py:1330-1411`에 있다. 논문의
33 °C와는 14 °C 차이고, 20 °C scale은 논문 PSN보다 약 `1.84×` 낙관적이다. 더구나
example config는 temperature/density dephasing을 모두 0으로 둔다
(`analysis/rydberg_cell_heating/example_config.json:5-27`).

이 결과는 **현재 extra view가 잘못 계산한다는 단순 수치 버그**라기보다, opacity,
dephasing, Zeeman participation, detector reference가 아직 논문 curve에 검증되지 않은
조건부 모델임을 뜻한다. 따라서 이 sweep를 온도 최적화의 절대 reference로 쓰면 안 된다.

가장 싼 개선은 현재 sweep에 논문 기준점과 `UNVALIDATED / qualitative trend only` 상태를
함께 표시하는 것이다. 이미 푼 배열의 O(1) 비교라 추가 OBE solve가 없다. 실제 curve를
맞출 때는 paper/measurement의 linewidth·transmission·sensitivity를 동시에 사용해 기존
temperature/density dephasing fit을 calibration해야 하며
(`gabes/rydberg_experiment.py:593-795`), opacity 하나로 결과를 억지로 맞추면 안 된다.

### 4.2 새 근거: probe-power optimum도 외부 검증이 필요하다

0.5–20 µW diagnostic에서 현재 PSN sensitivity는 `4 µW`에서 `5.9120`으로 최소이고,
`6 µW`에서는 `6.0824 nV/cm/√Hz`였다. 방향성은 논문의 “저 power shot noise와 고 power
broadening 사이 최적점”을 재현하지만, 논문의 measured optimum `6 µW`, measured
`12.5(8)`, PSN `11.2 nV/cm/√Hz`와 scale 및 최적 위치가 다르다. 기존 테스트는 probe
power가 linewidth를 단조 증가시키는지만 확인한다 (`tests/test_rydberg_eit.py:236-246`).

권고는 exact paper target을 hard-code하는 것이 아니라, probe-power sensitivity curve에
reference overlay/status와 model/reference optimum ratio를 추가하는 것이다. 현재 finite-IF
sweep 결과의 O(N) 후처리이고, 실측 CSV가 들어오면 같은 validation contract를 재사용할 수 있다.

## 5. 기존 미완료 제안의 현재 비용과 우선순위

### P1. 숨은 `ls=0.001` effective opacity

`compute()`는 `ls=0.001`을 고정 반환한다 (`gabes/schemes/rydberg.py:625-643`). 이 값은
⁸⁵Rb D2 `F=3→F'=4`의 atomic line strength라기보다 velocity/Zeeman/polarization,
optical pumping, overlap, fit amplitude가 섞인 effective opacity anchor다. 2026-08-02
진단에서 `ls=0.0005/0.001/0.002`는 default sensitivity를 각각
`11.782/6.082/3.239 nV/cm/√Hz`로 바꿨다. 오늘 확인한 20 °C 논문 차이와 직접 연결된다.

`atom_participation_fraction`은 표시용 `N_eff`에만 들어가고 OBE/감도에는 들어가지 않는다
(`gabes/schemes/rydberg.py:732-746`). `ls`를 atomic `C_F²`와 분리된
`effective_opacity_participation`으로 공개하고 uncertainty sweep에 넣는 것은 solver 차원을
늘리지 않는다. atomic phasor를 cache하면 transmission/noise 후처리만 다시 하면 된다.

### P1. extra-view cache dependency가 여전히 잘못돼 있다

Scheme contract상 `recompute=False`는 heavy solve에 쓰이지 않아야 한다
(`gabes/schemes/base.py:7-17`). 앱은 extra view에 recompute item만 전달한다
(`streamlit_app.py:608-619`, `1455-1458`, `1507-1513`). 그러나 cell-heating extra view는
`if_khz`, detector QE/path/reference/RIN/electronics, RF dipole/angular factor, cell length를
실제 finite-IF sensitivity에 사용한다 (`gabes/schemes/rydberg.py:1330-1411`).

2026-08-02의 app-equivalent 진단에서는 사용자가 IF `200 kHz`, QE `0.2`, path `0.2`로
바꿔도 extra view가 default `6.0824`를 냈고, full params는 `21.6711`을 냈다. 현재 코드가
그 뒤 바뀌지 않았으므로 문제는 그대로다. `ExtraView`별 dependency key와 live metadata를
도입하고 atomic phasor/detector postprocess cache를 분리해야 한다. 물리식은 바뀌지 않으며,
오히려 사용자가 선택한 물리를 정확히 적용하는 수정이다.

### P1. 저투과·edge optimum에서도 숫자가 hero가 된다

`_readout()`은 EIT contrast를 계산한 뒤 버리고, slope/IF optimum을 전체 scan에서 찾는다
(`gabes/schemes/rydberg.py:952-973`, `1030-1051`, `1137-1159`). 이전 80 °C EIT 진단의
`T(0)=1.799×10⁻⁶`, contrast `1.799×10⁻⁶`, FWHM `0.307 MHz`, IF optimum `-8.96 MHz`
조건에서도 linewidth가 hero였다. contrast, detected power, samples-per-width, edge distance를
표시하고 central-feature window 밖 optimum을 status 처리하는 것은 O(Nscan) 후처리라
추가 solve가 없다.

### P2. reference arm이 detuning마다 이상적으로 재균형된다

각 probe detuning에서 `reference_power = signal_power × ratio`를 새로 만들고 동일 DC balance를
적용한다 (`gabes/schemes/rydberg.py:849-885`). 이는 한 개의 고정 reference arm이 아니라
detuning별 최적 detector envelope다. fixed reference power 또는 balance detuning을 입력받아
scan 전체에 유지해야 한다. detector loop만 바뀌므로 원자 solve overhead는 0이다.

### P2. AT peak와 sensitivity optimum 보간

AT peak/optimum은 grid point다 (`gabes/schemes/rydberg.py:662-666`, `894-906`,
`994-1005`). 3-point quadratic interpolation과 curvature/resolution status는 O(1)이며
추가 solve 없이 grid stair-step을 줄인다.

### P2. RF dipole provenance를 실행 가능하게 보존

`1326.257243... e a₀`는 ARC 3.10.2 stretched-state σ+ 결과라는 주석과 config 값만 있고,
그 값을 재생성하는 script/test는 없다 (`gabes/schemes/rydberg.py:29-33`;
`analysis/rydberg_cell_heating/example_config.json:88-96`). ARC는 fine/hyperfine-resolved
matrix-element API를 별도로 제공하므로, 사용한 isotope/state/`mJ`/polarization 호출,
ARC version과 결과를 작은 provenance script로 고정하는 편이 좋다
([ARC matrix-element documentation](https://arc-alkali-rydberg-calculator.readthedocs.io/en/latest/generated/arc.alkali_atom_functions.AlkaliAtom.getDipoleMatrixElementHFS.html)).
앱 runtime overhead는 0이다.

## 6. 기존 제안의 계산비용과 물리 보존 판단

| 항목 | 상태 | 계산비용 | 실험물리 판단 |
|---|---|---:|---|
| power/diameter→Rabi, transit, scalar temperature/density dephasing | 완료 | O(1) 전처리 | 같은 4준위 OBE를 유지하며 knob 방향성을 유용하게 보존한다. anchor임은 계속 밝혀야 한다. |
| finite-IF weak-SIG + detector noise | 완료 | default AT readout 약 `203 ms` | full time-domain보다 훨씬 싸고 실제 linear response다. negligible overhead는 아니지만 interactive 사용에 충분하다. |
| heater/effective/cold-spot 분리 | 완료 | O(1) | 실험 의미가 크고 solve 차원을 늘리지 않는다. |
| axial Beer–Lambert / `N_eff` | 완료, approximate | O(Nz·Nscan), local OBE 재해 없음 | density-rescaled diagnostic으로만 적절하다. spatial sensitivity라 부르면 안 된다. |
| SAM field/uncertainty | 완료 | O(1) | far-field scalar calibration으로 유용하나 horn/cell standing wave는 보존하지 않는다. |
| opacity 공개, validation status, contrast/edge, peak interpolation | 미완료 | O(1)–O(Nscan) | 물리를 잃지 않고 거의 공짜로 신뢰도를 높이는 최우선 작업이다. |
| low-order polarization/Zeeman proxy | deferred GROUP B | scalar/소차원, 같은 solve 가능 | 측정으로 calibration하고 phenomenological uncertainty band로 표시해야 한다. |
| full Zeeman/polarization/ionization | deferred GROUP C | 상태수와 dense Liouvillian 급증 | 절대 calibration에는 중요하지만 negligible-overhead가 아니다. |
| full time-domain LO+SIG/lock-in | deferred GROUP C | time grid/filter/noise process | 별도 엔진과 실험 protocol 합의가 필요하다. |
| 3-D RF field / three-photon imaging | deferred GROUP C | 공간 격자·EM/이미지 모델 | 현재 scalar SAM/4-level scheme의 작은 확장이 아니다. |

## 7. 동작을 바꾸지 않는 순수 코드 최적화

### 7.1 finite-IF scan의 affine Liouvillian assembly

`_superheterodyne_readout()`은 801 detuning마다 Python에서 Hamiltonian과 Liouvillian을
새로 만든다 (`gabes/schemes/rydberg.py:811-818`). detuning dependence는 affine이므로
`L(s)=L(0)+s[L(1)-L(0)]`을 broadcasting할 수 있다.

오늘 재측정:

- 현재 list/build median: `103.30 ms`
- affine broadcast median: `3.856 ms`
- 약 `26.8×` 빠름
- `np.array_equal=True`, max absolute difference `0`

수치와 물리를 bit-identical하게 유지하면서 default AT readout의 큰 병목을 제거한다.

### 7.2 static steady state를 finite-IF response에 재사용

`weak_signal_response_from_liouvillian()`은 steady state를 다시 푼 뒤 두 sideband를 푼다
(`gabes/rydberg_electrometry.py:165-188`). 그러나 static spectrum은 같은 801개 Liouvillian의
steady state를 이미 계산했다. affine kernel이 optional full `rho₀`를 반환해 finite-IF 호출에
전달하면 된다.

오늘 재측정:

- 현재 response 전체: `35.834 ms`
- 미리 계산한 `rho₀` 전달: `25.284 ms`
- response 단계 약 `29.4%` 절감
- 두 sideband 배열 모두 `np.array_equal=True`

801×4×4 complex `rho₀`는 약 0.2 MiB라 memory cost도 작다. 단, 기본 static path에서
full density matrix를 항상 materialize하지 말고 finite-IF 요청에만 opt-in 반환하는 것이 좋다.

### 7.3 atomic phasor와 detector/SAM postprocess cache 분리

RF dipole, QE, path efficiency, reference balance, RIN, electronics, ENBW, SAM은
`L₀ ± iω_IF`를 바꾸지 않는다 (`gabes/schemes/rydberg.py:837-923`). OBE+IF-dependent complex
phasor 배열만 cache하고 detector loop를 다시 계산하면 calibration slider가 즉시 반응한다.
이는 기존 cache-correctness 수정과 같은 경계에서 처리해야 한다.

### 7.4 공통 affine Doppler loop hoist

`kernels._affine_scan_chi_real()`은 velocity마다 `base+sA+kvB`를 다시 채운다
(`gabes/kernels.py:426-449`). `base+sA`를 velocity loop 밖으로 옮기는 기존 2026-07-22
제안은 현재도 유효하다. 당시 `190→166 ms`, bit-identical이었다. Rydberg Doppler-on뿐 아니라
Lambda에도 공통 이득이다.

## 8. 문서·테스트·예제의 실험 reference 적합성

### 좋은 점

- README는 finite-IF, conditional sensitivity, temperature/cold spot, SAM, deferred scope를
  현재 코드와 대체로 일치하게 설명한다 (`README.md:94-116`).
- 분석 README는 evidence status와 approximation boundary를 명확히 분리한다
  (`analysis/rydberg_cell_heating/README.md:22-39`, `78-117`).
- 테스트는 기준 linewidth/AT split, linear-response 방정식, noise convention, 온도/cold spot,
  axial integration, SAM uncertainty, CSV units/provenance를 폭넓게 검사한다
  (`tests/test_rydberg_eit.py:29-439`; `tests/test_rydberg_electrometry.py:38-203`;
  `tests/test_rydberg_experiment.py:19-173`; `tests/test_rydberg_experimental_csv.py:14-109`).
- example config는 raw EIT/RF/temperature 입력을 `PENDING`, RF dipole을 `REFERENCE`,
  detector/geometry를 `ASSUMED`, 출력은 `PREDICTED`로 둔다
  (`analysis/rydberg_cell_heating/example_config.json:53-103`). 실제 데이터 없이 validation
  완료처럼 보이게 하지 않는다.

### 보완할 점

- 사용자 가이드는 여전히 static spectrum, linewidth, AT split, dispersion만 설명해
  finite-IF/noise/cold-spot/SAM 계층을 누락한다
  (`docs/Userguide/GABES_User_Guide_v2.html:534-537`, `595-605`).
- `tests/test_rydberg_eit.py:1-5`도 “public UI shows the static spectrum only”라고 해 현재와
  모순된다.
- bundled example의 raw EIT/RF/PSD/temperature 파일은 모두 `null/PENDING`이다
  (`analysis/rydberg_cell_heating/example_config.json:53-78`). pipeline demonstration이지
  experiment-model agreement 증거가 아니다.
- 현재 test는 논문의 33 °C temperature optimum, 6 µW sensitivity optimum, hidden opacity,
  fixed-reference detector, high-temperature trust status, Streamlit extra-view dependency를
  검증하지 않는다. 내부 수학 regression은 강하지만 외부 실험 validation contract가 약하다.

## 9. 검증

- Scheme 3 및 cell-heating 관련: `60 passed in 26.47 s`
- 저장소 전체 회귀: `378 passed in 143.79 s`
- benchmark와 physics diagnostic은 production code를 수정하지 않고 실행했다.

## 10. 결론과 권장 순서

현재 Scheme 3은 **실제 물리를 담은 유용한 반정량적 실험 설계 도구**다. compact 4-level
OBE, power/beam/transit dependence, AT field conversion, finite-IF complex response,
balanced-noise budget, cold-spot/axial/SAM/provenance 계층은 모두 실험자가 knob 방향과
검출 chain을 이해하는 데 가치가 있다.

그러나 **절대 atomic electrometry 또는 cell-heating 최적화 reference**로 인용하기에는 아직
부족하다. 가장 직접적인 근거는 현재 모델이 논문의 20 °C PSN scale을 약 1.84배 낙관적으로
예측하고, 내장 온도 optimum을 논문의 약 33 °C가 아닌 47 °C로 고른다는 점이다. 이 차이는
숨은 `ls=0.001`, zero-default density/temperature dephasing, stretched-state RF dipole,
idealized detector balance가 아직 외부 calibration curve로 검증되지 않았음을 보여 준다.

권장 순서는 다음과 같다.

1. **P1:** extra-view dependency/cache와 live metadata를 바로잡고 app-path 회귀 test를 넣는다.
2. **P1:** temperature/probe-power sweep에 paper comparison과 `UNVALIDATED` status를 추가한다.
3. **P1:** `ls=0.001`을 effective opacity calibration으로 공개하고 uncertainty를 감도에 전파한다.
4. **P1:** contrast/transmitted-power/edge/resolution status로 hero 숫자를 보호한다.
5. **P2:** fixed reference arm, AT/optimum interpolation, 실행 가능한 ARC provenance를 추가한다.
6. **속도:** affine finite-IF assembly, static `rho₀` 재사용, atomic-phasor cache를 적용한다.
7. 그 뒤 실제 EIT/RF/PSD/temperature trace와 detector/SAM calibration으로 저차 Zeeman proxy를
   검증하고, full Zeeman/time-domain/spatial GROUP C의 필요성을 판단한다.

1–6은 solver 차원을 늘리지 않거나 결과를 그대로 보존하는 저비용 작업이다. full Zeeman이나
full lock-in보다 먼저 처리해야 현재 구현된 좋은 물리가 과도한 정밀도나 잘못된 최적점으로
전달되지 않는다.
