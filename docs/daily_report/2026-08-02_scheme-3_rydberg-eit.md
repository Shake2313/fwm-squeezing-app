# 2026-08-02 Scheme 3 물리 검토: Rydberg-EIT electrometry

## 1. 오늘의 선택과 현재 5개 scheme 순서

- 로컬 날짜는 2026-08-02이고, `n = (2 mod 5) + 1 = 3`이다.
- 드롭다운 순서는 `gabes/schemes/__init__.py:19-25`의 `_SCHEMES`가 결정한다.
- 따라서 오늘의 대상은 **Scheme 3, `rydberg_eit` — Rydberg-EIT electrometry**이다
  (`gabes/schemes/rydberg.py:91-100`).

| 번호 | registry 객체 | 현재 표시 이름 | 근거 |
|---:|---|---|---|
| 1 | `SASScheme()` | Absorption spectroscopy (OD / SAS) | `gabes/schemes/sas.py:53-56` |
| 2 | `LambdaScheme()` | Λ coherence (EIT / AT / CPT) | `gabes/schemes/absorption.py:462-503` |
| 3 | `RydbergEITScheme()` | Rydberg-EIT electrometry | `gabes/schemes/rydberg.py:91-100` |
| 4 | `MagnetoScheme()` | Magneto-optics (Hanle/MOR) | `gabes/schemes/magneto.py:127-172` |
| 5 | `FWMScheme()` | Four-wave mixing (Squeezing / Biphoton) | `gabes/schemes/fwm.py:1800-1803` |

평가 대상은 `main`의 `a82bbf7` 위 현재 working tree다. Rydberg finite-IF,
cell-heating, CSV 분석 계층은 아직 커밋되지 않은 현재 파일까지 포함해 읽었으며, 이 보고서
외의 기존 변경은 수정하지 않았다.

## 2. 먼저 찾은 기존 제안, TODO, daily report, issue note

기존 개선안은 **있다**. 새 제안을 만들기 전에 다음을 먼저 대조했다.

- `docs/checklist.json:8-40`에는 power-to-Rabi, finite-IF electrometry,
  heated-cell model, experimental report, axial/SAM 계층이 모두 `done`으로 기록돼 있다.
- `docs/checklist.json:116-120`의 저차 polarization/Zeeman proxy는 `deferred`다.
- `docs/checklist.json:158-183`의 full Zeeman/polarization/field/ionization,
  time-domain superhet, 3-D RF field, three-photon imaging은 명시적인 GROUP C다.
- Scheme 3의 기존 일일 보고서는 `docs/daily_report/2026-06-22_scheme-3_rydberg-eit.md`,
  `2026-06-27`, `2026-07-02`, `2026-07-12`, `2026-07-22`를 확인했다.
- 특히 2026-07-22 보고서가 지적한 숨은 `line_strength=0.001`, 고온 저대비 hero,
  AT peak grid 정밀도 문제는 `docs/daily_report/2026-07-22_scheme-3_rydberg-eit.md:107-176`에
  이미 제안돼 있다. 오늘 코드에서도 세 항목은 아직 그대로다.
- 저장소 안에서는 별도 TODO/issue-note 파일을 찾지 못했다. 공개
  [GitHub Issues](https://github.com/Shake2313/fwm-squeezing-app/issues?q=is%3Aissue)도
  2026-08-02 확인 시 0건이다.

즉 “기존 제안이 없다”가 아니라, **과거 저비용 신뢰성 제안 일부는 미완료이고,
그 사이 더 큰 finite-IF/열·검출 계층이 구현된 상태**다.

## 3. 현재 구현이 담는 실제 물리

### 3.1 정적 원자·광학 모델

현재 원자 코어는 ⁸⁵Rb
`5S1/2 F=3 → 5P3/2 F'=4 → 40D5/2 → 39F7/2`의 축약 4준위 cascade다.
5P spontaneous decay, 작은 Rydberg decay, 5S–40D와 40D–39F coherence dephasing,
780 nm probe, 481 nm coupling, 37 GHz RF LO를 정적 OBE에 넣는다
(`gabes/schemes/rydberg.py:4-15`, `69-88`, `413-433`, `554-623`).

- probe와 coupling Rabi는 power와 1/e² beam diameter에 대해
  `Ω ∝ √P/d`로 변하며 reference point에 anchor된다
  (`gabes/schemes/rydberg.py:435-457`).
- beam diameter와 effective vapor temperature는 transit broadening에도 들어간다
  (`gabes/schemes/rydberg.py:493-501`).
- residual two-photon Doppler는 per-level wave-vector ratio로 표현하고, `doppler=on`이면
  Maxwell velocity class를 평균한다 (`gabes/schemes/rydberg.py:503-523`, `606-615`).
- optical coherence는 physical susceptibility와 Beer–Lambert transmission으로 변환된다
  (`gabes/schemes/rydberg.py:668-677`; `gabes/observables.py:393-418`).
- EIT FWHM은 half-height crossing을 선형 보간하고, AT는 local transparency maxima로
  분리를 읽는다 (`gabes/schemes/rydberg.py:925-950`, `981-1029`).

이 부분은 EIT linewidth, power broadening, RF Autler–Townes splitting, microwave detuning에
따른 center shift를 실험 전 knob-scan하는 데 실제로 유용하다. 1.6 MHz EIT linewidth와
37 GHz에서 12.5(8) nV/cm/√Hz를 보고한
[Ju et al., arXiv:2606.04354](https://arxiv.org/abs/2606.04354)의 기준점도 명시적으로
문서화돼 있다 (`gabes/schemes/rydberg.py:1495-1527`).

### 3.2 finite-IF weak-SIG와 detector chain

새 electrometry 계층은 단순 정적 slope proxy보다 한 단계 진전됐다.

- LO-dressed Liouvillian 주위에서 `L₀ ± iω_IF`의 두 sideband 선형응답을 푼다
  (`gabes/rydberg_electrometry.py:116-199`).
- 복소 velocity-class phasor를 magnitude가 아니라 coherent average한다
  (`gabes/rydberg_electrometry.py:69-113`, `216-232`).
- RF dipole, angular/polarization factor, peak/RMS convention을 SI field와 Rabi 사이에
  명시적으로 둔다 (`gabes/rydberg_electrometry.py:235-321`).
- signal/reference photodiode의 one-sided Schottky ASD, correlated RIN, electronics ASD를
  합치고 field ASD로 환산한다 (`gabes/rydberg_electrometry.py:324-488`).
- Scheme은 probe-detuning scan에서 response 자체가 아니라 total noise-equivalent field가
  최소인 점을 고른다 (`gabes/schemes/rydberg.py:769-923`).

전용 테스트는 sideband 방정식, Hermiticity, low-IF static derivative limit, high-IF roll-off,
peak/RMS round trip, balanced noise를 직접 검사한다
(`tests/test_rydberg_electrometry.py:38-186`). 따라서 **수학적 선형응답 primitive 자체는
실제 물리이고 내부 일관성도 좋다**. 다만 full LO+SIG waveform, mixer/lock-in phase/filter,
sampling과 colored technical noise는 구현하지 않았다는 경계를 코드와 문서가 솔직히 밝힌다
(`gabes/rydberg_electrometry.py:1-35`; `docs/checklist.json:165-169`).

### 3.3 온도, cold spot, axial density, SAM, 실험 CSV

실험가 관점에서 좋아진 부분이다.

- heater setpoint, sensor, effective vapor temperature, cold spot을 분리하고 출처 문자열을
  보존한다 (`gabes/rydberg_experiment.py:47-162`).
- sealed cell에서는 cold spot이 vapor pressure를 정하고, hot optical region의 local density를
  ideal-gas 비로 낮춘다 (`gabes/schemes/rydberg.py:459-491`).
- `T(z), n(z)`, column density, segmented Beer–Lambert, Gaussian overlap `N_eff`가 별도
  primitive로 존재한다 (`gabes/rydberg_experiment.py:165-443`).
- SAM은 `E_RMS = √(30PG)/r`, peak/RMS, 표준불확도, `r/(2D²/λ)` 경고를 명시한다
  (`gabes/rydberg_experiment.py:446-589`).
- EIT/RF/PSD CSV는 unit 변환, duplicate median, SHA-256 provenance를 보존하고 signal을
  임의 정규화하지 않는다 (`gabes/rydberg_experimental_csv.py:1-7`, `211-379`, `446-718`).
- batch workflow는 `MEASURED/PPT/REFERENCE/FITTED/PREDICTED/ASSUMED/PENDING`을 분리한다
  (`analysis/rydberg_cell_heating/README.md:22-39`). axial spectrum은 local OBE 재해가 아니라
  density-rescaled lumped line shape임을 명시한다 (`analysis/rydberg_cell_heating/README.md:86-104`).

이는 실험 notebook/report 인프라로서 좋은 방향이다. 특히 setpoint를 vapor temperature나
cold spot과 자동 동일시하지 않는 것은 실제 heated-cell 분석에서 중요하다.

## 4. 현재 기준점과 수치 진단

JIT warm 상태에서 현재 default를 직접 실행했다. 시간은 이 Windows 환경의 median이며
절대 benchmark라기보다 병목 비율을 보는 값이다.

| 조건 | 주요 결과 | 시간 |
|---|---|---:|
| EIT, Doppler off | `T(0)=0.94002`, contrast `0.13539`, FWHM `1.61306 MHz`, max slope `0.122/MHz` | static compute `6.07 ms`, readout `0.49 ms` |
| AT, Doppler off | `T(0)=0.80603`, split `3.50 MHz`, optimum `-2.220 MHz` | static compute `21.45 ms`, finite-IF readout `289.36 ms` |
| AT default sensitivity | total = PSN = `6.0824 nV/cm/√Hz` | RIN/electronics default가 0이므로 동일 |
| 10점 cell-heating extra view | 10개 EIT + 10개 AT/finite-IF point | `2.80 s` (default 입력) |

현재 값은 논문의 12.5 또는 11.2를 직접 대입한 결과는 아니다
(`tests/test_rydberg_eit.py:64-81`). 그러나 아래 hidden opacity와 ideal detector assumption에
강하게 의존하므로, “독립적인 절대 예측”보다 **조건부 모델 출력**으로 읽어야 한다.

## 5. 실험물리 reference로서 가장 중요한 신뢰 경계

### P1. 숨은 `ls=0.001`이 이제 절대 감도까지 지배한다

`compute()`는 여전히 `ls=0.001`을 raw result에 고정한다
(`gabes/schemes/rydberg.py:625-643`). 같은 저장소의 ⁸⁵Rb D2
`F=3→F'=4` AutoOD `C_F²`는 1.0이다 (`gabes/species.py:226-245`). 따라서 0.001은
원자 선강도가 아니라 stationary velocity fraction, mF/polarization participation,
optical pumping, mode overlap, fitted opacity가 섞인 **effective opacity anchor**다.

이 factor만 바꾼 진단은 다음과 같다.

| effective `ls` | AT `T(0)` | PSN/total sensitivity |
|---:|---:|---:|
| 0.0005 | 0.89779 | 11.782 nV/cm/√Hz |
| 0.0010 | 0.80603 | 6.082 nV/cm/√Hz |
| 0.0020 | 0.64968 | 3.239 nV/cm/√Hz |

즉 factor 2의 의미 선택만으로 감도가 약 factor 1.9 변한다. 그런데 이 값은 UI, derived
table, uncertainty budget, example config 어디에도 calibration input으로 드러나지 않는다.
현재 `atom_participation_fraction`은 오직 표시용 `N_eff`에만 들어가고 OBE/감도에는 들어가지
않는다 (`gabes/schemes/rydberg.py:732-746`).

**권고:** 기존 2026-07-22 제안대로 이름을 `effective_opacity_participation`처럼 바꾸고,
atomic `C_F²`와 분리해 table/help/config에 노출하며 uncertainty/sensitivity sweep에 포함한다.
이는 scalar scaling과 이미 계산한 spectrum 재매핑이므로 solver overhead가 사실상 0이다.
이 작업 전에는 6.08 nV/cm/√Hz를 절대 atomic reference로 인용하지 않는 편이 안전하다.

### P1. UI cache dependency가 extra-view 물리를 무시한다

Scheme contract는 `recompute=False` knob가 solve에 들어가지 않고 observables가 싸야 한다고
정한다 (`gabes/schemes/base.py:7-17`). app은 extra view에도 recompute item만 전달한다
(`streamlit_app.py:608-619`, `1417-1475`). 그러나 현재 Rydberg에서는:

- `if_khz`가 `recompute=False`지만 실제 `L₀ ± iω_IF` solve에 들어간다
  (`gabes/schemes/rydberg.py:253-256`, `788-818`).
- detector QE/path/reference/RIN/electronics, RF dipole/angular factor, cell length가 모두
  cell-heating sensitivity extra view 결과에 필요하지만 cache key와 extra compute 입력에서 빠진다
  (`gabes/schemes/rydberg.py:278-323`, `1330-1411`).
- `mw_frequency_ghz`는 “display only”인데 `compute()`의 raw에 저장한 뒤 table이 raw를 읽으므로,
  app solve-cache 경로에서는 사용자가 바꾼 표시값이 갱신되지 않는다
  (`gabes/schemes/rydberg.py:222-224`, `654-655`, `1225-1245`).

실제 app과 같은 `recompute_items`만 extra view에 넘긴 진단에서, 사용자가
`IF=200 kHz`, `QE=0.2`, `path=0.2`로 바꿔도 sweep 첫 점은 default `6.0824`를 냈다.
full params를 직접 넘기면 `21.6711 nV/cm/√Hz`였다. `mw_frequency=55 GHz`도 cached raw에는
37 GHz로 남았다. 테스트는 extra view에 full params를 직접 넘기므로 이 UI 경로를 잡지 못한다
(`tests/test_rydberg_eit.py:419-431`).

**권고:** `ExtraView`별 dependency key를 도입하고, finite-IF atomic phasor는 OBE+IF key로
별도 cache하며 detector/SAM scaling만 navigate postprocess로 분리한다. 표시 metadata는 live
params에서 읽는다. 물리식은 그대로이며, 올바른 설정을 쓰게 만드는 캐시 수정이다.

### P1. 저대비·불투명 조건에서도 linewidth가 hero다

80 °C EIT 진단은 현재도 다음 결과를 냈다.

- `T(0)=1.799×10⁻⁶`, central contrast `1.799×10⁻⁶`
- 보고 FWHM `0.3070 MHz`
- IF optimum `-8.96 MHz`, scan edge `-9 MHz`
- hero는 여전히 `EIT linewidth`, `Transmission at resonance`

`_readout()`이 `_eit_features()`의 contrast를 버리고 finite width만 hero 조건으로 쓰며,
slope/IF를 전체 scan에서 찾기 때문이다 (`gabes/schemes/rydberg.py:952-973`, `1030-1051`,
`1137-1159`). 이는 중앙 Rydberg feature가 아니라 background/edge를 계측 optimum으로
오인할 수 있다.

**권고:** contrast, detected probe power, samples-per-width를 표시하고, central feature window
안에서 slope/IF를 찾으며 edge/low-contrast status를 hero보다 우선한다. 모두 기존 배열의
O(Nscan) 후처리이므로 추가 OBE solve가 없다.

### P2. reference arm은 detuning마다 이상적으로 재균형된 envelope다

main sensitivity scan은 각 후보 detuning에서 `reference_power = signal_power × ratio`로 다시
만들고 `reference_weight=1/ratio`로 DC balance한다
(`gabes/schemes/rydberg.py:849-885`). 따라서 한 개의 고정 reference arm을 scan한 것이 아니라,
각 detuning마다 이상적으로 다시 맞춘 detector들의 sensitivity envelope다. default RIN=0,
electronics=0, dark current UI 없음도 낙관적이다.

**권고:** fixed reference optical power 또는 한 balance detuning을 입력받고 scan 전체에서
동일 reference를 사용하며, photodiode saturation/headroom과 dark current를 table에 표시한다.
detector loop만 바뀌므로 원자 solve overhead는 0이다. 이 개선은 full lock-in model보다 훨씬
싸면서 실제 balanced detection reference로서의 의미를 크게 높인다.

### P2. AT 위치와 optimum은 아직 grid 좌표다

AT maxima와 finite-IF optimum은 grid point를 그대로 쓴다
(`gabes/schemes/rydberg.py:661-666`, `994-1005`, `894-906`). 현재 default scan에서:

| 점수 | 간격 | AT split | sensitivity optimum |
|---:|---:|---:|---:|
| 401 | 0.0555 MHz | 3.44 MHz | -2.220 MHz |
| 801 | 0.02775 MHz | 3.50 MHz | -2.220 MHz |
| 1601 | 0.013875 MHz | 3.50 MHz | -2.206 MHz |

default 결과는 두 자리 표시에서 대체로 안정적이지만, 큰 LO에서는 기존 보고서가 보인
한-grid-step 오차가 더 커진다. 3-point quadratic interpolation과 resolution status는 추가
solve 없이 가능하다.

## 6. 기존 제안의 비용과 물리 보존 평가

| 기존 항목 | 현재 상태 | 계산비용 | 평가 |
|---|---|---:|---|
| power/waist→Rabi, transit, temperature/density scalar dephasing | 완료 | O(1) 전처리 | 같은 4준위 OBE를 유지하며 비용 대비 유용성이 높다. 다만 anchor/phenomenology임을 유지해야 한다. |
| finite-IF weak-SIG + noise chain | 완료 | default AT readout 약 289 ms | 실제 동적 선형응답을 보존한다. negligible overhead는 아니지만 full time-domain보다 훨씬 싸다. |
| heater/effective/cold-spot 분리 | 완료 | O(1) | 실험 의미가 크고 solve 차원을 늘리지 않는다. |
| axial Beer–Lambert / `N_eff` | 완료, approximate | O(Nz·Nscan), local OBE 재해 없음 | 현재 명시처럼 density-rescaled diagnostic으로는 적절하다. spatial sensitivity라고 부르면 안 된다. |
| SAM field/uncertainty | 완료 | O(1) | far-field 조건 안에서는 유용하다. near-field/standing wave를 보존하지는 못한다. |
| effective opacity 공개, contrast/edge status, peak interpolation | 미완료 | O(1) 또는 O(Nscan) | 물리를 잃지 않고 사실상 공짜로 신뢰도를 높이는 최우선 항목이다. |
| low-order polarization/Zeeman proxy | deferred GROUP B | scalar/소차원, 거의 동일 solve | full manifold보다 싸지만 데이터로 calibration하고 phenomenological band라고 표시해야 한다. |
| full Zeeman/polarization/ionization | deferred GROUP C | 상태수 증가, dense Liouvillian 급증 | 절대 calibration에는 중요하지만 negligible-overhead가 아니다. |
| full time-domain LO+SIG/lock-in | deferred GROUP C | time grid·filter·noise process 추가 | 목적과 detector protocol 합의가 필요한 별도 엔진이다. |
| 3-D RF field / three-photon imaging | deferred GROUP C | 공간 격자·EM/이미지 모델 | 현재 scalar SAM/4-level scheme에 작은 패치로 넣을 수 없다. |

## 7. 동작을 바꾸지 않는 순수 코드 최적화

### 7.1 finite-IF scan의 affine Liouvillian assembly

`_superheterodyne_readout()`은 801개 detuning마다 Python에서 Hamiltonian과 Liouvillian을
새로 만든다 (`gabes/schemes/rydberg.py:811-818`). 이 Hamiltonian은 scan detuning에 affine다.
따라서 static path처럼

`L(s) = L(0) + s·[L(1)-L(0)]`

을 한 번 만든 뒤 broadcasting하면 된다. 진단 결과:

- 현재 list/build: median `142.44 ms`
- affine broadcast: median `5.00 ms`
- max absolute difference `0`, `np.array_equal=True`

조립 단계만 약 28.5배 빨라지고, 전체 default AT readout 289 ms 중 상당 부분을 없앤다.
finite-IF solve 자체는 그대로라 behavior와 물리는 변하지 않는다.

### 7.2 atomic phasor와 detector postprocess cache 분리

RF dipole, QE, path efficiency, reference balance, RIN, electronics, ENBW, SAM은
`L₀ ± iω_IF`를 바꾸지 않는다. 현재는 이 navigate knob 하나를 바꿔도 observables가 801점
atomic solve를 다시 수행한다. OBE+IF-dependent complex phasor array를 별도 cache하고 detector
loop만 다시 계산하면 결과를 바꾸지 않으면서 실험 calibration slider가 즉시 반응한다.
이는 위 UI dependency 수정과 함께 하는 것이 안전하다.

### 7.3 기존 affine Doppler loop hoist는 아직 유효

`kernels._affine_scan_chi_real()`은 velocity마다 `base + s*A + kv*B` 전체를 다시 채운다
(`gabes/kernels.py:426-449`). `base+s*A`를 velocity loop 밖으로 옮기는 기존 제안은 현재도
미구현이다. 2026-07-22 진단의 `190→166 ms`, bit-identical 결과를 그대로 재검토할 가치가
있다 (`docs/daily_report/2026-07-22_scheme-3_rydberg-eit.md:232-239`). Doppler-on뿐 아니라
Lambda affine path에도 공통 이득이다.

### 7.4 작은 최적화

- 동일 dephasing의 `_atom()`/dissipator를 immutable bounded cache하면 probe-power sweep의
  반복 조립을 줄일 수 있다 (`gabes/schemes/rydberg.py:413-433`). 다만 ms 이하라 7.1보다
  우선순위가 낮다.
- CSV duplicate group median은 Python group loop다
  (`gabes/rydberg_experimental_csv.py:332-340`). 큰 duplicate-heavy trace에서는 vectorized
  grouped reduction 후보가 되지만 일반 Rydberg trace의 병목은 아니다.
- CSV는 sort 후 같은 x를 median merge하므로 bidirectional/hysteretic sweep를 한 곡선으로
  접을 수 있다. 최적화보다 먼저 pre-sort monotonic-run warning을 추가해야 하며 O(N)이다.

## 8. 문서·테스트·예제의 reference 적합성

### 잘된 점

- README는 finite-IF, conditional sensitivity, cold spot, SAM, deferred scope를 현재 코드와
  대체로 일치하게 설명한다 (`README.md:94-116`).
- 분석 README는 status와 approximation boundary를 명확히 분리한다
  (`analysis/rydberg_cell_heating/README.md:22-39`, `78-117`).
- tests는 reference defaults, 1.6 MHz linewidth, AT split, low-IF derivative, noise convention,
  temperature/cold spot, axial integration, SAM uncertainty, CSV units/provenance를 폭넓게 검사한다
  (`tests/test_rydberg_eit.py:29-439`; `tests/test_rydberg_electrometry.py:38-203`;
  `tests/test_rydberg_experiment.py:19-173`; `tests/test_rydberg_experimental_csv.py:14-109`).
- `analysis/rydberg_cell_heating/example_config.json:53-110`은 실제 measurement input을
  `PENDING`, dipole을 `REFERENCE`, detector/geometry를 `ASSUMED`, output을 `PREDICTED`로 둔다.
  실제 원시 데이터 없이 validation 완료처럼 보이게 하지 않는 점이 좋다.

### 보완할 점

- 사용자 가이드는 아직 “정적 스펙트럼, EIT linewidth, AT split, dispersion”만 설명해
  finite-IF/noise/cold-spot/SAM 계층을 반영하지 않는다
  (`docs/Userguide/GABES_User_Guide_v2.html:534-537`, `595-605`).
- `tests/test_rydberg_eit.py:1-5`의 module docstring도 “public UI shows the static spectrum only”로
  현재와 모순된다.
- 현재 test는 숨은 opacity 의미/uncertainty, high-temperature contrast/edge status,
  AT interpolation, fixed-reference detector, Streamlit recompute/extra-view dependency를 검증하지
  않는다. 내부 수학 테스트는 강하지만 **실험 calibration 신뢰성 test가 아직 부족하다**.
- example config의 raw EIT/RF/PSD/temperature inputs가 모두 `null/PENDING`이므로, 현재 bundled
  example은 pipeline demonstration이지 독립적인 experiment-model agreement 증거는 아니다
  (`analysis/rydberg_cell_heating/example_config.json:53-78`).

## 9. 검증

- Scheme 3 전용 및 cell-heating workflow 관련:
  `60 passed in 54.90 s`
- 저장소 정책의 전체 회귀:
  `232 passed in 94.29 s`
- benchmark/diagnostic은 production code를 수정하지 않고 실행했다.

## 10. 결론과 우선순위

현재 Scheme 3은 과거의 “정적 spectrum proxy”보다 분명히 강해졌다. compact 4-level OBE,
power/beam/transit dependence, AT field conversion, finite-IF complex response, balanced-noise budget,
cold-spot/axial/SAM/provenance 계층은 모두 실제 실험 설계와 데이터 정리에 유용하다.

그러나 **절대 atomic electrometry reference**로는 아직 부족하다. 가장 큰 이유는 full Zeeman이
없다는 사실 자체보다 먼저, `ls=0.001`이라는 실질적인 sensitivity calibration factor가
숨겨져 있고, ideal detector balance와 low-contrast/edge 조건이 hero 숫자에 반영되지 않으며,
UI extra-view cache가 사용자의 detector/IF 설정을 실제 계산에서 누락하기 때문이다.

권장 순서는 다음과 같다.

1. **P1:** extra-view/cache dependency와 live metadata를 바로잡고 회귀 test를 추가한다.
2. **P1:** `ls=0.001`을 effective opacity/participation calibration으로 공개하고 uncertainty를
   감도에 전달한다.
3. **P1:** contrast/transmitted-power/edge/resolution status로 hero를 보호한다.
4. **P2:** fixed reference arm/balance-point 모델과 AT/optimum interpolation을 넣는다.
5. **속도:** finite-IF Liouvillian affine assembly와 atomic-phasor cache를 적용한다.
6. 그 다음에 measured trace와 detector/SAM calibration이 확보되면 low-order Zeeman proxy를
   검증하고, full Zeeman/time-domain/spatial GROUP C의 필요성을 다시 판단한다.

이 순서의 1–5는 solver 차원을 늘리지 않거나 behavior를 보존하는 저비용 작업이다. full
Zeeman이나 full lock-in보다 먼저 처리해야, 현재 구현된 좋은 물리가 실험가에게 과도한
정밀도나 잘못된 calibration으로 전달되지 않는다.
