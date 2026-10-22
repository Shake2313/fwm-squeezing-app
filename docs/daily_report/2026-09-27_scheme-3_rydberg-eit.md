# GABES 일일 Scheme 물리 검토 — 2026-09-27

## 0. 선택·현재 순서

Asia/Seoul 현지 일자 27일. `n = (27 mod 5) + 1 = 3`.

현재 드롭다운 순서(`gabes/schemes/__init__.py:19-25`):

1. OD / SAS — `SASScheme`
2. Lambda coherence (EIT / AT / CPT) — `LambdaScheme`
3. **Rydberg-EIT electrometry — `RydbergEITScheme`**
4. Hanle / EIA / NMOR — `MagnetoScheme`
5. FWM — `FWMScheme`

`gabes/schemes/__init__.py:27-39`의 Lambda/Magneto 단일-regime alias는 드롭다운
항목 아님.

현재 HEAD `5b29c4f`. Rydberg 생산 코드, 전용 테스트, cell-heating workflow는
2026-09-22 검토 뒤 커밋·작업 트리 변화 없음. 저장소에는 다른 대규모 기존 dirty/staged
작업 있음. 보존함.

## 1. 선행 제안 검색 결과

먼저 `CLAUDE.md`, `README.md`, Scheme 3 일일 보고서 15개, `docs/checklist.json`,
Rydberg 소스·테스트·분석 workflow, 파일명 기반 issue/TODO/backlog/roadmap 검색 수행.

- Rydberg 코드·테스트의 `TODO/FIXME/XXX/HACK`: 없음.
- 로컬 `issue`, `todo`, `backlog`, `roadmap` 이름 파일: 없음.
- **기존 개선안 있음. 새 제안을 처음부터 만들 상황 아님.** 주요 항목:
  `docs/checklist.json:2936-3018`, `3357-3371`, `3540-3619`.
- 직전 상세 비용 검토: `docs/daily_report/2026-09-17_scheme-3_rydberg-eit.md:225-305`.
- 직전 Scheme 3 검토: `docs/daily_report/2026-09-22_scheme-3_rydberg-eit.md:42-65`.

상태 추적 주의: `rydberg-extra-view-dependency-correctness`는 checklist에서 여전히
`ready`. 그러나 2026-09-05 코드가 non-recompute dependency key와 finite-IF phasor
LRU 일부 구현함(`gabes/schemes/rydberg.py:821-865`, `1528-1556`;
`tests/test_rydberg_eit.py:45-56`, `355-427`). 완전 미구현 아님. 반면 detector-only
변경에도 온도 extra view가 20회 원자 `compute`를 다시 호출함. 원자/phasor/
Beer–Lambert/detector 네 계층 cache 완료 기준은 미충족. **부분 구현 + checklist 상태
갱신 지연**으로 판단.

## 2. 연구자 판정

**판정: 실제 원자물리 있음. 반정량적 실험 설계·내부 모델 교차검증·교육 참조로
유용. 절대 감도, 온도 최적점, 편광/Zeeman 응답, 공간 RF장을 위한 독립 실험 기준으로는
아직 부적합.**

유용 범위:

- ⁸⁵Rb `5S₁/₂ F=3 → 5P₃/₂ F′=4 → 40D₅/₂ → 39F₇/₂` 4준위 cascade
  EIT/AT 선형·비선형 경향.
- probe/coupling power, beam size, transit broadening, optical depth, RF dressing,
  finite IF, detector noise 사이 trade-off 탐색.
- AT splitting↔field convention, RMS/peak, RF dipole, balanced-detector PSN/RIN/
  electronics bookkeeping 확인.
- heater setpoint/effective vapor/cold spot 구분, axial density·Beer–Lambert,
  SAM uncertainty, CSV 단위·provenance를 포함한 실험 분석 골격.

부적합 범위:

- 논문용 절대 sensitivity/temperature optimum 인용.
- full Zeeman/polarization, optical pumping, stray-field, ionization/charge dynamics.
- horn/cell standing-wave RF map, full time-domain LO+SIG/lock-in transfer.
- held-out 측정 없이 현재 calibration과 같은 1.6 MHz EIT 폭을 독립 예측 성공으로 해석.

## 3. 구현된 물리

### 3.1 원자계·스펙트럼

- 4준위 topology와 범위: `gabes/schemes/rydberg.py:1-15`, `72-91`, `417-437`.
  5P decay, 작은 Rydberg decay, 5S–40D와 40D–39F dephasing 분리.
- probe/coupling Rabi는 `√P/d` scaling. 절대 dipole 예측 아님. 기준 운전점 anchor임
  (`:41-55`, `439-461`). 코드 자체가 bare D2 probe dipole이 lumped model에서 약
  100배 과대라고 명시.
- heater/effective/cold-spot 분리와 cold-spot pressure 기반 local density:
  `:463-505`, `558-590`. 실험 셀 해석에 유용한 구분.
- detuning마다 정상상태 Liouvillian을 풀어 probe coherence를 얻음. 정적 spectrum은
  이미 affine kernel 사용(`:507-527`). 801점 scan은 `:544-556`; Hamiltonian은
  `:600-608`.
- Doppler-off는 정적 한 class. Doppler-on은 2 m/s, 4σ 절단 Maxwell grid와 level별
  residual `k` ratio 사용(`:610-619`). 두 모드는 단순 정밀도 tier가 아니라 calibration
  의미가 다른 물리 경로.
- residual Zeeman은 full magnetic manifold가 아니라 χ의 Gaussian convolution
  (`:621-627`).
- Beer–Lambert 자체는 표준. 그러나 `ls=0.001`은 숨은 effective opacity/participation
  (`:629-643`, `672-681`). 공개 line strength 아님.

### 3.2 finite-IF electrometry·검출

- LO-dressed steady state 주위 weak-SIG 선형 응답. Hermitian perturbation의 minus
  sideband 하나를 풀고 plus sideband를 conjugate transpose로 복원
  (`gabes/rydberg_electrometry.py:116-194`). 물리 보존 최적화.
- Doppler-on phasor는 복소 coherent average
  (`gabes/schemes/rydberg.py:787-806`; `gabes/rydberg_electrometry.py:211-227`).
- RF dipole, angular/polarization factor, peak/RMS convention 분리
  (`gabes/rydberg_electrometry.py:230-316`). AT split helper도 detuning 제거와 line-pulling
  한계를 명시.
- balanced detector는 두 arm Schottky noise, correlated RIN, electronics를 분리
  (`:319-438`). field ASD 변환은 absolute detector slope가 있어야 가능
  (`:441-514`). 문헌의 11.2/12.5 숫자를 결과에 주입하지 않음.
- 현재 scan은 각 detuning에서 reference power를 signal power×ratio로 다시 잡음
  (`gabes/schemes/rydberg.py:902-927`). 실험의 고정 reference arm이 아니라 이상적
  pointwise-balance envelope.

### 3.3 실험 helper·예제

- 온도 상태와 axial profile: `gabes/rydberg_experiment.py:53-294`.
- segmented Beer–Lambert: `:297-343`.
- Gaussian overlap/effective atom number: `:346-448`.
- SAM field·불확도·far-field guard: `:452-594`.
- temperature/density dephasing fit과 identifiability: `:597-803`.
- EIT/RF/PSD CSV의 단위, duplicate, SHA-256 provenance:
  `gabes/rydberg_experimental_csv.py:35-100`, `441-691`.
- workflow 증거 상태 `MEASURED/PPT/REFERENCE/FITTED/PREDICTED/ASSUMED/PENDING`:
  `analysis/rydberg_cell_heating/README.md:22-39`.
- 예제 입력의 EIT/RF/temperature는 모두 `PENDING`
  (`analysis/rydberg_cell_heating/example_config.json:53-78`). RF dipole은 `REFERENCE`,
  electrometry 설정은 `ASSUMED` (`:88-102`). 검증 dataset 아님.

## 4. 오늘의 수치 재검증

동일 worktree. JIT warm-up 뒤 median. headless readout 사용.

| 조건 | 결과 | 해석 |
|---|---:|---|
| 기본 EIT compute | 3.915 ms | 빠른 compact 4-level solve |
| 기본 EIT `T(0)` / 대비 / FWHM | 0.940016615 / 0.135391492 / 1.613058583 MHz | 1.6 MHz anchor 규모 재현 |
| 기본 EIT samples/FWHM | 71.69 | 기본 폭은 scan 분해능 충분 |
| 기본 AT compute | 3.802 ms | 정적 spectrum만 |
| 기본 AT `T(0)` / splitting | 0.806029283 / 3.4965 MHz | 3.7 MHz LO와 같은 규모; line pulling 있음 |
| 기본 finite-IF 총·PSN ASD | 6.082387 / 6.082387 nV/cm/√Hz | 조건부 계산. 문헌 11.2/12.5보다 낙관적 |
| 기본 optimum detuning | −2.220 MHz | 801점 grid 결과 |
| finite-IF cold / cached headless | 73.91 / 22.35 ms | phasor solve와 scalar detector loop 비용 |
| `ls=0.0005/0.001/0.002` 총 ASD | 11.781990 / 6.082387 / 3.238776 nV/cm/√Hz | 숨은 opacity에 강한 의존 |
| 80 °C EIT `T(0)` / 대비 / FWHM | 1.798923×10⁻⁶ / 1.798923×10⁻⁶ / 0.307015 MHz | 거의 암흑인데 numeric width hero 유지 |
| 80 °C samples/FWHM | 13.65 | 수치 분해는 됨. **실험 관측성은 아님** |
| LO=0.1 MHz | AT unresolved, 총 ASD 10.50 nV/cm/√Hz hero | 미분해 feature에도 precise sensitivity 노출 |
| Doppler-on compute / 총 ASD | 159.40 ms / 1007.281 nV/cm/√Hz | off보다 40.7배 느림; calibration regime 격차 큼 |
| 기본 cell-heating extra view | 842.66 ms, 20 atomic computes | 최적 47 °C, 0.971731 nV/cm/√Hz. `PREDICTED` |

문헌 맥락 `references/rydberg_ju_experiment_context.md:20-22`, `94-100`, `140-155`:
PSN 11.2(4), 측정 12.5(8) nV/cm/√Hz, 논문 온도 모델 optimum 약 33 °C.
GABES 내장 sweep 47 °C와 불일치. 운전점 권고 금지.

### 기본 grid 수렴 확인

동일 범위 global refinement:

| 점 수 | 간격 | 총 ASD | optimum | AT split |
|---:|---:|---:|---:|---:|
| 801 | 27.75 kHz | 6.082387 | −2.220000 MHz | 3.496500 MHz |
| 1601 | 13.875 kHz | 6.078733 | −2.206125 MHz | 3.496500 MHz |
| 3201 | 6.9375 kHz | 6.078733 | −2.206125 MHz | 3.496500 MHz |

기본점 총 ASD의 801→수렴 변화 약 0.060%. 기본점 numerical resolution 양호.
그러나 저대비·scan-edge·강한 detuning 조건의 관측성 보장은 아님. 전 범위 refinement는
83→147→347 ms. 실제 guard는 optimum 주변 국소 몇 점 refinement가 더 적절.

## 5. 실험-물리 참조 한계

1. **anchor와 예측 혼동 위험.** transit factor와 probe/coupling Rabi anchor가 기준
   1.6 MHz에 맞춰짐(`gabes/schemes/rydberg.py:41-55`). 같은 폭 재현은 독립 검증 아님.
2. **숨은 opacity.** `ls=0.001`이 절대 transmission·shot noise·sensitivity 지배.
   provenance·불확도 없음.
3. **관측가능성 gate 부족.** `_eit_features`는 양의 대비와 crossing만 확인
   (`:979-1003`). 80 °C 반례처럼 대비 10⁻⁶에도 폭 hero 유지. weak LO 미분해 AT에도
   sensitivity hero 유지(`:1078-1082`, `1127-1207`; test가 현재 동작 고정:
   `tests/test_rydberg_eit.py:134-141`).
4. **receiver가 지나치게 이상적.** reference arm이 scan point마다 완전 balance. 고정 arm,
   balance point, residual RIN을 별도 모델링해야 실험 readout 참조 가능.
5. **Doppler-off/on 동등 calibration 아님.** on에서 1007 nV/cm/√Hz. finite-IF on은
   full detuning scan이 아니라 static-slope 한 점만 평가(`gabes/schemes/rydberg.py:787-819`).
6. **누락 물리 큼.** full Zeeman/polarization, optical pumping, spatial RF standing wave,
   measured transfer function, nonlinear/time-domain lock-in 없음(`:1-15`, `1559-1585`).
7. **held-out 실험 데이터 없음.** workflow 구조는 좋지만 예제 핵심 측정은 PENDING.

## 6. 기존 개선안: 비용·물리 보존

| 기존 항목 | 현재 상태·계산 비용 | 물리 보존 판단 |
|---|---|---|
| `rydberg-effective-opacity-and-sweep-trust` P1 ready | `ls`를 named input으로 노출하는 기본 비용 0. uncertainty sweep만 선택적 선형 비용. | 기본값 0.001 고정 시 기존 spectrum/readout `≤1e-12` 보존 가능. line strength와 participation 중복 금지. **최우선.** |
| `rydberg-readout-validity-guardrails` P1 needs_decision | contrast, detector power, samples/FWHM, edge margin은 기존 801점의 O(N) 후처리. 국소 refinement만 소수 추가 solve. | 원자 물리 불변. threshold 정책 필요. fixed reference arm은 receiver 정의를 현실화하며 기존 pointwise envelope를 별도 diagnostic으로 보존 가능. |
| `rydberg-extra-view-dependency-correctness` P1 ready | dependency key·phasor LRU는 일부 완료. detector electronics 5 pA 변경만으로도 temperature view가 20 atomic computes 재실행; 313.66 ms. | 네 계층 cache 분리 시 결과 불변. detector-only knob는 atomic steady state 재사용해야 함. checklist를 `in_progress` 성격으로 정정 필요. |
| `rydberg-experimental-validation-dataset` P2 blocked_external | 계산보다 측정·metadata·uncertainty·hash 비용 지배. | 내부 simulation으로 대체 불가. 완료 전 absolute/optimum은 `PREDICTED`/`UNSUPPORTED`. |
| `rydberg-doppler-on-pole-residue` P3 ready | 현재 on 159 ms, off 3.9 ms. row당 작은 eigendecomposition으로 가속 가능. | 현재 2 m/s·4σ **이산 측도**와 coherent phasor를 유지할 때 behavior-preserving. continuous Faddeeva는 별도 물리 변경. default off라 P3 타당. |
| `rydberg-zeeman-polarization-reference-mode` P2 parked | Hilbert/Liouville 차원, velocity/B grid, 상태·편광 provenance 필요. 무시 가능한 overhead 아님. | 특정 상태·편광·B dataset 뒤 opt-in reference로만 활성화 타당. |
| `rydberg-readout-transfer-function-validation` P2 blocked_external | 측정 complex `H(IF)`를 phasor에 곱하는 비용 O(N_IF), 사실상 작음. | 선형 원자 응답 보존. 측정 gain/phase 확보 뒤 full time-domain보다 먼저. |
| `rydberg-rf-field-map-import` P2 blocked_external | O(N_space×N_scan), mesh convergence 필요. 무시 가능한 overhead 아님. | provenance-tagged 외부 EM/측정 map import가 내부 3-D EM 구현보다 타당. uniform-map reduction gate 필요. |

## 7. 순수 코드 최적화

### 7.1 finite-IF Liouvillian affine 조립 — 우선

현재 Doppler-off finite-IF는 801 detuning마다 Python에서 Hamiltonian과 Liouvillian을
재조립(`gabes/schemes/rydberg.py:807-815`). detuning 의존성 affine.

- 현재 list 조립 median: **35.359 ms**.
- `L(scan)=L(0)+scan·[L(1)-L(0)]` broadcast: **2.889 ms**.
- **12.24×**, `np.array_equal=True`, max abs difference 0.

행렬 조립만 변경. solver, scan, 물리, 부동소수점 결과까지 동일. 가장 안전한 patch.

### 7.2 detector/noise scan 벡터화

현재 801점마다 dataclass 3개와 noise/sensitivity 객체 생성
(`gabes/schemes/rydberg.py:908-946`). `gabes/rydberg_electrometry.py:415-429`의
Schottky/RIN 식을 배열로 평가하고 최적점 하나만 객체화 가능.

- scalar loop: **17.541 ms**.
- vector 식: **0.0243 ms**, 약 **722×**.
- max relative difference `3.10×10⁻¹⁶`.

표시·물리 행동 불변. bitwise 동일은 아님; machine-roundoff gate와 scalar fixture 필요.

### 7.3 extra-view 네 계층 cache

1. atomic steady-state spectrum
2. finite-IF atomic phasor
3. cell/opacity Beer–Lambert
4. detector/RF dipole/SAM conversion

현재 2번만 명시적 LRU. 1번은 app cache이나 extra loop 내부 재사용 안 함. detector-only
변경이 20 atomic computes를 재실행. 계층 분리하면 속도 향상, 결과 정의 불변.

### 7.4 문자열 재파싱 제거

온도 sweep은 표시용 `"Max spectral slope"`, `"IF discriminator"` 문자열을 split해
float로 되돌림(`gabes/schemes/rydberg.py:1412-1423`). `_readout`에 numeric field를
직접 보관. 비용 작지만 포맷-계산 결합 제거. 동작 불변.

### 7.5 pole–residue 위치

Doppler-on 전용 병목 가속. default path 이득 없음. 먼저 affine 조립·detector 배열화·
cache layering 수행. pole port는 기존 이산 Maxwell 측도 보존 gate와 함께 P3 유지.

## 8. 문서·테스트·예제 적합성

좋음:

- `README.md:115-137`: finite-IF, 온도 분리, SAM, deferred 범위를 짧고 정확히 명시.
- `analysis/rydberg_cell_heating/README.md:22-39`, `72-117`: 증거 상태, provenance,
  lumped axial approximation, conditional electrometry를 명시.
- 저 IF→정적 LO 미분 일치: `tests/test_rydberg_eit.py:309-337`.
- sideband/Hermiticity/한 번 solve: `tests/test_rydberg_electrometry.py:39-94`.
- temperature/cold-spot 분리: `tests/test_rydberg_eit.py:456-475`.
- axial/SAM/dephasing identifiability: `tests/test_rydberg_experiment.py:19-200`.
- CSV 단위·hash·음수 PSD 거부: `tests/test_rydberg_experimental_csv.py:14-109`.

부족:

- User Guide는 정적 EIT/AT만 기술
  (`docs/Userguide/GABES_User_Guide_v2.html:535-538`, `617-625`). finite-IF 감도,
  conditional calibration, temperature/provenance, opacity 한계 누락.
- `tests/test_rydberg_eit.py:1-5`의 “public UI shows the static spectrum only”는 현재 UI와
  불일치.
- observability acceptance fixture 없음: 80 °C low-transmission, edge optimum,
  fixed-reference, local refinement를 회귀로 고정해야 함.
- opacity provenance와 Doppler mode별 calibration scope를 직접 고정하는 테스트 없음.

## 9. 검증 결과

- 전용 Scheme/helper/workflow/adapter:
  `69 passed in 13.67 s`.
- Rydberg 관련 affine kernel·scheme render 보강:
  `30 passed, 8 deselected in 7.54 s`.
- 필수 `python -m pytest -q`:
  24.36 s 후 collection 중단. `references/AutoOD/test_autood_gui_adjustments.py`가
  선택 의존성 `sympy`를 import하나 현재 환경에 없음. 생산 Rydberg failure 아님.
  전체 회귀 통과도 아님.

## 10. 권장 순서

1. **P1 opacity provenance 공개.** 비용 0. 절대값 신뢰도에 최대 효과.
2. **P1 observability/fixed-reference gate.** 80 °C·weak-LO 반례부터 fixture화.
3. **P1 성능 patch.** bit-identical affine 조립 → detector vectorization → four-layer cache.
4. **P1 상태·문서 정리.** extra-view 항목을 partial/in-progress로 갱신. User Guide와 stale
   test docstring 수정.
5. **P2 held-out 데이터.** EIT/RF/PSD/temperature와 고정 reference/transfer function.
6. **P3 Doppler pole.** Doppler-on이 실제 기본/sweep 병목이 된 뒤 이산 측도 보존 이식.
7. Zeeman/편광·RF map·time-domain은 특정 측정과 runtime budget 전까지 opt-in/blocked 유지.

핵심: 새 거대 물리 엔진보다 **현재 compact physics의 opacity provenance, 관측 가능성,
receiver 조건, held-out 검증**이 먼저임.
