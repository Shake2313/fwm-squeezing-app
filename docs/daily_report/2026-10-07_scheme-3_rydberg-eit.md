# GABES 일일 Scheme 물리 검토 — 2026-10-07

## 0. 선택·현재 순서

Asia/Seoul 현지 일자 7일.

`n = (7 mod 5) + 1 = 3`.

실제 드롭다운 순서(`gabes/schemes/__init__.py:19-24`):

1. OD / SAS — `SASScheme`
2. Lambda coherence (EIT / AT / CPT) — `LambdaScheme`
3. **Rydberg-EIT electrometry — `RydbergEITScheme`**
4. Hanle / EIA / NMOR — `MagnetoScheme`
5. FWM — `FWMScheme`

`gabes/schemes/__init__.py:27-39`의 Lambda/Magneto 단일-regime alias는
드롭다운 항목 아님.

현재 `main@fcb312d`. 대규모 기존 staged/unstaged/untracked 작업 있음. 보존함.
Rydberg 생산 코드·전용 테스트는 HEAD 대비 변경 없음. 핵심 파일 마지막 커밋은
`d98bdf4`(2026-09-05). 현재 `gabes/schemes/rydberg.py` SHA-256은
`8E658687771478839F4951A712C2AA4EDC4E1337074841A32515B654397BA877`.

## 1. 먼저 찾은 기존 제안·TODO·issue

먼저 다음 읽음.

- `CLAUDE.md`, `README.md` 전체.
- 기존 Scheme 3 보고서 16개(2026-06-22…2026-09-27)를 검색·대조하고 직전 보고서를
  상세 검토.
- `docs/checklist.json`의 Rydberg 완료·활성·parked 항목.
- `gabes/schemes/rydberg.py`, electrometry/experiment/CSV helper.
- 전용 테스트, cell-heating workflow 문서·예제, User Guide, 실험 맥락 문서.
- 파일명 `TODO/issue/backlog/roadmap`, 코드 `TODO/FIXME/XXX/HACK`.
- GitHub 공개 open issue API.

검색 결과:

- Rydberg 코드·테스트·분석 경로 inline `TODO/FIXME/XXX/HACK`: 없음.
- 별도 로컬 Rydberg TODO/issue-note 파일: 없음.
- 공개 GitHub open issue: 0건.
- **기존 개선안 많음. 새 물리 로드맵을 임의 추가할 상황 아님.** 활성 핵심은
  `docs/checklist.json:3145-3234`, `3566-3581`, `3749-3829`.

기존 항목:

- `rydberg-extra-view-dependency-correctness` — P1, `ready`.
- `rydberg-effective-opacity-and-sweep-trust` — P1, `ready`.
- `rydberg-readout-validity-guardrails` — P1, `needs_decision`.
- `rydberg-experimental-validation-dataset` — P2, `blocked_external`.
- `rydberg-doppler-on-pole-residue` — P3, `ready`.
- `rydberg-zeeman-polarization-reference-mode` — P2, `parked`.
- `rydberg-readout-transfer-function-validation` — P2, `blocked_external`.
- `rydberg-rf-field-map-import` — P2, `blocked_external`.

추가로 기존 실험 맥락 문서 자체가 status 개선을 이미 제안함.
`references/rydberg_ju_experiment_context.md:45-46`은 같은 장치 논문에 나온 값들을
`ASSUMED`가 아니라 `REFERENCE`로 분류하라 명시. 그러나 예제의 electrometry block은
문헌값과 이상적 detector 가정을 한 `ASSUMED` 상태로 묶음
(`analysis/rydberg_cell_heating/example_config.json:89-102`). blanket 승격은 부정확.
입력별 분리 필요:

- `REFERENCE`: 6 µW, 30 mW, 0.15 mm, 50 mm, 37 GHz, 약 95% QE처럼 논문에 나온
  동일 장치 값.
- `ASSUMED`: 0.5 path split, 완전 balance, RIN correlation=1, electronics=0,
  participation/overlap=1처럼 측정 없는 이상화.

이는 metadata/status 정리. 계산 비용 O(1), spectrum·감도 값 불변.

## 2. 연구자 판정

**판정: 실제 유용한 compact 원자물리 구현. 반정량적 실험 설계, convention 확인,
내부 모델 교차검증, 교육용 코드 참조에 적합. 절대 감도·온도 최적점·Zeeman/편광·
공간 RF장·논문 수치의 독립 기준으로는 아직 부적합.**

유용 범위:

- ⁸⁵Rb `5S₁/₂ F=3 → 5P₃/₂ F′=4 → 40D₅/₂ → 39F₇/₂` compact
  4준위 cascade EIT/AT.
- probe/coupling power, beam diameter, transit broadening, optical depth, RF dressing,
  finite IF, detector-noise trade-off 탐색.
- AT split↔field, angular/RMS/peak convention, PSN/RIN/electronics budget 검산.
- heater/effective/cold-spot 온도 분리, axial density·Beer–Lambert, SAM uncertainty,
  CSV unit·hash provenance를 포함한 실험 분석 골격.

부적합 범위:

- 현 상태 수치를 publication-grade 절대 sensitivity/temperature optimum으로 인용.
- full Zeeman/polarization, optical pumping, stray field, ionization/charge dynamics.
- horn/cell standing-wave RF map, nonlinear SIG, full time-domain LO+SIG/lock-in.
- held-out 측정 없이 1.6 MHz anchor 재현을 독립 예측 성공으로 해석.

## 3. 실제 구현된 물리

### 3.1 원자계·정적 spectrum

- topology·한계: `gabes/schemes/rydberg.py:1-15`, `73-91`, `417-437`.
  5P decay, 작은 Rydberg decay, 5S–40D와 40D–39F dephasing 분리.
- probe/coupling Rabi는 `√P/d` anchored scaling
  (`gabes/schemes/rydberg.py:439-461`). 절대 optical dipole 예측 아님.
- `TRANSIT_FACTOR=0.6`, `PROBE_RABI_REF_MHZ=2.0`은 기준 1.6 MHz 폭에 맞춘 anchor
  (`gabes/schemes/rydberg.py:41-55`). 같은 폭 재현은 validation 아님.
- heater/effective/cold-spot 분리와 cold-spot pressure 기반 local density:
  `gabes/schemes/rydberg.py:463-505`.
- detuning-affine steady-state scan: `gabes/schemes/rydberg.py:507-527`, `558-619`.
  기본 Doppler-off는 한 static class. opt-in Doppler-on은 2 m/s, 4σ Maxwell grid와
  per-level residual `k` ratio.
- residual Zeeman은 full manifold가 아니라 susceptibility Gaussian convolution
  (`gabes/schemes/rydberg.py:530-542`, `621-627`).
- Beer–Lambert는 표준. 그러나 `ls=0.001`은 숨은 effective opacity/participation
  (`gabes/schemes/rydberg.py:629-643`). 공개 atomic line strength 아님.

### 3.2 finite-IF electrometry·detector

- LO-dressed steady state 주위 1차 weak-SIG solve:
  `gabes/rydberg_electrometry.py:116-194`.
- minus sideband 한 번 solve. plus sideband는 Hermitian conjugate 복원
  (`gabes/rydberg_electrometry.py:185`). 물리 보존 최적화.
- Doppler phasor는 complex coherent average
  (`gabes/rydberg_electrometry.py:211-227`, `gabes/schemes/rydberg.py:773-819`).
- RF dipole, angular factor, peak/RMS, detuned split helper:
  `gabes/rydberg_electrometry.py:231-316`.
- 두 photodiode Schottky noise, correlated RIN, electronics와 field-ASD 변환:
  `gabes/rydberg_electrometry.py:358-514`.
- 문헌 11.2/12.5 nV/cm/√Hz를 결과에 강제 주입하지 않음
  (`tests/test_rydberg_eit.py:80-113`). 좋은 신뢰경계.
- 현재 reference arm power는 각 scan point의 signal power×ratio로 재설정
  (`gabes/schemes/rydberg.py:894-944`). 실제 고정 reference arm보다 이상적인
  pointwise-balance envelope.

### 3.3 실험 helper·예제

- 온도 상태: `gabes/rydberg_experiment.py:80-168`.
- axial density: `gabes/rydberg_experiment.py:191-294`.
- segmented Beer–Lambert: `gabes/rydberg_experiment.py:307-343`.
- Gaussian overlap/effective atom 수: `gabes/rydberg_experiment.py:346-449`.
- SAM field·불확도·far-field guard: `gabes/rydberg_experiment.py:485-595`.
- temperature/density dephasing fit과 identifiability:
  `gabes/rydberg_experiment.py:661-803`.
- EIT/RF/PSD CSV unit와 SHA-256 provenance:
  `gabes/rydberg_experimental_csv.py:35-100`, `441-691`.
- workflow 증거 상태 `MEASURED/PPT/REFERENCE/FITTED/PREDICTED/ASSUMED/PENDING`:
  `analysis/rydberg_cell_heating/README.md:22-39`.
- 예제 핵심 EIT/RF/temperature 입력은 모두 `PENDING`
  (`analysis/rydberg_cell_heating/example_config.json:55-77`). validation dataset 아님.

## 4. 오늘의 수치 재검증

동일 worktree. JIT warm-up 뒤 median. headless readout.

| 조건 | 결과 | 해석 |
|---|---:|---|
| 기본 EIT compute | 4.204 ms | compact 4-level 정적 solve |
| 기본 EIT `T(0)` / 대비 / FWHM | 0.940016615 / 0.135391492 / 1.613058583 MHz | anchor 규모 재현 |
| 기본 EIT samples/FWHM | 71.69 | 기본 폭의 grid 분해능 충분 |
| 기본 AT compute | 4.013 ms | 정적 spectrum만 |
| 기본 AT `T(0)` / splitting | 0.806029283 / 3.496500 MHz | 3.7 MHz LO와 같은 규모. line pulling 있음 |
| 기본 총·PSN ASD | 6.082387 / 6.082387 nV/cm/√Hz | 조건부 예측. 문헌 11.2/12.5보다 낙관적 |
| 기본 optimum detuning | −2.220 MHz | 801-point discrete optimum |
| finite-IF cold / phasor-cached readout | 78.34 / 21.71 ms | 조립·detector loop 비용 큼 |
| `ls=0.0005/0.001/0.002` 총 ASD | 11.781990 / 6.082387 / 3.238776 | 숨은 opacity가 절대 감도 지배 |
| 80 °C EIT `T(0)` / 대비 / FWHM | 1.798923×10⁻⁶ / 1.798923×10⁻⁶ / 0.307015 MHz | 사실상 암흑인데 numeric linewidth 유지 |
| LO=0.1 MHz | AT unresolved, 총 ASD 10.4950 | 미분해 feature에도 numeric sensitivity hero 유지 |
| Doppler-on compute / 총 ASD | 156.90 ms / 1007.281 nV/cm/√Hz | off보다 약 37배 느림. calibration regime 격차 큼 |
| 기본 cell-heating view | 758.78 ms, 원자 compute 20회 | 최적 47 °C / 0.971731. `PREDICTED` |
| detector-only 변경 뒤 temperature view | 305.26 ms, 원자 compute 20회 | cache 계층 미분리 증거 |

### 4.1 기본 grid 수렴

같은 ±11.1 MHz 범위:

| 점 수 | 간격 | 총 ASD | optimum | AT split |
|---:|---:|---:|---:|---:|
| 801 | 27.750 kHz | 6.082387 | −2.220000 MHz | 3.496500 MHz |
| 1601 | 13.875 kHz | 6.078733 | −2.206125 MHz | 3.496500 MHz |
| 3201 | 6.9375 kHz | 6.078733 | −2.206125 MHz | 3.496500 MHz |

기본 총 ASD의 801→수렴 변화 약 0.0601%. 기본 운전점 numerical grid 충분.
그러나 low-contrast, near-zero power, scan-edge, unresolved feature의 실험 관측성은
grid 수렴과 다른 문제.

### 4.2 온도 최적점 비교

저장소 실험 맥락은 논문 temperature-vs-PSN model optimum 약 33 °C를 명시
(`references/rydberg_ju_experiment_context.md:20-22`). 현재 내장 sweep은 47 °C.
측정 첫 scan도 38 °C에서 좋아졌다가 42.5/44 °C에서 나빠졌지만 측정 순서와 drift가
분리되지 않음(`references/rydberg_ju_experiment_context.md:108-126`).

따라서 47 °C를 운전점 권고로 쓰면 안 됨.

## 5. 실험-물리 참조 한계

1. **anchor와 prediction 혼동.** transit/probe/coupling anchor가 기준 폭에 맞음.
2. **숨은 opacity.** `ls=0.001`이 transmission, photon shot noise, sensitivity를 지배.
3. **관측가능성 gate 부족.** `_eit_features`는 양의 대비와 half-height crossing만 확인
   (`gabes/schemes/rydberg.py:979-1003`). 80 °C 반례의 대비 10⁻⁶도 hero 폭 유지.
4. **미분해 AT에서도 감도 hero.** 테스트가 AT status와 동시에 총 감도 hero를 고정
   (`tests/test_rydberg_eit.py:134-141`).
5. **receiver 이상화.** reference arm이 scan point마다 signal과 함께 변함. 실험은
   고정 reference power 또는 한 balance point 필요.
6. **Doppler-off/on은 precision tier 아님.** calibration 의미가 다른 경로. on은
   1007 nV/cm/√Hz. on finite-IF는 full scan 대신 static-slope 한 점만 계산.
7. **누락 물리 큼.** full Zeeman/polarization, spatial RF standing wave, measured transfer
   function, nonlinear/time-domain readout 없음.
8. **held-out 데이터 없음.** example의 raw EIT/RF/temperature 입력은 `PENDING`.

## 6. 기존 개선안: 비용·물리 보존 평가

| 기존 항목 | 비용 | 물리·행동 보존 판단 |
|---|---|---|
| `rydberg-effective-opacity-and-sweep-trust` | named input·provenance 표시는 O(1), 사실상 0. uncertainty sweep만 opt-in 선형 비용. | 기본값 0.001 고정 시 기존 배열·readout `≤1e-12` 보존 가능. line strength와 participation 분리 필요. **최우선.** |
| `rydberg-readout-validity-guardrails` | contrast, detector power, samples/FWHM, edge margin은 기존 801 배열의 O(N) 후처리. 소수 local solve refinement만 작은 추가 비용. | 원자 spectrum 불변. hero 공개 여부만 gate 가능. fixed-reference readout은 정의 변경이므로 현재 pointwise envelope도 별도 diagnostic으로 보존. |
| `rydberg-extra-view-dependency-correctness` | dependency key와 phasor LRU는 일부 구현(`gabes/schemes/rydberg.py:821-865`). 그러나 detector-only 변경도 20 atomic computes, 305 ms. | atomic/phaser/Beer–Lambert/detector 네 cache 분리 시 수치 불변. checklist `ready`는 실제로 partial/in-progress 성격. |
| `rydberg-experimental-validation-dataset` | 계산보다 측정·metadata·uncertainty·hash 비용 지배. | simulation으로 대체 불가. 완료 전 absolute/optimum은 `PREDICTED` 또는 `UNSUPPORTED`. |
| `rydberg-doppler-on-pole-residue` | 현재 on 156.9 ms, off 4.2 ms. row당 작은 eigendecomposition으로 큰 가속 가능. | 현재 2 m/s·4σ **이산 측도**와 coherent phasor를 유지하고 per-row guard 통과하면 behavior-preserving. default off라 P3 타당. |
| `rydberg-zeeman-polarization-reference-mode` | Hilbert/Liouville 차원, B·velocity grid, dipole provenance 필요. 무시 가능 overhead 아님. | 특정 state/polarization/B dataset 뒤 opt-in reference로만 타당. |
| `rydberg-readout-transfer-function-validation` | 측정 complex `H(IF)` 곱은 O(N_IF), 계산 비용 작음. | 원자 응답 보존. 측정 gain/phase 없어서 외부 blocked. full time-domain보다 먼저 할 일. |
| `rydberg-rf-field-map-import` | O(N_space×N_scan), mesh convergence 필요. 무시 가능 overhead 아님. | provenance-tagged external EM/측정 map이 내부 3-D EM보다 타당. uniform-map reduction gate 필요. |
| 문헌값/가정값 status 분리 | 입력별 metadata O(1). | 계산값 완전 불변. composite block 전체를 `REFERENCE`로 승격하지 말고 published/assumed 입력 분리. |

## 7. 순수 코드 최적화

### 7.1 BLAS controller 지연 캐시 — 가장 넓은 저위험 이득

현재 `core.blas_single_thread()`는 매 solve마다 `threadpool_limits(...)` context를 새로 만듦
(`gabes/core.py:23-34`). lazy `ThreadpoolController`를 한 번 만들고 limiter를 재사용한
prototype:

- 기본 EIT median: **3.7466 → 0.8530 ms**.
- **4.39×**.
- `chi_bar np.array_equal=True`.

원자물리·solver 불변. NumPy/BLAS load 뒤 lazy initialization, 병렬 호출 중첩, process-wide
thread-limit 복원을 회귀로 확인해야 함. Lambda 등 다른 작은-batch scheme에도 이득 가능.

### 7.2 finite-IF Liouvillian affine 조립

현재 Doppler-off finite-IF는 801 detuning마다 Python에서 Hamiltonian/Liouvillian 재조립
(`gabes/schemes/rydberg.py:807-815`). detuning 의존은 affine.

- list 조립: **34.718 ms**.
- `L(scan)=L(0)+scan·[L(1)-L(0)]` broadcast: **2.479 ms**.
- **14.0×**, `np.array_equal=True`, max abs difference 0.

행렬 조립만 변경. solve·scan·물리·부동소수점 결과까지 동일.

### 7.3 detector/noise scan 벡터화

현재 801점마다 channel/detector/noise/sensitivity dataclass 객체 생성
(`gabes/schemes/rydberg.py:910-946`). 식은 elementwise 배열화 가능.

- scalar loop: **17.256 ms**.
- vector 식: **0.0316 ms**.
- 약 **546×**.
- max relative difference `3.10×10⁻¹⁶`.

표시·물리 행동 불변. bitwise는 아님. machine-roundoff gate와 scalar fixture 필요.

### 7.4 extra-view 네 계층 cache

분리 대상:

1. atomic steady-state spectrum.
2. finite-IF atomic phasor.
3. cell/opacity Beer–Lambert.
4. detector/RF-dipole/SAM conversion.

현재 2번 일부만 명시적 LRU. detector electronics만 바꿔도 temperature view가 원자
`compute` 20회 호출. 1·3·4 경계 분리하면 결과 정의 없이 큰 속도 향상.

### 7.5 작은 정리

- temperature sweep이 표시 문자열 `"Max spectral slope"`, `"IF discriminator"`를
  split해 float로 재파싱(`gabes/schemes/rydberg.py:1412-1423`). `_readout` numeric field
  직접 사용. 포맷-계산 결합 제거.
- EIT temperature loop는 `_readout` 뒤 `_transmission`을 다시 호출. 기존 readout 배열
  재사용 가능. 이득 작지만 행동 불변.

## 8. 문서·테스트·예제 적합성

좋음:

- `README.md:117-137`: finite-IF, 온도 분리, SAM, deferred 범위를 정확히 제한.
- `analysis/rydberg_cell_heating/README.md:22-39`, `72-117`: evidence status,
  provenance, lumped axial 한계, conditional electrometry 명시.
- finite-IF sideband residual·Hermiticity·한 번 solve:
  `tests/test_rydberg_electrometry.py:39-94`.
- low-IF→static derivative: `tests/test_rydberg_eit.py:309-337`.
- detector-only phasor reuse: `tests/test_rydberg_eit.py:355-427`.
- temperature/cold spot, axial/SAM/dephasing identifiability:
  `tests/test_rydberg_eit.py:456-517`, `tests/test_rydberg_experiment.py:19-200`.
- CSV unit·duplicate·hash·음수 PSD 거부:
  `tests/test_rydberg_experimental_csv.py:14-109`.

부족:

- User Guide는 정적 EIT/AT만 기술
  (`docs/Userguide/GABES_User_Guide_v2.html:535-538`, `617-625`). finite-IF 감도,
  conditional calibration, temperature/provenance, hidden opacity 한계 누락.
- `tests/test_rydberg_eit.py:1-5`의 “public UI shows the static spectrum only”는 현재 UI와
  불일치.
- example의 paper-derived 값과 ideal detector 가정이 composite `ASSUMED` 아래 섞임.
- 80 °C low-transmission, fixed-reference, edge optimum, local refinement acceptance fixture 없음.
- opacity provenance와 Doppler mode별 calibration scope를 직접 고정하는 테스트 없음.

## 9. 검증

- 전용 Scheme/helper/workflow/adapter:
  `69 passed in 93.03 s`.
- Rydberg + affine kernel + scheme render 보강:
  `43 passed, 26 deselected in 9.16 s`.
- 필수 `python -m pytest -q`:
  `2140 passed, 4 skipped, 2 failed in 2644.39 s`.
  두 실패 모두 Rydberg 아님.
  - `tests/quantum/test_thermal_campaign_extend.py`의 fresh-capsule extraction이 Windows
    임시 `.git/grand-challenge-runs/...` 중첩 경로에서
    `ballistic_linear_channel.py` destination을 만들지 못해 `FileNotFoundError`.
  - `tests/test_fwm_excess_noise.py`의 Streamlit ScrubField-fallback rerun에서
    gain-proportional noise mode 전환 뒤 metric이 이전 constant 값 2.0을 유지
    (기대 0.04).

전체 회귀 통과 아님. 그러나 전용 Rydberg 69개와 관련 kernel/render 43개는 통과.

## 10. 권장 순서

1. **P1 opacity provenance 공개.** runtime 0, 절대값 신뢰도 효과 최대.
2. **P1 observability/fixed-reference gate.** 80 °C·weak-LO 반례 fixture부터.
3. **P1 순수 성능.** lazy cached BLAS controller → bit-identical affine 조립 → detector
   vectorization → four-layer cache.
4. **P1 provenance/status·문서 정리.** published input과 ideal assumption 분리. User Guide,
   stale test docstring, checklist partial 상태 정정.
5. **P2 held-out dataset.** EIT/RF/PSD/temperature, fixed reference, transfer function.
6. **P3 Doppler pole.** Doppler-on이 실제 sweep 병목이 된 뒤 현재 이산 측도 보존 이식.
7. Zeeman/편광·RF map·time-domain은 특정 측정·runtime budget 전까지 opt-in/blocked 유지.

핵심: 새 거대 물리 엔진보다 **현재 compact model의 opacity provenance, 관측 가능성,
receiver 조건, 입력별 evidence status, held-out 검증**이 먼저.
