# GABES 일일 Scheme 물리 검토 — 2026-10-12

## 0. 선택·현재 순서

Asia/Seoul 현지 일자 12일.

`n = (12 mod 5) + 1 = 3`.

실제 드롭다운 순서(`gabes/schemes/__init__.py:19-24`):

1. OD / SAS — `SASScheme`
2. Lambda coherence (EIT / AT / CPT) — `LambdaScheme`
3. **Rydberg-EIT electrometry — `RydbergEITScheme`**
4. Hanle / EIA / NMOR — `MagnetoScheme`
5. FWM — `FWMScheme`

`gabes/schemes/__init__.py:27-39`의 Lambda/Magneto 단일-regime alias는
드롭다운 항목 아님.

현재 `main@39be4ae`. 광범위한 기존 staged/unstaged/untracked 작업 있음. 보존함.
Rydberg 생산 코드·전용 테스트는 worktree 변경 없음. `gabes/schemes/rydberg.py`
마지막 커밋 `d98bdf4`(2026-09-05), SHA-256
`8E658687771478839F4951A712C2AA4EDC4E1337074841A32515B654397BA877`.
2026-10-07 Scheme 3 검토 뒤 생산 코드·전용 테스트 변화 없음.

## 1. 먼저 찾은 기존 제안·TODO·issue

먼저 읽고 검색함.

- `CLAUDE.md`, `README.md` 전체.
- Scheme 3 기존 일일 보고서 17개(2026-06-22…2026-10-07). 직전 보고서 상세 대조.
- `docs/checklist.json`의 Rydberg 완료·활성·parked 항목.
- 현재 scheme/electrometry/experiment/CSV 코드, 전용 테스트.
- cell-heating workflow 문서·예제·테스트, 실험 맥락, User Guide.
- 파일명 `TODO/issue/backlog/roadmap`, 관련 코드 `TODO/FIXME/XXX/HACK`.
- GitHub Issues REST API. 2026-10-12 05:05 KST 기준 공개 issue 0건.

검색 결과:

- Rydberg 경로 inline `TODO/FIXME/XXX/HACK`: 없음.
- 별도 로컬 Rydberg TODO/issue-note 파일: 없음.
- **기존 개선안 충분. 새 거대 물리 로드맵 추가 불필요.**
- 활성 핵심: `docs/checklist.json:3147-3234`, `3568-3581`, `3751-3829`.

기존 항목:

| ID | 우선도·상태 | 핵심 |
|---|---|---|
| `rydberg-extra-view-dependency-correctness` | P1 `ready` | atomic/finite-IF/Beer–Lambert/detector 캐시 경계 |
| `rydberg-effective-opacity-and-sweep-trust` | P1 `ready` | 숨은 `ls=0.001` 공개·provenance·불확도 |
| `rydberg-readout-validity-guardrails` | P1 `needs_decision` | 저대비·저투과·경계·미분해 hero 차단 |
| `rydberg-experimental-validation-dataset` | P2 `blocked_external` | EIT/RF/PSD/온도 held-out 측정 |
| `rydberg-doppler-on-pole-residue` | P3 `ready` | Doppler-on exact-discrete-measure 가속 |
| `rydberg-zeeman-polarization-reference-mode` | P2 `parked` | 제한된 검증형 Zeeman/편광 reference |
| `rydberg-readout-transfer-function-validation` | P2 `blocked_external` | 측정 complex `H(IF)` 검증 |
| `rydberg-rf-field-map-import` | P2 `blocked_external` | provenance 있는 외부 RF map 평균 |

실험 맥락의 기존 status 지적도 유효.
`references/rydberg_ju_experiment_context.md:26-46`의 동일 장치 출판값과
`analysis/rydberg_cell_heating/example_config.json:8-27`의 이상적 detector 가정이
하나의 base-parameter 묶음. source note는 RF dipole을 `REFERENCE`, 전체 electrometry를
`ASSUMED`로 구분하지만 각 입력 provenance는 아직 분리 안 됨
(`example_config.json:80-107`).

- `REFERENCE`: 6 µW, 30 mW, 0.15 mm, 50 mm, 37 GHz, 40 kHz, 약 95% QE.
- `ASSUMED`: 0.5 path split, 완전 balance, RIN correlation=1, electronics=0,
  participation/overlap=1.

입력별 metadata 분리 비용 O(1). spectrum·감도 값 불변.

## 2. 연구자 판정

**판정: 실제 원자물리 있음. compact 4준위 OBE와 finite-IF 검출 사슬은 반정량적
실험 설계·convention 검산·교육/코드 참조에 유용. 절대 sensitivity, 온도 optimum,
Zeeman/편광, 공간 RF장, lock-in 전체 사슬의 정량 기준으로는 부적합.**

유용 범위:

- ⁸⁵Rb `5S₁/₂ F=3 → 5P₃/₂ F′=4 → 40D₅/₂ → 39F₇/₂` cascade EIT/AT.
- probe/coupling power·beam diameter·transit broadening·density·RF dressing 탐색.
- LO-dressed steady state 주위 finite-IF weak-SIG 선형응답.
- AT split↔field, angular/peak/RMS convention, balanced PSN/RIN/electronics 예산 검산.
- heater/effective/cold-spot 분리, axial density·Beer–Lambert, SAM uncertainty,
  CSV unit/hash provenance 실험 분석 골격.

부적합 범위:

- 현재 숫자를 publication-grade 절대 감도 또는 47 °C 운전점으로 인용.
- full Zeeman/polarization·optical pumping·stray field·ionization/charge dynamics.
- horn/cell standing-wave RF map, nonlinear SIG, full time-domain LO+SIG/lock-in.
- 같은 1.6 MHz anchor 재현을 독립 예측 성공으로 해석.

## 3. 실제 구현 물리

### 3.1 원자계·정적 spectrum

- topology·Doppler ratio: `gabes/schemes/rydberg.py:60-91`.
- 4-level atom·5S–40D/Rydberg–Rydberg dephasing 분리:
  `gabes/schemes/rydberg.py:417-437`.
- `TRANSIT_FACTOR=0.6`, `PROBE_RABI_REF_MHZ=2.0`: 기준 1.6 MHz 폭에 맞춘 anchor
  (`gabes/schemes/rydberg.py:41-55`). validation 아님.
- probe/coupling Rabi `√P/d` anchored scaling:
  `gabes/schemes/rydberg.py:439-461`. 절대 optical dipole 예측 아님.
- heater/effective/cold-spot 분리, cold-spot pressure 기반 local density:
  `gabes/schemes/rydberg.py:463-505`.
- detuning-affine steady-state scan: `gabes/schemes/rydberg.py:507-527`,
  `558-619`.
- 기본 Doppler-off = static class. opt-in on = 2 m/s, 4σ Maxwell grid와 residual
  `k` ratio (`gabes/schemes/rydberg.py:614-619`).
- residual Zeeman = full manifold 아님. susceptibility Gaussian convolution
  (`gabes/schemes/rydberg.py:530-542`, `621-627`).
- Beer–Lambert 표준. 단 `ls=0.001` 숨은 effective opacity/participation
  (`gabes/schemes/rydberg.py:629-643`). 공개 atomic line strength 아님.

### 3.2 finite-IF electrometry·detector

- LO-dressed 상태 주위 1차 weak-SIG solve:
  `gabes/rydberg_electrometry.py:116-194`.
- minus sideband 한 번 solve. plus sideband Hermitian-conjugate 복원:
  `gabes/rydberg_electrometry.py:173-193`. 물리 보존 최적화.
- Doppler phasor complex coherent average:
  `gabes/rydberg_electrometry.py:211-227`, `gabes/schemes/rydberg.py:773-819`.
- RF dipole·angular factor·peak/RMS·detuned split:
  `gabes/rydberg_electrometry.py:231-316`.
- photodiode·Schottky noise·correlated RIN·electronics·field ASD:
  `gabes/rydberg_electrometry.py:319-514`.
- 문헌 11.2/12.5 nV/cm/√Hz 결과 강제 주입 없음
  (`tests/test_rydberg_eit.py:80-113`). 좋은 신뢰경계.
- 현재 reference arm power = 각 scan point의 signal power×ratio
  (`gabes/schemes/rydberg.py:902-946`). 실제 고정 reference arm보다 이상적
  pointwise-balance envelope.

### 3.3 실험 helper·예제

- 온도 상태: `gabes/rydberg_experiment.py:80-168`.
- axial density: `gabes/rydberg_experiment.py:191-294`.
- segmented Beer–Lambert: `gabes/rydberg_experiment.py:307-343`.
- Gaussian overlap/effective atom 수: `gabes/rydberg_experiment.py:346-449`.
- SAM field·uncertainty·far-field guard: `gabes/rydberg_experiment.py:485-595`.
- temperature/density dephasing fit·identifiability:
  `gabes/rydberg_experiment.py:598-803`.
- CSV unit·duplicate·SHA-256 provenance:
  `gabes/rydberg_experimental_csv.py:35-100`, `441-691`.
- workflow evidence status:
  `analysis/rydberg_cell_heating/README.md:22-39`.
- 예제 EIT/RF/temperature raw 입력 모두 `PENDING`:
  `analysis/rydberg_cell_heating/example_config.json:53-78`. validation dataset 아님.

## 4. 오늘의 수치 재검증

동일 worktree. headless. warm-up 뒤 median. wall time은 host/BLAS 부하 의존.

| 조건 | 결과 | 해석 |
|---|---:|---|
| 기본 EIT compute | 20.105 ms | compact 정적 4-level solve |
| 기본 EIT `T(0)` / 대비 / FWHM | 0.940016615 / 0.135391492 / 1.613058583 MHz | anchor 규모 재현 |
| 기본 EIT samples/FWHM | 71.69 | 기본 폭 grid 충분 |
| 기본 AT compute | 17.328 ms | 정적 spectrum만 |
| 기본 AT `T(0)` / split | 0.806029283 / 3.4965 MHz | 3.7 MHz LO와 같은 규모. line pulling 있음 |
| 기본 총·PSN ASD | 6.082387 / 6.082387 nV/cm/√Hz | 조건부 예측. 문헌값보다 낙관적 |
| 기본 optimum detuning | −2.220 MHz | 801-point discrete optimum |
| `ls=.0005/.001/.002` 총 ASD | 11.781990 / 6.082387 / 3.238776 | 숨은 opacity가 절대 감도 지배 |
| 80 °C `T(0)` / 대비 / FWHM | 1.798923×10⁻⁶ / 1.798923×10⁻⁶ / 0.307015 MHz | 사실상 암흑인데 numeric 폭 유지 |
| LO=0.1 MHz | AT unresolved / 10.495029 nV/cm/√Hz | 미분해인데 sensitivity hero 유지 |
| Doppler-on compute / 총 ASD | 480.168 ms / 1007.281 nV/cm/√Hz | off 대비 약 24×. calibration regime 격차 큼 |
| cell-heating view | 2.020 s, 원자 compute 20회 | optimum 47 °C / 0.971731. `PREDICTED` |

직전 2026-10-07과 물리 수치 동일. wall time만 크게 변함. 성능 결론은 절대 ms보다
같은 process 내 paired ratio 사용.

### 4.1 온도 optimum 신뢰경계

저장소 실험 맥락: 논문 temperature-vs-PSN optimum 약 33 °C
(`references/rydberg_ju_experiment_context.md:20-22`). 첫 측정 scan은 38 °C에서 좋아졌다가
42.5/44 °C에서 악화. 측정 순서·drift 미분리
(`references/rydberg_ju_experiment_context.md:108-128`).

현재 내장 optimum 47 °C. 운전점 권고 금지.

## 5. 실험-물리 참조 한계

1. **anchor와 prediction 혼동.** probe/transit/coupling 상수가 기준 폭에 맞음.
2. **숨은 opacity.** `ls=0.001`이 transmission·shot noise·sensitivity 지배.
3. **관측가능성 gate 부족.** `_eit_features`는 양의 대비와 half-height crossing만 확인
   (`gabes/schemes/rydberg.py:979-1003`). 80 °C 대비 10⁻⁶도 hero 폭 유지.
4. **미분해 AT에서도 감도 hero.** weak-LO 테스트가 이 동작을 고정
   (`tests/test_rydberg_eit.py:126-141`).
5. **receiver 이상화.** reference arm이 scan point마다 재균형. 고정 power 또는 명시한
   한 balance point 필요.
6. **Doppler off/on은 precision tier 아님.** on은 감도 1007 nV/cm/√Hz. on finite-IF는
   full scan 대신 static-slope 한 점만 계산.
7. **누락 물리 큼.** Zeeman/편광, 공간 RF, measured transfer function, nonlinear/time-domain.
8. **held-out 데이터 없음.** raw EIT/RF/temperature 입력 `PENDING`.

## 6. 기존 개선안: 비용·물리 보존

| 기존 항목 | 계산비용 | 물리·행동 보존 판단 |
|---|---|---|
| effective opacity 공개 | 이름·provenance·status O(1). uncertainty sweep만 opt-in 선형 비용. | 기본 0.001 고정 시 배열·readout `≤1e-12` 보존 가능. atomic line strength와 분리. **최우선.** |
| observability guard | 기존 801 배열의 대비·power·samples/width·edge O(N). 두세 local refinement만 작은 비용. | 원자 spectrum 불변. hero 공개 여부만 차단 가능. fixed-reference readout은 정의 변경이므로 기존 pointwise envelope도 diagnostic 보존. |
| extra-view cache correctness | 네 계층 key·LRU 관리 O(1). detector 변경 때 원자 20회 재계산 제거. | atomic/phasor/Beer–Lambert/detector 경계 분리 시 수치 불변. 현재 finite-IF phasor cache만 부분 구현. |
| held-out dataset | 계산보다 측정·metadata·uncertainty·hash 비용 지배. | simulation 대체 불가. 완료 전 absolute/optimum은 `PREDICTED`/`UNSUPPORTED`. |
| Doppler pole-residue | scan row당 작은 eigendecomposition. 현재 velocity-class solve 제거. | 2 m/s·4σ 이산 측도와 coherent phasor, per-row direct guard 유지 시 behavior-preserving. default off라 P3 타당. |
| Zeeman/polarization reference | Hilbert/Liouville 차원·B/velocity grid·dipole provenance. 무시 가능 overhead 아님. | 특정 state/polarization/B dataset 뒤 opt-in만 타당. |
| measured `H(IF)` | 측정 complex factor 적용 O(N_IF). 작음. | 원자 응답 보존. 외부 데이터 없어 blocked. full time-domain보다 먼저. |
| RF field-map import | O(N_space×N_scan), mesh convergence 필요. 큼. | uniform-map reduction·외부 uncertainty gate 필요. 내부 3-D EM solver는 불필요. |
| 입력별 evidence status | metadata O(1). | 계산 완전 불변. published 값과 ideal assumption 분리. |

기존 제안이 물리 공백과 trust boundary를 이미 포괄. 오늘 새 physics 항목 추가 안 함.

## 7. 순수 코드 최적화

### 7.1 BLAS controller 지연 캐시

현재 `core.blas_single_thread()`는 매 solve마다 `threadpool_limits(...)` 객체 생성
(`gabes/core.py:23-34`). process 내 `ThreadpoolController` 한 번 생성 뒤 limiter 재사용한
prototype:

- 기본 EIT median: **22.5093 → 4.3460 ms**.
- **5.18×**.
- `chi_bar np.array_equal=True`.

원자물리·solver·부동소수점 결과 불변. NumPy/BLAS 로드 뒤 lazy init, 중첩 context,
병렬 호출, process-wide thread limit 복원 회귀 필요.

### 7.2 finite-IF Liouvillian affine 조립

현재 Doppler-off finite-IF: 801 detuning마다 Python Hamiltonian/Liouvillian 재조립
(`gabes/schemes/rydberg.py:804-815`). detuning 의존 affine.

- list 조립 median: **220.5785 ms**.
- `L(scan)=L(0)+scan·[L(1)-L(0)]`: **20.3237 ms**.
- **10.85×**, `np.array_equal=True`, max abs difference 0.

조립 microbenchmark. 전체 solve 속도 향상률과 혼동 금지. 물리·scan·solver 불변.

### 7.3 detector/noise scan 벡터화

현재 801점마다 dataclass 생성·검증·noise/sensitivity 계산
(`gabes/schemes/rydberg.py:910-946`). 식은 elementwise.

- scalar loop median: **100.5659 ms**.
- NumPy vector 식: **0.0368 ms**.
- 약 **2733×**. max relative difference `3.10×10⁻¹⁶`.

역시 loop microbenchmark. 전체 finite-IF 속도 향상률 아님. scalar fixture와 machine-roundoff
gate 필요.

### 7.4 extra-view 계층 cache

분리:

1. atomic steady-state spectrum.
2. finite-IF atomic phasor.
3. cell/opacity Beer–Lambert.
4. detector/RF-dipole/SAM conversion.

현재 2번 일부만 LRU. temperature view는 10 온도×EIT/AT = 원자 compute 20회
(`gabes/schemes/rydberg.py:1384-1441`). detector electronics만 바뀌어도 전체 view 재계산.
계층 분리 시 수치 정의 불변.

### 7.5 작은 정리

- temperature sweep이 표시 문자열 `Max spectral slope`, `IF discriminator`를 split해 float
  재파싱(`gabes/schemes/rydberg.py:1412-1423`). `_readout` numeric field 직접 전달.
- EIT temperature loop는 `_readout` 뒤 `_transmission` 재호출. 이미 계산한 배열 재사용.

## 8. 문서·테스트·예제 적합성

좋음:

- `README.md:117-137`: finite-IF, 온도 분리, SAM, deferred 범위 제한.
- `analysis/rydberg_cell_heating/README.md:22-39`, `78-117`: evidence status,
  provenance, lumped axial 한계, conditional electrometry 명시.
- sideband residual·Hermiticity·1회 solve:
  `tests/test_rydberg_electrometry.py:39-94`.
- low-IF→static derivative: `tests/test_rydberg_eit.py:309-337`.
- detector-only phasor reuse: `tests/test_rydberg_eit.py:340-427`.
- temperature/cold spot·axial/SAM·dephasing identifiability:
  `tests/test_rydberg_eit.py:429-517`, `tests/test_rydberg_experiment.py:19-200`.
- CSV unit·duplicate·hash·음수 PSD 거부:
  `tests/test_rydberg_experimental_csv.py:14-109`.

부족:

- User Guide는 정적 EIT/AT만 설명
  (`docs/Userguide/GABES_User_Guide_v2.html:535-538`, `617-625`). finite-IF 감도,
  conditional calibration, temperature/provenance, 숨은 opacity 한계 누락.
- `tests/test_rydberg_eit.py:1-5`의 “public UI shows the static spectrum only”는 현재 UI와 불일치.
- 예제 base params의 published 값과 ideal detector 가정 provenance가 입력별로 분리 안 됨.
- 80 °C low-transmission, fixed-reference, scan-edge optimum, local-refinement acceptance fixture 없음.
- opacity provenance와 Doppler mode별 calibration scope 직접 고정 테스트 없음.

## 9. 검증

- 전용 Scheme/helper/workflow/adapter:
  `69 passed in 77.07 s`.
- 필수 `python -m pytest -q`:
  `2202 passed, 4 skipped, 2 failed in 5891.90 s`.
- 두 실패 모두 Rydberg 아님. 직전 전체 회귀와 동일:
  - `tests/quantum/test_thermal_campaign_extend.py`의 fresh-capsule extraction이 Windows
    임시 `.git/grand-challenge-runs/...` 중첩 경로에서
    `analysis/grand_challenge/reference/ballistic_linear_channel.py` destination을 만들지
    못해 `FileNotFoundError`.
  - `tests/test_fwm_excess_noise.py`의 Streamlit ScrubField fallback에서
    gain-proportional mode 전환 뒤 metric이 이전 constant 값 2.0 유지. 기대 0.04.

전체 회귀 통과 아님. 대상 Rydberg 전용 69개는 전부 통과.

생산 코드 수정 없음. 이번 보고서만 추가.

## 10. 권장 순서

1. **P1 opacity provenance 공개.** runtime 사실상 0. 절대값 신뢰도 효과 최대.
2. **P1 observability/fixed-reference gate.** 80 °C·weak-LO 반례 fixture부터.
3. **P1 순수 성능.** cached BLAS controller → affine Liouvillian 조립 → detector
   vectorization → four-layer cache.
4. **P1 provenance/status·문서.** published input과 ideal assumption 분리. User Guide,
   stale test docstring, checklist partial 상태 정정.
5. **P2 held-out dataset.** EIT/RF/PSD/temperature, fixed reference, transfer function.
6. **P3 Doppler pole.** 실제 Doppler-on sweep 병목 뒤 현재 이산 측도 보존 이식.
7. Zeeman/편광·RF map·time-domain은 측정·runtime budget 전까지 opt-in/blocked 유지.

핵심: 새 거대 엔진보다 **opacity provenance, 관측가능성, receiver 조건,
입력별 evidence status, held-out 검증** 먼저.
