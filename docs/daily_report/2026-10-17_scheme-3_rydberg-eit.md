# GABES 일일 Scheme 물리 검토 — 2026-10-17

## 0. 선택·현재 다섯 scheme

Asia/Seoul 현지 17일.

`n = (17 mod 5) + 1 = 3`.

실제 드롭다운 등록 순서(`gabes/schemes/__init__.py:19-39`):

1. OD / SAS — `SASScheme`
2. Lambda coherence (EIT / AT / CPT) — `LambdaScheme`
3. **Rydberg-EIT electrometry — `RydbergEITScheme`**
4. Hanle / EIA / NMOR — `MagnetoScheme`
5. FWM — `FWMScheme`

Lambda·Magneto 단일-regime 객체는 alias. 드롭다운 항목 아님.

검토 worktree: `main@39be4ae`. 기존 staged/unstaged/untracked 370건 있음. 보존함.
Rydberg 생산 코드·전용 테스트 diff 없음. 핵심 scheme SHA-256
`8E658687771478839F4951A712C2AA4EDC4E1337074841A32515B654397BA877`.
직전 2026-10-12 검토와 같음. 마지막 관련 생산 커밋 `d98bdf4`(2026-09-05).

## 1. 선행 제안부터 검색

먼저 확인함.

- `CLAUDE.md`, `README.md`.
- Head Manager memory·현행 blueprint 카드·`NEXT_SESSION.md`. 현행 Grand Challenge 카드는
  FWM gain 보정. 이번 Scheme 3 감사와 직접 의존성 없음.
- Scheme 3 과거 일일 보고서 18개(2026-06-22…2026-10-12).
- `docs/checklist.json`의 Rydberg 완료·활성·parked 항목.
- Rydberg scheme/electrometry/experiment/CSV 코드, 전용 테스트, heating workflow·예제,
  User Guide, 실험 맥락.
- 파일명 TODO/issue/backlog/roadmap, 관련 코드 `TODO/FIXME/XXX/HACK`.
- GitHub Issues API. 2026-10-17 KST 기준 open/closed issue 0건.

Rydberg inline TODO·별도 로컬 issue-note 없음. 기존 개선안 충분. 새 대형 물리 제안 불필요.

| 기존 ID | 현재 상태 | 위치 |
|---|---|---|
| `rydberg-extra-view-dependency-correctness` | P1 `ready`; 아래 fresh 검증상 문제 문구 일부 stale | `docs/checklist.json:3469-3488` |
| `rydberg-effective-opacity-and-sweep-trust` | P1 `ready` | `docs/checklist.json:3489-3511` |
| `rydberg-readout-validity-guardrails` | P1 `needs_decision` | `docs/checklist.json:3512-3536` |
| `rydberg-experimental-validation-dataset` | P2 `blocked_external` | `docs/checklist.json:3537-3561` |
| `rydberg-doppler-on-pole-residue` | P3 `ready` | `docs/checklist.json:3890-3907` |
| `rydberg-zeeman-polarization-reference-mode` | P2 `parked` | `docs/checklist.json:4073-4099` |
| `rydberg-readout-transfer-function-validation` | P2 `blocked_external` | `docs/checklist.json:4100-4126` |
| `rydberg-rf-field-map-import` | P2 `blocked_external` | `docs/checklist.json:4127-4153` |

## 2. 현재 구현 물리

### 2.1 원자·광학

- ⁸⁵Rb `5S₁/₂ F=3 → 5P₃/₂ F′=4 → 40D₅/₂ → 39F₇/₂` 4준위 cascade
  (`gabes/schemes/rydberg.py:4-15`, `73-91`).
- spontaneous decay와 5S–40D / 40D–39F dephasing 분리
  (`gabes/schemes/rydberg.py:417-437`).
- probe·coupling Rabi는 `√P/d` anchored scaling. 절대 optical dipole 예측 아님
  (`gabes/schemes/rydberg.py:41-55`, `439-461`).
- heater setpoint, effective vapor temperature, cold spot 분리. cold-spot 압력과 local
  ideal-gas density 사용(`gabes/schemes/rydberg.py:463-505`;
  `gabes/rydberg_experiment.py:52-167`).
- detuning-affine steady-state OBE. 기본 Doppler-off 한 class. opt-in on은 2 m/s·4σ
  Maxwell grid와 residual `(k_probe-k_coupling)v`
  (`gabes/schemes/rydberg.py:507-527`, `558-627`).
- residual Zeeman은 full manifold 아님. susceptibility Gaussian convolution
  (`gabes/schemes/rydberg.py:530-542`, `621-627`).
- Beer–Lambert. 단 line-strength 인수 `ls=0.001` 하드코딩
  (`gabes/schemes/rydberg.py:629-643`). 실제 atomic line strength보다 effective
  opacity/participation 보정으로 해석해야 함.

### 2.2 finite-IF electrometry·검출

- LO-dressed steady state 주위 1차 weak-SIG finite-IF Liouvillian 응답
  (`gabes/rydberg_electrometry.py:116-194`).
- minus sideband 한 번 solve, plus는 Hermitian conjugate 복원. 물리 보존 최적화
  (`gabes/rydberg_electrometry.py:173-193`;
  `tests/test_rydberg_electrometry.py:39-93`).
- Doppler class complex phasor coherent average
  (`gabes/rydberg_electrometry.py:211-227`; `gabes/schemes/rydberg.py:773-819`).
- RF dipole, angular factor, peak/RMS, detuned AT split 단위 명시
  (`gabes/rydberg_electrometry.py:231-316`).
- one-sided Schottky PSN, correlated RIN, electronics, field ASD
  (`gabes/rydberg_electrometry.py:319-483`). 문헌 11.2/12.5 nV/cm/√Hz 강제 주입 안 함
  (`tests/test_rydberg_eit.py:80-97`).
- 현재 reference arm은 각 probe-detuning 점의 transmitted signal에 비례해 다시 맞춤
  (`gabes/schemes/rydberg.py:902-946`). 실제 fixed reference arm보다 이상적인 envelope.

### 2.3 실험 helper·예제

- axial cold-spot-limited density와 segmented Beer–Lambert
  (`gabes/rydberg_experiment.py:191-343`).
- Gaussian overlap atom 수·SAM field/uncertainty·far-field guard
  (`gabes/rydberg_experiment.py:346-595`).
- temperature/density dephasing fit은 rank·condition으로 식별 불가를 표시
  (`gabes/rydberg_experiment.py:598-800`).
- EIT/RF/PSD CSV는 단위 변환, duplicate median, SHA-256, parsing assumptions 보존
  (`gabes/rydberg_experimental_csv.py:34-103`, `441-713`).
- workflow는 `MEASURED/PPT/REFERENCE/FITTED/PREDICTED/ASSUMED/PENDING` 분리
  (`analysis/rydberg_cell_heating/README.md:22-39`).
- 예제 raw EIT/RF/temperature 입력은 모두 `PENDING`
  (`analysis/rydberg_cell_heating/example_config.json:53-78`). held-out validation 아님.

## 3. ⟦AI⟧ 실험 원자물리 연구자 판정

**판정: 실제 물리 있음. compact 4준위 cascade EIT/AT, finite-IF 원자 mixer,
RF-field convention, balanced detector, heated-cell 분석 골격은 반정량적 실험 설계·단위
검산·코드 reference에 유용. 절대 sensitivity·온도 optimum·Zeeman/편광·공간 RF장·전체
lock-in의 정량 reference에는 아직 부적합.**

쓸 만한 범위:

- EIT/AT line-shape 추세, probe/coupling/waist/transit/dephasing scan.
- LO/SIG finite-IF 선형응답 convention 검산.
- AT↔field, peak/RMS, PSN/RIN/electronics 예산 sanity check.
- temperature/cold spot/axial density/SAM/CSV provenance 분석 골격.

쓰면 안 되는 범위:

- 현재 6.08 nV/cm/√Hz 또는 47 °C를 publication-grade 예측/운전점으로 인용.
- full Zeeman/polarization·optical pumping·stray field·ionization/charge 포함 주장.
- horn/cell standing-wave, nonlinear SIG, measured `H(IF)`, lock-in bandwidth 재현 주장.
- 1.6 MHz anchor 재현을 독립 물리 검증으로 해석.

## 4. fresh 수치 점검

현재 worktree, headless. 핵심 값은 2026-10-12와 동일.

| 조건 | 현재 결과 | 의미 |
|---|---:|---|
| 기본 EIT `T(0)` / 대비 / FWHM | 0.940016615 / 0.135391492 / 1.613058583 MHz | anchor 규모 재현 |
| 기본 AT split | 3.4965 MHz | 3.7 MHz LO와 같은 규모. optical line pulling 있음 |
| 기본 총 / PSN ASD | 6.082387 / 6.082387 nV/cm/√Hz | 이상적 detector 조건의 conditional prediction |
| 기본 superhet optimum detuning | −2.220 MHz | 801점 discrete optimum |
| `ls=.0005/.001/.002` 총 ASD | 11.781990 / 6.082387 / 3.238776 | 숨은 opacity가 절대 감도 지배 |
| 80 °C `T(0)` / 대비 / FWHM | 1.798923×10⁻⁶ / 1.798923×10⁻⁶ / 0.307015 MHz | 사실상 암흑인데 numeric width 유지 |
| LO=0.1 MHz | AT unresolved / 10.495029 nV/cm/√Hz | 미분해인데 sensitivity hero 유지 |
| Doppler-on 총 ASD | 1007.281 nV/cm/√Hz | off와 두 자릿수 이상 모델 격차 |
| 10점 heating sweep optimum | 47 °C / 0.971731 nV/cm/√Hz | 외부 paper-context optimum 약 33 °C와 불일치 |

실험 맥락은 paper model optimum 약 33 °C, 최초 heating scan은 38 °C 최선이나 측정 순서와
drift 미분리라고 기록(`references/rydberg_ju_experiment_context.md:20-22`, `108-128`).
현재 47 °C는 운전점 권고 아님.

## 5. ⟦AI⟧ 기존 개선안 비용·물리 보존 평가

| 기존 제안 | 계산비용 | 현재 판단 |
|---|---|---|
| effective opacity 공개 | 이름·provenance·status O(1); 불확도 sweep만 opt-in 선형 비용 | 기본 0.001 유지 시 수치 `≤1e-12` 보존 가능. 절대값 신뢰도 최우선. |
| observability guard | 기존 801 배열의 contrast·power·samples/width·edge O(N); 국소 refinement 소수 | 원자 solve 불변. hero 공개만 차단 가능. fixed-reference는 정의 변경이므로 기존 pointwise envelope도 diagnostic 보존. |
| extra-view cache correctness | key 검사는 O(1); 계층 분리는 medium 구현·실행비용 감소 | 아래 fresh 검증상 key 정확성은 현 코드에서 이미 충족. checklist의 3.5× 불일치 문구 정정 필요. 남은 것은 cache 계층 성능. |
| held-out dataset | 계산보다 측정·metadata·uncertainty 비용 지배 | simulation으로 대체 불가. 완료 전 absolute/optimum은 `PREDICTED`/`UNSUPPORTED`. |
| measured `H(IF)` | complex factor 적용 O(N_IF), 작음 | 원자 응답 보존. 측정 transfer function 없어 blocked. |
| Doppler pole-residue | scan row당 작은 eigendecomposition; velocity-class 반복 제거 | 현재 2 m/s·4σ 이산 측도와 direct guard 유지하면 behavior-preserving. default off라 P3 적절. |
| Zeeman/polarization reference | Hilbert/Liouville 차원·B/velocity grid 증가. 큼 | 특정 state/polarization/B dataset 전에는 parked 타당. |
| RF field map | O(N_space×N_scan), mesh refinement 필요. 큼 | uniform-map 1% reduction·외부 uncertainty 필요. 내부 3-D EM solver 불필요. |
| 입력별 evidence metadata·Guide 정정 | O(1) | 계산 불변. published 입력과 ideal detector assumption 분리 가능. |

### 5.1 stale extra-view 문제 재검증

실제 UI는 `recompute_keys ∪ ExtraView.param_keys`를 cache key와 compute 입력으로 사용
(`gabes_ui/plotcard.py:62-67`, `723-728`). Rydberg는 probe sweep에 `cell_mm`, heating
sweep에 IF·RF dipole·detector·participation 입력을 선언
(`gabes/schemes/rydberg.py:1528-1556`). 전용 test도 현재 목록을 고정
(`tests/test_rydberg_eit.py:45-56`).

모든 navigate-only readout 값을 크게 바꾼 한 fresh fixture에서:

- probe-power view: 전체 params 직접 계산 vs 실제 filtered cache params, diff 없음.
- heating view: 동일 비교, 모든 반환 numeric/string field diff 없음.

따라서 `docs/checklist.json:3473`의 “cached/direct 3.5× 차이”는 현 코드에서 재현 안 됨.
단, detector만 바뀌어도 extra-view 전체 함수가 다시 실행되어 heating loop의 EIT/AT
`self.compute` 20회가 반복됨(`gabes/schemes/rydberg.py:1403-1442`). correctness는 닫고,
“atomic spectrum → finite-IF phasor → Beer–Lambert → detector/SAM” 성능 계층화만 별도
항목으로 남기는 편이 정확함. 현재 test는 key 목록만 검사. filtered/direct 결과 parity도 추가 필요.

## 6. ⟦AI⟧ 동작 불변 코드 최적화

1. **BLAS controller 지연 캐시.** 현재 매 solve마다 `threadpool_limits` 객체 생성
   (`gabes/core.py:23-34`). 2026-10-12 paired prototype는 EIT median
   22.5093→4.3460 ms, 5.18×, `chi_bar array_equal=True`. NumPy/BLAS 로드 뒤 lazy init,
   nested context·병렬 호출·복원 test 필요
   (`docs/daily_report/2026-10-12_scheme-3_rydberg-eit.md:213-224`).
2. **finite-IF Liouvillian affine 조립.** 현재 801 detuning마다 Python Hamiltonian/
   Liouvillian 재조립(`gabes/schemes/rydberg.py:807-815`). 같은 scan을
   `L(s)=L(0)+s[L(1)-L(0)]`로 batch 조립 가능. 직전 microbenchmark 220.5785→20.3237 ms,
   10.85×, array-equal. solve 물리·grid 불변
   (`docs/daily_report/2026-10-12_scheme-3_rydberg-eit.md:226-235`).
3. **detector/noise scan 벡터화.** 801점마다 dataclass 생성·검증
   (`gabes/schemes/rydberg.py:908-946`). elementwise NumPy 식으로 직전
   100.5659→0.0368 ms, 최대 상대차 `3.10×10⁻¹⁶`. scalar fixture gate 필요
   (`docs/daily_report/2026-10-12_scheme-3_rydberg-eit.md:237-247`).
4. **extra-view 4계층 cache.** detector 변경 때 static spectrum 20회 재계산 제거.
   현재 finite-IF atomic phasor LRU만 존재(`gabes/schemes/rydberg.py:821-865`). 수치 정의 불변.
5. **작은 정리.** heating view가 formatted metric 문자열을 float로 다시 parsing하고 EIT
   transmission을 중복 계산(`gabes/schemes/rydberg.py:1407-1423`). `_readout` numeric field와
   이미 계산한 배열 직접 재사용. 물리 불변.

`BLAS cache → affine assembly → detector vectorization → 계층 cache` 순서 권장.

## 7. ⟦AI⟧ 문서·테스트·예제 적합성

좋음:

- README는 finite-IF, conditional absolute sensitivity, temperature/SAM, deferred 범위 명시
  (`README.md:139-161`).
- workflow는 evidence status와 lumped axial 한계 명시
  (`analysis/rydberg_cell_heating/README.md:22-39`, `78-117`).
- sideband residual·Hermiticity·한 번 solve, low-IF static derivative, detector-only phasor
  reuse 검증(`tests/test_rydberg_electrometry.py:39-119`;
  `tests/test_rydberg_eit.py:309-427`).
- temperature/cold spot·axial/SAM·dephasing identifiability·CSV 단위/hash/음수 거부 검증
  (`tests/test_rydberg_eit.py:456-531`; `tests/test_rydberg_experiment.py:19-200`;
  `tests/test_rydberg_experimental_csv.py:14-109`).

부족:

- User Guide는 정적 EIT/AT와 split만 설명. finite-IF·conditional sensitivity·temperature/
  opacity·provenance 누락(`docs/Userguide/GABES_User_Guide_v2.html:535-538`, `617-625`).
- test module docstring의 “public UI shows the static spectrum only”는 현재 기능과 불일치
  (`tests/test_rydberg_eit.py:1-5`).
- example config는 같은 `base_params` 안의 published apparatus 값과 ideal detector 가정을
  입력별 provenance로 분리하지 않음(`analysis/rydberg_cell_heating/example_config.json:8-27`,
  `80-107`).
- 80 °C 저투과, fixed-reference, scan-edge optimum, local refinement acceptance test 없음.
- extra-view test는 목록만 pin. 실제 cache-key/full-param parity 미검사.

## 8. 검증

- fresh physics fixture: 기본 EIT/AT, opacity 3점, 80 °C, weak-LO, Doppler-on,
  10점 heating sweep 완료. 위 표 수치.
- extra-view altered navigate-input filtered/direct parity: 두 view 모두 diff 없음.
- Rydberg scheme/electrometry/experiment/CSV/workflow/adapter: **69 passed in 24.06 s**.
- 필수 `python -m pytest -q`: **2434 passed, 4 skipped, 1 failed, 5 warnings in
  2664.92 s**.
- 유일 실패는 Rydberg 아님. 직전 전체 회귀와 같은 Windows nested fresh-capsule
  `FileNotFoundError`: `tests/quantum/test_thermal_campaign_extend.py:367` →
  `analysis/grand_challenge/reference/ballistic_linear_channel.py` destination directory 없음.
  warning 5건은 `gabes_ui/export.py:197`의 naive `datetime.utcnow()` deprecation.

전체 회귀 통과 아님. 대상 Rydberg 69개는 전부 통과.

생산 코드·Grand Challenge state 수정 없음. 오늘 보고서만 추가.

## 9. ⟦AI⟧ 권장 순서

1. P1 effective-opacity 이름·provenance·`PREDICTED` status 공개. 실행 overhead 사실상 0.
2. P1 observability/fixed-reference gate. 80 °C·weak-LO·edge fixture 먼저.
3. stale cache item 정리: dependency correctness 완료/검증, 남은 계층 cache는 성능 항목으로 분리.
4. 순수 최적화: BLAS controller, affine assembly, detector vectorization, 4계층 cache.
5. 입력별 evidence metadata와 User Guide/test docstring 동기화.
6. held-out EIT/RF/PSD/temperature dataset 뒤 transfer function·Zeeman·RF map 판단.
7. Doppler-on이 실제 sweep 병목일 때 pole-residue 이식.

핵심: 현재 장점은 실제 finite-IF 원자응답과 명시적 검출 convention. 현재 약점은 계산식보다
**effective opacity·관측가능성·receiver 조건·held-out 검증**. 새 거대 엔진보다 이 경계 먼저.
