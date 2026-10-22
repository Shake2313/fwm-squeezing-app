# GABES 일일 Scheme 물리 검토 — 2026-09-07

## 0. 오늘의 선택과 범위

- 로컬 날짜의 일(day) `7`에 대해 `n = (7 mod 5) + 1 = 3`이므로 **Scheme 3 — Rydberg-EIT electrometry**를 검토했다.
- 현재 registry 순서는 `gabes/schemes/__init__.py:19-25`의 ① OD/SAS, ② Lambda EIT/AT/CPT, ③ Rydberg-EIT electrometry, ④ Hanle/EIA/NMOR, ⑤ FWM이다. alias는 `gabes/schemes/__init__.py:27-39`에 따로 있어 이 순서를 바꾸지 않는다.
- 검토 기준은 현재 worktree의 Scheme 3 구현이다. 대상 Scheme 3 파일은 `HEAD`의 `d98bdf4`와 같았으며, worktree의 다른 staged/unstaged 변경은 수정하지 않았다.
- 먼저 `docs/daily_report`의 기존 Scheme 3 보고서 11개(2026-06-22…2026-09-02), `docs/checklist.json`, 코드의 TODO/FIXME, 분석 예제, GitHub issue를 확인했다. 2026-09-07 기준 공개 issue는 open/closed 모두 0개였고, 실제 미완료 항목은 checklist에 집중돼 있다.

## 1. 결론

**판정: 실험 설계·추세 탐색·데이터 정리용으로는 유용하지만, 절대 감도와 온도 최적점을 reference로 인용하기에는 아직 반정량적(semi-quantitative)이다.**

현재 구현은 다음 용도에는 충분히 강하다.

- ⁸⁵Rb `5S–5P–40D` cascade와 약 37 GHz `40D–39F` RF leg의 EIT/AT 형상을 빠르게 비교한다.
- LO-dressed 정상상태 주위의 유한 IF 약신호 응답과 balanced detector PSN/RIN/electronics budget을 한 단위계로 연결한다.
- heater setpoint, effective vapor temperature, cold-spot temperature를 구분하고 실험 CSV, axial profile, SAM 보정을 붙인다.
- 모든 입력의 provenance가 갖춰졌을 때 “어떤 항목이 감도를 제한하는가”를 조건부로 비교한다.

그러나 다음 이유로 **절대 감도 표준, 보정 곡선, 온도 optimum의 실험 reference**로 쓰면 안 된다.

1. `ls=0.001` 유효 opacity 계수가 공개 파라미터가 아닌 채 감도에 약 2배 수준의 영향을 준다 (`gabes/schemes/rydberg.py:48-51`, `gabes/schemes/rydberg.py:641-644`).
2. 저투과·저대비·scan-edge optimum을 차단하는 observability gate가 없다. 80 °C fixture는 contrast가 `1.80×10⁻⁶`인데도 유한 FWHM을 내고, IF optimum은 scan 끝에서 불과 0.04 MHz 안쪽이다 (`gabes/schemes/rydberg.py:978-1003`, `gabes/schemes/rydberg.py:1083-1094`).
3. balanced detector reference arm이 scan의 각 점에서 signal transmission을 따라 다시 설정된다. 이는 고정 reference 또는 한 detuning에서 balance한 실험 광학계와 다르다 (`gabes/schemes/rydberg.py:916-927`).
4. 측정 EIT/RF/PSD/temperature의 predeclared held-out dataset이 아직 없다 (`docs/checklist.json:583-603`; `analysis/rydberg_cell_heating/example_config.json:53-78`).
5. full Zeeman/polarization, 공간 RF field, 측정된 mixer/lock-in transfer는 의도적으로 빠져 있다 (`gabes/schemes/rydberg.py:1-15`; `README.md:135-137`).

따라서 UI와 보고서의 현재 “conditional prediction” 표현은 유지해야 하며, 숫자는 같은 모델 안의 상대 비교에 우선 사용해야 한다 (`README.md:117-129`).

## 2. 현재 물리 구현 검토

### 2.1 원자계와 optical response

- `RydbergEITScheme`은 4-level compact cascade를 만든다 (`gabes/schemes/rydberg.py:72-104`, `_cascade_skeleton`; `gabes/schemes/rydberg.py:417-437`, `_atom`).
- probe/coupling power와 beam diameter를 Rabi frequency로 변환하되, coupling은 문헌 조건 30 mW → 3 MHz에 anchored되어 있다 (`gabes/schemes/rydberg.py:439-461`). 이는 실험 장치 전체의 절대 beam-field calibration이 아니라 명시적 reference scaling이다.
- density는 cold spot, thermal motion과 temperature-dependent dephasing은 effective vapor temperature를 사용한다 (`gabes/schemes/rydberg.py:463-505`). 이 분리는 heated-cell 실험에서 물리적으로 적절하다.
- spectrum은 801점 고정 scan이며 affine detuning kernel 또는 pointwise fallback으로 susceptibility를 계산한다 (`gabes/schemes/rydberg.py:507-527`, `gabes/schemes/rydberg.py:544-556`).
- uncompensated magnetic broadening은 full Zeeman population/selection-rule solve가 아니라 Gaussian convolution proxy다 (`gabes/schemes/rydberg.py:529-542`, `gabes/schemes/rydberg.py:621-629`). 현재 설명처럼 비교용 proxy로만 해석해야 한다.
- Doppler-on finite-IF 경로는 최대 static slope 한 점의 velocity phasor를 coherent average한다 (`gabes/schemes/rydberg.py:773-807`). 이것은 빠른 first-order estimator이지 모든 detuning의 완전한 velocity-resolved RF readout은 아니다.

### 2.2 finite-IF electrometry와 detector

- 약신호 solver는 LO-dressed Liouvillian에서 정상상태를 구하고 `L₀+iω_IF I` 선형계를 푼다 (`gabes/rydberg_electrometry.py:116-185`).
- `d98bdf4`에서는 Hermitian RF perturbation과 density-matrix Hermiticity를 이용해 minus sideband 하나만 풀고 plus sideband를 conjugate transpose로 복원한다 (`gabes/rydberg_electrometry.py:173-185`). 물리식을 바꾸지 않는 타당한 최적화이며 관련 회귀 테스트가 통과한다.
- Doppler velocity class는 complex phasor 상태에서 평균하므로 magnitude를 먼저 평균하는 오류를 피한다 (`gabes/rydberg_electrometry.py:211-227`).
- RF transition dipole, angular factor, peak/RMS convention, detector QE/path/reference/RIN/electronics가 분리돼 있어 calibration ledger의 구조는 좋다 (`gabes/schemes/rydberg.py:683-750`, `gabes/schemes/rydberg.py:867-976`).
- 다만 detector reference optical power가 매 detuning에서 `signal_power × reference_ratio`로 계산된다 (`gabes/schemes/rydberg.py:916-927`). scan 동안 고정하거나 명시한 balance detuning에서 고정해야 실제 balanced receiver와 대응한다.

### 2.3 실험 분석 경로

- `analysis/rydberg_cell_heating` workflow는 clean CSV contract, temperature semantics, axial Beer–Lambert, effective atom number, SAM uncertainty를 분리한다 (`analysis/rydberg_cell_heating/README.md:22-39`, `analysis/rydberg_cell_heating/README.md:41-94`, `analysis/rydberg_cell_heating/README.md:96-117`).
- 문헌과 같은 apparatus의 공개 조건은 probe 6 µW, coupling 30 mW, diameter 0.15 mm, cell 50 mm, coupling Rabi 3 MHz, LO 3.7 MHz, IF 40 kHz, detector QE 95%다. 현재 defaults 테스트도 이를 고정한다 (`tests/test_rydberg_eit.py:31-42`).
- 같은 apparatus의 보고값은 PSN 11.2 nV/cm/√Hz, 측정 12.5(8) nV/cm/√Hz이고, B compensation은 linewidth 약 1.9 → 1.6 MHz와 peak 약 18% 증가를 보였다. 현재 기본 FWHM 1.613 MHz는 해당 scale과 맞지만, 감도 6.082 nV/cm/√Hz는 PSN reference보다 약 1.84배 낙관적이다. 이 차이는 calibration 완료가 아니라 미식별 opacity/detector/overlap 자유도의 증거로 봐야 한다.

## 3. 현재 코드 수치 점검

Windows/Python 3.12 로컬 실행에서 현재 worktree를 직접 계산했다. timing은 단일 프로세스 wall clock이며 장비 독립 성능 보장은 아니다.

| 점검 | 현재 결과 | 해석 |
|---|---:|---|
| 기본 EIT 중심 transmission | 0.940016615 | 정상 범위 |
| 기본 EIT contrast | 0.135391492 | 기본점은 관측 가능 |
| 기본 EIT FWHM | 1.613058583 MHz | 보상 후 실험 scale과 근접 |
| 기본 AT 중심 transmission | 0.806029283 | detector power가 충분한 기본점 |
| 기본 AT split | 3.4965 MHz | LO Rabi 3.7 MHz와 같은 scale |
| finite-IF optimum | -2.2200 MHz | 기본 scan 내부 |
| 기본 PSN = total sensitivity | 6.082387215 nV/cm/√Hz | 문헌 PSN 11.2보다 1.84배 낙관적 |
| `ls=0.0005` total sensitivity | 11.781990083 nV/cm/√Hz | hidden opacity의 큰 영향 |
| `ls=0.001` total sensitivity | 6.082387215 nV/cm/√Hz | 현재 기본값 |
| `ls=0.002` total sensitivity | 3.238775723 nV/cm/√Hz | 4× opacity 범위에서 3.64× 변화 |
| 80 °C EIT contrast | 1.798922847×10⁻⁶ | 사실상 관측 불가지만 차단되지 않음 |
| 80 °C EIT FWHM | 0.307015017 MHz | 정밀해 보이는 숫자가 계속 노출됨 |
| 80 °C IF optimum | -8.96 MHz | ±9 MHz scan edge에서 0.04 MHz |

성능 측정:

| 경로 | 시간 | 메모 |
|---|---:|---|
| 첫 EIT compute | 2250.7 ms | JIT cold start 포함 |
| warm EIT compute 중앙값, 7회 | 3.883 ms | 정적 interactive 사용에 충분 |
| warm AT compute 중앙값, 7회 | 4.000 ms | 정적 spectrum만 |
| cached finite-IF readout 중앙값, 7회 | 22.465 ms | atomic phasor cache hit 후 detector scan 포함 |
| heating extra view, 전체 입력 | 816.1 ms | 10개 온도에서 EIT+AT |
| detector noise만 변경한 heating view | 317.2 ms | atomic phasor는 재사용하지만 정적 spectrum은 다시 계산 |

추가 뷰 dependency correctness는 9월 2일 보고 당시보다 개선됐다. `ExtraView.param_keys`가 추가됐고 (`gabes/schemes/base.py:80-86`), 앱은 recompute key와 view별 live key의 합집합으로 extra cache를 만든다 (`streamlit_app.py:1982-1988`). Rydberg view는 cell/IF/RF dipole/detector/SAM 의존성을 선언한다 (`gabes/schemes/rydberg.py:1528-1555`), 테스트도 key와 atomic-cache 경계를 고정한다 (`tests/test_rydberg_eit.py:45-56`, `tests/test_rydberg_eit.py:355-427`).

다만 `docs/checklist.json:515-532`의 acceptance 전체가 완료된 것은 아니다. detector-only 변경에서 finite-IF atomic solve는 피하지만 20개의 정적 EIT/AT compute와 string metric 추출이 반복되므로, “atomic steady-state / finite-IF phasor / Beer–Lambert / detector-SAM”의 네 경계가 완전히 분리된 상태는 아니다. checklist의 단순 `ready` 표시는 현재 구현과 어긋나며 **partially implemented**로 갱신하는 편이 정확하다.

## 4. 기존 제안의 비용과 우선순위

| 기존 항목 | 물리/계산 비용 | 오늘의 판단 |
|---|---|---|
| extra-view dependency correctness (`docs/checklist.json:515-532`) | key bookkeeping은 O(1), cache miss의 원자계 계산은 그대로 | 핵심 correctness와 atomic phasor 재사용은 반영됨. 남은 것은 cache layer 분리와 checklist 상태 동기화. **P1 마무리** |
| effective opacity 공개 (`docs/checklist.json:535-555`) | named scalar와 uncertainty sweep; 기본 계산 비용 거의 0, uncertainty sweep만 선택적 선형 증가 | default를 그대로 두면 spectrum은 bit/roundoff 수준으로 보존 가능. 절대 감도 신뢰도에 가장 큰 즉시 효과. **P1** |
| observability/readout guardrail (`docs/checklist.json:558-580`) | contrast/power/edge/sample 수는 O(Nscan); 두 단계 local refinement도 작은 추가 scan | 원자 물리를 바꾸지 않고 과신만 막는다. fixed reference와 함께 **최우선 P1** |
| held-out validation dataset (`docs/checklist.json:583-603`) | runtime보다 실험 획득·provenance 비용이 지배 | 내부 코딩으로 대체 불가. 데이터가 오기 전 결과는 PREDICTED 유지. **P2 / blocked_external** |
| Zeeman/polarization reference mode (`docs/checklist.json:914-936`) | Hilbert/Liouville 차원이 커지고 velocity·field grid와 곱해져 고비용 | default에 넣지 말고 특정 state/polarization/B dataset이 있을 때 opt-in. **parked가 맞음** |
| measured readout transfer (`docs/checklist.json:939-961`) | measured complex transfer를 기존 phasor에 곱하는 것은 O(Nf), 거의 negligible | full time-domain engine보다 싸고 물리를 보존한다. 측정 transfer 확보 후 **P2** |
| RF field-map import (`docs/checklist.json:964-982`) | O(Nspace×Nscan), mesh convergence와 큰 데이터 필요 | 내부 3-D EM solver는 범위 과대. 외부 field map이 생길 때 opt-in. **blocked_external** |

새로운 대형 물리 모델 제안은 필요하지 않다. 이미 기록된 P1 세 항목(opacity, guardrail/fixed reference, cache layering)이 현재 신뢰도와 사용성을 가장 싸게 높인다.

## 5. 물리 결과를 보존하는 순수 코드 최적화

### 5.1 off-Doppler Liouvillian batch 조립 — 가장 값싼 다음 최적화

`_superheterodyne_atomic_phasors`의 Doppler-off 경로는 801개 Hamiltonian에 대해 Python list comprehension으로 `core.liouvillian`을 호출한다 (`gabes/schemes/rydberg.py:808-813`). 같은 Hamiltonian은 scan detuning에 affine하므로, detuning-independent Liouvillian과 coefficient를 한 번 만들고 NumPy broadcasting으로 batch를 조립할 수 있다.

- 현재 list 조립 중앙값: 35.552 ms
- affine batch 조립 중앙값: 2.758 ms
- speedup: 12.89×
- 비교: `np.array_equal=True`, 최대 절대차 `0`

이는 해를 근사하거나 grid를 줄이지 않고 **동일한 행렬을 더 싸게 구성**한다. cache miss와 온도 sweep을 줄이는 최우선 순수 성능 변경이다.

### 5.2 detector scan vectorization

cache hit 뒤에도 `_superheterodyne_readout`은 최대 801개 점에서 `PhotodiodeChannel`, `BalancedDetector`, noise, sensitivity 객체를 Python loop로 만든다 (`gabes/schemes/rydberg.py:908-946`). power, shot noise, correlated RIN, electronics, field ASD를 배열식으로 계산하고 최적점 한 개에만 dataclass를 만들면 식과 단위를 유지하면서 22 ms post-processing을 줄일 수 있다. scalar reference fixtures와 배열 전체 `allclose`를 동시에 고정해야 한다.

### 5.3 heating view의 네 단계 cache 분리

현재 temperature sweep은 각 온도에서 EIT compute/readout과 AT compute/readout을 모두 다시 부른다 (`gabes/schemes/rydberg.py:1384-1442`). 다음 경계를 별도로 cache하면 detector, cell, SAM 조정에 원자 spectrum을 다시 풀 필요가 없다.

1. temperature/probe/coupling/dephasing → atomic raw spectrum
2. IF/LO/RF dephasing → atomic response phasor
3. cell/opacity → Beer–Lambert transmission
4. QE/path/reference/noise/RF dipole/SAM → detector와 field conversion

현재 detector-only 두 번째 sweep도 317 ms이므로 사용자 체감 이득이 있다. 저장할 데이터는 온도 10개 × 801점 complex spectrum 정도로 작고, cache size를 제한하면 메모리 부담도 낮다.

### 5.4 formatted metric 재파싱 제거

temperature sweep은 `_readout`의 표시 문자열에서 숫자를 `split()[0]`로 다시 읽는다 (`gabes/schemes/rydberg.py:1412-1423`). 비용은 작지만 locale/format 변화에 약하다. `_readout`에 이미 계산된 numeric slope/discriminator를 구조화 필드로 반환하고 renderer만 format하면 동일 물리를 더 안전하게 재사용할 수 있다.

## 6. 문서·테스트·예제 검토

### 잘 된 부분

- README의 Scheme 표와 Rydberg 설명은 finite-IF, detector budget, cold spot, SAM, conditional sensitivity, 의도적 deferred scope를 현재 코드와 맞게 설명한다 (`README.md:14`, `README.md:115-137`).
- 테스트는 defaults, linewidth/AT split, sensitivity 비주입, IF/temp, detector noise, field convention, atomic cache, extra-view dependency, effective atoms/SAM/heating/CSV를 폭넓게 덮는다 (`tests/test_rydberg_eit.py:31-56`, `tests/test_rydberg_eit.py:59-124`, `tests/test_rydberg_eit.py:281-427`, `tests/test_rydberg_eit.py:456-556`).
- 분석 workflow는 PREDICTED/ASSUMED/REFERENCE/PENDING를 구분하고 clean CSV 계약을 제공한다 (`analysis/rydberg_cell_heating/README.md:22-43`).

### 수정이 필요한 부분

1. HTML User Guide는 Rydberg를 정적 EIT/AT와 dispersion만 있는 것처럼 설명한다 (`docs/Userguide/GABES_User_Guide_v2.html:535-538`, `docs/Userguide/GABES_User_Guide_v2.html:617-625`). finite-IF, detector budget, temperature/cold spot, SAM, conditional status를 추가하고 `userguide_assets/rydberg.png`도 현재 3-panel AT view로 갱신해야 한다.
2. `tests/test_rydberg_eit.py:1-5`의 module docstring은 “public UI shows the static spectrum only”라고 해 현재 기능과 모순된다.
3. `test_weak_lo_reports_unresolved_at_status`는 unresolved AT에서도 `Total sensitivity` hero가 계속 존재해야 한다고 고정한다 (`tests/test_rydberg_eit.py:134-141`). 이것은 observability guardrail 목표와 반대다. guardrail 구현 시 numeric hero가 suppression/status로 바뀌는 테스트로 수정해야 한다.
4. high-temperature tiny-contrast, scan-edge optimum, fixed reference arm, successive refinement, public opacity provenance에 대한 회귀 fixture가 없다. 이 네 fixture를 P1 acceptance와 직접 연결해야 한다.
5. `analysis/rydberg_cell_heating/example_config.json:8-27`에는 문헌 reference 값(probe/coupling/cell/LO/IF/QE)과 실제 미측정 가정(path efficiency, reference ratio, zero electronics/RIN)이 한 `base_params`에 섞여 있다. dipole은 REFERENCE로 잘 분리되어 있지만 (`:92-96`), electrometry 전체는 ASSUMED 하나로만 묶인다 (`:98-102`). 각 파라미터 또는 논리 그룹에 provenance status를 붙이는 것이 재현성에 낫다.

## 7. 검증 결과

- Scheme 3 관련 6개 테스트 파일: **69 passed in 14.23 s**
  - `tests/test_rydberg_eit.py`
  - `tests/test_rydberg_electrometry.py`
  - `tests/test_rydberg_experiment.py`
  - `tests/test_rydberg_experimental_csv.py`
  - `analysis/rydberg_cell_heating/test_workflow.py`
  - `analysis/rydberg_cell_heating/test_adapter_integration_regressions.py`
- 제품 tests + Rydberg analysis: **595 passed in 307.17 s**
- 루트 `python -m pytest -q`: 기존 `references/AutoOD/test_autood_gui_adjustments.py` 수집 중 optional `sympy`가 없어 중단. Scheme 3 실패는 아니다.
- 이번 실행은 코드나 checklist를 고치지 않았고 이 검토 보고서만 추가했다. 기존 dirty worktree 변경은 보존했다.

## 8. 다음 행동 순서

1. **P1 — observability guardrail + fixed reference:** 저대비/저투과/edge/unresolved에서 numeric hero를 숨기고, reference arm의 balance detuning을 명시한다.
2. **P1 — effective opacity 공개:** `ls=0.001`을 provenance·범위·uncertainty가 있는 `effective_opacity_participation`으로 분리하고 default `≤1e-12` 보존 테스트를 둔다.
3. **P1 — affine Liouvillian batch + cache layering:** exact matrix fixture로 12.9× assembly 개선을 고정하고 detector-only sweep에서 atomic solve가 0회인지 검증한다.
4. **P1 — checklist/docs/tests 동기화:** extra-view 항목을 partially implemented로 갱신하고 User Guide·stale docstring·weak-LO hero 테스트를 현재 계약에 맞춘다.
5. **P2 — held-out data:** measured EIT/RF/PSD/temperature와 complex readout transfer가 들어오기 전 절대 감도·온도 optimum은 PREDICTED로 유지한다.

현재 가장 중요한 것은 full Zeeman이나 time-domain 엔진을 늘리는 일이 아니다. **이미 계산된 compact physics를 관측 가능성·provenance·고정 receiver 조건으로 정직하게 gate하는 것**이 가장 작은 비용으로 Scheme 3의 실험 reference 가치를 높인다.
