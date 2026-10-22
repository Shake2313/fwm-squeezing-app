# GABES 일일 Scheme 물리 검토 — 2026-09-12

## 0. 오늘의 선택과 검토 범위

- 로컬 날짜의 일(day) `12`에 대해 `n = (12 mod 5) + 1 = 3`이므로 오늘은
  **Scheme 3 — Rydberg-EIT electrometry**를 검토했다.
- 실제 dropdown registry는 `gabes/schemes/__init__.py:19-25`에 다음 순서로 정의돼 있다.
  1. `SASScheme` — OD / SAS
  2. `LambdaScheme` — EIT / AT / CPT
  3. `RydbergEITScheme` — Rydberg-EIT electrometry
  4. `MagnetoScheme` — Hanle / EIA / NMOR
  5. `FWMScheme` — seeded FWM / biphoton
- `gabes/schemes/__init__.py:27-39`의 Lambda/Magneto alias는 직접 조회용이며 위 다섯
  dropdown 순서를 바꾸지 않는다.
- 먼저 `CLAUDE.md`, `README.md`, `docs/checklist.json`, `docs/daily_report`의 Scheme 3
  보고서 12개(2026-06-22…2026-09-07), 코드의 TODO/FIXME, 분석 workflow와 예제,
  현재 테스트를 확인했다.
- 공개 GitHub issue는 [open](https://github.com/Shake2313/fwm-squeezing-app/issues)과
  [closed](https://github.com/Shake2313/fwm-squeezing-app/issues?q=is%3Aissue+state%3Aclosed)
  모두 0건이다. 저장소 안의 별도 TODO/issue 문서는 없고, 실제 개선 제안은
  `docs/checklist.json`과 기존 일일 보고서에 모여 있다.

## 1. 이전 Scheme 3 검토 이후의 변화

직전 검토 기준 `d98bdf4`에서 현재 `HEAD c0c46f0`까지 Scheme 3 원자계 코드, helper,
테스트, 분석 예제와 문서는 변하지 않았다. 해당 구간의 커밋은 checklist 진행 기록,
FWM trust 문구, 공통 앱 export 기능이며, Rydberg 관련 파일 diff는 0이다.

새 공통 기능인 결과 export는 Scheme 3의 실험 reference 가치에는 긍정적이다.
`streamlit_app.py:854-890`은 Git revision, scheme/cache/default/readout version, 전체
파라미터, raw 결과, metric/table, plot 배열과 experimental-overlay provenance를 JSON에
넣고, `streamlit_app.py:893-947`은 JSON 및 CSV/ZIP 다운로드를 제공한다. 다만 현재
테스트에서 이 helper들을 직접 호출하지 않으며, lazily rendered extra view와 named
scheme arrays, round-trip 재현은 아직 남아 있다. checklist도 정확히 `in_progress`로
기록한다(`docs/checklist.json:2498-2519`). 따라서 재현성은 개선됐지만 아직 완결된
experiment record 형식은 아니다.

## 2. 결론

**판정: 실제 Rydberg-EIT/AT와 유한-IF 검출 물리를 구현한 유용한 반정량적 실험 설계
도구다. 그러나 절대 감도, 온도 optimum, RF field calibration의 독립 reference로
인용하기에는 아직 부족하다.**

잘 맞는 용도는 다음과 같다.

- ⁸⁵Rb `5S–5P–40D` cascade EIT와 약 37 GHz `40D–39F` RF dressing의 line shape,
  linewidth, AT splitting과 parameter trend 탐색.
- LO-dressed 정상상태 주위의 유한-IF 약신호 응답, balanced detector의 PSN/RIN/
  electronics budget, RF dipole 변환을 일관된 단위계로 비교.
- heater setpoint, effective vapor temperature, cold spot을 구분한 heated-cell 시나리오,
  axial Beer–Lambert, effective atom number, SAM 불확도, unit-aware CSV 분석.
- 각 입력의 provenance가 별도로 확보된 뒤 “어떤 항목이 현재 감도를 제한하는가”를
  조건부로 비교.

다음 용도에는 아직 부적합하다.

- `6.08 nV/cm/√Hz` 같은 현재 기본 sensitivity를 장치의 절대 예측 또는 calibration
  curve로 사용하는 일.
- 내부 sweep의 “최적 온도”를 측정 검증 없이 실험 운전점으로 채택하는 일.
- Zeeman/polarization-resolved spectrum, horn/cell standing wave, nonlinear/time-domain
  lock-in response를 현재 compact model이 대표한다고 가정하는 일.
- 수치적으로 존재하는 FWHM/최적점을 낮은 대비에서도 관측 가능한 feature로 해석하는 일.

## 3. 현재 구현이 담고 있는 실제 물리

### 3.1 원자계와 정적 optical spectrum

- `RydbergEITScheme`은 ⁸⁵Rb `5S1/2 F=3 → 5P3/2 F'=4 → 40D5/2 → 39F7/2`
  4준위 cascade를 명시한다(`gabes/schemes/rydberg.py:1-15`, `72-91`, `417-437`).
- probe와 481 nm coupling power는 beam diameter와 함께 `Ω ∝ √P/d`로 Rabi frequency를
  구동한다(`gabes/schemes/rydberg.py:439-461`). 다만 절대 dipole 계산이 아니라
  문헌 operating point에 맞춘 effective anchor이므로 optical pumping/Zeeman/hyperfine
  누락을 흡수한 반경험적 calibration이다(`gabes/schemes/rydberg.py:41-55`).
- cold spot은 sealed-cell vapor pressure를, effective temperature는 local ideal-gas
  density·thermal motion·transit/dephasing를 정한다(`gabes/schemes/rydberg.py:463-505`).
  heater metadata만 바꿔 원자 solve가 달라지지 않는 것도 테스트한다
  (`tests/test_rydberg_eit.py:456-475`). 이 분리는 heated-cell 실험에 유용하고 물리적으로
  타당하다.
- 정적 spectrum은 801점 scan에 대해 affine Liouvillian kernel을 사용하고, 필요하면
  pointwise steady-state fallback으로 같은 susceptibility를 계산한다
  (`gabes/schemes/rydberg.py:507-527`, `544-556`).
- counter-propagating 780/481 nm geometry의 residual two-photon Doppler를 level별
  wavevector ratio로 담고, 기본은 논문에 맞춘 Doppler-off transit-limited 모드다
  (`gabes/schemes/rydberg.py:72-91`, `610-619`).
- uncompensated magnetic broadening은 full Zeeman manifold가 아니라 compensated
  susceptibility의 Gaussian convolution proxy다(`gabes/schemes/rydberg.py:621-627`).
  B-compensation trend에는 유용하지만 polarization/B-dependent absolute spectrum은 아니다.
- Beer–Lambert 변환 자체는 표준적이지만 effective line-strength/opacity 역할의
  `ls=0.001`이 코드에 고정돼 있다(`gabes/schemes/rydberg.py:629-644`, `672-681`).
  이것이 절대 transmission과 sensitivity의 가장 큰 공개되지 않은 자유도다.

### 3.2 finite-IF electrometry와 detector chain

- 약한 SIG는 LO-dressed Liouvillian 주위에서 `L₀+iω_IF I` 선형계를 풀어 sideband
  density matrix를 얻는다(`gabes/rydberg_electrometry.py:116-194`). 이는 실제
  frequency-domain linear response이며 단순 정적 slope 재명명이 아니다.
- Hermitian RF perturbation을 이용해 minus sideband 하나만 풀고 plus sideband를
  conjugate transpose로 복원한다(`gabes/rydberg_electrometry.py:173-185`). 물리와
  Hermiticity를 보존하는 올바른 최적화이며 `tests/test_rydberg_electrometry.py:39-94`가
  식과 solve 횟수를 검증한다.
- Doppler-on 경로는 complex response phasor를 velocity weight로 평균한다
  (`gabes/schemes/rydberg.py:787-806`, `gabes/rydberg_electrometry.py:211-227`). magnitude를
  먼저 평균하지 않는 점은 중요하다. 다만 이 모드는 full detuning scan이 아니라 정적
  최대-slope 한 점만 푼다.
- RF dipole, angular/polarization factor, RMS/peak convention이 분리돼 있고 AT splitting의
  detuning 제거와 field 변환도 명시돼 있다(`gabes/rydberg_electrometry.py:231-316`).
- detector는 one-sided Schottky ASD, correlated RIN, electronics noise를 합성한 뒤 field
  responsivity로 나눈다(`gabes/rydberg_electrometry.py:403-483`). 문헌 sensitivity
  상수를 결과에 주입하지 않는 구조는 좋다.
- 문제는 각 detuning에서 reference optical power를 그 지점의 transmitted signal power에
  비례시켜 다시 만든다는 점이다(`gabes/schemes/rydberg.py:902-927`). 실제 fixed arm 또는
  한 operating detuning에서 맞춘 balance가 아니라 pointwise-balanced envelope이므로
  sensitivity scan을 낙관적으로 만들 수 있다.
- low transmission/contrast, scan-edge optimum, under-resolved doublet에 대한 통합 validity
  gate가 없다. EIT는 `contrast > 0`이면 FWHM을 내고(`gabes/schemes/rydberg.py:979-1003`),
  weak-LO AT에서도 `Total sensitivity`가 hero로 남는다
  (`gabes/schemes/rydberg.py:1127-1207`, `tests/test_rydberg_eit.py:134-141`).

### 3.3 실험 분석 helper와 예제

- `gabes/rydberg_experiment.py:53-95`, `171-307`, `368-448`은 온도 상태, axial
  profile/segmented Beer–Lambert, Gaussian-overlap atom number를 각각 분리한다.
- `gabes/rydberg_experiment.py:485-594`는 SAM field와 first-order standard uncertainty,
  `2D²/λ` far-field 경고를 제공한다. 내부 3-D horn/cell EM solver라고 과장하지 않는다.
- `gabes/rydberg_experimental_csv.py:35-100`, `441-691`은 EIT/RF/PSD trace의 단위,
  중복, non-finite row, SHA-256 provenance와 unit assumption을 보존한다.
- 분석 workflow는 `MEASURED/PPT/REFERENCE/FITTED/PREDICTED/ASSUMED/PENDING`을 분리하고
  (`analysis/rydberg_cell_heating/README.md:22-39`), axial model의 한계도 명시한다
  (`analysis/rydberg_cell_heating/README.md:78-117`). 이는 실험 reference code로서 좋은
  설계다.
- 그러나 `analysis/rydberg_cell_heating/example_config.json:8-27`의 한 `base_params` 안에는
  문헌 reference 조건과 미측정 detector/path/reference/noise 가정이 섞여 있다. raw EIT,
  sensitivity, temperature log도 모두 `PENDING`이다(`:53-78`). dipole은 별도
  `REFERENCE`지만(`:92-96`) 전체 electrometry는 하나의 `ASSUMED` status다(`:98-102`).
  예제는 pipeline demo이지 held-out validation 예제가 아니다.

## 4. 문헌 기준과 오늘의 수치 재검증

Ju et al.의 [원 논문](https://arxiv.org/abs/2606.04354)은 residual magnetic field와 laser
noise를 줄여 1.6 MHz EIT linewidth를 얻고, 37 GHz에서 측정 감도
`12.5(8) nV/cm/√Hz`와 계산 PSN limit `11.2 nV/cm/√Hz`를 보고한다. 현재 defaults는
probe 6 µW, coupling 30 mW, 0.15 mm beam, 50 mm cell, LO 3.7 MHz, IF 40 kHz를
고정한다(`gabes/schemes/rydberg.py:110-140`, `tests/test_rydberg_eit.py:31-42`).

현재 worktree에서 직접 계산한 결과는 다음과 같다.

| 점검 | 결과 | 해석 |
|---|---:|---|
| 기본 EIT 중심 transmission | 0.940016615 | 기본점은 충분히 관측 가능 |
| 기본 EIT contrast | 0.135391492 | 정상 범위 |
| 기본 EIT FWHM | 1.613058583 MHz | 논문의 1.6 MHz scale과 일치 |
| 기본 AT 중심 transmission | 0.806029283 | detector power가 남는 조건 |
| 기본 AT split | 3.496500 MHz | LO Rabi 3.7 MHz와 같은 scale |
| finite-IF optimum | −2.220000 MHz | 기본 scan 내부 |
| 기본 PSN = total sensitivity | 6.082387215 nV/cm/√Hz | 문헌 PSN보다 1.84배 낙관적 |
| `ls=0.0005` total sensitivity | 11.781990083 nV/cm/√Hz | opacity 영향이 큼 |
| `ls=0.001` total sensitivity | 6.082387215 nV/cm/√Hz | 현재 hard-coded 기본값 |
| `ls=0.002` total sensitivity | 3.238775723 nV/cm/√Hz | 4× opacity 범위에서 3.64× 변화 |
| 80 °C EIT transmission/contrast | 1.798922847×10⁻⁶ | 사실상 관측 불가 |
| 80 °C EIT FWHM | 0.307015017 MHz | 정밀 숫자가 계속 hero로 노출 |
| 80 °C static IF optimum | −8.96 MHz | ±9 MHz scan edge에서 0.04 MHz |

논문과 비슷한 기본 linewidth는 유용한 sanity check지만, probe/coupling Rabi와 opacity가
reference point에 맞춘 effective 값이므로 독립 예측 검증은 아니다. 특히 hard-coded opacity
하나가 sensitivity를 수배 바꾸고, 현재 기본 PSN이 논문보다 이미 1.84배 좋다는 사실은
`6.08`을 절대 reference로 사용하지 말아야 한다는 직접 증거다.

## 5. 기존 개선 제안과 계산비용 평가

기존 제안은 충분하며 새로운 대형 physics proposal은 필요하지 않다.

| 기존 checklist 항목 | 계산비용과 물리 보존성 | 오늘의 판단 |
|---|---|---|
| extra-view dependency correctness (`docs/checklist.json:1483-1500`, P1 ready) | key bookkeeping은 O(1). IF가 atomic phasor만, detector/cell/RF conversion이 후처리만 무효화되면 원자 물리는 동일하다. 현재 `ExtraView.param_keys`와 bounded phasor cache는 반영됐지만 네 단계 cache 경계가 완결되지는 않았다. | **partially implemented**로 상태를 동기화하고 cache layering을 마무리할 가치가 큼 |
| effective opacity / sweep trust (`docs/checklist.json:1503-1523`, P1 ready) | hard-coded scalar를 named input으로 옮기는 기본 경로 비용은 사실상 0. uncertainty sweep만 선택적으로 샘플 수에 선형 증가한다. default `0.001` 회귀를 두면 현재 spectrum을 `≤1e-12`로 보존 가능하다. | 절대값 신뢰도에 가장 큰 즉시 효과. **P1** |
| readout validity guardrails (`docs/checklist.json:1526-1548`, P1 needs_decision) | contrast, power, samples/FWHM, edge distance는 기존 801점 배열의 O(N) 후처리라 추가 solve가 없다. 두 번의 deterministic local refinement만 작은 선택적 비용이다. fixed reference arm은 원자 physics를 바꾸지 않고 detector model을 실제 장치에 맞춘다. | **최우선 P1**. threshold 정책 결정 필요 |
| held-out validation dataset (`docs/checklist.json:1551-1572`, P2 blocked_external) | runtime보다 측정·metadata·권한·hash 관리 비용이 지배한다. 내부 계산으로 대체할 수 없다. | 데이터가 오기 전 absolute/optimum은 계속 `PREDICTED` |
| Zeeman/polarization reference mode (`docs/checklist.json:1985-2007`, P2 parked) | Hilbert/Liouville 차원과 velocity/B grid가 함께 커져 기본 4-level solve보다 매우 비싸다. 선택한 state/polarization/B dataset 없이는 자유도만 늘어난다. | opt-in research mode로 parked가 맞음 |
| measured transfer function (`docs/checklist.json:2010-2032`, P2 blocked_external) | 측정 complex transfer를 현재 phasor에 곱하는 것은 O(Nf)로 거의 무시 가능. full time-domain solver보다 훨씬 싸며 linear response를 보존한다. | 측정 transfer 확보 후 우선 적용 |
| RF field-map import (`docs/checklist.json:2035-2058`, P2 blocked_external) | O(Nspace×Nscan)와 mesh convergence 비용. 내부 3-D EM보다 외부 map import가 타당하다. | named geometry/map이 생길 때 opt-in |

## 6. 동작을 바꾸지 않는 순수 코드 최적화

### 6.1 Doppler-off Liouvillian affine batch 조립

`_superheterodyne_atomic_phasors`는 801개 detuning마다 Python에서
`core.build_liouvillian`을 호출한다(`gabes/schemes/rydberg.py:808-815`). Hamiltonian과
Liouvillian은 probe detuning에 affine이므로 기준 행렬과 단위-detuning coefficient를 한 번
만든 뒤 broadcasting할 수 있다.

오늘 같은 process에서 직접 비교한 중앙값은 다음과 같다.

- 현재 list assembly: `79.819 ms`
- affine batch assembly: `3.111 ms`
- speedup: `25.66×`
- 검증: `np.array_equal=True`, 최대 절대차 `0`

grid 축소나 근사가 아니라 **동일한 행렬을 더 싸게 만드는 변경**이다. exact-array fixture를
추가하면 가장 안전한 첫 최적화다.

### 6.2 detector scan vectorization

atomic phasor cache hit 뒤에도 801개 점 각각에서 `PhotodiodeChannel`, `BalancedDetector`,
noise, sensitivity dataclass를 생성한다(`gabes/schemes/rydberg.py:908-946`). signal/reference
photocurrent, Schottky variance, correlated RIN, electronics, field ASD는 모두 배열식이다
(`gabes/rydberg_electrometry.py:403-429`, `463-483`). 전체 후보는 NumPy로 계산하고 최종
optimum 한 점에만 dataclass를 만들 수 있다. scalar helper 결과와 배열 전체를 동시에
`allclose`로 고정하면 동작을 바꾸지 않는다.

### 6.3 heating extra view의 네 계층 cache

온도 10점마다 EIT와 AT `compute/_readout`을 모두 반복한다
(`gabes/schemes/rydberg.py:1384-1465`). 오늘의 3회 중앙값에서는 기본 sweep과 detector-only
변경 sweep이 각각 약 `675 ms`, `759 ms`로 둘 다 수백 ms였다. detector noise만 바뀌어도
정적 EIT/AT spectrum을 다시 만드는 구조임을 보여 준다. 다음을 별도 bounded cache로 나누면
원자 physics는 그대로 유지된다.

1. temperature/probe/coupling/dephasing → atomic raw spectrum
2. IF/LO/RF dephasing → finite-IF atomic phasor
3. cell/opacity → Beer–Lambert transmission
4. QE/path/fixed reference/noise/RF dipole/SAM → detector와 field conversion

### 6.4 formatted metric 재파싱 제거

temperature sweep은 `_readout`이 표시용 문자열로 만든 slope/discriminator를 다시
`split()[0]`으로 숫자로 읽는다(`gabes/schemes/rydberg.py:1412-1423`). 큰 병목은 아니지만
locale/format 변화에 취약하다. numeric readout 구조를 먼저 반환하고 renderer만 format하면
동일 물리를 더 안전하고 조금 더 빠르게 재사용할 수 있다.

## 7. 문서·테스트·예제의 reference 적합성

### 좋은 점

- README는 finite-IF, detector budget, cold spot, SAM, conditional sensitivity와 deferred
  scope를 현재 코드에 맞게 설명한다(`README.md:115-137`).
- Scheme 테스트 31개와 helper/workflow 테스트는 defaults, line width, AT splitting,
  IF response, Hermiticity, RMS/peak, detector noise, cache, temperature/cold spot, axial/SAM,
  CSV provenance를 폭넓게 다룬다(`tests/test_rydberg_eit.py:31-540` 및 관련 helper tests).
- `test_superhet_low_if_matches_static_lo_derivative`는 finite-IF 구현의 저주파 한계를 독립
  finite difference와 비교한다(`tests/test_rydberg_eit.py:309-337`). 단순 숫자 golden보다
  강한 물리 검증이다.
- 분석 workflow는 measurement와 prediction의 evidence class를 구분하고, clean CSV와
  SHA-256 source manifest를 제공한다.
- 새 app export는 같은 render에 사용한 파라미터와 raw/plot/provenance를 외부로 가져갈 수
  있게 해 재현성을 높였다.

### 남은 문제

1. HTML User Guide는 Rydberg를 정적 EIT/AT와 dispersion만 있는 것처럼 설명한다
   (`docs/Userguide/GABES_User_Guide_v2.html:535-538`, `617-625`). finite-IF, detector
   budget, temperature/cold spot, SAM, conditional status와 현재 3-panel figure가 빠져 있다.
2. `tests/test_rydberg_eit.py:1-5`는 여전히 “public UI shows the static spectrum only”라고 해
   현재 공개 sensitivity와 모순된다.
3. weak-LO test는 unresolved AT에서도 `Total sensitivity` hero가 유지돼야 한다고 고정한다
   (`tests/test_rydberg_eit.py:134-141`). observability guardrail을 구현할 때 suppression/status
   계약으로 바꿔야 한다.
4. high-temperature tiny-contrast, scan-edge rejection, fixed reference arm, sub-grid/local
   refinement, public opacity provenance를 acceptance test가 아직 직접 고정하지 않는다.
5. 새 app export helper는 테스트에서 직접 참조되지 않는다. round-trip, complex raw array,
   CSV unit header, Rydberg experimental overlay와 lazy extra-view export를 검증해야 한다.
6. example config는 좋은 실행 예제지만 입력별 provenance와 held-out 측정이 없으므로
   experimental validation example로 오해하지 않도록 유지해야 한다.

## 8. 검증 결과

- Scheme 3 관련 6개 파일: **69 passed in 29.68 s**
  - `tests/test_rydberg_eit.py`
  - `tests/test_rydberg_electrometry.py`
  - `tests/test_rydberg_experiment.py`
  - `tests/test_rydberg_experimental_csv.py`
  - `analysis/rydberg_cell_heating/test_workflow.py`
  - `analysis/rydberg_cell_heating/test_adapter_integration_regressions.py`
- 필수 루트 `python -m pytest -q`: **collection 실패**. 60.70 s 뒤
  `references/AutoOD/test_autood_gui_adjustments.py`가 선언되지 않은 optional `sympy`를
  import해 `ModuleNotFoundError`가 발생했다. Scheme 3 테스트 실패는 아니다.
- 기존 dirty/staged 사용자 변경은 건드리지 않았다. production code, checklist, docs/example
  source는 수정하지 않았고 오늘의 보고서만 추가했다.

## 9. 권장 실행 순서

1. **P1 — observability/fixed-reference 정책 결정 및 구현:** low contrast/power, edge,
   unresolved, unconverged 결과에서 numeric hero를 억제하고 reference balance detuning을
   명시한다.
2. **P1 — effective opacity 공개:** `ls=0.001`을 provenance·범위·불확도가 있는 named
   parameter로 분리하고 default 보존 회귀를 둔다.
3. **P1 — exact affine batch + cache layering:** array-exact Liouvillian fixture와
   detector-only 변경 시 atomic solve 0회 회귀를 추가한다.
4. **P1/P2 — 상태·문서·테스트 동기화:** extra-view 항목을 partially implemented로,
   User Guide/docstring/weak-LO test/export round-trip을 현재 계약에 맞춘다.
5. **P2 — held-out measurement:** EIT/RF/PSD/temperature/complex transfer 데이터가 오기
   전 절대 sensitivity와 heating optimum은 `PREDICTED`로 유지한다.

현재 Scheme 3의 가장 큰 약점은 4준위 OBE가 가짜라는 것이 아니라, **실제 compact physics의
출력을 관측 가능성·opacity provenance·고정 receiver 조건으로 충분히 gate하지 않는 것**이다.
이 세 가지는 full Zeeman이나 3-D EM보다 훨씬 싸고, 실험자가 잘못 확신할 위험을 더 크게 줄인다.
