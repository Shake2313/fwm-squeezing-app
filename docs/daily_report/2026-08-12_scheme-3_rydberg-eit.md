# 2026-08-12 Scheme 3 물리 검토 — Rydberg-EIT electrometry

## 1. 선택 규칙과 현재 다섯 스킴

- 현지 날짜: 2026-08-12 (Asia/Seoul), 일자 `12`
- `n = (12 mod 5) + 1 = 3`
- 오늘의 대상: **Scheme 3 — Rydberg-EIT electrometry**

드롭다운의 실제 순서는 `gabes/schemes/__init__.py:19-24`의 `_SCHEMES`가 결정한다.
README의 기능 정의(`README.md:11-17`)와 일치한다.

| 번호 | 등록 클래스 | 현재 사용자-facing 정의 |
|---:|---|---|
| 1 | `SASScheme()` | 흡수 분광 OD / SAS |
| 2 | `LambdaScheme()` | Λ 결맞음 EIT / AT / CPT |
| **3** | `RydbergEITScheme()` | **Rydberg-EIT / microwave AT electrometry** |
| 4 | `MagnetoScheme()` | Hanle / EIA / NMOR |
| 5 | `FWMScheme()` | seeded FWM squeezing / generic SFWM biphoton |

## 2. 먼저 확인한 기존 제안

새 제안부터 만들 필요는 없었다. 2026-07-22, 08-02, 08-07 Scheme 3 보고서와
현재 체크리스트에 이미 구체적인 개선안이 있다.

- **P1 ready — extra-view cache correctness:** navigate-only IF, QE, path, reference arm,
  cell length, RF dipole가 extra-view의 실제 계산에는 들어가지만 app cache key에서는 빠진다
  (`docs/checklist.json:464-480`).
- **P1 ready — effective opacity와 sweep trust:** 숨은 `ls=0.001`을 명시적 calibration
  parameter로 만들고 온도/파워 sweep을 검증 전 `PREDICTED`로 한정한다
  (`docs/checklist.json:484-500`).
- **P1 needs_decision — readout validity:** 저투과·저대비·경계·저해상도 조건의 hero 수치를
  막고, fixed reference arm과 peak/optimum interpolation을 요구한다
  (`docs/checklist.json:507-528`).
- **P2 blocked_external — held-out dataset:** EIT/RF/PSD/온도 측정과 provenance가 없어
  absolute validation이 막혀 있다 (`docs/checklist.json:532-552`).
- full Zeeman/polarization reference mode는 구체적 상태·편광·B-field·검증 스펙트럼이 생길
  때까지 `parked`다 (`docs/checklist.json:850-868`). measured transfer function과 RF field map도
  각각 외부 데이터가 필요한 `blocked_external` 항목이다 (`docs/checklist.json:875-919`).

2026-08-07 보고서는 온도 최적점 47 °C 대 논문 약 33 °C, probe-power 최적점
4 µW 대 논문 6 µW를 이미 지적했고(`docs/daily_report/2026-08-07_scheme-3_rydberg-eit.md:122-159`),
숨은 opacity·cache·저투과 hero·fixed reference·grid 보간·ARC provenance와 순수 최적화까지
정리했다(`:165-282`). 오늘은 이 제안들이 현재 코드에서도 유효한지 재검증했다.

## 3. 현재 구현이 담고 있는 실제 물리

### 3.1 축약 4준위 cascade OBE

원자 코어는 ⁸⁵Rb
`5S1/2 F=3 → 5P3/2 F'=4 → 40D5/2 → 39F7/2`의 4준위 cascade다.
780 nm probe, 481 nm counter-propagating coupling, 약 37 GHz RF LO를 Hamiltonian에 넣고,
5P spontaneous decay, Rydberg decay, 5S–40D와 40D–39F dephasing을 분리한다
(`gabes/schemes/rydberg.py:4-15`, `69-88`, `413-433`, `554-623`).

- probe/coupling Rabi는 `Ω ∝ √P/d`로 reference operating point에 anchor된다
  (`gabes/schemes/rydberg.py:435-457`).
- effective temperature는 thermal velocity와 transit broadening에, cold spot은 sealed-cell
  vapor pressure에 들어간다 (`gabes/schemes/rydberg.py:459-501`, `559-580`).
- residual two-photon Doppler는 per-level wave-vector ratio와 Maxwell velocity average로
  표현된다 (`gabes/schemes/rydberg.py:69-88`, `503-523`, `606-615`).
- susceptibility는 absolute density/dipole와 Beer–Lambert 전파로 transmission에 매핑된다
  (`gabes/schemes/rydberg.py:625-677`).

이는 기준점 주변에서 EIT linewidth, power/waist/transit trade-off, RF AT splitting과 detuning
방향을 탐색하는 데 실제로 유용하다. 다만 probe/coupling Rabi가 bare dipole에서 독립적으로
계산된 값이 아니라 fitted anchor이고(`gabes/schemes/rydberg.py:39-53`), residual Zeeman은
full manifold가 아니라 Gaussian convolution이다(`gabes/schemes/rydberg.py:617-623`).

### 3.2 finite-IF weak-SIG와 detector chain

정적 slope proxy보다 강한 실제 선형응답이 구현돼 있다.

- LO-dressed steady state 주위에서 `L₀ ± iω_IF` sideband를 풀며, 응답의 복소 위상과
  Hermitian observable convention이 명시돼 있다
  (`gabes/rydberg_electrometry.py:1-35`, `116-199`).
- velocity-class phasor는 magnitude가 아니라 complex coherent average된다
  (`gabes/rydberg_electrometry.py:69-113`, `216-232`).
- RF dipole, angular/polarization factor, peak/RMS convention을 전기장과 Rabi 사이에 둔다
  (`gabes/rydberg_electrometry.py:236-321`).
- balanced detector는 독립 photodiode shot noise, correlated RIN, electronics ASD를 합성해
  field ASD로 환산한다 (`gabes/rydberg_electrometry.py:324-488`).
- Scheme은 801점 probe scan에서 total noise-equivalent field가 최소인 점을 고른다
  (`gabes/schemes/rydberg.py:769-923`).

따라서 finite-IF primitive는 수학적으로나 물리적으로 유용하다. 그러나 이는 full LO+SIG
time trace, mixer/lock-in transfer, sampling, colored technical noise를 포함한 실험 acquisition
모델은 아니다. 코드도 그 경계를 명시한다 (`gabes/rydberg_electrometry.py:3-7`).

### 3.3 heated-cell, axial, SAM, CSV 계층

실험 reference로서 특히 좋은 부분이다.

- setpoint, sensor, effective vapor temperature, cold spot을 구분하고 source 문자열을 보존한다
  (`gabes/rydberg_experiment.py:47-162`).
- cold-spot-limited `n(z)`, column density, segmented Beer–Lambert를 별도 primitive로 둔다
  (`gabes/rydberg_experiment.py:165-338`).
- Gaussian optical overlap의 effective atom number를 participation과 geometry factor로 분리한다
  (`gabes/rydberg_experiment.py:341-443`).
- SAM은 `E_RMS=√(30PG)/r`, 표준불확도, `r/(2D²/λ)` 경고를 구현한다
  (`gabes/rydberg_experiment.py:447-589`).
- EIT/RF/PSD CSV는 unit assumption, ignored/duplicate row, SHA-256을 보존하고 detector signal을
  임의 정규화하지 않는다 (`gabes/rydberg_experimental_csv.py:1-7`, `211-379`, `446-718`).

분석 workflow는 `MEASURED/PPT/REFERENCE/FITTED/PREDICTED/ASSUMED/PENDING`을 분리하며
(`analysis/rydberg_cell_heating/README.md:22-39`), axial 결과가 density-rescaled lumped line shape이고
spatial finite-IF solve가 아님을 명시한다. 이는 실험 notebook/report 인프라로 적절하다.

## 4. 실험물리 reference 적합성 판정

**판정: 실제 물리를 구현한 유용한 semi-quantitative experiment-planning/reference code다.
그러나 absolute Rydberg electrometry나 cell-heating optimum의 계측 기준으로는 아직 부적합하다.**

[Ju et al.](https://arxiv.org/abs/2606.04354)은 1.6 MHz EIT linewidth와 37 GHz에서
12.5(8) nV/cm/√Hz를 보고하며, 계산 PSN limit은 20 °C·6 µW에서 11.2 nV/cm/√Hz,
온도 optimum은 약 33 °C라고 기술한다. 현재 GABES는 이 실험의 topology와 기준 knob를
잘 반영하지만, linewidth/Rabi anchor와 숨은 opacity가 같은 논문 근처에 맞춰져 있어
독립적인 hold-out 검증은 아니다.

가장 중요한 trust boundary는 full Zeeman 부재 자체보다 먼저 다음 세 가지다.

1. absolute transmission과 sensitivity를 지배하는 `ls=0.001`이 calibration parameter로
   노출되지 않는다 (`gabes/schemes/rydberg.py:625-643`).
2. app extra-view가 사용자의 detector/IF/cell 입력 일부를 실제 계산에 전달하지 않는다
   (`streamlit_app.py:617-622`, `1555-1558`, `1607-1613`).
3. low-contrast/scan-edge 결과도 precise hero가 될 수 있다
   (`gabes/schemes/rydberg.py:926-1051`, `1137-1148`).

## 5. 오늘의 수치 재검증

수치는 JIT warm 상태의 현재 working tree에서 계산했다. 시간은 이 Windows 환경의 median으로,
절대 성능보다 병목 비율을 보는 값이다.

| 조건 | 현재 결과 |
|---|---|
| EIT default | `T(0)=0.9400166`, contrast `0.1353915`, FWHM `1.6130586 MHz`, max slope `0.122/MHz` |
| AT default | `T(0)=0.8060293`, split `3.4965 MHz`, optimum `-2.2200 MHz` |
| AT sensitivity | PSN = total = `6.082387 nV/cm/√Hz` (default RIN/electronics = 0) |
| warm timing | EIT static `20.1 ms`, AT static `4.27 ms`, AT readout `74.8 ms` |

### 5.1 숨은 opacity 의존성

동일 atomic response에서 `raw["ls"]`만 바꿨다. solver state는 그대로이고 optical propagation과
field sensitivity만 재계산했다.

| effective `ls` | AT `T(0)` | PSN/total sensitivity |
|---:|---:|---:|
| 0.0005 | 0.897791 | 11.78199 nV/cm/√Hz |
| **0.0010** | **0.806029** | **6.08239 nV/cm/√Hz** |
| 0.0020 | 0.649683 | 3.23878 nV/cm/√Hz |

factor 2의 의미 선택이 sensitivity를 약 1.9배 움직인다. 같은 저장소의 실제 D2 atomic line
strength와 구분되는 effective opacity/participation calibration으로 노출하기 전에는 default
감도를 절대 atomic prediction으로 인용하면 안 된다.

### 5.2 외부 optimum과 불일치

- 내장 10점 온도 sweep: best = **47 °C**, `0.971731 nV/cm/√Hz`; 20 °C는 `6.082387`.
- 논문: 계산 optimum 약 **33 °C**, 20 °C·6 µW PSN limit **11.2**.
- 현재 probe-power sensitivity sweep: best = **4 µW**, `5.911989`; 논문 optimum은 **6 µW**.

이는 논문 값을 hard-code하라는 뜻이 아니다. 반대로 현재 sweep를 `PREDICTED/UNVALIDATED`로
표시하고 opacity·dephasing·detector uncertainty band를 함께 내야 한다. 현재 example의 raw EIT,
RF, temperature inputs도 모두 `PENDING`이다
(`analysis/rydberg_cell_heating/example_config.json:53-78`).

### 5.3 low-contrast/edge hero 재현

80 °C EIT에서 `T(0)=contrast=1.799×10⁻⁶`인데도 FWHM `0.307015 MHz`와
`Transmission at resonance`가 hero다. IF optimum은 `-8.96 MHz`, scan edge는 `-9 MHz`다.
이는 중앙 feature의 실험 관측 가능성보다 scan background/edge를 정밀 숫자로 승격하는 사례다.

### 5.4 extra-view cache 불일치 재현

`IF=200 kHz`, `QE=0.2`, path efficiency `0.2`, cell `20 mm`, microwave metadata `55 GHz`를
선택한 뒤 app과 같은 recompute-only dictionary로 temperature extra-view를 호출했다.

| 경로 | 20 °C sensitivity | 선택된 best T |
|---|---:|---:|
| app cache-equivalent recompute-only | `6.08239` | 47 °C |
| direct full-parameter | `52.18463` | 55 °C |

현재 `ExtraView`에는 dependency declaration이 없고(`gabes/schemes/base.py:78-84`), app은
recompute key만 전달한다. 전용 test는 full params를 직접 전달하므로 이 UI 경로를 잡지 못한다
(`tests/test_rydberg_eit.py:419-431`). 또한 temperature sweep은 사용자가 낮은 LO를 선택해도
`max(user_LO, 3.7 MHz)`로 강제한다 (`gabes/schemes/rydberg.py:1374-1379`). 이 override도
입력 semantics와 status에 명시하거나 제거해야 한다.

### 5.5 grid 좌표 의존성

| scan points | Δx | AT split | sensitivity optimum | total sensitivity |
|---:|---:|---:|---:|---:|
| 401 | 0.05550 MHz | 3.4410 MHz | -2.2200 MHz | 6.08239 |
| 801 | 0.02775 MHz | 3.4965 MHz | -2.2200 MHz | 6.08239 |
| 1601 | 0.013875 MHz | 3.4965 MHz | -2.206125 MHz | 6.07873 |

default는 대체로 안정적이지만 peak와 optimum은 여전히 grid 좌표다
(`gabes/schemes/rydberg.py:662-666`, `894-906`, `994-1005`). 보간과 successive-refinement status가
없는 현재 표시 정밀도는 수치 격자보다 한 단계 과하다.

## 6. 기존 개선안의 비용과 물리 보존성

| 기존 항목 | 현재 상태 | 추가 계산비용 | 물리 보존 평가 |
|---|---|---:|---|
| ExtraView dependency/caches 분리 | P1 `ready`, medium | key 작성 O(k); detector-only 변경은 atomic re-solve를 제거 | 같은 식을 올바른 입력으로 계산하므로 결과 semantics를 복원한다. atomic steady state, finite-IF phasor, Beer–Lambert, detector/SAM cache 경계를 분리해야 한다. |
| `ls=0.001` 공개 + uncertainty sweep | P1 `ready`, medium | 현재 모델에서는 O(Nscan) transmission/noise remap, OBE solve 없음 | atomic line strength와 effective participation을 분리하면 현재 default를 1e-12 수준으로 보존할 수 있다. 향후 optical pumping에 넣을 때만 새 physics solve가 필요하다. |
| contrast/power/edge/resolution guardrail | P1 `needs_decision`, medium | O(Nscan), 새 OBE solve 없음 | 기존 배열에서 status와 local fit을 내므로 물리를 잃지 않는다. threshold policy는 detector calibration 유무에 따라 결정해야 한다. |
| fixed reference arm + peak/optimum interpolation | 위 P1에 통합 | detector loop O(Nscan), interpolation O(1) | 실제 한 reference arm의 의미를 복원한다. atomic response는 그대로다. |
| ARC regeneration/provenance fixture | 위 P1에 통합 | 평상시 O(1), optional validation script만 ARC 실행 | field scale의 state/polarization/version을 재현 가능하게 만든다. ARC API는 matrix element를 `e a₀` 단위로 정의한다([ARC 문서](https://arc-alkali-rydberg-calculator.readthedocs.io/en/latest/generated/arc.alkali_atom_functions.AlkaliAtom.getDipoleMatrixElementHFS.html)). |
| held-out EIT/RF/PSD dataset | P2 `blocked_external`, large | batch fit/validation은 중간; 핵심 비용은 측정·provenance | absolute reference에 필수다. fit/holdout를 사전 분리해야 calibration circularity를 피한다. |
| measured transfer function | P2 `blocked_external`, medium | 선형 complex filter O(Nf), atomic solve 없음 | full time-domain보다 훨씬 싸고, 측정된 필요가 있을 때만 확장하는 현재 checklist 판단이 적절하다. |
| full Zeeman/polarization reference | P2 `parked`, research | 상태수 증가로 dense Liouvillian 비용 급증 | absolute polarization/B dependence에는 중요하지만 negligible-overhead가 아니다. compact model 신뢰 경계를 먼저 고치는 순서가 맞다. |
| RF field-map import | P2 `blocked_external`, large | O(Nspace·Nscan) 또는 interpolation/cached reduction | uniform-map limit을 scalar SAM과 맞추는 별도 reference path가 적절하다. 내부 3-D horn solver를 만드는 것은 범위를 벗어난다. |

## 7. 동작을 바꾸지 않는 순수 코드 최적화

### 7.1 finite-IF Liouvillian affine assembly — 최우선

현재 801 detuning마다 Python에서 Hamiltonian과 Liouvillian을 새로 만든다
(`gabes/schemes/rydberg.py:811-818`). detuning dependence는 정확히 affine이므로
`L(s)=L(0)+s[L(1)-L(0)]`을 broadcasting할 수 있다.

- 현재 list assembly median: `34.254 ms`
- affine broadcast median: `2.347 ms`
- **14.60× faster**, `np.array_equal=True`, max absolute difference `0`

solve 방정식이나 부동소수 결과를 바꾸지 않는 가장 확실한 최적화다.

### 7.2 이미 구한 steady state 재사용

동일 Liouvillian batch에 대해 finite-IF helper가 `rho₀`를 다시 구한다
(`gabes/rydberg_electrometry.py:165-188`). optional `steady_state` 경로를 사용한 오늘의 진단은:

- fresh response: `14.093 ms`
- supplied `rho₀`: `9.542 ms`
- **1.48× faster**, `rho_minus/rho_plus` bit-identical

static affine kernel에서 필요한 `rho₀`를 선택적으로 반환하거나 phasor cache layer가 보존하면 된다.
메모리 증가를 피하려면 full scan density matrix를 영구 보관하지 말고 finite-IF 계산에 필요한
batch lifetime으로 제한하는 편이 안전하다.

### 7.3 atomic phasor와 detector/SAM postprocess cache 분리

QE, optical path, reference ratio, RIN, electronics, ENBW, RF dipole, SAM은 LO-dressed atomic
phasor를 바꾸지 않는다 (`gabes/schemes/rydberg.py:837-923`). phasor까지의 결과만 cache하고
detector/SAM을 다시 매핑하면 slider 응답이 빨라지며 값도 변하지 않는다. 이는 성능 최적화와
P1 cache correctness를 같은 dependency contract로 해결한다.

### 7.4 후순위

- Doppler affine kernel의 `base+sA`를 velocity loop 밖으로 hoist하는 기존 bit-identical 제안은
  현재도 유효하다 (`gabes/kernels.py:426-449`; 이전 진단은 190→166 ms).
- Rydberg CSV duplicate median의 Python group loop(`gabes/rydberg_experimental_csv.py:332-340`)는
  대형 duplicate-heavy trace에서 vectorization 후보지만 일반 trace의 병목은 아니다. 성능보다
  먼저 acquisition-order reversal/hysteresis를 sort 전에 경고해야 한다.

## 8. 문서·테스트·예제 평가

### 잘된 점

- README는 finite-IF, conditional sensitivity, cold spot, SAM, deferred scope를 현재 코드와
  대체로 일치하게 설명한다 (`README.md:94-116`).
- 전용 tests는 reference knobs/linewidth/AT split, finite-IF sideband/low-IF limit,
  peak/RMS, balanced noise, temperature/cold spot, axial integration, SAM uncertainty,
  CSV units/provenance를 폭넓게 검사한다
  (`tests/test_rydberg_eit.py:29-431`; `tests/test_rydberg_electrometry.py:38-203`;
  `tests/test_rydberg_experiment.py:19-173`; `tests/test_rydberg_experimental_csv.py:14-109`).
- `analysis/rydberg_cell_heating/example_config.json`은 measured input을 `PENDING`, ARC dipole을
  `REFERENCE`, detector/geometry를 `ASSUMED`, 결과를 `PREDICTED`로 둔다. raw data 없이
  validation 완료처럼 보이게 하지 않는다.

### 보완할 점

- 사용자 가이드는 여전히 Rydberg를 “정적 스펙트럼, EIT linewidth, AT split, dispersion”으로만
  설명해 finite-IF/noise/cold-spot/SAM을 누락한다
  (`docs/Userguide/GABES_User_Guide_v2.html:534-537`, `603-613`).
- `tests/test_rydberg_eit.py:4-5`도 “public UI shows the static spectrum only”라고 해 현재 코드와
  모순된다.
- 테스트는 내부 수학에는 강하지만 hidden opacity 의미/uncertainty, 논문의 33 °C·6 µW optimum,
  high-temperature hero suppression, fixed reference arm, sub-grid convergence, app-path extra-view
  dependency를 검증하지 않는다. 특히 현재 `test_at_heroes...`는 guardrail 없이 sensitivity를
  hero로 두는 동작을 고정한다 (`tests/test_rydberg_eit.py:99-105`).
- bundled example은 pipeline demonstration이지 measured experiment-model agreement 예제가 아니다.
  최소한 한 개의 provenance-tagged synthetic fixture와, 별도로 실제 held-out trace가 준비될 때의
  expected workflow를 구분해 제공하는 것이 좋다.

## 9. 검증

- Scheme 3 + finite-IF + heated-cell + CSV + workflow 관련:
  `60 passed in 11.90 s`
- 저장소 정책의 전체 회귀:
  `402 passed in 77.41 s`
- production code는 수정하지 않았다. 기존 dirty working tree를 보존하고 이 보고서만 추가했다.

## 10. 결론과 권장 순서

Scheme 3은 toy line shape를 넘어선다. compact cascade OBE, residual Doppler, power/beam/transit,
finite-IF complex response, balanced detector noise, cold spot/axial density, SAM, unit/provenance CSV는
실험 설계와 데이터 정리에 모두 쓸 만한 실제 물리다.

그러나 현재의 `6.08 nV/cm/√Hz`, 47 °C optimum, 4 µW optimum은 **조건부 모델 출력**이지
absolute reference가 아니다. 가장 큰 이유는 hidden opacity, 잘못된 extra-view dependency,
observability guardrail 부재와 held-out 측정 부재다.

권장 순서는 다음과 같다.

1. **P1:** ExtraView별 dependency contract와 app-path 회귀 test를 만들고 atomic/phasor/detector
   cache를 분리한다.
2. **P1:** `ls=0.001`을 provenance와 uncertainty가 있는 effective opacity/participation으로
   공개하고 sweep를 `PREDICTED/UNVALIDATED`로 표시한다.
3. **P1:** contrast, detected power, edge clearance, samples/FWHM, refinement convergence로 hero를
   gate하고 reference arm을 한 detuning에서만 balance한다.
4. **P1/P2:** AT/optimum interpolation과 executable ARC provenance를 추가한다.
5. **속도:** affine Liouvillian assembly와 `rho₀` 재사용을 적용한다.
6. 그 뒤 measured EIT/RF/PSD/temperature holdout을 확보해 compact model을 검증하고, 결과가 요구할
   때만 full Zeeman, measured transfer function, RF field-map reference mode를 활성화한다.

1–5는 solve 차원을 늘리지 않거나 결과를 그대로 보존하면서 신뢰성과 속도를 함께 높인다.
full Zeeman/time-domain/3-D 확장보다 먼저 처리하는 것이 실험 reference로서 가장 비용 효율적이다.
