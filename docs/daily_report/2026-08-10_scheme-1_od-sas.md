# 2026-08-10 Scheme 1 물리 리뷰 — OD / SAS

## 선택 규칙과 현재 다섯 scheme

- 로컬 날짜는 `2026-08-10`, day-of-month는 `10`이다.
- `n = (day mod 5) + 1 = (10 mod 5) + 1 = 1`이므로 오늘의 대상은
  **Scheme 1**이다.
- UI 순서는 `gabes/schemes/__init__.py`의 `_SCHEMES` 리스트가 직접 정한다
  (`gabes/schemes/__init__.py:12-25`, `gabes/schemes/__init__.py:37-39`).

| 번호 | 현재 scheme | registry 인스턴스 | 런타임 이름 / 제목 |
|---:|---|---|---|
| 1 | OD / SAS | `SASScheme()` | `sas` / `Absorption spectroscopy (OD / SAS)` |
| 2 | Lambda coherence (EIT / AT / CPT) | `LambdaScheme()` | `lambda` / `Λ coherence (EIT / AT / CPT)` |
| 3 | Rydberg-EIT electrometry | `RydbergEITScheme()` | `rydberg_eit` / `Rydberg-EIT electrometry` |
| 4 | Hanle / EIA / NMOR | `MagnetoScheme()` | `magneto` / `Magneto-optics (Hanle/MOR)` |
| 5 | FWM | `FWMScheme()` | `fwm` / `Four-wave mixing (Squeezing / Biphoton)` |

Scheme 1의 정의는 `gabes/schemes/sas.py:53-63`에 있다. README의 사용자용
설명도 같은 순서와 범위를 사용한다 (`README.md:8-16`).

## 먼저 검색한 기존 제안과 issue/TODO

물리 판단 전에 다음을 검색·대조했다.

- 저장소 지침과 개요: `CLAUDE.md`, `README.md`
- 현재 계획/TODO: `docs/checklist.json`, `gabes/constants.py:79-89`
- 과거 Scheme 1 보고서:
  `2026-06-25_scheme-1_od-sas.md`,
  `2026-07-10_scheme-1_od-sas.md`,
  `2026-07-20_scheme-1_od-sas.md`,
  `2026-07-30_scheme-1_od-sas.md`,
  `2026-08-05_scheme-1_od-sas.md`
- 구현: `gabes/schemes/sas.py`, `gabes/species.py`, `gabes/constants.py`,
  `gabes/core.py`, `gabes/doppler.py`, `gabes/lineshape.py`,
  `gabes/experimental_csv.py`, `streamlit_app.py`
- 테스트: `tests/test_sas.py`, `tests/test_experimental_csv.py`,
  `tests/test_absorption.py`, `tests/test_headless_observables.py`,
  `tests/test_schemes_render.py`
- 문서·예제: `README.md:83-99`,
  `docs/Userguide/GABES_User_Guide_v2.html:536-565`,
  `docs/Userguide/userguide_assets/od.png`,
  `docs/Userguide/userguide_assets/sas.png`,
  `references/AutoOD/ReferenceOD.csv`

기존 개선 제안은 **있다**. 현재 checklist는 과거 일일 보고서의 핵심을 다음처럼
정리했다.

1. 충돌 계수 provenance, pressure shift, pump/probe homogeneous-width 일관성:
   P1, `ready` (`docs/checklist.json:209-243`).
2. species/line별 Rabi 및 saturation convention:
   P1, `ready` (`docs/checklist.json:246-264`).
3. sub-Doppler linewidth/lock-slope 신뢰성:
   `done` (`docs/checklist.json:267-285`).
4. CSV 절대 보정과 sweep 방향 보존:
   `done` (`docs/checklist.json:288-305`).
5. full velocity-changing-collision(VCC) kernel:
   P2 연구 항목, `parked` (`docs/checklist.json:563-582`).

과거의 자유로운 low-order polarization/Zeeman proxy 제안은 독립 측정 없이
fit knob만 늘릴 위험 때문에 현재 checklist에서 `rejected`되었다
(`docs/checklist.json:38-42`). 검증 데이터 없는 phenomenological Dicke knob도
충돌 계수 항목의 1차 범위에서 명시적으로 제외되었다
(`docs/checklist.json:238-241`). 이 두 결정은 실험 reference의 식별 가능성을
지키는 방향으로 타당하다.

별도 로컬 issue/proposal/note 파일은 `docs/checklist.json` 외에는 발견되지 않았다.
공개 [GitHub Issues](https://github.com/Shake2313/fwm-squeezing-app/issues)는
2026-08-10 확인 시 open/closed 모두 0건이었다. 코드 TODO는 Ne broadening을
gas/species/line table과 pressure shift로 승격하라는 한 항목이 직접 관련된다
(`gabes/constants.py:79-89`).

## 결론 요약

**현재 OD/SAS는 실제 원자물리를 구현한 유용한 warm-vapor 실험 계획용
semi-quantitative reference다.** 실제 Rb/Cs hyperfine 상수, 질량, 자연 선폭,
Wigner-6j/3j 선세기, Maxwell 속도 평균, counter-propagating pump/probe,
CG-branched spontaneous decay, transit relaxation, hyperfine optical pumping을
사용한다 (`gabes/schemes/sas.py:1-37`, `gabes/species.py:127-182`,
`gabes/species.py:336-434`). 특히 pump-off 85Rb D1 절대 OD와 pump-on
enhanced/inverted crossover가 한 모델에서 연속적으로 나온다는 점은 교육용 toy를
넘어서는 강점이다.

이번 current working tree에서는 지난 리뷰의 큰 신뢰성 결함 두 가지가 실제로
고쳐졌다.

- under-resolved SAS에서는 numeric linewidth/lock slope를 hero로 올리지 않고
  `resolution-limited` 상태를 올린다 (`gabes/schemes/sas.py:371-390`,
  `gabes/lineshape.py:141-250`).
- CSV extrema scaling은 명시적으로 상대 정규화이며, 절대 투과는 dark/reference
  또는 gain/offset 증거가 있을 때만 허용된다. 원래 취득 순서의 forward/reverse
  branch도 보존하고, 병합 overlay의 hysteresis 손실을 경고한다
  (`gabes/experimental_csv.py:212-321`, `gabes/experimental_csv.py:441-520`,
  `streamlit_app.py:931-977`, `streamlit_app.py:1115-1140`).

그러나 다음 용도에는 아직 **정량 reference로 부적합**하다.

- buffer-cell SAS의 절대 선중심, 폭, saturation/hole-burning
- 고온 Cs 및 natural-Rb mixture의 self-broadened linewidth
- species/line/polarization 사이의 절대 saturation-power 비교
- mF·편광·Zeeman-resolved contrast/branching
- 변조·복조 전달함수를 포함한 실제 servo error signal
- 외부 측정으로 검증된 pumped-SAS 절대 contrast와 linewidth

따라서 pure Rb 셀의 line finding, pump-off OD sanity check, pump power·온도·transit
경향, 해상도 상태가 `resolved`인 경우의 lock 후보 탐색에는 적합하다. 반대로
정밀 linewidth·collision coefficient·absolute discriminator calibration에는 현재
숫자를 그대로 인용하면 안 된다.

## 현재 구현이 담는 실제 물리

### 1. pump-off OD와 pump-on SAS가 같은 모델이다

laser scan을 `Δ`, 원자 속도를 `v`라 하면 pump와 probe의 atom-frame detuning은
각각 `Δ + kv`, `Δ - kv`로 들어간다 (`gabes/schemes/sas.py:4-30`). Pump가 0이면
population factor가 정확히 1로 줄어 unit-area Lorentzian의 Maxwell 평균, 즉
multi-line Voigt OD가 된다. Pump를 켜면 velocity-selective saturation과
hyperfine optical pumping이 같은 스펙트럼 위에 Lamb dip과 crossover를 만든다
(`gabes/schemes/sas.py:211-253`).

실험 노브도 온도, cell length, isotope/natural abundance, D1/D2, pump power,
waist, transit relaxation, Ne pressure로 구성되어 있다
(`gabes/schemes/sas.py:65-107`). `cell_mm`와 `line_strength`는 post-process knob라
불필요한 OBE 재계산을 피한다 (`gabes/schemes/sas.py:82-99`,
`gabes/schemes/sas.py:324-331`).

### 2. hyperfine optical pumping은 실제로 구현되어 있다

Rb-85, Rb-87, Cs-133의 hyperfine A/B, centroid, mass, abundance와 자연 선폭이
명시돼 있다 (`gabes/species.py:127-182`). 허용 `Fg→Fe` 선세기는
Wigner-6j/3j에서 계산하고, 같은 angular factor가 absorption weight, relative
pump Rabi, spontaneous branching에 일관되게 들어간다
(`gabes/species.py:226-270`, `gabes/species.py:364-419`).

`build_manifold()`는 모든 허용 F-state와 CG-branched decay를 만들고, transit
relaxation으로 population을 thermal ground distribution으로 돌려보낸다
(`gabes/species.py:336-396`). 낮은 transit rate에서 crossover transmission이
커지는 테스트는 단순 Lorentzian 합이 아니라 optical pumping이 작동한다는
유용한 내부 검증이다 (`tests/test_sas.py:144-158`).

다만 각 F는 하나의 lumped state다. mF, pump/probe polarization, beam profile,
magnetic field, diffusion은 포함되지 않는다. 따라서 line centre와 qualitative
crossover 구조는 유용하지만 absolute pumped contrast와 branching은 실험별
reference가 아니다.

### 3. pump-off AutoOD anchor는 좋지만 regression 계약은 느슨하다

85Rb D1, 90 °C, 12.5 mm, pump off 조건에서 internal AutoOD primitive와 다시
비교한 결과는 다음과 같다.

| 비교량 | 현재 차이 |
|---|---:|
| integrated absorption area | `0.039317 %` |
| peak absorption | `0.012573 %` |

즉 현재 값 자체는 공개 문서의 `<0.1 %` 약속을 만족한다
(`gabes/schemes/sas.py:123-151`, `README.md:259-268`). 그러나 직접 보호하는
테스트는 여전히 area/peak 모두 `<1 %`만 요구한다
(`tests/test_sas.py:101-111`). 향후 0.5 % drift도 테스트를 통과할 수 있으므로,
같은 solve를 사용하는 assertion을 `1e-3` 이하로 강화해야 한다. 계산 overhead는
없다.

## 현재 정량 상태

### 1. sub-Doppler readout: 과신은 해결, 기본 grid는 여전히 제한적

현재 natural-Rb D1 기본값은 `P=0.5 mW`, `40 °C`, `75 mm`, 1401 points다.
headless readout은 다음을 낸다.

- hero: `SAS resolution = resolution-limited`
- provisional sub-Doppler FWHM: `21.56 MHz`
- samples/FWHM: `2.7`
- Peak OD: `0.54`
- calculated single-line Gaussian Doppler FWHM: `518.7 MHz`

따라서 기본 linewidth를 숨기지는 않지만 신뢰 가능한 hero로 승격하지 않는 현재
동작이 옳다. UI가 허용하는 4001 points에서는 같은 natural-Rb 기본값이
`resolved`, FWHM `20.36 MHz`, samples/FWHM `7.2`, lock slope
`0.0316 /MHz`가 된다.

checklist acceptance와 같은 predeclared 85Rb D1 fixed-window refinement를 다시
측정한 결과는 다음과 같다.

| points | warm compute | 상태 | interpolated FWHM | samples/FWHM | local slope |
|---:|---:|---|---:|---:|---:|
| 1401 | `0.2334 s` | resolution-limited | `18.7049 MHz` | `3.73` | provisional `0.0337 /MHz` |
| 2801 | `0.4416 s` | resolved | `17.8414 MHz` | `7.11` | `0.0379 /MHz` |
| 5601 | `1.0458 s` | resolved | `17.7409 MHz` | `14.14` | `0.0394 /MHz` |

2801→5601에서 FWHM 차이는 약 `0.56 %`, slope 차이는 약 `3.81 %`다.
이 결과는 checklist의 완료 주장과 일치한다
(`docs/checklist.json:267-284`, `tests/test_sas.py:244-294`). 5601 points는
UI 범위를 넘는 reference-only check이며, 매번 강제할 필요가 없다. 현행 O(N)
상태 계산은 추가 OBE solve 없이 물리를 보존한다.

### 2. collision model은 여전히 P1 경계다

`gamma_eff = gamma_nat + self + Ne broadening`은 계산하지만
(`gabes/schemes/sas.py:169-185`), pump Liouvillian은 natural-γ manifold로 먼저
만들어진다 (`gabes/species.py:347-396`, `gabes/schemes/sas.py:221-227`).
`gamma_eff`는 pump table 범위/간격과 probe Lorentzian에만 직접 들어간다
(`gabes/schemes/sas.py:225-252`).

예를 들어 85Rb D1, Ne 20 Torr에서는 현재 probe FWHM convention이
`5.75 + 3.91×20 = 83.95 MHz`인데 pump spontaneous/dephasing dynamics는 여전히
5.75 MHz natural manifold다. 따라서 pressure-broadened hole burning과
saturation은 내부적으로 일관되지 않는다. Pressure shift와 Dicke narrowing도
없다 (`gabes/constants.py:79-89`, `gabes/schemes/sas.py:84-87`).

현재 P1 checklist가 요구하듯 pressure shift와 elastic optical dephasing을
분리한 provenance-backed `(gas, species, line)` table을 만들고, pump/probe에
동일한 homogeneous-width convention을 적용해야 한다
(`docs/checklist.json:209-236`). 데이터/API/테스트 작업은 cross-scheme라
`effort=large`지만, 런타임은 table lookup과 같은 크기의 Liouvillian이므로
거의 늘지 않는다. Pure-cell 85Rb D1 anchor는 pressure=0 회귀로 보존할 수 있다.

### 3. Rb self-broadening 계수가 Cs에도 적용된다

`BETA_SELF` 주석은 Rb 계수라고 명시하지만 (`gabes/species.py:49-50`),
`self_broadened_gamma(iso, N)`는 `iso`를 사용하지 않고 `BETA_SELF*N`만 반환한다
(`gabes/species.py:211-213`). Scheme 1은 이 값을 Cs에도 적용한다
(`gabes/schemes/sas.py:175-185`).

| Cs 온도 | 현재 코드가 더하는 `β_Rb N_Cs / 2π` |
|---:|---:|
| 22 °C | `0.002186 MHz` |
| 30 °C | `0.004776 MHz` |
| 100 °C | `1.018599 MHz` |
| 150 °C | `15.175086 MHz` |
| 200 °C | `124.009810 MHz` |

기본 Cs preset에서는 작지만 허용 slider 상단에서는 자연 선폭보다 커진다. 또
natural Rb에서 `number_density()`는 elemental total density를 반환한 뒤 abundance를
곱하므로 (`gabes/species.py:193-207`, `gabes/schemes/sas.py:175-182`), self collision이
total-Rb density를 써야 하는지 isotope-partial density를 써야 하는지도 coefficient
정의와 함께 고정해야 한다. 종/선별 계수와 `unsupported/qualitative` 상태를 넣는
런타임 비용은 O(1)이다.

### 4. 모든 species/line이 하나의 Rb 기준 saturation intensity를 쓴다

`I_SAT = 4.484 mW/cm²` 한 값이 constants에 있고 (`gabes/constants.py:25-31`),
`pump_rabi_from_power()`가 모든 Rb/Cs D1/D2에서 이를 재사용한다
(`gabes/species.py:273-283`). 현재 기본 `0.5 mW`, `w=1 mm`에서는 모든 종/선이
동일하게 `Ω/γ ≈ 1.884`가 된다. Pump-power trend 자체는 유용하지만 종과 D-line
사이의 절대 saturation power 비교에는 근거가 부족하다.

P1 `species-line-rabi-and-saturation-provenance` 항목처럼 dipole, linewidth,
polarization, intensity, angular/cyclic-frequency convention을 기록한 lookup/derived
scale을 사용해야 한다 (`docs/checklist.json:246-263`). 계산비용은 O(1)이고
pump-off AutoOD anchor에는 영향을 주지 않는다. 독립 검증 없이 자유 보정 knob를
늘리는 방식은 피해야 한다.

## CSV 비교 경로의 현재 가치와 한계

CSV 개선은 과거 제안대로 구현됐다.

- extrema mode는 `Relative normalized transmission`으로 명시된다
  (`gabes/experimental_csv.py:24-29`, `gabes/experimental_csv.py:901-992`).
- dark/reference와 signed gain/offset 모드는 절대 투과를 복원하고, 0–1 밖의
  calibrated sample을 숨기지 않는다 (`gabes/experimental_csv.py:754-891`).
- sweep reversal은 sort/merge 전에 센다. Acquisition-order `SweepBranch`가 양방향
  trace를 보존한다 (`gabes/experimental_csv.py:267-312`,
  `gabes/experimental_csv.py:441-495`).
- UI는 merged compatibility overlay가 hysteresis를 잃는다고 명시적으로 경고한다
  (`streamlit_app.py:1115-1122`).

합성 절대 보정은 `rtol=1e-9`로 원래 transmission을 복원하고
(`tests/test_experimental_csv.py:117-159`), round-trip sweep의 5 forward/5 reverse
sample과 두 branch를 정확히 보존한다 (`tests/test_experimental_csv.py:223-247`).
이 개선은 O(N) import/post-process라 OBE overhead가 0이다.

실험 관점의 남은 작은 gap은 두 가지다.

1. 기본 relative mode는 여전히 absolute OD/contrast가 아니다. 사용자는 반드시
   measured dark/reference 또는 detector law를 제공해야 한다.
2. branch는 데이터 구조에 보존되지만 현재 overlay는 merged trace만 그린다
   (`streamlit_app.py:1149-1180`). Hysteresis를 실제 비교해야 한다면 forward/reverse
   branch selector/overlay를 추가할 수 있다. 이는 O(N) plotting이며 OBE solve는
   늘지 않는 P3 usability 개선이다.

README는 이 계약을 정확히 설명한다 (`README.md:83-99`). 반면 User Guide의 OD/SAS
절은 새 CSV calibration modes를 설명하지 않는다
(`docs/Userguide/GABES_User_Guide_v2.html:536-565`). 실험 reference로 배포하려면
README의 calibration/hysteresis 문단을 가이드에도 옮기는 것이 좋다.

## 문서와 정적 예제 평가

`od.png`와 `sas.png`는 같은 natural-Rb D1 조건에서 pump off/on 차이를 즉시 보여
교육적 가치가 높다. 하지만 두 PNG는 Transmission과 OD를 한 이미지에 쌓은 이전
layout이다. 현재 코드는 두 figure view를 분리한다
(`gabes/schemes/sas.py:342-369`). 또한 SAS PNG에는 현재의
`resolution-limited` hero, samples/FWHM, half-height edge가 보이지 않는다.

가이드 본문은 새 resolution contract를 정확히 설명하므로
(`docs/Userguide/GABES_User_Guide_v2.html:559-565`), PNG를 현재 UI에서 다시 만들고
조건·commit/cache version·resolution status를 caption에 기록하면 코드와 예제가
실험 reference로서 다시 일치한다. 계산비용이 아니라 문서 재생성 작업이다.

## 기존 개선안의 계산비용과 물리 보존성

| 기존 제안 | 현재 상태 | 계산비용 | 물리 보존성 / 판단 |
|---|---|---:|---|
| sub-Doppler interpolation/status/hero gating | 완료 | O(N) post-process | OBE 불변, 과신 방지. 유지 |
| CSV 상대/절대 calibration + sweep branch | 완료 | O(N) import | OBE 불변, 실험 의미 복구. 유지 |
| species/line/gas 충돌 계수 + pressure shift | P1 ready | O(1) lookup, 같은 solve 차원 | sourced coefficient와 zero-pressure 회귀가 있으면 pure-cell anchor 보존 가능 |
| added width의 pump optical dephasing 반영 | P1 ready의 일부 | 같은 Liouvillian 크기·solve 수 | 현재 buffer-cell inconsistency를 고치는 저-overhead physics change |
| species/line별 Rabi/saturation provenance | P1 ready | O(1) | pump-off anchor 불변; pumped 결과는 의도적으로 교정 |
| AutoOD `<0.1 %` 테스트 강화 | 과거 제안, 미완료 | 같은 solve | 문서 약속만 더 강하게 보호 |
| phenomenological Dicke knob | 1차 범위 제외 | scalar면 낮음 | 검증 데이터 없이는 bias/fit knob 위험. 현재 제외가 타당 |
| low-order polarization/Zeeman proxy | rejected | 낮을 수 있음 | 독립 측정 없이는 식별 불가. 재도입 비권고 |
| full VCC kernel | parked | velocity class 결합으로 높음 | 특정 cell/kernel/dataset/runtime이 정해진 opt-in reference mode만 타당 |
| full mF/polarization solve | 현행 범위 밖 | density-matrix 차원 급증 | absolute contrast가 실제 요구되고 검증 trace가 있을 때만 별도 설계 |

핵심은 **구현 노력과 계산 overhead를 구분하는 것**이다. 충돌 coefficient API는
cross-scheme provenance 작업이라 개발 effort는 크지만 interactive runtime은 거의
늘지 않는다. 반대로 full VCC는 물리적으로 좋을 수 있어도 velocity separability를
깨므로 negligible-overhead 제안이 아니다.

## 동작을 바꾸지 않는 순수 코드 최적화

### 1. pump-off analytic population fast path — 가장 큰 solver-side 이득

Pump가 0이면 각 transition의 population factor는 정확히 `w=1`인데 현재도
Liouvillian, `_pump_pops()` table, level별 interpolation array를 만든다
(`gabes/schemes/sas.py:221-252`, `gabes/schemes/sas.py:505-512`).

85Rb D1, 90 °C, 12.5 mm, 1401 points에서 같은 scan/probe Voigt를 직접 계산한
benchmark는 다음과 같다.

| 경로 | median runtime |
|---|---:|
| 현재 | `0.246875 s` |
| analytic `w=1` | `0.081636 s` |
| speedup | **`3.024×`** |

scan axis는 exact equal이고 max relative alpha 차이는 `7.95×10⁻16`이었다. 조건을
`pump_power_mw <= 0`으로 제한하고 AutoOD·natural-Rb·buffer/temperature 회귀로
보호하면 물리와 public behavior를 바꾸지 않는 최우선 최적화다.

### 2. running median vectorization — 현재 구현 완료

현재 `subdoppler_feature()`는 `sliding_window_view`와 vector median을 사용한다
(`gabes/lineshape.py:166-173`). 같은 1401-point window에서 이전 Python loop와
재비교한 결과는 `20.9345 ms → 0.12455 ms`, **168.1×**, max difference 0이었다.
과거 제안이 정확히 구현된 positive result이며 되돌릴 이유가 없다.

### 3. duplicate-heavy CSV median merge vectorization — 아직 유효

`_sort_and_merge()`는 duplicate group마다 Python에서 `np.median`을 호출한다
(`gabes/experimental_csv.py:506-520`). CSV 상한과 같은 500,000 rows,
각 x 5회 반복 synthetic trace에서 다음을 측정했다.

| 구현 | median runtime |
|---|---:|
| 현재 group loop | `1.666681 s` |
| `lexsort((y,x))` + vector median indices | `0.146120 s` |
| speedup | **`11.41×`** |

unique x는 exact equal, merged y max difference는 0이었다. Odd/even group,
signed zero, extreme finite values, stable ordering을 추가 회귀로 보호하면 public
behavior를 바꾸지 않는 좋은 importer 최적화다. 일반 실험 trace에서는 OBE가 더
큰 비용이지만 row-limit 근처 파일에는 직접 체감된다.

### 4. 낮은 우선순위

Pump-on 경로는 같은 `deff_grid`에 대해 level마다 `np.interp` search를 반복한다
(`gabes/schemes/sas.py:237-252`). 공통 index/weight geometry를 재사용할 수 있지만
strict equivalence와 peak memory benchmark가 먼저다. `build_manifold()`는 이미
`lru_cache(maxsize=64)`이므로 (`gabes/species.py:336-337`) angular-data cache 확대는
우선순위가 낮다.

## 테스트 결과와 해석

Scheme 1 관련 테스트:

```text
python -m pytest tests/test_sas.py tests/test_experimental_csv.py \
  tests/test_absorption.py tests/test_headless_observables.py \
  tests/test_schemes_render.py -q

93 passed in 15.80s
```

전체 저장소 테스트:

```text
MPLBACKEND=Agg python -m pytest -q

402 passed in 62.01s
```

테스트는 current working tree의 수치/렌더/API 계약이 일관됨을 보여 준다. 특히 새
resolution/CSV tests는 과거 trust 문제를 직접 보호한다
(`tests/test_sas.py:208-294`, `tests/test_experimental_csv.py:100-247`). 그러나
테스트 통과가 omitted mF/polarization/collision physics나 held-out 실험 검증을
대체하지는 않는다. Pumped-SAS 검증은 현재 sharp feature, total absorption 감소,
transit-rate trend 같은 내부·정성 검사 중심이다 (`tests/test_sas.py:133-158`).

## 최종 우선순위

1. **P1, runtime 거의 0:** `(gas, species, line)` 충돌 coefficient provenance와
   unsupported status를 만들고, pressure shift와 optical dephasing을 분리해
   pump/probe width convention을 일치시킨다.
2. **P1, O(1):** species/line/polarization별 Rabi/saturation scale을 근거와 함께
   도입한다. Pump-off AutoOD는 frozen anchor로 유지한다.
3. **P2, 추가 solve 0:** AutoOD 공개 `<0.1 %` 약속에 맞춰 regression tolerance를
   `1e-3` 이하로 강화한다.
4. **P2/P3, 계산비용 없음:** current UI로 OD/SAS PNG를 재생성하고 User Guide에
   CSV calibration/hysteresis contract를 넣는다.
5. **순수 성능:** pump-off analytic path를 먼저, duplicate-heavy CSV vector merge를
   다음으로 구현한다. Running median 최적화는 이미 완료됐다.
6. **필요할 때만 heavy:** full VCC와 full mF/polarization solver는 특정 실험 cell,
   sourced kernel, measured trace, held-out validation, runtime budget이 정해진 뒤
   opt-in reference mode로 설계한다.

종합하면 현재 Scheme 1의 가장 좋은 정체성은 **AutoOD-anchored pump-off OD와
hyperfine-pumping SAS를 연결한 빠른 실험 계획 도구**다. 새 resolution과 CSV
계약 덕분에 연구자가 숫자를 잘못 신뢰할 위험은 크게 줄었다. 다음 단계는 solver를
무겁게 만드는 것이 아니라, 충돌·Rabi provenance와 공개 검증 계약을 먼저 닫아
어디까지 정량 reference인지 더 명확하게 만드는 것이다.
