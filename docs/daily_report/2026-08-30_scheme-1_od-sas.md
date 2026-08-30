# 2026-08-30 Scheme 1 물리 리뷰 — OD / SAS

## 선택 규칙과 현재 다섯 scheme

서울 현지 날짜의 일(day)은 30이므로

```text
n = (30 mod 5) + 1 = 1
```

오늘의 대상은 **Scheme 1, Absorption spectroscopy (OD / SAS)** 다. 현재 UI
등록 순서는 `gabes/schemes/__init__.py:19-25`에 다음과 같이 고정돼 있다.

| 번호 | 등록 인스턴스 | 사용자-facing scheme |
|---:|---|---|
| 1 | `SASScheme()` | OD / SAS |
| 2 | `LambdaScheme()` | Lambda coherence (EIT / AT / CPT) |
| 3 | `RydbergEITScheme()` | Rydberg-EIT electrometry |
| 4 | `MagnetoScheme()` | Hanle / EIA / NMOR |
| 5 | `FWMScheme()` | seeded FWM gain diagnostic / generic SFWM biphoton |

`ODScheme`은 현재 dropdown scheme이 아니라 2-level 내부 검증 primitive로만 남아
있다 (`gabes/schemes/__init__.py:12-18`).

## 먼저 검색한 기존 제안·TODO·issue note

다음 자료를 먼저 확인했다.

- 저장소 지침과 현재 설명: `CLAUDE.md`, `README.md`
- 작업 registry: `docs/checklist.json`
- Scheme 1 과거 보고서 일곱 건:
  `2026-06-25`, `2026-07-10`, `2026-07-20`, `2026-07-30`, `2026-08-05`,
  `2026-08-10`, `2026-08-25`
- 구현: `gabes/schemes/sas.py`, `gabes/species.py`, `gabes/constants.py`,
  `gabes/core.py`, `gabes/doppler.py`, `gabes/lineshape.py`,
  `gabes/experimental_csv.py`
- 테스트: `tests/test_sas.py`, `tests/test_experimental_csv.py`,
  `tests/test_absorption.py`, `tests/test_headless_observables.py`,
  `tests/test_schemes_render.py`
- 문서·정적 예제: `README.md:8-16`, `README.md:87-103`,
  `README.md:280-317`, `docs/Userguide/GABES_User_Guide_v2.html:520-584`,
  `docs/Userguide/userguide_assets/od.png`,
  `docs/Userguide/userguide_assets/sas.png`, `references/AutoOD/ReferenceOD.csv`

기존 개선 제안은 **충분히 존재한다**. 현재 checklist의 Scheme 1 관련 항목은
다음과 같다.

1. 충돌 계수 provenance, pressure shift, pump/probe homogeneous-width 일관성:
   P1, `ready`, effort `large` (`docs/checklist.json:209-243`).
2. species/line별 Rabi와 saturation convention:
   P1, `ready`, effort `medium` (`docs/checklist.json:246-264`).
3. resolution-aware sub-Doppler readout:
   P1, `done` (`docs/checklist.json:267-285`).
4. OD/SAS paraffin population memory MVP:
   P2, `done` (`docs/checklist.json:288-301`).
5. cell/beam geometry와 dark-flight를 포함한 coated-cell 확장:
   P2, research, `parked` (`docs/checklist.json:304-325`).
6. CSV 절대 보정과 sweep 방향 보존:
   P1, `done` (`docs/checklist.json:328-350`).
7. full velocity-changing-collision(VCC) kernel:
   P2, research, `parked` (`docs/checklist.json:603-625`).

코드의 직접 TODO는 `gabes/constants.py:79-89`의 Ne broadening scalar를
gas/species/line coefficient table과 pressure shift로 승격하라는 항목 하나다. 별도
로컬 issue/proposal/note 파일이나 `.github` issue 문서는 발견되지 않았다. 연결된
`Shake2313/fwm-squeezing-app` GitHub issue도 오늘 확인 시 open 0건, closed 0건이었다.

8월 25일 Scheme 1 리뷰 뒤에는 실제 변경이 있다. 8월 28일의 `5cc7eee`가 opt-in
paraffin population-memory 경로를 추가했고, 뒤의 validation 정리도 이 코드를 현재
형태로 문서화했다. 따라서 이번 보고서는 이전 결론을 단순 반복하지 않고 이 새 경로를
중점적으로 평가한다.

## 결론 요약

**비코팅 OD/SAS 기본 경로는 실제 alkali hyperfine data와 광펌핑을 담은 유용한
warm-vapor 실험 계획용 정성/반정량 reference다.** Pump-off 85Rb D1 절대 OD,
isotope/hyperfine line finding, pump와 transit-rate trend, 초기 lock 후보 탐색에는
단순 Lorentzian 합보다 분명히 가치가 크다.

새 **paraffin 옵션도 물리적으로 무의미한 transmission multiplier가 아니다.** 각
입사 ground-F 상태에 대한 조건부 OBE, 자연붕괴 뒤 ground-F transfer, 벽 완화,
velocity rethermalization을 연결하므로 “장시간 hyperfine population memory가 SAS를
어떻게 바꿀 수 있는가”를 보는 정성적 가설 도구로 유용하다
(`gabes/schemes/sas.py:238-306`, `gabes/schemes/sas.py:582-701`).

그러나 이를 **실제 paraffin cell의 절대 스펙트럼 reference로 쓰면 안 된다.** 현재
`T1=25.1 ms`는 한 87Rb cell의 300 K 측정값이고, 동일 `transit_khz`가 bright-region
exit와 dark return cadence를 모두 맡는다 (`gabes/schemes/sas.py:147-153`,
`gabes/schemes/sas.py:671-688`). Cell/beam volume ratio, dark-flight distribution,
scan rate, Zeeman coherence, coating/adsorption shift, wall-collision statistics,
species·temperature·ageing별 T1이 없다. 오늘 수치에서도 이 가정이 거의 완전한
dark-hyperfine pumping을 만들어 결과가 parameterization에 매우 민감했다.

따라서 전체 Scheme 1을 가장 정확히 부르면 다음과 같다.

> **AutoOD-anchored pump-off OD와 hyperfine-pumping SAS를 연결한 빠른 실험 계획
> 모델. Paraffin mode는 calibrated coated-cell solver가 아니라 quasi-static
> population-memory what-if reference.**

다음 용도에는 여전히 정량 reference로 부적합하다.

- buffer-cell SAS의 절대 line centre, linewidth, Dicke/VCC narrowing
- 고온 Cs 또는 natural-Rb mixture의 정확한 self broadening
- species/D-line/polarization 사이의 절대 saturation-power 비교
- mF·편광·Zeeman-resolved contrast와 branching
- 변조·복조 전달함수를 포함한 servo error signal
- coated-cell scan hysteresis, Ramsey bright/dark evolution, wall shift
- held-out 측정으로 검증된 pumped-SAS absolute contrast/linewidth

## 현재 구현이 담는 실제 물리

### 1. Pump-off OD와 pump-on SAS가 같은 모델에서 이어진다

레이저 detuning을 `Δ`, 원자 속도를 `v`라 하면 pump/probe는 각각 `Δ+kv`와
`Δ-kv`를 본다 (`gabes/schemes/sas.py:4-30`). Pump가 꺼지면 population factor가
정확히 1이 되어 unit-area Lorentzian의 Maxwell 평균, 즉 multi-line Voigt OD가
된다. Pump가 켜지면 velocity-selective saturation과 hyperfine optical pumping이
같은 흡수 스펙트럼 위에 Lamb dip과 crossover를 만든다
(`gabes/schemes/sas.py:248-305`).

Rb-85, Rb-87, Cs-133의 hyperfine A/B, centroid, mass, abundance, 자연 선폭은
`gabes/species.py:127-182`에 있고, Wigner-6j/3j 기반 세기는 absorption weight,
relative pump Rabi, spontaneous branching에 일관되게 들어간다
(`gabes/species.py:226-270`, `gabes/species.py:366-422`). Transit reset은 thermal
ground-F distribution으로 population을 되돌린다 (`gabes/species.py:337-399`).

이는 실험적으로 유용한 실제 물리다. 다만 각 F는 하나의 lumped state라 mF,
polarization, magnetic field, beam profile, spatial diffusion은 표현하지 않는다.

### 2. 새 paraffin 경로는 population transfer를 실제로 푼다

코팅 옵션은 Advanced, 기본 off, heavy-cache key에 포함되고 Generic toy에서는
숨겨진다 (`gabes/schemes/sas.py:85-92`, `tests/test_sas.py:109-115`). Pump가 0이면
legacy OD 경로를 그대로 사용해 coated/uncoated array가 bitwise identical이다
(`gabes/schemes/sas.py:187-190`, `tests/test_sas.py:228-235`).

Pump-on coated path는 다음 순서다.

1. 각 incoming ground-F basis state로 reload되는 조건부 bright-region steady state를
   푼다 (`gabes/schemes/sas.py:582-618`).
2. Beam을 떠날 때의 excited population을 natural spontaneous branching으로 ground
   reservoir에 투영한다 (`gabes/schemes/sas.py:621-668`).
3. `cycle_rate(I-T)`와 `1/T1` wall relaxation을 합친 2×2 정상상태를 푼다
   (`gabes/schemes/sas.py:671-693`).
4. Pump가 준비한 population과 probe line profile을 **같은 velocity class**에서
   곱한 뒤 Maxwell 평균한다 (`gabes/schemes/sas.py:273-305`,
   `tests/test_sas.py:204-211`).

확률 보존, basis-state 선형성, far-wing axis coverage, thermal identity limit는
테스트된다 (`tests/test_sas.py:118-201`). 이는 구현 수준의 좋은 검증이다. 그러나
현재 end-to-end physics test는 coated 결과가 finite/nonnegative이며 uncoated와
`0.1%` 이상 다르다는 것만 요구한다 (`tests/test_sas.py:238-252`). 실제 coated-cell
linewidth/contrast나 scan dynamics를 고정하는 독립 reference는 아니다.

### 3. CSV 비교 경로는 의미 경계를 잘 지킨다

CSV extrema scaling은 `Relative normalized transmission`으로 표시되고, measured
dark/reference 또는 detector gain/offset이 있을 때만 absolute transmission을 만든다
(`README.md:87-103`, `gabes/experimental_csv.py:754-891`). Acquisition-order
forward/reverse branch도 정렬 전에 보존한다. 이는 “모양 맞춤”과 절대 투과 검증을
구분하는 실험 reference로서 중요한 장점이다.

## 오늘 재확인한 정량 상태

### 1. 비코팅 기본값은 여전히 resolution-limited다

현재 natural-Rb D1 기본값(`P=0.5 mW`, `40 °C`, `75 mm`, 1401 points)의 warm
median compute time은 오늘 환경에서 `0.715 s`였다.

| readout | 현재 값 |
|---|---:|
| status | `resolution-limited` |
| provisional sub-Doppler FWHM | `21.56 MHz` |
| samples/FWHM | `2.7` |
| peak OD | `0.54` |
| calculated single-line Gaussian Doppler FWHM | `518.7 MHz` |

현재 코드는 이 상태에서 numeric lock slope를 hero로 내지 않는 것이 옳다
(`gabes/schemes/sas.py:424-443`, `tests/test_sas.py:357-384`).

### 2. Pump-off AutoOD anchor는 좋지만 공개 계약보다 테스트가 느슨하다

85Rb D1, 90 °C, 12.5 mm, pump off에서 내부 `ODScheme` primitive와 비교한 오늘
결과는 다음과 같다.

| 비교량 | 상대 차이 |
|---|---:|
| integrated absorption area | `-0.0393169 %` |
| peak absorption | `-0.0125734 %` |

따라서 README와 in-app 설명의 `<0.1 %` 약속은 현재 만족한다
(`gabes/schemes/sas.py:132-146`, `README.md:289-297`). 그러나 직접 회귀는 아직
`<1 %`만 요구한다 (`tests/test_sas.py:215-225`). 같은 계산의 assertion만
`1e-3` 이하로 강화하면 추가 solve 없이 공개 계약을 실제로 보호할 수 있다.

### 3. Paraffin 결과는 강하고, 가정에 매우 민감하다

새 경로를 보기 위해 pure 85Rb D1, 45 °C, 75 mm, 0.5 mW, 401 points를 비교했다.
시간은 warm median 3회다.

| 항목 | uncoated | coated |
|---|---:|---:|
| compute | `0.0981 s` | `0.3102 s` |
| peak `alpha` | `35.9987 m⁻¹` | `0.0659012 m⁻¹` |
| minimum transmission | `0.06721` | `0.99507` |
| integrated absorption ratio | 1 | `0.0050601` |

즉 이 fixture에서 코팅 memory는 적분 흡수를 비코팅의 약 **0.506%**로 줄인다.
Reservoir의 thermal ground-F population은 `[0.4167, 0.5833]`이지만 scan 전체에서
각 component는 약 `0.000737`부터 `0.999188`까지 움직였다. 기본 wall/cycle rate
ratio는 `6.34×10⁻5`라 수천~수만 번의 effective return이 벽 완화보다 먼저
누적된다.

같은 transfer map에서 T1만 바꾼 진단은 대표 두 hyperfine line에서 ground-0
population을 다음처럼 바꿨다. 이는 production knob sweep가 아니라 현재 고정값의
민감도를 확인한 것이다.

| T1 | line A | line B |
|---:|---:|---:|
| `1 µs` | `0.3822` | `0.4720` |
| `100 µs` | `0.0423` | `0.9486` |
| `25.1 ms` | `0.00111` | `0.99907` |

강한 dark-state pumping 자체는 충분히 가능한 물리 방향이다. 하지만 이 정도의
결과는 literature T1 하나와 return cadence 가정에 의해 지배되므로, cell-specific
측정 없이 절대 contrast를 인용하면 안 된다.

수치 해상도도 주의해야 한다. 같은 coated fixture를 401/801/1401 points로 올렸을
때 적분 면적은 빠르게 수렴했지만 401-point curve를 1401 grid에 보간한 최대 차이는
1401 peak의 `17.7%`였다. 1401 points에서도 narrow feature는 `4.69 samples/FWHM`로
`resolution-limited`다. 현재 status gate는 유효하지만, coated spectrum을 정량
비교하려면 별도 fixed-window refinement가 필요하다.

### 4. 기존 collision/Rabi trust 경계는 변하지 않았다

`gamma_eff = gamma_nat + self + Ne broadening`은 probe Lorentzian과 pump-state
detuning grid에는 들어가지만, species pump manifold의 spontaneous/dephasing은
natural Γ로 지어진다 (`gabes/schemes/sas.py:196-209`,
`gabes/schemes/sas.py:248-268`, `gabes/species.py:337-399`). 예를 들어 85Rb D1,
Ne 20 Torr의 probe FWHM convention은 `5.75 + 3.91×20 = 83.95 MHz`인 반면 pump
dynamics는 natural 5.75 MHz다. Pressure shift, Dicke narrowing, VCC도 없다
(`gabes/constants.py:79-89`).

`self_broadened_gamma(iso, N)`는 `iso`를 사용하지 않고 Rb `BETA_SELF`를 모든
species에 적용한다 (`gabes/species.py:211-213`). 이 때문에 Cs에 잘못 더해지는
폭은 100/150/200 °C에서 각각 `1.01860/15.17509/124.00981 MHz`다.

또 `pump_rabi_from_power()`는 모든 Rb/Cs D1/D2에 한 Rb-85 D1 `I_sat`를 쓴다
(`gabes/constants.py:22-37`, `gabes/species.py:273-283`). 0.5 mW, 1 mm에서 오늘
확인한 여섯 조합의 `Omega/Gamma`는 모두 `1.883984`였다. 따라서 power trend는
유용하지만 species/line 사이의 absolute saturation 비교에는 부적합하다.

## 기존 개선안의 계산비용과 물리 보존성

| 기존 제안 | 상태 | 계산비용 | 물리 보존성 / 판단 |
|---|---|---:|---|
| sub-Doppler interpolation/status gate | done | O(N) post-process | OBE 불변, under-resolved 숫자 과신 방지. 유지 |
| CSV relative/absolute calibration + sweep branch | done | O(N) import | OBE 불변, 실험 의미 복구. 유지 |
| paraffin population-memory MVP | done, opt-in | 오늘 401점에서 약 `3.16×` | 실제 ground-F memory를 추가하지만 기본 off라 legacy physics 보존 |
| gas/species/line 충돌 coefficient + pressure shift | P1 ready, effort large | O(1) lookup, 동일 solve 차원 | zero-pressure anchor를 보존하면서 잘못된 fallback과 pump/probe 불일치를 교정 가능 |
| added width의 pump optical dephasing 반영 | 위 P1 일부 | 동일 Liouvillian 크기·solve 수 | buffer-cell 내부 일관성을 고치는 negligible-overhead physics change |
| species/line별 Rabi/saturation provenance | P1 ready | O(1) | pump-off anchor 불변, pumped 결과는 의도적으로 교정 |
| AutoOD `<0.1 %` regression 강화 | 과거 제안, 미완료 | 추가 solve 0 | 문서 약속만 더 강하게 보호 |
| calibrated coated-cell transport/dark-flight | parked, research | scan×velocity/time quadrature 가능 | cell geometry·측정 trace 없이는 무거운 자유도만 늘어남 |
| full VCC kernel | parked, research | velocity classes 결합으로 높음 | 특정 cell/kernel/held-out dataset/runtime budget가 있을 때만 opt-in reference mode가 타당 |
| full mF/polarization solve | 현행 범위 밖 | density-matrix 차원 급증 | absolute contrast 요구와 held-out trace가 있을 때 별도 설계 |

핵심은 개발 effort와 runtime overhead를 분리하는 것이다. 충돌 coefficient API는
cross-scheme 작업이라 개발 범위는 크지만 실제 계산은 lookup과 같은 크기의
Liouvillian이므로 거의 공짜다. 반대로 coated transport와 full VCC는 가정을
명시하고 검증 데이터를 확보하기 전에는 “더 정교해 보이는” 것 이상의 가치가 없다.

## 오늘 확인한 저비용 trust 개선

기존 제안이 충분하므로 새 물리 knob를 임의로 늘릴 필요는 없다. 다만 새 paraffin
경로에서 계산비용 없이 바로 닫을 수 있는 presentation/provenance gap은 있다.

1. `raw`에는 `paraffin_coated`와 `paraffin_t1_s`가 있지만
   (`gabes/schemes/sas.py:230-236`), figure title과 metric/table은 이를 전혀 쓰지
   않는다 (`gabes/schemes/sas.py:386-464`). Exported figure에 `coated,
   population-only, T1 reference` status를 넣는 것은 O(1)이며 physics를 바꾸지 않는다.
2. User Guide의 OD/SAS 절은 새 paraffin 옵션을 설명하지 않고 6월 8일 uncoated
   stacked PNG만 보여 준다 (`docs/Userguide/GABES_User_Guide_v2.html:552-584`).
   README의 quasi-static 한계와 CSV calibration/hysteresis 계약을 옮기고 current UI
   figure를 재생성해야 한다. Interactive overhead는 0이다.
3. Coated end-to-end test에 최소한 model/provenance label, fixed-window grid
   convergence, T1/cycle-rate sensitivity ledger를 추가할 수 있다. 이는 물리 검증을
   대신하지 않지만 silent regression과 과도한 claim을 막는다. Cell-specific held-out
   spectrum이 생기기 전에는 `semi-quantitative` status를 올리지 않아야 한다.

## 동작을 바꾸지 않는 순수 코드 최적화

### 1. Pump-off exact analytic population fast path — 여전히 최우선

Pump가 0이면 모든 transition의 population factor는 수학적으로 정확히 `w=1`인데,
현재도 Liouvillian과 `_pump_pops()` table을 만든다
(`gabes/schemes/sas.py:248-305`, `gabes/schemes/sas.py:704-711`).

85Rb D1, 90 °C, 12.5 mm, 1401 points에서 오늘 5회 warm median을 비교했다.

| 경로 | runtime |
|---|---:|
| 현재 | `0.26489 s` |
| analytic `w=1` | `0.08671 s` |
| speedup | **`3.05×`** |

Scan은 bitwise exact였고 최대 상대 alpha 차이는 `5.63×10⁻16`이었다.
`pump_power_mw <= 0`에만 엄격히 적용하고 AutoOD/natural-Rb/buffer/temperature
회귀로 보호하면 public behavior를 바꾸지 않는 안전한 최적화다. Coating-on +
pump-off도 이미 legacy path를 쓰므로 그대로 이득을 받는다.

### 2. Duplicate-heavy CSV median merge vectorization

`_sort_and_merge()`는 duplicate group마다 Python에서 `np.median`을 호출한다
(`gabes/experimental_csv.py:506-520`). 500,000 rows, 각 x가 5회 반복되는 fixture에서
`(x,y)` lexicographic sort 뒤 중앙 index를 벡터로 읽는 후보를 비교했다.

| 경로 | runtime |
|---|---:|
| 현재 group loop | `1.98035 s` |
| vector median indices | `0.17748 s` |
| speedup | **`11.16×`** |

Unique x와 merged y는 exact equal이었다. Odd/even group, signed zero, stable ordering,
extreme finite value 회귀를 추가한 뒤 적용할 가치가 높다.

### 3. Coated basis-conditioned solve의 shared operator / multi-RHS

각 source ground-F의 reset dissipator는 공통 `L_base - gamma_t I`에 source-dependent
trace-one injection만 다르게 쓸 수 있다. 두 Liouvillian을 따로 만들고 푸는 현재
`_basis_reset_pump_pops()` (`gabes/schemes/sas.py:609-618`) 대신 공통 real-basis
matrix에 두 RHS를 한 번에 풀어 본 후보는 해당 stage를 `0.12159 -> 0.09901 s`
(**1.23×**)로 줄였다. Population 최대 차이는 `3.65×10⁻12`였다.

전체 compute 이득은 작으므로 3순위다. Trace-row 처리와 source injection 항등식을
독립 테스트로 고정한 뒤에만 적용해야 한다.

### 4. 적용하지 말아야 할 후보

Coated transfer와 모든 level population을 한 대형 gathered array로 보간해 공통
search를 재사용하는 순진한 vectorization은 오늘 stage를 `0.14355 -> 0.17779 s`로
오히려 느리게 했고 메모리도 늘렸다. Chunked 설계와 end-to-end benchmark 없이
적용할 이유가 없다. Running-median vectorization은 이미 구현돼 있다
(`gabes/lineshape.py:166-173`).

## 문서와 정적 예제 평가

README는 pump-off limit, AutoOD normalization, paraffin population-only/T1 경계,
CSV relative/absolute calibration을 현재 코드와 대체로 정확히 설명한다
(`README.md:87-103`, `README.md:280-317`). In-app `info()`도 paraffin model이
calibrated geometry나 coherence model이 아님을 명시한다
(`gabes/schemes/sas.py:132-169`).

반면 `od.png`와 `sas.png`는 2026-06-08 생성본이며 Transmission과 OD를 한 이미지에
쌓은 이전 layout이다. 현재 코드는 두 개의 single-panel `figure_views`를 만든다
(`gabes/schemes/sas.py:395-422`). 정적 SAS 예제에는 resolution status,
samples/FWHM, coating state, cache/git/model provenance가 없다. User Guide는 “모두
GABES가 직접 계산한 실제 스펙트럼”이라고 말하지만 (`docs/Userguide/
GABES_User_Guide_v2.html:520-524`), 현재 code version의 재생성 가능성을 제공하지
않는다.

따라서 정적 예제는 **uncoated qualitative line-shape illustration**로는 유효하지만,
현재 UI나 coated-cell reference의 재현 가능한 예제로는 부적합하다.

## 테스트 결과와 검증 범위

Scheme 1 관련 테스트:

```text
python -m pytest -q tests/test_sas.py tests/test_experimental_csv.py \
  tests/test_absorption.py tests/test_headless_observables.py \
  tests/test_schemes_render.py

101 passed in 25.18s
```

전체 저장소 회귀:

```text
MPLBACKEND=Agg python -m pytest -q

487 passed in 315.65s (0:05:15)
```

테스트는 atomic-data algebra, AutoOD anchor, 49/25 manifold ratio, CG branching,
hyperfine-pumping trend, resolution gate, CSV calibration/sweep branch, coated transfer
probability와 legacy invariance를 잘 보호한다. 다만 테스트 통과는 omitted
mF/polarization/collision/transport physics 또는 외부 실험 검증을 뜻하지 않는다.

## 최종 우선순위

1. **P1, runtime 거의 0:** `(gas, species, line)` 충돌 coefficient provenance와
   unsupported status를 만들고 pressure shift와 elastic optical dephasing을 분리해
   pump/probe width convention을 일치시킨다.
2. **P1, O(1):** species/line/polarization별 Rabi/saturation scale을 독립 근거와
   함께 도입한다. Pump-off 85Rb D1 AutoOD는 frozen anchor로 유지한다.
3. **P2, solve 추가 0:** 공개 `<0.1 %` AutoOD 약속에 맞춰 regression tolerance를
   강화하고, coated state/T1/`population-only` status를 figure와 metric에 남긴다.
4. **순수 성능:** pump-off analytic fast path를 먼저, duplicate-heavy CSV vector
   merge를 다음으로 구현한다. Coated multi-RHS solve는 낮은 우선순위다.
5. **P2/P3, runtime 0:** current UI로 OD/SAS PNG를 재생성하고 User Guide에 paraffin
   한계와 CSV calibration/hysteresis 계약, 생성 조건·git hash를 기록한다.
6. **검증 데이터가 생길 때만 heavy:** coated transport/dark-flight, full VCC,
   full mF/polarization은 cell geometry, sourced rates, scan protocol, calibration/holdout
   trace, runtime budget이 정해진 뒤 opt-in reference mode로 설계한다.

종합하면 Scheme 1은 실제 실험에 유용한 물리를 이미 많이 갖췄다. 현재 가장 큰
개선은 solver를 무겁게 하는 것이 아니라 **충돌·Rabi provenance, coated-result
provenance, 공개 검증 계약**을 닫고 exact fast path로 속도를 회수하는 것이다.
