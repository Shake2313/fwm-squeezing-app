# 2026-08-25 Scheme 1 물리 리뷰 — OD / SAS

## 선택 규칙과 현재 다섯 scheme

- 로컬 날짜는 `2026-08-25`, day-of-month는 `25`이다.
- `n = (day mod 5) + 1 = (25 mod 5) + 1 = 1`이므로 오늘의 대상은
  **Scheme 1**이다.
- UI 순서는 `gabes/schemes/__init__.py`의 `_SCHEMES` 리스트가 직접 정한다
  (`gabes/schemes/__init__.py:12-25`, `gabes/schemes/__init__.py:37-39`).

| 번호 | 현재 scheme | registry 인스턴스 | 런타임 이름 / 제목 |
|---:|---|---|---|
| 1 | OD / SAS | `SASScheme()` | `sas` / `Absorption spectroscopy (OD / SAS)` |
| 2 | Lambda coherence (EIT / AT / CPT) | `LambdaScheme()` | `lambda` / `Λ coherence (EIT / AT / CPT)` |
| 3 | Rydberg-EIT electrometry | `RydbergEITScheme()` | `rydberg_eit` / `Rydberg-EIT electrometry` |
| 4 | Hanle / EIA / NMOR | `MagnetoScheme()` | `magneto` / `Magneto-optics (Hanle/MOR)` |
| 5 | FWM | `FWMScheme()` | `fwm` / `Four-wave mixing (Gain diagnostic / Biphoton)` |

Scheme 1의 사용자용 정의는 README의 첫 행과 일치한다
(`README.md:8-16`, `gabes/schemes/sas.py:53-63`).

## 먼저 검색한 기존 제안·TODO·issue note

물리 판단 전에 다음을 검색하고 현재 코드와 다시 대조했다.

- 저장소 지침과 개요: `AGENTS.md`, `CLAUDE.md`, `README.md`
- 현재 작업 registry: `docs/checklist.json`
- 코드 TODO/FIXME/XXX: `gabes/constants.py:79-89`에 OD/SAS와 직접 관련된 TODO
  하나가 있다. Ne broadening 상수를 gas/species/line 계수표와 pressure shift로
  승격하라는 내용이다.
- 과거 Scheme 1 보고서 여섯 건:
  `2026-06-25`, `2026-07-10`, `2026-07-20`, `2026-07-30`, `2026-08-05`,
  `2026-08-10`
- 구현: `gabes/schemes/sas.py`, `gabes/species.py`, `gabes/constants.py`,
  `gabes/core.py`, `gabes/doppler.py`, `gabes/lineshape.py`,
  `gabes/experimental_csv.py`, `streamlit_app.py`
- 테스트: `tests/test_sas.py`, `tests/test_experimental_csv.py`,
  `tests/test_absorption.py`, `tests/test_headless_observables.py`,
  `tests/test_schemes_render.py`
- 문서·예제: `README.md:8-16`, `README.md:87-103`, `README.md:280-307`,
  `docs/Userguide/GABES_User_Guide_v2.html:520-584`,
  `docs/Userguide/userguide_assets/od.png`,
  `docs/Userguide/userguide_assets/sas.png`,
  `references/AutoOD/ReferenceOD.csv`

기존 개선 제안은 **있다**. 현재 checklist가 Scheme 1 관련 제안을 다음처럼
정리하고 있으므로 새 물리 knob를 임의로 제안할 상황은 아니다.

1. 충돌 계수 provenance, pressure shift, pump/probe homogeneous-width 일관성:
   P1, `ready`, 개발 effort `large` (`docs/checklist.json:209-243`).
2. species/line별 Rabi와 saturation convention:
   P1, `ready`, effort `medium` (`docs/checklist.json:246-264`).
3. sub-Doppler linewidth/lock-slope 신뢰성:
   `done` (`docs/checklist.json:267-285`).
4. CSV 절대 보정과 sweep 방향 보존:
   `done` (`docs/checklist.json:288-305`).
5. full velocity-changing-collision(VCC) kernel:
   P2 연구 항목, `parked` (`docs/checklist.json:563-582`).

`rg --files`로 검색한 결과 별도 로컬 issue/proposal/note 파일이나 `.github`
issue template은 발견되지 않았다. 과거의 자유로운 low-order polarization/Zeeman
proxy는 측정 없이 fit 자유도만 늘릴 위험 때문에 checklist에서 기각됐고, 검증
데이터 없는 phenomenological Dicke knob도 충돌 계수 P1의 범위 밖이다. 이 판단은
실험 reference의 식별 가능성을 지키는 방향으로 타당하다.

## 결론 요약

**현재 OD/SAS는 실제 원자물리를 담은 유용한 warm-vapor 실험 계획용
semi-quantitative reference다.** 실제 Rb/Cs hyperfine 상수, 질량, 자연 선폭,
Wigner-6j/3j 선세기, Maxwell 속도 평균, counter-propagating pump/probe,
CG-branched spontaneous decay, transit reset, hyperfine optical pumping을 사용한다
(`gabes/schemes/sas.py:1-37`, `gabes/species.py:127-182`,
`gabes/species.py:336-437`). Pump-off 85Rb D1 절대 OD와 pump-on
enhanced/inverted crossover가 같은 모델에서 연속적으로 나오는 점은 단순 교육용
Lorentzian 합보다 훨씬 낫다.

적합한 용도는 다음과 같다.

- pure-Rb 또는 저압 셀의 hyperfine line finding과 isotope assignment
- pump-off OD scale sanity check
- pump power, 온도, cell length, transit rate 변화의 방향성 확인
- resolution status가 `resolved`인 조건에서의 lock 후보 위치와 local DC slope 탐색
- 절대 보정 근거가 있는 CSV와의 transmission 비교

반대로 다음에는 아직 정량 reference로 부적합하다.

- buffer-cell SAS의 절대 선중심·선폭·saturation/hole-burning
- 고온 Cs와 natural-Rb mixture의 self-broadened linewidth
- species/line/polarization 사이의 절대 saturation-power 비교
- mF·편광·Zeeman-resolved contrast와 branching
- 변조·복조 전달함수를 포함한 실제 servo error signal
- held-out 측정으로 검증된 pumped-SAS 절대 contrast와 linewidth

따라서 연구자가 이 코드를 인용할 때 가장 정확한 표현은
**“AutoOD-anchored pump-off OD와 hyperfine-pumping SAS를 연결한 빠른 정성/반정량
실험 계획 모델”**이다. 정밀 collision metrology나 lock electronics calibration
reference로 부르는 것은 과하다.

## 현재 구현이 담는 실제 물리

### 1. pump-off OD와 pump-on SAS가 하나의 모델이다

laser scan을 `Δ`, 원자 속도를 `v`라 하면 pump와 probe detuning은 각각
`Δ + kv`, `Δ - kv`로 들어간다 (`gabes/schemes/sas.py:4-30`). Pump가 0이면
population factor가 1로 줄어 unit-area Lorentzian의 Maxwell 평균, 즉 multi-line
Voigt OD가 된다. Pump를 켜면 velocity-selective saturation과 hyperfine optical
pumping이 같은 스펙트럼 위에 Lamb dip과 crossover를 만든다
(`gabes/schemes/sas.py:211-253`).

실험 노브는 온도, cell length, isotope/natural abundance, D1/D2, pump power,
waist, transit relaxation, Ne pressure다 (`gabes/schemes/sas.py:65-107`).
`cell_mm`와 `line_strength`는 navigate-only 후처리 knob여서 매번 OBE를 다시 풀지
않는다 (`gabes/schemes/sas.py:82-99`, `gabes/schemes/sas.py:324-331`).

### 2. hyperfine optical pumping은 실제로 구현되어 있다

Rb-85, Rb-87, Cs-133의 hyperfine A/B, centroid, mass, abundance, 자연 선폭이
명시돼 있다 (`gabes/species.py:127-182`). 허용 `Fg→Fe` 선세기는 Wigner-6j/3j에서
계산하고, 같은 angular factor가 absorption weight, relative pump Rabi,
spontaneous branching에 일관되게 들어간다 (`gabes/species.py:226-270`,
`gabes/species.py:366-414`).

`build_manifold()`는 허용 F-state 전체와 CG-branched decay를 만들고 transit
relaxation으로 population을 thermal ground distribution으로 되돌린다
(`gabes/species.py:337-422`). 낮은 transit rate에서 crossover transmission이
커지는 테스트는 단순 선형 line sum이 아니라 optical pumping이 작동한다는 유용한
내부 검증이다 (`tests/test_sas.py:144-158`).

다만 각 F는 하나의 lumped state다. mF, pump/probe polarization, magnetic field,
beam profile, diffusion은 없다. 따라서 hyperfine 중심과 crossover의 질적 구조는
유용하지만 absolute pumped contrast와 branching은 실험별 기준값이 아니다.

### 3. experimental CSV 경로는 의미 경계를 잘 지킨다

CSV 기본 extrema scaling은 명시적으로 `Relative normalized transmission`이며,
dark/reference 또는 gain/offset이 제공될 때만 `Absolute transmission`을 만든다
(`gabes/experimental_csv.py:1-29`, `gabes/experimental_csv.py:754-891`).
절대 보정에서는 0–1 밖 sample도 숨기지 않고 경고와 함께 보존한다
(`gabes/experimental_csv.py:838-890`). Sweep reversal은 정렬 전에 세고 원래
forward/reverse branch를 보존하며, 호환 merged trace가 hysteresis를 잃는다고
명시한다 (`gabes/experimental_csv.py:212-321`,
`gabes/experimental_csv.py:441-520`).

합성 절대 보정은 원래 transmission을 `rtol=1e-9`로 복원하고, round-trip fixture는
5 forward/5 reverse sample을 정확히 보존한다
(`tests/test_experimental_csv.py:100-159`,
`tests/test_experimental_csv.py:223-247`). 이는 실험 데이터와 모델을 비교할 때
“상대 모양 맞춤”과 “절대 투과 검증”을 구분하게 해 주는 중요한 장점이다.

## 오늘 재확인한 정량 상태

### 1. 기본 sub-Doppler readout은 정직하지만 under-resolved다

현재 natural-Rb D1 SAS 기본값은 `P=0.5 mW`, `40 °C`, `75 mm`, 1401 points다.
오늘 headless compute는 다음을 냈다.

- hero: `SAS resolution = resolution-limited`
- provisional sub-Doppler FWHM: `21.56 MHz`
- samples/FWHM: `2.7`
- Peak OD: `0.54`
- calculated single-line Gaussian Doppler FWHM: `518.7 MHz`

즉 숫자를 아예 숨기지는 않지만 신뢰 가능한 linewidth/lock slope로 승격하지 않는
현재 동작이 옳다 (`gabes/schemes/sas.py:371-390`,
`gabes/lineshape.py:141-250`).

checklist와 같은 85Rb D1 fixed-window refinement를 재현한 결과는 다음과 같다.
시간은 오늘 환경의 단일 warm 실행값이며 상대 비교용이다.

| points | compute | 상태 | interpolated FWHM | samples/FWHM | local slope |
|---:|---:|---|---:|---:|---:|
| 1401 | `0.535 s` | resolution-limited | `18.7049 MHz` | `3.73` | 숨김 |
| 2801 | `1.001 s` | resolved | `17.8414 MHz` | `7.11` | `0.0379 /MHz` |
| 5601 | `1.926 s` | resolved | `17.7409 MHz` | `14.14` | `0.0394 /MHz` |

2801→5601에서 FWHM 차이는 약 `0.56 %`, slope 차이는 약 `3.81 %`다.
Reference-only 5601점 계산을 interactive default로 강제할 필요는 없고, 현재처럼
O(N) status와 opt-in refinement로 신뢰도를 구분하는 편이 비용 대비 좋다
(`tests/test_sas.py:244-294`).

### 2. pump-off AutoOD 값은 좋지만 regression 계약이 문서보다 느슨하다

85Rb D1, 90 °C, 12.5 mm, pump off에서 internal AutoOD primitive와 다시 비교했다.

| 비교량 | 현재 상대 차이 |
|---|---:|
| integrated absorption area | `-0.039317 %` |
| peak absorption | `-0.012573 %` |

현재 값은 README와 User Guide의 `<0.1 %` 약속을 만족한다
(`gabes/schemes/sas.py:123-151`, `README.md:280-297`,
`docs/Userguide/GABES_User_Guide_v2.html:409-412`). 그러나 직접 보호하는
`test_pump_off_reproduces_autood_85rb_d1()`은 area/peak 모두 `<1 %`만 요구한다
(`tests/test_sas.py:101-111`). 같은 solve의 assertion을 `1e-3` 이하로 강화하면
추가 계산비용 없이 공개 계약을 실제로 보호한다.

### 3. collision model은 여전히 P1 경계다

`gamma_eff = gamma_nat + self + Ne broadening`은 계산하지만
(`gabes/schemes/sas.py:163-185`), realistic species pump Liouvillian은 natural-Γ로
만든 manifold를 그대로 사용한다 (`gabes/species.py:337-399`,
`gabes/schemes/sas.py:221-227`). `gamma_eff`는 pump population table의 범위/간격과
probe Lorentzian에만 직접 들어간다 (`gabes/schemes/sas.py:225-252`).

예를 들어 85Rb D1, Ne 20 Torr에서 probe FWHM convention은
`5.75 + 3.91×20 = 83.95 MHz`지만 pump spontaneous/dephasing dynamics는 natural
5.75 MHz manifold다. 따라서 pressure-broadened hole burning과 saturation이 내부적으로
일관되지 않는다. Pressure shift, Dicke narrowing, VCC도 없다
(`gabes/constants.py:79-89`, `gabes/schemes/sas.py:84-87`).

또 `BETA_SELF`는 Rb 계수라고 주석에 명시되어 있지만
(`gabes/species.py:49-50`), `self_broadened_gamma(iso, N)`는 `iso`를 사용하지 않고
모든 원소에 같은 `BETA_SELF*N`을 적용한다 (`gabes/species.py:211-213`). Cs slider
범위에서 잘못 더해지는 `β_Rb N_Cs / 2π`는 다음과 같다.

| Cs 온도 | 현재 추가 폭 |
|---:|---:|
| 22 °C | `0.002186 MHz` |
| 30 °C | `0.004776 MHz` |
| 100 °C | `1.018599 MHz` |
| 150 °C | `15.175086 MHz` |
| 200 °C | `124.009810 MHz` |

기본 Cs preset에서는 작지만 허용 slider 상단에서는 자연 선폭보다 훨씬 커진다.
species/line/gas별 sourced coefficient와 unsupported status가 필요하다.

### 4. 모든 species/line이 하나의 Rb 기준 saturation intensity를 쓴다

`I_SAT = 4.484 mW/cm²` 한 값이 constants에 있고 (`gabes/constants.py:22-37`),
`pump_rabi_from_power()`가 모든 Rb/Cs D1/D2에서 이를 재사용한다
(`gabes/species.py:273-283`). `P=0.5 mW`, `w=1 mm`에서 확인한 `Ω/Γ`는 여섯
species/line 조합 모두 정확히 `1.883984`였다.

따라서 pump-power trend는 유용하지만 종·D-line·편광 사이의 absolute saturation
power 비교에는 근거가 부족하다. Checklist의 P1대로 dipole, linewidth,
polarization, intensity, angular/cyclic-frequency convention을 기록한 lookup/derived
scale을 써야 한다 (`docs/checklist.json:246-264`).

## 기존 개선안의 계산비용과 물리 보존성

| 기존 제안 | 현재 상태 | 런타임 비용 | 물리 보존성 / 판단 |
|---|---|---:|---|
| sub-Doppler interpolation/status/hero gate | 완료 | O(N) post-process | OBE 불변, 과신 방지. 유지 |
| CSV 상대/절대 calibration + sweep branch | 완료 | O(N) import | OBE 불변, 실험 의미 복구. 유지 |
| gas/species/line 충돌 계수 + pressure shift | P1 ready, 개발 effort large | O(1) lookup, 같은 solve 차원·횟수 | sourced coefficient와 zero-pressure 회귀가 있으면 pure-cell anchor를 보존하면서 잘못된 fallback을 교정 |
| added width를 pump optical dephasing에 반영 | 위 P1의 일부 | 같은 Liouvillian 크기·solve 수 | 현재 buffer-cell 내부 불일치를 고치는 negligible-overhead physics change |
| species/line별 Rabi/saturation provenance | P1 ready, effort medium | O(1) | pump-off AutoOD 불변; pumped 결과는 의도적으로 교정 |
| AutoOD `<0.1 %` 테스트 강화 | 과거 제안, 미완료 | 추가 solve 0 | 문서의 현재 약속만 더 강하게 보호 |
| phenomenological Dicke knob | 1차 범위 제외 | scalar면 낮음 | 검증 데이터 없이는 non-identifiable fit knob. 현재 제외가 타당 |
| low-order polarization/Zeeman proxy | rejected | 낮을 수 있음 | 독립 측정 없이는 대비 왜곡을 자유롭게 흡수. 재도입 비권고 |
| full VCC kernel | parked | velocity-class coupling으로 높음 | 특정 cell/kernel/dataset/runtime이 정해진 opt-in reference mode만 타당 |
| full mF/polarization solve | 현행 범위 밖 | density-matrix 차원 급증 | absolute contrast 요구와 held-out trace가 있을 때 별도 설계 |

핵심은 **개발 effort와 interactive runtime overhead를 구분하는 것**이다. 충돌
coefficient API는 cross-scheme provenance 작업이라 구현 범위는 크지만, 계산 자체는
lookup과 기존 크기의 Liouvillian이므로 거의 공짜다. 반대로 full VCC는 속도 클래스의
독립 Maxwell 평균을 깨므로 negligible-overhead 개선이 아니다.

## 동작을 바꾸지 않는 순수 코드 최적화

### 1. pump-off analytic population fast path — solver 측 최우선

Pump가 0이면 각 transition의 population factor는 정확히 `w=1`인데 현재도
Liouvillian, `_pump_pops()` table, level별 interpolation array를 만든다
(`gabes/schemes/sas.py:221-252`, `gabes/schemes/sas.py:505-512`).

85Rb D1, 90 °C, 12.5 mm, 1401 points에서 같은 scan과 analytic `w=1` Voigt 합을
warm median 5회로 비교했다.

| 경로 | median runtime |
|---|---:|
| 현재 | `0.54647 s` |
| analytic `w=1` | `0.15648 s` |
| speedup | **`3.49×`** |

scan axis는 bitwise exact였고, 최대 상대 alpha 차이는 `5.63×10⁻16`이었다.
조건을 `pump_power_mw <= 0`으로 엄격히 제한하고 AutoOD, natural-Rb, buffer,
temperature 회귀로 보호하면 public behavior를 바꾸지 않는 안전한 최적화다.

### 2. duplicate-heavy CSV median merge vectorization — importer 최우선

`_sort_and_merge()`는 duplicate x group마다 Python에서 `np.median`을 호출한다
(`gabes/experimental_csv.py:506-520`). 500,000 rows, 각 x가 5회 반복된 synthetic
trace에서 `(x, y)` 정렬 후 group 중앙 index를 벡터로 뽑는 후보를 비교했다.

| 경로 | median runtime |
|---|---:|
| 현재 group loop | `5.50935 s` |
| vector median indices | `0.08122 s` |
| speedup | **`67.83×`** |

unique x는 exact equal, merged y 최대 차이는 0이었다. CSV 상한 근처에서는 실제
체감이 크다. Odd/even group, signed zero, extreme finite value, stable sorting을
회귀로 고정한 뒤 적용할 가치가 높다.

### 3. 이미 완료됐거나 낮은 우선순위인 항목

- `subdoppler_feature()`의 running median은 이미
  `sliding_window_view + np.median(axis=1)`로 vectorized돼 있다
  (`gabes/lineshape.py:166-173`). 과거 최적화가 완료된 상태다.
- Pump-on 경로는 같은 uniform `deff` 축에 대해 각 level마다 `np.interp` search를
  반복하고, `(scan, velocity)` 크기의 `pop_at` 배열을 level별로 보관한다
  (`gabes/schemes/sas.py:237-252`). 공통 interpolation index/weight를 재사용할 수
  있지만, `np.interp`와의 strict equivalence와 peak-memory benchmark가 먼저다.
- `build_manifold()`는 이미 `lru_cache(maxsize=64)`다
  (`gabes/species.py:337-338`). angular-data caching 확대는 우선순위가 낮다.

## 문서와 정적 예제 평가

README는 OD/SAS의 pump-off limit, AutoOD normalization, line-weight convention,
CSV relative/absolute calibration, sweep-direction 계약을 현재 코드와 대체로 정확히
설명한다 (`README.md:87-103`, `README.md:280-307`). User Guide도 Gaussian 폭과
measured Voigt/envelope 폭을 구분하고, `resolved`일 때만 lock slope를 보고한다고
정확히 적는다 (`docs/Userguide/GABES_User_Guide_v2.html:577-583`).

남은 문서 gap은 세 가지다.

1. `od.png`와 `sas.png`는 2026-06-08 생성본으로 Transmission과 OD를 한 이미지에
   쌓은 과거 layout이다. 현재 코드는 두 개의 `figure_views`로 분리한다
   (`gabes/schemes/sas.py:342-369`). SAS PNG에는 현재의 resolution status,
   samples/FWHM, half-height edge도 없다.
2. User Guide OD/SAS 절은 새 CSV relative/absolute calibration과
   forward/reverse branch/hysteresis 계약을 설명하지 않는다. 이 정보는 README에만
   충분히 들어 있다 (`README.md:87-103`).
3. 문서의 AutoOD `<0.1 %` 주장과 `tests/test_sas.py`의 `<1 %` 회귀 허용치가
   일치하지 않는다. 현재 값은 약속을 만족하지만 guardrail은 약하다.

현재 UI로 PNG를 재생성하고 condition, git hash/cache version, resolution status를
caption에 기록하며, README의 CSV 계약을 Guide에 옮기면 추가 solve 없이 실험
reference로서 문서와 실행 화면을 일치시킬 수 있다.

## 테스트 결과와 검증 범위

Scheme 1 관련 테스트:

```text
python -m pytest -q tests/test_sas.py tests/test_experimental_csv.py \
  tests/test_absorption.py tests/test_headless_observables.py \
  tests/test_schemes_render.py

93 passed in 51.08s
```

전체 저장소 테스트:

```text
MPLBACKEND=Agg python -m pytest -q

476 passed in 524.16s (0:08:44)
```

관련 테스트는 atomic-data algebra, AutoOD anchor, 49/25 manifold ratio,
CG branching, hyperfine-pumping trend, resolution gate, absolute CSV calibration,
sweep branch 보존, headless/render contract를 잘 보호한다
(`tests/test_sas.py:58-175`, `tests/test_sas.py:178-362`,
`tests/test_experimental_csv.py:83-270`).

그러나 테스트 통과는 omitted mF/polarization/collision physics나 외부 실험 검증을
대체하지 않는다. Pumped-SAS 검증은 현재 sharp feature, integrated absorption 감소,
transit-rate trend 같은 내부·정성 검사 중심이다 (`tests/test_sas.py:133-158`).
측정 스펙트럼의 absolute contrast/linewidth를 fit subset과 held-out subset으로 나눈
검증은 아직 없다.

## 최종 우선순위

1. **P1, 런타임 거의 0:** `(gas, species, line)` 충돌 coefficient provenance와
   unsupported status를 만들고 pressure shift와 elastic optical dephasing을 분리해
   pump/probe width convention을 일치시킨다.
2. **P1, O(1):** species/line/polarization별 Rabi/saturation scale을 근거와 함께
   도입한다. Pump-off 85Rb D1 AutoOD는 frozen anchor로 유지한다.
3. **P2, 추가 solve 0:** AutoOD 공개 `<0.1 %` 약속에 맞춰 regression tolerance를
   `1e-3` 이하로 강화한다.
4. **순수 성능:** pump-off exact analytic path를 먼저, duplicate-heavy CSV vector
   merge를 다음으로 구현한다. 오늘 환경에서 각각 3.49×와 67.83× 후보 이득이었다.
5. **P2/P3, 계산비용 없음:** current UI로 OD/SAS PNG를 재생성하고 User Guide에 CSV
   calibration/hysteresis 계약을 넣는다.
6. **필요할 때만 heavy:** full VCC와 full mF/polarization solver는 특정 실험 cell,
   sourced kernel, measured trace, held-out validation, runtime budget이 정해진 뒤
   opt-in reference mode로 설계한다.

종합하면 현재 Scheme 1의 가장 강한 부분은 **실제 hyperfine data와 optical pumping을
사용하면서도 under-resolved 결과와 상대 CSV 보정을 과신하지 않게 막는 것**이다.
다음 단계는 solver를 무겁게 만드는 것이 아니라 충돌·Rabi provenance와 공개 검증
계약을 먼저 닫고, exact fast path로 속도를 회수하는 것이다.
