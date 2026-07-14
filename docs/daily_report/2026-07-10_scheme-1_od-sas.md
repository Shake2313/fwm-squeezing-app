# 2026-07-10 Scheme 1 Review: OD / SAS

## 선택 기준과 현재 스킴 순서

- 오늘 현지 날짜는 `2026-07-10`이고, day-of-month는 `10`이다.
- 따라서 `n = (day mod 5) + 1 = (10 mod 5) + 1 = 1`이다.
- 현재 등록 순서는 `OD / SAS -> Lambda coherence -> Rydberg-EIT electrometry -> Hanle / EIA / NMOR -> FWM`이다 (`gabes/schemes/__init__.py:19-24`, `README.md:10-16`).
- 오늘 대상은 1번째 스킴인 `SASScheme`, UI 제목 `Absorption spectroscopy (OD / SAS)`이다 (`gabes/schemes/sas.py:53-62`).

## 검토한 근거

- 스킴 등록과 순서: `gabes/schemes/__init__.py:19-24`
- OD/SAS 파라미터 표면: `gabes/schemes/sas.py:67-106`
- species 기반 OD/SAS 계산 경로: `gabes/schemes/sas.py:157-246`
- generic hole-burning fallback: `gabes/schemes/sas.py:248-312`
- headless/figure observables와 lock slope metric: `gabes/schemes/sas.py:318-401`, `gabes/schemes/sas.py:429-445`
- isotope/line data와 hyperfine manifold: `gabes/species.py:182-300`, `gabes/species.py:338-370`
- Ne buffer broadening의 현재 구현과 TODO: `gabes/constants.py:81-89`, `gabes/schemes/sas.py:84-87`
- 회귀/물리 테스트: `tests/test_sas.py:95-179`, `tests/test_absorption.py:57-105`, `tests/test_headless_observables.py`
- 기존 개선안: `docs/checklist.json:22-28`, `docs/checklist.json:29-36`, `docs/checklist.json:88-95`, `docs/checklist.json:102-108`
- 이전 OD/SAS 리뷰: `docs/daily_report/2026-06-25_scheme-1_od-sas.md`

## 현재 정의의 물리 구조

이 스킴은 pump power 하나로 같은 매질 위에서 OD와 SAS를 연결한다. `pump_power_mw = 0`이면 weak-probe linear absorption이고, pump를 올리면 같은 hyperfine manifold에서 velocity-selective hole burning과 hyperfine optical pumping이 생긴다 (`gabes/schemes/sas.py:67-70`, `gabes/schemes/sas.py:157-179`).

현실적인 species 모드에서는 `species.build_manifold()`가 isotope/line별 `{Fg}->{Fe}` hyperfine manifold를 만들고, CG-branched decay와 transit relaxation을 포함한다 (`gabes/species.py:338-370`). 각 isotope component마다 vapor density, self-broadening, Ne broadening, pump Rabi가 계산되고 (`gabes/schemes/sas.py:169-179`), pump steady-state population table을 만든 뒤 probe Lorentzian과 Maxwell velocity weight로 absorption coefficient를 합산한다 (`gabes/schemes/sas.py:205-246`).

이 점은 실험물리학자 관점에서 중요하다. OD와 SAS가 서로 다른 toy plot 두 개가 아니라, 같은 atomic data와 같은 cell parameter에서 pump를 껐다 켜며 이어진다. 85Rb D1 pump-off limit가 기존 AutoOD 검증용 `ODScheme`과 integrated/peak scale에서 1% 안에 맞도록 테스트되어 있고 (`tests/test_sas.py:95-105`), 85Rb D1 ground-manifold strength ratio가 `49/25` 근처임도 테스트된다 (`tests/test_sas.py:108-113`).

## 오늘 실행한 확인

직접 실행한 관련 테스트:

```text
python -m pytest tests/test_sas.py tests/test_absorption.py tests/test_headless_observables.py -q
28 passed in 9.87s
```

추가로 몇 가지 대표 조건을 직접 계산했다.

| 조건 | compute | headless observables | 주요 결과 |
|---|---:|---:|---|
| 85Rb D1 OD, 90 C, 12.5 mm | 0.265 s | 0.038 s | Peak OD 6.58, sub-Doppler 없음 |
| 85Rb D2 SAS, 30 C, 50 mm, 1.5 mW | 0.454 s | 0.023 s | Peak OD 0.45, narrowest sub-Doppler 14.4 MHz |
| 같은 SAS, transit 2000 kHz | 0.460 s | 0.023 s | crossover `T(1.719 GHz)=0.698` |
| 같은 SAS, transit 100 kHz | 0.439 s | 0.023 s | crossover `T(1.719 GHz)=0.891` |
| 같은 SAS, transit 20 kHz | 0.480 s | 0.024 s | crossover `T(1.719 GHz)=0.964` |
| 85Rb D1 OD, 50 C, Ne 0 Torr | 0.254 s | 0.023 s | Peak OD 0.32 |
| 85Rb D1 OD, 50 C, Ne 50 Torr | 0.193 s | 0.023 s | Peak OD 0.25 |

이 결과는 기대한 물리 방향과 맞다. Pump-off에서는 sub-Doppler readout이 사라지고, pump-on에서는 10 MHz대 Doppler-free feature가 나타난다. Transit relaxation을 낮추면 atoms leaving/entering rate가 줄어 optical pumping이 더 강하게 누적되어 crossover transmission이 증가한다. Ne pressure는 현재 모델에서 homogeneous broadening만 키우므로 peak OD가 낮아지고 lock slope가 완만해진다.

## 실험물리학적 평가

현재 OD/SAS 스킴은 **실제 분광 실험을 계획하고 해석하는 semi-quantitative reference**로 쓸 만하다. 특히 다음 항목들은 실험자가 바로 믿고 사용할 수 있는 쪽에 가깝다.

- isotope/line 선택이 실제 Rb/Cs D-line hyperfine data와 연결되어 있다 (`gabes/schemes/sas.py:72-79`, `gabes/species.py:182-249`).
- pump power와 waist가 pump Rabi로 연결된다 (`gabes/schemes/sas.py:88-90`, `gabes/species.py:273-280`).
- pump-off OD absolute scale이 AutoOD 검증 경로와 같은 정규화를 따른다 (`gabes/schemes/sas.py:205-245`, `gabes/species.py:249-270`, `tests/test_sas.py:95-105`).
- SAS에서 단순 Lamb dip뿐 아니라 hyperfine optical pumping이 crossover enhancement/inversion을 만든다는 핵심 서명이 들어 있다 (`tests/test_sas.py:127-141`, `README.md:187-198`).
- `Lock slope`와 `Lock detuning` metric이 있어 laser-lock discriminator 관점의 1차 실험 readout도 제공한다 (`gabes/schemes/sas.py:429-445`, `tests/test_sas.py:161-168`).

다만 **정량 metrology reference**라고 부르기에는 아직 제한이 분명하다.

- Buffer gas는 `neon_buffer_broadening()` 하나로 homogeneous FWHM만 더한다. Pressure shift, gas/species/line별 coefficient, Dicke narrowing은 없다 (`gabes/constants.py:81-89`, `gabes/schemes/sas.py:84-87`).
- Full velocity-changing collision은 없다. 현재는 velocity classes를 Maxwell weight로 독립 평균하므로, collision-rich buffer cell에서 velocity redistribution까지 맞추는 모델은 아니다 (`gabes/schemes/sas.py:218-245`).
- Polarization, Zeeman sublevel imbalance, optical pumping anisotropy, pump/probe alignment asymmetry, etalon/RAM/background 같은 실제 SAS setup distortion은 모델 밖이다. `docs/checklist.json`도 이를 full manifold가 아닌 low-order proxy 후보로 남겨 두고 있다 (`docs/checklist.json:88-95`).
- Transit relaxation은 실험적으로 유용한 phenomenological knob지만, beam profile, diffusion, wall collision을 분해한 model은 아니다 (`gabes/schemes/sas.py:91-94`, `gabes/species.py:338-356`).

따라서 연구자가 "이 isotope/line/cell에서 어느 hyperfine feature가 보이고, pump를 올리면 crossover가 어느 방향으로 바뀌며, lock slope 후보가 어디쯤인가"를 보는 용도로는 좋다. 반면 buffer-gas cell의 absolute line center, pressure-shifted lock point, Dicke-narrowed linewidth budget을 publication-grade fit parameter로 바로 쓰기에는 부족하다.

## 기존 개선안과 계산 부하 분석

### 1. Buffer-gas pressure shift + Dicke narrowing

체크리스트의 `buffer-gas-pressure-shift`가 OD/SAS에 가장 직접적으로 해당한다 (`docs/checklist.json:22-28`). 현재 `ne_pressure_torr`는 `buffer_gamma`로만 들어가고 (`gabes/schemes/sas.py:163`, `gabes/schemes/sas.py:258`), help text도 pressure shift와 Dicke narrowing이 없다고 명시한다 (`gabes/schemes/sas.py:84-87`).

계산 부하는 구현 깊이에 따라 크게 다르다.

- Pressure shift만 넣으면 부하는 거의 없다. Transition center 또는 isotope/line reference offset에 pressure-dependent scalar를 더하면 되므로 solve dimension과 velocity grid가 그대로다.
- Phenomenological Dicke narrowing도 낮은 부하다. `gamma_eff` 또는 Doppler width에 low-order effective correction을 주는 방식이면 현재 table/interpolation 구조를 유지한다.
- Gas/species/line coefficient table도 런타임 비용은 사실상 lookup뿐이다. 코드/테스트/문서 작업은 필요하지만 계산량은 증가하지 않는다.
- Full velocity-changing collision까지 가면 다른 문제다. velocity classes가 결합되므로 현재 separable Maxwell average 구조가 깨지고, 큰 block solve나 iterative redistribution이 필요해질 수 있다. 체크리스트도 별도 GROUP C 항목으로 분리해 둔 판단이 맞다 (`docs/checklist.json:102-108`).

결론: **pressure shift + coefficient table + phenomenological Dicke narrowing**은 물리적 실용성을 꽤 올리면서 interactive runtime을 거의 유지할 수 있는 좋은 다음 단계다. Full VCC는 별도 설계가 필요하다.

### 2. Low-order polarization / Zeeman proxies

`low-order-polarization-zeeman-proxies`도 OD/SAS와 Rydberg에 걸친 남은 제안이다 (`docs/checklist.json:88-95`). OD/SAS에서는 full Zeeman manifold를 풀지 않고도 polarization purity, effective branching selector, asymmetry broadening, line-dependent contrast scale 같은 proxy를 둘 수 있다.

부하는 대체로 작다. 각 transition weight나 effective pumping rate에 scalar modifier를 곱하는 정도라면 현재 `A_t`, `w=(pop_g-pop_e)/p_ground`, `gamma_eff` 조립에 낮은 차수의 보정만 추가하면 된다 (`gabes/schemes/sas.py:234-245`). 다만 모델 해석이 조심스럽다. 실제 polarization physics처럼 보이게 만들수록 validation 없이 fit knob가 늘어날 위험이 있다. 실험 레퍼런스 가치를 올리려면 "Zeeman-resolved solver"가 아니라 "setup imperfection proxy"라고 UI/문서에 명확히 표시하는 편이 좋다.

### 3. Figureless observables와 lock readout

이전 OD/SAS 리뷰에서 코딩 최적화 후보였던 figure 분리와 lock-point readout은 이미 반영되어 있다. `observables(..., include_figures=False)`가 figure 생성을 건너뛰고 (`gabes/schemes/sas.py:318-401`), `tests/test_headless_observables.py`가 이를 확인한다. 오늘 측정에서도 headless observables는 약 `0.02-0.04 s`로 가볍다.

이 개선은 batch review/report 자동화에 실제로 유용하다. Matplotlib figure가 필요 없는 metric/table 경로에서는 현재 구조가 충분히 빠르다.

## 순수 코딩 측면의 속도 개선 후보

물리 기능을 건드리지 않는 범위에서 보면, 큰 병목은 species-mode `compute()`의 pump population table과 transition별 interpolation이다 (`gabes/schemes/sas.py:205-245`).

우선순위가 높은 후보는 다음이다.

1. `species.build_manifold(iso, line, transit_rate=gt)` 결과 중 transit-rate와 무관한 skeleton을 캐시한다. 현재는 isotope/line/transit 조합마다 manifold를 재조립한다 (`gabes/species.py:338-370`, `gabes/schemes/sas.py:169-174`). 다만 transit decay가 atom decay list에 들어가므로 cache boundary를 조심해야 한다.
2. 같은 component 안에서 `pop_at` interpolation을 level별로 미리 만들고 재사용하는 현재 구조는 괜찮다 (`gabes/schemes/sas.py:229-232`). 추가 개선은 `deff_grid.ravel()`과 `probe_base` 같은 큰 배열을 component별로 한 번만 만들고, transitions loop 내부 allocation을 더 줄이는 방향이다.
3. Natural Rb처럼 isotope component가 여러 개인 경우 marker string/table은 isotope/line 기준으로 거의 정적이다 (`gabes/schemes/sas.py:187-202`, `gabes/schemes/sas.py:374-379`). compute cache가 이미 Streamlit에 있지만, headless batch에서 반복 호출이 많다면 marker/table building을 작은 pure helper로 캐시할 수 있다.
4. `scan_points`가 1401이고 compute가 0.2-0.5 s 수준이라 일반 UI에서는 급하지 않다. 최적화가 필요해지는 곳은 parameter sweep이다. 이 경우 figureless path를 쓰고, pressure/cell/line_strength처럼 solve를 다시 할 필요가 없는 knobs를 recompute tier에서 빼는 현재 설계를 유지하는 것이 가장 큰 이득이다 (`CLAUDE.md:16-18`, `streamlit_app.py:460-477`).

## 종합 판단

OD/SAS는 현재 GABES 안에서 가장 "실험실 냄새가 나는" 스킴 중 하나다. Pump-off OD가 검증된 absolute absorption scale에 닿아 있고, pump-on SAS가 hyperfine optical pumping과 crossover enhancement를 구현한다. 실험자가 saturated absorption setup을 조정하거나 lock 후보를 찾는 데는 충분히 유용한 semi-quantitative reference다.

다음 물리 개선은 full VCC나 full Zeeman으로 바로 뛰기보다, 계산 부하가 거의 없는 **buffer gas coefficient table + pressure shift + phenomenological Dicke narrowing**이 가장 효율적이다. 그다음에는 명확히 phenomenological이라고 표시한 low-order polarization/Zeeman proxy가 실험 setup mismatch를 설명하는 데 도움이 될 수 있다.
