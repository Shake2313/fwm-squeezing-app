# 2026-07-20 Scheme 1 물리 리뷰 — OD / SAS

## 선택 규칙과 현재 다섯 scheme

- 현지 날짜는 `2026-07-20`이고 day-of-month는 `20`이다.
- `n = (day mod 5) + 1 = (20 mod 5) + 1 = 1`이므로 오늘 대상은 첫 번째 scheme이다.
- 실제 UI 순서는 registry의 `_SCHEMES` 리스트가 정하며, README의 표와 일치한다 (`gabes/schemes/__init__.py:19-24`, `README.md:8-16`).

| 번호 | 현재 scheme | 코드 정의 | 핵심 출력 |
|---:|---|---|---|
| 1 | OD / SAS | `SASScheme` | pump-off Doppler OD, pump-on Lamb dip/crossover |
| 2 | Lambda coherence (EIT / AT / CPT) | `LambdaScheme` | transparency, AT splitting, CPT |
| 3 | Rydberg-EIT electrometry | `RydbergEITScheme` | cascade EIT, microwave AT |
| 4 | Hanle / EIA / NMOR | `MagnetoScheme` | zero-field transmission/rotation |
| 5 | FWM | `FWMScheme` | seeded gain/squeezing, SFWM biphoton |

따라서 검토 대상은 `SASScheme`, UI 제목 `Absorption spectroscopy (OD / SAS)`이다 (`gabes/schemes/sas.py:53-63`).

## 조사 범위와 기존 제안

다음 자료를 먼저 검색한 뒤 현재 코드를 검토했다.

- 프로젝트 지침과 개요: `CLAUDE.md`, `README.md`
- 이전 Scheme 1 보고서: `docs/daily_report/2026-06-25_scheme-1_od-sas.md`, `docs/daily_report/2026-07-10_scheme-1_od-sas.md`
- 계획/TODO: `docs/checklist.json`, `gabes/constants.py:79-89`
- 핵심 구현: `gabes/schemes/sas.py`, `gabes/species.py`, `gabes/core.py`, `gabes/doppler.py`, `gabes/lineshape.py`
- 물리/회귀/UI 테스트: `tests/test_sas.py`, `tests/test_absorption.py`, `tests/test_headless_observables.py`, `tests/test_schemes_render.py`
- 예제: 별도 `examples/` 디렉터리는 없다. 사용자 가이드의 실제 계산 OD/SAS 그림과 기본 app 실행이 현재 예제 역할을 한다 (`docs/GABES_User_Guide_v2.html:518-573`).
- 로컬 파일명/본문에서 issue note를 찾았지만 별도 issue 문서는 없었다. GitHub `Shake2313/fwm-squeezing-app`의 open/closed Issues도 2026-07-20에 확인했으며 양쪽 모두 0건이었다.

기존 개선 제안은 **있다**. 체크리스트에는 (1) gas/species/line별 buffer-gas 계수, pressure shift, phenomenological Dicke narrowing (`docs/checklist.json:22-26`), (2) low-order polarization/Zeeman proxy (`docs/checklist.json:88-92`), (3) full velocity-changing collision(VCC) (`docs/checklist.json:102-106`)이 남아 있다. 반면 figureless observables는 완료 상태다 (`docs/checklist.json:29-33`). 2026-07-10 보고서가 제안했던 lock readout도 현재 구현되었고, 최근에는 검출된 sub-Doppler feature의 한 FWHM 안에서 slope를 찾도록 개선되었다 (`gabes/schemes/sas.py:471-514`, `tests/test_sas.py:185-230`).

## 결론 요약

**OD/SAS는 실제 물리를 담은 유용한 scheme이며, pure/warm alkali cell의 선 배정·OD 규모·펌프 의존성·crossover 경향을 보는 semi-quantitative 실험 레퍼런스로 적합하다.** Pump-off OD와 pump-on SAS가 같은 atomic manifold와 cell parameter 위에서 연속적으로 이어지고, 실제 Rb/Cs hyperfine 자료, Doppler 평균, CG-branched decay, transit relaxation, hyperfine optical pumping을 사용한다.

그러나 지금 상태를 publication-grade linewidth/lock-point 또는 buffer-cell metrology 레퍼런스로 쓰면 안 된다. 오늘 확인한 가장 중요한 이유는 두 가지다.

1. `gamma_eff`의 Ne/self broadening이 probe Lorentzian에는 들어가지만 species-mode pump OBE의 optical dephasing에는 들어가지 않는다.
2. 기본 global scan에서 hero metric인 가장 좁은 sub-Doppler FWHM과 lock slope가 충분히 수렴하지 않는다.

둘 다 solver 차원을 키우지 않고 개선할 수 있다. 첫 항목은 같은 Liouvillian 크기에 optical-coherence dephasing을 더하는 일이고, 둘째는 해상도 경고와 half-height edge interpolation만으로도 우선 정직하게 다룰 수 있다.

## 현재 구현이 담는 실제 물리

### 강점

- `pump_power_mw = 0`이면 선형 Doppler OD, pump가 켜지면 같은 스펙트럼에 velocity-selective hole burning과 optical pumping이 생긴다 (`gabes/schemes/sas.py:67-95`, `gabes/schemes/sas.py:151-179`). OD와 SAS를 서로 무관한 toy curve로 붙인 구조가 아니다.
- Rb-85, Rb-87, Cs-133의 D1/D2 centroid, hyperfine A/B, 질량, 자연선폭과 자연동위원소비가 명시돼 있다 (`gabes/species.py:127-182`).
- `build_manifold()`는 허용된 모든 `Fg→Fe` 전이, 상대 Rabi, CG-branched spontaneous decay, thermal ground distribution으로의 transit reset을 만든다 (`gabes/species.py:337-422`). 이 decay redistribution이 실제 alkali SAS에서 중요한 enhanced/inverted crossover를 만든다.
- pump power와 1/e² waist를 Rabi frequency로 변환하고, pump OBE의 steady-state population을 velocity/detuning 축에서 푼다 (`gabes/species.py:273-283`, `gabes/schemes/sas.py:205-246`). Probe는 약한 선형 readout으로 취급된다.
- pump-off Rb-85 D1의 적분 흡수와 peak가 AutoOD validation primitive와 1% 이내이며, ground-manifold 적분 세기비 `49/25`도 고정돼 있다 (`tests/test_sas.py:101-121`, `tests/test_absorption.py:81-100`). 이는 단순한 모양 맞추기보다 훨씬 강한 기준점이다.
- pump-on sharp feature, 총흡수 감소, transit rate가 낮아질수록 강해지는 hyperfine-pumping crossover, natural-Rb isotope overlay가 테스트된다 (`tests/test_sas.py:132-173`).
- 사용자 가이드의 pump-off/pump-on 예제 그림은 같은 natural-Rb D1 조건에서 넓은 Doppler envelope와 날카로운 SAS 구조가 함께 나타남을 보여 준다 (`docs/GABES_User_Guide_v2.html:551-572`).

### 레퍼런스로서의 범위

신뢰도가 높은 사용은 다음과 같다.

- isotope와 D1/D2 hyperfine group/line의 대략적 배정
- 온도·셀 길이에 따른 OD 규모와 Doppler envelope 경향
- pump power·waist·transit relaxation에 따른 saturation/optical-pumping 방향성
- pure-cell 조건에서 lock 후보와 crossover를 찾는 1차 설계

주의가 필요한 사용은 다음과 같다.

- buffer-gas cell의 절대 line centre, pressure-shifted lock point, Dicke/VCC linewidth
- 편광과 Zeeman sublevel 분포가 중요한 feature contrast
- 실제 modulation spectroscopy error signal의 절대 discriminator slope
- 서로 다른 species/line에서 pump mW를 직접 비교하는 정량 saturation calibration

현재 lock metric은 직접 계산한 transmission의 `|dT/dΔ|` 최대값이지, EOM/FM 변조·demodulation phase·검출기 transfer function을 포함한 실험 error signal은 아니다 (`gabes/schemes/sas.py:471-514`). 따라서 “잠글 만한 국소 slope proxy”로는 유용하지만 servo gain 기준값으로 해석하면 안 된다.

## 오늘 확인한 물리·수치 정확도 위험

### 1. Species-mode에서 homogeneous broadening이 pump OBE에 빠짐

Species path는 `gamma_eff = gamma_nat + self broadening + Ne broadening`을 계산한다 (`gabes/schemes/sas.py:163-179`). 하지만 pump Liouvillian은 자연선폭으로 만들어진 `man.atom`을 그대로 사용한다 (`gabes/schemes/sas.py:215-221`, `gabes/species.py:349-351`, `gabes/species.py:394-399`). `gamma_eff`는 이후 Δ_eff 표 간격/범위와 probe Lorentzian HWHM에만 쓰인다 (`gabes/schemes/sas.py:219-246`). 반대로 generic path는 `atoms.sas_atom(..., gamma=gamma_eff)`를 만들어 pump와 probe 양쪽에 넓힘을 넣는다 (`gabes/schemes/sas.py:258-288`).

즉 species-mode에서 pressure를 바꿔도 pump population Liouvillian 자체는 동일하다. 이는 collisional homogeneous broadening이 pump saturation과 hole-burning response에도 영향을 주어야 하는 실험 상황과 맞지 않는다. 오늘 진단용으로 기존 added FWHM의 절반을 모든 optical coherence의 pure-dephasing rate로 넣어 동일한 6-level/36-density-variable solve를 돌렸을 때, Rb-85 D2, 1.5 mW, 30 °C에서 현재 결과 대비 `alpha`의 최대 상대 변화가 Ne 5 Torr에서 0.279, 20 Torr에서 0.388이었다. 이 숫자는 보정 모델의 정답이 아니라 **현재 누락에 대한 민감도 경고**다. 실행시간은 각각 약 0.29→0.30 s, 0.30→0.30 s로 사실상 같았다.

개선은 체크리스트의 buffer-gas 항목에 포함해 처리하는 것이 좋다. Spontaneous population decay는 자연 Γ로 유지하고, added homogeneous FWHM의 절반을 optical-coherence dephasing에 넣으면 solver dimension과 solve count가 변하지 않는다. Pressure shift와 gas/species/line coefficient table 역시 scalar lookup/offset이므로 런타임 증가는 무시할 수 있다. 정확한 coefficient와 dephasing convention은 문헌/셀 데이터로 검증해야 한다.

### 2. 기본 sub-Doppler FWHM/lock slope가 해상도 제한을 받음

`narrowest_subdoppler()`는 running-median residual의 half-height를 넘는 첫 sample 두 개로 폭을 반환하며 edge interpolation을 하지 않는다 (`gabes/lineshape.py:107-130`). 이 값은 현재 SAS의 hero metric이고, lock slope 탐색 범위도 그 FWHM에 의존한다 (`gabes/schemes/sas.py:354-385`, `gabes/schemes/sas.py:471-500`).

오늘 동일 조건에서 `scan_points`만 바꾼 결과는 다음과 같다.

| 조건 | points | sample 간격 | FWHM | lock slope | lock detuning |
|---|---:|---:|---:|---:|---:|
| natural Rb D1 기본 | 401 | 28.200 MHz | 56.4 MHz | 0.0067/MHz | -1143.8 MHz |
| 〃 | 1401 (기본) | 8.057 MHz | 32.2 MHz | 0.0241/MHz | -1304.9 MHz |
| 〃 | 4001 | 2.820 MHz | 22.6 MHz | 0.0316/MHz | -1301.7 MHz |
| Rb-85 D2, 1.5 mW | 401 | 16.848 MHz | 33.7 MHz | 0.0074/MHz | -1147.0 MHz |
| 〃 | 1401 | 4.814 MHz | 14.4 MHz | 0.0170/MHz | -1171.0 MHz |
| 〃 | 4001 | 1.685 MHz | 15.2 MHz | 0.0222/MHz | -1160.4 MHz |

Natural-Rb 기본 FWHM은 1401→4001점에서 30% 변하고, 기본 1401점 폭은 네 sample 정도뿐이다. Rb-85 D2 FWHM은 더 안정적이지만 lock slope는 26% 변한다. 현재 테스트는 검출 feature 근처에 lock point가 있는지는 확인하지만 scan-resolution convergence는 고정하지 않는다 (`tests/test_sas.py:185-230`).

가장 싼 개선은 (1) `FWHM / sample_step`을 함께 계산해 sample 수가 부족하면 `resolution-limited` 상태를 내고, (2) half-height 양쪽 edge를 선형 보간하며, (3) 1401/4001점 대표 회귀 테스트를 추가하는 것이다. 이들은 추가 OBE solve 없이 가능하다. Publication-grade 숫자가 필요할 때만 선택적 local refinement를 켜는 편이 좋다.

### 3. 모든 species/line에 고정된 Rb-85 D1 `I_sat`

`I_SAT = 4.484 mW/cm²`는 Rb-85 D1 상수 묶음에 한 번 정의돼 있고 (`gabes/constants.py:22-31`), `pump_rabi_from_power()`는 모든 isotope와 D1/D2에 같은 `I_sat`을 쓴다 (`gabes/species.py:273-283`). Line별 자연 Γ는 바뀌지만 `Ω/Γ = sqrt(I/(2 I_sat))`는 동일하다. 실제로 0.5 mW, 1 mm 조건에서 코드의 `Ω/Γ`는 Rb-85/Rb-87/Cs 및 D1/D2 모두 `1.883984`였다.

따라서 pump-power trend는 유용하지만 species/line 사이의 절대 saturation power 비교는 아직 정량 기준이 아니다. 명시적인 polarization convention과 함께 per-species/line saturation-intensity 또는 reduced-dipole calibration table을 넣으면 lookup 비용만 들고, full Zeeman solve 없이도 reference 가치가 올라간다. 다만 검증 자료 없이 새 fit knob만 늘려서는 안 된다.

## 기존 개선안의 계산비용 평가

| 기존 항목 | 물리 가치 | 계산비용 | 판단 |
|---|---|---|---|
| Gas/species/line coefficient + pressure shift | buffer cell의 centre/폭 해석 | table lookup와 scalar offset | **거의 0, 우선 적용 가능** |
| Phenomenological Dicke narrowing | 현재의 단순 broadening 편향 완화 | 기존 width algebra 유지 시 거의 0 | **저비용**, 적용 범위 명시 필요 |
| Added broadening을 pump optical dephasing에도 반영 | saturation/hole-burning 일관성 | 동일 Liouvillian 크기·동일 solve 수 | **저비용이면서 correctness 우선순위 높음** |
| Low-order polarization/Zeeman proxy | 편광 불순도·비대칭·contrast 설명 | transition weight/rate scalar 보정이면 낮음 | 저비용이나 calibration 없이는 fit knob 위험 |
| Full VCC | buffer-rich cell의 velocity redistribution | velocity classes 결합, block/iterative solve 가능성 | **고비용**, 별도 설계와 사용자 범위 결정 필요 |
| Figureless observables | batch/report 속도 | 이미 완료 | 유지; `tests/test_headless_observables.py:38-89`가 보호 |

## 순수 코드 최적화 후보

### 1. Pump-off analytic population fast path — 가장 확실함

현재 `pump_power_mw = 0`에서도 pump Hamiltonian/Liouvillian, `_pump_pops()` 표, level별 `np.interp`를 모두 실행한다 (`gabes/schemes/sas.py:215-246`, `gabes/schemes/sas.py:461-468`). 하지만 이 극한에서는 ground population이 `p_ground`, excited population이 0으로 detuning과 무관하므로 `w=1`을 바로 쓸 수 있다.

오늘 Rb-85 D1, 90 °C, 12.5 mm에서 이 경로를 진단용으로 대체했을 때:

- 현재: 0.2573 s
- analytic pump-off path: 0.0823 s
- speedup: 3.13×
- 최대 상대 `alpha` 차이: `5.63×10⁻16`

물리와 동작을 바꾸지 않는 가장 좋은 최적화다. Pump-off AutoOD 회귀와 natural-Rb isotope overlay 테스트로 보호하면 된다.

### 2. Running median vectorization — headless readout 개선

`narrowest_subdoppler()`는 1401개 위치마다 Python list comprehension으로 median을 다시 계산한다 (`gabes/lineshape.py:113-119`). 같은 edge padding과 window를 `numpy.lib.stride_tricks.sliding_window_view` + axis median으로 바꾼 진단에서:

- 현재 median 단계: 21.507 ms
- vectorized 단계: 0.218 ms
- speedup: 98.6×
- 결과 최대 절대 차이: 0

이는 전체 OBE compute가 아니라 observables의 해당 단계만의 수치다. 그래도 현재 headless readout이 약 0.04–0.06 s이므로 체감 가능한 순수 최적화다.

### 3. Pump-on interpolation/accumulation fused kernel — 그 다음 후보

Natural-Rb D1 기본 compute의 cProfile 결과 0.753 s 중 `_component_alpha()`가 0.736 s, 8회의 `np.interp`가 0.230 s, `_pump_pops()`가 0.133 s였다 (`gabes/schemes/sas.py:205-246`). Level별 `pop_at` 대형 배열을 모두 만든 뒤 transition별 Lorentzian을 누적하는 구조가 핵심 병목이다.

단순히 uniform-grid index arithmetic을 NumPy로 바꾼 실험은 오히려 0.57×로 느려졌으므로 권장하지 않는다. 가치가 있는 방향은 기존 optional Numba kernel 패턴을 따라 interpolation과 transition/velocity 누적을 한 kernel에서 처리하고 중간 배열과 메모리 traffic을 줄이는 것이다. 먼저 bitwise/엄격 tolerance 회귀와 401/1401/4001점 benchmark가 필요하다.

이전 보고서의 manifold skeleton/marker 캐시는 우선순위가 낮다. `build_manifold()`는 이미 `lru_cache(maxsize=64)`를 쓴다 (`gabes/species.py:337-338`), 오늘 profile에서 line-strength/Wigner/marker 관련 비용은 수 ms 이하였다.

## 검증 결과

관련 테스트:

```text
python -m pytest tests/test_sas.py tests/test_absorption.py tests/test_headless_observables.py -q
31 passed in 21.13s
```

전체 저장소 테스트:

```text
MPLBACKEND=Agg python -m pytest -q
162 passed in 32.16s
```

현재 테스트는 atomic constants/strengths, AutoOD scale, pressure-broadened cold OD, pump-on feature, hyperfine-pumping crossover, figure/headless contract, lock target 위치를 잘 보호한다. 추가로 필요한 회귀는 species-mode pump OBE의 added-dephasing 반영과 FWHM/lock-slope 해상도 수렴이다.

## 최종 우선순위

1. **즉시·무부하:** sub-Doppler FWHM에 samples-per-width와 `resolution-limited` 상태를 붙이고 half-height edge를 보간한다.
2. **저부하 correctness:** 기존 buffer-gas 개선안에 species-mode pump optical dephasing 반영을 포함하고 pressure shift/coefficient table을 함께 검증한다.
3. **저부하 calibration:** species/line/polarization convention별 saturation-intensity 또는 dipole calibration을 도입한다.
4. **순수 성능:** pump-off analytic fast path, running-median vectorization, 이후에만 fused pump-on kernel을 검토한다.
5. **필요할 때만 고비용:** full VCC와 full Zeeman/polarization solver는 실제 대상 셀·편광·정확도 요구가 정해진 뒤 범위를 합의한다.

종합하면 OD/SAS는 이미 “실제 실험에 쓸 만한 물리”를 구현한다. 특히 pure-cell hyperfine spectroscopy의 구조와 pump-off 절대 OD 기준은 강하다. 다음 단계는 무거운 full model보다, 현재 출력의 수치 해상도를 정직하게 표시하고 이미 계산한 homogeneous broadening을 pump OBE에도 일관되게 전달하는 것이다. 두 개선 모두 interactive 성능을 사실상 보존하면서 실험 레퍼런스 신뢰도를 크게 높인다.
