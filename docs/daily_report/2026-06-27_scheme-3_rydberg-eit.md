# 2026-06-27 Scheme 3 Review: Rydberg-EIT electrometry

## 오늘의 선택 스킴

- 오늘의 현지 날짜는 `2026-06-27`이며 규칙 `n = (day mod 5) + 1`에 따라 `n = (27 mod 5) + 1 = 3`이다.
- 현재 등록된 다섯 개 스킴 순서는 `OD / SAS -> Lambda coherence -> Rydberg-EIT electrometry -> Hanle / EIA / NMOR -> FWM`이다 (`gabes/schemes/__init__.py:19-24`, `README.md:8-16`).
- 따라서 오늘 검토 대상은 3번째 스킴 `RydbergEITScheme`가 담당하는 `Rydberg-EIT electrometry`이다 (`gabes/schemes/rydberg.py:80-88`).

## 이번 검토에서 읽은 근거

- 스킴 등록/순서: `gabes/schemes/__init__.py:19-24`
- 스킴 본체: `gabes/schemes/rydberg.py:80-577`
- residual Doppler를 위한 per-level `k` ratio 정의: `gabes/atoms.py:21-48`, `gabes/atoms.py:84-98`
- 증기 밀도/Beer-Lambert 연결: `gabes/species.py:192-207`, `gabes/observables.py`의 `absorption_coefficient()` / `transmission()` 호출 위치 `gabes/schemes/rydberg.py:385-395`
- 기준 테스트: `tests/test_rydberg_eit.py:29-200`, `tests/test_kernels.py:116-136`
- 기존 개선안 메모: `docs/checklist.json:4-10`
- 이전 일일 보고서: `docs/daily_report/2026-06-22_scheme-3_rydberg-eit.md`

## 현재 정의와 물리 모델

이 스킴은 85Rb의 정적 4준위 cascade ladder `5S1/2 F=3 -> 5P3/2 F'=4 -> 40D5/2 -> 39F7/2`를 푼다. RF leg는 `40D5/2 -> 39F7/2`이며, UI는 `EIT`와 `AT electrometry` 두 regime를 한 스킴 안에서 바꿔 보여 준다 (`gabes/schemes/rydberg.py:114-118`, `gabes/schemes/rydberg.py:191-211`).

중요한 점은 지금 HEAD에서 beam/power knob가 더 이상 단순 display-only가 아니라는 것이다.

- `probe_power_uw`는 `_probe_rabi()`를 통해 약한 probe Rabi로 들어간다 (`gabes/schemes/rydberg.py:119-123`, `gabes/schemes/rydberg.py:226-234`, `gabes/schemes/rydberg.py:316`).
- `coupling_power_mw`와 `beam_diameter_mm`는 `_coupling_rabi()`로부터 유효 `Omega_c`를 정하고 (`gabes/schemes/rydberg.py:124-142`, `gabes/schemes/rydberg.py:213-224`, `gabes/schemes/rydberg.py:317`), 같은 `beam_diameter_mm`가 transit broadening도 정한다 (`gabes/schemes/rydberg.py:129-137`, `gabes/schemes/rydberg.py:236-243`, `gabes/schemes/rydberg.py:309`).
- dephasing은 이제 `5S-40D`와 `40D-39F` 두 채널로 분리되어 있다 (`gabes/schemes/rydberg.py:151-165`, `gabes/schemes/rydberg.py:191-210`, `gabes/schemes/rydberg.py:307-314`).
- `doppler="on"`일 때는 residual `(k_probe - k_coupling) v` mismatch를 Maxwell average로 태우고, 이를 위해 `AtomModel.doppler_ratios`가 실제로 사용된다 (`gabes/atoms.py:31-36`, `gabes/atoms.py:84-98`, `gabes/schemes/rydberg.py:168-175`, `gabes/schemes/rydberg.py:333-349`).
- readout은 transmission 기반이며 EIT에서는 linewidth, AT에서는 splitting과 center shift를 준다 (`gabes/schemes/rydberg.py:423-459`).

즉 현재 모델은 “정적 optical spectrum + microwave dressing”에 초점을 맞춘 4준위 lumped ladder이고, time-domain superheterodyne demodulation은 여전히 범위 밖이다 (`gabes/schemes/rydberg.py:8-11`, `gabes/schemes/rydberg.py:560-576`).

## 직접 실행과 수치 확인

실행:

```bash
python -m pytest tests/test_rydberg_eit.py tests/test_kernels.py -q
```

결과:

- `26 passed in 12.14s`
- 따라서 현재 HEAD에서는 reference-default test, probe-power broadening, AT center shift, residual Doppler broadening, affine kernel fast path 일치성이 모두 통과한다 (`tests/test_rydberg_eit.py:43-127`, `tests/test_rydberg_eit.py:150-200`, `tests/test_kernels.py:123-136`).

추가로 직접 수치를 뽑아 보면:

- 기본 `EIT`에서 `EIT linewidth ~= 1.61 MHz`, `Transmission at resonance ~= 0.940`였다.
- 기본 `AT electrometry`에서 `RF AT splitting ~= 3.50 MHz`로 기본 `lo_rabi_mhz = 3.7`를 약간 밑도는 정도였고, 이는 test가 요구하는 `0.88 <= split / lo <= 1.0` 범위와 맞는다 (`tests/test_rydberg_eit.py:52-61`).
- 워밍업 후 평균 시간은 `EIT compute ~= 9.23 ms`, `AT compute ~= 5.01 ms`, `_readout ~= 1 ms`, `observables() ~= 173-176 ms`였다. 즉 지금도 solver보다 Matplotlib figure 생성이 체감 병목이다 (`gabes/schemes/rydberg.py:423-513`).
- `doppler=off` 대비 `doppler=on` 평균 compute 시간은 약 `11.8 ms -> 518.8 ms`였다. residual Doppler는 물리적으로 의미가 있지만 “거의 공짜”는 아니다.

파라미터 반응도 물리적으로 꽤 납득 가능했다.

- `beam_diameter_mm: 0.30 mm`로 넓히면 `Omega_c`와 probe drive가 함께 줄고 transit broadening이 `0.076 MHz`로 감소해 `EIT linewidth ~= 0.68 MHz`까지 좁아졌다.
- 반대로 `beam_diameter_mm: 0.10 mm`로 줄이면 `Omega_c`, probe drive, transit broadening이 모두 커져 `EIT linewidth ~= 2.67 MHz`가 되었다.
- `probe_power_uw: 1 -> 10 uW`는 `probe_rabi_mhz: 0.816 -> 2.582`로 바뀌며 linewidth가 `1.42 -> 1.75 MHz`로 넓어졌다. Ju 논문의 probe-power broadening 방향성과 맞다 (`gabes/schemes/rydberg.py:226-234`, `tests/test_rydberg_eit.py:174-184`).
- `mw_detuning_mhz: +/-4`에서는 `AT center shift ~= -/+1.88 MHz`가 나와 detuned microwave가 dressed-doublet 중심을 끄는 효과가 metric으로 읽힌다 (`gabes/schemes/rydberg.py:433-450`, `tests/test_rydberg_eit.py:93-115`).
- `temp_c: 20 -> 60 C`에서는 EIT linewidth는 거의 그대로인데 공진 transmission이 크게 떨어졌다. 현재 temperature가 주로 number density와 transit term을 통해 들어가고, collision-enhanced dephasing이나 pressure shift는 직접 모델링하지 않기 때문이다 (`gabes/species.py:192-207`, `gabes/schemes/rydberg.py:236-243`, `gabes/schemes/rydberg.py:301-315`).

## 실험물리학 연구자 관점 평가

### 어디까지는 실험에 유용한가

이 스킴은 **정적 Rydberg-EIT/AT 스펙트럼을 빠르게 탐색하는 실험 planning reference**로는 상당히 유용하다.

- 실험자가 실제로 만지는 probe power, coupling power, beam diameter가 이제 스펙트럼에 직접 연결된다. 이 변화는 lab-facing 가치가 크다 (`gabes/schemes/rydberg.py:119-142`, `gabes/schemes/rydberg.py:213-243`, `gabes/schemes/rydberg.py:316-318`).
- beam diameter가 “밝기 상승”과 “transit broadening 증가”를 동시에 일으키도록 묶어 둔 점이 좋다. 좁은 빔이 무조건 좋지 않다는 warm-vapor 실험 감각을 그대로 전달한다 (`gabes/schemes/rydberg.py:129-137`, `gabes/schemes/rydberg.py:236-243`).
- compensated / uncompensated EIT를 함께 보여 주고, probe-power sweep extra view까지 둔 것은 실제 정렬과 자장 보상 작업에 매우 실용적이다 (`gabes/schemes/rydberg.py:467-483`, `gabes/schemes/rydberg.py:515-558`).
- affine-scan kernel 검증이 별도 테스트로 묶여 있어, 빠른 경로를 쓴다고 물리 결과가 틀어질 가능성이 낮다 (`tests/test_kernels.py:123-136`).

그래서 오늘 기준 이 코드는 **정성적을 넘어서 반정량적 semi-quantitative Rydberg-EIT reference**로는 충분히 참고할 만하다. 적어도 “어떤 knob가 linewidth, contrast, splitting을 어느 방향으로 움직이는가”를 읽는 데에는 꽤 강하다.

### 아직 최종 실험 레퍼런스로 보기 어려운 지점

다만 **절대적인 electrometry calibration reference**라고 부르기에는 아직 선이 있다.

- 모델은 여전히 Zeeman sublevel, polarization selection rule, stray E/B mixing, state-selective ionization, optical pumping redistribution을 해상하지 않는 lumped 4준위다 (`gabes/schemes/rydberg.py:191-211`).
- `if_khz`는 아직 solve나 observables에 들어가지 않으므로, 실제 superheterodyne readout sensitivity는 계산하지 못한다. 현재 `Max spectral slope`는 좋은 proxy이지만 detector chain 자체는 아니다 (`gabes/schemes/rydberg.py:166-167`, `gabes/schemes/rydberg.py:454-458`, `gabes/schemes/rydberg.py:560-576`).
- temperature는 density와 transit term에는 들어가지만 collision-enhanced dephasing, pressure shift, Dicke narrowing, charge/field-noise broadening으로 연결되지는 않는다 (`gabes/species.py:192-207`, `gabes/schemes/rydberg.py:236-243`, `gabes/schemes/rydberg.py:301-315`).
- `doppler=on`이 residual mismatch를 보여 주긴 하지만, 현재 help 문구가 직접 말하듯 calibrated reference가 아니라 “what-if” 용이다 (`gabes/schemes/rydberg.py:168-175`).

정리하면 이 스킴은 **정적 스펙트럼 실험 설계용으로는 꽤 좋은 반정량 레퍼런스**이지만, **최종 감도 예측이나 절대적인 전계 캘리브레이션 레퍼런스**로 쓰기에는 아직 물리 채널이 모자란다.

## 기존 문서에서 찾은 개선안과 현재 반영 상태

이번 검토에서는 “기존에 제안된 개선안이 있는가”에 대해 답을 분명히 할 수 있었다. 있다. 다만 상당수는 이미 반영되었다.

### 1. checklist의 `rydberg-power-to-rabi`

- `docs/checklist.json:4-10`의 `rydberg-power-to-rabi` 항목이 대표적인 기존 제안이다.
- 그 내용은 coupling power/beam geometry를 `Omega_c`에 연결하고, residual Doppler, 2-channel dephasing, AT center shift, cached constants까지 함께 넣자는 것이다.
- 현재 HEAD는 이 항목을 사실상 이미 반영했다. 대응 위치는 `_coupling_rabi()` (`gabes/schemes/rydberg.py:213-224`), residual Doppler option (`gabes/schemes/rydberg.py:168-175`, `gabes/schemes/rydberg.py:333-349`), 2-channel dephasing (`gabes/schemes/rydberg.py:151-165`, `gabes/schemes/rydberg.py:191-210`), `AT center shift` metric (`gabes/schemes/rydberg.py:433-450`), cached probe-line constants (`gabes/schemes/rydberg.py:45-78`)이다.
- 계산 부하는 대부분 매우 작다. `_coupling_rabi()`, `_probe_rabi()`, `_transit_rate_mhz()`는 전부 스칼라 전처리라 solver의 상태 수나 scan 길이를 바꾸지 않는다 (`gabes/schemes/rydberg.py:213-243`).
- 실제 측정에서도 warm `compute()`는 여전히 `5-12 ms` 수준이다. 즉 이 개선안은 **연산량 자체는 거의 늘리지 않고 물리 유용성은 크게 올린 사례**로 보인다.
- 다만 UI 차원에서는 `probe_power_uw`, `coupling_power_mw`, `beam_diameter_mm`가 이제 solve knob가 되었으므로, 예전보다 “슬라이더를 움직일 때 재계산이 더 자주 일어나는” 비용은 생겼다 (`gabes/schemes/rydberg.py:119-137`).

### 2. 이전 일일 보고서의 제안

`docs/daily_report/2026-06-22_scheme-3_rydberg-eit.md`에는 power-to-Rabi, residual Doppler, 두 채널 dephasing, `AT center shift`, metric-only fast path 같은 제안이 남아 있다. 그런데 현재 코드를 보면 그중 상당수가 이미 구현되었다.

- `_readout()` headless path가 이미 들어와 있어 metric/table 계산과 figure build를 분리한다 (`gabes/schemes/rydberg.py:423-459`).
- `AT center shift`도 이미 metric으로 들어가 있다 (`gabes/schemes/rydberg.py:441-450`).
- residual Doppler와 two-channel dephasing도 위에서 본 대로 구현되었다.

오히려 현재는 **기존 문서가 최신 코드를 못 따라온 부분**이 보인다.

- 예전 보고서에는 `probe_power_uw`, `coupling_power_mw`, `beam_diameter_mm`가 display-only라고 적혀 있으나 (`docs/daily_report/2026-06-22_scheme-3_rydberg-eit.md`), 현재 코드는 그렇지 않다 (`gabes/schemes/rydberg.py:119-142`, `gabes/schemes/rydberg.py:213-243`, `gabes/schemes/rydberg.py:316-318`).
- `docs/checklist.json:8`의 설명도 `probe_power_uw stays display-only`라고 되어 있는데, 현재 구현은 probe power까지 OBE drive에 연결했다 (`gabes/schemes/rydberg.py:226-234`, `tests/test_rydberg_eit.py:174-184`).

따라서 “기존 개선안을 반영하면 부하가 얼마나 드는가”라는 질문에 대한 오늘의 답은:

- **power-to-Rabi / center-shift / 2채널 dephasing / cached constants는 거의 공짜에 가깝다.**
- **residual Doppler는 물리 가치는 있지만 명백히 무겁다.** `doppler=off ~ 11.8 ms`, `doppler=on ~ 518.8 ms`.

## 남아 있는 조심스러운 후보 개선안

기존 개선안이 아예 없는 것은 아니었고 상당수는 이미 반영되었으므로, 오늘 새 제안은 “남은 빈 곳”에만 조심스럽게 한정한다.

### 1. IF / superheterodyne readout proxy를 한 단계 더 실험가 친화적으로 확장

- 지금은 `if_khz`가 메타데이터이고, readout은 정적 transmission slope까지만 본다 (`gabes/schemes/rydberg.py:166-167`, `gabes/schemes/rydberg.py:454-458`, `gabes/schemes/rydberg.py:560-576`).
- 실제 vapor-cell electrometry에서는 slope 자체보다 “IF 체인에서 어느 detuning이 가장 큰 discriminator를 주는가”가 중요할 때가 많다.
- full time-domain demodulation까지 가지 않더라도, lock-in/superhet proxy를 추가하는 정도는 계산 부하를 거의 늘리지 않고 실험적 유용성을 올릴 수 있다. 예를 들면 small-signal FM dithering에 대한 `dT/dnu`, `d²T/dnu²`, 또는 finite `if_khz`에서의 차분 readout proxy를 metric으로 둘 수 있다.

### 2. temperature-linked phenomenological dephasing / pressure-shift 항

- 현재 temperature는 density와 transit에 비해 coherence broadening 채널에 약하게만 연결된다 (`gabes/species.py:192-207`, `gabes/schemes/rydberg.py:236-243`, `gabes/schemes/rydberg.py:301-315`).
- warm vapor 실험에서는 셀 온도 상승이 충돌성 broadening, charge noise, stray-field sensitivity를 함께 키우는 경우가 많다.
- 따라서 `rydberg_dephasing_mhz = gamma0 + a * (T - T_ref)` 같은 낮은 차수의 phenomenological 항을 추가하면, 계산 부하를 사실상 거의 늘리지 않으면서 실험 감각을 더 잘 따라갈 수 있다.

### 3. polarization / Zeeman 효과의 “저차 버전”

- full Zeeman manifold까지 가면 모델과 계산이 모두 무거워진다.
- 그러나 완전한 다준위 확장 대신, effective polarization selector나 uncompensated broadening asymmetry 같은 1-2개 스칼라 knob만 추가해도 stray polarization mismatch나 optical-pumping sensitivity를 조금 더 잘 흉내낼 수 있다.
- 이 경로는 residual Doppler full average보다 계산 부하가 훨씬 작고, 실험가가 느끼는 “왜 오늘 피크가 이렇게 찌그러졌지?”를 설명하는 데 더 직접적일 수 있다.

## 물리와 무관한 순수 코딩 최적화 관점

### 1. `observables()` 병목은 여전히 figure 생성이다

- `_readout()`는 대략 `1 ms`인데 `observables()`는 `173-176 ms`였다 (`gabes/schemes/rydberg.py:423-459`, `gabes/schemes/rydberg.py:461-513`).
- 현재도 `Scheme.cache_observables = True`라 unchanged rerun은 캐시되지만 (`gabes/schemes/base.py:73-77`, `streamlit_app.py:756-761`), solve knob가 바뀌면 figure는 다시 만든다.
- 따라서 자동 보고서나 batch sweep 경로에서 figure가 필요 없을 때 `_readout()`만 쓰는 public/headless 진입점을 더 노골적으로 노출하면 속도 이득이 있다. 물리는 전혀 바뀌지 않는다.

### 2. `_atom()`의 Liouvillian 재조립을 더 줄일 여지

- 지금도 topology skeleton은 `_cascade_skeleton()`으로 캐시해 두었지만 (`gabes/schemes/rydberg.py:58-78`), `_atom()`은 호출마다 새 `AtomModel`을 만들고 그 안에서 Lindblad와 `S_v`를 다시 구성한다 (`gabes/atoms.py:45-48`, `gabes/schemes/rydberg.py:191-211`).
- dephasing이 바뀌는 항은 사실 소수의 diagonal damping뿐이므로, base dissipator를 캐시하고 dephasing delta만 더하는 식으로 가면 compute를 조금 더 깎을 수 있다.

### 3. Doppler-on 경로의 velocity grid 캐시

- `doppler=on`의 큰 비용 중 일부는 velocity-grid 생성과 다중 velocity solve에서 온다 (`gabes/schemes/rydberg.py:337-343`).
- 동일한 `temp_c`, `dv`, `cutoff_sigma` 조합에서 `doppler.velocity_grid()` 결과를 캐시하면 특히 반복 스캔에서 약간의 이득이 있다.
- 물론 핵심 비용은 여전히 다중 class solve라서, 이 최적화는 “조금 더 나아짐” 정도이지 게임체인저는 아니다.

## 결론

- 오늘 기준 `RydbergEITScheme`은 **정적 Rydberg-EIT/AT 실험 설계와 knob sensitivity 파악용으로는 꽤 강한 semi-quantitative reference**다.
- 특히 예전 제안이던 power-to-Rabi, beam-size/transit coupling, residual Doppler option, 2채널 dephasing, `AT center shift`, headless readout이 이미 구현되어 있어서, 6월 22일 평가보다 실험가 친화성이 확실히 올라갔다.
- 반면 **최종 electrometry calibration reference**로 쓰기에는 IF/demodulation chain, 충돌성/압력성 broadening, polarization/Zeeman 분해가 아직 부족하다.
- 비용 대비 가장 좋은 다음 한 걸음은 **IF 기반 readout proxy**와 **temperature-linked phenomenological dephasing**처럼 계산 부하를 거의 늘리지 않는 실험가용 readout 보강이다.
