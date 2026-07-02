# 2026-07-02 Scheme 3 Review: Rydberg-EIT electrometry

## 선택 결과

- 오늘 현지 날짜는 `2026-07-02`이고, 규칙 `n = (day mod 5) + 1`에 따라 `n = (2 mod 5) + 1 = 3`이다.
- 현재 등록된 다섯 개 스킴 순서는 `OD / SAS -> Lambda coherence -> Rydberg-EIT electrometry -> Hanle / EIA / NMOR -> FWM`이다 (`gabes/schemes/__init__.py:19-24`, `README.md:8-16`).
- 따라서 오늘 검토 대상은 3번째 스킴인 `RydbergEITScheme`, 즉 `Rydberg-EIT electrometry`이다 (`gabes/schemes/rydberg.py:81-114`).

## 읽은 근거

- 스킴 등록/순서: `gabes/schemes/__init__.py:19-24`
- Rydberg-EIT 본체: `gabes/schemes/rydberg.py:81-616`
- 공통 headless observables 계약: `gabes/schemes/base.py:68-106`
- UI 캐시 경로: `streamlit_app.py:465`, `streamlit_app.py:757-758`
- Rydberg 전용 테스트: `tests/test_rydberg_eit.py:29-232`
- residual Doppler per-level ratio: `gabes/atoms.py:32`, `gabes/atoms.py:88-98`
- 기존 개선안 메모: `docs/checklist.json:7`, `docs/checklist.json:28`, `docs/checklist.json:63`, `docs/checklist.json:105`, `docs/checklist.json:112`
- 이전 Rydberg 일일 리뷰: `docs/daily_report/2026-06-22_scheme-3_rydberg-eit.md`, `docs/daily_report/2026-06-27_scheme-3_rydberg-eit.md`

## 현재 물리 구현 평가

현재 `RydbergEITScheme`은 85Rb `5S F=3 -> 5P F'=4 -> 40D -> 39F`의 축약 4준위 cascade OBE다. 780 nm probe, 481 nm coupling, 37 GHz RF leg를 lumped ladder로 묶고, optical susceptibility를 Beer-Lambert transmission으로 바꾼다 (`gabes/schemes/rydberg.py:1-14`, `gabes/schemes/rydberg.py:200-219`, `gabes/schemes/rydberg.py:307-394`, `gabes/schemes/rydberg.py:402-411`).

실험가 관점에서 좋은 점은 lab knob가 이미 상당히 들어왔다는 것이다.

- `probe_power_uw`는 `_probe_rabi()`를 통해 probe Rabi로 들어가며 probe-power broadening을 만든다 (`gabes/schemes/rydberg.py:119-123`, `gabes/schemes/rydberg.py:235-244`, `gabes/schemes/rydberg.py:344`).
- `coupling_power_mw`와 `beam_diameter_mm`는 `_coupling_rabi()`로 coupling Rabi를 정하고, 같은 beam diameter가 transit broadening도 정한다 (`gabes/schemes/rydberg.py:124-142`, `gabes/schemes/rydberg.py:222-233`, `gabes/schemes/rydberg.py:246-253`, `gabes/schemes/rydberg.py:345`).
- `residual_zeeman_mhz`는 compensated/uncompensated EIT curve 차이로 보이며, `temp_dephasing_mhz_per_c`는 기본값 0인 opt-in phenomenological broadening이다 (`gabes/schemes/rydberg.py:154-165`, `gabes/schemes/rydberg.py:319-326`, `tests/test_rydberg_eit.py:79-91`).
- `if_khz`는 이제 단순 metadata가 아니라 finite-difference IF discriminator metric에 들어간다 (`gabes/schemes/rydberg.py:175-176`, `gabes/schemes/rydberg.py:448-460`).
- AT regime에서는 `RF AT splitting`뿐 아니라 detuned microwave의 dressed-doublet center shift도 metric으로 나온다 (`gabes/schemes/rydberg.py:462-484`, `tests/test_rydberg_eit.py:108-130`).

직접 실행한 기본 수치도 물리 방향성은 꽤 설득력 있다.

- EIT 기본값: `EIT linewidth = 1.61 MHz`, `Transmission at resonance = 0.940`, `Max spectral slope = 0.122 /MHz`, `IF discriminator = 0.122 /MHz`.
- AT 기본값: `RF AT splitting = 3.50 MHz`, `AT center shift = +0.00 MHz`, `Transmission at resonance = 0.806`.
- `probe_power_uw: 1 -> 6 -> 10 uW`에서 linewidth가 `1.417 -> 1.613 -> 1.754 MHz`로 증가했다.
- `beam_diameter_mm: 0.10 -> 0.15 -> 0.30 mm`에서 linewidth가 `2.674 -> 1.613 -> 0.680 MHz`로 줄었다. 작은 beam이 intensity를 키우는 동시에 transit broadening도 키운다는 감각을 잘 보여 준다.

따라서 이 스킴은 **정적 Rydberg-EIT/AT spectrum planning, knob sensitivity, Ju et al. 식 기준 조건 재현 확인용 semi-quantitative reference**로는 충분히 유용하다. 다만 최종 계측 장비 calibration reference로 쓰기에는 아직 좁다. Zeeman sublevel, polarization selection, stray E/B field mixing, optical-pumping redistribution, ionization/charge-noise channel은 lumped 4준위 안에 들어 있지 않다 (`gabes/schemes/rydberg.py:200-219`, `docs/checklist.json:105`). 또한 full time-domain lock-in/superheterodyne chain은 아직 static discriminator 수준이다 (`gabes/schemes/rydberg.py:440-460`, `gabes/schemes/rydberg.py:603-616`, `docs/checklist.json:112`).

## 기존 개선안 반영 상태와 부하

기존 개선안은 명확히 존재한다. 가장 중요한 항목은 `docs/checklist.json`의 `rydberg-power-to-rabi`였고, 현재는 `done`으로 정리되어 있다 (`docs/checklist.json:7`). 내용상 power-to-Rabi, residual Doppler option, two-channel dephasing, AT center-shift, cached constants, finite-IF proxy, opt-in temperature dephasing이 함께 들어온 상태다.

부하 평가는 두 층으로 나뉜다.

- 거의 공짜에 가까운 개선: `_coupling_rabi()`, `_probe_rabi()`, `_transit_rate_mhz()`, `temp_dephasing_mhz_per_c`, `AT center shift`, `IF discriminator`는 스칼라 전처리 또는 이미 계산된 spectrum의 후처리다 (`gabes/schemes/rydberg.py:222-253`, `gabes/schemes/rydberg.py:448-484`). solver 차원이나 scan 길이를 늘리지 않으므로 물리를 따라가는 비용 대비 효과가 매우 좋다.
- 무거운 개선: `doppler="on"`은 residual two-photon Doppler를 Maxwell velocity grid로 평균하므로 비용이 크다 (`gabes/schemes/rydberg.py:177-185`, `gabes/schemes/rydberg.py:350-357`). 오늘 환경에서는 `doppler=off` EIT compute 평균이 약 `29 ms`였고, `doppler=on`은 약 `5.3 s`였다. 이 옵션은 calibrated default라기보다 what-if knob로 유지하는 현재 선택이 타당하다.

이전 보고서에서 남아 있던 `IF / superheterodyne readout proxy`와 `temperature-linked dephasing`은 현재 코드에 부분 반영되었다 (`gabes/schemes/rydberg.py:160-176`, `gabes/schemes/rydberg.py:319-326`, `gabes/schemes/rydberg.py:448-460`). 따라서 오늘의 새 제안은 더 좁게 잡는 편이 맞다.

## 조심스러운 물리 개선 후보

1. **IF proxy를 noise-aware sensitivity proxy로 확장**

현재 IF metric은 `T(nu +/- IF)`의 finite difference로 static discriminator를 뽑는다 (`gabes/schemes/rydberg.py:448-460`). 여기에 probe power 기반 shot-noise scale, detector bandwidth, optional technical noise floor를 곱해 `nV/cm/sqrt(Hz)`가 아니라도 normalized sensitivity proxy를 만들 수 있다. spectrum을 새로 풀 필요가 없고 후처리만 늘어나므로 계산 부하는 거의 없다. 다만 공개 숫자가 Ju et al. reference sensitivity처럼 보이지 않도록 "relative / proxy"임을 UI와 테스트에서 분명히 해야 한다 (`tests/test_rydberg_eit.py:64-77`).

2. **low-order polarization / Zeeman proxy**

`docs/checklist.json:63`에는 OD/SAS와 Rydberg에 effective polarization selector 또는 uncompensated asymmetry broadening 같은 저차 proxy를 넣자는 deferred 항목이 있다. full manifold는 `docs/checklist.json:105`의 GROUP C에 해당해 무겁지만, asymmetry/contrast-reduction scalar는 spectrum 후처리나 dephasing budget에만 걸 수 있다. 계산 부하는 작고, 실제 실험자가 흔히 보는 편광 불순도와 미보상 자기장 왜곡을 설명하는 데 도움이 된다.

3. **temperature dephasing의 기본 경험식은 아직 보류**

`temp_dephasing_mhz_per_c`는 이미 opt-in으로 있다 (`gabes/schemes/rydberg.py:160-165`). 이를 자동 기본값으로 켜면 실험 친화성은 올라가지만, vapor cell 조건과 stray-field 환경에 따라 계수가 크게 달라진다. 지금처럼 기본 0으로 두고, 문서/프리셋에 "warm-cell stress test"만 추가하는 쪽이 물리적으로 더 정직해 보인다. 계산 부하는 없지만 calibration 책임이 생긴다.

## 순수 코딩 최적화 후보

- 오늘 측정에서 `_readout()`은 평균 약 `2.2 ms`, `observables(include_figures=False)`는 약 `0.35 ms`, figure 포함 `observables()`는 평균 약 `1.2 s`였다. 즉 자동 보고서, 테스트, sweep, batch scan은 가능하면 `headless_observables()` 또는 `observables(..., include_figures=False)`를 호출해야 한다 (`gabes/schemes/base.py:92-106`, `gabes/schemes/rydberg.py:496-556`).
- normal UI는 이미 `cache_observables=True`와 `_cached_observables()`를 사용한다 (`gabes/schemes/base.py:68-76`, `streamlit_app.py:465`, `streamlit_app.py:757-758`). 추가로 extra view나 batch report 쪽에서 figure 생성을 opt-in으로 더 엄격히 분리하면 원래 기능을 침해하지 않고 체감 속도를 올릴 수 있다.
- `doppler=on` 반복 스캔에서는 `doppler.velocity_grid(T, mass, dv, cutoff_sigma)` 결과 캐시가 작게 도움이 될 수 있다 (`gabes/doppler.py:15`, `gabes/schemes/rydberg.py:353-355`). 다만 병목의 대부분은 velocity-class solve라서 grid 캐시만으로는 큰 개선은 어렵다.
- `_atom()`은 topology skeleton을 캐시하지만 dephasing이 바뀔 때마다 `AtomModel`을 새로 만들고 dissipator를 다시 구성한다 (`gabes/schemes/rydberg.py:58-78`, `gabes/schemes/rydberg.py:200-219`, `gabes/atoms.py:88-98`). base dissipator와 dephasing delta를 분리 캐시하는 최적화는 가능하지만, 현재 default compute가 수십 ms 수준이라 우선순위는 figure/headless 경로보다 낮다.

## 검증

실행:

```bash
python -m pytest tests/test_rydberg_eit.py -q
python -m pytest tests/test_kernels.py -q
```

결과:

- `tests/test_rydberg_eit.py`: `13 passed in 29.43s`
- `tests/test_kernels.py`: `14 passed in 43.51s`

## 결론

오늘 기준 `RydbergEITScheme`은 이전 6월 22일 리뷰 때의 "display-only knob가 많다"는 약점을 상당히 벗어났다. probe/coupling power, beam diameter, transit broadening, finite IF discriminator, opt-in temperature dephasing까지 연결되어 있어 실험 설계자가 knob 방향성을 확인하는 reference로는 꽤 쓸 만하다.

다만 절대 electrometry sensitivity 또는 최종 calibration reference로 부르기에는 아직 static/lumped 모델이다. 다음으로 비용 대비 효과가 좋은 방향은 full Zeeman이나 full demodulation으로 바로 가는 것이 아니라, headless/batch 경로를 적극 쓰면서 IF proxy와 low-order polarization/Zeeman distortion을 조심스럽게 확장하는 것이다.
