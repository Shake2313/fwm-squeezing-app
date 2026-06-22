# 2026-06-22 Scheme 3 Review: Rydberg-EIT electrometry

## 선택 결과

- 오늘의 현지 날짜는 `2026-06-22`이므로 `n = (22 mod 5) + 1 = 3`이다.
- 현재 등록된 다섯 개 스킴 순서는 `OD / SAS -> Lambda coherence -> Rydberg-EIT electrometry -> Hanle / EIA / NMOR -> FWM`이다 (`gabes/schemes/__init__.py:20-24`, `README.md:12-16`).
- 따라서 오늘 검토 대상은 3번째 스킴 `RydbergEITScheme`이다 (`gabes/schemes/rydberg.py:24-31`).

## 이번 검토에서 읽은 근거

- 스킴 등록 및 순서: `gabes/schemes/__init__.py:12-24`
- 스킴 본체: `gabes/schemes/rydberg.py:24-289`
- 흡수/투과 읽기: `gabes/observables.py:287-312`
- 85Rb 증기 밀도식: `gabes/species.py:192-207`
- affine fast path: `gabes/kernels.py:299-350`
- 기준 테스트: `tests/test_rydberg_eit.py:28-72`, `tests/test_kernels.py:103-136`
- 기존 개선안 메모: `docs/checklist.json:4-10`

## 현재 물리 정의와 구현 범위

- 이 스킴은 정적 4준위 cascade ladder `5S1/2 F=3 -> 5P3/2 F'=4 -> 40D5/2 -> 39F7/2`를 푼다 (`gabes/schemes/rydberg.py:4-11`, `111-127`).
- 프로브는 항상 약한 선형 probe로 고정되고 `PROBE_RABI = 1e-3`가 `gamma_e`에 곱해진다 (`gabes/schemes/rydberg.py:20-21`, `151`).
- 실제 OBE를 움직이는 주된 노브는 `coupling_rabi_mhz`, `lo_rabi_mhz`, `mw_detuning_mhz`, `rydberg_dephasing_mhz`다 (`gabes/schemes/rydberg.py:77-90`, `148-156`).
- `probe_power_uw`, `coupling_power_mw`, `beam_diameter_mm`, `mw_frequency_ghz`, `if_khz`는 코드와 help 문구에서 명시적으로 display-only다 (`gabes/schemes/rydberg.py:60-72`, `85-95`, `202-206`, `269-277`).
- 도플러 처리는 현재 `off`만 허용된다. 즉 residual Doppler mismatch, velocity-class pumping, transit spread는 이 스킴에 직접 들어오지 않는다 (`gabes/schemes/rydberg.py:93-95`).
- 출력은 optical transmission/dispersive spectrum과 RF AT splitting까지만 다룬다. time-domain superheterodyne demodulation은 범위 밖이라고 코드가 직접 말한다 (`gabes/schemes/rydberg.py:8-11`, `281-288`).

## 직접 실행한 확인

실행 명령:

```bash
python -m pytest tests/test_rydberg_eit.py tests/test_kernels.py -q
```

결과:

- `18 passed in 8.48s`
- 따라서 기준값 테스트와 numba affine fast path 일치성은 현재 HEAD에서 통과한다 (`tests/test_rydberg_eit.py:42-72`, `tests/test_kernels.py:119-136`).

추가로 스킴을 직접 호출해 확인한 핵심 수치는 다음과 같다.

- EIT default에서 warm `compute()`는 약 `5-10 ms`, `observables()`는 약 `140-170 ms`였다. 첫 호출은 JIT/초기화 영향으로 `compute()` 약 `291 ms`, `observables()` 약 `821 ms`까지 올라갔다 (`gabes/schemes/rydberg.py:170-179`, `225-279`, `gabes/kernels.py:299-350`).
- EIT default linewidth는 약 `1.665 MHz`였고, 이는 테스트가 요구하는 `1.3-1.9 MHz` 범위 안이다 (`tests/test_rydberg_eit.py:42-48`).
- AT default에서 RF splitting은 약 `3.663 MHz`였고, 기본 `lo_rabi_mhz = 3.7`과 거의 일치했다 (`gabes/schemes/rydberg.py:45`, `257-259`, `tests/test_rydberg_eit.py:51-57`).

`chi_bar` 자체가 어떤 노브에 반응하는지도 따로 확인했다.

- `probe_power_uw`, `coupling_power_mw`, `beam_diameter_mm`, `mw_frequency_ghz`, `if_khz`를 바꿔도 `max|Δchi_bar| = 0.0`이었다. 코드 정의와 정확히 일치한다 (`gabes/schemes/rydberg.py:60-72`, `85-95`).
- `temp_c`와 `cell_mm`를 바꿔도 `chi_bar`는 변하지 않았다. 이는 `temp_c`가 `species.number_density()`를 통해 원자수밀도 `N`만 바꾸고 (`gabes/schemes/rydberg.py:146-147`, `gabes/species.py:192-207`), `cell_mm`가 observables 단계의 Beer-Lambert 길이 `L`만 바꾸기 때문이다 (`gabes/schemes/rydberg.py:193`, `231-235`, `gabes/observables.py:304-312`).
- 실제로 AT default에서 `temp_c: 20 -> 40 C`로 올리면 공진 투과가 `0.755 -> 0.105`, 최대 기울기가 `0.0507 -> 0.1307 /MHz`로 크게 바뀌지만 splitting은 `3.663 MHz`로 고정됐다. 즉 "증기 밀도에 따른 흡수 깊이"는 잡지만, ladder coherence 자체의 thermal broadening은 거의 안 잡는다.
- `cell_mm: 50 -> 70 mm`로 바꾸면 공진 투과는 `0.755 -> 0.674`로 바뀌지만 splitting은 역시 그대로다. 이것도 Beer-Lambert 후처리와 일치한다 (`gabes/schemes/rydberg.py:231-235`).
- 반면 `lo_rabi_mhz: 3.7 -> 6.0`은 splitting을 `3.663 -> 5.85 MHz`로 키웠고, `coupling_rabi_mhz: 3.0 -> 5.0`은 transparency 구조와 splitting을 함께 움직였다. 이 부분은 실제로 실험자가 기대하는 정적 AT 직관과 잘 맞는다 (`gabes/schemes/rydberg.py:152-165`, `255-259`).

## 실험물리 관점 평가

### 어디까지는 유용한가

- 이 스킴은 "정적 Rydberg-EIT/AT 스펙트럼을 빠르게 탐색하는 도구"로는 충분히 쓸 만하다. 4준위 OBE에 optical susceptibility, density, Beer-Lambert transmission을 일관되게 연결하고 있기 때문이다 (`gabes/schemes/rydberg.py:139-207`, `225-279`, `gabes/observables.py:287-312`).
- 특히 `lo_rabi_mhz`가 AT splitting으로 거의 1:1 대응하는지 바로 확인할 수 있어, RF dressing이 EIT window를 얼마나 찢는지 직관적으로 파악하는 실험 planning tool로 좋다 (`tests/test_rydberg_eit.py:51-57`).
- fast path가 affine 구조를 잘 이용하고 있어서 interactive scan 용도로 매우 빠르다 (`gabes/schemes/rydberg.py:170-179`, `gabes/kernels.py:299-350`).

### 어디서부터 레퍼런스로는 조심해야 하나

- 정량 electrometry 레퍼런스로는 아직 "반정량적 reference" 수준이 맞다. 가장 큰 이유는 실험자가 실제로 돌리는 beam power와 waist가 스펙트럼을 직접 움직이지 않기 때문이다 (`gabes/schemes/rydberg.py:60-72`, `77-82`; `docs/checklist.json:5-8`).
- `if_khz`와 superheterodyne demodulation 체인이 빠져 있으므로, 현재 metric의 `Max spectral slope`는 좋은 내부 proxy일 뿐 실제 readout sensitivity 자체는 아니다 (`gabes/schemes/rydberg.py:91-95`, `263-266`, `281-288`; `tests/test_rydberg_eit.py:60-72`).
- `temp_c`가 실제로는 수밀도만 바꾸고 coherence kernel에는 안 들어가므로, 실제 셀 온도 상승에서 흔한 residual Doppler washout, collision-enhanced dephasing, transit redistribution을 직접 추적하기 어렵다 (`gabes/schemes/rydberg.py:146-147`, `170-179`; `gabes/species.py:192-207`).
- Zeeman sublevel, polarization selection rules, stray electric/magnetic mixing, state-specific ionization/dephasing 분해가 없는 4준위 lumped ladder이므로, 절대선폭과 절대감도 피팅에 그대로 쓰기에는 정보가 부족하다 (`gabes/schemes/rydberg.py:111-127`).

정리하면:

- "정적 spectrum shape와 splitting trend를 보는 실험 설계 도구"로는 유용하다.
- "장비 power knob까지 포함한 정량적 electrometry calibration reference"로 쓰기에는 아직 한 단계 부족하다.

## 저장소에 이미 있던 개선안과 계산 부하

이번 검토에서 Rydberg 관련 기존 개선안으로 명시적으로 찾은 것은 하나였다.

- `rydberg-power-to-rabi`: beam power/geometry를 intensity와 Rabi frequency로 연결하자는 제안 (`docs/checklist.json:4-10`)

이 개선안을 반영할 때의 계산 부하는 다음처럼 본다.

- 핵심 solver 차원은 그대로 4준위다. 새 상태를 늘리거나 새 속도축을 도입하지 않으므로 FLOP 증가는 거의 없다.
- 추가되는 일은 power, waist, dipole에서 `Omega_c`를 계산하는 전처리 algebra가 대부분이므로, 순수 연산량 증가는 미미하다.
- 다만 UI 체감 속도는 조금 달라진다. 지금은 `coupling_power_mw`와 `beam_diameter_mm`가 `recompute=False`여서 slider 이동이 solve를 다시 안 돌지만 (`gabes/schemes/rydberg.py:65-72`), power-to-Rabi를 연결하면 이 노브들이 solve knob가 된다.
- 그래도 warm solve 자체가 이미 `5-10 ms` 수준이라, figure 생성 비용을 제외하면 사용성 악화는 크지 않을 가능성이 높다.

즉, 이 기존 개선안은 **물리 해석력 상승 대비 계산 부하가 매우 작은 편**이다. 오늘 기준 최우선 후보로 여전히 타당하다.

## 기존 제안이 커버하지 않는 조심스러운 후속 후보

기존 메모에 없는 항목 중에서는 다음 정도가 현실적이다.

- `rydberg_dephasing_mhz`를 단일 숫자 대신 `ground-40D`, `40D-39F` 두 채널 정도로 나누기.
  - 행렬 차원은 그대로라 계산 부하는 작다.
  - 실제 셀 상태 변화나 RF 잡음의 원인을 좀 더 실험적으로 분리해서 볼 수 있다.
- metric에 "AT center shift"를 추가하기.
  - 지금은 splitting만 바로 보이고 detuned microwave가 중심을 얼마나 끄는지는 직접 읽기 어렵다 (`gabes/schemes/rydberg.py:209-223`, `255-267`).
  - 계산 부하는 사실상 없다.
- residual Doppler를 옵션으로 추가하기.
  - 이건 기존 두 제안보다 무겁다. velocity class 수에 거의 비례해 비용이 늘고, 현재의 `kv=[0], weights=[1]` 구조가 깨진다 (`gabes/kernels.py:320-323`).
  - 그래도 커널 구조 자체는 이미 있어 완전히 새 solver를 만들 필요는 없다.

## 물리와 무관한 순수 코딩 최적화 관점

- 현재 Rydberg 스킴에서 큰 병목은 warm solve보다 `observables()`의 figure 생성이다. 직접 측정에서도 compute는 `5-10 ms`, figure 포함 observables는 `140-170 ms`였다 (`gabes/schemes/rydberg.py:225-279`).
- 따라서 원래 기능을 해치지 않고 속도를 올리려면 "숫자만 필요한 경로"와 "figure까지 필요한 경로"를 분리하는 것이 가장 효과적이다.
- 미세 최적화로는 `species.RB85.line("D2")`, `lam`, `k_vec`, `dipole` 같은 불변 상수를 모듈 또는 클래스 레벨에 캐시할 수 있다 (`gabes/schemes/rydberg.py:140-145`).
- `_atom()`은 dephasing 값만 바뀌는데 매번 동일 topology를 다시 조립한다. topology skeleton을 재사용하고 dephasing만 주입하는 방식도 가능하다 (`gabes/schemes/rydberg.py:111-127`).

## 결론

- 오늘 기준 GABES의 `RydbergEITScheme`은 **빠르고 안정적인 정적 Rydberg-EIT/AT 스펙트럼 탐색기**로는 충분히 실용적이다.
- 그러나 **실험자가 실제로 돌리는 power/waist/IF knob까지 반영하는 정량 electrometry 레퍼런스**로 보기에는 아직 제한이 뚜렷하다.
- 저장소에 이미 적혀 있던 `power-to-Rabi` 연결은 물리적 유용성을 가장 크게 올리면서도 계산 부하는 거의 늘리지 않는, 우선순위가 좋은 개선안이다 (`docs/checklist.json:4-10`).
