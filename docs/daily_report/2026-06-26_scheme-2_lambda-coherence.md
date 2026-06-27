# 2026-06-26 Scheme 2 Review: Lambda coherence

## 오늘 선택된 스킴

- 오늘 현지 날짜는 `2026-06-26`이고 day-of-month는 `26`이다.
- 규칙 `n = (day mod 5) + 1`에 따라 `n = (26 mod 5) + 1 = 2`이므로 오늘 검토 대상은 2번째 스킴이다.
- 현재 등록된 다섯 개 스킴 순서는 `OD / SAS -> Lambda coherence -> Rydberg-EIT electrometry -> Hanle / EIA / NMOR -> FWM`이다 (`gabes/schemes/__init__.py:19-24`, `README.md:8-16`).
- 따라서 오늘 검토 대상은 `LambdaScheme`이 담당하는 `Lambda coherence (EIT / AT / CPT)`이다 (`gabes/schemes/absorption.py:447-556`).

## 이번 검토에서 읽은 범위

- 스킴 등록/순서: `gabes/schemes/__init__.py:19-24`
- Lambda 스킴 본체: `gabes/schemes/absorption.py:447-678`
- Doppler 평균 및 affine kernel 경로: `gabes/schemes/absorption.py:44-140`
- 매질 상수/밀도/파수 선택: `gabes/schemes/absorption.py:407-432`
- 3준위 원자 모델: `gabes/atoms.py:135-149`
- 흡수/투과/군지수 readout: `gabes/observables.py:364-405`
- 물리 테스트: `tests/test_absorption.py:109-158`
- numba affine fast path 일치성 테스트: `tests/test_kernels.py:102-136`

## 구현이 실제로 무엇을 하고 있는가

- 이 스킴은 hyperfine-resolved manifold가 아니라 **대칭 branching을 가진 축약 3준위 Lambda 모델** 하나를 EIT/AT/CPT 세 regime로 보여 준다. excited state는 하나이고, decay는 두 바닥상태로 `Gamma/2`씩 나뉘며, ground coherence decay는 `buffer_ground_relax_khz` 하나로 대표된다 (`gabes/atoms.py:135-149`).
- scan variable은 probe detuning이며, Hamiltonian에서 `H[1,1] = Dc - s`로 들어가므로 coupling detuning을 중심으로 한 two-photon resonance 구조를 읽는다 (`gabes/schemes/absorption.py:592-599`).
- Doppler-on일 때도 모델은 **단일 excited-state shift + co-propagating geometry**를 가정해 two-photon resonance를 정확히 Doppler-free로 둔다. 즉 velocity averaging은 broad optical pedestal에는 들어가지만 dark resonance 자체는 residual `k_p-k_c` mismatch 없이 계산된다 (`gabes/schemes/absorption.py:9-17`, `59-72`, `117-140`).
- 온도는 매질 density, Maxwell velocity width, optical wavenumber 선택에 들어가지만, Lambda 스킴 자체에는 SAS 쪽처럼 self-broadening이나 pressure shift 모델이 없다 (`gabes/schemes/absorption.py:407-432`, `578-608`).
- `cell_mm`는 `recompute=False`라서 무거운 `chi_bar` solve를 다시 돌리지 않고 Beer-Lambert transmission만 다시 그린다. 실험자가 셀 길이만 바꿔 optical depth 감을 볼 때는 이 선택이 좋다 (`gabes/schemes/absorption.py:523-524`, `617-619`).

## 테스트와 직접 실행 결과

- `python -m pytest tests/test_absorption.py tests/test_kernels.py -q`를 실행했고 `24 passed in 8.61s`로 통과했다.
- 저장소 테스트는 다음 정도를 보장한다.
  - AT splitting이 coupling Rabi와 거의 1:1로 맞는지 (`tests/test_absorption.py:109-117`)
  - cold EIT에서 공명 투과가 충분히 커지는지 (`tests/test_absorption.py:120-128`)
  - CPT dark resonance가 natural linewidth보다 좁은지 (`tests/test_absorption.py:130-145`)
  - numba affine kernel이 NumPy reference와 일치하는지 (`tests/test_kernels.py:107-136`)
- 직접 돌린 기본 동작점 결과는 다음과 같았다.
  - `EIT` 기본값: `Transmission at resonance ~= 0.014`, `Window FWHM ~= 0.46 MHz`, `Group index ~= 1.0e5`
  - `AT` 기본값: `AT splitting ~= 46.0 MHz`, expected `Omega_c ~= 46.0 MHz`, center transmission `~= 0.989`
  - `CPT` 기본값: `Transmission at resonance ~= 0.706`, `Window FWHM ~= 923 kHz`
- warm 상태에서 `EIT compute()`는 약 `70 ms`였고 `cell_mm = 3, 15, 60 mm`로 바꿔도 거의 같았다. 반면 `observables()`는 약 `550 ms`로 figure 생성이 더 큰 비용이었다. 즉 인터랙티브 체감의 주병목은 solver보다 plotting 쪽이다.
- 한편 EIT 기본값에서 `temp_c: 50 -> 90 C`로 올려도 창 폭은 거의 `0.46 MHz`로 유지되고 공명 투과만 크게 떨어졌다. 코드 구조상 temperature가 density와 Doppler pedestal에는 강하게 들어가지만 dark-state coherence linewidth budget에는 거의 안 들어간다는 뜻이다 (`gabes/schemes/absorption.py:565-567`, `579-608`).

## 실험 원자광학 연구자 관점의 평가

- **장점은 분명하다.** textbook Lambda physics를 빠르게 확인하는 용도에는 매우 좋다. coupling Rabi를 올리면 EIT에서 AT로 자연스럽게 넘어가고, AT splitting이 `Omega_c`를 따라가는지 test까지 걸려 있다 (`gabes/schemes/absorption.py:650-665`, `tests/test_absorption.py:109-117`).
- `Transmission`, `Re chi`, `group index`를 한 자리에서 같이 보여 주는 점도 실험가에게 유용하다. “투명창이 생기는가”뿐 아니라 “분산 기울기가 실제로 가파른가”를 바로 볼 수 있기 때문이다 (`gabes/schemes/absorption.py:621-647`, `gabes/observables.py:397-405`).
- `cell_mm`를 navigate-only knob로 둔 선택도 실험 planning에는 실용적이다. 같은 susceptibility에서 셀 길이만 바꿔 transmission/OD 감도를 빨리 스캔할 수 있다 (`gabes/schemes/absorption.py:523-524`, `617-619`).

## 어디까지 실험물리학적 레퍼런스로 쓸 수 있나

- 오늘 기준 이 코드는 **정성적 내지 반정량적 Lambda spectroscopy reference**로는 충분히 쓸 만하다.
- 특히 다음 질문에는 잘 답한다.
  - coupling Rabi를 이 정도 주면 AT 분리가 몇 MHz쯤 나와야 하는가
  - ground-coherence relaxation을 올리면 EIT/CPT 창이 얼마나 무뎌지는가
  - 같은 susceptibility에서 셀 길이를 바꾸면 transmission이 얼마나 급해지는가
  - slow-light sign과 relative steepness가 어느 regime에서 가장 큰가
- 그러나 **실험의 절대 linewidth, absolute contrast, isotope-resolved line assignment를 바로 믿는 1차 기준 모델**로 쓰기에는 아직 부족하다.
- 이유는 다음과 같다.
  - 원자 구조가 hyperfine/Zeeman-resolved가 아니라 대칭 3준위 lumped model이다 (`gabes/atoms.py:135-149`).
  - Doppler 처리가 co-propagating exact-cancellation을 가정하므로 residual angle mismatch, counter-propagating Raman, unequal `k`에 따른 two-photon washout이 없다 (`gabes/schemes/absorption.py:9-17`, `117-140`).
  - temperature가 density와 Doppler pedestal을 바꾸지만, 실제 warm vapor EIT/CPT에서 흔한 collision-enhanced dephasing, pressure shift, Dicke narrowing, optical-pumping redistribution은 직접 들어가지 않는다 (`gabes/schemes/absorption.py:407-432`, `578-608`).
  - 테스트도 paper-anchored warm-vapor EIT linewidth/contrast를 맞추는 수준까지는 아직 가지 않고, 주로 내부 일관성과 textbook invariant를 본다 (`tests/test_absorption.py:109-158`).

## 저장소에 이미 있던 개선안과 계산 부하

- 이번 검토에서는 Lambda 스킴에 대해 `docs/checklist.json`, 기존 `docs/daily_report`, `README.md`, `CLAUDE.md` 안에서 **별도로 명시된 기존 개선안은 찾지 못했다**.
- 즉 오늘 보고서의 Lambda 개선 제안은 새로 제안하는 후보이며, 기존 저장소 TODO를 이어받는 항목은 아니다.

## 조심스럽게 제안할 수 있는 물리 개선 후보

### 1. `coupling_rabi_mhz`에 대응하는 실험 노브를 추가하기

- 지금 스킴은 coupling을 바로 `MHz`로 받는다 (`gabes/schemes/absorption.py:511-516`).
- 실험가는 실제로는 power와 beam size를 조절하므로, Rydberg 스킴처럼 `coupling_power_mw`와 waist/diameter에서 `Omega_c`를 유도하는 보조 입력을 두면 lab-facing reference 가치가 커진다.
- 계산 부하는 거의 없다. 대부분은 solve 전에 `Omega_c`를 산술 변환하는 단계라 FLOP 증가가 미미하다.
- 다만 UI상 그 노브들이 `recompute=True` solve knob가 되므로, slider 체감은 지금보다 조금 무거워질 수 있다.

### 2. buffer-gas / pressure-shift / Dicke narrowing의 저차 모델을 Lambda에도 넣기

- 현재 Lambda에는 buffer pressure knob 자체가 없고, optical linewidth도 사실상 species natural linewidth만 쓴다 (`gabes/schemes/absorption.py:407-432`).
- warm vapor EIT/CPT 실험에서는 pressure broadening, pressure shift, Dicke narrowing이 line center와 dark-feature width를 함께 바꾼다.
- 이 항목은 **저부하로 넣기 좋다**. line center offset, homogeneous width correction, phenomenological Dicke narrowing factor 정도는 scalar correction이라 solve 차수는 그대로다.
- full velocity-changing-collision kernel까지 가면 다른 문제지만, 그 전 단계는 속도 손실이 거의 없으면서 물리를 상당히 따라갈 수 있다.

### 3. residual two-photon Doppler / beam-angle mismatch knob 추가

- 현재 구현의 핵심 가정은 “two-photon resonance는 Doppler-free”이다 (`gabes/schemes/absorption.py:117-123`).
- 실제 셋업에서 probe/coupling angle이 완전 0이 아니거나 파장이 충분히 다르면 dark resonance가 퍼진다.
- `delta_k_eff` 또는 beam angle 기반 residual mismatch를 옵션으로 두면, warm-vapor EIT/CPT의 실험적 민감도를 훨씬 현실적으로 줄 수 있다.
- 이 역시 Hamiltonian과 velocity shift에 affine하게만 들어가게 유지하면 계산 부하는 크지 않을 가능성이 높다.

### 4. hyperfine-resolved Lambda manifold + optical pumping

- 이건 물리적으로는 가장 강력한 업그레이드다. 특히 natural Rb, Cs에서 “어느 hyperfine pair를 보고 있는가”, “optical pumping 때문에 contrast가 왜 이렇게 바뀌는가”를 더 직접 설명할 수 있다.
- 하지만 계산 부하는 작지 않다. 지금 3준위는 density-matrix 차원이 `3^2 = 9`인데, 예를 들어 8준위만 돼도 `8^2 = 64`가 된다 (`gabes/atoms.py:22-54`, `135-149`).
- batched linear solve의 비용은 대략 차원 세제곱에 민감하므로, 인터랙티브 속도는 여러 배에서 한두 자릿수 이상까지 느려질 수 있다.
- 따라서 이 항목은 “부하를 거의 만들지 않으면서” 넣는 후보로는 보기 어렵다.

## 순수 코딩 측면의 속도 개선 후보

### 1. metric/table 경로와 Matplotlib figure 경로를 분리하기

- 현재 `observables()`는 항상 figure 2장을 만든다 (`gabes/schemes/absorption.py:610-647`).
- 실측으로는 EIT warm solve가 약 `70 ms`인데 `observables()`가 약 `550 ms`였다. 이 스킴의 병목은 solver보다 plotting이다.
- metric/table만 필요할 때 figure 생성을 건너뛰거나, Streamlit에서 collapsed 상태의 그림은 lazy-render하게 하면 원래 기능을 해치지 않고 체감 속도를 크게 올릴 수 있다.

### 2. `atoms.lambda3(...)`와 매질 파생상수 캐시

- 현재 `compute()`마다 `atoms.lambda3(...)`를 새로 만들고 (`gabes/schemes/absorption.py:587`), `_medium_from_params()`도 species/line/temp에서 `dipole`, `k_vec`, `omega0`, `mass`를 다시 계산한다 (`gabes/schemes/absorption.py:407-432`).
- 이 비용이 지금의 주병목은 아니지만, 배치 스캔이나 자동 보고서 루프에서는 누적된다.
- `(gamma, gamma_gg)`와 `(species, line, temp_c)` 키 기반의 작은 cache를 두면 원래 물리를 전혀 바꾸지 않고 미세한 속도 이득을 얻을 수 있다.

### 3. `chi_bar -> alpha/xphys` 중간결과 재사용

- `cell_mm`는 navigate-only knob인데, 그때도 `observables.absorption_coefficient()`가 매번 다시 호출된다 (`gabes/schemes/absorption.py:614-618`, `gabes/observables.py:381-384`).
- `alpha`와 `xphys`를 raw 또는 별도 cache에 유지하면 cell-length 슬라이더 이동 시 다시 계산할 것이 Beer-Lambert transmission뿐이어서 더 가벼워질 수 있다.
- 이 역시 물리 모델을 건드리지 않는 안전한 최적화다.

## 종합 판단

- 오늘 기준 `LambdaScheme`은 **교과서적 EIT/AT/CPT를 빠르게 스캔하고 실험 조건 변화의 방향성을 읽는 용도**로는 충분히 유용하다.
- 특히 AT splitting 대 `Omega_c`, CPT 협폭성, group-index sign/steepness 같은 핵심 직관을 빠르게 주는 점은 강점이다.
- 다만 실제 warm-vapor 실험의 절대 linewidth, residual Doppler washout, pressure-shifted line center, hyperfine optical-pumping contrast까지 그대로 믿고 가져갈 정도의 최종 레퍼런스는 아직 아니다.
- 기존 저장소에 Lambda 전용 개선안은 아직 보이지 않았고, 오늘 기준 가장 비용 대비 효과가 좋은 새 후보는 **lab knob 기반 `Omega_c` 입력**과 **저차 buffer-gas / residual Doppler 보정**이다. 이 둘은 계산 부하를 거의 늘리지 않으면서 실험가가 체감하는 현실성을 가장 많이 올릴 가능성이 크다.
