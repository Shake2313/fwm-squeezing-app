# 2026-06-23 Scheme 4 Review: Hanle / EIA / NMOR

## 선택 결과

- 현재 등록된 다섯 개 스킴 순서는 `OD / SAS -> Lambda coherence -> Rydberg-EIT electrometry -> Hanle / EIA / NMOR -> FWM`이다 (`gabes/schemes/__init__.py:19-24`, `README.md:10-16`).
- 오늘 현지 날짜는 `2026-06-23`이므로 `n = (23 mod 5) + 1 = 4`이고, 검토 대상은 4번째 스킴 `MagnetoScheme`이다 (`gabes/schemes/magneto.py:127-172`).

## 이번 검토에서 읽은 근거

- 스킴 등록/순서: `gabes/schemes/__init__.py:12-24`
- 스킴 설명과 참고문헌: `gabes/schemes/magneto.py:350-369`
- 파라미터/기본 regime 정의: `gabes/schemes/magneto.py:186-348`
- 계산 본체: `gabes/schemes/magneto.py:372-597`
- 관측량/metric/derived table: `gabes/schemes/magneto.py:599-697`
- Zeeman manifold 및 TOC(transfer of coherence) 구현 이유: `gabes/zeeman.py:71-103`
- kernel 가속 경로: `gabes/kernels.py:200-280`
- 물리 테스트: `tests/test_magneto.py:111-230`, `tests/verify_hanle_eit_eia.py:119-206`
- kernel parity 테스트: `tests/test_kernels.py:62-96`
- 기존 개선안 메모: `docs/checklist.json:18-23`

## 현재 정의와 물리 내용

- 이 스킴은 `87Rb D1`의 근영자기장 Zeeman manifold를 대상으로, transmission 기반 Hanle/EIA와 rotation 기반 NMOR를 하나의 엔진으로 묶는다 (`gabes/schemes/magneto.py:127-172`, `350-369`).
- 기본 readout은 `Transmission` 또는 `NMOR rotation`이고, preset/regime는 `EIT dip`, `EIA peak`, `Buffer Hanle`, `Buffer LCA`, `NMOR` 다섯 개다 (`gabes/schemes/magneto.py:284-322`).
- 파라핀 셀은 `light region <-> dark region` 2영역 교환 모델을 쓰고, `wall_coherence_ms`, `transit_relax_khz`, `dark_return_khz`로 Ramsey형 협폭 구조를 만든다 (`gabes/schemes/magneto.py:401-468`, `512-555`).
- buffer-gas 셀은 단일 영역 모델이며, `buffer_ground_relax_khz`와 `collisional_depol_khz`를 통해 ground-state relaxation을 넣는다 (`gabes/schemes/magneto.py:393-400`, `558-580`).
- 편광은 QWP 각도를 원형기저 가중치로 바꿔 `σ+ / σ-` 구동 비를 조절한다 (`gabes/schemes/magneto.py:113-124`, `431`, `490-509`).
- 중요한 물리 포인트는 spontaneous emission을 polarization-grouped jump operator `Σ_q`로 묶어 TOC를 살렸다는 점이다. 이 덕분에 `Fe = Fg + 1`에서 intrinsic EIA가 나오는 경로를 코드가 실제로 가진다 (`gabes/zeeman.py:94-101`, `tests/test_magneto.py:170-183`, `tests/verify_hanle_eit_eia.py:121-173`).

## 로컬 검증 및 간단 실험

- `python -m pytest tests/test_magneto.py tests/test_kernels.py -q`를 실행했고 `28 passed`였다.
- `tests/test_magneto.py`는
  - 선형 편광 파라핀 셀에서 zero-field dip,
  - 원편광 쪽으로 갈 때 peak 전환,
  - wall coherence 증가 시 협폭화,
  - buffer mode에서 broad Hanle 유지,
  - `Fe=Fg+1`에서 intrinsic EIA,
  - transverse field가 있어야 circular-light LCA가 나타남,
  - NMOR zero crossing
  을 확인한다 (`tests/test_magneto.py:123-214`).
- `tests/test_kernels.py`는 paraffin, buffer, NMOR 케이스에서 numba kernel과 NumPy 기준 해가 일치함을 확인한다 (`tests/test_kernels.py:62-96`).

직접 돌려 본 quick experiment 결과:

- 파라핀 Hanle default 근처에서 warm `compute()`는 대체로 `219-259 ms`, `observables()`는 `115-135 ms`였다. figure 생성까지 합치면 대략 `0.34-0.39 s` 수준이다 (`gabes/schemes/magneto.py:599-697`).
- 같은 조건에서 `doppler off / velocity_classes=1`이면 `compute() ~82 ms`, `doppler on / 5 class`면 `~156 ms`, `9 class`면 `~215 ms`, `scan_points 201`이면 `~552 ms`였다. 즉 현재 비용의 주축은 `(B scan points) x (velocity classes)`와 2영역 solve다 (`gabes/schemes/magneto.py:413-476`, `512-580`, `gabes/kernels.py:200-280`).
- buffer Hanle는 단일 영역이라 같은 `121 x 9` 조건에서 `compute() ~54 ms`로 paraffin보다 빨랐다 (`gabes/schemes/magneto.py:469-471`, `558-580`).
- 파라핀 `F=2 -> F'=1`에서 QWP `0 deg -> 45 deg`로 바꾸면 zero-field amplitude가 `-0.020 -> +0.073`으로 바뀌어 dip에서 peak로 뒤집혔다. 이는 코드 의도와 테스트 결과와 일치한다 (`gabes/schemes/magneto.py:113-124`, `tests/test_magneto.py:123-138`).
- 같은 전이에서 `wall_coherence_ms = 0.05 -> 10`으로 늘리면 중앙 폭이 `0.473 -> 0.151 uT`로 줄어들어, anti-relaxation coating 품질이 협폭 구조를 만든다는 실험적 직관과 잘 맞았다 (`gabes/schemes/magneto.py:405`, `512-555`, `tests/test_magneto.py:141-151`).
- buffer 셀에서 `ne_pressure_torr = 0 -> 20`으로 바꾸면 optical broadening은 `0 -> 78.2 MHz`가 되지만, ground relaxation은 pressure로부터 자동 계산되지 않고 여전히 사용자가 `buffer_ground_relax_khz`, `collisional_depol_khz`를 직접 잡아야 한다 (`gabes/constants.py:50-60`, `gabes/schemes/magneto.py:390-396`).

## 실험물리 관점 평가

### 잘 구현된 점

- 이 스킴은 단순한 장난감 3준위 Hanle 근사가 아니라, 실제 `87Rb D1` hyperfine 전이 선택, Zeeman manifold, `σ±` 편광 성분, branching decay, TOC까지 넣은 compact warm-vapor 모델이다 (`gabes/schemes/magneto.py:206-215`, `490-509`; `gabes/zeeman.py:71-103`).
- 파라핀 셀을 2영역 light/dark exchange로 나눈 것은 실험적으로 매우 유용하다. wall-preserved coherence가 broad transit pedestal 위에 narrow central feature를 얹는 구조는 코팅 셀 Hanle/NMOR 데이터를 읽을 때 핵심이다 (`gabes/schemes/magneto.py:358-360`, `401-468`, `512-555`).
- `Fe=Fg+1` intrinsic EIA와 circular-light LCA를 test로 명시해 둔 점은 좋다. 즉 코드가 “곡선이 예뻐 보인다” 수준이 아니라, 적어도 Moon/Yu/Lezama 계열의 부호 전환 physics를 의식적으로 맞추고 있다 (`gabes/schemes/magneto.py:361-369`, `tests/test_magneto.py:170-202`, `tests/verify_hanle_eit_eia.py:121-206`).
- NMOR도 단순히 transmission을 재포장한 것이 아니라 `χ_+ - χ_-`의 실수부 차이로 rotation을 만들어 slope를 readout metric으로 주고 있다 (`gabes/schemes/magneto.py:616-655`).

### 어디까지 레퍼런스로 쓸 수 있는가

- **정성적-반정량적 실험 레퍼런스**로는 꽤 쓸 만하다.
  - 어떤 전이에서 dip/peak가 나와야 하는지,
  - 파라핀 셀에서 wall coherence가 협폭폭에 어떻게 먹는지,
  - QWP/잔류 transverse field가 line shape를 어떻게 뒤집는지,
  - buffer cell과 paraffin cell이 왜 다른 폭과 contrast를 보이는지
  를 빠르게 점검하기에 좋다.
- 특히 “자기장 차폐가 충분한가”, “원편광 성분이 너무 커서 EIA 쪽으로 갔는가”, “wall coherence를 이 정도로 봐도 되는가” 같은 실험실의 1차 sanity check 용도로는 유용하다.

### 아직 조심해야 하는 한계

- **정량 magnetometer calibration reference**로 쓰기에는 아직 이르다.
- 가장 큰 이유는 buffer-gas 관련 physics가 아직 너무 압축돼 있기 때문이다. 현재 `ne_pressure_torr`는 optical homogeneous broadening만 바꾸고 (`gabes/schemes/magneto.py:390-391`, `gabes/constants.py:50-60`), 실제로 pressure가 바꾸는 diffusion, ground relaxation, pressure shift, Dicke narrowing은 직접 들어가지 않는다.
- 다시 말해 buffer 셀에서 실험자가 실제로 바꾸는 “압력”은 코드 안에서 곧바로 `buffer_ground_relax_khz`와 `collisional_depol_khz`로 이어지지 않는다 (`gabes/schemes/magneto.py:393-396`, `241-250`). 이 상태에서는 “20 Torr 셀이라면 linewidth가 얼마”를 정량 예측하기 어렵다.
- 파라핀 셀도 `wall_coherence_ms`, `transit_relax_khz`, `dark_return_khz`가 모두 phenomenological knob라서, absolute fit은 가능해도 재료/셀 geometry로부터 predictive하게 계산되는 단계는 아니다 (`gabes/schemes/magneto.py:231-259`, `401-405`).
- 따라서 현재 수준을 한 문장으로 요약하면, **실험적 현상 분류와 파라미터 감도 탐색에는 강하고, 절대적 자력계 성능 예측이나 buffer-cell 정밀 피팅에는 아직 보수적으로 써야 하는 코드**다.

## 저장소에 이미 있던 개선안과 계산 부하

이번 검토에서 magneto에 직접 연결되는 기존 개선안으로 명시적으로 찾은 것은 magneto 전용 TODO가 아니라 공통 항목 하나였다.

- `buffer-gas-pressure-shift`: Ne broadening을 gas/species/line table로 올리고 pressure shift와 Dicke narrowing까지 포함하자는 항목 (`docs/checklist.json:18-23`, `gabes/constants.py:50-60`)

이 개선안을 magneto 관점에서 나눠 보면 계산 부하는 다르게 볼 수 있다.

### 1. pressure shift만 넣는 경우

- `laser_detuning_mhz` 또는 내부 `dL`에 pressure-dependent offset을 한 번 더하는 수준이면 된다 (`gabes/schemes/magneto.py:432`, `447`).
- FLOP 증가는 사실상 무시 가능하다. 스캔당 giant solve 수는 그대로다.
- 실험물리적 가치는 꽤 있다. 특히 buffer 셀에서 line center alignment를 할 때 “왜 최적 detuning이 셀마다 조금씩 다르지?”를 설명하는 데 도움 된다.

### 2. Dicke narrowing을 phenomenological effective-width로 넣는 경우

- coarse Doppler 모델이나 `doppler_scale`, 혹은 effective optical width `gamma_opt` 쪽에 pressure-dependent correction을 넣는 방식이면 계산 부하는 매우 작다 (`gabes/schemes/magneto.py:390-391`, `440-443`).
- 완전한 미시적 VCC(velocity-changing collisions) 모델은 아니지만, 현재 코드 구조와 인터랙티브 속도를 거의 유지하면서도 실제 buffer-cell line shape를 훨씬 더 따라갈 수 있다.

### 3. velocity-class mixing을 가진 full collision kernel로 가는 경우

- 이 경우는 부하가 크다. 지금은 각 velocity class가 `deff = dL - kv`로 분리되어 batched solve가 가능하지만 (`gabes/schemes/magneto.py:447-476`), VCC가 들어오면 velocity classes끼리 결합되어 현재의 separable grid가 깨진다.
- 그러면 행렬 차원이 실질적으로 커지거나 iterative solve가 필요해져, interactive tool로서의 장점이 크게 줄 수 있다.

정리하면, **기존 checklist의 buffer-gas 개선안은 pressure shift + phenomenological Dicke narrowing 정도까지는 부하를 거의 늘리지 않으면서 물리를 상당히 개선할 수 있고, full VCC까지 가면 속도 희생이 커진다.**

## 추가로 조심스럽게 제안할 수 있는 물리 개선 후보

기존 checklist에 magneto 전용 항목은 보이지 않았지만, 현재 코드 구조를 보면 다음 두 가지는 우선순위가 높아 보인다.

### 1. `ne_pressure_torr -> buffer_ground_relax_khz / collisional_depol_khz`의 기본 매핑

- 지금은 pressure와 relaxation이 분리된 독립 노브다 (`gabes/schemes/magneto.py:227-250`, `393-396`).
- 실험자 입장에서는 “20 Torr면 대략 이 정도 diffusion-limited relaxation” 같은 기본값이 있어야 훨씬 직접적이다.
- 단순 경험식/테이블 기반 기본 매핑을 두고 Advanced에서 override만 허용하면, 계산 부하는 거의 0에 가깝고 실험적 사용성은 크게 오른다.

### 2. longitudinal offset field `B0`의 explicit knob

- 현재는 scan이 항상 `-Bmax ... +Bmax` 대칭이고, systematic field는 transverse residual만 따로 둔다 (`gabes/schemes/magneto.py:413-425`).
- 실제 실험에서는 영점 오프셋이 line shape 분류만큼 자주 문제를 만든다. `B0`를 넣어 zero crossing/peak center가 어긋나는 상황을 직접 맞추게 하면 실험 비교가 쉬워진다.
- 이것도 내부적으로는 `b_z` 축을 한 번 shift하는 일이라 계산 부하는 사실상 없다.

## 순수 코딩 측면의 속도 개선 후보

- 현재 warm path에서도 `observables()`가 `~115-135 ms`로 적지 않다. Matplotlib figure 생성이 매번 들어가기 때문이다 (`gabes/schemes/magneto.py:599-697`). metric/table 계산과 figure 생성 경로를 분리하면, figure가 필요 없는 자동 검증/배치 분석에서 바로 이득을 볼 수 있다.
- `_hamiltonian()`은 호출마다 `angular_momentum_matrices()`를 다시 만든다 (`gabes/schemes/magneto.py:490-499`). 이는 `(Fg, Fe)`별로 캐시해도 원래 물리를 바꾸지 않는다.
- `zeeman.zeeman_manifold()`의 구조적 부분(coupling topology, grouped emission skeleton)은 relaxation 값과 분리 가능한 부분이 있다 (`gabes/zeeman.py:71-103`). 이를 template화해서 재사용하면 repeated compute에서 파이썬 레벨 오버헤드를 줄일 수 있다.
- paraffin path는 2영역 `2M x 2M` solve라 본질적으로 비싸다 (`gabes/schemes/magneto.py:512-555`, `gabes/kernels.py:237-280`). 따라서 여기서는 미세한 Python 최적화보다 `scan_points`, `velocity_classes`, figure 생성, cached small-matrix ingredients` 쪽을 다듬는 편이 효과적이다.

## 최종 판단

- 오늘 기준 `MagnetoScheme`은 **실제 원자광학 실험물리학자가 Hanle/EIA/NMOR 실험을 해석할 때, 현상 판별과 파라미터 감도 탐색에 충분히 도움 되는 반정량적 레퍼런스**다.
- 특히 파라핀 셀의 협폭 Hanle/NMOR, TOC 기반 intrinsic EIA, circular-light LCA까지 한 UI 안에서 빠르게 비교할 수 있다는 점은 강점이다.
- 다만 buffer-gas physics가 아직 optical broadening 중심의 compact model이므로, **절대 linewidth/offset/sensitivity를 실험의 최종 숫자로 피팅하는 1차 레퍼런스**로 쓰기에는 보수적으로 접근하는 편이 맞다.
- 저장소에 이미 적혀 있던 `buffer-gas-pressure-shift` 개선안은 magneto에도 직접 가치가 크며, pressure shift와 phenomenological Dicke narrowing 정도는 연산 부하를 거의 늘리지 않고도 물리를 한 단계 끌어올릴 수 있다.
