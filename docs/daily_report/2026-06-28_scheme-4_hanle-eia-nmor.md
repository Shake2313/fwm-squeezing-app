# 2026-06-28 Scheme 4 Review: Hanle / EIA / NMOR

## 오늘의 선택

- 현재 등록된 다섯 개 스킴 순서는 `OD / SAS -> Lambda coherence -> Rydberg-EIT electrometry -> Hanle / EIA / NMOR -> FWM`이다 (`gabes/schemes/__init__.py`, `README.md:8-16`).
- 오늘의 현지 날짜는 `2026-06-28`이므로 `n = (28 mod 5) + 1 = 4`이고, 검토 대상은 4번째 스킴 `MagnetoScheme`이다 (`gabes/schemes/__init__.py`).

## 이번 검토에서 읽은 근거

- 스킴 등록 순서와 개요: `gabes/schemes/__init__.py`, `README.md:8-16`
- Magneto 파라미터/기본 regime/UI 선언: `gabes/schemes/magneto.py:186-322`
- Magneto 물리 설명 및 참고문헌: `gabes/schemes/magneto.py:357-378`
- Magneto 계산 본체: `gabes/schemes/magneto.py:394-517`
- Magneto readout/metric/derived table: `gabes/schemes/magneto.py:607-699`
- 버퍼가스 broadening TODO: `gabes/constants.py:50-60`
- 저장소의 기존 개선안 메모: `docs/checklist.json:20-39`
- 기존 Magneto 일일 보고서: `docs/daily_report/2026-06-23_scheme-4_hanle-eia-nmor.md`
- 물리/회귀 테스트: `tests/test_magneto.py:154-242`, `tests/test_kernels.py:62-96`

## 현재 정의와 구현 수준

- 이 스킴은 `87Rb D1`의 near-zero-field Zeeman manifold를 대상으로, transmission 기반 Hanle/EIA와 rotation 기반 NMOR를 하나의 엔진에 묶는다 (`gabes/schemes/magneto.py:199-205`, `357-366`).
- regime는 `EIT dip`, `EIA peak`, `Buffer Hanle`, `Buffer LCA`, `NMOR` 다섯 개이며, merged entry 안에서 readout과 cell model까지 함께 바꾸는 구조다 (`gabes/schemes/magneto.py:174-197`, `292-327`, `342-352`).
- paraffin cell은 `light region <-> dark region` 2영역 교환 모델을 쓰며 `wall_coherence_ms`, `transit_relax_khz`, `dark_return_khz`가 협폭 Ramsey형 중심 구조를 결정한다 (`gabes/schemes/magneto.py:238-260`, `406-416`, `468-475`, `520-553`).
- buffer gas cell은 단일 영역 모델이고, `ne_pressure_torr`는 현재 optical homogeneous broadening으로만 직접 들어가며 ground-state relaxation은 `buffer_ground_relax_khz`, `collisional_depol_khz` override로 따로 받는다 (`gabes/schemes/magneto.py:232-255`, `394-405`; `gabes/constants.py:50-60`).
- probe polarization은 `qwp_deg`를 circular-basis drive weight로 바꿔 `sigma+ / sigma-` 구동 비를 조절한다. 그래서 EIT-like dip과 EIA/MIA-like peak의 전환이 “readout 선택”이 아니라 실제 line shape 결과로 나온다 (`gabes/schemes/magneto.py:107-124`, `213-215`, `436-438`, `634-640`).
- 중요한 물리적 장점은 spontaneous emission을 polarization-grouped jump operator로 묶는 Zeeman manifold를 써서 TOC(transfer of coherence) 경로를 살린 점이다. 이 때문에 `Fe = Fg + 1` 전이에서 intrinsic EIA가 실제로 나온다 (`gabes/zeeman.py:81-101`, `tests/test_magneto.py:179-192`).

## 지난 Magneto 보고서 이후 확인되는 변화

- 2026-06-23 보고서에서 조심스럽게 제안했던 longitudinal zero-offset 성격의 knob가 이제 실제로 구현됐다.
- 현재 `b_offset_ut`가 advanced parameter로 들어가고 (`gabes/schemes/magneto.py:219-223`), 계산 시 `b_physical_ut = b_ut + b_offset_ut`로 실제 물리 축에 더해지며 (`gabes/schemes/magneto.py:418-422`), 결과 table에도 노출된다 (`gabes/schemes/magneto.py:694-695`).
- 이 동작은 테스트로도 고정돼 있다 (`tests/test_magneto.py:170-176`).
- 실험물리학 관점에서 이 추가는 꽤 유용하다. 실제 Hanle/NMOR 셋업에서는 “스캔 중심은 0으로 명령했지만 shielding/coils offset 때문에 공진 중심이 밀린다”가 매우 흔한데, 이제 그 상황을 코드 안에서 직접 표현할 수 있다.
- 계산 부하는 사실상 0에 가깝다. `b_offset_ut`는 단지 `b_z`를 만들기 전 스칼라 shift 한 번 더하는 수준이라 solver 차원, velocity class 수, scan length를 바꾸지 않는다 (`gabes/schemes/magneto.py:418-422`, `456-478`).

## 테스트와 로컬 확인

- `python -m pytest tests/test_magneto.py tests/test_kernels.py -q`를 실행했고 `29 passed in 20.35s`였다.
- 테스트는 buffer mode broadening, intrinsic EIA, LCA의 transverse-field 의존성, NMOR zero crossing, kernel parity를 모두 확인한다 (`tests/test_magneto.py:154-242`, `tests/test_kernels.py:62-96`).
- 로컬 quick check에서도
  - default는 `Transmission at B=0 = 0.998`, `Zero-field feature = EIT-like dip`,
  - `qwp_deg = 45`는 `EIA/MIA-like peak`,
  - `b_offset_ut = 0.25`는 `crossover`로 바뀌고 `mean(b_physical_ut - b_ut) = 0.25`,
  - buffer 20 Torr preset류는 넓은 `EIT-like dip`
  를 보여, 모델이 실제로 “편광/offset/cell model”에 반응하고 있음을 다시 확인했다.

## 실험물리학자의 관점에서 본 평가

### 무엇이 실제 물리에 유용한가

- 이 코드는 장난감 3준위 Hanle 곡선이 아니라, 실제 `87Rb D1` hyperfine 전이 선택, Zeeman manifold, polarization selection, branching decay, TOC를 넣은 compact warm-vapor 모델이다 (`gabes/schemes/magneto.py:206-209`, `394-416`, `499-516`; `gabes/zeeman.py:81-101`).
- paraffin cell을 2영역 light/dark exchange로 나눈 점이 특히 좋다. 코팅 셀 Hanle/NMOR에서 broad pedestal 위에 narrow central feature가 어떻게 생기는지 설명하기에 적절하다 (`gabes/schemes/magneto.py:406-416`, `468-475`, `520-553`).
- `qwp_deg`, `residual_transverse_b_ut`, `transition`을 함께 바꾸며 EIT/EIA/LCA 부호 전환을 보는 것은 실험 정렬과 편광 점검에 매우 실용적이다 (`gabes/schemes/magneto.py:213-223`, `261-270`, `292-327`; `tests/test_magneto.py:179-211`).
- NMOR가 단순히 transmission 재포장이 아니라 `chi_+ - chi_-`의 실수부 차이로 rotation을 만들고 zero-crossing slope를 metric으로 주는 점도 실험 readout과 잘 맞는다 (`gabes/schemes/magneto.py:624-631`, `648-663`).

### 어디까지 레퍼런스로 쓸 수 있는가

- 오늘 기준 이 스킴은 **Hanle/EIA/NMOR 현상 분류, 부호 전환 이해, 셀/편광/잔류장 민감도 탐색용 semi-quantitative reference**로는 충분히 쓸 만하다.
- 특히 “지금 내 신호가 EIT 쪽인가 EIA 쪽인가”, “잔류 transverse field가 LCA를 만들고 있나”, “wall coherence를 어느 정도로 잡아야 이 정도 협폭이 나오나” 같은 실험실의 1차 sanity check에는 유용하다.
- 반면 **정량적인 magnetometer calibration reference**로 보기에는 아직 이르다.
- 가장 큰 이유는 buffer-cell physics가 아직 압축돼 있기 때문이다. `ne_pressure_torr`는 optical broadening만 직접 바꾸고 (`gabes/schemes/magneto.py:394-396`; `gabes/constants.py:50-60`), pressure가 실제로 바꾸는 diffusion-limited ground relaxation, pressure shift, Dicke narrowing, velocity-changing collisions는 직접 모델링하지 않는다.
- paraffin cell 쪽도 `wall_coherence_ms`, `transit_relax_khz`, `dark_return_khz`가 현상론적 knob라서 absolute fit은 가능해도, coating 재질/beam geometry로부터 predictive하게 바로 계산하는 단계는 아니다 (`gabes/schemes/magneto.py:238-260`, `406-410`).

## 저장소에 이미 있던 개선안과 계산 부하

이번에는 기존 제안이 분명히 있다. Magneto와 직접 연결되는 항목은 최소 세 가지다.

### 1. `buffer-gas-pressure-shift`

- 위치: `docs/checklist.json:21-25`
- 내용: `neon_buffer_broadening()` 중심의 현재 단순 모델을 gas/species/line coefficient table로 확장하고, pressure shift와 phenomenological Dicke narrowing을 OD/SAS, Magneto, Lambda에 넣자는 제안이다.
- Magneto에 반영할 때의 계산 부하는 대체로 작다.
- pressure shift는 `laser_detuning_mhz` 또는 내부 `dL`에 pressure-dependent offset을 더하는 전처리 수준으로 넣을 수 있어 거의 공짜다 (`gabes/schemes/magneto.py:439`, `454-478`).
- Dicke narrowing도 low-order effective correction을 `gamma_opt` 또는 Doppler scalar 쪽에 넣는다면 solver 차원이나 velocity grid coupling을 건드리지 않아 부하가 매우 작다 (`gabes/schemes/magneto.py:395-396`, `447-450`).
- 반대로 full VCC까지 가면 현재의 separable `(B x velocity)` 구조가 깨지므로 부하가 커질 가능성이 높다 (`gabes/schemes/magneto.py:454-478`, `520-580`).
- 따라서 이 기존 개선안은 **pressure shift + phenomenological Dicke narrowing 정도까지는 거의 부하를 늘리지 않으면서 물리를 유의미하게 올릴 수 있는 좋은 후보**다.

### 2. `magneto-buffer-relaxation-map`

- 위치: `docs/checklist.json:35-39`
- 내용: `ne_pressure_torr`가 현재는 mainly optical broadening에만 영향을 주고, `buffer_ground_relax_khz`와 `collisional_depol_khz`는 독립 override로 남아 있는 문제를 완화하기 위해, pressure 기반 기본 매핑을 넣자는 제안이다.
- 이것도 계산 부하는 거의 없다.
- 구현은 pressure/temperature/cell 길이 기반의 경험식 또는 간단한 table lookup으로 기본값을 채우고, Advanced override를 그대로 두면 된다. 즉 heavy solve 전의 파라미터 변환 단계에서 끝난다.
- 실험물리적 가치는 크다. 지금은 “20 Torr 셀”을 고르면 linewidth budget의 핵심 일부를 사용자가 다시 손으로 넣어야 하는데, 이 개선이 들어가면 buffer 셀 preset의 현실감이 크게 오른다.

### 3. `figureless-observables-paths`

- 위치: `docs/checklist.json:28-32`
- 이건 직접적인 물리 개선안은 아니지만, 저장소에 이미 적혀 있는 공통 성능 개선안이다.
- Magneto에서는 특히 중요하다. 현재 `observables()`는 항상 Matplotlib figure를 만든다 (`gabes/schemes/magneto.py:607-679`).
- solver보다 plotting이 느린 경우가 반복 보고됐고, 이 항목은 물리를 바꾸지 않으면서 batch sweep, 테스트, 보고서 자동화의 wall-clock을 줄일 가능성이 크다.

## 오늘 기준으로 조심스럽게 제안할 수 있는 추가 물리 개선

기존 checklist가 이미 상당히 핵심을 짚고 있어서, 새 제안은 많이 덧붙이기보다 빈 곳만 적는다.

### 1. buffer pressure가 ground relaxation으로 자동 이어지는 기본 경로 강화

- 사실상 `magneto-buffer-relaxation-map`의 구체화다.
- `ne_pressure_torr -> buffer_ground_relax_khz / collisional_depol_khz` 기본값을 넣고, 고급 사용자는 override하게 두는 방식이 현재 구조와 가장 잘 맞는다 (`gabes/schemes/magneto.py:232-255`, `394-405`).
- 계산 부하는 거의 0이고, 실험적 유용성은 크다.

### 2. pressure shift + low-order Dicke narrowing의 저부하 도입

- 이것도 기존 checklist를 Magneto 관점에서 재확인한 것이다.
- full collision kernel 대신 `gamma_opt` 또는 Doppler scalar correction 수준에서 시작하면 interactive speed를 지키면서 buffer-cell reference 가치를 높일 수 있다 (`gabes/constants.py:50-60`, `gabes/schemes/magneto.py:395-396`, `447-450`).

## 순수 코딩 측면의 속도 개선 후보

- 가장 우선순위가 높은 것은 여전히 **figureless/headless observables 경로**다. 지금 `observables()`는 invalid branch를 포함해 항상 figure 생성 쪽으로 들어간다 (`gabes/schemes/magneto.py:607-679`). 기능을 침해하지 않고 배치 분석 속도를 올리기 좋다.
- `_hamiltonian()`은 호출마다 `zeeman.angular_momentum_matrices()`를 다시 만든다 (`gabes/schemes/magneto.py:499-505`). `(Fg, Fe)` 또는 블록 크기별 캐시를 두면 물리를 건드리지 않고 미세한 반복 비용을 줄일 수 있다.
- 본질적인 큰 비용은 여전히 `(scan_points x velocity_classes)`와 paraffin 2영역 solve다 (`gabes/schemes/magneto.py:418-478`, `520-580`). 따라서 Python 미세 최적화보다
  - figure 생성 분리,
  - 작은 행렬 ingredient 캐시,
  - preset/default에서 과도한 `scan_points`, `velocity_classes`를 피하는 쪽
  이 실효성이 더 크다.

## 최종 판단

- 오늘 기준 `MagnetoScheme`은 **실험실에서 Hanle/EIA/NMOR 데이터를 읽고 정렬/편광/잔류장/셀 완화의 방향성을 판단하는 데 꽤 유용한 semi-quantitative reference**다.
- 특히 paraffin two-region model, TOC 기반 intrinsic EIA, NMOR slope readout, 그리고 이제 들어간 `b_offset_ut`까지 합치면 “실험가가 실제로 부딪히는 near-zero-field 문제”를 꽤 잘 담고 있다.
- 다만 buffer-gas cell의 pressure-dependent relaxation, pressure shift, Dicke narrowing이 아직 얕기 때문에, **절대적 linewidth/offset/sensitivity를 바로 믿는 정량 기준 모델**로 쓰기에는 아직 조심해야 한다.
- 저장소에 이미 적혀 있는 개선안 중에서는 `buffer-gas-pressure-shift`와 `magneto-buffer-relaxation-map`이 물리 가치 대비 계산 부하가 가장 좋고, 순수 성능 면에서는 `figureless-observables-paths`가 가장 즉효다.
