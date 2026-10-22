# 2026-09-23 Scheme 4 물리 검토 — Hanle / EIA / NMOR

## 1. 선택과 현재 다섯 scheme

- 로컬 날짜: **2026-09-23**, Asia/Seoul.
- 계산: **n = (23 mod 5) + 1 = 4**.
- 대상: **Scheme 4, Magneto-optics (Hanle / EIA / NMOR)**.
- 현재 HEAD: **4cb011f**.
- 기준: 현재 working tree. 이미 있던 staged/unstaged/untracked 사용자 변경 보존.

실행 registry 순서:

| n | key | 정의 |
|---:|---|---|
| 1 | `sas` | Absorption spectroscopy. pump-off OD + pump-on hyperfine-pumping SAS |
| 2 | `lambda` | 3-level Λ coherence. EIT / AT / CPT |
| 3 | `rydberg_eit` | cascade Rydberg-EIT + microwave AT electrometry |
| 4 | `magneto` | ⁸⁷Rb D1 Zeeman Hanle / EIA / NMOR |
| 5 | `fwm` | seeded ⁸⁵Rb D1 double-Λ mean-field + generic SFWM biphoton |

근거: `gabes/schemes/__init__.py:14-24`, `README.md:8-16`. `all_schemes()`도 이 순서 반환.

## 2. 먼저 찾은 기존 제안, TODO, issue

검색 범위:

- Scheme 4 기존 일일 보고서 12개. 최신은
  `docs/daily_report/2026-09-18_scheme-4_hanle-eia-nmor.md`.
- `docs/checklist.json`의 Magneto, collision, Doppler, paraffin, headless 항목.
- `docs`, `analysis`, `gabes`, `tests` 안 `TODO`, `FIXME`, `issue`, `proposal`,
  `Hanle`, `EIA`, `NMOR` 검색.
- `.github`, 독립 TODO/issue 파일 검색.
- public GitHub Issues API. 현재 open/closed 합계 0.

별도 로컬 issue/TODO 파일 없음. 직접 관련 코드 TODO는 고정 Ne broadening을
gas/species/line coefficient table과 pressure shift로 승격하라는
`gabes/constants.py:79-89`.

**기존 개선 제안 있음. 새 제안부터 만들 상황 아님.** 현재 핵심:

| 항목 | 현재 상태 | 위치 |
|---|---|---|
| collisional coefficient provenance + pressure shift | P1, ready, large | `docs/checklist.json:2630-2662` |
| observable convergence / physical zero / trust label | P1, ready, medium | `docs/checklist.json:2848-2867` |
| geometry-aware buffer relaxation | P2, needs_decision, large | `docs/checklist.json:2870-2900` |
| full velocity-changing collision | P2, parked, research | `docs/checklist.json:3035-3056` |
| pole-residue Doppler reuse | P3, done | `docs/checklist.json:3263-3281` |
| exact-Doppler 뒤 EIT-dip preset 재검토 | P2, needs_decision, small | `docs/checklist.json:3299-3318` |
| paraffin transport/coating calibration | P2, parked, research | `docs/checklist.json:2728-2760` |

이전 보고서에만 남은 구체 제안도 존재:

- paraffin transit exchange와 in-beam pure dephasing 분리.
- circular/near-circular 입력에서 “polarization-plane rotation” hero gate.
- 생성 Hanle–FWM 예제 provenance와 현행 코드 재생성.

Scheme 4 production 코드와 테스트는 2026-09-14 이후 수정 시각 변화 없음.
9월 18일 검토 뒤 관련 commit도 없음. 이번 결과는 새 구현 평가가 아니라 현행 물리와
남은 신뢰 경계 재검증.

## 3. 현재 구현 물리

### 3.1 원자, 편광, 자기장

- ⁸⁷Rb D1, `Fg ∈ {1,2}`, `Fe ∈ {1,2}` 네 transition. line strength와 low-field
  hyperfine Landé `gF` 계산 (`gabes/schemes/magneto.py:21-30,55-104,384-400`).
- 선택한 hyperfine pair의 전 `mF` 상태 사용. σ⁺/σ⁻/π CG 결합과 CG-branched spontaneous
  emission 생성 (`gabes/zeeman.py:17-45,69-123`).
- 방출 jump를 편광별 한 연산자로 묶음. excited coherence → ground coherence TOC 보존.
  intrinsic EIA에 필요한 실제 물리 (`gabes/zeeman.py:101-108`;
  `tests/test_magneto.py:37-63,458-471`).
- QWP 각도 → complex circular-basis drive weight. longitudinal field, zero offset,
  residual transverse field를 같은 Zeeman Hamiltonian에 넣음
  (`gabes/schemes/magneto.py:107-118,208-219,430-450,555-574`).

### 3.2 완화와 셀

- buffer cell: natural emission 고정. Ne optical broadening은 elastic optical dephasing.
  reservoir reload와 ground pure dephasing 별도 (`gabes/schemes/magneto.py:401-417`;
  `gabes/zeeman.py:125-149`).
- paraffin cell: illuminated/dark 두 영역 교환. dark 영역에 wall coherence lifetime.
  light block은 trace-one conditional state, occupation은 diagnostic
  (`gabes/schemes/magneto.py:418-428,576-635`).
- susceptibility에 local vapor density 한 번만 적용
  (`gabes/schemes/magneto.py:530-552`; `gabes/observables.py:556-570`).

### 3.3 Doppler와 readout

- Doppler-on: 각 B row를 real Hermitian-basis affine shift-invert system으로 변환.
  pole-residue + Faddeeva로 **절단 없는 Maxwell 평균**. direct `v=0` guard 실패 row만
  ±6σ 수렴 grid fallback (`gabes/schemes/magneto.py:45-52,496-516,662-747`).
- 독립 옛 solver ±6σ, 1 m/s reference 대비 χ 최대 상대차 `1.6e-7`, guard 최대
  `4.9e-12`, fallback 0, 최소 71× 가속 기록
  (`analysis/pole_reuse/DEVLOG.md:136-158`; `analysis/pole_reuse/GATE.md:42-63`).
- 관측량: `α = k Imχ`, `T = exp(-αL)`,
  `θ = kL Re(χ+ - χ-)/4` (`gabes/schemes/magneto.py:788-823`).

즉 단순 Lorentzian 장난감 아님. selection rule, optical pumping, TOC, Zeeman precession,
transverse-field switching, two-region Ramsey core, Beer–Lambert/NMOR readout 연결.

## 4. 실험 원자물리 reference 판정

### 쓸 만한 범위

- linear/circular/elliptical drive에 따른 Hanle dip/peak/crossover 방향.
- transverse residual field, longitudinal zero offset 영향.
- paraffin wall coherence, dark return, buffer ground relaxation의 상대 추세.
- detuning/intensity/cell length/QWP/B-range 초기 설계.
- 측정 trace sign/topology sanity check와 반정량 fit seed.

코드 참고문헌도 Lee–Moon paraffin MIA, Moon–Kim Ramsey EIA/MIT, Yu–Lee–Moon
buffer LCA를 명시 (`gabes/schemes/magneto.py:362-381`).

### 절대 reference로 못 쓰는 이유

1. 선택한 hyperfine pair 안의 닫힌 모델. 실제 D1 open transition의 다른 ground F 누출,
   repump, 전체 excited-hyperfine interference 없음 (`gabes/zeeman.py:69-130`).
2. Ne collision은 `3.91 MHz/Torr` broadening scalar 하나. pressure shift, gas/line/T 의존,
   diffusion, spin destruction, Dicke narrowing, VCC 없음 (`gabes/constants.py:79-89`).
3. Gaussian beam, field gradient, cell geometry, diffusion 공간 평균 없음.
4. NMOR 식은 thin/weak-medium circular birefringence readout. 임의 ellipticity의 Jones/Stokes
   propagation 아님. 원편광에서는 “편광면” 자체가 정의되지 않음
   (`gabes/schemes/magneto.py:107-118,788-823`).
5. 절대 contrast/rotation/slope는 코드도 `external validation required` 표시
   (`gabes/schemes/magneto.py:41-43,825-857`).
6. sensitivity는 이 원자 모델만으로 안 나옴. optical power, detector/noise, measured calibration 필요.

판정: **실제 물리를 가진 정성·반정량 experimental-planning/reference**. publication-grade
absolute linewidth, rotation, contrast, pT/√Hz reference 아님.

## 5. 2026-09-23 수치 재검증

### 5.1 built-in 121-point exact-Maxwell

동일 process warm 5회 median.

| regime | feature | α(0)-wing [m⁻¹] | 중앙 폭 [µT] | T(Bcmd=0) | 3-point θ slope [mrad/µT] | median |
|---|---|---:|---:|---:|---:|---:|
| EIT dip | **crossover** | +0.357285 | 0.981271 | 0.787789 | -102.569 | 61.94 ms |
| EIA peak | peak | +14.930043 | 0.363514 | 0.856642 | -23.566 | 58.72 ms |
| Buffer Hanle | dip | -4.670913 | 4.308938 | 0.757571 | -4.394 | 30.08 ms |
| Buffer LCA | crossover | +0.012284 | 1.693135 | 0.629190 | -0.041 | 70.69 ms |
| NMOR | dip | -2.816494 | 0.113283 | 0.827167 | -252.147 | 58.09 ms |

모든 preset guard `≤ 3.0e-12`, fallback 0. exact Maxwell 수치 경로 정상. 다만 “EIT dip”
버튼은 실제 수렴 결과에서 crossover. checklist의 preset 재검토 항목 재현.

### 5.2 B-grid와 NMOR slope

현행 hero slope는 `np.gradient` 중앙 3점 (`gabes/schemes/magneto.py:800,817-823`).

| points | NMOR 중앙 폭 [µT] | 3-point [mrad/µT] | 5-point [mrad/µT] |
|---:|---:|---:|---:|
| 121 | 0.113283 | -252.147 | -285.285 |
| 201 | 0.110669 | -293.582 | -314.868 |
| 401 | 0.109904 | -314.602 | -321.608 |
| 801 | 0.109482 | -320.119 | -321.958 |

기본 hero는 801-point 5-point 값보다 절댓값 **21.7% 작음**. 401-point 5-point는
801-point 5-point 대비 약 0.11%. `magneto-observable-convergence-and-trust` 그대로 P1.
5-point postprocess 비용 O(1). 401-point solve도 현재 수백 ms 이하. 기본 전체 grid를 무조건
키우기보다 local high-order slope + convergence/status 병기가 안전.

### 5.3 paraffin transit 의미 중복

같은 `transit_relax_khz`가:

1. light→dark exchange `γout`,
2. light-region `gamma_gg` pure dephasing

둘 다 됨 (`gabes/schemes/magneto.py:419-424`). exchange kernel도 light block에
`-γout I` 적용 (`gabes/schemes/magneto.py:621-625`).

201-point EIT preset. light `gamma_gg`만 0으로 둔 비저장 진단:

| | 현행 | exchange-only 진단 |
|---|---:|---:|
| feature [m⁻¹] | +0.357285, crossover | -4.428289, dip |
| 중앙 폭 [µT] | 0.980818 | 0.123955 |
| T(B=0) | 0.787789 | 0.857372 |

χ 최대차 / 현행 χ 최대값 약 **35.0%**. “EIT dip” 문제와 직접 연결. 두 rate 의미 분리는
행렬 차원·solve 수 불변, runtime overhead 0. 그러나 어느 rate를 얼마로 둘지는 측정/geometry
근거 필요. 숫자만 제거하면 안 됨.

### 5.4 command zero와 physical zero

Hamiltonian은 `Bphysical = Bcommand + Boffset` 올바르게 사용
(`gabes/schemes/magneto.py:430-434`). hero는 command-axis 0 index 사용
(`gabes/schemes/magneto.py:788-823`).

`b_offset = +0.25 µT`, 401-point EIT:

- command zero, physical `+0.25 µT`: `α=23.906932 m⁻¹`, `T=0.787360`.
- physical zero, command `-0.25 µT`: `α=23.852542 m⁻¹`, `T=0.787789`.

수치차 작아도 의미 다름. physical-zero interpolation과 두 기준 병기는 O(nB), 재 solve 없음.

### 5.5 circular input의 NMOR hero

QWP `45°` drive는 한 circular component가 정확히 0
(`gabes/schemes/magneto.py:107-118`). 그런데 현행 NMOR hero는 여전히

- slope `-32.550 mrad/µT`,
- peak `5.904 mrad`

보고. pure circular input에는 linear polarization plane 없음. 최소 수정은 원형/근원형에서
plane-rotation hero를 unsupported로 gate. O(1), atomic solve 불변. arbitrary ellipticity까지
정량화하려면 Jones/Stokes propagation과 detector model 필요. 후자는 별도 검증 과제.

### 5.6 문헌 비교 문구

`python tests/verify_hanle_eit_eia.py` 현재 결과:

- paraffin linear `0.950 mG` 대 paper `0.12 mG`: 약 7.9×.
- paraffin circular `2.305 mG` 대 paper `0.20 mG`: 약 11.5×.
- buffer low-power `3.349 mG` 대 reference `2.4 mG`: 약 39.5% 큼.

스크립트는 sign switch를 `MATCH`, 앞 두 폭을 `same sub-mG order`로 기술
(`tests/verify_hanle_eit_eia.py:243-271`). sign/trend와 absolute width 판정을 분리해야 함.
추가 계산비용 0.

## 6. 기존 개선안: 비용과 물리 보존

| 기존 제안 | 비용 | 물리 보존 판단 | 권고 |
|---|---:|---|---|
| exact pole-residue Doppler | 완료. 30–71 ms/preset | 옛 coarse grid보다 정확. 절단 없는 같은 Maxwell 평균 | 유지. guard/fallback 보존 |
| observable convergence/trust | slope O(1), physical zero/status O(nB) | atomic model 불변. numerical/readout 신뢰만 개선 | **가장 싼 P1** |
| EIT-dip preset 재검토 | relabel 0; retune은 offline sweep | transit 의미 먼저 결정해야 함 | transit 분리 뒤 retune/relabel |
| sourced pressure shift/coefficient | lookup, scalar shift, 같은 dissipator 차원. runtime 미미 | source/units/sign/uncertainty 명시 시 compact physics 개선 | P1. 구현 연구량은 large |
| geometry-aware relaxation | geometry→scalar rate 자체 미미 | geometry/diffusion/wall 자료 없으면 false precision | decision 유지 |
| full VCC | velocity classes 결합. 시간·메모리 큼 | separable exact-Maxwell 경로 깨짐 | opt-in research만 |
| transit exchange/dephasing 분리 | solve overhead 0 | 큰 물리 영향. 별도 rate 의미와 데이터 필요 | 독립 acceptance로 승격 |
| circular NMOR gate | O(1) | 잘못된 observable claim만 제거 | 즉시 가능 |
| full Jones/Stokes propagation | 최소 O(nB), self-consistent면 더 큼 | 검증된 입력/검출 정의 필요 | 별도 P2/reference |

negligible-overhead 최우선은 새 거대 물리 아님. slope/zero/status, transit semantics,
truthful preset, sourced scalar pressure shift.

## 7. 동작 불변 순수 코드 최적화

### 7.1 exact path B-stack 지연 생성

현재 exact Doppler도 먼저 `L0_all`, `Ld_all` 전체 complex B-stack 생성
(`gabes/schemes/magneto.py:479-485`). 정상 exact path는 `C_xy`, `C_z`, `Cd_xy`만 사용.
stack은 guard fallback의 `solve()`에서만 필요 (`gabes/schemes/magneto.py:487-512`).

기본 paraffin F=2→F'=1, 121 points fresh 측정:

- 불필요 stack 두 개: **15.125 MiB**.
- broadcast 조립 median **9.534 ms**, p95 16.43 ms.
- built-in 5 regime fallback 모두 0.

`solve(rows, deff)` 안에서 requested rows만 같은 식으로 생성하면 exact normal path는 allocation 0.
Doppler-off는 전 row, forced fallback은 failed row만 생성. 결과 불변. exact/off/forced-fallback parity
test 필요. 현재 최고 순수 최적화.

### 7.2 angular-momentum matrix cache

`_hamiltonian()`이 호출마다 `angular_momentum_matrices(Fg/Fe)` 재생성
(`gabes/schemes/magneto.py:555-574`; `gabes/zeeman.py:48-66`). fresh pair median `38.2 µs`.
`F`별 immutable cache면 동작 불변. 전체 30–71 ms compute에 비해 작음. 낮은 우선순위.

### 7.3 이미 끝난 최적화는 재제안 금지

- headless observables 완료 (`docs/checklist.json:2774-2803`).
- exact path는 real-basis affine generator와 세 linear readout만 사용
  (`gabes/schemes/magneto.py:681-747`).
- `zeeman_manifold()`와 Hermitian basis 이미 cache (`gabes/zeeman.py:69-70`,
  `gabes/core.py:60-93`).

과거 dense-grid weighted-coherence contraction은 이제 rare fallback/off 경로만 대상. 우선순위 낮음.

## 8. 문서, 테스트, 예제

### 테스트

강점:

- `tests/test_magneto.py:28-526`: CG, emission normalization, TOC, elastic dephasing,
  relaxation 분리, positivity, conditional normalization, density-once, offset, feature sign,
  LCA, NMOR zero crossing.
- `tests/test_magneto_doppler.py:54-114`: exact/reference, guard/fallback, off-path,
  provenance, velocity knob 제거.

빈 곳:

- exact Doppler 후 preset label truth.
- 5-point/local-fit slope와 B-grid convergence.
- physical-zero hero.
- circular-input NMOR gate.
- open-hyperfine/repump scope status.
- literature absolute-width tolerance.
- generated example의 source/config/model hash 일치.

### 사용자 문서

- User Guide runtime은 `~0.35 s`, 설명은 `B scan × velocity`
  (`docs/Userguide/GABES_User_Guide_v2.html:778-784`). 현재 exact path는 30–71 ms,
  normal path velocity-grid 없음. 낡음.
- Magneto 한계는 “실험별 조정 필요”만 표기
  (`docs/Userguide/GABES_User_Guide_v2.html:852-854`). open hyperfine, absolute scale,
  circular NMOR, collision 한계 더 명시 필요.
- “절대 스케일, 선폭, NMOR 영점 교차 자동 확인” 문구는 실제 보장보다 강함
  (`docs/Userguide/GABES_User_Guide_v2.html:935-937`).

### 분석 예제

실질 예제: `analysis/squeezing/resonant_hanle_squeezing_reference.py`.

- 제거된 `velocity_classes`를 계속 config와 params에 넣음
  (`analysis/squeezing/resonant_hanle_experiment_config.example.json:15-18`;
  `analysis/squeezing/resonant_hanle_squeezing_reference.py:370-397`). 현재 무시됨.
- transmission slope도 3-point `np.gradient`
  (`analysis/squeezing/resonant_hanle_squeezing_reference.py:401-413`).
- 저장 Markdown: B=-0.04 µT에서 `T=0.9962`, `|dT/dB|=0.002533/µT`
  (`analysis/squeezing/resonant_hanle_squeezing_reference.md:20-25`).
- 오늘 같은 Hanle 입력 직접 실행: `T=0.642560`, `|dT/dB|=0.0679216/µT`.
- 저장 보고서는 measured Hanle CSV 없음도 명시하면서 구체 sensitivity 제시
  (`analysis/squeezing/resonant_hanle_squeezing_reference.md:41-60`).
- 생성 시각, git dirty/commit, Magneto cache version, config/source hash 없음
  (`analysis/squeezing/resonant_hanle_squeezing_reference.py:1020-1171`).

따라서 저장 example은 현행 experimental reference 아님. 재생성 전 sensitivity 수치 인용 금지.
provenance 추가 비용 O(file size), OBE 대비 무시 가능.

## 9. 권고 순서

1. **P1 trust:** 5-point/local odd-fit slope, B-grid convergence status, command/physical zero 분리,
   literature trend/absolute 판정 분리.
2. **P1 semantics:** paraffin transit exchange와 light pure dephasing 분리. 독립 knob/acceptance/test.
3. **P1 truthful preset:** 2번 뒤 EIT preset을 문헌 조건으로 retune하거나 crossover로 relabel.
4. **P1 observable scope:** circular/near-circular NMOR plane-rotation hero unsupported gate.
5. **P1 reproducibility:** stale Hanle–FWM Markdown/User Guide 재생성. source/config/git/model hash 포함.
6. **P1 collision data:** gas/species/line별 sourced broadening/shift table.
7. **P1 pure speed:** exact path `L0_all/Ld_all` lazy rows.
8. **보류:** full VCC, geometry-derived relaxation, open-hyperfine/repump, full polarization
   propagation은 dataset/runtime budget 있는 opt-in reference.

## 10. 검증 기록

- registry와 현재 순서 직접 실행 확인.
- 5 built-in regime exact-Maxwell 재실행, warm timing, guard/fallback 확인.
- 121/201/401/801-point central width와 3/5-point slope 확인.
- paraffin light dephasing 제거 비저장 진단.
- command zero / physical zero 비교.
- QWP 0/45° NMOR readout 비교.
- current Hanle example input 직접 실행.
- `python tests/verify_hanle_eit_eia.py` 실행.
- focused pytest:
  `python -m pytest -q tests/test_magneto.py tests/test_magneto_doppler.py tests/test_kernels.py tests/test_headless_observables.py tests/test_schemes_render.py tests/test_resonant_hanle_reference.py tests/test_docs_consistency.py`
  → **81 passed, 1 failed in 16.74 s**.
- focused 유일 실패: 기존 working tree에서 삭제된
  `docs/FWM physics and analytic reconstruction/FWM_physics.tex`를 docs consistency test가 찾음.
  Scheme 4 실패 아님.
- AGENTS.md 필수 root `python -m pytest -q` 실행. **31.86 s 뒤 collection 중단**.
  `references/AutoOD/rb_od_calculator_optimized_final.py`가 optional `sympy`를 import하나 미설치.
- 기존 사용자 변경 보존. production 코드 수정 없음. 오늘 보고서만 추가.

## 11. 최종 판정

Scheme 4는 **유용한 실제 원자물리** 구현. Zeeman CG, optical pumping, grouped-emission TOC,
transverse-field switching, coated-cell two-region coherence, buffer relaxation, exact Maxwell 평균,
Beer–Lambert/NMOR readout 모두 실험 설계와 초기 trace 해석에 가치 있음.

그러나 현행 “EIT dip”은 crossover. 기본 NMOR slope는 21.7% 미수렴. paraffin transit rate는
exchange와 pure dephasing에 중복 적용. physical zero와 hero zero 다름. pure circular에도
plane-rotation 수치 출력. open hyperfine/repump, collision provenance, 공간 평균, absolute scale
부족. 생성 sensitivity example도 현행 코드와 불일치.

따라서 등급: **정성·반정량 experimental atomic-physics reference / planning tool**.
**absolute magnetometry, publication-grade linewidth/rotation/sensitivity reference로 부적합**.
최고 비용효율 작업은 큰 새 모델 아님. transit 의미, slope/zero/status, observable gate,
provenance 먼저.
