# 2026-09-28 Scheme 4 물리 검토 — Hanle / EIA / NMOR

## 1. 선택과 현재 다섯 scheme

- 로컬 날짜: **2026-09-28**, Asia/Seoul.
- 계산: **n = (28 mod 5) + 1 = 4**.
- 대상: **Scheme 4, Magneto-optics (Hanle / EIA / NMOR)**.
- 현재 HEAD: **`5b29c4f`**. 기준은 현재 working tree. 기존 staged/unstaged/untracked 변경 보존.

실행 registry 순서:

| n | key | 현재 정의 |
|---:|---|---|
| 1 | `sas` | Absorption spectroscopy: pump-off OD + pump-on hyperfine-pumping SAS |
| 2 | `lambda` | 3-level Lambda coherence: EIT / AT / CPT |
| 3 | `rydberg_eit` | cascade Rydberg-EIT + microwave AT electrometry |
| 4 | `magneto` | ⁸⁷Rb D1 Zeeman Hanle / EIA / NMOR |
| 5 | `fwm` | seeded ⁸⁵Rb D1 double-Lambda mean-field + generic SFWM biphoton |

근거: `gabes/schemes/__init__.py:12-24`, `README.md:8-16`. `all_schemes()`도 이 순서 반환.

## 2. 먼저 검색한 기존 제안, TODO, issue

검색 범위:

- Scheme 4 기존 일일 보고서 **13개**. 최신은
  `docs/daily_report/2026-09-23_scheme-4_hanle-eia-nmor.md`.
- `docs/checklist.json`의 Magneto, collision, paraffin, Doppler, VCC 항목.
- `gabes`, `tests`, `analysis`, `docs`의 `TODO`, `FIXME`, `issue`, `proposal`,
  `Hanle`, `EIA`, `NMOR`.
- 독립 TODO/issue note 파일과 공개 GitHub Issues API.

결과:

- 별도 로컬 issue/TODO 파일 없음.
- 공개 `Shake2313/fwm-squeezing-app` issue: open/closed 합계 **0**.
- 직접 관련 코드 TODO 하나: 고정 Ne broadening을 gas/species/line 계수표와 pressure
  shift로 승격 (`gabes/constants.py:79-89`).
- **기존 개선안 많음. 새 거대 모델을 먼저 제안할 상황 아님.**

현재 핵심 항목:

| 항목 | 상태 | 위치 |
|---|---|---|
| collision coefficient provenance + pressure shift | P1, ready, large | `docs/checklist.json:2630-2663` |
| Magneto dynamics semantics | done | `docs/checklist.json:2805-2823` |
| light-region normalization | done | `docs/checklist.json:2826-2845` |
| observable convergence / physical zero / trust | P1, ready, medium | `docs/checklist.json:2848-2867` |
| geometry-aware buffer relaxation | P2, needs_decision, large | `docs/checklist.json:2870-2899` |
| full velocity-changing collision | P2, parked, research | `docs/checklist.json:3035-3053` |
| exact pole-residue Doppler | done | `docs/checklist.json:3263-3279` |
| exact-Doppler 뒤 EIT-dip preset 재검토 | P2, needs_decision, small | `docs/checklist.json:3299-3317` |
| paraffin transport/coating calibration | P2, parked, research | `docs/checklist.json:2728-2751` |

과거 Scheme 4 보고서에만 남은 구체 제안:

- paraffin transit exchange와 in-beam pure dephasing 분리.
- circular/near-circular 입력에서 polarization-plane-rotation hero gate.
- Hanle–FWM 생성 예제 provenance와 현행 코드 재생성.

`4cb011f..5b29c4f`의 Scheme 4 관련 경로에는 commit 변화 없음. 현재
`gabes/schemes/magneto.py` SHA-256도 9월 18일 checkpoint의
`3FF0032E...`와 같음. 이번 검토는 새 구현 평가가 아니라 현행 결론 재검증.

## 3. 현재 구현하는 실제 물리

### 3.1 원자와 광장

- ⁸⁷Rb D1. `Fg, Fe ∈ {1,2}` 네 hyperfine pair. Wigner line strength, transition
  dipole, low-field hyperfine Landé `gF` 사용
  (`gabes/schemes/magneto.py:21-30,55-104,384-400`).
- 선택한 pair의 전 `mF` 상태 사용. σ⁺/σ⁻/π CG 결합과 CG-branched spontaneous
  emission 생성 (`gabes/zeeman.py:69-123`).
- spontaneous emission을 편광별 jump operator로 묶어 excited coherence에서 ground
  coherence로 transfer of coherence 보존. EIA/LCA에 필요한 실제 물리
  (`gabes/zeeman.py:101-108`; `tests/test_magneto.py:37-63,458-490`).
- QWP angle을 complex circular-basis drive로 변환. longitudinal scan, B-zero offset,
  residual transverse field를 Zeeman Hamiltonian에 함께 적용
  (`gabes/schemes/magneto.py:107-118,208-219,430-450,555-574`).

### 3.2 셀과 완화

- buffer cell: natural emission/TOC는 `GAMMA_D1`. Ne broadening은 elastic optical
  dephasing. reservoir reload와 collisional ground dephasing 별도
  (`gabes/schemes/magneto.py:401-417`; `gabes/zeeman.py:125-149`).
- paraffin cell: illuminated/dark 두 영역 교환. dark 영역 wall-coherence lifetime.
  optical state는 trace-one conditional state, light-region occupation은 진단값
  (`gabes/schemes/magneto.py:418-428,576-635`).
- local vapor density를 susceptibility에 한 번만 적용
  (`gabes/schemes/magneto.py:530-552`; `gabes/observables.py:556-570`).

### 3.3 Doppler와 관측량

- Doppler-on: 각 B row를 real Hermitian-basis affine shift-invert system으로 변환.
  pole-residue + Faddeeva로 **절단 없는 Maxwell 평균**. direct `v=0` guard 실패 row만
  ±6σ 수렴 grid fallback (`gabes/schemes/magneto.py:45-52,496-516,663-747`).
- 테스트가 exact 평균과 독립 fallback grid를 비교하고 method provenance, off path를 고정
  (`tests/test_magneto_doppler.py:54-114`).
- 관측량: `alpha = k Im(chi)`, `T = exp(-alpha L)`,
  `theta = kL Re(chi+ - chi-)/4` (`gabes/schemes/magneto.py:788-823`).

따라서 단순 Lorentzian 장난감 아님. selection rule, optical pumping, TOC, Zeeman
precession, transverse-field switching, two-region Ramsey core, Beer–Lambert/NMOR
readout을 연결.

## 4. 실험 원자물리 reference 판정

### 유용한 범위

- Hanle dip/peak/crossover 부호와 topology sanity check.
- linear/circular/elliptical drive, residual transverse field, longitudinal offset 추세.
- wall coherence, dark return, buffer reload/dephasing의 상대 민감도.
- detuning, intensity, cell length, B range의 초기 실험 설계.
- 측정 trace의 반정량 fit seed.

코드도 Lee–Moon paraffin MIA, Moon–Kim Ramsey EIA/MIT, Yu–Lee–Moon buffer LCA를
명시 (`gabes/schemes/magneto.py:362-381`).

### 절대 reference로 부적합한 이유

1. 선택한 hyperfine pair만 닫힌 모델. 실제 D1 open transition의 다른 ground F 누출,
   repump, 전체 excited-hyperfine interference 없음 (`gabes/zeeman.py:69-130`).
2. Ne collision은 `3.91 MHz/Torr` broadening scalar. pressure shift, gas/line/T 의존,
   diffusion, spin destruction coefficient provenance, Dicke narrowing, VCC 없음
   (`gabes/constants.py:79-89`).
3. Gaussian beam/spatial intensity, magnetic-field gradient, cell geometry, diffusion 평균 없음.
4. NMOR은 thin/weak-medium circular birefringence 식. arbitrary ellipticity의 Jones/Stokes
   propagation 아님. pure circular input에는 polarization plane 자체가 없음
   (`gabes/schemes/magneto.py:107-118,788-823`).
5. `line_strength`가 0.01–2.0 calibration knob. 절대 contrast/rotation/slope는 코드도
   `external validation required` 표시 (`gabes/schemes/magneto.py:41-43,280-284,825-857`).
6. absolute sensitivity에는 measured optical power, detector noise, transfer calibration 필요.
   production scheme은 이를 제공하지 않음.

판정: **실제 물리를 가진 정성·반정량 experimental-planning/reference**.
**absolute magnetometry, publication-grade linewidth/rotation/sensitivity reference 아님.**

## 5. 2026-09-28 현행 수치 재검증

### 5.1 built-in 121-point exact-Maxwell

동일 process warm 7회 median. 모든 guard `<= 3.0e-12`, grid fallback 0.

| regime | 현행 feature | `alpha(0)-wing` [m⁻¹] | central width [µT] | T(command B=0) | median |
|---|---|---:|---:|---:|---:|
| EIT dip | **crossover** | +0.357285 | 0.981271 | 0.787789 | 59.05 ms |
| EIA peak | peak | +14.930043 | 0.363514 | 0.856642 | 57.14 ms |
| Buffer Hanle | dip | -4.670913 | 4.308938 | 0.757571 | 30.01 ms |
| Buffer LCA | crossover | +0.012284 | 1.693135 | 0.629190 | 68.82 ms |
| NMOR | dip | -2.816494 | 0.113283 | 0.827167 | 58.24 ms |

exact Maxwell 수치 경로 정상. 그러나 “EIT dip” 버튼은 이름과 다른 crossover.
`magneto-eit-dip-preset-exact-doppler` 그대로 재현.

### 5.2 B-grid와 NMOR slope

현행 hero는 중앙 3-point `np.gradient` (`gabes/schemes/magneto.py:800,817-823`).

| B points | central width [µT] | 3-point slope [mrad/µT] | 5-point slope [mrad/µT] | warm median |
|---:|---:|---:|---:|---:|
| 121 | 0.113283 | -252.147 | -285.285 | 49.44 ms |
| 201 | 0.110669 | -293.582 | -314.868 | 80.33 ms |
| 401 | 0.109904 | -314.602 | -321.608 | 144.37 ms |
| 801 | 0.109482 | -320.119 | -321.958 | 281.26 ms |

기본 hero는 801-point 5-point보다 절댓값 **21.7% 작음**. 401-point 5-point는
801-point 5-point와 0.11% 차이. 전체 기본 grid를 무조건 키우면 약 3–6배 비용.
local high-order stencil/odd fit + convergence status가 더 효율적.

### 5.3 paraffin transit 의미 중복

`transit_relax_khz` 하나가 light→dark exchange `gamma_out`과 light-region
`gamma_gg` pure dephasing 둘 다 됨 (`gabes/schemes/magneto.py:419-424`). exchange
kernel도 light block에 `-gamma_out I` 적용 (`gabes/schemes/magneto.py:621-625`).

201-point EIT preset. light `gamma_gg`만 0인 비저장 진단:

| | 현행 | exchange-only 진단 |
|---|---:|---:|
| feature | crossover | dip |
| `alpha(0)-wing` [m⁻¹] | +0.357285 | -4.428289 |
| central width [µT] | 0.980818 | 0.123955 |
| T(B=0) | 0.787789 | 0.857372 |

`chi_probe` 최대차 / 현행 최대값 **35.0%**. rate 의미 분리는 matrix 차원과 solve 수를
바꾸지 않아 runtime overhead 0. 그러나 별도 rate 값은 geometry/measurement 근거 필요.
무근거로 현행 dephasing을 삭제하면 안 됨.

### 5.4 command zero, physical zero

Hamiltonian은 `Bphysical = Bcommand + Boffset` 사용 (`gabes/schemes/magneto.py:430-434`).
hero는 command-axis 0 index 사용 (`gabes/schemes/magneto.py:788-823`).

`b_offset=+0.25 µT`, 401-point EIT:

- command zero = physical `+0.25 µT`: `alpha=23.906932 m⁻¹`, `T=0.787360`.
- physical zero = command `-0.25 µT`: `alpha=23.852542 m⁻¹`, `T=0.787789`.

차이는 작아도 의미 다름. physical-zero interpolation과 두 기준 표시는 O(nB), re-solve 없음.

### 5.5 pure circular의 NMOR hero

QWP 45°는 한 circular component가 정확히 0 (`gabes/schemes/magneto.py:107-118`).
401-point 현행 NMOR은 그래도 slope `-41.37 mrad/µT`, peak `5.944 mrad`를
polarization rotation으로 출력. pure circular에는 linear polarization plane 없음.
최소 gate는 O(1), atomic solve 불변. full Jones/Stokes propagation은 별도 검증 과제.

### 5.6 문헌 비교 스크립트

`python tests/verify_hanle_eit_eia.py` 현행 결과:

- paraffin linear: `0.950 mG` 대 paper `0.12 mG` → 약 **7.9배**.
- paraffin circular: `2.305 mG` 대 paper `0.20 mG` → 약 **11.5배**.
- buffer low-power: `3.349 mG` 대 reference `2.4 mG` → **39.5% 큼**.

그런데 스크립트는 sign switch를 `MATCH`, 앞 두 폭을 `same sub-mG order`로 표현
(`tests/verify_hanle_eit_eia.py:243-271`). sign/trend PASS와 absolute-width CHECK/FAIL을
분리해야 함. 추가 계산비용 0.

## 6. 기존 개선안의 비용과 물리 보존성

| 기존 제안 | 계산비용 | 물리 보존 판단 | 권고 |
|---|---:|---|---|
| exact pole-residue Maxwell | 완료. 30–69 ms/preset | 옛 coarse grid보다 정확. 같은 Maxwell 물리 | guard/fallback 유지 |
| observable convergence/trust | postprocess O(nB); 401 full scan은 144 ms | atomic model 불변. 수치/readout 신뢰 개선 | **최우선 P1** |
| physical-zero hero 분리 | O(nB), re-solve 0 | 의미 오류만 제거 | 즉시 가능 |
| literature trust label | 0 | 과장 문구만 제거 | 즉시 가능 |
| EIT preset 재검토 | relabel 0; retune은 offline sweep | transit semantics 결정 전 숫자 tuning 위험 | transit 분리 뒤 retune/relabel |
| sourced collision coefficients/pressure shift | lookup + scalar detuning shift; solve 차원 동일 | source/units/sign/uncertainty 있으면 compact physics 개선 | runtime 미미, 조사/API 작업 large |
| transit exchange/dephasing 분리 | solve overhead 0 | 물리 영향 큼. 서로 다른 과정 | 별도 knob/acceptance/test 필요 |
| circular NMOR gate | O(1) | 정의 안 되는 observable claim만 차단 | 즉시 가능 |
| geometry-aware relaxation | geometry→rate 계산 자체 작음 | 데이터 없으면 false precision | needs_decision 유지 |
| calibrated coated-cell transport | bright/dark distribution·spatial model이면 큼 | 실제 Ramsey transport 개선 | dataset 있는 opt-in research |
| full VCC | velocity classes 결합. 시간·메모리 큼 | separable exact-Maxwell 경로 깨짐 | parked reference solver만 |

negligible-overhead 우선순위: slope/zero/status, transit 의미 분리, truthful preset,
sourced scalar pressure shift. full VCC나 geometry transport를 interactive default로 넣지 말 것.

## 7. 동작 불변 순수 코드 최적화

### 7.1 exact path B-stack 지연 생성 — 최우선

현행 exact Doppler도 먼저 `L0_all`, `Ld_all` 전체 complex B-stack 생성
(`gabes/schemes/magneto.py:479-485`). 정상 exact path는 `C_xy`, `C_z`, `Cd_xy`만 사용.
stack은 guard fallback의 `solve()`에서만 필요 (`gabes/schemes/magneto.py:487-512`).

현재 121-point paraffin F=2→F'=1 측정:

- 불필요 stack 두 개: **15.125 MiB**.
- 조립 median **9.95 ms**, p95 **10.56 ms**.
- built-in 5 regime fallback 모두 0.

`solve(rows,deff)` 안에서 요청 row만 동일 식으로 생성하면 결과/물리 불변. exact 정상 경로는
allocation 0. Doppler-off는 전 row, fallback은 failed row만 생성. exact/off/forced-fallback
parity test 필요.

### 7.2 angular-momentum matrix cache — 낮은 우선순위

`_hamiltonian()`이 같은 Fg/Fe의 `Fx,Fy,Fz`를 반복 생성
(`gabes/schemes/magneto.py:555-574`; `gabes/zeeman.py:48-66`). immutable `F`별 cache면
결과 불변. 선행 측정 약 38 µs/pair. 전체 30–69 ms 대비 작음.

### 7.3 BLAS controller 재사용 — 낮은 우선순위

`core.blas_single_thread()`는 매번 새 `threadpool_limits` context 생성
(`gabes/core.py:17-34`). 선행 측정 steady 0.406 ms 대 cached controller 0.002 ms.
late-loaded BLAS와 thread-safety 검증 필요. 7.1보다 이득 작음.

이미 끝난 것 재제안 금지:

- headless observables 완료 (`docs/checklist.json:2774-2778`).
- real-basis pole readout 완료 (`gabes/schemes/magneto.py:681-747`).
- `zeeman_manifold()`와 Hermitian basis cache 완료
  (`gabes/zeeman.py:69-70`; `gabes/core.py:61-94`).

## 8. 문서, 테스트, 예제

### 테스트 강점

- `tests/test_magneto.py:28-526`: CG, emission normalization, TOC, elastic dephasing,
  reload/dephasing 분리, positivity, conditional normalization, density-once, offset,
  feature sign, LCA, NMOR zero crossing.
- `tests/test_magneto_doppler.py:54-114`: exact/reference, guard/fallback, off path,
  provenance, obsolete velocity knob 제거.

남은 test gap:

- exact Maxwell 뒤 preset label truth.
- B-grid convergence와 local odd/5-point slope.
- physical-zero hero.
- circular-input NMOR gate.
- open-hyperfine/repump scope status.
- literature absolute-width tolerance.
- generated example의 config/source/model/git hash.

### 사용자 문서

- User Guide Magneto runtime `~0.35 s`, 비용 설명 `B scan × velocity`
  (`docs/Userguide/GABES_User_Guide_v2.html:778-790`). 현재 exact path는 30–69 ms이고
  normal path velocity grid 없음. 낡음.
- Magneto 한계는 “실험별 조정 필요”만 표기
  (`docs/Userguide/GABES_User_Guide_v2.html:852-854`). open hyperfine, absolute scale,
  circular NMOR, collision 한계 누락.
- “절대 스케일, 선폭, NMOR 영점 교차 자동 확인”은 실제 보장보다 강함
  (`docs/Userguide/GABES_User_Guide_v2.html:935-937`).

### 분석 예제

실질 예제: `analysis/squeezing/resonant_hanle_squeezing_reference.py`.

- 삭제된 `velocity_classes`를 config/CLI/params에 계속 넣음
  (`analysis/squeezing/resonant_hanle_experiment_config.example.json:15-18`;
  `analysis/squeezing/resonant_hanle_squeezing_reference.py:370-397,1647-1651`). 현재 무시됨.
- transmission slope도 3-point `np.gradient`
  (`analysis/squeezing/resonant_hanle_squeezing_reference.py:401-413`).
- 저장 Markdown은 `T=0.9962`, `|dT/dB|=0.002533/µT`
  (`analysis/squeezing/resonant_hanle_squeezing_reference.md:20-25`). 9월 23일 현행 직접 실행은
  같은 Hanle 입력에서 `T=0.642560`, `|dT/dB|=0.0679216/µT`.
- measured Hanle CSV 없음도 명시하면서 구체 pT/√Hz 수치 제시
  (`analysis/squeezing/resonant_hanle_squeezing_reference.md:36-60`).
- 생성 시각, git dirty/commit, Magneto cache version, config/source hash 없음.
  report writer도 command dict만 저장 (`analysis/squeezing/resonant_hanle_squeezing_reference.py:1019-1171`).

따라서 저장 example은 현행 experimental reference 아님. 재생성 전 sensitivity 수치 인용 금지.
provenance 추가 비용은 O(file size), OBE 대비 무시 가능.

## 9. 권고 순서

1. **P1 trust:** local converged slope/odd fit, B-grid status, command/physical zero 분리,
   literature trend/absolute 판정 분리.
2. **P1 semantics:** paraffin transit exchange와 light pure dephasing 분리. 독립 rate 이름,
   acceptance, limit test 추가.
3. **P1 truthful preset:** 2번 뒤 EIT preset을 문헌 조건으로 retune하거나 crossover로 relabel.
4. **P1 observable scope:** circular/near-circular NMOR plane-rotation hero unsupported gate.
5. **P1 reproducibility:** stale Hanle–FWM Markdown/User Guide 재생성. source/config/git/model
   hash 포함.
6. **P1 collision data:** gas/species/line별 sourced broadening/shift table.
7. **P1 pure speed:** exact path `L0_all/Ld_all` lazy rows.
8. **보류:** full VCC, geometry-derived relaxation, open-hyperfine/repump, full polarization
   propagation은 dataset/runtime budget 있는 opt-in reference.

## 10. 검증 기록

- registry와 현재 순서 확인.
- 기존 Scheme 4 보고서 13개, checklist, TODO/FIXME, local/public issue 검색.
- five built-in regime exact-Maxwell 재실행, warm timing, guard/fallback 확인.
- 121/201/401/801-point central width, 3/5-point slope와 timing 확인.
- paraffin light dephasing 제거 비저장 진단.
- command zero / physical zero, QWP 45° NMOR 확인.
- B-stack allocation/time microbenchmark.
- `python tests/verify_hanle_eit_eia.py` 실행.
- focused pytest:
  `python -m pytest -q tests/test_magneto.py tests/test_magneto_doppler.py tests/test_kernels.py tests/test_headless_observables.py tests/test_schemes_render.py tests/test_resonant_hanle_reference.py tests/test_docs_consistency.py`
  → **81 passed, 1 failed in 16.71 s**.
- focused 유일 실패: 기존 working tree에서 삭제된
  `docs/FWM physics and analytic reconstruction/FWM_physics.tex`를 docs consistency test가 찾음.
  Scheme 4 실패 아님.
- AGENTS.md 필수 root `python -m pytest -q` 실행. **23.82 s 뒤 collection 중단**.
  `references/AutoOD/rb_od_calculator_optimized_final.py`의 optional `sympy` 미설치.
- production 코드 수정 없음. 오늘 보고서만 추가.

## 11. 최종 판정

Scheme 4는 **유용한 실제 원자물리** 구현. Zeeman CG, optical pumping, grouped-emission TOC,
transverse-field switching, two-region coated-cell coherence, buffer relaxation, exact Maxwell 평균,
Beer–Lambert/NMOR readout은 실험 설계와 초기 trace 해석에 가치 있음.

하지만 현행 “EIT dip”은 crossover. 기본 NMOR slope는 21.7% 미수렴. paraffin transit는
exchange와 pure dephasing에 중복 적용. physical zero와 hero zero 다름. pure circular에도
plane-rotation 출력. open hyperfine/repump, collision provenance, spatial averaging, absolute scale
부족. 생성 sensitivity example도 현행 코드와 불일치.

등급: **정성·반정량 experimental atomic-physics reference / planning tool**.
**absolute magnetometry, publication-grade linewidth/rotation/sensitivity reference로 부적합**.
최고 비용효율 작업: 큰 새 모델 아님. transit 의미, slope/zero/status, observable gate,
provenance, lazy B-stack 먼저.
