# 2026-10-03 Scheme 4 물리 검토 — Hanle / EIA / NMOR

## 1. 오늘 선택과 현재 다섯 scheme

- 로컬 날짜: **2026-10-03**, Asia/Seoul.
- 계산: **n = (3 mod 5) + 1 = 4**.
- 대상: **Scheme 4, Magneto-optics (Hanle/MOR)**. 물리 범위는 Hanle / EIA / NMOR.
- 현재 HEAD: **`fcb312d`**. 판정 기준은 현재 working tree.
- Scheme 4 관련 working-tree 변경 존재. `gabes/schemes/magneto.py`,
  `gabes/pole_doppler.py`, `docs/checklist.json` 수정 상태. exact-Doppler 테스트와
  개발 기록은 untracked. 기존 사용자 변경 보존.

실행 registry 순서:

| n | key | 현재 정의 |
|---:|---|---|
| 1 | `sas` | Absorption spectroscopy (OD / SAS) |
| 2 | `lambda` | Λ coherence (EIT / AT / CPT) |
| 3 | `rydberg_eit` | Rydberg-EIT electrometry |
| 4 | `magneto` | Magneto-optics (Hanle/MOR) |
| 5 | `fwm` | Four-wave mixing (Squeezing / Biphoton) |

근거: registry 인스턴스 순서 `gabes/schemes/__init__.py:19-24`; `all_schemes()`가
그대로 복사 반환 `gabes/schemes/__init__.py:37-39`; UI가 이 순서의 title list 사용
`streamlit_app.py:131-149`. 실제 title 정의는 `gabes/schemes/sas.py:178`,
`gabes/schemes/absorption.py:466-507`, `gabes/schemes/rydberg.py:94-97`,
`gabes/schemes/magneto.py:151-167`, `gabes/schemes/fwm.py:3626-3629`.

## 2. 기존 제안부터 검색

검색 범위:

- 기존 Scheme 4 일일 보고서 **14개**. 최신:
  `docs/daily_report/2026-09-28_scheme-4_hanle-eia-nmor.md`.
- `docs/checklist.json`, `analysis/pole_reuse/DEVLOG.md`, `GATE.md`.
- `gabes`, `tests`, `analysis`, `docs`, `README.md`의 `TODO`, `FIXME`, `issue`,
  `proposal`, `Hanle`, `EIA`, `NMOR`.
- issue/TODO/proposal 이름의 로컬 파일.

결과:

- **기존 개선안 많음.** 새 거대 모델을 먼저 제안할 상황 아님.
- 별도 Scheme 4 로컬 issue-note 파일 없음.
- 직접 관련 코드 TODO 하나: Ne broadening scalar를 gas/species/line 계수표와 pressure
  shift로 승격 `gabes/constants.py:79-89`.
- exact pole-residue Maxwell 이식은 checklist상 완료
  `docs/checklist.json:3472-3488`. 그러나 현재 HEAD에는 없고 working tree에만 있음.
  clean checkout 재현성 없음. 실험 reference provenance 관점에서 먼저 정리 필요.
- 9월 28일 뒤 Scheme 4 생산 코드 hash 변화 없음. `magneto.py` SHA-256 시작값
  `3FF0032E...` 동일. 이후 두 commit은 FWM/Grand Challenge 범위.

현재 핵심 기존 항목:

| 항목 | 상태 | 위치 |
|---|---|---|
| collision coefficient provenance + pressure shift | P1, ready, large | `docs/checklist.json:2838-2876` |
| Magneto dynamics semantics | done | `docs/checklist.json:3013-3032` |
| light-region normalization | done | `docs/checklist.json:3034-3054` |
| observable convergence / physical zero / trust | P1, ready, medium | `docs/checklist.json:3056-3076` |
| geometry-aware buffer relaxation | P2, needs_decision, large | `docs/checklist.json:3078-3108` |
| coated-cell transport 확장 | P2, parked, research | `docs/checklist.json:2936-2960` |
| exact-Doppler 뒤 EIT-dip preset 재검토 | P2, needs_decision, small | `docs/checklist.json:3507-3526` |

과거 보고서에 반복된 아직 미완 항목:

- paraffin light→dark exchange와 in-beam pure dephasing 분리.
- command zero와 physical zero 분리.
- NMOR slope 국소 수렴.
- pure/near-circular 입력에서 polarization-plane rotation hero 차단.
- open-hyperfine/repump scope 표시.
- 문헌 sign/trend와 absolute width 판정 분리.
- Hanle–FWM 예제와 User Guide 재생성·provenance.

## 3. 현재 코드가 구현하는 실제 물리

### 3.1 원자와 광장

- ⁸⁷Rb D1. `Fg, Fe ∈ {1,2}` 네 선택 transition. Wigner line strength,
  transition dipole, low-field hyperfine Landé `gF` 사용
  (`gabes/schemes/magneto.py:21-30,55-104,384-400`).
- 선택 pair의 모든 `mF` 상태. σ⁺/σ⁻/π CG 결합. 편광별 grouped spontaneous-emission
  jump operator가 excited coherence→ground coherence transfer, 즉 EIA/LCA 핵심 TOC 보존
  (`gabes/zeeman.py:69-123`; `tests/test_magneto.py:37-63`).
- QWP angle을 complex circular drive로 변환. longitudinal scan, B offset, residual transverse
  field 모두 Zeeman Hamiltonian에 적용
  (`gabes/schemes/magneto.py:107-118,208-219,430-450,555-574`).

### 3.2 셀과 완화

- buffer cell: natural emission/TOC는 `GAMMA_D1`. Ne broadening은 elastic optical
  dephasing. reservoir reload와 collisional ground dephasing 별도
  (`gabes/schemes/magneto.py:401-417`; `gabes/zeeman.py:125-149`).
- paraffin cell: illuminated/dark 두 영역 교환. dark region에 wall-coherence lifetime.
  susceptibility는 trace-one conditional light state와 local density 한 번 사용
  (`gabes/schemes/magneto.py:418-428,576-635`; `gabes/observables.py:556-570`).
- normalization, density-once, rate-swap, positivity 회귀 존재
  (`tests/test_magneto.py:233-302,352-446`).

### 3.3 Doppler와 관측량

- Doppler-on: 각 B row를 real Hermitian affine system으로 변환. shift-invert pole-residue와
  Faddeeva로 **untruncated Maxwell average**. direct `v=0` guard 실패 row만 수렴 grid fallback
  (`gabes/schemes/magneto.py:45-52,496-516,663-747`;
  `gabes/pole_doppler.py:82-122,253-337`).
- exact/reference grid, forced fallback, Doppler-off, provenance 테스트
  (`tests/test_magneto_doppler.py:54-114`).
- 관측량: `α = k Im χ`, `T = exp(-αL)`,
  `θ = kL Re(χ+ - χ-)/4` (`gabes/schemes/magneto.py:788-823`).

따라서 단순 Lorentzian toy 아님. selection rule, optical pumping, TOC, Zeeman precession,
transverse-field switching, two-region Ramsey core, exact thermal average, Beer–Lambert/NMOR
readout 연결.

## 4. 실험 원자물리 reference 판정

### 유용

- Hanle dip/peak/crossover 부호와 topology sanity check.
- linear/circular/elliptical drive, residual transverse field, B offset 추세.
- wall coherence, dark return, buffer reload/dephasing의 상대 민감도.
- detuning, intensity, cell length, B range 초기 설계.
- 측정 trace 반정량 fit seed.

코드도 Lee–Moon paraffin MIA, Moon–Kim Ramsey EIA/MIT, Yu–Lee–Moon buffer LCA를 명시
`gabes/schemes/magneto.py:362-381`.

### 절대 reference로 부적합

1. 선택한 `Fg↔Fe` pair만 닫힘. 실제 D1 open transition의 다른 ground F 누출, repump,
   전체 excited-hyperfine interference 없음 (`gabes/zeeman.py:69-130`).
2. Ne collision = 고정 `3.91 MHz/Torr` broadening. pressure shift, gas/line/T 의존,
   diffusion, spin-destruction provenance, Dicke narrowing, VCC 없음
   (`gabes/constants.py:79-89`).
3. Gaussian beam/spatial intensity, field gradient, finite cell geometry, diffusion 평균 없음.
4. NMOR 식은 thin/weak-medium circular birefringence. arbitrary ellipticity Jones/Stokes
   propagation 아님. pure circular에는 polarization plane 없음
   (`gabes/schemes/magneto.py:107-118,788-823`).
5. `line_strength` 0.01–2.0는 calibration knob. 절대 contrast/rotation/slope는 코드도
   `external validation required` 표시
   (`gabes/schemes/magneto.py:41-43,280-284,825-857`).
6. absolute sensitivity에 필요한 measured optical power, detector noise, transfer calibration
   없음.
7. 오늘 exact-Doppler 생산 경로와 test가 HEAD 밖 working tree에만 있음. 숫자 provenance 약함.

판정: **실제 물리를 가진 정성·반정량 experimental-planning/reference**.
**absolute magnetometry, publication-grade linewidth/rotation/sensitivity reference 아님.**

## 5. 2026-10-03 수치 재검증

### 5.1 내장 121-point exact-Maxwell

동일 process warm 5회 median. guard 모두 `≤ 3.0e-12`, fallback 0.

| regime | 현행 feature | `α(0)-wing` [m⁻¹] | central width [µT] | T(command B=0) | median |
|---|---|---:|---:|---:|---:|
| EIT dip | **crossover** | +0.357285 | 0.981271 | 0.787789 | 56.12 ms |
| EIA peak | peak | +14.930043 | 0.363515 | 0.856642 | 57.29 ms |
| Buffer Hanle | dip | −4.670913 | 4.308938 | 0.757571 | 29.96 ms |
| Buffer LCA | crossover | +0.012284 | 1.693135 | 0.629190 | 68.18 ms |
| NMOR | dip | −2.816494 | 0.113283 | 0.827167 | 57.36 ms |

exact thermal path 정상. 그러나 “EIT dip” 버튼은 여전히 crossover.
`magneto-eit-dip-preset-exact-doppler` 미결 그대로.

### 5.2 B-grid와 NMOR slope

현행 hero = 중앙 3-point `np.gradient` (`gabes/schemes/magneto.py:800,817-823`).

| B points | central width [µT] | 3-point slope [mrad/µT] | 5-point slope [mrad/µT] | 실행 |
|---:|---:|---:|---:|---:|
| 121 | 0.113283 | −252.147 | −285.285 | 58.56 ms |
| 201 | 0.110669 | −293.582 | −314.868 | 88.64 ms |
| 401 | 0.109904 | −314.602 | −321.608 | 175.39 ms |
| 801 | 0.109482 | −320.119 | −321.958 | 322.52 ms |

기본 hero는 801-point 5-point보다 절댓값 **21.7% 작음**. 401-point 5-point는
801-point 5-point와 0.11% 차이. full grid 상시 증가보다 local higher-order/odd fit와
convergence status가 싸고 물리 보존.

### 5.3 paraffin transit 의미 중복

`transit_relax_khz`가 light→dark exchange `gamma_out`과 light-region `gamma_gg`
pure dephasing 둘 다 됨 (`gabes/schemes/magneto.py:419-424`). exchange kernel도 light block에
`-gamma_out I` 적용 (`gabes/schemes/magneto.py:621-625`).

201-point EIT preset. light `gamma_gg`만 0으로 둔 비저장 진단:

| | 현행 | exchange-only 진단 |
|---|---:|---:|
| feature | crossover | dip |
| `α(0)-wing` [m⁻¹] | +0.357285 | −4.428289 |
| central width [µT] | 0.980818 | 0.123955 |
| T(B=0) | 0.787789 | 0.857372 |

`χ_probe` 최대차 / 현행 최대값 **35.0%**. rate 의미 분리는 matrix 크기와 solve 수를
안 바꿈. runtime overhead 사실상 0. 단, 새 rate 값은 geometry/measurement 근거 필요.

### 5.4 zero와 observable 정의

- `b_offset=+0.25 µT`, 401-point EIT:
  - command zero = physical `+0.25 µT`: `α=23.906932 m⁻¹`, `T=0.787360`.
  - physical zero = command `−0.25 µT`: `α=23.852542 m⁻¹`, `T=0.787789`.
- Hamiltonian은 physical axis 사용 `gabes/schemes/magneto.py:430-434`; hero는 command-axis
  중앙 사용 `gabes/schemes/magneto.py:788-823`. 두 기준 표시와 interpolation = O(nB),
  atomic re-solve 0.
- QWP 45° NMOR, 401-point: slope `−41.37 mrad/µT`, peak `5.944 mrad`를 여전히
  plane rotation으로 출력. pure circular에는 정의 불가. 최소 gate O(1).

### 5.5 문헌 cross-check

`python tests/verify_hanle_eit_eia.py` 오늘 실행:

- paraffin linear: `0.950 mG` 대 paper `0.12 mG` → **7.9배**.
- paraffin circular: `2.305 mG` 대 paper `0.20 mG` → **11.5배**.
- buffer low-power: `3.349 mG` 대 reference `2.4 mG` → **39.5% 큼**.

스크립트는 sign switch를 `MATCH`, 앞 두 폭을 `same sub-mG order`로 출력
`tests/verify_hanle_eit_eia.py:243-271`. 2.305 mG는 sub-mG 아님. trend PASS와
absolute-width CHECK/FAIL 분리 필요. 계산비용 0.

## 6. 기존 개선안: 비용과 물리 보존

| 기존 제안 | 계산비용 | 물리 보존 판단 | 권고 |
|---|---:|---|---|
| exact pole-residue Maxwell | 완료된 working-tree 경로. 30–68 ms/preset | 옛 coarse grid보다 정확. 같은 Maxwell 물리 | 먼저 commit/provenance 고정. guard/fallback 유지 |
| observable convergence/trust | postprocess O(nB); 국소 fit은 re-solve 0 | atomic model 불변. 숫자 신뢰 개선 | **최우선 P1** |
| physical-zero hero 분리 | O(nB), re-solve 0 | 의미 오류만 제거 | 즉시 가능 |
| literature trust label | 0 | 과장 문구만 제거 | 즉시 가능 |
| EIT preset 재검토 | relabel 0; retune은 offline sweep | transit semantics 결정 전 tuning 위험 | transit 분리 뒤 retune/relabel |
| transit exchange/dephasing 분리 | solve 차원·횟수 불변 | 서로 다른 과정. 영향 큼 | 독립 rate·limit test 추가 |
| circular NMOR gate | O(1) | 정의 안 되는 claim 차단 | 즉시 가능 |
| sourced collision table/pressure shift | lookup + scalar detuning shift | source/units/sign/uncertainty 있으면 compact physics 개선 | runtime 무시 가능, 조사/API effort 큼 |
| geometry-aware buffer budget | rate 계산 자체 작음 | 입력/계수 없으면 false precision | needs_decision 유지 |
| calibrated coated transport | bright/dark distribution·spatial quadrature 필요 | 실제 Ramsey transport 개선 | research 비용 큼 |
| full open-hyperfine/repump | 현재 F2→F1 `n=8,M=64`; full D1 약 `n=16,M=256` | 실제 leakage/interference 회복 | dense solve 대략 64×, memory 16×. opt-in reference |
| full VCC | velocity classes 상호 결합 | exact separable pole 경로 깨짐 | interactive default 부적합 |

negligible-overhead 우선: slope/zero/status, transit 의미 분리, truthful label,
sourced scalar pressure shift. full VCC/transport/full-D1을 기본 경로에 바로 넣지 말 것.

## 7. 동작 불변 순수 코드 최적화

### 7.1 exact path B-stack 지연 생성 — 최우선

현행 exact path도 `L0_all`, `Ld_all` complex B-stack을 먼저 생성
`gabes/schemes/magneto.py:479-485`. 정상 exact path는 `C_xy`, `C_z`, `Cd_xy`만 사용.
stack은 Doppler-off 또는 guard fallback의 `solve()`에서만 필요
`gabes/schemes/magneto.py:487-512`.

오늘 121-point paraffin F=2→F'=1:

- 불필요 stack 두 개: **15.125 MiB**.
- 조립 median **9.51 ms**, p95 **10.77 ms**.
- 요청 row 식으로 만든 배열과 eager stack slice: `np.array_equal=True`.

`solve(rows,deff)` 안에서 요청 row만 같은 broadcast 식으로 생성하면 물리·부동소수 연산식
불변. exact 정상 path allocation 0. off는 전 row, fallback은 failed row만 생성.

### 7.2 BLAS controller 재사용 — 중간 우선순위

`core.blas_single_thread()`가 매번 새 threadpool query/context 생성
`gabes/core.py:17-34`.

- 빈 context: 현재 `0.672 ms`, cached `ThreadpoolController` `0.0049 ms`.
- EIT default median: `50.81 → 48.10 ms`, **1.056×**.
- `chi_probe`, `chi_p`, `chi_m` bitwise 동일.

late-loaded BLAS, nested/concurrent context, thread safety 회귀 뒤 적용. 물리 불변.

### 7.3 angular-momentum matrix cache — 낮은 우선순위

`angular_momentum_matrices(F)`는 매 호출 배열 재생성 `gabes/zeeman.py:48-66`.
`_hamiltonian()`이 같은 Fg/Fe에 반복 호출 `gabes/schemes/magneto.py:555-574`.
오늘 F=2/F=1 pair: `35.74 µs → 0.18 µs`. 전체 30–68 ms 대비 작음.
read-only cached matrices + mutation test면 안전.

이미 완료된 것을 재제안하지 않음:

- `zeeman_manifold()` cache `gabes/zeeman.py:69-70`.
- Hermitian basis cache `gabes/core.py:61-94`.
- headless observables `docs/checklist.json:2982-2987`.
- exact real-basis pole readout `gabes/schemes/magneto.py:681-747`.

## 8. 문서, 테스트, 예제

### 테스트 강점

- `tests/test_magneto.py`: CG, emission normalization, TOC, elastic dephasing,
  reload/dephasing 분리, positivity, conditional normalization, density-once, offset,
  sign, LCA, NMOR zero crossing.
- `tests/test_magneto_doppler.py:54-114`: exact/reference, forced fallback, off path,
  provenance, obsolete velocity knob 제거.

남은 gap:

- exact Maxwell 뒤 preset label truth.
- B-grid convergence와 local 5-point/odd-fit slope.
- physical-zero hero.
- circular-input NMOR gate.
- open-hyperfine/repump scope status.
- literature absolute-width tolerance.
- generated example source/config/model/git hash.

낡은 테스트/스크립트 입력도 남음. `velocity_classes`는 schema에서 제거됐지만
`tests/test_magneto.py:166-169,412-421`, `tests/verify_hanle_eit_eia.py:71-81,98-107,157-187`가
계속 전달. Python dict라 조용히 무시. 테스트가 실제 knob를 검증한다고 오해 가능.

### User Guide

- Magneto `~0.35 s`, 비용 `B scan × velocity` 표기
  `docs/Userguide/GABES_User_Guide_v2.html:778-790`. 현행 normal path 30–68 ms,
  velocity grid 없음.
- Magneto 한계는 “실험별 조정”뿐 `...v2.html:852-854`. open hyperfine,
  absolute scale, circular NMOR, collision 한계 누락.
- “절대 스케일, 선폭, NMOR 영점 교차 자동 확인” 문구
  `...v2.html:935-937`는 실제 검증보다 강함.

### 분석 예제

실질 예제: `analysis/squeezing/resonant_hanle_squeezing_reference.py`.

- 삭제된 `velocity_classes`를 config/CLI/params에 계속 저장·전달
  (`analysis/squeezing/resonant_hanle_experiment_config.example.json:15-18`;
  `analysis/squeezing/resonant_hanle_squeezing_reference.py:370-397,1647-1651`).
- transmission slope도 3-point `np.gradient`
  `analysis/squeezing/resonant_hanle_squeezing_reference.py:401-413`.
- 저장 Markdown은 `T=0.9962`, `|dT/dB|=0.002533/µT`
  `analysis/squeezing/resonant_hanle_squeezing_reference.md:20-25`.
  오늘 현행 direct default는 `T=0.642560`, `|dT/dB|=0.0679216/µT`.
- measured Hanle CSV 없음인데 구체 pT/√Hz 출력
  `analysis/squeezing/resonant_hanle_squeezing_reference.md:36-60`.
- writer는 command dict만 보존. 생성시각, commit/dirty, Magneto cache version,
  source/config hash 없음 `...reference.py:1019-1170`.
- unit test는 frequency anchor, shot-noise scaling, CSV fit을 검사하지만 저장 report의
  현행-code 재현성은 검사하지 않음 `tests/test_resonant_hanle_reference.py:14-131`.

저장 예제는 현행 experimental reference 아님. sensitivity 숫자 인용 금지.
provenance 추가 비용 O(file size), OBE 대비 무시 가능.

## 9. 권고 순서

1. **재현성 P0:** exact-Doppler 코드·테스트·DEVLOG를 review/commit. HEAD에서 같은 숫자 재현.
2. **trust P1:** local converged slope/odd fit, B-grid status, physical/command zero 분리,
   literature trend/absolute 판정 분리.
3. **semantics P1:** paraffin exchange와 light pure dephasing 분리. 독립 이름·limit test.
4. **truthful preset P1:** 3번 뒤 EIT preset을 문헌 조건으로 retune 또는 crossover relabel.
5. **observable scope P1:** circular/near-circular NMOR plane-rotation hero gate.
6. **docs/example P1:** User Guide와 Hanle–FWM Markdown 재생성. source/config/git/model hash.
7. **collision P1:** gas/species/line별 sourced broadening/shift table.
8. **pure speed P1:** exact path B-stack lazy rows. 그 뒤 BLAS controller cache.
9. **보류:** full VCC, geometry transport, full open-hyperfine/repump, full polarization
   propagation은 dataset/runtime budget 있는 opt-in reference.

## 10. 검증 기록

- registry, title, UI 사용 순서 확인.
- Scheme 4 기존 보고서 14개, checklist, pole-reuse 기록, TODO/FIXME/local issue-note 검색.
- current source/docs/tests/example 읽음.
- 내장 5 regime 재실행, warm timing, guard/fallback 확인.
- 121/201/401/801 B-grid width, 3/5-point NMOR slope 재측정.
- transit dephasing 제거 비저장 진단, offset, QWP 45° NMOR 확인.
- B-stack, BLAS context, angular-matrix cache microbenchmark.
- `python tests/verify_hanle_eit_eia.py` 실행.
- focused pytest:
  `python -m pytest -q tests/test_magneto.py tests/test_magneto_doppler.py tests/test_kernels.py tests/test_headless_observables.py tests/test_schemes_render.py tests/test_resonant_hanle_reference.py tests/test_docs_consistency.py`
  → **85 passed in 16.39 s**.
- 필수 root `python -m pytest -q` 실행. **25.83 s 뒤 collection 중단**.
  `references/AutoOD/rb_od_calculator_optimized_final.py`의 optional `sympy` 미설치.
  Scheme 4 회귀 실패 아님.
- production 코드 수정 없음. 오늘 보고서만 추가.

## 11. 최종 판정

Scheme 4는 유용한 실제 원자물리 구현. Zeeman CG, optical pumping, grouped-emission TOC,
transverse-field switching, coated-cell two-region coherence, buffer relaxation, exact Maxwell
평균, Beer–Lambert/NMOR readout은 실험 설계와 초기 trace 해석에 가치 있음.

하지만 built-in EIT label 불일치, 21.7% 미수렴 NMOR hero slope, transit 의미 중복,
physical-zero 혼동, pure-circular plane-rotation 출력, 문헌 폭 7.9–11.5배 불일치,
open-hyperfine/repump·collision·spatial 평균·absolute calibration 부족. 저장 sensitivity example도
현행 코드와 불일치. exact-Doppler가 HEAD 밖이라는 재현성 문제 추가.

등급: **정성·반정량 experimental atomic-physics reference / planning tool**.
**absolute magnetometry, publication-grade linewidth/rotation/sensitivity reference로 부적합**.
최고 비용효율 작업: 새 대형 모델 아님. 먼저 provenance, transit 의미, slope/zero/status,
observable gate, lazy B-stack.
