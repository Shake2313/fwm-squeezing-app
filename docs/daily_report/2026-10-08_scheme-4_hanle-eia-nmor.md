# 2026-10-08 Scheme 4 물리 검토 — Hanle / EIA / NMOR

## 1. 선택

- 로컬 날짜: **2026-10-08**, Asia/Seoul.
- 계산: **n = (8 mod 5) + 1 = 4**.
- 대상: **Scheme 4, Magneto-optics (Hanle/MOR)**.
- 검토 범위: 현재 working tree. HEAD `39be4ae`.
- 저장소 dirty. 기존 사용자 변경 보존. 생산 코드 수정 안 함.

현재 실행 순서:

| n | key | 정의 |
|---:|---|---|
| 1 | `sas` | Absorption spectroscopy (OD / SAS) |
| 2 | `lambda` | Λ coherence (EIT / AT / CPT) |
| 3 | `rydberg_eit` | Rydberg-EIT electrometry |
| 4 | `magneto` | Magneto-optics (Hanle/MOR) |
| 5 | `fwm` | Four-wave mixing (Squeezing / Biphoton) |

근거:

- registry 실제 인스턴스 순서 `gabes/schemes/__init__.py:19-24`.
- `all_schemes()`가 이 목록 복사 반환 `gabes/schemes/__init__.py:37-39`.
- UI가 반환 순서 그대로 title 목록과 selectbox 구성 `streamlit_app.py:130-150`.
- Magneto 실제 key/title `gabes/schemes/magneto.py:151-167`.
- README의 사용자-facing 다섯 묶음 `README.md:8-16`.

## 2. 변경 확인

- 직전 Scheme 4 검토: `docs/daily_report/2026-10-03_scheme-4_hanle-eia-nmor.md`.
- 그 뒤 HEAD 변화는 `39be4ae feat(od): improve CSV fitting and Rb ratios` 한 건. Scheme 4 생산 파일 변경 없음.
- 현재 `gabes/schemes/magneto.py` SHA-256:
  `3FF0032E5F80003D651DFE27E0F5FA158D1E66549C8EA5661FD0F3E9E87E0F0E`.
  10월 3일과 같음.
- `gabes/zeeman.py` SHA-256:
  `543AE7E8888B8D6C491D8E5884A698BF639F300C9D624ED2EB66B8FE345923B2`.
- exact-Doppler 생산 변경은 여전히 HEAD 밖 working-tree 수정.
  `tests/test_magneto_doppler.py`, `analysis/pole_reuse/`도 untracked.
- clean HEAD의 Magneto는 아직 `velocity_classes=9` + scalar `doppler_scale` 경로.
  현재 working tree는 exact pole-residue Maxwell 경로
  `gabes/schemes/magneto.py:45-52,496-516,663-747`.
- 따라서 오늘 수치는 현재 작업트리에서 재현. clean checkout 재현성은 아직 없음.

## 3. 기존 제안 먼저 검색

검색:

- 기존 Scheme 4 일일 보고서 **15개**.
- `docs/checklist.json` 전체.
- `analysis/pole_reuse/DEVLOG.md`, `GATE.md`, `gate_report.json`.
- README, User Guide, `gabes`, `tests`, `analysis`, `docs`의
  `TODO`, `FIXME`, `issue`, `proposal`, `Hanle`, `EIA`, `NMOR`.
- issue/TODO/proposal 이름의 로컬 파일.
- GitHub 공개 저장소 open issue API. **0개**.

결과: **기존 개선안 있음. 충분함. 새 대형 물리 roadmap 항목 추가 불필요.**

직접 항목:

| 항목 | 상태 | 위치 |
|---|---|---|
| collision coefficient provenance + pressure shift | P1, ready, large | `docs/checklist.json:2841-2879` |
| coated-cell transport 확장 | P2, parked, research | `docs/checklist.json:2939-2962` |
| Magneto dynamics semantics | done | `docs/checklist.json:3016-3034` |
| light-region normalization | done | `docs/checklist.json:3037-3056` |
| observable convergence / physical zero / trust | P1, ready, medium | `docs/checklist.json:3059-3078` |
| geometry-aware buffer relaxation | P2, needs_decision, large | `docs/checklist.json:3081-3108` |
| full VCC reference | P2, parked, research | `docs/checklist.json:3246-3265` |
| pole-residue Doppler reuse | done | `docs/checklist.json:3474-3490` |
| exact-Doppler 뒤 EIT preset 재검토 | P2, needs_decision, small | `docs/checklist.json:3510-3527` |

코드 TODO도 기존 항목과 같음. Ne를 고정 broadening scalar가 아닌
gas/species/line 계수표와 pressure shift로 승격
`gabes/constants.py:79-89`.

별도 로컬 Scheme 4 issue-note 파일 없음.

## 4. 현재 구현 물리

### 4.1 실제 물리 있음

- 원자: ⁸⁷Rb D1. `Fg, Fe ∈ {1,2}` 중 한 hyperfine pair 선택.
  transition strength, reduced dipole, low-field hyperfine Landé `gF` 사용
  `gabes/schemes/magneto.py:21-30,55-104,384-400`.
- 선택 pair의 모든 `mF` 상태. σ⁺/σ⁻/π CG 결합
  `gabes/zeeman.py:69-123`.
- spontaneous emission은 편광별 grouped jump operator `Σq`.
  excited coherence가 ground coherence로 전달됨. EIA/LCA에 필요한 TOC 살아 있음
  `gabes/zeeman.py:101-123`; 회귀 `tests/test_magneto.py:37-63`.
- QWP 각도를 복소 circular-basis drive로 변환
  `gabes/schemes/magneto.py:107-118`.
- longitudinal B scan, longitudinal offset, residual transverse B, azimuth 모두 Hamiltonian 입력
  `gabes/schemes/magneto.py:430-450,555-574`.
- buffer cell: natural emission/TOC, elastic Ne optical dephasing, isotropic reservoir reload,
  ground depolarization 분리
  `gabes/schemes/magneto.py:401-417`; `gabes/zeeman.py:125-149`.
- paraffin cell: illuminated/dark 두 영역 교환. dark ground coherence에 wall lifetime
  `gabes/schemes/magneto.py:418-428,576-635`.
- two-region light state는 trace-one conditional state. local vapor density 한 번 적용.
  occupation은 진단값
  `gabes/schemes/magneto.py:576-635`; 회귀 `tests/test_magneto.py:233-302`.
- Doppler-on: affine real-Hermitian system → shift-invert pole/residue → Faddeeva
  untruncated Maxwell average. 각 B row direct guard. 실패 row만 수렴 grid fallback
  `gabes/schemes/magneto.py:663-747`; `gabes/pole_doppler.py:253-337`.
- 관측량:
  `α = k Imχ`, `T = exp(-αL)`, `θ = kL Re(χ+−χ−)/4`
  `gabes/schemes/magneto.py:788-823`.

따라서 toy Lorentzian 아님. Zeeman selection rule, optical pumping, TOC,
transverse-field switching, two-region Ramsey core, thermal average, Beer–Lambert/NMOR 연결 존재.

### 4.2 실험 연구자에게 유용

- zero-field feature 부호: dip / peak / crossover sanity check.
- linear↔elliptical↔circular drive 추세.
- residual transverse field, longitudinal offset 영향.
- wall coherence, dark return, buffer reload/dephasing 상대 민감도.
- detuning, intensity, cell length, B scan 범위 초기 설계.
- 측정 trace의 반정량 fit seed.

내장 참고 문헌도 Lee–Moon paraffin MIA, Moon–Kim Ramsey EIA/MIT,
Yu–Lee–Moon buffer LCA 명시 `gabes/schemes/magneto.py:362-381`.

### 4.3 절대 experimental reference 아님

1. 한 `Fg↔Fe` pair만 닫힘. 실제 D1 open transition의 다른 ground F 누출,
   repump, 두 excited hyperfine 간섭 없음 `gabes/zeeman.py:69-130`.
2. Ne 충돌은 고정 `3.91 MHz/Torr` broadening뿐. pressure shift, gas/line/T 의존,
   diffusion, spin destruction provenance, Dicke narrowing, VCC 없음
   `gabes/constants.py:79-89`.
3. Gaussian beam, spatial intensity, field gradient, finite cell geometry,
   diffusion 평균 없음.
4. NMOR는 circular birefringence thin/weak-medium 식. 일반 Jones/Stokes 전파 아님.
   pure circular에는 polarization plane 자체 없음
   `gabes/schemes/magneto.py:107-118,788-823`.
5. `line_strength`는 0.01–2 calibration knob. 절대 contrast/rotation/slope는 코드도
   `external validation required` 표시
   `gabes/schemes/magneto.py:41-43,280-284,825-857`.
6. absolute sensitivity에 필요한 measured optical power, detector transfer/noise,
   calibration/holdout trace 없음.
7. exact-Doppler 생산 경로가 HEAD 밖. clean-checkout provenance 약함.

판정: **실제 물리를 가진 정성·반정량 experimental-planning/reference**.
**absolute magnetometry, publication-grade linewidth/rotation/sensitivity reference 아님.**

## 5. 오늘 수치 재검증

### 5.1 내장 121-point exact Maxwell

동일 process warm 5회 median. guard 최대 `2.97e-12`. fallback 모두 0.

| regime | 현재 feature | `α(0)-wing` [m⁻¹] | central width [µT] | T(command 0) | median |
|---|---|---:|---:|---:|---:|
| EIT dip | **crossover** | +0.357285 | 0.981271 | 0.787789 | 66.21 ms |
| EIA peak | peak | +14.930043 | 0.363514 | 0.856642 | 62.13 ms |
| Buffer Hanle | dip | −4.670913 | 4.308938 | 0.757571 | 32.18 ms |
| Buffer LCA | crossover | +0.012284 | 1.693135 | 0.629190 | 73.46 ms |
| NMOR | absorption dip | −2.816494 | 0.113283 | 0.827167 | 60.91 ms |

핵심: `EIT dip` 버튼은 여전히 crossover. 기존 checklist 항목 정확.

### 5.2 NMOR B-grid 수렴

현재 hero는 중앙 `np.gradient` 한 번
`gabes/schemes/magneto.py:800,817-823`.

| B points | width [µT] | 현행 3-point slope [mrad/µT] | 5-point slope [mrad/µT] | 실행 |
|---:|---:|---:|---:|---:|
| 121 | 0.113283 | −252.147 | −285.285 | 61.62 ms |
| 201 | 0.110669 | −293.582 | −314.868 | 97.12 ms |
| 401 | 0.109904 | −314.602 | −321.608 | 188.89 ms |
| 801 | 0.109482 | −320.119 | −321.958 | 367.52 ms |

- 기본 hero 절댓값은 801-point 5-point 기준보다 **21.7% 작음**.
- 401-point 5-point는 801-point 5-point와 **0.11%** 차이.
- 기존 121점의 5-point만 써도 오차 약 11.4%. acceptance `<5%`에는 부족.
- 해법: 전체 401점을 기본으로 강제하지 말고 physical zero 근처에 별도 촘촘한 5점 추가.
  추가 pole row 4개면 row 수 +3.3% 수준. 실제 end-to-end 비용은 구현 때 측정 필요.

### 5.3 paraffin transit 의미 중복

`transit_relax_khz`가 두 역할:

- light→dark exchange `gamma_out`.
- light-region ground pure dephasing `gamma_gg`.

근거 `gabes/schemes/magneto.py:419-424`; exchange block에도 `-gamma_out I`
`gabes/schemes/magneto.py:621-625`.

201-point EIT preset. light `gamma_gg`만 0으로 둔 비저장 진단:

| | 현행 | exchange-only 진단 |
|---|---:|---:|
| feature | crossover | dip |
| `α(0)-wing` [m⁻¹] | +0.357285 | −4.428289 |
| central width [µT] | 0.980818 | 0.123955 |
| T(B=0) | 0.787789 | 0.857372 |

`χ_probe` 최대 상대 변화 **35.04%**. 의미 중복 영향 큼.
rate 분리는 matrix 크기와 solve 수를 안 바꿈. runtime overhead 사실상 0.
단, 별도 rate 기본값은 beam/cell geometry나 측정값 필요.

### 5.4 zero와 observable 의미

`b_offset = +0.25 µT`, 401-point EIT:

- command zero = physical `+0.25 µT`: `α=23.906932 m⁻¹`, `T=0.787360`.
- physical zero = command `−0.25 µT`: `α=23.852542 m⁻¹`, `T=0.787789`.

Hamiltonian은 physical axis 사용 `gabes/schemes/magneto.py:430-434`.
hero는 command-axis 중앙 사용 `gabes/schemes/magneto.py:788-823`.
두 값 표시와 interpolation은 O(nB), atomic re-solve 0.

QWP 45° NMOR, 401 points:

- slope `−41.369 mrad/µT`.
- peak `5.944 mrad`.

여전히 polarization-plane rotation으로 출력. pure circular에는 정의 불가.
observable gate 비용 O(1).

### 5.5 문헌 검증 스크립트

`python tests/verify_hanle_eit_eia.py`:

- paraffin linear: `0.950 mG` 대 paper `0.12 mG` → **7.9배**.
- paraffin circular: `2.305 mG` 대 paper `0.20 mG` → **11.5배**.
- buffer low-power: `3.349 mG` 대 reference `2.4 mG` → **39.5% 큼**.

스크립트는 sign switch를 `MATCH`, 앞 두 폭을 `same sub-mG order`로 출력
`tests/verify_hanle_eit_eia.py:243-271`.
`2.305 mG`는 sub-mG 아님. trend PASS와 absolute-width CHECK/FAIL 분리 필요.
계산비용 0.

## 6. 기존 개선안 비용 평가

| 기존 제안 | 계산비용 | 물리 보존 / 판단 |
|---|---:|---|
| exact pole-residue Maxwell | 현재 32–73 ms/preset | 옛 9-class보다 정확. 같은 Maxwell 물리. 먼저 review/commit해 provenance 고정 |
| B-grid/slope convergence | 기존 점 5-point O(1), 그러나 불충분. local 4 row 소폭 추가 | atomic model 불변. hero 신뢰 개선. 최우선 |
| physical/command zero 분리 | O(nB), 재풀이 0 | 의미 오류 제거. 즉시 가능 |
| literature trust label | 0 | 과장 제거. 즉시 가능 |
| open-hyperfine/repump scope status | O(1) | 계산 안 한 물리 명시. 즉시 가능 |
| circular NMOR gate | O(1) | 정의 불가 claim 차단. 즉시 가능 |
| EIT preset relabel | 0 | 현재 계산을 정직하게 표시. retune은 transit 의미 분리 뒤 |
| transit exchange/dephasing 분리 | 차원·solve 수 동일 | 서로 다른 과정 분리. 영향 35%. 독립 rate/limit test 필요 |
| sourced collision table + pressure shift | lookup + scalar shift, 사실상 O(1) | compact physics 개선. 조사/API effort 큼 |
| geometry-aware buffer budget | rate 계산은 작음 | 입력·계수 없으면 false precision. needs_decision 유지 |
| calibrated coated transport | dark-time/spatial quadrature 필요 | 실물리 개선. research 비용 큼 |
| full selected-D1 open manifold | 현행 F2→F1 `n=8, M=64`; full D1 약 `n=16, M=256` | dense LU 대략 64× flops, 16× matrix memory. opt-in reference 적합 |
| full VCC | velocity classes 상호 결합 | separable exact-pole 경로 깨짐. interactive default 부적합 |

결론: 물리 보존 + 무시 가능한 overhead 우선순위는
**zero/status/gate → local slope → transit 의미 분리 → sourced scalar collision shift**.

## 7. 동작 불변 순수 코드 최적화

### 7.1 exact path B-stack 지연 생성 — 최우선

현재 exact path도 `L0_all`, `Ld_all` 전체 complex stack 먼저 생성
`gabes/schemes/magneto.py:479-485`.
정상 exact path는 `C_xy`, `C_z`, `Cd_xy`만 사용
`gabes/schemes/magneto.py:509-512,663-747`.
stack은 Doppler-off 또는 guard fallback `solve()`에서만 필요
`gabes/schemes/magneto.py:487-507,517-528`.

오늘 121-point paraffin F2→F'1:

- 불필요 두 stack: **15.125 MiB**.
- 조립 median **9.969 ms**, p95 **20.727 ms**.
- 요청 row 직접 broadcast와 eager slice: `np.array_equal=True`.

`solve(rows,deff)` 안에서 요청 row만 같은 식으로 조립하면 결과 불변.
normal exact path allocation/조립 제거. off는 전 row, fallback은 실패 row만 생성.

### 7.2 BLAS controller 재사용 — 중간 우선순위

`core.blas_single_thread()`가 매번 새 threadpool limit context 생성
`gabes/core.py:17-34`.

오늘 cached `ThreadpoolController` 비저장 prototype:

- 빈 context: `2.7923 → 0.0176 ms`.
- EIT compute median: `52.9884 → 49.8232 ms`, 약 **6.0% 감소**.
- `chi_probe`, `chi_p`, `chi_m` bitwise 동일.

late-loaded BLAS, nested/concurrent context, thread safety 회귀 뒤 적용.

### 7.3 angular-momentum matrix cache — 낮은 우선순위

`angular_momentum_matrices(F)` 매번 배열 생성
`gabes/zeeman.py:48-66`.

- F=2/F=1 pair 호출 median: `39.50 → 0.40 µs`.
- 전체 30–73 ms 대비 작음.
- read-only `lru_cache` + mutation test면 안전.

이미 된 것 재제안 안 함:

- `zeeman_manifold()` cache `gabes/zeeman.py:69-70`.
- Hermitian basis cache `gabes/core.py:61-94`.
- exact real-basis pole readout `gabes/schemes/magneto.py:681-747`.
- figureless observables checklist 완료.

## 8. 테스트·문서·예제

### 테스트 강점

- `tests/test_magneto.py`: CG, emission normalization, TOC, elastic dephasing,
  reload/dephasing 분리, positivity, conditional normalization, density-once,
  B offset, sign, LCA, NMOR zero crossing.
- `tests/test_magneto_doppler.py:54-114`: exact/reference grid, forced fallback,
  off path, method provenance, obsolete velocity knob 제거.
- 오늘 집중 회귀 **85 passed**.

남은 gap:

- exact Maxwell 뒤 preset label truth.
- B-grid/local slope convergence.
- physical-zero hero.
- circular-input NMOR gate.
- open-hyperfine/repump scope status.
- literature absolute-width tolerance.
- generated example source/config/model/git hash.

낡은 입력:

- schema에서 `velocity_classes` 제거됐지만
  `tests/test_magneto.py:166-169,412-421`,
  `tests/verify_hanle_eit_eia.py:71-81,98-107,157-187`,
  `analysis/squeezing/resonant_hanle_experiment_config.example.json:17`,
  `analysis/squeezing/resonant_hanle_squeezing_reference.py:395`가 계속 전달.
- dict extra key라 조용히 무시. 실제 knob 검증처럼 보일 위험.

### User Guide

- Magneto 비용 `~0.35 s`, 구조 `B scan × velocity`
  `docs/Userguide/GABES_User_Guide_v2.html:778-790`.
  현행 정상 경로 32–73 ms, velocity grid 없음.
- 한계는 “실험별 조정” 한 줄뿐
  `docs/Userguide/GABES_User_Guide_v2.html:852-854`.
  open hyperfine, absolute scale, circular NMOR, collision 한계 누락.
- “절대 스케일, 선폭, NMOR 영점 교차 자동 확인”
  `docs/Userguide/GABES_User_Guide_v2.html:935-937`.
  실제 absolute 검증보다 강함.

### Hanle–FWM 예제

`analysis/squeezing/resonant_hanle_squeezing_reference.py`:

- obsolete `velocity_classes` 보존·전달.
- transmission slope도 `np.gradient` 3-point
  `analysis/squeezing/resonant_hanle_squeezing_reference.py:401-413`.
- 저장 Markdown은 measured Hanle CSV 없음을 밝히지만
  `analysis/squeezing/resonant_hanle_squeezing_reference.md:47`,
  구체 pT/√Hz sensitivity 계속 출력.
- 저장 명령은 있으나 생성시각, commit/dirty, Magneto cache version,
  source/config hash 없음 `...reference.py:1019-1170`.
- unit test는 frequency anchor, shot-noise scaling, CSV fit 검사.
  저장 report의 현행-code 재현성은 검사 안 함
  `tests/test_resonant_hanle_reference.py:14-131`.

저장 sensitivity 숫자는 현행 experimental reference로 인용 금지.

## 9. 권고 순서

1. **P0 재현성:** exact-Doppler 생산 코드, test, DEVLOG review/commit.
2. **P1 trust:** physical/command zero 분리, local-converged NMOR slope,
   B-grid status, 문헌 trend/absolute 판정 분리.
3. **P1 semantics:** paraffin exchange와 light pure dephasing 분리.
4. **P1 preset truth:** 3번 뒤 EIT preset 문헌 근거 retune 또는 crossover relabel.
5. **P1 observable scope:** circular/near-circular NMOR plane-rotation gate.
6. **P1 docs/example:** User Guide와 Hanle–FWM report 재생성. provenance hash 추가.
7. **P1 collision:** sourced gas/species/line broadening + pressure shift.
8. **P1 pure speed:** lazy B-stack. 그 뒤 BLAS controller cache.
9. **보류:** full VCC, calibrated spatial transport, full open-D1/repump,
   full Jones/Stokes propagation. dataset + runtime budget 있는 opt-in reference로만.

## 10. 검증 기록

- `CLAUDE.md`, `README.md` 읽음.
- registry/UI 순서 확인.
- 기존 Scheme 4 보고서 15개, checklist, pole-reuse 기록, TODO/FIXME/local issue-note 검색.
- GitHub open issue 0 확인.
- 현재 source, tests, User Guide, Hanle–FWM example 읽음.
- 5 regime warm timing, pole guard/fallback 확인.
- 121/201/401/801 NMOR slope/width 재측정.
- transit dephasing 제거 비저장 진단, offset, QWP 45° NMOR 재측정.
- lazy B-stack, cached BLAS controller, angular-matrix cache 비저장 benchmark.
- `python tests/verify_hanle_eit_eia.py` 실행.
- focused pytest:
  `python -m pytest -q tests/test_magneto.py tests/test_magneto_doppler.py tests/test_kernels.py tests/test_headless_observables.py tests/test_schemes_render.py tests/test_resonant_hanle_reference.py tests/test_docs_consistency.py`
  → **85 passed in 16.93 s**.
- root pytest: `python -m pytest -q`
  → **2202 passed, 4 skipped, 2 failed in 3208.79 s (53:28)**.
  두 실패 모두 Scheme 4 밖:
  - `tests/quantum/test_thermal_campaign_extend.py::test_empty_campaign_extension_passes_both_real_capsules_without_any_native_solve`:
    Windows Grand Challenge fresh-capsule destination `FileNotFoundError`.
  - `tests/test_fwm_excess_noise.py::test_hidden_noise_sliders_retain_values_across_streamlit_reruns`:
    FWM Streamlit fallback slider가 gain-proportional 전환 뒤 이전 constant metric 유지.
- 생산 코드 수정 없음. 이 보고서만 추가.

## 11. 최종 판정

Scheme 4 실제 물리 있음. Zeeman CG, optical pumping, grouped-emission TOC,
transverse-field switching, two-region coated-cell coherence, buffer relaxation,
exact Maxwell average, Beer–Lambert/NMOR readout 유용.

실험 설계, feature 부호, 상대 추세, 초기 fit seed에 좋음.

절대 reference 아님. EIT label 불일치, 기본 NMOR slope 21.7% 미수렴,
transit 의미 중복, physical-zero 혼동, circular plane-rotation 출력,
문헌 폭 7.9–11.5배 불일치, open-hyperfine/repump·collision·spatial 평균·absolute calibration 부족.
exact-Doppler도 HEAD 밖.

등급: **정성·반정량 experimental atomic-physics reference / planning tool**.
**absolute magnetometry, publication-grade linewidth/rotation/sensitivity reference로 부적합**.

최고 비용효율 작업: 새 대형 모델 아님.
먼저 provenance, zero/slope/status, transit 의미, observable gate, lazy B-stack.
