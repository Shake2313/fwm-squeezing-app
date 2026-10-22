# 2026-10-13 Scheme 4 물리 검토 — Hanle / EIA / NMOR

## 1. 선택과 현재 다섯 scheme 순서

- 로컬 날짜: **2026-10-13**, Asia/Seoul.
- 계산: **n = (13 mod 5) + 1 = 4**.
- 대상: **Scheme 4, Magneto-optics (Hanle/MOR)**.
- 검토 기준: 현재 working tree, branch `main`, HEAD `39be4ae`.
- 저장소 dirty. 기존 사용자 변경 보존. 생산 코드 수정 없음.

현재 UI 순서:

| n | key | 현재 정의 |
|---:|---|---|
| 1 | `sas` | Absorption spectroscopy (OD / SAS) |
| 2 | `lambda` | Λ coherence (EIT / AT / CPT) |
| 3 | `rydberg_eit` | Rydberg-EIT electrometry |
| 4 | `magneto` | Magneto-optics (Hanle/MOR) |
| 5 | `fwm` | Four-wave mixing (Squeezing / Biphoton) |

근거:

- 실제 registry 인스턴스 순서: `gabes/schemes/__init__.py:12-25`.
- `all_schemes()`가 이 목록 복사 반환: `gabes/schemes/__init__.py:37-39`.
- UI가 같은 순서의 title 목록과 selectbox 사용: `streamlit_app.py:130-148`.
- 사용자-facing 정의: `README.md:8-16`.
- Scheme 4 key/title/caption: `gabes/schemes/magneto.py:151-167`.

## 2. 재현 범위와 변경 확인

- 직전 Scheme 4 검토: `docs/daily_report/2026-10-08_scheme-4_hanle-eia-nmor.md`.
- 그 뒤 HEAD 변화 없음. HEAD는 계속 `39be4ae`.
- 현재 생산 파일과 테스트 hash도 직전 검토와 같음.

| 파일 | SHA-256 |
|---|---|
| `gabes/schemes/magneto.py` | `3FF0032E5F80003D651DFE27E0F5FA158D1E66549C8EA5661FD0F3E9E87E0F0E` |
| `gabes/zeeman.py` | `543AE7E8888B8D6C491D8E5884A698BF639F300C9D624ED2EB66B8FE345923B2` |
| `gabes/pole_doppler.py` | `06A1255AF017A82A4948D507A4E6A6D2112DC95880F168099FE4571F15724F5A` |
| `tests/test_magneto.py` | `1049940A553C343F87376563DC5802BEF8A0B440468089776BB916233827C32A` |
| `tests/test_magneto_doppler.py` | `EC0EF9C355B1AE3EA91EE092DD82D016A30B190DDF3F0F8D8B986784A2AB19E1` |

중요 재현성 경계:

- exact-Doppler 경로는 현재 working tree에 존재. `gabes/schemes/magneto.py`와
  `gabes/pole_doppler.py`는 modified, `tests/test_magneto_doppler.py`와
  `analysis/pole_reuse/`는 untracked.
- `git diff HEAD --numstat` 기준 Magneto는 `+177/-67`, pole engine은 `+98/-0`.
- clean HEAD의 Magneto는 아직 9-class coarse velocity grid + scalar dilution.
- 현재 README와 실행 결과는 working-tree exact pole–residue 구현 기준.
- 따라서 오늘 수치 재현에는 현재 dirty working tree 필요. **clean checkout 재현성 미확보**.

## 3. 기존 제안·보고서·TODO·issue 먼저 검색

검색 범위:

- 기존 Scheme 4 일일 보고서 **16개**: 2026-06-23부터 2026-10-08.
- `docs/checklist.json` 전체와 retired/superseded 기록.
- `analysis/pole_reuse/DEVLOG.md`, `GATE.md`, `gate_report.json`.
- README, User Guide, code, tests, analysis의 `TODO`, `FIXME`, `issue`, `proposal`,
  `Hanle`, `EIA`, `NMOR`.
- issue/TODO/proposal 이름의 로컬 파일: 별도 파일 없음.
- GitHub repository Issues API, open/closed 합계: **0개**.

결과: **기존 개선안 있음. 충분함. 새 대형 roadmap 항목 불필요.**

현재 직접 관련 항목:

| 항목 | 현재 상태 | 위치 | 오늘 판단 |
|---|---|---|---|
| collision coefficient provenance + pressure shift | P1, ready, large | `docs/checklist.json:2841-2879` | lookup/shift 실행비는 작음. 계수 조사·API·검증 비용 큼 |
| Magneto dynamics semantics | done | `docs/checklist.json:3016-3034` | buffer natural lifetime/elastic dephasing/reload/depolarization 분리 완료 |
| light-region normalization | done | `docs/checklist.json:3037-3056` | conditional trace-one + local density one-use 완료 |
| observable convergence/trust | P1, ready, medium | `docs/checklist.json:3059-3078` | 오늘도 미완료. 가장 직접적인 신뢰성 병목 |
| geometry-aware buffer relaxation | P2, needs_decision, large | `docs/checklist.json:3081-3109` | pressure 하나로 자동 rate 생성 금지. geometry/data 필요 |
| full VCC reference | P2, parked, research | `docs/checklist.json:3246-3265` | interactive default 부적합 |
| pole-residue Doppler reuse | done | `docs/checklist.json:3474-3490` | working tree 구현은 정확·빠름. commit provenance 남음 |
| exact-Doppler 뒤 EIT preset 재검토 | P2, needs_decision, small | `docs/checklist.json:3510-3527` | 현재도 label과 결과 불일치 |

이전 `buffer-gas-pressure-shift`는 sourced coefficient 항목으로 대체됨
`docs/checklist.json:2647-2651`. 이전 `magneto-buffer-relaxation-map`도 pressure-only
자동 매핑을 버리고 geometry-aware budget으로 대체됨 `docs/checklist.json:2653-2657`.
이 결정이 맞음. 계산이 싸도 근거 없는 rate 자동 생성은 물리 개선 아님.

## 4. 현재 Scheme 4가 실제로 푸는 물리

### 4.1 원자와 광장

- 원자: ⁸⁷Rb D1. `Fg, Fe ∈ {1,2}` 중 한 hyperfine pair 선택
  `gabes/schemes/magneto.py:21-30,60-104,384-400`.
- 선택 pair의 모든 `mF` state 포함. σ⁺/σ⁻/π Clebsch–Gordan coupling 구성
  `gabes/zeeman.py:69-123`.
- spontaneous emission을 편광별 grouped jump operator `Σq`로 구성.
  excited-state coherence가 ground-state coherence로 전달. EIA/LCA에 필요한
  transfer of coherence(TOC) 포함 `gabes/zeeman.py:101-123`.
- QWP 각도를 complex circular-basis drive weight로 변환
  `gabes/schemes/magneto.py:107-118`.
- longitudinal B scan, coil/shield offset, residual transverse B, transverse azimuth를
  Zeeman Hamiltonian에 직접 반영 `gabes/schemes/magneto.py:430-450,556-574`.
- intensity를 Steck D1 saturation scale로 Rabi frequency에 변환
  `gabes/schemes/magneto.py:27-30,446-450`.

### 4.2 완화와 cell

- buffer cell: natural emission/TOC는 `GAMMA_D1`에 고정. Ne broadening 절반을
  elastic optical-coherence dephasing으로 넣음. isotropic transit reservoir reload와
  collisional ground dephasing 분리 `gabes/schemes/magneto.py:401-417`,
  `gabes/zeeman.py:125-149`.
- paraffin cell: illuminated region + dark region coupled steady state. dark region은
  unlit. wall coherence lifetime, bright exit, dark return 포함
  `gabes/schemes/magneto.py:418-428,586-635`.
- two-region light block을 trace-one conditional optical state로 정규화.
  light-region occupation은 diagnostic. local vapor density는 susceptibility에 한 번만 적용
  `gabes/schemes/magneto.py:577-635,548-552`.

### 4.3 Doppler와 관측량

- Doppler-on: B-dependent affine real-Hermitian Liouvillian → shift-invert poles/residues
  → Faddeeva Gaussian mean. velocity cutoff/grid 없는 Maxwell average
  `gabes/schemes/magneto.py:663-747`, `gabes/pole_doppler.py:253-337`.
- 각 B row를 direct solve로 guard. 실패 row만 ±6σ converged grid fallback
  `gabes/schemes/magneto.py:676-747`.
- transmission: `α = k Imχ`, `T = exp(-αL)`.
- rotation: `θ = kL Re(χ+−χ−)/4`
  `gabes/schemes/magneto.py:788-823`.

결론: toy Lorentzian 아님. Zeeman selection rule, optical pumping, spontaneous
coherence transfer, transverse-field switching, two-region Ramsey core, thermal average,
Beer–Lambert/NMOR 연결이 실제로 있음.

## 5. 실험 원자물리 연구자 관점 평가

### 5.1 유용한 용도

- zero-field feature의 dip/peak/crossover 부호 sanity check.
- linear→elliptical→circular polarization 추세.
- residual transverse field와 longitudinal zero offset 영향 확인.
- wall coherence, dark return, buffer reload/dephasing의 상대 민감도.
- B scan 범위, detuning, intensity, cell length의 초기 설계.
- measured trace fitting 전 parameter seed와 모델 누락 진단.
- exact Maxwell 평균 덕분에 옛 9-class aliasing은 제거됨. `analysis/pole_reuse/DEVLOG.md:136-158`.

### 5.2 원 논문과의 일치 범위

내장 reference 세 개의 적용 범위가 서로 다름 `gabes/schemes/magneto.py:373-381`.

- Lee–Moon 2013: ⁸⁷Rb D1 `F=2→F'=1`, paraffin Hanle. 현 compact model과 직접
  비교 가능한 가장 가까운 anchor. 논문은 linear CPT 폭 약 0.12 mG, circular MIA
  폭 약 0.20 mG 보고
  ([Optica, DOI 10.1364/JOSAB.30.002301](https://doi.org/10.1364/JOSAB.30.002301)).
- Yu et al. 2010: ⁸⁷Rb D1 + Ne buffer gas, circular LCA, residual transverse field와
  diffusive interaction time 중요. 2.4 mG 보고
  ([APS, DOI 10.1103/PhysRevA.81.023416](https://doi.org/10.1103/PhysRevA.81.023416)).
  현재 buffer LCA topology 검증에는 적절. diffusion을 생략한 절대 폭 anchor로는 부족.
- Moon–Kim 2014: paraffin Ramsey EIA→MIT이지만 **⁸⁷Rb D2 `F=2→F'=3`**.
  현재 Scheme 4는 D1, `Fe∈{1,2}`만 지원. 따라서 coated-cell switching의 qualitative
  reference일 뿐 D1 preset의 숫자 검증 자료가 아님
  ([DOI 10.1364/OE.22.018604](https://doi.org/10.1364/OE.22.018604)).

이 source/line/transition 적용 범위를 UI와 검증 출력에 같이 표시해야 함. 계산비용 O(1).
기존 `magneto-observable-convergence-and-trust` 범위 안에서 해결 가능.

### 5.3 절대 experimental reference가 아닌 이유

1. **선택한 hyperfine pair만 닫힘.** 실제 D1 open transition의 다른 ground F 누출,
   repump, 두 excited-hyperfine interference 없음 `gabes/zeeman.py:69-130`.
2. **충돌 provenance 부족.** Ne는 고정 broadening scalar. pressure shift, gas/line/T
   의존성, diffusion, spin-destruction uncertainty, Dicke narrowing, VCC 없음
   `gabes/constants.py:79-89`.
3. **공간 평균 없음.** Gaussian transverse intensity, diffusion, field gradient, finite
   cell geometry, beam/cell volume ratio 없음.
4. **paraffin transit 의미 중복.** `transit_relax_khz`가 light→dark exchange와 light-region
   ground pure dephasing에 동시에 사용 `gabes/schemes/magneto.py:419-424,621-625`.
5. **rotation 식 적용 범위 좁음.** weak/thin circular birefringence 식. 일반 Jones/Stokes
   propagation 아님. pure circular input에는 polarization plane 자체가 없음.
6. **절대 scale 미검증.** `line_strength`가 0.01–2 calibration knob. 코드도 absolute
   contrast/rotation/slope를 `external validation required`로 표시
   `gabes/schemes/magneto.py:41-43,280-284,825-857`.
7. **magnetometer sensitivity 미제공.** measured power, detector transfer/noise,
   calibrated coil conversion, held-out trace 없음.
8. **clean-checkout provenance 없음.** exact-Doppler production/test/DEVLOG이 HEAD 밖.

판정: **실제 물리를 가진 정성·반정량 experimental-planning/reference**.
**absolute magnetometry, publication-grade linewidth/rotation/sensitivity reference 아님.**

## 6. 오늘 수치 재검증

동일 process warm 5회 median. 현재 working-tree exact Maxwell 경로.
모든 preset에서 pole guard ≤ `2.97e-12`, fallback row 0.

| regime | 현재 feature | `α(0)-wing` [m⁻¹] | central width [µT] | T(command 0) | median |
|---|---|---:|---:|---:|---:|
| EIT dip | **crossover** | +0.357285 | 0.981271 | 0.787789 | 73.24 ms |
| EIA peak | peak | +14.930043 | 0.363514 | 0.856642 | 73.34 ms |
| Buffer Hanle | dip | −4.670913 | 4.308938 | 0.757571 | 39.01 ms |
| Buffer LCA | crossover | +0.012284 | 1.693135 | 0.629190 | 89.28 ms |
| NMOR | absorption dip | −2.816494 | 0.113283 | 0.827167 | 76.82 ms |

### 6.1 EIT preset truth

`EIT dip` button은 exact Maxwell 평균에서 계속 crossover. checklist 결정 항목 정확
`docs/checklist.json:3510-3527`. 현재 label 그대로 두려면 model result와 모순.
무근거 retune보다 우선 **crossover로 relabel** 또는 explicit `not a dip under current model`
status가 안전.

`Buffer LCA`도 현 built-in point에서는 `α(0)-wing=+0.012284 m⁻¹`뿐이고 코드 분류는
**crossover**. LCA topology 자체는 transverse-field test로 검증되지만, 기본 버튼이
명료한 absorption peak를 보인다는 뜻은 아님. 같은 preset-truth 점검에 포함할 가치 있음.

### 6.2 NMOR B-grid 수렴

현재 hero slope는 중앙 `np.gradient`, 사실상 3-point
`gabes/schemes/magneto.py:800,817-823`.

| B points | width [µT] | 현행 3-point [mrad/µT] | 5-point [mrad/µT] | runtime |
|---:|---:|---:|---:|---:|
| 121 | 0.113283 | −252.147 | −285.285 | 78.79 ms |
| 201 | 0.110669 | −293.582 | −314.868 | 112.57 ms |
| 401 | 0.109904 | −314.602 | −321.608 | 214.88 ms |
| 801 | 0.109482 | −320.119 | −321.958 | 472.17 ms |

- 기본 hero 절댓값은 801-point 5-point 기준보다 **21.7% 작음**.
- 401-point 5-point는 801-point 5-point와 **0.11%** 차이.
- 121-point 5-point도 약 11.4% 작음.
- 전체 기본 grid를 401점으로 늘리면 runtime 약 2.7배. negligible overhead 아님.
- physical zero 주위에 별도 4 pole row를 추가하면 row 수 +3.3%. 기존 checklist의
  locally converged derivative가 가장 좋은 비용/정확도 절충.

### 6.3 paraffin transit 의미 중복

201-point EIT preset. light-region `gamma_gg`만 0으로 둔 비저장 진단.
exchange rate와 matrix dimension/solve count는 유지.

| | 현행 | exchange-only 진단 |
|---|---:|---:|
| feature | crossover | dip |
| `α(0)-wing` [m⁻¹] | +0.357285 | −4.428289 |
| central width [µT] | 0.980818 | 0.123955 |
| T(B=0) | 0.787789 | 0.857372 |

`χ_probe` 최대 상대 변화 **35.04%**. runtime overhead는 사실상 0이지만 결과 변화 큼.
따라서 단순 코드 정리 아님. bright exchange와 independent light-region dephasing의
기본값/provenance 결정 필요. 물리 채널 분리는 강하게 권고.

### 6.4 command zero와 physical zero

`b_offset=+0.25 µT`, 401-point EIT:

- command zero = physical `+0.25 µT`: `α=23.906932 m⁻¹`, `T=0.787360`.
- physical zero = command `−0.25 µT`: `α=23.852542 m⁻¹`, `T=0.787789`.

Hamiltonian은 physical axis 사용 `gabes/schemes/magneto.py:430-434`.
hero는 command-axis 중앙 사용 `gabes/schemes/magneto.py:788-823`.
두 값 표시/보간은 O(nB), atomic re-solve 0.

### 6.5 circular NMOR gate

QWP 45°, 401 points:

- 현행 3-point slope: `−41.369 mrad/µT`.
- 5-point slope: `−42.382 mrad/µT`.
- peak `|rotation| = 5.944 mrad`.

하지만 pure circular input에는 polarization-plane angle 정의 없음. 숫자를 hero로 내보내면
비물리적 의미 부여. linear/near-linear gate 또는 Stokes/Jones readout 필요. gate 비용 O(1).

### 6.6 문헌 검증 스크립트

`python tests/verify_hanle_eit_eia.py` fresh output:

- paraffin linear: `0.950 mG` vs Lee–Moon `0.12 mG` → **7.9배**.
- paraffin circular: `2.305 mG` vs `0.20 mG` → **11.5배**.
- buffer low-power circular LCA: `3.349 mG` vs Yu `2.4 mG` → **39.5% 큼**.

스크립트는 sign switch를 `MATCH`라 하고, 0.95–2.30 mG를 논문 0.12–0.20 mG와
`same sub-mG order`라고 출력 `tests/verify_hanle_eit_eia.py:243-271`.
2.305 mG는 sub-mG 아님. **trend/sign PASS와 absolute-width CHECK/FAIL 분리 필요.**
계산비용 0.

## 7. 기존 물리 개선안의 비용과 우선순위

| 기존 제안 | 계산비용 | 물리 보존 / 판단 |
|---|---:|---|
| exact pole-residue Maxwell | 현재 39–89 ms/preset | 옛 9-class보다 정확·빠름. 먼저 review/commit해 provenance 고정 |
| physical/command zero 분리 | O(nB), 재풀이 0 | 의미 오류 제거. 즉시 가능 |
| literature/source applicability label | O(1) | D1/D2·transition mismatch까지 표시. 즉시 가능 |
| open-hyperfine/repump scope status | O(1) | 계산하지 않은 물리 명시. 즉시 가능 |
| circular NMOR gate | O(1) | 정의 불가 claim 차단. 즉시 가능 |
| local NMOR slope | 추가 B row 약 4개, 기본 대비 약 +3.3% | atomic model 불변. hero 수치 신뢰 개선 |
| EIT preset relabel | 0 | 현재 계산을 정직하게 표시. retune은 rate 의미 분리/데이터 뒤 |
| transit exchange/dephasing 분리 | 차원·solve 수 동일 | computational overhead 거의 0. 결과 35% 변화. independent rate 필요 |
| sourced collision table + pressure shift | lookup + scalar detuning shift, O(1) | compact physics 개선. source/API/validation effort large |
| geometry-aware buffer budget | rate 계산 자체 작음 | geometry·diffusion·wall·uncertainty 입력 없으면 false precision |
| full selected D1 manifold | 현재 F2→F'1: `n=8, M=64`; full D1 약 `n=16, M=256` | dense LU 대략 64× flops, 16× matrix memory. opt-in reference만 |
| calibrated coated transport | dark-time/spatial quadrature 필요 | research 비용. 실측 geometry/trace 먼저 |
| full VCC | velocity classes 상호 결합 | exact separable pole 경로 깨짐. interactive default 부적합 |

무시 가능한 overhead로 물리를 보존/개선하는 순서:

1. zero/status/source-scope/circular gate.
2. local-converged NMOR derivative.
3. transit 채널 의미 분리. 단, 기본값은 데이터/명시적 가정 필요.
4. sourced pressure shift와 elastic-dephasing coefficient.

## 8. 동작 불변 순수 코드 최적화

### 8.1 exact path B-stack 지연 생성 — 최우선

현재 `L0_all`, `Ld_all` complex stack을 Doppler mode와 무관하게 먼저 생성
`gabes/schemes/magneto.py:469-485`.

정상 exact path는 `C_xy`, `C_z`, `Cd_xy`만 사용
`gabes/schemes/magneto.py:509-512,663-747`.
stack은 Doppler-off 또는 failed-guard grid fallback의 `solve()`에서만 필요
`gabes/schemes/magneto.py:487-507,517-528`.

오늘 121-point paraffin F2→F'1 benchmark:

- 두 stack: **15.125 MiB**.
- eager assembly median **21.03 ms**, p95 **35.77 ms**.
- 요청 row broadcast와 full eager result: `np.array_equal=True`.
- fallback 3 row만 lazy assembly: median **0.129 ms**.

`solve(rows, deff)` 안에서 필요한 row만 같은 식으로 조립하면 정상 exact path에서
allocation/assembly 전부 제거. 수학·합 순서·출력 불변.

### 8.2 BLAS threadpool controller 재사용 — 중간 우선순위

`core.blas_single_thread()`가 매 호출 새 threadpool limit context 생성
`gabes/core.py:17-34`.

cached `ThreadpoolController` 비저장 prototype:

- 빈 context median: `2.692 → 0.0113 ms`.
- EIT compute median: `68.453 → 64.712 ms`, 약 **5.5% 감소**.
- `chi_probe`, `chi_p`, `chi_m`: bitwise 동일.

적용 전 late-loaded BLAS, nested context, concurrent call, process-fork/thread safety 회귀 필요.

### 8.3 angular-momentum matrix cache — 낮은 우선순위

`angular_momentum_matrices(F)`가 매번 배열 생성 `gabes/zeeman.py:48-66`.

- F=2/F=1 pair median: `51.8 → 0.4 µs`.
- 전체 compute 대비 매우 작음.
- read-only `lru_cache` + mutation regression이면 안전. 우선순위 낮음.

이미 된 최적화 재제안 안 함:

- `zeeman_manifold()` cache `gabes/zeeman.py:69-70`.
- Hermitian basis cache `gabes/core.py:61-94`.
- exact real-basis pole readout `gabes/schemes/magneto.py:681-747`.
- figureless observables path.

## 9. 코드·테스트·문서·예제 평가

### 9.1 테스트 강점

- CG known values, emission normalization, TOC source
  `tests/test_magneto.py:28-63`.
- elastic optical dephasing이 lifetime/branching/TOC를 보존
  `tests/test_magneto.py:66-103`.
- transit reload와 collisional dephasing 채널 분리, positivity/trace/Hermiticity
  `tests/test_magneto.py:105-159,352-446`.
- conditional light-state normalization, occupation, local-density one-use
  `tests/test_magneto.py:233-302`.
- paraffin dip/peak, wall narrowing, buffer broadening, intrinsic EIA, LCA transverse-field
  requirement, NMOR odd crossing `tests/test_magneto.py:305-503`.
- exact Maxwell vs converged fallback grid, off-path isolation, method provenance
  `tests/test_magneto_doppler.py:54-114`.
- Numba grid kernel vs NumPy parity `tests/test_kernels.py:62-142`.

남은 test gap:

- 121→401/801 B-grid/local slope convergence.
- physical-zero hero와 command-zero secondary readout.
- circular/near-circular NMOR gate.
- open-hyperfine/repump scope status.
- literature sign/trend와 absolute-width tolerance 분리.
- exact-Doppler test/DEVLOG의 committed clean-checkout 실행.
- generated example source/config/model/git hash.

낡은 인자:

- production schema에서 `velocity_classes` 제거됨 `gabes/schemes/magneto.py:285-293`.
- 하지만 `tests/test_magneto.py:166-169,412-421`,
  `tests/verify_hanle_eit_eia.py:71-81,98-107,157-187`,
  `analysis/squeezing/resonant_hanle_experiment_config.example.json:15-18`,
  `analysis/squeezing/resonant_hanle_squeezing_reference.py:370-397,1648-1651`가 계속 전달.
- extra dict key라 조용히 무시. 수렴 knob 검증처럼 오해 가능. 제거 필요.

### 9.2 README와 User Guide

- README는 exact Maxwell 경로를 설명 `README.md:519-524`. 현재 working tree와 일치.
  그러나 clean HEAD와 불일치.
- User Guide는 Magneto를 `~0.35 s`, `B scan × velocity`로 설명
  `docs/Userguide/GABES_User_Guide_v2.html:778-790`. 현 정상 경로는 pole 평균, velocity grid 없음.
- Magneto 한계는 “실험별 조정” 한 줄뿐
  `docs/Userguide/GABES_User_Guide_v2.html:852-854`.
  open hyperfine, absolute scale, circular NMOR, collision provenance, B-grid 미수렴 누락.
- “절대 스케일, 선폭, NMOR 영점 교차 자동 확인” 문구
  `docs/Userguide/GABES_User_Guide_v2.html:935-937`는 실제 검증 범위보다 강함.
- static `docs/Userguide/userguide_assets/hanle.png`는 2026-06-08 생성.
  현재와 같은 표면 조건(F2→F'1, QWP 0°, 0.8 mW/cm², ±2 µT)이지만
  예전 sharp transmission peak `T≈0.9985` 표시. 현재 exact 경로는 `T(0)=0.787789`,
  crossover. 이미지 stale.

### 9.3 Hanle–FWM analysis example

`analysis/squeezing/resonant_hanle_squeezing_reference.py`와 example JSON이 유일한
실질적 Scheme 4 batch example. standalone Magneto-only example 없음.

장점:

- beam power/waist→intensity, Hanle transmission slope, balanced detector shot/electronic
  noise, optional measured CSV affine calibration을 한 workflow에 연결
  `...reference.py:357-420,665-735`.
- frequency anchor, responsivity, CSV unit parsing/fit에 unit tests 있음
  `tests/test_resonant_hanle_reference.py:14-141`.

한계:

- `velocity_classes` obsolete.
- transmission slope가 coarse `np.gradient` 3-point
  `analysis/squeezing/resonant_hanle_squeezing_reference.py:401-413`.
- 기본 `source_anchor=paper`; OE gain/IDS를 paper anchor로 주입
  `...reference.py:637-656,1629-1633`.
- 저장 Markdown은 measured Hanle CSV 없음 명시
  `analysis/squeezing/resonant_hanle_squeezing_reference.md:45-60`.
  그런데 ideal/zero electronics defaults 아래 구체 pT/√Hz 계속 출력.
- example report에 생성시각, git commit/dirty, Magneto cache version, source/config hash 없음
  `...reference.py:1019-1173`.

저장된 sensitivity 숫자는 현재 Scheme 4의 experimental validation 아님.
조건부 system study로만 사용.

## 10. 권고 순서

1. **P0 재현성:** exact-Doppler production code, test, pole-reuse evidence review/commit.
2. **P1 trust:** physical/command zero 분리, local-converged NMOR slope,
   B-grid status, 문헌 trend/absolute 판정 분리.
3. **P1 source scope:** 각 reference에 isotope/line/transition/cell/readout applicability 표시.
4. **P1 semantics:** paraffin exchange와 light pure dephasing 분리.
5. **P1 preset truth:** 4번 뒤 EIT/Buffer-LCA preset 근거 retune 또는 crossover relabel.
6. **P1 observable scope:** circular/near-circular NMOR plane-rotation gate.
7. **P1 docs/example:** User Guide/image와 Hanle–FWM report 재생성. provenance hash 추가.
8. **P1 collision:** sourced gas/species/line broadening + pressure shift.
9. **P1 pure speed:** lazy B-stack. 그 뒤 BLAS controller reuse.
10. **보류:** full VCC, calibrated spatial transport, full open D1/repump,
    full Jones/Stokes propagation. dataset + runtime budget 있는 opt-in reference만.

## 11. 검증 기록

- `CLAUDE.md`, `README.md` 전체 읽음.
- registry/UI 순서 확인.
- 기존 Scheme 4 보고서 16개, checklist, pole-reuse 기록, TODO/FIXME/local issue-note 검색.
- GitHub open/closed issue 0 확인.
- current Magneto/Zeeman/pole code, tests, User Guide, static example image,
  Hanle–FWM source/config/report 읽음.
- primary literature metadata/abstract 재확인.
- 5 regime warm timing, pole guard/fallback 확인.
- 121/201/401/801 NMOR slope/width 재측정.
- transit dephasing 제거 비저장 진단, offset, QWP 45° NMOR 재측정.
- lazy B-stack, cached BLAS controller, angular-matrix cache 비저장 benchmark.
- `python tests/verify_hanle_eit_eia.py` 실행.
- focused pytest:
  `python -m pytest -q tests/test_magneto.py tests/test_magneto_doppler.py tests/test_kernels.py tests/test_headless_observables.py tests/test_schemes_render.py tests/test_resonant_hanle_reference.py tests/test_docs_consistency.py`
  → **85 passed in 45.46 s**. 전체 pytest와 동시 실행해 직전 단독 실행보다 느림.
- required root pytest: `python -m pytest -q`
  → **2184 passed, 4 skipped, 20 failed in 4332.63 s**. 실패는 Scheme 4 집중
  테스트가 아니라 Grand Challenge thermal-campaign 격리 capsule/grid harness 19건
  (그중 1건은 기존 Windows nested-destination `FileNotFoundError`)과 기존 FWM
  Streamlit hidden-noise-slider 상태 유지 1건이다. 즉 Magneto 회귀는 위 85건에서
  모두 통과했지만, 현재 dirty worktree의 저장소 전체 green 상태를 주장할 수는 없다.
- 생산 코드 수정 없음. 이 보고서와 automation memory만 갱신 예정.

## 12. 최종 판정

Scheme 4는 유용한 실제 원자물리 구현. Zeeman CG, optical pumping, grouped-emission
TOC, transverse-field switching, two-region coated-cell coherence, buffer relaxation,
exact Maxwell average, Beer–Lambert/NMOR readout이 있음.

실험 설계, feature 부호, 상대 추세, 초기 fit seed에 적합.

절대 reference로는 부적합. 기본 EIT/Buffer-LCA label 불일치, NMOR slope 21.7% 미수렴,
transit 의미 중복, physical-zero 혼동, circular plane-rotation 출력, paraffin 문헌 폭
7.9–11.5배 불일치, open-hyperfine/repump·collision·spatial 평균·absolute calibration 부족.
Moon–Kim 2014 D2 reference도 현 D1 scheme의 숫자 anchor가 아님. exact-Doppler는 HEAD 밖.

등급: **정성·반정량 experimental atomic-physics reference / planning tool**.
**absolute magnetometry, publication-grade linewidth/rotation/sensitivity reference로 부적합**.

최고 비용효율 작업은 새 대형 모델 아님.
먼저 provenance, zero/slope/status, source applicability, transit 의미, observable gate,
lazy B-stack.
