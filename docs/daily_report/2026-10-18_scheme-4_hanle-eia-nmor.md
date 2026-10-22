# 2026-10-18 Scheme 4 물리 검토 — Hanle / EIA / NMOR

## 1. 선택과 현재 다섯 scheme 순서

- 로컬 날짜: **2026-10-18**, Asia/Seoul.
- 계산: **n = (18 mod 5) + 1 = 4**.
- 대상: **Scheme 4, Magneto-optics (Hanle/MOR)**.
- 기준: 현재 working tree, branch `main`, HEAD `39be4ae`.
- 저장소 dirty. 기존 변경 보존. 생산 코드 수정 없음.

현재 UI 순서:

| n | key | 정의 |
|---:|---|---|
| 1 | `sas` | Absorption spectroscopy (OD / SAS) |
| 2 | `lambda` | Λ coherence (EIT / AT / CPT) |
| 3 | `rydberg_eit` | Rydberg-EIT electrometry |
| 4 | `magneto` | Magneto-optics (Hanle/MOR) |
| 5 | `fwm` | Four-wave mixing (Squeezing / Biphoton) |

근거:

- registry 순서: `gabes/schemes/__init__.py:19-25`.
- `all_schemes()`가 그 목록 복사 반환: `gabes/schemes/__init__.py:37-39`.
- UI가 같은 title 배열과 selectbox 사용: `streamlit_app.py:132-150`.
- 사용자-facing 정의: `README.md:8-16`.
- Scheme 4 key/title/caption: `gabes/schemes/magneto.py:151-167`.

## 2. 기존 기록 먼저 검색한 결과

검색 범위:

- Scheme 4 일일 보고서 **17개**: 2026-06-23부터 2026-10-13.
- `docs/checklist.json`의 active/done/superseded/parked 항목.
- `analysis/pole_reuse/DEVLOG.md`와 gate 기록.
- README, User Guide, tests, analysis 예제의 `TODO`, `FIXME`, `issue`,
  `proposal`, `Hanle`, `EIA`, `NMOR`.
- 로컬 TODO/issue 이름 파일. Magneto 전용 별도 issue 파일 없음.

기존 개선안 있음. 새 대형 roadmap 불필요.

현재 직접 관련 registry:

| 항목 | 상태 | 근거 |
|---|---|---|
| provenance-backed collision coefficient + pressure shift | P1, ready, large | `docs/checklist.json:3163-3197` |
| Magneto dynamics semantics | done | `docs/checklist.json:3338-3357` |
| light-region normalization | done | `docs/checklist.json:3359-3379` |
| observable convergence/trust | P1, ready, medium | `docs/checklist.json:3381-3400` |
| geometry-aware buffer relaxation | P2, needs decision, large | `docs/checklist.json:3403-3432` |
| full VCC reference | P2, parked, research | `docs/checklist.json:3568-3588` |
| exact-Doppler 뒤 EIT preset 검토 | P2, needs decision, small | `docs/checklist.json:3832-3849` |

옛 `magneto-buffer-relaxation-map`은 pressure 하나로 relaxation을 자동 생성하지 말라는
이유로 geometry-aware budget에 대체됨 `docs/checklist.json:2885-2888,3403-3432`.
계산이 싸도 입력 근거 없으면 물리 개선 아님.

## 3. 변경·재현 범위

직전 Scheme 4 보고서: `docs/daily_report/2026-10-13_scheme-4_hanle-eia-nmor.md`.
그 뒤 핵심 Magneto 코드와 전용 테스트 byte 동일.

| 파일 | 현재 SHA-256 |
|---|---|
| `gabes/schemes/magneto.py` | `3FF0032E5F80003D651DFE27E0F5FA158D1E66549C8EA5661FD0F3E9E87E0F0E` |
| `gabes/zeeman.py` | `543AE7E8888B8D6C491D8E5884A698BF639F300C9D624ED2EB66B8FE345923B2` |
| `gabes/pole_doppler.py` | `06A1255AF017A82A4948D507A4E6A6D2112DC95880F168099FE4571F15724F5A` |
| `tests/test_magneto.py` | `1049940A553C343F87376563DC5802BEF8A0B440468089776BB916233827C32A` |
| `tests/test_magneto_doppler.py` | `EC0EF9C355B1AE3EA91EE092DD82D016A30B190DDF3F0F8D8B986784A2AB19E1` |

재현성 경계:

- exact-Doppler production은 HEAD 대비 `magneto.py +177/-67`,
  `pole_doppler.py +98/-0`.
- `tests/test_magneto_doppler.py`, `analysis/pole_reuse/`는 untracked.
- clean HEAD는 옛 9-class velocity grid + scalar dilution.
- 오늘 수치는 dirty working tree 필요. clean-checkout 재현성 아직 없음.

## 4. 현재 Scheme 4가 푸는 물리

### 4.1 원자·광장

- ⁸⁷Rb D1. `Fg, Fe ∈ {1,2}` 중 한 hyperfine pair 선택
  `gabes/schemes/magneto.py:21-30,60-104,384-400`.
- 선택 pair의 모든 `mF` 상태. σ⁺/σ⁻/π CG coupling
  `gabes/zeeman.py:69-123`.
- 편광별 grouped spontaneous-emission jump. excited coherence의 ground
  coherence transfer 포함. EIA/LCA 핵심 `gabes/zeeman.py:101-156`.
- QWP를 circular-basis complex drive weight로 변환
  `gabes/schemes/magneto.py:107-118`.
- longitudinal scan, B-zero offset, residual transverse field, azimuth을 Zeeman
  Hamiltonian에 반영 `gabes/schemes/magneto.py:214-219,268-279,432-450`.
- intensity를 ⁸⁷Rb D1 saturation scale로 Rabi frequency 변환
  `gabes/schemes/magneto.py:27-30,446-450`.

### 4.2 완화·셀

- buffer cell: 자연 방출/TOC는 `GAMMA_D1`. Ne broadening은 elastic optical
  dephasing. transit reservoir reload와 collisional ground dephasing 분리
  `gabes/schemes/magneto.py:402-417`.
- paraffin cell: illuminated region + dark region coupled steady state.
  bright exit, dark return, wall-coherence lifetime 포함
  `gabes/schemes/magneto.py:418-428,586-635`.
- light block은 trace-one conditional state. occupation 별도 diagnostic.
  local vapor density는 susceptibility에 한 번 적용
  `gabes/schemes/magneto.py:532-552,586-635`.

### 4.3 Doppler·관측량

- Doppler-on: affine real-Hermitian Liouvillian → shift-invert poles/residues →
  Faddeeva Gaussian mean. velocity cutoff/grid 없는 Maxwell average
  `gabes/schemes/magneto.py:663-746`, `gabes/pole_doppler.py:119-153,253-337`.
- 각 B row direct guard. 실패 row만 ±6σ grid fallback
  `gabes/schemes/magneto.py:711-746`.
- transmission: `α = k Imχ`, `T = exp(-αL)`.
- rotation: `θ = kL Re(χ+−χ−)/4`
  `gabes/schemes/magneto.py:789-823`.

Toy Lorentzian 아님. Zeeman selection rule, optical pumping, TOC, transverse-field
switching, two-region Ramsey core, exact thermal average, Beer–Lambert/NMOR 연결 있음.

## 5. 오늘 수치 확인

현재 working tree. 동일 process 3회 median. 시간은 동시 전체 테스트 부하 영향 있음.
모든 preset에서 pole fallback 0.

| regime | 현재 feature | `α(0)-wing` [m⁻¹] | `T(command 0)` | pole guard max | median |
|---|---|---:|---:|---:|---:|
| EIT dip | **crossover** | +0.357285 | 0.787789 | 2.56e-12 | 163.7 ms |
| EIA peak | peak | +14.930043 | 0.856642 | 6.19e-14 | 119.8 ms |
| Buffer Hanle | dip | −4.670913 | 0.757571 | 1.47e-13 | 50.7 ms |
| Buffer LCA | **crossover** | +0.012284 | 0.629190 | 2.25e-12 | 132.1 ms |
| NMOR | absorption dip | −2.816494 | 0.827167 | 2.97e-12 | 136.5 ms |

NMOR hero: slope `−252.15 mrad/µT`, peak rotation `10.94 mrad`.
직전 보고의 해시 동일 수렴 감사 재사용:

- 기본 121-point 3-point slope는 801-point 5-point `−321.958 mrad/µT`보다
  절댓값 **21.7% 작음** `docs/daily_report/2026-10-13_scheme-4_hanle-eia-nmor.md:205-228`.
- paraffin light-region `gamma_gg` 제거 진단은 같은 solve 차원에서 EIT preset을
  crossover `+0.357285 m⁻¹`에서 dip `−4.428289 m⁻¹`로 바꾸고 χ를 35.04% 변경
  `docs/daily_report/2026-10-13_scheme-4_hanle-eia-nmor.md:229-243`.
- QWP 45° NMOR도 plane-rotation hero를 출력. pure circular input에는 편광면 각도 정의 없음
  `docs/daily_report/2026-10-13_scheme-4_hanle-eia-nmor.md:255-265`.

`python tests/verify_hanle_eit_eia.py` fresh 결과:

- paraffin linear: `0.950 mG` vs 문헌 `0.12 mG` → 7.9배.
- paraffin circular: `2.305 mG` vs 문헌 `0.20 mG` → 11.5배.
- buffer low-power circular LCA: `3.349 mG` vs 문헌 `2.4 mG` → 39.5% 큼.
- 스크립트는 sign switch를 `MATCH`, 폭 0.95–2.30 mG를 `same sub-mG order`라 출력
  `tests/verify_hanle_eit_eia.py:248-271`. 2.305 mG는 sub-mG 아님.

## 6. ⟦AI⟧ 실험 원자물리 평가

### 6.1 유용한 범위

- zero-field dip/peak/crossover 부호 sanity check.
- linear→elliptical→circular polarization 추세.
- residual transverse field와 longitudinal offset 영향.
- wall coherence, dark return, buffer reload/dephasing 상대 민감도.
- B scan, detuning, intensity, cell length의 초기 설계.
- measured trace fitting 전 parameter seed와 누락 물리 진단.
- exact Maxwell 평균으로 옛 9-class aliasing 제거
  `analysis/pole_reuse/DEVLOG.md:20-32,136-158`.

### 6.2 absolute reference가 아닌 이유

1. 선택 hyperfine pair만 닫힘. 실제 D1 open transition의 다른 ground F 누출,
   repump, excited-hyperfine interference 없음 `gabes/zeeman.py:69-156`.
2. Ne는 fixed broadening scalar. pressure shift, gas/line/T provenance, diffusion,
   spin-destruction uncertainty, Dicke narrowing, VCC 없음
   `gabes/constants.py:86-89`, `docs/checklist.json:3163-3197`.
3. Gaussian transverse intensity, diffusion, field gradient, finite cell geometry,
   beam/cell volume 평균 없음.
4. paraffin `transit_relax_khz`가 bright→dark exchange와 light-region ground
   dephasing에 함께 들어감 `gabes/schemes/magneto.py:418-426`.
5. rotation은 weak/thin circular-birefringence 식. 일반 Jones/Stokes propagation 아님.
6. `line_strength`는 0.01–2 calibration knob. 코드가 absolute contrast/rotation/slope를
   `external validation required`로 표시 `gabes/schemes/magneto.py:43,280-284,825-832`.
7. measured power, detector transfer/noise, coil calibration, held-out trace 없음.
8. exact-Doppler 구현·전용 테스트가 clean HEAD 밖.

판정: **실제 물리를 가진 정성·반정량 experimental-planning/reference**.
**absolute magnetometry, publication-grade linewidth/rotation/sensitivity reference 아님.**

## 7. ⟦AI⟧ 기존 물리 개선안 비용 평가

| 기존 제안 | 계산비용 | 물리 보존·효과 |
|---|---:|---|
| exact pole-residue Maxwell | 현재 0 fallback, 대화형 속도 | 옛 grid보다 정확·빠름. 우선 clean-checkout provenance 고정 |
| command/physical zero 분리 | O(nB), atomic solve 0 | 의미 오류 제거. 즉시 가능 |
| source/applicability·open-hyperfine status | O(1) | 계산하지 않은 범위 명시. 즉시 가능 |
| circular NMOR gate | O(1) | 정의 불가 plane-angle claim 차단 |
| local-converged NMOR slope | 약 4개 local B row, 기본 대비 약 +3.3% | 원자 모델 불변. hero 신뢰 개선 |
| EIT/Buffer-LCA preset relabel | 0 | 현재 converged 결과와 label 일치. 무근거 retune보다 안전 |
| bright exchange / light dephasing 분리 | matrix 차원·solve 수 동일 | 연산 overhead 거의 0. 결과는 35% 변함. 독립 rate/provenance 필요 |
| sourced pressure shift·elastic coefficient | lookup + scalar detuning shift, O(1) | compact physics 개선. 문헌/API/validation 노력은 큼 |
| geometry-aware buffer budget | rate 계산 자체 작음 | geometry·diffusion·wall·uncertainty 없으면 false precision |
| full selected D1 manifold | n≈8→16, dense LU 대략 64× flops | default에 비쌈. opt-in reference만 |
| calibrated coated transport | dark-time/spatial quadrature | 연구 비용 큼. geometry/trace 먼저 |
| full VCC | velocity classes 상호 결합 | exact separable pole 경로 깨짐. negligible 아님 |

무시 가능한 runtime으로 먼저 할 것:

1. zero/status/source-scope/circular gate.
2. local-converged NMOR derivative.
3. transit channel 의미 분리. 단, 기본값은 independent data 또는 명시적 가정 필요.
4. sourced pressure shift와 elastic-dephasing coefficient.

## 8. ⟦AI⟧ 동작 불변 순수 코드 최적화

### 8.1 exact path B-stack 지연 생성 — 최우선 미구현

현재 `L0_all`, `Ld_all` full complex stack을 먼저 생성
`gabes/schemes/magneto.py:479-485`.
정상 Doppler exact path는 `C_xy`, `C_z`, `Cd_xy`만 사용. stack은 Doppler-off 또는
failed-guard grid fallback의 `solve()`에서만 필요
`gabes/schemes/magneto.py:487-512,663-746`.

해시 동일 직전 benchmark:

- paraffin 121-point 두 stack 15.125 MiB.
- eager assembly median 21.03 ms.
- fallback 3 row lazy assembly 0.129 ms.
- requested-row broadcast와 eager result array-identical.

`solve(rows, deff)` 안에서 필요한 row만 같은 식으로 조립. 수학·합 순서·출력 불변.

### 8.2 BLAS threadpool controller 재사용 — 중간 우선순위

`core.blas_single_thread()`가 호출마다 context 생성 `gabes/core.py:23-34`.
직전 해시 동일 prototype은 empty context 2.692→0.0113 ms,
EIT compute 68.453→64.712 ms, χ arrays bitwise 동일.
late-loaded BLAS, nested/concurrent calls, fork 안전 회귀 필요.

### 8.3 angular-momentum matrix cache — 낮은 우선순위

`angular_momentum_matrices(F)`는 매번 배열 생성 `gabes/zeeman.py:48-66`.
직전 prototype F=2/F=1 pair 51.8→0.4 µs. 전체 compute 영향 작음.
read-only cache + mutation regression이면 안전. `zeeman_manifold()` 자체는 이미 cached
`gabes/zeeman.py:69-70`.

### 8.4 이미 있는 example 최적화 — 재제안 금지

Hanle operating scan은 `cell_mm`가 readout-only인 점을 써서 T/intensity당 atomic raw를
한 번 계산하고 여러 cell length에 재사용 `analysis/squeezing/resonant_hanle_squeezing_reference.py:908-943`.

오늘 4개 cell length `[5,10,25,50] mm` 비교:

- 각 length 직접 solve vs raw reuse: transmission/slope array-exact, bias index 동일.
- median `1.1090 s → 0.2905 s`, **3.82×**.

현재 working tree에 이미 있음. 새 후보 아님. 전용 reuse regression은 아직 없음.

## 9. ⟦AI⟧ 코드·테스트·문서·예제 평가

### 9.1 테스트

강점:

- CG known values, emission normalization, TOC `tests/test_magneto.py:28-63`.
- elastic dephasing의 lifetime/branching/TOC 보존 `tests/test_magneto.py:66-103`.
- transit reload와 collisional dephasing 분리, positivity/trace/Hermiticity
  `tests/test_magneto.py:105-159,352-446`.
- conditional normalization와 density one-use `tests/test_magneto.py:233-302`.
- paraffin dip/peak, wall narrowing, buffer broadening, intrinsic EIA, LCA, NMOR
  `tests/test_magneto.py:305-503`.
- exact Maxwell vs converged fallback, off-path isolation, provenance
  `tests/test_magneto_doppler.py:54-114`.

남은 gap:

- 121→401/801 B-grid/local slope convergence.
- physical-zero hero와 command-zero secondary.
- circular/near-circular NMOR gate.
- open-hyperfine/repump scope status.
- literature sign/trend와 absolute-width 판정 분리.
- clean-checkout exact-Doppler 실행.
- example raw-reuse parity와 provenance hash.

낡은 인자:

- production schema는 `velocity_classes` 제거 `gabes/schemes/magneto.py:285-293`.
- tests/example/config는 계속 전달
  `tests/test_magneto.py:168,421`,
  `tests/verify_hanle_eit_eia.py:80,106,163,187`,
  `analysis/squeezing/resonant_hanle_squeezing_reference.py:370-397,1649-1651`,
  `analysis/squeezing/resonant_hanle_experiment_config.example.json:17`.
- extra dict key라 조용히 무시. 수렴 knob처럼 오해 가능. 제거 필요.

### 9.2 README·User Guide

- README는 exact Maxwell을 설명. current working tree와 일치 `README.md:511-524`.
- User Guide는 Magneto를 `~0.35 s`, `B scan × velocity`로 설명
  `docs/Userguide/GABES_User_Guide_v2.html:778-790`. 현 exact path에 velocity grid 없음.
- Magneto 한계는 실험별 조정 한 줄뿐
  `docs/Userguide/GABES_User_Guide_v2.html:849-854`.
- “절대 스케일, 선폭, NMOR 영점 교차 자동 확인” 문구는 실제 범위보다 강함
  `docs/Userguide/GABES_User_Guide_v2.html:934-937`.
- `docs/Userguide/userguide_assets/hanle.png`는 옛 sharp transmission peak
  `T≈0.99825`. 현재 같은 표면 조건은 `T(0)=0.787789`, crossover.
  Guide caption도 “좁은 Ramsey형 투과 피크”라 stale
  `docs/Userguide/GABES_User_Guide_v2.html:628-639`.

### 9.3 Hanle–FWM example

`analysis/squeezing/resonant_hanle_squeezing_reference.py`가 유일한 실질 batch example.
standalone Magneto-only example 없음.

장점:

- beam power/waist→intensity, transmission slope, balanced detector noise,
  optional measured CSV calibration 연결 `...reference.py:357-420,665-735`.
- cell-length raw reuse 정확·빠름 `...reference.py:908-943`.
- frequency anchor, responsivity, CSV unit/fit 테스트 있음
  `tests/test_resonant_hanle_reference.py:14-141`.

한계:

- obsolete `velocity_classes` 계속 노출.
- transmission slope는 coarse `np.gradient` 3-point `...reference.py:401-413`.
- 기본 `source_anchor=paper`; OE gain/IDS anchor 주입 `...reference.py:1631-1633`.
- saved report는 measured Hanle CSV 없음 표시하지만 ideal/zero-electronics 아래 구체
  pT/√Hz 출력 `...reference.py:1022-1163`.
- report에 생성시각, git commit/dirty, Magneto cache version, source/config hash 없음.

저장 sensitivity 숫자는 Scheme 4 experimental validation 아님. 조건부 system study.

## 10. 검증 기록

- `CLAUDE.md`, `README.md`, Grand Challenge Head Manager memory/dashboard/handoff 확인.
- Scheme registry/UI 순서 실행 확인.
- 기존 Scheme 4 보고서 17개, checklist, DEVLOG, local TODO/issue notes 검색.
- current Magneto/Zeeman/pole code, tests, User Guide, static Hanle image,
  Hanle–FWM source/config/report 확인.
- 5 regime fresh 계산. 모든 pole guard ≤2.97e-12, fallback 0.
- 문헌 비교 스크립트 fresh 실행.
- example 4-cell direct/reuse parity와 3.82× timing 확인.
- focused pytest:
  `python -m pytest -q tests/test_magneto.py tests/test_magneto_doppler.py tests/test_kernels.py tests/test_headless_observables.py tests/test_schemes_render.py tests/test_resonant_hanle_reference.py tests/test_docs_consistency.py`
  → **85 passed in 48.57 s**.
- required root pytest: `python -m pytest -q`
  → **2434 passed, 4 skipped, 1 failed, 5 warnings in 3433.82 s**.
  유일한 실패는 기존 Windows Grand Challenge nested fresh-capsule 경로의
  `FileNotFoundError`:
  `tests/quantum/test_thermal_campaign_extend.py:367`에서 시작해
  `analysis/grand_challenge/thermal_campaign.py:328`의 중첩 destination open 실패.
  Scheme 4 집중 회귀와 무관. 경고 5건은 `gabes_ui/export.py:197`의
  `datetime.utcnow()` deprecation.
- 생산 코드와 Grand Challenge core/state 수정 없음. 오늘 보고서와 automation memory만 갱신.

## 11. ⟦AI⟧ 최종 판정과 순서

Scheme 4는 유용한 실제 원자물리 구현. Zeeman CG, optical pumping, grouped-emission
TOC, transverse-field switching, two-region coated-cell coherence, buffer relaxation,
exact Maxwell average, Beer–Lambert/NMOR readout 있음.

적합: feature 부호, 상대 추세, 실험 설계, 초기 fit seed, 누락 물리 진단.

부적합: absolute magnetometry, publication-grade linewidth/rotation/sensitivity.
원인: EIT/Buffer-LCA label 불일치, NMOR slope 21.7% 미수렴, transit 의미 중복,
physical-zero 혼동, circular plane-rotation 출력, paraffin 폭 7.9–11.5배 불일치,
open-hyperfine/repump·collision·spatial 평균·absolute calibration 부족, HEAD 밖 provenance.

권고 순서:

1. exact-Doppler production/test/evidence clean-checkout 재현성 고정.
2. physical/command zero, local slope, B-grid status, 문헌 trend/absolute 판정 분리.
3. reference 적용 범위와 open-hyperfine/repump status.
4. paraffin exchange/light-dephasing 의미 분리.
5. EIT/Buffer-LCA relabel 또는 근거 있는 retune.
6. circular NMOR gate.
7. User Guide/image/example provenance 갱신.
8. sourced collision coefficient + pressure shift.
9. lazy B-stack. 그 뒤 cached BLAS controller.
10. full VCC, calibrated spatial transport, full open D1/repump, Jones/Stokes는
    dataset·runtime budget 있는 opt-in reference만.

Head Manager closure:

`⟦AI⟧ Scheme 4의 장기 경로에 새 변경이 있는가? → 핵심 source hash와 수치 결과 동일,
새 Grand Challenge dependency 없음 → 기존 봉인 수치 재사용, 재현성·trust gap만 재확인 →
다음 행동은 P0 clean-checkout provenance와 P1 zero/slope/status.`
