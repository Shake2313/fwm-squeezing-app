# 2026-10-16 Scheme 2 — Λ coherence (EIT / AT / CPT) 물리 검토

## 1. 선택·재현 범위

- 로컬 날짜: `2026-10-16`, Asia/Seoul.
- `n = (16 mod 5) + 1 = 2`.
- 선택: **Scheme 2, Λ coherence (EIT / AT / CPT)**.
- 기준: 현재 `main`, `HEAD 39be4ae`.
- 작업 트리에는 검토 전부터 staged/unstaged/untracked 항목이 `1154`개 있음. 기존 변경 건드리지 않음.
- 오늘 생산 코드 변경 없음. 현재 `gabes/schemes/absorption.py` SHA-256은
  `F368A878...F356`; 2026-10-01 manifest의 값과 같음:
  [source manifest](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/analysis/fwm_gain_hotfix/tpd_gain_diagnostic_20261001/source_manifest_before.json:536>).
- 현재 `tests/test_absorption.py` SHA-256은 `D042FC4B...F73`; 같은 manifest와 같음:
  [test manifest](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/analysis/fwm_gain_hotfix/tpd_gain_diagnostic_20261001/source_manifest_before.json:637>).
  따라서 2026-10-11 검토 뒤 Lambda 생산 코드·핵심 시험은 바뀌지 않음.

현재 dropdown 순서. registry의 실제 생성 순서 기준:
[registry](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/schemes/__init__.py:19>),
[README 요약](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/README.md:8>).

| 번호 | key | 현행 정의 | 구현 |
|---:|---|---|---|
| 1 | `sas` | Absorption OD / SAS | `SASScheme` |
| **2** | **`lambda`** | **Λ coherence: EIT / AT / CPT** | `LambdaScheme` |
| 3 | `rydberg_eit` | Rydberg-EIT electrometry | `RydbergEITScheme` |
| 4 | `magneto` | Hanle / EIA / NMOR | `MagnetoScheme` |
| 5 | `fwm` | Squeezing / Biphoton FWM | `FWMScheme` |

`eit`, `at`, `cpt`는 시험·직접 호출 alias. dropdown 별도 scheme 아님:
[aliases](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/schemes/__init__.py:27>).
`ODScheme`도 등록되지 않은 2준위 검증 primitive.

## 2. 기존 개선 제안부터 검색

검색 범위:

- `CLAUDE.md`, `README.md`, `docs/checklist.json`.
- 기존 Scheme 2 일일 보고서 **18개**, 2026-06-26~2026-10-11.
- 저장소의 `TODO`, `FIXME`, `XXX`, `HACK`, issue/proposal/roadmap/DEVLOG 파일·문자열.
- 로컬 DEVLOG 3개. Lambda 관련 별도 메모는 pole 이식 항목 하나:
  [pole reuse DEVLOG](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/analysis/pole_reuse/DEVLOG.md:168>).
- 공개 GitHub Issues API, 2026-10-16 22:18 KST: 응답 `[]`, open/closed issue 0개.
  [repository Issues](https://github.com/Shake2313/fwm-squeezing-app/issues).

결과:

- Lambda 생산 코드 안 TODO/FIXME 없음.
- Lambda 전용 local issue-note 없음.
- **개선안 부재 아님.** checklist와 선행 보고서에 충분히 존재. 새 physics roadmap 추가 불필요.
- 2026-10-11 보고서의 수치·판정도 현재 source hash와 맞음:
  [직전 Scheme 2 보고서](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/docs/daily_report/2026-10-11_scheme-2_lambda-coherence.md:25>).

### 기존 제안의 비용·물리 보존성

| 기존 항목 | 상태 | 계산비 | 물리 보존·변경 판단 |
|---|---|---:|---|
| `lambda-spectrum-validity` | P1 `ready`; [checklist](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/docs/checklist.json:3442>) | 같은 601 scan점·같은 3준위/9×9 solve. 오늘 warm median `37.045→42.157 ms`, **+13.8%** | scan 연산자와 velocity-only Raman 연산자 분리, 중심집중 격자, 보간, resolution/edge/visibility gate. compact OBE 보존. 잘못된 의미와 미수렴 readout 교정. 최우선. |
| `collisional-coefficient-provenance-and-pressure-shift` | P1 `ready`; [checklist](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/docs/checklist.json:3163>) | coefficient lookup·scalar 조립 `O(1)`; solve 차원 동일 | zero-pressure 기본값 보존 가능. sourced pressure shift와 elastic optical dephasing은 실제 line center/폭 변경. runtime보다 종·선·기체별 출처·부호·단위·불확도·validity가 큰 비용. |
| `lambda-pole-residue-doppler` | P3 `ready`, P1 의존; [checklist](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/docs/checklist.json:3870>) | 선행 prototype 약 `1.8×`; 절대 절감 수십 ms | 현재 `dv=2 m/s`, `±4σ` 이산 measure와 row guard를 그대로 쓰면 behavior-preserving 가속. continuous Gaussian/Faddeeva로 바꾸면 measure 변경이라 순수 최적화 아님. P1 뒤 적용. |
| `paraffin-coated-cell-extensions` | P2 `parked`; [checklist](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/docs/checklist.json:3261>) | 새 bright/dark 상태·시간/비행 통계 필요 | OD/SAS population reservoir 복사 불가. Lambda ground coherence, Ramsey phase, wall/return 통계 필요. negligible overhead 아님. checklist의 `schemes/lambda_system.py` 경로는 현재 존재하지 않음. |
| `full-velocity-changing-collision-kernel` | P2 `parked`; [checklist](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/docs/checklist.json:3568>) | velocity classes 결합; 큰 선형계/iterative solve | 현 독립 Maxwell 평균 구조 변경. kernel·검증 spectrum·runtime budget 선택 전 기본 탑재 부적합. |
| `lambda-hyperfine-resolved-manifold` | P2 `parked`; [checklist](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/docs/checklist.json:3590>) | Hilbert/Liouville 차원·매개변수 증가 | 특정 isotope/transition/polarization와 calibration/holdout spectrum 필요. compact 기본값 대신 opt-in reference가 맞음. |

정확한 비공선 2-D 공동 속도 평균도 무료 아님. 공용 기하식은 이미 존재:
[2-D Raman shift](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/doppler.py:106>).
선행 61 scan점 직접합은 현 1-D `3–5 ms` 대비 `1.2–1.3 s`. 수백 배.
기본 경로 즉시 교체보다 opt-in reference와 두 속도축 수렴 gate가 안전.

## 3. 현재 구현이 푸는 물리

가짜 선모양 fit 아님. **정상상태 scalar 3준위 Λ optical Bloch equation**.

- 상태 `g1`, `g2`, `e`. `e→g1`, `e→g2` 각각 `Γ/2`; ground Raman coherence에
  `γgg`: [atoms.lambda3](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/atoms.py:168>).
- probe 고정 `Ωp=10⁻³Γ`: [constant](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/schemes/absorption.py:29>).
- Hamiltonian은 probe, coupling `Ωc`, coupling detuning, two-photon scan 포함:
  [Lambda compute](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/schemes/absorption.py:637>).
- EIT·AT·CPT는 같은 엔진의 세 regime/default:
  [regime defaults](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/schemes/absorption.py:472>).
- coupling power·1/e² diameter는 실험자 Rabi anchor에 `√P/d` 적용:
  [beam helper](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/beam.py:9>),
  [scheme call](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/schemes/absorption.py:587>).
- Doppler-on: 601 scan점 × `dv=2 m/s`, `±4σ` 1-D Maxwell class. affine real-basis batched LU:
  [Doppler path](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/schemes/absorption.py:118>),
  [kernel](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/kernels.py:558>).
- 한 complex susceptibility에서 `α=k Imχ`, `T=exp(-αL)`, `Reχ`, `n_g` 계산:
  [observables](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/observables.py:565>),
  [group index](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/observables.py:616>).
- `cell_mm`, `line_strength`는 navigate-only. 원자 재풀이 없이 readout만 변경:
  [schema](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/schemes/absorption.py:559>),
  [test](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/tests/test_absorption.py:110>).

### 실험가에게 실제로 유용한 범위

- 2광자 공명 EIT dark feature와 정상분산의 정성/반정량 확인.
- **resolved** AT splitting과 `Ωc` 관계.
- coupling power/diameter, `γgg`, opacity 변화 방향.
- 코드 읽기 쉬운 3준위 OBE reference, 교육, 초기 sanity check.

### 모델 경계

- scalar Λ. `Fg,mF`, polarization selection, Zeeman 구조, 비대칭 branching, hyperfine optical pumping, repump 없음.
- `Rb (natural)`은 ⁸⁵Rb/⁸⁷Rb hyperfine susceptibility 합 아님. 우세 isotope의 line constant와
  abundance-weighted scalar density를 한 Λ에 사용:
  [medium construction](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/schemes/absorption.py:434>).
- Rabi는 anchor. CG·편광·mode overlap에서 절대 `Ωc` 예측 아님.
- `buffer_ground_relax_khz`는 출처 있는 pressure law가 아니라 하나의 generic Raman-coherence decay:
  [control](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/schemes/absorption.py:553>).
- species/line/gas-specific optical collision width, pressure shift, Dicke narrowing, VCC 없음.
- CPT는 별도 bichromatic clock-CPT/광펌핑 엔진 아님. 같은 weak-probe Λ의 좁은 preset.
- `group_index()`는 sampled `Reχ`의 local gradient. finite-pulse bandwidth, group delay propagation,
  distortion, detector SNR 없음.

## 4. 오늘 수치 재검증

현재 source로 새 계산. 동일 process에서 Numba warm-up 후 수행.

### 기본 preset

| regime | 중심 투과 | 현재 UI hero | 해석 |
|---|---:|---|---|
| EIT, Doppler on | `0.013558281` | FWHM `0.459680 MHz`, `n_g=100305` | 투과 1.36%. 선폭·기울기 grid 미수렴. |
| AT, Doppler off | `0.989089759` | split `45.968 MHz` | `Ωc/2π=45.968 MHz`와 일치. resolved AT reference로 유용. |
| CPT, Doppler off | `0.706327628` | FWHM `0.923317 MHz`, `n_g=95255` | 보간 FWHM은 `0.805733 MHz`; UI 14.6% 큼. |

기본값은 EIT `(Ωc/2π,γgg/2π,T,L)=(17.238 MHz,57.46 kHz,50 °C,15 mm)`,
AT `(45.968 MHz,57.46 kHz,25 °C,3 mm)`, CPT `(5.746 MHz,28.73 kHz,25 °C,3 mm)`.

### EIT grid 수렴

같은 OBE. scan sample 위치/readout만 비교.

| 계산 | 보간 FWHM | 중심 `n_g` | warm median |
|---|---:|---:|---:|
| 현 uniform 601 | `0.229898639 MHz` | `100305.329` | `37.045 ms` |
| `sinh(5u)` 중심집중 601 | `0.159141369 MHz` | `112707.749` | `42.157 ms` |
| uniform 9,601 reference | `0.159102406 MHz` | `112692.820` | reference용 |

- 현재 `window_fwhm()`은 half-height crossing 보간 없이 바깥 sample 좌표 차 반환:
  [implementation](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/lineshape.py:52>).
- UI `0.459680 MHz`는 9,601점 수렴값의 **2.89배**.
- uniform 601에 보간만 해도 `0.229899 MHz`; 아직 수렴값보다 **44.5% 큼**.
- 중심집중 601은 FWHM `+0.0245%`, `n_g +0.0132%`로 9,601점 기준 재현.
- solve 수·행렬 차원 동일. 이번 PC 시간 `+13.8%`; checklist `≤20%` gate 안.
- 수렴 `n_g`를 15 mm에 단순 환산한 delay scale은 약 `5.64 µs`. 하지만 투과 1.36%,
  finite pulse·distortion·SNR 미계산. slow-light 성능값 아님.

### beam-angle cold-limit 결함

`_affine_scan_coeffs()`는 `A=dL/ds-S_v`로 scan shift와 velocity shift를 묶음:
[affine coefficients](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/schemes/absorption.py:45>).

Doppler off, angle `0→10 mrad`:

```text
max|Δχ| / max|χ₀| = 0.0294203819635
```

`v=0`이면 residual Doppler는 0. 2.942% 변화는 물리 효과 아님.
Optical scan operator와 velocity-only Raman operator 분리 필요. 행렬 크기·scan점·solve 수 증가 없음.

또 `_two_photon_k_ratio()`는 probe와 coupling에 같은 `medium["k_vec"]`를 넘김:
[current ratio](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/schemes/absorption.py:597>).
실제 co-propagating longitudinal `k_p-k_c`를 0으로 버림. 비공선 transverse `v_x`와 axial `v_z`도
현재 한 속도축에 결박. 공용 정확식은 두 성분을 분리:
[noncollinear formula](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/doppler.py:106>).

source 불변이라 선행 2-D 직접합은 재사용함:
[2026-10-11 결과](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/docs/daily_report/2026-10-11_scheme-2_lambda-coherence.md:147>).
1 mrad에서 현 1-D `T(0)=0.034796`, 독립 2-D `0.131721`; **3.79배 차이**.
현재 beam-angle slider는 정량 alignment tolerance reference 아님.

## 5. ⟦AI⟧ 실험 reference 판정

**제한부 합격.** Scheme 2는 유용한 실제 compact 3준위 정상상태 OBE.
resolved AT, Rabi/power/diameter/dephasing trend, 교육·초기 설계의 반정량 reference로 적합.

다음 용도는 현재 부적합:

- 절대 EIT/CPT linewidth·contrast.
- usable slow-light delay·delay-bandwidth·pulse distortion.
- beam-angle 정렬 tolerance.
- clock-CPT, buffer/coated-cell metrology.
- 특정 isotope/hyperfine/Zeeman/polarization spectrum.

핵심 위험은 solver가 아니라 **의미·sampling·claim gate**. 정밀 소수 출력이 현재 모델 신뢰도보다 앞섬.
특히 기본 EIT의 큰 `n_g`를 낮은 투과와 미수렴 grid 없이 단독 인용하면 안 됨.

## 6. 코드·시험·문서·예제 품질

### 코드

- 3준위/9 real-coordinate affine solve. 작고 추적 쉬움.
- complex susceptibility 하나에서 흡수와 분산을 함께 계산. 무관한 curve를 사후 부착하지 않음.
- Numba/NumPy reference parity 존재:
  [kernel parity test](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/tests/test_kernels.py:145>).
- navigate-only opacity/cell remap과 headless path 좋음.
- 단점: 잘못 결합된 scan/velocity operator, uniform grid 의존 hero, 실험 visibility 미평가.

### 시험

현재 잡는 것:

- AT split≈`Ωc`: [test](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/tests/test_absorption.py:135>).
- Rabi `√P/d`: [test](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/tests/test_absorption.py:146>).
- warm-angle broadening trend: [test](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/tests/test_absorption.py:157>).
- cold EIT, sub-natural CPT: [EIT](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/tests/test_absorption.py:176>),
  [CPT](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/tests/test_absorption.py:186>).
- regime defaults·unresolved status: [test](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/tests/test_absorption.py:203>).
- 세 preset과 edge rendering: [render cases](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/tests/test_schemes_render.py:97>).

빈칸:

- Doppler-off angle invariance.
- 실제 `k_p-k_c` 크기·부호.
- 2-D noncollinear joint average.
- clustered/high-resolution EIT/CPT FWHM·`n_g` 회귀.
- samples/FWHM, edge clearance, experimental visibility gate.

현재 warm-angle test는 잘못된 1-D 경로의 상대 trend를 고정. P1에서 cold-limit·analytic RMS·2-D reference로 보강 필요.

### 문서·예제

- README는 scheme 한 줄, architecture, 시험 진입점 제공:
  [scheme table](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/README.md:8>),
  [test entry](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/README.md:163>).
  scalar Λ, 각도 한계, grid/visibility 경계 설명 없음.
- User Guide는 3준위 Λ, 투과, `Reχ`, FWHM, 군지수 설명:
  [guide](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/docs/Userguide/GABES_User_Guide_v2.html:595>).
- “느린 빛의 서명” 문구에 grid 미수렴, 기본 투과 1.36%, finite-pulse 미검증 경고 없음:
  [caption](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/docs/Userguide/GABES_User_Guide_v2.html:603>).
- EIT asset 제목은 `Ωc=6.00 MHz`, `γgg=10.0 kHz`; 현재 EIT 기본은 `17.238 MHz`, `57.46 kHz`.
  AT asset은 현재 기본과 일치. CPT 그림 없음.
- Guide Lambda 성능 `~0.6 s`는 오늘 headless `~0.04 s`와 오래된 차이:
  [timing table](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/docs/Userguide/GABES_User_Guide_v2.html:778>).
- standalone Lambda example/config/CSV 없음. 실제 예제 역할은 시험 파일과 정적 EIT/AT 그림뿐.

재현 가능한 headless sweep 예제 하나면 reference 품질 좋아짐: power/diameter, angle, `γgg`, cell length;
입력 JSON, commit/source hash, grid/visibility status, CSV 출력. 새 physics·추가 solve 요구 아님.

## 7. ⟦AI⟧ 동작을 바꾸지 않는 순수 코드 최적화

1. **BLAS `ThreadpoolController` 재사용.**
   현재 매 solve마다 `threadpool_limits(...)` context 생성:
   [core](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/core.py:23>).
   선행 측정은 기본 EIT `39.207→32.848 ms`, `chi_bar` array-equal, **1.19×**.
   module-level controller 재사용 후보. process/thread safety·중첩 context 회귀 필요.

2. **Affine 행렬 희소 update.**
   현재 모든 `(scan,velocity)`에서 9×9 전체 `base+sA+kvB` 재작성:
   [inner loop](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/kernels.py:571>).
   Lambda `A_coef`, `B_coef` 비영 원소는 각각 81개 중 4개.
   `base+sA`를 velocity loop 밖에서 만들고 B 위치만 갱신 가능. LU·합산 순서·weight 유지, 엄격 parity gate 필요.

3. **Pole-residue port.**
   P1 operator 교정 뒤 현재 이산 Maxwell measure 그대로 port.
   선행 prototype은 상대차 `8.37e-14`, 약 `1.8×`:
   [직전 평가](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/docs/daily_report/2026-10-11_scheme-2_lambda-coherence.md:240>).

4. **낮은 우선순위 immutable cache.**
   `_medium_from_params()`, `atoms.lambda3()`, affine coefficient의 scalar-key cache 가능.
   선행 측정 `_medium_from_params()` `17.25 µs/call`, 전체의 약 `0.05%`:
   [측정 기록](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/docs/daily_report/2026-09-16_scheme-2_lambda-coherence.md:241>).
   mutable array·float-key 오염 위험 대비 이득 작음.

## 8. 검증

집중 회귀:

```text
python -m pytest -q \
  tests/test_absorption.py tests/test_kernels.py \
  tests/test_headless_observables.py tests/test_schemes_render.py \
  tests/test_beam.py tests/test_fwm_angular_doppler.py

64 passed in 41.39 s
```

전체 회귀:

```text
python -m pytest -q

2434 passed, 4 skipped, 1 failed, 5 warnings in 3229.90 s (0:53:49)
```

실패 1건은 Lambda 무관. 이전 Scheme 2 검토와 같은 Windows nested fresh-capsule 경로 문제:

- `tests/quantum/test_thermal_campaign_extend.py::`
  `test_empty_campaign_extension_passes_both_real_capsules_without_any_native_solve`
- 호출 위치: [test](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/tests/quantum/test_thermal_campaign_extend.py:367>).
- nested `.git/grand-challenge-runs/.../.git/grand-challenge-runs/...` 안에서
  `analysis/grand_challenge/reference/ballistic_linear_channel.py`의 부모 경로가 만들어지기 전에
  `destination.open('xb')`가 실행되어 `FileNotFoundError`:
  [extract_bundle](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/analysis/grand_challenge/thermal_campaign.py:328>).
- warning 5건은 `gabes_ui/export.py`의 naive `datetime.utcnow()` deprecation. Lambda 수치와 무관:
  [export timestamp](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes_ui/export.py:197>).

판정: Scheme 2 집중 회귀 `64/64` 통과. 전체 실패는 기존 플랫폼/격리 경로 문제이며 이번 보고서 변경이 만든 회귀 아님.

## 9. ⟦AI⟧ 결론·우선순위

1. `lambda-spectrum-validity`: cold operator 분리, 실제 `k_p-k_c`, 중심집중 601, 보간,
   resolution/edge/visibility status. **거의 같은 비용, 가장 큰 신뢰도 개선.**
2. 정확한 2-D 비공선 평균은 opt-in reference + 직접합 수렴 gate. negligible overhead 아님.
3. sourced collision coefficient·pressure shift 연결. runtime 거의 0; provenance가 본 작업.
4. 순수 속도: cached BLAS controller. P1 뒤 pole port와 희소 affine 조립.
5. Guide·EIT asset·CPT 예제 갱신. 낮은 transmission과 finite-pulse/group-delay 경계 명시.
6. Hyperfine/Zeeman, coated-cell Ramsey, full VCC는 독립 데이터·runtime budget 있을 때만 opt-in research mode.

Grand Challenge 장기 경로 영향: 이번 일일 검토는 production Lambda의 reference 경계 재확인.
현 FWM Grand Challenge 카드·공식 `0/4`·봉인 결과·다음 질문 변경 없음. 새 batch 없음.

`⟦AI⟧ 질문 → 증거 → 장기 효과 → 다음 행동`:
Scheme 2가 실험 reference로 충분한가 → real 3-level OBE는 유용하나 grid와 angle 의미 결함 재현 →
독립 apparatus→absolute FWM 경로는 그대로, cross-scheme 신뢰성 debt만 재확인 → 기존 P1부터 처리.
