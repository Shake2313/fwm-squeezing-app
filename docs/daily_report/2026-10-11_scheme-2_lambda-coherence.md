# 2026-10-11 Scheme 2 — Λ coherence (EIT / AT / CPT) 물리 검토

## 1. 선택·현재 순서

- 로컬 날짜: `2026-10-11`, Asia/Seoul.
- `n = (11 mod 5) + 1 = 2`.
- 선택: **Scheme 2, Λ coherence (EIT / AT / CPT)**.
- 기준: 현재 `main`, `HEAD 39be4ae`.
- 시작 전부터 광범위한 staged/unstaged/untracked 변경 존재. 건드리지 않음.
- Lambda 생산 코드·핵심 시험·Guide 그림은 2026-10-06 검토 뒤 변화 없음. 관련 경로의 이후 commit도 없음.

실제 dropdown 순서. 근거: `gabes/schemes/__init__.py:12-25,37-39`.

| 번호 | key | 표시·물리 | 구현 |
|---:|---|---|---|
| 1 | `sas` | Absorption OD / SAS | `gabes/schemes/sas.py` |
| **2** | **`lambda`** | **Λ coherence: EIT / AT / CPT** | `gabes/schemes/absorption.py:460-800` |
| 3 | `rydberg_eit` | Rydberg-EIT electrometry | `gabes/schemes/rydberg.py` |
| 4 | `magneto` | Hanle / EIA / NMOR | `gabes/schemes/magneto.py` |
| 5 | `fwm` | Squeezing / Biphoton FWM | `gabes/schemes/fwm.py` |

`eit`, `at`, `cpt`는 시험·직접 호출 alias. dropdown 별도 항목 아님:
`gabes/schemes/__init__.py:27-34`. `ODScheme`도 등록되지 않은 2준위 검증 primitive.

## 2. 기존 제안부터 검색

먼저 다음을 검색·대조함.

- `CLAUDE.md`, `README.md`, `docs/checklist.json`.
- 기존 Scheme 2 일일 보고서 **17개**, 2026-06-26~2026-10-06.
- `docs`, `analysis`, `gabes`, `tests`의 `TODO`, `FIXME`, `XXX`, `HACK`,
  issue/proposal/roadmap/DEVLOG 이름·문자열.
- GitHub repository Issues API. open/closed 실제 issue **0건**.

결과:

- Lambda 코드 안 TODO/FIXME 없음.
- Lambda 전용 local issue-note/DEVLOG 없음.
- 개선안 부재 아님. checklist와 기존 일일 보고서에 충분히 존재.
- 새 physics roadmap 항목 추가 불필요. 기존 항목의 비용·보존성 재평가가 맞음.

| 기존 항목 | 상태·근거 | 계산비·물리 보존 판단 |
|---|---|---|
| `lambda-spectrum-validity` | P1 `ready`; `docs/checklist.json:3120-3144` | 최우선. scan/velocity operator 분리, 중심 집중 601점, half-height 보간, 해상도·edge·visibility gate. 행렬 차원과 solve 수 그대로. 오늘 paired median `34.581→34.669 ms`, **+0.25%**. 같은 compact 물리 보존. |
| `collisional-coefficient-provenance-and-pressure-shift` | P1 `ready`; `docs/checklist.json:2841-2878` | sourced pressure shift·elastic optical dephasing은 scalar lookup/조립. solve 차원 동일, 실행비 사실상 0. 진짜 비용은 종·선·기체별 출처, 부호, 단위, 불확도, 검증자료. |
| `lambda-pole-residue-doppler` | P3 `ready`, 위 P1 의존; `docs/checklist.json:3548-3565` | 현 1-D `dv=2 m/s`, `±4σ` measure를 보존하면 순수 가속. operator 의미 먼저 수정해야 함. |
| `paraffin-coated-cell-extensions` | P2 `parked`; `docs/checklist.json:2939-2962` | Lambda는 OD/SAS population reservoir 복사 불가. ground coherence, Ramsey phase, bright/dark flight·return 통계 필요. negligible overhead 아님. checklist area의 `schemes/lambda_system.py`는 존재하지 않음. 실제 파일은 `schemes/absorption.py`. |
| `full-velocity-changing-collision-kernel` | P2 `parked`; `docs/checklist.json:3246-3265` | velocity class 결합. 현 독립 Maxwell 평균 파괴. kernel·dataset·runtime budget 전에는 opt-in reference만 타당. |
| `lambda-hyperfine-resolved-manifold` | P2 `parked`; `docs/checklist.json:3268-3286` | isotope/전이/polarization/held-out spectrum 먼저 필요. Hilbert/Liouville 차원·식별 매개변수 증가. compact 기본값 교체 부적합. |

정확한 비공선 2-D 공동 속도 평균도 물리적으로 필요하지만 무료 아님. 선행 검토의
61 scan점 직접합은 현 1-D `3–5 ms` 대비 `1.2–1.3 s`; 수백 배. 생산 기본 즉시 교체보다
opt-in reference, 두 velocity 축 수렴 gate, 이후 가속이 안전.

## 3. 현재 구현: 무엇을 실제로 푸나

가짜 선모양 fit 아님. **정상상태 3준위 Λ optical Bloch equation**.

- 상태 `g1`, `g2`, `e`. `e→g1`, `e→g2` 각각 `Γ/2`; ground Raman coherence에
  `γgg`: `gabes/atoms.py:168-190`.
- probe 고정 `Ωp=10⁻³Γ`: `gabes/schemes/absorption.py:29-31`.
- coupling `Ωc`, coupling detuning, two-photon detuning을 넣은 Hamiltonian:
  `gabes/schemes/absorption.py:637-667`.
- EIT·AT·CPT는 같은 엔진의 세 regime/default: `gabes/schemes/absorption.py:472-493,574-611`.
- coupling power·1/e² diameter는 실험자 Rabi anchor에 `√P/d` 적용:
  `gabes/beam.py:9-20`, `gabes/schemes/absorption.py:587-595`.
- Doppler-on Lambda: 601 scan점, `dv=2 m/s`, `±4σ` 1-D Maxwell class,
  affine real-basis batched LU: `gabes/schemes/absorption.py:118-140,613-667`,
  `gabes/kernels.py:558-625`.
- 동일 complex susceptibility에서 `α=k Imχ`, `T=exp(-αL)`, `Reχ`, `n_g` 계산:
  `gabes/observables.py:556-581,589-615`.
- `cell_mm`, `line_strength`는 navigate-only. 원자 재풀이 없이 Beer–Lambert/readout만 갱신:
  `gabes/schemes/absorption.py:559-562,681-726`; 회귀 `tests/test_absorption.py:110-132`.

유효 용도:

- 2광자 공명 EIT dark feature·정상분산의 정성/반정량 확인.
- **resolved** AT splitting과 `Ωc` 관계.
- coupling power, diameter, ground-dephasing, opacity 변화 방향.
- 짧은 parameter sanity check, 교육, 초기 실험 설계.

## 4. 실험 reference 경계

- scalar Λ. `Fg,mF`, polarization selection, Zeeman 구조, 비대칭 branching,
  hyperfine optical pumping, repump 없음.
- `Rb (natural)`은 isotope별 susceptibility 합 아님. 우세 isotope의 선 상수와
  abundance-weighted scalar density를 한 Λ에 넣음: `gabes/schemes/absorption.py:422-447`.
- Rabi는 **anchor**. CG·편광·모드 중첩에서 절대 `Ωc`를 예측하지 않음.
- `buffer_ground_relax_khz`는 출처 있는 buffer-pressure law가 아니라 하나의 generic
  Raman-coherence decay: `gabes/schemes/absorption.py:553-557`, `gabes/atoms.py:168-187`.
- Lambda에는 species/line/gas-specific optical collisional width, pressure shift, Dicke
  narrowing, VCC 없음.
- CPT는 별도 bichromatic clock-CPT 엔진 아님. 같은 weak-probe Λ의 좁은 preset:
  `gabes/schemes/absorption.py:482-493,613-635`.
- `group_index()`는 샘플된 `Reχ`의 centered gradient. finite pulse bandwidth, group delay,
  distortion, detector SNR 없음: `gabes/observables.py:607-615`.

**판정:** 유용한 실제 compact 3준위 원자물리. resolved AT와 parameter trend의
반정량 reference로 적합. 절대 EIT/CPT linewidth, slow-light delay, beam-angle tolerance,
clock-CPT, buffer/coated-cell, 특정 hyperfine/polarization reference로 부적합.

## 5. 오늘 재계산

동일 process Numba warm-up 뒤 headless median. 이 PC 참고값.

| 기본 regime | warm median | 중심 투과 | 현재 UI readout |
|---|---:|---:|---|
| EIT, Doppler on | `35.782 ms` | `0.013558281` | FWHM `0.459680 MHz`, `n_g=100305` |
| AT, Doppler off | `3.693 ms` | `0.989089759` | split `45.968 MHz` = `Ωc/2π` |
| CPT, Doppler off | `4.334 ms` | `0.706327628` | FWHM `0.923317 MHz`, `n_g=95255` |

기본값:

- EIT: `(Ωc/2π, γgg/2π, T, L)=(17.238 MHz, 57.46 kHz, 50 °C, 15 mm)`.
- AT: `(45.968 MHz, 57.46 kHz, 25 °C, 3 mm)`.
- CPT: `(5.746 MHz, 28.73 kHz, 25 °C, 3 mm)`.

### 5.1 EIT/CPT 선폭·군지수 미수렴

동일 OBE. scan grid와 readout만 비교.

| 계산 | EIT 보간 FWHM | EIT `n_g` | CPT 보간 FWHM | CPT `n_g` |
|---|---:|---:|---:|---:|
| 현 uniform 601 | `0.229899 MHz` | `100305` | `0.805733 MHz` | `95255` |
| `sinh(5u)` 중심 집중 601 | `0.159141 MHz` | `112708` | `0.803540 MHz` | `95272` |
| uniform 9,601 reference | `0.159102 MHz` | `112693` | `0.803539 MHz` | `95272` |

현재 `window_fwhm()`은 half-height crossing 보간 없이 바깥 sample 좌표 차를 반환:
`gabes/lineshape.py:52-70`. Scheme은 양의 유한 폭이면 hero 표시:
`gabes/schemes/absorption.py:776-800`.

- UI EIT `0.459680 MHz`: 수렴값의 **2.89배**. 중앙 step `0.229840 MHz`가 실제
  FWHM보다 큼. 한 폭 안에 1 sample interval도 안 됨.
- 601점에 선형보간만 적용해도 EIT `0.229899 MHz`. 아직 수렴값보다 44.5% 큼.
- UI CPT `0.923317 MHz`: 수렴값보다 **14.9% 큼**.
- EIT `n_g`: 수렴값보다 약 **11.0% 낮음**.
- 중심 집중 601의 중앙 step `0.01549 MHz`, 약 10.3 samples/FWHM.
- 실행비: uniform/clustered median `34.581/34.669 ms`, **+0.25%**.

따라서 brute-force 9,601점 불필요. 기존 P1의 clustered 601 + 보간 + 해상도/edge/
visibility status가 물리를 바꾸지 않고 문제 해결.

수렴 `n_g`를 15 mm에 단순 환산하면 group delay scale은 약 `5.64 µs`. 하지만 투과가
`1.36%`이고 pulse propagation·bandwidth·distortion·SNR을 풀지 않음. slow-light 성능값으로
인용 금지.

### 5.2 beam-angle 경로의 cold-limit 오류

`_affine_scan_coeffs()`가 `A=dL/ds-S_v`로 scan shift와 velocity shift를 묶음:
`gabes/schemes/absorption.py:45-57`. 오늘 Doppler-off angle `0→10 mrad` 재계산:

```text
max|Δχ| / max|χ₀| = 0.02942038196
```

`v=0`이면 residual Doppler는 0. 2.942% 변화는 물리 효과 아님. Optical scan operator와
velocity-only Raman operator 분리 필요. 행렬 크기·scan점·solve 수 증가 없음.

또한 `_two_photon_k_ratio()`는 probe와 coupling에 같은 `medium["k_vec"]`를 넘김:
`gabes/schemes/absorption.py:597-601`. 실제 co-propagating longitudinal
`k_p-k_c`를 0으로 버림. 비공선일 때도 `v_z`, `v_x` 독립인데 현재 한 velocity 축에 결박:
`gabes/atoms.py:95-109,168-189`; 정확한 공용 기하식은 `gabes/doppler.py:106-172`.

선행 동일-core 2-D 직접 평균의 50 °C, 3 mm 공명 결과:

| angle | 현 1-D `T(0)` | 독립 2-D `T(0)` |
|---:|---:|---:|
| 0 mrad | `0.423098` | `0.423098` |
| 1 mrad | `0.034796` | `0.131721` |
| 10 mrad | `0.002679` | `0.004026` |

1 mrad에서 3.79배 차이. 정량 alignment reference 아님. operator cold-limit 수정은 무료;
정확한 2-D 평균은 무료 아님. 둘을 한 개선으로 뭉개면 안 됨.

## 6. 코드·문서·시험·예제 품질

### 코드

- 3준위/9 real-coordinate affine solve. 작고 읽기 쉬움. Numba와 NumPy reference parity 있음.
- complex susceptibility 한 개에서 흡수와 분산을 같이 냄. 사후에 무관한 curve를 붙이지 않음.
- navigate-only opacity/cell remap과 headless path 좋음.
- 문제는 solver보다 의미·sampling·claim gate. 정확한 소수점 출력이 모델 신뢰도보다 앞섬.

### 시험

현재 잡는 것:

- AT split≈`Ωc`: `tests/test_absorption.py:135-143`.
- Rabi `√P/d`: `tests/test_absorption.py:146-154`.
- warm-angle broadening trend: `tests/test_absorption.py:157-173`.
- cold EIT, sub-natural CPT: `tests/test_absorption.py:176-200`.
- regime defaults·unresolved status: `tests/test_absorption.py:203-241`.
- Numba/NumPy affine parity: `tests/test_kernels.py:145-179`.
- 세 preset과 edge rendering: `tests/test_schemes_render.py:97-107`.

빈칸:

- Doppler-off angle invariance.
- 실제 `k_p-k_c`와 부호.
- 2-D non-collinear joint average.
- clustered/high-resolution FWHM·`n_g` 회귀.
- samples/FWHM, edge clearance, experimental visibility gate.

현재 warm-angle test는 잘못된 1-D 경로의 상대 trend를 고정. P1 때 대체/보강 필요.

### 문서·예제

- README는 기능 한 줄과 시험 진입점만 제공: `README.md:8-16,67-75,163-180`.
  scalar Λ, 각도, 선폭, low-transmission slow-light 경계 없음.
- User Guide는 3준위 Λ, 투과, `Reχ`, FWHM, 군지수 설명:
  `docs/Userguide/GABES_User_Guide_v2.html:519-549,595-613`.
- “느린 빛의 서명” 문구에 grid 미수렴, 기본 EIT 투과 1.36%, finite-pulse 미검증 경고 없음:
  `docs/Userguide/GABES_User_Guide_v2.html:603-606`.
- EIT 그림 제목은 `Ωc=6.00 MHz`, `γgg=10.0 kHz`; 현재 기본은
  `17.238 MHz`, `57.46 kHz`. AT 그림은 현재 기본과 일치. CPT 그림 없음.
  assets: `docs/Userguide/userguide_assets/eit.png`, `at.png`.
- Guide Lambda 성능 `~0.6 s`: `docs/Userguide/GABES_User_Guide_v2.html:778-784`.
  오늘 headless solve `35.8 ms`; 오래된 수치. 그림/하드웨어 조건 manifest도 없음.
- 독립 `examples/` 또는 Lambda script/config/CSV 예제 없음. 실제 예제 역할은 위 시험 파일.

실험 reference 품질에는 headless power/angle/dephasing sweep 예제 하나가 유용. 입력 JSON,
commit, grid/visibility status, CSV 출력 포함. 새 physics·추가 solve 없이 재현성만 개선.

## 7. 동작을 바꾸지 않는 순수 코드 최적화

1. **BLAS `ThreadpoolController` 재사용.**
   - 현 `core.blas_single_thread()`는 매번 `threadpool_limits(...)` 생성:
     `gabes/core.py:23-34`.
   - 오늘 EIT median `39.207→32.848 ms`, **1.19×**. `chi_bar`
     `np.array_equal=True`, 최대차 0.
   - 공용 저위험 후보. process/thread safety·중첩 context 회귀 필요.

2. **Affine 행렬 조립 희소 update.**
   - `_affine_scan_chi_real()`은 모든 `(scan,velocity)`에서 9×9 전체를 다시 씀:
     `gabes/kernels.py:571-595`.
   - 기본 Lambda의 `A_coef`, `B_coef` 비영 원소는 각각 81개 중 4개.
   - `base+sA`를 velocity loop 밖에서 만들고 B 위치만 갱신 가능. 동일 LU, 순서,
     weight 유지. 엄격 parity 필요.

3. **기존 pole-residue port.**
   - checklist `docs/checklist.json:3548-3565`.
   - 선행 nonpersistent prototype은 현 601×711 이산 measure 상대차 `8.37e-14`,
     약 `1.8×`. 현 1-D 행동의 순수 가속.
   - `lambda-spectrum-validity`의 operator 의미 수정 뒤 적용. 2-D 물리 교정 아님.

4. **낮은 우선순위 cache.**
   - `_medium_from_params()`, `atoms.lambda3()`, affine coefficient의 scalar-key cache 가능.
   - 현재 35 ms 중 작은 몫. mutable array·float-key 오염 방지 필요.

## 8. 검증

관련 회귀:

```text
python -m pytest -q \
  tests/test_absorption.py tests/test_kernels.py \
  tests/test_headless_observables.py tests/test_schemes_render.py \
  tests/test_beam.py tests/test_fwm_angular_doppler.py

64 passed in 12.08 s
```

AGENTS.md 요구 전체 회귀:

```text
python -m pytest -q

2202 passed, 4 skipped, 2 failed in 4828.96 s (1:20:28)
```

실패 2건. 둘 다 Lambda 무관. 직전 일일 검토와 동일한 기존 실패.

1. `tests/quantum/test_thermal_campaign_extend.py::`
   `test_empty_campaign_extension_passes_both_real_capsules_without_any_native_solve`
   - Windows nested fresh-capsule 경로에서
     `analysis/grand_challenge/reference/ballistic_linear_channel.py` 생성 전
     `FileNotFoundError`.
2. `tests/test_fwm_excess_noise.py::`
   `test_hidden_noise_sliders_retain_values_across_streamlit_reruns`
   - ScrubField 미등록 fallback 환경에서 gain-proportional mode가
     `0.02×(3−1)=0.04` 대신 이전 constant noise `2.0` 유지.

관련 Lambda 집중 회귀는 64/64 통과. 오늘 보고서 외 생산 코드 변경 없음.

## 9. 결론·우선순위

Scheme 2는 **실제 정상상태 3준위 Λ OBE를 구현한 유용한 반정량 실험 계획 reference**.
resolved AT, Rabi/power/diameter/dephasing trend에 좋음.

현재 EIT/CPT linewidth와 group index는 기본 grid 미수렴. beam-angle 경로는 cold-limit와
joint-velocity geometry 위반. 절대 linewidth, slow-light delay, 정렬 tolerance로 인용 금지.

권장 순서:

1. 기존 P1: cold operator 분리, 중심 집중 601점, 보간, resolution/edge/visibility status.
   거의 같은 비용.
2. 실제 coupling wavevector와 2-D joint Maxwell은 opt-in reference + 직접합 수렴 gate부터.
   무시 가능한 비용 아님.
3. sourced collision coefficient·pressure shift 연결. 적용 runtime 거의 0; provenance가 본 작업.
4. 순수 속도: cached BLAS controller. P1 뒤 pole port와 희소 affine 조립.
5. Guide·EIT asset·CPT 예제 갱신. 현재 기본값과 claim boundary 명시.
6. Hyperfine/Zeeman, coated-cell Ramsey, full VCC는 독립 데이터·runtime budget이 생길 때만
   opt-in research mode.
