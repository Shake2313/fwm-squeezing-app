# 2026-10-06 Scheme 2 — Λ coherence (EIT / AT / CPT) 물리 검토

## 1. 선택, 범위, 현재 순서

- 로컬 날짜: `2026-10-06` (Asia/Seoul).
- `n = (6 mod 5) + 1 = 2`.
- 선택: **Scheme 2, Λ coherence (EIT / AT / CPT)**.
- 기준: 현재 `main`, `HEAD fcb312d`.
- 검토 시작 전부터 광범위한 staged/unstaged 변경 존재. 기존 변경 건드리지 않음.
- Lambda 생산 코드·핵심 시험은 HEAD와 작업 트리 모두 2026-09-26 검토 뒤 변화 없음. `gabes/schemes/absorption.py` 마지막 관련 커밋은 `d98bdf4`.

실제 dropdown 순서. 근거: `gabes/schemes/__init__.py:19-25,37-39`.

| 번호 | key | 표시·물리 | 구현 |
|---:|---|---|---|
| 1 | `sas` | OD / SAS | `gabes/schemes/sas.py` |
| **2** | **`lambda`** | **Λ coherence: EIT / AT / CPT** | `gabes/schemes/absorption.py:466-800` |
| 3 | `rydberg_eit` | Rydberg-EIT electrometry | `gabes/schemes/rydberg.py` |
| 4 | `magneto` | Hanle / EIA / NMOR | `gabes/schemes/magneto.py` |
| 5 | `fwm` | Squeezing / Biphoton FWM | `gabes/schemes/fwm.py` |

`eit`, `at`, `cpt`는 직접 호출·시험용 alias. dropdown 별도 항목 아님: `gabes/schemes/__init__.py:27-34`. `ODScheme`도 등록되지 않은 2준위 검증 primitive.

## 2. 기존 제안, TODO, issue 먼저 검색

검색 범위:

- `CLAUDE.md`, `README.md`, `docs/checklist.json`.
- 기존 Scheme 2 일일 보고서 16개: 2026-06-26부터 2026-09-26.
- `docs`, `analysis`, `gabes`, `tests`의 `TODO`, `FIXME`, `XXX`, `HACK`, `issue`, `proposal`.
- TODO/issue/note 이름의 로컬 파일.
- 공개 GitHub Issues API: 2026-10-06 open/closed 합계 **0건**.

결과:

- Lambda 코드 안 별도 TODO/FIXME 없음.
- 별도 Lambda local issue-note 파일 없음.
- 개선안 부재 아님. 구체적 physics·신뢰성·성능 제안은 checklist와 선행 보고서에 이미 충분히 존재.
- 오늘 새 physics roadmap 항목 만들지 않음. 기존 항목의 비용과 물리 보존 여부 재평가.

| 기존 항목 | 상태·근거 | 계산비·물리 보존 판단 |
|---|---|---|
| `lambda-spectrum-validity` | P1 `ready`; `docs/checklist.json:3118-3142` | 최우선. scan/velocity operator 분리, 중심 집중 601점, half-height 보간, resolution/edge/visibility gate는 행렬 차원·solve 수 그대로. 거의 무료. 정확한 비공선 2-D 공동 속도 평균은 별도. 무료 아님. |
| `collisional-coefficient-provenance-and-pressure-shift` | P1 `ready`; `docs/checklist.json:2839-2876` | pressure shift·optical dephasing 적용 자체는 scalar 조립. 런타임 증가 무시 가능. 큰 비용은 종/선/기체별 출처·부호·단위·불확도·holdout 데이터. |
| `lambda-pole-residue-doppler` | P3 `ready`, P1 의존; `docs/checklist.json:3546-3563` | 현 1-D `dv=2 m/s`, `±4σ` 이산 Maxwell measure를 유지하면 순수 가속. P1에서 operator 의미 고정 뒤 적용. |
| `paraffin-coated-cell-extensions` | P2 `parked`; `docs/checklist.json:2937-2959` | Λ는 population reservoir가 아니라 ground coherence, Ramsey phase, dark-flight/return 통계 필요. 새 상태·적분·실험 자료 필요. negligible overhead 아님. |
| `full-velocity-changing-collision-kernel` | P2 `parked`; `docs/checklist.json:3244-3263` | velocity class 결합. 현 독립 평균 파괴. 선택한 kernel·dataset·runtime budget 전에는 opt-in reference만 타당. |
| `lambda-hyperfine-resolved-manifold` | P2 `parked`; `docs/checklist.json:3266-3284` | isotope/transition/polarization/검증 spectrum 필요. Liouville 차원과 식별 매개변수 증가. compact 기본값 교체 부적합. |

## 3. 현재 구현: 실제 물리인가

### 3.1 실제로 푸는 것

가짜 fit curve 아님. 실제 compact 정상상태 OBE.

- `atoms.lambda3()`는 `g1`, `g2`, `e` 3준위 Lindblad 모델. `e→g1`, `e→g2` 각각 `Γ/2`; ground Raman coherence에 `γgg`: `gabes/atoms.py:168-190`.
- probe 고정 `Ωp = 10⁻³Γ`: `gabes/schemes/absorption.py:29-31`.
- coupling `Ωc`, two-photon detuning을 넣은 Hamiltonian 뒤 정상상태 밀도행렬 풀이: `gabes/schemes/absorption.py:637-667`.
- 동일 엔진의 세 regime: EIT, 강결합 AT, 좁은 CPT scan. 기본값: `gabes/schemes/absorption.py:486-493,574-585`.
- Doppler-on Lambda는 601 scan점 × `dv=2 m/s`, `±4σ` 1-D Maxwell class: `gabes/schemes/absorption.py:118-140,613-635`.
- readout: `χ̄ → χphys → α=k Imχ → T=exp(-αL)`. `Reχ`도 같은 susceptibility에서 나옴: `gabes/observables.py:556-581`.
- coupling power·1/e² diameter는 측정자가 정한 Rabi anchor에 `sqrt(P)/d` 적용: `gabes/beam.py:9-20`, `gabes/schemes/absorption.py:587-595`.

따라서 다음 용도 유효.

- two-photon resonance EIT 암흑선·정상분산의 qualitative/semi-quantitative 확인.
- resolved AT splitting과 `Ωc` 관계.
- power·beam diameter·ground dephasing·opacity 변화 방향성.
- 짧은 시간 안 parameter sanity check, 교육·설계 초기값 탐색.

### 3.2 실험 reference 경계

- scalar Λ. `Fg,mF`, polarization selection, Zeeman 구조, 비대칭 branching, hyperfine optical pumping, repump 없음.
- Rb/Cs D1/D2 선택은 선폭·파장·질량·밀도·dipole을 바꿈. 하지만 `Rb (natural)`은 동위원소별 susceptibility 합이 아님. 우세 isotope 선 상수와 abundance-weighted scalar density를 한 Λ에 넣음: `gabes/schemes/absorption.py:422-447`.
- Rabi knob는 **anchor**. CG·편광·모드 중첩으로 절대 `Ωc`를 예측하지 않음.
- `buffer_ground_relax_khz` 한 항이 transit, spin destruction, diffusion, VCC를 구분하지 않음: `gabes/schemes/absorption.py:553-557`.
- pressure shift와 sourced optical collisional dephasing 없음.
- CPT는 별도 bichromatic clock-CPT 엔진 아님. 같은 weak-probe Λ의 좁은 preset: `gabes/schemes/absorption.py:482-493,613-635`.
- `group_index()`는 샘플된 `Reχ`의 수치 gradient: `gabes/observables.py:607-615`. pulse bandwidth, group delay, pulse distortion, 검출 SNR 없음.

**실험가 판정:** 유용한 실제 3준위 원자물리. resolved AT와 knob trend의 반정량 reference로 적합. 절대 EIT/CPT linewidth·slow-light delay·clock-CPT·buffer/coated-cell·특정 hyperfine/polarization·beam-angle tolerance의 정량 기준으로 부적합.

## 4. 현재 작업 트리 재계산

동일 프로세스 Numba warm-up 뒤 headless 계산. 시간은 이 PC 참고값.

| 기본 regime | warm median | 중심 투과 | 현재 UI 결과 |
|---|---:|---:|---|
| EIT, Doppler on | `33.74 ms` | `0.013558281` | FWHM `0.459680 MHz`, `n_g=100305` |
| AT, Doppler off | `3.25 ms` | `0.989089759` | split `45.968 MHz` = `Ωc/2π` |
| CPT, Doppler off | `3.28 ms` | `0.706327628` | FWHM `0.923317 MHz`, `n_g=95255` |

기본값:

- EIT: `(Ωc/2π, γgg/2π, T, L) = (17.238 MHz, 57.46 kHz, 50 °C, 15 mm)`.
- AT: `(45.968 MHz, 57.46 kHz, 25 °C, 3 mm)`.
- CPT: `(5.746 MHz, 28.73 kHz, 25 °C, 3 mm)`.

### 4.1 EIT/CPT 선폭·군지수 미수렴

| 계산 | EIT 보간 FWHM | EIT `n_g` | CPT 보간 FWHM | CPT `n_g` |
|---|---:|---:|---:|---:|
| 현 uniform 601점 | `0.229899 MHz` | `100305` | `0.805733 MHz` | `95255` |
| `sinh(5u)` 중심 집중 601점 | `0.159141 MHz` | `112708` | `0.803540 MHz` | `95272` |
| uniform 9,601점 reference | `0.159102 MHz` | `112693` | `0.803539 MHz` | `95272` |

현재 `window_fwhm()`은 half-height crossing 보간 없이 바깥 sample 좌표를 뺌: `gabes/lineshape.py:52-70`. Scheme은 양의 유한 폭이면 hero로 표시: `gabes/schemes/absorption.py:776-799`.

- UI EIT `0.459680 MHz`는 수렴 보간값의 **2.89배**.
- UI CPT `0.923317 MHz`는 수렴값보다 **14.9% 큼**.
- EIT 중심 step `0.229840 MHz`; 실제 수렴 FWHM `0.159102 MHz`보다 큼.
- EIT `n_g`는 수렴값보다 약 **11.0% 낮음**.
- 같은 601점을 중심 집중했을 때 EIT 단일 paired 실행 `32.32→32.90 ms` (`+1.8%`, 잡음 수준). 9,601점은 `474 ms`.

결론: 기존 P1의 clustered grid + 보간 + resolution/edge/visibility gate는 물리 모델과 solve 수 보존. 거의 무료. 9,601점 brute force 필요 없음.

### 4.2 각도 경로: 두 독립 문제

첫째. cold-limit 오류.

`_affine_scan_coeffs()`가 `A = dL/ds - S_v`로 scan shift와 velocity shift를 묶음: `gabes/schemes/absorption.py:45-57`. Doppler-off인데 angle `0→10 mrad`가

`max|Δχ| / max|χ₀| = 0.02942038`

변화. `v=0`이면 residual Doppler는 0이어야 함. optical scan operator와 velocity-only Raman operator 분리 필요. 행렬 크기·solve 수 증가 없음.

둘째. non-collinear joint velocity 오류.

현재 `_two_photon_k_ratio()`는 `k_probe = k_coupling`을 넣어 `|Δk|/k`만 계산: `gabes/schemes/absorption.py:597-601`, `gabes/beam.py:42-55`. 이후 optical Doppler와 Raman Doppler를 같은 한 축 속도에 결박: `gabes/atoms.py:95-109,168-189`.

실제 기하:

```text
optical Doppler = k_p v_z
Raman Doppler = (k_p - k_c cosθ) v_z - k_c sinθ v_x
```

`v_x`, `v_z` 독립. 공용 정확식 이미 존재: `gabes/doppler.py:106-172`. Lambda 생산 경로는 미사용.

50 °C, 3 mm, 동일 compact OBE 공명점:

| 각도 | 현 1-D `T(0)` | 독립 2-D Maxwell `T(0)` | 2-D / 1-D |
|---:|---:|---:|---:|
| 0 mrad | `0.423098` | `0.423098` | `1.000` |
| 1 mrad | `0.034796` | **`0.131721`** | **`3.79`** |
| 10 mrad | `0.002679` | `0.004026` | `1.50` |

2-D 확인: `v_z` 711 × `v_x` 357, center 한 점 약 `0.30 s`. 현 전체 601점 1-D EIT 약 `0.034 s`. 정확한 2-D 공동 평균은 negligible overhead 아님. 생산 기본 즉시 교체보다 opt-in reference, 축별 수렴 gate, 이후 pole/적분 가속이 안전.

## 5. 코드·문서·시험·예제 품질

### 코드

- 작은 9×9 실수 Hermitian-basis solve. 구조 명료. Numba affine kernel과 headless path 있음.
- `cell_mm`, `line_strength` navigate-only: 재풀이 없이 Beer–Lambert/readout 갱신: `gabes/schemes/absorption.py:559-562,681-726`.
- 물리 한계보다 UI precision이 앞섬. `buffer ground relaxation`, angle, FWHM, group index가 정확한 실험량처럼 보일 위험.

### 시험

현재 시험이 잡는 것:

- AT split ≈ `Ωc`: `tests/test_absorption.py:135-143`.
- Rabi `sqrt(P)/d`: `tests/test_absorption.py:146-154`.
- warm-angle broadening trend: `tests/test_absorption.py:157-173`.
- cold EIT, sub-natural CPT: `tests/test_absorption.py:176-200`.
- regime defaults·unresolved status: `tests/test_absorption.py:203-241`.
- Numba/NumPy affine parity: `tests/test_kernels.py:145-179`.
- 세 preset·edge rendering: `tests/test_schemes_render.py:97-107`.

빈칸:

- Doppler-off angle invariance.
- 실제 `k_p-k_c`와 부호.
- 2-D non-collinear joint average.
- clustered/high-resolution FWHM·`n_g` 비교.
- samples/FWHM, edge clearance, experimental visibility gate.

현재 warm-angle test는 잘못된 1-D 경로의 상대 trend를 회귀로 고정. P1 수정 때 교체 필요.

### 문서·예제

- README 기능 요약과 시험 진입점은 맞음: `README.md:13,68-72,143`.
- scalar Λ, angle, linewidth 신뢰 경계는 README에 없음.
- User Guide는 3준위 Λ, 투과, `Reχ`, FWHM, 군지수 설명: `docs/Userguide/GABES_User_Guide_v2.html:531-534,595-613`.
- Guide의 “느린 빛의 서명”은 기본 `T=1.36%`, 미수렴 `n_g`, pulse 미검증 경고 없음: `:604-606`.
- `eit.png` 제목은 `Ωc=6.00 MHz`, `γgg=10.0 kHz`. 현재 기본 `17.238 MHz`, `57.46 kHz`와 불일치. `at.png`는 현재 기본과 일치. CPT 그림 없음.
- 그림 입력 manifest·생성 script·commit 기록 없음.
- Guide Lambda 성능 `~0.6 s`: `:778-784`. 현재 warm median `0.0337 s`; stale.
- 독립 `examples/` 디렉터리와 Lambda script/config/CSV 예제 없음. 실행 예제 역할을 시험 파일이 대신함.

실험 reference 품질 개선에는 headless power/angle/dephasing sweep 예제 하나 필요. 입력 JSON, commit, grid/visibility status, CSV 출력 포함. 새 physics 없이 재현성만 올림. 기존 P1 문서화 범위로 처리 가능.

## 6. 순수 코드 최적화

물리·현재 이산 평균·출력 행동을 바꾸지 않는 후보만 분리.

1. **BLAS thread controller 재사용 — 우선순위 높음.**
   - `core.blas_single_thread()`는 호출마다 `threadpool_limits(...)` 생성: `gabes/core.py:23-34`.
   - 오늘 empty context median `2.438→0.012 ms`.
   - 기본 EIT median `33.266→30.711 ms`: **7.7% 단축**.
   - `chi_bar` `np.array_equal=True`, 최대 절대차 `0`.
   - 공용 저위험 최적화. controller를 module-level로 재사용하고 process/thread safety 회귀 필요.

2. **Affine 행렬 조립 희소화 — 낮은 위험, 중간 가치.**
   - 기본 EIT real-basis `A_coef`, `B_coef`: 각각 81개 중 nonzero 4개.
   - 현재 모든 `(scan, velocity)`에서 81개 원소 재작성: `gabes/kernels.py:571-595`.
   - `base+sA`를 velocity loop 밖에서 만들고 B의 4개 위치만 갱신 가능.
   - 동일 LU, 동일 velocity 순서, 동일 weighted sum 유지. bitwise 또는 엄격 parity test 필요.

3. **1-D pole-residue port — P1 뒤 적용.**
   - checklist: `docs/checklist.json:3546-3563`.
   - `ShiftInvertSystem`: `gabes/pole_doppler.py:253-337`.
   - 2026-09-26 비영구 prototype: 현 601×711 measure 상대차 `8.37e-14`, median 약 `1.8×`.
   - 현 1-D 동작의 순수 가속. 2-D angle physics 교정 아님.
   - operator 의미가 잘못된 상태로 먼저 최적화하면 오류를 빠르게 만들 뿐. P1 의존 타당.

4. **낮은 우선순위.**
   - `_medium_from_params()`·affine coefficient scalar cache 가능. 총비용 작음.
   - 현재 navigate-only remap과 `headless_observables()` 구조 유지가 더 중요.

## 7. 검증

관련 회귀:

```text
python -m pytest -q \
  tests/test_absorption.py tests/test_kernels.py \
  tests/test_headless_observables.py tests/test_schemes_render.py \
  tests/test_beam.py tests/test_fwm_angular_doppler.py

64 passed in 12.07 s
```

AGENTS.md 요구 전체 시험:

```text
python -m pytest -q

collection 중단: 1 error in 24.29 s
references/AutoOD/test_autood_gui_adjustments.py
→ references/AutoOD/rb_od_calculator_optimized_final.py:5
→ ModuleNotFoundError: No module named 'sympy'
```

Lambda 회귀 실패 아님. optional AutoOD reference dependency 문제. 생산 파일 변경 없음.

## 8. 결론·우선순위

Scheme 2는 **실제 정상상태 3준위 Λ OBE를 구현한 유용한 반정량 실험 계획 reference**. resolved AT, Rabi/power/beam-size/dephasing trend에 좋음.

하지만 현재 EIT/CPT linewidth와 group index는 기본 격자 미수렴. beam-angle 경로는 cold-limit와 joint-velocity geometry 둘 다 위반. 이 숫자를 절대 linewidth, slow-light delay, alignment tolerance로 인용 금지.

권장 순서:

1. 기존 P1: cold operator 분리, 중심 집중 601점, 보간, resolution/edge/visibility status. 거의 같은 비용.
2. 실제 coupling wavevector와 2-D joint Maxwell 평균은 opt-in reference + 직접합 수렴 gate부터. 무시 가능한 비용 아님.
3. sourced collision coefficient·pressure shift 연결. 적용 runtime 거의 0, provenance가 본 작업.
4. 순수 속도: cached BLAS controller 먼저. P1 뒤 `ShiftInvertSystem` pole port.
5. Hyperfine/Zeeman, coated-cell Ramsey, full VCC는 독립 데이터·runtime budget 있을 때만 research mode.

