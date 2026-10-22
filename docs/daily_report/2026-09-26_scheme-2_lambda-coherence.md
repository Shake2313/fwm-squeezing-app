# 2026-09-26 Scheme 2 — Λ coherence (EIT / AT / CPT) 물리 검토

## 선택과 현재 스킴 순서

- 로컬 날짜: `2026-09-26` (Asia/Seoul).
- `n = (26 mod 5) + 1 = 2`.
- 기준 작업본: `main`, `HEAD 5b29c4f`. 검토 전부터 광범위한 staged/unstaged 변경이 있었으며 이 보고서는 이를 수정하지 않았다.
- 실제 dropdown 순서는 `gabes/schemes/__init__.py:19-25,37-39`의 `_SCHEMES`와 `all_schemes()`를 런타임에서 다시 확인했다.

| 번호 | key | 현재 표시·물리 | 구현 |
|---:|---|---|---|
| 1 | `sas` | OD / SAS | `gabes/schemes/sas.py` |
| **2** | **`lambda`** | **Λ coherence: EIT / AT / CPT** | `gabes/schemes/absorption.py:466-800` |
| 3 | `rydberg_eit` | Rydberg-EIT electrometry | `gabes/schemes/rydberg.py` |
| 4 | `magneto` | Hanle / EIA / NMOR | `gabes/schemes/magneto.py` |
| 5 | `fwm` | Squeezing / Biphoton FWM | `gabes/schemes/fwm.py` |

`eit`, `at`, `cpt`는 직접 호출·시험용 alias이고 dropdown의 별도 항목은 아니다 (`gabes/schemes/__init__.py:27-34`). `ODScheme`도 등록되지 않은 2준위 검증 primitive다.

## 기존 제안·TODO·issue 기록을 먼저 확인한 결과

`CLAUDE.md`, `README.md`, `docs/checklist.json`, 2026-06-26부터 2026-09-21까지의 Scheme 2 일일 보고서, `docs/`, `analysis/`, `gabes/`, `tests/`의 TODO/FIXME/issue/proposal 문자열과 해당 이름의 파일을 검색했다. Lambda 코드 안의 별도 TODO/FIXME나 로컬 issue-note 파일은 없지만, **구체적인 개선 제안은 이미 충분히 존재한다.** 따라서 새 항목을 발명하지 않고 현재 코드에 대해 비용과 우선순위를 다시 평가했다.

| 기존 항목 | 현재 상태 | 오늘의 비용·물리 보존 판단 |
|---|---|---|
| `lambda-spectrum-validity` | P1 `ready`, `docs/checklist.json:2909-2933` | 최우선. cold scan/velocity operator 분리, 중심 집중 601점, half-height 보간과 해상도 gate는 행렬 차원·solve 수를 늘리지 않아 거의 무료다. 정확한 비공선 공동 속도 평균은 별도이며 무료가 아니다. |
| `collisional-coefficient-provenance-and-pressure-shift` | P1 `ready`, `:2630-2667` | pressure shift와 optical dephasing 적용은 scalar 조립이라 런타임이 거의 늘지 않는다. 큰 비용은 종·선·기체별 계수의 출처, 부호, 단위, 불확도와 검증 데이터다. |
| `lambda-pole-residue-doppler` | P3 `ready`, `:3337-3354` | 현 1차원 2 m/s·±4σ 이산 Maxwell measure를 보존하는 순수 가속으로 가능하다. 다만 P1의 operator 의미를 먼저 고정해야 한다. |
| `paraffin-coated-cell-extensions` | P2 `parked`, `:2728-2751` | Λ에는 population reservoir가 아니라 ground coherence, Ramsey phase, dark-flight/return 통계가 필요하다. 데이터·상태공간이 추가되어 negligible overhead가 아니다. |
| `full-velocity-changing-collision-kernel` | P2 `parked`, `:3035-3054` | velocity class를 결합해 현재 독립 평균을 깨므로 무겁다. 특정 kernel·dataset·runtime budget이 생길 때 opt-in reference로만 타당하다. |
| `lambda-hyperfine-resolved-manifold` | P2 `parked`, `:3057-3075` | 특정 isotope/transition/polarization과 holdout spectrum을 먼저 정해야 한다. Liouville 차원과 식별할 매개변수가 커지므로 compact default 교체에는 부적합하다. |

직전 Scheme 2 검토(2026-09-21) 뒤 `absorption.py`, `atoms.py`, `beam.py`, `test_absorption.py`에는 committed/working-tree physics 변경이 없다. 이 묶음의 마지막 관련 커밋은 `d98bdf4`(2026-09-05)다.

## 현재 구현이 실제로 푸는 물리

### 유효한 핵심

`atoms.lambda3()`는 `g1`, `g2`, `e`의 3준위 Lindblad 모델이다. `e→g1/g2`에 각각 `Γ/2`의 자발붕괴를 두고, ground Raman coherence에 `gamma_gg`를 적용한다 (`gabes/atoms.py:168-190`). Scheme은 고정 weak probe `Ωp=10⁻³Γ`, coupling `Ωc`, two-photon detuning을 가진 정상상태 Liouvillian을 푼다 (`gabes/schemes/absorption.py:29-31,637-667`). 따라서 다음은 임의의 fit curve가 아니라 실제 compact OBE 응답이다.

- two-photon resonance의 EIT 암흑선과 정상분산,
- 강결합에서의 resolved Autler–Townes doublet,
- ground-coherence decay가 정하는 dark-resonance floor,
- `χ → α=k Imχ → T=exp(-αL)` 및 `Reχ`의 분산 readout (`gabes/observables.py:556-615`).

Coupling power와 1/e² diameter는 측정자가 정한 Rabi anchor에 `sqrt(P)/d`를 곱한다 (`gabes/beam.py:9-20`; `absorption.py:587-595`). 실험 knob의 **방향성과 상대 스케일**을 보는 데 유용하지만, 편광·CG 계수·모드 중첩에서 절대 `Ωc`를 예측하는 모델은 아니다.

### 실험 reference로서의 경계

- Rb/Cs D1/D2 선택은 파장·자연선폭·dipole·질량·밀도를 바꾸지만 (`absorption.py:422-447,525-530`), `Rb (natural)`은 동위원소별 susceptibility 합이 아니라 우세 isotope의 선 상수와 abundance-weighted scalar density를 하나의 Λ에 넣는다.
- `Fg,mF`, polarization selection, 비대칭 branching, hyperfine optical pumping, repumping이 없다. 따라서 특정 실험의 절대 contrast/line center를 맞추는 모델이 아니다.
- CPT는 별도 bichromatic clock-CPT 엔진이 아니라 같은 weak-probe Λ의 좁은 scan preset이다 (`absorption.py:482-493,613-635`). 균형 광장, light shift, clock transition, transient pumping의 정량 reference로 쓰면 안 된다.
- `buffer_ground_relax_khz`는 ground coherence 한 항으로 transit, spin destruction, diffusion, VCC를 구분하지 않는다 (`atoms.py:171-187`; `absorption.py:553-557`). Pressure shift와 optical collisional width도 Lambda 본체에 없다.
- `group_index()`는 샘플된 `Reχ`의 수치 gradient다 (`gabes/observables.py:607-615`). 유한 pulse bandwidth, group delay, pulse distortion, 검출 SNR은 평가하지 않는다.

**실험가 관점의 판정:** resolved AT splitting, Rabi/power/beam-size trend, ground-dephasing 및 opacity 변화의 1차 sanity check에는 유용한 실제 원자물리 코드다. 그러나 현재 EIT/CPT linewidth, nonzero beam-angle curve, clock-CPT, buffer/coated-cell, 특정 hyperfine/polarization의 절대 수치에는 정량 실험 reference가 아니다.

## 오늘의 재계산

Numba warm-up 뒤 동일 프로세스에서 headless 경로를 사용했다. 시간은 이 PC의 참고값이며 병렬 스케줄링에 민감하다.

| 기본 regime | warm median compute | 중심 투과 | 현재 UI 주 결과 |
|---|---:|---:|---|
| EIT, Doppler on | 약 `0.10–0.12 s` | `0.013558281` | FWHM `0.459680 MHz`, `n_g=100305` |
| AT, Doppler off | `4.14 ms` | `0.989089759` | split `46.0 MHz` = anchor `Ωc/2π=45.968 MHz` |
| CPT, Doppler off | `3.46 ms` | `0.706327628` | FWHM `0.923317 MHz`, `n_g=95255` |

현재 기본값은 EIT `(Ωc/2π, γgg/2π, T, L)=(17.238 MHz, 57.46 kHz, 50 °C, 15 mm)`, AT `(45.968 MHz, 57.46 kHz, 25 °C, 3 mm)`, CPT `(5.746 MHz, 28.73 kHz, 25 °C, 3 mm)`다 (`absorption.py:486-493,574-585`).

### EIT/CPT 선폭과 군지수는 기본 격자에서 미수렴

| 계산 | EIT 보간 FWHM | EIT `n_g` | CPT 보간 FWHM | CPT `n_g` |
|---|---:|---:|---:|---:|
| 현 uniform 601점 | `0.229899 MHz` | `100305` | `0.805733 MHz` | `95255` |
| `sinh(5u)` 중심 집중 601점 | `0.159141 MHz` | `112708` | `0.803540 MHz` | `95272` |
| uniform 9,601점 reference | `0.159102 MHz` | `112693` | `0.803539 MHz` | `95272` |

현재 `window_fwhm()`은 half-height 바깥 **샘플 좌표**를 빼며 보간·해상도·edge 판정을 하지 않는다 (`gabes/lineshape.py:52-70`). 그래서 UI의 EIT `0.459680 MHz`는 수렴 보간값의 **2.89배**, CPT `0.923317 MHz`는 **14.9%** 크다. EIT 중심 step은 `0.229840 MHz`라 실제 수렴 FWHM 안에 한 샘플도 온전히 들어가지 않는다. 같은 601점을 중심 집중했을 때 EIT wall time은 `0.1170→0.1148 s`로 측정 오차 안이었고, 고해상도 `1.52 s`와 달리 solve 수가 동일하다. 즉 기존 P1의 격자·보간·status 교정은 물리를 보존하면서 사실상 무료다.

### 각도 모델은 여전히 정량적으로 부적합

`_affine_scan_coeffs()`는 `A=dL/ds-S_v`로 scan shift와 velocity shift를 묶는다 (`absorption.py:45-57`). 그 결과 Doppler-off인데도 angle `0→10 mrad`가 `max|Δχ|/max|χ₀| = 0.02942038`을 만든다. 속도 0에서는 residual Doppler가 0이어야 하므로 명백한 cold-limit 오류다. Optical scan operator와 velocity-only Raman operator를 분리하면 차원·solve 수는 그대로다.

더 근본적으로 `_two_photon_k_ratio()`는 `k_probe=k_coupling`으로 `|Δk|/k`만 계산하고 (`absorption.py:597-601`; `gabes/beam.py:42-55`), 이를 optical Doppler와 같은 단일 `v_z`에 곱한다 (`atoms.py:95-109,168-189`). 실제 비공선 기하는

```text
optical Doppler = k_p v_z
Raman Doppler = (k_p - k_c cosθ) v_z - k_c sinθ v_x
```

이며 `v_x`, `v_z`는 독립이다. 공용 모듈에는 이미 이 식과 RMS budget이 있다 (`gabes/doppler.py:106-172`)고 FWM용 독립 시험도 있다 (`tests/test_fwm_angular_doppler.py:67-139`). 그러나 Lambda 생산 경로는 사용하지 않는다.

현재 1D 모델의 50 °C, 3 mm 공명 투과는 0/1/2/5/10 mrad에서 각각 `0.423098 / 0.034796 / 0.011022 / 0.004028 / 0.002679`였다. 직전 보고서의 동일 OBE 2D 직접 평균은 1 mrad `0.131721`, 10 mrad `0.004026`이었고 Lambda core가 그 뒤 바뀌지 않았으므로 반례는 그대로 적용된다. 61 scan점 직접합은 1D `3–5 ms` 대비 2D `1.2–1.3 s`로 수백 배였으므로 정확한 공동 평균을 “negligible overhead”로 분류할 수 없다. P1 acceptance에 2D joint benchmark, `θ→0`, 두 velocity 축 수렴과 runtime budget을 명시하고, 생산 기본이 아니라 opt-in reference/가속 경로부터 두는 편이 안전하다.

## 문서·테스트·예제 품질

- `README.md:13,68,143`은 기능과 시험 진입점을 정확히 요약하지만, compact scalar Λ의 한계와 angle/linewidth 신뢰 경계는 적지 않는다.
- User Guide는 3준위 Λ, EIT/AT/CPT, `Reχ`, FWHM, 군지수를 설명한다 (`docs/Userguide/GABES_User_Guide_v2.html:531-534,595-613`). 그러나 “느린 빛의 서명”이라는 문구 옆에 현재 grid 미수렴, 매우 낮은 기본 EIT 투과, pulse propagation 미검증 경고가 없다.
- Guide의 EIT 그림 제목은 `Ωc=6.00 MHz`, `γgg=10.0 kHz`로 현재 기본 `17.238 MHz`, `57.46 kHz`와 다르다 (`docs/Userguide/userguide_assets/eit.png`). AT 그림은 현재 `45.968 MHz`, `57.46 kHz`와 일치한다. CPT 그림은 없다. 두 그림 모두 입력 manifest·생성 스크립트·commit 기록이 없다.
- Guide 성능표의 Lambda `~0.6 s` (`GABES_User_Guide_v2.html:777-784`)는 오늘 EIT `~0.10–0.12 s`와 맞지 않는다.
- 별도 `examples/` 디렉터리나 Lambda 전용 script/config/CSV 예제는 없다. 실제 실행 예제 역할은 `tests/test_absorption.py:135-241`의 AT/Rabi/angle/EIT/CPT와 `tests/test_schemes_render.py:97-107`의 세 preset smoke case가 대신한다.

테스트는 AT split, `sqrt(P)/d`, warm angle broadening, cold EIT, sub-natural CPT, unresolved status, Numba/NumPy parity를 다룬다. 그러나 Doppler-off angle invariance, 실제 `k_p-k_c`, 2D joint average, clustered/high-resolution FWHM·`n_g`, samples/FWHM·edge·visibility gate는 없다. 현재 warm-angle 테스트 (`tests/test_absorption.py:157-173`)는 잘못된 1D 모델의 상대 경향만 고정한다.

## 순수 코드 최적화 후보

물리 의미나 이산 평균을 바꾸지 않는 후보만 적었다.

1. **BLAS thread controller 재사용.** `core.blas_single_thread()`는 호출 때마다 `threadpool_limits(...)`를 새로 만들어 Windows DLL 목록을 다시 훑는다 (`gabes/core.py:23-34`). 빈 context median은 현 `1.015 ms`, 미리 만든 `ThreadpoolController` 사용 시 `0.0062 ms`; EIT compute median은 `101.85→94.56 ms`(**7.2%**)였고 `chi_bar`는 bit-identical이었다. 공용 저위험 최적화다.
2. **기존 `lambda-pole-residue-doppler` 구현 시 `ShiftInvertSystem` 사용.** Lambda의 shift-null block은 singular라 단순 `AffineShiftSystem` 소거가 맞지 않지만, regular detuning에서 확장하는 `ShiftInvertSystem` (`gabes/pole_doppler.py:253-337`)은 현재 601×711 이산 평균을 `max|Δχ|/max|χ|=8.37e-14`로 재현했다. NumPy mean fallback을 포함한 비영구 prototype median은 현 compute `0.0994 s` 대비 `0.0550 s`(약 **1.8×**)였다. 이것은 현 1D 동작의 순수 가속이지 2D angle physics 교정이 아니다.
3. **Affine kernel의 희소 update.** 기본 EIT의 real-basis `A_coef`, `B_coef`는 각각 81개 중 유효 원소가 4개인데, `_affine_scan_chi_real()`은 모든 `(scan,velocity)`에서 81개 원소를 다시 쓴다 (`gabes/kernels.py:571-595`). `base+sA`를 velocity loop 밖에서 만들고 `B`의 비영 원소만 갱신하면 동일 9×9 LU와 Maxwell weight를 유지할 수 있다. Pole port보다 우선순위는 낮고 parity/roundoff 회귀가 필요하다.
4. `cell_mm`·`line_strength`의 navigate-only remap과 `headless_observables()`는 이미 좋은 최적화다 (`absorption.py:559-562,681-726`). `_medium_from_params()` scalar cache는 절감 가능량이 작아 우선순위가 낮다.

## 검증

```text
python -m pytest -q \
  tests/test_absorption.py tests/test_kernels.py \
  tests/test_headless_observables.py tests/test_schemes_render.py \
  tests/test_beam.py tests/test_fwm_angular_doppler.py

64 passed in 20.44 s
```

AGENTS.md가 요구한 전체 시험도 실행했다.

```text
python -m pytest -q

collection 중단, 1 error in 36.72 s
references/AutoOD/test_autood_gui_adjustments.py가
references/AutoOD/rb_od_calculator_optimized_final.py:5의 optional sympy를 찾지 못함
```

이는 오늘 검토한 Lambda 경로의 실패가 아니라 저장소의 optional reference dependency 문제다.

## 결론과 우선순위

Scheme 2는 **실제 정상상태 3준위 Λ OBE를 구현한 유용한 반정량 실험 계획 reference**다. 특히 resolved AT splitting과 Rabi·power·beam-size·dephasing trend는 연구실의 빠른 sanity check에 적합하다. 반면 현재 기본 EIT/CPT linewidth와 group index는 격자에 미수렴하고, beam-angle path는 cold-limit와 joint-velocity geometry를 위반한다. 이 숫자들을 절대 linewidth, slow-light delay, 정렬 tolerance로 인용하면 안 된다.

권장 순서는 다음과 같다.

1. 기존 P1에서 cold operator 분리, 중심 집중 601점, 보간·resolution/edge/visibility status를 먼저 적용한다. 비용은 사실상 동일하다.
2. 실제 `k_p-k_c`와 2D 공동 평균은 공용 `doppler.py` 식을 재사용하되 opt-in reference와 직접합 수렴 gate를 먼저 만든다. 이것은 무시 가능한 비용이 아니다.
3. 그 뒤 출처 있는 collision coefficient/pressure shift를 같은 solve 차원에 연결한다.
4. 순수 속도 개선은 cached thread controller를 먼저, P1 의미가 고정된 뒤 `ShiftInvertSystem` pole port를 적용한다.
5. Hyperfine/Zeeman, coated-cell Ramsey, full VCC는 독립 실험 데이터와 runtime budget이 생길 때만 opt-in research mode로 진행한다.
