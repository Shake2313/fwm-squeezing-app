# 2026-09-21 Scheme 2 — Λ coherence (EIT / AT / CPT) 물리 검토

## 선택·범위

- 서울 현지 날짜 `2026-09-21`; `n = (21 mod 5) + 1 = 2`.
- 기준: 현재 작업 트리, `main`의 `HEAD 2998dc1`. 검토 시작 전부터 다수 파일이 staged/unstaged로 수정돼 있음. 기존 변경은 건드리지 않음.
- 실제 dropdown 순서: `gabes/schemes/__init__.py:19-25,37-39`의 `all_schemes()`를 런타임에서도 확인.

| 번호 | key | 현재 표시·물리 | 구현 |
|---:|---|---|---|
| 1 | `sas` | OD / SAS | `gabes/schemes/sas.py` |
| **2** | **`lambda`** | **Λ coherence: EIT / AT / CPT** | `gabes/schemes/absorption.py:466-800` |
| 3 | `rydberg_eit` | Rydberg-EIT electrometry | `gabes/schemes/rydberg.py` |
| 4 | `magneto` | Hanle / EIA / NMOR | `gabes/schemes/magneto.py` |
| 5 | `fwm` | Squeezing / Biphoton FWM | `gabes/schemes/fwm.py` |

`eit`, `at`, `cpt`는 직접 API·시험용 alias. dropdown 항목 3개가 아님 (`gabes/schemes/__init__.py:27-34`). `ODScheme`도 별도 검증 primitive (`:18`).

## 기존 개선 제안 먼저 검색

`CLAUDE.md`, `README.md`, `docs/checklist.json`, Scheme 2 선행 일일 보고서 14개(2026-06-26~09-16), `docs/`, `analysis/`, `gabes/`, `tests/`의 TODO/FIXME/issue·proposal 기록, 로컬 issue/TODO 이름 파일을 확인. 별도 Lambda issue-note 파일 없음. 공개 [GitHub Issues](https://github.com/Shake2313/fwm-squeezing-app/issues)는 현재 열린 항목 0개. **개선 제안은 이미 존재**하므로 새 항목을 만들기보다 그 범위와 비용을 재평가했다.

| 기존 항목 | 상태·근거 | 현 판단 |
|---|---|---|
| `lambda-spectrum-validity` | P1 ready, `docs/checklist.json:2909-2933` | 냉각 한계 오류, `Δk`, 선폭 격자·표시 신뢰성. 오늘 확인한 **광학/Raman 속도 상관 오류**도 이 항목에서 다뤄야 함. |
| `collisional-coefficient-provenance-and-pressure-shift` | P1 ready, `docs/checklist.json:2630-2667` | 적용 산술은 작지만 종·선·기체별 계수, 부호, 단위, 적용 범위와 검증값 필요. |
| `lambda-pole-residue-doppler` | P3 ready, P1 의존, `docs/checklist.json:3337-3354` | 일치하는 1차원 이산 Maxwell measure에는 순수 가속 후보. 2차원 각도 평균을 저절로 해결하지 않음. |
| `paraffin-coated-cell-extensions` | P2 parked, `docs/checklist.json:2728-2751` | Λ에는 암흑 구간의 **ground coherence**·Ramsey 위상·귀환 통계 필요. OD/SAS population reservoir 재사용 불가. |
| `full-velocity-changing-collision-kernel` | P2 parked, `docs/checklist.json:3035-3054` | 속도 class 결합, 독립 평균 상실. 실험 자료·kernel·runtime budget 없이는 기본 경로 부적합. |
| `lambda-hyperfine-resolved-manifold` | P2 parked, `docs/checklist.json:3057-3075` | 동위원소·전이·편광·독립 검증 spectrum 정해진 뒤 opt-in. Liouville 크기와 매개변수 수 증가. |

## 실제 구현과 실험 적합성

`atoms.lambda3()`는 `g1,g2,e` 3준위, `e→g1/g2` 각각 `Γ/2`, ground coherence 감쇠 `gamma_gg`를 둔다 (`gabes/atoms.py:168-190`). 고정 weak probe `Ωp=10⁻³Γ`, coupling `Ωc`, two-photon scan으로 정상상태 밀도행렬을 푼다 (`gabes/schemes/absorption.py:29-31,637-667`). 이 OBE는 EIT 암흑선, 분산, 강결합 AT 이중선을 실제 원자 응답으로 만든다. 실험용 power/직경 trend는 측정·설정한 Rabi anchor에 `√P/d`를 곱한다 (`gabes/beam.py:9-20`; `absorption.py:587-595`). 절대 `Ωc`를 편광·CG·모드 중첩에서 예측하는 모델은 아니다.

선택한 Rb/Cs D1/D2로 밀도, 질량, 파장, 자연선폭, reduced dipole을 정한다 (`gabes/schemes/absorption.py:422-447,525-530`). 다만 `Rb (natural)`은 동위원소별 susceptibility 합이 아니라 우세 isotope의 선 상수와 abundance 가중 밀도·질량을 **하나의 scalar Λ**에 쓴다 (`:434-446`). `Fg,mF`, 편광, 비대칭 branching, repumping 없음. `CPT` 역시 독립 clock-CPT 계산이 아니라 같은 weak-probe Λ의 좁은 scan preset (`:486-493,613-635`). 균형 두 광장, 시간 의존 pumping, light/magnetic shift 기준으로 쓰기 어렵다.

Doppler-on은 2 m/s, ±4σ **1차원** Maxwell class마다 601개 scan을 affine Liouvillian으로 푼다 (`gabes/schemes/absorption.py:118-140,613-635`; `gabes/kernels.py:559-640`). 검출량은 `χphys→α=k Imχ→T=e⁻ᵅᴸ`; `n_g`는 샘플된 `Reχ`의 중심 차분 (`gabes/observables.py:556-615`). 따라서 resolved AT 간격, Rabi·온도·셀 길이 경향을 보는 반정량 실험 계획에는 유용하다. 절대 선폭·contrast, buffer/coated-cell, 특정 hyperfine·편광, clock-CPT, 유한 pulse delay·SNR의 실험 기준으로는 부족하다. [EIT 정리 논문](https://journals.aps.org/rmp/abstract/10.1103/RevModPhys.77.633)과 [EIT/AT 판별 연구](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.107.163604)는 해석 범위의 물리 배경이다. 이름만으로 두 현상을 실험적으로 판별했다고 주장하지 않는다.

## 오늘 재계산: 선폭·냉각 한계

동일 프로세스 warm JIT. 시간은 이 PC의 참고값이며 부하에 민감하다. `headless_observables()` 사용; 고해상도 비교에서는 `_scan`만 임시 치환하고 생산 코드는 수정하지 않았다.

| preset | 기본 공명 투과율 | 현재 표시 선폭 | 보간된 601점 선폭 | 고해상도 9,601점 보간 선폭 | 현재 `n_g` |
|---|---:|---:|---:|---:|---:|
| EIT, Doppler on | 0.0135583 | 0.459680 MHz | 0.229899 MHz | **0.159102 MHz** | 100,305 |
| AT, Doppler off | 0.989090 | 분리 46.0 MHz | — | — | — |
| CPT, Doppler off | 0.706328 | 0.923317 MHz | 0.805733 MHz | **0.803539 MHz** | 95,255 |

AT 분리는 같은 preset `Ωc/2π=46.0 MHz`와 일치. EIT 기본 중심 step `0.229840 MHz`: 보간 선폭 기준 겨우 **1.00 sample/FWHM**. 현재 `window_fwhm()`은 half-height 바깥 샘플 좌표를 빼고 해상도·edge 판정을 안 한다 (`gabes/lineshape.py:52-70`); scheme은 양의 유한 폭이면 hero로 낸다 (`gabes/schemes/absorption.py:776-799`). 표시 EIT 선폭은 수렴값의 **2.89배**. CPT 표시값은 수렴값보다 **14.9%** 큼. EIT `n_g`는 수렴값 약 112,693보다 **11.0% 낮음**. `sinh(5u)`로 같은 601점을 중심에 모으면 EIT 보간 선폭 0.159141 MHz, `n_g=112,708` (고해상도와 각각 약 0.025%, 0.013% 차이); CPT 0.803540 MHz, `n_g=95,272`. 601점 수 동일, paired median `33.473→33.614 ms` (+0.42%). 보간·resolution/visibility status는 선폭 표시의 정확도 향상이지 새 물리 추가가 아니다.

Doppler-off에서 각도 `0→10 mrad`: `max|Δχ̄|/max|χ̄₀|=0.02942038`. 속도 0에서 `Δk·v=0`이어야 한다. `_affine_scan_coeffs()`의 `A=dL/ds−S_v`가 ground residual Doppler까지 scan 항에 섞는다 (`gabes/schemes/absorption.py:45-57,86-95`; `gabes/atoms.py:176-189`). 이 부분은 optical scan operator와 velocity-only operator를 분리하면 행렬 크기·solve 수 그대로 고칠 수 있다. `tests/test_absorption.py:157-173`은 warm angle broadening만 확인하며 냉각 불변성은 검증하지 않는다.

## 오늘 새 확인: 각도에서 **두 속도의 상관관계**

현재 `_two_photon_k_ratio()`는 `k_probe=k_coupling`을 전달하고 (`gabes/schemes/absorption.py:597-601`), `beam.angled_residual_k_ratio()`는 `|Δk|/k`만 반환한다 (`gabes/beam.py:42-55`). `AtomModel.S_v`는 이를 광학 Doppler와 같은 `v_z`에 곱한다 (`gabes/atoms.py:95-109`; `absorption.py:125-140`). 따라서 광학 detuning과 Raman detuning의 속도 상관을 강제로 1로 만든다. 실제 비평행 기하에서는

```text
광학 Doppler = k_p v_z
Raman Doppler = (k_p − k_c cosθ) v_z − k_c sinθ v_x
```

`v_x`와 `v_z`는 독립 Maxwell 성분. 같은 `k_p=k_c` 근사에서도 실제 상관계수는 작은 각도에서 약 `θ/2` (`10 mrad`에서는 0.005), 현재 1차원 계산은 1. **`|Δk|`의 크기만 맞춰도 joint spectrum은 보존되지 않는다.** 별도의 longitudinal `k_p−k_c` 누락 문제도 있다 (`absorption.py:597-601`). 실제 `k_c`를 정하려면 현재 scalar 매질에 없는 특정 두 ground hyperfine 준위·coupling 전이 선택이 필요하다. 두 문제는 구별해야 한다.

동일 3준위 Liouvillian을 종·횡속도에 각각 affine 이동시킨 비영구 2차원 직접 합으로 공명점만 대조했다. `S_e=comm_super(diag(0,0,1))`, `S_g=comm_super(diag(0,1,0))`로 두고 `L=L₀+s(dL₀/ds−S_e)+k_p v_z S_e+[(k_p−k_c cosθ)v_z−k_c sinθ v_x]S_g`를 풀었다. 50 °C, 3 mm, 기본 EIT, `k_p=k_c`, `v_z` 2 m/s·±4σ, `v_x` 4 m/s·±4σ, 동일 `χ→T` 정규화. 이 계산은 별도 구현·실험 검증이 아닌 **현 scalar 모델 안의 기하 반례**다.

| 각도 | 현재 1D `T(0)` | 2D `T(0)` | 2D / 현재 |
|---:|---:|---:|---:|
| 0 mrad | 0.423098 | 0.423098 | 1.000 |
| 1 mrad | 0.034796 | **0.131721** | **3.79** |
| 2 mrad | 0.011022 | **0.037093** | **3.37** |
| 5 mrad | 0.004028 | 0.008236 | 2.04 |
| 10 mrad | 0.002679 | 0.004026 | 1.50 |

횡격자 8→4→2 m/s에서 1 mrad `T=0.13171985→0.13172121→0.13172334`; 10 mrad `0.00402625→0.00402616→0.00402618`. 0 mrad의 2D 경로는 현재값과 기계 정밀도로 일치했다. `s=0` 비교이므로 위 scan-operator 냉각 오류는 이 차이를 설명하지 못한다. 실제 Cs 증기에서도 수 mrad 비평행 EIT 선폭 변화가 측정됐다 ([Carvalho *et al.*, PRA 70, 063818](https://journals.aps.org/pra/abstract/10.1103/PhysRevA.70.063818)); 이는 각도 민감성의 독립 근거이지 위 GABES 숫자의 실험 검증은 아니다. 충돌성 Dicke 좁힘은 별도 물리이며 ([Shuker *et al.*, PRA 76, 023813](https://journals.aps.org/pra/abstract/10.1103/PhysRevA.76.023813)), 현재 `gamma_gg`만으로 대체할 수 없다.

**비용 판정 수정:** cold operator 분리, `k_p−k_c` 부호·주파수, 601점 중심 격자·표시 gate는 추가 solve 없이 가능. 하지만 정확한 비평행 **공동** 속도 평균은 현재 `scan × v_z` 경로에서 `scan × v_z × v_x`로 커진다. 같은 61개 scan점·711개 `v_z`·357개 `v_x`의 warm 직접합 두 차례: 1D `3.09–5.26 ms`, 2D `1.19–1.30 s` (해당 구현에서 226–419배). 전체 601점의 성능 보장은 아니지만 **무시 가능한 overhead는 아님**. 단순 `|Δk|` 치환을 “거의 무료인 물리 교정”으로 인정하면 안 된다. `lambda-spectrum-validity`의 각도 합격 조건에 2D joint benchmark, `θ→0`·`v_x`/`v_z` 수렴, runtime budget을 명시하는 편이 타당하다. 그전에는 beam-angle 수치 곡선을 **정량 실험 기준으로 쓰지 말아야 한다**. [기존 체크리스트](../checklist.json)의 full-vector 문구가 이 상관 조건까지 보장하는지 재검토 필요.

## 개선 비용·순수 가속

| 제안 | 물리/동작 보존 | 비용 판단 |
|---|---|---|
| cold scan/velocity operator 분리 | compact OBE 유지, 비물리 각도 의존 교정 | 추가 solve 0; 회귀로 `doppler=off` 불변성 요구 |
| 중심 집중 601점 + half-height 보간 + 해상도·edge·visibility 상태 | OBE 동일, 정확한 샘플링·보고 | 오늘 EIT paired +0.42%; 선형 후처리 |
| 실제 광주파수 차이 `k_p−k_c` 반영 | co-propagating Raman 잔차 교정 | scalar 산술은 작음. 특정 hyperfine 전이 선택·출처 필요; 비평행 joint 평균 문제까지 해결하지는 않음 |
| 비평행 2D Maxwell 평균 | 기하의 공동 분포 보존 | 순진한 class 곱은 비싸다. opt-in reference·pole/적분 가속을 2D 직접합으로 검증해야 함 |
| 출처 있는 pressure shift/광학 dephasing | 실험 물리 확장 | 적용 `O(1)`; 계수 provenance·독립 데이터가 큰 비용 |
| full VCC, coated-cell Ramsey, hyperfine/Zeeman | 새로운 물리 | 계산·데이터 비용 큼; parked 적절 |

순수 코드 최적화: 현재 `line_strength`와 `cell_mm`는 navigate-only (`gabes/schemes/absorption.py:559-562,668-690`), headless remap median `0.126 ms` vs warm EIT compute `33.5 ms`; 유지할 구조다. `kernels.affine_scan_chi()`는 scan×velocity 행렬의 affine 결합을 매 inner-loop에서 수행 (`gabes/kernels.py:571-595`); 희소한 비영 원소만 갱신하는 specialization은 **동일 행렬·동일 이산 가중치**로 비교하면 행동 보존 가속 후보. 그러나 9×9 solve 자체가 이미 빠르므로 프로파일과 회귀 후 채택할 가치가 있다. 기존 P3 pole-residue port는 1D 평균에서 같은 2 m/s·±4σ measure와 direct-row guard를 유지하면 순수 수치 가속 가능; continuous Gaussian으로 바꾸면 계산 measure도 바뀌므로 별도 정확도 검증 필요 (`docs/checklist.json:3337-3354`). Immutable `(species,line)` 상수 cache는 가능하지만 `_medium_from_params()` (`absorption.py:422-447`)가 작은 비용이라 우선순위 낮다. `headless_observables()` 지원 (`absorption.py:470`; `gabes/schemes/base.py`)은 이미 그림 없는 batch에 적합하다.

## 문서·시험·예제와 검증

`README.md:13,68,143`은 기능·시험 진입점만 요약. User Guide는 EIT/AT 그림·`Reχ`와 Λ의 scalar 한계를 소개 (`docs/Userguide/GABES_User_Guide_v2.html:595-613,849-851`). 다만 guide의 Λ 시간 `~0.6 s` (`:781`)는 오늘 warm 계산 `~0.0335 s`와 맞지 않는다. “느린 빛의 서명” (`:604-606`) 옆에 `T=1.36%`, 미수렴 선폭, pulse bandwidth·visibility 미평가를 적어야 한다. `eit.png`, `at.png`는 재현 가능한 입력 manifest/생성 스크립트가 없고 CPT 그림도 없다. 저장소에 `examples/` 디렉터리나 Lambda 전용 script/config/CSV 예제 없음. 현재 실행 가능한 예제 역할은 `tests/test_absorption.py:135-241`의 AT·power·angle·EIT·CPT와 `tests/test_schemes_render.py:97-107`의 preset 렌더가 담당한다. 실험 참고 품질에는 입력 JSON·commit·격자/visibility status를 저장하는 headless angle/power sweep 예제가 유익하다.

핵심 시험 공백: Doppler-off angle invariance, full optical `Δk`, **2D 비평행 공동 평균**, clustered vs high-resolution FWHM/`n_g`, sample/FWHM·edge·visibility gate. 현재 warm-angle 시험 (`tests/test_absorption.py:157-173`)은 상대 경향만 검사한다. affine kernel parity는 `tests/test_kernels.py:145-191`에 있다.

```text
python -m pytest -q tests/test_absorption.py tests/test_kernels.py tests/test_headless_observables.py tests/test_schemes_render.py
55 passed (12.35 s)

python -m pytest -q
collection 중단 (30.51 s): references/AutoOD/test_autood_gui_adjustments.py →
rb_od_calculator_optimized_final.py:5의 optional sympy import 불가
```

결론: Scheme 2는 실제 정상상태 Λ OBE에 기반한 **반정량 실험 계획 reference**다. resolved AT·Rabi trend는 쓸 만하다. 현재 EIT 선폭과 유한 각도 곡선은 실험 수치로 인용하기 어렵다. 우선 기존 P1에서 냉각 불변성·선폭 표시를 거의 같은 비용으로 바로잡고, 각도는 **2D 공동 속도 평균의 정확도·계산비를 별도 판정**해야 한다. 그 뒤 충돌 계수와 1D pole 가속, 마지막으로 자료가 마련된 확장 모델을 검토하는 순서가 적절하다.
