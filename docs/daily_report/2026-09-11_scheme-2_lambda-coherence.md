# 2026-09-11 Scheme 2 — Lambda coherence (EIT / AT / CPT) 물리 검토

## 1. 선택 규칙과 현재 다섯 scheme

- 로컬 날짜: `2026-09-11` (Asia/Seoul)
- 계산: `n = (11 mod 5) + 1 = 2`
- 오늘 대상: **Scheme 2, `LambdaScheme` — Λ coherence (EIT / AT / CPT)**
- 검토 기준: `main`, `HEAD c0c46f0`, 현재의 staged/unstaged 작업 트리 포함

런타임 registry를 직접 조회한 순서는 다음과 같다
(`gabes/schemes/__init__.py:12-25,37-39`).

| 번호 | registry key | 표시 이름 | 구현 |
|---:|---|---|---|
| 1 | `sas` | Absorption spectroscopy (OD / SAS) | `gabes/schemes/sas.py` |
| 2 | `lambda` | Λ coherence (EIT / AT / CPT) | `gabes/schemes/absorption.py:466-800` |
| 3 | `rydberg_eit` | Rydberg-EIT electrometry | `gabes/schemes/rydberg.py` |
| 4 | `magneto` | Magneto-optics (Hanle/MOR) | `gabes/schemes/magneto.py` |
| 5 | `fwm` | Four-wave mixing (Squeezing / Biphoton) | `gabes/schemes/fwm.py` |

`eit`, `at`, `cpt`는 테스트와 direct API를 위한 alias이고 dropdown의 별도 scheme은
아니다 (`gabes/schemes/__init__.py:27-34`). `ODScheme`도 현재 dropdown에 없는 2-level
검증 primitive다.

작업 트리에는 검토 전부터 많은 사용자 변경이 있었다. 특히 오늘 검토한
`gabes/schemes/absorption.py`, `tests/test_absorption.py`, README/User Guide/checklist도
이미 수정 상태다. 이를 고치거나 되돌리지 않았고, 오늘은 이 보고서만 추가했다.

## 2. 기존 개선 제안부터 확인한 결과

검색 범위는 `CLAUDE.md`, `README.md`, `docs/checklist.json`, Scheme 2 선행 일일 보고서
12개(2026-06-26~2026-09-06), 코드의 TODO/FIXME/issue/proposal 문자열, 그리고
TODO/issue/proposal 이름의 로컬 파일이다. 별도의 로컬 issue note 파일은 없었지만,
**Scheme 2 개선안 자체는 이미 충분히 존재한다.** 따라서 새 기능을 임의로 늘리기보다
기존 제안의 현재 유효성과 비용을 재검토했다.

| 기존 제안 | 현재 상태 | 오늘의 판단 |
|---|---|---|
| `lambda-spectrum-validity` | P1, `ready` | cold-limit 오류, EIT/CPT 격자와 hero 신뢰성 문제를 함께 다루는 최우선 항목 (`docs/checklist.json:1125-1147`) |
| `collisional-coefficient-provenance-and-pressure-shift` | P1, `ready` | 계산 추가량은 작지만 계수의 출처·부호·단위·불확도 확보가 핵심 (`docs/checklist.json:854-889`) |
| `paraffin-coated-cell-extensions` | P2, `parked` | Lambda에는 population reservoir 복사가 아니라 ground coherence와 bright/dark Ramsey 동역학이 필요 (`docs/checklist.json:950-970`) |
| `full-velocity-changing-collision-kernel` | P2, `parked` | 검증 데이터와 runtime budget이 있는 opt-in 연구 solver만 타당 (`docs/checklist.json:1249-1268`) |
| `lambda-hyperfine-resolved-manifold` | P2, `parked` | isotope/transition/polarization/held-out spectrum이 정해지기 전에는 일반화하지 않는 판단이 맞음 (`docs/checklist.json:1271-1288`) |

직전 Scheme 2 검토의 기준 `d98bdf4`에서 현재 `c0c46f0`까지 committed 변경은
checklist, FWM, Streamlit 영역뿐이다. 현재 Lambda 계산식, 601점 scan, residual-k 전달,
문서 자산에는 9월 6일 이후 물리 변경이 없다. 당시 완료된 `line_strength` navigate-only
변경은 현재 작업 트리에도 유지돼 있다 (`gabes/schemes/absorption.py:561-562,681-688`).

## 3. 현재 구현이 푸는 실제 물리

### 3.1 정상상태 3-level OBE

`atoms.lambda3()`는 `g1`, `g2`, `e`의 세 준위를 구성하고 excited state에서 두 ground
state로 각각 `Γ/2` 붕괴를 둔다. `gamma_gg`는 population exchange가 아니라
`rho_g1,g2`에 직접 작용하는 phenomenological ground-coherence `T2` 감쇠율이다
(`gabes/atoms.py:168-190`).

Hamiltonian에는 고정된 weak probe `Ω_p = 10^-3 Γ`, coupling field, coupling detuning,
two-photon scan이 들어간다 (`gabes/schemes/absorption.py:29-31,637-667`). 따라서 EIT
투명창, 정상분산, strong-coupling dressed-state AT doublet, ground-coherence-limited dark
resonance는 임의 curve가 아니라 steady-state density-matrix 해다. 이 핵심은 실험 전
sanity check에 쓸 수 있는 실제 원자물리다.

다만 CPT는 balanced bichromatic clock-CPT 모델이 아니라 같은 weak-probe Lambda의
좁은 scan preset이다 (`gabes/schemes/absorption.py:486-493,613-635`). 실제 hyperfine
branching, optical pumping, Zeeman/polarization, clock light shift, pumping transient는 없다.

### 3.2 실험 노브와 매질

Rb/Cs 및 D1/D2 선택은 밀도, 질량, 파장, 자연선폭, lifetime-derived reduced dipole을
바꾼다 (`gabes/schemes/absorption.py:422-447`). Coupling power와 1/e² beam diameter는
사용자가 정한 Rabi anchor를 `sqrt(P)/d`로 스케일한다
(`gabes/beam.py:9-20`; `gabes/schemes/absorption.py:587-595`). 따라서 power/diameter
trend에는 유용하지만 polarization, 실제 transition dipole, mode overlap로부터 절대
`Ω_c`를 예측하는 calibration은 아니다.

`Rb (natural)`은 85Rb/87Rb susceptibility를 각각 풀어 더하지 않고, 우세 isotope의
line constants와 abundance-weighted density/mass를 하나의 scalar Lambda에 넣는다
(`gabes/schemes/absorption.py:434-446`). Pure isotope도 특정 `Fg→Fe`, `mF`, polarization을
고르지 않는다. 매질 이름은 hyperfine-resolved 실험선이 아니라 **effective scalar
D-line medium**으로 읽어야 한다.

### 3.3 Doppler 평균, Beer–Lambert, group index

Doppler-on Lambda는 `dv=2 m/s`, `±4σ` Maxwell grid에서 scan×velocity steady-state를
계산한다 (`gabes/schemes/absorption.py:118-140`). 기본 50 °C EIT는 601 scan point와
약 711 velocity class를 사용한다. 결과는

```text
chi_phys = -2 N |d|^2 line_strength chi_bar / (epsilon_0 hbar)
alpha = k Im(chi_phys)
T = exp(-alpha L)
```

로 변환된다 (`gabes/observables.py:517-542`). `n_g`는 dilute-medium
`n ≈ 1 + Re(chi)/2`와 sampled local gradient로 계산한다
(`gabes/observables.py:550-576`). 이는 CW local-dispersion 진단이지 finite pulse의
delay, distortion, bandwidth, 검출 SNR까지 포함한 slow-light propagation 결과는 아니다.

## 4. 실험 원자물리 reference 적합성

| 질문 | 적합성 | 이유 |
|---|---|---|
| `Ω_c` 증가에 따른 EIT→AT와 split sanity check | 좋음 | 3-level OBE와 resolved absorption maxima를 직접 사용 |
| coupling power/diameter 및 Raman `T2` trend | 유용 | anchored `sqrt(P)/d`와 단일 coherence-dephasing rate |
| warm-cell beam misalignment sensitivity | **수정 후 유용** | residual-Doppler 방향은 맞지만 현재 cold invariant를 위반 |
| 특정 isotope/hyperfine/polarization의 contrast·width | 낮음 | scalar 3-level, symmetric branching, Zeeman/pumping 부재 |
| clock-CPT light shift/자기장 민감도 | 부적합 | balanced two-field clock model과 time dynamics 부재 |
| finite-pulse slow-light delay/왜곡/효율 | 부적합 | local `n_g`만 있고 visibility/bandwidth gate가 없음 |
| buffer/coated-cell EIT absolute reference | 부적합 | collision, diffusion, dark flight, Ramsey phase가 단일 `T2`에 축약됨 |

결론적으로 Scheme 2는 **정상상태 EIT/AT/CPT line-shape와 knob sensitivity를 보는
반정량적 experiment-planning reference**다. 특히 resolved AT는 유용하다. 특정 셀의
절대 선폭·contrast·clock shift 또는 usable slow-light delay를 인용하는 기준으로는 아직
부족하다.

## 5. 현재 작업 트리 수치 재검증

Numba warm-up 뒤 같은 프로세스에서 반복 측정했다. 절대 wall time은 Windows 부하에
따라 변하므로 물리값과 상대비를 중심으로 본다.

| regime | 기본 계산 median | 현재 출력 |
|---|---:|---|
| EIT | `36.16 ms` | `T_res=0.0135583`, 표시 FWHM `0.459680 MHz`, `n_g=100,305` |
| AT | `3.246 ms` | split `46.0 MHz`, expected `46.0 MHz`, `T_center=0.989090` |
| CPT | `3.218 ms` | `T_res=0.706328`, 표시 FWHM `0.923317 MHz`, `n_g=95,255` |

### 5.1 EIT/CPT 601점 격자는 hero 선폭에 부족하다

모든 regime이 601점 uniform scan을 사용한다 (`gabes/schemes/absorption.py:613-635`).
`window_fwhm()`은 half-height를 둘러싼 바깥 sample의 간격을 그대로 반환하며 보간,
samples/FWHM, edge clearance를 보고하지 않는다 (`gabes/lineshape.py:52-70`).

| 경우 | 중앙 step | 표시 방식 width | 선형보간 FWHM | samples/FWHM | 중앙 `n_g` |
|---|---:|---:|---:|---:|---:|
| EIT uniform 601 | `0.229840 MHz` | `0.459680 MHz` | `0.229899 MHz` | `1.0003` | `100,305` |
| EIT `sinh(5u)` clustered 601 | `0.015488 MHz` | `0.186156 MHz` | `0.159141 MHz` | `10.28` | `112,708` |
| EIT uniform 9,601 | `0.014365 MHz` | `0.172380 MHz` | `0.159102 MHz` | `11.08` | `112,693` |
| CPT uniform 601 | `0.230829 MHz` | `0.923317 MHz` | `0.805733 MHz` | `3.49` | `95,255` |
| CPT clustered 601 | `0.015555 MHz` | `0.834352 MHz` | `0.803540 MHz` | `51.66` | `95,272` |
| CPT uniform 9,601 | `0.014427 MHz` | `0.807903 MHz` | `0.803539 MHz` | `55.70` | `95,272` |

기본 EIT hero FWHM은 converged/interpolated `0.159102 MHz`보다 약 `2.89x` 크고,
`n_g`는 약 `11.0%` 작다. CPT 표시 width도 converged 값보다 약 `14.9%` 크다. 현재
`width_valid`는 finite positive만 검사하므로 이 값들을 hero로 계속 표시한다
(`gabes/schemes/absorption.py:776-799`).

같은 601점을 중앙에 모은 paired benchmark의 median 비는 clustered/uniform
`1.0020`, 즉 약 `+0.2%`였다. solve 수를 늘리지 않고 high-resolution 기준에 EIT FWHM
`0.025%`, `n_g 0.013%` 이내로 맞췄다. 반면 9,601점 uniform EIT는 이 실행에서
`552 ms`로 약 15배 비쌌다. 기존 `<=20%` budget 안에서 clustered 601과 readout gate가
가장 합리적이다.

### 5.2 residual Doppler cold-limit 오류도 재현된다

현재 affine 조립은 scan coefficient에서 `S_v`를 빼고 velocity coefficient로 다시
사용한다 (`gabes/schemes/absorption.py:45-57`). Angle-dependent ground shift가 optical
scan operator와 섞이므로 원자 속도가 0인 Doppler-off 계산도 angle에 따라 변한다.

- Doppler off, 0→10 mrad: `max|Delta chi_bar| / max|chi_bar| = 2.942038%`
- 10 mrad에서 `|Delta k|/k = 0.009999958`
- Doppler on, 50 °C, 3 mm: 공명 transmission `0.423098 → 0.00267946`

실제 residual `Delta k·v`는 `v=0`에서 0이어야 한다. scan operator와 velocity-only
operator를 분리하면 Hilbert dimension, scan 수, velocity class, solve 수가 모두 그대로라
물리를 바로잡으면서 계산 오버헤드는 사실상 없다.

또 `_two_photon_k_ratio()`는 두 field에 같은 `medium["k_vec"]`를 넘겨 0 mrad에서
residual을 정확히 0으로 만든다 (`gabes/schemes/absorption.py:597-601`;
`gabes/beam.py:36-55`). 실제 co-propagating fields는 ground-state splitting 때문에
작은 longitudinal `Delta k`가 있다. Pure-isotope 경로에서 실제 두 wavevector를 쓰는 것은
scalar coefficient 변경이며 추가 solve가 필요 없다.

기본 EIT의 `n_g=100,305`를 15 mm에 단순 적용하면 약 `5.02 µs`지만 같은 점의
transmission은 `1.36%`이고 linewidth도 미수렴이다. 따라서 이 수치를 usable pulse
delay로 승격하면 안 된다.

## 6. 기존 물리 개선안의 비용과 물리 보존성

| 개선 | 물리 보존성 | 예상 계산비용 | 권고 |
|---|---|---:|---|
| scan/velocity operator 분리 | 잘못 섞인 Doppler 항만 분리, 3-level model 유지 | 추가 solve 없음 | 즉시 P1 |
| clustered 601 + interpolated width/status | 같은 OBE/같은 점 수, 해석 신뢰도 개선 | 측정 `+0.2%`, postprocess `O(N)` | 즉시 P1 |
| samples/edge/contrast/visibility gate와 `T`, delay 병기 | 모델을 바꾸지 않고 과신 방지 | `O(N)`, 추가 solve 없음 | P1과 함께 |
| pure-isotope collinear `Delta k` | 실제 두 optical frequency의 momentum bookkeeping | scalar 산술 | operator 수정 뒤 P2 |
| sourced pressure shift/elastic width | line center와 homogeneous width 개선 | `O(1)` lookup, solve 차원 동일 | provenance 확보 후 P1 |
| natural-Rb isotope별 incoherent sum | scalar mixture보다 실제 isotope 구조에 가까움 | isotope별 solve로 대략 2배 | opt-in 또는 pure isotope 우선 |
| coated-cell ground-coherence/bright-dark model | Ramsey/dark-flight 물리 보존 | trajectory quadrature로 중간~높음 | geometry/data 전까지 parked |
| full VCC/diffusion | velocity redistribution/Dicke physics | velocity class 결합으로 높음 | 한 scheme·한 dataset 연구 mode |
| hyperfine/Zeeman Lambda | line assignment, polarization, pumping | Liouville 차원 급증 | 특정 transition opt-in |

Pressure shift/elastic width의 병목은 runtime이 아니라 신뢰할 계수와 불확도다. 출처 없는
빠른 phenomenological knob는 experimental reference의 정확도를 높이지 않는다. Coated
cell도 단일 exponential `T2`만 하나 더 넣으면 싸지만 bright/dark return과 Ramsey phase를
잃으므로 물리 보존 개선으로 볼 수 없다.

## 7. 동작을 바꾸지 않는 순수 코드 최적화

### 완료된 고효율 항목

`line_strength`는 heavy cache key에서 빠져 있고 저장된 `chi_bar`를 observables에서 live
scale한다 (`gabes/schemes/absorption.py:561-562,668-688`). 현재 측정에서 full EIT
recompute `35.73 ms`, 같은 raw의 headless remap `0.263 ms`, 약 `136x` 차이다. 물리 배열은
`line_strength`와 독립이고 회귀가 이 cache contract를 고정한다
(`tests/test_absorption.py:110-132`). 이 최적화는 완료됐으며 다시 구현할 필요가 없다.

### 남은 exact 후보는 낮은 우선순위

기본 EIT에서 다음 immutable construction median을 측정했다.

| 작업 | median |
|---|---:|
| `_medium_from_params()` | `27.1 µs` |
| `atoms.lambda3()` | `214.8 µs` |
| `doppler.velocity_grid()` | `27.5 µs` |
| 합계 | `269.4 µs` |

exact parameter tuple의 작은 LRU와 read-only array를 쓰면 behavior-preserving이지만,
이 합계는 35.73 ms EIT compute의 약 `0.75%`에 불과하다. cold repeated-scan workload에서
실제 병목으로 확인될 때만 적용할 가치가 있다. 현재 affine kernel, real-coordinate solve,
figureless observables가 이미 큰 비용을 줄이고 있어, 오늘 검사에서는 그 밖의 고가치
behavior-preserving 최적화를 찾지 못했다.

## 8. 문서·테스트·예제의 reference 품질

### 문서

- README의 한 줄 요약은 기능은 맞지만 weak-probe/scalar/Rabi-anchor/natural-Rb
  effective-medium 한계를 밝히지 않는다 (`README.md:13`).
- 모듈 docstring은 knob가 `Γ` 단위라고 쓰지만 UI는 이미 MHz/kHz의
  `physical-units-v3`다 (`gabes/schemes/absorption.py:16-18,468-469`).
- `buffer_ground_relax_khz` help는 buffer-gas collision처럼 읽히지만 구현은 직접 지정한
  ground-coherence dephasing rate다 (`gabes/schemes/absorption.py:553-557`;
  `gabes/atoms.py:168-190`).
- User Guide는 EIT의 가파른 분산을 “느린 빛의 서명”으로 설명하지만 낮은 transmission,
  격자 미수렴, finite-pulse 한계를 붙이지 않는다
  (`docs/Userguide/GABES_User_Guide_v2.html:595-611`). Lambda 한계 표도 full
  polarization/Zeeman 부재만 기록한다 (`...:849-851`).
- Guide는 각 scheme의 Reference/about 패널에 1차 문헌이 있다고 말하지만
  (`...:915-925`), `LambdaScheme`은 `info()`를 override하지 않는다
  (`gabes/schemes/base.py:116-117`; `gabes/schemes/absorption.py:466-800`). Lambda 전용
  source/provenance가 없다.
- coated-cell checklist의 area가 존재하지 않는 `schemes/lambda_system.py`를 가리킨다
  (`docs/checklist.json:950-953`). 실제 구현 파일은 `schemes/absorption.py`다.

### 테스트

현재 테스트는 AT split≈`Ω_c`, power/diameter scaling, warm-angle broadening, cold EIT,
sub-natural CPT, line-strength live remap, unresolved AT/EIT status를 검사한다
(`tests/test_absorption.py:110-241`). Affine kernel parity, headless, render contract도 별도
테스트가 있다 (`tests/test_kernels.py:103-136`;
`tests/test_headless_observables.py:38-89`; `tests/test_schemes_render.py:97-107`).

그러나 다음 회귀가 없다.

- Doppler-off angle invariance
- EIT/CPT default FWHM 및 `n_g`의 high-resolution convergence
- interpolated edges, samples/FWHM, scan-edge/contrast/visibility gate
- co-propagating hyperfine `Delta k`
- natural-Rb effective-medium provenance/status
- transmission·pulse bandwidth·group delay의 동시 validity gate

Warm-angle test는 의도한 warm effect만 검사하므로 cold operator 오류를 잡지 못한다.
현재 unresolved test도 crossing 자체가 없을 때만 막아 기본 EIT의 1 sample/FWHM을 hero로
남긴다.

### 예제와 그림

Lambda 전용 실행 script/config/CSV example은 없다. 가장 가까운 실행 예제는
`tests/test_absorption.py:135-241`의 direct API 호출이다. User Guide의 EIT/AT PNG를
직접 확인했으며 둘 다 2026-06-08 생성본이다.

- `docs/Userguide/userguide_assets/eit.png`: `Omega_c=6.00 MHz`,
  `gamma_gg=10.0 kHz`; 현재 EIT default와 다름
- `docs/Userguide/userguide_assets/at.png`: `Omega_c=45.97 MHz`,
  `gamma_gg=57.5 kHz`; 현재 AT scale과 비슷함
- CPT 그림은 없고, 두 PNG 모두 source script, full parameter set, commit, grid/status,
  scalar-model badge가 없다.

Coupling power/diameter, angle, dephasing, cell-length sweep을 headless CSV/Markdown로
내보내고 parameter JSON, commit, grid/status, scalar-model badge를 묶는 재현 가능한
Lambda 예제는 solver 변경 없이 reference 품질을 크게 높일 수 있다.

## 9. 테스트 결과

```powershell
python -m pytest -q tests/test_absorption.py tests/test_kernels.py `
  tests/test_headless_observables.py tests/test_schemes_render.py `
  tests/test_docs_consistency.py
# Lambda/kernel/headless/render 53 passed;
# docs consistency 1 failed: 기존 작업 트리에서 삭제된 FWM_physics.tex 부재

python -m pytest -q
# collection error after 37.50 s:
# references/AutoOD가 undeclared optional dependency sympy를 요구
```

첫 실패는 Scheme 2와 무관한 삭제 상태
`docs/FWM physics and analytic reconstruction/FWM_physics.tex` 때문이다. 루트 수집 중단도
`references/AutoOD/rb_od_calculator_optimized_final.py:5`의 `sympy.physics.wigner` import
때문이다. 관련 Lambda 테스트에서 물리/수치 회귀 실패는 없었다. 기존 사용자 파일은
복원하거나 수정하지 않았다.

## 10. 최종 판단과 권고 순서

Scheme 2의 compact 3-level OBE, Maxwell average, anchored Rabi scaling,
Beer–Lambert readout은 **EIT/AT/CPT의 방향·형태·강결합 split을 빠르게 확인하는 유용한
반정량적 도구**다. 하지만 기본 EIT linewidth는 약 2.89배 과대 표시되고, CPT width도
14.9% 높으며, angle residual이 cold-limit를 위반한다. 큰 `n_g`도 낮은 transmission과
격자 상태 없이 slow-light 성능으로 읽으면 안 된다.

권고 순서는 다음과 같다.

1. `lambda-spectrum-validity` P1을 완료한다: operator 분리, clustered 601,
   interpolated FWHM, samples/edge/contrast/visibility status, transmission·delay 병기.
2. Doppler-off angle-invariance와 high-resolution EIT/CPT 회귀를 먼저 추가한다.
3. Pure-isotope co-propagating `Delta k`, natural-Rb effective-medium 경고,
   `buffer_ground_relax` 명칭/도움을 정리한다.
4. 출처와 uncertainty가 있는 pressure shift/elastic width만 low-order model에 넣는다.
5. provenance manifest가 있는 Lambda sweep example을 만들고 EIT/AT 그림을 갱신하며
   CPT 예제를 추가한다.
6. Coated-cell, VCC, hyperfine/Zeeman은 geometry, sourced rate, calibration/holdout data,
   runtime budget이 준비될 때 opt-in research mode로 시작한다.

가장 중요한 개선은 새 물리를 크게 추가하는 일이 아니라, 이미 계산한 compact physics를
거의 같은 비용으로 올바르게 샘플링하고 과신하지 않도록 상태를 표시하는 일이다.
