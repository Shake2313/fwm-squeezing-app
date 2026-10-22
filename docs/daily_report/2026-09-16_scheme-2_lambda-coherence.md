# 2026-09-16 Scheme 2 — Lambda coherence (EIT / AT / CPT) 물리 검토

## 1. 선택과 현재 순서

- 로컬 날짜: `2026-09-16`, Asia/Seoul.
- 계산: `n = (16 mod 5) + 1 = 2`.
- 대상: **Scheme 2, `LambdaScheme` — Λ coherence (EIT / AT / CPT)**.
- 기준: `main`, `HEAD 08268ab`, 현재 staged/unstaged 작업 트리 포함.

런타임 `all_schemes()` 직접 조회 결과. 정의·순서 근거:
`gabes/schemes/__init__.py:12-26, 37-39`.

| 번호 | registry key | 표시 이름 | 구현 |
|---:|---|---|---|
| 1 | `sas` | Absorption spectroscopy (OD / SAS) | `gabes/schemes/sas.py` |
| 2 | `lambda` | Λ coherence (EIT / AT / CPT) | `gabes/schemes/absorption.py:466-800` |
| 3 | `rydberg_eit` | Rydberg-EIT electrometry | `gabes/schemes/rydberg.py` |
| 4 | `magneto` | Magneto-optics (Hanle/MOR) | `gabes/schemes/magneto.py` |
| 5 | `fwm` | Four-wave mixing (Squeezing / Biphoton) | `gabes/schemes/fwm.py` |

`eit`, `at`, `cpt`: direct API·시험용 alias. dropdown 별도 scheme 아님
(`gabes/schemes/__init__.py:28-34`). `ODScheme`: 같은 파일 안 2-level 검증
primitive. dropdown 미등록.

작업 트리 검토 전부터 대량 수정 상태. Scheme 2 핵심 파일도 staged/unstaged가
겹침: `gabes/schemes/absorption.py`, `tests/test_absorption.py`, README, User Guide,
checklist. 기존 변경 복원·수정 안 함. 오늘 새 보고서만 추가.

## 2. 기존 제안 먼저 검색

검색 범위:

- `CLAUDE.md`, `README.md`, `docs/checklist.json`.
- Scheme 2 선행 보고서 13개, `2026-06-26`~`2026-09-11`.
- `docs`, `analysis`, `gabes`, `tests` 안 TODO/FIXME/issue/proposal/improvement 문자열.
- TODO/issue/proposal 이름 파일.

별도 로컬 issue-note/TODO 파일 없음. **개선 제안은 이미 충분함.** 새 physics
proposal 추가 안 함. 기존 항목 비용·물리 보존성 재평가.

| 기존 항목 | 상태 | 현재 판단 |
|---|---|---|
| `lambda-spectrum-validity` | P1, `ready` | 최우선. cold-limit 오류, residual `Δk`, EIT/CPT 격자, hero 신뢰성 함께 수정. `docs/checklist.json:2137-2161` |
| `collisional-coefficient-provenance-and-pressure-shift` | P1, `ready` | 산술비 작음. 어려움은 계수 출처·적용범위·부호·불확도. `docs/checklist.json:1858-1895` |
| `lambda-pole-residue-doppler` | P3, `ready`, P1 의존 | 9월 14일 추가. 현재 2 m/s Maxwell 합 가속. 절대 절감 작음. `docs/checklist.json:2565-2582` |
| `paraffin-coated-cell-extensions` | P2, `parked` | Lambda는 population reservoir 복사 불가. ground coherence, bright/dark flight, Ramsey phase 필요 |
| `full-velocity-changing-collision-kernel` | P2, `parked` | dataset·kernel·runtime budget 없는 기본 탑재 위험. `docs/checklist.json:2263-2282` |
| `lambda-hyperfine-resolved-manifold` | P2, `parked` | isotope/transition/polarization/holdout spectrum 선택 전 일반화 금지. `docs/checklist.json:2285-2303` |

직전 검토 뒤 committed Lambda source/test 변경 없음. 이 파일들의 마지막 commit:
`d98bdf4` (`2026-09-05`, `perf(core): streamline solves and readouts`). 현재 실제
working-tree source도 9월 6일 이후 수정 없음. 변화는 주로 FWM와 checklist. 새 Scheme 2
항목은 `lambda-pole-residue-doppler`.

## 3. 현재 구현: 무슨 물리 푸나

### 3.1 정상상태 3-level OBE

준위: `g1`, `g2`, `e`. excited decay: `e→g1`, `e→g2` 각각 `Γ/2`. Ground
coherence에 phenomenological `gamma_gg` 적용
(`gabes/atoms.py:168-190`). 이것은 population exchange, transit, diffusion 각각을
분리한 모델 아님.

Hamiltonian:

```text
H_g2g2 = Δc - s
H_g1e = Ωp/2
H_g2e = Ωc/2
Ωp = 10^-3 Γ
```

근거: `gabes/schemes/absorption.py:29-31, 637-667`. Steady-state density matrix
직접 풀이. 따라서 EIT transparency, normal dispersion, strong-field AT doublet,
ground-coherence-limited dark resonance: 임의 fit curve 아님. 실제 원자물리 핵심 있음.

CPT 해석 경계 중요. 현재 CPT = 같은 weak-probe Λ 모델의 좁은 two-photon scan preset
(`gabes/schemes/absorption.py:486-493, 613-635`). Balanced bichromatic clock-CPT,
optical-pumping transient, Zeeman clock shift 모델 아님.

### 3.2 실험 노브와 매질

지원 선택: natural Rb, ⁸⁵Rb, ⁸⁷Rb, ¹³³Cs; D1/D2
(`gabes/species.py:181`; `gabes/schemes/absorption.py:525-530`). 선택에 따라 밀도,
질량, 파장, 자연선폭, lifetime-derived reduced dipole 변경
(`gabes/schemes/absorption.py:422-447`; `gabes/species.py:192-248`).

Coupling power와 1/e² beam diameter:

```text
Ωc = Ωanchor sqrt(P/Pref) (dref/d)
```

근거: `gabes/beam.py:9-20`; `gabes/schemes/absorption.py:587-595`. Power/waist
trend 유용. 단, 실제 polarization·hyperfine CG·mode overlap에서 절대 Ωc 예측 아님.
사용자가 정한 Rabi anchor 중심.

`Rb (natural)`: 85/87 isotope별 hyperfine susceptibility 합 아님. 우세 isotope line
constant + abundance-weighted density/mass를 하나의 scalar Λ에 사용
(`gabes/schemes/absorption.py:434-446`). Pure isotope도 특정 `Fg→Fe`, `mF`, polarization
미선택. 이름은 **effective scalar D-line medium**으로 읽어야 함.

### 3.3 Doppler 평균과 readout

Doppler-on: scan 601점 × `dv=2 m/s`, `±4σ` Maxwell class. Affine Liouvillian kernel로
batched solve (`gabes/schemes/absorption.py:45-57, 118-140`; `gabes/kernels.py`,
`tests/test_kernels.py:145-191`). Doppler-off: velocity 0 한 class.

물리 susceptibility·흡수·투과:

```text
χphys = -2 N |d|² line_strength χbar / (ε0 ħ)
α = k Im(χphys)
T = exp(-αL)
```

근거: `gabes/observables.py:539-568`; scheme 연결
`gabes/schemes/absorption.py:681-690`.

`n_g`: sampled local derivative 기반 dilute-medium group index
(`gabes/observables.py:590-607`; `gabes/schemes/absorption.py:776-793`). CW local
dispersion 진단. Finite pulse delay·왜곡·효율·검출 SNR 아님.

## 4. 실험 연구자 관점 적합성

| 실험 질문 | 적합성 | 판단 |
|---|---|---|
| Ωc 증가 시 EIT→AT, resolved split scale | 좋음 | 3-level OBE, 실제 absorption maxima, split≈Ωc 시험 있음 |
| power/diameter, `T2`, 온도 trend | 유용 | anchored `sqrt(P)/d`, Maxwell 평균, Beer–Lambert |
| warm-cell beam misalignment trend | **수정 후 유용** | residual Doppler 방향 맞음. 현재 cold-limit 위반 |
| 특정 isotope/hyperfine/polarization contrast·width | 낮음 | scalar line, symmetric branching, Zeeman·optical pumping 없음 |
| clock-CPT light shift·magnetic sensitivity | 부적합 | balanced clock fields·시간동역학 없음 |
| slow-light usable delay·bandwidth·pulse distortion | 부적합 | local `n_g`만 있음. transmission/visibility/bandwidth gate 없음 |
| buffer/coated-cell absolute linewidth | 부적합 | collisions·diffusion·dark-flight를 단일 `gamma_gg`로 축약 |

판정: **실제 EIT/AT/CPT 정상상태 물리와 knob 방향을 주는 반정량적
experiment-planning reference**. 특히 resolved AT splitting, coupling scaling sanity check에
좋음. 특정 실험의 절대 linewidth·contrast·clock shift·pulse delay 기준으로는 부족.

## 5. 현재 수치 재검증

Windows, 같은 프로세스 warm JIT 뒤 median. 절대 시간은 PC 부하 의존.

| regime | 기본 compute median | 현재 hero/readout |
|---|---:|---|
| EIT | `33.568 ms` | `T_res=0.0135583`, 표시 FWHM `0.459680 MHz`, `n_g=100,305` |
| AT | `3.316 ms` | split `46.0 MHz`, expected `46.0 MHz`, `T_center=0.989090` |
| CPT | `3.286 ms` | `T_res=0.706328`, 표시 FWHM `0.923317 MHz`, `n_g=95,255` |

### 5.1 EIT/CPT linewidth 미수렴

현재 모든 regime uniform 601점 (`gabes/schemes/absorption.py:613-635`).
`window_fwhm()`은 half-height 바깥 sample 간격 반환. edge interpolation·sample/FWHM·edge
clearance 없음 (`gabes/lineshape.py:52-70`). Scheme은 finite positive면 hero 허용
(`gabes/schemes/absorption.py:776-799`).

독립 보간 재계산:

| 경우 | 중앙 step | 보간 FWHM | samples/FWHM | `n_g` |
|---|---:|---:|---:|---:|
| EIT uniform 601 | `0.229840 MHz` | `0.229899 MHz` | `1.0003` | `100,305` |
| EIT `sinh(5u)` clustered 601 | `0.015488 MHz` | `0.159141 MHz` | `10.28` | `112,708` |
| EIT uniform 9,601 | `0.014365 MHz` | `0.159102 MHz` | `11.08` | `112,693` |
| CPT uniform 601 | `0.230829 MHz` | `0.805733 MHz` | `3.49` | `95,255` |
| CPT clustered 601 | `0.015555 MHz` | `0.803540 MHz` | `51.66` | `95,272` |
| CPT uniform 9,601 | `0.014427 MHz` | `0.803539 MHz` | `55.70` | `95,272` |

현재 **표시** EIT FWHM `0.459680 MHz`: converged/interpolated `0.159102 MHz`의
`2.89×`. 현재 `n_g`: 약 `11.0%` 낮음. CPT 표시값: converged 값보다 약 `14.9%`
높음. Clustered 601은 solve 수 그대로. 직전 paired benchmark 비용비 `1.0020`; 이번
수치도 9,601 reference와 EIT FWHM `0.025%`, group index `0.013%` 이내.

### 5.2 residual Doppler cold-limit 위반

현재 `_affine_scan_coeffs()`가 scan coefficient에서 `S_v`를 빼고 velocity coefficient로
다시 둠 (`gabes/schemes/absorption.py:45-57`). Angle-dependent ground shift가 optical
scan operator와 섞임. 결과:

- Doppler off, angle `0→10 mrad`: `max|Δχbar|/max|χbar| = 2.942038%`.
- 10 mrad: residual ratio `|Δk|/k = 0.009999958`.
- Doppler on, 50 °C, 3 mm: resonance transmission `0.423098→0.00267946`.

실제 velocity-only `Δk·v`: `v=0`이면 0. Doppler-off angle dependence 없어야 함.
Optical scan-shift operator와 velocity-only residual operator 분리 필요. Dimension,
scan 수, velocity class, solve 수 불변. 계산 overhead 사실상 0.

또 `_two_photon_k_ratio()`는 probe/coupling에 같은 `medium["k_vec"]` 전달
(`gabes/schemes/absorption.py:597-601`; `gabes/beam.py:50-69`). 0 mrad residual을 정확히
0으로 만듦. 실제 co-propagating Λ 두 광장은 ground splitting만큼 주파수 다름. 지원
isotope에서 각 optical frequency 사용: scalar 산술만 추가.

기본 EIT `n_g≈1.00×10^5`를 15 mm에 단순 대입하면 delay 약 5 µs. 하지만 `T=1.36%`,
linewidth 미수렴. usable slow-light delay로 인용 금지.

## 6. 기존 physics 개선안: 비용과 보존성

| 개선 | 물리 보존/변경 | 계산비 | 권고 |
|---|---|---:|---|
| scan/velocity operator 분리 | 잘못 섞인 항만 교정. compact 3-level OBE 유지 | 추가 solve 0 | P1 즉시 |
| clustered 601 + interpolated FWHM/status | 같은 OBE·같은 601점. sampling만 신뢰 가능하게 | paired `+0.2%`, postprocess `O(N)` | P1 즉시 |
| samples/edge/visibility gate + T·delay 병기 | physics 출력 불변. 과신 방지 | `O(N)`, 거의 0 | P1 함께 |
| full probe-minus-coupling `Δk` | momentum bookkeeping 교정 | scalar 산술 | operator 분리 뒤 |
| sourced collision coefficient/pressure shift | 실제 line center·homogeneous width 개선 | lookup `O(1)` | 출처/범위 확보 후 |
| coated-cell bright/dark Ramsey | 새 실제 physics | trajectory/phase 평균, 중~대 | opt-in research |
| VCC kernel | velocity class 결합. separable average 깨짐 | 큰 행렬/iterative solve | dataset 선택 뒤 |
| hyperfine/Zeeman Λ | branching·polarization·pumping 추가 | Liouville 차원 제곱 증가 | held-out data 뒤 |

P1은 physics 보존이 아니라 **오류 수정 + 같은 compact 모델의 해상도 보장**. 비용 거의
없음. Heavy proposal은 실제 physics 늘리지만 overhead 무시 불가. 현 parked 판단 타당.

## 7. 순수 코드 최적화

### 7.1 이미 완료: `line_strength` navigate-only

현재 `line_strength`는 recompute key 아님
(`gabes/schemes/absorption.py:561-562, 668-688`; 시험
`tests/test_absorption.py:110-132`). 같은 `chi_bar`에 readout scale만 재적용.

- headless remap median: `0.1185 ms`.
- 전체 EIT compute median: `33.7423 ms`.
- 비용비: 약 `285×` 가벼움.

행동·물리 정확히 보존. 이미 좋은 최적화. 다만 파일 staged/unstaged 겹침 상태이므로 merge
시 회귀 주의.

### 7.2 새 기존 제안: pole-residue Doppler

`lambda-pole-residue-doppler`: affine velocity dependence를 pole-residue로 바꿔 711개
속도 class 직접 solve 제거 (`docs/checklist.json:2565-2582`).

- 현재 default EIT `33.57 ms`; Doppler-off `≈3.3 ms`. 이 PC에서 이론상 절대 절감 상한
  대략 30 ms/compute.
- 3-level, 9×9 Liouvillian이라 현재 compiled batched solve도 이미 빠름. 절대 이득 작음.
- 현재 discrete 2 m/s Maxwell measure를 `PoleResidues.discrete_mean`로 평가하고 per-row direct
  guard 두면 physics/behavior를 tolerance 안에서 보존 가능.
- Continuous Faddeeva/Gaussian mean 채택 시 4σ cutoff measure 변경. 이것은 순수 최적화 아님.
- 먼저 `lambda-spectrum-validity`가 operator 정의를 바꿈. 그 뒤 port해야 중복 작업 없음.

결론: P3 순서 타당. Parameter sweeps·extra views 추가 전에는 급하지 않음.

### 7.3 더 작은 기회

`_medium_from_params()`는 매 compute마다 species line constant, dipole, mass 재구성
(`gabes/schemes/absorption.py:422-447`). 이번 측정 `17.25 µs/call`; EIT compute의
약 `0.05%`. `(species,line)` immutable 부분 cache 가능. Behavior 보존되지만 우선순위 낮음.

Figure 없는 batch는 이미 `headless_observables()` 사용 가능
(`gabes/schemes/base.py:101-112`; Lambda opt-in `absorption.py:470`). 이것이 plotting skip의
정답. 새 numerical shortcut 불필요.

## 8. 문서·예제·시험 적합성

### 8.1 문서

README는 scheme 존재·핵심 출력·시험 진입점만 짧게 기록
(`README.md:13, 68, 143, 341`). User Guide는 3-level Λ, EIT/AT/CPT, Re χ, FWHM,
group index를 설명 (`docs/Userguide/GABES_User_Guide_v2.html:531-534, 595-613,
849-851`).

부족/낡음:

- User Guide runtime `~0.6 s` (`:781`): 현재 warm EIT `~0.034 s`와 약 18배 차이.
- EIT caption `느린 빛의 서명` (`:604-606`): transmission·grid validity·pulse bandwidth
  경고 없음.
- Lambda limit는 “full polarization/Zeeman 아님” 한 줄뿐 (`:849-851`). Natural-Rb scalar
  effective-medium, symmetric branching, anchored Ωc, linewidth 미수렴, cold-angle bug,
  visibility 미평가 누락.
- `eit.png`, `at.png`: 2026-06-08 생성. CPT 그림 없음. source script·parameter manifest·commit·
  grid status 없음.

### 8.2 예제

`examples/` 디렉터리 없음. Lambda 전용 script/config/CSV example 없음. 가장 가까운 실행
예제: `tests/test_absorption.py:135-241`, `tests/test_schemes_render.py:97-107`.

실험 reference 품질 개선에는 새 solver보다 재현 가능한 headless sweep 예제가 값 큼:
power/diameter, angle, `gamma_gg`, cell length sweep; params JSON, commit, grid/visibility status,
scalar-medium badge 같이 저장. 기존 proposal 범위 안 문서 작업. 새 독립 proposal 불필요.

### 8.3 시험

좋은 범위:

- AT split≈Ωc: `tests/test_absorption.py:135-143`.
- power/diameter Rabi scaling: `:146-154`.
- warm angle broadening: `:157-173`.
- EIT transparency, CPT sub-natural: `:176-200`.
- regime defaults·unresolved status: `:203-241`.
- affine compiled kernel vs NumPy reference: `tests/test_kernels.py:145-191`.
- render cases: `tests/test_schemes_render.py:97-107`.

빈 핵심 시험:

- Doppler-off angle invariance.
- Full optical-frequency `Δk` analytic thermal width.
- Clustered grid vs high-resolution EIT/CPT FWHM·group index.
- samples/FWHM·edge clearance·visibility status.
- natural-Rb scalar-medium qualifier/provenance.

## 9. 검증 결과

```powershell
python -m pytest -q tests/test_absorption.py tests/test_kernels.py `
  tests/test_headless_observables.py tests/test_schemes_render.py `
  tests/test_docs_consistency.py
# 55 passed, 1 failed in 12.29 s
```

실패 1개: Scheme 2 아님. 기존 작업 트리에서 삭제된
`docs/FWM physics and analytic reconstruction/FWM_physics.tex`를
`tests/test_docs_consistency.py:20`이 읽음 → `FileNotFoundError`.

저장소 필수 전체 시험:

```powershell
python -m pytest -q
# collection 중단: 1 error in 21.15 s
# ModuleNotFoundError: No module named 'sympy'
```

원인: ignored reference GUI test가
`references/AutoOD/rb_od_calculator_optimized_final.py:5`의 optional `sympy.physics.wigner`
import 요구. Scheme 2 회귀 실패 아님. 관련 55개 동작 시험 통과.

## 10. 결론·권고 순서

Scheme 2: 장난감 curve 아님. Compact 3-level steady-state OBE, Maxwell average, anchored
Rabi, Beer–Lambert, dispersion 있음. **Resolved AT와 knob trend에는 유용한 실제 물리
reference.** 하지만 현재 EIT hero linewidth 2.89배 과대, cold residual-Doppler 불변성
위반 2.942%, 큰 `n_g`에 낮은 transmission/visibility 경고 없음. 절대 실험 reference 아님.

순서:

1. `lambda-spectrum-validity` 완료: operator 분리, full `Δk`, clustered 601,
   interpolated FWHM, resolution/edge/visibility status, transmission·delay 병기.
2. 위 빈 시험 먼저 추가. 같은 solve 수, overhead 거의 0.
3. User Guide runtime·경고·그림 갱신. CPT + provenance sweep example 추가.
4. Collision coefficient: 출처·부호·단위·validity domain 확보된 값만 `O(1)` 적용.
5. P1 안정화 뒤 `lambda-pole-residue-doppler`. Discrete measure + guard면 순수 가속;
   continuous Gaussian 변경은 별도 physics gate.
6. Coated cell, VCC, hyperfine/Zeeman: geometry·dataset·holdout·runtime budget 준비 뒤 opt-in.

오늘 새 physics proposal 없음. 가장 큰 개선: 기존 compact physics를 같은 비용으로 제대로
샘플링하고, 실험자가 숫자 신뢰 범위를 즉시 보게 만드는 것.
