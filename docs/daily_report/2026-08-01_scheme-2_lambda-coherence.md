# 2026-08-01 Scheme 2 물리 검토 — Lambda coherence

## 1. 오늘의 선택과 현재 다섯 scheme

서울 현지 날짜의 day-of-month는 `1`이므로

```text
n = (1 mod 5) + 1 = 2
```

이다. 실제 dropdown 순서는 `gabes/schemes/__init__.py:19-25`의 `_SCHEMES`가
결정하며, 현재 정의는 다음과 같다. `eit`, `at`, `cpt`는 직접 호출·테스트용
alias일 뿐 dropdown 항목을 추가하지 않는다
(`gabes/schemes/__init__.py:27-39`).

| 순번 | 등록 인스턴스 | 사용자 표시 | 핵심 물리 |
|---:|---|---|---|
| 1 | `SASScheme()` | Absorption spectroscopy (OD / SAS) | pump-off OD, pump-on SAS와 hyperfine optical pumping |
| 2 | `LambdaScheme()` | Lambda coherence (EIT / AT / CPT) | 축약 3준위 Λ 결맞음 |
| 3 | `RydbergEITScheme()` | Rydberg-EIT electrometry | cascade EIT, microwave AT, electrometry |
| 4 | `MagnetoScheme()` | Hanle / EIA / NMOR | Zeeman manifold의 투과·회전 |
| 5 | `FWMScheme()` | FWM (Squeezing / Biphoton) | seeded FWM과 SFWM source estimate |

README의 표도 같은 다섯 항목을 같은 순서로 요약한다
(`README.md:8-16`). 오늘 검토 대상은 2번 `LambdaScheme`이다. 기준 HEAD는
`a82bbf7`, branch는 `main`이다.

## 2. 선행 제안·TODO·issue note 검색 결과

먼저 다음을 검색·대조했다.

- Scheme 2 선행 보고서:
  `docs/daily_report/2026-06-26_scheme-2_lambda-coherence.md`,
  `2026-07-01_scheme-2_lambda-coherence.md`,
  `2026-07-06_scheme-2_lambda-coherence.md`,
  `2026-07-16_scheme-2_lambda-coherence.md`,
  `2026-07-31_scheme-2_lambda-coherence.md`
- 계획/TODO: `docs/checklist.json`, `gabes/constants.py:79-89`
- 구현: `gabes/schemes/absorption.py`, `gabes/atoms.py`, `gabes/beam.py`,
  `gabes/core.py`, `gabes/kernels.py`, `gabes/observables.py`,
  `gabes/lineshape.py`
- 테스트: `tests/test_absorption.py`, `tests/test_kernels.py`,
  `tests/test_headless_observables.py`, `tests/test_schemes_render.py`
- 문서·예제: `README.md`, `docs/Userguide/GABES_User_Guide_v2.html`,
  `docs/Userguide/userguide_assets/eit.png`, `at.png`

별도의 로컬 issue/proposal 파일은 없다. 저장소 내부에서 이 역할은
`docs/checklist.json`, `gabes/constants.py`의 TODO, 일일 보고서가 맡고 있다.
따라서 “기존 개선안이 없음”이 아니라, 이미 다음 제안들이 축적되어 있다.

| 기존 항목 | 현재 상태 | 물리 보존·계산 비용 판단 |
|---|---|---|
| lab-facing `sqrt(P)/d` Rabi | 구현 완료 | solve 전 scalar 환산이라 차원·solve 수가 늘지 않는다 (`gabes/schemes/absorption.py:529-548,583-591`). |
| residual two-photon Doppler / beam angle | 구현 완료 | 기존 Maxwell 평균의 계수만 바꿔 추가 solve가 없다. 다만 아래 cold-limit 오류가 남아 있다 (`docs/checklist.json:109-113`). |
| figureless/headless readout | 구현 완료 | 물리를 바꾸지 않고 Matplotlib만 생략한다 (`docs/checklist.json:57-61`; `tests/test_headless_observables.py:38-67`). |
| gas/species/line별 pressure shift·broadening·Dicke proxy | deferred, GROUP B | 검증 계수표와 scalar shift/유효 폭으로 구현하면 solver 차원은 그대로다 (`docs/checklist.json:49-54`; `gabes/constants.py:79-89`). |
| full velocity-changing collision(VCC) | deferred, GROUP C | velocity class를 서로 결합해 현재의 독립 Maxwell 평균을 깨므로 무겁다 (`docs/checklist.json:129-134`). |
| hyperfine/Zeeman-resolved Λ + optical pumping | deferred, GROUP C | 실험 충실도는 크게 오르지만 Liouvillian 차원과 preset 의미가 바뀐다 (`docs/checklist.json:136-141`). |

7월 31일 검토 이후 HEAD와 Scheme 2 관련 소스·테스트는 바뀌지 않았다. 따라서
오늘의 핵심은 새 기능 평가보다 기존 제안의 유효성 재검증과 문서 진단의 정정이다.

## 3. 현재 구현하는 실제 물리

### 3.1 3준위 정상상태 OBE

`atoms.lambda3()`는 `g1`, `g2`, `e`의 세 상태, 들뜬 상태에서 두 바닥상태로
각각 `Γ/2`의 자발붕괴, 그리고 `ρ_g1g2`의 Raman 결맞음 감쇠를 구성한다
(`gabes/atoms.py:157-179`). `buffer_ground_relax_khz`는 population reload나
확산 모델이 아니라 이 단일 coherence-dephasing rate이다
(`gabes/atoms.py:175-176`).

Hamiltonian은 고정 약한 probe `Ω_p = 10^-3 Γ`, coupling `Ω_c`, coupling
detuning과 two-photon detuning을 담는다
(`gabes/schemes/absorption.py:29-31,633-656`). 그러므로 다음 현상은 실제
Lindblad 정상상태 OBE에서 나온다.

- EIT의 destructive interference와 정상분산
- 강한 coupling에서의 Autler–Townes dressing과 이중선
- 좁은 two-photon dark resonance
- coupling power·beam diameter·ground-coherence `T2`에 대한 올바른 방향성

다만 CPT는 별도의 bichromatic clock-CPT 엔진이 아니라 같은 weak-probe Λ의 좁은
scan preset이다. 두 광장의 균형, 실제 optical pumping, light shift, clock contrast를
절대적으로 예측하지 않는다.

### 3.2 scalar D-line 매질과 실험 노브

species/line 선택은 온도 밀도, 자연선폭, 파장, 질량, reduced dipole을 바꾼다
(`gabes/schemes/absorption.py:418-443,521-526`). coupling은 보정된 Rabi anchor를
`sqrt(power)/diameter`로 옮긴다 (`gabes/beam.py:9-20`). 실험에서 측정한
`Ω_c`를 power/beam-size sweep으로 확장하는 데는 유용하지만 polarization, CG,
mode overlap을 first-principles로 예측하는 환산은 아니다.

특히 `Rb (natural)`은 두 isotope의 분리된 susceptibility 합이 아니다. 가장 풍부한
isotope의 line constants를 택하고 natural-Rb 전체 밀도를 한 scalar Λ에 넣는다
(`gabes/schemes/absorption.py:430-442`). 따라서 isotope-specific Raman resonance나
절대 participating density의 reference로 읽으면 안 된다.

### 3.3 Doppler 평균과 관측량

Doppler-on Λ는 `4σ`, `dv = 2 m/s` Maxwell grid에서 scan×velocity 정상상태를 풀고,
velocity-weighted probe coherence만 축약한다
(`gabes/schemes/absorption.py:118-140`; `gabes/kernels.py:454-495`). beam angle은
`|Δk|/k`를 ground-level Doppler ratio에 넣는다
(`gabes/beam.py:42-55`; `gabes/atoms.py:165-178`).

`χ_bar`는 density와 dipole을 거쳐 `χ_phys`, `α = k Im(χ)`,
`T = exp(-αL)`로 변환된다 (`gabes/observables.py:393-418`). group index는
`Re χ`의 국소 수치 미분이다 (`gabes/observables.py:426-434`). 이는 실제
선형응답 물리이지만 finite-pulse propagation이나 검출 SNR 계산은 아니다.

## 4. 실험 원자물리 reference로서의 평가

### 유용한 범위

현재 구현은 **실제 물리를 가진 semi-quantitative 실험 계획 reference**로 유용하다.

- coupling power/diameter를 바꿀 때 EIT에서 AT로 넘어가는 scale
- AT split이 보정된 `Ω_c`를 따르는지 확인
- Raman `T2`가 EIT/CPT 창을 흐리는 방향과 규모
- mrad 정렬 오차가 warm-vapor transparency를 씻어내는 경향
- 같은 susceptibility에서 cell length·유효 line strength를 바꾼 transmission

현 HEAD의 기본 출력과 warm solve median은 다음과 같다.

| regime | compute median | 현재 hero/readout |
|---|---:|---|
| EIT | `0.174 s` | `T_res = 0.014`, FWHM `0.46 MHz`, `n_g = 1.003e5` |
| AT | `0.0046 s` | split `46.0 MHz`, expected `46.0 MHz`, `T_center = 0.989` |
| CPT | `0.0042 s` | `T_res = 0.706`, FWHM `923.32 kHz`, `n_g = 9.525e4` |

AT split, `sqrt(P)/d`, warm-angle broadening, cold EIT, sub-natural CPT,
Numba/NumPy parity, unresolved-feature status는 테스트되어 있다
(`tests/test_absorption.py:110-216`; `tests/test_kernels.py:102-136`).

### 절대값 reference로 부적합한 범위

다음에는 아직 논문급 절대 reference가 아니다.

- natural-Rb isotope/hyperfine line assignment와 participating density
- polarization·Zeeman·CG branching·optical pumping에 따른 contrast
- buffer-cell pressure shift, homogeneous broadening, diffusion, VCC
- clock-CPT light shift·균형 광장·절대 linewidth/contrast
- finite-pulse delay·distortion·bandwidth와 실제 detector SNR

따라서 linewidth·contrast·group index의 경향과 order of magnitude에는 쓸 수 있지만,
실험 피팅 파라미터나 불확도 없는 절대 예측치로 인용해서는 안 된다.

## 5. 재확인된 핵심 문제와 저비용 개선

### P0 — EIT FWHM과 group index의 scan-grid 미수렴

모든 regime의 `_scan()`은 601점을 사용한다
(`gabes/schemes/absorption.py:609-631`). `window_fwhm()`은 half-height 경계를
sample 단위로 걷고 보간하지 않는다 (`gabes/lineshape.py:19-37`). 기본 EIT의
표시 FWHM `0.45968 MHz`는 격자의 정확히 두 칸이다.

| uniform points | step | sample FWHM | edge-interpolated FWHM | center `n_g` | 오늘 1회 solve |
|---:|---:|---:|---:|---:|---:|
| 601 | 0.22984 MHz | 0.45968 MHz | 0.22990 MHz | 100,305 | 0.097 s |
| 1,201 | 0.11492 MHz | 0.22984 MHz | 0.14169 MHz | 116,233 | 0.373 s |
| 2,401 | 0.05746 MHz | 0.22984 MHz | 0.16174 MHz | 113,989 | 0.692 s |
| 4,801 | 0.02873 MHz | 0.17238 MHz | 0.15884 MHz | 112,966 | 1.427 s |
| 9,601 | 0.01436 MHz | 0.17238 MHz | 0.15910 MHz | 112,693 | 2.035 s |

현재 hero FWHM은 수렴 보간값보다 약 `2.9x` 크고 `n_g`는 약 `11%` 작다. CPT도
현재 sample 값 `923.32 kHz`에 대해 9,601점 보간 수렴값은 `803.54 kHz`이다.

개선 비용은 매우 작게 유지할 수 있다.

1. `samples_per_width`, `resolution-limited` status와 half-height edge interpolation은
   기존 배열의 `O(N)` 후처리라 추가 OBE solve가 없다.
2. EIT/CPT에 601점 center-clustered grid를 쓰면 solve 수는 그대로다. 오늘의
   `sinh(5u)` 진단은 보간 FWHM `0.15914 MHz`, `n_g = 112,708`로 9,601점 결과와
   각각 약 `0.03%`, `0.01%` 이내였다. 넓은 wing 회귀는 추가해야 한다.
3. 9,601점 uniform scan은 약 2초여서 publication/reference용 opt-in으로 두는 편이
   낫다.

### P0 — residual Doppler가 cold spectrum도 바꾸는 계수 오류

현재 `_affine_scan_coeffs()`는

```text
A_coef = dL/ds - S_v
B_coef = S_v
```

를 만든다 (`gabes/schemes/absorption.py:45-57`). residual ratio가 ground level까지
`S_v`에 들어가므로 velocity가 0이어도 scan coefficient에 angle-dependent 항이
섞인다. 물리적인 residual Doppler는 `Δk·v`이므로 `doppler=off`에서 angle은
스펙트럼을 바꾸지 않아야 한다.

현 HEAD의 cold EIT 0 대 10 mrad 재검증 결과는 다음과 같다.

- `max |Δχ| / max |χ| = 2.9420e-2`
- `cell_mm = 3`에서 최대 transmission 차이 `2.9769e-3`

scan operator와 velocity operator를 분리해 residual ratio를 `kv` 항에만 적용하면
solver 차원·velocity 수·solve 횟수가 모두 같다. 즉 물리를 바로잡으면서 runtime
증가는 사실상 없다. `doppler=off` angle invariance 회귀 테스트가 반드시 필요하다.

### P1 — 수치 feature와 실험 검출 가능성을 구분하지 않음

현재 readout은 half-height crossing 또는 두 local maximum이 아예 없을 때만
`unresolved`를 낸다 (`gabes/schemes/absorption.py:723-795`). 절대 contrast나
detector noise floor는 보지 않는다.

허용 범위의 `20 °C`, `0.5 mm`, `line_strength = 0.01`, **`doppler=on`** EIT에서

```text
T_min = 0.9994266174
T_max = 0.9999261407
modeled contrast = 4.99523e-4
```

인데도 hero는 `T_res = 1.000`, `FWHM = 0.92 MHz`다. 7월 31일 보고서가 이 사례를
`doppler=off`라고 적은 것은 진단 조건 표기 오류이며, 오늘 바로잡는다. 실제
`doppler=off`에서는 넓은 cold absorption 구조 때문에 전역 transmission range가
약 `3.27%`이다.

modeled contrast, samples-per-width, optional detector threshold를 기존 배열에서
계산하면 `O(N)` 후처리일 뿐 추가 solve가 없다. noise 정보를 모르면
`numerically resolved; experimental visibility not assessed`라고 표시하는 것만으로도
오용을 줄일 수 있다.

### P1 — group index가 loss 없는 slow-light 예측처럼 보임

`group_index()`는 local dilute-medium derivative이고 pulse 전파·왜곡은 풀지 않는다
(`gabes/observables.py:426-434`). 기본 EIT의 수렴값은 `n_g ≈ 1.127e5`지만 resonance
transmission은 약 `1.36%`에 불과하다. 큰 `n_g`만으로 좋은 delay line이라 판단할 수
없다.

추가 solve 없이 group delay, delay-bandwidth product, `T_res`를 함께 보여 주고
finite-pulse 예측이 아님을 명시해야 한다. 문서상 `README.md:18-20`은 group-index
readout을 아직 roadmap으로 적지만, 실제 hero는 이미 이를 출력하므로 정정이 필요하다.

### P1 — rate와 매질 의미가 과도하게 축약됨

- `buffer_ground_relax_khz`라는 이름은 buffer collision budget처럼 보이지만 실제로는
  단일 Raman `T2` dephasing이다 (`gabes/schemes/absorption.py:549-553`;
  `gabes/atoms.py:175-176`). `ground-coherence dephasing / 2π`로 이름을 바꾸고,
  transit/diffusion/spin-destruction budget은 검증 계수가 있을 때 분리해야 한다.
  같은 3준위 dissipator의 scalar rate 합이므로 solve 차원은 늘지 않는다.
- natural-Rb가 effective scalar medium임을 UI/Derived table에 경고하는 것은 runtime 0의
  신뢰도 개선이다. 실제 isotope 합을 원하면 hyperfine/Zeeman reference mode가 필요하다.
- pressure shift·추가 homogeneous width·제한적 Dicke proxy는 검증된 coefficient table을
  전제로 같은 차원에서 구현 가능하다. 반면 full VCC와 full manifold는 저비용 수정이
  아니므로 별도 opt-in heavy solver가 맞다.

## 6. 동작을 바꾸지 않는 순수 코드 최적화

### 6.1 `line_strength`를 navigate-only로 이동

`line_strength`는 현재 기본 `recompute=True`를 상속한다
(`gabes/schemes/absorption.py:557-558`). 그러나 Hamiltonian이나 `chi_bar`에는 들어가지
않고, solve 뒤 `χ_phys`를 곱할 때만 쓰인다
(`gabes/schemes/absorption.py:664-668,677-684`;
`gabes/observables.py:393-413`).

오늘 기본 EIT에서 line-strength-only 재 solve median은 `0.172 s`, 같은 raw에서
전체 headless metric/table remap은 `0.193 ms`였다. 약 `890x` 차이다.

`ParamSpec(..., recompute=False)`로 바꾸고 `observables()`가 live
`params["line_strength"]`를 읽게 하면 된다. 여러 값에서 `chi_bar` bitwise 불변,
`alpha` 선형 scaling, transmission equivalence를 회귀하면 behavior-preserving이다.

### 6.2 affine kernel의 scan-constant matrix hoist

현재 kernel은 모든 `(scan, velocity)` 쌍에서 9×9 전체
`base + s*A_coef + kv*B_coef`를 다시 조립한 뒤 trace row를 덮어쓴다
(`gabes/kernels.py:426-449`). `base + s*A_coef`를 scan당 한 번 만들고 inner loop에서는
copy와 `kv*B_coef`만 더하며, 어차피 덮는 trace row는 조립하지 않을 수 있다.

7월 31일의 동일 HEAD 진단은 `0.1344 -> 0.1193 s`, `1.13x`, `chi` max absolute
difference `0`이었다. 기존 Numba/NumPy `<1e-11` parity와 Doppler on/off 회귀를
유지하면 물리를 바꾸지 않는 최적화다.

### 6.3 우선순위가 낮은 항목

- `_medium_from_params()`와 immutable `atoms.lambda3()` template cache는 반복 batch에는
  도움이 되지만 각각 수십~수백 μs 규모라 0.1초대 warm EIT solve보다 작다.
- headless readout은 이미 구현되었다. 자동 sweep·보고서는 계속
  `headless_observables()`를 사용해야 한다.
- figure 생성은 계산 물리와 별개로 느리므로 batch에서 Matplotlib을 호출하지 않는 것이
  object micro-cache보다 효과가 크다.

## 7. 문서·테스트·예제의 reference 품질

- README의 Scheme 2 한 줄은 기능을 정확히 요약하지만 scalar/weak-probe/natural-Rb
  effective-density 한계를 말하지 않는다 (`README.md:13`). 또한 이미 구현된 group-index를
  roadmap으로 적는다 (`README.md:18-20`).
- 사용자 가이드는 3준위 Λ와 full polarization/Zeeman 부재를 밝힌다
  (`docs/Userguide/GABES_User_Guide_v2.html:530-533,801-802`). 다만 group index를
  “느린 빛의 서명”으로만 설명하고 transmission loss·finite-pulse 한계를 함께 보이지
  않는다 (`docs/Userguide/GABES_User_Guide_v2.html:577-591`).
- `eit.png`는 제목에 `Ω_c = 6.00 MHz`, `gamma_gg = 10.0 kHz`를 표시하지만 현재 EIT
  default는 `17.238 MHz`, `57.46 kHz`이다
  (`gabes/schemes/absorption.py:482-490`). 의도적인 비기본 예제일 수 있으나 온도·셀 길이와
  생성 script/config가 없어 exact reproduction provenance가 없다. `at.png`의
  `45.97 MHz`, `57.5 kHz`는 현재 AT default와 일치한다.
- 실행 가능한 Lambda 전용 example은 없다. 현재 closest example은
  `tests/test_absorption.py:110-216`의 direct API 호출이다. 짧은 headless
  power/diameter·angle·dephasing sweep과 CSV/Markdown provenance 예제를 추가하면 runtime
  물리를 건드리지 않고 실험 reference 재현성을 높일 수 있다.
- 테스트는 textbook invariant와 UI contract에는 좋다. 그러나 cold-angle invariance,
  FWHM/grid convergence, natural-Rb isotope semantics, contrast/detectability threshold,
  group-delay/loss status는 아직 고정하지 않는다.

## 8. 검증

오늘 현 HEAD에서 직접 실행했다.

- 현재 registry와 EIT/AT/CPT 기본 solve·headless metrics
- 601–9,601점 uniform resolution sweep과 601점 center-clustered 진단
- cold 0/10 mrad angle invariance 진단
- 허용 범위 저대비 warm-EIT 진단
- `line_strength` 재 solve 대 headless remap benchmark
- 관련 테스트와 전체 test suite

```powershell
python -m pytest tests/test_absorption.py tests/test_kernels.py tests/test_headless_observables.py tests/test_schemes_render.py -q
# 52 passed in 47.14 s

python -m pytest -q
# 232 passed in 125.36 s
```

초기에 짧은 command timeout으로 중단된 두 pytest 시도는 완전한 명령으로 다시 실행했고,
위 두 최종 실행은 모두 통과했다.

## 9. 최종 판단과 우선순위

Scheme 2는 실제 3준위 Lindblad OBE, Maxwell velocity average, lab-facing Rabi scaling,
Beer–Lambert propagation을 결합하므로 장난감 곡선 생성기가 아니다. **EIT/AT/CPT의
scale과 knob sensitivity를 빠르게 보는 반정량 실험 reference**로는 충분히 가치가 있다.

그러나 현재 가장 먼저 고칠 것은 더 큰 atomic manifold가 아니라 다음 두 P0이다.

1. 같은 601 solve를 center-clustered grid에 배치하고 edge interpolation·resolution
   status를 추가해 FWHM과 `n_g`의 수렴성을 회복한다.
2. residual Doppler의 scan/velocity coefficient를 분리하고 cold-angle invariance를
   회귀한다.

둘 다 solver 차원과 solve 수를 늘리지 않거나 그대로 유지한다. 그 다음은 modeled
contrast와 experimental visibility의 구분, natural-Rb scalar 경고, Raman-dephasing rate
명확화, group-index와 loss의 동시 표시가 우선이다. full VCC와 hyperfine/Zeeman Λ는
실험적 가치가 크지만 별도 검증 자료와 opt-in heavy-mode 설계가 필요한 장기 항목이다.
