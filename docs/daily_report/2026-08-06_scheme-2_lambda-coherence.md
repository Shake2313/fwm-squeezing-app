# 2026-08-06 Scheme 2 물리 검토 — Lambda coherence

## 1. 선택 규칙과 현재 다섯 scheme

현지 날짜는 2026-08-06이다.

```text
n = (day mod 5) + 1 = (6 mod 5) + 1 = 2
```

따라서 오늘의 대상은 **Scheme 2, `LambdaScheme()` — Λ coherence
(EIT / AT / CPT)** 이다. 드롭다운 순서는 `gabes/schemes/__init__.py:12-24`의
`_SCHEMES`가 결정하며 현재 정의는 다음과 같다.

| 순번 | 등록 클래스 | UI scheme | 핵심 물리 |
|---:|---|---|---|
| 1 | `SASScheme()` | OD / SAS | Doppler 흡수, 포화흡수, hyperfine pumping |
| 2 | `LambdaScheme()` | Λ coherence (EIT / AT / CPT) | 축약 3준위 Λ 결맞음 |
| 3 | `RydbergEITScheme()` | Rydberg-EIT electrometry | cascade EIT, microwave AT, electrometry |
| 4 | `MagnetoScheme()` | Hanle / EIA / NMOR | Zeeman OBE, 투과·회전 readout |
| 5 | `FWMScheme()` | FWM | seeded squeezing, SFWM biphoton |

`eit`, `at`, `cpt` 인스턴스는 직접 호출과 테스트를 위한 alias이며 드롭다운을
늘리지 않는다 (`gabes/schemes/__init__.py:27-34`). README의 사용자-facing 표도
같은 다섯 엔진을 설명한다 (`README.md:8-16`). 기준 HEAD는
`f3898a57629cb3b7d38f42e5cbbe0e819d9edaaf`이다. 작업 트리는 기존 사용자
변경을 포함해 dirty이지만, Scheme 2 소스와 관련 테스트에는 2026-08-01 검토 후
커밋 또는 로컬 diff가 없다.

## 2. 선행 제안·TODO·issue note 검색

먼저 다음 위치를 검색했다.

- 기존 Scheme 2 보고서:
  `docs/daily_report/2026-06-26_scheme-2_lambda-coherence.md`,
  `2026-07-01_scheme-2_lambda-coherence.md`,
  `2026-07-06_scheme-2_lambda-coherence.md`,
  `2026-07-16_scheme-2_lambda-coherence.md`,
  `2026-07-31_scheme-2_lambda-coherence.md`,
  `2026-08-01_scheme-2_lambda-coherence.md`
- 계획/TODO: `docs/checklist.json`, `gabes/constants.py:79-89`
- 현재 설명: `README.md`, `CLAUDE.md`,
  `docs/Userguide/GABES_User_Guide_v2.html`
- 테스트와 예제 후보: `tests/`, `analysis/`, `docs/Userguide/userguide_assets/`

별도 로컬 issue/proposal 파일이나 Lambda 전용 실행 예제는 없다. 이 저장소에서는
`docs/checklist.json`, 코드 TODO, 일일 보고서가 그 역할을 한다. 즉 기존 개선안이
없는 것이 아니라 다음 제안이 이미 축적되어 있다.

| 기존 항목 | 상태 | 계산비용·물리 보존 판단 |
|---|---|---|
| lab-facing `sqrt(P)/d` Rabi | 구현 완료 | scalar Rabi 환산만 추가하며 solve 차원은 같다 (`absorption.py:529-548,583-591`). |
| residual two-photon Doppler / beam angle | 부분적으로 유효 | 같은 Maxwell 계수만 바꾸므로 추가 solve가 없다. 다만 cold-limit 오류와 아래의 종방향 `Δk` 누락이 남는다 (`docs/checklist.json:109-113`). |
| figureless/headless readout | 구현 완료 | Matplotlib만 생략하며 물리는 동일하다 (`docs/checklist.json:57-61`; `tests/test_headless_observables.py:38-67`). |
| gas/species/line pressure shift·broadening·Dicke proxy | deferred, GROUP B | 검증된 scalar shift/width로 제한하면 solver 차원은 그대로다 (`docs/checklist.json:49-54`; `constants.py:79-89`). |
| full velocity-changing collision(VCC) | deferred, GROUP C | velocity class를 서로 결합하므로 현재의 독립 Maxwell 평균보다 훨씬 무겁다 (`docs/checklist.json:129-134`). |
| hyperfine/Zeeman-resolved Λ + optical pumping | deferred, GROUP C | 상태 수와 Liouvillian 차원을 크게 늘리고 preset 의미도 바꾼다 (`docs/checklist.json:137-141`). |

2026-08-01 보고서의 수치·제안은 현재 소스에서도 대부분 재현됐다. 아래에서는
우선순위를 현 checklist 정책의 P1/P2/P3로 표기한다.

## 3. 현재 구현하는 실제 물리

### 3.1 3준위 정상상태 Lindblad OBE

`atoms.lambda3()`는 `g1`, `g2`, `e` 세 상태를 만들고 들뜬 상태가 두
바닥상태로 각각 `Γ/2`로 붕괴하게 한다. `gamma_gg`는 두 바닥상태 사이의
Raman coherence를 직접 감쇠시킨다 (`gabes/atoms.py:157-179`). 이는 장난감
선모양 합이 아니라 trace-normalized 정상상태 density matrix solve이다.

Hamiltonian에는 고정된 약한 probe `Ω_p = 10^-3 Γ`, coupling `Ω_c`, coupling
detuning과 probe/two-photon scan이 들어간다
(`gabes/schemes/absorption.py:29-31,633-663`). 따라서 다음 현상은 실제 OBE에서
나온다.

- EIT destructive interference와 정상분산
- 강한 coupling에서의 Autler–Townes dressing과 이중선
- Raman coherence가 만드는 CPT dark resonance
- coupling power, beam diameter, ground-coherence dephasing의 올바른 경향

다만 CPT는 별도의 균형 bichromatic clock-CPT 엔진이 아니라 같은 weak-probe
Λ의 좁은 scan preset이다. 실제 clock contrast, light shift, 편광 선택규칙,
광펌핑을 절대적으로 예측하지 않는다.

### 3.2 scalar D-line 매질과 실험 노브

Rb/Cs 및 D1/D2 선택은 온도별 원자수밀도, 자연선폭, 질량, 파장, reduced
dipole을 바꾼다 (`absorption.py:418-443,521-526`). coupling은 보정된 Rabi
anchor를 `sqrt(power)/diameter`로 스케일한다 (`absorption.py:583-591`;
`gabes/beam.py:9-20`). 실험에서 한 점의 측정 `Ω_c`를 power/beam-size sweep으로
확장하기에는 유용하지만, CG 계수·편광·공간 mode overlap을 first-principles로
계산하는 절대 power-to-Rabi 모델은 아니다.

특히 `Rb (natural)`은 두 isotope susceptibility를 합산하지 않는다. 가장 큰
abundance의 isotope가 line constant를 제공하고, natural-Rb 전체 밀도를 한
scalar Λ에 넣는다 (`absorption.py:430-442`). 따라서 특정 isotope Raman 선의
절대 participating density로 읽으면 안 된다.

### 3.3 Doppler 평균과 관측량

Doppler-on Λ는 `4σ`, `dv=2 m/s` Maxwell velocity grid에서 scan×velocity
정상상태를 풀고, velocity-weighted probe coherence만 반환한다
(`absorption.py:118-140`; `gabes/kernels.py:454-495`). beam angle은
`|Δk|/k`를 ground-level velocity shift에 추가한다
(`gabes/beam.py:42-55`; `gabes/atoms.py:165-178`).

그 뒤 `χ_bar`는 density와 dipole을 포함한 `χ_phys`,
`α = k Im(χ)`, `T = exp(-αL)`로 변환된다
(`gabes/observables.py:393-418`). group index는 `Re χ`의 국소 수치 미분이다
(`gabes/observables.py:426-434`). 이는 정상상태 선형응답 물리이지만 finite-pulse
propagation, distortion, detector SNR을 풀지는 않는다.

## 4. 2026-08-06 직접 재현값

Numba warm-up 뒤 현재 기본값을 3회 계산한 median과 headless readout은 다음과
같다. 절대 시간은 이 실행 환경의 참고값이다.

| regime | compute median | 현재 주요 readout |
|---|---:|---|
| EIT | `0.1106 s` | `T_res=0.01356`, 표시 FWHM `0.45968 MHz`, `n_g=100,305` |
| AT | `0.00705 s` | split `46.0 MHz`, expected `≈46.0 MHz`, `T_center=0.98909` |
| CPT | `0.00932 s` | `T_res=0.70633`, 표시 FWHM `923.32 kHz`, `n_g=95,255` |

AT split, `sqrt(P)/d`, warm angle broadening, cold EIT, sub-natural CPT,
Numba/NumPy parity, unresolved-feature status가 테스트된다
(`tests/test_absorption.py:110-216`; `tests/test_kernels.py:102-136`).

## 5. 실험 원자물리 reference로서의 판단

### 유용한 범위

현재 Scheme 2는 **실제 물리를 가진 반정량적 실험 계획 reference**이다.

- coupling 세기를 바꿀 때 EIT에서 AT로 넘어가는 scale
- AT splitting이 보정한 `Ω_c`를 따르는지 확인
- Raman `T2`가 EIT/CPT window를 흐리는 방향과 대략적 크기
- mrad 정렬 오차가 warm-vapor transparency를 씻어내는 경향
- 같은 susceptibility에서 cell length와 유효 line strength를 바꾼 transmission

이 용도에서는 compact 3-level 모델의 단순성이 오히려 장점이다. 노브와 결과의
인과관계가 분명하고 한 sweep이 빠르다.

### 절대 reference로 부적합한 범위

다음에는 아직 논문급 절대 기준으로 쓰면 안 된다.

- natural-Rb isotope/hyperfine line assignment와 participating density
- polarization, Zeeman, CG branching, optical pumping에 따른 contrast
- buffer-cell pressure shift, 추가 homogeneous broadening, diffusion, VCC
- clock-CPT light shift, 균형 광장, 절대 linewidth/contrast
- finite-pulse delay·distortion·bandwidth와 실제 detector visibility

따라서 현재 linewidth, contrast, group index는 경향과 order of magnitude에는
쓸 수 있지만 실험 세팅의 불확도 없는 예측치로 인용해서는 안 된다.

## 6. 핵심 문제와 저비용 개선

### P1 — 기본 EIT FWHM과 group index가 scan grid에 미수렴

모든 regime은 601점 uniform scan을 쓴다 (`absorption.py:609-631`).
`window_fwhm()`은 half-height 경계를 sample 단위로 걷고 보간하지 않는다
(`gabes/lineshape.py:19-37`). 기본 EIT 표시값 `0.45968 MHz`는 격자의 정확히
두 칸이다.

오늘의 재진단 결과:

| scan | 중앙 step | sample FWHM | edge-interpolated FWHM | center `n_g` | solve |
|---|---:|---:|---:|---:|---:|
| uniform 601 | `0.22984 MHz` | `0.45968 MHz` | `0.22990 MHz` | `100,305` | `0.126 s` |
| center-clustered 601, `sinh(5u)` | `0.01549 MHz` | `0.18616 MHz` | `0.15914 MHz` | `112,708` | `0.122 s` |
| uniform 9,601 | `0.01437 MHz` | `0.17238 MHz` | `0.15910 MHz` | `112,693` | `1.700 s` |

즉 현재 hero FWHM은 수렴 보간값보다 약 `2.89×` 크고 `n_g`는 약 `11%`
작다. CPT도 601점의 sample 값 `923.32 kHz`가 같은 배열의 edge 보간만으로
`805.73 kHz`가 된다.

물리를 보존하면서 비용을 거의 0으로 유지할 수 있다.

1. `samples_per_width`, `resolution-limited` status와 half-height edge
   interpolation은 기존 배열의 `O(N)` 후처리라 추가 OBE solve가 없다.
2. EIT/CPT는 같은 601점을 중앙에 모으면 solve 수와 오늘 측정 runtime이
   사실상 같다. broad wing·off-center detuning 회귀를 추가해야 한다.
3. 9,601점 uniform scan은 느린 opt-in validation/reference mode로만 두는 편이 낫다.

### P1 — residual Doppler 계수가 cold spectrum까지 바꾸는 오류

현재 affine 조립은

```text
A_coef = dL/ds - S_v
B_coef = S_v
L(s,kv) = base + s*A_coef + kv*B_coef
```

를 쓴다 (`absorption.py:45-57`; `kernels.py:454-482`). residual ratio가
ground-level `S_v`에 들어가므로 velocity가 0인 Doppler-off solve에서도
scan coefficient가 angle-dependent가 된다. 물리적인 residual Doppler는
`Δk·v`이므로 `v=0`에서는 angle이 spectrum을 바꾸면 안 된다.

오늘 cold EIT, 0 대 10 mrad 비교에서 다시 확인했다.

- `max |Δχ_bar| / max |χ_bar| = 2.9420×10^-2`
- 3 mm cell에서 `max |ΔT| = 7.8849×10^-3`

scan optical-shift operator와 velocity-only residual operator를 분리하면 solver
차원, scan 점수, velocity solve 수가 모두 같다. 따라서 정확성을 회복하면서
runtime을 늘리지 않는 P1 수정이다. `doppler=off` angle invariance와
Numba/NumPy parity 회귀 테스트가 필요하다.

### P2 — 같은 진행 방향의 probe–coupling 주파수 차 `Δk`가 0으로 고정됨

기존 beam-angle 제안의 범위에서 새로 확인한 누락이다. 현재
`_two_photon_k_ratio()`는 probe와 coupling에 같은 `medium["k_vec"]`를 넘긴다
(`absorption.py:594-597`). 따라서 `angle=0`이면 helper가 정확히 0을 반환한다
(`gabes/beam.py:36-55`). 그러나 hyperfine Λ의 두 광장은 같은 D line 안에서도
ground-state splitting만큼 주파수가 달라 종방향 `Δk`가 남는다.

`species.py:154-180`의 hyperfine A와 질량으로 계산한 thermal residual Gaussian
FWHM은 다음과 같다.

| isotope | ground split | `|Δk|/k` | 25 °C residual FWHM |
|---|---:|---:|---:|
| ⁸⁵Rb | `3035.732 MHz` | `8.05×10^-6` | `4.07 kHz` |
| ⁸⁷Rb | `6834.683 MHz` | `1.81×10^-5` | `9.07 kHz` |
| ¹³³Cs | `9192.632 MHz` | `2.74×10^-5` | `9.86 kHz` |

현재 기본 EIT/CPT 폭에는 작은 보정이지만 저전력·장수명 clock-CPT에서는
무시하기 어렵다. pure-isotope mode에서는 line frequency에 ground splitting을
더해 실제 두 wavevector를 만들고, natural-Rb scalar mode에는 isotope ambiguity
status를 표시할 수 있다. 기존 velocity coefficient에 상수만 더하므로 추가 solve는
없다. 단, 앞의 cold-limit coefficient 분리를 먼저 해야 이 보정이 Doppler-off
spectrum을 오염시키지 않는다.

### P2 — 수치적으로 보이는 feature와 실험 검출 가능성을 구분하지 않음

현재 readout은 half-height crossing 또는 두 local maximum이 없을 때만
`unresolved`를 표시한다 (`absorption.py:723-795`). modeled contrast나 detector
noise floor는 보지 않는다. 2026-08-01의 허용 범위 저대비 예에서는
`T_max-T_min ≈ 5.0×10^-4`인데도 FWHM이 hero가 됐다.

modeled contrast, samples-per-width, optional detector threshold를 `O(N)`
후처리로 추가할 수 있다. noise 정보가 없으면 “numerically resolved;
experimental visibility not assessed”라고 표시하는 것만으로도 과신을 줄인다.

### P2 — group index가 loss 없는 slow-light 예측처럼 보임

`group_index()`는 dilute-medium local derivative이고 pulse propagation·distortion을
풀지 않는다 (`observables.py:426-434`). 수렴된 기본 EIT는 `n_g≈1.127×10^5`로
크지만 resonance transmission은 `1.36%`뿐이다. `n_g`만으로 좋은 delay line이라고
판단할 수 없다.

추가 solve 없이 `T_res`, group delay `(n_g-1)L/c`, linewidth 기반의 제한적
delay-bandwidth proxy를 함께 보이고 finite-pulse 예측이 아님을 명시해야 한다.
실제 pulse propagation은 별도 opt-in 계산이어야 한다.

### P2 — rate와 매질 의미가 지나치게 압축됨

- `buffer_ground_relax_khz`는 이름과 help상 buffer relaxation budget처럼 보이지만
  실제로는 population reload가 없는 단일 Raman-coherence dephasing이다
  (`absorption.py:549-553`; `atoms.py:175-176`). 이름을
  `ground-coherence dephasing / 2π`로 명확히 하는 것은 0-cost 개선이다.
- natural-Rb가 effective scalar medium임을 UI와 derived table에 경고하는 것도
  runtime 0의 신뢰도 개선이다.
- pressure shift, 추가 homogeneous width, 검증 범위가 명시된 Dicke proxy는 같은
  solve 차원에서 가능하다. full VCC는 분리된 heavy reference solver가 맞다.

## 7. 동작을 바꾸지 않는 순수 코드 최적화

### 7.1 `line_strength`를 navigate-only로 이동

`line_strength`는 현재 기본 `recompute=True`를 상속한다
(`absorption.py:557-558`). 그러나 Hamiltonian과 `chi_bar`에는 들어가지 않고,
solve 뒤 `χ_phys`를 스케일할 때만 쓰인다
(`absorption.py:664-668,677-684`; `observables.py:393-413`).

오늘 기본 EIT에서:

- line-strength-only 재 solve median: `0.12289 s`
- 같은 raw의 headless remap median: `0.1783 ms`
- 차이: 약 `689×`

현재는 observables가 `raw["ls"]`를 읽기 때문에 params만 바꾸면 결과가 전혀
바뀌지 않는 것도 확인했다. 따라서 다음 두 변경을 함께 해야 한다.

1. `ParamSpec(..., recompute=False)`
2. `observables()`에서 live `params["line_strength"]` 사용

여러 line-strength 값에서 `chi_bar` bitwise 불변, `alpha` 선형 scaling,
transmission/metric 동등성을 테스트하면 물리와 사용자-visible 동작을 그대로
보존할 수 있다.

### 7.2 affine kernel의 scan-constant matrix hoist

현재 kernel inner velocity loop는 매번 전체
`base + s*A_coef + kv*B_coef`를 조립한다 (`gabes/kernels.py:426-449`).
`base + s*A_coef`를 scan마다 한 번 계산하고, velocity loop에서는 copy와
`kv*B_coef`만 더할 수 있다. 2026-08-01 동일 소스 진단은 `1.13×` 개선과
`chi` max absolute difference `0`을 보였다. 현재 관련 소스가 변하지 않았으므로
여전히 적용 가능하다. Numba/NumPy parity와 Doppler on/off 회귀가 안전장치다.

### 7.3 낮은 우선순위

- `_medium_from_params()`와 immutable `atoms.lambda3()` template cache는 반복
  batch에는 도움이 되지만 생성비가 warm EIT solve보다 훨씬 작다.
- headless readout은 이미 구현됐다. 자동 sweep과 보고서는 계속
  `headless_observables()`를 써야 한다.
- figure 생성을 생략하는 효과가 object micro-cache보다 크므로 Matplotlib을
  batch path에서 호출하지 않는 현재 구조를 유지해야 한다.

## 8. 문서·테스트·예제 평가

### 문서

- README의 Scheme 2 한 줄은 기능을 정확히 요약하지만 scalar/weak-probe/
  natural-Rb effective-density 한계를 말하지 않는다 (`README.md:13`).
- 이미 구현된 group-index readout이 아직 roadmap에 있다 (`README.md:18-20`;
  실제 출력은 `absorption.py:771-787`).
- 사용자 가이드는 3준위 Λ와 full polarization/Zeeman 모델의 차이를 명시한다
  (`docs/Userguide/GABES_User_Guide_v2.html:530-533,800-802`). 반면 EIT 그림을
  “느린 빛의 서명”이라고 설명하면서 큰 loss와 finite-pulse 한계를 함께 보여주지
  않는다 (`...:575-591`).
- module docstring은 노브가 Γ 단위라고 적지만 UI는 `physical-units-v3`의
  MHz/kHz 노브다 (`absorption.py:16-18,462-465`). P3 문서 정리 대상이다.

### 예제

`docs/Userguide/userguide_assets/eit.png`와 `at.png`는 유용한 시각 예제다.
그러나 Lambda 전용 생성 script/config가 없어 exact reproduction provenance가
없다. `eit.png`에 표시된 `Ω_c=6.00 MHz`, `gamma_gg=10.0 kHz`는 현재 EIT
default `≈17.24 MHz`, `≈57.46 kHz`와 다르다. 실행 가능한 가장 가까운 예는
`tests/test_absorption.py:110-216`의 direct API 호출이다.

낮은 비용으로 headless power/diameter, angle, dephasing, cell-length sweep을
CSV/Markdown으로 저장하는 Lambda 예제를 추가하면 실험 reference 재현성이
크게 좋아진다. 각 output에는 parameter manifest, 코드 commit, grid/status를
포함해야 한다.

### 테스트

현재 테스트는 textbook invariant와 UI contract에는 강하다. 다만 다음은 아직
고정하지 않는다.

- Doppler-off angle invariance
- FWHM/grid convergence와 edge interpolation
- collinear hyperfine `Δk`
- natural-Rb scalar semantics
- modeled contrast/detectability status
- group delay와 loss의 동시 표시
- Lambda example artifact의 provenance

## 9. 검증

오늘 current working tree에서 직접 실행했다.

```powershell
python -m pytest tests/test_absorption.py tests/test_kernels.py tests/test_headless_observables.py tests/test_schemes_render.py -q
# 52 passed in 23.41s

python -m pytest -q
# 313 passed in 156.79s
```

추가로 기본 EIT/AT/CPT solve, 601–9,601점 resolution 진단, 601점
center-clustered scan, cold 0/10 mrad 비교, line-strength remap benchmark,
isotope별 종방향 residual-Doppler 폭을 실행했다.

## 10. 최종 판단과 우선순위

Scheme 2는 실제 3준위 Lindblad OBE, Maxwell velocity average, lab-facing Rabi
scaling, Beer–Lambert propagation을 결합한다. 따라서 **EIT/AT/CPT의 scale과
knob sensitivity를 빠르게 보는 반정량적 실험 reference**로는 충분히 유용하다.

그러나 현재 가장 먼저 고칠 것은 더 큰 atomic manifold가 아니라 다음 두 P1이다.

1. 같은 601 solve를 center-clustered grid에 배치하고 edge interpolation,
   samples-per-width, resolution status로 linewidth와 `n_g`의 수렴성을 회복한다.
2. scan optical-shift와 velocity residual operator를 분리해 Doppler-off angle
   invariance를 복원한다.

그 다음에는 종방향 hyperfine `Δk`, modeled contrast/experimental visibility,
natural-Rb scalar 경고, Raman-dephasing 명칭, group-index와 loss의 동시 표시가
비용 대비 효과가 크다. 순수 성능 면에서는 `line_strength`의 navigate-only 이동이
약 `689×`의 불필요한 재계산을 없애므로 가장 명확하다. full VCC와
hyperfine/Zeeman Λ는 가치가 있지만 별도 검증 자료와 opt-in heavy solver 설계가
필요한 장기 항목이다.
