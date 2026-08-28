# 2026-08-13 Scheme 4 물리 검토 — Hanle / EIA / NMOR

## 1. 오늘의 선택과 현재 scheme 순서

- 로컬 날짜: 2026-08-13 (Asia/Seoul)
- 계산: `n = (13 mod 5) + 1 = 4`
- 선택: **Scheme 4, `magneto` — Magneto-optics (Hanle/MOR)**

현재 드롭다운 순서는 실제 registry의 `_SCHEMES` 배열 기준으로 다음과 같다
(`gabes/schemes/__init__.py:19-25`). README의 사용자 설명도 같은 다섯 항목을 같은
순서로 제시한다 (`README.md:8-16`).

| 순서 | registry 이름 | 사용자-facing scheme | 핵심 출력 |
|---:|---|---|---|
| 1 | `sas` | OD / SAS | Doppler OD, hyperfine-pumping SAS |
| 2 | `lambda` | Lambda coherence | EIT / AT / CPT |
| 3 | `rydberg_eit` | Rydberg-EIT electrometry | cascade EIT, microwave AT, finite-IF electrometry |
| 4 | `magneto` | Magneto-optics | Hanle/EIA transmission, NMOR rotation |
| 5 | `fwm` | Four-wave mixing | seeded squeezing, SFWM biphoton |

## 2. 먼저 확인한 기존 문서·제안·TODO·issue

새 제안을 만들기 전에 다음을 검색하고 읽었다.

- 구현: `gabes/schemes/magneto.py:1-812`, `gabes/zeeman.py:1-162`,
  `gabes/atoms.py:22-109`, `gabes/kernels.py:276-410`,
  `gabes/observables.py:393-418`
- 사용자 문서: `README.md:8-16`,
  `docs/Userguide/GABES_User_Guide_v2.html:538-541`, `616-627`, `811-813`
- 테스트·검증 스크립트: `tests/test_magneto.py:28-535`,
  `tests/test_kernels.py:62-99`, `tests/test_schemes_render.py:118-126`,
  `tests/test_resonant_hanle_reference.py:14-131`,
  `tests/verify_hanle_eit_eia.py:120-206`, `240-271`
- 실험 예제: `analysis/squeezing/resonant_hanle_squeezing_reference.py:352-412`,
  `527-626`, `657-727`, `829-899`, `1006-1158` 및
  `analysis/squeezing/resonant_hanle_experiment_config.example.json:1-47`
- 생성된 실험 문서: `analysis/squeezing/resonant_hanle_squeezing_reference.md:18-60`,
  `docs/Resonant Hanle analysis/resonant_hanle_squeezing_reference_ko.tex:49-64`,
  `369-402`
- 계획: `docs/checklist.json:209-239`, `339-430`, `563-582`
- 과거 Scheme 4 보고서: 2026-06-23, 06-28, 07-03, 07-13, 07-23,
  07-28, 08-03, 08-08의 `docs/daily_report/*scheme-4_hanle-eia-nmor.md`

`docs/checklist.json:3`은 승격 검토가 2026-08-07 보고서까지만 동기화됐고 08-08
보고서는 cutoff 밖이라고 명시한다. 따라서 08-08 보고서의 NMOR 편광 범위 제안은 현재
체크리스트에 아직 승격되지 않았다. 코드 안의 관련 TODO는 Ne 계수표/pressure shift 승격
메모 하나다 (`gabes/constants.py:79-89`). 공개 GitHub Issues API도 현재 issue 0개를
반환했다.

중요하게도 08-08 보고서 뒤 현재 작업 트리에는 `magneto.py`, `zeeman.py`,
`tests/test_magneto.py`의 미커밋 수정이 있다. 아래 평가는 마지막 commit이 아니라 이 **현재
작업 트리**를 기준으로 한다.

## 3. 현재 구현이 담고 있는 실제 원자물리

이 scheme은 임의 Lorentzian을 그리는 toy model은 아니다.

1. 선택한 ⁸⁷Rb D1 `Fg -> Fe`의 모든 `mF` 상태를 만들고, `q=-1,0,+1`의
   Clebsch–Gordan 결합을 계산한다 (`gabes/zeeman.py:69-123`).
2. 자발방출을 편광별 `Sigma_q` jump operator로 묶어 excited-state Zeeman coherence가
   ground-state coherence로 전달되는 TOC를 보존한다 (`gabes/zeeman.py:101-123`).
   정규화와 TOC matrix element는 테스트가 직접 고정한다
   (`tests/test_magneto.py:28-64`).
3. QWP 각도를 복소 `sigma+/-` 구동으로 바꾸고 (`gabes/schemes/magneto.py:117-128`),
   longitudinal scan, B-zero offset, residual transverse field를 ground/excited Zeeman
   Hamiltonian에 넣는다 (`gabes/schemes/magneto.py:440-455`, `537-556`).
4. 파라핀 셀은 illuminated/dark 두 density-matrix block의 교환으로 Ramsey narrowing을
   나타내며 (`gabes/schemes/magneto.py:497-507`, `558-617`), buffer cell은 단일 영역
   steady state를 푼다 (`gabes/schemes/magneto.py:508-510`, `619-642`).
5. 현재 작업 트리는 자연방출과 Ne elastic optical dephasing을 분리한다. 자연 수명·branching·
   TOC는 `GAMMA_D1`에 남고, Ne FWHM의 절반만 optical-coherence dephasing으로 들어간다
   (`gabes/schemes/magneto.py:411-425`; `gabes/zeeman.py:132-156`). 이는 충돌 폭을
   자발방출률로 잘못 쓰던 이전 구현보다 물리적으로 타당하다.
6. 파라핀 light block은 이제 trace-one conditional state로 정규화되고, light-region
   occupation은 별도 진단값이다. local vapor density는 감수율에 정확히 한 번만 들어간다
   (`gabes/schemes/magneto.py:468-474`, `502-535`, `558-565`).
7. 흡수는 `alpha=k Im(chi)`, 투과는 Beer–Lambert, 회전은
   `theta=kL Re(chi_+-chi_-)/4`로 계산한다 (`gabes/schemes/magneto.py:683-705`).
8. 별도 분석 예제는 measured trace가 있으면
   `signal = offset + scale * model(B_scale*B_measured+B_offset)`를 fit하고
   (`analysis/squeezing/resonant_hanle_squeezing_reference.py:527-626`), 검출 shot/electronic
   noise를 field sensitivity로 환산한다
   (`analysis/squeezing/resonant_hanle_squeezing_reference.py:657-727`).

따라서 **부호 전환, transverse-field 의존성, Ramsey core, TOC 기반 EIA, NMOR 분산형
신호, measured-trace 초기 fit**을 연구하는 데 실제로 유용하다.

## 4. 실험 연구자 관점의 현재 판정

### 유용한 용도

- QWP, 잔류장, wall lifetime, buffer relaxation을 바꿨을 때 dip/peak 부호와 폭의 방향성을
  확인하는 실험 설계
- TOC 기반 intrinsic EIA와 circular-light LCA가 나타나는 조건의 sanity check
- shield/coil offset과 유효 relaxation의 초기 탐색
- measured Hanle trace의 B축 scale/offset 및 affine detector scale fit
- headless parameter sweep와 정성·반정량 trend 비교

### 아직 부적합한 용도

- 절대 linewidth, contrast, NMOR slope, `pT/sqrt(Hz)`를 독립 표준값처럼 인용
- 단일 Ne 압력으로 pressure shift, diffusion, spin destruction, Dicke narrowing을 예측
- 다른 ground hyperfine manifold로의 실제 누출과 repump가 중요한 open transition의
  절대 population dynamics
- 임의 타원/원편광에서 현재 `Re(chi_+-chi_-)`를 곧바로 편광면 회전이라고 해석

최종 평가는 **실제 원자물리를 담은 유용한 정성·반정량 실험 reference**이지만,
**absolute magnetometry/line-shape reference는 아니다**. 현재 테스트는 내부 불변량과 trend에
강하지만 held-out measured spectrum으로 절대 scale과 선폭을 닫지 않았다.

## 5. 08-08 이후 구현된 기존 P1의 평가

| 기존 항목 | 현재 상태 | 계산비용 | 물리 보존 평가 |
|---|---|---:|---|
| `magneto-dynamics-semantics` (`docs/checklist.json:339-357`) | `done` | 같은 Hilbert/Liouville 차원과 같은 `(B,v)` solve. dissipator 행렬만 달라져 runtime overhead는 사실상 0 | 자연방출/TOC와 elastic broadening 분리는 올바른 수정이다. rate-swap, trace, Hermiticity, positivity 테스트도 추가됐다 (`tests/test_magneto.py:66-159`, `352-446`). 다만 아래 7.2의 `collisional_depol` 의미는 아직 완전하지 않다. |
| `magneto-light-region-normalization` (`docs/checklist.json:360-379`) | `done` | solve 후 trace 정규화 `O(n_B n_v n_level^2)`, 새 solve 없음 | local density convention을 일관되게 만든 올바른 수정이다. occupation=1/81과 trace-one, density-once, gamma_out=0 limit이 고정돼 있다 (`tests/test_magneto.py:208-302`). |

정규화 수정은 absolute scale을 크게 바꾸므로 단순 표시 변경이 아니다. 현재 기본 121×9
EIT에서 `T(B=0)=0.8684`이며, 수정 전 08-08 보고서의 `0.9983`과 직접 비교할 수 없다.

## 6. 현재 수치와 수렴성

### 6.1 built-in 121×9 기본값

| regime | 중앙 absorption feature [m⁻¹] | 중앙 폭 [µT] | `T(B_cmd=0)` | NMOR slope [mrad/µT] |
|---|---:|---:|---:|---:|
| EIT dip | -1.58463 | 0.16626 | 0.86841 | -81.41 |
| EIA peak | +5.87981 | 0.22327 | 0.94260 | -15.76 |
| Buffer Hanle | -6.87430 | 4.59180 | 0.78051 | -5.65 |
| Buffer LCA | +0.01455 | 1.77694 | 0.66276 | -0.031 |
| NMOR | -4.11793 | 0.13202 | 0.90948 | -129.42 |

여기서 feature는 absorption 중앙값과 양 끝 평균의 차이다. NMOR 행의 폭은 동일 raw state의
absorption 폭이고, slope가 실제 NMOR 핵심 진단값이다.

### 6.2 고해상도 reference schedule

| regime | `(B,v)` | feature [m⁻¹] | 중앙 폭 [µT] | slope [mrad/µT] | runtime | 반환 `rho` |
|---|---:|---:|---:|---:|---:|---:|
| EIT | 401×321 | +0.31909 | 1.02995 | -114.25 | 11.51 s | 125.7 MiB |
| EIT | 401×641 | +0.35718 | 0.98055 | -116.79 | 23.77 s | 251.0 MiB |
| NMOR | 401×321 | -2.77992 | 0.11034 | -307.82 | 13.99 s | 125.7 MiB |
| NMOR | 401×641 | -2.81604 | 0.10990 | -314.55 | 27.89 s | 251.0 MiB |

- 기본 EIT는 641-class에서 **dip에서 peak로 부호가 뒤집힌다**.
- 321→641에서 EIT 폭 변화는 약 4.8%, NMOR slope 변화는 약 2.2%이지만 EIT amplitude는
  약 11.9% 변한다.
- 기본 NMOR slope 크기는 641-class 값보다 약 59% 작다.
- 따라서 `magneto-observable-convergence-and-trust`의 기본 `coarse/nonconverged` 상태 표시는
  여전히 필요하다 (`docs/checklist.json:382-400`). 고해상도 solve를 UI 기본값으로 올리는 것은
  비용상 부적합하고, opt-in reference 또는 느린 검증 경로가 맞다.

## 7. 남은 기존 제안과 이번에 확인한 새 문제

### 7.1 기존 `magneto-observable-convergence-and-trust` — P1, 그대로 유효

대부분의 신뢰성 수정은 새 solve 없이 가능하다.

- `b_offset_ut=0.25 µT` 진단에서 command zero의 feature는 `-0.15793 m⁻¹`이지만
  physical zero 샘플은 `-1.44453 m⁻¹`이다. 현재 view는 command zero를 기준으로
  `crossover`, FWHM `2.355 µT`를 hero로 낸다. Hamiltonian에는 physical axis가 올바르게
  들어가지만 (`gabes/schemes/magneto.py:440-455`), observables는 여전히 `raw["b_ut"]`의
  중앙을 쓴다 (`gabes/schemes/magneto.py:683-705`). physical/command zero 분리는
  `O(n_B)` postprocess다.
- NMOR slope는 아직 중앙의 단일 `np.gradient` 값이다
  (`gabes/schemes/magneto.py:712-727`). local odd fit 또는 고차 중앙 미분은 `O(n_B)`다.
- `tests/verify_hanle_eit_eia.py`의 현재 출력은 paraffin `1.193/1.737 mG` 대 문헌
  `0.12/0.20 mG`, 즉 약 `9.94×/8.69×` 차이인데도 `MATCH`, `same sub-mG order`라고
  쓴다 (`tests/verify_hanle_eit_eia.py:243-271`). Buffer 저전력점도 현재 `3.389 mG`
  대 `2.4 mG`인데 “reaches”라고 출력한다 (`tests/verify_hanle_eit_eia.py:174-206`).
  trend PASS와 absolute-width CHECK/FAIL을 분리하는 것은 계산비용 0이다.
- 실제 ⁸⁷Rb D1 open transition은 다른 ground hyperfine manifold로 decay하지만 현재 builder는
  선택한 `Fg` 안으로만 자발방출을 되돌린다 (`gabes/zeeman.py:94-130`). 먼저 model-scope
  status를 `O(1)`로 표시해야 한다. full open-hyperfine/repump solve는 상태 수와 검증 표면을
  키우므로 opt-in reference로 분리해야 한다.

### 7.2 새 P1: 파라핀 `transit_relax_khz`의 이중 적용

UI는 이 값을 “atoms leaving the illuminated region”으로 정의한다
(`gabes/schemes/magneto.py:267-272`). 생산 코드에서는 같은 `gamma_light`가

1. `gamma_out` light→dark 교환률이 되고 (`gabes/schemes/magneto.py:429-432`),
2. light atom의 `gamma_gg` projective dephasing에도 다시 들어간다
   (`gabes/schemes/magneto.py:433-438`).

교환 행렬 자체가 light-block coherence를 `-gamma_out`으로 제거하고 dark block으로 전달한다
(`gabes/schemes/magneto.py:599-607`; `gabes/kernels.py:354-363`). 따라서 하나의 “beam을
떠나는 rate”가 교환과 별도 순수 탈위상으로 두 번 작용한다. 201×41 기본 EIT 진단에서 별도
light `gamma_gg`만 제거하고 교환은 유지하면

- 최대 `chi_probe` 상대 차이: **48.8%**
- feature: `-1.44274 -> -4.84764 m⁻¹`
- 중앙 폭: `0.15206 -> 0.14481 µT`

로 바뀌었다. 이는 무시할 수 없는 의미 중복이다. 권고는 **transit은 두-region exchange로만
표현**하고, 실제 in-beam pure dephasing이 필요하면 별도 이름·근거·rate로 추가하는 것이다.
행렬 크기나 solve 횟수는 변하지 않으므로 runtime overhead는 0이다. production construction을
직접 고정하는 회귀 테스트가 필요하다. 현재 gamma_out→0 fixture는 별도 `gamma_gg`가 없는
2-level fixture라 이 중복을 잡지 못한다 (`tests/test_magneto.py:223-258`).

### 7.3 새/재개 P1: `collisional_depol_khz`는 population depolarization이 아니다

UI help는 이를 spin-destruction/depolarization으로 설명한다
(`gabes/schemes/magneto.py:261-266`). 그러나 구현은 ground projector들의 pure-dephasing
channel이라 모든 ground diagonal population을 그대로 둔다 (`gabes/zeeman.py:132-149`).
테스트도 population column 변화가 정확히 0임을 요구한다 (`tests/test_magneto.py:105-159`).
완전 편극된 `mF=+F` population에 이 dissipator만 적용한 진단에서도 모든 population derivative가
0이었다. 즉 orientation/alignment의 population 성분을 파괴하는 spin-destruction 모델은 아니다.

두 가지 안전한 선택이 있다.

- 의도가 coherence-only라면 이름을 `ground_pure_dephasing_khz`로 바꾸고 범위를 명시한다.
- 실제 depolarization이면 ground block을 trace-preserving isotropic state로 완화하는 완전양의
  channel을 구현하고 orientation/alignment decay fixture를 추가한다.

둘 다 같은 Liouville 차원에서 dissipator만 바꾸므로 계산 overhead는 사실상 0이다.

### 7.4 08-08 기존 제안: 임의 QWP에서 NMOR plane-rotation hero를 제한

현재 QWP=45°는 `|E_+|=0`, `|E_-|=sqrt(2)`의 순수 원편광이다
(`gabes/schemes/magneto.py:117-128`). 그런데 `_coherences()`는 각 circular susceptibility를
해당 입력 성분이 아니라 공통 total `Omega`로 나누고 (`gabes/schemes/magneto.py:644-659`),
모든 QWP에서 같은 plane-rotation 식을 적용한다 (`gabes/schemes/magneto.py:683-695`).
현재 QWP=45° NMOR view는 여전히

- `Rotation at B=0 = 0.74 mrad`
- `Slope = -5.15 mrad/µT`
- `Peak |rotation| = 0.74 mrad`
- “zero crossing near B=0” note

를 낸다. 순수 원편광에는 선형 편광면의 major-axis angle이 정의되지 않으므로 이 hero는 물리적
범위를 벗어난다. 먼저 비선형 입력에서 `unsupported/diagnostic` gate를 거는 것은 `O(1)`이다.
그 뒤 thin/undepleted 조건에서 입력 성분별 정규화와 Jones/Stokes postprocess를 추가하면
ellipse orientation과 ellipticity를 `O(n_B)`로 계산할 수 있다. self-consistent polarization
propagation은 별도의 고비용 연구 문제다.

### 7.5 새 P1 reference 문제: 생성된 FWM–Hanle 산출물이 현재 물리와 불일치

생성된 Markdown은 probe-resonant 행에 `T=0.9962`,
`|dT/dB|=2.533e-3 /µT`, `29 pT/sqrt(Hz)`를 기록한다
(`analysis/squeezing/resonant_hanle_squeezing_reference.md:20-25`). 한글 TeX/PDF도 같은
`29 pT/sqrt(Hz)`와 `5.25 pT/sqrt(Hz)`를 요약한다
(`docs/Resonant Hanle analysis/resonant_hanle_squeezing_reference_ko.tex:49-64`,
`386-397`). 그러나 현재 코드로 `compute_hanle(detuning=0)`만 재실행하면 같은 bias
`-0.04 µT`에서

- `T = 0.73733`
- `|dT/dB| = 0.14899 /µT`

이다 (`analysis/squeezing/resonant_hanle_squeezing_reference.py:396-412`). 정규화 수정 뒤
slope가 약 58.8배 달라졌으므로 기존 절대 sensitivity 표를 현재 코드의 reference로 사용하면 안 된다.

권고는 physics 변경 뒤 산출물을 재생성하고, 보고서에 생성 시각, git commit, dirty 상태,
`MagnetoScheme.cache_version`, 주요 입력과 source hash를 기록하는 것이다. 메타데이터 비용은
무시 가능하다. measured CSV가 없다는 현재 문서의 명시
(`analysis/squeezing/resonant_hanle_squeezing_reference.md:45-60`)도 유지해야 한다.

## 8. 나머지 기존 개선안의 비용과 우선순위

| 항목 | 상태 | 계산비용 | 권고 |
|---|---|---:|---|
| `magneto-observable-convergence-and-trust` (`docs/checklist.json:382-400`) | P1 ready | status/physical-zero/local fit/literature label은 `O(n_B)` 또는 0. 641-class solve만 24–28 s, 251 MiB | 싼 trust/readout 수정은 기본값에 즉시 적용하고, 고해상도 schedule은 opt-in/slow test로 둔다. |
| `collisional-coefficient-provenance-and-pressure-shift` (`docs/checklist.json:209-239`) | P1 ready, large | table lookup, scalar detuning shift, 기존 dephasing matrix이므로 runtime overhead 거의 0 | gas/species/line/온도 범위·단위·부호·불확도를 소싱한 뒤 적용한다. 구현 작업은 크지만 solve 비용은 늘지 않는다. |
| `geometry-aware-buffer-relaxation-budget` (`docs/checklist.json:404-429`) | P2 needs_decision, large | geometry에서 scalar rate를 계산하므로 solve overhead 거의 0 | geometry, wall condition, diffusion/spin-destruction 자료를 먼저 선택한다. pressure만으로 자동 rate를 만들면 false precision이다. |
| `full-velocity-changing-collision-kernel` (`docs/checklist.json:563-582`) | P2 parked, research | velocity classes를 결합해 현재 separable batch solve를 깨므로 메모리·시간 증가가 큼 | 특정 kernel·dataset·runtime budget이 정해질 때만 opt-in reference로 구현한다. interactive 기본 경로에는 부적합하다. |

## 9. 동작을 바꾸지 않는 순수 코드 최적화

### 9.1 가장 우선: kernel-side weighted coherence contraction

현재 kernel은 `(B,v,n,n)` 전체 `rho`를 반환하고
(`gabes/kernels.py:310-410`), Python이 필요한 세 coherence를 뽑아 Maxwell 가중합한다
(`gabes/schemes/magneto.py:502-515`, `644-659`). 401×641, 8-level에서 반환 배열만
251.0 MiB다.

kernel 내부 velocity loop에서 `chi_+`, `chi_-`, `chi_probe`의 가중합만 누적해 반환하면
출력은 세 개의 401-point complex array, 약 19 KiB 수준이다. LU solve와 관측량 정의는 그대로이므로
물리는 변하지 않는다. 현재 NumPy fallback 및 kernel parity 기준
(`tests/test_kernels.py:82-99`)으로 세 susceptibility의 상대오차를 고정해야 한다.

### 9.2 affine real-basis assembly를 먼저 수행

현재는 complex `C_xy + B*C_z` stack을 만든 뒤 kernel wrapper가 각 B stack에
`U†LU`를 적용한다 (`gabes/schemes/magneto.py:485-501`;
`gabes/kernels.py:321-330`, `399-410`). `C_xy`와 `C_z`만 먼저 real Hermitian basis로
바꾸고 그 basis에서 B-affine stack을 만들 수 있다.

401×64×64 light stack microbenchmark:

- 현재 조립+변환: **68.81 ms**
- 선변환 affine 조립: **8.20 ms**
- 속도 향상: **8.39×**
- generator 최대 절대 차이: `5.82e-11` (행렬 scale 대비 roundoff)
- complex 중간 stack: 25.06 MiB; real stack: 12.53 MiB

두-region은 light/dark 두 stack 모두 이 이득을 받는다. 전체 runtime은 LU solve가 지배하므로
8.39×가 전체 solve 속도 향상을 뜻하지는 않지만, allocation과 basis-transform 비용을 안전하게
줄인다.

### 9.3 작은 최적화: angular-momentum matrix cache

`_hamiltonian()`이 호출될 때 `Fx,Fy,Fz`를 다시 만든다
(`gabes/schemes/magneto.py:537-545`; `gabes/zeeman.py:48-66`). `(Fg,Fe)`별 작은 LRU cache는
동작을 바꾸지 않지만, 호출 횟수가 계산당 몇 번뿐이므로 위 두 항목보다 우선순위가 낮다.

`zeeman_manifold()`와 Hermitian basis는 이미 cache돼 있으므로 중복 제안하지 않는다
(`gabes/zeeman.py:69-70`; `gabes/core.py:60-93`). Headless observables도 이미 구현됐다
(`gabes/schemes/magneto.py:136`, `661-812`).

## 10. 권고 실행 순서

1. **P1 semantics:** 파라핀 transit exchange와 light pure-dephasing 중복을 제거하고,
   `collisional_depol`을 실제 depolarization 또는 명시적 pure dephasing으로 확정한다.
2. **P1 trust:** physical zero/command zero 분리, coarse/nonconverged 및 open-hyperfine scope status,
   NMOR local slope fit, 문헌 PASS/CHECK/FAIL 라벨을 구현한다.
3. **P1 polarization:** 선형 편광면이 정의되지 않는 QWP 범위에서 NMOR hero를 gate하고,
   검증된 Jones/Stokes 후처리를 추가한다.
4. **P1 reproducibility:** 현재 물리로 FWM–Hanle Markdown/NPZ/PNG/TeX/PDF를 재생성하고
   모델 버전·git 상태·source hash를 기록한다. 재생성 전 기존 pT 수치는 stale로 표시한다.
5. **P1/P2 coefficients:** sourced pressure shift/elastic coefficient table을 추가하고,
   geometry-aware relaxation은 입력·자료가 정해진 뒤 진행한다.
6. **성능:** kernel-side coherence contraction을 먼저, real-basis affine assembly를 다음으로 한다.
7. **Parked:** full VCC와 full open-hyperfine reference는 dataset과 runtime budget 없이 기본 경로로
   넣지 않는다.

## 11. 검증 기록

- registry 직접 확인: `sas -> lambda -> rydberg_eit -> magneto -> fwm`
- 공개 GitHub Issues API: issue 0개
- targeted:
  `python -m pytest -q tests/test_magneto.py tests/test_kernels.py tests/test_headless_observables.py tests/test_schemes_render.py tests/test_resonant_hanle_reference.py`
  -> **72 passed in 16.09 s**
- `python tests/verify_hanle_eit_eia.py`:
  paraffin `1.193/1.737 mG` 대 문헌 `0.12/0.20 mG`; buffer 저전력점 `3.389 mG`
  대 문헌 `2.4 mG`; 스크립트의 정성 label은 여전히 과장됨
- full suite: `python -m pytest -q` -> **402 passed in 64.73 s**
- 기존 dirty README/checklist/Rydberg/SAS/SABES/분석 작업은 보존했으며, 이번 자동화는 이
  보고서와 자동화 memory만 변경한다.
