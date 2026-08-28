# 2026-08-23 Scheme 4 물리 검토 — Hanle / EIA / NMOR

## 1. 오늘의 scheme 선택과 현재 순서

- 로컬 날짜: **2026-08-23** (`Asia/Seoul`)
- 계산: `n = (23 mod 5) + 1 = 4`
- 오늘의 대상: **Scheme 4 — Magneto-optics (Hanle / EIA / NMOR)**

현재 dropdown registry는 다음 순서다.

| n | registry key | 현재 정의 |
|---:|---|---|
| 1 | `sas` | Absorption OD/SAS: pump off의 Doppler OD와 pump on의 hyperfine-pumping SAS |
| 2 | `lambda` | 3-level Λ coherence: EIT / AT / CPT |
| 3 | `rydberg_eit` | 5S–5P–40D cascade EIT와 40D–39F microwave AT electrometry |
| 4 | `magneto` | ⁸⁷Rb D1 Zeeman manifold의 Hanle / EIA / NMOR |
| 5 | `fwm` | seeded ⁸⁵Rb D1 double-Λ gain 진단과 generic SFWM biphoton |

근거는 `gabes/schemes/__init__.py:12-25`의 병합 scheme 설명과 `_SCHEMES` 순서이며,
사용자 문서의 같은 정의는 `README.md:8-16`에 있다. 직접 실행한 registry도
`sas -> lambda -> rydberg_eit -> magneto -> fwm`이었다.

## 2. 먼저 검색한 기존 제안·문서·TODO·issue

검토 전에 다음을 검색했다.

- 전체 Scheme 4 일일 보고서 9개:
  `docs/daily_report/2026-06-23_scheme-4_hanle-eia-nmor.md`부터
  `docs/daily_report/2026-08-13_scheme-4_hanle-eia-nmor.md`까지
- 통합 작업 목록: `docs/checklist.json`
- 소스의 `TODO/FIXME/issue/proposal` 표식
- `README.md`, `docs/Userguide/GABES_User_Guide_v2.html`
- Scheme 4 코드와 공용 엔진: `gabes/schemes/magneto.py`, `gabes/zeeman.py`,
  `gabes/kernels.py`, `gabes/core.py`, `gabes/observables.py`, `gabes/constants.py`
- 테스트: `tests/test_magneto.py`, `tests/test_kernels.py`,
  `tests/test_headless_observables.py`, `tests/test_schemes_render.py`,
  `tests/test_resonant_hanle_reference.py`, `tests/verify_hanle_eit_eia.py`
- 예제/분석:
  `analysis/squeezing/resonant_hanle_squeezing_reference.py`,
  `analysis/squeezing/resonant_hanle_experiment_config.example.json`, 생성 Markdown과
  `docs/Resonant Hanle analysis/` 산출물
- 연결된 공개 GitHub 저장소 `Shake2313/fwm-squeezing-app`의 전체 issue: **0개**

별도 로컬 issue/TODO 파일은 없다. 개선안은 과거 일일 보고서와 `docs/checklist.json`에
집약돼 있다. Scheme 4 관련 핵심 registry 항목은 다음과 같다.

| 기존 항목 | 상태 | 위치 |
|---|---|---|
| lifetime / elastic dephasing / transit / depolarization 의미 분리 | `done`, P1 | `docs/checklist.json:339-357` |
| two-region light-state 정규화와 density-once | `done`, P1 | `docs/checklist.json:360-379` |
| observable 수렴·physical-zero·문헌 비교 신뢰 경계 | `ready`, P1 | `docs/checklist.json:382-401` |
| collision coefficient provenance와 pressure shift | `ready`, P1 | `docs/checklist.json:209-240` |
| geometry-aware buffer relaxation budget | `needs_decision`, P2 | `docs/checklist.json:404-429` |
| full velocity-changing-collision reference | `parked`, P2 | `docs/checklist.json:563-582` |

8월 13일 보고서에 새로 기록된 파라핀 transit 중복, population depolarization 명칭,
원편광 NMOR 범위, stale FWM–Hanle 산출물은 아직 독립 checklist item으로 승격되지 않았다
(`docs/daily_report/2026-08-13_scheme-4_hanle-eia-nmor.md:177-256`). 따라서 “기존 개선
제안이 없다”가 아니라, **구체적이고 재현 가능한 제안이 이미 여러 개 있으며 일부만 통합 목록에
반영된 상태**다.

Scheme 4 관련 마지막 commit은 문서 이동을 포함한 `a82bbf7`(2026-07-23)이고, 8월 13일 이후
관련 경로의 새 commit은 없다. 다만 현재 작업 트리에는 사용자 소유의 큰 미커밋 변경이 있으며,
오늘 수치는 그 현행 working tree를 평가한 것이다.

## 3. 현재 구현하는 물리

### 3.1 실제로 구현된 부분

1. ⁸⁷Rb D1의 선택한 `F_g ↔ F'_e` Zeeman sublevel을 명시적으로 만든다. 허용 전이와
   level 수를 검사하고 hyperfine Landé `g_F`를 계산한다
   (`gabes/schemes/magneto.py:394-410`).
2. `zeeman_manifold()`는 `m_F` 상태, σ±/π Clebsch–Gordan 결합, 편광별 grouped
   spontaneous-emission jump operator `Σ_q`를 구성한다
   (`gabes/zeeman.py:70-130`). 같은 편광로 방출되는 excited coherence가 ground coherence로
   전달되는 TOC(transfer of coherence)가 있어 `F_e=F_g+1` intrinsic EIA를 만들 수 있다.
3. Ne optical broadening은 자연방출 수명·branching·TOC와 분리된 완전양의 elastic optical
   dephasing channel로 들어간다 (`gabes/schemes/magneto.py:411-425`;
   `gabes/zeeman.py:132-156`). 이는 Ne broadened FWHM을 자발방출률로 쓰던 이전 구현보다
   물리적으로 타당하다.
4. magnetic Hamiltonian은 ground/excited Zeeman 항, longitudinal scan, residual transverse
   field, QWP가 만든 σ± drive를 포함한다 (`gabes/schemes/magneto.py:440-460`, `538-556`).
5. 파라핀 셀은 light/dark 두 영역의 교환 OBE와 wall coherence를 풀고, buffer 셀은
   single-region OBE와 isotropic reservoir reload를 쓴다
   (`gabes/schemes/magneto.py:468-510`, `568-642`).
6. 파라핀 light block은 trace-one conditional state로 정규화하고 light-region occupation을
   별도 진단값으로 둔다. vapor density `N`은 감수율에 정확히 한 번 들어간다
   (`gabes/schemes/magneto.py:502-534`, `558-565`).
7. Maxwell velocity class를 평균하고, 작은 coarse grid의 optical opacity는 fine scalar Voigt
   average로 보정한다 (`gabes/schemes/magneto.py:72-104`, `476-515`).
8. 관측량은 `α = k Imχ`, `T = exp(-αL)`,
   `θ = kL Re(χ_+ - χ_-)/4`이다 (`gabes/schemes/magneto.py:683-705`).
9. 별도 FWM–Hanle 예제는 measured CSV가 있으면
   `signal = offset + scale*model(B_scale*B_measured+B_offset)`를 fit하고
   (`analysis/squeezing/resonant_hanle_squeezing_reference.py:469-626`), shot/dark/electronic
   current noise를 field sensitivity로 환산한다
   (`analysis/squeezing/resonant_hanle_squeezing_reference.py:657-727`).

이 구조는 단순 Lorentzian 장난감이 아니다. CG selection rule, optical pumping, TOC, Zeeman
precession, transverse-field-induced sign change, wall-preserved Ramsey core, Beer–Lambert readout을
실제 밀도행렬 계산으로 연결한다. 문헌에서도 파라핀 ⁸⁷Rb D1에서 linear CPT가 circular
magnetic-field-induced absorption으로 바뀌며 0.12/0.20 mG 폭을 보였고
([Lee & Moon, JOSA B 30, 2301 (2013)](https://doi.org/10.1364/JOSAB.30.002301)),
Ne buffer cell에서는 circular polarization과 transverse field가 관련된 2.4 mG LCA가 보고됐다
([Yu et al., PRA 81, 023416 (2010)](https://doi.org/10.1103/PhysRevA.81.023416)).
EIA에서 excited-to-ground coherence transfer의 역할도 독립 실험으로 검증된 메커니즘이다
([Failache et al., arXiv:quant-ph/0211065](https://arxiv.org/abs/quant-ph/0211065)).

### 3.2 연구 reference로서의 범위 제한

- 선택한 단일 `F_g ↔ F'_e` manifold 안에서만 decay가 닫힌다. 실제 D1 open transition에서
  다른 ground hyperfine manifold로의 누출과 repump가 없다 (`gabes/zeeman.py:94-130`).
- `NEON_BUFFER_BROADENING_MHZ_PER_TORR = 3.91` 하나만 있고 gas/species/line/temperature
  provenance, pressure shift, diffusion, spin destruction, Dicke narrowing이 없다
  (`gabes/constants.py:80-89`).
- 기본 Doppler quadrature는 9 classes이고 UI 최대도 41 classes다
  (`gabes/schemes/magneto.py:295-302`). fine scalar opacity 보정은 전체 B-dependent velocity
  dynamics의 수렴을 대신하지 못한다.
- `line_strength`는 실험변수가 아닌 effective-OD calibration knob다
  (`gabes/schemes/magneto.py:290-294`). measured absolute transmission 없이 조절하면 절대
  contrast의 예측력을 만들지 못한다.
- NMOR의 `χ_+ - χ_-` 식은 thin/undepleted linear-polarization readout에는 유용하지만, 임의
  ellipticity에서 출력 ellipse를 Jones/Stokes로 전파하는 모델은 아니다. 정밀 NMOR는 atomic
  alignment의 생성·검출과 편광 전파 범위를 명시해야 한다
  ([Budker et al., RMP 74, 1153 (2002)](https://doi.org/10.1103/RevModPhys.74.1153)).
- absolute contrast/rotation/slope는 코드도 `external validation required`로 표시한다
  (`gabes/schemes/magneto.py:41-43`, `720-752`).

따라서 이 scheme은 **실제 물리를 포함하는 유용한 정성·반정량 실험 설계/초기 fit reference**다.
QWP, residual field, cell relaxation, detuning, intensity 변화의 방향과 sign topology를 탐색하는
데 적합하다. 그러나 **절대 linewidth, contrast, NMOR slope 또는 pT/√Hz를 독립 표준값처럼
인용하는 experimental-physics reference는 아니다**.

## 4. 현재 수치와 수렴성

### 4.1 built-in 121×9 결과

아래 feature는 중앙 absorption과 양쪽 scan wing 평균의 차이다. 음수는 EIT/CPT-like
absorption dip, 양수는 EIA/MIA-like absorption peak다.

| regime 순서 | feature [m⁻¹] | 중앙 폭 [µT] | `T(B_cmd=0)` | NMOR slope [mrad/µT] |
|---|---:|---:|---:|---:|
| EIT dip | -1.58463 | 0.16626 | 0.86841 | -81.41 |
| EIA peak | +5.87981 | 0.22327 | 0.94260 | -15.76 |
| Buffer Hanle | -6.87430 | 4.59180 | 0.78051 | -5.65 |
| Buffer LCA | +0.01455 | 1.77694 | 0.66276 | -0.031 |
| NMOR | -4.11793 | 0.13202 | 0.90948 | -129.42 |

이는 8월 13일 보고서 수치와 일치한다. 기본 첫 실행은 Numba warm-up을 포함해 약 0.63 s,
이후 각 regime은 약 0.09–0.14 s였다.

### 4.2 401×321→401×641 reference schedule

| regime | `(B,v)` | feature [m⁻¹] | 중앙 폭 [µT] | slope [mrad/µT] | runtime | 반환 `rho` |
|---|---:|---:|---:|---:|---:|---:|
| EIT | 401×321 | +0.31909 | 1.02995 | -114.25 | 11.14 s | 125.7 MiB |
| EIT | 401×641 | +0.35718 | 0.98055 | -116.79 | 21.85 s | 251.0 MiB |
| NMOR | 401×321 | -2.77992 | 0.11034 | -307.82 | 11.24 s | 125.7 MiB |
| NMOR | 401×641 | -2.81604 | 0.10990 | -314.55 | 22.93 s | 251.0 MiB |

- 기본 EIT는 641-class reference에서 **dip에서 peak로 부호가 바뀐다**.
- 321→641에서 EIT 폭은 약 4.8%, NMOR slope는 약 2.2% 변하지만 EIT amplitude는 약
  11.9% 변한다.
- 기본 NMOR slope `-129.42`는 641-class `-314.55 mrad/µT`보다 절댓값이 약 59% 작다.
- 그러므로 고해상도 solve를 UI 기본값으로 강제하는 것은 비용상 나쁘지만, 기본 hero에
  `coarse/nonconverged` status를 붙이는 것은 필수다. 이것이
  `magneto-observable-convergence-and-trust`의 핵심 acceptance다
  (`docs/checklist.json:395-400`).

## 5. 기존 개선안의 재평가와 계산비용

| 항목 | 현행 판정 | 추가 계산비용 | 물리 보존/권고 |
|---|---|---:|---|
| natural lifetime / elastic dephasing 분리 | 구현 완료 | 같은 `(B,v)` solve, 사실상 0 | 올바른 수정. emission branching/TOC가 Ne 압력에 오염되지 않는다. `tests/test_magneto.py:66-102`가 고정한다. |
| buffer reload와 ground dephasing 분리 | 구현 완료, 이름은 불완전 | dissipator만 달라져 0 | rate swap이 서로 다른 `χ`를 만들고 positivity fixture가 있다 (`tests/test_magneto.py:105-159`, `352-446`). 다만 현재 `collisional_depol`은 population orientation을 isotropic하게 만드는 channel이 아니라 projective pure dephasing이다. 아래 5.2 참조. |
| conditional light-state / density-once | 구현 완료 | solve 후 정규화 `O(n_B n_v n²)`, 실질적으로 미미 | trace-one, occupation 분리, density 2배→χ 2배가 테스트된다 (`tests/test_magneto.py:208-302`). |
| observable convergence/trust | 아직 `ready` | status·physical-zero·local fit은 `O(n_B)`; 641 solve만 22–23 s/251 MiB | 싼 신뢰 경계는 기본에, 고해상도는 opt-in reference에 둔다. |
| collisional coefficient/pressure shift | 아직 `ready` | table lookup·scalar shift·기존 dephasing matrix이므로 runtime overhead 거의 0 | 구현 작업은 large지만 물리 계산비용은 미미하다. source/units/sign/uncertainty와 unsupported 조합을 먼저 정의한다 (`docs/checklist.json:209-240`). |
| geometry-aware relaxation | `needs_decision` | geometry→scalar rates이므로 solve overhead 거의 0 | beam radius, cell geometry, wall, diffusion/spin-destruction source가 없으면 pressure-only false precision이다 (`docs/checklist.json:404-429`). |
| full VCC | `parked` | velocity classes를 결합해 시간·메모리 크게 증가 | 특정 kernel, dataset, runtime budget이 정해질 때만 opt-in reference로 활성화한다 (`docs/checklist.json:563-582`). |

### 5.1 P1: physical zero와 command zero를 분리

`b_offset_ut=0.25 µT`, 401×41에서

- command-zero feature: `-0.15793 m⁻¹`
- physical-zero feature: `-1.44453 m⁻¹`

Hamiltonian은 `b_physical = b_command + b_offset`을 올바르게 사용한다
(`gabes/schemes/magneto.py:440-444`). 그러나 observables는 여전히 `raw["b_ut"]`에서 가장
가까운 command-zero index를 중앙으로 잡는다 (`gabes/schemes/magneto.py:683-705`).
physical-zero interpolation과 두 기준의 동시 표시는 `O(n_B)` postprocess이며 새 OBE solve가 없다.

### 5.2 P1: 파라핀 transit 중복과 `collisional_depol` 명칭

파라핀의 하나뿐인 `transit_relax_khz`가

1. light→dark 교환 `gamma_out`으로 들어가고 (`gabes/schemes/magneto.py:429-432`),
2. light atom의 `gamma_gg` projective dephasing에도 다시 들어간다
   (`gabes/schemes/magneto.py:433-438`).

교환 행렬 자체가 light block을 `-gamma_out`으로 소거한다
(`gabes/schemes/magneto.py:599-607`; `gabes/kernels.py:377-386`). 별도 light pure dephasing만
제거하고 exchange를 유지한 201×41 진단은 다음과 같다.

| | 현행 | exchange-only |
|---|---:|---:|
| feature [m⁻¹] | -1.44274 | -4.84764 |
| 중앙 폭 [µT] | 0.15206 | 0.14481 |
| `T(B=0)` | 0.86326 | 0.93403 |
| NMOR slope [mrad/µT] | -88.19 | -230.72 |

`chi_probe` 최대 상대차는 **48.8%**다. 하나의 “beam을 떠나는 rate”를 교환과 순수
dephasing에 두 번 쓰지 말고 transit은 two-region exchange가 전담해야 한다. 필요한 in-beam
dephasing은 별도 이름·근거·rate로 추가해야 한다. Liouville 차원과 solve 수가 같으므로 runtime
overhead는 0이다.

또한 `collisional_depol_khz`는 UI에서 spin-destruction/depolarization으로 불리지만
(`gabes/schemes/magneto.py:261-266`), 구현은 모든 ground population을 보존하고 coherence만
감쇠한다 (`gabes/zeeman.py:132-149`; `tests/test_magneto.py:142-152`). 의도가 현재 동작이면
`ground_pure_dephasing_khz`로 이름을 바꾸는 것이 가장 안전하고 계산비용은 0이다. 실제
spin-destruction이면 ground block을 trace-preserving isotropic state로 완화하는 CP channel과
orientation/alignment decay test를 추가해야 하며, 이 역시 같은 Liouville 차원이라 solve 비용은
늘지 않는다.

### 5.3 P1: NMOR readout의 편광 범위

`_qwp_drive_weights(45°)`는 한 circular component가 정확히 0인 순수 원편광을 만든다
(`gabes/schemes/magneto.py:117-128`). 그런데 `_coherences()`는 각 circular coherence를 입력
성분이 아니라 공통 `Ω`로 나누고 (`gabes/schemes/magneto.py:645-659`), 모든 QWP에서 같은
plane-rotation 식을 적용한다 (`gabes/schemes/magneto.py:683-695`). NMOR preset에서 QWP만
45°로 바꾸면 여전히 `rotation(B=0)=0.74 mrad`, slope `-5.15 mrad/µT`와 “zero crossing” note를
낸다.

순수 원편광에는 선형 편광면의 major-axis angle이 정의되지 않는다. 우선 해당 범위를
`unsupported/diagnostic`으로 gate하는 것은 `O(1)`이다. 그다음 입력 σ± amplitude별 정규화와
thin-medium Jones/Stokes postprocess를 `O(n_B)`로 추가할 수 있다. self-consistent polarization
propagation은 별도의 고비용 연구 문제로 남겨야 한다.

### 5.4 P1: 문헌 비교 라벨과 생성 산출물

`python tests/verify_hanle_eit_eia.py`의 현행 출력은

- paraffin linear/circular: `1.193/1.737 mG` 대 논문 `0.12/0.20 mG`, 약 `9.94×/8.69×`
- buffer low-power LCA: `3.389 mG` 대 논문 `2.4 mG`, 약 `1.41×`

인데도 `MATCH`, “same sub-mG order”, “reaches”라고 쓴다
(`tests/verify_hanle_eit_eia.py:243-271`, `263-271`). 첫 두 계산값은 sub-mG조차 아니다.
sign/trend PASS와 absolute-width CHECK/FAIL을 분리하는 데 계산비용은 없다.

생성 Markdown은 probe-resonant 행에 `T=0.9962`, `|dT/dB|=0.002533 /µT`,
`29 pT/√Hz`를 기록한다
(`analysis/squeezing/resonant_hanle_squeezing_reference.md:20-25`). 현재
`compute_hanle(detuning=0)`를 다시 실행하면 같은 `B=-0.04 µT`에서
`T=0.737333`, `|dT/dB|=0.148987 /µT`다
(`analysis/squeezing/resonant_hanle_squeezing_reference.py:396-412`). 생성 문서는 measured CSV가
없고 compact/semi-quantitative라고 올바르게 경고하지만
(`analysis/squeezing/resonant_hanle_squeezing_reference.md:45-60`), 생성 시각, git commit/dirty,
`MagnetoScheme.cache_version`, source hash가 없다. 현재 산출물을 재생성하고 이 provenance를
추가하기 전에는 기존 `29`와 `5.25 pT/√Hz`를 reference 수치로 사용하면 안 된다.

### 5.5 문서/example의 신뢰 경계

- `README.md:15`, `38-39`, `61`은 scheme의 큰 범위를 간결하게 맞게 설명한다.
- 사용자 가이드는 Hanle/EIA/NMOR의 목적과 제한을 설명한다
  (`docs/Userguide/GABES_User_Guide_v2.html:619-630`, `829-833`). 그러나
  `userguide_assets/hanle.png`는 정규화 수정 전 `T≈0.998` scale의 그림이며 현재 기본
  `T(B=0)=0.868`과 불일치한다.
- 가이드의 “절대 스케일, 선폭, … NMOR 영점 교차 등을 자동 확인” 문구
  (`docs/Userguide/GABES_User_Guide_v2.html:911-914`)는 현재 테스트보다 강하다. 테스트는
  CG/CP/trace/positivity/sign/trend/internal normalization에는 강하지만 held-out measured
  absolute scale이나 문헌 선폭 tolerance를 pass/fail하지 않는다.
- 예제 config는 measured waist, detector noise, calibration CSV를 받을 수 있어 실험 workflow
  출발점으로 유용하다
  (`analysis/squeezing/resonant_hanle_experiment_config.example.json:10-46`). 그러나 실제 CSV가
  없는 기본 생성물을 absolute sensitivity validation으로 읽으면 안 된다.

## 6. 동작을 바꾸지 않는 순수 코드 최적화

### 6.1 최우선: kernel-side weighted coherence contraction

현재 Numba kernel은 `(B,v,n,n)` 전체 `rho`를 반환하고
(`gabes/kernels.py:299-353`, `356-433`), Python이 세 coherence를 추출해 velocity average한다
(`gabes/schemes/magneto.py:512-515`, `645-659`). 401×641, 8-level 반환 배열은 **251.0 MiB**다.

kernel에서 conditional trace 정규화 후 필요한 `chi_+`, `chi_-`, `chi_probe`의 Maxwell
가중합만 누적해 반환하면 세 401-point complex array, 약 **18.8 KiB**면 된다. LU solve와
관측량 정의는 그대로다. 현재 NumPy fallback/kernel parity test
(`tests/test_kernels.py:82-99`)를 세 susceptibility에 그대로 적용하면 behavior-preserving임을
검증할 수 있다.

### 6.2 affine generator를 real basis에서 조립

현재는 complex `C_xy + B*C_z` stack을 만든 뒤 wrapper가 전체 stack에 `U†LU`를 적용한다
(`gabes/schemes/magneto.py:485-501`; `gabes/kernels.py:339-352`, `415-432`). 실제 default
8-level generator의 401-point microbenchmark는 다음과 같다.

| 방식 | median |
|---|---:|
| complex stack 조립 후 real-basis 변환 | 82.53 ms |
| `C_xy`, `C_z`를 먼저 변환한 뒤 real affine stack 조립 | 8.17 ms |

- 조립 단계 속도 향상: **10.10×**
- generator 최대 절대차: `5.82e-11` (floating-point roundoff)
- complex stack: 25.06 MiB, real stack: 12.53 MiB

전체 runtime은 LU solve가 지배하므로 전체 10× 가속을 뜻하지 않지만, 두-region의 light/dark
stack allocation과 basis-transform 비용을 줄이는 안전한 최적화다. parity tolerance는 현재
kernel test에 고정해야 한다.

### 6.3 낮은 우선순위: angular-momentum matrix cache

`_hamiltonian()` 호출마다 같은 `(F_g,F_e)`의 `Fx,Fy,Fz`를 다시 만든다
(`gabes/schemes/magneto.py:538-545`; `gabes/zeeman.py:48-66`). 현재 장비에서 F=2/F=1 한 쌍은
호출당 약 36 µs다. 작은 `lru_cache`로 동작을 바꾸지 않고 없앨 수 있지만, 계산당 호출 횟수가
적어 앞의 두 최적화보다 우선순위가 낮다. `zeeman_manifold()`와 Hermitian basis는 이미 cache돼
있으므로 중복 제안하지 않는다 (`gabes/zeeman.py:69-70`; `gabes/core.py:62-93`).

## 7. 권고 실행 순서

1. **P1 semantics:** 파라핀 transit exchange와 light pure-dephasing 중복을 제거하고,
   `collisional_depol`을 명시적 ground pure dephasing 또는 실제 isotropic depolarization으로
   확정한다.
2. **P1 cheap trust:** physical-zero/command-zero 분리, default coarse/nonconverged status,
   open-hyperfine scope status, local odd-fit NMOR slope, 문헌 PASS/CHECK/FAIL을 구현한다.
3. **P1 polarization:** linear-polarization plane이 정의되지 않는 QWP 범위의 NMOR hero를 gate하고,
   검증된 Jones/Stokes thin-medium postprocess를 추가한다.
4. **P1 reproducibility/docs:** FWM–Hanle 산출물과 사용자 가이드 그림을 현재 물리로 재생성하고,
   model/git/dirty/source-hash provenance를 기록한다. 테스트 설명은 실제 검증 범위로 낮춘다.
5. **P1 coefficients:** gas/species/line별 sourced elastic broadening과 pressure shift를 추가한다.
   geometry-aware relaxation은 입력과 자료가 결정된 뒤 진행한다.
6. **성능:** kernel-side coherence contraction을 먼저, real-basis affine assembly를 다음으로 한다.
7. **Parked:** full VCC와 full open-hyperfine/repump reference는 held-out dataset과 runtime budget이
   생길 때만 opt-in으로 활성화한다.

## 8. 검증 기록

- registry 직접 확인: `sas -> lambda -> rydberg_eit -> magneto -> fwm`
- 공개 GitHub issues: 0개
- targeted:
  `python -m pytest -q tests/test_magneto.py tests/test_kernels.py tests/test_headless_observables.py tests/test_schemes_render.py tests/test_resonant_hanle_reference.py`
  → **72 passed in 15.69 s**
- `python tests/verify_hanle_eit_eia.py` 재실행:
  sign 전환과 TOC trend는 재현하지만 paraffin width는 문헌보다 약 9배 크고, buffer low-power
  width도 문헌보다 약 41% 크며 출력 label은 이를 과장한다.
- 고해상도 수렴, physical-zero, transit 중복, stale Hanle 예제와 optimization microbenchmark를
  오늘 working tree에서 재현했다.
- full suite: `python -m pytest -q` → **476 passed in 247.98 s (0:04:07)**
- 기존 dirty 작업은 보존했으며, 이번 자동화는 오늘 보고서와 automation memory만 변경한다.

## 9. 최종 판정

Scheme 4는 **유용한 실제 물리**를 구현한다. 특히 Zeeman CG 구조, TOC 기반 intrinsic EIA,
residual transverse field가 만드는 polarization sign switch, two-region Ramsey narrowing,
buffer-cell LCA, NMOR dispersion을 하나의 OBE 틀에서 연결한 점은 실험 조건 탐색과 초기 measured
trace fit에 가치가 있다.

다만 현재 기본 grid에서 물리 feature의 부호와 NMOR slope가 reference schedule과 수렴하지 않고,
relaxation channel 하나가 중복 적용되며, open hyperfine/repump·collision provenance·편광 전파가
빠져 있다. 생성 absolute sensitivity 산출물도 현재 코드와 불일치한다. 따라서 현 단계의 올바른
표현은 **정성·반정량 experimental planning/reference**이며, **absolute atomic magnetometry 또는
publication-grade linewidth/sensitivity reference는 아니다**.
