# 2026-09-08 Scheme 4 물리 검토 — Hanle / EIA / NMOR

## 1. 오늘의 선택과 현재 다섯 scheme

- 로컬 날짜: **2026-09-08** (`Asia/Seoul`)
- 계산: `n = (8 mod 5) + 1 = 4`
- 검토 대상: **Scheme 4 — Magneto-optics (Hanle / EIA / NMOR)**
- 현재 HEAD: `c0c46f0`

현재 dropdown registry 순서는 다음과 같다.

| n | registry key | 현재 정의 |
|---:|---|---|
| 1 | `sas` | Absorption OD/SAS: pump-off Doppler OD와 pump-on hyperfine-pumping SAS |
| 2 | `lambda` | 3-level Λ coherence: EIT / AT / CPT |
| 3 | `rydberg_eit` | 5S–5P–40D cascade EIT와 40D–39F microwave AT electrometry |
| 4 | `magneto` | ⁸⁷Rb D1 Zeeman manifold의 Hanle / EIA / NMOR |
| 5 | `fwm` | seeded ⁸⁵Rb D1 double-Λ mean-field 진단과 generic SFWM biphoton |

근거는 `gabes/schemes/__init__.py:12-25`의 설명과 `_SCHEMES` 순서다. 직접 실행한
`all_schemes()`도 `sas -> lambda -> rydberg_eit -> magneto -> fwm`을 반환했다.

## 2. 먼저 검색한 기존 제안, TODO, issue note

검토 전에 다음을 확인했다.

- Scheme 4 일일 보고서 10개:
  `docs/daily_report/2026-06-23_scheme-4_hanle-eia-nmor.md`부터
  `docs/daily_report/2026-08-23_scheme-4_hanle-eia-nmor.md`까지
- 통합 작업 목록 `docs/checklist.json`
- `docs`, `gabes`, `tests`, `analysis`의 `TODO`, `FIXME`, `issue`, `proposal`,
  `magneto`, `Hanle`, `EIA`, `NMOR` 검색
- `README.md`, `docs/Userguide/GABES_User_Guide_v2.html`
- 구현 `gabes/schemes/magneto.py`, `gabes/zeeman.py`, `gabes/kernels.py`,
  `gabes/core.py`, `gabes/constants.py`
- 테스트 `tests/test_magneto.py`, `tests/test_kernels.py`,
  `tests/test_headless_observables.py`, `tests/test_schemes_render.py`,
  `tests/test_resonant_hanle_reference.py`, `tests/verify_hanle_eit_eia.py`
- 실험 workflow 예제
  `analysis/squeezing/resonant_hanle_squeezing_reference.py`,
  `analysis/squeezing/resonant_hanle_experiment_config.example.json`, 생성 Markdown

별도 `.github` 디렉터리나 독립적인 로컬 issue/TODO 파일은 없었다. 코드 안의 직접 관련 TODO는
고정 Ne broadening scalar를 gas/species/line coefficient table과 pressure shift로 확장하라는
`gabes/constants.py:79-89` 하나다. 개선안이 없는 것이 아니라, 과거 보고서와 checklist에 이미
상당히 구체적으로 정리돼 있다.

현재 핵심 항목은 다음과 같다.

| 기존 항목 | 상태 | 위치 |
|---|---|---|
| lifetime / elastic dephasing / transit / depolarization 분리 | `done`, P1 | `docs/checklist.json:390-408` |
| two-region light-state 정규화와 density-once | `done`, P1 | `docs/checklist.json:411-430` |
| observable 수렴·physical zero·문헌 비교 신뢰 경계 | `ready`, P1 | `docs/checklist.json:433-452` |
| collision coefficient provenance와 pressure shift | `ready`, P1 | `docs/checklist.json:219-252` |
| geometry-aware buffer relaxation budget | `needs_decision`, P2 | `docs/checklist.json:455-480` |
| full velocity-changing-collision reference | `parked`, P2 | `docs/checklist.json:614-633` |

8월 23일 이후 commit `7244f25`가 Magneto claim-hardening을 저장소에 반영했다. 그러나 그
내용은 8월 23일 보고서가 당시 working tree에서 이미 평가한 natural/elastic 분리, 완전양의
dephasing, conditional light-state 정규화, density-once, external-validation status다. 따라서
오늘 확인한 현행 물리와 당시 보고서의 수치는 동일하다. 이후 Magneto 본체·Zeeman·테스트에
추가 변경은 없다.

## 3. 무엇을 실제로 구현하는가

### 3.1 유용한 실제 물리

1. ⁸⁷Rb D1의 선택한 `F_g ↔ F'_e`에 대해 허용 transition과 level 수를 검사하고 hyperfine
   Landé `g_F`를 계산한다 (`gabes/schemes/magneto.py:394-410`).
2. `zeeman_manifold()`는 `m_F` 상태, σ±/π Clebsch–Gordan 결합, 편광별 grouped
   spontaneous-emission operator `Σ_q`를 만든다 (`gabes/zeeman.py:69-123`). 같은 방출
   편광을 공유하는 excited coherence가 ground coherence로 전달되는 TOC가 살아 있어
   `F_e=F_g+1` intrinsic EIA 경로가 있다 (`tests/test_magneto.py:50-63,458-471`). 이는
   EIA의 coherence-transfer 메커니즘을 검증한
   [Failache et al.](https://arxiv.org/abs/quant-ph/0211065)의 물리와 방향이 맞다.
3. Ne optical broadening은 자연수명·branching·TOC를 바꾸지 않는 elastic optical
   dephasing으로 분리됐다 (`gabes/schemes/magneto.py:411-425`;
   `gabes/zeeman.py:132-149`). 관련 불변식은 `tests/test_magneto.py:66-102`가 고정한다.
4. longitudinal scan, longitudinal offset, residual transverse field와 QWP가 만든 σ± drive를
   Zeeman Hamiltonian에 함께 넣는다 (`gabes/schemes/magneto.py:440-460,538-556`).
5. paraffin 셀은 illuminated/dark 두 영역의 교환 OBE, buffer 셀은 isotropic reservoir reload를
   포함한 single-region OBE를 쓴다 (`gabes/schemes/magneto.py:468-510,568-642`).
6. paraffin light block을 trace-one conditional state로 정규화하고 occupation은 별도 diagnostic로
   둔다. local vapor density `N`은 감수율에 한 번만 들어간다
   (`gabes/schemes/magneto.py:502-534,558-565`;
   `tests/test_magneto.py:208-302`).
7. Maxwell velocity class를 평균하고, coarse optical opacity는 fine scalar Voigt average로
   보정한다 (`gabes/schemes/magneto.py:72-104,476-515`).
8. 관측량은 `α = k Imχ`, `T = exp(-αL)`,
   `θ = kL Re(χ_+ - χ_-)/4`로 계산한다 (`gabes/schemes/magneto.py:683-705`).

이것은 단순 Lorentzian toy curve가 아니다. selection rule, optical pumping, TOC, Zeeman
precession, transverse-field sign switching, wall-preserved Ramsey core와 Beer–Lambert readout을
실제 밀도행렬 계산으로 연결한다. 특히 ⁸⁷Rb D1 `F=2→F'=1` paraffin 셀에서 linear CPT가
circular magnetic-field-induced absorption으로 바뀌고 0.12/0.20 mG 폭이 관측됐다는
[Lee & Moon의 실험](https://opg.optica.org/josab/abstract.cfm?uri=josab-30-8-2301)과
비교 가능한 topology를 구현한다. NMOR를 atomic-vapor magnetometry와 연결하는 일반적 근거는
[Budker et al.의 review](https://journals.aps.org/rmp/abstract/10.1103/RevModPhys.74.1153)와
일치한다.

### 3.2 experimental-physics reference로서의 한계

- 선택한 한 `F_g ↔ F'_e` manifold 안에서 decay가 닫힌다. 실제 D1 open transition의 다른
  ground hyperfine manifold 누출과 repump가 없다 (`gabes/zeeman.py:94-130`).
- Ne는 `3.91 MHz/Torr` broadening scalar 하나뿐이며 gas/species/line/temperature provenance,
  pressure shift, diffusion, spin destruction, Dicke narrowing이 없다
  (`gabes/constants.py:79-89`).
- 기본 Doppler grid는 9 class, UI 최대도 41 class다
  (`gabes/schemes/magneto.py:295-302`). scalar opacity 보정은 B-dependent velocity dynamics의
  수렴을 보장하지 않는다.
- `line_strength`는 실험변수가 아니라 effective-OD calibration knob다
  (`gabes/schemes/magneto.py:290-294`). 독립 absolute transmission 없이 조정하면 예측력이
  생기지 않는다.
- 현행 NMOR 식은 약한 매질/선형 입력의 circular birefringence readout에는 유용하지만 임의
  ellipticity의 출력 ellipse를 Jones/Stokes로 전파하지 않는다.
- 코드도 absolute contrast/rotation/slope를 `external validation required`로 표시한다
  (`gabes/schemes/magneto.py:41-43,720-752`).

따라서 올바른 판정은 **실제 물리를 포함하는 유용한 정성·반정량 실험 설계 및 초기 fit
reference**다. QWP, 잔류장, 셀 완화, detuning, intensity 변화의 방향과 sign topology를
탐색하기에는 가치가 있다. 반면 absolute linewidth, contrast, NMOR slope 또는 pT/√Hz를
publication-grade reference로 인용하기에는 아직 부족하다.

## 4. 현행 수치 재검증

### 4.1 built-in 121×9 regime

feature는 중앙 absorption과 양쪽 scan wing 평균의 차이다. 음수는 dip, 양수는 peak다.

| regime | feature [m⁻¹] | 중앙 폭 [µT] | `T(B_cmd=0)` | NMOR slope [mrad/µT] |
|---|---:|---:|---:|---:|
| EIT dip | -1.584625 | 0.166258 | 0.868408 | -81.411 |
| EIA peak | +5.879812 | 0.223274 | 0.942600 | -15.757 |
| Buffer Hanle | -6.874297 | 4.591797 | 0.780505 | -5.654 |
| Buffer LCA | +0.014555 | 1.776939 | 0.662757 | -0.031 |
| NMOR | -4.117929 | 0.132017 | 0.909484 | -129.421 |

첫 실행은 Numba warm-up 포함 3.56 s, 이후 regime은 약 0.11–0.17 s였다.

### 4.2 deterministic 401×321→401×641 reference schedule

| regime | `(B,v)` | feature [m⁻¹] | 중앙 폭 [µT] | slope [mrad/µT] | runtime | 반환 `rho` |
|---|---:|---:|---:|---:|---:|---:|
| EIT | 401×321 | +0.319093 | 1.029954 | -114.254 | 11.70 s | 125.7 MiB |
| EIT | 401×641 | +0.357179 | 0.980545 | -116.791 | 24.45 s | 251.0 MiB |
| NMOR | 401×321 | -2.779918 | 0.110337 | -307.818 | 13.00 s | 125.7 MiB |
| NMOR | 401×641 | -2.816044 | 0.109904 | -314.549 | 26.58 s | 251.0 MiB |

- 기본 EIT는 9 class에서 **dip**이지만 321/641 class에서 **peak**로 부호가 뒤집힌다.
- 321→641에서 EIT 폭은 약 4.8%, slope는 약 2.2% 바뀌지만 feature amplitude는 약 11.9%
  바뀐다.
- 기본 NMOR slope `-129.4`는 641-class `-314.5 mrad/µT`보다 절댓값이 약 59% 작다.
- 고해상도를 interactive default로 강제하는 것은 비용상 부적절하다. 대신 기본 결과에
  `coarse/unconverged`를 표시하고 opt-in reference schedule을 제공하는 것이 맞다. 이는
  `magneto-observable-convergence-and-trust`의 acceptance와 정확히 일치한다
  (`docs/checklist.json:447-451`).

## 5. 남아 있는 P1 신뢰 문제

### 5.1 command zero와 physical zero 혼동

`b_offset_ut=0.25 µT`, 401×41에서 재검증한 값은 다음과 같다.

- command-zero feature: `-0.157925 m⁻¹`
- physical-zero feature: `-1.444534 m⁻¹`

Hamiltonian은 `B_physical = B_command + B_offset`을 올바르게 사용한다
(`gabes/schemes/magneto.py:440-444`). 그러나 observables는 displayed command axis
`raw["b_ut"]`에서 0에 가장 가까운 index를 hero로 쓴다
(`gabes/schemes/magneto.py:683-705`). physical-zero interpolation과 두 기준 병기는
`O(n_B)` postprocess이고 OBE 재풀이가 필요 없다.

### 5.2 paraffin transit rate가 exchange와 pure dephasing에 이중 적용

동일한 `transit_relax_khz`가

1. light→dark exchange `gamma_out`이 되고 (`gabes/schemes/magneto.py:429-432`),
2. light-region `gamma_gg` projective dephasing에도 들어간다
   (`gabes/schemes/magneto.py:433-435`).

exchange kernel 자체가 light block에 `-gamma_out I`를 적용한다
(`gabes/kernels.py:377-386`). light pure-dephasing만 제거하고 exchange를 유지한 201×41 진단은
다음과 같다.

| | 현행 | exchange-only 진단 |
|---|---:|---:|
| feature [m⁻¹] | -1.442743 | -4.847639 |
| 중앙 폭 [µT] | 0.152061 | 0.144810 |
| `T(B=0)` | 0.863261 | 0.934032 |
| slope [mrad/µT] | -88.192 | -230.719 |

`chi_probe` 최대 상대차는 **48.8%**다. transit은 two-region exchange가 전담하고 추가 in-beam
pure dephasing은 별도 이름과 근거를 가져야 한다. Liouville 차원과 solve 수는 그대로이므로
runtime overhead는 0이다.

`collisional_depol_khz`도 명칭상 spin destruction/depolarization처럼 보이지만 구현은 모든
ground population을 보존하는 projective pure dephasing이다
(`gabes/zeeman.py:132-141`; `tests/test_magneto.py:142-152`). 현재 동작이 의도라면
`ground_pure_dephasing_khz`로 이름을 바꾸는 것이 계산비용 0의 안전한 정리다. 실제
spin-destruction을 원한다면 orientation/alignment relaxation channel과 multipole decay test가
필요하다. 같은 Hilbert 차원 안에서 구현 가능하므로 solve 차수는 늘지 않는다.

### 5.3 순수 원편광에도 NMOR “편광면 회전”을 보고

`_qwp_drive_weights(45°)`는 한 circular component가 정확히 0인 순수 원편광을 만든다
(`gabes/schemes/magneto.py:117-128`). 그런데 `_coherences()`는 각 circular coherence를 입력
성분이 아니라 공통 `Ω`로 나누고 (`gabes/schemes/magneto.py:645-659`), observables는 모든
QWP에서 같은 plane-rotation 식을 쓴다 (`gabes/schemes/magneto.py:683-695`).

QWP 45°, 401×41 진단에서도 `rotation(B=0)=5.451 mrad`,
`slope=-21.284 mrad/µT`가 출력됐다. 순수 원편광에는 선형 편광면의 major-axis angle이
정의되지 않는다. 우선 해당 범위를 `unsupported/diagnostic`으로 gate하는 것은 `O(1)`이다.
다음 단계의 입력 σ± amplitude 정규화와 thin-medium Jones/Stokes postprocess도 `O(n_B)`로
가능하다. self-consistent polarization propagation은 별도 고비용 연구 문제로 남기는 것이 맞다.

### 5.4 문헌 비교 스크립트가 오차를 과장해 통과시킴

`python tests/verify_hanle_eit_eia.py`의 현행 출력은 다음과 같다.

- paraffin linear/circular: `1.193/1.737 mG` 대 논문 `0.12/0.20 mG`
  (약 `9.94×/8.69×`)
- buffer low-power LCA: `3.389 mG` 대 논문 `2.4 mG` (약 `1.41×`)

그런데 출력은 첫 비교를 “same sub-mG order”, sign switch를 `MATCH`, buffer 폭은 reference에
“reaches”라고 쓴다 (`tests/verify_hanle_eit_eia.py:243-271`). 1.19–1.74 mG는 sub-mG도 아니다.
sign/trend PASS와 absolute-width CHECK/FAIL을 분리하는 데 추가 계산비용은 없다.

### 5.5 생성 예제 보고서가 현행 코드와 불일치

생성 Markdown은 probe-resonant 행에 `T=0.9962`,
`|dT/dB|=0.002533 /µT`, `29 pT/√Hz`를 기록한다
(`analysis/squeezing/resonant_hanle_squeezing_reference.md:20-25`). 현재
`compute_hanle(detuning=0)` 재실행은 같은 `B=-0.04 µT`에서
`T=0.737333`, `|dT/dB|=0.148987 /µT`다. 계산 경로는
`analysis/squeezing/resonant_hanle_squeezing_reference.py:365-415`에 있다.

예제는 measured CSV가 없음을 명시하지만 생성 시각, git commit/dirty,
`MagnetoScheme.cache_version`, source/config hash를 기록하지 않는다. 현재 산출물을 재생성하고
portable provenance를 추가하기 전에는 기존 `29`와 `5.25 pT/√Hz`를 reference 수치로 쓰면 안
된다. provenance 기록은 hash/문자열 생성의 `O(file size)`이며 OBE 비용과 비교해 무시할 만하다.

## 6. 기존 개선안의 계산비용과 물리 보존성

| 기존 제안 | 현재 판정 | 계산비용 | 물리 보존/권고 |
|---|---|---:|---|
| natural lifetime / elastic dephasing 분리 | 구현 완료 | 같은 `(B,v)` solve; 실질 0 | 올바른 수정. emission/TOC가 Ne pressure에 오염되지 않는다. |
| conditional light state / density-once | 구현 완료 | solve 후 trace normalize `O(n_B n_v n²)` | 유의미한 solve 증가 없이 normalization ambiguity를 제거한다. |
| physical zero, model-scope, local slope, trust status | 아직 `ready` | `O(1)`–`O(n_B)` | OBE를 바꾸지 않는 metadata/postprocess이므로 즉시 적용할 가치가 높다. |
| deterministic 641-class reference | 아직 `ready` | 약 24–27 s, 반환 `rho` 251 MiB | interactive default가 아니라 opt-in certification에만 둔다. |
| paraffin transit 중복 제거 | 일일 보고서에만 존재 | solve 차원·횟수 동일 | 물리 영향은 크지만 계산 overhead 0. checklist에 독립 acceptance로 승격할 필요가 있다. |
| NMOR polarization gate/Jones–Stokes | 일일 보고서에만 존재 | gate `O(1)`, thin postprocess `O(n_B)` | 먼저 unsupported gate, 이후 검증된 thin-medium readout을 추가한다. |
| collisional coefficient/pressure shift | `ready` | table lookup·scalar shift; 미미 | source/units/sign/uncertainty와 unsupported 조합을 정의하면 기존 compact solver를 보존한다. |
| geometry-aware relaxation budget | `needs_decision` | geometry→scalar rates 자체는 미미 | 입력·계수·경계조건이 없으면 false precision이다. 결정 후 opt-in derived mode로 구현한다. |
| full VCC | `parked` | velocity class 결합으로 시간·메모리 크게 증가 | kernel·dataset·runtime budget이 정해질 때만 reference mode로 활성화한다. |
| full open-hyperfine/repump | 아직 별도 선정 필요 | Hilbert/Liouville 차원 증가, 검증면 큼 | held-out spectrum이 있는 opt-in reference로만 추진한다. |

## 7. 동작을 바꾸지 않는 순수 코드 최적화

### 7.1 최우선: kernel-side weighted coherence contraction

현재 kernel은 `(B,v,n,n)` 전체 `rho`를 반환하고
(`gabes/kernels.py:299-353,356-433`), Python이 `chi_+`, `chi_-`, `chi_probe`를 추출해 Maxwell
평균한다 (`gabes/schemes/magneto.py:512-515,645-659`). 401×641, 8-level에서 반환 complex
`rho`는 **251.0 MiB**다.

kernel 내부에서 trace normalization 뒤 필요한 세 coherence의 velocity-weighted sum만 누적하면
최종 반환은 세 401-point complex array, 약 **18.8 KiB**다. LU 방정식과 감수율 정의는
그대로라 physics/behavior를 보존할 수 있다. `kernels.affine_scan_chi()`가 이미 real-basis solve
후 coherence를 kernel 안에서 contraction하는 선례를 제공한다
(`gabes/kernels.py:436-518`). 현재 NumPy/kernel parity를 세 susceptibility에 고정해야 한다.

### 7.2 affine generator를 real basis에서 조립

현행은 complex `C_xy + B C_z` stack을 만든 뒤 전체 stack을 `U†LU`로 변환한다
(`gabes/schemes/magneto.py:485-501`; `gabes/kernels.py:424-432`). 8-level, 401-point 재벤치:

| 방식 | 평균 |
|---|---:|
| complex stack 조립 후 real-basis 변환 | 83.04 ms |
| `C_xy`, `C_z`를 한 번 변환 후 real affine stack 조립 | 8.26 ms |

- 조립 단계 약 **10.05×** 단축
- generator 최대 절대차 `5.82e-11` (부동소수 roundoff)
- stack 메모리 `25.06 MiB → 12.53 MiB`

전체 runtime은 LU가 지배하므로 전체 10× 가속을 뜻하지 않는다. 그러나 기존 real solve와
정확히 맞는 low-risk 최적화다.

### 7.3 낮은 우선순위: angular-momentum matrix cache

`_hamiltonian()` 호출마다 같은 `(F_g,F_e)`의 `Fx,Fy,Fz`를 다시 만든다
(`gabes/schemes/magneto.py:538-545`; `gabes/zeeman.py:48-66`). 현재 한 F=2/F=1 pair 생성은
약 **37.2 µs**다. `lru_cache`로 동작을 바꾸지 않고 없앨 수 있지만 계산당 호출 수가 적어
앞의 두 항목보다 우선순위가 낮다. `zeeman_manifold()`와 Hermitian basis는 이미 cache된다
(`gabes/zeeman.py:69-70`; `gabes/core.py:61-94`).

## 8. 문서·테스트·예제 적합성

- `tests/test_magneto.py`는 CG, emission normalization, TOC, elastic dephasing, CP
  ground dephasing, trace/positivity, light-state normalization, density-once, sign/trend와 NMOR
  zero crossing을 폭넓게 검증한다 (`tests/test_magneto.py:28-158,208-502`). 내부 물리 불변식
  테스트는 강하다.
- 그러나 641-class convergence, physical-zero hero, QWP 45° NMOR invalidity, open-hyperfine
  scope, 문헌 linewidth tolerance와 stale generated output은 자동 test가 아니다.
- 사용자 가이드는 큰 범위를 맞게 설명하지만 Magneto 한계는 “실험별 조정 필요” 한 줄로 너무
  약하다 (`docs/Userguide/GABES_User_Guide_v2.html:852-854`). 검증 섹션의 “절대 스케일,
  선폭, … NMOR 영점 교차 등을 자동으로 확인”은 실제 검증 범위보다 강하다
  (`docs/Userguide/GABES_User_Guide_v2.html:935-937`).
- `userguide_assets/hanle.png`는 2026-06-08 생성물로 현행 normalization/scale 이전 그림이다.
  오늘 기본값과 함께 재생성하고 commit/hash를 기록해야 한다.
- `analysis/squeezing/resonant_hanle_experiment_config.example.json`은 measured waist, detector
  noise, CSV calibration을 받을 수 있어 workflow 출발점으로 유용하다. 다만 기본 noise가 모두
  0이고 measured CSV도 비어 있으므로 absolute sensitivity validation 예제는 아니다.

## 9. 권고 우선순위

1. **P1 trust UI:** physical/command zero 분리, default coarse/unconverged status,
   open-hyperfine scope, local odd-fit slope, 문헌 PASS/CHECK/FAIL을 먼저 구현한다.
2. **P1 relaxation semantics:** paraffin transit exchange와 light pure dephasing의 이중 적용을
   제거하고 `collisional_depol`의 실제 의미를 pure dephasing 또는 isotropic depolarization으로
   확정한다.
3. **P1 polarization:** 순수 원편광과 near-circular input에서 plane-rotation hero를 gate한 뒤,
   검증된 thin-medium Jones/Stokes readout을 추가한다.
4. **P1 reproducibility:** FWM–Hanle Markdown/figures와 user-guide Hanle 그림을 현재 코드로
   재생성하고 model/git/dirty/config/source hash provenance를 저장한다.
5. **P1 coefficients:** gas/species/line별 sourced elastic broadening과 pressure shift를 추가한다.
6. **성능:** kernel-side weighted contraction, real-basis affine assembly 순으로 진행한다.
7. **보류:** full VCC, geometry-derived relaxation, open-hyperfine/repump는 데이터와 runtime
   budget이 정해진 opt-in reference로만 활성화한다.

## 10. 검증 기록

- registry 직접 확인: `sas -> lambda -> rydberg_eit -> magneto -> fwm`
- built-in 5 regime, 401×321/641 convergence, physical zero, circular-NMOR,
  paraffin transit 중복과 생성 예제 mismatch를 현행 HEAD에서 재현
- `python tests/verify_hanle_eit_eia.py`: sign/trend는 재현하지만 paraffin linewidth는 문헌보다
  약 9배 크고 출력 라벨이 이를 과장함
- targeted:
  `python -m pytest -q tests/test_magneto.py tests/test_kernels.py tests/test_headless_observables.py tests/test_schemes_render.py tests/test_resonant_hanle_reference.py tests/test_docs_consistency.py`
  → **73 passed in 18.59 s**
- product/analysis:
  `python -m pytest -q tests analysis/rydberg_cell_heating`
  → **595 passed in 306.18 s**
- root `python -m pytest -q`는
  `references/AutoOD/test_autood_gui_adjustments.py` 수집 중 선택적·미선언 `sympy`가 없어 중단
- 기존 dirty/staged 사용자 작업은 보존했다. 이번 자동화는 오늘 보고서와 automation memory만
  변경한다.

## 11. 최종 판정

Scheme 4는 **유용한 실제 물리**를 구현한다. Zeeman CG 구조, grouped-emission TOC,
residual transverse field에 따른 polarization sign switching, two-region Ramsey narrowing,
buffer-cell LCA와 NMOR dispersion을 하나의 OBE 틀로 연결한 점은 실험 조건 탐색과 초기 trace
fit에 충분한 가치가 있다.

그러나 기본 grid에서 feature 부호와 NMOR slope가 reference schedule에 수렴하지 않고,
paraffin transit rate가 두 완화 경로에 중복 적용되며, physical zero·open hyperfine·편광 전파·
collision provenance가 불완전하다. 생성 absolute-sensitivity 산출물도 현재 코드와 맞지 않는다.
따라서 현 단계의 정확한 표현은 **정성·반정량 experimental planning/reference**이며,
**absolute atomic magnetometry 또는 publication-grade linewidth/sensitivity reference는 아니다**.
