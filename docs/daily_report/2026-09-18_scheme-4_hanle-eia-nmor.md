# 2026-09-18 Scheme 4 물리 검토 — Hanle / EIA / NMOR

## 1. 오늘의 선택과 현재 다섯 scheme

- 로컬 날짜: **2026-09-18**, Asia/Seoul
- 계산: **n = (18 mod 5) + 1 = 4**
- 검토 대상: **Scheme 4 — Magneto-optics (Hanle / EIA / NMOR)**
- 현재 HEAD: **2998dc1**
- 검토 기준: 현재 working tree. 기존 staged/unstaged/untracked 사용자 변경은 보존했다.

현재 실행 registry 순서는 다음과 같다.

| n | registry key | 현재 정의 |
|---:|---|---|
| 1 | sas | Absorption OD/SAS: pump-off Doppler OD와 pump-on hyperfine-pumping SAS |
| 2 | lambda | 3-level Λ coherence: EIT / AT / CPT |
| 3 | rydberg_eit | cascade Rydberg-EIT와 microwave AT electrometry |
| 4 | magneto | ⁸⁷Rb D1 Zeeman manifold의 Hanle / EIA / NMOR |
| 5 | fwm | seeded ⁸⁵Rb D1 double-Λ mean-field와 generic SFWM biphoton |

근거는 gabes/schemes/__init__.py:14-24의 설명과 생성 순서다. README.md:8-16도 같은
제품 순서를 쓴다.

## 2. 먼저 검색한 기존 개선안, TODO, issue note

코드를 평가하기 전에 다음을 확인했다.

- 기존 Scheme 4 일일 보고서와 최근 2026-09-08 보고서
- docs/checklist.json의 Magneto, collision, Doppler 관련 항목
- docs, gabes, tests, analysis 안의 TODO, FIXME, issue, proposal, Hanle, EIA, NMOR 검색
- README와 사용자 가이드
- 구현 gabes/schemes/magneto.py, gabes/zeeman.py, gabes/pole_doppler.py,
  gabes/observables.py, gabes/constants.py
- 테스트 tests/test_magneto.py, tests/test_magneto_doppler.py,
  tests/test_kernels.py, tests/test_resonant_hanle_reference.py,
  tests/verify_hanle_eit_eia.py
- 예제 analysis/squeezing/resonant_hanle_squeezing_reference.py와 config/생성 Markdown
- pole-residue 개발 기록 analysis/pole_reuse/DEVLOG.md와 GATE.md

별도 .github 디렉터리나 독립 로컬 issue 파일은 없다. 공개 GitHub API로 현재 저장소의
open/closed issue도 0건이었다. 직접 관련 코드 TODO는 고정 Ne 계수를 gas/species/line
table과 pressure shift로 확장하라는 gabes/constants.py:80-89다.

따라서 **기존 제안이 없는 상태가 아니다.** 이미 다음 항목이 구체적으로 기록돼 있다.

| 기존 항목 | 상태 | 근거 |
|---|---|---|
| collision coefficient provenance와 pressure shift | P1, ready, large | docs/checklist.json:2267-2300 |
| observable 수렴·physical zero·신뢰 표시 | P1, ready, medium | docs/checklist.json:2485-2503 |
| geometry-aware buffer relaxation budget | P2, needs_decision, large | docs/checklist.json:2507-2532 |
| full velocity-changing-collision kernel | P2, parked, research | docs/checklist.json:2672-2689 |
| pole-residue Doppler 재사용 | P3, done | docs/checklist.json:2900-2915 |
| exact-Doppler 뒤 EIT-dip preset 재검토 | P2, needs_decision, small | docs/checklist.json:2936-2953 |

이번 검토는 새 거대 모델을 제안하기보다 이 항목들이 현재 코드와 수치에서 여전히 맞는지,
어떤 것은 이미 해결됐는지, 물리를 유지하면서 싸게 고칠 수 있는지를 재평가했다.

## 3. 현행 코드가 구현하는 실제 물리

### 3.1 원자와 광학 결합

1. ⁸⁷Rb D1의 Fg ∈ {1,2}, Fe ∈ {1,2} 네 transition을 지원하고 line strength와
   저자기장 hyperfine Landé gF를 계산한다
   (gabes/schemes/magneto.py:21-30,55-104,384-400).
2. zeeman_manifold()는 mF 상태, σ±/π Clebsch–Gordan 결합과 편광별 grouped
   spontaneous-emission operator를 만든다
   (gabes/zeeman.py:17-45,69-123).
3. grouped emission은 excited coherence를 ground coherence로 옮기는 TOC를 보존한다.
   intrinsic EIA 경로는 tests/test_magneto.py:37-63,458-471에서 직접 검증된다.
4. QWP 각도를 실제 complex σ± drive weight로 바꾸고, longitudinal scan/offset와 residual
   transverse field를 같은 Zeeman Hamiltonian에 넣는다
   (gabes/schemes/magneto.py:107-118,430-450,555-574).

### 3.2 완화와 셀 모델

- buffer cell: 자연 방출은 고정하고 Ne broadening을 elastic optical dephasing으로,
  reservoir reload와 ground pure dephasing을 별도 채널로 둔다
  (gabes/schemes/magneto.py:401-417; gabes/zeeman.py:125-149).
- paraffin cell: illuminated/dark 두 영역을 exchange시키고 wall coherence lifetime을 dark
  Liouvillian에 둔다. light block은 trace-one conditional state로 정규화하며 occupation은
  별도 diagnostic이다
  (gabes/schemes/magneto.py:418-428,576-635).
- local vapor density는 susceptibility에 한 번만 적용한다
  (gabes/schemes/magneto.py:530-552; gabes/observables.py:556-570).

### 3.3 Doppler와 관측량

현재 Doppler-on은 더 이상 9-class 근사가 아니다. 각 B 행을 real Hermitian basis의 affine
shift-invert 문제로 만들고 pole-residue와 Faddeeva 함수로 **절단 없는 Maxwell 평균**을
계산한다. direct v=0 solve guard에 실패한 행만 수렴된 grid로 fallback한다
(gabes/schemes/magneto.py:45-52,662-747).

개발 gate는 다음을 기록한다.

- 수렴된 ±6σ, 1 m/s reference 대비 최대 χ 상대차 1.6×10⁻⁷
- FWHM/중앙폭 0.5% 이내, NMOR slope/peak 0.1% 이내
- guard 최대 4.9×10⁻¹², fallback 0
- reference 대비 최소 71배 가속

근거: analysis/pole_reuse/DEVLOG.md:136-156,
analysis/pole_reuse/GATE.md:42-63. 현재 tests/test_magneto_doppler.py:54-114도 exact/reference,
off-path, fallback, provenance를 회귀 고정한다.

관측량은

- α = k Imχ
- T = exp(−αL)
- θ = kL Re(χ+ − χ−)/4

로 계산한다 (gabes/schemes/magneto.py:788-812). 단순 Lorentzian 그림이 아니라 selection rule,
optical pumping, TOC, Zeeman precession, transverse-field sign switching, two-region Ramsey
core와 Beer–Lambert readout을 하나의 OBE에 연결한다.

## 4. experimental atomic-physics reference 판정

### 4.1 유용한 범위

이 모델은 다음 용도에 실제로 유용하다.

- linear/circular/elliptical drive에 따른 Hanle dip, EIA/MIA peak, crossover의 방향 판정
- 잔류 transverse field와 longitudinal offset의 영향 탐색
- paraffin wall coherence와 buffer-cell ground relaxation의 상대 추세
- detuning, intensity, cell length, QWP 각도에 대한 초기 실험 설계
- 실험 trace의 초기 sign/topology 확인과 반정량 fit seed

코드 자체도 Lee–Moon paraffin-cell magnetic-field-induced absorption, Moon–Kim Ramsey
EIA/MIT, Yu–Lee–Moon buffer-gas LCA를 참고문헌으로 명시한다
(gabes/schemes/magneto.py:362-381).

### 4.2 publication-grade reference가 아닌 이유

1. **한 hyperfine pair 안에서 닫힌 모델이다.** 실제 D1 open transition의 다른 ground F로
   누출, repump, 전체 excited hyperfine interference가 없다
   (gabes/zeeman.py:69-130).
2. **collision 계수가 빈약하다.** Ne는 3.91 MHz/Torr scalar 하나며 pressure shift,
   gas/line/temperature provenance, diffusion, spin destruction, Dicke narrowing이 없다
   (gabes/constants.py:80-89).
3. **absolute scale은 외부 검증 전이다.** 코드도 transmission contrast와 NMOR
   rotation/slope에 external validation required를 표시한다
   (gabes/schemes/magneto.py:41-46,825-829).
4. **NMOR 식의 편광 적용 범위가 좁다.** 현재 식은 약한 매질과 선형 입력의 circular
   birefringence readout에는 유용하지만 임의 ellipticity를 Jones/Stokes로 전파하지 않는다.
   QWP=45°는 한 circular 성분이 0인데도 같은 plane-rotation hero를 낸다
   (gabes/schemes/magneto.py:107-118,750-764,788-823).
5. **공간 평균이 없다.** Gaussian beam intensity, field gradient, cell geometry와 diffusion을
   직접 풀지 않는다.
6. **sensitivity는 원자 모델 단독 결과가 아니다.** 별도 squeezing workflow의 optical power,
   detector efficiency/noise, calibration CSV가 필요하다.

따라서 정확한 표현은 **실제 물리를 포함한 정성·반정량 experimental-planning/reference**다.
절대 linewidth, contrast, 회전각, pT/√Hz를 독립 publication reference로 인용할 단계는 아니다.

## 5. 오늘의 수치 재검증

### 5.1 built-in 121-point, exact-Maxwell regimes

warm median은 같은 process에서 여러 번 실행한 현재 PC 값이다.

| regime | 현재 feature | α 중앙−wing [m⁻¹] | 중앙 폭 [µT] | T(Bcmd=0) | NMOR slope [mrad/µT] | warm median |
|---|---|---:|---:|---:|---:|---:|
| EIT dip | **crossover** | +0.357285 | 0.981271 | 0.787789 | −102.569 | 58.37 ms |
| EIA peak | peak | +14.9300 | 0.363515 | 0.856642 | −23.566 | 58.44 ms |
| Buffer Hanle | dip | −4.670913 | 4.30894 | 0.757571 | −4.394 | 30.30 ms |
| Buffer LCA | crossover | +0.012284 | 1.69314 | 0.629190 | −0.041 | 73.10 ms |
| NMOR | dip | −2.816494 | 0.113283 | 0.827167 | −252.147 | 58.63 ms |

모든 행의 exact pole guard는 최대 2.6×10⁻¹²였고 fallback은 0이었다. 정확 Doppler 평균은
성공했다. 반면 “EIT dip” 버튼은 이름과 달리 수렴된 모델에서 crossover다. 이는 이미
docs/checklist.json:2936-2953에 기록된 현상이다.

### 5.2 B-grid와 NMOR slope 수렴

현재 slope는 중앙점의 3-point np.gradient다
(gabes/schemes/magneto.py:800,817-823).

| scan points | EIT 중앙 폭 [µT] | EIT slope [mrad/µT] | NMOR 중앙 폭 [µT] | NMOR 3-point slope [mrad/µT] |
|---:|---:|---:|---:|---:|
| 121 | 0.981271 | −102.569 | 0.113283 | −252.147 |
| 201 | 0.980818 | −112.214 | 0.110669 | −293.582 |
| 401 | 0.980476 | −116.811 | 0.109904 | −314.602 |
| 801 진단 | 0.980472 | −118.005 | 0.109482 | −320.119 |

NMOR 5-point central derivative는 121/201/401/1601점에서 각각
−285.285/−314.868/−321.608/−321.972 mrad/µT였다. 즉 기본 121점 3-point slope는
고해상도 5-point 값보다 절댓값이 약 21.7% 작다. 401점 5-point는 1601점 대비 약 0.11%다.

따라서 기존 magneto-observable-convergence-and-trust 제안의
“locally converged odd fit 또는 higher-order derivative”가 정확하다. 5-point postprocess는
O(1)이고 사실상 무료다. 201/401점으로 scan을 늘리면 solve 비용은 약 82/153 ms로 선형
증가하지만 여전히 interactive 범위다. 기본값을 무조건 늘리기보다 slope convergence status와
reference preset을 병기하는 편이 안전하다.

### 5.3 paraffin transit rate 이중 적용

동일한 transit_relax_khz가

1. light→dark exchange γout이 되고
2. light-region atom의 γgg pure dephasing에도 들어간다

(gabes/schemes/magneto.py:419-424). exchange kernel 자체도 light block에 −γout I를 적용한다
(gabes/schemes/magneto.py:621-625).

201-point exact-Maxwell EIT preset에서 light γgg만 0으로 만든 비저장 진단:

| | 현행 | exchange-only 진단 |
|---|---:|---:|
| feature [m⁻¹] | +0.357285, crossover | −4.428289, dip |
| 중앙 폭 [µT] | 0.980818 | 0.123955 |
| T(B=0) | 0.787789 | 0.857372 |

χ의 최대 상대 변화는 **37.6%**다. 현행 EIT-preset 문제는 단순 label 문제가 아니라
transit semantics와 얽혀 있다. 먼저 transit exchange와 추가 in-beam pure dephasing의 물리
의미를 분리해야 한다. 행렬 차원과 solve 수는 같으므로 실행시간 overhead는 0이다.

### 5.4 command zero와 physical zero

Hamiltonian은 Bphysical = Bcommand + Boffset를 올바르게 사용한다
(gabes/schemes/magneto.py:430-444). 그러나 hero는 command axis의 0 index다
(gabes/schemes/magneto.py:788-812).

b_offset=0.25 µT, 401점 EIT 진단에서:

- command zero, 즉 Bphysical=+0.25 µT: α=23.90693 m⁻¹, T=0.787360
- physical zero, 즉 Bcommand=−0.25 µT: α=23.85254 m⁻¹, T=0.787789

이 조건에서는 수치 차이가 작지만 의미는 다르다. physical-zero interpolation과 두 기준
병기는 O(nB) postprocess이며 OBE 재풀이가 없다.

### 5.5 문헌 비교 스크립트의 과한 판정

python tests/verify_hanle_eit_eia.py의 오늘 출력:

- paraffin linear width 0.950 mG 대 논문 0.12 mG: 약 7.9배
- paraffin circular width 2.305 mG 대 논문 0.20 mG: 약 11.5배
- buffer low-power width 3.349 mG 대 reference 2.4 mG: 약 39.5% 큼

그런데 스크립트는 sign switch를 MATCH라 하고 첫 두 폭을 same sub-mG order라고 쓴다
(tests/verify_hanle_eit_eia.py:243-271). sign/trend PASS와 absolute-width CHECK/FAIL을
분리해야 한다. 계산비용은 0이다.

## 6. 기존 개선안의 비용과 물리 보존성

| 기존 제안 | 오늘의 판정 | 계산비용 | 물리 보존/권고 |
|---|---|---:|---|
| exact pole-residue Doppler | **완료·검증됨** | 기본 30–73 ms, reference 대비 최소 71× | 절단 없는 같은 Maxwell 물리를 더 정확하고 빠르게 계산한다. |
| observable convergence/trust | **P1 유지** | slope/zero/status O(1)–O(nB); 401점은 약 153–181 ms | 모델을 바꾸지 않는 numerical/trust 수정. 가장 값싸다. |
| EIT-dip preset 재검토 | **P2 needs_decision 유지** | relabel 0, retune 탐색만 추가 | transit semantics를 먼저 결정한 뒤 문헌 조건으로 retune해야 한다. |
| collision provenance/pressure shift | **P1 유지** | lookup·scalar shift·같은 dephasing 행렬: 미미 | source, units, sign, uncertainty와 unsupported 조합을 명시하면 compact physics를 보존한다. |
| geometry-aware relaxation | **P2 결정 필요** | geometry→rate 계산 자체는 미미 | geometry, diffusion, wall, spin-destruction 데이터 없이 pressure-only 자동값을 만들면 false precision이다. |
| full VCC | **parked 유지** | velocity class 결합으로 시간·메모리 크게 증가 | separable exact-Maxwell 경로를 깨므로 opt-in reference만 가능하다. |
| open hyperfine/repump | 큰 별도 과제 | Hilbert/Liouville 차원과 검증면 증가 | held-out spectrum이 있는 reference mode 전용이 적절하다. |
| transit exchange/dephasing 분리 | 과거 일일 보고서 제안, 아직 독립 checklist 없음 | solve 차원·횟수 동일 | 물리 영향 37.6%, overhead 0. 독립 acceptance로 승격할 가치가 높다. |
| circular-input NMOR gate/Jones–Stokes | 과거 일일 보고서 제안 | gate O(1), thin-medium postprocess O(nB) | 먼저 unsupported gate, 이후 검증된 Jones/Stokes readout. |

## 7. 동작을 바꾸지 않는 순수 코드 최적화

### 7.1 최우선: exact path에서 B-stack을 지연 생성

compute()는 exact Doppler를 선택해도 먼저 L0_all과 Ld_all 전체 complex B-stack을 만든다
(gabes/schemes/magneto.py:468-485). 그러나 정상 exact path는 C_xy, C_z, Cd_xy만 사용하고
이 stack은 guard fallback이 있을 때만 grid_average가 필요로 한다
(gabes/schemes/magneto.py:487-516,662-747).

121-point paraffin F=2→F′=1에서 이 불필요한 두 stack은:

- median 조립시간 약 **9.63 ms**
- 메모리 약 **15.13 MiB**

F=2→F′=2이면 약 **36.9 MiB**다. 오늘 모든 built-in regime은 fallback 0이었다. Doppler-off
또는 실제 failed rows가 생길 때만 동일 수식으로 필요한 row를 생성하면 결과는 동일하고
default runtime의 약 10 ms와 peak memory를 줄일 수 있다. exact-path parity, forced-fallback
parity, off-path bitwise test로 고정하면 된다.

### 7.2 angular-momentum matrix cache

_hamiltonian()은 같은 Fg/Fe의 Fx,Fy,Fz를 호출 때마다 다시 만든다
(gabes/schemes/magneto.py:555-574; gabes/zeeman.py:48-66). F-pair당 lru_cache를 적용하면
동작은 동일하다. 이전 측정은 F=2/F=1 pair당 약 37 µs라 우선순위는 낮다.

### 7.3 BLAS thread-limit controller 재사용

core.blas_single_thread()는 매번 threadpool_limits context를 새로 만든다
(gabes/core.py:18-34). 오늘 빈 context의 steady median은 약 0.406 ms, 이미 만든 controller의
limit context는 약 0.002 ms였다. 작은 비용이므로 낮은 우선순위다. 라이브러리 late-load와
thread-safety를 확인한 뒤 적용해야 한다.

### 7.4 과거 최적화 제안의 현재 상태

과거의 kernel-side velocity-weighted coherence contraction과 real-basis affine assembly는
옛 dense velocity grid에서 중요했다. 현재 default exact path는 이미 세 readout만 반환하고
generator를 real basis로 한 번 변환한다
(gabes/schemes/magneto.py:681-727). 따라서 기본 경로에서는 사실상 구현 완료다. 남은
dense contraction 최적화는 rare fallback/off 경로에만 해당해 우선순위가 내려갔다.

## 8. 문서·테스트·예제 적합성

### 8.1 테스트

tests/test_magneto.py는 CG, emission normalization, TOC, elastic dephasing, relaxation-channel
분리, positivity, conditional normalization, density-once, field offset, Hanle/EIA/LCA sign,
NMOR zero crossing을 폭넓게 검증한다
(tests/test_magneto.py:28-158,208-526).

tests/test_magneto_doppler.py는 exact average와 수렴 grid, guard/fallback, Doppler-off,
provenance를 검증한다 (tests/test_magneto_doppler.py:54-114). 수치 kernel 내부 검증은 강하다.

아직 자동 고정되지 않은 항목:

- exact Doppler에서 모든 preset label이 실제 feature와 일치하는지
- 5-point/odd-fit NMOR slope와 B-grid convergence
- physical-zero hero
- QWP near-circular에서 rotation의 지원 범위
- open-hyperfine/repump scope status
- 문헌 absolute linewidth tolerance
- 생성 example의 code/config/provenance 일치

### 8.2 사용자 문서

- User Guide의 Magneto 설명은 실제 큰 구조와 맞지만 한계는 “실험별 조정 필요” 수준으로
  너무 약하다 (docs/Userguide/GABES_User_Guide_v2.html:628-640,852-854).
- runtime 표는 여전히 약 0.35 s와 “B scan × velocity”라고 쓴다. exact pole path의 현재
  30–73 ms와 velocity-grid-free 구조로 갱신해야 한다
  (docs/Userguide/GABES_User_Guide_v2.html:782-784).
- “절대 스케일, 선폭, NMOR 영점 교차를 자동 확인”한다는 문구는 실제 보장보다 강하다
  (docs/Userguide/GABES_User_Guide_v2.html:935-937).

### 8.3 실제 분석 예제

독립 examples 디렉터리는 없고, 실질적인 예제는
analysis/squeezing/resonant_hanle_squeezing_reference.py다.

- hanle_params()는 현재 scheme에 더 이상 존재하지 않는 velocity_classes를 계속 넘긴다
  (analysis/squeezing/resonant_hanle_squeezing_reference.py:370-398;
  config example:17). 무시되므로 결과는 같지만 재현성 문구를 흐린다.
- _hanle_readout()도 transmission에 3-point np.gradient를 쓴다
  (analysis/squeezing/resonant_hanle_squeezing_reference.py:401-413).
- 저장된 Markdown은 B=−0.04 µT에서 T=0.9962,
  |dT/dB|=0.002533/µT, 29 pT/√Hz를 기록한다
  (analysis/squeezing/resonant_hanle_squeezing_reference.md:22-25).
- 오늘 같은 기본 workflow를 재실행한 값은 T=0.642560,
  |dT/dB|=0.0679216/µT였다. 저장 산출물은 현행 exact-Doppler 코드의 reference가 아니다.
- measured Hanle CSV가 없다고 명시하면서도 29와 5.25 pT/√Hz를 구체적으로 제시한다
  (analysis/squeezing/resonant_hanle_squeezing_reference.md:41,47).
- tests/test_resonant_hanle_reference.py:14-132는 frequency anchor, shot-noise scaling,
  noise, CSV fit을 검증하지만 생성 보고서가 현재 model/config/hash와 일치하는지는 검증하지 않는다.

생성 시각, git commit/dirty, Magneto cache_version, config/source hash와 모델 scope를 저장하고
현재 코드로 다시 생성하기 전에는 기존 sensitivity 숫자를 실험 reference로 쓰면 안 된다.
이 provenance 비용은 O(file size)이며 OBE에 비해 무시 가능하다.

## 9. 권고 우선순위

1. **P1 trust/numerics:** 5-point 또는 local odd-fit slope, B-grid convergence status,
   physical/command zero 분리, 문헌 PASS/CHECK/FAIL 분리를 구현한다.
2. **P1 relaxation semantics:** paraffin transit exchange와 light pure dephasing을 분리하고
   독립 checklist acceptance와 회귀 test를 만든다.
3. **P1 truthful preset:** 2번 결정 뒤 EIT preset을 문헌 근거로 retune하거나 현재 crossover로
   relabel한다.
4. **P1 reproducibility:** stale Hanle/FWM 생성 Markdown과 User Guide 그림/성능 문구를
   현행 코드로 갱신하고 model/git/config/source provenance를 넣는다.
5. **P1 coefficients:** sourced gas/species/line별 elastic broadening과 pressure shift를 추가한다.
6. **P1 low-risk speed:** exact path의 L0_all/Ld_all 지연 생성을 적용한다.
7. **P2 polarization:** circular/near-circular plane-rotation hero를 gate하고 검증된 thin-medium
   Jones/Stokes readout을 추가한다.
8. **보류:** full VCC, geometry-derived relaxation, open-hyperfine/repump는 데이터와 runtime
   budget을 먼저 정한 opt-in reference로 둔다.

## 10. 검증 기록

- 다섯 scheme registry와 현재 순서 직접 확인
- built-in 5 regime exact-Maxwell 출력과 guard/fallback 재실행
- 121/201/401/801-point B 수렴 및 5-point derivative 진단
- paraffin transit 중복 제거 monkeypatch 진단
- b_offset command-zero/physical-zero 비교
- resonant Hanle example의 현재 출력 재실행
- python tests/verify_hanle_eit_eia.py 실행
- focused pytest:
  python -m pytest -q tests/test_magneto.py tests/test_magneto_doppler.py
  tests/test_kernels.py tests/test_headless_observables.py tests/test_schemes_render.py
  tests/test_resonant_hanle_reference.py tests/test_docs_consistency.py
  → **81 passed, 1 failed in 17.43 s**
- 유일한 focused failure는 기존 working tree에서 삭제된
  docs/FWM physics and analytic reconstruction/FWM_physics.tex를 찾는 docs consistency test다.
  Scheme 4 실패가 아니다.
- AGENTS.md가 요구한 root python -m pytest -q는 25.50 s 뒤 collection에서 중단됐다.
  ignored optional reference
  references/AutoOD/rb_od_calculator_optimized_final.py가 설치되지 않은 sympy를 import하기
  때문이다.

## 11. 최종 판정

Scheme 4는 **유용한 실제 원자물리**를 구현한다. Zeeman CG 구조, optical pumping,
grouped-emission TOC, transverse-field sign switching, coated-cell two-region coherence,
buffer-cell relaxation, 정확 Maxwell 평균과 Beer–Lambert/NMOR readout은 실험 조건 탐색과
초기 trace 해석에 충분히 가치 있다. 특히 2026-09-14의 pole-residue 전환은 옛 coarse
velocity-class 오차를 제거하면서 runtime도 크게 낮춘 좋은 개선이다.

하지만 현재 “EIT dip” preset은 crossover이고, default NMOR slope는 B-grid에 수렴하지 않으며,
paraffin transit rate가 exchange와 pure dephasing에 중복 적용된다. physical zero, 임의 편광의
readout 의미, open hyperfine/repump, collision provenance, absolute scale도 아직 제한된다.
생성 sensitivity 예제는 현행 코드와 크게 어긋난다.

따라서 현재 적합도는 **정성·반정량 experimental atomic-physics reference 및 설계 도구**다.
**절대 magnetometry, publication-grade linewidth/rotation/sensitivity reference로는 아직
부적합**하다. 다음 투자 대비 효과가 가장 큰 작업은 새 거대 물리 모델이 아니라
transit semantics, slope/zero convergence와 정직한 status/provenance를 먼저 고치는 것이다.
