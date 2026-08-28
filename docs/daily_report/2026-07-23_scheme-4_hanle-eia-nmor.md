# 2026-07-23 Scheme 4 물리 검토: Hanle / EIA / NMOR

## 1. 오늘의 선택과 현재 다섯 스킴 순서

Asia/Seoul 로컬 날짜의 일(day)은 23이므로

`n = (23 mod 5) + 1 = 4`

이다. 드롭다운 순서는 `gabes/schemes/__init__.py:18-24`의 `_SCHEMES`가 결정하며,
현재 정의와 출력은 다음과 같다.

| 순번 | 레지스트리 인스턴스 / 이름 | 현재 표시 제목 | 주 출력 |
|---:|---|---|---|
| 1 | `SASScheme()` / `sas` | Absorption spectroscopy (OD / SAS) | 펌프-off OD와 펌프-on 포화흡수 스펙트럼 |
| 2 | `LambdaScheme()` / `lambda` | Λ coherence (EIT / AT / CPT) | 3준위 Λ 투과·분산, EIT/AT/CPT |
| 3 | `RydbergEITScheme()` / `rydberg_eit` | Rydberg-EIT electrometry | cascade EIT와 microwave AT |
| 4 | `MagnetoScheme()` / `magneto` | Magneto-optics (Hanle/MOR) | Hanle/EIA 투과와 NMOR 회전 |
| 5 | `FWMScheme()` / `fwm` | Four-wave mixing (Squeezing / Biphoton) | seeded gain/squeezing과 SFWM biphoton 추정 |

README의 사용자 관점 순서도 이 레지스트리와 일치한다 (`README.md:8-17`). 따라서 오늘의 대상은
`MagnetoScheme`이다 (`gabes/schemes/magneto.py:127-173`).

## 2. 선행 제안, 문서, 테스트, 예제, issue 노트 선검색

먼저 새 제안을 만들기 전에 다음을 조사했다.

- 현재 구현: `gabes/schemes/magneto.py:1-746`
- Zeeman/CG/TOC 기반 원자 모델: `gabes/zeeman.py:17-140`, `gabes/atoms.py:21-98`
- 수치 커널과 공통 관측량: `gabes/kernels.py:276-410`, `gabes/core.py:96-121`,
  `gabes/observables.py`
- 사용자 문서: `README.md:8-17`, `docs/Userguide/GABES_User_Guide_v2.html:538-541`,
  `610-619`, `799-805`
- 현재 물리 테스트: `tests/test_magneto.py:26-256`
- 문헌 대조 스크립트: `tests/verify_hanle_eit_eia.py:1-277`
- 실제 분석 예제: `analysis/squeezing/resonant_hanle_squeezing_reference.py:352-412`,
  `527-625`, `657-725`, `900-990`; 예제 설정은
  `analysis/squeezing/resonant_hanle_experiment_config.example.json`
- 계획/TODO: `docs/checklist.json:22-26`, `60-64`, `102-106` 및
  `gabes/constants.py:79-89`
- 이전 Scheme 4 보고서: 2026-06-23, 06-28, 07-03, 07-13의
  `docs/daily_report/*scheme-4_hanle-eia-nmor.md`

별도 `examples/` 디렉터리는 없지만, `analysis/squeezing/resonant_hanle_squeezing_reference.py`가
FWM probe를 87Rb D1 Hanle cell에 넣고 balanced detection·절대 감도·측정 CSV 보정까지 수행하는
실행 가능한 실험 예제다. 이 예제 자체도 현재 모델을 “compact and semi-quantitative”라고 명시하고
측정 trace 보정을 권한다 (`analysis/squeezing/resonant_hanle_squeezing_reference.py:1147-1150`).

로컬에는 별도 issue-note 파일이 없었다. 2026-07-23에 원격
[GitHub Issues](https://github.com/Shake2313/fwm-squeezing-app/issues?q=is%3Aissue)도 확인했으며
open/closed 모두 0건이다. 대신 아래와 같은 기존 개선안이 이미 명시돼 있었다.

| 기존 항목 | 상태 | 계산 비용 판단 | 오늘의 판단 |
|---|---|---|---|
| `buffer-gas-pressure-shift` | deferred | pressure shift는 scalar detuning 보정, 저차 Dicke 보정도 동일 차원이라 거의 0 | 여전히 물리 가치/비용비가 좋다. gas·species·line별 실험계수 근거가 필요하다. |
| `magneto-buffer-relaxation-map` | deferred | solve 전에 기본 scalar rate를 정하므로 사실상 0 | Ne 압력·온도·beam geometry에서 diffusion/relaxation 기본값을 만들되 measured override를 보존해야 한다. |
| `full-velocity-changing-collision-kernel` | deferred GROUP C | 속도군을 서로 결합해 현재 separable `(B, v)` solve를 깨므로 큼 | 정밀 buffer-cell fitting에는 필요할 수 있으나 interactive 기본 경로에는 무겁다. |
| `figureless-observables-paths` | done | figure 생성을 건너뜀 | 현재 headless 관측량은 대표 측정에서 약 0.18 ms로 충분히 가볍다 (`gabes/schemes/magneto.py:620-746`). |

이전 보고서가 제안했던 `b_offset_ut`는 이미 구현됐다 (`gabes/schemes/magneto.py:220-225`,
`431-435`). `zeeman_manifold()` template도 LRU cache를 사용한다 (`gabes/zeeman.py:69-70`).

## 3. 현재 구현이 담고 있는 실제 원자물리

### 실험적으로 유용한 핵심

이 스킴은 단순 Lorentzian toy curve보다 훨씬 실제 원자물리에 가깝다.

1. 선택한 87Rb D1 `Fg -> Fe`의 모든 `mF` 상태를 만들고, `q=-1,0,+1` 편광별 CG 결합을
   계산한다 (`gabes/zeeman.py:89-118`). 네 개 D1 hyperfine 전이를 모두 고를 수 있다
   (`gabes/schemes/magneto.py:47-54`, `207-210`).
2. 자발방출을 전이 하나씩 독립 jump로 쪼개지 않고 편광별 `Sigma_q`로 묶어 excited coherence가
   ground coherence로 전달되는 transfer of coherence(TOC)를 보존한다
   (`gabes/zeeman.py:96-118`). 이 덕분에 `Fe=Fg+1`의 intrinsic EIA와 `Fe<=Fg`의 EIT 부호를
   테스트할 수 있다 (`tests/test_magneto.py:179-192`).
3. QWP 각도를 실제 `sigma+/-` 복소 구동 진폭으로 바꾸고 (`gabes/schemes/magneto.py:113-124`),
   종방향 scan field와 잔류 횡장을 같은 Zeeman Hamiltonian에 넣는다
   (`gabes/schemes/magneto.py:431-479`, `511-530`). 이 구조가 linear CPT/Hanle dip,
   circular MIA/EIA, transverse-field LCA 전환을 만든다 (`tests/test_magneto.py:123-138`,
   `195-211`).
4. 파라핀 셀은 illuminated/dark 두 영역의 밀도행렬을 교환해 wall-preserved coherence와 Ramsey
   narrowing을 표현한다 (`gabes/schemes/magneto.py:419-429`, `481-488`, `532-576`). Buffer 셀은
   단일 영역에 ground relaxation과 collisional depolarization을 넣는다
   (`gabes/schemes/magneto.py:411-418`, `578-601`).
5. 투과는 Beer-Lambert `exp(-alpha L)`로, NMOR은
   `theta = kL Re(chi_+ - chi_-)/4`로 계산한다 (`gabes/schemes/magneto.py:642-653`). 즉 NMOR은
   Hanle 투과의 이름만 바꾼 readout이 아니라 편광 분산차에서 나온 실제 회전 관측량이다.
6. 최근 분석 예제는 측정한 Hanle CSV에 대해 B scale/offset과 detector affine scale을 fitting하고
   (`analysis/squeezing/resonant_hanle_squeezing_reference.py:536-625`), 광검출 shot/electronic noise를
   field sensitivity로 변환한다 (`analysis/squeezing/resonant_hanle_squeezing_reference.py:657-725`).
   실험 trace가 있다면 compact model을 현장 전달함수로 보정하는 실용적인 경로다.

문헌상 87Rb `F=2 -> F'=1` 파라핀 셀에서 linear polarization은 0.12 mG CPT, circular
polarization은 residual transverse field와 coherence preservation에 의한 0.20 mG absorption으로
전환된다. 저장소가 인용한 [Lee & Moon, JOSA B 30, 2301 (2013)](https://doi.org/10.1364/JOSAB.30.002301)의
실험 대상과 현재 manifold/QWP/two-region 구조는 정성적으로 잘 맞는다.

다만 현재 `tests/verify_hanle_eit_eia.py`의 실제 출력은 linear 1.195 mG, circular 1.739 mG다.
문헌값보다 각각 약 10.0배, 8.7배 넓고 첫 값은 sub-mG도 아니다. 따라서 스크립트 마지막의
“same sub-mG order” 문구 (`tests/verify_hanle_eit_eia.py:269-271`)는 수치 판정으로는 틀리며,
이전 Scheme 4 보고서가 지적한 exact-linewidth 불일치가 계속 남아 있다.

### 실험 레퍼런스로서의 등급

현재 코드는 다음 목적에는 유용하다.

- 전이, 편광, 횡장, wall lifetime, ground relaxation을 바꿨을 때 dip/peak/zero-crossing이 어떻게
  바뀌는지 확인하는 실험 설계와 sanity check
- Hanle/NMOR 동작점과 parameter sensitivity 탐색
- 측정 Hanle trace를 넣기 전의 compact forward model 및 보정 초기값

반면 절대 linewidth·절대 NMOR slope·절대 자력계 감도의 독립 표준으로는 아직 부족하다. 한 번에
하나의 `Fg -> Fe` manifold만 풀고 다른 hyperfine manifold, repump, 실제 빔 단면과 diffusion,
wall collision 분포, magnetic-field inhomogeneity를 포함하지 않는다 (`gabes/schemes/magneto.py:390-429`).
파라핀의 네 rate knob와 buffer ground-rate knob도 측정에서 직접 정해지는 effective parameter다
(`gabes/schemes/magneto.py:241-273`). 아래의 새 수치 진단까지 고려하면 **정성적/반정량적 실험
참조에는 좋지만, 측정 trace 또는 별도 수렴 검증 없이 논문 수치나 pT/sqrt(Hz) 성능을 인용하면 안
되는 수준**이다.

## 4. 오늘 확인한 새 물리·수치 신뢰 경계

### 4.1 Ne pressure broadening이 자발방출/TOC rate로 들어간다

현재 buffer 경로는

`gamma_opt = GAMMA_D1 + buffer_gamma`

를 만든 뒤 이를 `zeeman_manifold(..., gamma=gamma_opt)`에 넘긴다
(`gabes/schemes/magneto.py:407-415`). 그런데 `zeeman_manifold()`의 `gamma`는 편광별 방출 jump의
진폭 `sqrt(gamma)`가 되어 population decay와 TOC refeeding을 동시에 정한다
(`gabes/zeeman.py:96-118`). 따라서 Ne의 homogeneous pressure broadening이 단순 optical
dephasing이 아니라 들뜬상태 수명 단축 및 추가 자발방출/TOC처럼 작동한다.

저장소의 다른 원자 모델은 자연방출 `Gamma`와 추가 collisional optical dephasing을 분리하고,
추가 FWHM의 절반을 ground-excited coherence에 넣는 올바른 구조를 이미 설명한다
(`gabes/atoms.py:101-112`, `124-129`). 같은 원칙으로 진단용 대체 dissipator를 만들어 비교했다.
이는 완전한 buffer collision 모델이 아니라, **자연방출과 elastic optical dephasing을 분리했을 때의
민감도 진단**이다.

| 조건 | 현재 구현 | 자연방출 + optical pure-dephasing 진단 | 변화 |
|---|---:|---:|---:|
| 기본 Buffer Hanle, 20 Torr, 0.8 mW/cm²: `T(B=0)` | 0.75650 | 0.75742 | +0.12% |
| 같은 조건: absorption feature amplitude | -4.046 | -4.006 | 0.99% 작음 |
| 허용 전이 `F=1 -> F'=2`, 5 Torr, 20 mW/cm², 25 °C, 5 mm: `T(B=0)` | 0.48097 | 0.65145 | +35.4% |
| 같은 조건: absorption feature amplitude [m⁻¹] | 1.3274 | 0.6337 | 52.3% 작음 |
| 같은 조건: central width [µT] | 106.57 | 99.97 | 6.2% 작음 |

기본 약한 probe preset에서는 우연히 차이가 작지만, UI가 허용하는 강한 구동과 intrinsic-EIA
전이에서는 배경 투과와 대비가 크게 달라진다. 이 수정은 dissipator의 항만 바꾸고 Hilbert/Liouville
차원과 `(B, v)` grid를 그대로 유지하므로 solve overhead는 사실상 없다. 먼저 `gamma_natural`과
`gamma_optical_dephasing`을 API에서 분리하고, 0 Torr에서 byte/수치 동등성, pressure-on에서 excited
population lifetime과 TOC가 자연 `Gamma`에 고정되는 테스트를 추가하는 것이 가장 값싼 물리
정확도 개선이다. 단, 실제 Ne collision이 만드는 excited-state depolarization은 별도 실험계수로
추가해야 하며 optical broadening에 묵시적으로 섞어서는 안 된다.

### 4.2 9개 속도군 + scalar dilution은 비선형 Hanle Doppler 평균을 수렴시키지 못한다

현재는 `[-3 sigma_v, 3 sigma_v]` 균일 grid를 기본 9점으로 적분하고
(`gabes/schemes/magneto.py:68-77`, `291-298`), 같은 grid의 scalar Lorentzian 오차를
fine-grid Voigt 값으로 보정한 `doppler_scale`을 전체 susceptibility에 곱한다
(`gabes/schemes/magneto.py:80-100`, `460-467`, `498-507`). 이 방식은 선형 line-centre OD 크기를
빠르게 맞출 수 있지만, velocity마다 달라지는 optical pumping/Zeeman coherence의 부호와 모양은
scalar 하나로 복구할 수 없다.

기본 EIT-dip preset을 같은 B grid에서 속도군만 늘려 비교한 결과는 다음과 같다.

| velocity classes | `doppler_scale` | absorption feature amplitude [m⁻¹] | central width [µT] |
|---:|---:|---:|---:|
| 9 (UI 기본) | 0.0560 | -0.02005 (dip) | 0.1666 |
| 81 | 0.5092 | -0.01333 (dip) | 0.1257 |
| 321 | 0.9760 | +0.003545 (crossover/peak) | 1.0834 |
| 641 | 0.9972 | +0.004031 (crossover/peak) | 1.0253 |

즉 기본 “EIT dip”의 부호가 촘촘한 Maxwell 적분에서 뒤집힌다. 기본 EIA preset은 부호는 유지하지만
9→641점에서 amplitude가 0.07268→0.18322 m⁻¹, width가 0.2232→0.3634 µT로 변했다. 현재 단위
테스트는 빠른 `doppler=off`, 한 속도군 조건으로 dip/peak를 검사한다
(`tests/test_magneto.py:69-72`, `123-223`). 문헌 검증 스크립트도 5개 속도군을 사용한다
(`tests/verify_hanle_eit_eia.py:71-84`, `98-107`). 따라서 그 스크립트의 “MATCH”는 연속 Maxwell
평균 검증이 아니다.

이 문제의 완전한 수정은 공짜가 아니다. 121개 B점에서 641개 균일 속도군은 대표 실행에서 약
9-14 s가 걸려, 기본 9개 경로의 수백 ms보다 수십 배 무겁다. 다음 두 단계가 합리적이다.

1. **무추가-solve 안전장치:** 이미 계산한 `doppler_scale`, velocity spacing, optical width로
   “coarse Doppler quadrature / qualitative only” 상태를 표시하고, 논문/분석 예제의 hero linewidth와
   feature sign을 수렴값처럼 내보내지 않는다. 비용은 사실상 0이다.
2. **reference fidelity 경로:** 공명 속도 `v_res=Delta/k` 근처를 조밀하게 두는 composite/adaptive
   quadrature와 결과 tolerance 검사를 사용한다. 진단용 37-153점 비균일 grid는 641점보다 빨랐지만
   amplitude 수렴이 느렸으므로, 고정된 적은 점수를 물리 보존으로 단정하지 말고 자동 refinement와
   sign/width convergence를 검사해야 한다. 이 경로는 interactive 기본보다 비싸지만 full VCC처럼
   속도군끼리 결합하지는 않는다.

### 4.3 기본 B grid의 NMOR 영점 미분은 아직 수렴하지 않는다

NMOR hero slope는 중앙 표본에서 `np.gradient(rotation, x)`로 계산한다
(`gabes/schemes/magneto.py:671-680`). 기본 ±2 µT/121점, 9 velocity-class 조건에서 grid만 바꾸면:

| scan points | ΔB [µT] | 중앙차분 slope [mrad/µT] |
|---:|---:|---:|
| 51 | 0.0800 | -1.2606 |
| 121 (기본) | 0.0333 | -2.5525 |
| 201 | 0.0200 | -2.9511 |
| 401 | 0.0100 | -3.1503 |

기본값은 401점 값보다 약 19% 작다. 같은 기존 121점의 5-point derivative는 -2.8871 mrad/µT로
개선되지만, 401점 5-point 값 -3.2166 mrad/µT보다 여전히 약 10% 작다. 따라서 5-point/local odd
polynomial derivative와 `samples per central width` 상태 표시는 추가 solve 없이 정확도를 높일 수
있고, 정량 magnetometry 모드에서는 최소 201점 또는 자동 B-grid refinement가 필요하다. 이 개선은
full physics 변경이 아니라 관측량 추출 정확도 개선이다.

## 5. 계산 비용과 순수 코드 최적화

warm default paraffin 조건의 한 대표 profile은 총 0.182 s 중 0.168 s가 two-region steady-state
경로였고 (`gabes/schemes/magneto.py:532-576`, `gabes/kernels.py:379-410`), headless observables는
약 0.18 ms였다. 따라서 Matplotlib나 작은 Python helper보다 solve 결과의 메모리 흐름을 줄이는 것이
우선이다.

1. **커널 안에서 Doppler-weighted coherence를 직접 축약한다.** 현재 Numba 커널은 모든
   `(B,v)`의 full `rho`를 복원해 Python으로 돌려주고 (`gabes/kernels.py:399-410`), 이후 Python이
   `chi_+`, `chi_-`, `chi_probe`를 순회·평균한다 (`gabes/schemes/magneto.py:493-496`, `604-618`).
   커널에서 필요한 coherence 세 개만 weighted sum하면 결과를 바꾸지 않고 full `rho`
   materialization을 없앨 수 있다. 특히 121×641×8×8 complex `rho`는 약 76 MiB이므로 reference
   Doppler 경로에서 가치가 커진다. 기존 NumPy fallback/Numba parity test를 유지해야 한다.
2. **복소 B-stack을 만든 뒤 basis-transform하지 말고 affine real 계수를 전달한다.** 현재
   `C_xy + B*C_z` 복소 stack을 만들고 (`gabes/schemes/magneto.py:469-485`), 커널 wrapper가 다시
   `U^dagger L U`를 batched 수행한다 (`gabes/kernels.py:399-405`). `C_xy`, `C_z`, dark coefficient를
   한 번만 real Hermitian basis로 바꾸고 커널에서 B-affine 조립하면 동일 행렬을 더 적은 메모리로
   만든다. 물리는 변하지 않지만 구현·동등성 검증 비용은 1번보다 크다.
3. **작은 immutable helper만 제한적으로 cache한다.** `_hamiltonian()`은 호출마다
   `angular_momentum_matrices()`를 다시 만든다 (`gabes/schemes/magneto.py:511-520`,
   `gabes/zeeman.py:48-66`). 한 호출은 약 22 µs라 전체 solve에 비하면 작지만, read-only LRU는
   반복 sweep에서 안전한 미세 최적화다. `zeeman_manifold()` 자체는 이미 cache돼 있으므로 중복
   template cache는 필요 없다.

`threadpoolctl` controller 재사용도 미세 benchmark에서는 context 생성비를 줄였지만, 전체
Magneto compute 20회 비교에서는 변동폭보다 유의한 개선이 없었다. 현재 우선순위로 권하지 않는다.

## 6. 권고 우선순위와 최종 판단

1. **P0 물리 의미 수정:** Ne optical broadening을 natural spontaneous emission/TOC와 분리하고
   optical pure dephasing으로 넣는다. 행렬 차원과 solve 횟수는 그대로다.
2. **P0 신뢰 경계:** Doppler quadrature 수렴 상태를 hero/readout에 노출하고, 현재 5/9-class
   문헌 대조를 정량 MATCH로 취급하지 않는다. Reference fidelity에는 adaptive refinement를 둔다.
3. **P1 관측량 정확도:** NMOR slope를 고차 중앙미분/local fit으로 바꾸고 B samples-per-width 상태를
   표시한다. 비용은 거의 0이며, 정량 모드만 더 촘촘한 B grid를 쓴다.
4. **P1 기존 저부하 물리:** pressure shift, 저차 Dicke correction, pressure-to-relaxation mapping을
   실험계수 근거와 override를 유지한 채 구현한다.
5. **P2 성능:** 커널 내부 coherence contraction을 먼저 하고, 그 다음 affine real-basis 조립을
   검토한다. 이는 reference velocity refinement의 메모리 비용도 낮춘다.

`MagnetoScheme`은 Zeeman manifold, CG, TOC, QWP 편광, 횡장, two-region Ramsey exchange,
`chi_+-chi_-` NMOR을 한 엔진에 모은 실질적인 원자물리 모델이다. 실험가가 현상 부호와 knob
민감도를 이해하고 measured-trace fit의 초기 모델로 쓰기에는 매우 유용하다. 그러나 오늘의
수렴 검사는 기본 Doppler grid가 paraffin “EIT dip”의 부호까지 바꿀 수 있고, 기본 NMOR slope도
grid에 약 20% 민감함을 보였다. 그러므로 현재 버전은 **정성적/반정량적 실험 참조**로 평가하며,
절대 linewidth·대비·자력계 감도 레퍼런스로 승격하려면 dissipator 의미 분리, Doppler/B-grid 수렴
상태, 측정 trace calibration이 먼저 필요하다.

## 7. 검증 기록

보고서 작성 전 수치 진단은 production 파일을 수정하지 않고 runtime monkeypatch/독립 계산으로
수행했다. 최종 테스트 결과는 아래에 기록한다.

- `python tests/verify_hanle_eit_eia.py`: 51.6 s에 완료. 편광 부호 전환은 재현했으나 계산 폭
  1.195/1.739 mG 대 문헌 0.12/0.20 mG의 불일치를 재확인했다.
- `python -m pytest tests/test_magneto.py tests/test_kernels.py tests/test_headless_observables.py tests/test_schemes_render.py tests/test_resonant_hanle_reference.py -q`:
  **63 passed in 24.04 s**.
- AGENTS.md 필수 전체 회귀 `python -m pytest -q`: **189 passed in 59.68 s**.
