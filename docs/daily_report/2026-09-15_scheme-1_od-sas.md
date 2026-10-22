# GABES 일일 Scheme 물리 검토 — 2026-09-15

## 1. 선택·판정

- 현지 날짜: 2026-09-15, Asia/Seoul. `day = 15`.
- \(n=(15\bmod5)+1=1\).
- 대상: **Scheme 1 — Absorption spectroscopy (OD / SAS)**, registry key `sas`.
- 검토 HEAD: `51d5b74`. 기존 staged/unstaged 작업 다수 있음. 이번 실행은 production
  코드·기존 문서·테스트를 수정하지 않고 이 보고서만 추가함.

**판정:** 실제 원자·열·전파 물리 있음. **반정량 실험 계획·초기 해석 reference로
유용.** 특히 pump-off Rb 약한-probe 흡수, hyperfine 위치·상대 세기, Doppler
envelope, velocity-correlated hyperfine pumping, 약한-probe 분산, 유한-probe
Gaussian 전파 유용. 그러나 **절대 pumped-SAS, power-to-Rabi, buffer/coated-cell,
고온 Cs, 편광/자기장, 정량 선폭 fitting reference 아님.**

2026-09-10 Scheme 1 검토 뒤 scheme source/test/example 변화 없음. 다만 2026-09-14
checklist에 dense velocity arithmetic 성능 항목 새로 들어옴. 이전 물리 판정 유지,
성능 권고만 더 구체화됨.

## 2. 현재 다섯 scheme 정의·순서

드롭다운 순서의 단일 근거: [`_SCHEMES`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/schemes/__init__.py:19>).

| n | registry key | 표시명 | 정의 |
|---:|---|---|---|
| 1 | `sas` | Absorption spectroscopy (OD / SAS) | [`SASScheme`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/schemes/sas.py:175>) |
| 2 | `lambda` | Λ coherence (EIT / AT / CPT) | [`LambdaScheme`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/schemes/absorption.py:466>) |
| 3 | `rydberg_eit` | Rydberg-EIT electrometry | [`RydbergEITScheme`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/schemes/rydberg.py:94>) |
| 4 | `magneto` | Magneto-optics (Hanle/MOR) | [`MagnetoScheme`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/schemes/magneto.py:121>) |
| 5 | `fwm` | Four-wave mixing (Squeezing / Biphoton) | [`FWMScheme`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/schemes/fwm.py:3547>) |

`ODScheme`는 단일 two-level/85Rb-D1 내부 검증 primitive. UI 미등록. 실제 OD는
Scheme 1의 `pump_power_mw = 0` 극한
([`__init__.py:9-24`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/schemes/__init__.py:9>)).

## 3. 기존 개선 제안 먼저 검색

확인 범위:

- Scheme 1 선행 일일 보고서 10건: 2026-06-25, 07-10, 07-20, 07-30,
  08-05, 08-10, 08-25, 08-30, 09-05, 09-10. 최신:
  [`2026-09-10_scheme-1_od-sas.md`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/docs/daily_report/2026-09-10_scheme-1_od-sas.md:1>).
- 통합 작업 등록부:
  [`docs/checklist.json`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/docs/checklist.json:1858>).
- 직접 코드 TODO: Ne scalar를 gas/species/line coefficient table과 pressure shift로
  승격
  ([`constants.py:79-89`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/constants.py:79>)).
- 별도 로컬 `issue`/`TODO`/`proposal` 파일: 없음.
- 공개 GitHub issue: 2026-09-15 API 조회 결과 open/closed 모두 0.

따라서 **기존 제안 없는 상태 아님.** 이미 충분히 구체적. 새 physics wishlist보다
기존 P1 신뢰 계약과 안전한 성능 이식 먼저.

### 현재 관련 항목

| 항목 | 상태 | 판단 |
|---|---|---|
| collision coefficient provenance + pressure shift | P1 / large / `ready` | 여전히 최우선 physics 교정 |
| species/line Rabi + saturation provenance | P1 / medium / `ready` | pump-power 축 정량화에 필수 |
| sub-Doppler readout trust | P1 / small / `done` | gate 정상. 기본 grid 자체는 미수렴 |
| paraffin population memory | P2 / medium / `done` | opt-in 반정량 모델로 적절 |
| dense scan×velocity arithmetic 제거 | P3 / medium / `ready` | 09-14 신규. 가장 직접적인 속도 과제 |
| experimental fitting | P2 / large / `ready` | 명시적 실행·bounds·holdout 필요 |
| isotope fraction | P2 / medium / `needs_decision` | enriched-cell 수요 결정 먼저 |
| segmented cell | P2 / large / `ready` | pump depletion에 유용, 비쌈 |
| weak-probe Zeeman / full VCC | P2 / research / `parked` | 측정 target·dataset·runtime budget 먼저 |
| result export/provenance | P2 / medium / `in_progress` | round-trip·named arrays·extra views 남음 |

근거: [`collision`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/docs/checklist.json:1858>),
[`Rabi`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/docs/checklist.json:1898>),
[`dense arithmetic`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/docs/checklist.json:2548>),
[`fitting`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/docs/checklist.json:3160>),
[`segmented cell`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/docs/checklist.json:3260>),
[`export`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/docs/checklist.json:3287>).

## 4. 구현된 실제 물리

### 4.1 원자 데이터·절대 약한-probe 흡수

- 85Rb, 87Rb, 133Cs D1/D2. hyperfine A/B, centroid, 질량, 자연선폭,
  abundance 명시
  ([`species.py:130-184`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/species.py:130>)).
- Casimir hyperfine energy와 Wigner-6j/3j strength. manifold는 allowed Fg→Fe,
  CG-branched decay, thermal transit reset 구성
  ([`species.py:94-118`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/species.py:94>),
  [`337-416`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/species.py:337>)).
- 전이 적분 흡수:

  \[
  \int\alpha_t\,d\delta=
  \frac{\pi k C_{F,t}^{2}|d|^{2}N}
       {\epsilon_0\hbar\,2(2I+1)}.
  \]

  `C_F²`가 mF 합을 이미 포함. thermal sublevel population 한 번만 계수
  ([`species.py:248-269`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/species.py:248>)).
  85Rb D1 두 ground manifold 적분비 7/5. 전체 D1/D2 oscillator sum과 축퇴도 비를
  3 isotope 모두 독립 항등식으로 검사
  ([`test_sas.py:91-109`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/tests/test_sas.py:91>)).
- density temperature와 Doppler temperature 분리. tie 기본 on이라 legacy/default
  array exact 유지
  ([`sas.py:56-75`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/schemes/sas.py:56>),
  [`test_sas.py:698-739`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/tests/test_sas.py:698>)).

이 부분: Rb pump-off weak-probe scale·shape reference로 가장 강함. 단 “lab calculator와
software match”이지 held-out 실험 검증 아님.

### 4.2 SAS: 중요한 velocity correlation·hyperfine pumping

Pump atom-frame detuning `scan + kv`, counter-propagating probe `scan - kv`.
Pump가 만든 population difference와 probe Lorentzian을 **같은 velocity class**에서 곱한
뒤 Maxwell 평균
([`sas.py:590-631`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/schemes/sas.py:590>)).
두 평균을 독립 factorization하지 않음. Lamb dip/crossover에 필요한 핵심 물리.

CG-branched spontaneous decay + transit reload가 hyperfine optical pumping 생성.
Transit 감소 시 crossover 증강 테스트
([`test_sas.py:339-354`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/tests/test_sas.py:339>)).
단순 two-level hole-burning 그림보다 실험적 의미 큼.

### 4.3 분산·유한 probe·실험 CSV

- 약한 probe: 각 전이가 하나의 complex profile 제공. Im χ는 흡수, Re χ는 같은
  population·Doppler measure의 Kramers-Kronig partner
  ([`sas.py:539-553`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/schemes/sas.py:539>)).
  독립 Faddeeva Voigt와 0.1% 안 검사
  ([`test_sas.py:594-605`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/tests/test_sas.py:594>)).
- 유한 probe: Gaussian peak intensity, line-specific closed-two-level `I_sat`, power
  broadening, `dI/dz=-alpha(I)I`, 32점 radial quadrature
  ([`probe_saturation.py:34-82`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/probe_saturation.py:34>),
  [`166-260`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/probe_saturation.py:166>)).
- Rb D2만 supplied software-reference scope. D1/Cs는 generalized qualitative.
  finite pump+probe는 frozen pump population 위의 순차 근사. 음의 흡수 발생 시
  passivity limiter 보고
  ([`sas.py:113-168`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/schemes/sas.py:113>),
  [`466-501`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/schemes/sas.py:466>)).
- CSV: relative extrema normalization과 evidence-backed absolute calibration 분리.
  acquisition-order forward/reverse branch 보존
  ([`experimental_csv.py:214-314`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/experimental_csv.py:214>),
  [`443-497`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/experimental_csv.py:443>),
  [`743-881`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/experimental_csv.py:743>)).

## 5. 실험 reference로서 핵심 한계

### 5.1 `Optical density` 이름과 log convention 불일치

현재 공통 함수:

```text
Optical density = -log10(T) = alpha*L/ln(10)
```

([`observables.py:545-547`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/observables.py:545>)).
반면 supplied AutoOD 기술 문서: `OD=-ln(T)=alpha*L`
([`AutoOD_Natural_Rb_D2_Model.tex:403-405`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/references/AutoOD-Rb/docs/AutoOD_Natural_Rb_D2_Model.tex:403>)).

오늘 직접 계산, 85Rb D1, 50 °C, 12.5 mm:

| 같은 spectrum | peak |
|---|---:|
| natural optical depth `tau=-ln(T)` | `1.2556477622` |
| 현재 GABES `-log10(T)` | `0.5453208943` |
| ratio | `2.302585093 = ln(10)` |

Raw α/T 오류 아님. 명명·export에서 정확히 2.3026배 혼동 가능. **P1 trust issue.**
`Optical depth tau=-ln T`와 `Decadic absorbance A10=-log10 T` 동시 제공 권고.
추가 solve 없음, O(Nscan) scale만 필요.

### 5.2 pump power는 아직 정밀 Rabi calibration 아님

`pump_rabi_from_power()`가 모든 isotope/line에 하나의 Rb85-D1
`I_SAT=4.484 mW/cm²` 사용
([`species.py:272-282`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/species.py:272>),
[`constants.py:22-31`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/constants.py:22>)).

0.5 mW, 1 mm waist 직접 확인: Rb85/87, Cs, D1/D2 모두
`Omega/Gamma = 1.8839842329`. line/species/polarization 차이 사라짐. 따라서 power
scan은 trend knob. 절대 pump intensity 비교에는 부적합. `rabi_rel`도 lumped
`sqrt((2Fg+1)S/Tmax)` prescription이며 individual mF matrix element 아님
([`species.py:365-407`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/species.py:365>)).

### 5.3 collisional physics가 pump/probe에 일관되지 않음

- `self_broadened_gamma(iso,N)`가 `iso`를 받지만 single Rb coefficient만 반환
  ([`species.py:210-212`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/species.py:210>)).
- Cs fallback self width 직접 확인: 100/150/200 °C에서
  `1.0186 / 15.1751 / 124.0098 MHz`. 고온에서 지배적. 검증되지 않은 외삽.
- Ne는 `3.91 MHz/Torr` broadening-only. pressure shift·temperature dependence·Dicke/VCC
  없음
  ([`constants.py:79-89`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/constants.py:79>)).
- `gamma_eff`는 probe line profile·pump table spacing에 들어감. 그러나 species path의
  pump Liouvillian은 natural decay/transit manifold로 먼저 만들어져 self/buffer elastic
  optical dephasing 미포함
  ([`sas.py:393-413`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/schemes/sas.py:393>),
  [`568-585`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/schemes/sas.py:568>)).

결과: probe width와 pumped population dynamics가 서로 다른 homogeneous-width 계약 사용.

### 5.4 기본 선폭·공간·상태 모델 한계

- 기본 1401점 sub-Doppler FWHM은 2.34 samples뿐. UI가
  `resolution-limited`로 표시하고 Lock Slope를 숨기는 것은 좋음
  ([`sas.py:804-835`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/schemes/sas.py:804>),
  [`test_sas.py:406-423`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/tests/test_sas.py:406>)).
  하지만 숫자 자체는 grid-converged 아님.
- pump Rabi는 entrance power 한 번. z-direction pump depletion·temperature/field gradient
  없음. uniform slab.
- finite pump+probe는 joint coherent solve 아님. probe-induced hyperfine repumping 없음.
- hyperfine-lumped scalar model. mF, input polarization, Zeeman splitting,
  Faraday/Voigt dichroism·rotation 없음.
- paraffin: fixed `T1=25.1 ms`, one-zone return cadence. 한 87Rb cell reference를 모든
  cell에 쓰는 반정량 population-memory MVP
  ([`sas.py:44`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/schemes/sas.py:44>),
  [`1090-1113`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/schemes/sas.py:1090>)).
- normalization 변경 뒤에도 `cache_version = "7"`
  ([`sas.py:179`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/schemes/sas.py:179>)).
  앱 cache/export provenance 계약과 불일치 위험. O(1) version bump 필요.

## 6. 오늘의 수치 spot check

동일 머신, warm median 중심. `OD10`은 현재 UI base-10 값. `tau`는 natural optical
depth.

| 확인 | 결과 | 해석 |
|---|---:|---|
| 기본 SAS compute, 1401점 | `0.613 s` median | natural Rb D1, 0.5 mW, 40 °C, 75 mm |
| Tmin / peak tau / peak OD10 | `0.11805 / 2.13665 / 0.92794` | log convention 필수 |
| 기본 sub-Doppler FWHM | `18.8684 MHz`, `2.342 samples/FWHM` | `resolution-limited` |
| 4001점 | `17.6561 MHz`, `6.261 samples/FWHM`, `1.61 s` | gate는 `resolved`; 11201점 대비 2.0% 큼 |
| 5601 / 11201점 | `17.3596 / 17.3108 MHz`, `2.17 / 4.27 s` | 마지막 둘 0.282% 차 |
| Gaussian Doppler FWHM | `518.694 MHz` | 단일 optical line 계산 폭 |
| peak phase | `688.263 mrad` | single-pass phase, group delay 아님 |
| pump-off 1401점 | `0.299 s`, Tmin `0.10510` | pump population solve 생략 |
| velocity classes | Rb85 `1403`, Rb87 `1387` | 1401×Nv float64 배열 하나가 약 15 MiB/isotope |

기본 선폭은 11201점보다 9.0% 큼. `resolved` gate는 최소 sampling gate이지 convergence
certificate 아님. 정량 linewidth/lock export에는 fixed-window refinement와 수렴값 필요.

## 7. 문서·테스트·예제 평가

### 7.1 문서

README는 0–100 µW finite probe, Rb-D2 reference tier, D1/Cs generalized tier,
population-once, 7/5, two temperatures, dispersion을 현재 코드에 가깝게 설명
([`README.md:12`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/README.md:12>),
[`314-386`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/README.md:314>)).

남은 문제:

- README와 in-app info의 `validated`, `AutoOD ... <0.1%` 표현은 software/reference
  reproduction을 held-out experiment validation처럼 읽히게 함
  ([`README.md:12`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/README.md:12>),
  [`sas.py:296-306`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/schemes/sas.py:296>)).
- HTML User Guide는 “finite probe는 natural Rb D2에서만”이라고 씀. 실제 UI는 모든
  built-in species/line에 제공하고 tier로 제한을 설명
  ([`GABES_User_Guide_v2.html:551-584`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/docs/Userguide/GABES_User_Guide_v2.html:551>),
  [`sas.py:188-213`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/schemes/sas.py:188>)).
- Guide는 probe 기본값을 `0 mW`로 씀. 실제 UI `0 µW`.
- Guide Scheme 1 주요 출력 표에 weak-probe Dispersion/phase 없음. split temperature와
  Scheme-1 paraffin memory 설명도 없음
  ([`GABES_User_Guide_v2.html:527-530`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/docs/Userguide/GABES_User_Guide_v2.html:527>)).
- OD 자연로그/base-10 구분 없음.

### 7.2 예제

Scheme 1 독립 runnable script/config 없음. 실제 예제는 User Guide의
[`od.png`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/docs/Userguide/userguide_assets/od.png>)와
[`sas.png`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/docs/Userguide/userguide_assets/sas.png>)뿐.
둘 다 2026-06-08 생성.

시각 확인 결과:

- Transmission과 base-10 `Optical density`를 한 figure에 수직 배치. 현재 앱의
  Transmission/OD/Dispersion carousel과 다름.
- OD 그림 Tmin 약 0.27. 오늘 기본값 Tmin 0.118. thermal normalization 교정 전 숫자.
- finite-probe qualifier, resolution status, two-temperature/paraffin provenance,
  dispersion 없음.

개념 그림으로만 유효. current numerical example/figure regression 아님. current params,
`cache_version`, git revision, OD convention, resolution status를 담은 재생성 예제 필요.

### 7.3 테스트

오늘 실행:

- `python -m pytest -q tests/test_sas.py tests/test_absorption.py tests/test_probe_saturation.py tests/test_experimental_csv.py`
  → **132 passed in 65.64 s**.
- `python -m pytest -q references/AutoOD-NatRbD2/test_beam_saturation.py`
  → **29 passed in 6.82 s**.
- render/headless/style/docs 묶음 → **30 passed, 1 failed**. 실패는 Scheme 1 아님.
  기존 삭제 상태인 `docs/FWM physics and analytic reconstruction/FWM_physics.tex`를
  docs consistency test가 읽어 `FileNotFoundError`.
- 저장소 지침의 `python -m pytest -q` → collection 중단. ignored AutoOD GUI reference가
  기본 requirements에 없는 optional `sympy` import. **1 collection error in 27.98 s**.

좋은 coverage: oscillator sum/7:5, pump-off skip, same-velocity correlation, crossover,
coating invariance, resolution gate, Voigt dispersion sign, split temperature, Rb-D2
finite-probe anchors, passivity, CSV absolute/relative calibration·sweep branches.

남은 trust gap:

- `test_pump_off_reproduces_autood_85rb_d1`은 같은 repository의 `ODScheme`와 비교.
  완전 독립 실험 검증 아님
  ([`test_sas.py:260-271`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/tests/test_sas.py:260>)).
- independent AutoOD source/baseline은 `.gitignore`의 `references/*` 아래. 일반 clone에서
  보존되지 않음
  ([`.gitignore:19`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/.gitignore:19>)).
- held-out measured spectrum 없음.

## 8. 기존 개선안: 계산비용·physics 보존성

| 기존 제안 | 계산비용 | physics를 거의 무비용으로 보존? | 판단 |
|---|---|---|---|
| OD/log convention 분리 | O(Nscan) scale/alias, solve 0회 | **예, exact** | P1 trust |
| cache/provenance version bump | O(1) | **예, exact** | 즉시 |
| species/line/polarization Rabi | lookup/폐형식 O(1) | **예.** pump-off exact 불변 | P1 |
| collision/pressure shift/dephasing | lookup O(1), 같은 Liouvillian 차원·solve 수 | **예.** zero-pressure/Rb85-D1 anchor 고정 가능 | P1 |
| dense scan×velocity fused arithmetic | O(Nscan Nv) 동일, 임시 배열·메모리 traffic 감소 | **예, roundoff tolerance 안.** 같은 grid/interpolation 유지 조건 | P3 성능 |
| analytic Voigt pump-off | O(Nscan Nline), 현 O(Nscan Nv Nline)보다 작음 | **수학 모델은 보존**, 현 4σ discrete measure와 bitwise 동일 아님 | 독립 gate 필요 |
| pole-residue pumped SAS | 작은 pole solve + closed/discrete average | **조건부.** 현 Γ/8 pump-population interpolation 경로를 바꿈 | reference guard·domain gate 뒤 |
| result export | O(Nscan) serialization, solve 0회 | **예** | round-trip test 남음 |
| isotope fraction | natural-Rb는 기존 두 component weight만 변경; pure→mixture는 두 isotope 비용 | defaults exact 가능 | 수요 결정 뒤 |
| fitting | 대략 `N_eval × compute`; pump-on 50회면 오늘 default 기준 약 31 s + optimizer | interactive 기본에는 아니오 | explicit action/cache/holdout |
| segmented cell | 대략 Nsegment 선형, self-consistency면 더 큼 | 1 segment exact 가능, 다중은 비쌈 | opt-in convergence |
| full VCC / mF-Zeeman / calibrated coating | velocity/state/geometry 차원 증가 | **아니오** | parked 유지 |

가장 중요한 구분: coefficient/provenance 교정은 거의 공짜. full VCC·Zeeman·segmented
depletion은 새 검증 문제. 하나의 roadmap 항목처럼 취급하면 안 됨.

## 9. 동작 바꾸지 않는 순수 코드 최적화

### 9.1 pump-off scalar average fast path — 가장 안전

현재 helper:

```python
return (prepared_population * probe_profile) @ weights
```

([`sas.py:996-998`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/schemes/sas.py:996>)).
Pump off일 때 `prepared_population == 1.0`. `probe_profile @ weights`면 큰 임시
scan×velocity 곱 배열 없음.

오늘 warm benchmark:

| 경로 | 현재 | scalar fast | speedup | 결과 |
|---|---:|---:|---:|---|
| pump off, probe 0 µW | `0.299 s` | `0.280 s` | `1.067x` | α, Re χ array exact |
| pump off, probe 100 µW | `0.426 s` | `0.403 s` | `1.057x` | α와 finite table array exact |

이 PC에서는 이전 보고서보다 이득 작음. 그래도 risk 거의 0. 신규 checklist의 3x 목표를
혼자 달성하지 못함.

### 9.2 pump-off finite-probe thermal matvec 재사용

Pump off에서 `transition_alpha`와 `thermal_alpha`가 같은 `A_t*(Lp@wt)`인데 연속 두 번
계산
([`sas.py:629-638`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/schemes/sas.py:629>)).
`thermal_alpha = transition_alpha` 사용 가능. array exact. transition당 matvec 1회 제거.

### 9.3 fused compiled velocity loop — 가장 큰 잠재 이득

현재 isotope당 `probe_base`, `deff_grid`, level별 `pop_at`, 전이별 `arg/denom/Lp/Dp`
등 1401×약 1400 배열 반복 생성. 배열 하나 약 15 MiB. 자연 Rb는 isotope 두 개를
순차 처리. 병목은 solve보다 memory allocation/traffic
([`sas.py:590-631`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/schemes/sas.py:590>),
[`checklist:2548`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/docs/checklist.json:2548>)).

권고: optional Numba fused loop. 한 `(scan, velocity)` pass에서 population interpolation,
Lorentzian/dispersion, weighted reduction 동시 수행. 기존 `deff` grid, `np.interp` 선형
규칙, ±4σ weights 그대로. physics·measure 유지. 허용 차이는 summation-order roundoff만.
Numba off NumPy fallback 유지. Pump-off·pump-on·coated·두 quadrature·finite-probe residual
전체 domain gate 필요. “3x”는 구현 전 약속 아님. benchmark/peak RSS와 함께 판정.

### 9.4 32-point Legendre nodes/weights cache

Navigate-only `cell_mm`/`line_strength` 변화마다 `leggauss(32)` 재생성
([`probe_saturation.py:195`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/probe_saturation.py:195>)).
Read-only module constant 또는 `lru_cache` 가능.

오늘 headless finite-probe observables median `20.28 → 19.39 ms`, `1.046x`. 결과 exact.
Heavy compute에는 작지만 반복 interaction에 안전.

### 9.5 낮은 우선순위·주의

- `(iso,line)` Wigner/dipole coefficient cache: exact 가능. 전체 solve 대비 매우 작음.
- scan-row block 처리: peak RSS 감소 가능. Python block overhead 때문에 속도는 반드시
  측정. 결과는 summation order 유지 시 거의 exact.
- analytic Voigt와 pole path는 **순수 behavior optimization으로 분류하지 말 것.**
  continuous Gaussian, cutoff, interpolation semantics가 바뀔 수 있음. 별도 physics gate.

## 10. 권고 순서

1. `Optical depth=-ln(T)`와 `Decadic absorbance=-log10(T)` 분리. export exact conversion
   test 추가. `cache_version` 올림.
2. `validated` 표현을 `normalization-audited/software-reference matched`로 낮춤.
   ignored reference의 source hash/baseline을 versioned fixture로 보존.
3. 기존 P1 Rabi/collision 계약 구현. lookup/같은 matrix dimension이라 성능비용 미미.
4. fused compiled velocity loop. 현 quadrature/interpolation 고정, 기존 solver guard,
   max error·RSS·numba on/off·전 domain benchmark 요구.
5. 기본 1401점 linewidth를 정량값으로 export하지 않음. fixed-window local refinement와
   convergence status 저장.
6. User Guide와 OD/SAS 예제 그림 재생성. 모든 finite-probe tier, Dispersion,
   two-temperature, paraffin, log convention, version/provenance 동기화.
7. fitting/segmented cell/VCC/full Zeeman은 opt-in 또는 parked 유지. 실제 calibration과
   held-out data 먼저.

## 11. 결론

Scheme 1은 실제 실험 물리 참고 가치 있음. 특히 **pump-off Rb weak-probe absorption**과
**qualitative/semi-quantitative SAS 구조** 강함. thermal normalization, same-velocity
pump/probe correlation, hyperfine pumping, shared absorption/dispersion profile, finite-probe
Gaussian propagation, calibrated CSV import 모두 실질적.

그러나 실험값 인용 때 반드시 함께 기록:

- `OD`가 \(-\ln T\)인지 \(-\log_{10}T\)인지,
- pump Rabi와 collision coefficient의 species/line/polarization convention,
- sub-Doppler feature가 grid-converged인지 단지 resolution gate만 통과했는지,
- finite probe가 Rb-D2 reference scope인지 generalized/sequential approximation인지.

가장 좋은 다음 작업: 새 거대 physics 추가 아님. **OD 명명·cache/provenance 수정 →
Rabi/collision P1 → 같은 measure를 유지하는 fused velocity arithmetic.** 신뢰도와 속도
동시에 개선, 기본 물리 거의 무비용 보존 가능.
