# GABES 일일 Scheme 물리 검토 — 2026-10-05

## 선택·요약

- 현지 날짜: 2026-10-05, Asia/Seoul.
- `n = (5 mod 5) + 1 = 1`.
- 대상: **Scheme 1 — Absorption spectroscopy (OD / SAS)**. registry key `sas`.
- 기준: HEAD `fcb312d` + 현재 dirty working tree. 기존 변경 보존. production 코드 수정 없음.
- 2026-09-25 검토 뒤 Scheme 1 핵심 코드·테스트 커밋 변화 없음. 이후 두 커밋은 FWM/Grand Challenge 중심. README·checklist만 관련 파일 목록에 포함.
- **판정:** 실제 원자물리 있음. 약한-probe OD, hyperfine 위치·세기, Doppler/Voigt 폭, Beer–Lambert 전파, 같은 속도군 pump–probe 상관, hyperfine optical pumping SAS 구현. line assignment, scan 범위, pump-off OD, 상대 knob trend에 유용. 절대 pump Rabi, pumped-SAS contrast·linewidth·lock slope, buffer/coated cell, Cs, 편광·자기장 정량 기준으로는 부족. **조건 탐색용 semi-quantitative experimental reference**. 계측 표준·절대 fit model 아님.

## 현재 다섯 scheme 정의·순서

실제 dropdown 순서: [`gabes/schemes/__init__.py:19`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/schemes/__init__.py:19>). 내부 `ODScheme`은 등록 안 된 2-level 검증 primitive: [`absorption.py:211`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/schemes/absorption.py:211>).

| n | key | 정의 | 코드 |
|---:|---|---|---|
| 1 | `sas` | Absorption spectroscopy. pump off OD, pump on SAS | [`sas.py:175`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/schemes/sas.py:175>) |
| 2 | `lambda` | Λ coherence. EIT / AT / CPT | [`absorption.py:466`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/schemes/absorption.py:466>) |
| 3 | `rydberg_eit` | Rydberg-EIT electrometry / microwave AT | [`rydberg.py:94`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/schemes/rydberg.py:94>) |
| 4 | `magneto` | Magneto-optics. Hanle / EIA / MOR/NMOR | [`magneto.py:121`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/schemes/magneto.py:121>) |
| 5 | `fwm` | Four-wave mixing. Squeezing / Biphoton | [`fwm.py:3626`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/schemes/fwm.py:3626>) |

README 표와 일치: [`README.md:8`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/README.md:8>).

## 선행 개선안·TODO·issue 검색

검색 범위:

- 이전 Scheme 1 일일 보고서 13개. 최신: [`2026-09-25_scheme-1_od-sas.md`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/docs/daily_report/2026-09-25_scheme-1_od-sas.md>).
- 현재 `docs/checklist.json`, README, User Guide, 코드·테스트 `TODO/FIXME`, local issue/proposal/roadmap 파일명.
- 공개 GitHub Issues API. 2026-10-05 open/closed 합계 **0건**: [repository issues](https://github.com/Shake2313/fwm-squeezing-app/issues).

결과:

- **기존 개선안 있음.** 새 roadmap 부재 상태 아님.
- 별도 local issue-note 없음.
- 관련 코드 TODO 한 건: Ne coefficient table·pressure shift: [`constants.py:79`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/constants.py:79>).
- 2026-09-25 뒤 관련 checklist 상태 변화 없음.

| 기존 항목 | 상태 | 계산비용 | 물리 보존·가치 |
|---|---|---|---|
| collision coefficient provenance + pressure shift | P1 `ready` ([`checklist:2839`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/docs/checklist.json:2839>)) | O(1) table lookup·detuning shift. 기존 Liouvillian 차원·solve 수 유지 | **negligible overhead.** zero-pressure anchor exact 보존 가능. unsupported fallback 제거, pump/probe linewidth 계약 통일 |
| species/line Rabi + saturation provenance | P1 `ready` ([`checklist:2879`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/docs/checklist.json:2879>)) | O(1) coefficient 선택 | **negligible overhead.** pump-off OD 불변 가능. pump-on 축 의도적 교정 |
| dense scan×velocity arithmetic 제거 | P3 `ready` ([`checklist:3529`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/docs/checklist.json:3529>)) | 같은 O(`Nscan Nv`) 또는 pole reduction | fused loop면 물리 exact, 합산 허용오차만 검증. pole형은 기존 discrete measure 차이 별도 검증 필요 |
| experimental spectrum fitting | P2 `ready` ([`checklist:4141`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/docs/checklist.json:4141>)) | 최소 `N_eval × compute`; SAS는 매 평가 pump table solve | negligible 아님. uncertainty·identifiability·holdout 필수 |
| isotope fraction / element coverage | P2 `needs_decision` ([`checklist:4177`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/docs/checklist.json:4177>)) | forward weight만이면 작음; fit이면 `N_eval` 증가 | 측정 cell 필요성 결정 전 보류 적절 |
| weak-probe Zeeman/polarization reference | P2 research `parked` ([`checklist:4207`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/docs/checklist.json:4207>)) | Hamiltonian diagonalization + polarization propagation | 값 큰 개선. 무비용 아님. 측정 geometry·holdout 필요 |
| segmented cell / pump depletion | P2 `ready` ([`checklist:4241`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/docs/checklist.json:4241>)) | 대략 `Nz × atomic response`, cache/surrogate 없으면 큼 | uniform slab 물리 변경. 수렴 parameter 필요 |
| full VCC | P2 research `parked` ([`checklist:3244`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/docs/checklist.json:3244>)) | velocity classes 결합. 현재 separability 소실 | negligible 아님. opt-in reference가 맞음 |
| coated-cell transport 확장 | P2 research `parked` ([`checklist:2937`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/docs/checklist.json:2937>)) | dark-time/geometry/Zeeman-coherence 상태 추가 | cell별 데이터 없으면 식별 불가 |

신규 physics 항목 추가 안 함. 기존 P1 두 건 먼저. 계산비 거의 0, 신뢰 경계 개선 큼.

## 구현 물리 감사

### 실제로 유용한 물리

1. **알칼리 D선·hyperfine manifold.** ⁸⁵Rb/⁸⁷Rb/¹³³Cs D1/D2 상수, 허용 F→F′, CG 세기·붕괴 분기, thermal ground population: [`species.py:127`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/species.py:127>), [`species.py:336`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/species.py:336>).
2. **약한-probe 절대 면적.** `C_F²`에 magnetic sublevel 합 포함. 적분 흡수 `A_t ∝ C_F² |d|² N/[2(2I+1)]`; thermal population 이중계수 방지: [`species.py:248`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/species.py:248>). oscillator sum·manifold ratio 테스트: [`test_sas.py:91`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/tests/test_sas.py:91>).
3. **SAS 핵심 kinematics.** pump `Δ+kv`, probe `Δ−kv`; 같은 velocity class의 prepared population과 probe Lorentzian 곱 평균: [`sas.py:539`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/schemes/sas.py:539>), [`sas.py:590`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/schemes/sas.py:590>). factorized average 아님: [`test_sas.py:249`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/tests/test_sas.py:249>).
4. **Hyperfine optical pumping.** pump steady state, CG-branched decay, transit reload. crossover 증강·반전 가능: [`species.py:373`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/species.py:373>), [`sas.py:563`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/schemes/sas.py:563>). transit 감소 시 crossover 강화 테스트: [`test_sas.py:339`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/tests/test_sas.py:339>).
5. **흡수·분산 한 complex profile.** Im/Re χ 같은 population difference·Doppler quadrature 공유: [`sas.py:539`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/schemes/sas.py:539>). 독립 Voigt·부호 fixture: [`test_sas.py:594`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/tests/test_sas.py:594>), [`test_sas.py:607`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/tests/test_sas.py:607>).
6. **유한 probe.** line별 closed-two-level saturation, power broadening, longitudinal depletion, Gaussian radial average: [`probe_saturation.py:46`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/probe_saturation.py:46>), [`probe_saturation.py:95`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/probe_saturation.py:95>), [`probe_saturation.py:166`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/probe_saturation.py:166>). Rb D2 software anchor·passivity·grid convergence 테스트: [`test_probe_saturation.py:212`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/tests/test_probe_saturation.py:212>), [`test_probe_saturation.py:415`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/tests/test_probe_saturation.py:415>), [`test_probe_saturation.py:529`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/tests/test_probe_saturation.py:529>).
7. **실험 친화 기능.** density/Doppler temperature 분리, paraffin population memory, CSV absolute/relative calibration·sweep branch 보존. cold spot, heated body, detector trace 의미 구분에 도움.

### 정량 reference 한계

| 용도 | 판단 |
|---|---|
| isotope/hyperfine line 식별, scan span, 온도·길이 경향 | 좋음 |
| pump-off weak-probe Rb OD, Doppler/Voigt scale | 유용. software reference·원자상수 closure 중심 |
| Rb D2 finite-probe 경향·Gaussian depletion | 계산 reference로 유용. 독립 실험 검증 아님 |
| pumped SAS contrast, absolute lock slope, pump-power calibration | 정성/반정량 |
| high-pressure buffer cell, coated-cell absolute linewidth | 부적합 |
| Cs absolute saturation/self-broadening | 부적합. Rb coefficient fallback |
| magnetic field, polarization, Faraday/dichroism | 모델 밖 |

## 오늘 수치 재현

기본값: natural Rb D1, 40 °C, 75 mm, pump 0.5 mW, probe 0 µW, waist 1 mm, transit 100 kHz, 1401 points. warm compute 3회 median. wall time은 장비·부하 의존.

| 경로 | 결과 |
|---|---:|
| 기본 SAS compute | 0.8794 s median |
| SAS `Tmin` | 0.1180495075 |
| SAS peak `αL` | 2.136651188 |
| 표시 sub-Doppler FWHM | 18.8684 MHz |
| samples/FWHM | 2.3418 |
| resolution status | `resolution-limited` |
| pump-off OD compute | 0.4166 s median |
| OD `Tmin` | 0.1050956546 |
| OD peak `αL` | 2.252884348 |

수치 자체는 2026-09-25 결과와 동일. runtime만 변동. 현재 기본 FWHM은 실험 hero 값으로 인용 금지. 코드도 6 samples/FWHM 전에는 lock slope를 숨김: [`sas.py:1254`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/schemes/sas.py:1254>), [`test_sas.py:406`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/tests/test_sas.py:406>). 4×/8× 고정-window 수렴 fixture 별도 존재: [`test_sas.py:442`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/tests/test_sas.py:442>).

### pump Rabi 규약: 이전 보고서 비율 표 수정

현재 pump 변환은 모든 종·선에 `I_sat = 4.484 mW/cm²`: [`constants.py:30`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/constants.py:30>), [`species.py:272`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/species.py:272>). finite-probe helper는 line별 `I_sat = πhcΓ/(3λ³)`: [`probe_saturation.py:46`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/probe_saturation.py:46>).

| line | helper `I_sat` [mW/cm²] | 현재 pump Ω / helper식 Ω | helper식 Ω / 현재 pump Ω |
|---|---:|---:|---:|
| Rb D1 | 1.49584 | 0.5776 | 1.731 |
| Rb D2 | 1.66934 | 0.6102 | 1.639 |
| Cs D1 | 0.83270 | 0.4309 | 2.321 |
| Cs D2 | 1.10238 | 0.4958 | 2.017 |

2026-09-25 보고서는 마지막 열 숫자를 “현재/helper”로 잘못 표기. 실제로 현재 pump Ω가 helper식보다 **작음**. 여기서 helper식이 open hyperfine·mF·편광 pump의 정답이라는 뜻 아님. 방향 오류만 수정. 올바른 교정은 기존 P1처럼 transition·polarization 정의와 독립 reference를 함께 고정해야 함.

### collision/linewidth 계약

- probe `gamma_eff = natural + self + buffer`: [`sas.py:400`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/schemes/sas.py:400>).
- pump Liouvillian은 자연 붕괴·transit 사용. 같은 `gamma_eff` elastic optical dephasing 없음: [`sas.py:393`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/schemes/sas.py:393>), [`sas.py:568`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/schemes/sas.py:568>).
- Cs도 단일 Rb self-broadening β 사용: [`species.py:210`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/species.py:210>).

결과: pressure↑ 조건에서 probe 폭과 pump-prepared population의 optical width 계약 불일치. 현재 SAS hole/crossover 정량성 낮음. P1 collision item이 정확한 해결 위치.

## 동작 불변 code optimization

오늘 cProfile, 기본 SAS 한 번 0.979 s:

- `_component_response`: 0.932 s cumulative.
- `np.interp`: 0.235 s.
- `_velocity_correlated_average`: 0.133 s.
- `_pump_pops`: 0.176 s.
- dense 한 배열: ⁸⁵Rb `1401×1403 = 14.996 MiB`, ⁸⁷Rb `1401×1387 = 14.825 MiB`. `probe_base`, `deff_grid`, `pop_at`, `arg`, `denom`, `Lp`, `Dp` 반복.

근거 hot path: [`sas.py:539`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/schemes/sas.py:539>), [`sas.py:590`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/schemes/sas.py:590>).

1. **pump-off scalar fast path.** `_velocity_correlated_average(1.0, profile, weights)`가 `1.0 * profile` 전체 복사 뒤 matvec: [`sas.py:996`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/schemes/sas.py:996>). synthetic `1401×1403`, 9회 median: **7.6795 → 0.6484 ms, 11.84×**. `np.array_equal=True`, max diff 0. pump-off profile cProfile에서 이 helper 16회, 0.115 s cumulative. `np.isscalar(prepared_population)` 분기 뒤 `probe_profile @ weights` 사용. 물리·배열 exact.
2. **pump-off finite-probe 중복 matvec 제거.** `transition_alpha` 계산 직후 pump-off `thermal_alpha = A_t*(Lp@wt)` 다시 계산: [`sas.py:629`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/schemes/sas.py:629>), [`sas.py:633`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/schemes/sas.py:633>). scalar fast path와 함께 `thermal_alpha = transition_alpha` 재사용 가능. 동작 exact.
3. **기존 P3 fused compiled loop.** pump-on은 `np.interp`와 dense temporaries가 주 병목. velocity index 순서, 1 m/s·4σ discrete weights, interpolation rule, transition 합산 순서 고정한 Numba loop. 물리 exact, floating 허용오차만 gate. 3× 목표는 checklist에 이미 있음.
4. **32점 Legendre node cache.** finite-probe observable 호출마다 `leggauss(32)` 재생성: [`probe_saturation.py:195`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/probe_saturation.py:195>). order-key immutable cache 가능. 작은 이득. 반환 배열 mutation 방지 필요.
5. **finite-probe per-frequency inversion compiled loop.** `propagate_gaussian_spectrum`은 frequency마다 `np.interp` 두 번: [`probe_saturation.py:236`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/probe_saturation.py:236>). monotone column search를 compiled loop로 옮기면 수식·quadrature 유지. finite-probe batch에서만 가치. 먼저 profile 필요.

Wigner/manifold는 이미 `build_manifold` LRU cache: [`species.py:336`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/species.py:336>). 추가 cache 우선순위 낮음.

## 문서·예제·테스트

### 문서·예제

- README는 finite-probe 범위, thermal normalization, coating, 두 온도, dispersion 한계를 비교적 정확히 기록: [`README.md:386`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/README.md:386>).
- 그러나 scheme About은 ⁸⁵Rb D1이 AutoOD와 `<0.1%`라고 표현: [`sas.py:301`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/schemes/sas.py:301>). 현재 D1 핵심 test는 같은 저장소 내부 `ODScheme`과 SAS 경로 closure: [`test_sas.py:260`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/tests/test_sas.py:260>). 독립 실험 validation처럼 읽히면 안 됨.
- User Guide도 0.1% 문구 반복: [`GABES_User_Guide_v2.html:409`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/docs/Userguide/GABES_User_Guide_v2.html:409>), [`:563`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/docs/Userguide/GABES_User_Guide_v2.html:563>), [`:922`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/docs/Userguide/GABES_User_Guide_v2.html:922>).
- Guide는 probe knob를 natural Rb D2에서만 보인다고 쓰지만, 현재 UI는 모든 built-in species/line에 노출하고 D1/Cs를 generalized로 label: [`Guide:556`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/docs/Userguide/GABES_User_Guide_v2.html:556>), [`sas.py:195`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/schemes/sas.py:195>). Guide 기본 단위도 `0 mW`; 현재 knob는 `0 µW`: [`Guide:571`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/docs/Userguide/GABES_User_Guide_v2.html:571>).
- standalone `examples/` 디렉터리 없음. Scheme 1 실행 예제 없음. 현재 예제 역할은 tests와 User Guide 정적 그림.
- Guide `od.png`, `sas.png` 생성 시각 2026-06-08. 현재 normalization 뒤 수치와 불일치 가능. 2026-09-25에서 확인한 그림 `Tmin≈0.27/0.28` 대 현재 `0.105/0.118`. 재생성 필요.

### 테스트

집중 명령:

```text
python -m pytest -q tests/test_sas.py tests/test_absorption.py tests/test_probe_saturation.py tests/test_experimental_csv.py tests/test_lineshape.py
133 passed in 92.89s
```

강점: oscillator sum, hyperfine split, same-velocity correlation, pump-off limit, crossover pumping, coating opt-in/inert limit, dispersion sign, two-temperature contract, finite-probe anchors/passivity/convergence, CSV calibration·sweep direction.

빈칸: independent measured OD/SAS holdout, absolute pump Rabi/polarization, pressure shift·pump/probe common dephasing, Cs collision data, fitted residual·uncertainty.

저장소 필수 명령:

```text
python -m pytest -q
ERROR collecting references/AutoOD/test_autood_gui_adjustments.py
ModuleNotFoundError: No module named 'sympy'
1 error during collection in 37.30s
```

Scheme 1 집중 실패 아님. optional reference GUI가 test collection 때 `sympy`를 hard-import하는 환경 문제: [`rb_od_calculator_optimized_final.py:5`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/archive/AutoOD/rb_od_calculator_optimized_final.py:5>).

## 최종 권고

1. 현재 출력: line assignment, scan span, pump-off OD scale, 상대 knob trend에 사용.
2. 기본 SAS FWHM·lock slope: `resolution-limited`. 정량 인용 금지. 고정 feature window 수렴값 사용.
3. 기존 P1 Rabi·collision 계약 먼저. runtime 거의 그대로. 오늘 수정된 비율 방향상 현재 pump Ω는 closed-two-level helper보다 0.431–0.610배.
4. 문서의 “validated / 0.1%”를 “normalization-corrected software-reference agreement”로 통일. Guide knob 범위·단위 수정. 그림 재생성. 재현 가능한 Scheme 1 example config 추가.
5. pure speed: pump-off scalar 분기 즉시 적용 후보. 다음 fused scan×velocity loop. exact-discrete-measure regression 필수.
6. fitting, segmentation, VCC, Zeeman/polarization: 무비용 개선 아님. 측정 데이터·holdout·수렴 budget 뒤 opt-in reference로 진행.
