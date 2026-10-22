# GABES 일일 Scheme 물리 검토 — 2026-09-25

## 선택·결론

- Asia/Seoul 현지 날짜: 2026-09-25.
- `n = (25 mod 5) + 1 = 1`.
- 대상: **Scheme 1 — Absorption spectroscopy (OD / SAS)**. registry key `sas`.
- 검토 기준: HEAD `5b29c4f` + 현재 dirty working tree. 기존 변경 보존. production 코드 수정 없음.
- **판정:** 실제 원자물리 있음. 약한 probe OD, hyperfine 위치·적분 세기, Doppler 폭, Beer–Lambert 전파, 같은 속도군의 pump/probe 상관, hyperfine optical pumping 기반 SAS는 실험 계획에 유용. 그러나 절대 pump power 축, pumped-SAS contrast·linewidth·lock slope, buffer/coated cell, Cs, 편광·자기장 조건은 정량 실험 기준으로 부적합. **조건 탐색용 semi-quantitative reference**로 적절. 계측 표준·절대 fit model 아님.

## 현재 다섯 scheme 정의·순서

실제 등록 순서: [`gabes/schemes/__init__.py:19`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/schemes/__init__.py:19>). 내부 `ODScheme`은 등록 안 된 2-level 검증 primitive.

| n | key | 현재 정의 |
|---:|---|---|
| 1 | `sas` | OD / SAS. pump off: Doppler 흡수. pump on: Lamb dip·crossover·hyperfine pumping |
| 2 | `lambda` | Λ coherence. EIT / AT / CPT |
| 3 | `rydberg_eit` | Rydberg-EIT / microwave AT electrometry |
| 4 | `magneto` | Hanle / EIA / NMOR |
| 5 | `fwm` | seeded FWM gain 지표 / SFWM biphoton |

README 표와 일치: [`README.md:8`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/README.md:8>).

## 선행 개선안 검색

검색 범위: 12개 이전 Scheme 1 일일 보고서, `docs/checklist.json`, README/User Guide, `analysis`, 코드 `TODO/FIXME`, 로컬 issue/proposal/roadmap 파일명.

- **기존 개선안 있음.** 새 아이디어가 없는 상태 아님.
- 전용 로컬 issue-note 파일 없음. 관련 코드 TODO는 Ne broadening coefficient table·pressure shift: [`gabes/constants.py:79`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/constants.py:79>).
- P1 ready:
  - collision coefficient·pressure shift·공통 linewidth 계약: [`docs/checklist.json:2630`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/docs/checklist.json:2630>).
  - species/line/polarization별 pump Rabi·saturation provenance: [`docs/checklist.json:2670`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/docs/checklist.json:2670>).
- P3 ready: dense scan×velocity arithmetic 제거: [`docs/checklist.json:3320`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/docs/checklist.json:3320>).
- P2 ready/decision:
  - uncertainty-aware spectrum fitting: [`docs/checklist.json:3932`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/docs/checklist.json:3932>).
  - isotope fraction/원소 범위: [`docs/checklist.json:3968`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/docs/checklist.json:3968>).
  - segmented cell·pump depletion: [`docs/checklist.json:4032`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/docs/checklist.json:4032>).
- research parked:
  - full VCC: [`docs/checklist.json:3035`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/docs/checklist.json:3035>).
  - weak-probe Zeeman/polarization reference: [`docs/checklist.json:3998`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/docs/checklist.json:3998>).
  - calibrated coated-cell transport 확장: [`docs/checklist.json:2728`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/docs/checklist.json:2728>).
- 완료 기록: resolution-aware sub-Doppler gate, paraffin population-memory MVP, dispersive quadrature, density/Doppler 온도 분리: [`checklist:2691`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/docs/checklist.json:2691>), [`2712`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/docs/checklist.json:2712>), [`3876`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/docs/checklist.json:3876>), [`3904`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/docs/checklist.json:3904>).

따라서 이번 보고서에서 별도 신규 roadmap 항목 제안 안 함. 기존 P1 trust 계약 먼저.

## 구현 물리

### 유용한 실제 물리

1. **알칼리 D선 자료·hyperfine manifold.** ⁸⁵Rb/⁸⁷Rb/¹³³Cs D1/D2의 A/B, 중심 주파수, 질량, 자연선폭, abundance 사용: [`gabes/species.py:130`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/species.py:130>). 허용 F→F′, CG 기반 붕괴 분기, transit reset 구성: [`gabes/species.py:337`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/species.py:337>).
2. **약한 probe 절대 면적.** `C_F²`가 mF 합 포함. 열적 sublevel population 한 번만 계수. 적분 흡수 `C_F²/[2(2I+1)]`: [`gabes/species.py:248`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/species.py:248>). oscillator sum·manifold degeneracy ratio 테스트: [`tests/test_sas.py:91`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/tests/test_sas.py:91>).
3. **진짜 SAS 핵심.** 역행 pump `Δ+kv`, probe `Δ−kv`. 같은 velocity class에서 prepared population×probe Lorentzian 평균: [`gabes/schemes/sas.py:539`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/schemes/sas.py:539>), [`sas.py:590`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/schemes/sas.py:590>). factorized average가 아님을 테스트: [`tests/test_sas.py:249`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/tests/test_sas.py:249>). transit 감소 시 crossover 강화도 검증: [`tests/test_sas.py:339`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/tests/test_sas.py:339>).
4. **흡수·분산 공통 profile.** 같은 population difference와 Doppler quadrature로 Im/Re χ 생성: [`gabes/schemes/sas.py:539`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/schemes/sas.py:539>). 약한 probe에서 Kramers–Kronig 부호·Voigt reference 테스트 존재: [`tests/test_sas.py:594`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/tests/test_sas.py:594>).
5. **유한 probe 계산.** transition별 2-level saturation, power broadening, Gaussian radial average, `dI/dz=-α(I)I`: [`gabes/probe_saturation.py:34`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/probe_saturation.py:34>), [`probe_saturation.py:95`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/probe_saturation.py:95>), [`probe_saturation.py:166`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/probe_saturation.py:166>). Rb D2 software anchor·passivity·grid convergence 테스트 강함: [`tests/test_probe_saturation.py:212`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/tests/test_probe_saturation.py:212>), [`415`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/tests/test_probe_saturation.py:415>), [`529`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/tests/test_probe_saturation.py:529>).
6. **실험 친화 기능.** density temperature와 Doppler temperature 분리, coating population memory, CSV calibration/sweep branch 보존. 실제 cold-spot·heated body·상대/절대 detector trace 구분에 도움.

### 실험 reference 적합성

| 용도 | 판단 |
|---|---|
| isotope/hyperfine line 식별, scan 범위, 온도·길이 경향 | 좋음 |
| pump-off 약한-probe Rb OD, Doppler/Voigt scale | 유용. 단 software reference와 원자상수 closure 중심 |
| Rb D2 finite-probe 경향·Gaussian depletion | 유용한 계산 reference. 독립 실험 검증 아님 |
| pumped SAS contrast, absolute lock slope, pump power calibration | 정성/반정량 |
| high-pressure buffer cell, coated cell absolute linewidth | 부적합 |
| Cs absolute saturation/self-broadening | 부적합. Rb coefficient fallback |
| magnetic field, polarization, Faraday/dichroism | 모델 밖 |

## 오늘의 정량 감사

기본값: natural Rb D1, 40 °C, 75 mm, pump 0.5 mW, probe 0 µW, 1401 points. 3회 median. 하드웨어 의존 wall time.

| 경로 | 결과 |
|---|---:|
| 기본 SAS compute | 0.6103 s |
| SAS Tmin | 0.1180495075 |
| SAS peak αL | 2.136651188 |
| SAS peak −log₁₀T | 0.927935821 |
| 표시 sub-Doppler FWHM | 18.87 MHz |
| FWHM당 sample | 2.3 → `resolution-limited` |
| pump-off OD compute | 0.2971 s |
| pump-off Tmin / peak αL / peak −log₁₀T | 0.105095655 / 2.252884348 / 0.978415241 |

해상도 gate는 안전하게 lock slope hero를 숨김: [`tests/test_sas.py:406`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/tests/test_sas.py:406>). 고정 window 4×/8× 수렴 테스트는 별도 존재: [`tests/test_sas.py:442`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/tests/test_sas.py:442>). 기본 1401점 숫자 자체가 수렴값은 아님.

### pump Rabi 규약 감사

현재 pump 변환은 모든 종/선에 `I_sat=4.484 mW/cm²` 사용: [`gabes/species.py:272`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/species.py:272>). 반면 finite-probe helper의 같은 closed-two-level 식 `πhcΓ/(3λ³)`은 line-specific: [`gabes/probe_saturation.py:46`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/probe_saturation.py:46>).

| line | helper I_sat [mW/cm²] | 현재 pump Ω / helper식 Ω |
|---|---:|---:|
| Rb D1 | 1.49584 | 1.731 |
| Rb D2 | 1.66934 | 1.639 |
| Cs D1 | 0.83270 | 2.321 |
| Cs D2 | 1.10238 | 2.017 |

이 표가 helper식의 편광·open-transition 적용을 자동 정당화하지는 않음. 오히려 **현재 pump와 probe convention이 닫히지 않았음**을 수치화. 기존 P1 `species-line-rabi-and-saturation-provenance` 정확히 해당. correct 값은 polarization·mF·transition 정의와 독립 reference를 함께 고정한 뒤 결정.

### collision/linewidth 계약

- probe profile: `gamma_eff = gamma_natural + self + buffer`: [`gabes/schemes/sas.py:400`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/schemes/sas.py:400>).
- pump state Liouvillian: `build_manifold`의 자연 붕괴·transit 사용. `gamma_eff` optical dephasing 미주입: [`sas.py:393`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/schemes/sas.py:393>), [`sas.py:568`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/schemes/sas.py:568>).
- Cs도 `species.self_broadened_gamma`의 단일 Rb β 사용: [`gabes/species.py:210`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/species.py:210>).

따라서 buffer/self broadening이 probe width만 넓히고 pump-prepared population의 optical coherence width에는 같은 방식으로 안 들어감. 압력↑ 조건의 SAS hole/crossover 정량성 낮음.

## 기존 제안: 계산비용·물리 보존

| 기존 제안 | runtime 비용 | 물리 보존 판단 |
|---|---|---|
| species/line Rabi·saturation provenance | O(1) 계수 조회. matrix 크기·solve 수 동일 | **무시 가능한 overhead.** pump-on 물리는 의도적으로 교정. pump-off OD exact 유지 가능 |
| collision table·pressure shift·공통 linewidth API | O(1) lookup·detuning shift. dephasing을 기존 Liouvillian에 넣어 solve 차원 동일 | **대체로 negligible overhead.** 값은 바뀌지만 더 일관된 실제 물리. zero-pressure anchor 보존 가능 |
| OD/absorbance 명칭·software/experiment provenance·cache metadata | solve 0 | 물리 exact. trust 개선 |
| dense velocity fused arithmetic | O(Nscan Nv) 동일. 임시배열·memory traffic 감소 | 동일 1 m/s, 4σ weights·`np.interp` 순서를 허용오차 내 보존하면 pure optimization |
| isotope fraction | isotope 2개 합산 구조 유지. single run은 거의 동일 | 값만 선택한 조성대로 바뀜. fit까지 하면 N_eval배 |
| experimental fit | 적어도 N_eval×compute; identifiability·uncertainty 필요 | negligible 아님. 독립 데이터 필요 |
| segmented pump depletion | 대략 Nz×atomic response 또는 surrogate/caching 필요 | negligible 아님. uniform-slab physics 변경 |
| full VCC / mF-Zeeman / coating transport | velocity/state/geometry 차원 증가 | negligible 아님. opt-in research reference가 맞음 |

우선순위: 기존 두 P1은 실행비용이 거의 없고 실험 신뢰 이득 큼. dense arithmetic은 그 다음. fitting/segmentation/full VCC는 데이터·수렴 budget 확보 뒤.

## 동작 불변 speed 개선

오늘 기본 SAS cProfile 0.630 s:

- 두 isotope `_component_response`: cumulative 0.614 s.
- `np.interp`: 0.164 s.
- `_velocity_correlated_average`: 0.094 s.
- pump population solve: 0.071 s.

근거 병목: [`gabes/schemes/sas.py:539`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/schemes/sas.py:539>).

1. **pump-off scalar fast path.** helper는 `prepared_population=1.0`에도 `(1.0 * profile) @ weights`, 즉 전체 복사성 곱 수행: [`sas.py:996`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/schemes/sas.py:996>). 1401×1403 microbench median 5.047 ms → 직접 `profile @ weights` 0.564 ms, **8.94×**, SHA-256 동일. pump-off transition마다 적용 가능.
2. **pump-off finite-probe 중복 matvec 제거.** `transition_alpha` 계산 직후 pump-off이면 `thermal_alpha`가 같은 `A_t*(Lp@wt)`인데 다시 수행: [`sas.py:629`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/schemes/sas.py:629>). `thermal_alpha=transition_alpha` 재사용. array exact 가능.
3. **기존 P3 fused loop.** `probe_base`, `deff_grid`, level별 `pop_at`, `Lp`, `Dp`의 1401×~1400 임시배열 반복. 동일 interpolation·합산 순서 기준 compiled fused kernel 검증. 현재 가장 큰 전체 compute 후보.
4. **Legendre node cache.** 매 finite-probe readout마다 `leggauss(32)` 재생성: [`gabes/probe_saturation.py:195`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/probe_saturation.py:195>). `(order→nodes,weights)` immutable cache. 작은 이득.
5. `line_integrated_alpha`/Wigner cache는 profile에서 합계 약 0.002 s. 우선순위 낮음.

## 문서·예제·테스트 품질

### 문서·예제

- README는 finite-probe 제한, thermal normalization, coating, 두 온도를 비교적 정확히 기술: [`README.md:376`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/README.md:376>).
- 그러나 in-app About은 ⁸⁵Rb D1이 AutoOD에 `<0.1%`라고 계속 표현: [`gabes/schemes/sas.py:296`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/schemes/sas.py:296>). 현재 D1 테스트는 같은 저장소 `ODScheme`와 SAS 경로 비교: [`tests/test_sas.py:260`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/tests/test_sas.py:260>). independent experimental validation처럼 읽히면 안 됨.
- User Guide도 0.1% 문구 반복: [`docs/Userguide/GABES_User_Guide_v2.html:409`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/docs/Userguide/GABES_User_Guide_v2.html:409>), [`563`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/docs/Userguide/GABES_User_Guide_v2.html:563>), [`922`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/docs/Userguide/GABES_User_Guide_v2.html:922>).
- Guide의 finite-probe 설명은 natural Rb D2만 언급: [`User Guide:553`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/docs/Userguide/GABES_User_Guide_v2.html:553>). 실제 UI는 모든 built-in species/line에 노출하되 Rb D2 reference-scope, D1/Cs generalized.
- 독립 실행 가능한 Scheme 1 example script/config 없음. 예제는 Guide의 정적 [`od.png`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/docs/Userguide/userguide_assets/od.png>)·[`sas.png`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/docs/Userguide/userguide_assets/sas.png>) 중심. 2026-06-08 생성. 그림 OD Tmin≈0.27, SAS Tmin≈0.28. 오늘 current 값 0.105/0.118과 불일치. normalization 변경 뒤 재생성 필요.
- `references/*`는 gitignore 대상: [`.gitignore:19`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/.gitignore:19>). 일반 clone에서 AutoOD reference 재현성 약함. pinned test number는 남지만 source provenance bundle 아님.

### 테스트

- Scheme 1 집중: `test_sas.py`, `test_absorption.py`, `test_probe_saturation.py`, `test_experimental_csv.py`, `test_lineshape.py` → **133 passed / 62.94 s**.
- 강점: oscillator sum, hyperfine split, same-velocity correlation, pump-off limit, crossover pumping, coating opt-in/inert limit, dispersion sign, two-temperature contract, finite-probe anchors/passivity/convergence, CSV calibration·sweep direction.
- 빈칸: independent measured OD/SAS spectrum holdout, absolute pump Rabi/polarization, buffer-pressure shift/dephasing, Cs collision data, full experimental fit residual.
- repository 요구 `python -m pytest -q`: **collection error / 24.08 s**. `references/AutoOD/test_autood_gui_adjustments.py` import 중 optional `sympy` 없음. Scheme 1 집중 테스트 실패 아님.

## 연구자용 최종 권고

1. 현재 출력은 line assignment, scan span, pump-off OD scale, 상대 knob trend에 사용.
2. 기본 SAS FWHM·lock slope는 `resolution-limited`; 논문/실험 기준값으로 인용 금지. 고정 feature window grid 수렴 필요.
3. 기존 P1 Rabi와 collision 계약 먼저. 둘 다 runtime 거의 그대로. 특히 pump/probe가 같은 linewidth·frequency convention을 쓰게 할 것.
4. 문서의 “validated / 0.1%”를 “normalization-corrected software reference agreement”로 통일. 현재 그림 재생성. 재현 가능한 Scheme 1 example config 추가.
5. 그 뒤 P3 fused velocity arithmetic. 동일 quadrature·interpolation·array tolerance gate 필수.
6. segmentation, VCC, Zeeman/polarization은 무비용 개선 아님. 측정 데이터와 opt-in reference budget 준비 전 기본 경로에 넣지 말 것.

