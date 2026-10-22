# GABES 일일 Scheme 물리 검토 — 2026-10-15

## 선택·현재 순서

- 현지 시각: 2026-10-15, Asia/Seoul.
- `n = (15 mod 5) + 1 = 1`.
- 대상: **Scheme 1 — Absorption spectroscopy (OD / SAS)**, registry key `sas`.
- 기준: `main` HEAD `39be4ae`와 검토 시점 dirty working tree. 기존 변경 보존. production 코드 수정 안 함.

실제 UI 순서는 [`_SCHEMES`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/schemes/__init__.py:19>) 배열 순서다. `ODScheme`은 dropdown 미등록 2-level 검증 primitive다.

| n | key | 현재 정의 | 구현 |
|---:|---|---|---|
| 1 | `sas` | OD / SAS: pump off 선형 흡수, pump on 포화흡수 | [`SASScheme`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/schemes/sas.py:216>) |
| 2 | `lambda` | Λ coherence: EIT / AT / CPT | [`LambdaScheme`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/schemes/absorption.py:466>) |
| 3 | `rydberg_eit` | Rydberg-EIT / microwave AT electrometry | [`RydbergEITScheme`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/schemes/rydberg.py:94>) |
| 4 | `magneto` | Magneto-optics: Hanle / EIA / MOR·NMOR | [`MagnetoScheme`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/schemes/magneto.py:121>) |
| 5 | `fwm` | FWM: mean-field Squeezing indicator / Biphoton | [`FWMScheme`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/schemes/fwm.py:3626>) |

## 현재 코드 정체성과 재현성

- Scheme 1 주 파일 SHA-256: `sas.py = 48FAAEE8…`, `species.py = 932A47C8…`, `probe_saturation.py = B8495B3E…`.
- 세 파일 수정 시각과 blob은 2026-10-10 검토 뒤 불변. 이전 기본 수치·profile을 다시 계산하지 않고 해당 봉인 결과를 재사용한다.
- 그러나 현재 working-tree blob(`sas.py` `4968d9b…`, `species.py` `dea1a94…`)은 clean HEAD blob(`ee4376f…`, `5b9a7a8…`)과 다르다.
- 차이는 단순 문서가 아니다. 현재 작업 트리는 `C_F²`가 자기부준위를 이미 합했다는 점을 반영해 ground-manifold population 중복 계수를 제거한다: [`line_integrated_alpha`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/species.py:248>), [`ODScheme._hyperfine_alpha`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/schemes/absorption.py:147>).
- 현재 ⁸⁵Rb D1 ground-group 적분비는 7/5. clean HEAD는 옛 49/25 경로. 독립 lifetime sum rule과 모든 Rb/Cs D1/D2 degeneracy ratio 검사가 현재 경로를 고정한다: [`test_thermal_absorption_oscillator_sum_and_manifold_ratio`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/tests/test_sas.py:91>).

## 선행 개선안·TODO·issue 검색

Scheme 1 일일 보고서 15개(2026-06-25~2026-10-10), 현행 [`docs/checklist.json`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/docs/checklist.json>), README, User Guide, 코드 `TODO/FIXME`, issue/proposal/DEVLOG 파일명을 먼저 검색했다.

- **기존 개선안 충분.** 새 physics roadmap 부재 상태 아님.
- Scheme 1 전용 local issue/DEVLOG 없음. `.github` 디렉터리도 없음. 현행 issue ledger는 checklist와 이전 일일 보고서.
- 직접 코드 TODO는 Ne broadening scalar를 gas/species/line coefficient table과 pressure shift로 승격하라는 항목: [`constants.py:79`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/constants.py:79>).
- 2026-10-10 뒤 Scheme 1 코드·전용 보고서 변경 없음. checklist 전역 파일은 바뀌었으나 아래 Scheme 1 항목의 상태·내용은 그대로.
- ledger 불일치 지속. `od-sas-isotope-fraction-and-element-coverage`는 아직 fixed menu / `needs_decision`라 적지만 custom Rb fraction은 구현·테스트 완료: [`checklist:4363`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/docs/checklist.json:4363>), [`_isotope_composition`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/schemes/sas.py:94>).
- 2026-10-10 공개 GitHub issue 0건 기록은 이전 보고서에 있음. 이번 환경은 `gh` 미설치, GitHub issue endpoint 접근 실패. 현 원격 issue 수를 새로 단정하지 않음.

## ⟦AI⟧ 기존 제안 비용·물리 보존 판정

| 기존 항목 | 상태·근거 | 계산비 | 판정 |
|---|---|---:|---|
| collision coefficient provenance, pressure shift, 공통 optical dephasing | P1 `ready`; [`checklist:3025`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/docs/checklist.json:3025>) | coefficient lookup·폭/중심 갱신 O(1), 같은 matrix 차원·solve 수 | **runtime negligible.** zero-pressure anchor exact 보존 가능. Cs가 Rb self-broadening을 상속하는 오류와 pump/probe linewidth 불일치 교정. 문헌·단위·부호 검증 effort 큼. |
| species/line/polarization별 pump Rabi·saturation | P1 `ready`; [`checklist:3065`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/docs/checklist.json:3065>) | O(1) lookup/폐형식 | **runtime negligible.** pump-off OD exact 불변. pump-on power 축 의도적 교정. open transition·편광 convention 필요. |
| custom Rb isotope fraction | 코드는 완료, ledger stale | natural Rb와 같은 최대 두 isotope solve. endpoint는 한 component 생략 | 총 elemental density·고정 ⁸⁵Rb centroid·legacy natural 배열 보존. 새 원소와 ledger 상태는 분리 필요. |
| affine CSV alignment | 구현 완료 | 고정 이론 뒤 bounded 후처리, OBE solve 0 | 유용한 display match. 물리 파라미터 fit·주파수 calibration·검증 아님: [`fit_trace_alignment`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/experimental_alignment.py:199>). |
| physical spectrum fitting | P2 `ready`; [`checklist:4327`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/docs/checklist.json:4327>) | 약 `N_eval × compute`, pumped SAS면 매 평가 pump table | negligible 아님. uncertainty·identifiability·fit/holdout 분리 필요. |
| dense scan×velocity 제거 | P3 `ready`; [`checklist:3715`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/docs/checklist.json:3715>) | fused loop는 같은 O(`Nscan Nv`), 메모리 traffic 감소 | 같은 discrete measure·보간·합산 규칙이면 behavior 보존 가능. Pole/연속 Maxwell 변경은 별도 수치모델. |
| weak-probe Zeeman/polarization reference | P2 research `parked`; [`checklist:4393`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/docs/checklist.json:4393>) | hyperfine-Zeeman diagonalization + mode/Jones propagation | 실험 가치 큼. 무부하 아님. 측정 geometry·polarization·holdout 선행. |
| segmented cell / pump depletion | P2 `ready`; [`checklist:4427`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/docs/checklist.json:4427>) | 대략 `Nz × atomic response` | high-OD pumped SAS에 필요. 수렴 parameter를 가진 opt-in reference가 맞음. |
| full VCC / coated transport | research parked | velocity·history·geometry 상태 결합 | negligible 아님. cell별 자료 없으면 식별 불가. |
| export / provenance | P2 `in_progress`; [`checklist:4454`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/docs/checklist.json:4454>) | O(`Nscan`) 직렬화, solve 0 | physics exact. named arrays·round trip·lazy view 범위 남음. |

가성비 순서: **collision/Rabi provenance → 문서·ledger 동기화 → exact arithmetic 최적화**. Full VCC·Zeeman·segmentation은 같은 범주 아님.

## 구현된 실제 물리

1. **원자·hyperfine 구조.** ⁸⁵Rb/⁸⁷Rb/¹³³Cs D1/D2의 A/B, centroid, 질량, 자연선폭: [`Isotope`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/species.py:130>). Wigner-6j/3j line strength와 CG-branched decay: [`cf2`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/species.py:225>), [`build_manifold`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/species.py:336>).
2. **약한-probe 절대 면적.** 현재 식은
   `A_t = π k C_F² |d|² N / [ε₀ ħ 2(2I+1)]`.
   `C_F²`의 자기부준위 합과 thermal population을 한 번만 센다: [`line_integrated_alpha`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/species.py:248>). AutoOD-Rb도 2026-09-06 같은 오류를 제거하고 ⁸⁵Rb 7/5, ⁸⁷Rb 5/3을 기록: [`AutoOD-Rb README:114`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/references/AutoOD-Rb/README.txt:114>).
3. **실제 SAS kinematics.** 같은 velocity class에 pump `Δ+kv`, counter-propagating probe `Δ−kv` 적용. factorized 평균 아님: [`_component_response`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/schemes/sas.py:600>), [`same-class test`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/tests/test_sas.py:249>).
4. **Hyperfine optical pumping.** pump steady state, CG-branched decay, thermal transit reload가 population transfer와 crossover 증강·반전을 만듦: [`pump_hamiltonian`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/species.py:419>), [`pump table`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/schemes/sas.py:624>), [`crossover test`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/tests/test_sas.py:339>).
5. **흡수·분산 한 complex profile.** 같은 population difference와 Maxwell weight로 Im/Re χ 계산. `α=k Imχ`, Kramers–Kronig 부호와 zero crossing 검사 있음: [`sas.py:603`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/schemes/sas.py:603>), [`test_sas.py:594`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/tests/test_sas.py:594>).
6. **유한 probe.** line별 two-level saturation, power broadening, longitudinal depletion, Gaussian radial average: [`power_broadened_profile_table`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/probe_saturation.py:95>), [`propagate_gaussian_spectrum`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/probe_saturation.py:166>). Rb D2는 software-reference scope. D1/Cs는 generalized estimate. Pump와 probe 동시 유한일 때 frozen-pump separable approximation; 비수동이면 residual 축소를 명시: [`sas.py:522`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/schemes/sas.py:522>).
7. **실험 workflow.** density/Doppler temperature 분리, paraffin population memory, absolute/relative CSV calibration, sweep branch 보존, custom isotope ratio, fixed-theory alignment/export.

## ⟦AI⟧ 실험 원자물리 연구자 판정

**실제 물리 있음. 조건 탐색용 반정량 reference로 유용. 절대 SAS metrology reference 아님.**

| 실험 질문 | 적합성 |
|---|---|
| isotope/hyperfine line assignment, scan span, 온도·길이·조성 경향 | 좋음 |
| pump-off weak-probe Rb OD, Doppler/Voigt scale | 유용. lifetime sum·소프트웨어 reference closure 강함. 독립 measured holdout 아님 |
| Rb D2 finite-probe trend·Gaussian depletion | 계산 reference로 유용. 실험 포화 검증 아님 |
| pumped SAS feature 위치·상대 경향 | 유용한 반정량 지도 |
| pumped SAS contrast/FWHM/lock slope·절대 pump power | 정량 인용 부적합. 기본 feature도 resolution-limited |
| high-pressure buffer cell·pressure-shifted lock point | 부적합 |
| Cs 절대 saturation/self-broadening | 부적합. Rb scalar self-broadening 상속: [`self_broadened_gamma`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/species.py:210>) |
| polarization·mF pumping·magnetic-field asymmetry·Faraday/dichroism | 모델 밖 |
| 데이터에서 T·L·조성 추정 | 현재 alignment 범위 밖. uncertainty-aware physical fit 없음 |

주요 정량 한계:

- Pump Ω는 모든 종·선에 고정 ⁸⁵Rb D1 `I_sat=4.484 mW/cm²` 사용: [`constants.py:30`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/constants.py:30>), [`pump_rabi_from_power`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/species.py:272>).
- Probe `Γ_eff = Γ_nat + self + Ne`와 달리 pump Liouvillian은 같은 elastic optical dephasing을 일관되게 받지 않음: [`sas.py:442`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/schemes/sas.py:442>), [`sas.py:456`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/schemes/sas.py:456>), [`sas.py:629`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/schemes/sas.py:629>).
- Pump는 입사 power로 한 번 정해짐. z별 pump depletion·온도/자기장 gradient 없음.
- mF·polarization 미해결. `rabi_rel`은 실제 개별 Zeeman matrix element가 아니라 lumped manifold prescription: [`species.py:365`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/species.py:365>).
- 기본 SAS FWHM은 FWHM당 약 2.34 sample. 코드가 6-sample gate 전에 lock slope를 숨기는 것은 올바른 trust behavior: [`_subdoppler_metrics`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/schemes/sas.py:1328>), [`resolution test`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/tests/test_sas.py:406>).

## 현재 수치 snapshot

입력·코드 hash가 2026-10-10과 같아 이전 측정을 재사용한다. natural Rb D1, 40 °C, 75 mm, pump 0.5 mW, probe 0 µW, waist 1 mm, `γ_t/2π=100 kHz`, 1,401 points.

| 경로 | 결과 |
|---|---:|
| 기본 SAS compute | 0.7825 s median |
| SAS `Tmin`, peak `αL` | 0.1180495075, 2.136651188 |
| 표시 FWHM | 18.86836 MHz, 2.3418 samples/FWHM, `resolution-limited` |
| pump-off OD compute | 0.3707 s median |
| pump-off `Tmin`, peak `αL` | 0.1050956546, 2.252884348 |
| custom ⁸⁵Rb fraction 0.35 | 0.7025 s median |

원 기록: [`2026-10-10 report:103`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/docs/daily_report/2026-10-10_scheme-1_od-sas.md:103>). 기본 FWHM을 실험 hero 값으로 인용 금지.

## ⟦AI⟧ 동작 불변 순수 코드 최적화

1. **Pump-off scalar fast path.** 현재 `(1.0 * Lp) @ wt`: [`_velocity_correlated_average`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/schemes/sas.py:1070>). `Lp @ wt` 직접 분기. 2026-10-10 current-hash benchmark `6.4255 → 0.8587 ms`, 7.48×, `array_equal=True`.
2. **Pump-off finite-probe 중복 matvec 제거.** `transition_alpha` 계산 직후 같은 `Lp @ wt`를 `thermal_alpha`로 재계산: [`sas.py:690`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/schemes/sas.py:690>), [`sas.py:694`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/schemes/sas.py:694>). Pump off면 `thermal_alpha = transition_alpha` exact.
3. **Dense fused kernel.** `deff_grid`, level별 `pop_at`, `arg/denom/Lp/Dp`, `np.interp`를 velocity-order 고정 compiled loop로 융합. 1 m/s·4σ discrete measure, Γ/8 pump table, interpolation, transition sum order를 regression gate로 고정. Roundoff tolerance 안 behavior 보존 가능. Pole/continuous Maxwell 변환은 pure optimization으로 분류 금지.
4. **Legendre node cache.** 매 finite-probe readout의 `leggauss(32)` 재생성: [`probe_saturation.py:195`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/probe_saturation.py:195>). order-key immutable cache. 수식·node 불변. 작은 readout 최적화.
5. **Frequency inversion loop compile.** frequency마다 `np.interp`를 도는 Python loop: [`probe_saturation.py:236`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/probe_saturation.py:236>). 같은 monotone table·연산 순서의 compiled loop. finite-probe batch에서 가치 큼.

`build_manifold`는 이미 bounded LRU cache: [`species.py:336`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/species.py:336>). 추가 Wigner cache 우선순위 낮음. Alignment coarse search 축소는 acceptance behavior를 바꿀 수 있어 profile·gate 없이 pure optimization 아님.

## 문서·예제 감사

좋은 점:

- README는 finite-probe scope, custom isotope fraction, alignment가 physical parameter fit이 아님을 설명: [`README:102`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/README.md:102>), [`README:129`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/README.md:129>).
- Tests가 실행 가능한 예제 역할. natural/custom fraction exactness: [`test_od_isotope_ratio.py:21`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/tests/test_od_isotope_ratio.py:21>). seconds-axis CSV fit와 physics 불변: [`test_od_comparison_ui.py:137`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/tests/test_od_comparison_ui.py:137>).

남은 문제:

- `SASScheme.info()`는 아직 “lab AutoOD calculator to <0.1%”라 표현. 현재 소프트웨어 closure와 독립 실험 검증 경계가 흐림: [`sas.py:352`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/schemes/sas.py:352>).
- User Guide도 같은 0.1% 표현 반복: [`Guide:409`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/docs/Userguide/GABES_User_Guide_v2.html:409>), [`Guide:922`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/docs/Userguide/GABES_User_Guide_v2.html:922>).
- Guide는 finite-probe를 natural Rb D2만 가능하다고 쓰고 기본을 `0 mW`로 표기. 실제 UI는 모든 built-in D1/D2·Rb/Cs에 0–100 µW를 노출하고 tier를 표시: [`Guide:556`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/docs/Userguide/GABES_User_Guide_v2.html:556>), [`Guide:573`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/docs/Userguide/GABES_User_Guide_v2.html:573>), [`sas.py:236`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/schemes/sas.py:236>).
- Guide에 custom isotope fraction·Fit alignment workflow 없음.
- 정적 `od.png`, `sas.png`는 2026-06-08 생성. 현 normalization·controls·readout 미반영.
- standalone `examples/` 디렉터리 없음. 재현 가능한 Scheme 1 script/config 없음. README·tests가 de facto examples.

문서 동기화는 solve 0. 실험 reference 신뢰를 즉시 올리는 negligible-overhead 작업.

## 테스트

### Scheme 1 집중 회귀

```text
python -m pytest -q tests/test_sas.py tests/test_absorption.py \
  tests/test_probe_saturation.py tests/test_experimental_csv.py \
  tests/test_experimental_alignment.py tests/test_od_isotope_ratio.py \
  tests/test_lineshape.py tests/test_schemes_render.py --tb=short
208 passed, 1 skipped in 94.92 s

python -m pytest -q tests/test_od_comparison_ui.py --tb=short
8 passed in 18.94 s
```

합계 **216 passed, 1 skipped**. Scheme 1 assertion 실패 없음.

Skip 1건은 [`test_reference_od_full_file_is_imported_and_smoothed`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/tests/test_experimental_csv.py:473>). 테스트는 없는 `references/AutoOD/ReferenceOD.csv`를 찾고, 실제 파일은 [`archive/AutoOD/ReferenceOD.csv`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/archive/AutoOD/ReferenceOD.csv>)에 있음. 물리 실패는 아니나 reference CSV importer coverage가 조용히 비활성. 권위 source 경로를 결정해 skip을 복구해야 함. 추가 solve·runtime 없음.

### 저장소 전체 회귀

```text
python -m pytest -q
2255 passed, 4 skipped, 2 failed in 2524.85 s (0:42:04)
```

실패 2건은 Scheme 1 집중 회귀 밖이다.

- [`test_empty_campaign_extension_passes_both_real_capsules_without_any_native_solve`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/tests/quantum/test_thermal_campaign_extend.py:357>): Windows 임시 경로 아래 Grand Challenge fresh-capsule 추출 중 중첩 목적지의 부모 디렉터리가 없어 `FileNotFoundError`.
- [`test_hidden_noise_sliders_retain_values_across_streamlit_reruns`](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/tests/test_fwm_excess_noise.py:470>): Scheme 5 FWM의 숨겨진 noise slider를 다시 보일 때 metric이 기대값 `0.04` 대신 이전 constant-mode 값 `2.0`을 유지.

따라서 전체 suite는 green이 아니지만, traceback과 위 집중 회귀 결과상 두 실패 모두 오늘 검토한 Scheme 1 물리 계산·CSV 비교 경로와 독립이다.

## ⟦AI⟧ 최종 우선순위

1. **P0 재현성:** 현재 thermal normalization correction과 독립 tests를 clean checkout에서 재현 가능하게 정리. 현재 live 결과와 HEAD 결과가 다름. 이동 뒤 끊긴 `ReferenceOD.csv` test 경로도 복구.
2. **P1, runtime 거의 0:** `(gas, species, line)` collision/pressure-shift/elastic-dephasing table. Pump/probe homogeneous-width convention 통일.
3. **P1, O(1):** species/line/polarization별 pump Rabi·saturation provenance. Pump-off OD exact 보존.
4. **P1 trust, solve 0:** isotope checklist를 “Rb fraction done / 새 원소 needs decision”로 분리. `SASScheme.info()`·Guide·그림·단위·workflow 동기화. standalone example 추가.
5. **성능:** scalar fast path와 중복 matvec부터. 다음 discrete-measure-preserving fused kernel.
6. **별도 opt-in research:** fitting, segmented pump depletion, full VCC, Zeeman/polarization. 측정 target·holdout·수렴 budget 전 기본 경로에 혼합 금지.

## 결론

OD/SAS는 실제 원자물리를 구현한다. Hyperfine 에너지·세기, thermal Doppler, same-class pump/probe 상관, optical pumping, Beer–Lambert/Gaussian 전파가 연결됨. Line assignment, cell 조건 탐색, 스펙트럼 해석에 유용한 반정량 reference다.

가장 큰 즉시 문제는 solver 복잡도 부족보다 **재현성·계수 provenance**다. 현재 올바른 thermal normalization은 dirty tree에만 있고, clean HEAD는 다른 물리를 계산한다. 그다음 위험은 fixed Rb-D1 pump Rabi, pressure shift 부재, pump/probe dephasing 불일치, mF/polarization 부재다. 이 중 collision·Rabi·문서 교정은 계산시간을 거의 늘리지 않고 신뢰도를 크게 올릴 수 있다.

⟦AI⟧ Head Manager closure: 이번 작업은 Scheme 1 유지보수 감사. Grand Challenge의 현행 FWM thermal 카드·공식 milestone 0/4·다음 질문을 바꾸는 새 증거 없음. Blueprint·registry의 Grand Challenge 상태는 수정하지 않음.
