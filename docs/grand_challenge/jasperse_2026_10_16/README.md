# Jasperse 식·NF·seed 잡음·효율 검수 — 2026-10-16

## ⟦AI⟧ 사용자 직접 요청 요약

원문 요청은 [decision ledger](../head_manager/decisions.jsonl)의
`REQ-2026-10-16-JASPERSE`, `GAIN-ROLLBACK`, `CONJUGATE-EFFICIENCY`에 분리 보존했다.

첨부 Jasperse 학위논문/논문으로 식·기호를 재검증한다. 선형 잡음비는 NF,
로그 표시는 Squeezing Degree [dB]로 통일한다. 고정 N을 제거하고 seed의
고전 잡음 a와 계산된 shot noise b, 계산된 양팔 total 효율을 표시한다.
잘못된 자료에 의한 semi-empirical gain correction은 rollback하되 개발 방법은
보존한다. 다음 턴 사용자가 새 데이터를 제공할 때 이 방법을 다시 사용한다.
순증폭 때문에 conjugate 효율을 표시하지 않는 선택지 대신, Jasperse 효율에
대응할 물리량을 GABES에서 추출할 방법을 검토하라는 추가 요청을 반영했다.

## ⟦AI⟧ 검증 결과와 적용 범위

원문 기준은 [Jasperse·Turner·Scholten, Optics Express 19, 3765–3774 (2011)](https://doi.org/10.1364/OE.19.003765)과
첨부 2010 학위논문이다. [저자 preprint](https://arxiv.org/abs/1012.3482)는
2010 제출·2011 개정이며, 학술지 출판 연도와 구분했다.

**흡수와 방출을 분리하면 conjugate 효율을 계산할 수 있다.** 기존 retarded
감수율은 두 과정의 차이만 제공한다. 같은 reduced Lindblad 모델의 finite-seed
주기 상태에 ordered QRT를 적용하여 양의 흡수율 α_abs와 방출률 α_em [m⁻¹]을
별도로 얻었다. `T_abs=exp(−α_abs L)`를 cell 효율로 사용하고,
`T_net=exp[−(α_abs−α_em)L]`는 별도 진단으로 보존한다. 순증폭을 1로 자르거나
투과율의 역수로 바꾸지 않는다.

조건부 장치점(blue Δ=1 GHz, 113 °C, pump 294.8 mW, 가정 seed 3.918 µW,
L=12.5 mm, lab δ=−23.7 MHz)에서 probe 흡수/방출은
126.97683/102.10175 m⁻¹, conjugate는 15.375923/19.193586 m⁻¹이다.
계산된 흡수 투과율은 각각 0.2044953, 0.8251426이다. 후단과 QE의 곱을
각 값에 곱해 sidebar의 total 효율을 표시한다.

이 대응은 **흡수 채널의 survival 효율**이다. 방출이 존재하는 셀 전체를
Jasperse의 순수 손실 채널 하나로 대체하는 것은 아니다. 현재 QRT 진단은
균일한 지정 pump/seed, zero coherent conjugate의 reduced four-level 상태와
절단 1D Maxwell 평균이다. 전체 spatial/depletion/Zeeman/RF covariance 전파는
포함하지 않는다. 미시 NF와 실제 실험 squeezing의 인증은 열린 과제로 유지한다.
QRT는 공통 probe-seeded 상태의 미분 응답이다. 서로 다른 seed 상태의 유한 χ
비율로 구성한 기존 transfer와 완전히 같은 행렬이라고 가정하지 않는다.
기존 weak-field normalization의 `0.74/12`는 유지한 조건부 입력이며,
이번에 rollback한 `0.5594938027` coupling multiplier·local affine fit과 구분한다.
따라서 이번 rollback 자체가 완전한 parameter-free 예측을 인증하지는 않는다.
[독립 QRT 유도](qrt_absorption_derivation.md), [구현 검수](qrt_absorption_review.md).

| 항목 | Jasperse와의 대조·현재 규약 |
|---|---|
| NF | 2011 Eq.3의 `Var(N_probe−N_conjugate)/mean(N_probe+N_conjugate)` |
| Squeezing Degree [dB] | `10 log10(NF)`; SQL=0 dB, 음수는 SQL 이하 |
| intrinsic gain | Jasperse `G=cosh²S`; 이 S는 squeeze parameter이며 NF를 뜻하지 않음 |
| a, b | UI의 명시적 규약: one-sided classical RIN a, shot RIN `b=2ℏω_seed/P_seed`, 둘 다 Hz⁻¹ |
| seed Fano factor | `F=1+a/b`; optical PSD는 `aP²`, `2ℏωP` [W²/Hz]; 적분 분산에는 대역폭 필요 |
| Jasperse η_a,η_b | 셀 이후 손실; GABES는 이를 after-cell optical transmission×QE로 전달 |
| cell 및 total 효율 | `η_in cell=T_abs`; `η_probe/conjugate=T_abs×η_after cell×η_QE` |

Jasperse의 원문은 lowercase a,b를 이 noise 규약으로 유일하게 정의하지 않는다.
학위논문의 검출기 보정 `AP²+BP+C`와 구분해 PSD/단위를 명시했다.
비대칭 후단 손실은 2011 Eq.5와 독립 공분산으로 검산했다. 기존 L1은
분포 probe 손실·conjugate 무손실일 때 Eq.16과 24개 독립 조건에서
최대 3.55×10⁻¹⁵ 차이로 일치한다. 일반 conjugate 분포 잡음은 Eq.17의
creation-reservoir 항이 필요하므로 L1과 동일하다고 주장하지 않는다.

2010 학위논문에는 covariance 재배열 부호(p.20), cosh/sinh 급수 부호(p.27),
Eq.3.45의 항/계수(p.31)에 차이가 있다. 원문을 시각 대조했으며, 독립 재유도와
일치하는 2011 Eq.16/17을 참조 구현으로 사용한다. 공식 erratum 확인 주장은 아니다.
[원문·기호·수식 상세 대조](math_reference.md).

셀 응답은 NF에 사용하는 원자 transfer에 이미 들어 있다. sidebar의 total 효율을
NF에 다시 곱하지 않고, 후단 효율과 QE만 별도 적용한다. 검출 효율이 0이면
분모인 detected SQL도 0이므로 NF는 unavailable이다. 형식적인 η→0 한계값과
실제 무신호 판독을 구분한다.

## ⟦AI⟧ 코드·화면·GC 검수

- 현재 화면에서 a 아래 b, after-cell/QE 아래 total 효율을 조작 불가능한
  슬라이더로 표시한다. 계산값은 solver 입력/cache key에 포함하지 않는다.
  NF는 선형 값에서 보간한 후 dB로 변환한다. export와 SABES도 현재 NF를 사용한다.
  [실제 화면](../../../output/playwright/jasperse_20261016_qe80.png),
  [값 변경·조작 불가 검증](ui_validation.json).
- 현행 모든 tier에서 기존 coupling multiplier와 local gain fit을 적용하지 않는다.
  역사 재현은 `legacy_gain_calibration_fixture=True`와 기존 enable flag를 함께
  지정할 때만 가능하다. 과거 방법·자료·fixture·허용 오차는 보존했다.
  **다음 턴 새 데이터가 필요하다는 사용자 요청을 인계에 기록했다.** 새 fit은 0회다.
- 독립 ordered-engine/fast QRT 9점 최대 차이 6.95×10⁻¹²,
  two-level analytic control 2.20×10⁻¹⁵, 생산 dissipator parity 1.44×10⁻¹⁶.
  Floquet refinement와 적응 Maxwell 오차도 별도로 검사했다.
- GC 21진입점·34 FWM dependency binding·8 observable node·94 수치 source,
  4 native 보고서의 전체 수치 payload와 37봉인 artifact를 대조했다.
  마지막 표시 경계 수정까지 반영한 [영향 검수](release_gc_impact_comparison.json)는
  PASS다. 새 thermal transport 계산 0, 공식 milestone 0/4 유지.
- 전체 `python -m pytest -q`: **2392통과·4skip·0실패**, 868.00초.
  [검증 기록](code_validation.json), [최종 receipt](release_full_pytest_receipt.json).
  Python263개 실행 전후 불변, [37봉인 자료·GC 동일 소스 최종 대조](post_validation_integrity.json) PASS.
  최초 전체 실행의8실패는 이전 UI/검출계약 기대값이어서 두 테스트 파일만 이관했다.
  최초 로그는 보존했고, 생산 Python 변경 없이 최종 전체 회귀를 통과했다.
  새 현재 UI 의존성이 된 SciPy도 requirements에 명시했다.

과거 결과와 GC의 미수렴 v2는 보존한다. 이번 결과는 평균장 절대 이득의 기존
잔차를 해결한 것으로 승격하지 않는다. [이전 이득 기준면 검토](../mean_gain_2026_10_15/README.md).
