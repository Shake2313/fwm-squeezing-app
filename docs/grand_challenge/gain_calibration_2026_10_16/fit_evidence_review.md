# ⟦AI⟧ 새 ratio 교정 후보·기존 방법 한계 검수

2026-10-16. [사전 명세](fit_spec.md), [원시·최초 적합](fit.json), [spline pivot](spline_fit.json), [원자 baseline](uncorrected_baseline.json), [기준면·부호 감사](independent_plane_sign_review.md).

## 보존·원자 baseline

실측23점: Sheet1 두 팔 전압/35mV. old Gold $C_{mix}$·옛 TPD affine 모두 OFF. 기본 장치 조건과113°C·294.804mW 사용. Baseline은 actual23점과 internalδ±60MHz401점, 2scans; production Python source 실행 전후 SHA 일치. 원자 비교 실행 약0.73초. 추가 spline pivot의 atomic solve0.

비교용 seed proxy3.921545933µW는 $P_{off}/(0.8046×0.9709)$로 산출한 외부 입사 추정치를 원자 입력처럼 사용한 **조건부 baseline**. 0.8046을 atomic $A_0$라고 채택한 결과 아님. 실제 vapor 입력/TPD별 투과 미확정.

명목 원자 peak: probe236.27244@EOM3012.032439MHz, conjugate242.61058@3011.732439MHz. 측정3022-3044MHz 범위 밖. 측정 probe870/35=24.85714는3034/3033 두 점 plateau, conjugate890/35=25.42857는3033MHz.

**옛 gain-only monotone affine의 한계:** 고정 양의 기준면 환산을 가정한 순서 역전 하한은 probe99.8025%, conjugate98.5499% 최대 상대오차. 새 baseline에서 동일 LP hinge 방법을 직접 적합해도 각각 그 하한과 일치. 따라서 계수만 갱신해 이 RF 형상을 설명할 수 없음. TPD 의존 $A_0/\tau$가 미확정이므로 이 하한을 장치의 절대 원자 모델 오차로 주장하지 않음.

## 사전 고정 cohort와 판정

OOXML 순서3044→3022의 0-based even12점 train(양끝 포함), odd11점 development validation. 두 팔 모두 validation 최대 상대오차≤20%가 선택 관문. 모든 자료는 이미 검수에 사용되어 untouched holdout이 아님. 모두 자연로그 ratio minimax LP; noise·원자 계수 fit 없음.

| 후보 | 팔별 coefficient 수 | Probe validation max | Conjugate validation max | 양팔 관문 |
|---|---:|---:|---:|---|
| Global Chebyshev degree3 | 4 | 98.803% | 90.008% | FAIL |
| Global Chebyshev degree5 | 6 | 33.381% | 29.065% | FAIL |
| Global Chebyshev degree7 | 8 | 33.156% | 19.141% | FAIL |
| Cubic spline, uniform interior knot2 | 6 | 28.187% | 24.145% | FAIL |
| Cubic spline, uniform interior knot3 | 7 | 25.813% | 18.055% | FAIL |
| Cubic spline, uniform interior knot4 | 8 | 30.482% | 18.897% | FAIL |

최초 spec은 모두 실패 시 최소 joint validation 후보를 제한적으로 남기는 규칙이라 degree7을 `acceptance_passed=false`로 기록. 이후 사전 선언한 spline pivot은 모두 실패 시 미선택 규칙이므로 `selected_candidate=null`. 실패를 full-data training 잔차로 뒤집지 않음. 두 파일 모두 원문 유지.

## 설명용 full23 refit

| 후보 | Probe max/RMS relative | Conjugate max/RMS relative | Continuous peak EOM (P/C) |
|---|---|---|---|
| Chebyshev7 | 16.266% /11.762% | 8.843% /6.520% | 3032.132 /3032.298MHz |
| Cubic spline3 interior | 14.834% /10.626% | 11.765% /8.308% | 3032.565 /3033.046MHz |

전점 refit는 training 결과이며 validation 통과·실험 불확도·peak 위치 독립 검증이 아님. Coefficient는 RF 전구간 smooth function, 점별 lookup 아님. Domain 밖 extrapolation 금지.

실패를 구별하는 점: Probe validation3025MHz의 측정4.6mV에 Chebyshev7은6.125mV(+33.156%), spline3은5.787mV(+25.813%). 다음3023MHz의1.17mV는 각각−25.053%/−15.748%. 나머지 probe validation 오차는 Chebyshev7 약3.46% 이하, spline3 약5.76% 이하. 국소 basis는 최악 tail 오차를 줄였지만20% 관문에는 부족. 독립 전자 offset·반복 측정 불확도 없이 해당 두 점을 outlier로 삭제하거나 threshold를 완화할 근거 없음.

판정: **기존 gain-only 교정의 표현 실패 확인; 새 empirical RF 후보는 사전 development 관문 미통과.** 제한적 표시/추가 모델 채택은 상위 판단 사항. Atomic Gs/Gc, canonical M, NF, dressed absorption, SABES/Grand Challenge no-fit 경로를 바꿀 근거 아님. Raw experimental ratio와 empirical curve를 명시적으로 분리하는 계약 유지.

후속 AI 총괄 판단: 전역 예측 fit은 채택하지 않고23점 log-PCHIP 측정자료 보간으로 용도를 전환했다.
이후 [보간 기록](interpolation_record.json)과 [최종 검수](README.md)가 현행 채택 범위를 소유한다.
