# ⟦AI⟧ 실측 이득 참조 구현 독립 검수

2026-10-16. 채택 상태: 상위 작업자가 채택한 **23점 실측 참조 보간**의 구현 검수. 독립 예측 이득 모델 또는 원자·잡음 모델의 재보정 승인이 아니다. 검수자는 Python 소스를 수정하지 않았다.

## 판정

검수 범위에서 남은 조치 사항 없음. 새 참조는 두 팔 검출기 전압을 공통 pump-off seed 35 mV로 나눈 비율이다. 측정 기준면에서 `.8046`·QE를 다시 곱하지 않는다. `.8046`은 창을 포함한 전체 셀 투과율 추정치이므로 atomic $A_0$로 채택하지 않는다. 비교 preset의 seed **3.921545933233986 µW**는 $P_{off}/(0.8046\times0.9709)$ 외부 입사 proxy를 원자 입력으로 가정한 값이며, 독립 측정된 vapor 입력이 아니다.

- **부호:** lower seed 규약에서 $f_{EOM}=\nu_{HF}-\delta_{internal}$; 정확한 $\nu_{HF}=3035.732439$ MHz 사용. Spreadsheet의 $TPD=f_{EOM}-3036$과 $\delta_{internal}=-TPD-0.267561$ MHz로 연결된다. preset 3034 MHz는 internal $+1.732439$ MHz.
- **적용 범위:** 기본 OFF. 명시적 proxy 가정·명목 장치 조건·lower branch가 일치해야 적용된다. RF 3022–3044 MHz 밖, 다른 branch, 옛 gain fixture 조합은 사용 불가. 숫자 비교 허용폭은 단위 환산 부동소수점 오차용이며 실험 불확도가 아니다.
- **경계 수정 확인:** 실제 광 주파수 축 $D-f_{EOM}/1000$을 역변환하면 3044.0000000000005 MHz가 되어 끝점이 배제되던 문제가 있었다. 통합부의 양끝 16 ULP 이내 정규화 후 23점 모두 포함, 양팔 node 일치, scalar/array endpoint 일치가 확인됐다. 3044+$10^{-8}$ MHz는 계속 제외한다. 독립 helper의 엄격한 support·NaN 정책은 유지된다.
- **원자·잡음 분리:** `fwm.py`가 원자 $G_s,G_c$, canonical transfer, NF 계산 뒤 별도 `measured_gain_reference`를 추가한다. NF에는 원자 canonical map과 seed $a/b$, 셀 이후 손실·QE를 사용한다. 실측 참조 비율 또는 `.8046`을 NF에 대입하지 않아 셀 흡수를 중복 반영하지 않는다. SABES/Grand Challenge no-fit 계산으로 보간을 전달하는 경로도 추가하지 않았다.
- **표시·주장:** 자연로그 비율 PCHIP는 양팔 23점 재현과 구간 shape/no-overshoot 검사를 통과한다. Probe peak는 3033–3034 MHz plateau, conjugate는 3033 MHz이다. 전점 사용 후의 재현은 독립 validation이 아니다. 고정 교대 cohort의 사후 보간 진단 최대 상대오차는 probe **19.92247%**, conjugate **15.54499%**이며, `development_acceptance_passed=False`와 독립 예측 검증 없음이 함께 명시된다. 실패한 Chebyshev·spline 적합 기록은 보존된다.

## 검증과 동결 기준

직접 실행: `python -B -m pytest -q tests/test_fwm_measured_gain.py tests/test_fwm_measured_gain_integration.py` → **43 passed**, 기존 UTC export deprecation warning 1건. 이 검수는 전체 회귀를 대체하지 않으며 전체 회귀·최종 통합 판정은 상위 작업자 담당이다. UI 담당자가 별도로 157 targeted tests와 실제 Streamlit preset/checkbox rerun에서 proxy slider full precision 및 `context_applied=True`를 확인했다고 보고했다. AppTest의 기존 segmented-mode serializer 문제 때문에 해당 모드 위젯만 harness에서 제외된 범위다.

채택 참조 기록: `interpolation_record.json`, SHA256 `0a52cc123edff9749a9d1bcde0649e5e9cd6961961149304e9c165aab7a50b64`.

| 동결 파일 | SHA256 |
|---|---|
| `gabes/fwm_measured_gain.py` | `df35e67a034cb028bb69a4b81f7867fb2e36f6db84bfdf02947dd9e5e55b553c` |
| `gabes/schemes/fwm.py` | `d3e1d2399befe967357d0400a25ad4fd3a6b395a35c0bc48d91c617f396e5e98` |
| `tests/test_fwm_measured_gain.py` | `1a9f7bfe8c9aefccc00775a4fc9bb087964f81bdad6b7b3d28c6d754a24c461e` |
| `tests/test_fwm_measured_gain_integration.py` | `0f97434c835d91df496a29c47b64d0955e85e25b17bc94bb2eb87407af666553` |

이 SHA 집합에 대해 독립 검수 완료·소스 동결 확인. 이후 Python 변경이 있으면 해당 변경만 추가 검수한다.

## 분석 실행기 식별

`interpolation_record.json`의 `runner_sha256`는 최초 baseline/적합 실행기 `fit.py`가 아니라, 최종 보간 실행기 `interpolate.py`의 바이트를 가리킨다. 현재 파일을 직접 해시하여 기록값과 일치함을 확인했다. 사후 소스 변경에 따른 해시 불일치가 아니다. 보간 실행기는 `fit.py`의 공통 오차 지표 함수만 import하며, 새 원자 계산 없이 보존된 raw ratio에서 참조를 만든다.

| 실행기 | 현재 SHA256 |
|---|---|
| `analysis/fwm_gain_hotfix/calibration_20261016/fit.py` | `13eec2b4554572a890d1b935c62eca28d075b67744d1266009ba93fc71d57aef` |
| `analysis/fwm_gain_hotfix/calibration_20261016/spline_pivot.py` | `5a1610d609dee94da7e534c1c745357ff3e15fa6b7bbac06a86732db5aaffa8a` |
| `analysis/fwm_gain_hotfix/calibration_20261016/interpolate.py` | `4ae5fd6b63d24351776247fdbbcb5fa4fb072e64fe51244e4b4b576cb695b814` |
