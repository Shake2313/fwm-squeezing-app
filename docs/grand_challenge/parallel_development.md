# Grand Challenge 병렬 개발 — RF 대역과 입력 불확도

2026-09-11. 두 위상 원자 모델을 실제 probe/conjugate optical modes에 연결하는
작업과 독립적으로 구현한 downstream 연구 도구다. 현재 branch에서 신규 파일만
추가했고, 기존 spatial/kinetic/periodic solver나 production UI를 변경하지 않았다.
Grand Challenge 전체와 absolute experimental prediction은 계속 미완료다.

## 구현한 범위

- [`spectrum_analysis.py`](../../gabes/quantum/spectrum_analysis.py): 명시한 RF 축의
  선형 PSD/SQL 비율에서 SQL 또는 target-dB 아래인 모든 연결 구간을 찾는다.
  서로 떨어진 대역, threshold equality, invalid/missing/unconverged 구간을 구분한다.
  양쪽 threshold crossing이 있는 구간에만 유한 `width_hz`를 제공한다. 계산 범위
  끝이나 누락 구간에서 잘린 대역은 `width_hz=None`이며 관측된 span을 별도 제공한다.
  같은 연속 spectrum provider를 다시 평가하여 crossing bracket을 세분화할 수 있다.
- [`uncertainty.py`](../../gabes/quantum/uncertainty.py): 입력의 값·단위·표준불확도와
  명시한 joint correlation/covariance를 받아 상관 Gaussian 표본과 다변량 1차
  불확도 전파를 제공한다. 입력들이 서로 다른 SI 크기를 가져도 dimensionless
  correlation에서 검증한다. 실패한 표본은 버리거나 다시 뽑지 않는다.
- [`downstream_audit.py`](../../analysis/grand_challenge/downstream_audit.py): 독립
  해석식으로 crossing과 출력 covariance를 확인하고, 기존
  `power_normalized_readout`의 gain·S₋ 계산에 입력 불확도를 연결한다.

이 작업은 checklist의 `fwm-rf-squeezing-spectrum-and-bandwidth` 중 독립적으로
허용된 분석·검증 부분과 Grand Challenge milestone 2의 입력 불확도 부분을
진전시킨다. 실측 detector/SQL calibration, 일반적인 band topology 변화의 통계,
UI 연결, held-out 실험 검증까지 완료한 것은 아니다. 다른 task가 갱신 중인
`checklist.json`과 `research_log.md`는 이 작업에서 수정하지 않는다.

## 검증 결과

최종 불변 보고서: [`s2_downstream_report_v2.json`](s2_downstream_report_v2.json).
확인용 그림: [`s2_downstream_v2.png`](s2_downstream_v2.png).
초기 `s2_downstream_report.json`과 `s2_downstream.png`도 보존한다. v2는 독립
검토 후 극단적인 PSD 비율의 보간 overflow와 완전히 상관된 입력의 numerical
null-mode 처리를 보완한 source snapshot이다. 이전 보고서를 덮어쓰지 않았다.

해석 spectrum은 `R(f)=0.25+((f-c)/b)^2`, `c=2.1 MHz`, `b=1.3 MHz`다.
SQL crossing은 정확히 `c ± b sqrt(0.75)`다. 81→321→1281점 RF grid와 별도로
local bisection을 검사했다. 0.01 Hz bracket tolerance에서 최대 crossing 오차는
약 **0.00285 Hz**였다. 이는 이 연속 해석 fixture의 root-finding 정확도다.
좁은 미관측 spectral feature나 atomic covariance 수렴을 보장하지 않는다.

`c,b`의 표준불확도를 20 kHz, 50 kHz, correlation을 0.6으로 선언했다.
독립적인 해석 Jacobian과 비교한 대역 경계/폭 covariance 상대 오차는
**8.56×10⁻⁹**다. 고정 seed의 512회 Gaussian 전파는 해석 covariance와 약
**0.62%** 차이를 보이며 동일 환경의 재실행 결과가 일치한다. 이 finite-sample
비교는 일반적인 Monte Carlo convergence 또는 coverage 보증이 아니다.

기존 reduced atom에는 pump power 600 mW와 waist 530 µm를 사용했다.
이 두 입력에만 각각 1% 표준불확도와 correlation 0.6을 가정했고, 나머지
소비 입력은 보고서의 ledger 값에 고정했다. 해당 불확도 값은 측정치가 아니다.

| 출력 | nominal | 조건부 표준불확도 |
|---|---:|---:|
| Probe power gain | 1.11065587 | 0.00081043 |
| Conjugate power gain | 0.11356221 | 0.00080062 |
| S₋ / SQL, 0.1 MHz | 0.84519784 | 0.00093337 |
| S₋ / SQL, 1 MHz | 0.84522156 | 0.00093386 |
| S₋ / SQL, 4 MHz | 0.84558329 | 0.00094133 |

미분 step을 절반으로 줄였을 때 output covariance 상대 변화는 **6.72×10⁻⁷**다.
이 결과는 single-velocity, prescribed pump, weak-field/bright approximation이며
진행 중인 일반 optical-mode 모델의 새 결과로 해석하지 않는다. RF 세 점은
모두 SQL 아래지만 양끝 crossing을 보지 못했으므로 intrinsic bandwidth는
`None`이다. Source와 detector의 bandwidth는 서로 다른 측정 평면으로 기록된다.

동일 electronic low-pass를 numerator PSD와 SQL에 함께 적용한 대조군은
정규화 비율과 대역폭을 보존한다. Flat sub-SQL 대조군은 전체 계산 범위의
span만 반환하며 유한 intrinsic width를 만들지 않는다. Target dataset을 이용한
입력·correlation evidence, 누락 evidence, 비물리 covariance/PSD, 다른 provider를
사용한 crossing refinement와 artifact 덮어쓰기는 거부 또는 실패로 남긴다.

## 재현과 연결

저장소 root에서 실행한다. 출력 경로는 새 이름이어야 한다.

```powershell
python -m pytest -q tests/quantum/test_spectrum_analysis.py tests/quantum/test_uncertainty.py tests/quantum/test_downstream_audit.py
python -m analysis.grand_challenge.downstream_audit --output NEW.json --plot NEW.png
python -m pytest -q
```

기존 source/readout 결과의 연결 예:

```python
from gabes.quantum.spectrum_analysis import SpectralTrace, analyze_squeezing_bands

trace = SpectralTrace.from_intensity_difference(result.spectrum, contribution="quantum")
sql_bands = analyze_squeezing_bands(trace)
target_bands = analyze_squeezing_bands(trace, threshold_db=-3.0)
```

`layer="source"`는 해당 spectrum 자체가 cell-output 측정 평면에서 계산되었을
때만 지정한다. Detector loss를 적용한 spectrum의 label을 바꿔 source spectrum을
얻을 수 없다. Uncertainty callback은 같은 input draw로 모든 gain·RF 출력값을
함께 반환해야 상관관계가 보존된다. 이번 band-edge uncertainty 예제는 하나의
유한 대역이 유지되는 조건이며 censored 또는 분리된 대역이 나오면 계산을 중단한다.

다변량 covariance 및 joint-distribution propagation의 배경은
[JCGM 102:2011](https://www.bipm.org/en/doi/10.59161/jcgm102-2011)이다.
Gaussian 입력분포, correlation과 fixed-input 가정의 실험적 적절성은 별도 근거가
필요하다. 구현 상세는 [RF 대역 유도](spectrum_analysis_derivation.md)와
[입력 불확도 유도](uncertainty_derivation.md)에 둔다.

## 자동 검사

신규 검사 **50 passed**: RF 대역 16개, 입력 불확도 32개, 통합 보고서 2개다.
극단적인 finite PSD/SQL 비율 `1e308`에서도 교차점은 bracket 내부에 머문다.
완전히 상관된 5개 입력의 공통 오차는 동일 표본을 공유하며, 차이 신호에
`10^9`의 Jacobian을 적용해도 가짜 분산이 생기지 않는다. `ρ=1−10⁻¹²`처럼
정확히 1이 아닌 상관의 실제 양의 분산은 유지한다.

v2 보고서의 선언한 controls와 실행 중 source stability가 모두 통과했으며,
기록한 source SHA-256 **25개를 재확인**했다. 그림도 렌더링하여 확인했다.

최종 필수 `python -m pytest -q`: **855 passed, 1 failed (204.71 s)**.
동시에 개발 중인 optical-mode task의 검사도 이 전체 실행에 포함되어 있다.
유일한 실패는 작업 전부터 존재한
`tests/test_docs_consistency.py::test_repository_visibility_wording_is_consistently_public`의
`docs/FWM physics and analytic reconstruction/FWM_physics.tex` 누락이다.
해당 삭제나 문서 검사는 변경하지 않았다.
