# 실험용 FWM gain 보정 — 2026-10-22

사용자 요청: 10월13일 측정 기반 보정토글, 피크위치 보정, 측정범위밖·운전조건변경에도 계속계산. GC와분리.

## ⟦AI⟧ 방법·적용범위

원시23점: `gabes.fwm_measured_gain.observations()`. 펌프ON양팔전압/펌프OFFseed35mV.
기준면은검출기앞비율. 원자매질경계gain·잡음공분산의보정으로자동해석하지않음.
10/13측정·10/16전달 workbook SHA256: `2653fe73abca75550b2b8665da8c8fff16f5c28f155e029f1e6da155a9409ed0`.

명목계산: OPD+1GHz,113°C,pump294.804mW,seedproxy3.921545933µW,L12.5mm,waists530/330µm,angle0.32°,Balanced401점.
계산피크δ24.75MHz; 실측probeplateauδ1.732439–2.732439MHz,conjugateδ2.732439MHz.

팔별주파수이동/폭·양의아핀gain식부터log잔차최소화. 아핀식만의최대상대오차probe25.55%/conjugate28.77%.
잔여shape차이는동일23점의boundedlog보정계수로반영. 런타임식:

\[
x_j=c_j+w_j(\delta-d_j),\qquad
R_j(\delta,p)=G_j^{model}(x_j,p)C_j(\delta),\quad
C_j(\delta_i)=R_{j,i}^{measured}/G_j^{nominal}(x_{j,i}).
\]

`C_j`는측정영역안log-PCHIP,밖은경계보정계수고정. **gain값고정이아님**: 현재운전조건모델을계속곱함.
아핀offset도명목계산과현재계산의세기비로전달하여,원래conjugate가0이면보정뒤에도0.
원자계산scan의밖에서는양의지수연장으로probe1/conjugate0에접근. 이추가연장범위도별도표시.
다른온도·pump·seed·OPD·waist·길이·각도·branch는상태만extrapolation으로바꾸며토글차단없음.
검출효율·고전잡음조절로재피팅하지않음. canonicalM·원자raw배열·GC계산수정없음.

## ⟦AI⟧ 검증·재현

- 명목모델의23점양팔비율:최대상대오차3.4×10⁻¹⁶. 동일자료fit이며독립예측검증아님.
- 0.005MHz격자보정피크:probeδ2.170MHz/25.0164,conjugateδ2.880MHz/25.4719. 측정점사이작은peak차이는보간·모델곱의결과.
- `tests/test_fwm_experiment_gain.py`:17통과. 범위밖±10⁶MHz유한성,전계수상태유지,운전조건변경9종,모델응답변화,무신호conjugate보존,OFF정확복원·raw불변·현행명목/변경조건실제원자solve확인.
- [재현script](fit_gain_profile.py), [전체fit기록·명목모델](gain_profile.json).
  `python analysis/fwm_experiment_20261022/fit_gain_profile.py`로기록재생성.
- embedded런타임계수는`gabes/fwm_experiment_gain.py`에저장. 실제운전중refit·분석폴더파일읽기없음.
