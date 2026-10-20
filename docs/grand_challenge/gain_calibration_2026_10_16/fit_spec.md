# ⟦AI⟧ 새 측정 ratio 교정 사전 명세

선언: 2026-10-16, 새 적합 실행 전. 사용자 승인 범위: 새 workbook에 맞춘 재교정, 기존 개발 방법 보존, 두 팔·기준면·부호·Grand Challenge 영향 검수. 본 실행은 offline analysis만 소유하며 생산 코드/원자 모델/잡음 경로를 변경하지 않음.

입력: `D:/20261013_TPD_vs_Gain 및 스퀴징 특성 추정.xlsx`, SHA-256 `2653fe73abca75550b2b8665da8c8fff16f5c28f155e029f1e6da155a9409ed0`.
관측: Sheet1 EOM3022-3044MHz의 두 팔 전압23점 /35mV. 같은 converter·after-cell0.9709 사용은 사용자 확인. T113°C, L12.5mm, waists530/330µm, angle0.32° 확인. OPD+1.0GHz는 근사 운전값이며 정확도 별도. 0.8046은 windows 포함 full-cell transmission으로 해석; atomic $A_0$ 채택하지 않음.

## 모델과 데이터 분리

1. 현행 보정 OFF 원자 모델을 실제23개 주파수와 별도 집중 격자에서 계산. $\delta_{internal}=3035.732439-f_{EOM}$MHz. baseline 상대 순서/peak를 점검하고 old monotone hinge $F(G)$의 순서 역전 하한을 두 팔 각각 계산.
2. 비교용 seed proxy: $P_{off,ref}/(0.8046\times0.9709)$. 이는 **외부 입사면 추정치**를 atomic input처럼 사용한 조건부 비교이며 독립 vapor-input 교정이 아님. 이 proxy나 full-cell factor로 원자 이득을 측정 ratio라고 인증하지 않음. baseline fitting으로 물리 파라미터 변경하지 않음.
3. 새 교정은 측정 ratio 전용 global RF 모델:

   $$x=(f_{EOM}[MHz]-3033)/11,\qquad\log R_j=\sum_{k=0}^{d}c_{jk}T_k(x),\quad d\in\{3,5,7\}.$$

   degree별 finite global coefficient 모델. 점별 lookup/배율/원자 transfer 수정 아님. 양팔 계수 따로, 선택 기준 공통. Domain3022-3044MHz 밖으로 적용/검증 주장하지 않음.
4. 원래 OOXML row 순서(3044→3022)의 0-based even indices12개(양끝 포함)를 development train, odd indices11개를 development validation으로 고정. 새 데이터23점은 이미 검수에 모두 사용되어 **untouched holdout이 아님**.
5. 이전 방법의 선형계획 minimax 절차 보존. Log-ratio Chebyshev basis의 최대 absolute log residual을 LP로 최소화. 새 결과는 실제 ratio 상대오차 $|\hat R/R-1|$로 평가. 이는 log-linear global model의 multiplicative minimax이며 옛 원자 gain-only hinge의 수치 계수를 재사용하지 않음.
6. 양팔 모두 development validation 최대 상대오차≤20%인 가장 낮은 degree 선택. 어느 것도 통과하지 않으면 모든 후보 실패를 보존하고 양팔 joint validation 최대오차가 가장 작은 degree를 선택하되 acceptance 실패 scope flag 유지. 동일 점수이면 낮은 degree.
7. 선택된 degree를23점 전체에 offline refit. 전점 상대 max/RMS·ratio-unit RMSE·측정 peak와 연속 domain 모델 peak를 두 팔 모두 기록. 이 refit 결과는 training 진단이며 confidence interval/독립 검증 아님.

## 고정 제약·출력

- 원시 workbook·기존 Gold/TPD calibration 증거 bytes 보존. 이전 `--fit/--capture-before` runner 호출 금지.
- 새 runner `analysis/fwm_gain_hotfix/calibration_20261016/fit.py`; 결과 이 폴더의 `fit.json`, baseline/candidate 근거.
- 새 fit는 $R_j=V_{j,on}/V_{seed,off}$에만 해당. Atomic Gs/Gc, canonical M, NF, dressed absorption 및 SABES/GC no-fit 경로로 계수 주입 금지.
- 원자 baseline의 source SHA를 실행 전후 기록하고 일치 확인. 변경 중이면 baseline을 재사용하지 않고 원인 보고.
- 새 atomic ensemble solve는 baseline 비교용2 scans(실측23점+집중 격자)로 제한. Chebyshev/LP 후보는 원자 solve 추가0.
- 다음 독립 관문: 새 장치 조건·반복 측정과 실제 seed 입력·전압 응답/기준면 확인. 이번 global measured-ratio fit는 그 관문을 대신하지 않음.

## ⟦AI⟧ 사전 선언한 bounded basis pivot

최초 Chebyshev3/5/7의 양팔 joint development validation 최대오차는98.803%/33.381%/33.156%. 모두 사전20% 관문 실패. 최초 `fit.json`과 계수·실패 판정은 역사 증거로 유지하고 덮어쓰지 않음. 이 절은 새 spline fitting 실행 전에 추가.

새 가설: 전구간 고차 다항식의 형상 제약이 steep asymmetric tail을 충분히 표현하지 못함. 더 높은 polynomial degree·점별 lookup·peak 강제 없이 **piecewise cubic B-spline log-ratio**로 basis만 바꿔 판별.

- 후보3개: domain3022-3044MHz에 균일 간격의 **interior knot2/3/4개**, clamped cubic(끝점4중 반복). 각각 coefficient6/7/8개. 측정점마다 knot를 두지 않음.
- 기존 alternating-row train12/validation11과20% 양팔 최대 상대오차 관문 그대로. 팔별 natural-log minimax LP, 같은 source 데이터와 단위 사용. Holdout 이동·강제 peak 제약·추가 noise fitting 없음.
- 관문 통과 후보 중 interior knot 수가 가장 적은 모델 선택. 모두 실패하면 `selected_candidate=null`, 실패를 보고하고 생산 채택 근거로 삼지 않음. 후보별23점 refit는 설명용 max/RMS·연속 peak 비교에만 사용; 최종 선택을 training 잔차로 뒤집지 않음.
- 새 결과 `spline_fit.json`, 새 model ID. 최초 Chebyshev 파일은 그대로. 기존 atomic baseline 재사용, **새 atomic solve0**. 원래 spec의 hash는 최초 fit 기록에 보존하고 이번 파일은 append된 현행 명세로 관리.

## ⟦AI⟧ 사전 선언한 sampled reference 전환

Spline2/3/4 interior knot의 양팔 development validation 최대오차는28.187%/25.813%/30.482%로 모두20% 관문 실패. `spline_fit.json`의 미선택·실패 결과 보존. 새로운 예측 모델을 통과했다고 주장하지 않고, 상위 승인에 따라 **23개 원시 관측점의 결정론적 log-PCHIP interpolation**으로 용도 전환. 이 절은 새 interpolation record 계산 전에 선언.

- EOM3022-3044MHz의 각 팔 양의 ratio23개 모두 명시적 source point. Natural-log ratio에 PCHIP, 결과를 exponential로 환산. Free/fitted coefficient·forced peak·새 atomic solve 없음.
- 지원 범위는 관측 RF 구간만. 외삽은 unavailable. 원자 Gs/Gc·canonical M·NF를 대체하지 않고 측정 ratio의 sampled calibration/reference로 따로 표시.
- Exact-node 재현, 구간별 endpoint 범위/no-overshoot, 단조 구간 derivative sign, 양수, support 밖NaN을 확인. 연속 peak는 노드 최대값과 plateau 구간까지 보고; 별도 peak 위치 fitting 없음.
- 동일 고정 even12/odd11 development interpolation 진단도 보존. 실패해도 점 삭제·cohort 이동·예측 관문 완화 없음. 이 진단의 성패와 독립적으로 **predictive development acceptance는 이전 실패 상태** 유지. Sampled interpolation의 유효성만 별도로 판정.
- 새 `interpolation_record.json`과 scientific comparison/residual PNG, 새 model ID. 원래 workbook·두 실패 fit 기록·baseline bytes 그대로 유지. 구현 helper는 원시 source point와 source SHA를 고정하며 runtime fitting하지 않음.
