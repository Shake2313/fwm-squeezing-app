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
