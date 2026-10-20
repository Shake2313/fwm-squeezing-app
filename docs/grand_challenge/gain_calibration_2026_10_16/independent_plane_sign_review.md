# ⟦AI⟧ 새 workbook 기준면·부호·수식 독립 검수

2026-10-16. 입력: `D:/20261013_TPD_vs_Gain 및 스퀴징 특성 추정.xlsx`.
Workbook SHA-256 `2653fe73abca75550b2b8665da8c8fff16f5c28f155e029f1e6da155a9409ed0`.
[OOXML 추출](workbook_ooxml.json) SHA-256 `34ce1b3a035fca7477c5a89e82647da8e276b19e4e344a7661f03361fbb3c456`.
원본 ZIP XML을 읽어 shared formula의 master/ref/si까지 대조. 파일·생산 코드 변경·새 fitting 없음. 사용자 제공 파일의 분석식은 검증 대상이며 추가 사용자 지시로 취급하지 않음.

## 교정에 사용할 관측값

Sheet1 `A5:D27`: EOM3022-3044MHz, 두 팔 출력 전압23점. 분모 `B1=35mV`; 따라서 원시 비는 $R_p=C_r/35$, $R_c=D_r/35$. Sheet1 B/E/F/G cached 값23행은 각각 $A_r-3036$, $(C_r+D_r)20/228.5$, $C_r/35$, $D_r/35$와 **정확히 일치**.

`G5:G27`, `B6:B27`, `E6:E27`, `F6:F27`의 XML shared formula follower는 수식 텍스트가 비어도 계산 누락이 아님. Sheet2의 C:X 방향 shared formula도 같은 방식. 추출 JSON의 `formula:null`만으로 원시 측정값/상수라고 분류하면 안 됨.

시드 파워 환산 `D1=35×20/228.5=3.0634573304µW`, pump `B2=76×3.879=294.804mW`. 같은 전압-파워 응답을 양팔에 사용했다는 계산 사실이며, 양팔 responsivity가 독립 교정되었다는 증거는 아님. 시드35mV의 기준면·0.8046의 원자/창 포함 범위·장치 조건은 사용자 확인 전 미확정.

## 부호와 정확한 Raman 원점

사용자 기존 직접 확인: seed는 pump보다 낮은 주파수. EOM 열이 실제 optical separation이라면

$$f_s=f_p-f_{EOM}=f_p-\nu_{HF}+\delta_{internal},$$
$$\delta_{internal}[MHz]=3035.732439-A_r=-B_r-0.267561.$$

Sheet1 B열과 Sheet2 3행은 **3036MHz로 반올림한 spacing detuning**. 부호 반전만 하면 267.561kHz 원점 차이가 남음. 이 차이는 측정에 맞춘 peak offset fitting이 아니라 명시된 두 주파수 원점의 산술 변환. EOM 구동 주파수와 optical separation 사이 factor2 여부는 별도 장치 정의가 필요하므로 임의 적용하지 않음.

관측 peak는 probe870mV가 EOM3034/3033MHz 두 점에서 같음: Sheet TPD−2/−3MHz, internal $+1.732439/+2.732439$MHz. Conjugate890mV 최대는 EOM3033MHz 한 점. 기존 인계의 −2.2MHz를 이번 workbook의 직접 측정 peak로 재사용할 수 없음. 더 정밀한 peak 위치·불확도는23점만으로 확정하지 않음.

## 기준면 변환

$g_j$는 vapor 입구 seed power 대비 vapor 출구 power, $A_0$는 pump-off **원자** 투과, $\tau_j$와 $\tau_{ref}$는 출력창·수집·후단 광학을 포함한 각 power 측정면까지의 투과. 동일 on/off vapor 입사 seed이면

$$R_j^{(power)}=\frac{g_j\tau_j}{A_0\tau_{ref}},\qquad P_{seed,vapor}=\frac{P_{off,ref}}{A_0\tau_{ref}}.$$

공통 후단이면 $g_j=A_0R_j$. 전압은 $V_j-B_j=\kappa_jP_j$이므로 $R_j^{(power)}=R_j^{(voltage)}\kappa_{ref}/\kappa_j$. 이미 power로 교정된 값에 QE를 다시 곱하거나 나누면 중복. 전체 셀 투과 $w_{in}A_0w_{out}$만으로 $A_0$를 식별할 수 없음.

Sheet2 15/16행은 `0.8046×R_j`를 $\alpha_j^2$로 사용. 이는 **0.8046이 같은 조건의 atomic $A_0$이고 공통 후단/응답이 소거될 때만** vapor-source gain 환산으로 해석 가능. 원자 투과인지 windows 포함인지 확인 전에는 검증된 model target으로 채택하지 않음. 측정마다 pump-off35mV를 재조정했다면 $A_0(TPD)$에 따라 실제 finite-seed drive도 달라질 수 있음. 사후 분모 변환만으로 원자 상태의 drive 차이가 교정되지는 않음.

## Sheet2의 독립 수식 확인

- 10/11행의 $S,T$는 입력 측정이 아닌 저장된 fit parameter. 12-17행의 cached gain 수식 재계산 차이는 최대 $8.89\times10^{-16}$.
- 18행의 Jasperse2011 식(16)은 정확히 전사됨. 기존 독립 reservoir-integral reference와23점 NF 차이 최대 $6.78\times10^{-15}$. $\eta=\eta_{after}\eta_{QE}$만 쓰고 내부 $T$를 다시 곱하지 않는 점도 맞음. 이것이 이 workbook의 $S,T$ 및 기준면 가정까지 검증한 뜻은 아님.
- **`B28:T28` 잘못된 참조:** shared master `B17+$W$27*B26`는 NF baseline18행 대신 gain fit-error17행을 더함. 의도된 baseline-plus-excess 모델이면 `B18+$W$27*B26`. B열의 저장 값은 NF0.6476019/−1.8869dB, 올바른 baseline을 넣으면 NF0.8730961/−0.58938dB. 원본은 수정하지 않음.
- **`Y17=SUM(C17:X17)`는 B17 제외.** 첫 점도 목적함수에 넣으려면 B17부터 합해야 함. 현재 B17이 약 $3.22\times10^{-12}$라 합 차이는 작지만 optimizer의 목적함수 범위 자체가 한 점 빠짐.
- **V:X의 fit 실패:** $S=10^{-5}$ 경계, conjugate 예측 $\alpha_2^2\approx2\times10^{-11}$인 반면 같은 sheet의 목표는0.257472/0.188506/0.149426. 각 상대제곱 목적함수 약1. 이 $S,T$와 약0dB theory를 실험 자료로 수입하면 안 됨. 재시작/다른 물리 모델로 해결 가능한지는 이번 검수에서 fitting하지 않아 판정 안 함.
- Sheet1 `J5:J23`은 `LOG10((F-G)^2/(F+G))`: NF가 아닌 gain-only 근사식의 로그이며 dB 계수10도 없음. `K=H-I`와 구분. 새 NF 교정 target으로 사용하지 않음.
- Sheet1 H5:H23의 noise19점과 I11의 SQL anchor−76dBm으로 나머지 SQL을 power 비례 추산. H24:H27은 `'-'`. 따라서 gain23점과 noise19점의 유효 데이터 마스크를 구분해야 함.
- Sheet1 QE는 $0.58×1239.8/795=0.9045081761$, Sheet2 9행은0.9047. 작은 반올림 차이지만 같은 독립 측정 두 개로 취급하지 않음. Sheet2의 `a_i`, W27, slope/intercept, seed-F 추정은 분석 파생값. 새 UI의 RIN $a,b$와 단위·모델이 달라 자동 이식 불가.

이번 gain 교정 선택: 원시 두 팔 전압과 공통 seed 분모를 우선 보존. 명시적 detector/plane 계약에 따른 변환만 추가. Sheet2의 fitted $S,T$, NF, excess-noise 회귀값은 calibration 관측값에서 제외.

## 새 monotone fit 전 필수 판별

아직 새 fitting 안 함. 아래는 이전 continuous two-piece affine 방법의 실행 전 수학적 점검 순서.

1. 정확한 sign/origin 변환과 apparatus 입력을 고정하고 **현행 보정 OFF** 원자 모델을23개 실제 detuning에서 직접 평가. 옛 Gold $C_{mix}$·옛 affine·옛 frozen baseline 숫자는 재사용하지 않음.
2. 팔별 model gain $x_i$와 동일 기준면 target $y_i>0$를 구성. 시드 입력이 가정이면 그 조건을 함께 기록. 광손실·QE·noise 계수로 gain 오차를 숨기지 않음.
3. $x_i\le x_j$인데 $y_i>y_j$이면 모든 비감소 $F$는

   $$\max_k\left|F(x_k)/y_k-1\right|\ge\frac{y_i-y_j}{y_i+y_j}$$

   를 만족. 모든 역전 pair의 최대를 팔별로, 두 팔의 최대를 공통 acceptance 하한으로 보고. 거의 같은 $x$는 모델 수치 오차를 고려해 판정. 동일 $x$의 다른 target도 같은 식 적용.
4. 비감소 $F(G)$는 고유 peak를 이동시키지 못함. 새 target의 peak 위치/팔별 ordering과 불일치가 크면 knot·계수만 늘리기 전에 이 방법의 표현 한계를 판정. 별도 empirical detuning registration을 도입하면 물리적 부호/원점 변환과 구분하고 독립 parameter·scope로 기록.
5. 방법이 표현 가능한 경우 고정 knot의 비음수 양쪽 slope·상대 minimax LP를 두 팔 모두 평가. Knot 선택도 같은 데이터에 맞추면 training임을 기록. In-sample 잔차는 confidence interval/holdout 아님.

출력 교정은 **측정 ratio readout**에만 적용. Canonical transfer·원자 response·dressed absorption·seed RIN/NF·GC no-fit 경로 보존. Corrected gain의 제곱근으로 새 quantum $M$을 만들지 않음. 옛 `apply_power_gain`의 source-plane pump cap/absorbing-probe identity guard도 measured-ratio 경계로 자동 전용하지 않음. 역사 workbook/fit/source capsule bytes 유지, 새 결과는 새 model ID와 날짜 경로 사용.
