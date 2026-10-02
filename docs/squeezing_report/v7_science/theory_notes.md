# 독립 가우시안 증폭기 이론: v7 재현 자료

이 스크립트는 원자 GABES 엔진이 아니다. McCormick et al., PRA 78, 043816 (2008),
https://arxiv.org/abs/quant-ph/0703111 의 분포 이득/흡수 모형을 연속 극한으로 계산한다.
강한 coherent seed, conjugate vacuum, 펌프 비고갈, 완전 위상정합, 한 공간모드,
주파수에 무관한 계수, probe만의 분포 흡수 및 진공 유입을 가정한다.

## 정의와 계산

차원 없는 위치 z/L에서 r=arcosh(sqrt(G)), a=-ln(T),
Mx=[[-a/2,r],[r,0]], Mp=[[-a/2,-r],[-r,0]], A=diag(Mx,Mp).
G는 흡수를 끈 intrinsic gain, T는 이득을 끈 probe 투과도이다.
Vvac=I/2, D=diag(a/2,0,a/2,0), E=exp(A),
Vout=E Vvac E^T + integral_0^1 exp(A s) D exp(A^T s) ds.
입력 seed 진폭을 1로 두면 출력 평균 진폭은 E[:2,0]이다.
외부 효율 eta_p,eta_c는 Vdet=L Vout L^T+(I-L L^T)/2에 적용한다.
L=diag(sqrt(eta_p),sqrt(eta_c),sqrt(eta_p),sqrt(eta_c)).
q=(sqrt(eta_p)*Ap, -w*sqrt(eta_c)*Ac, 0, 0)일 때
R=2 q^T Vdet q / (eta_p Ap^2 + w^2 eta_c Ac^2), S=10log10(R).
광손실과 전자 가중치는 다른 연산이며 SNL에도 w^2가 필요하다.

## Gold 평균 출력만으로 제약한 예측

8 uW seed, 111/109 uW 출력은 Gp=13.875, Gc=13.625다.
평균 출력 두 개만 적합한 결과 G=16.025936244, T=0.821164524.
스퀴징 데이터는 적합하지 않았다. eta=.8694에서는 -7.818903 dB,
논문 총손실 13.5%를 그대로 쓴 eta=.865에서는
-7.709262 dB를 예측한다. 측정값은 -7.8 dB다.
이는 평균광학량으로 제약된 단순 양자모형의 일관성 검사이며,
원자 파라미터에서 출발한 GABES의 독립적인 전방 예측 검증이 아니다.
반올림된 광출력, 손실 기준면 및 효율 정의의 차이가 첫 예측의 작은 잔차보다 크다.
흡수의 공간 분포/공간 모드/기술잡음이 달라지면 같은 평균 출력에도 잡음은 달라진다.

## 실제 유한 최적점과 허용범위

적합 T를 고정하면 G=20.252754에서 -7.845024 dB (검색 G=1..200).
적합 G를 고정하면 T=0.930920에서 -8.096600 dB (검색 T=.3..1).
분포 흡수로 유입된 진공도 뒤의 이득에서 증폭되어 큰 G의 불이익이 생긴다.
공통 광학깊이 스케일 x가 r과 a를 함께 곱할 때 x=1.003382가 최적이다.
이 x는 추가 가정 없이는 온도/P/셀길이의 물리적 허용범위가 아니다.
표시 2D 영역 G=1..50,T=.65..1의 최솟값은 G 경계에 있으므로 전역 최적점이라 부르지 않는다.
허용범위는 theory_results.json의 distributed_plus0p5db_tolerances에 있다.

## 이상적 검출 허용범위

기준 G=15, eta=.8694, w=1에서 R0=0.160579310, S0=-7.943104 dB.
K=2G-1, 입력 seed 초과 Fano factor xi, 출력 SNL 기준 가산잡음 epsilon에 대해
R=1-eta+eta(1+xi)/K+epsilon.
독립 coherent 배경의 전체 SNL 중 비율 f가 있으면 Rtotal=(1-f)R+f.
동일한 잡음 여유를 공유하므로 표의 단일변수 허용범위를 동시에 최대까지 쓸 수 없다.
eta 개선, gain 개선, seed 잡음 억제, 배경 억제의 예산과 전자 가중치 허용구간은 JSON에 있다.
가중치 최적값 sqrt(G/(G-1))는 DC 균형 G/(G-1)과 다르다.
w 최적화 결과를 w=1의 통상 IDS와 동일한 측정값으로 비교하지 않는다.

## 유한 seed의 정확한 공분산과 다변수 함수

이상적 비축퇴 증폭기에서 입력 seed 평균 광자수 n, 분산 n+E를 두고 H=G-1이면,
mu_p=G n+H, mu_c=H(n+1),
V_p=G(2G-1)n+GH+G^2 E,
V_c=H(2G-1)n+GH+H^2 E,
C_pc=GH(2n+1)+GH E.
손실 이후 mu'_j=eta_j mu_j,
V'_j=eta_j^2 V_j+eta_j(1-eta_j)mu_j, C'=eta_p eta_c C.
R_w=(V'_p+w^2 V'_c-2w C')/(mu'_p+w^2 mu'_c).
독립 coherent 배경 Bp,Bc는 분자와 분모에 Bp+w^2 Bc를 함께 더하며,
독립 전자잡음은 분자에만 더한다. 기술잡음이 있는 배경은 자체 분산을 사용해야 한다.
이 공분산 관계는 본 스크립트의 bright-seed 근사보다 일반적이다.

K=2G-1, D0=eta_p G+w^2 eta_c(G-1)인 bright limit에서는
R_w=1+[2G(G-1)(eta_p-w eta_c)^2-2w^2 eta_c^2(G-1)]/D0
    +xi[eta_p G-w eta_c(G-1)]^2/D0.
따라서 양팔 광손실, 입력 seed 잡음, 전자 가중치는 서로 결합되어 있다.
G=15, eta_p=eta_c=.8694 기준 +0.5 dB 범위에서 다른 팔을 고정한 추가
probe 광손실 허용량은 5.9730%, conjugate 광손실은 1.3582%이다.
이는 양팔 공통 효율 허용치와 다른 조건이다.

## 실험 제어변수의 결합과 허용영역

Gaussian pump의 Omega는 sqrt(P)/w에 비례하므로,
dln(Omega)=0.5 dln(P)-dln(w). 같은 Rabi 주파수를 유지하려면
dP/P=2 dw/w가 필요하다. 그래도 transit와 mode overlap은 유지되지 않는다.
코드의 Rb 증기밀도식은 121 C에서 dln(N)/dT=.0583603 /K를 준다.
고정 NL 방향은 dL/L=-.0583603 dT이다. 1 K 상승은 밀도를 약 5.993% 증가시킨다.
이것은 optical-depth 보상 관계이며 스퀴징 허용 온도폭을 직접 뜻하지 않는다.
full Maxwell Raman-Doppler rms는 sigma_v |k_p-k_s|/(2pi),
vector magnitude를 사용하며 작은 각도에서는 sigma_v theta/lambda이다.
121 C, .32 deg에서는 약 1.380 MHz이고 .01 deg 증가는 약 43.13 kHz다.
분포 rms는 관측 gain FWHM이나 squeezing bandwidth와 동일하지 않다.

임의 실험변수 q에 대해 R=1-eta+eta/(2G-1)+epsilon의 기울기는
dR=-(1-1/K)deta-2eta dG/K^2+depsilon이다.
유한 최적점에서는 이득 증가의 이점이 손실/추가잡음 증가와 상쇄된다.
실제 stationary point의 dB Hessian H에서 +b dB 허용영역은
delta_q^T H delta_q <=2b이다. 다른 변수 고정 시 폭은 sqrt(2b/Hii),
다른 변수를 재최적화할 때 유효 곡률은
Hii-Hi,-i inv(H-i,-i) H-i,i (Schur complement)이다.
기울기가 0이 아닌 점/검색 경계에서는 이 대칭 타원 공식을 쓰지 않는다.
통계적 공분산 Cq가 있으면 일차 분산은 grad(S)^T Cq grad(S)이며,
stationary point의 평균 성능 저하는 0.5 trace(H Cq)부터 시작한다.

잔여 채널 지연 tau가 있고 이상적 공분산이 주파수 대역에서 일정하다면,
R(f,tau)=R0+[4eta G(G-1)/(2G-1)] [1-cos(2pi f tau)].
G15, eta.8694의 +0.5 dB 지연 허용값은 500 kHz에서 12.557 ns,
3.5 MHz에서 1.794 ns이다. 잡음 스펙트럼이 바뀌면 해당 주파수의 C(f)를 써야 한다.
검출 대역/두 광파의 군지연 보상은 함께 지정해야 한다.

## 검증

무손실 R=1-eta+eta/(2G-1), 순수 수동손실 coherent R=1,
전체 4x4 공분산의 최소 symplectic eigenvalue >=1/2를 검증했다.
독립 gain/loss 교대 분할 32/128/512단계가 정확한 연속 공분산으로 수렴한다.
상세 오차는 JSON validation에 기록한다.
