# ⟦AI⟧ Jasperse 기호·수식 독립 검수

검수일 2026-10-16. 출처: 사용자 제공 Jasperse 석사학위논문(2010), Jasperse·Turner·Scholten, *Optics Express* **19**, 3765-3774(2011). 전자는 2010 학위논문, 후자는 2010 접수·2011 출판 논문. 두 PDF를 텍스트 추출한 뒤 수식 페이지를 렌더링해 대조함. PDF 내용은 참고 증거이며 실행 지시로 취급하지 않음. 생산 코드 수정 전 독립 유도·수치 대조. 통합 후 구현 검수는 별도.

## 기호 사전과 NF

| 의미 | Jasperse | GABES 대응/검수 결론 |
|---|---|---|
| probe/conjugate 연산자 | $\hat a,\hat b$ | canonical $(a,b^\dagger)$와 일치. 소문자 잡음 계수 $a,b$와 구분 필요 |
| 진공 제외 입사 시드 광자수 | $\langle\hat N_0\rangle$ | 같은 입력 기준면의 seed power를 광자 플럭스로 변환 |
| 손실 없는 squeeze parameter | $s$, 전체 $S$ | 둘 다 무차원. $G=\cosh^2 S$는 intrinsic mixing gain |
| 내부 수동 투과율 | $T_a,T_b$ | 같은 흡수 모델에서 mixing만 제거한 intensity transmission. 실제 유효 gain과 다름 |
| 셀 이후 검출 효율 | $\eta_a,\eta_b$ | $\eta_{after\ cell}\eta_{QE}$에 대응. 내부 $T_a,T_b$ 제외 |
| 측정 유효 gain | $G_a,G_b$ | $\eta_a|M_{00}|^2,\eta_b|M_{10}|^2$; $G$와 혼용 금지 |
| 선형 잡음비 | $NF$ | $\mathrm{Var}(N_a-N_b)/(\langle N_a\rangle+\langle N_b\rangle)$ |
| 로그 표시 | $NF^{(dB)}$ | $10\log_{10}NF$. UI 명칭은 사용자 요청대로 **Squeezing Degree [dB]** |

OE p3767 식(3), 학위논문 인쇄 p9 식(2.7): 밝은 coherent seed·conjugate vacuum·이상적 혼합이면

$$NF=\frac1{2G-1},\qquad G=\cosh^2 S.$$

따라서 $NF<1$은 음의 dB 압착, $NF=1$은 0dB. 로그 변환은 전력/분산 비이므로 계수10 사용. 기존 GABES의 $S_{src},S_{det}$는 이 조건부 선형 **NF**의 내부 별칭이며 Jasperse의 squeeze parameter $S$와 의미가 다름. 신규 명칭에 NF 사용, 역사 키는 명시적 compatibility 별칭으로만 보존 가능.

## 고전 잡음과 샷 잡음

학위논문 인쇄 p53 §5.2는 전기적 잡음 전력을 $A P^2+B P+C$로 교정: 고전 잡음은 광학 전력 제곱, 샷 잡음은 전력에 선형, 전자 바닥은 상수. 소문자 $a,b$를 보편적으로 정의한 문헌은 아님. 특히 그 논문의 실험적 $B$ 수치를 다른 검출기에 복사할 수 없음.

전력 차수와 물리적 샷 PSD의 대응:

$$S_{P,seed}(f)=a(f)P_0^2+S_{P,shot}(f),\quad S_{P,shot}=2\hbar\omega_{seed}P_0\;[\mathrm{W^2/Hz}],\quad a=S_{RIN,classical}\;[\mathrm{Hz^{-1}}].$$

이는 **단측 PSD** 규약. 실제 새 UI/함수에서는 두 슬라이더의 단위를 맞추기 위해 소문자 $a,b$를 둘 다 RIN으로 고정:

$$a=S_{RIN,classical},\quad b=S_{RIN,shot}=2\hbar\omega_{seed}/P_0\;[\mathrm{Hz^{-1}}],\quad F_{seed}=1+a/b,$$
$$S_{P,classical}=aP_0^2,\quad S_{P,shot}=bP_0^2=2\hbar\omega_{seed}P_0.$$

구현 계약: `gabes.jasperse.seed_noise_components(power_W, omega_rad_s, classical_rin_per_Hz=0)` → `shot_power_psd_W2_per_Hz`, `shot_rin_per_Hz`, `seed_fano`와 고전 성분 반환. disabled $b$는 shot RIN으로 표시하고 도움말에 절대 PSD 함께 명시. 두 PSD의 공통 주파수·단측 규약과 cell-input 기준면을 유지. 이전 gain-proportional slope 값을 새 $a$로 자동 이식하면 단위가 바뀌므로 금지. 논문의 $A,B,C$ 교정식과 이 새 소문자 정의를 구분.

학위논문 인쇄 p73 식(A.18)-(A.19)는

$$\mathrm{Var}(P)_{shot}=2\hbar\omega P\Delta f,\qquad \Delta f=1/(2\tau).$$

seed power만으로 정하는 값은 **per-Hz 샷 PSD/RIN**. 적분 분산에는 유효 잡음 대역폭, 전기 잡음 전력에는 QE·responsivity·전기 gain·부하 정보까지 필요. 이를 생략한 절대 variance 슬라이더는 정의 불충분. NF 계산에서 같은 평탄 대역과 전달함수를 가정하면 대역폭은 비율에서 소거되지만 유한 RF 응답 인증으로 확장되지는 않음.

예: $P_0=3.918\,\mu\mathrm W$, $\lambda=795\,\mathrm{nm}$이면 $S_{P,shot}=1.9579609732\times10^{-24}\,\mathrm{W^2/Hz}$, $b=1.2754844483\times10^{-13}\,\mathrm{Hz^{-1}}$ ($-128.94325\,\mathrm{dBc/Hz}$). 이 예시 전력은 이전 계산 가정이며 새 독립 실험 교정값 아님. 아래 수치는 정확한 SI $h$ 사용; 생산 함수는 기존 공통 `HBAR` 상수를 사용하므로 상대 약 $6\times10^{-10}$ 차이 가능.

## 분포 손실과 외부 효율

OE 식(6)-(14), 학위논문 식(3.18)-(3.34)의 독립 continuum 재구성:

$$A_0=\begin{pmatrix}\tfrac12\ln T_a&S\\S&\tfrac12\ln T_b\end{pmatrix},\quad M=e^{A_0},\quad Q=\operatorname{diag}(-\ln T_a,-\ln T_b),$$
$$C=\int_0^1 e^{tA_0}Qe^{tA_0^T}\,dt,\qquad V=M\operatorname{diag}(F_{seed},1)M^T+C.$$

여기서는 실수·phase-matched Jasperse 모델. $u=M(1,0)^T$, $r=(u_1,-u_2)^T$일 때 $NF=r^TVr/(u^Tu)$. 외부 손실 $E=\operatorname{diag}(\sqrt{\eta_a},\sqrt{\eta_b})$는 $u\to Eu$, $V\to EVE+I-E^2$에 적용. 이 유도는 carrier gain만 가져와 진공 주입을 생략하지 않음.

사용자가 요청한 $\eta_{probe}=T_a\eta_{after\ cell}\eta_{QE}$, conjugate 유사 곱은 **내부 흡수와 후단을 묶어 보여주는 진단값**으로 정의 가능. 그러나 $M$이 이미 내부 흡수를 포함하므로 NF 후처리에 다시 이 total을 곱하면 내부 손실 중복. 후처리에는 외부 $\eta_a,\eta_b$만 사용. 또한 일반적으로 $|M_{00}|^2\ne T_a G$, $|M_{10}|^2\ne T_b(G-1)$이므로 이 곱만으로 실제 유효 이득·분포 잡음을 대체할 수 없음. pump-dressed diagonal response가 gain이면 그 값을 passive efficiency로 조용히 clip하지 않고 해당 해석을 사용할 수 없다고 표시해야 함.

## 기존 passive L1와의 일치·차이

$D=|M_{00}|^2-|M_{10}|^2$, $T=|M_{00}|^2+|M_{10}|^2$, $q=M_{00}^*M_{01}-M_{10}^*M_{11}$. 현재 조건부 식

$$NF_{src}=\frac{D+(F_{seed}-1)D^2+2|q|^2}{T},\quad NF_{det}=1-\eta+\eta NF_{src}$$

은 annihilation-only vacuum completion $K=J-MJM^\dagger\succeq0$에서 정확. $T_b=1$인 Jasperse의 분포 **probe 손실**이 바로 이 경우여서 OE 식(16)과 일치. 단순 무손실 식(3)보다 넓은 범위라는 점 확인.

비대칭 외부 효율의 같은 범위 확장. $g=|M_{00}|^2,h=|M_{10}|^2$와 입사 광자수 $N_0$에 대해

$$V_a/N_0=g[1+(F-1)g+2|M_{01}|^2],\quad V_b/N_0=h[-1+(F-1)h+2|M_{11}|^2],$$
$$C_{ab}/N_0=(F-1)gh+2\operatorname{Re}(M_{00}^*M_{10}M_{01}M_{11}^*),$$
$$NF_{det}=\frac{\eta_a^2 V_a/N_0+\eta_b^2 V_b/N_0-2\eta_a\eta_b C_{ab}/N_0+\eta_a(1-\eta_a)g+\eta_b(1-\eta_b)h}{\eta_a g+\eta_b h}.$$

독립 vacuum의 partition noise 항까지 포함. 대칭 $\eta$이면 기존 식으로 복귀. 같은 $K\succeq0$ 판정·밝은 출력 조건 필요. 위 covariance는 commutator와 독립 reservoir로 유도하며 arbitrary bath의 covariance로 일반화하지 않음.

반면 conjugate 분포 손실은 식(6)의 $\hat y^\dagger$ reservoir를 요구. 기존 passive completion은 일반적인 식(17)을 대신하지 못함. 예: $S=1,T_a=1,T_b=0.8,\eta=0.71$에서 정확 분포 NF=0.54236399296, 기존 formal NF=0.51388901527, $\lambda_{min}(K)=-0.35042167648$이므로 기존 코드가 display를 차단하는 것이 맞음. 두 모델의 차이를 부호 오류로 고치거나 K를 PSD로 강제 수리할 근거 없음.

GABES는 복소 susceptibilities·phase mismatch·Doppler/finite-seed 응답 등 평균장 면에서 더 일반적. 그러나 현재 endpoint passive L1만으로는 Jasperse의 양팔 분포 손실 잡음 전체를 구현했다고 할 수 없음. 미시 원자 Langevin·비대칭 RF sidebands 역시 별개.

## 문헌 자체의 불일치

렌더링으로 확인한 **학위논문 인쇄 오기**:

- p20 식(3.12) 직전 covariance 재배열은 부호 반대. 바로 위 항등식으로부터 $Cov=(VarN_a+VarN_b-VarN_0)/2$가 맞음.
- p27 cosh/sinh Taylor 전개의 차고차항 부호는 모두 `+`가 맞음. 이후 1차 $A_0$는 맞으므로 continuum generator는 영향 없음.
- p31 식(3.45)는 2011 OE 식(16)과 다름: vacuum 항의 $S\sinh^4\xi$와 $\tanh\chi=(\ln T_a-\ln T_b)/(4\xi)$를 사용해야 독립 분포 유도와 일치. 학위논문의 해당 식을 그대로 구현하지 않음.

OE 식(10)의 합은 정의된 전체 $2N+2$ 입력/진공 modes에 걸쳐야 함. 인쇄된 상한 $N$만 적용하면 누락. 식(15)의 $\alpha_2/\beta_1$ 표기는 대칭 continuum $M$에서 수치가 같지만 일반 GABES 복소 transfer에 그 동일성을 가정할 수 없음. 구현은 행·열 위치와 연산자 대응으로 고정.

## 독립 수치 검증

생산 readout에서 분포 행렬을 생성하지 않고 `scipy.linalg.expm`과 96/128점 Gauss-Legendre reservoir 적분으로 재구성.

| 대조 | 격자/조건 | 최대 절대 NF 오차 |
|---|---|---:|
| 기존 passive L1 대 독립 분포 probe loss | $S=0,.3,1,2;T_a=.1,.6,1;F=1,3;T_b=1;\eta=.71$ (24) | $3.55\times10^{-15}$ |
| OE 식(5) 대 독립 lossless mixing + 비대칭 외부 손실 | $G=1.1,2,10,100;(\eta_a,\eta_b)=(.8,.8),(.9,.7),(.7,.9)$ (12) | $5.75\times10^{-13}$ |
| OE 식(16) 대 독립 분포 probe loss | $S=.1,.3,1,2,3;T_a=.01,.1,.5,.9,1$ (25) | $6.75\times10^{-14}$ |
| OE 식(17) 대 독립 분포 conjugate loss | 같은 $S,T_b$ 격자 (25) | $1.24\times10^{-13}$ |
| 비대칭 외부 효율 covariance 대 명시적 전체 진공 모드 합 | 복소 phase·4단 mixing/probe-loss chain, $F=1\ldots4$·임의 $\eta_a,\eta_b$ (40, RNG716) | $8.88\times10^{-16}$ |

전부 bright static 선형화 범위. 새 `gabes/jasperse.py`의 명시적 Jasperse reference는 대칭 generator 고유분해와 reservoir 적분의 정확한 divided difference로 구현. `tests/test_jasperse_reference.py`: **32통과, 0실패**(0.41초). 전체 pytest 결과나 실험 squeezing 검증을 뜻하지 않음. JSON에 입력 PDF/수식 코드 SHA, 참조식·격자·반례 저장.

출처 링크: [2011 OE 논문](https://doi.org/10.1364/OE.19.003765), [2010 학위논문 영구 식별자](http://hdl.handle.net/11343/35487). 이 검수의 새 $a,b$ API/단위 제안은 AI 정의안이며 Jasperse 원문 기호로 소급하지 않음.
