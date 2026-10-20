# ⟦AI⟧ dressed dipole의 양의 흡수·방출률: 독립 QRT 유도

2026-10-16. 사용자 요구: conjugate에도 Jasperse 내부 투과에 대응하는 물리량 계산. 아래는 같은 GABES 4준위 GKSL 모델에서 흡수/방출을 따로 정의하는 유도·작은 반례 검증. $M$의 net gain에서 양의 잡음을 역추정하지 않음. [재현 코드](qrt_absorption_derivation.py), [수치](qrt_absorption_derivation.json). 다른 에이전트의 결과를 받기 전에 정적 유도·파일럿 완료, 이후 독립 경로를 교차 대조함.

## 1. 정적 ordered QRT와 단위

원자 lowering dipole $D_j=\sum_e c_{g_je}|g_j\rangle\langle e|$. Pump frame의 각 주파수 $\omega_j$는 optical beat이며 검출 RF와 다름. $R_\omega=-(\mathcal L+i\omega)^{-1}$를 trace-zero 부분공간에서 사용. $\rho$는 **같은 generator의 정상 상태**.

$$S_{abs,j}(\omega)=2\operatorname{Re}\operatorname{Tr}\{D_jR_\omega[D_j^\dagger\rho-\langle D_j^\dagger\rangle\rho]\},$$
$$S_{em,j}(-\omega)=2\operatorname{Re}\operatorname{Tr}\{D_jR_\omega[\rho D_j^\dagger-\langle D_j^\dagger\rangle\rho]\}.$$

단위 seconds. 앞 식은 $\int e^{i\omega t}\langle\delta D(t)\delta D^\dagger(0)\rangle dt$, 뒤는 동일 주파수의 반대 ordering. 물리적 stationary GKSL과 연결 상관함수의 Wiener-Khinchin 성질로 각각 비음수. 원자의 excited population만 이용한 bare Lorentzian 분해는 dressed coherence·Raman 경로를 놓침.

통상적인 Rabi 입력 $H_{test}/\hbar=(\Omega D^\dagger e^{-i\omega t}+h.c.)/2$의 coherence susceptibility $\chi_{raw}=\delta\langle D\rangle/\Omega$에 대해

$$\chi_{raw}=-\frac i2\operatorname{Tr}D R_\omega[D^\dagger,\rho],\qquad S_{abs}-S_{em}=-4\operatorname{Im}\chi_{raw}.$$

계수4는 full two-sided ordered spectrum의 계수2와 Rabi drive의 $1/2$에서 발생. GABES physical susceptibility는 $\chi_{phys}=-2nd_{eff}^2\chi_{raw}/(\epsilon_0\hbar)$. 따라서

$$\alpha_{abs,j}=\frac{n\omega_{opt,j}d_{eff}^2}{2\epsilon_0\hbar c}\,S_{abs,j},\qquad\alpha_{em,j}=\frac{n\omega_{opt,j}d_{eff}^2}{2\epsilon_0\hbar c}\,S_{em,j},$$
$$\alpha_{net,j}=\alpha_{abs,j}-\alpha_{em,j}=k_j\operatorname{Im}\chi_{phys,j}=-2\operatorname{Re}A_{jj}.$$

동일 결과: photon-flux mode $E=Qa$, $Q^2=2\hbar\omega_{opt}/(\epsilon_0cA)$, $g=d_{eff}Q/(2\hbar)$, 선밀도 $\lambda=nA$이면 $\alpha=\lambda g^2S$. 임의 면적 $A$는 소거. Atomic beat $\omega$와 양의 optical carrier $\omega_{opt}$를 혼동하면 안 됨.

**Nambu ordering 주의:** `reduced_readout_operators()`는 $(D_p,D_c^\dagger)$이므로 creator conjugate 좌표의 greater/lesser 순서는 physical lowering 좌표와 반대. `LocalNambuGenerator`에서 probe는 $\alpha_{abs}=D^>_{pp},\alpha_{em}=D^<_{pp}$; conjugate는 $\alpha_{abs}=D^<_{cc},\alpha_{em}=D^>_{cc}$. 두 팔 모두 $\alpha_{net}=-2\operatorname{Re}A_{jj}$ 확인 가능.

## 2. Jasperse $T_b$와의 정확한 관계

독립 흡수 채널의 투과는

$$T_{abs,j}=\exp[-\int_0^L\alpha_{abs,j}(z)\,dz],\qquad0<T_{abs,j}\le1.$$

공간 균일 prescribed state라면 $e^{-\alpha_{abs,j}L}$. 이는 물리적 absorption-channel survival을 정의하므로 conjugate net gain이 있어도 계산 가능. 반면 $e^{-\alpha_{net}L}$은 net diagonal transfer로 1을 넘을 수 있음. 두 값은 다른 물리량.

Jasperse의 $T_a,T_b$는 **수동 손실 reservoir만 있는** 분포 모델의 transmission. $\alpha_{em}=0$이면 위 $T_{abs}$가 바로 그 $T$와 일치. $\alpha_{em}>0$인 GABES에서는 양의 흡수 채널을 표시할 수 있지만, 같은 $T$ 하나만 Jasperse 식(16)/(17)에 넣어 전체 잡음이 같다고 주장할 수 없음. 방출률과 양팔 교차 상관·진공 주입도 유지해야 함. Total efficiency를 $T_{abs}\eta_{after}\eta_{QE}$로 표시해도 현 전체 $M$ 뒤에 그 내부 손실을 또 적용하면 중복.

## 3. 실제 production dissipator를 양의 jump로 재표현

추가 ground coherence decay $\gamma_{gg}$, optical coherence decay $\gamma_{opt}$의 정확한 동일 표현:

$$J_1=\sqrt{\gamma_{gg}/2}\,\operatorname{diag}(1,-1,0,0),$$
$$J_2=\sqrt{2(\gamma_{opt}-\gamma_{gg}/4)}\,\operatorname{diag}(0,0,1,1).$$

조건 $\gamma_{opt}\ge\gamma_{gg}/4$. 기존 radiative jumps·thermal replacement jumps와 함께 사용. 예시 $T=386.15K$에서 $\gamma_{gg}=11218.9920048\,s^{-1}$, $\gamma_{opt}=2916782.47158\,s^{-1}$, dissipator 상대 차이 $1.44\times10^{-16}$. 기존 `ExplicitReservoirs.from_atom`의 nonzero legacy-dephasing 거부를 우회한 근거 없는 import가 아니라, 각 coherence rate를 재현하는 explicit jump 증명. 조건이 깨지는 다른 입력에는 이 분해를 강제 적용할 수 없음.

## 4. finite seed 상태: periodic QRT가 필요

Production `chi_matrix_table`은 probe-only와 conjugate-only의 **서로 다른 finite drive 상태 두 개**에서 coherence/drive 비를 추출. 이는 한 공통 운전점의 infinitesimal Jacobian이 아님. 따라서 pump-only QRT를 그대로 finite-seed 물리량이라고 부르거나, QRT net rate와 finite quotient가 반드시 같아야 한다고 요구할 수 없음.

현재 입력 probe와 pump가 주어진 입구 조건의 정확한 공통 상태를 사용하려면

$$H(t)/\hbar=H_0+V e^{-i\nu t}+V^\dagger e^{i\nu t},\qquad V=\Omega_s D_p^\dagger/2,\quad\nu=-\Omega_{HF}+\delta,$$

에서 전체 $\rho_n$을 유지. Conjugate 입력은 vacuum. Probe lowering spectrum은 physical beat $\nu$, conjugate는 $-\nu$에서 평가. Floquet quasifrequency 0으로 접으면 각각 harmonic $+1,-1$. 동위원소/velocity/drive/dissipator·전이 계수와 단위를 같게 유지해야 함. 전파 중에는 생성된 conjugate carrier와 pump/seed 변화에 따라 local state를 재계산해야 하므로 입구의 prescribed-state 값을 전파가 해결된 값으로 승격하지 않음.

기존 `PeriodicAtomicNoise`가 이 상태의 drift·jump diffusion과 양 ordering을 이미 제공함. **$\rho_0$만 평균하여 정적 QRT에 넣으면 안 됨.** 본 파일럿은 전체 periodic state와 covariance를 사용.

## 5. 빠른 독립 연결 QRT 알고리즘

완전 traceless Hermitian basis $F_i$에 대해 $\mu_i(t)=\operatorname{Tr}[F_i\rho(t)]$,

$$C_{ij}(t)=\operatorname{Tr}[F_iF_j\rho(t)]-\mu_i(t)\mu_j(t),\quad w_i=\operatorname{Tr}(DF_i).$$

$\rho_n$의 차수가 $N$이면 mean-product convolution은 최대 $2N$까지 계산. Absorption source는 $C_nw^\dagger$, emission source는 $C_n^Tw^\dagger$. 같은 periodic traceless drift $A_n$에 대해

$$[A_0+i(\omega+n\nu)I]y_n+A_+y_{n-1}+A_-y_{n+1}=-s_n,$$
$$S(\omega)=2\operatorname{Re}(wy_0).$$

15×15 block Thomas로 두 ordering RHS를 함께 풂. 전체 $RDR^\dagger$ 행렬 생성·큰 행렬의 모든 RHS inversion을 피함. 이는 QRT의 다른 계산식이며 drift에서 noise를 정하는 최소 PSD completion이 아님. 정확한 GKSL·물리적 periodic state에 대한 full spectrum positivity가 근거. 유한 mean/response truncation은 별도 수렴 필요; 음의 수치를 clip하면 안 됨. Explicit thermal reset이 양수이면 traceless evolution에 해당 decay가 있어 연속 spectrum의 안정성 근거가 됨.

## 6. 독립 대조 결과

조건: blue OPD+1GHz, internal $\delta=+23.7MHz$, $T=386.15K$, pump294.8mW/530µm, probe3.918µW/330µm, 길이12.5mm, 단일 velocity0. Weak-field 환산에는 **기존** $d_{eff}^2=d^2(0.74/12)$ 사용. 이는 기존 normalization의 활용이며 새 no-fit absolute reduction 증명 아님.

| 상태 | probe $\alpha_{abs}/\alpha_{em}$ [m⁻¹] | conjugate $\alpha_{abs}/\alpha_{em}$ [m⁻¹] |
|---|---:|---:|
| Pump-only | 약22.36/6.93 | 약0.717/3.006 |
| 실제 finite seed | 53.964863/38.564247 | 5.025683/7.309750 |

Finite conjugate에서 $\alpha_{net}=-2.284067m^{-1}$이지만 $T_{abs,c}=0.93911152$로 물리적 흡수 채널 존재. 위 단일 velocity 값을 Maxwell 평균으로 오인하면 안 됨.

- 직접 연결 QRT 대 기존 full periodic $RDR^\dagger$: seed0, $10^{-6}$ power, 실제 power의 3조건에서 상대 약 $10^{-12}$ 일치.
- 두 수준 원자의 ground/incoherently excited 상태 Lorentzian 9조건: 최대 상대 오차 $2.20\times10^{-15}$. $S_{abs}=2p_g\Gamma_2/[\Gamma_2^2+(\omega-\omega_0)^2]$, $S_{em}$은 $p_g\to p_e$.
- 실제 finite seed의 pump-frame/production seed-frame 상태를 gauge 변환하여 64phase 대조: 최대 행렬 norm 차이 $1.30\times10^{-14}$.
- velocity $-200,0,+200m/s$에서 mean/response NF3→4 최대 상대 차이 $5.95\times10^{-15}$. 해당 셋은 Doppler 수렴의 증명이 아님.
- 실제 velocity0의 response NF2→3 약 $10^{-9}$, 3→4는 약 $4\times10^{-16}$ 이하. 계산 source action 약1.3-1.5ms/두 팔·velocity (이미 계산된 $\rho_n$ 기준). 전체 UI latency 주장은 아님.
- 같은 finite 상태의 infinitesimal conjugate $\operatorname{Im}\chi=2.524855708\times10^{-13}s$; production의 별도 conjugate-only quotient는 $2.530053938\times10^{-13}s$. 약0.21% 차이는 서로 다른 상태/미분 정의. 강제로 net rate를 맞춰 보정하지 않음.

판정: **양의 conjugate absorption-channel transmission 계산의 물리적 정의·동일 finite 입구 상태·독립 계산 경로 확보**. Maxwell 평균/전파·전체 미시 NF 통합은 별도 검수 범위. 공식 Grand Challenge 절대 gain/squeezing/thermal acceptance를 닫는 결과는 아님.

후속 구현 적대적 검토: `gabes/dressed_cell_absorption.py`의 연결 QRT가 처음에는 net 값을 같은 두 ordered source의 차이로 정의해 Kubo residual 검사가 항등식이 되는 점 지적. 구현 담당자가 $\rho_n$에서 직접 만든 **제3 commutator RHS**로 수정·재검증했다고 보고함. 이것이 독립 비교이며 covariance source 두 개와 같은 식을 반복하는 검사는 아님. 담당자의 별도 Maxwell 교차 검토: full engine 210 velocity 평가와 약 $2\times10^{-13}$ 일치, adaptive 상대 tolerance0.005→0.001의 rate 변화 $10^{-10}$ 미만. 이는 해당 담당자의 독립 실행 보고이며 본 파일럿의 3 velocity 결과와 구분. Quadrature error는 수치 적분기가 보고한 추정치로 표기하고 엄밀한 전구간 bound로 부르지 않음.
