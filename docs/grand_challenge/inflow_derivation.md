# Thermal box inflow — 경계 flux, 실제 chord와 공통 입사 위상

2026-09-14. [구현](../../gabes/quantum/inflow.py), [검사](../../tests/quantum/test_inflow.py), [audit](../../analysis/grand_challenge/inflow_audit.py), [수치 보고서](inflow_report_v1.json), [원자 characteristic 유도](transport_derivation.md).

유한 직육면체를 통과하는 isotropic Maxwell gas의 경계 유입률을 밀도·온도·질량·기하에서 정하고, 각 원자의 실제 first-exit chord를 기존 `BallisticPath`에 전달한다. 이어서 경로와 **하나의 공통 RF 입사 위상**별로 제공된 finite atomic wavepacket을 marked-Poisson 합으로 연결한다. 계산에 맞추는 유입률이나 transit linewidth는 없다.

이는 충돌 없는 독립 원자의 prescribed boundary ensemble이다. 실제 벽의 산란·흡착·재유입 상관, 내부 충돌, beam aperture, self-consistent Maxwell propagation 또는 검출기 모델은 구현하지 않는다. 입력값에 실험 출처를 제공할 수 있지만, 아래 synthetic controls와 수학적 항등식만으로 실측의 no-fit squeezing 예측이 검증되지는 않는다.

## 경계의 Maxwell flux

Box bounds를 \(\mathbf l,\mathbf h\), 길이를 \(L_i=h_i-l_i>0\), 부피를 \(V=L_xL_yL_z\)라 하자. 명시적인 number density \(n>0\), temperature \(T>0\), mass \(m>0\)에서

\[
\sigma^2=\frac{k_BT}{m},\qquad
f(\mathbf v)=\frac{\exp[-|\mathbf v|^2/(2\sigma^2)]}
                  {(2\pi\sigma^2)^{3/2}}
\]

이며 \(f\)는 셀 내부의 부피 속도 분포다. Outward normal \(\hat{\mathbf n}\)을 갖는 면에 들어오는 원자의 rate measure는

\[
dJ=n f(\mathbf v)|\mathbf v\cdot\hat{\mathbf n}|\,
   \mathbf1_{\mathbf v\cdot\hat{\mathbf n}<0}\,dA\,d^3v.
\]

면에 수직인 미소 길이 \(|\mathbf v\cdot\hat{\mathbf n}|dt\) 안의 원자가 \(dt\) 동안 면을 통과하므로 이 속도 인자가 필요하다. 입사 사건에서 관측하는 속도 분포는 부피 속도 분포와 다르다.

Inward normal speed를 \(u=-\mathbf v\cdot\hat{\mathbf n}>0\)라 하면

\[
\int_0^\infty
u\frac{e^{-u^2/(2\sigma^2)}}{\sqrt{2\pi}\sigma}\,du
=\frac{\sigma}{\sqrt{2\pi}},\qquad
J_f=\frac{nA_f\sigma}{\sqrt{2\pi}}.
\]

따라서 면별 rate를 제외하고 정규화한 입사 mark 분포는

\[
p_f(\mathbf r_{\rm in},u,v_1,v_2)=
\frac1{A_f}
\underbrace{\frac{u}{\sigma^2}e^{-u^2/(2\sigma^2)}}_{\text{Rayleigh normal speed}}
\prod_{j=1}^{2}
\underbrace{\frac{e^{-v_j^2/(2\sigma^2)}}{\sqrt{2\pi}\sigma}}_{\text{Gaussian tangential velocity}}.
\]

면 위치 두 좌표는 uniform, inward normal speed는 Rayleigh, tangential 두 속도는 독립 Gaussian이다. 여섯 면을 모두 더하면

\[
J_{\rm total}
=\frac{2n\sigma}{\sqrt{2\pi}}(L_yL_z+L_xL_z+L_xL_y)
=\frac{nA_{\rm total}\langle|\mathbf v|\rangle}{4},
\quad
\langle|\mathbf v|\rangle=\sigma\sqrt{\frac8\pi}.
\]

Isotropic gas의 기본 식 \(\Phi=n\langle v\rangle/4\)는 [University of Texas의 molecular-flux 유도](https://farside.ph.utexas.edu/teaching/355/Surveyhtml/node213.html)와 일치한다. 여기서 box-specific conditional distribution과 수치 검사식은 위 Maxwell 적분으로 직접 유도했다. 이 참고자료가 GABES의 구현이나 quantum-optical 결과를 검증하는 것은 아니다.

## 직접 6면 quadrature와 first-exit chord

`maxwell_box_inflow(lower_corner_m, upper_corner_m, *, temperature_K, mass_kg, density_m3, points_per_face_power, seed, source)`는 각 면에 \(N=2^p\)개의 scrambled Sobol 5D 노드를 만든다. 다섯 uniform 좌표를 다음과 같이 변환한다.

\[
r_1=l_1+L_1U_1,\quad r_2=l_2+L_2U_2,\quad
u=\sigma\sqrt{-2\log(1-U_3)},\quad
v_1=\sigma\Phi^{-1}(U_4),\quad v_2=\sigma\Phi^{-1}(U_5).
\]

입사면 normal 좌표는 해당 boundary에 정확히 고정하고, normal 속도의 부호는 inward로 정한다. 모든 노드의 rate는 양수이며 같은 면 안에서 \(w_j=J_f/N\)이다. `rate_s_inverse`에는 이미 density, area, flux weighting이 들어 있다. 후속 합에서 \(|v_n|\)이나 면적을 다시 곱하지 않는다. Face IDs는 `x-, x+, y-, y+, z-, z+` 순서다.

각 축의 앞쪽 boundary까지의 시간은

\[
t_i=
\begin{cases}
(h_i-r_{{\rm in},i})/v_i,&v_i>0,\\
(l_i-r_{{\rm in},i})/v_i,&v_i<0,\\
\infty,&v_i=0.
\end{cases}
\qquad
\tau=\min_i t_i.
\]

이 \(\tau\)가 실제 첫 출구의 residence time이다. Normal 방향 두 면 사이의 거리만 사용하면 tangential 방향으로 먼저 빠져나가는 원자를 과도하게 오래 남긴다. `ThermalInflowQuadrature`는 entry face, inward velocity, 모든 면의 노드 수, 양의 rate, 실제 chord residence를 검사하며, `path(index)`가 해당 `BallisticPath`를 반환한다.

Sobol은 30-bit dyadic grid를 사용하고 \(2^{-31}\)을 더해 각 좌표를 open cube의 midpoint로 옮긴다. 표현 가능한 범위는

\[
2^{-31}\le U_i\le1-2^{-31}
\]

이므로 Gaussian inverse CDF와 Rayleigh 변환의 무한 endpoint가 발생하지 않는다. 이는 유한 수치 quadrature의 처리이며 물리적 velocity cutoff를 새 입력으로 둔 것이 아니다. 각 uniform coordinate의 표현 구간 밖 tail probability는 한쪽당 \(2^{-31}\)이지만, 이 작은 확률 자체가 residence나 고차 모멘트의 적분 오차 상한은 아니다. 무한 Maxwell tail, 매우 느린 원자와 짧은 chord를 포함하는 연속 적분에 대한 정확한 전체 오차 보장은 제공하지 않는다.

`seed`와 `points_per_face_power`는 물리 파라미터가 아니라 수치 제어다. 여섯 면은 한 seed에서 재현 가능한 별도 child sequences를 사용한다. 해상도와 독립 seed를 함께 바꾸어 수렴을 확인해야 하며, seed별 오차의 단조 감소나 독립 노드에 근거한 Monte Carlo error bar를 가정하지 않는다.

## 점유수와 phase-space 모멘트의 독립 검사

입사점에서 age \(a\in[0,\tau]\)를 따라 내부점으로 옮기는 좌표변환의 Jacobian은

\[
d^3r=|\mathbf v\cdot\hat{\mathbf n}|\,dA\,da
\]

이다. Convex box에서 거의 모든 내부점과 nonzero velocity는 하나의 과거 입사점과 하나의 age를 갖는다. 따라서 임의의 적분 가능한 \(g(\mathbf r,\mathbf v)\)에 대해

\[
\sum_f\int_{\partial V_f,\ v_n<0}\!
n f(\mathbf v)|v_n|\,dA\,d^3v
\int_0^{\tau(\mathbf r_{\rm in},\mathbf v)}
g(\mathbf r_{\rm in}+\mathbf va,\mathbf v)\,da
=n\int_Vd^3r\int_{\mathbb R^3}f(\mathbf v)g(\mathbf r,\mathbf v)\,d^3v.
\]

이 phase-space identity에서 \(g=1\)을 택하면

\[
\sum_j w_j\tau_j\longrightarrow nV,
\qquad
\overline{\tau}_{\rm arrivals}
=\frac{nV}{J_{\rm total}}
=\frac{4V}{A_{\rm total}\langle|\mathbf v|\rangle}.
\]

구현의 `mean_occupancy`는 실제 \(\sum_jw_j\tau_j\)다. `equilibrium_atom_number`는 독립 analytic 값 \(nV\)다. 앞의 합을 뒤의 값에 맞추어 rate를 정규화하지 않는다.

추가로 centered position \(\mathbf x=\mathbf r-(\mathbf l+\mathbf h)/2\)를 사용하면 다음이 성립한다.

\[
\begin{aligned}
\sum_jw_j\tau_jv_{j,i}&\to0,\\
\sum_jw_j\tau_jv_{j,i}v_{j,k}&\to nV\sigma^2\delta_{ik},\\
\sum_jw_j\tau_j|\mathbf v_j|&\to nV\sigma\sqrt{8/\pi},\\
\sum_jw_j\tau_j|\mathbf v_j|^4&\to15nV\sigma^4,\\
\sum_jw_j\int_0^{\tau_j}x_i(a)\,da&\to0,\\
\sum_jw_j\int_0^{\tau_j}x_i(a)x_k(a)\,da
&\to nV(L_i^2/12)\delta_{ik}.
\end{aligned}
\]

Chord midpoint를 \(\mathbf c=\mathbf x_{\rm in}+\mathbf v\tau/2\), displacement를 \(\mathbf d=\mathbf v\tau\)라 쓰면 spatial second moment의 age 평균은 정확히 \(\mathbf c\mathbf c^T+\mathbf d\mathbf d^T/12\)다. 공간 적분의 추가 sample error 없이 chord의 위치·길이를 검사한다. `occupation_moments()`의 반환값들은 수치 점유 합이 아닌 analytic \(nV\)로 나누므로 잘못된 점유수를 숨기지 않는다.

면별 총 rate는 노드 가중치 정의로 정확하므로, 그것만으로 속도 분포나 residence가 옳다고 주장할 수 없다. 부피를 유지하며 box의 aspect ratio를 바꾸면 \(nV\)는 같고 면적과 \(J_{\rm total}\)은 달라진다. 이 경우에도 점유수와 각 shape의 spatial moment가 함께 수렴해야 한다. 전체 크기를 바꾸는 검사에서는 rate가 길이의 제곱, residence가 길이, 점유수가 길이의 세제곱에 비례한다.

내부점을 부피 분포에서 먼저 선택한 뒤 backtrace하고 \(1/\tau\)로 rate를 정하는 다른 중요도 적분도 수학적으로 구성할 수 있다. 그러나 그 구성에서는 점유 identity가 정의상 자동으로 성립할 수 있고 작은 chord에서 reciprocal weights가 날카로워진다. 이번 구현은 직접 경계 flux를 표본화하여 \(nV\)를 독립 검사로 남긴다.

## 하나의 공통 lab RF 입사 위상

Lab-time drive의 beat frequency를 \(\omega_{\rm RF}\)라 하면

\[
\phi_{\rm in}=\omega_{\rm RF}t_{\rm in}+\phi_0.
\]

각 원자의 입사 위치는 path quadrature가 정한다. `entry_phase_factors(path, common_phase_rad, wavevectors_rad_m, temporal_harmonics, *, offsets_rad)`는 각 구동항에

\[
\exp\left[i\left(
\mathbf k_\ell\cdot\mathbf r_{\rm in}
+h_\ell\phi_{\rm in}+\theta_{\ell,0}\right)\right]
\]

를 반환한다. \(h_\ell\)는 같은 RF beat의 signed integer harmonic이고 \(\theta_{\ell,0}\)는 고정된 입력 offset이다. 경로에 따른 spatial phase를 따로 uniform randomize하지 않으며, 두 구동항에 서로 독립적인 temporal phase를 주지 않는다.

입사 후의 \(\mathbf k_\ell\cdot\mathbf v a+h_\ell\omega_{\rm RF}a\) age evolution은 호출자의 Hamiltonian에서 처리한다. 이 helper가 velocity-dependent drive evolution이나 원자 내부 dynamics를 대신 계산하지 않는다. `uniform_entry_phases(count, offset_rad=...)`는 공통 \(\phi\)에 대한 periodic trapezoid rule이고, `EntryPhaseQuadrature`는 명시적인 양의 확률과 출처를 갖는 규칙도 받을 수 있다. Period-averaged 결과로 해석하려면 이 규칙이 uniform lab-time phase integral을 근사해야 한다.

예를 들어 고정 entry position에서 conditional mean vector가

\[
\mathbf m(\phi)=
\begin{pmatrix}\cos(\phi+\mathbf q\cdot\mathbf r_{\rm in})\\
\cos(\phi+\mathbf Q\cdot\mathbf r_{\rm in})\end{pmatrix}
\]

이면

\[
\langle\mathbf m\mathbf m^T\rangle_\phi
=\frac12
\begin{pmatrix}
1&\cos[(\mathbf q-\mathbf Q)\cdot\mathbf r_{\rm in}]\\
\cos[(\mathbf q-\mathbf Q)\cdot\mathbf r_{\rm in}]&1
\end{pmatrix}.
\]

두 위상을 독립적으로 평균하면 이 비대각 상관을 잘못 없앤다. 공통 phase rule의 refinement와 offset 변화로 이런 상관과 phase-dependent wavepacket의 수렴을 검사한다.

## Marked-Poisson 합: raw moment를 먼저 평균한다

`average_marked_poisson(inflow, phases, wavepacket_factory, *, source)`의 callback은 `wavepacket_factory(BallisticPath, scalar_common_phase_rad)`다. 기존 `transport.wavepacket`과 같은 다음 schema가 필요하다.

| Key | 의미와 shape |
| --- | --- |
| `frequency_axis` | frame을 선언한 `GeneratorFrequencyAxis` |
| `greater`, `lesser` | 각 ordering의 connected finite-pulse covariance `[frequency, readout, readout]` |
| `mean_pulse` | finite-pulse conditional mean `[frequency, readout]` |
| `residence_time_s` | 전달된 path의 실제 chord residence |

모든 packet은 같은 frame, 주파수 grid, readout 순서·단위·coupling convention을 사용해야 한다. 구현은 grid, frame, shape, residence와 covariance의 Hermiticity·PSD를 검사한다. Readout의 물리적 단위와 coupling 일치는 호출자가 `source`로 설명해야 한다.

경로 \(j\)와 위상 \(k\)의 rate·확률을 각각 \(w_j,p_k\), conditional pulse mean을 \(\mathbf m_{jk}\), connected covariance를 \(C_{jk}^{>,<}\)라 쓰면 raw ordered pulse moment는

\[
R_{jk}^{>,<}=C_{jk}^{>,<}+\mathbf m_{jk}\mathbf m_{jk}^\dagger
\]

이고 결과는

\[
S_0^{>,<}(\Omega)
=\sum_j w_j\sum_k p_kR_{jk}^{>,<}(\Omega).
\]

Poisson의 서로 다른 사건에 해당하는 product of means는 stream의 connected covariance를 구성할 때 취소되지만, 같은 사건의 conditional mean outer product는 남는다. 기본 Poisson shot-noise의 mean·variance 관계는 [Illinois의 Lecture 4, §4.5 Campbell's theorem](https://maxim.ece.illinois.edu/teaching/spring16/notes/lec04.pdf)에 제시되어 있다. 위 식은 이 독립 사건 합을 complex vector, 내부 ordered covariance와 path/phase marks에 적용한 유도다.

반환값은 `internal_greater`, `internal_lesser`, `poisson_number`를 분리하며, `greater`와 `lesser`는 각각 내부 기여와 number 기여의 합이다. Number 기여는 양 ordering에 같아 ordering 차에는 기여하지 않는다. `arrival_weighted_conditional_mean_pulse`는 진단용 mark 평균이다. 그 outer product를 `poisson_number` 대신 쓰거나 총 spectrum에서 다시 빼면 안 된다.

고전적인 identity readout \(O=I\)만 있어도 한 원자의 내부 connected covariance는 0인 반면

\[
m_j(\Omega)=\int_0^{\tau_j}e^{i\Omega a}da
=\tau_j e^{i\Omega\tau_j/2}
\operatorname{sinc}\left(\frac{\Omega\tau_j}{2\pi}\right),
\qquad \operatorname{sinc}(x)=\frac{\sin(\pi x)}{\pi x}
\]

이므로 \(S_0=\sum_jw_j|m_j|^2\)는 일반적으로 양수다. 공통 phase-dependent readout에서도 \(\langle\mathbf m\rangle_\phi=0\)이면서 \(\langle\mathbf m\mathbf m^\dagger\rangle_\phi\ne0\)일 수 있다. 평균을 먼저 취하거나 Poisson mean 항을 버리는 대조군은 이 두 경우를 통과하지 못해야 한다.

주기적으로 구동되는 실제 lab process에서 \(\phi_{\rm in}\)은 도착 시각과 연결된다. Uniform entry-phase 평균은 deterministic periodic stream mean을 뺀 covariance의 **zero-cyclic, period-averaged frequency-diagonal spectrum**을 계산한다. 이는 모든 원자에 독립적인 random drive를 부여하는 물리 가정도, 전체 과정의 full stationarity 증명도 아니다. 서로 \(\ell\omega_{\rm RF}\)만큼 다른 주파수 사이의 cyclic spectra와 coherent periodic mean의 delta lines는 반환하지 않는다.

Fourier convention은 기존 finite packet과 같이 \(Y(\Omega)=\int e^{i\Omega a}O(a)\,da\), \(S(\Omega)=\int e^{i\Omega\tau}C(\tau)\,d\tau\)다. Input readout의 단위가 \(U\)라면 pulse mean은 \(U\,s\), covariance는 \(U^2s^2\), rate-weighted stream spectrum은 \(U^2s\)다. One-sided detector PSD, optical SQL 또는 squeezing dB로 환산하지 않는다. 양의 가중 합이 PSD인 것과 quantum-optical canonical commutator가 보존되는 것은 별개 검사다.

## Audit controls, 재현성과 수치 결과

Audit는 코드의 algebraic checks, analytic phase-space limits와 해상도·독립 seed 수렴을 구분한다. 다음 대조군은 오류를 의도적으로 만들어 올바른 구현과 구분되어야 한다.

1. **부피 속도를 입사 속도로 오용:** 정확한 \(J_f\)는 유지하되 normal Rayleigh를 half-normal \(u=\sigma\Phi^{-1}((1+U)/2)\)로 바꾼다. 면별 rate는 여전히 정확하지만 \(\sum w\tau\)와 Maxwell occupation moments가 잘못된다.
2. **Poisson conditional mean 삭제:** \(C=0\)인 identity 또는 phase-dependent mean pulse에서 올바른 유한 number spectrum을 0으로 잘못 만든다.
3. **독립 두 위상 평균:** 고정 spatial phase 관계를 가진 두 readout의 공통 RF 상관을 잘못 없앤다.

열평형 점유 모멘트의 수렴은 geometry·velocity ensemble의 검사다. 실제 Rb wavepacket이 entry position, phase 또는 velocity에서 더 빠르게 변하면 해당 spectrum 자체의 path·phase quadrature를 다시 수렴시켜야 한다. Synthetic covariance와 identity pulse 검사는 실제 Rb 내부 dynamics의 독립 QRT 비교를 대체하지 않는다.

Versioned [보고서](inflow_report_v1.json)와 [plot](inflow_v1.png)은 실행 조건, 수치 제어, 통과 기준, 대조군과 source hash manifest를 보존한다. Source/test hash를 실행 전후에 대조하고 해당 수치가 어느 코드에 속하는지 기록해야 한다. 이미 생성된 report/plot은 덮어쓰지 않고 새 version으로 재실행한다. 공유 작업 중 코드가 바뀌었다면 과거 manifest를 현재 코드의 검증으로 취급하지 않는다.

Audit는 \(T=373\) K, \(n=10^{17}\,\mathrm{m}^{-3}\), repository의 \(^{85}\mathrm{Rb}\) mass를 **가정한 입력값**으로 사용한다. 같은 \(2\,\mathrm{mm}^3\) 부피의 cube와 \(4\times2\times0.25\,\mathrm{mm}^3\) thin box를 비교하고, 면당 \(2^{10},2^{12},2^{14},2^{16},2^{18}\)개 노드와 seeds 11, 211, 811을 사용한다. 최종 각 quadrature는 총 1,572,864 paths다. Phase-space tolerance는 dimensionless 0.002이며, identity-pulse spectrum은 마지막 refinement 변화와 독립 seed spread 모두 1.5% 이내여야 한다. 별도로 \(H=0\), jump가 없는 실제 transport identity-readout ODE를 여섯 path에서 적분하여 finite-pulse analytic 식과 비교한다.

v1 audit의 모든 선언된 controls가 통과했다. 다음 값은 최종 면당 \(2^{18}\)개 노드에서 세 seed 중 가장 큰 오차다. Zero-target mean과 비대각 모멘트도 포함하므로, 모멘트 오차는 analytic 0으로 나눈 상대오차가 아니라 \(\sigma\), \(\sigma^2\), \(L_i/\sqrt{12}\) 등 해당 자연 scale로 나눈 무차원 오차다.

| Geometry | 최대 모멘트 무차원 오차 | 최대 점유수 상대오차 | \(2^{16}\to2^{18}\) identity spectrum 변화 | 최종 spectrum seed spread |
| --- | ---: | ---: | ---: | ---: |
| Cube | \(1.9172\times10^{-4}\) | \(5.2121\times10^{-5}\) | 0.6230% | 0.05425% |
| Thin box | \(3.2278\times10^{-4}\) | \(1.2150\times10^{-4}\) | 0.8771% | 0.5957% |

Spectrum 수렴은 DC, 50, 150, 500 kHz의 identity readout을 별도로 비교한 것이다. Seed spread는 seed별 값의 range를 seed 평균으로 나눈 값이다. 이 표는 검사한 형상·주파수·수치 규칙에서 관측된 오차이며, 연속 변수 전체에 대한 보장된 상한이나 실험 불확도가 아니다.

같은 analytic 점유수 \(nV=2\times10^8\)에서 thin box/cube의 면적 비와 총 유입률 비는 각각 1.994874995667로 일치했다. 두 형상의 수치 점유수 비는 0.9999076574다. 경계 면적에 따라 유입률과 평균 residence가 함께 바뀌어 같은 volume occupation을 회복한다.

잘못된 half-normal normal 속도는 정확한 면별 rate를 유지하면서도 점유수를 cube에서 \(1.14143194\,nV\), thin box에서 \(1.69658593\,nV\)로 만들었다. 이 대조군은 두 형상 모두 기각되었다. 서로 다른 면별 면적 가중치만 맞추는 것으로 flux-weighted velocity distribution을 검증할 수 없음을 보여 준다.

Synthetic 공통 phase control에서 1, 2개 phase는 불충분했고 4, 8, 16개 phase의 raw ordered matrices는 analytic phase integral과 최대 \(1.14\times10^{-15}\) 상대오차로 일치했다. Ordering difference의 상대오차는 \(5.48\times10^{-15}\) 미만이다. 평균을 먼저 취하고 outer product를 만든 경우, 두 phase를 독립 평균한 경우, Poisson mean 항을 버린 경우의 matrix relative error는 각각 1.0, 0.67545, 0.81057이다.

실제 `transport.wavepacket`으로 적분한 identity-readout 대조군은 여섯 characteristic, DC와 50 kHz에서 analytic pulse와 \(3.98\times10^{-16}\) 미만의 상대오차로 일치했고, 내부 connected 기여는 정확히 0이었다. 이는 기존 finite-atom packet과 신규 rate 합의 interface 검사이며 실제 Rb optical readout의 검증은 아니다.

보고서의 11개 source/test SHA-256은 audit 실행 전후에 같았다. 생성된 `inflow_report_v1.json`의 SHA-256은 다음과 같다.

~~~text
262ecbfe648ad45f45845e2b31f33396fdc972c59ef7181b5a46ec55c37d63bc
~~~

신규 출력 이름으로 audit 실행:

~~~powershell
python -m analysis.grand_challenge.inflow_audit --output docs/grand_challenge/inflow_report_v1.json --plot docs/grand_challenge/inflow_v1.png
~~~

이미 v1 출력이 있으면 두 출력 이름을 함께 v2 등 새로운 version으로 바꾼다.

Targeted 검사:

~~~powershell
python -m pytest -q tests/quantum/test_inflow.py
~~~

Repository 통합 검사:

~~~powershell
python -m pytest -q
~~~
