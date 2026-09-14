# Quantum foundations: gabes-quantum-v1

2026-09-09. S0 foundations와 S1 atomic/field/readout 구현에 적용하는 규약이다. Explicit reservoirs, input metadata audits, atomic diffusion/QRT, reduced field M/D와 prescribed-segment propagation, four-sideband Gaussian channel 및 conditional bright-carrier intensity-difference/SQL이 구현되어 있다. Atomic 유도는 [derivation.md](derivation.md), field 유도는 [field_derivation.md](field_derivation.md), readout과 제한은 [readout_derivation.md](readout_derivation.md)에 있다. Independent-input absolute hot-vapor prediction과 실험 검증은 남아 있다.

**원자와 주파수**

- Atomic state는 `Tr rho=1`이며 density vector는 기존 core와 동일한 row-major `vec(rho)[n*i+j]=rho[i,j]`이다.
- `ExplicitReservoirs.generator(H)`의 H는 에너지 Hamiltonian을 ℏ로 나눈 값이며 rad/s 단위다. Evolution은 `d rho/dt = -i[H,rho] + D[rho]`이다.
- `CollapseChannel.operator`는 이미 sqrt(rate)를 포함하는 jump operator로 s⁻¹ᐟ² 단위다. Population/decay rates는 s⁻¹이다. 주파수 Hz 값을 decay rate에 무조건 2π배 하지 않는다.
- `OpticalDetunings`는 한 operating point의 one/two-photon detuning을 rad/s로 보관한다. `AnalysisFrequencyAxis`는 독립 RF 좌표이며 이 객체와 교환할 수 없다.
- RF 축은 유한하고 중복 없는 증가순 1-D vector다. 음·양 Ω와 DC를 허용한다. `AnalysisFrequencyAxis.from_hz`에서 Ω=2πf를 한 번 변환한다.
- Fourier convention은 `f(t)=∫dΩ/(2π) f(Ω) exp(-iΩt)`이다. 시간 미분은 −iΩ, field response의 보존 subspace는 기존 trace-zero core 방식으로 처리한다.
- Optical carrier와 Doppler-shifted atomic detuning, RF Ω를 별도 변수로 유지한다. `GeneratorFrequencyAxis`는 generator-frame rad/s vector와 frame label을 요구하며 lab RF type과 교환할 수 없다. S1의 minus-sector adapter는 `Omega_atom=-OMEGA_HF+delta+Omega_RF`를 사용한다. Companion은 반대 optical beat에 같은 RF를 더한다. 동일 explicit reservoirs를 제공한 기존 pump-only reference의 state와 photon-flux Maxwell transfer parity를 검증했다. Readout은 positive RF와 reflected partner를 함께 사용하고 coherent carrier는 별도의 명시적 RF=0 row에서만 계산한다.

**Single-atom ordered noise**

Complete traceless Hermitian basis는 `Tr(F_i F_j)=delta_ij`이며 identity는 stationary means에서 처리한다. `AtomicNoiseModel`의 real A는 atomic drift다. `C_ij=Tr(rho_ss Delta F_i Delta F_j)`와 reservoir별 `D_ij=Tr(rho_ss [L†,F_i][F_j,L])`는 ordered covariance/diffusion이다. D는 A나 C에서 역산하지 않고 explicit jumps에서 직접 얻는다. Atomic A,D의 단위는 s⁻¹이고 C는 dimensionless다.

```text
R(Omega_atom) = (-i Omega_atom I - A)^-1
S_ordered = R D R†
S_symmetrized = R ((D+D.T)/2) R†
C = integral dOmega_atom/(2 pi) S_ordered(Omega_atom)
```

Atomic spectrum의 단위는 seconds이며 symmetrization에는 ±Ω의 ordering이 들어간다. `S_sym(Omega)=(S_ordered(Omega)+S_ordered(-Omega).T)/2`이고, 일반적으로 단일 frequency의 실수부를 취하는 것과 다르다. 이 atomic spectrum은 detector squeezing이 아니다.

`A C+C A.T+D=0`, `K=C-C.T`에 대한 `A K+K A.T+D-D.T=0` 및 C,D,S의 Hermiticity/positivity를 검사한다. Atomic K는 state-dependent expected commutator다. 아래 finite-mode channel의 canonical J와 혼용하지 않는다. Strictly decaying traceless modes가 없으면 elastic terms가 필요한 경우로 거부한다. 유도, independent adjoint/QRT reference, analytic Lorentzian 기준과 numerical tolerances는 [derivation.md](derivation.md)를 따른다.

**정규화와 readout 계약**

Photon-flux field는 `E=Q a`, `Q_j=sqrt(2 ℏ ω_j/(ε₀ c A_j))`이며 기존 `gabes.observables.photon_flux_mode_matrix`의 convention을 사용한다. Gaussian intensity 1/e² radius w에 대해 기존 effective area는 `A=πw²/2`다. 새 field bridge는 공통 uniform area를 별도 입력으로 받는다. 이 constant-mode derivation에 서로 다른 Gaussian waists를 단순 대입해 transverse overlap까지 계산한 것으로 취급하지 않는다. Slice-average atomic diffusion의 `D_atom/(n A dz)`에서 field normalization을 유도한다. S0 channel fixture에 원자수 인자를 임의로 곱하지 않는다.

**Complex Nambu field channels**

`LocalNambuGenerator`는 frequency frame, ordered physical mode labels, annihilator/creator signs와 reservoir별 greater/lesser D를 보관한다. Main sector는 `(a_probe,a_conjugate†)`, signs는 `(1,−1)`이며 companion의 signs는 반대다. `g=d_eff Q/(2 hbar)`를 weak-field driving과 emission 양쪽에 쓰고 `M=n A C0 R B`, `D_>=n A C0 R D_atom R† C0†`, `D_<=n A C0 R D_atom.T R† C0†`를 각각 계산한다. D를 M으로부터 역산하지 않는다. M,D의 단위는 m⁻¹이다.

Local check는 `M J+J M†+D_>−D_<=0`, propagated check는 `T J T†+N_>−N_<=J`다. 각 noise ordering의 PSD/Hermiticity도 검사한다. Constant-segment noise는 `integral exp(M s) D exp(M† s) ds`이며 `compose_segments(first,second)`의 순서는 first→second다. Mode labels, axes, frames 또는 signs가 다르면 명시적 변환 없이 합성하지 않는다.

`sum_independent_classes`는 class density를 이미 적용한 M과 covariance를 합한다. Atomic bath/class independence를 깨는 collisions는 이 함수로 표현되지 않는다. 모든 source matrices, signed residuals와 고유값을 유지한다. Roundoff tolerance와 model assumptions는 [field_derivation.md](field_derivation.md)에 고정했다.

이 complex spectral path는 아래 real finite-temporal-mode `GaussianChannel`과 다른 객체다. `sideband_channels`는 main/companion과 reflected RF의 일관성을 검사한 뒤 positive RF마다 `(probe:+, conjugate:-, probe:-, conjugate:+)` 네 모드의 real Gaussian channel을 만든다. 현재 phase-selected closure에서 없는 cross-sector moments를 추가 생성하지 않는다. Prescribed-segment propagation은 pump depletion을 푸는 기능도 아니다.

`G_probe=P_probe,out/P_seed,in`, `G_conjugate=P_conjugate,out/P_seed,in`이라는 power 정의와 photon-flux gain을 구분한다. Conjugate conversion에는 ω_conjugate/ω_probe가 필요하다. Weak-field transfer로부터 얻는 gain과 finite-seed nonlinear mean-field ratio가 같다는 가정도 검증 대상이다.

Readout의 계약은 `S_minus(f)=PSD[i_p-g*i_c]/PSD_SQL`이다. `DetectorResponse`는 두 arm의 intensity transmission×QE, complex current response h(f) [A/A], 선언한 balance g, independent additive electronics PSD와 source를 보관한다. Means에는 sqrt(eta), covariance에는 동일한 vacuum-loss channel을 적용한다. SQL과 quantum PSD에는 같은 h,g,mean currents를 사용한다. Detector RF axis와 covariance RF axis가 다르면 거부한다.

`intensity_difference_spectrum`은 bright-carrier approximation에서 positive-frequency **one-sided A²/Hz**를 계산한다. `S_Hz,two(f)=S_omega(2pi f)`이고 one-sided는 그 두 배다. `SQL=2 e (|h_p|² I_p+g²|h_c|² I_c)`이며 ratio는 dimensionless, dB는 `10 log10(ratio)`다. Electronics는 따로 저장하고 quantum PSD에 더한 total ratio를 반환한다. 암묵적 subtraction은 없다. RF=0 detector PSD, quadratic fluctuation photocurrent, spontaneous mean powers와 실제 RBW/VBW 장치 응답은 이 구현에 포함되지 않는다.

**Finite temporal-mode channel**

`GaussianChannel`은 지정한 normalized temporal modes에 대한 real quadrature channel이다.

```text
r = (x_probe, p_probe, x_conjugate, p_conjugate, ...)
x = (a+a†)/sqrt(2), p = (a-a†)/(i*sqrt(2))
[r_i,r_j] = i J_ij
J = direct_sum([[0,1],[-1,0]])
V_ij = <{Delta r_i,Delta r_j}>/2
V_vac = I/2
V_out = X V_in X.T + Y
```

`input_modes`와 `output_modes`는 quadrature 순서를 결정하는 이름 tuple이다. 같은 이름의 순서가 다른 channels는 명시적 변환 없이 합성할 수 없다. 관측하지 않는 모드를 버리는 rectangular X도 허용한다. 복소 RF/Nambu matrix를 이 class에 직접 넣으면 거부한다. 현재 two-mode closure의 ±RF adapter는 `sidebands.py`에 구현했으며 arbitrary Floquet/multimode conversion은 아직 아니다.

`TopHatBand`는 positive center f₀, width B<2f₀, Gauss quadrature order를 명시한다. Frequency filter는 h=1/sqrt(B), `integral |h|² dOmega/(2pi)=1`이다. Flat matched filters와 vacuum orthogonal inputs에서 `average_spectral_channels`는 `Xbar=<X>`, `Yeff=<Y>+0.5<(X-Xbar)(X-Xbar).T>`를 사용한다. 추가 항은 spectral mode leakage의 vacuum noise다. 이는 normalized L² temporal mode이며 strictly finite-time window라는 뜻은 아니다. Band-integrated current variance에는 별도로 B를 곱한다. 유도와 convergence는 [readout_derivation.md](readout_derivation.md)에 기록한다.

채널 합성의 순서는 `compose_channels(first, second)`에서 first → second다.

```text
X_total = X_second X_first
Y_total = X_second Y_first X_second.T + Y_second
```

X는 요동의 전달행렬이다. Mean field를 이 객체에서 자동 계산하지 않는다. Saturation/depletion을 포함한 평균장은 다음 단계의 별도 nonlinear solution으로 유지한다.

양자 일관성 검사는 다음과 같다.

```text
Y >= 0
Y + (i/2)(J_out-X J_in X.T) >= 0
V + (i/2)J >= 0
```

Generator/channel 감사는 signed 최소 고유값과 numerical tolerance를 반환하고 eigenvalue를 clipping하지 않는다. Covariance 적용과 channel 합성은 invalid channel을 거부한다. `vacuum_attenuator`의 transmission은 amplitude가 아닌 intensity transmission이며 `X=diag(sqrt(eta))`, `Y=diag((1-eta)/2)`다. 이 예제는 analytic/device vacuum channel로, pumped-vapor diffusion의 유도가 아니다. 관련 formalism은 [Gaussian Quantum Information](https://arxiv.org/abs/1110.3234)을 따른다.

**Explicit atomic reservoirs**

`CollapseChannel`은 이름, jump operator, kind와 source를 요구한다. `ExplicitReservoirs`는 duplicate names와 dimension mismatch를 거부한다. 기존 `AtomModel`에서 가져올 때는 nonzero legacy coherence dephasing을 거부하고, 명시적 channels의 합과 실제 `atom.lindblad` 사이 assembly parity를 검사한다. 숨은 superoperator 수정도 이 단계에서 드러난다.

새 radiation 기준 모델은 다음처럼 생성한다.

```python
from gabes import atoms
from gabes.quantum.reservoirs import ExplicitReservoirs

model = ExplicitReservoirs.from_atom(
    atoms.double_lambda_rb85(gamma_gg=0),
    source="GABES D1 CF2 branching and natural decay data",
)
```

Thermal replacement는 `J_ij=sqrt(gamma*p_i)|i><j|`를 사용한다. Σp=1인 diagonal target에 대해 `ΣD[J_ij](rho)=gamma*(rho_th*Tr(rho)-rho)`이며, population뿐 아니라 coherence와 trace가 1이 아닌 operator에도 선형적으로 성립한다. 이를 기존 thermal-reset superoperator 및 직접 operator action과 대조한다. 이 동등성은 transport rate의 실험적 정확성을 뜻하지 않는다.

`audit_generator`는 trace preservation, Hermiticity preservation과 conditional complete positivity를 검사한다. `J(L)=Σ_ij |i><j|⊗L(|i><j|)`인 unnormalized Choi와 `P=I-|Omega><Omega|`, `|Omega>=vec(I)/sqrt(n)`에 대해 `P J(L) P>=0`를 검사한다. [CCP criterion](https://arxiv.org/html/2409.17072v2)에 대응한다. CCP tolerance는 projected dissipative scale과 floating-point cancellation floor를 사용한다. 큰 Hamiltonian norm 자체에 상대허용치를 곱해 작은 dissipative negativity를 숨기지 않는다.

총 generator의 CP 통과와 각 reservoir의 물리적 의미는 별도다. 121 °C의 기존 combined generator가 통과하더라도 nonzero legacy dephasing을 가진 모델의 explicit import는 거부된다. 새 collision model은 물리적 provenance와 모든 coherence에 대한 dynamics를 제시해야 한다.

**독립 입력 metadata**

`ParameterEvidence`에는 value, unit, standard uncertainty, source ID, estimation method, 사용한 dataset/observable, applicability, covariance-group 이름과 상태를 기록한다. 상태는 `independent`, `assumed`, `target_fitted`, `unknown` 중 하나다.

`audit_independent_inputs`는 required input IDs를 요구하며 누락·중복·target dataset 재사용·불완전한 provenance를 검출한다. 독립 beam-profile 데이터에 Gaussian fit을 한 waist는 허용하지만 FWM target dataset으로 정한 waist는 거부한다. 0 uncertainty는 명시된 exact input에 사용할 수 있다. Calibration의 수학적 fitting 여부와 target 데이터에 대한 의존성을 혼동하지 않는다.

이 검사는 제출된 metadata의 일관성만 확인한다. source 내용의 진실성·장치 적용 가능성을 자동으로 입증하거나 전체 파라미터 의존성 목록을 만들어 주지 않는다. 향후 physical prediction pipeline은 모든 소비 입력을 required IDs로 연결해야 한다. S1 atomic audit의 nominal pump/reset inputs는 구현 검증용 조건부 fixture로 표시하며 독립 실험 입력으로 승인하지 않는다. Uncertainty propagation, correlated parameter sampling, experiment-validation badge는 후속 구현이다.

**실행과 현재 검증 범위**

```powershell
python -m pytest -q tests/quantum
python -m analysis.grand_challenge.foundation_audit
python -m analysis.grand_challenge.foundation_audit --output docs/grand_challenge/s0_foundation_report.json
python -m analysis.grand_challenge.atomic_noise_audit --output NEW_ATOMIC_REPORT.json
python -m analysis.grand_challenge.field_noise_audit --output NEW_FIELD_REPORT.json
python -m analysis.grand_challenge.readout_audit --output NEW_READOUT.json --plot NEW_READOUT.png
```

`--output`은 기존 파일을 덮어쓰지 않는다. 재실행 저장에는 새 이름을 사용한다. 보고서는 source hashes, environment, signed generator/channel residuals 및 scope를 보존한다. Legacy default가 거부되는 것은 의도한 negative control이며 보고서의 `expected_controls_passed`는 모든 generator가 유효하다는 뜻이 아니다.

현재 tests는 reservoir/operator action, atomic Einstein/QRT, local/global field commutator, Maxwell parity, passive/amplifier limits와 covariance propagation을 검증한다. Readout에는 four-sideband CP/uncertainty, ideal bright-seed `1/(2G-1)`, equal/unequal detector losses, complex response/balance, same-current SQL, 직접 Nambu PSD, finite-band vacuum leakage와 quadrature convergence를 추가했다. SciPy dependency는 `requirements-quantum.txt`에 선언한다. Conditional spectrum 계산이 가능하지만 실험 squeezing이나 bandwidth를 검증한 것은 아니다.
