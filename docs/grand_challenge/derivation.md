# S1: explicit reservoirs에서 단일 원자 diffusion과 spectrum까지

2026-09-09. 구현: `gabes/quantum/diffusion.py`, reduced adapter: `gabes/fwm_quantum/model.py`. 독립 검증: `analysis/grand_challenge/reference/atomic_qrt.py`. 규약: [conventions.md](conventions.md). 수치 기록: [s1_atomic_noise_report.json](s1_atomic_noise_report.json).

후속 구현: 이 문서의 single-atom 결과를 이용한 conditional local field M/D, reciprocal coupling과 prescribed-segment propagation은 [field_derivation.md](field_derivation.md)에 기록했다. 아래 정리와 atomic snapshot은 원래 범위를 유지한다.

**이번에 닫은 문제는 명시한 Markov Lindbladian의 단일 원자 stationary two-point spectrum이다.** Drift와 jump diffusion을 독립적으로 계산하고, 별도의 density-operator QRT 적분과 대조한다. 아래 정리는 finite-dimensional 원자의 second moments에 대한 것이다. Optical propagation matrix M, effective field diffusion, intensity-difference squeezing은 아직 이 결과에 포함되지 않는다.

**가정과 좌표**

유한한 n-level 원자에 대해 시간에 무관한 generator를 선언한다. H는 에너지를 ℏ로 나눈 rad/s 값이며, 각 jump Lₖ에는 이미 sqrt(rate)가 들어 있다.

\[
\mathcal L\rho=-i[H,\rho]+\sum_k\left(L_k\rho L_k^\dagger-
\tfrac12\{L_k^\dagger L_k,\rho\}\right).
\]

입력은 trace-one positive stationary state ρₛₛ를 갖고, traceless subspace의 모든 모드는 엄격히 감쇠해야 한다. 구현은 이 조건의 수치적 분해 가능성도 검사한다. 고립된 비감쇠 coherence, 여러 steady states 또는 분해하지 못한 좁은 모드가 있으면 일반적인 연속 PSD를 반환하지 않는다. 이런 경우에는 initial-state dependence와 elastic delta terms를 별도로 다루어야 한다.

Hermitian traceless operator Fᵢ, i=1,…,n²−1를 `Tr(FᵢFⱼ)=δᵢⱼ`로 정규화한다. 첫 n−1개는 generalized diagonal generators이고, 나머지는 기존 `core.hermitian_basis`의 symmetric/antisymmetric off-diagonals다. Report에 실제 operator arrays를 저장하므로 이름만으로 순서를 추측할 필요가 없다. Row-major vec(Fᵢ)를 열로 놓은 B에 대해

\[
A=B^\dagger\mathcal L B,\qquad
b_i=\operatorname{Tr}[F_i\mathcal L(I/n)],\qquad
\frac{d\langle F\rangle}{dt}=A\langle F\rangle+b.
\]

A는 Hermiticity preservation 때문에 실수다. Complete basis이므로 이 closure에는 operator truncation이 없다. μᵢ=Tr(ρₛₛFᵢ), ΔFᵢ=Fᵢ−μᵢI로 중심화하면 adjoint evolution은 정확히 `L*(ΔFᵢ)=ΣⱼAᵢⱼΔFⱼ`가 된다. Stationary means와 A를 별도로 계산하며, 전체 finite atom을 Gaussian random variable이라고 가정하지 않는다.

**정리: ordered diffusion의 독립 구성과 stationary identities**

위 가정에서 다음 두 행렬을 정의한다.

\[
C_{ij}=\operatorname{Tr}(\rho_{ss}\Delta F_i\Delta F_j),\qquad
D_{ij}=\sum_k\operatorname{Tr}\left(\rho_{ss}
[L_k^\dagger,F_i][F_j,L_k]\right).
\]

각 reservoir의 diffusion D⁽ᵏ⁾와 C는 Hermitian positive semidefinite이며,

\[
AC+CA^T+D=0.
\]

증명은 jump-level product rule에서 시작한다. Hamiltonian commutator는 Leibniz rule을 만족하므로 아래 차이에서 상쇄된다. 각 dissipator를 전개하면

\[
\mathcal L^*(F_iF_j)-(\mathcal L^*F_i)F_j-F_i(\mathcal L^*F_j)
=\sum_k[L_k^\dagger,F_i][F_j,L_k].
\]

Kᵢ⁽ᵏ⁾=[Fᵢ,Lₖ]로 두면 우변의 기대값은 `Tr(ρₛₛ Kᵢ†Kⱼ)`다. 임의 복소 vector z에 대해 `z†D⁽ᵏ⁾z=Tr(ρₛₛ Q†Q)≥0`, `Q=ΣⱼzⱼKⱼ`이므로 Gram positivity가 따른다. 같은 논리가 ΔFᵢ에 적용되어 C≥0를 준다. Stationarity에서 `Tr[ρₛₛ L*(ΔFᵢΔFⱼ)]=0`이고, centered closure와 product rule을 대입하면 Lyapunov identity가 나온다.

Production code는 D를 jump commutators에서 직접 계산한다. `D=−AC−CAᵀ`는 검증에만 사용한다. 별도 reference는 full Liouville adjoint를 product에 작용시켜 Einstein 식의 세 항을 평가한다. 이 reference는 jump diffusion 구현을 import하지 않으며, 비정상상태의 임의 physical ρ에서도 두 product-rule 계산을 비교한다. 그 비교에는 stationarity가 필요하지 않다.

Atomic expected commutator matrix를 `K=C−Cᵀ`, 즉 `Kᵢⱼ=⟨[Fᵢ,Fⱼ]⟩`로 두면

\[
AK+KA^T+D-D^T=0.
\]

이 식을 별도로 검사하며 imaginary cross correlations를 버리지 않는다. 이것은 **stationary expected atomic commutator balance**다. Fᵢ에는 상태와 무관한 canonical bosonic J를 부여하지 않는다. Field의 input/output commutator preservation을 증명하는 일은 후속 단계다.

**유한 주파수 spectrum과 ordering**

시간 Fourier inverse에 exp(−iΩt)를 사용한다. Centered correlation과 two-sided angular-frequency spectrum의 정의는

\[
C_{ij}(\tau)=\langle\Delta F_i(\tau)\Delta F_j(0)\rangle_{ss},
\qquad S_{ij}(\Omega)=\int_{-\infty}^{\infty}d\tau\,
e^{+i\Omega\tau}C_{ij}(\tau).
\]

QRT에 의해 τ≥0에서 `C(τ)=exp(Aτ)C`, τ≤0에서는 `C(τ)=C exp(−Aᵀτ)`다. A가 strictly stable이므로 양쪽 적분이 수렴한다. `R(Ω)=(−iΩI−A)⁻¹`로 두면

\[
S_{\rm ord}(\Omega)=R(\Omega)C+C R(\Omega)^\dagger
=R(\Omega)D R(\Omega)^\dagger.
\]

마지막 등식은 Lyapunov identity에 양쪽 resolvents를 곱하면 얻는다. D≥0이므로 S_ord≥0이다. 이는 field diffusion D_field(Ω)를 유도한 식이 아니다. 여기서 atomic D는 white Markov diffusion이며 frequency dependence는 atomic response R에서 생긴다.

Symmetrized spectrum은 시간상 anticommutator에 대응한다.

\[
S_{\rm sym}(\Omega)=\frac{S_{\rm ord}(\Omega)+S_{\rm ord}(-\Omega)^T}{2}
=R(\Omega)\frac{D+D^T}{2}R(\Omega)^\dagger.
\]

단일 주파수에서 S_ord의 실수부를 취하는 방법은 이 ordering 변환과 일반적으로 다르다. 코드와 테스트는 복소 cross spectrum과 ±Ω를 함께 유지한다.

F, C는 dimensionless, A와 D는 s⁻¹, R와 S는 seconds다. 정규화는 `C=∫S_ord(Ω)dΩ/(2π)`이다. 이 atomic PSD를 곧바로 detector의 one-sided A²/Hz 또는 SQL-normalized S₋라고 부를 수 없다.

**독립 QRT reference와 해석적 기준**

Reference는 위 atomic A, D, R을 사용하지 않는다. Full n²-dimensional Liouvillian에서 source `ΔFⱼρₛₛ`를 만들고

\[
T_{ij}(\Omega)=\operatorname{Tr}\left[
\Delta F_i(-\mathcal L-i\Omega)^{-1}_{\operatorname{Tr}=0}
(\Delta F_j\rho_{ss})\right],\qquad S=T+T^\dagger
\]

를 직접 계산한다. 역행렬은 trace constraint를 포함한 bordered linear solve로 평가하므로 DC의 stationary pole도 처리한다. Centered source의 trace roundoff만 ρₛₛ 방향으로 제거한다. 독립이라는 말은 계산 경로가 다르다는 뜻이며, physical generator·state·operators 자체는 공유한다. 따라서 두 계산의 일치만으로 Hamiltonian의 실험적 정확성을 증명하지는 않는다. Density-operator QRT의 기본 정의는 [QuTiP correlation documentation](https://qutip.readthedocs.io/en/stable/guide/guide-correlation.html)에 설명되어 있다. 해당 문서의 `spectrum`은 exp(−iωτ)를 쓰므로 우리 Ω와 비교할 때 ω=−Ω로 변환해야 한다. 이 저장소 reference는 QuTiP에 의존하지 않는다.

주파수 부호와 normalization의 별도 기준은 ground-state two-level atom이다. `H=ω₀|e⟩⟨e|`, `L=√γ|g⟩⟨e|`, `Fₓ=(|g⟩⟨e|+|e⟩⟨g|)/√2`에 대해

\[
C_{xx}=\tfrac12,\qquad
S_{xx}(\Omega)=\frac{\gamma/2}{(\gamma/2)^2+(\Omega-\omega_0)^2}.
\]

테스트는 +ω₀의 peak, width, absolute magnitude, symmetrized ±ω₀ peaks, 수치 적분과 해석적 tails의 합을 확인한다. Driven two-level cases는 DC 및 복소 cross spectra를 포함해 QRT와 비교한다.

**GABES reduced four-level 연결**

`reduced_pump_noise`는 기존 `atoms.double_lambda_rb85(gamma_gg=0)`의 radiative channels와 기존 pump Hamiltonian을 사용한다.

```python
H = fwm.pump_hamiltonian_at_Deff_zero(pump_rabi, pump_rabi)
H -= effective_one_photon_detuning * diag(0, 0, 1, 1)
```

Optional thermal replacement는 `Jᵢⱼ=√(γₜpᵢ)|i⟩⟨j|`, `p=(5/12,7/12,0,0)`로 선언한다. Pump saturation에 의한 stationary population/coherence 변화는 이 원자 모델 안에서 풀지만 pump는 주어진 classical field이며 seed back-action과 depletion은 포함하지 않는다. Report의 reset rate `constants.GAMMA_GG`는 조건부 fixture다. 그 값을 독립 측정한 transport law나 검증된 collision rate라고 승인한 것이 아니다.

기존 `pump_only_weak_response_reference`에도 **동일한 explicit reservoirs**를 제공하여 ρₛₛ를 비교한다. 이 검사는 production default의 모든 dephasing/collision physics를 그대로 채택했다는 뜻이 아니다. Legacy coherence-only dephasing의 import 거부는 유지한다.

Lab RF와 static-pump atomic frequency는 별도 type이다. Standard minus branch의 atomic sector에서는

\[
\Omega_{\rm atom}=-\omega_{\rm hf}+\delta+\Omega_{\rm RF}.
\]

`minus_branch_atomic_frequencies`는 `OpticalDetunings`와 `AnalysisFrequencyAxis`를 받아 명시적 `GeneratorFrequencyAxis`를 반환한다. Positive lab RF가 negative generator frequency로 옮겨지는 것은 frame 변환이다. 이것만으로 ±Ω optical sideband/Nambu blocks, field readout 또는 measured S₋가 구성되지는 않는다.

**수치 관문과 재현 결과**

Signed 최소 고유값과 residual을 저장하고 covariance/PSD clipping을 하지 않는다. Stability gap은 `64 ε n² ||L||₂`보다 커야 한다. Stationarity relative residual은 10⁻¹⁰ 이하, Lyapunov 및 expected commutator residual은 `||D||`로 나누어 10⁻⁸ 이하를 요구한다. Diffusion Hermiticity와 상대 PSD tolerance는 10⁻¹⁰이다. C positivity의 absolute tolerance는 normalized basis에서 10⁻¹⁰이다. Spectrum은 각 Ω의 matrix norm을 기준으로 Hermiticity/positivity를 검사한다. 이 허용치는 물리적 오차 범위나 모든 operating point에서의 정확도 보증을 뜻하지 않는다.

Snapshot은 Δ/2π=0.9 GHz, δ/2π=−8 MHz, nominal pump 600 mW 및 radius 530 μm를 사용한다. 각 case는 generator Ω/2π=−4,−1,0,1,4 MHz와 lab RF 0.1–4 MHz의 40개 shifted frequencies에서 검증했다. 아래 QRT 오차는 각 주파수의 전체 matrix Frobenius relative error 중 최댓값이다.

| 조건 | Jump/adjoint Einstein 오차 | Spectrum/QRT 최대 오차 | Lyapunov residual |
|---|---:|---:|---:|
| Pump off + declared reset | 2.09×10⁻¹⁶ | 3.50×10⁻¹³ | 2.46×10⁻¹⁶ |
| Pump on, radiation only | 3.54×10⁻¹⁴ | 6.59×10⁻¹² | 8.97×10⁻¹⁴ |
| Pump on + declared reset | 3.35×10⁻¹⁴ | 6.65×10⁻¹² | 1.40×10⁻¹³ |

모든 case가 선언한 numerical gates를 통과했다. Semidefinite cases의 작은 음의 고유값도 원래 값으로 report에 남긴다. Pump-reference state 차이는 Frobenius norm으로 최대 3.68×10⁻¹⁴다. 이 수치는 세 조건·두 축에서의 구현 비교이며 실험 오차나 exhaustive domain validation이 아니다.

```powershell
python -m pytest -q tests/quantum
python -m analysis.grand_challenge.atomic_noise_audit --output NEW_REPORT.json
```

보고서는 state, basis, A, C, reservoir별 D, ordered/symmetrized spectra, frequency frames, signed diagnostics, source hashes와 environment를 저장한다. 기존 파일을 덮어쓰지 않는다. 이전 S0 report는 당시 source snapshot으로 유지한다.

**Atomic 단계에서 정의한 후속 연결**

1. 동일한 explicit model의 weak-field driving/readout matrices를 photon-flux convention으로 유도하고 기존 pump-only response와 비교한다.
2. 원자 밀도·coarse-graining volume·mode area를 명시하여 atomic Langevin sources를 collective polarization과 field noise로 옮긴다. Velocity-class covariance와 density weights를 이 단계에서 도출한다.
3. ±Ω companion blocks를 포함해 local field M(Ω), D_field(Ω)와 reservoir commutator를 검증한 뒤 distributed propagation을 합성한다.
4. Mean fields, shot-noise normalization, collection/detection을 연결한 이후에 gain과 S₋를 출력한다. 현재 Grand Challenge milestone 1은 계속 진행 중이다.
