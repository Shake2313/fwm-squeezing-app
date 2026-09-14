# S1: 원자 diffusion에서 photon-flux field M/D와 구간 전파까지

2026-09-09. 구현: `gabes/quantum/traveling.py`, `gabes/fwm_quantum/field.py`. 전제인 단일 원자 정리는 [derivation.md](derivation.md), 수치 기록은 [s1_field_noise_report.json](s1_field_noise_report.json)에 있다.

후속 구현: conditional four-sideband Gaussian channel, normalized temporal mode와 bright-carrier intensity-difference/SQL은 [readout_derivation.md](readout_derivation.md)에 연결했다. 아래 field snapshot과 당시 적용 범위는 원래 기록으로 유지한다.

후속 정규화 감사에서 기존 weak-field의 공통 1/12가 uniform-manifold absorption의 두 F 평균을 동시에 재현하지 못함을 확인했다. 이전의 structural convention은 이론적으로 검증된 full-Zeeman reduction이 아니다. Pump/weak dipole을 함께 변경하는 명시적 variant와 독립 dipole reference는 [normalization_derivation.md](normalization_derivation.md)에 있다.

**이번 단계는 선택한 reduced double-Λ 광장 sector에서 M(Ω), 두 ordered noise matrices, distributed source-noise integral을 연결한다.** Independent atomic jump diffusion을 사용하며, M에서 필요한 noise를 역산하지 않는다. Local 및 propagated field commutator를 검사한다. 아직 검출기의 intensity-difference spectrum이나 실험 no-fit 예측을 출력하지 않는다.

**적용 범위와 근사**

원자에는 앞 단계의 finite-dimensional stationary Markov model을 사용한다. Optical field에 연결할 때는 주어진 classical pump를 기준으로 weak quantum fields에 대해 선형화하고, 원자 response/noise 계수를 pump-only state에서 평가한다. 원자 자체의 complete-basis two-point result와 달리, 이것은 전체 interacting atom-light system의 exact solution이라는 주장이 아니다. Generated fields의 back-action이 pump-only state를 바꾸는 영역은 적용 범위 밖이다.

현재 spatial model은 다음 조건을 명시한다.

- Independent atoms, 한 velocity class, 공통의 균일한 transverse area A, prescribed pump/state/density.
- 기존 reduced reference와 같은 선택된 probe/conjugate transitions 및 phase-selected two-mode sectors. 추가 optical modes, 다른 transition paths와 cross-sector scattering의 크기는 아직 평가하지 않았다.
- Retarded time에서 전파하며, 분석 RF 대역은 optical beat보다 매우 작다. Companion sector를 함께 계산한다. Bare geometric mismatch는 선택적으로 diagonal phase에만 더한다.
- 원자 reservoir는 선택된 traveling input fields와 독립인 Markov baths로 처리한다. Recurrent scattering, bath-mediated atom correlations, radiation trapping은 포함하지 않는다.

이것은 hot-vapor angular-Doppler·full Zeeman·Gaussian transverse-mode model의 완성이 아니다. 특히 공통 uniform area를 바꾸는 수학적 불변성과, 실제 Gaussian beam waist를 바꾸는 실험은 다르다. Collective atom-field Hamiltonian, propagation, microscopic stochastic noise를 잇는 선행 접근은 [Florez의 microscopic multimode FWM 연구](https://arxiv.org/html/2512.15051v1)에도 있다. 아래에서는 그 논문의 numerical factors를 가져오지 않고 GABES의 photon-flux 및 Rabi/2 규약에서 계수를 유도한다.

**Slice 평균과 원자수 정규화**

Number density를 n [m⁻³], common area를 A [m²], linear density를 λ=nA [m⁻¹]라 두자. 길이 dz의 slice에는 Nₛ=λ dz개의 원자가 있고,

\[
\overline F_i=\frac{1}{N_s}\sum_{a=1}^{N_s}F_i^{(a)}.
\]

독립 원자의 connected correlations에서는 a≠b 항이 사라지므로 slice-average covariance와 diffusion은 각각 `C_atom/Nₛ`, `D_atom/Nₛ`다. Continuum에서 atomic Langevin source ξ는

\[
\langle\xi_i(z,t)\xi_j(z',t')\rangle
=\frac{(D_{\rm atom})_{ij}}{\lambda}\delta(z-z')\delta(t-t').
\]

따라서 전파 구간을 세분화해도 결과에 임의의 dz 또는 총 cell atom number가 남아서는 안 된다. Atomic class v의 density가 nᵥ=n wᵥ라면 먼저 λᵥ=nA wᵥ를 사용해 class별 M,D를 만들고 **covariance를 합한다**. Noise amplitude를 평균한 뒤 제곱하지 않는다. `sum_independent_classes`는 이미 class density가 적용된 결과만 합하며 reservoir labels를 보존한다. 실제 velocity quadrature 생성과 그 수렴은 후속 공정이다.

**Reciprocal atom-field coupling**

광장 b의 각 성분은 annihilator 또는 creator다. J=diag(sⱼ), sⱼ=+1/−1로 두고 `[b_i(t),b_j†(t')]=J_ij δ(t-t')`를 사용한다. Main sector의 순서는 `(a_probe, a_conjugate†)`이며 physical mode labels와 signs를 객체에 저장한다.

Peak electric envelope와 photon flux의 기존 변환은

\[
E_j=Q_jb_j,\qquad Q_j=\sqrt{\frac{2\hbar\omega_j}{\epsilon_0 c A}},
\qquad g_j=\frac{d_{\rm eff}Q_j}{2\hbar}.
\]

gⱼ의 단위는 s⁻¹ᐟ²이다. `H_int/hbar=Σⱼ gⱼ bⱼ Oⱼ†+h.c.`로 쓰며, main sector의 O는 probe lowering과 conjugate raising이다. Transition별 dimensionless dipole scales는 O 안에 들어 있다. GABES의 Hamiltonian이 `Rabi/2`를 사용하는 것이 g 정의의 factor 2를 결정한다.

Complete atomic basis Fᵢ와 stationary ρ에 대해

\[
W_{ji}=\operatorname{Tr}(O_jF_i),\qquad
B_{ij}=-ig_j\operatorname{Tr}\{F_i[O_j^\dagger,\rho]\},\qquad
C_0=-iJ\operatorname{diag}(g)W.
\]

B는 density-operator driving commutator에서 계산한다. C₀는 동일 Hamiltonian에서 광장의 Heisenberg equation으로 얻는 방출 계수다. 공통 면적의 continuum atom-field equations는

\[
\partial_t\Delta\overline F=A_{\rm atom}\Delta\overline F+B b+\xi,
\qquad \partial_z b=\lambda C_0\Delta\overline F.
\]

평균장이 있는 경우 여기서 b는 그 평균 주변의 weak fluctuation 좌표다. Strong classical pump의 mean evolution은 이 식에서 자동으로 풀리지 않는다.

**독립 elimination과 field diffusion**

\[
R(\omega)=(-i\omega I-A_{\rm atom})^{-1},\qquad
M(\omega)=\lambda C_0R(\omega)B,
\]

\[
D_>(\omega)=\lambda C_0R D_{\rm atom}R^\dagger C_0^\dagger,
\qquad
D_<(\omega)=\lambda C_0R D_{\rm atom}^T R^\dagger C_0^\dagger.
\]

M,D₍>,<₎의 단위는 m⁻¹이다. Greater ordering은 `⟨fᵢ(ω) fⱼ(ω′)†⟩`, lesser ordering은 `⟨fⱼ(ω′)† fᵢ(ω)⟩`에서 `2π δ(ω−ω′) δ(z−z′)`를 제외한 행렬이다. 여기서 †는 Fourier-transformed operator의 adjoint를 뜻한다. 각 explicit reservoir의 두 D를 별도로 계산·저장한다.

Slice source의 `D_atom/λ`에 양쪽 field coefficient λC₀가 작용하므로 최종 계수는 λ다. Fixed n, pump Rabi, uniform geometry에서 `λgᵢgⱼ`의 area 의존성은 상쇄된다. 실제 pump power를 고정한 채 waist를 바꾸면 atomic state와 mode overlaps가 달라지므로 같은 결론을 사용할 수 없다.

Lesser ordering은 greater의 행렬 원소를 같은 주파수에서 transpose하는 것과 일반적으로 다르다. Independent QRT 검증에서는 `S_atom,<(ω)=S_atom,>(−ω)ᵀ`를 직접 계산하고 C₀로 projection한다. Atomic diffusion/R 구현을 공유하지 않는 앞 단계의 full-Liouville QRT solver를 사용한다.

**Local field commutator 정리**

원자의 expected commutator를 K=C_atom−C_atomᵀ로 두면 reciprocal coupling에서 `B J=−K C₀†`가 성립한다. 이 관계와 앞 단계의 `A_atom K+K A_atomᵀ+D_atom−D_atomᵀ=0`를 사용하면

\[
M J+J M^\dagger+D_>-D_<=0.
\]

실제로 `MJ+JM†=−λ C₀(RK+KR†)C₀†`이고, atomic identity를 resolvents로 감싸면 `RK+KR†=R(D_atom−D_atomᵀ)R†`이므로 위 등식이 따른다. D₍>,<₎는 각각 PSD다. 이 증명은 drift가 정해지면 arbitrary minimum noise를 선택하는 절차가 아니다. 이미 explicit jumps에서 독립 계산한 noise가 reciprocal field coupling과 맞는지 확인하는 정리다.

Implementation은 이 식의 signed numerical residual 및 reservoir별 positivity/Hermiticity를 검사한다. Negative control로 D₍>,<₎에만 추가로 1/12를 곱하면 검사가 실패한다. Source covariance를 clipping하거나 M에 맞추어 보정하지 않는다.

**기존 1/12와 reduced dipole의 의미**

기존 Maxwell convention은 `chi_phys=−2n d² s chi_bar/(epsilon_0 hbar)`, structural s=1/12다. 새로운 경로가 같은 weak-field normalization을 갖도록 report에서 `d_eff=d/sqrt(12)`를 **구동 B와 방출 C₀ 양쪽**에 선언한다. 따라서 M과 D 모두 같은 d_eff² scale을 갖는다. Transition strengths와 stationary population을 다시 곱하지 않는다. 기존 residual 0.74 및 squeezing efficiency는 사용하지 않는다.

이 선택은 기존 reduced response와의 연결 규약이다. Full Zeeman dynamics의 완전한 coarse-graining이나 실제 장치의 absolute dipole calibration을 증명하지 않는다. 특히 기존 `rabi_freq(power, waist)`로 정한 pump Rabi와 이 weak-field dipole을 하나의 microscopic power-to-Hamiltonian 규약으로 닫는 일은 남아 있다. Report는 pump Rabi를 조건부 입력으로 고정하며, 실측 입력이 준비되었다고 표시하지 않는다.

동일 explicit reservoirs와 pump Rabi를 기존 `pump_only_weak_response_reference`에 넣어 얻은 chi를 기존 `gain_from_chi`와 photon-flux 변환에 통과시킨 결과와 새 transfer를 대조한다. 이 parity에는 noise 구현을 공유하지 않는다. 동일 reduced physics에서의 software parity이며 실험적 타당성의 증거는 아니다.

**Main/companion과 주파수**

Lab RF를 Ω, minus-branch optical beat를 ω_b=−ω_hf+δ라고 두면

\[
\omega_{\rm main}=\omega_b+\Omega,\qquad
\omega_{\rm companion}=-\omega_b+\Omega.
\]

Companion에서는 O→O†, J→−J를 적용하고 별도로 계산한다. RF axis를 반사하면 `M_comp(Ω)=M_main(−Ω)*`, `D_comp,>(Ω)=D_main,<(−Ω)*`이고 반대 ordering도 같은 관계를 갖는다. Report는 −4…4 MHz의 81개 RF points와 두 generator axes를 모두 저장한다. Symmetric mismatch terms는 main에서 `(−iΔk/2,+iΔk/2)`, companion에서는 반대다. Refractive response는 M의 atomic diagonal에 이미 들어 있으므로 mismatch에 다시 더하지 않는다.

두 sector를 계산했다고 arbitrary broad-band real quadrature covariance가 자동으로 만들어지는 것은 아니다. Normalized temporal filters와 ±RF optical modes를 지정하는 adapter가 다음 단계다.

**구간 전파와 distributed noise**

길이 ℓ 안에서 M,D가 일정한 경우

\[
T=e^{M\ell},\qquad N_{>,<}=\int_0^\ell ds\,e^{Ms}D_{>,<}e^{M^\dagger s}.
\]

`constant_segment`는 row-major covariance equation의 Kronecker generator `M⊗I+I⊗M*`에 두 source columns를 붙인 matrix exponential로 이 적분을 평가한다. 이는 loss를 cell 뒤에 한 번 적용하는 근사가 아니다. 여러 구간을 first→second 순서로 합성하면

\[
T_{21}=T_2T_1,\qquad N_{21}=T_2N_1T_2^\dagger+N_2.
\]

각 propagation과 합성 뒤에

\[
TJT^\dagger+N_>-N_<=J
\]

및 두 added-noise matrices의 positivity를 검사한다. `vacuum_output`은 input greater/lesser `diag((1±sⱼ)/2)`를 전파한 ordered spectral covariances다. Detector PSD, finite-mode covariance 또는 squeezing dB를 반환하는 API가 아니다. Segment별 state/pump가 주어지면 전파할 수 있지만 self-consistent pump depletion은 아직 풀지 않는다.

**수치 검증과 재현**

Local/global commutator tolerance는 각 식의 항별 matrix norms 합으로 나눈 10⁻⁸이다. PSD/Hermiticity는 해당 ordering norm에 10⁻⁸을 곱한 값에 paired-ordering norm의 `64 epsilon` roundoff allowance를 더해 판정한다. 이는 vacuum에서 정확히 0인 ordering을 10⁻³³ 수준 projection roundoff로 나누는 문제를 피한다. 원래 signed eigenvalues와 matrices는 그대로 저장한다. Physical uncertainties와 이 numerical tolerances는 별개다.

해석적 two-level tests에서는 emission-only ground state가 `D_< = 0`, `N_> = 1−|T|²`를 주어 vacuum을 보존한다. 완전 inverted reset은 `D_> = 0`, `N_< = |T|²−1`을 주는 quantum-limited amplifier가 된다. 이 inverted fixture는 실제 hot Rb medium에 대한 모델 주장이 아니다.

그 밖에 density/dipole/common-area scaling, independent-class splitting, ±RF companion, no-atom/zero-length, 직접 numerical quadrature와 covariance integral, constant-segment semigroup, noncommuting prescribed-segment order 및 잘못된 noise coefficient의 거부를 검사했다.

Snapshot 조건은 nominal Δ/2π=0.9 GHz, δ/2π=−8 MHz, n=10¹⁸ m⁻³, uniform A=1.2×10⁻⁷ m², ℓ=12.5 mm, Δk=0이다. Pump-on cases는 기존 600 mW/530 μm helper의 Rabi를 조건부 값으로 사용한다. Reset도 앞 단계의 declared fixture다.

| 조건 | 기존 Maxwell transfer 상대 차이 | Field noise/QRT 최대 상대 차이 | Local commutator residual | Global commutator residual |
|---|---:|---:|---:|---:|
| Pump off + reset | 7.08×10⁻¹⁶ | 2.28×10⁻¹⁵ | 1.80×10⁻¹⁸ | 1.42×10⁻¹⁶ |
| Pump on, radiation only | 2.35×10⁻¹⁴ | 5.58×10⁻¹¹ | 7.52×10⁻¹⁴ | 2.40×10⁻¹⁴ |
| Pump on + reset | 2.75×10⁻¹⁴ | 3.26×10⁻¹¹ | 1.09×10⁻¹³ | 2.85×10⁻¹⁴ |

81 frequency points에서 평가한 matrix relative errors다. Companion 관계의 최대 오차는 3.81×10⁻¹³, 구간 분할의 최대 상대 차이는 5.60×10⁻¹⁶다. 세 경우 모두 declared gates를 통과했고 noise-only factor negative control은 거부되었다. 이것은 exhaustive operating-domain scan이 아니다.

```powershell
# New environments: install the optional research dependency set.
python -m pip install -r requirements-quantum.txt
python -m pytest -q tests/quantum
python -m analysis.grand_challenge.field_noise_audit --output NEW_FIELD_REPORT.json
python -m pytest -q
```

Matrix exponential에는 SciPy를 사용하며 version을 report에 저장한다. 기존 S0/atomic-only 보고서는 당시 source snapshot으로 유지한다. 새 report도 기존 파일을 덮어쓰지 않는다.

**Field 단계에서 정의한 후속 공정**

1. Normalized temporal-mode/±RF quadrature adapter, quantum uncertainty 검사와 mean-field-dependent intensity-difference/SQL readout을 연결한다.
2. Pump Rabi, weak-field dipole, optical carriers, density와 reset/transport 입력의 physical provenance를 통합한다. 이를 해결하기 전에는 이번 결과를 실험 no-fit absolute prediction으로 올리지 않는다.
3. Angular-Doppler/velocity quadrature와 transverse mode overlaps를 유도하고 독립 convergence를 검사한다. 그 뒤 mean-field back-action과 self-consistent segmentwise depletion을 추가한다.
