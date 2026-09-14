# S1 — bright-carrier readout의 생략 항과 유한 seed 진단

2026-09-09. [정규화 감사](normalization_derivation.md)의 uniform-Zeeman-RMS 조건부 모델을 사용한다. **현재 8 μW seed에서 Gaussian quadratic readout 보정은 최대 1.033×10⁻⁵ dB, 국소 polarization의 finite-seed 변화는 최대 9.738×10⁻⁵이다.** 이 결과는 아래에 선언한 광학 수집 대역과 근사에 대한 진단이다. 전체 hot-vapor 오차나 실험 −7.8 dB를 검증한 결과는 아니다.

최종 [기계 판독 보고서](s1_seed_validity_report_v2.json)와 [검사한 그림](s1_seed_validity_v2.png). 첫 [보고서](s1_seed_validity_report.json)는 격자 수렴 기준을 통과하지 못한 이력으로 남긴다.

**1. 무엇을 보충했는가**

출력 광장을 aⱼ(t)=βⱼ+dⱼ(t), ⟨dⱼ⟩=0으로 쓰면 photon flux는

\[
\hat I_j=|\beta_j|^2+\beta_j^*d_j+\beta_jd_j^\dagger+d_j^\dagger d_j.
\]

이전 readout은 가운데 두 항의 noise와 coherent mean |β|²를 사용했다. 이번에는 다음을 추가한다.

- Φⱼ=⟨dⱼ†dⱼ⟩에 의한 spontaneous mean, mean photocurrent와 SQL 증가.
- dⱼ†dⱼ−Φⱼ의 quadratic photocurrent spectrum과 두 beam 간 intensity correlation.
- Pump-only 상태가 finite seed에서도 적절한지 확인하는 국소 periodic mean-state 계산.

앞의 두 항은 **Gaussian 광장이라는 추가 가정** 아래 계산한다. Gaussian state에서는 네 연산자 평균이 두 연산자 contractions로 정해진다. 선형 항과 quadratic 항 사이의 covariance는 zero-mean 세 연산자 평균이므로 0이다. Direct detection에서 normal/symmetric ordering을 구분해야 한다는 배경은 [Najafabadi et al., Intensity correlations in the Wigner representation](https://pmc.ncbi.nlm.nih.gov/articles/PMC11667585/)을 참고한다. 아래 두 beam의 연속 주파수 convolution은 이 저장소의 Fourier convention에서 직접 유도하고 별도 Fock-space 계산과 대조했다.

Atomic diffusion에서 exact second-order correlations를 얻었다고 해서 field의 connected fourth cumulant까지 0이 되는 것은 아니다. 이번 Wick closure는 non-Gaussian atomic fluorescence나 field fourth cumulant의 크기를 증명하지 않는다.

**2. Optical offset과 photocurrent RF를 분리한다**

\[
d_j(t)=\int d\nu\,d_j[\nu]e^{-i2\pi\nu t},\qquad
[d_j[\nu],d_k^\dagger[\nu']]=\delta_{jk}\delta(\nu-\nu').
\]

ν는 각 optical carrier 주위의 envelope offset [Hz]이고, Ω/2π=f는 검출 current의 RF [Hz]이다. 기존 phase-selected pair에서 유지하는 moments는

\[
\langle d_j^\dagger[\nu]d_j[\nu']\rangle=n_j(\nu)\delta(\nu-\nu'),
\quad
\langle d_p[\nu]d_c[\nu']\rangle=m(\nu)\delta(\nu+\nu').
\]

Normal interbeam와 same-beam anomalous moments는 이 sector closure에 없다. m_cp(ν)=m_pc(−ν)이다. n과 m는 dimensionless이고 Φⱼ=∫nⱼ(ν)dν의 단위는 photons/s다.

`PairedOpticalBins`는 이 함수를 동일 폭의 piecewise-constant bins로 표현한다. 각 두-mode pair에서 두 ordered covariance의 positivity를 검사하고, nonnegative occupations를 요구한다. `paired_optical_bins`는 main/companion symmetry를 검증한 뒤

\[
n_p(\nu)=C^<_{00}(\nu),\quad
n_c(\nu)=C^>_{11}(-\nu),\quad
m(\nu)=C^>_{01}(\nu)
\]

를 사용한다. Occupation을 구하기 위해 큰 vacuum identity를 빼지 않는다. 기존 four-sideband quadrature covariance에서 계산한 occupation과의 parity도 검사한다.

**3. Quadratic current의 spectrum**

Mean product를 뺀 photon-flux covariance의 two-sided Fourier spectrum을 Cⱼₖ(f)라 하자. Wick contraction에서

\[
C_{jj}(f)=\Phi_j+\int d\nu\,n_j(\nu)n_j(\nu+f),
\]
\[
C_{pc}(f)=\int d\nu\,m^*(\nu)m(\nu+f),\qquad
C_{cp}(f)=C_{pc}^*(f).
\]

Φ 항은 spontaneous photon 자체의 shot noise다. Anomalous convolution의 복소 위상은 상대적인 beam delay에 필요하므로 버리지 않는다. 정상상태의 평균 DC delta peak는 위 식에 포함하지 않는다.

검출 loss η를 먼저 nⱼ→ηⱼnⱼ, m→√(ηₚη꜀)m로 적용하고, w=(hₚ,−g h꜀)로 current를 합한다. Positive one-sided PSD는

\[
S_{\rm quad}^{(1)}(f)=2e^2\,w(f)C(f)w^\dagger(f),\qquad
\Delta S_{\rm SQL}^{(1)}(f)=2e^2\sum_j|w_j|^2\Phi_{j,\rm det}.
\]

단위는 A²/Hz이며 추가 2π는 없다. 기존 bright-carrier 결과와 합칠 때 numerator에 quadratic PSD를 더하고 denominator에도 spontaneous mean의 SQL을 더한다. Electronics는 한 번만 더하며 subtraction이나 spectrum fitting은 없다.

\[
S_-^{\rm Gaussian}(f)=
\frac{S_{\rm bright}^{(1)}(f)+S_{\rm quad}^{(1)}(f)}
{S_{\rm SQL,coh}^{(1)}(f)+\Delta S_{\rm SQL}^{(1)}(f)}.
\]

Pair-correlated spontaneous photon은 numerator와 denominator에 다르게 기여하므로 보정이 항상 squeezing을 악화시키는 것은 아니다. 낮은 seed에서 그림의 값이 더 음수가 되는 이유는 이 조건부 covariance와 mean SQL의 조합이다. 이를 실험 squeezing에 맞춘 결과로 사용하지 않는다.

**4. Finite optical band의 의미와 수렴**

각 carrier 주위에 동일한 ideal optical top-hat collection |ν|<B를 선언한다. 그 밖에는 canonical vacuum port가 있다. 따라서 band 밖의 occupation은 없지만 vacuum commutator는 살아 있고, spontaneous shot term Φ는 RF에 대해 white다. Optical filter의 손실 port를 제거한 band-limited commutator로 shot noise를 계산하면 잘못된 결과가 된다.

구현은 zero-padded bins의 겹치는 길이를 정확히 계산한다. RF가 bin 폭의 정수배가 아니면 인접 두 lag의 overlap을 사용한다. FFT의 circular wraparound는 없다. Coherent carrier와 ±RF sidebands가 이 optical band 안에 있는 경우에만 이전 bright readout과 합친다.

RF는 0.1–4 MHz, Gaussian peak pump는 600 mW/530 μm, n=10¹⁸ m⁻³, A=1.2×10⁻⁷ m², L=12.5 mm, Δ/2π=0.9 GHz, δ/2π=−8 MHz, declared reset=2π×100 kHz, η=0.85이다. 아래는 0.5 MHz bin에서 optical collection을 바꾼 대조다.

| Optical half-width B | Probe spontaneous flux [s⁻¹] | Conjugate spontaneous flux [s⁻¹] |
|---|---:|---:|
| 8 MHz | 1.863743×10⁶ | 1.849742×10⁶ |
| 32 MHz | 8.358611×10⁶ | 8.297646×10⁶ |
| 128 MHz | 4.822724×10⁷ | 4.551670×10⁷ |
| 512 MHz | 6.145551×10⁷ | 5.852211×10⁷ |
| 1024 MHz | 6.444965×10⁷ | 6.145760×10⁷ |

B=512→1024 MHz에서도 collected flux가 최대 5.016% 증가한다. 따라서 **unfiltered fluorescence tail의 상한은 아직 없다.** 기존 RF band의 photon occupation 적분을 전체 spontaneous flux로 사용하는 것도 적절하지 않다. 이 optical band는 아직 독립적으로 측정한 실험 filter가 아니다.

반면 같은 B=1024 MHz에서 bin width를 1, 0.5, 0.25, 0.125, 0.0625, 0.03125 MHz로 줄이는 것은 수치 수렴 검사다. 처음 0.5→0.25 MHz 비교는 quadratic PSD의 최대 상대 변화 0.00172155로 사전 기준 10⁻⁴를 실패했다. 마지막 0.0625→0.03125 MHz 변화는 **3.18534×10⁻⁵**로 통과했다. 기준을 완화하지 않았으며 실패 보고서도 보존했다. 최종 65,536 bins와 complex m 값을 v2 JSON에 저장했다.

최종 source flux는 (6.4450819828×10⁷, 6.1458773214×10⁷) photons/s다. 이는 선언한 두 collected bands의 값이며 모든 방향·편광·모드의 fluorescence가 아니다.

**5. Gaussian closure와 별개인 finite-seed mean 검사**

Pump-only generator L₀를 유지한 채, weak propagation이 예측한 z=0,L/2,L의 carrier amplitudes를 각각 고정한다. 같은 reciprocal dipole로

\[
H(t)/\hbar=H_0/\hbar+V e^{-i\omega_b t}+V^\dagger e^{i\omega_b t},
\quad \omega_b=-\omega_{\rm hf}+\delta,
\]
\[
V=g_p\beta_pO_p^\dagger+g_c\beta_c^*O_c^\dagger
\]

를 구성한다. Cp=−i[V,·], Cm=−i[V†,·]이고 ρ(t)=Σₙρₙexp(−inωbt)로 전개하면

\[
(L_0+in\omega_b)\rho_n+C_p\rho_{n-1}+C_m\rho_{n+1}=0.
\]

Floquet orders 2/3/4를 비교하고 96개 phase에서 trace/Hermiticity/positivity를 검사한다. 경계 밖 harmonic residual도 저장한다. Zero-harmonic의 forcing이 seed power에 비례해 작아지는 영역에서는 cancellation floor를 포함한 residual 검사를 병행한다. 기존 forcing-relative residual의 값 자체도 숨기지 않고 기록한다.

독립 dense Floquet matrix와의 baseline 최대 원소 차이는 5.15×10⁻¹⁴ 이하이다. 별도의 trace-bordered 1차 Liouvillian solve로 ρ₁^(weak)를 얻고, 이 polarization이 기존 microscopic field drift M b와 일치하는지 검사했다. Global drive phase를 회전했을 때 ρₙ의 phase가 exp(inθ)로 변하는지도 검사한다.

이것은 prescribed carrier에 대한 **국소 back-action 진단**이다. 세 위치 사이의 오차 상한이나 self-consistent propagation을 증명하지 않는다. 주기적으로 변하는 atomic state를 단순 평균해 기존 stationary D에 넣지 않았으며, finite-seed Floquet quantum diffusion은 아직 연결하지 않았다.

**6. 현재 8 μW 결과와 seed 범위**

최종 B=1024 MHz 조건에서 8 μW seed의 결과는 다음과 같다.

| 진단량 | 값 |
|---|---:|
| 검출 spontaneous/coherent SQL 비 | 3.212414×10⁻⁶ |
| Quadratic/bright PSD의 최대 비 | 9.844130×10⁻⁷ |
| 0.1–4 MHz의 최대 readout 보정 | 1.032969×10⁻⁵ dB |
| 1 MHz 기존 S₋ | −0.7302943383 dB |
| 1 MHz Gaussian 보정 후 S₋ | −0.7303045589 dB |
| 세 z 위치의 최대 평균 원자상태 trace distance | 5.104738×10⁻⁵ |
| 세 z 위치의 최대 polarization 상대 변화 | 9.737180×10⁻⁵ ≈ 0.00974% |

0.01 pW부터 10 mW까지 seed를 바꿨다. `readout 변화 <0.01 dB`, `spontaneous/coherent SQL <1%`, `국소 polarization 변화 <1%`라는 선언된 진단 기준을 사용한다. 이는 coefficient fitting이나 엄밀한 전체 모델 오차 상한이 아니다.

- 1 nW에서는 readout 변화가 약 0.0813 dB라 이 기준을 실패한다.
- 검사한 10 nW, 100 nW, 1 μW, 8 μW, 10 μW, 100 μW에서는 세 기준을 통과한다.
- 1 mW에서는 국소 polarization 변화가 약 1.20%, 10 mW에서는 약 10.54%라 실패한다.

위 값은 샘플링한 점들에 대한 판정이다. 정확한 연속 구간의 임계 seed나 실제 장치의 운용 범위를 인증하지 않는다. Gaussian plotting curve는 원자상태 변화가 큰 seed에서도 비교용으로 표시되며, 그 영역의 corrected spectrum이 finite-seed quantum 예측은 아니다.

**7. 검증·남은 과제**

Thermal top-hat의 triangular bunching과 white shot limit, pure twin-beam의 intensity correlation, equal loss, unequal losses/complex detector weights, vacuum, fractional-bin overlap을 검사했다. 네 optical modes의 truncated Fock state에서 bilinear photocurrent operators를 직접 적용한 결과와 Gaussian formula의 차이는 8.72×10⁻¹⁵ 이하이다. 이 reference는 Wick convolution 코드를 사용하지 않는다.

추가 13개를 포함해 quantum 테스트 83개가 통과했다. 보고서의 source hashes 28개를 확인했다. 전체 suite 결과는 [research_log.md](research_log.md)에 기록한다.

```powershell
python -m pytest -q tests/quantum
python -m analysis.grand_challenge.seed_validity_audit --output NEW.json --plot NEW.png
python -m pytest -q
```

다음 연결점은 periodic finite-seed state와 일관된 drift/diffusion이다. Stationary D에 평균 ρ만 대입하는 방식 대신 Floquet noise와 harmonic correlations를 유지해야 한다. Full Zeeman의 signed CG·실제 편광, non-Gaussian fourth cumulants, 실제 collection/filter inputs와 대역 밖 기여, angular-Doppler, depletion 및 held-out 실험 검증도 남아 있다. Grand Challenge와 milestone 1은 계속 진행 중이다.

**2026-09-09 후속 구현:** [Periodic atomic noise 단계](periodic_noise_derivation.md)에서 주기적 atomic A/D와 harmonic correlations를 구현하고 독립 시간영역 QRT로 검증했다. 같은 8 μW prescribed exit에서 특정 probe atomic-noise 성분은 base 1 MHz에서 약 40% 증가했다. 위의 작은 mean-polarization 및 Gaussian readout 진단을 전체 finite-seed noise 오차 상한으로 사용할 수 없다는 결과다. 본 문서의 기존 S₋ 수치는 그대로 조건부 pump-state/Gaussian 결과이며, 새 atomic noise를 finite-seed field propagation과 검출 readout에 연결하는 작업은 남아 있다.
