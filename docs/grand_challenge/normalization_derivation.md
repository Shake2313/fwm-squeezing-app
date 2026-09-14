# S1 — pump/weak-field dipole과 manifold population 정규화 감사

2026-09-09. 결과: pump와 weak field에 같은 dipole 행렬을 쓰는 조건부 경로를 구현했다. 독립적인 자기부준위 계산으로 기존 공통 1/12의 약한 흡수 한계 불일치를 확인했다. **강한 pump에서 24준위를 4준위로 줄이는 정리는 아직 증명하지 않았다.**

재현 자료: [JSON 보고서](s1_normalization_report.json), [비교 그림](s1_normalization_comparison.png). 이전 readout 보고서는 당시 규약을 기록한 이력으로 보존한다. 이번 변경은 `gabes/fwm_quantum` 연구 경로에만 적용한다.

**1. 전기장과 Rabi의 같은 규약**

원형 Gaussian beam의 전체 power P, intensity 1/e² radius w를 사용한다.

\[
I_0=\frac{2P}{\pi w^2},\qquad E_0=\sqrt{\frac{2I_0}{\epsilon_0c}},
\qquad \Omega_{ge}=\frac{d_{ge}E_0}{\hbar},\qquad
(H/\hbar)_{ge}=\Omega_{ge}/2.
\]

여기서 E(t)=E₀ cos(ωt)인 real peak field다. Weak traveling field는 photon flux amplitude a에 대해 E₀=Qa, Q=√(2ℏω/ε₀cA), g=dQ/(2ℏ)를 쓴다. A=πw²/2인 동일 local intensity에서 `2g√(P/ℏω)=Ω`가 성립하는지 자동 검사한다. 일반적인 Gaussian pump와 collected weak mode의 transverse overlap을 이 등식으로 대신하지 않는다. 이번 cell은 pump 중심의 E₀를 균일하게 샘플링한다.

Steck의 reduced dipole dJ 규약에서는 D1의 Jg=Je=1/2에 대해 Γ=ω₀³dJ²/(3πε₀ℏc³)다. 다른 Wigner–Eckart 규약의 √(2Jg+1)을 추가하면 안 된다. Rabi, saturation intensity, angular strengths의 원전은 [Steck, Eqs. 34–45, 49–50와 Tables 4, 7–8](https://steck.us/alkalidata/rubidium85numbers.pdf)이다.

기존 `rabi_freq`를 역산하면

\[
d_{\rm pump,implicit}=\hbar\Gamma\sqrt{\frac{\epsilon_0c}{4I_{\rm sat}}}
=1.4646827595\times10^{-29}\ {\rm C\,m}.
\]

이는 dJ/√3와 거의 같으며, 기존 quantum fixture의 weak-field base dipole dJ/√12=7.3257088906×10⁻³⁰ C m의 **1.9993734141배**다. 따라서 이전 fixture는 weak-field coupling과 noise 내부에서 commutator를 보존하더라도 pump와 weak field가 같은 절대 dipole을 사용한 것은 아니다. Noise에만 상수를 곱해 해결할 수 없다.

새 RMS 경로는 현재 Lindbladian의 Γ=2π×5.746 MHz와 centroid에서 dJ=2.5367896004×10⁻²⁹ C m를 유도한다. 저장된 dJ=2.5377×10⁻²⁹를 사용해 역산한 Γ와의 차이는 +0.0717886%다. 이는 현재 상수 표의 내부 일관성 차이이며 실험적으로 유의한 불일치라는 주장은 아니다. `species`의 5.7500 MHz 값으로 Γ를 암묵적으로 바꾸지 않는다.

**2. Sublevel sum과 manifold mean의 차이**

`rho[F,F]`가 전체 F manifold population pF일 때, 균일한 m population의 약한 흡수는

\[
W_{FF'}=p_F\,\frac{1}{2F+1}\sum_{m,m'}
|d^{(q)}_{Fm,F'm'}|^2
=p_F d_J^2\frac{S_{FF'}}3.
\]

GABES의 CF² 표와 독립적인 I⊗J dipole reference를 비교하면

\[
3C_F^2=(2F+1)S_{FF'}/3.
\]

왼쪽은 한 polarization의 sublevel **합**이다. 기존 4준위의 pF와 함께 여기에 공통 1/12를 적용하면 목표 mean의 (2F+1)/12배만 남는다.

| F→F′ | S/3: 필요한 평균 d²/dJ² | 기존 3CF²/12 | 기존/필요 |
|---|---:|---:|---:|
| 2→2 | 2/27 | 10/324 | 5/12 |
| 2→3 | 7/27 | 35/324 | 5/12 |
| 3→2 | 5/27 | 35/324 | 7/12 |
| 3→3 | 4/27 | 28/324 | 7/12 |

즉 합산된 `3CF²` 뒤의 scalar는 F=2에서 1/5, F=3에서 1/7이어야 한다. **하나의 공통 scalar로 두 manifold를 동시에 맞출 수 없다.** Thermal pF=(2F+1)/12가 이미 density matrix에 있으면 각 m의 population은 pF/(2F+1)=1/12이다. 전체 pF에 다시 1/12를 곱하는 것은 이 연산과 다르다.

이 결론은 B=0, incoherent uniform sublevels, 약한 absorption 조건에서의 수학적 대조다. Pump가 m population과 coherence를 바꾸는 경우에는 이 평균 자체를 다시 계산해야 한다.

독립 reference `analysis/grand_challenge/reference/uncoupled_d1.py`는 CF²나 6j, 기존 Zeeman builder를 사용하지 않는다. I=5/2와 J=1/2의 닫힌 형태 spinor로 coupled basis를 만들고, 전자 Pauli spherical tensor에 nuclear identity를 곱해 변환한다. 각 σ⁻/π/σ⁺의 평균을 비교하고, 전체 흡수·방출 closure와 basis unitarity를 검사한다. 평균 세기 차이는 최대 1.12×10⁻¹⁶, closure 잔차는 2.23×10⁻¹⁶ 이하이다. 이는 dipole reference이며 full-atom FWM solver는 아니다.

**3. 명시적인 두 reciprocal variant**

- `legacy-reciprocal`: 이전 weak-field 행렬 dJ√(3CF²)/√12를 pump에도 적용한다. Reciprocity 변경의 영향을 분리하는 기준이다.
- `uniform-zeeman-rms`: 자연폭에서 유도한 dJ와 `d_ge=dJ√(S_FF′/3)`를 pump Hamiltonian, weak drive, polarization readout에 함께 사용한다. Uniform-manifold absorption sum rule을 만족한다.

두 번째 variant의 양의 RMS amplitude는 합산된 **세기**를 대표한다. 실제 CG 경로의 부호, loop phase, Zeeman coherence, 서로 다른 pump/probe 편광의 interference를 복원하지 못한다. 따라서 RMS 결과를 full atom에서 유도한 강한 pump FWM 값이라고 해석하지 않는다. 모든 variant에서 pump state부터 M, 두 noise ordering, covariance와 current를 다시 계산한다. Radiative collapse branching과 reset 모델은 동일하다.

**4. Optical carrier의 기준점**

현재 Hamiltonian에서 Δ=0은 F=2→F′=3이다. `OMEGA_D1`은 fine-structure centroid이므로 carrier를 만들 때 다음 offset이 필요하다.

\[
\omega_{23}=\omega_{D1}+\frac{7\omega_{\rm hf,g}+5\omega_{\rm hf,e}}{12},\qquad
\omega_P=\omega_{23}+\Delta,
\]
\[
\omega_{\rm pr}=\omega_P-\omega_{\rm hf,g}+\delta,\qquad
\omega_{\rm co}=\omega_P+\omega_{\rm hf,g}-\delta.
\]

Offset/2π=1.921502256083 GHz다. Atomic Hamiltonian의 Δ를 이만큼 옮기는 수정이 아니다. 이미 그 transition에 대해 정의된 Δ는 유지하고, 광자 에너지·field normalization에 사용하는 optical carrier만 고친다. Carrier만 수정한 대조군에서 probe gain은 1.0517029257→1.0517034724로 변한다. 또한 `species.A_P12×3`와 Hamiltonian excited splitting에는 1 kHz의 rounding 차이가 있어, 이 adapter는 Hamiltonian에 쓰인 값을 일관되게 따른다.

**5. 같은 조건부 입력에서의 비교**

P=600 mW, w=530 μm, n=10¹⁸ m⁻³, weak area=1.2×10⁻⁷ m², L=12.5 mm, Δ/2π=0.9 GHz, δ/2π=−8 MHz, declared reset=2π×100 kHz, seed=8 μW. Detector는 η=0.85, g=1, flat current response, electronics=0이다. 이 값들은 실험의 독립 측정 자료로 승인되지 않았다.

| Variant | Probe power gain | Conjugate power gain | 검출 S₋, 0.1–4 MHz [dB] |
|---|---:|---:|---:|
| 이전 mixed fixture | 1.05170293 | 0.05398004 | −0.365852…−0.365321 |
| Carrier만 수정 | 1.05170347 | 0.05398059 | −0.365856…−0.365325 |
| 같은 historical dipole | 1.02316796 | 0.02568712 | −0.166100…−0.163718 |
| 같은 Zeeman RMS dipole | 1.11065587 | 0.11356221 | −0.730416…−0.728436 |

Independent ordered-Nambu current reference와 PSD의 최대 상대 차이는 3.37×10⁻¹⁵이다. 네 조건의 channel CP, commutator와 covariance 검사가 통과한다. 내부 검증 통과만으로 어느 variant가 실제 −7.8 dB를 설명하는지 결정하지 않는다.

**6. 입력과 증거를 실제 계산 값에 결합**

`ReducedPowerInputs`와 `power_normalized_readout`은 같은 dipole contract에서 pump부터 readout까지 연결한다. `consumed_inputs`는 pump/seed power, waist, density, area, detuning, reset, 길이, phase, atomic frequency/decay와 detector의 주파수별 response/electronics 값을 기록한다. Density와 reset는 직접 입력이며, 온도로부터 암묵적으로 추정하지 않는다. Photon-flux normalization의 carrier와 Rabi, dipole 행렬은 파생 값으로 저장한다.

`audit_consumed_inputs`는 실제 소비한 260개 scalar의 value/unit과 `ParameterEvidence`를 대조한다. 작은 dipole도 무시하지 않도록 단위 변환 후 exact equality를 사용한다. 누락·다른 단위·다른 값·assumed status·target-dataset leakage는 통과하지 않는다. 독립 실험의 beam-profile fit은 허용하지만 target squeezing을 이용한 input 추정은 독립 증거가 아니다. 이 보고서의 입력은 모두 assumed로 기록해 **independent-input audit가 예상대로 실패**한다. Scalar metadata가 모두 갖춰져도 geometry와 reduced-model 가정이 실험적으로 증명되는 것은 아니다.

```powershell
python -m pytest -q tests/quantum
python -m analysis.grand_challenge.normalization_audit --output NEW.json --plot NEW.png
python -m pytest -q
```

Quantum 검사는 추가 11개를 포함해 70개가 통과했다. 전체 검사 기록은 [research_log.md](research_log.md)에 남긴다. CLI는 기존 report/plot을 덮어쓰지 않는다.

다음 구현 대상은 bright/weak-seed approximation의 정량적 오차이다. 대역 밖 spontaneous mean과 quadratic current, finite-seed atomic back-action을 계산해 gain/S₋의 적용 범위를 정해야 한다. Full Zeeman 확장에서는 이번 독립 dipole reference를 사용해 polarization과 signed CG 경로를 복원하고 RMS surrogate와 비교한다. 독립 입력 측정·불확도, Doppler/modes, depletion, held-out 검증은 계속 남아 있다.

후속 구현: [seed_validity_derivation.md](seed_validity_derivation.md)에 finite optical-band Gaussian correction과 국소 finite-seed mean 진단을 연결했다. 현재 8 μW fixture의 보정은 작지만 full tail/비-Gaussian cumulant 상한과 finite-seed quantum diffusion은 아직 검증하지 않았다.
