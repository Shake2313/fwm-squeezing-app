# Pump-only transport ensemble — 공통 lab 주파수와 수렴 증거

2026-09-14. [구현](../../gabes/fwm_quantum/transport_ensemble.py), [검사](../../tests/quantum/test_transport_ensemble.py), [audit](../../analysis/grand_challenge/transport_ensemble_audit.py), [보고서](transport_ensemble_report_v1.json), [plot](transport_ensemble_v1.png).

검증 근거와 함께 제공된 단일 원자의 pump-only finite wavepacket을 열적 경계 flux로 합산하는 adapter다. [열적 유입](inflow_derivation.md)의 양의 경로별 rate와 실제 first-exit residence를 사용하고, [moving reduced Rb packet](rb_transport_derivation.md)의 공통 lab RF·Nambu convention을 명시적인 입력 계약으로 연결한다. 하나의 공통 입사 위상 적분은 pump-only symmetry로 수행하므로 경로마다 위상을 바꾸어 원자를 다시 풀지 않는다.

**구현 대조군은 통과했지만, 최종 toy ensemble은 마지막 refinement 변화가 6.515%로 선언한 5% 기준을 넘어 미수렴·미인증이다.** 단일 path evidence와 위상·source·number 합산의 정확성을 실제 ensemble quadrature 수렴과 구분한다. 실제 열적 Rb polarization ensemble, smooth Gaussian Hamiltonian, nonlocal Maxwell 전파 또는 실험 squeezing이 이 결과로 검증된 것은 아니다. `certified`는 선언한 모델·수치 예산·제공된 비교 증거에 한정된 상태다.

## 단일 원자 packet의 명시적 convention

`PumpOnlyConvention`에는 공통 `AnalysisFrequencyAxis`, number density, reciprocal coupling scales, dimensionless readout scales, ordered reservoir source names, 정지 원자의 carrier offsets, 각 port의 상대 wavevector, model ID·scope·출처가 들어간다. Nambu 순서와 phase charges는 다음으로 고정한다.

\[
\mathbf O=(O_p,O_c,O_p^\dagger,O_c^\dagger)^T,\qquad
\mathbf s=(1,-1,-1,1).
\]

정지 carrier offsets는 하나의 scalar beat \(b\)로부터 정확히 \(\boldsymbol\nu^{(0)}=b(1,-1,-1,1)\)이어야 한다. Dagger pairing만 맞는 임의의 두 정지 carrier를 허용하지 않는다. 상대 wavevectors는 dagger 쌍에서 부호가 반대여야 한다. Coupling과 readout scales는 양수이고 dagger 쌍에서 같아야 한다. Source ledger는 `atomic_inflow`로 시작하며 이름이 중복될 수 없다. 이 첫 source는 한 원자의 boundary state 내부 불확도다. 독립 원자의 유입 수 변동은 뒤에서 별도 `poisson_number`로 계산한다.

Model scope는 `smooth`, `prescribed_segments`, `analytic_fixture` 중 하나를 명시한다. 마지막 두 scope에서 통과한 구간 모델이나 대조군을 smooth 모델의 수렴 증거로 바꾸지 않는다. `model_id`, `provenance`, `phase_provenance`는 구체적인 모델과 출처를 기록해야 한다.

Factory는 `SuppliedPathPacket(packet, convention, evidence)`를 반환한다. 주요 packet schema는 다음과 같다.

| Key | 내용 |
| --- | --- |
| `analysis_axis` | 모든 경로가 공유하는 lab RF axis |
| `frequencies_rad_s` | 각 경로의 Doppler offset을 포함한 `[frequency, 4]` 적분 주파수 |
| `source_names` | Convention과 같은 순서의 boundary·reservoir 이름 |
| `greater`, `lesser` | Connected atomic pulse covariance `[frequency, 4, 4]` |
| `greater_by_source`, `lesser_by_source` | 완전 source 분해 `[source, frequency, 4, 4]` |
| `mean_pulse` | Conditional pulse mean `[frequency, 4]` |
| `retarded_response` | Reciprocal weak-drive response `[frequency, 4, 4]` |
| `residence_time_s` | 공급된 경로의 실제 first-exit residence |
| `metadata` | Mode ordering, coupling, entry phase, carrier offsets와 모델 범위 |
| `audit` | 단일 packet 내부 검사 상태 |

Adapter는 shape, finite values, mode/source 순서, 실제 residence, coupling과 carrier geometry를 검사한다. Connected covariance는 양 ordering과 source별로 Hermitian·PSD여야 하며 각 ordering의 합은 source ledger 합과 일치해야 한다. Source closure의 상대 scale은 각 RF row에서 따로 계산하여, 강한 다른 RF row가 약한 RF row의 누락 source를 숨기지 않도록 한다. Readout scale로 정한 작은 SI roundoff tolerance는 허용 오차이며 covariance에 더하는 noise나 eigenvalue repair가 아니다.

기존 native atomic packet에서 density는 사용하지 않는 입력일 수 있다. 그래서 convention이 density를 명시하고 inflow와 정확히 일치시킨다. Packet metadata에도 density가 있다면 그것도 같아야 한다. Optional readout units·scales가 metadata에 제공된 경우에는 convention과 비교하고, 그렇지 않으면 convention의 선언을 사용한다. 이런 선언과 일치 검사는 호출자의 Hamiltonian이나 readout이 실제로 그 단위를 구현했는지 자동으로 증명하지 않는다.

## 경로별 Doppler offset과 공통 lab RF

Lab RF row를 \(\Omega_\alpha\), 정지 carrier offset을 \(\nu_j^{(0)}\), pump에 대한 상대 wavevector를 \(\mathbf q_j\)라 하자. 속도 \(\mathbf v_\ell\)인 경로의 offset과 적분 주파수는

\[
\nu_{\ell j}=\nu_j^{(0)}-\mathbf q_j\cdot\mathbf v_\ell,
\qquad
w_{\ell,\alpha j}=\Omega_\alpha+\nu_{\ell j}.
\]

코드의 `carrier_offsets_at_rest_rad_s`와 `port_wavevectors_rad_m`가 각각 \(\nu_j^{(0)}\), \(\mathbf q_j\)다. Packet metadata의 offset과 실제 적분 주파수를 이 식과 대조한다. 같은 output row에 모이는 원자들은 모두 같은 **lab \(\Omega_\alpha\)**를 뜻해야 한다. 서로 다른 원자의 generator-frame 주파수를 그대로 같은 row로 합산하지 않는다.

Dagger pairing은 \(\nu_{\ell,p^\dagger}=-\nu_{\ell,p}\), \(\nu_{\ell,c^\dagger}=-\nu_{\ell,c}\)다. 정지 carrier의 공통 beat 관계를 만족해도 **움직이는 원자의** probe와 conjugate offsets가 서로 반대여야 한다는 뜻은 아니다. 비공선 geometry에서는

\[
\nu_{\ell,p}+\nu_{\ell,c}
=-(\mathbf k_p+\mathbf k_c-2\mathbf k_0)\cdot\mathbf v_\ell
\]

가 남을 수 있다. 이 convective loop mismatch를 평균으로 없애거나 scalar compensation으로 닫지 않는다. 실제 entry position에서 정해진 \(-\mathbf q_j\cdot\mathbf r_{\rm in}\) spatial phase도 packet에 유지된다.

Native `smooth`·`prescribed_segments` packet은 metadata의 `entry_optical_demodulation_phases_rad`를 반드시 제공해야 한다. 두 lowering ports의 값이

\[
\begin{pmatrix}\theta_p\\\theta_c\end{pmatrix}
=\begin{pmatrix}\phi\\-\phi\end{pmatrix}
-\begin{pmatrix}\mathbf q_p^T\\\mathbf q_c^T\end{pmatrix}\mathbf r_{\rm in}
\pmod{2\pi}
\]

인지 실제 path entry position으로 검사한다. Wrapped phase 차의 허용 오차는 \(10^{-10}\) rad다. `analytic_fixture`에서는 이 metadata를 생략할 수 있지만, 제공했다면 같은 검사를 적용한다. 이 검사는 entry geometry 선언의 일치를 확인하며 caller의 Hamiltonian 구현 전체를 검증하지는 않는다.

## Pump-only 공통 위상 적분의 정확한 mask

Pump-only Hamiltonian, jumps와 boundary state가 common lab entry phase \(\phi\)에 의존하지 않는다고 가정한다. Readout은 charge \(s_j\), reciprocal input drive는 반대 charge \(-s_j\)를 가진다. 따라서

\[
D(\phi)=\operatorname{diag}(e^{is_j\phi}),\quad
\mathbf m(\phi)=D(\phi)\mathbf m(0),\quad
C^{>,<}(\phi)=D(\phi)C^{>,<}(0)D^\dagger(\phi).
\]

Source별 covariance도 같은 법칙을 따른다. Response의 output row는 \(e^{is_j\phi}\), input drive column은 \(e^{-is_k\phi}\)를 얻으므로

\[
\chi(\phi)=D(\phi)\chi(0)D^\dagger(\phi)
\]

다. Entry phase의 uniform 적분은

\[
\frac1{2\pi}\int_0^{2\pi}e^{i(s_j-s_k)\phi}d\phi
=\delta_{s_j,s_k}.
\]

이를 Hadamard product에 쓰는 mask로 쓰면

\[
M=
\begin{pmatrix}
1&0&0&1\\
0&1&1&0\\
0&1&1&0\\
1&0&0&1
\end{pmatrix},\qquad
\langle X\rangle_\phi=M\odot X
\]

이며 \(X\)는 total/source covariance, response 또는 **이미 계산한 mean outer product**다. Charge가 같은 포트 사이의 상관은 남는다. 각 포트를 독립 위상 평균하면 이 상관을 잘못 없앤다.

모든 charge가 nonzero이므로 \(\langle\mathbf m\rangle_\phi=0\)이지만

\[
\langle\mathbf m\mathbf m^\dagger\rangle_\phi
=M\odot(\mathbf m\mathbf m^\dagger)
\]

는 일반적으로 0이 아니다. Mean을 먼저 평균하고 outer product를 계산하면 Poisson number contribution을 잃는다. Mask는 어느 하나의 고정 entry phase에서 계산한 packet에도 적용할 수 있다. 살아남는 same-charge 성분은 그 공통 위상에 의존하지 않기 때문이다.

이 결과는 **pump-only symmetry를 선언하고 실제 모델에서 확인했다는 조건**에 따른다. Finite seed saturation이나 phase-dependent boundary/Hamiltonian이 들어가면 위 증명을 그대로 쓸 수 없다. `phase_provenance`와 mode charges가 맞는다는 형식 검사만으로 실제 공급된 Hamiltonian의 phase independence가 검증되지는 않는다.

주기적 lab process에서 이 uniform phase average는 deterministic periodic mean을 제외한 zero-cyclic, period-averaged connected spectrum에 대응한다. Coherent periodic mean의 delta lines, 다른 cyclic sectors 또는 full stationarity를 반환·증명하지 않는다.

## 경계 rate를 한 번만 적용하는 atomic stream

경로 \(\ell\)의 경계 rate를 \(J_\ell>0\)라 쓰자. 이미 [thermal inflow](inflow_derivation.md)에서 \(n f(\mathbf v)|\mathbf v\cdot\hat{\mathbf n}|\,dA\,d^3v\)를 적분한 값이며 density와 area가 포함되어 있다. 따라서 별도의 density, velocity Jacobian 또는 residence multiplier를 다시 넣지 않는다.

\[
\begin{aligned}
S_r^{>,<}(\Omega)&=\sum_\ell J_\ell\,M\odot C_{\ell r}^{>,<}(\Omega),\\
S_{\rm number}(\Omega)&=\sum_\ell J_\ell\,M\odot
  [\mathbf m_\ell(\Omega)\mathbf m_\ell^\dagger(\Omega)],\\
S^{>,<}(\Omega)&=\sum_rS_r^{>,<}(\Omega)+S_{\rm number}(\Omega),\\
\chi_{\rm stream}(\Omega)&=\sum_\ell J_\ell\,M\odot\chi_\ell(\Omega).
\end{aligned}
\]

Source 합에서 boundary-state noise를 한 번 포함하고, number contribution을 별도로 한 번 더한다. 두 항은 같은 현상이 아니며 어느 하나를 다른 항에 흡수하지 않는다. Number term은 두 ordering에 동일하므로 ordering 차에 기여하지 않는다.

`readout_units=('1',)*4`인 dimensionless atomic readout에서 single-pulse mean은 s, covariance는 s² 단위다. Weak optical input amplitude의 단위는 \(\sqrt{\text{photon}/s}\)이고 photon count를 dimensionless로 취급하면 reciprocal coupling은 \(s^{-1}/\sqrt{\text{photon}/s}=s^{-1/2}\)다. Packet response는 pulse mean의 input-amplitude 미분이므로 \(s^{3/2}\) 단위이고, rate 합 후에는 다음이 된다.

| Quantity | Single path | Rate-weighted stream |
| --- | --- | --- |
| Connected/raw covariance | s² | s |
| Retarded response | \(s^{3/2}\) | \(\sqrt{s}\) |

Adapter는 upstream packet에 포함된 coupling을 다시 곱하지 않는다. 이런 atomic units는 traveling optical channel, photon-flux covariance 또는 SQL normalization이 아니다. 광학장으로의 실제 비국소 coupling과 검출기 단위 변환은 별도 단계다.

Aggregate `retarded_response`는 각 path의 전체 age 구간에서 수행한 causal double integral을 rate로 가중한 값이다. 이는 선언된 port-drive 패턴에 대한 formal atomic response이며, uniform 또는 canonical optical propagation matrix \(M(\Omega)\)가 아니다. 공간에 따라 달라지는 실제 Maxwell driving을 연결하려면 서로 다른 readout·drive 위치를 보존하는 response kernel \(\mathcal R(\mathbf r,\mathbf r';\Omega)\)가 필요하다. 이미 두 age를 모두 적분한 행렬을 그 spatial kernel 대신 사용하면 안 된다.

## 경로별 numerical evidence 계약

`packet_digest(path, packet, convention)`는 소비할 numerical arrays, port frequencies, 실제 entry/velocity/residence, mode·scale·density·model 선언과 출처를 묶는 SHA-256 identifier다. Mean outer product는 mean array로부터 결정되므로 mean이 바뀌면 digest도 바뀐다. `ConvergenceEvidence`는 이 exact target digest와 model scope에 연결된다.

비교 하나는 `NumericalComparison`으로 기록한다. `errors`와 `tolerances`는 같은 이름의 finite numerical 값이어야 하며, error norm과 normalization scale, zero-signal 처리 방법을 `error_definition`에 명시한다. `reference_id`, `candidate_id`와 서로 다른 `control_values`가 비교 계산을 식별한다. Refined control은 coarse에서 fine으로 증가해야 한다. `independent_reference`·`independent_scramble`에서는 reference ID와 candidate ID가 같을 수 없다. Candidate 자체를 독립 reference로 인용하는 것을 거절하는 검사이며, 서로 다른 ID만으로 계산의 수학적 독립성이 증명되지는 않는다. 단순 `passed=True`나 미리 존재한 packet audit boolean만으로 convergence를 인정하지 않는다.

Per-path evidence가 다루어야 할 metrics는 `greater`, `lesser`, 두 `*_by_source`, `mean_pulse`, `mean_outer`, `retarded_response`다. 작은 mean 자체의 오차와 mean outer product의 오차를 구분한다. 내부 covariance만 통과시키고 Poisson mean이나 response를 생략할 수 없다.

각 path는 refinement와 independent reference를 모두 요구한다.

| Model scope | 요구하는 refinement | 추가 독립 비교 |
| --- | --- | --- |
| `smooth` | 세 해상도의 연속된 두 비교 | Fine candidate의 independent reference |
| `prescribed_segments` | Coarse/fine 비교 하나 | Candidate의 independent reference |
| `analytic_fixture` | Coarse/fine 비교 하나 | Candidate의 independent reference |

Smooth 비교는 첫 비교의 fine control·candidate ID가 다음 비교의 coarse control·reference ID와 이어져야 한다. 각 종류의 마지막 candidate는 실제 packet digest에 연결되어야 하고, 필요한 모든 metrics가 선언된 예산을 통과해야 한다. Native packet의 내부 `audit.passed`도 별도로 true여야 한다.

이 계약은 **전달받은 observed error ledger의 구조·예산·대상 일치**를 검증한다. 참조 계산을 자동으로 실행하거나, caller가 제공한 error가 실제 계산으로 얻어졌는지, independent reference가 구현상 독립인지, provenance 문장이 사실인지까지 증명하지 않는다. Audit artifact와 재현 가능한 reference 계산이 그 추가 근거를 제공해야 한다. Reference와 refinement 결과가 모두 같은 누락 물리를 공유할 가능성도 실험 검증과 별도로 남는다.

현재 별도의 single-path raw QRT 검증은 이 계약이 요구하는 **독립 source별 covariance와 retarded response**를 제공하지 않는다. 따라서 해당 단일 경로 보고서를 adapter의 complete `independent_reference` evidence로 자동 승격할 수 없다. Total raw/connected moments의 독립 일치는 그 측정량의 근거이며, 빠진 source·response error ledger를 대신하지 않는다. 필요한 증거가 없는 frozen Rb packet을 받아도 certification을 거절하는 상태를 유지한다.

## 완전한 stream 이후 actual-integrand 수렴

`aggregate_pump_only_stream(inflow, convention, packet_factory)`는 각 경로에 factory를 정확히 한 번 호출하고 즉시 rate-weighted 항을 누적한다. Native solver가 필요하면 호출자가 factory에서 수행한다. Adapter는 모든 path/phase packet을 보관하지 않으며, phase를 바꾼 재계산이나 별도의 원자 solve를 시작하지 않는다.

잘못된 schema·convention이나 계산 예외는 즉시 거절한다. Packet audit 또는 path evidence가 실패하면 반환값은 `spectra=None`, `certified=False`이며, 실패한 path index·출처·rate와 `consumed_path_count`를 남긴다. 이후 경로를 계속 풀거나 실패한 경로를 빼고 남은 rate를 정규화하지 않는다. 미완성 누적합을 physical ensemble 결과로 제공하지 않는다.

모든 경로가 통과하면 source closure를 유지한 spectra와 `candidate_digest`를 반환하지만 `certified=False`다. 이 상태는 **모든 단일 path의 제공된 증거를 통과한 diagnostic ensemble candidate**다. \(\sum_\ell J_\ell\tau_\ell\simeq nV\)나 Maxwell moments의 수렴만으로 actual wavepacket 적분까지 수렴했다고 간주하지 않는다.

`certify_pump_only_stream(candidate, ensemble_evidence)`는 실제 stream output에 대해 다음 두 종류의 추가 비교를 요구한다.

1. `ensemble_refinement`: 경로 quadrature 해상도를 바꾼 actual-integrand 합의 비교.
2. `independent_scramble`: 독립 quadrature scramble로 다시 계산한 actual-integrand 합의 비교.

둘 다 raw `greater`·`lesser`, source별 두 covariance, `poisson_number`, `retarded_response`의 오차를 포함해야 한다. 이 비교 대상은 thermal occupation controls가 아니라 최종 atomic spectrum·response 자체다. Spectrum digest에는 convention, 전체 inflow, 소비된 packet digest chain, sources, 모든 누적 수치와 rate·density ledger가 포함된다. 누적 이후 outputs나 계산 정체성이 바뀌면 certification을 거절한다.

제공된 evidence가 모두 통과해 `certified=True`가 되어도 선언한 model scope와 수치 예산에서의 조건부 상태다. 실제 열적 Rb에 대해 이 상태를 보고하려면 실제 Rb path packets의 적절한 refinement·독립 reference와 실제 Rb integrand의 ensemble 비교가 필요하다.

## Audit, 결과와 재현

Audit는 exact toy/constant controls를 사용하여 common-phase mask, source·number 분해, 공통 lab frequency, rate의 density counting, 실패 경로 처리와 evidence 계약을 검사한다. Analytic fixture를 통해 adapter가 올바른 수치를 합산하고 잘못된 증거를 거절하는지를 확인하는 것이며, thermal Rb dynamics를 대신 검증하지 않는다.

Fixture는 \(H=0\), jump가 없는 2-level atom과 명시적인 coherent mixed boundary state, identity 성분을 포함하는 네 readout을 사용한다. Carrier offsets와 상대 wavevectors는 모두 0이며, \(\phi\)는 \(D(\phi)\)로 **부여한 공통 phase mark**다. 실제 nonzero Rb beat의 coherence dynamics나 그 lab-time phase를 계산한 것이 아니다. Source에는 `atomic_inflow`만 있으므로 nonzero microscopic jump-source ensemble의 수렴도 이 대조군으로 검증하지 않는다.

모든 port의 적분 주파수가 같은 이 fixture에서

\[
F_\tau(\Omega)=\int_0^\tau e^{i\Omega a}da
=\tau e^{i\Omega\tau/2}\operatorname{sinc}\!\left(\frac{\Omega\tau}{2\pi}\right),
\quad
T_\tau(\Omega)=\int_0^\tau(\tau-u)e^{i\Omega u}du
=\tau^2\frac{e^z-1-z}{z^2},\quad z=i\Omega\tau
\]

이며 \(T_\tau(0)=\tau^2/2\)다. \(\operatorname{sinc}(x)=\sin(\pi x)/(\pi x)\) convention을 사용한다. Packet covariance는 \(|F_\tau|^2C_0\), mean은 \(F_\tau\boldsymbol\mu_0\), response는 \(T_\tau\chi_0\)로 독립 계산할 수 있다. 이 causal-triangle 식을 unequal-port frequencies의 일반식으로 재사용하지 않는다.

각 실제 chord에서 24-point와 48-point Gauss–Legendre age 적분을 다시 계산하고, fine 결과를 위 closed form과 비교하여 모든 path metrics의 observed error를 만든다. 별도로 실제 `transport.wavepacket` ODE의 일곱 metrics를 같은 constant model의 analytic 결과와 비교한다. Production mask를 사용하지 않는 네 literal phase의 직접 행렬합도 독립 reference로 계산한다.

최종 immutable [보고서](transport_ensemble_report_v1.json)는 `implementation_controls_passed=true`, `expected_controls_passed=false`, `final_candidate.certified=false`다. 구현 검사 통과를 전체 ensemble 인증 통과로 바꾸지 않는다.

가정한 box는 \(0.4\times0.2\times0.2\,\mathrm{mm}^3\), \(T=373\) K, \(n=3\times10^{16}\,\mathrm m^{-3}\)이며 mass만 repository의 \(^{85}\mathrm{Rb}\) 값을 사용한다. Lab RF samples는 DC, 10, 30 kHz다. 면당 \(2^4,2^6,2^8,2^{10}\)개 노드와 seeds 11, 211의 모든 grid에서 실제 fixture integrands를 다시 계산했다. 전체 callback 16,320회가 모두 해당 경로 한 번씩만 소비했고, 최종 nominal grid는 6,144 paths다. Adapter의 공통 phase 처리 때문에 발생한 추가 atomic solve는 없다.

Path error budget은 \(10^{-9}\), ensemble budget은 5%로 사전 선언했다. Error는 각 RF row의 Frobenius difference를 reference-row norm으로 나눈 뒤 RF별 최댓값을 사용한다. Zero-signal 근처에서는 명시적인 단위 일치 scale의 \(10^{-12}\)배를 분모의 바닥값으로 사용하며, source blocks는 각 RF에서 하나의 Frobenius norm으로 비교한다. 이 수치는 관측된 비교 오차이며 전체 연속 적분의 엄밀한 상한이나 실험 불확도가 아니다.

| Fine 면당 노드 수 | Total paths | 이전 grid 대비 최대 output 변화 | 독립 scramble 최대 output 차이 | Ensemble 인증 |
| --- | ---: | ---: | ---: | --- |
| \(2^6\) | 384 | 1.7063% | 7.8597% | 실패 |
| \(2^8\) | 1,536 | 7.2242% | 9.7873% | 실패 |
| \(2^{10}\) | 6,144 | 6.5152% | 4.0642% | 실패 |

Final independent-scramble 비교는 5% 안에 들지만 마지막 refinement는 기준을 넘는다. 오차는 비단조이고, 두 조건을 동시에 만족하지 못해 최종 spectrum은 diagnostic candidate로 남는다. 해당 nominal grid의 점유수는 \(1.00110193117\,nV\), 즉 \(nV\) 오차가 약 0.110193%다. 작은 occupation error만으로 covariance·Poisson·response의 실제 적분 수렴을 인증할 수 없음을 이 실행에서도 확인했다.

| 별도 검사 | 저장된 최대 상대오차 또는 결과 |
| --- | ---: |
| 전체 path Gauss 24→48 refinement | \(1.512\times10^{-15}\) 미만 |
| 전체 path Gauss 48 versus closed form | \(1.844\times10^{-15}\) 미만 |
| 직접 네 phase 합 versus production mask | \(3.116\times10^{-15}\) 미만 |
| 실제 constant-atom ODE의 일곱 packet metrics | \(9.754\times10^{-16}\) 미만 |
| Mean outer 삭제에 따른 raw greater 변화 | 29.8626% |
| 평균을 먼저 취한 뒤 outer product를 계산한 number error | 100% |

Actual ODE control은 한 constant-atom chord에서 수행한 별도 일곱 metrics 검사다. Phase control은 모든 grid에서 네 literal phase를 직접 평균한 결과다. 이러한 내부 검사와 source closure의 통과가 위 표의 ensemble 미수렴을 해소하지는 않는다.

Evidence 음성 대조군에서도 audit boolean만 있는 packet, path evidence만 있고 ensemble evidence가 없는 결과, occupation-only evidence가 모두 인증을 받지 못했다. 세 번째 경로에 실패한 mean-outer error를 의도적으로 넣으면 callback 세 번 뒤 중단하고 `spectra=None`을 반환했다. 전체 선언된 rate를 유지하고 부분합을 재정규화하지 않았다.

별도로 실제 native reduced Rb packet을 1 MHz, 한 frozen segment에서 계산하여 mode·coupling·공간 위상·Doppler schema를 확인했다. 내부 audit는 통과했지만 path convergence evidence가 없으므로 첫 경로만 소비한 뒤 `certified=false`, `spectra=None`으로 거절했다. 이 진단은 실제 Rb 열적 ensemble 수렴 검증이 아니다.

실행 전후 source/test/reference SHA-256 63개가 같았고 최종 파일과도 대조했다. 저장된 보고서의 SHA-256은 다음과 같다.

~~~text
b2571c4777456b0a53ce528bfefecbc99c9dc013f475e79719f29948b9d9141c
~~~

Targeted 검사는 신규 66개를 포함해 126개가 통과했다. 전체 `python -m pytest -q`는 **1,185 passed, 1 failed**였다. 유일한 실패는 기존에 누락된 `docs/FWM physics and analytic reconstruction/FWM_physics.tex`를 읽는 `tests/test_docs_consistency.py::test_repository_visibility_wording_is_consistently_public`이다. 이 문서 작업에서 해당 파일이나 검사 조건을 바꾸지 않았다.

신규 artifact 경로로 실행:

~~~powershell
python -m analysis.grand_challenge.transport_ensemble_audit --output docs/grand_challenge/transport_ensemble_report_v1.json --plot docs/grand_challenge/transport_ensemble_v1.png
~~~

v1 outputs가 이미 있으면 report와 plot을 모두 새로운 version 이름으로 바꾼다. 이전 수치나 source manifest를 덮어쓰지 않는다.

Targeted 검사와 repository 통합 검사:

~~~powershell
python -m pytest -q tests/quantum/test_transport_ensemble.py
python -m pytest -q
~~~
