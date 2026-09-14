# Grand Challenge research log

**2026-09-09 — S0 첫 구현**

사용자가 청사진의 첫 구현을 요청했다. 현재 branch에서 `gabes/quantum`의 frequency/input contracts, explicit reservoirs/CCP audit, finite temporal-mode channel composition을 추가했다. 기존 reduced model의 nonzero coherence-only dephasing은 자동 변환할 수 없는 입력으로 처리한다.

검증한 내용:

- Explicit radiative channels는 기존 `double_lambda_rb85(gamma_gg=0)`의 dissipator를 정확히 재현한다.
- Thermal replacement jump construction은 기존 reset superoperator 및 arbitrary complex operator에 대한 직접 계산과 일치한다.
- Default four-level generator의 complete-positivity 반례와 121 °C combined generator의 통과를 별도로 재현한다. 후자의 통과만으로 개별 reservoir provenance를 승인하지 않는다.
- Vacuum loss, ideal two-mode squeezing과 collection의 covariance/commutator/uncertainty를 검증했다. Gain/loss 합성 순서에 따른 noise 차이도 확인한다.
- 독립 계측 fitting과 목표 FWM 데이터 fitting을 dataset provenance로 구분한다.

재현 명령과 규약: `conventions.md`. 이번 snapshot: `s0_foundation_report.json`. 이전 설계 당시의 `foundation_probe.json`은 원래 기록을 유지한다.

검사 결과: `python -m pytest -q tests/quantum`은 23 passed. `python -m pytest -q`는 674 passed, 1 failed (154.08 s). 실패는 이번 작업 전에도 존재한 `tests/test_docs_consistency.py::test_repository_visibility_wording_is_consistently_public`의 `FWM_physics.tex` 누락이다. 새 quantum tests는 모두 통과했고 실패한 기존 검사를 skip/xfail로 숨기거나 삭제된 문서를 임의 복원하지 않았다.

현재 범위: S0 foundations의 첫 패키지. 전체 operating domain의 reservoir 감사, mean-field result contract와 coupled uncertainty ledger까지 완료한 것은 아니다. Microscopic D, RF field response, physical S_minus 및 experimental validation은 아직 구현하지 않았다.

**S0 종료 당시 다음 미해결 작업**

1. S1의 최소 Hamiltonian/reference state를 standard minus branch의 기존 pump-only reference와 연결하고 basis/정규화 parity를 고정한다.
2. 명시한 radiative/reset model에서 complete atomic operator basis의 ordered Einstein diffusion을 유도한다. 같은 source를 사용하되 D 구현을 공유하지 않는 direct two-time/QRT reference와 비교한다.
3. Pump-only weak-seed reference의 보존 subspace와 선택된 reservoir를 포함한 stationary-state 유일성·stability를 확인한다.
4. Atomic noise를 photon-flux-normalized field noise로 변환할 때 density/volume/mode-area 인자를 유도하고, 그 뒤 finite-RF M/D와 전파를 연결한다.
5. 실제 collision/dephasing 설정은 출처가 명시된 GKSL channels와 coherence 변화의 검증이 필요하다. Legacy fitting coefficients를 새 모델에 옮겨 해결하지 않는다.

체크리스트의 Grand Challenge 및 첫 milestone을 진행 중으로 표시한다. 나머지 실험·full-atom milestone의 완료 조건은 유지한다.

**2026-09-09 — S1 atomic diffusion과 독립 QRT reference**

사용자의 다음 단계 요청으로 S0 기록의 1–3번에 해당하는 single-atom foundation을 구현했다. `gabes/quantum/diffusion.py`는 complete traceless Hermitian basis에서 atomic A를 구하고 explicit jump commutators로 reservoir별 ordered D를 독립 계산한다. C와 D의 positivity, stationary Lyapunov identity 및 expected atomic commutator balance를 자동 검사한다. 이 commutator 검사는 상태에 의존하는 atomic algebra에 대한 것이며, 아직 field input/output preservation을 검증한 것은 아니다.

독립 reference는 `analysis/grand_challenge/reference/atomic_qrt.py`다. 하나는 full Liouville adjoint의 Einstein product rule, 다른 하나는 trace-bordered density-operator QRT 적분이다. Production diffusion이나 atomic drift/resolvent 구현을 공유하지 않는다. Physical generator/state/operators는 공유하므로 이 일치는 구현 검증이며 새로운 실험적 증거가 아니다.

Reduced adapter는 radiation-only four-level atom에 명시한 optional thermal-reset jumps만 더한다. 기존 pump-only reference에도 같은 reservoir를 넣어 stationary-state parity를 확인했다. Lab RF Ω와 static-pump generator frequency를 다른 type으로 유지한다. Pump-off/reset, pumped/radiation-only, pumped/reset 세 경우를 generator DC 부근과 minus-sector frequencies에서 비교했다. 비감쇠 atomic modes는 continuous stationary PSD로 오인하지 않도록 거부한다.

근거 문서: [derivation.md](derivation.md). Snapshot: [s1_atomic_noise_report.json](s1_atomic_noise_report.json). 전체 state, operator basis, drift, reservoir별 diffusion, ordered/symmetrized spectra, frequency axes, source hashes와 signed diagnostics를 보관했다. S0 보고서는 당시 snapshot으로 유지하며 source hash를 현재 버전에 맞추려고 덮어쓰지 않았다.

관측한 수치:

- Jump/adjoint Einstein relative error의 최댓값은 3.54×10⁻¹⁴.
- 두 주파수 축·세 조건에서 spectrum/QRT matrix relative error의 최댓값은 6.66×10⁻¹².
- Stationary Lyapunov relative residual은 최대 1.41×10⁻¹³, expected commutator residual은 최대 1.14×10⁻¹³.
- 기존 pump reference와 stationary-state Frobenius 차이는 최대 3.68×10⁻¹⁴.
- 모든 선언된 spectral PSD/Hermiticity gates가 통과했다. Semidefinite 경우의 작은 음의 고유값은 clipping하지 않고 저장했다.

검사 결과: `python -m pytest -q tests/quantum`은 37 passed. 이번에 추가한 14개에는 analytic two-level Lorentzian의 sign/width/normalization 및 integral, driven QRT, four-level parity, nonstationary/nondecaying 거부와 report 보존 검사가 포함된다. 필수 전체 실행 `python -m pytest -q`는 **688 passed, 1 failed (152.54 s)**. 실패는 S0에서도 확인한 `tests/test_docs_consistency.py::test_repository_visibility_wording_is_consistently_public`이며, 작업 전 삭제 상태였던 `docs/FWM physics and analytic reconstruction/FWM_physics.tex`가 없어서 발생했다. 해당 기존 검사와 문서 삭제 상태는 이번 작업에서 변경하지 않았다.

재현 명령:

```powershell
python -m pytest -q tests/quantum
python -m analysis.grand_challenge.atomic_noise_audit --output NEW_ATOMIC_REPORT.json
python -m pytest -q
```

현재 범위는 static classical pump, single atom/velocity, declared radiative/reset Markov model의 second-order spectrum이다. Nominal pump와 reset은 조건부 reference inputs다. No-fit apparatus-input ledger, collision law, Doppler ensemble, field-noise normalization, propagation 및 measured S₋ 예측은 완료되지 않았다. Grand Challenge와 reduced four-level milestone은 계속 진행 중이다.

**Atomic 단계 종료 당시 다음 미해결 작업**

1. Explicit reservoir model에서 weak-field driving/readout matrices를 유도하고 기존 photon-flux response와 basis·주파수·정규화 parity를 검증한다.
2. Density/coarse-graining volume/mode-area factors를 도출하여 atomic diffusion을 collective polarization과 field D(Ω)로 옮긴다. 독립 velocity classes의 covariance weight를 함께 고정한다.
3. ±Ω companion blocks와 local field commutator를 검증하고, 그 이후 segmentwise M/D propagation 및 finite-temporal-mode channel로 연결한다.
4. Collision/dephasing reservoir의 physical provenance와 모든 consumed inputs의 measurement ledger를 마련한다. 이후 mean fields, detected SQL 및 uncertainty를 연결하여 gain과 S₋를 출력한다.

**2026-09-09 — S1 conditional field M/D와 prescribed-segment propagation**

사용자의 후속 요청으로 atomic→field 연결을 구현했다. `gabes/quantum/traveling.py`는 동일 weak-field Hamiltonian에서 density-operator driving B와 radiation C₀를 만들고, 앞 단계의 atomic jump diffusion으로 greater/lesser field noise를 각각 계산한다. Slice-average noise `D_atom/(n A dz)`에서 density/area normalization을 도출했다. Independent classes는 각 class density가 적용된 covariance를 합하는 인터페이스를 제공한다. 아직 Doppler quadrature나 angular geometry를 생성·검증한 것은 아니다.

`gabes/fwm_quantum/field.py`는 기존 reduced transition readouts와 main/companion generator frequencies를 연결한다. Field drift와 noise의 local commutator identity 및 propagated identity를 검사한다. 공통 uniform transverse area와 explicit effective dipole을 요구하며, 서로 다른 Gaussian waists나 phenomenological efficiency를 입력으로 받아 보정하는 경로는 만들지 않았다.

중요한 정규화 결정: 기존 structural 1/12를 재현하려면 `d_eff=d/sqrt(12)`를 weak-field driving과 emission 양쪽에 적용해야 한다. M은 그대로 두고 noise에만 1/12를 추가하면 commutator가 실패한다. 이 negative control을 보존했다. 반면 기존 pump-power→Rabi와 weak-field dipole의 완전한 microscopic parity는 아직 검증되지 않았으므로 pump Rabi 및 reduced dipole을 conditional inputs로 표시한다. 0.74 residual이나 fitted loss는 사용하지 않았다.

Propagation은 constant segment의 matrix exponential과 distributed covariance integral을 계산한다. 여러 segment의 앞쪽 noise를 뒤쪽 transfer로 전달한 뒤 뒤쪽 source noise를 더한다. Segment의 pump/state/density는 prescribed inputs다. 이를 self-consistent pump depletion으로 표시하지 않는다.

유도와 적용 범위: [field_derivation.md](field_derivation.md). Snapshot: [s1_field_noise_report.json](s1_field_noise_report.json). Report는 −4…4 MHz의 81 RF points, 두 generator axes, 주파수별 M과 reservoir별 두 noise orderings, propagated added noise와 vacuum-output spectral covariances, signed diagnostics, code hashes와 NumPy/SciPy versions를 보관한다. 이전 S0와 atomic report는 원본 snapshot으로 유지했다.

세 조건(pump off/reset, pumped/radiative-only, pumped/reset)의 관측 결과:

- 기존 동일-reservoir pump reference 및 Maxwell photon-flux transfer와의 최대 상대 차이: 2.75×10⁻¹⁴.
- 독립 full-Liouville QRT projection과 greater/lesser field noise의 최대 상대 차이: 5.58×10⁻¹¹.
- Local field commutator relative residual: 최대 1.09×10⁻¹³. Propagated commutator residual: 최대 2.85×10⁻¹⁴.
- Companion symmetry 상대 차이: 최대 3.81×10⁻¹³. Constant medium을 두 구간으로 나눈 propagation 상대 차이: 최대 5.60×10⁻¹⁶.
- 모든 local/global PSD·Hermiticity gates 통과. 잘못된 noise-only coefficient는 세 조건 모두 거부되었다.

검증 과정에서 passive vacuum의 정확히 0인 ordering을 projection roundoff(~10⁻³³) 자체로 나누면 Hermiticity relative check가 거짓 실패하는 경우를 확인했다. 해당 ordering norm에 비례하는 tolerance에 paired reservoir scale의 `64 epsilon` allowance를 더했다. 실제 matrices나 고유값은 수정하지 않았으며, tolerance와 signed results를 문서화했다. Analytic ground-state absorber와 inverted-reset amplifier가 각각 필요한 vacuum noise 및 `G−1` added noise를 재현한다.

검사 결과: `python -m pytest -q tests/quantum`은 **47 passed**. 이번에 추가한 10개 검사는 analytic attenuation/amplification, Maxwell/QRT parity, area/density/dipole scaling, class splitting, companion, numerical integration와 semigroup, noncommuting segment order, invalid coefficient/coordinate 거부 및 immutable report를 포함한다. 필수 전체 실행 `python -m pytest -q`는 **698 passed, 1 failed (229.41 s)**. 실패는 이전과 동일한 `tests/test_docs_consistency.py::test_repository_visibility_wording_is_consistently_public`의 기존 삭제된 `FWM_physics.tex` 누락이다. 해당 테스트와 삭제 상태는 변경하지 않았다.

새 field propagation의 matrix exponential은 SciPy를 사용한다. Optional dependency set은 `requirements-quantum.txt`에 선언했다. 기존 환경에는 이미 SciPy 1.17.0이 있어 추가 설치 없이 검증했다.

```powershell
python -m pytest -q tests/quantum
python -m analysis.grand_challenge.field_noise_audit --output NEW_FIELD_REPORT.json
python -m pytest -q
```

**Field 단계 종료 당시 다음 미해결 작업**

1. Normalized temporal-mode/±RF quadrature 변환, quantum uncertainty 및 mean-field-dependent intensity-difference/SQL readout. 현재 출력된 Nambu spectral covariances를 실험 S₋로 부르지 않는다.
2. Pump-power/Rabi, weak-field dipole, optical carrier, density/reset/transport의 physical provenance와 consumed-input ledger 통합. Conditional numerical fixtures를 독립 계측으로 승격하지 않는다.
3. Angular-Doppler/velocity quadrature, transverse mode overlap과 독립 수렴. Independent-class sum이 collision-induced class correlations까지 표현하지는 않는다.
4. Generated-field back-action과 self-consistent segmentwise pump depletion, 이후 full-atom 및 held-out experimental validation.

Grand Challenge와 reduced-model milestone은 계속 진행 중이다. 조건부 field M/D와 전파는 구현했지만 네 완료 조건의 전체 달성을 뜻하지 않는다.

**2026-09-09 — S1 conditional gain, four-sideband와 intensity-difference/SQL**

사용자의 다음 작업 요청으로 `gabes/quantum/sidebands.py`, `gabes/quantum/readout.py`, `gabes/fwm_quantum/readout.py`를 추가했다. 같은 atomic model에서 coherent carrier transfer와 finite-RF noise를 얻어 gain과 conditional linearized S₋(f)를 함께 출력한다. Four physical optical sidebands를 real Gaussian channel로 변환하고 channel CP와 covariance quantum uncertainty를 검사한다. DC는 명시적 RF=0 carrier row에만 사용한다.

Finite spectral band는 `TopHatBand`의 normalized flat filter와 Gauss quadrature로 정의한다. Spectrum 안에서 X가 변할 때 unobserved orthogonal input modes에서 들어오는 vacuum term을 포함한다. 단순히 평균 X와 평균 Y만 쓰면 phase-varying passive fixture의 CP가 깨지는 negative control을 검증했다. Arbitrary complex filters나 strictly finite-time windows를 구현한 것은 아니다.

Detector model은 두 η, 선언한 balancing g, frequency-resolved complex current response와 additive electronics PSD를 받는다. Means와 sideband covariance에 같은 loss를 적용하고, 같은 detected currents·g·response로 SQL을 계산한다. 출력은 positive-frequency one-sided A²/Hz, quantum/total SQL ratios와 dB다. Additive electronics는 별도로 더하며 암묵적 subtraction이나 squeezing fitting은 없다.

독립 검증 `analysis/grand_challenge/reference/direct_readout.py`는 quadrature conversion을 사용하지 않고 ordered Nambu covariances at ±RF에 current weights를 직접 적용한다. Coherent SQL, ideal bright-seed `1/(2G−1)`, equal loss `1−eta+eta*S_source`, unequal losses/balance와 complex detector phase를 해석적 기준과 비교했다.

유도: [readout_derivation.md](readout_derivation.md). 최종 보고서: [s1_readout_report_v2.json](s1_readout_report_v2.json). 확인한 그림: [s1_readout_spectrum_v2.png](s1_readout_spectrum_v2.png). 초기 그림에서 shared plot style이 curve와 legend의 색을 다르게 만드는 것을 발견해 명시적 palette로 수정했다. 초기 report/plot은 그대로 두고 source hashes가 맞는 v2를 최종본으로 저장했다. 데이터 계산식은 이 그림 수정에서 바뀌지 않았다.

기존 field fixture의 n=10¹⁸ m⁻³, A=1.2×10⁻⁷ m², ℓ=12.5 mm, Δ/2π=0.9 GHz, δ/2π=−8 MHz, Δk=0를 유지했다. Seed=8 μW, source η=1 및 declared detector η=0.85, g=1, h=1, electronics=0를 사용한다. Pump-on은 기존 600 mW/530 μm helper의 Rabi이며 physical normalization은 여전히 conditional이다.

Pumped/reset baseline 결과:

- Coherent probe power gain 1.05170293, conjugate power gain 0.05398004. Photon-flux gain은 별도 값으로 저장했다.
- 0.1–4 MHz source S₋는 −0.433740…−0.433105 dB. Declared 85% detector에서는 −0.365852…−0.365321 dB.
- 같은 current의 one-sided source SQL은 약 1.817408×10⁻²⁴ A²/Hz, detector SQL은 약 1.544797×10⁻²⁴ A²/Hz.
- 세 조건 전체에서 direct Nambu PSD와의 최대 상대 차이 3.35×10⁻¹⁵, equal-loss ratio identity 차이 최대 4.45×10⁻¹⁶.
- 0.8–1.2 MHz normalized band의 orders 4/8/16이 통과했다. Direct weighted PSD와의 상대 차이는 최대 2.59×10⁻¹⁶, order 8→16 ratio 변화는 1.12×10⁻¹⁶.

보고서는 sampled-band spontaneous flux도 부분 진단으로 보관한다. 이는 전체 fluorescence나 quadratic photocurrent noise의 bound가 아니다. `delta-a-dagger delta-a`를 생략한 bright-carrier approximation, weak atomic back-action 및 spontaneous mean contribution의 전체 검증은 별도 작업이다. 실제 −7.8 dB를 설명하거나 실험 uncertainty 안의 절대 예측을 달성했다고 표시하지 않는다.

검사 결과: `python -m pytest -q tests/quantum`은 **59 passed**. 이번 추가 12개는 sidebands/CP/uncertainty, analytic detector/readout limits, independent contraction, normalized filter/leakage와 artifact preservation을 검증한다. 필수 전체 실행 `python -m pytest -q`는 **710 passed, 1 failed (148.81 s)**. 실패는 앞 단계와 동일한 `tests/test_docs_consistency.py::test_repository_visibility_wording_is_consistently_public`의 기존 `FWM_physics.tex` 누락이다. 기존 문서 삭제나 해당 테스트는 변경하지 않았다.

```powershell
python -m pytest -q tests/quantum
python -m analysis.grand_challenge.readout_audit --output NEW_READOUT.json --plot NEW_READOUT.png
python -m pytest -q
```

**현재 다음 미해결 작업**

1. Pump-power/Rabi, weak-field dipole, optical carrier·transition conventions 및 density/reset를 하나의 physical input ledger로 감사한다. Conditional fixture와 independent measurement의 구분을 유지한다.
2. Bright/weak-seed validity와 대역 밖 fluorescence, spontaneous mean, quadratic current correction의 크기를 평가한다. 이후 필요한 finite-seed back-action을 연결한다.
3. Angular-Doppler/velocity와 transverse-mode integration 및 convergence, self-consistent segmentwise depletion과 full atom.
4. 실제 detection/SQL/RBW/VBW inputs와 uncertainty를 적용하고 held-out gain/spectrum을 검증한다.

Grand Challenge 및 milestone 1은 계속 진행 중이다. Conditional gain/S₋ 계산 경로가 열렸으며 experimental no-fit 완료 상태로 바꾸지는 않았다.

## 2026-09-09 — S1 pump/weak dipole과 manifold 평균 감사

Power→Rabi와 weak-field dipole을 같은 규약으로 연결하는 다음 단계를 구현했다. [유도와 적용 범위](normalization_derivation.md), [기계 판독 보고서](s1_normalization_report.json), [확인한 비교 그림](s1_normalization_comparison.png)에 결과를 보존했다.

핵심 발견은 두 가지다. 기존 saturation-intensity pump helper의 implicit dipole은 이전 quantum weak-field base dipole의 **1.9993734141배**다. 또 CF²의 `3CF²`가 한 polarization에 대한 sublevel sum임을 독립 I⊗J dipole reference로 확인했다. 이를 전체 manifold population과 결합할 때 평균은 F=2에서 1/5, F=3에서 1/7이므로 공통 1/12는 두 weak-absorption weights를 동시에 맞출 수 없다. 이전 structural-convention 해석의 한계를 field/readout 유도 문서에 추가했다. 이 결론을 강한 pump의 full-atom reduction 정리로 확대하지 않는다.

`gabes/fwm_quantum/normalization.py`는 Gaussian peak field, 자연폭에서의 Steck-convention dipole, 고정된 두 reciprocal variant와 F2→F′3 carrier anchor를 제공한다. `legacy-reciprocal`은 기존 weak dipole을 pump에 같이 적용한다. `uniform-zeeman-rms`는 `d_ge=dJ sqrt(S_FFprime/3)`를 양의 4준위 amplitude로 사용해 unpolarized weak-absorption mean을 재현한다. 후자의 nonlinear loop phases/Zeeman optical pumping은 아직 보존되지 않는다. 기존 model/field/readout에 명시적 transition-scale 입력을 추가했으며 기본 호출은 이전 결과를 유지한다.

Optical carrier 기준은 centroid보다 1.921502256083 GHz 높은 F2→F′3로 연결했다. Hamiltonian의 Δ는 유지하며 optical photon-energy normalization만 고친다. `species`와 Hamiltonian의 excited splitting 차이 1 kHz와 stored dipole/decay의 +0.0717886% 내부 차이도 숨기지 않고 기록했다. Natural decay는 기존 2π×5.746 MHz를 사용한다.

`gabes/fwm_quantum/inputs.py`의 `ReducedPowerInputs`와 `power_normalized_readout`은 같은 dipole map에서 pump state, M, D, coherent gain과 S₋를 함께 계산한다. 실제 소비하는 scalar 260개의 SI 값과 단위를 evidence에 결합하고, 작은 dipole 값의 오차까지 잡는다. Fixture evidence는 assumed이며 독립 입력 gate가 예상대로 실패한다. Density, reset와 geometry도 아직 독립 측정치가 아니다.

동일 600 mW/530 μm, n=10¹⁸ m⁻³, A=1.2×10⁻⁷ m², 12.5 mm, Δ/2π=0.9 GHz, δ/2π=−8 MHz, seed=8 μW, η=0.85의 비교:

| Variant | Probe power gain | 검출 S₋, 0.1–4 MHz [dB] |
|---|---:|---:|
| 이전 mixed fixture | 1.05170293 | −0.365852…−0.365321 |
| Optical carrier만 수정 | 1.05170347 | −0.365856…−0.365325 |
| 같은 historical dipole | 1.02316796 | −0.166100…−0.163718 |
| 같은 Zeeman RMS dipole | 1.11065587 | −0.730416…−0.728436 |

이는 normalization sensitivity이며 실제 −7.8 dB의 설명이나 fitted improvement가 아니다. Variants마다 state부터 전부 재계산했고 보상 coefficient는 없다. Dipole reference strength 차이 최대 1.12×10⁻¹⁶, ordered-Nambu PSD reference 상대 차이 최대 3.37×10⁻¹⁵, propagated commutator 상대 잔차 최대 2.85×10⁻¹⁴다. Source covariance uncertainty 최소 고유값은 네 조건 전체에서 5.89×10⁻⁴ 이상이다. 보고서의 source hashes 26개가 현재 코드와 일치함을 확인했다.

새 검사 11개는 dipole/linewidth 규약, Gaussian intensity, σ⁻/π/σ⁺의 독립 평균, population-counting 반례, photon-flux reciprocity, Casimir carrier anchor, evidence의 실제 value/unit 및 held-out leakage, no-atom SQL와 immutable artifact를 검증한다. `python -m pytest -q tests/quantum`: **70 passed**. 필수 `python -m pytest -q`: **721 passed, 1 failed (255.76 s)**. 실패는 이전과 같은 `tests/test_docs_consistency.py::test_repository_visibility_wording_is_consistently_public`의 기존 `FWM_physics.tex` 누락이다. 해당 문서 삭제나 테스트를 바꾸지 않았다.

**현재 다음 미해결 작업**

1. Bright/weak-seed approximation의 오차를 계산한다. 대역 밖 spontaneous mean, quadratic photocurrent correction 및 finite-seed atomic back-action을 평가한다.
2. 독립 dipole reference를 full Zeeman Hamiltonian/회전된 실제 편광/collapse로 확장해 RMS surrogate의 강한 pump 한계를 검증한다.
3. 독립 실측 입력·불확도, angular-Doppler와 transverse convergence, self-consistent depletion, held-out gain/S₋ 검증을 연결한다.

체크리스트의 `normalization_implementation`에 발견·구현·검사를 추가했다. Grand Challenge와 milestone 1은 진행 중이며 experimental validation 상태는 false다.

## 2026-09-09 — S1 quadratic photocurrent와 finite-seed 진단

이번 단계는 bright-carrier readout에서 생략한 spontaneous mean/SQL과 Gaussian quadratic photocurrent를 구현하고, 별도 local Floquet mean으로 finite-seed back-action을 평가했다. [유도](seed_validity_derivation.md), [최종 보고서](s1_seed_validity_report_v2.json), [확인한 그림](s1_seed_validity_v2.png)에 수치와 한계를 기록했다.

`gabes/quantum/photocounting.py`는 optical offset bins의 n_p,n_c,m_pc를 ordered Nambu covariances에서 추출하고, 두-beam physical moment 조건과 main/companion reflection을 검사한다. Gaussian Wick factorization으로 normal/anomalous convolutions와 spontaneous shot noise를 계산한다. 두 detector losses, complex response와 balancing을 적용하고 기존 bright PSD 및 coherent SQL에 각각 quadratic PSD와 spontaneous SQL을 더한다. Fractional RF lag는 zero-padded piecewise-constant bins의 겹치는 길이로 계산한다. Optical band 밖의 vacuum port를 유지하며, coherent beat sidebands가 band 안에 있는지 확인한다.

이는 exact atomic second-order theorem을 Gaussian fourth moments로 확장하는 추가 가정이다. Non-Gaussian connected fourth cumulants의 크기나 unfiltered fluorescence의 상한을 얻은 것은 아니다. 이전 bright-carrier API의 계산식은 유지한다.

정규화 단계의 같은 uniform-Zeeman-RMS fixture로 B=8/32/128/512/1024 MHz optical half-width를 비교했다. B=512→1024 MHz에서도 collected spontaneous flux가 최대 5.016% 증가한다. 실제 optical collection/filter 입력이 필요하며, RF=0.1–4 MHz 적분만을 전체 photon flux로 해석하지 않는다.

처음 B=1024 MHz에서 bin width=0.5→0.25 MHz의 최대 quadratic PSD 상대 변화가 0.00172155로 사전 수렴 기준 10⁻⁴를 실패했다. [초기 보고서](s1_seed_validity_report.json)와 그림을 그대로 보존하고, 0.125/0.0625/0.03125 MHz까지 격자를 세분화해 다시 계산했다. 최종 adjacent-grid 변화는 3.18534×10⁻⁵로 통과한다. 기준을 완화하지 않았다. 최종 65,536 bins의 n,m와 source hashes 28개를 JSON에 보존하고 해시 일치를 확인했다.

`gabes/fwm_quantum/seed_validity.py`는 기존 pump-only generator에 reciprocal seed/conjugate Fourier Hamiltonian을 추가해 periodic mean을 계산한다. Weak propagation에서 예측한 carrier를 z=0,L/2,L에서 고정하고 Floquet orders 2/3/4를 비교한다. Density matrix의 trace/Hermiticity/96-phase positivity, ladder 밖 residual을 기록했다. Tiny-seed forcing이 작아질 때 rounding을 과대평가하지 않도록 cancellation-aware residual을 사용하며 기존 relative residual도 보존한다. 별도 trace-bordered first-order response가 microscopic field drift M b와 일치하고, dense extended Floquet matrix 및 harmonic phase rotation과의 parity를 검사했다.

현재 8 μW, B=1024 MHz에서:

- Source spontaneous flux는 probe 6.4450819828×10⁷, conjugate 6.1458773214×10⁷ photons/s.
- Detector spontaneous/coherent SQL 비는 3.212414×10⁻⁶, quadratic/bright PSD의 최대 비는 9.844130×10⁻⁷.
- 0.1–4 MHz의 Gaussian readout 보정은 최대 1.032969×10⁻⁵ dB. 1 MHz는 −0.7302943383→−0.7303045589 dB.
- 세 z 위치의 최대 mean-state trace distance 5.104738×10⁻⁵, polarization 상대 변화 9.737180×10⁻⁵(약 0.00974%).
- Dense Floquet reference와 최대 원소 차이 5.15×10⁻¹⁴ 이하. 독립 four-mode Fock current reference와 Gaussian fourth-moment covariance 차이 8.72×10⁻¹⁵ 이하.

Seed sweep은 0.01 pW–10 mW를 다룬다. Readout 변화<0.01 dB, spontaneous/coherent SQL<1%, 국소 polarization 변화<1%를 선언한 진단 기준으로 사용했다. 검사한 10 nW/100 nW/1 μW/8 μW/10 μW/100 μW는 통과하고, 1 nW는 readout 보정 때문에, 1 mW 이상은 polarization 변화 때문에 실패한다. 연속 구간의 임계값이나 실제 장치 운용 범위를 인증하지 않는다. Mean Floquet back-action을 새로운 stationary diffusion으로 대체하거나 quantum spectrum에 적용하지 않았다.

추가 13개 검사는 thermal triangular bunching/white shot, pure paired-state/loss, fractional overlap과 complex phase, vacuum, unequal detector/SQL/electronics, independent Fock current, sideband occupation parity, dense Floquet/order/weak-power/phase/field-drift limits와 artifact preservation을 다룬다. `python -m pytest -q tests/quantum`: **83 passed**. 필수 `python -m pytest -q`: **734 passed, 1 failed (150.44 s)**. 실패는 이전과 같은 `tests/test_docs_consistency.py::test_repository_visibility_wording_is_consistently_public`의 기존 `FWM_physics.tex` 누락이다. 해당 삭제나 테스트는 변경하지 않았다.

**현재 다음 미해결 작업**

1. Finite-seed periodic state에서 quantum drift/diffusion과 harmonic noise correlations를 일관되게 구성한다. 평균 ρ만 stationary D에 넣는 우회는 하지 않는다.
2. Non-Gaussian fourth cumulant와 band/out-of-sector fluorescence의 크기를 실제 collection/filter inputs와 함께 평가한다.
3. Full Zeeman signed CG·편광, angular-Doppler/modes, independent input uncertainties, self-consistent depletion 및 held-out gain/S₋ 검증을 연결한다.

체크리스트의 `seed_validity_implementation`에 기록했다. Grand Challenge/milestone 1은 진행 중이며 absolute hot-vapor prediction, unfiltered tail bound, finite-seed quantum noise 및 experimental validation은 아직 완료가 아니다.

## 2026-09-09 — S1 periodic atomic quantum noise

Finite-seed periodic mean에 대응하는 complete-basis atomic drift와 microscopic diffusion을 구현했다. [유도](periodic_noise_derivation.md), [최종 보고서](s1_periodic_noise_report_v2.json), [확인한 그림](s1_periodic_noise_v2.png)에 convention·검증·범위를 기록했다. 이번 단계는 원자층이며 finite-seed field M/D·propagation·S₋는 아직 연결하지 않는다.

`gabes/quantum/periodic.py`는 H(t)/ℏ=H₀+Ve⁻ⁱᵛᵗ+V†eⁱᵛᵗ와 explicit constant jumps를 받는다. Complete traceless Hermitian basis의 A_q는 generator에서, reservoir별 D_q는 Trρ_q[L†,F_i][F_j,L]에서 각각 구한다. Periodic covariance의 시간 미분을 포함한 Ċ=AC+CAᵀ+D와 K̇=AK+KAᵀ+D−Dᵀ를 검사한다. K(t)는 state-dependent atomic commutator이며 canonical field J를 대신하지 않는다.

Harmonic ladder의 drift block은 A_(h−k)+ihνδ_hk, diffusion block은 D_(h−k)다. Lesser ordering은 atomic indices만 transpose한다. Quasifrequency와 physical harmonic frequency를 구분하고, reservoir별 Toeplitz PSD, lifted stability와 전체 ordered spectrum PSD를 확인한다. Nonzero Fourier coefficient 하나의 PSD를 요구하지 않는다. Mean와 response cutoff를 분리해 고정된 내부 h=−1,0,1 blocks를 비교했다. Finite ladder 경계에서 무한계의 covariance identity가 정확히 성립한다고 주장하지 않는다.

공통 `reduced_pump_system`으로 기존 stationary/새 periodic adapter의 Hamiltonian과 reservoirs를 공유한다. 기존 stationary power/readout baseline 테스트가 유지된다. 새로운 adapter는 prescribed local probe/conjugate를 받아 단일 velocity, uniform-Zeeman-RMS reduced atom을 계산한다.

독립 `reference/periodic_qrt.py`는 production A/D/Floquet mean을 사용하지 않는다. Full Liouvillian을 시간 적분해 periodic fixed point와 one-period QRT integral을 얻고, trace-zero 공간에서 모든 이후 주기를 geometric resolvent로 합한다. 임의 finite-time tail cutoff를 사용하지 않으며, periodic coherent mean을 뺀 connected spectrum을 비교한다. General periodic-correlation 배경은 [Navarrete-Benlloch 등 (2021)](https://arxiv.org/abs/2005.08249)을 참고했다.

같은 600 mW pump, 530 μm waist, Δ/2π=0.9 GHz, δ/2π=−8 MHz fixture의 zero seed, 8 μW entrance, 8 μW prescribed exit, 1 mW prescribed exit를 검사했다. 마지막 case는 weak-field exit carrier를 √125배한 stress test이고 self-consistent output이 아니다.

- 모든 15 atomic operators와 3×3 내부 harmonic blocks의 시간영역 QRT 상대 차이는 5.83×10⁻¹³–2.86×10⁻¹²이다. 비교 base frequency는 0, 1 MHz다.
- Periodic state reference의 최대 원소 차이는 3.88×10⁻¹³ 이하, 독립 adjoint product-rule diffusion과 차이는 1.90×10⁻¹⁴ 이하이다.
- Dynamic Einstein/commutator 상대 잔차는 각각 1.50×10⁻¹⁶, 6.34×10⁻¹⁷ 이하이다. Phase-sampled state/ordered covariance minimum은 5.60×10⁻⁴보다 크다. Diffusion의 미소 음수 최소 고유값은 norm 대비 약 10⁻¹⁶이며 clipping하지 않았다.
- Mean order 4→5 변화는 4.00×10⁻¹⁸, response order 4→5의 내부 spectrum 상대 변화는 8.16×10⁻¹⁵ 이하이다. 1 mW에서 response order 1→2 변화는 3.37×10⁻⁵로 실제 cutoff 오차가 드러난다.
- Zero-seed stationary limit과 차이는 3.90×10⁻¹⁵ 이하, reference phase 16→32 변화는 5.63×10⁻¹⁴, ODE tolerance 변화는 1.01×10⁻¹²이다.
- 같은 periodic A에서 D₀만 남기면 내부 spectrum Frobenius 상대 변화가 8 μW exit에서 0.02259%, 1 mW에서 2.4113%다. 이 ablation은 microscopic noise의 대체 모델로 채택하지 않는다.

선택한 probe h=1 lesser atomic-noise spectrum은 8 μW prescribed exit/base 1 MHz에서 zero seed보다 40.02% 증가한다. Base 0.1 MHz에서는 약 3.88배다. 이전 mean polarization 변화 0.00974%와 다른 진단량이다. 작은 mean-state back-action이 noise 변화를 제한한다는 결론을 내릴 수 없으며, 이 atomic component를 검출 S₋의 같은 변화율로 해석하지 않는다. 기존 −0.73 dB 조건부 readout의 finite-seed 보정값을 얻은 것은 아니다.

초기 보고서/그림은 보존했다. v2는 scalar input metadata, 상대 PSD 최소 고유값, 그림의 log ratio 및 heatmap contrast를 보강한 실행이다. 최종 source hashes 31개가 일치한다. 추가 10개 검사는 dynamic jump/product-rule identity, signed beat, zero-seed shifted blocks, time-origin phase, greater/lesser reflection, 시간영역 QRT와 cutoff 수렴, D-average 반례, 비감쇠/부족한 cutoff/잘못된 frequency 거부, Rb report와 immutable artifacts를 다룬다.

`python -m pytest -q tests/quantum`: **93 passed (23.91 s)**. 필수 `python -m pytest -q`: **744 passed, 1 failed (164.48 s)**. 실패는 이전과 같은 `tests/test_docs_consistency.py::test_repository_visibility_wording_is_consistently_public`의 기존 `FWM_physics.tex` 누락이다. 해당 삭제나 테스트는 변경하지 않았다.

**현재 다음 미해결 작업**

1. Periodic atom의 독립 field-drive/readout couplings와 full Nambu·harmonic/mode mapping을 유도한다. 두 sector 사이 및 추가 harmonics를 포함하고 canonical field commutator와 physical-mode uncertainty를 검사한다.
2. Finite-seed field propagation/readout을 연결해 이번 atomic-noise 변화가 S₋에 미치는 영향을 계산한다. Self-consistent carrier/pump depletion과 prescribed-carrier 진단을 구분한다.
3. Full Zeeman signed CG·편광, angular-Doppler/modes, 실제 collection inputs와 non-Gaussian fourth cumulants, 실측 불확도 및 held-out gain/S₋ 검증을 이어간다.

체크리스트의 `periodic_noise_implementation`에 기록했다. Periodic atomic quantum noise는 구현되었으며 finite-seed field quantum noise, absolute hot-vapor prediction과 experimental validation은 false다. Grand Challenge/milestone 1은 진행 중이다.

## 2026-09-10 — S1 finite-seed two-band microscopic field propagation

[유도](periodic_field_derivation.md), [최종 보고서](s1_periodic_field_report_v2.json), [확인한 그림](s1_periodic_field_v2.png)에 기록했다. Periodic atomic A/D를 두 optical carrier band의 full Nambu field로 연결하고, nonlinear probe/conjugate mean과 같은 상태에서 gain 및 bright intensity-difference spectrum을 계산했다. Classically prescribed pump, single velocity, uniform area, zero phase mismatch와 기존 reduced RMS dipoles를 유지한다. 추가 optical ports와 pump depletion은 포함하지 않는다.

`gabes/quantum/periodic_field.py`의 `PeriodicFieldPorts`는 physical label, integer carrier harmonic, atomic lowering operator와 coupling을 명시한다. Full field order는 (p,c,p†,c†), atomic harmonics는 (+1,−1,−1,+1), canonical J는 diag(+1,+1,−1,−1)이다. 내부 atomic harmonic ladder와 physical optical ports를 구분한다. Generic extra-port test는 구현되었지만 그 port를 Rb 장치의 실제 mode로 가정하지 않는다.

Driving B_(h,j)=−ig_j TrF[O_j†,ρ_(h−h_j)]는 density commutator에서, C₀는 Maxwell readout에서 독립 유도한다. M=λC₀RB와 D_field=λC₀R D_atom R†C₀†를 계산하고 reservoir별 양 ordering을 유지한다. Atomic dynamic K의 lifted identity 및 BJ=−K C₀†에서 canonical local relation을 유도한다. Noise를 commutator defect로 역산하지 않는다.

Response order 1은 8 μW midpoint에서 field commutator 상대 잔차 6.94×10⁻⁷로 10⁻⁸ 기준을 실패하며 반환을 거부한다. D의 PSD는 통과하는 사례다. Order 2/3/4는 통과하고 마지막 M/D 변화는 각각 1.07×10⁻¹⁴, 5.40×10⁻¹⁴다. Finite-ladder boundary defect를 새로운 noise로 수선하지 않는다.

`gabes/fwm_quantum/periodic_cell.py`는 β(z)에 따른 periodic mean을 매번 계산하여 nonlinear mean ODE를 적분한다. DC quantum transfer는 mean map의 tangent이지 mean을 직접 전달하는 행렬이 아니다. 네 real input 방향의 독립 finite difference와 DC transfer는 최대 5.69×10⁻¹⁰ 상대 오차로 일치한다. T(0)β_in을 output mean으로 잘못 쓰면 8 μW에서 상대 차이 5.54×10⁻⁵가 생긴다.

처음에는 constant-midpoint quantum segments 4/8/16을 비교했지만 8→16 S₋ 변화 5.19234×10⁻⁶ dB가 선언한 10⁻⁷ dB 기준을 실패했다. [초기 보고서](s1_periodic_field_report.json)와 그림을 보존했다. Dense mean trajectory를 따라 T,N⁾,N⁼를 직접 적분하는 adaptive DOP853을 추가했다. rtol 2×10⁻⁹→2×10⁻¹¹ 비교는 8 μW에서 4.00×10⁻¹⁵ dB, 1 mW에서 2.42×10⁻¹² dB로 통과한다. Midpoint 4/8/16과 adaptive 결과의 차이가 약 4배씩 줄어드는 것도 확인했다. Constant generator에 대한 adaptive covariance는 기존 exact augmented-exponential integral과 대조했다.

Full ±RF Nambu transfer에서 four physical sideband U,V와 added covariance를 구성한다. 이전 두 sector 사이에 없던 normal/anomalous correlations를 포함하며 reflection, channel CP, physical covariance uncertainty를 검사한다. 독립 ordered current 계산과 검출 PSD의 상대 차이는 3.87×10⁻¹⁴ 이하이다. Unequal detector losses, complex response 및 balancing도 검사한다. 기존 stationary Gaussian quadratic correction을 periodic readout에 재사용하지 않는다.

동일한 600 mW/530 μm pump, 12.5 mm, Δ/2π=0.9 GHz, δ/2π=−8 MHz, η=0.85 조건:

| Seed | Probe gain | S₋(1 MHz) [dB] | Pump-state 대비 최대 변화 [dB] |
|---|---:|---:|---:|
| 10 nW | 1.1106558457 | −0.7302942150 | 1.23715×10⁻⁷ |
| 8 μW | 1.1106391761 | −0.7301956922 | 9.89590×10⁻⁵ |
| 1 mW | 1.1086057895 | −0.7181562313 | 0.0121771 |

이전 약 40% 증가는 선택한 atomic component의 변화였다. Full readout에서는 상관이 중요하다. Main/companion cross blocks를 M와 D에서 함께 제거한 진단은 8 μW/0.1 MHz에서 S₋를 −0.73031724→−0.69092746 dB로 바꾼다. 1 mW에서 같은 진단은 +2.56099 dB, full 계산은 −0.71823912 dB다. Block pinching은 비교용 별도 모델이며 fitted correction이나 채택한 microscopic prediction이 아니다.

독립 full-Liouville harmonic source solve와 M 상대 차이는 5.52×10⁻¹⁴, 시간영역 QRT를 readout에 투영한 D와 차이는 4.63×10⁻¹²이다. Adaptive evaluation의 local commutator 잔차는 1.30×10⁻¹³ 이하, global은 5.96×10⁻¹⁴ 이하이다. 최소 channel CP 고유값은 9.79×10⁻⁵, 검출 covariance uncertainty minimum은 0.01223 이상이다. Source hashes 33개가 일치한다.

추가 10개 검사는 zero-seed static sector/channel parity, independent full-Liouville M/time-QRT D, nonlinear mean tangent, full ordered photocurrent, no-atom/zero-length, time-origin phase와 explicit extra port, 부족한 cutoff/잘못된 좌표의 거부, adaptive/exact constant propagation 및 immutable scientific artifact를 다룬다. `python -m pytest -q tests/quantum`: **103 passed (72.34 s)**.

필수 `python -m pytest -q`: **754 passed, 1 failed (215.64 s)**. 실패는 이전과 같은 `tests/test_docs_consistency.py::test_repository_visibility_wording_is_consistently_public`의 기존 `FWM_physics.tex` 누락이다. 해당 삭제나 테스트는 변경하지 않았다.

**현재 다음 미해결 작업**

1. Non-collinear spatial phases와 optical-mode selection을 명시한다. Atomic harmonic 수렴과 추가 physical modes를 버리는 근사의 검증을 구분한다.
2. Velocity-class별 mean/M/D를 같은 geometry에서 연결하고, 독립 class의 noise covariance를 density weights로 합산한다. Angular-Doppler와 transverse quadrature를 독립 수렴시킨다.
3. Full Zeeman signed CG·실제 편광, pump depletion/quantum ports, 실제 collection과 non-Gaussian cumulants, 독립 입력 불확도 및 held-out gain/S₋를 연결한다.

`periodic_field_implementation`에 두-band finite-seed field 완료 범위를 기록한다. Grand Challenge/milestone 1은 진행 중이며 absolute hot-vapor prediction과 experimental validation은 false다.

## 2026-09-11 — S1 kinetic local mean/noise와 thermal pump-state cell

[유도](kinetic_derivation.md)에 moving-atom carrier phases, class별 microscopic noise 합산과 적용 한계를 기록했다. `gabes/fwm_quantum/kinetic.py`는 velocity별 finite-seed 국소 mean/M/D를 계산하고, pump-state weak-field 한계에서는 Maxwell 평균 뒤 cell gain과 bright S₋까지 전파한다. General noncollinear geometry, nonlinear finite-seed thermal trajectory와 slow-envelope atomic transport는 아직 구현하지 않았다.

새로 확인한 조건은 νₚ,ᵥ+ν꜀,ᵥ=(2k₀−kₚ−k꜀)·v다. 기존 single space-time phase solver는 q꜀=−qₚ인 geometry에서 닫힌다. 동일한 반대 vacuum 각도 ±5 mrad에서는 Δk=(0.63791649,0,197.590726) rad/m, 12.5 mm의 axial phase=2.469884 rad, 394.15 K의 loop Doppler rms=38817.87 rad/s가 남는다. 이 입력은 자동 거부하고 mismatch를 없애는 보정값을 넣지 않는다.

각도에 따른 응답을 검사하는 별도 synthetic fixture는 kc=2k₀−kₚ를 명시한다. Conjugate k의 크기는 vacuum 값의 1.00002499903배다. 이는 측정한 refractive index 또는 Maxwell dispersion 해가 아니다. Common z-plane photon flux에서 g_j∝1/√cosθ_j를 mean와 readout에 함께 적용하며, refractive energy normalization은 미검증이다. Lab optical photon energies와 RF Ω는 Doppler 적분 중 고정한다.

독립 class의 M와 D⁾/D⁼는 λᵥ=n A_z wᵥ를 한 번만 사용해 합한다. Streaming과 class별 저장 경로를 대조했다. 같은 velocity class를 0.25/0.75로 나누어도 M와 D는 변하지 않는다. √weight를 먼저 더한 잘못된 noise amplitude sum은 PSD를 유지하지만 canonical commutator 상대 잔차 0.00304713으로 거부된다. Commutator defect에서 D를 역산하지 않는다.

Full 16-component Liouville forced response 및 connected QRT reference와의 상대 차이는 M 1.37×10⁻¹⁴, 양 ordering의 D 각각 1.96×10⁻¹² 이하다. 비교에는 세 discrete velocities와 vz=700 m/s의 near-resonant 원자가 포함된다. Finite-seed 국소 평균 polarization의 independent finite difference와 DC drift tangent는 1.12×10⁻⁹ 상대 차이다. 이 discrete reference 검사를 finite-seed Maxwell 전체의 수렴으로 해석하지 않는다.

Velocity factory는 coplanar geometry에서 vy를 해석적으로 적분하고 vx,vz에 독립 Maxwell Gauss–Legendre grid를 사용한다. Retained square의 확률을 정규화하고 omitted Maxwell probability를 저장한다. 온도 394.15 K는 kinetic 입력이며 density와 reset rate를 자동 추정하지 않는다. Narrow resonance 때문에 작은 확률 tail과 실제 D의 적분 오차는 다를 수 있다.

[초기 보고서](s1_kinetic_report.json)와 [그림](s1_kinetic.png)은 수렴 실패를 보존한다. N_z=768→1024, N_x=32, cutoff=7σ에서 S₋ 변화는 2.61×10⁻⁹ dB지만 D의 최대 RF별 상대 변화는 1.10446×10⁻⁵로 선언한 10⁻⁵를 실패했다. Local/global commutator와 PSD는 모두 통과한다. S₋만 비교했다면 이 미수렴을 놓쳤을 것이다.

`kinetic_refinement_audit`는 parent source hashes 35개를 검증하고 N_z=1536, N_x=48, cutoff=7/8σ 계산으로 종방향 및 tail 비교를 보강한다. Initial report를 덮어쓰지 않고 parent artifact의 hash도 남긴다. Transverse 및 6→7σ 비교는 기존 N_z=1024의 검증이라는 provenance를 유지한다. 허용 오차를 바꾸지 않는다.

다음 유도 좌표는 q=qₚ, Q=qₚ+q꜀인 (h,ℓ) harmonic lattice다. Convective shift는 hν−(hq+ℓQ)·v이며 probe drive=(1,0), conjugate drive=(−1,1)이다. Q=0에서 같은 physical spatial phase를 합쳐 기존 모델로 환원해야 한다. 이 좌표에서 mean/D와 optical projection을 함께 유도하고 직접 이동 궤적의 Liouville/QRT와 비교하는 것을 다음 검증 단위로 남겼다.

Slow-envelope transport와 inter-slice noise correlations, velocity-changing collisions, full Zeeman, 실제 transverse collection, pump depletion, higher-order photocurrent, independent measured input uncertainty와 held-out gain/S₋ 검증은 남아 있다. 체크리스트의 `kinetic_implementation`에 구현과 미구현 범위를 구분하며 Grand Challenge와 milestone 1은 계속 진행 중이다.

최종 [v2 보고서](s1_kinetic_report_v2.json)와 [확인한 그림](s1_kinetic_v2.png)은 모든 선언한 controls를 통과했다. N_z=1024→1536의 최대 D 상대 변화 8.4203×10⁻⁷, 1536 격자의 7→8σ 변화 9.7244×10⁻⁸이다. 기존 transverse와 6→7σ 검사는 각각 1.9748×10⁻¹⁰, 8.4347×10⁻⁷이다. 기준은 그대로 10⁻⁵다. Source hashes 36개와 parent artifact hash를 확인했다.

선택한 1536×48/7σ synthetic closed fixture의 probe gain=1.110091613, conjugate gain=0.113166637, S₋(0.1/1/4 MHz)=−0.728325262/−0.728207393/−0.726408844 dB. Local/global commutator 상대 잔차는 각각 4.53×10⁻¹⁵/2.04×10⁻¹⁵ 이하, 최소 channel CP=1.00765×10⁻⁴, detected uncertainty minimum=0.0122276이다. 이 결과는 fixed-pump weak-field/bright approximation의 조건부 계산이며 실험 −7.8 dB를 재현한 것이 아니다.

새 kinetic 검사 **10 passed (2.71 s)**. 최종 필수 `python -m pytest -q`: **764 passed, 1 failed (261.94 s)**, quantum 113개 모두 포함. 실패는 이전과 같은 `tests/test_docs_consistency.py::test_repository_visibility_wording_is_consistently_public`의 기존 `FWM_physics.tex` 누락이다. 해당 삭제나 테스트는 변경하지 않았다. Refinement command 추가 후 전체 검사를 다시 실행했으며, artifact 덮어쓰기·변경된 source 재사용·refinement 결과의 자기 재비교 거부도 별도로 확인했다.

## 2026-09-11 — S1 two-phase convective atomic mean / noise

Ultra의 정확한 대칭 가속 패치 뒤 사용자의 연구 재개 요청에 따라 진행했다. [유도](spatial_derivation.md), [불변 보고서](s1_spatial_report.json), [확인한 그림](s1_spatial.png)에 결과를 보존했다. 이번 구현은 일반 nonclosed plane-wave geometry의 **국소 원자층**이다. Optical field M/D, thermal finite-seed propagation, gain과 S₋는 아직 연결하지 않는다.

`gabes/quantum/spatial.py`는 (h,ℓ) lattice에서 convective trace-zero mean을 풀고 Liouvillian drift와 jump-product D를 독립 구성한다. Sampled density/covariance positivity, dynamic Einstein/atomic-commutator identities, finite response stability와 reservoir별 양 ordering PSD를 검사한다. 응답은 요청한 resolvent rows만 정확히 풀어 얻으며 내부 source columns나 harmonics를 생략하지 않는다. 이 대수적 계산량 절감의 근거도 주석에 적었다.

`gabes/fwm_quantum/spatial.py`는 실제 k0,kp,kc를 유지하고 q=kp−k0, Q=kp+kc−2k0, drive labels (1,0)/(−1,1)을 사용한다. Q=0이면 Vp+Vc†를 먼저 합친 기존 periodic model로 quotient한다. Q·v=0만으로 Q≠0의 spatial grating을 지우지 않는다. 기존 kinetic field API의 closure guard와 Ultra kernel은 변경하지 않았다.

독립 `reference/torus_qrt.py`는 full-density trajectory relaxation과 connected two-time regression을 실제 두 moving phase로 적분한다. Toy frequencies (2.3,√0.5) rad/s에서 lattice mean과 최대 원소 차이는 3.281×10⁻¹⁰, atomic spectrum의 상대 차이는 6.686×10⁻¹³이다. Mean (5,5)→(6,6), response h와 ℓ, reference phase/history/ODE tolerance를 독립 refinement했다. Reference controls의 최대 변화는 1.249×10⁻¹¹이다. 작은 mean order는 dynamic consistency를 실패하여 거부된다.

Rb v=0, Q≠0에서는 spatial phase마다 독립 full-Liouville periodic QRT를 계산하고 geometric tail을 정확히 합한 뒤 그 결과를 공간 Fourier project했다. Spectrum 상대 오차 8.097×10⁻¹³, mean 최대 원소 오차 5.689×10⁻¹³이다. Spatial/temporal reference phase grids 4×16→8×16→8×24가 수렴한다. 이 reference를 움직이는 모든 Rb velocity의 장시간 QRT 검증으로 확대하지 않는다.

기존 600 mW/530 μm, Δ/2π=0.9 GHz, δ/2π=−8 MHz, reduced RMS dipoles와 reset을 유지했다. Local probe는 8 μW에 대응하고 conjugate는 βc=0.25i βp로 지정했다. 이는 self-consistent exit carrier가 아니다. Vacuum ±5 mrad geometry의 Δk를 그대로 유지해 세 velocities (0,0,0), (150,0,380), (0,0,700) m/s를 계산했다. 두 moving cases의 loop convection은 75,180.163 / 138,313.508 rad/s다.

Loop convection을 0으로 바꾼 뒤 state/A/D를 모두 재계산한 진단에서 atomic spectrum의 양 ordering norm 변화는 최대 4.738×10⁻⁷, mean 최대 원소 변화는 1.945×10⁻⁸이다. 이는 해당 atomic blocks의 sensitivity이며 optical S₋ 오차 상한이나 loop 항 생략의 일반적 근거가 아니다. Mean (3,3)→(4,4) 변화는 2.03×10⁻¹⁵ 이하, 마지막 h/ℓ response refinement는 1.03×10⁻¹⁴ 이하이다.

보고서의 선언한 controls가 모두 통과했고 source hashes 36개와 predecessor hash를 확인했다. 이전 kinetic v2의 수치를 현재 source가 재생산했다고 주장하지 않고 historical context로만 연결했다. 신규 spatial 검사 13개를 포함한 최종 필수 `python -m pytest -q`는 **791 passed, 1 failed (199.08 s)**다. Quantum 126개와 이전 Ultra 가속 검사가 통과했다. 유일한 실패는 이전부터 누락된 `FWM_physics.tex`를 읽는 기존 문서 검사이며 해당 삭제/검사는 변경하지 않았다.

**다음 미해결 작업**

1. Auxiliary atomic (h,ℓ) coordinates를 실제 optical spatial ports에 투영하고, 동일 lab RF에서 B/C/M/D와 공간 noise correlations를 유도한다. 모든 lattice label을 독립 photon mode로 세지 않는다.
2. Q→0 coherent quotient와 기존 field parity, canonical optical commutator/uncertainty 및 physical-mode truncation을 검사한 뒤 general-geometry gain/S₋로 연결한다.
3. Finite-seed thermal mean/noise propagation, velocity/space/z convergence와 slow-envelope/inter-slice transport를 검증한다.
4. Full Zeeman/편광, measured collection, velocity-changing collisions, pump depletion, higher-order current, independent measured inputs와 held-out experimental validation을 잇는다.

체크리스트의 `spatial_implementation`에 완료 범위와 다음 단계를 기록했다. Grand Challenge와 reduced hot-vapor milestone, absolute experimental prediction은 계속 진행 중이다.
## 2026-09-11 — S1 유한 면적의 stationary-center microscopic field

[유도와 결과](spatial_field_derivation.md), [최종 보고서](s1_spatial_field_report_v2.json), [그림](s1_spatial_field_v2.png). 이전 two-phase atom에서 실제 probe/conjugate profiles로 이어지는 첫 canonical field projection을 구현했다. 정지 원자 중심, fixed uniform pump, 두 fixed normalized transverse profiles의 Galerkin model이다. Nonclosed vacuum geometry를 허용하지만 moving thermal field를 완성한 것은 아니다.

`gabes/quantum/spatial_modes.py`는 positive area weights와 두 complex profiles의 각각의 ∫|u|²dA=1을 요구한다. `gabes/fwm_quantum/spatial_cell.py`는 같은 profile을 Hamiltonian drive와 Maxwell readout에 reciprocal하게 적용하고, independent atoms의 M 및 양 ordering D를 n dA로 한 번만 합산한다. E=diag(η,η*)가 J와 가환하므로 projected canonical commutator가 보존된다. Profile projection 자체를 loss channel로 부르거나 commutator defect에서 noise를 만들지 않는다.

공통 time-origin phase를 atom/drive/readout에서 함께 제거해 rapid optical phase별 solve를 없앴으며 full loop phase Q·r는 유지한다. η가 정확히 같은 quadrature rows는 H/B/C/D도 같으므로 area weights를 먼저 합산한다. Zero density에서는 모든 atomic source가 정확히 사라진다. 세 최적화 모두 수학적 동등성 주석을 남겨 불필요한 반복 solve가 다시 들어오지 않게 했다. 기존 Ultra kernel/physical grids/cutoffs는 변경하지 않았다.

`analysis/grand_challenge/reference/spatial_field.py`는 original probe/conjugate phases와 complex unequal profiles를 사용한다. Independent dense full-Liouville forcing의 M 오차는 1.58×10⁻¹³, time-domain QRT 양 ordering D 오차는 최대 3.86×10⁻¹²다. Reference phase 16→24, mean/response (4,3)→(5,3)→(5,4), spatial Nx=2/4/8/12를 따로 비교했다. New adapter의 stationary section은 앞선 two-phase Rb atom의 state/D와도 일치한다.

Nonlinear β(z) mean과 같은 trajectory의 full Nambu quantum tangent/noise를 적분하고 ±RF four-sideband CP, covariance uncertainty 및 direct ordered-current PSD를 검사한다. DC transfer와 nonlinear mean map의 독립 미분은 5.77×10⁻¹⁰ 상대 오차로 일치한다. Canonical reference area A₀를 7배 바꾸어도 gη가 정확히 같아 결과가 보존된다. Fixed-aperture Q→0, collinear old mean/M/D, area splitting, vacuum limits와 unequal complex detector response도 검사한다.

조건은 600 mW/530 µm, n=10¹⁸ m⁻³, L=12.5 mm, seed=8 µW, Δ/2π=0.9 GHz, δ/2π=−8 MHz, explicit reset=2π×100 kHz, η_p=η_c=.85와 400×300 µm top-hat이다. Reset은 이 stationary-center fixture의 declared Markov replacement이며 실제 transit을 예측한 값이 아니다. RF=0.1/1/4 MHz 세 점을 계산했다.

| Probe/conjugate 각도 | Probe gain | Conjugate gain | S₋(1 MHz) [dB] |
|---|---:|---:|---:|
| 0/0 mrad | 1.1106391761 | .1135452814 | −.7301956922 |
| +5/−5 mrad | 1.0620337919 | .0650877602 | −.4328093597 |
| +5/−4 mrad, Nx=12 | 1.0277727399 | .0306514955 | −.2040641891 |

Unequal-angle Nx=4→8의 S 변화 3.5784×10⁻⁶ dB와 gain 변화 5.19×10⁻⁷은 선언한 10⁻⁶ dB / 10⁻⁷ 기준을 실패했다. Nx=8→12는 각각 3.93×10⁻¹³ dB / 2.56×10⁻¹⁴로 통과한다. Quantum ODE rtol 2×10⁻⁹→2×10⁻¹¹의 S 변화는 9.84×10⁻¹² dB, mean ODE refinement의 trajectory 상대 변화는 7.30×10⁻¹²다. 광학 mode completeness는 이 quadrature 수렴으로 입증되지 않는다.

Moving-grating negative control은 v=(150,0,380) m/s, +5/−4 mrad, 강한 prescribed β=√125×(5×10⁶,.5i×10⁶)를 쓴다. Moving torus를 푼 뒤 grating phase만 고정하면 ordered D는 PSD지만 field commutator 상대 잔차는 2.96091344×10⁻⁵로 실패한다. Mean (4,4)→(5,5), response 3→4에서도 잔차가 약 4.14×10⁻¹⁵만 변한다. Dynamic K identity에서 ω₂∂θ₂K를 잃은 잘못된 reduction이다. Adapter는 모든 nonzero velocity를 거부하며, 기존 kinetic closure guard도 그대로 남긴다. 이 반례는 모든 moving-medium reduction의 불가능성 증명이 아니다.

초기 `s1_spatial_field_report.json`도 해당 버전의 controls를 통과했다. V2는 cutoff-refined negative control, consumed scalar/derived ledger, 더 완전한 dependency manifest와 실행 전후 source stability를 추가해 전체 계산을 다시 실행한 결과다. 초기 evidence는 덮어쓰지 않았다.

V2의 모든 controls 및 실행 전후 source stability가 통과했고 SHA-256 **41개**와 이전 atomic evidence hash를 재확인했다. Ultra verification의 source hashes도 일치한다. 그림을 렌더링하여 확인했다. 신규 spatial-field 검사 **15 passed**; 최종 `python -m pytest -q`는 **856 passed, 1 failed (231.76 s)**다. 별도 downstream 신규50개도 포함한다. 실패는 이전부터 누락된 `docs/FWM physics and analytic reconstruction/FWM_physics.tex`를 읽는 문서 consistency 검사 하나이며 해당 삭제 상태나 테스트를 수정하지 않았다. 재현: `python -m analysis.grand_challenge.spatial_field_audit --output NEW.json --plot NEW.png`.

다음 구현은 moving-atom characteristic transport의 작은 독립 검증 문제부터 시작한다. `(∂t+v·∇)δF=AδF+Bδb+f`에서 envelope/profile의 궤적상 변화, aperture inflow state와 reservoir를 함께 정하고, retarded field-response kernel 및 서로 다른 z의 noise correlations를 도출한다. Atomic inflow의 number statistics도 명시해야 한다. 시간영역 직접 전파와 비교하고 v→0에서 이번 stationary-center field로 환원되는지 검사한 뒤 moving Rb/velocity quadrature로 확장한다. Local z-white noise를 가정한 채 grating만 고정하는 경로로 되돌아가지 않는다.

Full signed hyperfine/Zeeman·편광, diffraction/walkoff/추가 optical modes, collisions, pump depletion, higher-order photocurrent, 독립 입력 계측과 held-out gain/RF spectra는 계속 남는다. 체크리스트 `spatial_field_implementation`에 기록했으며 reduced hot-vapor milestone 및 Grand Challenge는 진행 중이다.

## 2026-09-11 — 별도 downstream RF 대역 / 상관 입력 불확도

사용자가 별도로 허용한 병렬 작업의 [설명](parallel_development.md), [v2 보고서](s2_downstream_report_v2.json), [그림](s2_downstream_v2.png)을 체크리스트 `downstream_analysis_implementation`에 연결했다. 신규50 tests와 보고서 controls가 통과하고 source hashes25개를 재확인했다. 위 spatial solver와 독립적인 분석 도구이며 새 spatial field의 불확도를 실행한 결과는 아니다.

SQL/target-dB connected bands, crossing refinement, censored/missing/unconverged edges와 intrinsic width 부재를 처리한다. Correlated Gaussian samples 및 J C Jᵀ propagation을 구현하고, 완전 상관 null modes와 거의 완전 상관의 실제 양의 분산을 구분한다. 기존 reduced pump-state model에는 pump power/waist 각각1%와 correlation .6을 **가정**한 joint gain/S₋ sensitivity만 연결했다. 실측 불확도나 일반 band topology의 통계, 실험 validation으로 승격하지 않는다.
## 2026-09-11 — S1 characteristic transport와 boundary / cross-position noise

사용자가 다음 단계 및 기존 병렬 작업의 추가 진행을 요청했다. [수송 유도](transport_derivation.md), [최종 보고서](s1_transport_report_v2.json), [그림](s1_transport_v2.png)에 prescribed finite atomic characteristics와 별도 passive canonical optical bridge를 기록했다. Main write set은 신규 `gabes/quantum/transport.py`, `analysis/grand_challenge/transport_fixtures.py`, `reference/transport_qrt.py`, `transport_audit.py` 및 신규 tests다. 기존 spatial/kinetic guard나 production Ultra 코드를 변경하지 않았다.

`BallisticPath`는 실제 entry/velocity/residence time을 명시한다. ρ_in에서 출발한 time-dependent Lindblad mean을 complete Hermitian basis의 A와 독립 jump-product D에 연결한다. Covariance 초기조건 C_in과 각 internal reservoir를 분리해 전파하며, 같은 원자가 다음 수치 구간에 들어갈 때 초기화하지 않는다. Source 합을 density matrix에서 직접 계산한 C/K와 대조한다.

`two_time_covariance`는 같은 원자의 두 위치가 공유하는 inflow와 jump history를 유지한다. `retarded_kernel`은 causal C_out(a)R(a,b)B(b)를 반환한다. Finite readout wavepacket의 전체 double-time integral은 Y/τ를 포함한 augmented covariance ODE로 계산하며, SI microsecond 단위의 작은 covariance를 absolute tolerance가 놓치지 않게 한다. Exact propagator composition, constant jump-product cache와 Y/τ coordinate change의 근거를 주석에 남겼다. 감쇠한 fundamental matrix의 inverse를 사용하지 않는다.

Independent full-density QRT는 production의 A,D나 covariance를 받지 않는다. Actual Liouvillian에서 positive-time triangle, reversed ordering 및 driven density response를 직접 적분한다. 2준위 moving fixture의 QRT 오차는 양 ordering 최대 **1.77×10⁻¹³**, retarded response **3.42×10⁻¹⁴**, DC nonlinear mean finite difference **4.03×10⁻¹⁰**이다. Source 합 C/K의 오차는 **1.25×10⁻¹⁴ / 3.80×10⁻¹⁶**이며 3준위 transient 및 SI/time-rescaling 검사가 추가되었다.

조건은 synthetic 2준위, r_in=(−300,0,0) µm, v=(150,0,100) m/s, τ=4 µs다. Spatially varying Gaussian drive와 phase, 두 complex readouts를 사용하고 −.5…+.5 MHz의21 samples를 계산한다. Reference 비교는 −.3/0/+.3 MHz다. Rb/실험 파라미터가 아니며 declared reset bath는 ballistic renewal과 별개다.

Sampled full covariance의 직접 double integral은 age samples17/33/65/129에서 상대 오차 .005931/.001472/.0003672/.00009176으로 수렴했다. 같은 kernel의 cross-position blocks를 지우면 올바른 norm의 .2548/.1274/.06372/.03186만 남으며 grid를 늘릴수록 잘못 0에 접근한다. Inflow를 지우는 다른 대조군은 나머지 noise가 PSD여도 atomic commutator defect **.20196**, spectrum norm 변화 **14.59%**를 만든다. 해당 atomic fixture의 수치이며 hot-Rb squeezing 오차율은 아니다.

Explicit independent Poisson entry의 beam PSD는 J(C_Y+mean_pulse mean_pulse†)다. Mean-pulse 항은 number statistics로 정해지고 별도 noise fit coefficient가 아니다. Identity readout은 internal noise=0이지만 count noise가 남는 해석 대조군과 일치한다. 모든 원자가 같은 stationary age protocol을 겪는다는 제한이 있으며, CW double-Λ에는 entry-phase average가 추가로 필요하다. Constant atom/fixed length/v→0의 finite-residence spectrum은 기존 stationary atomic spectrum에 수렴하지만 nonlinear spatial FWM 전체의 limit를 검증한 것은 아니다. Nonzero mean의 DC transit peak는 delta limit일 수 있다.

병렬 검증 에이전트는 [독립 passive transport bridge](ballistic_linear_channel.md)를 신규 reference/test/doc로 구현했다. Optical/atomic inflow와 distributed reservoir의 commutator 기여는 각각 .2945585873/.3979254850/.3075159276이며 합은1, full residual **2.23×10⁻¹⁶**이다. Independent Gramian/ODE와 nonconstant-profile refinement 및 nonlocal source-kernel 적분을 확인했다. Inflow 제거는 CP minimum −.1989627425, reservoir cross-z covariance 제거는 commutator deficit 약 .30568을 만든다. Passive oscillator testbed로서 finite Rb atom과 구별한다.

[2026-08 Liouville–transport 선행 연구](https://arxiv.org/html/2608.15130v1)도 확인했다. Finite mode와 boundary renewal을 함께 다루는 Cs D2 noise 연구이며 common scale factor와 longitudinal propagation 생략을 명시한다. 이를 no-fit hot-Rb FWM squeezing의 완료로 인용하거나 수송 연구 자체의 선행 부재를 주장하지 않는다.

초기 transport report/plot은 보존했다. Final v2는 추가 SI/3준위 검사를 포함한 source snapshot으로 다시 계산했고, 모든 controls 및 source stability가 통과했다. **12개 source hashes**와 predecessor artifact hash를 재확인했으며 그림을 확인했다. 신규 atomic16 +passive24 tests 통과. 필수 전체 `python -m pytest -q`는 **961 passed, 1 failed (216.03 s)**다. 별도 spatial uncertainty33 tests도 포함한다. 실패는 기존 `FWM_physics.tex` 누락 한 건이며 삭제나 검사를 변경하지 않았다.

다음은 flux-weighted boundary/path/velocity measures 및 entry-phase average를 정의하고, moving periodic Rb의 빠른 내부 dynamics와 느린 envelope transport를 결합하는 것이다. 같은 characteristics에서 nonlocal Maxwell response/noise를 구성한 후 optical commutator/CP와 stationary spatial-field limit를 검증해야 한다. Full atom, optical modes, collisions, pump depletion 및 독립 실측/holdout은 남아 있다. `transport_implementation`에 기록했고 Grand Challenge 전체는 진행 중이다.

## 2026-09-11 — 병렬 후속: stationary spatial field의 joint uncertainty

기존 병렬 작업에 추가로 [새 공간 모델 불확도](spatial_uncertainty_derivation.md)를 맡겼으며 [보고서](spatial_uncertainty_report_v1.json), [그림](spatial_uncertainty_v1.png)을 확인해 `spatial_uncertainty_implementation`에 통합했다. 신규 analysis runner/test/doc/artifacts만 추가되었고 기존 uncertainty 또는 solver는 수정하지 않았다.

σP=.006 W, σθ_c=10 µrad를 **가정**하고 correlation=0 및 ±.6을 비교한다. Angle마다 actual wavevectors와 같은 microscopic mean/M/D를 다시 계산한다. Nx8 coarse/fine finite differences와 Nx12 fine의 **14 unique cell solves**가 모두 transfer/CP/uncertainty 관문을 통과했다. Step refinement의 output-covariance 상대 변화는 **7.46×10⁻⁷**, spatial refinement는 **9.59×10⁻¹⁰**이다.

Nx12에서 (Gp,Gc,R_.1MHz,R_1MHz,R_4MHz)=(1.0277727399,.0306514955,.9540881696,.9540993085,.9542703241), standard uncertainties=(.00060807,.00060870,.00092742,.00092727,.00092503)이다. 1 MHz의 표준불확도는 correlation −.6/0/+.6에서 .00106715/.00092727/.00076214로 바뀐다. 동일 Jacobian의 signed cross term 효과이며 correlation을 data에 맞추지 않는다. 세 RF samples의 intrinsic bandwidth는 모두 None이다.

Parent stationary-field 보고서 이후 kernels/doppler/fwm의3 sources가 달라진 점을 명시했다. 새 계산의 nominal gain/dB는 parent와 정확히 같지만 parent source parity를 주장하지 않는다. Current source **34개** hashes와 stability를 확인했고 신규 **33 tests**도 전체961pass에 포함되었다. Perturbed points마다 harmonic/ODE refinement를 재실행한 것은 아니며 실측 입력 불확도, nonlinear probability distribution, moving thermal field 또는 experimental no-fit validation으로 승격하지 않는다.


## 2026-09-14 — Moving reduced Rb finite segments와 thermal inflow

사용자가 다음 단계 및 기존 병렬 작업의 추가 진행을 요청했다. 본 작업은
`gabes/quantum/segmented_transport.py`, `gabes/fwm_quantum/transport.py`,
`analysis/grand_challenge/rb_transport_audit.py`와 신규 테스트를 추가했다.
독립 subagent는 full-density raw QRT reference/test/doc를 작성했고,
기존 병렬 작업은 `gabes/quantum/inflow.py`, inflow audit/test/doc를 추가했다.
현재 branch에서 기존 사용자 변경을 보존했으며 production Ultra는 수정하지 않았다.

### Exact frozen-segment characteristic

Prescribed moving pump-only reduced Rb D1의 실제 entry density, four radiative
jumps, Gaussian pump profile midpoint, nonclosed vacuum geometry와 모든 carrier
Doppler phases를 보존한다. Boundary renewal과 중복되는 Markov transit reset은
거절한다. Finite seed power/phase, density와 cell length는 single-path solve에서
미사용으로 명시한다. Weak Nambu ports의 drive/readout coupling은 같은 dipole과
z-plane flux convention을 사용한다. Maxwell field channel은 아직 없다.

Complete-basis covariance에 rotated finite-pulse accumulators를 붙이고,
각 microscopic jump-product D가 actual rho에 선형인 성질로 homogeneous block
exponential을 구성했다. Both orderings와 boundary/internal sources를 따로
전파하고 cross-segment covariance를 유지한다. Mean과 retarded response는 별도
작은 linear lift다. Density-generator zero modes나 decaying map의 역행렬을
요구하지 않으며 PSD/commutator defect로 D를 만들지 않는다.

동일 H/O/V/dt에서는 entering state와 무관하게 propagator가 같으므로
exponential을 재사용하되 실제 state는 계속 전파한다. Seven identical segments
control은 14개 중 12개의 중복 exponential을 제거한다. Pump-only entry phase는
s=(1,-1,-1,1)의 정확한 row/column phase이므로 uniform common-phase 평균을
Kronecker delta(s_j,s_k)로 수행한다. 이 두 수학적 근거와 finite-seed 제한은 코드
주석/테스트에 남겼다. 여러 RF와 여러 구간을 동시에 다루는 2-/3-level 독립 QRT
regression도 추가하여 같은 segment entry density가 모든 RF에 사용되도록 검증한다.

### Independent evidence and failed envelope gate

Independent `reference/segmented_qrt.py`는 A/D/propagated covariance를 사용하지
않고 raw O-dagger*rho 또는 rho*O-dagger sources의 full-density QRT triangle을
전파한다. Global mean outer product는 끝에서 한 번만 뺀다.

Immutable evidence:
- `docs/grand_challenge/rb_transport_report_v1.json`
- `docs/grand_challenge/rb_transport_v1.png`
- `docs/grand_challenge/rb_transport_derivation.md`
- `docs/grand_challenge/segmented_qrt_reference.md`

2 microseconds, v=(150,10,100) m/s, entry=(-100,20,0) micrometers,
P=.6 W, w=530 micrometers, vacuum angles .006/-.005 rad의 declared control이다.
RF .1/1/4 MHz와 Gaussian segments8/16/32/64/128을 계산했다. Constant physical Rb
QRT max relative error1.183e-11; constant1->4 refinement1.198e-12; moving8-segment
QRT3.481e-12; actual phase re-solve vs theorem1.327e-13. Direct density covariance
vs source sum max4.302e-13. Boundary omission changes wavepacket norm18.06% and
leaves atomic exit K defect2.896e-5; deleting cross-segment covariance changes
norm333.46% despite positive remaining noise. These are atomic fixture norms,
not optical squeezing errors.

Smooth Gaussian envelope is UNCONVERGED. At64->128, relative changes are
[greater .00399397, lesser .00372007, mean .02955849, retarded .00010738,
exit .00012832]. The declared criterion1e-3 requires two consecutive refinements
for all five quantities. The report separately records exact controls=true,
smooth convergence=false, all certification gates=false. Numerical variation is
not yet attributed to a specific missing physical mechanism.

### Parallel thermal inflow

`inflow.py` integrates n*f(v)*|v.normal|dA dv over all six box faces with incoming
Rayleigh normal speeds, Gaussian tangential speeds, positive weights and the
true first-exit chord. Occupation is not rescaled to nV. One common laboratory
entry phase retains fixed spatial offsets. Marked-Poisson averaging uses
conditional connected moments plus the conditional mean outer product BEFORE
averaging; lab-periodic output is the zero-cyclic connected spectrum, without
coherent mean delta lines. Actual Rb ensemble integration remains open.

Immutable `inflow_report_v1.json`/`inflow_v1.png` and `inflow_derivation.md` record
T373K, n1e17 m^-3 and equal-volume cube/thin-box controls. Three seeds and
six faces times2^10..2^18 points give final maximum occupancy error0.01215%,
scaled phase-space moment error0.03228%, last identity-spectrum refinement0.8771%.
Wrong half-normal normal velocities produce occupancy excess14.14%/69.66%.
Wrong independent phases and deleting Poisson means also fail their controls.
This certifies the tested measure and identity observables, not a thermal Rb
polarization integral or an optical SQL-normalized field spectrum.

### Verification and next work

All current Rb24 source hashes and inflow11 source hashes match their immutable
run manifests. Both figures were rendered and inspected. Final full
`python -m pytest -q`: **1061 passed, 1 failed (323.85 s)**. New tests: transport17,
independent QRT31, inflow43, total91. Only failure is the pre-existing missing
`docs/FWM physics and analytic reconstruction/FWM_physics.tex` in the repository
visibility documentation check. Its deletion and the test remain unchanged.

Next: certify the smooth-envelope characteristic, connect the actual Rb
polarization wavepackets to thermal flux/path/phase averaging with independent
spectral convergence, add finite-seed periodic dynamics, and close nonlocal
Maxwell response/noise with optical commutator/CP and stationary-limit checks.
The no-fit experimental Grand Challenge remains in progress.


## 2026-09-14 — Continuous Gaussian characteristic and guarded ensemble adapter

The previous frozen Gaussian refinement left the selected physical Rb path
unconverged. This step preserves the same reduced pump-only Hamiltonian,
dipoles, four radiative jumps, moving path and carrier phases, and integrates
H(a)=H0+f(a)H1 continuously. No fitted loss, diffusion factor or experimental
squeezing target is used.

### Continuous microscopic D and exact reuse

`gabes/quantum/smooth_transport.py` assembles sparse affine moment generators
once and integrates the actual density, independently derived jump-product D,
source covariances, cross-time readout histories, mean pulse and retarded
response with DOP853. `gabes/fwm_quantum/smooth_transport.py` binds the same Rb
ports and Doppler/entry phases to a continuous transverse Gaussian envelope.

Equal-time atomic C_r is independent of analysis RF, and its other Hermitian
basis ordering is exactly C_r.T. One shared atomic solve per source replaces
2*nf identical solves. Exact accumulator coordinates eta_j=1+abs(w_j)*T are
undone in physical outputs. Mathematical comments and regressions preserve
these identities: 5461 complex variables instead of 11086 duplicated variables.
This is a variable-count reduction, not a measured production Ultra speedup.

Parallel independent `analysis/grand_challenge/reference/smooth_qrt.py` uses
full-density raw regression triangles, not A or D, and subtracts the global
mean outer product at the end. It verifies total ordered noise, mean and exit
state. Independent source-resolved and retarded-response references for this
physical 2 us Rb trajectory remain to be added; shorter/generic controls do not
supply that stronger coverage automatically.

### Failed reference refinement retained; unchanged criteria then passed

Immutable artifacts:

- `smooth_transport_report_v1.json`: initial reference refinement FAILED.
- `smooth_transport_report_v2.json` and `smooth_transport_v2.png`: tighter QRT
  and unchanged primary physics; all declared selected-path controls PASS.
- `smooth_transport_derivation.md` and `smooth_qrt_reference.md`: derivation,
  independence and scope.

The same entry=(-100,20,0) micrometers, v=(150,10,100) m/s, duration2 us,
P=.6W, waist530 micrometers and RF .1/1/4 MHz were used. Primary refinements
(3e-9,3e-12)->(1e-9,1e-12)->(3e-10,3e-13) give maximum relative changes
1.185e-7 and4.136e-8, both below1e-3. Initial QRT mean refinement2.036e-5
failed its2e-6 budget. Tightening QRT to(2e-12,2e-16) and(2e-13,2e-17)
reduced that change to2.013e-7. Primary versus finest independent QRT errors:
greater4.364e-10, lesser4.213e-10, mean4.640e-9, exit3.873e-12; all below5e-6.

The v2 runner verifies all27 original source/test hashes and parent primary
tolerances before reusing those identical primary results. It retains the
failed parent report SHA256 and old QRT cases; its own source makes28 hashes.
A changed source or existing output is refused. This avoids redundant solves
without removing independent evidence or changing acceptance criteria.

Finest primary atomic covariance closure residual2.446e-10, commutator
residual7.931e-15; density and source-resolved ordered PSD checks pass without
clipping. Against this continuous solution, the old128-segment mean differs
1.8371%; greater/lesser differ1.457e-4/1.448e-4. This identifies the numerical
envelope error for the same physical submodel; it does not identify omitted
experimental physics. The old frozen report remains unchanged and failed.

Primary wall times349.13/407.37/465.35s; tighter independent QRT119.78/159.39s
while other jobs ran. The primary finest required1,770,761 RHS evaluations.
These costs motivate further validated acceleration before thermal quadrature.
No approximation or new fast method has yet been certified across that domain.

### Verification

New tests: main continuous transport18, independent QRT40, parallel guarded
ensemble adapter66; total124. Final `python -m pytest -q`: **1185 passed,
1 failed in256.50s**. The sole failure still references the previously deleted
`docs/FWM physics and analytic reconstruction/FWM_physics.tex` in
`test_repository_visibility_wording_is_consistently_public`. The deletion and
test were preserved. The earlier full run also gave1185passed/1samefail317.49s.
The supplementary audit runner passed current-manifest validation, stale-source
refusal and immutable-output refusal checks. All28 current source hashes and
the immutable parent report hash match; the final figure was visually inspected.

Selected-path numerical convergence is now established within the stated
controls. Actual thermal Rb source/response reference coverage and integrand
convergence, finite-seed dynamics, nonlocal Maxwell optical commutator/CP,
full-atom/collision/mode/depletion ablations and out-of-sample experiment remain
open. The no-fit Grand Challenge and production squeezing remain incomplete.


### Parallel ensemble result and remaining quadrature failure

`gabes/fwm_quantum/transport_ensemble.py` now connects supplied pump-only packets
to positive thermal boundary rates. It checks actual common lab RF, path entry
phase/Doppler, Nambu modes, units, reciprocal coupling and density. Common phase
averaging is the exact charge mask; retain mask*(mean*mean-dagger), not the outer
product of the phase-averaged mean. One callback per consumed path avoids phase
re-solves. Hash-bound observed path errors and actual-integrand ensemble errors
are separate requirements. A failed/missing path returns no partial spectra;
no dropped-path reweighting is permitted. These checks validate evidence
contracts, not the truth of an arbitrary caller's claimed independence.

Final `transport_ensemble_report_v1.json`, `transport_ensemble_v1.png` and
`transport_ensemble_derivation.md` use a declared constant two-level analytic
fixture. All63 current source hashes match; source hashes were stable during
the run. Implementation, four-phase sums, actual ODE comparison, streaming and
negative controls PASS. Path Gauss24/48 versus closed-form maximum error1.84e-15;
actual generic transport error9.754e-16. Across powers4/6/8/10 and two seeds,
16,320 path callbacks were made; final6,144 paths per seed.

The final toy ensemble remains UNCONVERGED: last actual-integrand refinement
6.515211% exceeds the predeclared5% budget, despite an independent scramble
change4.064197% and occupancy error0.110193%. Thus
implementation_controls_passed=true, expected_controls_passed=false,
final_candidate.certified=false. The report and nonzero runner exit retain this
failure. It is not a thermal Rb ensemble or an optical squeezing calculation.

Deleting Poisson mean outer terms changes the raw greater norm29.86258%; taking
the phase average before the outer product deletes100% of that number term.
An actual native reduced-Rb packet passes its own internal quantum audit but is
correctly refused without explicit path evidence: spectra=None, certified=false,
one consumed path. The toy has no radiative jumps, so it cannot establish
nonzero jump-source ensemble convergence. Rate-integrated response also does
not supply the spatially nonlocal Maxwell M kernel.

The parallel figure was visually inspected and its gates/remaining work added
to checklist.json. Complete independent physical-Rb source/response coverage,
convergence-preserving acceleration and actual thermal integrand refinement are
the next actionable work, before finite-seed/nonlocal optical closure.


## 2026-09-14 — 잡음원별 adjoint 검증·앙상블 수렴 확장

### 독립 역방향 계산

추가: `reference/adjoint_transport.py`, `adjoint_transport_audit.py`.
전방 full-density 1회; 역방향 관측자 K·source 적분·복소 응답 전파.
기존 forward drift/D/covariance 입력 없음. 실제 H·jump 정의만 공유.

Lindblad product 항등식으로 boundary + 각 bath source 분리.
모든 source에 전체 jump가 켜진 동일 density/관측자 사용.
Jump 제거 후 차감하는 ablation과 구분. Mean outer는 boundary에서 한 번 차감.

정확한 좌표: eta=1+abs(w)*T, K=eta*exp(-i*w*t)*A/T.
복원에서 eta·T 모두 취소. 잡음·효율 fitting 없음.
RF/source와 무관한 density solve 재사용 근거 주석 포함.
실제 Rb: forward16 + backward624 complex variables. 기존 forward lift5461.
변수 수 비교만; 생산 Ultra speedup 주장 없음.

24 reference tests + 4 evidence tests 통과.
해석적 ground-state boundary/vacuum noise, 복소 응답 부호·주파수,
물리적 cosine-drive finite difference, jump phase/splitting/permutation,
identity shift, seconds 재척도, 다중 source/RF, 변경 source·output 거부 포함.

기존 smooth v2 primary3해: 28source hash·입력·metadata 일치 후 재사용.
새 source별/RF별 norm에서 두 successive 최대 변화2.146e-7/7.526e-8.
mean_outer 포함. 둘 다1e-3 기준 통과. 물리적 adjoint 비교 수치는 아래 추가.

### 병렬: 이전 toy ensemble 실패 해결

신규 `transport_ensemble_refinement.py`, tests51개, 별도 한국어 정리 MD.
Immutable `transport_ensemble_refinement_report_v1.json` 및 PNG.
기존 실패 `transport_ensemble_report_v1.json` 보존.

물리 고정: H=0, no jumps, fixed state/readout/drive, zero offsets.
F48/T48 scalar 합으로 정확히 분리. 증명 주석·기존 callback·4phase parity.
이 분리식의 실제 Rb 적용 금지. Nonzero jump-source thermal 검증 아님.

사전계획 p10/12/14/16 × seeds11/211/811.
총1,566,720 actual boundary paths. 최종393,216/seed.
각 path24/48/closed-form 7metric1e-9 유지.
모든 grid closed-form stream6metric2e-12 유지.
각 seed 최근두grid edge + 마지막두grid 모든 directed seed pair: 모두5% 요구.

p14 gate: seed211의 p10→12 변화6.709%로 실패.
p16 gate: 최근두edge 최대2.3116449%, scramble 최대0.88354265% → 통과.
마지막p14→16 edge 최대0.38575%.
모든 path 오차최대4.5284e-14; closed-form stream3.2950e-16.
기존p10 report 재현2.943e-15. Output·density·mean outer 보정 없음.
전체 sweep9.56s; fixture 전용 batch 계산의 실행시간. 실제 Rb 속도 배수 아님.

17source hash·parent SHA256 일치. Figure 직접 검수.
`constant_atom_ensemble_converged=true`.
`production_adapter_certified=false`: 고해상도 callback packet-digest chain 미생성.
`physical_Rb_ensemble_certified=false`: 실제 Rb 경로 계산 범위 밖.

후속: 실제 Rb source/response 경로 evidence + thermal integrand 수렴,
검증 오차 보존 가속, 위치별 Maxwell kernel, finite seed.


### 전체 검사

`python -m pytest -q`: **1264 passed, 1 failed, 232.46s**.
신규79: adjoint24 + evidence4 + 병렬refinement51.
유일 실패: 기존 삭제된 `docs/FWM physics and analytic reconstruction/FWM_physics.tex`
참조 `test_repository_visibility_wording_is_consistently_public`.
삭제 상태·테스트 유지. 이번 변경과 별개.


### 실제 Rb source·response 결과

`adjoint_transport_report_v1.json`: complete path evidence 및 전체선언gate 통과.
Source34개·parentreport SHA256 일치. 기존2us/3RF 조건 유지.
Primary vs adjoint: source별greater7.804e-9, lesser6.928e-9,
mean2.703e-8, mean_outer4.069e-8, retarded1.022e-11, exit5.528e-9.
Independent 자체refinement 최대5.562e-8 <2e-6.
Primary 독립일치 최대4.069e-8 <5e-6.

실제5source(입사경계+4radiative jumps)의 signed PSD 통과.
Fine최소greater2.29358e-22/lesser2.29649e-22s²;
anti-Hermitian residual최대1.488e-34s², roundoff floor3.790e-26s².
출력 clipping·보정 없음. Figure 직접검수.

총QRT-only evidence·변경response의stale digest 거부 확인.
이제 선택된 경로에는 앙상블 계약이 요구하는 source/response coverage 확보.
다른 경로·전체 thermalRb ensemble·광장 Maxwell 인증으로 확대 금지.

Reference coarse478.85s/fine553.92s. Fineforward1,161,104 RHS,
backward2,016,269 RHS. 현재 병목: 광학/carrier 진동을 직접 시간분해하는 비용.
다음 작업: 같은 검증 오차를 유지하는 가속 후 실제 thermal integrand 수렴.
No-fit Grand Challenge 계속 진행 중; 실험−7.8dB·absolute spectrum 미검증.

## 2026-09-14 — 같은 오차 기준의 경로 가속·thermal Rb 준비

한국어 caveman 정리 유지. 새 solver·audit만 추가; 기존 report/source 불변.
Spectral convolution: density16 + augmented19 eigenmodes.
489×489 source lift, 96×96 first-moment lift의 정확한 분리식.
RF/source가 공유하는 density 및 eigenproblem 중복 제거. 증명 주석 유지.
불안정 basis는 기존 block exponential. Noise fitting·clipping 없음.

Gaussian 시간 적분은 전체 affine lift CF4. 작은 matrix exponent만 정확하다고
연속 경로 수렴으로 인정하지 않음. Coarse N128에서는 quantum audit 통과해도
mean_outer 약1.97% 오차. N4096도 약0.080%로 독립5e-6 기준 실패.
해결: 같은 물리식을 직접 refinement. 별도 mean recentering/hybrid 미도입.

History 기본17개 macro endpoint. 모든 내부 substep audit는 streaming으로 유지.
고해상도 covariance history의 GB 단위 저장 회피; 검증점 축소 없음.

신규 solver73 + evidence18 검사 통과. 공식3해 audit 진행 후 결과 추가.
병렬 작업: Maxwell boundary의 실제6경로 × reference2tol,
입력·source hash에 묶인 재사용 cache, actual-rate 재결합 gate.
열적 integrand grid/scramble 수렴은 별도. 유리한 경로만 선택 금지.

### CF4 공식 결과

`exponential_transport_report_v1.json`: 전체선언gate 통과.
Source38개·독립 adjoint report hash 일치. 실행 중 소스 불변.
N4096/8192/16384 =63.146/129.600/250.310s.
두 successive edge 최대9.41853e-4,1.61380e-4 <1e-3.
첫 edge 기준 근방; 다른 thermal path에는 별도 실제 비교 필수.

최종 독립 비교: greater8.00745e-10, lesser8.41531e-10,
source greater2.65216e-9/lesser2.48047e-9,
mean5.73269e-7, mean_outer8.69080e-7,
complex response8.89983e-12, exit4.10696e-8. 모두5e-6 이내.
Complete selected-path evidence 승인; 변경 response digest 거부.
독립 수식·코드검토 blocker0. 수렴 figure 검수.

동일 physical midpoint16 구간 직접비교:
old14.89924s, spectral0.155581s →95.7649배.
8metric 차이 최대1.76453e-11. 같은 coarse 이산화·프로세스·BLAS 조건.
생산 Ultra·전체 연속 경로·전체 thermal ensemble 속도 배수 아님.
당시 actual BLAS 수 미기록인 과거 adaptive report와 공정 속도비 미주장.
신규 timing은 NumPy/SciPy OpenBLAS 모두 실제1 thread 확인.

### 합동 전체 검사

`python -m pytest -q`: **1487 passed, 1 failed, 338.88s**.
신규223 = solver73 + selected evidence18 + thermal/cache83 + thermal audit49.
유일 실패: 기존 삭제된 `FWM_physics.tex`를 읽는
`test_repository_visibility_wording_is_consistently_public`.
기존 삭제·테스트 유지. 새 코드 검증 실패 없음.

### Thermal Rb v1: 실패 보존, 같은 물리로 refinement

고정 p0/seed11 경로6개 × independent reference2개 완료.
CF4 primary3개씩 계산 후 `rb_thermal_ensemble_report_v1.json` 생성.
Source62개 불변, elapsed531.39s. Path0/3/5 통과, path1/2/4 실패.
모든 finest-vs-independent는 통과; 문제는 coarse successive edge.
Path1 mean_outer2.25747e-3 >1e-3.
Path4 두 edge1.01347e-2,2.55107e-3 >1e-3.

허용오차·경로·물리계수 그대로. 공통 maxdt를 단계별 절반으로 줄임.
기존 두 수준의 동일 계산은 검증된 cache 재사용, 각 추가 수준6개만 새 계산.
v2도 기존 path4 edge 때문에 미달 예상. v2 결과 보존 후 v3까지 검증.
코드·test 불변이므로 전체검사 반복 불필요. 최종 observed gate는 결과로 판정.

`rb_thermal_ensemble_report_v2.json`: path4의 첫 edge2.55107e-3만 실패.
나머지5경로 통과. 새 두 번째 edge는 전6경로 최대2.34323e-6 <1e-3.
최종 독립 비교 최대1.03484e-7 <5e-6.
Source62개 불변, elapsed660.885s, primary18개 중12cache hit·6개 새 계산.
v1/v2 실패 그대로 보존. v3 finest6개는 workers4로 병렬 선계산.

### Thermal Rb v3: 6경로 증거·한 격자 진단 완료

`rb_thermal_ensemble_report_v3.json`: 고정6경로 모두 통과.
SHA256 `c5cfc2142020b38648ab8ccbbc7a68ad777690647057cdf33c6ee13848f2737f`.
Final maxdt =1.220703125e-10,6.103515625e-11,3.0517578125e-11s.
각 경로 N=ceil(τ/maxdt). τ·경로·물리 입력·허용오차 불변.

두 successive edge의 전6경로 최대2.3432344e-6/9.4625991e-8 <1e-3.
독립 비교 최대1.8223767e-8 <5e-6.
Independent 자체 refinement 최대3.4601050e-8 <2e-6.
Source62개·12개 원본 reference file hash 직접 일치 확인.

누적 실제 계산: primary30 + reference12 =42개 cache.
최종 finest6개만 workers4·각 BLAS1로 추가. 해당 실행627.447s.
별도 `rb_thermal_finest_jobs_v1.json`에 실행·source·spec 기록.
v3 집계의 primary18개 전부 cache hit; grid callback6개,cache hit30개.
집계 시간은 실제 ODE 계산 시간이 아님.

병렬 독립 가중합 검수: source별 두 ordering, Poisson number,
복소 response, total 등8항목 모두 bitwise 일치. 금지 phase 항0.
추가 nested grid의 동일6경로 raw8metric·rate 반감·digest 재결합도 검증.
원자 계산 반복 없이12/12 cache hit. Code/test 변경 없음.

출력: 한 p0/seed11 grid의 진단용 atomic stream.
6/6 경로 모두 합산, `spectra` 존재. `certified=false`.
Occupancy/nV=0.6176750449. 다른11개 grid/seed 조합 `NOT_EVALUATED`.
전체 thermal ensemble·nonlocal Maxwell·finite seed·absolute squeezing 미인증.
v1/v2 실패 report 보존. 새 ensemble 수렴 주장으로 소급 변경 금지.

공유 README·checklist 갱신. 정리 MD 한국어 caveman 적용.
합동 전체검사1487통과/기존삭제문서1실패 결과 유지.
