# PSR 모델 한계와 개선 구현 계획

사용자 요청, 2026-11-02: “PSR 레퍼런스 DB를 바탕으로 현재 PSR 스킴의 모델 한계와 개선점을 정밀 검토하고, 기존 원자·편광 계산을 최대한 재사용하는 방향으로 개선 우선순위와 구현 계획을 정리”.

## ⟦AI⟧ 검토와 제안

검토 기준: commit `387549130d5d2926759c667826a6c650a5b0fe95`의 PSR 구현과 [문헌 DB 27편](../../references/psr_experiments/README.md). 아래 우선순위·인터페이스·완료 기준은 구현 제안. 사용자 채택 또는 구현 완료를 뜻하지 않음. 이번 작업은 문헌·코드 검토와 제한된 읽기 전용 수치 확인.

### 1. 핵심 판단

현재 PSR은 **다준위 원자의 고전적 self-rotation 계산에 조건부 Gaussian 잡음 모형을 연결한 스킴**. 원자 구조·광펌핑·편광 규약·검출 손실·수치 진단은 재사용 가치가 높음. 반면 실제 squeezing/anti-squeezing의 주파수 스펙트럼을 예측하려면 **동일한 구동 원자 모형에서 편광 응답과 Langevin 잡음을 함께 계산**해야 함. 고전적 회전 곡선 또는 채널의 완전양성 검증만으로 실험적 squeezing을 결정할 수 없음. [Hsu 2006, pp.7–9](../../references/psr_experiments/papers/hsu_2006_pra.json), [Lezama 2008, pp.3–6](../../references/psr_experiments/papers/lezama_2008_pra.json).

과학적 최우선은 **유한 분석 주파수의 원자 잡음**. 다만 실험의 집속된 mW 파워 영역을 다루기 위한 구현 선행 조건은 **광학 포화를 포함한 다준위 정상상태와 명시적 dissipative operators**. 기존 해석·수치 코어를 재사용하고, PSR 원자–장 결합과 homodyne 검출 부분을 추가하는 구성이 적절함.

첫 개선 목표는 **Rb87 D1, 선형 입사 carrier + 직교 vacuum, B=0, 비코팅·무완충가스 셀의 단일 통과**. 이 제한된 영역에서 최소·최대 잡음, 최적 위상, 원자 잡음의 증가 영역을 함께 검증. 현재 조건부 모드는 비교 기준으로 유지. B 제어·seed·코팅 셀·공간 다중 모드는 후속 확장.

### 2. 현재 계산과 한계

[원자 응답](../../gabes/psr.py#L64)은 0, ±h, ±h/2의 QWP 설정에서 `multilevel_response`를 호출. Jones/Stokes 변화율로

\[
g=\partial_\epsilon(d\psi/dz),\quad
\alpha=-S_0^{-1}dS_0/dz,\quad \kappa=gL,\quad a=\alpha L
\]

계산. [잡음 모형](../../gabes/psr_noise.py#L20)은 정규화된 길이 좌표에서

\[
\frac{dV}{dt}=MV+VM^T+D,\quad
M=\begin{pmatrix}-a/2&\kappa\\0&-a/2\end{pmatrix},\quad
D=\frac{a(1+n_{\rm ex})}{2}I
\]

를 정확히 적분. 입력 vacuum은 `V=I/2`; 외부 검출 효율은 한 번 적용. 내보낸 `covariance`는 `2V`, 즉 SQL=I. `source_added_covariance` 등 채널 행렬은 canonical I/2 규약을 유지함. 이 두 규약의 변환은 보존해야 함.

| 현재 가정·구현 | 실제 한계와 영향 | 문헌 근거·개선 방향 |
|---|---|---|
| 정적 회전 기울기 → 일정한 shear | DC 광펌핑 응답이 MHz quadrature 응답과 같다는 보장 없음. 현재 analysis frequency 입력·RF 스펙트럼 없음 | Lezama pp.3–6: 주파수별 두 편광 응답과 잡음 계산. 주파수 슬라이더만 붙이는 방식으로 해결 불가 |
| `n_ex` 하나로 등방 잡음 가정 | fluorescence·spin noise·편광 상관·주파수 의존성 미계산. `n_ex=0`은 잡음의 생략이며 실제 squeezing의 엄밀한 상한도 아님 | Hsu pp.7–9; Mikhailov 2009 pp.3–5: squeezing 영역과 주변 excess-noise 영역을 함께 설명해야 함 |
| 약한 광학 여기에서 excited state 제거 | 바닥상태 광펌핑은 포함하지만 full optical saturation은 없음. 실험 파워를 UI에 넣는 것만으로 모델 범위 확장되지 않음 | Hsu p.8; Barreiro 2011 pp.2–3. 다준위 full-OBE 기준 계산 필요 |
| carrier의 α를 직교 vacuum 손실에도 적용 | 두 편광의 차등 흡수·분산과 conjugate coupling 미도출 | 전체 편광 fluctuation Jacobian과 같은 원자 상태의 잡음으로 교체 |
| 길이 방향 일정한 pump·계수 | local pumping·pump depletion·편광 변화에 따른 재계산 없음. 같은 OD의 밀도/길이 변경을 과도하게 동일시 | Rochester 2001 pp.6–7; Zhang 2017 pp.2–5: 평균장과 잡음 전파 동시 갱신 |
| 독립 Maxwell 속도군의 χ 평균 | 현재 고전 응답에는 적절한 근사이나, 잡음 평균까지 정해 주지는 않음 | Mikhailov 2009 Fig.5: 멀리 detune된 속도군도 저주파 잡음에 기여. 응답·잡음 각각 수렴 확인 |
| 한 개 유효 공간 모드·국소 intensity | 논문의 power/waist/focus를 직접 대응하기 어려움. LO와 aperture의 복소 모드 중첩 미계산 | Zhang 2013·2016: 단일 손실 η로 설명되지 않는 공간 효과 |
| 양의 carrier·throughput만 확인 | Stokes 선형화에 필요한 밝은 carrier 대비 직교 모드 fluctuation 크기의 정량 검증 없음 | beam area·주파수 대역을 정하고 carrier flux와 같은 관측 대역의 fluctuation을 비교해야 함 |
| B=0·직교 vacuum·등방 reset | 자기장, coherent seed, buffer gas·벽 충돌·원자 이동의 응답/잡음 없음 | Horrom·Lima·Zhao·Wang의 결과는 별도 확장 조건. 기존 인자 전달만으로 해당 물리 구현되지 않음 |

코드상 확인된 누락은 [현재 적용 범위](README.md)에 이미 대부분 명시돼 있음. 이번 검토에서 현재 선언한 모형과 모순되는 새 구현 오류는 확인되지 않음. 개선 대상은 주로 모형의 물리적 범위.

**수치적으로 안정된 조건부 값과 실험 정확도는 별개.** 현재 기본 선택점(Rb87 D1, 70°C, 75 mm, 0.1 mW/cm², reset 1000 s⁻¹, Δ=−900 MHz, 외부 η=0.9025, `n_ex=0`)의 직접 계산:

| 속도 적분 차수 / 실제 속도군 수 | κ | OD | 최소 / 최대 잡음 dB |
|---|---:|---:|---:|
| 4 / 116 | 0.4231674 | 0.2857264 | −1.3782537 / +1.4792840 |
| 6 / 174 | 0.4231621 | 0.2857261 | −1.3782376 / +1.4792654 |
| 10 / 290 | 0.4231620 | 0.2857261 | −1.3782375 / +1.4792652 |

이 한 점의 수치 일치는 확인됨. 차수 6의 virtual-excitation bound≈0.01114, ellipticity drift×L≈−0.00878. 같은 조건에서 intensity 1·10 mW/cm²로 올리면 bound≈0.1114·1.1139로 현재 0.05 screen을 초과. 수치 정상상태 검사는 여전히 통과하므로 두 상태를 구분해야 함. 이 bound는 보수적 engineering screen이며 실제 오차율이 아님.

기존 `|OD|≤0.3`, `|drift×L|≤0.1`, derivative 1% 등의 screen을 일반적인 물리 경계로 승격하지 않음. 또한 Quick/Balanced/Fine 광학 격자 간격은 75/37.5/18.75 MHz. 직접 선택점 계산은 선택점 오류를 막지만 격자 사이의 좁은 구조까지 보증하지 않음.

### 3. 재사용 범위와 필요한 연결부

| 기존 코드 | 재사용할 부분 | 제한·추가 작업 |
|---|---|---|
| [`species.level_structure`](../../gabes/species.py#L338), `angular` | 준위 에너지, signed dipole, Γ, isotope mass, angular coefficients | PSR full-OBE의 원자 데이터 원본. F-lumped SAS용 `build_manifold`로 대체하지 않음 |
| [`magneto_multilevel`](../../gabes/magneto_multilevel.py#L65) | 약한 여기 기준 계산, Jones 규약, 밀도·공명 중심 속도 적분 | 현재 반환 χ·진단만으로 원자 확산을 복원할 수 없음. 내부 상태·generator·jump operators를 별도 명시적 인터페이스로 공급 |
| [`AtomModel`](../../gabes/atoms.py#L23), [`core`](../../gabes/core.py#L580) | row-major Liouvillian, 정상상태, trace-zero 선형 응답, 선형대수 | 새로운 독립 OBE solver를 만들 필요 없음. optical detuning·Hamiltonian은 rad/s, reset은 s⁻¹ 규약 유지 |
| [`ExplicitReservoirs`](../../gabes/quantum/reservoirs.py#L40) | 실제 jump operators와 generator 일치 검사, 명시적 thermal reset | 방출·reset·dephasing의 물리 출처/단위/종류 구분. 단순 감쇠율 추가 뒤 잡음 생략 금지 |
| [`stationary_atomic_noise`](../../gabes/quantum/diffusion.py#L176) | jump-derived ordered diffusion, QRT·2차 모멘트 검사 | 작은 기준 계산부터 사용. 모든 jump별 큰 행렬 저장 방식은 열평균 적용 전 최적화 필요 |
| [`eliminate_atomic_noise`](../../gabes/quantum/traveling.py#L106), `sum_independent_classes` | 원자 fluctuation→장 drift, greater/lesser noise, 독립 속도군 합산 | PSR의 편광·±Ω 모드와 photon-flux 결합 정규화는 새로 정의 |
| [`constant_segment`](../../gabes/quantum/traveling.py#L219), `compose_segments`, [`GaussianChannel`](../../gabes/quantum/channels.py#L59) | 국소 전파, 분산 잡음 적분·합성, 검출 손실 | segment 도구 자체는 pump depletion을 풀지 않음. PSR 쪽에서 평균장을 전파하고 국소 상태 재계산 |
| [`sidebands`](../../gabes/quantum/sidebands.py), [`readout`](../../gabes/quantum/readout.py) | RF 대역·채널 평균의 규약과 해석 테스트 | 현재 주요 builder는 probe/conjugate 4측파대 IDS 전용. PSR용 homodyne/Stokes 측파대 adapter 필요 |
| [`probe_saturation.gaussian_peak_intensity`](../../gabes/probe_saturation.py#L34), [`spatial_modes`](../../gabes/quantum/spatial_modes.py) | 명시적 Gaussian 단위 환산, reciprocal mode-projection 원칙 | SAS 흡수 saturation law를 PSR 양자 잡음에 적용하지 않음. 현재 spatial builder의 FWM 모드 수·이름도 그대로 사용 불가 |
| [`PSRScheme`](../../gabes/schemes/psr.py#L52) | scan axis, 0 dB, 두 orthogonal phase traces, 직접 선택점·cache·export | 새 backend와 RF 관측량 추가 시에도 기존 상태·실패 gap·selected-point 일치 보존 |

재사용 시 특히 주의할 다섯 항목:

1. **QWP 미분과 횡자기장.** 현재 QWP 변화는 β≈0에서 ψ≈β, ε≈−β를 함께 변화시킴. B=0 회전 대칭에서는 현재 slope에 문제가 없지만, 횡자기장·비등방성이 있으면 두 미분이 섞임. 독립적인 Jones real/imag perturbation과 전체 편광 Jacobian을 먼저 정의해야 함. [`_polarization`](../../gabes/magneto_multilevel.py#L111), [`_stokes_rates`](../../gabes/psr.py#L37).
2. **복소 주파수 응답.** [`GroundCoordinates`](../../gabes/magneto_coordinates.py#L28)의 Hermitian 정적 상태용 pack/projection에 복소 Fourier 응답을 그대로 넣으면 필요한 성분을 잃음. basis는 재사용하되 복소 응답 좌표를 유지.
3. **reservoir 출처.** `AtomModel.with_reset`의 reset을 `decay`로 가져오면 `ExplicitReservoirs.from_atom`이 방출로 분류할 수 있음. 수학적 generator 일치와 별도로 방출·transit noise의 이름/출처를 정확히 구성. [`thermal_reset_channels`](../../gabes/quantum/reservoirs.py#L114)를 별도로 합성하는 방향 권장.
4. **다른 작업의 미추적 코드.** 로컬 `gabes/fwm_zeeman_response.py`는 현재 미추적이며 Rb85 D1·FWM 주파수·π-pump·기존 Rabi 정규화에 특화. 구현 패턴 참고는 가능하지만 PSR의 새 의존성으로 삼거나 Claude 작업 중 이동·리팩터링하지 않음. 공유된 준위 데이터와 코어에서 구성.
5. **약한 여기/full-OBE 비교.** 현재 weak backend는 ground-F coherence를 생략함. 두 계산의 ground-F secularization, reservoir/reset, detuning·dipole/Rabi 규약을 맞춘 영역에서 비교. full backend가 추가 coherence를 유지한다면 그 근사 차이를 solver 오류와 구별하고 무조건적인 수치 일치를 요구하지 않음.

### 4. 새 계산 경로와 규약

아래 API 이름은 제안이며 아직 존재하지 않음. 최초 구현은 PSR 전용 모듈을 추가해 기존 공유 엔진의 기본 동작을 보존하는 방향.

```text
prepare_psr_atom(atomic_inputs)
  → 준위·Hamiltonian·명시적 reservoir·편광 basis·mode normalization
solve_psr_pump(prepared, detuning, local_jones_field, velocity_nodes)
  → 각 속도군 rho_ss·mean polarization·excited populations·진단
psr_local_generator(pump_state, rf_axis, mode_spec)
  → 편광/측파대 drift와 두 순서의 Langevin noise
propagate_psr(local_states, geometry)
  → 평균 carrier와 주파수별 장 전달·잡음
psr_homodyne_readout(transfer, carrier, detector, phase_axis, rf_filter)
  → SQL 기준 phase/min/max noise·검출 기준면·대역 정보
```

주파수 응답은 예를 들어 `e^(−iΩt)` 규약에서

\[
(-i\Omega I-A_v)\,\delta x_v=B_v\,\delta a+F_v
\]

로 구성. `A_v`, coupling, diffusion은 **동일한** 구동 정상상태·reservoir에서 유도. 이것은 설계식이며 고전 χ에 임의의 잡음을 붙이는 처방이 아님. optical Δ[MHz]와 analysis f[Hz]를 분리하고 `Ω=2πf`를 한 번 적용. 정적 reset의 s⁻¹에는 2π를 적용하지 않음.

- **상관·모드:** 두 광학 편광과 ±Ω, normal/anomalous correlation을 유지. 선형 입사·B=0 대칭으로 직교 모드 쌍만 남길 수 있는지는 full reference와 비교해 확인. FWM의 beat offset·mode label을 이식하지 않음.
- **속도·밀도:** coherent response와 독립 원자의 ordered noise kernel을 각각 평균한 뒤 같은 장을 전파. 각 class에 밀도/확률을 한 번 적용. noise에 Maxwell 확률의 제곱을 쓰거나, 최종 dB·각 class의 전체 출력 채널을 단순 평균하지 않음. [Lezama Eqs.21–29](https://arxiv.org/abs/0710.0609), [Mikhailov Fig.5](https://arxiv.org/abs/0903.3156).
- **손실:** 도출된 편광 drift와 확산에 들어간 매질 손실은 별도 vacuum-loss channel로 다시 추가하지 않음. 외부 광학 손실·QE는 현재처럼 뒤에서 한 번 적용.
- **전파:** 국소 Jones 평균장으로 pump state를 다시 풀고 동일한 상태에서 response/noise 산출. 평균장 depletion을 무시하는 첫 backend에는 그 범위를 유지. Lezama의 유한주파수 모형도 원문 p.4에서 평균 pump depletion을 생략함.
- **검출 대역:** RBW 숫자 자체를 임의 mode duration으로 취급하지 않음. RF transfer/filter, PSD의 단면/양면 규약과 SQL 기준을 먼저 정의. 유한 대역 출력은 PSD/선형 variance에서 적분한 뒤 dB 변환.
- **공간 검출:** Gaussian power 평균은 LO 투영과 다름. `b_LO = ∫u_LO*(r)b(r)d²r`에 대응하는 공분산·교차상관으로 검출 모드를 계산. local squeezing dB의 반경 평균으로 대체하지 않음.
- **기술잡음:** 검출 전자잡음·LO 검출 위상 흔들림·검출기 imbalance에 의한 공통 모드 누설은 downstream에서 처리. 반면 입사 pump/seed의 RIN·위상/주파수·편광 요동은 원자–장 응답을 통과하는 입력 기술잡음. 입력 covariance/forcing으로 전파하며 임의의 출력 가산 잡음으로 대체하지 않음. 검출 위상의 Gaussian jitter는 covariance를 위상 분포로 평균하는 방식부터 가능. 전자잡음 차감 전후 데이터와 검출기 포화는 별도 조건.

현재 model ID `effective_psr_uniform_shear_distributed_loss_v1` 유지. full-OBE 평균장만 추가한 중간 단계도 기존 가정된 diffusion을 사용한다면 여전히 조건부 모형으로 표시. 미시적 잡음 backend에 별도 version·가정·입출력 metadata 부여.

### 5. 개선 우선순위와 구현 단위

과학적 영향 순위: **유한주파수 원자 잡음 → 포화·편광별 전파 → 공간/LO → 실험 제어·수송**. 아래는 실제 의존성을 반영한 구현 순서. 새 승인 절차가 아니라 작업 분할과 검증 경계의 제안.

| 단계 | 구현 범위·산출물 | 재사용·의존성 | 단계 완료 기준 | 규모 |
|---|---|---|---|---|
| S0 비교·입출력 규약 | 기존 conditional 결과 고정. optical 기준점, power/beam 정의, 분석 주파수·대역, detector plane, input seed 계약. 문헌 benchmark 레코드 | PSR export와 DB의 quantity/source/basis 필드 | 서로 다른 운전점·보정값을 합치지 않음. 단위·부호·SQL·phase 규약 테스트. 미보고는 미보고로 유지 | 작음 |
| S1 다준위 정상상태 기준 | Rb87 D1·B=0·단일 속도 full-OBE. 두 excited hyperfine manifold 유지. 명시적 방출·reset, 독립 Jones perturbation | `level_structure`, `core`, `ExplicitReservoirs`. S0 필요 | trace/Hermiticity/positivity/residual, dipole branching, secular/reset 규약을 맞춘 약한 여기 overlap 확인. 파워에 따른 정상상태 포화 확인 | 중간 |
| S2 유한주파수 잡음·검출 | 같은 상태에서 diffusion·복소 response, ±Ω 편광 coupling, uniform thin-cell 전파, PSR homodyne adapter | `diffusion`, `traveling`, `channels`. S1 필요 | QRT/확산 identity, commutator·CP·uncertainty, vacuum·passive loss·ideal shear, ±Ω/위상 규약 확인. 주파수에 따른 min/max/최적 위상 출력 | 큼 |
| S3 열평균·계산 효율 | drift와 noise를 모두 속도 적분. 공명·잡음 구조에 따른 refinement. 직접 기준 계산과 최적화 비교 | 기존 velocity quadrature와 S2 | 회전뿐 아니라 min/max·noise kernel의 속도 수렴. 잘린 Maxwell tail 영향과 메모리·시간 측정. thermal scan 전 성능 기준 수립 | 큼 |
| S4 평균장·잡음의 길이 전파 | pump power·Jones 편광을 z 방향 갱신. 각 구간 상태·편광별 손실·잡음 재계산 | `constant_segment`, `compose_segments`. S1–S3 | thin/uniform 극한 복원, z 분할 수렴, 매질 vacuum noise 중복 없음, 동일 OD/다른 밀도·길이 비교 | 큼 |
| S5 측정·첫 공간 보정 | 실제 power/waist/focus 입력, 선언된 Gaussian geometry, LO overlap. 검출 위상/전자잡음과 매질을 통과하는 입력 기술잡음 분리 | 검출부는 S0–S2와 병행 가능. 입력 잡음은 S2 응답, 공간 source는 S4 필요 | source/detector plane 일치, phase ellipse와 반대 quadrature 함께 비교. 광학 손실과 검출기 포화를 구분 | 중간 |
| S6 실험별 확장 | 필요에 따라 횡/종 B, orthogonal seed, D2·Rb85, 다중 공간 모드·다중 통과, buffer gas·coating·sensor cells | 해당 물리의 정상상태·reservoir·관측량을 각각 확장 | 각 확장에 해당하는 문헌과 음성 결과 검증. 한 기능의 성공을 다른 regime 검증으로 전용하지 않음 | 기능별 큼 |

**첫 구현 묶음 권장: S0 + S1.** 실험에 연결할 단위·편광 계약과 포화 가능한 원자 상태를 먼저 확보. 이어 S2에서 단일 속도·균일 pump의 제한된 finite-RF spectrum을 완성. 전체 열분포·긴 셀·모든 isotope를 처음부터 한 번에 확대하지 않음.

낮은 비용의 병행 개선: 기존 광학 scan의 adaptive sampling, 직접 선택점·경고 정보 보존, 단위가 명확한 power/waist 환산, 검출 손실·위상 jitter. 이는 사용성과 비교 가능성을 높이지만 미시적 잡음 계산의 대체물이 아님.

실험적 근사를 빠르게 활용할 필요가 있으면 measured covariance/PSD에 근거한 주파수별 추가 잡음 모드를 별도로 둘 수 있음. fitting 데이터·공식·측정 기준면·외삽 조건·예측 한계 명시. 최소 squeezing만 맞추거나 `n_ex` 하나에 saturation·RIN·mode mismatch를 모두 흡수하지 않음. 실험기반 보정의 성공과 원자 잡음 유도의 성공을 구분.

### 6. 성능과 cache 설계

현재는 Rb87/Rb85의 retained ground density가 34/74개 좌표이며, B=0 계산은 더 작은 even-Δm sector를 활용. full-OBE는 다른 규모:

| 원자 | 전체 Hilbert 차원 | 전체 density 좌표 수 |
|---|---:|---:|
| Rb87 D1 | 16 | 256 |
| Rb85 D1 | 24 | 576 |
| Rb87 D2 | 24 | 576 |
| Rb85 D2 | 36 | 1296 |

현재 `stationary_atomic_noise`의 `channel_diffusion[jump,operator,operator]`를 그대로 확장하면 **thermal-reset jump 저장만** complex128 기준 Rb87 D1 약 0.124 GiB, Rb85 D2 약 10.8 GiB/속도군. 실제 차원으로 계산한 allocation 추정이며 실행시간 측정치는 아님. 기본 선택점 차수 6도 174개 속도군이므로 전체 배열을 보관하는 설계는 부적절.

기준 계산 후 검토할 최적화: 물리적 reservoir별 diffusion 합산·streaming, 필요한 optical observables로의 QRT/resolvent 투영, 대칭으로 닫힌 coherence sector. jump provenance는 유지. 일반 ellipticity·횡 B에서 block closure가 깨질 수 있으므로 full reference와 비교하고 fallback 마련. 평균 응답이 빠르다는 이유로 잡음 계산도 빠르다고 가정하지 않음.

cache는 세 층으로 분리:

- **원자 평균장:** isotope/line, Δ, pump Jones field·intensity, B, 온도·속도군, reset/충돌, geometry와 local pump 상태.
- **스펙트럼:** 평균장 cache + RF axis·field coupling·noise model. RF 변경은 정상상태를 재사용하지만 새로운 resolvent/noise 계산 필요.
- **검출:** θ, QE, 외부 손실, 검출 LO 위상 jitter·전자잡음·기존 RF 구간 안의 filter. 가능한 한 원자 재계산 없음. 작은 입력 기술잡음은 같은 평균장의 선형 전달을 재사용하되 입력 covariance/forcing에서 반영. 광학 seed나 큰 pump 요동이 원자 평균장을 바꾸면 첫 층부터 재평가.

최적화 전 목표 격자의 시간·peak memory를 측정해 UI의 Quick/Balanced/Fine 범위 설정. 이번 검토는 전체 full-OBE 성능을 측정하지 않았으므로 초 단위 완료시간이나 현재 81점 속도의 유지 약속은 하지 않음.

### 7. 문헌 비교와 회귀 검증

| 비교 단계 | 우선 문헌·원문 위치 | 평가할 항목·주의점 |
|---|---|---|
| 이상적·수학적 기준 | [Matsko 2002](../../references/psr_experiments/papers/matsko_2002_pra.json), pp.2–4; Lezama pp.3–6·Appendix | 무손실 shear의 covariance와 주파수 응답·확산 identity. Matsko Eq.29의 amplitude ratio와 특수 dB 표기는 통상 noise-power dB와 직접 비교하지 않음 |
| 첫 thermal 경향 | [Mikhailov 2009](../../references/psr_experiments/papers/mikhailov_2009_jmo.json), Fig.3–5 | 1.2 MHz density/power map의 squeezing 영역과 주변 excess noise, 두 excited hyperfine의 영향. 현재 DB에는 digitized map이 없으므로 즉시 점별 오차 점수를 주장하지 않음 |
| 첫 phase/Stokes 경향 | [Barreiro 2011](../../references/psr_experiments/papers/barreiro_2011_pra.json), pp.2–3·Figs.2–3; [Valente 2015](../../references/psr_experiments/papers/valente_2015_pra.json), Eqs.3–12 | detuning별 최적 phase, covariance/양 quadrature. Barreiro: 30 mW·73°C·50 mm·2.7 MHz에서 측정 2.9±0.1 dB, 손실 보정 4.7±0.2 dB. 200 µm waist는 intensity e⁻¹ 정의. >40 mW 감소는 검출 증폭기 포화 가능성도 언급됨 |
| 필수 음성 결과 | [Hsu 2006](../../references/psr_experiments/papers/hsu_2006_pra.json), pp.6–9; [Agha 2010](../../references/psr_experiments/papers/agha_2010_oe.json), p.3; [Horrom 2011 냉각 원자](../../references/psr_experiments/papers/horrom_2011_jmo.json), p.7 | strong rotation·D2·cold atoms가 자동으로 squeezing을 뜻하지 않음. Hsu 결과도 모든 PSR 영역의 불가능성 정리가 아님 |
| 확장된 주파수·매질 | [Ries 2003](../../references/psr_experiments/papers/ries_2003_pra.json), pp.3–4; [Mikhailov 2008](../../references/psr_experiments/papers/mikhailov_2008_ol.json), pp.2–3 | Ries D2: 5 MHz에서 측정 0.85±0.05 dB와 보정 1.23±0.07 dB 분리. Mikhailov 저주파 결과는 Ne buffer gas이므로 최초 무완충가스 절대값 비교에서 분리 |
| 공간·전파 | [Zhang 2013](../../references/psr_experiments/papers/zhang_2013_ol.json), [2016](../../references/psr_experiments/papers/zhang_2016_pra.json), [2017](../../references/psr_experiments/papers/zhang_2017_pra.json); [Lanning 2017](../../references/psr_experiments/papers/lanning_2017_pra.json) | aperture/LO phase·다중 통과·밀도/길이 비등가성. Lanning은 고전적 Born 근사, Zhang의 illustrative mode model은 보정·원자 잡음 생략 포함. 완성된 noise backend로 취급하지 않음 |
| 제어·응용 | [Zhao 2025](../../references/psr_experiments/papers/zhao_2025_sciadv.json), [Lima 2024](../../references/psr_experiments/papers/lima_2024_pra.json), [Novikova 2015](../../references/psr_experiments/papers/novikova_2015_pra.json), [Monsa 2026](../../references/psr_experiments/papers/monsa_2026_prapplied.json), [Wang 2026](../../references/psr_experiments/papers/wang_2026_lsa.json) | B·seed·기술잡음·센싱·코팅 수송 분리. 최저 noise phase와 최적 sensing SNR 불일치 가능. Monsa source/sensor 셀 분리, Wang은 squeezing-level 동기화 결과이며 직접 dynamic phase locking 증거는 아님 |

검증 결과는 네 층으로 구분:

1. **기존 회귀:** [`test_psr.py`](../../tests/test_psr.py), [`test_psr_noise.py`](../../tests/test_psr_noise.py), [`test_psr_scheme.py`](../../tests/test_psr_scheme.py)의 부호·밀도/길이 scaling·실패 gap·loss·selected-point cache 보존. 기존 모드와 새 모드의 물리적 예측이 항상 같을 필요는 없음.
2. **새 수학·수치 검증:** 약한 여기 overlap, trace-zero response/finite difference, QRT 일치, reservoir noise, ±Ω 공액 관계, spectral/검출 SQL 정규화, vacuum·loss·ideal shear. 대칭 축 회전, velocity·z·RF 적분 수렴과 직접/최적화 결과 비교.
3. **문헌 경향·절대값:** 먼저 phase·frequency·power/density 추세와 음성 결과. 단위·검출 기준면·입력이 충분한 경우에만 절대값 오차 계산. missing beam/η/reset을 조용히 조정해 가장 좋은 dB에 맞추지 않음.
4. **실험적 fitting·외삽:** 수행한다면 apparatus group 단위 분할, 학습 파라미터·범위·공통 측정 장치 명시. 같은 scan의 행을 무작위로 나누어 독립 검증으로 부르지 않음. 이미 열람한 27편 DB 자체는 미열람 blind validation이 아님.

공유 코어를 수정하지 않는 adapter 단계에는 관련 PSR·원자·quantum primitive 테스트를 선택. 공유 solver·reservoir·normalization을 변경하면 그 호출자와 교차 스킴 회귀까지 확대. 기존 Grand Challenge의 봉인된 캠페인을 이 계획 때문에 재실행하지 않으며, 그 프로젝트의 no-fit·완전성 조건을 PSR 구현의 별도 승인 요건으로 만들지 않음.

### 8. 이번 검토의 확인 범위와 다음 작업

27편 DB record를 대조하고 핵심 이론·실험의 보관 원문을 재확인. Lezama p.4 및 Barreiro p.3 수식·그림을 직접 렌더링해 확인. 코드·문헌·재사용 구조를 독립적으로 나누어 검토했고, 기본 선택점의 적분 차수/강도 확인은 위 표의 제한된 조건에서 수행.

이번 결과는 문서·계획. 생산 solver/UI 변경, 실험 fitting, 전체 수치 회귀, 새로운 미시적 squeezing 예측은 수행하지 않음. 다음 구현 단위는 **S0 + S1: 측정·편광 규약 고정과 Rb87 D1 단일 속도 full-OBE 기준 backend**. 이후 S2의 유한주파수 response/diffusion에 같은 상태·reservoir를 연결하는 방향 권장.
