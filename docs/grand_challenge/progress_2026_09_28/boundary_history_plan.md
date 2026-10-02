# Thermal 이후 경계 이력 검증 계획 — 2026-09-28

**설계 문서. 새 물리 계산·수렴 인증 없음.** 현재 thermal 캠페인 완료·수렴 여부와 독립. 다음 반증 대상: “고정된 collection mode 안의 원자 stream에서, 가상 옆 경계의 즉시 비편극 재설정이 5% 이내 영향인가?”

09-18 정지 원자 probe는 옆면 중점의 고정 Hamiltonian 아래 입사 전 pump 노출에 따른 상태 변화만 확인. 이동 경로·열적 평균·source noise·실제 셀 경계 검증 아님. 원본 [probe](../../../analysis/grand_challenge/checkpoint_2026_09_18/physics_probe.txt), [점검](../checkpoint_2026_09_18.md) 보존.

## 1. 고정할 대상과 현재 코드의 함정

| 항목 | 고정값·주의 |
|---|---|
| Collection column Q | x,y ∈ [−a,a], a=√(1.2×10⁻⁷)/2 m=173.205 μm; z ∈ [−6.25,6.25] mm |
| Photon-flux 면적 | A=1.2×10⁻⁷ m² 유지. 물리 영역 확대와 함께 A를 늘리면 coupling도 달라짐 |
| Pump | 0.6 W, 1/e² intensity radius w=530 μm; amplitude exp[−(x²+y²)/w²] |
| 원자·환경 | 동일 reduced RMS ⁸⁵Rb, T=373 K, 독립 고정 n=10¹⁸ m⁻³, radiative jumps; transit reset 0 |
| Detuning·geometry | Δ=2π×0.9 GHz, δ=−2π×8 MHz, probe/conjugate angle 0.006/−0.005 rad |
| RF | 동일 laboratory 0.1/1/4 MHz; optical detuning과 구분 |
| 바꿀 항목 | 원자가 Q에 도착하기 전 pump 이력의 시작 위치·그에 따른 조건부 입사 상태 |

`reduced_ballistic_problem`의 readout/drive는 현재 경로 전체에서 일정한 transverse projection. `smooth_rb_problem`은 pump envelope만 연속화. **현재 box를 넓혀 같은 함수를 그대로 호출하면 collection support까지 넓어짐.** 경계 이력 하나의 효과로 해석 불가.

`ThermalRbModel.boundary_state`도 모든 경로에 같은 단일 4×4 상태. 새 이력 상태는 위치·속도별로 다르므로 이 필드 하나 교체로 thermal 모델 구현 불가. Solver 수정 전 독립 진단 runner에서 경로별 상태·이력 hash를 명시하고 새 결과 경로 사용. 이력이 달라진 계산에 기존 thermal cache·certificate 전용 금지. 정확히 같은 physical request의 재사용은 원본 hash·동일성 근거를 따로 보존. 구현 시 모델 identity에 이력 규칙·물리 경계·collection support 포함.

## 2. 첫 실험: 같은 Q 경로에 입사 전 구간 연결

물리 계산 영역을 B_R=[−R,R]²×[−L/2,L/2], L=12.5 mm로 정의. R=a, w, 2w, 3w 순서. Q·A·pump·n·T 고정. R=a는 기존 계산 복원 control.

Q의 각 입사점 rᵢ와 속도 v에 대해 **역방향 직선이 B_R 경계에 처음 닿는 시간** τ를 계산:

\[
\tau=\min_k\begin{cases}
(r_{i,k}-\ell_k)/v_k,&v_k>0,\\
(h_k-r_{i,k})/(-v_k),&v_k<0,\\
\infty,&v_k=0,
\end{cases}\quad
r_b=r_i-v\tau.
\]

입사 전 r(t)=r_b+vt, 0≤t≤τ. 바깥 경계의 가정 상태 ρ_b=diag(5/12,7/12,0,0)를 같은 pump-only Liouvillian으로 전파해 ρᵢ 산출. 이어 기존 Q chord만 측정. τ=0이면 ODE 생략·ρᵢ=ρ_b. Convex Q·충돌 없는 직선에서는 방문 한 번; 재입사 모델 없음.

**실제 축 방향 끝 고정.** z=±6.25 mm를 넘어 pump history 연장 금지. 역추적이 먼저 셀 끝에 닿으면 거기서 종료; axial 입사 경로는 τ=0일 수 있음. 이 L은 현재 conditional fixture의 셀 길이이며 실제 window·wall의 상태 재설정은 여전히 미검증. 알려진 물리 측벽이 R보다 안쪽이면 그 벽에서 종료하고 wall model 필요. 실제 셀 반경·coating·충돌 자료 없이 3w 영역을 실측 셀로 명명하지 않음.

R=3w 옆면 중점 intensity fraction exp(−18)≈1.52×10⁻⁸. 이것만으로 이력 오차 상계 증명 불가. 느린 원자의 긴 노출, axial 입사, 경로별 ∫|Ω_p(t)|²dt와 τ 분포도 기록. R=2w→3w 변화가 크면 더 큰 R 또는 wall model로 후속; cutoff 성공 선언 금지.

## 3. Hamiltonian·위상 일치

전체 pre-entry와 Q 구간에서 같은 **static pump frame** 사용. [현재 Hamiltonian](../../../gabes/schemes/fwm.py)의 순서 (g₁,g₂,e₂,e₃), angular-frequency 단위:

\[
H(t)=\operatorname{diag}(0,\omega_{HF},-\omega_{eHF}-\Delta_v,-\Delta_v)
 +e^{-r_\perp(t)^2/w^2}H_{pump},\quad
\Delta_v=\Delta-k_0\cdot v.
\]

H_pump는 `reduced_dipoles('uniform-zeeman-rms')`와 동일 power-to-Rabi·transition scales 사용. Radiative dissipator 한 번만 적용. Finite-seed Floquet frame 또는 δ를 ground diagonal로 쓰는 frame으로 중간 교체 금지.

같은 전역 pump gauge에서 연속 전파하면 ρᵢ는 그대로 이어받음. 다른 gauge를 쓰는 구현은 접점에서 두 frame의 관계 W=U_inner U_pre†를 명시하여 ρᵢ=Wρ_preW†, 모든 readout·drive·source covariance도 동일 변환. Population만 넘기면 pump-induced coherence 손실.

Beat b=−ω_HF+δ, s=(1,−1,−1,1), q=(k_p−k₀,k_c−k₀,−k_p+k₀,−k_c+k₀). Pre-entry 시작 laboratory beat phase φ_b를 사용하면 Q 입사에서

\[
\phi_i=\phi_b+b\tau,\qquad
\theta_j=s_j\phi_i-q_j\cdot r_i,\qquad
\nu_j=s_jb-q_j\cdot v,\quad w_j=\Omega+\nu_j.
\]

이는 θᵢ=θ_b+ντ와 일치. ν_p+ν_c=−(k_p+k_c−2k₀)·v 보존; 임의 momentum closure 금지. 두 구간에서 entry phase를 독립 난수로 뽑지 않음.

접점 시간 원점을 바꿀 때 raw pulse에도 Fourier 위상 필요. 고정된 같은 수치 readout/drive를 쓰는 일반 좌표에서 drive의 각주파수를 ξ_k로 별도 정의; V_k(t)=V_k exp(−iξ_kt). 따라서 Y_global=D_wY_local, C_global=D_wC_localD_w†, R_global=D_wR_localD_ξ†; D_w=diag exp(iw_jτ), D_ξ=diag exp(iξ_kτ). Response 적분의 exp(iw_jt−iξ_kt′)에 t=τ+a, t′=τ+a′를 넣으면 exp[i(w_j−ξ_k)τ]가 나오는 정확한 좌표 변환. ξ_k는 속도 성분 v_k와 다름. Readout 자체에 θᵢ를 이미 넣었다면 ντ를 재곱하지 않음. 해당 경우 pulse의 남는 공통 exp(iΩτ)만 처리. 직접 전 경로 적분과 대조해 이중 위상 적용 방지.

Pump-only state는 공통 beat phase와 무관. 따라서 기존 M_jk=δ(s_j,s_k) 평균 계속 가능. 유한 seed·위상 의존 경계 상태 추가 시 이 전제부터 재검증.

## 4. 상태 변경과 전체 noise history 구분

첫 산출물: 경로별 Δρᵢ, exit state, mean pulse, **합계** greater/lesser, mean outer, 복소 retarded response의 변화. 명칭은 `conditional pre-entry state-history sensitivity`.

Pre-entry에서 collection·weak drive가 정확히 0이고 원자가 Markov GKSL을 따른다는 선언 아래, Q 입사 시 **완전한 ρᵢ**는 이후 Q 안의 전체 원자 상관 계산에 충분. 그러나 현재 solver는 C(ρᵢ)를 모두 `atomic_inflow`에 배정하므로 **바깥 boundary와 pre-entry radiative jump의 기원별 noise는 복원하지 못함**. 바뀐 `atomic_inflow`를 기존 외부 boundary source와 동일 물리량으로 비교하지 않음. 코드: `smooth_transport.py` 135–136행에서 `_moments(rho0,f)`의 `cin`을 첫 source block에만 저장.

실제 Hermitian basis F에서 `_moments`는 C_ab=Tr(ρF_aF_b)−⟨F_a⟩⟨F_b⟩. 반대 ordering은 C_ba, 즉 **Cᵀ**. C†로 바꾸면 ordering 교환 사라짐. 현재 solver도 `transpose_columns` 및 `samples.swapaxes(-1,-2)` 사용. 출발·접점·종료에서 같은 규약 유지.

완전한 source ledger를 유지하려면 바깥 경계에서 C_b(0)=C(ρ_b), jump별 C_r(0)=0으로 시작:

\[
\dot C_b=AC_b+C_bA^T,\qquad
\dot C_r=AC_r+C_rA^T+D_r(\rho),\qquad
D_{r,ab}=\operatorname{Tr}\rho[L_r^\dagger,F_a][F_b,L_r].
\]

Q 입사에서 각 C_r(τ)를 넘김. 합 Σ_r C_r(τ)=C(ρᵢ) 검사. `outer_inflow`, `pre_entry:jump:*`, `in_collection:jump:*`를 구분한 ledger로 reference 대조. Current `smooth_wavepacket`은 boundary density만 받아 C₀=C(ρᵢ), 나머지 source=0으로 초기화하므로 이 전달에는 새 명시적 API·contract 필요.

Collection/readout 또는 weak drive가 pre-entry에도 존재하면 추가로 누적 mean·response·atomic–readout cross blocks·readout covariance를 접점에 전달. 두 구간 covariance를 별개로 계산 후 더하면 cross terms 손실. Q 밖의 유한 Gaussian collection, 재입사·충돌·벽 memory까지 검증했다고 주장 금지.

## 5. Flux·density 회계

첫 실험은 **기존 Q inflow quadrature와 동일 속도·rate** 유지. 이력에 따라 원자 수·flux를 재가중하지 않음. Face별 N=2ᵖ, σ=√(k_BT/m)에서

\[
\lambda_j=\frac{nA_{face}\sigma}{\sqrt{2\pi}N},\qquad
S^{\gtrless}=\sum_j\lambda_j M\circ
\left[C_j^{\gtrless}+\mu_j\mu_j^\dagger\right].
\]

λ에 density 포함. 추가 n·ground population·nV/계산 occupancy·outer/inner area ratio 곱셈 없음. Poisson mean outer 한 번 포함. μ의 phase average가 0이어도 μμ† 삭제 금지. 이력 없는 외부 경로를 Q 입사 경로에 다시 합산하지 않음.

후속 독립 geometry control: B_R의 실제 incoming flux를 적분하고 fixed Q에 도달하지 않는 경로의 pulse를 0으로 유지. Q 도달 경로만 남겨 rate 재정규화 금지. Outer-flux 결과와 inner-flux 역추적 결과의 합계 stream 비교. Collection support가 Q라는 동일 정의 아래 수렴 시 일치해야 함.

## 6. 사전 판정·독립 검산

1. 기존 thermal p2→p3→p4·세 seed gate 판정 먼저 보존. 미수렴이면 이력 테스트는 paired 경로 진단만 보고; ensemble 경계 오차 인증 보류.
2. 대표 경로는 결과 보기 전 선택: 각 입사면, ±v_z, 긴/짧은 τ, 강/약 pre-entry exposure 포함. State는 full-density DOP853와 별도 refined midpoint exponential 또는 CF4로 대조. Trace·Hermiticity·양의성 잔차 보존; clipping 금지.
3. 접점 누락 검출 control: R=a 복원, pump=0, τ=0, arbitrary common entry phase, 다른 시간 원점. Constant-H exact expm과 전 경로 QRT를 사용. Nonzero mismatch·coherence가 있는 경로 필수.
4. Fixed-Q split 진단과 전 경로 적분 비교. 후자는 Q 밖 readout/drive=0인 정확한 window를 사용하며 접점에서 적분 분할. 현재 adjoint/QRT도 fixed readout만 받으므로 window 전달을 명시적으로 구현한 별도 reference 필요. Source 분해는 backward observable·Lindblad-product reference, 총 ordered covariance는 raw full-density QRT로 확인.
5. 기존 path 기준 유지: primary 연속 두 refinement <10⁻³, independent reference 자체 <2×10⁻⁶, primary-reference <5×10⁻⁶. RF별·source별 complex Frobenius 비교와 기존 SI dark floor 사용. Gate 완화 금지.
6. 새 history 모델마다 thermal quadrature 수렴 재검사. 기존 5% gate는 각 모델의 수치 수렴 판단. 경계 민감도는 같은 Q·동일 quadrature의 paired 변화와 독립 seed spread 별도 표기. 두 값의 단순 차이만으로 엄밀한 5% 오차 상계라 부르지 않음. 수치 오차가 effect 판별 여유보다 작지 않으면 `INCONCLUSIVE`.
7. 반증 조건: 수렴한 두 history 가정에서 하나 이상의 declared total-stream metric 변화 >5%이며 수치 불확실성으로 설명되지 않음. 이때 “즉시 reset의 영향 ≤5%” 거부. 여기 5%는 **탐색적 경계 민감도 판정선**; 전체 물리 모델의 오차 허용치·상계 아님. 모두 ≤5%여도 시험한 R·모델·RF 범위에서만 안정성 관측. 실제 wall model의 정확도 증명 아님.

허용 결론: 고정 collection 아래 경계 상태 가정 민감도, 특정 조건의 source-history bookkeeping 검산, 수렴한 conditional atomic stream 변화. **허용 안 됨:** 전체 thermal physics 인증, Maxwell M(Ω)·광장 CP 인증, 실제 gain·squeezing dB 예측, 실제 셀의 boundary-history 오차 상계. Nonlocal Maxwell·finite seed·depletion·full Zeeman·실험 holdout은 별도 미완료.

코드 근거: [Rb path](../../../gabes/fwm_quantum/transport.py), [smooth adapter](../../../gabes/fwm_quantum/smooth_transport.py), [moment solver](../../../gabes/quantum/smooth_transport.py), [thermal model](../../../analysis/grand_challenge/rb_thermal_ensemble.py), [inflow](../../../gabes/quantum/inflow.py), [independent adjoint](../../../analysis/grand_challenge/reference/adjoint_transport.py).
