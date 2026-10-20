# Dressed 흡수·방출 분리 — 2026-10-16

## 사용자 요구

순 대각 이득을 투과율로 표시하거나 계산 불가로 두는 대신, Jasperse의
conjugate 손실 계수에 대응하는 물리량을 현 GABES에서 추출할 것.

## ⟦AI⟧ 결과

가능. 현 4준위 Lindblad 모델의 **양의 흡수 스펙트럼**과 **양의 유도
방출 스펙트럼**을 기존 quantum-regression 엔진으로 분리. 순 감쇠만
보면 잃는 두 양의 항을 ordered noise에서 각각 계산. 임의 계수·피팅 없음.

물리적 annihilation port의 local ordered field noise를 `D⁾`, `D⁼`라 할 때

\[
\alpha_{\rm abs,j}=D^{>}_{jj},\qquad
\alpha_{\rm em,j}=D^{<}_{jj},\qquad
\alpha_{\rm net,j}=\alpha_{\rm abs,j}-\alpha_{\rm em,j}
                  =-2\Re A_{jj}.
\]

단위 m⁻¹. Creator Nambu 좌표에서는 greater/lesser 순서 반대.
모든 속도 클래스에서 각각 평균한 뒤

\[
\eta_{{\rm in\ cell},j}^{\rm abs}
=\exp[-L\langle\alpha_{\rm abs,j}\rangle_v]
\]

로 표시. `exp[−L〈α_net〉]`와 구분. Jasperse의 순수 손실 `T_b`와
같은 흡수 항이며, 방출·위상 민감 잡음이 함께 존재하면 `T_b` 하나만으로
전체 GABES quantum channel을 재현하지는 못함.

### 현 dissipator 그대로 재구성

현재 ground collision `γ_g`, optical collision `γ_o`에 대해

\[
J_g=\sqrt{\gamma_g/2}\,\mathrm{diag}(1,-1,0,0),\qquad
J_o=\sqrt{2(\gamma_o-\gamma_g/4)}\,\mathrm{diag}(0,0,1,1).
\]

`γ_o≥γ_g/4`이면 완전 양의 두 jump가 기존 coherence 감쇠를 정확히
재현. 기존 radiative branching과 thermal-reset jump 추가.
예시 운전점의 전체 dissipator 상대 차이 **1.44×10⁻¹⁶**.
추가 damping·온도·noise scale 도입 없음. 부등식 불만족 시 계산 거부.

### 계산 경로·독립 검증

`local_dressed_absorption(method="ordered-diffusion")`은 기존
`periodic_atomic_noise` → `eliminate_periodic_atom`을 그대로 연결.
ordered PSD·CCR·동적 Einstein 검사를 통과한 `R D R†` 기준 경로.

빠른 경로는 같은 주기적 상태 `ρ_n`에서 connected covariance
`C_n=〈F_i F_j〉_n−(〈F_i〉〈F_j〉)_n` 구성. 평균의 convolution을
`|n|≤2N_mean`까지 보존. `C_n`, `C_nᵀ` 두 초기원으로 직접 QRT 계산.
별도로 `−i[D†,ρ_n]/2` commutator 원을 세 번째 RHS로 조립하여
`α_abs−α_em`과 retarded 응답의 일치 검수. 15차원 block-Thomas 풀이.
평균장 상태가 양·trace 1·Hermitian인지 검사. 양의 thermal reset이
trace-zero sector 수축을 보장하는 조건에서만 빠른 경로 사용.

| 독립 대조 | 결과 |
|---|---:|
| pump-only χ → Maxwell와 ordered-noise drift | 상대 차이 3.57×10⁻¹⁴ |
| zero-seed periodic → stationary | 상대 차이 2.69×10⁻¹⁵ |
| 9개 `(δ,v)` direct QRT → `R D R†` | 최대 상대 차이 6.95×10⁻¹² |
| `δ=0,23.7,52.269382 MHz`, `v=−200,0,+200 m/s` | 위 9점 범위 |
| actual seed에서 harmonic 3 → 4, v=0·resonance root | 최대 상대 차이 3.51×10⁻¹⁴ |

Response harmonic 1은 CCR 검사를 통과하지 못해 채택하지 않음.
독립 검토 경로: [직접 QRT 유도](qrt_absorption_derivation.py),
[9점 대조](qrt_connected_controls.json),
[기존 engine 기준](qrt_efficiency_probe.json).

### 조건부 실제 설정 예시

`T=386.15 K`, pump `294.8 mW`, `w_p=530 µm`, `w_s=330 µm`,
OPD `+1 GHz` blue, seed 아래쪽, internal `δ=+23.7 MHz`,
`L=12.5 mm`, 가정한 vapor-entry seed `3.917987 µW`, 기존 weak-field
normalization `0.74/12`. 동일한 **probe-only 유한 seed 상태**에서 계산.

| 위치·항 | Probe | Conjugate |
|---|---:|---:|
| v=0 `α_abs` [m⁻¹] | 53.96476 | 5.025669 |
| v=0 `α_em` [m⁻¹] | 38.56414 | 7.309736 |
| v=0 흡수만의 투과율 | 0.509381 | 0.939112 |
| Maxwell 평균 `α_abs` [m⁻¹] | 126.97683 | 15.375923 |
| Maxwell 평균 `α_em` [m⁻¹] | 102.10175 | 19.193586 |
| Maxwell 평균 흡수만의 투과율 | **0.20449529** | **0.82514261** |
| Maxwell 평균 순 대각 강도비 | 0.73275893 | 1.04887775 |

Conjugate가 순 이득이더라도 흡수 계수는 양수이며 투과율 ≤1.
흡수만의 효율과 유도 방출을 함께 보존하므로 `max(α_net,0)`로
만든 값과 다름.

속도는 기존 production과 같은 collinear `Δ_eff=Δ−K_VEC v`와
±4σ normalized Maxwell. Dressed 에너지 차로 찾은 Raman root
`v=−635.907803 m/s`를 적분 breakpoint에 포함. Adaptive Gauss–Kronrod
210회 평가, NF3/4 control 포함 214개 local 계산. Numeric 약 **0.77 s**;
cache hit 약 **10 ms**. Full ordered engine 기준 동일 평균 약 6.65 s.

0.5%→0.1% 적분 설정 변경의 값 차이 <10⁻¹⁰ 상대.
빠른 wrapper가 보고한 흡수율의 **추정** 상대 적분 오차는
probe 2.70×10⁻⁵, conjugate 2.08×10⁻⁵. 엄밀한 오차 상계나
Maxwell tail 수렴 증명은 아님. 후자는 4σ 규약으로 한정.
원자료: [Doppler 기준](qrt_doppler_probe.json),
[현 wrapper](qrt_selected_wrapper.json).

### 재사용 API·해석 범위

[`selected_dressed_cell_absorption`](../../../gabes/dressed_cell_absorption.py)
입력: `D_GHz,T,P_pump,P_probe,w_pump,w_probe,delta_mhz,L,
weak_line_strength,transit_rate`. `weak_line_strength`에는 현행
`effective_line_strength` 입력. 입력 seed 기준면은 atomic entrance.
반환: `probe/conjugate`별 `alpha_abs_per_m`, `alpha_em_per_m`,
`alpha_net_per_m`, `in_cell_absorption_transmission`, 별도 net transmission,
적분·harmonic 검수 정보. 마지막 32개 운전점 cache.

현 유한 seed χ는 probe-only·conjugate-only **서로 다른 상태**의
coherent quotient를 조합. 새 QRT는 **하나의 probe-only 상태**에서
두 팔 fluctuation을 계산하므로 동일한 유한 seed differential map으로
간주할 수 없음. 균일 입력 pump·seed, coherent conjugate 0을 셀 전체에
가정. 생성광과 depletion을 따라 적분한 평균장 경로·Zeeman reduction·
각도 Doppler·실험 투과율·physical squeezing 검증으로 승격하지 않음.

`η_after cell η_QE`와 곱한 total efficiency는 요청한 loss bookkeeping.
원래 atomic transfer에 이 셀 손실을 또 적용하지 않음. 방출 및 두 팔
상관 잡음은 별도 채널로 남으며, 본 변경이 GC의 전체 quantum propagation
문제를 닫지는 않음. 기존 microscopic GC 함수·봉인 자료 수정 없음.

### 검사

`python -m pytest -q tests/test_dressed_cell_absorption.py` → **16 passed**.
정확한 dissipator, independent `R D R†`, zero-seed susceptibility,
활성 conjugate의 양의 흡수, Raman pole harmonic 수렴, Maxwell 적분,
cache 격리, 순수 손실 한계, 잘못된 입력 검수. 전체 회귀는 메인 담당.
