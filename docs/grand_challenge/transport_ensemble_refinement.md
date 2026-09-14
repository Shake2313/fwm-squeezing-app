# Constant-atom ensemble refinement — 같은 적분, 강화한 검사

2026-09-14. [계산](../../analysis/grand_challenge/transport_ensemble_refinement.py), [검사](../../tests/quantum/test_transport_ensemble_refinement.py), [새 보고서](transport_ensemble_refinement_report_v1.json), [plot](transport_ensemble_refinement_v1.png).

**고정 fixture 수렴 통과.** p16에서 최근 두 refinement 최대 2.312%, 마지막 두 grid의 directed seed 비교 최대 0.884%. 기존 5% 이내. Production adapter·physical Rb 인증은 false 유지.

이전 [adapter 유도](transport_ensemble_derivation.md)·[v1 보고서](transport_ensemble_report_v1.json)의 **같은 constant-atom fixture** 재계산. Boundary·integrand·path budget \(10^{-9}\)·stream budget 5% 유지. 경로 수와 독립 seed 수 증가. 정확한 factorization으로 계산량 절감. 실패 path 삭제·density 재조정·허용 오차 완화 없음.

Historical v1: 마지막 refinement 6.5152% > 5%, `certified=false`. 파일·실패 상태 유지. 새 계산의 수렴 여부는 별도 기록. Production callback 인증과 physical Rb 인증도 별도.

## 고정한 모델과 입력

| 항목 | 유지한 값·조건 |
| --- | --- |
| Atom | 2-level, \(H=0\), jump 없음 |
| Boundary state | 기존 coherent mixed density matrix |
| Readout·reciprocal drive | 기존 네 행렬, 모든 path에서 동일 |
| Phase | \(s=(1,-1,-1,1)\)의 uniform common mark |
| Carrier offsets·relative wavevectors | 모두 0 |
| Box | \(0.4\times0.2\times0.2\,\mathrm{mm}^3\) |
| Gas | 373 K, \(n=3\times10^{16}\,\mathrm m^{-3}\), 기존 \(^{85}\mathrm{Rb}\) mass |
| Lab RF | DC, 10, 30 kHz |
| Boundary quadrature | 기존 여섯 면 Maxwell flux, 실제 first-exit chord |

모든 물리·기하·상태 값은 가정한 fixture 입력. Rb mass 사용만으로 Rb optical model이 되지 않음. Common phase도 부여한 mark. 실제 nonzero lab beat dynamics 계산 아님. Source는 `atomic_inflow` 하나. Nonzero microscopic jump-source ensemble은 미검증.

## F·T factorization 유도

고정 density matrix \(\rho_0\), readout \(O_j\), drive \(V_k\). \(H=0\), jump 없음 → unperturbed state와 operators 일정. 모든 port에서 같은 \(\Omega\). 경로별 변수는 residence \(\tau\).

\[
\mu_j=\operatorname{Tr}(O_j\rho_0),\quad
K_{jk}=\operatorname{Tr}\{O_j[-i(V_k\rho_0-\rho_0V_k)]\}.
\]

Connected ordered matrices \(C_0^{>,<}\)도 일정. 유한 pulse와 causal triangle:

\[
\begin{aligned}
F_\tau(\Omega)
&=\int_0^\tau e^{i\Omega a}da
=\tau e^{i\Omega\tau/2}
  \operatorname{sinc}\left(\frac{\Omega\tau}{2\pi}\right),\\
T_\tau(\Omega)
&=\int_0^\tau da\int_0^a db\,e^{i\Omega(a-b)}
=\int_0^\tau(\tau-u)e^{i\Omega u}du
=\tau^2\frac{e^z-1-z}{z^2},\quad z=i\Omega\tau.
\end{aligned}
\]

\(\operatorname{sinc}(x)=\sin(\pi x)/(\pi x)\). DC limits: \(F_\tau(0)=\tau\), \(T_\tau(0)=\tau^2/2\). Small \(z\)에서 analytic Taylor series 사용. 취소 오차 방지용 동일 함수 표현.

\[
\mathbf m_\tau=F_\tau\boldsymbol\mu,\qquad
C_\tau^{>,<}=|F_\tau|^2C_0^{>,<},\qquad
\chi_\tau=T_\tau K.
\]

Common-phase mask \(M_{jk}=\delta_{s_j,s_k}\). 실제 path rate \(J_\ell\)로

\[
A(\Omega)=\sum_\ell J_\ell|F_{\tau_\ell}(\Omega)|^2,\qquad
B(\Omega)=\sum_\ell J_\ell T_{\tau_\ell}(\Omega).
\]

따라서 전체 행렬 복원:

\[
\begin{aligned}
S_{\rm internal}^{>,<}
&=A\,M\odot C_0^{>,<},\\
S_{\rm number}
&=A\,M\odot(\boldsymbol\mu\boldsymbol\mu^\dagger),\\
S^{>,<}
&=S_{\rm internal}^{>,<}+S_{\rm number},\\
\chi_{\rm stream}
&=B\,M\odot K.
\end{aligned}
\]

Source별 배열도 동일 복원. 이 fixture에서는 boundary source 하나가 internal 전체. Mean을 먼저 평균하지 않음. \(\langle\mathbf m\rangle_\phi=0\)이어도 \(M\odot\mathbf m\mathbf m^\dagger\) 유지.

\(B\)는 **complex** 가중합. 실수부만 남기면 response phase 손실. DC에서는 \(B=A/2\). 일반 RF에서는 그 관계 사용 불가.

\(J_\ell\)에 density·면적·Maxwell flux 이미 포함. 추가 \(n\), \(\tau\), velocity weight 곱하지 않음. \(A,B\) 단위 s. \(K\)에 reciprocal coupling의 \(s^{-1/2}\) 포함 → stream covariance s, response \(\sqrt{s}\). Response는 전체 age double integral. 공간 Maxwell kernel \(\mathcal R(\mathbf r,\mathbf r';\Omega)\), optical propagation \(M(\Omega)\), SQL spectrum 아님.

**정확 분리 조건:** 기존 고정 atom·state·readout·drive, \(H=0\), jump 없음, 모든 port offset=0, spatial wavevector=0, 동일 common phase rule. 위치·속도별 Hamiltonian, Gaussian pump, unequal port frequency, finite seed가 생기면 현재 분리 재사용 불가. `factorized_grid`는 공급 atom의 digest가 기존 fixture와 다르면 거절.

## Candidate는 Gauss 48, analytic 식은 독립 비교

모든 실제 path에서 세 계산 수행:

1. Gauss–Legendre 24-point: coarse \(F,T\).
2. Gauss–Legendre 48-point: candidate \(F,T\).
3. Sinc/exponential·DC Taylor 식: analytic reference.

최종 \(A,B\)는 **Gauss 48 값**으로 누적. Reference 값을 candidate에 대입하여 오차를 0으로 만들지 않음. 모든 boundary point·velocity·rate·face·residence 유지. Chunk는 작업 배열의 묶음 크기일 뿐. 대표 path 치환·표본 삭제 없음.

이전 일곱 path metrics 유지: `greater`, `lesser`, 두 `*_by_source`, `mean_pulse`, `mean_outer`, `retarded_response`. 행렬이 scalar×고정행렬이므로

\[
\|(c-r)X\|_F=|c-r|\,\|X\|_F
\]

사용. 기존 RF별 Frobenius error와 동일한 수치 기준. Vector mean에도 같은 norm identity 적용. 분모의 기존 단위 일치 zero-signal floor 유지. 24→48·48→analytic 두 비교를 **각 path의 모든 metrics**에서 수행. 예산 \(10^{-9}\).

추가로 모든 grid의 reconstructed stream 여섯 metrics를 analytic 가중합과 비교. 기존 `PHASE_BUDGET=2e-12` 유지. Path 검사만으로 누적 구현의 오차를 숨기지 않음.

완전 boundary arrays와 세 종류 \(F,T\) arrays의 SHA-256 chain 기록. Grid의 spectra·counts·오차·출처도 content digest로 결합. 실패 path 수가 0이 아니거나 counts·digest·metrics가 불완전하면 수렴 gate 거절. 실패 path를 빼고 rate를 재정규화하지 않음.

## 강화한 ensemble gate

사전 grid: 면당 \(2^{10},2^{12},2^{14},2^{16}\). 독립 seeds 11, 211, 811. 적어도 세 연속 grid 필요. 처음 통과 가능 지점 p14. 통과 못하면 p16까지. Path evidence 실패 시 인증 불가.

각 grid에서 모든 seed의 actual covariance·source·number·complex response를 다시 계산. 비교 metrics 여섯:
`greater`, `lesser`, 두 `*_by_source`, `poisson_number`, `retarded_response`.

최종 gate 구성:

| 비교 | 요구량 |
| --- | ---: |
| 최근 두 adjacent refinement × 세 seeds | 6 |
| 마지막 두 grid × 모든 directed seed pairs | 12 |
| 총 actual-output comparisons | 18 |

각 비교에서 여섯 metrics 모두 기존 RF별 최대오차 ≤5% 필요. Seed pair는 양방향 계산. Reference norm이 분모라 방향을 바꾸면 상대오차도 달라짐. 한 방향 통과만으로 다른 방향 통과를 추정하지 않음.

Grid·seed 누락, 중복, 미선언 조합 거절. 모든 grid의 density는 finite·positive이며 서로 정확히 같아야 함. Path·closed-form error도 finite·nonnegative여야 함. 최적 seed 선택 없음. 최근 한 edge나 마지막 grid만 우연히 통과한 경우도 불충분. Gate는 저장된 실제 행렬에서 error 재계산. \(nV\)·Maxwell moments는 진단값이며 인증 기준 대신 사용하지 않음.

이 강화 gate도 관측된 수치 수렴 검사. 연속 적분 전체의 엄밀한 오차 상한·통계적 confidence interval·실험 불확도 아님.

## 기존 adapter와의 연결, 인증 경계

낮은 두 grid에서 실제 production callback 결과와 factorized 결과 비교. p0/seed11: 6 paths. p2/seed211: 24 paths. Production mask와 별개인 네 literal phase 직접 합도 대조. 기존 \(2\times10^{-12}\) 예산 유지.

별도 실제 constant-atom ODE의 일곱 packet metrics 검사 유지. Historical v1의 p10/seed11 최종 spectra도 새 factorization으로 재현. Parent 보고서 SHA-256 고정:

~~~text
b2571c4777456b0a53ce528bfefecbc99c9dc013f475e79719f29948b9d9141c
~~~

고해상도 계산의 `packet_callback_count=0`. 실제 path 적분은 전부 수행하지만 production `SuppliedPathPacket`·per-path packet-digest chain은 만들지 않음. 따라서 상태 구분:

| 상태 | 의미 |
| --- | --- |
| `constant_atom_ensemble_converged` | 이 고정 fixture의 factorized actual-integrand 수렴 |
| `production_adapter_certified=false` | 고해상도 production callback 인증 근거 미구성 |
| `physical_Rb_ensemble_certified=false` | 실제 Rb thermal ensemble 검증 없음 |

같은-grid 일치·정확 분리 유도·전경로 evidence가 fixture 계산의 근거. 이를 production certification으로 소급 치환하지 않음. Historical v1의 미인증 candidate도 그대로. 실제 Rb에는 독립 path/source/response 검증과 실제 Rb integrand의 ensemble 수렴 필요.

## 실행 결과와 재현

최종 [refinement 보고서](transport_ensemble_refinement_report_v1.json): `implementation_controls_passed=true`, `expected_controls_passed=true`. `constant_atom_ensemble_converged=true`. `production_adapter_certified=false`, `physical_Rb_ensemble_certified=false`.

p10·p12에서는 세 연속 grid 조건 미충족. p14에서 p10→p12/seed211 변화 6.7092%로 실패. 선언한 다음 grid p16까지 계산. 최종 gate는 p12→p14→p16의 두 edge와 p14·p16의 모든 directed seed pairs. 18개 비교 모두 통과.

| 최종 비교 | 비교 수 | 최대 matrix error |
| --- | ---: | ---: |
| p12→p14→p16, 세 seeds | 6 | 2.311645% |
| 그중 p14→p16 | 3 | 0.385750% |
| p14·p16의 모든 directed seed pairs | 12 | 0.883543% |

두 번째 행은 첫 행의 부분집합. 추가 gate로 중복 계산하지 않음. Stream budget 5% 그대로. 세 seed 중 유리한 결과만 선택하지 않음.

Grid별 path 평가수 합 1,566,720. 최종 seed당 393,216 paths. 모두 Gauss24·Gauss48·analytic 적분 수행. 고해상도 packet callback 0회. 최종 세 seeds의 최대 점유수 상대오차 약 0.01879%; nominal \(N/(nV)=1.000187815015\). 점유수는 진단값. 위 actual-matrix gate가 수렴 판단 근거.

| 독립·구현 검사 | 최대 상대오차 | 기존 예산 |
| --- | ---: | ---: |
| 전경로 일곱 metrics, 두 비교 종류 | \(4.529\times10^{-14}\) 미만 | \(10^{-9}\) |
| 모든 grid의 Gauss48 stream versus closed form | \(3.296\times10^{-16}\) 미만 | \(2\times10^{-12}\) |
| 같은 낮은 grid의 production callback parity | \(3.057\times10^{-16}\) 미만 | \(2\times10^{-12}\) |
| 독립 네 literal phase 직접 합 | \(4.704\times10^{-16}\) 미만 | \(2\times10^{-12}\) |
| Historical v1 p10/seed11 재현 | \(2.943\times10^{-15}\) 미만 | \(2\times10^{-12}\) |
| 실제 constant-atom ODE의 일곱 metrics | \(9.754\times10^{-16}\) 미만 | \(2\times10^{-8}\) |

Historical v1의 실패 flag·bytes·SHA-256 유지. 동일 p10 결과 재현 확인 후 더 높은 grid 수렴만 새로 기록. 과거 production candidate를 소급 인증하지 않음.

관련 tests **160 passed (3.35 s)**: 신규 refinement 51 + adapter 66 + inflow 43. 실제 audit **9.564 s**, BLAS 한 thread. 수치 비교·집계·artifact 실행 기록이며 일반 성능 보장 아님.

주작업의 최종 `python -m pytest -q`: **1264 passed, 1 failed (232.46 s)**. 실패는 기존 `docs/FWM physics and analytic reconstruction/FWM_physics.tex` 누락. 이번 병렬 변경의 신규 실패 없음. 메인 adjoint 28 + refinement 51 = 신규 79 tests 포함.

최종 독립 artifact 검사: 12개 grid의 실제 complex matrices에서 18개 비교·108개 metric errors 재계산. 저장값과 차이 최대 \(1.74\times10^{-18}\). Row content digests·source 합·최종 spectra·17개 source hashes·historical parent hash 일치.

실제로 import한 project modules·실행 driver·소비한 tests의 source hash 17개: 실행 전후 동일, 최종 파일과 일치. 무관한 병렬 reference 파일 전체를 glob하여 포함하지 않음. 새 보고서 SHA-256:

~~~text
82bf1f805257145e06091a51babe62ba4fd865eda61aef3b4fa47aca8f89c9bf
~~~

PNG SHA-256:

~~~text
4d23d9aa2323f7a296439df4556e32e33f38eac8032b9d658a004f2e1147410c
~~~

~~~powershell
python -m analysis.grand_challenge.transport_ensemble_refinement --output docs/grand_challenge/transport_ensemble_refinement_report_v1.json --plot docs/grand_challenge/transport_ensemble_refinement_v1.png
python -m pytest -q tests/quantum/test_transport_ensemble_refinement.py
python -m pytest -q
~~~

Report/plot 기존 경로에 덮어쓰기 거절. 재실행은 두 출력 모두 새 version 이름 사용. Source hash 실행 전후 불일치 시 artifact 저장 중단.
