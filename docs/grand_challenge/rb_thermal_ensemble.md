# 열적 Rb stream: 경로 재사용·수렴 증거

최종 v3: **p0/seed 11의 고정된 6경로 모두 검증 통과**, 진단용 atomic spectrum 생성.
독립 reference 12개 완료. V1·v2 실패 산출물도 불변 보존.
**전체 thermal ensemble·광학 예측은 미인증.** 다른 11개 grid/seed 미평가.

기반: [연속 Gaussian 경로](smooth_transport_derivation.md),
[독립 adjoint 검증](adjoint_transport_derivation.md),
[pump-only ensemble 계약](transport_ensemble_derivation.md).
실행기: `analysis/grand_challenge/rb_thermal_ensemble.py`.
Pilot audit: `analysis/grand_challenge/rb_thermal_ensemble_audit.py`.

**물리 범위.** Undepleted pump-only reduced ⁸⁵Rb D1.
실제 signed wavevector·경로별 Doppler·공간 입사 위상 유지.
Gaussian pump는 경로 전체에서 연속 평가. 궤적·체류시간 자르기, 속도 평균,
pump 평균, phenomenological transit reset 없음.
Finite seed saturation·full Zeeman·nonlocal Maxwell·검출기 SQL·실험 squeezing 인증은 범위 밖.

**새 경계 가정.** 기존 단일 prescribed path에서 열적 경계로 확장.
기본 영역은 측면 길이 $\sqrt{A_z}$, 길이 $L$인 정사각 열린 기둥.
$A_z=1.2\times10^{-7}\,\mathrm{m^2}$, $L=12.5\,\mathrm{mm}$.
기존 optical normalization 면적을 물질 유입 경계 면적으로 사용하는 **새 가정**.
측정한 실제 cell 벽·beam aperture로 해석하지 않음.

여섯 면에서 Maxwell incoming flux 유입. 매 입사 상태
$\rho_{\rm in}=\operatorname{diag}(5/12,7/12,0,0)$.
외부 reservoir에서 새 원자 공급, 첫 출구에서 경로 종료.
경계 밖 pumping 이력·재진입 상관·벽 충돌·코팅 memory 없음.
기본 $T=373\,\mathrm K$, $n=10^{18}\,\mathrm{m^{-3}}$는 별도 지정값.
온도로 밀도를 추정·맞춤하지 않음.

Pump $P=0.6\,\mathrm W$, $w=530\,\mathrm{\mu m}$;
진폭 envelope $f=\exp[-r_\perp^2/w^2]$.
기본 **측면** 경계에서 $f\simeq0.808$–$0.899$,
광강도 $I/I_0=f^2\simeq0.652$–$0.808$.
$z$ 입·출구 면 중심에서는 $f=1$.
따라서 이 경계는 pump가 사라진 먼 영역이 아님.
밖에서도 pumping된 원자가 들어오는 실제 cell에는 현재 입사 상태 가정 추가 검증 필요.

**물리 경로 캐시.** 저장 대상: boundary rate 곱하기 전 원자 packet.
키에는 정확한 입사 위치·속도·체류시간·공통 입사 위상, 전체 SI 입력,
carrier geometry·RF axis·입사 밀도행렬·pump 중심·규약 포함.
Solver family·해상도·세부 parameter, 의존 소스 SHA-256,
Python·NumPy·SciPy 버전도 포함. 좌표 반올림으로 다른 경로 합치지 않음.

같은 seed에서 Sobol power 증가 → 각 face의 기존 경로 재등장.
그러나 face별 블록 길이가 바뀌므로 전역 node 번호와 `path.source`는 달라질 수 있음.
캐시 키에서 grid power·전역 index·`path.source`·node rate 제외.
원본 source와 cache payload hash는 이력에 보존.

실제 재사용 검사: p0/seed 11의 6경로 ↔ p1/seed 11의 indices 0·2·4·6·8·10.
입사 위치·속도·$\tau$ 완전 동일, rate 정확히 0.5배.
V1 finest spec으로 두 node 집합 조회: **12/12 cache hit**.
각 대응 pair의 key·payload·8개 raw metric 완전 동일.
Faces 1…5는 새 source에 따라 packet digest 갱신;
face 0은 source도 같으므로 같은 digest 유지가 정상.
소스 62개·cache 파일 6개 불변 확인, provider·ODE 호출 0회.
이 검사는 기존 6경로 재사용 확인. P1의 나머지 경로·전체 stream spectrum은 미평가.

현재 경로에서 packet shape·RF·Doppler·공간 위상·source 합·Hermiticity·PSD·출구 상태 재검증.
Primary 각 해의 `packet_digest`를 현재 `path.source`로 다시 계산하고,
저장 배열에서 비교 오차 재산출. 현재 candidate/reference 연결로 evidence 재구성.
이전 node의 인증 객체를 그대로 옮기지 않음.

입력 배열은 불변 복사, solver의 중첩 parameter도 동결.
Source manifest·계산 identity·payload hash로 저장값 확인.
Cache hit/miss 모두 소스 해시 전후 검사.
기존 JSON은 배타적 생성으로 덮어쓰기 금지; 부분·손상 record는 거부.
의존 소스 manifest는 실제 provider와 그 의존성을 포함해야 함.
해시 일치는 계산 출처·내용 확인이며 물리 모델의 타당성 증명은 아님.

독립 reference는 기존 35개 소스 manifest·불변 job identity 검증 후 연결.
같은 입력·경로·$H_0,H_1,\rho_{\rm in},O,V$·RF·물리 metadata 대조.
Cache hit에서도 반환된 8개 metric·metadata·numerics를 검증된 원본 job과 다시 대조.
원본 refinement와 경로 evidence가 서로 다른 reference를 쓰는 경우 거부.

**경로 하나의 통과 조건.** Primary 3개 해상도 + 수학적으로 독립된 reference 2개 해상도.
계산값 총 5개 필요. 기존 cache에서 읽어도 아래 오차를 실제 배열로 다시 계산.

| 비교 | 최대 상대 오차 |
| --- | ---: |
| Primary의 연속 두 refinement 각각 | $10^{-3}$ |
| Finest primary vs finest independent reference | $5\times10^{-6}$ |
| Independent reference 두 해상도 | $2\times10^{-6}$ |

검사량: greater·lesser, 각각의 모든 named source,
mean pulse·mean outer·복소 retarded response·exit state.
Source별·RF별 Frobenius norm; mean은 RF별 Euclidean norm.
다른 source나 강한 RF row에 약한 항의 오차 숨기지 않음.

분모는 reference norm과 $128\epsilon_{\rm mach}$ 곱하기 선언된 SI scale 중 큰 값.
경로 covariance scale은 $\tau^2\|O\|^2$, mean은 $\tau\|O\|$,
response는 $\tau^2\|O\|\|V\|$, exit state는 1.
Floor는 비교 분모이며 스펙트럼에 더하는 잡음이 아님.
NaN·음수 오차·실패 audit·누락 source·누락 해상도는 통과 불가.

**수치 방법·검증.** Primary는 CF4 + fixed-generator spectral source 적분.
시간 의존 Hamiltonian의 근사 오차는 경로별 refinement와 독립 reference로 검사.
고정 generator의 spectral 적분 항등식만으로 전체 경로 정확도를 인증하지 않음.
Selected 2 µs 경로·fixed-map 공식 검증:
[exponential report v1](exponential_transport_report_v1.json) PASS.

| Selected-path 비교 | 최대 상대 오차 | 기준 |
| --- | ---: | ---: |
| CF4 4096 → 8192 | $9.4185276\times10^{-4}$ | $10^{-3}$ |
| CF4 8192 → 16384 | $1.6138011\times10^{-4}$ | $10^{-3}$ |
| Finest CF4 vs independent reference | $8.6908012\times10^{-7}$ | $5\times10^{-6}$ |

이 결과로 thermal primary 실행 시작. 새 열적 경로의 5개 해 비교는 각각 별도 수행.

열적 pilot **v1**의 세 maximum step:

\[
\Delta t_{\max}=\left(
 \frac{2\,\mathrm{\mu s}}{4096},
 \frac{2\,\mathrm{\mu s}}{8192},
 \frac{2\,\mathrm{\mu s}}{16384}\right),\qquad
N_j=\left\lceil\frac{\tau_j}{\Delta t_{\max}}\right\rceil,
\qquad \Delta t_j=\frac{\tau_j}{N_j}\leq\Delta t_{\max}.
\]

`step_plan`으로 선언. 긴 경로는 필요한 segment 수 증가.
2 µs는 step을 정하는 수치 기준. 각 경로의 실제 $\tau_j$·밀도·영역은 그대로.
Reference 두 tolerance는 `(rtol, atol)=(1e-9, 1e-12)`, `(3e-10, 3e-13)`.
각 reference의 `max_step_s=τ/64`도 계획과 원본 numerics에 명시.

V1 전체 결과: paths **1·2·4 FAIL**, paths 0·3·5 PASS.
실패 경로의 두 edge 최대 오차는 모두 mean outer. 기준은 각각 $10^{-3}$.

| Path | 첫 primary edge | 다음 primary edge | 경로 gate |
| --- | ---: | ---: | --- |
| 1 | $2.2574663\times10^{-3}$ | $3.9706150\times10^{-4}$ | FAIL |
| 2 | $1.6099793\times10^{-3}$ | $3.9118190\times10^{-4}$ | FAIL |
| 4 | $1.0134718\times10^{-2}$ | $2.5510686\times10^{-3}$ | FAIL |

Finest vs independent reference는 모든 경로 통과.
최대 $2.4465310\times10^{-6}<5\times10^{-6}$, path 1 mean outer.
두 연속 edge가 모두 필요한 계약이므로 뒤 비교 통과로 첫 실패를 지우지 않음.
Path 0·1의 저장 배열·hash·원본 reference 연결은 별도 읽기 전용 검사에서도 확인.

V1 실패 report·plot 생성 완료, 불변 보존.
V2는 **모든 6경로**에서 v1 maximum step을 각각 절반으로 축소:
$(2\,\mathrm{\mu s}/8192,\ 2\,\mathrm{\mu s}/16384,\ 2\,\mathrm{\mu s}/32768)$.
V2 첫 edge는 v1 두 번째 edge 재사용.
V2 완료 결과: **path 4 첫 edge만 FAIL**, 오차 $2.5510686\times10^{-3}$.

**V3 최종 검증.** 모든 6경로에 같은 maximum-step 계획 적용:
$(1.220703125\times10^{-10},\ 6.103515625\times10^{-11},\ 3.0517578125\times10^{-11})\,\mathrm s$.
같은 밀도·경계·입사 상태·실제 체류시간·오차 기준 유지.
실제 macro-segment 수 $N_j=\lceil\tau_j/\Delta t_{\max}\rceil$:

| Path | 실제 $\tau$ (µs) | Coarse $N$ | Middle $N$ | Fine $N$ |
| --- | ---: | ---: | ---: | ---: |
| 0 | 1.064721504 | 8723 | 17445 | 34889 |
| 1 | 0.380078373 | 3114 | 6228 | 12455 |
| 2 | 0.996237328 | 8162 | 16323 | 32645 |
| 3 | 0.330242817 | 2706 | 5411 | 10822 |
| 4 | 1.073989402 | 8799 | 17597 | 35193 |
| 5 | 0.194352308 | 1593 | 3185 | 6369 |

표의 $\tau$만 표시용 반올림. 실제 계산·키는 원래 값 사용.
6경로·모든 source·RF·8개 metric 중 최악 오차:

| 비교 | 최대 상대 오차 | 기준 | 결과 |
| --- | ---: | ---: | --- |
| Primary coarse → middle | $2.3432344\times10^{-6}$ | $10^{-3}$ | PASS |
| Primary middle → fine | $9.4625991\times10^{-8}$ | $10^{-3}$ | PASS |
| Fine vs independent reference | $1.8223767\times10^{-8}$ | $5\times10^{-6}$ | PASS |
| Independent reference refinement | $3.4601050\times10^{-8}$ | $2\times10^{-6}$ | PASS |

Primary 두 최악값: path 1 mean outer. Independent 두 최악값: path 4 mean outer.
새 finest 6개 계산은 기존 worker로 병렬 실행 완료.
[실제 job manifest](rb_thermal_finest_jobs_v1.json): workers 4, 각 BLAS 1 thread.
2026-09-14 UTC **07:40:40.565933 → 07:51:08.012928**, elapsed **627.4470053 s**.
각 job의 실제 환경·시간·cache key 보존.
Manifest SHA-256: `d7090ee4e396d63d3d1c0aa2604c4ac06fdf874b3bfbd3f6ca761909713a449e`.

Thermal primary 신규 계산 총 **30개**: v1 18개 + v2 추가 6개 + v3 추가 6개.
독립 reference 12개 포함, 원자 해 cache 총 **42개**.
V3 audit의 primary 요청 18개는 모두 cache hit.
최종 grid callback 6회 × 경로당 5해 = **30 cache hit**, cache miss 0.
V3 report elapsed **18.3182300 s**는 검증·집계 시간.
42개 원자 해를 새로 계산한 시간이나 solver 가속 배수로 해석하지 않음.

**실행 환경·검증 기록.** 현재 `threadpoolctl 3.6.0` 설치.
새 reference worker는 실행 context 안에서 모든 실제 BLAS runtime의
`num_threads == 1` 확인하고 `runtime_blas_threads` 저장. 확인 실패 시 계산 거부.
Driver `rb_thermal_reference_jobs.py` SHA-256:
`8435761140df87cf4cd71771b871113600a0d007632707343abbe7a4a9afe480`.
이는 driver 한 파일의 hash. 35개 소스 manifest 전체 digest로 해석하지 않음.

과거 보고서의 `blas_single_thread` 호출·BLAS 1 label만으로 실제 thread 수 입증 불가.
당시 패키지 부재 시 helper가 no-op였을 수 있음.
과거 immutable 보고서는 보존. 새 실행과 같은 thread 조건이었다는 성능 비교 전제 없음.
이번 문서에서 과거 대비 가속 배수 주장 없음.

신규 thermal 계층 132개 + CF4 관련 91개, 대상 테스트 **223개 통과**.
합동 `python -m pytest -q`: **1487 passed, 1 failed, 338.88 s**.
기존 실패: `test_docs_consistency.py::test_repository_visibility_wording_is_consistently_public`.
기존 삭제 파일 `docs/FWM physics and analytic reconstruction/FWM_physics.tex`를
읽다가 `FileNotFoundError`. 전체 테스트가 전부 통과한 상태는 아님.
이 수치는 통합 작업의 실행 결과. Thermal ODE 완료 수·ensemble 수렴과 구분.
이후 계산 동안 source·test 동결 유지. V3 소스 62개 실행 전후 일치.

**Rate는 한 번.** Face별 $2^p$ node,
$\sigma=\sqrt{k_BT/m}$일 때 node rate

\[
\lambda_j=\frac{nA_{\rm face}\sigma}{\sqrt{2\pi}\,2^p}.
\]

Power가 1 증가하면 재사용 경로의 rate는 절반.
항상 현재 grid rate로 전 경로 합산. Density 추가 곱셈·occupancy 재정규화 없음.
계산된 $\sum_j\lambda_j\tau_j$와 해석적 $nV$의 차이도 남김.

하나의 공통 입사 beat phase를 평균.
$s=(1,-1,-1,1)$, $M_{ab}=\delta_{s_a,s_b}$이면

\[
S^{\gtrless}(\Omega)=\sum_j\lambda_j\,
 M\circ\left[\sum_r C^{\gtrless}_{j,r}(\Omega)
                  +\mu_j(\Omega)\mu_j^\dagger(\Omega)\right],
\qquad
R(\Omega)=\sum_j\lambda_j\,M\circ R_j(\Omega).
\]

$r$: `atomic_inflow`와 모든 실제 jump source.
Mean outer는 Poisson 원자 수 잡음으로 한 번 추가.
평균 위상이 0이어도 이 항은 제거하지 않음.
결과는 period-averaged zero-cyclic atomic stream.
Deterministic mean line·다른 cyclic spectrum 제외.
Retarded response의 복소부 유지; 공간 Maxwell 행렬로 해석하지 않음.

V3 저장 raw packet에서 위 rate·phase mask·source 합·mean outer를 직접 재합산.
보고서의 총 **8개 행렬 배열과 `np.array_equal` 완전 일치**.
위상 평균에서 제거되는 성분은 정확히 0.
복소 response의 $\|\operatorname{Im}R\|_F=0.006023184794457149$;
허수부 보존 확인. 보고서 수치 단위: covariance s, response $\sqrt{\mathrm s}$.

**Grid 진단과 ensemble 인증.** 한 grid의 모든 경로 통과 → 진단용 spectrum 생성.
이 단계의 `certified=False`. $nV$·속도 moment 일치만으로 인증하지 않음.
최종 v3도 이 상태: `path_evidence_passed=True`, `certified=False`.
실제 occupancy/$nV$ = **0.617675044895484**. 계산된 차이 그대로 유지.

Ensemble gate: 같은 모델·수치 계획에서 최소 3개 grid, seeds 11·211·811.
최근 두 refinement를 각 seed에서 검사: 6개 비교.
마지막 두 grid에서는 모든 directed seed pair 검사: 12개 비교.
총 18개 비교 모두 5% 이내 필요.
Greater·lesser·모든 source·Poisson number·복소 response 각각 RF별 검사.
같은 candidate digest를 서로 다른 seed의 독립 계산으로 인정하지 않음.
행렬·metadata digest 및 각 경로 ledger의 실제 8개 오차도 재검사.
관측된 수치 수렴 조건이며 통계적 신뢰구간은 아님.
V3 `ensemble_gate=False`: 필요한 추가 grid와 독립 seed가 없음.
Thermal ensemble·production adapter·physical optical prediction 인증 모두 `False`.

한 경로라도 실패 → `spectra=None`, 실패 사유·경로·현재 rate 기록.
부분 spectrum 숨김; 누락 경로 제외 후 재정규화 금지.
마지막 packet 검증·최종 합산 실패도 별도 진단에 남김.

**미평가 범위.** 기본 모델·물리 경로로 생성한 workload는 p0…p3, 세 seed.
전체 평가 시 grid 경로 항목 합계 270개.
중첩 경로의 계산 재사용 가능; 270은 고유 경로 수나 ODE 완료 수가 아님.

| Power | Face당 node | Grid당 경로 | 세 seed 경로 항목 |
| --- | ---: | ---: | ---: |
| 0 | 1 | 6 | 18 |
| 1 | 2 | 12 | 36 |
| 2 | 4 | 24 | 72 |
| 3 | 8 | 48 | 144 |
| 합계 | — | — | 270 |

첫 실제 pilot: p0, seed 11, 고정된 6개 경로.
각 경로의 두 adjoint reference, **12개 모두 완료·검증**.
[불변 reference index](rb_thermal_reference_jobs_v1/index.json) SHA-256:
`8921eb4ad2c2456018d539f861aa4e561d903efb012ec622d43125d9ec0b95f2`.
독립 재계산: 모든 경로·8개 metric refinement 통과.
최대 $3.4601050\times10^{-8}<2\times10^{-6}$;
path 4, mean outer, RF 4 MHz. 소스 35개·원본 identity·record hash도 대조 완료.
이 reference와 v3 primary 3개 해상도로 6경로 각각 통과.
완료 범위는 한 grid의 진단. Ensemble 인증 없음.
나머지 p0…p3 × 세 seed 조합은 `NOT_EVALUATED` 유지.

**불변 산출물.**

| 버전 | Report | 그림 | 경로 검증 |
| --- | --- | --- | --- |
| v1 | [JSON](rb_thermal_ensemble_report_v1.json) | [PNG](rb_thermal_ensemble_v1.png) | Paths 1·2·4 실패, `spectra=None` |
| v2 | [JSON](rb_thermal_ensemble_report_v2.json) | [PNG](rb_thermal_ensemble_v2.png) | Path 4 실패, `spectra=None` |
| v3 | [JSON](rb_thermal_ensemble_report_v3.json) | [PNG](rb_thermal_ensemble_v3.png) | 6경로 통과, 진단 spectrum |

최종 v3 report SHA-256:
`c5cfc2142020b38648ab8ccbbc7a68ad777690647057cdf33c6ee13848f2737f`.
V1·v2 실패 기록 보존. 같은 물리·기준에서 수치 refinement를 추가한 이력.
세 PNG는 동일: 공통 geometry·residence·reference refinement만 표시.
Primary refinement 판정은 각 JSON과 위 표에 있음.

**재현.** 저장소 root에서 실행. 기존 reference·cache는 원본 해시 검증 후 재사용.
Reference driver는 workers 2 고정; 별도 `--workers` 옵션 없음.

```powershell
python -m analysis.grand_challenge.rb_thermal_reference_jobs --output-dir docs/grand_challenge/rb_thermal_reference_jobs_v1
```

Audit은 새 output·plot 이름 필수. 아래는 같은 cache·최종 step·workers 4 사용.
현재 동결 소스와 완성 cache에서는 primary 18개 모두 hit.

```powershell
python -m analysis.grand_challenge.rb_thermal_ensemble_audit `
  --reference-dir docs/grand_challenge/rb_thermal_reference_jobs_v1 `
  --cache-dir docs/grand_challenge/rb_thermal_path_cache_v1 `
  --max-steps-s 1.220703125e-10 6.103515625e-11 3.0517578125e-11 `
  --compute-primary --workers 4 `
  --output docs/grand_challenge/rb_thermal_ensemble_report_v3_recheck.json `
  --plot docs/grand_challenge/rb_thermal_ensemble_v3_recheck.png
```

예시 산출물이 이미 있으면 새 이름 지정. 실행 성공도 한 grid 경로 검증 범위에 한정.
