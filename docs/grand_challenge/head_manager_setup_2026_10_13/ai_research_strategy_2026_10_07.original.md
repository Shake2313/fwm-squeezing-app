# Grand Challenge 개선을 위한 AI 연구 기법 검토

검토일: 2026-10-07. 과학 근거: 2026-09-28 봉인 자료. 기존 관리·코드 감사: 2026-10-02.

**권고: 다음 계산이 병목 원인을 구분하도록 설계.** 우선 기존 캐시로 분산 감소 후보를 검토하고, 현재 p3/seed811 pilot에서 원래 estimator를 평가. 경계 이력과 작은 nonlocal 광장 연결은 독립 소과제로 검증. 대규모 토큰·비공개 모델·추가 학습은 실행 전제로 삼지 않음.

이 문서는 외부 사례를 평가한 **개선 제안**. 새 원자 계산·수렴 인증·전략 채택 기록 아님. 현행 실행 권한과 카드는 같은 [blueprint](blueprint.md)에 있음. 공식 완료 **0/4**, no-fit 계약·5% 수렴 gate·봉인된 v2 실패 유지.

## 목적과 실제 진행 상황

목적: 독립 측정 입력으로 hot ⁸⁵Rb D1 double-Λ FWM의 절대 probe/conjugate gain과 전체 intensity-difference spectrum S₋(Ω)을 동시에 예측. 평균장·응답·원자 잡음·공간 전파·검출이 같은 모델과 정규화를 사용해야 함. Parameter-free는 입력이 없다는 뜻이 아니라 목표 gain/squeezing에 맞춘 자유 계수가 없다는 뜻. [공식 계약](../checklist.json), [총괄 설계](blueprint.md).

| 판단 대상 | 현재 증거와 제한 |
|---|---|
| GC-1 reduced theorem/model | 진행 중. 기초·단일 원자·조건부 field/readout·선택 경로 등 여섯 scoped WP 종결. 실제 thermal/nonlocal 연결 미종결 |
| GC-2 독립 입력 절대 예측 | 대기. GC-1과 독립 입력 ledger 필요 |
| GC-3 새 조건 검증 | 외부 raw gain/RF/SQL/calibration과 untouched 조건 필요 |
| GC-4 full atom·physics ablation | 검증된 reduced 경로와 잔차 필요 |
| Frozen v2 | 96고유 경로·480계산 기록·4/9격자 감사. 경로·합산 감사 통과, 앙상블 수렴 실패 |

P2→P3/seed11 변화: greater 10.5141%, lesser 7.8477%, source별 13.7886%/12.3309%, Poisson 3.6776%, response 36.1694%. 여섯 중 다섯 metric이 5% 초과. 개별 경로 refinement는 최대 약 5.47×10⁻⁶. 현재 주 병목은 ODE 허용오차보다 앙상블 적분·경계 이력·광장 연결. 이 수치 차이는 참 적분 오차의 엄밀한 상계가 아님. [봉인 진행 보고](progress_2026_09_28/README.md).

현재 카드: p3/seed811, 120기존 계산 재사용·재검증 + 새120계산 상한. p2→p3/811와 p3/11↔811 양방향 비교. 전 RF/source·여섯 metric 평가. 하나라도 5% 초과면 전략 검토; 미완료/경로 실패면 INCONCLUSIVE. 통과해도 새 window 설계의 부분 근거. 자동 p4/p5 확대·실패한 v2의 남은960계산 일괄 실행 금지. [인계](NEXT_SESSION.md).

어려운 연구 과제라는 평가는 타당한 성격의 문제. 다만 수학 conjecture와 난이도 순위를 비교할 객관적 근거는 없음. 여기에는 새로운 유도뿐 아니라 충분한 장치 입력·수치 수렴·실험 검증이 함께 필요. 계산량으로 누락 측정을 대체할 수 없음.

## 최근 발표의 정확한 범위

| 사례와 날짜 | 확인한 주장 및 검증 상태 | 방법을 해석할 때의 제한 |
|---|---|---|
| OpenAI non-sofic group, 2026-08-01 | 존재 구성 논문과 Lean artifact 공개. Andreas Thom의 기술 해설과 Fournier-Facio의 후속 연구 존재 | 기존 Kun 및 Kun–Thom 정리 사이의 조건 불일치를 연결한 명제가 핵심. 모든 group-theory 문제를 푸는 일반 기법은 아님 |
| OpenAI Navier–Stokes, 2026-09-08 | 매끄러운 외력이 있는 3차원 비압축성 계의 유한시간 singularity, Clay C/D를 충족한다고 발표. 논문·Lean 공개 | 무외력 문제의 blowup 증명으로 확대 금지. Clay 09-11은 해결된 것으로 보인다고 반응하면서 평가·귀속 절차 진행 명시 |
| OpenAI quasi-Riemann, 2026-09-30 논문·10-06 공개 묶음 | Re(s)>7/8의 고정 zero-free half-plane. 제타·Dirichlet L 및 지정 Hecke L 범위의 formalization 명세 공개 | 전체 RH의 1/2 critical line과 다른 목표. formalization에서 후속 응용은 제외. 10-05 대안 11/12는 더 약한 영역 |
| Anthropic RH 관련, 2026-08-10 | 전체 RH 시도는 실패. Critical line의 simple zeros 비율 하한을 기존41.6%에서 약67.25%로 개선. 회사 발표는67.2%로 표기 | 전체 RH 해결 아님. 독립 Lamzouri 논문은 더 간단한 재증명과 추가 결과 제공 |

Non-sofic 근거: [OpenAI 발표](https://openai.com/index/ten-advances-in-mathematics/), [논문 Chapter 3](https://cdn.openai.com/pdf/ten-proofs-oai.pdf), [Thom의 직접 해설](https://mathoverflow.net/questions/513866/what-are-the-key-new-ideas-in-the-proof-of-nonsoficity-of-groups-in-openai-s-con/), [Fournier-Facio 후속 논문](https://arxiv.org/abs/2608.02025). 공개 Lean을 이번 검토에서 직접 빌드한 것은 아님.

Navier–Stokes 근거: [OpenAI 발표와 실행 방법](https://openai.com/index/navier-stokes-solution/), [Clay 09-11 발표](https://www.claymath.org/news/navier-stokes-announcement/). [Quanta의 Fefferman·Córdoba·Martínez-Zoroa 인터뷰](https://www.quantamagazine.org/ai-has-solved-one-of-maths-1-million-millennium-prize-problems-20260908/)는 선행 인간 연구와 Lean 명제가 원래 문제에 해당하는지 확인할 필요를 설명.

Quasi-Riemann 근거: [정확한 formalization 범위](https://raw.githubusercontent.com/openai/math/main/lean/docs/003.md), [09-30 논문](https://github.com/openai/math/blob/main/preprints/The-Quasi-Riemann-Hypothesis-September-30-2026/paper.pdf). 공개 명세에 근거한 범위 확인이며 전체 증명의 독립 인증은 아님.

Anthropic 근거: [발표·방법](https://www.anthropic.com/research/riemann-zeta), [Lamzouri 독립 재증명](https://arxiv.org/abs/2609.02882). 비율 하한 개선과 zero-free 영역은 서로 다른 결과.

10-06 OpenAI 묶음은 722원고·372결과군, 약4,000문제를 제시한 평가에서 선정. 검증 단계가 다양하고 모든 원고가 Lean formalization을 갖춘 것은 아님. 평균 약3시간 ChatGPT Pro 상당 계산이라는 설명도 **비공개 내부 모델** 기준. 제타 zero-free 연구는 고정 절차의 예외. 이를 개인 계정의 성공률·필요 시간으로 읽으면 안 됨. [공식 공개 설명](https://openai.com/index/sharing-ai-progress-in-mathematics/), [저장소의 방법·제한](https://github.com/openai/math).

## 개인 사용자가 가져올 수 있는 기법

아래 응용은 외부 자료를 바탕으로 한 **프로젝트별 제안**. 해당 기업이 GABES에서 효과를 보였다는 뜻 아님.

| 실제 문서에 나타난 방법 | 가져올 실행 형태 | Grand Challenge 적용 |
|---|---|---|
| 기존 정리의 출력과 다음 정리의 가정 사이 bridge 찾기 | 한 세션에 하나의 빠진 연결 명제 | Moving-atom retarded response/noise → nonlocal Maxwell → 검출 mode 연결 |
| 쉬운 보조문제부터 풀고 유망 경로에 자원 재배분 | 같은 병목을 보존하는 작은 fixture | Thermal 전체 대신 한 경로·두 공간점·한 optical mode의 end-to-end 검산 |
| 접근 다양화 후 중간 결과 교환 | 제안자·독립 유도자·반례 검토자, 최대3역할 | 서로 다른 적분 표현·경계 이력·광장 연결만 비교. 단순 seed 반복과 구분 |
| 실패 분류와 재시도 금지 목록 | 실패 가정·최소 반례·재개 조건을 짧게 저장 | 수렴 실패, 물리 모델 실패, 자료 부족, 이미 아는 항등식 분리 |
| 독립 재증명·의도적으로 틀린 세계의 negative control | 원래 풀이를 보지 않는 검토와 실패해야 하는 fixture | Density 중복·source cross-term 제거·위상 부호 오류·비물리 noise를 검출하는 감사 |
| AI가 코드 후보 생성, 실행 evaluator가 선별 | 후보2개·고정 비용·기준 통과 후 속도 평가 | Quadrature/control variate 후보의 정확도·비용 평가 |
| 복잡한 성공을 짧은 메커니즘으로 다시 설명 | 최소 수식·가정·오차·반례 요약 | 긴 계산에서 재사용할 reduced lemma 추출 |

Bridge 예: Thom은 여러 expander를 주는 정리와 하나의 expander가 필요한 정리의 불일치를 핵심으로 설명. 새 명제가 그 간격을 메움. 변수를 bounded quantity로 바꾸는 장치도 사용. 여기서 가져올 것은 **필요 가정을 충족시키는 연결 문제를 분리하는 사고법**. GABES의 metric 분모를 바꿔 통과시키는 허가는 아님. [Thom 직접 해설](https://mathoverflow.net/questions/513866/what-are-the-key-new-ideas-in-the-proof-of-nonsoficity-of-groups-in-openai-s-con/).

실패·검토 예: Anthropic 공개 부록은 알려진 결과·RH와 동치인 주장·유한 수치 관찰·동어반복을 걸러내고, 서로의 답을 못 보는 검토자와 RH가 거짓인 합성 사례로 숨은 가정을 검사한 과정을 설명. 같은 모델의 합의는 인증이 아니라고 명시. [과정 부록 pp.8–12](https://www-cdn.anthropic.com/d7f3ecf1d01392d887f8bc974ca187e2a121b1ed.pdf).

Code search 예: Google DeepMind AlphaEvolve는 LLM 코드 제안·자동 실행 평가·후보 보존의 반복을 사용. 개인은 작은 후보 집합과 기존 참조 구현으로 이 구조를 실행 가능. 전체 플랫폼·대규모 진화 탐색 재현은 전제하지 않음. [DeepMind 공식 방법](https://deepmind.google/blog/alphaevolve-a-gemini-powered-coding-agent-for-designing-advanced-algorithms/).

전문가 관점: Tao의 공개 대담 해설은 문헌·계산·코드·접근 실패 판별에 AI를 사용하고, 문제 선택·workflow·검증이 중요해졌다고 설명. [대담 해설](https://forum.openai.com/public/blogs/terence-tao-ai-is-ready-for-primetime-in-math-and-theoretical-physics-2026-03-10). [The Verge의 Kun·Fournier-Facio·Maynard 인터뷰](https://www.theverge.com/ai-artificial-intelligence/977273/the-ai-takeover-of-mathematics-has-begun)는 선행 연구 귀속과 성공 사례 비용의 해석 문제를 지적. 검증 가능한 메커니즘과 실제 실패 비용을 함께 기록할 이유.

## 제외할 기업 자원과 과장된 해석

- Navier–Stokes의 약1만 동시 에이전트·해당 문제 약1,300억 출력 토큰·진행 중 내부 모델 학습. [실행 공개](https://openai.com/index/navier-stokes-solution/).
- Anthropic의 3,100만 출력 토큰·약60 subagent를 개인의 기본 실행법으로 복제. [방법 공개](https://www.anthropic.com/research/riemann-zeta).
- 비공개 모델과 전용 강화학습을 구독 계정에서 그대로 사용할 수 있다는 가정.
- 성공한10결과 약2,000달러 환산 또는 평균3시간 상당을 전체 탐색비용·재현 비용·성공 보장으로 사용. 실패 시도와 선정 효과를 분리해야 함.
- “계속해/자신을 믿어”를 검증된 독립 기법으로 취급. 격려의 효과를 분리한 실험은 제시되지 않음.
- OpenAI discovery walkthrough를 원본 전체 실행 로그로 취급. 문서 서두상 다른 모델이 원래 사고 기록과 최종 논문을 읽고 재구성한 서술. [Walkthrough 서두](https://cdn.openai.com/pdf/reasoning-walkthroughs.pdf).

## 이미 갖춘 절차와 실제 추가점

이미 있음: 사전 질문·gate·계산 상한, 독립 QRT/reference, analytic limit, 실패 보존, source-sealed cache, 소과제 종결, 반복 중단, 입력 provenance, 총책임자 검토. 새 framework를 만드는 것만으로 연구 진전 없음.

추가점은 **진단 분할을 새로운 estimator 후보로 바꾸기**, **blind 유도와 악의적 negative control**, **원자→광장 bridge를 작은 문제로 종결**, **수치적으로 최적화하는 계수와 물리 입력의 명확한 분리**. 성공 여부는 같은 비용의 독립 평가로 판단.

## 개선안 1 앙상블 적분의 분산 감소

새 solve 없이 기존 source-bound raw packet과 geometry로 개발용 후보를 검토. 현재 face·residence·velocity 분할은 signed 합산 진단이며, 아직 새 표본 배분 규칙은 아님.

긴 경로만 추가하는 접근은 부족. Response 0.1MHz 변화36.1694%에는 체류2–4μs +19.9332%p, 1–2μs +11.6899%p 기여. Source별 lesser의 최악 atomic_inflow는 0.5μs 미만 경로가 중요. 서로 다른 분할의 기여 합산 금지. [전체 진단](progress_2026_09_28/README.md).

후보 A: 입사면×짧은/중간/긴 체류 또는 transverse-speed의 소수 strata에 표본 배분. 원래 incoming-flux measure의 weight/Jacobian 유지. Geometry bin 빈도는 정확한 층 확률이 아니므로, 해석적 확률 또는 독립 수렴한 확률 적분과 그 오차가 필요. 더 작은 첫 후보는 기존 Sobol 좌표의 고정 구간을 분할해 알려진 mass를 사용하는 것. 중요층만 남기거나 density를 다시 맞추지 않음. 분산 감소의 성공은 미검증.

후보 B: 해석적 control variate. 균일 정상 밀도·ballistic first-exit·닫힌 collection volume Q·동일 incoming-flux measure에서 다음 항등식을 독립 유도:

\[
\int \tau(\ell)\,d\Lambda(\ell)=nV.
\]

각 복소 RF/source 행렬의 unweighted path functional F에 대해 개발 자료에서 고정한 행렬 C를 쓰면:

\[
I[F]=I[F-C\tau]+C\,nV,\qquad
\widehat I_C=\sum_\ell J_\ell(F_\ell-C\tau_\ell)+C\,nV.
\]

C의 단위는 F/time. C는 적분 오차를 줄이는 **수치 계수**, 목표 gain/squeezing에 맞추는 physical coefficient 아님. Unbiasedness는 표본 규칙의 무편향성과 평가 표본에서 독립적으로 고정한 C를 수학적으로 보인 경우에만 주장. 독립 scramble의 spread 감소는 진단이며 미지의 참 적분값·참 오차를 인증하지 않음. Deterministic quadrature에서는 unbiasedness를 자동 주장하지 않으며 원래 적분으로 수렴함을 별도 검사.

이는 nV/추정점유수로 spectrum을 곱하는 사후 재정규화와 다른 식. Density·arrival rate·원래 packet은 그대로. Residual F−Cτ가 PSD일 필요는 없지만 최종 물리 covariance/source closure/uncertainty는 검사 필요. Occupancy 항등식은 해당 ballistic 모델에서만 사용; collision·inhomogeneous density로 자동 확장 금지. 기초 measure·점유수 조건: [inflow 유도](inflow_derivation.md), [ensemble 합산 계약](transport_ensemble_derivation.md).

C도 원래 source 합산·ordering·Hermiticity·sideband conjugacy·phase-selection 항등식을 만족하도록 묶음. Total-source 계수는 source별 계수의 합과 같아야 함. 행렬 성분별 독립 최적화로 관계를 깨뜨리지 않음. 수정된 유한격자 결과에 물리적으로 허용되지 않는 음의 고유값이 나타나면 FAIL/INCONCLUSIVE; clipping 금지.

원래 p3/811 결과는 원래 규약으로 먼저 감사. 후보 B는 별도 numerical contract·출력의 shadow diagnostic으로만 비교. C의 개발 split과 평가 split을 관찰 전에 고정. P3/811은 p2/811의 모든 경로를 포함하므로 p2/811로 C를 정했으면 독립 CV 평가가 아님. 독립 estimator 검증에는 C 고정 뒤 별도 새 seed를 선언. 이미 본 p2/p3 자료는 development이며 untouched로 재명명 금지. 큰 행렬 상쇄를 보존하고 전체 여섯 metric·RF/source를 평가. 캐시 결과에서 개선되어도 새 평가의 5% gate가 자동 통과하는 것은 아님.

## 개선안 2 경계 이력의 작은 반증 실험

현재 artificial boundary의 반폭 약173.205μm, pump waist530μm. 입사면은 pump가 꺼진 먼 곳이 아니며, fresh thermal state reset은 실제 셀에서 검증할 물리 가정. [현행 Rb scope](rb_thermal_ensemble.md), [기존 미실행 history 계획](progress_2026_09_28/boundary_history_plan.md).

대표6경로를 기존 기하 strata로 사전 선정. Collection Q·mode area·density·velocity·arrival weight를 고정하고 pre-entry pumping history만 바꿔 paired 비교. 기존 계획의 R=a,w,2w,3w 계열을 사용하되 zero-history 복원과 pump-off control부터 확인. 상태·mean·총 covariance·복소 response를 비교; noise source ancestry는 history 연장에 맞춰 별도 추적.

복원/독립 reference 실패는 구현 질문; 그 검사가 통과한 뒤 history 차이가 수치 오차보다 크면 **선택 경로의 conditional history sensitivity**. 사전 정의한 history-independence 명제와 threshold를 위반한 경우에만 해당 명제의 반례. Thermal ensemble이나 실험 모델 오류의 크기로 일반화 금지. 차이가 작아도 선정한 경로만 허용.

대표 경로·history별5해상도와 재사용 가능성을 preflight에서 집계. 새90solve 첫 상한 후보는 6경로×3새history(w,2w,3w)×5spec이며, R=a의30기록을 exact-spec으로 재사용할 수 있을 때만 성립. Source/model 변경으로 재사용이 안 되거나 복원·control/reference 요구가 상한을 넘으면 범위를 축소하거나 비용을 다시 선언. 정확한 wall time은 미측정. 다른 모델 identity의 packet을 역사 v2에 재결합하지 않음.

## 개선안 3 원자에서 광장으로 가는 연결 문제 분리

Non-sofic 사례와 가장 직접적으로 닮은 질문: 원자 코드가 주는 retarded response/source-resolved correlation을 nonlocal Maxwell과 collection-mode solver가 요구하는 kernel로 **어떤 가정 아래 변환할 수 있는가?**

한 경로 또는 고정 소수 경로·작은 spatial basis로 다음 전달 계약을 유도: 평균 polarization, causal response kernel K(z,z′;Ω), noise C(z,z′;Ω), ±Ω pairing, photon-flux normalization. Mean response만으로 noise를 결정하지 않음. 같은 원자가 두 위치를 지나는 correlation과 reservoir 중복 계상 검사.

구현자는 kernel discretization, 독립 검토자는 직접 integral equation으로 검산. No-atom·zero-length·passive limit와 물리 temporal mode의 CP/uncertainty 검사. 작은 fixture의 end-to-end PASS는 연결 방식의 범위 있는 증거; 실제 thermal optical 인증은 WP-THERM/BOUNDARY 후속 조건 유지. 기존 [WP-OPTICAL](blueprint.md)의 구체화이며 새 master plan 아님.

## 개선안 4 독립 검증과 실패 기억을 강화

최대3역할: 제안자, 독립 유도/구현자, 반례 검토자. 같은 모델·사람이 시간차로 수행해도 됨. 독립 유도자는 최종 공식을 먼저 보지 않고 입력·관측량·가정만 받음. 결과 공유는 초기 유도 뒤 수행. 역할 분리는 통계적 독립성이나 전문가 인증을 뜻하지 않음.

검토 질문: “결과가 맞는가?”에 더해 “틀린 입력에도 PASS하는가?”

- Density를 두 번 곱한 estimator, Poisson mean outer를 평균 뒤 계산한 estimator, 필요한 cross correlation을 지운 estimator가 negative-control에서 실패해야 함.
- 한 RF/source만 좋은 값으로 숨긴 결과, 임의 PSD clipping, 관측량/sideband 부호를 바꾼 결과를 거절해야 함.
- 내부 불변량을 통과하는 모델도 독립 QRT·조건별 관측량 검사를 통과하는지 확인. 같은 모델의 reservoir consistency와 실제 장치 적합성 구분.

이미 검사되는 경우 재작성하지 않고 coverage 확인 후 빠진 control만 추가. 실패 기록은 접근·가정·최소 반례·재개 조건·결정 영향으로 요약. “공식과 동치인 재서술”과 “새 입력/검증을 닫은 성과”를 구분.

## 개선안 5 입력 부족과 모델 식별성을 먼저 밝힘

WP-INPUT ledger를 짧은 독립 산출물로 종결: 입력별 측정/문헌/가정/target-fitted/unknown, 단위·공분산·사용 위치·획득 경로. 자료 소유자가 없으면 unknown 유지. GC-2를 돌리기 전에 어떤 입력이 없어서 절대 예측이 conditional인지 표시.

같은 주어진 입력으로 허용되는 두 reservoir/history/collection 모델이 서로 다른 S₋를 만들 수 있는지 작은 반례 탐색. 가능하면 해당 입력 집합으로 유일한 절대 예측이 식별되지 않는 것. 더 긴 추론 대신 필요한 독립 측정 한 가지를 지정. 평균 gain 일치만으로 noise law가 결정된다는 가정 금지.

Sensitivity는 목표 trace fitting 없이 독립 입력 uncertainty 범위에서 수행. Numerical error·input uncertainty·model discrepancy를 분리. Holdout 조건과 gain/RF/SQL 비교 정의는 자료 보기 전에 고정. 이미 탐색한 operating point를 untouched로 바꾸지 않음.

## 제안하는 실행 순서와 중단 기준

비용은 **초기 범위 제안**. 기존 card의120상한 외에는 measured runtime 또는 채택된 실행 권한이 아님. 전체 병렬 실행은 전제하지 않음.

| 순서 | 닫을 질문 | 계산 범위 | 판정과 다음 행동 |
|---|---|---|---|
| 1 | 정확한 measure를 유지하면서 predictable part를 뺄 수 있는가? | 기존 cache/geometry, 새 atomic solve0. CV identity 독립 유도와 개발용 비교 | Raw sum·units·source·density 보존 실패면 폐기. 개선도 exploratory로 기록 |
| 2 | 현재 pilot이 quadrature 확대를 지지하는가? | 기존 p3/811 card의120reused+새120상한. 원래 estimator 감사 | 여섯 metric과 양방향 비교 모두≤5%일 때 새 window 설계. 하나라도 초과면 PIVOT_REVIEW |
| 3 | 적분 전략 전환 시 동일 비용에서 새 estimator가 낫나? | 최대2후보·별도 선언·평가 seed·공동120새solve 상한 후보. 표본 수/5spec 비용 preflight 필요 | Raw 관측량 보존·독립 평가와 모든 gate 확인. 보편 가속·참 오차 상계 주장 금지 |
| 병행 A | 기존 fresh boundary가 선택 경로에서 영향이 큰가? | 대표6경로·history 비교·새90solve 첫 상한 후보 | 복원/reference 선행. 유의한 차이는 경계 가정의 제한 반례; 아니면 범위 있는 음성 결과 |
| 병행 B | 원자 packet에서 광장/검출로 계약이 닫히나? | One-path/small basis fixture, 독립2경로 | Limit·normalization·correlation·CP 검산을 통과한 작은 연결 명제 |
| 병행 C | 필요한 독립 입력과 acquisition이 정의됐나? | Atomic solve0, ledger 작성 | Unknown 보존·담당/경로/재개 조건 기록으로 WP-INPUT inventory 종결 |

매 세션 주 질문 하나·후보 최대2개·검토 역할 최대3개를 기본 범위로 제안. 나머지는 순차 진행 가능. 실제 제공되는 usage limit 안에서 사고·wall time·새/reused solve를 기록하고 예산 도달 시 checkpoint. Token 무제한 또는 추가 API 지출이 선행조건은 아님.

성공 metric: 봉인된 독립 근거로 빠진 연결 하나가 닫힘, 가정 하나가 반증됨, 다음 실행의 선택지가 줄어듦. 토큰·코드·보고서 개수와 구분. 공식 네 milestone의 acceptance는 그대로.

## 총책임자 판정

새 근거: 외부 사례의 범위·비용·실제 workflow를 교차 확인; 현행 병목과 연결. 특히 bridge 분리, cheap falsification, blind 재유도, matrix control variate가 구체적 제안으로 도출됨.

판정: 출처·전략 검토 완료. 제안 기법의 GABES 성능 **INCONCLUSIVE**, 새 thermal/optical/실험 인증 없음. 기존 v2 FAIL·공식0/4 유지. 연구 절차의 많은 부분은 이미 있으며, 다음 투자는 작은 판별 실험의 실행에 집중.

두 독립 검토 반영: 외부 결과 범위·반올림, CV의 표본 독립성·계수 제약, 층 확률의 오차, history 비용·반례 범위. [문서 검증 기록](ai_research_strategy_2026_10_07_verification.json): 로컬 링크·기존 기술 본문·연구 일지 원래 prefix·목표/no-fit/milestone·선택 봉인 증거 보존 확인. 코드 검사 결과는 이번 검토의 증거로 새로 주장하지 않음.

다음 단일 과학 질문: **p3/seed811의 원래 estimator 결과가 새 적분 window 투자를 지지하는가, quadrature/경계 전략 전환을 요구하는가?** 저비용 캐시 검토와 입력·광장 연결 설계는 이 질문과 병행 가능. 제안된 새 estimator·history 실험은 각 별도 card를 고정한 뒤 실행.
