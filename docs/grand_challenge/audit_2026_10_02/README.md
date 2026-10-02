# Grand Challenge 총괄 감사 — 2026-10-02

**판정: 실제 연구 진전 있음. 무의미한 반복으로 단정할 근거 없음. 종결 단위·의사결정·범위 관리 보강 필요.**
공식 milestone0/4 유지. 이번 감사 원자 solve0회. 현행 설계는 별도 청사진이 아니라
[기존 blueprint](../blueprint.md)의 현황판·작업 단위·현재 결정·총책임자 검토로 편집.

## 감사 방식과 근거

실제 기술/R&D 관리의 아래 방식을 연구 규모에 맞춰 적용. 공식 인증·TRL 등급 부여 아님.

| 원전 | 채택 방식 | 여기 적용 |
|---|---|---|
| [NASA Systems Engineering Handbook, Technical Assessment](https://www.nasa.gov/reference/6-7-technical-assessment/) | 기술 지표·WBS별 평가, review 진입/성공 기준, 근거 있는 교정 결정 | 종결 가능한 WP, 턴 종료 PI 검토, failed gate 뒤 다음 투자 판단 |
| [GAO Technology Readiness Assessment Guide](https://www.gao.gov/products/gao-20-48g) | 성숙도와 이를 입증한 증거를 함께 평가 | 내부 계산·수치 수렴·실제 장치 검증의 별도 상태 |
| [Scrum Guide 2020, Definition of Done / Sprint Review](https://scrumguides.org/scrum-guide.html) | 작은 검증 결과, 명시적 Done, 결과 검토로 계획 수정 | 한 턴 한 질문·검증된 종결·다음 결정; 전체 Scrum 도입 아님 |

세 독립 검토 축: **성과/반복 연혁**, **수치·물리 critical path**, **문서/종결 governance**.
같은 증거를 참조한 검토이므로 독립 실험 재현과 다름. 관리 검토 분담과 독립 QRT 구현의 의미 구별.
읽기 범위: repository orientation, 현행 checklist/handoff, 최초 blueprint, 필요한 유도/연혁,
봉인 ensemble·comparison·path-union 보고의 관련 필드. 전체 원자 cache의 새 native 감사 없음.

## 핵심 발견과 조치

| 발견 | 근거 | 이번 조치 |
|---|---|---|
| 작은 성과가 공식0/4에 가려짐 | implementation record23개, milestone4개·acceptance22항목. 최초 S1 single velocity와 공식 GC-1 hot/noncollinear 수렴의 범위 차이 | 소과제6개 종결 인정; failed declaration 판정1개 종결. GC0/4·원래 기준 유지 |
| master가 초기 초안에 머묾 | blueprint·registry `design_draft_2026-09-09`, 최초 착수표 | 같은 blueprint 앞부분을 현행 총괄 설계로 전환; 기술 본문 보존 |
| v2 PASS 위한 추가 계산은 목표 달성 불가 | `[2,3,4]`의 두 edge 필수; p2→p3/11 실패 | 남은960계산 일괄 실행 대신 작은 진단 카드 |
| 실제 새 경로·독립 seed 결과 있음 | [09-28 보고](../progress_2026_09_28/README.md), [union](../progress_2026_09_28/batch-path-union-verification.json):96paths/480records; p3/11 새120·reused120 | 과거 진전 인정; 경로 PASS를 ensemble PASS로 확대하지 않음 |
| ODE 정확도보다 quadrature 차이가 큼 | primary maxima5.47×10⁻⁶/2.19×10⁻⁷ vs10⁻³; independent refinement1.88×10⁻⁷ vs2×10⁻⁶; ensemble response36.17% | 새 실패 없는 ODE 정밀도 강화 보류; ensemble·경계·광장 연결 집중 |
| 큰 window에 긴 미관측 경로 발생 | 아래 frozen geometry 진단 | p3/811만 다음 pilot; automatic p4/p5 금지 |
| boundary·optical·input/holdout이 계속 후속 | [history 계획](../progress_2026_09_28/boundary_history_plan.md)은 설계만 완료; nonlocal Maxwell·독립 입력 미인증 | 열린 의존성 명시; input ledger·작은 optical 설계 병행 가능 |
| 역사 `next`가 현행으로 읽힐 위험 | gain 인계의 p1/11 나머지11경로 문구 | 역사 날짜·현행 참조 명시 |

필수 실패 수치, reference p2/11 → candidate p3/11:

| Metric | 최대 변화 | 5% gate |
|---|---:|---|
| Greater / Lesser | 10.514126% / 7.847664% | 둘 다 실패 |
| Source-resolved greater / lesser | 13.788644% / 12.330948% | 둘 다 실패 |
| Poisson number | 3.677599% | 통과 |
| Retarded response | 36.169445% | 실패 |

[공동 ensemble 감사](../thermal_campaign_v2/ensemble-p2-three-seeds-p3-s11.json):4/9격자 감사,
평가한 refinement1개 실패, p2 directed scramble6개 모두 실패. 다른5격자는 미계산.
Integrity/산술 성공은 thermal 수렴과 별도. 상대 변화는 참 적분 오차의 엄밀한 상계·confidence interval 아님.

## 새 cheap geometry 진단

[실행 소스](geometry_probe.py), [봉인 결과](geometry_probe.json).
현행 working-tree 물리 구현 대신 **v2 sealed sources.zip의 inflow 코드 직접 import**.
108개 raw source byte·plan seal·환경 대조·declared9grid의 path/rate 목록 일치 확인.
ZIP 추출/수정·원자 cache 읽기 없음. geometry0.194초, 출력 전 전체1.391초; 단일 실행 시간.

| 비교 | 관측 |
|---|---|
| 완료 p3/seed11 | ≥4 μs 경로0개; 최대3.99232 μs |
| 미계산 p3/seed811 | ≥4 μs 경로2개; index20/45,6.01058/5.08527 μs |
| p5 seed11/211/811 | 최대 체류6.27116 /15.25765 /23.89051 μs |
| p14 세 seed geometry 추정 | ≥4 μs arrival fraction1.6656–1.6869%, occupancy/nV 기여8.9362–9.0293% |

범위: inflow quadrature의 표본 기하·가중치·체류시간. **atomic noise share·누락 tail 오차·thermal 수렴 인증 아님**.
Occupancy/nV가1에서 벗어나도 arrival rate 재정규화·density 재곱셈 금지.
긴 경로만으로 실패 설명 금지. 기존 층별 진단은 짧은 경로 inflow noise와 signed 상쇄도 보임.

최초 실행:

```powershell
python -B docs/grand_challenge/audit_2026_10_02/geometry_probe.py
```

재현: `--output`으로 **이 감사 폴더 안의 새 JSON 이름** 지정. 기존 결과 덮어쓰기 거부.
예: `--output docs/grand_challenge/audit_2026_10_02/geometry_recheck.json`.

## 총책임자 결정

다음 질문: **p3/811의 새 독립 경로에서 refinement·scramble 차이가5% 안으로 줄어드는가?**
새120계산 상한, 기존120기록 재검증·재결합. 전 RF/source·6metric, p2→p3/811과 p3/11↔811 양방향 비교.
Candidate-SI-floor ensemble/diagnosis로 카드 판정; standalone coarse-floor 비교는 별도 보존.
Path gate 유지. 통과는 새 window 설계의 부분 근거; 실패는 quadrature/경계 전략 검토; 미완료는 `INCONCLUSIVE`.
새 선언·다음 atomic pilot은 이번에 시작하지 않음. 자세한 카드·종결/반복 규칙은 blueprint 소유.

검증 overhead도 존재: 역사 sampled cache load2.756초 중 source 검사2.441초,
native reads12→7회 개선. 전체 속도·누적 token 절감으로 환산할 자료 없음.
정확한 provenance 검사 유지; 성공한 full audit·render·발행을 이유 없이 반복하지 않도록 운영 규칙 보강.

## 검증과 보존

- [시작 hash 목록](before.json):257개 파일의 선택 baseline. 기술 본문·일지 prefix·checklist의 변경 허용 필드 밖 semantic hash 포함.
- [보존·링크·registry 검증](verification.json): 선택 범위 보존 결과. 전체 cache tree/저장소 inventory 아님.
- [현재 회귀 결과](test_summary.json): `python -m pytest -q` → **2091 passed /3 skipped /0 failed**,803.24초. 현재 dirty working tree 결과; geometry helper는 별도 실행·봉인 검증.
- 과거6개 보고의 gate·명시한 범위로 소과제 종결 인정. 과거 실행의 새 native 재인증 아님.
- 기존 physical source·tests·선택 sealed evidence·공식 milestone와 acceptance 변경 없음. geometry helper·관리 문서 추가/수정.
- 기존 branch `main` 유지. 사용자 staged/unstaged 변경을 포함한 작업트리에 commit·stash·reset 수행 없음.

과거 token 사용 복원 불가. token 총량·진척률·실제 장치 정확도를 추정 수치로 만들지 않음.
