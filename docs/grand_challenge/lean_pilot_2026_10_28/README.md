# Lean 단회 시범 결과 — 2026-10-28

사용자 요청: 간단한 사례로 1회 가동한 뒤 **Grand Challenge 정식 도입 전에 종료·보고**.
원문·사전 범위: [protocol.json](protocol.json).

## ⟦AI⟧ 판정과 검증 범위

**시범 PASS. 정식 도입 미결정.** 유리수 공분산 채널의 일반 항등식 6개와 정확한 수치 사례·반례 7개를 Lean 4.35.0으로 증명했다.
각 정리의 공리 의존성은 `propext, Classical.choice, Quot.sound`만 포함한다.
미완성 증명, 추가 공리, native 평가 공리 없음. 두 거짓 명제의 증명 시도는 예상대로 실패했다.
판정 근거: [최종 receipt](run/receipt.json), [정리·공리 출력](run/proof_audit.log), [잡음 단순합 거부](run/negative_naive.log), [잡음 누락 거부](run/negative_missing.log).

모형은 `C(g,y)(v)=g*v+y`, 합성은 `(g2*g1, g2*y1+y2)`.
g는 **공분산/파워 배율**, 전기장 진폭 배율이 아니다.
무차원 quadrature, hbar=1, 진공 분산 v=1/2. 입력–저장고 및 구간별 저장고의 독립성을 가정한다.
물리적 pure loss는 0≤eta≤1, y=(1−eta)/2. 일반 대수 항등식은 모든 유리수 계수에서 성립하며 물리적 허용 영역을 증명한 것은 아니다.

| 사례 | 추가 잡음 y | 진공 출력 v | 결과 |
|---|---:|---:|---|
| 50% 손실 두 번, 올바른 합성 | 3/8 | 1/2 | 정확한 등식 증명 |
| 앞 구간 잡음을 전달하지 않고 단순합 | 1/2 | 5/8 | 1/2와 다름을 증명; 거짓 등식 거부 |
| 첫 저장고 잡음 누락 | 1/4 | 3/8 | 1/2와 다름을 증명; 거짓 등식 거부 |

추가로 gain=2인 이상 증폭과 50% 손실의 순서를 바꾸면 양쪽 총 gain은 1이지만 추가 잡음이 달라짐을 증명했다.
동일 gain만으로 noise가 정해지지 않는 작은 반례다.
모든 정리는 [Main.lean](../../../analysis/grand_challenge/lean_pilot/Main.lean)에 있다.

## ⟦AI⟧ 실제 비용과 시행착오

| 항목 | 관측값 |
|---|---|
| 공식 Windows 배포본 | 902,262,322 bytes, 약 860.5 MiB |
| 다운로드 | curl 진행 기록 약 92초 |
| 압축 해제 | 시작 후 24.15초 이내 완료 확인; 정밀 종료 시간 아님 |
| 최초 빌드 | 24.12초, 증명 구문 3곳 미완료로 실패 |
| 수정 후 빌드 | 1.68초, PASS |
| 최종 재현 스크립트 | 약 5.22초, PASS |
| 그 안의 정리 재검사 | 1.579초 |
| 오류 대조 2개 | 1.428초 + 1.381초 |

첫 빌드에서 기본 simp의 Rat 덧셈 처리와 True로 단순화된 목표의 rfl 사용을 수정했다.
초기 오류 대조는 decide의 reduction 실패였으므로 **수학적 거부 성공으로 세지 않았다**.
최종 대조는 증명된 witness로 목표를 정규화한 뒤 grind 진단의 정확한 거짓 명제까지 검사한다.
최초 빌드·대조 실패도 [initial_build.log](initial_build.log), [initial_NaiveNoise.log](initial_NaiveNoise.log), [initial_MissingNoise.log](initial_MissingNoise.log)에 보존했다.
보조 중간 로그는 개발 진단이며 최종 판정은 run/ 아래 receipt와 로그가 소유한다.

공식 릴리스 SHA-256 일치 확인. 버전·해시·설치 경로: [toolchain.json](toolchain.json).
도구는 사용자 임시 폴더에 휴대형으로 설치되어 OS 정리 시 사라질 수 있다. 아래 setup 스크립트로 같은 버전을 복원한다.
Mathlib 및 외부 Lean 패키지 의존성은 0개. 전역 PATH·루트 의존성·CI·production 코드는 변경하지 않았다.

## ⟦AI⟧ 효용 판단과 한계

작은 연결 항등식을 일반 변수에 대해 검사하고, 잡음 전달 누락의 구체적 반례를 남기는 보조 도구로 효용 확인.
표준 라이브러리만으로 시작할 수 있었고 준비 후 검사 비용은 작았다.
형식화·오류 진단의 개발 비용은 별도로 발생했다. 이 1회 사례로 복잡한 GC 전체의 생산성 개선까지 확정할 수는 없다.

이번에는 정확한 **유리수·스칼라** 대수만 검증했다.
실수/복소수 행렬, multimode 상관, 완전양성(CP), 원자→광장 미시 연결, Python 구현과의 대응은 검증하지 않았다.
열운동 적분 수렴이나 실험 예측의 타당성도 판정하지 않았다.
따라서 기존 GC milestone·과학 결과는 바뀌지 않는다.

총책임자 검토: 작은 proof/반례 검증 가능성을 확인하여 방법 평가 질문을 종료했다.
정식 도입은 보류하며 추가 실행은 시작하지 않는다.
향후 사용자가 도입을 검토할 때의 다음 질문은 “실제 GC 행렬식 한 개에 연결하는 비용이 얻는 검증 가치보다 작은가?”이다.
이 질문과 다음 후보는 AI 제안이며 채택된 계획이 아니다.

## 재현

저장소 루트의 PowerShell에서:

```powershell
.\analysis\grand_challenge\lean_pilot\setup.ps1
.\analysis\grand_challenge\lean_pilot\run_pilot.ps1
```

setup은 고정된 공식 배포본을 내려받아 SHA-256을 확인하고 휴대형으로 푼다.
run은 build → 직접 증명 재검사/13개 공리 감사 → 거짓 명제 2개 거부 확인 순서다.
의도된 실패는 negative/에만 있으며 기본 빌드에는 포함되지 않는다.
재실행 증거를 보존하려면 run에 새 `-OutputDir`를 지정한다.

변경 영향별 검증: Lean 독립 프로젝트와 기록 링크/해시 검사만 실행.
공유 엔진·수치 입력·pytest discovery 변경이 없어 전체 pytest·GC numerical campaigns는 의도적으로 생략했다.
이번 범위의 미실행 필수 검증은 없다.
