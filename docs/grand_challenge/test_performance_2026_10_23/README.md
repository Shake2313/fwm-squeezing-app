# 검사 속도 개선 — 2026-10-23

사용자 요청: “반복 FWM 계산과 GC 자료 복사·검증, G: 파일 접근 중 병목이 큰 곳부터 확인해서, 검증 범위와 정확도를 유지하면서 줄일 수 있는 부분을 최적화하고 전후 시간을 비교해주세요.”

## ⟦AI⟧ 구현·측정 결과

| 범위 | 변경 전 | 변경 후 | 조건 |
|---|---:|---:|---|
| SABES bridge 전체 22검사 | 379.44s | 105.31s | 같은 G: 경로·Python 3.12·BLAS/Numba 8 threads |
| 소스123파일 반복검증 | 3.812s | 0.157s | 동일 C: 사본, 5회 중앙값, manifest 완전 동일 |
| GC native exact-transfer 검사 본체 | 92.62s | 16.08s | 같은 C: 사본, setup 제외; 아래 조건 참고 |

- `tests/test_sabes_bridge.py`: 검출기 조건만 다른 4검사의 동일 Ultra 계산을 모듈 fixture로 공유. 전체 입력 dict 일치 강제, 반환값 deepcopy. Bridge·geometry·readout은 매번 실행. 원래 22개 test 이름·assertion AST 동일; grid·허용오차 변경 없음. 입력이 달라지면 독립 solve 사용을 요구하며 실패.
- `analysis/grand_challenge/source_provenance.py`: 정확한 source bytes+파일명 기준, 최대256개 순수 검증결과 재사용. 매번 실제 파일 resolve·alias 검사·읽기 유지. 동일 내용의 encoding/compile/hash 계산만 재사용; mtime/size만 믿는 캐시 아님. 새 manifest·row는 매번 생성. 변경 전 profile의 compile은 2.430/2.874s. 변경검출·봉인·독립 capsule 검증 유지.
- 검증: SABES22통과, provenance96통과/기존 Windows link권한3skip, thermal-transfer37통과. 합계 **155통과·3skip**. 같은 크기·mtime의 변경 검출과 manifest 변조 격리 추가 확인. 원본 G: provenance도96통과/3skip(5.82s). 전체 회귀 실행·전체 실행시간 측정은 하지 않음.

## ⟦AI⟧ G: 접근·복사 판단

컴파일 캐시가 채워진 상태에서 같은123파일의 **읽기 포함 전체 검증**을3회 비교:
C: 중앙값0.174s, G:4.069s. 매회 manifest 동일. 첫 G:32.144s로 변동 큼.
123파일 로컬 복사118.998s; FWM용103파일 사본174.788s(의존성이 완비되지 않아 검사에는 사용하지 않은 사본). 짧은 검사마다 복사하면 오히려 느려질 수 있어 자동 전체복사 기능은 추가하지 않음.

장시간 반복 검증은 byte대조한 로컬 사본이 유리할 수 있음. 필요한 코드·데이터·설정과 collection 범위를 모두 보존하고, 원본 절대경로 계약 검사는 원본에서 실행해야 함. 이번 변경은 필수 source reread·해시계약·독립 verifier·artifact 복사를 생략하지 않음. 복사 자체의 안전한 대폭 축소는 이번 측정으로 입증되지 않았음.

## ⟦AI⟧ 측정 한계·재현

다른 채팅의 전체 검사와 일부 측정이 동시 실행됨. SABES는 전후 각1회로, OS/Numba cache와 부하가 전체 시간차에 포함됨. 3개의 중복 원자 solve 제거는 구조적으로 확정이나72%를 항상 보장하지 않음.
GC transfer 전측정은 동일 사본의 host provenance 함수를 이전 구현으로 주입한 단일 case, 후측정은37case 모듈 내 해당 case. Fresh capsule과 setup 비용을 통제한 전체 전후 비교는 아님. 반복 소스검증24배를 전체 회귀 배속으로 환산하지 않음.

- [측정값·명령·소스 SHA256](measurements.json).
- 원시 로그·JUnit·assertion 대조: `output/test_performance_20261023/`.
- 소스검증 재현: `python output/test_performance_20261023/gc/benchmark.py <source-root>`; 같은 폴더의 before/after helper 사용, 지정 root의 `gabes/`·`analysis/grand_challenge/` 검사.
- SABES: `python -m pytest -q tests/test_sabes_bridge.py --durations=0`.
- GC: `python -m pytest -q tests/quantum/test_source_provenance.py tests/quantum/test_thermal_campaign_extend.py --durations=20`.

총괄 검토: 반복 검증 비용 감소, 물리 solver·공식GC milestone·과거 봉인자료 변경 없음. 다음 병목은 독립 capsule의 cold start와 G: 실제 읽기이며, 새 profile 없이 추가 캐시/검증생략을 확대하지 않음.
