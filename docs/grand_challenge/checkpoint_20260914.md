# 2026-09-14 중간 저장

- 기반: `f556af17aec10bd47dbc6144534d0f2642f1e0ae`.
- Claude Fast/Balanced: `0462e98`, `f556af1`. 기존 독립 커밋 유지.
- 이번 분리: Ultra Hermitian solve 최적화 / microscopic FWM 연구 체크포인트.
- 연구 범위: 원자 Langevin 잡음, 광장·검출·공간 전파, characteristic transport, CF4, thermal pilot 6경로.
- 전체 thermal ensemble·nonlocal Maxwell·실험 절대 squeezing 인증: 미완료.
- 최신 수치·범위: [연구 기록](research_log.md), [thermal pilot](rb_thermal_ensemble.md).

## 분리본 검증

- 커밋 후보만 별도 폴더에 추출. `python -m pytest -q` 실행.
- 결과: **1388 passed, 11 failed, 25 errors, 1 skipped**. 244.28초.
- 실패·오류 36건 모두 `ValueError: parent source changed: rerun the primary/reference audit`.
- 영향: `test_adjoint_transport_audit.py` 1건, `test_exponential_transport_audit.py` 10건, `test_rb_thermal_ensemble_audit.py` fixture 25건.
- 원인: 기존 보고서가 실행 당시 전체 소스의 원시 바이트에 결합. 분리본은 미커밋 OD/SAS 변경을 제외하고, 기존 추적 소스의 혼합 줄바꿈도 재현하지 못함.
- 신규 연구·Ultra 증거 파일 269개: Git blob과 작업 원본 바이트 일치. `.gitattributes`로 줄바꿈 변환 방지. `git diff --check` 통과.
- 커밋은 연구 체크포인트. 깨끗한 checkout의 전체 테스트 통과 인증 아님.

## 재개 지점

- 소스 정규화·출처 계약 정리 후 primary/reference audit 재실행. 새 보고서 생성.
- 기존 보고서·해시·성공 플래그 덮어쓰기 금지. 실패를 숨기는 검증 완화 금지.
- 경로·줄바꿈이 다른 checkout에서도 검증 가능한 출처 형식 검토. 물리적 입력·연산 동일성 별도 확인.
- Ultra의 `Q_-h = R_h*` 증명·일반 solver fallback 유지. 출처 해시 문제는 음의 harmonic 체인 solve를 복원할 근거가 아님.
