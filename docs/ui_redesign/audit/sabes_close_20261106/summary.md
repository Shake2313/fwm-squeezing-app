# SABES 잔여 작업 커밋 · 2026-11-06

직접 사용자 요청: 작업 트리에 남아 있는 SABES 미커밋 변경을 일괄 커밋·푸시해 현재 작업 상태 마감.

## ⟦AI⟧ 범위·검증

- Part1/Part2 배치·glyph·라벨·검출기 두 입력·Fit/확대/이동·키보드/터치·iframe 재진입 이벤트, 관련 테스트와 디자인 규칙/계획/로그 포함.
- 장치 조사 문서·Drive/Notion/채팅 근거·인벤토리·재생성 스크립트와 기존 테이블 검수 기록 함께 보존. 오래된 audit 결과는 해당 날짜의 역사 기록.
- 실제 파일이 HEAD와 같은 `sabes/noise.py`·`tests/test_sabes_eom_noise.py`는 stale staged deletion만 해제. 파일 삭제 없음.
- 임시 PID 기록 두 개와 브라우저 테스트 다운로드 한 개는 로컬 유지·정확한 경로로 ignore. 다른 작업의 staged/미커밋 변경 보존.
- 현재 작업 트리: `python -m pytest -q tests/test_sabes*.py --durations=10`의 13개 파일을 Python glob으로 전달, **246 통과**, 28.08s.
- 분리된 커밋 후보(HEAD의 공통 엔진/UI + 이번 SABES 변경): 같은 13개 파일 **246 통과**, 38.26s. [검증 기록](validation.json). 문서 링크43개·frontend JavaScript 구문 검사 통과.
- 전체 수치 회귀·GC 검사는 계산식·공유 solver·입력 변경이 없어 생략. 문서 JSON·링크와 소스 구문을 함께 확인.
- 공용 UI 후속 커밋을 보존해 재배치한 뒤 페이지·NF·에디터 연동 37개 추가 통과. 기존 246개 검증 범위와 재배치 기준은 validation.json에 구분.
