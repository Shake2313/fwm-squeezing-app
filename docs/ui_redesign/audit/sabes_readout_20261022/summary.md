# SABES 상단 요약 · 2026-10-22

직접 사용자 요청: 잡음 수준·Gain G_s/G_c·Detector clearance 세 카드. 나머지는 상세 탭·작은 상태 배지로 이동, 잘못된 라벨 수정.

## ⟦AI⟧ 반영·화면 검수

| 항목 | 결과 |
|---|---|
| 상단 카드 | 9→3. Noise level [dB], Gain G_s / G_c, Detector clearance |
| 모델 상태 | 제목 옆 공유 console 상태 배지. Detection에도 값 보존 |
| Detection | 계산 NF [dB]·NF 선형값·algebraic floor·RIN 진단과 해석 범위 |
| Power budget | 잔류 carrier : seed, quantized 입력 기준과 적용 범위 |
| 라벨 | G_s는 probe/seed gain. Pump gain·Squeezing Degree 상단 라벨 제거 |
| 1440×900 | 세 카드 한 줄. 카드 높이 76.55px |
| 390×844·320×780 | 세 카드 세로 배치. 카드 잘림·문서 가로 overflow 0 |
| 720×450 | 좁은 viewport 검수. 실제 브라우저 zoom 조작 아님 |
| 접힌 사이드바 | 열기 버튼과 제목 분리, 상태 배지 줄바꿈 확인 |

Lab 표시값 −8.16 dB·13.49/12.68·7.0 dB, NF 0.15266·floor −9.16 dB·carrier 0.051%. 사용자 제공 기존 화면과 동일.
Sim 프로필에서 −7.82 dB·372.96/375.04·24.6 dB로 갱신, Lab 복귀 확인. Sim의 더 긴 Gain 값도 1024/900/800/720/390/320px에서 카드 overflow 0.

브라우저 console 오류 0. 기존 Streamlit 기능·테마 경고 존재. 계산식·입력·solver 경로 변경 없음.

- [화면 수치](metrics.json), [Sim 프로필](paper_profile.json), [Sim 폭별 확인](paper_widths.json)
- 화면 PNG: `output/playwright/sabes_readout_20261022/`의 `desktop.png`, `mobile.png`, `small.png`, `readout.png`, `detection.png`, `power_budget.png`. PNG는 기존 git 제외 정책 유지.
- 최종 관련 검사: `python -m pytest -q tests/test_sabes_page.py tests/test_sabes_current_nf.py tests/test_docs_consistency.py` → 31 통과, 4.34s.
- 전체 `python -m pytest -q`: **2707 통과·3 skip·13 warnings**, 1191.07s. [Receipt](full_pytest_receipt.json)·[로그](full_pytest.log).

전체 검사 시작 후 사이드바 열기 버튼과 제목을 분리하는 CSS 여백만 추가. 위 31개 검사·브라우저 검수는 최종 CSS 기준. 소스 변경 내역은 receipt로 구분.
검사 종료코드 0. 수집 스크립트의 종료 후 UTF-8 읽기만 Windows cp949 로그에서 실패했으며, 로그를 UTF-8로 정규화해 보존.

## ⟦AI⟧ 2026-10-23 추가 검수

상단 변경의 커밋 전 리뷰·focused 검증은 [추가 검수 기록](../sabes_readout_review_20261023/summary.md) 참조. 위 전체 회귀 수치는 10/22 당시 결과.
