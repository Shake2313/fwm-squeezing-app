# ⟦AI⟧ PSR 문헌 조사 기록

조사일 2026-10-30. 검색 cutoff 동일. 사용자 요청은 광범위한 PSR 개선 근거 수집; 체계적 문헌고찰의 완전성·모든 최신 논문 포괄을 주장하지 않음.

## 검색·선정

검색어군: polarization self-rotation / vacuum squeezing / polarization squeezing; atomic excess noise / Heisenberg-Langevin / Doppler; cold atoms / buffer gas / wall coating; spatial modes / self-defocusing / orbital angular momentum / multipass; Stokes tomography / entanglement / EIT filtering; low-frequency / magnetometry / magnetic field / seeded / artificial intelligence control / squeezed-light array.

초기 논문과 최근 원문의 참고문헌을 앞뒤로 추적. 검색엔진·arXiv·출판사·저자 대학 저장소 사용. 수치와 모델 주장은 실제 원문 PDF 또는 전체 XML에서 확인. 연구실 소개·검색 snippet은 발견용이며 canonical 수치 출처로 사용하지 않음. 서지 정보는 논문 표제, 출판사 및 arXiv journal reference와 대조; 출판사 등록 메타데이터도 보조 사용.

원문 기반 27편 선정. classical SR, noise-inclusive theory, 실제 squeezing, null result, 공간·검출·응용을 함께 포함. OPO·비축 FWM·spin squeezing 등 별개 메커니즘은 원문 인용에 등장했다는 이유만으로 PSR 실험 행에 추가하지 않음. 비교 가치가 있는 cavity X-system, 일반 공간 전파, NMOR·coherence-array는 `context_only`로 수록.

## 원문 접근

| 논문 | 직접 접근 문제 | 요청한 브라우저 사용 결과 |
|---|---|---|
| Zhao 2025 | 출판사 PDF 403; PMC 자동 수집이 challenge 응답 | PMC 본문·방법·Fig. 1–3 설명 확인. PDF 다운로드는 시간 초과. 공개 Europe PMC 전체 JATS XML 확보·해시 저장. |
| Monsa 2026 | APS 직접 PDF 401; preprint와 출판본 값 다름 | 기관 구독 브라우저로 최종 본문 및 PDF 확보. 출판본 4±0.25 dB / 보정 4.5±0.35 dB로 기록. 이전 preprint 따로 보존. |
| Rochester 2001 | APS 직접 PDF 401 | 기관 구독 브라우저에서 출판본 PDF 확보. classical response·전파 방법 확인. |

각 논문의 `sources`가 실제 열람 판본과 경로의 최종 근거. arXiv PDF의 자동 조판 날짜가 출판 연도와 달라도 논문 연도로 오인하지 않음. 저널 원제와 preprint 원제가 다른 Hsu·Lanning·Josse 등은 원문 판본을 함께 표시.

세부 로그: [실험](research_experiments.json), [공간·잡음·응용](research_modes.json), [기초 문헌 원문 파일](research_foundations_sources.json). 추가 후보를 찾더라도 초록만으로 운전점 수치를 채우지 않음.
