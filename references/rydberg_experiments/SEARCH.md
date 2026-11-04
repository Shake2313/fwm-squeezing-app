# ⟦AI⟧ Rydberg-EIT 문헌 탐색·접근 기록

조사일·cutoff 2026-11-03. 사용자 요청은 실험 reference DB. 검색 시작 뒤 Ju et al. 2026을 당분간 Gold로 삼는 사용자 지시 반영. 검색 범위의 완전성·전 분야 최신 기록 순위는 평가하지 않음.

## 검색·선정

주요 검색어군은 Rydberg EIT electrometry, microwave electric field, bright atomic resonance, Autler-Townes, atomic superheterodyne, homodyne, phase/polarization, photon shot noise, SQL, laser frequency noise, optical depth, temperature, repumping, transit time, flat-top, frequency comb, instantaneous bandwidth, MEMS, off-resonant RF.

2007년 EIT 기초와 2012년 밝은 공명 실험부터 Gold 논문의 인용 문헌·후속 2025–2026 원문을 앞뒤로 추적. Rb·Cs 열증기 중심. 두 감도 이론·OD 광학 최적화·냉각 원자 SQL 비교는 `context_only`. 통신만을 다룬 논문, 형광 전용 판독, EIT가 없는 beam·tweezer 센서는 이번 core에서 제외.

검색엔진은 발견에 사용. Canonical 수치·finding은 확보한 PDF의 해당 위치 확인 후 기록. Crossref DOI 레코드와 arXiv citation metadata는 제목·저자·출판일·판본 연결에 사용. 검색 snippet·기존 내부 실험 요약만으로 실측 데이터를 만들지 않음.

## 원문 접근

- Sedlacek 2012의 Nature 페이지를 브라우저로 열어 구독 제한 확인. 공개 arXiv:1205.4461 원고 확보. 출판 제목과 preprint 제목 차이는 source version에 기록.
- Yang 2024 PRA의 HTTP 원문 요청 실패 뒤 브라우저에서 기관 접근 확인. 출판사 PDF 링크로 내려받아 `sources/yang_2024_pra.pdf` 보관. PDF 9쪽·SHA-256 기록. 기관 로그인 추가 조작 없음.
- Mohapatra 2007, Sedlacek 2013, Kumar 2017 OE, Su 2022는 출판사 직접 파일 접근이 제한되어 공개 arXiv 원고 사용.
- Jing 2020은 저자 대학의 공개 PDF. NIST mixer·repump·AT·comb 원문은 NIST 연구 저장소 공개 PDF. 2025–2026 공개 출판본은 Nature·MDPI 경로와 공개 arXiv 사용.
- 잘못 연결된 arXiv:2402.00718은 fluorescence/transmission 비교 논문으로 확인. 대역폭 논문의 올바른 arXiv:2402.17942와 PDF 제목을 대조해 서지 연결 수정. 형광 비교 논문은 DB에 포함하지 않음.
- 최초 임시 수집 키의 이름은 source 파일명에 남을 수 있음. Canonical 논문 ID·저자·출판연도는 실제 서지로 확인. Dixon 2023, Artusio-Glimpse 2024, Zhang 2023, Yao 2026, Ma 2026, Viray 2026 기록에 반영.

26편 원문 모두 로컬 확보. `sources/`에 PDF와 페이지 구분 추출문 보관. `research_sources.json`에 실제 응답·대체 URL·메타데이터·해시 기록. 다운로드 성공은 자동 원문 검수 판정이 아니며 canonical JSON은 수작업 검수 뒤 생성.

## 수록 판단

보고된 최저 감도만 모아 순위를 만들지 않음. Repump·광학 homodyne·RF LO·비공명 Stark·다중 광자·냉각 원자는 각 장치 조건에서 저장. Carrier 조정 범위와 순간 IF 대역폭 구분. 동일 장치 계열은 잠정 `apparatus_group`으로 연결.

수치로 읽지 않은 그림 점은 수치 표본으로 추가하지 않음. 공동 스캔에 실험 수 증가를 위한 Cartesian 확장 없음. Metadata와 원문이 상충하면 출처 판본·위치 유지. Cai의 직접 응답/잡음 추정 및 Viray의 표/고정 신호 평가 차이는 기록별 정의 보존.
