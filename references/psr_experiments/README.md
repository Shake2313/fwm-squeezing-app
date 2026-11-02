# PSR · Polarization Squeezing 문헌 DB

사용자 요청일: 2026-10-30. PSR 스킴 개선에 필요한 원문을 폭넓게 확인하고, 기존 [FWM DB](../fwm_experiments/README.md) 구조를 참조해 저장.

**[검색·비교 열람본](index.html)** · [논문 목록 CSV](papers.csv) · [실험·계산 조건 CSV](experiments.csv) · [수치·출처 CSV](observations.csv) · [주요 결과 CSV](findings.csv) · [전체 JSON](database.json)

## ⟦AI⟧ 수록 범위와 열람

2001–2026년 **27편**. 고전적 self-rotation, 초기 squeezing 이론·실험, 재현 실패와 원자 잡음, 냉각 원자, 유한 주파수, 공간 다중 모드·OAM, 다중 통과, Stokes 상태 복원, EIT 전파, 자력계, seed·자기장·장시간 안정화, 코팅 셀의 다중 빔 배열 포함. 관련 메커니즘은 `context_only`; 자유공간 PSR와 cavity·seeded Kerr·확산 결합 배열을 동일 실험으로 취급하지 않음.

논문별 원본은 [papers/](papers/). 42개 실험·계산·장치 기록, 249개 수치·미보고 항목, 53개 주요 결과를 연결. 이 개수는 독립 실험 표본 수가 아님. 동일 장치 그룹·공통 논문의 스캔·이론 사례 포함.

27편 모두 원문 확인. 로컬 보관: 26편 PDF, Zhao 2025는 공개 원문 JATS XML. Monsa의 이전 preprint PDF도 별도 보존. `sources/`의 원문 파일마다 판본·URL·SHA-256·PDF 페이지 수 기록. Zhao의 추가 PMC 브라우저 출처는 로컬 snapshot 없음 경고가 있으나, 같은 논문의 전체 XML은 `main` 출처로 저장됨. 원문 파일은 기존 FWM DB와 동일하게 로컬 보관하고 Git에서는 제외; JSON·CSV·HTML·검증 도구는 추적 가능.

`index.html`은 데이터 내장형. 파일을 직접 열어 오프라인 검색·필터·근거 펼치기 가능. 논문 출처의 페이지 링크는 로컬 PDF로 이동. 웹 DOI·원문 링크만 인터넷 필요. CSV는 Excel 호환 UTF-8 BOM.

## ⟦AI⟧ 읽을 때 지킬 구분

- `squeezing_measured` 등은 **SQL 아래 감소량을 양수 dB**로 기록. GABES 그래프의 음수 잡음 dB와 부호가 다름. Anti-squeezing은 SQL 위 증가량. 원문 고유 정의는 `basis` 확인.
- 직접 측정, 전자잡음 차감, 검출손실 역보정은 별도 quantity. 검출 기준면·분석 주파수·RBW/VBW·LO 조건 보존.
- 다른 온도·파워·주파수의 최적값을 한 운전점으로 합치지 않음. 공동 스캔을 Cartesian 측정점으로 확장하지 않음. 그림 숫자 임의 추출 없음.
- 미보고는 `null`/`not_reported`. 빔 diameter·FWHM·waist·1/e² 반경을 임의 환산하지 않음. 출판본과 preprint 값이 다르면 출처 판본별 기록.
- `findings.claim`은 논문 내용의 재서술. `model_implication`과 아래 개선 방향은 **AI 제안**, 구현 완료나 사용자 채택을 뜻하지 않음.

[스키마](SCHEMA.md) · [검색 범위·접근 기록](SEARCH.md) · [교차 검수](REVIEW.md).

## ⟦AI⟧ 현재 GABES 모델에 대한 개선 방향 — 제안, 미구현

현재 구현은 [모델 문서](../../docs/psr/README.md)의 정적 회전 기울기, 일정한 shear·공통 흡수, 가정한 등방 excess noise에 기반. 아래는 문헌을 현재 근사와 대조한 우선순위이며 새로운 완료 기준이 아님.

| 순서 | 개선 대상 | 원문 근거 | 재사용·검증 방향 |
|---|---|---|---|
| 1 | 분석 주파수와 원자 Langevin 잡음 | [Lezama 2008](papers/lezama_2008_pra.json), [Mikhailov 2009](papers/mikhailov_2009_jmo.json), [Hsu 2006](papers/hsu_2006_pra.json) | 기존 원자 연산자·Zeeman 구조를 활용해 주파수별 응답과 확산 공분산 구성. squeezing·anti-squeezing·관측 실패 영역 함께 대조. 정적 χ만으로 양자 잡음을 결정하지 않음. |
| 2 | 광학 포화·펌프 전파·검출 편광의 손실 | [Rochester 2001](papers/rochester_2001_pra.json), [Ries 2003](papers/ries_2003_pra.json), [Zhang 2017](papers/zhang_2017_pra.json) | 국소 원자 계산을 길이 방향으로 재사용. 밀도 증가와 길이 증가를 같은 OD 효과로 축약하지 않음. 실제 실험 파워를 weak-excitation 근사와 비교. |
| 3 | 빔 모드·초점·LO 중첩 | [Zhang 2013](papers/zhang_2013_ol.json), [Zhang 2016](papers/zhang_2016_pra.json), [Lanning 2017](papers/lanning_2017_pra.json) | 먼저 Gaussian 평균·초점 위치·모드 중첩, 이후 다중 모드 공분산 검토. 단일 scalar η로 모든 공간 효과를 설명하지 않음. |
| 4 | 검출잡음·자기장·seed·안정화 | [Horrom 2012](papers/horrom_2012_pra.json), [Valente 2015](papers/valente_2015_pra.json), [Lima 2024](papers/lima_2024_pra.json), [Zhao 2025](papers/zhao_2025_sciadv.json) | 위상 잡음·전자잡음·RIN·검출손실 분리. seed 유무 및 외부 자기장에 따라 모델 범위 구분. 실험적 보정도 허용하되 측정 근거와 적용 조건 유지. |
| 5 | 원자 이동·벽 코팅·응용 셀 | [Novikova 2015](papers/novikova_2015_pra.json), [Monsa 2026](papers/monsa_2026_prapplied.json), [Wang 2026](papers/wang_2026_lsa.json) | 생성 셀과 센싱 셀의 조건 분리. 코팅·확산·다중 빔 결합은 별도 확장. 최소 잡음 위상과 최대 센싱 SNR이 일치하는지 확인. |

첫 비교 대상은 단순한 단일 통과 D1 실험과 잡음 포함 이론. cavity·seeded·coated-array 결과의 최대 dB를 한 보정 목표로 섞지 않음. 이 자료는 이미 열람한 문헌 풀; 미열람 독립 예측 검증이나 GABES 예측 정확도 인증이 아님.

## ⟦AI⟧ 재생성·검증

저장소 루트에서 실행:

```powershell
python references/psr_experiments/build_db.py
python references/psr_experiments/build_db.py --check
python references/psr_experiments/verify_db.py
```

원본 JSON 수정 후 생성. `--check`는 파일을 쓰지 않고 schema·출처·생성물 일치 확인. [validation.json](validation.json), [validation_detail.json](validation_detail.json), [manifest.json](manifest.json)에 결과·해시 기록. `_curate_foundations.py`는 초기 수작업 큐레이션 재현용이며 일상 갱신은 canonical JSON 편집 후 `build_db.py`만 사용.

검증 범위는 문헌 DB·열람본·CSV 이전 경로. 물리 엔진 변경 없음. 수치 회귀·GC 캠페인은 실행 대상 아님. 레거시 FWM CSV는 [archive 원본](../../archive/fwm_squeezing_paper_parameters.csv)과 [복구 manifest](../../archive/fwm_squeezing_csv_migration_2026-10-30.json) 참조.
