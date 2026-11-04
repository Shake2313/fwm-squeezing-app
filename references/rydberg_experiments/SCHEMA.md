# Rydberg-EIT 문헌 DB 구조 · version 1

사용자 요청은 실험 레퍼런스 DB 구축. 구조는 FWM·PSR의 논문별 JSON → CSV·통합 JSON·오프라인 HTML 경로 계승. Ju 2026 Gold 채택은 직접 사용자 결정.

## ⟦AI⟧ 원본과 생성물

- 원본 `papers/<paper_id>.json`. UTF-8, 파일명·ID 일치.
- 생성물 `database.json`, `papers.csv`, `experiments.csv`, `observations.csv`, `findings.csv`, `index.html`, `validation.json`, `manifest.json`.
- CSV는 UTF-8 BOM, JSON·HTML은 UTF-8. HTML 내부 데이터로 오프라인 검색 가능.
- `schema.json`은 JSON Schema 2020-12. PSR의 공통 필드 유지, squeezing 전용 검사는 전기장 감도 단위·양수 수치 검사로 교체.

논문 필드는 `schema_version`, `paper_id`, `title`, `authors`, `year`, `publication`, `scope`, `study_type`, `topics`, `sources`, `experiments`, `findings`. `study_type`은 `experiment`, `theory`, `mixed`, `review`. `scope.inclusion`은 `primary`, `context_only`, `excluded`.

`experiments[]`는 `experiment_id`, `label`, `kind`, `conditions`, `results` 포함. `kind`로 운전점·공동 스캔·장치·잡음 측정·이론 계산 구분. 그림만 있는 스캔은 finding으로 연결 가능. 독립 숫자나 실험 표본을 임의 생성하지 않음.

`species`, `line`, `configuration`, `measurement`, `seed_mode`, `reference_role`, `apparatus_group`은 실험 재정의 → 논문 scope 순으로 상속. 원문에서 동위원소를 주소 지정했다는 사실과 셀 동위원소 조성 구분.

## ⟦AI⟧ Gold reference

`reference_policy.json`에 `gold_reference_paper_id`, `origin`, `adopted_on`, `status`, `scope`, 직접 요청문 보존. 현재 `origin=direct_user_request`, `status=active_temporary`. 같은 내용을 Gold 논문의 `reference_designation`과 통합 JSON에 저장. `scope.reference_role=gold_reference`는 단 한 논문에 부여. 검증기에서 일치 검사.

열람본 상단 배너·논문 badge·기본 정렬에 반영. 새 Gold 지정은 사용자 결정 기록과 함께 갱신.

## ⟦AI⟧ 출처

`sources[]`는 `source_id`, `url`, `access`, `version`. 로컬 원문에 `path`, `sha256`, PDF `pages` 포함. `path`는 저장소 내부 POSIX 상대경로. `../` 원문 재사용 허용, 저장소 밖 경로 금지.

`access`는 `full_text`, `metadata_only`, `abstract_only`, `unavailable`. 브라우저에서만 원문을 읽으면 `path=null`, `access_method`, `access_note` 필수. 로컬 snapshot 부재 경고. 초록 수치는 `verification=metadata_only`를 명시. 원문에서 확인한 수치·finding은 해당 원문 `source_id`와 PDF 페이지·절·표·그림 `locator`를 연결.

출판 서지는 `publication`, 수치 검수 판본은 `sources.version`으로 구분. 출판 DOI가 확인되어도 preprint 수치를 출판본 검수값으로 승격하지 않음. `research_sources.json`은 발견·접근 기록. Canonical evidence는 논문 JSON.

## ⟦AI⟧ 수치·단위

공통 필드 예시.

```json
{"quantity":"probe_power","value":6,"unit":"uW","value_type":"reported","basis":"Reported probe power","reference":null,"plane":"not_specified","source_id":"main","locator":"PDF p. 3","evidence":"Reported probe power","notes":""}
```

`value_type`은 `reported`, `approximate`, `range`, `lower_bound`, `upper_bound`, `derived`, `digitized`, `not_reported`. 미보고는 `value=null`. 구간은 `value=null`, `lower`·`upper`. 환산은 `derivation.formula`와 입력 근거. 그림 추출에는 방법·불확도 필요. 모든 수치는 유한값.

`field_sensitivity`, `psn_field_sensitivity`, `sql_field_sensitivity`는 양수 수치와 `/sqrt(Hz)` 포함 단위 필요. 최소 검출장은 `minimum_detectable_field`. PSD를 보고한 신호 감도, 직접 변조 감도, atom SQL과 photon PSN은 `basis`·quantity·실험 기록으로 구분. 같은 논문에서 평가 방법이 다르면 별도 기록.

단위는 원문의 물리적 배율 유지. `normalized_unit`은 `µW/uW`, `µm/um`, `°C/degC`, `%/percent` 표기 별칭만 통일. 예를 들어 79 µV/m/√Hz를 79 nV/cm/√Hz로 바꾸지 않음. `uncertainty`는 인쇄값 보존. 신뢰수준 미보고는 추정하지 않음. `basis`에 peak/RMS·RBW·취득 시간·장 기준면 정의 보존.

상태 명칭·room-temperature 등 범주값은 문자열 가능. 수치 미보고를 0으로 대체하지 않음. `probe_detuning`, `rf_detuning`, `intermediate_frequency`는 서로 다른 quantity. `rf_tuning_range`, optical comb span과 `instantaneous_bandwidth`도 구분.

`findings.claim`은 논문 내용 재서술. `model_implication`은 ⟦AI⟧ 해석을 넣을 때 출처 관측과 분리하는 필드. 본 구축에서는 빈 문자열 유지.

## ⟦AI⟧ 검증

`build_db.py --check`는 schema·ID·원문 SHA-256·PDF 페이지 수·범위·출처 연결 검증 후 생성물의 예상 바이트 비교. 파일을 쓰지 않음. `verify_db.py`는 CSV 모든 필드 왕복, JSON 동일성, stale view 탐지, 잘못된 근거·단위·환산식 거부, Gold 정책 일치 확인.

`manifest.json`은 canonical JSON·문서·생성기·생성물·원문 해시 기록. 자기 자신과 검증 영수증은 순환 해시에서 제외. 구조 PASS와 물리 모델 예측 검증은 별도 범위.
