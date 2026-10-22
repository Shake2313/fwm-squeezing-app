# SABES 장치 반영 · 2026-10-21

## 사용자 요청·확정

10/16 장치 상태 기준, 이후 날짜·출처 보존하며 갱신. 조사·사용자 Q1–Q26 답변·제공 SQL/FWM 도식 코드 반영. 사소한 미확인 사항 보류; 중요한 입력도 관련 작업에 필요할 때만 재질문.

근거: [사용자 확인 JSON](../apparatus_research_2026_10_21/apparatus_confirmed_2026_10_16.json), [확인·검수 기록](../apparatus_research_2026_10_21/user_confirmations_review.md), [Notion README](https://app.notion.com/p/3fa6cba14fee81c9b447c1d6e36cabae).

## ⟦AI⟧ 구현·적용 범위

### 장치 snapshot과 조회

[apparatus.py](../../../sabes/apparatus.py)와 [apparatus_20261016.json](../../../sabes/apparatus_20261016.json): 실행용 장치 기록. `baseline_date=2026-10-16`, `confirmed_on=2026-10-21`. 장비·광학 노드·RF/Monitor/제어/온도 연결·운용 사실·출처·보류 질문 분리.

| API | 용도 |
|---|---|
| `default_apparatus()` | 기준 snapshot 로드 |
| `get_equipment(id)`, `get_fact(id)` | 확정·추정·미확인 상태 포함 조회 |
| `equipment_for_node(node_id)` | optic popup 연결. 관련 장치 여러 개일 수 있어 tuple 반환 |
| `equipment_rows()`, `connection_rows(category)`, `fact_rows(category)` | 장치·배선·사실 표 |
| `validate_layouts(layouts)` | 누락·오래된 노드 매핑 검출 |
| `as_dict()`, `to_json()`, `from_dict()`, `load()` | 출처·상태·명시적 null 보존한 저장·갱신 |
| `deferred_questions(context)` | 해당 작업에 필요한 보류 질문만 조회 |

현재 LDC501·DSOX1204A·DSA875·TC200·CCM1-PBS252/M·교체한 EOSPACE PM EOM·PDB450A/S3883 반영. 일반 APC–APC fiber의 길이·MFD·PM·AR는 null 유지. 실제 장비 제어·핀 매핑을 기종 기록만으로 생성하지 않음.

화면은 `?app=sabes`로 접근. 새 세션은 Lab 프로필로 시작하며, `Apparatus & wiring` 탭에서 장비·연결·사실·출처를 조회하고 현재 세션 JSON을 내려받는다. 그림의 optic 또는 `Inspect hardware`에서 부품 정보를 확인한다. 프로필 변경·Reset은 선택한 프로필의 입력·전송계수로 복원하며, optic 편집값은 창을 닫아도 유지한다.

### 10/16 실험 프로필과 논문 프로필

`presets.py`의 `Lab reference · 2026-10-16`: 장치 기준일10/16, 운전점은10/13 FWM 기준 T113°C·pump294.804mW·셀 입사 seed proxy약3.92µW 유지. seed proxy는 검출기 기준35mV를 기존 계수로 환산한 값에 pump-OFF cell 투과0.8046·셀 후 효율0.9709를 보정한 계산값. fEOM3.039GHz·pump/seed 반경530/330µm·QE90.5%·셀 후 광학효율 약97.09% 적용. 10/16 현재 온도의 신규 실측이라는 뜻 아님. 반경은1/e² intensity radius; 정확한 waist z 재측정 아님. `Sim 2025 model reference` 별도 유지. 프로필 출처·선택 상태 결과 내보내기에 보존.

채택 RF 표시11.8dBm은9/03 발생기 기록. EOM 포트 보정 없으므로 effective RF 모델 입력과 별도 유지.

확인된 집계 효율과 per-optic 추정 곱 구분. 현행 fiber 기종 변경이 개별 전달손실의 신규 실측을 뜻하지 않음. TA 전류·온도도 고정 실험값으로 승격하지 않고 측정 기반 휴리스틱 유지.

기록 반경 적용 중에는 telescope 배율 대신530/330µm 입력을 쓴다. 배율을 탐색하려면 두 번째 telescope 렌즈의 `Use recorded cell mode radii`를 해제한다. 개별 검출 광학계 손실을 탐색하려면 BPD의 `Use recorded post-cell efficiency and QE`를 해제한다. 활성 입력·계수·출처·unknown은 `session_export.py`의 strict JSON에 보존된다.

### 광로·검출 기하

- Noise Eater: etalon3→GT 이후, seed output fiber 이전. 삽입·ON·미최적화. 이 위치에서 투과율만 적용; noise suppression 추가하지 않음.
- SQL: pump 우회→HWP→두 번째 PBS(CCM1-PBS252/M)→두 검출 입력 전환 광로 추가. **실제 절차는 같은 BPD에서 SQL/FWM 전환·총 검출 광파워 정합, 이미 사용자 확인됨.** 모델의50:50 coherent 분기는 별도 가정. 우회 경로는 attenuation·readout switching을 아직 구현하지 않아 정합된 SQL 기준을 계산하지 않음; 이를 실제 실험의 정합 미확인으로 해석하지 않음.
- 약23inch=584.2mm는 **셀 출사면 끝→D-mirror**. 출사면에서 겹침 확인하고 D-mirror의 probe–pump 간격으로 각도 계산. 셀 중심 거리·두 twin 전체 간격으로 바꾸지 않음.
- Gaussian 전파의 mode origin은 명목상 셀 중심. 검출기 전파거리는 `d_mode = d_exit + L_cell/2`; 실제 교차 z·waist z는 미측정. 도식 좌표는 schematic이며 실측 좌표 아님.

### 전기·열 기록과 미확인 범위

E8257D→N–SMA adapter→SMA cable→EOM. 11.8 dBm은9/03 발생기 표시값; EOM 포트 RF power·cable loss 미등록. BPD gain switch10⁵ V/A, RF→SA, Monitor+→scope 주로probe, Monitor−→scope 주로conjugate. Scope 채널·종단 미확인; Monitor 전압과 RF 차동 gain 구분.

SAS lock OFF, SAS 육안 수동 정렬. 테이블 미접지; SA를 같은 테이블로 이동한 뒤 큰 종전 전기잡음 부재는 관찰 보고. 그린세이프 멀티탭·홈메이드 filter 확인, 현 배선은 과거 도면과 다름.

Heater는 자체 제작. 셀 윗면 센서, FWM 무자기차폐, cold spot 미측정. 회로도 내부 AD603/LT1210/Pt1000 경로는 `schematic`; 현재 보드 부품 실장 검증 아님. 10/07 FWM 기록의 설정2.262V·TP1약2.227V·OD fitting약113°C는 별도 값. Hanle 같은 보드의 유사 특성은 `expectation`; Hanle 독립 V–T 실측 교정으로 승계하지 않음.

Waveform generator 모델 null·`deferred`. 제어·배선·coil·주파수 응답 작업에 실제 필요할 때 이전 보류 사실 알리고 질문. 자동 시간 알림이나 현 작업 선행 질문 없음.

### 이후 갱신·검증

후속 실측은 원본 출처·측정일·적용 장치·단위·기준면과 함께 추가. `confirmed`·`recorded`·`tentative`·`schematic`·`expectation`·`unknown`·`deferred` 구분 유지. 미확인 값은0이나 과거 기종 사양 대신 null. 장치 재배치 후 `validate_layouts()`와 관련 경로 검증 수행.

**검증 완료:** 전체 `python -m pytest -q` **2570통과·3skip·0실패**,902.49초. 실행 중 코드·설정238파일 해시 불변([receipt](full_pytest_receipt.json), [log](full_pytest.log)). 원자 계산 코드와 기준 artifact74파일도 작업 전 해시와 동일([검사](physics_source_integrity.json)).

[실제 화면 검증](ui_validation.json): Lab/Paper 전환, SQL PBS/flip 도식, 장치 popup, 효율0.9709→0.92 편집·닫기·재개방·JSON 저장, RF 포트/검출기 종단 null 보존 확인. 자동 AppTest에서도 editor 표시/숨김·checkbox·Reset 후 상태 보존과 warning 부재 확인. 투과율0의 dark detector, QE0 거부, post-cell profiler/검출 기하 일치, SQL의 원자 NF 미적용을 회귀검사에 포함. 기존 원자 모델·물리 claim gate·공식 GC 완료 상태는 유지.
