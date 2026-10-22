# 사용자 확인 반영과 Notion 추가 검수

날짜: 2026-10-21. 기본 snapshot: 2026-10-16, 이후 항목별 지속 갱신. 직접 출처: 현재 사용자 Q1–Q26 답변 및 첨부 도식. SABES 구현 변경 전 Notion README 갱신 요청.

| 범위 | 확정·방침 | 보류·조건 |
|---|---|---|
| NEL03A/M | etalon3→GT 뒤, 삽입·ON | 미최적화. 당장은 투과율만 반영; suppression 없음 |
| EOM/RF | PM-0S5-10-PFA-PFA-780-UL로 실제 교체. E8257D→N–SMA adapter→SMA cable→EOM | 9/03 RF 기록을 최신으로 사용. 교체일·포트 실파워 미등록 |
| 기기 | ECDL/TA SRS LDC501, scope DSOX1204A, SA DSA875, etalon TC200, PBS 전부CCM1-PBS252/M | Waveform generator는 관련 작업 시작 때 다시 알려 질문 |
| 광학/기하 | 일반APC–APC fiber로 변경,렌즈Thorlabs plano-convex,셀 내부 교차 | 정확한교차z 추후 측정;렌즈/셀 상세 추후. PM/MFD/길이/AR를 과거OZ에서 승계하지 않음 |
| 거리 | 약23inch=584.2mm는셀출사면끝→D-mirror. 출사면중첩·Dmirror의probe–pump간격으로각도측정 | 셀중심/양twins전체간격으로 바꾸지 않음 |
| TA | 매번전류/온도조절,측정기반휴리스틱허용. CL LJ1477L1추정 | chip/asphere중요도낮음,Sim자료우선. 고정TEC25°C현행화 금지 |
| SAS | lock안함,SAS육안수동정렬 | lockservo구현가정제외 |
| BPD | 원PDB450A회로,eye만S3883. gain1e5;RF→SA,mon+scope주로probe/−주로conj | 종단미상 유지.4MHz는명목,새대역측정아님 |
| 전원 | 테이블미접지,SA같은테이블이동후큰noise없음.그린세이프/홈메이드필터 | 배선기존도면과다름. 접지완료/과거전선그대로 추정금지 |
| Heater | 자체제작 Docs회로.센서셀윗면,FWM자기차폐없음,coldspot미측정 | V–T는FWM실측.Hanle같은보드·비슷한특성은사용자예상,Hanle실측으로승격금지 |
| Hanle power | 앞27.17/뒤26.53µW로순서정정 | 기존반대기록은이력. 새증폭주장없음 |
| 효율 | 기존확인값유지,변경시사용자알림 | Q17낮은우선도.보류사항은필요할때질문 |

## ⟦AI⟧ 1회 추가 검수

- 독립 대조 결과: 현재/역사/미등록 기종 구분,셀출사면거리기준,Heater TP1센서/설정분리,FWM대Hanle교정범위,개조BPD명목대역,4/03·4/28 coil 기록 구분 필요.
- Hanle 최신F2→F′1·2사이 탐색 유지. FWM무차폐를Hanle µ-metal제거로해석하지 않음.
- 10/16직접사용자NF선택도 기존README와대조: 9/17경험식은 역사로보존하고NF/10logNF를현재표기로반영. 새이론유도없음.
- NF근거: `1016 랩미팅 자료 초안 작성`,thread01a13fc0-4a20-7503-8704-0b08ca48da4c,turn01a14025-25db-73d2-9335-b4b75f60b44e,2026-10-16 00:18:58 KST: “F는 일반적인 Fano factor가 아니라 Jaspers 논문의 NF입니다. 저도 NF라고 부르겠습니다.”
- 추가NF근거: thread01a13d78-c95a-70a1-9694-0c6f878a82a9,turn01a1433f-ed5b-7733-a17b-77cfc87e39f3,2026-10-16 14:47:05 KST: “스퀴징 정도를 나타내는 부호는 Jasperse에 따라 NF로 표현하기로 하고, 이 값의 로그 값으로 스퀴징 정도를 표현해주세요. 로그스케일 용어는 Squeezing Degree [dB]로 그대로입니다.”

정밀좌표·지연질문에새마감/알림자동화없음. 해당작업을재개할때확인하도록프로젝트메모리에보존.
