# SABES 장치 질문의 응답 상태

2026-10-21 사용자 Q1–Q26 답변 반영. **기본 snapshot은 2026-10-16**, 이후 변경만 날짜·출처를 붙여 계속 갱신.

기준문서: [Notion README_Squeezed_Light_Magnetometry](https://app.notion.com/p/3fa6cba14fee81c9b447c1d6e36cabae).
출처·추가 검수: [직접 확인 기록](user_confirmations_review.md). 이번 답변으로 해결된 내용을 다시 질문하지 않음. 사소한 세부·기존 자료로 찾을 항목은 보류하고, 중요한 입력만 실제 구현·분석에 필요할 때 질문.

| ID | 상태 | 사용자 답변·이후 처리 |
|---|---|---|
| Q01 | 확정 | 10/16 기준, 지속 갱신 |
| Q02 | 추후 자료 | 최신 테이블 사진 제공 예정. 이번 도식은 SQL/FWM 배치 기준이며 실측 좌표 아님 |
| Q03 | 확정 | NEL03A/M은 etalon3→GT 뒤, 삽입·ON·미최적화. 투과율만 모델에 반영 |
| Q04 | 확정 | PM-0S5-10-PFA-PFA-780-UL로 실제 교체. E8257D→N–SMA connector→SMA cable→EOM. 최신 RF는9/03 기록 |
| Q05 | 대부분 확정 | SRS LDC501, Keysight DSOX1204A, Rigol DSA875, Thorlabs TC200. Heater는 Notion Docs의 자체 제작 회로 |
| Q06 | 위치만 추후 측정 | 교차점은 셀 내부, 정확한 z는 모름. 사용자 추가 측정 예정. 기존530/330µm 반경 유지 |
| Q07 | 확정 | 셀 출사면 끝→D-mirror 약23inch. 출사면에서 겹침 확인·D-mirror에서 probe–pump 간격으로 각도 측정 |
| Q08 | 보류 | Thorlabs plano-convex, 상세 품번은 나중 |
| Q09 | 낮은 우선도 | CL=LJ1477L1 추정. chip/asphere는 Sim 자료 우선. TA 전류·온도 변동, 측정 기반 휴리스틱 허용 |
| Q10 | 확정 | OZ 고출력 fiber 사용 안 함. 현재 일반 APC–APC. 과거 길이/MFD/PM/AR 사양 승계 금지 |
| Q11 | 일부 확정·보류 | PBS 전부 CCM1-PBS252/M. 나머지 Thorlabs, 상세 나중 |
| Q12 | 기존 자료 우선 | 현재 정밀 자료 없음. Gisung Sim 자료 참조. etalon controller TC200 확인 |
| Q13 | 보류 | FWM cell 상세는 나중 |
| Q14 | 확정 | 사용자 제공 이미지대로 SQL/FWM 전환 배치. Notion에 원본 첨부 완료 |
| Q15 | 확정 범위 | 현재 SAS lock 안 함, 사용자가 SAS를 육안으로 보고 수동 정렬. controller LDC501 |
| Q16 | 확정 | PDB450A 회로 동일·eye만S3883. gain10⁵. RF→SA, Monitor+→scope 주로probe, −→scope 주로conjugate |
| Q17 | 낮은 우선도 | 현 단계에서 중요하지 않다는 사용자 판단. 추가 port/환산 질문 보류 |
| Q18 | 확인 | 거의 항상 기존 SA setting 사용. USB 저장, 파일명 정밀 구분 없어 과거 데이터 검색 어려움 |
| Q19 | 확정 | 현재 테이블 미접지. SA 같은 테이블 이동 후 큰 전기잡음 없음. 그린세이프·홈메이드 filter, 과거와 다른 배선 |
| Q20 | 확정 범위 | 자체 제작 Docs 회로. 현장 V–T는 FWM 보드 측정. Hanle 같은 보드이며 비슷할 것이라는 사용자 예상; Hanle 실측 교정 아님 |
| Q21 | 확정 | 센서는 셀 윗면, FWM 자기 차폐 없음, cold spot 미측정. Hanle 기존µ-metal과 구분 |
| Q22 | 보류 | Hanle 상세는 나중 |
| Q23 | 기존 자료 우선 | coil 측정 기록은 랩미팅 자료 우선.4/03과4/28 실드 전후 범위 구분 |
| Q24 | 보류 | Hanle 상세 scan/주파수는 나중. 최신F2→F′1·2사이 방향 유지 |
| Q25 | 정정 완료 | 앞뒤가 바뀌었음. 올바른 값: **앞27.17µW, 뒤26.53µW** |
| Q26 | 유지 | 기존 효율·power 기준 유지. 변경되면 사용자가 알리기로 함 |

## 관련 작업 시작 때 다시 알릴 항목

**Waveform generator 기종:** 그 기기를 제어하거나 coil/주파수 응답·배선을 구체화할 때, “기종은 나중에 확인하기로 했던 항목”임을 알리고 질문. 현재 별도 질문이나 시간 기반 알림을 만들지 않음.

## ⟦AI⟧ 다음 작업

노션 장치 기준 갱신 완료. 후속 SABES 패치는 위 확정 입력으로 시작 가능. 기존 solver의 물리 인증·gain 보정 상태와 구분하며, 정확한 기하나 전기 전달함수를 실제로 필요로 하는 작업에서만 남은 중요 입력을 확인.
