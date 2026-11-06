# Drive 실험 배치 근거

조사일: 2026-10-21 KST. 원본 읽기 전용. PPTX는 ZIP/XML 본문·발표 노트 추출, HTML은 script/style 제외 텍스트 추출. 전체 목록·SHA-256: [manifest](lab_sources_manifest.json). 검색용 원문: [text](lab_sources_text.json). 클라우드 원본 ID: [Drive metadata](drive_source_metadata.json).

주 폴더: [Lab meeting](https://drive.google.com/drive/folders/17a64-NZC7E0gvAhmruBFtUeAT1H_XjZw). 로컬 경로 `G:/내 드라이브/원자광학연구실/Lab meeting`.
추가 이력: `G:/내 드라이브/Docs/Gisung Sim`의 2022–2024 PPTX 38개.
ECDL 제작 참고: `G:/내 드라이브/Docs/ECDL laser 구성요소 및 제작방법_2024.08.02.pptx` 15페이지. 총 색인79파일=사용자PPTX27+HTML4+기준도식1+Sim38+ECDL제작1+Markdown8. 본문/노트2020구획. 텍스트 전량 색인, 관련 항목 선별 대조. 모든 이미지·OLE 데이터의 전량 시각 검수는 아님.

## 현재 장치에 가까운 기록

| 항목 | 기록된 내용 | 근거·시점 | 적용 범위 |
|---|---|---|---|
| 광원 | Littrow ECDL + TA/MOPA. TA 뒤 비구면·실린더 렌즈로 모드 성형 | 20260814 PPTX p.3, 20260821 pp.3–4 | 제작 방식·광로. LD·TA chip 모델은 미기재 |
| ECDL 운전 이력 | 25.550°C, PZT 5.42, 전류143.43mA, λvac794.973±0.001nm | 20260814 p.4, Notion08/11 인용 | 해당일 복원값. 현행 고정 설정 아님 |
| TA 출력 결합 | 정방향 렌즈 뒤 거리5/38/51/68cm의 빔 크기5.460/3.498/2.018/2.809mm. 역방향10/38/51cm에서1.000/0.982/0.984mm | 20260821 p.4 | 크기 정의·축 미명시. 셀 waist로 대입 불가 |
| 콜리메이터 | F260(d3.33mm) 결합 효율 낮았음. F220(d2.1mm) 사용 검토 | 20260821 p.4 | 검토 기록. 교체 완료로 처리 불가 |
| RF 발생기 | Agilent E8257D Analog Signal Generator | [20260904 PPTX](https://docs.google.com/presentation/d/1BLzoITnLuACuMvAMAGGYEadvaoczjhHO/edit) pp.4–5, 본문 실험일09/03 | 원문 모델명 |
| Fiber EOM | EOSPACE PM-0S5-10-PFA-PFA-780-UL | 같은 덱 pp.4–5 | Sim2025의 AZ-0S5 표기와 다름. 교체/오기 여부 미해결 |
| RF 경로 변경 | 출력단 Alcatel/Ferrocom60A74-01, 5.9–8.5GHz 부품 제거. RF18→11.8dBm에서 seed 최적화 | 같은 덱 p.5 | 발생기 표시값. EOM 커넥터의 실측 전달 파워 아님 |
| EOM·필터 power chain | Input10.05mW → EOM0.374mW → etalon1 160µW → etalon2+3 8µW | 20260619 HTML, 실험06/18 | 당시 정렬 상태. etalon2·3 인접 설치로 중간 power 측정 불가 |
| 광학 손실 | FWM셀 전26.0µW, detector20.0×232.0/228.5=20.306µW, ηOpt≈0.781 | 10/13 사용자 정정이 수록된 LAB_CONTEXT.md, [10/12 수정 덱](https://docs.google.com/presentation/d/1zT3oIiFfqsuDCq67ZqH-Rf-9kI5-F6jV/edit) p.3 | Pump OFF, 셀 포함 경로. after-cell만의 효율과 구분 |
| 검출기 | Hamamatsu S3883, 사용자 채택 R795=0.58A/W, QE0.905, 위 경로와 결합η≈0.707 | LAB_CONTEXT.md 10/13 사용자 정정 | 확인된 분석 기준. 다시 미확인으로 되돌리지 않음 |
| SQL/FWM | 같은 셋업·flip mirror 교차, 동일 검출 파워 맞춤, 공통 전기 구성·electronic noise | LAB_CONTEXT.md 10/01 사용자 확인 | 이미 확인된 절차 |
| 가열 | set2.262V=표상120°C, 최종TP1 약2.227V, OD 수동 fitting 약113°C | 10/12 수정 덱 pp.3,10–11; [10/16 덱](https://docs.google.com/presentation/d/1bM-8levUoo28_ZMuILzn776YzeW1h5Ax/edit) pp.3,18 | Setpoint·센서 환산·OD 유효온도 별개 |
| OD 조건 | 9.0µW, 빔 직경1.06mm | 10/12 수정 덱 p.10 | 해당 OD 측정의 입력·모델 범위 |
| Magnetic cell 상태 | 10/06 제거 상태. 10/16 기록에는 재삽입 완료·히터 연결·미가열 | LAB_CONTEXT10/13의10/06 기록, 10/16 덱 p.14·직접 사용자 채팅 | 날짜에 따른 변경. 모순으로 제거하지 않음 |
| Magnetic cell 통과 power | 앞26.53µW, 뒤27.17µW; EOM3.039GHz, seed35mV, 센서1.722V | 10/16 덱 p.14 | 원문값 그대로. 비율>1을 수동 셀 증폭·투과율로 채택하지 않음 |
| Magnetic cell 교체 이력 | 길이19.0→51.5mm, Ne4→50Torr, 새 셀 동위원소 명시 없음 | 20260521 HTML p.1 | 5월 상태. 10월 재삽입 셀과 동일 물체인지는 별도 확인 |

## 직접 확인한 그림과 사진

- `Media & Reference/setup.png`: Part1 도식. ECDL/OI/HWP·SAS 분기 → MOPA/CL/OI/HWP → PBS에서 pump/seed 분기. Seed는 polarizer/HWP/fiber EOM → etalon1 → HWP/두 etalon → GT → Noise Eater → delivery fiber. 이 그림의 Noise Eater는 필터 이후에 있음.
- `Media & Reference/20260619/20260619 스퀴징 셋업 part 1.jpg`: 실제 테이블 사진. 구멍 격자, 여러 회전 마운트·미러·파이버 결합부·enclosure·RF 배선 육안 확인. 촬영각 원근, 가림, 표찰 해상도로 절대좌표·전체 part number 확정 불가.
- `Media & Reference/20260521/셋업_펑.png`: 위 Part1 도식 일부에 사고 위치를 강조한 별표. 실측 도면 아님.
- `Media & Reference/Templetes_PPTX/Intensity Difference Squeezing Setup.pptx` pp.1–2: 두 부분의 편집 가능한 기준 도식. p.2에는 별도 pump-SQL bypass HWP/BS/flip mirror2개, PBS/GTP/beam blocker, iris5724(1.6mm)·5725(2.0mm), 400/150/300/150/100/75mm 렌즈 표기.

⟦AI⟧ 도식은 광로·상대 배치의 근거. 이미지 픽셀·PPT 좌표를 mm 또는 실측 테이블 좌표로 승격하지 않음. 사진의 격자 피치·원점·카메라 보정이 없으므로 좌표 추출은 보류.

## Sim 선배 실험 이력에서 보완한 사항

| 자료 | 확인 내용 | 현행 적용 제한 |
|---|---|---|
| 220923 p.12 | 당시 연구실 SA=Rigol DSA875, KAIST 비교기=Keysight N9020B | 현 SA 모델 확정 근거로 사용 불가 |
| 221118 pp.6,17;221125 pp.6,17 | Multitab/UPS, SA,120V·30V DC supply, Homodyne detector, inductor/noise filter 구성 | 2022 homodyne 전원망. PDB450A 또는 현 heater 전원으로 합치지 않음 |
| 230303 pp.7–9 | EMI filter+inductor+ferrite core noise filter,120V supply 잡음 억제 시도 | 실장·접지 구조의 역사 |
| 230330 pp.2–7 | F220APC 쌍, phase modulator, SA210-5B 10GHz Fabry–Perot(535–820nm), RF3GHz amplitude1V 표기 | 장치 점검 이력. 1V의 peak/RMS·종단 불명 |
| 230428 p.7 | 85Rb pure12.5mm AR-coated, Newport10BC16PC.6(PBS 표기), mu-metal; 요구 waist650/350µm | 초기 설계. 현 실물 모델/waist 확정 아님 |
| 230814 p.23;230901 p.10 | PDB450A의 Si detector→S3883 개선 논의. R795<0.58A/W와 QE~93%가 함께 표기됨 | 현 사용자 채택0.905와 구분. 약93%를 재사용하지 않음 |
| 231006 p.5 | 평행광·4f 변경 후 squeezing 차이 거의 없음. F220 평행광 품질 양호, right-angle mirror 일부 clipping 영향 작다고 기록 | 당시 결과. 모든 clipping 손실의 일반 결론 아님 |
| 240129 p.3 | source branch의 Noise Eater 추가 도식 | 역사적 설치 맥락 |
| 240304 p.2 | T~120°C,angle~0.3°,waist~560/340µm;렌즈400/150/300/150/100/75mm,flip mirror2개 | 최종 논문530/330µm와 날짜·상태 구분 |
| 240412 pp.10–14 | three-successive etalon filtering의 이전/당시 source 도식 변경,filter 위치 비교 | 그림 숫자·최종 통과 개수는 특정 판본/경로와 매칭 필요 |
| 240419 p.4 | 최적 etalon 위치에서 squeezing,filtering peak~0.9GHz 표기 | tuning 상태 이력 |

## ⟦AI⟧ 남는 공백

ECDL 제작자료 pp.4,6–13에는 PZT=Thorlabs PC4QM(BNC), 온도센서=TH10K, TEC15–9pin 케이블 TE+1/TE−3/Sensor7·8, LD9pin Interlock1·2/Cathode5/Anode9, driver LDC501 예시가 있음. Collimator C230220P-B/LT230220P-B 예시 및 2024/07 A230TM-B+ring 사용 기록도 있음. 이는 연구실 제작 기준·예시이며 현재795nm 장치의 실장 부품·핀맵을 확인한 결과로 옮기지 않음. 특히 동자료의1560nm grating600/mm·27.9° 예시를795nm에 대입하지 않음.

실측(x,y,z)·테이블 원점/피치·미러 법선·마운트 높이·각도·렌즈 vertex 간격·파이버 길이/MFD·편광축·RF포트별 케이블/감쇠·현 전원/접지망은 원자료만으로 완결되지 않음. 최신 전경/격자 기준 및 장비 표찰을 기존 도식에 대응시키면 상당 부분 해소 가능. 이는 에이전트의 자료 공백이며, 사용자의 미측정/미실행으로 해석하지 않음.
