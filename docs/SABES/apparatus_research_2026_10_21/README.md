# SABES 실제 실험 장치 복원 조사

조사일: **2026-10-21 KST**. 사용자 요청: 옵틱 위치·빔 기하·기종·전기 구성을 실제 테이블에 가깝게 반영하기 위한 자료 조사, 확인 정보 정리, 부족 정보 질문화.

**후속 사용자 확인 반영(10/21):** 기본 snapshot은10/16. [Notion 기준문서](https://app.notion.com/p/3fa6cba14fee81c9b447c1d6e36cabae) 추가 검수·갱신 완료. 현재 기종·파이버·RF 배선·Noise Eater·거리 기준·power 정정은 [직접 확인 기록](user_confirmations_review.md), 기존26질문의 상태는 [응답 상태표](questions.md)가 소유. 아래 최초 조사 표에서 ‘미확인/질문’인 항목은 이 후속 확인이 우선하며, 당시 근거·코드 스냅샷은 역사로 보존.

⟦AI⟧ 후속 요청에 따른 [SABES 코드 반영·검증](../apparatus_patch_2026_10_21/README.md)은 별도 기록. 아래 개수·코드 가정은 패치 전 조사 시점의 스냅샷이다.

**광원→EOM/에탈론→FWM→SQL 전환/검출→Hanle의 연결과 주요 기종은 상당 부분 복원됨. 절대 배치 좌표·현행 RF/전원 포트 연결·보드 revision은 미완.** 최신 직접 장치 상태는 주로 10/16, 최신 외부 피드백은 10/19. 조사일의 실물 상태로 일괄 확정하지 않음. 이 문서는 조사·근거 목록이며, 후속 코드 변경 이력은 위 구현 기록이 소유한다.

- **질문 응답:** [26개 응답 상태](questions.md). 해결·기존 자료 우선·사용자 보류를 구분. 지금 추가 질문 없음.
- **후속 구현 입력:** [10/16 기준 사용자 확인 JSON](apparatus_confirmed_2026_10_16.json). 코드 가정 스냅샷과 분리.
- **상세 근거:** [Drive·사진·이력](drive_evidence.md), [Notion](notion_evidence.md), [사용자 채팅](chat_evidence.md), [SABES·논문·매뉴얼](sabes_paper_evidence.md).
- **조사 시점 구현 전체 목록:** [63노드·61링크·85계수 JSON](sabes_implementation_inventory.json), [개수·출처 요약](sabes_inventory_summary.json). 실제 장치표와 분리된 코드 스냅샷.
- **원본 접근:** [Lab meeting](https://drive.google.com/drive/folders/17a64-NZC7E0gvAhmruBFtUeAT1H_XjZw), [Notion 기준문서](https://app.notion.com/p/3fa6cba14fee81c9b447c1d6e36cabae), [Sim 2025](https://doi.org/10.1038/s41598-025-86479-w).

## 자료 범위와 출처 구분

| 범위 | 실제 열람 |
|---|---|
| Drive 로컬 원문 | **79파일**: 사용자 PPTX 27·HTML 4(15개 회차), 기준 도식 1, Sim 선배 PPTX 38, ECDL 제작 PPTX 1, Markdown 8. 본문·발표 노트 **2020구획** 색인. 관련 항목 선별 대조 |
| 시각 자료 | 실제 20260619 테이블 사진, Part 1/2 기준 도식, 사고 표시 도식, Sim 논문 그림·매뉴얼 관련 페이지 |
| Notion | 프로젝트 관계 Experiment Log **42페이지** 전부 fetch. 별도 10/19 피드백·장치·히터·전원·Tasks·기준문서 추가 확인 |
| 채팅 | 관련 **22개 채팅** 열람. 최근·보관 목록 및 로컬 세션 673개에서 주제 검색. 직접 사용자 발언과 AI 답변 분리 |
| 논문·제품문서 | Sim 2025 전 9쪽, PDB450A 매뉴얼, S3883 데이터시트. Allen 2026은 타 실험실 장비 혼입 여부 대조 |
| 현행 구현 | sabes/·설정·layout·검출·계측기와 nominal 경로 실행. 원자 solver 재계산 없이 장치 입력 비교 |

Notion 첨부 사진·PDF·영상은 파일블록 참조까지만 접근되어 픽셀·영상 본문 미열람. 모든 PPT 이미지·Origin OLE·전체 ChatGPT 과거 대화를 전량 조사한 것은 아님. 원문·SHA-256은 [manifest](lab_sources_manifest.json), [검색용 추출문](lab_sources_text.json), [Drive 원본 ID](drive_source_metadata.json)에 보존.

출처 등급: **사용자 직접 확인 / 날짜 있는 실험기록 / 제작·제품 사양 / 논문 / 역사·재고 / 코드 가정·fit / 미확인**. 노션 기준문서는 Claude의 2차 종합이므로 원문의 ‘확정’ 태그와 원로그를 함께 추적. AI 해석·업데이트 제안은 ⟦AI⟧ 절에 배치.

## 최초 조사에서 복원한 광로와 기종

상세 출처·쪽수·Notion URL·채팅 turn ID는 연결된 근거 문서에 수록.

| 블록 | 확보 정보 | 근거·상태 |
|---|---|---|
| ECDL 내부 | 자작 795 nm Littrow. 사용자 설명의 거울 수 3→**2개** 정정. LD·collimation·grating·PZT/TEC 구성 | 채팅 CHAT01–03, 제작자료. ‘거울’ 2개의 개별 역할은 도식 대응 필요 |
| Laser control | **SRS LDC-502**(3/26). 8/11 ECDL 25.550°C/143.43 mA/PZT 5.42/794.973 nm(vac) | Notion. 날짜별 운전점 |
| TA/MOPA | 앞·뒤 비구면 렌즈, 출력 CL 원통렌즈. CL=f70 mm, **LJ1477L1은 추정**. TA 25°C, 1 A→355 mW·2 A→약 1.34 W(8/11) | 채팅 CHAT02·Notion 8/11–12. TA chip·비구면 품번 미확보 |
| Pump fiber 입력 | 8/12 F230APC-780→**8/26 F260APC-780 교체** | Notion 8/24 페이지의 8/26 추가. 원래 도식 F230과 날짜 구분 |
| Pump fiber | **OZ Optics SMJ-A3AHPC3A-780-5/125-3AS-2 AR2**, SM/non-PM, MFD 5 µm, 2 m, 양단 AR, 입력 HPC air-gap FC/APC | Notion 8/12. non-PM 편광 drift 기록 |
| Pump coupling | 8/26 1083→320 mW=29.5%. 7월 1.4 W→230 mW는 이전 상태 | Notion·직접 채팅. 코드 기본 70%와 구분 |
| Seed RF | **Agilent/Keysight E8257D PSG**, 9/03 RF 18→11.8 dBm 최적화 | [9/04 덱 pp.4–5](https://docs.google.com/presentation/d/1BLzoITnLuACuMvAMAGGYEadvaoczjhHO/edit) |
| Fiber EOM | 현 계보 슬라이드 **EOSPACE PM-0S5-10-PFA-PFA-780-UL**, 논문 **AZ-0S5-10-PFA-PFA-780-UL** | 모델 표기 충돌. 명판·교체 여부 질문 Q04 |
| RF 출력단 변경 | **Alcatel/Ferrocom 60A74-01(5.9–8.5 GHz)** 제거 후 11.8 dBm 최적화 | 같은 덱 p.5. 발생기 표시와 EOM 포트 실파워 별개 |
| Seed filtering | etalon 3매 직렬. 논문 FSR 17 GHz/FWHM 950 MHz/peak T85%, 사용 기록 약 0.8 | 제조사·품번·개별 실측 스캔 미확보. PID P200/I5/D6(5/21) |
| Noise Eater | **Thorlabs NEL03A/M**, 도식상 3 etalon/GT 뒤·delivery fiber 앞. 1/07 T92.8%, 5/28 제외, 9/09 상태 불충분 | 현 삽입/가동 여부 미확인. 2 MHz 봉우리 제거 효과도 확정 안 됨 |
| Part 2 pump | collimator 5725(2.0 mm)→QWP→HWP→f300/f150→GTP→합류 PBS | 도식·기준문서. 5725는 전체 품번 아님 |
| Part 2 seed | collimator 5724(1.6 mm)→QWP→HWP→f400/f150→GTP→합류 PBS | 같은 근거. 5724는 전체 품번 아님 |
| FWM cell | **Pure 85Rb, L12.5 mm**, 미량 87Rb 가능성(~0.1%)은 사용자 추정. 직교 선편광 | 직접 채팅·기준문서. 인증 순도 수치 아님 |
| 셀 기하 | **pump 530 µm/seed 330 µm=1/e² 반경**, 교차각 실험 약 0.3°/기본 0.32° | 사용자 확인·당시 GABES default. 좌표·waist 위치·오차 별도 |
| 주파수 | OPD≈+1 GHz, 85Rb D1 F2→F′3 기준 blue. **seed=pump−fEOM** | 최신 직접 확인. 실험 TPD=fEOM−3036 MHz |
| Pump 제거/검출 | 셀 후 PBS→beam blocker, D-mirror/iris로 twins 분리, probe f75/conj f100→BPD | 기준 도식. 미러·iris·렌즈 품번/좌표 미확보 |
| SQL 경로 | pump 우회→HWP/BS→두 flip mirror→같은 BPD 양 eye | 도식+프로토콜. 동일 총 검출 파워·공통 전기 구성 확정 |
| BPD | **Thorlabs PDB450A, 두 PD 모두 Hamamatsu S3883 교체** | Notion 10/13. 순정 NEP/대역을 개조품 성능으로 자동 채택 불가 |
| 보조 계측 | 파장계 **Advantest Q8326**, HI780 fiber 입력, 공용·상시 미연결 | Notion 8/11. Keysight laser spectrum analyzer는 모델 불명 |
| SA/scope 후보 | 2025 재고표 **RIGOL RSA5065N / Keysight DSOX1204A / EDU33212A** | 역사적 재고. 현재 EOM 셋업 각 기기와 동일한지는 미확인 |

각 블록의 전후 연결은 확보. 모든 steering mirror와 실물 ID의 1:1 대응은 미완. 자세한 source branch 및 ECDL 제작·RF 재고는 근거 문서에 보존.

## 전기·제어·검출 연결

| 계통 | 확보 연결·설정 | 남은 범위 |
|---|---|---|
| RF 생성 | E8257D→(9/03 부적합 5.9–8.5 GHz 소자 제거)→EOSPACE EOM | 현재 cable·attenuator·isolator·filter·종단·포트별 power |
| Laser scan | SAS 피크 위치로 OPD 맞춤. 9/30 scope CH1=SAS, CH2=trigger, CH3/4=FWM 출력 | 현재 lock/manual 상태·servo/scan controller·trigger 분기 |
| BPD | RF 출력→SA, 개별 monitor→scope 계통 | 실제 gain·두 monitor 선·50 Ω/high-Z·cable loss |
| SA 운전 | RBW 30 kHz, VBW 300 Hz, sweep 2 s, attenuation 10 dB, span 0–3/4 MHz가 기준문서 설정 | 실험별 detector·평균·ENBW·trace와 설정 대응 |
| 접지 이력 | SA·detector가 다른 테이블. 두 테이블 wire 연결 후 저주파 잡음 소멸(7월 직접 발언) | 현재 접점·chassis/BNC shield·보호접지·전원 분배 |
| Heater 도면 | 1 MHz 발진기 U5/J1 외부 입력→AD603→LT1210→DC 차단 2 µF→16.5 Ω heater | 도면과 현재 교체 보드의 동일성 |
| Heater feedback | Pt1000→TLC2272(U3B)→TP1→U3A/설정 TP2→AD603 GPOS, ±16 V→±5 V | 센서 실부착·실장 R/C·발진원·전원·board ID |
| Heater 환산 | 도면 120°C=TP1 1.566 V, 현장표 120°C=2.262 V, 최종 약 2.227 V | 서로 다른 환산. OD 약 113°C와 구분 |
| 역사 RF 재고 | ZFBT-4R2GW+·ZX73-2500+·FK-3000+·ZHL-1-2W-S+·ZHL-42W+·JFW 50DR-129·MACOM HH-108-BNC | SU(1,1)/AOM 계보 재고. 현재 EOM 회로 전체 장착으로 처리 불가 |

[Heater 회로도 페이지](https://app.notion.com/p/24f6cba14fee804eb94eea4d11b3652e)에 원본 PDF/PNG 참조. 첨부 픽셀은 미열람이며 위 회로는 페이지 텍스트에 근거. ECDL 제작자료의 PC4QM PZT·TH10K·LDC501 예시·LD/TEC 핀맵도 확보했으나, 현재 LDC-502 장치 실장표와 구분.

## Hanle와 날짜별 변경

**87Rb magnetic cell은 51.5 mm·Ne 50 Torr가 사용자 확인 사양.** 주문제작 원통 µ-metal 3겹(약 지름 19 cm×길이 40 cm), 3축 Helmholtz의 z=광축 기록. 4/28 실드 입구 z0, 코일 응답 최대 z20–23 cm. 코일 반경·권수·간격과 현행 전류→B 보정행렬 미확보.

| 시점 | 상태 |
|---|---|
| 5/21 | 구 19 mm/Ne 4 Torr→신 51.5 mm/Ne 50 Torr 교체 이력 |
| 10/06 | FWM power chain 측정 때 Hanle 셀 제거 |
| 10/14·16 사용자 | 탐색 전이 **87Rb D1 F2→F′1·2 사이**로 변경. 오래된 Task의 F1 목표와 구분 |
| 10/16 | 셀 재삽입·히터 연결 완료, **미가열**. EOM 3.039 GHz, seed 35 mV, 온도계 1.722 V |
| 10/16 power | 셀 앞/뒤 26.53/27.17 µW 원문. 순서·변동·offset 질문 필요. 투과>1 모델로 채택 안 함 |
| 10/19 방문 | RF 필터·셀 출사면 가까운 교차·셀 교체/이동은 타 연구실 운용/추천. 사용자 테이블 채택 사실 미확보 |

## 운전점과 효율의 기준면

| 항목 | 보존할 값·뜻 |
|---|---|
| Sim 2025 기준 | 85Rb 121°C, pump 600 mW, seed 8 µW(최고 IDS)/10 µW(setup), 530/330 µm, 0.32°. 현재 운전점과 구분 |
| 10/13 gain run | 113°C, OPD 약 +1 GHz, seed 35 mV **pump OFF·셀 후 검출기 기준**. 양팔 같은 전압 환산·기준면 사용자 확인 |
| 35 mV 환산 | 3.06346 µW는 기존 계수 환산. 창 포함 pump OFF 투과 0.8046과 after-cell 0.9709로 나눈 입사 proxy 약 3.92 µW도 계산값 |
| 사용자 선택 QE | R795=0.58 A/W, ηQE=0.905. 제품표의 0.58은 780 nm 점사양. 795 nm 값은 사용자 채택 근사/기준 |
| 셀 포함 optical 효율 | 26.0→20.306 µW, ηOpt≈0.781, QE 결합 η≈0.707. 10/06 power chain·10/13 확인 |
| 셀 후만의 효율 | 20.0/20.6≈0.9709. 위 0.781과 대체 관계 아님. 양팔 공통 적용 사용자 확인 |
| BPD 부하 | 매뉴얼 10⁵ V/A는 high-Z, 50 Ω 종단에서는 절반. monitor 전압 환산과 RF 차동 gain 별개 |
| 논문 loss | optical 5.5(2)%+system 8.0(2)%=total 13.5(2)% 서술. PDB450A·S3883·QE 90.47% 직접 기재 없음 |

## ⟦AI⟧ SABES 업데이트에 반영할 차이

채택 상태: **후속 구현 제안**. 실제 장치에 대한 새 단정이나 추가 사용자 승인 조건 아님.

| 순서 | 구현에서 확인한 차이 | 필요한 업데이트 |
|---|---|---|
| 1 | 63노드·61링크, 도식 좌표. 85계수: nominal47/fitted13/lab10/paper8/datasheet7 | 실측 mm(x,y,z) 분리, 기종·serial·기준면·날짜·출처·채택 상태 추가 |
| 2 | Part 1의 36링크 전부 물리 길이 없음. Part 2는 25중 14링크에 길이 | 부족 길이 null 보존. design/nominal을 실측으로 승격하지 않음 |
| 3 | 도식의 SQL bypass HWP·BS·flip mirror 2개가 현 63노드에 없음 | 확정 전환 프로토콜을 광로 분기로 표시 |
| 4 | Noise Eater 도식 위치는 etalon/GT 뒤, 코드 손실은 EOM 앞. 렌더 노드 연결도 없음 | 위치 정합, inserted/active/bypass 분리. 최신 상태는 Q03 연결 |
| 5 | CL 설명의 출입 방향 불일치, TA·편광·fiber 전달 단순화 | TA 앞뒤 비구면/CL 등록, non-PM 편광·측정 profile 보존 |
| 6 | RF/DC/monitor/trigger/temperature 배선 graph 없음 | 기기별 port, impedance·termination·cable loss·channel·ground 구분 |
| 7 | 기본 RF 18 dBm, Vπ 4.4 V 가정. 9/03 표시 11.8 dBm 이력 | 발생기 표시/전달 RF/변조 깊이 분리. PM/AZ 표기 충돌 보존 |
| 8 | 명목 경로 waist 실행값 540/324 µm, 문헌·사용자 기준 530/330 µm와 차이 | 공통 1.08 보정과 실측 w0,z0,M² 분리. 기존 확인 취소 아님 |
| 9 | lab EOM loss 5.3 dB는 역산, etalon leak fit 0.076497, verified에 fit 포함 | 독립 측정/fit/문헌/사양과 검증 범위 분리 |
| 10 | 후단 94.5%·QE 90.45%, residual 1.004074로 문헌 loss 정합 | 실험 snapshot별 기준면 효율 분리. 원자 셀 손실 중복 방지 |
| 11 | PD 자동 focus, SA 4 MHz 단일 pole 합성, stock NEP 재사용 | 실제 PD 거리·spot·개조 revision·측정 전달함수·SA setting 반영 |
| 12 | Hanle 현장 사양·전이·코일/실드·heater 상태를 담는 장치층 부족 | probe→Hanle와 conjugate reference 분기, 온도·B 교정·mode 상태 추가 |

권장 장치 필드: id, category, manufacturer, model, serial, installed, active, position_mm, orientation, ports, parameters, units, reference_plane, measured_at, source, origin, uncertainty, adoption_status.

광선: from/to port, path_mm, k_vector, frequency, polarization, power_plane, qx/qy, beam_size_definition, aperture. 전기 연결: from/to port, connector, impedance, termination, cable, loss/transfer, ground/shield.

절대좌표가 없어도 부품·기종·광로·전기 port·출처를 먼저 등록 가능. 좌표는 빈값으로 두고 기준 도식 보기와 실측 보기 분리. 문헌·과거 구성·현장 구성을 서로 덮어쓰지 않음.

## ⟦AI⟧ 총책임자 검토와 다음 질문

실제 장치 입력→같은 원자/잡음 모델→광전파→검출 경로에서 이번 조사는 **WP-INPUT/WP-OPTICAL의 입력 근거**를 보완. 새 gain/noise 물리 인증·공식 milestone 완료 판정 없음. 기존 봉인 계산 반복 불필요.

10/21 후속 확인으로 기본 snapshot은 **10/16**이며 주요 기종·RF 직접 연결·Noise Eater 상태·거리 기준면은 해결. [응답 상태](questions.md)와 [확인 JSON](apparatus_confirmed_2026_10_16.json)을 후속 SABES 입력으로 사용. 세부 좌표·기종은 기존 자료 우선/보류, 중요한 값만 필요할 때 질문. Waveform generator 기종은 관련 작업 시작 때 사용자에게 다시 알림.

검증: [receipt](validation.json). 원본79파일과 SABES 입력 소스6개의 SHA-256 변화 없음. JSON·문서 링크 점검 통과. 저장소 필수 `python -m pytest -q`: **2494통과·3skip**, 903.39초. 11개 경고는 기존 `datetime.utcnow()` 폐기 예고. 장치 실측·시뮬레이션 물리 인증과 별도인 소프트웨어 회귀검사.
