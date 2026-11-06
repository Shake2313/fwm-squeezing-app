# ⟦AI⟧ Notion 셋업 근거 — 2026-10-21

요청: 현행 SABES 를 실제 테이블의 옵틱 위치·빔 기하·기종·전기 구성에 맞춰 확장하기 위한 자료 조사. 이 파일은 Notion 근거 편집본. 새 실험 결과·사용자 채택 결정 아님. 근거의 실험 날짜와 문서 수정일 구분.

## 자료 범위·판독 한계

- [Squeezed Light Magnetometry](https://app.notion.com/p/2da6cba14fee804d9d04f343edca4f73) 관계 필드의 42개 Experiment Log 페이지 전부 fetch. 2025-11-28–2026-10-14 기록, 템플릿 포함. 최신 별도 검색으로 2026-10-19 방문 피드백 추가.
- [README_Squeezed_Light_Magnetometry](https://app.notion.com/p/3fa6cba14fee81c9b447c1d6e36cabae): 문서 자체 최종 갱신 2026-10-16, 도구 수정일 2026-10-15T19:11:21Z. Claude 가 41개 로그·Tasks 2개·Lab memo·Drive·대화를 합친 **2차 기준문서**. 원문의 ‘확정’은 사용자 확인을 옮겼다는 태그. 원문에 없는 인간 저작·추가 채택을 추정하지 않음.
- SU(1,1) Interferometry → [중요] 셋업 → 구성요소·전원/필터·영상·PID 문서, 2025년 AOM/locking·전자 구성 관련 로그도 확인.
- 검색어 묶음: SABES; squeezing/FWM/Rb85/optical experimental setup; 레이저/광학/셋업/EOM/AOM; 202610 실험/에탈론; 장비/모델/PDA/DBM/RF amplifier; E8257D/EOSPACE/NEL03A; LDC-502/Moku/SAS/PID/locking; RF dBm; 히터/코일/회로도/전원; F260APC; 20260812/20260826/20260428; probe/렌즈/PBS/GTP/초점거리. AI 검색의 순위·상위 결과 제약 있어 워크스페이스 전수조사로 표현하지 않음.
- Notion fetch 의 사진·PDF·영상은 `notion-file-block://` 참조만 확인. 첨부 픽셀/영상 자체는 이 조사에서 미열람. 이미지에서만 읽을 수 있는 라벨·위치·회로 핀은 보류. 외부 제조사 링크 사양을 이번 조사에서 독립 검증하지 않음.

## 1. 현재 FWM 광로 — 문서상 복원

출처: [기준문서 §§0,2](https://app.notion.com/p/3fa6cba14fee81c9b447c1d6e36cabae). 옵틱 거리·초점거리는 Drive `calibration_sim2025.json` lab 항목에서 옮긴 값. 최종 좌표는 문서에 없음. 아래는 해당 2차 문서의 기능상 경로 요약. TA 이후 분기 위치·PBS 순서·noise eater 상태를 현 테이블에서 직접 추적한 결과는 아님. 실제 부품별 광로는 Drive 도안·사진 및 최신 채팅과 대조해 채택.

```text
ECDL(795 nm, 자작 Littrow)
→ OI → HWP → glass(SAS pickoff) → M → HWP → M → M
→ TA → 출력 비구면 → 원통렌즈 CL(f70 mm)
→ F260APC-780 → pump SM fiber(OZ Optics, 2 m, non-PM)
→ collimator(1/e² 직경 2.0 mm) → QWP → HWP
→ telescope(f300/f150) → GTP → PBS(시드 합류)
→ 85Rb FWM cell(L12.5 mm, pump/seed 직교 선편광, 교차각 약0.3°)
→ PBS(pump 제거) → D-mirror/iris(probe·conjugate 분리)
→ probe lens f75 / conjugate lens f100
→ flip-mirror 경로 선택 → PDB450A(S3883 개조)
→ RF 출력: SA / monitor 출력: oscilloscope

Seed:
pump PBS pickoff → 선형편광자 → EOSPACE fiber EOM
→ etalon1 → etalon2 → etalon3 → GTP
→ NEL03A/M → 출력 fiber
→ collimator(1/e² 직경1.6 mm) → QWP → HWP
→ telescope(f400/f150) → GTP → 합류 PBS
→ 같은 FWM cell

SQL:
첫 PBS 에서 분기한 pump coherent 광(원자셀 미통과)
→ 두 번째 PBS 로 두 팔 분리 → flip mirror 전환 → 같은 detector
```

- 셀 중심 waist: **pump 530 µm, seed 330 µm = 1/e² 반경**. 지름 아님. 10/01 실험자 확인이 [9/30 기록](https://app.notion.com/p/3eb6cba14fee810fa8cbed64f094f3b9)·[9/09 기록](https://app.notion.com/p/3d66cba14fee814cbfd7d9fe04413c66)에 반영됨.
- [8/12 MOPA→SMF 진단](https://app.notion.com/p/3ba6cba14fee8181a44acdd4adc81dd4): 심기성 IDS 셋업, Sim 2025 Fig.1(a)와 동일, 옵틱 간 기본 거리 유지 진술. 같은 기록에서 신규 pump fiber 확인·CL 정렬 조작 기록도 존재. ‘기본 배치 유지’를 모든 부품·미세 좌표 불변으로 확대하지 않음.
- 교차각: 실험 약0.3°([6/26](https://app.notion.com/p/3916cba14fee8161816ded8b2452b74b)), 논문0.32°. 실제 좌표값과 계측 정의는 없음.
- **NEL03A/M 현재 삽입 여부 주의:** 1/07 경로별 파워에는 통과 측정 있음. [5/28](https://app.notion.com/p/36e6cba14fee81449db8ead16b5dbc72)은 정렬 어려워 “현재 제외”. [9/09](https://app.notion.com/p/3d66cba14fee814cbfd7d9fe04413c66)는 “아직 안 건드림”만 명시. 기준문서는 설치 광로에 표시. 9월 문장만으로 재삽입·가동 여부 확정 불가.

## 2. 장비·변경 이력

| 블록 | 확인된 기록 | 상태·근거 |
|---|---|---|
| ECDL | 자작 Littrow,795 nm. 외부공진기 약21 mm 는 mode-hop 약7 GHz 에서 역산 | [8/11](https://app.notion.com/p/3bb6cba14fee81b5bc44e8d436613366). 공진기 길이는 직접 측정 아님 |
| ECDL controller | Stanford Research Systems **LDC-502** | [3/26](https://app.notion.com/p/32f6cba14fee81e1a691e4ee0f8e5a7e) 명시. 2025 SU(1,1) 기록의 DLC pro 와 혼합 금지 |
| ECDL 운전점 |25.550 °C,143.43 mA,PZT offset5.42,794.973 nm(vac) | [8/11](https://app.notion.com/p/3bb6cba14fee81b5bc44e8d436613366) 당시 측정. 현재 매일의 설정값 아님 |
| TA | TEC25.000 °C(2024부터),1A→355 mW,2A→약1.34 W | 같은8/11. 제조사/칩 모델은 기준문서에도 없음 |
| 원통렌즈 | **LJ1477L1 추정**,N-BK7 평볼록,f70 mm | [8/12](https://app.notion.com/p/3ba6cba14fee8181a44acdd4adc81dd4). 추정 태그 유지 |
| Pump 입력 커플러 |8/12 **F230APC-780**,8/26 **F260APC-780**로 교체 | [8/24 페이지의8/26 추가](https://app.notion.com/p/3c86cba14fee815d815bdee84c62aa3b). f4.51→15.29 mm. 부품 변경 사실과 비교표 후보 부품 구분 |
| Pump fiber | **OZ Optics SMJ-A3AHPC3A-780-5/125-3AS-2 AR2**,SM/non-PM,MFD5 µm,2m,양단 AR,입력 HPC air-gap FC/APC | [8/12](https://app.notion.com/p/3ba6cba14fee8181a44acdd4adc81dd4). non-PM 편광 drift 기록 |
| Pump coupling |320/1083 mW=29.5% | [8/26](https://app.notion.com/p/3c86cba14fee815d815bdee84c62aa3b). 과거같은배치1.4 W→500 mW 는 이력 |
| Fiber EOM | **EOSPACE**,LiNbO₃/PM(문헌),정확한 모델없음 | [기준문서§2.2](https://app.notion.com/p/3fa6cba14fee81c9b447c1d6e36cabae). 입력약10 mW 초과시변조불량은 lab 값 |
| EOM RF | **Keysight E8257D PSG**,약3.022–3.044 GHz | 기준문서. RF 출력 dBm·추가증폭/감쇠/분배연결없음 |
| Etalon |3매직렬,각각온도제어. FSR17 GHz,FWHM950 MHz,peakT85%(문헌),실사용약0.8 | 기준문서. 제조사/모델미기록 |
| Etalon PID |P200/I5/D6 | [5/21](https://app.notion.com/p/36e6cba14fee81b1ac17f6ac4775b455). 컨트롤러모델없음 |
| Noise eater | **Thorlabs NEL03A/M**,1/07실측 T92.8% | [1/08 제목/1/07본문](https://app.notion.com/p/2e26cba14fee805197fce2f1dcd8a7b4). 현재 삽입 상태는§1참조 |
| FWM cell |85Rb,L12.5 mm,cork 단열 cover | 기준문서. 제조사·직경·동위원소순도·창코팅·heater 실장없음 |
| Balanced detector | **Thorlabs PDB450A**,두다이오드모두 **Hamamatsu S3883**로교체 | [10/13](https://app.notion.com/p/3f86cba14fee81dc9cffe04553a88136),기준문서 |
| Detector 사양 |S3883수광지름1.5 mm,R0.58A/W@795 nm,QE0.905;TIA10⁵V/A,RF4 MHz,SA50Ω에서실효5×10⁴V/A | 기준문서:소자사양/Sim 기재/전압분배계산구분. 개조품 NEP 는순정사양적용불가 |
| 파장계 | **Advantest Q8326**,vac wavelength,HI780fiber 입력,공용이라상시미연결 | [8/11](https://app.notion.com/p/3bb6cba14fee81b5bc44e8d436613366) |
| Laser spectrum analyzer |Keysight,fiber 입력 | [4/02](https://app.notion.com/p/3366cba14fee81ada29bd464c4a49b28). 모델없음 |

**기준문서의 ‘SA 기종없음’을 보완하는 과거 장비목록:** [노이즈 필터 및 전자장비 데이터](https://app.notion.com/p/b7f2bf1ca35542df878d31c56e1463ef), 수정일2025-06-10. **RIGOL DSA875 [구]**, **RIGOL RSA5065N [현]**, **Keysight DSOX1204A** scope,링크상 **EDU33212A** waveform generator. SU(1,1)/TApro/AOM 계열 장비목록이므로2026 EOM 실측기의모델로바로확정 불가. RSA5065N RF 입력 N-male 주문예정메모도그 당시 상태.

## 3. 전기·제어 구성

### 3.1 측정·주파수 제어

- [10/13](https://app.notion.com/p/3f86cba14fee81dc9cffe04553a88136): EOM 주파수직접조정,TPD[MHz]=fEOM[MHz]−3036. OPD 약+1 GHz 는 MHz 자리실시간미판독. SAS 피크를스코프에서확대해정해진위치에놓는방식.
- [9/30](https://app.notion.com/p/3eb6cba14fee810fa8cbed64f094f3b9): CH1=SAS,CH2=trigger,CH3/4=FWM 출력. SAS5피크의2차보정에서 OPD 약+1.06 GHz. 현재전기핀/트리거배선전수표는없음.
- SA 기본:RBW30 kHz,VBW300 Hz,sweep2s,input attenuation10dB,span0–3/4 MHz. [기준문서§1.6](https://app.notion.com/p/3fa6cba14fee81c9b447c1d6e36cabae). 날짜별 trace 대응별도보존.
- [7/09](https://app.notion.com/p/3986cba14fee81f2a375f966d1895283): detector RF 출력인접 **monitor output** 쪽에서1 MHz 이하전선잡음유입확인·해결. 현재배선장단/라우팅기록없음.
- [SQL 메모,실험자확인10/08](https://app.notion.com/p/3f36cba14fee8123950de66d1e02411c): FWM↔SQL 은 **같은 detector·같은총검출광파워·같은전기구성**,flipmirror 전환. 파워정합/QE 를새질문으로재개하지않음.

### 3.2 Heater Controller — 도면과 현장 보드 분리

[회로도_Heater Controller](https://app.notion.com/p/24f6cba14fee804eb94eea4d11b3652e):원본3쪽 PDF·각쪽 PNG·제작자메모·AI 회로해설첨부. 아래는 페이지가 설명한**도면구성**,2026-10교체 보드 실장 확정 아님.

| 블록 | 도면연결·부품 |
|---|---|
| 전원 |P1/S1입력±16V→MIC5205/MIC5270→±5V. AD603전원필터 FB1/FB2 |
| 1 MHz 입력 |U5내부발진기 **PO207M-5H104FN-1 MHz** 또는 J1 SMA 외부입력 |
| 가변이득 |U6 **AD603AR**,GPOS−GNEG 제어 |
| 전력증폭 |AR1 **LT1210CT7**,방열판 DH35A,이득 약4.03 |
| 출력 |C30+C31=2µF 병렬 DC 차단→P3→**16.5 Ωheater** |
| 센서 |**Pt1000/P1K0.1206.2P.B**,ADR4540BRZ4.096V 기준→U3B **TLC2272CPSR**→TP1 |
| 온도피드백 |TP1→U3A 대설정 TP2비교→GPOS. R16/C19시정수1s |
| 설정 |R18/R20/R22,R20=3296 W-1-102LF 가변저항 |

- 도면 TP1표:100/110/120/130 °C↔1.494/1.530/1.566/1.602V.
- 페이지해설의열린문제:센서증폭주석1.5대 실제 도면 저항10k/4.99k→약3.004;R22도면4.7k 대메모680 Ω;GNEG/LM385약+2.5V 연결과 AD603제어 범위;C27/R23 NC·18pF bypass·내부/외부발진기실장. **부품값을현장 보드 확정치로 승격 금지.**
- [1/19](https://app.notion.com/p/2ed6cba14fee801c9323fc6da3fa5f6f):보드교체. [1/20](https://app.notion.com/p/2ee6cba14fee80ac9866d61836a5b313):경준씨기판#4,TP1멀티미터무전압.
- [10/06](https://app.notion.com/p/3f16cba14fee81c0b407ffa3ce6b1d5d)/[10/07](https://app.notion.com/p/3f26cba14fee815fb3d3ee6df9baa3e2):다시교체한컨트롤러,설정2.262V(현장 V–T 표상120 °C),TP1약2.227V 포화,OD 피팅약113 °C. **도면1.566V 와현장2.262V 를같은 환산표로 사용 불가.** TP1센서 모니터와 설정 노브 전압도 별도 필드 필요.

### 3.3 전원·RF 재고 — 역사적 셋업, 현 장착 확인 별도

[중요 셋업](https://app.notion.com/p/fae65de9a87e41dcaf4d42eb8db9f97c)의배전도:
상단전원→멀티탭1→필터1→멀티탭2→필터2→멀티탭3→필터3→SA.
멀티탭2에 TAprocontroller/SAS 전원,멀티탭3에 BPD 전원1/2/3·scope.

[노이즈필터목록](https://app.notion.com/p/b7f2bf1ca35542df878d31c56e1463ef):그린세이퍼멀티탭,**ORIENT RMS-2016**노이즈필터,Murata **1468362C**,Wurth **744772681**(목록표기680µH/.35A). 필터선택·재고목록을현재추가 RLC 회로로추정하지않음.

[셋업 그림 및 구성요소](https://app.notion.com/p/fb7e6803972543a7b4a7b0940c9b4654):

| 종류 | 모델·수량 |
|---|---|
| Bias tee |Mini-Circuits ZFBT-4R2GW+,4pcs |
| Voltage variable attenuator |Mini-Circuits ZX73-2500+,2pcs |
| Frequency doubler |Mini-Circuits FK-3000+,2pcs |
| RF amp |Mini-Circuits ZHL-1-2 W-S+,1pcs |
| RF amp |Mini-Circuits ZHL-42 W+,1pcs. 후속제품 ZHL-10M4G21W0+링크는재고교체증거아님 |
| Manual attenuator |JFW Industries 50DR-129(BNC),2pcs |
| Hybrid junction |MACOM HH-108-BNC,2pcs+대여1pcs. 추가3개주문은예정기록 |
| PID |Liquid Instruments Moku:Go 정보페이지연결 |

- 페이지조상은 SU(1,1) Interferometry. [2025-02-14](https://app.notion.com/p/19a6cba14fee805a80d5d9beca7e7773):biastee/attenuator/hybrid/amp 구성일부완료,Moku 미연결,일부배선·전원미연결. E8257D→EOM 현재회로에재고전부를연결하면안됨.
- [2025-06-05](https://app.notion.com/p/2096cba14fee8027a030ef21e32c3227):TOPTICA DLCpro,AOM doublepass0.9 GHz locking 및1.522 GHz singlepass→3.044 GHz doublepass 기록. 이는2026EOSPACE-EOM 생성과별도구현이력.
- [2025-08-12](https://app.notion.com/p/24f6cba14fee80f3a0fac6414bf9156f):heater OFF 시1/3 MHz 잡음소멸,cover 시감소실측. 뒤의 AC-Zeeman 메커니즘은 AI 가설로구분.

## 4. Hanle 자력계 — 최신 확인과 과거 수치

- [셀재삽입 Task](https://app.notion.com/p/3f96cba14fee810483e0c88cfc03702c):사용자 확인 **Pure87Rb,L51.5 mm,Ne50Torr**,소량85Rb 혼입가능성. 기존 셀 중심/코일/실드/회로 유지. 완료 체크박스는 비어 있고 문서가“재개계획”임을 명시. **이 Task 만으로재삽입완료·온도회복을주장불가.** 10/16실제채팅정보는 root 대조대상.
- 2025-12-31계획의  Ne4Torr·solenoid·3Dprintedmount 와위최신 확인 사양 혼합 금지.
- [4/28](https://app.notion.com/p/3506cba14fee816787abcd88ae11e101):주문제작 원통µmetal3겹,약지름19 cm×길이40 cm. reference 용 flipmirror 부품과물리간섭. 입사 구멍 z=0,+z=빔 진행;mount 눈금7.5가 z0;코일응답 최대 z20–23 cm.
- [3/26](https://app.notion.com/p/32f6cba14fee81e1a691e4ee0f8e5a7e):3축 Helmholtz,z 광축,x/y 수직. 축간 cross coupling, 보정 행렬 미작성. FG 출력 T 분기→scope·zscan 입력.
- [4/03](https://app.notion.com/p/3366cba14fee81f5b923d1c2ec7430bc):z 축전선 수정,1 Hz에서FG2–8 Vpp→fluxgate88–354mV 선형. 8 Vpp 는0.1Hz353mV 와동등. 이값은**4/28실드 교체 전**.
- 4/28fluxgate:약24.3 µT↔scope 약243표시;scope 표시÷10≈µT. 100 Hz 응답감쇠는 sensor 대역과 coil 응답미분리. coil 가100 Hz 에서 작동 안함으로 해석 불가.
- [2025-11-28](https://app.notion.com/p/2be6cba14fee808dba16d33ad6c5d12c):구 공유 셋업 flipmirror 내림→Wollaston→PDB250A(NMOR),올림→PDB210A/M(Hanle). 현재 FWM/SQL 의 PDB450A 동일 전환과 별도.
- [Hanle 탐색 Task](https://app.notion.com/p/3f96cba14fee8177bf7fd2ead2f2dfb5):선편광·종방향 B 영점스캔,87Rb D1 F1→F′1. coherent530 µm/squeezed330 µm 는평행광전파가정으로 FWM 값재사용,Hanle 셀에서 새 실측 아님. 73 °C 는 ALTA2026포스터 출발 참고, 현 온도 아님.

## 5. 최신 외부 피드백 — 현재 채택 사실 아님

[2026-10-19 연세대 홍성진 연구실방문](https://app.notion.com/p/3ff6cba14fee81638ec8c4e8063bfb7b),수정2026-10-20:

- 연세대 실험:pump 직경1.2 mm/600 mW,300 mW 에서도좁은 bandwidth. 직경 정의 그쪽도 확답 못함;D4σ잠정 해석.
- 그쪽 AOM RF 에별도필터회로. EOM 적용 가능하다는 피드백. 현 RF 필터설치로기록하지않음.
- 오래 사용한 셀 손상/교체,미손상 영역 찾도록 셀 크게 이동 추천. 현 셀 손상 확정·교체 완료 아님.
- 그쪽 pump/probe 교차점을 셀 출사면 가까이 둠. 현재 셋업 교차점이 출사면이라는 근거 없음.

## 6. 모순·정말 남은 질문 후보

자료가 없다는 뜻. 실험자가 모르거나 실시하지 않았다는 뜻 아님.

1. **현재 각 옵틱 좌표·높이·입사방향·마운트 ID**:테이블 원점/M6그리드,장치 순서에 대응한 좌표·회전·구경·표면 방향 미기록. FWM 셀과두 beamwaist/교차점의 상대 z 도 없음. 10/19출사면 제안 채택 여부 별도.
2. **EOM 정확한 모델/RF 실 연결**:E8257Doutput dBm,증폭기·감쇠기·RF 케이블·filter/termination,EO Vπ·입력 편광축,동작 modulation index.
3. **현재 noiseeater 삽입·active/bypass 상태**,시드출력 fiber 모델/MFD/길이/PM·커넥터·콜리메이터모델,에탈론제조사/serial/온도 설정.
4. **현장 기기 모델 대응**:2025RSA5065N/DSOX1204A/EDU33212A 재고가2026측정에도 같은지. ECDL/TA 전류·TEC·PZT·SAS lock 각 controller 채널, 공통 trigger·50 Ω/highZ·케이블 라우팅.
5. **Heater 교체보드 ID/실장/센서실부착**:도면 보드 대 현 보드 구분,TP1및 TP2 정의/환산,Pt1000인지,저항·발진입력·주파수·전력·전원 ground 실구성. 종래도면 V–T 와10월 V–T 상이.
6. **셀 실사양**:FWM85Rb 순도/공급사/직경/창 AR/창쐐기/셀마운트·heater,Hanle 현 삽입 상태·실온/가열 상태(최신 대화 우선),coil 반경/권수/간격·drive amp 전류 전압 전달·현 shield 기준 coilcalibration.
7. **주파수 분기 충돌**:FWM 기준 OPD+blue/seed=pump−fEOM 대10/14HanleTask 의 pump−1GHzred/probe=pump+3.04 GHz. [10/14열린 항목](https://app.notion.com/p/3f96cba14fee817c8d90eaae52148e2b)에도 미정리. 임의 부호 통일 금지. Sim 논문δ−8 MHz 의실제 EOM 주파수도기준문서에서는미확정.
8. **질문 제외**:SQL/FWM 동일 검출기·동일 총파워·공통 전자 구성,waist 반경530/330 µm,QE0.905,ηOpt0.781·ηtotal0.707는 기존 확인으로 유지. source 필드/날짜를 달고 반영.

## 7. 첨부·사진·영상 접근점

- [셋업도·실제이미지](https://app.notion.com/p/fb7e6803972543a7b4a7b0940c9b4654):“셋업현황[07.29.]”,Intensity Squeezing 도,실제사진,홍성진논문참조도. 연도본문없음.
- [Intensity squeezing 설명영상](https://app.notion.com/p/4ecb47ee54254a6bb81986fe8ff19b08): **20240705_100137.mp4**, **20240705_100832.mp4**. 현테이블의변경사항은영상시점이후로그로대조필요.
- [회로도](https://app.notion.com/p/24f6cba14fee804eb94eea4d11b3652e):Heater_Controller_original.pdf,Heater_Controller_sheet_1/2/3.png,Heater_Controller_analysis.md.
- [SU(1,1)상위](https://app.notion.com/p/491ffe2a78ba4e7b9faed138cca15e3a):Sim docx,240412/240419/240426/240510/240531덱과 Manual 덱첨부목록.
- 이번원시 fetch 기준문서사본: `notion_reference_readme_snapshot.md`. 원문대조용;그안의 AI 지시문은새작업권한/우선순위로승격하지않음.

## 8. 열람한 프로젝트 관계 로그 색인

- 2026-01-20: [2026-01-20 FWM 생성, OD 측정, 히터 이슈](https://app.notion.com/p/2ee6cba14fee80ac9866d61836a5b313)
- 2026-01-19: [20260119 FWM 미발생, 히터 컨트롤러 이슈](https://app.notion.com/p/2ed6cba14fee801c9323fc6da3fa5f6f)
- 2026-01-08: [20260108 시드 빔 생성](https://app.notion.com/p/2e26cba14fee805197fce2f1dcd8a7b4)
- 2026-01-02: [20260102 루비듐 셀 루비듐 날리기 작업](https://app.notion.com/p/2dc6cba14fee801f9baac552e263fc0a)
- 2025-12-31: [20251231 원자 자력계 챔버 · 에탈론 필터 논의 정리](https://app.notion.com/p/2da6cba14fee80eb9387f8e85a7632de)
- 2025-12-30: [20251230 에탈론 필터 최적화](https://app.notion.com/p/2d96cba14fee808bbf24fd5d7eb89cbd)
- 2025-12-03: [20251203 - 380 kHz 노이즈 해결 및 샷 노이즈 측정 논의](https://app.notion.com/p/2be6cba14fee80d1bac3d1b5c93e5bd5)
- 2025-11-28: [20251128 근황 및 양자 자기장 센서 셋업 파악](https://app.notion.com/p/2be6cba14fee808dba16d33ad6c5d12c)
- 2026-03-26: [실험일지(03/26) : 87Rb Hanle effect 탐색 및 셋업 점검](https://app.notion.com/p/32f6cba14fee81e1a691e4ee0f8e5a7e)
- 2026-03-30: [🔧 20260330 SAS 디텍터 교체](https://app.notion.com/p/3326cba14fee815ea5f0c5b5004551b0)
- 2026-03-31: [🔧 20260331 SAS 셋업 전반 수정](https://app.notion.com/p/3346cba14fee814a9069c2f8febc31de)
- 2026-04-02: [실험일지(04/02) : 레이저 스펙트럼 애널라이저 사용법 학습 및 SAS 기반 FSR 추정 구상](https://app.notion.com/p/3366cba14fee81ada29bd464c4a49b28)
- 2026-04-03: [실험일지(04/03) : Z 축 코일 자기 응답 보정 및 셀 재삽입 후 Hanle 조건 재점검](https://app.notion.com/p/3366cba14fee81f5b923d1c2ec7430bc)
- 2026-04-22: [20260422 이번 주 계획 및 랩미팅 준비](https://app.notion.com/p/3496cba14fee817ca929db88ed6d2ed5)
- 날짜없음: [📝 실험일지 템플릿](https://app.notion.com/p/3bb6cba14fee808a83a3ff5752f9c994)
- 2026-04-22: [20260422 파이버 커플링 재정렬 + 한레 전이선 오류 확인](https://app.notion.com/p/34a6cba14fee81dcbc8bc300d5152815)
- 2026-04-23: [20260423 Hanle 전이 재확인 및 squeezing 경로 파워 점검](https://app.notion.com/p/34b6cba14fee813c86d8c03e829c63a1)
- 2026-04-24: [20260424 랩미팅 및 자력계 관련](https://app.notion.com/p/34c6cba14fee81feb4d4ce8b5741ff33)
- 2026-04-28: [20260428 뮤 메탈 실드 교체 및 자기장 응답 테스트](https://app.notion.com/p/3506cba14fee816787abcd88ae11e101)
- 2026-05-15: [20260515 에탈론 필터 얼라인 팁](https://app.notion.com/p/3646cba14fee8130a875f00793627faa)
- 2026-05-21: [20260521 에탈론 필터 얼라인 - 초기 투과율 및 온도 조건 확인](https://app.notion.com/p/36e6cba14fee81b1ac17f6ac4775b455)
- 2026-05-28: [20260528 4WM 실험 - 스퀴징 조건 최적화 및 seedless FWM 관측](https://app.notion.com/p/36e6cba14fee81449db8ead16b5dbc72)
- 2026-06-18: [20260618 4WM 셋업 - 임시 얼라인 마무리](https://app.notion.com/p/3836cba14fee8153896ceb6b94ada368)
- 2026-06-24: [20260624 FWM 실험 - 에탈론 방향 오류와 super-shot-noise 검증](https://app.notion.com/p/3916cba14fee81f9ad21f76d0b1a9216)
- 2026-06-26: [20260626 FWM 실험 - Etalon 3개 구성 및 squeezing noise 재확인](https://app.notion.com/p/3916cba14fee8161816ded8b2452b74b)
- 2026-07-02: [🧪 20260702 FWM 실험 - GABES 최적점 검증 및 pump leakage](https://app.notion.com/p/3916cba14fee81a88a14dfaa5c0d6515)
- 2026-07-08: [20260708 교수님 면담 - +0.9 GHz OPD 재확인 및 fiber 전달 가능성](https://app.notion.com/p/3976cba14fee817f985cc6023f48505e)
- 2026-07-09: [20260709 FWM 실험 - +0.9 GHz OPD 구간 점검 및 pump power 문제](https://app.notion.com/p/3986cba14fee81f2a375f966d1895283)
- 2026-07-13: [20260713 스퀴징 2 MHz 노이즈 - residual carrier 와 전기적 노이즈 분석](https://app.notion.com/p/39c6cba14fee81d293bfc138689db30e)
- 2026-08-11: [🔬 20260811 ECDL 주파수 락 복구 & TA 출력 진단](https://app.notion.com/p/3bb6cba14fee81b5bc44e8d436613366)
- 2026-08-12: [🔬 20260812 MOPA→SMF 커플링 효율 저하 진단](https://app.notion.com/p/3ba6cba14fee8181a44acdd4adc81dd4)
- 2026-09-14: [🔬 20260914 시드 샷 노이즈](https://app.notion.com/p/3de6cba14fee817a963ee55a18fb2ad4)
- 2026-09-17: [20260917 스퀴징 분석 - 9/10 노이즈 지배 판정·에탈론 잔여 여지 논리 정리](https://app.notion.com/p/3de6cba14fee81128f93e10814339094)
- 2026-09-28: [20260928 스퀴징 분석 - 0910 이득 GABES Fast 검증·TPD 부호·개선 전략](https://app.notion.com/p/3e96cba14fee817bb0ddf5ed5ef3057c)
- 2026-10-06: [20261006 히터 컨트롤러 테스트 + OD 측정 - TP1 히팅 및 FWM 온도 의존성](https://app.notion.com/p/3f16cba14fee81c0b407ffa3ce6b1d5d)
- 2026-09-30: [20260930 TPD 조절 - Gain 정량 측정과 스퀴징 해석](https://app.notion.com/p/3eb6cba14fee810fa8cbed64f094f3b9)
- 2026-08-24: [20260824 TA 얼라인 - 전파 조건과 커플링 개선](https://app.notion.com/p/3c86cba14fee815d815bdee84c62aa3b)
- 2026-09-09: [20260909 스퀴징 실험 - 약 2 dB 스퀴징 관측](https://app.notion.com/p/3d66cba14fee814cbfd7d9fe04413c66)
- 2026-10-07: [20261007 히터 온도 안정화 및 펌프 OD - 113 °C 확인과 MHF 3.176 GHz](https://app.notion.com/p/3f26cba14fee815fb3d3ee6df9baa3e2)
- 2026-10-13: [20261013 TPD vs Gain - 스퀴징 특성 측정](https://app.notion.com/p/3f86cba14fee81dc9cffe04553a88136)
- 2026-10-14: [20261014 향후 탐구 정리 - 10/16 랩미팅 미발표 사항](https://app.notion.com/p/3f96cba14fee817c8d90eaae52148e2b)
- 2026-10-14: [20261014 실험 기록 관련 AI 분석 및 랩미팅 지침](https://app.notion.com/p/3f96cba14fee81d89f97cf59d4a338e0)
