# SABES 현행 구현·Sim 논문 대조

조사일: 2026-10-21. 작성: AI. 채택 상태: 조사 자료·구현 변경 제안. 코드·기본값 변경 없음.

범위: 현행 `sabes/`, 로컬 Sim 2025 논문 전 9쪽, 셋업 원본 PNG 2개, PDB450A 매뉴얼, S3883 데이터시트. Allen 2026 논문은 다른 실험실 장비가 섞였는지 대조하는 용도로만 확인. 현재 랩미팅의 EOM 기종은 2026-09-04 슬라이드 4·5에서 좁게 교차 확인. 전체 랩미팅·Notion 조사 결과는 상위 조사 문서 참조.

기계 판독 원장: [sabes_current_inventory.json](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/docs/SABES/apparatus_research_2026_10_21/sabes_current_inventory.json>). 63노드·61광로·85계수의 원시 필드,정의 파일·줄,기존 provenance,추가 근거상태,설정값,기본 기하 계산 결과 포함. 기존 `lab`/`fitted`/`measured` 꼬리표를 보존하면서 실제 측정과 구분. 원본 코드 SHA-256 함께 저장.

## 출처

| ID | 출처·정확한 위치 | 근거 성격 |
|---|---|---|
| P25 | [Sim, Kim, Moon, Scientific Reports 15, 7727 (2025)](https://doi.org/10.1038/s41598-025-86479-w), [로컬 PDF](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/references/Intensity-difference squeezing from four-wave mixing in hot 85Rb and 87Rb atoms in electro-optic-modulation-based single diode laser pumping system.pdf>) | 해당 논문 장치·운전점. 오늘 장치의 자동 증거 아님 |
| D1 | [IDS setup part 1.png](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/references/IDS setup part 1.png>) | 기존 저장 셋업 그림. 축척·작성일·좌표 원점 없음 |
| D2 | [IDS setup part 2.png](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/references/IDS setup part 2.png>) | 기존 저장 셋업 그림. SQL 우회광로 포함 |
| C | [calibration_sim2025.json](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/sabes/calibration_sim2025.json>) | 현행 설정·추정·맞춤 계수. 실측 여부 개별 판별 필요 |
| B | [beamline.py:54](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/sabes/beamline.py:54>) | 설정 및 광원 전달 모델 |
| L | [layout/parts.py:56](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/sabes/layout/parts.py:56>) | 노드·광로·길이·손실 경로 |
| T | [layout/table.py:6](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/sabes/layout/table.py:6>) | 표시 좌표·물리 길이 분리 계약 |
| DET | [detection.py:56](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/sabes/detection.py:56>) | 후단 검출 모델 |
| E | [instruments/electrical.py:130](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/sabes/instruments/electrical.py:130>) | 합성 SA·오실로스코프 출력 |
| M | [PDB450A-Manual.pdf](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/references/PDB450A-Manual.pdf>) | Thorlabs 2019 매뉴얼. PDF쪽/인쇄쪽 구분 |
| H | [s3883 광다이오드.pdf](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/references/s3883 광다이오드.pdf>), [Hamamatsu 제품 페이지](https://www.hamamatsu.com/us/en/product/optical-sensors/photodiodes/si-photodiodes/S3883.html) | 소자 명목 사양. 개조 BPD 전체 성능 아님 |
| A26 | [Allen 외, OE 34, 19382 (2026)](https://doi.org/10.1364/OE.600911), [로컬 PDF](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/references/2026 OE Long-term stabilization ids.pdf>) | Oak Ridge/JILA 실험. 사용자 장치로 이전 금지 |
| LM0904 | [박수인 랩미팅 20260904.pptx](<G:/내 드라이브/원자광학연구실/Lab meeting/박수인 랩미팅 20260904.pptx>) 슬라이드 4·5, `ppt/slides/slide4.xml`, `slide5.xml` | 현재 계보의 실험기록. 본문 날짜 09/03 |

PDF쪽은 파일 첫 페이지를 1로 계산. PNG는 투명 배경 때문에 검은 글씨가 안 보일 수 있음. [Part 1 백색 배경 렌더](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/docs/SABES/apparatus_research_2026_10_21/source_preview/IDS_setup_part_1_white.png>)·[Part 2 백색 배경 렌더](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/docs/SABES/apparatus_research_2026_10_21/source_preview/IDS_setup_part_2_white.png>)는 원본을 흰 PDF 페이지에 배치한 열람용 렌더. 원본 미변경.

## Sim 2025에서 직접 확인되는 장치·운전값

| 항목 | 논문에 실제 기재된 정보 | 위치·조건 |
|---|---|---|
| 광원 | 자체 제작 ECDL, 795 nm, Rb D1 부근 | P25 p.2 Experimental schematic |
| 증폭 | 자체 제작 MOPA·tapered amplifier; SMF 공간 모드 정제 언급 | P25 p.2. TA 칩·드라이버 모델, 전류·온도·출력 곡선 없음 |
| 분기 | 증폭광 대부분 pump, PBS로 일부를 fiber EOM에 분기 | P25 p.2·Fig.1(a) |
| EOM 기종 | **EOSPACE AZ-0S5-10-PFA-PFA-780-UL** | P25 p.2. 문자 그대로 기록. 현행 랩미팅 PM-0S5 표기와 다름 |
| EOM 역할 | pump보다 red인 probe-seed 생성; 원자종별 약 3.036/6.835 GHz 간격; carrier·다른 sideband 동반 | P25 p.2 |
| 필터 | etalon 3개 직렬 | P25 p.2·Fig.1(a) |
| etalon 명목 사양 | FWHM 950 MHz, FSR 17 GHz, peak transmission 85% | P25 p.2. 제조사·모델·온도·각도·개별 실측 스캔 없음 |
| 셀 앞 빔 형상 | pump·seed 각각 두 렌즈 4F 계통, 셀에서 느슨한 초점 | P25 pp.2–3, Fig.1(a): L1–L2·L3–L4. 초점거리 수치 없음 |
| 셀 종류·길이 | isotopically pure 85Rb 및 87Rb 셀, 각각 12.5 mm | P25 p.3. 제조사·모델·순도 수치·지름·창 코팅·buffer gas 언급 없음 |
| 셀 온도 | 85Rb 약121 °C; 87Rb 약131 °C | P25 p.3·p.4. 온도 센서 위치·heater/PID 계통 없음 |
| 셀 중심 pump 모드 | **530 µm, 1/e² intensity radius**, 600 mW | P25 p.3. 지름 아님 |
| 셀 중심 seed 모드 | **330 µm, 1/e² intensity radius**; setup 설명10 µW | P25 p.3. 최고 IDS 측정에서는8 µW(p.4) |
| pump–seed 각 | 약0.32°, 최고 squeezing에서 선택 | P25 p.3. 3차원 기울기·교차 위치·거리 측정 없음 |
| conjugate 방향 | pump 반대편 약0.32° | P25 p.3·Fig.1(b). probe–conjugate 총각 약0.64°는 기하 유도 |
| 편광·결합 | PBS 결합·분리; pump와 twins 분리 | P25 Fig.1·p.3. 실물 H/V 축·마운트 각도 없음 |
| pump 제거 | PBS 약3000:1 extinction, 다음 D-shaped mirror 공간 분리 | P25 p.3. 산란이 주요 한계. 셀→DM 거리 없음 |
| BPD | transimpedance **10⁵ V/A** | P25 p.3. 제조사·PDB450A 명칭·S3883 명칭 없음 |
| SA | RBW **30 kHz**, VBW **300 Hz** | P25 p.3. 제조사·모델·종단·검파·평균·cable 정보 없음 |
| SQL 광로 | 셀 앞 pump를 우회시켜 50:50 BS로 coherent beams 두 개 생성, intensity difference 측정 | P25 p.3. 실제 flip 동작·출력 맞춤 조건은 사용자 기록 별도 |
| 85Rb 최고 IDS 운전점 | T121 °C, Δ≈+0.9 GHz, δ=−8 MHz, pump600 mW, seed8 µW | P25 pp.3–4. Δ 기준85Rb F=2→F′=3 Fig.2 caption |
| 85Rb 출력 | probe111 µW, conjugate109 µW, 합220 µW; gain≈15 | P25 p.4. 논문 반올림 수치:111/8=13.875와 gain≈15를 같은 정확도로 강제 일치시키지 않음 |
| 85Rb 결과 | 직접 IDS−7.8 dB, 주영역40–500 kHz; 대역 약3.5 MHz | P25 p.4 Fig.3 |
| 잡음 특징 | 약2 MHz 봉우리: laser-source 전자 잡음으로 서술 | P25 p.4. 특정 전원장치·회로까지 특정하지 않음 |
| power slope | 분석주파수150 kHz; FWM/SQL 기울기비0.134≈−8.7 dB | P25 pp.4–5 Fig.3 |
| loss 서술 | optical5.5(2)%, system8.0(2)%, total detection13.5(2)% | P25 p.5. 논문은8%를 QE 손실이라고 직접 명명하지 않음 |
| 87Rb 운전점·결과 | Δ≈1.3 GHz, δ=−11 MHz, T131 °C; pump600 mW; IDS−5.0 dB | P25 pp.4–5. 85Rb trace와 혼합 금지 |
| carrier 억제 | 2 etalons: carrier/seed≈0.5%; 3 etalons: 훨씬0.5% 미만 | P25 pp.5–6. 3개 정확한 비율 없음 |
| 데이터 가용성 | 대응 저자에게 합리적 요청 시 제공 | P25 p.7. 원자료·부록 공개 제공 없음 |

P25 p.6의 “0.8 GHz squeezing bandwidth reduction”은 같은 논문의 MHz 대역과 일관되지 않는 단위 표기. ⟦AI⟧ 추정 단위로 고쳐 수치 입력하지 말 것. 모델 근거로는 정성적 감소만 사용.

## 기존 셋업 그림에서 확인되는 실제 연결 순서

### Part 1

그림 D1 관찰. 위치는 화면상의 상대 배치만 의미. 빔 진행방향 따라 기록.

- ECDL → OI → HWP → 편광 분기 → MOPA/TA 입력. 별도 가지: f=300 렌즈 → 거울 → natural-Rb reference cell → PBS → PD.
- TA 출력 → CL → steering mirrors → OI → HWP → PBS 분기. **CL은 그림상 TA 출력 쪽**. SABES note는 “amplifier input”이라고 적어 출입 방향 대조 필요.
- pump 분기 → F230 → delivery fiber.
- seed 분기 → Polarizer → steering mirror → HWP → F220 → fiber EOM → 출력 collimator → etalon1 → mirror → HWP → mirror → etalon2 → etalon3 → GT → **Noise Eater** → mirror → F220 delivery fiber.
- 따라서 D1의 Noise Eater는 **3-etalon/GT 뒤·마지막 fiber 앞**. 현행 SABES `ROUTE_SEED_TO_EOM`의 upstream 배치는 원본과 다름.

### Part 2

그림 D2 관찰. “5724”, “5725”는 그림 표기만 확보; 제조사·전체 품번 미확정.

- seed fiber 출력 collimator **5724 (1.6 mm)** → QWP → HWP → f400 → f150 → steering mirror → GTP → steering mirror → 입구 PBS → Rb cell.
- pump fiber 출력 collimator **5725 (2.0 mm)** → QWP → HWP → f300 → f150 → steering mirror → GTP → 입구 PBS → Rb cell.
- 셀 후 PBS가 pump를 beam blocker 방향으로 분리. probe·conjugate는 벌어져 각 mirror → iris → flip mirror 통과/반사 경로 → 마지막 반사계 → 검출 렌즈 → BPD.
- probe 최종 검출 렌즈 f75, conjugate f100. D2에 별도 f75/f100 표기 확인.
- 입구 PBS에서 빠지는 pump 우회광 → mirror → **HWP → BS → 두 flip mirrors** → 같은 BPD 두 포트. SQL 경로가 그림에 명시됨.
- 그림의 pump–seed 각 표기는 약0.3°. P25의 정량0.32°와 구분 보관.

## 현행 SABES가 저장·계산하는 값

조사 시점 `default_calibration()`·`nominal_settings()`·`build_source_chain()`·검출 geometry를 직접 실행. 원자 solver 실행 없음.

| 분류 | 값·상태 | 정확한 근거 |
|---|---|---|
| calibration 개수 | **85개**: nominal47, lab10, fitted13, paper8, datasheet7 | C 전체. 기존 UI plan의42개·placeholder25개 서술은 과거 스냅샷 |
| provenance 구현 | `fitted`, `datasheet`, `paper`를 모두 `verified=True`로 취급 | [calibration.py:32](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/sabes/calibration.py:32>) |
| display 노드 | Part1 38개·Part2 25개 =63개 | L. 노드의 x/y는 임의 도식좌표 |
| link | Part1 36개·Part2 25개 =61개 | L. Part1 36개 모두 `path_m=None`; Part2 14개에 길이 존재 |
| 광원 UI 기본값 | ECDL40 mW; TA2.5 A; HWP_TA0° | B:54. 실험기록으로부터 검증된 운전점 아님 |
| ECDL 조절 기록 | 온도·전류·PZT 기본 NaN, TA 온도 NaN | B:89. laser 주파수·lock 동역학 계산 안 함 |
| pump 분기 HWP | nominal solver8.099631° | [beamline.py:533](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/sabes/beamline.py:533>), 목표600 mW 역산. 마운트 실측각 아님 |
| EOM 앞 polarizer | nominal solver19.160651°, 목표 입력7.0 mW | 같은 함수. 실험각 아님 |
| seed trim HWP | nominal solver37.633311°, 목표8 µW | 같은 함수. 실험각 아님 |
| EOM RF | 18 dBm; `f_EOM=νHF+8 MHz=3043.732439 MHz` | B:68–70. δ=νHF−f_EOM=−8 MHz 모델 규약 |
| EOM Vπ | 4.4 V | C:80 `nominal`. exact PM 또는 AZ 기종의 실측 사양 아님 |
| EOM optical loss | 5.3 dB | C:62 `lab`지만 주석상−1 sideband 효율 약10%를 맞춘 역산. 독립 측정 아님 |
| EOM optical limit |10 mW | C `lab`. 논문에 이 rating 없음. 현재 기록/제품 정확 사양과 구분 |
| etalon peak | **80%** | C:50 `lab`. P25의 명목85%와 구분 |
| etalon leakage | 단계당0.076497 | C:44 `fitted`. 두-etalon carrier/seed≈0.5% 한 숫자에 맞춘 phenomenological floor |
| nominal 3-etalon 결과 | carrier/wanted seed≈0.000595410=**0.0595410%** | 위 코드 실행. 독립 측정 스펙트럼 아님 |
| collimator·f | pump 지름2.0 mm, f300→150; seed 지름1.6 mm, f400→150 | C `lab`, D2와 수치 일치. 지름 정의는 그림 자체에 없음 |
| waist 보정 | 모든 모드 반지름×1.08 | C:212 `fitted`. 소자 아님·실제 빔측정 대체 아님 |
| 실제 nominal 셀 waist | pump**540 µm**, seed**324 µm** | 1000×0.5×1.08;800×0.375×1.08. 논문530/330 µm와 정확 일치하지 않음 |
| pump fiber 효율 |70% | C:128 `nominal`. TA 출력의 M² penalty와 다른 광손실도 별도 곱함 |
| seed fiber 효율 | 입력60% nominal, 마지막50% lab | C. 서로 다른 fiber 구간 기준면 확인 필요 |
| TA 계수 | threshold0.5 A, slope1.05 W/A, saturation20 mW, M²1.6, ASE3% | C 모두 nominal. 합리적 형태를 둔 모형 가정 |
| 셀·geometry 기본 |121 °C, Δ0.9 GHz, cell12.5 mm | B·C. 논문 기준값 |
| separation knob |3.2628 mm, 셀→DM0.5842 m에서 arctan해0.3199978° | B:107, L:31. pump↔한 twin 기준. twins 전체 간격6.5256 mm는 모델 유도 |
| 검출 렌즈 | probe75 mm, conjugate100 mm | DET:70, D2 |
| PD 위치 | 렌즈 뒤 Gaussian waist로 자동 배치, defocus0 mm | DET:199. 실측 렌즈–PD 거리 입력 아님 |
| nominal spot | probe 반지름26.315 µm, conjugate35.954 µm | 검출 모델 실행값. 실측 프로파일 아님 |
| cell 후 광학효율 |94.500033% | nominal per-optic 손실과 residual1.004074를 맞춰 P25의5.5% 손실 재현 |
| 검출 QE |90.453881% | DET:236, `R=0.58 A/W`, λ795 nm로 계산. 고정 responsivity 가정 |
| nominal 합효율 |약85.47895%=94.500033%×90.453881% | ⟦AI⟧ 위 두 모델 입력의 곱. 현재 실험 효율 아님 |
| pump leakage |−60 dBm | DET:69. 설명은 observed 입력이지만 **기본 숫자 자체는 실측 근거 없음** |
| unwanted-mode RIN |−100 dBc/Hz, 분석대역100 kHz | DET:75–78. 주석이 명시적으로 simulation setting임을 밝힘 |
| SA 기본 |0.1–5 MHz,401 points,RBW30 kHz | E:138. VBW·averaging·detector mode 없음 |
| SA transfer |단일 pole4 MHz roll-off,50 Ω로 noise power 변환 | E:65·186. 실제 PDB 개조 전달함수·SA 교정 아님 |

### 현재 물리 거리의 출처 등급

| 구간 | 현재 값 | 코드 등급·의미 |
|---|---|---|
| 셀→D mirrors |584.2 mm=23 inch | L:31 `measured` 선언. C/P25/D1/D2에는 측정 원문·날짜 없음. 현재 장치에서 재확인할 역사값 |
| seed: HWP→L1,L1→L2,L2→mirror |400,550,150 mm | L:321 `design`. 앞/뒤 초점면·afocal 가정으로 둔 값. 현장 자/좌표 측정 아님 |
| pump: HWP→L1,L1→L2,L2→mirror |300,450,150 mm | L:335 `design`. 같은 범위 |
| 입구 PBS→셀 |50 mm | L:345 `nominal` |
| 셀→출구 PBS |100 mm | L:347 `nominal` |
| 각 DM→iris |100 mm | L:359 `nominal` |
| 각 iris→검출 렌즈 |215.8 mm | L:363 `nominal`; 셀→렌즈 전체900 mm를 만들도록 둠 |
| 각 검출 렌즈→PD |저장 길이 없음 | L:367 link에 `None`; DET가 focus 위치로 계산 |

`measured`인 PBS→DM484.2 mm도 “셀→DM584.2 mm에서 nominal 셀→PBS100 mm를 뺀 값”. 전체23 inch와 하위 구간484.2 mm의 측정 신뢰도는 다름. `Node.x/y/angle`·`Link.via`를 mm 좌표로 환산할 근거 없음.

## BPD·photodiode·전자 연결: 사양과 실물 구분

| 항목 | 확인 사실 | 적용 범위·출처 |
|---|---|---|
| PDB450A 기본 소자 | stock Si PIN, 직경0.8 mm, typical 최대 responsivity0.53 A/W | M PDFp.21/인쇄p.19. S3883 개조기와 다름 |
| gain별 bandwidth |10³/10⁴/10⁵/10⁶/10⁷ V/A →150/45/4/0.3/0.1 MHz | M PDFp.21/인쇄p.19. loaded 전달함수 실측 아님 |
| 명목 NEP |gain10⁵,DC–4 MHz에서3.3 pW/√Hz | M PDFp.21. stock 전체 detector 사양; S3883 교체 후 유지 보장 없음 |
| gain·load |표 gain은 high-Z; **50 Ω load에서는½** | M PDFp.21 각주. 현재 SA 변환은 이 factor를 명시적으로 반영하지 않음 |
| RF 출력 |SMA,출력 임피던스50 Ω, swing±4.6 V(gain10³–10⁶),±10 V(10⁷), high-Z 기준 | M PDFp.16/인쇄p.14. 50 Ω load에서 swing½ |
| monitor 출력 |SMA,+/− 개별 포트,220 Ω,10 V/mW@peak,DC–1 MHz;10 V high-Z/1.55 V50 Ω | M PDFp.16. RF OUTPUT gain과 동일시 금지 |
| 입력·전원 |FC removable;±12 V250 mA linear supply | M PDFpp.8·16. 실물 cable/adapter 계통 미확인 |
| 포화 조건 |CW9 mW@820 nm 표값. 균형 동작에서는 **두 입력 차이**가 포화 기준 미만 | M PDFpp.9–10·21. 현행 DET:295는 각 arm power>9 mW 개별 경고 |
| 선형성·광손상 |최대 선형성4 W/cm² 이하,수광면을 넓게 조사;20 mW damage threshold | M PDFpp.7·9. 숫자 그대로 stock 장비 지침 |
| S3883 수광면 |직경1.5 mm,유효1.7 mm²,TO-5 | H p.1 |
| S3883 명목 감도 |780 nm에서0.58 A/W;830 nm0.60 A/W | H p.1 표. **795 nm0.58은 사용자가 채택한 근사/보정값**이지 표에795 nm 점값으로 명시된 숫자 아님 |
| S3883 명목 응답 |VR20 V에서300 MHz,capacitance6 pF@1 MHz,NEP6.7×10⁻¹⁵ W/√Hz@peak | H p.1. photodiode 단독. BPD 전체4 MHz/NEP와 혼합 금지 |
| 사용자 장치 QE |사용자 확정값0.905,0.58 A/W@795 nm 채택 | 2026-10-13 `lab-context/LAB_CONTEXT.md`·`FWM_NOTATION_AND_LOGIC.md`; 상위 조사에서 직접 사용자 확인 계보 보관 |

### 논문과 구현 사이의 검출 효율 충돌

1. P25 p.5: optical loss5.5(2)%, system loss8.0(2)%, total13.5(2)%로 **가산 서술**.
2. README/기존 paper CSV: historical QE92%,after-cell94.5% → total86.94%. 이는 논문의 직접 QE 문장이 아님.
3. 현행 SABES: S3883 responsivity0.58를 사용 → QE90.4539%,after-cell94.5% → total약85.479%.
4. C:114의 “90.47 percent quoted in the paper”는 **확인한 Sim PDF에 없음**. PDB450A/S3883 명칭도 Sim 본문 전체에 없음.
5. A26 p.3는 PDB450A-SP2+S3883을 명시하나 **다른 저자·다른 실험**. 이것을 Sim2025 또는 사용자 장비의 근거로 쓰지 않음.

⟦AI⟧ 수정 후보: 위 세 효율을 서로 다른 reference/snapshot으로 보관. 사용자가 확인한 QE0.905는 재질문하지 않음. 현장 `ηOpt`가 셀 포함인지·후단만인지·어느 측정일인지 보존하고, atomic cell loss와 검출 loss 이중 계수를 피하도록 기준면 별 필드 필요.

## ⟦AI⟧ 원본 그림·논문과 현행 SABES의 수정 후보

채택 미정. 현재 구현 변경은 하지 않음.

| 우선도 | 관찰한 차이 | 업데이트에 필요한 정보 |
|---|---|---|
| 높음 | SQL bypass HWP·BS·2 flip mirrors가 D2에 있지만 SABES 전체63노드에 없음 | 원본 광로·확정 SQL/FWM 프로토콜 그대로 모드 스위치/노드로 추가. 원래 프로토콜 재질문 불필요 |
| 높음 | D1 Noise Eater는 etalon/GT 뒤. L:404는 EOM 앞 `ROUTE_SEED_TO_EOM`에 손실만 곱함. 렌더 node는 있지만 link 자체 없음 | 현장 최신 순서,서보 모델,구동기·PD·reference·대역 등. 도식 원본 순서 우선 보존 |
| 높음 | D1 CL은 MOPA 출력쪽인데 L:89 note는 입력. 모델은 단일M²로 대체 | input asphere/output asphere/cylindrical lens의 정확한 출입 순서·축·실물 기종·거리 |
| 높음 | L에 수광·광학 노드만 있음. 전기적 graph,RF cable·감쇠·isolator·bias·servo가 없음 | RF·DC·monitor·trigger·temperature 계통 별 포트 연결,임피던스·gain·대역 |
| 높음 | EOM paper AZ-0S5와 LM0904 PM-0S5 기종표기 상충 | 명판/구매문서로 동일 물건의 표기오차인지 교체인지 확인. 지금은 두 원문 보존 |
| 높음 | RF18 dBm 기본값,반파전압4.4 V nominal. LM0904는11.8 dBm optimum과 RF 경로 변경 기록 | 발생기 표시 power·cable/isolator loss·EOM port 실제 power 분리.12 dBm 주변 측정점 등록 |
| 높음 | nominal optical output이 논문 waist530/330과 다르게540/324 µm | 두 빔·두 축의 w0,z0,M²/astigmatism 측정. 공통1.08가 정의 차이 때문이라는 추측을 실험사실로 취급하지 않음 |
| 중간 | telescope 간격·초점위치는 ideal4F 가정. `path_m`와 source ABCD 전달이 별개로 구성됨 | 모든 렌즈 위치와 q parameter 연결. 렌즈 이동이 실제 전파·셀 모드에 반영되게 후속 설계 |
| 중간 | GT after filters node는 `seed_gtp_deg`를 소유하지만 polarization operation은 delivery trim 뒤 단 한 번 적용(B:427) | 실제 각각의 GT/HWP/QWP에 독립 Jones 전달·축·누설·loss 구분 |
| 중간 | pump QWP/HWP/GTP는 삽입손실만. 입력 polarization/amplifier behavior는 단순화 | 각 arm 측정 Stokes,부품 retardance·angle·zero reference |
| 중간 | pump scatter −60 dBm은 실측값처럼 설명되나 기본 숫자 출처 없음 | detector plane별 pump-only dark-subtracted DC와 RF PSD. 현재 가정값 표시 보존 |
| 중간 | probe와 conjugate가 항상 같은 seed transverse mode를 복제(DET:176) | 두 출력 beam profiles·divergence·센터·색수차·collection clipping 측정 |
| 중간 | detector가 이상 Gaussian focus에 자동 위치; nominal26/36 µm radius는 수광면 대비 매우 작음 | 실제 PD defocus·거리·spot radius·alignment/CMRR. focus 계산값을 장치 배치로 취급하지 않음 |
| 중간 | SA 한-pole,white/RIN scalar 합성;VBW·ENBW·detector/log average·AC/DC coupling 미구현 | 실제 SA 설정 snapshot·BPD transfer·dark/SQL spectral traces·50 Ω termination |
| 중간 | stock NEP·bandwidth·9 mW 포화값을 S3883 개조기에 재사용 | 개조 revision·bias·feedback capacitance·dark PSD·gain/load별 calibration |
| 중간 | per-optic 손실 nominal 합계보다 높은 총투과를 맞추는 residual1.004074>1 | passive optic transmission으로 쓰지 말고 aggregate calibration correction으로 별도 표기 |
| 중간 | calibration `verified`가 fitted 식별 가능성과 datasheet의 실물 일치까지 보장하는 인상 | 출처 문서·원시 측정·날짜·기종·오차·기준면·fit 방법·adoption 필드 분리 |
| 낮음 | ECDL frequency/control,thermal drift,locks/noise eater dynamics 없음 | 실제 수행 필요한 경우에만 closed-loop/control 모델 범위 확정 |

## ⟦AI⟧ 추가 확인 질문 후보

상위 조사에서 답을 찾으면 삭제. **자료에서 못 찾은 항목**이며, 실험자가 모른다거나 실험이 미수행이라는 뜻 아님.

1. 최종 갱신 기준일/실험 구성은 어느 날인가. 원본2025 도식,RF09/03,최근10월 셋업 중 기본 실행 snapshot을 명시할 필요.
2. EOM 명판의 전체 품번은 PM-0S5인가 AZ-0S5인가. 논문·슬라이드 차이가 오기/교체 중 무엇인가.
3. 최신 실측 테이블 치수·홀 pitch·원점,각 optic 중심(x,y,z),빔 높이,mirror normal,cell center와 lens distance 자료가 있는가.
4. 5724/5725 collimator 제조사·전체 품번·파장·fiber type·MFD·APC/PC·serial은 무엇인가. 그림의1.6/2.0 mm가 어떤 diameter 정의인가.
5. 각 렌즈 전체 품번·방향·유효구경·coating·focal tolerance·cylindrical 축과 각 측정 기준면은 무엇인가.
6. etalon 제조사·품번·개별scan·온도/tilt coefficient·heater/PID·sensor/lock 기록이 있는가.3개가 같은 사양·동일 leakage인가.
7. 최신 Noise Eater 위치·모델·actuator·servo PD·setpoint·loop bandwidth와 실제 켜짐 상태는 무엇인가.
8. BPD 개조 도면·실제gain·50 Ω/high-Z 종단·bias·두 monitor 연결·scope/SA port·RF cable loss·CMRR/dark transfer 측정은 어디에 있는가.
9. 셀 제조사·purity·window coating/tilt·heater/sensor 배치·cold-spot·자기차폐/코일 적용 상태는 무엇인가.
10. 셀 전·후·검출면별 power/waist/uncertainty의 최신 측정 표 또는 원데이터 경로는 어디인가. 이미 확정한 QE0.905·효율기준면·SQL 프로토콜은 유지.

## ⟦AI⟧ 검증 범위

P25 전9쪽 text extraction,Fig.1 시각 대조,D1/D2 흰 배경 시각 대조,M PDFp.21 표 시각 확인 완료. 실제 기본 광원/기하 계산 실행해 출력 숫자 확인. 코드 수정·원자 solver·시험 suite 실행 없음. 별도 실물 측정·장치 동일성 인증 없음. Allen2026 장비·수치를 사용자 장치에 가져오지 않음.
