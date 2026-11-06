<callout icon="🤖" color="gray_bg">
<b>용도.</b> 이 셋업(85Rb FWM 강도차 스퀴징 광원 + 87Rb Hanle 자력계)에서 잘 바뀌지 않는 조건·장비 기종·계수·측정 규약·잘 나왔던 조건을 모은 기준 문서. 이 프로젝트를 다루는 AI 에이전트는 실험일지보다 이 문서를 먼저 읽는다.

<b>우선순위.</b> 사용자 최신 지시 → 이 문서의 <code>확정</code> 항목 → 실험일지 최신 기록 → 과거 AI 분석. 기준 논문(Sim 2025)은 설계 기준선, 실험일지는 현재 구현 상태와 편차.

<b>태그.</b> <code>확정</code> 사용자가 확인한 조건 · <code>측정 MM-DD</code> 실험일지의 실측값(2026년) · <code>사양</code> 데이터시트 · <code>문헌</code> 논문·선배 자료 · <code>AI</code> 에이전트 계산·해석(실험자 의견 아님).

<b>범위 밖.</b> 최근 논의, 분석 진행 상황, 다음 실험 계획은 넣지 않는다. 그것은 아래 실험일지 DB와 Tasks에 있다. 값이 바뀌면 해당 줄만 고치고 근거 날짜를 갱신한다.

기본 장치 구성 기준일: <mention-date start="2026-10-16"/>. 이후 변경은 항목별 날짜·근거를 붙여 계속 갱신.

최종 갱신: <mention-date start="2026-10-21"/> — 사용자 Q1–Q26 직접 답변·제공 도식 반영.
</callout>

<table-of-contents />

## 0. 셋업 한눈에

- <b>목표:</b> 열증기 85Rb D1 이중-Λ FWM으로 probe·conjugate 강도차 스퀴징 광을 만들고, 이를 87Rb Hanle(zero-field level crossing, CPT) 자력계의 probe로 써서 같은 셀 입사 probe 파워의 coherent 광 대비 감도 개선을 보인다. quantum enhancement 주장은 동일 파워 개선이 검증된 경우에만 한다.
- **셋업 출처:** 심기성 선배 IDS 셋업과 Sim, Kim, Moon, Sci. Rep. 15, 7727 (2025) Fig. 1(a)를 계승. 08-12에는 기본 배치를 유지했다고 확인했으나, 이후 fiber EOM 교체·펌프 파이버 변경·계측기 이동을 반영. 논문은 설계 기준선, 아래 현재 기종은 실제 운용 구성. `확정 10-21` (사용자 Q1–Q26 답변)
- <b>주파수 생성:</b> 펌프와 시드는 같은 레이저에서 나온다. 펌프 일부를 fiber EOM으로 변조하고 에탈론 3매로 적색 −1차 측대역만 골라 시드로 쓴다. 별도 probe 레이저 없음.

```txt
[광원]  ECDL(795 nm, 자작 Littrow) - OI - HWP - glass(SAS 픽오프) - M - HWP - M - M - TA - 비구면 - CL(원통, f=70 mm)
        -> F260APC-780 -> 일반 APC–APC 파이버(현재 길이·MFD·PM 여부 미등록) -> 콜리메이터(1/e2 지름 2.0 mm)
        -> QWP, HWP -> 망원경 f300/f150 -> GTP -> PBS(시드와 합침)          [펌프 기준 waist 530 um, 정확한 z 미등록]
[시드]  펌프에서 PBS 픽오프 -> 선형편광자 -> fiber EOM(EOSPACE PM-0S5-10-PFA-PFA-780-UL, RF: E8257D)
        -> 에탈론 1 -> 에탈론 2 -> 에탈론 3 -> GTP -> 노이즈 이터(NEL03A/M, 삽입·ON·미최적화) -> 출력 파이버
        -> 콜리메이터(1/e2 지름 1.6 mm) -> QWP, HWP -> 망원경 f400/f150 -> GTP -> PBS   [시드 기준 waist 330 um, 정확한 z 미등록]
[셀]    펌프와 시드는 직교 선편광 -> 85Rb 셀(12.5 mm), 교차각 약 0.3 deg
        -> PBS(펌프 제거) -> D-mirror, iris로 probe/conjugate 분리 -> 렌즈(probe f75 / conj f100)
[검출]  플립 미러로 FWM(probe+conj) <-> SQL 기준광 전환 -> PDB450A(S3883 개조) -> RF: RIGOL DSA875 / Monitor +·-: Keysight DSOX1204A
[자력계] Hanle 셀(Pure 87Rb, 51.5 mm, Ne 50 Torr) + 원통형 mu-metal 3겹 + 3축 Helmholtz 코일
```

- 망원경·콜리메이터 지름·렌즈 초점거리는 셋업 도면 기준 값(Drive <code>calibration_sim2025.json</code>의 lab 항목). waist 530/330 µm는 1/e² 반경이며 지름이 아니다 <code>확정 10-01</code>.

- **위치·각도 기준:** 빔 교차점은 FWM 셀 내부이며 정확한 종방향 위치는 추가 측정 후 갱신. 약 23 inch(584.2 mm)는 **셀 출사면 끝→D-mirror** 거리. 출사면에서 빔 겹침을 확인하고, D-mirror 위치에서 probe–pump 간격을 재어 삼각형으로 각도를 구함. 셀 중심 기준 거리나 probe–conjugate 전체 간격으로 해석하지 않음. `확정 10-21` (사용자 Q1–Q26 답변)
- **SQL/FWM 전환 배치:** 아래 도식이 사용자 확인 기준. 파선 pump 우회 경로와 두 flip mirror로 같은 BPD 양 eye에 SQL을 넣음. 도식 좌표는 실측 mm 좌표가 아님. `확정 10-21` (사용자 Q1–Q26 답변)

<image source="file-upload://4006cba1-4fee-8120-95d2-00b2fdc21c5e" />

사용자 제공 FWM·SQL 전환 도식. 장치 기준 2026-10-16, 확인·제공 2026-10-21.

## 1. 측정 규약 — 다시 묻거나 재확인 과제로 만들지 않는다

### 1.1 SQL·FWM 비교 <code>확정 10-01, 10-08, 10-13</code>

- 샷 노이즈(SQL)와 FWM은 같은 셋업에서 <b>플립 미러로 전환</b>해 <b>같은 디텍터</b>로 잰다. 두 측정의 <b>총 검출 광파워를 항상 같게</b> 맞춘다. 전기적 셋업도 같으므로 전자잡음은 공통 배경으로 본다.
- SQL 기준광: 펌프를 첫 PBS에서 분기한 coherent 빔. 원자 셀을 통과하지 않는다. 두 번째 PBS에서 다시 나눠 두 eye에 넣는다. SQL 각 팔 파워는 probe·conjugate 파워의 평균 수준(현재 둘은 거의 같음).
- 이 파워 정합을 '확인 필요', '미확인', 'SQL calibration mismatch 가능성', 새 체크리스트로 다시 만들지 않는다. 기록에서 설명이 생략돼 있어도 미실시로 읽지 않는다. 상세: <mention url="https://app.notion.com/p/3f36cba14fee8123950de66d1e02411c" />
- SQL을 한 점에서만 실측한 스캔은 검출 광파워 비례로 환산한다: $L_{SQL}(P)=L_{SQL}(P_0)+10\log_{10}(P/P_0)$. 10/13 기준점: $P_0=81.3$ µW에서 −76 dBm(RBW 30 kHz) <code>측정 10-13</code>.
- 다른 비교 조건과 섞지 않는다: Hanle의 coherent/squeezed 감도 비교에서 맞추는 것은 <b>셀 입사 probe 파워</b>. 셀 전후 시드 잡음 비교, 과거 homodyne LO 조건, NMOR/Hanle 모드별 디텍터 차이도 각자의 맥락.

### 1.2 기호와 정의 — 10/16 NF 선택 반영

- **현재 잡음 표기:** Jasperse의 **NF**로 스퀴징 정도를 쓰고, 로그 표시는 **10 log₁₀(NF), Squeezing Degree [dB]**. 사용자 직접 요청 10-16. 9/17–10/13 기록의 정규화 잡음비 **S**와 아래 과거 경험식은 당시 기호로 보존. Jasperse의 intrinsic squeezing parameter S와 혼동하지 않음. 별도 표기가 없으면 전자잡음 미차감(raw).
- 전자잡음 차감은 선형 파워에서: $S_{opt}=(P_{FWM}-P_{Elec})/(P_{SQL}-P_{Elec})$. dBm 값끼리 빼지 않는다.
- <b>Gain 정의:</b> $G_s=P_{probe,\,pump\,on}/P_{seed,\,pump\,off}$, $G_c=P_{conj,\,pump\,on}/P_{seed,\,pump\,off}$. 모두 <b>셀 통과 후 디텍터 앞 같은 위치</b>에서 잰 값. 셀 입사 시드 기준이 아니다. 10/13부터는 디텍터 전압비로 계산(시드 35 mV 기준).
- **과거 경험식(9/17 랩미팅 8쪽, 10/13까지의 해석):** 아래 식과 수치는 이력으로 보존. 10/16부터 사용자 선택은 Jasperse NF 해석이며, 아래 a를 새 모델의 seed 고전 잡음 a와 동일한 계수로 취급하지 않음.

$$
S=(1-\eta)+\eta\left[\frac{(G_s-G_c)^2}{G_s+G_c}+a\max(G_s-1,0)\right]
$$

- **과거 경험식의 항:** N_ex=a(G−1)은 당시 seed 경로 초과잡음 항. 비공개 TeX 원고의 단일 이득식과 구분하며, 위 식의 두 항은 모두 η 괄호 안에 있음.
- **과거 점별 a 역산(10/13):** 최적점 부근 0.02, 저이득 0.17, Gc 우세 분기 0.14→28. 당시 경험식에서도 전 구간을 단일 a로 설명하지 못함. 현재 Jasperse 해석은 seed 고전 잡음 a와 shot-noise b를 구분. 실측 양팔 gain에서 역산한 intrinsic S·effective transmission Ta는 측정 기반 모델 해석값.
- 기호·논리 상세: Drive <code>원자광학연구실/Lab meeting/lab-context/FWM_NOTATION_AND_LOGIC.md</code>.

### 1.3 검출 효율 <code>확정 10-13</code>

<table fit-page-width="true" header-row="true">
<tr>
<td>
항목
</td>
<td>
값
</td>
<td>
근거
</td>
</tr>
<tr>
<td>
$\eta_{Opt}$
</td>
<td>
<b>0.781</b>
</td>
<td>
셀 앞 26.0 µW → 디텍터 20.306 µW(펌프 차단). 20.306 = 20.0 × 232.0/228.5(셀 틸팅 후 환산) <code>측정 10-06</code>
</td>
</tr>
<tr>
<td>
$\eta_{QE}$
</td>
<td>
<b>0.905</b>
</td>
<td>
Hamamatsu S3883, 0.58 A/W @795 nm <code>사양</code>
</td>
</tr>
<tr>
<td>
$\eta=\eta_{QE}\eta_{Opt}$
</td>
<td>
<b>0.707</b>
</td>
<td>
위 두 값의 곱
</td>
</tr>
<tr>
<td>
셀 뒤 광학계 투과
</td>
<td>
0.971
</td>
<td>
20.0/20.6 µW <code>측정 10-06</code>
</td>
</tr>
<tr>
<td>
펌프 OFF 시드의 셀 투과(창 포함)
</td>
<td>
≈0.805
</td>
<td>
20.92/26.0 µW, TPD +3, 113 °C. probe도 같다고 둠 <code>사용자 추정 10-15</code>
</td>
</tr>
<tr>
<td>
conjugate 셀 내 효율
</td>
<td>
1로 가정
</td>
<td>
85Rb D1 F=2→F′에서 blue 쪽으로 약 4 GHz 떨어짐 <code>사용자 가정 10-15</code>
</td>
</tr>
</table>

- 이 값들을 '미확정'으로 쓰거나 셀 뒤만의 효율로 바꾸지 않는다. 골드 레퍼런스의 검출 손실은 13.5%, 셀 뒤 손실 5.5% <code>문헌</code>.

### 1.4 파워 환산 계수

<table fit-page-width="true" header-row="true">
<tr>
<td>
항목
</td>
<td>
값
</td>
<td>
근거·주의
</td>
</tr>
<tr>
<td>
디텍터 전압 → 광파워
</td>
<td>
<b>0.08753 µW/mV</b> (20.0 µW ↔ 228.5 mV)
</td>
<td>
디텍터 출력을 오실로스코프로 읽은 값 <code>측정 10-06, 10-15</code>. 두 eye의 µW/mV 비율 같음 → seed·probe·conj에 같은 계수 <code>확정 10-13</code>
</td>
</tr>
<tr>
<td>
ND06A @795 nm
</td>
<td>
T = 0.2578 → <b>×3.879</b>
</td>
<td>
<code>측정 10-06</code> (08-11 교정값 0.2586). 펌프는 셀 직전에서 ND06A 거쳐 측정
</td>
</tr>
<tr>
<td>
ND03A @795 nm
</td>
<td>
T = 0.544
</td>
<td>
<code>측정 08-11</code>. 흡수형 ND는 공칭보다 NIR 투과가 높다 → 사용 파장에서 실측
</td>
</tr>
<tr>
<td>
Si 파워미터 헤드(공칭 500 mW)
</td>
<td>
300 mW 이상 비선형
</td>
<td>
355 mW 입사 → 310.9 mW 표시(−12.5%). W급은 ND로 표시값 200 mW 이하로 낮춰 읽는다 <code>측정 08-11</code>
</td>
</tr>
<tr>
<td>
파워미터 파장 설정
</td>
<td>
795 nm
</td>
<td>
780/755 nm 설정 상태로 잘못 읽은 사례 있음(06-26, 07-02)
</td>
</tr>
</table>

- 9/09–9/10 기록의 펌프 "86.2 × 4 mW"는 ND06A 공칭 계수로 환산한 값으로 보인다. 10/06 이후는 ×3.879 <code>AI 해석</code>.

### 1.5 주파수 규약

- <b>OPD(Δ):</b> 85Rb D1 <b>F=2→F′=3</b> 기준, <b>+ = blue</b> <code>확정 10-01</code>. 일상적으로 쓰는 "+0.9 / +1.0 / +1.06 GHz"는 모두 이 기준. GABES의 <code>opd</code>도 같은 정의.
- **OPD 설정·판독:** 현재 **SAS에 전자 lock을 걸지 않고**, 사용자가 오실로스코프의 SAS 피크를 눈으로 보고 위치를 맞춤. `확정 10-21` (사용자 Q1–Q26 답변) 현재 MHz 자리까지 실시간 판독하지 않음(10/13). 사후 SAS 피크 5개의 시간축 2차 보정에서 9/30 OPD 약 +1.06 GHz, 모형별 차이 ±0.03 GHz 수준; 선형 보정은 약 90 MHz 과대.
- <b>TPD(δ):</b> $\mathrm{TPD\,[MHz]}=f_{EOM}\,[\mathrm{MHz}]-3036$ <code>확정 10-13</code>. 85Rb 바닥상태 초미세 분리 $\nu_{HF}=3035.732$ MHz이므로 GABES의 lab TPD($f_{EOM}-\nu_{HF}$)와 0.268 MHz 차이. GABES 내부 δ = −(lab TPD) <code>AI 09-28</code>.
- 시드 = 펌프 − $f_{EOM}$ (적색 −1차 측대역).
- **과거 Task 표기:** 10/14 Task의 ‘pump red / probe=pump+3.04 GHz’는 위 현재 주파수 규약과 다른 과거 표기. 현재는 seed=pump−fEOM과 blue OPD를 기준으로 보관. Hanle의 상세 scan 기준·과거 Task 대응은 실제 탐색/구현에 필요할 때만 확인하고, 지금 반복 질문하지 않음.

<toggle>
Rb D1 SAS 기준 주파수 (85Rb F=2→F′=3 = 0, Steck rev 2.3.4)

| SAS 피크 | 상대 주파수 (MHz) |
| --- | --- |
| 85Rb F=2→F′=2 | −361.58 |
| 85Rb F=2→F′=(2,3) crossover | −180.79 |
| 85Rb F=2→F′=3 | 0 |
| 87Rb F=1→F′=1 | +1918.69 |
| 87Rb F=1→F′=2 | +2733.19 |

- D1 여기상태 분리: 85Rb 361.58 MHz, 87Rb 814.5 MHz(삼중선 폭이 약 2.25배 넓은 쪽이 87Rb). D1 동위원소 이동 77.583 MHz.
- D1 선 순서(주파수 증가 방향): 87Rb F=2 — 1.22 GHz — 85Rb F=3 — 3.04 GHz — 85Rb F=2 — 2.58 GHz — 87Rb F=1.
</toggle>

### 1.6 SA·데이터 기록 규약

- <b>SA 기본 설정:</b> span 0–3 또는 0–4 MHz, <b>RBW 30 kHz, VBW 300 Hz</b>, sweep 2 s, Input Att. 10 dB <code>측정 07-02, 09-30</code>. (05-28은 VBW 100 Hz.) 골드 레퍼런스의 대표 분석 주파수 150 kHz.
- <b>Trace 대응은 날마다 다를 수 있다.</b> 9/30: Tr1 = FWM, Tr2 = SQL, Tr3 = Elec <code>확정 09-30</code>. 파일마다 확인한다.
- 분석 주파수: 지금까지는 조건마다 수백 kHz 대역의 SA 최저점을 읽었다(10/13). 대역값은 dBm을 선형으로 바꿔 합의 비로 계산한다.
- <b>CSV 함정:</b> SA CSV의 주파수 열 머리글이 MHz여도 실제 값은 Hz로 저장돼 있다(09-09, 09-14 파일).
- **저장:** 대개 위 SA 설정으로 측정하고 USB에 저장. 정밀한 파일명/trace 식별 체계를 사용하지 않아 과거 파일 검색이 어려울 수 있음. 필요한 분석 시 실제 파일의 헤더·trace를 확인하며 이름만으로 날짜/조건을 추정하지 않음. `확정 10-21` (사용자 Q1–Q26 답변)
- 화면상 잡음 요동 약 ±0.5 dB(10/13).
- 그래프 색: FWM 검정 / Elec. 빨강 / SQL 파랑. 그래프 양식은 Drive <code>lab-context/GRAPH_STYLE_GUIDE.md</code>.

## 2. 장비 기종과 사양

### 2.1 광원 (ECDL–TA MOPA, 펌프 파이버)

<table fit-page-width="true" header-row="true">
<tr>
<td>
항목
</td>
<td>
기종·값
</td>
<td>
근거·주의
</td>
</tr>
<tr>
<td>ECDL</td>
<td>자작 Littrow ECDL, 795 nm. 외부 공진기 약 21 mm는 mode-hop 간격 약 7 GHz로 역산한 값.</td>
<td>연속 튜닝 PZT·전류에 따른 mode 선택은 08-11 기록. 현재 SAS lock 없이 수동 정렬. 내부 거울 수는 사용자 정정 기준 2개(08-13); 개별 부품 역할은 도식과 대조.</td>
</tr>
<tr>
<td>ECDL·TA controller</td>
<td>**Stanford Research Systems LDC501**</td>
<td>`확정 10-21` (사용자 Q1–Q26 답변). 과거 자료의 LDC-502를 현 기종으로 사용하지 않음.</td>
</tr>
<tr>
<td>
ECDL 동작점 (85Rb D1 F=2→F′)
</td>
<td>
25.550 °C / 143.43 mA / PZT offset 5.42 / 794.973 nm(vac)
</td>
<td>
<code>측정 08-11</code>. 그레이팅 조작으로 그 이전 모드 맵은 무효. 85Rb F=3→F′ 경유점: 25.50 °C / 145.65 mA / 6.52 / 794.983 nm
</td>
</tr>
<tr>
<td>
SAS 스캔 방향
</td>
<td>
스코프 화면 좌→우 = 주파수 증가
</td>
<td>
왼쪽 클러스터 85Rb F=2, 오른쪽 87Rb F=1. error signal 기울기 부호의 기준 <code>측정 08-11</code>
</td>
</tr>
<tr>
<td>
Mode-hop-free 범위
</td>
<td>
약 3.176 GHz, 85Rb F=2 전이선 왼쪽에서 mode hop
</td>
<td>
<code>측정 10-07</code>. OD·SAS 스캔에서 85Rb F=3 블록이 사라지거나 평평한 구간이 보이면 mode hop부터 의심
</td>
</tr>
<tr>
<td>
OI (ECDL 뒤)
</td>
<td>
투과 72% (기대 85–90%)
</td>
<td>
앞에 HWP가 없어 보정 불가, 미해결 <code>측정 08-11</code>. ECDL out 50.03 mW → TA 입력 29.2 mW
</td>
</tr>
<tr>
<td>TA (MOPA)</td>
<td>전류·온도는 매번 조금씩 변경. 실제 입출력 측정 기반의 휴리스틱 모델 사용.</td>
<td>앞·뒤 asphere 존재. chip·asphere 품번은 낮은 우선도의 미등록 정보. 08-11 이력: TEC25.000°C,1 A→355 mW,2 A→1.34 W 이상,무시드 ASE11.57 mW. 이 날짜의 운전점을 현재 고정값으로 사용하지 않음. `확정 10-21` (사용자 Q1–Q26 답변)</td>
</tr>
<tr>
<td>원통렌즈 CL</td>
<td>Thorlabs LJ1477L1 **추정**, f70 mm.</td>
<td>기종은 사용자도 추정으로 유지. TA 출력 측 비구면·CL의 상세 모델/배치는 기존 Sim 자료를 우선 활용. 회전 오차별 5%/31% 손실 수치는 과거 ⟦AI⟧ 계산이며 실측 아님.</td>
</tr>
<tr>
<td>
펌프 파이버 커플러
</td>
<td>
Thorlabs F260APC-780 (f 15.29 mm)
</td>
<td>
08-26에 F230APC-780에서 교체(빔 지름 약 3.3 mm에 정합) <code>측정 08-26</code>
</td>
</tr>
<tr>
<td>펌프 파이버</td>
<td>**현재 일반 APC–APC 파이버. OZ 고출력 파이버는 사용하지 않음.**</td>
<td>`확정 10-21` (사용자 Q1–Q26 답변). 현재 길이·MFD·PM 여부·AR 사양은 미등록. 08-12의 OZ SMJ-A3AHPC3A-780-5/125-3AS-2 AR2,2 m,MFD5µm,non-PM,air-gap/양단AR 사양은 과거 사용 이력으로만 보존.</td>
</tr>
<tr>
<td>
펌프 커플링 효율
</td>
<td>
320/1083 mW = 29.5% (목표 50%)
</td>
<td>
<code>측정 08-26</code>. 스퀴징 펌프로는 충분. 같은 배치의 과거 실적: 1.4 W → 500 mW
</td>
</tr>
<tr>
<td>
파장계
</td>
<td>
Advantest Q8326
</td>
<td>
진공 파장 표시. HI780 파이버 입력, 70–80 µW면 충분. RANGE 480–1000 nm, FUNCTION CHOP/LASER. 공용 장비라 상시 연결 아님 <code>측정 08-11</code>
</td>
</tr>
<tr>
<td>
레이저 스펙트럼 애널라이저
</td>
<td>
Keysight (파이버 입력)
</td>
<td>
analysis 기능으로 중심 파장·선폭. overexposure 표시되면 입력 감쇠 <code>04-02</code>
</td>
</tr>
</table>

### 2.2 시드 생성 (fiber EOM + 에탈론 3매)

<table fit-page-width="true" header-row="true">
<tr>
<td>
항목
</td>
<td>
기종·값
</td>
<td>
근거·주의
</td>
</tr>
<tr>
<td>Fiber EOM</td>
<td>**EOSPACE PM-0S5-10-PFA-PFA-780-UL**</td>
<td>논문의 AZ-0S5-10-PFA-PFA-780-UL에서 실제 기기를 교체함. `확정 10-21` (사용자 Q1–Q26 답변). 교체 날짜는 미등록. 광 입력 약10mW 운용 제한은 lab 기록이며 제조사 rating으로 승격하지 않음.</td>
</tr>
<tr>
<td>RF 소스</td>
<td>**Agilent/Keysight E8257D PSG**</td>
<td>현재 연결: **N–SMA 변환 커넥터→SMA 케이블→EOM**. `확정 10-21` (사용자 Q1–Q26 답변). 최신 확보 RF 설정은 09-03 실험/09-04 랩미팅 기록: 60A74-01(5.9–8.5GHz) 제거 후 발생기 표시11.8dBm 최적화. EOM 포트의 실측 RF power와 구분.</td>
</tr>
<tr>
<td>에탈론 필터</td>
<td>3매 직렬. Sim2025 명목 FSR17GHz/FWHM950MHz/peak T85%.</td>
<td>현재 개별 제조사·기종·정밀 실측 자료는 없음. 사용 기록 peak 약0.8,2매 carrier:seed 약0.5%,3매는 그보다 훨씬 작음. 세부사항은 Gisung Sim 자료로 보완, 필수 작업 때만 추가 확인.</td>
</tr>
<tr>
<td>
에탈론 온도 튜닝
</td>
<td>
약 1 GHz/°C
</td>
<td>
0.1 °C 단위로는 거의 안 움직임, 1 °C 단위로 조정 <code>측정 01-19</code>. 제어 정밀도 0.1 °C, 설정 대비 약 0.1 °C 오프셋 <code>12-30</code>
</td>
</tr>
<tr>
<td>에탈론 온도 컨트롤러 PID</td>
<td>**Thorlabs TC200**, 과거 PID P200/I5/D6</td>
<td>기종 `확정 10-21` (사용자 Q1–Q26 답변). PID 숫자는 05-21 기록이며 현재 설정을 새로 측정한 값은 아님.</td>
</tr>
<tr>
<td>
편광자
</td>
<td>
Glan-Thompson(GTP)
</td>
<td>
에탈론 3 뒤
</td>
</tr>
<tr>
<td>노이즈 이터</td>
<td>**Thorlabs NEL03A/M — 삽입·전원 ON·미최적화**</td>
<td>에탈론3→GT 뒤, 출력 파이버 앞. `확정 10-21` (사용자 Q1–Q26 답변). 현 SABES 적용 범위는 **투과율만**; noise suppression/feedback 성능을 가정하지 않음. 약93% 투과는 01-07 이력이며 현재 재측정값 아님.</td>
</tr>
<tr>
<td>
전송용 예비 파이버
</td>
<td>
APC–APC, AR 코팅, 5 m
</td>
<td>
스퀴징 광을 다른 셋업으로 보낼 때 사용 예정 <code>07-08</code>
</td>
</tr>
</table>

- <b>파이버끼리 직접 커플링(fiber-to-fiber) 금지</b> <code>확정 07-08</code>.

<toggle>
시드 경로 파워 체인 실측 사례 (정상 범위 감 잡기용)

| 단계 | 01-07 | 06-18 | 06-26 |
| --- | --- | --- | --- |
| EOM 입력 | 9.91 mW | 10.05 mW | 7.3 mW |
| EOM 출력 | 1.94 mW (19.6%) | 0.374 mW (3.7%) | 0.365 mW (5.0%) |
| 에탈론 1 후 | 0.210 mW (10.8%) | 160 µW | 77 µW |
| 에탈론 2 후 | 121.6 µW (57.9%) | — | 49 µW |
| 에탈론 3 후 | 67.9 µW (55.8%) | 8 µW (2+3 합산) | (2매 운용) |
| GTP 후 | 52.7 µW | — | 49 µW |
| 노이즈 이터 후 | 48.9 µW | — | — |
| 출력 파이버 후 (최종 시드) | 8.3 µW (17%) | — | 21 µW (42.9%) |

- 경험치: 에탈론 0차(EOM OFF) 최대 투과 약 0.8, 1차(EOM ON) 투과는 EOM 출력의 약 0.1. 시드는 셀 앞 약 8 µW면 사용 가능, 최근 운용은 3–6 µW.
</toggle>

### 2.3 FWM 셀과 가열

<table fit-page-width="true" header-row="true">
<tr>
<td>
항목
</td>
<td>
값
</td>
<td>
근거·주의
</td>
</tr>
<tr>
<td>FWM 셀</td>
<td>Pure85Rb 셀, 길이12.5mm. 미량87Rb(~0.1%) 가능성은 사용자 추정.</td>
<td>길이·종류는 기존 확인 유지. 공급사·직경·창·코팅 상세는 사용자 보류. 코르크 단열 cover는 과거 기록; **현재 FWM 자기 차폐 없음**. `확정 10-21` (사용자 Q1–Q26 답변)</td>
</tr>
<tr>
<td>빔 waist (1/e² 반경)</td>
<td>펌프530µm / 시드330µm, **1/e² 반경**.</td>
<td>`확정 10-01` 기준값 유지. 교차점은 셀 내부, 정확한 z 위치는 미정이며 추후 측정 자료 반영. waist/교차점을 정확히 셀 중앙으로 새 확정하지 않음. `확정 10-21` (사용자 Q1–Q26 답변)</td>
</tr>
<tr>
<td>교차각</td>
<td>약0.3°(논문0.32°), 셀 출사면에서 빔 겹침 확인.</td>
<td>출사면 끝→D-mirror 약23inch=584.2mm, D-mirror 위치의 **probe–pump 간격**으로 삼각형 각도 측정. `확정 10-21` (사용자 Q1–Q26 답변). 거리·겹침은 근사 계측. probe–conjugate 총간격과 구분.</td>
</tr>
<tr>
<td>펌프 제거</td>
<td>셀 뒤 PBS + D-mirror + iris. **PBS 전부 Thorlabs CCM1-PBS252/M**.</td>
<td>소광비3000:1은 Sim2025 서술, 현재 각 부품의 독립 실측값 아님. +0.9–1.0GHz에서 pump leakage가 작았던 기록 보존. PBS 기종 `확정 10-21` (사용자 Q1–Q26 답변)</td>
</tr>
<tr>
<td>렌즈·기타 편광 옵틱</td>
<td>렌즈는 Thorlabs plano-convex. 그 외 waveplate·편광 옵틱도 Thorlabs.</td>
<td>정확한 품번·각도·좌표는 사용자 보류. 제공 도식을 배치 기준으로 유지. `확정 10-21` (사용자 Q1–Q26 답변)</td>
</tr>
<tr>
<td>히터 컨트롤러</td>
<td>**자체 제작**, Notion Docs의 회로도_Heater Controller와 연결.</td>
<td>회로 근거: <mention-page url="https://app.notion.com/p/24f6cba14fee804eb94eea4d11b3652e"/>. **TP1은 센서 온도 전압**, 설정 노브 전압과 분리. 10월 초 보드 교체 이력. 사용자 확인상 FWM/Hanle는 같은 보드 형식이나 실장·보정 동일성을 별도 실측한 것은 아님. `확정 10-21` (사용자 Q1–Q26 답변)</td>
</tr>
<tr>
<td>센서·열 상태</td>
<td>센서는 셀 윗면 부착. cold spot 미측정.</td>
<td>FWM 자기 차폐 없음. `확정 10-21` (사용자 Q1–Q26 답변). 회로도의 센서는 Pt1000; 실제 보드의 실장 상세가 별도 확인된 것으로 확대하지 않음.</td>
</tr>
<tr>
<td>V–T 환산 기준</td>
<td>10/06–07 V–T 표·OD 자료는 **FWM 보드에서 측정**.</td>
<td>회로도 주석120°C=1.566V와 현장표120°C=2.262V는 서로 다른 보정 이력. 현장표는 10/12·10/16 랩미팅 부록 및 10/06 로그 참조. Hanle도 비슷할 것이라는 설명은 사용자 예상이며 별도 Hanle 실측 교정 아님.</td>
</tr>
<tr>
<td>
현 운용 온도
</td>
<td>
약 113 °C
</td>
<td>
설정 2.262 V(표상 120 °C)에서 TP1 2.227 V에 포화. OD 피팅 113.0 °C와 표 보간 112.9 °C 일치 <code>측정 10-07</code>. 설정을 올리면 더 오를 것으로 봄
</td>
</tr>
</table>

- <b>OD 측정법:</b> 약한 펌프 빔으로 스캔, 오실로스코프 raw CSV를 GABES OD 시뮬레이션과 수동 피팅해 온도 역산(10/07: 9.0 µW, 빔 지름 1.06 mm). 스캔에 mode hop이 섞이면 OD가 비정상으로 나온다. 포화·광펌핑을 피하려면 약 1 µW 이하가 이상적 <code>AI 권장</code>.

### 2.4 검출계

<table fit-page-width="true" header-row="true">
<tr>
<td>
항목
</td>
<td>
값
</td>
<td>
근거·주의
</td>
</tr>
<tr>
<td>Balanced detector</td>
<td>**Thorlabs PDB450A, 두 eye만 Hamamatsu S3883으로 교체. 회로는 기존 PDB450A와 동일.**</td>
<td>`확정 10-21` (사용자 Q1–Q26 답변). 수광 지름1.5mm. R795=0.58A/W·QE0.905는 사용자 채택 기준. 제품 데이터시트의0.58A/W 표값은780nm이며795nm 실측 사양으로 단정하지 않음.</td>
</tr>
<tr>
<td>TIA gain / 대역</td>
<td>**현재 gain 스위치10⁵ V/A**. 순정 매뉴얼에서 해당 설정 RF 대역4MHz.</td>
<td>gain은 사용자 직접 확인10-21. 매뉴얼의 high-Z gain은50Ω load에서 절반(실효5×10⁴V/A). 150/45/4/0.3/0.1MHz는 순정 gain별 명목 대역이며 개조품의 새 전달함수 측정 아님.</td>
</tr>
<tr>
<td>출력 포트 연결</td>
<td>RF→SA, Monitor +→오실로스코프(지금은 주로 probe), Monitor −→오실로스코프(지금은 주로 conjugate).</td>
<td>`확정 10-21` (사용자 Q1–Q26 답변). 양 monitor의 광학 arm 대응은 현 배치 기준. 실제 load/종단은 자료에 없는 값을 가정하지 않음.</td>
</tr>
<tr>
<td>오실로스코프</td>
<td>**Keysight DSOX1204A**</td>
<td>`확정 10-21` (사용자 Q1–Q26 답변)</td>
</tr>
<tr>
<td>스펙트럼 애널라이저</td>
<td>**RIGOL DSA875**</td>
<td>`확정 10-21` (사용자 Q1–Q26 답변). 2025 다른 셋업 재고 RSA5065N과 구분. 최근 SA setting은1.6절.</td>
</tr>
<tr>
<td>
포화·선형 한계
</td>
<td>
CW 포화 약 9 mW, 4 W/cm²
</td>
<td>
<code>사양</code>(순정 기준). 빔이 수광면을 넉넉히 덮게, 강하게 집속하지 않는다
</td>
</tr>
<tr>
<td>
NEP
</td>
<td>
사양값 적용 불가
</td>
<td>
사양 NEP는 순정 다이오드 기준. 개조품은 접합 용량이 달라 전자잡음은 반드시 빔 차단 실측 trace 사용
</td>
</tr>
<tr>
<td>
SQL 레벨 기준점
</td>
<td>
81.3 µW → −76 dBm (RBW 30 kHz)
</td>
<td>
<code>측정 10-13</code>. 이론값 $2eI\,R^2\,\mathrm{RBW}/50\,\Omega$(R = 5×10⁴ V/A) = −76.4 dBm으로 일치 → 검출 이득·QE 값이 서로 정합 <code>AI 검산</code>
</td>
</tr>
<tr>
<td>
전자잡음(빔 차단)
</td>
<td>
RBW 30 kHz에서 대략 −85…−95 dBm, 날·대역에 따라 다름
</td>
<td>
09-14: −85.9 dBm(0.5–1 MHz) · 09-30: −88.4 dBm(0.1–0.3 MHz) · 10-13 메모: 약 −95 dBm 스케일. 고주파로 갈수록 바닥이 올라와 스퀴징이 나빠 보임. SQL–Elec 간격 사례: 09-09 12.6 dB, 09-10 7.5 dB
</td>
</tr>
<tr>
<td>
자력계 쪽 디텍터(구 공용 셋업)
</td>
<td>
Hanle 모드 PDB210A/M, NMOR 모드 PDB250A/M (둘 다 0.56 A/W)
</td>
<td>
<code>2025-11-28</code>. 기종마다 TIA gain이 달라 샷 노이즈 비교 시 주의. 현재 FWM↔SQL 전환과는 무관
</td>
</tr>
</table>

### 2.5 자력계(Hanle) 쪽

<table fit-page-width="true" header-row="true">
<tr>
<td>
항목
</td>
<td>
값
</td>
<td>
근거·주의
</td>
</tr>
<tr>
<td>
Hanle 셀
</td>
<td>
Pure 87Rb, 길이 51.5 mm, Ne 50 Torr
</td>
<td>
<code>확정 10-14</code>. 소량 85Rb 혼입 가능성은 남김. 셀 중심 위치·코일·실드·구동 회로는 별도 언급 전까지 유지. 셀 중심 = 코일 응답 최대이고 외부장 영향이 평탄한 구간
</td>
</tr>
<tr>
<td>목표 전이</td>
<td>**87Rb D1 F=2→F′=1·2 사이**</td>
<td>사용자10-14·10-16 직접 방향 변경. 과거04-22의 F1→F′1 CPT 계획과 F1→F′2의 EIA/광펌핑 해석은 과거 전이에 대한 참고이며 현재 목표로 사용하지 않음. 상세 scan/주파수 대응은 사용자 보류.</td>
</tr>
<tr>
<td>기준일 삽입·가열 상태</td>
<td>**10/16 셀 삽입·히터 연결 완료, 미가열 상온.**</td>
<td>EOM3.039GHz,seed35mV,온도계1.722V. 10/21 사용자Q1에서10/16을 기준 snapshot으로 선택.</td>
</tr>
<tr>
<td>10/16 셀 전후 파워</td>
<td>**앞27.17µW / 뒤26.53µW**</td>
<td>사용자가10/21에 기존 앞뒤 순서가 바뀌었다고 정정. 파워미터 측정. 이전26.53/27.17 표기는 정정 전 값으로만 보존.</td>
</tr>
<tr>
<td>
운전 방침
</td>
<td>
스퀴징 운전점을 유지한 probe로 Hanle 탐색
</td>
<td>
Hanle만을 위해 레이저를 다른 전이로 먼저 옮기지 않는다 <code>확정 10-14</code>. 선편광, 종방향 자기장을 영점 부근에서 스캔
</td>
</tr>
<tr>
<td>
Hanle 셀에서의 빔 크기
</td>
<td>
coherent probe 530 µm / squeezed probe 330 µm
</td>
<td>
평행광이라 전파 중 퍼짐 무시하고 FWM 셀 값을 그대로 쓴 가정. 새로 잰 값 아님 <code>확정 10-14</code>
</td>
</tr>
<tr>
<td>
µ-metal 실드
</td>
<td>
주문제작 원통형 3겹, 지름 약 19 cm × 길이 약 40 cm
</td>
<td>
<code>04-28</code>부터 사용. reference 빔용 플립 미러 부품과 간섭 있음
</td>
</tr>
<tr>
<td>
코일
</td>
<td>
3축 Helmholtz
</td>
<td>
축이 완전히 직교하지 않아 cross-coupling 있음, 보정 행렬은 미작성 <code>03-26</code>. 빔 입사 구멍을 z = 0으로 둘 때 코일 응답 최대 z 약 20–23 cm <code>04-28</code>
</td>
</tr>
<tr>
<td>자기장 교정</td>
<td>4/28 기록: 셀 중심 fluxgate, scope 표시÷10≈µT.</td>
<td>**이전4/03**: 1Hz에서FG2–8Vpp→88–354mV 선형(4/28 실드 교체 전). **4/28**: 약24.3µT↔scope약243표시.100Hz 응답 감소는 센서·coil 기여가 분리되지 않음. 두 날짜의 결과를 같은 최신 교정으로 합치지 않음. 원기록: <mention-page url="https://app.notion.com/p/3366cba14fee81f5b923d1c2ec7430bc"/> · <mention-page url="https://app.notion.com/p/3506cba14fee816787abcd88ae11e101"/>.</td>
</tr>
<tr>
<td>
셀 온도 참고
</td>
<td>
ALTA 2026 포스터 Hanle 측정 73 °C
</td>
<td>
재개 출발값(실제 사용 온도는 따로 기록). 120 °C는 선폭이 너무 넓어짐
</td>
</tr>
<tr>
<td>
문헌 기준
</td>
<td>
Yu et al., PRA 81, 023416 (2010)
</td>
<td>
Ne 50 Torr, 87Rb F=2→F′=1 선폭 1.2 mG. F=1→F′=1 CPT 투과 피크 contrast 약 0.5%
</td>
</tr>
</table>

- 실드 교체(04-28) 이전의 코일 영점 값(X 0.045 / Y 0.077 / Z 0.095 A, Z offset 0.52 V 등)은 현 실드에서 다시 확인해야 한다.

### 2.6 전원·테이블·미등록 기기

- **현재 테이블 접지 없음.** SA를 detector와 **같은 테이블로 옮긴 후**, 과거의 큰 전기적 잡음은 들어오지 않음. 접지를 완료한 결과로 해석하지 않음. `확정 10-21` (사용자 Q1–Q26 답변)
- 멀티탭 **그린세이프**, 전원 필터 **자체 제작**. 배선은 과거 SU(1,1) 도면과 일부 다르며 상세 port/접점은 미등록. 기존 재고·다른 셋업 배전도를 현재 회로에 자동 적용하지 않음. `확정 10-21` (사용자 Q1–Q26 답변)
- **Waveform generator 모델 미등록.** 해당 기기의 제어·연결·주파수 응답을 구현/분석할 때 사용자에게 다시 알리고 질문. 지금 추가 확인 과제로 재요청하지 않음. `확정 10-21` (사용자 Q1–Q26 답변)

## 3. 실험적으로 잘 나왔던 조건

<table fit-page-width="true" header-row="true">
<tr>
<td>
날짜
</td>
<td>
OPD
</td>
<td>
TPD (f\_EOM)
</td>
<td>
펌프 (셀 앞)
</td>
<td>
시드
</td>
<td>
셀 온도
</td>
<td>
$G_s$
</td>
<td>
결과
</td>
</tr>
<tr>
<td>
Sim 2025 (골드)
</td>
<td>
+0.9 GHz
</td>
<td>
−8 MHz (논문 표기)
</td>
<td>
600 mW
</td>
<td>
8 µW
</td>
<td>
121 °C
</td>
<td>
약 15
</td>
<td>
−7.8 dB 직접 (40–500 kHz), probe/conj 111/109 µW
</td>
</tr>
<tr>
<td>
05-28
</td>
<td>
+0.9 GHz
</td>
<td>
0 (3.036 GHz)
</td>
<td>
400 mW
</td>
<td>
27.8 µW
</td>
<td>
113 °C
</td>
<td>
13.5
</td>
<td>
bright FWM (스퀴징 미측정). 같은 날 seedless bright FWM 약 30 µW 관측
</td>
</tr>
<tr>
<td>
09-10
</td>
<td>
+1.0 GHz
</td>
<td>
+3.5 MHz
</td>
<td>
344.8 mW
</td>
<td>
4.0 µW
</td>
<td>
118 °C
</td>
<td>
7.5
</td>
<td>
raw −3.45 dB / Elec 차감 −4.52 dB (0.13–0.30 MHz). 에탈론 투과 중심을 FWM 이득 지점으로 옮긴 뒤
</td>
</tr>
<tr>
<td>
09-30
</td>
<td>
+1.06 GHz
</td>
<td>
+3 (3.039 GHz)
</td>
<td>
—
</td>
<td>
3.7 µW
</td>
<td>
—
</td>
<td>
13.4
</td>
<td>
−0.66 dB (0.1–0.3 MHz). EOM을 내리면 gain은 118.7까지 오르나 잡음은 SQL 위 +10 dB
</td>
</tr>
<tr>
<td>
10-13
</td>
<td>
약 +1.0 GHz
</td>
<td>
+2 (3.038) / +3 (3.039)
</td>
<td>
294.8 mW
</td>
<td>
3.06 µW (35 mV)
</td>
<td>
약 113 °C
</td>
<td>
13.6 / 11.9
</td>
<td>
−3.0 dB 직접 측정 / −3.42 dB(SQL 환산 추산)
</td>
</tr>
</table>

- 시드 파워는 디텍터 앞 펌프 OFF 값. 골드 레퍼런스 TPD −8 MHz가 랩 TPD 규약에서 어느 부호인지는 미확정(심기성 선배 240419 덱은 'TPD 8 MHz'로 표기).

<b>반복해서 확인된 경향</b>

- <b>Gain 최대점 ≠ 스퀴징 최적점.</b> 스퀴징은 $G_s\gtrsim G_c$ 분기(TPD ≥ 0)에서만 SQL 아래. $G_c$ 우세 분기(TPD ≤ −3)는 같은 gain에서도 5–13 dB 나쁘다(10/13).
- <b>+0.9–1.06 GHz(blue)가 실험적 최적 영역</b>(심기성 선배가 여러 OPD를 실험으로 탐색). GABES 최적점 −2.1 GHz에서는 gain은 크게 나오지만 펌프 누설 때문에 스퀴징 측정 불가(07-02).
- <b>FWM이 서는 조건:</b> 셀 약 112–120 °C 이상, 펌프 약 300 mW 이상(선배 권장 500 mW), 에탈론이 올바른 −1차 측대역을 통과. 125 °C 이상은 초과잡음 의심(07-02).
- <b>날마다 다시 맞출 것:</b> 셀 산란(셀 틸팅·이동으로 투과 최대), 에탈론 온도(심기성 선배 240419 덱 '날마다 바뀌는 주요 요인').
- OPD를 MHz 단위로만 바꿔도 FWM 잡음이 크게 달라질 수 있다(10/13 관찰).
- 시드는 크게 쓸 필요 없다. 이상적인 경우 seed/probe를 작게 쓰는 편(06-26 메모).

## 4. 알려진 잡음원과 함정

<table fit-page-width="true" header-row="true">
<tr>
<td>
현상
</td>
<td>
원인·상태
</td>
<td>
기록
</td>
</tr>
<tr>
<td>
380 kHz 피크
</td>
<td>
실험실 열화상 CCTV의 IR이 디텍터로 유입(광학적 잡음). 해결
</td>
<td>
2025-12-03
</td>
</tr>
<tr>
<td>
1 MHz 이하 자글자글한 전선 잡음
</td>
<td>
디텍터 monitor output 단자(RF output 인접) 쪽에서 유입. 해결 → 배선 바꿀 때 재발 주의
</td>
<td>
07-09
</td>
</tr>
<tr>
<td>0.2–0.6 MHz 정체불명 잡음</td>
<td>당시 전원 노이즈 필터3단에서도 남았고 케이블·커넥터 교체로 일부 개선.</td>
<td>06-26 역사 기록. 필터 효과가 없었다는 관측만으로 접지 기원을 배제하지 않음. 현재 SA 이동 후의 큰 잡음 감소는 아래 별도 상태와 구분.</td>
</tr>
<tr>
<td>
약 1.9–2 MHz 피크(+pedestal)
</td>
<td>
펌프에는 없고 시드에 있음. 셀 통과 시드 단일 eye에서 1.92 MHz 피크. 광원·전자계 기원 잡음을 에탈론 투과 기울기가 FM→AM으로 바꾸는 경로가 유력: PZT로 레이저 주파수를 에탈론 중심에 두면 시드의 피크가 거의 사라짐. 저주파 초과잡음과는 별개. 교수님·선배 의견: 어쨌든 제거하는 게 좋다
</td>
<td>
07-13, 09-14
</td>
</tr>
<tr>
<td>
셀 통과 시드의 넓은 대역 초과잡음
</td>
<td>
단일 eye 12 µW에서 shot noise보다 약 7 dB 위($F_{det}$ 약 9.5, 0.1–0.5 MHz)
</td>
<td>
09-14
</td>
</tr>
<tr>
<td>
펌프 누설
</td>
<td>
공명에서 먼 OPD(−2.1 GHz)에서 PBS를 지나 디텍터로 유입. +0.9 GHz에서는 거의 없음. 판별: 시드 차단·펌프만 켠 상태의 잡음
</td>
<td>
07-02, 06-26
</td>
</tr>
<tr>
<td>
시드 파워는 나오는데 FWM이 안 생김
</td>
<td>
에탈론이 원치 않는 측대역을 통과. 파워미터로는 구분 불가 → 에탈론 방향부터 확인
</td>
<td>
06-24
</td>
</tr>
<tr>
<td>
시드 파워 드리프트
</td>
<td>
fiber EOM 입력 편광 문제로 추정(직전 HWP 미조정). 매 측정 35 mV로 재조정
</td>
<td>
10-13
</td>
</tr>
<tr>
<td>
Gain의 수 분 단위 드리프트
</td>
<td>
펌프 파워·OPD는 불변. 원인 미확정
</td>
<td>
09-30
</td>
</tr>
<tr>
<td>
오실로스코프와 파워미터의 probe/conj 비 불일치
</td>
<td>
795 nm 외 빛 혼입 가능성(미확정)
</td>
<td>
06-26
</td>
</tr>
<tr>
<td>
TA 출력이 낮게 읽힘
</td>
<td>
Si 파워미터 헤드 포화. ND로 낮춰 읽으면 정상
</td>
<td>
08-11
</td>
</tr>
<tr>
<td>
TA 출력 실제 저하(1.4 W → 700 mW)
</td>
<td>
mode hop 대응 중 ECDL 내부 거울을 바꾼 여파. 동작점 복구로 회복. TA 소자 열화 아님
</td>
<td>
08-11
</td>
</tr>
</table>

## 5. 표준 작업 절차 요약

### 5.1 에탈론 정렬

1. 에탈론 입력 파워를 먼저 기록한다(local maximum 판정 기준).
2. 각도: 한 축씩 노브를 크게 돌려 크기가 거의 같은 큰 피크 2개를 찾는다. 그 사이의 아주 작은 피크(큰 피크의 약 1/100)가 그 축의 0(수직 입사). x·y를 번갈아 반복해 x = y = 0. 이 단계 투과 약 50%.
3. 온도: 투과가 오르는 방향으로 조정. 최적이 실온 아래면 10 °C 이상 올려 FSR을 한 칸 넘긴다. 잘 맞으면 80–90%.
4. EOM OFF에서 0차 최대를 잡은 뒤, EOM ON에서 온도로 1차 피크로 옮긴다. 1차에 맞으면 EOM을 껐을 때 투과가 거의 사라진다.
5. +1/−1차 구분은 파장계. 실용적으로는 뒤쪽 정렬이 정상일 때 FWM이 바로 보이면 맞는 측대역.
6. 마지막으로 에탈론 투과 중심을 실제 FWM 이득이 나는 지점에 맞춘다(투과 최대가 곧 FM→AM 변환 0은 아님).

- 출처: 정한솔 선배 팁(05-15), 12-30, 01-07, 06-24, 09-09 실험일지.

### 5.2 FWM이 안 나올 때 점검 순서

펌프·시드 공간 중첩 → 시드 주파수(SAS와 함께 보며 에탈론 투과 위치가 펌프에서 3.036 GHz인지) → 에탈론 측대역 방향 → 셀 온도 → 펌프 파워 (01-19, 06-24).

### 5.3 파이버 커플링 (김희우 선배 조언, 08-24)

- 한 위치의 세기 분포만이 아니라 전파 조건(평행광인지 수렴인지)까지 맞춘다. Thorlabs 콜리메이터에는 평행광을 넣는다.
- TA 입력 방향 역방향 ASE도 평행광이어야 한다(타원인 것은 정상). 역 ASE가 맺는 초점거리 D로 비구면 축 위치 오차 $\delta\approx f_a^2/D$를 읽는다.
- '쓸 만한' 수준에서 멈추지 말고 최대 효율점에서 잠근다. 최대점 근처는 흔들려도 손실이 작다.

### 5.4 루비듐 셀 창 정리 (01-02)

- 토치로 셀 양 옆면을 번갈아 1초 이내로 짧게 가열해 루비듐을 중앙으로 모은다. 창에 열을 집중하지 않고, 끝나면 충분히 식힌다.

## 6. 소프트웨어와 참고 문서 위치

- <b>GABES</b>(Generic Atomic Bloch Equation Solver): GitHub <code>Shake2313/fwm-squeezing-app</code>, 로컬 <code>G:\내 드라이브\Coding\GABES - Generic Atomic Bloch Equation Solver</code>. 장비 계수 파일 <code>calibration_sim2025.json</code>(Drive)은 항목마다 provenance(datasheet / paper / lab / nominal / fitted)가 있고, nominal은 자리표시값이다. 사용자 판단(10/15): 모델이 아직 완전하지 않아 실험 레퍼런스 근거로 주장하기엔 빈약.
- <b>Drive</b> <code>원자광학연구실/Lab meeting/</code>: <code>AGENTS.md</code>(랩미팅 자료 작성 규칙), <code>README.md</code>, <code>lab-context/LAB_CONTEXT.md</code>, <code>FWM_NOTATION_AND_LOGIC.md</code>, <code>NOTION_SOURCE_INDEX.md</code>, <code>SIM_LAB_MEETING_ARCHIVE.md</code>, <code>GRAPH_STYLE_GUIDE.md</code>. 랩미팅 자료 규칙은 이 문서에 중복하지 않는다.
- **심기성 선배 자료:** Drive `G:/내 드라이브/Docs/Gisung Sim/`의220729–240627 랩미팅38개. 최종 최적 조건은240419덱. 논문: Sim, Kim, Moon, Sci. Rep.15,7727(2025),DOI10.1038/s41598-025-86479-w. TA chip/asphere·etalon 상세는 이 자료에서 먼저 보완. 부품이 확인되지 않으면 TA는 사용자 선택에 따라 측정 기반 휴리스틱 범위를 유지.
- <b>데이터시트(Drive):</b> <code>PDB450A-Manual.pdf</code>, <code>s3883 광다이오드.pdf</code>.
- <b>Notion:</b> 실험일지는 아래 Experiment Log DB(이 프로젝트로 필터된 보기), 계획은 Tasks DB, 반복 측정 방법은 Docs › Lab memo.

## 7. 보류 정보와 재질문 시점

사용자 방침10-21: 사소한 세부나 추가 자료 조사로 확인 가능한 항목은 보류하고, **정말 중요한 정보만 필요한 작업 시점에 다시 질문**. 미등록을 미실시 실험이나 실험자 무지로 해석하지 않음.

- **사진·기하:** 최신 테이블 사진은 추후 제공 예정. 교차점은 셀 내부지만 정확한 위치는 추후 측정. 제공된 SQL/FWM 도식은 현재 배치 근거이며 실측 좌표 아님.
- **Waveform generator:** 모델이 필요한 제어/코일/주파수 응답 작업을 시작할 때 다시 알리고 질문.
- **일단 보류:** 렌즈/기타 Thorlabs 부품 상세, FWM/Hanle 셀 공급사·창·mount, Hanle 상세 scan/기준 주파수. chip/asphere·etalon·coil 측정 기록은 Gisung Sim/랩미팅/Notion 기존 자료 우선.
- **낮은 우선도:** Q17의 전압 판독 port/추가 환산 세부는 현 단계에서 더 캐묻지 않음. USB 파일은 엄격한 이름 체계가 없어 필요 분석 때 실제 파일로 대응.
- **이미 해결:** 기준일10/16, PM-EOM 실제 교체, EOM 연결·RF 최신 기록, NEL 삽입ON/투과율만 반영, LDC501·DSOX1204A·DSA875·TC200·CCM1-PBS252/M, 현재APC–APC fiber, manual SAS, BPD gain/출력 연결, SA 같은 테이블·미접지, 센서 윗면·cold spot 미측정, 셀 전후 power 정정.
- **기존 확인 유지:** SQL/FWM 파워 정합·공통 전기 구성,530/330µm 반경,QE0.905,ηOpt0.781/ηtotal0.707/after-cell0.9709. 값이 바뀌면 사용자가 알리기로 함. 골드 레퍼런스의 당시 RF/TPD는 문헌 재현에 필요할 때만 원자료로 대응.

## 변경 이력

- **2026-10-21 사용자 확인 반영:**10/16기준의 지속 갱신,장비/배선/파이버/Noise Eater/기하/온도 기준 보정,제공 도식 추가,Hanle power 순서 정정. 출처는 현재 사용자 Q1–Q26 답변 및10/16직접 NF 요청. 위 장치·방법 갱신은 사용자 확인, 미등록 세부는 필요 시점 재질문으로 보류.
- **⟦AI⟧ 추가 검수:**기존 조사·직접 사용자 발언과1회 대조.584.2mm 기준면·FWM/Hanle 온도 교정·문헌/개조 detector 사양·4/03/4/28 coil 자료의 범위를 분리. 새로운 suppression 성능·Hanle 교정·실물 좌표를 추정값으로 채우지 않음.

- 2026-10-16 초판 (Claude). 출처: 이 프로젝트 실험일지 41건(2025-11-28 – 2026-10-14), Tasks 2건, Lab memo 'SQL 측정 방법에 대해', Drive <code>lab-context</code> 문서·<code>calibration_sim2025.json</code>·데이터시트, 이전 대화 기록.

<p />
