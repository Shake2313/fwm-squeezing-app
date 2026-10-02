"""Build the Korean report from retained diagnostic results, without model edits."""
import csv
import json
from pathlib import Path

OUT = Path(__file__).resolve().parent
ROOT = OUT.parent.parent


def read(name):
    return json.loads((OUT/name).read_text(encoding="utf-8"))


def link(name, path, line=None):
    target = (ROOT/path).as_posix()
    if line is not None:
        target += f":{line}"
    return f"[{name}](<{target}>)"


def main():
    data, summary, extra = read("results.json"),read("summary.json"),read("additional_checks.json")
    mix = read("mixing_sensitivity.json") if (OUT/"mixing_sensitivity.json").exists() else None
    manifest=read("source_manifest_before.json")
    table=list(csv.DictReader((OUT/"comparison.csv").open(encoding="utf-8-sig")))
    c=summary["fits"]["exact_corrected"]
    u=summary["fits"]["exact_uncorrected"]
    def figure(name):
        return f"![{name}](<{(OUT/(name+'.png')).as_posix()}>)"
    sections=[]
    sections.append(fr"""# TPD–Gain 실측·GABES 비교 진단

작성일: 2026-10-01, Asia/Seoul. **진단·비교·보고만 수행. 기존 코드·설정·보정식·계수·기준 데이터 유지.**

최종 지정 조건: **118 °C, Pump 380 mW, OPD +1.06 GHz, 셀 입력 Seed 3.7 µW**. 나머지 조건만 현 GABES Default. 최초 첨부의 OPD/Seed 해석은 사용자의 후속 정정으로 대체. +0.9 GHz·8 µW를 본 비교에 사용하지 않음.

## 1. 판단

- 동일 EOM 10점을 직접 계산한 현행 보정 Gain: **6.351→11.275**. 실측: **13.432→118.649**. 증가 방향 일치, 증가비 **1.775배 대 8.833배**. 절대값과 상승폭 불일치 유지.
- **단일 출력 배율로 설명 부족**: 전체 10점 최소제곱 적합 RMSE **25.696**. 실측/모델 비가 **2.115→10.523**으로 변함.
- **Affine는 상당 부분 근사 가능**: \(G_{{meas}}≈23.027G_{{sim}}-141.568\), RMSE **4.177**, \(R^2=0.9870\). 따라서 “곡선이 비선형이므로 모든 선형 보정 실패”라는 주장은 성립하지 않음. 다만 저 Gain 점의 상대 잔차와 구조적 잔차 유지, 큰 음의 offset은 전역 물리 보정에 부적합. 측정 분산이 없어 통계적 기각 불가.
- **공통 온도 ±5 °C만으로 크기·형태·피크 동시 설명 불가**. 113–123 °C에서 3.030 GHz Gain **4.456–34.654**, 증가비 최대 **1.784배**. ±5 °C 민감도 범위 밖 탐색에서 마지막 점을 맞춘 **127.886 °C**도 첫 점 **85.467**을 예측해 전체 곡선과 불일치.
- **피크 일치 전제도 재확인 필요**. 현행 보정, 118 °C의 밀집 계산 피크 **EOM 3.0220 GHz, Gain 14.576**. 사용자의 약 3.030 GHz 관측과 약 **8 MHz** 차이. 실측 표는 3.030에서 끝나므로 감소 쪽 정량 검증 불가.
- 현행 보정은 Gain을 낮추고 이 구간의 증가비도 줄임. **보정 해제만으로 해결되지 않음**: 해제 Gain **80.962→219.501**, 원값 RMSE **81.819**로 더 큼. 여기서 해제는 `gain_closure_enabled=False`만 뜻함. 기존 residual 0.74는 유지.
- **즉시 계수·함수 변경 근거 부족**. 단일 내부 계수 재적합도 잔차와 피크를 해결하지 못함. 먼저 측정면·시간 drift·밀도/실효 온도·기준점 부호 확인. 별도로 Fast 적응 보간 문제 확인됨. 이번 10점에서는 영향 작지만 확장 곡선에 큰 오류 가능.

주요 원인 후보: 검증되지 않은 보정의 조건 외삽, 생산 모델의 1차원 Doppler/축약 준위·모드 가정, 실효 온도/밀도와 시간 drift, 불완전한 pump-off 기준면 교정. 이번 자료로 기여도 분리 불가. Efficiency·Noise를 Gain 맞춤에 사용하지 않음.

## 2. 사용 조건·버전

| 항목 | 사용값 | 근거·규약 |
|---|---:|---|
| 표시 온도 | 118 °C = 391.15 K | 사용자 지정; 실효 온도 미교정 |
| Pump | 380 mW | 사용자 지정; 과거 320/344.8 mW 미사용 |
| OPD | +1.06 GHz | 사용자 후속 지정; ⁸⁵Rb D1 F=2→F′=3 기준, blue positive |
| 셀 입력 Seed | 3.7 µW | 사용자 후속 지정; 모든 본 계산에서 공통 |
| 셀 길이 | 12.5 mm | Default |
| Pump waist | 530 µm | 1/e² **강도 radius**; diameter 1.06 mm |
| Seed waist | 330 µm | 1/e² 강도 radius; diameter 0.66 mm |
| Conjugate collection waist | 330 µm | 별도 입력 없음 → Seed와 같다고 가정 |
| 교차각 | 0.32° = 5.585 mrad | Default; 전파 mismatch/axial crossing에 적용 |
| 원자·준위 | pure ⁸⁵Rb D1, 축약 4준위 double-Λ | 생산 유한 Seed Floquet 모델 |
| Raman branch | −1, red sideband | \(\nu_{{seed}}=\nu_{{pump}}-f_{{EOM}}\) |
| density | 1.80740×10¹⁹ m⁻³ | CRC pressure/\(k_BT\), 118 °C |
| transit reset | γₜ/2π = 100 kHz | 독립 교정 없는 inherited Default |
| ground collision dephasing | γcoll/2π ≈ 2.414 kHz | 현재 density 기반 계산 |
| residual line strength | 0.74 | 과거 기준 residual; 아래 mixing 보정과 별개 |
| macroscopic normalization | 1/12 | 해당 branch 코드값; residual 포함 coupling 0.74/12 |
| overlap/polarization/Zeeman 상대 penalty | 각각 1 | Default; 추가 적합 없음 |
| Semi-empirical mixing factor | 0.5594938027, ON | 현행 Fast/Balanced |
| detection efficiency | 86.94% | Default; QE 92%, legacy post-cell loss 5.5% |
| Excess noise | constant N=0, slope=0 | Default; 내부 pump-scatter diagnostic κ=0.1 |
| residual EOM carrier/other sidebands | 각각 0 µW | Default; 생산 원자 응답에는 unapplied |
| 실제 전파 | Ultra propagation, 64 segments | 모든 현재 seeded tier 공통; axial Gaussian crossing + 근사 pump budget + 최종 saturation |
| UI Default solver | Fast | 401점 ±550 MHz, 적응 pole–residue |
| 본 표 solver | Balanced와 동일 pole–residue, 모든 10점 직접 solve | 출력/TPD 보간 없음; 현행 보정 유지 |
| Floquet | N_F=3, 인접 N_F=2 감사 | 추가 N_F=4 교차 검증 |
| Doppler | 1D, Δv=1 m/s, ±4σ | 118 °C에서 1567 velocity classes |

Default 자체는 온도121 °C·Pump600 mW·OPD0.9 GHz·Seed8 µW. 위 네 사용자 지정값으로만 대체. 저장된 `phase_detail=Balanced` 값 대신 `resolution` tier가 실제 전파 `ultra`를 지정함. bare `compute_spectrum()`의 module defaults는 UI와 다르므로 본 스크립트는 `FWMScheme().defaults()`를 읽고 물리값을 명시 전달.

Git branch **main**, HEAD **`{manifest['head']}`**. 실행 당시 다수 staged/unstaged·untracked 변경 존재. 결과는 깨끗한 HEAD가 아니라 **현재 작업 트리**의 계산. 기존 파일 678개 SHA-256, 초기 Git 상태, 전체 defaults/사용 params는 {link('실행 manifest','analysis/tpd_gain_diagnostic_20261001/source_manifest_before.json')}에 보존. Python 3.14.2, NumPy 2.4.1, SciPy 1.17.0. 본 캠페인은 NumPy pole backend와 guard를 사용했고 compiled Numba 경로로 같은 10점 추가 검증.

조건·변환 근거: {link('UI schema','gabes/schemes/fwm.py',3638)}, {link('실행 mapping','gabes/schemes/fwm.py',4029)}, {link('beam intensity/radius','gabes/constants.py',92)}, {link('density','gabes/hyperfine.py',78)}.

## 3. Gain 정의·주파수 대응

실측은 같은 검출기 앞에서 **pump-on Probe / pump-off Seed**. 표의 Seed **3.7/3.8 µW를 그대로 분모**로 사용해 반올림 이전 Gain 재계산. 모델 입력 Seed는 사용자 지정 **3.7 µW**로 공통 유지. 3.8 µW 분모를 3.7로 바꾸거나 지점별 입력 온도·배율을 적합하지 않음.

모델 `G_s`: 셀 출구 Probe power / 셀 입력 Seed power. `G_c`: 셀 출구 Conjugate power / 같은 셀 입력 Seed. 작은 신호 canonical gain과 pump-budget/saturation 적용 후 출력 Gain 구분; 본 표는 실제 표시용 `G_s`,`G_c` 사용. pump-off transmission으로 자동 나누지 않음. 같은 downstream Probe 경로·입력 모드/파워가 유지되면

\[
G_{{det}}=G_{{source}}/T_{{off}}.
\]

경로가 달라지면 추가 collection/path 비가 필요. 총 η를 이 식에 곱해 Gain을 맞출 근거 없음. 따라서 완전히 같은 측정면으로 변환한 비교는 현재 미완료이며, 아래 잔차는 **제공된 pump-on/off 비와 source-gain의 조건부 비교**.

별도 4.7→4.6 µW 구간의 \(T=0.978723\)를 전 점에 가정해도 source-equivalent 실측은 2.13% 감소, 마지막 Gain **116.124**. 관측 불일치를 없애지 못함. 이 값은 각 EOM의 전체 셀/검출 경로 transmission이나 pump-on 흡수, 양팔 η 측정이 아님. 같은 경로 공식만으로 현행 Gain을 맞추려면 T_off가 약 **47.3%→9.50%**여야 하나, 이를 지지하는 측정 없음. 입력3.7과 검출기 앞3.7/3.8도 별도 측정 불확실성 없이 정밀 transmission 교정으로 사용 불가.

기준 hyperfine **νHF=3.035732439 GHz**. red sideband에서

\[
\delta_{{GABES}}/(2\pi)\,[\mathrm{{MHz}}]=(3.035732439-f_{{EOM}}[\mathrm{{GHz}}])\,1000.
\]

따라서 3.039→3.030 GHz는 δ **−3.267561→+5.732439 MHz**. 실험 기록의 반올림된 “3.039=+3 MHz” 표시는 GABES와 부호·원점 차이 존재. 본 비교는 EOM 절대 주파수를 사용해 정합. OPD를 바꾸어도 동일 EOM의 δ는 동일.

근거: {link('hyperfine constant','gabes/constants.py',26)}, {link('branch center/TPD','gabes/schemes/fwm.py',2479)}, {link('EOM red-sideband mapping','sabes/beamline.py',96)}, {link('gain convention','gabes/schemes/fwm.py',4166)}.

## 4. 실측 10점·잔차

잔차 \(r=G_{{meas}}-G_{{sim}}\). 상대오차는 \(100(G_{{sim}}-G_{{meas}})/G_{{meas}}\). 두 부호 정의를 구분.

| EOM GHz | δ MHz | 실측 Gain | 보정 ON G_s | 보정 OFF G_s | r | 상대오차 % | 실측/ON |
|---:|---:|---:|---:|---:|---:|---:|---:|""")
    for row in table:
        sections.append("| " + " | ".join(f"{float(row[k]):.{digits}f}" for k,digits in
            (("EOM_GHz",3),("TPD_GABES_MHz",6),("gain_measured",4),
             ("gain_exact_corrected",4),("gain_exact_uncorrected",4),
             ("residual_measured_minus_sim",4),("simulation_minus_measured_pct",2),
             ("measured_over_sim",4)))+" |")
    sections.append(fr"""
Probe/Conjugate 원자료, 현행 Conjugate Gain, UI Fast/Balanced 표시값 포함 전 정밀도 자료: {link('comparison.csv','analysis/tpd_gain_diagnostic_20261001/comparison.csv')}.

Conjugate 예측 5.477→10.524, 실측 Conjugate/Seed 12.838→121.622. Probe와 Conjugate의 수기 동일 기재는 무한 정밀도 균형 측정이 아님. 이를 lossless `G_s−G_c=1` 또는 양팔 calibration 증거로 사용하지 않음.

{figure('gain_comparison')}

실측 총 증가 **105.216**, 모델 **4.924**. 인접 1 MHz 이동의 실측 증가량은 1.757,9.405,7.248,15.455,10.811,17.568,17.219,15.213,10.541; 모델은 0.454–0.609. 실측 후반 증가세 둔화는 보이지만 끝 이후 감소 데이터 없음. \(G/G_{{3.030}}\) 정규화는 배율만 제거하며 offset 효과는 남김. offset까지 제거한 비교는 아래 affine 진단 그래프에 따로 표시.

## 5. 단일 배율·Affine

10점 전체 **동일 가중 최소제곱**. 측정 표준편차가 없어 inverse-variance 가중, 오차막대, confidence interval, p-value 생성하지 않음.

| 입력 곡선 | 적합 | a | b | RMSE | R² | 최대 상대 잔차 % |
|---|---|---:|---:|---:|---:|---:|""")
    for label, fits in (("현행 ON",c),("OFF",u)):
        for name in ("identity","scale","affine"):
            v=fits[name]
            sections.append(f"| {label} | {name} | {v['a']:.6f} | {v['b']:.6f} | {v['rmse']:.4f} | {v['r2']:.6f} | {v['max_abs_relative_pct']:.2f} |")
    sections.append(fr"""
단일 배율 ON \(a={c['scale']['a']:.6f}\): 앞부분 과대, 뒤쪽 과소; 잔차 **−32.69,−34.23,…,+30.65,+36.77**. 한 공통 배율로 상승폭을 설명하지 못함.

Affine ON: 잔차 **+8.757,+0.068,−1.610,−6.071,−2.919,−4.959,−0.718,+2.797,+4.065,+0.590**. 작은 절대 RMSE와 동시에 저점 65.19% 오차·중앙 음/양끝 양의 구조 존재. 이 잔차가 실제 측정 분산보다 큰지는 미확정. R²는 통계적 검증·외삽 타당성의 증거가 아님.

\(23.027G-141.568\)은 G≈6.148 이하에서 음수; G=1이면 **−118.54**. 따라서 이번 좁은 상승 구간의 진단 근사로만 유효. 최종 Gain에 바로 적용할 물리 보정으로 제안하지 않음. 일반 affine는 양의 a일 때 모델 피크 위치를 옮기지도 못함. Gain 곡선 자체 비선형성과 두 Gain 사이 affine 관계의 가능성을 구분.

{figure('affine_diagnostic')}

## 6. 공통 온도·피크

온도 센서는 cell body, 과거 voltage→temperature 관계 사용. OD 기반 density 확인 없음. 당일 정밀도·온도 시간 로그 없음. **±1,±3,±5 °C는 단계적 민감도 범위**이며 장치 보장 정밀도·통계 신뢰구간 아님. UI 입력 범위 60–150 °C 안에서 계산했으나 UI 범위가 모델의 실험적 유효 범위를 보장하지 않음. 온도는 density, collision, Maxwell velocity weights를 함께 바꿈.

113–123 °C를 **0.5 °C 간격**으로 조사. 각 곡선의 온도는 전체 EOM에 공통. 주요 점:

| 공통 T °C | G(3.039) | G(3.030) | 마지막/첫째 | 원값 RMSE | 피크 EOM GHz | 피크 Gain |
|---:|---:|---:|---:|---:|---:|---:|""")
    for t in (113,115,117,118,119,121,123):
        tag=f"temperature_{t:.2f}"
        g=data[tag]["G_s"]
        pk=summary["peaks"][f"peak_{t:.2f}"]
        sections.append(f"| {t} | {g[0]:.4f} | {g[-1]:.4f} | {g[-1]/g[0]:.4f} | {summary['fits'][tag]['identity']['rmse']:.4f} | {pk['eom_GHz']:.5f} | {pk['gain']:.4f} |")
    best=extra["exploratory_fits"]["exploratory_best_temperature"]["identity"]
    sections.append(fr"""
温度 +1/+3/+5 °C에서도 마지막 Gain **13.923/21.689/34.654**. ±5 °C에서 원값 RMSE 최소는 상한123 °C **44.498**, 아직 충분한 설명 아님. 공통 온도에 추가 공통 출력 배율을 허용해도 최저 RMSE 약25.7이며 증가 형태 불일치 유지. Affine가 추가되면 온도별 a,b가 크게 바뀌어 온도·보정의 원인을 식별할 수 없음.

추가 탐색 **113–132 °C**는 허용 오차 확대가 아닌 원인 분리용. 고정 현행 보정에서 전체 RMSE 최소 온도 **{extra['exploratory_best_common_temperature_C']:.3f} °C**, RMSE **{best['rmse']:.3f}**, 첫/끝 Gain **46.244/71.368**. 마지막 점만 일치시키는 **{extra['temperature_matching_last_point_C']:.3f} °C**에서는 첫 Gain **85.467**, 피크 약 **3.0310 GHz/G119.668**. 크기와 피크를 가까이 가져오면서 앞쪽 형태를 잃음. 이 탐색을 실제 셀 온도 추정으로 제시하지 않음.

{figure('temperature_sensitivity')}

피크 스캔: **EOM 2.990–3.050 GHz**, 0.25 MHz 간격, 모든 점 직접 solve. 118 °C ON **3.0220 GHz/G14.576**, OFF **3.0105 GHz/G1216.881**. 피크 값은 이 구간 grid 최대이며 0.25 MHz는 계산 간격이지 실측 정밀도·오차막대가 아님. 구간 밖 모든 extrema의 전역 탐색은 아님. OPD0.9 GHz 대조에서도 현행 피크3.0210/G15.219로 나타남. 따라서 OPD 후속 정정만으로 피크·규모 정합 회복 안 됨.

{figure('peak_sensitivity')}

실측 표의 maximum은 endpoint118.649. 사용자의 **3.030 부근 turnover 관측은 존중하되, 감소 구간 정량 검증과 분리**. 노션·두 채팅·원자료 폴더에서 추가 광파워 감소점 찾지 못함.

## 7. Semi-empirical 보정 감사

현식:

\[
C_{{mix}}=0.5594938027,\quad
(\bar\chi_{{ss}},\bar\chi_{{cs}},\bar\chi_{{sc}},\bar\chi_{{cc}})
\mapsto(\bar\chi_{{ss}},C_{{mix}}\bar\chi_{{cs}},C_{{mix}}\bar\chi_{{sc}},\bar\chi_{{cc}}).
\]

유한 Seed atomic solve → Doppler average → 두 off-diagonal scaling → Maxwell matrix → 64-segment propagation/pump budget → final saturation. diagonal absorption/dispersion 및 atomic state 변화 없음. `line_strength=0.74`는 별도 macroscopic residual. 보정은 **내부 계수의 일정 배율**이며 \(G_{{ON}}=C_{{mix}}G_{{OFF}}\)가 아님. 이번 ON/OFF 비도 첫점0.07844→끝점0.05137로 변함.

{link('상수·provenance','gabes/fwm_gain_closure.py',12)}, {link('apply 함수','gabes/fwm_gain_closure.py',104)}, {link('적용 위치','gabes/schemes/fwm.py',3141)}, {link('Maxwell matrix','gabes/observables.py',48)}. Fast/Balanced에서만 eligible. Ultra tier는 보정 OFF, 따라서 ON Balanced 대 OFF Ultra 차이를 순수 solver 오차로 해석하면 안 됨.

Calibration: Sim2025, OPD0.9 GHz, 코드 δ−8 MHz, T121 °C, Pump600 mW, Seed8 µW, 동일 L/waist/angle. 목표 Gain **15.5**는 문헌 약15–16의 대표값. 별도 반올림 광파워111/8=13.875와 구분. 기존 uncorrected Gold394.299→corrected15.5. η/Noise를 fit에 사용하지 않음. 현재 조건은 OPD+0.16 GHz, T−3 °C, Pump−220 mW, Seed−4.3 µW. 단일 기준점 외삽 정확도 미검증; 저장 held-out Liu6.445 vs8, McCormick1.926 vs9도 혼합 결과.

Gaussian projected pump participation **0.720626**은 코드에 기록된 **미적용 진단값**. 실제 계수0.559494의 온도·공간·Zeeman·편광 기여 분해는 미확정. 이 진단값으로 현 계수를 대체하지 않음.

### 내부 공통 계수의 진단 민감도

기존 상수·함수 수정 없이 untouched χ에 진단 인자 C를 넣어 기존 propagation 함수 호출. C현행과 C=1 결과가 각각 production ON/OFF에 수치적으로 일치하는지 확인. 이는 가상 민감도 평가이며 기본 계수 교체·변경안 적용 없음.

118 °C에서 전체10점 최소제곱을 줄이는 C≈**0.879072**, RMSE **19.588**, 첫/끝 **39.678→96.579**, 증가비 **2.434배**. 피크 약 **3.014 GHz/G288**로 실측 관측에서 더 멀어짐. 113/123 °C의 같은 진단에서 계수≈1.199/0.658로 달라짐. 온도와 계수 분리 불가; 113 °C의 C>1은 attenuation형 participation 해석과도 맞지 않음. 이 결과는 새 계수 채택 근거 아님.

이 C를 전역 교체한다고 가정하면 frozen Gold Gain **15.5→162.017**, 약 **10.45배**. 현재 기준점 재현을 크게 훼손. 조건마다 다른 계수를 넣어 해결하면 predictivity를 잃고 과적합 위험 증가. 실제 교체 없음.

진단 재현·기준점 영향: {link('mixing sensitivity','analysis/tpd_gain_diagnostic_20261001/mixing_sensitivity.json')}. 전체 온도/계수 함수 공간에 대한 불가능성 증명은 수행하지 않음. 시험한 단일계수와 공통온도 가설에서 남는 구조적 잔차를 보고.

### 기준점 TPD 부호의 미확정 사항

[Sim2025 원문](https://www.nature.com/articles/s41598-025-86479-w)은 red seed와 δ−8 MHz를 명시. [Fig.1(b)](https://www.nature.com/articles/s41598-025-86479-w/figures/1)의 에너지 그림은 paper δ=(pump−seed)−HF를 시사하지만, 명시적 signed equation·Gold 실제 EOM 없음. **그림에서 추론한 부호 차이 후보이며 확정 bug 아님**. 기존 ledger는 paper−8을 GABES−8로 입력. 코드 기준점 EOM **3.043732439 GHz**; 반대 부호 해석이면 **3.027732439 GHz/GABES+8 MHz**. 재보정 전에 실제 기준 EOM 또는 원저자의 부호 정의 확인 필요.

## 8. 계산 모드·수렴·η/Noise

| 교차 검증 | 10점 G_s 최대 상대 차이 |
|---|---:|
| ON pole vs ON direct-grid control | 2.46×10⁻¹⁰ % |
| OFF pole vs OFF Ultra direct grid | 2.61×10⁻¹⁰ % |
| N_F=3 vs4 | 1.40×10⁻¹⁰ % |
| NumPy vs compiled Numba pole/guard | 5.55×10⁻¹¹ % |
| Δv1→0.5 m/s, ±4σ | 0.00001779 % |
| Δv0.5 m/s, cutoff±5σ vs baseline | 0.004227 % |
| UI Default Fast readout vs exact rows | 0.3485 % |
| UI Balanced readout vs exact rows | 0.3499 % |
| exact rows η/内部 noise change | 0 % |
| input Seed3.7→3.8 µW control | 0.001728 % |
| input Seed3.7→8 µW control | 0.07426 % |

全10点 Floquet3→2/4→3 **CONVERGED**、pole guard max relative约1.22×10⁻¹²、fallback0。小於實測差異。401점 focused UI의 **2.75 MHz grid + `np.interp`**는 1 MHz EOM 자료보다 거침. 직접 점 계산을 본 결론에 사용. Fast의 spline 보간과 최종 UI linear interpolation은 별개.

**수치 구현 문제 후보 확인:** η는 정확한 physical Gain에 들어가지 않지만 Fast의 refinement score에는 noise를 통해 입력. η86.94→50%에서 solved rows **69→63**, 실측10점 readout 최대 차이는 **0.00002938%**로 작음. 그러나 확장 full401 curve의 **EOM3.010982439 GHz**에서 Gain **2.442647→1.182485**, **51.59%** 변화. Balanced는 η 변경에 수치 roundoff 이내 불변. 두 모드의 Floquet convergence badge는 이 adaptive artifact를 잡지 못함.

따라서 η/Noise로 Gain을 맞추지 않음. 확장 피크·온도 분석은 모든 점 solve한 Balanced/direct rows 사용. 수치 개선 후보는 η-independent Gain 검증과 적응 spline의 직접점 비교·fallback. 실제 변경 미실행.

TPD·사용자 noise는 readout-only. Temp/Pump/OPD/Seed/η/closure는 compute cache key; 직접 스크립트는 Streamlit cache 사용 안 함. 소스 변경 전 생성된 실행 session/cache를 실제 default 재현 증거로 삼지 않음. {link('cache','streamlit_app.py',82)}, {link('Fast scoring','gabes/schemes/fwm.py',3055)}, {link('operating point interpolation','gabes/schemes/fwm.py',3466)}.

## 9. 판단별 다음 단계

| 후보 | 현재 판단 | 변경 전 검증 |
|---|---|---|
| 현 보정 유지+온도만 변경 | ±5 °C에서 곡선 설명 부족; 확장 탐색도 전체 곡선 실패 | OD density/실효 온도, 시간 안정성 |
| 한 출력 배율 | 전체10점 잔차 큼 | 반복측정 이후에도 residual trend 유지 여부 |
| 한 내부 mixing 계수 | 진단 최적값도 형태·피크 실패 | 부호 확정, 독립 기준·held-out 조건 |
| Affine 출력 보정 | 좁은 구간 근사 양호, 큰 음 offset·피크 이동 불가 | 측정 background/zero 검증; 전역 적용 금지 |
| 현 함수 형태/적용 범위 | 외삽 정확도 부족 후보; 새 함수의 필요성 확정 못함 | 물리 participation/geometry/temperature와 구분 가능한 다조건 데이터 |
| Gain 면·TPD 단위/부호 | 현10점 EOM 정합 확인; plane transmission·Gold sign 미확정 | matched pump-off paths, Gold 실제 EOM |
| 수치 계산 | 본10점 수렴 양호; 확장 Fast adaptive artifact 확인 | direct Balanced 대조, η 바꾸어도 same G 검증 |
| 부족한 물리 | 1D Doppler·4준위·유한 reference·unmeasured conjugate mode 후보 | density/beam/angle/polarization 교정 후 별도 reference 비교 |

0.32°·118 °C 조건의 별도 analytic Raman-Doppler budget은 transverse RMS 약 **1.375 MHz**. 생산 응답은 transverse Raman velocity average를 포함하지 않음. 이는 형태 민감도 후보이며, 8 MHz peak 차이의 원인으로 단정하지 않음. full Zeeman은 CG-sum bookkeeping으로만 존재. microscopic quantum noise 부재가 이번 mean optical gain 잔차의 직접 원인이라고 단정하지도 않음.

수정한다면 검토 위치: `gabes/schemes/fwm.py` adaptive response 검증·실험 Gain 기준면 비교 인터페이스; 기준점 부호 검증 후 `gabes/fwm_gain_closure.py` calibration provenance/적용 범위. **고차 다항식·점별 보정·η/Noise fit 제안 없음.** 특정 TPD에서만 맞춘 계수는 과적합·Gold 훼손·다른 Pump/온도 외삽 위험. 승인 이후에도 frozen Gold+기존 held-out+별도 실측 구간으로 검증하고, amplitude·shape·peak 세 지표를 함께 평가해야 함.

## 10. 최소 추가 측정·자료 한계

1. **짧은 randomized/interleaved scan:** EOM3.039/3.035/3.030을 바꾸되 매번 한 기준 EOM으로 복귀. 시각, 표시 온도, Pump/OPD, pump-off Seed와 pump-on Probe/Conjugate 함께 기록. 순차 drift와 TPD dependence 분리.
2. **피크 양쪽:** 3.032/3.031/3.030/3.029/3.028 GHz 최소 점을 반복, 필요시 범위 확장. 감소 branch·peak location 확인. 같은3.030에서 수 분 반복해 시간 변화도 확보.
3. **OD 또는 독립 density/실효 온도:** 별도 heater 추정과 분리; beam path와 cold spot 관계 점검. density와 Doppler-width 변화에 대한 독립 근거 확보.
4. **맞춘 기준면:** 각 EOM에서 입력/검출기 앞 pump-off Seed, pump-on background·모드/collection 기록. 사용된3.7 µW input과 detector 기준 분모 관계 확인. 4.7→4.6 한 기록을 전체 transmission으로 확대하지 않음.
5. **보정 재교정 전 Gold 실제 EOM·부호 확인**, 이후 최소 한 독립 Pump 또는 온도 곡선을 held-out으로 사용. 변경 함수 평가용 데이터와 coefficient fit용 데이터 분리.

SA CSV는 frequency/dBm 잡음 스펙트럼으로 분류, Gain 추출 미사용. Scope seconds/Volts는 SAS OPD 교정 참고; 미교정 channel을 광파워 독립 검증으로 미사용. +1.06 GHz는 저장 SAS output maximum과 이번 사용자 지정값을 구분해 기록. raw 폴더의 과거344.8 mW·9/09/9/10 조건 합치지 않음.

노션 두 페이지·참고 채팅 둘 확인. 추가 감소점·반복 분산·당일 온도/시간 로그 없음. pump와 OPD가 drift 중 안정적이었다는 사용자 기록 있으나 온도는 당시 확인하지 않았음. 따라서 drift를 온도 원인으로 확정 못함. 오차막대·신뢰구간·통계 유의성 미생성.

참고: [9/30 측정 기록](https://www.notion.so/3eb6cba14fee810fa8cbed64f094f3b9), [9/28 분석 기록](https://www.notion.so/3e96cba14fee817bb0ddf5ed5ef3057c), [TPD Gain 정량 확인](https://chatgpt.com/c/6abcb00b-7e7c-83ee-b84e-64ac3c05dc50), [실험 기록 분석 제안](https://chatgpt.com/c/6abb4445-2ac8-83ee-85ab-5672f9828e4b). 과거 AI의 수치·해석은 새로운 측정 증거로 사용하지 않음.

## 11. 재현·보존·검증 기록

신규 출력은 이 `analysis/tpd_gain_diagnostic_20261001/` 폴더만 사용. diagnostic callbacks의 입력을 바꾼 가상 계산과 production coefficient 수정은 구분. 기존 상수·기준 파일 수정, commit, stage, branch 변경 없음.

- {link('기본 진단 스크립트','analysis/tpd_gain_diagnostic_20261001/run_diagnostic.py')}: 직접10점, UI 모드, 수렴, ±5 °C, 밀집 피크·그래프.
- {link('추가 진단','analysis/tpd_gain_diagnostic_20261001/additional_checks.py')}: 확대 온도 가설, compiled 대조, η/Noise/Seed controls.
- {link('mixing 진단','analysis/tpd_gain_diagnostic_20261001/mixing_sensitivity.py')}: 상수 수정 없는 C 민감도·Gold 영향.
- {link('전체 계산 결과','analysis/tpd_gain_diagnostic_20261001/results.json')}, {link('fit/수렴 summary','analysis/tpd_gain_diagnostic_20261001/summary.json')}, {link('추가 수치 감사','analysis/tpd_gain_diagnostic_20261001/additional_checks.json')}, {link('온도 fit CSV','analysis/tpd_gain_diagnostic_20261001/temperature_fits.csv')}.

`python -m pytest -q` 전체 실행: **2077 passed, 3 skipped**, 877.17초. 캐시는 신규 진단 폴더로 지정. 기본·추가·mixing 진단 스크립트 실행 완료, 네 그래프 시각 확인 완료. 계산 전후 기존 파일 **678개 SHA-256 동일**, HEAD/main 유지.

검증 상태는 {link('verification.json','analysis/tpd_gain_diagnostic_20261001/verification.json')}에 최종 기록. 모델 수치 수렴·repository tests와 실제 실험 agreement는 서로 다른 검증 범위.
""")
    text="\n".join(sections)
    # Korean report; preserve math/identifiers and ASCII-only plot axes.
    text=text.replace("温度", "온도").replace(
        "全10点 Floquet3→2/4→3 **CONVERGED**、pole guard max relative约1.22×10⁻¹²、fallback0。小於實測差異。",
        "전체10점 Floquet3→2/4→3 **CONVERGED**, pole guard 최대 상대차 약1.22×10⁻¹², fallback0. 실측 차이보다 작음. ")
    (OUT/"REPORT_ko.md").write_text(text,encoding="utf-8")
    print("Report saved", OUT/"REPORT_ko.md")


if __name__=="__main__":
    main()
