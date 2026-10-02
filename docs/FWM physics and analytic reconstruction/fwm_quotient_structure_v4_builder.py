"""Build the self-contained FWM atlas; standard library only.

Edit the node/edge data here, then run this file. The delivered HTML needs no
builder, companion assets, network, or JavaScript to read its full content.
"""
from pathlib import Path
from html import escape
import json

ROOT = Path(__file__).resolve().parent

PRESERVED_CONTEXT = r'''  <section class="sheet evidence" id="evidence" aria-labelledby="evidence-title">
    <div class="evidence-header"><h2 id="evidence-title">주장마다 다른 검증 근거</h2><small>선언된 참조 fixture의 보고값 · 이 페이지에서 재계산한 결과가 아닙니다.</small></div>
    <p>방정식 잔차, 구현 간 차이, 극한 수렴, 실험 비교를 별도로 읽습니다. 작은 잔차가 작은 물리적 예측 오차를 뜻하지는 않습니다.</p>
    <div class="table-scroll"><table>
      <caption class="sr-only">검증 대상, 보고된 결과, 성립 범위와 근거</caption>
      <thead><tr><th scope="col">검증 대상</th><th scope="col">보고된 결과</th><th scope="col">이 결과가 말하는 범위</th></tr></thead>
      <tbody>
        <tr><td>참조 응답의 방정식 잔차</td><td><b>7.32 × 10⁻¹⁷</b><p>정규화 최대 잔차</p></td><td>참조 응답이 자기 방정식을 만족하는 정도. 생산 경로와의 차이가 아닙니다.<p>Pump-only reference diagnostics · <code>tab:pump-reference-diagnostics</code></p></td></tr>
        <tr><td>pump 상태의 프레임 비교</td><td><b>8.37 × 10⁻¹⁵</b><p>정적 프레임 ↔ N<sub>F</sub>=3</p></td><td>minus의 pump-only 상태를 게이지 변환해 Floquet 상태와 비교. 유한 시드 전체 응답의 동등성을 뜻하지 않습니다.</td></tr>
        <tr><td>유한 시드 → 무한소 응답</td><td><b>2.503 × 10⁻³</b><p>시드 8 µW · N<sub>F</sub>=3</p></td><td>속도 평균한 네 복소 응답 성분을 참조 최대 응답으로 정규화한 최악 차이. 시드 감소에 따른 비교는 아래 표에 있습니다.</td></tr>
        <tr><td>Floquet 절단 수렴</td><td><span class="status limited">스캔마다 확인</span><p>N<sub>F</sub>=3/2 비교 경로 구현</p></td><td>현재 파라미터에서의 합격 여부는 실제 실행 결과로 판정합니다. 과거 N<sub>F</sub>=1의 16.8% 변화는 모든 차수의 실패를 뜻하지 않습니다.</td></tr>
        <tr><td>2D Raman–Doppler</td><td><span class="status limited">참조 범위 내 검증</span><p>rms ≈ 1.380 MHz</p></td><td>121 °C, 0.32°, minus 참조의 해석 rms·격자·cutoff 대조. 생산 스캔의 2D 승격이나 잡음 스펙트럼 검증이 아닙니다.</td></tr>
        <tr><td>실험 이득 비교</td><td><span class="status fail">큰 불일치</span></td><td>문헌 운전점의 절대 이득을 현재 축약 모델이 재현하지 못합니다. 수치 일관성과 실험 일치는 다른 주장입니다.</td></tr>
        <tr><td>물리적 스퀴징 스펙트럼</td><td><span class="status need">미완성</span></td><td>미시 M<sub>q</sub>(Ω), D(Ω), 입력 상태 및 검출 응답의 완결된 연결이 필요합니다. 이득 기반 숫자로 대체할 수 없습니다.</td></tr>
      </tbody>
    </table></div>
    <details><summary>시드 세기를 줄이면 참조 응답에 어떻게 접근하는가?</summary><div class="disclosure-content">
      <div class="table-scroll"><table><caption class="sr-only">원문 유한 시드 극한 비교 표</caption><thead><tr><th scope="col">Rabi 진폭 비율</th><th scope="col">시드 파워</th><th scope="col">최악 정규화 차이</th></tr></thead><tbody><tr><td>1.00</td><td>8.0000 µW</td><td>2.503 × 10⁻³</td></tr><tr><td>0.30</td><td>0.7200 µW</td><td>5.093 × 10⁻⁴</td></tr><tr><td>0.10</td><td>0.0800 µW</td><td>9.254 × 10⁻⁵</td></tr><tr><td>0.03</td><td>0.0072 µW</td><td>9.138 × 10⁻⁶</td></tr></tbody></table></div>
      <p>같은 1D 속도 평균과 다섯 Raman detuning에서 N<sub>F</sub>=3 응답을 비교한 보고값입니다. 마지막 세 점의 진폭에 대한 log–log 기울기는 1.749입니다. 이 유한 범위의 관찰로 정확한 2차 수렴 법칙을 주장하지 않습니다.</p>
      <p>응답 잔차 표의 별도 검사 범위: Ω<sub>SA</sub>/2π = 0, 0.1, 1, 4 MHz; δ/2π = −30, −8, 0, 12, 40 MHz. 서로 다른 오차량을 하나의 ‘인증 진폭’으로 합치지 않습니다.</p>
      <p><a href="squeezing_analytic_reconstruction_v3.tex">원문 TeX</a> · 검색 키: <code>tab:pump-reference-diagnostics</code>, <code>tab:finite-seed-reference-limit</code></p>
    </div></details>
  </section>

  <details class="math-details" id="math-notes">
    <summary>수학적 구조 더 읽기 · 합성, 조건, 순서와 잡음</summary>
    <div class="math-grid">
      <article><h3>01 · 검출은 두 경로의 합류</h3><p>이론적으로 평균장을 예측하면 D → F가 필요합니다. D(평균장), X(이득 진단), E(공분산), F(검출)에 D &lt; X, D &lt; F, E &lt; F만 있는 유도 부분순서는 N형입니다. 실제 의존관계를 보존하며 SP를 강제하지 않습니다.</p><p>평균장 진폭·위상을 실험 입력으로 대체하면 다른 조건부 모델이 됩니다. 그 입력은 숨기지 않고 검출 포트에 표시해야 합니다.</p></article>
      <article><h3>02 · Gaussian channel의 합성</h3><div class="equation-box formula">V ↦ XVXᵀ + Y<br>(X₂,Y₂) ∘ (X₁,Y₁)<br>= (X₂X₁, X₂Y₁X₂ᵀ + Y₂)</div><p>앞 구간의 잡음은 뒤 구간에서 다시 변환됩니다. 독립적인 구간 저장고를 가정하며, 구간 간 상관이 있으면 추가 항이 필요합니다. 여기의 X는 위 도식의 이득 진단 식별자 X와 다른 채널 행렬입니다.</p></article>
      <article><h3>03 · 이득 검사와 정준성은 다름</h3><div class="equation-box formula">𝒢<sub>flux</sub> = |T₁₁|² − |T₂₁|²<br>TJT† = J, &nbsp; J = diag(1,−1)</div><p>𝒢<sub>flux</sub>=1만으로 일반 전달행렬의 정준성이 보장되지 않습니다. T=diag(1,2)는 이 스칼라 검사를 통과하지만 TJT†=diag(1,−4)입니다. 이상 IDS는 입력 상태와 검출 조건까지 요구합니다.</p></article>
      <article><h3>04 · 허용 가능한 잡음은 하나가 아님</h3><div class="equation-box formula">Y + (i/2)(Ω<sub>R</sub> − XΩ<sub>R</sub>Xᵀ) ⪰ 0</div><p>정규화된 실수 quadrature 모드, V<sub>vac</sub>=I/2의 Gaussian channel 완전양성 조건입니다. Ω<sub>R</sub>은 symplectic 행렬이며 분석 주파수가 아닙니다. 주파수별 적용에는 companion sideband를 포함한 모드 조립이 선행합니다.</p><p>물리적으로 허용되는 채널이라고 실제 원자 잡음 모델인 것은 아닙니다. 고정된 전달·입력·검출 조건에서 W₂−W₁⪰0이면 두 번째 모델의 측정 잡음이 더 크지만, 일반적인 두 W는 부분순서상 비교 불가능할 수 있습니다.</p></article>
      <article><h3>05 · 평균과 전파는 교환 불가</h3><div class="equation-box formula">exp(L⟨M(v)⟩) ≠ ⟨exp(LM(v))⟩</div><p>일반적으로 성립하는 비동등성입니다. 같은 장을 구동하는 속도군의 응답을 평균하는 것과, 서로 다른 전달계를 평균하는 것은 다른 모델입니다. 독립 속도군의 확산은 잡음 진폭을 평균한 뒤 제곱하지 않고 공분산 기여로 더합니다.</p></article>
      <article><h3>06 · 규약과 근거의 역할</h3><p>규약을 공통 환경으로 모으는 것은 동치관계에 의한 수학적 몫과 다릅니다. 이 도식은 입력·출력 의존관계와 조건을 명시하며, 검증 연결을 정의되지 않은 자연변환이나 2-cell로 부르지 않습니다.</p><p>검증은 잔차·구현 비교·극한·실험으로 나눕니다. 수학적 정식화를 더한다면 객체·사상·합성과 오차의 합성 규칙을 별도로 정의해야 합니다.</p></article>
    </div>
    <div class="detail-links"><a href="https://arxiv.org/abs/1510.05468" target="_blank" rel="noopener noreferrer">과정과 string diagram</a><a href="https://arxiv.org/abs/1110.3234" target="_blank" rel="noopener noreferrer">Gaussian quantum information</a><a href="https://arxiv.org/abs/0909.0408" target="_blank" rel="noopener noreferrer">Gaussian channel 합성</a><a href="https://epubs.siam.org/doi/10.1137/0211023" target="_blank" rel="noopener noreferrer">SP 순서 구조</a></div>
  </details>
'''

NUMBERS = {'input':'01','atomic':'02','reference':'R','mean':'03','gain':'04',
           'micro':'05','covariance':'06','detector':'07','output':'08'}
CODE = {
 'input':('gabes/constants.py','rabi_freq; 1/e² intensity radius'),
 'atomic':('gabes/schemes/fwm.py','chi_matrix_table · compute_spectrum · assess_floquet_scan_convergence'),
 'reference':('gabes/schemes/fwm.py','pump_only_weak_response_reference · pump_only_weak_response_noncollinear_reference'),
 'mean':('gabes/observables.py','gain_from_chi · canonical_transfer_diagnostics'),
 'gain':('gabes/observables.py','canonical_transfer_diagnostics; ideal-limit comparisons'),
 'micro':('gabes/schemes/fwm.py','seeded_validation_claim_gate; microscopic diffusion unavailable'),
 'covariance':('analysis/squeezing/analytic_reconstruction/final_run.py','analytic identity fixtures; not a production quantum solver'),
 'detector':('sabes/detection.py','detector boundary; the microscopic FWM covariance is not supplied'),
 'output':('gabes/schemes/fwm.py','seeded_validation_claim_gate; mean-field Squeezing indicator only')}

GROUPS, EDGES = [], []
KINDS = {'ready': '생산 구현', 'reference': '별도 참조', 'analytic': '가정하의 해석식',
         'need': '전체 연결 미완료', 'test': '대수·수치 시험', 'input': '물리 입력'}
EDGE_KINDS = {'data': '입력 전달', 'alternative': '동일 문제의 대안',
              'limit': '가정하의 축약', 'validation': '검증 대상', 'analogy': '형태의 유사성'}

def group(id, title, subtitle, kind, part, brief, evidence):
    GROUPS.append(dict(id=id, title=title, subtitle=subtitle, kind=kind, part=part,
                       brief=brief, evidence=evidence, nodes=[]))

def node(id, title, out, xy, math, inputs, units, body, condition, refs,
         kind=None, evidence=None, compact=None):
    g=GROUPS[-1]
    g['nodes'].append(dict(id=id, group=g['id'], title=title, out=out, xy=xy,
        math=math, compact=compact or math.split('\n')[0], inputs=inputs, units=units,
        body=body, condition=condition, refs=refs.split(), kind=kind or g['kind'],
        evidence=evidence or g['evidence']))

def edge(source, target, label, kind='data', pending=False):
    EDGES.append(dict(source=source,target=target,label=label,kind=kind,pending=pending))

group('input','운전점과 입력','구동 · 원자 · 수집 모드','input','I · Definitions and conventions',
      '같은 입력도 쓰이는 위치가 다릅니다. 광학 구동, 원자 상태, 전파 형상, 입력 요동, 검출 보정을 각각 전달합니다.',
      '입력의 출처와 모델 가정입니다. 이 화면은 장치 파라미터를 측정하거나 실시간 계산하지 않습니다.')
node('i-levels','네 준위와 double-Λ','구동·읽기 전이의 정의',(1,1),
     'g₁=5S₁/₂(F=2), g₂=5S₁/₂(F=3)\ne₂=5P₁/₂(F′=2), e₃=5P₁/₂(F′=3)',
     '⁸⁵Rb D1 원자종·초미세 준위·전이 강도','4준위 ρ: 4×4; vectorized ρ: 16성분',
     'minus 모드쌍에서 probe는 g₂ 쪽, conjugate는 g₁ 쪽 optical coherence를 읽습니다. 펌프는 두 ground manifold를 구동하며 두 excited manifold의 진폭이 함께 기여합니다.',
     'Zeeman 평균을 사용한 축약 모델입니다. 전체 mF·편광 분해 모델과 동일하지 않습니다.','eq:pump-hamiltonian eq:age-normalization')
node('i-frequency','광학 carrier와 RF 축','Δ · δ · Ω_beat · Ω_SA',(2,1),
     'Δ=ω_pump−ω₁₃\nω_p=ω_pump−ω_HF+δ\nω_c=2ω_pump−ω_p\nΩ_beat=ω_HF−δ\nδa(t)=∫dΩ/(2π) δa(Ω)e^(−iΩt)',
     'pump/probe 주파수; 별도의 RF 분석 주파수 Ω_SA','식: rad/s; 수치 축: 주파수/(2π)',
     'Ω_beat는 Floquet 조화의 간격입니다. Ω_SA는 probe와 conjugate carrier 주위의 잡음 sideband를 지정합니다. p는 probe이고 pump의 약자가 아닙니다.',
     '평균장 δ 스캔을 RF 잡음 주파수 스캔으로 바꾸어 읽지 않습니다. 두 Raman branch는 별도 모드쌍입니다.','eq:frequency-definitions eq:fourier',compact='ω_c=2ω_pump−ω_p')
node('i-drive','입사 파워와 원자 환경','Rabi · N · 속도분포 · 완화율',(1,2),
     'I₀=2P/(πw²), Ω_ge=a_ge d E₀/ℏ\nf(v_z)∝exp[−v_z²/(2σ_v²)]\nσ_v=√(k_B T/m), N=N(T)',
     '셀에 입사한 pump/seed 파워, 1/e² 강도 반경 w, 온도 T','Ω·완화율: s⁻¹; N: m⁻³; σ_v: m/s',
     '파워와 빔 반경은 원자 구동을, 온도는 밀도·Doppler 가중치·충돌률을 바꿉니다. 생산 경로는 순수 ⁸⁵Rb 밀도를 사용합니다. EOM 전압과 측대역 변환은 셀 입사 입력을 만드는 실험 장치 쪽에 속합니다.',
     '모델의 transit 재장전과 ground/optical pure dephasing은 별개입니다. 현행 경험적 완화율에는 독립 측정에 따른 불확실성이 있습니다.','eq:rates-at-literature app:parameters')
node('i-modes','수집 모드와 입력 평균장','A_p · A_c · α_in · L · θ',(2,2),
     'P_j=ℏω_j|α_j|²\nα_in=(α_p,in, α_c,in*)ᵀ\nΩ=𝒢 α, M_can=𝒢⁻¹M_Ω𝒢',
     '수집 면적, probe/conjugate 모드, 셀 길이 L, 교차각 θ, 입력 carrier','A: m²; L: m; α: s⁻¹/²; 광자속: s⁻¹',
     'probe만 시드한 운전에서는 α_c,in=0입니다. 원자 코드의 Rabi 진폭과 정준 광자속 진폭은 변환이 필요합니다. 모드 면적은 실제 물리 입력이며 단위 규약만으로 정해지지 않습니다.',
     'conjugate 면적을 probe와 같다고 두면 조건부 모드 가정입니다. 파워만으로 복소 위상을 복원할 수 없습니다.','eq:flux-normalization eq:canonical-rescaling')
node('i-state','입력 양자 상태','V_in(Ω) · companion 입력',(1,3),
     'V_vac=I/2\ncoherent displacement: α≠0, V=I/2',
     '두 모드의 입력 요동과 상관, coherent displacement','정규화된 spectral covariance; 2×2 paired block / 4×4 quadrature',
     'coherent seed의 평균 진폭과 vacuum 크기의 요동을 따로 전달합니다. conjugate 진공도 독립 입력 포트입니다. 임의의 squeezed/thermal/correlated 입력에는 해당 V_in이 필요합니다.',
     '밝은 coherent seed+vacuum conjugate와 vacuum-input TMSV 시험은 서로 다른 초기 상태입니다.','eq:covariance eq:nambu-quadrature')
node('i-detection','셀 이후 측정 조건','η · H_p/c · w · S_el · SQL',(2,3),
     'η_j=η_path,j η_QE,j\nH_j=H_j(Ω), w=w(Ω)',
     '경로 투과·검출 QE·전자 응답·가중치·dark noise·RBW/VBW','η: 무차원; H: 보정된 전자 전달함수; S_el/SQL: 같은 PSD 단위',
     '외부 손실, 검출기의 주파수 응답, 전자 잡음은 별도 입력입니다. GABES의 총 효율 표시와 실험의 path/QE 분해는 같은 검출 경계를 기술합니다.',
     '셀 내부 흡수를 η에 다시 넣지 않습니다. 실험 SQL은 같은 검출 파워·가중치·분석 대역폭으로 보정해야 합니다.','eq:detection-loss eq:photocurrent tab:loss-ledger')

group('atomic','생산 원자 응답','두 기준장 Floquet 해 · 1D 평균','ready','II–III · Floquet harmonics and convergence',
      '유한 probe 기준장과 유한 conjugate 기준장으로 각각 풀어 응답행렬의 두 열을 만듭니다. 별도 pump-only 정상상태를 먼저 푸는 순서가 아닙니다.',
      '연분수·compiled kernel·dense block의 절단 문제 일치가 검사됩니다. 현재 스캔의 수렴 상태는 실제 실행 결과로만 판정합니다.')
node('a-assembly','구동과 산일자 조립','ℒ₀ · 𝒞₊ · 𝒞₋',(1,1),
     'ℒ[ρ]=−i[H/ℏ,ρ]+Σ D[L_μ]ρ\nℒ_reload[ρ]=γ_t(ρ_th Trρ−ρ)\nρ_th=diag(5/12,7/12,0,0)',
     '준위·가중 dipole·Rabi·Δ_eff·δ·branch·온도별 완화율','ℒ₀,𝒞±: 16×16, s⁻¹; ρ: trace-normalized',
     '자발방출 branching, optical/ground pure dephasing, 열적 재장전이 같은 Liouvillian에 들어갑니다. 재장전은 coherence만 감쇠시키는 채널이 아니라 population도 보충합니다.',
     '4준위·Markov·회전파 가정 안의 계산입니다. ground population은 ρ에서 한 번만 셉니다.','eq:pump-hamiltonian eq:rates-at-literature')
node('a-probe','probe 기준장 해','ρ₀ᵃ · ρ₊ᵃ',(1,2),
     '(ℒ₀+inΩ_beat)ρ_n+𝒞₊ρ_(n−1)+𝒞₋ρ_(n+1)=0\nprobe=Ω_p,ref; conjugate=0\nTrρ₀=1; Trρ_(n≠0)=0',
     '조립한 블록, 유한 probe 기준장, N_F','각 ρ_n: 4×4; n=−N_F,…,+N_F',
     'probe 입력에 대한 probe와 conjugate polarization을 같은 유한장 해에서 추출합니다. 기본 생산 차수는 N_F=3입니다.',
     '유한 기준장에 대한 비율은 무한소 도함수와 일반적으로 같지 않습니다.','eq:floquet-block',compact='probe 입력 → {ρ_n}ᵃ')
node('a-conjugate','conjugate 기준장 해','ρ₀ᵇ · ρ₊ᵇ',(2,2),
     '(ℒ₀+inΩ_beat)ρ_n+𝒞₊ρ_(n−1)+𝒞₋ρ_(n+1)=0\nprobe=0; conjugate=Ω_c,ref\nTrρ₀=1; Trρ_(n≠0)=0',
     '별도 기준장 블록, 유한 conjugate 기준장, 같은 N_F','각 ρ_n: 4×4',
     '응답행렬의 다른 열을 얻기 위한 두 번째 계산입니다. 실제 실험에서 conjugate를 함께 시드한다는 뜻은 아닙니다.',
     'probe 해와 응답을 합쳐 행렬을 조립하지만, 두 Raman branch의 χ를 합산하지 않습니다.','eq:floquet-block',compact='conjugate 입력 → {ρ_n}ᵇ')
node('a-readout','가중 polarization 읽기','4개 복소 reduced response',(1,3),
     'χ̄_pp=Pread_p(ρ₀ᵃ)/Ω_p,ref\nχ̄_cp=Pread_c(ρ₊ᵃ)/Ω_p,ref\nχ̄_pc=Pread_p(ρ₀ᵇ)/Ω_c,ref\nχ̄_cc=Pread_c(ρ₊ᵇ)/Ω_c,ref',
     '두 해의 ρ₀,ρ₊; ground별 CG-weighted readout','χ̄: s (Rabi angular rate로 나눈 coherence)',
     'bar는 reduced response를 뜻합니다. 속도 평균은 별도로 ⟨χ̄⟩_v라고 씁니다. 코드의 s 하첨자는 이 문서의 probe p에 대응합니다. Nambu 행렬 조립에서는 conjugation 규약을 적용합니다.',
     '구동과 polarization 각각의 진폭 가중치가 곱으로 전이 강도를 만듭니다. 읽기 가중치를 거시 변환에서 다시 적용하지 않습니다.','eq:susceptibility-dimension eq:chi-matrix eq:age-normalization',compact='ρ / Ω_ref → χ̄ [s]')
node('a-average','생산 1D 속도 평균','⟨χ̄_jk⟩_v',(1,4),
     '⟨χ̄_jk⟩_v=∫dv_z f(v_z)χ̄_jk(Δ−k_pump v_z,δ)\nexp(L⟨M⟩) ≠ ⟨exp(LM)⟩',
     '4개 복소 응답표, Δ_eff 보간, 1D Maxwell 가중치','복소 2×2 reduced response; 각 성분 s',
     '같은 광학장을 구동하는 속도군의 복소 응답을 먼저 평균합니다. 2D 비공선 Doppler는 별도 참조에 있으며 생산 평균을 자동 교체하지 않습니다.',
     '각 속도군을 통과한 전달행렬을 평균하는 모델과 다릅니다. 평균 응답은 저장고 diffusion을 공급하지 않습니다.','sec:doppler sec:pole-doppler',compact='χ̄(v_z) → ⟨χ̄⟩_v')
node('a-order','인접 차수 전체 경로','N_F=3 / 2의 응답·전파',(2,4),
     'N_F=3: χ̄₃ → T₃\nN_F=2: χ̄₂ → T₂',
     '같은 운전점·스캔·정규화에서 두 차수의 계산','동일 차원과 기저의 비교쌍',
     '각 차수에서 응답뿐 아니라 평균장 전파까지 계산합니다. 검증은 full scan의 모든 보고 성분과 feature에 붙습니다.',
     'trace·positivity 또는 단일 최적점의 일치만으로 절단 수렴을 선언하지 않습니다.','sec:floquet-convergence sec:physicality-blind',kind='test')
node('a-gate','전체 스캔 수렴 판정','CONVERGED / UNCONVERGED',(1,5),
     '|ΔT|≤10⁻¹⁰+0.01 max(|T_N|,|T_(N−1)|)\nχ̄의 절대 바닥값: 성분별 full-scan max의 10⁻⁸\ngain 변화 <1%; retained wrapped phase <0.5°\nprobe-gain 최적점 이동 ≤ 최종 격자 한 간격',
     '두 차수의 복소 응답·정준 T·gain·wrapped phase·gain 최적점','transfer: 무차원; response: s; phase: 각도',
     '차원이 다른 reduced response에 transfer의 절대 허용오차를 그대로 쓰지 않습니다. gate 실패는 UNCONVERGED이며 이 정적 HTML은 실제 운전점에 대한 합격을 주장하지 않습니다.',
     '인접 차수 일치는 선언된 수치 기준의 통과입니다. 원자 축약 가정이나 절대 gain의 실험 타당성을 보증하지 않습니다.','sec:floquet-convergence sec:physicality-blind',kind='test',compact='응답 + T + gain + phase + 최적점')

group('reference','pump-only 약응답 참조','정적 상태 · pole/residue · 2D 선택','reference','II–III · Pump-only state and weak-field response',
      '공유 원자 조립 위의 별도 실행 경로입니다. 생산 finite-seed와의 관계는 gauge 비교와 seed→0 극한으로 검사합니다. minus만 지원합니다.',
      '선언된 minus 참조 범위에서 수치 검증되었습니다. 7.32e−17은 자기 방정식 잔차이고, 유한 seed와의 차이는 별도로 보고합니다.')
node('r-frame','minus 정적 pump 프레임','H_pump · ℒ_pump',(1,1),
     'U₋(t)=exp(+iΩ_beat t|g₂⟩⟨g₂|)\nρ_pump=U₋†ρ_seed U₋\nH_pump/ℏ=ω_HF|g₂⟩⟨g₂|−Σ_e Δ_e(v)|e⟩⟨e|\n             +½Σ_ge(Ω_ge^pump|e⟩⟨g|+h.c.)',
     '같은 4준위 원자 조립, pump 구동, minus branch','H/ℏ: s⁻¹; ℒ_pump: 16×16',
     'g₂ 대각은 δ+Ω_beat=ω_HF가 되고 pump 결합이 정적이 됩니다. phase-invariant dissipator와 thermal reload는 유지됩니다.',
     'plus 참조는 gauge parity가 성립하지 않아 unsupported입니다. 생산 plus 응답 자체가 없다는 의미는 아닙니다.','eq:pump-frame-transform eq:pump-hamiltonian',compact='U₋ → ℒ_pump')
node('r-state','trace-one 정상상태','ρ_ss',(1,2),
     'ℒ_pump ρ_ss=0, Trρ_ss=1\n16×16 계의 한 행 → trace constraint',
     'ℒ_pump와 trace 구속','ρ_ss: 4×4, 무차원',
     'population과 coherence가 상호 피드백하는 정상상태입니다. 유한 seed가 없는 pump 상태에서만 참조 약응답을 정의합니다.',
     '현재 참조는 undepleted이고 pump 상태를 z 구간마다 다시 풀지 않습니다.','eq:pump-nullspace')
node('r-gauge','프레임·seed 극한 검증','서로 다른 비교량',(2,2),
     'g₂ row coherence → n=−1; column → n=+1\n나머지 → n=0; |n|≥2 → 0\ngauge parity: 8.37×10⁻¹⁵\nfinite seed 8 µW 응답 차이: 2.503×10⁻³',
     '정적 ρ_ss, pump-only Floquet 해, 별도의 finite-seed 응답','정규화 오차; 원문 보고값',
     'gauge 비교는 pump 상태끼리 비교합니다. seed 감소에 따른 네 복소 응답의 차이는 별도 극한 검증입니다. 마지막 세 점의 Rabi 진폭 기울기 1.749는 정확한 2차 법칙이 아닙니다.',
     '공유 원자 모델 안의 검사이며 독립 실험 검증은 아닙니다.','tab:pump-reference-diagnostics tab:finite-seed-reference-limit',kind='test')
node('r-source','trace-zero 구동과 주파수','−𝒱_kρ_ss · ω_k',(1,3),
     '𝒱_kρ_ss=−i[V_k,ρ_ss], Tr(𝒱_kρ_ss)=0\nω_k=−ω_HF+δ+Ω_SA',
     'ρ_ss, probe/conjugate sideband 섭동, 독립 Ω_SA','source: vectorized 16성분; ω_k: rad/s',
     'ω_k는 정적 pump 프레임에서의 resolvent 주파수입니다. optical δ와 RF Ω_SA가 독립적으로 들어갑니다.',
     'ω_k=0과 Ω_SA=0은 일반적으로 다릅니다. 정상 모드가 여러 개면 전체 stationary subspace를 처리합니다.','eq:linear-resolvent eq:pump-frame-response-frequency sec:zero-pole')
node('r-direct','구속 직접 선형해','δρ_k(Ω_SA)',(1,4),
     'δρ_k=−(ℒ_pump+iω_k)⁻¹_(Tr=0)𝒱_kρ_ss',
     'trace-zero source와 ℒ_pump','δρ: 4×4; 입력 Rabi에 대한 도함수는 s',
     'trace-zero 부분공간에서 선형계를 풉니다. 영 주파수의 stationary subspace는 projection 또는 Drazin 구성으로 다룹니다.',
     '유일한 trace mode의 residue가 0이어도 전체 역행렬의 0/0을 그대로 계산하지 않습니다.','eq:linear-resolvent sec:zero-pole')
node('r-poles','쌍직교 pole와 residue','λ_m · A_m^(jk)',(2,4),
     'ℒ|r_m⟩=λ_m|r_m⟩, ⟨l_m|r_n⟩=δ_mn\nχ̄_jk(ω_k)=Σ_m A_m^(jk)/(λ_m+iω_k)\nA_m^(jk)=−⟨w_j|r_m⟩⟨l_m|𝒱_k|ρ_ss⟩\nλ_m=−Γ_m+iω_m',
     '같은 ℒ_pump·source·polarization readout','λ: s⁻¹; A: 해당 reduced-response residue',
     'pole 중심은 resolvent 축의 −Imλ, HWHM은 −Reλ입니다. RF 축에는 ω_k 관계로 옮깁니다. residue가 작은 좁은 pole는 측정에서 보이지 않을 수 있습니다.',
     '단순극 합은 대각화 가능할 때의 같은 선형 문제 표현입니다. residue와 속도 평균·검출을 거치지 않은 pole 폭은 스퀴징 대역폭이 아닙니다.','eq:biorthogonal eq:chi-pole-sum eq:pole-parametrization',compact='χ̄=Σ A_m/(λ_m+iω_k)')
node('r-defective','결함 행렬의 대안','고차극 · 직접/등고선 해',(2,5),
     '(J_λ+iωI)⁻¹=Σ_(r=0)^(p−1) (−1)^r N^r/(λ+iω)^(r+1)',
     '결함 또는 ill-conditioned 고유기저','원래 응답과 같은 차원',
     'Jordan block은 고차극을 만듭니다. 수치 계산에는 Jordan-free 직접 선형해 또는 contour/Riesz-projector 구성이 적절합니다.',
     '고유값의 중복만으로 결함이라 판정하지 않습니다. 원문의 해석 대안이며 모든 경로가 생산에 구현되었다는 뜻은 아닙니다.','eq:jordan-resolvent',kind='analytic',compact='단순극 실패 → 안정한 평가')
node('r-response','Nambu 응답 조립','χ′(Ω_SA;v)',(1,6),
     '(𝒫_p,𝒫_c*)ᵀ=ε₀[χ_pp,χ_pc;χ_cp*,χ_cc*](E_p,E_c*)ᵀ\n참조 reduced derivative 순서: (p,c*)',
     '직접해 또는 유효한 spectral 표현의 polarization 도함수','reduced derivative: 2×2, s; physical χ: 무차원',
     '두 입력에 대한 복소 도함수를 조립합니다. 미시 요동 drift를 만드는 원자 재료를 제공하지만 B나 D를 반환하지 않습니다.',
     '생산 finite-seed 응답과 같은 실행 경로가 아닙니다. 참조와 실제 구동 상태의 대응을 정당화해야 합니다.','eq:chi-matrix eq:linear-resolvent')
node('r-geometry','속도와 비공선 detuning','1D 기본 / 2D opt-in',(2,6),
     'Δ_eff=Δ−k_pump v_z\nδ_eff=δ+(k_pump−k_p cosθ)v_z−k_p sinθ v_x\n⟨X⟩_v=∬dv_z dv_x f(v_z)f(v_x)X(Δ_eff,δ_eff)',
     'T, θ, bare k, 1D 또는 독립 (v_z,v_x) Maxwell 격자','δ_eff: rad/s; f(v)dv: 무차원',
     '2D wrapper는 lab Ω_beat와 Ω_SA를 고정하고 원자 detuning만 이동시킵니다. 121 °C·0.32°에서 Raman rms는 약 1.380 MHz입니다.',
     'δ_eff를 생산 beat 계산에 그대로 넣으면 다른 문제를 풉니다. 2D는 minus 참조에 한정됩니다.','eq:doppler-detunings eq:two-dimensional-velocity tab:angular-doppler-convergence')
node('r-average','복소 응답의 속도 구적','⟨χ′(Ω_SA)⟩_v',(1,7),
     '⟨χ′⟩_v=Σ_v w_v χ′(Ω_SA;v)',
     '속도별 full complex response, 정규화 Maxwell 가중치','2×2 reduced derivative, s',
     '속도별 population과 pole·residue 변화를 함께 평균합니다. 기본 1D와 선택 2D는 명시된 기하 안에서 별도로 수렴 검사합니다.',
     '광학 pole의 Doppler smear, Raman Lorentzian HWHM, 비공선 Gaussian rms를 구별합니다. 두 폭을 임의 제곱합해 측정 bandwidth라 하지 않습니다.','sec:pole-doppler tab:angular-doppler-convergence')
node('r-voigt','단일극 Faddeeva 극한','조건부 해석 평균',(2,7),
     '⟨1/(x₀−kv+iΓ_o)⟩=−i√(π/2)/(kσ_v) · w_F(z)\nz=(x₀+iΓ_o)/(√2 kσ_v), w_F(z)=e^(−z²)erfc(−iz)\nIm⟨1/(x₀−kv+iΓ_o)⟩=−πV(x₀;kσ_v,Γ_o)',
     '속도에 선형인 단일 pole 위치, 상수 residue, 한 속도 성분','base integral: s; V: 단위 정규화 Voigt, s',
     '일반 속도 구적의 특수한 해석 축약입니다. w_F는 Faddeeva 함수이고 검출 전자 가중치 w와 다릅니다.',
     '속도별 population·Raman self-energy·두 detuning 기하는 보통 이 가정을 깨므로 모든 χ가 Voigt인 것은 아닙니다.','eq:faddeeva eq:voigt-correct',kind='analytic',compact='단일 선형 pole → Faddeeva')

group('mean','평균장 전파','복소 진폭 · 위상 · 두 이득','ready','IV · Semiclassical propagation and gain',
      '복소 χ와 기하에서 Maxwell 전파를 구성합니다. 평균장 진폭은 이득 진단과 검출 선형화 양쪽으로 전달됩니다.',
      'Maxwell/Q 구현 비교는 수치 일관성 검사입니다. 문헌 운전점에서 절대 gain의 큰 불일치가 남아 있습니다.')
node('d-chi','무차원 χ로 변환','physical χ_jk',(1,1),
     'χ_jk=−2Nd²ℓ_norm⟨χ̄_jk⟩_v/(ε₀ℏ)\na_ge=√(3C_F²), C_F²=(10,35,35,28)/81\nℓ_norm=(1/12)×경험적 잔여 계수×추가 결합 인자',
     '평균 reduced response, pure-⁸⁵Rb N, reduced D1 dipole, 결합 정규화','χ: 무차원; reduced response: s',
     '구동과 readout의 a_ge가 전이 강도를 만듭니다. ρ가 pump-modified population을 이미 포함하므로 별도 평형 p_F를 다시 곱하지 않습니다. 현행 잔여 계수 0.74는 미시 유도값이 아닙니다.',
     '추가 mode/polarization/Zeeman 참여 penalty의 기본값은 1입니다. UI에서 숨긴 내부 인자를 독립 측정으로 확정된 값처럼 읽지 않습니다.','eq:susceptibility-dimension eq:age-normalization eq:ell-effective-corrected')
node('d-geometry','bare k와 기하 위상부정합','k_p · k_c · Δk_geom',(2,1),
     'Δk_geom=2k_pump,z−k_p,z−k_c,z\nk_j=ω_j/c\nQ=diag(√[2ℏω_p/(ε₀cA_p)], √[2ℏω_c/(ε₀cA_c)])',
     'carrier 주파수, 교차각, 실제 수집 면적','k·Δk: m⁻¹; Q: 광자속 진폭 → 전기장',
     'Option A는 vacuum carrier k와 기하만 mismatch에 씁니다. Reχ의 분산은 Maxwell 대각에서 이미 들어갑니다.',
     '굴절률로 k를 다시 보정해 분산을 두 번 넣지 않습니다. mode area와 주파수 정규화를 생략하지 않습니다.','eq:electric-propagation eq:geometric-mismatch')
node('d-maxwell','Maxwell 행렬과 정준 변환','M_cl',(1,2),
     'M_E = [ +ik_pχ_pp/2−iΔk_geom/2,   +ik_pχ_pc/2 ;\n          −ik_cχ_cp*/2,             −ik_cχ_cc*/2+iΔk_geom/2 ]\nM_cl=Q⁻¹M_EQ\nα̃=(α_p e^(−iΔkz/2), α_c* e^(+iΔkz/2))ᵀ',
     'χ의 네 성분 + bare k/mismatch + Q','M_E,M_cl: 2×2, m⁻¹; paired basis (p,c*)',
     '대각 χ의 허수부는 흡수/증폭을, 실수부는 분산을 반영합니다. conjugate 성분의 복소켤레와 부호는 phase-conjugate 기저에서 옵니다.',
     '수동 uncoupled probe에서 |α_p(L)|²=e^(−k_p Imχ_pp L)|α_p(0)|²가 되어야 합니다.','eq:classical-propagation eq:electric-propagation',compact='χ + Δk + Q → M_cl')
node('d-constant','상수 drift의 행렬지수','T_cl=e^(M_cl L)',(1,3),
     's=(a+d)/2, q²=(a−d)²/4+bc\nT_cl=e^(sL)[cosh(qL)I + sinh(qL)/q · (M_cl−sI)]\nq=0: T_cl=e^(sL)[I+L(M_cl−sI)]',
     '상수 M_cl=[a,b;c,d], L','T: 2×2 무차원; s,q: m⁻¹',
     'Cayley–Hamilton의 (M_cl−sI)²=q²I로 지수급수가 닫힙니다. sinh(qL)/q의 연속값은 L이며 고유벡터 역행렬이 필요 없습니다.',
     'undepleted pump·일정 밀도·온도·overlap·하나의 수집 모드 가정입니다.','eq:matrix-exponential eq:cayley-hamilton',compact='T_cl=exp(M_cl L)')
node('d-segment','z 의존 전파','순서 있는 T_K⋯T₁',(2,3),
     'T_k=exp(M_cl,k Δz_k)\nT_total=T_K⋯T₂T₁',
     '국소 계수·분절 길이·공간 overlap','각 T_k: 2×2 무차원',
     '서로 다른 구간의 행렬은 일반적으로 교환하지 않습니다. 실제 국소 계수에 대한 piecewise-constant 전파를 순서대로 합성합니다. Ultra는 공간 overlap과 근사 pump 예산을 반영합니다.',
     '현재 근사 pump-budget 제한은 매 구간 원자 상태를 다시 푸는 완전한 세 장 고갈 해가 아닙니다.','eq:matrix-exponential eq:segmented-covariance',compact='T_total=T_K⋯T₁')
node('d-amplitude','출력 복소 평균장','α_p,out · α_c,out → 검출',(1,4),
     'α̃_out=T_cl α̃_in\nα̃_in=(α_p,in,0)ᵀ → probe-seed-only',
     '선택한 전파자, 입력 carrier, phase gauge','α: s⁻¹/²; 파워: ℏω|α|²',
     '출력의 진폭과 위상이 검출 photocurrent의 방향을 정합니다. 측정한 carrier로 대체할 수도 있지만 그때는 측정값을 조건으로 하는 모델입니다.',
     'post-hoc power cap으로 얻은 gain만으로 이 복소 평균장을 유일하게 복원하지 않습니다. M_cl의 T를 일반적인 T_q(Ω)로 동일시하지 않습니다.','eq:mean-field-gains eq:photocurrent')
node('d-gains','두 이득과 정규화','G_p · G_c · T_small-signal',(2,4),
     'G_p=P_p,out/P₀=|T₁₁|²\nG_c=P_c,out/P₀=(ω_c/ω_p)|T₂₁|²',
     '정준 작은 신호 T와 같은 probe-seed 입력 P₀','gain: 무차원 power ratio',
     '두 gain은 독립적으로 계산합니다. photon-flux gain과 power gain을 구별합니다. 현행 seed-dependent 유효 power budget 출력은 별도 provenance를 갖습니다.',
     '정준 진단에는 작은 신호 선형 T를 사용합니다. capped/유효 gain을 그 선형 T의 원소처럼 역으로 해석하지 않습니다.','eq:mean-field-gains eq:gain-gap-observable')

group('gain','이득 진단','광자속 관계 · 조건부 이상 극한','analytic','IV / VII · Gain convention and ideal IDS',
      '평균장 gain에서 검사할 수 있는 것과 이상 증폭기의 비교식을 읽습니다. 공분산을 만드는 경로는 아닙니다.',
      '스칼라 gain 관계는 전체 정준성·입력 상태·검출 조건을 검증하지 않습니다. 현재 mean-field indicator는 물리적 스펙트럼이 아닙니다.')
node('x-flux','광자속 gap','𝒢_flux',(1,1),
     'g_c^(Φ)=|T₂₁|²=(ω_p/ω_c)G_c\n𝒢_flux=G_p−g_c^(Φ)\n𝒢_power=G_p−G_c',
     '두 gain과 광학 주파수','무차원',
     '두 gap은 출력입니다. near-degenerate 파장에서 수치적으로 가까워도 정준 metric에 대응하는 것은 flux gap입니다.',
     '𝒢_flux≠1이면 두 연산자만의 닫힌 정준 증폭기 해석에 추가 채널이 필요합니다.','eq:gain-gap-observable eq:gain-gap')
node('x-canonical','전체 정준성 조건','TJT†=J와 입력 가정',(2,1),
     'J=diag(1,−1), TJT†=J\n반례: T=diag(1,2)\n𝒢_flux=1 이지만 TJT†=diag(1,−4)',
     '두 gain뿐 아니라 full complex T, 선택한 모드 정규화','2×2 metric identity',
     'flux gap 하나가 맞는다고 정준성 전체가 성립하지 않습니다. 이상 IDS에는 밝은 coherent seed, vacuum conjugate와 지정 검출 조건도 필요합니다.',
     '열린 산일 셀의 교환관계는 저장고 항까지 포함해 검사합니다.','eq:commutator-integral eq:gain-gap',kind='test')
node('x-ideal','밝은 seed 이상 비교','S_ideal · 대칭 외부 손실',(1,2),
     'S_ideal=1/(G_p+g_c^(Φ))=1/(2G_p−1)\nS(η)=1−η+ηS_ideal',
     '무손실 정준 증폭기, 밝은 coherent seed, vacuum conjugate, 비가중 차동','SQL-normalized noise, 무차원',
     '두 번째 식은 대칭 외부 손실일 때입니다. 비대칭 효율·전자 응답·가중치에서는 full covariance와 같은 조건의 SQL을 사용합니다.',
     '현재 GABES gain 지표의 일반 계산식을 뜻하지 않습니다. flux gap 검사만 통과했다고 이 이상 가정들이 성립하는 것은 아닙니다.','eq:ids-ideal eq:ids-symmetric-loss eq:ids-balanced')

group('micro','미시 요동과 저장고','M_q · B · N_res · D','need','V · Quantum-Langevin propagation',
      '구동된 원자의 주파수 응답과 미시 저장고 상관이 필요합니다. 평균한 유한 seed 응답만으로 잡음을 정할 수 없습니다.',
      '현재 생산에 주파수 의존 미시 M_q와 원자 D의 완결된 연결은 없습니다. 참조 χ′는 재료이며 독립진공 완성은 시험용입니다.')
node('m-pair','Nambu pairing','Â(+Ω,−Ω) · J',(1,1),
     'Â(z,Ω)=(δâ_p(z,Ω),δâ_c†(z,−Ω))ᵀ\n[Â(Ω),Â†(Ω′)]=2πJ δ(Ω−Ω′)\nJ=diag(1,−1)',
     '수집 두 모드·정규화·Fourier 부호','two-component paired spectral block',
     'FWM은 반대 RF offset의 phase-conjugate 요동을 결합합니다. companion pairing을 함께 풀어 검출의 네 quadrature를 구성합니다.',
     '두 성분에 같은 creation/annihilation metric을 적용하지 않습니다. 모드 정규화가 바뀌면 metric을 사용한 진단도 재계산합니다.','eq:nambu eq:canonical-rescaling',kind='analytic')
node('m-drift','원자 요동 drift 구성','M_q(z,Ω)',(1,2),
     '∂_z Â=M_q(z,Ω)Â+B(z,Ω)F̂\nM_q(z,Ω) ≠ M_cl  (일반적인 경우)',
     '구동 원자 모델·상태, 해당 상태에서의 주파수 응답, 모드 결합','M_q: 2×2, m⁻¹',
     'pump-only χ′에서 실제 구동 상태의 요동 drift로 연결하려면 상태·약요동·전파 근사를 정당화해야 합니다. 생산 finite-seed 응답에서 자동으로 얻는 값이 아닙니다.',
     '현재의 원자 참조는 미시 양자 propagation solver 전체를 제공하지 않습니다.','eq:langevin eq:linear-resolvent')
node('m-reservoir','미시 저장고 상관','B · J_F · N_res',(2,2),
     '[F̂(z,Ω),F̂†(z′,Ω′)]=2πJ_F δ(z−z′)δ(Ω−Ω′)\n½⟨{F̂(z,Ω),F̂†(z′,Ω′)}⟩=2πN_res δ(z−z′)δ(Ω−Ω′)',
     '원자 Lindblad/Einstein 관계 또는 독립 근거가 있는 reservoir 모델','B: 2×r; J_F,N_res: r×r; 공간 white-noise 정규화',
     '교환관계와 symmetrized 상관은 다른 정보입니다. 실제 자발방출·원자 population·구동 상태의 상관이 이 입력을 정해야 합니다.',
     'K의 고유값만으로 symmetrized N_res를 추론하지 않습니다. 현재 미시 원자 diffusion은 미완성입니다.','eq:physical-realizability eq:covariance')
node('m-k','교환관계 채널 분류','K=BJ_FB†',(1,3),
     'K=−(M_qJ+JM_q†)\nM_qJ+JM_q†+BJ_FB†=0\nTJT†+∫T(L,z)BJ_FB†T(L,z)† dz=J',
     '정준 M_q와 J; reservoir commutator','K: Hermitian 2×2, m⁻¹',
     'K의 양 고유값은 loss-type, 음 고유값은 gain-type 채널을 요구합니다. 둘 다 증폭기의 물리적 reservoir 구성에 나타날 수 있습니다.',
     'K 부호는 채널 분류이며 합격/불합격이 아닙니다. 모드 면적이나 Rabi 변환이 달라지면 K도 바뀝니다.','eq:K-matrix eq:physical-realizability eq:commutator-integral',kind='analytic')
node('m-diffusion','symmetrized diffusion','D(z,Ω)',(2,3),
     'D=BN_resB†\n독립 속도군: D_total=Σ_v D_v  (가중치 포함)',
     '미시 B와 symmetrized N_res, 속도군별 공분산 기여','D: 2×2, covariance per unit length',
     '독립 속도군은 diffusion power로 더합니다. 속도평균한 잡음 진폭을 제곱하는 방식은 상관 구조를 바꿉니다.',
     'K는 antisymmetric reservoir constraint를 정하며 D를 유일하게 정하지 않습니다. 실제 미시 noise 입력이 필요합니다.','eq:covariance')
node('m-vacuum','독립진공 완성 시험','선택한 D_vac',(1,4),
     'K=UΛ_KU†, B_vac=U|Λ_K|^(1/2)\nJ_F=diag(signΛ_K)\nD_vac=½B_vacB_vac†',
     '상수 drift의 K; 독립 vacuum 채널이라는 추가 선택','D_vac: 2×2, covariance/m',
     '정해진 K를 만족시키는 대수 fixture입니다. 상수 drift와 정적 검출에 넣으면 평탄한 잡음 진단이 되며 bandwidth를 예측하지 않습니다.',
     '원자 자발방출 noise의 유도나 보편적 최저 잡음이 아닙니다. 이 시험의 통과로 미시 연결 미완성 상태를 해제하지 않습니다.','eq:minimal-dilation',kind='test')

group('covariance','공분산 전파','입력 요동 + 분포 잡음의 합류','analytic','V · Covariance and Lyapunov solution',
      '입력 요동을 전파한 항과 셀의 모든 위치에서 주입된 잡음을 더합니다. 이 계산 구조는 M_q와 D가 공급될 때 평가할 수 있습니다.',
      '원문의 대수 fixture에서 직접 구적·φ₁·Van Loan을 비교했습니다. 그 검증은 실제 원자 diffusion의 검증과 다릅니다.')
node('v-propagator','요동 전파자','T_q(z₂,z₁;Ω)',(1,1),
     'T_q(z₂,z₁;Ω)=𝒫 exp[∫_(z₁)^(z₂) M_q(z,Ω)dz]\n상수 M_q: T_q(L)=exp(M_qL)',
     '미시 M_q(z,Ω)와 L','T_q: 2×2 무차원',
     'T_q는 입력 요동과 각 위치의 reservoir contribution을 모두 전파합니다. 평균장 T_cl을 재사용하려면 두 drift가 같다는 별도 조건이 필요합니다.',
     'z 의존 계수는 path ordering이 필요합니다.','eq:langevin-solution')
node('v-input','입력 요동의 전파','T_q V_in T_q†',(1,2),
     'V_input,out=T_q V_in T_q†',
     'T_q(L,0;Ω), 외부 입력 상태 V_in(Ω)','2×2 paired covariance',
     'coherent displacement와 별개인 입력 covariance를 변환합니다. coherent+vacuum 입력의 경우 정규화된 V_in=I/2입니다.',
     '이 항만으로 열린 산일 셀의 출력 covariance가 완성되지는 않습니다.','eq:covariance')
node('v-integral','셀 내부 분포 잡음','W(L;Ω)',(2,2),
     'W=∫₀ᴸ T_q(L,z)D(z,Ω)T_q(L,z)† dz',
     '같은 T_q(L,z;Ω), 미시 D(z,Ω)','W: 2×2 paired covariance',
     '각 z에서 생성된 잡음은 남은 구간을 다시 통과합니다. 광학 흡수와 gain이 섞여 있으므로 전체 잡음을 출구 beam splitter 하나로 대체하지 않습니다.',
     '독립 local Markov 저장고를 가정합니다. 구간 사이 상관이 있으면 추가 cross term이 필요합니다.','eq:covariance eq:W-integral')
node('v-lyapunov','상수계수 Lyapunov 문제','𝓜 · W(0)=0',(1,3),
     '𝓜[X]=M_qX+XM_q†\ndW/dL=𝓜[W]+D, W(0)=0',
     'z에 독립인 M_q와 D','𝓜: 행렬공간 선형 map, m⁻¹',
     '왼쪽/오른쪽 곱이 교환하므로 exp(u𝓜)[D]=exp(uM_q)D exp(uM_q†)입니다. 아래 세 표현은 같은 W를 평가하는 대안입니다.',
     '세 결과를 순차적으로 통과하거나 합산하지 않습니다.','eq:lyapunov-superoperator eq:differential-lyapunov')
node('v-phi','좌표무관 φ₁ 해','W=Lφ₁(L𝓜)[D]',(1,4),
     'φ₁(x)=(eˣ−1)/x=Σ_(k≥0) x^k/(k+1)!\nφ₁(0)=1\nW=Lφ₁(L𝓜)[D]',
     '상수계수 Lyapunov 문제','W: 2×2 covariance',
     'φ₁은 entire function이므로 finite L에서 영 고유값이나 결함점의 대수적 특이성이 없습니다. 안정한 행렬함수 평가를 사용합니다.',
     '이 사실이 exceptional point의 물리적 mode coalescence나 민감도까지 없앤다는 뜻은 아닙니다.','eq:lyapunov-closed-form',compact='W=Lφ₁(L𝓜)[D]')
node('v-eigen','조건부 고유기저 전개','μ_m+μ_n* · residue',(2,4),
     'M_q=S diag(μ_m)S⁻¹, D̃=S⁻¹DS^(−†)\nW̃_mn=D̃_mn Lφ₁[(μ_m+μ_n*)L]\nW=S W̃ S†; spec𝓜={μ_m+μ_n*}\ns_q=tr(M_q)/2, q_q²=(a_q−d_q)²/4+b_qc_q\nμ_±=s_q±q_q',
     '대각화 가능하고 조건수가 양호한 M_q, D','μ,s_q,q_q: m⁻¹',
     '각 채널이 확산을 어떻게 증폭·간섭시키는지 해석하는 식입니다. 허수 교차 exponent와 실수 gain exponent가 함께 나타납니다.',
     'cond(S)를 확인합니다. 결함·준결함점에서는 φ₁/블록지수로 평가합니다. classical (s,q)와 값이 같다는 가정은 하지 않습니다.','eq:lyapunov-eigenbasis eq:lyapunov-spectrum eq:Mq-eigenvalues',compact='W̃_mn=D̃_mn Lφ₁[(μ_m+μ_n*)L]')
node('v-vanloan','Van Loan 블록지수','고유벡터 없이 W 평가',(1,5),
     'Z=[−M_q, D; 0, M_q†]\nexp(ZL)=[F₁₁,F₁₂;0,F₂₂]\nW=F₂₂†F₁₂',
     '상수 M_q와 D','Z: 4×4 block matrix; W: 2×2',
     '한 번의 4×4 행렬지수로 같은 분포 잡음 적분을 평가합니다. 원문에서 권장하는 구현 방식입니다.',
     '고유기저의 불안정성을 피하지만 미시 D의 부재를 해결하는 알고리즘은 아닙니다.','eq:van-loan-block eq:van-loan-result',compact='W=F₂₂†F₁₂')
node('v-segments','z 의존 국소 채널 합성','V^(k+1)',(2,5),
     'T_k=exp(M_k Δz_k), W_k=W(Δz_k;M_k,D_k)\nV^(k+1)=T_k V^(k) T_k†+W_k\n(T₂,W₂)∘(T₁,W₁)=(T₂T₁,T₂W₁T₂†+W₂)',
     '실제 국소 M_k,D_k, 구간 길이, 초기 V_in','각 T_k: 무차원; W_k,V: covariance',
     '각 구간에서 평가한 W_k를 다음 구간에 다시 전파합니다. 상수계수 전체 W를 얻은 뒤 잡음을 추가로 한 번 더 넣는 과정이 아닙니다.',
     'piecewise-constant 계수와 독립 구간 저장고 조건에서 합성합니다. 평균 drift의 지수로 임의 대체하지 않습니다.','eq:segmented-covariance')
node('v-output','두 항의 합류','2×2 V_out',(1,6),
     'V_out=T_q V_in T_q†+W\nW=W†; D⪰0 ⇒ W⪰0; L=0 ⇒ W=0',
     '입력 전파와 W, 또는 z-분절 경로의 최종 V','2×2 paired covariance',
     '이 합류가 원자 내부 noise를 포함한 한 Nambu block의 결과입니다. 형식의 일치·Hermiticity·PSD는 수치 invariants입니다.',
     '2×2 V_out을 4×4 detector 행렬에 곧바로 곱하지 않습니다.','eq:covariance eq:differential-lyapunov')
node('v-quadrature','companion과 quadrature 조립','4×4 V_R,out',(1,7),
     'Ξ(Ω)=(a_p(Ω),a_c(Ω),a_p†(−Ω),a_c†(−Ω))ᵀ\nC=(1/√2)[1,0,1,0; −i,0,i,0; 0,1,0,1; 0,−i,0,i]\nV_R=C V_Ξ C†; R=(X_p,Y_p,X_c,Y_c)ᵀ',
     '두 paired block, companion 입력/저장고, ±Ω 상관 조립','V_Ξ,V_R: 4×4; vacuum I₄/2',
     '반대 sideband companion과 필요한 상관을 조립해 검출 quadrature로 변환합니다. 네 벡터 성분은 두 carrier의 annihilation/creation bookkeeping입니다.',
     '네 독립 광학 carrier로 해석하지 않습니다. real-mode Gaussian CP 조건의 Ω_R도 RF 주파수 Ω와 다른 symplectic 행렬입니다.','eq:nambu-quadrature',compact='paired blocks → V_Ξ → V_R')

group('detector','외부 손실과 차동 검출','평균장 + 공분산 + 전자 응답','analytic','VI · External detection model',
      '공분산과 carrier에 같은 외부 손실을 적용한 뒤, 검출된 평균장으로 측정 벡터와 SQL을 함께 만듭니다. 아래의 TMSV 손실 통계는 별도 이상 시험입니다.',
      '검출 대수는 주어진 상태·효율·전자 응답 조건에서 성립합니다. 필요한 미시 covariance와 실험별 검출 보정은 별도 입력입니다.')
node('f-loss','공분산의 외부 순수 손실','V_R,det',(1,1),
     'K_η=diag(√η_p,√η_p,√η_c,√η_c)\nK_ℓ=diag(√(1−η_p),√(1−η_p),√(1−η_c),√(1−η_c))\nV_R,det=K_η V_R,out K_ηᵀ+K_ℓ(I₄/2)K_ℓᵀ',
     '4×4 셀 출력 V_R,out, 외부 η_p/c, 독립 vacuum ports','검출 V: 4×4 quadrature covariance',
     'cell 이후의 path loss와 detector QE는 vacuum beam splitter로 나타냅니다. 내부 absorption/spontaneous emission은 이미 분포 drift와 noise에 속합니다.',
     'exp(−OD_cell)를 여기서 다시 적용하지 않습니다. 실제 thermal 환경이면 vacuum I/2 대신 그 환경 covariance가 필요합니다.','eq:detection-loss eq:pure-loss-operator',compact='V_det=K_ηV_outK_ηᵀ+vacuum')
node('f-carrier','검출 carrier와 위상','α_p,det · α_c,det',(2,1),
     'α_j,det=√η_j α_j,out\nP_j,det=η_jℏω_j|α_j,out|²',
     'mean 그룹의 복소 carrier 또는 명시적으로 측정한 carrier, 같은 외부 η','α_det: s⁻¹/²; power: W',
     '공분산 손실과 같은 검출 경계를 적용합니다. 출력 평균장을 직접 측정값으로 대체하면 그 측정값에 조건부인 예측이 됩니다.',
     '전자 식에서 α_det를 쓸 때 η를 다시 곱하지 않습니다. mean-field power cap은 carrier 위상을 제공하지 않습니다.','eq:flux-normalization eq:photocurrent')
node('f-current','광전류 선형화','측정 벡터 m(Ω)',(1,2),
     'δi_j∝H_j(Ω)[α_j,det* δa_j,det(Ω)+α_j,det δa_j,det†(−Ω)]\nδi_−=δi_p−w(Ω)δi_c\nδi_−=c(Ω)ᵀ R_det, m=c*\nc=√2(H_p Reα_p, H_p Imα_p, −wH_c Reα_c, −wH_c Imα_c)ᵀ',
     '검출 carrier 진폭·위상, H_p/H_c, w; α는 모두 검출값','c,m: 4성분; 전자 gain 보정은 동일하게 적용',
     '측정 벡터는 상수가 아니라 carrier와 전자 응답에서 구성됩니다. R=(X_p,Y_p,X_c,Y_c), δa=(X+iY)/√2 규약으로 전개했습니다.',
     '밝은 carrier 주위의 선형화입니다. cᵀV c*=m†V m이므로 복소 주파수 응답의 conjugation을 유지합니다.','eq:photocurrent eq:spectrum',compact='α_det + H + w → m')
node('f-sql','같은 검출 조건의 SQL','S_SQL(Ω)',(2,2),
     'S_SQL∝|H_p|²|α_p,det|²+|w|²|H_c|²|α_c,det|²\n실험 SQL: 동일 power, H, w, RBW, VBW',
     '같은 검출 carrier·전자 전달·가중치·실험 보정','전자 PSD; S_el과 같은 단위',
     '두 팔의 독립 coherent shot noise를 같은 차동 측정 조건으로 정규화합니다. w가 복소수라면 |w|²이고 실수일 때만 w²입니다.',
     'DC power 균형 가중치가 항상 normalized noise의 최소 가중치는 아닙니다. dark noise 포함/차감 방식을 명시해야 합니다.','eq:photocurrent eq:spectrum')
node('f-measure','검출 covariance 투영','S_raw(Ω) · SQL',(1,3),
     'S_raw(Ω)=m†V_R,det(Ω)m+S_el(Ω)',
     'V_R,det, m, 보정된 전자 PSD, 같은 조건의 SQL','분자와 분모: 같은 electronic PSD 단위',
     '셀 출력의 covariance와 복소 carrier를 여기서 함께 사용합니다. mean-field gain만으로 이 contraction을 만들 수 없습니다.',
     '미시 covariance가 없으면 물리적 스펙트럼은 미완성입니다. pump scatter 등 기술 잡음은 위치와 보정이 확인된 입력으로 다룹니다.','eq:spectrum tab:loss-ledger')
node('f-tmsv','이상 TMSV의 공유 pair 수','thermal marginal',(1,5),
     '|ψ⟩=√(1−λ)Σ_n λ^(n/2)|n,n⟩\nλ=tanh²r, n̄=sinh²r\nP(N=n)=n̄ⁿ/(1+n̄)^(n+1)\nVarN=n̄(n̄+1), g²=2',
     '각 팔 한 시공간 모드, vacuum-input 이상 two-mode squeezer','N: photon count; n̄: 단일 모드 평균 점유',
     '전체 순수 상태라도 한 팔만 보면 thermal marginal입니다. 이 local uncertainty는 공유 pair number의 불확실성이며 독립 thermal field를 주입한 것이 아닙니다.',
     '별도 이상 광자수 시험입니다. 밝은 seeded beam에 g²=2를 고정하거나 선형 photocurrent에 zero-carrier 상태를 넣지 않습니다.','eq:tmsv-fock eq:tmsv-marginal',kind='test',compact='공유 N → 두 thermal marginal')
node('f-thinning','두 팔의 독립 순수 손실','K_j|N ~ Binomial(N,η_j)',(1,6),
     'a_j,det=√η_j a_j,out+√(1−η_j)v_j\nE[K_j|N]=η_jN\nVar(K_j|N)=η_j(1−η_j)N\nVarK_j=η_j²VarN+η_j(1−η_j)E[N]',
     '공유 N과 두 독립 vacuum loss port','K_j: 검출 photon count',
     '첫 분산 항은 줄어든 source fluctuation, 두 번째는 binomial partition noise입니다. vacuum이 Poisson photon population을 주입한다는 뜻은 아닙니다.',
     '순수 손실 후 단일 모드 marginal은 평균 ηn̄의 thermal 상태입니다. 평균이 0이 아니면 g²=2는 유지됩니다.','eq:pure-loss-operator eq:loss-thinning eq:loss-total-variance',kind='test',compact='VarK=η²VarN+η(1−η)E[N]')
node('f-difference','차동의 열적 요동 상쇄','대칭 / 비대칭 NRF',(1,7),
     'Cov(K_p,K_c)=η_pη_c n̄(n̄+1)\nVar(K_p−wK_c)=n̄[η_p(1−η_p)+w²η_c(1−η_c)]\n                  +(η_p−wη_c)² n̄(n̄+1)\nη_p=η_c=η, w=1: NRF=Var(K_p−K_c)/E[K_p+K_c]=1−η',
     'TMSV count moments, 효율, 실수 count 가중치 w','NRF: 무차원',
     '대칭 차동은 공유 pair fluctuation을 상쇄합니다. 비대칭이면 남는 열적 항이 생깁니다. 이 시험에서 w=η_p/η_c는 그 항을 지웁니다.',
     'count 가중치와 밝은 carrier의 RF 최적 검출 가중치는 같은 문제가 아닙니다. X의 bright-seed 식과 초기 상태가 다릅니다.','eq:two-arm-loss-statistics eq:tmsv-weighted-difference eq:tmsv-symmetric-nrf',kind='test',compact='대칭 순수 손실: NRF=1−η')
node('f-states','seeded 상태와 thermal 환경','서로 다른 물리 입력',(2,6),
     'seeded marginal: ρ_j=D(α_j)ρ_th(n̄_th)D†(α_j)\ng²=1+[n̄_th²+2n̄_th|α|²]/(n̄_th+|α|²)²\nthermal loss: V_det=ηV_out+(1−η)(n_E+½)I\nn̄_det=ηn̄_out+(1−η)n_E',
     'coherent displacement 또는 실제 환경 점유 n_E','n_E,n̄_th: mode occupancy',
     '밝은 시드의 marginal은 displaced thermal이며 carrier가 우세하면 g²→1입니다. 실제 thermal 환경은 n_E>0일 때 광자를 더합니다.',
     '일반 광학 순수 손실의 포트는 보통 n_E≈0입니다. 이 출구 loss model로 분포하는 원자 Langevin reservoir를 대체하지 않습니다.','eq:seeded-displaced-thermal eq:thermal-attenuator',kind='test',compact='displaced thermal ≠ vacuum loss port')

group('output','검출 스펙트럼','RF Ω_SA별 SQL 정규화','need','VI–VII · Assembled observable and validation',
      '최종 목표는 조건이 명시된 S_−(Ω_SA)입니다. 현재 mean-field Squeezing indicator의 δ 스캔과 이 물리적 RF 스펙트럼을 구별합니다.',
      '주파수 의존 미시 drift·diffusion과 실제 detector calibration이 필요합니다. 물리적 스퀴징과 bandwidth 예측은 아직 unavailable입니다.')
node('o-spectrum','정규화와 dB','S_−(Ω_SA)',(1,1),
     'S_−(Ω)=[m†V_R,det(Ω)m+S_el(Ω)]/S_SQL(Ω)\nS_dB(Ω)=10 log₁₀ S_−(Ω)',
     '검출 covariance 투영, 전자 noise, 동일 조건 SQL','S: 무차원; S_dB: dB; RF 축: Ω_SA/(2π)',
     'S<1은 정의한 측정·SQL 조건에서의 sub-SQL noise입니다. 단일 scalar gain이나 원자 pole 폭으로 이 주파수 함수를 대신하지 않습니다.',
     '입력 상태·효율·가중치·분석 대역폭·dark 처리·구동 조건을 함께 명시해야 합니다.','eq:spectrum')
node('o-claim','예측의 성립 범위','필요 입력과 검증 경계',(2,1),
     '원자 M_q(Ω), D(Ω) + 입력 V_in + 검출 H,w,η,SQL\n→ 물리적 spectrum 비교에 필요한 입력',
     '모델 provenance, 수치 수렴, 미시 noise, 독립 실험 비교','상태·검증 주장; 계산 숫자가 아님',
     '대수적 일치, 내부 방정식 잔차, 절단 수렴, 실험적 예측 오차를 따로 읽습니다. 준위 축약과 결합 정규화의 불확실성도 절대 gain에 영향을 줍니다.',
     '현재 절대 gain은 실험과 불일치하고 미시 noise는 미완성입니다. 정적 vacuum fixture의 평탄한 곡선은 bandwidth 예측이 아닙니다.','tab:validation-hierarchy tab:model-layer-comparison')

# Edges are the single source for SVG wires, mobile ports, text ledgers and JSON.
for s,t,l in [
 ('i-levels','a-assembly','준위·전이 강도'),('i-frequency','a-assembly','Δ·δ·branch·beat'),
 ('i-drive','a-assembly','Rabi·N·완화율'),('i-drive','a-average','속도 가중치'),
 ('i-frequency','d-geometry','bare carrier 주파수'),('i-modes','d-geometry','L·θ·면적'),
 ('i-modes','d-amplitude','입력 복소 carrier'),('i-state','v-input','V_in'),
 ('i-state','v-segments','초기 V_in'),('i-state','v-quadrature','companion 입력 상태'),
 ('i-detection','f-loss','η·vacuum port'),('i-detection','f-carrier','동일 η'),
 ('i-detection','f-current','H·w'),('i-detection','f-sql','SQL 보정'),('i-detection','f-measure','S_el'),
 ('a-assembly','a-probe','probe 기준장 블록'),('a-assembly','a-conjugate','conjugate 기준장 블록'),
 ('a-probe','a-readout','ρ₀ᵃ·ρ₊ᵃ'),('a-conjugate','a-readout','ρ₀ᵇ·ρ₊ᵇ'),
 ('a-readout','a-average','네 복소 χ̄'),('a-average','d-chi','평균 reduced response'),
 ('i-levels','r-frame','공유 원자 조립'),('i-drive','r-frame','pump·완화율'),
 ('r-frame','r-state','ℒ_pump'),('r-state','r-source','ρ_ss'),('i-frequency','r-source','δ·독립 Ω_SA'),
 ('r-source','r-direct','trace-zero 구동'),('r-direct','r-response','직접 polarization'),
 ('r-response','r-average','속도별 복소 도함수'),('r-geometry','r-average','1D / opt-in 2D 가중치'),
 ('i-drive','r-geometry','T·σ_v'),('i-modes','r-geometry','θ·k'),
 ('d-chi','d-maxwell','physical χ'),('d-geometry','d-maxwell','k·Δk·Q'),
 ('d-maxwell','d-constant','상수 M_cl'),('d-constant','d-amplitude','T_cl'),
 ('d-constant','d-gains','정준 T'),('d-segment','d-amplitude','순서 있는 전파자'),
 ('d-segment','d-gains','정규화한 작은 신호 T'),('d-gains','x-flux','G_p·G_c'),
 ('d-gains','x-canonical','full complex T'),('d-amplitude','f-carrier','복소 출력 carrier'),
 ('m-pair','m-drift','정준 paired basis'),('m-drift','m-k','M_q·J'),
 ('m-reservoir','m-k','J_F·B 교환관계'),('m-reservoir','m-diffusion','B·N_res'),
 ('m-k','m-vacuum','K와 추가 vacuum 선택'),('v-propagator','v-input','T_q(L,0)'),
 ('v-propagator','v-integral','T_q(L,z)'),('v-input','v-output','입력 전파 항'),
 ('v-output','v-quadrature','paired output'),('f-loss','f-measure','V_R,det'),
 ('f-carrier','f-current','검출 진폭·위상'),('f-carrier','f-sql','동일 검출 파워'),
 ('f-current','f-measure','측정 벡터 m'),('f-sql','f-measure','같은 SQL'),
 ('f-measure','o-spectrum','PSD·SQL'),('o-spectrum','o-claim','주장 범위'),
 ('f-tmsv','f-thinning','공유 pair number'),('f-thinning','f-difference','count moments')]: edge(s,t,l)
for s,t,l in [('a-average','a-order','두 차수 응답'),('d-gains','a-order','두 차수 정준 전파'),
 ('a-order','a-gate','full-scan 비교량'),('r-state','r-gauge','gauge map 대상'),
 ('a-average','r-gauge','별도 finite-seed 극한 비교'),('r-direct','r-gauge','참조 응답 비교')]:
    edge(s,t,l,'validation')
for s,t,l in [('r-source','r-poles','대각화 가능시'),('r-poles','r-response','같은 응답의 spectral 표현'),
 ('r-defective','r-response','안정한 응답 평가'),('v-lyapunov','v-phi','좌표무관 행렬함수'),
 ('v-lyapunov','v-eigen','조건 좋은 고유기저'),('v-lyapunov','v-vanloan','권장 블록지수'),
 ('v-phi','v-output','같은 W'),('v-eigen','v-output','같은 W · 조건수 guard'),
 ('v-vanloan','v-output','같은 W'),('v-segments','v-output','z 의존 전체 V')]:edge(s,t,l,'alternative')
for s,t,l in [('r-source','r-defective','결함 행렬 조건'),('r-response','r-voigt','단일 선형 pole·상수 residue'),
 ('d-maxwell','d-segment','z 의존 계수'),('x-flux','x-ideal','이상 가정 추가'),
 ('x-canonical','x-ideal','정준성·입력·검출 조건 추가'),('v-integral','v-lyapunov','상수 M_q·D')]:edge(s,t,l,'limit')
for s,t,l in [('i-levels','m-drift','구동 원자 모델'),('r-average','m-drift','주파수 응답 재료·상태 대응 필요'),
 ('m-drift','v-propagator','미시 M_q 공급'),('m-diffusion','v-integral','미시 D 공급'),
 ('m-drift','v-segments','국소 M_k'),('m-diffusion','v-segments','국소 D_k'),
 ('m-reservoir','v-quadrature','companion reservoir 상관'),('v-quadrature','f-loss','물리적 V_R,out')]:edge(s,t,l,pending=True)
edge('d-constant','v-eigen','trace–discriminant 형태만 유사','analogy')
edge('f-states','f-tmsv','입력 상태와 환경 비교','analogy')
edge('f-thinning','f-loss','순수 손실의 count / covariance 표현','analogy')
edge('i-levels','i-drive','가중 dipole·원자종')
edge('i-frequency','i-modes','carrier 에너지·정규화')
edge('i-drive','m-reservoir','실제 구동 상태와 reservoir 모델',pending=True)
edge('a-probe','m-drift','유한 시드 상태 주위의 요동 선형화 필요',pending=True)
edge('a-probe','m-reservoir','실제 구동 원자 상관의 미시 연결',pending=True)

CSS = r'''
:root{--paper:#f3f4f0;--card:#fff;--ink:#182d36;--muted:#51636c;--line:#d5dddc;--teal:#086b64;--wash:#e9f4ef;--amber:#805010;--amber-bg:#fff5df;--red:#923c32;--blue:#365f8a;--purple:#644e79}
*{box-sizing:border-box}html{scroll-behavior:smooth;scroll-padding-top:115px}body{margin:0;background:var(--paper);color:var(--ink);font:15px/1.65 'Malgun Gothic',system-ui,sans-serif;word-break:keep-all;overflow-wrap:anywhere}a{color:var(--teal);text-underline-offset:4px}button{font:inherit;color:inherit;cursor:pointer}button,a,summary{touch-action:manipulation}a:focus-visible,button:focus-visible,summary:focus-visible{outline:3px solid var(--blue);outline-offset:4px}button{border:1px solid var(--line);background:white;border-radius:7px;padding:8px 13px}h1,h2,h3,h4,p{margin:0}h1{font-size:clamp(32px,4vw,54px);line-height:1.22;letter-spacing:-.05em}h2{font-size:25px;line-height:1.4;letter-spacing:-.025em}h3{font-size:19px}h4{font-size:18px}p+p{margin-top:10px}small{font-size:12px}sub,sup{line-height:0}[hidden]{display:none!important}.shell{max-width:1530px;margin:auto;padding:42px 36px 50px}.eyebrow{font-size:11px;font-weight:700;letter-spacing:.15em;color:var(--teal);margin-bottom:13px}.masthead{display:flex;gap:40px;justify-content:space-between;align-items:end}.intro{max-width:740px;color:var(--muted);margin-top:16px;font-size:16px}.scope-note{max-width:310px;padding:20px;border-left:3px solid #b89145;background:#fffaf0;font-size:13px}.scope-note b{display:block;margin-bottom:6px}.header-meta{display:flex;gap:8px;flex-wrap:wrap;margin-top:23px}.pill,.status{display:inline-block;padding:3px 8px;border-radius:5px;font-size:11px;font-weight:650;line-height:1.7}.pill{background:#e5ebe6;color:#3d555c}.status.ready,.status.input{background:var(--wash);color:var(--teal)}.status.reference,.status.limited,.status.analytic{background:#edf1f7;color:var(--blue)}.status.need{background:var(--amber-bg);color:var(--amber)}.status.test{background:#f0edf4;color:var(--purple)}.status.fail{background:#fbefea;color:var(--red)}.toolbar{position:sticky;top:0;z-index:10;background:#f3f4f0f5;border-bottom:1px solid var(--line);padding:14px 0;margin:23px 0 22px;display:flex;gap:15px;align-items:center;justify-content:space-between}.views{display:flex;gap:4px;background:#e6eae6;padding:4px;border-radius:9px}.views button{border:0;background:transparent;font-size:13px;font-weight:700}.views button[aria-pressed=true]{background:var(--ink);color:white}.tools{display:flex;gap:8px;align-items:center;font-size:12px}.tools>a{margin-right:8px}.js-only{display:none}.enhanced .js-only{display:flex}.noscript{background:var(--amber-bg);padding:16px;border-radius:8px;margin:18px 0}.sheet{background:white;border:1px solid var(--line);border-radius:13px;padding:24px}.section-heading{display:flex;justify-content:space-between;align-items:baseline;gap:18px;margin-bottom:14px}.section-heading p{color:var(--muted);font-size:12px}.overview-layout{display:grid;grid-template-columns:minmax(0,1fr) 285px;gap:24px;align-items:start}.overview{padding:24px 38px}.overview-graph{position:relative;display:grid;grid-template-columns:1fr 1fr;gap:54px 92px;padding:12px 30px}.overview-node{position:relative;z-index:1;display:block;background:white;border:1px solid #aebfc0;border-radius:10px;padding:14px 16px;text-decoration:none;color:var(--ink);min-width:0}.overview-node.wide{grid-column:1 / -1;width:65%;justify-self:center}.overview-node h3{font-size:17px}.overview-node p{font-size:12px;color:var(--muted);margin:5px 0 8px}.overview-node:hover,.overview-node[aria-current=true]{border-color:var(--teal);background:var(--wash);box-shadow:0 0 0 2px #086b6415}.overview-node .number{float:right;color:var(--muted);font-size:11px}.jump{display:block;color:var(--teal);font-size:11px;margin-top:8px}.overview-aside{display:grid;gap:16px}.overview-aside p{font-size:13px;color:var(--muted);margin-top:9px}.overview-aside h3{font-size:16px}.frequency-axis{display:flex;justify-content:space-between;gap:8px;border-bottom:1px solid #9badb1;padding:14px 0 8px;font-family:Cambria,serif;font-size:18px}.legend{display:flex;gap:14px;flex-wrap:wrap;margin-top:18px;font-size:11px;color:var(--muted)}.legend span{display:inline-flex;align-items:center;gap:5px}.sample{width:24px;border-top:2px solid #687f83}.sample.alternative{border-color:var(--blue);border-top-style:dashed}.sample.limit{border-color:var(--purple);border-top-style:dotted}.sample.validation{border-color:var(--amber);border-top-style:dashed}.sample.analogy{border-color:#85818b;border-top-style:dotted}.sample.pending{border-color:var(--amber);border-top-style:double;border-top-width:4px}.svg-wires{position:absolute;inset:0;width:100%;height:100%;pointer-events:none;overflow:visible}.wire{fill:none;stroke:#74898d;stroke-width:1.5}.wire.alternative{stroke:var(--blue);stroke-dasharray:7 4}.wire.limit{stroke:var(--purple);stroke-dasharray:2 4}.wire.validation{stroke:var(--amber);stroke-dasharray:9 3 2 3}.wire.analogy{stroke:#85818b;stroke-dasharray:1 5}.wire.pending{stroke:var(--amber)}.wire.active{stroke:var(--teal);stroke-width:2.5}.wire-label{font:11px 'Malgun Gothic',sans-serif;fill:#435961;paint-order:stroke;stroke:white;stroke-width:6px;stroke-linejoin:round}.overview-links{font-size:12px;margin-top:20px}.route-text{display:grid;gap:7px;padding-left:22px;margin:13px 0;color:var(--muted);font-size:12px}.route-text a{font-weight:600}.route-kind{font-size:10px;color:var(--muted);margin-right:5px}.toc{display:flex;gap:7px;flex-wrap:wrap;margin:15px 0 24px}.toc a{border:1px solid var(--line);background:#fff;border-radius:6px;padding:6px 11px;text-decoration:none;font-size:12px}.toc a[aria-current=true]{background:var(--ink);color:white;border-color:var(--ink)}.detail-heading{margin-top:38px}.detail-layout{display:grid;grid-template-columns:minmax(0,1fr) 340px;gap:22px;align-items:start}.groups{min-width:0}.group{border:1px solid var(--line);border-radius:12px;background:#fff;margin-bottom:14px;scroll-margin-top:110px;overflow:clip}.group>summary{list-style:none;padding:19px 23px;cursor:pointer;display:flex;align-items:center;gap:13px}.group>summary::-webkit-details-marker{display:none}.group>summary:before{content:'+';color:var(--teal);font:22px Cambria,serif}.group[open]>summary:before{content:'−'}.group>summary .summary-title{flex:1;min-width:0}.group>summary strong{font-size:19px;display:block}.group>summary small{display:block;color:var(--muted);margin-top:2px}.group .group-number{font-size:11px;font-weight:700;color:var(--teal)}.group[open]>summary{border-bottom:1px solid var(--line);background:#f7faf7}.group-intro{padding:22px 24px 8px;color:var(--muted);font-size:14px}.part-label{font-size:11px;color:var(--teal);margin-bottom:8px;letter-spacing:.04em}.group-nav{display:flex;flex-wrap:wrap;gap:8px;margin-top:12px;font-size:11px}.group-nav a{border:1px solid var(--line);padding:3px 8px;border-radius:5px;text-decoration:none}.local-graph{position:relative;display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:64px 65px;padding:28px 45px 26px;margin:8px 0;align-items:stretch}.graph-node{position:relative;z-index:1;grid-column:var(--col);grid-row:var(--row);border:1px solid #afc0c2;border-radius:9px;padding:13px 14px;text-decoration:none;color:var(--ink);background:white;min-width:0;align-self:stretch}.graph-node .node-label{display:flex;justify-content:space-between;gap:6px;align-items:center;font-size:10px;color:var(--muted);margin-bottom:7px}.graph-node h4{font-size:15px;line-height:1.5}.graph-node p{font-size:12px;margin-top:6px;line-height:1.65;color:var(--muted)}.graph-node .status{font-size:10px;padding:1px 5px}.graph-node:hover,.graph-node[aria-current=true]{border:1px solid var(--teal);background:var(--wash);box-shadow:0 0 0 2px #086b6420}.graph-node.kind-need{border-color:#b69b64;background:#fffcf4}.graph-node.kind-test{border-left:3px solid #887197}.graph-node.kind-reference{border-left:3px solid #6085a5}.card-math,.card-evidence{display:none}.formula,.card-math{font-family:'Cambria Math',Cambria,serif;word-break:normal;overflow-wrap:anywhere}.card-math{color:var(--ink)!important}.mode-math .card-flow,.mode-evidence .card-flow{display:none}.mode-math .card-math,.mode-evidence .card-evidence{display:block}.mobile-ports{display:none}.route-ledger{border-top:1px solid var(--line);padding:16px 24px;font-size:12px}.route-ledger>summary{cursor:pointer;font-weight:650}.route-ledger ol{columns:2;column-gap:30px}.route-ledger li{break-inside:avoid;margin-bottom:9px}.reference-records{border-top:1px solid var(--line);padding:19px 24px}.reference-records>h3{font-size:14px;margin-bottom:12px;color:var(--muted)}.node-record{border-bottom:1px solid var(--line);padding:10px 0;scroll-margin-top:115px}.node-record:last-child{border-bottom:0}.node-record>summary{cursor:pointer;font-size:13px;display:list-item;list-style-position:inside}.node-record>summary .status{margin-left:8px}.node-content{padding:18px 0 5px;font-size:13px}.node-content>h4{margin-bottom:8px}.node-content .meaning{margin:12px 0 17px;color:var(--muted)}.facts{display:grid;grid-template-columns:72px minmax(0,1fr);gap:9px 12px;margin:14px 0}.facts dt{font-size:11px;font-weight:700;color:var(--muted)}.facts dd{margin:0;min-width:0}.equation-box{padding:14px;background:#f1f5f3;border-radius:7px;white-space:pre-wrap;font-size:16px;line-height:1.9;margin:13px 0;overflow-wrap:anywhere}.condition{border-left:3px solid #b79b68;background:#fffaf0;padding:11px 13px;font-size:12px}.condition b{display:block;color:var(--amber);font-size:11px;margin-bottom:4px}.evidence-text{margin:13px 0;font-size:12px;color:var(--muted)}.source-refs{font-size:11px;color:var(--muted)}.source-refs code{font-family:Consolas,monospace;background:#f2f4f2;padding:1px 3px;overflow-wrap:anywhere}.inspector{position:sticky;top:95px;background:#fff;border:1px solid var(--line);border-radius:12px;max-height:calc(100vh - 112px);overflow:auto;padding:22px;scrollbar-width:thin;scroll-margin-top:115px}.inspector h3{font-size:20px;line-height:1.4}.inspector .node-content{padding-top:5px}.inspector .node-content>h4{display:none}.inspector .facts{grid-template-columns:1fr;gap:4px}.inspector .facts dd{margin-bottom:7px}.inspector .equation-box{font-size:15px}.inspector .group-nav{margin-bottom:8px}.inspector .source-refs{font-size:10px}.inspector .selection-tip{font-size:13px;color:var(--muted);margin-top:12px}.inspector-empty{display:none}.all-expanded .detail-layout{grid-template-columns:1fr}.all-expanded .inspector{display:none}.all-expanded .local-graph{max-width:1050px;margin:auto}.all-expanded .reference-records{display:grid;grid-template-columns:1fr 1fr;gap:14px 28px}.all-expanded .reference-records>h3{grid-column:1/-1}.all-expanded .node-record{align-self:start}.all-expanded .graph-node{max-width:none}.no-wire-note{font-size:11px;color:var(--muted);padding:0 24px 12px}.group-break{grid-column:1/-1;grid-row:4;font-size:12px;color:var(--purple);border-top:1px dashed #b4a9bf;padding-top:16px}.evidence{margin-top:30px}.evidence-header{display:flex;justify-content:space-between;gap:20px;align-items:baseline;margin-bottom:12px}.evidence-header small{color:var(--muted)}.evidence>p{font-size:13px;color:var(--muted);margin:10px 0}.table-scroll{overflow-x:auto}table{width:100%;border-collapse:collapse;text-align:left;font-size:13px}caption{text-align:left;color:var(--muted);font-size:12px;margin-bottom:8px}th,td{padding:12px 14px;border-bottom:1px solid var(--line);vertical-align:top}th{font-size:11px;color:var(--muted)}td p{font-size:12px;color:var(--muted);margin-top:5px}details>summary{cursor:pointer}.evidence details,.math-details{margin-top:18px}.disclosure-content{padding-top:15px;font-size:13px}.disclosure-content p+p{margin-top:10px}.math-details{background:#fff;border:1px solid var(--line);border-radius:12px;padding:22px 24px}.math-details>summary{font-weight:700}.math-grid{display:grid;grid-template-columns:1fr 1fr;gap:24px;margin-top:20px}.math-grid article{border-top:2px solid var(--line);padding-top:14px}.math-grid h3{font-size:16px}.math-grid p{font-size:13px;color:var(--muted);margin-top:10px}.detail-links{display:flex;flex-wrap:wrap;gap:13px;font-size:12px;margin-top:18px}.glossary{margin-top:28px}.glossary summary{font-weight:700}.glossary table{margin-top:18px}.footer{border-top:1px solid var(--line);margin-top:28px;padding-top:18px;display:flex;justify-content:space-between;gap:20px;font-size:12px;color:var(--muted)}.footer nav{display:flex;gap:18px;flex-wrap:wrap}.sr-only{position:absolute;width:1px;height:1px;padding:0;margin:-1px;overflow:hidden;clip:rect(0,0,0,0);white-space:nowrap;border:0}.skip:focus{position:fixed;width:auto;height:auto;clip:auto;background:#fff;z-index:20;padding:12px;left:12px;top:12px}.top-link{font-size:12px;display:inline-block;margin-top:16px}
@media(max-width:1180px){.shell{padding:30px 22px}.overview-layout{grid-template-columns:1fr 235px;gap:16px}.overview{padding:20px}.overview-graph{gap:50px 52px;padding:10px 15px}.detail-layout{grid-template-columns:minmax(0,1fr) 300px;gap:15px}.local-graph{gap:62px 40px;padding:25px 32px}.inspector{padding:18px}.graph-node{padding:12px}.masthead{gap:25px}.scope-note{max-width:265px}}
@media(max-width:960px){.detail-layout{grid-template-columns:1fr}.inspector{position:static;max-height:none}.masthead{display:block}.scope-note{max-width:none;margin-top:20px}.overview-layout{grid-template-columns:1fr}.overview-aside{grid-template-columns:1fr 1fr}.local-graph{gap:64px 75px;padding:25px 50px}.toolbar{top:0}.inspector .facts{grid-template-columns:75px 1fr;gap:8px}.inspector .facts dd{margin-bottom:0}.inspector .source-refs{font-size:11px}}
@media(max-width:600px){html{scroll-padding-top:140px}.shell{padding:23px 14px 35px}.masthead h1{font-size:34px}.intro{font-size:14px}.toolbar{align-items:stretch;flex-direction:column;gap:8px;padding:10px 0;margin-top:18px}.views button{flex:1;padding:8px}.tools{justify-content:space-between}.tools button{padding:6px 9px;font-size:12px}.overview{padding:18px 14px}.overview-graph{display:flex;flex-direction:column;gap:12px;padding:0}.overview-node.wide{width:100%}.overview-node{padding:13px}.overview-graph .svg-wires,.local-graph .svg-wires{display:none}.overview-node .jump{display:none}.overview-aside{grid-template-columns:1fr}.mobile-ports{display:block;border-top:1px solid var(--line);margin-top:9px;padding-top:8px;font-size:10px;color:var(--teal)}.legend{gap:10px;font-size:10px}.section-heading{display:block}.section-heading p{margin-top:6px}.detail-heading{margin-top:28px}.toc{gap:5px}.toc a{font-size:11px;padding:5px 8px}.group>summary{padding:15px 14px;gap:9px}.group>summary strong{font-size:17px}.group>summary>.status{display:none}.group-intro{padding:18px 16px 5px;font-size:13px}.local-graph{display:flex;flex-direction:column;gap:14px;padding:16px}.graph-node{padding:14px}.graph-node h4{font-size:16px}.graph-node p{font-size:13px}.graph-node .node-label{font-size:10px}.group-break{margin-top:8px}.route-ledger,.reference-records{padding:16px}.route-ledger ol{columns:1}.all-expanded .reference-records{display:block}.node-content h4{font-size:17px}.equation-box{font-size:15px;padding:12px}.facts{grid-template-columns:1fr;gap:4px}.facts dd{margin-bottom:8px}.inspector .facts{grid-template-columns:1fr}.math-grid{grid-template-columns:1fr}.sheet,.math-details{padding:19px 16px}.evidence-header{display:block}.evidence-header small{display:block;margin-top:8px}h2{font-size:22px}.table-scroll table{min-width:570px}th,td{padding:10px}.footer{flex-direction:column;gap:12px}.no-wire-note{padding:0 16px 12px}.group,.node-record,.inspector{scroll-margin-top:140px}}
@media(prefers-reduced-motion:reduce){html{scroll-behavior:auto}}
@media print{html{scroll-behavior:auto}body{background:#fff;font-size:10pt}.shell{max-width:none;padding:0}.toolbar,.inspector,.toc,.js-only,.top-link,.jump,.overview-aside{display:none!important}.masthead{display:block}.scope-note{max-width:none}.overview-layout,.detail-layout{display:block}.overview-graph{gap:15px 30px;padding:0}.overview-node{break-inside:avoid}.overview .svg-wires,.local-graph,.no-wire-note,.legend{display:none}.overview-links{display:block}.group{overflow:visible;break-inside:auto;border:0;border-top:2px solid #888;margin-top:24px}.group>summary{background:#fff!important;padding:12px 0}.group-intro,.reference-records,.route-ledger{padding:12px 0}.group>summary:before{content:''}.node-record{break-inside:avoid;padding:12px 0}.node-content{padding-top:8px}.all-expanded .reference-records{display:block}.equation-box{white-space:pre-wrap;font-size:11pt}.math-grid{display:block}.math-grid article{margin-bottom:18px;break-inside:avoid}details> :not(summary){display:block!important}.facts{display:grid!important}.table-scroll{overflow:visible}.table-scroll table{min-width:0}.source-refs{font-size:8pt}.footer{font-size:8pt}.status{print-color-adjust:exact}a{color:inherit;text-decoration:none}.route-ledger ol{columns:1}.sheet{border:0;padding:12px 0}.evidence table{font-size:9pt}tr{break-inside:avoid}}
'''

SCRIPT = r'''
"use strict";
const model=JSON.parse(document.getElementById('atlas-data').textContent);
const groups=model.groups, edges=model.edges, nodeIndex=new Map(groups.flatMap(g=>g.nodes.map(n=>[n.id,n])));
const groupEls=new Map(groups.map(g=>[g.id,document.getElementById('group-'+g.id)]));
const inspector=document.getElementById('inspector');
let selectedGroup='atomic',selectedNode='a-assembly',allMode=false,printState=null,drawFrame=0;
const motion=()=>matchMedia('(prefers-reduced-motion:reduce)').matches?'instant':'smooth';
function scheduleDraw(){cancelAnimationFrame(drawFrame);drawFrame=requestAnimationFrame(drawAll);}
function setCurrent(){
 document.querySelectorAll('[data-group-link]').forEach(a=>a.setAttribute('aria-current',String(a.dataset.groupLink===selectedGroup)));
 document.querySelectorAll('.graph-node').forEach(a=>a.setAttribute('aria-current',String(a.dataset.node===selectedNode)));
}
function inspectNode(id){
 const n=nodeIndex.get(id);if(!n)return;
 selectedNode=id;selectedGroup=n.group;
 const source=document.getElementById(id).querySelector('.node-content');
 const clone=source.cloneNode(true);clone.querySelectorAll('[id]').forEach(x=>x.removeAttribute('id'));
 inspector.replaceChildren();
 const eyebrow=document.createElement('p');eyebrow.className='eyebrow';eyebrow.textContent=groups.find(g=>g.id===n.group).title+' / '+id;
 const title=document.createElement('h3');title.textContent=n.title;title.tabIndex=-1;
 inspector.append(eyebrow,title,clone);
 const back=document.createElement('a');back.href='#group-'+n.group;back.textContent='↑ 내부 도식으로 돌아가기';back.className='top-link';inspector.append(back);
 inspector.scrollTop=0;setCurrent();scheduleDraw();
 document.getElementById('announcement').textContent=n.title+' — '+n.out;
}
function activateGroup(id,scroll=false){
 if(!groupEls.has(id))return;
 selectedGroup=id;
 if(!allMode)groupEls.forEach((el,key)=>{el.open=key===id;});else groupEls.get(id).open=true;
 if(nodeIndex.get(selectedNode)?.group!==id)selectedNode=groups.find(g=>g.id===id).nodes[0].id;
 inspectNode(selectedNode);
 if(scroll)groupEls.get(id).scrollIntoView({block:'start',behavior:motion()});
}
function showNode(id,scroll=false){
 const n=nodeIndex.get(id);if(!n)return;
 activateGroup(n.group,false);inspectNode(id);
 if(allMode){const record=document.getElementById(id);record.open=true;if(scroll)record.scrollIntoView({block:'start',behavior:motion()});}
 else if(scroll&&matchMedia('(max-width:960px)').matches){inspector.scrollIntoView({block:'start',behavior:motion()});inspector.querySelector('h3').focus({preventScroll:true});}
 else if(scroll){const card=document.querySelector('.graph-node[data-node="'+id+'"]');card?.scrollIntoView({block:'nearest',behavior:motion()});}
}
function applyHash(scroll=false){
 if(document.body.classList.contains('explorer-active'))return;
 let id;try{id=decodeURIComponent(location.hash.slice(1));}catch{return;}
 if(id.startsWith('group-'))activateGroup(id.slice(6),scroll);
 else if(nodeIndex.has(id))showNode(id,scroll);
 else if(id){const el=document.getElementById(id);if(el){for(let p=el;p;p=p.parentElement)if(p.tagName==='DETAILS')p.open=true;if(scroll)el.scrollIntoView({block:'start',behavior:motion()});}}
}
function changeHash(id){if(location.hash!=='#'+id)history.pushState(null,'','#'+id);}
function expandAll(value){
 allMode=value;document.body.classList.toggle('all-expanded',value);
 document.getElementById('expand-all').setAttribute('aria-pressed',String(value));
 document.getElementById('expand-all').textContent=value?'선택 그룹으로 접기':'전체 상세 펼치기';
 groupEls.forEach((el,id)=>{el.open=value||id===selectedGroup;});
 document.querySelectorAll('.node-record,.route-ledger,.node-connections,.glossary,.math-details,.evidence details').forEach(el=>{el.open=value;});
 if(!value)groupEls.get(selectedGroup).scrollIntoView({block:'start',behavior:motion()});
 scheduleDraw();
}
const NS='http://www.w3.org/2000/svg';
function svgTag(name,attrs){const el=document.createElementNS(NS,name);for(const [k,v] of Object.entries(attrs))el.setAttribute(k,v);return el;}
function drawGraph(container,links,selector,prefix){
 const svg=container.querySelector('.svg-wires');if(!svg)return;svg.replaceChildren();
 if(matchMedia('(max-width:600px)').matches||container.getBoundingClientRect().height===0)return;
 const rect=container.getBoundingClientRect();svg.setAttribute('viewBox',`0 0 ${rect.width} ${rect.height}`);
 const defs=svgTag('defs',{}),marker=svgTag('marker',{id:prefix+'-arrow',viewBox:'0 0 10 10',refX:9,refY:5,markerWidth:6,markerHeight:6,orient:'auto'});
 marker.append(svgTag('path',{d:'M1 1 L9 5 L1 9',fill:'none',stroke:'#74898d','stroke-width':1.6}));defs.append(marker);svg.append(defs);
 const cards=new Map([...container.querySelectorAll(selector)].map(el=>[el.dataset.node||el.dataset.groupLink,el]));
 const pos=el=>{const r=el.getBoundingClientRect();return{x:r.x-rect.x,y:r.y-rect.y,w:r.width,h:r.height};};
 links.forEach((e,i)=>{
  if(!cards.has(e.source)||!cards.has(e.target))return;
  const a=pos(cards.get(e.source)),b=pos(cards.get(e.target));let x1=a.x+a.w/2,y1=a.y+a.h+2,x2=b.x+b.w/2,y2=b.y-4;
  let path,lx,ly;const down=y2-y1;
  if(down>115||down<0){
   const lane=12+(i%3)*4;const useRight=(a.x+b.x)>rect.width;
   const side=useRight?rect.width-lane:lane;
   x1=useRight?a.x+a.w+2:a.x-2;y1=a.y+a.h*.62;x2=useRight?b.x+b.w+3:b.x-3;y2=b.y+b.h*.35;
   path=`M${x1},${y1} H${side} V${y2} H${x2}`;lx=side;ly=(y1+y2)/2;
  }else{const mid=(y1+y2)/2;path=`M${x1},${y1} V${mid} H${x2} V${y2}`;lx=(x1+x2)/2;ly=mid-5;}
  const active=e.source===selectedNode||e.target===selectedNode;
  const p=svgTag('path',{d:path,class:`wire ${e.kind||'data'}${e.pending?' pending':''}${active?' active':''}`,'marker-end':`url(#${prefix}-arrow)`});
  const title=svgTag('title',{});title.textContent=e.label;p.append(title);svg.append(p);
  const label=svgTag('text',{x:lx,y:ly,class:'wire-label','text-anchor':'middle'});
  label.textContent=prefix==='overview'?e.label:String(i+1).padStart(2,'0')+(e.pending?' ◇':'');svg.append(label);
 });
}
function drawAll(){
 drawGraph(document.getElementById('overview-graph'),model.spine,'.overview-node','overview');
 groups.forEach(g=>{const el=groupEls.get(g.id);if(!el.open)return;
  const local=edges.filter(e=>nodeIndex.get(e.source).group===g.id&&nodeIndex.get(e.target).group===g.id);
  drawGraph(el.querySelector('.local-graph'),local,'.graph-node',g.id);
 });
}
document.documentElement.classList.add('enhanced');
document.querySelectorAll('.node-record,.route-ledger,.node-connections,.glossary,.overview-links').forEach(el=>{el.open=false;});
document.querySelectorAll('.group>summary').forEach(summary=>summary.addEventListener('click',event=>{
 if(document.body.classList.contains('explorer-active'))return;
 event.preventDefault();const id=summary.parentElement.dataset.group;
 if(allMode){summary.parentElement.open=!summary.parentElement.open;scheduleDraw();return;}
 changeHash('group-'+id);activateGroup(id,false);
}));
document.addEventListener('click',event=>{
 if(document.body.classList.contains('explorer-active'))return;
 const a=event.target.closest('a[href^="#"]');if(!a)return;
 const id=a.getAttribute('href').slice(1);
 if(nodeIndex.has(id)){event.preventDefault();changeHash(id);showNode(id,true);}
 else if(id.startsWith('group-')&&groupEls.has(id.slice(6))){event.preventDefault();changeHash(id);activateGroup(id.slice(6),true);}
});
document.querySelectorAll('.node-record').forEach(record=>record.addEventListener('toggle',()=>{
 if(document.body.classList.contains('explorer-active'))return;
 if(record.open&&!allMode&&!printState){inspectNode(record.id);groupEls.get(nodeIndex.get(record.id).group).open=true;}
}));
document.querySelectorAll('[data-view]').forEach(button=>button.addEventListener('click',()=>{
 document.body.classList.remove('mode-math','mode-evidence');
 if(button.dataset.view!=='flow')document.body.classList.add('mode-'+button.dataset.view);
 document.querySelectorAll('[data-view]').forEach(b=>b.setAttribute('aria-pressed',String(b===button)));
 scheduleDraw();
}));
document.getElementById('expand-all').addEventListener('click',()=>expandAll(!allMode));
document.getElementById('print-page').addEventListener('click',()=>window.print());
window.addEventListener('hashchange',()=>applyHash(true));window.addEventListener('popstate',()=>applyHash(true));
window.addEventListener('beforeprint',()=>{if(printState)return;printState=[...document.querySelectorAll('details')].map(el=>[el,el.open]);printState.forEach(([el])=>{el.open=true;});});
window.addEventListener('afterprint',()=>{const saved=printState;if(!saved)return;saved.forEach(([el,open])=>{el.open=open;});setTimeout(()=>{printState=null;scheduleDraw();},0);});
const resizeObserver=new ResizeObserver(scheduleDraw);
document.querySelectorAll('.local-graph,#overview-graph').forEach(el=>resizeObserver.observe(el));
activateGroup(selectedGroup);applyHash(false);document.fonts.ready.then(scheduleDraw);
'''

SPINE = [dict(source=s,target=t,label=l,kind='data',pending=p) for s,t,l,p in [
 ('input','atomic','운전점',False),('atomic','mean','복소 응답',False),
 ('input','micro','구동 원자 모델',True),('mean','gain','두 이득',False),
 ('micro','covariance','M_q · D',True),('mean','detector','복소 평균장',False),
 ('covariance','detector','4×4 공분산',True),('detector','output','측정 · SQL',True)]]

GLOSSARY = [
 ('pump / p / c','강한 pump / probe(시드가 입사하는 모드) / conjugate. p는 pump의 약자가 아니다.'),
 ('Δ / δ','pump의 one-photon detuning / Raman two-photon detuning. 식에서는 rad/s.'),
 ('Ω_beat / Ω_SA','pump–probe optical beat / 독립 RF 분석 주파수. Floquet 조화와 측정 sideband의 역할이 다르다.'),
 ('ω_k / λ_m','정적 pump-frame resolvent 주파수 / 원자 Liouvillian 고유값. pole 중심은 ω_k 축에서 −Imλ_m.'),
 ('ρ_n / ρ_ss','생산 finite-seed Floquet harmonic / 별도 pump-only 정상상태. 같은 상태를 부르는 두 이름이 아니다.'),
 ('χ̄ / ⟨χ̄⟩_v / χ','시간 차원 reduced response / 그것의 속도 평균 / 거시 무차원 physical susceptibility.'),
 ('ℒ / M_cl / M_q','원자 density matrix의 Liouvillian(s⁻¹) / 평균장 drift(m⁻¹) / 주파수 의존 요동 drift(m⁻¹).'),
 ('T_cl / T_q','평균장 / 양자 요동 전파자. 둘 다 무차원 행렬이지만 일반적으로 서로 다르다.'),
 ('Q / 𝒢','정준 photon-flux 진폭 → 전기장 / 정준 진폭 → Rabi 진폭 변환. 모드 면적과 carrier 에너지가 필요하다.'),
 ('J / J_F / K','paired 광학 commutator metric / reservoir metric / drift가 요구하는 commutator-defect matrix.'),
 ('B / N_res / D','reservoir coupling / symmetrized reservoir 상관 / 광학 diffusion D=BN_resB†. D는 원자 산일자와 다른 양이다.'),
 ('T V_in T† / W','입력 요동의 전파 / 셀 안에서 분포 생성되어 전파된 noise. 두 covariance 기여를 합산한다.'),
 ('V_Ξ / V_R','creation/annihilation sideband 순서의 4×4 covariance / (X_p,Y_p,X_c,Y_c) quadrature covariance.'),
 ('φ₁ / w_F / w','entire Lyapunov 함수 / Faddeeva 함수 / 검출 전자 가중치. 서로 다른 함수다.'),
 ('G_p / G_c / g_c^(Φ)','동일 probe 입력 P₀에 대한 두 power gain / conjugate photon-flux gain. G_c=(ω_c/ω_p)g_c^(Φ).'),
 ('η / H / S_el / SQL','셀 이후 효율 / 전자 주파수 응답 / 전자 PSD / 동일 측정 조건의 shot-noise 기준.'),
 ('TMSV / displaced thermal','진공 입력 이상 두 모드 squeezed vacuum / coherent seed가 있는 이상 출력의 한 팔 상태.'),
 ('NRF / S_−(Ω)','정의한 모드의 photon-count 차동 정규화 / RF 분석 주파수별 검출 PSD 정규화. 입력·측정 조건을 함께 봐야 한다.'),
 ('Ω_R / Ω_SA','real-quadrature symplectic 행렬 / RF 주파수. Gaussian complete positivity 식에서 같은 기호로 읽지 않는다.')]

LEVEL_DIAGRAM = '''<figure class="level-figure"><svg viewBox="0 0 600 335" role="img" aria-label="대표 double-Lambda 회로: pump가 두 ground manifold를 구동하고 probe는 g2, conjugate는 g1에서 excited manifold를 연결합니다.">
<g fill="none" stroke="#51636c" stroke-width="3"><path d="M35 90H180 M400 50H560 M35 270H180 M400 240H560"/></g>
<g fill="#182d36" font-family="Cambria,serif" font-size="25"><text x="35" y="72">e₂ (F′=2)</text><text x="400" y="32">e₃ (F′=3)</text><text x="35" y="307">g₁ (F=2)</text><text x="400" y="278">g₂ (F=3)</text></g>
<g fill="none" stroke="#086b64" stroke-width="3"><path d="M130 265L430 58 M135 251L130 265L146 263 M416 60L430 58L424 73"/><path d="M430 235L130 95 M144 96L130 95L138 108 M417 235L430 235L424 223"/></g>
<g fill="none" stroke="#365f8a" stroke-width="3"><path d="M490 235V55 M483 69L490 55L497 69"/></g>
<g fill="none" stroke="#805010" stroke-width="3"><path d="M78 265V95 M71 109L78 95L85 109"/></g>
<g font-family="Cambria,serif" font-size="25"><text x="244" y="130" fill="#086b64">pump</text><text x="250" y="235" fill="#086b64">pump</text><text x="500" y="170" fill="#365f8a">probe</text><text x="10" y="190" fill="#805010">c</text></g>
</svg><figcaption>minus 모드쌍의 대표 Λ 회로 · 에너지 간격은 비례 축이 아닙니다. c는 conjugate입니다. 실제 축약 계산은 각 ground에서 두 excited manifold로의 가중 결합·readout을 모두 포함합니다.</figcaption></figure>'''

def build():
    import re
    from fwm_quotient_structure_v4_science import extend, INTRO, SOURCE_CATALOG
    extend(GROUPS, EDGES, KINDS, EDGE_KINDS, NUMBERS, CODE)
    nodes={n['id']:n for g in GROUPS for n in g['nodes']}
    gs={g['id']:g for g in GROUPS}
    assert len(nodes)==sum(len(g['nodes']) for g in GROUPS)
    assert all(e['source'] in nodes and e['target'] in nodes for e in EDGES)
    tex_path=ROOT/'squeezing_analytic_reconstruction_v3.tex'
    tex=tex_path.read_text(encoding='utf-8') if tex_path.exists() else ''
    sections=json.loads((ROOT/'fwm_quotient_structure_v4_sources.json').read_text(encoding='utf-8'))
    current='Definitions and conventions'
    for line in tex.splitlines():
        m=re.match(r'\\(?:(?:sub)*section|part|chapter)\*?\{(.*)\}',line)
        if m:current=m.group(1)
        for label in re.findall(r'\\label\{([^}]+)\}',line):sections[label]=current
    sections.update({k:v[0] for k,v in SOURCE_CATALOG.items()})
    used=sorted({r for n in nodes.values() for r in n['refs']})
    assert set(used)<=sections.keys(),set(used)-sections.keys()
    if tex:
        actual_labels=set(re.findall(r'\\label\{([^}]+)\}',tex))
        missing_labels=set(used)-SOURCE_CATALOG.keys()-actual_labels
        assert not missing_labels, ('Referenced labels absent from current theory source', missing_labels)
    h=escape
    def badge(kind):return '<span class="status '+kind+'">'+h(KINDS[kind])+'</span>'
    def refid(key):return 'source-'+key.replace(':','-')
    def link(id):return '<a href="#'+id+'">'+h(nodes[id]['title'])+'</a>'
    def route(e,number=None):
        mark='<span class="route-kind">'+h(EDGE_KINDS[e['kind']])+(' · ◇ 미완성' if e['pending'] else '')+'</span>'
        return '<li data-source="'+e['source']+'" data-target="'+e['target']+'">'+mark+link(e['source'])+' → '+link(e['target'])+'<br>'+h(e['label'])+'</li>'
    def groupnav(g):
        ins=sorted({nodes[e['source']]['group'] for e in EDGES if nodes[e['target']]['group']==g['id'] and nodes[e['source']]['group']!=g['id']})
        outs=sorted({nodes[e['target']]['group'] for e in EDGES if nodes[e['source']]['group']==g['id'] and nodes[e['target']]['group']!=g['id']})
        return '<div class="group-nav">'+''.join('<a href="#group-'+x+'">'+lead+h(gs[x]['title'])+'</a>' for lead,values in [('입력·비교 ← ',ins),('출력·비교 → ',outs)] for x in values)+'</div>'
    def noderecord(n):
        incoming=[e for e in EDGES if e['target']==n['id']]
        outgoing=[e for e in EDGES if e['source']==n['id']]
        refs=' · '.join('<a href="#'+refid(r)+'"><code>'+h(r)+'</code></a>' for r in n['refs'])
        return f'''<details class="node-record" id="{n['id']}" open>
<summary><strong>{h(n['title'])}</strong> {badge(n['kind'])}</summary>
<div class="node-content"><h4>{h(n['title'])}</h4>{badge(n['kind'])}
<p class="meaning">{h(n['body'])}</p>
<dl class="facts"><dt>입력</dt><dd>{h(n['inputs'])}</dd><dt>출력</dt><dd>{h(n['out'])}</dd><dt>차원·기저</dt><dd>{h(n['units'])}</dd></dl>
<div class="equation-box formula">{h(n['math'])}</div>
<p class="condition"><b>성립 조건과 경계</b>{h(n['condition'])}</p>
<p class="evidence-text"><b>검증 범위</b><br>{h(n['evidence'])}</p>
<p class="source-refs">원문 {h(gs[n['group']]['part'])}<br>{refs}</p>
<details open class="node-connections"><summary>입력 출처와 출력 목적지</summary>
<ul class="route-text">{''.join(route(e) for e in incoming+outgoing) or '<li>독립 외부 입력 / 명시된 상태 시험</li>'}</ul></details>
</div></details>'''
    def graphcard(n):
        incoming=[e for e in EDGES if e['target']==n['id']]
        outgoing=[e for e in EDGES if e['source']==n['id']]
        ports='입력 ← '+(' · '.join(nodes[e['source']]['title'] for e in incoming) or '명시된 외부 입력')+'<br>출력 → '+(' · '.join(nodes[e['target']]['title'] for e in outgoing) or '이 그룹의 출력·조건')
        return f'''<a class="graph-node kind-{n['kind']}" data-node="{n['id']}" href="#{n['id']}" style="--col:{n['xy'][0]};--row:{n['xy'][1]}">
<span class="node-label">{h(n['id'])}{badge(n['kind'])}</span><h4>{h(n['title'])}</h4>
<p class="card-flow">{h(n['out'])}</p><p class="card-math">{h(n['compact'])}</p><p class="card-evidence">{h(n['evidence'])}</p>
<span class="mobile-ports">{ports}</span></a>'''
    grouphtml=[]
    for g in GROUPS:
        local=[e for e in EDGES if nodes[e['source']]['group']==g['id'] and nodes[e['target']]['group']==g['id']]
        cross=[e for e in EDGES if (nodes[e['source']]['group']==g['id']) != (nodes[e['target']]['group']==g['id'])]
        codepath,codewords=CODE[g['id']]
        grouphtml.append(f'''<details class="group" data-group="{g['id']}" id="group-{g['id']}" open>
<summary><span class="group-number">{NUMBERS[g['id']]}</span><span class="summary-title"><strong>{h(g['title'])}</strong><small>{h(g['subtitle'])} · 내부 {len(g['nodes'])}개 노드</small></span>{badge(g['kind'])}</summary>
<div class="group-intro"><p class="part-label">{h(g['part'])}</p><p>{h(g['brief'])}</p>{groupnav(g)}{LEVEL_DIAGRAM if g['id']=='input' else ''}</div>
<div class="local-graph" aria-label="{h(g['title'])} 내부 연결도"><svg class="svg-wires" aria-hidden="true"></svg>
{''.join(graphcard(n) for n in g['nodes'])}
{'<p class="group-break">별도 이상 시험 · TMSV와 손실 통계 / 위 검출 경로의 미시 covariance를 대체하지 않습니다.</p>' if g['id']=='detector' else ''}</div>
<p class="no-wire-note">선의 번호는 아래 연결 목록과 대응합니다. 계산·대안·축약·검증은 서로 다른 연결입니다.</p>
<details class="route-ledger" open><summary>내부 연결 {len(local)}개 · 입력·출력 포트 {len(cross)}개</summary>
<h4>내부 연결 / 선 번호</h4><ol class="route-text local-routes">{''.join(route(e,i+1) for i,e in enumerate(local))}</ol>
<h4>다른 그룹과의 연결</h4><ul class="route-text">{''.join(route(e) for e in cross)}</ul></details>
<section class="reference-records"><h3>모든 노드의 수식·조건·근거</h3>{''.join(noderecord(n) for n in g['nodes'])}</section>
<div class="route-ledger source-refs">현재 구현·검사 대응: <a href="../../{codepath}">{codepath}</a><br>{h(codewords)}<br>구현 링크는 추가 추적용입니다. 이 HTML의 설명은 파일 하나에 모두 포함되어 있습니다.</div>
</details>''')
    overview=[]
    for id in ['input','atomic','mean','micro','gain','covariance','detector','output']:
        g=gs[id];wide=id in ['input','atomic','detector','output']
        parents=[e for e in SPINE if e['target']==id]
        overview.append(f'''<a href="#group-{id}" class="overview-node{' wide' if wide else ''}" data-group-link="{id}">
<span class="number">{NUMBERS[id]}</span><h3>{g['title']}</h3><p>{g['subtitle']}</p>{badge(g['kind'])}
<span class="jump">내부 구조 {len(g['nodes'])}개 노드 ↗</span><span class="mobile-ports">입력 ← {h(' + '.join(gs[e['source']]['title'] for e in parents) or '실험 파라미터')}<br>출력 → {h(' + '.join(gs[e['target']]['title'] for e in SPINE if e['source']==id) or '조건이 명시된 RF 스펙트럼')}</span></a>''')
    legend='<div class="legend">'+''.join('<span><i class="sample '+k+'"></i>'+v+'</span>' for k,v in EDGE_KINDS.items())+'<span><i class="sample pending"></i>◇ 미완성 연결</span></div>'
    # Retain the original evidence and mathematical explanations; add node targets.
    from fwm_quotient_structure_v4_science import update_context
    context=update_context(PRESERVED_CONTEXT)
    for title,target in [('참조 응답의 방정식 잔차','r-direct'),('pump 상태의 프레임 비교','r-gauge'),
      ('유한 시드 → 무한소 응답','r-gauge'),('Floquet 절단 수렴','a-gate'),('2D Raman–Doppler','r-geometry'),
      ('실험 이득 비교','d-gains'),('물리적 스퀴징 스펙트럼','o-claim')]:
        context=context.replace('<td>'+title+'</td>','<td><a href="#'+target+'">'+title+'</a></td>')
    context=context.replace('과거 N<sub>F</sub>=1의 16.8% 변화는 모든 차수의 실패를 뜻하지 않습니다.','명시된 fixture에서 N<sub>F</sub>=1→2는 G<sub>p</sub>를 16.8% 바꾸고 계수 갭 부호를 뒤집습니다. trace/positivity 검사만으로 절단 수렴을 판정하지 않습니다.')
    catalog=''.join('<tr id="'+refid(r)+'"><td><code>'+h(r)+'</code></td><td>'+h(sections[r])+(' · <a href="'+h(SOURCE_CATALOG[r][1])+'">증거와 유도</a>' if r in SOURCE_CATALOG else '')+'</td></tr>' for r in used)
    toc=''.join('<a data-group-link="'+g['id']+'" href="#group-'+g['id']+'">'+NUMBERS[g['id']]+' '+h(g['title'])+'</a>' for g in GROUPS)
    from fwm_quotient_structure_v4_science import EXTRA_SPINE
    model=dict(version='v4',groups=GROUPS,edges=EDGES,spine=SPINE + EXTRA_SPINE,
               source='squeezing_analytic_reconstruction_v3.tex',code=CODE)
    payload=json.dumps(model,ensure_ascii=False,separators=(',',':')).replace('<','\\u003c')
    sourcehash=__import__('hashlib').sha256(tex_path.read_bytes()).hexdigest() if tex else 'catalog-only; scientific content embedded'
    page=f'''<!doctype html>
<html lang="ko"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<meta name="color-scheme" content="light"><meta name="description" content="GABES FWM: 유한 시드 원자 응답에서 평균장, 미시 요동, 공분산과 검출 스펙트럼까지의 해석 지도.">
<title>FWM · 원자 응답에서 검출 잡음까지 — 논리와 물리</title>
<!-- Generated by fwm_quotient_structure_v4_builder.py. All content and graph data are embedded. -->
<!-- Source TeX SHA256: {sourcehash} -->
<style>{CSS}
.level-figure{{max-width:620px;margin:20px auto 8px;background:#f7f9f6;padding:12px;border-radius:9px}}.level-figure svg{{width:100%;height:auto;display:block}}.level-figure figcaption{{font-size:11px;color:var(--muted);margin-top:8px}}@media print{{.level-figure{{max-width:450px;break-inside:avoid}}}}
html:not(.enhanced) .detail-layout{{display:block}}html:not(.enhanced) .inspector{{display:none}}html:not(.enhanced) .svg-wires{{display:none}}
</style></head><body>
<a href="#atlas" class="sr-only skip">전체 물리 흐름으로 건너뛰기</a>
<main class="shell"><header class="masthead"><div><p class="eyebrow">GABES / FWM ANALYTIC ATLAS · LOGIC & PHYSICS</p>
<h1>FWM 이론의 논리와<br>검출 예측의 조건</h1>
<p class="intro">가정, 정확한 동치, 정보의 손실, 수치 근사와 실험 검증을 구별합니다.<br>각 결론이 필요로 하는 입력·수식·증거를 같은 지도에서 읽습니다.</p>
<div class="header-meta"><span class="pill">⁸⁵Rb D1 · double-Λ</span><span class="pill">{len(GROUPS)}개 그룹 · {len(nodes)}개 내부 노드</span><span class="pill">단일 HTML · 오프라인</span></div></div>
<aside class="scope-note"><b>조건부 미시 이론과 장치 예측을 구분</b>명시한 GKSL 원자·Gaussian 광장에서는 drift, noise, 검출 연결이 검산되어 있습니다. 움직이는 열원자 전체, nonlocal Maxwell, 독립 실측 입력의 연결은 미완료입니다.</aside></header>
<nav class="toolbar" aria-label="지도 탐색"><div class="views js-only" role="group" aria-label="정보 보기">
<button data-view="flow" aria-pressed="true">물리 흐름</button><button data-view="math" aria-pressed="false">수식과 가정</button><button data-view="evidence" aria-pressed="false">검증 근거</button></div>
<div class="tools"><a href="#atlas">전체 지도</a><a href="#evidence">검증 표</a><span class="js-only"><button id="expand-all" aria-pressed="false">전체 상세 펼치기</button></span><span class="js-only"><button id="print-page">인쇄</button></span></div></nav>
<noscript><p class="noscript">모든 설명과 연결 목록이 펼쳐져 있습니다. JavaScript 없이도 노드·수식·조건을 읽을 수 있습니다. 각 카드의 링크는 상세 본문으로 이동합니다.</p></noscript>
{INTRO}<section id="atlas"><div class="section-heading"><h2>전체 물리 흐름</h2><p>화살표는 의존관계 · 카드를 선택하면 그룹 내부로 이동</p></div>
<div class="overview-layout"><div class="sheet overview"><div id="overview-graph" class="overview-graph"><svg class="svg-wires" aria-hidden="true"></svg>{''.join(overview)}</div>{legend}
<details class="overview-links" open><summary>전체 연결을 텍스트로 읽기</summary><ol class="route-text">{''.join('<li><a href="#group-'+e['source']+'">'+gs[e['source']]['title']+'</a> → <a href="#group-'+e['target']+'">'+gs[e['target']]['title']+'</a> · '+e['label']+(' / ◇ 미완성' if e['pending'] else '')+'</li>' for e in SPINE + EXTRA_SPINE)}</ol></details></div>
<aside class="overview-aside"><section class="sheet"><p class="eyebrow">REFERENCE LANE</p><h3>pump-only 약응답 참조</h3>{badge('reference')}<p>정적 pump 상태 → trace-zero Nambu 응답 → 1D / opt-in 2D 평균</p><p>원자 조립을 공유하는 별도 계산입니다. minus 한정이며 미시 diffusion을 공급하지 않습니다.</p><p><a href="#group-reference" data-group-link="reference">R · 내부 구조 {len(gs['reference']['nodes'])}개 노드 ↗</a></p><p>응답 재료 <a href="#m-drift">→ 미시 drift 구성</a><br>비교 근거 <a href="#r-gauge">→ 프레임·seed 극한</a></p></section>
<section class="sheet"><h3>두 주파수의 역할</h3><div class="frequency-axis"><span>ω_p</span><span>ω_pump</span><span>ω_c</span></div><p>광학 carrier의 개념도 / 비례 축 아님</p><p class="formula">Ω_beat=ω_pump−ω_p</p><p class="formula">δa_p(+Ω_SA) ↔ δa_c†(−Ω_SA)</p><p>RF sideband를 companion block과 조립합니다. 네 성분은 네 독립 carrier가 아닙니다.</p><a class="top-link" href="#i-frequency">주파수 규약 읽기 →</a></section>
<section class="sheet"><h3>현재 출력의 의미</h3><p>GABES의 Squeezing indicator는 gain 기반 평균장 지표입니다. 여기의 최종 RF 스펙트럼과 같은 출력으로 해석하지 않습니다.</p><p><a href="#o-claim">구현 범위와 필요한 입력 →</a></p></section></aside></div></section>
<section aria-labelledby="detail-title"><div class="section-heading detail-heading"><h2 id="detail-title">그룹 내부 구조</h2><p>카드 → 노드 상세 · 원문 없이도 모든 내용을 읽을 수 있습니다.</p></div><nav class="toc" aria-label="그룹 목차">{toc}</nav>
<div class="detail-layout"><div class="groups">{''.join(grouphtml)}</div><aside class="inspector" id="inspector" aria-label="선택한 노드 상세" tabindex="-1"><p class="eyebrow">NODE INSPECTOR</p><h3>내부 노드를 선택하세요</h3><p class="selection-tip">입력과 출력, 수식, 성립 조건, 검증 범위를 이곳에서 함께 읽습니다.</p></aside></div></section>
<details class="sheet glossary" id="glossary" open><summary>기호·용어·단위 사전</summary><div class="table-scroll"><table><caption>같은 문자라도 물리적 객체와 사용 기저를 구별합니다.</caption><thead><tr><th>기호</th><th>뜻과 읽는 법</th></tr></thead><tbody>{''.join('<tr><td>'+h(k)+'</td><td>'+h(v)+'</td></tr>' for k,v in GLOSSARY)}</tbody></table></div></details>
{context}
<details class="sheet glossary" id="sources" open><summary>원문 근거 찾아보기 · {len(used)}개 검색 키</summary><p class="evidence-text">수식 키는 이론 TeX의 section·equation·table을, logic:·gc: 키는 이 문서의 정의와 별도 유도·증거를 가리킵니다. 핵심 식과 조건은 위 노드 본문에 내장되어 있습니다. <a href="squeezing_analytic_reconstruction_v3.tex">원문 TeX</a>와 <a href="squeezing_analytic_reconstruction_v3.pdf">PDF</a>는 추가 검증용입니다.</p><div class="table-scroll"><table><thead><tr><th>원문 검색 키</th><th>원문 section</th></tr></thead><tbody>{catalog}</tbody></table></div></details>
<footer class="footer"><p>FWM 이론의 순서·동치·조건부 예측<br>근거 범위: 2026-09-18. 배지는 실시간 실행 상태가 아닙니다.</p><nav><a href="#group-structure">순서와 동치</a><a href="squeezing_analytic_reconstruction_v3.tex">이론 TeX</a><a href="squeezing_analytic_reconstruction_v3.pdf">이론 PDF</a><a href="../squeezing_report/squeezing_report_v8.pdf">검증 보고서</a></nav></footer>
<p id="announcement" class="sr-only" aria-live="polite" aria-atomic="true"></p></main>
<script type="application/json" id="atlas-data">{payload}</script>
<script>{SCRIPT}</script></body></html>'''
    # Author the explorer separately, but deliver a single offline HTML file.
    view_css=(ROOT/'fwm_quotient_structure_v4_view.css').read_text(encoding='utf-8')
    view_js=(ROOT/'fwm_quotient_structure_v4_view.js').read_text(encoding='utf-8')
    page=page.replace('</head>', '<style>'+view_css+'</style></head>')
    page=page.replace('</body>', '<script>'+view_js+'</script></body>')
    # Scientific content carries a dated evidence scope, without a change log.
    output=ROOT/'fwm_quotient_structure_v4.html'
    output.write_text(page,encoding='utf-8',newline='\n')
    print(f'Built {output.name}: {len(GROUPS)} groups, {len(nodes)} nodes, {len(EDGES)} edges, {len(used)} source keys.')

if __name__=='__main__':
    build()
