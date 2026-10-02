"""Scientific definitions, conditional submodels, and evidence for the FWM atlas.

This module is authored independently of earlier atlases. The builder embeds all
records in its standalone HTML; no Python or network is needed by the reader.
"""

SOURCE_CATALOG = {
    'logic:order': ('정보 준순서와 동치류의 부분순서 · 본문 정의', '#q-information'),
    'logic:quotient': ('정확한 동치·몫·합성에 필요한 조건 · 본문 명제', '#q-quotient'),
    'logic:assumptions': ('가정의 함의 순서와 조건부 추론 · 본문 정의', '#q-assumptions'),
    'logic:noise': ('고정 관측에서의 Loewner 단조성·최소 잡음 반례', '#q-noise'),
    'logic:counterexamples': ('스칼라 점수·rate 근사·N형의 적용 경계', '#q-counterexamples'),
    'gc:foundations': ('명시적 reservoir와 GKSL 검사', '../grand_challenge/conventions.md'),
    'gc:diffusion': ('Atomic Einstein diffusion와 독립 QRT 유도', '../grand_challenge/derivation.md'),
    'gc:normalization': ('Pump/weak dipole와 manifold 평균', '../grand_challenge/normalization_derivation.md'),
    'gc:periodic': ('Finite-seed periodic drift·noise·두-band 광장', '../grand_challenge/periodic_field_derivation.md'),
    'gc:spatial': ('정지 원자 공간 projection과 nonclosed geometry', '../grand_challenge/spatial_field_derivation.md'),
    'gc:uncertainty': ('공동 입력 불확도·stationary spatial channel', '../grand_challenge/spatial_uncertainty_derivation.md'),
    'gc:transport': ('Characteristic transport·boundary와 cross-position noise', '../grand_challenge/transport_derivation.md'),
    'gc:rb-transport': ('Moving reduced Rb 원자의 finite segment 모델', '../grand_challenge/rb_transport_derivation.md'),
    'gc:inflow': ('Maxwell boundary flux·marked Poisson 합산', '../grand_challenge/inflow_derivation.md'),
    'gc:adjoint': ('Backward observable의 독립 source·response 적분', '../grand_challenge/adjoint_transport_derivation.md'),
    'gc:cf4': ('정확 spectral source 적분과 CF4 경로 평가', '../grand_challenge/exponential_transport_derivation.md'),
    'gc:thermal': ('실제 thermal 경로·경계 조건·재사용 계약', '../grand_challenge/rb_thermal_ensemble.md'),
    'gc:ensemble': ('세 seed의 p2 복소 원자 spectrum·수렴 판정', '../grand_challenge/thermal_campaign_v2/ensemble-p2-three-seeds.json'),
    'gc:checkpoint': ('독립 수학 probe·경계 이력·적용 범위', '../grand_challenge/checkpoint_2026_09_18.md'),
    'gc:experiment': ('독립 입력·out-of-sample·full-physics 완료 조건', '../checklist.json'),
}

INTRO = '''<section class="sheet logic-primer" id="logic-primer" aria-labelledby="logic-title">
<p class="eyebrow">RELATIONS HAVE TYPES</p><h2 id="logic-title">한 줄의 정확도 사다리로 정렬하지 않는다</h2>
<p>이 지도는 <strong>계산 의존성</strong>, <strong>정확히 보존되는 정보</strong>, <strong>가정에 따른 축약</strong>, <strong>수치 해상도</strong>, <strong>실험 근거</strong>를 분리한다. 같은 출력값 하나, 작은 방정식 잔차, 양의 공분산만으로 두 물리 모델을 동치로 놓지 않는다.</p>
<div class="relation-grid">
<a href="#q-information"><b>정보 순서 A ≼ B</b><span>B의 출력으로 A의 출력을 정확하게 복원할 수 있음. 같은 입력 도메인과 관측 계약 필요.</span></a>
<a href="#q-equivalence"><b>표현의 동치 A ≃ B</b><span>허용된 프레임·기저 변환 아래 전체 관련 상태·관측·noise·검출이 보존됨.</span></a>
<a href="#q-numerics"><b>해상도 n ⊑ n′</b><span>같은 물리 모델의 더 세밀한 계산. 관측 오차의 단조 감소나 실험 정확도를 뜻하지 않음.</span></a>
<a href="#q-evidence"><b>검증 E ⊢ C | Γ</b><span>명시한 가정 Γ와 영역에서 증거 E가 주장 C를 지지함. 모든 주장 사이의 전역 순서 아님.</span></a>
</div><p class="condition"><b>현재 성립하는 범위</b>명시된 reduced GKSL 원자와 조건부 Gaussian 광장에서 microscopic diffusion·QRT·검출 계산이 연결되어 있다. 이동하는 열원자 전체의 quadrature는 미수렴이며, nonlocal Maxwell·full atom·독립 실측 입력을 합친 절대 squeezing 예측은 성립하지 않았다.</p>
<p><a href="#group-structure">정의와 증명 읽기 →</a> · <a href="#group-microscopic">조건부 quantum 이론 →</a> · <a href="#group-transport">열원자 transport와 남은 연결 →</a></p>
</section>'''

EXTRA_SPINE = [
    dict(source='input', target='microscopic', label='명시한 원자·reservoir', kind='data', pending=False),
    dict(source='microscopic', target='micro', label='조건부 M·D', kind='limit', pending=False),
    dict(source='microscopic', target='transport', label='같은 H·jump', kind='data', pending=False),
    dict(source='transport', target='micro', label='nonlocal 광장 연결', kind='data', pending=True),
]


def update_context(context):
    context=context.replace('<tbody>', '''<tbody>
        <tr><td><a href="#r-poles">Spectral resolvent 항등식</a></td><td><b>1.0 × 10⁻¹²</b><p>상대 Frobenius 차이</p></td><td>16×16 원자 Liouvillian의 RF 0.1–4 MHz에서 spectral expansion과 direct solve 비교. 전체 이론의 보편적 오차가 아닙니다. <code>tab:analytic-verification</code></td></tr>
        <tr><td><a href="#a-gate">Floquet·SI 전파 구현 대조</a></td><td><b>≤ 1.64 × 10⁻¹⁴</b><p>N<sub>F</sub>=3 fixture 최대 차이</p></td><td>별도 finite-block·literal SI Maxwell/Q 경로와 생산 계산의 대조. 원자 Hamiltonian·산일자 가정은 공유하므로 독립 장치 검증이 아닙니다.</td></tr>''',1)
    return (context.replace('현재 축약 모델이 재현하지 못합니다.',
        'no-fit 축약 모델이 일반적으로 재현하지 못합니다. 단일 운전점 gain에 맞춘 C_mix=0.5594938027은 독립 검증이 아닙니다.')
        .replace('미시 M<sub>q</sub>(Ω), D(Ω), 입력 상태 및 검출 응답의 완결된 연결이 필요합니다.',
        '조건부 원자·광장에서는 미시 M<sub>q</sub>(Ω), D(Ω)와 검출이 연결되어 있습니다. 실제 열원자 장치의 수렴된 nonlocal 광장·입력·검출 연결은 남아 있습니다.')
        .replace('수학적 정식화를 더한다면 객체·사상·합성과 오차의 합성 규칙을 별도로 정의해야 합니다.',
        '이 문서의 정보 준순서, 가정 함의, 정확한 동치류는 아래의 정의된 영역에서만 사용합니다. 수치 오차의 근접성은 동치관계가 아니며 자동으로 합성되지 않습니다.')
        .replace('이론적으로 평균장을 예측하면 D → F가 필요합니다.',
        'D를 평균장, X를 gain 진단, E를 공분산, F를 검출이라 부르자. 평균장으로부터 검출을 예측하면 D → F가 필요합니다.'))


def extend(groups, edges, kinds, edge_kinds, numbers, code):
    """Append complete records and refine stale scope statements in place."""
    if any(g['id'] == 'structure' for g in groups):
        return
    kinds.update(conditional='조건부 연구 구현', definition='정의·명제')
    edge_kinds.update(logic='논리 전제', equivalence='정확한 표현 동치',
                      refinement='동일 모델의 수치 세분', evidence='범위가 있는 증거')
    numbers.update(structure='L', microscopic='Q', transport='T')
    code.update(structure=('docs/FWM physics and analytic reconstruction/fwm_quotient_structure_v4_science.py',
                          '명시적 준순서·동치류·조건부 추론; 물리 solver 아님'),
                microscopic=('gabes/quantum/diffusion.py', 'explicit jump-products, atomic drift and independent QRT'),
                transport=('gabes/fwm_quantum/transport_ensemble.py', 'characteristic packets, current path evidence and conditional ensemble gate'))

    lookup = {n['id']: n for g in groups for n in g['nodes']}
    def amend(id, **updates):
        lookup[id].update(updates)

    amend('m-drift', body='같은 구동 원자 상태의 주파수 응답을 정준 광학 mode로 사상한다. Single-velocity·two-band 및 stationary spatial 조건부 경로는 구현되어 있다. 생산 finite-seed 비율과 일반적으로 같은 drift는 아니다.',
          condition='움직이는 thermal 원자의 nonlocal response를 실제 광장에 연결하는 일반 이론은 별도 미완료다. 조건부 국소 M_q의 존재가 이 연결을 대신하지 않는다.')
    amend('m-reservoir', body='교환관계와 symmetrized 상관은 다른 정보다. 명시적 jump operators에서 source별 atomic ordered diffusion을 유도하고 별도 QRT와 대조한다. Optical elimination과 collective normalization을 거쳐야 field D가 된다.',
          condition='K만으로 N_res를 추론하지 않는다. 원자 Markov diffusion과 원자 제거 후 주파수 의존 field diffusion을 구별한다.')
    amend('o-claim', condition='실제 열원자 장치의 절대 gain·S₋는 아직 미검증이다. 조건부 microscopic readout은 존재하지만 thermal quadrature, nonlocal Maxwell, full atom과 실험 입력의 검증을 대체하지 않는다.')
    amend('d-chi', condition='이 식은 생산 평균장 convention이다. 연구의 uniform-Zeeman-RMS 모델은 pump/weak field에 같은 dipole을 사용하며 manifold별 평균을 다시 정한다. 1/12와 RMS convention을 섞지 않는다. 추가 penalty 기본값 1과 경험적 0.74는 독립 측정값이 아니다.')
    amend('i-drive', body='생산 경로는 온도로 pure-⁸⁵Rb 밀도·Doppler·일부 완화율을 정한다. 조건부 microscopic/thermal 경로는 density와 temperature를 별도 입력으로 선언한다. 파워와 1/e² 반경은 Rabi 구동을 정하며 beam area·transit law는 별도 물리 가정이다.')
    lookup['r-geometry']['body'] += ' 이 fixture의 해석 rms 1.380173564 MHz와 40×40·5σ_v 구적 1.380163304 MHz의 차이는 −7.43×10⁻⁴%다.'
    lookup['r-poles']['refs'].append('tab:analytic-verification')
    for g in groups:
        # Section names are stable scientific references; part numbers are not.
        if ' · ' in g['part']:
            prefix, rest = g['part'].split(' · ', 1)
            if prefix and all(c in 'IVX' for c in prefix):
                g['part'] = rest
        if g['id'] == 'micro':
            g['brief'] = '선언된 국소 모델에서 구현된 미시 경로와 실제 움직이는 열원자 장치에 필요한 연결을 구분한다. 동일 pumped generator가 response와 source별 noise를 함께 정한다.'
            g['evidence'] = '조건부 single-velocity·finite-seed two-band·stationary spatial M/D는 독립 검사되었다. 일반 thermal nonlocal field와 실험 절대 squeezing은 미완료다.'
        if g['id'] == 'output':
            g['evidence'] = '조건부 microscopic S₋는 계산 가능하다. 생산 앱의 physical squeezing, full thermal 예측과 실험 bandwidth는 아직 미인증이다.'
        if g['id'] == 'mean':
            g['evidence'] = '복소 전파의 구현 일관성과 실험 예측은 별도다. 단일 gain에 fitted C_mix를 적용한 결과는 no-fit 타당성의 근거가 아니다.'

    def group(id, title, subtitle, kind, brief, evidence):
        g = dict(id=id, title=title, subtitle=subtitle, kind=kind,
                 part={'structure': '수학적 관계와 조건부 추론',
                       'microscopic': '원자 reservoir에서 조건부 광장까지',
                       'transport': 'Moving atoms, open boundary and field closure'}[id],
                 brief=brief, evidence=evidence, nodes=[])
        groups.append(g)
        return g

    def node(g, id, title, out, math, inputs, units, body, condition, refs,
             kind=None, evidence=None):
        i = len(g['nodes'])
        n = dict(id=id, group=g['id'], title=title, out=out,
                 xy=(1+i % 2, 1+i//2), math=math, compact=math.split('\n')[0],
                 inputs=inputs, units=units, body=body, condition=condition,
                 refs=refs.split(), kind=kind or g['kind'], evidence=evidence or g['evidence'])
        g['nodes'].append(n)
        lookup[id] = n

    def edge(s, t, label, kind='data', pending=False):
        edges.append(dict(source=s, target=t, label=label, kind=kind, pending=pending))

    g = group('structure', '순서·동치·조건부 추론', '비교할 객체와 관계를 먼저 정의', 'definition',
              '계산 도식의 화살표는 모두 하나의 순서관계가 아니다. 물리 specification, 표현, 수치 실행과 증거를 구분한 뒤 각각의 관계만 합성한다.',
              '아래는 정의와 명제다. 현재 코드의 모든 관계가 정리의 가정을 만족한다는 자동 인증은 아니다.')
    node(g, 'q-objects', '비교의 객체와 관측 계약', 'M · Γ · n · E',
         '비교 context 𝒦=(D, Γ, 관측 interfaces, 허용 Π) 고정\nM=(H, {L_r}, boundary, optical modes, P_M)\nP_M : D → O ; n=수치 설정 ; E=증거 묶음',
         '입력 도메인 D, 허용 가정 Γ, state/noise/field model, 관측공간 O',
         'P는 지정 단위·기저·검출 조건의 예측 map; n과 E는 M 자체가 아님',
         '예측 대상에 평균 파워만 둘지, ±Ω의 full covariance와 SQL까지 둘지 먼저 정한다. 상태가 달라도 한 점 gain이 우연히 같을 수 있다. 비교 가능한 관측과 입력 식별을 선언하기 전에는 우열이나 동치를 말하지 않는다.',
         '다른 isotope·입력 domain·detector convention은 대응 map 또는 공통 제한 영역 없이 비교하지 않는다. 물리 상수·area·noise를 단위 규약으로 숨기지 않는다.',
         'logic:order eq:frequency-definitions eq:flux-normalization')
    node(g, 'q-information', '정확한 정보 회수의 준순서', 'A ≼ B ⇔ P_A=π∘P_B',
         'A ≼ B ⇔ ∃ π : O_B→O_A, P_A=π∘P_B on D\nA ≼ A (π=id)\nA≼B, B≼C ⇒ A≼C (π_AC=π_AB∘π_BC)',
         '고정 context 𝒦의 공통 D·Γ·관측 interfaces, 항등·합성에 닫힌 허용 Π의 π',
         '준순서(preorder): 반사적·추이적; total order 아님',
         'B가 보존한 출력 정보로 A를 정확히 회수한다는 뜻이다. 허용 π는 입력에 따른 임의 보상이나 target fitting을 포함하지 않는다. 같은 상태의 full complex field로 power를 구할 수 있지만 power만으로 phase를 유일하게 복원하지 못한다.',
         '노드의 데이터 전달 화살표와 ≼를 같게 두지 않는다. 더 많은 물리 가정을 넣은 다른 예측 P_B가 기존 P_A를 정확히 포함한다고 자동 가정하지 않는다.',
         'logic:order')
    node(g, 'q-assumptions', '가정의 함의 순서', '강한 가정은 좁은 적용 영역',
         'Γ₁ ⪯_Γ Γ₂ ⇔ Mod(Γ₂) ⊆ Mod(Γ₁)\nΓ₂ ⊢ Γ₁ : Γ₂가 Γ₁보다 강한 제약\nΓ ⊢ C : 선언한 조건에서만 결론 C',
         '공통 물리 해석을 가진 가정들의 집합과 만족하는 모델 클래스',
         '가정식의 준순서; Mod가 같은 논리적 동치류에서는 부분순서',
         '상수 drift, 독립 Markov reservoirs, bright coherent seed, symmetric loss는 각각 적용 영역을 제한한다. 조건을 더하면 닫힌 식이 생길 수 있지만 실제 장치가 조건을 만족한다는 근거도 필요하다.',
         '공집합인 Mod(Γ)는 모순된 가정이며 유효한 물리 모델이 아니다. 임의의 두 허용 모델의 join/meet가 존재한다고 주장하지 않는다. union/intersection을 물리적 격자로 부르려면 폐쇄성과 일관성을 증명해야 한다.',
         'logic:assumptions')
    node(g, 'q-equivalence', '표현의 정확한 동치', '프레임·기저·관측을 함께 변환',
         'ρ′=UρU† ; O′=UOU† ⇒ Tr(O′ρ′)=Tr(Oρ)\nR′=SR ; V′=SVSᵀ ; m′=S⁻ᵀm ⇒ m′ᵀV′m′=mᵀVm\ncanonical S : SΩ_RSᵀ=Ω_R',
         'invertible 표현 변환, H·jump·state·noise·readout의 일관된 사상',
         '같은 physical observable을 보존하는 표현; 단위 변경도 명시',
         '변환과 역변환 아래 동역학과 허용 관측이 일치하면 같은 물리 문제의 표현이다. 시간 의존 U는 Hamiltonian의 gauge 미분항도 변환한다. Power-only agreement나 mean-only gauge parity로 전체 noise·detector 동치를 선언하지 않는다.',
         'pump-only minus의 정적/Floquet 비교는 대응되는 harmonic과 상태 범위에서만 유효하다. finite seed·plus branch·누락 optical ports까지 동치를 확장하지 않는다.',
         'logic:quotient eq:pump-frame-transform tab:pump-reference-diagnostics')
    node(g, 'q-quotient', '몫의 정의와 보존되는 합성', '[A] ≤ [B]를 정당화',
         'A ∼ B ⇔ A≼B and B≼A\n[A] ≤ [B] ⇔ A≼B\n반사·추이: ≼에서 상속; 양방향이면 [A]=[B] (반대칭)\nπ₂∘π₁ 허용 ⇒ 대표를 바꾸어도 관계가 동일',
         '고정 𝒦=(D,Γ,관측 interfaces,Π) 위의 준순서; Π는 항등·합성에 닫힘',
         '동치류의 부분순서; 모든 표현을 하나로 합치는 연산 아님',
         'A′∼A와 B∼B′이면 A′≼A≼B≼B′이므로 순서가 대표 선택에 무관하다. 이 동치는 고정한 관측 계약의 정보 회수에 관한 것이며, 계약에 없는 모든 물리 관측의 동치를 뜻하지 않는다. Gaussian channel을 표현별로 몫내어 합성하려면 입출력 interface가 맞고 해당 동치가 합성과 양립하는 congruence인지 따로 확인한다.',
         '규약을 옆 패널로 이동하는 것은 수학적 quotient가 아니다. 임의 ε-근접은 추이적이지 않으므로 같은 증명으로 몫을 만들 수 없다. 계산 DAG의 그룹 접기도 별개의 시각화다.',
         'logic:order logic:quotient')
    node(g, 'q-numerics', '수치 세분과 오차', 'M 고정, n만 변경',
         'n=(N_F, velocity rule/cutoff, Δt, RF grid, z step, tolerances)\nn ⊑ n′ : 선언한 각 해상도가 n′에서 더 세밀함\n||x_n−x_n′||≤ε ≠ ||x_n−x_exact||≤ε',
         '동일 H·reservoir·boundary·입력·observable, 해상도별 원시 결과',
         '벡터 해상도의 부분순서; 여러 축이 반대로 변하면 비교 불가',
         '해상도 순서는 계산 설정을 정렬한다. 관측량의 실제 오차가 반드시 단조 감소하는 것은 아니다. 인접 결과 일치는 선언한 검사의 통과이지 참 해에 대한 오차 상계가 아니다. Richardson 평가에도 수렴 차수·점근 영역 등의 조건이 필요하다.',
         '4→24준위, 1D→일반 비공선 transport, pump-only→finite seed, density law 변경은 물리 축이다. 같은 모델의 더 작은 time step과 섞지 않는다. 서로 다른 Sobol seed는 nested refinement가 아니다.',
         'logic:order sec:floquet-convergence gc:ensemble')
    node(g, 'q-evidence', '주장마다 다른 증거', 'E ⊢ C | Γ,D,ε',
         'C_num: 같은 방정식을 주어진 오차 기준으로 계산\nC_CP: 선언한 quantum channel이 허용됨\nC_QRT: 같은 GKSL model의 독립 대수 경로 일치\nC_exp: untouched 조건의 gain+S₋가 실험 불확도 안에서 예측됨',
         '검증 대상, domain, model/input/source identity, tolerance, 실제 결과',
         '증거의 종류와 범위; 물리 정확도 %의 전역 척도 아님',
         '잔차, 구현 비교, 물리 불변량, 수치 수렴, 독립 실험은 서로 다른 주장을 검사한다. 같은 generator를 쓰는 QRT와 diffusion 구현의 일치는 공유한 실제 장치 가정의 오류를 배제하지 않는다. provenance hash는 출처를 확인하며 물리 타당성을 증명하지 않는다.',
         'CP/QRT 통과 ⇒ 실험 정확도라는 함의는 없다. 데이터 공개·설정 동결·불확도·no-refit 조건을 갖춘 실험 gate를 따로 요구한다.',
         'gc:checkpoint gc:experiment tab:validation-hierarchy')
    node(g, 'q-noise', '잡음의 부분순서와 한계', '고정 관측에서만 단조',
         'W₁ ⪯ W₂ ⇔ W₂−W₁ ⪰ 0\n동일 m, SQL>0: S₂−S₁=m†(W₂−W₁)m/SQL ≥ 0\nX=0, Y_s=diag(s,s⁻¹)/2, s>0: CP 허용\ns≠t이면 Y_s−Y_t의 두 대각 부호가 반대',
         '같은 transfer·input·mode·carrier·detector 조건의 두 covariance',
         'Loewner 부분순서; 작은 trace나 한 측정값의 순서와 다름',
         '두 허용 covariance가 일반적으로 비교 가능하지 않다. X=0인 prepare channel은 모든 순수 squeezed covariance Y_s를 허용한다. 서로 다른 s의 잡음은 비교 불가능하고, 모든 s보다 작은 PSD Y는 0이어야 하나 0은 uncertainty를 위반한다. 따라서 이 허용 집합에는 보편적인 least noise가 없다.',
         '최소 원소(least)와 더 줄일 수 없는 원소(minimal)를 구별한다. 특정 scalar 관측의 최적화는 별도 문제다. K의 defect를 채운 독립 vacuum D를 실제 원자 noise의 보편적 하한이라 하지 않는다.',
         'logic:noise eq:physical-realizability eq:minimal-dilation')
    node(g, 'q-and', '결론은 전제들의 합류', 'AND와 대안 OR를 구별',
         '{same mean, full covariance, η,H,w,SQL} ⇒ S₋(Ω)\n{M,D constant} ⇒ W=integral=φ₁=Van Loan\n{mean gain only} ⇏ S₋ ; {CP} ⇏ experimental validation',
         '한 결론이 실제로 소비하는 모든 전제와 조건',
         'hyperedge/전제 집합; 일반 선 하나는 전체 충분조건 아님',
         '여러 입력 화살표는 서로 바꿔 쓸 수 있는 증거 목록이 아니라 동시 필요 입력이다. φ₁·고유기저·Van Loan은 조건이 맞을 때 같은 문제를 계산하는 대안이며 세 noise를 더하는 경로가 아니다. 원자와 검출의 필요 입력 중 하나라도 없으면 결론의 범위를 제한한다.',
         '추이적 축약 그림에서 지운 선도 원시 값 전달에 필요할 수 있다. reachability만 보존한 Hasse 그림이 계산의 모든 인자·공유 입력을 보존한다는 뜻은 아니다.',
         'logic:assumptions eq:covariance eq:spectrum')
    node(g, 'q-counterexamples', '축약과 단일 점수의 반례', 'N형 · rate 근사 · gain 점수',
         'D<X, D<F, E<F만인 4점 순서: induced N형\nscore=clip_[0,1][(G_p−g_c)²/(G_p+g_c)]\n(G_p,g_c)=(0.2,0.1) ⇒ score=0.0333\npassive coherent attenuation + vacuum ⇒ matched SQL ratio=1',
         '의존관계·회수 가능량·축약 scale hierarchy의 명시',
         '반례; 세 공식을 하나의 물리 모델로 동일시하지 않음',
         '평균장 D와 covariance E가 검출 F에서 합류하고 gain 진단 X가 별도로 있을 때, 네 점의 해당 제한 순서에는 N형이 남을 수 있다. 전체 구조의 SP/N-free를 가정하지 않는다. Gain 점수의 수동 영역 sub-SQL 값은 microscopic covariance가 아니다. Rate+Sylvester 제거는 eliminated coherence가 유지 변수보다 빠르다는 scale hierarchy를 요구한다. Ω_pump/(2π)=0.708 GHz와 excited splitting 361.58 MHz인 알려진 조건은 약한 혼합 제거를 정당화하지 못한다.',
         'N형 명제는 선언한 네 점 관계에 한정한다. 특정 식의 실패를 모든 rate model의 실패로 확대하지 않는다. 이 반례들은 실험의 squeezing 크기를 정하지 않는다.',
         'logic:counterexamples eq:gain-gap-observable sec:physicality-blind')

    g = group('microscopic', '조건부 microscopic 이론', '같은 generator의 response와 noise', 'conditional',
              'Reduced GKSL atom의 미시 noise와 field readout을 수학적으로 연결한다. 각 확장은 어떤 원자·mode·경계·선형화 조건을 사용하는지 남긴다.',
              '명시된 fixture의 독립 QRT·null/passive·CP·수렴 검사. 실제 hot-vapor apparatus의 완전성 또는 실험 정확도 인증이 아니다.')
    node(g, 'q-gksl', '명시적 원자 generator', 'H와 source별 L_r',
         'ℒρ=−i[H/ℏ,ρ]+Σ_r(L_rρL_r†−{L_r†L_r,ρ}/2)\nL_r: √(s⁻¹), H/ℏ: s⁻¹ ; Trρ=1, ρ⪰0',
         '원자 basis, 출처와 적용 영역이 있는 radiative/reset/dephasing reservoirs',
         'finite-dimensional Markov generator; state와 bath의 provenance',
         '직접 coherence decay만 더하는 수정은 총 CP 한 점의 통과만으로 독립 reservoir가 되지 않는다. 각 noise source는 명시한 collapse 또는 PSD Kossakowski matrix로 추적한다. H와 dissipator, state를 같은 모델에서 조립한다.',
         'GKSL 허용성은 선언한 Markov 모델의 일관성이다. 실제 collision memory·velocity correlation·Zeeman structure가 무시 가능하다는 실험 증명은 아니다.',
         'gc:foundations gc:diffusion')
    node(g, 'q-diffusion', 'Jump products의 ordered diffusion', 'A · Dᵍ · Dˡ · reservoir 합',
         'A_ij=Tr[F_i ℒ(F_j)], Tr(F_iF_j)=δ_ij\nD⁾_ij=Σ_r Tr[ρ [L_r†,F_i][F_j,L_r]]\nĊ=AC+CAᵀ+D⁾ ; K=C−Cᵀ\nD⁽=(D⁾)ᵀ (같은 equal-time Hermitian basis)',
         'complete traceless Hermitian atomic basis, 같은 ρ와 L_r',
         'A: s⁻¹; D atomic covariance/s; K atomic commutator',
         'Drift는 generator에서, diffusion은 reservoir jump products에서 각각 계산한다. Complete atomic basis의 adjoint evolution은 affine-linear이므로 exact second moment와 QRT를 닫을 수 있다. Atomic Gaussian fourth-cumulant 가정을 여기서 요구하지 않는다.',
         'Atomic K는 traveling optical metric J가 아니다. 원자 Markov D와 elimination 뒤의 field D(Ω), ordered와 symmetrized covariance를 혼용하지 않는다. D를 A의 commutator defect만으로 맞추지 않는다.',
         'gc:diffusion')
    node(g, 'q-qrt', '독립 QRT와 극한 검사', '동일 모델의 다른 대수 경로',
         't≥0: ⟨O_i(t)O_j(0)⟩=Tr[O_i e^(ℒt)(O_jρ_ss)]\nconnected correlation: 평균 곱 차감\nS(ω)=∫dt e^(iωt) C(t)',
         '정상/periodic 상태, 같은 generator, 적절한 stationary projection',
         '원자 spectral correlation; 같은 ordering·Fourier 규약',
         'Jump-product/resolvent 계산을 full-Liouville two-time correlation과 비교한다. 2·3·4준위 random GKSL 30모델에서 모델당 generator 주파수 5점(−13,−3,0,2,11 rad/s), 총150회 비교 최대 상대차 9.985×10⁻¹⁶. No-atom·zero-length·passive와 channel 조건은 별도 검사다.',
         '이 수치는 synthetic generator의 대수 검사다. 150개의 실험 RF bin 검사나 실제 원자의 누락 physics 검증이 아니다. 같은 H·state의 오류는 두 방식이 공유할 수 있다.',
         'gc:checkpoint gc:diffusion', kind='test')
    node(g, 'q-rms', 'Pump·weak dipole의 일관성', 'reciprocity와 manifold mean',
         'I₀=2P/(πw²), Ω=dE₀/ℏ\nW_FF′=p_F d_J² S_FF′/3\nd_RMS,FF′=d_J√(S_FF′/3)\n3C_F²=(2F+1)S_FF′/3',
         '같은 dipole table, spontaneous branching, state의 manifold population',
         'RMS amplitude는 polarization-resolved signed CG amplitude가 아님',
         'Pump/weak implicit dipole 비 1.9993734141인 불일치는 내부 commutator가 맞아도 가능하다. Sublevel sum을 manifold mean으로 바꾸는 분모는 F=2에서5, F=3에서7이며 공통1/12로 두 값을 동시에 맞출 수 없다. Independent uncoupled I⊗J reference의 strength 대조 잔차는1.12×10⁻¹⁶ 이하.',
         '균일한 비편극 weak absorption의 평균을 맞춘 양의 RMS는 강한 pump의 Zeeman population·CG 부호·loop interference를 복원하지 않는다. 24→4준위 full strong-pump 정리는 증명되지 않았다.',
         'gc:normalization')
    node(g, 'q-periodic', 'Finite-seed periodic noise', 'harmonic 사이 상관 유지',
         'H(t)=H₀+Ve^(−iνt)+V†e^(iνt)\nD_q⁾, A_q, C_q ; Ċ=AC+CAᵀ+D⁾\n⟨δF_h(ω)δF_k†(ω′)⟩=2πδ(ω−ω′)S_hk⁾(ω)',
         '같은 finite-seed periodic state, explicit baths, harmonic cutoffs',
         'ν=−ω_HF+δ; independent RF frequency retained',
         'Periodic atomic drift와 diffusion을 독립적으로 만들고 harmonic cross-correlation을 유지한다. 두 optical bands에 대한 full Nambu elimination, nonlinear probe/conjugate mean, adaptive covariance propagation과 bright readout이 조건부로 연결되어 있다.',
         'Pump는 classical·fixed다. Additional physical optical ports와 full depleted thermal cell은 별도다. 작은 mean 변화가 각 noise 성분의 작은 변화를 보장하지 않는다. q≠0 Fourier block 하나는 density matrix나 PSD matrix일 필요가 없다.',
         'gc:periodic')
    node(g, 'q-spatial', 'Stationary spatial mode projection', '선언한 geometry의 광장',
         'mode-projected response + source covariance\nlocal M_q(z,Ω), D(z,Ω) → T_q,W → detected S₋',
         '정지 atomic center, reduced atom, 실제 signed wavevectors, 선언한 mode profiles',
         '광학 mode와 auxiliary phase coordinate를 구별',
         'Stationary centers에서 nonclosed wavevector geometry의 물리 optical projection과 두-band mean/noise/readout을 연결한다. 이는 phase coordinate 자체를 추가 광학 mode로 부르는 것을 피한다. Synthetic top-hat profiles와 fixed pump 조건에서 검산되었다.',
         'Moving thermal transport와 diffraction·walkoff·mode completeness는 포함되지 않는다. 세 RF sample의 sub-SQL 값으로 intrinsic bandwidth를 정할 수 없다.',
         'gc:spatial')
    node(g, 'q-gaussian', 'Gaussian field와 photocurrent', '2차 noise와 4차 closure 구별',
         'a=β+d ; I=|β|²+β*d+βd†+d†d\nGaussian Wick: ⟨dddd⟩=pair contractions\nS_quad와 spontaneous mean의 SQL를 함께 보정',
         'field covariance, declared optical collection band, bright carrier or full quadratic model',
         'Gaussian field는 atomic GKSL exact second moments보다 추가 가정',
         'Bright photocurrent는 carrier에 선형인 항을 사용한다. Finite-band Gaussian quadratic readout은 spontaneous photon mean·shot noise·convolution을 추가하고 같은 SQL도 바꾼다. Atomic second-moment 정확성은 connected field fourth cumulant가0임을 뜻하지 않는다.',
         'Unfiltered fluorescence tail와 비Gaussian cumulant의 full bound는 열려 있다. Vacuum-input TMSV exact count 식과 bright-seed 선형화의 적용 영역을 구별한다.',
         'gc:periodic gc:checkpoint eq:photocurrent')
    node(g, 'q-uncertainty', '입력·오차·검출의 공동 조건', '조건부 sensitivity와 실험 예측',
         'Σ_y≈J_θ Σ_θ J_θᵀ (local linear uncertainty)\nS₋=PSD(i_p−wi_c)/SQL(same means,H,w)',
         '독립 input ledger, correlations, model validity, numerical error, detector calibration',
         '입력 불확도·수치 오차·model discrepancy는 별도',
         'Stationary conditional channel의 공동 입력 sensitivity와 uncertainty propagation이 구현되어 있다. 독립 beam-profile fit은 허용되지만 target gain/squeezing을 맞춘 coupling/noise/loss는 no-fit 입력이 아니다. 같은 평균과 효율의 검출 SQL를 사용한다.',
         '가정한 uncertainty는 실측 uncertainty가 아니다. 국소 선형 전파는 nonlinear distribution과 model discrepancy를 검증하지 않는다. 실제 장치 절대 예측에는 provenance·measurement applicability·heldout가 필요하다.',
         'gc:uncertainty gc:experiment')

    g = group('transport', '열원자 transport와 광장 연결', '경로 해 · ensemble · nonlocal 관측', 'conditional',
              '움직이는 원자의 입사 상태, 경로 전체의 구동과 다른 위치 사이 상관을 유지한다. 정밀한 개별 경로와 수렴된 열적 ensemble은 별도 관문이다.',
              '조건부 open-column reduced Rb stream. Thermal quadrature 미수렴; moving-atom Maxwell·실험 squeezing 인증 없음.')
    node(g, 't-characteristic', '유한 atomic characteristic', '입구부터 출구까지의 H(t)',
         'r(t)=r_in+vt, 0≤t≤τ\nρ̇(t)=ℒ[r(t),v]ρ(t), ρ(0)=ρ_in\npump envelope f(r)=exp(−r_⊥²/w²)',
         'signed k, v, entry phase/state, full residence τ, beam profile',
         '단일 원자 경로, 시간 s·위치 m; uniform transit reset과 다른 모델',
         'Gaussian pump를 경로 전체에서 연속 평가한다. 움직이는 reduced Rb atom의 원자 response와 source별 noise를 계산한다. 체류시간 자르기·대표속도·평균 pump로 치환하지 않는다.',
         '현재 thermal path는 pump-only reduced RMS atom이다. 실제 cell wall collision, re-entry history, finite seed saturation은 별도 물리 입력/모델이다.',
         'gc:transport gc:thermal')
    node(g, 't-boundary', '입사 상태와 경계 noise', 'boundary + distributed + cross terms',
         'C(t,t′)=propagated C_in + reservoir-generated correlations\n총 wavepacket covariance에 cross-position terms 유지',
         'ρ_in, entry flux statistics, finite-time propagators, source correlations',
         '원자 wavepacket coordinates는 canonical optical modes가 아님',
         '처음 들어오는 상태의 불확실성과 서로 다른 구간의 원자 상관이 출력에 영향을 준다. 특정 Rb fixture에서 boundary contribution을 지우면 greater norm18.06%, cross-segment covariance를 지우면333.46%가 바뀐다. 이것은 해당 모델에서 생략 항의 크기를 드러내는 대조다.',
         '위 수치는 실험 squeezing 오차가 아니다. 각 위치에서 새 독립 원자를 시작하는 local white-noise 대체가 타당한지는 별도 근사 검증이 필요하다.',
         'gc:rb-transport gc:transport')
    node(g, 't-poisson', 'Maxwell 유입과 Poisson 합', '독립 path covariance의 가중 합',
         'dλ=n f_M(v)(−v·n_face)₊ dA d³v\nstream spectrum = Σ_paths λ_path × unweighted packet\nnumber fluctuations include mean-pulse outer terms',
         'open boundary measure, phase law, independent marks, density once',
         'λ: s⁻¹; microscopic packet과 arrival rate를 구별',
         '여섯 면의 Maxwell incoming flux로 경로를 고른다. 독립 원자의 noise는 covariance 기여로 더하며 noise amplitude의 평균을 제곱하지 않는다. Marked Poisson number fluctuations는 connected atomic noise와 별도 mean-pulse term을 요구한다.',
         'Independent classes/marks 가정은 velocity-changing collision과 re-entry correlations가 있을 때 재검토한다. Density를 boundary rate와 optical normalization에 중복 곱하지 않는다. Occupancy를 nV에 맞추는 재정규화는 사용하지 않는다.',
         'gc:inflow gc:thermal')
    node(g, 't-solvers', '경로 해와 독립 참조', 'CF4 ↔ adjoint',
         'primary: 세 Δt ; independent: 두 adjoint tolerances\nmax ||CF4−adjoint||_rel = 6.080×10⁻⁸\n72/72 paths meet source/RF/metric gates',
         '같은 physical path·full τ·H·jumps·state·readout·source definitions',
         '8 metrics, RF0.1/1/4MHz; raw matrices, source-resolved errors',
         'Exact spectral source integrals와 CF4 path propagation을 backward-observable adjoint reference와 대조한다. Primary refinement budgets10⁻³, independent comparison5×10⁻⁶, independent refinement2×10⁻⁶를 적용한다. 독립 참조는 같은 numerical method의 copy가 아니다.',
         '이러한 작은 경로 잔차는 thermal quadrature 또는 boundary assumptions의 오류를 제한하지 않는다. Shared physical generator에 대한 독립 수치 비교다.',
         'gc:adjoint gc:cf4 gc:ensemble', kind='test')
    node(g, 't-thermal', '조건부 열린 Rb 기둥', 'T와 n을 별도 지정',
         'T=373 K, n=10¹⁸ m⁻³, L=12.5mm\nA=1.2×10⁻⁷m², side=346.41μm, pump w=530μm\nρ_in=diag(5/12,7/12,0,0)\nside-boundary I/I₀=0.6523…0.8077',
         '가상 square column, Gaussian undepleted pump, unpolarized inflow',
         'temperature/density/boundary는 measured apparatus defaults가 아님',
         '물질 유입 경계를 optical normalization area와 같게 둔 선언된 모델이다. 옆 경계에도 상당한 pump가 남아 있으므로 실제 cell에서 들어오기 전 pumping 이력을 무시할 수 있는지 확인해야 한다. Side midpoint의 정지 원자 probe에서는1μs 노출로 F=3 population0.58333→0.86975가 된다.',
         '정지 원자의 이력 probe는 thermal ensemble 오차 상계가 아니다. Sim121°C 장치의 wall·aperture·incoming state를 재현한 것이 아니다. Collection mode를 고정하고 물리 경계/입사 전 경로를 바꾸는 sensitivity가 필요하다.',
         'gc:thermal gc:checkpoint')
    node(g, 't-convergence', 'Ensemble 수렴의 별도 관문', '경로 PASS, ensemble 미수렴',
         'p2: 3scrambles ×24paths ×5evaluations =360evaluations\n경로별 primary3 + adjoint2 ; 총72고유 경로\n6 directed seed pairs: 0/6 meet 5%\n36 metric maxima=5.70…43.61%',
         'source/RF별 complex spectra, consecutive refinement, independent scrambles',
         '72고유 경로의 p2 quadrature; 서로 다른 Sobol seed는 동일 seed의 해상도 수렴을 입증하지 않음',
         '각 방향의 최악 metric은26.80–43.61%다. 5.70–43.61%는36개 metric 최대값 전체의 범위다. Actual cache·source·path gate·직접 합산은 검증되었지만 아직 thermal integral이5% 기준에 수렴하지 않았다.',
         '더 큰 p3/p4 quadrature를 평가해도 통과를 보장하지 않는다. 인접/seed 차이는 참 적분 오차의 엄밀한 상계가 아니다. 독립 Sobol seed는 실험 Δ·δ·T·pump 조건의 holdout이 아니다.',
         'gc:ensemble', kind='test')
    node(g, 't-nonlocal', 'Moving-atom Maxwell의 조건', '위치쌍 response/noise 필요',
         'δP(r,Ω)=∫dr′ χ(r,r′;Ω)δE(r′,Ω)+F_P(r,Ω)\n⟨F_P(r,Ω)F_P†(r′,Ω′)⟩ → spatial covariance kernel',
         '서로 다른 source/readout 위치, mode projection, canonical normalization, mean fields',
         'χ(r,r′) kernel과 이미 적분한 atomic stream matrix는 서로 다른 객체',
         '움직이는 원자는 한 위치의 광장에 반응하고 다른 위치에서 읽힐 수 있다. Path-integrated retarded response만으로 일반적인 local M(Ω)를 유일하게 정할 수 없다. 필요한 source/readout 위치 정보를 보존한 field elimination과 boundary-value propagation을 닫아야 한다.',
         'Local approximation은 correlation length와 field variation scale 등의 조건을 검증한 뒤 사용한다. 현재 thermal stream의 검사 통과를 optical commutator, SQL, S₋ 인증으로 넘기지 않는다.',
         'gc:transport gc:thermal', kind='need')
    node(g, 't-prediction', '절대 예측과 실험 판정', '조건을 동시에 충족',
         '{converged transport, same mean/M/D, mode collection,\n independent inputs ±uncertainty, detector/SQL} ⇒ conditional prediction\nuntouched experiments + no-refit metrics ⇒ empirical validation',
         '수렴된 이론 bundle, 단일 provenance ledger, gain과 전체 RF spectrum',
         'finite seed·self-consistent three-field depletion·full atom는 각각 범위 명시',
         '실험 목표는 독립 입력으로 두 gain과 S₋(Ω)를 함께 예측하는 것이다. Target-fitted0.74, C_mix, noise/loss 계수는 no-fit 입력으로 승격하지 않는다. 알려진−7.8dB 한 점은 blind holdout이 아니다. 중요한 physics의 원인은 no-refit add-one/leave-one-out 및 새 holdout로 확인한다.',
         'Thermal quadrature, nonlocal Maxwell, full-atom 적용성과 독립 실험 입력의 조건이 아직 함께 충족되지 않았다. 실패/반례는 적용 경계를 밝히지만 장치의 절대 예측을 입증하지 않는다. 데이터가 없으면 조건부 sensitivity까지만 주장한다.',
         'gc:experiment', kind='need')

    for s, t, label, kind in [
        ('q-objects','q-information','공통 domain·관측·허용 회수 map','logic'),
        ('q-objects','q-assumptions','물리 가정의 해석','logic'),
        ('q-information','q-quotient','양방향 정확 회수의 동치류','logic'),
        ('q-equivalence','q-quotient','표현 교체 시 관측과 합성 보존','logic'),
        ('q-assumptions','q-and','필요 조건의 동시 충족','logic'),
        ('q-numerics','q-evidence','고정 모델에서의 측정된 수치 오차','evidence'),
        ('q-noise','q-and','같은 m·SQL에서만 잡음 비교','logic'),
        ('q-counterexamples','q-information','gain 정보와 covariance 정보 구별','evidence'),
        ('q-counterexamples','q-and','모든 실제 입력선 보존','logic'),
        ('q-gksl','q-diffusion','같은 state와 jump operators','data'),
        ('q-gksl','q-qrt','같은 ℒ의 독립 Liouville 경로','data'),
        ('q-diffusion','q-qrt','ordered spectrum 비교 대상','validation'),
        ('q-rms','q-gksl','reciprocal Hamiltonian·readout','data'),
        ('q-gksl','q-periodic','periodic state·explicit bath','data'),
        ('q-diffusion','q-periodic','source별 harmonic noise','data'),
        ('q-periodic','q-gaussian','두-band field covariance와 밝은 carrier','limit'),
        ('q-spatial','q-uncertainty','조건부 spatial prediction Jacobian','data'),
        ('q-gaussian','q-uncertainty','같은 출력과 SQL의 sensitivity','data'),
        ('q-diffusion','m-diffusion','원자 제거·mode 사상 조건 필요','limit'),
        ('q-periodic','m-drift','two-band fixed-pump 조건부 drift','limit'),
        ('q-spatial','m-drift','stationary-center 조건부 drift','limit'),
        ('q-gaussian','f-current','bright approximation 또는 quadratic 보충','limit'),
        ('q-uncertainty','o-claim','독립 입력의 적용 범위','evidence'),
        ('q-qrt','q-evidence','동일 GKSL의 대수 비교','evidence'),
        ('r-frame','q-equivalence','minus pump-only gauge와 domain','evidence'),
        ('r-direct','r-poles','대각화 가능한 동일 resolvent의 표현','equivalence'),
        ('v-phi','v-vanloan','동일 상수 M,D,L에서 같은 W','equivalence'),
        ('m-vacuum','q-noise','허용 noise의 한 선택','evidence'),
        ('i-modes','q-objects','mode·normalization·observable contract','logic'),
        ('q-gksl','t-characteristic','위치·시간 의존 H와 같은 jumps','data'),
        ('t-characteristic','t-boundary','full finite-time evolution','data'),
        ('t-boundary','t-poisson','connected packet와 mean pulse','data'),
        ('t-characteristic','t-solvers','동일 경로의 독립 계산 대상','validation'),
        ('t-thermal','t-characteristic','boundary·state·beam·signed k','data'),
        ('t-thermal','t-poisson','Maxwell incoming measure','data'),
        ('t-poisson','t-convergence','grid별 source-resolved spectra','data'),
        ('t-solvers','t-convergence','경로 수치 gate는 필요조건','evidence'),
        ('q-numerics','t-convergence','consecutive grids·independent seeds','logic'),
        ('t-boundary','t-nonlocal','cross-position source/readout 정보','data'),
        ('t-convergence','t-prediction','수렴한 thermal 자료 필요','evidence'),
        ('t-nonlocal','t-prediction','field closure와 canonical covariance','logic'),
        ('q-uncertainty','t-prediction','independent ledger·detector uncertainty','logic'),
        ('q-and','t-prediction','모든 전제의 동시 충족','logic'),
        ('t-convergence','q-evidence','observed refinement/scramble 오류','evidence'),
    ]:
        edge(s, t, label, kind, pending=(s, t) in {
            ('t-boundary','t-nonlocal'), ('t-convergence','t-prediction'),
            ('t-nonlocal','t-prediction'), ('q-uncertainty','t-prediction')})

    theory_refs = {
        'q-gksl': 'eq:gc-gksl eq:gc-reset',
        'q-diffusion': 'eq:gc-atomic-drift eq:gc-einstein eq:gc-covariance-identity',
        'q-qrt': 'eq:gc-qrt eq:gc-ordered-spectra',
        'q-rms': 'eq:gc-si-drive eq:gc-manifold-strength',
        'q-periodic': 'eq:gc-lifted-generator eq:gc-period-qrt',
        'q-spatial': 'sec:gc-kinetic tab:gc-spatial',
        'q-gaussian': 'sec:gc-readout eq:gc-sql',
        'q-uncertainty': 'sec:gc-contract',
        't-characteristic': 'eq:gc-moving-phases sec:gc-transport',
        't-boundary': 'eq:gc-transport-covariance eq:gc-integrated-output',
        't-poisson': 'eq:gc-arrival-noise eq:gc-thermal-flux',
        't-solvers': 'eq:gc-adjoint-output eq:gc-cf4',
        't-thermal': 'sec:gc-ensemble',
        't-convergence': 'tab:gc-thermal-convergence',
        't-nonlocal': 'eq:gc-nonlocal-maxwell',
        't-prediction': 'sec:gc-contract',
    }
    for id, refs in theory_refs.items():
        lookup[id]['refs'].extend(refs.split())

    # Scientific distinctions are present in both graph data and static records.
    required = {'body','inputs','out','units','math','condition','evidence','refs'}
    assert all(required <= n.keys() and n['refs'] for g in groups for n in g['nodes'])
    assert len(lookup) == sum(len(g['nodes']) for g in groups)
