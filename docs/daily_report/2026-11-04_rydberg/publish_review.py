"""⟦AI⟧ Publish Korean Rydberg review with artifact checks and atomic state merge."""
from collections import Counter
from datetime import datetime, timezone, timedelta
import hashlib
import json
import os
from pathlib import Path
import re
import xml.etree.ElementTree as ET

ROOT=Path(__file__).resolve().parents[3]
OUT=Path(__file__).resolve().parent
REPORT='docs/daily_report/2026-11-04_scheme-3_rydberg-eit.md'
EVIDENCE=REPORT[:-3]+'_evidence.json'
PREVIOUS='docs/daily_report/2026-10-28_scheme-3_rydberg-eit.md'
STATE=ROOT/'docs/daily_report/review_state.json'
TOPIC='rydberg_eit'

def now(): return datetime.now(timezone(timedelta(hours=9))).isoformat()
def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p): return json.loads(Path(p).read_text('utf-8-sig'))
def write(p,d): Path(p).write_text(json.dumps(d,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8')
def link(label,path,line=None): return f'[{label}](<{(ROOT/path).as_posix()}{":"+str(line) if line else ""}>)'
def test_summary(path):
    root=ET.parse(path).getroot()
    cases=list(root.iter('testcase'))
    return dict(cases=len(cases),failures=sum(len(list(c.iter('failure'))) for c in cases),
        errors=sum(len(list(c.iter('error'))) for c in cases),
        skipped=sum(len(list(c.iter('skipped'))) for c in cases),
        modules=dict(Counter(c.get('classname') for c in cases)))

assert not (ROOT/REPORT).exists() and not (ROOT/EVIDENCE).exists(), 'Archive previous publication before a retry.'
snapshot=read(OUT/'snapshot.json'); probe=read(OUT/'probe.json')
checks=read(OUT/'checks.json'); integration=read(OUT/'integration_probe.json')
integration_checks=read(OUT/'integration_checks.json'); git=read(OUT/'git.json')
old=read(ROOT/(PREVIOUS[:-3]+'_evidence.json'))
initial_state_sha=sha(STATE); initial_state=read(STATE)
canonical={x['id']:x for x in read(ROOT/'docs/checklist.json')['items']}
hashes=snapshot['hashes']; closure={p:sha(ROOT/p) for p in hashes}
changed_during=[p for p,h in hashes.items() if closure[p]!=h]
assert not changed_during, changed_during
assert not snapshot['changes_since_initial']
assert checks['returncode']==integration_checks['returncode']==0
test_results=test_summary(OUT/'pytest.xml'); adapter_tests=test_summary(OUT/'integration_pytest.xml')
assert test_results['cases']==100 and adapter_tests['cases']==7
assert all(v==0 for s in (test_results,adapter_tests) for k,v in s.items() if k in ('failures','errors','skipped'))

assessment=[]
for case in integration['cases']:
    errors=[r['relative_error'] for r in case['comparisons']]
    coordinate_ok=abs(case['scheme_mode']['batch_optimum_probe_detuning_mhz']-
        case['scheme_mode']['actual_optimum_probe_detuning_mhz'])<1e-10
    resonance_ok=abs(case['resonance']['batch_transmission']-case['resonance']['native_transmission'])<1e-10
    assessment.append(dict(case=case['label'],phasor_parity_tolerance=1e-10,
        phasor_parity_pass=max(errors)<1e-10,coordinate_identity_pass=coordinate_ok,
        resonance_identity_pass=resonance_ok,center_relative_phasor_error=errors[2]))
assert all(assessment[0][k] for k in ('phasor_parity_pass','coordinate_identity_pass','resonance_identity_pass'))
assert not assessment[1]['phasor_parity_pass'] and not assessment[2]['coordinate_identity_pass']
write(OUT/'integration_assessment.json',dict(origin='⟦AI⟧',status='failed',
    criteria='Same physical Hamiltonian and scalar output must agree; software identity tolerance only, not experimental acceptance.',
    cases=assessment))

start='2026-11-04T04:55:27+09:00'; completed=now()
new_id='review-2026-11-04-rydberg-batch-optical-axis'
old_by_id={f['id']:f for f in old['findings']}
findings=[]
def finding(fid,disposition,next_action,evidence,scope,reopen=None):
    previous=old_by_id.get(fid,{})
    index=initial_state['finding_index'].get(fid,{})
    task=canonical.get(fid)
    f=dict(id=fid,checklist_id=fid if task else None,
        canonical_status_observed=task.get('status') if task else None,
        disposition=disposition,first_seen=previous.get('first_seen', '2026-08-12' if fid=='rydberg-doppler-on-pole-residue' else '2026-11-04'),
        first_seen_source=previous.get('first_seen_source',index.get('first_report',REPORT)),
        last_confirmed='2026-11-04',evidence=evidence,next_action=next_action,
        reopen_when=reopen or previous.get('reopen_when','관련 구현·입력·사용 범위 변경 또는 수리 뒤 같은 계약 실패 시 재검토.'),
        origin='⟦AI⟧',adoption=previous.get('adoption',
            'Existing checklist scope/status observed; no checklist edit or new implementation approval.' if task else
            'New provisional AI review observation; not yet adopted into the checklist or approved for implementation.'),
        confirmation_scope=scope)
    findings.append(f)
    return f

finding(new_id,'new',
    '배치 static/configured/scheme 경로가 공통 detuning Hamiltonian·resonance/readout·좌표 export를 사용하도록 정합화. 양쪽 축·비영점 counterpart·기본 parity 검증 후 사용.',
    ['docs/daily_report/2026-11-04_rydberg/integration_probe.json','docs/daily_report/2026-11-04_rydberg/integration_assessment.json','analysis/rydberg_cell_heating/adapters.py'],
    'Fresh deterministic three-fixture physical-point/metadata comparisons; default passes, shifted scans fail. Review priority P1; no product change.')['review_priority']='P1'
finding('review-2026-10-22-rydberg_eit-selective-dephasing-positivity','resolution_candidate',
    '현행 CP 수학적 반례 해소 후보 유지. 기존 negative-state fixture·채널 양성/trace 검증 근거 보존; 새 rate 의미와 실측 power curve 검증은 별도.',
    ['gabes/schemes/rydberg.py','tests/test_rydberg_model_physics.py','docs/daily_report/2026-11-04_rydberg/pytest.xml'],
    '21 fresh model tests include conditional CP and finite channel positivity, invalid rates, legacy negative stationary state, optical-axis parity and exports; not measured laser/transit validation.')
finding('rydberg-readout-validity-guardrails','reconfirmed',
    '기존 P1 관측량별 optical-grid/edge 수렴 기준·고정/단일점 balance reference 계약 검증. 속도축 수렴과 별도.',
    ['docs/daily_report/2026-11-04_rydberg/probe.json','gabes/schemes/rydberg.py','docs/checklist.json'],
    'Fresh four default readouts and receiver instrumentation; numerical_convergence/experimental_visibility NOT_ASSESSED. No repeated optical-grid sweep.')
finding('rydberg-experimental-validation-dataset','reconfirmed',
    'Gold 문헌 조건과 별도로 보유 측정 EIT/RF/PSD·온도/receiver calibration·owner/provenance 확보; fit/holdout 사전 구분.',
    ['references/rydberg_experiments/reference_policy.json','references/rydberg_experiments/papers/ju_2026_preprint.json','analysis/rydberg_cell_heating/example_config.json','docs/checklist.json'],
    'Current local policy/Gold record/example inspected; literature database is not a completed held-out apparatus dataset.')
finding('rydberg-csv-acquisition-order','reconfirmed',
    '기존 EIT/RF 정렬 전 방향·branch 또는 median 병합 손실 경고 보존. PSD 정렬과 구분.',
    ['gabes/rydberg_experimental_csv.py','docs/daily_report/2026-10-28_scheme-3_rydberg-eit_probe_results.json'],
    'Static source identity of CSV loader and numeric helper; prior synthetic 2026-10-28 result historical, not replayed.')
finding('rydberg-effective-opacity-and-sweep-trust','reconfirmed',
    'done 유지. CP 변경 후 현재 원자응답에 opacity 1회 적용·zero response·시나리오·export 계약 검증됨; 옛 원자 배열과 동등성 주장 금지.',
    ['tests/test_rydberg_opacity.py','docs/daily_report/2026-11-04_rydberg/pytest.xml'],
    '24 fresh opacity tests plus integration opacity tests; current propagation contract, not old non-CP spectral parity or measured opacity calibration.')
finding('rydberg-extra-view-dependency-correctness','reconfirmed',
    'done 유지. 기존 extra-view 입력 의존 회귀와 잔여 계층 cache 성능 과제 분리.',
    ['tests/test_rydberg_eit.py','docs/daily_report/2026-11-04_rydberg/pytest.xml','gabes/schemes/rydberg.py'],
    'Fresh declared dependency and detector/cache/sweep tests; no end-to-end browser replay.')
finding('rydberg-layered-sweep-cache','reconfirmed',
    '기존 ready 유지. 원자/IF/propagation/detector 계층 solve 수·exact-key·출력 parity·메모리 경계를 검증하는 기존 과제.',
    ['gabes/schemes/rydberg.py','docs/checklist.json','docs/daily_report/2026-11-04_rydberg/pytest.xml'],
    'Current source/checklist and focused existing cache tests; no new performance benchmark.')
finding('rydberg-finite-if-affine-postprocess','reconfirmed',
    'ready 유지. 새 detuning helper와 CP 채널을 보존하는 조립·detector 배열화 후보; 채택 시 출력 parity와 end-to-end 비용 계측.',
    ['gabes/schemes/rydberg.py','docs/checklist.json'],
    'Source-only bottleneck scope; 2026-10-27 candidate timings not current evidence.')
finding('scheme-current-documentation-and-examples','changed',
    'Rydberg static-only docstring 지적 해소. CP rate·새 optical axes·현재 source/config를 Guide/예제/저장 그림에 반영하는 잔여 범위; 다른 scheme 문서 과제 보존.',
    ['tests/test_rydberg_eit.py','docs/Userguide/GABES_User_Guide_v2.html','gabes/schemes/rydberg.py','docs/checklist.json'],
    'Current docstring correction confirmed; Guide static image already labelled historical/non-validation. Overall checklist ready unchanged.')
finding('rydberg-doppler-on-pole-residue','changed',
    'done 유지. Exact Maxwell adoption/finite-IF resolving-grid 완료 범위 보존; 현 CP 버전에 fresh focal parity. 이전 stage timings 재측정 안 함.',
    ['docs/checklist.json','docs/daily_report/2026-11-04_rydberg/git.json','gabes/affine_doppler.py','docs/daily_report/2026-11-04_rydberg/pytest.xml'],
    'Observed canonical done; fresh affine/grid/guard/fallback and finite-IF refinement tests. 92,532-row and timing claims remain historical.')
findings[-1]['first_seen_source']='docs/daily_report/2026-08-12_scheme-3_rydberg-eit.md'

report=f'''# GABES 일일 리뷰 — 2026-11-04 Rydberg-EIT

## ⟦AI⟧ 선택·변경·판정

- 실제 시작 **2026-11-04 수요일 04:55:27 KST**. 설정상 예정 05:00보다 이른 전달. 실제 요일대로 `rydberg_eit`, `_SCHEMES` 등록 확인. `weekday-v1`; mod 5 미사용.
- 출발점: {link('10-28 같은 주제 보고서',PREVIOUS)}·{link('상태', 'docs/daily_report/review_state.json')}·{link('공식 체크리스트','docs/checklist.json')}. 10-19 취합 cutoff 10-18 유지.
- `main@676a7a61ee693dbde5ca4940c9500f24f965f8cb`. 이전 42개 source 중 **16개 변경**. 현재 관련 **56개 hash 종료까지 동일**. CP/새 광학축 구현은 기존 **미커밋 작업 트리**에 존재; HEAD만의 결과 아님. 사용자 편집 보존.
- **새 지적 1건, 기존 양성성 지적 해결 후보 1건. 리뷰 완료 / 검증 판정 failed.** 기존 시험 **107 passed**이나 새 배치 광학축 계약 probe 실패. scalar 모델의 조건부 실험 reference 가치는 유지; 비영점 detuning 배치 결과는 아래 수정 전 실험 판단에 사용 보류.

근거: {link('동반 evidence',EVIDENCE)}·{link('소스/환경','docs/daily_report/2026-11-04_rydberg/snapshot.json')}·{link('git 한정 조회','docs/daily_report/2026-11-04_rydberg/git.json')}.

## ⟦AI⟧ 주요 관측·다음 행동

### 1. `{new_id}` — new · 잠정 P1 · 체크리스트 미등록

새 probe/coupling 독립 스캔이 {link('본체 공통 Hamiltonian','gabes/schemes/rydberg.py',841)}에는 반영됐으나 {link('배치 adapter','analysis/rydberg_cell_heating/adapters.py',436)}는 기존 `diag(H)/2π = [0,−s,−s,−s−ΔMW]`를 재조립. 실제 식은 `[0,−Δp,−(Δp+Δc),−(Δp+Δc+ΔMW)]`. 양쪽 고정 detuning 0인 기본 fixture는 일치; 비영점 counterpart·coupling 스캔에서 다른 원자 응답을 올바른 static transmission과 조합함.

- **Probe 스캔, 고정 Δc=+12 MHz:** 중심 Δp=−12 MHz에서 configured atomic phasor `−1.49340e−12 + 2.05160e−14 i`, 본체와 같은 물리점 직접 응답 `6.14097e−11 − 9.57783e−11 i` (단위 `(rad/s)⁻¹`, 상대 복소오차 **100.73%**). RF/detector 환산 전 불일치; 그 차이로 감도 오차율을 직접 단정하지 않음.
- **Coupling 스캔, 고정 Δp=+12 MHz:** scheme 모드 배치 export가 최적 probe를 **−10.446 MHz**로 기록. 실제 probe는 **+12 MHz**, −10.446은 최적 coupling 좌표. {link('좌표 복사 위치','analysis/rydberg_cell_heating/adapters.py',367)}.
- {link('Static 요약','analysis/rydberg_cell_heating/adapters.py',165)}도 resonance를 `argmin(abs(x))`로 찾음. 위 probe 스캔에서 `transmission_at_resonance` **0.787767** 대 본체 **0.980231**. 실제 두 광자 중심은 −12 MHz.

재현: {link('probe 소스','docs/daily_report/2026-11-04_rydberg/integration_probe.py')}·{link('세 fixture 원시 결과','docs/daily_report/2026-11-04_rydberg/integration_probe.json')}·{link('명시 계약 판정','docs/daily_report/2026-11-04_rydberg/integration_assessment.json')}. 5개 검색점만 사용한 identity 진단; 수렴된 감도 최적화 아님.

**다음/완료 조건:** 배치 static/configured/scheme 경로를 공통 Hamiltonian·resonance·물리 detuning export와 연결. 기본 parity 유지하며 양쪽 축·비영점 counterpart·필드 단위를 검증. 같은 원인 3증상을 1건으로 추적. 이번 리뷰에서 제품/시험/체크리스트 수정 안 함; 구현 승인은 별도.

### 2. `review-2026-10-22-rydberg_eit-selective-dephasing-positivity` — resolution_candidate

{link('현행 채널','gabes/schemes/rydberg.py',464)}은 `Lopt=√(2γg) diag(0,0,1,1)`, `LRF=√(2γRF) diag(0,0,0,1)`. 따라서 `γ02=γg`, `γ23=γRF`, `γ03=γg+γRF`. 독립 pair 감쇠와 **다른 물리 모델**. 새 {link('21개 모델 시험','tests/test_rydberg_model_physics.py')}에서 conditional CP·유한 채널 양성/trace·이전 음의 정상상태 fixture·양쪽 optical axis/velocity·cache/export 통과. 과거 반례를 현행 결함으로 반복하지 않음. 실측 laser phase/transit rate와의 동일성은 별도 미검증. 공식 checklist 없는 기존 AI 관측의 해결 후보이며 자동 과제 폐쇄 아님.

### 3. `rydberg-readout-validity-guardrails` — reconfirmed · P1 needs_decision

현행 801점·20 °C·probe 6 µW·coupling 30 mW·직경 0.15 mm·셀 50 mm·opacity 0.001에서:

| 모드 | 현재 조건부 출력 | 범위 |
|---|---|---|
| EIT Doppler off / on | 선폭 **1.603494 / 2.878587 MHz** | optical-grid 수렴 인증 아님 |
| AT Doppler off / on, LO 3.7 MHz·IF 40 kHz | total=PSN **13.2413 / 2229.8020 nV/cm/√Hz** | RIN·electronics=0, QE 0.95·path 0.5·RMS dipole 규약 |
| AT Doppler on | static AT **미분해**, finite-IF 감도는 PREDICTED | static split과 finite-IF 유효성 별도 |

{link('원시 입력·출력·receiver 계측','docs/daily_report/2026-11-04_rydberg/probe.json')}. Doppler-off AT 기준광은 **2.519709–2.930825 µW**, 801점마다 signal과 정확히 balance. 고정 기준광 장치 검증 아님. `numerical_convergence`·`experimental_visibility` 모두 `NOT_ASSESSED`. 속도축 refinement PASS가 optical-grid/edge·측정 SNR을 대신하지 않음.

**다음:** 기존 acceptance대로 관측량별 수렴·edge 기준과 고정/단일점 balance reference 계약 검증. 임의 투과율 cutoff로 실험 불가 판정하거나 static AT 미분해만으로 finite-IF 감도를 일괄 기각하지 않음.

## ⟦AI⟧ 완료 유지·정정·잔여 항목

| ID | 이번 판정 |
|---|---|
| `rydberg-doppler-on-pole-residue` | **changed, done 유지**. 11-03 `e251527`에서 exact Maxwell 평균·finite-IF resolving grid 채택. 현행 affine/guard/fallback·velocity refinement 시험 PASS. 기존 완료를 미완료로 재지적 안 함. |
| `rydberg-effective-opacity-and-sweep-trust` | **reconfirmed, done**. 현재 CP 원자응답에서 opacity 1회 적용·zero response·시나리오·표시/export 시험 PASS. 옛 non-CP 원자 배열 불변 주장은 하지 않음. |
| `rydberg-extra-view-dependency-correctness` | **reconfirmed, done**. 입력 선언·cache·sweep 시험 통과. 잔여 계층 성능 과제와 구분. |
| `scheme-current-documentation-and-examples` | **changed, ready 유지**. 시험 docstring의 static-only 설명 이미 수정. Guide 그림도 정적 예시/finite-IF 검증 아님을 명시. CP·새 축·현행 source/config 예제 동기화는 별도 잔여 범위. |
| `rydberg-experimental-validation-dataset` | **reconfirmed, blocked_external**. {link('11-03 사용자 Gold 지정','references/rydberg_experiments/reference_policy.json')}과 {link('26편 문헌 DB','references/rydberg_experiments/README.md')}는 새 근거. 보유 측정 EIT/RF/PSD·receiver·온도 metadata를 갖춘 held-out 자료 완료와 다름. 예제 실측 입력 PENDING 유지. |
| `rydberg-csv-acquisition-order` | **reconfirmed, ready**. loader·numeric helper hash 동일; branch/방향 손실 구조 유지. 10-28 synthetic 숫자는 historical, 재실행 안 함. |
| `rydberg-layered-sweep-cache` / `rydberg-finite-if-affine-postprocess` | **reconfirmed, ready**. 기존 계층 재사용·조립/배열화 범위 유지. 새 benchmark·속도 후보 없음. |

중요 정정: exact-Maxwell 완료 기록의 **Doppler-on 1557.57 nV/cm/√Hz**, “Doppler-off bitwise”는 **CP 변경 전 소스 비교**에 한정. 현재 작업 트리는 `rydberg-eit-v8-cp-optical-detunings`; 오늘 **2229.80**을 해당 과거 수치와 동일시하지 않음. 이번 probe는 변경별 감도 기여를 분리한 실험이 아님. Gold 감도와 가까운 Doppler-off **13.24**도 zero technical noise·미보정 opacity·Rabi/transit 가정의 조건부 결과; 독립 실험 재현 근거로 승격하지 않음.

## ⟦AI⟧ 검증 범위·생략

- **Fresh 100 passed, 23.02 s, 2 deprecation warnings:** `python -m pytest -q tests/test_rydberg_model_physics.py tests/test_affine_doppler.py tests/test_rydberg_eit.py tests/test_rydberg_opacity.py --durations=8` (+ JUnit). 모듈별 21/18/37/24. {link('명령/실행시간','docs/daily_report/2026-11-04_rydberg/checks.json')}·{link('로그','docs/daily_report/2026-11-04_rydberg/pytest.txt')}. 경고는 export의 `datetime.utcnow()`; 물리 실패로 분류 안 함.
- **Fresh 7 passed, 28.02 s:** adapter integration 전체 + workflow의 calibrated finite-IF 시험 1개. {link('명령','docs/daily_report/2026-11-04_rydberg/integration_checks.json')}·{link('로그','docs/daily_report/2026-11-04_rydberg/integration_pytest.txt')}. **기존 시험은 새 비영점 광학축 probe 실패를 검출하지 못함.** 전체 판정 failed를 green으로 덮지 않음.
- 환경: 원본 G: 작업공간, Python 3.12.6·NumPy 2.2.0·SciPy 1.18.0·Numba 0.65.1, `MPLBACKEND=Agg`. threading/backend는 snapshot에 저장. 시간은 시험 수행시간이며 cold/warm 성능 benchmark 아님.
- **Reused:** source 동일한 CSV 방향 손실의 정적 구현 관측. **Historical:** 10-27 optical-grid/positivity/속도 수치, 10-28 receiver 숫자, exact-Maxwell commit의 92,532행·시간·모든 slider 검증. 날짜·소스 범위 유지; 오늘 broad PASS로 재사용 안 함.
- 전체 pytest·GC audit/campaign·sealed 연구·대형 sweep/fit·ARC 재생성·원격 issue/신규 문헌 검색·브라우저 검토 생략. 이번 세션 제품 변경 없음; 새 CP·광학축과 영향받는 배치 caller를 집중 확인. 전체 엔진/실험/주간 UI PASS 주장 없음. 한정 TODO/FIXME 검색 결과 없음은 원격 issue 부재 근거 아님.
- 재현 뒤 제품 수리는 수행하지 않음. **미해결 필수 확인:** 새 batch optical-axis 계약은 수정 후 실패 fixture·기존 parity 재검증 필요. 기존 관측량별 수렴·실측 검증 과제도 열림. 이번 리뷰의 링크/hash·상태 게시 외 추가 필수 검증 보류 없음.

다음 우선 행동: **배치 새 광학축 연결 수정**, 이후 기존 readout 수렴·receiver 계약. 미채택 AI 검토 제안; 이번 리뷰로 사용자 목표·우선순위·구현 승인 추가 안 함. 다음 실제 목요일은 `magneto`; 누락이 생겨도 배정 이동 없음.
'''
(ROOT/REPORT).write_text(report,encoding='utf-8')

link_errors=[]; link_count=0
for target in re.findall(r'\]\(<([^>]+)>\)',report):
    path_text=re.sub(r':\d+$','',target)
    path=Path(path_text)
    link_count+=1
    if not path.is_file() and path!=ROOT/EVIDENCE: link_errors.append(target)
assert not link_errors,link_errors
prior_hashes=old['source_state']['hashes_at_closure']
changed_since=[p for p,h in prior_hashes.items() if hashes.get(p)!=h]
assert len(changed_since)==16
evidence=dict(schema_version=1,origin='⟦AI⟧',workflow='weekday-v1',local_date='2026-11-04',
    timezone='Asia/Seoul',weekday='Wednesday',mode='scheme',topic=TOPIC,
    started_at=start,completed_at=completed,scheduled_for='2026-11-04T05:00:00+09:00',trigger='scheduled',
    timing_note='Actual delivery began 4m33s before configured 05:00 slot; actual Wednesday routing retained.',
    outcome='completed',validation_status='failed',
    validation_scope='107 fresh existing focused tests passed; new deterministic batch optical-axis physical-point/metadata contract probe failed. No full-suite, experimental calibration, global convergence, browser, benchmark or GC PASS.',
    previous_report=PREVIOUS,report_path=REPORT,report_sha256=sha(ROOT/REPORT),
    source_state=dict(head=git['head']['stdout'].strip(),branch=git['branch']['stdout'].strip(),
        hashes_at_start=hashes,hashes_at_closure=closure,changed_during_review=changed_during,
        previous_manifest_count=len(prior_hashes),changed_since_previous_manifest=changed_since,
        environment=snapshot['environment'],registry=snapshot['registry'],git_evidence='docs/daily_report/2026-11-04_rydberg/git.json',
        limitations=['Bounded tracked diff confirms existing uncommitted Rydberg code/test and README/checklist edits; whole-worktree clean/staged/untracked inventory not claimed.',
            'New model-physics test is absent from git ls-files at review time; actual bytes hashed and tested.',
            'Prior broad numerical and performance evidence not promoted across changed sources/runtime scope.']),
    changes=['Shared affine exact-Maxwell and resolving finite-IF velocity grid completed since prior same-topic review.',
             'Existing uncommitted CP nested dephasing, independent optical scan axes/fixed counterpart detunings and live probe-Rabi anchor.',
             'User temporarily adopted Ju Gold and local literature DB on 2026-11-03; no held-out apparatus validation inferred.',
             'Static-only test docstring corrected; current opacity tests now cover propagation with CP physics.'],
    findings=findings,checks=[
        dict(id='focused-existing-tests',evidence_kind='fresh',original_date='2026-11-04',
            command=checks['command'],path='docs/daily_report/2026-11-04_rydberg/pytest.txt',
            result='100 passed, 0 skipped/failed, 2 deprecation warnings, 23.02 s',status='passed',
            scope='21 model,18 affine,37 scheme,24 opacity software fixtures; no measured-data acceptance.'),
        dict(id='existing-adapter-tests',evidence_kind='fresh',original_date='2026-11-04',command=integration_checks['command'],
            path='docs/daily_report/2026-11-04_rydberg/integration_pytest.txt',result='7 passed, 28.02 s',status='passed',
            scope='Existing default thermal/dephasing/opacity handoff plus one configured workflow; lacks nonzero-counterpart optical-axis contract.'),
        dict(id='new-batch-optical-axis-contract',evidence_kind='fresh',original_date='2026-11-04',status='failed',
            command='python docs/daily_report/2026-11-04_rydberg/integration_probe.py',
            path='docs/daily_report/2026-11-04_rydberg/integration_probe.json',
            parameters='20 C, AT electrometry, Doppler off, five detuning search points; default, probe scan/fixed coupling12 MHz, coupling scan/fixed probe12 MHz.',
            observable='atomic phasor same physical-point parity, optimum physical-probe metadata, resonance transmission identity',
            units='phasor (rad/s)^-1; detuning MHz; transmission dimensionless',
            result=assessment,scope='Default matches; changed-axis/fixed-counterpart contract fails before detector scaling. Not a converged sensitivity/error-percentage or experimental benchmark.'),
        dict(id='current-default-readouts',evidence_kind='fresh',original_date='2026-11-04',status='observed',
            command='python docs/daily_report/2026-11-04_rydberg/run_review.py --probe-only',path='docs/daily_report/2026-11-04_rydberg/probe.json',
            parameters='Four default EIT/AT x Doppler cases; full parameters in JSON.',
            observable='linewidth,conditional field ASD,opacity/visibility/convergence status,receiver powers',
            units='MHz,nV/cm/sqrt(Hz),W,dimensionless transmission',
            result='Current CP defaults, pointwise balanced receiver and NOT_ASSESSED convergence/visibility recorded.',
            scope='No optical-grid refinement or measured apparatus validation.'),
        dict(id='csv-source-observation',evidence_kind='reused',original_date='2026-10-28',
            path='docs/daily_report/2026-10-28_scheme-3_rydberg-eit.md',
            parameters='Unchanged rydberg_experimental_csv.py and _csv_numeric.py source.',
            result='Static direction/branch-loss implementation remains; original synthetic quantitative result historical.',
            scope='Source-only observation; no historical numerical/benchmark PASS adopted.'),
        dict(id='prior-scans-and-campaigns',evidence_kind='historical',original_date='2026-10-27',
            path='docs/daily_report/2026-10-27_scheme-3_rydberg-eit_evidence.json',
            result='Old grid,non-CP transient and candidate timing evidence retains original source/assumption scope.',
            scope='Not reused as current validation.'),
        dict(id='exact-maxwell-completion',evidence_kind='historical',original_date='2026-11-03',
            path='docs/daily_report/2026-11-04_rydberg/git.json',parameters='e251527 commit and checklist scope',
            result='Done/explicit user exact-Maxwell decision observed; 1557.57 ASD and off-bitwise were pre-CP comparison only.',
            scope='Current tests confirm bounded parity; no full slider sweep or stage benchmark rerun.')],
    intentional_omissions=['Full pytest and unrelated GC audit/sealed campaigns: review makes no product/shared-engine edits and makes no global-release claim.',
        'Large sweeps/fits/ARC regeneration/current benchmark/browser/remote issues and fresh literature search omitted; topic is bounded Wednesday code/evidence review.'],
    deferred_required_validation=[dict(id=new_id,scope='Repair then replay both shifted axes, configured atomic phasors, static resonance and scheme physical-detuning exports with default parity.',status='unresolved_existing_product_contract_failure')],
    open_project_validation=['Observable-specific optical-grid/edge and fixed-reference receiver acceptance.',
        'Held-out EIT/RF/PSD/thermal/detector calibration; literature conditions do not replace it.'],
    new_actionable_findings=1,resolution_candidates=1,
    rejected_interpretations=['Old selective dephasing counterexample is not an unchanged current-model defect.',
        'Doppler completed optimization is not an open task; pre-CP bitwise/value evidence is not current CP spectral parity.',
        'Static-only test docstring and unlabeled finite-IF Guide criticism no longer applies.',
        '107 existing test passes do not negate the fresh batch optical-axis failure.'],
    probe_execution_note='Initial review harness JSON serialization failed on DetectorNoiseBudget after the 100 tests passed. Dataclass serialization corrected in review script; only --probe-only rerun. Product/test bytes unchanged; no duplicate test run.',
    state_update=dict(pre_write_sha256=initial_state_sha,method='Reread and merge latest state; compare immediate pre-replace SHA, write sibling temp then os.replace; preserve unknown fields and unrelated topic pointers.'),
    artifact_validation=dict(evidence_kind='fresh',local_link_count=link_count,local_link_errors=[],strict_json=True,
        report_sha256_verified=True,source_hash_count=len(closure),source_hash_errors=[]))
artifacts=[p for p in OUT.iterdir() if p.is_file() and p.name not in ('publication.json',)]
evidence['artifacts']={p.relative_to(ROOT).as_posix():sha(p) for p in artifacts}
write(ROOT/EVIDENCE,evidence)
assert read(ROOT/EVIDENCE)['report_sha256']==sha(ROOT/REPORT)
assert (ROOT/EVIDENCE).is_file()

entry={k:evidence[k] for k in ('local_date','topic','mode','workflow','outcome','validation_status','validation_scope','report_path','report_sha256','completed_at')}
entry.update(evidence_path=EVIDENCE,evidence_sha256=sha(ROOT/EVIDENCE))
for attempt in range(5):
    expected=sha(STATE); current=read(STATE)
    current['latest_completed'][TOPIC]=entry
    current['latest_attempt'][TOPIC]=entry.copy()
    for f in findings:
        prior=current['finding_index'].get(f['id'],{})
        index={**prior,'checklist_id':f['checklist_id'],
            'first_report':prior.get('first_report',f['first_seen_source']),
            'latest_report':REPORT,'latest_observation_report':REPORT,
            'disposition':f['disposition'],'next_action':f['next_action'],'reopen_when':f['reopen_when'],
            'origin':'⟦AI⟧','adoption':f['adoption'],'confirmation_scope':f['confirmation_scope']}
        if f['id']=='scheme-current-documentation-and-examples':
            index['next_action']='각 scheme 문서·그림·실행 예제의 현재 모델/입력/provenance 정합화; scheme별 관측 범위 참조.'
            index.setdefault('topic_observations',{})['lambda']=dict(report=prior.get('latest_observation_report'),next_action=prior.get('next_action'))
            index['topic_observations']['rydberg_eit']=dict(report=REPORT,next_action=f['next_action'])
        current['finding_index'][f['id']]=index
    current['updated_at']=completed
    tmp=STATE.with_name('review_state.2026-11-04-rydberg.tmp')
    write(tmp,current)
    if sha(STATE)!=expected:
        tmp.unlink(); continue
    os.replace(tmp,STATE)
    break
else: raise RuntimeError('Concurrent state updates prevented safe publication.')
saved=read(STATE)
assert saved['latest_completed'][TOPIC]==entry
for topic,value in initial_state['latest_completed'].items():
    if topic!=TOPIC: assert saved['latest_completed'][topic]==value
final_hashes={p:sha(ROOT/p) for p in hashes}
assert final_hashes==closure,'Source changed during publication; revise affected claims before completion.'
receipt=dict(origin='⟦AI⟧',published_at=now(),report_path=REPORT,report_sha256=sha(ROOT/REPORT),
    evidence_path=EVIDENCE,evidence_sha256=sha(ROOT/EVIDENCE),state_sha256=sha(STATE),
    numerical_validation_status='failed',existing_tests_passed=107,new_actionable_findings=1,
    report_links_checked=link_count,source_hashes_stable=len(final_hashes),
    unrelated_completed_topic_pointers_preserved=True,unknown_state_fields_preserved=True,
    artifact_validation_status='passed')
write(OUT/'publication.json',receipt)
print(json.dumps(receipt,ensure_ascii=False,indent=2))
