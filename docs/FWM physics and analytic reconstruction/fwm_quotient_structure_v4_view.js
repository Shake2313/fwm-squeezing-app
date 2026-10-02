/* Progressive enhancement: all scientific records live in the static document. */
(() => {
'use strict';
const data = JSON.parse(document.getElementById('atlas-data').textContent);
const groupById = new Map(data.groups.map(g => [g.id, g]));
const nodeById = new Map(data.groups.flatMap(g => g.nodes.map(n => [n.id, n])));
const kinds = {input:'물리 입력',ready:'생산 구현',reference:'별도 참조',analytic:'조건부 해석',test:'대수·수치 시험',need:'전체 연결 미완료',conditional:'조건부 연구 구현',definition:'정의·명제'};
const edgeKinds = {data:'입력 전달',alternative:'동일 문제의 대안',limit:'가정하의 축약',validation:'검증 대상',analogy:'형태의 유사성',logic:'논리 전제',equivalence:'정확한 표현 동치',refinement:'수치 세분',evidence:'범위가 있는 증거'};
const numbers = {input:'01',atomic:'02',mean:'03',gain:'04',micro:'05',covariance:'06',detector:'07',output:'08',reference:'R',structure:'L',microscopic:'Q',transport:'T'};
const groupFormula = {input:'Δ · δ · Ω_SA · V_in · η',atomic:'{ρ_n} → ⟨χ̄⟩ᵥ',mean:'α_out = T_cl α_in',gain:'G_p · G_c · TJT†',reference:'Tr δρ = 0 · χ_ref(Ω)',micro:'M_q(Ω) · D(Ω)',covariance:'V_out = T_q V_in T_q† + W',detector:'V_R,det · α_det → PSD / SQL',output:'S_−(Ω) · 10 log₁₀ S_−',structure:'A ≼ B · A ∼ B · [A] ≤ [B]',microscopic:'ℒ → A, D → QRT → field',transport:'path → ensemble → nonlocal field'};
const esc = x => String(x).replace(/[&<>"']/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
const svgNS='http://www.w3.org/2000/svg';
const makeSVG=(tag,attrs={})=>{const el=document.createElementNS(svgNS,tag);Object.entries(attrs).forEach(([k,v])=>el.setAttribute(k,v));return el;};
const app=document.createElement('section'); app.id='explorer'; app.setAttribute('aria-label','FWM 계산 구조 탐색');
app.innerHTML=`<header class="atlas-header"><div class="atlas-brand"><span class="atlas-mark" aria-hidden="true">χ</span><div><strong>GABES <span class="atlas-version">FWM</span></strong><small>THEORY · ORDER · EVIDENCE</small></div></div><nav class="atlas-actions" aria-label="보기와 내보내기"><button id="atlas-reader" aria-pressed="false">문서로 읽기</button><button id="atlas-print">인쇄</button><button id="atlas-download" class="primary">HTML 저장 ↓</button></nav></header>
<div class="atlas-status"><span class="status-dot" aria-hidden="true"></span><span><strong>조건부 microscopic 이론 구현 · 실제 장치 예측 미인증</strong> · 정확한 동치, 수치 수렴, 실험 검증을 구별합니다.</span><a href="#o-claim">구현 범위 ↗</a></div>
<p class="atlas-reader-note">모든 그룹의 수식·조건·연결을 펼친 문서입니다. ‘다이어그램으로’ 버튼으로 탐색 화면에 돌아갑니다.</p>
<div class="atlas-workspace"><nav class="atlas-nav" aria-label="계산 그룹"><p class="atlas-nav-heading">CALCULATION MAP</p><button data-open-group="overview" aria-current="true"><span class="nav-num">◈</span><span>전체 흐름</span></button><div class="atlas-nav-divider"></div>${data.groups.filter(g=>g.id!=='reference').map(g=>`<button data-open-group="${g.id}"><span class="nav-num">${numbers[g.id]}</span><span>${esc(g.title)}</span></button>`).join('')}<div class="atlas-nav-divider"></div><button data-open-group="reference"><span class="nav-num">R</span><span>pump-only 참조</span></button><div class="atlas-nav-foot">${nodeById.size} nodes · ${data.edges.length} relations<br>물리·논리·수치·실험의 구분<br>배지는 실행 결과가 아닙니다.</div></nav>
<div class="atlas-main"><div class="canvas-head"><div class="canvas-breadcrumb"><button data-open-group="overview">전체 지도</button><span id="atlas-crumb">/ 계산의 두 경로</span></div><div class="canvas-heading"><h1 id="atlas-title">FWM 이론의 구조와 성립 조건</h1><small id="atlas-count">${data.groups.length} GROUPS</small></div><p class="canvas-subtitle" id="atlas-subtitle">그룹을 열어 계산의 내부 구조와 성립 조건을 읽습니다.</p></div>
<div class="canvas-tools"><div class="atlas-segment" role="group" aria-label="노드 표시"><button data-atlas-view="flow" aria-pressed="true">물리 흐름</button><button data-atlas-view="math" aria-pressed="false">수식</button><button data-atlas-view="evidence" aria-pressed="false">검증</button></div><div class="zoom-tools"><button id="zoom-out" aria-label="축소">−</button><output class="zoom-value" id="zoom-value" aria-label="확대 비율">100%</output><button id="zoom-in" aria-label="확대">+</button><button id="zoom-fit" title="전체 연결도 맞춤 (0)">전체 맞춤</button></div></div>
<div class="atlas-viewport" id="atlas-viewport" tabindex="0" role="region" aria-label="연결도. 드래그 또는 방향키로 이동, 더하기와 빼기로 확대 및 축소, 0으로 전체 맞춤."><div class="atlas-scene" id="atlas-scene"></div><svg class="atlas-minimap" id="atlas-minimap" aria-hidden="true"></svg></div>
<div class="canvas-bottom"><span><i></i>입력 전달</span><span><i class="pending"></i>◇ 미완성 연결</span><span><i class="validation"></i>검증·대안·축약</span><span class="canvas-hint">드래그 이동 · 휠 확대</span></div></div>
<aside class="atlas-inspector" id="atlas-inspector" aria-label="선택한 노드의 수식과 조건" tabindex="-1"></aside></div>`;
document.body.insertBefore(app,document.body.firstChild);
document.body.classList.add('explorer-active','atlas-overview');
const viewport=document.getElementById('atlas-viewport'),scene=document.getElementById('atlas-scene'),panel=document.getElementById('atlas-inspector'),minimap=document.getElementById('atlas-minimap');
const state={group:'overview',node:null,view:'flow',x:0,y:0,scale:1,width:1140,height:730,positions:new Map(),links:[],reader:false};
function hash(id){if(location.hash!=='#'+id)history.pushState(null,'','#'+id);}
function resetReader(value){
 state.reader=value;document.body.classList.toggle('reader-mode',value);
 const b=document.getElementById('atlas-reader');b.textContent=value?'다이어그램으로':'문서로 읽기';b.setAttribute('aria-pressed',String(value));
 if(value){document.querySelectorAll('.shell details').forEach(el=>el.open=true);document.body.classList.add('all-expanded');}
 else{document.body.classList.remove('all-expanded');requestAnimationFrame(()=>fit(false));}
}
function paintCamera(){
 const low=Math.min(.25,Math.min(viewport.clientWidth/state.width,viewport.clientHeight/state.height)*.85);
 state.scale=Math.max(low,Math.min(2.25,state.scale));
 const width=state.width*state.scale,height=state.height*state.scale;
 state.x=width<viewport.clientWidth?(viewport.clientWidth-width)/2:Math.max(viewport.clientWidth-width-70,Math.min(70,state.x));
 state.y=height<viewport.clientHeight?(viewport.clientHeight-height)/2:Math.max(viewport.clientHeight-height-70,Math.min(70,state.y));
 scene.style.transform=`translate(${state.x}px,${state.y}px) scale(${state.scale})`;
 document.getElementById('zoom-value').value=Math.round(state.scale*100)+'%';
 drawMinimap();
}
function fit(all=true){
 const w=viewport.clientWidth||800,h=viewport.clientHeight||550;
 const full=Math.min((w-40)/state.width,(h-32)/state.height,1.15);
 state.scale=all?full:Math.min((w-35)/state.width,Math.max(full,.78),1.05);
 state.x=(w-state.width*state.scale)/2;state.y=all?(h-state.height*state.scale)/2:20;paintCamera();
}
function zoom(factor,cx=viewport.clientWidth/2,cy=viewport.clientHeight/2){
 const old=state.scale;state.scale=Math.max(.15,Math.min(2.25,state.scale*factor));
 state.x=cx-(cx-state.x)*state.scale/old;state.y=cy-(cy-state.y)*state.scale/old;paintCamera();
}
function drawMinimap(){
 minimap.replaceChildren();minimap.setAttribute('viewBox',`0 0 ${state.width} ${state.height}`);
 state.positions.forEach((p,id)=>minimap.append(makeSVG('rect',{x:p.x,y:p.y,width:p.w,height:p.h,class:'map-node'+(id===state.node?' selected':'')})));
 minimap.append(makeSVG('rect',{x:-state.x/state.scale,y:-state.y/state.scale,width:viewport.clientWidth/state.scale,height:viewport.clientHeight/state.scale,class:'map-window'}));
}
const overviewXY={input:[35,78],atomic:[305,78],mean:[575,78],gain:[845,78],reference:[305,330],micro:[575,330],covariance:[845,330],detector:[845,576],output:[575,576],structure:[35,830],microscopic:[305,830],transport:[575,830]};
const localXY={
 input:[[0,0],[1,0],[0,1],[1,1],[0,2],[1,2]],
 atomic:[[0,0],[0,1],[1,1],[0,2],[0,3],[1,3],[0,4]],
 reference:[[0,0],[0,1],[1,1],[0,2],[0,3],[1,3],[1,4],[0,4],[1,5],[0,5],[0,6]],
 mean:[[0,0],[1,0],[0,1],[0,2],[1,2],[0,3],[1,3]],
 gain:[[0,0],[1,0],[0,1]],
 micro:[[0,0],[0,1],[1,1],[0,2],[1,2],[0,3]],
 covariance:[[0,0],[0,1],[1,1],[0,2],[0,3],[1,3],[0,4],[1,4],[0,5],[0,6]],
 detector:[[0,0],[1,0],[0,1],[1,1],[0,2],[0,4],[0,5],[0,6],[1,5]],
 output:[[0,0],[1,0]]};
function layout(group){
 state.positions.clear();
 if(group==='overview'){
  state.width=1140;state.height=1030;
  Object.entries(overviewXY).forEach(([id,[x,y]])=>state.positions.set(id,{x,y,w:230,h:140}));
  state.links=data.spine.map(e=>({...e}));
  // Reference links summarize real cross-group edges without inventing a
  // production prerequisite. Full typed relations remain in node records.
  for(const [s,t,kind] of [['atomic','reference','data'],['reference','micro','data']]){
   const actual=data.edges.filter(e=>nodeById.get(e.source).group===s&&nodeById.get(e.target).group===t&&e.kind===kind);
   if(actual.length)state.links.push({source:s,target:t,kind,label:s==='atomic'?'공통 원자 조립':'응답 재료',pending:actual.some(e=>e.pending)});
  }
 }else{
  const g=groupById.get(group),coords=localXY[group]||g.nodes.map(n=>[n.xy[0]-1,n.xy[1]-1]);state.width=740;
  g.nodes.forEach((n,i)=>{const [x,y]=coords[i];state.positions.set(n.id,{x:55+x*355,y:60+y*235,w:275,h:170});});
  state.height=Math.max(...[...state.positions.values()].map(p=>p.y+p.h))+70;
  state.links=data.edges.filter(e=>nodeById.get(e.source).group===group&&nodeById.get(e.target).group===group);
 }
}
function routePath(a,b,index,overview=false){
 const x1=a.x+a.w/2,x2=b.x+b.w/2,y1=a.y+a.h/2,y2=b.y+b.h/2;
 if(overview){
  if(Math.abs(a.y-b.y)<30){const right=b.x>a.x;return {d:`M${right?a.x+a.w:a.x} ${y1} H${right?b.x:b.x+b.w}`,x:(right?a.x+a.w+b.x:a.x+b.x+b.w)/2,y:y1-12};}
  if(a.x===b.x)return {d:`M${x1} ${a.y+a.h} V${b.y}`,x:x1+8,y:(a.y+a.h+b.y)/2};
  if(a.x===575&&a.y===78&&b.x===845){const lane=1120;return{d:`M${x1} ${a.y+a.h} V246 H${lane} V${y2} H${b.x+b.w}`,x:lane-10,y:510,anchor:'end'};}
  const mid=(a.y+a.h+b.y)/2;return{d:`M${x1} ${a.y+a.h} V${mid} H${x2} V${b.y}`,x:(x1+x2)/2,y:mid-10};
 }
 const down=b.y-(a.y+a.h);
 if(down>=0&&down<140){const mid=(a.y+a.h+b.y)/2+((index%3)-1)*10;return{d:`M${x1} ${a.y+a.h} V${mid} H${x2} V${b.y}`,x:(x1+x2)/2,y:mid-9};}
 if(Math.abs(a.y-b.y)<30){const right=b.x>a.x;return{d:`M${right?a.x+a.w:a.x} ${y1} H${right?b.x:b.x+b.w}`,x:(right?a.x+a.w+b.x:a.x+b.x+b.w)/2,y:y1-10};}
 const right=(a.x+b.x)>state.width/2,side=right?state.width-12-(index%4)*7:12+(index%4)*7;
 const ax=right?a.x+a.w:a.x,bx=right?b.x+b.w:b.x;
 return{d:`M${ax} ${y1+20} H${side} V${y2-20} H${bx}`,x:side+(right?-6:6),y:(y1+y2)/2,anchor:right?'end':'start'};
}
function renderGraph(){
 scene.replaceChildren();scene.style.width=state.width+'px';scene.style.height=state.height+'px';
 const svg=makeSVG('svg',{class:'atlas-edges',width:state.width,height:state.height,'aria-hidden':'true'});
 const defs=makeSVG('defs');
 for(const [id,color] of [['normal','#8298b4'],['pending','#b28948'],['active','#315fa8'],['validation','#8776ac']]){
  const marker=makeSVG('marker',{id:'atlas-arrow-'+id,viewBox:'0 0 10 10',refX:9,refY:5,markerWidth:6,markerHeight:6,orient:'auto'});marker.append(makeSVG('path',{d:'M1 1L9 5L1 9',fill:'none',stroke:color,'stroke-width':1.6}));defs.append(marker);
 }svg.append(defs);scene.append(svg);
 if(state.group==='overview'){
  [['MEAN FIELD · 현재 생산 경로',35,25,''],['QUANTUM NOISE · 조건부 모델',575,273,'국소 M_q·D 구현 / full thermal 연결 미완료'],['DETECTION · 두 경로의 합류',575,525,''],['LOGIC · MICROSCOPIC · TRANSPORT',35,777,'서로 다른 관계의 정의와 적용 범위']].forEach(([label,x,y,small])=>{const l=document.createElement('div');l.className='canvas-lane';l.style.left=x+'px';l.style.top=y+'px';l.innerHTML=esc(label)+(small?'<small>'+esc(small)+'</small>':'');scene.append(l);});
 }else if(state.group==='detector'){
  const l=document.createElement('div');l.className='canvas-lane rule';l.style.cssText='left:55px;top:875px;width:630px';l.textContent='별도 이상 시험 · TMSV 광자수 통계';scene.append(l);
 }
 state.links.forEach((e,i)=>{
  const a=state.positions.get(e.source),b=state.positions.get(e.target);if(!a||!b)return;
  const r=routePath(a,b,i,state.group==='overview');
  const p=makeSVG('path',{d:r.d,class:`atlas-edge ${e.kind}${e.pending?' pending':''}`,'data-edge':i,'marker-end':`url(#atlas-arrow-${e.pending?'pending':e.kind==='validation'?'validation':'normal'})`});
  const title=makeSVG('title');title.textContent=edgeKinds[e.kind]+': '+e.label+(e.pending?' · 미완성':'');p.append(title);svg.append(p);
  const label=makeSVG('text',{x:r.x,y:r.y,'text-anchor':r.anchor||'middle',class:'atlas-edge-label'+(e.pending?' pending':''),'data-edge-label':i});
  label.textContent=state.group==='overview'?(e.pending?'◇ ':'')+e.label:String(i+1).padStart(2,'0')+(e.pending?' ◇':'');svg.append(label);
 });
 state.positions.forEach((p,id)=>{
  const isGroup=state.group==='overview',n=isGroup?groupById.get(id):nodeById.get(id),button=document.createElement('button');
  button.className='atlas-node kind-'+n.kind;button.dataset.canvasNode=id;button.setAttribute('aria-pressed',String(id===state.node));
  button.style.cssText=`left:${p.x}px;top:${p.y}px;width:${p.w}px;height:${p.h}px`;
  const description=isGroup?(state.view==='math'?groupFormula[id]:state.view==='evidence'?n.evidence:n.subtitle):state.view==='math'?n.compact:state.view==='evidence'?n.evidence:n.out;
  button.innerHTML=`<span class="node-top"><span class="node-key">${isGroup?numbers[id]:esc(id)}</span><span class="node-status">${esc(kinds[n.kind])}</span></span><strong>${esc(n.title)}</strong><span class="node-output">${esc(description)}</span>${isGroup?`<span class="node-foot">${n.nodes.length}개 내부 노드 <span aria-hidden="true">↗</span></span>`:''}`;
  button.addEventListener('click',()=>isGroup?openGroup(id):selectNode(id,true));scene.append(button);
 });
 highlight();drawMinimap();
}
function highlight(){
 const connected=new Set(state.node?[state.node]:[]),activeEdges=new Set();
 if(state.node){
  // Only one-hop incident relations are emphasized. Validation and analogy do
  // not become transitive computational dependencies.
  state.links.forEach((e,i)=>{if(e.source===state.node||e.target===state.node){connected.add(e.source);connected.add(e.target);activeEdges.add(i);}});
 }
 scene.querySelectorAll('[data-canvas-node]').forEach(el=>{el.classList.toggle('is-dim',!!state.node&&!connected.has(el.dataset.canvasNode));el.setAttribute('aria-pressed',String(el.dataset.canvasNode===state.node));});
 scene.querySelectorAll('[data-edge]').forEach(el=>{const i=Number(el.dataset.edge),e=state.links[i],active=activeEdges.has(i);el.classList.toggle('is-active',active);el.classList.toggle('is-dim',!!state.node&&!active);el.setAttribute('marker-end',`url(#atlas-arrow-${e.pending?'pending':e.kind==='validation'?'validation':active?'active':'normal'})`);});
 scene.querySelectorAll('[data-edge-label]').forEach(el=>{const active=activeEdges.has(Number(el.dataset.edgeLabel));el.classList.toggle('is-active',active);el.classList.toggle('is-dim',!!state.node&&!active);});
}
function portButtons(g){
 const records=data.edges.filter(e=>(nodeById.get(e.source).group===g.id)!==(nodeById.get(e.target).group===g.id));
 const seen=new Set();return records.map(e=>{
  const incoming=nodeById.get(e.target).group===g.id,other=nodeById.get(incoming?e.source:e.target),key=other.group+incoming+e.kind;
  if(seen.has(key))return '';seen.add(key);
  return `<button data-open-node="${other.id}" class="${e.pending?'pending':''}" title="${esc(edgeKinds[e.kind]+': '+e.label)}">${incoming?'← ':'→ '}${esc(groupById.get(other.group).title)}${e.pending?' ◇':''}<br><small>${esc(edgeKinds[e.kind])}</small></button>`;
 }).join('');
}
function groupPanel(g){
 const local=state.links.length;
 panel.innerHTML=`<p class="inspector-kicker">GROUP ${numbers[g.id]} / ${esc(g.part)}</p><h2>${esc(g.title)}</h2><span class="status ${g.kind}">${esc(kinds[g.kind])}</span><p class="inspector-summary">${esc(g.brief)}</p><div class="inspector-stats"><div><span class="inspector-number">${g.nodes.length}</span><small>내부 노드</small></div><div><span class="inspector-number">${local}</span><small>내부 연결</small></div></div><div class="inspector-section"><h3>다른 그룹과의 관계</h3><div class="atlas-ports">${portButtons(g)}</div></div><div class="inspector-section"><h3>내부 노드 바로가기</h3><div class="atlas-ports">${g.nodes.map(n=>`<button data-open-node="${n.id}">${esc(n.title)}</button>`).join('')}</div></div><div class="inspector-section"><h3>검증 범위</h3><p class="inspector-summary">${esc(g.evidence)}</p></div><button class="view-record" data-read-record="group-${g.id}">이 그룹을 문서로 읽기 ↗</button>`;
 panel.scrollTop=0;
}
// Native MathML keeps representative equations typeset and entirely offline.
// The original full equation text remains available directly below each one.
const mi=x=>`<mi>${x}</mi>`,mo=x=>`<mo>${x}</mo>`,mn=x=>`<mn>${x}</mn>`,row=(...x)=>`<mrow>${x.join('')}</mrow>`,sub=(x,s)=>`<msub>${mi(x)}<mtext>${s}</mtext></msub>`,sup=(x,s)=>`<msup>${x}${String(s).startsWith('<')?s:mo(s)}</msup>`,frac=(a,b)=>`<mfrac>${a}${b}</mfrac>`;
const alpha=sub('α','out'),tc=sub('T','cl'),tq=sub('T','q'),vin=sub('V','in'),vout=sub('V','out'),mq=sub('M','q');
const matrix=(values)=>`<mrow><mo>[</mo><mtable>${values.map(r=>'<mtr>'+r.map(c=>'<mtd>'+c+'</mtd>').join('')+'</mtr>').join('')}</mtable><mo>]</mo></mrow>`;
const representative={
 'i-state':row(sub('V','vac'),mo('='),frac(mi('I'),mn(2))),
 'i-drive':row(sub('I','0'),mo('='),frac(row(mn(2),mi('P')),row(mi('π'),sup(mi('w'),'2')))),
 'd-amplitude':row(alpha,mo('='),tc,sub('α','in')),
 'd-constant':row(tc,mo('='),sup(mi('e'),row(sub('M','cl'),mi('L')))),
 'd-maxwell':row(frac(row(mo('∂'),mi('α')),row(mo('∂'),mi('z'))),mo('='),sub('M','cl'),mi('α')),
 'd-gains':row(sub('G','p'),mo('='),frac(sub('P','p,out'),sub('P','0')),mo(','),sub('G','c'),mo('='),frac(sub('P','c,out'),sub('P','0'))),
 'm-k':row(mi('K'),mo('='),mo('−'),mq,mi('J'),mo('−'),mi('J'),sup(mq,'†')),
 'm-diffusion':row(mi('D'),mo('='),mi('B'),sub('N','res'),sup(mi('B'),'†')),
 'v-input':row(sub('V','input,out'),mo('='),tq,vin,sup(tq,'†')),
 'v-integral':row(mi('W'),mo('='),`<msubsup><mo>∫</mo><mn>0</mn><mi>L</mi></msubsup>`,tq,mo('('),mi('L'),mo(','),mi('z'),mo(')'),mi('D'),mo('('),mi('z'),mo(','),mi('Ω'),mo(')'),sup(tq,'†'),mo('('),mi('L'),mo(','),mi('z'),mo(')'),mi('d'),mi('z')),
 'v-output':row(vout,mo('='),tq,vin,sup(tq,'†'),mo('+'),mi('W')),
 'v-phi':row(sub('φ','1'),mo('('),mi('z'),mo(')'),mo('='),frac(row(sup(mi('e'),'z'),mo('−'),mn(1)),mi('z')),mo(','),sub('φ','1'),mo('('),mn(0),mo(')'),mo('='),mn(1)),
 'f-loss':row(sub('V','R,det'),mo('='),sub('K','η'),sub('V','R,out'),sup(sub('K','η'),'T'),mo('+'),sub('K','ℓ'),frac(sub('I','4'),mn(2)),sup(sub('K','ℓ'),'T')),
 'x-canonical':row(mi('T'),mi('J'),sup(mi('T'),'†'),mo('='),mi('J'))
};
function mathFor(id){const inner=representative[id];return inner?`<div class="math-display"><math xmlns="http://www.w3.org/1998/Math/MathML" display="block">${inner}</math></div>`:'';}
function selectNode(id,update=true){
 const n=nodeById.get(id);if(!n)return;
 if(state.group!==n.group){openGroup(n.group,false);}
 state.node=id;if(update)hash(id);
 const source=document.getElementById(id).querySelector('.node-content').cloneNode(true);
 source.querySelectorAll('[id]').forEach(x=>x.removeAttribute('id'));
 panel.innerHTML=`<p class="inspector-kicker">${numbers[n.group]} / ${esc(id)}</p><h2>${esc(n.title)}</h2>`;panel.append(source);
 const rendered=mathFor(id),equation=source.querySelector('.equation-box');
 if(rendered&&equation){const holder=document.createElement('div');holder.innerHTML=rendered;equation.before(holder.firstChild);const details=document.createElement('details');details.className='raw-equation';const summary=document.createElement('summary');summary.textContent='전체 식과 표기 읽기';details.append(summary);equation.before(details);details.append(equation);}
 const connections=source.querySelector('.node-connections');if(connections)connections.open=true;
 const local=state.links.map((e,i)=>({e,i})).filter(({e})=>e.source===id||e.target===id);
 if(local.length){const section=document.createElement('div');section.className='inspector-section';section.innerHTML='<h3>강조된 내부 연결 / 선 번호</h3><ol class="route-text">'+local.map(({e,i})=>`<li value="${i+1}"><span class="route-kind">${esc(edgeKinds[e.kind])}${e.pending?' · ◇ 미완성':''}</span>${esc(e.label)}</li>`).join('')+'</ol>';panel.append(section);}
 const b=document.createElement('button');b.className='view-record';b.dataset.readRecord=id;b.textContent='전체 문서에서 이 노드 읽기 ↗';panel.append(b);panel.scrollTop=0;highlight();drawMinimap();
 document.getElementById('announcement').textContent=n.title+' — '+n.out;
 if(matchMedia('(max-width:720px)').matches&&update){panel.scrollIntoView({block:'start',behavior:matchMedia('(prefers-reduced-motion:reduce)').matches?'instant':'smooth'});panel.focus({preventScroll:true});}
}
function centerNode(id){const p=state.positions.get(id);if(!p)return;state.x=viewport.clientWidth/2-(p.x+p.w/2)*state.scale;state.y=viewport.clientHeight/2-(p.y+p.h/2)*state.scale;paintCamera();}
function openGroup(id,update=true){
 if(id!=='overview'&&!groupById.has(id))return;
 resetReader(false);state.group=id;state.node=null;document.body.classList.toggle('atlas-overview',id==='overview');
 app.querySelectorAll('[data-open-group]').forEach(b=>b.setAttribute('aria-current',String(b.dataset.openGroup===id)));
 if(update)hash(id==='overview'?'atlas':'group-'+id);
 const g=groupById.get(id);
 document.getElementById('atlas-title').textContent=g?g.title:'FWM 이론의 구조와 성립 조건';
 document.getElementById('atlas-crumb').textContent=g?'/ '+numbers[id]+' · '+g.subtitle:'/ 계산의 두 경로';
 document.getElementById('atlas-subtitle').textContent=g?'노드를 선택하면 수식·조건과 직접 연결된 입력·출력을 함께 읽습니다.':'평균장과 공분산의 합류를 보존하고, 아래 세 그룹에서 논리적 관계와 조건부 미시 결과를 읽습니다.';
 document.getElementById('atlas-count').textContent=g?g.nodes.length+' NODES':data.groups.length+' GROUPS';
 layout(id);renderGraph();if(g)groupPanel(g);requestAnimationFrame(()=>fit(id==='overview'));
}
function readRecord(id){resetReader(true);const el=document.getElementById(id);if(el){for(let p=el;p;p=p.parentElement)if(p.tagName==='DETAILS')p.open=true;hash(id);requestAnimationFrame(()=>el.scrollIntoView({block:'start'}));}}
function handleHash(){let id;try{id=decodeURIComponent(location.hash.slice(1));}catch{return;}
 if(nodeById.has(id)){selectNode(id,false);requestAnimationFrame(()=>centerNode(id));}
 else if(id.startsWith('group-')&&groupById.has(id.slice(6)))openGroup(id.slice(6),false);
 else if(id==='atlas'||!id)openGroup('overview',false);
 else if(document.getElementById(id))readRecord(id);
}
app.addEventListener('click',e=>{
 const g=e.target.closest('[data-open-group]'),n=e.target.closest('[data-open-node]'),r=e.target.closest('[data-read-record]'),v=e.target.closest('[data-atlas-view]');
 if(g){openGroup(g.dataset.openGroup);return;}
 if(n){selectNode(n.dataset.openNode,true);requestAnimationFrame(()=>centerNode(n.dataset.openNode));return;}
 if(r){readRecord(r.dataset.readRecord);return;}
 if(v){state.view=v.dataset.atlasView;app.querySelectorAll('[data-atlas-view]').forEach(b=>b.setAttribute('aria-pressed',String(b===v)));renderGraph();return;}
});
// Capture internal links in both inspector and reading view before static
// atlas navigation; this avoids competing scroll and selection handlers.
document.addEventListener('click',e=>{
 const a=e.target.closest('a[href^="#"]');if(!a||!document.body.classList.contains('explorer-active'))return;
 const id=a.getAttribute('href').slice(1);
 if(nodeById.has(id)||id.startsWith('group-')||id==='atlas'){
  e.preventDefault();e.stopImmediatePropagation();hash(id);
  if(state.reader){readRecord(id);return;}
  if(nodeById.has(id)){selectNode(id,false);requestAnimationFrame(()=>centerNode(id));}else openGroup(id==='atlas'?'overview':id.slice(6),false);
 }else if(document.getElementById(id)){e.preventDefault();e.stopImmediatePropagation();readRecord(id);}
},true);
document.getElementById('atlas-reader').addEventListener('click',()=>{resetReader(!state.reader);window.scrollTo({top:0});});
document.getElementById('atlas-print').addEventListener('click',()=>window.print());
document.getElementById('atlas-download').addEventListener('click',()=>{
 const clone=document.documentElement.cloneNode(true);clone.querySelector('#explorer')?.remove();clone.classList.remove('enhanced');
 const body=clone.querySelector('body');body.classList.remove('explorer-active','atlas-overview','reader-mode','all-expanded','mode-math','mode-evidence');
 clone.querySelectorAll('details').forEach(d=>d.setAttribute('open',''));clone.querySelectorAll('.svg-wires').forEach(s=>s.replaceChildren());
 clone.querySelector('#inspector')?.replaceChildren();clone.querySelectorAll('[aria-current]').forEach(el=>el.removeAttribute('aria-current'));
 // Source and code links are preserved locally. Hosted output rewrites these
 // to the embedded provenance catalog so its download is self-contained too.
 const blob=new Blob(['<!doctype html>\n'+clone.outerHTML],{type:'text/html;charset=utf-8'}),url=URL.createObjectURL(blob),a=document.createElement('a');
 a.href=url;a.download='fwm_quotient_structure_v4.html';document.body.append(a);a.click();a.remove();setTimeout(()=>URL.revokeObjectURL(url),1000);
 document.getElementById('announcement').textContent=nodeById.size+'개 노드와 모든 수식·조건·연결이 포함된 HTML을 저장합니다.';
});
document.getElementById('zoom-fit').addEventListener('click',()=>fit(true));document.getElementById('zoom-in').addEventListener('click',()=>zoom(1.2));document.getElementById('zoom-out').addEventListener('click',()=>zoom(1/1.2));
viewport.addEventListener('wheel',e=>{e.preventDefault();const r=viewport.getBoundingClientRect();zoom(Math.exp(-Math.max(-120,Math.min(120,e.deltaY))*.003),e.clientX-r.left,e.clientY-r.top);},{passive:false});
const pointers=new Map();let gesture=null,moved=false;
function gestureBase(){const p=[...pointers.values()];if(p.length>=2){gesture={type:'pinch',cx:(p[0].x+p[1].x)/2,cy:(p[0].y+p[1].y)/2,d:Math.hypot(p[0].x-p[1].x,p[0].y-p[1].y),scale:state.scale,x:state.x,y:state.y};}else if(p.length===1){gesture={type:'pan',cx:p[0].x,cy:p[0].y,x:state.x,y:state.y};}else gesture=null;}
viewport.addEventListener('pointerdown',e=>{if(e.button!==0)return;const r=viewport.getBoundingClientRect();pointers.set(e.pointerId,{x:e.clientX-r.left,y:e.clientY-r.top});gestureBase();moved=false;if(!e.target.closest('button'))viewport.setPointerCapture(e.pointerId);});
viewport.addEventListener('pointermove',e=>{if(!pointers.has(e.pointerId)||!gesture)return;const r=viewport.getBoundingClientRect();pointers.set(e.pointerId,{x:e.clientX-r.left,y:e.clientY-r.top});const p=[...pointers.values()];if(gesture.type==='pinch'&&p.length>=2){const cx=(p[0].x+p[1].x)/2,cy=(p[0].y+p[1].y)/2,d=Math.hypot(p[0].x-p[1].x,p[0].y-p[1].y);state.scale=Math.min(2.25,Math.max(.15,gesture.scale*d/Math.max(1,gesture.d)));state.x=cx-(gesture.cx-gesture.x)*state.scale/gesture.scale;state.y=cy-(gesture.cy-gesture.y)*state.scale/gesture.scale;moved=true;}else{const dx=p[0].x-gesture.cx,dy=p[0].y-gesture.cy;if(Math.hypot(dx,dy)>5)moved=true;if(moved){state.x=gesture.x+dx;state.y=gesture.y+dy;}}
 if(moved){viewport.classList.add('dragging');paintCamera();}});
const endPointer=e=>{pointers.delete(e.pointerId);if(viewport.hasPointerCapture(e.pointerId))viewport.releasePointerCapture(e.pointerId);gestureBase();if(!pointers.size)viewport.classList.remove('dragging');};
viewport.addEventListener('pointerup',endPointer);viewport.addEventListener('pointercancel',endPointer);viewport.addEventListener('lostpointercapture',e=>{pointers.delete(e.pointerId);gestureBase();});
viewport.addEventListener('click',e=>{if(moved){e.preventDefault();e.stopImmediatePropagation();moved=false;}},true);
viewport.addEventListener('keydown',e=>{if(e.target!==viewport)return;if(['ArrowLeft','ArrowRight','ArrowUp','ArrowDown','+','=','-','0','Escape'].includes(e.key)){e.preventDefault();if(e.key==='0')fit(true);else if(e.key==='+'||e.key==='=')zoom(1.2);else if(e.key==='-')zoom(1/1.2);else if(e.key==='Escape'){state.node=null;if(groupById.has(state.group))groupPanel(groupById.get(state.group));highlight();}else{state.x+=e.key==='ArrowLeft'?70:e.key==='ArrowRight'?-70:0;state.y+=e.key==='ArrowUp'?70:e.key==='ArrowDown'?-70:0;paintCamera();}}});
window.addEventListener('hashchange',handleHash);window.addEventListener('popstate',handleHash);
new ResizeObserver(()=>{if(!state.reader)paintCamera();}).observe(viewport);
openGroup('overview',false);handleHash();
})();
