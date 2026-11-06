// Scientific views use the same SVG, nearest-sample inspector and pointer lifecycle as OD/SAS.
const spectralSessions = new Map();
const conditionColors = ['#07655e','#a54c27','#565095','#426b92','#6b7240','#9b3e65','#70604a','#334f68'];
const conditionDashes = ['', '7 4', '2 4', '8 3 2 3', '4 2', '10 3', '1 3', '6 2 1 2'];

export function conditionNavigationGrid(dataset, id, fallback) {
  const native=dataset.series.find(row=>row.id===id)?.x.filter(Number.isFinite);
  return native?.length ? native : fallback;
}

function plottedCurves(block, panel) {
  return panel.curves.map((definition, index) => {
    const condition = block.series.find(c => c.id === definition.series_id);
    const colorIndex=block.series.indexOf(condition);
    return {x:condition.observable_x?.[definition.observable] || condition.x, y:condition.observables[definition.observable],
      sample_status:condition.sample_status,sample_reasons:condition.sample_reasons,
      observable_status:condition.observable_status?.[definition.observable],
      unavailable_reason:condition.unavailable_observables?.[definition.observable],
      shortLabel:condition.parameter ? String(condition.parameter.value)+(condition.parameter.unit || '') : condition.label,
      label:definition.label || condition.label, color:conditionColors[colorIndex % conditionColors.length],
      cyclic_period:panel.cyclic_period,points_only:panel.sampled_only,
      condition_id:condition.id,condition_ids:condition.condition_ids,
      dash:definition.dash ?? conditionDashes[colorIndex % conditionDashes.length]};
  });
}

export function panelCurves(block,panel,enabled=[],selectedId=null) {
  const filtered={...panel,curves:panel.curves.filter(c=>(c.optional_overlays || []).every(key=>enabled.includes(key)))};
  const curves=plottedCurves(block,filtered),row=block.series.find(row=>row.id===selectedId);
  for(const key of enabled){const overlay=row?.overlays?.[key];if(overlay?.panel===panel.key){
    const base=curves.find(c=>c.condition_id===selectedId);
    curves.push({x:overlay.x,y:overlay.y,label:overlay.label,color:base?.color || conditionColors[0],dash:'5 4',condition_id:selectedId});
  }}
  return curves;
}

function spectralRange(curves, window, reference, logarithmic=false) {
  let low=Infinity, high=-Infinity;
  for(const curve of curves) for(let i=0;i<curve.x.length;i++) {
    if(curve.x[i]<window[0] || curve.x[i]>window[1]) continue;
    const y=curve.y[i];
    if(Number.isFinite(y) && (!logarithmic || y>0)) {low=Math.min(low,y);high=Math.max(high,y);}
  }
  if(!Number.isFinite(low)) return logarithmic ? [1,10] : [0,1];
  if(logarithmic){const a=Math.log10(low),b=Math.log10(high),pad=(b-a || .2)*.08;return [10**(a-pad),10**(b+pad)];}
  if(Number.isFinite(reference)){low=Math.min(low,reference);high=Math.max(high,reference);}
  const pad=(high-low || Math.abs(high)*.1 || 1)*.08;
  return [low-pad,high+pad];
}

export function stackedCurves(curves, window) {
  const [lo,hi]=spectralRange(curves,window),span=hi-lo || 1,step=span*1.12;
  const axisTicks=ticks(lo,hi,3),scaleStep=axisTicks.length>1?axisTicks[1]-axisTicks[0]:span/3;
  return {series:curves.map((curve,index)=>({...curve,display_offset:index*step})),
    ylim:[lo,hi+Math.max(0,curves.length-1)*step],scaleStep};
}

function readable(value, unit='', status='') {
  const text=typeof value==='number' ? (Number.isFinite(value) && value!==0 && (Math.abs(value)>=1e6 || Math.abs(value)<1e-4) ? value.toExponential(5).replaceAll('-','−') : format(value,6))
    : value==null ? (status || 'unavailable') : String(value);
  return text+(value!=null && unit ? ' '+unit : '')+(status && value!=null && status!=='valid' && status!=='resolved' ? ' · '+status : '');
}

export default function mountSpectral({parentElement,data,key,setStateValue}) {
  const root=parentElement.querySelector('.sc');
  if(!root) return;
  const d=data.dataset, grid=data.grid, prior=spectralSessions.get(key);
  if(prior?.root===root && prior.context===data.context && prior.revision===data.revision) return;
  prior?.cleanup();
  const selection={...(prior?.context===data.context ? prior.selection : data.selection)};
  selection.window=[...selection.window];
  if(!d.series.some(row=>row.id===selection.series_id))selection.series_id=data.selection.series_id || d.series[0].id;
  if(!d.panels.some(panel=>panel.key===selection.panel_key))selection.panel_key=d.panels[0].key;
  selection.display_mode=selection.display_mode || 'curves';
  const domain=[grid[0],grid[grid.length-1]], pointer={drag:null};
  let separate=prior?.separate ?? true;
  root.innerHTML=`
    <div class="sc-heading"><h3></h3><p class="mc-frequency-reference mc-muted"></p></div>
    <div class="mc-legend sc-legend" aria-label="Conditions compared"></div>
    <div class="sc-profiles sc-view-switch" role="group" aria-label="Analysis"></div>
    <p class="sc-profile-note mc-muted" hidden></p>
    <div class="sc-explore-controls">
      <div class="sc-view-switch sc-display-switch" role="group" aria-label="Comparison display"><button type="button" data-view="curves">Curves</button><button type="button" data-view="map">Map</button></div>
      <label>Quantity <select class="sc-panel-choice"></select></label>
      <label>Selected condition <select class="sc-condition-choice"></select></label>
      <label class="sc-axis-label" hidden>Map vertical axis <select class="sc-axis-choice"></select></label>
    </div>
    <div class="sc-options"></div><button type="button" class="sc-optimum" hidden>Inspect sampled optimum</button>
    <label class="sc-separate" hidden><input type="checkbox"> Separate transmission curves</label>
    <div class="mc-layout"><div class="mc-graphs">
      <div class="mc-overview"></div>
      <div class="sc-map" hidden></div>
      <details class="mc-range"><summary>Zoom window</summary><p>Drag the shaded region or its edges. Drag the dotted cursor inside the region. Click outside to recenter.</p><div class="sc-range">
        <label>From <input type="number" step="any"></label>
        <label>To <input type="number" step="any"></label>
        <button type="button" class="sc-full">Full scan</button></div></details>
      <p class="sc-scope-note mc-muted" hidden>Shaded intervals include samples outside model checks or unavailable samples. Inspect the selected sample's scope below.</p>
      <div class="sc-panels"></div><div class="sc-summary"></div>
    </div><aside class="mc-reader" aria-label="Selected samples and derived readouts">
      <h3>Selected sample</h3><p class="mc-loc"></p>
      <div class="mc-cursor"><label><span class="sc-xlabel"></span><input type="number" step="any"></label>
        <button type="button" class="sc-prev" aria-label="Previous sample">←</button>
        <button type="button" class="sc-next" aria-label="Next sample">→</button></div>
      <div class="sc-scroll"><table class="mc-values sc-samples"><caption>Nearest native sample for each condition</caption>
        <thead><tr><th>Condition</th><th class="sc-sample-x"></th><th>Observable</th><th>Value</th></tr></thead><tbody></tbody></table></div>
      <label class="sc-all-values"><input type="checkbox"> Show all observables</label>
      <div class="sc-sample-reasons mc-muted"></div>
      <p class="mc-muted sc-native-note">Each curve retains its own solved grid. Rows report that condition's actual sampled coordinate.</p>
      <h3 class="sc-readout-title">Condition readouts</h3><p class="mc-muted sc-readout-description">Derived from each complete scan; independent of the inspection cursor.</p>
      <div class="sc-readouts"></div><div class="sc-notes"></div>
      <details class="mc-model"><summary>Calculation conditions and scope</summary><p class="sc-fixed"></p><p class="sc-source"></p><div class="sc-all-notes"></div><div class="sc-scopes"></div></details>
      <button type="button" class="mc-export">Download comparison JSON</button>
      <output class="mc-sr" aria-live="polite"></output>
    </aside></div>`;
  const q=s=>root.querySelector(s);
  function profile(){return d.profiles?.find(p=>p.key===selection.profile) || d.profiles?.[0];}
  function visiblePanels(){const keys=profile()?.panels;return keys?d.panels.filter(p=>keys.includes(p.key)):d.panels;}
  for(const p of d.profiles || []){const button=element('button',p.label);button.type='button';button.dataset.profile=p.key;button.onclick=()=>{selection.profile=p.key;selection.panel_key=p.default_panel;refresh();publish();};q('.sc-profiles').append(button);}
  for(const option of d.optional_overlays || []){const label=element('label'),check=element('input');check.type='checkbox';check.dataset.overlay=option.key;
    check.onchange=()=>{selection.overlays=[...q('.sc-options').querySelectorAll('input:checked')].map(input=>input.dataset.overlay);refresh();publish();};label.append(check,document.createTextNode(option.label));q('.sc-options').append(label);}
  q('.sc-optimum').onclick=()=>{const point=d.series.find(row=>row.id===selection.series_id)?.optimum;if(Number.isFinite(point?.x)){choose(point.x);refresh();publish();}};
  const panelChoice=q('.sc-panel-choice'),conditionChoice=q('.sc-condition-choice');
  const axisChoice=q('.sc-axis-choice');
  for(const axis of d.condition_axes || []){const option=element('option',axis.label+' ['+axis.unit+']');option.value=axis.key;axisChoice.append(option);}
  axisChoice.onchange=()=>{selection.condition_axis=axisChoice.value;refresh();publish();};
  function quantityOptions(){const panels=visiblePanels(),signature=panels.map(p=>p.key).join('|');
    if(panelChoice.dataset.panels!==signature){panelChoice.replaceChildren();for(const panel of panels){const option=element('option',panel.y_label+(panel.unit && panel.unit!=='1'?' ['+panel.unit+']':''));option.value=panel.key;panelChoice.append(option);}panelChoice.dataset.panels=signature;}
    if(!panels.some(panel=>panel.key===selection.panel_key))selection.panel_key=profile()?.default_panel || panels[0].key;
  }
  for(const row of d.series){const option=element('option',row.label);option.value=row.id;conditionChoice.append(option);}
  panelChoice.onchange=()=>{selection.panel_key=panelChoice.value;refresh();publish();};
  conditionChoice.onchange=()=>{selectCondition(conditionChoice.value,selection.x);refresh();publish();};
  for(const button of q('.sc-display-switch').querySelectorAll('button'))button.onclick=()=>{selection.display_mode=button.dataset.view;refresh();publish();};
  q('.sc-all-values input').onchange=()=>refresh();
  q('.sc-heading h3').textContent=d.title || (d.mode ? d.mode+' · condition comparison' : 'Condition comparison');
  q('.mc-frequency-reference').textContent=d.axis.label+' ['+d.axis.unit+'] · '+d.axis.reference+(d.axis.scale==='log'?' · logarithmic axis':'');
  q('.sc-xlabel').textContent=d.axis.label+' ('+d.axis.unit+')';
  q('.sc-sample-x').textContent='Sample ('+d.axis.unit+')';
  const from=q('.sc-range label:first-child input'),to=q('.sc-range label:nth-child(2) input');
  from.setAttribute('aria-label','From ('+d.axis.unit+')');to.setAttribute('aria-label','To ('+d.axis.unit+')');
  const freq=q('.mc-cursor input');
  const firstCurves=plottedCurves(d,d.panels[0]);
  q('.sc-separate').hidden=!d.panels.some(panel=>panel.stackable);
  q('.sc-separate input').checked=separate;
  q('.sc-separate input').onchange=e=>{separate=e.target.checked;entry.separate=separate;refresh();};
  q('.sc-scope-note').hidden=!d.series.some(c=>c.sample_status?.some(status=>status!=='available'));
  const legend=q('.sc-legend');
  function conditionLegend(){
    legend.replaceChildren();
    const compact=selection.display_mode==='map' || firstCurves.length>8;
    const visible=compact?firstCurves.filter(curve=>curve.condition_id===selection.series_id):firstCurves;
    if(compact)legend.append(element('span',d.series.length+' calculated conditions'));
    visible.forEach(c=>{
    const item=element('span'),line=element('i');
    line.style.borderTopColor=c.color;line.style.borderTopStyle=c.dash?'dashed':'solid';
    item.append(line,document.createTextNode(c.label));legend.append(item);
    });
  }
  const notes=q('.sc-notes');
  (d.notes || []).forEach(note=>q('.sc-all-notes').append(element('p',note,'mc-model')));
  q('.sc-fixed').textContent=Object.entries(d.fixed || {}).map(([k,v])=>k+' = '+String(v)).join(' · ');
  q('.sc-source').textContent=Object.entries(d.source || {}).filter(([,v])=>typeof v!=='object').map(([k,v])=>k+': '+String(v)).join(' · ');
  for(const c of d.series)if(c.scope){
    const details=element('details'),summary=element('summary',c.label+' · model and validity');
    details.append(summary,element('pre',JSON.stringify(c.scope,null,2)));q('.sc-scopes').append(details);
  }
  for(const marker of d.markers||[])q('.sc-notes').append(element('p',marker.label+' · '+readable(marker.x,d.axis.unit),'mc-model'));
  function conditionReadouts(){
    const readouts=q('.sc-readouts');readouts.replaceChildren();
    const row=d.series.find(row=>row.id===selection.series_id);
    q('.sc-readout-description').textContent=profile()?.readout_description || 'Derived from each complete scan; independent of the inspection cursor.';
    const readingKeys=profile()?.readouts;
    const readings=readingKeys?readingKeys.map(key=>row?.readouts?.find(item=>item.key===key)).filter(Boolean):(row?.readouts || []);
    for(const item of readings){
      const section=element('div',undefined,'sc-metric');
      section.append(element('h4',item.label),element('p',readable(item.value,item.unit==='1'?'':item.unit,item.status)));
      if(item.reason)section.append(element('small',item.reason));readouts.append(section);
    }
    q('.sc-readout-title').hidden=!row?.readouts?.length;
  }
  let frame=0, width=0;
  const listeners=[];
  function publish() {
    selection.condition=d.series.find(row=>row.id===selection.series_id)?.parameter || null;
    const value={...selection,window:[...selection.window],context:data.context};
    const stamp=JSON.stringify(value);
    if(stamp===entry.lastSent)return;
    entry.lastSent=stamp;setStateValue('selection',value);
  }
  function choose(frequency) {
    const [a,b]=selection.window;
    const native=conditionNavigationGrid(d,selection.series_id,grid);
    let i=nearestSample(native,Math.min(b,Math.max(a,frequency)));
    if(native[i]<a && i+1<native.length)i++;
    if(native[i]>b && i && native[i-1]>=a)i--;
    if(native[i]<a || native[i]>b)selection.window=[Math.min(a,native[i]),Math.max(b,native[i])];
    selection.x=native[i];
  }
  function selectCondition(id,frequency){
    const row=d.series.find(row=>row.id===id);
    if(!row)return;
    selection.series_id=id;
    const sample=curveSample({x:row.x,y:row.observables[d.panels[0].curves.find(c=>c.series_id===id)?.observable] || []},grid,frequency);
    choose(sample.x==null?frequency:sample.x);
  }
  function apply(result) {
    selection.window=[...result.window];choose(result.cursor);
  }
  function shape(block,panel,host,window,overview=false,summary=false) {
    let curves=panelCurves(block,panel,selection.overlays || [],summary?null:selection.series_id);
    if(!overview && !summary && selection.display_mode==='map')curves=curves.filter(curve=>curve.condition_id===selection.series_id);
    if(summary){
      const selected=[];
      for(const curve of curves){const i=curve.condition_ids?.indexOf(selection.series_id);if(i>=0 && Number.isFinite(curve.y[i]))selected.push({...curve,x:[curve.x[i]],y:[curve.y[i]],label:'Selected condition · '+curve.label,markers:true,point_kind:'selected',legend:false});}
      const objective=curves.find(c=>c.condition_id===panel.optimum_series_id);
      if(objective && panel.goal){const valid=objective.y.map((y,i)=>({x:objective.x[i],y})).filter(p=>Number.isFinite(p.y));
        if(valid.length){const best=valid.reduce((a,b)=>(panel.goal==='min'?b.y<a.y:b.y>a.y)?b:a);selected.push({...objective,x:[best.x],y:[best.y],label:'Best sampled condition · '+format(best.x,6)+' '+block.axis.unit,markers:true,point_kind:'optimum'});}}
      curves=[...curves,...selected];
    }
    const stacking=!overview && panel.stackable && separate ? stackedCurves(curves,window) : null;
    if(stacking)curves=stacking.series;
    const referenceGrid=summary ? [...new Set(block.series.flatMap(s=>s.x))].sort((a,b)=>a-b) : grid;
    const spec={title:panel.title,grid:referenceGrid,index:nearestSample(referenceGrid,summary?(d.series.find(row=>row.id===selection.series_id)?.parameter?.value ?? selection.x):selection.x),
      unit:block.axis.unit,divisor:1,xLabel:block.axis.label,xlim:window,
      xScale:block.axis.scale,tickLabels:block.axis.tick_labels,discrete:Boolean(block.axis.tick_labels),
      marks:summary?[]:[...(d.markers||[]),...(d.series.find(row=>row.id===selection.series_id)?.markers || [])].filter(m=>!m.profiles || m.profiles.includes(profile()?.key)).map(m=>({...m,label:m.display_label || m.label})),cursor:!summary,
      ylim:stacking?.ylim || spectralRange(curves,window,panel.reference,panel.scale==='log'),series:curves,
      yScale:panel.scale,yUnit:panel.unit,stacked:stacking?{scaleStep:stacking.scaleStep,scaleUnit:panel.stack_scale_unit,label:panel.stack_label}:null,
      height:stacking?Math.max(250,curves.length*34+65):undefined,
      yLabel:stacking?(panel.stack_label || 'Condition')+' · '+panel.y_label+' (offset)':panel.y_label+(panel.unit && panel.unit!=='1'?' ['+panel.unit+']':'')+(panel.scale==='log'?' · log':''),
      reference:panel.reference,referenceLabel:panel.reference_label,
      keyboardInspector:!overview,
      description:panel.title+'. Nearest solved samples are shown for each condition.'};
    if(overview) spec.band=selection.window;
    spectrum(host,spec);
    if(!overview && (panel.note || stacking))host.append(element('p',stacking?'Vertical offsets only; all curves share the same transmission scale. Hover and the readout table show the original values.':panel.note,'sc-panel-note mc-muted'));
    if(!overview && (summary || curves.some((c,i)=>c.label!==firstCurves[i]?.label))){
      const legend=element('div',undefined,'mc-legend sc-panel-legend');
      curves.filter(c=>c.legend!==false).forEach(c=>{const item=element('span'),line=element('i');line.style.borderTopColor=c.color;line.style.borderTopStyle=c.dash?'dashed':'solid';item.append(line,document.createTextNode(c.label));legend.append(item);});
      host.append(legend);
    }
    if(summary)host.querySelector('svg')._navigation.summary=true;
  }
  function drawCuts(){
    const panels=q('.sc-panels');panels.replaceChildren();
    // Mount every grid child before measuring: :only-child otherwise makes the
    // first SVG measure a full row and then shrink after the second is added.
    const available=visiblePanels(),active=available.find(panel=>panel.key===selection.panel_key),groups=new Map();
    for(const panel of available)if(panel.cut_group && !groups.has(panel.cut_group))groups.set(panel.cut_group,panel.key);
    if(active?.cut_group)groups.set(active.cut_group,active.key);
    const visible=available.filter(panel=>!panel.cut_group || groups.get(panel.cut_group)===panel.key);
    const ordered=visible[0].wide ? visible : [...visible.filter(panel=>!panel.wide),...visible.filter(panel=>panel.wide)];
    const narrow=ordered.filter(panel=>!panel.wide),lastNarrow=narrow.length%2?narrow.at(-1):null;
    const hosts=ordered.map(panel=>{const host=element('div');if(panel.wide || panel===lastNarrow)host.className='sc-panel-wide';return host;});panels.append(...hosts);
    ordered.forEach((panel,i)=>shape(d,panel,hosts[i],selection.window));
  }
  function drawSummary(){
    const summary=q('.sc-summary');summary.replaceChildren();
    const panels=(d.summary?.panels || []).filter(p=>!profile()?.summaries || profile().summaries.includes(p.key));
    summary.classList.toggle('sc-summary-grid',panels.length>1);
    const hosts=panels.map(()=>element('div'));summary.append(...hosts);
    if(d.summary)for(const [i,panel] of panels.entries()){
      const host=hosts[i];
      const xs=d.summary.series.flatMap(c=>c.x).filter(Number.isFinite), a=Math.min(...xs),b=Math.max(...xs);
      const pad=a===b?(Math.abs(a)*.05 || 1):0;
      shape(d.summary,panel,host,[a-pad,b+pad],false,true);
    }
  }
  function draw() {
    for(const stop of listeners.splice(0))stop();
    const overview=q('.mc-overview');overview.replaceChildren();
    shape(d,{...d.panels[0],title:'Full range · '+d.panels[0].title},overview,domain,true);
    const host=q('.sc-map');host.replaceChildren();host.hidden=selection.display_mode!=='map';
    const panel=d.panels.find(panel=>panel.key===selection.panel_key);
    if(!host.hidden)spectrumMap(host,{dataset:d,panel,window:selection.window,selection,onKeyboardSelect:(id,x)=>{entry.restoreFocus=true;selectCondition(id,x);refresh();publish();}});
    drawCuts();drawSummary();
    if(entry.restoreFocus)host.querySelector('svg')?.focus({preventScroll:true});
  }
  function refresh(preview=false) {
    quantityOptions();
    for(const button of q('.sc-profiles').querySelectorAll('button'))button.setAttribute('aria-pressed',String(button.dataset.profile===profile()?.key));
    q('.sc-profiles').hidden=!d.profiles?.length;
    const note=q('.sc-profile-note');note.textContent=profile()?.note || '';note.hidden=!note.textContent;
    for(const option of d.optional_overlays || []){const input=q('input[data-overlay="'+option.key+'"]');input.checked=(selection.overlays || []).includes(option.key);input.parentElement.hidden=Boolean(option.profiles && !option.profiles.includes(profile()?.key));}
    q('.sc-optimum').hidden=!profile()?.inspect_optimum || !Number.isFinite(d.series.find(row=>row.id===selection.series_id)?.optimum?.x);
    conditionLegend();
    panelChoice.value=selection.panel_key;conditionChoice.value=selection.series_id;
    axisChoice.value=selection.condition_axis || 'parameter';
    q('.sc-axis-label').hidden=selection.display_mode!=='map' || !d.condition_axes?.length;
    for(const button of q('.sc-display-switch').querySelectorAll('button')){button.setAttribute('aria-pressed',String(button.dataset.view===selection.display_mode));if(button.dataset.view==='map')button.disabled=d.series.length<2;}
    freq.value=String(selection.x);freq.min=domain[0];freq.max=domain[1];
    from.value=selection.window[0];to.value=selection.window[1];
    q('.mc-loc').textContent=d.axis.tick_labels?.[Math.round(selection.x)] || format(selection.x,8)+' '+d.axis.unit;
    const body=q('.sc-samples tbody');body.replaceChildren();
    const reasons=new Set();
    for(const panel of visiblePanels().filter(panel=>q('.sc-all-values input').checked || panel.key===selection.panel_key))for(const curve of panelCurves(d,panel,selection.overlays || [],selection.series_id)){
      const row=element('tr'),point=curveSample(curve,grid,selection.x);
      for(const text of [curve.label,point.x==null?(curve.unavailable_reason?'—':'outside scan'):d.axis.tick_labels?.[Math.round(point.x)] || format(point.x,8),panel.y_label,readable(point.y,panel.unit==='1'?'':panel.unit,point.y==null?curve.unavailable_reason:'')])
        row.append(element('td',text));
      const scope=sampleScope(curve,grid,point.x);
      if(scope.status && scope.status!=='available' && !curve.unavailable_reason)row.lastElementChild.append(element('small',scope.status,'sc-sample-status'));
      scope.reasons.forEach(reason=>reasons.add(reason));
      if(curve.condition_id===selection.series_id)row.className='sc-selected-row';body.append(row);
    }
    const scopeNotes=q('.sc-sample-reasons');scopeNotes.replaceChildren();
    for(const reason of reasons)scopeNotes.append(element('p',reason));
    q('output').textContent='Selected '+format(selection.x,8)+' '+d.axis.unit+'. '+body.textContent;
    conditionReadouts();
    if(preview){
      q('.sc-map svg')?._navigation.mark(selection);
      if(selection.display_mode==='map')drawCuts();drawSummary();
    }else draw();
  }
  const minWidth=Math.max(Number.EPSILON,Math.min(domain[1]-domain[0],grid.length>1?grid[1]-grid[0]:1));
  function editWindow(){
    const a=from.valueAsNumber,b=to.valueAsNumber;
    if(!Number.isFinite(a)||!Number.isFinite(b)||a<domain[0]||b>domain[1]||b-a<minWidth){
      to.setCustomValidity('Choose an ordered window inside the calculated range.');return;
    }
    to.setCustomValidity('');apply({window:[a,b],cursor:selection.x});refresh();publish();
  }
  from.onchange=editWindow;to.onchange=editWindow;
  for(const input of [from,to])input.onkeydown=e=>{if(e.key==='Enter'){e.preventDefault();editWindow();}if(e.key==='Escape'){input.value=input===from?selection.window[0]:selection.window[1];}};
  freq.onchange=()=>{if(Number.isFinite(freq.valueAsNumber)){choose(freq.valueAsNumber);refresh();publish();}};
  function stepSample(delta){const native=conditionNavigationGrid(d,selection.series_id,grid);choose(native[Math.max(0,Math.min(native.length-1,nearestSample(native,selection.x)+delta))]);refresh();publish();}
  q('.sc-prev').onclick=()=>stepSample(-1);
  q('.sc-next').onclick=()=>stepSample(1);
  q('.sc-full').onclick=()=>{apply({window:domain,cursor:selection.x});refresh();publish();};
  q('.mc-export').onclick=()=>{
    const blob=new Blob([JSON.stringify({...d,selection:{...selection,window:[...selection.window]}},null,2)],{type:'application/json'});
    const url=URL.createObjectURL(blob),a=document.createElement('a');a.href=url;a.download=d.scheme+'-comparison.json';a.click();setTimeout(()=>URL.revokeObjectURL(url),1000);
  };
  const graphs=q('.mc-graphs');
  const stop=bindSpectrumGesture(graphs,pointer,{
    begin(e){
      const svg=e.composedPath().find(n=>n._navigation && !n._navigation.summary);
      if(!svg)return null;
      const n=svg._navigation,{x,y}=spectrumPointer(svg,e,n.w,n.h);
      if(x<n.left||x>n.right||y<n.top||y>n.bottom)return null;
      if(n.kind==='condition-map')return {svg,action:'condition',window:[...selection.window],cursor:selection.x,series_id:selection.series_id};
      const frequency=n.toFrequency(x),[a,b]=selection.window;
      let action='cursor';
      if(n.overview){
        const tol=Math.min(7,Math.max(2,(n.sx(b)-n.sx(a))/5));
        if(Math.abs(x-n.sx(selection.x))<=tol)action='cursor';
        else if(Math.abs(x-n.sx(a))<=tol)action='left';
        else if(Math.abs(x-n.sx(b))<=tol)action='right';
        else action=frequency>=a&&frequency<=b?'pan':'center';
      }
      return {svg,action,window:[a,b],cursor:selection.x,x:frequency,startX:e.clientX};
    },
    preview(drag,e){
      const n=drag.svg._navigation,{x}=spectrumPointer(drag.svg,e,n.w,n.h);
      if(drag.action==='condition'){
        const p=spectrumPointer(drag.svg,e,n.w,n.h),sample=conditionMapSample(n.rows,n.toFrequency(p.x),n.toCondition(p.y),Math.abs(n.toFrequency(p.x+5)-n.toFrequency(p.x-5))/2);
        if(sample?.row.id && sample.point.x!=null)selectCondition(sample.row.id,sample.point.x);
        refresh(true);root.dataset.gesture='preview';return;
      }
      apply(navigateSpectrumWindow(drag.action,drag,n.toFrequency(x),domain,minWidth));
      const index=nearestSample(grid,selection.x);
      for(const svg of graphs.querySelectorAll('svg'))if(svg._navigation && !svg._navigation.summary && svg._navigation.kind!=='condition-map')previewMarkers(svg,selection.window,index);
      q('.sc-map svg')?._navigation.mark(selection);
      root.dataset.gesture='preview';q('.mc-loc').textContent='Preview · '+format(selection.x,8)+' '+d.axis.unit;
    },
    commit(drag,e){
      if(drag.action==='condition')entry.restoreFocus=true;
      if(drag.action==='pan'&&Math.abs(e.clientX-drag.startX)<2)choose(drag.x);
      root.dataset.gesture='ready';refresh();publish();
    },
    cancel(drag){if(drag.series_id)selection.series_id=drag.series_id;apply({window:drag.window,cursor:drag.cursor});root.dataset.gesture='ready';refresh();}
  });
  const observer=new ResizeObserver(records=>{const next=records[0].contentRect.width;if(Math.abs(next-width)<.5)return;width=next;cancelAnimationFrame(frame);frame=requestAnimationFrame(()=>refresh());});
  const entry={root,context:data.context,revision:data.revision,selection,separate,restoreFocus:prior?.restoreFocus || false,lastSent:JSON.stringify({...selection,context:data.context}),cleanup:()=>{observer.disconnect();stop();cancelAnimationFrame(frame);}};
  spectralSessions.set(key,entry);refresh();observer.observe(root);
  return ()=>{entry.cleanup();if(spectralSessions.get(key)===entry)spectralSessions.set(key,{context:data.context,separate:entry.separate,restoreFocus:entry.restoreFocus,selection:{...selection,window:[...selection.window]},cleanup:()=>{}});};
}
