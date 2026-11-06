// Scientific map cells own their native x coordinates. No interpolation, shared
// array index, display offset or normalization enters the sampled values.
export function conditionMapRows(dataset, panel, conditionAxis='parameter') {
  const axis=dataset.condition_axes?.find(axis=>axis.key===conditionAxis);
  const coordinate=value=>axis?.transform ? axis.transform.factor*Math.pow(value,axis.transform.exponent) : value;
  const definitions=new Map(panel.curves.map(c=>[c.series_id,c]));
  const rows=dataset.series.filter(c=>definitions.has(c.id)).map(condition=>{
    const definition=definitions.get(condition.id),parameter=condition.parameter || {};
    return {id:condition.id,label:condition.label,value:parameter.value,displayValue:condition.coordinates?.[conditionAxis] ?? coordinate(parameter.value),condition,
      curve:{x:condition.observable_x?.[definition.observable] || condition.x,y:condition.observables[definition.observable],points_only:panel.sampled_only,
        unavailable_reason:condition.unavailable_observables?.[definition.observable],
        sample_status:condition.sample_status,sample_reasons:condition.sample_reasons,
        observable_status:condition.observable_status?.[definition.observable]}};
  });
  for(const value of dataset.sweep?.requested_values || [])if(!rows.some(row=>row.value===value))
    rows.push({id:null,value,displayValue:coordinate(value),label:String(value)+(dataset.varied?.unit?' '+dataset.varied.unit:''),condition:null,curve:{x:[],y:[]}});
  const numeric=rows.length>1 && rows.every(row=>typeof row.displayValue==='number' && Number.isFinite(row.displayValue)) && new Set(rows.map(row=>row.displayValue)).size===rows.length;
  if(numeric)rows.sort((a,b)=>a.displayValue-b.displayValue);
  rows.forEach((row,i)=>{row.coordinate=numeric?row.displayValue:i;});
  return {rows,numeric,axis:axis || dataset.varied};
}

export function mapNativeCells(curve) {
  return curve.x.map((x,i)=>({x,y:curve.y[i],
    low:i?(curve.x[i-1]+x)/2:x,
    high:i+1<curve.x.length?(x+curve.x[i+1])/2:x}));
}

export function conditionMapSample(rows, x, coordinate, tolerance=0) {
  if(!rows.length)return null;
  let row=rows[0];
  for(const candidate of rows)if(Math.abs(candidate.coordinate-coordinate)<Math.abs(row.coordinate-coordinate))row=candidate;
  if(!row.condition)return {row,point:{x:null,y:null},scope:{status:'not calculated',reasons:[]}};
  if(row.curve.unavailable_reason)return {row,point:{x:null,y:null},scope:{status:row.curve.unavailable_reason,reasons:[]}};
  if(row.curve.points_only && row.curve.x.length){
    const native=row.curve.x[nearestSample(row.curve.x,x)];
    if(Math.abs(native-x)>tolerance)return {row,point:{x:null,y:null},scope:{status:'not evaluated',reasons:[]}};
    x=native;
  }
  const point=curveSample(row.curve,[],x);
  return {row,point,scope:sampleScope(row.curve,[],point.x)};
}

export function mapColorLimits(rows, logarithmic=false) {
  const values=rows.flatMap(row=>row.curve.y).filter(v=>Number.isFinite(v) && (!logarithmic || v>0));
  if(!values.length)return logarithmic?[1,10]:[0,1];
  let low=Infinity,high=-Infinity;for(const value of values){low=Math.min(low,value);high=Math.max(high,value);}
  if(low===high){const pad=Math.abs(low)*.05 || 1;low=logarithmic?low/2:low-pad;high=logarithmic?high*2:high+pad;}
  return [low,high];
}

function mapColor(value, limits, logarithmic, diverging, period) {
  const [low,high]=limits,clamp=x=>Math.max(0,Math.min(1,x));
  if(period){const amount=((value+period/2)%period+period)%period/period*4;
    const palette=[[48,100,161],[145,77,145],[237,166,57],[29,132,115],[48,100,161]],i=Math.floor(amount),t=amount-i;
    return 'rgb('+palette[i].map((v,k)=>Math.round(v+(palette[i+1][k]-v)*t)).join(',')+')';}
  let amount=logarithmic?(Math.log(value)-Math.log(low))/(Math.log(high)-Math.log(low)):(value-low)/(high-low);
  let a=[239,245,243],b=[7,101,94];
  if(diverging){const extent=Math.max(Math.abs(low),Math.abs(high));amount=Math.abs(value)/extent;b=value<0?[48,100,161]:[174,53,61];}
  return 'rgb('+a.map((v,i)=>Math.round(v+(b[i]-v)*clamp(amount))).join(',')+')';
}

function spectrumMap(host, {dataset,panel,window:range,selection,onKeyboardSelect}) {
  const fig=element('figure',undefined,'mc-figure sc-map-figure');
  const title=element('figcaption',panel.title+' · condition map');
  const surface=element('div',undefined,'sc-map-surface'),canvas=element('canvas');
  canvas.setAttribute('aria-hidden','true');
  const svg=svgNode('svg',{role:'application',tabindex:0,
    'aria-label':panel.title+' condition map. Arrow keys select native samples and calculated conditions.'});
  surface.append(canvas,svg);fig.append(title,surface);host.append(fig);
  const {rows,numeric,axis}=conditionMapRows(dataset,panel,selection.condition_axis);
  const w=Math.max(170,Math.round(host.clientWidth || 360)),h=285;
  const left=w<300?49:65,right=w-15,top=26,bottom=h-39;
  const sx=spectrumScale(range[0],range[1],left,right,dataset.axis.scale==='log' && range[0]>0);
  const ylo=numeric?rows[0].coordinate:-.5,yhi=numeric?rows.at(-1).coordinate:rows.length-.5;
  const sy=spectrumScale(ylo,yhi,bottom,top);
  const logarithmic=panel.scale==='log',period=panel.cyclic_period,limits=period?[-period/2,period/2]:mapColorLimits(rows,logarithmic),diverging=panel.reference===0 && !logarithmic && !period;
  const hasValues=rows.some(row=>row.curve.y.some(y=>Number.isFinite(y) && (!logarithmic || y>0)));
  const dpr=Math.min(2,window.devicePixelRatio || 1);
  canvas.width=Math.round(w*dpr);canvas.height=Math.round(h*dpr);
  const ctx=canvas.getContext('2d');ctx.scale(dpr,dpr);ctx.clearRect(0,0,w,h);
  ctx.fillStyle='#edf1f1';ctx.fillRect(left,top,right-left,bottom-top);
  const sampledPoints=[];
  for(let i=0;i<rows.length;i++){
    const row=rows[i],a=numeric?(i?(rows[i-1].coordinate+row.coordinate)/2:ylo):i-.5;
    const b=numeric?(i+1<rows.length?(row.coordinate+rows[i+1].coordinate)/2:yhi):i+.5;
    row.bounds=[a,b];
    if(panel.sampled_only || row.curve.x.length===1){
      row.curve.x.forEach((x,j)=>{const y=row.curve.y[j];if(x>=range[0] && x<=range[1] && Number.isFinite(y) && (!logarithmic || y>0))sampledPoints.push({x:sx.forward(x),y:sy.forward(row.coordinate),fill:mapColor(y,limits,logarithmic,diverging,period)});});
      continue;
    }
    for(const cell of mapNativeCells(row.curve)){
      if(!Number.isFinite(cell.y) || (logarithmic && cell.y<=0) || cell.high<range[0] || cell.low>range[1])continue;
      // Shared boundaries must land on the same device pixel. Fractional canvas
      // rectangles otherwise leave pale seams that resemble spectral structure.
      const pixel=value=>Math.round(value*dpr)/dpr;
      const x1=pixel(sx.forward(Math.max(cell.low,range[0]))),x2=pixel(sx.forward(Math.min(cell.high,range[1])));
      if(!Number.isFinite(x1)||!Number.isFinite(x2))continue;
      const y1=pixel(sy.forward(b)),y2=pixel(sy.forward(a));
      if(x2<=x1 || y2<=y1)continue;
      ctx.fillStyle=mapColor(cell.y,limits,logarithmic,diverging,period);
      ctx.fillRect(x1,y1,x2-x1,y2-y1);
    }
  }
  svg.setAttribute('viewBox',[0,0,w,h].join(' '));
  const add=(tag,attrs,text)=>{const node=svgNode(tag,attrs,text);svg.append(node);return node;};
  for(const point of sampledPoints)add('circle',{class:'sc-map-native-point',cx:point.x,cy:point.y,r:3.5,fill:point.fill,stroke:'#203d3a','stroke-width':.7});
  add('path',{d:`M${left},${top}V${bottom}H${right}`,fill:'none',stroke:'#8b9f99'});
  for(const tick of spectrumXTicks(range[0],range[1],Math.max(2,Math.floor((right-left)/75)),dataset.axis.scale==='log')){
    const px=sx.forward(tick);add('line',{x1:px,x2:px,y1:bottom,y2:bottom+3,stroke:'#738986'});
    add('text',{x:px,y:bottom+16,'text-anchor':'middle'},format(tick,5));
  }
  const stride=Math.max(1,Math.ceil(rows.length/9));
  const labelRows=rows.filter((row,i)=>!(i%stride) || row.id===selection.series_id || i===rows.length-1);
  const priority=row=>row.id===selection.series_id?0:(row===rows[0] || row===rows.at(-1)?1:2);
  const labelPositions=[];
  labelRows.sort((a,b)=>priority(a)-priority(b)).forEach(row=>{
    const py=sy.forward(row.coordinate);
    if(labelPositions.some(y=>Math.abs(y-py)<14))return;
    labelPositions.push(py);
    const label=numeric?format(row.displayValue,5):String(row.value ?? row.label);
    add('text',{x:left-5,y:sy.forward(row.coordinate)+3,'text-anchor':'end'},label);
  });
  add('text',{x:4,y:12,'font-weight':600},(axis?.label || 'Condition')+(axis?.unit?' ['+axis.unit+']':''));
  add('text',{x:right,y:h-2,'text-anchor':'end'},dataset.axis.label+' ['+dataset.axis.unit+']');
  for(const guide of dataset.map_guides || []){
    // Declared physical-coordinate guide, clipped without changing native samples.
    const points=[ylo,yhi].map(y=>({x:guide.x_per_y*y+guide.x_offset,y}));
    if(guide.x_per_y)for(const x of range){const y=(x-guide.x_offset)/guide.x_per_y;if(y>=ylo && y<=yhi)points.push({x,y});}
    const inside=points.filter(p=>p.x>=range[0] && p.x<=range[1]).sort((a,b)=>a.y-b.y);
    if(inside.length>=2){const a=inside[0],b=inside.at(-1);const line=add('line',{class:'sc-map-guide',x1:sx.forward(a.x),x2:sx.forward(b.x),y1:sy.forward(a.y),y2:sy.forward(b.y),stroke:'#203d3a','stroke-dasharray':'5 4','pointer-events':'none'});line.append(svgNode('title',{},guide.label));}
  }
  const band=add('rect',{class:'sc-map-selection',x:left,width:right-left,fill:'none',stroke:'#254c48','stroke-width':1.5,'pointer-events':'none'});
  const ruler=add('line',{class:'sc-map-cursor',y1:top,y2:bottom,stroke:'#203d3a','stroke-dasharray':'2 3','pointer-events':'none'});
  function mark(current){
    const row=rows.find(r=>r.id===current.series_id);
    band.setAttribute('visibility',row?'visible':'hidden');
    if(row){band.setAttribute('y',sy.forward(row.bounds[1]));band.setAttribute('height',sy.forward(row.bounds[0])-sy.forward(row.bounds[1]));}
    const px=sx.forward(current.x);ruler.setAttribute('x1',px);ruler.setAttribute('x2',px);
    ruler.setAttribute('visibility',px>=left && px<=right?'visible':'hidden');
  }
  mark(selection);
  const legend=element('div',undefined,'sc-map-key');
  const low=element('span',format(diverging?-Math.max(...limits.map(Math.abs)):limits[0],5));
  const high=element('span',format(diverging?Math.max(...limits.map(Math.abs)):limits[1],5));
  const ramp=element('i');ramp.style.background=period?'linear-gradient(to right,#3064a1,#914d91,#eda639,#1d8473,#3064a1)':diverging?'linear-gradient(to right,#3064a1,#eff5f3,#ae353d)':'linear-gradient(to right,#eff5f3,#07655e)';
  legend.append(element('span',panel.y_label+(panel.unit && panel.unit!=='1'?' ['+panel.unit+']':'')+(logarithmic?' · log':'')));
  if(hasValues)legend.append(low,ramp,high);else legend.append(element('span','No finite samples for this quantity.'));
  fig.append(legend);
  for(const guide of dataset.map_guides || [])fig.append(element('p','Dashed guide: '+guide.label,'sc-panel-note mc-muted'));
  fig.append(element('p',panel.sampled_only?'Only calculated sample points are shown. Blank regions have not been evaluated.':'Each band is a calculated condition. Blank cells are missing, outside the native scan, or not calculated.','sc-panel-note mc-muted'));
  const tip=element('div',undefined,'mc-tooltip');tip.hidden=true;fig.append(tip);
  svg.addEventListener('pointermove',event=>{
    if(event.buttons){tip.hidden=true;return;}
    const p=spectrumPointer(svg,event,w,h);
    if(p.x<left||p.x>right||p.y<top||p.y>bottom){tip.hidden=true;return;}
    const sample=conditionMapSample(rows,sx.inverse(p.x),sy.inverse(p.y),Math.abs(sx.inverse(p.x+5)-sx.inverse(p.x-5))/2);
    tip.replaceChildren(element('strong',sample.row.label),element('div',sample.point.x==null?sample.scope.status:
      format(sample.point.x,7)+' '+dataset.axis.unit+' · '+panel.y_label+' '+readable(sample.point.y,panel.unit==='1'?'':panel.unit)));
    if(sample.scope.status && sample.scope.status!=='available')tip.append(element('small',sample.scope.status));
    sample.scope.reasons.forEach(reason=>tip.append(element('small',reason)));
    tip.hidden=false;
    const box=surface.getBoundingClientRect(),tipWidth=tip.offsetWidth;
    tip.style.left=Math.max(0,Math.min(w-tipWidth,p.x+12))+'px';
    tip.style.top=(Math.max(35,Math.min(h-70,(event.clientY-box.top)+20)))+'px';
  });
  svg.addEventListener('pointerleave',()=>{tip.hidden=true;});
  svg.addEventListener('keydown',event=>{
    if(!['ArrowLeft','ArrowRight','ArrowUp','ArrowDown','Home','End'].includes(event.key))return;
    event.preventDefault();
    const available=rows.filter(row=>row.condition);
    let row=available.find(row=>row.id===selection.series_id) || available[0];
    if(!row)return;
    if(event.key==='ArrowUp'||event.key==='ArrowDown')row=available[Math.max(0,Math.min(available.length-1,available.indexOf(row)+(event.key==='ArrowUp'?1:-1)))];
    const visible=row.curve.x.filter(x=>x>=range[0] && x<=range[1]);
    if(!visible.length)return;
    let index=nearestSample(visible,selection.x);
    if(event.key==='ArrowLeft')index=Math.max(0,index-1);
    if(event.key==='ArrowRight')index=Math.min(visible.length-1,index+1);
    if(event.key==='Home')index=0;if(event.key==='End')index=visible.length-1;
    onKeyboardSelect(row.id,visible[index]);
  });
  svg._navigation={kind:'condition-map',rows,w,h,left,right,top,bottom,toFrequency:sx.inverse,toCondition:sy.inverse,mark};
  svg.style.touchAction='none';
  return svg;
}
