// All interactions read the solver's original transmission and optical-depth samples.
const SVG_NS = 'http://www.w3.org/2000/svg';
const ink = {off:'#4c606c', on:'#07655e', cursor:'#192f39'};

export function nearestSample(grid, value) {
  let lo = 0, hi = grid.length - 1;
  while (lo < hi) {
    const mid = (lo + hi) >> 1;
    if (grid[mid] < value) lo = mid + 1; else hi = mid;
  }
  return lo && Math.abs(value - grid[lo-1]) <= Math.abs(grid[lo] - value) ? lo-1 : lo;
}

function format(value, digits=5) {
  if (!Number.isFinite(value)) return 'unavailable';
  if (value === 0) return '0';
  return Number(value.toPrecision(digits)).toString().replace('-', '−');
}

function element(tag, text, cls) {
  const node = document.createElement(tag);
  if (text !== undefined) node.textContent = text;
  if (cls) node.className = cls;
  return node;
}

function svgNode(tag, attributes, text) {
  const node = document.createElementNS(SVG_NS, tag);
  for (const [name, value] of Object.entries(attributes)) node.setAttribute(name, value);
  if (text !== undefined) node.textContent = text;
  return node;
}

function ticks(lo, hi, target) {
  const span = hi - lo;
  if (span <= 0) return [lo];
  const rough = span / target, scale = 10 ** Math.floor(Math.log10(rough));
  const step = [1, 2, 2.5, 5, 10].find(n => n >= rough / scale) * scale;
  const values = [];
  for (let i = Math.ceil(lo / step - 1e-9); i <= Math.floor(hi / step + 1e-9); i++)
    values.push(i * step);
  return values;
}

function valueRange(curves, key, first, last, initialMax) {
  let lo=0, hi=initialMax;
  for (const curve of curves) {
    for (let i=first; i<=last; i++) {
      const y=curve[key][i];
      if (Number.isFinite(y)) {lo=Math.min(lo,y);hi=Math.max(hi,y);}
    }
  }
  return [lo, hi+(hi-lo || 1)*.04];
}

function curveSample(curve, grid, frequency) {
  const x = curve.x || grid;
  if (!x.length || frequency < x[0] || frequency > x[x.length-1]) return {x:null,y:null};
  const index = nearestSample(x, frequency);
  return {x:x[index],y:curve.y[index]};
}

export function sampleScope(curve, grid, x) {
  if(x==null)return {status:'outside scan',reasons:[]};
  const index=nearestSample(curve.x || grid,x);
  return {status:curve.sample_status?.[index] || curve.observable_status?.[index] || '',reasons:curve.sample_reasons?.[index] || []};
}

export function displaySample(curve, point) {
  return {x:point.x,y:Number.isFinite(point.y)?point.y+(curve.display_offset || 0):point.y};
}

export function sampleInsidePlot(point, geometry) {
  if(point.x==null || !Number.isFinite(point.y))return false;
  const x=geometry.sx(point.x),y=geometry.sy(point.y);
  return x>=geometry.left && x<=geometry.right && y>=geometry.top && y<=geometry.bottom;
}

export function spectrumXTicks(lo, hi, target, logarithmic=false) {
  if(!logarithmic || lo<=0)return ticks(lo,hi,target);
  const result=[];
  for(let decade=Math.floor(Math.log10(lo));decade<=Math.ceil(Math.log10(hi));decade++)
    for(const multiplier of [1,2,5]){
      const value=multiplier*10**decade;
      if(value>=lo && value<=hi)result.push(value);
    }
  if(result.length>Math.max(6,target*2)){
    const major=result.filter(value=>Math.abs(Math.log10(value)-Math.round(Math.log10(value)))<1e-9);
    if(major.length>=2)return major.filter((_,i)=>i%Math.max(1,Math.ceil(major.length/(target+1)))===0);
  }
  return result.length>=2 ? result : ticks(lo,hi,target);
}

export function sampleWarningIntervals(curves) {
  const ranges=[];
  for(const curve of curves)for(let i=0;i<(curve.sample_status?.length || 0);i++){
    if(!curve.sample_status[i] || curve.sample_status[i]==='available')continue;
    const x=curve.x;
    ranges.push([i?(x[i-1]+x[i])/2:x[i],i+1<x.length?(x[i]+x[i+1])/2:x[i]]);
  }
  ranges.sort((a,b)=>a[0]-b[0]);
  const merged=[];
  for(const range of ranges){const prior=merged.at(-1);if(prior && range[0]<=prior[1])prior[1]=Math.max(prior[1],range[1]);else merged.push([...range]);}
  return merged;
}

export function spectrumScale(lo, hi, left, right, logarithmic=false) {
  const transform=logarithmic ? Math.log10 : x=>x;
  const inverse=logarithmic ? x=>10**x : x=>x;
  const a=transform(lo),span=transform(hi)-a || 1;
  return {forward:x=>left+(transform(x)-a)/span*(right-left),
    inverse:pixel=>inverse(a+(pixel-left)/(right-left)*span)};
}

function sampleInspector(fig, svg, spec, geometry) {
  const {left, right, top, bottom, sx, sy, w, h, toFrequency} = geometry;
  const overlay = svgNode('g', {'aria-hidden':'true', class:'mc-hover-overlay', visibility:'hidden'});
  const ruler = svgNode('line', {y1:top,y2:bottom,stroke:'#506d78','stroke-dasharray':'4 3','stroke-width':1});
  overlay.append(ruler);
  const dots = spec.series.map(curve => {
    const dot = svgNode('circle', {r:3.8,fill:curve.color,stroke:'white','stroke-width':1.3});
    overlay.append(dot);
    return dot;
  });
  svg.append(overlay);
  const tip = element('div', undefined, 'mc-tooltip');
  tip.hidden = true;
  tip.setAttribute('aria-hidden','true');
  const title = element('strong');
  tip.append(title);
  const values = spec.series.map(curve => {
    const row = element('div', undefined, 'mc-tip-row');
    const sample = element('i', undefined, curve.dash ? 'off' : 'on');
    sample.style.borderColor=curve.color;
    const value = element('b');
    row.append(sample,element('span',curve.label),value); tip.append(row);
    return value;
  });
  const spoken = element('span', undefined, 'mc-sr');
  spoken.setAttribute('role','status'); fig.append(tip, spoken);
  let index = spec.index;
  let prior = -1;
  const [lo,hi] = spec.xlim;
  let first = nearestSample(spec.grid, lo), last = nearestSample(spec.grid, hi);
  if (first<last && spec.grid[first]<lo) first++;
  if (last>first && spec.grid[last]>hi) last--;
  function hide() { tip.hidden=true; overlay.setAttribute('visibility','hidden'); prior=-1; }
  function show(i, pointerY, announce=false) {
    index = Math.max(first,Math.min(last,i));
    const x = spec.grid[index], px = sx(x);
    if (index !== prior || tip.hidden) {
      const unitScale = spec.divisor ?? (spec.unit === 'GHz' ? 1000 : 1);
      const label=spec.tickLabels?.[Math.round(x)];
      title.textContent = (label || format(x/unitScale,8)+' '+spec.unit)+' · '+spec.yLabel;
      ruler.setAttribute('x1',px);ruler.setAttribute('x2',px);
      spec.series.forEach((curve,k)=>{
        const point=curveSample(curve,spec.grid,x), y=point.y,shown=displaySample(curve,point);
        const scope=sampleScope(curve,spec.grid,point.x);
        values[k].textContent=point.x==null ? 'outside scan' : format(y,6)+(Number.isFinite(y) && spec.yUnit && spec.yUnit!=='1' ? ' '+spec.yUnit : '')+(curve.x ? ' @ '+format(point.x,7)+' '+spec.unit : '')+(scope.status && scope.status!=='available'?' · '+scope.status:'');
        values[k].parentElement.title=scope.reasons.join('; ');
        const inside=sampleInsidePlot(shown,geometry);
        dots[k].setAttribute('visibility',inside?'inherit':'hidden');
        if (inside) {dots[k].setAttribute('cx',sx(point.x));dots[k].setAttribute('cy',sy(shown.y));}
      });
      prior=index;
    }
    tip.hidden=false;overlay.setAttribute('visibility','visible');
    const sr=svg.getBoundingClientRect(), fr=fig.getBoundingClientRect();
    const ax=sr.left-fr.left+px*sr.width/w;
    const ay=pointerY == null ? sr.top-fr.top+top*sr.height/h : pointerY-fr.top;
    const tw=tip.offsetWidth, th=tip.offsetHeight;
    const xLeft=ax+12+tw<=fr.width?ax+12:ax-tw-12;
    tip.style.left=Math.max(0,Math.min(fr.width-tw,xLeft))+'px';
    tip.style.top=Math.max(sr.top-fr.top,Math.min(fr.height-th,ay-th/2))+'px';
    if (announce) spoken.textContent = title.textContent+'. '+spec.series.map((c,k)=>c.label+' '+values[k].textContent).join(', ');
  }
  svg.addEventListener('pointermove',event=>{
    if(event.buttons || event.pointerType==='touch') {hide(); return;}
    const box=svg.getBoundingClientRect(), px=(event.clientX-box.left)*w/box.width, py=(event.clientY-box.top)*h/box.height;
    if(px<left || px>right || py<top || py>bottom) {hide();return;}
    show(nearestSample(spec.grid,toFrequency(px)),event.clientY);
  });
  svg.addEventListener('pointerdown',hide);
  svg.addEventListener('pointerleave',hide);
  if(spec.keyboardInspector) {
    svg.setAttribute('tabindex','0');
    svg.setAttribute('aria-label',spec.description+' Use Left or Right to inspect a sample.');
    svg.addEventListener('focus',()=>{if(svg.matches(':focus-visible')) show(spec.index,null,true);});
    svg.addEventListener('blur',hide);
    svg.addEventListener('keydown',event=>{
      if(event.key==='Escape') {hide();return;}
      if(!['ArrowLeft','ArrowRight','Home','End'].includes(event.key)) return;
      event.preventDefault();
      const next=event.key==='Home'?first:event.key==='End'?last:index+(event.key==='ArrowLeft'?-1:1);
      show(next,null,true);
    });
  }
}

function spectrum(host, spec) {
  const fig = element('figure', undefined, 'mc-figure');
  const caption = element('figcaption', spec.title);
  const svg = svgNode('svg', {role:'img', 'aria-label':spec.description});
  fig.append(caption, svg); host.append(fig);
  const w = Math.max(170, Math.round(host.clientWidth || 360)), h = spec.height || (spec.tickLabels ? 226 : spec.compact ? 145 : 181);
  const left = spec.stacked?58:spec.tickLabels ? Math.min(92,w/3) : w < 300 ? 44 : 52, right = w - 13, top = 23, bottom = h - (spec.tickLabels ? 80 : 36);
  const [xlo,xhi] = spec.xlim, [ylo,yhi] = spec.ylim;
  const scale=spectrumScale(xlo,xhi,left,right,spec.xScale==='log' && xlo>0);
  const sx = scale.forward;
  const yScale=spectrumScale(ylo,yhi,bottom,top,spec.yScale==='log' && ylo>0);
  const sy=yScale.forward;
  const add = (tag, attrs, text) => svg.append(svgNode(tag, attrs, text));
  svg.setAttribute('viewBox', [0,0,w,h].join(' '));
  if (spec.band) {
    const a = Math.max(left, sx(spec.band[0])), b = Math.min(right, sx(spec.band[1]));
    add('rect', {class:'mc-nav-window',x:a,y:top,width:b-a,height:bottom-top,fill:'#dceee4'});
    for (const [side,x] of [['left',a],['right',b]])
      add('path',{class:'mc-nav-edge', 'data-side':side,
        d:`M0,${top}v${bottom-top}m-2,0h4`,transform:`translate(${x} 0)`,
        fill:'none',stroke:'#769886','stroke-width':1});
  }
  for(const [a,b] of sampleWarningIntervals(spec.series)){
    const low=Math.max(a,xlo),high=Math.min(b,xhi);
    if(low<high)add('rect',{class:'mc-sample-warning',x:sx(low),y:top,width:sx(high)-sx(low),height:bottom-top,fill:'#ad7023','fill-opacity':.16});
  }
  for (const y of spec.stacked?[]:spectrumXTicks(ylo,yhi,3,spec.yScale==='log')) {
    const py = sy(y);
    add('line', {x1:left, x2:right, y1:py, y2:py, stroke:y === 0 ? '#879894' : '#dbe3e1'});
    add('text', {x:left-5, y:py+3, 'text-anchor':'end'}, format(y,4));
  }
  add('path', {d:'M'+left+' '+top+'V'+bottom+'H'+right, fill:'none', stroke:'#9cadaa'});
  const divisor = spec.divisor ?? (spec.unit === 'GHz' ? 1000 : 1);
  const xticks=spec.tickLabels ? spec.grid.filter(x=>x>=xlo && x<=xhi).filter((_,i)=>i%Math.max(1,Math.ceil((xhi-xlo+1)/((right-left)/35)))===0)
    : spectrumXTicks(xlo/divisor,xhi/divisor,Math.max(2,Math.floor((right-left)/72)),spec.xScale==='log');
  for (const tick of xticks) {
    const px = sx(tick*divisor);
    add('line',{x1:px,x2:px,y1:bottom,y2:bottom+3,stroke:'#738986'});
    add('text',{x:px,y:bottom+15,'text-anchor':spec.tickLabels?'end':'middle',
      ...(spec.tickLabels?{transform:`rotate(-45 ${px} ${bottom+15})`}:{})},spec.tickLabels?.[Math.round(tick)] || format(tick,5));
  }
  add('text',{x:4,y:12,'font-weight':'600'},spec.yLabel);
  add('text',{x:right,y:h-2,'text-anchor':'end'},(spec.xLabel || 'Relative frequency')+' ('+spec.unit+')');
  if(spec.stacked){
    for(const curve of spec.series){
      const curveGrid=curve.x || spec.grid;
      const baseline=curve.y.find((value,index)=>Number.isFinite(value) && curveGrid[index]>=xlo && curveGrid[index]<=xhi);
      if(baseline===undefined)continue;
      const py=sy(baseline+(curve.display_offset || 0));
      add('line',{x1:left,x2:right,y1:py,y2:py,stroke:'#dbe3e1','stroke-width':.7});
      add('text',{x:left-5,y:py+3,'text-anchor':'end',fill:curve.color},curve.shortLabel || curve.label);
    }
    const amount=spec.stacked.scaleStep,py=sy(ylo+amount),base=sy(ylo);
    add('path',{d:`M${right-6},${py}h6m-3,0v${base-py}m-3,0h6`,fill:'none',stroke:'#526b70'});
    add('text',{x:right-12,y:(py+base)/2+3,'text-anchor':'end'},format(amount,3)+' '+(spec.stacked.scaleUnit || ''));
  }
  if (Number.isFinite(spec.reference) && spec.reference>=ylo && spec.reference<=yhi) {
    add('line',{x1:left,x2:right,y1:sy(spec.reference),y2:sy(spec.reference),stroke:'#4b6068','stroke-dasharray':'4 3'});
    add('text',{x:right-2,y:sy(spec.reference)-4,'text-anchor':'end'},spec.referenceLabel || format(spec.reference)+' '+(spec.yUnit||''));
  }
  const groupedMarks=new Map();
  for(const m of spec.marks || []){
    const x=m.x ?? m.x_mhz,prior=groupedMarks.get(x);
    groupedMarks.set(x,prior?{...prior,label:prior.label+' / '+m.label}:{...m,x_mhz:x});
  }
  const marks = [...groupedMarks.values()].filter(m => m.x_mhz >= xlo && m.x_mhz <= xhi);
  let rightmostLabel = -Infinity;
  for (const m of marks) {
    const px = sx(m.x_mhz), co = m.kind === 'crossover_candidate';
    const line = svgNode('line',{x1:px,x2:px,y1:top+4,y2:bottom,stroke:co?'#8b6e48':'#7c8d88','stroke-dasharray':co?'2 4':'4 4','stroke-width':.8});
    line.append(svgNode('title',{},m.label)); svg.append(line);
    const abbreviation = co ? 'CO?' : String(m.label).replace(/^.*?\s(?=\d+\s*[→-])/, '');
    const labelWidth = Math.max(19,abbreviation.length*6);
    if (marks.length <= 5 && px - labelWidth / 2 > rightmostLabel + 2 && px + labelWidth/2 < w) {
      add('text',{x:px,y:top-2,'text-anchor':'middle','font-size':9},abbreviation);
      rightmostLabel = px + labelWidth/2;
    }
  }
  const grid = spec.grid;
  for (const [curveIndex,c] of spec.series.entries()) {
    const curveGrid=c.x || grid;
    let path = '', pen = false, previousY=null;
    for (let i = nearestSample(curveGrid,xlo); i < curveGrid.length && curveGrid[i] <= xhi; i++) {
      if (curveGrid[i] < xlo) continue;
      if(c.points_only || (c.cyclic_period && previousY!==null && Math.abs(c.y[i]-previousY)>c.cyclic_period/2))pen=false;
      previousY=c.y[i];
      const y=c.y[i];
      if (!Number.isFinite(y) || (spec.yScale==='log' && y<=0)) {pen=false; continue;}
      const shown=displaySample(c,{x:curveGrid[i],y});
      path += (pen?'L':'M')+sx(curveGrid[i]).toFixed(2)+','+sy(shown.y).toFixed(2);
      pen=true;
      if(spec.discrete || c.markers || curveGrid.length<=12){
        const px=sx(curveGrid[i]),py=sy(shown.y);
        if(c.point_kind==='optimum')add('path',{d:`M${px},${py-5}L${px+5},${py}L${px},${py+5}L${px-5},${py}Z`,fill:'white',stroke:c.color,'stroke-width':1.7});
        else add('circle',{cx:px,cy:py,r:c.point_kind==='selected'?5:spec.discrete?5-Math.min(curveIndex,3):3,
          fill:c.point_kind==='selected'?'white':spec.discrete?'none':c.color,
          ...((spec.discrete || c.point_kind==='selected')?{stroke:c.color,'stroke-width':c.point_kind==='selected'?2:1.3}:{})});
      }
    }
    if(!spec.discrete)add('path',{d:path,fill:'none',stroke:c.color,'stroke-width':1.9,'stroke-dasharray':c.dash || 'none','stroke-linejoin':'round'});
  }
  const px = sx(grid[spec.index]);
  if (spec.cursor !== false && px >= left && px <= right) {
    add('line',{class:'mc-nav-cursor',x1:px,x2:px,y1:top,y2:bottom,stroke:ink.cursor,'stroke-dasharray':'1.5 3'});
    for (const c of spec.series) {
      const point=curveSample(c,grid,grid[spec.index]),shown=displaySample(c,point);
      const sampleX=point.x==null?left:sx(point.x),sampleY=Number.isFinite(shown.y)?sy(shown.y):top;
      const dot=svgNode('circle',{class:'mc-nav-dot',cx:Number.isFinite(sampleX)?sampleX:left,cy:Number.isFinite(sampleY)?sampleY:top,r:2.8,fill:c.color,stroke:'white','stroke-width':.9});
      dot.setAttribute('visibility',sampleInsidePlot(shown,{left,right,top,bottom,sx,sy})?'visible':'hidden');svg.append(dot);
    }
  }
  sampleInspector(fig, svg, spec, {left,right,top,bottom,sx,sy,w,h,toFrequency:scale.inverse});
  svg._navigation={w,h,left,right,top,bottom,sx,sy,grid,series:spec.series,
    overview:Boolean(spec.band),toFrequency:scale.inverse};
  svg.style.touchAction='none';
}

function previewMarkers(svg, window, index) {
  const n=svg._navigation, px=n.sx(n.grid[index]);
  const visible=px>=n.left && px<=n.right;
  const line=svg.querySelector('.mc-nav-cursor');
  if(line) {
    line.setAttribute('x1',px);line.setAttribute('x2',px);
    line.setAttribute('visibility',visible?'visible':'hidden');
  }
  svg.querySelectorAll('.mc-nav-dot').forEach((dot,j)=>{
    const point=curveSample(n.series[j],n.grid,n.grid[index]),shown=displaySample(n.series[j],point),y=shown.y;
    const sampleX=point.x==null?n.left:n.sx(point.x),inside=visible&&sampleInsidePlot(shown,n);
    dot.setAttribute('cx',Number.isFinite(sampleX)?sampleX:n.left);
    dot.setAttribute('visibility',inside?'visible':'hidden');
    if(inside){dot.setAttribute('cy',n.sy(y));}
  });
  if(!n.overview) return;
  const a=n.sx(window[0]), b=n.sx(window[1]), band=svg.querySelector('.mc-nav-window');
  band.setAttribute('x',a);band.setAttribute('width',b-a);
  svg.querySelectorAll('.mc-nav-edge').forEach(edge=>edge.setAttribute('transform',
    `translate(${edge.dataset.side==='left'?a:b} 0)`));
}
