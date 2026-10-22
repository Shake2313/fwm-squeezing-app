// Keep the preview's interpolation/affine rules independently testable.
export function interpolateCurve(c,x,linear=false) {
  if (x<c.x[0]-1e-9 || x>c.x.at(-1)+1e-9) return NaN;
  let lo=0, hi=c.x.length-1;
  while (hi-lo>1) {const mid=(lo+hi)>>1; if(c.x[mid]<=x) lo=mid; else hi=mid;}
  if(Math.abs(x-c.x[lo])<1e-9) return c.y[lo] ?? NaN;
  if(Math.abs(x-c.x[hi])<1e-9) return c.y[hi] ?? NaN;
  if(c.y[lo]==null || c.y[hi]==null) return NaN;
  const t=(x-c.x[lo])/(c.x[hi]-c.x[lo]);
  return linear ? 10*Math.log10((1-t)*10**(c.y[lo]/10)+t*10**(c.y[hi]/10))
                : (1-t)*c.y[lo]+t*c.y[hi];
}

export function alignPoint(x,y,pivot,original,values) {
  return [pivot+(x-pivot-original.shift)*values.multiplier/original.multiplier+values.shift,
          y==null ? NaN : (y-original.signal_offset)/original.signal_scale*values.signal_scale+values.signal_offset];
}

export default function(component) {
  const {parentElement, data: d, setStateValue} = component;
  const root = parentElement.querySelector('.lp');
  if (!root) return;
  root.removeAttribute('aria-busy');
  const stage = root.querySelector('.lp-stage'), svg = root.querySelector('svg');
  const status = root.querySelector('output'), controls = root.querySelector('.lp-controls');
  root.querySelector('img').src = d.image;
  // Match the plot card's viewport cap, with room for the preview readout.
  // Constrain the stage, not just the image, so SVG and raster always agree.
  stage.style.width = `min(100%, calc((100vh - 18rem) * ${d.width/d.height}))`;
  svg.setAttribute('viewBox', `0 0 ${d.width} ${d.height}`);
  const values = {...d.original}, ns = 'http://www.w3.org/2000/svg';
  const state = root._live || (root._live = {seq: 0});
  // Rerenders/unmounts must not accumulate window listeners or animation frames.
  state.cleanup?.();
  svg.replaceChildren(); controls.replaceChildren();
  const make = (name, attrs, parent=svg) => {
    const el = document.createElementNS(ns, name);
    Object.entries(attrs).forEach(([k,v]) => el.setAttribute(k, v));
    parent.append(el); return el;
  };
  const px = (a,x) => a.left + (x-a.xlim[0])/(a.xlim[1]-a.xlim[0])*a.width;
  const py = (a,y) => {
    const f = a.log ? Math.log10 : (v) => v;
    if (y == null || !Number.isFinite(y) || (a.log && y <= 0)) return NaN;
    return a.top + (f(a.ylim[1])-f(y))/(f(a.ylim[1])-f(a.ylim[0]))*a.height;
  };
  const defs=make('defs',{}), groups=d.axes.map((a,i)=> {
    const id=`clip-${i}`, clip=make('clipPath',{id},defs);
    make('rect',{x:a.left,y:a.top,width:a.width,height:a.height},clip);
    return make('g',{'clip-path':`url(#${id})`});
  });
  const inputs = {};
  let pending=false, frame=0, drag=null;
  const fmt=(v,digits=6)=> Number.isFinite(v) ? Number(v.toPrecision(digits)).toString() : 'unavailable';
  const cursors = d.kind==='cursor' ? groups.map((g)=> ({
    line:make('line',{stroke:'crimson','stroke-width':1.5,'stroke-dasharray':'6 4'},g),
    point:make('circle',{r:4,fill:'crimson'},g)
  })) : [];
  const paths = d.kind==='alignment' ? d.traces.map((t)=>make('path',{
    fill:'none',stroke:t.color,'stroke-width':t.dashed?1.3:1.8,
    'stroke-opacity':t.alpha,'stroke-dasharray':t.dashed?'2 3':'none',
    'vector-effect':'non-scaling-stroke'
  },groups[0])) : [];
  function paint(preview=false) {
    if(d.kind==='alignment') {
      d.traces.forEach((t,j)=> {
        let path='', start=true;
        t.x.forEach((x,i)=> {
          const [xx,yy]=alignPoint(x,t.y[i],d.pivot,d.original,values);
          const u=px(d.axes[0],xx), v=t.y[i]==null ? NaN : py(d.axes[0],yy);
          if(!Number.isFinite(u)||!Number.isFinite(v)){start=true;return;}
          path+=(start?'M':'L')+u.toFixed(3)+' '+v.toFixed(3); start=false;
        });
        paths[j].setAttribute('d',path);
      });
      status.textContent=(preview?'Preview · ':'')+`X shift ${fmt(values.shift)} · X scale ×${fmt(values.multiplier)}`;
    } else {
      const ys=d.curves.map((c,i)=>interpolateCurve(c,values.tpd,i===1&&d.linear_noise));
      cursors.forEach((c,i)=> {
        const a=d.axes[i], x=px(a,values.tpd), y=py(a,ys[i]);
        Object.entries({x1:x,x2:x,y1:a.top,y2:a.top+a.height}).forEach(([k,v])=>c.line.setAttribute(k,v));
        c.point.setAttribute('cx',x); c.point.setAttribute('cy',Number.isFinite(y)?y:0);
        c.point.style.display=Number.isFinite(y)?'':'none';
      });
      status.textContent=(preview?'Preview · ':'')+`δ ${fmt(values.tpd)} MHz · Gain ${d.gain_digits===3?'≈ ':''}${fmt(ys[0],d.gain_digits||6)} · Noise ${fmt(ys[1])} dB`;
    }
  }
  function schedule() {
    if(!frame) frame=requestAnimationFrame(()=>{
      frame=0;paint(Object.keys(values).some(k=>values[k]!==d.original[k]));
    });
  }
  function commit() {
    if(frame){cancelAnimationFrame(frame);frame=0;}
    if(!Object.keys(values).some(k=>values[k]!==d.original[k])) {paint();return;}
    paint(true); pending=true;
    root.setAttribute('aria-busy','true');
    for(const pair of Object.values(inputs)) pair.forEach(el=>el.disabled=true);
    status.textContent+=' · Synchronizing…';
    state.seq++;
    setStateValue('gesture',JSON.stringify({n:state.seq,context:d.context,values}));
  }
  function cancel() {
    Object.assign(values,d.original);
    for(const [key,pair] of Object.entries(inputs)) pair.forEach(el=>el.value=values[key]);
    if(frame){cancelAnimationFrame(frame);frame=0;}
    paint();
  }
  for(const spec of d.controls) {
    const label=document.createElement('label'), text=document.createElement('span');
    text.textContent=spec.label; label.append(text);
    const number=document.createElement('input'), range=document.createElement('input');
    for(const el of [number,range]) {
      el.type=el===number?'number':'range'; el.min=spec.min; el.max=spec.max; el.step=spec.step;
      el.value=values[spec.name]; el.setAttribute('aria-label',spec.label);
      el.oninput=()=> {
        const v=el.valueAsNumber;
        if(pending || !Number.isFinite(v)||v<spec.min||v>spec.max) return;
        values[spec.name]=v; (el===number?range:number).value=v; schedule();
      };
      el.onchange=commit;
      el.onkeydown=e=>{if(e.key==='Escape'){e.preventDefault();cancel();el.blur();}};
      el.onpointercancel=cancel;
    }
    inputs[spec.name]=[number,range]; label.append(number,range); controls.append(label);
  }
  const inAxis=(e)=> {
    const r=svg.getBoundingClientRect(), x=(e.clientX-r.left)*d.width/r.width, y=(e.clientY-r.top)*d.height/r.height;
    return d.axes.find(a=>x>=a.left&&x<=a.left+a.width&&y>=a.top&&y<=a.top+a.height);
  };
  function move(e) {
    if(!drag || pending) return;
    const r=svg.getBoundingClientRect(), a=drag.axis;
    const dx=(e.clientX-drag.x)*d.width/r.width/a.width*(a.xlim[1]-a.xlim[0]);
    if(d.kind==='alignment') {
      const spec=d.controls.find(s=>s.name==='shift');
      values.shift=Math.min(spec.max,Math.max(spec.min,drag.value+dx));
      inputs.shift.forEach(el=>el.value=values.shift);
    } else {
      const x=(e.clientX-r.left)*d.width/r.width;
      values.tpd=Math.min(a.xlim[1],Math.max(a.xlim[0],a.xlim[0]+(x-a.left)/a.width*(a.xlim[1]-a.xlim[0])));
    }
    schedule();
  }
  svg.onpointerdown=e=> {
    if(e.button!==0||pending) return;
    const a=inAxis(e); if(!a)return;
    svg.setPointerCapture(e.pointerId);
    drag={x:e.clientX,value:values.shift,axis:a};
    if(d.kind==='cursor')move(e);
  };
  svg.onpointermove=move;
  svg.onpointerup=e=> {
    if(!drag)return; move(e);drag=null;
    if(svg.hasPointerCapture(e.pointerId))svg.releasePointerCapture(e.pointerId);
    commit();
  };
  svg.onpointercancel=()=>{drag=null;cancel();};
  const onPreview=e=> {
    if(e.detail?.key!==d.event_key||pending)return;
    values.tpd=e.detail.value;
    if(e.detail.cancel){cancel();return;}
    schedule();
  };
  window.addEventListener('gabes-preview',onPreview);
  const cleanup=()=>{window.removeEventListener('gabes-preview',onPreview);if(frame)cancelAnimationFrame(frame);};
  state.cleanup=cleanup;
  paint();
  return cleanup;
}
