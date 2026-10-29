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

// One request is in flight at a time; newer completed gestures replace the
// queued request. The displayed preview may be newer still (an active drag).
export function sameValues(a,b) {
  return Object.keys(a).length===Object.keys(b).length &&
    Object.keys(a).every(k=>a[k]===b[k]);
}

export function receiveRender(state,data) {
  if(!state || state.context!==data.context) {
    return {context:data.context,seq:data.ack||0,values:{...data.original},
      data,inflight:null,queued:null,previewing:false,drag:null,frame:0};
  }
  state.data=data;
  state.seq=Math.max(state.seq,data.ack||0);
  if(state.inflight && data.ack>=state.inflight.n) state.inflight=null;
  if(!state.inflight && state.queued && sameValues(state.queued,data.original)) state.queued=null;
  if(!state.inflight && !state.queued && !state.previewing) state.values={...data.original};
  return state;
}

export function releaseGesture(state,intent=0) {
  state.previewing=false;
  state.queued={...state.values};
  state.queuedIntent=intent;
  return nextGesture(state);
}

export function nextGesture(state) {
  if(state.inflight || !state.queued) return null;
  if(sameValues(state.queued,state.data.original)) {state.queued=null;return null;}
  state.inflight={n:++state.seq,context:state.context,values:{...state.queued},intent:state.queuedIntent||0};
  state.queued=null;
  return state.inflight;
}

// Streamlit may replace the HTML wrapper while retaining its component key.
// Its JS module is cached by content hash; keep interactions outside that DOM.
const mountedInteractions=new Map();

export function interactionForComponent(key,data) {
  const state=receiveRender(mountedInteractions.get(key),data);
  mountedInteractions.set(key,state);
  return state;
}

export default function(component) {
  const {parentElement, data:d, setStateValue}=component;
  const root=parentElement.querySelector('.lp');
  if(!root)return;
  const mountKey=component.key;
  const previous=mountedInteractions.get(mountKey);
  previous?.cleanup?.();
  const state=root._live=interactionForComponent(mountKey,d);
  state.root=root;
  const stage=root.querySelector('.lp-stage'), svg=root.querySelector('svg');
  const status=root.querySelector('output'), controls=root.querySelector('.lp-controls');
  const img=root.querySelector('img');
  if(img.getAttribute('src')!==d.image)img.src=d.image;
  stage.style.width=`min(100%, calc(var(--g-plot-height, calc(100dvh - 22rem)) * ${d.width/d.height}))`;
  stage.style.aspectRatio=`${d.width} / ${d.height}`;
  svg.setAttribute('viewBox',`0 0 ${d.width} ${d.height}`);
  // The SVG element remains in place, retaining pointer capture across reruns.
  // Only its drawing children change. Active state lives in the module store.
  svg.replaceChildren();
  const ns='http://www.w3.org/2000/svg';
  const make=(name,attrs,parent=svg)=> {
    const el=document.createElementNS(ns,name);
    Object.entries(attrs).forEach(([k,v])=>el.setAttribute(k,v));
    parent.append(el);return el;
  };
  const px=(a,x)=>a.left+(x-a.xlim[0])/(a.xlim[1]-a.xlim[0])*a.width;
  const py=(a,y)=> {
    const f=a.log?Math.log10:(v)=>v;
    if(y==null||!Number.isFinite(y)||(a.log&&y<=0))return NaN;
    return a.top+(f(a.ylim[1])-f(y))/(f(a.ylim[1])-f(a.ylim[0]))*a.height;
  };
  const defs=make('defs',{}), groups=d.axes.map((a,i)=> {
    const id=`lp-${d.scheme||d.kind}-${i}`,clip=make('clipPath',{id},defs);
    make('rect',{x:a.left,y:a.top,width:a.width,height:a.height},clip);
    return make('g',{'clip-path':`url(#${id})`});
  });
  const inputs={};
  const fmt=(v,digits=6)=>Number.isFinite(v)?Number(v.toPrecision(digits)).toString():'unavailable';
  const cursors=d.kind==='cursor'?groups.map(g=>({
    line:make('line',{stroke:'crimson','stroke-width':1.5,'stroke-dasharray':'6 4'},g),
    point:make('circle',{r:4,fill:'crimson'},g)
  })):[];
  const paths=d.kind==='alignment'?d.traces.map(t=>make('path',{
    fill:'none',stroke:t.color,'stroke-width':t.dashed?1.3:1.8,
    'stroke-opacity':t.alpha,'stroke-dasharray':t.dashed?'2 3':'none',
    'vector-effect':'non-scaling-stroke'
  },groups[0])):[];

  function paint() {
    const values=state.values;
    const phase=state.previewing?'preview':state.inflight||state.queued?'syncing':'ready';
    let text;
    if(d.kind==='alignment') {
      d.traces.forEach((t,j)=> {
        let path='',start=true;
        t.x.forEach((x,i)=> {
          const [xx,yy]=alignPoint(x,t.y[i],d.pivot,d.original,values);
          const u=px(d.axes[0],xx),v=t.y[i]==null?NaN:py(d.axes[0],yy);
          if(!Number.isFinite(u)||!Number.isFinite(v)){start=true;return;}
          path+=(start?'M':'L')+u.toFixed(3)+' '+v.toFixed(3);start=false;
        });
        paths[j].setAttribute('d',path);
      });
      text=`X shift ${fmt(values.shift)} · X scale ×${fmt(values.multiplier)}`;
    } else {
      const ys=d.curves.map((c,i)=>interpolateCurve(c,values.tpd,i===1&&d.linear_noise));
      cursors.forEach((c,i)=> {
        const a=d.axes[i],x=px(a,values.tpd),y=py(a,ys[i]);
        Object.entries({x1:x,x2:x,y1:a.top,y2:a.top+a.height}).forEach(([k,v])=>c.line.setAttribute(k,v));
        c.point.setAttribute('cx',x);c.point.setAttribute('cy',Number.isFinite(y)?y:0);
        c.point.style.display=Number.isFinite(y)?'':'none';
      });
      text=`δ ${fmt(values.tpd)} MHz · Gain ${d.gain_digits===3?'≈ ':''}${fmt(ys[0],d.gain_digits||6)} · Noise ${fmt(ys[1])} dB`;
    }
    const prefix=phase==='preview'?'Preview · ':phase==='syncing'?'Updating · ':'';
    status.textContent=prefix+text+(phase!=='ready'?' · Result cards show the last confirmed values.':'');
    root.setAttribute('aria-busy',String(phase==='syncing'));
    root.dataset.phase=phase;
    window.dispatchEvent(new CustomEvent('gabes-plot-status',{
      detail:{scheme:d.scheme,phase,text:status.textContent,owner:root,
        eventKey:d.event_key,values:{...state.values}}
    }));
  }
  function schedule() {
    if(!state.frame)state.frame=requestAnimationFrame(()=>{state.frame=0;paint();});
  }
  function send(payload) {
    if(payload)setStateValue('gesture',JSON.stringify(payload));
    paint();
  }
  function commit(intent) {
    if(state.frame){cancelAnimationFrame(state.frame);state.frame=0;}
    if(!Number.isSafeInteger(intent)) {
      intent=Math.max(Date.now(),(window.__gabesInputIntent||0)+1);
      window.__gabesInputIntent=intent;
    }
    send(releaseGesture(state,intent));
  }
  function cancel() {
    state.previewing=false;
    state.values={...(state.queued||state.inflight?.values||d.original)};
    for(const [key,pair] of Object.entries(inputs))pair.forEach(el=>el.value=state.values[key]);
    if(state.frame){cancelAnimationFrame(state.frame);state.frame=0;}
    paint();
  }
  // Preserve focused numeric inputs and native range drags across replies.
  const schema=JSON.stringify(d.controls.map(c=>[c.name,c.label]));
  if(controls.dataset.schema!==schema) {
    controls.replaceChildren();controls.dataset.schema=schema;
    for(const spec of d.controls) {
      const label=document.createElement('label'),text=document.createElement('span');
      label.dataset.control=spec.name;text.textContent=spec.label;label.append(text);
      const number=document.createElement('input'),range=document.createElement('input');
      number.type='number';range.type='range';label.append(number,range);controls.append(label);
    }
  }
  for(const spec of d.controls) {
    const label=Array.from(controls.children).find(el=>el.dataset.control===spec.name);
    const [number,range]=label.querySelectorAll('input');
    for(const el of [number,range]) {
      el.min=spec.min;el.max=spec.max;el.step=spec.step;
      if(parentElement.activeElement!==el)el.value=state.values[spec.name];
      el.disabled=false;el.setAttribute('aria-label',spec.label);
      el.oninput=()=> {
        const v=el.valueAsNumber;
        if(!Number.isFinite(v)||v<spec.min||v>spec.max)return;
        state.previewing=true;state.values[spec.name]=v;
        (el===number?range:number).value=v;schedule();
      };
      el.onchange=commit;
      el.onkeydown=e=>{if(e.key==='Escape'){e.preventDefault();cancel();el.blur();}};
      el.onpointercancel=cancel;
    }
    inputs[spec.name]=[number,range];
  }
  const inAxis=e=> {
    const r=svg.getBoundingClientRect(),x=(e.clientX-r.left)*d.width/r.width,y=(e.clientY-r.top)*d.height/r.height;
    return d.axes.findIndex(a=>x>=a.left&&x<=a.left+a.width&&y>=a.top&&y<=a.top+a.height);
  };
  function move(e) {
    if(!state.drag || e.pointerId!==state.drag.pointerId)return;
    const r=svg.getBoundingClientRect(),a=state.drag.axis;
    const dx=(e.clientX-state.drag.x)*d.width/r.width/a.width*(a.xlim[1]-a.xlim[0]);
    if(d.kind==='alignment') {
      const spec=d.controls.find(s=>s.name==='shift');
      state.values.shift=Math.min(spec.max,Math.max(spec.min,state.drag.value+dx));
      inputs.shift.forEach(el=>el.value=state.values.shift);
    } else {
      const x=(e.clientX-r.left)*d.width/r.width;
      state.values.tpd=Math.min(a.xlim[1],Math.max(a.xlim[0],a.xlim[0]+(x-a.left)/a.width*(a.xlim[1]-a.xlim[0])));
    }
    state.previewing=true;schedule();
  }
  svg.onpointerdown=e=> {
    if(e.button!==0)return;
    const index=inAxis(e);if(index<0)return;
    svg.setPointerCapture(e.pointerId);
    state.drag={x:e.clientX,value:state.values.shift,axis:d.axes[index],pointerId:e.pointerId};
    state.previewing=true;
    if(d.kind==='cursor')move(e);
  };
  const pointerUp=e=> {
    if(!state.drag || e.pointerId!==state.drag.pointerId)return;
    move(e);state.drag=null;
    if(svg.hasPointerCapture(e.pointerId))svg.releasePointerCapture(e.pointerId);
    commit();
  };
  const pointerCancel=e=> {
    if(!state.drag || e.pointerId!==state.drag.pointerId)return;
    state.drag=null;cancel();
  };
  // A component HTML replacement can remove the SVG holding pointer capture.
  // Window listeners let an already active gesture finish on its new wrapper.
  window.addEventListener('pointermove',move,true);
  window.addEventListener('pointerup',pointerUp,true);
  window.addEventListener('pointercancel',pointerCancel,true);
  const onPreview=e=> {
    if(e.detail?.key!==d.event_key)return;
    if(e.detail.discard) {
      // A direct text/reset edit uses the knob's Python unit parser. Drop any
      // unsent graph request so it cannot overwrite that newer explicit edit.
      state.queued=null;state.inflight=null;state.drag=null;cancel();return;
    }
    if(e.detail.cancel){cancel();return;}
    if(!Number.isFinite(e.detail.value))return;
    state.values.tpd=e.detail.value;
    if(e.detail.commit){e.detail.handled=true;commit(e.detail.intent);}
    else {state.previewing=true;schedule();}
  };
  window.addEventListener('gabes-preview',onPreview);
  const cleanup=()=> {
    window.removeEventListener('gabes-preview',onPreview);
    window.removeEventListener('pointermove',move,true);
    window.removeEventListener('pointerup',pointerUp,true);
    window.removeEventListener('pointercancel',pointerCancel,true);
    if(state.frame){cancelAnimationFrame(state.frame);state.frame=0;}
  };
  state.cleanup=cleanup;
  send(nextGesture(state));
  return ()=> {
    cleanup();
    queueMicrotask(()=> {
      const host=parentElement.host||parentElement;
      if(state.root===root && !host.isConnected && mountedInteractions.get(mountKey)===state)
        mountedInteractions.delete(mountKey);
    });
  };
}
