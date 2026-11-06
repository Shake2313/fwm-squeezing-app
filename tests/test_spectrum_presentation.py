"""Display offsets and logarithmic SVG scales preserve native scientific samples."""
import base64
from pathlib import Path
import shutil
import subprocess

import pytest

from gabes_ui.component_runtime import spectrum_component_js
from gabes_ui.spectral_compare import comparison_plots


DOM = r"""
import assert from 'node:assert/strict';
class Node {
  constructor(tag){
    this.tag=tag;this.children=[];this.attrs={};this.style={};this.dataset={};
    this.clientWidth=600;this.events={};this.offsetWidth=100;this.offsetHeight=24;
  }
  setAttribute(key,value){this.attrs[key]=String(value);}
  append(...nodes){for(const node of nodes){node.parentElement=this;this.children.push(node);}}
  addEventListener(name,callback){this.events[name]=callback;}
  getBoundingClientRect(){return {left:0,top:0,width:600,height:400};}
  querySelectorAll(selector){
    const name=selector.slice(1);
    return descendants(this).filter(node=>(node.attrs.class || node.className || '').split(' ').includes(name));
  }
  querySelector(selector){return this.querySelectorAll(selector)[0] || null;}
  matches(){return true;}
}
const descendants=node=>[node,...node.children.flatMap(descendants)];
globalThis.document={createElement:tag=>new Node(tag),createElementNS:(_,tag)=>new Node(tag)};
const {spectrum,curveSample,spectralRange,plottedCurves,stackedCurves,displaySample,
       previewMarkers,spectrumScale,sampleInsidePlot,conditionMapRows,mapNativeCells,
       conditionMapSample,mapColorLimits,spectrumMap}=await import('MODULE');
const near=(actual,expected)=>assert.ok(Math.abs(actual-expected)<=1e-9*Math.max(1,Math.abs(expected)),`${actual} != ${expected}`);
const drawing=spec=>{
  const host=new Node('div');
  spectrum(host,{title:'Test',description:'Native test samples',unit:'MHz',yLabel:'Value',...spec});
  return {host,svg:host.children[0].children[1]};
};
const finiteSvg=svg=>{
  for(const node of descendants(svg))for(const [name,value] of Object.entries(node.attrs)){
    assert.ok(!/NaN|Infinity/.test(value),`${node.tag}.${name} contains ${value}`);
    if(['x','y','cx','cy','x1','x2','y1','y2','width','height'].includes(name))
      assert.ok(Number.isFinite(Number(value)),`${node.tag}.${name} is not finite: ${value}`);
  }
};
"""


@pytest.fixture(scope='module')
def module_url():
    # Import the exact concatenated runtime module. Its default mount function
    # stays uncalled; exporting private pure/drawing helpers adds no substitute.
    source = spectrum_component_js('spectral_compare.js') + (
        '\nexport {spectrum,curveSample,spectralRange,plottedCurves,previewMarkers,spectrumMap};\n')
    return 'data:text/javascript;base64,' + base64.b64encode(source.encode('utf-8')).decode('ascii')


def _node(module_url, body):
    node=shutil.which('node')
    if not node:
        pytest.skip('Node.js required for shared SVG presentation tests')
    script=DOM.replace('MODULE',module_url)+'\n'+body
    result=subprocess.run([node,'--input-type=module'],input=script,
                          text=True,encoding='utf-8',capture_output=True,timeout=20)
    assert result.returncode == 0, result.stderr[-5000:]


def test_offsets_keep_raw_sampling_hover_values_and_dataset_unchanged(module_url):
    _node(module_url,r"""
const block={series:[
  {id:'a',label:'A',parameter:{value:0,unit:'deg'},x:[0,1,2],observables:{T:[31,32,33]}},
  {id:'b',label:'B',parameter:{value:45,unit:'deg'},x:[.1,1.1,2.1],observables:{T:[41,42,43]}}
]};
const before=JSON.stringify(block),panel={curves:[
  {series_id:'a',observable:'T'},{series_id:'b',observable:'T'}
]};
const curves=plottedCurves(block,panel),stack=stackedCurves(curves,[0,2.1]);
assert.equal(JSON.stringify(block),before);
assert.strictEqual(stack.series[1].x,block.series[1].x);
assert.strictEqual(stack.series[1].y,block.series[1].observables.T);
assert.ok(stack.series[1].display_offset>0);
const raw=curveSample(stack.series[1],[0,1,2],1);
assert.deepEqual(raw,{x:1.1,y:42});
assert.deepEqual(displaySample(stack.series[1],raw),{x:1.1,y:42+stack.series[1].display_offset});
assert.deepEqual(curveSample(stack.series[1],[],1),raw);
const {host,svg}=drawing({grid:[0,1,2],index:1,xlim:[0,2.1],ylim:stack.ylim,
  series:stack.series,yUnit:'%',stacked:{scaleStep:stack.scaleStep,scaleUnit:'percentage points T'}});
const nav=svg._navigation;
svg.events.pointermove({clientX:nav.sx(1),clientY:60,buttons:0,pointerType:'mouse'});
const shownText=descendants(host).filter(node=>node.tag==='b').map(node=>node.textContent);
assert.ok(shownText.some(text=>text.startsWith('42 % @ 1.1 MHz')),
          `Hover must show native 42 rather than its display offset: ${shownText}`);
assert.equal(JSON.stringify(block),before);
finiteSvg(svg);
""")


def test_stacked_limits_cover_all_native_paths_in_the_current_window(module_url):
    _node(module_url,r"""
const curves=[
  {x:[-10,-1,0,1,10],y:[1e6,5,6,7,-1e6],color:'curve-a',label:'A'},
  {x:[-2,-.5,.5,2],y:[-1e6,-8,9,1e6],color:'curve-b',label:'B'},
  {x:[-1,0,1],y:[100,null,110],color:'curve-c',label:'C'},
  {x:[-1,0,1],y:[null,null,null],color:'curve-d',label:'D'}
];
const window=[-1,1],stack=stackedCurves(curves,window);
for(const curve of stack.series)for(let i=0;i<curve.x.length;i++){
  if(curve.x[i]>=window[0] && curve.x[i]<=window[1] && Number.isFinite(curve.y[i])){
    const shown=displaySample(curve,{x:curve.x[i],y:curve.y[i]});
    assert.ok(shown.y>stack.ylim[0] && shown.y<stack.ylim[1]);
  }
}
const {svg}=drawing({grid:[-1,-.5,0,.5,1],index:2,xlim:window,ylim:stack.ylim,
  series:stack.series,stacked:{scaleStep:stack.scaleStep},height:300});
const nav=svg._navigation;
for(const node of descendants(svg).filter(node=>node.tag==='path' && curves.some(c=>c.color===node.attrs.stroke))){
  for(const [,x,y] of node.attrs.d.matchAll(/[ML](-?[\d.]+),(-?[\d.]+)/g)){
    assert.ok(Number(x)>=nav.left-.01 && Number(x)<=nav.right+.01);
    assert.ok(Number(y)>=nav.top-.01 && Number(y)<=nav.bottom+.01,
              `Visible path escaped stacked limits: ${node.attrs.d}`);
  }
}
finiteSvg(svg);
""")


def test_stacked_export_arrays_are_the_unshifted_native_values():
    dataset=dict(axis=dict(label='Field',unit='uT',reference='Command zero'),
                 series=[dict(id='a',label='A',x=[0.,1.,2.],observables={'T':[31.,None,33.]}),
                         dict(id='b',label='B',x=[.1,1.1,2.1],observables={'T':[41.,42.,43.]})],
                 panels=[dict(title='Separated transmission',y_label='Transmission',unit='%',
                              stackable=True,curves=[dict(series_id='a',observable='T'),
                                                     dict(series_id='b',observable='T')])])
    traces=comparison_plots(dataset)[0]['traces']
    for trace,condition in zip(traces,dataset['series']):
        assert trace['x_data'] is condition['x']
        assert trace['y_data'] is condition['observables']['T']
    assert traces[1]['y_data'] == [41.,42.,43.]
    assert traces[0]['y_data'][1] is None


def test_log_y_uses_positive_limits_equal_decades_and_gaps_for_nonpositive_samples(module_url):
    _node(module_url,r"""
const grid=[0,1,2,3,4,5,6,7],y=[1,-1,0,null,10,100,NaN,1000];
const curve={x:grid,y,color:'curve-positive',label:'ASD'};
const range=spectralRange([curve],[0,7],undefined,true);
assert.ok(range[0]>0 && range[0]<1 && range[1]>1000 && range.every(Number.isFinite));
const scale=spectrumScale(range[0],range[1],145,23,true);
for(const [a,b,c] of [[1,10,100],[10,100,1000]])near(scale.forward(a)-scale.forward(b),scale.forward(b)-scale.forward(c));
for(const value of [1,2.5,10,100,1000])near(scale.inverse(scale.forward(value)),value);
const {svg}=drawing({grid,index:0,xlim:[0,7],ylim:range,yScale:'log',series:[curve]});
const path=descendants(svg).find(node=>node.tag==='path' && node.attrs.stroke===curve.color);
assert.equal((path.attrs.d.match(/M/g)||[]).length,3);
assert.equal((path.attrs.d.match(/[ML]/g)||[]).length,4);
const nav=svg._navigation;
near(nav.sy(1)-nav.sy(10),nav.sy(10)-nav.sy(100));
assert.ok(nav.sy(1)>nav.sy(10) && nav.sy(10)>nav.sy(100));
finiteSvg(svg);
for(const missing of [[0,-1,null,NaN],[]]){
  const fallback=spectralRange([{x:missing.map((_,i)=>i),y:missing}],[0,4],undefined,true);
  assert.ok(fallback[0]>0 && fallback[0]<fallback[1] && fallback.every(Number.isFinite));
}
""")


@pytest.mark.parametrize('invalid', ['0','-1','null','NaN'])
def test_log_y_cursor_hover_and_drag_keep_finite_svg_coordinates_for_gaps(module_url,invalid):
    _node(module_url,r"""
const grid=[0,1,2],curve={x:grid,y:[1,INVALID,100],color:'curve-gap',label:'ASD'};
const {host,svg}=drawing({grid,index:1,xlim:[0,2],ylim:[.5,200],yScale:'log',series:[curve]});
finiteSvg(svg);
const dot=svg.querySelector('.mc-nav-dot');
assert.equal(dot.attrs.visibility,'hidden');
const nav=svg._navigation;
previewMarkers(svg,[0,2],1);
finiteSvg(svg);
svg.events.pointermove({clientX:nav.sx(1),clientY:60,buttons:0,pointerType:'mouse'});
finiteSvg(svg);
assert.equal(sampleInsidePlot({x:1,y:curve.y[1]},nav),false);
""".replace('INVALID',invalid))


def test_stacked_zoom_labels_and_baselines_use_visible_native_samples(module_url):
    _node(module_url,r"""
const curves=[
  {x:[-10,-1,0,1,10],y:[1e6,5,6,7,-1e6],color:'curve-a',label:'A'},
  {x:[-10,-.5,.5,10],y:[-1e6,8,9,1e6],color:'curve-b',label:'B'},
  {x:[-10,10],y:[-1e6,1e6],color:'curve-empty',label:'No sample in zoom'}
];
const window=[-1,1],stack=stackedCurves(curves,window);
const {svg}=drawing({grid:[-1,-.5,0,.5,1],index:2,xlim:window,ylim:stack.ylim,
  series:stack.series,stacked:{scaleStep:stack.scaleStep},height:300});
const nav=svg._navigation;
for(const curve of stack.series){
  const index=curve.x.findIndex((x,i)=>x>=window[0] && x<=window[1] && Number.isFinite(curve.y[i]));
  const label=descendants(svg).find(node=>node.tag==='text' && node.textContent===curve.label);
  if(index<0){assert.equal(label,undefined);continue;}
  assert.ok(label,`Visible curve ${curve.label} needs its angle label`);
  const expected=nav.sy(curve.y[index]+curve.display_offset);
  near(Number(label.attrs.y),expected+3);
  assert.ok(Number(label.attrs.y)>=nav.top && Number(label.attrs.y)<=nav.bottom+3,
            `${curve.label} label escaped zoom: ${label.attrs.y}`);
}
const baselines=descendants(svg).filter(node=>node.tag==='line' && Number(node.attrs['stroke-width'])===.7);
assert.equal(baselines.length,2);
for(const line of baselines){
  const y=Number(line.attrs.y1);
  assert.ok(y>=nav.top && y<=nav.bottom,`Baseline escaped zoom: ${y}`);
  assert.equal(line.attrs.y1,line.attrs.y2);
}
finiteSvg(svg);
""")
