"""Actual bundled map sampler and drawing retain native coordinates and gaps."""
from test_spectrum_presentation import module_url, _node


def test_native_nonuniform_coordinates_and_blank_cells(module_url):
    _node(module_url,r"""
const block={varied:{key:'power',label:'Power',unit:'mW'},series:[
  {id:'a',label:'A',parameter:{value:1},x:[-3,0,3],observables:{T:[.2,.8,.2]}},
  {id:'b',label:'B',parameter:{value:4},x:[-2,.25,2],observables:{T:[.1,null,.1]}}
],sweep:{requested_values:[1,2,4]}};
const panel={curves:[{series_id:'a',observable:'T'},{series_id:'b',observable:'T'}]};
const before=JSON.stringify(block),data=conditionMapRows(block,panel);
assert.deepEqual(data.rows.map(row=>row.value),[1,2,4]);assert.equal(data.numeric,true);
const native=conditionMapSample(data.rows,.2,4);
assert.equal(native.point.x,.25);assert.equal(native.point.y,null);
assert.equal(conditionMapSample(data.rows,2.5,4).point.x,null);
assert.equal(conditionMapSample(data.rows,0,2).scope.status,'not calculated');
assert.equal(conditionMapSample(data.rows,.1,1).point.y,.8);
assert.deepEqual(mapNativeCells(data.rows[2].curve).map(cell=>[cell.low,cell.high]),[[-2,-.875],[-.875,1.125],[1.125,2]]);
assert.equal(JSON.stringify(block),before);
assert.deepEqual(mapColorLimits(data.rows),[.1,.8]);
const categorical=conditionMapRows({series:block.series.map((row,i)=>({...row,parameter:{value:['D1','D2'][i]}}))},panel);
assert.equal(categorical.numeric,false);assert.deepEqual(categorical.rows.map(row=>row.coordinate),[0,1]);
""")


def test_canvas_and_svg_use_native_cells_shared_scale_and_keyboard_samples(module_url):
    _node(module_url,r"""
globalThis.window={devicePixelRatio:1};
const painted=[];
Node.prototype.getContext=function(){return {scale(){},clearRect(){},fillRect(...args){painted.push({fill:this.fillStyle,args});}};};
const dataset={axis:{label:'Detuning',unit:'MHz'},varied:{label:'Power',unit:'mW'},series:[
  {id:'a',label:'A',parameter:{value:1},x:[-3,0,3],observables:{T:[.1,.9,.1]}},
  {id:'b',label:'B',parameter:{value:4},x:[-2,.25,2],observables:{T:[.2,null,.2]}}
]};
const panel={title:'Transmission',y_label:'T',unit:'1',curves:[{series_id:'a',observable:'T'},{series_id:'b',observable:'T'}]};
const selection={series_id:'a',x:0},chosen=[],host=new Node('div');
const svg=spectrumMap(host,{dataset,panel,window:[-3,3],selection,onKeyboardSelect:(...args)=>chosen.push(args)});
finiteSvg(svg);const n=svg._navigation;
assert.equal(painted.length,6); // background + three A cells + two finite B cells
const colored=painted.slice(1);assert.ok(colored.every(cell=>cell.args.every(Number.isFinite)));
assert.ok(colored.every(cell=>cell.args.every(Number.isInteger))); // no fractional raster seams
assert.equal(colored[0].args[0]+colored[0].args[2],colored[1].args[0]);
assert.equal(colored[1].args[0]+colored[1].args[2],colored[2].args[0]);
const bCells=colored.slice(3);
assert.ok(bCells.every(cell=>cell.args[0]>n.left)); // B cannot extend to A's -3
svg.events.keydown({key:'ArrowUp',preventDefault(){}});
assert.deepEqual(chosen.at(-1),['b',.25]);
svg.events.keydown({key:'ArrowRight',preventDefault(){}});
assert.deepEqual(chosen.at(-1),['a',3]);
n.mark({series_id:'b',x:.25});finiteSvg(svg);
const rect=svg.querySelector('.sc-map-selection');assert.equal(rect.attrs.visibility,'visible');
assert.deepEqual(mapColorLimits([{curve:{y:[0,-1,null,1,100]}}],true),[1,100]);
const large={curve:{y:Array.from({length:200001},(_,i)=>i)}};
assert.deepEqual(mapColorLimits([large]),[0,200000]);
""")


def test_alternate_rabi_axis_keeps_native_parameter_and_pending_identity(module_url):
    _node(module_url,r"""
const dataset={varied:{label:'Coupling power',unit:'mW'},condition_axes:[
  {key:'parameter',label:'Coupling power',unit:'mW'},
  {key:'coupling_rabi',label:'Effective Rabi',unit:'MHz',transform:{factor:3,exponent:.5}}
],series:[{id:'p1',label:'P1',parameter:{value:1},coordinates:{coupling_rabi:3},x:[-1,0,1],observables:{T:[.2,.8,.2]}},
{id:'p4',label:'P4',parameter:{value:4},coordinates:{coupling_rabi:6},x:[-2,0,2],observables:{T:[.2,.9,.2]}}],sweep:{requested_values:[1,2,4]}};
const before=JSON.stringify(dataset),panel={curves:[{series_id:'p1',observable:'T'},{series_id:'p4',observable:'T'}]};
const data=conditionMapRows(dataset,panel,'coupling_rabi');
assert.deepEqual(data.rows.map(row=>row.coordinate),[3,3*Math.sqrt(2),6]);
assert.equal(data.rows[2].value,4);assert.equal(data.rows[2].id,'p4');
assert.equal(conditionMapSample(data.rows,0,6).row.value,4);
assert.equal(data.axis.unit,'MHz');assert.equal(JSON.stringify(dataset),before);
""")


def test_selected_condition_navigation_uses_its_grid_not_the_union(module_url):
    _node(module_url,r"""
const {conditionNavigationGrid}=await import('MODULE');
const data={series:[{id:'a',x:[-3,0,3]},{id:'b',x:[-1,.1,1]}]},union=[-3,-1,0,.1,1,3];
assert.deepEqual(conditionNavigationGrid(data,'a',union),[-3,0,3]);
assert.deepEqual(conditionNavigationGrid(data,'b',union),[-1,.1,1]);
assert.strictEqual(conditionNavigationGrid(data,'missing',union),union);
""".replace('MODULE',module_url))


def test_sparse_if_map_paints_only_native_points_and_keeps_unsolved_space_blank(module_url):
    _node(module_url,r"""
globalThis.window={devicePixelRatio:1};const painted=[];
Node.prototype.getContext=function(){return {scale(){},clearRect(){},fillRect(...args){painted.push(args);}};};
const block={axis:{label:'Optical frequency',unit:'MHz'},varied:{label:'RF',unit:'MHz'},series:[
 {id:'off',label:'RF off',parameter:{value:0},x:[-4,0,4],observable_x:{asd:[]},observables:{asd:[]}},
 {id:'on',label:'RF on',parameter:{value:4},x:[-4,0,4],observable_x:{asd:[-1,1]},observables:{asd:[12,14]}}
]};
const panel={key:'asd',title:'ASD',y_label:'Field ASD',unit:'nV/cm/√Hz',scale:'log',sampled_only:true,
 curves:[{series_id:'off',observable:'asd'},{series_id:'on',observable:'asd'}]};
const before=JSON.stringify(block),rows=conditionMapRows(block,panel).rows;
assert.deepEqual(rows[1].curve.x,[-1,1]);
assert.equal(conditionMapSample(rows,0,4,.01).scope.status,'not evaluated');
assert.deepEqual(conditionMapSample(rows,1.001,4,.01).point,{x:1,y:14});
const svg=spectrumMap(new Node('div'),{dataset:block,panel,window:[-4,4],selection:{series_id:'on',x:1},onKeyboardSelect(){}});
assert.equal(painted.length,1);
assert.equal(svg.querySelectorAll('.sc-map-native-point').length,2);finiteSvg(svg);
assert.equal(JSON.stringify(block),before);
""")


def test_optional_curves_and_each_observable_grid_remain_native(module_url):
    _node(module_url,r"""
const {panelCurves}=await import('MODULE');
const block={series:[{id:'nominal',label:'Nominal',x:[-4,0,4],observable_x:{response:[-.2,.7]},
observables:{T:[.2,.8,.2],response:[2,3]},overlays:{residual:{panel:'T',label:'Residual',x:[-4,0,4],y:[.3,.7,.3]}}},
{id:'scenario',label:'Assumed scenario',x:[-4,0,4],observables:{T:[.1,.9,.1]}}]};
const before=JSON.stringify(block);
const panel={key:'T',curves:[{series_id:'nominal',observable:'T'},
{series_id:'scenario',observable:'T',optional_overlays:['opacity']}]};
assert.equal(panelCurves(block,panel).length,1);
const curves=panelCurves(block,panel,['opacity','residual'],'nominal');assert.equal(curves.length,3);
assert.strictEqual(curves[0].y,block.series[0].observables.T);assert.equal(curves[2].label,'Residual');
const response=panelCurves(block,{curves:[{series_id:'nominal',observable:'response'}]});
assert.deepEqual(response[0].x,[-.2,.7]);assert.equal(response[0].color,curves[0].color);
assert.equal(JSON.stringify(block),before);
""".replace('MODULE',module_url))


def test_wrapped_phase_and_sparse_paths_never_draw_false_connecting_lines(module_url):
    _node(module_url,r"""
const x=[0,1,2],phase=[179,-179,-170],before=[...phase];
const {svg}=drawing({grid:x,index:1,xlim:[0,2],ylim:[-180,180],cursor:false,
 series:[{x,y:phase,color:'phase',label:'Wrapped phase',cyclic_period:360}]});
const path=descendants(svg).find(n=>n.tag==='path' && n.attrs.stroke==='phase');
assert.equal((path.attrs.d.match(/M/g)||[]).length,2); // break at the ±180° wrap
assert.equal((path.attrs.d.match(/L/g)||[]).length,1);assert.deepEqual(phase,before);
const sparse=drawing({grid:x,index:1,xlim:[0,2],ylim:[0,4],cursor:false,
 series:[{x,y:[1,2,3],color:'sparse',label:'Sparse IF samples',points_only:true}]});
const points=descendants(sparse.svg).find(n=>n.tag==='path' && n.attrs.stroke==='sparse');
assert.ok(!points.attrs.d.includes('L'));finiteSvg(sparse.svg);
""")


def test_selected_and_optimum_markers_are_distinct_without_changing_samples(module_url):
    _node(module_url,r"""
const x=[1],y=[2],before=JSON.stringify([x,y]);
const selected=drawing({grid:x,index:0,xlim:[0,2],ylim:[0,3],cursor:false,
series:[{x,y,color:'native',label:'Selected',markers:true,point_kind:'selected'}]});
assert.ok(descendants(selected.svg).some(n=>n.tag==='circle' && n.attrs.r==='5' && n.attrs.fill==='white'));
const optimum=drawing({grid:x,index:0,xlim:[0,2],ylim:[0,3],cursor:false,
series:[{x,y,color:'native',label:'Sampled optimum',markers:true,point_kind:'optimum'}]});
assert.ok(descendants(optimum.svg).some(n=>n.tag==='path' && n.attrs.fill==='white' && n.attrs.d.endsWith('Z')));
finiteSvg(optimum.svg);assert.equal(JSON.stringify([x,y]),before);
""")


def test_empty_map_has_no_invented_numeric_color_range(module_url):
    _node(module_url,r"""
globalThis.window={devicePixelRatio:1};
Node.prototype.getContext=function(){return {scale(){},clearRect(){},fillRect(){}};};
const data={axis:{label:'Optical frequency',unit:'MHz'},varied:{label:'IF',unit:'kHz'},series:[
{id:'empty',label:'RF off',parameter:{value:40},x:[-1,0,1],observable_x:{asd:[]},observables:{asd:[]},unavailable_observables:{asd:'RF off'}}]};
const panel={title:'ASD',y_label:'Field ASD',unit:'nV/cm/√Hz',scale:'log',curves:[{series_id:'empty',observable:'asd'}]},host=new Node('div');
spectrumMap(host,{dataset:data,panel,window:[-1,1],selection:{series_id:'empty',x:0},onKeyboardSelect(){}});
const legend=host.querySelector('.sc-map-key');
assert.ok(legend.children.some(node=>node.textContent==='No finite samples for this quantity.'));
assert.ok(!legend.children.some(node=>node.tag==='i'));
assert.equal(conditionMapSample(conditionMapRows(data,panel).rows,0,40).scope.status,'RF off');
""")
