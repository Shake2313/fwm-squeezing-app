"""Logarithmic zooms retain visible, correctly positioned frequency ticks."""
import base64
from pathlib import Path
import shutil
import subprocess

import pytest


def test_short_log_windows_render_frequency_ticks():
    node = shutil.which("node")
    if not node:
        pytest.skip("Node.js required for spectrum tick tests")
    source = (Path(__file__).resolve().parents[1] /
              "assets/ui/spectrum_plot.js").read_text(encoding="utf-8")
    module = "data:text/javascript;base64," + base64.b64encode(
        (source + "\nexport {spectrum};\n").encode()).decode()
    script = r"""
import assert from 'node:assert/strict';
class Node {
  constructor(tag){this.tag=tag;this.children=[];this.attrs={};this.style={};this.clientWidth=600;}
  setAttribute(key,value){this.attrs[key]=String(value);}
  append(...nodes){this.children.push(...nodes);}
  addEventListener(){}
}
globalThis.document={createElement:tag=>new Node(tag),createElementNS:(_,tag)=>new Node(tag)};
const {spectrum}=await import('MODULE');
const descendants=node=>[node,...node.children.flatMap(descendants)];
for(const [lo,hi] of [[2,5],[2.1,2.2],[2,50],[.001,1000]]) {
  const host=new Node('div'),grid=[lo,Math.sqrt(lo*hi),hi];
  spectrum(host,{grid,index:1,xlim:[lo,hi],ylim:[0,1],unit:'MHz',xScale:'log',
    series:[{x:grid,y:[.2,.4,.6],color:'green',label:'RF'}],title:'RF',yLabel:'T'});
  const svg=host.children[0].children[1],n=svg._navigation;
  const ticks=descendants(svg).filter(el=>el.tag==='text'&&Number(el.attrs.y)===n.bottom+15);
  assert.ok(ticks.length>=2,`Window ${lo} to ${hi} needs frequency tick labels`);
  for(const tick of ticks){
    const value=Number(tick.textContent),pixel=Number(tick.attrs.x);
    assert.ok(Number.isFinite(value)&&value>=lo&&value<=hi);
    assert.ok(pixel>=n.left-1e-8&&pixel<=n.right+1e-8);
    assert.ok(Math.abs(pixel-n.sx(value))<1e-8);
  }
}
""".replace("MODULE", module)
    result = subprocess.run([node, "--input-type=module"], input=script,
                            text=True, encoding="utf-8", capture_output=True,
                            timeout=20)
    assert result.returncode == 0, result.stderr[-4000:]
