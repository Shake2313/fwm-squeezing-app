"""Pan, cursor confinement and one commit per gesture are shared chart contracts."""
import base64
from pathlib import Path
import shutil
import subprocess

import pytest


def test_pointer_preview_commit_cancel_and_window_geometry():
    node=shutil.which("node")
    if not node:
        pytest.skip("Node.js required for spectrum gesture tests")
    src=Path(__file__).resolve().parents[1] / "assets/ui/spectrum_gestures.js"
    url="data:text/javascript;base64,"+base64.b64encode(src.read_bytes()).decode()
    script=r"""
import assert from 'node:assert/strict';
import {bindSpectrumGesture,navigateSpectrumWindow} from 'MODULE';
const domain=[0,100], initial={window:[20,50],cursor:44,x:25};
const adjust=(action,x)=>navigateSpectrumWindow(action,initial,x,domain,2);
assert.deepEqual(adjust('pan',100),{window:[70,100],cursor:94});
assert.deepEqual(adjust('pan',-50),{window:[0,30],cursor:24});
assert.deepEqual(adjust('cursor',75),{window:[20,50],cursor:50});
assert.deepEqual(adjust('center',99),{window:[70,100],cursor:99});
assert.deepEqual(adjust('left',49),{window:[48,50],cursor:48});
assert.deepEqual(adjust('right',18),{window:[20,22],cursor:22});
const events=[], frames=new Map();let seq=0;
globalThis.window=new EventTarget();
globalThis.requestAnimationFrame=cb=>{frames.set(++seq,cb);return seq;};
globalThis.cancelAnimationFrame=id=>frames.delete(id);
const flush=()=>{const batch=[...frames.values()];frames.clear();for(const cb of batch)cb();};
const surface={},state={drag:null};
const bind=()=>bindSpectrumGesture(surface,state,{
  begin:()=>({}),preview:(_,e)=>events.push('p'+e.clientX),
  commit:(_,e)=>events.push('c'+e.clientX),cancel:()=>events.push('cancel')
});
function pointer(type,x,id=3){const e=new Event(type,{cancelable:true});
  Object.assign(e,{pointerId:id,clientX:x,button:0});return e;}
let cleanup=bind();
const down=pointer('pointerdown',2);surface.onpointerdown(down);
assert.equal(down.defaultPrevented,true); // dragging a plot must not select its tick text
window.dispatchEvent(pointer('pointermove',21));window.dispatchEvent(pointer('pointermove',22));
window.dispatchEvent(pointer('pointerup',25,7));flush(); // second finger has no effect
assert.deepEqual(events,['p2','p22']);
cleanup();cleanup=bind(); // replacement component owns the same gesture state
window.dispatchEvent(pointer('pointermove',30));window.dispatchEvent(pointer('pointerup',41));flush();
assert.deepEqual(events,['p2','p22','p41','c41']);
surface.onpointerdown(pointer('pointerdown',2));window.dispatchEvent(pointer('pointermove',90));
const esc=new Event('keydown');Object.assign(esc,{key:'Escape'});window.dispatchEvent(esc);
window.dispatchEvent(pointer('pointerup',95));flush();cleanup();
assert.deepEqual(events.slice(-2),['p2','cancel']);assert.equal(state.drag,null);
""".replace("MODULE",url)
    subprocess.run([node,"--input-type=module"],input=script,text=True,
                   capture_output=True,check=True,timeout=20)
