"""Dropdown keyboard policy blocks editing while preserving navigation and number entry."""
import base64
import json
from pathlib import Path
import shutil
import subprocess

import pytest


def test_native_dropdown_policy_keeps_navigation_and_numeric_text_edits():
    node = shutil.which("node")
    if not node:
        pytest.skip("Node.js is required for browser keyboard policy")
    source = Path(__file__).resolve().parents[1] / "assets/ui/dropdown_policy.js"
    url = "data:text/javascript;base64," + base64.b64encode(source.read_bytes()).decode()
    script = f"""
import mount from {json.dumps(url)};
const listeners=new Map(), attrs={{}};
const dropdown={{readOnly:false,setAttribute:(k,v)=>attrs[k]=v,matches:()=>true}};
const numeric={{matches:()=>false}};
const marker={{}};
globalThis.document={{body:{{}},querySelectorAll:()=>[dropdown],
  addEventListener:(name,fn)=>listeners.set(name,fn),
  removeEventListener:name=>listeners.delete(name)}};
globalThis.MutationObserver=class {{observe(){{}}disconnect(){{}}}};
const cleanup=mount({{parentElement:{{querySelector:()=>marker}}}});
const press=(key,target=dropdown,ctrlKey=false)=>{{
  const event={{key,target,ctrlKey,preventDefault(){{this.blocked=true}},stopPropagation(){{}}}};
  listeners.get('keydown')(event);return !!event.blocked;
}};
const edited=['a','4','Backspace','Delete'].map(key=>press(key));
const navigation=['ArrowDown','ArrowUp','Home','End','Enter','Escape','Tab',' '].map(key=>press(key));
const numericEdit=press('4',numeric);
const paste=press('v',dropdown,true);
cleanup();
console.log(JSON.stringify({{readonly:dropdown.readOnly,attrs,edited,navigation,numericEdit,paste,remaining:listeners.size}}));
"""
    result = subprocess.run([node, "--input-type=module"], input=script, text=True,
                            encoding="utf-8", capture_output=True, check=True, timeout=15)
    state = json.loads(result.stdout)
    assert state["readonly"] and state["attrs"]["inputmode"] == "none"
    assert state["attrs"]["aria-autocomplete"] == "none"
    assert all(state["edited"]) and not any(state["navigation"])
    assert not state["numericEdit"] and state["paste"]
    assert state["remaining"] == 0
