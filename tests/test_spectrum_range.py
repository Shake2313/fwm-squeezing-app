"""Direct SAS window edits reject invalid limits before mutating navigation."""
import base64
from pathlib import Path
import re
import shutil
import subprocess

import pytest


def test_sas_window_input_validation_and_recovery():
    node = shutil.which("node")
    if not node:
        pytest.skip("Node.js required for spectrum range tests")
    source = (Path(__file__).resolve().parents[1] /
              "assets/ui/model_compare.js").read_text(encoding="utf-8")
    body = re.search(r"function applyRange\(\)\{(.*?)\n  \}", source, re.S).group(1)
    encoded = base64.b64encode(body.encode()).decode()
    script = r"""
import assert from 'node:assert/strict';
const body=Buffer.from('BODY','base64').toString('utf8');
const from={valueAsNumber:20},to={valueAsNumber:80,validity:'',reports:0,
  setCustomValidity(message){this.validity=message;},
  reportValidity(){this.reports++;}};
const state={sample_index:1},d={x:[0,50,100]},full=[0,100];
let applied=null,refreshes=0,published=0;
const edit=new Function('from','to','minWidth','format','applyWindow','d','state',
  'refresh','publish','full',body);
const apply=()=>edit(from,to,2,String,value=>{applied=value;},d,state,
  ()=>refreshes++,()=>published++,full);
for(const [a,b] of [[NaN,80],[20,NaN],[Infinity,80],[20,-Infinity],
                    [-1,80],[20,101],[101,110],[-10,-1],[80,20],[20,21]]) {
  from.valueAsNumber=a;to.valueAsNumber=b;
  apply();
  assert.equal(applied,null);
  assert.equal(refreshes,0);
  assert.equal(published,0);
  assert.ok(to.validity.length>0,`${a}, ${b} must explain its rejection`);
}
assert.equal(to.reports,10);
// A valid edit after an invalid one clears the error and commits exactly once.
from.valueAsNumber=0;to.valueAsNumber=100;
apply();
assert.equal(to.validity,'');
assert.deepEqual(applied,{window:[0,100],cursor:50});
assert.equal(refreshes,1);
assert.equal(published,1);
assert.equal(to.reports,10);
""".replace("BODY", encoded)
    result = subprocess.run([node, "--input-type=module"], input=script,
                            text=True, encoding="utf-8", capture_output=True,
                            timeout=20)
    assert result.returncode == 0, result.stderr[-4000:]
