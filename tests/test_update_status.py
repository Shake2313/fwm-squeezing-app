"""Status never claims synchronization while a browser gesture is pending."""
import base64
import json
from pathlib import Path
import shutil
import subprocess

import pytest


def test_browser_update_states_cover_pending_preview_and_failure():
    node = shutil.which('node')
    if not node:
        pytest.skip('Node.js required for browser status logic')
    source = Path(__file__).resolve().parents[1] / 'assets/ui/update_status.js'
    url = 'data:text/javascript;base64,' + base64.b64encode(source.read_bytes()).decode()
    script = f'''import {{updateMessage}} from {json.dumps(url)};
console.log(JSON.stringify([
updateMessage(false,'syncing'),updateMessage(true,'ready'),
updateMessage(false,'preview'),updateMessage(false,'ready'),
updateMessage(false,'ready',true)]));'''
    result = subprocess.run([node, '--input-type=module'], input=script, text=True, encoding='utf-8',
                            capture_output=True, check=True, timeout=15)
    messages = json.loads(result.stdout)
    assert [m['phase'] for m in messages] == ['syncing', 'syncing', 'preview', 'ready', 'error']
    assert 'committed' in messages[2]['text']
    assert messages[3]['text'] == ''
    assert 'failed' in messages[-1]['text']


def test_removed_plot_does_not_leave_update_status_stuck():
    node = shutil.which('node')
    if not node:
        pytest.skip('Node.js required for browser status lifecycle')
    source = Path(__file__).resolve().parents[1] / 'assets/ui/update_status.js'
    url = 'data:text/javascript;base64,' + base64.b64encode(source.read_bytes()).decode()
    script = f'''import mount from {json.dumps(url)};
const listeners=new Map();let running=false,observe;
globalThis.window={{addEventListener:(n,f)=>listeners.set(n,f),removeEventListener:n=>listeners.delete(n)}};
globalThis.MutationObserver=class {{constructor(f){{observe=f;}}observe(){{}}disconnect(){{}}}};
const root={{dataset:{{}},textContent:''}},app={{getAttribute:()=>running?'running':'notRunning'}};
globalThis.document={{querySelector:s=>s.includes('stException')?null:app}};
const cleanup=mount({{parentElement:{{querySelector:()=>root}},data:{{scheme:'fwm'}}}});
const owner={{isConnected:true}},phases=[];
listeners.get('gabes-plot-status')({{detail:{{scheme:'fwm',phase:'syncing',owner}}}});
running=true;observe();running=false;observe();phases.push(root.dataset.phase);
running=true;observe();owner.isConnected=false;running=false;observe();phases.push(root.dataset.phase);
cleanup();console.log(JSON.stringify({{phases,text:root.textContent,listeners:listeners.size}}));'''
    result = subprocess.run([node, '--input-type=module'], input=script, text=True,
                            encoding='utf-8', capture_output=True, check=True, timeout=15)
    assert json.loads(result.stdout) == {'phases': ['syncing', 'ready'], 'text': '', 'listeners': 0}
