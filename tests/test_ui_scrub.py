"""ScrubField gesture handling (gabes_ui.scrub)."""
import pytest

from gabes_ui.scrub import decode_gesture, resolve_gesture, state_key

KNOB = dict(unit="mW", vmin=0.0, vmax=2.0, integer=False, default=0.5)


def test_component_keys_drop_the_reserved_double_underscore():
    # Streamlit refuses "__" inside a bidirectional component's id.
    assert state_key("sas__pump_power_mw") == "scrub-sas-pump_power_mw"
    assert "__" not in state_key("fwm__excess_noise_slope")


def test_a_gesture_arrives_as_json_text():
    assert decode_gesture('{"n": 3, "v": 1.25}') == {"n": 3, "v": 1.25}
    assert decode_gesture({"n": 1, "r": True}) == {"n": 1, "r": True}
    assert decode_gesture(None) is None
    assert decode_gesture("[]") is None


def test_a_dragged_value_is_taken_as_it_comes():
    assert resolve_gesture({"n": 1, "v": 1.25}, **KNOB) == (1.25, None)


def test_typed_text_is_parsed_in_the_knobs_unit():
    assert resolve_gesture({"n": 1, "t": "1.25"}, **KNOB) == (1.25, None)
    assert resolve_gesture({"n": 1, "t": "1200 uW"}, **KNOB) == (1.2, None)


def test_out_of_range_entry_lands_on_the_boundary_and_says_so():
    value, note = resolve_gesture({"n": 1, "t": "1.25 W"}, **KNOB)
    assert value == 2.0
    assert note == "clamped to 2"


def test_an_entry_in_the_wrong_quantity_leaves_the_knob_alone():
    value, note = resolve_gesture({"n": 1, "t": "3 MHz"}, **KNOB)
    assert value is None
    assert note == "expected mW"


def test_reset_restores_the_default():
    assert resolve_gesture({"n": 1, "r": True}, **KNOB) == (0.5, None)


def test_integer_knobs_stay_integers():
    knob = dict(unit="", vmin=201, vmax=4001, integer=True, default=1401)
    value, _ = resolve_gesture({"n": 1, "v": 1400.6}, **knob)
    assert value == 1401 and isinstance(value, int)


def test_nothing_to_do_is_not_an_error():
    assert resolve_gesture(None, **KNOB) == (None, None)
    assert resolve_gesture({"n": 1}, **KNOB) == (None, None)


@pytest.mark.parametrize("text,expected", [("2", 2.0), ("0", 0.0), ("-1", 0.0)])
def test_the_range_is_always_respected(text, expected):
    value, _ = resolve_gesture({"n": 1, "t": text}, **KNOB)
    assert value == expected


def test_late_scrub_callback_cannot_overwrite_newer_input(monkeypatch):
    import json
    from gabes_ui import scrub
    key = 'sas__pump_power_mw'
    component = state_key(key)
    state = {key: 0.5, component: {'value': json.dumps({'n': 2, 'v': 1.4})}}
    monkeypatch.setattr(scrub.st, 'session_state', state)
    commit = lambda: scrub._commit(key, component, 'mW', 0., 2., False, 0.5, None)
    commit()
    state[component]['value'] = json.dumps({'n': 1, 'v': 0.8})
    commit()
    assert state[key] == 1.4
    assert state[component + '-ack'] == 2


def test_invalid_typed_value_is_acknowledged_without_changing_physics(monkeypatch):
    import json
    from gabes_ui import scrub
    key = 'sas__pump_power_mw'
    component = state_key(key)
    state = {key: 0.5, component: {'value': json.dumps({'n': 3, 't': 'bad'})}}
    monkeypatch.setattr(scrub.st, 'session_state', state)
    scrub._commit(key, component, 'mW', 0., 2., False, 0.5, None)
    assert state[key] == 0.5
    assert state[component + '-ack'] == 3
    assert state[key + scrub.NOTE_SUFFIX]


def test_browser_draft_pointer_and_keyboard_survive_replaced_html():
    """Exercise the real handlers with replacement DOM and an older reply."""
    import base64
    import json
    import shutil
    import subprocess
    from gabes_ui import scrub
    node = shutil.which("node")
    if not node:
        pytest.skip("Node.js is needed to exercise the browser input lifecycle")
    url = "data:text/javascript;base64," + base64.b64encode(scrub._JS.encode()).decode()
    script = r'''
import assert from 'node:assert/strict';
import render, {scrubStateForComponent} from MODULE_URL;
globalThis.window = new EventTarget();
globalThis.CustomEvent = class extends Event {
  constructor(name, options={}) { super(name); this.detail=options.detail; }
};
const sent = [];
const observers = new Set();
const app = {getAttribute(){return 'notRunning'}};
globalThis.document = {querySelector(){return app}};
globalThis.MutationObserver = class {
  constructor(callback){this.callback=callback}
  observe(){observers.add(this)}
  disconnect(){observers.delete(this)}
};
function element() {
  return { hidden:false, textContent:'', value:'', style:{}, selectionStart:0, selectionEnd:0,
    classList:{toggle(){},add(){},remove(){}}, setAttribute(){},removeAttribute(){},
    setPointerCapture(){},hasPointerCapture(){return false},releasePointerCapture(){},
    getBoundingClientRect(){return {left:0,width:100}},
    focus(){this.focused=true},select(){},
    setSelectionRange(a,b){this.selectionStart=a;this.selectionEnd=b} };
}
const defaults = {value:0.5,text:'0.5',min:0,max:2,step:0.1,decimals:1,default:0.5,
  unit:'mW',label:'Power',ack:0};
function mount(key, data={}) {
  const children = new Map();
  const root = {...element(),isConnected:true,querySelector(name){
    if(!children.has(name))children.set(name,element());return children.get(name);
  }};
  const dispose = render({key,data:{...defaults,...data},
    parentElement:{querySelector(){return root}},
    setStateValue(name,payload){sent.push({key,root,...JSON.parse(payload)});} });
  return {root,children,dispose};
}
const event = (values={}) => ({button:0,pointerId:7,clientX:50,preventDefault(){},...values});
function pointer(type,clientX) {
  const e = new Event(type);Object.assign(e,{pointerId:7,clientX});window.dispatchEvent(e);
}
const wait = ms => new Promise(resolve=>setTimeout(resolve,ms));
// An uncommitted unit-bearing draft and caret survive a different input's reply.
let first=mount('draft');
first.children.get('.sf-value').onkeydown(event({key:'Enter'}));
const input=first.children.get('.sf-input');
input.value='1200 uW';input.selectionStart=4;input.selectionEnd=4;input.oninput();
first.dispose();first.root.isConnected=false;
let next=mount('draft',{value:0.6,text:'0.6'});
await Promise.resolve();
assert.equal(next.children.get('.sf-input').value,'1200 uW');
assert.equal(next.children.get('.sf-input').focused,true);
assert.equal(next.children.get('.sf-input').selectionStart,4);
assert.equal(sent.length,0);
next.children.get('.sf-input').onkeydown(event({key:'Enter'}));
assert.equal(sent.at(-1).t,'1200 uW');
next.dispose();next.root.isConnected=false;
// A rejected text response still acknowledges its request and restores the knob.
let rejected=mount('draft',{ack:1,value:0.6,text:'0.6',note:'expected mW'});
assert.equal(rejected.root._sf.pending,0);
assert.equal(rejected.children.get('.sf-num').textContent,'0.6');
rejected.dispose();rejected.root.isConnected=false;
// Release reaches the new render after the element that captured it is removed.
first=mount('pointer');
first.children.get('.sf-track').onpointerdown(event());
pointer('pointermove',75);
first.dispose();first.root.isConnected=false;
next=mount('pointer');
assert.equal(next.root._sf.shown,1.5);
pointer('pointerup',90);
assert.equal(sent.at(-1).v,1.8);
assert.equal(sent.at(-1).root,next.root);
assert.equal(sent.filter(item=>item.key==='pointer').length,1);
// An older reply cannot replace a pending value; its own acknowledgement can.
next.dispose();next.root.isConnected=false;
let oldReply=mount('pointer');
assert.equal(oldReply.root._sf.shown,1.8);
oldReply.dispose();oldReply.root.isConnected=false;
let acknowledged=mount('pointer',{ack:1,value:1.8,text:'1.8'});
assert.equal(acknowledged.root._sf.pending,0);
acknowledged.dispose();acknowledged.root.isConnected=false;
// The replacement owns the debounce: exactly one commit, from its callback.
first=mount('keyboard');
first.children.get('.sf-value').onkeydown(event({key:'ArrowUp'}));
first.dispose();first.root.isConnected=false;
next=mount('keyboard');
await wait(200);
assert.equal(sent.filter(item=>item.key==='keyboard').length,1);
assert.equal(sent.at(-1).v,0.6);
assert.equal(sent.at(-1).root,next.root);
next.dispose();next.root.isConnected=false;
// Graph-owned updates keep the knob synchronized through intermediate replies.
let owner={isConnected:true};
const plot = (phase,value) => window.dispatchEvent(new CustomEvent('gabes-plot-status',{
  detail:{eventKey:'fwm__tpd',phase,values:{tpd:value},owner}
}));
first=mount('graph-knob',{previewKey:'fwm__tpd'});
plot('syncing',1.4);
assert.equal(first.root._sf.shown,1.4);
first.dispose();first.root.isConnected=false;
next=mount('graph-knob',{previewKey:'fwm__tpd'});
assert.equal(next.root._sf.shown,1.4);
plot('ready',1.4);
assert.equal(next.root._sf.previewPending,false);
plot('preview',0.8);
next.children.get('.sf-value').onkeydown(event({key:'Enter'}));
next.children.get('.sf-input').value='1.2';
next.children.get('.sf-input').oninput();
plot('syncing',0.9);
assert.equal(next.children.get('.sf-input').value,'1.2');
next.children.get('.sf-input').onkeydown(event({key:'Escape'}));
owner.isConnected=false;
for(const observer of observers)observer.callback();
assert.equal(next.root._sf.previewPending,false);
assert.equal(next.root._sf.shown,0.5);
next.dispose();next.root.isConnected=false;
// True unmounts retire the input rather than resurrecting a later stale draft.
first=mount('removed');
first.children.get('.sf-value').onkeydown(event({key:'Enter'}));
first.dispose();first.root.isConnected=false;
await wait(120);
assert.equal(scrubStateForComponent('removed',defaults).editing,false);
console.log('ok');
'''.replace("MODULE_URL", json.dumps(url))
    result = subprocess.run([node, "--input-type=module"], input=script, text=True,
                            capture_output=True, check=True, timeout=15)
    assert result.stdout.strip() == "ok"
