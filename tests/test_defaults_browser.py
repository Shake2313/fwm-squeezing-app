"""Browser defaults transport boundaries (scheme validation lives separately)."""
import json
import shutil
import subprocess
from types import SimpleNamespace

import pytest

from gabes_ui import defaults_browser as browser


def event(**updates):
    value = {"id": "request-1", "action": "save", "records": [], "error": None}
    return dict(value, **updates)


@pytest.mark.parametrize("value", [None, [], "{bad json", event(id=""),
    event(action="execute"), event(records={}), event(records=[None]),
    event(error={}), "x" * (3 * 1024 * 1024 + 1)],
    ids=["none", "array", "json", "id", "action", "records", "record", "error", "size"])
def test_rejects_malformed_events(value):
    assert browser.decode_event(value) is None


def test_accepts_json_acknowledgement_and_storage_failure():
    saved = event(records=[{"name": "My experiment"}])
    assert browser.decode_event(json.dumps(saved)) == saved
    failed = event(error="Browser storage is disabled.")
    assert browser.decode_event(failed) == failed


def test_mount_receives_ack_before_rerun_and_uses_stable_load_token(monkeypatch):
    state, calls, received = {}, [], []
    monkeypatch.setattr(browser, "st", SimpleNamespace(session_state=state))

    def component(**kwargs):
        calls.append(kwargs)
        state[kwargs["key"]] = {"event": json.dumps(event())}
        kwargs["on_event_change"]()

    monkeypatch.setattr(browser, "_COMPONENT", component)
    command = {"id": "request-1", "action": "save", "record": {"name": "A"}}
    assert browser.render_browser_store(None, "fwm", command, received.append) == event()
    browser.render_browser_store(None, "fwm")
    assert received == [event()]
    assert calls[0]["data"]["command"] == command
    assert calls[0]["data"]["readToken"] == calls[1]["data"]["readToken"]
    assert calls[0]["height"] == 0
    browser.render_browser_store(None, "sas")
    assert calls[2]["data"]["readToken"] != calls[0]["data"]["readToken"]


def test_unsupported_runtime_keeps_local_import_available(monkeypatch):
    monkeypatch.setattr(browser, "_COMPONENT", None)
    assert browser.render_browser_store(None, "fwm") is None


def test_component_mount_fallback(monkeypatch):
    monkeypatch.setattr(browser, "st", SimpleNamespace(session_state={}))

    def unavailable(**kwargs):
        raise AttributeError("AppTest component registry is unavailable")

    monkeypatch.setattr(browser, "_COMPONENT", unavailable)
    assert browser.render_browser_store(None, "fwm") is None


def test_unrelated_type_errors_remain_visible(monkeypatch):
    monkeypatch.setattr(browser, "st", SimpleNamespace(session_state={}))

    def broken(**kwargs):
        raise TypeError("a real integration problem")

    monkeypatch.setattr(browser, "_COMPONENT", broken)
    with pytest.raises(TypeError, match="real integration problem"):
        browser.render_browser_store(None, "fwm")


def test_apptest_component_registry_falls_back():
    from streamlit.testing.v1 import AppTest
    app = AppTest.from_string(
        'from gabes_ui.defaults_browser import render_browser_store\n'
        'import streamlit as st\n'
        'render_browser_store(st,"fwm")\n'
        'st.write("Local import remains available")\n'
    ).run()
    assert not app.exception
    assert app.markdown[0].value == "Local import remains available"


@pytest.mark.skipif(shutil.which("node") is None, reason="Node unavailable for JS bridge check")
def test_lost_acknowledgements_replay_after_remount_without_rewriting():
    script = browser._JS.replace(
        "export default function (component)", "function bridge(component)") + r'''
const events = [], stored = new Map();
let writes = 0;
globalThis.window = {localStorage: {
  getItem: key => stored.get(key) ?? null,
  setItem: (key, value) => {writes += 1; stored.set(key, value);}
}};
const run = data => bridge({data, setStateValue:(key,value)=>events.push(JSON.parse(value))});
const record = {name:'Lab A',scheme:'fwm',params:{x:1}};
const pending = {scheme:'fwm',readToken:'read1',command:{id:'save1',action:'save',record}};
run(pending);
const lostSave = events.pop();
run({scheme:'sas',readToken:'read2',command:null});
run(pending); // A new component instance returns to the pending FWM command.
if (JSON.stringify(events.at(-1)) !== JSON.stringify(lostSave)) throw Error('lost success acknowledgement');
if (writes !== 1) throw Error('replayed save wrote again');
let failures = 0;
window.localStorage.setItem = () => {failures += 1;const error=new Error('full');error.name='QuotaExceededError';throw error;};
const failed = {...pending,command:{id:'save2',action:'save',record:{...record,name:'Lab B'}}};
run(failed);
const lostError = events.pop();
run(failed);
if (!lostError.error || JSON.stringify(events.at(-1)) !== JSON.stringify(lostError)) throw Error('lost error acknowledgement');
if (failures !== 1) throw Error('replayed error retried write');
if (events.at(-1).records[0].name !== 'Lab A') throw Error('failed save changed stored list');
'''
    result = subprocess.run(
        [shutil.which("node"), "--input-type=module"], input=script,
        text=True, capture_output=True, timeout=15,
    )
    assert result.returncode == 0, result.stderr


def test_unregistered_real_runtime_component_error_remains_visible(monkeypatch):
    """The same real Streamlit registry exception is not an app fallback."""
    import importlib
    from streamlit.runtime import Runtime

    if browser._COMPONENT is None:
        pytest.skip("Streamlit components v2 unavailable")
    bidi = importlib.import_module("streamlit.components.v2.bidi_component.main")
    context = SimpleNamespace(_enqueue=lambda message: None)
    runtime = SimpleNamespace(bidi_component_registry={})
    monkeypatch.setattr(bidi, "get_script_run_ctx", lambda: context)
    monkeypatch.setattr(Runtime, "instance", lambda: runtime)
    monkeypatch.setattr(browser, "st", SimpleNamespace(session_state={}))
    with pytest.raises(ValueError, match="Component 'gabes_defaults_store' is not registered"):
        browser.render_browser_store(None, "fwm")


@pytest.mark.parametrize("failure_source", ["component", "callback"])
def test_apptest_real_value_errors_remain_visible(monkeypatch, failure_source):
    """Even this exact message cannot hide a component or event callback bug."""
    from streamlit.testing.v1 import AppTest

    message = "Component 'gabes_defaults_store' is not registered"

    def component(**kwargs):
        if failure_source == "component":
            raise ValueError(message)
        browser.st.session_state[kwargs["key"]] = {"event": json.dumps(event())}
        kwargs["on_event_change"]()

    monkeypatch.setattr(browser, "_COMPONENT", component)
    app = AppTest.from_string(
        'from gabes_ui import defaults_browser as browser\n'
        'import streamlit as st\n'
        'def fail(event):\n'
        '    raise ValueError("Component \'gabes_defaults_store\' is not registered")\n'
        'browser.render_browser_store(st, "fwm", on_event=fail)\n'
    ).run(timeout=10)
    assert len(app.exception) == 1
    assert app.exception[0].message == message
