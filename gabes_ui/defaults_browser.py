"""Small browser-local storage bridge for named parameter defaults.

The native controls own the interface and validate records against the scheme.
This component only transports JSON to and from localStorage on the app origin.
No parameter values are applied here, and storage failures are reported back to
the caller without replacing the existing stored document.
"""
import contextlib
import json
import uuid

import streamlit as st
from .component_runtime import missing_apptest_component


_HTML = '<span class="gabes-defaults-store" aria-hidden="true"></span>'
_CSS = ".gabes-defaults-store{display:none;}"
_JS = r"""
export default function (component) {
  const { data, setStateValue } = component;
  const key = "gabes.defaults.v1." + data.scheme;
  // Components can remount during a rerun. Keep acknowledgement state for this
  // page, rather than asking Python to rerun again for every remount.
  const stores = window.__gabesDefaultsStores || (window.__gabesDefaultsStores = new Map());
  const state = stores.get(key) || { reads: new Set(), commands: new Map() };
  // Retain acknowledgements as well as IDs: an unmount can discard an event
  // after localStorage was updated but before Python received the response.
  if (!(state.commands instanceof Map)) state.commands = new Map();
  stores.set(key, state);

  const read = () => {
    const text = window.localStorage.getItem(key);
    if (text == null) return Object.create(null);
    if (text.length > 2097152) throw new Error("Saved defaults are too large to load.");
    const value = JSON.parse(text);
    if (!value || typeof value !== "object" || Array.isArray(value)) {
      throw new Error("Saved defaults have an invalid storage format.");
    }
    const records = Object.create(null);
    for (const name of Object.keys(value)) {
      const record = value[name];
      if (!record || typeof record !== "object" || Array.isArray(record) ||
          record.name !== name || record.scheme !== data.scheme) {
        throw new Error("A saved default is invalid; existing browser data was kept.");
      }
      records[name] = record;
    }
    return records;
  };
  const send = (id, action, records, error) => {
    const acknowledgement = {
      id, action, records: Object.values(records), error: error || null
    };
    if (action !== "load") state.commands.set(id, acknowledgement);
    setStateValue("event", JSON.stringify(acknowledgement));
  };
  const errorText = (error) => {
    if (error && error.name === "QuotaExceededError") {
      return "Browser storage is full. Save locally or remove a saved default.";
    }
    if (error && error.name === "SecurityError") {
      return "Browser storage is disabled. Save locally instead.";
    }
    return error && error.message ? error.message : "Browser storage is unavailable.";
  };
  const command = data.command;
  if (command && state.commands.has(command.id)) {
    // Python clears a pending command after its acknowledgement. Replay until
    // it does so; repeating the identical state value causes no extra rerun.
    setStateValue("event", JSON.stringify(state.commands.get(command.id)));
    return;
  }
  if (command && !state.commands.has(command.id)) {
    state.reads.add(data.readToken);
    let records = Object.create(null);
    try {
      records = read();
      const next = Object.assign(Object.create(null), records);
      if (command.action === "save") {
        const record = command.record;
        if (!record || typeof record.name !== "string" || !record.name.trim() ||
            record.scheme !== data.scheme) throw new Error("Invalid default to save.");
        next[record.name] = record;
      } else if (command.action === "delete") {
        if (typeof command.name !== "string") throw new Error("Invalid default to remove.");
        delete next[command.name];
      } else {
        throw new Error("Unknown browser storage action.");
      }
      const text = JSON.stringify(next);
      if (text.length > 2097152) throw new Error("Saved defaults exceed the browser storage limit.");
      // Publish only after a successful write. A failed save retains the list
      // actually stored in the browser, never the proposed replacement.
      window.localStorage.setItem(key, text);
      records = next;
      send(command.id, command.action, records, null);
    } catch (error) {
      send(command.id, command.action, records, errorText(error));
    }
  } else if (!state.reads.has(data.readToken)) {
    state.reads.add(data.readToken);
    try {
      send("load-" + data.readToken, "load", read(), null);
    } catch (error) {
      send("load-" + data.readToken, "load", Object.create(null), errorText(error));
    }
  }
}
"""


try:
    _COMPONENT = st.components.v2.component(
        "gabes_defaults_store", html=_HTML, css=_CSS, js=_JS)
except (AttributeError, RuntimeError):  # older Streamlit without components v2
    _COMPONENT = None


def _component_key(scheme_name):
    return "defaults-store-" + scheme_name.replace("__", "-")


def decode_event(value):
    """Parse a component event without trusting malformed browser payloads."""
    if isinstance(value, str):
        if len(value) > 3 * 1024 * 1024:
            return None
        try:
            value = json.loads(value)
        except (ValueError, TypeError):
            return None
    if not isinstance(value, dict):
        return None
    if not isinstance(value.get("id"), str) or not value["id"]:
        return None
    if value.get("action") not in ("load", "save", "delete"):
        return None
    if not isinstance(value.get("records"), list):
        return None
    if any(not isinstance(record, dict) for record in value["records"]):
        return None
    if value.get("error") is not None and not isinstance(value["error"], str):
        return None
    return value


def _state_event(comp_key):
    state = st.session_state.get(comp_key)
    if state is None:
        return None
    try:
        return decode_event(state["event"])
    except (KeyError, TypeError):
        return None


def _unsupported_test_registry(error):
    """Recognize the two known AppTest-only component registry failures."""
    if isinstance(error, ValueError):
        return missing_apptest_component(error, "gabes_defaults_store")
    frame = error.__traceback__
    while frame is not None:
        if frame.tb_frame.f_code.co_name == "_bidi_component":
            local = frame.tb_frame.f_locals
            definition = local.get("component_def")
            javascript = getattr(definition, "js_content", None)
            return javascript is not None and not isinstance(javascript, str)
        frame = frame.tb_next
    return False


def render_browser_store(container, scheme_name, command=None, on_event=None):
    """Mount the storage bridge; ``on_event(event)`` runs before the rerun.

    ``command`` is None or a dict with a unique ``id``, ``action`` of ``save``
    (with ``record``) or ``delete`` (with ``name``). Every acknowledgement
    includes the resulting ``records`` list and an optional ``error`` string.
    Loading storage never applies a default. AppTest and old runtimes may not
    mount v2 components; the local JSON save/import path remains available.
    """
    if _COMPONENT is None:
        return None
    comp_key = _component_key(scheme_name)
    token_key = "_" + comp_key + "-read-token"
    if token_key not in st.session_state:
        st.session_state[token_key] = uuid.uuid4().hex

    def receive():
        event = _state_event(comp_key)
        if event is not None and on_event is not None:
            on_event(event)

    ctx = contextlib.nullcontext() if container is None or container is st else container
    try:
        with ctx:
            _COMPONENT(
                key=comp_key,
                data={"scheme": scheme_name, "readToken": st.session_state[token_key],
                      "command": command},
                default={"event": None}, height=0, width="stretch",
                on_event_change=receive,
            )
    except (AttributeError, RuntimeError):
        # AppTest has no component registry. Restrict this fallback to mounting
        # errors; errors inside the event callback must remain visible.
        return None
    except (TypeError, ValueError) as error:
        if not _unsupported_test_registry(error):
            raise
        return None
    return _state_event(comp_key)
