"""Named parameter snapshots: built-ins, browser storage and portable JSON."""
from uuid import uuid4

import streamlit as st

from gabes_ui.controls import current_params
from gabes_ui.defaults_browser import render_browser_store
from gabes_ui.defaults_store import (
    DefaultValidationError, create_default, default_filename, dump_default,
    parse_default,
)


def _key(scheme, name):
    return f"{scheme.name}___{name}"


def _params(scheme):
    # Defaults render before the rail. Required modes must already be restored
    # here so the first save/download snapshot cannot retain a stale None.
    return current_params(scheme.name, scheme)


def builtin_defaults(scheme):
    """Complete snapshots prevent settings from the previous preset leaking in."""
    presets = scheme.presets()
    if not presets:
        return {"builtin:0": create_default(scheme, "Default settings", scheme.defaults())}
    return {f"builtin:{i}": create_default(
                scheme, preset.name, dict(scheme.defaults(), **preset.values))
            for i, preset in enumerate(presets)}


def _saved(scheme):
    return st.session_state.setdefault(_key(scheme, "saved_defaults"), {})


def _notice(scheme, message, error=False):
    st.session_state[_key(scheme, "defaults_notice")] = (message, error)


def _apply(scheme, record):
    for name, value in record["params"].items():
        st.session_state[f"{scheme.name}__{name}"] = value


def selector_state_key(scheme):
    revision = st.session_state.get(_key(scheme, "defaults_revision"), 0)
    base = _key(scheme, "named_default")
    return f"{base}:{revision}" if revision else base


def _activate(scheme, ident, refresh=False):
    # Keep the chosen name outside widget state. In a browser acknowledgement
    # rerun, the selectbox can still carry the previous client-side selection.
    if refresh:
        revision_key = _key(scheme, "defaults_revision")
        st.session_state[revision_key] = st.session_state.get(revision_key, 0) + 1
    st.session_state[_key(scheme, "active_default")] = ident
    st.session_state[selector_state_key(scheme)] = ident


def _receive(scheme, event):
    """Storage acknowledgement arrives as a widget callback, before controls."""
    pending_key = _key(scheme, "defaults_command")
    command = st.session_state.get(pending_key)
    acknowledged = command and command["id"] == event.get("id")
    if acknowledged:
        st.session_state.pop(pending_key, None)
    if event.get("error"):
        _notice(scheme, "Browser storage unavailable: " + str(event["error"]), True)
        return
    builtins = {record["name"] for record in builtin_defaults(scheme).values()}
    invalid = 0
    for raw in event.get("records", []):
        try:
            record = parse_default(scheme, raw)
            if record["name"] in builtins:
                raise DefaultValidationError("A built-in default cannot be replaced.")
            _saved(scheme)[record["name"]] = record
        except DefaultValidationError:
            invalid += 1
    if acknowledged:
        _activate(scheme, "saved:" + command["record"]["name"], refresh=True)
        _notice(scheme, f'Saved “{command["record"]["name"]}” in this browser.')
    elif invalid:
        _notice(scheme, f"Skipped {invalid} incompatible browser default(s).", True)


def _named_record(scheme, name, builtins):
    record = create_default(scheme, name, _params(scheme))
    if record["name"] in {r["name"] for r in builtins.values()}:
        raise DefaultValidationError("Choose a name different from a built-in default.")
    return record


def _save_browser(scheme, name_key, builtins):
    try:
        record = _named_record(scheme, st.session_state[name_key], builtins)
    except DefaultValidationError as exc:
        _notice(scheme, str(exc), True)
        return
    st.session_state[_key(scheme, "defaults_command")] = {
        "id": uuid4().hex, "action": "save", "record": record,
    }


def _import(scheme, upload_key, builtins):
    upload = st.session_state.get(upload_key)
    if upload is None:
        return
    try:
        record = parse_default(scheme, upload.getvalue())
        if record["name"] in {r["name"] for r in builtins.values()}:
            raise DefaultValidationError("The file uses a reserved built-in name.")
    except DefaultValidationError as exc:
        _notice(scheme, "Cannot import default: " + str(exc), True)
        return
    _saved(scheme)[record["name"]] = record
    _apply(scheme, record)
    _activate(scheme, "saved:" + record["name"], refresh=True)
    st.session_state[_key(scheme, "default_name")] = record["name"]
    _notice(scheme, f'Imported “{record["name"]}”. Save in browser to keep it after reload.')


def render_defaults(container, scheme):
    """Render before any parameter widgets so loading also restores mode safely."""
    builtins = builtin_defaults(scheme)
    render_browser_store(
        container, scheme.name,
        command=st.session_state.get(_key(scheme, "defaults_command")),
        on_event=lambda event: _receive(scheme, event),
    )
    records = dict(builtins)
    records.update({"saved:" + name: record for name, record in _saved(scheme).items()})
    current = _params(scheme)
    selector_key = selector_state_key(scheme)
    active_key = _key(scheme, "active_default")
    selected = st.session_state.get(active_key, st.session_state.get(selector_key))
    if selected not in records or records[selected]["params"] != current:
        selected = next((ident for ident, record in records.items()
                         if record["params"] == current), None)
    st.session_state[active_key] = selected
    st.session_state[selector_key] = selected

    def _selected():
        if selector_key != selector_state_key(scheme):
            return  # A response from the selector replaced by save/import.
        record = records.get(st.session_state.get(selector_key))
        if record:
            _activate(scheme, st.session_state[selector_key])
            _apply(scheme, record)

    container.markdown("<div class='gabes-group-header'>DEFAULTS</div>",
                       unsafe_allow_html=True)
    container.selectbox(
        "Defaults", list(records), key=selector_key, index=None,
        format_func=lambda ident: records[ident]["name"],
        placeholder="Current settings", label_visibility="collapsed",
        on_change=_selected, help="Load a complete named set of settings.",
    )
    actions = container.container(key="gabes_defaults_actions")
    left, right = actions.columns([3, 2], gap="small")
    with left.popover("Save as default", use_container_width=True):
        name_key = _key(scheme, "default_name")
        st.text_input("Default name", value="My default", key=name_key, max_chars=80)
        try:
            record = _named_record(scheme, st.session_state[name_key], builtins)
            error = None
        except DefaultValidationError as exc:
            record, error = None, str(exc)
        replacing = record and record["name"] in _saved(scheme)
        st.button(
            "Replace in browser" if replacing else "Save in browser",
            key=_key(scheme, "save_browser"), use_container_width=True,
            disabled=record is None or bool(st.session_state.get(_key(scheme, "defaults_command"))),
            on_click=_save_browser, args=(scheme, name_key, builtins),
        )
        st.caption("Saved in this browser’s cache. Clearing site data may remove it.")
        st.download_button(
            "Save locally", dump_default(record) if record else "",
            file_name=default_filename(record["name"]) if record else "default.json",
            mime="application/json", key=_key(scheme, "download_default"),
            disabled=record is None, use_container_width=True, on_click="ignore",
            help="Download the current settings as a named JSON file.",
        )
        if error:
            st.caption(error)
    with right.popover("Import", use_container_width=True):
        upload_key = _key(scheme, "import_default")
        st.file_uploader(
            "Import local default", type=["json"], key=upload_key,
            max_upload_size=1, on_change=_import, args=(scheme, upload_key, builtins),
            help="Load a JSON default exported from this scheme.",
        )
    notice = st.session_state.get(_key(scheme, "defaults_notice"))
    if notice:
        message, error = notice
        if error:
            container.warning(message)
        else:
            container.caption(message)
