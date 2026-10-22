"""Recognize a missing v2 registry only in Streamlit's test runner."""
import sys


def missing_apptest_component(error, component_name):
    frame = error.__traceback__
    while frame is not None and frame.tb_next is not None:
        frame = frame.tb_next
    if frame is None or not isinstance(error, ValueError):
        return False
    local = frame.tb_frame.f_locals
    ctx = local.get("ctx")
    runner = getattr(getattr(ctx, "_enqueue", None), "__self__", None)
    testing = sys.modules.get("streamlit.testing.v1.local_script_runner")
    runner_type = getattr(testing, "LocalScriptRunner", None)
    return (
        frame.tb_frame.f_code.co_name == "_bidi_component"
        and frame.tb_frame.f_globals.get("__name__")
        == "streamlit.components.v2.bidi_component.main"
        and local.get("component_def") is None
        and local.get("component_name") == component_name
        and str(error) == f"Component '{component_name}' is not registered"
        and isinstance(runner_type, type)
        and isinstance(runner, runner_type)
    )
