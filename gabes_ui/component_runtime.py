"""Recognize a missing v2 registry only in Streamlit's test runner."""
import sys
from pathlib import Path


def spectrum_component_js(asset_name):
    """Bundle common pointer handling with a Streamlit v2 inline ES module."""
    assets = Path(__file__).resolve().parents[1] / "assets" / "ui"
    names = ["spectrum_gestures.js"]
    if asset_name in {"model_compare.js", "spectral_compare.js"}:
        names.append("spectrum_plot.js")
    if asset_name == "spectral_compare.js":
        names.append("spectrum_map.js")
    names.append(asset_name)
    return "\n".join((assets / name).read_text(encoding="utf-8") for name in names)


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
