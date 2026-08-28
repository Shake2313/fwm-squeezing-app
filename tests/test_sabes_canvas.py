"""
Canvas component checks (Stage B of the optical-table plan).

    python tests/test_sabes_canvas.py   # or: pytest tests/test_sabes_canvas.py

What can be pinned here is the Python side and the shipped frontend contract:
spec construction, package data being present, and the postMessage protocol the
hand-written frontend must implement. The click round-trip itself is a browser
behaviour and was verified interactively against the running app -- see
docs/SABES/optical_table_ui_plan.md.
"""
import re
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from sabes import components  # noqa: E402
from sabes.components import canvas as canvas_module  # noqa: E402

FRONTEND = ROOT / "sabes" / "components" / "frontend" / "index.html"


# --------------------------------------------------------------- packaging

def test_frontend_ships_as_package_data():
    """No build step: the component is a static file next to its module."""
    assert FRONTEND.is_file()
    assert canvas_module._FRONTEND_DIR == FRONTEND.parent


def test_no_build_tooling_crept_in():
    """The whole point of hand-writing the protocol is having no toolchain."""
    directory = FRONTEND.parent
    for forbidden in ("package.json", "package-lock.json", "node_modules",
                      "yarn.lock", "vite.config.js"):
        assert not (directory / forbidden).exists(), forbidden


def test_importing_the_package_does_not_shadow_the_module():
    """`canvas` is a module; the drawing function is `svg_canvas`.

    Binding the function to the package attribute `canvas` would make
    `from sabes.components import canvas` silently return a function, which is
    exactly the bug this naming avoids.
    """
    import types
    assert isinstance(canvas_module, types.ModuleType)
    assert callable(components.svg_canvas)
    assert not hasattr(components, "canvas") or isinstance(
        components.canvas, types.ModuleType)


# --------------------------------------------------------------- protocol

def test_frontend_implements_the_three_protocol_messages():
    source = FRONTEND.read_text(encoding="utf-8")
    for message in ("streamlit:componentReady",
                    "streamlit:setFrameHeight",
                    "streamlit:setComponentValue"):
        assert message in source, message
    assert "streamlit:render" in source
    assert "isStreamlitMessage" in source


def test_frontend_guards_against_swallowed_repeat_clicks():
    """Streamlit only reruns when a component value CHANGES.

    Without a monotonic sequence number, clicking the same optic twice would do
    nothing the second time and read as a broken UI.
    """
    source = FRONTEND.read_text(encoding="utf-8")
    assert re.search(r"seq\s*\+=\s*1", source)
    assert "payload.seq = seq" in source


def test_frontend_converts_clicks_to_spec_coordinates():
    """A click must mean the same thing regardless of how the SVG was scaled."""
    source = FRONTEND.read_text(encoding="utf-8")
    assert "getScreenCTM" in source and "inverse()" in source


def test_frontend_reports_height_after_parent_resize():
    """The production table must recover when a hidden tab becomes visible."""
    source = FRONTEND.read_text(encoding="utf-8")
    assert 'window.addEventListener("resize", reportHeight)' in source
    assert "ResizeObserver(reportHeight)" in source


# --------------------------------------------------------------- spec model

def test_make_spec_has_the_keys_the_frontend_reads():
    spec = components.make_spec([], width=800, height=250)
    assert set(spec) == {"width", "height", "background", "shapes", "probe"}
    assert spec["width"] == 800.0 and spec["height"] == 250.0
    assert spec["probe"] is None


def test_rect_and_circle_are_clickable_by_default_and_beams_never_are():
    assert components.rect("pbs", 1, 2, 3, 4)["clickable"] is True
    assert components.circle("lens", 1, 2, 3)["clickable"] is True
    # Clicking a beam means clicking the table under it -- that is how an
    # instrument head gets placed on the beam.
    assert components.beam([(0, 0), (1, 1)])["clickable"] is False


def test_decorative_optics_can_opt_out_of_being_clickable():
    """Steering mirrors are geometry, not parameters; the table must say so."""
    mirror = components.rect("m1", 0, 0, 10, 3, clickable=False)
    assert mirror["clickable"] is False


def test_spec_is_json_serialisable():
    """It crosses the iframe boundary as JSON, so numpy scalars would break it."""
    import json
    spec = components.make_spec(
        [components.beam([(0, 1), (2, 3)]),
         components.rect("cell", 5, 6, 7, 8, angle=45, label="Rb"),
         components.circle("iris", 9, 10, 2)],
        probe=(4, 5))
    restored = json.loads(json.dumps(spec))
    assert restored["shapes"][1]["angle"] == 45.0
    assert restored["probe"] == [4.0, 5.0]


def test_geometry_values_are_coerced_to_float():
    shape = components.rect("x", 1, 2, 3, 4)
    assert all(isinstance(shape[k], float) for k in ("x", "y", "w", "h", "angle"))
    path = components.beam([(0, 1), (2, 3)])
    assert path["points"] == [[0.0, 1.0], [2.0, 3.0]]


def test_missing_frontend_is_reported_clearly(monkeypatch):
    monkeypatch.setattr(canvas_module, "_component", None)
    monkeypatch.setattr(canvas_module, "_FRONTEND_DIR", ROOT / "does-not-exist")
    with pytest.raises(FileNotFoundError, match="package data"):
        canvas_module._instance()


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-q"]))
