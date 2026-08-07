"""
SABES page wiring: control coverage, the three compute tiers, and routing.

    python tests/test_sabes_page.py   # or: pytest tests/test_sabes_page.py

The tier claim is tested through `bridge.solve_key`, not through wall-clock
timing: a knob belongs to tier A or C exactly when it leaves the solve key alone,
and that is a property of the code rather than of the machine it runs on.
"""
import sys
from dataclasses import fields, replace
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from gabes.schemes import fwm  # noqa: E402
from sabes import bridge  # noqa: E402
from sabes.beamline import SetupSettings  # noqa: E402
from sabes.detection import DetectionSettings  # noqa: E402

import sabes_page  # noqa: E402


# ------------------------------------------------------------------ controls

def test_every_control_names_a_real_setting_field():
    known = {f.name for f in fields(SetupSettings)}
    known |= {f.name for f in fields(DetectionSettings)}
    known |= {"eom_offset_mhz"} | set(sabes_page.ETALON_KEYS)
    for name, (label, unit, low, high, step, help_text) in \
            sabes_page.CONTROLS.items():
        assert name in known, f"{name} is not a settings field"
        assert label and help_text, f"{name} is missing label or help"
        assert low < high and step > 0


def test_defaults_cover_every_control_and_land_on_the_paper_point():
    defaults = sabes_page._defaults()
    for name in sabes_page.CONTROLS:
        assert name in defaults, name
    # The generator control is an offset from the hyperfine splitting, and the
    # default offset is what produces the -8 MHz operating point.
    assert defaults["eom_offset_mhz"] == pytest.approx(8.0, abs=1e-6)
    assert defaults["resolution"] in sabes_page.FIDELITIES
    # The shared etalon knob is gone; each stage carries its own.
    assert "etalon_detune_ghz" not in defaults
    for name in sabes_page.ETALON_KEYS:
        assert name in defaults


def test_every_control_is_owned_by_an_optic_on_the_drawing():
    """Stage D's whole point: the sidebar no longer holds the parameters.

    A knob no optic claims would be unreachable once the sidebar stops listing
    them, so this is the check that keeps the drawing complete.
    """
    orphans = set(sabes_page.CONTROLS) - sabes_page._owned_parameters()
    assert not orphans, f"no optic owns: {sorted(orphans)}"


def test_the_sidebar_no_longer_carries_the_parameter_wall():
    """Parameters live on their optic; the sidebar keeps only global controls.

    Checked as a property rather than a widget count: the sidebar may render a
    parameter *only* as the orphan fallback, which a healthy build never enters.
    """
    source = (ROOT / "sabes_page.py").read_text(encoding="utf-8")
    sidebar = source[source.index("def _render_sidebar"):
                     source.index("def _owned_parameters")]
    assert "st.number_input" not in sidebar
    assert sidebar.count("_render_control") == 1        # the fallback only
    assert sidebar.index("if orphans:") < sidebar.index("_render_control")
    assert "Click it in the" in sidebar          # points at the drawing instead
    # And the fallback is empty in this build, so nothing is actually listed.
    assert not set(sabes_page.CONTROLS) - sabes_page._owned_parameters()


# --------------------------------------------------------------- the tiers

def _params(settings=None, detection=None, **overrides):
    return bridge.run(settings or SetupSettings(),
                      detection=detection or DetectionSettings(),
                      solve=False, **overrides).params


def test_detection_knobs_never_touch_the_solve_key():
    """Tier C: lenses, defocus and iris must not trigger a Bloch solve.

    The iris does change `loss_pct`, which IS a solve knob -- but only when it
    actually clips. Wide open, moving it is free.
    """
    base = bridge.solve_key(_params())
    for detection in (DetectionSettings(probe_lens_focal_mm=400.0),
                      DetectionSettings(conjugate_lens_focal_mm=250.0),
                      DetectionSettings(pd_defocus_mm=30.0),
                      DetectionSettings(pump_leakage_dbm=-40.0)):
        assert bridge.solve_key(_params(detection=detection)) == base, detection


def test_pump_leakage_is_a_readout_input_not_a_solve_input():
    """It changes what the analyser shows, never what the atoms do."""
    base = bridge.solve_key(_params())
    leaky = bridge.solve_key(
        _params(detection=DetectionSettings(pump_leakage_dbm=-30.0)))
    assert leaky == base


def test_cell_facing_knobs_do_change_the_solve_key():
    base = bridge.solve_key(_params())
    for settings in (replace(SetupSettings(), cell_temp_c=125.0),
                     replace(SetupSettings(), opd_ghz=1.2),
                     replace(SetupSettings(), hwp_split_deg=14.0),
                     replace(SetupSettings(), seed_telescope_f1_mm=250.0),
                     replace(SetupSettings(), dmirror_separation_mm=5.0)):
        assert bridge.solve_key(_params(settings)) != base, settings


def test_the_generator_is_navigate_only():
    """TPD moves the readout along an already-solved curve, so it must not
    re-solve -- the same two-tier split GABES uses for its own TPD slider."""
    base = bridge.solve_key(_params())
    moved = replace(SetupSettings(),
                    eom_frequency_hz=SetupSettings().eom_frequency_hz + 3e6)
    assert bridge.solve_key(_params(moved)) == base
    assert _params(moved)["tpd"] != _params()["tpd"]


def test_quantisation_absorbs_knob_noise_without_hiding_real_change():
    """A nudge far below any physical significance is a cache hit."""
    base = bridge.solve_key(_params())
    nudged = replace(SetupSettings(),
                     hwp_split_deg=SetupSettings().hwp_split_deg + 1e-7)
    assert bridge.solve_key(_params(nudged)) == base
    real = replace(SetupSettings(),
                   hwp_split_deg=SetupSettings().hwp_split_deg + 0.05)
    assert bridge.solve_key(_params(real)) != base


def test_quantised_parameters_are_what_actually_gets_solved():
    """Displayed and solved values must be the same number, not near-equal."""
    params = _params()
    for key, step in bridge.PARAM_GRID.items():
        value = params[key]
        assert abs(value / step - round(value / step)) < 1e-6, key


def test_solve_key_holds_only_recompute_knobs():
    params = _params()
    keys = {k for k, _ in bridge.solve_key(params)}
    recompute = set(fwm.FWMScheme().recompute_keys())
    assert keys <= recompute | {"mode"}
    assert "tpd" not in keys
    for essential in ("opd", "temp_c", "pump_mw", "pump_waist_um", "qe_pct"):
        assert essential in keys


# ------------------------------------------ widget state vs. real state
#
# Streamlit widget state does not outlive the widget, and every editable knob on
# this page lives inside a dialog that is torn down as soon as it closes. These
# pin the separation that keeps a closing dialog from rewriting the setup.

def _numeric_input_calls():
    """Every `*.number_input(...)` call in the page, with its enclosing def."""
    import ast

    tree = ast.parse((ROOT / "sabes_page.py").read_text(encoding="utf-8"))
    found = []
    for parent in ast.walk(tree):
        if not isinstance(parent, (ast.FunctionDef, ast.AsyncFunctionDef)):
            continue
        for node in ast.walk(parent):
            if (isinstance(node, ast.Call)
                    and isinstance(node.func, ast.Attribute)
                    and node.func.attr == "number_input"):
                found.append((parent.name, node))
    return found


def test_no_numeric_input_is_keyed_on_a_setting():
    """The bug this exists to prevent, in one sentence.

    A `number_input` keyed straight on a setting loses that setting the moment
    its dialog closes: Streamlit leaves the key holding the widget's default,
    which for these inputs is `min_value`. Opening an optic and closing it again
    without touching anything used to set every knob on it to its minimum --
    silently, and a zeroed transmission then took the Bloch solve down with it.
    Every input therefore goes through `_edit_number`, which keys on the edit
    namespace and copies into the real state in `on_change`.
    """
    calls = _numeric_input_calls()
    assert calls, "expected the page to render numeric inputs"
    owners = {owner for owner, _ in calls}
    assert owners == {"_edit_number"}, (
        f"numeric inputs outside _edit_number: {sorted(owners - {'_edit_number'})}")


def test_the_edit_namespace_is_invisible_to_everything_that_reads_state():
    """A widget key must not look like a setting or a calibration coefficient.

    `_session_calibration` scans session state by prefix, so an edit-namespace
    key that matched would feed a torn-down widget's value into the physics.
    """
    for name in ("t_m_mopa_out", "ecdl_power_mw"):
        widget = sabes_page._widget_key(sabes_page._key(name))
        assert widget != sabes_page._key(name)
        assert not widget.startswith(sabes_page.TRANSMISSION_PREFIX)
        assert widget.startswith(sabes_page.WIDGET_PREFIX)


def test_every_optic_transmission_has_a_coefficient_to_edit():
    """`_edit_number` edits stored state, so the state has to exist.

    A node naming a `transmission_key` the shipped calibration does not carry
    would have nothing seeded for its dialog to point at.
    """
    from sabes import layout
    from sabes.calibration import default_calibration

    calibration = default_calibration()
    orphans = sorted({node.transmission_key
                      for part in layout.LAYOUTS for node in part.nodes
                      if node.transmission_key
                      and node.transmission_key not in calibration})
    assert not orphans, f"no calibration coefficient for: {orphans}"


def test_commit_copies_the_widget_value_into_the_real_state():
    state_key = sabes_page._key("ecdl_power_mw")
    session = {state_key: 40.0, sabes_page._widget_key(state_key): 12.5}
    original = sabes_page.st.session_state
    sabes_page.st.session_state = session
    try:
        sabes_page._commit_edit(state_key)
    finally:
        sabes_page.st.session_state = original
    assert session[state_key] == 12.5


# ------------------------------------------------- clicking an optic

def test_a_canvas_click_does_not_cost_a_second_script_run():
    """The drawing applies a click on its own side, so Python need not rerun.

    Selecting an optic used to run the whole page twice -- once to receive the
    click and once so the panels below could see it -- and each run redrew every
    tab. Re-reading the selection after the canvas is what removed the second.
    """
    source = (ROOT / "sabes_page.py").read_text(encoding="utf-8")
    body = source[source.index("def _render_table_tab"):
                  source.index("def _markdown_table")]
    handled = body.index("_handle_canvas_event")
    panels = body.index("Selected optic")
    assert "st.rerun()" not in body[handled:panels], \
        "handling a canvas click must not force a rerun"
    # ...and the panels below must be given the click's own value.
    assert body.count('st.session_state.get(_key("selected"))') >= 2

    frontend = (ROOT / "sabes" / "components" / "frontend" / "index.html").read_text(
        encoding="utf-8")
    # Both click kinds paint themselves before telling Python, which is what
    # makes skipping the rerun invisible to the user.
    assert "lastSpec.probe = [point.x, point.y]" in frontend
    assert frontend.count("render(lastSpec);") >= 2


def test_the_open_dialog_flag_is_consumed_not_read():
    """Otherwise a dialog dismissed with ESC pops back open on the next rerun."""
    source = (ROOT / "sabes_page.py").read_text(encoding="utf-8")
    assert 'st.session_state.pop(_key("open_optic"), None)' in source
    assert 'st.session_state.get(_key("open_optic"))' not in source


# ------------------------------------------------------- figure caching

def test_figures_are_cached_as_bytes_under_streamlit_image_width():
    """`st.pyplot` re-rasterises on every rerun; that was most of a click.

    Two separate costs: drawing the figure, and Streamlit re-opening anything
    wider than its content width to resize and re-encode it. Caching PNG bytes
    that already fit removes both.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from PIL import Image
    from io import BytesIO

    figure, axis = plt.subplots(figsize=(12.0, 4.0))     # wider than the cap
    axis.plot([0, 1], [0, 1])
    data = sabes_page._figure_png(figure)
    assert isinstance(data, bytes)
    image = Image.open(BytesIO(data))
    assert image.width <= sabes_page.MAX_IMAGE_WIDTH_PX, image.size

    source = (ROOT / "sabes_page.py").read_text(encoding="utf-8")
    assert "st.pyplot(" not in source, "figures must go through the PNG cache"


# ------------------------------------------------------------------ routing

def test_the_router_constant_matches_what_the_page_expects():
    source = (ROOT / "streamlit_app.py").read_text(encoding="utf-8")
    assert 'SABES_QUERY_VALUE = "sabes"' in source
    assert 'st.query_params.get("app") == SABES_QUERY_VALUE' in source
    # The router must sit before GABES touches the sidebar, or both pages draw.
    assert source.index("sabes_page.render") < source.index("st.sidebar.image")


def test_sabes_is_not_reachable_from_the_scheme_registry():
    """SABES is a separate site, not a scheme: it must stay out of the dropdown."""
    from gabes import schemes
    titles = " ".join(s.title for s in schemes.all_schemes()).lower()
    assert "sabes" not in titles
    assert "sabes" not in " ".join(schemes.REGISTRY).lower()


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-q"]))


# ------------------------------------------------- page/model drift

def test_every_model_attribute_the_page_reads_actually_exists():
    """Catches the class of bug a browser would find and a unit test would not.

    Renaming a field on the geometry or readout leaves the page compiling and
    the suite green, and blows up only when someone opens that tab. This walks
    the page's attribute accesses on the model objects and checks them for real.
    """
    import ast

    from sabes import bridge
    from sabes.beamline import build_source_chain
    from sabes.detection import DetectionSettings, geometry, readout

    settings = SetupSettings()
    detection_settings = DetectionSettings()
    chain = build_source_chain(settings)
    geom = geometry(chain, settings, detection_settings)
    read = readout(geom, chain, (14.0, 13.0), settings, detection_settings)
    live = {"geom": geom, "readout": read, "chain": chain,
            "detection": detection_settings, "settings": settings}

    source = (ROOT / "sabes_page.py").read_text(encoding="utf-8")
    missing = []
    for node in ast.walk(ast.parse(source)):
        if not isinstance(node, ast.Attribute):
            continue
        if not isinstance(node.value, ast.Name) or node.value.id not in live:
            continue
        target = live[node.value.id]
        if not hasattr(target, node.attr):
            missing.append(f"{node.value.id}.{node.attr} (line {node.lineno})")
    assert not missing, "page reads attributes the model does not have: " + \
        ", ".join(missing)


def test_arm_readouts_expose_what_the_detection_tab_prints():
    from sabes.beamline import build_source_chain
    from sabes.detection import DetectionSettings, geometry, readout

    settings = SetupSettings()
    chain = build_source_chain(settings)
    geom = geometry(chain, settings, DetectionSettings())
    read = readout(geom, chain, (14.0, 13.0), settings, DetectionSettings())
    for attribute in ("power_w", "spot_radius_m", "power_density_w_per_cm2",
                      "residual_pump_w"):
        assert hasattr(read.arms[0], attribute), attribute
    for attribute in ("twin_radius_at_dmirror_m", "pump_separation_m",
                      "separation_margin", "optics_transmission"):
        assert hasattr(geom, attribute), attribute
