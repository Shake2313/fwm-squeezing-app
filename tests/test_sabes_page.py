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
from types import SimpleNamespace

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from gabes.schemes import fwm  # noqa: E402
from sabes import bridge  # noqa: E402
from sabes.beamline import SetupSettings, build_source_chain  # noqa: E402
from sabes.calibration import default_calibration  # noqa: E402
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


def test_eom_and_etalon_throughput_have_no_duplicate_transmission_control():
    by_id = {node.id: node for part in sabes_page.layout.LAYOUTS
             for node in part.nodes}
    assert by_id["eom"].transmission_key == ""
    for name in ("etalon_1", "etalon_2", "etalon_3"):
        assert by_id[name].transmission_key == ""
        assert by_id[name].params


def test_session_calibration_changes_only_editable_overrides(monkeypatch):
    base = default_calibration()
    editable = next(name for name in sabes_page._editable_transmission_keys()
                    if base.get(name).verified)
    state = {
        sabes_page.SESSION_PREFIX + editable: base.value(editable),
        sabes_page.SESSION_PREFIX + "t_lens_probe": 0.75,
        # Hidden duplicate from an older session must not affect the EOM model.
        sabes_page.SESSION_PREFIX + "t_eom": 0.10,
    }
    monkeypatch.setattr(sabes_page.st, "session_state", state)
    calibration = sabes_page._session_calibration()
    assert calibration.get(editable).provenance == base.get(editable).provenance
    assert calibration.value("t_lens_probe") == pytest.approx(0.75)
    assert calibration.get("t_lens_probe").provenance == "hand"
    assert calibration.value("t_eom") == base.value("t_eom")


def test_solve_knobs_uses_the_session_calibration(monkeypatch):
    sentinel = object()
    seen = []
    monkeypatch.setattr(sabes_page.st, "session_state", {})
    monkeypatch.setattr(sabes_page, "_current",
                        lambda: (SetupSettings(), DetectionSettings()))
    monkeypatch.setattr(sabes_page, "_session_calibration", lambda: sentinel)

    def solved(value):
        def solve(target, settings, calibration):
            seen.append(calibration)
            return value
        return solve

    monkeypatch.setattr(sabes_page, "solve_split_angle_deg", solved(10.0))
    monkeypatch.setattr(sabes_page, "solve_seed_polarizer_deg", solved(20.0))
    monkeypatch.setattr(sabes_page, "solve_seed_trim_deg", solved(30.0))
    sabes_page._apply_solved_knobs()
    assert seen == [sentinel, sentinel, sentinel]


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


# ---------------------------------------------------------- production canvas

def test_canvas_events_round_trip_through_the_real_table_state(monkeypatch):
    """Production optic and table clicks survive without the retired spike."""
    state = {}
    monkeypatch.setattr(sabes_page.st, "session_state", state)

    optic = {"seq": 1, "kind": sabes_page.canvas_module.KIND_OPTIC,
             "id": "cell"}
    assert sabes_page._handle_canvas_event(optic, "part2") is True
    assert state[sabes_page._key("selected")] == "cell"
    assert state[sabes_page._key("open_optic")] == "cell"
    assert sabes_page._handle_canvas_event(optic, "part2") is False

    point = {"seq": 2, "kind": sabes_page.canvas_module.KIND_POINT,
             "x": 12.5, "y": 34.0}
    assert sabes_page._handle_canvas_event(point, "part2") is True
    assert state[sabes_page._probe_key("part2")] == [12.5, 34.0]
    assert state[sabes_page._key("selected")] is None


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
                      DetectionSettings(pump_leakage_dbm=-40.0),
                      DetectionSettings(eom_unwanted_rin_db_per_hz=-80.0)):
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
    assert keys <= recompute | {"mode"} | set(bridge.SOLVE_AUDIT_KEYS)
    assert "tpd" not in keys
    for essential in ("opd", "temp_c", "pump_mw", "pump_waist_um", "qe_pct"):
        assert essential in keys


def test_cached_solve_key_preserves_eom_audit_metadata():
    params = _params()
    reconstructed = fwm.FWMScheme().defaults()
    reconstructed.update(dict(bridge.solve_key(params)))
    raw = fwm.FWMScheme().compute(reconstructed)
    ledger = raw["eom_seed_spectrum"]
    assert ledger["provenance"] == params["eom_seed_spectrum_provenance"]
    assert ledger["status"] == "unsupported"
    assert ledger["application"] == "unapplied"


def test_page_routes_global_noise_only_to_the_post_cell_balanced_analyser():
    source = (ROOT / "sabes_page.py").read_text(encoding="utf-8")
    block = source[source.index('if choice == "Spectrum analyser"'):
                   source.index("raise KeyError(choice)")]
    assert 'reading.link.beam in ("probe", "conjugate")' in block
    assert "pump_leakage_dbm=(detection.pump_leakage_dbm" in block
    assert "eom_noise=(readout.eom_noise" in block


def test_scope_scan_uses_the_selected_twin_gain_and_current_path_power(
        monkeypatch):
    chain = build_source_chain()
    raw = {
        "probe_axis_GHz": np.array([-1.0, 0.0, 1.0]),
        "G_s": np.array([50.0, 50.0, 50.0]),
        "G_c": np.array([1.0, 2.0, 3.0]),
    }
    result = SimpleNamespace(
        raw=raw, gains=(15.0, 2.0), chain=chain,
        geometry=SimpleNamespace(optics_transmission=1.0, arms=()),
    )
    reading = SimpleNamespace(
        link=SimpleNamespace(beam="conjugate"),
        beam=sabes_page._twin_states(result)["conjugate"],
    )
    detector = sabes_page.instruments.Photodiode(
        responsivity_a_per_w=0.5, transimpedance_v_per_a=2.0)
    signal = detector.convert(reading.beam)
    monkeypatch.setattr(
        sabes_page.st, "session_state",
        {sabes_page._key("scope_mode"): "Scan"})

    trace = sabes_page._scope_reading(signal, result, reading).trace
    wanted = reading.beam.power_at(-chain.seed_offset_hz)
    background = reading.beam.total_power_w - wanted
    expected_power = raw["G_c"] * wanted / 2.0 + background
    assert np.allclose(trace.series["signal"], expected_power)


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
    for attribute in ("rin_db_per_hz", "fractional_intensity_rms",
                      "unwanted_detector_power_w", "classical_rin_excess_sql"):
        assert hasattr(read.eom_noise, attribute), attribute
    for attribute in ("twin_radius_at_dmirror_m", "pump_separation_m",
                      "separation_margin", "optics_transmission"):
        assert hasattr(geom, attribute), attribute


def test_post_cell_twin_states_amplify_only_the_wanted_mode():
    from sabes.beamline import build_source_chain

    chain = build_source_chain()
    result = SimpleNamespace(
        raw={}, gains=(15.0, 14.0), chain=chain,
        geometry=SimpleNamespace(
            optics_transmission=0.85,
            arms=(SimpleNamespace(name="probe", optics_transmission=0.9),
                  SimpleNamespace(name="conjugate",
                                  optics_transmission=0.8))))
    twin = sabes_page._twin_states(result)

    assert twin["probe"].power_at(chain.seed_offset_hz) == pytest.approx(
        chain.wanted_seed_sideband_power_w * 15.0 * 0.9)
    assert twin["conjugate"].power_at(-chain.seed_offset_hz) == pytest.approx(
        chain.wanted_seed_sideband_power_w * 14.0 * 0.8)
    assert twin["conjugate"].power_at(chain.seed_offset_hz) == 0.0
    conjugate_frequency = sabes_page.instruments.Wavemeter(
        wanted_offset_hz=-chain.seed_offset_hz).measure(twin["conjugate"])
    assert conjugate_frequency.quantity(
        "Dominant-line offset").value == pytest.approx(
            -chain.seed_offset_hz / 1e9)
    assert twin["probe"].power_at(0.0) == pytest.approx(
        chain.eom_residual_carrier_power_w * 0.9)
    assert twin["conjugate"].power_at(0.0) == 0.0
