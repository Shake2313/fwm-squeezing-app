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
    known |= {"eom_offset_mhz", "etalon_detune_ghz", "resolution"}
    for group, controls in sabes_page.CONTROL_GROUPS:
        assert controls, group
        for name, label, unit, low, high, step, help_text in controls:
            assert name in known, f"{group}/{name} is not a settings field"
            assert help_text, f"{group}/{name} has no help text"
            assert low < high and step > 0


def test_defaults_cover_every_control_and_land_on_the_paper_point():
    defaults = sabes_page._defaults()
    for _group, controls in sabes_page.CONTROL_GROUPS:
        for name, *_ in controls:
            assert name in defaults, name
    # The generator control is an offset from the hyperfine splitting, and the
    # default offset is what produces the -8 MHz operating point.
    assert defaults["eom_offset_mhz"] == pytest.approx(8.0, abs=1e-6)
    assert defaults["resolution"] in sabes_page.FIDELITIES


def test_control_names_are_unique_across_groups():
    seen = []
    for _group, controls in sabes_page.CONTROL_GROUPS:
        seen += [name for name, *_ in controls]
    assert len(seen) == len(set(seen))


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
                      DetectionSettings(iris_radius_mm=3.0)):
        assert bridge.solve_key(_params(detection=detection)) == base, detection


def test_a_clipping_iris_does_change_the_solve_key():
    base = bridge.solve_key(_params())
    clipped = bridge.solve_key(
        _params(detection=DetectionSettings(iris_radius_mm=0.4)))
    assert clipped != base


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
