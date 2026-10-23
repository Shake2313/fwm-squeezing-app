"""Physical-field readout and sampled-grid trust, independent of atomic dynamics."""
import importlib.util
from pathlib import Path
from unittest.mock import patch

import numpy as np
import pytest

from gabes import observables, schemes
from gabes.schemes import magneto
from gabes_ui.readout import partition_readout


@pytest.fixture(scope="module")
def template():
    sc = schemes.get("magneto")
    params = dict(sc.defaults(), doppler="off", scan_points=21)
    return sc, params, sc.compute(params)


def _synthetic_raw(template, offset=0.257, points=401):
    sc, params, raw = template
    params = dict(params, b_offset_ut=offset)
    command = np.linspace(-1.0, 1.0, points)
    physical = command + offset
    alpha = 5.0 - 2.0 / (1.0 + (physical / 0.2) ** 2)
    rotation = 0.003 + 0.01 * physical + 0.02 * physical ** 3
    coupling = observables.chi_phys(
        1.0, raw["N_eff"], dipole=raw["dipole"],
        line_strength=params["line_strength"])
    chi_probe = 1j * alpha / (raw["k_vec"] * coupling)
    chi_p = rotation / (0.25 * raw["k_vec"] * params["cell_mm"] * 1e-3 * coupling)
    changed = dict(raw, b_ut=command, b_physical_ut=physical,
                   b_offset_ut=offset, chi_probe=chi_probe,
                   chi_p=chi_p.astype(complex), chi_m=np.zeros_like(chi_probe))
    return sc, params, changed, alpha, rotation


def test_physical_zero_and_command_zero_are_distinct(template):
    sc, params, raw, alpha, rotation = _synthetic_raw(template)
    view = sc.headless_observables(raw, params)
    data = view["magneto_readout"]
    transmission = np.exp(-alpha * params["cell_mm"] * 1e-3)
    assert data["default_reference"] == "physical zero"
    assert data["physical_zero_command_ut"] == -0.257
    assert data["physical_zero_transmission"] == pytest.approx(
        np.interp(-0.257, raw["b_ut"], transmission))
    assert abs(data["physical_zero_transmission"] - data["command_zero_transmission"]) > .005
    assert data["physical_zero_rotation_rad"] == pytest.approx(0.003, abs=1e-7)
    assert data["command_zero_rotation_rad"] == pytest.approx(
        0.003 + .01 * .257 + .02 * .257 ** 3)
    assert data["slope_rad_per_ut"] == pytest.approx(.01, abs=1e-12)
    assert data["slope_status"] == "sample-grid consistent"
    metrics = {m["label"]: m for m in view["metrics"]}
    assert "Transmission at physical B=0" in metrics
    assert "Transmission at command B=0" not in metrics
    assert "Peak shot-noise SNR" in metrics
    derived = next(table["markdown"] for table in view["tables"]
                   if table["title"] == "Derived quantities")
    assert "Transmission at command B=0" in derived
    assert "closed hyperfine pair" in metrics["Model scope"]["value"]
    assert "repump omitted" in metrics["Model scope"]["value"]


@pytest.mark.parametrize("signal", [magneto.SIGNAL_TRANSMISSION, magneto.SIGNAL_NMOR])
def test_physical_zero_outside_scan_is_unavailable(template, signal):
    sc, params, raw, _, _ = _synthetic_raw(template, offset=3.0)
    params["signal_type"] = signal
    view = sc.headless_observables(raw, params)
    data = view["magneto_readout"]
    assert not data["physical_zero_in_scan"]
    assert data["physical_zero_transmission"] is None
    assert data["physical_zero_rotation_rad"] is None
    assert data["central_width_ut"] is None
    assert data["slope_rad_per_ut"] is None
    assert data["width_status"] == data["slope_status"] == "physical zero outside scan"
    assert data["command_zero_transmission"] is not None
    heroes = [m for m in view["metrics"] if m.get("tier") == "hero"]
    assert heroes and all(m.get("kind") == "status" or m["value"] == "unavailable"
                          for m in heroes)
    displayed_heroes, secondary, _, _ = partition_readout(view["metrics"], view["hero_count"])
    assert all(m["value"] != "unavailable" for m in displayed_heroes + secondary)


def test_central_width_uses_both_asymmetric_crossings():
    # Deliberately omit an exact center sample. The left/right half widths are
    # 0.2 and 0.5, so doubling the right half would give the wrong answer 1.0.
    x = np.linspace(-1.0, 1.0, 1000)
    signal = np.maximum(0.0, np.where(x < 0.0, 2.0 * (1.0 + x / .4),
                                    2.0 * (1.0 - x)))
    readout = magneto._central_feature(x, signal)
    assert readout["left"] == pytest.approx(-.2, abs=.002)
    assert readout["right"] == pytest.approx(.5, abs=.002)
    assert readout["width"] == pytest.approx(.7, abs=.003)


def test_local_derivative_handles_zero_between_samples():
    x = np.linspace(-1.0, 1.0, 101) + .013
    rotation = .1 + .2 * x + .5 * x ** 2 - .7 * x ** 3 + .8 * x ** 4
    assert magneto._local_derivative(x, rotation) == pytest.approx(.2, abs=1e-12)
    assert magneto._local_derivative(x, rotation, stride=2) == pytest.approx(.2, abs=1e-12)


def test_readout_never_solves_or_mutates_raw_and_render_matches(template):
    sc, params, raw, _, _ = _synthetic_raw(template)
    saved = {k: v.copy() for k, v in raw.items() if isinstance(v, np.ndarray)}
    fail = AssertionError("readout attempted an atomic or native linear solve")
    with patch.object(sc, "compute", side_effect=fail), \
         patch("gabes.core.steady_state_from_liouvillian", side_effect=fail), \
         patch("gabes.core.steady_state_batched", side_effect=fail), \
         patch("numpy.linalg.solve", side_effect=fail), \
         patch("numpy.linalg.lstsq", side_effect=fail):
        headless = sc.headless_observables(raw, params)
    with patch.object(sc, "compute", side_effect=fail):
        rendered = sc.observables(raw, params)
    assert rendered["magneto_readout"] == headless["magneto_readout"]
    for key, before in saved.items():
        np.testing.assert_array_equal(raw[key], before)
    axis = rendered["figure"].axes[0]
    assert axis.get_xlabel() == "Commanded magnetic field B  [uT]"
    references = {line.get_label(): line.get_xdata() for line in axis.lines}
    assert np.allclose(references["Physical zero"], -.257)
    import matplotlib.pyplot as plt
    plt.close(rendered["figure"])


@pytest.fixture(scope="module")
def preset_schedule():
    sc = schemes.get("magneto")
    results = {}
    for name, preset in sc.recommended_defaults(sc.defaults()).items():
        if "EIT dip" not in name and "NMOR" not in name:
            continue
        key = "NMOR" if "NMOR" in name else "EIT"
        results[key] = {}
        for count in (121, 201, 401, 801):
            params = dict(sc.defaults(), **preset)
            params["scan_points"] = count
            results[key][count] = sc.headless_observables(
                sc.compute(params), params)["magneto_readout"]
    return results


@pytest.mark.parametrize("preset", ["EIT", "NMOR"])
def test_deterministic_preset_schedule_preserves_sign_width_and_slope(preset_schedule, preset):
    schedule = preset_schedule[preset]
    for count in (201, 401, 801):
        current = schedule[count]
        assert current["central_width_ut"] > 0.0
        assert current["slope_rad_per_ut"] < 0.0
    fine, finer = schedule[401], schedule[801]
    assert fine["feature_amplitude_per_m"] * finer["feature_amplitude_per_m"] > 0.0
    for field in ("central_width_ut", "slope_rad_per_ut"):
        assert abs(fine[field] - finer[field]) / abs(finer[field]) <= .05
    assert fine["width_status"] == finer["width_status"] == "sample-grid consistent"
    assert fine["slope_status"] == finer["slope_status"] == "sample-grid consistent"


def test_default_nmor_slope_is_withheld_when_nested_stencils_disagree(preset_schedule):
    current = preset_schedule["NMOR"][121]
    assert current["slope_status"] == "nonconverged on sampled grid"
    assert current["slope_relative_change"] > .3
    sc = schemes.get("magneto")
    params = next(dict(sc.defaults(), **preset) for name, preset in
                  sc.recommended_defaults(sc.defaults()).items() if "NMOR" in name)
    view = sc.headless_observables(sc.compute(params), params)
    metrics = {m["label"]: m for m in view["metrics"]}
    assert metrics["Slope dtheta/dB"]["value"] == "unavailable"
    assert "tier" not in metrics["Slope dtheta/dB"]
    assert metrics["Slope status"]["tier"] == "hero"
    assert view["hero_count"] == 1
    assert "not a full convergence certificate" in metrics["Slope status"]["help"]


def test_ground_pure_dephasing_keeps_the_existing_api_key(template):
    sc, _, _ = template
    spec = next(s for s in sc.param_schema() if s.name == "collisional_depol_khz")
    assert spec.label == "Ground pure dephasing"
    assert spec.default == .5 and spec.unit == "kHz" and spec.recompute
    assert "without redistribution of Zeeman populations" in spec.help
    # The population-preserving channel is independently pinned by
    # test_transit_reload_and_collisional_dephasing_are_distinct_channels.


@pytest.fixture(scope="module")
def verifier():
    path = Path(__file__).with_name("verify_hanle_eit_eia.py")
    spec = importlib.util.spec_from_file_location("magneto_verifier", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_literature_width_comparisons_do_not_promote_sign_agreement(verifier, monkeypatch, capsys):
    def fake_run(qwp, **kwargs):
        return dict(amp=-1.0 if qwp == 0.0 else 1.0,
                    fwhm_ut=.095 if qwp == 0.0 else .2305)

    monkeypatch.setattr(verifier, "run", fake_run)
    monkeypatch.setattr(verifier, "transition_dependence", lambda: None)
    verifier.main()
    output = capsys.readouterr().out
    assert "Sign/trend PASS" in output
    assert output.count("Absolute-width FAIL") == 2
    assert "same sub-mG order" not in output and "MATCH" not in output
    assert "CHECK is not experimental validation" in output
    assert verifier.width_comparison(1.1, 1.0)["status"] == "CHECK"
    assert verifier.width_comparison(float("nan"), 1.0)["status"] == "CHECK"
    assert verifier.width_comparison(3.349, 2.4)["status"] == "FAIL"
