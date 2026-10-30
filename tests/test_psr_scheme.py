"""The PSR scheme exposes conditional Stokes noise without changing atomic data."""
from copy import deepcopy
import json

import matplotlib
import numpy as np
import pytest

matplotlib.use("Agg")

from gabes import psr as atomic_psr, schemes
from gabes.schemes import psr as psr_scheme
from gabes_ui.export import build_export_payload, collect_view_plot_payloads, to_json_primitive


PHASE_VIEW = "Quadrature phase scan"
OPTICAL_VIEW = "Optical detuning"


@pytest.fixture(scope="module")
def atomic_case():
    """Three actual Doppler-averaged atomic points, without a full UI scan."""
    return atomic_psr.compute_psr_response([-1000.0, -950.0, -900.0])


def parameters(**changes):
    values = psr_scheme.PSRScheme().defaults()
    values.update(detuning_mhz=-950.0, view=PHASE_VIEW)
    values.update(changes)
    return values


def direct_response(raw, selected):
    """Independent backend call at the selected frequency, with frozen inputs."""
    c = raw["conditions"]
    return atomic_psr.compute_psr_response(
        [selected], isotope=c["isotope"], line=c["line"], intensity_mw_cm2=c["intensity_mw_cm2"],
        temp_c=c["temp_c"], cell_length_m=c["cell_length_m"],
        ground_relaxation_hz=c["ground_relaxation_hz"], doppler=c["doppler"],
        velocity_order=c["velocity_order"], ellipticity_step_deg=c["ellipticity_step_deg"],
        density_m3=c["density_m3"])


def test_psr_is_a_separate_registered_scheme_with_readout_only_navigation():
    scheme = psr_scheme.PSRScheme()
    assert isinstance(schemes.get("psr"), psr_scheme.PSRScheme)
    assert "psr" in [item.name for item in schemes.all_schemes()]
    assert schemes.get("psr") is not schemes.get("fwm")
    assert scheme.supports_headless_observables
    specs = {item.name: item for item in scheme.param_schema()}
    assert set(specs["view"].choices) == {PHASE_VIEW, OPTICAL_VIEW}
    readout = {"view", "detuning_mhz", "after_cell_eff_pct", "qe_pct", "excess_noise_per_od"}
    assert readout.isdisjoint(scheme.recompute_keys())
    assert all(not specs[name].recompute for name in readout)
    assert {"isotope", "line", "intensity_mw_cm2", "temp_c", "cell_mm",
            "relaxation_hz", "doppler", "resolution", "finite_difference_deg"}.issubset(scheme.recompute_keys())
    assert scheme.defaults()["excess_noise_per_od"] == 0.0


def test_navigation_does_not_enter_atomic_compute_calls(monkeypatch, atomic_case):
    calls = []

    def capture(*args, **kwargs):
        calls.append((args, kwargs))
        return deepcopy(atomic_case)

    monkeypatch.setattr(atomic_psr, "compute_psr_response", capture)
    # Cover either an imported function or a module-qualified backend call.
    monkeypatch.setattr(psr_scheme, "compute_psr_response", capture, raising=False)
    scheme = psr_scheme.PSRScheme()
    scheme.compute(parameters())
    scheme.compute(parameters(view=OPTICAL_VIEW, detuning_mhz=-925.0,
                              after_cell_eff_pct=60.0, qe_pct=80.0, excess_noise_per_od=2.0))
    assert len(calls) == 2
    np.testing.assert_array_equal(calls[0][0][0], calls[1][0][0])
    assert calls[0][1] == calls[1][1]


def test_headless_readout_is_pure_and_has_two_explicit_noise_heroes(monkeypatch, atomic_case):
    import matplotlib.pyplot as plt

    before = json.dumps(to_json_primitive(atomic_case), sort_keys=True)
    monkeypatch.setattr(plt, "subplots", lambda *a, **kw: pytest.fail("headless readout made a figure"))
    view = psr_scheme.PSRScheme().headless_observables(atomic_case, parameters())
    assert view["figure"] is None
    point = view["psr_readout"]["selected_point"]
    assert point["available"] and point["approximation_valid"]
    assert len(view["metrics"]) >= 2
    assert all(metric["tier"] == "hero" for metric in view["metrics"][:2])
    assert all(metric.get("tier") != "hero" for metric in view["metrics"][2:])
    assert view["metrics"][0]["value"] == f"{point['minimum_db']:.2f} dB"
    assert view["metrics"][1]["value"] == f"{point['maximum_db']:.2f} dB"
    assert json.dumps(to_json_primitive(atomic_case), sort_keys=True) == before


def test_selected_point_solves_atomic_shear_and_od_without_extrapolation(atomic_case):
    scheme = psr_scheme.PSRScheme()
    selected = -925.0
    view = scheme.headless_observables(atomic_case, parameters(detuning_mhz=selected))
    readout = view["psr_readout"]
    point = readout["selected_point"]
    assert point["available"]
    assert point["detuning_mhz"] == selected
    direct = direct_response(atomic_case, selected)
    for name in ("shear", "optical_depth"):
        assert point[name] == pytest.approx(direct[name][0], rel=1e-12)
    assert point["transmission"] == pytest.approx(np.exp(-point["optical_depth"]))
    assert "direct atomic solve" in point["selection_method"]
    assert point["atomic_diagnostics"]["health"]["virtual_excitation_bound"] == pytest.approx(
        direct["atomic_health"]["virtual_excitation_bound"][0])
    outside = scheme.headless_observables(atomic_case, parameters(detuning_mhz=-899.0))
    point = outside["psr_readout"]["selected_point"]
    assert not point["available"]
    assert all(item["value"] == "unavailable" for item in outside["metrics"][:2])


@pytest.mark.parametrize("resolution", ["Quick", "Balanced", "Fine"])
def test_doppler_free_selected_point_resolves_the_off_grid_resonance(resolution):
    from gabes.psr_noise import psr_gaussian_readout

    scheme = psr_scheme.PSRScheme()
    params = parameters(isotope="Rb85", doppler=False, temp_c=20., cell_mm=100.,
                        intensity_mw_cm2=.1, relaxation_hz=1000., detuning_mhz=200.,
                        resolution=resolution)
    raw = scheme.compute(params)
    direct = direct_response(raw, 200.)
    expected = psr_gaussian_readout(direct["shear"][0], direct["optical_depth"][0],
                                  detection_efficiency=.9025)
    point = scheme.headless_observables(raw, params)["psr_readout"]["selected_point"]
    assert point["available"] and point["approximation_valid"]
    assert point["shear"] < 0  # Quick interpolation formerly reversed this sign.
    assert point["shear"] == pytest.approx(-.3509031433, abs=1e-9)
    assert point["minimum_db"] == pytest.approx(-1.237065736, abs=1e-9)
    assert point["min_angle_deg"] == pytest.approx(49.85245354, abs=1e-7)
    for key in ("minimum_db", "maximum_db", "min_angle_deg"):
        assert point[key] == pytest.approx(expected[key], rel=1e-12)


def test_selected_atomic_cache_ignores_readout_controls_and_uses_raw_conditions(monkeypatch, atomic_case):
    psr_scheme._cached_selected_response.cache_clear()
    calls = []
    original = psr_scheme.compute_psr_response

    def capture(*args, **kwargs):
        calls.append((args, kwargs))
        return original(*args, **kwargs)

    monkeypatch.setattr(psr_scheme, "compute_psr_response", capture)
    scheme = psr_scheme.PSRScheme()
    try:
        p = parameters(detuning_mhz=-925., isotope="Rb85", temp_c=99., intensity_mw_cm2=8.)
        first = scheme.headless_observables(atomic_case, p)
        assert len(calls) == 1
        assert calls[0][1]["isotope"] == atomic_case["conditions"]["isotope"]
        assert calls[0][1]["temp_c"] == atomic_case["conditions"]["temp_c"]
        assert calls[0][1]["density_m3"] == atomic_case["conditions"]["density_m3"]
        point = first["psr_readout"]["selected_point"]
        health = point["atomic_diagnostics"]["health"].copy()
        point["atomic_diagnostics"]["health"]["trace_error"] = 123.
        point["atomic_diagnostics"]["conditions"]["temp_c"] = -100.
        again = scheme.headless_observables(atomic_case, dict(p, qe_pct=40., after_cell_eff_pct=60.,
                                                            excess_noise_per_od=2., view=OPTICAL_VIEW))
        assert len(calls) == 1
        diagnostics = again["psr_readout"]["selected_point"]["atomic_diagnostics"]
        assert diagnostics["health"] == health
        assert diagnostics["conditions"]["temp_c"] == atomic_case["conditions"]["temp_c"]
        scheme.headless_observables(atomic_case, dict(p, detuning_mhz=-975.))
        assert len(calls) == 2
        changed = deepcopy(atomic_case)
        changed["conditions"]["ground_relaxation_hz"] *= 2
        scheme.headless_observables(changed, p)
        assert len(calls) == 3
        assert calls[-1][1]["ground_relaxation_hz"] == changed["conditions"]["ground_relaxation_hz"]
    finally:
        psr_scheme._cached_selected_response.cache_clear()


@pytest.mark.parametrize("mode", [PHASE_VIEW, OPTICAL_VIEW])
def test_both_plot_modes_match_the_headless_noise_data(mode, atomic_case):
    import matplotlib.pyplot as plt

    scheme = psr_scheme.PSRScheme()
    params = parameters(view=mode)
    headless = scheme.headless_observables(atomic_case, params)
    view = scheme.observables(atomic_case, params)
    try:
        assert view["metrics"] == headless["metrics"]
        axis = view["figure"].axes[0]
        readout = view["psr_readout"]
        if mode == PHASE_VIEW:
            scan = readout["phase_scan"]
            assert "phase" in axis.get_xlabel().lower()
            for index, key in enumerate(("noise_a_db", "noise_b_db")):
                np.testing.assert_array_equal(axis.lines[index].get_xdata(), scan["phase_degrees"])
                np.testing.assert_array_equal(axis.lines[index].get_ydata(), scan[key])
            assert any(np.all(np.asarray(line.get_ydata()) == 0.0) for line in axis.lines[2:])
        else:
            spectrum = readout["spectrum"]
            assert "detuning" in axis.get_xlabel().lower()
            for index, key in enumerate(("minimum_db", "maximum_db")):
                np.testing.assert_array_equal(axis.lines[index].get_xdata(), spectrum["detunings_mhz"])
                np.testing.assert_allclose(axis.lines[index].get_ydata(), spectrum[key], equal_nan=True)
        view["figure"].canvas.draw()
    finally:
        plt.close(view["figure"])


def test_detection_loss_and_assumed_excess_noise_update_readout_only(atomic_case):
    scheme = psr_scheme.PSRScheme()
    original = deepcopy(atomic_case)
    clean = scheme.headless_observables(atomic_case, parameters(after_cell_eff_pct=100.0, qe_pct=100.0))
    lossy = scheme.headless_observables(atomic_case, parameters(after_cell_eff_pct=50.0, qe_pct=80.0))
    noisy = scheme.headless_observables(atomic_case, parameters(after_cell_eff_pct=100.0, qe_pct=100.0,
                                                             excess_noise_per_od=2.0))
    clean_covariance = clean["psr_readout"]["phase_scan"]["covariance"]
    np.testing.assert_allclose(lossy["psr_readout"]["phase_scan"]["covariance"],
                               0.4 * clean_covariance + 0.6 * np.eye(2), atol=1e-12)
    assert noisy["psr_readout"]["selected_point"]["minimum_db"] > clean["psr_readout"]["selected_point"]["minimum_db"]
    for key in ("shear", "optical_depth", "transmission", "g_per_m", "alpha_per_m"):
        np.testing.assert_array_equal(atomic_case[key], original[key])


def test_zero_density_preserves_vacuum_normalization():
    raw = atomic_psr.compute_psr_response([-1000.0, -950.0, -900.0], density_m3=0.0, doppler=False)
    view = psr_scheme.PSRScheme().headless_observables(raw, parameters())
    point = view["psr_readout"]["selected_point"]
    assert point["available"]
    assert point["minimum_db"] == pytest.approx(0.0, abs=1e-12)
    assert point["maximum_db"] == pytest.approx(0.0, abs=1e-12)
    np.testing.assert_allclose(view["psr_readout"]["phase_scan"]["covariance"], np.eye(2))


def test_zero_detection_efficiency_does_not_claim_a_measurable_stokes_noise(atomic_case):
    view = psr_scheme.PSRScheme().headless_observables(atomic_case, parameters(qe_pct=0.0))
    point = view["psr_readout"]["selected_point"]
    assert not point["available"]
    assert point["status"]
    assert all(item["value"] == "unavailable" for item in view["metrics"][:2])


def test_absent_input_carrier_cannot_define_a_polarization_shot_noise_reference():
    raw = atomic_psr.compute_psr_response([-1000.0, -950.0, -900.0],
                                          intensity_mw_cm2=0.0, doppler=False)
    view = psr_scheme.PSRScheme().headless_observables(raw, parameters(intensity_mw_cm2=0.0))
    assert not view["psr_readout"]["selected_point"]["available"]
    assert all(item["value"] == "unavailable" for item in view["metrics"][:2])


def test_applicability_warning_stays_local_and_is_not_noise_validation(atomic_case):
    raw = deepcopy(atomic_case)
    raw["approximation_valid"][1] = False
    raw["conditional_model_applicable"][1] = False
    raw["approximation_reasons"][1] = ["uniform-pump screen exceeded"]
    raw["approximation_status"][1] = "outside declared approximation screens"
    scheme = psr_scheme.PSRScheme()
    view = scheme.headless_observables(raw, parameters())
    point = view["psr_readout"]["selected_point"]
    assert point["available"]  # An explicitly conditional estimate remains useful.
    assert not point["approximation_valid"]
    assert "uniform-pump screen exceeded" in point["approximation_reasons"]
    table = "\n".join(item["markdown"] for item in view["tables"])
    assert "uniform-pump screen exceeded" in table
    unaffected = scheme.headless_observables(raw, parameters(detuning_mhz=-1000.0))
    assert unaffected["psr_readout"]["selected_point"]["approximation_valid"]


def test_isolated_outside_sample_shades_both_interpolation_intervals(atomic_case):
    import matplotlib.pyplot as plt

    raw = deepcopy(atomic_case)
    raw["approximation_valid"][1] = False
    raw["approximation_reasons"][1] = ["uniform-pump screen exceeded"]
    view = psr_scheme.PSRScheme().observables(raw, parameters(view=OPTICAL_VIEW))
    try:
        shading = view["figure"].axes[0].collections[0]
        vertices = np.concatenate([path.vertices for path in shading.get_paths()])
        assert vertices[:, 0].min() == raw["detunings_mhz"][0]
        assert vertices[:, 0].max() == raw["detunings_mhz"][-1]
        assert vertices[:, 1].max() == 1.0
    finally:
        plt.close(view["figure"])


def test_invalid_neighbors_leave_gaps_but_do_not_poison_a_valid_direct_point(atomic_case):
    scheme = psr_scheme.PSRScheme()
    for key, reason_key in (("numerical_valid", "numerical_reasons"),
                            ("passive_noise_valid", "passive_noise_reasons")):
        raw = deepcopy(atomic_case)
        raw[key][1] = False
        raw[reason_key][1] = ["synthetic invalid atomic point"]
        view = scheme.headless_observables(raw, parameters(detuning_mhz=-925.0))
        assert view["psr_readout"]["selected_point"]["available"]
        assert view["psr_readout"]["selected_point"]["approximation_valid"]
        assert not view["psr_readout"]["spectrum"]["available"][1]
        assert np.isnan(view["psr_readout"]["spectrum"]["minimum_db"][1])
        assert all(item["value"] != "unavailable" for item in view["metrics"][:2])


def test_failed_direct_point_is_unavailable_despite_valid_scan_neighbors(monkeypatch, atomic_case):
    import matplotlib.pyplot as plt

    selected = -925.
    failed = direct_response(atomic_case, selected)
    failed["numerical_valid"][:] = False
    failed["passive_noise_valid"][:] = False
    failed["approximation_valid"][:] = False
    failed["numerical_reasons"] = [["atomic solve failed: LinAlgError"]]
    failed["passive_noise_reasons"] = [["numerical response unavailable"]]
    for key in ("shear", "optical_depth", "transmission"):
        failed[key][:] = np.nan
    psr_scheme._cached_selected_response.cache_clear()
    monkeypatch.setattr(psr_scheme, "compute_psr_response", lambda *a, **kw: deepcopy(failed))
    view = None
    try:
        view = psr_scheme.PSRScheme().observables(atomic_case, parameters(detuning_mhz=selected, view=OPTICAL_VIEW))
        point = view["psr_readout"]["selected_point"]
        assert not point["available"]
        assert not point["approximation_valid"]
        assert "atomic solve failed" in point["status"]
        assert point["shear"] is None and point["optical_depth"] is None
        assert np.all(view["psr_readout"]["spectrum"]["available"])
        for line in view["figure"].axes[0].lines[:2]:
            plotted = np.asarray(line.get_ydata())[np.asarray(line.get_xdata()) == selected]
            assert plotted.size == 1 and np.isnan(plotted[0])
    finally:
        if view is not None:
            plt.close(view["figure"])
        psr_scheme._cached_selected_response.cache_clear()


@pytest.mark.parametrize("mode", [PHASE_VIEW, OPTICAL_VIEW])
def test_failed_nan_atomic_point_has_unavailable_values_and_correct_failure_reason(mode, atomic_case):
    import matplotlib.pyplot as plt

    raw = deepcopy(atomic_case)
    for key in ("shear", "optical_depth", "transmission", "g_per_m", "alpha_per_m"):
        raw[key][1] = np.nan
    for key in ("numerical_valid", "passive_noise_valid", "approximation_valid"):
        raw[key][1] = False
    raw["numerical_reasons"][1] = ["atomic solve failed: LinAlgError"]
    raw["passive_noise_reasons"][1] = ["numerical response unavailable"]
    view = psr_scheme.PSRScheme().observables(raw, parameters(view=mode))
    try:
        point = view["psr_readout"]["selected_point"]
        assert not point["available"]
        assert point["shear"] is None
        assert point["optical_depth"] is None
        assert point["transmission"] is None
        assert "atomic solve failed" in point["status"]
        assert "no detected bright carrier" not in point["status"]
        shown = [str(item["value"]) for item in view["metrics"]]
        shown.extend(item["markdown"] for item in view["tables"])
        shown.extend(text.get_text() for axis in view["figure"].axes for text in axis.texts)
        assert all("nan" not in text.lower() for text in shown)
        assert all("no detected bright carrier" not in text.lower() for text in shown)
        metrics = {item["label"]: item["value"] for item in view["metrics"]}
        assert metrics["Cell transmission"] == "unavailable"
        assert metrics["PSR shear κ"] == "unavailable"
        spectrum = view["psr_readout"]["spectrum"]
        assert not spectrum["available"][1]
        assert np.isnan(spectrum["minimum_db"][1])
        view["figure"].canvas.draw()
    finally:
        plt.close(view["figure"])


@pytest.mark.parametrize("selected", [-975., -925.])
def test_off_grid_optical_plot_evaluates_selected_coefficients_and_preserves_raw_spectrum(selected, atomic_case):
    import matplotlib.pyplot as plt
    from gabes.psr_noise import psr_gaussian_readout

    original = json.dumps(to_json_primitive(atomic_case), sort_keys=True)
    params = parameters(view=OPTICAL_VIEW, detuning_mhz=selected)
    view = psr_scheme.PSRScheme().observables(atomic_case, params)
    try:
        readout = view["psr_readout"]
        point, spectrum = readout["selected_point"], readout["spectrum"]
        raw_axis = atomic_case["detunings_mhz"]
        direct = direct_response(atomic_case, selected)
        shear, depth = direct["shear"][0], direct["optical_depth"][0]
        expected = psr_gaussian_readout(
            shear, depth, excess_noise_per_od=params["excess_noise_per_od"],
            detection_efficiency=params["after_cell_eff_pct"]*params["qe_pct"]/10000.)
        assert point["available"]
        np.testing.assert_array_equal(spectrum["detunings_mhz"], raw_axis)
        assert spectrum["minimum_db"].shape == raw_axis.shape
        for index, key in enumerate(("minimum_db", "maximum_db")):
            line = view["figure"].axes[0].lines[index]
            x, y = np.asarray(line.get_xdata()), np.asarray(line.get_ydata())
            selected_indices = np.flatnonzero(x == selected)
            assert len(x) == len(raw_axis)+1
            assert selected_indices.size == 1
            assert y[selected_indices[0]] == pytest.approx(expected[key])
            assert point[key] == pytest.approx(expected[key])
            assert view["metrics"][index]["value"] == f"{y[selected_indices[0]]:.2f} dB"
            np.testing.assert_array_equal(np.delete(x, selected_indices), raw_axis)
            np.testing.assert_allclose(np.delete(y, selected_indices), spectrum[key])
        # A direct solve replaces interpolation, including interpolation of
        # linear coefficients and of the sampled dB trace.
        assert abs(point["minimum_db"]-np.interp(selected, raw_axis, spectrum["minimum_db"])) > 1e-6
        plots = collect_view_plot_payloads(view, view["figure"])
        assert selected in plots[0]["traces"][0]["x_data"]
        assert json.dumps(to_json_primitive(atomic_case), sort_keys=True) == original
    finally:
        plt.close(view["figure"])


def test_export_retains_classical_response_noise_assumptions_and_scan(atomic_case):
    import matplotlib.pyplot as plt

    scheme = psr_scheme.PSRScheme()
    params = parameters()
    view = scheme.observables(atomic_case, params)
    try:
        plots = collect_view_plot_payloads(view, view["figure"])
        payload = build_export_payload(
            scheme=scheme, raw=atomic_case, params=params, view=view, plot_payloads=plots,
            comparison_payload=None, cache_version=scheme.cache_version, readout_cache_version="psr-test")
        decoded = json.loads(json.dumps(payload, allow_nan=False))
        assert decoded["scheme"]["name"] == "psr"
        data = decoded["view"]["psr_readout"]
        assert data["selected_point"]["available"]
        np.testing.assert_array_equal(data["spectrum"]["detunings_mhz"], atomic_case["detunings_mhz"])
        assert data["phase_scan"]["scope"]
        assert "SQL=I" in data["phase_scan"]["convention"]
        assert data["phase_scan"]["excess_noise_per_od"] == 0.0
        assert decoded["raw"]["model_scope"] == atomic_case["model_scope"]
        assert decoded["plots"][0]["traces"][0]["x_unit"] == "deg"
    finally:
        plt.close(view["figure"])
