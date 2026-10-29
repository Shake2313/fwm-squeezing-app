"""Historical detuning adapter contracts; public unified routing is tested separately."""
import json
from unittest.mock import patch

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pytest

from gabes.schemes import magneto_budker as budker
from gabes.schemes.magneto import SinglePairMagnetoReference
from gabes_ui.export import collect_view_plot_payloads


def _params(**changes):
    params = {spec.name: spec.default for spec in budker.param_schema()}
    params.update(calculation=budker.BUDKER, intensity_mw_cm2=0.8, temp_c=20.0,
                  beam_diameter_mm=2.0, detection_efficiency=0.8,
                  noise_bandwidth_hz=1.0, magnetic_signal_ft=10.0)
    params.update(changes)
    return params


def _fake_backend(detunings_mhz, *, b_ut, **conditions):
    delta = np.asarray(detunings_mhz)
    fields = np.asarray(b_ut)
    slope = 1.0e6 * (1 + np.exp(-((delta - 300.0) / 200.0) ** 2))
    transmission = np.full((fields.size, delta.size), 0.96)
    return dict(
        detunings_mhz=delta, b_ut=fields, transmission=transmission,
        optical_depth=-np.log(transmission), rotation_rad=fields[:, None] * 1e-6 * slope,
        frequency_hz=377e12 if conditions["line"] == "D1" else 384e12,
        isotope=conditions["isotope"], line=conditions["line"], ground_f=conditions["ground_f"],
        density_m3=conditions["density_m3"] if conditions["density_m3"] is not None else 2e15,
        reference_centroid_mhz=123.0,
        detuning_reference="line-strength centroid of selected ground-F optical component",
        model_scope="synthetic multilevel backend for readout tests",
        diagnostics=dict(optical_excitation_status="weak excitation", excited_fraction_max=0.001,
                         solver_status="physical steady state"),
        velocity_quadrature=dict(convergence_status="not assessed; compare interval orders",
                                 interval_order=conditions["velocity_points"]),
    )


@pytest.fixture
def synthetic_backend(monkeypatch):
    calls = []

    def capture(delta, **conditions):
        calls.append((np.array(delta), conditions))
        return _fake_backend(delta, **conditions)

    monkeypatch.setattr(budker.magneto_multilevel, "multilevel_response", capture)
    return calls


def test_budker_schema_keeps_physical_controls_in_compute_and_mode_scope():
    specs = budker.param_schema()
    assert all(spec.recompute for spec in specs)
    assert all(spec.name.startswith("budker_") for spec in specs)
    assert all(spec.visible_if == {"calculation": budker.BUDKER} for spec in specs)
    defaults = _params()
    assert defaults["budker_isotope"] == "Rb85"
    assert defaults["budker_density_cm3"] == 4.5e9
    assert defaults["budker_velocity_mixing"] == "complete"


def test_unregistered_reference_retains_historical_opt_in_adapter(synthetic_backend):
    scheme = SinglePairMagnetoReference()
    params = scheme.defaults()
    assert params["calculation"] == budker.LEGACY
    params["calculation"] = budker.BUDKER
    assert all(spec.name in scheme.recompute_keys() for spec in budker.param_schema())
    assert {"beam_diameter_mm", "detection_efficiency", "noise_bandwidth_hz", "magnetic_signal_ft"}.isdisjoint(scheme.recompute_keys())
    raw = scheme.compute(params)
    assert raw["budker_scan"] is True
    assert scheme.headless_observables(raw, params)["magneto_readout"]["calculation"] == budker.BUDKER
    assert len(synthetic_backend) == 1


def test_large_declared_field_withholds_snr_but_preserves_asd(synthetic_backend):
    params = _params()
    raw = budker.compute(params)
    small = budker.observables(raw, params, include_figures=False)["magneto_readout"]
    large = budker.observables(raw, dict(params, magnetic_signal_ft=1e6),
                               include_figures=False)["magneto_readout"]
    assert large["optimum_field_asd_t_per_sqrt_hz"] == small["optimum_field_asd_t_per_sqrt_hz"]
    assert large["optimum_snr"] is None
    assert not any(large["sensitivity"]["snr_valid"])
    assert "exceeds tested" in large["sensitivity"]["signal_linearity_status"]
    assert len(synthetic_backend) == 1


def test_budker_compute_routes_units_centroid_and_five_field_stencil(synthetic_backend):
    params = _params(budker_isotope="Rb87", budker_line="D2", budker_ground_center="lower",
                     budker_bias_ut=0.01, budker_detuning_points=4)
    raw = budker.compute(params)
    delta, call = synthetic_backend[0]
    assert np.count_nonzero(delta == 0) == 1
    assert delta.size == 5
    np.testing.assert_allclose(call["b_ut"], 0.01 + 1e-5 * np.array([-1, -.5, 0, .5, 1]))
    assert call["ground_f"] == 1
    assert call["density_m3"] == 4.5e15
    assert call["cell_length_m"] == pytest.approx(0.1)
    assert call["gamma_out_hz"] == 5e4
    assert call["qwp_deg"] == 0.0
    assert call["line"] == "D2"
    assert all(raw["slope_valid"])
    expected = 1e6 * (1 + np.exp(-((delta - 300) / 200) ** 2))
    np.testing.assert_allclose(raw["slope_rad_per_t"], expected, rtol=1e-11)


def test_budker_detector_changes_reuse_raw_and_preserve_asd_enbw(synthetic_backend):
    params = _params()
    raw = budker.compute(params)
    saved = raw["slope_rad_per_t"].copy()
    with patch.object(budker.magneto_multilevel, "multilevel_response", side_effect=AssertionError("readout solved atoms")), \
         patch("numpy.linalg.solve", side_effect=AssertionError("readout solved a matrix")):
        base = budker.observables(raw, params, include_figures=False)["magneto_readout"]
        wider = budker.observables(raw, dict(params, noise_bandwidth_hz=4.0), include_figures=False)["magneto_readout"]
        lossy = budker.observables(raw, dict(params, detection_efficiency=0.2), include_figures=False)["magneto_readout"]
    assert len(synthetic_backend) == 1
    np.testing.assert_array_equal(raw["slope_rad_per_t"], saved)
    assert wider["sensitivity"]["field_asd_t_per_sqrt_hz"] == base["sensitivity"]["field_asd_t_per_sqrt_hz"]
    assert wider["optimum_snr"] == pytest.approx(base["optimum_snr"] / 2)
    assert lossy["optimum_field_asd_t_per_sqrt_hz"] == pytest.approx(2 * base["optimum_field_asd_t_per_sqrt_hz"])
    assert base["optimum_detuning_mhz"] == 300.0
    assert base["snr_improvement_over_centroid"] > 1
    assert base["reference_centroid_mhz"] == 123.0
    assert "not assessed" in base["velocity_quadrature"]["convergence_status"]
    json.dumps(base, allow_nan=False)


def test_budker_compute_ignores_detection_keys(synthetic_backend):
    budker.compute(_params())
    budker.compute(_params(beam_diameter_mm=8.0, noise_bandwidth_hz=900.0,
                           detection_efficiency=0.1, magnetic_signal_ft=100.0))
    first, second = synthetic_backend
    np.testing.assert_array_equal(first[0], second[0])
    assert first[1].keys() == second[1].keys()
    for key in first[1]:
        np.testing.assert_equal(first[1][key], second[1][key])
    assert "beam_diameter_mm" not in first[1]


def test_budker_unresolved_derivative_has_no_sensitivity_or_displayed_slope(monkeypatch):
    def nonlinear(delta, **kw):
        result = _fake_backend(delta, **kw)
        field_t = np.asarray(kw["b_ut"]) * 1e-6
        h = np.max(abs(field_t))
        result["rotation_rad"] *= 1 + (field_t[:, None] / h) ** 2
        return result

    monkeypatch.setattr(budker.magneto_multilevel, "multilevel_response", nonlinear)
    params = _params()
    raw = budker.compute(params)
    assert not np.any(raw["slope_valid"])
    view = budker.observables(raw, params)
    try:
        data = view["magneto_readout"]
        assert data["optimum_index"] is None
        assert data["slope_rad_per_t"] == [None] * len(raw["detunings_mhz"])
        assert data["sensitivity"]["slope_rad_per_t"] == data["slope_rad_per_t"]
        assert all(metric.get("kind") == "status" for metric in view["metrics"] if metric.get("tier") == "hero")
        assert np.all(np.isnan(view["figure"].axes[0].lines[0].get_ydata()))
        json.dumps(data, allow_nan=False)
    finally:
        plt.close(view["figure"])


@pytest.mark.parametrize("shape, expected", [("center", 0.0), ("edge", 1200.0), ("zero", None)])
def test_budker_optimum_is_not_forced_off_resonance_and_edges_are_flagged(monkeypatch, shape, expected):
    def changed(delta, **kw):
        result = _fake_backend(delta, **kw)
        delta = np.asarray(delta)
        slopes = (1e6 * (1 + np.exp(-(delta / 200) ** 2)) if shape == "center" else
                  1e6 * (2 + delta / 1200) if shape == "edge" else np.zeros(delta.shape))
        result["rotation_rad"] = np.asarray(kw["b_ut"])[:, None] * 1e-6 * slopes
        return result

    monkeypatch.setattr(budker.magneto_multilevel, "multilevel_response", changed)
    params = _params()
    data = budker.observables(budker.compute(params), params, include_figures=False)["magneto_readout"]
    assert data["optimum_detuning_mhz"] == expected
    assert data["optimum_at_scan_edge"] == (shape == "edge")
    if shape == "center":
        assert data["snr_improvement_over_centroid"] == 1.0
    if shape == "edge":
        assert "expand optical scan" in data["optimum_status"]
    json.dumps(data, allow_nan=False)


def test_budker_figures_and_export_use_the_same_physical_curves(synthetic_backend):
    params = _params()
    raw = budker.compute(params)
    view = budker.observables(raw, params)
    try:
        headless = budker.observables(raw, params, include_figures=False)
        assert headless["magneto_readout"] == view["magneto_readout"]
        assert headless["figure"] is None
        assert len(view["figure"].axes) == 3
        view["figure"].canvas.draw()
        payloads = collect_view_plot_payloads(view, view["figure"])
        curves = [curve for payload in payloads for curve in payload["traces"]
                  if curve.get("quantity_metadata")]
        assert len(curves) == 3
        assert all(curve["quantity_metadata"]["detuning_reference"] == raw["backend"]["detuning_reference"]
                   for curve in curves)
        np.testing.assert_allclose(view["figure"].axes[2].lines[0].get_ydata(),
                                   np.asarray(view["magneto_readout"]["sensitivity"]["field_asd_t_per_sqrt_hz"]) * 1e15)
    finally:
        plt.close(view["figure"])


def test_budker_high_od_remains_conditional_and_is_flagged(synthetic_backend):
    params = _params()
    raw = budker.compute(params)
    raw["backend"]["optical_depth"][:] = 0.5
    raw["backend"]["transmission"][:] = np.exp(-0.5)
    data = budker.observables(raw, params, include_figures=False)["magneto_readout"]
    assert data["optimum_index"] is not None
    assert data["max_optical_depth"] == 0.5
    assert "high-OD quantitative use unsupported" in data["propagation_status"]


def test_budker_vapor_mode_does_not_apply_manual_density(synthetic_backend):
    raw = budker.compute(_params(budker_density_mode="vapor", budker_density_cm3=9e10))
    assert synthetic_backend[0][1]["density_m3"] is None
    assert raw["backend"]["density_m3"] == 2e15


def test_budker_failed_atomic_state_does_not_produce_sensitivity(monkeypatch):
    def failed(delta, **kw):
        result = _fake_backend(delta, **kw)
        result["diagnostics"]["solver_status"] = "steady-state validation failed"
        return result

    monkeypatch.setattr(budker.magneto_multilevel, "multilevel_response", failed)
    params = _params()
    data = budker.observables(budker.compute(params), params, include_figures=False)["magneto_readout"]
    assert data["optimum_index"] is None
    assert not any(data["sensitivity"]["valid"])
    assert data["diagnostics"]["solver_status"] == "steady-state validation failed"


def test_budker_minimal_real_backend_integration():
    params = _params(budker_isotope="Rb87", budker_doppler="off",
                     budker_velocity_points=2, budker_detuning_points=3,
                     budker_detuning_max_mhz=100.0, intensity_mw_cm2=0.01)
    raw = budker.compute(params)
    data = budker.observables(raw, params, include_figures=False)["magneto_readout"]
    assert len(data["transmission"]) == 3
    assert all(0 < value <= 1 for value in data["transmission"])
    assert data["reference_centroid_mhz"] == raw["backend"]["reference_centroid_mhz"]
    assert data["conditions"]["qwp_deg"] == 0.0
    json.dumps(data, allow_nan=False)
