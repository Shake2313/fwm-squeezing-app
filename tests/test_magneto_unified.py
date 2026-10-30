"""The public magneto scheme has one atomic model for every axis and preset."""
import json
from unittest.mock import patch

import numpy as np
import pytest

from gabes.schemes.magneto_unified import (
    UnifiedMagnetoScheme, FIELD_SCAN, OPTICAL_SCAN, TRANSMISSION, ROTATION,
    PARAFFIN, BUFFER,
)
from gabes.schemes import magneto_unified as unified
from gabes_ui.export import collect_view_plot_payloads


def fake_response(delta, *, b_ut, **conditions):
    delta, b = np.asarray(delta, float), np.asarray(b_ut, float)
    slope = 1e7*(1+np.exp(-((delta-300)/200)**2))
    shape = (b.size, delta.size)
    T = .9 + b[:, None]*1e-6*1000 + np.zeros(shape)
    linear = conditions["qwp_deg"] == 0
    theta = b[:, None]*1e-6*slope if linear else np.full(shape, np.nan)
    return dict(
        transmission=T, rotation_rad=theta, optical_depth=-np.log(T),
        rotation_slope_rad_per_t=np.broadcast_to(slope, shape) if linear else np.full(shape, np.nan),
        transmission_slope_per_t=np.full(shape, 1000.), slope_valid=np.ones(shape, bool),
        slope_relative_residual=np.full(shape, 1e-14),
        frequency_hz=377e12, reference_centroid_mhz=-1000.,
        detuning_reference="line-strength centroid of selected ground-F optical component",
        isotope=conditions["isotope"], line=conditions["line"], ground_f=conditions["ground_f"],
        density_m3=conditions["density_m3"] or 2e15,
        model_scope="synthetic multihyperfine test model", rotation_status=(
            "linear-input circular-birefringence approximation" if linear else "unsupported polarization"),
        diagnostics=dict(solver_status="physical steady state", optical_excitation_status="weak excitation",
                         magnetic_response_status="analytic DC response"),
        velocity_quadrature=dict(convergence_status="not assessed; compare interval orders"),
    )


@pytest.fixture
def backend(monkeypatch):
    calls = []

    def capture(delta, **conditions):
        calls.append((np.array(delta), conditions))
        return fake_response(delta, **conditions)

    monkeypatch.setattr(unified.magneto_multilevel, "multilevel_response", capture)
    return calls


def test_unified_schema_has_no_model_switch_and_only_physical_condition_presets():
    sc = UnifiedMagnetoScheme()
    specs = {s.name: s for s in sc.param_schema()}
    assert "calculation" not in specs and not any(name.startswith("budker_") for name in specs)
    assert sc.defaults()["scan_axis"] == FIELD_SCAN
    assert sc.defaults()["isotope"] == "Rb85"
    assert specs["scan_axis"].choices == (FIELD_SCAN, OPTICAL_SCAN)
    assert {"signal_type", "detection_efficiency", "beam_diameter_mm", "noise_bandwidth_hz", "show_shot_noise"}.isdisjoint(sc.recompute_keys())
    assert {"cell_mm", "density_mode", "density_cm3"}.isdisjoint(sc.recompute_keys())
    assert {"isotope", "line", "qwp_deg", "b_offset_ut", "laser_detuning_mhz"}.issubset(sc.recompute_keys())
    for preset in sc.presets():
        assert not {"calculation", "model"}.intersection(preset.values)
        assert preset.values.get("scan_axis", FIELD_SCAN) == FIELD_SCAN


def test_regime_defaults_are_sidebar_snapshots_without_topbar_choices(monkeypatch):
    from gabes_ui import defaults
    sc = UnifiedMagnetoScheme()
    records = defaults.builtin_defaults(sc)
    by_name = {record["name"]: record for record in records.values()}
    expected = {"Hanle", "NMOR"}
    assert set(by_name) == expected
    assert [record["name"] for record in records.values()] == ["Hanle", "NMOR"]
    assert sc.recommended_defaults(sc.defaults()) is None
    assert not any(spec.applies_defaults for spec in sc.param_schema())
    state = {}
    monkeypatch.setattr(defaults.st, "session_state", state)
    for name in expected:
        # A prior optical/circular/buffer selection must not leak into the
        # complete snapshot loaded by the sidebar dropdown.
        state.update(magneto__scan_axis=OPTICAL_SCAN, magneto__qwp_deg=17.,
                     magneto__cell_type=BUFFER, magneto__laser_detuning_mhz=123.)
        defaults._apply(sc, by_name[name])
        loaded = defaults._params(sc)
        assert loaded == by_name[name]["params"]
        assert loaded["scan_axis"] == FIELD_SCAN
        assert loaded["signal_type"] == (ROTATION if name == "NMOR" else TRANSMISSION)
        assert loaded["cell_type"] == (BUFFER if name == "Hanle" else PARAFFIN)
        assert loaded["qwp_deg"] == 0.


def test_both_readout_defaults_use_the_unified_backend(backend):
    sc = UnifiedMagnetoScheme()
    names = {"Hanle", "NMOR"}
    for preset in sc.presets():
        if preset.name not in names:
            continue
        p = dict(sc.defaults(), **preset.values)
        raw = sc.compute(p)
        assert raw["model"] == "Budker multihyperfine"
        assert backend[-1][1]["magnetic_derivative"] is True
    assert len(backend) == 2


@pytest.mark.parametrize("mode,signal", [(None, ROTATION), ("hanle", TRANSMISSION),
                                        ("eia", TRANSMISSION), ("nmor", ROTATION)])
def test_every_public_alias_uses_the_same_backend(mode, signal, backend):
    sc = UnifiedMagnetoScheme(mode)
    params = dict(sc.defaults(), doppler="off", scan_points=11)
    raw = sc.compute(params)
    assert raw["model"] == "Budker multihyperfine"
    assert params["signal_type"] == signal
    assert backend[-1][1]["magnetic_derivative"] is True
    assert sc.name == (mode or "magneto")


def test_both_axes_share_conditions_and_identical_operating_point(backend):
    sc = UnifiedMagnetoScheme()
    p = dict(sc.defaults(), scan_points=11, laser_detuning_mhz=600., b_offset_ut=.0002,
             residual_transverse_b_ut=.001, transverse_field_angle_deg=90.,
             gamma_out_hz=4321., gamma_in_hz=8., pure_dephasing_hz=13.)
    braw = sc.compute(p)
    bview = sc.headless_observables(braw, p)["magneto_readout"]
    op = bview["operating_index"]
    optical = dict(p, scan_axis=OPTICAL_SCAN, detuning_points=5)
    oraw = sc.compute(optical)
    oview = sc.headless_observables(oraw, optical)["magneto_readout"]
    for key in braw["conditions"]:
        assert braw["conditions"][key] == oraw["conditions"][key]
    selected = list(oview["x"]).index(600.)
    assert bview["transmission"][op] == oview["transmission"][selected]
    assert bview["rotation_slope_rad_per_t"][op] == oview["rotation_slope_rad_per_t"][selected]
    assert bview["sensitivity"]["field_asd_t_per_sqrt_hz"][op] == oview["sensitivity"]["field_asd_t_per_sqrt_hz"][selected]
    assert backend[-1][1]["b_x_ut"] == pytest.approx(0., abs=1e-18)
    assert backend[-1][1]["b_y_ut"] == pytest.approx(.001)


def test_buffer_switch_selects_single_region_without_ignored_pressure_knob(backend):
    sc = UnifiedMagnetoScheme()
    p = dict(sc.defaults(), cell_type=BUFFER, ground_relaxation_hz=4000.)
    raw = sc.compute(p)
    assert raw["conditions"]["gamma_out_hz"] == raw["conditions"]["gamma_in_hz"] == 0
    assert raw["conditions"]["ground_relaxation_hz"] == 4000.
    assert "ne_pressure_torr" not in sc.defaults()
    scope = sc.headless_observables(raw, p)["magneto_readout"]["cell_scope"]
    assert "pressure broadening" in scope and "omitted" in scope


def test_optical_axis_uses_analytic_slopes_without_finite_difference_gaps(backend):
    sc = UnifiedMagnetoScheme()
    p = dict(sc.defaults(), scan_axis=OPTICAL_SCAN, detuning_points=25, show_shot_noise=False)
    raw = sc.compute(p)
    view = sc.observables(raw, p)
    data = view["magneto_readout"]
    assert all(data["slope_valid"])
    assert all(value is not None for value in data["rotation_slope_rad_per_t"])
    assert all(value is not None for value in data["sensitivity"]["field_asd_t_per_sqrt_hz"])
    assert np.asarray(backend[0][1]["b_ut"]).size == 3
    assert "analytic" in data["slope_method"]
    for ax in view["figure"].axes:
        assert np.isfinite(ax.lines[0].get_ydata()).all()
    view["figure"].canvas.draw()
    payload = collect_view_plot_payloads(view, view["figure"])
    assert len(payload[0]["traces"]) >= 3
    json.dumps(data, allow_nan=False)
    import matplotlib.pyplot as plt
    plt.close(view["figure"])


def test_true_zero_response_remains_explicitly_unavailable(monkeypatch):
    def zero(delta, **conditions):
        result = fake_response(delta, **conditions)
        result["rotation_slope_rad_per_t"] = np.zeros_like(result["rotation_slope_rad_per_t"])
        return result
    monkeypatch.setattr(unified.magneto_multilevel, "multilevel_response", zero)
    sc = UnifiedMagnetoScheme()
    p = dict(sc.defaults(), scan_axis=OPTICAL_SCAN)
    data = sc.headless_observables(sc.compute(p), p)["magneto_readout"]
    assert data["optimum_index"] is None
    assert all(v is None for v in data["sensitivity"]["field_asd_t_per_sqrt_hz"])
    assert all("zero magnetic slope" in value for value in data["sensitivity"]["status"])
    assert all(value == 0 for value in data["rotation_slope_rad_per_t"])


def test_circular_input_retains_transmission_response_and_suppresses_rotation(backend):
    sc = UnifiedMagnetoScheme()
    p = dict(sc.defaults(), qwp_deg=45., signal_type=TRANSMISSION, scan_points=11)
    raw = sc.compute(p)
    data = sc.headless_observables(raw, p)["magneto_readout"]
    assert all(data["sensitivity"]["valid"])
    assert not any(data["rotation_sensitivity"]["valid"])
    assert all(value is not None for value in data["transmission_slope_per_t"])
    assert len(backend) == 1


def test_detector_changes_reuse_the_atomic_state(backend):
    sc = UnifiedMagnetoScheme()
    p = dict(sc.defaults(), scan_axis=OPTICAL_SCAN)
    raw = sc.compute(p)
    before = sc.headless_observables(raw, p)["magneto_readout"]
    with patch.object(unified.magneto_multilevel, "multilevel_response", side_effect=AssertionError("readout solved")):
        after = sc.headless_observables(raw, dict(p, detection_efficiency=.2, noise_bandwidth_hz=4.))["magneto_readout"]
    assert after["optimum_field_asd_t_per_sqrt_hz"] == pytest.approx(2*before["optimum_field_asd_t_per_sqrt_hz"])
    assert after["optimum_snr"] == pytest.approx(before["optimum_snr"]/4)


def test_large_signal_and_local_nonlinearity_never_mask_analytic_asd(backend, monkeypatch):
    sc = UnifiedMagnetoScheme()
    p = dict(sc.defaults(), scan_axis=OPTICAL_SCAN)
    raw = sc.compute(p)
    ordinary = sc.headless_observables(raw, p)["magneto_readout"]
    large = sc.headless_observables(raw, dict(p, magnetic_signal_ft=1e6))["magneto_readout"]
    assert large["optimum_field_asd_t_per_sqrt_hz"] == ordinary["optimum_field_asd_t_per_sqrt_hz"]
    assert large["optimum_snr"] is None
    assert not any(large["sensitivity"]["snr_valid"])
    assert "exceeds tested" in large["sensitivity"]["signal_linearity_status"]

    def curved(delta, *, b_ut, **conditions):
        result = fake_response(delta, b_ut=b_ut, **conditions)
        # Keep this deliberately outside the 5% boundary under readout scaling.
        result["rotation_rad"] += 1.1e30*(np.asarray(b_ut)[:, None]*1e-6)**3
        return result
    monkeypatch.setattr(unified.magneto_multilevel, "multilevel_response", curved)
    curved_raw = sc.compute(p)
    curved_view = sc.headless_observables(curved_raw, p)["magneto_readout"]
    assert all(curved_view["slope_valid"])
    assert all(v is not None for v in curved_view["sensitivity"]["field_asd_t_per_sqrt_hz"])
    assert curved_view["optimum_snr"] is None


@pytest.mark.parametrize("axis", [FIELD_SCAN, OPTICAL_SCAN])
def test_real_small_multilevel_pipeline_for_both_axes(axis):
    sc = UnifiedMagnetoScheme()
    p = dict(sc.defaults(), isotope="Rb87", scan_axis=axis, doppler="off", scan_points=11,
             detuning_points=5, intensity_mw_cm2=.01, show_shot_noise=False)
    raw = sc.compute(p)
    view = sc.headless_observables(raw, p)
    assert raw["backend"]["diagnostics"]["solver_status"] == "physical steady state"
    assert raw["model"] == "Budker multihyperfine"
    assert np.isfinite(view["magneto_readout"]["transmission"]).all()
    assert all(view["magneto_readout"]["slope_valid"])
