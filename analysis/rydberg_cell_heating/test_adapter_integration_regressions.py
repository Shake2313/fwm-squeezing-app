"""Regression tests for thermal/dephasing hand-off into report adapters."""

from __future__ import annotations

import numpy as np
import pytest
import json

from analysis.rydberg_cell_heating.adapters import (
    ElectrometryAdapter,
    StaticRydbergAdapter,
)
from gabes import species
from gabes import observables


def test_static_adapter_forwards_effective_and_cold_spot_temperatures():
    static = StaticRydbergAdapter()
    point = static.temperature_point(
        50.0,
        "EIT",
        {},
        include_spectrum=False,
        heater_setpoint_c=55.0,
        cold_spot_temp_c=40.0,
    )

    expected_density = (
        species.number_density(species.RB85, 40.0 + 273.15)
        * (40.0 + 273.15)
        / (50.0 + 273.15)
    )
    assert point.params["temperature_model"] == "Separated"
    assert point.params["heater_setpoint_c"] == 55.0
    assert point.params["effective_temp_c"] == 50.0
    assert point.params["cold_spot_temp_c"] == 40.0
    assert point.raw["N"] == pytest.approx(expected_density)


def test_configured_electrometry_rebuild_includes_density_dephasing(monkeypatch):
    static = StaticRydbergAdapter()
    point = static.temperature_point(
        50.0,
        "AT electrometry",
        {"density_dephasing_mhz_per_1e16_m3": 0.02},
        include_spectrum=False,
    )
    assert point.raw["density_dephasing_mhz"] > 0.0

    captured_ground_dephasing: list[float] = []
    original_atom = static.scheme._atom

    def recording_atom(ground_dephasing: float, rf_dephasing: float):
        captured_ground_dephasing.append(float(ground_dephasing))
        return original_atom(ground_dephasing, rf_dephasing)

    monkeypatch.setattr(static.scheme, "_atom", recording_atom)
    adapter = ElectrometryAdapter(
        static,
        {
            "enabled": True,
            "mode": "configured",
            "status": "ASSUMED",
            "transition_dipole_c_m": 2.0e-27,
            "angular_factor": 0.5,
            "field_amplitude_convention": "peak",
            "quantum_efficiency": 0.8,
            "signal_probe_power_scale": 1.0,
            "probe_detuning_search_points": 3,
        },
    )
    result = adapter.evaluate(point)

    expected_mhz = sum(
        point.raw[name]
        for name in (
            "rydberg_dephasing_mhz",
            "temperature_dephasing_mhz",
            "density_dephasing_mhz",
            "transit_mhz",
        )
    )
    assert result["available"]
    assert captured_ground_dephasing
    assert captured_ground_dephasing[0] / (2.0 * np.pi * 1.0e6) == pytest.approx(
        expected_mhz
    )


@pytest.fixture(scope="module")
def opacity_point():
    static = StaticRydbergAdapter()
    point = static.temperature_point(40.0, "AT electrometry", {}, include_spectrum=False)
    return static, point


def _configured_adapter(static):
    return ElectrometryAdapter(static, {
        "enabled": True, "mode": "configured", "status": "ASSUMED",
        "transition_dipole_c_m": 2e-27, "angular_factor": .5,
        "field_amplitude_convention": "peak", "quantum_efficiency": .8,
        "signal_probe_power_scale": 1.0, "probe_detuning_search_points": 11,
    })


def test_default_opacity_keeps_legacy_transmission(opacity_point):
    _, point = opacity_point
    raw = point.raw
    alpha, _ = observables.absorption_coefficient(
        raw["chi_bar"], raw["k_vec"], raw["N"], dipole=raw["dipole"], line_strength=.001)
    np.testing.assert_array_equal(point.readout["T_trans"],
                                  np.exp(-alpha * (point.params["cell_mm"] * 1e-3)))
    assert point.result["effective_opacity"]["value"] == .001
    assert point.result["effective_opacity"]["calibration_status"] == "UNCALIBRATED"
    assert "spectrum" not in point.result


def test_configured_opacity_scales_once_and_reuses_atomic_phasors(opacity_point, monkeypatch):
    static, point = opacity_point
    adapter = _configured_adapter(static)
    nominal = adapter.evaluate(point)
    assert nominal["available"]
    params = dict(point.params, effective_opacity_participation=.002)
    # A deliberately poisoned deprecated alias catches accidental reuse/doubling.
    raw = dict(point.raw, ls=12345.0)
    with monkeypatch.context() as check:
        check.setattr(static.scheme, "compute", lambda *a, **k: pytest.fail("extra static solve"))
        check.setattr("analysis.rydberg_cell_heating.adapters.weak_signal_response",
                      lambda *a, **k: pytest.fail("extra configured finite-IF solve"))
        changed = static.point_from_raw(
            raw, params, 40.0, False,
            atomic_response=point.readout["superhet"]["atomic_response"],
            atomic_response_cache=point.atomic_response_cache)
        result = adapter.evaluate(changed)
    assert result["available"]
    for old, new in zip(nominal["detuning_scan"], result["detuning_scan"]):
        assert old["probe_detuning_mhz"] == new["probe_detuning_mhz"]
        # dT/dE is proportional to opacity*T; this independent ratio detects
        # both missing opacity and applying it twice to a fixed atomic phasor.
        expected_ratio = 2.0 * new["transmission"] / old["transmission"]
        assert new["current_responsivity_a_per_v_m"] == pytest.approx(
            old["current_responsivity_a_per_v_m"] * expected_ratio, rel=1e-12)


@pytest.mark.parametrize("mode", ["configured", "scheme"])
def test_zero_opacity_has_no_sensitivity_optimum(opacity_point, mode):
    static, point = opacity_point
    zero = static.point_from_raw(
        point.raw, dict(point.params, effective_opacity_participation=0.0), 40.0, False,
        atomic_response_cache=point.atomic_response_cache)
    adapter = (_configured_adapter(static) if mode == "configured" else
               ElectrometryAdapter(static, {"enabled": True, "mode": "scheme"}))
    result = adapter.evaluate(zero)
    assert not result["available"]
    assert result["status"] == "PREDICTED"
    assert "optimum_probe_detuning_mhz" not in result
    np.testing.assert_array_equal(zero.readout["T_trans"], np.ones_like(zero.readout["x"]))
    assert zero.result["metrics"]["max_spectral_slope_per_mhz"]["value"] == 0.0
    json.dumps(zero.result, allow_nan=False)
    json.dumps(result, allow_nan=False)
