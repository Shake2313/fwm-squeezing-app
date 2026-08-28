"""Regression tests for thermal/dephasing hand-off into report adapters."""

from __future__ import annotations

import numpy as np
import pytest

from analysis.rydberg_cell_heating.adapters import (
    ElectrometryAdapter,
    StaticRydbergAdapter,
)
from gabes import species


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
