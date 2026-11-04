"""Regression tests for thermal/dephasing hand-off into report adapters."""

from __future__ import annotations

import numpy as np
import pytest
import json

from analysis.rydberg_cell_heating.adapters import (
    ElectrometryAdapter,
    StaticRydbergAdapter,
    _finite_or_none,
)
from gabes import species
from gabes import observables
from gabes.rydberg_electrometry import weak_signal_response

MHZ = 2.0 * np.pi * 1.0e6


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


# Optical scan axes: the scanned laser takes the coordinate, the other keeps
# its fixed detuning, and the bare two-photon resonance is Δp + Δc = 0.
_AXIS_CASES = {
    "default": {},
    "probe_scan_fixed_coupling": {"coupling_detuning_mhz": 12.0, "mw_detuning_mhz": 1.5},
    "coupling_scan_fixed_probe": {"scan_axis": "Coupling detuning", "probe_detuning_mhz": 12.0},
    "coupling_scan_resonant_probe": {"scan_axis": "Coupling detuning"},
}


@pytest.fixture(scope="module")
def axis_points():
    static = StaticRydbergAdapter()
    return static, {
        (name, view): static.temperature_point(20.0, view, overrides, include_spectrum=True)
        for name, overrides in _AXIS_CASES.items()
        for view in ("EIT", "AT electrometry")
    }


def _independent_hamiltonian(raw, probe_mhz, coupling_mhz, microwave_mhz):
    """diag(H) = [0, -Δp, -(Δp+Δc), -(Δp+Δc+Δ_MW)], built from physical detunings."""
    dp, dc, dmw = probe_mhz * MHZ, coupling_mhz * MHZ, microwave_mhz * MHZ
    hamiltonian = np.zeros((4, 4), dtype=complex)
    hamiltonian[1, 1] = -dp
    hamiltonian[2, 2] = -(dp + dc)
    hamiltonian[3, 3] = -(dp + dc) - dmw
    hamiltonian[0, 1] = hamiltonian[1, 0] = raw["probe_rabi_mhz"] * MHZ / 2.0
    hamiltonian[1, 2] = hamiltonian[2, 1] = raw["coupling_rabi_mhz"] * MHZ / 2.0
    hamiltonian[2, 3] = hamiltonian[3, 2] = raw["lo_rabi_mhz"] * MHZ / 2.0
    return hamiltonian


@pytest.mark.parametrize("case", list(_AXIS_CASES))
def test_configured_phasors_solve_the_physical_optical_point(axis_points, case):
    static, points = axis_points
    point = points[(case, "AT electrometry")]
    raw, params = point.raw, point.params
    adapter = ElectrometryAdapter(static, dict(
        _configured_adapter(static).config, probe_detuning_search_points=5))
    result = adapter.evaluate(point)
    assert result["available"]
    axis = params["scan_axis"]
    assert result["scan_axis"] == axis
    fixed_probe = float(params["probe_detuning_mhz"])
    fixed_coupling = float(params["coupling_detuning_mhz"])
    atom = static.scheme._atom(sum(raw[name] for name in (
        "rydberg_dephasing_mhz", "temperature_dephasing_mhz",
        "density_dephasing_mhz", "transit_mhz")) * MHZ, raw["rf_dephasing_mhz"] * MHZ)
    for row in result["detuning_scan"]:
        scan = row["scan_detuning_mhz"]
        if axis == "Probe detuning":
            assert (row["probe_detuning_mhz"], row["coupling_detuning_mhz"]) == (scan, fixed_coupling)
        else:
            assert (row["probe_detuning_mhz"], row["coupling_detuning_mhz"]) == (fixed_probe, scan)
        np.testing.assert_allclose(
            np.asarray(static.optical_detunings_mhz(params, scan)) * MHZ,
            static.scheme._optical_detunings_at_scan(params, scan * MHZ), rtol=1e-15)
        response = weak_signal_response(
            atom, _independent_hamiltonian(raw, row["probe_detuning_mhz"],
                                           row["coupling_detuning_mhz"],
                                           float(params["mw_detuning_mhz"])),
            2.0 * np.pi * float(params["if_khz"]) * 1.0e3, signal_transition=(2, 3))
        expected = response.real_observable_phasor_per_angular_rabi(
            static.scheme._absorption_quadrature_operator())
        actual = complex(row["atomic_phasor_per_angular_rabi_real"],
                         row["atomic_phasor_per_angular_rabi_imag"])
        assert abs(actual - expected) <= 1e-12 * abs(expected)
    optimum = min(result["detuning_scan"], key=lambda row: row["total_field_asd_nv_cm_per_sqrt_hz"])
    assert result["optimum_scan_detuning_mhz"] == optimum["scan_detuning_mhz"]
    assert result["optimum_probe_detuning_mhz"] == optimum["probe_detuning_mhz"]
    assert result["optimum_coupling_detuning_mhz"] == optimum["coupling_detuning_mhz"]


@pytest.mark.parametrize("case", list(_AXIS_CASES))
def test_scheme_mode_exports_the_physical_optimum(axis_points, case):
    static, points = axis_points
    point = points[(case, "AT electrometry")]
    superhet = point.readout["superhet"]
    result = ElectrometryAdapter(static, {"enabled": True, "mode": "scheme"}).evaluate(point)
    assert result["available"]
    assert result["scan_axis"] == point.params["scan_axis"]
    assert result["optimum_scan_detuning_mhz"] == superhet["optimum_scan_detuning_mhz"]
    assert result["optimum_probe_detuning_mhz"] == superhet["optimum_probe_detuning_mhz"]
    assert result["optimum_coupling_detuning_mhz"] == superhet["optimum_coupling_detuning_mhz"]
    if point.params["scan_axis"] == "Coupling detuning":
        # The probe never moves in a coupling scan.
        assert result["optimum_probe_detuning_mhz"] == pytest.approx(
            float(point.params["probe_detuning_mhz"]), abs=1e-12)


@pytest.mark.parametrize("view", ["EIT", "AT electrometry"])
@pytest.mark.parametrize("case", list(_AXIS_CASES))
def test_static_summary_reads_the_two_photon_resonance(axis_points, case, view):
    static, points = axis_points
    point = points[(case, view)]
    readout, metrics = point.readout, point.result["metrics"]
    scan = point.result["optical_scan"]
    fixed = ("coupling_detuning_mhz" if point.params["scan_axis"] == "Probe detuning"
             else "probe_detuning_mhz")
    assert scan["axis"] == point.params["scan_axis"]
    assert scan["two_photon_resonance_scan_mhz"] == -float(point.params[fixed])
    assert point.result["spectrum"]["scan_axis"] == point.params["scan_axis"]
    assert metrics["transmission_at_resonance"]["value"] == (
        readout["rydberg_readout"]["resonance_transmission"])
    if view == "EIT":
        assert metrics["eit_linewidth_mhz"]["value"] == _finite_or_none(readout["width"])
