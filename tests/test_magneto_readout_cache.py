"""Density and path length reuse the same microscopic magneto response."""
import copy
import json

import numpy as np
import pytest

from gabes import constants, species
from gabes.schemes import magneto_unified as unified
from gabes.schemes.magneto_unified import (
    UnifiedMagnetoScheme, FIELD_SCAN, OPTICAL_SCAN, TRANSMISSION,
    RESPONSE_DENSITY_M3, RESPONSE_LENGTH_M,
)


def _physical_response(delta, *, b_ut, **conditions):
    """A local Beer-Lambert response with a resolvable +/-1 pT curvature."""
    delta, field = np.asarray(delta), np.asarray(b_ut)*1e-6
    shape = (field.size, delta.size)
    frequency = 377e12
    k = 2*np.pi*frequency/constants.C_LIGHT
    density = conditions["density_m3"]
    length = conditions["cell_length_m"]
    density_ratio = density/RESPONSE_DENSITY_M3
    slope = np.broadcast_to(2e7*(1+np.exp(-((delta-300)/200)**2)), shape)
    theta = field[:, None]*slope
    od = np.broadcast_to(1+4e10*field[:, None], shape)
    probe = 1j*od/(k*RESPONSE_LENGTH_M)*density_ratio
    cp = probe+2*theta/(k*RESPONSE_LENGTH_M)*density_ratio
    cm = probe-2*theta/(k*RESPONSE_LENGTH_M)*density_ratio
    dp = np.full(shape, 1j*4e10/(k*RESPONSE_LENGTH_M))*density_ratio
    dcp = dp+2*slope/(k*RESPONSE_LENGTH_M)*density_ratio
    dcm = dp-2*slope/(k*RESPONSE_LENGTH_M)*density_ratio
    optical_depth = k*length*probe.imag
    transmission = np.exp(-optical_depth)
    return dict(
        transmission=transmission, optical_depth=optical_depth,
        rotation_rad=.25*k*length*(cp-cm).real,
        rotation_slope_rad_per_t=.25*k*length*(dcp-dcm).real,
        transmission_slope_per_t=-transmission*k*length*dp.imag,
        chi_plus=cp, chi_minus=cm, chi_probe=probe,
        chi_plus_slope_per_t=dcp, chi_minus_slope_per_t=dcm,
        chi_probe_slope_per_t=dp, slope_valid=np.ones(shape, bool),
        slope_relative_residual=np.zeros(shape), frequency_hz=frequency,
        density_m3=density, cell_length_m=length,
        reference_centroid_mhz=-1000., isotope=conditions["isotope"],
        line=conditions["line"], ground_f=conditions["ground_f"],
        detuning_reference="test component centroid", model_scope="test local response",
        rotation_status="linear-input circular-birefringence approximation",
        diagnostics=dict(solver_status="physical steady state",
                         optical_excitation_status="weak excitation",
                         magnetic_response_status="analytic DC response"),
        velocity_quadrature=dict(convergence_status="not assessed"),
    )


def _compute_params(scheme, params):
    return {name: params[name] for name in scheme.recompute_keys()}


@pytest.mark.parametrize("axis", [FIELD_SCAN, OPTICAL_SCAN])
def test_density_and_length_are_readout_controls_without_atomic_solves(monkeypatch, axis):
    calls = []

    def response(*args, **kwargs):
        calls.append(kwargs)
        return _physical_response(*args, **kwargs)

    monkeypatch.setattr(unified.magneto_multilevel, "multilevel_response", response)
    scheme = UnifiedMagnetoScheme()
    params = dict(scheme.defaults(), scan_axis=axis, scan_points=11, detuning_points=5)
    cheap = {"density_mode", "density_cm3", "cell_mm"}
    assert cheap.isdisjoint(scheme.recompute_keys())
    assert "temp_c" in scheme.recompute_keys()
    raw = scheme.compute(_compute_params(scheme, params))
    before = copy.deepcopy(raw)
    for changed in (dict(params, density_cm3=3e10, cell_mm=50.),
                    dict(params, density_mode="vapor", cell_mm=20.),
                    dict(params, density_cm3=0.), dict(params, cell_mm=0.)):
        assert _compute_params(scheme, changed) == _compute_params(scheme, params)
        readout = scheme.headless_observables(raw, changed)["magneto_readout"]
        expected_density = (species.number_density(species.RB85, params["temp_c"]+273.15)
                            if changed["density_mode"] == "vapor" else changed["density_cm3"]*1e6)
        assert readout["density_m3"] == expected_density
        assert readout["density_mode"] == changed["density_mode"]
        assert readout["conditions"]["cell_length_m"] == changed["cell_mm"]*1e-3
        assert readout["conditions"]["density_m3"] == (
            None if changed["density_mode"] == "vapor" else expected_density)
        json.dumps(readout, allow_nan=False)
    assert len(calls) == 1
    assert calls[0]["density_m3"] == RESPONSE_DENSITY_M3
    assert calls[0]["cell_length_m"] == RESPONSE_LENGTH_M
    assert raw["conditions"] == before["conditions"]
    for name, value in raw["backend"].items():
        if isinstance(value, np.ndarray):
            np.testing.assert_array_equal(value, before["backend"][name])


def test_transmission_linearity_noise_and_status_follow_display_density(monkeypatch):
    monkeypatch.setattr(unified.magneto_multilevel, "multilevel_response", _physical_response)
    scheme = UnifiedMagnetoScheme()
    params = dict(scheme.defaults(), scan_axis=OPTICAL_SCAN, detuning_points=5,
                  signal_type=TRANSMISSION, density_cm3=1e10)
    raw = scheme.compute(_compute_params(scheme, params))
    ordinary = scheme.headless_observables(raw, params)["magneto_readout"]
    thick = scheme.headless_observables(raw, dict(params, density_cm3=3e10))["magneto_readout"]
    empty = scheme.headless_observables(raw, dict(params, density_cm3=0.))["magneto_readout"]
    assert all(ordinary["sensitivity"]["snr_valid"])
    assert not any(thick["sensitivity"]["snr_valid"])
    assert all(thick["sensitivity"]["valid"])
    assert thick["max_optical_depth"] == pytest.approx(3*ordinary["max_optical_depth"])
    assert np.asarray(thick["sensitivity"]["signal_linearity_relative_error"]).min() > .05
    assert empty["max_optical_depth"] == 0.
    assert "OD > 0.1" not in empty["propagation_status"]
    assert empty["transmission"] == [1.]*len(empty["x"])
    assert not any(empty["sensitivity"]["valid"])
    assert empty["optimum_index"] is None
    assert empty["optimum_snr"] is None
    assert thick["shot_noise"]["detected_photon_flux_per_s"] != ordinary["shot_noise"]["detected_photon_flux_per_s"]


def test_direct_compute_also_uses_an_independent_nonzero_response_reference(monkeypatch):
    monkeypatch.setattr(unified.magneto_multilevel, "multilevel_response", _physical_response)
    scheme = UnifiedMagnetoScheme()
    params = dict(scheme.defaults(), scan_axis=OPTICAL_SCAN, detuning_points=3)
    before = scheme.compute(params)
    after = scheme.compute(dict(params, density_mode="vapor", density_cm3=0., cell_mm=0.))
    assert before["conditions"] == after["conditions"]
    assert before["response_reference"] == after["response_reference"]
    assert after["backend"]["density_m3"] == RESPONSE_DENSITY_M3
    assert after["backend"]["cell_length_m"] == RESPONSE_LENGTH_M
    np.testing.assert_array_equal(before["backend"]["chi_probe"], after["backend"]["chi_probe"])


@pytest.mark.parametrize("axis", [FIELD_SCAN, OPTICAL_SCAN])
@pytest.mark.parametrize("readout_params", [
    dict(density_mode="manual", density_cm3=2e9, cell_mm=25.),
    dict(density_mode="vapor", cell_mm=150.),
    dict(density_mode="manual", density_cm3=0., cell_mm=100.),
    dict(density_mode="manual", density_cm3=2e9, cell_mm=0.),
])
def test_reused_numerical_response_matches_fresh_density_and_length(axis, readout_params):
    scheme = UnifiedMagnetoScheme()
    params = dict(scheme.defaults(), isotope="Rb87", scan_axis=axis, doppler="off",
                  scan_points=11, detuning_points=3, intensity_mw_cm2=.01,
                  laser_detuning_mhz=100., b_offset_ut=.0002, **readout_params)
    raw = scheme.compute(_compute_params(scheme, params))
    readout = scheme.headless_observables(raw, params)["magneto_readout"]
    conditions = dict(raw["conditions"], cell_length_m=params["cell_mm"]*1e-3,
                      density_m3=params["density_cm3"]*1e6 if params["density_mode"] == "manual" else None)
    fresh = unified.magneto_multilevel.multilevel_response(
        raw["detunings_mhz"], b_ut=raw["physical_b_ut"], **conditions)
    select = (raw["operating_field_index"], slice(None)) if axis == OPTICAL_SCAN else (slice(None), 0)
    for name in ("transmission", "rotation_rad", "rotation_slope_rad_per_t", "transmission_slope_per_t"):
        expected = fresh[name][select]
        # Equivalent density-prefactor ordering changes cancellation-level
        # slopes, so include a curve-scale absolute tolerance.
        np.testing.assert_allclose(readout[name], expected, rtol=1e-9,
                                   atol=1e-11*np.max(np.abs(expected))+1e-15)
    assert readout["density_m3"] == fresh["density_m3"]
    _, expected_linearity = unified._finite_amplitude_linearity(
        fresh, raw["physical_b_ut"], params["b_offset_ut"], axis)
    for name in ("rotation", "transmission"):
        budget = readout[f"{name}_sensitivity"]
        np.testing.assert_array_equal(budget["signal_linearity_valid"], expected_linearity[name]["valid"])
        np.testing.assert_array_equal(readout["slope_valid"], fresh["slope_valid"][select])


@pytest.mark.parametrize("change", [dict(density_cm3=-1.), dict(density_cm3=float("nan")),
                                     dict(cell_mm=-1.), dict(cell_mm=float("inf")),
                                     dict(density_mode="invalid")])
def test_readout_rejects_invalid_density_or_length(monkeypatch, change):
    monkeypatch.setattr(unified.magneto_multilevel, "multilevel_response", _physical_response)
    scheme = UnifiedMagnetoScheme()
    params = scheme.defaults()
    raw = scheme.compute(_compute_params(scheme, params))
    with pytest.raises(ValueError):
        scheme.headless_observables(raw, dict(params, **change))
