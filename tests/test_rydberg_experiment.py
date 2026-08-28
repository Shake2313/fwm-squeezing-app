import numpy as np
import pytest

from gabes import species
from gabes.rydberg_experiment import (
    TemperatureDensityDephasingModel,
    axial_cell_profile,
    effective_atom_number,
    fit_effective_temperature_calibration,
    fit_temperature_density_dephasing,
    gaussian_overlap_area_m2,
    integrate_beer_lambert,
    linear_axial_cell_profile,
    resolve_cell_temperature,
    sam_field_calibration,
)


def test_temperature_state_keeps_setpoint_sensor_effective_and_cold_spot_separate():
    state = resolve_cell_temperature(
        50.0,
        [44.0, 42.0],
        sensor_offsets_c=[1.0, -1.0],
        sensor_weights=[1.0, 3.0],
    )
    assert state.heater_setpoint_c == 50.0
    assert state.sensor_temperatures_c == (44.0, 42.0)
    assert state.corrected_sensor_temperatures_c == (45.0, 41.0)
    assert state.effective_vapor_temp_c == pytest.approx(42.0)
    assert state.cold_spot_temp_c == 41.0
    assert state.effective_source == "weighted_sensor_mean"
    assert state.cold_spot_source == "coldest_sensor"
    assert state.saturated_cold_spot_density_m3() == pytest.approx(
        species.number_density(species.RB85, 314.15)
    )


def test_axial_cold_spot_pressure_and_segmented_beer_lambert():
    profile = linear_axial_cell_profile(
        0.01, 40.0, 50.0, points=5, density_mode="cold_spot_limited"
    )
    assert profile.length_m == pytest.approx(0.01)
    assert profile.cold_spot_temp_c == 40.0
    assert profile.density_m3[0] == pytest.approx(profile.cold_spot_density_m3)
    assert profile.density_m3[-1] == pytest.approx(
        profile.cold_spot_density_m3 * 313.15 / 323.15
    )
    assert profile.column_density_m2 == pytest.approx(
        np.trapezoid(profile.density_m3, profile.z_m)
    )

    external_cold = linear_axial_cell_profile(
        0.01, 50.0, 50.0, points=3,
        density_mode="cold_spot_limited",
        pressure_cold_spot_temp_c=30.0,
    )
    expected_cold_density = species.number_density(species.RB85, 303.15)
    assert external_cold.profile_min_temp_c == 50.0
    assert external_cold.cold_spot_temp_c == 30.0
    assert np.allclose(
        external_cold.density_m3,
        expected_cold_density * 303.15 / 323.15,
    )

    z = np.array([0.0, 0.5, 1.0])
    alpha = np.array([[1.0, 2.0], [1.0, 4.0], [1.0, 6.0]])
    result = integrate_beer_lambert(z, alpha)
    assert np.allclose(result.optical_depth, [1.0, 4.0])
    assert np.allclose(result.transmission, np.exp([-1.0, -4.0]))
    assert np.allclose(np.prod(result.segment_transmission, axis=0), result.transmission)


def test_local_saturation_is_diagnostic_upper_envelope_on_hot_side():
    z = [0.0, 0.5, 1.0]
    temperatures = [40.0, 45.0, 50.0]
    cold_limited = axial_cell_profile(z, temperatures, density_mode="cold_spot_limited")
    local = axial_cell_profile(z, temperatures, density_mode="local_saturation")
    assert local.density_m3[0] == pytest.approx(cold_limited.density_m3[0])
    assert local.density_m3[-1] > cold_limited.density_m3[-1]


def test_gaussian_effective_atom_number_states_geometry_and_participation():
    density = 1.0e16
    length = 0.01
    waist = 1.0e-3
    single_area = np.pi * waist**2 / 2.0
    assert gaussian_overlap_area_m2(waist) == pytest.approx(single_area)
    assert gaussian_overlap_area_m2(waist, waist) == pytest.approx(single_area / 2.0)

    estimate = effective_atom_number(
        density,
        length_m=length,
        probe_radius_m=waist,
        participation_fraction=0.5,
        overlap_efficiency=0.8,
    )
    assert estimate.column_density_m2 == pytest.approx(density * length)
    assert estimate.atoms == pytest.approx(density * length * single_area * 0.4)


def test_sam_field_and_log_uncertainty_with_far_field_warning():
    calibration = sam_field_calibration(
        0.0,
        10.0,
        1.0,
        source_power_std_db=0.2,
        antenna_gain_std_db=0.3,
        distance_std_m=0.01,
        frequency_hz=10.0e9,
        antenna_max_dimension_m=0.2,
    )
    assert calibration.power_at_antenna_w == pytest.approx(1.0e-3)
    assert calibration.field_v_m == pytest.approx(np.sqrt(30.0e-3 * 10.0))
    expected_relative = np.sqrt(
        (np.log(10.0) / 20.0) ** 2 * (0.2**2 + 0.3**2) + 0.01**2
    )
    assert calibration.relative_standard_uncertainty == pytest.approx(expected_relative)
    assert calibration.standard_uncertainty_v_m == pytest.approx(
        calibration.field_v_m * expected_relative
    )
    assert calibration.far_field_ratio < 1.0
    assert calibration.warning is not None
    assert calibration.field_n_v_cm == pytest.approx(calibration.field_v_m * 1e7)


def test_temperature_density_dephasing_fit_recovers_independent_slopes():
    temperature = np.array([20.0, 25.0, 30.0, 35.0, 40.0, 45.0])
    density = np.array([1.0, 4.0, 2.0, 6.0, 3.0, 5.0]) * 1e16
    measured = 0.20 + 0.010 * (temperature - 30.0) + 0.050 * (
        density - 3.0e16
    ) / 1e16
    result = fit_temperature_density_dephasing(
        temperature,
        density,
        measured,
        reference_temp_c=30.0,
        reference_density_m3=3.0e16,
    )
    assert result.identifiable
    assert result.rmse_mhz < 1e-12
    assert result.model.baseline_mhz == pytest.approx(0.20)
    assert result.model.temperature_slope_mhz_per_c == pytest.approx(0.010)
    assert result.model.density_slope_mhz_per_1e16_m3 == pytest.approx(0.050)
    assert result.model.evaluate(35.0, 4e16) == pytest.approx(0.30)


def test_dephasing_fit_flags_temperature_density_collinearity():
    temperature = np.array([20.0, 25.0, 30.0, 35.0])
    density = (temperature - 20.0) * 1e16
    measured = np.array([0.1, 0.2, 0.3, 0.4])
    result = fit_temperature_density_dephasing(
        temperature, density, measured, reference_temp_c=20.0, reference_density_m3=0.0
    )
    assert not result.identifiable
    assert result.warning is not None


def test_dephasing_model_clips_only_when_requested():
    clipped = TemperatureDensityDephasingModel(0.1, -1.0, 0.0, 20.0, 0.0)
    raw = TemperatureDensityDephasingModel(
        0.1, -1.0, 0.0, 20.0, 0.0, clip_nonnegative=False
    )
    assert clipped.evaluate(21.0, 0.0) == 0.0
    assert raw.evaluate(21.0, 0.0) == pytest.approx(-0.9)


def test_effective_temperature_calibration_fit():
    calibration = fit_effective_temperature_calibration(
        [20.0, 30.0, 40.0, 50.0], [19.0, 27.0, 35.0, 43.0]
    )
    assert calibration.intercept_c == pytest.approx(3.0)
    assert calibration.slope == pytest.approx(0.8)
    assert calibration.evaluate(45.0) == pytest.approx(39.0)
