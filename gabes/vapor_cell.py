"""Heated vapor cell: temperatures, beam-path density and axial Beer--Lambert.

The Bloch-equation solvers take one vapor temperature (Doppler width, collisional
rates) and a number density.  A laboratory cell instead exposes heater set points
and contact sensors, and in a sealed cell the coldest reservoir fixes the vapor
pressure, so the density along the beam is ``n_sat(T_cold)·T_cold/T_path`` rather
than the saturated density at the beam-path temperature.  This module keeps those
layers explicit for every scheme: the temperature state and the heater-to-vapor
calibration, the cold-spot-limited density, axial temperature/density profiles,
and the trapezoidal Beer--Lambert integral along ``z``.

Temperatures ending in ``_c`` are °C, ``_k`` are K; densities are m⁻³.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

import numpy as np

from . import species
from ._arrays import readonly_float_array as _readonly


def _finite_scalar(value: float, name: str) -> float:
    result = float(value)
    if not np.isfinite(result):
        raise ValueError(f"{name} must be finite.")
    return result


def _temperature_c(value: float, name: str) -> float:
    result = _finite_scalar(value, name)
    if result <= -273.15:
        raise ValueError(f"{name} must be above absolute zero.")
    return result


def _least_squares_covariance(weighted_design: np.ndarray) -> np.ndarray:
    inverse_design = np.linalg.pinv(weighted_design)
    return inverse_design @ inverse_design.T


def cold_spot_limited_density(cold_spot_density_m3, cold_spot_k, path_k):
    """Vapor density at ``path_k`` [K] when a colder reservoir at ``cold_spot_k``
    [K] fixes the pressure.

    The reservoir holds the saturated pressure P = n_cold k_B T_cold; the vapor
    along the beam is an ideal gas at that pressure, so n = n_cold·T_cold/T_path.
    ``cold_spot_density_m3`` is the saturated density at the reservoir (any
    isotope weighting already applied); ``path_k`` may be an array.
    """

    return cold_spot_density_m3 * cold_spot_k / path_k


@dataclass(frozen=True)
class CellTemperatureState:
    """One temperature observation with laboratory and vapor values separated.

    ``sensor_temperatures_c`` contains the reported readings, while
    ``corrected_sensor_temperatures_c`` includes the supplied sensor offsets.
    ``effective_vapor_temp_c`` is the temperature used for Doppler/dephasing
    calculations.  ``cold_spot_temp_c`` controls the equilibrium vapor pressure.
    The two need not be equal in a heated cell.
    """

    heater_setpoint_c: float
    sensor_temperatures_c: tuple[float, ...]
    sensor_offsets_c: tuple[float, ...]
    corrected_sensor_temperatures_c: tuple[float, ...]
    effective_vapor_temp_c: float
    cold_spot_temp_c: float
    effective_source: str
    cold_spot_source: str

    def saturated_cold_spot_density_m3(self, isotope=species.RB85) -> float:
        """Saturated density at the cold spot (not a spatial average)."""

        return float(
            species.number_density(isotope, self.cold_spot_temp_c + 273.15)
        )

    def beam_path_density_m3(self, isotope=species.RB85) -> float:
        """Density at the effective vapor temperature with the cold spot fixing
        the pressure: n = n_sat(T_cold)·T_cold/T_eff (see
        :func:`cold_spot_limited_density`)."""

        return cold_spot_limited_density(
            self.saturated_cold_spot_density_m3(isotope),
            self.cold_spot_temp_c + 273.15,
            self.effective_vapor_temp_c + 273.15,
        )


def resolve_cell_temperature(
    heater_setpoint_c: float,
    sensor_temperatures_c: Sequence[float] = (),
    *,
    sensor_offsets_c: Sequence[float] | None = None,
    sensor_weights: Sequence[float] | None = None,
    effective_vapor_temp_c: float | None = None,
    cold_spot_temp_c: float | None = None,
) -> CellTemperatureState:
    """Resolve explicitly labelled heater, sensor, effective, and cold-spot values.

    Missing effective temperature defaults to the weighted mean of corrected
    sensors (or the heater set point when there is no sensor).  Missing cold-spot
    temperature defaults to the coldest corrected sensor (or the set point).
    The returned ``*_source`` strings make those fallbacks visible to reports.

    Sensor offsets follow ``corrected = reported + offset``.
    """

    setpoint = _temperature_c(heater_setpoint_c, "heater_setpoint_c")
    sensors = tuple(
        _temperature_c(value, f"sensor_temperatures_c[{index}]")
        for index, value in enumerate(sensor_temperatures_c)
    )
    if sensor_offsets_c is None:
        offsets = (0.0,) * len(sensors)
    else:
        offsets = tuple(
            _finite_scalar(value, f"sensor_offsets_c[{index}]")
            for index, value in enumerate(sensor_offsets_c)
        )
        if len(offsets) != len(sensors):
            raise ValueError(
                "sensor_offsets_c must have the same length as "
                "sensor_temperatures_c."
            )
    corrected = tuple(
        _temperature_c(value + offset, f"corrected sensor {index}")
        for index, (value, offset) in enumerate(zip(sensors, offsets))
    )

    if sensor_weights is None:
        weights = None
    else:
        weights_array = np.asarray(sensor_weights, dtype=float)
        if weights_array.shape != (len(sensors),):
            raise ValueError(
                "sensor_weights must have the same length as "
                "sensor_temperatures_c."
            )
        if not np.all(np.isfinite(weights_array)) or np.any(weights_array < 0.0):
            raise ValueError("sensor_weights must be finite and non-negative.")
        if float(np.sum(weights_array)) <= 0.0:
            raise ValueError("sensor_weights must contain at least one positive value.")
        weights = weights_array

    if effective_vapor_temp_c is not None:
        effective = _temperature_c(
            effective_vapor_temp_c, "effective_vapor_temp_c"
        )
        effective_source = "explicit"
    elif corrected:
        effective = float(np.average(np.asarray(corrected), weights=weights))
        effective_source = "weighted_sensor_mean" if weights is not None else "sensor_mean"
    else:
        effective = setpoint
        effective_source = "heater_setpoint_fallback"

    if cold_spot_temp_c is not None:
        cold_spot = _temperature_c(cold_spot_temp_c, "cold_spot_temp_c")
        cold_source = "explicit"
    elif corrected:
        cold_spot = min(corrected)
        cold_source = "coldest_sensor"
    else:
        cold_spot = setpoint
        cold_source = "heater_setpoint_fallback"

    return CellTemperatureState(
        heater_setpoint_c=setpoint,
        sensor_temperatures_c=sensors,
        sensor_offsets_c=offsets,
        corrected_sensor_temperatures_c=corrected,
        effective_vapor_temp_c=effective,
        cold_spot_temp_c=cold_spot,
        effective_source=effective_source,
        cold_spot_source=cold_source,
    )


@dataclass(frozen=True)
class AxialCellProfile:
    """Temperature and vapor density sampled along the optical path."""

    z_m: np.ndarray
    temperature_c: np.ndarray
    density_m3: np.ndarray
    density_mode: str
    cold_spot_temp_c: float
    profile_min_temp_c: float
    cold_spot_density_m3: float

    @property
    def length_m(self) -> float:
        return float(self.z_m[-1] - self.z_m[0])

    @property
    def column_density_m2(self) -> float:
        return float(np.trapezoid(self.density_m3, self.z_m))


def axial_cell_profile(
    z_m: Sequence[float],
    temperature_c: Sequence[float],
    *,
    isotope=species.RB85,
    density_mode: str = "cold_spot_limited",
    pressure_cold_spot_temp_c: float | None = None,
) -> AxialCellProfile:
    """Build an axial thermal/density profile.

    Supported density modes are:

    ``cold_spot_limited``
        The cold spot fixes a uniform vapor pressure.  Local density then scales
        as ``P_cold / (k_B T(z))``.  This is the recommended sealed-cell model.
    ``uniform_cold_spot``
        Use the saturated cold-spot density everywhere, a common lumped model.
    ``local_saturation``
        Evaluate saturated vapor density independently at every location.  This
        is useful as an upper-envelope diagnostic but is not pressure-equilibrated.

    ``pressure_cold_spot_temp_c`` may name a colder reservoir outside the sampled
    optical path.  When omitted, the minimum sampled temperature sets pressure;
    a declared pressure cold spot warmer than any sampled point is rejected.
    """

    z = np.asarray(z_m, dtype=float)
    temp = np.asarray(temperature_c, dtype=float)
    if z.ndim != 1 or temp.ndim != 1 or z.size != temp.size or z.size < 2:
        raise ValueError(
            "z_m and temperature_c must be one-dimensional arrays of equal "
            "length with at least two points."
        )
    if not np.all(np.isfinite(z)) or not np.all(np.isfinite(temp)):
        raise ValueError("z_m and temperature_c must contain only finite values.")
    if np.any(np.diff(z) <= 0.0):
        raise ValueError("z_m must be strictly increasing.")
    if np.any(temp <= -273.15):
        raise ValueError("temperature_c must be above absolute zero.")

    profile_min_c = float(np.min(temp))
    cold_c = (
        profile_min_c
        if pressure_cold_spot_temp_c is None
        else _temperature_c(
            pressure_cold_spot_temp_c, "pressure_cold_spot_temp_c")
    )
    if cold_c > profile_min_c + 1.0e-12:
        raise ValueError(
            "pressure_cold_spot_temp_c cannot be warmer than the minimum "
            "sampled axial temperature")
    cold_k = cold_c + 273.15
    cold_density = float(species.number_density(isotope, cold_k))
    mode = str(density_mode).strip().lower()
    if mode == "cold_spot_limited":
        density = cold_spot_limited_density(cold_density, cold_k, temp + 273.15)
    elif mode == "uniform_cold_spot":
        density = np.full_like(temp, cold_density)
    elif mode == "local_saturation":
        density = np.asarray(
            [species.number_density(isotope, value + 273.15) for value in temp],
            dtype=float,
        )
    else:
        raise ValueError(
            "density_mode must be 'cold_spot_limited', 'uniform_cold_spot', "
            "or 'local_saturation'."
        )
    return AxialCellProfile(
        z_m=_readonly(z),
        temperature_c=_readonly(temp),
        density_m3=_readonly(density),
        density_mode=mode,
        cold_spot_temp_c=cold_c,
        profile_min_temp_c=profile_min_c,
        cold_spot_density_m3=cold_density,
    )


def linear_axial_cell_profile(
    length_m: float,
    left_temp_c: float,
    right_temp_c: float,
    *,
    points: int = 101,
    isotope=species.RB85,
    density_mode: str = "cold_spot_limited",
    pressure_cold_spot_temp_c: float | None = None,
) -> AxialCellProfile:
    """Convenience builder for a linearly varying cell temperature."""

    length = _finite_scalar(length_m, "length_m")
    if length <= 0.0:
        raise ValueError("length_m must be positive.")
    if isinstance(points, bool) or int(points) != points or int(points) < 2:
        raise ValueError("points must be an integer of at least 2.")
    left = _temperature_c(left_temp_c, "left_temp_c")
    right = _temperature_c(right_temp_c, "right_temp_c")
    z = np.linspace(0.0, length, int(points))
    temperature = np.linspace(left, right, int(points))
    return axial_cell_profile(
        z, temperature, isotope=isotope, density_mode=density_mode,
        pressure_cold_spot_temp_c=pressure_cold_spot_temp_c,
    )


@dataclass(frozen=True)
class SegmentedBeerLambert:
    """Segment-wise and total Beer--Lambert attenuation."""

    segment_optical_depth: np.ndarray
    segment_transmission: np.ndarray
    optical_depth: np.ndarray | float
    transmission: np.ndarray | float


def integrate_beer_lambert(
    z_m: Sequence[float], alpha_per_m: np.ndarray | Sequence[float]
) -> SegmentedBeerLambert:
    """Integrate local absorption coefficients along ``z`` by trapezoids.

    The first axis of ``alpha_per_m`` is axial.  Trailing axes (for example a
    detuning scan) are preserved.  Negative coefficients are allowed so the same
    helper can represent gain.
    """

    z = np.asarray(z_m, dtype=float)
    alpha = np.asarray(alpha_per_m, dtype=float)
    if z.ndim != 1 or z.size < 2 or np.any(np.diff(z) <= 0.0):
        raise ValueError("z_m must be a strictly increasing 1-D array.")
    if alpha.ndim < 1 or alpha.shape[0] != z.size:
        raise ValueError("alpha_per_m first dimension must match z_m.")
    if not np.all(np.isfinite(z)) or not np.all(np.isfinite(alpha)):
        raise ValueError("z_m and alpha_per_m must contain only finite values.")
    dz_shape = (z.size - 1,) + (1,) * (alpha.ndim - 1)
    segment_od = 0.5 * (alpha[:-1] + alpha[1:]) * np.diff(z).reshape(dz_shape)
    total_od_array = np.sum(segment_od, axis=0)
    segment_t = np.exp(-segment_od)
    total_t_array = np.exp(-total_od_array)
    total_od: np.ndarray | float
    total_t: np.ndarray | float
    if total_od_array.ndim == 0:
        total_od = float(total_od_array)
        total_t = float(total_t_array)
    else:
        total_od = _readonly(total_od_array)
        total_t = _readonly(total_t_array)
    return SegmentedBeerLambert(
        segment_optical_depth=_readonly(segment_od),
        segment_transmission=_readonly(segment_t),
        optical_depth=total_od,
        transmission=total_t,
    )


@dataclass(frozen=True)
class LinearTemperatureCalibration:
    """Linear map from heater set point to effective vapor temperature."""

    intercept_c: float
    slope: float
    covariance: np.ndarray
    rmse_c: float

    def evaluate(self, heater_setpoint_c: float | Sequence[float]) -> np.ndarray | float:
        value = np.asarray(heater_setpoint_c, dtype=float)
        if not np.all(np.isfinite(value)):
            raise ValueError("heater_setpoint_c must be finite.")
        result = self.intercept_c + self.slope * value
        if result.ndim == 0:
            return float(result)
        return result


def fit_effective_temperature_calibration(
    heater_setpoint_c: Sequence[float],
    effective_vapor_temp_c: Sequence[float],
    *,
    effective_temp_std_c: Sequence[float] | None = None,
) -> LinearTemperatureCalibration:
    """Fit ``T_effective = intercept + slope * T_setpoint``."""

    setpoint = np.asarray(heater_setpoint_c, dtype=float)
    effective = np.asarray(effective_vapor_temp_c, dtype=float)
    if (
        setpoint.ndim != 1
        or effective.shape != setpoint.shape
        or setpoint.size < 2
        or not np.all(np.isfinite(setpoint))
        or not np.all(np.isfinite(effective))
    ):
        raise ValueError("Temperature calibration needs two equal finite 1-D arrays.")
    design = np.column_stack((np.ones_like(setpoint), setpoint))
    if effective_temp_std_c is None:
        weighted_design = design
        weighted_y = effective
        sigma = None
    else:
        sigma = np.asarray(effective_temp_std_c, dtype=float)
        if (
            sigma.shape != effective.shape
            or not np.all(np.isfinite(sigma))
            or np.any(sigma <= 0.0)
        ):
            raise ValueError("effective_temp_std_c must be finite, positive, and match data.")
        weighted_design = design / sigma[:, None]
        weighted_y = effective / sigma
    coefficients, _, rank, _ = np.linalg.lstsq(
        weighted_design, weighted_y, rcond=None
    )
    if rank < 2:
        raise ValueError("Heater set points must contain at least two distinct values.")
    residuals = effective - design @ coefficients
    covariance = _least_squares_covariance(weighted_design)
    if sigma is None and setpoint.size > 2:
        covariance *= float(np.sum(residuals**2) / (setpoint.size - 2))
    return LinearTemperatureCalibration(
        intercept_c=float(coefficients[0]),
        slope=float(coefficients[1]),
        covariance=_readonly(covariance),
        rmse_c=float(np.sqrt(np.mean(residuals**2))),
    )


__all__ = [
    "AxialCellProfile",
    "CellTemperatureState",
    "LinearTemperatureCalibration",
    "SegmentedBeerLambert",
    "axial_cell_profile",
    "cold_spot_limited_density",
    "fit_effective_temperature_calibration",
    "integrate_beer_lambert",
    "linear_axial_cell_profile",
    "resolve_cell_temperature",
]
