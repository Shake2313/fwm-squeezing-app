"""Experiment-facing helpers for heated-cell Rydberg electrometry.

The Bloch-equation solver deliberately works with an atomic-vapor temperature
and a local absorption coefficient.  A laboratory setup instead exposes heater
set points, one or more contact sensors, beam geometry, and calibrated RF source
power.  This module keeps those layers explicit and supplies small, independently
testable conversions between them.

No function in this module assumes a full Zeeman, time-domain heterodyne, or horn
near-field model.  Rates ending in ``_mhz`` are ordinary-frequency linewidths
(``rate / 2 pi``), matching :mod:`gabes.schemes.rydberg`.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

import numpy as np

from . import constants, species


_DENSITY_SCALE_M3 = 1.0e16


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


def _readonly(values: np.ndarray | Sequence[float]) -> np.ndarray:
    result = np.array(values, dtype=float, copy=True)
    result.setflags(write=False)
    return result


def _least_squares_covariance(weighted_design: np.ndarray) -> np.ndarray:
    inverse_design = np.linalg.pinv(weighted_design)
    return inverse_design @ inverse_design.T


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
        density = cold_density * cold_k / (temp + 273.15)
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


def gaussian_overlap_area_m2(
    probe_radius_m: float, coupling_radius_m: float | None = None
) -> float:
    """Effective transverse area for normalized Gaussian intensity weights.

    Radii are 1/e^2 *intensity* radii.  With one beam the weight is its
    intensity, giving ``pi*w^2/2``.  With two beams the weight is the product of
    their normalized intensities.
    """

    probe = _finite_scalar(probe_radius_m, "probe_radius_m")
    if probe <= 0.0:
        raise ValueError("probe_radius_m must be positive.")
    if coupling_radius_m is None:
        return float(np.pi * probe**2 / 2.0)
    coupling = _finite_scalar(coupling_radius_m, "coupling_radius_m")
    if coupling <= 0.0:
        raise ValueError("coupling_radius_m must be positive.")
    return float(np.pi / (2.0 * (probe**-2 + coupling**-2)))


@dataclass(frozen=True)
class EffectiveAtomNumber:
    """Weighted participating-atom estimate and its geometric factors."""

    atoms: float
    column_density_m2: float
    effective_area_m2: float
    participation_fraction: float
    overlap_efficiency: float


def effective_atom_number(
    density_m3: float | Sequence[float] | AxialCellProfile,
    *,
    length_m: float | None = None,
    z_m: Sequence[float] | None = None,
    probe_radius_m: float,
    coupling_radius_m: float | None = None,
    participation_fraction: float = 1.0,
    overlap_efficiency: float = 1.0,
) -> EffectiveAtomNumber:
    """Estimate ``integral n(r) W(r) dV`` for Gaussian beam overlap.

    A scalar density requires ``length_m``.  An axial density array requires
    ``z_m``.  Passing :class:`AxialCellProfile` supplies both automatically.
    ``participation_fraction`` can represent isotope/velocity/Zeeman selection;
    ``overlap_efficiency`` is a separately reported one-sided geometry penalty.
    """

    participation = _finite_scalar(participation_fraction, "participation_fraction")
    overlap = _finite_scalar(overlap_efficiency, "overlap_efficiency")
    if not 0.0 <= participation <= 1.0:
        raise ValueError("participation_fraction must lie in [0, 1].")
    if not 0.0 <= overlap <= 1.0:
        raise ValueError("overlap_efficiency must lie in [0, 1].")

    if isinstance(density_m3, AxialCellProfile):
        if z_m is not None or length_m is not None:
            raise ValueError(
                "Do not supply z_m or length_m with an AxialCellProfile."
            )
        density_values = density_m3.density_m3
        z_values = density_m3.z_m
        column = float(np.trapezoid(density_values, z_values))
    else:
        density_values = np.asarray(density_m3, dtype=float)
        if not np.all(np.isfinite(density_values)) or np.any(density_values < 0.0):
            raise ValueError("density_m3 must be finite and non-negative.")
        if density_values.ndim == 0:
            if length_m is None or z_m is not None:
                raise ValueError(
                    "A scalar density requires length_m and does not accept z_m."
                )
            length = _finite_scalar(length_m, "length_m")
            if length <= 0.0:
                raise ValueError("length_m must be positive.")
            column = float(density_values) * length
        else:
            if density_values.ndim != 1 or z_m is None or length_m is not None:
                raise ValueError(
                    "An axial density array requires z_m and does not accept length_m."
                )
            z_values = np.asarray(z_m, dtype=float)
            if (
                z_values.shape != density_values.shape
                or not np.all(np.isfinite(z_values))
                or np.any(np.diff(z_values) <= 0.0)
            ):
                raise ValueError(
                    "z_m must be finite, strictly increasing, and match density_m3."
                )
            column = float(np.trapezoid(density_values, z_values))

    area = gaussian_overlap_area_m2(probe_radius_m, coupling_radius_m)
    atoms = column * area * participation * overlap
    return EffectiveAtomNumber(
        atoms=float(atoms),
        column_density_m2=column,
        effective_area_m2=area,
        participation_fraction=participation,
        overlap_efficiency=overlap,
    )


@dataclass(frozen=True)
class SAMFieldCalibration:
    """Far-field standard-antenna-method electric-field calibration."""

    field_v_m: float
    standard_uncertainty_v_m: float
    relative_standard_uncertainty: float
    power_at_antenna_w: float
    net_gain_dbi: float
    amplitude_convention: str
    far_field_ratio: float | None
    warning: str | None

    @property
    def field_n_v_cm(self) -> float:
        """Field in nV/cm."""

        return self.field_v_m * 1.0e7

    @property
    def field_nv_cm(self) -> float:
        """Field in nV/cm (compact spelling)."""

        return self.field_n_v_cm

    @property
    def standard_uncertainty_n_v_cm(self) -> float:
        return self.standard_uncertainty_v_m * 1.0e7

    @property
    def standard_uncertainty_nv_cm(self) -> float:
        return self.standard_uncertainty_n_v_cm


def sam_field_calibration(
    source_power_dbm: float,
    antenna_gain_dbi: float,
    distance_m: float,
    *,
    cable_loss_db: float = 0.0,
    additional_loss_db: float = 0.0,
    field_correction: float = 1.0,
    source_power_std_db: float = 0.0,
    antenna_gain_std_db: float = 0.0,
    cable_loss_std_db: float = 0.0,
    additional_loss_std_db: float = 0.0,
    distance_std_m: float = 0.0,
    field_correction_std: float = 0.0,
    amplitude_convention: str = "rms",
    frequency_hz: float | None = None,
    antenna_max_dimension_m: float | None = None,
) -> SAMFieldCalibration:
    """Convert generator power to a SAM free-space field with uncertainty.

    For RMS amplitude, ``E = sqrt(30 P G) / r``.  Peak amplitude is larger by
    ``sqrt(2)``.  Cable/additional losses are positive dB quantities subtracted
    from source power.  Independent one-standard-deviation inputs are propagated
    by first-order log uncertainty.

    When frequency and maximum antenna dimension are supplied, ``far_field_ratio``
    is ``r / (2 D^2/lambda)`` and a warning is emitted below unity.
    """

    power_dbm = _finite_scalar(source_power_dbm, "source_power_dbm")
    gain_dbi = _finite_scalar(antenna_gain_dbi, "antenna_gain_dbi")
    distance = _finite_scalar(distance_m, "distance_m")
    cable_loss = _finite_scalar(cable_loss_db, "cable_loss_db")
    additional_loss = _finite_scalar(additional_loss_db, "additional_loss_db")
    correction = _finite_scalar(field_correction, "field_correction")
    if distance <= 0.0 or correction <= 0.0:
        raise ValueError("distance_m and field_correction must be positive.")
    if cable_loss < 0.0 or additional_loss < 0.0:
        raise ValueError("Losses must be non-negative dB values.")
    uncertainty_values = {
        "source_power_std_db": source_power_std_db,
        "antenna_gain_std_db": antenna_gain_std_db,
        "cable_loss_std_db": cable_loss_std_db,
        "additional_loss_std_db": additional_loss_std_db,
        "distance_std_m": distance_std_m,
        "field_correction_std": field_correction_std,
    }
    uncertainty = {}
    for name, value in uncertainty_values.items():
        parsed = _finite_scalar(value, name)
        if parsed < 0.0:
            raise ValueError(f"{name} must be non-negative.")
        uncertainty[name] = parsed

    convention = str(amplitude_convention).strip().lower()
    if convention not in {"rms", "peak"}:
        raise ValueError("amplitude_convention must be 'rms' or 'peak'.")
    power_at_antenna_w = 1.0e-3 * 10.0 ** (
        (power_dbm - cable_loss - additional_loss) / 10.0
    )
    gain_linear = 10.0 ** (gain_dbi / 10.0)
    coefficient = 30.0 if convention == "rms" else 60.0
    field = correction * np.sqrt(coefficient * power_at_antenna_w * gain_linear) / distance

    db_to_log_field = np.log(10.0) / 20.0
    relative_variance = db_to_log_field**2 * (
        uncertainty["source_power_std_db"] ** 2
        + uncertainty["antenna_gain_std_db"] ** 2
        + uncertainty["cable_loss_std_db"] ** 2
        + uncertainty["additional_loss_std_db"] ** 2
    )
    relative_variance += (uncertainty["distance_std_m"] / distance) ** 2
    relative_variance += (
        uncertainty["field_correction_std"] / correction
    ) ** 2
    relative_std = float(np.sqrt(relative_variance))

    far_field_ratio = None
    warning = None
    if (frequency_hz is None) != (antenna_max_dimension_m is None):
        raise ValueError(
            "frequency_hz and antenna_max_dimension_m must be supplied together."
        )
    if frequency_hz is not None and antenna_max_dimension_m is not None:
        frequency = _finite_scalar(frequency_hz, "frequency_hz")
        dimension = _finite_scalar(
            antenna_max_dimension_m, "antenna_max_dimension_m"
        )
        if frequency <= 0.0 or dimension <= 0.0:
            raise ValueError(
                "frequency_hz and antenna_max_dimension_m must be positive."
            )
        far_field_distance = 2.0 * dimension**2 * frequency / constants.C_LIGHT
        far_field_ratio = float(distance / far_field_distance)
        if far_field_ratio < 1.0:
            warning = (
                "Observation point is inside the 2D^2/lambda far-field "
                "criterion; the SAM free-space formula may be biased."
            )

    return SAMFieldCalibration(
        field_v_m=float(field),
        standard_uncertainty_v_m=float(field * relative_std),
        relative_standard_uncertainty=relative_std,
        power_at_antenna_w=float(power_at_antenna_w),
        net_gain_dbi=float(gain_dbi - cable_loss - additional_loss),
        amplitude_convention=convention,
        far_field_ratio=far_field_ratio,
        warning=warning,
    )


@dataclass(frozen=True)
class TemperatureDensityDephasingModel:
    """Linear phenomenological Rydberg-dephasing model."""

    baseline_mhz: float
    temperature_slope_mhz_per_c: float
    density_slope_mhz_per_1e16_m3: float
    reference_temp_c: float
    reference_density_m3: float
    clip_nonnegative: bool = True

    def evaluate(
        self,
        temperature_c: float | Sequence[float],
        density_m3: float | Sequence[float],
    ) -> np.ndarray | float:
        temperature = np.asarray(temperature_c, dtype=float)
        density = np.asarray(density_m3, dtype=float)
        try:
            temperature, density = np.broadcast_arrays(temperature, density)
        except ValueError as exc:
            raise ValueError("temperature_c and density_m3 are not broadcastable.") from exc
        if (
            not np.all(np.isfinite(temperature))
            or not np.all(np.isfinite(density))
            or np.any(temperature <= -273.15)
            or np.any(density < 0.0)
        ):
            raise ValueError(
                "temperature_c must be finite and above absolute zero; "
                "density_m3 must be finite and non-negative."
            )
        result = (
            float(self.baseline_mhz)
            + float(self.temperature_slope_mhz_per_c)
            * (temperature - float(self.reference_temp_c))
            + float(self.density_slope_mhz_per_1e16_m3)
            * (density - float(self.reference_density_m3))
            / _DENSITY_SCALE_M3
        )
        if self.clip_nonnegative:
            result = np.maximum(result, 0.0)
        if result.ndim == 0:
            return float(result)
        return result


@dataclass(frozen=True)
class DephasingFitResult:
    """Weighted least-squares result for a temperature/density model."""

    model: TemperatureDensityDephasingModel
    covariance: np.ndarray
    fitted_mhz: np.ndarray
    residuals_mhz: np.ndarray
    rmse_mhz: float
    reduced_chi_squared: float | None
    rank: int
    condition_number: float
    identifiable: bool
    parameter_names: tuple[str, ...]
    warning: str | None


def fit_temperature_density_dephasing(
    temperature_c: Sequence[float],
    density_m3: Sequence[float],
    dephasing_mhz: Sequence[float],
    *,
    dephasing_std_mhz: Sequence[float] | None = None,
    reference_temp_c: float | None = None,
    reference_density_m3: float | None = None,
    fit_temperature: bool = True,
    fit_density: bool = True,
    clip_nonnegative: bool = True,
) -> DephasingFitResult:
    """Fit a linear phenomenological dephasing budget with NumPy only.

    Temperature and vapor density are often strongly correlated in a heating
    sweep.  ``identifiable`` is therefore false when the design is rank deficient
    or has condition number above ``1e10``; the numerical fit is still returned
    so reports can expose, rather than hide, that degeneracy.  Disable one slope
    when the available data cannot identify both.
    """

    temperature = np.asarray(temperature_c, dtype=float)
    density = np.asarray(density_m3, dtype=float)
    measured = np.asarray(dephasing_mhz, dtype=float)
    if (
        temperature.ndim != 1
        or density.shape != temperature.shape
        or measured.shape != temperature.shape
        or temperature.size < 2
    ):
        raise ValueError(
            "temperature_c, density_m3, and dephasing_mhz must be equal-length "
            "1-D arrays with at least two points."
        )
    if (
        not np.all(np.isfinite(temperature))
        or not np.all(np.isfinite(density))
        or not np.all(np.isfinite(measured))
        or np.any(temperature <= -273.15)
        or np.any(density < 0.0)
        or np.any(measured < 0.0)
    ):
        raise ValueError("Fit inputs contain an invalid temperature, density, or rate.")
    ref_temp = (
        float(np.mean(temperature))
        if reference_temp_c is None
        else _temperature_c(reference_temp_c, "reference_temp_c")
    )
    ref_density = (
        float(np.mean(density))
        if reference_density_m3 is None
        else _finite_scalar(reference_density_m3, "reference_density_m3")
    )
    if ref_density < 0.0:
        raise ValueError("reference_density_m3 must be non-negative.")

    columns = [np.ones_like(temperature)]
    names = ["baseline_mhz"]
    if fit_temperature:
        columns.append(temperature - ref_temp)
        names.append("temperature_slope_mhz_per_c")
    if fit_density:
        columns.append((density - ref_density) / _DENSITY_SCALE_M3)
        names.append("density_slope_mhz_per_1e16_m3")
    design = np.column_stack(columns)

    if dephasing_std_mhz is None:
        sigma = None
        weighted_design = design
        weighted_y = measured
    else:
        sigma = np.asarray(dephasing_std_mhz, dtype=float)
        if (
            sigma.shape != measured.shape
            or not np.all(np.isfinite(sigma))
            or np.any(sigma <= 0.0)
        ):
            raise ValueError(
                "dephasing_std_mhz must be finite, positive, and match the data."
            )
        weighted_design = design / sigma[:, None]
        weighted_y = measured / sigma

    coefficients, _, rank, singular_values = np.linalg.lstsq(
        weighted_design, weighted_y, rcond=None
    )
    fitted = design @ coefficients
    residuals = measured - fitted
    dof = measured.size - int(rank)
    if singular_values.size == 0 or singular_values[-1] == 0.0:
        condition = float("inf")
    else:
        condition = float(singular_values[0] / singular_values[-1])

    covariance_basis = _least_squares_covariance(weighted_design)
    if sigma is None:
        residual_variance = (
            float(np.sum(residuals**2) / dof) if dof > 0 else float("nan")
        )
        covariance = covariance_basis * residual_variance
        reduced_chi_squared = None
    else:
        covariance = covariance_basis
        reduced_chi_squared = (
            float(np.sum((residuals / sigma) ** 2) / dof) if dof > 0 else None
        )

    values = dict(
        baseline_mhz=float(coefficients[0]),
        temperature_slope_mhz_per_c=0.0,
        density_slope_mhz_per_1e16_m3=0.0,
    )
    for name, value in zip(names, coefficients):
        values[name] = float(value)
    model = TemperatureDensityDephasingModel(
        **values,
        reference_temp_c=ref_temp,
        reference_density_m3=ref_density,
        clip_nonnegative=bool(clip_nonnegative),
    )
    identifiable = bool(rank == design.shape[1] and condition < 1.0e10)
    warning = None
    if not identifiable:
        warning = (
            "Temperature and density coefficients are not independently "
            "identifiable from this data; fix one slope or add an independent "
            "density/temperature measurement."
        )
    return DephasingFitResult(
        model=model,
        covariance=_readonly(covariance),
        fitted_mhz=_readonly(fitted),
        residuals_mhz=_readonly(residuals),
        rmse_mhz=float(np.sqrt(np.mean(residuals**2))),
        reduced_chi_squared=reduced_chi_squared,
        rank=int(rank),
        condition_number=condition,
        identifiable=identifiable,
        parameter_names=tuple(names),
        warning=warning,
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
    "DephasingFitResult",
    "EffectiveAtomNumber",
    "LinearTemperatureCalibration",
    "SAMFieldCalibration",
    "SegmentedBeerLambert",
    "TemperatureDensityDephasingModel",
    "axial_cell_profile",
    "effective_atom_number",
    "fit_effective_temperature_calibration",
    "fit_temperature_density_dephasing",
    "gaussian_overlap_area_m2",
    "integrate_beer_lambert",
    "linear_axial_cell_profile",
    "resolve_cell_temperature",
    "sam_field_calibration",
]
