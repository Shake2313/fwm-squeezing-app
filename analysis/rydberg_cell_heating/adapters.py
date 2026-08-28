"""Adapters between the report workflow and GABES Rydberg implementations.

The current :class:`RydbergEITScheme` exposes scalar observables publicly, but
keeps its transmission arrays behind private helpers.  That dependency is
deliberately confined to :class:`StaticRydbergAdapter`; the report generator and
CSV layer do not depend on scheme internals.

Optional readout modules can participate through a small, versioned protocol:

``REPORT_ADAPTER_API = 1``
``def report_temperature_point(*, params, raw, readout, helper_config): ...``

The callable must return a JSON-compatible mapping.  In particular, an absolute
sensitivity implementation may return ``field_sensitivity_nv_cm_sqrt_hz`` and
``psn_limit_nv_cm_sqrt_hz``.  Unknown or absent helpers remain ``PENDING``; the
static IF discriminator is never relabelled as an absolute sensitivity.
"""
from __future__ import annotations

import importlib
import importlib.util
import math
import re
from dataclasses import dataclass
from typing import Any, Mapping

import numpy as np

from gabes import constants, observables, schemes
from gabes.rydberg_electrometry import (
    BalancedDetector,
    PhotodiodeChannel,
    RFDipoleCoupling,
    balanced_detector_noise,
    current_responsivity_from_atomic_phasor,
    electrometry_sensitivity,
    weak_signal_response,
)


ALLOWED_STATUSES = frozenset({
    "MEASURED",
    "PPT",
    "REFERENCE",
    "FITTED",
    "PREDICTED",
    "ASSUMED",
    "PENDING",
})


def validate_status(value: str) -> str:
    """Return a canonical evidence status or raise a useful error."""
    status = str(value).strip().upper()
    if status not in ALLOWED_STATUSES:
        allowed = ", ".join(sorted(ALLOWED_STATUSES))
        raise ValueError(f"unknown evidence status {value!r}; expected one of {allowed}")
    return status


def _finite_or_none(value: Any) -> float | None:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


def _metric_number(metrics: list[dict[str, Any]], label: str) -> float | None:
    """Extract the leading number from a scheme display metric."""
    for metric in metrics:
        if metric.get("label") != label:
            continue
        match = re.search(r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?",
                          str(metric.get("value", "")))
        return _finite_or_none(match.group(0)) if match else None
    return None


@dataclass
class StaticPoint:
    """One model point plus the non-serializable objects helpers may need."""

    result: dict[str, Any]
    params: dict[str, Any]
    raw: dict[str, Any]
    readout: dict[str, Any]


class StaticRydbergAdapter:
    """Thin adapter around the current registered ``rydberg_eit`` scheme."""

    def __init__(self) -> None:
        self.scheme = schemes.get("rydberg_eit")
        required = ("compute", "recommended_defaults", "_readout", "_eit_features")
        missing = [name for name in required if not hasattr(self.scheme, name)]
        if missing:
            raise RuntimeError(
                "the registered Rydberg scheme is incompatible with the report "
                f"adapter; missing {', '.join(missing)}")

    @property
    def identity(self) -> dict[str, str]:
        return {
            "name": self.scheme.name,
            "class": type(self.scheme).__name__,
            "cache_version": str(self.scheme.cache_version),
            "defaults_version": str(self.scheme.defaults_version),
            "adapter_note": (
                "Transmission arrays currently use the scheme's private _readout "
                "API; this dependency is isolated in adapters.py."
            ),
        }

    def parameters(self, view: str, overrides: Mapping[str, Any]) -> dict[str, Any]:
        defaults = self.scheme.recommended_defaults(self.scheme.defaults())
        if not isinstance(defaults, Mapping) or view not in defaults:
            raise ValueError(f"unknown Rydberg view {view!r}")
        params = dict(defaults[view])
        params.update(dict(overrides))
        params["view"] = view
        return params

    def temperature_point(
        self,
        temperature_c: float,
        view: str,
        overrides: Mapping[str, Any],
        include_spectrum: bool = True,
        *,
        heater_setpoint_c: float | None = None,
        cold_spot_temp_c: float | None = None,
    ) -> StaticPoint:
        params = self.parameters(view, overrides)
        if heater_setpoint_c is None and cold_spot_temp_c is None:
            params["temperature_model"] = "Linked"
            params["temp_c"] = float(temperature_c)
        else:
            params["temperature_model"] = "Separated"
            params["heater_setpoint_c"] = float(
                temperature_c if heater_setpoint_c is None else heater_setpoint_c)
            params["effective_temp_c"] = float(temperature_c)
            params["cold_spot_temp_c"] = float(
                temperature_c if cold_spot_temp_c is None else cold_spot_temp_c)
        raw = self.scheme.compute(params)
        readout = self.scheme._readout(raw, params)
        x = np.asarray(readout["x"], dtype=float)
        transmission = np.asarray(readout["T_trans"], dtype=float)
        width, contrast = self.scheme._eit_features(x, transmission)
        slope = float(np.nanmax(np.abs(np.gradient(transmission, x))))
        ic = int(np.argmin(np.abs(x)))
        metrics = readout.get("metrics", [])

        result: dict[str, Any] = {
            "temperature_c": float(temperature_c),
            "view": view,
            "status": "PREDICTED",
            "parameters": params,
            "metrics": {
                "number_density_m3": {
                    "value": _finite_or_none(raw.get("N")),
                    "unit": "m^-3",
                    "status": "PREDICTED",
                },
                "eit_linewidth_mhz": {
                    "value": _finite_or_none(width),
                    "unit": "MHz",
                    "status": "PREDICTED",
                },
                "eit_peak_contrast": {
                    "value": _finite_or_none(contrast),
                    "unit": "1",
                    "status": "PREDICTED",
                },
                "max_spectral_slope_per_mhz": {
                    "value": _finite_or_none(slope),
                    "unit": "MHz^-1",
                    "status": "PREDICTED",
                },
                "if_discriminator_per_mhz": {
                    "value": _metric_number(metrics, "IF discriminator"),
                    "unit": "MHz^-1",
                    "status": "PREDICTED",
                    "qualifier": "static finite-difference proxy; not absolute sensitivity",
                },
                "if_optimum_detuning_mhz": {
                    "value": _metric_number(metrics, "IF optimum detuning"),
                    "unit": "MHz",
                    "status": "PREDICTED",
                },
                "transmission_at_resonance": {
                    "value": _finite_or_none(transmission[ic]),
                    "unit": "1",
                    "status": "PREDICTED",
                },
                "at_splitting_mhz": {
                    "value": _metric_number(metrics, "RF AT splitting"),
                    "unit": "MHz",
                    "status": "PREDICTED",
                },
                "at_inferred_field_v_m": {
                    "value": _finite_or_none(readout.get("at_field_v_m")),
                    "unit": "V/m",
                    "status": "PREDICTED",
                    "qualifier": (
                        "converted from the full-spectrum model RF Rabi parameter; "
                        "conditional on the declared RF dipole/angular factor"
                    ),
                },
                "at_split_field_estimate_v_m": {
                    "value": _finite_or_none(
                        readout.get("at_split_field_estimate_v_m")),
                    "unit": "V/m",
                    "status": "PREDICTED",
                    "qualifier": (
                        "ideal generalized two-level split estimate; full-spectrum "
                        "fit is preferred because of optical line pulling"
                    ),
                },
            },
        }
        if include_spectrum:
            result["spectrum"] = {
                "detuning_mhz": x.tolist(),
                "transmission": transmission.tolist(),
                "status": "PREDICTED",
            }
        return StaticPoint(result=result, params=params, raw=raw, readout=readout)


class ElectrometryAdapter:
    """Compose finite-IF atomic response, RF coupling, and detector noise.

    No transition dipole or detector calibration is supplied implicitly.  A
    disabled/incomplete config therefore remains ``PENDING``.  When enabled,
    the adapter searches probe detuning and reports the minimum predicted total
    noise-equivalent field together with PSN and technical contributions.
    """

    def __init__(self, static: StaticRydbergAdapter, config: Mapping[str, Any]) -> None:
        self.static = static
        self.config = dict(config)
        self.enabled = bool(self.config.get("enabled", False))
        self.mode = str(self.config.get("mode", "configured")).strip().lower()
        if self.mode not in {"configured", "scheme"}:
            raise ValueError("electrometry.mode must be 'configured' or 'scheme'")
        self.input_status = validate_status(self.config.get("status", "PENDING"))

    def summary(self) -> dict[str, Any]:
        if not self.enabled:
            return {
                "available": False,
                "status": "PENDING",
                "message": str(self.config.get(
                    "note", "absolute electrometry calibration is disabled")),
            }
        if self.mode == "scheme":
            return {
                "available": True,
                "status": "PREDICTED",
                "input_status": self.input_status,
                "mode": "scheme",
                "message": (
                    "using the scheme-native finite-IF response and declared "
                    "RF/detector parameters"
                ),
            }
        required = (
            "transition_dipole_c_m", "angular_factor", "field_amplitude_convention",
            "quantum_efficiency", "signal_probe_power_scale",
        )
        missing = [name for name in required if self.config.get(name) is None]
        if missing:
            return {
                "available": False,
                "status": "PENDING",
                "message": f"missing electrometry calibration fields: {', '.join(missing)}",
            }
        return {
            "available": True,
            "status": "PREDICTED",
            "input_status": self.input_status,
            "mode": "configured",
            "message": "finite-IF linear response and detector/noise chain ready",
        }

    def _detector(self, transmitted_probe_power_w: float) -> BalancedDetector:
        cfg = self.config
        signal = PhotodiodeChannel.from_quantum_efficiency(
            transmitted_probe_power_w,
            float(cfg["quantum_efficiency"]),
            float(cfg.get("probe_wavelength_m", 780.0e-9)),
            dark_current_a=float(cfg.get("signal_dark_current_a", 0.0)),
        )
        reference_power = cfg.get("reference_optical_power_w")
        reference = None
        if reference_power is not None:
            reference = PhotodiodeChannel.from_quantum_efficiency(
                float(reference_power),
                float(cfg.get("reference_quantum_efficiency", cfg["quantum_efficiency"])),
                float(cfg.get("probe_wavelength_m", 780.0e-9)),
                dark_current_a=float(cfg.get("reference_dark_current_a", 0.0)),
            )
        return BalancedDetector(
            signal=signal,
            reference=reference,
            reference_weight=float(cfg.get("reference_weight", 1.0)),
            electronic_noise_current_asd_a_per_sqrt_hz=float(
                cfg.get("electronic_noise_current_asd_a_per_sqrt_hz", 0.0)),
            relative_intensity_noise_per_sqrt_hz=float(
                cfg.get("relative_intensity_noise_per_sqrt_hz", 0.0)),
            rin_correlation=float(cfg.get("rin_correlation", 1.0)),
        )

    def evaluate(self, point: StaticPoint) -> dict[str, Any]:
        state = self.summary()
        if not state["available"]:
            return state
        if point.params.get("view") != "AT electrometry" or point.raw["lo_rabi_mhz"] <= 0:
            return {
                "available": False,
                "status": "PENDING",
                "message": "finite-IF absolute sensitivity requires the LO-dressed AT view",
            }
        if self.mode == "scheme":
            superhet = point.readout.get("superhet")
            if superhet is None:
                return {
                    "available": False,
                    "status": "PENDING",
                    "message": "scheme-native finite-IF readout was not produced",
                }
            sensitivity = superhet["sensitivity"]
            rf_coupling = superhet["rf_coupling"]
            return {
                "available": True,
                "status": "PREDICTED",
                "input_status": self.input_status,
                "method": "scheme-native finite-IF weak-signal Liouvillian response",
                "if_hz": float(point.params.get("if_khz", 40.0)) * 1.0e3,
                "field_amplitude_convention": rf_coupling.field_amplitude_convention,
                "effective_rf_dipole_c_m": rf_coupling.effective_dipole_c_m,
                "superheterodyne_responsivity": superhet[
                    "current_responsivity_a_per_v_m"],
                "field_sensitivity_nv_cm_sqrt_hz": (
                    sensitivity.total_field_asd_nv_cm_per_sqrt_hz),
                "psn_limit_nv_cm_sqrt_hz": (
                    sensitivity.psn_field_asd_nv_cm_per_sqrt_hz),
                "technical_limit_nv_cm_sqrt_hz": (
                    sensitivity.technical_field_asd_nv_cm_per_sqrt_hz),
                "optimum_probe_detuning_mhz": superhet["optimum_detuning_mhz"],
                "optimization_scope": superhet["optimization_scope"],
                "measurement_enbw_hz": superhet["enbw_hz"],
            }
        if point.params.get("doppler", "off") != "off":
            return {
                "available": False,
                "status": "PENDING",
                "message": (
                    "the report adapter does not yet coherently average finite-IF "
                    "phasors over Doppler classes"
                ),
            }

        cfg = self.config
        scheme = self.static.scheme
        mhz = 2.0 * np.pi * 1.0e6
        raw = point.raw
        params = point.params
        ground_dephasing = (
            raw["rydberg_dephasing_mhz"]
            + raw["temperature_dephasing_mhz"]
            + raw["density_dephasing_mhz"]
            + raw["transit_mhz"]
        ) * mhz
        atom = scheme._atom(ground_dephasing, raw["rf_dephasing_mhz"] * mhz)
        probe = raw["probe_rabi_mhz"] * mhz
        coupling = raw["coupling_rabi_mhz"] * mhz
        lo = raw["lo_rabi_mhz"] * mhz
        microwave_detuning = float(params.get("mw_detuning_mhz", 0.0)) * mhz
        if_omega = 2.0 * np.pi * float(params.get("if_khz", 40.0)) * 1.0e3
        rf_coupling = RFDipoleCoupling(
            float(cfg["transition_dipole_c_m"]),
            angular_factor=float(cfg["angular_factor"]),
            field_amplitude_convention=str(cfg["field_amplitude_convention"]),
        )

        x_all = np.asarray(point.readout["x"], dtype=float)
        max_points = max(3, int(cfg.get("probe_detuning_search_points", 161)))
        if x_all.size > max_points:
            detunings = np.linspace(float(x_all[0]), float(x_all[-1]), max_points)
        else:
            detunings = x_all
        transmission = np.asarray(point.readout["T_trans"], dtype=float)
        input_probe_power_w = float(params["probe_power_uw"]) * 1.0e-6
        signal_power_scale = float(cfg["signal_probe_power_scale"])

        # Hermitian observable with Tr(O rho) = Im(rho_10).
        absorption_operator = np.zeros((4, 4), dtype=complex)
        absorption_operator[0, 1] = -0.5j
        absorption_operator[1, 0] = +0.5j
        susceptibility_scale = (
            -2.0 * raw["N"] * raw["ls"] * raw["dipole"] ** 2
            / (constants.EPS_0 * constants.HBAR)
        )
        propagation_scale = (
            -raw["L"] * raw["k_vec"] * susceptibility_scale / probe
        )

        records: list[dict[str, Any]] = []
        for detuning_mhz in detunings:
            s = float(detuning_mhz) * mhz
            hamiltonian = np.zeros((4, 4), dtype=complex)
            hamiltonian[1, 1] = -s
            hamiltonian[2, 2] = -s
            hamiltonian[3, 3] = -s - microwave_detuning
            hamiltonian[0, 1] = hamiltonian[1, 0] = probe / 2.0
            hamiltonian[1, 2] = hamiltonian[2, 1] = coupling / 2.0
            hamiltonian[2, 3] = hamiltonian[3, 2] = lo / 2.0
            response = weak_signal_response(
                atom,
                hamiltonian,
                if_omega,
                signal_transition=(2, 3),
                signal_phase_rad=float(cfg.get("signal_phase_rad", 0.0)),
            )
            atomic_phasor = response.real_observable_phasor_per_angular_rabi(
                absorption_operator)
            t0 = float(np.interp(detuning_mhz, x_all, transmission))
            transmitted_power = input_probe_power_w * signal_power_scale * t0
            detector = self._detector(transmitted_power)
            current_per_atomic = (
                detector.signal.responsivity_a_per_w
                * input_probe_power_w
                * signal_power_scale
                * t0
                * propagation_scale
            )
            current_phasor = current_responsivity_from_atomic_phasor(
                atomic_phasor, current_per_atomic, rf_coupling)
            responsivity = abs(current_phasor)
            if not math.isfinite(responsivity) or responsivity <= 0.0:
                continue
            noise = balanced_detector_noise(detector)
            sensitivity = electrometry_sensitivity(noise, responsivity)
            records.append({
                "probe_detuning_mhz": float(detuning_mhz),
                "transmission": t0,
                "current_responsivity_a_per_v_m": sensitivity.current_responsivity_a_per_v_m,
                "psn_field_asd_nv_cm_per_sqrt_hz": (
                    sensitivity.psn_field_asd_nv_cm_per_sqrt_hz),
                "technical_field_asd_nv_cm_per_sqrt_hz": (
                    sensitivity.technical_field_asd_nv_cm_per_sqrt_hz),
                "total_field_asd_nv_cm_per_sqrt_hz": (
                    sensitivity.total_field_asd_nv_cm_per_sqrt_hz),
                "atomic_phasor_per_angular_rabi_real": float(np.real(atomic_phasor)),
                "atomic_phasor_per_angular_rabi_imag": float(np.imag(atomic_phasor)),
            })
        if not records:
            return {
                "available": False,
                "status": "PENDING",
                "message": "finite-IF scan produced no non-zero calibrated responsivity",
            }
        optimum = min(records, key=lambda row: row["total_field_asd_nv_cm_per_sqrt_hz"])
        return {
            "available": True,
            "status": "PREDICTED",
            "input_status": self.input_status,
            "method": "finite-IF weak-signal Liouvillian response",
            "if_hz": if_omega / (2.0 * np.pi),
            "field_amplitude_convention": rf_coupling.field_amplitude_convention,
            "effective_rf_dipole_c_m": rf_coupling.effective_dipole_c_m,
            "superheterodyne_responsivity": optimum[
                "current_responsivity_a_per_v_m"],
            "field_sensitivity_nv_cm_sqrt_hz": optimum[
                "total_field_asd_nv_cm_per_sqrt_hz"],
            "psn_limit_nv_cm_sqrt_hz": optimum[
                "psn_field_asd_nv_cm_per_sqrt_hz"],
            "technical_limit_nv_cm_sqrt_hz": optimum[
                "technical_field_asd_nv_cm_per_sqrt_hz"],
            "optimum_probe_detuning_mhz": optimum["probe_detuning_mhz"],
            "detuning_scan": records,
        }


class OptionalHelperAdapter:
    """Discover and invoke a report-aware optional readout helper safely."""

    API_VERSION = 1
    DEFAULT_CALLABLE = "report_temperature_point"

    def __init__(self, spec: Mapping[str, Any]) -> None:
        self.spec = dict(spec)
        self.module_name = str(self.spec.get("module", "")).strip()
        self.callable_name = str(
            self.spec.get("callable", self.DEFAULT_CALLABLE)).strip()
        self.required = bool(self.spec.get("required", False))
        if not self.module_name:
            raise ValueError("optional helper requires a non-empty 'module'")
        self._module: Any | None = None
        self._callable: Any | None = None
        self.state = "PENDING"
        self.message = "not discovered"

    def discover(self) -> dict[str, Any]:
        if importlib.util.find_spec(self.module_name) is None:
            self.message = "module not installed"
            if self.required:
                raise RuntimeError(f"required helper module {self.module_name!r} not found")
            return self.summary()
        module = importlib.import_module(self.module_name)
        version = getattr(module, "REPORT_ADAPTER_API", None)
        if version != self.API_VERSION:
            self.message = (
                f"module found but REPORT_ADAPTER_API={version!r}; "
                f"version {self.API_VERSION} required"
            )
            if self.required:
                raise RuntimeError(self.message)
            return self.summary()
        target = getattr(module, self.callable_name, None)
        if not callable(target):
            self.message = f"module found but callable {self.callable_name!r} is absent"
            if self.required:
                raise RuntimeError(self.message)
            return self.summary()
        self._module = module
        self._callable = target
        self.state = "AVAILABLE"
        self.message = "versioned report adapter ready"
        return self.summary()

    def summary(self) -> dict[str, Any]:
        return {
            "module": self.module_name,
            "callable": self.callable_name,
            "required": self.required,
            "state": self.state,
            "status": "PREDICTED" if self.state == "AVAILABLE" else "PENDING",
            "message": self.message,
        }

    def evaluate(self, point: StaticPoint) -> dict[str, Any]:
        if self._callable is None:
            return {"status": "PENDING", "message": self.message}
        result = self._callable(
            params=dict(point.params),
            raw=point.raw,
            readout=point.readout,
            helper_config=dict(self.spec.get("config", {})),
        )
        if not isinstance(result, Mapping):
            raise TypeError(
                f"{self.module_name}.{self.callable_name} must return a mapping")
        normalized = dict(result)
        normalized["status"] = validate_status(normalized.get("status", "PREDICTED"))
        return normalized
