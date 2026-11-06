"""OD/SAS pump-power comparisons from the native optical readout."""

from itertools import combinations
import re

import numpy as np

from .. import __version__
from ..model_comparison import SCHEMA, condition_views, same_grid, sample_values


# The solver emits one marker per allowed transition as "<isotope> <Fg>→<Fe>′".
# A generic manifold has model detunings but no named isotope/hyperfine groups.
_TRANSITION = re.compile(r"(?P<isotope>.+?) (?P<ground>\d+(?:\.\d+)?)→(?P<excited>\d+(?:\.\d+)?)′")


def _transition_groups(raw, frequency_mhz):
    by_ground = {}
    for x_ghz, label in raw.get("markers", ()):
        match = _TRANSITION.fullmatch(label)
        if match is None:
            continue
        key = match["isotope"], match["ground"]
        by_ground.setdefault(key, []).append((float(x_ghz) * 1000, match["excited"]))

    groups = []
    step = float(frequency_mhz[1] - frequency_mhz[0])
    for (isotope, ground), lines in by_ground.items():
        lines.sort()
        first, last = lines[0][0], lines[-1][0]
        margin = max((last - first) * 0.22, 4 * step)
        window = [max(first - margin, float(frequency_mhz[0])),
                  min(last + margin, float(frequency_mhz[-1]))]
        if window[0] >= window[1]:
            continue
        markers = [
            {"x_mhz": x, "label": f"{ground}→{excited}′", "kind": "transition"}
            for x, excited in lines
        ]
        markers.extend(
            {"x_mhz": (x1 + x2) / 2, "kind": "crossover_candidate",
             "label": f"CO candidate ({e1}′, {e2}′)"}
            for (x1, e1), (x2, e2) in combinations(lines, 2) if x1 != x2
        )
        groups.append({"id": f"{isotope}-F{ground}", "label": f"{isotope} F={ground}",
                       "window_mhz": window,
                       "markers": sorted(markers, key=lambda m: m["x_mhz"])})
    groups.sort(key=lambda g: g["window_mhz"][0])
    for group in groups:
        lo, hi = group["window_mhz"]
        group["overlaps"] = [
            other["id"] for other in groups if other is not group
            and max(lo, other["window_mhz"][0]) < min(hi, other["window_mhz"][1])
        ]
    return groups


def pump_power_comparison(scheme, params, *, resolve_raw, current_raw=None, current_view=None):
    """Compare the entered pump power with the same cell and probe at zero pump power."""
    key = scheme.model_comparison_parameter
    current = float(params[key])
    if not np.isfinite(current) or current < 0:
        raise ValueError("Pump power must be finite and non-negative.")
    powers = [0., current] if current else [0.]
    conditions = condition_views(scheme, params, key, powers, resolve_raw=resolve_raw,
                                 current_raw=current_raw, current_view=current_view)
    series = []
    for condition, raw, view in conditions:
        readout = view["optical_readout"]
        if "resolution" not in readout:
            # A running app may have cached a pre-comparison optical view.
            readout = scheme.headless_observables(raw, condition)["optical_readout"]
        x = np.asarray(readout["x"], dtype=float)
        transmission = np.asarray(readout["transmission"], dtype=float)
        optical_depth = np.asarray(readout["optical_depth"], dtype=float)
        if transmission.shape != x.shape or optical_depth.shape != x.shape or x.ndim != 1 or len(x) < 2:
            raise ValueError("Incomplete OD/SAS optical samples for model comparison.")
        if not series:
            native_unit = readout["x_unit"]
            reference = (f"{raw['frequency_reference']} {params['line']} line centroid"
                         if raw["mode"] == "species" else "Generic model zero detuning")
            scale = {"GHz": 1e3, "MHz": 1.0}[native_unit]
            native_axis = x
            if not np.all(np.isfinite(x)) or not np.all(np.diff(x) > 0):
                raise ValueError("OD/SAS requires a finite, increasing frequency axis.")
            frequency_mhz = x * scale
            groups = _transition_groups(raw, frequency_mhz)
            baseline = transmission
            mode = raw["mode"]
            composition = raw.get("isotope_fractions", {})
        else:
            if readout["x_unit"] != native_unit:
                raise ValueError("Inconsistent native frequency units in OD/SAS comparison.")
            same_grid(native_axis, x)
        probe_scope = next((m for m in view["metrics"]
                            if m.get("label") == "Probe model scope"), {})
        series.append({
            "power_mw": condition[key],
            "transmission": sample_values(transmission),
            "optical_depth": sample_values(optical_depth),
            "effective": bool(readout["effective"]),
            "delta": sample_values(transmission - baseline),
            "resolution": readout["resolution"],
            "attenuation_scope": readout["attenuation_scope"],
            "probe_model_scope": probe_scope.get("value"),
            "probe_model_note": probe_scope.get("help"),
        })
    return {
        "schema": SCHEMA, "scheme": scheme.name,
        "x": frequency_mhz.tolist(),
        "axis": {"label": "Relative frequency", "unit": "MHz",
                 "reference": reference, "overview_unit": native_unit},
        "series": series, "groups": groups,
        "default_power_mw": current,
        "fixed": {k: v for k, v in params.items() if k != key},
        "varied": {"key": key, "unit": "mW", "values": powers},
        "observables": {
            "transmission": {"label": "Transmission", "unit": "1"},
            "optical_depth": {"label": "Optical depth", "symbol": "τ", "unit": "1",
                              "definition": readout["conventions"]["optical_depth"]},
            "delta": {"label": "Transmission difference", "unit": "1",
                      "definition": "T(pump on) - T(0 mW) at the same sample"},
        },
        "source": {
            "model": type(scheme).__name__, "kind": mode,
            "package_version": __version__,
            "method": "compute + optical_readout", "cache_version": scheme.cache_version,
            "native_x_unit": native_unit, "native_frequency_reference": reference,
            "native_to_mhz": scale, "samples_per_condition": len(native_axis),
            "isotope_fractions": composition, "interpolation": "none",
            "provenance": "Numerical model", "calibration": "none",
            "baseline_power_mw": 0.0, "missing_value": "null",
            "crossover_basis": "Midpoint of allowed lines from the same ground F; not a detected peak.",
            "group_scope": "Group windows show total transmission, including other transitions in overlapping windows.",
        },
    }
