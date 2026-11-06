"""Lambda view declarations over existing native solves and readouts."""
import numpy as np

from .. import __version__
from ..model_comparison import condition_views, sample_values
from ..lineshape import level_crossing

_COMPARE_KEYS = {"Coupling power": "coupling_power_mw",
                 "Ground relaxation": "buffer_ground_relax_khz",
                 "Beam angle": "beam_angle_mrad", "Detuning map": "coupling_detuning_mhz"}


def comparison_key(params):
    return _COMPARE_KEYS.get(params.get("lambda_compare"), "coupling_power_mw")


def _powers(current, bounds):
    """Three in-domain powers including the exact current value."""
    if not np.isfinite(current) or current < 0:
        raise ValueError("Coupling power must be finite and non-negative.")
    low, high = bounds
    if current == 0:
        return [0., float(np.clip(.25, low, high)), float(np.clip(1., low, high))]
    values = {current}
    for value in (current / 4, current * 4, current / 16, current * 16, .25, 1.):
        values.add(float(np.clip(value, low, high)))
        if len(values) == 3:
            break
    return sorted(values)


def _conditions(scheme, params, key):
    spec = next(sp for sp in scheme.param_schema() if sp.name == key)
    current = float(params[key])
    if key == "coupling_power_mw":
        values = _powers(current, (spec.vmin, spec.vmax))
        if scheme._mode(params) != "cpt" and values[0] > 0:
            values = sorted(set(np.geomspace(values[0], values[-1], 9).tolist() + [current]))
    elif key == "coupling_detuning_mhz":
        half = max(2., min(20., float(params.get("coupling_rabi_mhz", 3.))))
        values = sorted(set(np.linspace(max(spec.vmin, current-half), min(spec.vmax, current+half), 9).tolist()+[current]))
    elif key == "buffer_ground_relax_khz":
        values = sorted({current, *[float(np.clip(v, spec.vmin, spec.vmax))
                                   for v in (current/4, current*4, 1., 10.)]})
    else:
        values = sorted({0., 1., 5., 10., current})
    # Uniform/log grids can round their central anchor by an ulp. Keep one
    # exact current condition instead of an almost-identical extra band.
    values = sorted({current if np.isclose(v, current, rtol=1e-13, atol=1e-14) else float(v) for v in values})
    return spec, values


def _readouts(view, readout, mode):
    metrics = {m["label"]: m for m in view["metrics"]}
    central = metrics["Transmission at center" if mode == "at" else "Transmission at resonance"]
    rows = [dict(key="coupling_rabi", label="Effective Ωc/2π", value=readout["coupling_rabi_mhz"], unit="MHz"),
            dict(key="central_transmission", label="Resonance transmission", value=central["numeric_value"], unit="1")]
    if mode == "at":
        metric = metrics["AT splitting"]
        rows.append(dict(key="doublet_spacing", label="Absorption-peak spacing", value=metric["numeric_value"],
                         unit="MHz", status="resolved" if metric["numeric_value"] is not None else "unresolved",
                         reason="Two absorption maxima (transmission minima); central contrast must exceed 1%."))
        return rows
    width_metric = metrics["Window FWHM"]
    width = width_metric["numeric_value"]
    scale = 1000. if width_metric["unit"] == "kHz" else 1.
    x, transmission = np.asarray(readout["x_mhz"]), np.asarray(readout["transmission"])
    center = int(np.argmin(np.abs(x-readout["center_mhz"])))
    floor = float(np.nanmin(transmission))
    contrast = float(transmission[center]-floor)
    half = (transmission[center]+floor)/2
    left = right = np.nan
    if np.isfinite(half) and contrast > 0:
        left, _ = level_crossing(x, transmission, center, -1, half)
        right, _ = level_crossing(x, transmission, center, 1, half)
    bracketed = np.isfinite([left, right]).all() and contrast > 0
    intervals = width/scale/float(np.median(np.diff(x))) if width is not None else None
    status = ("unresolved" if width is None else "scan-edge-limited" if not bracketed else
              "resolution-limited" if intervals < 6 else "resolved")
    # Keep the original metric untouched; the comparison withholds clipped widths.
    shown_width = width if bracketed else None
    rows += [dict(key="contrast", label="Transparency contrast ΔT", value=contrast if np.isfinite(contrast) else None,
                  unit="1", reason="T(δ = 0) − min T of this complete native scan; absolute transmission difference."),
             dict(key="window_width", label="Transparency-window FWHM", value=shown_width,
                  unit=width_metric["unit"], status=status,
                  reason="Both half-height crossings must lie inside the scan. Fewer than six native intervals is approximate."),
             dict(key="contrast_per_width", label="Contrast / FWHM",
                  value=contrast/shown_width if shown_width is not None and shown_width > 0 else None,
                  unit=width_metric["unit"]+"⁻¹", status=status,
                  reason="Spectral quality indicator; not clock stability or Allan deviation."),
             dict(key="samples_per_fwhm", label="Samples per FWHM", value=intervals, unit="1")]
    return rows


def _summary(curves, mode, key, spec):
    effective = [next(r["value"] for r in c["readouts"] if r["key"] == "coupling_rabi") for c in curves]
    rabi_axis = mode == "at" and key == "coupling_power_mw"
    axis = (dict(label="Effective Ωc/2π", unit="MHz", reference="Fixed beam diameter and Rabi anchor") if rabi_axis else
            dict(label=spec.label, unit=spec.unit, reference="Actual varied control; all other conditions fixed"))
    xs = effective if rabi_axis else [c["parameter"]["value"] for c in curves]
    keys = ["doublet_spacing"] if mode == "at" else ["contrast", "window_width"] + (["contrast_per_width"] if mode == "cpt" else [])
    series, panels = [], []
    for i, key_name in enumerate(keys):
        entries = [next(r for r in c["readouts"] if r["key"] == key_name) for c in curves]
        entry = entries[0]
        sid = "extracted" if i == 0 else key_name
        series.append(dict(id=sid, label=entry["label"], x=xs, observables=dict(value=[r["value"] for r in entries]),
                           sample_status=[r.get("status", "available") if r.get("status") != "resolved" else "available" for r in entries],
                           condition_ids=[c["id"] for c in curves]))
        panels.append(dict(key=key_name, title=entry["label"], y_label=entry["label"], unit=entry["unit"],
                           note=entry.get("reason"), curves=[dict(series_id=sid, observable="value")]))
    if mode == "at":
        series.append(dict(id="input", label="Input Ωc/2π reference", x=xs, observables=dict(value=effective),
                           condition_ids=[c["id"] for c in curves]))
        panels[0]["curves"].append(dict(series_id="input", observable="value"))
    return dict(axis=axis, series=series, panels=panels)


def coupling_comparison(scheme, params, *, resolve_raw, current_raw=None, current_view=None):
    key, mode = comparison_key(params), scheme._mode(params)
    detuning_map = key == "coupling_detuning_mhz"
    unit = "kHz" if mode == "cpt" and not detuning_map else "MHz"
    axis = dict(label="Two-photon detuning Δp − Δc" if unit == "kHz" else "Probe detuning Δp",
                unit=unit, reference="Two-photon resonance" if unit == "kHz" else "Bare probe transition")
    spec, values = _conditions(scheme, params, key)
    curves = []
    for i, (condition, raw, view) in enumerate(condition_views(
            scheme, params, key, values, resolve_raw=resolve_raw,
            current_raw=current_raw, current_view=current_view)):
        if "lambda_readout" not in view:
            view = scheme.headless_observables(raw, condition)
        readout = view["lambda_readout"]
        x = np.asarray(readout["x_mhz"])
        shown_x = (x-readout["center_mhz"])*1000 if unit == "kHz" else x
        value = float(condition[key])
        rabi = float(readout["coupling_rabi_mhz"])
        label = f"{value:.6g} {spec.unit}" + (f" · Ωc/2π {rabi:.5g} MHz" if key == "coupling_power_mw" else "")
        metrics = {m["label"]: m for m in view["metrics"]}
        markers = [dict(x=0. if unit == "kHz" else readout["center_mhz"], label="δ = 0")]
        if mode == "at":
            markers += [dict(x=float(position), label=f"Absorption peak {j+1}")
                        for j, position in enumerate(metrics["AT splitting"].get("sample_positions", []))]
        markers.sort(key=lambda marker: marker["x"])
        curves.append(dict(id=f"condition-{i}", label=label,
                           parameter=dict(key=key, label=spec.label, value=value, unit=spec.unit),
                           coordinates=dict(coupling_rabi=rabi), markers=markers,
                           x=sample_values(shown_x),
                           observables=dict(transmission=sample_values(readout["transmission"]),
                                            optical_depth=sample_values(np.asarray(readout["absorption_per_m"])*float(condition["cell_mm"])*1e-3),
                                            dispersion=sample_values(readout["dispersion"])),
                           readouts=_readouts(view, readout, mode)))
    panels = [dict(key="transmission", title="Transmission", y_label="Transmission T", unit="1", cut_group="intensity"),
              dict(key="optical_depth", title="Optical depth", y_label="Optical depth OD", unit="1", cut_group="intensity"),
              dict(key="dispersion", title="Dispersion from the same solve", y_label="Re χ", unit="1", reference=0)]
    for panel in panels:
        panel["curves"] = [dict(series_id=c["id"], observable=panel["key"]) for c in curves]
    summary = _summary(curves, mode, key, spec)
    display = dict(default_mode="curves" if mode == "cpt" and not detuning_map else "map")
    if mode in {"eit", "cpt"} and not detuning_map:
        selected = next(c for c in curves if c["parameter"]["value"] == params[key])
        width = next(r["value"] for r in selected["readouts"] if r["key"] == "window_width")
        half = max((width or 0)*2, 10*float(np.median(np.diff(selected["x"]))))
        center = 0. if unit == "kHz" else float(params["coupling_detuning_mhz"])
        display["window"] = [center-half, center+half]
    if detuning_map:
        detunings = [c["parameter"]["value"] for c in curves]
        margin = max(2., (max(detunings)-min(detunings))*.2)
        if mode == "at":
            margin = max(margin, max((abs(m["x"]-c["parameter"]["value"]) for c in curves
                                     for m in c["markers"]), default=0.))
        display["window"] = [min(detunings)-margin, max(detunings)+margin]
    result = dict(schema="gabes-spectral-comparison-v1", scheme=scheme.name, mode=mode, axis=axis,
                  title=mode.upper()+" · "+("Detuning map" if detuning_map else spec.label+" comparison"),
                  series=curves, panels=panels, summary=summary, display=display,
                  fixed={k: v for k, v in params.items() if k != key},
                  varied=dict(key=key, label=spec.label, unit=spec.unit, values=[c["parameter"]["value"] for c in curves]),
                  source=dict(model=type(scheme).__name__, cache_version=scheme.cache_version,
                              package_version=__version__, method="compute + headless lambda_readout",
                              provenance="Numerical model", calibration="none", native_probe_unit="MHz",
                              interpolation="none", default_parameter=float(params[key]),
                              contrast_definition="T at two-photon resonance minus the minimum T of the complete native scan",
                              width_definition="Existing sample-bracketed half-height width; clipped widths withheld"),
                  notes=["Calculated scalar three-level Λ response; a selected regime is a preset, not an EIT/ATS classifier.",
                         "Each condition retains its own 601-point frequency grid. Zoom does not improve native resolution.",
                         "Contrast and FWHM use the minimum transmission of that condition's complete scan as background."])
    if key == "coupling_power_mw":
        # Supplied by the scheme's existing beam conversion, not a new UI physics model.
        medium_gamma = float(current_raw.get("gamma_mhz", 1.)) if current_raw is not None else 1.
        unit_rabi = scheme._coupling_rabi_from_params({**params, key: 1.}, medium_gamma, scheme._DEF[mode]["oc_mhz"])
        result["condition_axes"] = [dict(key="parameter", label="Coupling power", unit="mW"),
                                    dict(key="coupling_rabi", label="Effective Ωc/2π", unit="MHz",
                                         transform=dict(exponent=.5, factor=unit_rabi))]
        display["condition_axis"] = "coupling_rabi" if mode == "at" else "parameter"
    if detuning_map:
        result["map_guides"] = [dict(label="Two-photon resonance Δp = Δc", x_per_y=1., x_offset=0.)]
        result["notes"].append("Δp × Δc: diagonal δ = 0. Rows retain shifted native optical grids; uncovered regions stay blank.")
    return result
