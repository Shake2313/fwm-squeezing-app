"""
85Rb FWM source -> 87Rb D1 Hanle magnetometer reference analysis.

This script connects three pieces that are separate in the Streamlit app:

1. Seeded 85Rb D1 double-Lambda FWM (`gabes.schemes.fwm`).
2. 87Rb D1 Hanle transmission (`gabes.schemes.magneto`).
3. Unequal-arm probe/reference balanced detection
   (`gabes.observables.balanced_twin_beam_noise`).

The frequency lock is computed from the isotope hyperfine constants instead of
assuming the sign convention of a paper's two-photon detuning.  The GABES FWM
probe axis is referenced to the 85Rb D1 F=2 -> F'=3 transition.

Typical use:

    python analysis/squeezing/resonant_hanle_squeezing_reference.py --fidelity fine

Outputs:
    analysis/squeezing/resonant_hanle_squeezing_reference.npz
    analysis/squeezing/resonant_hanle_squeezing_reference.png
    analysis/squeezing/resonant_hanle_squeezing_reference.md
"""
from __future__ import annotations

import argparse
import json
import math
import os
import sys
import time
from dataclasses import dataclass
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from gabes import constants, observables, species  # noqa: E402
from gabes.schemes import fwm, magneto  # noqa: E402
from gabes.core import blas_single_thread  # noqa: E402


OUT_STEM = "resonant_hanle_squeezing_reference"
RB85_REF = (species.RB85, 2, 3)       # GABES FWM frequency origin
RB87_PROBE = (species.RB87, 2, 2)     # 87Rb D1 F=2 -> F'=2
RB87_CONJ = (species.RB87, 1, 1)      # 87Rb D1 F=1 -> F'=1
PLANCK = 2.0 * math.pi * constants.HBAR
B_UNIT_TO_UT = {
    "ut": 1.0,
    "microtesla": 1.0,
    "microteslas": 1.0,
    "nt": 1e-3,
    "nanotesla": 1e-3,
    "mt": 1e3,
    "millitesla": 1e3,
    "t": 1e6,
    "tesla": 1e6,
    "g": 100.0,
    "gauss": 100.0,
    "mg": 0.1,
    "milligauss": 0.1,
}


@dataclass(frozen=True)
class SourceConfig:
    name: str
    label: str
    D_GHz: float
    delta_code_mhz: float
    temp_c: float
    pump_w: float
    seed_w: float
    cell_m: float
    w_pump_m: float
    w_probe_m: float
    angle_deg: float
    line_strength: float
    loss_frac: float
    qe: float
    target: str
    paper_note: str
    anchor_gain: float | None = None
    anchor_source_s_db: float | None = None
    anchor_direct_s_db: float | None = None


def transition_frequency_hz(iso: species.Isotope, Fg: int, Fe: int,
                            line: str = "D1") -> float:
    Je, nu0, _gamma_mhz, A_e, B_e = iso.line(line)
    Eg = species.hf_energy_mhz(iso.A_S, 0.0, iso.I, iso.Jg, Fg)
    Ee = species.hf_energy_mhz(A_e, B_e, iso.I, Je, Fe)
    return nu0 + (Ee - Eg) * 1e6


def offset_from_fwm_origin_ghz(iso: species.Isotope, Fg: int, Fe: int) -> float:
    ref_iso, ref_Fg, ref_Fe = RB85_REF
    return (transition_frequency_hz(iso, Fg, Fe)
            - transition_frequency_hz(ref_iso, ref_Fg, ref_Fe)) / 1e9


def source_probe_offset_ghz(D_GHz: float, delta_code_mhz: float,
                            branch: int = -1) -> float:
    return fwm.branch_center_GHz(D_GHz, branch) + delta_code_mhz * 1e-3


def source_conjugate_offset_ghz(D_GHz: float, delta_code_mhz: float,
                                branch: int = -1) -> float:
    probe = source_probe_offset_ghz(D_GHz, delta_code_mhz, branch)
    return 2.0 * D_GHz - probe


def D_for_probe_lock(probe_offset_ghz: float, delta_code_mhz: float,
                     branch: int = -1) -> float:
    # probe = D + branch*nu_HF + delta
    return probe_offset_ghz - branch * constants.NU_HF / 1e9 - delta_code_mhz * 1e-3


def D_for_conjugate_lock(conj_offset_ghz: float, delta_code_mhz: float,
                         branch: int = -1) -> float:
    # conjugate = 2D - probe = D - branch*nu_HF - delta
    return conj_offset_ghz + branch * constants.NU_HF / 1e9 + delta_code_mhz * 1e-3


def D_delta_for_double_lock(probe_offset_ghz: float, conj_offset_ghz: float,
                            branch: int = -1) -> tuple[float, float]:
    D = 0.5 * (probe_offset_ghz + conj_offset_ghz)
    delta_ghz = probe_offset_ghz - D - branch * constants.NU_HF / 1e9
    return D, delta_ghz * 1e3


def make_source_configs() -> list[SourceConfig]:
    rb87_probe = offset_from_fwm_origin_ghz(*RB87_PROBE)
    rb87_conj = offset_from_fwm_origin_ghz(*RB87_CONJ)
    oe_delta_code = -16.0  # OE paper's +16 MHz in the opposite beat-note sign.
    D_probe = D_for_probe_lock(rb87_probe, oe_delta_code)
    D_conj = D_for_conjugate_lock(rb87_conj, oe_delta_code)
    D_double, delta_double = D_delta_for_double_lock(rb87_probe, rb87_conj)
    return [
        SourceConfig(
            name="sim_2025_ids",
            label="Sim 2025 IDS optimum",
            D_GHz=0.9,
            delta_code_mhz=-8.0,
            temp_c=121.0,
            pump_w=0.6,
            seed_w=8e-6,
            cell_m=12.5e-3,
            w_pump_m=530e-6,
            w_probe_m=330e-6,
            angle_deg=0.32,
            line_strength=1.0,
            loss_frac=0.055,
            qe=fwm.QE_DETECTOR,
            target="off_resonant_ids",
            paper_note="Sim et al. 85Rb IDS optimum; not 87Rb-resonant.",
        ),
        SourceConfig(
            name="oe_probe_locked",
            label="OE probe resonant",
            D_GHz=D_probe,
            delta_code_mhz=oe_delta_code,
            temp_c=89.0,
            pump_w=1.0,
            seed_w=8e-6,
            cell_m=12.0e-3,
            w_pump_m=0.95e-3,     # OE gives 1/e^2 diameter = 1.9 mm.
            w_probe_m=0.30e-3,    # OE gives 1/e^2 diameter = 0.6 mm.
            angle_deg=0.45,
            line_strength=1.0,
            loss_frac=0.10,       # approx output-window + optics + detector path.
            qe=0.95,
            target="probe_87rb_F2_Fp2",
            paper_note="Frequency-locked to 87Rb F=2 -> F'=2 probe resonance.",
            anchor_gain=5.6,
            anchor_source_s_db=-8.3,
            anchor_direct_s_db=-5.4,
        ),
        SourceConfig(
            name="oe_conjugate_locked",
            label="OE conjugate resonant",
            D_GHz=D_conj,
            delta_code_mhz=oe_delta_code,
            temp_c=89.0,
            pump_w=1.0,
            seed_w=8e-6,
            cell_m=12.0e-3,
            w_pump_m=0.95e-3,
            w_probe_m=0.30e-3,
            angle_deg=0.45,
            line_strength=1.0,
            loss_frac=0.10,
            qe=0.95,
            target="conjugate_87rb_F1_Fp1",
            paper_note="Frequency-locked to 87Rb F=1 -> F'=1 conjugate resonance.",
            anchor_gain=3.8,
            anchor_source_s_db=-8.1,
            anchor_direct_s_db=-5.0,
        ),
        SourceConfig(
            name="oe_double_locked",
            label="OE double resonant",
            D_GHz=D_double,
            delta_code_mhz=delta_double,
            temp_c=91.0,
            pump_w=1.0,
            seed_w=8e-6,
            cell_m=12.0e-3,
            w_pump_m=0.95e-3,
            w_probe_m=0.30e-3,
            angle_deg=0.50,
            line_strength=1.0,
            loss_frac=0.10,
            qe=0.95,
            target="probe_and_conjugate_87rb",
            paper_note="Both modes frequency-locked to the two 87Rb D1 transitions.",
            anchor_gain=1.9,
            anchor_source_s_db=-4.7,
            anchor_direct_s_db=-3.5,
        ),
    ]


def fidelity_settings(name: str) -> dict:
    name = name.lower()
    if name == "quick":
        return dict(coarse_points=81, velocity_step=8.0,
                    velocity_cutoff=2.5, phase_detail=fwm.PHASE_BALANCED,
                    window_ghz=0.22)
    if name == "balanced":
        return dict(coarse_points=161, velocity_step=5.0,
                    velocity_cutoff=3.0, phase_detail=fwm.PHASE_BALANCED,
                    window_ghz=0.30)
    if name == "ultra":
        return dict(coarse_points=181, velocity_step=4.0,
                    velocity_cutoff=3.0, phase_detail=fwm.PHASE_ULTRA,
                    window_ghz=0.26)
    return dict(coarse_points=241, velocity_step=3.0,
                velocity_cutoff=3.0, phase_detail=fwm.PHASE_FINE,
                window_ghz=0.35)


def parse_float_axis(text: str) -> np.ndarray:
    vals = [float(part.strip()) for part in str(text).split(",")
            if part.strip()]
    if not vals:
        raise ValueError("axis list must contain at least one number")
    return np.asarray(vals, dtype=float)


def args_with(args: argparse.Namespace, **updates) -> argparse.Namespace:
    data = vars(args).copy()
    data.update(updates)
    return argparse.Namespace(**data)


def gaussian_peak_intensity_mw_cm2(power_w: float, waist_m: float) -> float:
    if power_w <= 0.0 or waist_m <= 0.0:
        return 0.0
    return 0.1 * 2.0 * float(power_w) / (math.pi * float(waist_m) ** 2)


def gaussian_waist_um_for_peak_intensity(power_w: float,
                                         intensity_mw_cm2: float) -> float:
    if power_w <= 0.0 or intensity_mw_cm2 <= 0.0:
        return float("nan")
    waist_m = math.sqrt(0.1 * 2.0 * float(power_w)
                        / (math.pi * float(intensity_mw_cm2)))
    return 1e6 * waist_m


def compute_source(cfg: SourceConfig, settings: dict) -> dict:
    probe_target = source_probe_offset_ghz(cfg.D_GHz, cfg.delta_code_mhz)
    half = settings["window_ghz"]
    t0 = time.time()
    with blas_single_thread():
        spec = fwm.compute_spectrum(
            cfg.D_GHz,
            T=cfg.temp_c + 273.15,
            P_pump=cfg.pump_w,
            P_probe=cfg.seed_w,
            w_pump=cfg.w_pump_m,
            w_probe=cfg.w_probe_m,
            line_strength=cfg.line_strength,
            loss_frac=cfg.loss_frac,
            qe=cfg.qe,
            L=cfg.cell_m,
            coarse_points=settings["coarse_points"],
            fine_points=0,
            scan_min=probe_target - half,
            scan_max=probe_target + half,
            velocity_step=settings["velocity_step"],
            velocity_cutoff=settings["velocity_cutoff"],
            branch=-1,
            phase_detail=settings["phase_detail"],
            pump_probe_angle_deg=cfg.angle_deg,
            model_fidelity=settings["phase_detail"],
        )
    op = fwm.operating_point(spec, cfg.delta_code_mhz, branch=-1)
    rb87_probe = offset_from_fwm_origin_ghz(*RB87_PROBE)
    rb87_conj = offset_from_fwm_origin_ghz(*RB87_CONJ)
    probe_offset = float(op["probe_GHz"])
    conj_offset = 2.0 * cfg.D_GHz - probe_offset
    return dict(
        config=cfg,
        spectrum=spec,
        op=op,
        runtime_s=time.time() - t0,
        probe_offset_ghz=probe_offset,
        conjugate_offset_ghz=conj_offset,
        probe_detuning_from_87_mhz=(probe_offset - rb87_probe) * 1e3,
        conjugate_detuning_from_87_mhz=(conj_offset - rb87_conj) * 1e3,
        paper_delta_equiv_mhz=-cfg.delta_code_mhz,
    )


def available_arm_power_w(source: dict, arm: str,
                          args: argparse.Namespace) -> float:
    cfg = source["config"]
    src = source_readout_values(source, args)
    gain = src["G_s"] if arm == "probe" else src["G_c"]
    return max(gain * cfg.seed_w, 0.0)


def effective_pre_atten(source: dict, arm: str,
                        args: argparse.Namespace) -> float:
    base = available_arm_power_w(source, arm, args)
    if arm == "probe":
        requested_uw = float(getattr(args, "probe_input_power_uw", 0.0))
        fallback = float(args.probe_pre_atten)
    else:
        requested_uw = float(getattr(args, "reference_input_power_uw", 0.0))
        fallback = float(args.reference_pre_atten)
    if requested_uw <= 0.0:
        return fallback
    requested = requested_uw * 1e-6
    if requested > base * (1.0 + 1e-9):
        raise ValueError(
            f"requested {arm} input power {requested_uw:g} uW exceeds "
            f"available FWM {arm} power {base * 1e6:g} uW")
    return requested / max(base, 1e-300)


def arm_input_power_w(source: dict, arm: str,
                      args: argparse.Namespace) -> float:
    return available_arm_power_w(source, arm, args) \
        * effective_pre_atten(source, arm, args)


def resolved_hanle_intensity_mw_cm2(
        args: argparse.Namespace, probe_power_w: float | None = None) -> float:
    mode = getattr(args, "hanle_intensity_mode", "explicit")
    if mode == "waist":
        if probe_power_w is None:
            raise ValueError("probe_power_w is required in waist intensity mode")
        waist_um = float(args.hanle_probe_waist_um)
        if waist_um <= 0.0:
            raise ValueError("--hanle-probe-waist-um must be > 0 in waist mode")
        return gaussian_peak_intensity_mw_cm2(probe_power_w, waist_um * 1e-6)
    return float(args.hanle_intensity)


def hanle_params(detuning_mhz: float, args: argparse.Namespace,
                 probe_power_w: float | None = None) -> dict:
    cell_type = (magneto.CELL_BUFFER if args.cell == "buffer"
                 else magneto.CELL_PARAFFIN)
    intensity = resolved_hanle_intensity_mw_cm2(args, probe_power_w)
    params = {
        "signal_type": magneto.SIGNAL_TRANSMISSION,
        "cell_type": cell_type,
        "transition": "F=2 -> F'=2",
        "intensity_mw_cm2": intensity,
        "qwp_deg": args.qwp_deg,
        "b_max_ut": args.b_max_ut,
        "laser_detuning_mhz": detuning_mhz,
        "temp_c": args.hanle_temp_c,
        "cell_mm": args.hanle_cell_mm,
        "ne_pressure_torr": args.ne_torr,
        "wall_coherence_ms": args.wall_ms,
        "buffer_ground_relax_khz": args.buffer_ground_khz,
        "collisional_depol_khz": args.collisional_depol_khz,
        "transit_relax_khz": args.transit_khz,
        "dark_return_khz": args.dark_return_khz,
        "residual_transverse_b_ut": args.residual_b_ut,
        "transverse_field_angle_deg": args.transverse_angle_deg,
        "line_strength": args.hanle_line_strength,
        "doppler": "on",
        "velocity_classes": args.velocity_classes,
        "scan_points": args.scan_points,
    }
    return params


def compute_hanle(detuning_mhz: float, args: argparse.Namespace,
                  probe_power_w: float | None = None) -> dict:
    scheme = magneto.MagnetoScheme()
    params = hanle_params(detuning_mhz, args, probe_power_w)
    raw = scheme.compute(params)
    xprobe = observables.chi_phys(
        raw["chi_probe"], raw["N_eff"], dipole=raw["dipole"],
        line_strength=params["line_strength"])
    b = raw["b_ut"]
    L = params["cell_mm"] * 1e-3
    alpha = raw["k_vec"] * np.imag(xprobe)
    trans = np.exp(-alpha * L)
    slope = np.gradient(trans, b)  # per microtesla.
    finite = np.isfinite(slope) & np.isfinite(trans)
    idx = int(np.nanargmax(np.where(finite, np.abs(slope), np.nan)))
    return dict(params=params, raw=raw, b_ut=b, transmission=trans,
                slope_per_ut=slope, bias_idx=idx)


def _split_csv_line(line: str) -> list[str]:
    if "," in line:
        return [p.strip() for p in line.split(",")]
    if "\t" in line:
        return [p.strip() for p in line.split("\t")]
    return line.split()


def _is_numeric_row(parts: list[str]) -> bool:
    if not parts:
        return False
    try:
        [float(p) for p in parts]
    except ValueError:
        return False
    return True


def _delimiter_from_line(line: str):
    if "," in line:
        return ","
    if "\t" in line:
        return "\t"
    return None


def _resolve_column(names: tuple[str, ...], requested: str,
                    default_index: int) -> str:
    if not names:
        raise ValueError("no named columns available")
    if requested:
        if requested in names:
            return requested
        lowered = {n.lower(): n for n in names}
        key = requested.lower()
        if key in lowered:
            return lowered[key]
        raise ValueError(f"column {requested!r} not found in {names}")
    if default_index >= len(names):
        raise ValueError("not enough columns in measured Hanle CSV")
    return names[default_index]


def _resolve_column_index(n_columns: int, requested: str,
                          default_index: int) -> int:
    if requested:
        idx = int(requested)
    else:
        idx = default_index
    if idx < 0 or idx >= n_columns:
        raise ValueError(f"column index {idx} outside 0..{n_columns - 1}")
    return idx


def load_hanle_csv(path: str | Path, *, b_column: str = "",
                   signal_column: str = "", b_unit: str = "uT") -> dict:
    csv_path = Path(path).expanduser()
    if not csv_path.is_absolute():
        csv_path = (Path.cwd() / csv_path).resolve()
    lines = csv_path.read_text(encoding="utf-8-sig").splitlines()
    first = next((ln for ln in lines
                  if ln.strip() and not ln.lstrip().startswith("#")), "")
    if not first:
        raise ValueError(f"{csv_path} contains no data rows")
    delimiter = _delimiter_from_line(first)
    has_header = not _is_numeric_row(_split_csv_line(first))
    unit_key = str(b_unit).strip().lower().replace("µ", "u").replace("μ", "u")
    if unit_key not in B_UNIT_TO_UT:
        raise ValueError(f"unknown B unit {b_unit!r}")
    factor = B_UNIT_TO_UT[unit_key]

    if has_header:
        data = np.genfromtxt(csv_path, delimiter=delimiter, names=True,
                             dtype=float, comments="#", encoding="utf-8-sig")
        data = np.atleast_1d(data)
        names = data.dtype.names or ()
        b_name = _resolve_column(names, b_column, 0)
        y_name = _resolve_column(names, signal_column, 1)
        b_raw = np.asarray(data[b_name], dtype=float)
        signal = np.asarray(data[y_name], dtype=float)
        b_label, y_label = b_name, y_name
    else:
        data = np.genfromtxt(csv_path, delimiter=delimiter, dtype=float,
                             comments="#", encoding="utf-8-sig")
        data = np.asarray(data, dtype=float)
        if data.ndim == 1:
            data = data[None, :]
        if data.shape[1] < 2:
            raise ValueError("measured Hanle CSV needs at least two columns")
        b_idx = _resolve_column_index(data.shape[1], b_column, 0)
        y_idx = _resolve_column_index(data.shape[1], signal_column, 1)
        b_raw = data[:, b_idx]
        signal = data[:, y_idx]
        b_label, y_label = str(b_idx), str(y_idx)

    b_ut = np.asarray(b_raw, dtype=float) * factor
    signal = np.asarray(signal, dtype=float)
    finite = np.isfinite(b_ut) & np.isfinite(signal)
    if np.count_nonzero(finite) < 3:
        raise ValueError("measured Hanle CSV needs at least three finite rows")
    order = np.argsort(b_ut[finite])
    return dict(
        path=str(csv_path),
        b_ut=b_ut[finite][order],
        signal=signal[finite][order],
        b_column=b_label,
        signal_column=y_label,
        b_unit=b_unit,
        n_points=int(np.count_nonzero(finite)),
    )


def _fit_affine_model(model_y: np.ndarray,
                      measured_y: np.ndarray) -> tuple[float, float, float]:
    A = np.column_stack([np.ones_like(model_y), model_y])
    offset, scale = np.linalg.lstsq(A, measured_y, rcond=None)[0]
    fit = offset + scale * model_y
    rms = float(np.sqrt(np.mean((measured_y - fit) ** 2)))
    return float(offset), float(scale), rms


def fit_hanle_calibration(data_b_ut: np.ndarray, data_signal: np.ndarray,
                          model_b_ut: np.ndarray,
                          model_signal: np.ndarray, *,
                          offset_span_ut: float = 0.5,
                          offset_points: int = 101,
                          scale_span: float = 0.10,
                          scale_points: int = 41) -> dict:
    data_b_ut = np.asarray(data_b_ut, dtype=float)
    data_signal = np.asarray(data_signal, dtype=float)
    model_b_ut = np.asarray(model_b_ut, dtype=float)
    model_signal = np.asarray(model_signal, dtype=float)
    finite = (np.isfinite(data_b_ut) & np.isfinite(data_signal)
              & np.isfinite(model_b_ut).all()
              & np.isfinite(model_signal).all())
    if np.count_nonzero(finite) < 3:
        raise ValueError("not enough finite measured Hanle points")
    data_b_ut = data_b_ut[finite]
    data_signal = data_signal[finite]
    order = np.argsort(model_b_ut)
    model_b_ut = model_b_ut[order]
    model_signal = model_signal[order]
    offsets = np.linspace(-abs(offset_span_ut), abs(offset_span_ut),
                          max(int(offset_points), 1))
    scales = np.linspace(1.0 - abs(scale_span), 1.0 + abs(scale_span),
                         max(int(scale_points), 1))
    best = None
    best_model = None
    for scale in scales:
        for offset in offsets:
            x = scale * data_b_ut + offset
            model_at_data = np.interp(x, model_b_ut, model_signal,
                                      left=np.nan, right=np.nan)
            valid = np.isfinite(model_at_data)
            if np.count_nonzero(valid) < 3:
                continue
            y0, yscale, rms = _fit_affine_model(
                model_at_data[valid], data_signal[valid])
            if best is None or rms < best["rms"]:
                fitted = y0 + yscale * model_at_data
                best = dict(
                    b_scale=float(scale),
                    b_offset_ut=float(offset),
                    signal_offset=float(y0),
                    signal_scale=float(yscale),
                    rms=float(rms),
                    n_fit_points=int(np.count_nonzero(valid)),
                )
                best_model = fitted
    if best is None:
        raise ValueError("measured Hanle data do not overlap model B scan")
    resid = data_signal - best_model
    y_centered = data_signal - float(np.mean(data_signal))
    sst = float(np.sum(y_centered * y_centered))
    sse = float(np.nansum(resid * resid))
    best["r2"] = float(1.0 - sse / sst) if sst > 0.0 else float("nan")
    best["normalized_rms"] = float(best["rms"]
                                   / max(np.ptp(data_signal), 1e-300))
    best["fitted_signal"] = best_model
    best["residual"] = resid
    best["model_b_at_data_ut"] = best["b_scale"] * data_b_ut \
        + best["b_offset_ut"]
    return best


def maybe_calibrate_hanle(results: list[dict],
                          args: argparse.Namespace) -> dict | None:
    if not args.hanle_calibration_csv:
        return None
    data = load_hanle_csv(
        args.hanle_calibration_csv,
        b_column=args.hanle_calibration_b_column,
        signal_column=args.hanle_calibration_signal_column,
        b_unit=args.hanle_calibration_b_unit)
    target = args.hanle_calibration_source
    result = None
    for r in results:
        if r["source"]["config"].name == target:
            result = r
            break
    if result is None:
        usable = [r for r in results if r["balanced"]["shot"]["usable"]]
        result = usable[0] if usable else results[0]
    fit = fit_hanle_calibration(
        data["b_ut"], data["signal"], result["hanle"]["b_ut"],
        result["hanle"]["transmission"],
        offset_span_ut=args.hanle_calibration_offset_span_ut,
        offset_points=args.hanle_calibration_offset_points,
        scale_span=args.hanle_calibration_scale_span,
        scale_points=args.hanle_calibration_scale_points)
    return dict(data=data, fit=fit, source_name=result["source"]["config"].name,
                source_label=result["source"]["config"].label)


def source_readout_values(source: dict, args: argparse.Namespace) -> dict:
    cfg = source["config"]
    use_anchor = (args.source_anchor == "paper"
                  and cfg.anchor_gain is not None
                  and cfg.anchor_source_s_db is not None)
    if use_anchor:
        Gs = float(cfg.anchor_gain)
        Gc = max(Gs - 1.0, 0.0)
        S_dB = float(cfg.anchor_source_s_db)
        S_linear = 10.0 ** (S_dB / 10.0)
        label = "paper_anchor"
    else:
        op = source["op"]
        Gs = float(op["G_s"])
        Gc = float(op["G_c"])
        S_dB = float(op["S_dB"])
        S_linear = None
        label = "gabes_model"
    return dict(G_s=Gs, G_c=Gc, S_dB=S_dB,
                source_noise_linear=S_linear, label=label)


def detector_responsivity_aw(args: argparse.Namespace) -> float:
    if args.responsivity_aw > 0:
        return float(args.responsivity_aw)
    return args.detector_qe * args.probe_wavelength_nm / 1239.841984


def absolute_field_sensitivity(
        S, weight, P_probe0, P_ref0, T, dT_per_ut,
        args: argparse.Namespace, *, probe_path_eff=None,
        reference_path_eff=None) -> dict:
    """Shot-noise and detector-noise equivalent magnetic sensitivity.

    Returns uT/sqrt(Hz), pT/sqrt(Hz), and coherent-reference values for the same
    optical powers and electronic weight. `S` is the normalized balanced noise
    from `balanced_twin_beam_noise`; `dT_per_ut` is the Hanle transmission slope.
    """
    R = detector_responsivity_aw(args)
    probe_path_eff = (args.probe_path_eff if probe_path_eff is None
                      else float(probe_path_eff))
    reference_path_eff = (args.reference_path_eff if reference_path_eff is None
                          else float(reference_path_eff))
    detector_dark_current_A = 1e-9 * float(
        getattr(args, "detector_dark_current_na", 0.0))
    detector_current_noise = 1e-12 * float(
        getattr(args, "detector_current_noise_pa_sqrt_hz", 0.0))
    balanced_current_noise = 1e-12 * float(
        getattr(args, "balanced_electronic_noise_pa_sqrt_hz", 0.0))
    measurement_bw = max(float(getattr(args, "measurement_bandwidth_hz", 1.0)),
                         0.0)
    P_s_diode = P_probe0 * probe_path_eff * T
    P_c_diode = P_ref0 * reference_path_eff
    I_s = R * P_s_diode
    I_c = R * P_c_diode
    optical_shot_i = np.sqrt(2.0 * constants.ELEMENTARY_CHARGE
                             * np.maximum(I_s + weight * weight * I_c,
                                          1e-300))
    dark_i = np.sqrt(2.0 * constants.ELEMENTARY_CHARGE
                     * detector_dark_current_A
                     * np.maximum(1.0 + weight * weight, 0.0))
    detector_i = np.sqrt(detector_current_noise * detector_current_noise
                         * np.maximum(1.0 + weight * weight, 0.0)
                         + balanced_current_noise * balanced_current_noise)
    squeezed_optical_i = np.sqrt(np.maximum(S, 1e-300)) * optical_shot_i
    squeezed_total_i = np.sqrt(squeezed_optical_i * squeezed_optical_i
                              + dark_i * dark_i
                              + detector_i * detector_i)
    coherent_total_i = np.sqrt(optical_shot_i * optical_shot_i
                               + dark_i * dark_i
                               + detector_i * detector_i)
    signal_slope_A_per_ut = np.abs(R * P_probe0 * probe_path_eff * dT_per_ut)
    shot_limited_uT = squeezed_optical_i \
        / np.maximum(signal_slope_A_per_ut, 1e-300)
    coherent_shot_uT = optical_shot_i \
        / np.maximum(signal_slope_A_per_ut, 1e-300)
    sens_uT = squeezed_total_i / np.maximum(signal_slope_A_per_ut, 1e-300)
    coh_uT = coherent_total_i / np.maximum(signal_slope_A_per_ut, 1e-300)
    return dict(
        responsivity_aw=R,
        photocurrent_probe_A=I_s,
        photocurrent_reference_A=I_c,
        optical_shot_current_A_per_sqrtHz=optical_shot_i,
        squeezed_optical_current_A_per_sqrtHz=squeezed_optical_i,
        detector_dark_current_A=detector_dark_current_A,
        detector_current_noise_A_per_sqrtHz=detector_i,
        total_squeezed_current_A_per_sqrtHz=squeezed_total_i,
        total_coherent_current_A_per_sqrtHz=coherent_total_i,
        shot_limited_uT_per_sqrtHz=shot_limited_uT,
        shot_limited_pT_per_sqrtHz=shot_limited_uT * 1e6,
        coherent_shot_uT_per_sqrtHz=coherent_shot_uT,
        coherent_shot_pT_per_sqrtHz=coherent_shot_uT * 1e6,
        sensitivity_uT_per_sqrtHz=sens_uT,
        sensitivity_pT_per_sqrtHz=sens_uT * 1e6,
        coherent_uT_per_sqrtHz=coh_uT,
        coherent_pT_per_sqrtHz=coh_uT * 1e6,
        sensitivity_rms_pT= sens_uT * 1e6 * math.sqrt(measurement_bw),
        coherent_rms_pT=coh_uT * 1e6 * math.sqrt(measurement_bw),
    )


def point_noise_metrics(source: dict, args: argparse.Namespace, T: float,
                        dT_per_ut: float, *, probe_path_eff: float,
                        reference_path_eff: float,
                        weight_scale: float = 1.0) -> dict:
    cfg = source["config"]
    src = source_readout_values(source, args)
    Gs, Gc = src["G_s"], src["G_c"]
    probe_pre = effective_pre_atten(source, "probe", args)
    reference_pre = effective_pre_atten(source, "reference", args)
    P_probe0 = max(Gs * cfg.seed_w * probe_pre, 0.0)
    P_ref0 = max(Gc * cfg.seed_w * reference_pre, 0.0)
    eta_s = np.clip(
        probe_pre * probe_path_eff * args.detector_qe * T,
        0.0, 1.0)
    eta_c = np.clip(
        reference_pre * reference_path_eff * args.detector_qe,
        0.0, 1.0)
    mean_s = eta_s * Gs
    mean_c = eta_c * Gc
    w = np.sqrt(mean_s / max(mean_c, 1e-30)) * float(weight_scale)
    S = observables.balanced_twin_beam_noise(
        Gs, Gc, eta_s=eta_s, eta_c=eta_c, reference_weight=w,
        source_noise=src["source_noise_linear"],
        seed_excess_noise=args.seed_excess_noise,
        reference_excess_noise=args.reference_excess_noise)
    abs_sens = absolute_field_sensitivity(
        S, w, P_probe0, P_ref0, T, dT_per_ut, args,
        probe_path_eff=probe_path_eff,
        reference_path_eff=reference_path_eff)
    return dict(
        S=float(S),
        S_dB=float(10.0 * np.log10(max(float(S), 1e-300))),
        weight=float(w),
        sensitivity_pT_per_sqrtHz=float(abs_sens["sensitivity_pT_per_sqrtHz"]),
        coherent_pT_per_sqrtHz=float(abs_sens["coherent_pT_per_sqrtHz"]),
        shot_limited_pT_per_sqrtHz=float(abs_sens["shot_limited_pT_per_sqrtHz"]),
        coherent_shot_pT_per_sqrtHz=float(abs_sens["coherent_shot_pT_per_sqrtHz"]),
    )


def tolerance_sweeps(source: dict, hanle: dict, best_idx: int,
                     args: argparse.Namespace) -> dict:
    T0 = float(hanle["transmission"][best_idx])
    dT0 = float(hanle["slope_per_ut"][best_idx])
    pe = np.linspace(args.loss_sweep_min, 1.0, args.loss_sweep_points)
    re = np.linspace(args.loss_sweep_min, 1.0, args.loss_sweep_points)
    loss_noise = np.full((pe.size, re.size), np.nan)
    loss_sens = np.full_like(loss_noise, np.nan)
    for i, ep in enumerate(pe):
        for j, er in enumerate(re):
            m = point_noise_metrics(source, args, T0, dT0,
                                    probe_path_eff=float(ep),
                                    reference_path_eff=float(er))
            loss_noise[i, j] = m["S_dB"]
            loss_sens[i, j] = m["sensitivity_pT_per_sqrtHz"]

    weight_error = np.linspace(-args.balance_error_max,
                               args.balance_error_max,
                               args.balance_error_points)
    balance_noise = np.full(weight_error.size, np.nan)
    balance_sens = np.full_like(balance_noise, np.nan)
    for k, err in enumerate(weight_error):
        m = point_noise_metrics(source, args, T0, dT0,
                                probe_path_eff=args.probe_path_eff,
                                reference_path_eff=args.reference_path_eff,
                                weight_scale=1.0 + float(err))
        balance_noise[k] = m["S_dB"]
        balance_sens[k] = m["sensitivity_pT_per_sqrtHz"]

    best_loss_idx = np.unravel_index(int(np.nanargmin(loss_sens)),
                                     loss_sens.shape)
    keeps_3db = loss_noise <= -3.0
    probe_margin = float(np.nanmin(pe[np.any(keeps_3db, axis=1)])) \
        if np.any(keeps_3db) else float("nan")
    ref_margin = float(np.nanmin(re[np.any(keeps_3db, axis=0)])) \
        if np.any(keeps_3db) else float("nan")
    ok_balance = np.where(balance_noise <= -3.0)[0]
    if ok_balance.size:
        balance_window = (float(weight_error[ok_balance[0]]),
                          float(weight_error[ok_balance[-1]]))
    else:
        balance_window = (float("nan"), float("nan"))
    return dict(
        probe_eff_axis=pe,
        reference_eff_axis=re,
        loss_noise_dB=loss_noise,
        loss_sensitivity_pT=loss_sens,
        best_loss_probe_eff=float(pe[best_loss_idx[0]]),
        best_loss_reference_eff=float(re[best_loss_idx[1]]),
        best_loss_sensitivity_pT=float(loss_sens[best_loss_idx]),
        min_probe_eff_for_3dB=probe_margin,
        min_reference_eff_for_3dB=ref_margin,
        balance_error_axis=weight_error,
        balance_noise_dB=balance_noise,
        balance_sensitivity_pT=balance_sens,
        balance_error_window_3dB=balance_window,
    )


def weighted_noise_over_hanle(source: dict, hanle: dict,
                              args: argparse.Namespace,
                              balance_mode: str) -> dict:
    cfg = source["config"]
    src = source_readout_values(source, args)
    Gs = src["G_s"]
    Gc = src["G_c"]
    probe_pre = effective_pre_atten(source, "probe", args)
    reference_pre = effective_pre_atten(source, "reference", args)
    P_probe0 = max(Gs * cfg.seed_w * probe_pre, 0.0)
    P_ref0 = max(Gc * cfg.seed_w * reference_pre, 0.0)
    T = np.clip(hanle["transmission"], 0.0, 1.5)
    dT = hanle["slope_per_ut"]
    eta_s = np.clip(
        probe_pre * args.probe_path_eff * args.detector_qe * T,
        0.0, 1.0)
    eta_c = np.clip(
        reference_pre * args.reference_path_eff * args.detector_qe,
        0.0, 1.0)
    P_s = cfg.seed_w * Gs * eta_s
    P_c = cfg.seed_w * Gc * eta_c
    S = observables.balanced_twin_beam_noise(
        Gs, Gc, eta_s=eta_s, eta_c=eta_c, reference_weight=balance_mode,
        source_noise=src["source_noise_linear"],
        seed_excess_noise=args.seed_excess_noise,
        reference_excess_noise=args.reference_excess_noise)
    # Reconstruct the electronic weight used by the observable for sensitivity
    # estimates in optical-power units.
    if balance_mode == "raw":
        w = np.ones_like(P_s)
    elif balance_mode == "dc":
        w = P_s / np.maximum(P_c, 1e-30)
    elif balance_mode == "shot":
        w = np.sqrt(P_s / np.maximum(P_c, 1e-30))
    else:
        raise ValueError(balance_mode)
    shot_power = P_s + w * w * P_c
    signal_slope = np.abs(P_probe0 * args.probe_path_eff
                          * args.detector_qe * dT)
    rel_field_noise = np.sqrt(np.maximum(S * shot_power, 1e-300)) \
        / np.maximum(signal_slope, 1e-300)
    coh_field_noise = np.sqrt(np.maximum(shot_power, 1e-300)) \
        / np.maximum(signal_slope, 1e-300)
    slope_valid = np.abs(dT) >= args.min_hanle_slope_per_ut
    abs_sens = absolute_field_sensitivity(
        S, w, P_probe0, P_ref0, T, dT, args)
    total_field_noise = abs_sens["sensitivity_pT_per_sqrtHz"]
    valid = np.isfinite(total_field_noise) & (signal_slope > 0) & slope_valid
    usable = bool(np.any(valid))
    if usable:
        idx = int(np.nanargmin(np.where(valid, total_field_noise, np.nan)))
    else:
        idx = int(np.nanargmax(np.abs(dT)))
    return dict(
        balance_mode=balance_mode,
        S=S,
        S_dB=10.0 * np.log10(np.maximum(S, 1e-300)),
        weight=w,
        P_s=P_s,
        P_c=P_c,
        rel_field_noise=rel_field_noise,
        coh_field_noise=coh_field_noise,
        absolute=abs_sens,
        best_idx=idx,
        usable=usable,
        source_readout=src,
        best_noise_reduction_dB=-10.0 * math.log10(max(float(S[idx]), 1e-300)),
        best_sensitivity_gain_dB=-5.0 * math.log10(max(float(S[idx]), 1e-300)),
    )


def hanle_operating_scan(source: dict, args: argparse.Namespace) -> dict:
    temp_axis = parse_float_axis(args.hanle_temp_axis_c)
    intensity_axis = parse_float_axis(args.hanle_intensity_axis)
    cell_axis = parse_float_axis(args.hanle_cell_axis_mm)
    shape = (temp_axis.size, intensity_axis.size, cell_axis.size)
    sensitivity = np.full(shape, np.nan)
    shot_limited = np.full(shape, np.nan)
    coherent = np.full(shape, np.nan)
    noise = np.full(shape, np.nan)
    best_b = np.full(shape, np.nan)
    transmission = np.full(shape, np.nan)
    slope = np.full(shape, np.nan)
    weight = np.full(shape, np.nan)
    equiv_waist_um = np.full(shape, np.nan)
    usable = np.zeros(shape, dtype=bool)
    best = None
    probe_power_w = arm_input_power_w(source, "probe", args)

    for ti, temp_c in enumerate(temp_axis):
        for ii, intensity in enumerate(intensity_axis):
            for li, cell_mm in enumerate(cell_axis):
                scan_args = args_with(
                    args,
                    hanle_temp_c=float(temp_c),
                    hanle_intensity=float(intensity),
                    hanle_intensity_mode="explicit",
                    hanle_cell_mm=float(cell_mm),
                )
                hanle = compute_hanle(
                    source["probe_detuning_from_87_mhz"], scan_args,
                    probe_power_w=probe_power_w)
                q = weighted_noise_over_hanle(source, hanle, scan_args, "shot")
                if not q["usable"]:
                    continue
                idx = q["best_idx"]
                value = float(
                    q["absolute"]["sensitivity_pT_per_sqrtHz"][idx])
                if not np.isfinite(value):
                    continue
                key = (ti, ii, li)
                usable[key] = True
                sensitivity[key] = value
                shot_limited[key] = float(
                    q["absolute"]["shot_limited_pT_per_sqrtHz"][idx])
                coherent[key] = float(
                    q["absolute"]["coherent_pT_per_sqrtHz"][idx])
                noise[key] = float(q["S_dB"][idx])
                best_b[key] = float(hanle["b_ut"][idx])
                transmission[key] = float(hanle["transmission"][idx])
                slope[key] = float(hanle["slope_per_ut"][idx])
                weight[key] = float(q["weight"][idx])
                equiv_waist_um[key] = gaussian_waist_um_for_peak_intensity(
                    probe_power_w, float(hanle["params"]["intensity_mw_cm2"]))
                if best is None or value < best["sensitivity_pT"]:
                    best = dict(
                        temp_idx=ti,
                        intensity_idx=ii,
                        cell_idx=li,
                        temp_c=float(temp_c),
                        intensity_mw_cm2=float(intensity),
                        cell_mm=float(cell_mm),
                        sensitivity_pT=value,
                        shot_limited_pT=shot_limited[key],
                        coherent_pT=coherent[key],
                        noise_dB=noise[key],
                        bias_ut=best_b[key],
                        transmission=transmission[key],
                        slope_per_ut=slope[key],
                        weight=weight[key],
                        equivalent_waist_um=equiv_waist_um[key],
                    )

    return dict(
        temp_axis_c=temp_axis,
        intensity_axis_mw_cm2=intensity_axis,
        cell_axis_mm=cell_axis,
        sensitivity_pT=sensitivity,
        shot_limited_pT=shot_limited,
        coherent_pT=coherent,
        noise_dB=noise,
        best_B_ut=best_b,
        transmission=transmission,
        slope_per_ut=slope,
        weight=weight,
        equivalent_waist_um=equiv_waist_um,
        usable=usable,
        best=best,
    )


def summarize_frequency_table() -> str:
    rows = [
        ("85Rb F=2 -> F'=3", 0.0),
        ("85Rb F=2 -> F'=2", offset_from_fwm_origin_ghz(species.RB85, 2, 2)),
        ("85Rb F=3 -> F'=3", offset_from_fwm_origin_ghz(species.RB85, 3, 3)),
        ("85Rb F=3 -> F'=2", offset_from_fwm_origin_ghz(species.RB85, 3, 2)),
        ("87Rb F=2 -> F'=2", offset_from_fwm_origin_ghz(species.RB87, 2, 2)),
        ("87Rb F=1 -> F'=1", offset_from_fwm_origin_ghz(species.RB87, 1, 1)),
    ]
    out = ["| Transition | offset from 85Rb F=2 -> F'=3 [GHz] |",
           "|---|---:|"]
    for label, val in rows:
        out.append(f"| {label} | {val:+.6f} |")
    return "\n".join(out)


def write_report(results: list[dict], args: argparse.Namespace,
                 out_md: Path, calibration: dict | None = None) -> None:
    lines = [
        "# Resonant 85Rb-FWM / 87Rb-Hanle squeezing reference",
        "",
        "This report is generated by `analysis/squeezing/resonant_hanle_squeezing_reference.py`.",
        "The FWM probe frequency is referenced to 85Rb D1 F=2 -> F'=3.",
        "The Hanle cell is 87Rb D1 F=2 -> F'=2 transmission.",
        "",
        "## Frequency anchors",
        "",
        summarize_frequency_table(),
        "",
        "## Source and magnetometer summary",
        "",
        "| Source | source used | D [GHz] | delta_code [MHz] | paper delta equiv [MHz] | probe detuning from 87Rb F=2->F'=2 [MHz] | conjugate detuning from 87Rb F=1->F'=1 [MHz] | model Gs | used Gs | used Gc | source S used [dB] | probe into Hanle [uW] | Hanle I [mW/cm^2] | equiv waist [um] | Hanle bias [uT] | T(bias) | abs(dT/dB) [1/uT] | balanced noise [dB] | noise reduction [dB] | squeezed total [pT/sqrtHz] | squeezed shot-only [pT/sqrtHz] | coherent total [pT/sqrtHz] |",
        "|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for r in results:
        cfg = r["source"]["config"]
        op = r["source"]["op"]
        h = r["hanle"]
        q = r["balanced"]["shot"]
        src = q["source_readout"]
        i = q["best_idx"]
        best_noise = f"{q['S_dB'][i]:.2f}" if q["usable"] else "n/a"
        noise_reduction = (f"{q['best_noise_reduction_dB']:.2f}"
                           if q["usable"] else "n/a")
        sens = (f"{q['absolute']['sensitivity_pT_per_sqrtHz'][i]:.3g}"
                if q["usable"] else "n/a")
        shot_sens = (f"{q['absolute']['shot_limited_pT_per_sqrtHz'][i]:.3g}"
                     if q["usable"] else "n/a")
        coh = (f"{q['absolute']['coherent_pT_per_sqrtHz'][i]:.3g}"
               if q["usable"] else "n/a")
        probe_power_uw = arm_input_power_w(r["source"], "probe", args) * 1e6
        hanle_intensity = float(h["params"]["intensity_mw_cm2"])
        equiv_waist = gaussian_waist_um_for_peak_intensity(
            probe_power_uw * 1e-6, hanle_intensity)
        lines.append(
            f"| {cfg.label} | {src['label']} | {cfg.D_GHz:+.6f} "
            f"| {cfg.delta_code_mhz:+.2f} "
            f"| {r['source']['paper_delta_equiv_mhz']:+.2f} "
            f"| {r['source']['probe_detuning_from_87_mhz']:+.2f} "
            f"| {r['source']['conjugate_detuning_from_87_mhz']:+.2f} "
            f"| {op['G_s']:.3g} | {src['G_s']:.3g} | {src['G_c']:.3g} "
            f"| {src['S_dB']:.2f} "
            f"| {probe_power_uw:.3g} | {hanle_intensity:.3g} "
            f"| {equiv_waist:.3g} "
            f"| {h['b_ut'][i]:+.4f} | {h['transmission'][i]:.4f} "
            f"| {abs(h['slope_per_ut'][i]):.3e} "
            f"| {best_noise} | {noise_reduction} | {sens} "
            f"| {shot_sens} | {coh} |"
        )
    lines += [
        "",
        "## Tolerance summary",
        "",
        "| Source | min probe path eff for >3 dB | min reference path eff for >3 dB | balance error window for >3 dB | best sensitivity in loss grid [pT/sqrtHz] | best loss-grid efficiencies probe/ref |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for r in results:
        cfg = r["source"]["config"]
        tol = r.get("tolerance")
        q = r["balanced"]["shot"]
        if tol is None or not q["usable"]:
            lines.append(f"| {cfg.label} | n/a | n/a | n/a | n/a | n/a |")
            continue
        lo, hi = tol["balance_error_window_3dB"]
        win = "n/a" if not np.isfinite(lo) else f"{100*lo:+.1f}% to {100*hi:+.1f}%"
        pmin = tol["min_probe_eff_for_3dB"]
        rmin = tol["min_reference_eff_for_3dB"]
        lines.append(
            f"| {cfg.label} | "
            f"{pmin:.2f}" if np.isfinite(pmin) else f"| {cfg.label} | n/a"
        )
        lines[-1] += (
            f" | {rmin:.2f}" if np.isfinite(rmin) else " | n/a"
        )
        lines[-1] += (
            f" | {win} | {tol['best_loss_sensitivity_pT']:.3g} "
            f"| {tol['best_loss_probe_eff']:.2f}/{tol['best_loss_reference_eff']:.2f} |"
        )
    lines += [
        "",
        "## Hanle operating-point scan",
        "",
        "| Source | best 87Rb temp [C] | intensity [mW/cm^2] | equiv waist [um] | cell length [mm] | Hanle bias [uT] | T(bias) | abs(dT/dB) [1/uT] | balanced noise [dB] | squeezed total [pT/sqrtHz] | squeezed shot-only [pT/sqrtHz] | coherent total [pT/sqrtHz] |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for r in results:
        cfg = r["source"]["config"]
        opt = r.get("optimization")
        best = opt.get("best") if opt is not None else None
        if best is None:
            lines.append(
                f"| {cfg.label} | n/a | n/a | n/a | n/a | n/a | n/a | n/a | n/a | n/a | n/a | n/a |")
            continue
        lines.append(
            f"| {cfg.label} | {best['temp_c']:.1f} "
            f"| {best['intensity_mw_cm2']:.3g} "
            f"| {best['equivalent_waist_um']:.3g} "
            f"| {best['cell_mm']:.3g} "
            f"| {best['bias_ut']:+.4f} "
            f"| {best['transmission']:.4f} "
            f"| {abs(best['slope_per_ut']):.3e} "
            f"| {best['noise_dB']:.2f} "
            f"| {best['sensitivity_pT']:.3g} "
            f"| {best['shot_limited_pT']:.3g} "
            f"| {best['coherent_pT']:.3g} |"
        )
    lines += [
        "",
        "## Measured Hanle calibration",
        "",
    ]
    if calibration is None:
        lines.append("No measured Hanle CSV was provided.")
    else:
        data = calibration["data"]
        fit = calibration["fit"]
        lines += [
            f"Calibration source: `{calibration['source_name']}` ({calibration['source_label']}).",
            "",
            "| data file | B column | signal column | points | B scale | B offset [uT] | signal offset | signal scale | RMS | normalized RMS | R^2 |",
            "|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|",
            f"| {Path(data['path']).name} | {data['b_column']} | {data['signal_column']} "
            f"| {fit['n_fit_points']} | {fit['b_scale']:.6g} "
            f"| {fit['b_offset_ut']:+.6g} | {fit['signal_offset']:.6g} "
            f"| {fit['signal_scale']:.6g} | {fit['rms']:.6g} "
            f"| {fit['normalized_rms']:.6g} | {fit['r2']:.6g} |",
        ]
    lines += [
        "",
        "Notes:",
        "- `delta_code` is the GABES branch=-1 two-photon-detuning convention.",
        "- `paper delta equiv` is shown with the opposite sign for comparison with the OE beat-note convention.",
        "- For OE resonant configurations, `paper_anchor` uses the reported/inferred OE gain and source IDS because the compact GABES 4-level model is Sim-anchored and under-predicts this resonant-FWM gain without a separate OE calibration.",
        f"- Rows with `n/a` do not meet the minimum Hanle slope threshold, abs(dT/dB) >= {args.min_hanle_slope_per_ut:.1e} 1/uT.",
        "- `best balanced noise` uses electronic shot-weight balancing, w=sqrt(Is/Ic), at the best field-sensitivity point.",
        "- Absolute sensitivity uses shot current noise plus optional dark-current and electronic-current noise floors. If `--responsivity-aw` is left at 0, it is computed from detector QE and wavelength.",
        "- Tolerance rows use the best Hanle bias point and scan path efficiencies plus electronic balance gain error around the shot-weight optimum.",
        "- `equiv waist` is the Gaussian 1/e^2 field-radius waist that would give the listed peak intensity at the simulated probe power. Use `--hanle-intensity-mode waist` with a measured waist to make this coupling explicit.",
        "- Measured Hanle calibration fits `signal = offset + scale * model(B_scale * B_measured + B_offset)`. This separates coil/current-axis calibration from arbitrary detector voltage normalization.",
        "- The Hanle operating-point scan varies the susceptibility-model intensity, 87Rb cell temperature, and cell length. Optical powers are set by the FWM source powers, measured-power overrides, and configured pre-attenuations.",
        "- The Hanle model is compact and semi-quantitative; buffer-gas diffusion, pressure shift and Dicke narrowing are still phenomenological unless calibrated.",
        "",
        "## Command parameters",
        "",
        "```text",
        str(vars(args)),
        "```",
    ]
    out_md.write_text("\n".join(lines), encoding="utf-8")


def plot_results(results: list[dict], out_png: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 2, figsize=(13.5, 9.0))
    ax0, ax1, ax2, ax3 = axes.ravel()
    rb87_probe = offset_from_fwm_origin_ghz(*RB87_PROBE)

    for r in results:
        cfg = r["source"]["config"]
        spec = r["source"]["spectrum"]
        h = r["hanle"]
        q = r["balanced"]["shot"]
        x = (spec["probe_axis_GHz"] - rb87_probe) * 1e3
        ax0.plot(x, spec["G_s"], lw=1.3, label=cfg.label)
        ax1.plot(x, spec["S_dB"], lw=1.3, label=cfg.label)
        ax2.plot(h["b_ut"], h["transmission"], lw=1.3, label=cfg.label)
        ax3.plot(h["b_ut"], q["S_dB"], lw=1.3, label=cfg.label)
        i = q["best_idx"]
        ax2.scatter([h["b_ut"][i]], [h["transmission"][i]], s=22)
        ax3.scatter([h["b_ut"][i]], [q["S_dB"][i]], s=22)

    ax0.axvline(0, color="black", lw=0.8, ls=":")
    ax1.axvline(0, color="black", lw=0.8, ls=":")
    ax0.set_xlabel("Probe detuning from 87Rb F=2 -> F'=2 [MHz]")
    ax0.set_ylabel("Probe gain Gs")
    ax0.set_title("85Rb FWM gain near 87Rb resonance")
    ax0.set_yscale("log")
    ax1.set_xlabel("Probe detuning from 87Rb F=2 -> F'=2 [MHz]")
    ax1.set_ylabel("IDS [dB]")
    ax1.set_title("Source intensity-difference squeezing")
    ax2.set_xlabel("B [uT]")
    ax2.set_ylabel("87Rb Hanle transmission")
    ax2.set_title("Hanle transfer function seen by probe")
    ax3.set_xlabel("B [uT]")
    ax3.set_ylabel("Balanced noise [dB]")
    ax3.set_title("After Hanle probe loss + external conjugate reference")
    for ax in axes.ravel():
        ax.grid(alpha=0.25)
    ax0.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(out_png, dpi=140)


def plot_tolerance(results: list[dict], out_png: Path,
                   args: argparse.Namespace) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    usable = [r for r in results if r["balanced"]["shot"]["usable"]
              and r.get("tolerance") is not None]
    if not usable:
        return
    best = min(
        usable,
        key=lambda r: r["balanced"]["shot"]["absolute"]["sensitivity_pT_per_sqrtHz"]
        [r["balanced"]["shot"]["best_idx"]])
    cfg = best["source"]["config"]
    tol = best["tolerance"]
    q = best["balanced"]["shot"]
    h = best["hanle"]
    i = q["best_idx"]

    fig, (ax0, ax1) = plt.subplots(1, 2, figsize=(13.0, 4.8))
    im = ax0.imshow(
        tol["loss_noise_dB"].T, origin="lower", aspect="auto",
        extent=[tol["probe_eff_axis"][0], tol["probe_eff_axis"][-1],
                tol["reference_eff_axis"][0], tol["reference_eff_axis"][-1]],
        cmap="viridis_r")
    fig.colorbar(im, ax=ax0, label="balanced noise [dB]")
    ax0.contour(tol["probe_eff_axis"], tol["reference_eff_axis"],
                tol["loss_noise_dB"].T, levels=[-6.0, -3.0, -1.0],
                colors="white", linewidths=0.9)
    ax0.scatter([args.probe_path_eff], [args.reference_path_eff],
                color="crimson", s=36, label="default")
    ax0.set_xlabel("probe path efficiency")
    ax0.set_ylabel("reference path efficiency")
    ax0.set_title(f"{cfg.label}: path-loss tolerance")
    ax0.legend(fontsize=8)

    ax1.plot(100.0 * tol["balance_error_axis"], tol["balance_noise_dB"],
             color="#1f77b4", lw=1.8)
    ax1.axhline(-3.0, color="crimson", ls="--", lw=1.0)
    ax1.axvline(0, color="black", ls=":", lw=0.8)
    ax1.set_xlabel("electronic balance error from shot-weight optimum [%]")
    ax1.set_ylabel("balanced noise [dB]")
    ax1.set_title(f"Balance tolerance at B={h['b_ut'][i]:+.3f} uT")
    ax1.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(out_png, dpi=140)


def plot_hanle_optimization(results: list[dict], out_png: Path) -> bool:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    usable = [r for r in results
              if r.get("optimization") is not None
              and r["optimization"].get("best") is not None]
    if not usable:
        return False
    best_result = min(
        usable,
        key=lambda r: r["optimization"]["best"]["sensitivity_pT"])
    cfg = best_result["source"]["config"]
    opt = best_result["optimization"]
    best = opt["best"]
    ti = best["temp_idx"]
    ii = best["intensity_idx"]
    li = best["cell_idx"]
    sens = opt["sensitivity_pT"][ti]
    noise = opt["noise_dB"][ti]
    intensity_axis = opt["intensity_axis_mw_cm2"]
    cell_axis = opt["cell_axis_mm"]

    fig, (ax0, ax1) = plt.subplots(1, 2, figsize=(13.0, 4.8))
    im0 = ax0.imshow(sens.T, origin="lower", aspect="auto",
                     cmap="magma_r")
    fig.colorbar(im0, ax=ax0, label="squeezed sensitivity [pT/sqrtHz]")
    ax0.scatter([ii], [li], color="cyan", edgecolors="black", s=44,
                label="best")
    ax0.set_title(f"{cfg.label}: best temp {best['temp_c']:.1f} C")
    ax0.set_xlabel("Hanle probe intensity [mW/cm^2]")
    ax0.set_ylabel("87Rb cell length [mm]")
    ax0.set_xticks(np.arange(intensity_axis.size))
    ax0.set_xticklabels([f"{x:g}" for x in intensity_axis], rotation=35)
    ax0.set_yticks(np.arange(cell_axis.size))
    ax0.set_yticklabels([f"{x:g}" for x in cell_axis])
    ax0.legend(fontsize=8)

    im1 = ax1.imshow(noise.T, origin="lower", aspect="auto",
                     cmap="viridis_r")
    fig.colorbar(im1, ax=ax1, label="balanced noise [dB]")
    ax1.scatter([ii], [li], color="cyan", edgecolors="black", s=44)
    ax1.set_title("Balanced noise at optimum B bias")
    ax1.set_xlabel("Hanle probe intensity [mW/cm^2]")
    ax1.set_ylabel("87Rb cell length [mm]")
    ax1.set_xticks(np.arange(intensity_axis.size))
    ax1.set_xticklabels([f"{x:g}" for x in intensity_axis], rotation=35)
    ax1.set_yticks(np.arange(cell_axis.size))
    ax1.set_yticklabels([f"{x:g}" for x in cell_axis])
    fig.tight_layout()
    fig.savefig(out_png, dpi=140)
    return True


def plot_hanle_calibration(calibration: dict | None,
                           out_png: Path) -> bool:
    if calibration is None:
        return False
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    data = calibration["data"]
    fit = calibration["fit"]
    fig, (ax0, ax1) = plt.subplots(2, 1, figsize=(8.8, 7.0),
                                   sharex=True,
                                   gridspec_kw={"height_ratios": [3, 1]})
    ax0.plot(data["b_ut"], data["signal"], "o", ms=3.5,
             label="measured")
    ax0.plot(data["b_ut"], fit["fitted_signal"], "-", lw=1.8,
             label="fit")
    ax0.set_ylabel("signal [arb.]")
    ax0.set_title(
        f"Measured Hanle calibration: {calibration['source_label']}")
    ax0.grid(alpha=0.25)
    ax0.legend(fontsize=8)

    ax1.axhline(0.0, color="black", lw=0.8, ls=":")
    ax1.plot(data["b_ut"], fit["residual"], "o-", ms=3.0, lw=1.0)
    ax1.set_xlabel("measured B axis [uT]")
    ax1.set_ylabel("residual")
    ax1.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(out_png, dpi=140)
    return True


def save_npz(results: list[dict], out_npz: Path,
             args: argparse.Namespace,
             calibration: dict | None = None) -> None:
    names = np.array([r["source"]["config"].name for r in results])
    used = [r["balanced"]["shot"]["source_readout"] for r in results]
    D = np.array([r["source"]["config"].D_GHz for r in results])
    delta = np.array([r["source"]["config"].delta_code_mhz for r in results])
    probe_det = np.array([r["source"]["probe_detuning_from_87_mhz"]
                          for r in results])
    conj_det = np.array([r["source"]["conjugate_detuning_from_87_mhz"]
                         for r in results])
    model_Gs = np.array([r["source"]["op"]["G_s"] for r in results])
    model_Gc = np.array([r["source"]["op"]["G_c"] for r in results])
    model_source_S = np.array([r["source"]["op"]["S_dB"] for r in results])
    Gs = np.array([u["G_s"] for u in used])
    Gc = np.array([u["G_c"] for u in used])
    source_S = np.array([u["S_dB"] for u in used])
    source_labels = np.array([u["label"] for u in used])
    best_noise = np.array([r["balanced"]["shot"]["S_dB"]
                           [r["balanced"]["shot"]["best_idx"]]
                           for r in results])
    best_B = np.array([r["hanle"]["b_ut"][r["balanced"]["shot"]["best_idx"]]
                       for r in results])
    usable = np.array([r["balanced"]["shot"]["usable"] for r in results],
                      dtype=bool)
    sens = np.array([
        r["balanced"]["shot"]["absolute"]["sensitivity_pT_per_sqrtHz"]
        [r["balanced"]["shot"]["best_idx"]]
        for r in results
    ])
    shot_sens = np.array([
        r["balanced"]["shot"]["absolute"]["shot_limited_pT_per_sqrtHz"]
        [r["balanced"]["shot"]["best_idx"]]
        for r in results
    ])
    coh_sens = np.array([
        r["balanced"]["shot"]["absolute"]["coherent_pT_per_sqrtHz"]
        [r["balanced"]["shot"]["best_idx"]]
        for r in results
    ])
    coh_shot_sens = np.array([
        r["balanced"]["shot"]["absolute"]["coherent_shot_pT_per_sqrtHz"]
        [r["balanced"]["shot"]["best_idx"]]
        for r in results
    ])
    probe_power_uw = np.array([
        arm_input_power_w(r["source"], "probe", args) * 1e6
        for r in results
    ])
    reference_power_uw = np.array([
        arm_input_power_w(r["source"], "reference", args) * 1e6
        for r in results
    ])
    hanle_intensity = np.array([
        float(r["hanle"]["params"]["intensity_mw_cm2"])
        for r in results
    ])
    equivalent_waist_um = np.array([
        gaussian_waist_um_for_peak_intensity(
            arm_input_power_w(r["source"], "probe", args),
            float(r["hanle"]["params"]["intensity_mw_cm2"]))
        for r in results
    ])
    if any(r.get("tolerance") is not None for r in results):
        first_tol = next(r["tolerance"] for r in results
                         if r.get("tolerance") is not None)
        loss_noise = np.stack([
            r["tolerance"]["loss_noise_dB"]
            if r.get("tolerance") is not None
            else np.full_like(first_tol["loss_noise_dB"], np.nan)
            for r in results
        ])
        balance_noise = np.stack([
            r["tolerance"]["balance_noise_dB"]
            if r.get("tolerance") is not None
            else np.full_like(first_tol["balance_noise_dB"], np.nan)
            for r in results
        ])
        probe_eff_axis = first_tol["probe_eff_axis"]
        reference_eff_axis = first_tol["reference_eff_axis"]
        balance_error_axis = first_tol["balance_error_axis"]
    else:
        loss_noise = np.empty((len(results), 0, 0))
        balance_noise = np.empty((len(results), 0))
        probe_eff_axis = reference_eff_axis = balance_error_axis = np.array([])
    if any(r.get("optimization") is not None for r in results):
        first_opt = next(r["optimization"] for r in results
                         if r.get("optimization") is not None)
        opt_temp_axis = first_opt["temp_axis_c"]
        opt_intensity_axis = first_opt["intensity_axis_mw_cm2"]
        opt_cell_axis = first_opt["cell_axis_mm"]
        opt_shape = first_opt["sensitivity_pT"].shape
        opt_sensitivity = np.stack([
            r["optimization"]["sensitivity_pT"]
            if r.get("optimization") is not None
            else np.full(opt_shape, np.nan)
            for r in results
        ])
        opt_shot = np.stack([
            r["optimization"]["shot_limited_pT"]
            if r.get("optimization") is not None
            else np.full(opt_shape, np.nan)
            for r in results
        ])
        opt_coherent = np.stack([
            r["optimization"]["coherent_pT"]
            if r.get("optimization") is not None
            else np.full(opt_shape, np.nan)
            for r in results
        ])
        opt_noise = np.stack([
            r["optimization"]["noise_dB"]
            if r.get("optimization") is not None
            else np.full(opt_shape, np.nan)
            for r in results
        ])
        opt_best_B = np.stack([
            r["optimization"]["best_B_ut"]
            if r.get("optimization") is not None
            else np.full(opt_shape, np.nan)
            for r in results
        ])
        opt_equiv_waist = np.stack([
            r["optimization"]["equivalent_waist_um"]
            if r.get("optimization") is not None
            else np.full(opt_shape, np.nan)
            for r in results
        ])
        opt_usable = np.stack([
            r["optimization"]["usable"]
            if r.get("optimization") is not None
            else np.zeros(opt_shape, dtype=bool)
            for r in results
        ])

        def opt_best_value(result: dict, key: str) -> float:
            opt = result.get("optimization")
            best = opt.get("best") if opt is not None else None
            return float(best[key]) if best is not None else float("nan")

        opt_best_sens = np.array([
            opt_best_value(r, "sensitivity_pT") for r in results])
        opt_best_shot = np.array([
            opt_best_value(r, "shot_limited_pT") for r in results])
        opt_best_coh = np.array([
            opt_best_value(r, "coherent_pT") for r in results])
        opt_best_noise = np.array([
            opt_best_value(r, "noise_dB") for r in results])
        opt_best_temp = np.array([
            opt_best_value(r, "temp_c") for r in results])
        opt_best_intensity = np.array([
            opt_best_value(r, "intensity_mw_cm2") for r in results])
        opt_best_cell = np.array([
            opt_best_value(r, "cell_mm") for r in results])
        opt_best_bias = np.array([
            opt_best_value(r, "bias_ut") for r in results])
        opt_best_equiv_waist = np.array([
            opt_best_value(r, "equivalent_waist_um") for r in results])
    else:
        opt_temp_axis = opt_intensity_axis = opt_cell_axis = np.array([])
        empty_opt = np.empty((len(results), 0, 0, 0))
        opt_sensitivity = opt_shot = opt_coherent = empty_opt
        opt_noise = opt_best_B = opt_equiv_waist = empty_opt
        opt_usable = np.empty((len(results), 0, 0, 0), dtype=bool)
        opt_best_sens = opt_best_shot = opt_best_coh = opt_best_noise = np.full(
            len(results), np.nan)
        opt_best_temp = opt_best_intensity = opt_best_cell = opt_best_bias = (
            np.full(len(results), np.nan))
        opt_best_equiv_waist = np.full(len(results), np.nan)
    if calibration is None:
        cal_source = np.array([""])
        cal_data_b = np.array([])
        cal_data_signal = np.array([])
        cal_fit_signal = np.array([])
        cal_residual = np.array([])
        cal_model_b_at_data = np.array([])
        cal_params = np.full(7, np.nan)
    else:
        cal_data = calibration["data"]
        cal_fit = calibration["fit"]
        cal_source = np.array([calibration["source_name"]])
        cal_data_b = cal_data["b_ut"]
        cal_data_signal = cal_data["signal"]
        cal_fit_signal = cal_fit["fitted_signal"]
        cal_residual = cal_fit["residual"]
        cal_model_b_at_data = cal_fit["model_b_at_data_ut"]
        cal_params = np.array([
            cal_fit["b_scale"],
            cal_fit["b_offset_ut"],
            cal_fit["signal_offset"],
            cal_fit["signal_scale"],
            cal_fit["rms"],
            cal_fit["normalized_rms"],
            cal_fit["r2"],
        ])
    np.savez(out_npz, names=names, D_GHz=D, delta_code_mhz=delta,
             probe_detuning_mhz=probe_det, conjugate_detuning_mhz=conj_det,
             source_label=source_labels,
             model_G_s=model_Gs, model_G_c=model_Gc,
             model_source_S_dB=model_source_S,
             G_s=Gs, G_c=Gc, source_S_dB=source_S,
             balanced_S_dB=best_noise, best_B_ut=best_B,
             usable_hanle_reference=usable,
             sensitivity_pT_per_sqrtHz=sens,
             shot_limited_pT_per_sqrtHz=shot_sens,
             coherent_pT_per_sqrtHz=coh_sens,
             coherent_shot_pT_per_sqrtHz=coh_shot_sens,
             probe_input_power_uw=probe_power_uw,
             reference_input_power_uw=reference_power_uw,
             hanle_intensity_mw_cm2=hanle_intensity,
             hanle_equivalent_waist_um=equivalent_waist_um,
             tolerance_probe_eff_axis=probe_eff_axis,
             tolerance_reference_eff_axis=reference_eff_axis,
             tolerance_balance_error_axis=balance_error_axis,
             tolerance_loss_noise_dB=loss_noise,
             tolerance_balance_noise_dB=balance_noise,
             hanle_opt_temp_axis_c=opt_temp_axis,
             hanle_opt_intensity_axis_mw_cm2=opt_intensity_axis,
             hanle_opt_cell_axis_mm=opt_cell_axis,
             hanle_opt_sensitivity_pT=opt_sensitivity,
             hanle_opt_shot_limited_pT=opt_shot,
             hanle_opt_coherent_pT=opt_coherent,
             hanle_opt_noise_dB=opt_noise,
             hanle_opt_best_B_ut_grid=opt_best_B,
             hanle_opt_equivalent_waist_um=opt_equiv_waist,
             hanle_opt_usable=opt_usable,
             hanle_opt_best_sensitivity_pT=opt_best_sens,
             hanle_opt_best_shot_limited_pT=opt_best_shot,
             hanle_opt_best_coherent_pT=opt_best_coh,
             hanle_opt_best_noise_dB=opt_best_noise,
             hanle_opt_best_temp_c=opt_best_temp,
             hanle_opt_best_intensity_mw_cm2=opt_best_intensity,
             hanle_opt_best_equivalent_waist_um=opt_best_equiv_waist,
             hanle_opt_best_cell_mm=opt_best_cell,
             hanle_opt_best_bias_ut=opt_best_bias,
             hanle_calibration_source=cal_source,
             hanle_calibration_params=cal_params,
             hanle_calibration_param_labels=np.array([
                 "b_scale", "b_offset_ut", "signal_offset",
                 "signal_scale", "rms", "normalized_rms", "r2"]),
             hanle_calibration_b_data_ut=cal_data_b,
             hanle_calibration_signal_data=cal_data_signal,
             hanle_calibration_signal_fit=cal_fit_signal,
             hanle_calibration_residual=cal_residual,
             hanle_calibration_model_b_at_data_ut=cal_model_b_at_data)


def apply_json_config(args: argparse.Namespace) -> argparse.Namespace:
    path = getattr(args, "config_json", "")
    if not path:
        return args
    cfg_path = Path(path).expanduser()
    if not cfg_path.is_absolute():
        cfg_path = (Path.cwd() / cfg_path).resolve()
    with cfg_path.open("r", encoding="utf-8") as f:
        data = json.load(f)
    if not isinstance(data, dict):
        raise ValueError("--config-json must point to a JSON object")
    valid = set(vars(args))
    normalized = {str(k).replace("-", "_"): v for k, v in data.items()}
    unknown = sorted(set(normalized) - valid)
    if unknown:
        raise ValueError(f"unknown config key(s): {unknown}")
    for key, value in normalized.items():
        setattr(args, key, value)
    args.config_json = str(cfg_path)
    return args


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--config-json", default="",
                   help="Optional JSON file with argument-name keys to override CLI/defaults.")
    p.add_argument("--fidelity", choices=("quick", "balanced", "fine", "ultra"),
                   default="fine")
    p.add_argument("--source-anchor", choices=("paper", "model"),
                   default="paper",
                   help="Use OE paper gain/IDS anchors when available, or raw GABES model output.")
    p.add_argument("--only", nargs="*", default=None,
                   help="Optional source names to run.")
    p.add_argument("--cell", choices=("paraffin", "buffer"), default="paraffin")
    p.add_argument("--hanle-temp-c", type=float, default=25.0)
    p.add_argument("--hanle-cell-mm", type=float, default=10.0)
    p.add_argument("--hanle-intensity", type=float, default=0.8,
                   help="Probe intensity in the 87Rb Hanle cell [mW/cm^2].")
    p.add_argument("--hanle-intensity-mode", choices=("explicit", "waist"),
                   default="explicit",
                   help="Use --hanle-intensity directly, or compute it from probe power and --hanle-probe-waist-um.")
    p.add_argument("--hanle-probe-waist-um", type=float, default=0.0,
                   help="Gaussian 1/e^2 field-radius waist used when --hanle-intensity-mode=waist.")
    p.add_argument("--hanle-line-strength", type=float, default=1.0)
    p.add_argument("--qwp-deg", type=float, default=0.0)
    p.add_argument("--b-max-ut", type=float, default=2.0)
    p.add_argument("--scan-points", type=int, default=201)
    p.add_argument("--velocity-classes", type=int, default=13)
    p.add_argument("--wall-ms", type=float, default=3.0)
    p.add_argument("--transit-khz", type=float, default=80.0)
    p.add_argument("--dark-return-khz", type=float, default=1.0)
    p.add_argument("--residual-b-ut", type=float, default=0.05)
    p.add_argument("--transverse-angle-deg", type=float, default=0.0)
    p.add_argument("--min-hanle-slope-per-ut", type=float, default=1e-6,
                   help="Rows below this |dT/dB| are marked unusable for Hanle sensing.")
    p.add_argument("--ne-torr", type=float, default=20.0)
    p.add_argument("--buffer-ground-khz", type=float, default=20.0)
    p.add_argument("--collisional-depol-khz", type=float, default=2.0)
    p.add_argument("--probe-path-eff", type=float, default=0.96)
    p.add_argument("--reference-path-eff", type=float, default=0.96)
    p.add_argument("--detector-qe", type=float, default=0.95)
    p.add_argument("--responsivity-aw", type=float, default=0.0,
                   help="Photodiode responsivity. 0 computes QE*lambda/1239.84.")
    p.add_argument("--probe-wavelength-nm", type=float, default=794.98)
    p.add_argument("--probe-pre-atten", type=float, default=1.0)
    p.add_argument("--reference-pre-atten", type=float, default=1.0)
    p.add_argument("--probe-input-power-uw", type=float, default=0.0,
                   help="Measured probe power entering the Hanle cell. 0 uses FWM gain*seed*pre-atten.")
    p.add_argument("--reference-input-power-uw", type=float, default=0.0,
                   help="Measured conjugate/reference power before the reference path. 0 uses FWM gain*seed*pre-atten.")
    p.add_argument("--detector-dark-current-na", type=float, default=0.0,
                   help="Dark current per photodiode, added as unsqueezed shot noise.")
    p.add_argument("--detector-current-noise-pa-sqrt-hz", type=float, default=0.0,
                   help="Input-referred current-noise density per photodiode.")
    p.add_argument("--balanced-electronic-noise-pa-sqrt-hz", type=float,
                   default=0.0,
                   help="Additional differential electronics current-noise density.")
    p.add_argument("--measurement-bandwidth-hz", type=float, default=1.0,
                   help="Bandwidth used for RMS pT columns stored in the NPZ.")
    p.add_argument("--seed-excess-noise", type=float, default=0.0)
    p.add_argument("--reference-excess-noise", type=float, default=0.0)
    p.add_argument("--loss-sweep-min", type=float, default=0.50)
    p.add_argument("--loss-sweep-points", type=int, default=51)
    p.add_argument("--balance-error-max", type=float, default=0.20)
    p.add_argument("--balance-error-points", type=int, default=81)
    p.add_argument("--skip-hanle-optimize", action="store_true",
                   help="Skip the Hanle intensity/temperature/cell-length operating-point scan.")
    p.add_argument("--hanle-temp-axis-c", default="25",
                   help="Comma-separated 87Rb cell temperatures for the operating scan.")
    p.add_argument("--hanle-intensity-axis", default="0.2,0.5,0.8,1.5,3.0",
                   help="Comma-separated Hanle probe intensities [mW/cm^2] for the operating scan.")
    p.add_argument("--hanle-cell-axis-mm", default="5,10,25,50",
                   help="Comma-separated 87Rb cell lengths [mm] for the operating scan.")
    p.add_argument("--hanle-calibration-csv", default="",
                   help="Optional measured Hanle scan CSV for model-vs-data calibration.")
    p.add_argument("--hanle-calibration-source", default="oe_probe_locked",
                   help="Source configuration whose Hanle trace should be fit to the CSV.")
    p.add_argument("--hanle-calibration-b-column", default="",
                   help="Measured B column name, or zero-based index for headerless CSV.")
    p.add_argument("--hanle-calibration-signal-column", default="",
                   help="Measured signal column name, or zero-based index for headerless CSV.")
    p.add_argument("--hanle-calibration-b-unit", default="uT",
                   help="Measured B-axis unit: uT, nT, mT, T, G, or mG.")
    p.add_argument("--hanle-calibration-offset-span-ut", type=float,
                   default=0.5,
                   help="Fit search half-span for B offset [uT].")
    p.add_argument("--hanle-calibration-offset-points", type=int,
                   default=101)
    p.add_argument("--hanle-calibration-scale-span", type=float,
                   default=0.10,
                   help="Fit search half-span around B scale=1.")
    p.add_argument("--hanle-calibration-scale-points", type=int,
                   default=41)
    return apply_json_config(p.parse_args())


def main() -> None:
    os.environ.setdefault("OMP_NUM_THREADS", str(os.cpu_count() or 1))
    args = parse_args()
    settings = fidelity_settings(args.fidelity)
    configs = make_source_configs()
    if args.only:
        keep = set(args.only)
        configs = [cfg for cfg in configs if cfg.name in keep]
        missing = keep - {cfg.name for cfg in configs}
        if missing:
            raise SystemExit(f"unknown --only source(s): {sorted(missing)}")

    results = []
    print(f"Running {len(configs)} source configs at fidelity={args.fidelity}")
    for cfg in configs:
        print(f"  FWM: {cfg.name} ...", flush=True)
        source = compute_source(cfg, settings)
        src = source_readout_values(source, args)
        print(
            f"    model Gs={source['op']['G_s']:.3g}, "
            f"used {src['label']} Gs={src['G_s']:.3g}, "
            f"S={src['S_dB']:.2f} dB, "
            f"probe det87={source['probe_detuning_from_87_mhz']:+.1f} MHz "
            f"({source['runtime_s']:.1f}s)",
            flush=True)
        probe_power_w = arm_input_power_w(source, "probe", args)
        hanle = compute_hanle(source["probe_detuning_from_87_mhz"], args,
                              probe_power_w=probe_power_w)
        balanced = {
            mode: weighted_noise_over_hanle(source, hanle, args, mode)
            for mode in ("raw", "dc", "shot")
        }
        tolerance = (tolerance_sweeps(source, hanle,
                                      balanced["shot"]["best_idx"], args)
                     if balanced["shot"]["usable"] else None)
        optimization = None
        if not args.skip_hanle_optimize:
            print("    Hanle operating scan ...", flush=True)
            optimization = hanle_operating_scan(source, args)
            best = optimization["best"]
            if best is None:
                print("    no usable Hanle operating point in scan", flush=True)
            else:
                print(
                    f"    best Hanle: T={best['temp_c']:.1f} C, "
                    f"I={best['intensity_mw_cm2']:.3g} mW/cm^2, "
                    f"L={best['cell_mm']:.3g} mm, "
                    f"{best['sensitivity_pT']:.3g} pT/sqrtHz",
                    flush=True)
        results.append(dict(source=source, hanle=hanle, balanced=balanced,
                            tolerance=tolerance,
                            optimization=optimization))

    calibration = maybe_calibrate_hanle(results, args)
    if calibration is not None:
        fit = calibration["fit"]
        print(
            f"Hanle calibration: source={calibration['source_name']}, "
            f"B scale={fit['b_scale']:.6g}, "
            f"B offset={fit['b_offset_ut']:+.4g} uT, "
            f"norm RMS={fit['normalized_rms']:.3g}",
            flush=True)

    out = ROOT / "analysis" / "squeezing"
    out_npz = out / f"{OUT_STEM}.npz"
    out_png = out / f"{OUT_STEM}.png"
    out_tol_png = out / f"{OUT_STEM}_tolerance.png"
    out_opt_png = out / f"{OUT_STEM}_hanle_optimize.png"
    out_cal_png = out / f"{OUT_STEM}_hanle_calibration.png"
    out_md = out / f"{OUT_STEM}.md"
    save_npz(results, out_npz, args, calibration)
    plot_results(results, out_png)
    plot_tolerance(results, out_tol_png, args)
    wrote_opt = plot_hanle_optimization(results, out_opt_png)
    wrote_cal = plot_hanle_calibration(calibration, out_cal_png)
    if not wrote_cal and out_cal_png.exists():
        out_cal_png.unlink()
    write_report(results, args, out_md, calibration)
    print(f"saved {out_npz}")
    print(f"saved {out_png}")
    print(f"saved {out_tol_png}")
    if wrote_opt:
        print(f"saved {out_opt_png}")
    if wrote_cal:
        print(f"saved {out_cal_png}")
    print(f"saved {out_md}")


if __name__ == "__main__":
    main()
